#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from mpi4py import MPI
from os import cpu_count
import numpy as np
import sys
import h5py
from contextlib import redirect_stderr
from io import StringIO

from .impescoordinates import ImpesCoordinates
from .outputstream import OutputStream
from .veloxchemlib import mpi_master
from .errorhandler import assert_msg_critical
from .inputparser import (parse_input, print_keywords, print_attributes)

with redirect_stderr(StringIO()) as fg_err:
    import geometric


class ImpesDriver():
    """
    Implements the potential energy surface interpolation driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - z_matrix: the Z-matrix of indices defining the internal
          coordinates.
        - impes_coordinate: the new molecular geometry at which the energy
          is to be calculated by interpolation.
        - energy: the energy determined by interpolation.
        - gradient: the Cartesian gradient determined by interpolation.
        - exponent_p: exponent required in the interpolation procedure (both
          simple and Shepard interpolation).
        - exponent_q: exponent required in the Shepard interpolation procedure.
        - confidence_radius: confidence radius required by the Shepard
          interpolation procedure.
    """

    def __init__(self, z_matrix, comm=None, ostream=None):
        """
        Initializes the IMPES driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # MPI information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()
        ostream = OutputStream(sys.stdout)
        self.ostream = ostream

        # Create ImpesCoordinate, z-matrix required!
        self.impes_coordinate = ImpesCoordinates(z_matrix)#, self.comm,self.ostream)

        self.dihedral_function = None
        # simple or Shepard interpolation
        self.interpolation_type = 'shepard'
        self.exponent_p = None
        self.scaling_time = False
        self.exponent_q = None
        self.confidence_radius = 0.5
        
        self.distance_thrsh = 0.5
        self.z_matrix = z_matrix

        self.old_gradient = None

        # kpoint file with QM data
        self.checkpoint_file_name = None
        self.qm_data_points = None
        # Name lables for the QM data points
        self.labels = None

        self._input_keywords = {
            'im_settings': {
                'interpolation_type':
                    ('str', 'type of interpolation (simple/Shepard)'),
                'exponent_p': ('int', 'the main exponent'),
                'exponent_q': ('int', 'the additional exponent (Shepard IM)'),
                'confidence_radius': ('float', 'the confidence radius'),
                'checkpoint_file_name':
                    ('str', 'the name of the chk file with QM data'),
                'labels': ('seq_fixed_str', 'the list of QM data point labels'),
            }
        }

    def print_keywords(self):
        """
        Prints the input keywords of the ImpesDriver.
        """

        print_keywords(self._input_keywords, self.ostream)

    def print_attributes(self):
        """
        Prints the attributes of the ImpesDriver.
        """

        print_attributes(self._input_keywords, self.ostream)

    def update_settings(self, impes_dict=None):
        """
        Updates settings in the ImpesDriver.

        :param impes_dict:
            The input dictionary of settings for IMPES.
        """

        if impes_dict is None:
            impes_dict = {}

        self.impes_dict = dict(impes_dict)
        
        self.impes_coordinate.update_settings(impes_dict)

        im_keywords = {
            key: val[0]
            for key, val in self._input_keywords['im_settings'].items()
        }

        parse_input(self, im_keywords, impes_dict)

        if self.interpolation_type == 'shepard' and self.exponent_q is None:
            self.exponent_q = self.exponent_p / 2.0

    def read_labels(self):
        """
        Read data point labels from checkpoint file.

        :returns a list of data point labels.

        """

        assert_msg_critical(
            self.checkpoint_file_name is not None,
            'ImpesDriver: Please provide a chekpoint file name.')

        h5f = h5py.File(self.checkpoint_file_name, 'r')

        keys = h5f.keys()

        remove_from_label = ""
        if self.impes_coordinate.use_inverse_bond_length:
            remove_from_label += "_rinv"
        else:
            remove_from_label += "_r"
        if self.impes_coordinate.use_cosine_dihedral:
            remove_from_label += "_cosine"
        else:
            remove_from_label += "_dihedral"
        remove_from_label += "_energy"

        labels = []
        for key in keys:
            if remove_from_label in key:
                label = key.replace(remove_from_label, "")
                if label not in labels:
                    labels.append(label)

        h5f.close()
        return labels

    def define_impes_coordinate(self, coordinates):
        """Defines the current coordinate based on the coordinates array.
           The energy and gradient of this new configuration
           will be determined by interpolation.

            :param coordinates:
                a numpy array of Cartesian coordinates.
        """
        self.impes_coordinate.reset_coordinates_impes_driver(coordinates)

    def compute(self, molecule, qm_data_points=None, chk_file=None, labels=None, NACs=False):
        """Computes the energy and gradient by interpolation
           between pre-defined points.

            :param coordinates:
                a numpy array of Cartesian coordinates.
            :param fname:
                the name of the checkpoint file containing the pre-defined
                data points.
            :param labels:
                a list of data point labels (optional).
                if this parameter is None, all datapoints from fname will be
                used.
        """

        self.NACs = NACs
        self.molecule = molecule

        self.define_impes_coordinate(molecule.get_coordinates_in_bohr())

        if labels:
            self.labels = labels
        if chk_file:
            self.checkpoint_file_name = chk_file
        if qm_data_points is None:
            qm_data_points = self.read_qm_data_points()
        if self.interpolation_type == 'simple':
            self.simple_interpolation(qm_data_points)
        elif self.interpolation_type == 'shepard':
            self.shepard_interpolation(qm_data_points)
        else:
            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)

    def simple_interpolation(self, qm_data_points):
        """Performs a simple interpolation.

        :param qm_data_points:
            A list of ImpesCoordinates corresponding to configurations
            to be used for interpolation (which have been calculated
            quantum mechanically).
        """

        natms = self.impes_coordinate.cartesian_coordinates.shape[0]

        sum_weights = 0.0
        potentials = []
        gradients = []
        weights = []
        weight_gradients = []
        sum_weight_gradients = np.zeros((natms, 3))

        n_points = len(qm_data_points)

        # Determine weights, weight gradients,
        # energy values, and energy gradients for interpolation
        for data_point in qm_data_points:
            distance, weight_gradient = self.cartesian_distance(data_point)
            weight = 1.0 / (distance**(2 * self.exponent_p))
            sum_weights += weight
            sum_weight_gradients += weight_gradient
            potential = self.compute_potential(data_point)
            gradient = self.compute_gradient(data_point)
            gradients.append(gradient)
            potentials.append(potential)
            weights.append(weight)
            weight_gradients.append(weight_gradient)

        # Perform the interpolation and save the result in
        # self.impes_coordinate.energy and self.impes_coordinate.gradient
        self.impes_coordinate.energy = 0
        self.impes_coordinate.gradient = np.zeros((natms, 3))

        weights = np.array(weights) / sum_weights

        for i in range(n_points):
            self.impes_coordinate.energy += weights[i] * potentials[i]
            self.impes_coordinate.gradient += (
                weights[i] * gradients[i] +
                potentials[i] * weight_gradients[i] / sum_weights -
                potentials[i] * weights[i] * sum_weight_gradients / sum_weights)
        
    
    def shepard_interpolation(self, qm_data_points):
        """Performs a simple interpolation.

        :param qm_data_points:
            A list of ImpesCoordinates corresponding to configurations
            to be used for interpolation (which have been calculated
            quantum mechanically).
        """
        
        natms = self.impes_coordinate.cartesian_coordinates.shape[0]

        sum_weights = 0.0
        potentials = []
        NACs = []
        gradients = []
        weights = []
        weight_gradients = []
        sum_weight_gradients = np.zeros((natms, 3))


        distances_and_gradients = []
        min_distance = float('inf')
        self.time_step_reducer = False

        for i, data_point in enumerate(qm_data_points):
            
            distance, weight_gradient, _ = self.cartesian_distance(data_point)
            print('Distance in IMpes', distance)
            if abs(distance) < min_distance:
                min_distance = abs(distance)

            distances_and_gradients.append((distance, i, weight_gradient))   
        
        close_distances = None
        close_distances = [
            (qm_data_points[index], distance, wg) 
            for distance, index, wg in distances_and_gradients 
            if abs(distance) <= min_distance + self.distance_thrsh]
            
        for qm_data_point, distance, weight_grad in close_distances:
            denominator = (
                (distance / self.confidence_radius)**(2 * self.exponent_p) +
                (distance / self.confidence_radius)**(2 * self.exponent_q))
            weight = 1.0 / denominator

            sum_weights += weight
            sum_weight_gradients += weight_grad
            potential = self.compute_potential(qm_data_point)
            gradient = self.compute_gradient(qm_data_point)
            gradients.append(gradient)
            potentials.append(potential)
            weights.append(weight)
            weight_gradients.append(weight_grad)
            if self.NACs:
                nac = self.compute_NAC(qm_data_point)
                NACs.append(nac)
            
        n_points = len(close_distances)
        self.impes_coordinate.energy = 0
        self.impes_coordinate.gradient = np.zeros((natms, 3))
        self.impes_coordinate.NAC = np.zeros((natms, 3))
        weights = np.array(weights) / sum_weights
        

        for i in range(n_points):

            if self.NACs:
                self.impes_coordinate.NAC += weights[i] * NACs[i]
            self.impes_coordinate.energy += weights[i] * potentials[i]

            self.impes_coordinate.gradient += (
                weights[i] * gradients[i] +
                potentials[i] * weight_gradients[i] / sum_weights -
                potentials[i] * weights[i] * sum_weight_gradients / sum_weights)
        
        print('Energy and Gradient', self.impes_coordinate.energy, self.impes_coordinate.gradient, weight_gradients, sum_weight_gradients, sum_weights, weights)

    def read_qm_data_points(self):
        """ Reads the QM data points to be used for interpolation
            from a checkpoint file and saves them to self.qm_data_points
            as a list of objects.
        """

        qm_data_points = []

        assert_msg_critical(
            self.checkpoint_file_name is not None,
            'ImpesDriver: Please provide a chekpoint file name.')
       
        if not self.labels:
            labels = self.read_labels()
            self.labels = labels  
        z_matrix = self.impes_coordinate.z_matrix
       
        for label in self.labels:
            data_point = ImpesCoordinates(z_matrix)#, self.comm, self.ostream)
            data_point.update_settings(self.impes_dict)
            data_point.read_hdf5(self.checkpoint_file_name, label)
            qm_data_points.append(data_point)

        self.qm_data_points = qm_data_points

        return qm_data_points

    def principal_angle(self, angle_rad):
        return (angle_rad + np.pi) % (2.0 * np.pi) - np.pi

    def compute_potential(self, data_point):
        """Calculates the potential energy surface at self.impes_coordinate
           based on the energy, gradient and Hessian of data_point.

           :param data_point:
                ImpesCoordinates object.
        """

        energy = data_point.energy
        grad = data_point.internal_gradient
        hessian = data_point.internal_hessian
        dist_check = (self.impes_coordinate.internal_coordinates_values - data_point.internal_coordinates_values)

        # implement sin for keeping the structure
        # change the code so it can take on any function that is being used for the dihedral
        # TODO: check if the bond angle can be described by a function as well?
        for i, element in enumerate(self.impes_coordinate.z_matrix): 
            if len(element) == 4:
                    dist_check[i] = np.sin(dist_check[i])
                    # dist_check[i] = abs(self.principal_angle(dist_check[i]))

        pes = (energy + np.matmul(dist_check.T, grad) +
               0.5 * np.linalg.multi_dot([dist_check.T, hessian, dist_check]))
        
        return pes

    
    def compute_gradient(self, data_point):
        """Calculates part of the cartesian gradient
           for self.impes_coordinate based on the gradient
           and Hessian of data_point.
           J ( g_a + delta_a H_a )
           eq. (4) JCTC 12, 5235-5246

           :param data_point:
                ImpesCoordinates object.
        """

        natm = data_point.cartesian_coordinates.shape[0]
        grad = data_point.internal_gradient.copy()
        hessian = data_point.internal_hessian
        im_b_matrix = self.impes_coordinate.b_matrix                
        dist_org = (self.impes_coordinate.internal_coordinates_values - data_point.internal_coordinates_values)
        dist_check_hessian = (self.impes_coordinate.internal_coordinates_values - data_point.internal_coordinates_values)
        
        for i, element in enumerate(self.impes_coordinate.z_matrix):
            
            if len(element) == 4:
                dist_check_hessian[i] = np.sin(dist_check_hessian[i])

        dist_hessian_cos = np.matmul(dist_check_hessian.T, hessian)

        internal_gradient_hess_cos = grad + dist_hessian_cos
        
        for i, element in enumerate(self.impes_coordinate.z_matrix):
            if len(element) == 4:
                internal_gradient_hess_cos[i] *= np.cos(dist_org[i])
            if len(element) == 2:
                internal_gradient_hess_cos[i] *= -(self.impes_coordinate.internal_coordinates_values[i])**2
            else:
                internal_gradient_hess_cos[i] *= 1.0

        gradient = (np.matmul(im_b_matrix.T, internal_gradient_hess_cos)).reshape(natm, 3)

        self.old_gradient = gradient
        return gradient

    def simple_weight_gradient(self, distance_vector, distance,
                               rotation_matrix):
        """ Returns the derivative of an unormalized simple interpolation
            weight with respect to the Cartesian coordinates in
            self.impes_coordinate

            :param distance_vector:
                The Cartesian distance vector between
                the current data_point and self.impes_coordinate.
            :param distance:
                The norm of the distance vector * sqrt(N), N number of atoms.
        """
        weight_gradient = (-2 * self.exponent_p * distance_vector /
                           (distance**(2 * self.exponent_p + 2)))

        return weight_gradient

    def shepard_weight_gradient(self, distance_vector, distance):
        """ Returns the derivative of an unormalized Shepard interpolation
            weight with respect to the Cartesian coordinates in
            self.impes_coordinate

            :param distance_vector:
                The Cartesian distance vector between
                the current data_point and self.impes_coordinate.
            :param distance:
                The norm of the distance vector * sqrt(N), N number of atoms.
        """

        denominator = (
            (distance / self.confidence_radius)**(2 * self.exponent_p) +
            (distance / self.confidence_radius)**(2 * self.exponent_q))
        derivative_p = (self.exponent_p * distance_vector *
                        distance**(2 * self.exponent_p - 2) /
                        self.confidence_radius**(2 * self.exponent_p))
        derivative_q = (self.exponent_q * distance_vector *
                        distance**(2 * self.exponent_q - 2) /
                        self.confidence_radius**(2 * self.exponent_q))

        weight_gradient = (-1.0 * (2 * (derivative_p + derivative_q)) *
                           (1.0 / (denominator**2)))

        return weight_gradient
    
    
    def cartesian_distance(self, data_point):
        """Calculates and returns the cartesian distance between
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.

           :param data_point:
                ImpesCoordinates object
        """
        # First, translate the cartesian coordinates to zero
        target_coordinates = data_point.calculate_translation_coordinates()
        reference_coordinates = (
            self.impes_coordinate.calculate_translation_coordinates())

        # Then, determine the rotation matrix which
        # aligns data_point (target_coordinates)
        # to self.impes_coordinate (reference_coordinates)
        rotation_matrix = geometric.rotate.get_rot(target_coordinates,
                                                   reference_coordinates)

        # Rotate the data point
        rotated_coordinates = np.dot(rotation_matrix, target_coordinates.T).T
        # rotation_weights = geometric.rotate.get_rot(rotated_coordinates,
        #                                             target_coordinates)

        # Calculate the Cartesian distance
        distance = (np.linalg.norm(rotated_coordinates - reference_coordinates))

        # Calculate the gradient of the interpolation weights
        # (required for energy gradient interpolation)
        distance_vector = (reference_coordinates - rotated_coordinates)

        distance_vector_norm = np.zeros(reference_coordinates.shape[0])
        for i in range(len(distance_vector_norm)):
            distance_vector_norm[i] += np.linalg.norm(distance_vector[i])

        if distance < 1e-10:
            distance = 1e-7
            distance_vector[:] = 0
        if self.interpolation_type == 'shepard':
            weight_gradient = self.shepard_weight_gradient(
                distance_vector, distance)
        elif self.interpolation_type == 'simple':
            weight_gradient = self.simple_weight_gradient(
                distance_vector, distance, None)
        else:
            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)

        return distance, weight_gradient, distance_vector_norm


    def get_energy(self):
        """ Returns the potential energy obtained by interpolation.
        """
        return self.impes_coordinate.energy

    def get_gradient(self):
        """ Returns the gradient obtained by interpolation.
        """
        return self.impes_coordinate.gradient
