#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
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
import multiprocessing as mp
from collections import Counter
import os
import numpy as np
import math
import random
from scipy.optimize import linear_sum_assignment
import sys
from .profiler import Profiler
import h5py
from contextlib import redirect_stderr
from io import StringIO

from .interpolationdatapoint import InterpolationDatapoint
from .atommapper import AtomMapper
from .molecule import Molecule
from .outputstream import OutputStream
from .veloxchemlib import mpi_master
from. veloxchemlib import hartree_in_kcalpermol, bohr_in_angstrom
from .errorhandler import assert_msg_critical
from .inputparser import (parse_input, print_keywords, print_attributes)

with redirect_stderr(StringIO()) as fg_err:
    import geometric

def process_data_point(args):
    index, data_point, cartesian_distance = args

    distance, weight_gradient, _, swapped = cartesian_distance(data_point)
    return abs(distance), (distance, index, weight_gradient, swapped)

class InterpolationDriver():
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
        - imforcefieldfile: File containing the necessary information to construct the interpolation forcefield
    """

    def __init__(self, z_matrix=None, comm=None, ostream=None):
        """
        Initializes the IMPES driver.
        """

        # if comm is None:
        #     comm = MPI.COMM_WORLD

        # if ostream is None:
        #     if comm.Get_rank() == mpi_master():
        #         ostream = OutputStream(sys.stdout)
        #     else:
        #         ostream = OutputStream(None)

        # # MPI information
        # self.comm = comm
        # self.rank = self.comm.Get_rank()
        # self.nodes = self.comm.Get_size()
        # ostream = OutputStream(sys.stdout)
        # self.ostream = ostream

        # Create InterpolationDatapoint, z-matrix required!
        self.impes_coordinate = InterpolationDatapoint(z_matrix)#, self.comm,self.ostream)

        self.dihedral_function = None
        # simple or Shepard interpolation
        self.interpolation_type = 'shepard'
        self.exponent_p = None
        self.scaling_time = False
        self.exponent_q = None
        self.confidence_radius = 0.5
        
        self.distance_thrsh = 2.0
        self.distance_thrsh = 2.0
        self.z_matrix = z_matrix
        self.impes_dict = None
        self.sum_of_weights = None
        self.print = False
        self.store_weights = False
        self.perform_mapping = True
        self.symmetry_sub_groups = None
        self.weights = []
        self.int_coord_weights = {}
        self.potentials = []

        # kpoint file with QM data
        self.imforcefieldfile = None
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
                'imforcefieldfile':
                    ('str', 'the name of the chk file with QM data'),
                'labels': ('seq_fixed_str', 'the list of QM data point labels'),
            }
        }

    # def print_keywords(self):
    #     """
    #     Prints the input keywords of the ImpesDriver.
    #     """

    #     print_keywords(self._input_keywords, self.ostream)

    # def print_attributes(self):
    #     """
    #     Prints the attributes of the ImpesDriver.
    #     """

    #     print_attributes(self._input_keywords, self.ostream)

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
            self.imforcefieldfile is not None,
            'ImpesDriver: Please provide a chekpoint file name.')

        h5f = h5py.File(self.imforcefieldfile, 'r')

        keys = h5f.keys()

        remove_from_label = ""
        z_matrix_label = ''
        if self.impes_coordinate.use_inverse_bond_length:
            remove_from_label += "_rinv"
            z_matrix_label += '_rinv'
        else:
            remove_from_label += "_r"
            z_matrix_label += '_r'
        remove_from_label += "_dihedral"
        remove_from_label += "_energy"

        z_matrix_label += '_dihedral'
        z_matrix_bonds = z_matrix_label + '_bonds'
        z_matrix_angles = z_matrix_label + '_angles'
        z_matrix_dihedrals = z_matrix_label + '_dihedrals'
        z_matrix_labels = [z_matrix_bonds, z_matrix_angles, z_matrix_dihedrals]
        z_matrix = []
            
        labels = []
        counter = 0
        for key in keys:
            if remove_from_label in key:

                label = key.replace(remove_from_label, "")
                if label not in labels:
                    labels.append(label)
                
                if counter == 0:
                    for label_obj in z_matrix_labels:

                        z_label = label + label_obj
                        current_z_list = [z_list.tolist() for z_list in list(h5f.get(z_label))]
                        z_matrix.extend(current_z_list)
                counter = 1


        h5f.close()
        return labels, z_matrix

    def define_impes_coordinate(self, coordinates):
        """Defines the current coordinate based on the coordinates array.
           The energy and gradient of this new configuration
           will be determined by interpolation.

            :param coordinates:
                a numpy array of Cartesian coordinates.
        """
        self.impes_coordinate.reset_coordinates_impes_driver(coordinates)

    def compute(self, molecule, qm_data_points=None, chk_file=None, labels=None):
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
        
        self.molecule = molecule
        self.define_impes_coordinate(molecule.get_coordinates_in_bohr())

        if labels:
            self.labels = labels
        if chk_file:
            self.imforcefield_file = chk_file

        if self.qm_data_points is None:
            self.qm_data_points = self.read_qm_data_points()

        if self.interpolation_type == 'simple':
            self.simple_interpolation()
        elif self.interpolation_type == 'shepard':
            self.shepard_interpolation()
        else:
            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)

    def simple_interpolation(self):
        """Performs a simple interpolation.

        :param qm_data_points:
            A list of InterpolationDatapoint corresponding to configurations
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

        distances_and_gradients = []
        min_distance = float('inf')
        self.time_step_reducer = False
        
        org_internal_coordinates = self.impes_coordinate.internal_coordinates_values.copy()
        org_b_matrix = self.impes_coordinate.b_matrix.copy()
        org_mask = np.arange(len(self.impes_coordinate.z_matrix))

        z_matrix_dict = { tuple(sorted(element)): i 
                  for i, element in enumerate(self.impes_coordinate.z_matrix) }

        # Determine weights, weight gradients,
        # energy values, and energy gradients for interpolation
        for i, data_point in enumerate(self.qm_data_points):
            
            distance, weight_gradient, _, swapped = self.cartesian_distance(data_point)

            if abs(distance) < min_distance:
                min_distance = abs(distance)

            distances_and_gradients.append((distance, i, weight_gradient, swapped))

        # min_distance, distances_and_gradients = self.process_in_parallel()

        used_labels = []
        close_distances = None
        close_distances = [
            (self.qm_data_points[index], distance, wg, index, swapped_tuple) 
            for distance, index, wg, swapped_tuple in distances_and_gradients 
            if abs(distance) <= min_distance + self.distance_thrsh]
        distances_from_points = [] 
        internal_distance_dict = []
        counter_idx = 0
        for qm_data_point, distance, weight_grad, label_idx, swapped_tuple_obj in close_distances:
            confidence_radius = 1.0
            if qm_data_point.confidence_radius is not None:
                confidence_radius = qm_data_point.confidence_radius
            new_coordintes = False
            reordered_int_coord_values = org_internal_coordinates.copy()
            mask = []
            org_b_matrix_cp = org_b_matrix.copy()
            if swapped_tuple_obj is not None:
                for i, element in enumerate(self.impes_coordinate.z_matrix):
                    # Otherwise, reorder the element
                    reordered_element = [swapped_tuple_obj[3].get(x, x) for x in element]
                    key = tuple(sorted(reordered_element))
                    z_mat_index = z_matrix_dict.get(key)
                    mask.append(z_mat_index)
                    # print(z_mat_index, key, element)
                    reordered_int_coord_values[i] = (float(org_internal_coordinates[z_mat_index]))
                
                org_b_matrix_cp[org_mask] = org_b_matrix[mask]
                org_b_matrix_cp = org_b_matrix_cp.reshape(len(self.impes_coordinate.z_matrix), len(self.molecule.get_labels()), 3)
                org_b_matrix_cp = org_b_matrix_cp[:, swapped_tuple_obj[1], :]
                org_b_matrix_cp = org_b_matrix_cp.reshape(len(self.impes_coordinate.z_matrix), len(self.molecule.get_labels()) * 3)
                # print('reordered internal coordinates', reordered_int_coord_values)
                new_coordintes = True
            weight = 1.0 / (distance**(2 * self.exponent_p))
            used_labels.append(label_idx)
            sum_weights += weight
            sum_weight_gradients += weight_grad
            cuurent_int_dist = []
            for i, element in enumerate(self.impes_coordinate.z_matrix):
                current_diff = qm_data_point.internal_coordinates_values[i] - reordered_int_coord_values[i]
                if len(element) == 4:
                    current_diff = np.sin(current_diff)
                cuurent_int_dist.append(current_diff)
            internal_distance_dict.append(np.array(cuurent_int_dist))

            
            potential = self.compute_potential(qm_data_point, reordered_int_coord_values)
            gradient = self.compute_gradient(qm_data_point, swapped_tuple_obj, reordered_int_coord_values, org_b_matrix_cp)
            gradients.append(gradient)
            potentials.append(potential)
            weights.append(weight)
            self.weights.append(weight)
            self.potentials.append(potential)
            weight_gradients.append(weight_grad)
            distances_from_points.append(distance)
            counter_idx += 1
            if new_coordintes:
                # print('B-matrix', self.impes_coordinate.internal_coordinates_values, '\n org B-Matrix \n', org_internal_coordinates)
                # print('\n difference B-matrix', np.linalg.norm(self.impes_coordinate.b_matrix[:, :] - org_b_matrix_cp[:, :]))
                # print('\n ', np.linalg.norm(np.array(reordered_int_coord_values) - np.array(self.impes_coordinate.internal_coordinates_values)))
                self.impes_coordinate.internal_coordinates_values = org_internal_coordinates
                self.impes_coordinate.b_matrix = org_b_matrix

        # Perform the interpolation and save the result in
        # self.impes_coordinate.energy and self.impes_coordinate.gradient
        self.impes_coordinate.energy = 0
        self.impes_coordinate.gradient = np.zeros((natms, 3))
        n_points = len(close_distances)
        weights = np.array(weights) / sum_weights

        for i in range(n_points):
            self.int_coord_weight[self.labels[used_labels[i]]] = weights[i]
            self.impes_coordinate.energy += weights[i] * potentials[i]
            self.impes_coordinate.gradient += (
                weights[i] * gradients[i] +
                potentials[i] * weight_gradients[i] / sum_weights -
                potentials[i] * weights[i] * sum_weight_gradients / sum_weights)



        
    def process_in_parallel(self):
        # Prepare the data for parallel processing.
        # Each task is a tuple: (index, data_point, reference to self.cartesian_distance)
    
        tasks = [(i, data_point, self.cartesian_distance) 
                 for i, data_point in enumerate(self.qm_data_points)]
        print('initialize tasks', tasks)
        # Use a Pool to map the tasks.
        with mp.Pool() as pool:
            results = pool.map(process_data_point, tasks)
        
        # Initialize your variables.
        min_distance = float('inf')
        distances_and_gradients = []
        
        # Combine results from each process.
        for abs_distance, info in results:
            if abs_distance < min_distance:
                min_distance = abs_distance
            distances_and_gradients.append(info)
        
        return min_distance, distances_and_gradients

    def shepard_interpolation(self):
        """Performs a simple interpolation.

        :param qm_data_points:
            A list of InterpolationDatapoint corresponding to configurations
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


        distances_and_gradients = []
        min_distance = float('inf')
        self.time_step_reducer = False
        
        org_internal_coordinates = self.impes_coordinate.internal_coordinates_values.copy()
        org_b_matrix = self.impes_coordinate.b_matrix.copy()
        org_mask = np.arange(len(self.impes_coordinate.z_matrix))

        z_matrix_dict = { tuple(sorted(element)): i 
                  for i, element in enumerate(self.impes_coordinate.z_matrix) }
        for i, data_point in enumerate(self.qm_data_points):
            
            distance, weight_gradient, _, swapped = self.cartesian_distance(data_point)

            if abs(distance) < min_distance:
                min_distance = abs(distance)

            distances_and_gradients.append((distance, i, weight_gradient, swapped))

        # min_distance, distances_and_gradients = self.process_in_parallel()

        used_labels = []
        close_distances = None
        close_distances = [
            (self.qm_data_points[index], distance, wg, index, swapped_tuple) 
            for distance, index, wg, swapped_tuple in distances_and_gradients 
            if abs(distance) <= min_distance + self.distance_thrsh + 10000]
        distances_from_points = [] 
        internal_distance_dict = []
        counter_idx = 0
        for qm_data_point, distance, weight_grad, label_idx, swapped_tuple_obj in close_distances:
            
            # confidence_radius = 1.0
            # if qm_data_point.confidence_radius is not None:
            #     confidence_radius = qm_data_point.confidence_radius
            new_coordintes = False
            reordered_int_coord_values = org_internal_coordinates.copy()
            mask = []
            org_b_matrix_cp = org_b_matrix.copy()
            if swapped_tuple_obj is not None:
                for i, element in enumerate(self.impes_coordinate.z_matrix):
                    # Otherwise, reorder the element
                    reordered_element = [swapped_tuple_obj[3].get(x, x) for x in element]
                    key = tuple(sorted(reordered_element))
                    z_mat_index = z_matrix_dict.get(key)
                    mask.append(z_mat_index)
                    # print(z_mat_index, key, element)
                    reordered_int_coord_values[i] = (float(org_internal_coordinates[z_mat_index]))
                
                org_b_matrix_cp[org_mask] = org_b_matrix[mask]
                org_b_matrix_cp = org_b_matrix_cp.reshape(len(self.impes_coordinate.z_matrix), len(self.molecule.get_labels()), 3)
                org_b_matrix_cp = org_b_matrix_cp[:, swapped_tuple_obj[1], :]
                org_b_matrix_cp = org_b_matrix_cp.reshape(len(self.impes_coordinate.z_matrix), len(self.molecule.get_labels()) * 3)
                # print('reordered internal coordinates', reordered_int_coord_values)
                new_coordintes = True
            denominator = (
                (distance / qm_data_point.confidence_radius)**(2 * self.exponent_p) +
                (distance / qm_data_point.confidence_radius)**(2 * self.exponent_q))
            weight = 1.0 / denominator
            used_labels.append(label_idx)
            sum_weights += weight
            sum_weight_gradients += weight_grad
            cuurent_int_dist = []
            for i, element in enumerate(self.impes_coordinate.z_matrix):
                current_diff = qm_data_point.internal_coordinates_values[i] - reordered_int_coord_values[i]
                if len(element) == 4:
                    current_diff = np.sin(current_diff)
                cuurent_int_dist.append(current_diff)
            internal_distance_dict.append(np.array(cuurent_int_dist))

            
            potential = self.compute_potential(qm_data_point, reordered_int_coord_values)
            gradient = self.compute_gradient(qm_data_point, swapped_tuple_obj, reordered_int_coord_values, org_b_matrix_cp)
            gradients.append(gradient)
            potentials.append(potential)
            weights.append(weight)
            self.weights.append(weight)
            self.potentials.append(potential)
            weight_gradients.append(weight_grad)
            distances_from_points.append(distance)
            counter_idx += 1
            if new_coordintes:
                # print('B-matrix', self.impes_coordinate.internal_coordinates_values, '\n org B-Matrix \n', org_internal_coordinates)
                # print('\n difference B-matrix', np.linalg.norm(self.impes_coordinate.b_matrix[:, :] - org_b_matrix_cp[:, :]))
                # print('\n ', np.linalg.norm(np.array(reordered_int_coord_values) - np.array(self.impes_coordinate.internal_coordinates_values)))
                self.impes_coordinate.internal_coordinates_values = org_internal_coordinates
                self.impes_coordinate.b_matrix = org_b_matrix
        
        n_points = len(close_distances)
        if self.print:
            for i in range(n_points - 1):
                print('\nint distane \n', i, distances_from_points[i], np.linalg.norm(internal_distance_dict[i]), '\n',
                      i + 1, distances_from_points[i+ 1], np.linalg.norm(internal_distance_dict[i+ 1]), '\n\n', internal_distance_dict[i], '\n\n', internal_distance_dict[i + 1], '\n\n', internal_distance_dict[i] - internal_distance_dict[i + 1])
        
                
        
        self.impes_coordinate.energy = 0
        self.impes_coordinate.gradient = np.zeros((natms, 3))
        self.impes_coordinate.NAC = np.zeros((natms, 3))
        weights = np.array(weights) / sum_weights
        self.sum_of_weights = sum_weights
        
        for i in range(n_points):
            
            if self.store_weights:
                self.int_coord_weights[self.labels[used_labels[i]]] = weights[i]
            self.impes_coordinate.energy += weights[i] * potentials[i]

            self.impes_coordinate.gradient += (
                weights[i] * gradients[i] +
                potentials[i] * weight_gradients[i] / sum_weights -
                potentials[i] * weights[i] * sum_weight_gradients / sum_weights)
        
        # print('Energy, Gradient', self.impes_coordinate.energy, self.impes_coordinate.gradient)
        
    def read_qm_data_points(self):
        """ Reads the QM data points to be used for interpolation
            from a checkpoint file and saves them to self.qm_data_points
            as a list of objects.
        """

        qm_data_points = []

        assert_msg_critical(
            self.imforcefieldfile is not None,
            'ImpesDriver: Please provide a chekpoint file name.')
       
        if not self.labels:
            labels, z_matrix = self.read_labels()
            self.labels = labels  
            self.z_matrix = z_matrix
        z_matrix = self.impes_coordinate.z_matrix
       
        for label in self.labels:
            data_point = InterpolationDatapoint(z_matrix)#, self.comm, self.ostream)
            data_point.update_settings(self.impes_dict)
            data_point.read_hdf5(self.imforcefieldfile, label)
            qm_data_points.append(data_point)

        self.qm_data_points = qm_data_points

        return qm_data_points

    def principal_angle(self, angle_rad):
        return (angle_rad + np.pi) % (2.0 * np.pi) - np.pi

    def compute_potential(self, data_point, current_internal_coordinates_values):
        """Calculates the potential energy surface at self.impes_coordinate
           based on the energy, gradient and Hessian of data_point.

           :param data_point:
                InterpolationDatapoint object.
        """

        energy = data_point.energy
        grad = data_point.internal_gradient
        hessian = data_point.internal_hessian
        dist_check = (current_internal_coordinates_values - data_point.internal_coordinates_values)

        for i, element in enumerate(self.impes_coordinate.z_matrix): 
            if len(element) == 4:
                    dist_check[i] = np.sin(dist_check[i])

        pes = (energy + np.matmul(dist_check.T, grad) +
               0.5 * np.linalg.multi_dot([dist_check.T, hessian, dist_check]))
        
        return pes

    def compute_gradient(self, data_point, mapping_list, current_internal_coordinates_values, current_b_matrix):
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
             

        dist_check = (current_internal_coordinates_values - data_point.internal_coordinates_values)
        dist_org = (current_internal_coordinates_values - data_point.internal_coordinates_values)
        for i, element in enumerate(self.impes_coordinate.z_matrix):
            if len(element) == 4:    
                dist_check[i] = np.sin(dist_check[i])
                grad[i] *= np.cos(dist_org[i])

        
        dist_hessian = np.matmul(dist_check.T, hessian)

        for i, element in enumerate(self.impes_coordinate.z_matrix):
            if len(element) == 4:
                dist_hessian[i] *= np.cos(dist_org[i])
            
        internal_gradient_hess_cos = grad + dist_hessian

        gradient = (np.matmul(current_b_matrix.T, internal_gradient_hess_cos)).reshape(natm, 3)

        if mapping_list is not None:

            gradient[mapping_list[1]] = gradient[mapping_list[0]]


        return gradient

    def simple_weight_gradient(self, distance_vector, distance):
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
    
    def trust_radius_weight_gradient(self, datapoint):
        
        confidence_radius = datapoint.confidence_radius
        distance, _, _, _ = self.cartesian_distance(datapoint)
        denominator = (
                (distance / confidence_radius)**(2 * self.exponent_p) +
                (distance / confidence_radius)**(2 * self.exponent_q))**2
        trust_radius_weight_gradient = -2.0 * (( -1.0 * (self.exponent_p * ((distance)**(2 * self.exponent_p) * (1/confidence_radius)**(2 * self.exponent_p)) / confidence_radius) - 
                                             (self.exponent_q * ((distance)**(2 * self.exponent_q) * (1/confidence_radius)**(2 * self.exponent_q)) / confidence_radius))) / denominator
        return trust_radius_weight_gradient

    def shepard_weight_gradient(self, distance_vector, distance, confidence_radius):
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
            (distance / confidence_radius)**(2 * self.exponent_p) +
            (distance / confidence_radius)**(2 * self.exponent_q))
        derivative_p = (self.exponent_p * distance_vector *
                        distance**(2 * self.exponent_p - 2) /
                        confidence_radius**(2 * self.exponent_p))
        derivative_q = (self.exponent_q * distance_vector *
                        distance**(2 * self.exponent_q - 2) /
                        confidence_radius**(2 * self.exponent_q))

        weight_gradient = (-1.0 * (2 * (derivative_p + derivative_q)) *
                           (1.0 / (denominator**2)))

        return weight_gradient
    
    def cartesian_distance(self, data_point):
        """Calculates and returns the cartesian distance between
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.

           :param data_point:
                InterpolationDatapoint object
        """
        # First, translate the cartesian coordinates to zero
        target_coordinates = data_point.calculate_translation_coordinates()
        reference_coordinates = (
            self.impes_coordinate.calculate_translation_coordinates())
        
        target_coordinates_core = target_coordinates.copy()
        reference_coordinates_core = reference_coordinates.copy()

        core_detected = False
        if self.perform_mapping and len(self.symmetry_sub_groups[1]) != 0 and (len(self.symmetry_sub_groups[0]) - len(self.symmetry_sub_groups[1][0]) > 1):
            
            target_coordinates_core = np.delete(target_coordinates, self.symmetry_sub_groups[1][0], axis=0)
            reference_coordinates_core = np.delete(reference_coordinates, self.symmetry_sub_groups[1][0], axis=0)
            core_detected = True

        # Then, determine the rotation matrix which
        # aligns data_point (target_coordinates)
        # to self.impes_coordinate (reference_coordinates)

        rotation_matrix = geometric.rotate.get_rot(target_coordinates_core,
                                                   reference_coordinates_core)

        # Rotate the data point
        rotated_coordinates = np.dot(rotation_matrix, target_coordinates.T).T
        # rotation_weights = geometric.rotate.get_rot(rotated_coordinates,
        #                                             target_coordinates)

        org_reference_coord = reference_coordinates.copy()
        org_mapping_list = [i for i in range(reference_coordinates.shape[0])]
        mapping_list = [[i for i in range(reference_coordinates.shape[0])]]
        swapped = None

        if core_detected and 1 == 2:
    
            non_group_atoms = []
            for i in range(len(self.molecule.get_labels())):
                if i not in self.symmetry_sub_groups[1][0]:
                    non_group_atoms.append(i)

            current_molecule_trunc_coords = np.delete(reference_coordinates, non_group_atoms, axis=0)
            data_point_molecule_trun_coords = np.delete(rotated_coordinates, non_group_atoms, axis=0)

            mapping_list, mapping_dict, assigned = self.perform_symmetry_assignment(self.symmetry_sub_groups[0], self.symmetry_sub_groups[1][0], 
                                                                      current_molecule_trunc_coords, data_point_molecule_trun_coords)
            if assigned:
                swapped = (org_mapping_list, mapping_list[0], reference_coordinates, mapping_dict)
            # mapping_list = atom_mapper.perform(atom_map_1=self.symmetry_sub_groups[0], env_groups_1=self.symmetry_sub_groups[1])

        reference_coordinates[org_mapping_list] = org_reference_coord[mapping_list[0]]
        # Calculate the Cartesian distance
        distance = (np.linalg.norm(rotated_coordinates - reference_coordinates))
        # Calculate the gradient of the interpolation weights
        # (required for energy gradient interpolation)
        distance_vector = (reference_coordinates - rotated_coordinates)
        
        distance_vector_norm = np.zeros(reference_coordinates.shape[0])
        for i in range(len(distance_vector_norm)):
            distance_vector_norm[i] += np.linalg.norm(distance_vector[i])
        confidence_radius = 1.0
        if data_point.confidence_radius is not None:
            confidence_radius = data_point.confidence_radius

        if distance < 1e-7:
            distance = 1e-5
            distance_vector[:] = 0
        if self.interpolation_type == 'shepard':
            weight_gradient = self.shepard_weight_gradient(
                distance_vector, distance, confidence_radius)
        elif self.interpolation_type == 'simple':
            weight_gradient = self.simple_weight_gradient(
                distance_vector, distance)

            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)

        return distance, weight_gradient, distance_vector_norm, swapped


    def determine_important_internal_coordinates(self, qm_energy, molecule, z_matrix, datapoints):
        """Determines the most important internal coordinates
           leading to the large deviation.
        """
        constraints = []
        for datapoint in datapoints:
            _, _, _, swapped = self.cartesian_distance(datapoint)
            coordinates_swapped = False
            if swapped is not None:
                self.define_impes_coordinate(swapped[2])
                coordinates_swapped = True

            internal_coord_elem_distance = []
            for elem_idx, element in enumerate(z_matrix):
                if len(element) == 4:
                    internal_coord_elem_distance.append(np.sin(self.impes_coordinate.internal_coordinates_values[elem_idx] - datapoint.internal_coordinates_values[elem_idx]))
                else:
                    internal_coord_elem_distance.append(self.impes_coordinate.internal_coordinates_values[elem_idx] - datapoint.internal_coordinates_values[elem_idx])

            N = len(z_matrix)

            partial_energies = np.zeros(N)

            # First handle linear parts
            for i in range(N):
                partial_energies[i] += internal_coord_elem_distance[i] * datapoint.internal_gradient[i]

            # Now handle quadratic parts (diagonal + cross)
            for i in range(N):
                # Diagonal
                partial_energies[i] += 0.5 * internal_coord_elem_distance[i] * datapoint.internal_hessian[i,i] * internal_coord_elem_distance[i]

                # Cross-terms: j>i to avoid double counting in sum
                for j in range(i+1, N):
                    cross = internal_coord_elem_distance[i] * datapoint.internal_hessian[i,j] * internal_coord_elem_distance[j]
                    # Split cross equally between i and j
                    partial_energies[i] += 0.5 * cross
                    partial_energies[j] += 0.5 * cross

            # The sum of partial_energies should match the total second-order approx:
            # total_energy_diff = sum(partial_energies)
            pred_E = self.compute_potential(datapoint, current_internal_coordinates_values=self.impes_coordinate.internal_coordinates_values)
            delta_E = abs(qm_energy - pred_E)
            single_energy_error = []
            weights = []

            for individual_contrib in partial_energies:
                
                energy_error = delta_E * (abs(individual_contrib) / sum(abs(partial_energies)))
                single_weight = abs(individual_contrib) / sum(abs(partial_energies))
                weights.append(single_weight)
                single_energy_error.append(energy_error)

            # Pair each energy error with its corresponding partial energy and internal coordinate
            contributions = list(zip(partial_energies, single_energy_error, weights, z_matrix))

            # Sort contributions by the energy error in descending order
            sorted_contributions = sorted(contributions, key=lambda x: x[1], reverse=True)

            # Print the sorted contributions with internal coordinates
            print('Delta E:', delta_E * hartree_in_kcalpermol())
            for contrib, error, ind_weight, coord in sorted_contributions[:]:
                if len(coord) == 2 and ind_weight == max(weights):
                    constraints.append(tuple(int(x) for x in coord))
                elif len(coord) == 3 and ind_weight == max(weights):
                    constraints.append(tuple(int(x) for x in coord))
                elif len(coord) == 4 and ind_weight > max(weights) * 0.7:
                    constraints.append(tuple(int(x) for x in coord))
                print(f'Internal Coordinate: {tuple(int(x) for x in coord)}, distance {internal_coord_elem_distance[z_matrix.index(coord)]}, Contribution: {contrib}, weight {ind_weight}, Error: {error * hartree_in_kcalpermol()}')
            print('Sum of Weights', sum(weights), sum(single_energy_error) * hartree_in_kcalpermol())
            
            if coordinates_swapped:
                self.define_impes_coordinate(molecule.get_coordinates_in_bohr())


        return constraints

    def get_energy(self):
        """ Returns the potential energy obtained by interpolation.
        """
        return self.impes_coordinate.energy

    def get_gradient(self):
        """ Returns the gradient obtained by interpolation.
        """
        return self.impes_coordinate.gradient
    

    def perform_symmetry_assignment(self, atom_map, sym_group, reference_group, datapoint_group):
        """ Performs the atom mapping. """

        new_map = atom_map.copy()
        mapping_dict = {}
        # cost = self.get_dihedral_cost(atom_map, sym_group, non_group_atoms)
        cost = np.linalg.norm(datapoint_group[:, np.newaxis, :] - reference_group[np.newaxis, :, :], axis=2)
        row, col = linear_sum_assignment(cost)
        assigned = False

        if not np.equal(row, col).all():
            assigned = True
            
            # atom_maps = self.linear_assignment_solver(cost)
            reordred_arr = np.array(sym_group)[col]
            new_map[sym_group] = new_map[reordred_arr]
            mapping_dict = {org: new for org, new in zip(np.array(sym_group), reordred_arr)}
        
        return [new_map], mapping_dict, assigned