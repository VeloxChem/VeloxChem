#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

from mpi4py import MPI
import multiprocessing as mp
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
        - imforcefield_file: File containing the necessary information to construct the interpolation forcefield
    """

    def __init__(self, z_matrix=None, comm=None, ostream=None):
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
        self.z_matrix = z_matrix
        self.impes_dict = None

        self.store_weights = False
        self.perform_mapping = True
        self.symmetry_sub_groups = None
        self.weights = {}

        # kpoint file with QM data
        self.imforcefield_file = None
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
                'imforcefield_file':
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
            self.imforcefield_file is not None,
            'ImpesDriver: Please provide a chekpoint file name.')

        h5f = h5py.File(self.imforcefield_file, 'r')

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

        n_points = len(self.qm_data_points)

        # Determine weights, weight gradients,
        # energy values, and energy gradients for interpolation
        for data_point in self.qm_data_points:
            distance, weight_gradient, _ = self.cartesian_distance(data_point)
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

        for i, data_point in enumerate(self.qm_data_points):
            
            distance, weight_gradient, _, swapped = self.cartesian_distance(data_point)

            if abs(distance) < min_distance:
                min_distance = abs(distance)

            distances_and_gradients.append((distance, i, weight_gradient, swapped))   
        
        used_labels = []
        close_distances = None
        close_distances = [
            (self.qm_data_points[index], distance, wg, index, swapped_tuple) 
            for distance, index, wg, swapped_tuple in distances_and_gradients 
            if abs(distance) <= min_distance + self.distance_thrsh]
        distances_from_points = [] 
        for qm_data_point, distance, weight_grad, label_idx, swapped_tuple_obj in close_distances:
            confidence_radius = 1.0
            if qm_data_point.confidence_radius is not None:
                confidence_radius = qm_data_point.confidence_radius
            coordiantes_swapped = False
            if swapped_tuple_obj is not None:
                self.define_impes_coordinate(swapped_tuple_obj[2])
                coordiantes_swapped = True
            denominator = (
                (distance / confidence_radius)**(2 * self.exponent_p) +
                (distance / confidence_radius)**(2 * self.exponent_q))
            weight = 1.0 / denominator
            used_labels.append(label_idx)
            sum_weights += weight
            sum_weight_gradients += weight_grad

            potential = self.compute_potential(qm_data_point)
            gradient = self.compute_gradient(qm_data_point, swapped_tuple_obj)
            gradients.append(gradient)
            potentials.append(potential)
            weights.append(weight)
            weight_gradients.append(weight_grad)
            distances_from_points.append(distance)
            if coordiantes_swapped:
                self.define_impes_coordinate(self.molecule.get_coordinates_in_bohr())
            
        n_points = len(close_distances)
        self.impes_coordinate.energy = 0
        self.impes_coordinate.gradient = np.zeros((natms, 3))
        self.impes_coordinate.NAC = np.zeros((natms, 3))
        weights = np.array(weights) / sum_weights
        

        for i in range(n_points):
            if self.store_weights:
                self.weights[self.labels[used_labels[i]]] = weights[i]
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
            self.imforcefield_file is not None,
            'ImpesDriver: Please provide a chekpoint file name.')
       
        if not self.labels:
            labels, z_matrix = self.read_labels()
            self.labels = labels  
            self.z_matrix = z_matrix
        z_matrix = self.impes_coordinate.z_matrix
       
        for label in self.labels:
            data_point = InterpolationDatapoint(z_matrix)#, self.comm, self.ostream)
            data_point.update_settings(self.impes_dict)
            data_point.read_hdf5(self.imforcefield_file, label)
            qm_data_points.append(data_point)

        self.qm_data_points = qm_data_points

        return qm_data_points

    def principal_angle(self, angle_rad):
        return (angle_rad + np.pi) % (2.0 * np.pi) - np.pi

    def compute_potential(self, data_point):
        """Calculates the potential energy surface at self.impes_coordinate
           based on the energy, gradient and Hessian of data_point.

           :param data_point:
                InterpolationDatapoint object.
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

    def compute_gradient(self, data_point, mapping_list):
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

        dist_check = (self.impes_coordinate.internal_coordinates_values - data_point.internal_coordinates_values)
        dist_org = (self.impes_coordinate.internal_coordinates_values - data_point.internal_coordinates_values)
        for i, element in enumerate(self.impes_coordinate.z_matrix):
            if len(element) == 4:    
                dist_check[i] = np.sin(dist_check[i])
                grad[i] *= np.cos(dist_org[i])

        
        dist_hessian = np.matmul(dist_check.T, hessian)

        for i, element in enumerate(self.impes_coordinate.z_matrix):
            if len(element) == 4:
                dist_hessian[i] *= np.cos(dist_org[i])
            
        internal_gradient_hess_cos = grad + dist_hessian

        gradient = (np.matmul(im_b_matrix.T, internal_gradient_hess_cos)).reshape(natm, 3)

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

        if self.perform_mapping and len(self.symmetry_sub_groups) != 0:
            
            target_coordinates_core = np.delete(target_coordinates, self.symmetry_sub_groups[0], axis=0)
            reference_coordinates_core = np.delete(reference_coordinates, self.symmetry_sub_groups[0], axis=0)

        # Then, determine the rotation matrix which
        # aligns data_point (target_coordinates)
        # to self.impes_coordinate (reference_coordinates)

        rotation_matrix = geometric.rotate.get_rot(target_coordinates_core,
                                                   reference_coordinates_core)

        # Rotate the data point
        rotated_coordinates = np.dot(rotation_matrix, target_coordinates.T).T
        # rotation_weights = geometric.rotate.get_rot(rotated_coordinates,
        #                                             target_coordinates)

        rotated_coordinates_cp = rotated_coordinates.copy()
        org_reference_coord = reference_coordinates.copy()
        org_mapping_list = [i for i in range(reference_coordinates.shape[0])]
        mapping_list = [[i for i in range(reference_coordinates.shape[0])]]
        if len(self.symmetry_sub_groups) != 0:
            current_molecule = Molecule(self.molecule.get_labels(), reference_coordinates, 'bohr')
            data_point_molecule = Molecule(self.molecule.get_labels(), rotated_coordinates_cp, 'bohr')
            

            atom_mapper = AtomMapper(data_point_molecule, current_molecule)

            mapping_list = atom_mapper.perform()

        reference_coordinates[org_mapping_list] = org_reference_coord[mapping_list[0]]
        # Calculate the Cartesian distance
        distance = (np.linalg.norm(rotated_coordinates - reference_coordinates))
        # Calculate the gradient of the interpolation weights
        # (required for energy gradient interpolation)
        distance_vector = (reference_coordinates - rotated_coordinates)
        swapped = None
        if list(org_mapping_list) != list(mapping_list[0]) :
            swapped = (org_mapping_list, mapping_list[0], reference_coordinates)
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
                distance_vector, distance, confidence_radius)
        else:
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
            pred_E = self.compute_potential(datapoint)
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
                if ind_weight > max(weights) * 0.7:
                    constraints.append(tuple(int(x) for x in coord))
                print(f'Internal Coordinate: {tuple(int(x) for x in coord)}, Contribution: {contrib}, weight {ind_weight}, Error: {error * hartree_in_kcalpermol()}')
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
    
