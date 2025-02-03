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
import random
from scipy.optimize import linear_sum_assignment
import sys
from .profiler import Profiler
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

        # Create ImpesCoordinate, z-matrix required!
        self.impes_coordinate = ImpesCoordinates(z_matrix)#, self.comm,self.ostream)

        self.dihedral_function = None
        # simple or Shepard interpolation
        self.interpolation_type = 'shepard'
        self.exponent_p = None
        self.scaling_time = False
        self.exponent_q = None
        self.confidence_radius = 0.5
        self.available_resources = None

        self.individual_weights = []
        self.important_weights = None
        self.distance_thrsh = 0.5
        self.z_matrix = z_matrix

        # kpoint file with QM data
        self.checkpoint_file_name = None
        self.qm_data_points = None
        # Name lables for the QM data points
        self.labels = None
        self.NACs = None
        self.swapping_map_bool = False
        self.swapped_indices_map = None

        self.core_structure = None
        self.non_core_structure = None
        self.non_core_symmetry_group = None
        self.non_core_symmetry_group_mapping = None
        self.indicies_swapping = {}
        self.swapped_indices_list = None

        self.data_points_gradient = None
        self.data_points_hessian = None

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

    def compute(self, molecule, qm_data_points=None, chk_file=None, labels=None, NACs=False, available_resources=None):
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

        self.available_resources = available_resources or mp.cpu_count()
        self.NACs = NACs
        self.molecule = molecule
        self.define_impes_coordinate(molecule.get_coordinates_in_bohr())
        qm_data_points_1 = qm_data_points
        if labels:
            self.labels = labels
        if chk_file:
            self.checkpoint_file_name = chk_file
        if qm_data_points_1 is None:
            self.read_qm_data_points()
            qm_data_points_1 = self.qm_data_points
        if self.interpolation_type == 'simple':
            self.simple_interpolation(qm_data_points_1)
        elif self.interpolation_type == 'shepard':
            if len(qm_data_points_1) > 800:
                self.shepard_interpolation(qm_data_points_1)
            else:
                self.shepard_interpolation_dsitance_check(qm_data_points_1)
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

    def guassian_process_regression(self, qm_data_points):
        '''This function uses a Gaussian approach in order
           to interpolate between different datapoints
        '''

        def kernel(self, data_point):
            '''defines the Kernel function which is covarince and
               allow for a smooth interpolation between the individual
               datapoints which is govern by the distance and the hyper-
               paramters which can be optimally paramaetrized
            '''
            
            sigma = 1.0
            l = 1.0

            distance, _, _ = self.cartesian_distance(data_point)

            kernel_function = sigma**2 * (1 + ((np.sqrt(5) * distance)/(l)) + ((5 * distance)**2/(3 * l**2))) * np.exp(-((np.sqrt(5) * distance)/(l)))


            return kernel_function
        
        
        natms = self.impes_coordinate.cartesian_coordinates.shape[0]
        potentials = []
        NACs = []
        gradients = []

        n_points = len(qm_data_points)

    def computation_of_distances(self, data_point):

        distance, weight_grad, _ = self.cartesian_distance(data_point)

        return {
        'distance': distance,
        'weight_gradient': weight_grad
    }

    def computation_of_distances_analysis(self, data_point):

        distance, weight_grad, _, swapped_indices_map = self.cartesian_just_distance_analysis(self.impes_coordinate, data_point)

        return {
        'distance': distance,
        'weight_gradient': weight_grad,
        'swapped_indices_map':swapped_indices_map
    }

    def computation_of_interpolation_properties(self, args):
        
        distance, qm_data_point, weight_gradient = args
        
        denominator = (
            (distance / self.confidence_radius)**(2 * self.exponent_p) +
            (distance / self.confidence_radius)**(2 * self.exponent_q))
        
        weight = 1.0 / denominator
        #sum_weights += weight
        #sum_weight_gradients += weight_gradient
        potential = self.compute_potential(qm_data_point)
        #print('potential in impesdriver class', potential)
        gradient = self.compute_gradient(qm_data_point)
        # gradients.append(gradient)
        # potentials.append(potential)
        # weights.append(weight)
        # weight_gradients.append(weight_gradient)
        if self.NACs:
            nac = self.compute_NAC(qm_data_point)
            # NACs.append(nac)
        return {
        'weight': weight,
        'potential': potential,
        'gradient': gradient,
        'weight_grad': weight_gradient}
    
    def shepard_interpolation_paralell(self, qm_data_points):
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

        #n_points = len(qm_data_points)

        # Determine weights, weight gradients,
        # energy values, and energy gradients for interpolation

        counter = 0
        point_density = 0
        prev_gradient = 1000000
        prev_weight_gradient = -1e7
        distances_and_gradients = []
        min_distance = float('inf')
        self.time_step_reducer = False
        num_processes = min(self.available_resources, len(close_distances))
        
        with mp.Pool(processes=8) as internal_pool:
            distances_grad = internal_pool.map(self.computation_of_distances, qm_data_points)

        distances = [distance['distances'] for distance in distances_grad]
        weight_grads = [weight_grad['weight_gradient'] for weight_grad in distances_grad]
        for i, distance in enumerate(distances):
            if abs(distance) < min_distance:
                min_distance = abs(distance)
 
            
            distances_and_gradients.append((distance, i))
        
        close_distances = [
        (distance, qm_data_points[index], weight_grads[index]) 
        for distance, index in distances_and_gradients 
        if abs(distance) <= min_distance + 0.5  
        ]

        
        with mp.Pool(processes=8) as internal_pool:
            results = internal_pool.map(self.computation_of_interpolation_properties, close_distances)

        weights = [result['weight'] for result in results]
        potentials = [result['potential'] for result in results]
        gradients = [result['gradient'] for result in results]
        weight_gradients = [result['weight_gradient'] for result in results]
        for i, weight in enumerate(weights):
            sum_weights += weight
            sum_weight_gradients += weight_gradients[i]
        n_points = len(close_distances)
        self.impes_coordinate.energy = 0
        self.impes_coordinate.gradient = np.zeros((natms, 3))
        self.impes_coordinate.NAC = np.zeros((natms, 3))
        weights = np.array(weights) / sum_weights
        for i in range(n_points):

            if self.NACs:
                self.impes_coordinate.NAC += weights[i] * NACs[i]
            self.impes_coordinate.energy += weights[i] * potentials[i]
            #print(weights[i], self.labels[i])
            print('WEIGHT', weights[i])
            self.impes_coordinate.gradient += (
                weights[i] * gradients[i] +
                potentials[i] * weight_gradients[i] / sum_weights -
                potentials[i] * weights[i] * sum_weight_gradients / sum_weights)
        print('here is the gradient', potentials, np.linalg.norm(self.impes_coordinate.gradient), self.impes_coordinate.energy, '\n\n')

    
    def calculate_single_data_point_properties(self, args):
        
        current_datapoint, distance, weight_gradient, swapped_indices_map, current_gradient = args

        result = {
            'potential': 0,
            'gradient': 0,
            'weight': 0,
            'weight_gradient': np.zeros((self.impes_coordinate.cartesian_coordinates.shape[0], 3))
        }

        denominator = (
            (distance / (self.confidence_radius * current_gradient))**(2 * self.exponent_p / current_gradient) +
            (distance / (self.confidence_radius * current_gradient))**(2 * self.exponent_q / current_gradient))
        weight = 1.0 / denominator
        # print('self.swapped_indices', len(self.swapped_indices_map))
        # if self.swapped_indices_map is not None:
            # self.define_impes_coordinate(self.swapped_indices_map[current_datapoint.point_label][0])
        if swapped_indices_map:
            self.define_impes_coordinate(swapped_indices_map[0])
        potential = self.compute_potential(current_datapoint)
        gradient = self.compute_gradient(current_datapoint, swapped_indices_map=swapped_indices_map)
        result['weight'] = weight
        result['weight_gradient'] = weight_gradient
        result['potential'] = potential
        result['impes_potential'] = weight * potential
        result['impes_gradient_1'] =  weight * gradient
        result['impes_gradient_2'] = weight_gradient * potential
        result['impes_gradient_3'] = weight * potential

        return result
        
    
    def shepard_interpolation(self, qm_data_points, core=False):
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
        distances_and_gradients = []
        gradients = []
        weights = []
        weight_gradients = []
        sum_weight_gradients = np.zeros((natms, 3))

        distances_and_gradients = []
        min_distance = float('inf')

        n_points = len(qm_data_points)
        ostream = OutputStream(sys.stdout)
        self.profiling = False
        self.memory_profiling = False
        self.memory_tracing = False
        self.timing = True

        profiler = Profiler({
                'timing': self.timing,
                'profiling': self.profiling,
                'memory_profiling': self.memory_profiling,
                'memory_tracing': self.memory_tracing,
            })
        profiler.set_timing_key(f'Full interpolation time')
        profiler.start_timer('Full Loop interpolation')

        # Determine weights, weight gradients,
        # energy values, and energy gradients for interpolation
        close_distances = None
        if len(self.non_core_symmetry_group) > 0:
            core = True
        # core = False
        print('Number of qm datapoints is given here:', len(qm_data_points))
        num_workers = min(os.cpu_count(), len(qm_data_points))
        self.swapped_indices_map  = None
        if core is True:
            self.swapped_indices_map  = {}
        distances_grad = None

        if not core:
            with mp.Pool(processes=4) as internal_pool:
                distances_grad = internal_pool.map(self.computation_of_distances, [qm_datapoint for qm_datapoint in qm_data_points])
            # distance, weight_gradient, _ = self.cartesian_distance(current_datapoint)
            distances = [distance['distance'] for distance in distances_grad]
            weight_grads = [weight_grad['weight_gradient'] for weight_grad in distances_grad]
            for i, distance in enumerate(distances):
                if abs(distance) < min_distance:
                    min_distance = abs(distance)
                
                distances_and_gradients.append((distance, i))
            
            close_distances = [
            (qm_data_points[index], distance, weight_grads[index], None, 1.0) 
            for distance, index in distances_and_gradients 
            if abs(distance) <= min_distance + self.distance_thrsh
        ]
        else:
            profiler.start_timer('distance check')
            with mp.Pool(processes=4) as internal_pool:
                distances_grad = internal_pool.map(self.computation_of_distances_analysis, [qm_datapoint for qm_datapoint in qm_data_points])
            # distance, weight_gradient, _ = self.cartesian_just_distance_analysis(self.impes_coordinate, current_datapoint)
            profiler.stop_timer('distance check')
            distances = [distance['distance'] for distance in distances_grad]
            weight_grads = [weight_grad['weight_gradient'] for weight_grad in distances_grad]
            swapped_indices = [swapped_indices['swapped_indices_map'] for swapped_indices in distances_grad]
            for i, distance in enumerate(distances):
                if abs(distance) < min_distance:
                    min_distance = abs(distance)
                
                distances_and_gradients.append((distance, i))
            
            close_distances = [
            (qm_data_points[index], distance, weight_grads[index], swapped_indices[index], 1.0) 
            for distance, index in distances_and_gradients 
            if abs(distance) <= min_distance + self.distance_thrsh
            ]

        # if self.data_points_gradient is not None:
        #     indices_for_mapping = [self.data_points_gradient[i] for distance, i in distances_and_gradients if abs(distance) <= min_distance + self.distance_thrsh]

        #     # print('here are the indices', indices_for_mapping, '\n\n', self.data_points_gradient, 'MIN', min(self.data_points_gradient), '\n\n HESS max', self.data_points_hessian['max'], 'min', self.data_points_hessian['min'])

        #     confidence_radii = [0.3 * (self.data_points_gradient[i] /  min(self.data_points_gradient)) for distance, i in distances_and_gradients if abs(distance) <= min_distance + self.distance_thrsh]

            # print('confidence radii \n', confidence_radii)

        num_workers = max(1, int(min(os.cpu_count(), len(close_distances)) / 5))
        print('here is the number of qm_datapoints which are left', len(close_distances), num_workers)
        # Use the minimum of available CPU cores or number of data points
        profiler.start_timer('calc datapoint properties')
        with mp.Pool(processes=num_workers) as pool:
            results = pool.map(self.calculate_single_data_point_properties, [(close_distance) for close_distance in close_distances])
        profiler.stop_timer('calc datapoint properties')
        self.impes_coordinate.energy = 0
        gradient_1 = np.zeros((natms, 3))
        gradient_2 = np.zeros((natms, 3))
        gradient_3 = 0
        self.impes_coordinate.NAC = np.zeros((natms, 3))
        available_potentials = []
        for result in results:
            sum_weights += result['weight']
            sum_weight_gradients += result['weight_gradient']
            available_potentials.append((result['potential'], result['weight']))
            self.impes_coordinate.energy += result['impes_potential']
            gradient_1 += result['impes_gradient_1']
            gradient_2 += result['impes_gradient_2']
            gradient_3 += result['impes_gradient_3']
        

        self.impes_coordinate.energy /= sum_weights
        self.impes_coordinate.gradient = gradient_1 / sum_weights + gradient_2 / sum_weights - gradient_3 / sum_weights * sum_weight_gradients / sum_weights
        
        best_potential = np.inf

        # for i, pot in enumerate(available_potentials):
        #     print('Here is the potential', pot[0] * 627.509474, best_potential * 627.509474, pot[1]/sum_weights)
        #     if pot[0] < best_potential:
        #         best_potential = pot[0]
        #         best_index = i
        #         print('Here is potential', best_potential * 627.509474, best_index, pot[1] / sum_weights)
                
        
        profiler.stop_timer('Full Loop interpolation')
        # profiler.print_timing(ostream)
        # ostream.flush()
        # check_pot = self.impes_coordinate.energy
        # check_norm_gradient = np.linalg.norm(self.impes_coordinate.gradient)
        # print('For the comparison 1', self.impes_coordinate.energy, np.linalg.norm(self.impes_coordinate.gradient))
        
        # sum_weights = 0.0
        # potentials = []
        # NACs = []
        # distances_and_gradients = []
        # gradients = []
        # weights = []
        # weight_gradients = []
        # sum_weight_gradients = np.zeros((natms, 3))

        # distances_and_gradients = []
        # min_distance = float('inf')

        # n_points = len(qm_data_points)
        
        # for i, data_point in enumerate(qm_data_points):
            
        #     if core is False:
        #         distance, weight_gradient, _ = self.cartesian_distance(data_point)
        #         denominator = (
        #             (distance / self.confidence_radius)**(2 * self.exponent_p) +
        #             (distance / self.confidence_radius)**(2 * self.exponent_q))
        #         weight = 1.0 / denominator
        #         if self.swapped_indices_map is not None:
        #             self.define_impes_coordinate(self.swapped_indices_map[data_point.point_label][0])
        #         sum_weights += weight
        #         sum_weight_gradients += weight_gradient
        #         potential = self.compute_potential(data_point)
        #         gradient = self.compute_gradient(data_point)
        #         gradients.append(gradient)
        #         potentials.append(potential)
        #         weights.append(weight)
        #         weight_gradients.append(weight_gradient)
        #         # if self.NACs:
        #         #     nac = self.compute_NAC(qm_data_points[i])
        #         #     NACs.append(nac)

        #     else:
        #         distance, weight_gradient, _, _ = self.cartesian_just_distance_analysis(self.impes_coordinate, data_point)

        #         denominator = (
        #             (distance / self.confidence_radius)**(2 * self.exponent_p) +
        #             (distance / self.confidence_radius)**(2 * self.exponent_q))
        #         weight = 1.0 / denominator
        #         if self.swapped_indices_map is not None:
        #             self.define_impes_coordinate(self.swapped_indices_map[data_point.point_label][0])
        #         sum_weights += weight
        #         sum_weight_gradients += weight_gradient
        #         potential = self.compute_potential(data_point)
        #         gradient = self.compute_gradient(data_point)
        #         gradients.append(gradient)
        #         potentials.append(potential)
        #         weights.append(weight)
        #         weight_gradients.append(weight_gradient)
        #         # if self.NACs:
        #         #     nac = self.compute_NAC(qm_data_points[i])
        #         #     NACs.append(nac)
            

        # self.impes_coordinate.energy = 0
        # self.impes_coordinate.gradient = np.zeros((natms, 3))
        # self.impes_coordinate.NAC = np.zeros((natms, 3))
        # weights = np.array(weights) / sum_weights
        # self.individual_weights.append(weights)
        # for i in range(n_points):

        #     if self.NACs:
        #         self.impes_coordinate.NAC += weights[i] * NACs[i]
        #     self.impes_coordinate.energy += weights[i] * potentials[i]

        #     self.impes_coordinate.gradient += (
        #         weights[i] * gradients[i] +
        #         potentials[i] * weight_gradients[i] / sum_weights -
        #         potentials[i] * weights[i] * sum_weight_gradients / sum_weights)
        # print('For the comparison 2', self.impes_coordinate.energy, np.linalg.norm(self.impes_coordinate.gradient))
        # print('difference', abs(self.impes_coordinate.energy - check_pot), abs(np.linalg.norm(self.impes_coordinate.gradient) - check_norm_gradient))
        # exit()

    def shepard_interpolation_dsitance_check(self, qm_data_points, core=False):
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


        #n_points = len(qm_data_points)

        # Determine weights, weight gradients,
        # energy values, and energy gradients for interpolation

        distances_and_gradients = []
        min_distance = float('inf')
        self.time_step_reducer = False
        # if len(self.non_core_symmetry_group) > 0:
        #     core = True
        # core = False
        print('Number of qm datapoints is given here:', len(qm_data_points), core)
        self.swapped_indices_map  = None
        if core is True:
            self.swapped_indices_map  = {}
        for i, data_point in enumerate(qm_data_points):
            
            current_gradient = 1.0
            if self.data_points_gradient is not None and 1==2:   
                current_gradient = self.data_points_gradient[i]
            
            if core is False:
                distance, weight_gradient, _ = self.cartesian_distance(data_point, current_gradient)
                print('Distance in IMpes', distance)
                if abs(distance) < min_distance:
                    min_distance = abs(distance)

                distances_and_gradients.append((distance, i, weight_gradient, None))
            else:
                distance, weight_gradient, _, swapped_indices_map = self.cartesian_just_distance_analysis(self.impes_coordinate, data_point, current_gradient)
                if abs(distance) < min_distance:
                    min_distance = abs(distance)
                distances_and_gradients.append((distance, i, weight_gradient, swapped_indices_map))
        
        print('here is the distances', distances_and_gradients)
        # confidence_radii = 
        # if self.data_points_gradient is not None:
        #     # indices_for_mapping = [self.data_points_gradient[i] for distance, i in distances_and_gradients if abs(distance) <= min_distance + self.distance_thrsh]

        #     # print('here are the indices', indices_for_mapping, '\n\n', self.data_points_gradient, 'MIN', min(self.data_points_gradient), '\n\n HESS max', self.data_points_hessian['max'], 'min', self.data_points_hessian['min'])

        #     confidence_radii = [0.3 * (self.data_points_gradient[i] /  min(self.data_points_gradient)) for distance, i in distances_and_gradients if abs(distance) <= min_distance + self.distance_thrsh]

        #     # print('confidence radii \n', confidence_radii)
        
        close_distances = None
        if self.data_points_gradient is not None and 1==2:
            close_distances = [
                (qm_data_points[index], distance, wg, si, self.data_points_gradient[index]) 
                for distance, index, wg, si in distances_and_gradients 
                if abs(distance) <= min_distance + self.distance_thrsh 
                ]
            
        else:
            close_distances = [
                (qm_data_points[index], distance, wg, si, 1.0) 
                for distance, index, wg, si in distances_and_gradients 
                if abs(distance) <= min_distance + self.distance_thrsh 
                ]
        
        # close_distances = [
        # (distance, index, weight_grad) 
        # for distance, index, weight_grad in distances_and_gradients 
        # if abs(distance) <= min_distance + self.distance_thrsh
        # ]

        # num_workers = max(1, int(min(os.cpu_count(), len(close_distances) / 5)))
        # print('here is the number of qm_datapoints which are left', len(close_distances), num_workers)
        # # Use the minimum of available CPU cores or number of data points
        # # profiler.start_timer('calc datapoint properties')
        # with mp.Pool(processes=num_workers) as pool:
        #     results = pool.map(self.calculate_single_data_point_properties, [(close_distance) for close_distance in close_distances])
        # # profiler.stop_timer('calc datapoint properties')
        # self.impes_coordinate.energy = 0
        # gradient_1 = np.zeros((natms, 3))
        # gradient_2 = np.zeros((natms, 3))
        # gradient_3 = 0
        # self.impes_coordinate.NAC = np.zeros((natms, 3))
        # available_potentials = []
        # for result in results:
        #     sum_weights += result['weight']
        #     sum_weight_gradients += result['weight_gradient']
        #     available_potentials.append((result['potential'], result['weight']))
        #     self.impes_coordinate.energy += result['impes_potential']
        #     gradient_1 += result['impes_gradient_1']
        #     gradient_2 += result['impes_gradient_2']
        #     gradient_3 += result['impes_gradient_3']
        
        # for result in results:
        #     self.individual_weights.append(result['weight'] / sum_weights)
        

        # self.impes_coordinate.energy /= sum_weights
        # self.impes_coordinate.gradient = gradient_1 / sum_weights + gradient_2 / sum_weights - gradient_3 / sum_weights * sum_weight_gradients / sum_weights

        # self.important_weights = [(distance, index) for distance, index, _ in close_distances]
        #print('Here are the important weights', self.important_weights)
        for qm_data_point, distance, weight_grad, swapped_indices, current_gradient in close_distances:
            denominator = (
                (distance / self.confidence_radius)**(2 * self.exponent_p) +
                (distance / self.confidence_radius)**(2 * self.exponent_q))
            weight = 1.0 / denominator
            if swapped_indices is not None:
                print('swapped indices map', swapped_indices)
                self.define_impes_coordinate(swapped_indices[0])
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
        self.individual_weights.append(weights)
        for i in range(n_points):

            if self.NACs:
                self.impes_coordinate.NAC += weights[i] * NACs[i]
            self.impes_coordinate.energy += weights[i] * potentials[i]

            self.impes_coordinate.gradient += (
                weights[i] * gradients[i] +
                potentials[i] * weight_gradients[i] / sum_weights -
                potentials[i] * weights[i] * sum_weight_gradients / sum_weights)
        
        print('Energy IMpes', self.impes_coordinate.energy * 627.509474)

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

    def compute_potential(self, data_point):
        """Calculates the potential energy surface at self.impes_coordinate
           based on the energy, gradient and Hessian of data_point.

           :param data_point:
                ImpesCoordinates object.
        """

        energy = data_point.energy
        grad = data_point.internal_gradient
        hessian = data_point.internal_hessian
        eigenvalues, eigenvectors = np.linalg.eigh(hessian)
        if min(eigenvalues) < -1e-2 and np.linalg.norm(grad) > 0.1:
            # The local harmonic approximation might not work and the time_step needs to be adjusted
            self.scaling_time = True
           
        # damped_eigenvalues = self.adjust_eigenvalues(eigenvalues)

        dist_check = (self.impes_coordinate.internal_coordinates_values - data_point.internal_coordinates_values)

        # implement sin for keeping the structure
        # change the code so it can take on any function that is being used for the dihedral
        # TODO: check if the bond angle can be described by a function as well?
        for i, element in enumerate(self.impes_coordinate.z_matrix):
            if len(element) == 4:
                if sum(i in self.non_core_symmetry_group for i in element) == 1 and 1 == 2:
                    # print(element, np.sin(dist[i]), np.sin((3.0/2.0) * dist[i]))
                    #dist[i] = 0.5 * (np.sin(3.0 * dist[i])) * np.tan(1.5 * dist[i])
                    dist_check[i] = 0.5 * (np.sin(3 * dist_check[i] + (2.0 * np.pi)))
                    
                else:
                    #dist[i] = 0.5 * (np.sin(3.0 * dist[i]* dist[i] + (3.0 * np.pi / 2))) + 0.5
                    dist_check[i] = (np.sin(dist_check[i]))

        pes = (energy + np.matmul(dist_check.T, grad) +
               0.5 * np.linalg.multi_dot([dist_check.T, hessian, dist_check]))
        
        return pes

    def adjust_eigenvalues(self, eigenvalues, alpha=5.0):
        adjusted_eigenvalues = []
        for ev in eigenvalues:
            if ev >= -1e2:
                adjusted_eigenvalues.append(ev)
            else:
                scaling_factor = np.exp(alpha * ev)
                adjusted_eigenvalues.append(ev * scaling_factor)
        return np.array(adjusted_eigenvalues)
    
    def compute_gradient(self, data_point, swapped_indices_map=None):
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

        for i, element in enumerate(self.impes_coordinate.z_matrix):
            if len(element) == 4:
                if sum(i in self.non_core_symmetry_group for i in element) == 1 and 1 == 2:
                    dist_check[i] = 0.5 * (np.sin(3.0 * dist_check[i] + (2.0 * np.pi)))
                    dist_check[i] *= (3.0/2.0) * np.cos(3.0 * dist_check[i])
                    grad[i] *= (3.0/2.0) * np.cos((3.0) * dist_check[i])
                else:   
                    dist_check[i] = np.sin(dist_check[i])
                    dist_check[i] *= np.cos(dist_check[i])
                    grad[i] *= np.cos(dist_check[i])

        
        dist_hessian = np.matmul(dist_check.T, hessian)

        gradient_check = (np.matmul(im_b_matrix.T, grad) +
                    np.matmul(im_b_matrix.T, dist_hessian)).reshape(natm, 3)
        
        if swapped_indices_map is not None:
            reference_swap_ind = swapped_indices_map[1]
            swapping_ind = swapped_indices_map[2]
            gradient_check[reference_swap_ind] = gradient_check[swapping_ind]

        return gradient_check

    def database_extracter(self, datafile, molecule):
        
        im_driver = ImpesDriver(molecule.get_z_matrix()) # -> implemented Class in VeloxChem that is capable to perform interpolation calculations for a given molecule and provided z_matrix and database
        im_driver.checkpoint_file_name = datafile
        labels = im_driver.read_labels()
        sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))

        z_mat = molecule.get_z_matrix()
        mol_labels = molecule.get_labels()

        impes_coordinate = ImpesCoordinates(z_mat) # -> implemented Class in VeloxChem that handles all transformations and database changes concerning the interpolation
        data_point_molecules = []
        data_point_energy = []
        data_point_gradient = []
        data_point_hessian = {}

        data_point_hessian['min'] = []
        data_point_hessian['max'] = []
        data_point_hessian['trace'] = []


        for label in sorted_labels:
            impes_coordinate.read_hdf5(datafile, label) # -> read in function from the ImpesDriver object
            current_molecule = Molecule(mol_labels, impes_coordinate.cartesian_coordinates, 'bohr') # -> creates a VeloxChem Molecule object

            eigenvalues, _ = np.linalg.eigh(impes_coordinate.internal_hessian)

            min_eigenvalue = np.min(eigenvalues)
            max_eigenvalue = np.max(eigenvalues)

            trace_hess = np.trace(impes_coordinate.internal_hessian)
            
            data_point_molecules.append(current_molecule)
            data_point_energy.append(float(impes_coordinate.energy))
            data_point_gradient.append(float(np.linalg.norm(impes_coordinate.internal_gradient)))
            data_point_hessian['min'].append(min_eigenvalue)
            data_point_hessian['max'].append(max_eigenvalue)
            data_point_hessian['trace'].append(trace_hess)



        return data_point_molecules, data_point_energy, data_point_gradient, data_point_hessian
    

    def compute_NAC(self, data_point):
        """Calculates part of the cartesian non-adiabatic couplings
           for self.impes_coordinate based on the NACs of data_point.

           :param data_point:
                ImpesCoordinates object.
        """

        natm = data_point.cartesian_coordinates.shape[0]
        nac = data_point.internal_NAC
        print('internal NAC', nac)

        dist = (self.impes_coordinate.internal_coordinates_values -
                data_point.internal_coordinates_values)

        transformation_vec = np.zeros_like(nac)
        for i, element in enumerate(self.impes_coordinate.z_matrix):
            if len(element) == 4:
                dist[i] = np.sin(dist[i])
                dist[i] *= np.cos(dist[i])
                nac[i] *= np.cos(dist[i])
            else:
                transformation_vec[i] += 1
        
        im_b_matrix = self.impes_coordinate.b_matrix
        

        #print('shapes', dist_hessian.shape, dist.shape, hessian.shape)
        NACs = (np.matmul(im_b_matrix.T, nac)).reshape(natm, 3)


        return NACs

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

    def shepard_weight_gradient(self, distance_vector, distance, current_gradient):
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
            (distance / self.confidence_radius * current_gradient)**(2 * self.exponent_p / current_gradient) +
            (distance / self.confidence_radius * current_gradient)**(2 * self.exponent_q / current_gradient))
        derivative_p = (self.exponent_p / current_gradient * distance_vector *
                        distance**(2 * self.exponent_p / current_gradient - 2) /
                        (self.confidence_radius * current_gradient)**(2 * self.exponent_p / current_gradient))
        derivative_q = (self.exponent_q / current_gradient * distance_vector *
                        distance**(2 * self.exponent_q / current_gradient - 2) /
                        (self.confidence_radius * current_gradient)**(2 * self.exponent_q / current_gradient))

        weight_gradient = (-1.0 * (2 * (derivative_p + derivative_q)) *
                           (1.0 / (denominator**2)))

        return weight_gradient
    

    def assign_atoms_hungarian(self, x, y):
        
        cost_matrix = np.linalg.norm(x[:, np.newaxis, :] - y[np.newaxis, :, :], axis=2)
        
        # Apply the Hungarian algorithm
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
        
        return row_ind, col_ind, cost_matrix
    
    
    def cartesian_distance(self, data_point, current_gradient):
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
                distance_vector, distance, current_gradient)
        elif self.interpolation_type == 'simple':
            weight_gradient = self.simple_weight_gradient(
                distance_vector, distance, None)
        else:
            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)

        return distance, weight_gradient, distance_vector_norm
    


    def cartesian_just_distance(self, data_point):
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

        # Calculate the Cartesian distance
        distance = (np.linalg.norm(rotated_coordinates - reference_coordinates))

        return distance


    def get_energy(self):
        """ Returns the potential energy obtained by interpolation.
        """
        return self.impes_coordinate.energy

    def get_gradient(self):
        """ Returns the gradient obtained by interpolation.
        """
        return self.impes_coordinate.gradient


    def analysis_tool(self, number_of_points, z_matrix, given_list=None, compare_each_data_point=False):

        # What I want to test here is a function that can be called in order to asses if an alignment lead to a better agreement result based on the energy
        # for that a database is used and the structures within are evaluated based on the predivtive performance of points with different weighting schemes
        #
        # - call the function with a given database
        # - the function will extract structures with the labels
        # - check the distance of a reference and anoter point j with the two different methods
        # - compare the energy of the points and the interpolation result
        # - estimate the accuracy of the different methods and the best mapping scheme

        if self.checkpoint_file_name is None:
            assert_msg_critical(
            self.checkpoint_file_name is not None, 'ImpesDriver: Please provide a chekpoint file name.')

        if self.non_core_structure is None:
            assert_msg_critical(
            self.checkpoint_file_name is not None, 'ImpesDriver: Please provide non core atoms.')
        
        labels = self.read_labels()
        sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))

        qm_data_points = []
        for label in sorted_labels:
            impes_coordinate = ImpesCoordinates(z_matrix)
            impes_coordinate.read_hdf5(self.checkpoint_file_name, label)

            qm_data_points.append(impes_coordinate)

        
        
        already_used_labels = []
        reference_label_index = 0
        reference_label = random.choice(sorted_labels)
        for i in range(number_of_points):

            while reference_label in already_used_labels and reference_label_index > 0 and i != 0:
                
                reference_label = random.choice(sorted_labels)
            
                reference_label_index = sorted_labels.index(reference_label)        
            
            already_used_labels.append(reference_label)
            reference_label_index = sorted_labels.index(reference_label)  
            if reference_label_index == 0:
                continue
            ref_qm_data_point = qm_data_points[reference_label_index]
            min_distance_core = [np.inf, 0]
            min_distance_all = [np.inf, 0]
            # for i, data_point in enumerate(qm_data_points[:reference_label_index - 1]):
                
            #     distance_all, distance_core = self.cartesian_just_distance_analysis(ref_qm_data_point, data_point)

            #     if distance_all < min_distance_all[0]:
            #         min_distance_all = [distance_all, i]
                    
            #     if distance_core < min_distance_core[0]:
            #         min_distance_core = [distance_core, i]

            print(qm_data_points[:reference_label_index], reference_label_index)
            self.define_impes_coordinate(ref_qm_data_point.cartesian_coordinates)
            if compare_each_data_point:
                for i in range(reference_label_index):
                    self.shepard_interpolation(qm_data_points[:i + 1])
                    energy_all = self.impes_coordinate.energy
                    weights_all = self.important_weights

                    self.define_impes_coordinate(ref_qm_data_point.cartesian_coordinates)
                    self.shepard_interpolation(qm_data_points[:i + 1], core=True)
                    energy_core = self.impes_coordinate.energy

                    print('\n\n', reference_label, '\n', weights_all, '\n', self.important_weights, '\n', energy_all * 627.509474, energy_core * 627.509474, ref_qm_data_point.energy * 627.509474)

            else:

                self.shepard_interpolation(qm_data_points[:reference_label_index])
                energy_all = self.impes_coordinate.energy
                weights_all = self.important_weights
                self.define_impes_coordinate(ref_qm_data_point.cartesian_coordinates)
                self.shepard_interpolation(qm_data_points[:reference_label_index], core=True)
                energy_core = self.impes_coordinate.energy

                print('\n\n', reference_label, '\n', weights_all, '\n', self.important_weights, '\n', energy_all * 627.509474, energy_core * 627.509474, ref_qm_data_point.energy * 627.509474)

    def interpolation_performance_analysis(self, molecule_traj, z_matrix):

        labels = self.read_labels()
        sorted_labels = sorted(labels, key=lambda x: int(x.split('_')[1]))

        qm_data_points = []
        for label in sorted_labels:
            impes_coordinate = ImpesCoordinates(z_matrix, self.molecule.get_labels())
            impes_coordinate.read_hdf5(self.checkpoint_file_name, label)

            qm_data_points.append(impes_coordinate)
        
        self.single_point_interpolation = {}
        self.combination_interpolation = {}
        self.data_point_energies = {}

        for i in range(len(sorted_labels[:7])):
            self.single_point_interpolation[sorted_labels[i]] = ([], [], [], [], None, None)
            self.data_point_energies[sorted_labels[i]] = qm_data_points[i].energy
            if i > 0:
                combination_label = '-'.join(map(str, sorted_labels[:i + 1]))
                print('combination_label', combination_label)
                self.combination_interpolation[combination_label] = ([], [], [], [], None, None)
            for j, mol in enumerate(molecule_traj):

                print('Current datapoint', i, 'Current molecule ', j, '\n')
                self.define_impes_coordinate(mol.get_coordinates_in_bohr())
                self.shepard_interpolation_dsitance_check([qm_data_points[i]])
                energy_all = self.impes_coordinate.energy
                gradient_all = np.linalg.norm(self.impes_coordinate.gradient)
                weights_all = self.important_weights
                self.define_impes_coordinate(mol.get_coordinates_in_bohr())
                self.shepard_interpolation_dsitance_check([qm_data_points[i]], core=True)
                energy_core = self.impes_coordinate.energy
                gradient_core = np.linalg.norm(self.impes_coordinate.gradient)
                
                self.single_point_interpolation[sorted_labels[i]][0].append(energy_all)
                self.single_point_interpolation[sorted_labels[i]][1].append(energy_core)
                self.single_point_interpolation[sorted_labels[i]][2].append(gradient_all)
                self.single_point_interpolation[sorted_labels[i]][3].append(gradient_core)
                current_tuple = self.single_point_interpolation[sorted_labels[i]]
                self.single_point_interpolation[sorted_labels[i]] = (current_tuple[0], current_tuple[1], current_tuple[2], current_tuple[3], weights_all, self.important_weights)

                if i > 0:
                    print('here is the combination', combination_label, len(qm_data_points[:i + 1]))
                    
                    self.define_impes_coordinate(mol.get_coordinates_in_bohr())
                    self.shepard_interpolation_dsitance_check(qm_data_points[:i + 1])
                    
                    energy_all = self.impes_coordinate.energy
                    gradient_all = np.linalg.norm(self.impes_coordinate.gradient)   
                    weights_all = self.important_weights
                    self.define_impes_coordinate(mol.get_coordinates_in_bohr())
                    self.shepard_interpolation_dsitance_check(qm_data_points[:i + 1], core=True)
                    energy_core = self.impes_coordinate.energy
                    gradient_core = np.linalg.norm(self.impes_coordinate.gradient)
                    print('Weights symm', self.important_weights)
                    
                    self.combination_interpolation[combination_label][0].append(energy_all)
                    self.combination_interpolation[combination_label][1].append(energy_core)
                    self.combination_interpolation[combination_label][2].append(gradient_all)
                    self.combination_interpolation[combination_label][3].append(gradient_core)
                    current_tuple = self.combination_interpolation[combination_label]
                    self.combination_interpolation[combination_label] = (current_tuple[0], current_tuple[1], current_tuple[2], current_tuple[3], weights_all, self.important_weights)            



    def cartesian_just_distance_analysis(self, ref_data_point, data_point, current_gradient, data_point_idx=None):
        """Calculates and returns the cartesian distance between
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.

           :param data_point:
                ImpesCoordinates object
        """

        #print('QM_Datapoint internal coordinates \n\n')

        # First, translate the cartesian coordinates to zero
        target_coordinates = data_point.calculate_translation_coordinates()
        reference_coordinates = (
            ref_data_point.calculate_translation_coordinates())

        # remove the non_core atoms (currently just H)
        target_coordinates_core = np.delete(target_coordinates, self.non_core_symmetry_group, axis=0)
        reference_coordinates_core = np.delete(reference_coordinates, self.non_core_symmetry_group, axis=0)

        # Then, determine the rotation matrix which
        # aligns data_point (target_coordinates)
        # to self.impes_coordinate (reference_coordinates)
        
        weightning_dictionary = {'H':1.0, 'C': 1.00, 'F': 1.0, 'N': 1.0, 'O': 1.0, 'S': 1.00}

        mol_atom_labels = self.molecule.get_labels()

        weight_vector = np.array([weightning_dictionary[atom_label] for atom_label in mol_atom_labels ])
    

        rotation_matrix_core = geometric.rotate.get_rot(target_coordinates_core,
                                                   reference_coordinates_core)
        
        # Rotate the data point
        rotated_coordinates_core = np.dot(rotation_matrix_core, target_coordinates.T).T

        reference_group_idx = [idx for idx, _ in enumerate(self.molecule.get_labels()) if idx not in self.non_core_symmetry_group]
        reference_group = np.delete(np.array(reference_coordinates), reference_group_idx, axis=0)
        current_group = np.delete(np.array(rotated_coordinates_core), reference_group_idx, axis=0)

       
        
        row_ind, col_ind, cost_matrix = self.assign_atoms(current_group, reference_group)
        swapped_indices = {}
        swapped_elements_map = {}
        swapped_indices_list = []
        reference_indices_list = []
        indices_that_need_to_be_swapped = []
        element = 0
        current_element_list = []
        for i, j in zip(row_ind, col_ind):
            swapped_indices[self.non_core_symmetry_group[i]] = self.non_core_symmetry_group[j]
            swapped_indices_list.append(self.non_core_symmetry_group[i])
            current_element_list.append(self.non_core_symmetry_group[j])
            reference_indices_list.append(self.non_core_symmetry_group[j])
            indices_that_need_to_be_swapped.append(self.non_core_symmetry_group[j])
            if len(current_element_list) == 3:

                swapped_elements_map[element] = tuple(current_element_list)
                current_element_list = []
                element += 1

        y_assigned = reference_group[col_ind]
        ref_structure_check = reference_coordinates.copy()
        impes_new_coordinates = ref_data_point.cartesian_coordinates.copy()
        impes_assign = impes_new_coordinates[reference_indices_list]
        impes_new_coordinates[swapped_indices_list] = impes_assign
        ref_structure_check[swapped_indices_list] = y_assigned

        self.swapped_indices_map[data_point.point_label] = [impes_new_coordinates, reference_indices_list, swapped_indices_list]
        swapped_indices_map = [impes_new_coordinates, reference_indices_list, swapped_indices_list]
        
        # Calculate the Cartesian distance
        distance_core = (np.linalg.norm(rotated_coordinates_core - ref_structure_check))
        
        #print('Distance norm', distance_core)
        # Calculate the gradient of the interpolation weights
        # (required for energy gradient interpolation)
        distance_vector_core = (ref_structure_check - rotated_coordinates_core)
        distance_vector_core_norm = np.zeros(reference_coordinates_core.shape[0])

        for i in range(len(distance_vector_core_norm)):
            distance_vector_core_norm[i] += np.linalg.norm(distance_vector_core[i])

        if distance_core < 1e-10:
            distance_core = 1e-7
            distance_vector_core[:] = 0
        if self.interpolation_type == 'shepard':
            weight_gradient = self.shepard_weight_gradient(
                distance_vector_core, distance_core, current_gradient)
        elif self.interpolation_type == 'simple':
            weight_gradient = self.simple_weight_gradient(
                distance_vector_core, distance_core, None)
        else:
            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)

        return distance_core, weight_gradient, distance_vector_core_norm, swapped_indices_map


    def assign_atoms(self, x, y):
    
        cost_matrix = np.linalg.norm(x[:, np.newaxis, :] - y[np.newaxis, :, :], axis=2)
        
        # Apply the Hungarian algorithm
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
        
        return row_ind, col_ind, cost_matrix
        


    def invert_positions_through_center(self, positions, center):

        """
        Invert positions through the given center.

        Parameters:
        - positions: numpy.ndarray of shape (N, 3)
        - center: numpy.ndarray of shape (3,)

        Returns:
        - inverted_positions: numpy.ndarray of shape (N, 3)
        """
        return 2 * center - positions

    def find_correspondence_via_inversion(self, x, y, carbon_index, hydrogen_indices):
        """
        Find correspondence between hydrogens in x and y via inversion through the carbon atom.

        Parameters:
        - x: numpy.ndarray of shape (N, 3) for Structure A
        - y: numpy.ndarray of shape (N, 3) for Structure B
        - carbon_index: int, index of the carbon atom
        - hydrogen_indices: list of int, indices of the hydrogen atoms

        Returns:
        - correspondence: dict mapping indices in x to indices in y
        """
        correspondence = {}
        carbon_x = x[carbon_index]

        bond_atom = [x[index] for i, index in enumerate(self.molecule.get_connectivity_matrix()[carbon_index]) if index == 1 and i not in hydrogen_indices]
        
        bond_vector = carbon_x - bond_atom[0]

        reflection_point = self.find_perpendicular_point(carbon_x, bond_vector, x[hydrogen_indices[0]])

        print('bond vector', bond_vector, carbon_index)

        # Invert hydrogens in x through carbon_x
        x_inverted = x.copy()
        if np.linalg.norm(x[hydrogen_indices[0]] - y[hydrogen_indices[0]]) > 1e-4: 
            x_inverted[hydrogen_indices] = self.invert_positions_through_center(x[hydrogen_indices], reflection_point)
        
        else:
            x_inverted[hydrogen_indices] = x[hydrogen_indices]
        # For each inverted hydrogen in x, find the closest hydrogen in y
        new_mask = []
        for h_idx in hydrogen_indices:
            inverted_position = x_inverted[h_idx]
            position = x[h_idx]
            # Compute distances to hydrogens in y
            distances = np.linalg.norm(y[hydrogen_indices] - inverted_position, axis=1)
            distances_org = np.linalg.norm(y[hydrogen_indices] - position, axis=1)
            print('Distances to inverted', distances, distances_org)
            # Find the index of the closest hydrogen
            min_idx = hydrogen_indices[np.argmin(distances)]
            correspondence[h_idx] = min_idx
            new_mask.append(min_idx)

        return correspondence, new_mask, x_inverted


    def find_perpendicular_point(self, line_point, bond_vector, target_point):
        """
        Find the point on the bond line that is perpendicular to the target point.
        
        Parameters:
        line_point (array-like): A point on the line, e.g., the bond start point.
        bond_vector (array-like): The direction vector of the bond (bond vector).
        target_point (array-like): The point to which the perpendicular is drawn.
        
        Returns:
        numpy.ndarray: The point on the line that is perpendicular to the target point.
        """
        # Convert inputs to NumPy arrays
        line_point = np.array(line_point)
        bond_vector = np.array(bond_vector)
        target_point = np.array(target_point)
        
        # Compute the vector from the line point to the target point
        vector_to_target = target_point - line_point
        
        # Compute the projection scalar (t)
        t = np.dot(vector_to_target, bond_vector) / np.dot(bond_vector, bond_vector)
        
        # Find the perpendicular point on the line
        perpendicular_point = line_point + t * bond_vector
        
        return perpendicular_point
