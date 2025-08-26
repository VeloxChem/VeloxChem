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

from tkinter import E
from mpi4py import MPI
import multiprocessing as mp
from collections import Counter
import os
import numpy as np
import math
import random
from scipy.optimize import linear_sum_assignment
from scipy.linalg import cho_factor, cho_solve
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
        - imforcefield_file: File containing the necessary information to construct the interpolation forcefield
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
        self.z_matrix = z_matrix
        self.impes_dict = None
        self.sum_of_weights = None
        self.sum_of_weights_grad = None
        self.print = False
        self.store_weights = False
        self.perform_mapping = True
        self.symmetry_information = None
        self.use_symmetry = True
        self.weights = {}
        self.int_coord_weights = {}
        self.potentials = []
        self.gradients = []
        self.distance_molecule = None

        self.calc_hess = True        
        
        self.skip_first_check = 0
        rng = np.random.default_rng(1)
        self.A = rng.uniform(-1, 1, size=(3*4, 3*4))

        # kpoint file with QM data
        self.imforcefield_file = None
        self.qm_data_points = None
        self.symmetry_dihedral_lists = {}
        self.qm_symmetry_data_points = {}
        self.qm_symmetry_data_points_1 = None
        self.qm_symmetry_data_points_2 = None

        self.bond_rmsd = None
        self.angle_rmsd = None
        self.dihedral_rmsd = None
        self.averaged_int_dist = None
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
        
        self.bond_rmsd = []
        self.angle_rmsd = []
        self.dihedral_rmsd = []
        self.molecule = molecule

        self.distance_molecule = Molecule.from_xyz_string(molecule.get_xyz_string())
        
        # for dihedral in self.symmetry_information[7]:
        #     self.distance_molecule.set_dihedral([dihedral[0] + 1, dihedral[1] + 1, dihedral[2] + 1, dihedral[3] + 1], 0.0, 'degree')
        
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

        def numerical_gradient(f, X, dp, eps=1e-5):
            """
            Compute numerical gradient of scalar function f(X) w.r.t. coordinates X.
            
            Parameters
            ----------
            f : callable
                Function taking X (Nc,3) and returning scalar float.
            X : ndarray
                Coordinate array of shape (Nc,3).
            eps : float
                Finite difference step.
            
            Returns
            -------
            grad : ndarray
                Numerical gradient with same shape as X.
            """
            grad = np.zeros_like(X, dtype=float)
            Nc = X.shape[0]
            for a in range(Nc):
                for mu in range(3):
                    dX = np.zeros_like(X)
                    dX[a, mu] = eps
                    current_geometry_plus = X + dX
                    self.define_impes_coordinate(current_geometry_plus)
                    _, _, denominator_p, _, _, _ = f(dp)
                    f_plus = 1.0 / denominator_p
                    # f_plus, _, _  = f(dp, self.impes_coordinate.internal_coordinates_values)
                    current_geometry_minus = X - dX
                    self.define_impes_coordinate(current_geometry_minus)
                    # f_minus, _, _  = f(dp, self.impes_coordinate.internal_coordinates_values)
                    _, _, denominator_m, _, _, _  = f(dp)
                    f_minus = 1.0 / denominator_m
                    grad[a, mu] = (f_plus - f_minus) / (2 * eps)
            
            return grad
        
        natms = self.impes_coordinate.cartesian_coordinates.shape[0]

        sum_weights = 0.0
        sum_weights_cart = 0.0
        potentials = []
        gradients = []
        hessian_error = []
        self.weights = {}
        beysian_error = []
        weights_cart = []
        averaged_int_dists = []
        weight_gradients_cart = []
        used_labels = []

        sum_weight_gradients_cart = np.zeros((natms - len(self.symmetry_information[4]), 3))

        distances_and_gradients = []
        min_distance = float('inf')
        self.time_step_reducer = False
        

        if not self.use_symmetry and 1==2:
            for i, data_point in enumerate(self.qm_data_points[:]):
                
                distance, denominator, weight_gradient, distance_vector, dihedral_dist = self.cartesian_distance(data_point)

                if abs(distance) < min_distance:
                    min_distance = abs(distance)

                distances_and_gradients.append((distance, dihedral_dist, i, denominator, weight_gradient, distance_vector))
        
        else:
            for i, data_point in enumerate(self.qm_data_points[:]):
                
                distance, dihedral_dist, denominator, weight_gradient, distance_vec, _ = self.cartesian_distance_symmetry(data_point)
                
                # grad = numerical_gradient(self.cartesian_distance_symmetry, self.molecule.get_coordinates_in_bohr(), data_point)
                
                # self.define_impes_coordinate(self.molecule.get_coordinates_in_bohr())
                # print('Gradient difference ', weight_gradient - grad, weight_gradient[0,0], grad[0,0])
              

                if abs(distance) < min_distance:
                    min_distance = abs(distance)

                distances_and_gradients.append((distance, dihedral_dist, i, denominator, weight_gradient, distance_vec))

        # used_labels = []
        close_distances = None
        close_distances = [
            (self.qm_data_points[index], distance, dihedral_dist, denom, wg, distance_vec, index) 
            for distance, dihedral_dist, index, denom, wg, distance_vec in distances_and_gradients 
            if abs(distance) <= min_distance + self.distance_thrsh]
        # distances_from_points = [] 
        # internal_distance_dict = []
        for qm_data_point, distance, dihedral_dist, denominator_cart, weight_grad_cart, distance_vector, label_idx in close_distances:
            
            weight_cart = 1.0 / (denominator_cart)

            # used_labels.append(label_idx)
            sum_weights_cart += weight_cart
            sum_weight_gradients_cart += weight_grad_cart


            potential, gradient, r_i = self.compute_potential(qm_data_point, self.impes_coordinate.internal_coordinates_values)

            # grad = numerical_gradient(self.compute_potential, self.molecule.get_coordinates_in_bohr(), qm_data_point)
                
            # self.define_impes_coordinate(self.molecule.get_coordinates_in_bohr())
            # print('Gradient difference ', gradient - grad)
           
            gradients.append(gradient)
            # bey_int_coord[self.symmetry_information[5]] = qm_data_point.internal_coordinates_values[self.symmetry_information[5]]
            # _, var = qm_data_point.bayes_predict(bey_int_coord)
            averaged_int_dists.append(qm_data_point.internal_coordinates_values)
            # beysian_error.append(var)
            hessian_error.append(r_i)
            potentials.append(potential)
            used_labels.append(label_idx)
            weights_cart.append(weight_cart)
            self.potentials.append(potential)
            self.gradients.append(gradient)
            weight_gradients_cart.append(weight_grad_cart)
            # distances_from_points.append(distance)

        n_points = len(close_distances)
        # if self.print:
        #     for i in range(n_points - 1):
        #         print('\nint distane \n', i, distances_from_points[i], np.linalg.norm(internal_distance_dict[i]), '\n',
        #               i + 1, distances_from_points[i+ 1], np.linalg.norm(internal_distance_dict[i+ 1]), '\n\n', internal_distance_dict[i], '\n\n', internal_distance_dict[i + 1], '\n\n', internal_distance_dict[i] - internal_distance_dict[i + 1])
        
                
        
        # --- initialise accumulators -------------------------------------------------
        self.impes_coordinate.energy    = 0.0
        self.impes_coordinate.gradient  = np.zeros((natms, 3))
        self.impes_coordinate.NAC       = np.zeros((natms, 3))       # if you need it

        # --- 1.  raw (unnormalised) weights and their gradients ----------------------
        w_i          = np.array(weights_cart, dtype=np.float64)        # ← rename
        grad_w_i     = np.array(weight_gradients_cart, dtype=np.float64)   # shape (n_pts, natms, 3)

        S            = w_i.sum()                         # Σ wᵢ
        sum_grad_w   = grad_w_i.sum(axis=0)              # Σ ∇wᵢ      shape (natms, 3)
        # for lbl, wi in zip(used_labels, w_i):
        #     self.weights[lbl] = wi
        self.sum_of_weights = S
        self.sum_of_weights_grad = sum_grad_w

        # --- 2.  normalised weights and their gradients ------------------------------
        W_i          = w_i / S
        grad_W_i     = (grad_w_i * S - w_i[:, None, None] * sum_grad_w) / S**2

        # --- 3.  accumulate energy and gradient --------------------------------------
        potentials   = np.array(potentials, dtype=np.float64)        # Uᵢ
        gradients    = np.array(gradients,  dtype=np.float64)        # ∇Uᵢ  shape (n_pts, natms, 3)

        self.impes_coordinate.energy   = np.dot(W_i, potentials)     # Σ Wᵢ Uᵢ

        # ∇U = Σ Wᵢ ∇Uᵢ  +  Σ Uᵢ ∇Wᵢ
        if len(self.symmetry_information[3]) != natms:
            self.impes_coordinate.gradient = (np.tensordot(W_i, gradients, axes=1))
            # self.impes_coordinate.gradient[self.symmetry_information[4]] += (gradients[:, self.symmetry_information[4], :].sum(axis=0))
            # Add contributions only to the selected rows
            
            self.impes_coordinate.gradient[self.symmetry_information[3]] += np.tensordot(potentials, grad_W_i, axes=1)
        else:

            self.impes_coordinate.gradient = (np.tensordot(W_i, gradients, axes=1) + np.tensordot(potentials, grad_W_i, axes=1))


        # --- 4.  book-keeping (optional) ---------------------------------------------
        for lbl, Wi in zip(used_labels, W_i):
            self.weights[lbl] = Wi

        # self.sum_of_weights      = W_i.sum()          # if you really need it later
        self.averaged_int_dist   = np.tensordot(W_i, averaged_int_dists, axes=1)
        
            # self.impes_coordinate.energy += potentials[i]

            # self.impes_coordinate.gradient += ( gradients[i])

            # sigma_sq += weights_cart[i] * hessian_error[i]/0.5
        
        # print('Bey Error', np.sqrt(sigma_sq) * hartree_in_kcalpermol() )
        
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

    def te_weight_general(self, theta, periodicity=1):
        """
        General smooth gate function centered at 0 and 2π.
        
        Parameters
        ----------
        theta : float
            Input angle (in radians).
        periodicity : int or float
            Controls the angular width of the smoothing zone: width = π / periodicity.

        Returns
        -------
        weight : float
            A value between 0 and 1, smoothly falling from 1 to 0 over the defined window.
        """
        def _smoother(t):
            """Quintic smoother-step kernel (C² on 0 ≤ t ≤ 1)."""
            return 6*t**5 - 15*t**4 + 10*t**3

        width = np.pi / periodicity
        x = np.mod(theta, 2*np.pi)

        if x <= width:
            t = 1.0 - x / width
            return _smoother(t)
        elif x >= 2*np.pi - width:
            t = 1.0 - (2*np.pi - x) / width
            return _smoother(t)
        else:
            return 0.0
        
    def te_weight_gradient_general(self, theta, b_matrix_col, periodicity=1):
        """
        Gradient of the general te_weight function, scaled by the B-matrix column.
        
        Parameters
        ----------
        theta : float
            Input angle (in radians).
        b_matrix_col : np.ndarray
            The relevant column of the B-matrix (Cartesian direction).
        periodicity : int or float
            Controls the angular width of the smoothing zone: width = π / periodicity.

        Returns
        -------
        gradient : np.ndarray
            Weighted gradient vector (same shape as b_matrix_col).
        """
        def _smoother_prime(t):
            """Derivative of the quintic smoother-step kernel."""
            return 30*t**4 - 60*t**3 + 30*t**2

        width = np.pi / periodicity
        x = np.mod(theta, 2*np.pi)

        wprime = 0.0  # default if outside the taper zones

        if x <= width:
            t = 1.0 - x / width
            dt_dθ = -1.0 / width
            wprime = _smoother_prime(t) * dt_dθ

        elif x >= 2*np.pi - width:
            t = 1.0 - (2*np.pi - x) / width
            dt_dθ = 1.0 / width
            wprime = _smoother_prime(t) * dt_dθ

        return wprime * b_matrix_col
    
        


    def compute_potential(self, data_point, org_int_coords):
        """Calculates the potential energy surface at self.impes_coordinate
           based on the energy, gradient and Hessian of data_point.

           :param data_point:
                InterpolationDatapoint object.
        """
        def pinv_cut(mat, rcond=1e-8):
            """Moore–Penrose pseudo‑inverse with controllable cut‑off."""
            U, s, Vt = np.linalg.svd(mat, full_matrices=False)
            large = s > rcond
            Sinv = np.zeros_like(s)
            Sinv[large] = 1.0 / s[large]
            return (Vt.T * Sinv) @ U.T

        def periodic_delta(phi, phi_ref):
            """Return smallest signed angular difference in radians (-π, π]."""
            return ((phi - phi_ref + np.pi) % (2*np.pi)) - np.pi
        # implement sin for keeping the structure
        # change the code so it can take on any function that is being used for the dihedral
        # TODO: check if the bond angle can be described by a function as well?
        # print('Mapping', mapping_list[1])
        # cos_weight_func = 0.5 * (1+np.cos())
        
        
        pes = 0.0
        gradient = None
        natm = data_point.cartesian_coordinates.shape[0]
        # print(len(self.qm_symmetry_data_points))
        dp_label = data_point.point_label 
        if len(self.qm_symmetry_data_points[dp_label]) > 1 and self.use_symmetry:

            symmetry_data_points = self.qm_symmetry_data_points[dp_label]
  
            symmetry_weights = []
            potentials = []
            pot_gradients = []
            symmetry_weight_gradients = []
            hessian_error = []
            
            for i, symmetry_data_point in enumerate(symmetry_data_points):
                # if i == 0:
                #     symmetry_data_point.confidence_radius = 1.5
                energy = symmetry_data_point.energy
                masks = symmetry_data_point.mapping_masks

                for ww, mask in enumerate(masks):
                    
                    
                    grad = symmetry_data_point.internal_gradient.copy()
                    grad[masks[0]] = grad[mask]

                    hessian = symmetry_data_point.internal_hessian.copy()
                    hessian[np.ix_(masks[0], masks[0])] = hessian[np.ix_(mask, mask)]


                    dist_org = (org_int_coords.copy() - symmetry_data_point.internal_coordinates_values[mask])

                    # print('mask', (org_int_coords.copy(), symmetry_data_point.internal_coordinates_values[mask]))

                    dist_check = (org_int_coords.copy() - symmetry_data_point.internal_coordinates_values[mask])
                    dist_correlation = (org_int_coords.copy() - symmetry_data_point.internal_coordinates_values[mask])
                    
                    sum_sym_dihedral = 0.0
                    sum_sym_dihedral_prime = np.zeros_like(data_point.cartesian_coordinates.reshape(-1))
           
                    for i, element in enumerate(self.impes_coordinate.z_matrix[self.symmetry_information[-1][1]:], start=self.symmetry_information[-1][1]): 

                        if tuple(sorted(element)) in self.symmetry_information[7][3]:                                       

                            dist_check[i] = np.sin(dist_org[i])
                            sum_sym_dihedral += self.te_weight_general(dist_org[i], 3)
                            sum_sym_dihedral_prime += self.te_weight_gradient_general(dist_org[i], self.impes_coordinate.b_matrix[i,:], 3)
                            dist_correlation[i] = 0.0



                        elif tuple(sorted(element)) in self.symmetry_information[7][2]:
                            dist_check[i] = np.sin(dist_org[i])
                            sum_sym_dihedral += self.te_weight_general(dist_org[i], 2)
                            sum_sym_dihedral_prime += self.te_weight_gradient_general(dist_org[i], self.impes_coordinate.b_matrix[i,:], 2)      
                            dist_correlation[i] = 0.0


                        else:
                            
                            dist_check[i] = np.sin(dist_org[i])  
                            dist_correlation[i] = np.sin(dist_org[i])

                        # print('sum sym dihedral', sum_sym_dihedral)
                    
                    # print('\n\n')
                    self.bond_rmsd.append(np.sqrt(np.mean(np.sum((dist_org[:self.symmetry_information[-1][0]])**2))))
                    self.angle_rmsd.append(np.sqrt(np.mean(np.sum(dist_org[self.symmetry_information[-1][0]:self.symmetry_information[-1][1]]**2))))
                    self.dihedral_rmsd.append(np.sqrt(np.mean(np.sum(dist_correlation[self.symmetry_information[-1][1]:]**2))))

                    if sum_sym_dihedral == 0.0:
                        continue
                    
                    
                    symmetry_weights.append(sum_sym_dihedral**10)

                    symmetry_weight_gradients.append(10 * sum_sym_dihedral**9 * sum_sym_dihedral_prime.reshape(natm, 3))

                    pes = (energy + np.matmul(dist_check.T, grad) +
                        0.5 * np.linalg.multi_dot([dist_check.T, hessian, dist_check]))

                    # hessian_error.append(np.sqrt(abs(np.linalg.multi_dot([dist_check.T, pinv_cut(hessian), dist_check]))))

                    potentials.append(pes)
                    dist_hessian = np.matmul(dist_check.T, hessian)
                    for i, element in enumerate(self.impes_coordinate.z_matrix[self.symmetry_information[-1][1]:], start=self.symmetry_information[-1][1]):
                        
                        grad[i] *= np.cos(dist_org[i])
                        dist_hessian[i] *= np.cos(dist_org[i])
                    
                    pes_prime = (np.matmul(self.impes_coordinate.b_matrix.T, (grad + dist_hessian))).reshape(natm, 3)
                    pot_gradients.append(pes_prime)

                
            # ---------------------------------------------------------------------------
            # 1. raw (unnormalised) Shepard weights and their gradients
            # ---------------------------------------------------------------------------
            w_i       = np.asarray(symmetry_weights,          dtype=np.float64)        # shape (m,)
            grad_w_i  = np.asarray(symmetry_weight_gradients, dtype=np.float64)        # shape (m, natm, 3)

            S         = w_i.sum()                              # Σ wᵢ      (scalar)
            sum_grad_w = grad_w_i.sum(axis=0)                  # Σ ∇wᵢ     (natm, 3)

            # ---------------------------------------------------------------------------
            # 2. normalised weights and their gradients
            # ---------------------------------------------------------------------------
            W_i       = w_i / S                                                    # (m,)

            grad_W_i  = (grad_w_i * S - w_i[:, None, None] * sum_grad_w) / S**2    # (m, natm, 3)
            # print('internal wieghts', W_i)
            # ---------------------------------------------------------------------------
            # 3. energy and gradient of the interpolated PES
            # ---------------------------------------------------------------------------
            U_i       = np.asarray(potentials,        dtype=np.float64)            # energies     (m,)
            grad_U_i  = np.asarray(pot_gradients,     dtype=np.float64)            # gradients    (m, natm, 3)

            pes       = np.dot(W_i, U_i)                      # Σ Wᵢ Uᵢ

            grad_pes  = ( np.tensordot(W_i,      grad_U_i, axes=1) +   # Σ Wᵢ ∇Uᵢ
                        np.tensordot(U_i,   grad_W_i, axes=1) )       # Σ Uᵢ ∇Wᵢ
            # grad_pes shape: (natm, 3)

            # ---------------------------------------------------------------------------
            # 4. book-keeping exactly as before
            # ---------------------------------------------------------------------------
            for j, Wi in enumerate(W_i):
                self.bond_rmsd[j]     *= Wi
                self.angle_rmsd[j]    *= Wi
                self.dihedral_rmsd[j] *= Wi

            gradient = grad_pes           # this is what the caller expects
            return pes, gradient, hessian_error
        
        else:
            hessian_error = 0.0
            energy = data_point.energy
            grad = data_point.internal_gradient.copy()
            hessian = data_point.internal_hessian.copy()
            dist_org = (org_int_coords.copy() - data_point.internal_coordinates_values)
            dist_check = (org_int_coords.copy() - data_point.internal_coordinates_values)

            # for bond_idx, element_bond in enumerate(self.impes_coordinate.z_matrix[:self.symmetry_information[-1][0]]):
            #     dist_check[bond_idx] = (1.0 / org_int_coords[bond_idx]) - (1.0 / data_point.internal_coordinates_values[bond_idx])
            for i, element in enumerate(self.impes_coordinate.z_matrix[self.symmetry_information[-1][1]:], start=self.symmetry_information[-1][1]): 
                                       
                dist_check[i] = np.sin(dist_org[i])


            self.bond_rmsd.append(np.sqrt(np.mean(np.sum((dist_org[:self.symmetry_information[-1][0]])**2))))
            self.angle_rmsd.append(np.sqrt(np.mean(np.sum(dist_org[self.symmetry_information[-1][0]:self.symmetry_information[-1][1]]**2))))
            self.dihedral_rmsd.append(np.sqrt(np.mean(np.sum(dist_check[self.symmetry_information[-1][1]:]**2))))
            pes = (energy + np.matmul(dist_check.T, grad) +
                        0.5 * np.linalg.multi_dot([dist_check.T, hessian, dist_check]))
            
            dist_hessian = np.matmul(dist_check.T, hessian)
            
            # for bond_idx, element_bond in enumerate(self.impes_coordinate.z_matrix[:self.symmetry_information[-1][0]]):
            #     grad[bond_idx] *= -1.0 / (org_int_coords[bond_idx])**2
            #     dist_hessian[bond_idx] *= -1.0 / (org_int_coords[bond_idx])**2

            for i, element in enumerate(self.impes_coordinate.z_matrix[self.symmetry_information[-1][1]:], start=self.symmetry_information[-1][1]):
                
                
                grad[i] *= np.cos(dist_org[i])
                dist_hessian[i] *= np.cos(dist_org[i])
            
            # masses = self.molecule.get_masses().copy()
            # masses_cart = np.repeat(masses, 3)
            # inv_sqrt_masses = 1.0 / np.sqrt(masses_cart)
            # self.impes_coordinate.inv_sqrt_masses = inv_sqrt_masses
            
            pes_prime = (np.matmul(self.impes_coordinate.b_matrix.T, (grad + dist_hessian))).reshape(natm, 3)

            return pes, pes_prime, (grad + dist_hessian)


        


    
    def te_weight_gradient(self, theta, b_matrix_col):
        
        def _smoother_prime(t):
            """Derivative of the smoother-step kernel."""
            return 30*t**4 - 60*t**3 + 30*t**2
        
        x = np.mod(theta, 2*np.pi)          # fold into [0, 2π)

        # initialise to zero so the "else" branch is automatic
        wprime = 0.0

        if x <= np.pi/3:                    # first taper
            t       = 1.0 - x / (np.pi/3)
            dt_dθ   = -1.0 / (np.pi/3)      #  d(1-x/L)/dθ   with  L = π/3
            wprime  = _smoother_prime(t) * dt_dθ

        elif x >= 5*np.pi/3:                # second taper
            t       = 1.0 - (2*np.pi - x) / (np.pi/3)
            dt_dθ   =  1.0 / (np.pi/3)
            wprime  = _smoother_prime(t) * dt_dθ

        # multiply by the B-matrix column to obtain the Cartesian gradient
        return wprime * b_matrix_col        # same shape as b_matrix_col
            
    def compute_gradient(self, data_point, mapping_list, current_internal_coordinates_values, org_int_coords, current_b_matrix, org_b_matrix):
        """Calculates part of the cartesian gradient
           for self.impes_coordinate based on the gradient
           and Hessian of data_point.
           J ( g_a + delta_a H_a )
           eq. (4) JCTC 12, 5235-5246

           :param data_point:
                ImpesCoordinates object.
        """

        natm = data_point.cartesian_coordinates.shape[0]
        gradient = None
        if len(self.qm_symmetry_data_points) > 0:
            # if self.qm_symmetry_data_points_1 is not None:
            #     dp_label = data_point.point_label  
            #     symmetry_data_point = self.qm_symmetry_data_points[dp_label]
            #     symmetry_data_point_1 = self.qm_symmetry_data_points_1[data_point.point_label]
            #     symmetry_data_point_2 = self.qm_symmetry_data_points_2[symmetry_data_point.point_label]
                
            #     energy = data_point.energy
            #     grad = data_point.internal_gradient.copy()
            #     grad_11 = data_point.internal_gradient.copy()
            #     grad_12 = symmetry_data_point_1[0].internal_gradient.copy()
            #     grad_13 = symmetry_data_point_1[1].internal_gradient.copy()


            #     hessian_11 = data_point.internal_hessian
            #     hessian_12 = symmetry_data_point_1[0].internal_hessian
            #     hessian_13 = symmetry_data_point_1[1].internal_hessian


            #     dist_check = (org_int_coords - data_point.internal_coordinates_values)
            #     org_dist_check = (org_int_coords - data_point.internal_coordinates_values)
                
                
            #     sym_energy = symmetry_data_point.energy
            #     sym_grad = symmetry_data_point.internal_gradient.copy()
            #     grad_21 = symmetry_data_point.internal_gradient.copy()
            #     grad_22 = symmetry_data_point_2[0].internal_gradient.copy()
            #     grad_23 = symmetry_data_point_2[1].internal_gradient.copy()
                
            #     hessian_21 = symmetry_data_point.internal_hessian
            #     hessian_22 = symmetry_data_point_2[0].internal_hessian
            #     hessian_23 = symmetry_data_point_2[1].internal_hessian
            #     sym_dist_check = (org_int_coords - symmetry_data_point.internal_coordinates_values)
                
                
            #     # org_sym_dist_check = (org_int_coords - symmetry_data_point.internal_coordinates_values)
                

            #     dist_check_11 = (org_int_coords - data_point.internal_coordinates_values)
            #     dist_check_12 = (org_int_coords - symmetry_data_point_1[0].internal_coordinates_values)
            #     dist_check_13 = (org_int_coords - symmetry_data_point_1[1].internal_coordinates_values)

            #     dist_check_21 = (org_int_coords - symmetry_data_point.internal_coordinates_values)
            #     dist_check_22 = (org_int_coords - symmetry_data_point_2[0].internal_coordinates_values)
            #     dist_check_23 = (org_int_coords - symmetry_data_point_2[1].internal_coordinates_values)

            #     sum_sym_dihedral_11_prime = np.zeros_like(data_point.internal_gradient)
            #     sum_sym_dihedral_12_prime = sum_sym_dihedral_11_prime.copy()
            #     sum_sym_dihedral_13_prime = sum_sym_dihedral_11_prime.copy()
            #     sum_sym_dihedral_21_prime = sum_sym_dihedral_11_prime.copy()
            #     sum_sym_dihedral_22_prime = sum_sym_dihedral_11_prime.copy()
            #     sum_sym_dihedral_23_prime = sum_sym_dihedral_11_prime.copy()

            #     sum_sym_dihedral_11 = 0.0
            #     sum_sym_dihedral_12 = 0.0
            #     sum_sym_dihedral_13 = 0.0
            #     sum_sym_dihedral_21 = 0.0
            #     sum_sym_dihedral_22 = 0.0
            #     sum_sym_dihedral_23 = 0.0


                

                
            #     org_molecule = Molecule(self.molecule.get_labels(), data_point.cartesian_coordinates, 'bohr')
                
            #     rot_molecule1 = Molecule(self.molecule.get_labels(), data_point.cartesian_coordinates, 'bohr')
            #     rot_molecule1.set_dihedral([5,2,1,7], org_molecule.get_dihedral([5,2,1,7], 'radian') - 2.0 * np.pi/3.0, 'radian')

                
            #     rot_molecule_11 = Molecule(self.molecule.get_labels(), rot_molecule1.get_coordinates_in_bohr(), 'bohr')

            #     rot_molecule2 = Molecule(self.molecule.get_labels(), data_point.cartesian_coordinates, 'bohr')

            #     rot_molecule2.set_dihedral([5,2,1,7], org_molecule.get_dihedral([5,2,1,7], 'radian') - 4.0 * np.pi/3.0, 'radian')
            #     rot_molecule_12 = Molecule(self.molecule.get_labels(), rot_molecule2.get_coordinates_in_bohr(), 'bohr')
                

            #     mapping_list_11, mapping_dict_11, assigned = self.perform_symmetry_assignment(self.symmetry_information[0], self.symmetry_information[1][0], 
            #                                                         org_molecule.get_coordinates_in_bohr()[self.symmetry_information[1][0]], rot_molecule_11.get_coordinates_in_bohr()[self.symmetry_information[1][0]])
                
            #     mapping_list_12, mapping_dict_12, assigned = self.perform_symmetry_assignment(self.symmetry_information[0], self.symmetry_information[1][0], 
            #                                                         org_molecule.get_coordinates_in_bohr()[self.symmetry_information[1][0]], rot_molecule_12.get_coordinates_in_bohr()[self.symmetry_information[1][0]])
            #     print('mapping_list', mapping_list_11, mapping_list_12)

            #     z_matrix_dict = {tuple(sorted(element)): i 
            #         for i, element in enumerate(self.impes_coordinate.z_matrix)}
                
            #     reorded_int_coords = [org_int_coords]
            #     b_matrices = [org_b_matrix]
            #     mapping_lists = [mapping_list_11, mapping_list_12]
            #     mapping_dict = [mapping_dict_11, mapping_dict_12]
            #     reordered_int_coord_values = org_int_coords.copy()
            #     for ord in range(2):
            #         org_mask = [i for i in range(len(self.impes_coordinate.z_matrix))]
            #         org_b_matrix_cp = org_b_matrix.copy()
            #         mask = []
            #         for i, element in enumerate(self.impes_coordinate.z_matrix):
            #             # Otherwise, reorder the element
            #             if len(element) == 4 or 1 == 1:
            #                 reordered_element = [mapping_dict[ord].get(x, x) for x in element]
            #                 key = tuple(sorted(reordered_element))

            #                 z_mat_index = z_matrix_dict.get(key)
            #                 mask.append(z_mat_index)
            #                 reordered_int_coord_values[i] = (float(org_int_coords[z_mat_index]))

            #             else:
            #                 mask.append(i)
            #                 reordered_int_coord_values[i] = (float(org_int_coords[i]))

            #             # print(z_mat_index, key, element)
            #         org_b_matrix_cp[org_mask] = org_b_matrix[mask]
            #         org_b_matrix_cp = org_b_matrix_cp.reshape(len(self.impes_coordinate.z_matrix), len(self.molecule.get_labels()), 3)
            #         org_b_matrix_cp = org_b_matrix_cp[:, mapping_lists[ord], :]
            #         org_b_matrix_cp = org_b_matrix_cp.reshape(len(self.impes_coordinate.z_matrix), len(self.molecule.get_labels()) * 3)

            #         b_matrices.append(org_b_matrix_cp)
            #         reorded_int_coords.append(reordered_int_coord_values)
                    
            #         # print('B-matrix', self.symmetry_information[0], b_matrices[0] - b_matrices[1], '\n org B-Matrix \n')

            #     for i, element in enumerate(self.impes_coordinate.z_matrix): 
            #         if len(element) == 4:                                       
                                
            #                 # dist_check[i] = 1.0/2.0*np.sin(3.0*dist_check[i])
                            
            #                 # dist_check[i] = 0.5*(1+np.sin(3.0*dist_check[i] - np.pi/2.0))
            #                 # dist_check[i] = 3.0/2.0*np.sin(dist_check[i])**2
            #                 # org scheme
            #                 dist_check_11[i] = np.sin(dist_check[i])
            #                 grad_11[i] *= np.cos(dist_check[i])
            #                 sum_sym_dihedral_11 += self.te_weight(dist_check[i])
            #                 sum_sym_dihedral_11_prime += self.te_weight_gradient(dist_check[i], b_matrices[0][i,:])

            #                 dist_check_12[i] = np.sin(dist_check_12[i])
            #                 grad_12[i] *= np.cos(org_int_coords[i] - symmetry_data_point_1[0].internal_coordinates_values[i])
            #                 sum_sym_dihedral_12 += self.te_weight(1.0*(org_int_coords[i] - symmetry_data_point_1[0].internal_coordinates_values[i]))
            #                 sum_sym_dihedral_12_prime += self.te_weight_gradient(1.0*(org_int_coords[i] - symmetry_data_point_1[0].internal_coordinates_values[i]), b_matrices[0][i,:])

            #                 dist_check_13[i] = np.sin(dist_check_13[i])
            #                 grad_13[i] *= np.cos(org_int_coords[i] - symmetry_data_point_1[1].internal_coordinates_values[i])
            #                 sum_sym_dihedral_13 += self.te_weight(1.0*(org_int_coords[i] - symmetry_data_point_1[1].internal_coordinates_values[i]))
            #                 sum_sym_dihedral_13_prime += self.te_weight_gradient(1.0*(org_int_coords[i] - symmetry_data_point_1[1].internal_coordinates_values[i]), b_matrices[0][i,:])
                            

            #                 dist_check_21[i] = np.sin(sym_dist_check[i])
            #                 grad_21[i] *= np.cos(sym_dist_check[i])
            #                 sum_sym_dihedral_21 += self.te_weight(sym_dist_check[i])
            #                 sum_sym_dihedral_21_prime += self.te_weight_gradient(sym_dist_check[i], b_matrices[0][i,:])

            #                 dist_check_22[i] = np.sin(dist_check_22[i])
            #                 grad_22[i] *= np.cos(org_int_coords[i] - symmetry_data_point_2[0].internal_coordinates_values[i])
            #                 sum_sym_dihedral_22 += self.te_weight(1.0*(org_int_coords[i] - symmetry_data_point_2[0].internal_coordinates_values[i]))
            #                 sum_sym_dihedral_22_prime += self.te_weight_gradient(1.0*(org_int_coords[i] - symmetry_data_point_2[0].internal_coordinates_values[i]), b_matrices[0][i,:])


            #                 dist_check_23[i] = np.sin(dist_check_23[i])
            #                 grad_23[i] *= np.cos(org_int_coords[i] - symmetry_data_point_2[1].internal_coordinates_values[i])
            #                 sum_sym_dihedral_23 += self.te_weight(1.0*(org_int_coords[i] - symmetry_data_point_2[1].internal_coordinates_values[i]))
            #                 sum_sym_dihedral_23_prime += self.te_weight_gradient(1.0*(org_int_coords[i] - symmetry_data_point_2[1].internal_coordinates_values[i]), b_matrices[0][i,:])
                            
            #                 # dist_check_12[i] = np.sin(org_int_coords[i] - (data_point.internal_coordinates_values[i] + 2.0*np.pi/3.0))
            #                 # grad_12[i] *= np.cos(org_int_coords[i] - (data_point.internal_coordinates_values[i] + 2.0*np.pi/3.0))
            #                 # sum_sym_dihedral_12 += self.te_weight(-1.0*(org_int_coords[i] - (data_point.internal_coordinates_values[i] + 2.0*np.pi/3.0)))
            #                 # sum_sym_dihedral_12_prime += self.te_weight_gradient(-1.0*(org_int_coords[i] - (data_point.internal_coordinates_values[i] + 2.0*np.pi/3.0)), b_matrices[0][i,:])

            #                 # dist_check_13[i] = np.sin(org_int_coords[i] - (data_point.internal_coordinates_values[i] + 4.0*np.pi/3.0))
            #                 # grad_13[i] *= np.cos(org_int_coords[i] - (data_point.internal_coordinates_values[i] + 4.0*np.pi/3.0))
            #                 # sum_sym_dihedral_13 += self.te_weight(-1.0*(org_int_coords[i] - (data_point.internal_coordinates_values[i] + 4.0*np.pi/3.0)))
            #                 # sum_sym_dihedral_13_prime += self.te_weight_gradient(-1.0*(org_int_coords[i] - (data_point.internal_coordinates_values[i] + 4.0*np.pi/3.0)), b_matrices[0][i,:])
                            

            #                 # dist_check_21[i] = np.sin(sym_dist_check[i])
            #                 # grad_21[i] *= np.cos(sym_dist_check[i])
            #                 # sum_sym_dihedral_21 += self.te_weight(-1.0*(sym_dist_check[i]))
            #                 # sum_sym_dihedral_21_prime += self.te_weight_gradient(-1.0*(org_int_coords[i] - symmetry_data_point.internal_coordinates_values[i]), b_matrices[0][i,:])

            #                 # dist_check_22[i] = np.sin(org_int_coords[i] - (symmetry_data_point.internal_coordinates_values[i] + 2.0*np.pi/3.0))
            #                 # grad_22[i] *= np.cos(org_int_coords[i] - (symmetry_data_point.internal_coordinates_values[i] + 2.0*np.pi/3.0))
            #                 # sum_sym_dihedral_22 += self.te_weight(-1.0*(org_int_coords[i] - (symmetry_data_point.internal_coordinates_values[i] + 2.0*np.pi/3.0)))
            #                 # sum_sym_dihedral_22_prime += self.te_weight_gradient(-1.0*(org_int_coords[i] - (symmetry_data_point.internal_coordinates_values[i] + 2.0*np.pi/3.0)), b_matrices[0][i,:])


            #                 # dist_check_23[i] = np.sin(org_int_coords[i] - (symmetry_data_point.internal_coordinates_values[i] + 4.0*np.pi/3.0))
            #                 # grad_23[i] *= np.cos(org_int_coords[i] - (symmetry_data_point.internal_coordinates_values[i] + 4.0*np.pi/3.0))
            #                 # sum_sym_dihedral_23 += self.te_weight(-1.0*(org_int_coords[i] - (symmetry_data_point.internal_coordinates_values[i] + 4.0*np.pi/3.0)))
            #                 # sum_sym_dihedral_23_prime += self.te_weight_gradient(-1.0*(org_int_coords[i] - (symmetry_data_point.internal_coordinates_values[i] + 4.0*np.pi/3.0)), b_matrices[0][i,:])

                
                
            #     v11, v12, v13, v21, v22, v23 = sum_sym_dihedral_11, sum_sym_dihedral_12, sum_sym_dihedral_13, sum_sym_dihedral_21, sum_sym_dihedral_22, sum_sym_dihedral_23
                
            #     print('unnormalized weights', v11, v12, v13, v21, v22, v23)

            #     v11_prime, v12_prime, v13_prime, v21_prime, v22_prime, v23_prime = sum_sym_dihedral_11_prime.reshape(natm, 3), sum_sym_dihedral_12_prime.reshape(natm, 3), sum_sym_dihedral_13_prime.reshape(natm, 3), sum_sym_dihedral_21_prime.reshape(natm, 3), sum_sym_dihedral_22_prime.reshape(natm, 3), sum_sym_dihedral_23_prime.reshape(natm, 3)

            #     total_weight_sum = v11 + v12 + v13 + v21 + v22 + v23
            #     total_weight_grad_sum = v11_prime + v12_prime + v13_prime + v21_prime + v22_prime + v23_prime



            #     w11 = v11 / total_weight_sum
            #     w12 = v12 / total_weight_sum    
            #     w13 = v13 / total_weight_sum
            #     w21 = v21 / total_weight_sum
            #     w22 = v22 / total_weight_sum
            #     w23 = v23 / total_weight_sum

            #     w11_prime = (v11_prime * total_weight_sum - v11 * total_weight_grad_sum) / total_weight_sum**2
            #     w12_prime = (v12_prime * total_weight_sum - v12 * total_weight_grad_sum) / total_weight_sum**2
            #     w13_prime = (v13_prime * total_weight_sum - v13 * total_weight_grad_sum) / total_weight_sum**2
            #     w21_prime = (v21_prime * total_weight_sum - v21 * total_weight_grad_sum) / total_weight_sum**2
            #     w22_prime = (v22_prime * total_weight_sum - v22 * total_weight_grad_sum) / total_weight_sum**2
            #     w23_prime = (v23_prime * total_weight_sum - v23 * total_weight_grad_sum) / total_weight_sum**2

            #     sum_w_prime = (w11_prime + w12_prime + w13_prime + w21_prime + w22_prime + w23_prime).sum(axis=(0,1))
            #     print('weights', w11, w12, w13, w21, w22, w23)
            #     print('sum_w_prime', sum_w_prime, total_weight_sum)
            #     # exit()

            #     pes_11 = (energy + np.matmul(dist_check_11.T, grad_11) +
            #         0.5 * np.linalg.multi_dot([dist_check_11.T, hessian_11, dist_check_11]))
            #     pes_12 = (symmetry_data_point_1[0].energy + np.matmul(dist_check_12.T, grad_12) +
            #         0.5 * np.linalg.multi_dot([dist_check_12.T, hessian_12, dist_check_12]))
            #     pes_13 = (symmetry_data_point_1[1].energy + np.matmul(dist_check_13.T, grad_13) +
            #         0.5 * np.linalg.multi_dot([dist_check_13.T, hessian_13, dist_check_13]))
            #     pes_21 = (sym_energy + np.matmul(dist_check_21.T, grad_21) +
            #         0.5 * np.linalg.multi_dot([dist_check_21.T, hessian_21, dist_check_21]))
            #     pes_22 = (symmetry_data_point_2[0].energy + np.matmul(dist_check_22.T, grad_22) +
            #         0.5 * np.linalg.multi_dot([dist_check_22.T, hessian_22, dist_check_22]))
            #     pes_23 = (symmetry_data_point_2[1].energy + np.matmul(dist_check_23.T, grad_23) +
            #         0.5 * np.linalg.multi_dot([dist_check_23.T, hessian_23, dist_check_23]))
                
            #     dist_hessian_11 = np.matmul(dist_check_11.T, hessian_11)
            #     dist_hessian_12 = np.matmul(dist_check_12.T, hessian_12)
            #     dist_hessian_13 = np.matmul(dist_check_13.T, hessian_13)
            #     dist_hessian_21 = np.matmul(dist_check_21.T, hessian_21)
            #     dist_hessian_22 = np.matmul(dist_check_22.T, hessian_22)
            #     dist_hessian_23 = np.matmul(dist_check_23.T, hessian_23)
                
            #     for i, element in enumerate(self.impes_coordinate.z_matrix):
            #         if len(element) == 4:

            #             dist_hessian_11[i] *= np.cos(dist_check[i])
            #             dist_hessian_12[i] *= np.cos(org_int_coords[i] - symmetry_data_point_1[0].internal_coordinates_values[i])
            #             dist_hessian_13[i] *= np.cos(org_int_coords[i] - symmetry_data_point_1[1].internal_coordinates_values[i])
            #             dist_hessian_21[i] *= np.cos(sym_dist_check[i])
            #             dist_hessian_22[i] *= np.cos(org_int_coords[i] - symmetry_data_point_2[0].internal_coordinates_values[i])
            #             dist_hessian_23[i] *= np.cos(org_int_coords[i] - symmetry_data_point_2[1].internal_coordinates_values[i])
                
            #     pes_11_prime = (np.matmul(b_matrices[0].T, (grad_11 + dist_hessian_11))).reshape(natm, 3)
            #     pes_12_prime = (np.matmul(b_matrices[0].T, (grad_12 + dist_hessian_12))).reshape(natm, 3)
            #     pes_13_prime = (np.matmul(b_matrices[0].T, (grad_13 + dist_hessian_13))).reshape(natm, 3)
            #     pes_21_prime = (np.matmul(b_matrices[0].T, (grad_21 + dist_hessian_21))).reshape(natm, 3)
            #     pes_22_prime = (np.matmul(b_matrices[0].T, (grad_22 + dist_hessian_22))).reshape(natm, 3)
            #     pes_23_prime = (np.matmul(b_matrices[0].T, (grad_23 + dist_hessian_23))).reshape(natm, 3)
                
            #     if mapping_list is not None:
            #         print('mapping_list, org', mapping_list_11, self.symmetry_information[0])
            #         print('old_map_list', mapping_list[1], mapping_list[0])
                

            #      # self.impes_coordinate.gradient += (
            # #     weights_cart[i] * gradients[i] +
            # #     potentials[i] * weight_gradients_cart[i] / sum_weights_cart -
            # #     potentials[i] * weights_cart[i] * sum_weight_gradients_cart / sum_weights_cart)
            #     pes_weight_11_prime = w11_prime * pes_11 + w11 * pes_11_prime
            #     pes_weight_12_prime = w12_prime * pes_12 + w12 * pes_12_prime
            #     pes_weight_13_prime = w13_prime * pes_13 + w13 * pes_13_prime
            #     pes_weight_21_prime = w21_prime * pes_21 + w21 * pes_21_prime
            #     pes_weight_22_prime = w22_prime * pes_22 + w22 * pes_22_prime
            #     pes_weight_23_prime = w23_prime * pes_23 + w23 * pes_23_prime

            #     # pes_weight_12_prime[mapping_list_11] = pes_weight_12_prime[self.symmetry_information[0]]
            #     # pes_weight_13_prime[mapping_list_12] = pes_weight_13_prime[self.symmetry_information[0]]
            #     # pes_weight_22_prime[mapping_list_11] = pes_weight_22_prime[self.symmetry_information[0]]
            #     # pes_weight_23_prime[mapping_list_12] = pes_weight_23_prime[self.symmetry_information[0]]
                
            #     gradient = pes_weight_11_prime + pes_weight_12_prime + pes_weight_13_prime + pes_weight_21_prime + pes_weight_22_prime + pes_weight_23_prime

                
            #     # print('gradient', gradient)





            # gradient_2p = 0.0

            # dp_label = data_point.point_label  
            # symmetry_data_point = self.qm_symmetry_data_points[dp_label][1]
            
            # energy = data_point.energy

            # org_molecule = Molecule(self.molecule.get_labels(), data_point.cartesian_coordinates, 'bohr')
            
            # rot_molecule1 = Molecule(self.molecule.get_labels(), data_point.cartesian_coordinates.copy(), 'bohr')
            # rot_molecule1.set_dihedral([7,2,1,5], org_molecule.get_dihedral([7,2,1,5], 'radian') - 2.0 * np.pi/3.0, 'radian')

            # print('org_mol_1', org_molecule.get_dihedral([4,2,1,3], 'radian'))
            # rot_molecule_11 = Molecule(self.molecule.get_labels(), rot_molecule1.get_coordinates_in_bohr(), 'bohr')

            # rot_molecule2 = Molecule(self.molecule.get_labels(), data_point.cartesian_coordinates.copy(), 'bohr')

            # rot_molecule2.set_dihedral([7,2,1,5], org_molecule.get_dihedral([7,2,1,5], 'radian') - 4.0 * np.pi/3.0, 'radian')
            # print('org_mol_2', org_molecule.get_dihedral([4,2,1,3], 'radian'))
            
            # rot_molecule_12 = Molecule(self.molecule.get_labels(), rot_molecule2.get_coordinates_in_bohr(), 'bohr')
            
            # print(org_molecule.get_xyz_string())
            # print('dihedrals', np.sin(rot_molecule_11.get_dihedral([4,2,1,3], 'radian')), np.sin(rot_molecule_12.get_dihedral([4,2,1,3], 'radian')))

            # mapping_list_11, mapping_dict_11, assigned = self.perform_symmetry_assignment(self.symmetry_information[0], self.symmetry_information[1][0], 
            #                                                       org_molecule.get_coordinates_in_bohr()[self.symmetry_information[1][0]], rot_molecule_11.get_coordinates_in_bohr()[self.symmetry_information[1][0]])
            
            # mapping_list_12, mapping_dict_12, assigned = self.perform_symmetry_assignment(self.symmetry_information[0], self.symmetry_information[1][0], 
            #                                                       org_molecule.get_coordinates_in_bohr()[self.symmetry_information[1][0]], rot_molecule_12.get_coordinates_in_bohr()[self.symmetry_information[1][0]])
            # print('mapping_list', mapping_list_11, mapping_list_12)

            # z_matrix_dict = {tuple(sorted(element)): i 
            #       for i, element in enumerate(self.impes_coordinate.z_matrix)}
            
            # mapping_lists = [mapping_list_11, mapping_list_12]
            # mapping_dict = [mapping_dict_11, mapping_dict_12]
            
            # org_mask = [i for i in range(len(self.impes_coordinate.z_matrix))]
            # org_b_matrix_cp = org_b_matrix.copy()
            # reorded_int_coords_1 = [data_point.internal_coordinates_values.copy()]
            # reorded_int_coords_2 = [symmetry_data_point.internal_coordinates_values.copy()]
            # b_matrices = [org_b_matrix]
            # masks = [org_mask]
            # for ord in range(2):
            #     mask = []
            #     reorded_int_coord_1 = np.zeros_like(data_point.internal_coordinates_values)
            #     reorded_int_coord_2 = np.zeros_like(data_point.internal_coordinates_values)
            #     for i, element in enumerate(self.impes_coordinate.z_matrix):
            #         # Otherwise, reorder the element
            #         if len(element) == 4 or 1 == 1:
            #             reordered_element = [mapping_dict[ord].get(x, x) for x in element]
            #             key = tuple(sorted(reordered_element))

            #             z_mat_index = z_matrix_dict.get(key)
            #             mask.append(z_mat_index)
            #             reorded_int_coord_1[i] = (float(reorded_int_coords_1[0][z_mat_index]))
            #             reorded_int_coord_2[i] = (float(reorded_int_coords_2[0][z_mat_index]))

            #         else:
            #             mask.append(i)
            #             reordered_int_coord_values[i] = (float(org_int_coords[i]))

            #         # print(z_mat_index, key, element)
            #     org_b_matrix_cp[org_mask] = org_b_matrix[mask]
            #     org_b_matrix_cp = org_b_matrix_cp.reshape(len(self.impes_coordinate.z_matrix), len(self.molecule.get_labels()), 3)
            #     org_b_matrix_cp = org_b_matrix_cp[:, mapping_lists[ord], :]
            #     org_b_matrix_cp = org_b_matrix_cp.reshape(len(self.impes_coordinate.z_matrix), len(self.molecule.get_labels()) * 3)

            #     b_matrices.append(org_b_matrix_cp)
            #     reorded_int_coords_1.append(reorded_int_coord_1)
            #     reorded_int_coords_2.append(reorded_int_coord_2)

                
            #     masks.append(mask)



            # grad = data_point.internal_gradient.copy()
            # grad_11 = data_point.internal_gradient.copy()
            # grad_12 = data_point.internal_gradient.copy()
            # grad_13 = data_point.internal_gradient.copy()

            # hessian_11 = data_point.internal_hessian.copy()
            # hessian_12 = data_point.internal_hessian.copy()
            # hessian_13 = data_point.internal_hessian.copy()

            # dist_check = (org_int_coords - data_point.internal_coordinates_values)
            
            # # symmetry point section
            # sym_energy = symmetry_data_point.energy

            # grad_21 = symmetry_data_point.internal_gradient.copy()
            # grad_22 = symmetry_data_point.internal_gradient.copy()
            # grad_23 = symmetry_data_point.internal_gradient.copy()
            
            # hessian_21 = symmetry_data_point.internal_hessian.copy()
            # hessian_22 = symmetry_data_point.internal_hessian.copy()
            # hessian_23 = symmetry_data_point.internal_hessian.copy()



            # sym_dist_check = (org_int_coords - symmetry_data_point.internal_coordinates_values)
            


            # dist_check_11 = (org_int_coords - reorded_int_coords_1[0])
            # dist_check_12 = (org_int_coords - reorded_int_coords_1[1])
            # dist_check_13 = (org_int_coords - reorded_int_coords_1[2])

            # dist_check_21 = (org_int_coords - reorded_int_coords_2[0])
            # dist_check_22 = (org_int_coords - reorded_int_coords_2[1])
            # dist_check_23 = (org_int_coords - reorded_int_coords_2[2])

            # dist_org_11 = (org_int_coords - reorded_int_coords_1[0])
            # dist_org_12 = (org_int_coords - reorded_int_coords_1[1])
            # dist_org_13 = (org_int_coords - reorded_int_coords_1[2])

            # dist_org_21 = (org_int_coords - reorded_int_coords_2[0])
            # dist_org_22 = (org_int_coords - reorded_int_coords_2[1])
            # dist_org_23 = (org_int_coords - reorded_int_coords_2[2])

            # grad_12[masks[0]] = grad_12[masks[1]]
            # grad_13[masks[0]] = grad_13[masks[2]]
            # grad_22[masks[0]] = grad_22[masks[1]]
            # grad_23[masks[0]] = grad_23[masks[2]]

            
            # hessian_12[masks[0], masks[0]] = hessian_12[masks[1], masks[1]]
            # hessian_13[masks[0], masks[0]] = hessian_13[masks[2], masks[2]]
            # hessian_22[masks[0], masks[0]] = hessian_22[masks[1], masks[1]]
            # hessian_23[masks[0], masks[0]] = hessian_23[masks[2], masks[2]]



            # sum_sym_dihedral_11_prime = np.zeros_like(data_point.internal_gradient)
            # sum_sym_dihedral_12_prime = np.zeros_like(data_point.internal_gradient)
            # sum_sym_dihedral_13_prime = np.zeros_like(data_point.internal_gradient)
            # sum_sym_dihedral_21_prime = np.zeros_like(data_point.internal_gradient)
            # sum_sym_dihedral_22_prime = np.zeros_like(data_point.internal_gradient)
            # sum_sym_dihedral_23_prime = np.zeros_like(data_point.internal_gradient)

            # sum_sym_dihedral_11 = 0.0
            # sum_sym_dihedral_12 = 0.0
            # sum_sym_dihedral_13 = 0.0
            # sum_sym_dihedral_21 = 0.0
            # sum_sym_dihedral_22 = 0.0
            # sum_sym_dihedral_23 = 0.0

            # for i, element in enumerate(self.impes_coordinate.z_matrix): 
            #     if len(element) == 4:   

            #             dist_check_11[i] = np.sin(dist_org_11[i])
            #             grad_11[i] *= np.cos(dist_org_11[i])
            #             sum_sym_dihedral_11 += self.te_weight(dist_org_11[i])
            #             sum_sym_dihedral_11_prime += self.te_weight_gradient(dist_org_11[i], b_matrices[0][i,:])                                    
                        
            #             dist_check_12[i] = np.sin(dist_org_12[i])
            #             grad_12[i] *= np.cos(dist_org_12[i])
            #             sum_sym_dihedral_12 += self.te_weight(1.0*(dist_org_12[i]))
            #             sum_sym_dihedral_12_prime += self.te_weight_gradient(1.0*(dist_org_12[i]), b_matrices[0][i,:])

            #             dist_check_13[i] = np.sin(dist_org_13[i])
            #             grad_13[i] *= np.cos(dist_org_13[i])
            #             sum_sym_dihedral_13 += self.te_weight(1.0*(dist_org_13[i]))
            #             sum_sym_dihedral_13_prime += self.te_weight_gradient(1.0*(dist_org_13[i]), b_matrices[0][i,:])

            #             dist_check_21[i] = np.sin(dist_org_21[i])
            #             grad_21[i] *= np.cos(dist_org_21[i])
            #             sum_sym_dihedral_21 += self.te_weight(1.0*(dist_org_21[i]))
            #             sum_sym_dihedral_21_prime += self.te_weight_gradient(1.0*(dist_org_21[i]), b_matrices[0][i,:])

            #             dist_check_22[i] = np.sin(dist_org_22[i])
            #             grad_22[i] *= np.cos(dist_org_22[i])
            #             sum_sym_dihedral_22 += self.te_weight(1.0*(dist_org_22[i]))
            #             sum_sym_dihedral_22_prime += self.te_weight_gradient(1.0*(dist_org_22[i]), b_matrices[0][i,:])


            #             dist_check_23[i] = np.sin(dist_org_23[i])
            #             grad_23[i] *= np.cos(dist_org_23[i])
            #             sum_sym_dihedral_23 += self.te_weight(1.0*(dist_org_23[i]))
            #             sum_sym_dihedral_23_prime += self.te_weight_gradient(1.0*(dist_org_23[i]), b_matrices[0][i,:])

            #             print('distances 2p', dist_check_11[i], dist_check_12[i], dist_check_13[i], dist_check_21[i], dist_check_22[i], dist_check_23[i])
            #             print('sum_sym_dihedral_11', sum_sym_dihedral_11, sum_sym_dihedral_12, sum_sym_dihedral_13, sum_sym_dihedral_21, sum_sym_dihedral_22, sum_sym_dihedral_23)
            
            
            # v11, v12, v13, v21, v22, v23 = sum_sym_dihedral_11, sum_sym_dihedral_12, sum_sym_dihedral_13, sum_sym_dihedral_21, sum_sym_dihedral_22, sum_sym_dihedral_23
            
            # print('unnormalized weights', v11, v12, v13, v21, v22, v23)

            # v11_prime, v12_prime, v13_prime, v21_prime, v22_prime, v23_prime = sum_sym_dihedral_11_prime.reshape(natm, 3), sum_sym_dihedral_12_prime.reshape(natm, 3), sum_sym_dihedral_13_prime.reshape(natm, 3), sum_sym_dihedral_21_prime.reshape(natm, 3), sum_sym_dihedral_22_prime.reshape(natm, 3), sum_sym_dihedral_23_prime.reshape(natm, 3)

            # total_weight_sum = v11 + v12 + v13 + v21 + v22 + v23
            # total_weight_grad_sum = v11_prime + v12_prime + v13_prime + v21_prime + v22_prime + v23_prime



            # w11 =  v11 / total_weight_sum
            # w12 =  v12 / total_weight_sum    
            # w13 =  v13 / total_weight_sum
            # w21 =  v21 / total_weight_sum
            # w22 =  v22 / total_weight_sum
            # w23 =  v23 / total_weight_sum

            # w11_prime = (v11_prime * total_weight_sum - v11 * total_weight_grad_sum) / total_weight_sum**2
            # w12_prime = (v12_prime * total_weight_sum - v12 * total_weight_grad_sum) / total_weight_sum**2
            # w13_prime = (v13_prime * total_weight_sum - v13 * total_weight_grad_sum) / total_weight_sum**2
            # w21_prime = (v21_prime * total_weight_sum - v21 * total_weight_grad_sum) / total_weight_sum**2
            # w22_prime = (v22_prime * total_weight_sum - v22 * total_weight_grad_sum) / total_weight_sum**2
            # w23_prime = (v23_prime * total_weight_sum - v23 * total_weight_grad_sum) / total_weight_sum**2

            # sum_w_prime = (w11_prime + w12_prime + w13_prime + w21_prime + w22_prime + w23_prime).sum(axis=(0,1))
            # print('weights', w11, w12, w13, w21, w22, w23)
            # print('sum_w_prime', sum_w_prime, total_weight_sum)
            # # exit()

            # pes_11 = (energy + np.matmul(dist_check_11.T, grad_11) +
            #     0.5 * np.linalg.multi_dot([dist_check_11.T, hessian_11, dist_check_11]))       
            # pes_12 = (energy + np.matmul(dist_check_12.T, grad_12) +
            #     0.5 * np.linalg.multi_dot([dist_check_12.T, hessian_12, dist_check_12]))
            # pes_13 = (energy + np.matmul(dist_check_13.T, grad_13) +
            #     0.5 * np.linalg.multi_dot([dist_check_13.T, hessian_13, dist_check_13]))
            
            # pes_21 = (sym_energy + np.matmul(dist_check_21.T, grad_21) +
            #     0.5 * np.linalg.multi_dot([dist_check_21.T, hessian_21, dist_check_21]))
            # pes_22 = (sym_energy + np.matmul(dist_check_22.T, grad_22) +
            #     0.5 * np.linalg.multi_dot([dist_check_22.T, hessian_22, dist_check_22]))
            # pes_23 = (sym_energy + np.matmul(dist_check_23.T, grad_23) +
            #     0.5 * np.linalg.multi_dot([dist_check_23.T, hessian_23, dist_check_23]))
            
            # dist_hessian_11 = np.matmul(dist_check_11.T, hessian_11)
            # dist_hessian_12 = np.matmul(dist_check_12.T, hessian_12)
            # dist_hessian_13 = np.matmul(dist_check_13.T, hessian_13)
            # dist_hessian_21 = np.matmul(dist_check_21.T, hessian_21)
            # dist_hessian_22 = np.matmul(dist_check_22.T, hessian_22)
            # dist_hessian_23 = np.matmul(dist_check_23.T, hessian_23)
            
            # for i, element in enumerate(self.impes_coordinate.z_matrix):
            #     if len(element) == 4:

            #         dist_hessian_11[i] *= np.cos(dist_org_11[i])
            #         dist_hessian_12[i] *= np.cos(dist_org_12[i])
            #         dist_hessian_13[i] *= np.cos(dist_org_13[i])
            #         dist_hessian_21[i] *= np.cos(dist_org_21[i])
            #         dist_hessian_22[i] *= np.cos(dist_org_22[i])
            #         dist_hessian_23[i] *= np.cos(dist_org_23[i])
            
            # pes_11_prime = (np.matmul(b_matrices[0].T, (grad_11 + dist_hessian_11))).reshape(natm, 3)
            # pes_12_prime = (np.matmul(b_matrices[0].T, (grad_12 + dist_hessian_12))).reshape(natm, 3)
            # pes_13_prime = (np.matmul(b_matrices[0].T, (grad_13 + dist_hessian_13))).reshape(natm, 3)
            # pes_21_prime = (np.matmul(b_matrices[0].T, (grad_21 + dist_hessian_21))).reshape(natm, 3)
            # pes_22_prime = (np.matmul(b_matrices[0].T, (grad_22 + dist_hessian_22))).reshape(natm, 3)
            # pes_23_prime = (np.matmul(b_matrices[0].T, (grad_23 + dist_hessian_23))).reshape(natm, 3)
            
            
            # pes_weight_11_prime = w11_prime * pes_11 + w11 * pes_11_prime
            # pes_weight_12_prime = w12_prime * pes_12 + w12 * pes_12_prime
            # pes_weight_13_prime = w13_prime * pes_13 + w13 * pes_13_prime
            # pes_weight_21_prime = w21_prime * pes_21 + w21 * pes_21_prime
            # pes_weight_22_prime = w22_prime * pes_22 + w22 * pes_22_prime
            # pes_weight_23_prime = w23_prime * pes_23 + w23 * pes_23_prime

            # print('PESPES GRADIENT SECTION', w11 * pes_11 + w12 * pes_12 + w13 * pes_13 + w21 * pes_21 + w22 * pes_22 + w23 * pes_23)
            # gradient = pes_weight_11_prime + pes_weight_12_prime + pes_weight_13_prime + pes_weight_21_prime + pes_weight_22_prime + pes_weight_23_prime
            # print('Gradient_1', gradient)
            # # gradient[mapping_list[0]] = gradient[mapping_list[1]]
            # # print('gradient', gradient - gradient_2p)
        
            dp_label = data_point.point_label 

                
            symmetry_data_points = self.qm_symmetry_data_points[dp_label]
    
            symmetry_weights = []
            potentials = []
            pot_gradients = []
            symmetry_weight_gradients = []

            for symmetry_data_point in symmetry_data_points:
                    energy = symmetry_data_point.energy
                    masks = symmetry_data_point.mapping_masks

                    for mask in masks:
                        
                        grad = symmetry_data_point.internal_gradient.copy()
                        grad[masks[0]] = grad[mask]

                        hessian = symmetry_data_point.internal_hessian.copy()
                        hessian[masks[0], masks[0]] = hessian[mask, mask]

                        dist_org = (org_int_coords - symmetry_data_point.internal_coordinates_values[mask])
                        dist_check = (org_int_coords - symmetry_data_point.internal_coordinates_values[mask])
                        
                        sum_sym_dihedral_prime = np.zeros_like(symmetry_data_point.internal_gradient)
                        sum_sym_dihedral = 0.0

                        for i, element in enumerate(self.impes_coordinate.z_matrix): 
                            if len(element) == 4:   

                                    dist_check[i] = np.sin(dist_org[i])
                                    grad[i] *= np.cos(dist_org[i])
                                    sum_sym_dihedral += self.te_weight(dist_org[i])
                                    sum_sym_dihedral_prime += self.te_weight_gradient(dist_org[i], org_b_matrix[i,:])

                        symmetry_weights.append(sum_sym_dihedral)
                        symmetry_weight_gradients.append(sum_sym_dihedral_prime.reshape(natm, 3))

                        pes = (energy + np.matmul(dist_check.T, grad) +
                            0.5 * np.linalg.multi_dot([dist_check.T, hessian, dist_check]))
                        
                        potentials.append(pes)

                        dist_hessian = np.matmul(dist_check.T, hessian)

                        for i, element in enumerate(self.impes_coordinate.z_matrix):
                            if len(element) == 4:

                                dist_hessian[i] *= np.cos(dist_org[i])
                        
                        pes_prime = (np.matmul(org_b_matrix.T, (grad + dist_hessian))).reshape(natm, 3)

                        pot_gradients.append(pes_prime)
        

            total_weight_sum = sum(symmetry_weights)
            total_weight_grad_sum = sum(symmetry_weight_gradients)

            weights = np.array(symmetry_weights) / total_weight_sum
            potentials = np.array(potentials)
            potentials_prime = np.array(pot_gradients)
            weights_prime = np.array(symmetry_weight_gradients)
            
            gradient = np.zeros_like(potentials_prime[0])


            for i in range(len(weights)):
                gradient += potentials_prime[i] * weights[i] + potentials[i] * weights_prime[i] / total_weight_sum - potentials[i] * weights[i] * total_weight_grad_sum / total_weight_sum    


            # exit()
        
        
        grad = data_point.internal_gradient.copy()
        org_grad = data_point.internal_gradient.copy()
        hessian = data_point.internal_hessian
            
    
        dist_check = (current_internal_coordinates_values - data_point.internal_coordinates_values)
        dist_org = (current_internal_coordinates_values - data_point.internal_coordinates_values)
        gradient_1 = None
        
        for i, element in enumerate(self.impes_coordinate.z_matrix):
            if len(element) == 4:
                if mapping_list is not None and 1==2:
    
                    dist_check[i] = 1/3*np.sin(3.0*dist_check[i])
                    grad[i] *= 3.0/3.0 * np.cos(3.0*dist_org[i])
                    # old code:
                    # dist_check[i] = np.sin(dist_check[i])
                    # grad[i] *= np.cos(dist_org[i])
                else:
                    
                    # dist_check[i] = 1.0/2.0*np.sin(3.0*dist_check[i])
                    # grad[i] *= 3.0/2.0 * np.cos(3.0*dist_org[i])
                    # dist_check[i] = 0.5*(1+np.sin(3.0*dist_check[i] - np.pi/2.0))
                    # grad[i] *= 3.0/2.0 * np.cos(-3.0*dist_org[i] + np.pi/2.0)
                    dist_check[i] = np.sin(dist_check[i])
                    grad[i] *= np.cos(dist_org[i])
                    
        
        dist_hessian = np.matmul(dist_check.T, hessian)
        for i, element in enumerate(self.impes_coordinate.z_matrix):
            if len(element) == 4:
                if mapping_list is not None and 1 == 2:
                    dist_hessian[i] *= 3.0/3.0 * np.cos(3.0 * dist_org[i])
                    # org scheme
                    # dist_hessian[i] *= np.cos(dist_org[i])
                else:
                    # dist_hessian[i] *= 3.0/3.0 * np.cos(3.0 * dist_org[i])
                    # dist_hessian[i] *= 3.0/2.0 * np.cos(3.0 * dist_org[i])
                    # dist_hessian[i] *= 3.0/2.0 * np.cos(-3.0*dist_org[i] + np.pi/2.0)
                    # print('dist_hessian', dist_hessian[i])
                    # org scheme
                    dist_hessian[i] *= np.cos(dist_org[i])
        internal_gradient_hess_cos = grad + dist_hessian
        gradient_1 = (np.matmul(current_b_matrix.T, internal_gradient_hess_cos)).reshape(natm, 3)
        # print('gradient', gradient)
        # if mapping_list is not None:
        #     gradient_1[mapping_list[1]] = gradient_1[mapping_list[0]]
        
        # print('gradient 2', gradient_1)
            # exit()
        if gradient is None:
            gradient = gradient_1
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
    
    def trust_radius_weight_gradient_hessian(self, datapoint):
        
        confidence_radius = datapoint.confidence_radius

        distance, _, _, _, _, _ = self.cartesian_distance_symmetry(datapoint)
        denominator = (
                (distance / confidence_radius)**(2 * self.exponent_p) +
                (distance / confidence_radius)**(2 * self.exponent_q))
        trust_radius_weight_gradient = -1.0 * ((( -2.0 * self.exponent_p * ((distance / confidence_radius)**(2 * self.exponent_p)) / confidence_radius) - 
                                             (2.0 * self.exponent_q * ((distance / confidence_radius)**(2 * self.exponent_q) / confidence_radius))) / denominator**2)
        return trust_radius_weight_gradient

    def trust_radius_weight_gradient(self, datapoint):
        
        confidence_radius = datapoint.confidence_radius

        distance, _, _, _, _, _ = self.cartesian_distance_symmetry(datapoint)
        denominator = (
                (distance / confidence_radius)**(2 * self.exponent_p) +
                (distance / confidence_radius)**(2 * self.exponent_q))
        trust_radius_weight_gradient = -1.0 * ((( -2.0 * self.exponent_p * ((distance / confidence_radius)**(2 * self.exponent_p)) / confidence_radius) - 
                                             (2.0 * self.exponent_q * ((distance / confidence_radius)**(2 * self.exponent_q) / confidence_radius))) / denominator**2)
        return trust_radius_weight_gradient

    def trust_radius_weight_gradient_gradient(self, datapoint):
        
        confidence_radius = datapoint.confidence_radius
        distance, _, _, _, distance_vector, _ = self.cartesian_distance_symmetry(datapoint)
        denominator = (
                (distance / confidence_radius)**(2 * self.exponent_p) +
                (distance / confidence_radius)**(2 * self.exponent_q))
        
        
        # trust_radius_weight_gradient_gradient_1_nominator = (((-2.0 * self.exponent_p**2 * distance_vector * distance**(2.0 * self.exponent_p - 2.0)) 
        #                                         /(confidence_radius**(2.0 * self.exponent_p) * confidence_radius))
        #                                         - (2.0 * self.exponent_q**2.0 * distance_vector * distance**(2.0 * self.exponent_q - 2.0)) 
        #                                         /(confidence_radius**(2.0 * self.exponent_q) * confidence_radius)) 

        # trust_radius_weight_gradient_gradient_2_nominator = ( 4.0 * (((self.exponent_p * distance_vector * distance**(2.0 * self.exponent_p - 2.0)) 
        #                                         /(confidence_radius**(2.0 * self.exponent_p)))
        #                                         + ( self.exponent_q * distance_vector * distance**(2.0 * self.exponent_q - 2.0)) 
        #                                         /(confidence_radius**(2.0 * self.exponent_q))) * (((-2.0 * (distance/confidence_radius)**(2.0 * self.exponent_p) * self.exponent_p) / confidence_radius) 
        #                                                                                        - ((2.0 * (distance/confidence_radius)**(2.0 * self.exponent_q) * self.exponent_q) / confidence_radius))) 

        
        trust_radius_weight_gradient_gradient_1_nominator = (((-4.0 * self.exponent_p**2 * distance_vector * distance**(2.0 * self.exponent_p - 2.0)) 
                                                /(confidence_radius**(2.0 * self.exponent_p) * confidence_radius))
                                                - (4.0 * self.exponent_q**2.0 * distance_vector * distance**(2.0 * self.exponent_q - 2.0)) 
                                                /(confidence_radius**(2.0 * self.exponent_q) * confidence_radius)) 

        trust_radius_weight_gradient_gradient_2_nominator = ( 2.0 * (((2.0 * self.exponent_p * distance_vector * distance**(2.0 * self.exponent_p - 2.0)) /(confidence_radius**(2.0 * self.exponent_p)))
                                                + (2.0 * self.exponent_q * distance_vector * distance**(2.0 * self.exponent_q - 2.0)) / (confidence_radius**(2.0 * self.exponent_q))) 
                                                * (((-2.0 * (distance/confidence_radius)**(2.0 * self.exponent_p) * self.exponent_p) / confidence_radius) 
                                                - ((2.0 * (distance/confidence_radius)**(2.0 * self.exponent_q) * self.exponent_q) / confidence_radius))) 

        
        trust_radius_weight_gradient_gradient = -1.0 * trust_radius_weight_gradient_gradient_1_nominator / denominator**2 + trust_radius_weight_gradient_gradient_2_nominator / denominator**3


        return trust_radius_weight_gradient_gradient

    def trust_radius_weight_gradient_gradient_hessian(self, datapoint):
        
        confidence_radius = datapoint.confidence_radius
        distance, _, _, _, distance_vector, grad_s = self.cartesian_distance_symmetry(datapoint)
   
        denominator = (
                (distance / confidence_radius)**(2 * self.exponent_p) +
                (distance / confidence_radius)**(2 * self.exponent_q))
        
        # trust_radius_weight_gradient_gradient_1_nominator = (((-2.0 * self.exponent_p**2 * distance_vector * distance**(2.0 * self.exponent_p - 2.0)) 
        #                                         /(confidence_radius**(2.0 * self.exponent_p) * confidence_radius))
        #                                         - (2.0 * self.exponent_q**2.0 * distance_vector * distance**(2.0 * self.exponent_q - 2.0)) 
        #                                         /(confidence_radius**(2.0 * self.exponent_q) * confidence_radius)) 

        # trust_radius_weight_gradient_gradient_2_nominator = ( 4.0 * (((self.exponent_p * distance_vector * distance**(2.0 * self.exponent_p - 2.0)) 
        #                                         /(confidence_radius**(2.0 * self.exponent_p)))
        #                                         + ( self.exponent_q * distance_vector * distance**(2.0 * self.exponent_q - 2.0)) 
        #                                         /(confidence_radius**(2.0 * self.exponent_q))) * (((-2.0 * (distance/confidence_radius)**(2.0 * self.exponent_p) * self.exponent_p) / confidence_radius) 
        #                                                                                        - ((2.0 * (distance/confidence_radius)**(2.0 * self.exponent_q) * self.exponent_q) / confidence_radius))) 

        # trust_radius_weight_gradient_gradient = -2.0 * trust_radius_weight_gradient_gradient_1_nominator / denominator**2 + trust_radius_weight_gradient_gradient_2_nominator / denominator**3
        

        trust_radius_weight_gradient_gradient_nominator_1_1 = 2.0 * ((2.0 * self.exponent_p * (distance / confidence_radius)**(2 * self.exponent_p - 1) / confidence_radius) + (2.0 * self.exponent_q * (distance / confidence_radius)**(2 * self.exponent_q - 1) / confidence_radius))
        trust_radius_weight_gradient_gradient_nominator_1_2 = ((-2.0 * self.exponent_p * (distance / confidence_radius)**(2 * self.exponent_p) / confidence_radius) - ( 2.0 * self.exponent_q * (distance / confidence_radius)**(2 * self.exponent_q) / confidence_radius))

        trust_radius_weight_gradient_gradient_nominator_2 = ((-2.0 * self.exponent_p * (2 * self.exponent_p - 1) * (distance / confidence_radius)**(2 * self.exponent_p - 1) / confidence_radius**2) - (2.0 * self.exponent_p * (distance / confidence_radius)**(2 * self.exponent_p - 1) / confidence_radius**2) -
                                                                    (2.0 * self.exponent_q * (2 * self.exponent_q - 1) * (distance / confidence_radius)**(2 * self.exponent_q - 1) / confidence_radius**2) - (2.0 * self.exponent_q * (distance / confidence_radius)**(2 * self.exponent_q - 1) / confidence_radius**2))
        
        # trust_radius_weight_gradient_gradient_nominator_2 = ((-4.0 * self.exponent_p**2 * (distance / confidence_radius)**(2 * self.exponent_p) / (confidence_radius*distance)) -
        #                                                             (4.0 * self.exponent_q**2 * (distance / confidence_radius)**(2 * self.exponent_q) / (confidence_radius*distance)))
        

        trust_radius_weight_gradient_gradient = (((trust_radius_weight_gradient_gradient_nominator_1_1 * trust_radius_weight_gradient_gradient_nominator_1_2 ) / denominator**3) - trust_radius_weight_gradient_gradient_nominator_2 / denominator**2) * grad_s


        return trust_radius_weight_gradient_gradient
    
    def trust_radius_weight_gradient_gradient_test(self, confidence_radius, datapoint):
        
        distance, _, _, _, distance_vector, _ = self.cartesian_distance_symmetry(datapoint)
        denominator = (
                (distance / confidence_radius)**(2 * self.exponent_p) +
                (distance / confidence_radius)**(2 * self.exponent_q))
        
        
        trust_radius_weight_gradient_gradient_1_nominator = (((-2.0 * self.exponent_p**2 * distance_vector * distance**(2.0 * self.exponent_p - 2.0)) 
                                                /(confidence_radius**(2.0 * self.exponent_p) * confidence_radius))
                                                - (2.0 * self.exponent_q**2.0 * distance_vector * distance**(2.0 * self.exponent_q - 2.0)) 
                                                /(confidence_radius**(2.0 * self.exponent_q) * confidence_radius)) 

        trust_radius_weight_gradient_gradient_2_nominator = ( 4.0 * (((self.exponent_p * distance_vector * distance**(2.0 * self.exponent_p - 2.0)) 
                                                /(confidence_radius**(2.0 * self.exponent_p)))
                                                + ( self.exponent_q * distance_vector * distance**(2.0 * self.exponent_q - 2.0)) 
                                                /(confidence_radius**(2.0 * self.exponent_q))) * (((-2.0 * (distance/confidence_radius)**(2.0 * self.exponent_p) * self.exponent_p) / confidence_radius) 
                                                                                               - ((2.0 * (distance/confidence_radius)**(2.0 * self.exponent_q) * self.exponent_q) / confidence_radius))) 

        trust_radius_weight_gradient_gradient = -2.0 * trust_radius_weight_gradient_gradient_1_nominator / denominator**2 + trust_radius_weight_gradient_gradient_2_nominator / denominator**3



        return np.linalg.norm(trust_radius_weight_gradient_gradient)
    
    
    def shepard_weight_test(self, confidence_radius, dp):
        """ Returns the derivative of an unormalized Shepard interpolation
            weight with respect to the Cartesian coordinates in
            self.impes_coordinate

            :param distance_vector:
                The Cartesian distance vector between
                the current data_point and self.impes_coordinate.
            :param distance:
                The norm of the distance vector * sqrt(N), N number of atoms.
        """

        weight_gradient = None

        distance, _, _, _, distance_vector, _ = self.cartesian_distance_symmetry(dp)
        
    
        denominator = (
            (distance / confidence_radius)**(2 * self.exponent_p) +
            (distance / confidence_radius)**(2 * self.exponent_q))
        
        weight = 1.0 / denominator
        
        return weight
    
    def shepard_weight_gradient_test(self, dp, confidence_radius):
        """ Returns the derivative of an unormalized Shepard interpolation
            weight with respect to the Cartesian coordinates in
            self.impes_coordinate

            :param distance_vector:
                The Cartesian distance vector between
                the current data_point and self.impes_coordinate.
            :param distance:
                The norm of the distance vector * sqrt(N), N number of atoms.
        """

        weight_gradient = None

        distance, _, _, _, distance_vector, _ = self.cartesian_distance_symmetry(dp)
    
        denominator = (
            (distance / confidence_radius)**(2 * self.exponent_p) +
            (distance / confidence_radius)**(2 * self.exponent_q))
        derivative_p = (self.exponent_p * distance_vector *
                        distance**(2 * self.exponent_p - 2) /
                        confidence_radius**(2 * self.exponent_p))
        derivative_q = (self.exponent_q * distance_vector *
                        distance**(2 * self.exponent_q - 2) /
                        confidence_radius**(2 * self.exponent_q))
        
        weight_gradient = (-1.0 * (2.0 * (derivative_p + derivative_q)) *
                       (1.0 / (denominator**2)))
        # weight_gradient = (1.0 / (denominator**2))
        

        return weight_gradient
    
    def numerical_weight_gradient(self,
            X, Xi, A, exponent_p, exponent_q, confidence_radius,
            delta=1e-6, scheme="central"
        ):
            """
            Finite-difference gradient of the weight wrt X (shape (N,3)).
            A is held fixed (built from the dp's Hessian), so only X is perturbed.

            Parameters
            ----------
            delta : float
                Displacement magnitude in the same length units as X, Xi.
            scheme : "central" | "forward"
                Central difference is O(delta^2) and recommended.

            Returns
            -------
            grad : (N,3) array
                Numerical gradient d w / d X.
            """

            def shepard_weight_only(X, Xi, A, exponent_p, exponent_q, confidence_radius):
                """ Computes the Shepard weight for the given X and Xi. """
                distance_vector = X - Xi
                quad_distance = np.linalg.multi_dot([distance_vector.flatten().T, A, distance_vector.flatten()])
                denominator = (
                    ( quad_distance / confidence_radius)**(2 * exponent_p) +
                    (quad_distance / confidence_radius)**(2 * exponent_q))
               
                weights = 1.0 / denominator
                return weights
            
            X = np.asarray(X).copy()
            Xi = np.asarray(Xi)
            N = X.shape[0]
            grad = np.zeros_like(X)

            # base evaluation if using forward scheme
            if scheme == "forward":
                w0 = shepard_weight_only(X, Xi, A, exponent_p, exponent_q, confidence_radius)

            # loop over coordinates
            for i in range(N):
                for a in range(3):
                    if scheme == "central":
                        X[i, a] += delta
                        wp = shepard_weight_only(X, Xi, A, exponent_p, exponent_q, confidence_radius)
                        X[i, a] -= 2*delta
                        wm = shepard_weight_only(X, Xi, A, exponent_p, exponent_q, confidence_radius)
                        X[i, a] += delta  # restore
                        grad[i, a] = (wp - wm) / (2.0 * delta)
                    elif scheme == "forward":
                        X[i, a] += delta
                        wp = shepard_weight_only(X, Xi, A, exponent_p, exponent_q, confidence_radius)
                        X[i, a] -= delta   # restore
                        grad[i, a] = (wp - w0) / delta
                    else:
                        raise ValueError("scheme must be 'central' or 'forward'")

            return grad
    
    def shepard_weight_gradient_hessian(self, distance_vector, quad_distance, confidence_radius, grad_s):
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
            ( quad_distance / confidence_radius)**(2 * self.exponent_p) +
            (quad_distance / confidence_radius)**(2 * self.exponent_q))
        
        
        derivative_p = ( self.exponent_p / confidence_radius *
                        (quad_distance/confidence_radius)**(2 * self.exponent_p - 1))
        derivative_q = (self.exponent_q / confidence_radius *
                        (quad_distance / confidence_radius)**(2 * self.exponent_q - 1))
        
        weight_gradient = (-1.0 * (2.0 * (derivative_p + derivative_q) * (1.0 / (denominator**2)))) * grad_s
        # weight_gradient = weight_gradient.reshape(-1, 3)
        
                    
        return  denominator, weight_gradient

    def shepard_weight_gradient(self, distance_vector, distance, confidence_radius,):
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
        
        weight_gradient = (-1.0 * (2.0 * (derivative_p + derivative_q)) *
                       (1.0 / (denominator**2)))
        
        
   
        return  denominator, weight_gradient
    
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
       
        confidence_radius = data_point.confidence_radius

        # natm = data_point.cartesian_coordinates.shape[0]
        # deriv_int_denominator = 0.0
        # deriv_int_nominator = np.zeros_like(distance_vector)
        dihedral_dist = 0.0
        # for i, element in enumerate(self.impes_coordinate.z_matrix[self.symmetry_information[-1][1]:], start=self.symmetry_information[-1][1]): 

        #     # print('ELEMENT', sorted(element))
        #     if tuple(sorted(element)) in self.symmetry_information[6][3]:
        #         dihedral_dist =+ 1/2*(1 + np.sin(3.0*(self.impes_coordinate.internal_coordinates_values[i] - data_point.internal_coordinates_values[i]) - np.pi/2))
        #         deriv_int_nominator += (self.impes_coordinate.b_matrix[i,:] * 
        #                     3.0/2.0 * np.cos(3.0*(self.impes_coordinate.internal_coordinates_values[i] - data_point.internal_coordinates_values[i]) + np.pi/2.0)).reshape(natm, 3)
        #         deriv_int_denominator += dihedral_dist
        
        # if abs(deriv_int_denominator) < 1e-6:
        #     deriv_int_denominator = 1e-6
        #     deriv_int_nominator[:] = 0 
        # print(deriv_int_denominator.shape, distance_vector.shape)
        if distance < 1e-9:
            distance = 1e-8
            distance_vector[:] = 0
        if self.interpolation_type == 'shepard':
            denominator, weight_gradient = self.shepard_weight_gradient(
                distance_vector, distance, confidence_radius, calc_cart=True)
        elif self.interpolation_type == 'simple':
            weight_gradient = self.simple_weight_gradient(
                distance_vector, distance)
        else:
            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)
        
        # print('dihedral_dist', 1/dihedral_dist, 1/denominator)
        # print('deriv_int_denominator', deriv_int_denominator)


        return distance, denominator, weight_gradient, distance_vector, dihedral_dist
    def calculate_translation_coordinates(self, given_coordinates):
        """Center the molecule by translating its geometric center to (0, 0, 0)."""
        center = np.mean(given_coordinates, axis=0)
        translated_coordinates = given_coordinates - center

        return translated_coordinates

    def cartesian_distance_symmetry(self, data_point):
        """Calculates and returns the cartesian distance between
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.

           :param data_point:
                InterpolationDatapoint object
        """

        def apply_Jc(A):   # A: (Nc,3)
            return A - A.mean(axis=0, keepdims=True)

        def make_spd_hessian(
            Hc, masses=None, floor_abs=None, floor_rel=1e-6, cap_rel=None,
            ridge_delta=None, use_mass_weight=False
        ):
            """
            Returns an SPD version of the Cartesian Hessian with the same units as Hc.
            Hc: (3N,3N) in Eh/bohr^2 (or any consistent units).
            masses: (N,) in atomic mass units (optional, only if use_mass_weight=True).
            floor_abs: absolute floor (same units as Hc), e.g. 1e-6 Eh/bohr^2
            floor_rel: relative floor as tau * median(diag(K or Hs)) if floor_abs is None
            cap_rel:  optional cap as kappa * percentile95(eigs)
            ridge_delta: alternative to eigen-flooring: target min eigenvalue (>0)
            use_mass_weight: if True, regularize K = M^{-1/2} H M^{-1/2} and map back
            """
            # 1) Symmetrize
            Hs = 0.5 * (Hc + Hc.T)

            if use_mass_weight:
                assert masses is not None
                N = len(masses)
                M = np.repeat(masses, 3)
                Minv_sqrt = 1.0 / np.sqrt(M)
                # K = M^{-1/2} Hs M^{-1/2}
                K = (Minv_sqrt[:,None] * Hs) * Minv_sqrt[None,:]
                A = K
            else:
                A = Hs

            # 2) SPD via ridge or eigen-floor
            if ridge_delta is not None:
                # Ridge shift
                w = np.linalg.eigvalsh(A)
                lam_min = w[0]
                mu = max(ridge_delta - lam_min, 0.0)
                Aspd = A + mu * np.eye(A.shape[0])
            else:
                # Eigen-floor (and optional cap)
                w, U = np.linalg.eigh(A)
                if floor_abs is None:
                    # relative floor from median diagonal as a robust scale
                    base = np.median(np.diag(A))
                    floor_abs = abs(floor_rel * base)
                w = np.maximum(np.abs(w), floor_abs)
                if cap_rel is not None:
                    cap = cap_rel * np.percentile(w, 95)
                    w = np.minimum(w, cap)
                Aspd = (U * w) @ U.T   # U @ diag(w) @ U.T

            if use_mass_weight:
                # Map back: Hspd = M^{1/2} Aspd M^{1/2}
                Msqrt = np.sqrt(np.repeat(masses, 3))
                Hspd = (Msqrt[:,None] * Aspd) * Msqrt[None,:]
            else:
                Hspd = Aspd

            # Final symmetrization to remove numeric dust
            return 0.5 * (Hspd + Hspd.T)
        

        
        # First, translate the cartesian coordinates to zero
        target_coordinates = data_point.cartesian_coordinates.copy()
        reference_coordinates = self.impes_coordinate.cartesian_coordinates.copy()

        target_coordinates_core = np.delete(target_coordinates, self.symmetry_information[4], axis=0)
        reference_coordinates_core = np.delete(reference_coordinates, self.symmetry_information[4], axis=0)

        center_target_coordinates_core = self.calculate_translation_coordinates(target_coordinates_core)
        center_reference_coordinates_core = self.calculate_translation_coordinates(reference_coordinates_core)

        # Then, determine the rotation matrix which
        # aligns data_point (target_coordinates)
        # to self.impes_coordinate (reference_coordinates)

        # rotation_matrix = geometric.rotate.get_rot(center_target_coordinates_core,
        #                                            center_reference_coordinates_core)

        # # Rotate the data point
        # rotated_coordinates = np.dot(rotation_matrix, center_target_coordinates_core.T).T

        
        # # rotation_weights = geometric.rotate.get_rot(rotated_coordinates,
        # #                                             target_coordinates)

        # org_reference_coord = center_reference_coordinates_core.copy()
        # org_mapping_list = [i for i in range(reference_coordinates.shape[0])]
        # mapping_list = [[i for i in range(reference_coordinates.shape[0])]]
        # swapped = None
    
        # non_group_atoms = []
        # for i in range(len(self.molecule.get_labels())):
        #     if i not in self.symmetry_information[1][0]:
        #         non_group_atoms.append(i)
        # current_molecule_trunc_coords = np.delete(reference_coordinates, self.symmetry_information[3], axis=0)
        # data_point_molecule_trun_coords = np.delete(rotated_coordinates, self.symmetry_information[3], axis=0)
        # mapping_list, mapping_dict, assigned = self.perform_symmetry_assignment(self.symmetry_information[0], self.symmetry_information[3], 
        #                                                           current_molecule_trunc_coords, data_point_molecule_trun_coords)
        # org_reference_coord[org_mapping_list] = reference_coordinates[mapping_list[0]]

        # org_reference_coord[self.symmetry_information[3]] = rotated_coordinates[self.symmetry_information[3]]
        assigned = False
        
        # if assigned and 1 == 2:
            
        #     swapped = (org_mapping_list, mapping_list[0], reference_coordinates, mapping_dict)
        #     print(swapped)
            # exit()
        # mapping_list = atom_mapper.perform(atom_map_1=self.symmetry_information[0], env_groups_1=self.symmetry_information[1])

        # # # # # Calculate the Cartesian distance
        # rotation_matrix = geometric.rotate.get_rot(target_coordinates, org_reference_coord)
        # # Rotate the data point
        # rotated_coordinates = np.dot(rotation_matrix, target_coordinates.T).T

        # # Calculate the Cartesian distance
        # distance = (np.linalg.norm(rotated_coordinates - org_reference_coord))
        # # Calculate the gradient of the interpolation weights
        # # (required for energy gradient interpolation)
        # distance_vector = (org_reference_coord - rotated_coordinates)
        distance = 0
        distance_vector = 0
        grad_s = 1.0
        Xc = center_reference_coordinates_core     # variable (current geometry)
        Yc = center_target_coordinates_core        # fixed (datapoint)
        R  = geometric.rotate.get_rot(Yc, Xc)
        dU = geometric.rotate.get_rot_der(Yc, Xc)   

        rotated_current = (R @ Yc.T).T
        distance_vector =  Xc - rotated_current
        distance = np.linalg.norm(distance_vector)

    
        # Xc = center_reference_coordinates_core     # variable (current geometry)
        # Yc = center_target_coordinates_core        # fixed (datapoint)
        # R_x  = geometric.rotate.get_rot(Xc, Yc)
        # dU_x = geometric.rotate.get_rot_der(Xc, Yc)   
        
        # rotated_current = (R_x @ Xc.T).T
        # distance_vector =  rotated_current - Yc
        # distance = np.linalg.norm(distance_vector) 

        # Nc = Xc.shape[0]
        # J = np.eye(Nc) - np.ones((Nc, Nc))/Nc      # (Nc,Nc)
        # J3 = np.kron(J, np.eye(3))                 # (3Nc,3Nc)

        # Hc = J3.T @ data_point.hessian @ J3 
        # E_ref = 1e-3 / 0.001593 # in hartree to scale the distance making it unitless
        
        # Hc_ul = Hc * E_ref
        # H_spd = make_spd_hessian(Hc_ul, floor_abs=1e-6)
        # # distance (use a small epsilon!)
        # eps = 1e-8
        #      # <-- this is the “quad_distance”
 
        # distance_2     = (distance_vector.reshape(-1) @ ((1.0 * H_spd) @ distance_vector.reshape(-1)))
        # org_distance = (distance_vector.reshape(-1) @ ((1.0 * H_spd) @ distance_vector.reshape(-1)))
        # distance = np.sqrt(distance_2)
        # # print('distance check', distance, 1/distance)
        # # if 1/distance > 2000:
        # #     exit()
        # y = ((1.0 * H_spd) @ distance_vector.reshape(-1)).reshape(-1, 3)      # (Nc,3)

        # rigid = y @ R_x                           # (Nc,3)
        # S = y.T @ Xc                                    # (3,3)

        # resp = np.einsum('baij,ij->ba', dU_x, S, optimize=True)   # (Nc,3)

        # z_mat = rigid + resp                            # (Nc,3)
        # z_mat = apply_Jc(z_mat)                        # (Nc,3)

        # grad_s = 1.0 * z_mat / distance                            # (Nc,3) 

        # def s_of_X(X_core, Yc):
        #     Xc = apply_Jc(X_core)
        #     R  = geometric.rotate.get_rot(Xc, Yc)
        #     Delta = (R @ Xc.T).T - Yc
        #     return float(Delta.reshape(-1) @ (H @ Delta.reshape(-1)))

        # u = np.random.randn(*Xc.shape); u = apply_Jc(u)   # kill translations
        # u /= np.linalg.norm(u)
        # for eps in [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]:
        #     s_plus  = s_of_X(Xc + eps*u, Yc)
        #     s_minus = s_of_X(Xc - eps*u, Yc)
        #     lhs = (s_plus - s_minus) / (2*eps)
        #     rhs = (grad_s.reshape(Xc.shape) * u).sum()
        #     print(f"eps={eps:.0e}, FD={lhs}, grad={rhs}, diff={lhs-rhs}")
        # # exit()
     
        # rotation_matrix_org = geometric.rotate.get_rot(target_coordinates,
        #                                            org_reference_coord)

        # # Rotate the data point
        # rotated_coordinates_org = np.dot(rotation_matrix_org, target_coordinates.T).T
        
        

        # # rotation_weights = geometric.rotate.get_rot(rotated_coordinates,
        # #                                             target_coordinates)

        # # Calculate the Cartesian distance
        # distance_org = (np.linalg.norm(rotated_coordinates_org - org_reference_coord))
        # distance_vec_org = ((rotated_coordinates_org - org_reference_coord))



        # for i in range(len(distance_vector)):

        #     print('Distance vector', i, distance, distance_org, distance_vector[i], distance_vec_org[i])
        
        # # print('\n\n')
        # natm = data_point.cartesian_coordinates.shape[0]
        
        # # for i in range(15, len(self.z_matrix)):
        # deriv_int_denominator = 0.0
        # deriv_int_nominator = np.zeros_like(distance_vector)

        # for i in range(15, len(self.z_matrix)):
        #     dihedral_dist = 1/2*(1 + np.sin(3.0*(self.impes_coordinate.internal_coordinates_values[i] - data_point.internal_coordinates_values[i]) - np.pi/2))
        #     deriv_int_nominator += (self.impes_coordinate.b_matrix[i,:] * 
        #                 3.0/2.0 * np.cos(3.0*(self.impes_coordinate.internal_coordinates_values[i] - data_point.internal_coordinates_values[i]) + np.pi/2.0)).reshape(natm, 3)
        #     deriv_int_denominator += (dihedral_dist)**2

        # if abs(deriv_int_denominator) < 1e-6:
        #     deriv_int_denominator = 1e-6
        #     deriv_int_nominator[:] = 0 
        
        # distance_vector_norm = np.zeros(reference_coordinates.shape[0])
        # for i in range(len(distance_vector_norm)):
        #     distance_vector_norm[i] += np.linalg.norm(distance_vector[i])
        confidence_radius = data_point.confidence_radius
        # w, v = np.linalg.eigh(data_point.hessian)
        # w_reg = np.where(w > 1e-8, w, 1e-8)   # clip eigenvalues
        # H_inv = (v / w_reg**0.5) @ v.T
        # w, V = np.linalg.eigh(data_point.hessian)
        # # set floor relative to a robust positive scale if available
        # pos = w[w > 0]
        # lam_floor = (np.median(pos) if pos.size else 1.0) * 1e-6
        # w_clipped = np.clip(w, lam_floor, None)
        # # Hq^{-gamma} (gamma=1 -> inverse; gamma=0.5 -> inverse sqrt)
        # invpow = w_clipped**(0.5)
        # Mq = (V * invpow) @ V.T
        # H_inv = 0.5*(Mq + Mq.T)
        # H_iden = np.ones_like(data_point.hessian)
        # I = np.eye(data_point.cartesian_coordinates.shape[0])
        # Rbig = np.kron(I, rotation_matrix)            # (3N,3N)

        # rotated_hessian = np.linalg.multi_dot([Rbig, data_point.hessian, Rbig.T])
        # quad_distance = np.linalg.multi_dot([distance_vector.flatten().T, rotated_hessian, distance_vector.flatten()])
        # print('quad_distance', quad_distance, distance)

        if not isinstance(confidence_radius, (float, int)):
            confidence_radius = 0.5

        dihedral_dist = 0.0
        if distance < 1e-8:
            distance = 1e-8
            distance_vector[:] = 0
            s = 1e-8

        # if distance_x < 1e-8:
        #     distance_x = 1e-8
        #     distance_vector_x[:] = 0
        #     s = 1e-8

        if self.interpolation_type == 'shepard':
            # weight_grad_numerical = self.numerical_weight_gradient(org_reference_coord, rotated_coordinates, rotated_hessian, self.exponent_p, self.exponent_q, confidence_radius)
            # denominator, weight_gradient = self.shepard_weight_gradient_hessian(
            #     distance_vector, distance, confidence_radius, grad_s)
            denominator, weight_gradient = self.shepard_weight_gradient(
                distance_vector, distance, confidence_radius)
            
            # denominator_x, weight_gradient_x = self.shepard_weight_gradient(
            #     distance_vector_x, distance_x, confidence_radius)
            
            # weight_gradient_x = weight_gradient_x - weight_gradient_x.mean(axis=0, keepdims=True)

            # print('Difference in weight gradient', weight_gradient - weight_gradient_x)
            # exit()
            # print('Numerical weight gradient', weight_grad_numerical)
            # print('Analytical weight gradient', weight_gradient_hess)
            # print('difference', weight_grad_numerical - weight_gradient_hess)
        elif self.interpolation_type == 'simple':
            weight_gradient = self.simple_weight_gradient(
                distance_vector, distance)
        else:
            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)
        
        # print('weight gradient hess before the rotational derivative', weight_gradient_hess)
        # weight_gradient_hess = weight_gradient_hess * grad_s
        # print('Different distance metrics:', s, distance)
        # print('weight hess', weight_gradient_hess, 'orginal hessian', weight_gradient)
       
        # distance = s
        # exit()

        return distance, dihedral_dist, denominator, weight_gradient, distance_vector, grad_s


    def determine_important_internal_coordinates(self, qm_energy, qm_gradient, molecule, z_matrix, datapoints, error_source='energy'):
        """Determines the most important internal coordinates
           leading to the large deviation.
        """
        constraints = []
        for datapoint in datapoints:
            internal_coord_elem_distance = []
            org_interal_coord_elem_distance = []
            for elem_idx, element in enumerate(self.z_matrix): 
                if len(element) == 4 and tuple(sorted(element)) in self.symmetry_information[7][3]: 
                    # internal_coord_elem_distance.append(0.5 * ( 1.0 + np.cos(3.0 * (self.impes_coordinate.internal_coordinates_values[elem_idx] - datapoint.internal_coordinates_values[elem_idx]) + np.pi)))
                    internal_coord_elem_distance.append(0.0)
                    org_interal_coord_elem_distance.append(0.0)
                
                elif len(element) == 4 and tuple(sorted(element)) in self.symmetry_information[7][2]: 
                    # internal_coord_elem_distance.append(0.5 * ( 1.0 + np.cos(2.0 * (self.impes_coordinate.internal_coordinates_values[elem_idx] - datapoint.internal_coordinates_values[elem_idx]) + np.pi)))                
                    internal_coord_elem_distance.append(0.0)
                    org_interal_coord_elem_distance.append(0.0)
                elif len(element) == 4:
                    
                    internal_coord_elem_distance.append(np.sin(self.impes_coordinate.internal_coordinates_values[elem_idx] - datapoint.internal_coordinates_values[elem_idx]))
                    org_interal_coord_elem_distance.append(np.cos(self.impes_coordinate.internal_coordinates_values[elem_idx] - datapoint.internal_coordinates_values[elem_idx]))
                else:
                    org_interal_coord_elem_distance.append(1.0)
                    internal_coord_elem_distance.append(self.impes_coordinate.internal_coordinates_values[elem_idx] - datapoint.internal_coordinates_values[elem_idx])

            N = len(z_matrix)

            partial_energies = np.zeros(N)
            partial_gradient = np.zeros(N) 

        
            # First handle linear parts
            for i in range(N):
                partial_energies[i] += internal_coord_elem_distance[i] * datapoint.internal_gradient[i]
                partial_gradient[i] += org_interal_coord_elem_distance[i] * datapoint.internal_gradient[i]


            # Now handle quadratic parts (diagonal + cross)
            for i in range(N):
                # Diagonal
                partial_energies[i] += 0.5 * internal_coord_elem_distance[i] * datapoint.internal_hessian[i,i] * internal_coord_elem_distance[i]
                partial_gradient[i] += org_interal_coord_elem_distance[i] * datapoint.internal_hessian[i,i] * internal_coord_elem_distance[i]

                # Cross-terms: j>i to avoid double counting in sum
                for j in range(i+1, N):
                    cross = internal_coord_elem_distance[i] * datapoint.internal_hessian[i,j] * internal_coord_elem_distance[j]
                    cross_prime_ij = org_interal_coord_elem_distance[i] * datapoint.internal_hessian[i,j] * internal_coord_elem_distance[j]
                    cross_prime_ji = org_interal_coord_elem_distance[j] * datapoint.internal_hessian[j,i] * internal_coord_elem_distance[i]

                    # Split cross equally between i and j
                    partial_energies[i] += 0.5 * cross
                    partial_energies[j] += 0.5 * cross
                    partial_gradient[i] += cross_prime_ij
                    partial_gradient[j] += cross_prime_ji

            # The sum of partial_energies should match the total second-order approx:
            # total_energy_diff = sum(partial_energies)

            variable_part = sum(partial_energies)        # Hartree
            E_pred_check  = datapoint.energy + variable_part

            pred_E, pred_G, _ = self.compute_potential(datapoint, self.impes_coordinate.internal_coordinates_values)
            pred_im_G_int = self.transform_gradient_to_internal_coordinates(molecule, pred_G, self.impes_coordinate.b_matrix)
            pred_qm_G_int = self.transform_gradient_to_internal_coordinates(molecule, qm_gradient, self.impes_coordinate.b_matrix)
            print('Difference in the gradient', pred_im_G_int - pred_qm_G_int)

            print("max energy diff",
                abs((E_pred_check - pred_E)))
            
            g_pred_check = partial_gradient              # already Hartree/bohr
            max_grad_diff = np.max(np.abs(g_pred_check - pred_im_G_int))
            print("max grad diff", max_grad_diff)
    
            delta_E = abs(qm_energy - pred_E)
            delta_G = np.linalg.norm(pred_qm_G_int - pred_im_G_int)

            delta_g = pred_qm_G_int - pred_im_G_int

            single_energy_error = []
            weights = []
            single_gradient_error = []
            abs_pg       = np.abs(partial_gradient)
            sum_abs_pg   = abs_pg.sum()

            eps = 1e-12

            abs_pE       = np.abs(partial_energies)
            
            wE_raw  = abs_pE * np.abs(delta_g)            # or: np.abs(Delta_q_eff * delta_g)
            w_E_1    = wE_raw / (wE_raw.sum() + eps)
            single_E_err_1 = delta_E * w_E_1
            
            w_E_2          = abs_pE / (abs_pE.sum() + eps)     # weights  Σw=1
            single_E_err = delta_E * w_E_2                     # each coord’s share

            print(' \n\n Energy error with QM 1: ', single_E_err_1 * hartree_in_kcalpermol(), 
                  "sum_1_kcal =", single_E_err_1.sum() * hartree_in_kcalpermol(), '2', single_E_err * hartree_in_kcalpermol(),
                  "sum_1_kcal =", single_E_err.sum() * hartree_in_kcalpermol(), '\n\n')

            # ---- gradient branch (the bit you asked about) -----------------
            # Here:  weight_i  = |partial_G_i|  / Σ|partial_G|
            abs_pG       = np.abs(partial_gradient)
            sum_abs_pG   = abs_pG.sum()
            if sum_abs_pG < eps:                             # force already perfect
                w_G          = np.zeros_like(abs_pG)
                single_G_err = np.zeros_like(abs_pG)
            else:
                w_G          = abs_pG / sum_abs_pG
                single_G_err = delta_G * w_G

            assert np.allclose(single_E_err.sum(), delta_E, atol=1e-12)
            assert np.allclose(single_G_err.sum(), delta_G, atol=1e-12)

            print('Error with QM ', np.allclose(single_E_err.sum(), delta_E, atol=1e-12), np.allclose(single_G_err.sum(), delta_G, atol=1e-12))
            
            for grad_idx, individual_contrib in enumerate(partial_energies):
                
                energy_error = delta_E * (abs(individual_contrib) / sum(abs(partial_energies)))  # each coord’s share
                single_weight = abs(individual_contrib) / sum(abs(partial_energies))
                weights.append(single_weight)
                single_energy_error.append(energy_error * hartree_in_kcalpermol())

                gradient_error = delta_G * w_G[grad_idx]     # share of |Δg|
                single_gradient_error.append(gradient_error  * hartree_in_kcalpermol() * bohr_in_angstrom())
            
            # Pair each energy error with its corresponding partial energy and internal coordinate
            contributions_energy = list(zip(partial_energies, single_E_err_1 * hartree_in_kcalpermol(), weights, z_matrix))
            # contributions_gradient = list(zip(partial_gradient, single_gradient_error, w_G, z_matrix))


            # Sort contributions by the energy error in descending order
            sorted_contributions = sorted(contributions_energy, key=lambda x: x[1], reverse=True)
            print('Delta E:', delta_E * hartree_in_kcalpermol(), constraints, error_source, weights)

            # if error_source == 'gradient':
            #     print('Delta G:', delta_G * hartree_in_kcalpermol() * bohr_in_angstrom(), constraints)
            #     sorted_contributions = sorted(contributions_gradient, key=lambda x: x[1], reverse=True)

                # # Print the sorted contributions with internal coordinates
                # print('Delta E:', delta_G * hartree_in_kcalpermol() * bohr_in_angstrom(), constraints)
                # for contrib, error, ind_weight, coord in sorted_contributions_gradient[:]:        
                #     print('additional coord', abs(ind_weight - max(weights)) < 1e-8)
                #     if tuple(int(x) for x in coord) in constraints:
                #         continue
                #     if len(coord) == 2 and ind_weight > max(weights) * 0.7:
                #         constraints.append(tuple(int(x) for x in coord))
                #     elif len(coord) == 3 and ind_weight > max(weights) * 0.7:
                #         constraints.append(tuple(int(x) for x in coord))
                #     elif len(coord) == 4 and ind_weight > max(weights) * 0.7 and tuple(sorted(coord)) not in self.symmetry_information[7][3] and tuple(sorted(coord)) not in self.symmetry_information[7][2]:
                #         constraints.append(tuple(int(x) for x in coord))
                #     print(f'Internal Coordinate: {tuple(int(x) for x in coord)}, Error: {error * hartree_in_kcalpermol() * bohr_in_angstrom()} kcal/mol Å^-1')#, distance {internal_coord_elem_distance[z_matrix.index(coord)]}, Contribution: {contrib}, weight {ind_weight}, Error: {error * hartree_in_kcalpermol()} kcal/mol')
                # print('Sum of Weights', sum(weights), sum(single_gradient_error) * hartree_in_kcalpermol() * bohr_in_angstrom())


            # Print the sorted contributions with internal coordinates
            
            for contrib, error, ind_weight, coord in sorted_contributions[:]:        
                print('additional coord', abs(ind_weight - max(weights)) < 1e-8)
                if tuple(int(x) for x in coord) in constraints:
                    continue
                if len(coord) == 2 and ind_weight > max(weights) * 0.7:
                    constraints.append(tuple(int(x) for x in coord))
                elif len(coord) == 3 and ind_weight > max(weights) * 0.7:
                    constraints.append(tuple(int(x) for x in coord))
                elif len(coord) == 4 and ind_weight > max(weights) * 0.7 and tuple(sorted(coord)) not in self.symmetry_information[7][3] and tuple(sorted(coord)) not in self.symmetry_information[7][2]:
                    constraints.append(tuple(int(x) for x in coord))
                print(f'Internal Coordinate: {tuple(int(x) for x in coord)}, Error: {error} kcal/mol')#, distance {internal_coord_elem_distance[z_matrix.index(coord)]}, Contribution: {contrib}, weight {ind_weight}, Error: {error * hartree_in_kcalpermol()} kcal/mol')
            print('Sum of Weights', sum(weights), sum(single_energy_error))

        return constraints
    
    def transform_gradient_to_internal_coordinates(self, molecule, gradient, b_matrix, tol=1e-6):
        """
        Transforms the gradient from Cartesian to internal coordinates.
        :param tol:
            Tolerance for the singular values of the B matrix.
        """
        
        dimension = gradient.shape[0] * 3 - 6
        if gradient.shape[0] == 2:
            dimension += 1
        
        g_matrix = np.dot(b_matrix, b_matrix.T)
        
        U, s, Vt = np.linalg.svd(g_matrix)
        
        print('eigenvalues', s)
        # Make zero the values of s_inv that are smaller than tol
        s_inv = np.array([1 / s_i if s_i > tol else 0.0 for s_i in s])
        
        # Critical assertion, check that the remaining positive values are equal to dimension (3N-6)
        number_of_positive_values = np.count_nonzero(s_inv)
        
        # If more elements are zero than allowed, restore the largest ones
        if number_of_positive_values > dimension:
            print('InterpolationDatapoint: The number of positive singular values is not equal to the dimension of the Hessian., restoring the last biggest elements')
            # Get indices of elements originally set to zero
            zero_indices = np.where(s_inv == 0.0)[0]
            
            # Sort these indices based on their corresponding values in 's' (descending order)
            sorted_zero_indices = zero_indices[np.argsort(s[zero_indices])[::-1]]
            
            # Calculate how many elements need to be restored
            num_to_restore = abs(number_of_positive_values - dimension)
            
            # Restore the largest elements among those set to zero
            for idx in sorted_zero_indices[:num_to_restore]:
                s_inv[idx] = 1 / s[idx]
        g_minus_matrix = np.dot(U, np.dot(np.diag(s_inv), Vt))
        gradient_flat = gradient.flatten()
        internal_gradient = np.dot(g_minus_matrix, np.dot(b_matrix, gradient_flat))

        return internal_gradient

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
        print('row col', row, col)
        if not np.equal(row, col).all():
            assigned = True
            
            # atom_maps = self.linear_assignment_solver(cost)
            reordred_arr = np.array(sym_group)[col]
            new_map[sym_group] = new_map[reordred_arr]

            for org, new in zip(np.array(sym_group), reordred_arr):
                print('org, neww', org, new)
            mapping_dict = {org: new for org, new in zip(np.array(sym_group), reordred_arr)}
        
        return [new_map], mapping_dict, assigned