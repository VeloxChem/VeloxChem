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
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import math
import random
import torch
import gpytorch
from typing import List, Tuple, Dict, Any, Optional
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
from .gprinterpolationdriver import GPRInterpolationDriver

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
        self.gpr_intdriver = None
        
        self.dihedral_function = None
        # simple or Shepard interpolation
        self.interpolation_type = 'shepard'
        self.weightfunction_type = 'cartesian'
        self.exponent_p = None
        self.scaling_time = False
        self.exponent_q = None
        self.confidence_radius = 0.5
        self.alpha = 1.0
        
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

        self.external_weights = None

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

        # Confidence radii optimization information
        self.dw_dalpha_list = []
        self.dw_dX_dalpha_list = []
        self.calc_optim_trust_radius = False

        self.bond_rmsd = None
        self.angle_rmsd = None
        self.dihedral_rmsd = None
        self.averaged_int_dist = None
        # Name lables for the QM data points
        self.labels = None
        self.use_inverse_bond_length = True
        self.use_eq_bond_length = False
        self.use_cosine_dihedral = False
        self.use_tc_weights = True


        self._input_keywords = {
            'im_settings': {
                'interpolation_type':
                    ('str', 'type of interpolation (simple/Shepard)'),
                'weightfunction_type':
                    ('str', 'type of interpolation (cartesian/cartesian-hessian)'),
                'exponent_p': ('int', 'the main exponent'),
                'exponent_q': ('int', 'the additional exponent (Shepard IM)'),
                'confidence_radius': ('float', 'the confidence radius'),
                'imforcefield_file':
                    ('str', 'the name of the chk file with QM data'),
                    'use_inverse_bond_length': ('bool', 'whether to use inverse bond lengths in the Z-matrix'),
                    'use_eq_bond_length': ('bool', 'whether to use eq bond lengths in the Z-matrix'),
                    'use_cosine_dihedral':('bool', 'wether to use cosine and sin for the diehdral in the Z-matrix'),
                    'use_tc_weights':('bool', 'weither to use target coustomized weights'),
                'labels': ('seq_fixed_str', 'the list of QM data point labels'),
            }
        }


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
        
        elif self.impes_coordinate.use_eq_bond_length:
            remove_from_label += "_eq"
            z_matrix_label += '_eq'
        
        else:
            remove_from_label += "_r"
            z_matrix_label += '_r'
        
        if self.impes_coordinate.use_cosine_dihedral:
            remove_from_label += "_cosine"
            z_matrix_label += '_cosine'
        else:
            remove_from_label += "_dihedral"
            z_matrix_label += '_dihedral'
        
        remove_from_label += "_energy"

        
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
                        current_z_list = [tuple(z_list.tolist()) for z_list in list(h5f.get(z_label))]
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

    def compute(self, molecule):
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
   
        self.define_impes_coordinate(molecule.get_coordinates_in_bohr())

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
        for i, data_point in enumerate(self.qm_data_points[:2]):
            
            distance, weight_gradient, _, swapped = self.cartesian_distance(data_point)

            if abs(distance) < min_distance:
                min_distance = abs(distance)

            distances_and_gradients.append((distance, i, weight_gradient, swapped))

        # min_distance, distances_and_gradients = self.process_in_parallel()

        used_labels = []
        close_distances = [
            (self.qm_data_points[index], distance, wg, index, swapped_tuple) 
            for distance, index, wg, swapped_tuple in distances_and_gradients 
            if abs(distance) <= min_distance + self.distance_thrsh]
        for qm_data_point, distance, weight_grad, label_idx, swapped_tuple_obj in close_distances:
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

            potential = self.compute_potential(qm_data_point, reordered_int_coord_values)
            gradient = self.compute_gradient(qm_data_point, swapped_tuple_obj, reordered_int_coord_values, org_b_matrix_cp)
            gradients.append(gradient)
            potentials.append(potential)
            weights.append(weight)
            self.potentials.append(potential)
            weight_gradients.append(weight_grad)
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
        weights = np.array(weights, dtype=np.float64) / sum_weights
        potentials = np.array(potentials, dtype=np.float64)
        gradients = np.array(gradients, dtype=np.float64)
        weight_gradients = np.array(weight_gradients, dtype=np.float64)

        for i, weight_i in enumerate(weights):
            self.int_coord_weight[self.labels[used_labels[i]]] = weight_i

        self.impes_coordinate.energy = np.dot(weights, potentials)
        self.impes_coordinate.gradient = (
            np.tensordot(weights, gradients, axes=1) +
            np.tensordot(potentials, weight_gradients, axes=1) / sum_weights -
            self.impes_coordinate.energy * sum_weight_gradients / sum_weights)



    def shepard_interpolation(self):
        """Performs a simple interpolation.

        :param qm_data_points:
            A list of InterpolationDatapoint corresponding to configurations
            to be used for interpolation (which have been calculated
            quantum mechanically).
        """

        def numerical_gradient(f, X, dp, eps=1e-6):
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
            grad_da_dw = np.zeros_like(X, dtype=float)

            Nc = X.shape[0]
            for a in range(Nc):
                for mu in range(3):
                    dX = np.zeros_like(X)
                    dX[a, mu] = eps
                    current_geometry_plus = X + dX
                    self.define_impes_coordinate(current_geometry_plus)
                    _, _, denominator_p, _, _, _, dw_da_p, _ = f(dp)
                    f_plus = 1.0 / denominator_p

                    # f_plus, _, _  = f(dp, self.impes_coordinate.internal_coordinates_values)
                    current_geometry_minus = X - dX
                    self.define_impes_coordinate(current_geometry_minus)
                    # f_minus, _, _  = f(dp, self.impes_coordinate.internal_coordinates_values)
                    _, _, denominator_m, _, _, _, dw_da_m, _ = f(dp)
                    f_minus = 1.0 / denominator_m
                    grad[a, mu] = (f_plus - f_minus) / (2 * eps)
                    if self.calc_optim_trust_radius:
                        grad_da_dw[a, mu] = (dw_da_p - dw_da_m) / (2 * eps)
            
            return grad, grad_da_dw
        

        natms = self.impes_coordinate.cartesian_coordinates.shape[0]

        sum_weights_cart = 0.0
        potentials = []
        gradients = []
        hessian_error = []
        self.weights = {}
        weights_cart = []
        averaged_int_dists = []
        weight_gradients_cart = []
        used_labels = []
        self.potentials = []
        self.gradients = []
        self.dw_dalpha_list = []
        self.dw_dX_dalpha_list = []

        masses = self.molecule.get_masses().copy()
        masses_cart = np.repeat(masses, 3)
        sqrt_masses = np.sqrt(masses_cart)

        sum_weight_gradients_cart = np.zeros((natms, 3))

        distances_and_gradients = []
        min_distance = float('inf')
        self.time_step_reducer = False
        
        self.calc_optim_trust_radius = True
        
        if self.weightfunction_type == 'cartesian':
            for i, data_point in enumerate(self.qm_data_points[:]):
                
                if self.external_weights is not None:
                    data_point.confidence_radius = self.external_weights[i]
                distance, dihedral_dist, denominator, weight_gradient, distance_vec, _, _, dw_da_dx = self.cartesian_distance(data_point)
                # if i < 2:
                #     num_grad, num_grad_dw_da = numerical_gradient(self.cartesian_distance, self.molecule.get_coordinates_in_bohr(), data_point)

                #     print('Grad_diff: ', i, '\n', num_grad, weight_gradient, num_grad - weight_gradient)
                    
                #     print('dw da grad', num_grad_dw_da, dw_da_dx, num_grad_dw_da - dw_da_dx)
                
                
                if abs(distance) < min_distance:
                    min_distance = abs(distance)
                distances_and_gradients.append((distance, dihedral_dist, i, denominator, weight_gradient, distance_vec))

        elif self.weightfunction_type == 'cartesian-hessian':
            for i, data_point in enumerate(self.qm_data_points[:]):
                
                distance, dihedral_dist, denominator, weight_gradient, distance_vec, _, _, dw_da_dx = self.cartesian_hessian_distance(data_point)
                
                # if i == 1:
                    
                    # num_grad, num_grad_dw_da = numerical_gradient(self.cartesian_hessian_distance, self.molecule.get_coordinates_in_bohr(), data_point)

                    # print('Grad_diff: ', i, '\n', num_grad, weight_gradient, num_grad - weight_gradient)
                    
                    # print('dw da grad', num_grad_dw_da, dw_da_dx, num_grad_dw_da - dw_da_dx)
                if abs(distance) < min_distance:
                    min_distance = abs(distance)
                distances_and_gradients.append((distance, dihedral_dist, i, denominator, weight_gradient, distance_vec))
        
        close_distances = None
        close_distances = [
            (self.qm_data_points[index], distance, dihedral_dist, denom, wg, distance_vec, index) 
            for distance, dihedral_dist, index, denom, wg, distance_vec in distances_and_gradients 
            if abs(distance) <= min_distance + self.distance_thrsh + 1000]

        for qm_data_point, distance, dihedral_dist, denominator_cart, weight_grad_cart, distance_vector, label_idx in close_distances:
            
            weight_cart = 1.0 / (denominator_cart)
        
            sum_weights_cart += weight_cart
            sum_weight_gradients_cart += weight_grad_cart

            potential, gradient_mw, r_i = self.compute_potential(qm_data_point, self.impes_coordinate.internal_coordinates_values)
            gradient = sqrt_masses * gradient_mw.reshape(-1)
            gradient = gradient.reshape(gradient_mw.shape)
            gradients.append(gradient)
            averaged_int_dists.append(qm_data_point.internal_coordinates_values)

            hessian_error.append(r_i)
            potentials.append(potential)
            used_labels.append(label_idx)
            weights_cart.append(weight_cart)
            self.potentials.append(potential)
            self.gradients.append(gradient)
            weight_gradients_cart.append(weight_grad_cart)
           
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

        self.impes_coordinate.gradient = (np.tensordot(W_i, gradients, axes=1) + np.tensordot(potentials, grad_W_i, axes=1))

        # --- 4.  book-keeping (optional) ---------------------------------------------
        for lbl, Wi in zip(used_labels, W_i):
            self.weights[lbl] = Wi
        # print('weights', self.weights)
        # self.sum_of_weights      = W_i.sum()          # if you really need it later
        self.averaged_int_dist   = np.tensordot(W_i, averaged_int_dists, axes=1)

    def shepard_interpolation_paral(self, max_workers=None):
        """Performs Shepard interpolation with parallelized datapoint-distance evaluation.

        This is intended as a comparison path against ``shepard_interpolation``.
        """

        natms = self.impes_coordinate.cartesian_coordinates.shape[0]

        sum_weights_cart = 0.0
        potentials = []
        gradients = []
        hessian_error = []
        self.weights = {}
        weights_cart = []
        averaged_int_dists = []
        weight_gradients_cart = []
        used_labels = []
        self.potentials = []
        self.gradients = []
        self.dw_dalpha_list = []
        self.dw_dX_dalpha_list = []

        masses = self.molecule.get_masses().copy()
        masses_cart = np.repeat(masses, 3)
        sqrt_masses = np.sqrt(masses_cart)

        sum_weight_gradients_cart = np.zeros((natms, 3))
        self.time_step_reducer = False
        self.calc_optim_trust_radius = True

        if self.external_weights is not None:
            for i, data_point in enumerate(self.qm_data_points[:]):
                data_point.confidence_radius = self.external_weights[i]

        def evaluate_datapoint(index_data_point):
            i, data_point = index_data_point
            if self.weightfunction_type == 'cartesian':
                distance, dihedral_dist, denominator, weight_gradient, distance_vec, _, dw_dalpha_i, dw_dX_dalpha_i = self.cartesian_distance(
                    data_point, store_alpha_gradients=False)
            elif self.weightfunction_type == 'cartesian-hessian':
                distance, dihedral_dist, denominator, weight_gradient, distance_vec, _, dw_dalpha_i, dw_dX_dalpha_i = self.cartesian_hessian_distance(
                    data_point, store_alpha_gradients=False)
            else:
                errtxt = "Unrecognized weight function type: "
                errtxt += self.weightfunction_type
                raise ValueError(errtxt)

            return (i, distance, dihedral_dist, denominator, weight_gradient,
                    distance_vec, dw_dalpha_i, dw_dX_dalpha_i)

        def collect_distances_and_gradients():
            min_distance = float('inf')
            distances_and_gradients = []

            datapoints = list(enumerate(self.qm_data_points[:]))
            n_points = len(datapoints)

            workers = max_workers
            if workers is None:
                workers = os.cpu_count() or 1
            workers = max(1, min(workers, n_points))

            if workers == 1 or n_points < 2:
                results = [evaluate_datapoint(item) for item in datapoints]
            else:
                with ThreadPoolExecutor(max_workers=workers) as executor:
                    results = list(executor.map(evaluate_datapoint, datapoints))

            for i, distance, dihedral_dist, denominator, weight_gradient, distance_vec, dw_dalpha_i, dw_dX_dalpha_i in results:
                if abs(distance) < min_distance:
                    min_distance = abs(distance)
                distances_and_gradients.append((distance, dihedral_dist, i, denominator, weight_gradient, distance_vec))

                if self.calc_optim_trust_radius:
                    self.dw_dalpha_list.append(dw_dalpha_i)
                    if self.weightfunction_type == 'cartesian':
                        dw_dX_dalpha_i_full = np.zeros((natms, 3))
                        dw_dX_dalpha_i_full[self.symmetry_information[3]] = dw_dX_dalpha_i
                        self.dw_dX_dalpha_list.append(dw_dX_dalpha_i_full)
                    else:
                        self.dw_dX_dalpha_list.append(dw_dX_dalpha_i)

            return min_distance, distances_and_gradients

        min_distance, distances_and_gradients = collect_distances_and_gradients()

        close_distances = [
            (self.qm_data_points[index], distance, dihedral_dist, denom, wg, distance_vec, index)
            for distance, dihedral_dist, index, denom, wg, distance_vec in distances_and_gradients
            if abs(distance) <= min_distance + self.distance_thrsh + 1000]

        for qm_data_point, distance, dihedral_dist, denominator_cart, weight_grad_cart, distance_vector, label_idx in close_distances:
            weight_cart = 1.0 / (denominator_cart)

            sum_weights_cart += weight_cart
            sum_weight_gradients_cart += weight_grad_cart

            potential, gradient_mw, r_i = self.compute_potential_paral(
                qm_data_point,
                self.impes_coordinate.internal_coordinates_values,
                max_workers=max_workers)
            gradient = sqrt_masses * gradient_mw.reshape(-1)
            gradient = gradient.reshape(gradient_mw.shape)
            gradients.append(gradient)
            averaged_int_dists.append(qm_data_point.internal_coordinates_values)

            hessian_error.append(r_i)
            potentials.append(potential)
            used_labels.append(label_idx)
            weights_cart.append(weight_cart)
            self.potentials.append(potential)
            self.gradients.append(gradient)
            weight_gradients_cart.append(weight_grad_cart)

        self.impes_coordinate.energy = 0.0
        self.impes_coordinate.gradient = np.zeros((natms, 3))
        self.impes_coordinate.NAC = np.zeros((natms, 3))

        w_i = np.array(weights_cart, dtype=np.float64)
        grad_w_i = np.array(weight_gradients_cart, dtype=np.float64)

        S = w_i.sum()
        sum_grad_w = grad_w_i.sum(axis=0)
        self.sum_of_weights = S
        self.sum_of_weights_grad = sum_grad_w

        W_i = w_i / S
        grad_W_i = (grad_w_i * S - w_i[:, None, None] * sum_grad_w) / S**2

        potentials = np.array(potentials, dtype=np.float64)
        gradients = np.array(gradients, dtype=np.float64)

        self.impes_coordinate.energy = np.dot(W_i, potentials)
        self.impes_coordinate.gradient = (
            np.tensordot(W_i, gradients, axes=1) +
            np.tensordot(potentials, grad_W_i, axes=1))

        for lbl, Wi in zip(used_labels, W_i):
            self.weights[lbl] = Wi
        self.averaged_int_dist = np.tensordot(W_i, averaged_int_dists, axes=1)

    def validate_shepard_parallel(
        self,
        molecule=None,
        max_workers=None,
        atol=1.0e-10,
        rtol=1.0e-8,
        restore_serial=True,
    ):
        """Validate serial vs parallel Shepard interpolation outputs.

        Returns a dictionary with per-quantity differences and pass/fail flags.
        """

        def _copy_if_array(value):
            return value.copy() if isinstance(value, np.ndarray) else value

        def _copy_list_maybe_arrays(values):
            return [_copy_if_array(v) for v in values]

        if molecule is not None:
            self.molecule = molecule
            self.distance_molecule = Molecule.from_xyz_string(molecule.get_xyz_string())
            self.define_impes_coordinate(molecule.get_coordinates_in_bohr())

        if self.qm_data_points is None:
            self.qm_data_points = self.read_qm_data_points()

        if self.interpolation_type != 'shepard':
            raise ValueError("validate_shepard_parallel requires interpolation_type='shepard'.")

        self.shepard_interpolation()

        serial_energy = float(self.impes_coordinate.energy)
        serial_gradient = self.impes_coordinate.gradient.copy()
        serial_nac = self.impes_coordinate.NAC.copy() if self.impes_coordinate.NAC is not None else None
        serial_weights = dict(self.weights)
        serial_sum_weights = float(self.sum_of_weights)
        serial_sum_weights_grad = self.sum_of_weights_grad.copy()
        serial_averaged_int_dist = self.averaged_int_dist.copy()
        serial_potentials = _copy_list_maybe_arrays(self.potentials)
        serial_gradients = _copy_list_maybe_arrays(self.gradients)
        serial_dw_dalpha_list = _copy_list_maybe_arrays(self.dw_dalpha_list)
        serial_dw_dX_dalpha_list = _copy_list_maybe_arrays(self.dw_dX_dalpha_list)

        self.shepard_interpolation_paral(max_workers=max_workers)

        parallel_energy = float(self.impes_coordinate.energy)
        parallel_gradient = self.impes_coordinate.gradient.copy()
        parallel_weights = dict(self.weights)
        parallel_sum_weights = float(self.sum_of_weights)
        parallel_sum_weights_grad = self.sum_of_weights_grad.copy()
        parallel_averaged_int_dist = self.averaged_int_dist.copy()

        all_weight_keys = sorted(set(serial_weights.keys()) | set(parallel_weights.keys()))
        serial_weight_values = np.array([serial_weights.get(k, 0.0) for k in all_weight_keys], dtype=np.float64)
        parallel_weight_values = np.array([parallel_weights.get(k, 0.0) for k in all_weight_keys], dtype=np.float64)

        energy_abs_diff = abs(serial_energy - parallel_energy)
        gradient_max_abs_diff = np.max(np.abs(serial_gradient - parallel_gradient))
        sum_weights_abs_diff = abs(serial_sum_weights - parallel_sum_weights)
        sum_weights_grad_max_abs_diff = np.max(np.abs(serial_sum_weights_grad - parallel_sum_weights_grad))
        averaged_int_dist_max_abs_diff = np.max(np.abs(serial_averaged_int_dist - parallel_averaged_int_dist))
        weights_max_abs_diff = (
            np.max(np.abs(serial_weight_values - parallel_weight_values))
            if serial_weight_values.size > 0 else 0.0
        )

        validation = {
            'energy_close': bool(np.isclose(serial_energy, parallel_energy, atol=atol, rtol=rtol)),
            'gradient_close': bool(np.allclose(serial_gradient, parallel_gradient, atol=atol, rtol=rtol)),
            'sum_weights_close': bool(np.isclose(serial_sum_weights, parallel_sum_weights, atol=atol, rtol=rtol)),
            'sum_weights_grad_close': bool(np.allclose(serial_sum_weights_grad, parallel_sum_weights_grad, atol=atol, rtol=rtol)),
            'averaged_int_dist_close': bool(np.allclose(serial_averaged_int_dist, parallel_averaged_int_dist, atol=atol, rtol=rtol)),
            'weights_close': bool(np.allclose(serial_weight_values, parallel_weight_values, atol=atol, rtol=rtol)),
            'energy_abs_diff': float(energy_abs_diff),
            'gradient_max_abs_diff': float(gradient_max_abs_diff),
            'sum_weights_abs_diff': float(sum_weights_abs_diff),
            'sum_weights_grad_max_abs_diff': float(sum_weights_grad_max_abs_diff),
            'averaged_int_dist_max_abs_diff': float(averaged_int_dist_max_abs_diff),
            'weights_max_abs_diff': float(weights_max_abs_diff),
            'n_weight_entries': int(serial_weight_values.size),
        }
        validation['all_close'] = bool(
            validation['energy_close'] and
            validation['gradient_close'] and
            validation['sum_weights_close'] and
            validation['sum_weights_grad_close'] and
            validation['averaged_int_dist_close'] and
            validation['weights_close']
        )

        if restore_serial:
            self.impes_coordinate.energy = serial_energy
            self.impes_coordinate.gradient = serial_gradient
            self.impes_coordinate.NAC = serial_nac
            self.weights = serial_weights
            self.sum_of_weights = serial_sum_weights
            self.sum_of_weights_grad = serial_sum_weights_grad
            self.averaged_int_dist = serial_averaged_int_dist
            self.potentials = serial_potentials
            self.gradients = serial_gradients
            self.dw_dalpha_list = serial_dw_dalpha_list
            self.dw_dX_dalpha_list = serial_dw_dX_dalpha_list

        return validation


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

            symmetry_weights_cos = []
            symmetry_weight_gradients_cos = []
            r_cut = 1.5

            for idx, symmetry_data_point in enumerate(symmetry_data_points):
                # if i == 0:
                #     symmetry_data_point.confidence_radius = 1.5
                energy = symmetry_data_point.energy
                masks = symmetry_data_point.mapping_masks


                # dihedral_maps = []
                # for key in self.symmetry_information[7][3]:

                #     if tuple(sorted(key[1:3])) not in dihedral_maps:
                        
                #         dihedral_maps.append(tuple(sorted(key[1:3])))

                # print(self.symmetry_information[7][3])
                # exit()
                for ww, mask in enumerate(masks[:]):
                    
                    sum_sym_3_dihedral_cos = {key: 0.0 for key in self.symmetry_information[7][3].keys()}
                    sum_sym_3_dihedral_prime_cos = {key :np.zeros_like(data_point.cartesian_coordinates.reshape(-1)) for key in self.symmetry_information[7][3].keys()}

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

                    dihedral_start = self.symmetry_information[-1][1]
                    if not self.use_cosine_dihedral:

                        dist_check[dihedral_start:] = np.sin(dist_org[dihedral_start:])
                    
                    dist_correlation[dihedral_start:] = dist_check[dihedral_start:]

                    for i, element in enumerate(self.impes_coordinate.z_matrix[dihedral_start:], start=dihedral_start): 

                        if element[1:3] in self.symmetry_information[7][3]:                                       
            
                            sum_sym_dihedral += self.te_weight_general(dist_org[i], 3)
                            sum_sym_dihedral_prime += self.te_weight_gradient_general(dist_org[i], self.impes_coordinate.b_matrix[i,:], 3)
                            dist_correlation[i] = 0.0
                            sum_sym_3_dihedral_cos[element[1:3]] += 2.0 * (1.0 - np.cos(dist_org[i]))
                            sum_sym_3_dihedral_prime_cos[element[1:3]] += 2.0 * np.sin( dist_org[i]) * self.impes_coordinate.b_matrix[i, :]

                        elif tuple(sorted(element)) in self.symmetry_information[7][2]:
                            
                            sum_sym_dihedral += self.te_weight_general(dist_org[i], 2)
                            sum_sym_dihedral_prime += self.te_weight_gradient_general(dist_org[i], self.impes_coordinate.b_matrix[i,:], 2)      
                            dist_correlation[i] = 0.0

                        # print('sum sym dihedral', sum_sym_dihedral)
                    
                    # print('\n\n')
                    self.bond_rmsd.append(np.sqrt(np.mean(np.sum((dist_org[:self.symmetry_information[-1][0]])**2))))
                    self.angle_rmsd.append(np.sqrt(np.mean(np.sum(dist_org[self.symmetry_information[-1][0]:self.symmetry_information[-1][1]]**2))))
                    self.dihedral_rmsd.append(np.sqrt(np.mean(np.sum(dist_correlation[self.symmetry_information[-1][1]:]**2))))
                    
                    combined_weights = []
                    combined_weights_derivative = []

                    for key, entry in sum_sym_3_dihedral_cos.items():
                        if entry >= r_cut**2:
                            # Outside the trust radius: weight and force are EXACTLY zero
                            combined_weights.append(0.0)
                            
                            combined_weights_derivative.append(np.zeros((natm, 3)))
                        else:
                            # Inside the trust radius: evaluate the smooth polynomial dome
                            decay_term = 1.0 - (entry / r_cut**2)
                            
                            combined_weights.append(decay_term**3)
                            
                            grad_factor = (-3.0 / r_cut**2) * (decay_term**2)
                            grad_matrix = grad_factor * sum_sym_3_dihedral_prime_cos[key].reshape(natm, 3)
                            
                            combined_weights_derivative.append(grad_matrix)

                    
                    W_total = np.prod(combined_weights)
                 
                    grad_W_total = np.zeros((natm, 3))

                    for I in range(len(sum_sym_3_dihedral_cos.keys())):
                        # Calculate the product of all scalar weights EXCEPT the current one (J != I)
                        # Using a generator expression handles any number of rotors automatically
                        prod_other_w = np.prod([combined_weights[J] for J in range(len(sum_sym_3_dihedral_cos.keys())) if J != I])
                        
                        # Add this rotor's contribution to the total multidimensional gradient
                        grad_W_total += combined_weights_derivative[I] * prod_other_w

                    # Now append W_total and grad_W_total to your global lists for the normalization step
                    symmetry_weights_cos.append(W_total)


                    symmetry_weight_gradients_cos.append(grad_W_total)

                    # symmetry_weights_cos.append((1.0 - sum_sym_dihedral_cos/r_cut**2)**3)
                    # symmetry_weight_gradients_cos.append((-3.0/r_cut**2 * (1.0 - sum_sym_dihedral_cos/r_cut**2)**2) * sum_sym_dihedral_prime_cos.reshape(natm, 3))

                    # if sum_sym_dihedral == 0.0:
                    #     continue

                    symmetry_weights.append(sum_sym_dihedral**4)


                    symmetry_weight_gradients.append(4 * sum_sym_dihedral**3 * sum_sym_dihedral_prime.reshape(natm, 3))

                    pes = (energy + np.matmul(dist_check.T, grad) +
                        0.5 * np.linalg.multi_dot([dist_check.T, hessian, dist_check]))

                    # hessian_error.append(np.sqrt(abs(np.linalg.multi_dot([dist_check.T, pinv_cut(hessian), dist_check]))))

                    potentials.append(pes)
                    dist_hessian = np.matmul(dist_check.T, hessian)

                    if not self.use_cosine_dihedral:

                        cos_dist = np.cos(dist_org[dihedral_start:])
                        grad[dihedral_start:] *= cos_dist
                        dist_hessian[dihedral_start:] *= cos_dist
                    
                    pes_prime = (np.matmul(self.impes_coordinate.b_matrix.T, (grad + dist_hessian))).reshape(natm, 3)
                    pot_gradients.append(pes_prime)


            w_i       = np.asarray(symmetry_weights_cos,          dtype=np.float64)       
            grad_w_i  = np.asarray(symmetry_weight_gradients_cos, dtype=np.float64)       

            S         = w_i.sum()                              
            sum_grad_w = grad_w_i.sum(axis=0)                  

            W_i       = w_i / S                             
            grad_W_i  = (grad_w_i * S - w_i[:, None, None] * sum_grad_w) / S**2    # (m, natm, 3)


            U_i       = np.asarray(potentials,        dtype=np.float64)            # energies     (m,)
            grad_U_i  = np.asarray(pot_gradients,     dtype=np.float64)            # gradients    (m, natm, 3)

            pes       = np.dot(W_i, U_i)                      # Σ Wᵢ Uᵢ

            grad_pes  = ( np.tensordot(W_i,      grad_U_i, axes=1) +   # Σ Wᵢ ∇Uᵢ
                        np.tensordot(U_i,   grad_W_i, axes=1) )       # Σ Uᵢ ∇Wᵢ
            # grad_pes shape: (natm, 3)

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
            
            dihedral_start = self.symmetry_information[-1][1]
            if not self.use_cosine_dihedral:
                dist_check[dihedral_start:] = np.sin(dist_org[dihedral_start:])
            

      
            self.bond_rmsd.append(np.sqrt(np.mean(np.sum((dist_org[:self.symmetry_information[-1][0]])**2))))
            self.angle_rmsd.append(np.sqrt(np.mean(np.sum(dist_org[self.symmetry_information[-1][0]:self.symmetry_information[-1][1]]**2))))
            self.dihedral_rmsd.append(np.sqrt(np.mean(np.sum(dist_check[self.symmetry_information[-1][1]:]**2))))

            pes = (energy + np.matmul(dist_check.T, grad) +
                        0.5 * np.linalg.multi_dot([dist_check.T, hessian, dist_check]))

            dist_hessian_eff = np.matmul(dist_check.T, hessian)


            if not self.use_cosine_dihedral:
                cos_dist = np.cos(dist_org[dihedral_start:])
                grad[dihedral_start:] *= cos_dist
                dist_hessian_eff[dihedral_start:] *= cos_dist


            pes_prime = (np.matmul(self.impes_coordinate.b_matrix.T, (grad + dist_hessian_eff))).reshape(natm, 3)

            return pes, pes_prime, (grad + dist_hessian_eff)

    def compute_potential_paral(self, data_point, org_int_coords, max_workers=None):
        """Parallelized potential evaluation for symmetry-expanded datapoints.

        Falls back to ``compute_potential`` when symmetry expansion is not active.
        """

        pes = 0.0
        gradient = None
        natm = data_point.cartesian_coordinates.shape[0]
        dp_label = data_point.point_label

        if len(self.qm_symmetry_data_points[dp_label]) > 1 and self.use_symmetry:

            symmetry_data_points = self.qm_symmetry_data_points[dp_label]

            symmetry_weights = []
            potentials = []
            pot_gradients = []
            symmetry_weight_gradients = []
            hessian_error = []

            symmetry_weights_cos = []
            symmetry_weight_gradients_cos = []
            r_cut = 1.5

            dihedral_start = self.symmetry_information[-1][1]
            bond_end = self.symmetry_information[-1][0]
            angle_end = self.symmetry_information[-1][1]
            sym3_keys = tuple(self.symmetry_information[7][3].keys())
            sym3_dict = self.symmetry_information[7][3]
            sym2_set = set(self.symmetry_information[7][2])

            def evaluate_symmetry_mask_task(task):
                symmetry_data_point, mask0, mask = task

                energy = symmetry_data_point.energy
                sum_sym_3_dihedral_cos = {key: 0.0 for key in sym3_keys}
                sum_sym_3_dihedral_prime_cos = {
                    key: np.zeros_like(data_point.cartesian_coordinates.reshape(-1))
                    for key in sym3_keys
                }

                grad = symmetry_data_point.internal_gradient.copy()
                grad[mask0] = grad[mask]

                hessian = symmetry_data_point.internal_hessian.copy()
                hessian[np.ix_(mask0, mask0)] = hessian[np.ix_(mask, mask)]

                dist_org = (org_int_coords.copy() -
                            symmetry_data_point.internal_coordinates_values[mask])
                dist_check = (org_int_coords.copy() -
                              symmetry_data_point.internal_coordinates_values[mask])
                dist_correlation = (org_int_coords.copy() -
                                    symmetry_data_point.internal_coordinates_values[mask])

                sum_sym_dihedral = 0.0
                sum_sym_dihedral_prime = np.zeros_like(
                    data_point.cartesian_coordinates.reshape(-1))

                dist_check[dihedral_start:] = np.sin(dist_org[dihedral_start:])
                dist_correlation[dihedral_start:] = dist_check[dihedral_start:]

                for i, element in enumerate(self.impes_coordinate.z_matrix[dihedral_start:],
                                            start=dihedral_start):

                    if element[1:3] in sym3_dict:
                        sum_sym_dihedral += self.te_weight_general(dist_org[i], 3)
                        sum_sym_dihedral_prime += self.te_weight_gradient_general(
                            dist_org[i], self.impes_coordinate.b_matrix[i, :], 3)
                        dist_correlation[i] = 0.0
                        sum_sym_3_dihedral_cos[element[1:3]] += 2.0 * (1.0 - np.cos(dist_org[i]))
                        sum_sym_3_dihedral_prime_cos[element[1:3]] += (
                            2.0 * np.sin(dist_org[i]) * self.impes_coordinate.b_matrix[i, :])

                    elif tuple(sorted(element)) in sym2_set:
                        sum_sym_dihedral += self.te_weight_general(dist_org[i], 2)
                        sum_sym_dihedral_prime += self.te_weight_gradient_general(
                            dist_org[i], self.impes_coordinate.b_matrix[i, :], 2)
                        dist_correlation[i] = 0.0

                bond_rmsd = np.sqrt(np.mean(np.sum((dist_org[:bond_end])**2)))
                angle_rmsd = np.sqrt(np.mean(np.sum(dist_org[bond_end:angle_end]**2)))
                dihedral_rmsd = np.sqrt(np.mean(np.sum(dist_correlation[dihedral_start:]**2)))

                combined_weights = []
                combined_weights_derivative = []
                for key, entry in sum_sym_3_dihedral_cos.items():
                    if entry >= r_cut**2:
                        combined_weights.append(0.0)
                        combined_weights_derivative.append(np.zeros((natm, 3)))
                    else:
                        decay_term = 1.0 - (entry / r_cut**2)
                        combined_weights.append(decay_term**3)
                        grad_factor = (-3.0 / r_cut**2) * (decay_term**2)
                        grad_matrix = grad_factor * sum_sym_3_dihedral_prime_cos[key].reshape(natm, 3)
                        combined_weights_derivative.append(grad_matrix)

                W_total = np.prod(combined_weights)
                grad_W_total = np.zeros((natm, 3))
                for I in range(len(sum_sym_3_dihedral_cos.keys())):
                    prod_other_w = np.prod(
                        [combined_weights[J] for J in range(len(sum_sym_3_dihedral_cos.keys()))
                         if J != I])
                    grad_W_total += combined_weights_derivative[I] * prod_other_w

                symmetry_weight = sum_sym_dihedral**4
                symmetry_weight_gradient = (
                    4 * sum_sym_dihedral**3 * sum_sym_dihedral_prime.reshape(natm, 3))

                pes_local = (energy + np.matmul(dist_check.T, grad) +
                             0.5 * np.linalg.multi_dot([dist_check.T, hessian, dist_check]))

                dist_hessian = np.matmul(dist_check.T, hessian)
                cos_dist = np.cos(dist_org[dihedral_start:])
                grad[dihedral_start:] *= cos_dist
                dist_hessian[dihedral_start:] *= cos_dist

                pes_prime = (
                    np.matmul(self.impes_coordinate.b_matrix.T, (grad + dist_hessian))).reshape(natm, 3)

                return (W_total, grad_W_total, pes_local, pes_prime,
                        bond_rmsd, angle_rmsd, dihedral_rmsd,
                        symmetry_weight, symmetry_weight_gradient)

            def evaluate_all_symmetry_tasks():
                tasks = []
                for symmetry_data_point in symmetry_data_points:
                    masks = symmetry_data_point.mapping_masks
                    for mask in masks[:]:
                        tasks.append((symmetry_data_point, masks[0], mask))

                n_tasks = len(tasks)
                workers = max_workers
                if workers is None:
                    workers = os.cpu_count() or 1
                workers = max(1, min(workers, n_tasks))

                if workers == 1 or n_tasks < 2:
                    return [evaluate_symmetry_mask_task(task) for task in tasks]

                with ThreadPoolExecutor(max_workers=workers) as executor:
                    return list(executor.map(evaluate_symmetry_mask_task, tasks))

            task_results = evaluate_all_symmetry_tasks()

            for (W_total, grad_W_total, pes_local, pes_prime, bond_rmsd, angle_rmsd,
                 dihedral_rmsd, symmetry_weight, symmetry_weight_gradient) in task_results:
                symmetry_weights_cos.append(W_total)
                symmetry_weight_gradients_cos.append(grad_W_total)
                symmetry_weights.append(symmetry_weight)
                symmetry_weight_gradients.append(symmetry_weight_gradient)
                potentials.append(pes_local)
                pot_gradients.append(pes_prime)
                self.bond_rmsd.append(bond_rmsd)
                self.angle_rmsd.append(angle_rmsd)
                self.dihedral_rmsd.append(dihedral_rmsd)

            w_i = np.asarray(symmetry_weights_cos, dtype=np.float64)
            grad_w_i = np.asarray(symmetry_weight_gradients_cos, dtype=np.float64)

            S = w_i.sum()
            sum_grad_w = grad_w_i.sum(axis=0)

            W_i = w_i / S
            grad_W_i = (grad_w_i * S - w_i[:, None, None] * sum_grad_w) / S**2

            U_i = np.asarray(potentials, dtype=np.float64)
            grad_U_i = np.asarray(pot_gradients, dtype=np.float64)

            pes = np.dot(W_i, U_i)
            grad_pes = (np.tensordot(W_i, grad_U_i, axes=1) +
                        np.tensordot(U_i, grad_W_i, axes=1))

            for j, Wi in enumerate(W_i):
                self.bond_rmsd[j] *= Wi
                self.angle_rmsd[j] *= Wi
                self.dihedral_rmsd[j] *= Wi

            gradient = grad_pes
            return pes, gradient, hessian_error

        return self.compute_potential(data_point, org_int_coords)
            

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

    
    def trust_radius_weight_gradient(self, confidence_radius, distance):
        

        denominator = (
                (distance / confidence_radius)**(2 * self.exponent_p) +
                (distance / confidence_radius)**(2 * self.exponent_q))
        trust_radius_weight_gradient = -1.0 * ((( -2.0 * self.exponent_p * ((distance / confidence_radius)**(2 * self.exponent_p)) / confidence_radius) - 
                                             (2.0 * self.exponent_q * ((distance / confidence_radius)**(2 * self.exponent_q) / confidence_radius))) / denominator**2)
        return trust_radius_weight_gradient

    
    def trust_radius_weight_gradient_gradient(self, confidence_radius, distance, distance_vector):
        
        # confidence_radius = datapoint.confidence_radius
        # distance, _, _, _, distance_vector, _ = self.cartesian_distance(datapoint)
        denominator = (
                (distance / confidence_radius)**(2 * self.exponent_p) +
                (distance / confidence_radius)**(2 * self.exponent_q))
        
        
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
        
        weight_gradient = (-1.0 * (2.0 * (derivative_p + derivative_q)) *
                       (1.0 / (denominator**2)))
        
        
   
        return  denominator, weight_gradient
    

    def trust_radius_tc_weight_gradient(self, confidence_radius, distance, imp_int_coords_distance):
        
        combined_distance_sq = (distance**2 + imp_int_coords_distance)

        denominator = (
                (combined_distance_sq / confidence_radius**2)**(self.exponent_p) +
                (combined_distance_sq / confidence_radius**2)**(self.exponent_q))
        
        trust_radius_weight_gradient = -1.0 * (((-2.0 * self.exponent_p * ((combined_distance_sq / confidence_radius**2)**(self.exponent_p)) / confidence_radius) - 
                                             (2.0 * self.exponent_q * ((combined_distance_sq / confidence_radius**2)**(self.exponent_q) / confidence_radius))) / denominator**2)
        return trust_radius_weight_gradient

    
    def trust_radius_tc_weight_gradient_gradient(self, confidence_radius, distance, distance_vector, imp_int_coords_distance, imp_int_coords_dist_derivative):
        
        combined_distance_sq = (distance**2 + imp_int_coords_distance)
        
        denominator = (
                (combined_distance_sq / confidence_radius**2)**(self.exponent_p) +
                (combined_distance_sq / confidence_radius**2)**(self.exponent_q))
        
        
        trust_radius_weight_gradient_gradient_1_nominator = (((-2.0 * self.exponent_p**2 * (imp_int_coords_dist_derivative + 2.0 * distance_vector.reshape(-1)) * combined_distance_sq**(self.exponent_p - 1.0)) 
                                                /(confidence_radius**(2.0 * self.exponent_p) * confidence_radius))
                                                - (2.0 * self.exponent_q**2.0 * (imp_int_coords_dist_derivative + 2.0 * distance_vector.reshape(-1)) * combined_distance_sq**(self.exponent_q - 1.0)) 
                                                /(confidence_radius**(2.0 * self.exponent_q) * confidence_radius)).reshape(distance_vector.shape[0], 3)

        trust_radius_weight_gradient_gradient_2_nominator = 2.0 * ((((self.exponent_p * (imp_int_coords_dist_derivative + 2.0 * distance_vector.reshape(-1)) * combined_distance_sq**(self.exponent_p - 1.0)) /(confidence_radius**(2.0 * self.exponent_p)))
                                                + (self.exponent_q * (imp_int_coords_dist_derivative + 2.0 * distance_vector.reshape(-1)) * combined_distance_sq**(self.exponent_q - 1.0)) / (confidence_radius**(2.0 * self.exponent_q))) 
                                                * (((-2.0 * (combined_distance_sq/confidence_radius**2)**(self.exponent_p) * self.exponent_p) / confidence_radius) 
                                                - ((2.0 * (combined_distance_sq/confidence_radius**2)**(self.exponent_q) * self.exponent_q) / confidence_radius))).reshape(distance_vector.shape[0], 3)

        
        trust_radius_weight_gradient_gradient = -1.0 * trust_radius_weight_gradient_gradient_1_nominator / denominator**2 + trust_radius_weight_gradient_gradient_2_nominator / denominator**3


        return trust_radius_weight_gradient_gradient
    
    def YM_target_customized_shepard_weight_gradient(self, distance_vector, distance, confidence_radius, imp_int_coords_distance, imp_int_coordinate_derivative):
        """ Returns the derivative of an unormalized Shepard interpolation
            weight with respect to the Cartesian coordinates in
            self.impes_coordinate

            :param distance_vector:
                The Cartesian distance vector between
                the current data_point and self.impes_coordinate.
            :param distance:
                The norm of the distance vector * sqrt(N), N number of atoms.
        """

        combined_distance = np.sqrt(imp_int_coords_distance)
        combined_distance_sq = (distance**2 + imp_int_coords_distance)
        
        denominator = (
            ((combined_distance_sq) / confidence_radius**2)**(self.exponent_p) +
            ((combined_distance_sq) / confidence_radius**2)**(self.exponent_q))
        derivative_p = ((self.exponent_p * ((2.0 * distance_vector.reshape(-1) + imp_int_coordinate_derivative) ) *
                        combined_distance_sq**(self.exponent_p - 1) /
                        confidence_radius**(2 * self.exponent_p)))
        derivative_q = ((self.exponent_q * ((2.0 * distance_vector.reshape(-1) + imp_int_coordinate_derivative)) *
                        combined_distance_sq**(self.exponent_q - 1) /
                        confidence_radius**(2 * self.exponent_q)))
        

        
        weight_gradient = (-1.0 * ((derivative_p + derivative_q)) *
                       (1.0 / (denominator**2)))
        
   
        return  denominator, weight_gradient.reshape(distance_vector.shape[0], 3)
    
    def VL_target_customized_shepard_weight_gradient(self, distance_vector, distance, confidence_radius, dq, dq_dx):
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
        
        final_denominator = denominator * np.exp(self.alpha * dq) 

        final_weight_gradient = (weight_gradient.reshape(-1) * np.exp(-self.alpha * dq) - 1.0/denominator * self.alpha * np.exp(-self.alpha * dq) * dq_dx).reshape(-1, 3)

        return  final_denominator, final_weight_gradient
    
    def trust_radius_weight_gradient_hessian(self, confidence_radius, distance, imp_int_coords_distance):

        combined_distance_sq = (distance + imp_int_coords_distance)

        denominator = (
                (combined_distance_sq / confidence_radius**2)**(self.exponent_p) +
                (combined_distance_sq / confidence_radius**2)**(self.exponent_q))
        trust_radius_weight_gradient = -1.0 * ((( -2.0 * self.exponent_p * ((combined_distance_sq)**(self.exponent_p)) / (confidence_radius * confidence_radius**(2 * self.exponent_p))) - 
                                             (2.0 * self.exponent_q * ((combined_distance_sq)**(self.exponent_p)) / (confidence_radius * confidence_radius**(2 * self.exponent_p)))) / denominator**2)

        return trust_radius_weight_gradient
    
    def trust_radius_weight_gradient_gradient_hessian(self, confidence_radius, distance, grad_s, imp_int_coords_distance, imp_int_coords_distance_derivative):
        
        combined_distance_sq = (distance + imp_int_coords_distance)
        denominator = (
                (combined_distance_sq / confidence_radius**2)**(self.exponent_p) +
                (combined_distance_sq / confidence_radius**2)**(self.exponent_q))

        trust_radius_weight_gradient_gradient_nominator_1_1 = ((self.exponent_p * (combined_distance_sq)**(self.exponent_p - 1) / (confidence_radius**(2 * self.exponent_p))) + (self.exponent_q * (combined_distance_sq)**(self.exponent_q - 1) / (confidence_radius**(2 * self.exponent_q))))
        trust_radius_weight_gradient_gradient_nominator_1_2 = 2.0 * ((-2.0 * self.exponent_p * (combined_distance_sq)**(self.exponent_p) / (confidence_radius * confidence_radius**(2 * self.exponent_p))) - ( 2.0 * self.exponent_q * (combined_distance_sq)**(self.exponent_q) / (confidence_radius * confidence_radius**(2 * self.exponent_q))))

        trust_radius_weight_gradient_gradient_nominator_2 = ((-2.0 * self.exponent_p**2 * (combined_distance_sq)**(self.exponent_p - 1) / (confidence_radius * confidence_radius**(2 * self.exponent_p))) - (2.0 * self.exponent_q**2 * (combined_distance_sq)**(self.exponent_q - 1) / (confidence_radius * confidence_radius**(2 * self.exponent_q))))
        
        trust_radius_weight_gradient_gradient = (((trust_radius_weight_gradient_gradient_nominator_1_1 * trust_radius_weight_gradient_gradient_nominator_1_2 ) / denominator**3) - trust_radius_weight_gradient_gradient_nominator_2 / denominator**2) * (grad_s + imp_int_coords_distance_derivative.reshape(-1, 3))


        return trust_radius_weight_gradient_gradient
    
    def shepard_weight_gradient_hessian(self, quad_distance, confidence_radius, grad_s, imp_int_coords_distance=0, imp_int_coords_distance_derivative=0):
        """ Returns the derivative of an unormalized Shepard interpolation
            weight with respect to the Cartesian coordinates in
            self.impes_coordinate

            :param distance_vector:
                The Cartesian distance vector between
                the current data_point and self.impes_coordinate.
            :param distance:
                The norm of the distance vector * sqrt(N), N number of atoms.
        """

        combined_distance_sq = (quad_distance + imp_int_coords_distance)

        denominator = (
            ((combined_distance_sq) / confidence_radius**2)**(self.exponent_p) +
            ((combined_distance_sq) / confidence_radius**2)**(self.exponent_q))
        
        derivative_p = ((self.exponent_p * 
                        combined_distance_sq**(self.exponent_p - 1) /
                        confidence_radius**(2 * self.exponent_p)))
        derivative_q = ((self.exponent_q *
                        combined_distance_sq**(self.exponent_q - 1) /
                        confidence_radius**(2 * self.exponent_q)))
        

        
        weight_gradient = (-1.0 * (1.0 * (derivative_p + derivative_q) *
                       (1.0 / (denominator**2)))) * (grad_s + imp_int_coords_distance_derivative.reshape(-1, 3))
                    
        return  denominator, weight_gradient
    
    
    def calculate_translation_coordinates(self, given_coordinates, weights=None):
        """Center the molecule by translating its geometric center to (0, 0, 0)."""
        
        if weights is not None:
            w = np.asarray(w, dtype=float)
            W = np.sum(w)
            return given_coordinates - np.sum(given_coordinates * w[:, None], axis=0) / W

        center = np.mean(given_coordinates, axis=0)
        translated_coordinates = given_coordinates - center



        return translated_coordinates

    def cartesian_distance(self, data_point, store_alpha_gradients=True):
        """Calculates and returns the cartesian distance between
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.

           :param data_point:
                InterpolationDatapoint object
        """

        def numerical_gradient(f, X, dp, eps=1e-6):
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

                    f_plus = self.impes_coordinate.internal_coordinates_values[0]
                    # f_plus, _, _  = f(dp, self.impes_coordinate.internal_coordinates_values)
                    current_geometry_minus = X - dX
                    self.define_impes_coordinate(current_geometry_minus)
                    # f_minus, _, _  = f(dp, self.impes_coordinate.internal_coordinates_values)
                   
                    f_minus = self.impes_coordinate.internal_coordinates_values[0]

                    grad[a, mu] = (f_plus - f_minus) / (2 * eps) 
            
            return grad.reshape(-1)

        # First, translate the cartesian coordinates to zero
        target_coordinates = data_point.cartesian_coordinates.copy()
        reference_coordinates = self.impes_coordinate.cartesian_coordinates.copy()


        active_atoms = np.delete(np.arange(reference_coordinates.shape[0]), self.symmetry_information[4])

        target_coordinates_core = target_coordinates[active_atoms]
        reference_coordinates_core = reference_coordinates[active_atoms]

        center_target_coordinates_core = self.calculate_translation_coordinates(target_coordinates_core)
        center_reference_coordinates_core = self.calculate_translation_coordinates(reference_coordinates_core)

        distance = 0
        distance_vector_sub = 0
        grad_s = 1.0
        Xc = center_reference_coordinates_core     # variable (current geometry)
        Yc = center_target_coordinates_core        # fixed (datapoint)
        R_x  = geometric.rotate.get_rot(Yc, Xc)

        rotated_current = (R_x @ Yc.T).T
        distance_vector_sub =  Xc - rotated_current
        distance = np.linalg.norm(distance_vector_sub)

        # distance_vector_sub = np.zeros_like(reference_coordinates)   # (natms,3)
        
        # distance_vector_sub[self.symmetry_information[3]] = distance_vector_sub


        # grad_s = distance_vector_sub @ R_x
        # grad_s -= grad_s.mean(axis=0, keepdims=True)

        # distance_vector_sub = grad_s
        
        dihedral_dist = 0.0
        if distance < 1e-8:
            distance = 1e-8
            distance_vector_sub[:] = 0

        H_eff = data_point.internal_hessian
        eigval, V = np.linalg.eigh(H_eff)

        participation = V**2  # shape (n_coord, n_modes)

        # For each internal coordinate, find dominant eigenmode
        dominant_mode_per_coord = np.argmax(participation, axis=1)
        coord_curvature = np.sum(eigval * (V**2), axis=1)
        # print(self.z_matrix)
        # print('True decomposition', coord_curvature)
            
        # sum of the derivative necessary 
        imp_int_coord_derivative_contribution = np.zeros_like(distance_vector_sub).reshape(-1)
        imp_int_coord_distance = 0
        dq_dx = np.zeros_like(distance_vector_sub).reshape(-1)
        Dimp_sq = 0.0
        if self.use_tc_weights:
            for element in data_point.imp_int_coordinates:
                idx = self.z_matrix.index(element)
                org_imp_int_coord_distance = (self.impes_coordinate.internal_coordinates_values[idx] - data_point.internal_coordinates_values[idx])
                unmw_b_matrix = self.impes_coordinate.b_matrix[idx, :] / self.impes_coordinate.inv_sqrt_masses
                
                if len(element) == 4:
                    sin_delta = np.sin(org_imp_int_coord_distance)
                    cos_delta = np.cos(org_imp_int_coord_distance)
                    imp_int_coord_derivative_contribution += 2.0 * unmw_b_matrix * cos_delta * sin_delta
                    imp_int_coord_distance += sin_delta**2

                    dq_dx +=  coord_curvature[idx] * 2.0 * unmw_b_matrix * np.cos(coord_curvature[idx] * org_imp_int_coord_distance) * np.sin(coord_curvature[idx] * org_imp_int_coord_distance) 
                    Dimp_sq += np.sin(coord_curvature[idx] * org_imp_int_coord_distance)**2
                else:
                    imp_int_coord_distance += org_imp_int_coord_distance**2
                    imp_int_coord_derivative_contribution += 2.0 * unmw_b_matrix * org_imp_int_coord_distance

                    dq_dx += coord_curvature[idx]**2 * 2.0 * unmw_b_matrix * org_imp_int_coord_distance
                    Dimp_sq += (coord_curvature[idx] * org_imp_int_coord_distance)**2
            
            if imp_int_coord_distance < 1e-8:
                imp_int_coord_distance = 1e-8
                imp_int_coord_derivative_contribution[:] = 0
        
        
        dw_dalhpa_i = 0
        dw_dX_dalpha_i = 0

        if self.interpolation_type == 'shepard':
            # denominator, weight_gradient_sub = self.shepard_weight_gradient(
            #     distance_vector_sub, distance, data_point.confidence_radius)
            denominator_imp_coord, weight_gradient_sub_imp_coord = self.YM_target_customized_shepard_weight_gradient(
                distance_vector_sub, distance, data_point.confidence_radius, imp_int_coord_distance, imp_int_coord_derivative_contribution)
            
            if 1 == 1:
                denominator_imp_coord, weight_gradient_sub_imp_coord = self.VL_target_customized_shepard_weight_gradient(
                    distance_vector_sub, distance, data_point.confidence_radius, Dimp_sq, dq_dx)
            
            if self.calc_optim_trust_radius:
                # dw_dalhpa_i = self.trust_radius_weight_gradient(data_point.confidence_radius, distance)
                # dw_dX_dalpha_i = self.trust_radius_weight_gradient_gradient(data_point.confidence_radius, distance, distance_vector)

                dw_dalhpa_i = self.trust_radius_tc_weight_gradient(data_point.confidence_radius, distance, imp_int_coord_distance)
                dw_dX_dalpha_i = self.trust_radius_tc_weight_gradient_gradient(data_point.confidence_radius, distance, distance_vector_sub, imp_int_coord_distance, imp_int_coord_derivative_contribution)

                if store_alpha_gradients:
                    self.dw_dalpha_list.append(dw_dalhpa_i)
                    dw_dX_dalpha_i_full = np.zeros_like(reference_coordinates)
                    dw_dX_dalpha_i_full[self.symmetry_information[3]] = dw_dX_dalpha_i
                    self.dw_dX_dalpha_list.append(dw_dX_dalpha_i_full)
        elif self.interpolation_type == 'simple':
            weight_gradient = self.simple_weight_gradient(
                distance_vector_sub, distance)
        else:
            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)
        
        # print(weight_gradient_sub, weight_gradient_sub_imp_coord, denominator_imp_coord, denominator)
        weight_gradient = np.zeros_like(reference_coordinates)   # (natms,3)
        weight_gradient[self.symmetry_information[3]] = weight_gradient_sub_imp_coord

        return distance, dihedral_dist, denominator_imp_coord, weight_gradient, distance_vector_sub, grad_s, dw_dalhpa_i, dw_dX_dalpha_i


    
    def cartesian_hessian_distance(self, data_point, store_alpha_gradients=True):
        """Calculates and returns the cartesian distance between
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.

           :param data_point:
                InterpolationDatapoint object
        """

        def denominator_func(a, distance):
            denominator = (
                (distance / a**2)**(self.exponent_p) +
                (distance / a**2)**(self.exponent_q))
            return 1.0 / denominator

        def numerical_gradient(f, distance, dp, eps=1e-6):
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
            w_da = 0

            a = dp.confidence_radius
            da = eps
            a_p = a + da

            w_da_p = f(a_p, distance)
            # f_plus, _, _  = f(dp, self.impes_coordinate.internal_coordinates_values)
            a_m = a - da
 
            # f_minus, _, _  = f(dp, self.impes_coordinate.internal_coordinates_values)
            w_da_m = f(a_m, distance)

            w_da += (w_da_p - w_da_m) / (2 * eps)
            
            return w_da

        def rigid_body_projector(X, masses, tol=1e-12):
            """
            Projector P (3N,3N) that removes translations + rotations in MASS-WEIGHTED coordinates.

            X:      (N,3) coordinates (use datapoint geometry, already centered is fine)
            masses: (N,)  atomic masses for these N atoms
            """
            X = np.asarray(X, float)
            m = np.asarray(masses, float)
            N = X.shape[0]
            n = 3*N

            # Use COM as origin (if you already centered by centroid, that's still OK:
            # changing origin only mixes rotation vectors with translations, which remain in the same rigid subspace)
            com = (m[:, None] * X).sum(axis=0) / m.sum()
            R = X - com

            # Build rigid basis B in mass-weighted coordinate space
            B = np.zeros((n, 6), float)
            mw = np.sqrt(np.repeat(m, 3))

            # translations
            for a in range(3):
                v = np.zeros((N, 3), float)
                v[:, a] = 1.0
                B[:, a] = (v.reshape(-1) * mw)

            # rotations: omega x r_i
            for k, omega in enumerate(np.eye(3)):
                v = np.cross(omega, R).reshape(-1)
                B[:, 3+k] = v * mw

            # Orthonormalize and detect rank (linear geometries -> rank 5)
            U, s, _ = np.linalg.svd(B, full_matrices=False)
            r = int(np.sum(s > tol * s[0])) if s[0] > 0 else 0
            Q = U[:, :r]

            P = np.eye(n) - Q @ Q.T
            return 0.5*(P + P.T), Q, s, r


        def apply_Jc(A):   # A: (Nc,3)
            return A - A.mean(axis=0, keepdims=True)


        def make_psd_stiffness_metric_from_massweighted_hessian(K, P, r_rigid,
                                                       floor_rel=1e-6, floor_abs=None,
                                                       cap_rel=None):
            """
            K:       (3N,3N) mass-weighted Hessian block (symmetric)
            P:       rigid-body projector in mass-weighted space
            r_rigid: number of rigid DOFs removed (typically 6, sometimes 5)
            Returns: K_metric (3N,3N) PSD stiffness metric in mass-weighted space
            """
            K = 0.5*(K + K.T)

            # Project out rigid DOFs (keeps a nullspace of dimension r_rigid)
            Kp = P @ K @ P
            Kp = 0.5*(Kp + Kp.T)

            w, U = np.linalg.eigh(Kp)

            # Enforce the rigid nullspace explicitly: zero the smallest r_rigid eigenvalues by magnitude
            idx = np.argsort(np.abs(w))
            rigid_idx = idx[:r_rigid]
            vib_idx   = idx[r_rigid:]

            lam = np.abs(w)                 # stiffness magnitude
            lam[rigid_idx] = 0.0            # keep rigid at exactly zero

            # Floor only the vibrational subspace if desired
            if vib_idx.size > 0:
                if floor_abs is None:
                    # robust scale from vibrational magnitudes
                    base = np.median(lam[vib_idx])
                    floor_abs = max(floor_rel * base, 0.0)

                lam[vib_idx] = np.maximum(lam[vib_idx], floor_abs)

                if cap_rel is not None:
                    cap = cap_rel * np.percentile(lam[vib_idx], 95)
                    lam[vib_idx] = np.minimum(lam[vib_idx], cap)

            K_metric = (U * lam) @ U.T
            return 0.5*(K_metric + K_metric.T)

        
        # First, translate the cartesian coordinates to zero
        target_coordinates = data_point.cartesian_coordinates.copy()
        reference_coordinates = self.impes_coordinate.cartesian_coordinates.copy()

        active_atoms = np.delete(np.arange(reference_coordinates.shape[0]), self.symmetry_information[4])

        target_coordinates_core = target_coordinates[active_atoms]
        reference_coordinates_core = reference_coordinates[active_atoms]

        center_target_coordinates_core = self.calculate_translation_coordinates(target_coordinates_core)
        center_reference_coordinates_core = self.calculate_translation_coordinates(reference_coordinates_core)

        distance = 0
        distance_vector = 0
        grad_s = 1.0

    
        Xc = center_reference_coordinates_core     # variable (current geometry)
        Yc = center_target_coordinates_core        # fixed (datapoint)
        R_x  = geometric.rotate.get_rot(Xc, Yc)
        dU_x = geometric.rotate.get_rot_der(Xc, Yc)   
        
        rotated_current = (R_x @ Xc.T).T
        distance_vector =  rotated_current - Yc
        distance = np.linalg.norm(distance_vector) 
        grad_s = distance_vector @ R_x
        grad_s -= grad_s.mean(axis=0, keepdims=True)

        H_eff = data_point.hessian
    
        idx3 = np.concatenate([3*active_atoms[:,None] + np.array([0,1,2])], axis=1).ravel()
        H_SS = H_eff[np.ix_(idx3, idx3)] # (Nc,)
        Msqrt = np.sqrt(np.repeat(self.molecule.get_masses(), 3))               # (3Nc,)

        # Build rigid-body projector in MASS-WEIGHTED space (removes translations + rotations)
        Pmw, Q, s, r_rigid = rigid_body_projector(Yc, self.molecule.get_masses())

        
        H_spd = make_psd_stiffness_metric_from_massweighted_hessian(
            H_SS, Pmw, r_rigid=r_rigid,
            floor_rel=1e-6, floor_abs=None,
            cap_rel=None
        )
        # H_spd = make_spd_hessian(Hc_ul, floor_abs=1e-6)
        H_spd_nmw = (Msqrt[:,None] * H_spd) * Msqrt[None,:]
        w, U   = np.linalg.eigh(H_spd)
        w_nmw, U   = np.linalg.eigh(H_spd_nmw)

        alpha = 1.1       # Your scaling factor for the Hessian importance
        beta  = 1e-6      # Regularization: ensures non-zero distance for rigid modes/flat regions
        d_vec = distance_vector.reshape(-1)
        dim = d_vec.size
        Metric = (alpha * H_spd_nmw) + (beta * np.eye(dim))
        distance_2     = (d_vec @ ((Metric) @ d_vec))
        distance = np.sqrt(distance_2)

        if distance < 1e-8:
            distance = 1e-8
            distance_vector[:] = 0
        y = (Metric @ d_vec).reshape(-1, 3)      # (Nc,3)

        rigid = y @ R_x                           # (Nc,3)
        S = y.T @ Xc                                    # (3,3)

        resp = np.einsum('baij,ij->ba', dU_x, S, optimize=True)   # (Nc,3)

        z_mat = rigid + resp                            # (Nc,3)
        z_mat = apply_Jc(z_mat)                        # (Nc,3)

        grad_s_sub = 1.0 * z_mat / distance                            # (Nc,3) 

        grad_s = np.zeros_like(reference_coordinates)   # (natms,3)
        grad_s[self.symmetry_information[3]] = grad_s_sub

        # sum of the derivative necessary 
        imp_int_coord_derivative_contribution = np.zeros_like(distance_vector).reshape(-1)
        imp_int_coord_distance = 0
        dihedral_dist = 0
        
        H_eff = data_point.internal_hessian
        eigval, V = np.linalg.eigh(H_eff)

        participation = V**2  # shape (n_coord, n_modes)

        # For each internal coordinate, find dominant eigenmode
        dominant_mode_per_coord = np.argmax(participation, axis=1)
        coord_curvature = np.sum(eigval * (V**2), axis=1)
        
        dq = []
        dq_dx = np.zeros_like(distance_vector).reshape(-1)
        # if data_point.point_label == "point_3_rinv_dihedral":
        #     data_point.imp_int_coordinates = [(1,0,2,3)]
        
        if self.use_tc_weights:
            for element in data_point.imp_int_coordinates:
                idx = self.z_matrix.index(element)
                org_imp_int_coord_distance = (self.impes_coordinate.internal_coordinates_values[idx] - data_point.internal_coordinates_values[idx])
                unmw_b_matrix = self.impes_coordinate.b_matrix[idx, :] / self.impes_coordinate.inv_sqrt_masses
                
                if len(element) == 4:
                    dq_dx +=  coord_curvature[idx] * 2.0 * unmw_b_matrix * np.cos(coord_curvature[idx] * org_imp_int_coord_distance) * np.sin(coord_curvature[idx] * org_imp_int_coord_distance)
                    dq.append(np.sin(coord_curvature[idx] * org_imp_int_coord_distance))
                else:
                    dq_dx += coord_curvature[idx]**2 * 2.0 * unmw_b_matrix * org_imp_int_coord_distance
                    dq.append(coord_curvature[idx] * org_imp_int_coord_distance)

            

        Dimp_sq = float(np.dot(np.array(dq), np.array(dq)))

        if Dimp_sq < 1e-8:
            imp_int_coord_distance = 1e-8
            imp_int_coord_derivative_contribution[:] = 0
            Dimp_sq = 1e-8
            dq_dx[:] = 0

        dw_dalpha_i = 0
        dw_dX_dalpha_i = 0
        if self.interpolation_type == 'shepard':
            denominator, weight_gradient = self.shepard_weight_gradient_hessian(distance, data_point.confidence_radius, grad_s, imp_int_coord_distance, imp_int_coord_derivative_contribution)  
            
            
            
            weight_gradient = (weight_gradient.reshape(-1) * np.exp(-self.alpha * Dimp_sq) - 1.0 / denominator * self.alpha * np.exp(-self.alpha * Dimp_sq) * dq_dx).reshape(-1,3)
            denominator = denominator * np.exp(self.alpha * Dimp_sq)
            

            if self.calc_optim_trust_radius:
                dw_dalpha_i = self.trust_radius_weight_gradient_hessian(data_point.confidence_radius, distance, imp_int_coord_distance)
                dw_dX_dalpha_i = self.trust_radius_weight_gradient_gradient_hessian(data_point.confidence_radius, distance, grad_s, imp_int_coord_distance, imp_int_coord_derivative_contribution)

                if store_alpha_gradients:
                    self.dw_dalpha_list.append(dw_dalpha_i)
                    self.dw_dX_dalpha_list.append(dw_dX_dalpha_i)


            
        elif self.interpolation_type == 'simple':
            weight_gradient = self.simple_weight_gradient(distance_vector, distance)
        else:
            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)
        
        

        return distance, dihedral_dist, denominator, weight_gradient, distance_vector, grad_s, dw_dalpha_i, dw_dX_dalpha_i
    
    
    
    def determine_important_internal_coordinates(
        self,
        qm_energy,
        qm_gradient,
        molecule,
        z_matrix,
        datapoints,
        dihedral_diff_const: bool = True,
    ):
        """
        New acquisition-oriented selection scheme.

        Main idea:
        ----------
        1) Diagnose each important datapoint locally.
        2) Build per-coordinate scores from:
            - hybrid relevance score
            - source displacement score
            - error score
            - causal source score
        3) Aggregate across important datapoints.
        4) Build coupled candidate blocks instead of selecting only one coordinate.
        5) Rank blocks using:
            - coordinate score
            - pair/coupling score
            - support across datapoints
            - novelty / redundancy penalty
        6) Return:
            - primary_constraint: best anchor coordinate
            - candidate_constraints: coordinates in the best block
            - ranked_blocks: rich metadata for downstream use
        """

        # ------------------------------------------------------------------
        # Hyperparameters
        # ------------------------------------------------------------------
        selection_rule = "coverage"   # 'relative' | 'coverage' | 'topk'
        relative_threshold = 0.10
        coverage_mass = 0.80
        topk = None

        rho_cap = 5.0
        alpha_trust = 0.4
        use_hybrid_score = True

        # Aggregation weights for coordinate-level global score
        w_hybrid = 0.45
        w_source = 0.10
        w_error = 0.20
        w_causal = 0.25

        # Build blocks
        max_anchor_candidates = 8
        max_block_partners = 4
        min_partner_pair_frac = 0.15      # partner pair-score relative to best partner
        min_partner_coord_frac = 0.20     # partner coord-score relative to anchor

        # Block ranking
        pair_weight = 0.35
        support_weight = 0.40
        novelty_weight = 0.60
        redundancy_weight = 0.90

        novelty_sigma = 0.20              # subspace-distance scale (tune)
        support_top_frac = 0.20           # top fraction in a datapoint to count "support"
        top_err_for_causal = 4
        top_src_per_err = 8

        eps = 1e-12

        constraints_to_exclude = []
        per_dp_results = []

        # ------------------------------------------------------------------
        # Exclude near-linear angles from being chosen as direct constraints
        # ------------------------------------------------------------------
        for element in self.z_matrix:
            if len(element) != 3:
                continue

            current_angle = molecule.get_angle_in_degrees(
                (element[0] + 1, element[1] + 1, element[2] + 1)
            )
            if abs(current_angle) < 2.0 or abs(current_angle) > 178.0:
                constraints_to_exclude.append(tuple(int(x) for x in element))

        masses = molecule.get_masses().copy()
        masses_cart = np.repeat(masses, 3)
        sqrt_masses = 1.0 / np.sqrt(masses_cart)
        qm_gradient_mw = qm_gradient.reshape(-1) * sqrt_masses

        q_current = np.asarray(self.impes_coordinate.internal_coordinates_values, dtype=float)
        N = len(z_matrix)

        # ------------------------------------------------------------------
        # Per-datapoint diagnostic pass
        # ------------------------------------------------------------------
        for datapoint, total_weight_contribution in datapoints:
            q_dp = np.asarray(datapoint.internal_coordinates_values, dtype=float)

            # raw and effective internal differences
            dq_raw, dq_eff, chain = self._compute_internal_differences(
                q_current=q_current,
                q_ref=q_dp,
                z_matrix=z_matrix,
                symmetry_information=self.symmetry_information,
            )

            # model and qm gradients in internal coordinates
            pred_E, pred_G_mw, _ = self.compute_potential(
                datapoint, self.impes_coordinate.internal_coordinates_values
            )
            pred_im_G_int = self.transform_gradient_to_internal_coordinates(
                molecule, pred_G_mw, self.impes_coordinate.b_matrix
            )
            pred_qm_G_int = self.transform_gradient_to_internal_coordinates(
                molecule, qm_gradient_mw.reshape(qm_gradient.shape), self.impes_coordinate.b_matrix
            )

            # local diagnostic without printing
            diag_result = self.compute_internal_gradient_diagnostics(
                z_matrix=z_matrix,
                dq_raw=dq_raw,
                H=datapoint.internal_hessian,
                g0=datapoint.internal_gradient,
                g_qm=pred_qm_G_int,
            )

            rho = np.asarray(diag_result["rho"], dtype=float)
            abs_delta_g_err = np.abs(np.asarray(diag_result["delta_g_err"], dtype=float))
            abs_dq = np.abs(np.asarray(diag_result["dq_eff"], dtype=float))
            abs_Hdq = np.abs(np.asarray(diag_result["Hdq"], dtype=float))

            # bounded trust factor in [0,1)
            rho_capped = np.minimum(rho, rho_cap)
            trust_factor = rho_capped / (1.0 + rho_capped)

            # --------------------------------------------------------------
            # Your existing combined energy/gradient local importance score
            # --------------------------------------------------------------
            delta_E = abs(qm_energy - pred_E)
            delta_G = np.linalg.norm(pred_qm_G_int - pred_im_G_int)
            delta_g = pred_qm_G_int - pred_im_G_int

            partial_energies, partial_gradient = self._compute_partial_taylor_contributions(
                dq_eff=dq_eff,
                chain=chain,
                H=np.asarray(datapoint.internal_hessian, dtype=float),
                g0=np.asarray(datapoint.internal_gradient, dtype=float),
            )

            abs_pE = np.abs(partial_energies)
            abs_work = np.abs(delta_g) * np.abs(dq_eff)
            w_E_raw = abs_pE * abs_work
            w_E = w_E_raw / (w_E_raw.sum() + eps)
            e_i = delta_E * w_E

            abs_pG = np.abs(partial_gradient)
            sum_abs_pG = abs_pG.sum()
            if sum_abs_pG < eps:
                w_G = np.zeros_like(abs_pG)
                g_i = np.zeros_like(abs_pG)
            else:
                w_G = abs_pG / (sum_abs_pG + eps)
                g_i = delta_G * w_G

            E_kcal = e_i * hartree_in_kcalpermol()

            L_ref = 0.1
            G_as_energy = g_i * L_ref
            G_kcal = G_as_energy * hartree_in_kcalpermol()

            lambda_grad = 0.5
            score_i_kcal = (1.0 - lambda_grad) * E_kcal + lambda_grad * G_kcal
            score_i_kcal = np.maximum(score_i_kcal, 0.0)

            w_tot = score_i_kcal / (score_i_kcal.sum() + eps)

            if use_hybrid_score:
                hybrid_score = w_tot * (alpha_trust + (1.0 - alpha_trust) * (1.0 - trust_factor))
            else:
                hybrid_score = w_tot.copy()

            hybrid_score = self._normalize_nonnegative(hybrid_score)

            # --------------------------------------------------------------
            # Additional explicit source / error / causal scores
            # --------------------------------------------------------------
            source_score = self._normalize_nonnegative(abs_dq)
            error_score = self._normalize_nonnegative(abs_delta_g_err)

            causal_score, pair_score = self._compute_causal_scores(
                z_matrix=z_matrix,
                dq_eff=dq_eff,
                H=np.asarray(datapoint.internal_hessian, dtype=float),
                abs_err=abs_delta_g_err,
                top_err=top_err_for_causal,
                top_src=top_src_per_err,
            )

            # support mask: which coordinates are in the top part for this datapoint?
            dp_coord_score = (
                w_hybrid * hybrid_score
                + w_source * source_score
                + w_error * error_score
                + w_causal * causal_score
            )
            dp_coord_score = self._normalize_nonnegative(dp_coord_score)

            order_dp = np.argsort(-dp_coord_score)
            k_support = max(1, int(np.ceil(support_top_frac * N)))
            support_mask = np.zeros(N, dtype=float)
            support_mask[order_dp[:k_support]] = 1.0

            per_dp_results.append({
                "datapoint": datapoint,
                "dp_weight": float(total_weight_contribution),
                "dq_raw": dq_raw.copy(),
                "dq_eff": dq_eff.copy(),
                "chain": chain.copy(),
                "rho": rho.copy(),
                "hybrid_score": hybrid_score.copy(),
                "source_score": source_score.copy(),
                "error_score": error_score.copy(),
                "causal_score": causal_score.copy(),
                "pair_score": pair_score.copy(),
                "coord_score": dp_coord_score.copy(),
                "support_mask": support_mask.copy(),
                "coords": [tuple(int(x) for x in c) for c in z_matrix],
            })

        if not per_dp_results:
            return [], [], []

        # ------------------------------------------------------------------
        # Normalize datapoint weights and keep important datapoints only
        # ------------------------------------------------------------------
        dp_weights = np.array([r["dp_weight"] for r in per_dp_results], dtype=float)
        if dp_weights.sum() < eps:
            dp_weights[:] = 1.0
        dp_weights /= dp_weights.sum()

        order = np.argsort(-dp_weights)
        cum = np.cumsum(dp_weights[order])

        keep_mask = cum <= coverage_mass
        if keep_mask.size > 0:
            keep_mask[0] = True

        keep_indices = order[np.where(keep_mask)[0]]

        if keep_mask.size > 0 and cum[keep_mask][-1] < coverage_mass and len(keep_indices) < len(order):
            next_idx = order[len(keep_indices)]
            keep_indices = np.append(keep_indices, next_idx)

        # ------------------------------------------------------------------
        # Aggregate coordinate-level and pair-level information
        # ------------------------------------------------------------------
        global_hybrid = np.zeros(N, dtype=float)
        global_source = np.zeros(N, dtype=float)
        global_error = np.zeros(N, dtype=float)
        global_causal = np.zeros(N, dtype=float)
        global_support = np.zeros(N, dtype=float)
        global_rho = np.zeros(N, dtype=float)
        global_pair = np.zeros((N, N), dtype=float)

        for idx in keep_indices:
            wk = dp_weights[idx]
            r = per_dp_results[idx]

            global_hybrid += wk * r["hybrid_score"]
            global_source += wk * r["source_score"]
            global_error += wk * r["error_score"]
            global_causal += wk * r["causal_score"]
            global_support += wk * r["support_mask"]
            global_rho += wk * r["rho"]
            global_pair += wk * r["pair_score"]

        global_hybrid = self._normalize_nonnegative(global_hybrid)
        global_source = self._normalize_nonnegative(global_source)
        global_error = self._normalize_nonnegative(global_error)
        global_causal = self._normalize_nonnegative(global_causal)
        global_support = global_support / (global_support.max() + eps)

        global_coord_score = (
            w_hybrid * global_hybrid
            + w_source * global_source
            + w_error * global_error
            + w_causal * global_causal
        )
        global_coord_score = self._normalize_nonnegative(global_coord_score)

        # ------------------------------------------------------------------
        # Select anchor coordinates according to the chosen rule
        # ------------------------------------------------------------------
        sorted_idx = np.argsort(-global_coord_score)
        sorted_weights = global_coord_score[sorted_idx]

        selected_pos = []
        if selection_rule == "relative":
            if sorted_weights.size > 0:
                wmax = sorted_weights.max()
                keep = sorted_weights >= (relative_threshold * wmax)
                selected_pos = list(np.where(keep)[0])

        elif selection_rule == "coverage":
            cumw = np.cumsum(sorted_weights)
            keep = cumw <= coverage_mass
            if keep.size > 0:
                keep[0] = True
            selected_pos = list(np.where(keep)[0])

            while (
                selected_pos
                and np.sum(sorted_weights[selected_pos]) < coverage_mass
                and selected_pos[-1] + 1 < len(sorted_weights)
            ):
                selected_pos.append(selected_pos[-1] + 1)

        elif selection_rule == "topk":
            k = topk if (topk is not None) else max(1, int(0.25 * N))
            selected_pos = list(range(min(k, N)))

        else:
            raise ValueError(f"Unknown selection_rule: {selection_rule}")

        candidate_anchor_indices = [int(sorted_idx[pos]) for pos in selected_pos]
        candidate_anchor_indices = candidate_anchor_indices[:max_anchor_candidates]

        # ------------------------------------------------------------------
        # Build and rank coupled blocks
        # ------------------------------------------------------------------
        ranked_blocks = []
        seen_blocks = set()

        for anchor_idx in candidate_anchor_indices:
            anchor_coord = tuple(int(x) for x in z_matrix[anchor_idx])

            # skip excluded anchors
            if anchor_coord in constraints_to_exclude:
                continue

            # For symmetric torsions excluded upstream, also skip as anchors
            if (
                len(anchor_coord) == 4
                and anchor_coord[1:3] in self.symmetry_information[7][3]
            ):
                continue

            block_indices = self._build_candidate_block(
                anchor_idx=anchor_idx,
                global_coord_score=global_coord_score,
                global_pair=global_pair,
                z_matrix=z_matrix,
                constraints_to_exclude=constraints_to_exclude,
                max_partners=max_block_partners,
                min_partner_pair_frac=min_partner_pair_frac,
                min_partner_coord_frac=min_partner_coord_frac,
            )

            key = tuple(sorted(block_indices))
            if key in seen_blocks:
                continue
            seen_blocks.add(key)

            # block-level components
            coord_mass = float(np.sum(global_coord_score[block_indices]))
            support_mass = float(np.mean(global_support[block_indices]))

            pair_mass = 0.0
            for a in range(len(block_indices)):
                for b in range(a + 1, len(block_indices)):
                    i = block_indices[a]
                    j = block_indices[b]
                    pair_mass += global_pair[i, j]

            dmin = self._min_block_distance_to_existing_datapoints(
                block_indices=block_indices,
                q_current=q_current,
                datapoints=[r["datapoint"] for r in per_dp_results],
                z_matrix=z_matrix,
                symmetry_information=self.symmetry_information,
            )

            novelty = dmin / (dmin + novelty_sigma + eps)
            redundancy = np.exp(- (dmin / (novelty_sigma + eps)) ** 2)

            acq_score = (
                coord_mass
                + pair_weight * pair_mass
            ) * (
                1.0 + support_weight * support_mass
            ) * (
                1.0 + novelty_weight * novelty
            ) * (
                1.0 - redundancy_weight * redundancy
            )

            block_coords = [tuple(int(x) for x in z_matrix[i]) for i in block_indices]

            ranked_blocks.append({
                "anchor_idx": int(anchor_idx),
                "anchor_coord": anchor_coord,
                "anchor_type": self.coord_type(anchor_coord),
                "block_indices": [int(i) for i in block_indices],
                "block_coords": block_coords,
                "coord_mass": float(coord_mass),
                "pair_mass": float(pair_mass),
                "support_mass": float(support_mass),
                "dmin": float(dmin),
                "novelty": float(novelty),
                "redundancy": float(redundancy),
                "acq_score": float(acq_score),
            })

        ranked_blocks = sorted(ranked_blocks, key=lambda x: x["acq_score"], reverse=True)

        if not ranked_blocks:
            return [], [], []

        # ------------------------------------------------------------------
        # Best block -> primary constraint and candidate constraints
        # ------------------------------------------------------------------
        best_block = ranked_blocks[0]

        # Prefer a torsion anchor inside the best block if it is competitive
        best_anchor_idx = best_block["anchor_idx"]
        best_anchor_score = global_coord_score[best_anchor_idx]

        torsion_candidates = []
        for idx in best_block["block_indices"]:
            coord = tuple(int(x) for x in z_matrix[idx])
            if self.coord_type(coord) == "dihedral":
                torsion_candidates.append(idx)

        if torsion_candidates:
            torsion_candidates = sorted(
                torsion_candidates,
                key=lambda i: global_coord_score[i],
                reverse=True
            )
            torsion_idx = torsion_candidates[0]
            if global_coord_score[torsion_idx] >= 0.6 * best_anchor_score:
                best_anchor_idx = torsion_idx

        primary_constraint = [tuple(int(x) for x in z_matrix[best_anchor_idx])]
        candidate_constraints = [tuple(int(x) for x in c) for c in best_block["block_coords"]]

        # ------------------------------------------------------------------
        # Debug printout
        # ------------------------------------------------------------------
        print("\n=== Global coordinate scores ===")
        print(f"{'rank':>4} {'idx':>4} {'type':>9} {'coord':>18} "
            f"{'score':>12} {'support':>12} {'rho':>12}")
        print("-" * 82)
        for rank, idx in enumerate(np.argsort(-global_coord_score)[:12]):
            coord = tuple(int(x) for x in z_matrix[idx])
            print(
                f"{rank:4d} {idx:4d} {self.coord_type(coord):>9} {self.format_coord(coord):>18} "
                f"{global_coord_score[idx]:12.6e} {global_support[idx]:12.6e} {global_rho[idx]:12.6e}"
            )

        print("\n=== Ranked candidate blocks ===")
        print(f"{'rank':>4} {'anchor':>18} {'type':>9} {'|block|':>8} "
            f"{'coord':>12} {'pair':>12} {'support':>12} "
            f"{'dmin':>12} {'novelty':>12} {'acq':>12}")
        print("-" * 118)
        for rank, blk in enumerate(ranked_blocks[:10]):
            print(
                f"{rank:4d} "
                f"{self.format_coord(blk['anchor_coord']):>18} "
                f"{blk['anchor_type']:>9} "
                f"{len(blk['block_indices']):8d} "
                f"{blk['coord_mass']:12.6e} "
                f"{blk['pair_mass']:12.6e} "
                f"{blk['support_mass']:12.6e} "
                f"{blk['dmin']:12.6e} "
                f"{blk['novelty']:12.6e} "
                f"{blk['acq_score']:12.6e}"
            )
            print(f"      block coords: {blk['block_coords']}")

        print("\nPrimary constraint:", primary_constraint)
        print("Candidate constraints (best block):", candidate_constraints)

        return primary_constraint, candidate_constraints, ranked_blocks


    # ======================================================================
    # Helper methods
    # ======================================================================

    def _normalize_nonnegative(self, x: np.ndarray, eps: float = 1e-12) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        x = np.maximum(x, 0.0)
        s = x.sum()
        if s < eps:
            return np.zeros_like(x)
        return x / s


    def _principal_torsion_delta(self, delta: float) -> float:
        """Wrap torsional difference to (-pi, pi]."""
        return (delta + np.pi) % (2.0 * np.pi) - np.pi


    def _compute_internal_differences(
        self,
        q_current: np.ndarray,
        q_ref: np.ndarray,
        z_matrix: List[Tuple[int, ...]],
        symmetry_information,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Returns
        -------
        dq_raw : raw internal differences
        dq_eff : effective differences used in your interpolation
        chain  : derivative chain factors (1 for non-torsions, cos(dq_raw) for torsions)
        """
        q_current = np.asarray(q_current, dtype=float)
        q_ref = np.asarray(q_ref, dtype=float)

        N = len(z_matrix)
        dq_raw = np.zeros(N, dtype=float)
        dq_eff = np.zeros(N, dtype=float)
        chain = np.ones(N, dtype=float)

        sym_torsion_centers = set()
        try:
            sym_torsion_centers = set(tuple(x) for x in symmetry_information[7][3])
        except Exception:
            sym_torsion_centers = set()

        for i, coord in enumerate(z_matrix):
            coord = tuple(int(x) for x in coord)

            if len(coord) == 4:
                # excluded/symmetric torsions -> zero effective displacement
                if tuple(coord[1:3]) in sym_torsion_centers:
                    dq_raw[i] = 0.0
                    dq_eff[i] = 0.0
                    chain[i] = 0.0
                    continue

                d = self._principal_torsion_delta(q_current[i] - q_ref[i])
                dq_raw[i] = d
                dq_eff[i] = np.sin(d)
                chain[i] = np.cos(d)
            else:
                d = q_current[i] - q_ref[i]
                dq_raw[i] = d
                dq_eff[i] = d
                chain[i] = 1.0

        return dq_raw, dq_eff, chain


    def _compute_partial_taylor_contributions(
        self,
        dq_eff: np.ndarray,
        chain: np.ndarray,
        H: np.ndarray,
        g0: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute per-coordinate split of Taylor energy and gradient contributions.
        Mirrors the logic you already had.
        """
        dq_eff = np.asarray(dq_eff, dtype=float)
        chain = np.asarray(chain, dtype=float)
        H = np.asarray(H, dtype=float)
        g0 = np.asarray(g0, dtype=float)

        N = len(dq_eff)
        partial_energies = np.zeros(N, dtype=float)
        partial_gradient = np.zeros(N, dtype=float)

        for i in range(N):
            partial_energies[i] += dq_eff[i] * g0[i]
            partial_gradient[i] += chain[i] * g0[i]

        for i in range(N):
            partial_energies[i] += 0.5 * dq_eff[i] * H[i, i] * dq_eff[i]
            partial_gradient[i] += chain[i] * H[i, i] * dq_eff[i]

            for j in range(i + 1, N):
                cross = dq_eff[i] * H[i, j] * dq_eff[j]
                cross_prime_ij = chain[i] * H[i, j] * dq_eff[j]
                cross_prime_ji = chain[j] * H[j, i] * dq_eff[i]

                partial_energies[i] += 0.5 * cross
                partial_energies[j] += 0.5 * cross
                partial_gradient[i] += cross_prime_ij
                partial_gradient[j] += cross_prime_ji

        return partial_energies, partial_gradient


    def _compute_causal_scores(
        self,
        z_matrix: List[Tuple[int, ...]],
        dq_eff: np.ndarray,
        H: np.ndarray,
        abs_err: np.ndarray,
        top_err: int = 4,
        top_src: int = 8,
        eps: float = 1e-12,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Build:
        causal_score[j] ~ how much coordinate j drives top failing responses
        pair_score[i,j] ~ coupling graph between response i and source j

        We use source terms |H[i,j] * dq_eff[j]|, weighted by top response error size.
        """
        dq_eff = np.asarray(dq_eff, dtype=float)
        H = np.asarray(H, dtype=float)
        abs_err = np.asarray(abs_err, dtype=float)

        N = len(z_matrix)
        causal_score = np.zeros(N, dtype=float)
        pair_score = np.zeros((N, N), dtype=float)

        err_order = np.argsort(-abs_err)
        top_err_idx = err_order[:min(top_err, N)]

        if abs_err.sum() < eps:
            return causal_score, pair_score

        err_weights = abs_err / (abs_err.sum() + eps)

        for i in top_err_idx:
            source_terms = np.abs(H[i, :] * dq_eff)
            src_order = np.argsort(-source_terms)
            src_keep = src_order[:min(top_src, N)]

            for j in src_keep:
                val = source_terms[j] * err_weights[i]
                causal_score[j] += val

                # symmetrize into a graph-like pair matrix
                pair_score[i, j] += val
                pair_score[j, i] += val

        causal_score = self._normalize_nonnegative(causal_score)

        # normalize pair matrix by max for easier downstream use
        pmax = np.max(pair_score)
        if pmax > eps:
            pair_score = pair_score / pmax

        return causal_score, pair_score


    def _build_candidate_block(
        self,
        anchor_idx: int,
        global_coord_score: np.ndarray,
        global_pair: np.ndarray,
        z_matrix: List[Tuple[int, ...]],
        constraints_to_exclude: List[Tuple[int, ...]],
        max_partners: int = 4,
        min_partner_pair_frac: float = 0.15,
        min_partner_coord_frac: float = 0.20,
    ) -> List[int]:
        """
        Build a coupled block around an anchor using:
        - strong pair/causal connectivity
        - non-negligible global coordinate score
        """
        N = len(z_matrix)
        anchor_score = global_coord_score[anchor_idx]

        pair_row = np.asarray(global_pair[anchor_idx], dtype=float).copy()
        partner_order = np.argsort(-pair_row)

        block = [int(anchor_idx)]

        # best pair value excluding self
        best_pair = 0.0
        for j in partner_order:
            if j == anchor_idx:
                continue
            best_pair = pair_row[j]
            break

        if best_pair <= 0.0:
            return block

        for j in partner_order:
            if j == anchor_idx:
                continue

            coord_j = tuple(int(x) for x in z_matrix[j])

            if coord_j in constraints_to_exclude:
                continue

            # skip excluded symmetric torsions as direct block members
            if len(coord_j) == 4 and coord_j[1:3] in self.symmetry_information[7][3]:
                continue

            pair_ok = pair_row[j] >= min_partner_pair_frac * best_pair
            coord_ok = global_coord_score[j] >= min_partner_coord_frac * anchor_score

            if pair_ok and coord_ok:
                block.append(int(j))

            if len(block) >= 1 + max_partners:
                break

        # small cleanup: sort block by descending coordinate score, keep anchor first
        partners = [i for i in block if i != anchor_idx]
        partners = sorted(partners, key=lambda i: global_coord_score[i], reverse=True)
        block = [int(anchor_idx)] + partners

        return block


    def _block_subspace_distance(
        self,
        block_indices: List[int],
        q_a: np.ndarray,
        q_b: np.ndarray,
        z_matrix: List[Tuple[int, ...]],
        symmetry_information,
    ) -> float:
        """
        Distance in the subspace of block_indices using your effective displacement definition.
        """
        q_a = np.asarray(q_a, dtype=float)
        q_b = np.asarray(q_b, dtype=float)

        d2 = 0.0
        sym_torsion_centers = set()
        try:
            sym_torsion_centers = set(tuple(x) for x in symmetry_information[7][3])
        except Exception:
            sym_torsion_centers = set()

        for i in block_indices:
            coord = tuple(int(x) for x in z_matrix[i])

            if len(coord) == 4:
                if tuple(coord[1:3]) in sym_torsion_centers:
                    continue
                d = self._principal_torsion_delta(q_a[i] - q_b[i])
                x = np.sin(d)
            else:
                x = q_a[i] - q_b[i]

            d2 += x * x

        return float(np.sqrt(d2))


    def _min_block_distance_to_existing_datapoints(
        self,
        block_indices: List[int],
        q_current: np.ndarray,
        datapoints: List[Any],
        z_matrix: List[Tuple[int, ...]],
        symmetry_information,
    ) -> float:
        """
        Minimal subspace distance between current structure and any existing datapoint
        in the candidate block.
        """
        if len(datapoints) == 0:
            return 1.0

        dmin = np.inf
        for dp in datapoints:
            q_dp = np.asarray(dp.internal_coordinates_values, dtype=float)
            d = self._block_subspace_distance(
                block_indices=block_indices,
                q_a=q_current,
                q_b=q_dp,
                z_matrix=z_matrix,
                symmetry_information=symmetry_information,
            )
            if d < dmin:
                dmin = d

        if not np.isfinite(dmin):
            dmin = 1.0

        return float(dmin)

    def coord_type(self, coord: Tuple[int, ...]) -> str:
        """Classify coordinate type from tuple length."""
        n = len(coord)
        if n == 2:
            return "bond"
        elif n == 3:
            return "angle"
        elif n == 4:
            return "dihedral"
        return "unknown"


    def format_coord(self, coord: Tuple[int, ...]) -> str:
            return str(tuple(int(x) for x in coord))


    def compute_internal_gradient_diagnostics(
        self,
        z_matrix,
        dq_raw,
        H,
        g0,
        g_qm,
        eps: float = 1e-12,
    ):
        """
        Compute local Taylor diagnostics in internal coordinates.

        dq_raw: raw internal difference q - q0 (wrapped to (-pi, pi] for torsions upstream)
        For torsions:
            dq_eff = sin(dq_raw)
            chain  = cos(dq_raw)
        For non-torsions:
            dq_eff = dq_raw
            chain  = 1

        IMPORTANT:
        Your interpolated internal gradient uses the mapping
            g_phi ≈ chain * ( g0 + H @ dq_eff )
        Therefore, the inverse-response operator in raw-gradient space is NOT H,
        but approximately:  H_eff = diag(chain) @ H @ diag(chain)
        (optionally plus a small diagonal correction term; see below).
        """
        dq_raw = np.asarray(dq_raw, dtype=float).reshape(-1)
        H = np.asarray(H, dtype=float)
        g0 = np.asarray(g0, dtype=float).reshape(-1)
        g_qm = np.asarray(g_qm, dtype=float).reshape(-1)

        N = len(z_matrix)
        if dq_raw.shape[0] != N:
            raise ValueError("dq_raw has wrong length")
        if g0.shape[0] != N:
            raise ValueError("g0 has wrong length")
        if g_qm.shape[0] != N:
            raise ValueError("g_qm has wrong length")
        if H.shape != (N, N):
            raise ValueError("H has wrong shape")

        dq_eff = dq_raw.copy()
        chain = np.ones(N)

        for i, coord in enumerate(z_matrix):
            if len(coord) == 4:
                dq_eff[i] = np.sin(dq_raw[i])
                chain[i] = np.cos(dq_raw[i])

        # Forward (matches your PES code)
        Hdq_raw = H @ dq_eff
        diag_resp_raw = np.diag(H) * dq_eff

        g0_eff = g0 * chain
        Hdq = Hdq_raw * chain
        diag_resp = diag_resp_raw * chain
        offdiag_resp = Hdq - diag_resp

        g_model = g0_eff + Hdq
        delta_g_err = g_qm - g_model
        rho = np.abs(delta_g_err) / (np.abs(Hdq) + eps)

        # Energy decomposition (kept as you had it; it mirrors your energy quadratic form)
        linear_energy = g0_eff * dq_eff
        diag_quad_energy = 0.5 * np.diag(H) * dq_eff * dq_eff

        cross_energy_split = np.zeros(N)
        for i in range(N):
            for j in range(i + 1, N):
                cross = dq_eff[i] * H[i, j] * dq_eff[j]
                cross_energy_split[i] += 0.5 * cross
                cross_energy_split[j] += 0.5 * cross

        pE_total = linear_energy + diag_quad_energy + cross_energy_split


        rows = []
        for i, coord in enumerate(z_matrix):
            rows.append({
                "idx": i,
                "coord": tuple(coord),
                "type": self.coord_type(tuple(coord)),
                "dq_raw": dq_raw[i],
                "dq_eff": dq_eff[i],
                "abs_dq": abs(dq_eff[i]),
                "chain": chain[i],
                "g0_eff": g0_eff[i],
                "diag_resp": diag_resp[i],
                "abs_diag_resp": abs(diag_resp[i]),
                "offdiag_resp": offdiag_resp[i],
                "abs_offdiag_resp": abs(offdiag_resp[i]),
                "Hdq": Hdq[i],
                "abs_Hdq": abs(Hdq[i]),
                "g_model": g_model[i],
                "g_qm": g_qm[i],
                "delta_g_err": delta_g_err[i],
                "abs_delta_g_err": abs(delta_g_err[i]),
                "rho": rho[i],
                "pE_total": pE_total[i],
                "abs_pE_total": abs(pE_total[i]),

            })

        by_abs_dq = sorted(rows, key=lambda r: r["abs_dq"], reverse=True)
        by_abs_Hdq = sorted(rows, key=lambda r: r["abs_Hdq"], reverse=True)
        by_abs_err = sorted(rows, key=lambda r: r["abs_delta_g_err"], reverse=True)
        by_rho = sorted(rows, key=lambda r: r["rho"], reverse=True)
        by_abs_pE = sorted(rows, key=lambda r: r["abs_pE_total"], reverse=True)

        return {
            "rows": rows,
            "dq_eff": dq_eff,
            "chain": chain,
            "Hdq": Hdq,
            "diag_resp": diag_resp,
            "offdiag_resp": offdiag_resp,
            "g_model": g_model,
            "delta_g_err": delta_g_err,
            "rho": rho,
            "pE_total": pE_total,


            "by_abs_dq": by_abs_dq,
            "by_abs_Hdq": by_abs_Hdq,
            "by_abs_err": by_abs_err,
            "by_rho": by_rho,
            "by_abs_pE": by_abs_pE,
        }


    def print_ranked_table(self,
        rows: List[Dict[str, Any]],
        title: str,
        n: int = 10,
    ) -> None:
        """Pretty-print a ranked summary table."""
        print(f"\n{title}")
        print("-" * len(title))
        header = (
            f"{'rank':>4} {'idx':>4} {'type':>9} {'coord':>18} "
            f"{'|dq|':>12} {'|diag|':>12} {'|offdiag|':>12} "
            f"{'|Hdq|':>12} {'|err|':>12} {'rho':>12}"
        )
        print(header)
        print("-" * len(header))

        for rank, r in enumerate(rows[:n]):
            print(
                f"{rank:4d} "
                f"{r['idx']:4d} "
                f"{r['type']:>9} "
                f"{self.format_coord(r['coord']):>18} "
                f"{r['abs_dq']:12.6e} "
                f"{r['abs_diag_resp']:12.6e} "
                f"{r['abs_offdiag_resp']:12.6e} "
                f"{r['abs_Hdq']:12.6e} "
                f"{r['abs_delta_g_err']:12.6e} "
                f"{r['rho']:12.6e}"
            )


    def print_source_decomposition(self,
        z_matrix: List[Tuple[int, ...]],
        dq: np.ndarray,
        H: np.ndarray,
        target_idx: int,
        top_m: int = 8,
    ) -> None:
        """
        For one response coordinate i, print the largest source terms H[i,j] * dq[j].
        This identifies which displaced coordinates drive the response in coordinate i.
        """
        z_matrix = [tuple(c) for c in z_matrix]
        dq = np.asarray(dq, dtype=float).reshape(-1)
        H = np.asarray(H, dtype=float)

        N = len(z_matrix)
        if target_idx < 0 or target_idx >= N:
            raise IndexError(f"target_idx {target_idx} out of range 0..{N-1}")

        terms = []
        for j in range(N):
            contrib = H[target_idx, j] * dq[j]
            terms.append({
                "j": j,
                "coord_j": z_matrix[j],
                "type_j": self.coord_type(z_matrix[j]),
                "dq_j": dq[j],
                "H_ij": H[target_idx, j],
                "contrib": contrib,
                "abs_contrib": abs(contrib),
            })

        terms_sorted = sorted(terms, key=lambda t: t["abs_contrib"], reverse=True)

        target_coord = z_matrix[target_idx]
        total = np.sum([t["contrib"] for t in terms])

        print(f"\nSource decomposition for target idx={target_idx}, coord={self.format_coord(target_coord)}")
        print(f"Total response sum_j H[i,j] * dq[j] = {total:.6e}")
        print(
            f"{'rank':>4} {'j':>4} {'type':>9} {'coord_j':>18} "
            f"{'dq_j':>12} {'H[i,j]':>12} {'H[i,j]*dq_j':>15}"
        )
        print("-" * 86)

        for rank, t in enumerate(terms_sorted[:top_m]):
            print(
                f"{rank:4d} "
                f"{t['j']:4d} "
                f"{t['type_j']:>9} "
                f"{self.format_coord(t['coord_j']):>18} "
                f"{t['dq_j']:12.6e} "
                f"{t['H_ij']:12.6e} "
                f"{t['contrib']:15.6e}"
            )


    def print_type_summary(self, rows: List[Dict[str, Any]]) -> None:
        """Aggregate diagnostics by coordinate type."""
        buckets: Dict[str, Dict[str, float]] = {}
        for r in rows:
            t = r["type"]
            if t not in buckets:
                buckets[t] = {
                    "count": 0,
                    "sum_abs_dq": 0.0,
                    "sum_abs_Hdq": 0.0,
                    "sum_abs_err": 0.0,
                    "sum_abs_pE": 0.0,
                }
            buckets[t]["count"] += 1
            buckets[t]["sum_abs_dq"] += r["abs_dq"]
            buckets[t]["sum_abs_Hdq"] += r["abs_Hdq"]
            buckets[t]["sum_abs_err"] += r["abs_delta_g_err"]
            buckets[t]["sum_abs_pE"] += r["abs_pE_total"]

        print("\nSummary by coordinate type")
        print("--------------------------")
        print(
            f"{'type':>9} {'count':>7} {'sum|dq|':>14} "
            f"{'sum|Hdq|':>14} {'sum|err|':>14} {'sum|pE|':>14}"
        )
        print("-" * 76)
        for t, vals in buckets.items():
            print(
                f"{t:>9} "
                f"{vals['count']:7d} "
                f"{vals['sum_abs_dq']:14.6e} "
                f"{vals['sum_abs_Hdq']:14.6e} "
                f"{vals['sum_abs_err']:14.6e} "
                f"{vals['sum_abs_pE']:14.6e}"
            )


    def run_full_internal_diagnostic(self,
        z_matrix: List[Tuple[int, ...]],
        dq: np.ndarray,
        H: np.ndarray,
        g0: np.ndarray,
        g_qm: np.ndarray,
        g_im: Optional[np.ndarray] = None,
        top_n: int = 10,
        decompose_top_err: int = 3,
    ) -> Dict[str, Any]:
        """
        Run the full diagnostic workflow and print useful summaries.

        Suggested usage:
            result = run_full_internal_diagnostic(z_matrix, dq, H, g0, g_qm)
        """
        result = self.compute_internal_gradient_diagnostics(
            z_matrix=z_matrix,
            dq_raw=dq,
            H=H,
            g0=g0,
            g_qm=g_qm,
        )

        rows = result["rows"]

        self.print_ranked_table(result["by_abs_dq"], "Ranked by |dq| (source displacement candidates)", n=top_n)
        self.print_ranked_table(result["by_abs_Hdq"], "Ranked by |H dq| (response channels)", n=top_n)
        self.print_ranked_table(result["by_abs_err"], "Ranked by |g_qm - g_model| (actual model failure)", n=top_n)
        self.print_ranked_table(result["by_rho"], "Ranked by rho = |err| / (|H dq| + eps)", n=top_n)
        self.print_ranked_table(result["by_abs_pE"], "Ranked by |partial energy contribution|", n=top_n)

        self.print_type_summary(rows)


        # Decompose the top few error coordinates to find their source coordinates
        for r in result["by_abs_err"][:decompose_top_err]:
            self.print_source_decomposition(
                z_matrix=z_matrix,
                dq=dq,
                H=H,
                target_idx=r["idx"],
                top_m=8,
            )
        
        return result   
    


    
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

        # Make zero the values of s_inv that are smaller than tol
        s_inv = np.array([1 / s_i if s_i > tol else 0.0 for s_i in s])
        
        # Critical assertion, check that the remaining positive values are equal to dimension (3N-6)
        number_of_positive_values = np.count_nonzero(s_inv)
        print(number_of_positive_values, dimension)
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

    def compute_internal_coordinates_values(self, coordinates):
        """
        Creates an array with the values of the internal coordinates
        and saves it in self.
        """

        cartesian_coordinates = coordinates

        n_atoms = cartesian_coordinates.shape[0]
        coords = cartesian_coordinates.reshape((n_atoms * 3))

        int_coords = []

        internal_coordinates = []
        for z in self.z_matrix:
            
            if len(z) == 2:
                q = geometric.internal.Distance(*z)
            elif len(z) == 3:
                q = geometric.internal.Angle(*z)
            elif len(z) == 4:
                q = geometric.internal.Dihedral(*z)
            else:
                assert_msg_critical(False, 'InterpolationDatapoint: Invalid entry size in Z-matrix.')
            internal_coordinates.append(q)


        for element, q in enumerate(internal_coordinates):
            if (isinstance(q, geometric.internal.Distance)):
                int_coords.append(1.0 / q.value(coords))
            elif len(self.z_matrix[element]) == 4:
                int_coords.append(np.sin(q.value(coords)))
            else:
                int_coords.append(q.value(coords))

        X = np.array(int_coords)

        return X
    
    def compute_ic_and_groups(self, coordinates):
        """
        Returns:
        X : np.ndarray shape (D,) feature vector
        groups : dict[str, np.ndarray] index arrays for feature groups
                keys: 'bond', 'angle', 'torsion' (torsion contains both sin/cos indices)
        types : list[str] parallel to features, e.g. 'bond','angle','torsion_sin','torsion_cos'
        """
        cart = np.asarray(coordinates, dtype=np.float64)
        n_atoms = cart.shape[0]
        coords = cart.reshape(n_atoms * 3)

        # Build geometric objects only once
        internal_coordinates = []
        for z in self.z_matrix:
            if len(z) == 2:
                internal_coordinates.append(('bond',     geometric.internal.Distance(*z)))
            elif len(z) == 3:
                internal_coordinates.append(('angle',    geometric.internal.Angle(*z)))
            elif len(z) == 4:
                internal_coordinates.append(('torsion',  geometric.internal.Dihedral(*z)))
            else:
                assert_msg_critical(False, 'Invalid entry size in Z-matrix.')

        X = []
        types = []
        bond_idx, angle_idx, torsion_idx = [], [], []
        k = 0

        for kind, q in internal_coordinates:
            val = q.value(coords)  # distances in Å, angles/torsions in rad (assuming)

            if kind == 'bond':
                # choose ONE: r or 1/r or log r ; avoid 1/r if you expect very short bonds
                feat = 1.0 / val                                # or: feat = val; or feat = np.log(val)
                X.append(feat); types.append('bond'); bond_idx.append(k); k += 1

            elif kind == 'angle':
                X.append(val); types.append('angle'); angle_idx.append(k); k += 1

            elif kind == 'torsion':
                # continuous encoding: add *two* features that should share a lengthscale
                X.append(np.sin(val)); types.append('torsion_sin'); torsion_idx.append(k); k += 1
                X.append(np.cos(val)); types.append('torsion_cos'); torsion_idx.append(k); k += 1

        X = np.asarray(X, dtype=np.float64)

        groups = {
            'bond':    np.asarray(bond_idx, dtype=int),
            'angle':   np.asarray(angle_idx, dtype=int),
            'torsion': np.asarray(torsion_idx, dtype=int),  # includes both sin & cos indices
        }
        return X, groups, types

    
