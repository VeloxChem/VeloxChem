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
import torch
import gpytorch
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
        self.use_cosine_dihedral = True

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
                    'use_cosine_dihedral':('bool', 'wether to use cosine and sin for the diehdral in the Z-matrix'),
                'labels': ('seq_fixed_str', 'the list of QM data point labels'),
            }
        }


        self.grimme_qf_params = {'k': {'H': 1.755, 'C': 2.463, 'N':2.559, 'O':2.579}, 'ken':-0.164,
                                 'EN':{'H': 2.1, 'C': 2.5, 'N':3.0, 'O':3.5}}
        self.eq_bond_force_constants = None

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

        # check_mol = Molecule(self.molecule.get_labels(), self.qm_data_points[0].cartesian_coordinates, 'bohr')
        # self.define_impes_coordinate(check_mol.get_coordinates_in_bohr())

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
        
        if not self.use_symmetry and 1==2:
            for i, data_point in enumerate(self.qm_data_points[:]):
                
                distance, denominator, weight_gradient, distance_vector, dihedral_dist = self.cartesian_distance(data_point)

                if abs(distance) < min_distance:
                    min_distance = abs(distance)

                distances_and_gradients.append((distance, dihedral_dist, i, denominator, weight_gradient, distance_vector))
        
        elif self.weightfunction_type == 'cartesian-hessian':
            for i, data_point in enumerate(self.qm_data_points[:]):
                
                distance, dihedral_dist, denominator, weight_gradient, distance_vec, _ = self.cartesian_hessian_distance(data_point)
              
                if abs(distance) < min_distance:
                    min_distance = abs(distance)

                distances_and_gradients.append((distance, dihedral_dist, i, denominator, weight_gradient, distance_vec))

        else:
            for i, data_point in enumerate(self.qm_data_points[:]):

                distance, dihedral_dist, denominator, weight_gradient, distance_vec, _ = self.cartesian_distance(data_point)
              

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
  
        
        #TODO: Grimme QF correction implementation comparison later
        self.grimme_qf_params = {'k': {'H': 1.755, 'C': 2.463, 'N':2.559, 'O':2.579}, 'ken':-0.164,
                                 'EN':{'H': 2.1, 'C': 2.5, 'N':3.0, 'O':3.5}}
            
        def grimme_function(r, re, kstr, a):
            grimme_func = kstr + kstr * (re / r)**a - 2.0 * kstr * (re / r)**(a/2)
            grimme_func_deriv = ((-a * kstr * (re/r)**a) / (r)) + ((a * kstr * (re / r)**(a/2)) / r)
            grimme_func_sec_deriv = ((a**2 * kstr * (re / r)**a) / (r**2) + (a * kstr * (re / r)**a) / (r**2) 
                                           - (a**2 * kstr * (re / r)**(a/2)) / (2.0 * r**2) - (a * kstr * (re / r)**(a/2.0)) / r**2)
            return grimme_func, grimme_func_deriv, grimme_func_sec_deriv
        
        def grimme_func_correction(r, re, kstr, a):
            V, dV, ddV = grimme_function(r, re, kstr, a)
            V0, dV0, ddV0 = grimme_function(re, re, kstr, a)
            
            dr = r - re
            dE = V - (V0 + dV0 * dr + 0.5 * ddV0 * dr**2)
            dg = dV - (dV0 + ddV0 * dr)
            dE_deriv_bmat = dg
            print('grimme corr', dr, dE, dE_deriv_bmat)
            return dE, dE_deriv_bmat   

        
        
        # --- initialise accumulators -------------------------------------------------
        self.impes_coordinate.energy    = 0.0
        self.impes_coordinate.gradient  = np.zeros((natms, 3))
        self.impes_coordinate.NAC       = np.zeros((natms, 3))       # if you need it

        # if self.eq_bond_force_constants is not None:
        #     for bond_idx, element_bond in enumerate(self.impes_coordinate.z_matrix[:]):
        #         if len(element_bond) !=2:
        #             break
                
        #         a = (self.grimme_qf_params['k'][self.molecule.get_labels()[element_bond[0]]] * self.grimme_qf_params['k'][self.molecule.get_labels()[element_bond[1]]]
        #             + self.grimme_qf_params['ken'] * (self.grimme_qf_params['EN'][self.molecule.get_labels()[element_bond[0]]] - self.grimme_qf_params['EN'][self.molecule.get_labels()[element_bond[1]]])**2)
        #         dE, dg = grimme_func_correction(self.impes_coordinate.internal_coordinates_values[bond_idx], self.eq_bond_force_constants[tuple(element_bond)]['r_eq'], self.eq_bond_force_constants[tuple(element_bond)]['k_st'], a)
        #         self.impes_coordinate.energy += dE
        #         self.impes_coordinate.gradient += (dg * self.impes_coordinate.b_matrix[bond_idx, :]).reshape(natms,3)

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

        # self.sum_of_weights      = W_i.sum()          # if you really need it later
        self.averaged_int_dist   = np.tensordot(W_i, averaged_int_dists, axes=1)
        
    def greedy_interpolation(self):

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


        ## determine the best optimal set for the interpolation of the current structure
        ## to generate the best Potential and Gradient approximation possible



        
        if not self.use_symmetry and 1==2:
            for i, data_point in enumerate(self.qm_data_points[:]):
                
                distance, denominator, weight_gradient, distance_vector, dihedral_dist = self.cartesian_distance(data_point)

                if abs(distance) < min_distance:
                    min_distance = abs(distance)

                distances_and_gradients.append((distance, dihedral_dist, i, denominator, weight_gradient, distance_vector))
        
        elif self.weightfunction_type == 'cartesian-hessian':
            for i, data_point in enumerate(self.qm_data_points[:]):
                
                distance, dihedral_dist, denominator, weight_gradient, distance_vec, _ = self.cartesian_hessian_distance(data_point)
              
                if abs(distance) < min_distance:
                    min_distance = abs(distance)

                distances_and_gradients.append((distance, dihedral_dist, i, denominator, weight_gradient, distance_vec))

        else:
            for i, data_point in enumerate(self.qm_data_points[:]):

                distance, dihedral_dist, denominator, weight_gradient, distance_vec, _ = self.cartesian_distance(data_point)
              

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

        # ∇U = Σ Wᵢ ∇Uᵢ  +  Σ Uᵢ ∇Wᵢ
        # if len(self.symmetry_information[3]) != natms:
        #     self.impes_coordinate.gradient = (np.tensordot(W_i, gradients, axes=1))
        #     # self.impes_coordinate.gradient[self.symmetry_information[4]] += (gradients[:, self.symmetry_information[4], :].sum(axis=0))
        #     # Add contributions only to the selected rows
            
        #     self.impes_coordinate.gradient[self.symmetry_information[3]] = np.tensordot(potentials, grad_W_i, axes=1)
        # else:

        self.impes_coordinate.gradient = (np.tensordot(W_i, gradients, axes=1) + np.tensordot(potentials, grad_W_i, axes=1))


        # --- 4.  book-keeping (optional) ---------------------------------------------
        for lbl, Wi in zip(used_labels, W_i):
            self.weights[lbl] = Wi

        # self.sum_of_weights      = W_i.sum()          # if you really need it later
        self.averaged_int_dist   = np.tensordot(W_i, averaged_int_dists, axes=1)



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
            
            if not self.use_cosine_dihedral:
                for i, element in enumerate(self.impes_coordinate.z_matrix[self.symmetry_information[-1][1]:], start=self.symmetry_information[-1][1]): 

                    dist_check[i] = np.sin(dist_org[i])
      
            self.bond_rmsd.append(np.sqrt(np.mean(np.sum((dist_org[:self.symmetry_information[-1][0]])**2))))
            self.angle_rmsd.append(np.sqrt(np.mean(np.sum(dist_org[self.symmetry_information[-1][0]:self.symmetry_information[-1][1]]**2))))
            self.dihedral_rmsd.append(np.sqrt(np.mean(np.sum(dist_check[self.symmetry_information[-1][1]:]**2))))

            pes = (energy + np.matmul(dist_check.T, grad) +
                        0.5 * np.linalg.multi_dot([dist_check.T, hessian, dist_check]))

            dist_hessian_eff = np.matmul(dist_check.T, hessian)


            if not self.use_cosine_dihedral:
                for i, element in enumerate(self.impes_coordinate.z_matrix[self.symmetry_information[-1][1]:], start=self.symmetry_information[-1][1]):
                    
                    grad[i] *= np.cos(dist_org[i])
                    dist_hessian_eff[i] *= np.cos(dist_org[i])
            
            
            pes_prime = (np.matmul(self.impes_coordinate.b_matrix.T, (grad + dist_hessian_eff))).reshape(natm, 3)

            return pes, pes_prime, (grad + dist_hessian_eff)

    
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
    
    def trust_radius_weight_gradient_hessian(self, confidence_radius, distance):

        denominator = (
                (distance / confidence_radius)**(2 * self.exponent_p) +
                (distance / confidence_radius)**(2 * self.exponent_q))
        trust_radius_weight_gradient = -1.0 * ((( -2.0 * self.exponent_p * ((distance / confidence_radius)**(2 * self.exponent_p)) / confidence_radius) - 
                                             (2.0 * self.exponent_q * ((distance / confidence_radius)**(2 * self.exponent_q) / confidence_radius))) / denominator**2)
        return trust_radius_weight_gradient

    
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

    
    def trust_radius_weight_gradient_gradient_hessian(self, confidence_radius, distance, grad_s):
        
   
        denominator = (
                (distance / confidence_radius)**(2 * self.exponent_p) +
                (distance / confidence_radius)**(2 * self.exponent_q))

        trust_radius_weight_gradient_gradient_nominator_1_1 = 2.0 * ((2.0 * self.exponent_p * (distance / confidence_radius)**(2 * self.exponent_p - 1) / confidence_radius) + (2.0 * self.exponent_q * (distance / confidence_radius)**(2 * self.exponent_q - 1) / confidence_radius))
        trust_radius_weight_gradient_gradient_nominator_1_2 = ((-2.0 * self.exponent_p * (distance / confidence_radius)**(2 * self.exponent_p) / confidence_radius) - ( 2.0 * self.exponent_q * (distance / confidence_radius)**(2 * self.exponent_q) / confidence_radius))

        trust_radius_weight_gradient_gradient_nominator_2 = ((-2.0 * self.exponent_p * (2 * self.exponent_p - 1) * (distance / confidence_radius)**(2 * self.exponent_p - 1) / confidence_radius**2) - (2.0 * self.exponent_p * (distance / confidence_radius)**(2 * self.exponent_p - 1) / confidence_radius**2) -
                                                                    (2.0 * self.exponent_q * (2 * self.exponent_q - 1) * (distance / confidence_radius)**(2 * self.exponent_q - 1) / confidence_radius**2) - (2.0 * self.exponent_q * (distance / confidence_radius)**(2 * self.exponent_q - 1) / confidence_radius**2))
        
        trust_radius_weight_gradient_gradient = (((trust_radius_weight_gradient_gradient_nominator_1_1 * trust_radius_weight_gradient_gradient_nominator_1_2 ) / denominator**3) - trust_radius_weight_gradient_gradient_nominator_2 / denominator**2) * grad_s


        return trust_radius_weight_gradient_gradient
    
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
    
    def shepard_weight_gradient_hessian(self, quad_distance, confidence_radius, grad_s):
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
    
    
    def calculate_translation_coordinates(self, given_coordinates):
        """Center the molecule by translating its geometric center to (0, 0, 0)."""
        center = np.mean(given_coordinates, axis=0)
        translated_coordinates = given_coordinates - center

        return translated_coordinates

    def cartesian_distance(self, data_point):
        """Calculates and returns the cartesian distance between
           self.coordinates and data_point coordinates.
           Besides the distance, it also returns the weight gradient,
           which requires the distance vector to be computed.

           :param data_point:
                InterpolationDatapoint object
        """

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
        R  = geometric.rotate.get_rot(Yc, Xc)

        rotated_current = (R @ Yc.T).T
        distance_vector_sub =  Xc - rotated_current
        distance = np.linalg.norm(distance_vector_sub)

        distance_vector = np.zeros_like(reference_coordinates)   # (natms,3)
        distance_vector[self.symmetry_information[3]] = distance_vector_sub
  


        dihedral_dist = 0.0
        if distance < 1e-8:
            distance = 1e-8
            distance_vector_sub[:] = 0
            

        if self.interpolation_type == 'shepard':
            denominator, weight_gradient_sub = self.shepard_weight_gradient(
                distance_vector_sub, distance, data_point.confidence_radius)
            if self.calc_optim_trust_radius:
                dw_dalhpa_i = self.trust_radius_weight_gradient(data_point.confidence_radius, distance)
                dw_dX_dalpha_i = self.trust_radius_weight_gradient_gradient(data_point.confidence_radius, distance, distance_vector)
                self.dw_dalpha_list.append(dw_dalhpa_i)
                self.dw_dX_dalpha_list.append(dw_dX_dalpha_i)        
        elif self.interpolation_type == 'simple':
            weight_gradient = self.simple_weight_gradient(
                distance_vector_sub, distance)
        else:
            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)
        
        weight_gradient = np.zeros_like(reference_coordinates)   # (natms,3)
        weight_gradient[self.symmetry_information[3]] = weight_gradient_sub

        return distance, dihedral_dist, denominator, weight_gradient, distance_vector, grad_s

    def cartesian_hessian_distance(self, data_point):
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
        # print('carteisna distance', distance)
        H_full = data_point.hessian                          # (3N,3N)

        idx3 = np.concatenate([3*active_atoms[:,None] + np.array([0,1,2])], axis=1).ravel()
        H_SS = H_full[np.ix_(idx3, idx3)]

        Nc = Xc.shape[0]
        J = np.eye(Nc) - np.ones((Nc, Nc))/Nc      # (Nc,Nc)
        J3 = np.kron(J, np.eye(3))                 # (3Nc,3Nc)

        Hc = J3.T @ H_SS @ J3 
        E_ref = 1e-3 / 0.001593 # in hartree to scale the distance making it unitless
        
        Hc_ul = Hc * E_ref
        H_spd = make_spd_hessian(Hc_ul, floor_abs=1e-6)

        w, U   = np.linalg.eigh(H_spd)
        lam_med = np.median(w)
        lam_lo  = 0.05 * lam_med
        lam_hi  = 10.0 * lam_med
        w_clip  = np.clip(w, lam_lo, lam_hi)
        H_clip  = (U * w_clip) @ U.T
        # H_clip = H_clip
 
        distance_2     = (distance_vector.reshape(-1) @ ((1.0 * H_clip) @ distance_vector.reshape(-1)))
        distance = np.sqrt(distance_2)
        if distance < 1e-8:
            distance = 1e-8
            distance_vector[:] = 0
        y = ((1.0 * H_clip) @ distance_vector.reshape(-1)).reshape(-1, 3)      # (Nc,3)

        rigid = y @ R_x                           # (Nc,3)
        S = y.T @ Xc                                    # (3,3)

        resp = np.einsum('baij,ij->ba', dU_x, S, optimize=True)   # (Nc,3)

        z_mat = rigid + resp                            # (Nc,3)
        z_mat = apply_Jc(z_mat)                        # (Nc,3)

        grad_s_sub = 1.0 * z_mat / distance                            # (Nc,3) 

        grad_s = np.zeros_like(reference_coordinates)   # (natms,3)
        grad_s[self.symmetry_information[3]] = grad_s_sub

        # components in the eigenbasis
        c = U.T @ distance_vector.reshape(-1)      # projections (u_k · d)
        contrib = w_clip * (c**2)                   # λ_k * (u_k·d)^2
        # print("Top contributing modes (λ_k * c_k^2):")
        # idx = np.argsort(contrib)[::-1]
        # for j in idx[:5]:
        #     print(f"k={j:2d}  λ={w_clip[j]:.3e}  (u·d)^2={c[j]**2:.3e}  term={contrib[j]:.3e}")


        dihedral_dist = 0.0

        if self.interpolation_type == 'shepard':
            denominator, weight_gradient = self.shepard_weight_gradient_hessian(distance, data_point.confidence_radius, grad_s)  
            if self.calc_optim_trust_radius:
                dw_dalhpa_i = self.trust_radius_weight_gradient_hessian(data_point.confidence_radius, distance)
                dw_dX_dalpha_i = self.trust_radius_weight_gradient_gradient_hessian(data_point.confidence_radius, distance, grad_s)
                self.dw_dalpha_list.append(dw_dalhpa_i)
                self.dw_dX_dalpha_list.append(dw_dX_dalpha_i)
            
        elif self.interpolation_type == 'simple':
            weight_gradient = self.simple_weight_gradient(distance_vector, distance)
        else:
            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)
        

        

        return distance, dihedral_dist, denominator, weight_gradient, distance_vector, grad_s


    def determine_important_internal_coordinates(self, qm_energy, qm_gradient, molecule, z_matrix, datapoints, dihedral_diff_const=True):
        """Determines the most important internal coordinates
           leading to the large deviation.
        """

        selection_rule='relative' # 'relative' | 'coverage' | 'topk'
        relative_threshold=0.7
        coverage_mass=0.8
        topk=None
        

        masses = molecule.get_masses().copy()
        masses_cart = np.repeat(masses, 3)
        sqrt_masses = 1.0 / np.sqrt(masses_cart)
        
        qm_gradient_mw = qm_gradient.reshape(-1) * sqrt_masses    # mass-weighted gradient
        constraints = []
        print(len(datapoints), self.z_matrix, self.symmetry_information)
        for datapoint in datapoints:
            dihedral_difference = []
            internal_coord_elem_distance = []
            org_interal_coord_elem_distance = []
            for elem_idx, element in enumerate(self.z_matrix): 
                if len(element) == 4 and tuple(sorted(element)) in self.symmetry_information[7][3]: 
                    # internal_coord_elem_distance.append(0.5 * ( 1.0 + np.cos(3.0 * (self.impes_coordinate.internal_coordinates_values[elem_idx] - datapoint.internal_coordinates_values[elem_idx]) + np.pi)))
                    internal_coord_elem_distance.append(0.0)
                    org_interal_coord_elem_distance.append(0.0)
                
                # elif len(element) == 4 and tuple(sorted(element)) in self.symmetry_information[7][2]: 
                #     # internal_coord_elem_distance.append(0.5 * ( 1.0 + np.cos(2.0 * (self.impes_coordinate.internal_coordinates_values[elem_idx] - datapoint.internal_coordinates_values[elem_idx]) + np.pi)))                
                #     internal_coord_elem_distance.append(0.0)
                #     org_interal_coord_elem_distance.append(0.0)
                elif len(element) == 4:
                    
                    internal_coord_elem_distance.append(np.sin(self.impes_coordinate.internal_coordinates_values[elem_idx] - datapoint.internal_coordinates_values[elem_idx]))
                    org_interal_coord_elem_distance.append(np.cos(self.impes_coordinate.internal_coordinates_values[elem_idx] - datapoint.internal_coordinates_values[elem_idx]))
                    print(z_matrix[elem_idx], self.impes_coordinate.internal_coordinates_values[elem_idx], datapoint.internal_coordinates_values[elem_idx])
                    if abs(internal_coord_elem_distance[-1]) > 0.3 and dihedral_diff_const:
                        dihedral_difference.append(element)
                    
                else:
                    org_interal_coord_elem_distance.append(1.0)
                    internal_coord_elem_distance.append(self.impes_coordinate.internal_coordinates_values[elem_idx] - datapoint.internal_coordinates_values[elem_idx])
            
            check_molecule = Molecule(molecule.get_labels(), datapoint.cartesian_coordinates, 'bohr')
            print(check_molecule.get_xyz_string())
            print(dihedral_difference)
            

            N = len(z_matrix)

            partial_energies = np.zeros(N)
            partial_gradient = np.zeros(N) 

            
            
            # First handle linear parts
            for i in range(N):
                partial_energies[i] += internal_coord_elem_distance[i] * datapoint.internal_gradient[i]
                partial_gradient[i] += org_interal_coord_elem_distance[i] * datapoint.internal_gradient[i]

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

            pred_E, pred_G_mw, _ = self.compute_potential(datapoint, self.impes_coordinate.internal_coordinates_values)
            pred_im_G_int = self.transform_gradient_to_internal_coordinates(molecule, pred_G_mw, self.impes_coordinate.b_matrix)
            pred_qm_G_int = self.transform_gradient_to_internal_coordinates(molecule, qm_gradient_mw.reshape(qm_gradient.shape), self.impes_coordinate.b_matrix)

            print("max energy diff",
                abs((E_pred_check - pred_E)))
            

            
            
            g_pred_check = partial_gradient           # already Hartree/bohr
            max_grad_diff = np.max(np.abs(g_pred_check - pred_im_G_int))
            print("max grad diff", max_grad_diff)
            print(g_pred_check, pred_im_G_int)

         

            delta_E = abs(qm_energy - pred_E)
            delta_G = np.linalg.norm(pred_qm_G_int - pred_im_G_int)

            print(' \n\n Energy error with QM ', delta_E * hartree_in_kcalpermol(), delta_G)
   
            delta_g = pred_qm_G_int - pred_im_G_int
            print('delta g', delta_g)
            single_energy_error = []
            weights = []
            single_gradient_error = []
            abs_pg       = np.abs(partial_gradient)
            sum_abs_pg   = abs_pg.sum()

            # effective displacements already computed as internal_coord_elem_distance (Δq_eff)
            eps    = 1e-12
            dq_eff = np.array(internal_coord_elem_distance)

            # --- energy part, includes gradient mismatch via abs_work ---
            abs_pE   = np.abs(partial_energies)
            abs_work = np.abs(delta_g) * np.abs(dq_eff)          # |Δg_i| |Δq_i|
            w_E_raw  = abs_pE * abs_work
            w_E      = w_E_raw / (w_E_raw.sum() + eps)
            e_i      = delta_E * w_E                             # Σ e_i = ΔE

            # --- gradient part (your existing decomposition) ---
            abs_pG     = np.abs(partial_gradient)
            sum_abs_pG = abs_pG.sum()
            if sum_abs_pG < eps:
                w_G = np.zeros_like(abs_pG)
                g_i = np.zeros_like(abs_pG)
            else:
                w_G = abs_pG / (sum_abs_pG + eps)
                g_i = delta_G * w_G                              # Σ g_i = ΔG

            # --- put both on an energy-like scale ---
            E_kcal = e_i * hartree_in_kcalpermol()

            L_ref         = 0.1                                  # bohr, choose sensibly
            G_as_energy   = g_i * L_ref                          # Hartree
            G_kcal        = G_as_energy * hartree_in_kcalpermol() * bohr_in_angstrom()
            print('Gradient in kcal', G_kcal)
            lambda_grad   = 0.5                                  # 0…1, how much you care about gradient
            score_i_kcal  = (1.0 - lambda_grad) * E_kcal + lambda_grad * G_kcal
            score_i_kcal  = np.maximum(score_i_kcal, 0.0)

            w_tot = score_i_kcal / (score_i_kcal.sum() + eps)

            # for bookkeeping: check sums
            assert np.allclose(e_i.sum(), delta_E, atol=1e-12)
            assert np.allclose(g_i.sum(), delta_G, atol=1e-12)

            # --- final per-coordinate entry ---
            contributions = list(
                zip(
                    partial_energies,
                    E_kcal,              # energy-error share (kcal/mol)
                    G_kcal,              # gradient-error share in energy units (kcal/mol)
                    w_tot,               # combined weight
                    z_matrix,
                )
            )

            sorted_contributions = sorted(
                contributions,
                key=lambda x: x[1] + x[2],        # total error = energy + gradient part
                reverse=True,
            )
            sorted_weights = np.array([c[3] for c in sorted_contributions])
            sorted_coords  = [tuple(int(x) for x in c[4]) for c in sorted_contributions]

            selected_idx = []
            if selection_rule == 'relative':
                if sorted_weights.size > 0:
                    wmax = sorted_weights.max()
                    keep = sorted_weights >= (relative_threshold * wmax)
                    selected_idx = list(np.where(keep)[0])
            elif selection_rule == 'coverage':
                cum = np.cumsum(sorted_weights)
                keep = cum <= coverage_mass
                if keep.size > 0:
                    keep[0] = True
                selected_idx = list(np.where(keep)[0])
            elif selection_rule == 'topk':
                k = topk if (topk is not None) else max(1, int(0.25 * N))
                selected_idx = list(range(min(k, N)))
            else:
                raise ValueError(f"Unknown selection_rule: {selection_rule}")

            

            for idx in selected_idx:
                coord = sorted_coords[idx]
                if tuple(coord) in constraints:
                    continue
                # optional symmetry filtering for torsions
                if len(coord) == 4 and (tuple(sorted(coord)) in self.symmetry_information[7][3]):
                    continue

                constraints.append(tuple(coord))
            
            for dev_dihedral in dihedral_difference:
                if tuple(dev_dihedral) in constraints:
                    continue
                
                constraints.append(tuple(dev_dihedral))

            print('Top contributors (first 10):')
            for i, (pE, e_kcal, g_kcal, w_i, coord) in enumerate(sorted_contributions[:10]):
                picked = 'SELECTED' if i in selected_idx else ''
                print(f"{i:2d} {tuple(int(x) for x in coord)}:"
                    f" |pE|={abs(pE):.3e} Ha,"
                    f" E_share={e_kcal:.3f} kcal/mol,"
                    f" G_share={g_kcal:.3f} kcal/mol,"
                    f" w_tot={w_i:.3f} {picked}")
            print('Sum of energy-weights:', float(sorted_weights.sum()))
            print('Selected constraints so far:', constraints)

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

    