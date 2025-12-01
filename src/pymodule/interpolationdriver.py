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



    def shepard_interpolation(self):
        """Performs a simple interpolation.

        :param qm_data_points:
            A list of InterpolationDatapoint corresponding to configurations
            to be used for interpolation (which have been calculated
            quantum mechanically).
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

        distances_and_gradients = []
        min_distance = float('inf')
        self.time_step_reducer = False
        
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

    