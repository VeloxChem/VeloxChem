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
import numpy as np
from time import time
from typing import List, Tuple, Dict, Any, Optional
import sys

from .profiler import Profiler
import h5py

from contextlib import redirect_stderr
from io import StringIO
from .interpolationdatapoint import InterpolationDatapoint

from .outputstream import OutputStream
from .veloxchemlib import mpi_master, bohr_in_angstrom

from .errorhandler import assert_msg_critical
from .inputparser import (parse_input)



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

        # simple or Shepard interpolation
        self.interpolation_type = 'shepard'
        self.weightfunction_type = 'cartesian'
        self.exponent_p = None
        self.scaling_time = False
        self.exponent_q = None
        self.confidence_radius = 0.5
        self.alpha = 1.0
        
        self.z_matrix = z_matrix
        self.impes_dict = None
        self.sum_of_weights = None
        self.sum_of_weights_grad = None
        self.print = False
        self.store_weights = False
        # symmetry metadata tuple produced during database-point generation
        self.symmetry_information = None
        # global switch controlling symmetry-expanded candidate evaluation

        self.weights = {}
        self.int_coord_weights = {}
        self.potentials = []
        self.gradients = []

        # kpoint file with QM data
        self.imforcefield_file = None
        self.qm_data_points = None
    
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
        self.use_cosine_dihedral = False
        self.use_tc_weights = False

        # Optional runtime profiling for interpolation bottleneck analysis.
        self.runtime_profile_enabled = False
        self.runtime_profile_print = False
        self.runtime_profile_profiler = None
        self.runtime_profile_totals = {}
        self.runtime_profile_last = {}
        self.runtime_profile_calls = 0

        self.tc_imp_gate_lambda = 0.7

        # Coordinate tolerances defining the sphere of influence.
        self.tc_imp_bond_sigma_angstrom = 0.04
        self.tc_imp_angle_sigma_degrees = 5.0
        self.tc_imp_dihedral_sigma_degrees = 12.0
        self.tc_imp_improper_sigma_degrees = 5.0

        # Numerical guards.
        self.tc_imp_gate_exp_clip = 50.0
        self.tc_imp_min_sigma = 1.0e-12

        # Optional runtime cache for symmetry-expanded task evaluation.
        # enables reuse of flattened symmetry tasks across compute() calls
        self.runtime_data_cache_enabled = True
        # invalidation flag for rebuilt symmetry task cache payload
        self._runtime_data_cache_dirty = True
        # core-label keyed cache of (datapoint, mask0, mask) symmetry tasks
        self._symmetry_task_cache = {}



        self._input_keywords = {
            'im_settings': {
                'interpolation_type':
                    ('str', 'type of interpolation (simple/Shepard)'),
                'weightfunction_type':
                    ('str', 'type of interpolation (cartesian/internal/cartesian-hessian)'),
                'exponent_p': ('int', 'the main exponent'),
                'exponent_q': ('int', 'the additional exponent (Shepard IM)'),
                'confidence_radius': ('float', 'the confidence radius'),
                'imforcefield_file':
                    ('str', 'the name of the chk file with QM data'),
                    'use_inverse_bond_length': ('bool', 'whether to use inverse bond lengths in the Z-matrix'),
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

                
        if self.interpolation_type == 'simple':
            AssertionError("simple interpolation scheme is not supported as it is considered worse then shepard interpolation.")

        if self.interpolation_type == 'shepard' and self.exponent_q is None:
            self.exponent_q = self.exponent_p / 2.0

    def enable_runtime_profiling(self, enabled=True, reset=True, print_summary=False):
        """
        Enables/disables lightweight runtime profiling for compute() calls.

        :param enabled:
            If True, timing data are collected.
        :param reset:
            If True, clear previous timing accumulators.
        :param print_summary:
            If True, print a short summary from compute().
        """

        self.runtime_profile_enabled = bool(enabled)
        self.runtime_profile_print = bool(print_summary)

        if not self.runtime_profile_enabled:
            self.runtime_profile_profiler = None
            return

        if self.runtime_profile_profiler is None or reset:
            self.runtime_profile_profiler = Profiler({'timing': True})
            self.runtime_profile_profiler.set_timing_key('InterpolationDriver.compute')
            self.runtime_profile_totals = {}
            self.runtime_profile_last = {}
            self.runtime_profile_calls = 0

    def _add_runtime_timing(self, label, dt):
        """
        Adds timing information to interpolation runtime profile.
        """

        if not self.runtime_profile_enabled or self.runtime_profile_profiler is None:
            return

        self.runtime_profile_profiler.add_timing_info(label, dt)
        self.runtime_profile_totals[label] = self.runtime_profile_totals.get(label, 0.0) + dt
        self.runtime_profile_last[label] = self.runtime_profile_last.get(label, 0.0) + dt

    def get_runtime_profile_summary(self):
        """
        Returns interpolation runtime profiling data.
        """

        return {
            'enabled': self.runtime_profile_enabled,
            'calls': self.runtime_profile_calls,
            'totals': dict(self.runtime_profile_totals),
            'last': dict(self.runtime_profile_last),
        }


    def mark_runtime_data_cache_dirty(self):
        """
        Marks the flattened symmetry task cache as outdated.
        """

        self._runtime_data_cache_dirty = True
    
    def _get_datapoint_by_label(self, dp_label):
        for dp in (self.qm_data_points or []):
            if getattr(dp, 'point_label', None) == dp_label:
                return dp
        raise KeyError(
            f"InterpolationDriver: datapoint label '{dp_label}' was not found in qm_data_points.")

    def _get_internal_coordinate_partitions(self):
        """
        Returns canonical internal-coordinate index partitions.

        Layout:
            bonds:            [0, bond_end)
            angles:           [bond_end, angle_end)
            proper dihedrals: [dihedral_start, dihedral_end)
            impropers/tail:   [dihedral_end, n_ic)
        """
        zmat_dict = getattr(self.impes_coordinate, 'z_matrix_dict', None)
        if isinstance(zmat_dict, dict):
            n_bonds = len(zmat_dict.get('bonds', []))
            n_angles = len(zmat_dict.get('angles', []))
            n_dihedrals = len(zmat_dict.get('dihedrals', []))
            dihedral_start = n_bonds + n_angles
            dihedral_end = dihedral_start + n_dihedrals
            return {
                'bond_end': int(n_bonds),
                'angle_end': int(dihedral_start),
                'dihedral_start': int(dihedral_start),
                'dihedral_end': int(dihedral_end),
            }

        if self.symmetry_information is not None:
            dihedral_start = int(self.symmetry_information[-1][0])
            dihedral_end = int(self.symmetry_information[-1][1])
            return {
                'bond_end': int(dihedral_start),
                'angle_end': int(dihedral_start),
                'dihedral_start': int(dihedral_start),
                'dihedral_end': int(dihedral_end),
            }

        n_ic = 0
        if getattr(self.impes_coordinate, 'z_matrix', None) is not None:
            n_ic = len(self.impes_coordinate.z_matrix)
        return {
            'bond_end': 0,
            'angle_end': 0,
            'dihedral_start': 0,
            'dihedral_end': int(n_ic),
        }
    
    def _build_candidates_for_label(self, dp_label):
        dp = self._get_datapoint_by_label(dp_label)   # from qm_data_points
        n_ic = len(dp.internal_coordinates_values)
        ident = np.arange(n_ic, dtype=np.int64)
        return [(dp, ident, ident)]

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
        z_matrix_impropers = z_matrix_label + '_impropers'

        z_matrix_labels = [z_matrix_bonds, z_matrix_angles, z_matrix_dihedrals, z_matrix_impropers]
        z_matrix = {"bonds": [], "angles": [], "dihedrals": [], "impropers": []}
            
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
                        if 'bonds' in z_label:
                            z_matrix['bonds'] = current_z_list
                        elif 'angles' in z_label:
                            z_matrix['angles'] = current_z_list
                        elif 'dihedrals' in z_label:
                            z_matrix['dihedrals'] = current_z_list
                        elif 'impropers' in z_label:
                            z_matrix['impropers'] = current_z_list
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
        if self.runtime_profile_enabled:
            timing_info = {}
            self.impes_coordinate.reset_coordinates_impes_driver(
                coordinates,
                timing_info=timing_info)
            for key, dt in timing_info.items():
                self._add_runtime_timing(
                    f'compute.define_impes_coordinate.{key}', dt)
        else:
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
        compute_t0 = time()
        if self.runtime_profile_enabled:
            self.runtime_profile_last = {}

        prep_t0 = time()
        self.bond_rmsd = []
        self.angle_rmsd = []
        self.dihedral_rmsd = []
        
        self.molecule = molecule 

        self._add_runtime_timing('compute.prepare_state', time() - prep_t0)

        define_t0 = time()
        self.define_impes_coordinate(molecule.get_coordinates_in_bohr())
        self._add_runtime_timing('compute.define_impes_coordinate', time() - define_t0)

        load_t0 = time()
        if self.qm_data_points is None:
            self.qm_data_points = self.read_qm_data_points()
        self._add_runtime_timing('compute.load_qm_data_points', time() - load_t0)

        if self.runtime_data_cache_enabled and self._runtime_data_cache_dirty:
            cache_t0 = time()
            self._add_runtime_timing('compute.prepare_runtime_cache', time() - cache_t0)

        if self.interpolation_type == 'shepard':
            interp_t0 = time()
            self.shepard_interpolation()
            self._add_runtime_timing('compute.shepard_interpolation', time() - interp_t0)
        else:
            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)

        self._add_runtime_timing('compute.total', time() - compute_t0)
        if self.runtime_profile_enabled:
            self.runtime_profile_calls += 1

        if self.runtime_profile_enabled and self.runtime_profile_print:
            total = self.runtime_profile_last.get('compute.total', 0.0)
            shep = self.runtime_profile_last.get('compute.shepard_interpolation', 0.0)
            if total > 0.0:
                print(
                    f"[InterpolationDriver] compute.total={total:.6f}s "
                    f"(shepard={shep:.6f}s, {100.0 * shep / total:.1f}% )"
                )

    def shepard_interpolation(self):
        """Performs a simple interpolation.

        :param qm_data_points:
            A list of InterpolationDatapoint corresponding to configurations
            to be used for interpolation (which have been calculated
            quantum mechanically).
        """
        shepard_t0 = time()
    
        natms = self.impes_coordinate.cartesian_coordinates.shape[0]

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

        distances_and_gradients = []
        min_distance = float('inf')
        self.time_step_reducer = False
        
        # self.calc_optim_trust_radius = True
        
        distance_scan_t0 = time()

        if self.weightfunction_type == 'cartesian':
            for data_point in self.qm_data_points[:]:

                distance, dihedral_dist, denominator, weight_gradient, distance_vec, _, _, _ = self.cartesian_distance(data_point)
                
                
                if abs(distance) < min_distance:
                    min_distance = abs(distance)
                distances_and_gradients.append((distance, dihedral_dist, self.qm_data_points.index(data_point), denominator, weight_gradient, distance_vec))

        elif self.weightfunction_type == 'internal':
            for i, data_point in enumerate(self.qm_data_points[:]):
                distance, dihedral_dist, denominator, weight_gradient, distance_vec, _, _, _ = self.internal_distance(data_point)

                if abs(distance) < min_distance:
                    min_distance = abs(distance)
                distances_and_gradients.append((distance, dihedral_dist, self.qm_data_points.index(data_point), denominator, weight_gradient, distance_vec))

        else:
            errtxt = "Unrecognized weight function type: "
            errtxt += self.weightfunction_type
            raise ValueError(errtxt)

        self._add_runtime_timing('shepard.distance_scan', time() - distance_scan_t0)
        
        filter_t0 = time()
        close_distances = None
        close_distances = [
            (self.qm_data_points[index], distance, dihedral_dist, denom, wg, distance_vec, index) 
            for distance, dihedral_dist, index, denom, wg, distance_vec in distances_and_gradients]
        self._add_runtime_timing('shepard.filter_close_points', time() - filter_t0)

        eval_loop_t0 = time()
        compute_potential_acc = 0.0

        for seq, (qm_data_point, distance, dihedral_dist, denominator_cart, weight_grad_cart, _, label_idx) in enumerate(close_distances):

            weight_cart = 1.0 / (denominator_cart)

            # if outer_results is None:
            potential_t0 = time()
            potential, gradient_mw, r_i = self.compute_potential(
                qm_data_point,
                self.impes_coordinate.internal_coordinates_values,
            )
            compute_potential_acc += time() - potential_t0
        
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

        self._add_runtime_timing('shepard.point_eval_loop', time() - eval_loop_t0)
        self._add_runtime_timing('shepard.compute_potential', compute_potential_acc)
           
        assembly_t0 = time()
        # --- initialise accumulators -------------------------------------------------
        self.impes_coordinate.energy    = 0.0
        self.impes_coordinate.gradient  = np.zeros((natms, 3))
        self.impes_coordinate.NAC       = np.zeros((natms, 3))       # if you need it

        # --- 1.  raw (unnormalised) weights and their gradients ----------------------
        w_i          = np.array(weights_cart, dtype=np.float64)        # ← rename
        grad_w_i     = np.array(weight_gradients_cart, dtype=np.float64)   # shape (n_pts, natms, 3)
   
        S            = w_i.sum()                        # Σ wᵢ
        sum_grad_w   = grad_w_i.sum(axis=0)            # Σ ∇wᵢ      shape (natms, 3)

        eps = 1e-12

        if S > eps:
            W_i = w_i / S
            grad_W_i = (grad_w_i * S - w_i[:, None, None] * sum_grad_w) / S**2
        else:
            # fallback if all raw weights are numerically zero
            n_pts = len(w_i)
            W_i = np.full_like(w_i, 1.0 / n_pts)
            grad_W_i = np.zeros_like(grad_w_i)

        self.sum_of_weights = S
        self.sum_of_weights_grad = sum_grad_w
        
        # --- 3.  accumulate energy and gradient --------------------------------------
        potentials   = np.array(potentials, dtype=np.float64)        # Uᵢ
        gradients    = np.array(gradients,  dtype=np.float64)        # ∇Uᵢ  shape (n_pts, natms, 3)

        self.impes_coordinate.energy   = np.dot(W_i, potentials)     # Σ Wᵢ Uᵢ

        self.impes_coordinate.gradient = (np.tensordot(W_i, gradients, axes=1) + np.tensordot(potentials, grad_W_i, axes=1))
        
        # --- 4.  book-keeping (optional) ---------------------------------------------
        for lbl, Wi in zip(used_labels, W_i):
            self.weights[lbl] = Wi
        # print(self.weights)
        self.averaged_int_dist   = np.tensordot(W_i, averaged_int_dists, axes=1)
        self._add_runtime_timing('shepard.assembly', time() - assembly_t0)
        self._add_runtime_timing('shepard.total', time() - shepard_t0)


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
    
    def compute_potential(self, data_point, org_int_coords):
        """Calculates the potential energy surface at self.impes_coordinate
           based on the energy, gradient and Hessian of data_point.

           :param data_point:
                InterpolationDatapoint object.
        """
        def principal_angle(delta):
            return (delta + np.pi) % (2.0 * np.pi) - np.pi


        pes = 0.0
        natm = data_point.cartesian_coordinates.shape[0]
        # print(len(self.qm_symmetry_data_points))
        
        bounds = self._get_internal_coordinate_partitions()
        bond_end = int(bounds['bond_end'])
        angle_end = int(bounds['angle_end'])
        dihedral_start = int(bounds['dihedral_start'])
        dihedral_end = int(bounds['dihedral_end'])

        energy = data_point.energy
        grad = data_point.internal_gradient.copy()
        hessian = data_point.internal_hessian.copy()
        dist_org = (org_int_coords.copy() - data_point.internal_coordinates_values)
        dist_check = (org_int_coords.copy() - data_point.internal_coordinates_values)
        chain = np.ones_like(dist_org)
        if not self.use_cosine_dihedral:
            d_prop = self._principal_torsion_delta(dist_org[dihedral_start:dihedral_end])
            dist_check[dihedral_start:dihedral_end] = np.sin(d_prop)
            d_imp = self._principal_torsion_delta(
                dist_org[dihedral_end:]
            )
            chain[dihedral_start:dihedral_end] = np.cos(d_prop)
            imp_slice = slice(dihedral_end, len(dist_org))
            dist_check[imp_slice] = 2.0 * np.tan(0.5 * d_imp)
            chain[imp_slice] = 1.0 / np.maximum(np.cos(0.5 * d_imp)**2, 1.0e-12)

        self.bond_rmsd.append(np.sqrt(np.mean(np.sum((dist_org[:bond_end])**2))))
        self.angle_rmsd.append(np.sqrt(np.mean(np.sum(dist_org[bond_end:angle_end]**2))))
        self.dihedral_rmsd.append(np.sqrt(np.mean(np.sum(dist_check[dihedral_start:dihedral_end]**2))))
        
        pes = (
            energy
            + np.matmul(dist_check.T, grad)
            + 0.5 * np.linalg.multi_dot([dist_check.T, hessian, dist_check])
        )

        dist_hessian_eff = np.matmul(dist_check.T, hessian)
        if not self.use_cosine_dihedral:

            grad *= chain
            dist_hessian_eff *= chain


        pes_prime = np.matmul(self.impes_coordinate.b_matrix.T, (grad + dist_hessian_eff)).reshape(natm, 3)
        return pes, pes_prime, (grad + dist_hessian_eff)

    def _iter_imp_internal_coordinate_rows(self, data_point):
        """
        Yields important internal-coordinate entries for a datapoint.

        Returns
        -------
        section : str
            One of 'bonds', 'angles', 'dihedrals', 'impropers'.
        coord : tuple[int, ...]
            Atom tuple defining the coordinate.
        idx : int
            Row index in self.impes_coordinate.z_matrix.

        Notes
        -----
        This implementation assumes use_cosine_dihedral=False for proper
        dihedral gates. That matches the current default and avoids ambiguity
        caused by duplicated cosine/sine rows.
        """

        imp = getattr(data_point, "imp_int_coordinates", None)
        if not isinstance(imp, dict):
            return

        zmat = [tuple(int(x) for x in z) for z in self.impes_coordinate.z_matrix]

        for section in ("bonds", "angles", "dihedrals", "impropers"):
            for coord_raw in imp.get(section, []):
                coord = tuple(int(x) for x in coord_raw)

                try:
                    idx = zmat.index(coord)
                except ValueError:
                    # Important coordinate is absent from the active z-matrix.
                    # This should not happen for a consistent database, but skip
                    # defensively rather than killing interpolation.
                    continue

                yield section, coord, idx

            
    def _imp_coordinate_sigma(self, section, idx, q_ref):
        """
        Returns the sigma value in the same coordinate representation as the
        internal coordinate row idx.

        For inverse bonds, q = 1/r, so a physical bond tolerance sigma_r is
        converted into q-space via sigma_q = sigma_r / r_ref^2.
        """

        sigma_bond_bohr = (
            float(self.tc_imp_bond_sigma_angstrom) / bohr_in_angstrom()
        )

        if section == "bonds":
            if self.use_inverse_bond_length:
                # q_ref is 1/r in bohr^-1.
                r_ref = 1.0 / max(abs(float(q_ref[idx])), self.tc_imp_min_sigma)
                sigma = sigma_bond_bohr / max(r_ref * r_ref, self.tc_imp_min_sigma)
            else:
                sigma = sigma_bond_bohr

        elif section == "angles":
            sigma = np.deg2rad(float(self.tc_imp_angle_sigma_degrees))

        elif section == "dihedrals":
            sigma = np.sin(np.deg2rad(float(self.tc_imp_dihedral_sigma_degrees)))

        elif section == "impropers":
            sigma = np.sin(np.deg2rad(float(self.tc_imp_improper_sigma_degrees)))


        else:
            raise ValueError(f"Unknown important-coordinate section: {section}")

        return max(float(abs(sigma)), self.tc_imp_min_sigma)

    def _important_coordinate_gate_metric(self, data_point, active_dofs):
        """
        Computes the dimensionless important-coordinate distance A_imp and
        its Cartesian derivative.

        A_imp = sum_k beta_k z_k^2 / sum_k beta_k

        The derivative is with respect to the active Cartesian coordinates
        used by the Cartesian weight gradient.

        Returns
        -------
        A_imp : float
            Dimensionless mean squared important-coordinate displacement.

        grad_A_imp : ndarray, shape (n_active_dofs,)
            dA_imp/dX for active Cartesian degrees of freedom.
        """

        q_cur = np.asarray(
            self.impes_coordinate.internal_coordinates_values,
            dtype=np.float64,
        )
        q_ref = np.asarray(
            data_point.internal_coordinates_values,
            dtype=np.float64,
        )

        B = np.asarray(self.impes_coordinate.b_matrix, dtype=np.float64)

        inv_sqrt_masses = getattr(self.impes_coordinate, "inv_sqrt_masses", None)
        if inv_sqrt_masses is not None:
            inv_sqrt_active = np.asarray(inv_sqrt_masses, dtype=np.float64)[active_dofs]
        else:
            inv_sqrt_active = None

        A_num = 0.0
        beta_sum = 0.0
        grad_A_num = np.zeros(active_dofs.size, dtype=np.float64)

        for section, coord, idx in self._iter_imp_internal_coordinate_rows(data_point):
            beta = 1.0

            dq_raw = float(q_cur[idx] - q_ref[idx])
            sigma = self._imp_coordinate_sigma(section, idx, q_ref)

            # B row is mass-weighted if inv_sqrt_masses was set on the datapoint.
            # Convert it back to derivative with respect to physical Cartesian
            # coordinates, matching the existing target-customized code path.
            b_row = B[idx, active_dofs].copy()
            if inv_sqrt_active is not None:
                b_row = b_row / inv_sqrt_active

            if section == "dihedrals":
                # Periodic, smooth metric.
                d = self._principal_torsion_delta(dq_raw)
                z = np.sin(d) / sigma
                dz_dx = np.cos(d) * b_row / sigma

            elif section == "impropers":
                d = self._principal_torsion_delta(dq_raw)
                z = np.sin(d) / sigma
                dz_dx = np.cos(d) * b_row / sigma


            elif section == "angles":
                z = dq_raw / sigma
                dz_dx = b_row / sigma

            elif section == "bonds":
                z = dq_raw / sigma
                dz_dx = b_row / sigma

            else:
                continue

            A_num += beta * z * z
            grad_A_num += beta * 2.0 * z * dz_dx
            beta_sum += beta

        if beta_sum <= 0.0:
            return 0.0, np.zeros(active_dofs.size, dtype=np.float64)

        A_imp = A_num / beta_sum
        grad_A_imp = grad_A_num / beta_sum

        return float(A_imp), grad_A_imp
    
    def _apply_imp_coordinate_penalty_gate(
            self,
            denominator_base,
            raw_weight_gradient_base,
            A_imp,
            grad_A_imp,
        ):
            """
            Applies the important-coordinate penalty gate to a raw Shepard weight.

            Base:
                w = 1 / denominator_base

            Gate:
                G = exp(-lambda * A_imp)

            Modified raw weight:
                w_mod = w * G

            Gradient:
                grad(w_mod) = G grad(w) + w grad(G)
                            = G grad(w) - lambda w G grad(A_imp)

            Returns
            -------
            denominator_mod : float
                1 / w_mod, so the surrounding code can still use
                raw_weight = 1 / denominator.

            raw_weight_gradient_mod : ndarray
                Gradient of w_mod with respect to active Cartesian coordinates.
            """

            if not np.isfinite(denominator_base) or denominator_base <= 0.0:
                return np.inf, np.zeros_like(raw_weight_gradient_base)

            lam = float(self.tc_imp_gate_lambda)
            exponent = min(lam * float(A_imp), float(self.tc_imp_gate_exp_clip))

            gate = np.exp(-exponent)
            w_base = 1.0 / denominator_base

            grad_base_flat = np.asarray(raw_weight_gradient_base, dtype=np.float64).reshape(-1)
            grad_A_flat = np.asarray(grad_A_imp, dtype=np.float64).reshape(-1)

            grad_gate = -lam * gate * grad_A_flat

            w_mod = w_base * gate
            grad_w_mod = gate * grad_base_flat + w_base * grad_gate

            if w_mod <= 1.0e-300:
                return np.inf, np.zeros_like(raw_weight_gradient_base)

            denominator_mod = 1.0 / w_mod

            return denominator_mod, grad_w_mod.reshape(raw_weight_gradient_base.shape)




    def trust_radius_weight_gradient(self, confidence_radius, distance):
        
        R = confidence_radius
        D = distance**2
        
        p = float(self.exponent_p)
        q = float(self.exponent_q)
        
        # Dimensionless ratio u = (d/R)^2
        u = D / (R**2)
        
        # Early Exit: If the distance is massively outside the confidence radius, 
        # the trust weight gradient flattens to absolute zero.
        if np.isscalar(u) and u > 1e6:
            return 0.0
            
        denom_base = u**p + u**q
        denom_base = np.clip(denom_base, a_min=1e-300, a_max=None)
        
        # Algebraically simplified numerator
        numerator = p * u**p + q * u**q
        
        trust_radius_weight_gradient = (2.0 / R) * (numerator / denom_base**2)
        
        return trust_radius_weight_gradient

    
    def trust_radius_weight_gradient_gradient(self, confidence_radius, distance, distance_vector):
        
        R = confidence_radius
        D = distance**2
        v = distance_vector
        
        p = float(self.exponent_p)
        q = float(self.exponent_q)
        
        u = D / (R**2)
        
        # Early Exit
        if np.isscalar(u) and u > 1e6:
            return np.zeros_like(v, dtype=np.float64)
            
        denom_base = u**p + u**q
        denom_base = np.clip(denom_base, a_min=1e-300, a_max=None)
        
        # Factored numerator components (the constants and R factors are moved to the multiplier)
        term1 = (p**2 * u**(p - 1.0) + q**2 * u**(q - 1.0)) / (denom_base**2)
        term2 = 2.0 * ((p * u**(p - 1.0) + q * u**(q - 1.0)) * (p * u**p + q * u**q)) / (denom_base**3)
        
        # Common scaling factor analytically pulled out of the original expressions
        multiplier = (4.0 * v) / (R**3)
        
        trust_radius_weight_gradient_gradient = multiplier * (term1 - term2)

        return trust_radius_weight_gradient_gradient
    

    def shepard_weight_gradient(self, distance_vector, distance, confidence_radius):
        
        R = confidence_radius
        D = distance**2
        v = distance_vector
        
        p = float(self.exponent_p)
        q = float(self.exponent_q)
        
        u = D / (R**2)
        
        # Early Exit: Return np.inf for the denominator and an array of zeros for the gradient
        if np.isscalar(u) and u > 1e6:
            return np.inf, np.zeros_like(v, dtype=np.float64)
            
        denom_base = u**p + u**q
        denom_base = np.clip(denom_base, a_min=1e-300, a_max=None)
        
        numerator = p * u**(p - 1.0) + q * u**(q - 1.0)
        
        weight_gradient = -2.0 * (numerator / (R**2 * denom_base**2)) * v
        
        return denom_base, weight_gradient
    

    def trust_radius_tc_weight_gradient(self, confidence_radius, distance, imp_int_coords_distance):
        
        R = confidence_radius
        # Incorporating the internal coordinates distance
        D = distance**2 + imp_int_coords_distance
        
        p = float(self.exponent_p)
        q = float(self.exponent_q)
        
        u = D / (R**2)
        
        # Early Exit
        if np.isscalar(u) and u > 1e6:
            return 0.0
            
        denom_base = u**p + u**q
        denom_base = np.clip(denom_base, a_min=1e-300, a_max=None)
        
        numerator = p * u**p + q * u**q
        
        trust_radius_weight_gradient = (2.0 / R) * (numerator / denom_base**2)
        
        return trust_radius_weight_gradient

    def trust_radius_tc_weight_gradient_gradient(self, confidence_radius, distance, distance_vector, imp_int_coords_distance, imp_int_coords_dist_derivative):
            
            # 1. Define base variables
            R = confidence_radius
            D = distance**2 + imp_int_coords_distance
            v = imp_int_coords_dist_derivative + 2.0 * distance_vector.reshape(-1)
            
            p = float(self.exponent_p)
            q = float(self.exponent_q)
            
            # 2. Compute the dimensionless ratio
            u = D / (R**2)
            
            # --- THE FIX: Early Exit for massive ratios ---
            # If u is so large that u**(3 * max(p,q)) exceeds ~1e308, float64 overflows.
            # Physically, the weight gradient at this immense relative distance is practically 0.0.
            if np.isscalar(u) and u > 1e6:
                return np.zeros((distance_vector.shape[0], 3), dtype=np.float64)
            
            # 3. Base denominator
            denom_base = u**p + u**q
            denom_base = np.clip(denom_base, a_min=1e-300, a_max=None)
            
            # 4. Compute the grouped terms 
            term1 = (p**2 * u**(p - 1.0) + q**2 * u**(q - 1.0)) / (denom_base**2)
            term2 = 2.0 * ((p * u**(p - 1.0) + q * u**(q - 1.0)) * (p * u**p + q * u**q)) / (denom_base**3)
            
            # 5. Apply the common multiplier
            multiplier = (2.0 * v) / (R**3)
            
            trust_radius_weight_gradient_gradient = multiplier * (term1 - term2)

            return trust_radius_weight_gradient_gradient.reshape(distance_vector.shape[0], 3)

    def YM_target_customized_shepard_weight_gradient(self, distance_vector, distance, confidence_radius, imp_int_coords_distance, imp_int_coordinate_derivative):
        
        R = confidence_radius
        D = distance**2 + imp_int_coords_distance
        v = 2.0 * distance_vector.reshape(-1) + imp_int_coordinate_derivative
        
        p = float(self.exponent_p)
        q = float(self.exponent_q)
        
        u = D / (R**2)
        
        # Early Exit
        if np.isscalar(u) and u > 1e6:
            return np.inf, np.zeros((distance_vector.shape[0], 3), dtype=np.float64)
            
        denom_base = u**p + u**q
        denom_base = np.clip(denom_base, a_min=1e-300, a_max=None)
        
        numerator = p * u**(p - 1.0) + q * u**(q - 1.0)
        
        weight_gradient = -1.0 * (numerator / (R**2 * denom_base**2)) * v
   
        return denom_base, weight_gradient.reshape(distance_vector.shape[0], 3)
    

    def VL_target_customized_shepard_weight_gradient(self, distance_vector, distance, confidence_radius, dq, dq_dx):
        
        R = confidence_radius
        D = distance**2
        
        p = float(self.exponent_p)
        q = float(self.exponent_q)
        
        u = D / (R**2)

        # Early Exit
        if np.isscalar(u) and u > 1e8:
            return np.inf, np.zeros((distance_vector.shape[0], 3), dtype=np.float64)
            
        denom_base = u**p + u**q
        denom_base = np.clip(denom_base, a_min=1e-300, a_max=None)
        
        numerator = p * u**(p - 1.0) + q * u**(q - 1.0)
        
        # Original gradient logic simplified with factored terms
        weight_gradient = -2.0 * (numerator / (R**2 * denom_base**2)) * distance_vector
        
        final_denominator = denom_base * np.exp(self.alpha * dq) 
        
        # 1.0/denominator is handled safely here because denom_base is clipped above 0
        final_weight_gradient = (weight_gradient.reshape(-1) * np.exp(-self.alpha * dq) - 1.0/denom_base * self.alpha * np.exp(-self.alpha * dq) * dq_dx).reshape(-1, 3)

        return final_denominator, final_weight_gradient
    
    
    def calculate_translation_coordinates(self, given_coordinates):
        """Center the molecule by translating its geometric center to (0, 0, 0)."""

        center = np.mean(given_coordinates, axis=0)
        translated_coordinates = given_coordinates - center



        return translated_coordinates

    def internal_distance(self, data_point, store_alpha_gradients=True):
        """Computes a Shepard weight denominator using internal-coordinate distance.

        The structure is analogous to ``cartesian_distance``, but the base distance
        metric is constructed from internal-coordinate displacements and mapped to
        Cartesian derivatives via Wilson B-matrix rows.
        """

        reference_coordinates = self.impes_coordinate.cartesian_coordinates.copy()
        natms = reference_coordinates.shape[0]

        active_atoms = np.delete(np.arange(natms), self.symmetry_information[4])
        active_dofs = np.array([3 * a + i for a in active_atoms for i in range(3)], dtype=np.int64)

        q_current = np.asarray(self.impes_coordinate.internal_coordinates_values, dtype=np.float64)
        q_ref = np.asarray(data_point.internal_coordinates_values, dtype=np.float64)

        dq_raw = q_current - q_ref
        dq_eff = dq_raw.copy()
        chain = np.ones_like(dq_raw)

        bounds = self._get_internal_coordinate_partitions()
        dstart = int(bounds["dihedral_start"])
        dend = int(bounds["dihedral_end"])

        for idx in range(len(dq_raw)):
            if dstart <= idx < dend:  # proper dihedrals only
                d = self._principal_torsion_delta(dq_raw[idx])
                dq_eff[idx] = np.sin(d)
                chain[idx] = np.cos(d)
            else:  # bonds, angles, impropers -> linear
                dq_eff[idx] = dq_raw[idx]
                chain[idx] = 1.0

        grad_Dint = np.zeros(active_dofs.shape[0], dtype=np.float64)
        if self.impes_coordinate.inv_sqrt_masses is not None:
            inv_sqrt_masses = self.impes_coordinate.inv_sqrt_masses[active_dofs]
        else:
            inv_sqrt_masses = np.ones(active_dofs.shape[0], dtype=np.float64)

        for idx in range(len(dq_eff)):
            b_row = self.impes_coordinate.b_matrix[idx, active_dofs] / inv_sqrt_masses
            grad_Dint += 2.0 * dq_eff[idx] * chain[idx] * b_row

        Dint_sq = float(np.dot(dq_eff, dq_eff))
        if Dint_sq < 1e-8:
            Dint_sq = 1e-8
            grad_Dint[:] = 0.0

        distance = np.sqrt(Dint_sq)
        dihedral_dist = 0.0
        distance_vector_sub = np.zeros((active_atoms.shape[0], 3), dtype=np.float64)
        grad_s = np.zeros_like(distance_vector_sub)

        Dimp_sq = 0.0
        dq_dx = np.zeros_like(grad_Dint)

        if self.use_tc_weights and isinstance(data_point.imp_int_coordinates, dict):
            H_eff = data_point.internal_hessian
            eigval, V = np.linalg.eigh(H_eff)
            coord_curvature = np.sum(eigval * (V**2), axis=1)

            for _, entries in data_point.imp_int_coordinates.items():
                for element in entries:
                    idx = self.impes_coordinate.z_matrix.index(element)
                    delta_q = (
                        self.impes_coordinate.internal_coordinates_values[idx] -
                        data_point.internal_coordinates_values[idx]
                    )
                    b_row = self.impes_coordinate.b_matrix[idx, active_dofs] / inv_sqrt_masses

                    if len(element) == 4:
                        arg = coord_curvature[idx] * delta_q
                        dq_dx += (
                            coord_curvature[idx] *
                            2.0 *
                            b_row *
                            np.cos(arg) *
                            np.sin(arg)
                        )
                        Dimp_sq += np.sin(arg)**2
                    else:
                        dq_dx += coord_curvature[idx]**2 * 2.0 * b_row * delta_q
                        Dimp_sq += (coord_curvature[idx] * delta_q)**2

            if Dimp_sq < 1e-8:
                Dimp_sq = 1e-8
                dq_dx[:] = 0.0

        R = data_point.confidence_radius
        p = float(self.exponent_p)
        q = float(self.exponent_q)
        u = Dint_sq / (R**2)

        if np.isscalar(u) and u > 1e6:
            denom_base = np.inf
            base_weight_grad = np.zeros_like(grad_Dint, dtype=np.float64)
        else:
            denom_base = u**p + u**q
            denom_base = np.clip(denom_base, a_min=1e-300, a_max=None)
            numerator = p * u**(p - 1.0) + q * u**(q - 1.0)
            base_weight_grad = -1.0 * (numerator / (R**2 * denom_base**2)) * grad_Dint

        if self.use_tc_weights:
            exp_decay = np.exp(-self.alpha * Dimp_sq)
            denominator = denom_base * np.exp(self.alpha * Dimp_sq)
            weight_grad_sub = (
                base_weight_grad * exp_decay -
                (self.alpha * exp_decay / denom_base) * dq_dx
            )
        else:
            denominator = denom_base
            weight_grad_sub = base_weight_grad

        weight_gradient_sub = weight_grad_sub.reshape(active_atoms.shape[0], 3)
        weight_gradient = np.zeros_like(reference_coordinates)
        weight_gradient[self.symmetry_information[3]] = weight_gradient_sub

        dw_dalpha_i = 0.0
        dw_dX_dalpha_i = np.zeros_like(weight_gradient_sub)

        if self.calc_optim_trust_radius:
            if self.use_tc_weights:
                dw_dalpha_i = self.trust_radius_tc_weight_gradient(
                    data_point.confidence_radius,
                    distance,
                    Dimp_sq)
                dw_dX_dalpha_i = self.trust_radius_tc_weight_gradient_gradient(
                    data_point.confidence_radius,
                    distance,
                    distance_vector_sub,
                    Dimp_sq,
                    dq_dx)
            else:
                dw_dalpha_i = self.trust_radius_tc_weight_gradient(
                    data_point.confidence_radius,
                    0.0,
                    Dint_sq)
                dw_dX_dalpha_i = self.trust_radius_tc_weight_gradient_gradient(
                    data_point.confidence_radius,
                    0.0,
                    distance_vector_sub,
                    Dint_sq,
                    grad_Dint)

            if store_alpha_gradients:
                self.dw_dalpha_list.append(dw_dalpha_i)
                dw_dX_dalpha_i_full = np.zeros_like(reference_coordinates)
                dw_dX_dalpha_i_full[self.symmetry_information[3]] = dw_dX_dalpha_i
                self.dw_dX_dalpha_list.append(dw_dX_dalpha_i_full)

        return (
            distance,
            dihedral_dist,
            denominator,
            weight_gradient,
            distance_vector_sub,
            grad_s,
            dw_dalpha_i,
            dw_dX_dalpha_i,
        )

    def cartesian_distance(self, data_point, store_alpha_gradients=True):
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
        active_dofs = np.array([3*a + i for a in active_atoms for i in range(3)])

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

        
        dihedral_dist = 0.0
        if distance < 1e-8:
            distance = 1e-8
            distance_vector_sub[:] = 0
      
        dw_dalhpa_i = 0
        dw_dX_dalpha_i = 0

        if self.interpolation_type == 'shepard':
            
            denominator_base, weight_gradient_sub_base = self.shepard_weight_gradient(
                distance_vector_sub,
                distance,
                data_point.confidence_radius,
            )
            if self.use_tc_weights:
                A_imp, grad_A_imp = self._important_coordinate_gate_metric(
                    data_point,
                    active_dofs,
                )

                denominator_imp_coord, weight_gradient_sub_imp_coord = (
                    self._apply_imp_coordinate_penalty_gate(
                        denominator_base,
                        weight_gradient_sub_base,
                        A_imp,
                        grad_A_imp,
                    )
                )
            else:
                denominator_imp_coord = denominator_base
                weight_gradient_sub_imp_coord = weight_gradient_sub_base

            if self.calc_optim_trust_radius:
                if self.use_tc_weights:
                    A_imp, grad_A_imp = self._important_coordinate_gate_metric(
                        data_point,
                        active_dofs,
                    )

                    lam = float(self.tc_imp_gate_lambda)
                    gate = np.exp(-min(lam * A_imp, self.tc_imp_gate_exp_clip))
                    grad_gate = -lam * gate * grad_A_imp

                    dw_base_dR = self.trust_radius_weight_gradient(
                        data_point.confidence_radius,
                        distance,
                    )

                    dgradw_base_dR = self.trust_radius_weight_gradient_gradient(
                        data_point.confidence_radius,
                        distance,
                        distance_vector_sub,
                    )

                    dw_dalhpa_i = gate * dw_base_dR

                    dw_dX_dalpha_i = (
                        gate * dgradw_base_dR.reshape(-1)
                        + grad_gate * dw_base_dR
                    ).reshape(distance_vector_sub.shape)
                else:
                    dw_dalhpa_i = self.trust_radius_weight_gradient(
                        data_point.confidence_radius,
                        distance,
                    )

                    dw_dX_dalpha_i = self.trust_radius_weight_gradient_gradient(
                        data_point.confidence_radius,
                        distance,
                        distance_vector_sub,
                    )

                if store_alpha_gradients:
                    self.dw_dalpha_list.append(dw_dalhpa_i)
                    dw_dX_dalpha_i_full = np.zeros_like(reference_coordinates)
                    dw_dX_dalpha_i_full[self.symmetry_information[3]] = dw_dX_dalpha_i
                    self.dw_dX_dalpha_list.append(dw_dX_dalpha_i_full)
                    
        else:
            errtxt = "Unrecognized interpolation type: "
            errtxt += self.interpolation_type
            raise ValueError(errtxt)
                
        # print(weight_gradient_sub, weight_gradient_sub_imp_coord, denominator_imp_coord, denominator)
        weight_gradient = np.zeros_like(reference_coordinates)   # (natms,3)
        weight_gradient[self.symmetry_information[3]] = weight_gradient_sub_imp_coord

        return distance, dihedral_dist, denominator_imp_coord, weight_gradient, distance_vector_sub, grad_s, dw_dalhpa_i, dw_dX_dalpha_i

    def _effective_dimension(self, scores, eps=1e-12):
        scores = np.asarray(scores, dtype=float)
        scores = np.maximum(scores, 0.0)
        s = scores.sum()
        if s < eps:
            return 0.0
        p = scores / s
        return 1.0 / (np.sum(p * p) + eps)
    

    def determine_important_internal_coordinates(
        self,
        qm_gradient,
        molecule,
        z_matrix,
        datapoints,
    ):
        """
        Error-source-focused selection scheme.

        Philosophy
        ----------
        This version explicitly separates three questions:

        1) Where does the local model fail?
        -> response score from |g_qm - g_model|

        2) Which displaced coordinates drive that failure?
        -> source blame score from a leave-one-source-out test

        3) Which coordinates should be constrained together?
        -> build local coupled blocks from a blame-derived pair graph

        Proper dihedrals:
            wrapped and represented with dq_eff = sin(dq_raw), chain = cos(dq_raw)

        Impropers:
            wrapped to principal interval if needed, but kept linear:
            dq_eff = dq_raw, chain = 1
z
        """
        # ------------------------------------------------------------------
        # Hyperparameters
        # ------------------------------------------------------------------
        selection_rule = "coverage"   # 'relative' | 'coverage' | 'topk'
        relative_threshold = 0.10
        coverage_mass = 0.80
        topk = None

        # final global score = source + response + support - helpful
        w_source = 0.60
        w_response = 0.30
        w_support = 0.10
        helpful_penalty = 0.15

        # local support
        support_top_frac = 0.20

        # block construction / ranking
        max_anchor_candidates = 8
        pair_weight = 0.35
        support_weight = 0.25

        component_score_z = 2.0
        component_pair_z = 1.5
        component_coverage_mass = 0.90
        component_size_penalty = 0.01
        max_component_size = None

        # final selection from best block
        max_constraints_to_return = None
        singleton_dominance = 1.35
        singleton_source_frac = 0.85
        singleton_response_frac = 0.85
        secondary_score_frac = 0.25
        min_secondary_pair_frac = 0.15

        eps = 1e-12

        coords_flat, kinds_flat = self._flatten_internal_coordinates(z_matrix)
        N = len(coords_flat)
        if max_component_size is None:
            max_component_size = min(20, max(5, int(np.ceil(np.sqrt(N)))))

        if N == 0:
            return [], [], []

        constraints_to_exclude = []
        per_dp_results = []

        # ------------------------------------------------------------------
        # Exclude near-linear angles from direct constraint selection
        # (but they still remain in diagnostics)
        # ------------------------------------------------------------------
        for element in z_matrix.get("angles", []):
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

        # ------------------------------------------------------------------
        # Per-datapoint diagnostic pass
        # ------------------------------------------------------------------
        for datapoint, total_weight_contribution in datapoints:
            q_dp = np.asarray(datapoint.internal_coordinates_values, dtype=float)

            dq_raw, dq_eff, chain = self._compute_internal_differences(
                q_current=q_current,
                q_ref=q_dp,
                z_matrix=z_matrix,
                symmetry_information=self.symmetry_information,
            )

            pred_qm_G_int = self.transform_gradient_to_internal_coordinates(
                molecule, qm_gradient_mw.reshape(qm_gradient.shape), self.impes_coordinate.b_matrix
            )

            diag_result = self.compute_internal_gradient_diagnostics_runtime_model(
                datapoint=datapoint,
                molecule=molecule,
                z_matrix=z_matrix,
                q_current=q_current,
                q_ref=q_dp,
                g_qm_int=np.asarray(pred_qm_G_int, dtype=float),
            )

            response_score = np.asarray(diag_result["response_score"], dtype=float)
            source_blame_score = np.asarray(diag_result["source_blame_score"], dtype=float)
            helpful_source_score = np.asarray(diag_result["helpful_source_score"], dtype=float)
            pair_score = np.asarray(diag_result["pair_score"], dtype=float)
            rho = np.asarray(diag_result["rho"], dtype=float)

            # local ranking used only for support counting
            local_coord_score = self._normalize_nonnegative(
                0.65 * source_blame_score + 0.35 * response_score
            )

            norm_displacement = self._normalized_internal_displacement(
                dq_eff=dq_eff,
                q_ref=q_dp,
                coords_flat=coords_flat,
                kinds_flat=kinds_flat,
            )
            
            support_mask = np.zeros(N, dtype=float)
            if local_coord_score.sum() > eps:
                order_dp = np.argsort(-local_coord_score)
                k_support = max(1, int(np.ceil(support_top_frac * N)))
                support_mask[order_dp[:k_support]] = 1.0

            per_dp_results.append({
                "datapoint": datapoint,
                "dp_weight": float(total_weight_contribution),
                "dq_raw": dq_raw.copy(),
                "dq_eff": dq_eff.copy(),
                "chain": chain.copy(),
                "rho": rho.copy(),
                "response_score": response_score.copy(),
                "source_blame_score": source_blame_score.copy(),
                "helpful_source_score": helpful_source_score.copy(),
                "pair_score": pair_score.copy(),
                "local_coord_score": local_coord_score.copy(),
                "support_mask": support_mask.copy(),
                "rows": diag_result["rows"],
                "norm_displacement": norm_displacement.copy(),
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
        # Aggregate global response/source/pair information
        # ------------------------------------------------------------------
        global_response = np.zeros(N, dtype=float)
        global_source_blame = np.zeros(N, dtype=float)
        global_helpful = np.zeros(N, dtype=float)
        global_support = np.zeros(N, dtype=float)
        global_rho = np.zeros(N, dtype=float)
        global_pair = np.zeros((N, N), dtype=float)
        global_displacement = np.zeros(N, dtype=float)

        for idx in keep_indices:
            wk = dp_weights[idx]
            r = per_dp_results[idx]

            global_response += wk * r["response_score"]
            global_source_blame += wk * r["source_blame_score"]
            global_helpful += wk * r["helpful_source_score"]
            global_support += wk * r["support_mask"]
            global_rho += wk * r["rho"]
            global_pair += wk * r["pair_score"]
            global_displacement += wk * r["norm_displacement"]

        global_response = self._normalize_nonnegative(global_response)
        global_source_blame = self._normalize_nonnegative(global_source_blame)
        global_helpful = self._normalize_nonnegative(global_helpful)

        if global_support.max() > eps:
            global_support = global_support / (global_support.max() + eps)

        if np.max(global_pair) > eps:
            global_pair = global_pair / (np.max(global_pair) + eps)

        raw_global_score = (
            w_source * global_source_blame
            + w_response * global_response
            + w_support * global_support
            - helpful_penalty * global_helpful
        )
        raw_global_score = np.maximum(raw_global_score, 0.0)
        global_coord_score = self._normalize_nonnegative(raw_global_score)

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
        # Build and rank coupled repair blocks
        # ------------------------------------------------------------------
        N_coords = N

        d_eff = self._effective_dimension(global_coord_score)

        max_constraints_to_return = int(
            np.clip(
                np.ceil(1.25 * d_eff),
                1,
                min(12, max(3, int(np.ceil(np.sqrt(N_coords)))))
            )
        )

        ranked_blocks = []
        seen_blocks = set()

        components = self._build_faulty_components(
            global_coord_score=global_coord_score,
            global_pair=global_pair,
            coords_flat=coords_flat,
            kinds_flat=kinds_flat,
            constraints_to_exclude=constraints_to_exclude,
            score_z=component_score_z,
            pair_z=component_pair_z,
        )

        # Fallback: if the graph builder finds nothing, use the best coordinate.
        if not components:
            best_idx = int(np.argmax(global_coord_score))
            components = [[best_idx]]

        for component in components:
            # Defensive filtering.
            component = [
                int(i) for i in component
                if tuple(int(x) for x in coords_flat[i]) not in constraints_to_exclude
            ]

            component = [
                int(i) for i in component
            ]

            if not component:
                continue

            # Steering score chooses the anchor inside the component.
            steering_score = {
                i: float(0.5 * global_source_blame[i] + 0.5 * global_response[i])
                for i in component
            }

            anchor_idx = max(component, key=lambda i: steering_score[i])
            anchor_coord = tuple(int(x) for x in coords_flat[anchor_idx])

            # Sort component coordinates by a mixture of individual score and coupling to anchor.
            pair_to_anchor = np.asarray(global_pair[anchor_idx], dtype=float)

            component_sorted = sorted(
                component,
                key=lambda i: (
                    0.55 * float(global_coord_score[i])
                    + 0.30 * float(steering_score[i])
                    + 0.15 * float(pair_to_anchor[i])
                ),
                reverse=True,
            )
            component_sorted = [anchor_idx] + [
                int(i) for i in component_sorted if int(i) != int(anchor_idx)
            ]
            # Trim oversized components by score coverage.
            total_component_score = float(np.sum(global_coord_score[component_sorted]))
            total_component_score = float(np.sum(global_coord_score[component_sorted]))
            block_indices = [int(anchor_idx)]
            cumulative_score = float(global_coord_score[anchor_idx])

            for i in component_sorted:
                i = int(i)
                if i == anchor_idx:
                    continue

                if len(block_indices) >= max_component_size:
                    break

                if (
                    total_component_score > eps
                    and cumulative_score >= component_coverage_mass * total_component_score
                    and len(block_indices) >= 1
                ):
                    break

                block_indices.append(i)
                cumulative_score += float(global_coord_score[i])

            key = tuple(sorted(block_indices))
            if key in seen_blocks:
                continue
            seen_blocks.add(key)

            coord_mass = float(np.sum(global_coord_score[block_indices]))
            support_mass = float(np.mean(global_support[block_indices])) if block_indices else 0.0

            # Use pair density, not raw pair sum.
            # Raw pair sum grows approximately as |block|^2 and would bias toward large components.
            pair_sum = 0.0
            n_pairs = 0

            for a in range(len(block_indices)):
                for b in range(a + 1, len(block_indices)):
                    i = block_indices[a]
                    j = block_indices[b]
                    pair_sum += float(global_pair[i, j])
                    n_pairs += 1

            pair_density = pair_sum / max(1, n_pairs)

            acq_score = (
                coord_mass
                + pair_weight * pair_density
                + support_weight * support_mass
                - component_size_penalty * len(block_indices)
            )

            block_coords = [tuple(int(x) for x in coords_flat[i]) for i in block_indices]
            block_types = [kinds_flat[i] for i in block_indices]

            ranked_blocks.append({
                "anchor_idx": int(anchor_idx),
                "anchor_coord": anchor_coord,
                "anchor_type": kinds_flat[anchor_idx],
                "block_indices": [int(i) for i in block_indices],
                "block_coords": block_coords,
                "block_types": block_types,
                "coord_mass": float(coord_mass),
                "pair_mass": float(pair_density),   # now normalized pair density
                "pair_sum": float(pair_sum),
                "support_mass": float(support_mass),
                "acq_score": float(acq_score),
            })

        ranked_blocks = sorted(ranked_blocks, key=lambda x: x["acq_score"], reverse=True)

        if not ranked_blocks:
            return [], [], []

        best_block = ranked_blocks[0]

        candidate_constraints, primary_constraint = self._select_constraints_from_block(
            block_indices=best_block["block_indices"],
            global_coord_score=global_coord_score,
            global_source_blame=global_source_blame,
            global_response_score=global_response,
            global_pair=global_pair,
            z_matrix=coords_flat,
            max_constraints=max_constraints_to_return,
            singleton_dominance=singleton_dominance,
            singleton_source_frac=singleton_source_frac,
            singleton_response_frac=singleton_response_frac,
            secondary_score_frac=secondary_score_frac,
            min_secondary_pair_frac=min_secondary_pair_frac,
        )

        fallback_constraints = self._select_fallback_locality_constraints(
            block_indices=best_block["block_indices"],
            candidate_constraints=candidate_constraints,
            primary_constraint=primary_constraint,
            coords_flat=coords_flat,
            kinds_flat=kinds_flat,
            constraints_to_exclude=constraints_to_exclude,
            global_coord_score=global_coord_score,
            global_pair=global_pair,
            global_displacement=global_displacement,
            max_constraints=max(2, min(8, max_constraints_to_return)),
        )

        # ------------------------------------------------------------------
        # Debug printout
        # ------------------------------------------------------------------
        print("\n=== Global coordinate scores (error-source focused) ===")
        print(
            f"{'rank':>4} {'idx':>4} {'type':>10} {'coord':>18} "
            f"{'score':>12} {'source':>12} {'response':>12} "
            f"{'support':>12} {'rho':>12} {'helpful':>12}"
        )
        print("-" * 122)
        for rank, idx in enumerate(np.argsort(-global_coord_score)[:12]):
            coord = tuple(int(x) for x in coords_flat[idx])
            print(
                f"{rank:4d} {idx:4d} {kinds_flat[idx]:>10} {self.format_coord(coord):>18} "
                f"{global_coord_score[idx]:12.6e} "
                f"{global_source_blame[idx]:12.6e} "
                f"{global_response[idx]:12.6e} "
                f"{global_support[idx]:12.6e} "
                f"{global_rho[idx]:12.6e} "
                f"{global_helpful[idx]:12.6e}"
            )

        print("\n=== Ranked candidate repair blocks ===")
        print(
            f"{'rank':>4} {'anchor':>18} {'type':>10} {'|block|':>8} "
            f"{'coord':>12} {'pair':>12} {'support':>12} {'acq':>12}"
        )
        print("-" * 102)
        for rank, blk in enumerate(ranked_blocks[:10]):
            print(
                f"{rank:4d} "
                f"{self.format_coord(blk['anchor_coord']):>18} "
                f"{blk['anchor_type']:>10} "
                f"{len(blk['block_indices']):8d} "
                f"{blk['coord_mass']:12.6e} "
                f"{blk['pair_mass']:12.6e} "
                f"{blk['support_mass']:12.6e} "
                f"{blk['acq_score']:12.6e}"
            )
            print(f"      block coords: {blk['block_coords']}")
            print(f"      block types : {blk['block_types']}")

        print("\nPrimary constraint:", primary_constraint)
        print("Candidate constraints:", candidate_constraints)

        return primary_constraint, candidate_constraints, best_block["block_coords"], fallback_constraints

    # Helper methods

    def _robust_positive_threshold(
        self,
        values: np.ndarray,
        z: float = 2.0,
        eps: float = 1e-12,
    ) -> float:
        """
        Robust threshold for positive values using median + z * MAD.

        This is more size-consistent than fixed fractions because the threshold
        adapts to the actual score distribution.
        """
        values = np.asarray(values, dtype=float)
        positive = values[values > eps]

        if positive.size == 0:
            return np.inf

        med = np.median(positive)
        mad = np.median(np.abs(positive - med)) + eps

        return float(med + z * 1.4826 * mad)


    def _is_coord_excluded_for_constraint(
        self,
        idx: int,
        coords_flat: List[Tuple[int, ...]],
        kinds_flat: List[str],
        constraints_to_exclude: List[Tuple[int, ...]],
    ) -> bool:
        """
        Centralized exclusion logic for direct constraint coordinates.
        """
        coord = tuple(int(x) for x in coords_flat[idx])

        if coord in constraints_to_exclude:
            return True

        return False


    def _build_faulty_components(
        self,
        *,
        global_coord_score: np.ndarray,
        global_pair: np.ndarray,
        coords_flat: List[Tuple[int, ...]],
        kinds_flat: List[str],
        constraints_to_exclude: List[Tuple[int, ...]],
        score_z: float = 2.0,
        pair_z: float = 1.5,
        min_atom_overlap: float = 0.25,
        bridge_score_frac: float = 0.20,
        bridge_score_z: float = 0.5,
        eps: float = 1e-12,
    ) -> List[List[int]]:
        """
        Build connected faulty coordinate components using core and bridge nodes.

        Core nodes:
            coordinates with clearly high error-source score.

        Bridge nodes:
            coordinates with moderate score that are chemically or blame-coupled
            to a core node.

        This avoids missing coordinates that are necessary to lock the failed
        geometry during constrained optimization but do not have top direct
        gradient-error scores.
        """
        global_coord_score = np.asarray(global_coord_score, dtype=float)
        global_pair = np.asarray(global_pair, dtype=float)

        N = len(global_coord_score)

        if N == 0:
            return []

        score_max = float(np.max(global_coord_score))

        if score_max <= eps:
            return []

        # ------------------------------------------------------------
        # Pair threshold
        # ------------------------------------------------------------
        upper = global_pair[np.triu_indices(N, k=1)]
        pair_threshold = self._robust_positive_threshold(
            upper,
            z=pair_z,
            eps=eps,
        )

        # If the pair matrix is almost empty, avoid pair_threshold = inf
        # blocking all chemically local bridge connections.
        if not np.isfinite(pair_threshold):
            pair_threshold = np.inf

        # ------------------------------------------------------------
        # Core and bridge score thresholds
        # ------------------------------------------------------------
        core_threshold = self._robust_positive_threshold(
            global_coord_score,
            z=score_z,
            eps=eps,
        )

        if not np.isfinite(core_threshold):
            core_threshold = score_max

        bridge_threshold = max(
            bridge_score_frac * score_max,
            self._robust_positive_threshold(
                global_coord_score,
                z=bridge_score_z,
                eps=eps,
            ),
        )

        if not np.isfinite(bridge_threshold):
            bridge_threshold = bridge_score_frac * score_max

        # Make sure bridge threshold is lower than core threshold.
        bridge_threshold = min(bridge_threshold, 0.75 * core_threshold)

        # ------------------------------------------------------------
        # Collect core and bridge-candidate nodes
        # ------------------------------------------------------------
        core_nodes = []
        bridge_candidates = []

        for i in range(N):
            if self._is_coord_excluded_for_constraint(
                i,
                coords_flat,
                kinds_flat,
                constraints_to_exclude,
            ):
                continue

            score_i = float(global_coord_score[i])

            if score_i >= core_threshold:
                core_nodes.append(int(i))

            elif score_i >= bridge_threshold:
                bridge_candidates.append(int(i))

        # Always retain at least the best non-excluded coordinate as a core.
        if not core_nodes:
            order = np.argsort(-global_coord_score)

            for i in order:
                i = int(i)

                if self._is_coord_excluded_for_constraint(
                    i,
                    coords_flat,
                    kinds_flat,
                    constraints_to_exclude,
                ):
                    continue

                core_nodes = [i]
                break

        if not core_nodes:
            return []

        # ------------------------------------------------------------
        # Accept only bridge nodes connected to a core node
        # ------------------------------------------------------------
        accepted_bridges = []

        for b in bridge_candidates:
            atoms_b = set(coords_flat[b])

            for c in core_nodes:
                atoms_c = set(coords_flat[c])

                atom_overlap = len(atoms_b.intersection(atoms_c)) / max(
                    1,
                    len(atoms_b.union(atoms_c)),
                )

                pair_ok = float(global_pair[b, c]) >= pair_threshold

                # If the pair graph is sparse, atom overlap is still allowed
                # to keep chemically local bridge nodes.
                overlap_ok = atom_overlap >= min_atom_overlap

                if pair_ok or overlap_ok:
                    accepted_bridges.append(int(b))
                    break

        nodes = set(core_nodes + accepted_bridges)

        if not nodes:
            return []

        # ------------------------------------------------------------
        # Build graph between all accepted nodes
        # ------------------------------------------------------------
        adjacency = {i: set() for i in nodes}
        node_list = sorted(nodes)

        for a, i in enumerate(node_list):
            atoms_i = set(coords_flat[i])

            for j in node_list[a + 1:]:
                atoms_j = set(coords_flat[j])

                pair_ij = float(global_pair[i, j])

                atom_overlap = len(atoms_i.intersection(atoms_j)) / max(
                    1,
                    len(atoms_i.union(atoms_j)),
                )

                pair_ok = pair_ij >= pair_threshold

                overlap_ok = atom_overlap >= min_atom_overlap

                # Additional score compatibility condition:
                # avoid connecting very weak bridge nodes to everything.
                score_compat_ok = (
                    min(global_coord_score[i], global_coord_score[j])
                    >= bridge_threshold
                )

                if pair_ok or (overlap_ok and score_compat_ok):
                    adjacency[i].add(j)
                    adjacency[j].add(i)

        # ------------------------------------------------------------
        # Connected components
        # ------------------------------------------------------------
        components = []
        visited = set()

        for start in node_list:
            if start in visited:
                continue

            stack = [start]
            visited.add(start)
            component = []

            while stack:
                u = stack.pop()
                component.append(int(u))

                for v in adjacency[u]:
                    if v not in visited:
                        visited.add(v)
                        stack.append(v)

            components.append(sorted(component))

        # ------------------------------------------------------------
        # Rank components by score mass
        # ------------------------------------------------------------
        components = sorted(
            components,
            key=lambda comp: float(np.sum(global_coord_score[comp])),
            reverse=True,
        )

        return components

    def _flatten_internal_coordinates(
        self,
        z_matrix: Dict[str, List[Tuple[int, ...]]]
    ) -> Tuple[List[Tuple[int, ...]], List[str]]:
        """
        Flatten z_matrix into a single coordinate list with explicit kinds.

        Order is fixed and must remain consistent everywhere.
        """
        coords = []
        kinds = []

        for coord in z_matrix.get("bonds", []):
            coords.append(tuple(int(x) for x in coord))
            kinds.append("bond")

        for coord in z_matrix.get("angles", []):
            coords.append(tuple(int(x) for x in coord))
            kinds.append("angle")

        for coord in z_matrix.get("dihedrals", []):
            coords.append(tuple(int(x) for x in coord))
            kinds.append("dihedral")

        for coord in z_matrix.get("impropers", []):
            coords.append(tuple(int(x) for x in coord))
            kinds.append("improper")

        return coords, kinds


    def _normalize_nonnegative(self, x: np.ndarray, eps: float = 1e-12) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        x = np.maximum(x, 0.0)
        s = x.sum()
        if s < eps:
            return np.zeros_like(x)
        return x / s


    def _principal_torsion_delta(self, delta: float) -> float:
        """Wrap angular difference to (-pi, pi]."""
        return (delta + np.pi) % (2.0 * np.pi) - np.pi

    def _improper_dihedral_displacement(self, delta):
        """
        Improper dihedral displacement used in the local Taylor model:
            2 * tan(delta / 2)
        """
        d = self._principal_torsion_delta(delta)
        return 2.0 * np.tan(0.5 * d)


    def _improper_dihedral_chain(self, delta):
        """
        Chain factor for d[2*tan(delta/2)]/d(delta):
            sec^2(delta / 2)
        """
        d = self._principal_torsion_delta(delta)
        cos_half = np.cos(0.5 * d)
        return 1.0 / np.maximum(cos_half * cos_half, 1.0e-12)



    def coord_type(self, coord: Tuple[int, ...], kind: Optional[str] = None) -> str:
        """
        Prefer explicit kind if available. Fall back to tuple-length classification.
        """
        if kind is not None:
            return kind

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


    def _compute_internal_differences(
        self,
        q_current: np.ndarray,
        q_ref: np.ndarray,
        z_matrix: Dict[str, List[Tuple[int, ...]]],
        symmetry_information,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Proper dihedrals:
            dq_eff = sin(delta), chain = cos(delta)

        Impropers:
            dq_eff = 2*tan(delta/2), chain = sec^2(delta/2)

        Bonds/angles:
            dq_eff = raw delta, chain = 1
        """
        q_current = np.asarray(q_current, dtype=float)
        q_ref = np.asarray(q_ref, dtype=float)

        coords_flat, kinds_flat = self._flatten_internal_coordinates(z_matrix)
        N = len(coords_flat)

        dq_raw = np.zeros(N, dtype=float)
        dq_eff = np.zeros(N, dtype=float)
        chain = np.ones(N, dtype=float)

        for idx, (_, kind) in enumerate(zip(coords_flat, kinds_flat)):
            if kind == "dihedral":
                d = self._principal_torsion_delta(q_current[idx] - q_ref[idx])
                dq_raw[idx] = d
                dq_eff[idx] = np.sin(d)
                chain[idx] = np.cos(d)

            elif kind == "improper":
                d = self._principal_torsion_delta(q_current[idx] - q_ref[idx])
                dq_raw[idx] = d
                dq_eff[idx] = self._improper_dihedral_displacement(d)
                chain[idx] = self._improper_dihedral_chain(d)

            else:
                d = q_current[idx] - q_ref[idx]
                dq_raw[idx] = d
                dq_eff[idx] = d
                chain[idx] = 1.0

        return dq_raw, dq_eff, chain



    def _compute_partial_taylor_contributions(
        self,
        dq_eff: np.ndarray,
        chain: np.ndarray,
        H: np.ndarray,
        g0: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Per-coordinate split of Taylor energy and gradient contributions.

        Energy is decomposed in the transformed displacement dq_eff.
        Gradient follows the chain-weighted form used by the local model.
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

    def _predict_internal_gradient_runtime(self, datapoint, molecule, org_int_coords):
        """
        Internal gradient from the exact interpolation model used in runtime:
        includes rotor projectors, cluster assembly, mask-weighting, etc.
        """
        org_int_coords = np.asarray(org_int_coords, dtype=np.float64)

        # Avoid polluting outer diagnostics buffers.
        old_bond = self.bond_rmsd
        old_angle = self.angle_rmsd
        old_dihedral = self.dihedral_rmsd
        self.bond_rmsd = []
        self.angle_rmsd = []
        self.dihedral_rmsd = []

        try:
            _, pred_g_mw_cart, _ = self.compute_potential(datapoint, org_int_coords)
        finally:
            self.bond_rmsd = old_bond
            self.angle_rmsd = old_angle
            self.dihedral_rmsd = old_dihedral

        pred_g_int = self.transform_gradient_to_internal_coordinates(
            molecule,
            np.asarray(pred_g_mw_cart, dtype=np.float64),
            self.impes_coordinate.b_matrix,
        )
        return np.asarray(pred_g_int, dtype=np.float64).reshape(-1)

    def compute_internal_gradient_diagnostics_runtime_model(
        self,
        *,
        datapoint,
        molecule,
        z_matrix,
        q_current,
        q_ref,
        g_qm_int,
        max_sources=None,
        eps=1.0e-12,
    ):
        """
        Source-blame diagnostics that mirror the actual interpolation model.
        For each source j, remove its displacement and re-evaluate the full model.
        """
        coords_flat, kinds_flat = self._flatten_internal_coordinates(z_matrix)
        N = len(coords_flat)

        q_current = np.asarray(q_current, dtype=np.float64).reshape(-1)
        q_ref = np.asarray(q_ref, dtype=np.float64).reshape(-1)
        g_qm = np.asarray(g_qm_int, dtype=np.float64).reshape(-1)

        if q_current.size != N or q_ref.size != N or g_qm.size != N:
            raise ValueError("runtime diagnostics: size mismatch")

        dq_raw, dq_eff, chain = self._compute_internal_differences(
            q_current=q_current,
            q_ref=q_ref,
            z_matrix=z_matrix,
            symmetry_information=self.symmetry_information,
        )

        # Current model prediction (exact runtime path).
        g_model = self._predict_internal_gradient_runtime(datapoint, molecule, q_current)
        delta_g_err = g_qm - g_model
        abs_err = np.abs(delta_g_err)

        # Select which sources to test (optional speed cap).
        source_indices = np.arange(N, dtype=np.int64)
        if isinstance(max_sources, int) and 0 < max_sources < N:
            source_indices = np.argsort(-np.abs(dq_eff))[:max_sources]

        blame_matrix = np.zeros((N, N), dtype=np.float64)

        for j in source_indices:
            q_wo_j = q_current.copy()

            # Remove displacement of source j in the same coordinate chart.
            if kinds_flat[j] in ("dihedral", "improper"):
                d = self._principal_torsion_delta(q_current[j] - q_ref[j])
                q_wo_j[j] = q_current[j] - d
            else:
                q_wo_j[j] = q_ref[j]

            g_model_wo_j = self._predict_internal_gradient_runtime(datapoint, molecule, q_wo_j)
            err_wo_j = np.abs(g_qm - g_model_wo_j)

            # Positive => source j worsens row-i error.
            blame_matrix[:, j] = abs_err - err_wo_j

        harmful_matrix = np.maximum(blame_matrix, 0.0)
        helpful_matrix = np.maximum(-blame_matrix, 0.0)

        response_score = self._normalize_nonnegative(abs_err)

        err_weights = abs_err / (abs_err.sum() + eps)
        source_blame_score = self._normalize_nonnegative(harmful_matrix.T @ err_weights)
        helpful_source_score = self._normalize_nonnegative(helpful_matrix.T @ err_weights)

        pair_score = harmful_matrix + harmful_matrix.T
        np.fill_diagonal(pair_score, 0.0)
        if np.max(pair_score) > eps:
            pair_score /= (np.max(pair_score) + eps)

        rho = abs_err / (np.abs(g_model) + eps)

        rows = []
        for idx, (coord, kind) in enumerate(zip(coords_flat, kinds_flat)):
            rows.append({
                "idx": idx,
                "coord": tuple(int(x) for x in coord),
                "type": self.coord_type(tuple(coord), kind),
                "dq_raw": float(dq_raw[idx]),
                "dq_eff": float(dq_eff[idx]),
                "abs_dq": float(abs(dq_eff[idx])),
                "chain": float(chain[idx]),
                "g_model": float(g_model[idx]),
                "g_qm": float(g_qm[idx]),
                "delta_g_err": float(delta_g_err[idx]),
                "abs_delta_g_err": float(abs_err[idx]),
                "rho": float(rho[idx]),
                "source_blame_score": float(source_blame_score[idx]),
                "helpful_source_score": float(helpful_source_score[idx]),
                "response_score": float(response_score[idx]),
                # Keep compatibility keys expected by print helpers.
                "diag_resp": 0.0,
                "abs_diag_resp": 0.0,
                "offdiag_resp": 0.0,
                "abs_offdiag_resp": 0.0,
                "Hdq": 0.0,
                "abs_Hdq": 0.0,
                "pE_total": 0.0,
                "abs_pE_total": 0.0,
                "g0": 0.0,
            })

        by_abs_dq = sorted(rows, key=lambda r: r["abs_dq"], reverse=True)
        by_abs_Hdq = sorted(rows, key=lambda r: r["abs_Hdq"], reverse=True)
        by_abs_err = sorted(rows, key=lambda r: r["abs_delta_g_err"], reverse=True)
        by_rho = sorted(rows, key=lambda r: r["rho"], reverse=True)
        by_abs_pE = sorted(rows, key=lambda r: r["abs_pE_total"], reverse=True)
        by_source_blame = sorted(rows, key=lambda r: r["source_blame_score"], reverse=True)

        return {
            "rows": rows,
            "dq_eff": dq_eff,
            "chain": chain,
            "Hdq": np.zeros(N, dtype=np.float64),
            "diag_resp": np.zeros(N, dtype=np.float64),
            "offdiag_resp": np.zeros(N, dtype=np.float64),
            "g_model": g_model,
            "delta_g_err": delta_g_err,
            "rho": rho,
            "pE_total": np.zeros(N, dtype=np.float64),
            "source_contrib": blame_matrix,
            "blame_matrix": blame_matrix,
            "harmful_matrix": harmful_matrix,
            "helpful_matrix": helpful_matrix,
            "response_score": response_score,
            "source_blame_score": source_blame_score,
            "helpful_source_score": helpful_source_score,
            "pair_score": pair_score,
            "by_abs_dq": by_abs_dq,
            "by_abs_Hdq": by_abs_Hdq,
            "by_abs_err": by_abs_err,
            "by_rho": by_rho,
            "by_abs_pE": by_abs_pE,
            "by_source_blame": by_source_blame,
        }

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

        Core model:
            g_model = chain * (g0 + H @ dq_eff)

        Proper dihedrals:
            dq_eff = sin(dq_raw), chain = cos(dq_raw)

        Impropers:
            dq_eff = dq_raw, chain = 1

        Important new object:
            source_contrib[i,j] = chain[i] * H[i,j] * dq_eff[j]

        This is the modeled contribution from source coordinate j into response row i.

        Error-source score:
            err_i = g_qm[i] - g_model[i]
            blame_ij = |err_i| - |err_i + source_contrib[i,j]|

        If blame_ij > 0, removing source j would reduce the error in row i.
        """
        dq_raw = np.asarray(dq_raw, dtype=float).reshape(-1)
        H = np.asarray(H, dtype=float)
        g0 = np.asarray(g0, dtype=float).reshape(-1)
        g_qm = np.asarray(g_qm, dtype=float).reshape(-1)

        coords_flat, kinds_flat = self._flatten_internal_coordinates(z_matrix)
        N = len(coords_flat)

        if dq_raw.shape[0] != N:
            raise ValueError("dq_raw has wrong length")
        if g0.shape[0] != N:
            raise ValueError("g0 has wrong length")
        if g_qm.shape[0] != N:
            raise ValueError("g_qm has wrong length")
        if H.shape != (N, N):
            raise ValueError("H has wrong shape")

        dq_raw = np.array(dq_raw, dtype=float).reshape(-1)
        dq_eff = np.zeros(N, dtype=float)
        chain = np.ones(N, dtype=float)

        for idx, (_, kind) in enumerate(zip(coords_flat, kinds_flat)):
            if kind == "dihedral":
                d = self._principal_torsion_delta(dq_raw[idx])
                dq_raw[idx] = d
                dq_eff[idx] = np.sin(d)
                chain[idx] = np.cos(d)

            elif kind == "improper":
                d = self._principal_torsion_delta(dq_raw[idx])
                dq_raw[idx] = d
                dq_eff[idx] = self._improper_dihedral_displacement(d)
                chain[idx] = self._improper_dihedral_chain(d)

            else:
                dq_eff[idx] = dq_raw[idx]
                chain[idx] = 1.0


        Hdq_raw = H @ dq_eff
        Hdq = chain * Hdq_raw

        diag_resp_raw = np.diag(H) * dq_eff
        diag_resp = chain * diag_resp_raw
        offdiag_resp = Hdq - diag_resp

        g_model = chain * (g0 + Hdq_raw)
        delta_g_err = g_qm - g_model
        rho = np.abs(delta_g_err) / (np.abs(Hdq) + eps)

        partial_energies, partial_gradient = self._compute_partial_taylor_contributions(
            dq_eff=dq_eff,
            chain=chain,
            H=H,
            g0=g0,
        )


        # ------------------------------------------------------------------
        # New source-blame analysis
        # ------------------------------------------------------------------
        source_contrib = (chain[:, None] * H) * dq_eff[None, :]

        # removing source j changes row-i error from err_i to err_i + contrib_ij
        blame_matrix = np.abs(delta_g_err)[:, None] - np.abs(delta_g_err[:, None] + source_contrib)

        harmful_matrix = np.maximum(blame_matrix, 0.0)
        helpful_matrix = np.maximum(-blame_matrix, 0.0)

        response_score = self._normalize_nonnegative(np.abs(delta_g_err))

        err_weights = np.abs(delta_g_err)
        err_weights = err_weights / (err_weights.sum() + eps)

        source_blame_score = harmful_matrix.T @ err_weights
        helpful_source_score = helpful_matrix.T @ err_weights

        source_blame_score = self._normalize_nonnegative(source_blame_score)
        helpful_source_score = self._normalize_nonnegative(helpful_source_score)

        pair_score = harmful_matrix + harmful_matrix.T
        np.fill_diagonal(pair_score, 0.0)

        if np.max(pair_score) > eps:
            pair_score = pair_score / (np.max(pair_score) + eps)

        rows = []
        for idx, (coord, kind) in enumerate(zip(coords_flat, kinds_flat)):
            rows.append({
                "idx": idx,
                "coord": tuple(coord),
                "type": self.coord_type(tuple(coord), kind),
                "dq_raw": dq_raw[idx],
                "dq_eff": dq_eff[idx],
                "abs_dq": abs(dq_eff[idx]),
                "chain": chain[idx],
                "g0": g0[idx],
                "diag_resp": diag_resp[idx],
                "abs_diag_resp": abs(diag_resp[idx]),
                "offdiag_resp": offdiag_resp[idx],
                "abs_offdiag_resp": abs(offdiag_resp[idx]),
                "Hdq": Hdq[idx],
                "abs_Hdq": abs(Hdq[idx]),
                "g_model": g_model[idx],
                "g_qm": g_qm[idx],
                "delta_g_err": delta_g_err[idx],
                "abs_delta_g_err": abs(delta_g_err[idx]),
                "rho": rho[idx],
                "pE_total": partial_energies[idx],
                "abs_pE_total": abs(partial_energies[idx]),
                "source_blame_score": source_blame_score[idx],
                "helpful_source_score": helpful_source_score[idx],
                "response_score": response_score[idx],
            })

        by_abs_dq = sorted(rows, key=lambda r: r["abs_dq"], reverse=True)
        by_abs_Hdq = sorted(rows, key=lambda r: r["abs_Hdq"], reverse=True)
        by_abs_err = sorted(rows, key=lambda r: r["abs_delta_g_err"], reverse=True)
        by_rho = sorted(rows, key=lambda r: r["rho"], reverse=True)
        by_abs_pE = sorted(rows, key=lambda r: r["abs_pE_total"], reverse=True)
        by_source_blame = sorted(rows, key=lambda r: r["source_blame_score"], reverse=True)

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
            "pE_total": partial_energies,
            "source_contrib": source_contrib,
            "blame_matrix": blame_matrix,
            "harmful_matrix": harmful_matrix,
            "helpful_matrix": helpful_matrix,
            "response_score": response_score,
            "source_blame_score": source_blame_score,
            "helpful_source_score": helpful_source_score,
            "pair_score": pair_score,
            "by_abs_dq": by_abs_dq,
            "by_abs_Hdq": by_abs_Hdq,
            "by_abs_err": by_abs_err,
            "by_rho": by_rho,
            "by_abs_pE": by_abs_pE,
            "by_source_blame": by_source_blame,
        }


    def _build_candidate_block(
        self,
        anchor_idx: int,
        global_coord_score: np.ndarray,
        global_pair: np.ndarray,
        z_matrix: List[Tuple[int, ...]],
        coord_kinds: List[str],
        constraints_to_exclude: List[Tuple[int, ...]],
        max_partners: int = 4,
        min_partner_pair_frac: float = 0.15,
        min_partner_coord_frac: float = 0.20,
    ) -> List[int]:
        """
        Build a coupled block around an anchor using blame-derived pair coupling
        plus global coordinate relevance.

        No dihedral bonus is used here.
        """
        N = len(z_matrix)
        anchor_score = float(global_coord_score[anchor_idx])

        pair_row = np.asarray(global_pair[anchor_idx], dtype=float).copy()
        partner_order = np.argsort(-pair_row)

        block = [int(anchor_idx)]

        best_pair = 0.0
        for j in partner_order:
            if j == anchor_idx:
                continue
            best_pair = float(pair_row[j])
            break

        if best_pair <= 0.0:
            return block

        anchor_atoms = set(z_matrix[anchor_idx])

        candidates = []
        for j in range(N):
            if j == anchor_idx:
                continue

            coord_j = tuple(int(x) for x in z_matrix[j])

            if coord_j in constraints_to_exclude:
                continue

            pair_ok = pair_row[j] >= min_partner_pair_frac * best_pair
            coord_ok = global_coord_score[j] >= min_partner_coord_frac * anchor_score

            if not (pair_ok and coord_ok):
                continue

            atoms_j = set(coord_j)
            overlap = len(anchor_atoms.intersection(atoms_j)) / max(1, len(anchor_atoms.union(atoms_j)))

            local_value = 0.65 * float(pair_row[j]) + 0.30 * float(global_coord_score[j]) + 0.05 * float(overlap)
            candidates.append((int(j), local_value))

        candidates = sorted(candidates, key=lambda x: x[1], reverse=True)

        for j, _ in candidates[:max_partners]:
            block.append(int(j))

        partners = [i for i in block if i != anchor_idx]
        partners = sorted(partners, key=lambda i: global_coord_score[i], reverse=True)
        block = [int(anchor_idx)] + partners

        return block


    def _select_constraints_from_block(
        self,
        block_indices: List[int],
        global_coord_score: np.ndarray,
        global_source_blame: np.ndarray,
        global_response_score: np.ndarray,
        global_pair: np.ndarray,
        z_matrix: List[Tuple[int, ...]],
        max_constraints: int = 3,
        singleton_dominance: float = 1.35,
        singleton_source_frac: float = 0.85,
        singleton_response_frac: float = 0.85,
        secondary_score_frac: float = 0.25,
        min_secondary_pair_frac: float = 0.15,
    ):
        """
        Choose final constraints from the winning block.

        Critical behavior:
        ------------------
        If one coordinate clearly dominates as both
        - a source of error and
        - a response channel,
        then return only that coordinate.

        This means:
        if a dihedral is truly the problematic piece, it is returned alone.
        No artificial non-dihedral is added in that case.
        """
        block_indices = list(block_indices)
        if not block_indices:
            return [], []

        eps = 1e-12

        block_scores = {i: float(global_coord_score[i]) for i in block_indices}
        steering_score = {
            i: float(0.5 * global_source_blame[i] + 0.5 * global_response_score[i])
            for i in block_indices
        }

        primary_idx = max(block_indices, key=lambda i: steering_score[i])
        primary_coord = tuple(int(x) for x in z_matrix[primary_idx])
        primary_constraint = [primary_coord]

        selected = [primary_coord]
        remaining = [i for i in block_indices if i != primary_idx]

        second_best = max([block_scores[i] for i in remaining], default=0.0)
        block_max_source = max([float(global_source_blame[i]) for i in block_indices], default=0.0)
        block_max_response = max([float(global_response_score[i]) for i in block_indices], default=0.0)

        # singleton rule:
        # if one coordinate clearly dominates, return only that coordinate
        if (
            block_scores[primary_idx] >= singleton_dominance * max(second_best, eps)
            and float(global_source_blame[primary_idx]) >= singleton_source_frac * max(block_max_source, eps)
            and float(global_response_score[primary_idx]) >= singleton_response_frac * max(block_max_response, eps)
        ):
            return selected, primary_constraint

        pair_row = np.asarray(global_pair[primary_idx], dtype=float)
        best_pair = max([float(pair_row[i]) for i in remaining], default=0.0)

        def secondary_rank_key(i):
            return 0.70 * float(pair_row[i]) + 0.30 * float(global_coord_score[i])

        remaining_sorted = sorted(remaining, key=secondary_rank_key, reverse=True)

        for i in remaining_sorted:
            if len(selected) >= max_constraints:
                break

            if block_scores[i] < secondary_score_frac * max(block_scores[primary_idx], eps):
                continue

            if best_pair > eps:
                pair_ok = float(pair_row[i]) >= min_secondary_pair_frac * best_pair
                if not pair_ok and block_scores[i] < 0.60 * max(block_scores[primary_idx], eps):
                    continue

            coord_i = tuple(int(x) for x in z_matrix[i])
            if coord_i not in selected:
                selected.append(coord_i)

        return selected, primary_constraint
    
    def _section_from_kind(self, kind):
        return {
            "bond": "bonds",
            "angle": "angles",
            "dihedral": "dihedrals",
            "improper": "impropers",
        }.get(kind)


    def _normalized_internal_displacement(self, dq_eff, q_ref, coords_flat, kinds_flat):
        disp = np.zeros(len(coords_flat), dtype=float)

        for idx, kind in enumerate(kinds_flat):
            section = self._section_from_kind(kind)
            if section is None:
                continue

            sigma = self._imp_coordinate_sigma(section, idx, q_ref)
            disp[idx] = abs(float(dq_eff[idx])) / max(float(sigma), 1.0e-12)

        return disp


    def _select_fallback_locality_constraints(
        self,
        *,
        block_indices,
        candidate_constraints,
        primary_constraint,
        coords_flat,
        kinds_flat,
        constraints_to_exclude,
        global_coord_score,
        global_pair,
        global_displacement,
        max_constraints=6,
        min_displacement_sigma=1.0,
        min_score_frac=0.10,
        min_pair_frac=0.10,
    ):
        eps = 1.0e-12
        selected = {tuple(int(x) for x in c)
                    for c in list(candidate_constraints) + list(primary_constraint)}

        selected_atoms = {a for coord in selected for a in coord}
        block_set = set(int(i) for i in block_indices)

        primary_idx = None
        if primary_constraint:
            primary_coord = tuple(int(x) for x in primary_constraint[0])
            for idx, coord in enumerate(coords_flat):
                if tuple(int(x) for x in coord) == primary_coord:
                    primary_idx = idx
                    break

        pair_row = np.zeros(len(coords_flat), dtype=float)
        best_pair = 0.0
        if primary_idx is not None:
            pair_row = np.asarray(global_pair[primary_idx], dtype=float)
            best_pair = max([float(pair_row[i]) for i in block_set if i != primary_idx], default=0.0)

        disp_score = self._normalize_nonnegative(global_displacement)
        block_score_ref = max([float(global_coord_score[i]) for i in block_set], default=float(np.max(global_coord_score)))

        ranked = []
        for idx, coord_raw in enumerate(coords_flat):
            coord = tuple(int(x) for x in coord_raw)
            if coord in selected:
                continue
            if self._is_coord_excluded_for_constraint(idx, coords_flat, kinds_flat, constraints_to_exclude):
                continue

            atoms = set(coord)
            atom_overlap = bool(selected_atoms.intersection(atoms)) if selected_atoms else False
            pair_ok = best_pair > eps and float(pair_row[idx]) >= min_pair_frac * best_pair
            local_ok = idx in block_set or atom_overlap or pair_ok

            disp_ok = float(global_displacement[idx]) >= min_displacement_sigma
            score_ok = float(global_coord_score[idx]) >= min_score_frac * max(block_score_ref, eps)

            if not local_ok or not (disp_ok or score_ok):
                continue

            rank_value = (
                0.50 * float(disp_score[idx])
                + 0.30 * float(global_coord_score[idx])
                + 0.20 * float(pair_row[idx])
            )
            ranked.append((idx, rank_value))

        ranked = sorted(ranked, key=lambda x: x[1], reverse=True)
        return [tuple(int(x) for x in coords_flat[i]) for i, _ in ranked[:max_constraints]]



    def print_ranked_table(
        self,
        rows: List[Dict[str, Any]],
        title: str,
        n: int = 10,
    ) -> None:
        """
        Pretty-print a ranked summary table.
        """
        print(f"\n{title}")
        print("-" * len(title))
        header = (
            f"{'rank':>4} {'idx':>4} {'type':>10} {'coord':>18} "
            f"{'|dq|':>12} {'|diag|':>12} {'|offdiag|':>12} "
            f"{'|Hdq|':>12} {'|err|':>12} {'rho':>12} {'src_blame':>12}"
        )
        print(header)
        print("-" * len(header))

        for rank, r in enumerate(rows[:n]):
            print(
                f"{rank:4d} "
                f"{r['idx']:4d} "
                f"{r['type']:>10} "
                f"{self.format_coord(r['coord']):>18} "
                f"{r['abs_dq']:12.6e} "
                f"{r['abs_diag_resp']:12.6e} "
                f"{r['abs_offdiag_resp']:12.6e} "
                f"{r['abs_Hdq']:12.6e} "
                f"{r['abs_delta_g_err']:12.6e} "
                f"{r['rho']:12.6e} "
                f"{r['source_blame_score']:12.6e}"
            )


    def print_type_summary(self, rows: List[Dict[str, Any]]) -> None:
        """
        Aggregate diagnostics by coordinate type.
        """
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
                    "sum_src_blame": 0.0,
                }
            buckets[t]["count"] += 1
            buckets[t]["sum_abs_dq"] += r["abs_dq"]
            buckets[t]["sum_abs_Hdq"] += r["abs_Hdq"]
            buckets[t]["sum_abs_err"] += r["abs_delta_g_err"]
            buckets[t]["sum_abs_pE"] += r["abs_pE_total"]
            buckets[t]["sum_src_blame"] += r["source_blame_score"]

        print("\nSummary by coordinate type")
        print("--------------------------")
        print(
            f"{'type':>10} {'count':>7} {'sum|dq|':>14} "
            f"{'sum|Hdq|':>14} {'sum|err|':>14} {'sum|pE|':>14} {'sum src':>14}"
        )
        print("-" * 94)
        for t, vals in buckets.items():
            print(
                f"{t:>10} "
                f"{vals['count']:7d} "
                f"{vals['sum_abs_dq']:14.6e} "
                f"{vals['sum_abs_Hdq']:14.6e} "
                f"{vals['sum_abs_err']:14.6e} "
                f"{vals['sum_abs_pE']:14.6e} "
                f"{vals['sum_src_blame']:14.6e}"
            )
   
    def transform_gradient_to_internal_coordinates(self, molecule, gradient, b_matrix, tol=1e-6):
        grad = np.asarray(gradient, dtype=np.float64).reshape(-1)
        B = np.asarray(b_matrix, dtype=np.float64)

        natm = grad.size // 3
        target_rank = 3 * natm - 6
        if natm == 2:
            target_rank += 1
        target_rank = max(1, int(target_rank))

        G = B @ B.T
        U, s, Vt = np.linalg.svd(G, full_matrices=False)

        s_inv = np.zeros_like(s)
        nz = np.where(s > tol)[0]

        if nz.size >= target_rank:
            keep = nz[np.argsort(s[nz])[-target_rank:]]
        else:
            keep = np.argsort(s)[-target_rank:]

        s_inv[keep] = 1.0 / s[keep]
        G_pinv = (U * s_inv) @ Vt
        return G_pinv @ (B @ grad)

    def get_energy(self):
        """ Returns the potential energy obtained by interpolation.
        """
        return self.impes_coordinate.energy

    def get_gradient(self):
        """ Returns the gradient obtained by interpolation.
        """
        return self.impes_coordinate.gradient
    

    # def perform_symmetry_assignment(self, atom_map, sym_group, reference_group, datapoint_group):
    #     """ Performs the atom mapping. """

    #     new_map = atom_map.copy()
    #     mapping_dict = {}
    #     # cost = self.get_dihedral_cost(atom_map, sym_group, non_group_atoms)
    #     cost = np.linalg.norm(datapoint_group[:, np.newaxis, :] - reference_group[np.newaxis, :, :], axis=2)
    #     row, col = linear_sum_assignment(cost)
    #     assigned = False
    #     print('row col', row, col)
    #     if not np.equal(row, col).all():
    #         assigned = True
            
    #         # atom_maps = self.linear_assignment_solver(cost)
    #         reordred_arr = np.array(sym_group)[col]
    #         new_map[sym_group] = new_map[reordred_arr]

    #         for org, new in zip(np.array(sym_group), reordred_arr):
    #             print('org, neww', org, new)
    #         mapping_dict = {org: new for org, new in zip(np.array(sym_group), reordred_arr)}
        
    #     return [new_map], mapping_dict, assigned

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
    
