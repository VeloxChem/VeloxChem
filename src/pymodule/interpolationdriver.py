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
from scipy.optimize import linear_sum_assignment
import sys
import os
from .profiler import Profiler
import h5py
import math
import copy
from contextlib import redirect_stderr
from io import StringIO
from .interpolationdatapoint import InterpolationDatapoint

from .outputstream import OutputStream
from .veloxchemlib import mpi_master

from .errorhandler import assert_msg_critical
from .inputparser import (parse_input)

from .imrotorbuilder import build_runtime_cluster_set


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
        
        self.z_matrix = z_matrix
        self.impes_dict = None
        self.sum_of_weights = None
        self.sum_of_weights_grad = None
        self.print = False
        self.store_weights = False
        self.perform_mapping = True
        # symmetry metadata tuple produced during database-point generation
        self.symmetry_information = None
        # global switch controlling symmetry-expanded candidate evaluation
        self.use_symmetry = True
        self.weights = {}
        self.int_coord_weights = {}
        self.potentials = []
        self.gradients = []

        # kpoint file with QM data
        self.imforcefield_file = None
        self.qm_data_points = None
        # per-label torsion descriptors used by symmetry-aware distance metrics
        self.symmetry_dihedral_lists = {}
        # grouped datapoints: core label -> [core point, symmetry-expanded points...]
        self.qm_symmetry_data_points = {}
        # legacy containers kept for backward-compatible symmetry workflows
        self.qm_symmetry_data_points_1 = None
        self.qm_symmetry_data_points_2 = None

        # active rotor-cluster metadata for the current datapoint family
        self.rotor_cluster_information = None
        # cached runtime masks/projectors derived from rotor_cluster_information
        self.rotor_cluster_runtime = None
        # family_label -> {core point, cluster registry, expected states}
        self.qm_rotor_cluster_banks = {}

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
        self.eq_bond_symmetry_mode = "masked_exact"  # "masked_exact" | "symmetrized"
        self.use_cosine_dihedral = False
        self.use_tc_weights = True

        # Optional runtime profiling for interpolation bottleneck analysis.
        self.runtime_profile_enabled = False
        self.runtime_profile_print = False
        self.runtime_profile_profiler = None
        self.runtime_profile_totals = {}
        self.runtime_profile_last = {}
        self.runtime_profile_calls = 0

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
                    'use_eq_bond_length': ('bool', 'whether to use eq bond lengths in the Z-matrix'),
                    'eq_bond_symmetry_mode': ('str', 'eq-bond symmetry mode: masked_exact or symmetrized'),
                    'use_cosine_dihedral':('bool', 'wether to use cosine and sin for the diehdral in the Z-matrix'),
                    'use_tc_weights':('bool', 'weither to use target coustomized weights'),
                    # 'use_mpi_preload': ('bool', 'enable MPI preload backend for compute_potential'),
                    'use_symmetry': ('bool', 'enable symmetry-expanded candidate evaluation'),
                'labels': ('seq_fixed_str', 'the list of QM data point labels'),
            }
        }


    def _get_eq_symmetry_mode(self):
        mode = getattr(self.impes_coordinate, 'eq_bond_symmetry_mode', self.eq_bond_symmetry_mode)
        if mode is None:
            mode = self.eq_bond_symmetry_mode
        mode = str(mode).strip().lower()
        if mode not in ('masked_exact', 'symmetrized'):
            raise ValueError(
                f"InterpolationDriver: invalid eq_bond_symmetry_mode='{mode}'. "
                "Use 'masked_exact' or 'symmetrized'."
            )
        return mode

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

    def _cluster_weight_metric_worker(self, payload):
        cluster_id, datapoint, org_int_coords, base_cache = payload
        d2_i, grad_d2_i = self._cluster_state_signature_metric(
            cluster_id, datapoint, org_int_coords, base_cache=base_cache
        )
        return float(d2_i), np.asarray(grad_d2_i, dtype=np.float64)

    def mark_runtime_data_cache_dirty(self):
        """
        Marks the flattened symmetry task cache as outdated.
        """

        self._runtime_data_cache_dirty = True


    def prepare_runtime_data_cache(self, force=False):
        """
        Pre-builds flattened symmetry task entries from symmetry datapoints.

        Each entry is a tuple ``(symmetry_data_point, mask0, mask)`` used by
        potential evaluation. The expensive object/mask traversal is done once
        and reused across timesteps.
        """
       
        if not self.runtime_data_cache_enabled:
            self._symmetry_task_cache = {}
            self._runtime_data_cache_dirty = False
            return

        if (not force) and (not self._runtime_data_cache_dirty):
            return

        cache = {}

        for dp_label, symmetry_data_points in (self.qm_symmetry_data_points or {}).items():
            entries = []
            for symmetry_data_point in symmetry_data_points:
                masks = getattr(symmetry_data_point, 'mapping_masks', None)
                if masks is None:
                    continue

                masks_arr = np.asarray(masks, dtype=np.int64)
                if masks_arr.ndim == 1:
                    masks_arr = masks_arr.reshape(1, -1)
                if masks_arr.shape[0] == 0:
                    continue

                mask0 = np.ascontiguousarray(masks_arr[0], dtype=np.int64)
                for mask in masks_arr:
                    entries.append(
                        (symmetry_data_point, mask0,
                         np.ascontiguousarray(mask, dtype=np.int64)))
            cache[dp_label] = entries

        self._symmetry_task_cache = cache
        self._runtime_data_cache_dirty = False
  
    def _use_symmetry_for_label(self, dp_label):
        sym_list = (self.qm_symmetry_data_points or {}).get(dp_label, [])
        return bool(len(sym_list) > 1)
    
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
        if self._use_symmetry_for_label(dp_label):
            return self._get_symmetry_task_entries(dp_label)
        dp = self._get_datapoint_by_label(dp_label)   # from qm_data_points
        n_ic = len(dp.internal_coordinates_values)
        ident = np.arange(n_ic, dtype=np.int64)
        return [(dp, ident, ident)]

    def _get_symmetry_task_entries(self, dp_label):
        """
        Returns flattened symmetry task entries for one datapoint label.
        """

        if self.runtime_data_cache_enabled:
            self.prepare_runtime_data_cache()
            cached_entries = self._symmetry_task_cache.get(dp_label)
            if cached_entries is not None:
                return cached_entries

        entries = []
        symmetry_data_points = (self.qm_symmetry_data_points or {}).get(dp_label, [])

        for symmetry_data_point in symmetry_data_points:
            masks = getattr(symmetry_data_point, 'mapping_masks', None)
            if masks is None:
                continue

            masks_arr = np.asarray(masks, dtype=np.int64)
            if masks_arr.ndim == 1:
                masks_arr = masks_arr.reshape(1, -1)
            if masks_arr.shape[0] == 0:
                continue

            mask0 = np.ascontiguousarray(masks_arr[0], dtype=np.int64)
            for mask in masks_arr:
                entries.append(
                    (symmetry_data_point, mask0,
                     np.ascontiguousarray(mask, dtype=np.int64)))
  
        return entries

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

    def _ensure_rotor_cluster_runtime(self):
        if self.rotor_cluster_information is None:
            self.rotor_cluster_runtime = None
            return
        n_ic = len(self.impes_coordinate.internal_coordinates_values)

       
        if self.rotor_cluster_runtime is None or self.rotor_cluster_runtime.n_ic != n_ic:
            self.rotor_cluster_runtime = build_runtime_cluster_set(
                self.rotor_cluster_information, n_ic=n_ic
            )


    def _get_core_projector(self):
        self._ensure_rotor_cluster_runtime()

        return self.rotor_cluster_runtime.core_mask


    def _get_cluster_projector(self, cluster_id):
        self._ensure_rotor_cluster_runtime()
        return self.rotor_cluster_runtime.clusters[cluster_id].projector_mask


    def _get_cluster_signature_rows(self, cluster_id):
        self._ensure_rotor_cluster_runtime()
        return self.rotor_cluster_runtime.clusters[cluster_id].grouped_signature_rows
    
    def _build_cluster_signature(self, cluster_id, int_coords):
        rows_by_rotor = self._get_cluster_signature_rows(cluster_id)
        grouped = []
        for rows in rows_by_rotor:
            grouped.append(np.asarray(int_coords[list(rows)], dtype=np.float64).copy())
        return grouped

    def _cluster_state_signature_metric(self, cluster_id, datapoint, org_int_coords, base_cache=None):
        masks = getattr(datapoint, "mapping_masks", None)

        if masks is None or len(masks) == 0:
            cur = self._build_cluster_signature(cluster_id, np.asarray(org_int_coords, dtype=np.float64))
            ref = self._build_cluster_signature(cluster_id, np.asarray(datapoint.internal_coordinates_values, dtype=np.float64))
            return self._cluster_signature_distance_and_gradient(cluster_id, cur, ref)

        best_d2 = None
        best_grad = None
        for mask in np.asarray(masks, dtype=np.int64):
            d2, gd2 = self._cluster_mask_signature_distance_and_gradient(
                cluster_id=cluster_id,
                datapoint=datapoint,
                mask=mask,
                org_int_coords=org_int_coords,
                base_cache=base_cache,
            )
            if best_d2 is None or d2 < best_d2:
                best_d2, best_grad = d2, gd2
        return best_d2, best_grad

    def _cluster_signature_distance_and_gradient(self, cluster_id, current_signature, reference_signature, b_matrix=None):
        rows_by_rotor = self._get_cluster_signature_rows(cluster_id)
        if b_matrix is None:
            b_matrix = self.impes_coordinate.b_matrix
        b_matrix = np.asarray(b_matrix, dtype=np.float64)

        ncart = b_matrix.shape[1]
        d2_total = 0.0
        grad_total = np.zeros(ncart, dtype=np.float64)
        n_groups = len(rows_by_rotor)

        for ridx, rows in enumerate(rows_by_rotor):
            rows_arr = np.asarray(rows, dtype=np.int64)
            brows = b_matrix[rows_arr, :]
            sig_c = np.asarray(current_signature[ridx], dtype=np.float64)
            sig_r = np.asarray(reference_signature[ridx], dtype=np.float64)
            delta = sig_c - sig_r
            terms = 2.0 * (1.0 - np.cos(delta))
            d2_r = float(np.mean(terms))
            grad_r = np.mean(2.0 * np.sin(delta)[:, None] * brows, axis=0)
            d2_total += d2_r
            grad_total += grad_r

        if n_groups > 0:
            d2_total /= n_groups
            grad_total /= n_groups

        return d2_total, grad_total
    
    def _evaluate_masked_cluster_state(self, cluster_id, datapoint, projector, org_int_coords, base_cache=None):
        masks = datapoint.mapping_masks
        if masks is None or len(masks) == 0:
            masks = [np.arange(len(datapoint.internal_coordinates_values), dtype=np.int64)]

        masked_E = []
        masked_G = []
        raw = []
        raw_grad = []

        for mask in masks:
            mask_arr = np.asarray(mask, dtype=np.int64)

            org_eval_m, b_eval_m = self._candidate_current_chart(
                symmetry_data_point=datapoint,
                mask=mask_arr,
                org_int_coords=org_int_coords,
                b_matrix=self.impes_coordinate.b_matrix,
                base_cache=base_cache,
            )

            E_m, G_m = self._restricted_taylor_eval_for_mask(
                datapoint=datapoint,
                mask=mask_arr,
                projector=projector,
                org_int_coords=org_int_coords,
                base_cache=base_cache,
                precomputed_current=(org_eval_m, b_eval_m),
            )

            d2_m, grad_d2_m = self._cluster_mask_signature_distance_and_gradient(
                cluster_id=cluster_id,
                datapoint=datapoint,
                mask=mask_arr,
                org_int_coords=org_int_coords,
                b_matrix=b_eval_m,
                base_cache=base_cache,
                precomputed_current=(org_eval_m, b_eval_m),
            )

            w_m, grad_w_m = self._normalize_symmetry_candidate_weight(d2_m, grad_d2_m)

            masked_E.append(E_m)
            masked_G.append(G_m)
            raw.append(w_m)
            raw_grad.append(grad_w_m)

        weights, grad_weights = self._normalize_mask_weight_pool(raw, raw_grad)

        state_E = np.dot(weights, np.asarray(masked_E, dtype=np.float64))
        state_G = (
            np.tensordot(weights, np.asarray(masked_G, dtype=np.float64), axes=1)
            + np.tensordot(np.asarray(masked_E, dtype=np.float64), grad_weights, axes=1)
        )
        return state_E, state_G
    
    def _cluster_mask_signature_distance_and_gradient(
        self,
        *,
        cluster_id,
        datapoint,
        mask,
        org_int_coords,
        b_matrix=None,
        base_cache=None,
        precomputed_current=None,
    ):
        mask = np.asarray(mask, dtype=np.int64)

        if precomputed_current is None:
            org_eval, b_eval = self._candidate_current_chart(
                symmetry_data_point=datapoint,
                mask=mask,
                org_int_coords=org_int_coords,
                b_matrix=self.impes_coordinate.b_matrix if b_matrix is None else b_matrix,
                base_cache=base_cache,
            )
        else:
            org_eval, b_eval = precomputed_current

        cur_sig = self._build_cluster_signature(cluster_id, np.asarray(org_eval, dtype=np.float64))
        ref_masked = np.asarray(datapoint.internal_coordinates_values, dtype=np.float64)[mask]
        ref_sig = self._build_cluster_signature(cluster_id, ref_masked)

        return self._cluster_signature_distance_and_gradient(
            cluster_id, cur_sig, ref_sig, b_matrix=b_eval
        )

    def _normalize_symmetry_candidate_weight(
        self,
        d2,
        grad_d2,
        *,
        confidence_radius=1.0,
        power=2,
        eps=1.0e-12,
    ):
        """
        Raw inverse-distance-like weight for one symmetry-mask candidate.
        Returns scalar weight and Cartesian gradient (natm,3).
        """
        rho2 = max(float(confidence_radius) ** 2, eps)
        x = float(d2) / rho2 + eps
        w = 1.0 / (2.0 * (x ** power))
        dw_dd2 = -(power / (2.0 * rho2)) * (x ** (-(power + 1)))
        grad_w = dw_dd2 * np.asarray(grad_d2, dtype=np.float64)

        natm = self.impes_coordinate.cartesian_coordinates.shape[0]
        return float(w), grad_w.reshape(natm, 3)


    def _normalize_mask_weight_pool(self, raw, raw_grad, *, eps=1.0e-14):
        """
        Normalize a pool of raw mask-local weights and their gradients.
        """
        raw = np.asarray(raw, dtype=np.float64)
        raw_grad = np.asarray(raw_grad, dtype=np.float64)

        natm = self.impes_coordinate.cartesian_coordinates.shape[0]
        if raw.size == 0:
            return (
                np.zeros((0,), dtype=np.float64),
                np.zeros((0, natm, 3), dtype=np.float64),
            )

        S = raw.sum()
        if S < eps:
            idx = int(np.argmax(raw))
            weights = np.zeros_like(raw)
            weights[idx] = 1.0
            grad_weights = np.zeros_like(raw_grad)
            return weights, grad_weights

        sum_grad = raw_grad.sum(axis=0)
        weights = raw / S
        grad_weights = (raw_grad * S - raw[:, None, None] * sum_grad) / (S ** 2)
        return weights, grad_weights

    def _restricted_taylor_eval_for_mask(
        self,
        *,
        datapoint,
        mask,
        projector,
        org_int_coords,
        base_cache=None,
        precomputed_current=None,
    ):
        mask = np.asarray(mask, dtype=np.int64)
        n_ic = len(datapoint.internal_coordinates_values)
        mask0 = np.arange(n_ic, dtype=np.int64)

        grad_eff = np.asarray(datapoint.internal_gradient, dtype=np.float64).copy()
        hess_eff = np.asarray(datapoint.internal_hessian, dtype=np.float64).copy()
        ref_masked = np.asarray(datapoint.internal_coordinates_values, dtype=np.float64)[mask]

        grad_eff[mask0] = grad_eff[mask]
        hess_eff[np.ix_(mask0, mask0)] = hess_eff[np.ix_(mask, mask)]

        bounds = self._get_internal_coordinate_partitions()
        dstart = int(bounds["dihedral_start"])
        dend = int(bounds["dihedral_end"])
        bend = int(bounds['bond_end'])
        aend = int(bounds['angle_end'])

        if precomputed_current is None:
            org_eval, b_eval = self._candidate_current_chart(
                symmetry_data_point=datapoint,
                mask=mask,
                org_int_coords=np.asarray(org_int_coords, dtype=np.float64),
                b_matrix=np.asarray(self.impes_coordinate.b_matrix, dtype=np.float64),
                base_cache=base_cache,
            )
        else:
            org_eval, b_eval = precomputed_current

        return self._restricted_taylor_eval(
            energy=float(datapoint.energy),
            grad=grad_eff,
            hess=hess_eff,
            ref_coords=ref_masked,
            org_int_coords=org_eval,
            b_matrix=b_eval,
            projector=np.asarray(projector, dtype=np.float64),
            dihedral_start=dstart,
            dihedral_end=dend,
            angle_end=aend,
            bond_end=bend,
        )

    def _restricted_taylor_eval(
        self,
        *,
        energy,
        grad,
        hess,
        ref_coords,
        org_int_coords,
        b_matrix,
        projector,
        dihedral_start,
        dihedral_end,
        angle_end,
        bond_end
    ):
        if self.bond_rmsd is None:
            self.bond_rmsd = []
        if self.angle_rmsd is None:
            self.angle_rmsd = []
        if self.dihedral_rmsd is None:
            self.dihedral_rmsd = []
        delta_raw = np.asarray(org_int_coords, dtype=np.float64) - np.asarray(ref_coords, dtype=np.float64)
        delta_eff = delta_raw.copy()
        chain = np.ones_like(delta_raw)
        delta_eff[dihedral_start:dihedral_end] = np.sin(delta_raw[dihedral_start:dihedral_end])
        chain[dihedral_start:dihedral_end] = np.cos(delta_raw[dihedral_start:dihedral_end])
        
        self.bond_rmsd.append(
            np.sqrt(np.mean(np.sum((delta_eff[:bond_end])**2)))
        )
        self.angle_rmsd.append(
            np.sqrt(np.mean(np.sum(delta_eff[bond_end:angle_end]**2)))
        )
        self.dihedral_rmsd.append(
            np.sqrt(np.mean(np.sum(delta_eff[dihedral_start:dihedral_end]**2)))
        )

        delta_cluster = projector * delta_eff
        
        U = float(energy + grad @ delta_cluster + 0.5 * delta_cluster @ hess @ delta_cluster)
        grad_eff_cluster = projector * (grad + hess @ delta_cluster)
        grad_internal = chain * grad_eff_cluster
        natm = b_matrix.shape[1] // 3
        grad_cart = (b_matrix.T @ grad_internal).reshape(natm, 3)
        return U, grad_cart

    def _assemble_clustered_family_model(self, family_label, org_int_coords):
        family_bank = self.qm_rotor_cluster_banks[family_label]
        family_cluster_info = family_bank.get("cluster_info")

        if family_cluster_info is not None:
            if self.rotor_cluster_information is not family_cluster_info:
                self.rotor_cluster_information = family_cluster_info
                self.rotor_cluster_runtime = None

        core_dp = family_bank["core"]
        bmat = self.impes_coordinate.b_matrix
        bounds = self._get_internal_coordinate_partitions()
        bend = int(bounds['bond_end'])
        aend = int(bounds['angle_end'])
        dstart = int(bounds["dihedral_start"])
        dend = int(bounds["dihedral_end"])

        eq_mode = None
        candidate_base_cache = None
        if self.impes_coordinate.use_eq_bond_length:
            eq_mode = self._get_eq_symmetry_mode()
            if eq_mode == 'masked_exact':
                candidate_base_cache = self.impes_coordinate.prepare_eq_candidate_base_cache()
 
        P0 = self._get_core_projector()

        core_E, core_G = self._restricted_taylor_eval(
            energy=float(core_dp.energy),
            grad=np.asarray(core_dp.internal_gradient, dtype=np.float64),
            hess=np.asarray(core_dp.internal_hessian, dtype=np.float64),
            ref_coords=np.asarray(core_dp.internal_coordinates_values, dtype=np.float64),
            org_int_coords=np.asarray(org_int_coords, dtype=np.float64),
            b_matrix=bmat,
            projector=P0,
            dihedral_start=dstart,
            dihedral_end=dend,
            angle_end=aend,
            bond_end=bend,
        )

        local_E = core_E
        local_G = core_G.copy()
        org_coords = np.asarray(org_int_coords, dtype=np.float64)

        for cid, cluster_bank in family_bank["clusters"].items():
            projector = self._get_cluster_projector(cid)
  
            weights, grad_weights, states = self._compute_cluster_weights(
                cid,
                cluster_bank,
                org_coords,
                base_cache=candidate_base_cache,
                # executor=executor,
                eq_mode=eq_mode,
            )

            if len(states) == 0:
                continue

            dE_list = []
            dG_list = []

            t0 = time()
            for _, dp in states:
                E_state, G_state = self._evaluate_masked_cluster_state(
                    cluster_id=cid,
                    datapoint=dp,
                    projector=projector,
                    org_int_coords=org_coords,
                    base_cache=candidate_base_cache,
                )
                dE_list.append(E_state - core_E)
                dG_list.append(G_state - core_G)
            self._add_runtime_timing('compute.cluster.state_eval_serial', time() - t0)

            dE = np.asarray(dE_list, dtype=np.float64)
            dG = np.asarray(dG_list, dtype=np.float64)

            local_E += np.dot(weights, dE)
            local_G += (
                np.tensordot(weights, dG, axes=1)
                + np.tensordot(dE, grad_weights, axes=1)
            )

        return local_E, local_G

    
    def _compute_cluster_weights(
        self,
        cluster_id,
        cluster_bank,
        org_int_coords,
        base_cache=None,
        # executor=None,
        eq_mode=None,
    ):
        populated = [
            (sid, dp) for sid, dp in sorted(cluster_bank["expected_states"].items())
            if dp is not None
        ]

        natm = self.impes_coordinate.cartesian_coordinates.shape[0]
        if len(populated) == 0:
            return (
                np.zeros((0,), dtype=np.float64),
                np.zeros((0, natm, 3), dtype=np.float64),
                [],
            )
        t0 = time()
        metric_results = [
            self._cluster_weight_metric_worker((cluster_id, dp, org_int_coords, base_cache))
            for _, dp in populated
        ]
        self._add_runtime_timing('compute.cluster.weight_eval_serial', time() - t0)

        raw = []
        raw_grad = []

        # for (_, _dp), (d2_i, grad_d2_i) in zip(populated, metric_results):
        #     rho2 = 1.0
        #     x = d2_i / rho2 + 1.0e-12
        #     w_i = 1.0 / (2.0 * x**2)
        #     dw_dd2 = -(2.0 / (2.0 * rho2)) * x**(-3)
        #     grad_w_i = dw_dd2 * grad_d2_i

        #     raw.append(w_i)
        #     raw_grad.append(grad_w_i.reshape(natm, 3))

        for (_, dp), (d2_i, grad_d2_i) in zip(populated, metric_results):
            R = getattr(dp, "confidence_radius", 1.0)
            if R is None:
                R = 1.0
            R = float(np.asarray(R, dtype=np.float64).reshape(-1)[0])
            R = max(R, 1.0e-8)

            rho2 = R * R
            u = max(float(d2_i) / rho2, 1.0e-12)

            p = float(self.exponent_p if self.exponent_p is not None else 2.0)
            q = float(self.exponent_q if self.exponent_q is not None else p / 2.0)

            denom = u**p + u**q
            w_i = 1.0 / denom

            ddenom_du = p * u**(p - 1.0) + q * u**(q - 1.0)
            du_dd2 = 1.0 / rho2
            dw_dd2 = -(ddenom_du * du_dd2) / (denom * denom)

            grad_w_i = dw_dd2 * grad_d2_i

            raw.append(w_i)
            raw_grad.append(grad_w_i.reshape(natm, 3))

        raw = np.asarray(raw, dtype=np.float64)
        raw_grad = np.asarray(raw_grad, dtype=np.float64)
        
        S = raw.sum()
        if S < 1.0e-14:
            idx = int(np.argmax(raw))
            weights = np.zeros_like(raw)
            weights[idx] = 1.0
            grad_weights = np.zeros_like(raw_grad)
            return weights, grad_weights, populated

        sum_grad = raw_grad.sum(axis=0)
        weights = raw / S

        grad_weights = (raw_grad * S - raw[:, None, None] * sum_grad) / (S**2)
        return weights, grad_weights, populated
    
    def _compute_potential_clustered(self, data_point, org_int_coords):
        family_label = data_point.family_label or data_point.point_label
        family_bank = self.qm_rotor_cluster_banks.get(family_label)

        if not isinstance(family_bank, dict):
            return None
        if family_bank.get("cluster_info") is None:
            return None
        if family_bank.get("core") is None:
            return None
        if not isinstance(family_bank.get("clusters"), dict):
            return None
        pes, gradient = self._assemble_clustered_family_model(family_label, org_int_coords)

        return pes, gradient, 0.0
    
    ###### End of the new rotor set up #####################
    
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

        # self.distance_molecule = Molecule.from_xyz_string(xyz_string)

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
            self.prepare_runtime_data_cache()
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
            for i, data_point in enumerate(self.qm_data_points[:]):

                distance, dihedral_dist, denominator, weight_gradient, distance_vec, _, _, dw_da_dx = self.cartesian_distance(data_point)
                
                
                if abs(distance) < min_distance:
                    min_distance = abs(distance)
                distances_and_gradients.append((distance, dihedral_dist, i, denominator, weight_gradient, distance_vec))

        elif self.weightfunction_type == 'internal':
            for i, data_point in enumerate(self.qm_data_points[:]):
                distance, dihedral_dist, denominator, weight_gradient, distance_vec, _, _, dw_da_dx = self.internal_distance(data_point)

                if abs(distance) < min_distance:
                    min_distance = abs(distance)
                distances_and_gradients.append((distance, dihedral_dist, i, denominator, weight_gradient, distance_vec))

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


    def _get_symmetry_torsion_meta(self):
        """
        Prepares row-index metadata for vectorized torsion handling.

        :returns:
            Dict containing dihedral start index, symmetry-center mapping and
            row indices for periodicity-3 and periodicity-2 contributions.
        """

        bounds = self._get_internal_coordinate_partitions()
        dihedral_start = bounds['dihedral_start']

        if self.symmetry_information is None or len(self.symmetry_information) <= 7:
            return {
                'sym3_keys': (),
                'sym3_rows': np.array([], dtype=np.int64),
                'sym3_center_ids': np.array([], dtype=np.int64),
                'sym2_rows': np.array([], dtype=np.int64),
                'dihedral_start': int(bounds['dihedral_start']),
                'dihedral_end': int(bounds['dihedral_end']),
            }

        sym3_keys = tuple(self.symmetry_information[7][3].keys())
        sym3_key_to_idx = {key: idx for idx, key in enumerate(sym3_keys)}
        sym2_set = {tuple(sorted(entry)) for entry in self.symmetry_information[7][2]}

        sym3_rows = []
        sym3_center_ids = []
        sym2_rows = []

        for row_idx, element in enumerate(self.impes_coordinate.z_matrix_dict['dihedrals']):
            center = tuple(element[1:3])
            if center in sym3_key_to_idx:
                sym3_rows.append(dihedral_start + row_idx)
                sym3_center_ids.append(sym3_key_to_idx[center])
            elif tuple(sorted(element)) in sym2_set:
                sym2_rows.append(row_idx)

        return {
            'sym3_keys': sym3_keys,
            'sym3_rows': np.asarray(sym3_rows, dtype=np.int64),
            'sym3_center_ids': np.asarray(sym3_center_ids, dtype=np.int64),
            'sym2_rows': np.asarray(sym2_rows, dtype=np.int64),
            'dihedral_start': int(bounds['dihedral_start']),
            'dihedral_end': int(bounds['dihedral_end']),
        }
    
    def _build_symmetry_phase_candidates(self, dp_label, torsion_meta):
        """
        Build cheap phase-labeled symmetry candidates.

        Returns
        -------
        candidates : list of dict
            Each dict contains:
            - symmetry_data_point
            - mask0
            - mask
            - candidate_full_phases
        """
        symmetry_task_entries = self._get_symmetry_task_entries(dp_label)
        candidates = []

        for symmetry_data_point, mask0, mask in symmetry_task_entries:
            candidate_full_phases = self._compute_collective_methyl_phases_full(
                symmetry_data_point.internal_coordinates_values[mask],
                torsion_meta
            )

            candidates.append({
                "symmetry_data_point": symmetry_data_point,
                "mask0": mask0,
                "mask": mask,
                "candidate_full_phases": candidate_full_phases,
            })

        return candidates

    def _vectorized_symmetry_dihedral_terms(self, dist_org, dist_correlation, b_matrix, torsion_meta):
        """
        Computes symmetry dihedral weights and derivatives in vectorized form.

        :param dist_org:
            Internal-coordinate displacement before dihedral sin transform.
        :param dist_correlation:
            Correlation displacement array updated in-place for selected rows.
        :param b_matrix:
            Current B matrix.
        :param torsion_meta:
            Metadata prepared by ``_get_symmetry_torsion_meta``.

        :returns:
            Tuple of scalar/array accumulators required in compute_potential.
        """

        # phase construction section
        ncart = b_matrix.shape[1]
        n_sym3 = len(torsion_meta['sym3_keys'])
        sym3_rows = torsion_meta['sym3_rows']
        center_ids = torsion_meta['sym3_center_ids']
        
        phases = np.zeros(n_sym3, dtype=np.float64)
        if sym3_rows.size == 0:
         return phases

        theta = dist_org[sym3_rows]  # radians

        # Fold threefold symmetry into the phase average
        exp_3theta = np.exp(1j * 3.0 * theta)

        z_real = np.bincount(center_ids, weights=exp_3theta.real, minlength=n_sym3)
        z_imag = np.bincount(center_ids, weights=exp_3theta.imag, minlength=n_sym3)
        counts = np.bincount(center_ids, minlength=n_sym3)

        z_real /= np.maximum(counts, 1)
        z_imag /= np.maximum(counts, 1)

        phases_3 = np.arctan2(z_imag, z_real)
        phases = phases_3 / 3.0

        ncart = b_matrix.shape[1]

        n_sym3 = len(torsion_meta['sym3_keys'])
        sum_sym_3_dihedral_cos = np.zeros(n_sym3, dtype=np.float64)
        sum_sym_3_dihedral_prime_cos = np.zeros((n_sym3, ncart), dtype=np.float64)

        sym3_rows = torsion_meta['sym3_rows']
        if sym3_rows.size > 0:
            theta3 = dist_org[sym3_rows]

            dist_correlation[sym3_rows] = 0.0

            cos_terms = 2.0 * (1.0 - np.cos(theta3))
            sin_terms = 2.0 * np.sin(theta3)
            center_ids = torsion_meta['sym3_center_ids']

            sum_sym_3_dihedral_cos += np.bincount(
                center_ids, weights=cos_terms, minlength=n_sym3)
                

            sin_b = sin_terms[:, None] * b_matrix[sym3_rows, :]
            np.add.at(sum_sym_3_dihedral_prime_cos, center_ids, sin_b)

        return (
                sum_sym_3_dihedral_cos,
                sum_sym_3_dihedral_prime_cos)

    def _compute_collective_methyl_phases_sym3(self, int_coords, torsion_meta, s=0):
        """
        Compute one threefold-symmetry-adapted phase per methyl center.
        int_coords must contain absolute torsion values in radians.
        """
        sym3_rows = torsion_meta['sym3_rows']
        center_ids = torsion_meta['sym3_center_ids']
        n_sym3 = len(torsion_meta['sym3_keys'])

        phases = np.zeros(n_sym3, dtype=np.float64)

        if sym3_rows.size == 0:
            return phases

        theta = int_coords[sym3_rows]

        exp_3theta = np.exp(1j * theta)

        z_real = np.bincount(center_ids, weights=exp_3theta.real, minlength=n_sym3)
        z_imag = np.bincount(center_ids, weights=exp_3theta.imag, minlength=n_sym3)
        counts = np.bincount(center_ids, minlength=n_sym3)

        z_real /= np.maximum(counts, 1)
        z_imag /= np.maximum(counts, 1)

        phases_3 = np.arctan2(z_imag, z_real)
        phases = phases_3 + (2.0 * np.pi / 3.0) * s
        return phases

    def _build_grouped_symmetry_signatures(self, int_coords, torsion_meta):
        """
        Build one torsion signature block per methyl rotor.

        Parameters
        ----------
        int_coords : np.ndarray
            Absolute internal coordinates in radians.
        torsion_meta : dict
            Must contain:
            - 'sym3_rows'
            - 'sym3_center_ids'
            - 'sym3_keys'

        Returns
        -------
        grouped : list[np.ndarray]
            grouped[r] is the torsion signature block for methyl rotor r.
        """
        sym3_rows = torsion_meta['sym3_rows']
        center_ids = torsion_meta['sym3_center_ids']
        n_sym3 = len(torsion_meta['sym3_keys'])

        grouped = []
        for rid in range(n_sym3):
            rows_r = sym3_rows[center_ids == rid]
            grouped.append(np.asarray(int_coords[rows_r], dtype=np.float64).copy())

        return grouped
    
    def _build_grouped_symmetry_brows(self, torsion_meta):
        sym3_rows = torsion_meta['sym3_rows']
        center_ids = torsion_meta['sym3_center_ids']
        n_sym3 = len(torsion_meta['sym3_keys'])

        grouped = []
        for rid in range(n_sym3):
            rows_r = sym3_rows[center_ids == rid]
            grouped.append(self.impes_coordinate.b_matrix[rows_r, :].copy())

        return grouped
    
    def _grouped_signature_distance_and_gradient(
        self,
        grouped_current,
        grouped_candidate,
        grouped_brows,
        rotor_weights=None
    ):
        """
        Returns
        -------
        d2_total : float
            Weighted grouped signature distance.
        grad_d2 : np.ndarray, shape (ncart,)
            Cartesian gradient of d2_total.
        """
        n_rotors = len(grouped_current)
        ncart = grouped_brows[0].shape[1] if n_rotors > 0 else self.impes_coordinate.b_matrix.shape[1]

        if rotor_weights is None:
            rotor_weights = np.ones(n_rotors, dtype=np.float64)

        d2_total = 0.0
        grad_total = np.zeros(ncart, dtype=np.float64)
        w_total = 0.0

        for rid in range(n_rotors):
            sig_c = grouped_current[rid]
            sig_r = grouped_candidate[rid]
            brows = grouped_brows[rid]   # shape (m_r, ncart)

            delta = sig_c - sig_r
            terms = 2.0 * (1.0 - np.cos(delta))
            d2_r = np.mean(terms)

            grad_terms = 2.0 * np.sin(delta)[:, None] * brows
            grad_r = np.mean(grad_terms, axis=0)

            wr = rotor_weights[rid]
            d2_total += wr * d2_r
            grad_total += wr * grad_r
            w_total += wr
        
        d2_total /= max(w_total, 1e-12)
        grad_total /= max(w_total, 1e-12)

        return d2_total, grad_total
    
    def _signature_inverse_weight_and_gradient(
        self,
        d2,
        grad_d2,
        confidence_radius=1.5,
        power=2,
        eps=1e-12
    ):
        """
        Raw inverse-distance-like weight and its gradient.
        """
        rho2 = max(confidence_radius**2, eps)
        x = d2 / rho2 + eps

        w = 1.0 / (2.0 * (x ** power))
        dw_dd2 = -(power / (2.0 * rho2)) * (x ** (-(power + 1)))

        grad_w = dw_dd2 * grad_d2
        return w, grad_w
    
    def _normalized_signature_weights(self, d2_array, confidence_radius=1.0, power=1,
                                  eps=1e-12, exact_tol=1e-10):
        """
        Shepard-like normalized inverse-distance weights.
        """
        d2_array = np.asarray(d2_array, dtype=np.float64)
        idx_best = int(np.argmin(d2_array))
        d2_min = d2_array[idx_best]

        if d2_min < exact_tol:
            w = np.zeros_like(d2_array)
            w[idx_best] = 1.0
            return w

        x = d2_array / max(confidence_radius**2, eps)
        raw = 1.0 / ((x + eps) ** power + (x + eps) ** power)

        S = raw.sum()
        if S < eps:
            w = np.zeros_like(raw)
            w[idx_best] = 1.0
            return w

        return raw / S
    
    def _compute_two_rotor_phase_distance(self, current_phases, reference_phases, rotor_weights=None):
        """
        current_phases:   shape (2,)
        reference_phases: shape (2,)
        rotor_weights:    optional shape (2,)
        """
        delta = current_phases - reference_phases   # shape (2,)
        mism = 2.0 * (1.0 - np.cos(delta))                   # shape (2,)

        
        
        if rotor_weights is not None:
            return np.sum(rotor_weights * mism)
        
        if np.sum(mism) == 0.0:
            return 1e-8
        return np.sum(mism)
    
    def _candidate_current_chart(self, symmetry_data_point, mask, org_int_coords, b_matrix, base_cache=None):
        """
        Return candidate-local current chart (org_int_coords, B).
        masked_exact  -> candidate-local eq masking (exact).
        symmetrized   -> global chart reuse (fast).
        """
        org = np.asarray(org_int_coords, dtype=np.float64)
        bmat = np.asarray(b_matrix, dtype=np.float64)

        if not self.impes_coordinate.use_eq_bond_length:
            return org, bmat

        mode = self._get_eq_symmetry_mode()
        if mode == 'symmetrized':
            return org, bmat

        eq_ref = getattr(symmetry_data_point, 'eq_bond_lengths', None)
        assert_msg_critical(
            eq_ref is not None,
            f"InterpolationDriver: missing eq_bond_lengths for symmetry_data_point "
            f"{getattr(symmetry_data_point, 'point_label', None)}."
        )

        if base_cache is not None:
            return self.impes_coordinate.build_masked_current_chart_from_reference_eq_fast(
                reference_eq_bond_lengths=eq_ref,
                mask=mask,
                base_cache=base_cache,
            )

        # Safe fallback if cache was not prepared.
        return self.impes_coordinate.build_masked_current_chart_from_reference_eq(
            reference_eq_bond_lengths=eq_ref,
            mask=mask,
        )
    
    def _evaluate_candidate_taylor(
            self,
            *,
            energy: float,
            grad: np.ndarray,
            hess: np.ndarray,
            ref_coords: np.ndarray,
            mask0: np.ndarray,
            mask: np.ndarray,
            org_int_coords: np.ndarray,
            b_matrix: np.ndarray,
            use_cosine_dihedral: bool,
            dihedral_start: int,
            dihedral_end: int) -> Tuple[float, np.ndarray, np.ndarray, np.ndarray]:
        """
        Evaluates one masked Taylor candidate in a parity-safe way shared by
        serial and MPI preload paths.
        """
        grad_eff = np.asarray(grad, dtype=np.float64).copy()
        hess_eff = np.asarray(hess, dtype=np.float64).copy()
        ref = np.asarray(ref_coords, dtype=np.float64)
        mask0_arr = np.asarray(mask0, dtype=np.int64)
        mask_arr = np.asarray(mask, dtype=np.int64)
        org = np.asarray(org_int_coords, dtype=np.float64)
        bmat = np.asarray(b_matrix, dtype=np.float64)

        grad_eff[mask0_arr] = grad_eff[mask_arr]
        hess_eff[np.ix_(mask0_arr, mask0_arr)] = hess_eff[np.ix_(mask_arr, mask_arr)]

        dist_org = org - ref[mask_arr]
        dist_check = dist_org.copy()
        if not use_cosine_dihedral:
            dist_check[dihedral_start:dihedral_end] = np.sin(
                dist_org[dihedral_start:dihedral_end])

        dist_correlation = dist_check.copy()

        U_i = float(
            energy
            + dist_check @ grad_eff
            + 0.5 * np.linalg.multi_dot([dist_check, hess_eff, dist_check])
        )

        dist_hess = dist_check @ hess_eff
        if not use_cosine_dihedral:
            c = np.cos(dist_org[dihedral_start:dihedral_end])
            grad_eff[dihedral_start:dihedral_end] *= c
            dist_hess[dihedral_start:dihedral_end] *= c

        natm = bmat.shape[1] // 3
        dU_i = (bmat.T @ (grad_eff + dist_hess)).reshape(natm, 3)

        return U_i, dU_i, dist_org, dist_correlation

    def compute_potential(self, data_point, org_int_coords):
        """Calculates the potential energy surface at self.impes_coordinate
           based on the energy, gradient and Hessian of data_point.

           :param data_point:
                InterpolationDatapoint object.
        """
        
        clustered_out = self._compute_potential_clustered(data_point, org_int_coords)
        if clustered_out is not None:
            return clustered_out
        

        pes = 0.0
        gradient = None
        natm = data_point.cartesian_coordinates.shape[0]
        # print(len(self.qm_symmetry_data_points))
        dp_label = data_point.point_label 
        
        symmetry_data_points = (self.qm_symmetry_data_points or {}).get(dp_label, [])
        use_symmetry_branch = (len(symmetry_data_points) > 1)
        bounds = self._get_internal_coordinate_partitions()
        bond_end = int(bounds['bond_end'])
        angle_end = int(bounds['angle_end'])
        dihedral_start = int(bounds['dihedral_start'])
        dihedral_end = int(bounds['dihedral_end'])

        symmetry_task_entries = ()
        if use_symmetry_branch:
            symmetry_task_entries = self._get_symmetry_task_entries(dp_label)
            use_symmetry_branch = len(symmetry_task_entries) > 0

        if use_symmetry_branch:

            potentials = []
            pot_gradients = []
            hessian_error = []

            torsion_meta = self._get_symmetry_torsion_meta()
            b_matrix = self.impes_coordinate.b_matrix

            eq_mode = None
            candidate_base_cache = None
            if self.impes_coordinate.use_eq_bond_length:
                eq_mode = self._get_eq_symmetry_mode()
                if eq_mode == 'masked_exact':
                    candidate_base_cache = self.impes_coordinate.prepare_eq_candidate_base_cache()

            current_grouped = self._build_grouped_symmetry_signatures(org_int_coords.copy(), torsion_meta)
            grouped_brows = self._build_grouped_symmetry_brows(torsion_meta)
            
            cand_grouped_signitures = []
            candidate_labels = []
  

            # reset diagnostics if needed
            self.bond_rmsd = []
            self.angle_rmsd = []
            self.dihedral_rmsd = []

            for symmetry_data_point, mask0, mask in symmetry_task_entries:
                energy = symmetry_data_point.energy

                candidate_full_grouped = self._build_grouped_symmetry_signatures(
                    symmetry_data_point.internal_coordinates_values[mask],
                    torsion_meta
                )

                cand_grouped_signitures.append(candidate_full_grouped)
                candidate_labels.append((symmetry_data_point.point_label, tuple(mask)))

                org_eval_i, b_eval_i = self._candidate_current_chart(
                    symmetry_data_point=symmetry_data_point,
                    mask=mask,
                    org_int_coords=org_int_coords,
                    b_matrix=b_matrix,
                    base_cache=candidate_base_cache,
                )
                
                pes_i, pes_prime_i, dist_org, dist_correlation = self._evaluate_candidate_taylor(
                    energy=energy,
                    grad=symmetry_data_point.internal_gradient,
                    hess=symmetry_data_point.internal_hessian,
                    ref_coords=symmetry_data_point.internal_coordinates_values,
                    mask0=mask0,
                    mask=mask,
                    org_int_coords=org_eval_i,
                    b_matrix=b_eval_i,
                    use_cosine_dihedral=bool(self.use_cosine_dihedral),
                    dihedral_start=dihedral_start,
                    dihedral_end=dihedral_end,
                )

                self.bond_rmsd.append(
                    np.sqrt(np.mean(np.sum((dist_org[:bond_end])**2)))
                )
                self.angle_rmsd.append(
                    np.sqrt(np.mean(np.sum(dist_org[bond_end:angle_end]**2)))
                )
                self.dihedral_rmsd.append(
                    np.sqrt(np.mean(np.sum(dist_correlation[dihedral_start:dihedral_end]**2)))
                )

                potentials.append(pes_i)
                pot_gradients.append(pes_prime_i)

            raw_w = []
            raw_grad_w = []

            for cand_grouped in cand_grouped_signitures:
                d2_i, grad_d2_i = self._grouped_signature_distance_and_gradient(
                    current_grouped,
                    cand_grouped,
                    grouped_brows
                )
        
                w_i, grad_w_i = self._signature_inverse_weight_and_gradient(
                    d2_i,
                    grad_d2_i,
                    confidence_radius=1.0,
                    power=2,
                    eps=1e-12
                )

                raw_w.append(w_i)
                raw_grad_w.append(grad_w_i.reshape(natm, 3))

            w_i_sig = np.asarray(raw_w, dtype=np.float64)
            grad_w_i_sig = np.asarray(raw_grad_w, dtype=np.float64)

            S_sig = w_i_sig.sum()
            sum_grad_w_sig = grad_w_i_sig.sum(axis=0)

            W_i_sig = w_i_sig / S_sig
            grad_W_i_sig = (grad_w_i_sig * S_sig - w_i_sig[:, None, None] * sum_grad_w_sig) / (S_sig**2)

            U_i = np.asarray(potentials, dtype=np.float64)
            grad_U_i = np.asarray(pot_gradients, dtype=np.float64)

            pes = np.dot(W_i_sig, U_i)

            grad_pes = (
                np.tensordot(W_i_sig, grad_U_i, axes=1)
                + np.tensordot(U_i, grad_W_i_sig, axes=1)
            )

            gradient = grad_pes

            return pes, gradient, hessian_error

        else:
            hessian_error = 0.0
            energy = data_point.energy
            grad = data_point.internal_gradient.copy()
            hessian = data_point.internal_hessian.copy()
            dist_org = (org_int_coords.copy() - data_point.internal_coordinates_values)
            dist_check = (org_int_coords.copy() - data_point.internal_coordinates_values)

            if not self.use_cosine_dihedral:
                dist_check[dihedral_start:dihedral_end] = np.sin(dist_org[dihedral_start:dihedral_end])

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
                cos_dist = np.cos(dist_org[dihedral_start:dihedral_end])
                grad[dihedral_start:dihedral_end] *= cos_dist
                dist_hessian_eff[dihedral_start:dihedral_end] *= cos_dist

            pes_prime = np.matmul(self.impes_coordinate.b_matrix.T, (grad + dist_hessian_eff)).reshape(natm, 3)

            return pes, pes_prime, (grad + dist_hessian_eff)
            

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
    
    
    def calculate_translation_coordinates(self, given_coordinates, weights=None):
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

            
        # sum of the derivative necessary 
        imp_int_coord_derivative_contribution = np.zeros_like(distance_vector_sub).reshape(-1)
        imp_int_coord_distance = 0
        dq_dx = np.zeros_like(distance_vector_sub).reshape(-1)
        Dimp_sq = 0.0
        if self.use_tc_weights:
            for key, entries in data_point.imp_int_coordinates.items():
                for element in entries:
                    idx = self.impes_coordinate.z_matrix.index(element)
                    org_imp_int_coord_distance = (self.impes_coordinate.internal_coordinates_values[idx] - data_point.internal_coordinates_values[idx])
                    unmw_b_matrix = self.impes_coordinate.b_matrix[idx, active_dofs] / self.impes_coordinate.inv_sqrt_masses[active_dofs]
                    # print(key, org_imp_int_coord_distance)
                    if len(element) == 4:
                        sin_delta = np.sin(org_imp_int_coord_distance)
                        cos_delta = np.cos(org_imp_int_coord_distance)
                        imp_int_coord_derivative_contribution += 2.0 * unmw_b_matrix * cos_delta * sin_delta
                        imp_int_coord_distance += sin_delta**2

                        dq_dx +=  2.0 * unmw_b_matrix * np.cos( org_imp_int_coord_distance) * np.sin( org_imp_int_coord_distance) 
                        Dimp_sq += np.sin(org_imp_int_coord_distance)**2
                    else:
                        imp_int_coord_distance += org_imp_int_coord_distance**2
                        imp_int_coord_derivative_contribution += 2.0 * unmw_b_matrix * org_imp_int_coord_distance

                        dq_dx += 2.0 * unmw_b_matrix * org_imp_int_coord_distance
                        Dimp_sq += (org_imp_int_coord_distance)**2
                
            if imp_int_coord_distance < 1e-8:
                imp_int_coord_distance = 1e-8
                imp_int_coord_derivative_contribution[:] = 0
        
        
        dw_dalhpa_i = 0
        dw_dX_dalpha_i = 0

        if self.interpolation_type == 'shepard':
            
            denominator_imp_coord, weight_gradient_sub_imp_coord = self.VL_target_customized_shepard_weight_gradient(
                distance_vector_sub, distance, data_point.confidence_radius, Dimp_sq, dq_dx)
            
            if self.calc_optim_trust_radius:

                dw_dalhpa_i = self.trust_radius_tc_weight_gradient(
                    data_point.confidence_radius, distance, Dimp_sq)
                dw_dX_dalpha_i = self.trust_radius_tc_weight_gradient_gradient(
                    data_point.confidence_radius,
                    distance,
                    distance_vector_sub,
                    Dimp_sq,
                    dq_dx)

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
        max_block_partners = 4
        min_partner_pair_frac = 0.15
        min_partner_coord_frac = 0.20
        pair_weight = 0.35
        support_weight = 0.25

        # final selection from best block
        max_constraints_to_return = max(1, int(len(molecule.get_labels()) / 3))
        singleton_dominance = 1.35
        singleton_source_frac = 0.85
        singleton_response_frac = 0.85
        secondary_score_frac = 0.25
        min_secondary_pair_frac = 0.15

        eps = 1e-12

        coords_flat, kinds_flat = self._flatten_internal_coordinates(z_matrix)
        N = len(coords_flat)

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

        for idx in keep_indices:
            wk = dp_weights[idx]
            r = per_dp_results[idx]

            global_response += wk * r["response_score"]
            global_source_blame += wk * r["source_blame_score"]
            global_helpful += wk * r["helpful_source_score"]
            global_support += wk * r["support_mask"]
            global_rho += wk * r["rho"]
            global_pair += wk * r["pair_score"]

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
        ranked_blocks = []
        seen_blocks = set()

        for anchor_idx in candidate_anchor_indices:
            anchor_coord = tuple(int(x) for x in coords_flat[anchor_idx])

            if anchor_coord in constraints_to_exclude:
                continue

            if (
                kinds_flat[anchor_idx] == "dihedral"
                and tuple(anchor_coord[1:3]) in self._get_symmetric_dihedral_centers(self.symmetry_information)
            ):
                continue

            block_indices = self._build_candidate_block(
                anchor_idx=anchor_idx,
                global_coord_score=global_coord_score,
                global_pair=global_pair,
                z_matrix=coords_flat,
                coord_kinds=kinds_flat,
                constraints_to_exclude=constraints_to_exclude,
                max_partners=max_block_partners,
                min_partner_pair_frac=min_partner_pair_frac,
                min_partner_coord_frac=min_partner_coord_frac,
            )

            key = tuple(sorted(block_indices))
            if key in seen_blocks:
                continue
            seen_blocks.add(key)

            coord_mass = float(np.sum(global_coord_score[block_indices]))
            support_mass = float(np.mean(global_support[block_indices])) if block_indices else 0.0

            pair_mass = 0.0
            for a in range(len(block_indices)):
                for b in range(a + 1, len(block_indices)):
                    i = block_indices[a]
                    j = block_indices[b]
                    pair_mass += global_pair[i, j]

            acq_score = coord_mass + pair_weight * pair_mass + support_weight * support_mass

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
                "pair_mass": float(pair_mass),
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
            coord_kinds=kinds_flat,
            max_constraints=max_constraints_to_return,
            singleton_dominance=singleton_dominance,
            singleton_source_frac=singleton_source_frac,
            singleton_response_frac=singleton_response_frac,
            secondary_score_frac=secondary_score_frac,
            min_secondary_pair_frac=min_secondary_pair_frac,
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

        return primary_constraint, candidate_constraints, best_block["block_coords"]


    # ======================================================================
    # Helper methods
    # ======================================================================

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


    def _get_symmetric_dihedral_centers(self, symmetry_information) -> set:
        """
        Return central atom pairs of torsions that should be symmetry-excluded.
        """
        try:
            return set(tuple(int(y) for y in x) for x in symmetry_information[7][3])
        except Exception:
            return set()


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
        Returns
        -------
        dq_raw : raw internal differences
        dq_eff : effective differences used in the interpolation model
        chain  : derivative chain factors

        Proper dihedrals:
            dq_raw wrapped to (-pi, pi]
            dq_eff = sin(dq_raw)
            chain  = cos(dq_raw)

        Impropers:
            dq_raw wrapped to (-pi, pi]
            dq_eff = dq_raw
            chain  = 1

        Bonds/angles:
            dq_eff = dq_raw
            chain  = 1
        """
        q_current = np.asarray(q_current, dtype=float)
        q_ref = np.asarray(q_ref, dtype=float)

        coords_flat, kinds_flat = self._flatten_internal_coordinates(z_matrix)
        N = len(coords_flat)

        dq_raw = np.zeros(N, dtype=float)
        dq_eff = np.zeros(N, dtype=float)
        chain = np.ones(N, dtype=float)

        sym_torsion_centers = self._get_symmetric_dihedral_centers(symmetry_information)

        for idx, (coord, kind) in enumerate(zip(coords_flat, kinds_flat)):
            if kind == "dihedral":
                if tuple(coord[1:3]) in sym_torsion_centers:
                    dq_raw[idx] = 0.0
                    dq_eff[idx] = 0.0
                    chain[idx] = 0.0
                    continue

                d = self._principal_torsion_delta(q_current[idx] - q_ref[idx])
                dq_raw[idx] = d
                dq_eff[idx] = np.sin(d)
                chain[idx] = np.cos(d)

            elif kind == "improper":
                d = self._principal_torsion_delta(q_current[idx] - q_ref[idx])
                dq_raw[idx] = d
                dq_eff[idx] = d
                chain[idx] = 1.0

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
            if kinds_flat[j] == "dihedral":
                # Keep periodic consistency.
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
            "source_contrib": blame_matrix,  # semantic reuse for compatibility
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

        sym_torsion_centers = self._get_symmetric_dihedral_centers(self.symmetry_information)

        dq_eff = np.zeros(N, dtype=float)
        chain = np.ones(N, dtype=float)

        for idx, (coord, kind) in enumerate(zip(coords_flat, kinds_flat)):
            if kind == "dihedral":
                if tuple(coord[1:3]) in sym_torsion_centers:
                    dq_eff[idx] = 0.0
                    chain[idx] = 0.0
                else:
                    dq_eff[idx] = np.sin(dq_raw[idx])
                    chain[idx] = np.cos(dq_raw[idx])

            elif kind == "improper":
                dq_eff[idx] = dq_raw[idx]
                chain[idx] = 1.0

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
        sym_torsion_centers = self._get_symmetric_dihedral_centers(self.symmetry_information)

        candidates = []
        for j in range(N):
            if j == anchor_idx:
                continue

            coord_j = tuple(int(x) for x in z_matrix[j])

            if coord_j in constraints_to_exclude:
                continue

            if coord_kinds[j] == "dihedral" and tuple(coord_j[1:3]) in sym_torsion_centers:
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
        coord_kinds: List[str],
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
    
