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
from contextlib import redirect_stderr
from io import StringIO
import json
import numpy as np
import h5py
from time import time
from sys import stdout

from .errorhandler import assert_msg_critical
from .inputparser import (parse_input, print_keywords, print_attributes)
from .outputstream import OutputStream
from .veloxchemlib import mpi_master


with redirect_stderr(StringIO()) as fg_err:
    import geometric


class InterpolationDatapoint:
    """
    Implements routines for the coordinates required for potential energy
    surface interpolation.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - hessian: The Cartesian Hessian in Hartree per Bohr**2.
        - energy: the total energy of the groound or excited state in Hartree
        - gradient: The Cartesian gradient in Hartree per Bohr
        - z_matrix: The Z-matrix (a list of tuples with indices of the atoms
          which define bonds, bond angles and dihedral angles).
        - internal_gradient: The gradient in internal coordinates.
        - internal_hessian: The Hessian in internal coordinates.
        - b_matrix: The Wilson B-matrix.
        - b2_matrix: The matrix of second order deriatives of
          internal coordinates with respect to Cartesian coordinates.
        - use_inverse_bond_length: Flag for using the inverse bond length
          instead of the bond length.
        - internal_coordinates: A list of geomeTRIC internal coordinate objects.
        - internal_coordinates_values: A numpy array with the values of the
          internal coordinates: in angstroms (bond lengths), 1/angstroms
          (if inverse bond lengths are used), radians (bond angles), radians or
          cosine of radians (dihedral angles).
    """

    def __init__(self, z_matrix=None, atom_labels=None, comm=None, ostream=None):
        """
        Initializes the InterpolationDatapoint object.

        :param z_matrix: a list of tuples with atom indices (starting at 0)
        that define bonds, bond angles, and dihedral angles.
        """
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(stdout)
            else:
                ostream = OutputStream(None)

        self.comm = None
        self.ostream = ostream

        self.point_label = None
        self.metadata = {}

        self.energy = None
        self.gradient = None
        self.hessian = None
        self.z_matrix_dict = z_matrix
        if z_matrix is not None:
            z_matrix_flat = self.flatten_z_matrix(z_matrix)
            self.z_matrix = z_matrix_flat
        self.atom_labels = atom_labels
        self.internal_gradient = None
        self.internal_hessian = None
        self.inv_sqrt_masses = None

        self.b_matrix = None

        # copy of Willson B matrix in r and phi
        # (needed for converting the b2_matrix in 1/r and/or cos(phi)).
        self.original_b_matrix = None
        self.b2_matrix = None

        self.confidence_radius = None
        self.use_inverse_bond_length = True
        self.use_eq_bond_length = False
        self.bond_coordinate_mode = "legacy"
        self.angle_coordinate_mode = "raw"
        self.use_eq_angle_cosine = False
        self.eq_bond_symmetry_mode = "masked_exact"
        # self.use_cos_angle = False
        self.use_vectorized_b_matrix = True
        self.use_vectorized_internal_coordinates_values = True
        self.identify_imp_int_coord = True

        self.eq_bond_lengths = None
        self.eq_angle_values = None

        # internal_coordinates is a list of geomeTRIC objects which represent
        # different types of internal coordinates (distances, angles, dihedrals)
        self.internal_coordinates = None
        self.imp_int_coordinates = []
        self.mapping_masks = None

        self.family_label = None
        self.bank_role = "global"
        self.cluster_id = None
        self.cluster_type = None
        self.cluster_rotor_ids = None
        self.cluster_state_id = 0
        self.dihedrals_to_rotate = None
        self.phase_signature = None
        self.active_atoms = None
        self.active_rows = None
        self.canonical_subset_key = ()
        self.response_rows = ()
        self.relaxation_policy_id = None
        self.is_anchor = False
        self.anchor_state_ids = ()
        self.local_factor_combination_mode = "signed_full"
        self.local_factor_overlap_source = None

        # internal_coordinates_values is a numpy array with the values of the
        # geomeTRIC internal coordinates. It is used in InterpolationDriver to
        # determine the distance in internal coordinates between two
        # InterpolationDatapoint objects.
        self.internal_coordinates_values = None
        self.cartesian_coordinates = None
        self._b_matrix_row_cache = None

        self._input_keywords = {
            'im_settings': {
                'use_inverse_bond_length':
                    ('bool',
                     'use the inverse bond lengths'),
                'use_eq_bond_length':
                    ('bool',
                     'use the log bond lengths'),
                'bond_coordinate_mode':
                    ('str',
                     'bond coordinate mode: legacy, r, inverse, eq, '
                     'eq_no_factor, relative, log_relative, scaled_inverse'),
                'angle_coordinate_mode':
                    ('str',
                     'angle coordinate mode: raw, delta, sin_delta, eq_cosine'),
                'use_eq_angle_cosine':
                    ('bool',
                     'use equilibrium-anchored cosine angle coordinates'),
                'eq_bond_symmetry_mode':
                    ('str', 'eq-bond symmetry mode: masked_exact or symmetrized'),
                'use_vectorized_b_matrix':
                    ('bool', 'use vectorized Wilson B-matrix construction'
                     ),
                # 'use_cos_angle':('bool', 'use cos angles -- careful for linear angles'),
                'use_vectorized_internal_coordinates_values':
                    ('bool', 'use vectorized internal-coordinate value evaluation'
                     ),
            }
        }

    @staticmethod
    def flatten_z_matrix(zm_dict):
        """Flatten canonical dict to legacy ordered list representation."""
        return (
            list(zm_dict["bonds"]) +
            list(zm_dict["angles"]) +
            list(zm_dict["dihedrals"]) +
            list(zm_dict["impropers"])
        )

    def print_keywords(self):
        """
        Prints the input keywords of the InterpolationDatapoint.
        """

        print_keywords(self._input_keywords, self.ostream)

    def print_attributes(self):
        """
        Prints the attributes of the InterpolationDatapoint.
        """

        print_attributes(self._input_keywords, self.ostream)

    def update_settings(self, impes_dict=None):
        """
        Updates settings for InterpolationDatapoint.

        :param impes_dict:
            The input dictionary of settings for interpolated mechanics.
        """

        if impes_dict is None:
            impes_dict = {}

        im_keywords = {
            key: val[0]
            for key, val in self._input_keywords['im_settings'].items()
        }

        parse_input(self, im_keywords, impes_dict)
        mode = self._normalize_bond_coordinate_mode(
            getattr(self, "bond_coordinate_mode", "legacy")
        )
        if mode != "legacy":
            self.use_inverse_bond_length = (mode == "inverse")
            self.use_eq_bond_length = self._bond_mode_requires_eq_reference(mode)

    def define_internal_coordinates(self):
        """
        Defines the internal coordinates from the Z-matrix.
        """

        assert_msg_critical(self.z_matrix is not None, 'InterpolationDatapoint: No Z-matrix defined.')
        self.internal_coordinates = []
        for z in self.z_matrix:

            if len(z) == 2:
                q = geometric.internal.Distance(*z)
            elif len(z) == 3:
                q = geometric.internal.Angle(*z)
            elif len(z) == 4:
                q = geometric.internal.Dihedral(*z)
            else:
                assert_msg_critical(False, 'InterpolationDatapoint: Invalid entry size in Z-matrix.')
            self.internal_coordinates.append(q)

    def _get_b_matrix_row_cache(self):
        """
        Builds/caches static row metadata for vectorized B-matrix assembly.
        """

        if self._b_matrix_row_cache is not None:
            return self._b_matrix_row_cache

        n_rows = len(self.z_matrix)
        row_sizes = np.fromiter((len(z) for z in self.z_matrix),
                                dtype=np.int8,
                                count=n_rows)
        bond_rows = np.flatnonzero(row_sizes == 2)
        angle_rows = np.flatnonzero(row_sizes == 3)
        dihedral_rows = np.flatnonzero(row_sizes == 4)

        dihedral_first = np.zeros(dihedral_rows.shape[0], dtype=bool)
        prev_dihedral = None
        dihedral_counter = 0
        for z in self.z_matrix:
            if len(z) != 4:
                continue
            if prev_dihedral != z:
                prev_dihedral = z
                dihedral_first[dihedral_counter] = True
            dihedral_counter += 1

        self._b_matrix_row_cache = {
            'row_sizes': row_sizes,
            'bond_rows': bond_rows,
            'angle_rows':  angle_rows,
            'dihedral_rows': dihedral_rows,
            'dihedral_first': dihedral_first,
        }

        return self._b_matrix_row_cache

    def _get_eq_bond_lengths_array(self, n_bonds):
        """
        Returns equilibrium bond lengths with shape validation.
        """
        assert_msg_critical(
            self.eq_bond_lengths is not None,
            'InterpolationDatapoint: No equilibrium bond lengths are defined.'
        )
        eq_values = np.asarray(self.eq_bond_lengths, dtype=np.float64).reshape(-1)
        assert_msg_critical(
            eq_values.size == int(n_bonds),
            'InterpolationDatapoint: Equilibrium bond-length size mismatch '
            f'(expected {int(n_bonds)}, got {eq_values.size}).'
        )
        return eq_values

    def _normalize_bond_coordinate_mode(self, mode):
        aliases = {
            None: "legacy",
            "": "legacy",
            "legacy": "legacy",
            "r": "r",
            "raw": "r",
            "plain": "r",
            "plain_r": "r",
            "inverse": "inverse",
            "rinv": "inverse",
            "1/r": "inverse",
            "eq": "eq",
            "equilibrium": "eq",
            "eq_bond": "eq",
            "eq_no_factor": "eq_no_factor",
            "eqnofactor": "eq_no_factor",
            "inverse_shift": "eq_no_factor",
            "relative": "relative",
            "relative_displacement": "relative",
            "log_relative": "log_relative",
            "log": "log_relative",
            "scaled_inverse": "scaled_inverse",
            "req_over_r_minus_1": "scaled_inverse",
        }
        key = str(mode).strip().lower() if mode is not None else None
        if key not in aliases:
            raise ValueError(
                "InterpolationDatapoint: unknown bond_coordinate_mode="
                f"'{mode}'."
            )
        return aliases[key]

    def get_bond_coordinate_mode(self):
        """
        Return the active bond coordinate mode.

        The default ``legacy`` value preserves the historical boolean API:
        ``use_inverse_bond_length=True`` gives ``inverse``; otherwise
        ``use_eq_bond_length=True`` gives the current equilibrium-bond
        coordinate; otherwise raw ``r`` is used.
        """
        mode = self._normalize_bond_coordinate_mode(
            getattr(self, "bond_coordinate_mode", "legacy")
        )
        if mode == "legacy":
            if self.use_inverse_bond_length:
                return "inverse"
            if self.use_eq_bond_length:
                return "eq"
            return "r"
        return mode

    def _bond_mode_requires_eq_reference(self, mode=None):
        mode = self.get_bond_coordinate_mode() if mode is None else mode
        return mode in {
            "eq",
            "eq_no_factor",
            "relative",
            "log_relative",
            "scaled_inverse",
        }

    def coordinate_label_suffix(self):
        """
        Return the HDF5 suffix for the active bond coordinate description.

        Existing modes keep their original suffixes. New experimental modes get
        explicit suffixes so multiple databases can coexist.
        """
        mode = self.get_bond_coordinate_mode()
        suffixes = {
            "r": "_r",
            "inverse": "_rinv",
            "eq": "_eq",
            "eq_no_factor": "_eqnofactor",
            "relative": "_rel",
            "log_relative": "_logrel",
            "scaled_inverse": "_scaledinv",
        }
        return suffixes[mode]

    def _get_eq_angle_values_array(self, n_angles):
        """
        Returns reference angle values with shape validation.
        """
        assert_msg_critical(
            self.eq_angle_values is not None,
            'InterpolationDatapoint: No equilibrium angle values are defined.'
        )
        eq_values = np.asarray(self.eq_angle_values, dtype=np.float64).reshape(-1)
        assert_msg_critical(
            eq_values.size == int(n_angles),
            'InterpolationDatapoint: Equilibrium angle size mismatch '
            f'(expected {int(n_angles)}, got {eq_values.size}).'
        )
        return eq_values

    def _normalize_angle_coordinate_mode(self, mode):
        aliases = {
            None: "raw",
            "": "raw",
            "raw": "raw",
            "theta": "raw",
            "radian": "raw",
            "delta": "delta",
            "theta_minus_theta_eq": "delta",
            "sin_delta": "sin_delta",
            "sin(theta-theta_eq)": "sin_delta",
            "eq_cosine": "eq_cosine",
            "eq_angle": "eq_cosine",
            "cosine": "eq_cosine",
        }
        key = str(mode).strip().lower() if mode is not None else None
        if key not in aliases:
            raise ValueError(
                "InterpolationDatapoint: unknown angle_coordinate_mode="
                f"'{mode}'."
            )
        return aliases[key]

    def get_angle_coordinate_mode(self):
        """
        Return the active angle coordinate mode.

        ``use_eq_angle_cosine`` is an alias for the benchmark-local
        equilibrium-anchored cosine mode:

            q = (cos(theta_eq) - cos(theta)) / sin(theta_eq)

        which has dq/dtheta = sin(theta) / sin(theta_eq).
        """
        if getattr(self, "use_eq_angle_cosine", False):
            return "eq_cosine"
        return self._normalize_angle_coordinate_mode(
            getattr(self, "angle_coordinate_mode", "raw")
        )

    def _angle_mode_requires_eq_reference(self, mode=None):
        mode = self.get_angle_coordinate_mode() if mode is None else mode
        return mode in {"delta", "sin_delta", "eq_cosine"}

    def _bond_transform(self, r, r_eq=None):
        """
        Returns q(r), dq/dr, d2q/dr2 for the active bond coordinate.

        Modes:
            r              : q = r
            inverse        : q = 1/r
            eq             : q = -r_eq^2 * (1/r - 1/r_eq)
            eq_no_factor   : q = -(1/r - 1/r_eq)
            relative       : q = (r - r_eq) / r_eq
            log_relative   : q = log(r / r_eq)
            scaled_inverse : q = r_eq/r - 1
        """
        mode = self.get_bond_coordinate_mode()
        r_input = np.asarray(r, dtype=np.float64)
        scalar_input = (r_input.ndim == 0)
        r_arr = np.atleast_1d(r_input).astype(np.float64, copy=False)

        if not self._bond_mode_requires_eq_reference(mode):
            if mode == "r":
                q = r_arr.copy()
                dq_dr = np.ones_like(r_arr)
                d2q_dr2 = np.zeros_like(r_arr)
            elif mode == "inverse":
                q = 1.0 / r_arr
                dq_dr = -1.0 / np.square(r_arr)
                d2q_dr2 = 2.0 / np.power(r_arr, 3)
            else:
                raise ValueError(f"Unsupported bond coordinate mode: {mode}")
        else:
            assert_msg_critical(
                r_eq is not None,
                'InterpolationDatapoint: Bond coordinate mode requires r_eq.'
            )
            r_eq_input = np.asarray(r_eq, dtype=np.float64)
            scalar_input = scalar_input and (r_eq_input.ndim == 0)
            r_arr, r_eq_arr = np.broadcast_arrays(
                np.atleast_1d(r_input),
                np.atleast_1d(r_eq_input),
            )
            r_arr = r_arr.astype(np.float64, copy=False)
            r_eq_arr = r_eq_arr.astype(np.float64, copy=False)

            if mode == "eq":
                q = -np.square(r_eq_arr) * (1.0 / r_arr - 1.0 / r_eq_arr)
                dq_dr = np.square(r_eq_arr) / np.square(r_arr)
                d2q_dr2 = -2.0 * np.square(r_eq_arr) / np.power(r_arr, 3)
            elif mode == "eq_no_factor":
                q = -(1.0 / r_arr - 1.0 / r_eq_arr)
                dq_dr = 1.0 / np.square(r_arr)
                d2q_dr2 = -2.0 / np.power(r_arr, 3)
            elif mode == "relative":
                q = (r_arr - r_eq_arr) / r_eq_arr
                dq_dr = 1.0 / r_eq_arr
                d2q_dr2 = np.zeros_like(r_arr)
            elif mode == "log_relative":
                q = np.log(r_arr / r_eq_arr)
                dq_dr = 1.0 / r_arr
                d2q_dr2 = -1.0 / np.square(r_arr)
            elif mode == "scaled_inverse":
                q = r_eq_arr / r_arr - 1.0
                dq_dr = -r_eq_arr / np.square(r_arr)
                d2q_dr2 = 2.0 * r_eq_arr / np.power(r_arr, 3)
            else:
                raise ValueError(f"Unsupported bond coordinate mode: {mode}")

        if scalar_input:
            return float(q[0]), float(dq_dr[0]), float(d2q_dr2[0])

        return q, dq_dr, d2q_dr2

    def _angle_transform(self, theta, theta_eq=None):
        """
        Returns q(theta), dq/dtheta, d2q/dtheta2 for the active angle mode.

        Modes:
            raw       : q = theta
            delta     : q = theta - theta_eq
            sin_delta : q = sin(theta - theta_eq)
            eq_cosine : q = (cos(theta_eq) - cos(theta)) / sin(theta_eq)
        """
        mode = self.get_angle_coordinate_mode()
        theta_input = np.asarray(theta, dtype=np.float64)
        scalar_input = (theta_input.ndim == 0)
        theta_arr = np.atleast_1d(theta_input).astype(np.float64, copy=False)

        if not self._angle_mode_requires_eq_reference(mode):
            q = theta_arr.copy()
            dq = np.ones_like(theta_arr)
            d2q = np.zeros_like(theta_arr)
        else:
            assert_msg_critical(
                theta_eq is not None,
                'InterpolationDatapoint: Angle coordinate mode requires theta_eq.'
            )
            theta_eq_input = np.asarray(theta_eq, dtype=np.float64)
            scalar_input = scalar_input and (theta_eq_input.ndim == 0)
            theta_arr, theta_eq_arr = np.broadcast_arrays(
                np.atleast_1d(theta_input),
                np.atleast_1d(theta_eq_input),
            )
            theta_arr = theta_arr.astype(np.float64, copy=False)
            theta_eq_arr = theta_eq_arr.astype(np.float64, copy=False)

            if mode == "delta":
                q = theta_arr - theta_eq_arr
                dq = np.ones_like(theta_arr)
                d2q = np.zeros_like(theta_arr)
            elif mode == "sin_delta":
                delta = theta_arr - theta_eq_arr
                q = np.sin(delta)
                dq = np.cos(delta)
                d2q = -np.sin(delta)
            elif mode == "eq_cosine":
                denom = np.sin(theta_eq_arr)
                safe = np.abs(denom) > 1.0e-8
                q = theta_arr.copy()
                dq = np.ones_like(theta_arr)
                d2q = np.zeros_like(theta_arr)
                q[safe] = (
                    np.cos(theta_eq_arr[safe]) -
                    np.cos(theta_arr[safe])
                ) / denom[safe]
                dq[safe] = np.sin(theta_arr[safe]) / denom[safe]
                d2q[safe] = np.cos(theta_arr[safe]) / denom[safe]
            else:
                raise ValueError(f"Unsupported angle coordinate mode: {mode}")

        if scalar_input:
            return float(q[0]), float(dq[0]), float(d2q[0])

        return q, dq, d2q

    def _calculate_b_matrix_vectorized(self):
        """
        Vectorized Wilson B-matrix construction.
        """
        assert_msg_critical(self.internal_coordinates is not None, 'InterpolationDatapoint: No internal coordinates are defined.')
        assert_msg_critical(self.cartesian_coordinates is not None, 'InterpolationDatapoint: No cartesian coordinates are defined.')

        n_atoms = self.cartesian_coordinates.shape[0]
        n_cart = n_atoms * 3
        n_rows = len(self.z_matrix)
        coords = self.cartesian_coordinates.reshape((n_cart))

        derivatives = np.zeros((n_rows, n_cart))
        for i, q in enumerate(self.internal_coordinates):
            derivative = q.derivative(coords).reshape(-1)
            derivatives[i, :derivative.shape[0]] = derivative

        bond_mode = self.get_bond_coordinate_mode()
        angle_mode = self.get_angle_coordinate_mode()
        use_original = (
            bond_mode != "r" or
            angle_mode != "raw"
        )
        
        if use_original:
            self.original_b_matrix = derivatives.copy()
        else:
            self.original_b_matrix = None

        row_scale = np.ones(n_rows, dtype=np.float64)
        row_cache = self._get_b_matrix_row_cache()

        bond_rows = row_cache['bond_rows']
        angle_rows = row_cache['angle_rows']
        if bond_rows.size > 0 and bond_mode != "r":
            bond_values = np.array(
                [self.internal_coordinates[idx].value(coords) for idx in bond_rows],
                dtype=np.float64)
            eq_values = None
            if self._bond_mode_requires_eq_reference(bond_mode):
                eq_values = self._get_eq_bond_lengths_array(bond_rows.size)
            _, dq_dr, _ = self._bond_transform(bond_values, eq_values)
            row_scale[bond_rows] = dq_dr

        if angle_rows.size > 0 and angle_mode != "raw":
            angle_values = np.array(
                [self.internal_coordinates[idx].value(coords) for idx in angle_rows],
                dtype=np.float64)
            eq_angle_values = self._get_eq_angle_values_array(angle_rows.size)
            _, dq_dtheta, _ = self._angle_transform(
                angle_values,
                eq_angle_values,
            )
            row_scale[angle_rows] = dq_dtheta
        
        self.b_matrix = derivatives * row_scale[:, np.newaxis]

        if self.inv_sqrt_masses is not None:
            self.b_matrix = self.b_matrix * self.inv_sqrt_masses

    def calculate_b_matrix(self):
        """
        Calculates the Wilson B-matrix.
        """
        self._calculate_b_matrix_vectorized()

    def calculate_b2_matrix(self):
        """
        Calculates the second derivative matrix of the internal coordinates
        """

        assert_msg_critical(self.internal_coordinates is not None, 'InterpolationDatapoint: No internal coordinates are defined.')
        assert_msg_critical(self.cartesian_coordinates is not None, 'InterpolationDatapoint: No cartesian coordinates are defined.')

        n_atoms = self.cartesian_coordinates.shape[0]
        coords = self.cartesian_coordinates.reshape((n_atoms * 3))

        # self.b2_matrix = np.zeros((len(self.z_matrix), n_atoms * 3, n_atoms * 3))
        self.b2_matrix = np.zeros((len(self.z_matrix), n_atoms * 3, n_atoms * 3))

        bond_mode = self.get_bond_coordinate_mode()
        angle_mode = self.get_angle_coordinate_mode()

        eq_bond_values = None
        if self._bond_mode_requires_eq_reference(bond_mode):
            bond_rows = self._get_b_matrix_row_cache()['bond_rows']
            eq_bond_values = self._get_eq_bond_lengths_array(bond_rows.size)

        eq_angle_values = None
        if self._angle_mode_requires_eq_reference(angle_mode):
            angle_rows = self._get_b_matrix_row_cache()['angle_rows']
            eq_angle_values = self._get_eq_angle_values_array(angle_rows.size)

        # prev_dihedral = None
        bond_counter = 0
        angle_counter = 0
        for i, z in enumerate(self.z_matrix):
            q = self.internal_coordinates[i]
            second_derivative = q.second_derivative(coords).reshape(-1, n_atoms * 3)

            if len(z) == 2:
                if bond_mode != "r":
                    r = q.value(coords)
                    r_eq = (
                        eq_bond_values[bond_counter]
                        if eq_bond_values is not None
                        else None
                    )
                    _, dq_dr, d2q_dr2 = self._bond_transform(r, r_eq)
                    self.b2_matrix[i] = dq_dr * second_derivative
                    for m in range(n_atoms):
                        for n in range(n_atoms):
                            self.b2_matrix[
                                i, m*3:(m+1)*3, n*3:(n+1)*3
                            ] += d2q_dr2 * np.outer(
                                self.original_b_matrix[i, m*3:(m+1)*3],
                                self.original_b_matrix[i, n*3:(n+1)*3],
                            )

                else:
                    self.b2_matrix[i] = second_derivative

                bond_counter += 1

            elif len(z) == 3:
                if angle_mode != "raw":
                    theta = q.value(coords)
                    _, dq_dtheta, d2q_dtheta2 = self._angle_transform(
                        theta,
                        eq_angle_values[angle_counter],
                    )
                    self.b2_matrix[i] = dq_dtheta * second_derivative
                    for m in range(n_atoms):
                        for n in range(n_atoms):
                            self.b2_matrix[
                                i, m*3:(m+1)*3, n*3:(n+1)*3
                            ] += d2q_dtheta2 * np.outer(
                                self.original_b_matrix[i, m*3:(m+1)*3],
                                self.original_b_matrix[i, n*3:(n+1)*3],
                            )
                else:
                    self.b2_matrix[i] = second_derivative

                angle_counter += 1

            else:
                self.b2_matrix[i] = second_derivative

        if self.inv_sqrt_masses is not None:
            self.b2_matrix = (self.b2_matrix *
                              self.inv_sqrt_masses.reshape(1, -1, 1) *
                              self.inv_sqrt_masses.reshape(1, 1, -1))

    def transform_gradient_to_internal_coordinates(self, tol=1e-6):
        """
        Transforms the gradient from Cartesian to internal coordinates.

        :param tol:
            Tolerance for the singular values of the B matrix.

        """

        dimension = self.gradient.shape[0] * 3 - 6
        if self.gradient.shape[0] == 2:
            dimension += 1

        if self.b_matrix is None:
            self.calculate_b_matrix()

        assert_msg_critical(self.gradient is not None, 'InterpolationDatapoint: No gradient is defined.')

        g_matrix = np.dot(self.b_matrix, self.b_matrix.T)

        U, s, Vt = np.linalg.svd(g_matrix)

        # Make zero the values of s_inv that are smaller than tol
        s_inv = np.array([1 / s_i if s_i > tol else 0.0 for s_i in s])

        g_minus_matrix = np.dot(U, np.dot(np.diag(s_inv), Vt))

        gradient_flat = self.gradient.flatten()
        self.internal_gradient = np.dot(g_minus_matrix, np.dot(self.b_matrix, gradient_flat))

    def transform_hessian_to_internal_coordinates(self, tol=1e-6):
        """
        Transforms the Hessian from Cartesian to internal coordinates.

        :param tol:
            Tolerance for the singular values of the B matrix.

        """

        if self.internal_gradient is None:
            self.transform_gradient_to_internal_coordinates()

        if self.b2_matrix is None:
            self.calculate_b2_matrix()

        dimension = self.gradient.shape[0] * 3 - 6

        if self.gradient.shape[0] == 2:
            dimension += 1

        g_matrix = np.dot(self.b_matrix, self.b_matrix.T)

        U, s, Vt = np.linalg.svd(g_matrix)

        # Make zero the values of s_inv that are smaller than tol
        s_inv = np.array([1 / s_i if s_i > tol else 0.0 for s_i in s])

        # number_of_positive_values = np.count_nonzero(s_inv)

        g_minus_matrix = np.dot(U, np.dot(np.diag(s_inv), Vt))

        b2_gradient = np.einsum("qxy,q->xy", self.b2_matrix, self.internal_gradient)
        self.b2_gradient = b2_gradient
        intermediate_matrix = np.dot(self.b_matrix, self.hessian - b2_gradient)

        self.internal_hessian = np.linalg.multi_dot([
            g_minus_matrix, intermediate_matrix, self.b_matrix.T, g_minus_matrix.T
        ])
        self.internal_hessian = 0.5 * (self.internal_hessian + self.internal_hessian.T)

    def backtransform_internal_gradient_to_cartesian_coordinates(self):
        '''
        Performs the back-transformation of the internal gradient to the
        Cartesian space.

        :returns:
            The gradient in Cartesian coordinates.
        '''

        cartesian_gradient = np.linalg.multi_dot([self.b_matrix.T, self.internal_gradient]).reshape(self.cartesian_coordinates.shape[0], 3)

        return cartesian_gradient

    def backtransform_internal_hessian_to_cartesian_coordinates(self):
        '''
        Performs the back-transformation of the internal hessian to the
        Cartesian space.

        :returns:
            The hessian in Cartesian coordinates.

        '''

        cartesian_hessian = np.linalg.multi_dot([self.b_matrix.T, self.internal_hessian, self.b_matrix]) + self.b2_gradient

        return cartesian_hessian

    def backtransform_gradient_and_hessian(self):
        """
        Performs all steps required to transform from internal coordinates
        to the Cartesian coordinates.

        """
        assert_msg_critical(self.internal_gradient is not None, 'InterpolationDatapoint: No internal gradient is defined.')
        assert_msg_critical(self.internal_hessian is not None, 'InterpolationDatapoint: No internal hessian is defined.')

        # Transform the gradient and Hessian to internal coordinates
        self.backtransform_internal_gradient_to_cartesian_coordinates()
        self.backtransform_internal_hessian_to_cartesian_coordinates()

    def transform_gradient(self):
        """
        Performs all steps required to transform from Cartesian coordinates
        to the internal coordinates defined in the Z-matrix.

        """

        # Create the internal coordinates through geomeTRIC
        if self.internal_coordinates is None:
            self.define_internal_coordinates()

        # Transform the gradient and Hessian to internal coordinates
        self.transform_gradient_to_internal_coordinates()

        # Save the values of the internal coordinates as a numpy array
        self.compute_internal_coordinates_values()

    def transform_gradient_and_hessian(self):
        """
        Performs all steps required to transform from Cartesian coordinates
        to the internal coordinates defined in the Z-matrix.

        """

        # Create the internal coordinates through geomeTRIC
        if self.internal_coordinates is None:
            self.define_internal_coordinates()

        # Transform the gradient and Hessian to internal coordinates
        self.transform_gradient_to_internal_coordinates()
        self.transform_hessian_to_internal_coordinates()

        # Save the values of the internal coordinates as a numpy array
        self.compute_internal_coordinates_values()

    def reset_coordinates_impes_driver(self, cartesian_coordinates, timing_info=None):
        """
        Resets the Cartesian coordinates for the InterpolationDriver.

        :param cartesian_coordinates:
                a numpy array of the new Cartesian coordinates.
        """
        self.cartesian_coordinates = cartesian_coordinates

        # Since the Cartesian coordinates have been reset,
        # the internal coordinates and transformation matrices
        # must be redifined.
        self.internal_coordinates_values = None
        self.b_matrix = None
        self.g_minus = None

        if self.internal_coordinates is None:
            define_t0 = time() if timing_info is not None else None
            self.define_internal_coordinates()
            if timing_info is not None:
                timing_info['define_internal_coordinates'] = (
                    timing_info.get('define_internal_coordinates', 0.0) +
                    (time() - define_t0))

        # The B matrix is required for the transformation of the
        # gradient, while B2 is required for the transformation of
        # the Hessian from Cartesian to internal coordinates and
        # vice-versa.

        b_t0 = time() if timing_info is not None else None
        self.calculate_b_matrix()
        if timing_info is not None:
            timing_info['calculate_b_matrix'] = (
                timing_info.get('calculate_b_matrix', 0.0) +
                (time() - b_t0))

        int_t0 = time() if timing_info is not None else None
        self.compute_internal_coordinates_values()
        if timing_info is not None:
            timing_info['compute_internal_coordinates_values'] = (
                timing_info.get('compute_internal_coordinates_values', 0.0) +
                (time() - int_t0))

        # The gradient and Hessian are reset to None.
        # The gradient of the new configuration will be
        # calculated by interpolation.
        self.internal_gradient = None
        self.internal_hessian = None

    def reset_coordinates(self, cartesian_coordinates):
        """
        Resets the Cartesian coordinates.

        :param cartesian_coordinates:
                a numpy array of the new Cartesian coordinates.
        """

        self.cartesian_coordinates = cartesian_coordinates

        # Since the Cartesian coordinates have been reset,
        # the internal coordinates and transformation matrices
        # must be redifined.
        self.internal_coordinates_values = None
        self.b_matrix = None
        self.b2_matrix = None
        self.g_minus = None

        if self.internal_coordinates is None:
            self.define_internal_coordinates()
        # The B matrix is required for the transformation of the
        # gradient, while B2 is required for the transformation of
        # the Hessian from Cartesian to internal coordinates and
        # vice-versa.
        self.calculate_b_matrix()
        self.calculate_b2_matrix()
        self.compute_internal_coordinates_values()

        # The gradient and Hessian are reset to None.
        # The gradient of the new configuration will be
        # calculated by interpolation.
        self.internal_gradient = None
        self.internal_hessian = None

    def _compute_internal_coordinates_values_vectorized(self):
        """
        Vectorized internal-coordinate value evaluation.
        """

        assert_msg_critical(
            self.internal_coordinates is not None,
            'InterpolationDatapoint: No internal coordinates are defined.')

        assert_msg_critical(
            self.cartesian_coordinates is not None,
            'InterpolationDatapoint: No cartesian coordinates are defined.')

        n_atoms = self.cartesian_coordinates.shape[0]
        coords = self.cartesian_coordinates.reshape((n_atoms * 3))

        n_rows = len(self.internal_coordinates)
        base_values = np.fromiter(
            (q.value(coords) for q in self.internal_coordinates),
            dtype=np.float64,
            count=n_rows)

        int_coords = base_values.copy()
        row_cache = self._get_b_matrix_row_cache()

        bond_rows = row_cache['bond_rows']
        angle_rows = row_cache["angle_rows"]
        bond_mode = self.get_bond_coordinate_mode()
        angle_mode = self.get_angle_coordinate_mode()

        if bond_rows.size > 0:
            if bond_mode != "r":
                eq_values = None
                if self._bond_mode_requires_eq_reference(bond_mode):
                    eq_values = self._get_eq_bond_lengths_array(bond_rows.size)
                q_values, _, _ = self._bond_transform(
                    base_values[bond_rows],
                    eq_values)
                int_coords[bond_rows] = q_values

        if angle_rows.size > 0:
            if angle_mode != "raw":
                eq_angle_values = self._get_eq_angle_values_array(angle_rows.size)
                q_values, _, _ = self._angle_transform(
                    base_values[angle_rows],
                    eq_angle_values)
                int_coords[angle_rows] = q_values

        self.internal_coordinates_values = int_coords

    def get_bond_permutation_from_mask(self, mask):
        """Converts a full internal-coordinate mask to a bond permutation."""

        assert_msg_critical(
            self.z_matrix is not None,
            'InterpolationDatapoint: No Z-matrix defined.')
        mask_arr = np.asarray(mask, dtype=np.int64).reshape(-1)
        n_ic = len(self.z_matrix)
        assert_msg_critical(
            mask_arr.size == n_ic,
            'InterpolationDatapoint: mask size mismatch '
            f'(expected {n_ic}, got {mask_arr.size}).')

        bond_rows = np.asarray(
            self._get_b_matrix_row_cache()['bond_rows'], dtype=np.int64)
        row_to_bond = {
            int(row): index for index, row in enumerate(bond_rows.tolist())
        }
        permutation = np.empty(bond_rows.size, dtype=np.int64)

        for target_bond, target_row in enumerate(bond_rows):
            source_row = int(mask_arr[int(target_row)])
            assert_msg_critical(
                source_row in row_to_bond,
                'InterpolationDatapoint: mask maps bond row '
                f'{int(target_row)} to non-bond row {source_row}.')
            permutation[target_bond] = row_to_bond[source_row]

        return permutation

    def get_masked_eq_bond_lengths(self, reference_eq_bond_lengths, mask):
        """Applies an internal-coordinate mask to equilibrium bond lengths."""

        bond_rows = np.asarray(
            self._get_b_matrix_row_cache()['bond_rows'], dtype=np.int64)
        eq_reference = np.asarray(
            reference_eq_bond_lengths, dtype=np.float64).reshape(-1)
        assert_msg_critical(
            eq_reference.size == bond_rows.size,
            'InterpolationDatapoint: eq_bond_lengths size mismatch in masked '
            f'mapping (expected {bond_rows.size}, got {eq_reference.size}).')
        return eq_reference[self.get_bond_permutation_from_mask(mask)]

    def symmetrize_eq_bond_lengths_from_masks(self, masks):
        """Averages equilibrium lengths over bond permutation orbits."""

        if (
            not self._bond_mode_requires_eq_reference()
            or self.eq_bond_lengths is None
        ):
            return self.eq_bond_lengths

        masks_array = np.asarray(masks, dtype=np.int64)
        if masks_array.ndim == 1:
            masks_array = masks_array.reshape(1, -1)
        if masks_array.shape[0] == 0:
            return self.eq_bond_lengths

        bond_rows = np.asarray(
            self._get_b_matrix_row_cache()['bond_rows'], dtype=np.int64)
        n_bonds = int(bond_rows.size)
        eq_values = np.asarray(
            self.eq_bond_lengths, dtype=np.float64).reshape(-1).copy()
        assert_msg_critical(
            eq_values.size == n_bonds,
            'InterpolationDatapoint: eq_bond_lengths size mismatch '
            f'(expected {n_bonds}, got {eq_values.size}).')

        parent = np.arange(n_bonds, dtype=np.int64)

        def find(index):
            while parent[index] != index:
                parent[index] = parent[parent[index]]
                index = parent[index]
            return index

        def union(left, right):
            left_root = find(left)
            right_root = find(right)
            if left_root != right_root:
                parent[right_root] = left_root

        for mask in masks_array:
            permutation = self.get_bond_permutation_from_mask(mask)
            for index, mapped_index in enumerate(permutation):
                union(index, int(mapped_index))

        groups = {}
        for index in range(n_bonds):
            groups.setdefault(find(index), []).append(index)
        for indices in groups.values():
            eq_values[indices] = float(np.mean(eq_values[indices]))

        self.eq_bond_lengths = eq_values
        return self.eq_bond_lengths

    def prepare_eq_candidate_base_cache(self):
        """Caches geometry-only data for candidate-local eq transformations."""

        assert_msg_critical(
            self.use_eq_bond_length,
            'InterpolationDatapoint: eq candidate cache requires '
            'use_eq_bond_length=True.')
        assert_msg_critical(
            self.cartesian_coordinates is not None,
            'InterpolationDatapoint: No Cartesian coordinates defined.')
        assert_msg_critical(
            self.internal_coordinates is not None,
            'InterpolationDatapoint: Internal coordinates not defined.')
        assert_msg_critical(
            self.b_matrix is not None and self.original_b_matrix is not None,
            'InterpolationDatapoint: B matrices are not available.')

        bond_rows = np.asarray(
            self._get_b_matrix_row_cache()['bond_rows'], dtype=np.int64)
        coordinates = self.cartesian_coordinates.reshape(-1)
        raw_bond_values = np.array([
            self.internal_coordinates[row].value(coordinates)
            for row in bond_rows
        ], dtype=np.float64)

        return {
            'bond_rows': bond_rows,
            'base_internal': np.asarray(
                self.internal_coordinates_values, dtype=np.float64).copy(),
            'base_b_matrix': np.asarray(
                self.b_matrix, dtype=np.float64).copy(),
            'base_original_b_matrix': np.asarray(
                self.original_b_matrix, dtype=np.float64).copy(),
            'inv_sqrt_masses': (
                None if self.inv_sqrt_masses is None else np.asarray(
                    self.inv_sqrt_masses, dtype=np.float64).copy()),
            'raw_bond_values': raw_bond_values,
        }

    def build_current_chart_with_eq_fast(
            self, eq_bond_lengths_candidate, base_cache):
        """Builds current internal coordinates and B matrix for one eq vector."""

        bond_rows = np.asarray(base_cache['bond_rows'], dtype=np.int64)
        eq_candidate = np.asarray(
            eq_bond_lengths_candidate, dtype=np.float64).reshape(-1)
        assert_msg_critical(
            eq_candidate.size == bond_rows.size,
            'InterpolationDatapoint: candidate eq size mismatch '
            f'(expected {bond_rows.size}, got {eq_candidate.size}).')

        bond_values, derivatives, _ = self._bond_transform(
            base_cache['raw_bond_values'], eq_candidate)
        internal_values = base_cache['base_internal'].copy()
        internal_values[bond_rows] = bond_values

        b_matrix = base_cache['base_b_matrix'].copy()
        bond_block = (
            base_cache['base_original_b_matrix'][bond_rows, :]
            * derivatives[:, None]
        )
        if base_cache['inv_sqrt_masses'] is not None:
            bond_block *= base_cache['inv_sqrt_masses'][None, :]
        b_matrix[bond_rows, :] = bond_block
        return internal_values, b_matrix

    def build_masked_current_chart_from_reference_eq_fast(
            self, reference_eq_bond_lengths, mask, base_cache):
        eq_masked = self.get_masked_eq_bond_lengths(
            reference_eq_bond_lengths, mask)
        return self.build_current_chart_with_eq_fast(eq_masked, base_cache)

    def compute_internal_coordinates_values(self):
        """
        Creates an array with the values of the internal coordinates
        and saves it in self.
        """
        self._compute_internal_coordinates_values_vectorized()

    def get_z_matrix_as_np_arrays(self):

        zmat = self.z_matrix_dict
        return {
            "bonds": np.array(zmat["bonds"], dtype=np.int64),
            "angles": np.array(zmat["angles"], dtype=np.int64),
            "dihedrals": np.array(zmat["dihedrals"], dtype=np.int64),
            "impropers": np.array(zmat["impropers"], dtype=np.int64),
        }

    def get_imp_int_coord_as_np_arrays(self):
        """
        Returns a dictionary with the numpy arrays corresponding to the bonds,
        bond angles, and dihedral angles defined by the Z-matrix.
        """
        assert_msg_critical(self.imp_int_coordinates is not None,
                            'InterpolationDatapoint: No Z-matrix is defined.')

        imp_coords = self.imp_int_coordinates
        return {
            "imp_bonds": np.array(imp_coords["bonds"], dtype=np.int64),
            "imp_angles": np.array(imp_coords["angles"], dtype=np.int64),
            "imp_dihedrals": np.array(imp_coords["dihedrals"], dtype=np.int64),
            "imp_impropers": np.array(imp_coords["impropers"], dtype=np.int64),
        }

    def calculate_translation_coordinates(self):
        """Center the molecule by translating its geometric center to (0, 0, 0)."""
        center = np.mean(self.cartesian_coordinates, axis=0)
        translated_coordinates = self.cartesian_coordinates - center

        return translated_coordinates

    def _write_string_dataset(self, h5f, name, value):
        if value is None:
            return
        dt = h5py.string_dtype(encoding="utf-8")
        h5f.create_dataset(name, data=np.array(value, dtype=object), dtype=dt)

    def _write_optional_array(self, h5f, name, value):
        if value is None:
            return
        array_value = np.asarray(value)
        if array_value.ndim == 0 or array_value.size == 0:
            h5f.create_dataset(name, data=array_value)
        else:
            h5f.create_dataset(name, data=array_value, compression="gzip")

    def _write_optional_string_sequence(self, h5f, name, value):
        if value is None:
            return
        if isinstance(value, str):
            self._write_string_dataset(h5f, name, value)
            return
        dt = h5py.string_dtype(encoding="utf-8")
        h5f.create_dataset(
            name,
            data=np.array([str(item) for item in value], dtype=object),
            dtype=dt,
        )

    def _metadata_json_safe(self, value):
        if isinstance(value, dict):
            return {
                str(key): self._metadata_json_safe(val)
                for key, val in value.items()
            }
        if isinstance(value, (list, tuple)):
            return [self._metadata_json_safe(item) for item in value]
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        return value

    def write_hdf5(self, fname, label):
        """
        Writes the energy, internal coordinates, internal gradient, and
        internal Hessian to a checkpoint file.

        :param fname:
            Name of the checkpoint file to be written.
        :param label:
            A string describing the coordinate.
        :param write_zmat:
            flag to determine if the Z-matrix should be writen to the
            checkpoint file. False by default.
        """
        valid_checkpoint = (fname and isinstance(fname, str))

        if valid_checkpoint:
            try:
                h5f = h5py.File(fname, 'a')
            except IOError:
                h5f = h5py.File(fname, 'w')

            label += self.coordinate_label_suffix()
            label += "_dihedral"

            assert_msg_critical(self.energy is not None,
                                'InterpolationDatapoint: No energy is defined.')

            assert_msg_critical(
                self.internal_gradient is not None,
                'InterpolationDatapoint: No internal gradient is defined.')

            assert_msg_critical(
                self.internal_hessian is not None,
                'InterpolationDatapoint: No internal Hessian is defined.')

            if self.atom_labels is not None:
                full_label = label + '_atom_labels'
                string_type = h5py.string_dtype(encoding='utf-8')
                h5f.create_dataset(full_label, data=self.atom_labels, dtype=string_type)

            full_label = label + "_energy"
            h5f.create_dataset(full_label, data=np.float64(self.energy))

            # write gradient in internal coordinates
            full_label = label + "_gradient"
            h5f.create_dataset(full_label,
                               data=self.internal_gradient,
                               compression='gzip')

            full_label = label + "_cartesian_gradient"
            h5f.create_dataset(full_label,
                               data=self.gradient,
                               compression='gzip')
            # write Hessian in internal coordinates
            full_label = label + "_hessian"
            h5f.create_dataset(full_label,
                               data=self.internal_hessian,
                               compression='gzip')

            full_label = label + "_cartesian_hessian"
            h5f.create_dataset(full_label,
                               data=self.hessian,
                               compression='gzip')

            if self.internal_coordinates is None:
                assert_msg_critical(
                    self.internal_coordinates_values is not None,
                    'InterpolationDatapoint: No internal coordinates are defined.')

                # write internal coordinates
                full_label = label + "_internal_coordinates"
                h5f.create_dataset(full_label,
                                   data=self.internal_coordinates_values,
                                   compression='gzip')
            else:
                # write internal coordinates
                full_label = label + "_internal_coordinates"
                self.compute_internal_coordinates_values()
                h5f.create_dataset(full_label,
                                   data=self.internal_coordinates_values,
                                   compression='gzip')

            if not self.use_eq_bond_length and self.eq_bond_lengths is None:
                assert_msg_critical(
                    self.eq_bond_lengths is not None,
                    'InterpolationDatapoint: No eq bond lenths are defined.')

            else:
                # write internal coordinates
                full_label = label + "_eq_bond_lengths"
                h5f.create_dataset(full_label,
                                   data=self.eq_bond_lengths,
                                   compression='gzip')

            if self.eq_angle_values is not None:
                full_label = label + "_eq_angle_values"
                h5f.create_dataset(
                    full_label,
                    data=np.asarray(self.eq_angle_values, dtype=np.float64),
                    compression='gzip',
                )

            if self.inv_sqrt_masses is not None:
                full_label = label + "_inv_sqrt_masses"
                h5f.create_dataset(
                    full_label,
                    data=np.asarray(self.inv_sqrt_masses, dtype=np.float64),
                    compression='gzip',
                )

            assert_msg_critical(
                self.confidence_radius is not None,
                'InterpolationDatapoint: No confidence radius is defined.')
            # write internal coordinates
            full_label = label + "_confidence_radius"
            h5f.create_dataset(full_label, data=np.float64(self.confidence_radius))

            assert_msg_critical(
                self.cartesian_coordinates is not None,
                'InterpolationDatapoint: No Cartesian coordinates are defined.')
            # write Cartesian coordinates
            full_label = label + "_cartesian_coordinates"
            h5f.create_dataset(full_label,
                               data=self.cartesian_coordinates,
                               compression='gzip')

            assert_msg_critical(
                self.z_matrix is not None,
                'InterpolationDatapoint: No Z-matrix is defined.')
            # Write the bonds, angles and dihedrals defined in the Z-matrix.
            # Bonds are saved in an array with two atoms indices per row.
            # Angles are saved in an array with three atom indices per row.
            # Dihedrals are saved in an array with four atom indices per row.
            zmat_dict = self.get_z_matrix_as_np_arrays()

            for key in zmat_dict.keys():
                h5f.create_dataset(label + '_' + key,
                                   data=zmat_dict[key],
                                   compression='gzip')

            if self.identify_imp_int_coord:
                assert_msg_critical(
                    self.imp_int_coordinates is not None,
                    'InterpolationDatapoint: No important internal coordinates are defined.')
                # Write the bonds, angles and dihedrals defined in the Z-matrix.
                # Bonds are saved in an array with two atoms indices per row.
                # Angles are saved in an array with three atom indices per row.
                # Dihedrals are saved in an array with four atom indices per row.
                imp_int_coord_dict = self.get_imp_int_coord_as_np_arrays()
                for key in imp_int_coord_dict.keys():
                    h5f.create_dataset(label + '_' + key,
                                       data=imp_int_coord_dict[key],
                                       compression='gzip')

            write_local_metadata = (
                self.mapping_masks is not None
                or self.family_label is not None
                or self.bank_role != "global"
                or self.cluster_id is not None
                or self.active_atoms is not None
                or self.active_rows is not None
                or bool(self.canonical_subset_key)
                or bool(self.response_rows)
                or self.relaxation_policy_id is not None
                or self.is_anchor
                or bool(self.anchor_state_ids)
                or self.local_factor_combination_mode != "signed_full"
                or self.local_factor_overlap_source is not None
            )

            if write_local_metadata:
                self._write_optional_array(h5f, label + "_masks", self.mapping_masks)
                self._write_string_dataset(h5f, label + "_family_label", self.family_label)
                self._write_string_dataset(h5f, label + "_bank_role", self.bank_role)
                self._write_string_dataset(h5f, label + "_cluster_id", self.cluster_id)
                self._write_string_dataset(h5f, label + "_cluster_type", self.cluster_type)
                self._write_optional_string_sequence(
                    h5f, label + "_cluster_rotor_ids", self.cluster_rotor_ids)
                self._write_optional_array(
                    h5f, label + "_cluster_state_id", [self.cluster_state_id])
                self._write_optional_array(
                    h5f, label + "_dihedrals_to_rotate", self.dihedrals_to_rotate)
                self._write_optional_array(
                    h5f, label + "_phase_signature", self.phase_signature)
                self._write_optional_array(
                    h5f, label + "_active_atoms", self.active_atoms)
                self._write_optional_array(
                    h5f, label + "_active_rows", self.active_rows)
                self._write_optional_string_sequence(
                    h5f,
                    label + "_canonical_subset_key",
                    self.canonical_subset_key,
                )
                self._write_optional_array(
                    h5f, label + "_response_rows", self.response_rows)
                self._write_string_dataset(
                    h5f,
                    label + "_relaxation_policy_id",
                    self.relaxation_policy_id,
                )
                self._write_optional_array(
                    h5f, label + "_is_anchor", [int(bool(self.is_anchor))])
                self._write_optional_array(
                    h5f, label + "_anchor_state_ids", self.anchor_state_ids)
                self._write_string_dataset(
                    h5f,
                    label + "_local_factor_combination_mode",
                    self.local_factor_combination_mode,
                )
                self._write_string_dataset(
                    h5f,
                    label + "_local_factor_overlap_source",
                    self.local_factor_overlap_source,
                )

            if self.metadata:
                metadata_json = json.dumps(
                    self._metadata_json_safe(self.metadata),
                    sort_keys=True,
                )
                self._write_string_dataset(
                    h5f, label + "_metadata_json", metadata_json)

            h5f.close()

    def read_hdf5(self, fname, label):
        """
        Reads the energy, internal coordinates, gradient, and Hessian from
        a checkpoint file,

        :param fname:
            Name of the checkpoint file.
            The file must exist before calling this routine.
        :param label:
            The label of the selected coordinates.
        """

        def _read_scalar_string(ds):
            if ds is None:
                return None
            val = ds[()]
            if isinstance(val, bytes):
                return val.decode("utf-8")
            return str(val)

        def _read_optional_array(h5f, name):
            ds = h5f.get(name)
            if ds is None:
                return None
            return np.array(ds)

        def _read_optional_string_sequence(h5f, name):
            ds = h5f.get(name)
            if ds is None:
                return None
            val = ds[()]
            if isinstance(val, bytes):
                return (val.decode("utf-8"),)
            arr = np.asarray(val).reshape(-1)
            out = []
            for item in arr:
                if isinstance(item, bytes):
                    out.append(item.decode("utf-8"))
                else:
                    out.append(str(item))
            return tuple(out)

        valid_checkpoint = (fname and isinstance(fname, str))

        if valid_checkpoint:
            h5f = h5py.File(fname, 'r')

            label += self.coordinate_label_suffix()
            label += "_dihedral"

            energy_label = label + "_energy"
            gradient_label = label + "_gradient"
            hessian_label = label + "_hessian"
            cart_hess_label = label + "_cartesian_hessian"
            cart_grad_label = label + "_cartesian_gradient"
            coords_label = label + "_internal_coordinates"
            eq_bond_length_label = label + "_eq_bond_lengths"
            eq_angle_values_label = label + "_eq_angle_values"
            inv_sqrt_masses_label = label + "_inv_sqrt_masses"
            cart_coords_label = label + "_cartesian_coordinates"
            confidence_radius_label = label + "_confidence_radius"
            mapping_masks_label = label + "_masks"
            metadata_json_label = label + "_metadata_json"

            z_matrix_bonds = label + '_bonds'
            z_matrix_angles = label + '_angles'
            z_matrix_dihedrals = label + '_dihedrals'
            z_matrix_impropers = label + '_impropers'

            imp_int_coord_bonds = label + '_imp_bonds'
            imp_int_coord_angles = label + '_imp_angles'
            imp_int_coord_dihedrals = label + '_imp_dihedrals'
            imp_int_coord_impropers = label + '_imp_impropers'

            self.point_label = label
            self.energy = np.array(h5f.get(energy_label))
            self.internal_gradient = np.array(h5f.get(gradient_label))
            self.gradient = np.array(h5f.get(cart_grad_label))
            self.internal_hessian = np.array(h5f.get(hessian_label))
            self.hessian = np.array(h5f.get(cart_hess_label))
            self.internal_coordinates_values = np.array(h5f.get(coords_label))
            self.eq_bond_lengths = np.array(h5f.get(eq_bond_length_label))
            self.eq_angle_values = _read_optional_array(
                h5f, eq_angle_values_label)
            self.inv_sqrt_masses = _read_optional_array(
                h5f, inv_sqrt_masses_label)
            self.cartesian_coordinates = np.array(h5f.get(cart_coords_label))
            self.confidence_radius = np.array(h5f.get(confidence_radius_label))
            self.mapping_masks = _read_optional_array(h5f, mapping_masks_label)
            metadata_json = _read_scalar_string(h5f.get(metadata_json_label))
            self.metadata = json.loads(metadata_json) if metadata_json else {}
            self.family_label = _read_scalar_string(
                h5f.get(label + "_family_label"))
            self.bank_role = (
                _read_scalar_string(h5f.get(label + "_bank_role")) or "global"
            )
            self.cluster_id = _read_scalar_string(
                h5f.get(label + "_cluster_id"))
            self.cluster_type = _read_scalar_string(
                h5f.get(label + "_cluster_type"))
            self.cluster_rotor_ids = _read_optional_string_sequence(
                h5f, label + "_cluster_rotor_ids")
            cluster_state_id = _read_optional_array(
                h5f, label + "_cluster_state_id")
            if cluster_state_id is not None and cluster_state_id.size > 0:
                self.cluster_state_id = int(cluster_state_id.reshape(-1)[0])
            self.dihedrals_to_rotate = _read_optional_array(
                h5f, label + "_dihedrals_to_rotate")
            self.phase_signature = _read_optional_array(
                h5f, label + "_phase_signature")
            self.active_atoms = _read_optional_array(
                h5f, label + "_active_atoms")
            self.active_rows = _read_optional_array(
                h5f, label + "_active_rows")
            self.canonical_subset_key = (
                _read_optional_string_sequence(
                    h5f, label + "_canonical_subset_key") or ()
            )
            response_rows = _read_optional_array(
                h5f, label + "_response_rows")
            self.response_rows = (
                () if response_rows is None else tuple(
                    int(row) for row in response_rows.reshape(-1))
            )
            self.relaxation_policy_id = _read_scalar_string(
                h5f.get(label + "_relaxation_policy_id"))
            is_anchor = _read_optional_array(h5f, label + "_is_anchor")
            self.is_anchor = bool(
                is_anchor.reshape(-1)[0]) if is_anchor is not None else False
            anchor_state_ids = _read_optional_array(
                h5f, label + "_anchor_state_ids")
            self.anchor_state_ids = (
                () if anchor_state_ids is None else tuple(
                    int(state_id)
                    for state_id in anchor_state_ids.reshape(-1))
            )
            self.local_factor_combination_mode = (
                _read_scalar_string(
                    h5f.get(label + "_local_factor_combination_mode"))
                or "signed_full"
            )
            self.local_factor_overlap_source = _read_scalar_string(
                h5f.get(label + "_local_factor_overlap_source"))

            z_matrix_dict = {}
            for kname, key in [
                (z_matrix_bonds, "bonds"),
                (z_matrix_angles, "angles"),
                (z_matrix_dihedrals, "dihedrals"),
                (z_matrix_impropers, "impropers")
            ]:
                ds = h5f.get(kname)
                if ds is not None:
                    z_matrix_dict[key] = [tuple(x.tolist()) for x in ds]
            self.z_matrix = self.flatten_z_matrix(z_matrix_dict)
            self._b_matrix_row_cache = None

            if self.identify_imp_int_coord:
                self.imp_int_coordinates = {}
                for kname, key in [
                    (imp_int_coord_bonds, "bonds"),
                    (imp_int_coord_angles, "angles"),
                    (imp_int_coord_dihedrals, "dihedrals"),
                    (imp_int_coord_impropers, "impropers")
                ]:
                    ds = h5f.get(kname)
                    if ds is not None:
                        self.imp_int_coordinates[key] = [tuple(x.tolist()) for x in ds]

            h5f.close()

    def cartesian_distance_vector(self, data_point):
        """Calculates and returns the cartesian distance between
           self and data_point.

           :param data_point:
                InterpolationDatapoint object
        """
        # First, translate the cartesian coordinates to zero
        target_coordinates = data_point.calculate_translation_coordinates()
        reference_coordinates = self.calculate_translation_coordinates()

        # Then, determine the rotation matrix which
        # aligns data_point (target_coordinates)
        # to self (reference_coordinates)
        rotation_matrix = geometric.rotate.get_rot(target_coordinates,
                                                   reference_coordinates)

        # Rotate the data point:
        rotated_coordinates = np.dot(rotation_matrix, target_coordinates.T).T

        # Calculate the distance:
        distance_vector = (reference_coordinates - rotated_coordinates)

        return distance_vector

    def remove_point_from_hdf5(self, fname, label, use_inverse_bond_length=False, use_eq_bond_length=False):
        """
        Removes a point (i.e., corresponding datasets) from an HDF5 file based on label and internal settings.

        :param fname: HDF5 filename
        :param label: base label for the data
        :param use_inverse_bond_length: whether to append '_rinv' or '_r'
        """
        if use_inverse_bond_length:
            label += "_rinv"
        if use_eq_bond_length:
            label += "_eq"
        else:
            label += "_r"

        label += "_dihedral"

        keys_to_remove = [
            label + "_energy",
            label + "_gradient",
            label + "_hessian",
            label + "_internal_coordinates",
            label + "eq_bond_lengths",
            label + "_cartesian_coordinates",
            label + "_cartesian_gradient",
            label + "_cartesian_hessian",
            label + "_confidence_radius",
            label + "_bonds",
            label + "_angles",
            label + "_dihedrals",
            label + "_impropers",
            label + "_imp_bonds",
            label + "_imp_angles",
            label + "_imp_dihedrals",
            label + "_imp_impropers"
        ]

        with h5py.File(fname, 'r+') as h5f:
            for key in keys_to_remove:
                if key in h5f:
                    del h5f[key]
                    print(f"Deleted: {key}")
                else:
                    print(f"Key not found (skipped): {key}")

    def update_confidence_radius(self, fname, label, new_confidence_radius):
        """
        Updates the confidence radius in the HDF5 file.

        :param fname: HDF5 filename
        :param label: base label for dataset lookup
        :param new_confidence_radius: new value(s) to store
        """
        with h5py.File(fname, 'r+') as h5f:
            label += self.coordinate_label_suffix()
            label += "_dihedral"

            confidence_radius_label = label + "_confidence_radius"

            if confidence_radius_label in h5f:
                dataset = h5f[confidence_radius_label]

                # Check if shape matches before replacing
                if np.shape(dataset) == np.shape(new_confidence_radius):
                    dataset[...] = new_confidence_radius
                else:
                    raise ValueError(f"Shape mismatch: existing {dataset.shape}, new {np.shape(new_confidence_radius)}")
            else:
                raise KeyError(f"{confidence_radius_label} not found in file.")
