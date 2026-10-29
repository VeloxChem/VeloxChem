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
import os
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
    Implements a datapoint of the interpolated potential energy surface: a
    geometry with its energy, gradient and Hessian transformed to the
    internal coordinates of a z-matrix, and its storage in the HDF5
    database. InterpolationDriver uses the same class for the current
    geometry (impes_coordinate), where only the internal coordinates and
    the B-matrix are needed.

    Internal coordinates, one per z-matrix entry in the flat order bonds,
    angles, proper dihedrals, impropers:
        - bonds: 1/r in 1/bohr (use_inverse_bond_length, the default),
          r_eq - r_eq^2/r in bohr (use_eq_bond_length, with eq_bond_lengths
          r_eq in bohr), or r in bohr (both flags False). A settings
          dictionary that gives only one of the flags as True sets the other
          one to False; both True raises.
        - angles, proper and improper dihedrals: radians (geomeTRIC Angle and
          Dihedral; there are no cosine dihedrals).

    Metric: if inv_sqrt_masses (1/sqrt(m), one entry per Cartesian
    component) is set, b_matrix and b2_matrix are derivatives with respect
    to the mass-weighted coordinates Q = sqrt(m) x, and gradient and
    hessian must be given in the same metric (m^-1/2 g and
    m^-1/2 H m^-1/2), as IMForceFieldGenerator does for use_mass_weight=True.

    :param z_matrix:
        The z-matrix dictionary {'bonds', 'angles', 'dihedrals',
        'impropers'} of 0-based atom-index tuples, or None (read_hdf5 then
        takes the z-matrix from the file).
    :param atom_labels:
        Optional atom labels; written to and read from the HDF5 file.
    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - z_matrix_dict: the z-matrix dictionary; z_matrix: its flat list.
        - point_label: after read_hdf5 the dataset prefix of the point,
          including the bond-mode suffix (e.g. 'point_1_eq_dihedral').
        - energy: the energy in Hartree.
        - gradient, hessian: the Cartesian gradient (natoms, 3) in
          Hartree/bohr and Hessian (3N, 3N) in Hartree/bohr**2, in the build
          metric (mass-weighted if inv_sqrt_masses was set; the HDF5
          datasets record it in the attribute 'mass_weighted').
        - internal_gradient, internal_hessian: the gradient and Hessian in
          the internal coordinates.
        - internal_coordinates: the geomeTRIC primitives;
          internal_coordinates_values: their values (units above).
        - cartesian_coordinates: the geometry in bohr, (natoms, 3).
        - b_matrix: the Wilson B-matrix dq/dx (or dq/dQ); b2_matrix: the
          second derivatives of q; original_b_matrix: the unscaled rows of
          the last B-matrix (dr/dx for bonds; no longer used by this class,
          kept for subclasses).
        - inv_sqrt_masses: 1/sqrt(m) per Cartesian component, or None.
        - use_inverse_bond_length, use_eq_bond_length: the bond mode.
        - eq_bond_lengths: r_eq per bond in bohr (eq mode).
        - confidence_radius: the confidence radius R in bohr.
        - imp_int_coordinates: the important coordinates, a dictionary
          {'bonds', 'angles', 'dihedrals', 'impropers'} of 0-based tuples in
          z-matrix orientation, used by the target-customised weights.
        - identify_imp_int_coord: write and read imp_int_coordinates.
        - atom_labels: the atom labels, or None.
    """

    def __init__(self, z_matrix=None, atom_labels=None, comm=None, ostream=None):
        """
        Initializes the InterpolationDatapoint object.

        :param z_matrix:
            The z-matrix dictionary {'bonds', 'angles', 'dihedrals',
            'impropers'} of 0-based atom-index tuples, or None.
        :param atom_labels:
            Optional atom labels.
        :param comm:
            The MPI communicator (default MPI.COMM_WORLD).
        :param ostream:
            The output stream (default: stdout on the master rank).
        """
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.ostream = ostream

        self.point_label = None

        self.energy = None
        self.gradient = None
        self.hessian = None
        self.z_matrix_dict = z_matrix
        if z_matrix is not None:
            z_matrix_flat = self.flatten_z_matrix(z_matrix)
            self.z_matrix = z_matrix_flat
        else:
            self.z_matrix = None
        self.atom_labels = atom_labels
        self.internal_gradient = None
        self.internal_hessian = None
        self.inv_sqrt_masses = None

        self.b_matrix = None

        # unscaled rows of the last B-matrix (dr/dx for bonds); not used by
        # this class any more, kept for subclasses
        self.original_b_matrix = None
        self.b2_matrix = None

        self.confidence_radius = None
        self.use_inverse_bond_length = True
        self.use_eq_bond_length = False
        self.identify_imp_int_coord = True

        self.eq_bond_lengths = None

        # internal_coordinates is a list of geomeTRIC objects which represent
        # different types of internal coordinates (distances, angles, dihedrals)
        self.internal_coordinates = None
        self.imp_int_coordinates = {
            'bonds': [],
            'angles': [],
            'dihedrals': [],
            'impropers': [],
        }

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
                     'use the bond coordinate r_eq - r_eq^2/r'),
            }
        }

    @staticmethod
    def flatten_z_matrix(zm_dict):
        """
        Returns the flat z-matrix list (bonds, angles, dihedrals, impropers)
        of a z-matrix dictionary.
        """
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

        if self.use_eq_bond_length and not self.use_inverse_bond_length:
            self.use_inverse_bond_length = False
        elif self.use_inverse_bond_length and not self.use_eq_bond_length:
            self.use_eq_bond_length = False
        elif self.use_eq_bond_length and self.use_inverse_bond_length:
            assert_msg_critical(
                False, 'use_inverse_bond_length and use_eq_bond_length '
                'cannot both be True. Set exactly one of them (or neither, for '
                'plain bond lengths in bohr).')

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
        Returns the row metadata (the bond rows) of the z-matrix, cached
        until the z-matrix changes.
        """

        # The cache is only valid for the z-matrix it was built from.
        z_matrix_key = tuple(tuple(int(i) for i in z) for z in self.z_matrix)
        if (self._b_matrix_row_cache is not None and
                self._b_matrix_row_cache['z_matrix'] == z_matrix_key):
            return self._b_matrix_row_cache

        n_rows = len(self.z_matrix)
        row_sizes = np.fromiter((len(z) for z in self.z_matrix),
                                dtype=np.int8,
                                count=n_rows)
        bond_rows = np.flatnonzero(row_sizes == 2)

        self._b_matrix_row_cache = {
            'z_matrix': z_matrix_key,
            'bond_rows': bond_rows,
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

        use_original = (
            self.use_inverse_bond_length or
            self.use_eq_bond_length
        )
        if use_original:
            self.original_b_matrix = derivatives.copy()
        else:
            self.original_b_matrix = None

        row_scale = np.ones(n_rows, dtype=np.float64)
        row_cache = self._get_b_matrix_row_cache()

        bond_rows = row_cache['bond_rows']
        if bond_rows.size > 0 and (self.use_inverse_bond_length or self.use_eq_bond_length):
            bond_values = np.array(
                [self.internal_coordinates[idx].value(coords) for idx in bond_rows],
                dtype=np.float64)

            if self.use_inverse_bond_length:
                row_scale[bond_rows] = -1.0 / np.square(bond_values)
            elif self.use_eq_bond_length:
                eq_values = self._get_eq_bond_lengths_array(bond_rows.size)
                _, dq_dr, _ = self._bond_transform(
                    bond_values, eq_values)
                row_scale[bond_rows] = dq_dr

        self.b_matrix = derivatives * row_scale[:, np.newaxis]

        if self.inv_sqrt_masses is not None:
            self.b_matrix = self.b_matrix * self.inv_sqrt_masses

    def calculate_b_matrix(self):
        """
        Calculates the Wilson B-matrix at cartesian_coordinates (bond rows
        scaled by dq/dr of the bond mode; mass-weighted if inv_sqrt_masses
        is set).
        """
        self._calculate_b_matrix_vectorized()

    def calculate_b2_matrix(self):
        """
        Calculates the second derivatives of the internal coordinates with
        respect to the Cartesian (or mass-weighted) coordinates at
        cartesian_coordinates.
        """

        assert_msg_critical(self.internal_coordinates is not None, 'InterpolationDatapoint: No internal coordinates are defined.')
        assert_msg_critical(self.cartesian_coordinates is not None, 'InterpolationDatapoint: No cartesian coordinates are defined.')

        n_atoms = self.cartesian_coordinates.shape[0]
        coords = self.cartesian_coordinates.reshape((n_atoms * 3))

        self.b2_matrix = np.zeros((len(self.z_matrix), n_atoms * 3, n_atoms * 3))

        eq_bond_values = None
        if self.use_eq_bond_length:
            bond_rows = self._get_b_matrix_row_cache()['bond_rows']
            eq_bond_values = self._get_eq_bond_lengths_array(bond_rows.size)

        bond_counter = 0
        for i, z in enumerate(self.z_matrix):
            q = self.internal_coordinates[i]
            second_derivative = q.second_derivative(coords).reshape(-1, n_atoms * 3)

            if len(z) == 2:
                # dr/dx at the current geometry (original_b_matrix may stem
                # from an earlier calculate_b_matrix call at another geometry)
                dr_dx = q.derivative(coords).reshape(-1)
                if self.use_inverse_bond_length:

                    r = q.value(coords)
                    r_inv_2 = 1.0 / (r * r)
                    r_inv_3 = r_inv_2 / r
                    self.b2_matrix[i] = -r_inv_2 * second_derivative

                    for m in range(n_atoms):
                        for n in range(n_atoms):
                            self.b2_matrix[i, m*3:(m+1)*3, n*3:(n+1)*3] += 2 * r_inv_3 * np.outer(dr_dx[m*3:(m+1)*3], dr_dx[n*3:(n+1)*3])

                elif self.use_eq_bond_length:
                    r = q.value(coords)
                    _, dq_dr, d2q_dr2 = self._bond_transform(
                        r, eq_bond_values[bond_counter])
                    self.b2_matrix[i] = dq_dr * second_derivative
                    for m in range(n_atoms):
                        for n in range(n_atoms):
                            self.b2_matrix[i, m*3:(m+1)*3, n*3:(n+1)*3] += d2q_dr2 * np.outer(dr_dx[m*3:(m+1)*3], dr_dx[n*3:(n+1)*3])

                else:
                    self.b2_matrix[i] = second_derivative

                bond_counter += 1

            else:
                self.b2_matrix[i] = second_derivative

        if self.inv_sqrt_masses is not None:
            self.b2_matrix = (self.b2_matrix *
                              self.inv_sqrt_masses.reshape(1, -1, 1) *
                              self.inv_sqrt_masses.reshape(1, 1, -1))

    def transform_gradient_to_internal_coordinates(self, tol=1e-6):
        """
        Transforms the gradient to internal coordinates,
        g_q = G^- B g with G = B B^T and G^- its pseudo-inverse keeping the
        singular values above tol (absolute).

        :param tol:
            Tolerance for the singular values of G.
        """

        assert_msg_critical(self.gradient is not None, 'InterpolationDatapoint: No gradient is defined.')

        if self.b_matrix is None:
            self.calculate_b_matrix()

        g_matrix = np.dot(self.b_matrix, self.b_matrix.T)

        U, s, Vt = np.linalg.svd(g_matrix)

        # Make zero the values of s_inv that are smaller than tol
        s_inv = np.array([1 / s_i if s_i > tol else 0.0 for s_i in s])

        g_minus_matrix = np.dot(U, np.dot(np.diag(s_inv), Vt))

        gradient_flat = self.gradient.flatten()
        self.internal_gradient = np.dot(g_minus_matrix, np.dot(self.b_matrix, gradient_flat))

    def transform_hessian_to_internal_coordinates(self, tol=1e-6):
        """
        Transforms the Hessian to internal coordinates,
        H_q = G^- B (H - sum_k (g_q)_k B2_k) B^T G^-, symmetrised.

        :param tol:
            Tolerance for the singular values of G.
        """

        if self.internal_gradient is None:
            self.transform_gradient_to_internal_coordinates()

        if self.b2_matrix is None:
            self.calculate_b2_matrix()

        g_matrix = np.dot(self.b_matrix, self.b_matrix.T)

        U, s, Vt = np.linalg.svd(g_matrix)

        # Make zero the values of s_inv that are smaller than tol
        s_inv = np.array([1 / s_i if s_i > tol else 0.0 for s_i in s])

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

    def transform_gradient_and_hessian(self):
        """
        Transforms gradient and hessian to internal_gradient and
        internal_hessian and computes the internal coordinate values, at
        cartesian_coordinates.
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
        Sets new Cartesian coordinates for the geometry of an
        InterpolationDriver and recomputes the B-matrix and the internal
        coordinate values. B2 is not needed at runtime and is reset.

        :param cartesian_coordinates:
            The new Cartesian coordinates in bohr, (natoms, 3).
        :param timing_info:
            Optional dictionary; the time spent in
            define_internal_coordinates, calculate_b_matrix and
            compute_internal_coordinates_values is added to its entries.
        """
        self.cartesian_coordinates = cartesian_coordinates

        # Since the Cartesian coordinates have been reset,
        # the internal coordinates and transformation matrices
        # must be redifined.
        self.internal_coordinates_values = None
        self.b_matrix = None
        self.b2_matrix = None

        if self.internal_coordinates is None:
            define_t0 = time() if timing_info is not None else None
            self.define_internal_coordinates()
            if timing_info is not None:
                timing_info['define_internal_coordinates'] = (
                    timing_info.get('define_internal_coordinates', 0.0) +
                    (time() - define_t0))

        # Only B is needed here: the model gradient is B^T dT/dq. B2 is
        # needed only to transform a Hessian (see reset_coordinates).

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
        Sets new Cartesian coordinates and recomputes B, B2 and the internal
        coordinate values (used when a datapoint is built).

        :param cartesian_coordinates:
            The new Cartesian coordinates in bohr, (natoms, 3).
        """

        self.cartesian_coordinates = cartesian_coordinates

        # Since the Cartesian coordinates have been reset,
        # the internal coordinates and transformation matrices
        # must be redifined.
        self.internal_coordinates_values = None
        self.b_matrix = None
        self.b2_matrix = None

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

        if bond_rows.size > 0:
            if self.use_inverse_bond_length:
                int_coords[bond_rows] = 1.0 / base_values[bond_rows]
            elif self.use_eq_bond_length:

                eq_values = self._get_eq_bond_lengths_array(bond_rows.size)
                q_values, _, _ = self._bond_transform(
                    base_values[bond_rows],
                    eq_values)
                int_coords[bond_rows] = q_values

        self.internal_coordinates_values = int_coords

    def _bond_transform(self, r, r_eq):
        """
        Returns q(r), dq/dr, d2q/dr2 for the equilibrium-scaled inverse
        bond coordinate

            q(r) = -r_eq^2 (1/r - 1/r_eq) = r_eq - r_eq^2 / r.

        q is affine in 1/r everywhere (there is no switching function). It
        is normalized such that q(r_eq) = 0 and dq/dr(r_eq) = 1, so that
        q = (r - r_eq) + O((r - r_eq)^2) near r_eq. r_eq must be positive.

        :param r:
            The bond length(s) in bohr.
        :param r_eq:
            The equilibrium bond length(s) in bohr.
        """
        r_input = np.asarray(r, dtype=np.float64)
        r_eq_input = np.asarray(r_eq, dtype=np.float64)
        scalar_input = (r_input.ndim == 0 and r_eq_input.ndim == 0)

        r_arr, r_eq_arr = np.broadcast_arrays(np.atleast_1d(r_input),
                                              np.atleast_1d(r_eq_input))
        r_arr = r_arr.astype(np.float64, copy=False)
        r_eq_arr = r_eq_arr.astype(np.float64, copy=False)

        assert_msg_critical(
            np.all(r_eq_arr > 0.0),
            'InterpolationDatapoint: Equilibrium bond lengths must be '
            'positive for use_eq_bond_length.')

        R = - np.square(r_eq_arr) * (1.0 / r_arr - 1.0 / r_eq_arr)
        dR = np.square(r_eq_arr) / np.square(r_arr)
        d2R = -np.square(r_eq_arr) * 2.0 / np.power(r_arr, 3)

        q = R
        dq_dr = dR
        d2q_dr2 = d2R

        if scalar_input:
            return float(q[0]), float(dq_dr[0]), float(d2q_dr2[0])

        return q, dq_dr, d2q_dr2

    def compute_internal_coordinates_values(self):
        """
        Creates an array with the values of the internal coordinates
        and saves it in self.
        """
        self._compute_internal_coordinates_values_vectorized()

    def get_z_matrix_as_np_arrays(self):
        """
        Returns the z-matrix as integer arrays with the keys bonds, angles,
        dihedrals and impropers (2, 3, 4 and 4 atom indices per row).
        """

        zmat = self.z_matrix_dict
        return {
            "bonds": np.array(zmat["bonds"], dtype=np.int64),
            "angles": np.array(zmat["angles"], dtype=np.int64),
            "dihedrals": np.array(zmat["dihedrals"], dtype=np.int64),
            "impropers": np.array(zmat["impropers"], dtype=np.int64),
        }

    def get_imp_int_coord_as_np_arrays(self):
        """
        Returns the important internal coordinates as integer arrays with
        the keys imp_bonds, imp_angles, imp_dihedrals and imp_impropers (2, 3,
        4 and 4 atom indices per row). An empty list or missing keys give
        empty arrays.
        """
        assert_msg_critical(self.imp_int_coordinates is not None,
                            'InterpolationDatapoint: No important internal coordinates are defined.')

        imp_coords = self.imp_int_coordinates
        if not isinstance(imp_coords, dict):
            # an empty sequence (the former default) means "no important
            # coordinates"; anything else cannot be stored
            assert_msg_critical(
                isinstance(imp_coords, (list, tuple)) and len(imp_coords) == 0,
                'InterpolationDatapoint: imp_int_coordinates must be a '
                'dictionary with the keys bonds, angles, dihedrals and '
                'impropers.')
            imp_coords = {}

        return {
            "imp_bonds": np.array(imp_coords.get("bonds", []), dtype=np.int64),
            "imp_angles": np.array(imp_coords.get("angles", []), dtype=np.int64),
            "imp_dihedrals": np.array(imp_coords.get("dihedrals", []), dtype=np.int64),
            "imp_impropers": np.array(imp_coords.get("impropers", []), dtype=np.int64),
        }

    def calculate_translation_coordinates(self):
        """
        Returns cartesian_coordinates translated to their geometric centre.
        """
        center = np.mean(self.cartesian_coordinates, axis=0)
        translated_coordinates = self.cartesian_coordinates - center

        return translated_coordinates

    def _write_optional_array(self, h5f, name, value):
        """
        Writes value as dataset name (gzip-compressed unless it is a
        scalar) and returns the dataset; does nothing for None.
        """
        if value is None:
            return None
        array_value = np.asarray(value)
        if array_value.ndim == 0:
            return h5f.create_dataset(name, data=array_value)
        else:
            return h5f.create_dataset(name, data=array_value, compression="gzip")

    @staticmethod
    def _bond_mode_suffix_from_flags(use_inverse_bond_length, use_eq_bond_length):
        """
        Returns the label suffix of a bond-coordinate mode.
        """

        if use_inverse_bond_length:
            return '_rinv'
        elif use_eq_bond_length:
            return '_eq'
        return '_r'

    def _full_point_label(self, label):
        """
        Returns the HDF5 dataset prefix of a datapoint label for the current
        bond-coordinate mode, e.g. 'point_1' -> 'point_1_eq_dihedral'.
        """

        return label + self._bond_mode_suffix_from_flags(
            self.use_inverse_bond_length,
            self.use_eq_bond_length) + '_dihedral'

    @staticmethod
    def _checked_file_name(fname, caller):
        """
        Returns fname as str; raises for anything that is not a non-empty
        str or os.PathLike.
        """

        assert_msg_critical(
            isinstance(fname, (str, os.PathLike)),
            f'InterpolationDatapoint.{caller}: the file name must be a str '
            f'or os.PathLike, got {type(fname).__name__}.')
        fname = os.fspath(fname)
        assert_msg_critical(
            isinstance(fname, str) and fname != '',
            f'InterpolationDatapoint.{caller}: invalid file name {fname!r}.')
        return fname

    def _missing_point_message(self, h5f, fname, full_label):
        """
        Builds the error message for a datapoint that is not in the file,
        including the bond-mode suffixes that the file does contain.
        """

        found = set()
        tail = '_dihedral_energy'
        for key in h5f.keys():
            if key.endswith(tail):
                found.add('_' + key[:-len(tail)].rsplit('_', 1)[-1] +
                          '_dihedral')
        found_txt = ', '.join(sorted(found)) if found else 'none'
        suffix = self._bond_mode_suffix_from_flags(
            self.use_inverse_bond_length,
            self.use_eq_bond_length) + '_dihedral'

        msg = f'InterpolationDatapoint.read_hdf5: datapoint {full_label} '
        msg += f'is not in {fname}. With use_inverse_bond_length='
        msg += f'{self.use_inverse_bond_length} and use_eq_bond_length='
        msg += f'{self.use_eq_bond_length} the datasets must be named '
        msg += f'<label>{suffix}_*; the file contains datapoints with the '
        msg += f'suffix(es): {found_txt}.'
        return msg

    def write_hdf5(self, fname, label):
        """
        Writes the datapoint to an HDF5 database.

        The datasets are named <label><suffix>_dihedral_<name> with the
        bond-mode suffix '_rinv', '_eq' or '_r', and <name> one of energy,
        gradient and hessian (internal), cartesian_gradient and
        cartesian_hessian (in the build metric, recorded in the attributes
        'mass_weighted' and 'bond_mode'), internal_coordinates,
        eq_bond_lengths (if set), confidence_radius, cartesian_coordinates,
        the z-matrix (bonds, angles, dihedrals, impropers), the important
        coordinates (imp_bonds, ...) and atom_labels (if set).

        All inputs are checked before the file is opened. The file is
        created if it does not exist and is never truncated. A label that
        exists already raises, and a failing write removes the datasets it
        wrote, so the database stays readable.

        :param fname:
            The HDF5 file name (str or os.PathLike).
        :param label:
            The datapoint label without suffix, e.g. 'point_12'.
        """

        fname = self._checked_file_name(fname, 'write_hdf5')
        label = self._full_point_label(label)

        # Validate every input before the file is opened, so that a failed
        # write can not leave a partial datapoint behind.
        assert_msg_critical(self.energy is not None,
                            'InterpolationDatapoint: No energy is defined.')

        assert_msg_critical(
            self.internal_gradient is not None,
            'InterpolationDatapoint: No internal gradient is defined.')

        assert_msg_critical(
            self.internal_hessian is not None,
            'InterpolationDatapoint: No internal Hessian is defined.')

        assert_msg_critical(
            self.gradient is not None,
            'InterpolationDatapoint: No Cartesian gradient is defined.')

        assert_msg_critical(
            self.hessian is not None,
            'InterpolationDatapoint: No Cartesian Hessian is defined.')

        if self.use_eq_bond_length:
            assert_msg_critical(
                self.eq_bond_lengths is not None,
                'InterpolationDatapoint: No eq bond lengths are defined.')

        assert_msg_critical(
            self.confidence_radius is not None,
            'InterpolationDatapoint: No confidence radius is defined.')

        assert_msg_critical(
            self.cartesian_coordinates is not None,
            'InterpolationDatapoint: No Cartesian coordinates are defined.')

        assert_msg_critical(
            self.z_matrix is not None and self.z_matrix_dict is not None,
            'InterpolationDatapoint: No Z-matrix is defined.')

        if self.internal_coordinates is None:
            assert_msg_critical(
                self.internal_coordinates_values is not None,
                'InterpolationDatapoint: No internal coordinates are defined.')
        else:
            self.compute_internal_coordinates_values()

        # Bonds, angles, dihedrals and impropers of the Z-matrix, one array
        # per type with 2, 3 or 4 atom indices per row.
        zmat_dict = self.get_z_matrix_as_np_arrays()

        imp_int_coord_dict = {}
        if self.identify_imp_int_coord:
            assert_msg_critical(
                self.imp_int_coordinates is not None,
                'InterpolationDatapoint: No important internal coordinates are defined.')
            imp_int_coord_dict = self.get_imp_int_coord_as_np_arrays()

        # Metric of the stored Cartesian gradient and Hessian: they are
        # written as given, i.e. mass-weighted if the point was built with a
        # mass-weighted B-matrix (inv_sqrt_masses set).
        cartesian_attrs = {
            'mass_weighted': bool(self.inv_sqrt_masses is not None),
            'bond_mode': self._bond_mode_suffix_from_flags(
                self.use_inverse_bond_length, self.use_eq_bond_length)[1:],
        }

        datasets = [
            (label + '_energy', np.float64(self.energy), None),
            (label + '_gradient', self.internal_gradient, None),
            (label + '_cartesian_gradient', self.gradient, cartesian_attrs),
            (label + '_hessian', self.internal_hessian, None),
            (label + '_cartesian_hessian', self.hessian, cartesian_attrs),
            (label + '_internal_coordinates',
             self.internal_coordinates_values, None),
            (label + '_eq_bond_lengths', self.eq_bond_lengths, None),
            (label + '_confidence_radius', np.float64(self.confidence_radius),
             None),
            (label + '_cartesian_coordinates', self.cartesian_coordinates,
             None),
        ]
        for key, value in zmat_dict.items():
            datasets.append((label + '_' + key, value, None))
        for key, value in imp_int_coord_dict.items():
            datasets.append((label + '_' + key, value, None))

        names = [name for name, value, _ in datasets if value is not None]
        if self.atom_labels is not None:
            names.insert(0, label + '_atom_labels')

        with h5py.File(fname, 'a') as h5f:
            # Mode 'a' creates the file if it does not exist. There is no
            # fallback to mode 'w': that would truncate an existing database
            # that is locked, corrupted or not writable.

            existing = [name for name in names if name in h5f]
            assert_msg_critical(
                not existing,
                f'InterpolationDatapoint.write_hdf5: datapoint {label} '
                f'already exists in {fname} ({len(existing)} of its '
                'datasets are present).')

            written = []
            try:
                if self.atom_labels is not None:
                    string_type = h5py.string_dtype(encoding='utf-8')
                    h5f.create_dataset(label + '_atom_labels',
                                       data=self.atom_labels,
                                       dtype=string_type)
                    written.append(label + '_atom_labels')

                for name, value, attrs in datasets:
                    dataset = self._write_optional_array(h5f, name, value)
                    if dataset is None:
                        continue
                    written.append(name)
                    if attrs:
                        for key, attr_value in attrs.items():
                            dataset.attrs[key] = attr_value

            except BaseException:
                # remove what was written of this datapoint, so that the
                # database stays readable, then re-raise
                for name in written:
                    if name in h5f:
                        del h5f[name]
                raise

    def read_hdf5(self, fname, label):
        """
        Reads a datapoint written by write_hdf5.

        The bond-mode settings select the label suffix and must be those of
        the database; a missing label raises with the suffixes found in the
        file. Missing datasets give None. The stored z-matrix replaces
        z_matrix and z_matrix_dict, the important coordinates and the atom
        labels are read if present, and point_label is set to the dataset
        prefix.

        :param fname:
            The HDF5 file name (str or os.PathLike).
        :param label:
            The datapoint label without suffix, e.g. 'point_12'.
        """

        fname = self._checked_file_name(fname, 'read_hdf5')
        label = self._full_point_label(label)

        with h5py.File(fname, 'r') as h5f:

            energy_label = label + "_energy"
            if energy_label not in h5f:
                assert_msg_critical(
                    False, self._missing_point_message(h5f, fname, label))

            def _read_array(name):
                # a missing dataset gives None, not array(None, dtype=object)
                dataset = h5f.get(name)
                return None if dataset is None else np.array(dataset)

            self.point_label = label
            self.energy = _read_array(energy_label)
            self.internal_gradient = _read_array(label + "_gradient")
            self.gradient = _read_array(label + "_cartesian_gradient")
            self.internal_hessian = _read_array(label + "_hessian")
            self.hessian = _read_array(label + "_cartesian_hessian")
            self.internal_coordinates_values = _read_array(
                label + "_internal_coordinates")
            self.eq_bond_lengths = _read_array(label + "_eq_bond_lengths")
            self.cartesian_coordinates = _read_array(
                label + "_cartesian_coordinates")
            self.confidence_radius = _read_array(label + "_confidence_radius")

            atom_labels = h5f.get(label + '_atom_labels')
            if atom_labels is not None:
                self.atom_labels = [
                    x.decode('utf-8') if isinstance(x, bytes) else str(x)
                    for x in atom_labels[()]
                ]

            # The Z-matrix of the point replaces the constructor Z-matrix
            # (both the dictionary and the flat list).
            z_matrix_dict = {}
            for key in ("bonds", "angles", "dihedrals", "impropers"):
                ds = h5f.get(label + '_' + key)
                if ds is not None:
                    z_matrix_dict[key] = [tuple(x.tolist()) for x in ds]
            if z_matrix_dict:
                for key in ("bonds", "angles", "dihedrals", "impropers"):
                    z_matrix_dict.setdefault(key, [])
                self.z_matrix_dict = z_matrix_dict
                self.z_matrix = self.flatten_z_matrix(z_matrix_dict)
            self._b_matrix_row_cache = None

            if self.identify_imp_int_coord:
                self.imp_int_coordinates = {}
                for key in ("bonds", "angles", "dihedrals", "impropers"):
                    ds = h5f.get(label + '_imp_' + key)
                    if ds is not None:
                        self.imp_int_coordinates[key] = [tuple(x.tolist()) for x in ds]
                    else:
                        self.imp_int_coordinates[key] = []

    def cartesian_distance_vector(self, data_point):
        """
        Returns the displacement between this datapoint and data_point after
        centring and rotating data_point onto this geometry (Kabsch over all
        atoms; the driver weights use only the atoms of
        symmetry_information[3]).

        :param data_point:
            The InterpolationDatapoint.

        :return:
            The displacement (natoms, 3) in bohr.
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
        Removes all datasets of a datapoint from an HDF5 file.

        :param fname:
            The HDF5 file name (str or os.PathLike).
        :param label:
            The datapoint label without suffix, e.g. 'point_12'.
        :param use_inverse_bond_length:
            The bond mode of the label (suffix '_rinv'; takes precedence).
        :param use_eq_bond_length:
            The bond mode of the label (suffix '_eq'); both False: '_r'.
        """

        fname = self._checked_file_name(fname, 'remove_point_from_hdf5')
        label += self._bond_mode_suffix_from_flags(use_inverse_bond_length,
                                                   use_eq_bond_length)
        label += "_dihedral"

        keys_to_remove = [
            label + "_energy",
            label + "_gradient",
            label + "_hessian",
            label + "_internal_coordinates",
            label + "_eq_bond_lengths",
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
            label + "_imp_impropers",
            label + "_atom_labels"
        ]

        deleted = []
        with h5py.File(fname, 'r+') as h5f:
            for key in keys_to_remove:
                if key in h5f:
                    del h5f[key]
                    deleted.append(key)

        skipped = [key[len(label) + 1:] for key in keys_to_remove
                   if key not in deleted]
        msg = f'Removed {len(deleted)} datasets of {label} from {fname}'
        if skipped:
            msg += f' (not present: {", ".join(skipped)})'
        self.ostream.print_info(msg + '.')
        self.ostream.flush()

    def update_confidence_radius(self, fname, label, new_confidence_radius):
        """
        Updates the confidence radius of a datapoint in the HDF5 file (bond
        mode of this datapoint).

        :param fname:
            The HDF5 file name (str or os.PathLike).
        :param label:
            The datapoint label without suffix, e.g. 'point_12'.
        :param new_confidence_radius:
            The new radius in bohr (same shape as the stored value).
        """

        fname = self._checked_file_name(fname, 'update_confidence_radius')

        with h5py.File(fname, 'r+') as h5f:
            label = self._full_point_label(label)

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
