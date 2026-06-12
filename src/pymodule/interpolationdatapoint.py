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
        # self.use_cos_angle = False
        self.use_vectorized_b_matrix = True
        self.use_vectorized_internal_coordinates_values = True
        self.identify_imp_int_coord = True
        
        self.eq_bond_lengths = None
        
        # internal_coordinates is a list of geomeTRIC objects which represent
        # different types of internal coordinates (distances, angles, dihedrals)
        self.internal_coordinates = None
        self.imp_int_coordinates = []


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


        im_keywords: dict[str, str] = {
            key: val[0]
            for key, val in self._input_keywords['im_settings'].items()
        }

        parse_input(self, im_keywords, impes_dict)

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
            self.use_eq_bond_length #or 
            # self.use_cos_angle
        )
        
        if use_original:
            self.original_b_matrix = derivatives.copy()
        else:
            self.original_b_matrix = None

        row_scale = np.ones(n_rows, dtype=np.float64)
        row_cache = self._get_b_matrix_row_cache()

        bond_rows = row_cache['bond_rows']
        angle_rows = row_cache['angle_rows']
        if bond_rows.size > 0 and (self.use_inverse_bond_length or self.use_eq_bond_length):
            bond_values = np.array(
                [self.internal_coordinates[idx].value(coords) for idx in bond_rows],
                dtype=np.float64)

            if self.use_inverse_bond_length:
                row_scale[bond_rows] = -1.0 / np.square(bond_values)
            elif self.use_eq_bond_length:
                eq_values = self._get_eq_bond_lengths_array(bond_rows.size)
                _, dq_dr, _ = self._switched_bond_transform(
                    bond_values, eq_values, eps_inner=0.005, eps_outer=0.01)
                row_scale[bond_rows] = dq_dr

        # if angle_rows.size > 0 and self.use_cos_angle:
        #     angle_values = np.array(
        #     [self.internal_coordinates[idx].value(coords) for idx in angle_rows],
        #     dtype=np.float64)

        #     row_scale[angle_rows] = -np.sin(angle_values)
        
        self.b_matrix = derivatives * row_scale[:, np.newaxis]

        if self.inv_sqrt_masses is not None:
            self.b_matrix = self.b_matrix * self.inv_sqrt_masses
            if self.use_inverse_bond_length or self.use_eq_bond_length:
                self.original_b_matrix = self.original_b_matrix

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

        eq_bond_values = None
        if self.use_eq_bond_length:
            bond_rows = self._get_b_matrix_row_cache()['bond_rows']
            eq_bond_values = self._get_eq_bond_lengths_array(bond_rows.size)

        
        # prev_dihedral = None
        bond_counter = 0
        for i, z in enumerate(self.z_matrix):
            q = self.internal_coordinates[i]
            second_derivative = q.second_derivative(coords).reshape(-1, n_atoms * 3)

            if len(z) == 2:
                if self.use_inverse_bond_length:

                    r = q.value(coords)
                    r_inv_2 = 1.0 / (r * r)
                    r_inv_3 = r_inv_2 / r
                    self.b2_matrix[i] = -r_inv_2 * second_derivative

                    for m in range(n_atoms):
                        for n in range(n_atoms):
                            self.b2_matrix[i, m*3:(m+1)*3, n*3:(n+1)*3] += 2 * r_inv_3 * np.outer(self.original_b_matrix[i, m*3:(m+1)*3], self.original_b_matrix[i, n*3:(n+1)*3])

                elif self.use_eq_bond_length:
                    r = q.value(coords)
                    _, r_deriv_2, r_deriv_3 = self._switched_bond_transform(
                        r, eq_bond_values[bond_counter], eps_inner=0.005, eps_outer=0.01)
                    self.b2_matrix[i] = r_deriv_2 * second_derivative
                    for m in range(n_atoms):
                        for n in range(n_atoms):
                            self.b2_matrix[i, m*3:(m+1)*3, n*3:(n+1)*3] += r_deriv_3 * np.outer(self.original_b_matrix[i, m*3:(m+1)*3], self.original_b_matrix[i, n*3:(n+1)*3])

                else:
                    self.b2_matrix[i] = second_derivative

                bond_counter += 1
            
            # elif len(z) == 3 and self.use_cos_angle:

            #     a = q.value(coords)
            #     a_sin_2 = -np.sin(a)
            #     a_cos_3 = -np.cos(a)
            #     self.b2_matrix[i] = a_sin_2 * second_derivative

            #     for m in range(n_atoms):
            #         for n in range(n_atoms):
            #             self.b2_matrix[i, m*3:(m+1)*3, n*3:(n+1)*3] += a_cos_3 * np.outer(self.original_b_matrix[i, m*3:(m+1)*3], self.original_b_matrix[i, n*3:(n+1)*3])

            
            else:
                self.b2_matrix[i] = second_derivative

        if self.inv_sqrt_masses is not None:
            self.b2_matrix = (self.b2_matrix *
                self.inv_sqrt_masses.reshape(1, -1, 1) *
                self.inv_sqrt_masses.reshape(1,  1, -1))

               
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

        number_of_positive_values = np.count_nonzero(s_inv)

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

        if bond_rows.size > 0:
            if self.use_inverse_bond_length:
                int_coords[bond_rows] = 1.0 / base_values[bond_rows]
            elif self.use_eq_bond_length:
               
                eq_values = self._get_eq_bond_lengths_array(bond_rows.size)
                q_values, _, _ = self._switched_bond_transform(
                    base_values[bond_rows],
                    eq_values,
                    eps_inner=0.005,
                    eps_outer=0.01)
                int_coords[bond_rows] = q_values
        # if angle_rows.size > 0 and self.use_cos_angle:
        #     int_coords[angle_rows] = np.cos(base_values[angle_rows])

        self.internal_coordinates_values = int_coords
    
    def smoothstep5(self, t):
        return 1.0 - 10.0*t**3 + 15.0*t**4 - 6.0*t**5

    def dsmoothstep5_dt(self, t):
        return -30.0*t**2 + 60.0*t**3 - 30.0*t**4

    def d2smoothstep5_dt2(self, t):
        return -60.0*t + 180.0*t**2 - 120.0*t**3


    def _switched_bond_transform(self, r, r_eq, eps_inner=0.05, eps_outer=0.15):
        """
        Returns q(r), dq/dr, d2q/dr2 for a bond coordinate that is
        smoothly switched from local r behavior near equilibrium to 1/r outside.
        """
        r_input = np.asarray(r, dtype=np.float64)
        r_eq_input = np.asarray(r_eq, dtype=np.float64)
        scalar_input = (r_input.ndim == 0 and r_eq_input.ndim == 0)

        r_arr, r_eq_arr = np.broadcast_arrays(np.atleast_1d(r_input),
                                              np.atleast_1d(r_eq_input))
        r_arr = r_arr.astype(np.float64, copy=False)
        r_eq_arr = r_eq_arr.astype(np.float64, copy=False)

        L = r_arr - r_eq_arr
        dL = np.ones_like(r_arr)
        d2L = np.zeros_like(r_arr)


        # R = - np.square(r_eq_arr) * (1.0 / r_arr - 1.0 / r_eq_arr)
        # dR = np.square(r_eq_arr) / np.square(r_arr)
        # d2R = -np.square(r_eq_arr) * 2.0 / np.power(r_arr, 3)

        R = - np.square(r_eq_arr) * (1.0 / r_arr - 1.0 / r_eq_arr)
        dR = np.square(r_eq_arr) / np.square(r_arr)
        d2R = -np.square(r_eq_arr) * 2.0 / np.power(r_arr, 3)

        x = np.log(r_arr / r_eq_arr)
        y = np.square(x)

        y1 = np.log(1.0 + eps_inner)**2
        y2 = np.log(1.0 + eps_outer)**2
        denom = y2 - y1

        s = np.ones_like(y)
        ds = np.zeros_like(y)
        d2s = np.zeros_like(y)

        mask_inner = (y <= y1)
        mask_outer = (y >= y2)
        mask_transition = (~mask_inner) & (~mask_outer)

        if np.any(mask_transition):
            t = (y[mask_transition] - y1) / denom
            P = self.smoothstep5(t)
            dP = self.dsmoothstep5_dt(t)
            d2P = self.d2smoothstep5_dt2(t)

            dy_dr = 2.0 * x[mask_transition] / r_arr[mask_transition]
            d2y_dr2 = 2.0 * (1.0 - x[mask_transition]) / np.square(r_arr[mask_transition])
            dt_dr = dy_dr / denom
            d2t_dr2 = d2y_dr2 / denom

            s[mask_transition] = P
            ds[mask_transition] = dP * dt_dr
            d2s[mask_transition] = d2P * np.square(dt_dr) + dP * d2t_dr2

        s[mask_outer] = 0.0
        ds[mask_outer] = 0.0
        d2s[mask_outer] = 0.0
        
        q = s * L + (1.0 - s) * R
        
        q = R
        dq_dr = dR
        d2q_dr2 = d2R
        # here the switch derivativve is being assembled 
        # dq_dr = ds * (L - R) + s * dL + (1.0 - s) * dR
        # d2q_dr2 = d2s * (L - R) + 2.0 * ds * (dL - dR) + s * d2L + (1.0 - s) * d2R

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
        if array_value.ndim == 0:
            h5f.create_dataset(name, data=array_value)
        else:
            h5f.create_dataset(name, data=array_value, compression="gzip")
        

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

            if self.use_inverse_bond_length:
                label += "_rinv"
            elif self.use_eq_bond_length:
                label += "_eq"
            else:
                label += "_r"

            
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

        valid_checkpoint = (fname and isinstance(fname, str))

        if valid_checkpoint:
            h5f = h5py.File(fname, 'r')

            if self.use_inverse_bond_length:
                label += "_rinv"
            elif self.use_eq_bond_length:
                label += "_eq"
            else:
                label += "_r"

            label += "_dihedral"

            energy_label = label + "_energy"
            gradient_label = label + "_gradient"
            hessian_label = label + "_hessian"
            cart_hess_label = label + "_cartesian_hessian"
            cart_grad_label = label + "_cartesian_gradient"
            coords_label = label + "_internal_coordinates"
            eq_bond_length_label = label + "_eq_bond_lengths"
            cart_coords_label = label + "_cartesian_coordinates"
            confidence_radius_label = label + "_confidence_radius"
        

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
            self.cartesian_coordinates = np.array(h5f.get(cart_coords_label))
            self.confidence_radius = np.array(h5f.get(confidence_radius_label))


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
            if self.use_inverse_bond_length:
                label += "_rinv"
            elif self.use_eq_bond_length:
                label += "_eq"
            else:
                label += "_r"

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