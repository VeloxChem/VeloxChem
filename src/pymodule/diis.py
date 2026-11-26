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

import numpy as np

from .mathutils import safe_solve


class Diis:
    """
    Implements direct inversion of the iterative subspace.

    Instance variables
        - error_vectors: The list of error vectors.
        - fock_matrices: The list of stored Fock/Kohn-Sham matrices.
        - fock_matrices_proj: The list of stored projected Fock/Kohn-Sham
          matrices (used in restricted open-shell SCF).
        - max_err_vecs: The maximum number of error vectors.
        - diis_thresh: The DIIS switch-on threshold.
        - scf_type: The type of SCF calculation.
        - b_matrix: The B matrix of error-vector inner products.
    """

    def __init__(self, max_err_vecs, diis_thresh, scf_type):
        """
        Initializes the DIIS driver.

        :param max_err_vecs:
            The maximum number of error vectors.
        :param diis_thresh:
            The DIIS switch-on threshold.
        :param scf_type:
            The type of SCF calculation ('restricted', 'unrestricted' or
            'restricted_openshell').
        """

        self.error_vectors = []

        self.fock_matrices = []
        self.fock_matrices_proj = []

        self.max_err_vecs = max_err_vecs
        self.diis_thresh = diis_thresh
        self.scf_type = scf_type

        self.b_matrix = np.zeros((max_err_vecs, max_err_vecs))

    def clear(self):
        """
        Clears the stored error vectors, Fock/Kohn-Sham matrices and B matrix.
        """

        self.error_vectors.clear()

        self.fock_matrices.clear()
        self.fock_matrices_proj.clear()

        self.b_matrix = np.zeros((self.max_err_vecs, self.max_err_vecs))

    def store_diis_data(self, fock_mat, den_mat, ovl_mat, e_mat, e_grad):
        """
        Stores error vector and Fock/Kohn-Sham matrix for the current
        iteration and updates the B matrix. For restricted open-shell SCF,
        the projected Fock/Kohn-Sham matrix is also stored.

        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param den_mat:
            The density matrix.
        :param ovl_mat:
            The overlap matrix (used in ROSCF).
        :param e_mat:
            The error vector.
        :param e_grad:
            The electronic gradient.
        """

        if e_grad < self.diis_thresh:

            if len(self.error_vectors) == self.max_err_vecs:
                self.error_vectors.pop(0)
                self.fock_matrices.pop(0)
                if self.scf_type == 'restricted_openshell':
                    self.fock_matrices_proj.pop(0)
                sub_bmat = self.b_matrix[1:, 1:].copy()
                self.b_matrix[:-1, :-1] = sub_bmat[:, :]

            self.error_vectors.append(e_mat.copy())
            self.fock_matrices.append([x.copy() for x in fock_mat])
            if self.scf_type == 'restricted_openshell':
                fock_proj = self.get_projected_fock(
                    fock_mat[0], fock_mat[1], den_mat[0], den_mat[1], ovl_mat)
                # Note: append a list
                self.fock_matrices_proj.append([fock_proj])

            n_vecs = len(self.error_vectors)
            for i in range(n_vecs):
                fij = np.vdot(self.error_vectors[i],
                              self.error_vectors[n_vecs - 1])
                self.b_matrix[i, n_vecs - 1] = fij
                self.b_matrix[n_vecs - 1, i] = fij

    def get_effective_fock(self, fock_mat):
        """
        Computes effective Fock/Kohn-Sham matrix by DIIS extrapolation of
        stored Fock/Kohn-Sham matrices.

        :param fock_mat:
            The current Fock/Kohn-Sham matrices, used when no error vectors
            have been stored.

        :return:
            The effective Fock/Kohn-Sham matrices as a tuple.
        """

        n_vecs = len(self.error_vectors)

        if n_vecs == 0:
            # No error vectors stored, e.g. when the electronic gradient is
            # still above the DIIS threshold. Use the current Fock/Kohn-Sham
            # matrices.
            return tuple(fock_mat)

        if n_vecs == 1:
            if self.scf_type == 'restricted_openshell':
                return tuple(self.fock_matrices_proj[0])
            else:
                return tuple(self.fock_matrices[0])

        else:
            weights = self.compute_weights()

            if self.scf_type == 'restricted':
                fock_matrices_a = [m[0] for m in self.fock_matrices]
                effmat_a = self._weighted_sum(weights, fock_matrices_a)
                # Note: return a tuple
                return (effmat_a,)

            elif self.scf_type == 'unrestricted':
                fock_matrices_a = [m[0] for m in self.fock_matrices]
                fock_matrices_b = [m[1] for m in self.fock_matrices]
                effmat_a = self._weighted_sum(weights, fock_matrices_a)
                effmat_b = self._weighted_sum(weights, fock_matrices_b)
                return (effmat_a, effmat_b)

            else:
                eff_fock_matrices = [m[0] for m in self.fock_matrices_proj]
                effmat = self._weighted_sum(weights, eff_fock_matrices)
                # Note: return a tuple
                return (effmat,)

    def compute_weights(self):
        """
        Computes DIIS weights from error vectors.

        :return:
            The DIIS weights.
        """

        n_vecs = len(self.error_vectors)

        bmat = np.zeros((n_vecs + 1, n_vecs + 1))
        bmat[:n_vecs, :n_vecs] = self.b_matrix[:n_vecs, :n_vecs]
        bmat[n_vecs, :n_vecs] = -1.0
        bmat[:n_vecs, n_vecs] = -1.0
        bmat[n_vecs, n_vecs] = 0.0

        bvec = np.zeros(n_vecs + 1)
        bvec[:n_vecs] = 0.0
        bvec[n_vecs] = -1.0

        return safe_solve(bmat, bvec)[:n_vecs]

    @staticmethod
    def _weighted_sum(weights, matrices):
        """
        Computes the weighted sum of matrices.

        :param weights:
            The weights.
        :param matrices:
            The matrices.

        :return:
            The weighted sum of matrices.
        """

        return sum([w * mat for w, mat in zip(weights, matrices)])

    @staticmethod
    def get_projected_fock(fa, fb, da, db, s):
        """
        Generates projected Fock matrix.

        :param fa:
            The Fock matrix of alpha spin.
        :param fb:
            The Fock matrix of beta spin.
        :param da:
            The density matrix of alpha spin.
        :param db:
            The density matrix of beta spin.
        :param s:
            The overlap matrix.

        :return:
            The projected Fock matrix.
        """

        naos = s.shape[0]

        inactive = np.matmul(s, db)
        active = np.matmul(s, da - db)
        virtual = np.eye(naos) - np.matmul(s, da)

        #       occ   act   vir
        #     +----------------+
        # occ | f0    fb    f0 |
        # act | fb    f0    fa |
        # vir | f0    fa    f0 |
        #     +----------------+

        f0 = 0.5 * (fa + fb)

        fcorr = np.linalg.multi_dot([inactive, fb - f0, active.T])
        fcorr += np.linalg.multi_dot([active, fa - f0, virtual.T])
        fcorr += fcorr.T

        return f0 + fcorr
