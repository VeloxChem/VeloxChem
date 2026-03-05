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

from .errorhandler import safe_solve

import numpy as np


class Diis:
    """
    Conventional DIIS driver that solves the augmented linear system.

    Error vectors are built from (FDS - SDF) in the orthogonal AO basis so
    the residuals match those used by the SCF drivers. The system enforces
    the normalization constraint and fallbacks to a regularized least
    squares solve when the B-matrix is singular.
    """

    def __init__(self):
        """
        Initializes the internal error-vector buffer.
        """

        self.error_vectors = []

    def compute_error_vectors_restricted(self, fock_matrices, density_matrices,
                                         overlap_matrix, oao_matrix):
        """
        Computes restricted-spin error vectors.

        :param fock_matrices:
            The list of AO Fock/Kohn-Sham matrices.
        :param density_matrices:
            The corresponding AO density matrices.
        :param overlap_matrix:
            The AO overlap matrix.
        :param oao_matrix:
            The orthogonalization matrix.
        """

        self._build_error_vectors(fock_matrices, density_matrices, overlap_matrix,
                                  oao_matrix)

    def compute_error_vectors_restricted_openshell(self, fock_matrices,
                                                   fock_matrices_beta,
                                                   density_matrices,
                                                   density_matrices_beta,
                                                   overlap_matrix, oao_matrix):
        """
        Computes restricted open-shell error vectors.

        :param fock_matrices:
            The alpha-spin AO Fock matrices.
        :param fock_matrices_beta:
            The beta-spin AO Fock matrices.
        :param density_matrices:
            The alpha-spin density matrices.
        :param density_matrices_beta:
            The beta-spin density matrices.
        :param overlap_matrix:
            The AO overlap matrix.
        :param oao_matrix:
            The orthogonalization matrix.
        """

        alpha_err = self._collect_error_matrices(fock_matrices, density_matrices,
                                                 overlap_matrix, oao_matrix)
        beta_err = self._collect_error_matrices(fock_matrices_beta,
                                                density_matrices_beta,
                                                overlap_matrix, oao_matrix)

        self.error_vectors = [a + b for a, b in zip(alpha_err, beta_err)]

    def compute_error_vectors_unrestricted(self, fock_matrices,
                                           fock_matrices_beta, density_matrices,
                                           density_matrices_beta,
                                           overlap_matrix, oao_matrix):
        """
        Computes unrestricted error vectors.

        :param fock_matrices:
            The alpha-spin AO Fock matrices.
        :param fock_matrices_beta:
            The beta-spin AO Fock matrices.
        :param density_matrices:
            The alpha-spin density matrices.
        :param density_matrices_beta:
            The beta-spin density matrices.
        :param overlap_matrix:
            The AO overlap matrix.
        :param oao_matrix:
            The orthogonalization matrix.
        """

        alpha_err = self._collect_error_matrices(fock_matrices, density_matrices,
                                                 overlap_matrix, oao_matrix)
        beta_err = self._collect_error_matrices(fock_matrices_beta,
                                                density_matrices_beta,
                                                overlap_matrix, oao_matrix)

        self.error_vectors = [
            np.vstack((a, b)) for a, b in zip(alpha_err, beta_err)
        ]

    def compute_weights(self):
        """
        Computes DIIS weights from the stored error vectors.

        :return:
            The weights normalized to sum to unity.
        :raises ValueError:
            If no error vectors have been collected.
        """

        if len(self.error_vectors) == 0:
            raise ValueError('DIIS: no error vectors available')

        if len(self.error_vectors) == 1:
            return np.array([1.0], dtype='float64')

        bmat = self._build_bmatrix()
        dim = bmat.shape[0]
        aug = np.zeros((dim + 1, dim + 1), dtype='float64')
        aug[:dim, :dim] = bmat
        aug[:dim, dim] = -1.0
        aug[dim, :dim] = -1.0

        rhs = np.zeros(dim + 1, dtype='float64')
        rhs[dim] = -1.0

        sol = safe_solve(aug, rhs)

        weights = sol[:dim]
        return weights

    def _collect_error_matrices(self, fock_matrices, density_matrices,
                                overlap_matrix, oao_matrix):
        """
        Builds spin-channel error matrices.

        :param fock_matrices:
            The AO Fock matrices for one spin.
        :param density_matrices:
            The AO density matrices for the same spin.
        :param overlap_matrix:
            The AO overlap matrix.
        :param oao_matrix:
            The orthogonalization matrix.
        :return:
            The list of residual arrays.
        """

        smat = overlap_matrix
        tmat = oao_matrix
        errs = []
        for fmat, dmat in zip(fock_matrices, density_matrices):
            fds = np.matmul(fmat, np.matmul(dmat, smat))
            err = np.matmul(tmat.T, np.matmul(fds - fds.T, tmat))
            errs.append(err)
        return errs

    def _build_error_vectors(self, fock_matrices, density_matrices,
                             overlap_matrix, oao_matrix):
        """
        Helper that stores the collected residuals.
        """

        errs = self._collect_error_matrices(fock_matrices, density_matrices,
                                            overlap_matrix, oao_matrix)
        self.error_vectors = [err for err in errs]

    def _build_bmatrix(self):
        """
        Builds the symmetric B matrix from flattened error vectors.

        :return:
            The B-matrix for the DIIS augmented system.
        """

        dim = len(self.error_vectors)
        bmat = np.zeros((dim, dim), dtype='float64')
        flat = [vec.reshape(-1) for vec in self.error_vectors]
        for i in range(dim):
            for j in range(i, dim):
                value = np.vdot(flat[i], flat[j])
                bmat[i, j] = value
                bmat[j, i] = value
        return bmat
