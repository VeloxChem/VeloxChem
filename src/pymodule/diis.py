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

from .veloxchemlib import weighted_sum_gpu, dot_product_gpu
from .errorhandler import assert_msg_critical


class Diis:
    """
    Implements direct inversion of the iterative subspace.

    Instance variables
        - error_vectors: The list of error vectors.
    """

    def __init__(self, max_err_vecs, diis_thresh, scf_type):
        """
        Initializes iterative subspace by setting list of error vectors to
        empty list.
        """

        self.error_vectors = []

        self.fock_matrices = []

        self.max_err_vecs = max_err_vecs
        self.diis_thresh = diis_thresh
        self.scf_type = scf_type

        self.b_matrix = np.zeros((max_err_vecs, max_err_vecs))

    def clear(self):

        self.error_vectors.clear()

        self.fock_matrices.clear()

        self.max_err_vecs = None
        self.diis_thresh = None

        self.b_matrix = None

    def store_diis_data(self, fock_mat, e_mat, e_grad):

        if e_grad < self.diis_thresh:

            if len(self.error_vectors) == self.max_err_vecs:
                self.error_vectors.pop(0)
                self.fock_matrices.pop(0)
                sub_bmat = self.b_matrix[1:, 1:].copy()
                self.b_matrix[:-1, :-1] = sub_bmat[:, :]

            self.error_vectors.append(e_mat.copy())
            self.fock_matrices.append([x.copy() for x in fock_mat])

            n_vecs = len(self.error_vectors)
            for i in range(n_vecs):
                fij = dot_product_gpu(self.error_vectors[i],
                                      self.error_vectors[n_vecs - 1])
                self.b_matrix[i, n_vecs - 1] = fij
                self.b_matrix[n_vecs - 1, i] = fij

    def get_effective_fock(self, fock_mat):

        n_vecs = len(self.error_vectors)

        if n_vecs == 0:
            return fock_mat

        elif n_vecs == 1:
            return self.fock_matrices[0]

        else:
            weights = self.compute_weights()

            if self.scf_type == 'restricted':
                fock_matrices_a = [m[0] for m in self.fock_matrices]
                effmat_a = weighted_sum_gpu(weights, fock_matrices_a)
                # Note: return a tuple
                return (effmat_a,)

            elif self.scf_type == 'unrestricted':
                fock_matrices_a = [m[0] for m in self.fock_matrices]
                fock_matrices_b = [m[1] for m in self.fock_matrices]
                effmat_a = weighted_sum_gpu(weights, fock_matrices_a)
                effmat_b = weighted_sum_gpu(weights, fock_matrices_b)
                return (effmat_a, effmat_b)

            else:
                assert_msg_critical(
                    False,
                    'DIIS: restricted open-shell not yet supported')

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

        return np.linalg.solve(bmat, bvec)[:n_vecs]
