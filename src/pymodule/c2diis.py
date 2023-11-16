#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
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

import numpy as np


class CTwoDiis:
    """
    Implements direct inversion of the iterative subspace in C2 form proposed
    by H. Seller.

    Instance variables
        - error_vectors: The list of error vectors.
    """

    def __init__(self):
        """
        Initializes iterative subspace by setting list of error vectors to
        empty list.
        """

        self.error_vectors = []

    def compute_error_vectors_restricted(self, fock_matrices, density_matrices,
                                         overlap_matrix, oao_matrix):
        """
        Computes error vectors for list of AO Fock matrices using (FDS - SDF)
        in orthogonal AO basis.

        :param fock_matrices:
            The list of AO Fock matrices.
        :param density_matrices:
            The list of AO density matrices.
        :param overlap_matrix:
            The overlap matrix.
        :param oao_matrix:
            The orthogonalization matrix.
        """

        smat = overlap_matrix
        tmat = oao_matrix.to_numpy()

        self.error_vectors.clear()

        for fmat, dmat in zip(fock_matrices, density_matrices):

            fds = np.matmul(fmat, np.matmul(dmat, smat))

            self.error_vectors.append(
                np.matmul(tmat.T, np.matmul(fds - fds.T, tmat)))

    def compute_error_vectors_restricted_openshell(self, fock_matrices,
                                                   fock_matrices_beta,
                                                   density_matrices,
                                                   density_matrices_beta,
                                                   overlap_matrix, oao_matrix):
        """
        Computes error vectors for list of AO Fock matrices using (FDS - SDF)
        in orthogonal AO basis.

        :param fock_matrices:
            The list of AO Fock matrices (alpha spin).
        :param fock_matrices_beta:
            The list of AO Fock matrices (beta spin).
        :param density_matrices:
            The list of AO density matrices (alpha spin).
        :param density_matrices_beta:
            The list of AO density matrices (beta spin).
        :param overlap_matrix:
            The overlap matrix.
        :param oao_matrix:
            The orthogonalization matrix.
        """

        smat = overlap_matrix
        tmat = oao_matrix.to_numpy()

        self.error_vectors.clear()

        for fmat_a, fmat_b, dmat_a, dmat_b in zip(fock_matrices,
                                                  fock_matrices_beta,
                                                  density_matrices,
                                                  density_matrices_beta):

            fds_a = np.matmul(fmat_a, np.matmul(dmat_a, smat))
            fds_b = np.matmul(fmat_b, np.matmul(dmat_b, smat))

            err_a = np.matmul(tmat.T, np.matmul(fds_a - fds_a.T, tmat))
            err_b = np.matmul(tmat.T, np.matmul(fds_b - fds_b.T, tmat))

            err_tot = err_a + err_b
            self.error_vectors.append(err_tot)

    def compute_error_vectors_unrestricted(self, fock_matrices,
                                           fock_matrices_beta, density_matrices,
                                           density_matrices_beta,
                                           overlap_matrix, oao_matrix):
        """
        Computes error vectors for list of AO Fock matrices using (FDS - SDF)
        in orthogonal AO basis.

        :param fock_matrices:
            The list of AO Fock matrices (alpha spin).
        :param fock_matrices_beta:
            The list of AO Fock matrices (beta spin).
        :param density_matrices:
            The list of AO density matrices (alpha spin).
        :param density_matrices_beta:
            The list of AO density matrices (beta spin).
        :param overlap_matrix:
            The overlap matrix.
        :param oao_matrix:
            The orthogonalization matrix.
        """

        smat = overlap_matrix
        tmat = oao_matrix.to_numpy()

        self.error_vectors.clear()

        for fmat_a, fmat_b, dmat_a, dmat_b in zip(fock_matrices,
                                                  fock_matrices_beta,
                                                  density_matrices,
                                                  density_matrices_beta):

            fds_a = np.matmul(fmat_a, np.matmul(dmat_a, smat))
            fds_b = np.matmul(fmat_b, np.matmul(dmat_b, smat))

            err_a = np.matmul(tmat.T, np.matmul(fds_a - fds_a.T, tmat))
            err_b = np.matmul(tmat.T, np.matmul(fds_b - fds_b.T, tmat))

            err_tot = np.vstack((err_a, err_b))
            self.error_vectors.append(err_tot)

    def compute_weights(self):
        """
        Computes C2-DIIS weights from error vectors using H. Sellers method
        (Int. J. Quantum Chem., vol. 45, pp. 31-41, 1993).

        :return:
            The array of C2-DIIS weights with smallest residual error.
        """

        bmat = self._comp_bmatrix()

        beigs, bvecs = np.linalg.eigh(bmat)

        weights = self._pick_weights(self._norm_bvectors(bvecs))

        return weights

    def _comp_bmatrix(self):
        """
        Computes B-matrix of C2-DIIS method using error vectors.

        :return:
            The B-matrix.
        """

        bdim = len(self.error_vectors)

        bmat = np.zeros(shape=(bdim, bdim), dtype='float64')
        bidx = np.triu_indices(bdim)

        for i, j in zip(bidx[0], bidx[1]):

            fij = np.vdot(self.error_vectors[i], self.error_vectors[j])

            bmat[i][j] = fij
            bmat[j][i] = fij

        return bmat

    def _norm_bvectors(self, bvectors):
        """
        Normalizes B-matrix eigenvectors by rescaling them to 1.0.

        :param bvectors:
            The array of B-matrix eigenvectors.

        :return:
            The normalized B-matrix eigenvectors.
        """

        sum_vecs = np.sum(bvectors, axis=0)

        norm_vecs = []

        for i in range(len(sum_vecs)):
            if abs(sum_vecs[i]) > 1.0e-6:
                norm_vecs.append(bvectors[:, i] / sum_vecs[i])

        return norm_vecs

    def _pick_weights(self, weights):
        """
        Picks normalize B-matrix eigenvector with smallest residual error by
        computing residual error for all eigenvectors of B_matrix.

        :param bvectors:
            The array of B-matrix eigenvectors.

        :return:
            The normalized B-matrix eigenvector.
        """

        fmin = 1.0e+8
        wmin = weights[0]

        for w in weights:
            evec = np.zeros(self.error_vectors[0].shape, dtype='float64')

            for f, v in zip(w, self.error_vectors):
                evec = evec + f * v

            fact = np.vdot(evec, evec)

            if fmin > fact:
                fmin = fact
                wmin = w

        return wmin
