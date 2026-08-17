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

from collections import deque

import numpy as np

from .diis import compute_pulay_weights


class KDiis:
    """
    Kollmar DIIS orbital updater.

    The Pulay error is the occupied-virtual orbital-gradient block represented
    in a fixed orthonormal AO frame. The same weights extrapolate the AO Fock
    matrices, after which first-order perturbation theory and a unitary orbital
    rotation produce the next iterate without diagonalizing the full Fock
    matrix.
    """

    def __init__(self,
                 max_vectors=10,
                 min_denominator=1.0e-4,
                 max_rotation=0.2):
        """
        Initializes the KDIIS history and orbital-step safeguards.

        :param max_vectors:
            The maximum number of Fock matrices and gradients retained.
        :param min_denominator:
            The minimum absolute occupied-virtual energy denominator.
        :param max_rotation:
            The maximum absolute first-order rotation amplitude.
        """

        self.max_vectors = max_vectors
        self.min_denominator = min_denominator
        self.max_rotation = max_rotation
        self.error_vectors = deque()
        self.fock_matrices = deque()

    def clear(self):
        """Clears the Fock-matrix and orbital-gradient history."""

        self.error_vectors.clear()
        self.fock_matrices.clear()

    def update_restricted(self,
                          fock_matrix,
                          coefficients,
                          number_occupied,
                          overlap_matrix,
                          oao_matrix,
                          level_shift=0.0):
        """
        Produces the next restricted set of molecular orbitals.

        :param overlap_matrix:
            The AO overlap matrix.
        :param oao_matrix:
            The orthogonalization matrix defining the fixed gradient frame.
        :return:
            A tuple containing updated coefficients and orbital energies.
        """

        return self._update((fock_matrix,), (coefficients,), (number_occupied,),
                            overlap_matrix, oao_matrix, level_shift)

    def update_unrestricted(self,
                            fock_matrices,
                            coefficients,
                            numbers_occupied,
                            overlap_matrix,
                            oao_matrix,
                            level_shift=0.0):
        """
        Produces the next alpha and beta sets of molecular orbitals.

        A single Pulay system is built from the stacked spin gradients so the
        same coefficients extrapolate both spin Fock matrices.

        :param overlap_matrix:
            The AO overlap matrix.
        :param oao_matrix:
            The orthogonalization matrix defining the fixed gradient frame.
        :return:
            Tuples containing updated coefficients and orbital energies.
        """

        return self._update(tuple(fock_matrices), tuple(coefficients),
                            tuple(numbers_occupied), overlap_matrix, oao_matrix,
                            level_shift)

    def _update(self, fock_matrices, coefficients, numbers_occupied,
                overlap_matrix, oao_matrix, level_shift):
        """Updates one or two spin channels using a shared Pulay solve."""

        gradients = []
        for fock, coeff, nocc in zip(fock_matrices, coefficients,
                                     numbers_occupied):
            fock_mo = np.linalg.multi_dot([coeff.T, fock, coeff])
            gradient_mo = fock_mo[nocc:, :nocc]

            # Gradients from different iterations must be compared in one
            # common coordinate system. Transform the occupied-virtual block
            # to the fixed orthonormal AO basis. This representation is
            # invariant to arbitrary rotations within the occupied and
            # virtual MO subspaces, including semicanonicalization.
            mo_to_oao = np.linalg.multi_dot(
                [oao_matrix.T, overlap_matrix, coeff])
            occupied = mo_to_oao[:, :nocc]
            virtual = mo_to_oao[:, nocc:]
            gradients.append(
                np.linalg.multi_dot([virtual, gradient_mo, occupied.T]))

        error_vector = np.concatenate(
            [gradient.reshape(-1) for gradient in gradients])
        self._append_history(fock_matrices, error_vector)

        weights = compute_pulay_weights(self.error_vectors, 'KDIIS')
        effective_fock = []
        for spin in range(len(fock_matrices)):
            eff = np.zeros_like(fock_matrices[spin])
            for weight, stored_fock in zip(weights, self.fock_matrices):
                eff += weight * stored_fock[spin]
            effective_fock.append(eff)

        new_coefficients = []
        new_energies = []
        for fock, coeff, nocc in zip(effective_fock, coefficients,
                                     numbers_occupied):
            new_coeff, energies = self._first_order_update(
                fock, coeff, nocc, level_shift)
            new_coefficients.append(new_coeff)
            new_energies.append(energies)

        return tuple(new_coefficients), tuple(new_energies)

    def _append_history(self, fock_matrices, error_vector):
        """Appends one synchronized Fock-matrix and gradient entry."""

        if len(self.error_vectors) == self.max_vectors:
            self.error_vectors.popleft()
            self.fock_matrices.popleft()

        self.error_vectors.append(error_vector)
        self.fock_matrices.append(tuple(fock.copy() for fock in fock_matrices))

    def _first_order_update(self, fock_matrix, coefficients, number_occupied,
                            level_shift):
        """Builds and applies the first-order occupied-virtual rotation."""

        nmo = coefficients.shape[1]
        nocc = number_occupied
        nvir = nmo - nocc

        if nocc == 0 or nvir == 0:
            energies = np.diag(
                np.linalg.multi_dot([coefficients.T, fock_matrix,
                                     coefficients])).copy()
            return coefficients.copy(), energies

        fock_mo = np.linalg.multi_dot(
            [coefficients.T, fock_matrix, coefficients])

        # Canonicalize only within the occupied and virtual spaces. These
        # rotations do not alter the density and make the perturbative energy
        # denominators well defined without a full Fock diagonalization.
        occ_energies, occ_rotation = np.linalg.eigh(fock_mo[:nocc, :nocc])
        vir_energies, vir_rotation = np.linalg.eigh(fock_mo[nocc:, nocc:])
        block_rotation = np.zeros((nmo, nmo), dtype='float64')
        block_rotation[:nocc, :nocc] = occ_rotation
        block_rotation[nocc:, nocc:] = vir_rotation

        semicanonical_coeff = coefficients @ block_rotation
        semicanonical_fock = np.linalg.multi_dot(
            [semicanonical_coeff.T, fock_matrix, semicanonical_coeff])
        gradient = semicanonical_fock[nocc:, :nocc]

        denominators = (vir_energies[:, None] - occ_energies[None, :] +
                        max(0.0, level_shift))
        signs = np.where(denominators < 0.0, -1.0, 1.0)
        safe_denominators = signs * np.maximum(np.abs(denominators),
                                               self.min_denominator)
        amplitudes = -gradient / safe_denominators

        largest = np.max(np.abs(amplitudes))
        if largest > self.max_rotation:
            amplitudes *= self.max_rotation / largest

        generator = np.zeros((nmo, nmo), dtype='float64')
        generator[nocc:, :nocc] = amplitudes
        generator[:nocc, nocc:] = -amplitudes.T
        rotation = self._exponential_skew_symmetric(generator)

        new_coefficients = semicanonical_coeff @ rotation
        rotated_fock = np.linalg.multi_dot(
            [new_coefficients.T, fock_matrix, new_coefficients])
        energies = np.diag(rotated_fock).copy()

        return new_coefficients, energies

    @staticmethod
    def _exponential_skew_symmetric(generator):
        """Computes a stable real exponential of a skew-symmetric matrix."""

        eigenvalues, eigenvectors = np.linalg.eigh(1.0j * generator)
        rotation = np.linalg.multi_dot([
            eigenvectors,
            np.diag(np.exp(-1.0j * eigenvalues)),
            eigenvectors.conj().T
        ])
        return np.real_if_close(rotation, tol=1000).real
