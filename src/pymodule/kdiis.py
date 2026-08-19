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

    The DIIS errors are first-order, energy-weighted orbital rotations. The
    same weights extrapolate accumulated orbital-rotation angles, and the next
    iterate is generated directly on the orbital manifold without extrapolating
    or diagonalizing a Fock matrix.

    Rotation operators are stored in a fixed orthonormal AO frame. This makes
    the history invariant to MO phases and to rotations within the occupied and
    virtual subspaces.
    """

    def __init__(self,
                 max_vectors=10,
                 min_denominator=1.0e-4,
                 max_rotation=0.2):
        """
        Initializes the KDIIS history and orbital-step safeguards.

        :param max_vectors:
            The maximum number of rotations and angles retained.
        :param min_denominator:
            The minimum absolute occupied-virtual energy denominator.
        :param max_rotation:
            The maximum absolute first-order rotation amplitude.
        """

        self.max_vectors = max_vectors
        self.min_denominator = min_denominator
        self.max_rotation = max_rotation

        self.error_vectors = deque()
        self.rotation_matrices = deque()
        self.angle_matrices = deque()

        self.current_angles = None
        self._last_occupied_projectors = None

        self.last_weights = None
        self.last_min_denominators = None
        self.last_max_rotations = None

    def clear(self):
        """
        Clears the stored rotations and accumulated orbital angles.
        """

        self.error_vectors.clear()
        self.rotation_matrices.clear()
        self.angle_matrices.clear()

        self.current_angles = None
        self._last_occupied_projectors = None

        self.last_weights = None
        self.last_min_denominators = None
        self.last_max_rotations = None

    def update_restricted(self,
                          fock_matrix,
                          coefficients,
                          number_occupied,
                          overlap_matrix,
                          oao_matrix,
                          level_shift=0.0):
        """
        Produces the next restricted set of molecular orbitals.

        :param fock_matrix:
            The current AO Fock/Kohn-Sham matrix.
        :param coefficients:
            The current AO molecular-orbital coefficient matrix, with occupied
            columns followed by virtual columns.
        :param number_occupied:
            The number of occupied molecular orbitals.
        :param overlap_matrix:
            The AO overlap matrix.
        :param oao_matrix:
            The orthogonalization matrix defining the fixed KDIIS frame.
        :param level_shift:
            The nonnegative shift added to occupied-virtual energy
            denominators.
        :return:
            A pair of one-element tuples containing the updated coefficient
            matrix and orbital-energy array, respectively.
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

        A single Pulay system is built from the stacked spin rotations so the
        same coefficients extrapolate both spin-channel orbital angles.

        :param fock_matrices:
            The current alpha- and beta-spin AO Fock/Kohn-Sham matrices.
        :param coefficients:
            The current alpha- and beta-spin AO molecular-orbital coefficient
            matrices, with occupied columns followed by virtual columns.
        :param numbers_occupied:
            The numbers of occupied alpha- and beta-spin molecular orbitals.
        :param overlap_matrix:
            The AO overlap matrix.
        :param oao_matrix:
            The orthogonalization matrix defining the fixed KDIIS frame.
        :param level_shift:
            The nonnegative shift added to occupied-virtual energy
            denominators.
        :return:
            A pair of two-element tuples containing the updated alpha- and
            beta-spin coefficient matrices and orbital-energy arrays,
            respectively.
        """

        return self._update(tuple(fock_matrices), tuple(coefficients),
                            tuple(numbers_occupied), overlap_matrix, oao_matrix,
                            level_shift)

    def _update(self, fock_matrices, coefficients, numbers_occupied,
                overlap_matrix, oao_matrix, level_shift):
        """
        Updates one or two spin channels using a shared KDIIS solve.
        """

        mo_to_oao = tuple(
            np.linalg.multi_dot([oao_matrix.T, overlap_matrix, coeff])
            for coeff in coefficients)

        self._reset_if_orbitals_changed(mo_to_oao, numbers_occupied)

        if self.current_angles is None:
            self.current_angles = tuple(
                np.zeros((orbitals.shape[0], orbitals.shape[0]),
                         dtype='float64') for orbitals in mo_to_oao)

        rotations = []
        min_denominators = []
        max_rotations = []
        for fock, coeff, nocc, orbitals in zip(fock_matrices, coefficients,
                                               numbers_occupied, mo_to_oao):
            rotation, min_denom, max_rot = self._build_rotation(
                fock, coeff, nocc, orbitals, level_shift)
            rotations.append(rotation)
            min_denominators.append(min_denom)
            max_rotations.append(max_rot)

        error_vector = np.concatenate(
            [rotation.reshape(-1) for rotation in rotations])
        self._append_history(tuple(rotations), self.current_angles,
                             error_vector)

        weights = compute_pulay_weights(self.error_vectors, 'KDIIS')
        self.last_weights = weights.copy()
        self.last_min_denominators = tuple(min_denominators)
        self.last_max_rotations = tuple(max_rotations)

        extrapolated_rotations = self._weighted_matrices(
            weights, self.rotation_matrices)
        extrapolated_angles = self._weighted_matrices(weights,
                                                      self.angle_matrices)

        new_coefficients = []
        new_energies = []
        new_angles = []
        step_was_limited = False

        for fock, coeff, orbitals, current_angle, mean_rotation, mean_angle in zip(
                fock_matrices, coefficients, mo_to_oao, self.current_angles,
                extrapolated_rotations, extrapolated_angles):
            target_angle = mean_angle + mean_rotation
            delta_oao = target_angle - current_angle

            delta_mo = np.linalg.multi_dot([orbitals.T, delta_oao, orbitals])
            delta_mo = 0.5 * (delta_mo - delta_mo.T)
            delta_mo, was_limited = self._limit_generator(delta_mo)
            step_was_limited = step_was_limited or was_limited

            orbital_rotation = self._exponential_skew_symmetric(delta_mo)
            new_coeff = np.matmul(coeff, orbital_rotation)
            new_coefficients.append(new_coeff)

            rotated_fock = np.linalg.multi_dot([new_coeff.T, fock, new_coeff])
            new_energies.append(np.diag(rotated_fock).copy())

            applied_delta = np.linalg.multi_dot(
                [orbitals, delta_mo, orbitals.T])
            new_angle = current_angle + applied_delta
            new_angles.append(0.5 * (new_angle - new_angle.T))

        self.current_angles = tuple(new_angles)
        self._last_occupied_projectors = self._occupied_projectors(
            tuple(
                np.linalg.multi_dot([oao_matrix.T, overlap_matrix, coeff])
                for coeff in new_coefficients), numbers_occupied)

        # A restricted extrapolated step indicates that the current DIIS model
        # is too aggressive. Preserve the accumulated angles but rebuild the
        # subspace from the next physical Fock/gradient pair.
        if step_was_limited:
            self.error_vectors.clear()
            self.rotation_matrices.clear()
            self.angle_matrices.clear()

        return tuple(new_coefficients), tuple(new_energies)

    def _build_rotation(self, fock_matrix, coefficients, number_occupied,
                        mo_to_oao, level_shift):
        """
        Builds an energy-weighted orbital rotation in the fixed OAO frame.
        """

        nmo = coefficients.shape[1]
        nocc = number_occupied
        nvir = nmo - nocc

        if nocc == 0 or nvir == 0:
            return np.zeros((nmo, nmo), dtype='float64'), float('inf'), 0.0

        fock_mo = np.linalg.multi_dot(
            [coefficients.T, fock_matrix, coefficients])

        # Semicanonicalization is used only to define a gauge-invariant diagonal
        # Hessian preconditioner. The orbitals themselves are not changed here.
        occ_energies, occ_rotation = np.linalg.eigh(fock_mo[:nocc, :nocc])
        vir_energies, vir_rotation = np.linalg.eigh(fock_mo[nocc:, nocc:])
        block_rotation = np.zeros((nmo, nmo), dtype='float64')
        block_rotation[:nocc, :nocc] = occ_rotation
        block_rotation[nocc:, nocc:] = vir_rotation

        semicanonical_fock = np.linalg.multi_dot(
            [block_rotation.T, fock_mo, block_rotation])
        gradient = semicanonical_fock[nocc:, :nocc]

        denominators = vir_energies[:, None] - occ_energies[None, :]
        safe_denominators = (
            np.maximum(np.abs(denominators), self.min_denominator) +
            max(0.0, level_shift))
        amplitudes = -gradient / safe_denominators
        largest_amplitude = np.max(np.abs(amplitudes))
        if largest_amplitude > self.max_rotation:
            amplitudes *= self.max_rotation / largest_amplitude

        generator_semicanonical = np.zeros((nmo, nmo), dtype='float64')
        generator_semicanonical[nocc:, :nocc] = amplitudes
        generator_semicanonical[:nocc, nocc:] = -amplitudes.T

        generator_mo = np.linalg.multi_dot(
            [block_rotation, generator_semicanonical, block_rotation.T])
        generator_oao = np.linalg.multi_dot(
            [mo_to_oao, generator_mo, mo_to_oao.T])
        generator_oao = 0.5 * (generator_oao - generator_oao.T)

        return (generator_oao, float(np.min(np.abs(denominators))),
                float(np.max(np.abs(amplitudes))))

    def _append_history(self, rotations, angles, error_vector):
        """
        Appends synchronized orbital rotations and accumulated angles.
        """

        if len(self.error_vectors) == self.max_vectors:
            self.error_vectors.popleft()
            self.rotation_matrices.popleft()
            self.angle_matrices.popleft()

        self.error_vectors.append(error_vector.copy())
        self.rotation_matrices.append(
            tuple(rotation.copy() for rotation in rotations))
        self.angle_matrices.append(tuple(angle.copy() for angle in angles))

    @staticmethod
    def _weighted_matrices(weights, history):
        """
        Extrapolates each spin-channel matrix using shared DIIS weights.
        """

        number_of_spins = len(history[0])
        result = []
        for spin in range(number_of_spins):
            matrix = np.zeros_like(history[0][spin])
            for weight, entry in zip(weights, history):
                matrix += weight * entry[spin]
            result.append(matrix)
        return tuple(result)

    def _limit_generator(self, generator):
        """
        Restricts an extrapolated orbital step to the configured trust radius.
        """

        largest = np.max(np.abs(generator)) if generator.size else 0.0
        if largest <= self.max_rotation:
            return generator, False
        return generator * (self.max_rotation / largest), True

    def _reset_if_orbitals_changed(self, mo_to_oao, numbers_occupied):
        """
        Resets stale history after an external occupied-subspace change.
        """

        if self._last_occupied_projectors is None:
            return

        projectors = self._occupied_projectors(mo_to_oao, numbers_occupied)
        changed = any(
            old.shape != new.shape or np.linalg.norm(old - new) > 1.0e-8
            for old, new in zip(self._last_occupied_projectors, projectors))
        if changed:
            self.clear()

    @staticmethod
    def _occupied_projectors(mo_to_oao, numbers_occupied):
        """
        Builds occupied-space projectors in the fixed orthonormal AO frame.
        """

        return tuple(orbitals[:, :nocc] @ orbitals[:, :nocc].T
                     for orbitals, nocc in zip(mo_to_oao, numbers_occupied))

    @staticmethod
    def _exponential_skew_symmetric(generator):
        """
        Computes a stable real exponential of a skew-symmetric matrix.
        """

        eigenvalues, eigenvectors = np.linalg.eigh(1.0j * generator)
        rotation = np.linalg.multi_dot([
            eigenvectors,
            np.diag(np.exp(-1.0j * eigenvalues)),
            eigenvectors.conj().T
        ])
        return np.real_if_close(rotation, tol=1000).real
