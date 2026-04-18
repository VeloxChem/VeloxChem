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

from .mathutils import symmetric_matrix_function


class OrbitalLocalization:
    """
    Orbital localization driver

    Methods:
        - pipek_mezey(projector="lowdin" | "mulliken")
        - boys(dipole_integrals)
        - edmiston_ruedenberg(eri)  # placeholder
    """

    def __init__(self):
        """
        Initializes the orbital localization driver
        """

        self.max_iter = 100
        self.thresh = 1e-6

        # TODO: use ostream
        self.silent = False

    def _rotate(self, i, j, theta, method="boys"):
        """
        Rotate MO pair around theta
        """

        c = np.cos(theta)
        s = np.sin(theta)

        # rotate MOs
        Ci = self.C[:, i].copy()
        Cj = self.C[:, j].copy()

        self.C[:, i] = c * Ci + s * Cj
        self.C[:, j] = -s * Ci + c * Cj

        if method == "pm":
            # rotate transformed MOs
            Ci = self.C_eff[:, i].copy()
            Cj = self.C_eff[:, j].copy()

            self.C_eff[:, i] = c * Ci + s * Cj
            self.C_eff[:, j] = -s * Ci + c * Cj

        elif method == "boys":
            # rotate dipole integrals
            ri = self.r[:, i, :].copy()
            rj = self.r[:, j, :].copy()

            self.r[:, i, :] = c * ri + s * rj
            self.r[:, j, :] = -s * ri + c * rj

            ri = self.r[:, :, i].copy()
            rj = self.r[:, :, j].copy()

            self.r[:, :, i] = c * ri + s * rj
            self.r[:, :, j] = -s * ri + c * rj

        else:
            raise ValueError(f"Unknown method {method} provided in _rotate")

    def _build_atom_projector(self):
        """
        Compute AO to atom projector
        """
        n_atoms = len(self.atom_indices)
        n_basis = self.C.shape[0]

        P = np.zeros((n_atoms, n_basis))

        for a, mask in enumerate(self.atom_indices):
            P[a, mask] = 1.0

        self.P_atom = P

    def _pm_pair(self, i, j):
        """
        Compute prerequesites for analytical Jacobian
        """
        Ci = self.C_eff[:, i]
        Cj = self.C_eff[:, j]

        return (np.matmul(self.P_atom,
                          (Ci * Ci)), np.matmul(self.P_atom, (Cj * Cj)),
                np.matmul(self.P_atom, (Ci * Cj)))

    def _pm_optimal_theta(self, Qi, Qj, Pij):
        """
        Approximate analytical angle for PM rotation
        (neglecting quartic terms)

        DOI: 10.1021/ct401016x
        """

        dQ = Qi - Qj

        K2 = np.sum(dQ * dQ - 4 * Pij * Pij)
        L2 = 4 * np.sum(dQ * Pij)

        if abs(K2) < 1e-14 and abs(L2) < 1e-14:
            return 0.0

        return 0.5 * np.arctan2(L2, K2)

    def _boys_optimal_theta(self, ri, rj, dij):
        """
        Exact boys rotation angle

        DOI: 10.1063/1.1681683
        """

        rij = ri - rj

        g = 2 * np.dot(rij, dij)
        h = np.dot(rij, rij) - 4 * np.dot(dij, dij)

        if abs(g) < 1e-14 and abs(h) < 1e-14:
            return 0.0

        return 0.25 * np.arctan2(g, h)

    def pipek_mezey(self, C, S, atom_map, projector="lowdin"):
        """
        Pipek–Mezey orbital localization

        Parameters
        ----------
        C : np.ndarray
            MO coefficients

        projector : str
            "lowdin"
            "mulliken"

        S : np.ndarray, optional
            Overlap matrix

        atom_map : list[int]
            AO -> atom mapping

        Returns
        -------
        C : np.ndarray
            Localized MO coefficients
        """

        self.C = np.array(C, copy=True)
        self.nao, self.norb = self.C.shape
        self.S = np.array(S)

        # transform MOs according to projector
        eigs = np.linalg.eigvalsh(self.S)
        min_eig = np.min(eigs)
        # TODO: Double-check whether -1.0e-10 is the right cutoff for
        # treating negative overlap eigenvalues as a hard error here.
        if min_eig < -1.0e-10:
            raise ValueError("Overlap matrix is not positive semidefinite; "
                             f"minimum eigenvalue = {min_eig:.3e}")
        if projector == "lowdin":
            # TODO: Double-check whether 1.0e-8 is the right screening
            # threshold for the overlap matrix function in PM localization.
            self.X = symmetric_matrix_function(self.S,
                                               lambda x: 1.0 / np.sqrt(x),
                                               thresh=1.0e-8)
        elif projector == "mulliken":
            self.X = symmetric_matrix_function(self.S, np.sqrt, thresh=1.0e-8)
        else:
            raise NotImplementedError(
                f"Requested projector {projector} not implemented")

        self.C_eff = np.matmul(self.X.T, self.C)

        self.atom_map = np.array(atom_map)
        self.n_atoms = np.max(atom_map) + 1

        # masks for AO to atom map
        self.atom_indices = [
            np.where(self.atom_map == A)[0] for A in range(self.n_atoms)
        ]

        self._build_atom_projector()

        for it in range(self.max_iter):

            delta = 0.0

            for i in range(self.norb):
                for j in range(i + 1, self.norb):

                    Qi, Qj, Pij = self._pm_pair(i, j)

                    theta = self._pm_optimal_theta(Qi, Qj, Pij)

                    if abs(theta) < 1e-12:
                        continue

                    # safeguard
                    theta *= 0.5

                    self._rotate(i, j, theta, method="pm")

                    delta += abs(theta)

            if delta < self.thresh:
                if not self.silent:
                    print(f"PM converged after {it:3d}  iterations")
                break

            if it == self.max_iter - 1:
                # return the object anyway, since unlike for SCF,
                # reaching the convergence threshold is not necessarily
                # required, as the physics are unchanged.
                if not self.silent:
                    print(f"PM only converged to delta = {delta:.6e} , "
                          f"instead of {self.thresh:.2e}")

        return self.C

    def boys(self, C, dipole_integrals):
        """
        Boys orbital localization

        Parameters
        ----------
        C : np.ndarray
            MO coefficients

        dipole integrals : np.ndarray, tuple(np.ndarray), list(np.ndarray)

        Returns
        -------
        C : np.ndarray
            Localized MO coefficients
        """

        self.C = np.array(C, copy=True)
        self.nao, self.norb = self.C.shape
        mu = dipole_integrals

        self.r = np.array(
            [np.matmul(self.C.T, np.matmul(mu[k], self.C)) for k in range(3)])

        for it in range(self.max_iter):

            delta = 0.0

            for i in range(self.norb):
                for j in range(i + 1, self.norb):

                    ri = self.r[:, i, i]
                    rj = self.r[:, j, j]
                    dij = self.r[:, i, j]

                    theta = self._boys_optimal_theta(ri, rj, dij)

                    if abs(theta) < 1e-12:
                        continue

                    self._rotate(i, j, theta, method="boys")

                    delta += abs(theta)

            if delta < self.thresh:
                if not self.silent:
                    print(f"Boys converged after {it:3d}  iterations")
                break

            if it == self.max_iter - 1:
                # return the object anyway, since unlike for SCF,
                # reaching the convergence threshold is not necessarily
                # required, as the physics are unchanged.
                if not self.silent:
                    print(f"Boys only converged to delta = {delta:.6e} , "
                          f"instead of {self.thresh:.2e}")

        return self.C

    def edmiston_ruedenberg(self, eri):
        """
        Placeholder for ER localization scheme
        """

        raise NotImplementedError(
            "ER localization requires AO->MO 2e integral transformation " +
            "and is not implemented.")
