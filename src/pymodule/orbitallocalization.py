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


class OrbitalLocalization:
    """
    Orbital localization driver

    Methods:
        - pipek_mezey(projector="lowdin" | "mulliken")
        - boys(dipole_integrals)
        - edmiston_ruedenberg(eri)  # placeholder
    """

    def __init__(self, C):
        """
        Parameters
        ----------
        C : np.ndarray
            MO coefficients
        """

        self.C = np.array(C, copy=True)
        self.nao, self.norb = self.C.shape

    # ============================================================
    # Core rotation
    # ============================================================

    def _rotate(self, i, j, theta, projector=""):

        c = np.cos(theta)
        s = np.sin(theta)

        Ci = self.C[:, i].copy()
        Cj = self.C[:, j].copy()

        self.C[:, i] = c*Ci + s*Cj
        self.C[:, j] = -s*Ci + c*Cj

        if projector == "lowdin":
            Ci = self.C_ortho[:, i].copy()
            Cj = self.C_ortho[:, j].copy()

            self.C_ortho[:, i] = c*Ci + s*Cj
            self.C_ortho[:, j] = -s*Ci + c*Cj

    # ============================================================
    # ---------------- PIPEK–MEZEY -------------------------------
    # ============================================================

    def _pm_pair_lowdin(self, i, j):

        Ci = self.C_ortho[:, i]
        Cj = self.C_ortho[:, j]

        Qi = np.array([np.dot(Ci[idx], Ci[idx]) for idx in self.atom_indices])
        Qj = np.array([np.dot(Cj[idx], Cj[idx]) for idx in self.atom_indices])
        Pij = np.array([np.dot(Ci[idx], Cj[idx]) for idx in self.atom_indices])

        return Qi, Qj, Pij

    def _pm_pair_mulliken(self, i, j):

        Ci = self.C[:, i]
        Cj = self.C[:, j]

        Qi = np.zeros(self.n_atoms)
        Qj = np.zeros(self.n_atoms)
        Pij = np.zeros(self.n_atoms)

        for A, idx in enumerate(self.atom_indices):
            SAA = self.S[np.ix_(idx, idx)]

            CiA = Ci[idx]
            CjA = Cj[idx]

            Qi[A] = CiA @ SAA @ CiA
            Qj[A] = CjA @ SAA @ CjA
            Pij[A] = CiA @ SAA @ CjA

        return Qi, Qj, Pij

    def _pm_optimal_theta(self, Qi, Qj, Pij):
        # This is not the exact analytical expression,
        # but in order to obtain a compact form neglects
        # the "quartic" order terms
        # 10.1021/ct401016x

        dQ = Qi - Qj

        K2 = np.sum(dQ**2 - 4*Pij**2)
        L2 = 4*np.sum(dQ * Pij)

        if abs(K2) < 1e-14 and abs(L2) < 1e-14:
            return 0.0

        return 0.5 * np.arctan2(L2, K2)

    def pipek_mezey(self, S, atom_map, projector="lowdin", max_iter=50, tol=1e-6):
        """
        Pipek–Mezey orbital localization

        Parameters
        ----------
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

        pair_func = (
            self._pm_pair_lowdin
            if projector == "lowdin"
            else self._pm_pair_mulliken
        )

        self.S = np.array(S)
        self.atom_map = np.array(atom_map)
        self.n_atoms = np.max(atom_map) + 1

        if projector == "lowdin":
            # Löwdin orthogonalization
            eigs, U = np.linalg.eigh(self.S)
            self.X = U @ np.diag(1.0 / np.sqrt(eigs)) @ U.T

            self.C_ortho = self.X.T @ self.C

        # atom masks
        self.atom_indices = [np.where(self.atom_map == A)[0]
                             for A in range(self.n_atoms)]

        for it in range(max_iter):

            delta = 0.0

            for i in range(self.norb):
                for j in range(i+1, self.norb):

                    Qi, Qj, Pij = pair_func(i, j)

                    theta = self._pm_optimal_theta(Qi, Qj, Pij)

                    if abs(theta) < 1e-12:
                        continue

                    # loose monotonic safeguard
                    self._rotate(i, j, theta, projector=projector)

                    Qi2, Qj2, Pij2 = pair_func(i, j)
                    f_new = np.sum(Qi2**2 + Qj2**2)
                    f_old = np.sum(Qi**2 + Qj**2)

                    if f_new < f_old:
                        theta *= 0.5
                        self._rotate(i, j, -theta, projector=projector)

                    delta += abs(theta)

            print(f"PM Iter {it:3d}  delta = {delta:.6e}")

            if delta < tol:
                break

        return self.C

    # ============================================================
    # -------------------- BOYS ----------------------------------
    # ============================================================

    def boys(self, dipole_integrals, max_iter=50, tol=1e-6):
        """
        Boys orbital localization

        Parameters
        ----------
        dipole integrals : np.ndarray, tuple(np.ndarray), list(np.ndarray)

        Returns
        -------
        C : np.ndarray
            Localized MO coefficients
        """

        mu = dipole_integrals

        for it in range(max_iter):

            delta = 0.0

            for i in range(self.norb):
                for j in range(i+1, self.norb):

                    ri = np.array([
                        self.C[:, i] @ mu[k] @ self.C[:, i]
                        for k in range(3)
                    ])

                    rj = np.array([
                        self.C[:, j] @ mu[k] @ self.C[:, j]
                        for k in range(3)
                    ])

                    rij = np.array([
                        self.C[:, i] @ mu[k] @ self.C[:, j]
                        for k in range(3)
                    ])

                    A = np.dot(ri - rj, ri - rj)
                    B = 2*np.dot(ri - rj, rij)

                    if abs(B) < 1e-12:
                        continue

                    # exact analytical angle for Jacobian
                    # 10.1063/1.1681683
                    theta = 0.25 * np.arctan2(2*B, A)

                    self._rotate(i, j, theta)
                    delta += abs(theta)

            print(f"Boys Iter {it:3d}  delta = {delta:.6e}")

            if delta < tol:
                break

        return self.C

    # ============================================================
    # -------- EDMISTON–RUEDENBERG (STRUCTURE ONLY) --------------
    # ============================================================

    def edmiston_ruedenberg(self, eri, max_iter=20):
        """
        eri: (nao, nao, nao, nao)
        WARNING: placeholder (needs density fitting / AO→MO transform)
        """

        raise NotImplementedError(
            "ER localization requires AO->MO 2e integral transformation "
            "and is not implemented in this prototype."
        )
