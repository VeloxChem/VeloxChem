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
import sys

from .mathutils import symmetric_matrix_function
from .oneeints import compute_electric_dipole_integrals, compute_overlap_integrals
from .molecularorbitals import MolecularOrbitals, molorb
from .visualizationdriver import VisualizationDriver
from .outputstream import OutputStream


class OrbitalLocalizationDriver:
    """
    Implements Pipek–Mezey and Boys orbital localization.

    :param ostream:
        Output stream.

    Instance variables
        - max_iter: Maximum number of Jacobi sweeps.
        - thresh: Convergence threshold for rotation angle sum.
        - method: Localization scheme ("pm" or "boys").
        - pm_projector: Population projector for PM ("mulliken" or "lowdin").
    """

    def __init__(self, ostream=None):
        """
        Initializes the orbital localization driver
        """

        if ostream is None:
            ostream = OutputStream(sys.stdout)
        self.ostream = ostream

        self.max_iter = 100
        self.thresh = 1e-6

        self.method = "pm"
        self.pm_projector = "lowdin"

        # Working arrays populated by localization routines
        self.C = None
        self.SC = None
        self.C_eff = None
        self.r = None
        self.P_atom = None
        self.atom_indices = None

    def _rotate(self, i, j, theta, method="boys"):
        """
        Rotates MO pair (i, j) by angle theta and updates associated
        auxiliary arrays.

        :param i:
            Index of first MO.
        :param j:
            Index of second MO.
        :param theta:
            Rotation angle in radians.
        :param method:
            Localization method ("pm" or "boys"); determines which auxiliary
            arrays are rotated alongside the MO coefficients.
        """

        c = np.cos(theta)
        s = np.sin(theta)

        # rotate MOs
        Ci = self.C[:, i].copy()
        Cj = self.C[:, j].copy()

        self.C[:, i] = c * Ci + s * Cj
        self.C[:, j] = -s * Ci + c * Cj

        if method == "pm":
            if self.pm_projector == "lowdin":
                # rotate Lowdin-transformed MOs
                Ci = self.C_eff[:, i].copy()
                Cj = self.C_eff[:, j].copy()

                self.C_eff[:, i] = c * Ci + s * Cj
                self.C_eff[:, j] = -s * Ci + c * Cj
            elif self.pm_projector == "mulliken":
                # rotate overlap-weighted MOs
                SCi = self.SC[:, i].copy()
                SCj = self.SC[:, j].copy()

                self.SC[:, i] = c * SCi + s * SCj
                self.SC[:, j] = -s * SCi + c * SCj

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
        Builds the AO-to-atom projector from atom indices.

        Sets self.P_atom to a (n_atoms, n_basis) matrix where each row
        corresponds to one atom and entries are 1.0 for AOs belonging to
        that atom.
        """
        n_atoms = len(self.atom_indices)
        n_basis = self.C.shape[0]

        P = np.zeros((n_atoms, n_basis))

        for a, mask in enumerate(self.atom_indices):
            P[a, mask] = 1.0

        self.P_atom = P

    def _pm_pair(self, i, j):
        """
        Computes per-atom populations and transition population for
        MO pair (i, j).

        :param i:
            Index of first MO.
        :param j:
            Index of second MO.
        :return:
            Tuple (Qi, Qj, Pij) where Qi/Qj are per-atom population
            vectors and Pij is the symmetrized transition population.
        """
        if self.pm_projector == "lowdin":
            Ci = self.C_eff[:, i]
            Cj = self.C_eff[:, j]

            return (
                np.matmul(self.P_atom, (Ci * Ci)),
                np.matmul(self.P_atom, (Cj * Cj)),
                np.matmul(self.P_atom, (Ci * Cj)),
            )

        elif self.pm_projector == "mulliken":
            Ci = self.C[:, i]
            Cj = self.C[:, j]

            SCi = self.SC[:, i]
            SCj = self.SC[:, j]

            # Symmetrized transition population keeps the pair-rotation
            # formulas in _pm_optimal_theta in the same 2*sin*cos*Pij form.
            Pij = 0.5 * np.matmul(self.P_atom, (Ci * SCj + Cj * SCi))

            return (
                np.matmul(self.P_atom, (Ci * SCi)),
                np.matmul(self.P_atom, (Cj * SCj)),
                Pij,
            )

        else:
            raise NotImplementedError(
                f"Requested projector {self.pm_projector} not implemented")

    def _pm_optimal_theta(self, Qi, Qj, Pij):
        """
        Returns the approximate optimal Pipek–Mezey rotation angle,
        neglecting quartic terms (DOI: 10.1021/ct401016x).

        :param Qi:
            Per-atom population vector for MO i.
        :param Qj:
            Per-atom population vector for MO j.
        :param Pij:
            Symmetrized transition population vector.
        :return:
            Rotation angle in radians.
        """

        dQ = Qi - Qj

        K2 = np.sum(dQ * dQ - 4 * Pij * Pij)
        L2 = 4 * np.sum(dQ * Pij)

        if abs(K2) < 1e-14 and abs(L2) < 1e-14:
            return 0.0

        return 0.5 * np.arctan2(L2, K2)

    def _boys_optimal_theta(self, ri, rj, dij):
        """
        Returns the exact Boys rotation angle
        (DOI: 10.1063/1.1681683).

        :param ri:
            Dipole expectation value vector for MO i (length 3).
        :param rj:
            Dipole expectation value vector for MO j (length 3).
        :param dij:
            Transition dipole vector between MOs i and j (length 3).
        :return:
            Rotation angle in radians.
        """

        rij = ri - rj

        g = 2 * np.dot(rij, dij)
        h = np.dot(rij, rij) - 4 * np.dot(dij, dij)

        if abs(g) < 1e-14 and abs(h) < 1e-14:
            return 0.0

        return 0.25 * np.arctan2(g, h)

    def _mol_orbs_wrapper(self, loc_orbs, scf_res, mo_range):
        """
        Wraps localized orbitals into a MolecularOrbitals container.

        :param loc_orbs:
            List of localized orbital coefficient matrices (one per spin).
        :param scf_res:
            SCF result dictionary providing orbital shapes and occupations.
        :param mo_range:
            One-based inclusive range (start, end) of localized MOs,
            or None if all MOs were localized.
        :return:
            MolecularOrbitals instance with localized coefficients and zero
            orbital energies.
        """

        # get 0-based indices from 1-based inclusive mo_range;
        # when mo_range is None all MOs were localized
        if mo_range is None:
            mo_start = 0
            mo_end = scf_res["C_alpha"].shape[1]
        else:
            mo_start = mo_range[0] - 1
            mo_end = mo_range[1]

        zero_energies = np.zeros(scf_res["E_alpha"].shape)[mo_start:mo_end]

        C_alpha_loc = np.zeros(scf_res["C_alpha"].shape)
        C_alpha_loc[:, mo_start:mo_end] = loc_orbs[0][:, :]
        if len(loc_orbs) == 2:
            C_beta_loc = np.zeros(scf_res["C_beta"].shape)
            C_beta_loc[:, mo_start:mo_end] = loc_orbs[1][:, :]

        if scf_res["scf_type"] == "restricted":
            return MolecularOrbitals(
                orbs=[C_alpha_loc],
                enes=[zero_energies],
                occs=[scf_res["occ_alpha"]],
                orbs_type=molorb.rest,
            )

        elif scf_res["scf_type"] == "unrestricted":
            return MolecularOrbitals(
                orbs=[C_alpha_loc, C_beta_loc],
                enes=[zero_energies, zero_energies],
                occs=[scf_res["occ_alpha"], scf_res["occ_beta"]],
                orbs_type=molorb.unrest,
            )

        elif scf_res["scf_type"] == "restricted_openshell":
            return MolecularOrbitals(
                orbs=[C_alpha_loc],
                enes=[zero_energies],
                occs=[scf_res["occ_alpha"], scf_res["occ_beta"]],
                orbs_type=molorb.restopen,
            )

        else:
            raise ValueError(f"scf type {scf_res['scf_type']} unknown")

    def pipek_mezey(self, molecule, basis, mos, mo_range=None):
        """
        Performs Pipek–Mezey orbital localization using Jacobi sweeps.
        Uses the projector specified by self.pm_projector.

        :param molecule:
            Molecule object.
        :param basis:
            AO basis set.
        :param mos:
            Molecular orbital coefficient matrix to localize.
        :param mo_range:
            One-based inclusive range (start, end) of MOs to localize,
            or None to localize all MOs.
        :return:
            Localized MO coefficient matrix.
        """

        # localize only user defined MOs, using one-based inclusive mo_range
        if mo_range:
            self.C = mos[:, mo_range[0]-1:mo_range[1]].copy()
        else:
            self.C = mos.copy()

        norb = self.C.shape[1]
        S = compute_overlap_integrals(molecule, basis)

        # check the validity of the overlap matrix
        eigs = np.linalg.eigvalsh(S)
        min_eig = np.min(eigs)
        if min_eig < -1.0e-10:
            raise ValueError("Overlap matrix is not positive semidefinite; "
                             f"minimum eigenvalue = {min_eig:.3e}")

        # transform MOs according to projector
        self.C_eff = None
        self.SC = None
        if self.pm_projector == "lowdin":
            sqrt_S = symmetric_matrix_function(S, np.sqrt, thresh=1.0e-12)
            self.C_eff = np.matmul(sqrt_S, self.C)
        elif self.pm_projector == "mulliken":
            self.SC = np.matmul(S, self.C)
        else:
            raise NotImplementedError(
                f"Requested projector {self.pm_projector} not implemented")

        # Map atoms to AOs
        vis_drv = VisualizationDriver()
        self.atom_indices = vis_drv.map_atom_to_atomic_orbitals(molecule, basis)

        self._build_atom_projector()

        for it in range(self.max_iter):

            delta = 0.0

            for i in range(norb):
                for j in range(i + 1, norb):

                    Qi, Qj, Pij = self._pm_pair(i, j)

                    theta = self._pm_optimal_theta(Qi, Qj, Pij)

                    if abs(theta) < 1e-12:
                        continue

                    # safeguard
                    theta *= 0.5

                    self._rotate(i, j, theta, method="pm")

                    delta += abs(theta)

            if delta < self.thresh:
                self.ostream.print_info(
                    f"PM converged after {it:3d}  iterations")
                self.ostream.flush()
                break

            if it == self.max_iter - 1:
                # return the object anyway, since unlike for SCF,
                # reaching the convergence threshold is not necessarily
                # required, as the physics are unchanged.
                self.ostream.print_info(
                    f"PM only converged to delta = {delta:.6e} , "
                    f"instead of {self.thresh:.2e}")
                self.ostream.flush()

        return self.C.copy()

    def boys(self, molecule, basis, mos, mo_range=None):
        """
        Performs Boys (Foster–Boys) orbital localization using Jacobi sweeps.

        :param molecule:
            Molecule object.
        :param basis:
            AO basis set.
        :param mos:
            Molecular orbital coefficient matrix to localize.
        :param mo_range:
            One-based inclusive range (start, end) of MOs to localize,
            or None to localize all MOs.
        :return:
            Localized MO coefficient matrix.
        """

        # localize only user defined MOs, using one-based inclusive mo_range
        if mo_range:
            self.C = mos[:, mo_range[0]-1:mo_range[1]].copy()
        else:
            self.C = mos.copy()

        norb = self.C.shape[1]
        S = compute_overlap_integrals(molecule, basis)

        # check the validity of the overlap matrix
        eigs = np.linalg.eigvalsh(S)
        min_eig = np.min(eigs)
        if min_eig < -1.0e-10:
            raise ValueError("Overlap matrix is not positive semidefinite; "
                             f"minimum eigenvalue = {min_eig:.3e}")

        mu = compute_electric_dipole_integrals(molecule, basis)

        self.r = np.array(
            [np.matmul(self.C.T, np.matmul(mu[k], self.C)) for k in range(3)])

        for it in range(self.max_iter):

            delta = 0.0

            for i in range(norb):
                for j in range(i + 1, norb):

                    ri = self.r[:, i, i]
                    rj = self.r[:, j, j]
                    dij = self.r[:, i, j]

                    theta = self._boys_optimal_theta(ri, rj, dij)

                    if abs(theta) < 1e-12:
                        continue

                    self._rotate(i, j, theta, method="boys")

                    delta += abs(theta)

            if delta < self.thresh:
                self.ostream.print_info(
                    f"Boys converged after {it:3d}  iterations")
                self.ostream.flush()
                break

            if it == self.max_iter - 1:
                # return the object anyway, since unlike for SCF,
                # reaching the convergence threshold is not necessarily
                # required, as the physics are unchanged.
                self.ostream.print_info(
                    f"Boys only converged to delta = {delta:.6e} , "
                    f"instead of {self.thresh:.2e}")
                self.ostream.flush()

        return self.C.copy()

    def edmiston_ruedenberg(self, molecule, basis, scf_res, mo_range=None):
        """
        Placeholder for ER localization scheme
        """

        raise NotImplementedError(
            "ER localization requires AO->MO 2e integral transformation " +
            "and is not implemented.")

    def compute(self, molecule, basis, scf_res, mo_range=None):
        """
        Top-level entry point for orbital localization.

        Dispatches to the method indicated by self.method, handles
        unrestricted spin cases, and wraps results in a
        MolecularOrbitals container.

        :param molecule:
            Molecule object.
        :param basis:
            AO basis set.
        :param scf_res:
            SCF result dictionary with keys "C_alpha", "C_beta",
            "E_alpha", "occ_alpha", "occ_beta", "scf_type".
        :param mo_range:
            One-based inclusive range (start, end) of MOs to localize,
            or None to localize all MOs.
        :return:
            Dictionary {"loc_orbs": MolecularOrbitals instance}.
        """

        if self.method == "boys":
            compute_func = self.boys
        elif self.method == "pm":
            compute_func = self.pipek_mezey
            if self.pm_projector not in ["mulliken", "lowdin"]:
                raise NotImplementedError(
                    f"only mulliken and lowdin projectors are currently "
                    f"implemented, not {self.pm_projector}")
        else:
            raise NotImplementedError(
                f"only boys and pm are currently supported, "
                f"not {self.method}")

        if mo_range:
            if len(mo_range) != 2:
                raise ValueError("mo_range object should be of length 2")
            if mo_range[0] == 0:
                raise ValueError("mo_range starts counting from 1")

        if scf_res["scf_type"] == "unrestricted":
            # localize alpha and beta independently
            alpha = compute_func(
                molecule, basis, scf_res["C_alpha"],
                mo_range=mo_range,
            )
            beta = compute_func(
                molecule, basis, scf_res["C_beta"],
                mo_range=mo_range,
            )
            return {
                "loc_orbs": self._mol_orbs_wrapper([alpha, beta], scf_res, mo_range),
            }
        else:
            alpha = compute_func(
                molecule, basis, scf_res["C_alpha"],
                mo_range=mo_range,
            )
            return {
                "loc_orbs": self._mol_orbs_wrapper([alpha], scf_res, mo_range),
            }
