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
from copy import deepcopy
import numpy as np
import math
import sys

from .veloxchemlib import XCFunctional, MolecularGrid
from .veloxchemlib import mpi_master, boltzmann_in_hartreeperkelvin
from .molecularorbitals import MolecularOrbitals, molorb
from .outputstream import OutputStream
from .scfdriver import ScfDriver
from .c2diis import CTwoDiis


class ScfUnrestrictedDriver(ScfDriver):
    """
    Implements spin unrestricted open shell SCF method with C2-DIIS and
    two-level C2-DIIS convergence accelerators.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes spin unrestricted open shell SCF driver to default setup
        (convergence threshold, initial guess, etc) by calling base class
        constructor.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        super().__init__(comm, ostream)

        self._scf_type = 'unrestricted'

    def _comp_gradient(self, fock_mat, ovl_mat, den_mat, oao_mat):
        """
        Computes spin unrestricted open shell electronic gradient using
        Fock/Kohn-Sham matrix. Overloaded base class method.

        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param ovl_mat:
            The overlap matrix.
        :param den_mat:
            The density matrix.
        :param oao_mat:
            The orthogonalization matrix.

        :return:
            The electronic gradient.
        """

        if self.rank == mpi_master():
            smat = ovl_mat
            tmat = oao_mat

            dmat_a = den_mat[0]
            dmat_b = den_mat[1]

            fmat_a = fock_mat[0]
            fmat_b = fock_mat[1]

            fds_a = np.matmul(fmat_a, np.matmul(dmat_a, smat))
            fds_b = np.matmul(fmat_b, np.matmul(dmat_b, smat))

            e_mat_a = np.matmul(tmat.T, np.matmul(fds_a - fds_a.T, tmat))
            e_mat_b = np.matmul(tmat.T, np.matmul(fds_b - fds_b.T, tmat))

            e_grad = np.linalg.norm(e_mat_a) + np.linalg.norm(e_mat_b)
            max_grad = max(np.max(np.abs(e_mat_a)), np.max(np.abs(e_mat_b)))
        else:
            e_grad = 0.0
            max_grad = 0.0

        e_grad = self.comm.bcast(e_grad, root=mpi_master())
        max_grad = self.comm.bcast(max_grad, root=mpi_master())

        return e_grad, max_grad

    def _comp_density_change(self, den_mat, old_den_mat):
        """
        Computes norm of spin unrestricted open shell density change between
        two density matrices. Overloaded base class method.

        :param den_mat:
            The current density matrix.
        :param old_den_mat:
            The previous density matrix.

        :return:
            The norm of change between two density matrices.
        """

        if self.rank == mpi_master():
            ddmat_a = den_mat[0] - old_den_mat[0]
            ddmat_b = den_mat[1] - old_den_mat[1]

            diff_den_a = np.linalg.norm(ddmat_a)
            diff_den_b = np.linalg.norm(ddmat_b)

            diff_den = max(diff_den_a, diff_den_b)
        else:
            diff_den = 0.0

        diff_den = self.comm.bcast(diff_den, root=mpi_master())

        return diff_den

    def _store_diis_data(self, fock_mat, den_mat, ovl_mat, e_grad):
        """
        Stores spin unrestricted open shell Fock/Kohn-Sham and density matrices
        for current iteration. Overloaded base class method.

        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param den_mat:
            The density matrix.
        :param ovl_mat:
            The overlap matrix (used in ROSCF).
        :param e_grad:
            The electronic gradient.
        """

        if self.rank == mpi_master() and e_grad < self.diis_thresh:

            if len(self._fock_matrices_alpha) == self.max_err_vecs:
                self._fock_matrices_alpha.popleft()
                self._fock_matrices_beta.popleft()

                self._density_matrices_alpha.popleft()
                self._density_matrices_beta.popleft()

            self._fock_matrices_alpha.append(fock_mat[0].copy())
            self._fock_matrices_beta.append(fock_mat[1].copy())

            self._density_matrices_alpha.append(den_mat[0].copy())
            self._density_matrices_beta.append(den_mat[1].copy())

    def _get_effective_fock(self, fock_mat, ovl_mat, oao_mat):
        """
        Computes effective spin unrestricted open shell Fock/Kohn-Sham matrix
        in OAO basis by applying Lowdin or canonical orthogonalization to AO
        Fock/Kohn-Sham matrix. Overloaded base class method.

        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param ovl_mat:
            The overlap matrix.
        :param oao_mat:
            The orthogonalization matrix.

        :return:
            The effective Fock/Kohn-Sham matrices.
        """

        if self.rank == mpi_master():

            if len(self._fock_matrices_alpha) == 1:

                return (np.copy(self._fock_matrices_alpha[0]),
                        np.copy(self._fock_matrices_beta[0]))

            if len(self._fock_matrices_alpha) > 1:

                acc_diis = CTwoDiis()

                acc_diis.compute_error_vectors_unrestricted(
                    self._fock_matrices_alpha, self._fock_matrices_beta,
                    self._density_matrices_alpha, self._density_matrices_beta,
                    ovl_mat, oao_mat)

                weights = acc_diis.compute_weights()

                return self._get_scaled_fock(weights)

            return tuple(fock_mat)

        return (None, None)

    def _get_scaled_fock(self, weights):
        """
        Computes effective spin unrestricted open shell Fock/Kohn-Sham matrix
        by summing Fock/Kohn-Sham matrices scaled with weigths.

        :param weights:
            The weights of Fock/Kohn-Sham matrices.

        :return:
            The scaled Fock/Kohn-Sham matrices.
        """

        effmat_a = np.zeros(self._fock_matrices_alpha[0].shape)
        effmat_b = np.zeros(self._fock_matrices_beta[0].shape)

        for w, fa, fb in zip(weights, self._fock_matrices_alpha,
                             self._fock_matrices_beta):

            effmat_a += w * fa
            effmat_b += w * fb

        return (effmat_a, effmat_b)

    def _gen_molecular_orbitals(self, molecule, eff_fock_mat, oao_mat):
        """
        Generates spin unrestricted molecular orbital by diagonalizing
        spin unrestricted open shell Fock/Kohn-Sham matrix. Overloaded base
        class method.

        :param molecule:
            The molecule.
        :param eff_fock_mat:
            The effective Fock/Kohn-Sham matrix.
        :param oao_mat:
            The orthogonalization matrix.

        :return:
            The molecular orbitals.
        """

        if self.rank == mpi_master():
            tmat = oao_mat
            eigs_a, evecs_a = np.linalg.eigh(
                np.linalg.multi_dot([tmat.T, eff_fock_mat[0], tmat]))
            eigs_b, evecs_b = np.linalg.eigh(
                np.linalg.multi_dot([tmat.T, eff_fock_mat[1], tmat]))

            orb_coefs_a = np.matmul(tmat, evecs_a)
            orb_coefs_b = np.matmul(tmat, evecs_b)
            orb_coefs_a, eigs_a = self._delete_mos(orb_coefs_a, eigs_a)
            orb_coefs_b, eigs_b = self._delete_mos(orb_coefs_b, eigs_b)

            occa = molecule.get_aufbau_alpha_occupation(eigs_a.size)
            occb = molecule.get_aufbau_beta_occupation(eigs_b.size)

            if self.pfon and (self.pfon_temperature > 0):

                self.ostream.print_info(
                    f'Applying pseudo-FON (T={self.pfon_temperature:.0f}K)')

                kT = boltzmann_in_hartreeperkelvin() * self.pfon_temperature
                inv_kT = 1.0 / kT

                nocc_a = molecule.number_of_alpha_electrons()
                e_fermi_a = 0.5 * (eigs_a[nocc_a - 1] + eigs_a[nocc_a])
                idx_start_a = max(0, nocc_a - self.pfon_nocc)
                idx_end_a = min(eigs_a.size, nocc_a + self.pfon_nvir)
                pfon_a = {}
                sum_pfon_a = 0.0
                for idx in range(idx_start_a, idx_end_a):
                    try:
                        exp_ene_kT = math.exp(
                            (eigs_a[idx] - e_fermi_a) * inv_kT)
                    except OverflowError:
                        exp_ene_kT = float('inf')
                    pfon_a[idx] = 1.0 / (1.0 + exp_ene_kT)
                    sum_pfon_a += pfon_a[idx]
                pfon_scale_a = self.pfon_nocc / sum_pfon_a
                for idx in range(idx_start_a, idx_end_a):
                    pfon_a[idx] *= pfon_scale_a
                    occa[idx] = pfon_a[idx]

                nocc_b = molecule.number_of_beta_electrons()
                e_fermi_b = 0.5 * (eigs_b[nocc_b - 1] + eigs_b[nocc_b])
                idx_start_b = max(0, nocc_b - self.pfon_nocc)
                idx_end_b = min(eigs_b.size, nocc_b + self.pfon_nvir)
                pfon_b = {}
                sum_pfon_b = 0.0
                for idx in range(idx_start_b, idx_end_b):
                    try:
                        exp_ene_kT = math.exp(
                            (eigs_b[idx] - e_fermi_b) * inv_kT)
                    except OverflowError:
                        exp_ene_kT = float('inf')
                    pfon_b[idx] = 1.0 / (1.0 + exp_ene_kT)
                    sum_pfon_b += pfon_b[idx]
                pfon_scale_b = self.pfon_nocc / sum_pfon_b
                for idx in range(idx_start_b, idx_end_b):
                    pfon_b[idx] *= pfon_scale_b
                    occb[idx] = pfon_b[idx]

            return MolecularOrbitals([orb_coefs_a, orb_coefs_b],
                                     [eigs_a, eigs_b], [occa, occb],
                                     molorb.unrest)

        return MolecularOrbitals()

    def get_scf_type_str(self):
        """
        Gets string for spin unrestricted open shell SCF calculation.
        Overloaded base class method.

        :return:
            The string for spin unrestricted open shell SCF calculation.
        """

        if self.embedding is not None:
            emb_type = ' with ' + self.embedding['settings']['embedding_method']
        else:
            emb_type = ''

        if self._dft:
            return "Spin-Unrestricted Kohn-Sham" + emb_type

        return "Spin-Unrestricted Hartree-Fock" + emb_type

    def natural_orbitals(self, scf_tensors=None):
        """
        Compute the UHF natural orbitals

        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            The natural orbitals.
        """

        if scf_tensors is None:
            scf_tensors = self.scf_tensors

        if self.rank == mpi_master():
            # Get total density
            D_total = scf_tensors['D_alpha'] + scf_tensors['D_beta']

            # Get some MO coefficients and create C^-1
            C = scf_tensors['C_alpha']
            S = scf_tensors['S']
            C_inv = np.matmul(S, C)

            # Transform total density to MO basis
            D_MO = np.linalg.multi_dot([C_inv.T, D_total, C_inv])

            # Diagonalize
            occupations, eigenvectors = np.linalg.eigh(D_MO)

            # Create the final orbitals
            C_natural = np.matmul(C, eigenvectors)

            # Compute the orbital energy as expectation value of the averaged Fock
            # matrix (they are not eigenvalues!)
            F_alpha = scf_tensors['F_alpha']
            F_beta = scf_tensors['F_beta']
            F_avg = 0.5 * (F_alpha + F_beta)

            orbital_energies = np.diag(
                np.linalg.multi_dot([C_natural.T, F_avg, C_natural]))

            # Sort by orbital energies or by occupation numbers?
            # idx = orbital_energies.argsort() # Sort by orbital energies
            idx = occupations.argsort()[::-1]  # Sort by occupation numbers
            orbital_energies = orbital_energies[idx]
            occupations = occupations[idx]
            C_natural = C_natural[:, idx]

            # Create the MolecularOrbitals object and return
            natural_orbitals = MolecularOrbitals([C_natural],
                                                 [orbital_energies],
                                                 [occupations], molorb.rest)
        else:
            natural_orbitals = MolecularOrbitals()

        return natural_orbitals

    def __deepcopy__(self, memo):
        """
        Implements deepcopy.

        :param memo:
            The memo dictionary for deepcopy.

        :return:
            A deepcopy of self.
        """

        new_scf_drv = ScfUnrestrictedDriver(self.comm, self.ostream)

        for key, val in vars(self).items():
            if isinstance(val, (MPI.Intracomm, OutputStream)):
                pass
            elif isinstance(val, XCFunctional):
                new_scf_drv.key = XCFunctional(val)
            elif isinstance(val, MolecularGrid):
                new_scf_drv.key = MolecularGrid(val)
            else:
                new_scf_drv.key = deepcopy(val)

        return new_scf_drv
