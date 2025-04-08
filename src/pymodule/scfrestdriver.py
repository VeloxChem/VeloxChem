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


class ScfRestrictedDriver(ScfDriver):
    """
    Implements spin restricted closed shell SCF method with C2-DIIS and
    two-level C2-DIIS convergence accelerators.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes spin restricted closed shell SCF driver to default setup
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

        self._scf_type = 'restricted'

    def _comp_gradient(self, fock_mat, ovl_mat, den_mat, oao_mat):
        """
        Computes spin restricted closed shell electronic gradient using
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

            dmat = den_mat[0]
            fmat = fock_mat[0]

            fds = np.matmul(fmat, np.matmul(dmat, smat))

            e_mat = np.matmul(tmat.T, np.matmul(fds - fds.T, tmat))
            e_grad = 2.0 * np.linalg.norm(e_mat)
            max_grad = np.max(np.abs(e_mat))
        else:
            e_grad = 0.0
            max_grad = 0.0

        e_grad = self.comm.bcast(e_grad, root=mpi_master())
        max_grad = self.comm.bcast(max_grad, root=mpi_master())

        return e_grad, max_grad

    def _comp_density_change(self, den_mat, old_den_mat):
        """
        Computes norm of spin restricted closed shell density change between
        two density matrices. Overloaded base class method.

        :param den_mat:
            The current density matrix.
        :param old_den_mat:
            The previous density matrix.

        :return:
            The norm of change between two density matrices.
        """

        if self.rank == mpi_master():
            dmat = den_mat[0]
            old_dmat = old_den_mat[0]
            ddmat = dmat - old_dmat

            diff_den = np.linalg.norm(ddmat)
        else:
            diff_den = 0.0

        diff_den = self.comm.bcast(diff_den, root=mpi_master())

        return diff_den

    def _store_diis_data(self, fock_mat, den_mat, ovl_mat, e_grad):
        """
        Stores spin restricted closed shell Fock/Kohn-Sham and density matrices
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
                self._density_matrices_alpha.popleft()

            self._fock_matrices_alpha.append(fock_mat[0].copy())
            self._density_matrices_alpha.append(den_mat[0].copy())

    def _get_effective_fock(self, fock_mat, ovl_mat, oao_mat):
        """
        Computes effective spin restricted closed shell Fock/Kohn-Sham matrix
        in OAO basis by applying Lowdin or canonical orthogonalization to AO
        Fock/Kohn-Sham matrix. Overloaded base class method.

        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param ovl_mat:
            The overlap matrix.
        :param oao_mat:
            The orthogonalization matrix.

        :return:
            The effective Fock/Kohn-Sham matrix.
        """

        if self.rank == mpi_master():

            if len(self._fock_matrices_alpha) == 1:
                return (np.copy(self._fock_matrices_alpha[0]),)

            if len(self._fock_matrices_alpha) > 1:
                acc_diis = CTwoDiis()
                acc_diis.compute_error_vectors_restricted(
                    self._fock_matrices_alpha, self._density_matrices_alpha,
                    ovl_mat, oao_mat)
                weights = acc_diis.compute_weights()

                return self._get_scaled_fock(weights)

            return tuple(fock_mat)

        return (None,)

    def _get_scaled_fock(self, weights):
        """
        Computes effective spin restricted closed shell Fock/Kohn-Sham matrix
        by summing Fock/Kohn-Sham matrices scalwd with weigths.

        :param weights:
            The weights of Fock/Kohn-Sham matrices.

        :return:
            The scaled Fock/Kohn-Sham matrix.
        """

        effmat = np.zeros(self._fock_matrices_alpha[0].shape)

        for w, fmat in zip(weights, self._fock_matrices_alpha):
            effmat = effmat + w * fmat

        return (effmat,)

    def _gen_molecular_orbitals(self, molecule, eff_fock_mat, oao_mat):
        """
        Generates spin restricted molecular orbital by diagonalizing
        spin restricted closed shell Fock/Kohn-Sham matrix. Overloaded base
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
            eigs, evecs = np.linalg.eigh(
                np.linalg.multi_dot([tmat.T, eff_fock_mat[0], tmat]))

            orb_coefs = np.matmul(tmat, evecs)
            orb_coefs, eigs = self._delete_mos(orb_coefs, eigs)

            occa = molecule.get_aufbau_alpha_occupation(eigs.size)

            if self.pfon and (self.pfon_temperature > 0):

                self.ostream.print_info(
                    f'Applying pseudo-FON (T={self.pfon_temperature:.0f}K)')

                kT = boltzmann_in_hartreeperkelvin() * self.pfon_temperature
                inv_kT = 1.0 / kT

                nocc_a = molecule.number_of_alpha_electrons()
                e_fermi_a = 0.5 * (eigs[nocc_a - 1] + eigs[nocc_a])
                idx_start_a = max(0, nocc_a - self.pfon_nocc)
                idx_end_a = min(eigs.size, nocc_a + self.pfon_nvir)
                pfon_a = {}
                sum_pfon_a = 0.0
                for idx in range(idx_start_a, idx_end_a):
                    try:
                        exp_ene_kT = math.exp((eigs[idx] - e_fermi_a) * inv_kT)
                    except OverflowError:
                        exp_ene_kT = float('inf')
                    pfon_a[idx] = 1.0 / (1.0 + exp_ene_kT)
                    sum_pfon_a += pfon_a[idx]
                pfon_scale_a = self.pfon_nocc / sum_pfon_a
                for idx in range(idx_start_a, idx_end_a):
                    pfon_a[idx] *= pfon_scale_a
                    occa[idx] = pfon_a[idx]

            return MolecularOrbitals([orb_coefs], [eigs], [occa], molorb.rest)

        return MolecularOrbitals()

    def get_scf_type_str(self):
        """
        Gets string for spin restricted closed shell SCF calculation.
        Overloaded base class method.

        :return:
            The string for spin restricted closed shell SCF calculation.
        """

        if self.embedding is not None:
            emb_type = ' with ' + self.embedding['settings'][
                'embedding_method']
        else:
            emb_type = ''

        if self._dft:
            return "Spin-Restricted Kohn-Sham" + emb_type

        return "Spin-Restricted Hartree-Fock" + emb_type

    def __deepcopy__(self, memo):
        """
        Implements deepcopy.

        :param memo:
            The memo dictionary for deepcopy.

        :return:
            A deepcopy of self.
        """

        new_scf_drv = ScfRestrictedDriver(self.comm, self.ostream)

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
