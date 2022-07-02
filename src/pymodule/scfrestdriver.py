#
#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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

from mpi4py import MPI
import numpy as np
import sys

from .veloxchemlib import mpi_master
from .veloxchemlib import MolecularOrbitals
from .veloxchemlib import molorb
from .veloxchemlib import fockmat
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
            smat = ovl_mat.to_numpy()
            tmat = oao_mat.to_numpy()

            dmat = den_mat.alpha_to_numpy(0)
            fmat = fock_mat.to_numpy(0)

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
            diff_mat = den_mat.sub(old_den_mat)
            ddmat = diff_mat.alpha_to_numpy(0)

            diff_den = np.linalg.norm(ddmat)
        else:
            diff_den = 0.0

        diff_den = self.comm.bcast(diff_den, root=mpi_master())

        return diff_den

    def _store_diis_data(self, i, fock_mat, den_mat):
        """
        Stores spin restricted closed shell Fock/Kohn-Sham and density matrices
        for current iteration. Overloaded base class method.

        :param i:
            The number of current SCF iteration.
        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param den_mat:
            The density matrix.
        """

        if self.rank == mpi_master():

            if not self._skip_iter:

                if len(self._fock_matrices) == self.max_err_vecs:

                    self._fock_matrices.popleft()
                    self._den_matrices.popleft()

                self._fock_matrices.append(fock_mat.alpha_to_numpy(0))
                self._den_matrices.append(den_mat.alpha_to_numpy(0))

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

            if len(self._fock_matrices) == 1:
                return np.copy(self._fock_matrices[0])

            if len(self._fock_matrices) > 1:
                acc_diis = CTwoDiis()
                acc_diis.compute_error_vectors(self._fock_matrices,
                                               self._den_matrices, ovl_mat,
                                               oao_mat)
                weights = acc_diis.compute_weights()

                return self._get_scaled_fock(weights)

            return fock_mat.alpha_to_numpy(0)

        return None

    def _get_scaled_fock(self, weights):
        """
        Computes effective spin restricted closed shell Fock/Kohn-Sham matrix
        by summing Fock/Kohn-Sham matrices scalwd with weigths.

        :param weights:
            The weights of Fock/Kohn-Sham matrices.

        :return:
            The scaled Fock/Kohn-Sham matrix.
        """

        effmat = np.zeros(self._fock_matrices[0].shape, dtype=float)

        for w, fmat in zip(weights, self._fock_matrices):
            effmat = effmat + w * fmat

        return effmat

    def _gen_molecular_orbitals(self, molecule, fock_mat, oao_mat):
        """
        Generates spin restricted molecular orbital by diagonalizing
        spin restricted closed shell Fock/Kohn-Sham matrix. Overloaded base
        class method.

        :param molecule:
            The molecule.
        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param oao_mat:
            The orthogonalization matrix.

        :return:
            The molecular orbitals.
        """

        if self.rank == mpi_master():
            tmat = oao_mat.to_numpy()
            fmo = np.matmul(tmat.T, np.matmul(fock_mat, tmat))

            eigs, evecs = np.linalg.eigh(fmo)
            orb_coefs = np.matmul(tmat, evecs)
            orb_coefs, eigs = self._delete_mos(orb_coefs, eigs)

            occa = molecule.get_aufbau_occupation(eigs.size, 'restricted')

            return MolecularOrbitals([orb_coefs], [eigs], [occa], molorb.rest)

        return MolecularOrbitals()

    def get_scf_type(self):
        """
        Gets string for spin restricted closed shell SCF calculation.
        Overloaded base class method.

        :return:
            The string for spin restricted closed shell SCF calculation.
        """

        pe_type = " with PE" if self._pe else ""

        if self._dft:
            return "Spin-Restricted Kohn-Sham" + pe_type

        return "Spin-Restricted Hartree-Fock" + pe_type

    def _update_fock_type(self, fock_mat):
        """
        Updates Fock matrix to fit selected functional in Kohn-Sham
        calculations.

        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        """

        if self.xcfun.is_hybrid():
            fock_mat.set_fock_type(fockmat.restjkx, 0)
            fock_mat.set_scale_factor(self.xcfun.get_frac_exact_exchange(), 0)
        else:
            fock_mat.set_fock_type(fockmat.restj, 0)

        return
