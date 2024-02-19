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

from mpi4py import MPI
from copy import deepcopy
import numpy as np
import sys
import os

from .veloxchemlib import mpi_master
from .veloxchemlib import fock_t as fockmat
from .veloxchemlib import XCFunctional, MolecularGrid
from .veloxchemlib import matmul_gpu, eigh_gpu, dot_product_gpu
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
            fds = matmul_gpu(fock_mat, matmul_gpu(den_mat.alpha_to_numpy(0), ovl_mat))
            fds = fds - fds.T
            # TODO: allow transpose in matmul_gpu
            e_mat = matmul_gpu(oao_mat.T.copy(), matmul_gpu(fds, oao_mat))

            e_mat_shape = e_mat.shape
            e_grad = 2.0 * np.linalg.norm(e_mat)
            max_grad = np.max(np.abs(e_mat))
        else:
            e_mat_shape = None
            e_grad = None
            max_grad = None

        e_mat_shape = self.comm.bcast(e_mat_shape, root=mpi_master())
        e_grad = self.comm.bcast(e_grad, root=mpi_master())
        max_grad = self.comm.bcast(max_grad, root=mpi_master())

        if self.rank != mpi_master():
            e_mat = np.zeros(e_mat_shape)
        self.comm.Bcast(e_mat, root=mpi_master())

        return e_mat, e_grad, max_grad

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
            dmat = den_mat.alpha_to_numpy(0)
            old_dmat = old_den_mat.alpha_to_numpy(0)
            ddmat = dmat - old_dmat

            #diff_den = np.linalg.norm(ddmat)
            diff_den = np.sqrt(dot_product_gpu(ddmat, ddmat))
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

            self._fock_matrices_alpha.append(fock_mat)
            self._density_matrices_alpha.append(den_mat.alpha_to_numpy(0))

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
                return np.copy(self._fock_matrices_alpha[0])

            if len(self._fock_matrices_alpha) > 1:
                acc_diis = CTwoDiis()
                acc_diis.compute_error_vectors_restricted(
                    self._fock_matrices_alpha, self._density_matrices_alpha,
                    ovl_mat, oao_mat)
                weights = acc_diis.compute_weights()

                return self._get_scaled_fock(weights)

            return fock_mat

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

        effmat = np.zeros(self._fock_matrices_alpha[0].shape)

        for w, fmat in zip(weights, self._fock_matrices_alpha):
            effmat += w * fmat

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
            tmat = oao_mat
            #fmo = np.matmul(tmat.T, np.matmul(fock_mat, tmat))
            fmo = matmul_gpu(tmat.T.copy(), matmul_gpu(fock_mat, tmat))

            # TODO: double check the criteria
            if max(fmo.shape) < 4096:
                eigs, evecs = np.linalg.eigh(fmo)
            else:
                eigs, evecs = eigh_gpu(fmo)
                evecs = evecs.T.copy()

            #orb_coefs = np.matmul(tmat, evecs)
            orb_coefs = matmul_gpu(tmat, evecs)
            orb_coefs, eigs = self._delete_mos(orb_coefs, eigs)

            occa = molecule.get_aufbau_alpha_occupation(eigs.size)

            return MolecularOrbitals([orb_coefs], [eigs], [occa], molorb.rest)

        return MolecularOrbitals()

    def get_scf_type_str(self):
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
