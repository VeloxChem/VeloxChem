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
import time as tm
import sys
import os

from .veloxchemlib import mpi_master
from .veloxchemlib import fock_t as fockmat
from .veloxchemlib import XCFunctional, MolecularGrid
from .veloxchemlib import matmul_gpu, eigh_gpu, dot_product_gpu
from .veloxchemlib import compute_error_vector_gpu
from .veloxchemlib import transform_matrix_gpu
from .molecularorbitals import MolecularOrbitals, molorb
from .outputstream import OutputStream
from .scfdriver import ScfDriver


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
            t0 = tm.time()
            e_mat = compute_error_vector_gpu(oao_mat, fock_mat, den_mat.alpha_to_numpy(0), ovl_mat)
            t1 = tm.time()
            self.ostream.print_info(f'    X^T(FDS-SDF)X: {t1-t0:.2f} sec')

            e_mat_shape = e_mat.shape
            #e_grad = 2.0 * np.linalg.norm(e_mat)
            e_grad = 2.0 * np.sqrt(dot_product_gpu(e_mat, e_mat))
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

    def _gen_molecular_orbitals(self, molecule, fock_mat, oao_mat,
                                num_gpus_per_node):
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
        :param num_gpus_per_node:
            The number of GPUs per MPI process.

        :return:
            The molecular orbitals.
        """

        if self.rank == mpi_master():
            t0 = tm.time()
            fmo = transform_matrix_gpu(oao_mat, fock_mat)
            t1 = tm.time()
            eigs, evecs_T = eigh_gpu(fmo, num_gpus_per_node)
            evecs = evecs_T.T
            t2 = tm.time()
            orb_coefs = matmul_gpu(oao_mat, evecs)
            t3 = tm.time()
            orb_coefs, eigs = self._delete_mos(orb_coefs, eigs)
            t4 = tm.time()
            self.ostream.print_info(f'    X^T(F)X      : {t1-t0:.2f} sec')
            self.ostream.print_info(f'    eigh         : {t2-t1:.2f} sec')
            self.ostream.print_info(f'    MOs          : {t3-t2:.2f} sec')
            self.ostream.print_info(f'    deleteMOs    : {t4-t3:.2f} sec')

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
