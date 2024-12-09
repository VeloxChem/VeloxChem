#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

from .veloxchemlib import mpi_master
from .veloxchemlib import XCFunctional, MolecularGrid
from .molecularorbitals import MolecularOrbitals, molorb
from .outputstream import OutputStream
from .scfdriver import ScfDriver
from .c2diis import CTwoDiis


class ScfRestrictedOpenDriver(ScfDriver):
    """
    Implements spin restricted open shell SCF method with C2-DIIS and
    two-level C2-DIIS convergence accelerators.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes spin restricted open shell SCF driver to default setup
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

        self._scf_type = 'restricted_openshell'

    def _comp_gradient(self, fock_mat, ovl_mat, den_mat, oao_mat):
        """
        Computes spin restricted open shell electronic gradient using
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

            e_mat = e_mat_a + e_mat_b
            e_mat *= np.sqrt(2)

            e_grad = np.linalg.norm(e_mat)
            max_grad = np.max(np.abs(e_mat))
        else:
            e_grad = 0.0
            max_grad = 0.0

        e_grad = self.comm.bcast(e_grad, root=mpi_master())
        max_grad = self.comm.bcast(max_grad, root=mpi_master())

        return e_grad, max_grad

    def _comp_density_change(self, den_mat, old_den_mat):
        """
        Computes norm of spin restricted open shell density change between
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
        Stores spin restricted open shell Fock/Kohn-Sham and density matrices
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
                self._fock_matrices_proj.popleft()

                self._density_matrices_alpha.popleft()
                self._density_matrices_beta.popleft()

            self._fock_matrices_alpha.append(fock_mat[0].copy())
            self._fock_matrices_beta.append(fock_mat[1].copy())
            self._fock_matrices_proj.append(
                self.get_projected_fock(
                    fock_mat[0],
                    fock_mat[1],
                    den_mat[0],
                    den_mat[1],
                    ovl_mat,
                ))

            self._density_matrices_alpha.append(den_mat[0].copy())
            self._density_matrices_beta.append(den_mat[1].copy())

    def _get_effective_fock(self, fock_mat, ovl_mat, oao_mat):
        """
        Computes effective spin restricted open shell Fock/Kohn-Sham matrix
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
                return (self._fock_matrices_proj[0],)

            if len(self._fock_matrices_alpha) > 1:

                acc_diis = CTwoDiis()

                acc_diis.compute_error_vectors_restricted_openshell(
                    self._fock_matrices_alpha, self._fock_matrices_beta,
                    self._density_matrices_alpha, self._density_matrices_beta,
                    ovl_mat, oao_mat)

                weights = acc_diis.compute_weights()

                return self._get_scaled_fock(weights)

            return tuple(fock_mat)

        return (None,)

    def _get_scaled_fock(self, weights):
        """
        Computes effective spin restricted open shell Fock/Kohn-Sham matrix
        by summing Fock/Kohn-Sham matrices scalwd with weigths.

        :param weights:
            The weights of Fock/Kohn-Sham matrices.

        :return:
            The scaled Fock/Kohn-Sham matrix.
        """

        effmat = np.zeros(self._fock_matrices_proj[0].shape)

        for w, fmat in zip(weights, self._fock_matrices_proj):
            effmat = effmat + w * fmat

        return (effmat,)

    def _gen_molecular_orbitals(self, molecule, eff_fock_mat, oao_mat):
        """
        Generates spin restricted molecular orbital by diagonalizing
        spin restricted projected open shell Fock/Kohn-Sham matrix. Overloaded base
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
            occb = molecule.get_aufbau_beta_occupation(eigs.size)

            return MolecularOrbitals([orb_coefs], [eigs], [occa, occb],
                                     molorb.restopen)

        return MolecularOrbitals()

    def get_projected_fock(self, fa, fb, da, db, s):
        """
        Generates projected Fock matrix.

        :param fa:
            The Fock matrix of alpha spin.
        :param fb:
            The Fock matrix of beta spin.
        :param da:
            The density matrix of alpha spin.
        :param db:
            The density matrix of beta spin.
        :param s:
            The overlap matrix.

        :return:
            The projected Fock matrix.
        """

        naos = s.shape[0]

        inactive = np.matmul(s, db)
        active = np.matmul(s, da - db)
        virtual = np.eye(naos) - np.matmul(s, da)

        #       occ   act   vir
        #     +----------------+
        # occ | f0    fb    f0 |
        # act | fb    f0    fa |
        # vir | f0    fa    f0 |
        #     +----------------+

        f0 = 0.5 * (fa + fb)

        fcorr = np.linalg.multi_dot([inactive, fb - f0, active.T])
        fcorr += np.linalg.multi_dot([active, fa - f0, virtual.T])
        fcorr += fcorr.T

        return f0 + fcorr

    def get_scf_type_str(self):
        """
        Gets string for spin restricted open shell SCF calculation.
        Overloaded base class method.

        :return:
            The string for spin unrestricted open shell SCF calculation.
        """

        if self.embedding_options is not None:
            emb_type = ' with ' + self.embedding_options['settings'][
                'embedding_method']
        else:
            emb_type = ''

        if self._dft:
            return "Spin-Restricted Open-Shell Kohn-Sham" + emb_type

        return "Spin-Restricted Open-Shell Hartree-Fock" + emb_type

    def __deepcopy__(self, memo):
        """
        Implements deepcopy.

        :param memo:
            The memo dictionary for deepcopy.

        :return:
            A deepcopy of self.
        """

        new_scf_drv = ScfRestrictedOpenDriver(self.comm, self.ostream)

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
