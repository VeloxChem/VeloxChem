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
            smat = ovl_mat.to_numpy()
            tmat = oao_mat.to_numpy()

            dmat_a = den_mat.alpha_to_numpy(0)
            dmat_b = den_mat.beta_to_numpy(0)

            fmat_a = fock_mat.alpha_to_numpy(0)
            fmat_b = fock_mat.beta_to_numpy(0)

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
            diff_mat = den_mat.sub(old_den_mat)
            ddmat_a = diff_mat.alpha_to_numpy(0)
            ddmat_b = diff_mat.beta_to_numpy(0)

            diff_den_a = np.linalg.norm(ddmat_a)
            diff_den_b = np.linalg.norm(ddmat_b)

            diff_den = max(diff_den_a, diff_den_b)
        else:
            diff_den = 0.0

        diff_den = self.comm.bcast(diff_den, root=mpi_master())

        return diff_den

    def _store_diis_data(self, i, fock_mat, den_mat):
        """
        Stores spin restricted open shell Fock/Kohn-Sham and density matrices
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

                    self._fock_matrices_beta.popleft()
                    self._den_matrices_beta.popleft()

                    self._fock_matrices_proj.popleft()

                self._fock_matrices.append(fock_mat.alpha_to_numpy(0))
                self._den_matrices.append(den_mat.alpha_to_numpy(0))

                self._fock_matrices_beta.append(fock_mat.beta_to_numpy(0))
                self._den_matrices_beta.append(den_mat.beta_to_numpy(0))

                self._fock_matrices_proj.append(
                    self.get_projected_fock(
                        fock_mat.alpha_to_numpy(0),
                        fock_mat.beta_to_numpy(0),
                        den_mat.alpha_to_numpy(0),
                        den_mat.beta_to_numpy(0),
                        self.scf_tensors['S'],
                    ))

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

            if len(self._fock_matrices) == 1:
                return self._fock_matrices_proj[0]

            if len(self._fock_matrices) > 1:

                acc_diis = CTwoDiis()

                acc_diis.compute_error_vectors_restricted_openshell(
                    self._fock_matrices, self._fock_matrices_beta,
                    self._den_matrices, self._den_matrices_beta, ovl_mat,
                    oao_mat)

                weights = acc_diis.compute_weights()

                return self._get_scaled_fock(weights)

        return None, None

    def _get_scaled_fock(self, weights):
        """
        Computes effective spin restricted closed shell Fock/Kohn-Sham matrix
        by summing Fock/Kohn-Sham matrices scalwd with weigths.

        :param weights:
            The weights of Fock/Kohn-Sham matrices.

        :return:
            The scaled Fock/Kohn-Sham matrix.
        """

        effmat = np.zeros(self._fock_matrices_proj[0].shape, dtype=float)

        for w, fmat in zip(weights, self._fock_matrices_proj):
            effmat = effmat + w * fmat

        return effmat

    def _gen_molecular_orbitals(self, molecule, fock_mat, oao_mat):
        """
        Generates spin restricted molecular orbital by diagonalizing
        spin restricted projected open shell Fock/Kohn-Sham matrix. Overloaded base
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
            fock_mat = np.linalg.multi_dot([tmat.T, fock_mat, tmat])

            eigs, evecs = np.linalg.eigh(fock_mat)

            orb_coefs = np.linalg.multi_dot([tmat, evecs])
            orb_coefs, eigs = self._delete_mos(orb_coefs, eigs)

            occ = molecule.get_aufbau_occupation(eigs.size, 'restricted')

            return MolecularOrbitals([orb_coefs], [eigs], [occ], molorb.rest)

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

        pe_type = " with PE" if self._pe else ""

        if self._dft:
            return "Spin-Restricted Open-Shell Kohn-Sham" + pe_type

        return "Spin-Restricted Open-Shell Hartree-Fock" + pe_type

    def _update_fock_type(self, fock_mat):
        """
        Updates Fock matrix to fit selected functional in Kohn-Sham
        calculations.

        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        """

        if self.xcfun.is_hybrid():
            fock_mat.set_fock_type(fockmat.unrestjkx, 0, 'alpha')
            fock_mat.set_scale_factor(self.xcfun.get_frac_exact_exchange(), 0,
                                      'alpha')
            fock_mat.set_fock_type(fockmat.unrestjkx, 0, 'beta')
            fock_mat.set_scale_factor(self.xcfun.get_frac_exact_exchange(), 0,
                                      'beta')
        else:
            fock_mat.set_fock_type(fockmat.unrestj, 0, 'alpha')
            fock_mat.set_fock_type(fockmat.unrestj, 0, 'beta')

        return
