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

    def comp_gradient(self, fock_mat, ovl_mat, den_mat, oao_mat):
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

            e_grad = np.linalg.norm(e_mat_a) + np.linalg.norm(e_mat_b)
            max_grad = max(np.max(np.abs(e_mat_a)), np.max(np.abs(e_mat_b)))
        else:
            e_grad = 0.0
            max_grad = 0.0

        e_grad = self.comm.bcast(e_grad, root=mpi_master())
        max_grad = self.comm.bcast(max_grad, root=mpi_master())

        return e_grad, max_grad

    def comp_density_change(self, den_mat, old_den_mat):
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

    def store_diis_data(self, i, fock_mat, den_mat):
        """
        Stores spin unrestricted open shell Fock/Kohn-Sham and density matrices
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

                self._fock_matrices.append(fock_mat.alpha_to_numpy(0))
                self._den_matrices.append(den_mat.alpha_to_numpy(0))

                self._fock_matrices_beta.append(fock_mat.beta_to_numpy(0))
                self._den_matrices_beta.append(den_mat.beta_to_numpy(0))

    def get_effective_fock(self, fock_mat, ovl_mat, oao_mat):
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

            if len(self._fock_matrices) == 1:

                return (np.copy(self._fock_matrices[0]),
                        np.copy(self._fock_matrices_beta[0]))

            if len(self._fock_matrices) > 1:

                acc_diis = CTwoDiis()

                acc_diis.compute_unrestricted_error_vectors(
                    self._fock_matrices, self._fock_matrices_beta,
                    self._den_matrices, self._den_matrices_beta, ovl_mat,
                    oao_mat)

                weights = acc_diis.compute_weights()

                return self.get_scaled_fock(weights)

            return fock_mat.alpha_to_numpy(0), fock_mat.beta_to_numpy(0)

        return None, None

    def get_scaled_fock(self, weights):
        """
        Computes effective spin unrestricted open shell Fock/Kohn-Sham matrix
        by summing Fock/Kohn-Sham matrices scalwd with weigths.

        :param weights:
            The weights of Fock/Kohn-Sham matrices.

        :return:
            The scaled Fock/Kohn-Sham matrices.
        """

        effmat_a = np.zeros(self._fock_matrices[0].shape)
        effmat_b = np.zeros(self._fock_matrices_beta[0].shape)

        for w, fa, fb in zip(weights, self._fock_matrices,
                             self._fock_matrices_beta):

            effmat_a += w * fa
            effmat_b += w * fb

        return effmat_a, effmat_b

    def gen_molecular_orbitals(self, molecule, fock_mat, oao_mat):
        """
        Generates spin unrestricted molecular orbital by diagonalizing
        spin unrestricted open shell Fock/Kohn-Sham matrix. Overloaded base
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

            fmo_a = np.matmul(tmat.T, np.matmul(fock_mat[0], tmat))
            fmo_b = np.matmul(tmat.T, np.matmul(fock_mat[1], tmat))

            eigs_a, evecs_a = np.linalg.eigh(fmo_a)
            eigs_b, evecs_b = np.linalg.eigh(fmo_b)

            orb_coefs_a = np.matmul(tmat, evecs_a)
            orb_coefs_b = np.matmul(tmat, evecs_b)

            orb_coefs_a, eigs_a = self.delete_mos(orb_coefs_a, eigs_a)
            orb_coefs_b, eigs_b = self.delete_mos(orb_coefs_b, eigs_b)

            occa, occb = molecule.get_aufbau_occupation(eigs_a.size,
                                                        'unrestricted')

            return MolecularOrbitals([orb_coefs_a, orb_coefs_b],
                                     [eigs_a, eigs_b], [occa, occb],
                                     molorb.unrest)

        return MolecularOrbitals()

    def get_scf_type(self):
        """
        Gets string for spin unrestricted open shell SCF calculation.
        Overloaded base class method.

        :return:
            The string for spin unrestricted open shell SCF calculation.
        """

        pe_type = " with PE" if self._pe else ""

        if self._dft:
            return "Spin-Unrestricted Kohn-Sham" + pe_type

        return "Spin-Unrestricted Hartree-Fock" + pe_type

    def update_fock_type(self, fock_mat):
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

    def natural_orbitals(self):
        """
        Compute the UHF natural orbitals

        :return:
            The natural orbitals.
        """

        # Get total density
        D_total = self.scf_tensors['D_alpha'] + self.scf_tensors['D_beta']

        # Get some MO coefficients and create C^-1
        C = self.scf_tensors['C_alpha']
        S = self.scf_tensors['S']
        C_inv = np.matmul(S, C)

        # Transform total density to MO basis
        D_MO = np.linalg.multi_dot([C_inv.T, D_total, C_inv])

        # Diagonalize
        occupations, eigenvectors = np.linalg.eigh(D_MO)

        # Create the final orbitals
        C_natural = np.matmul(C, eigenvectors)

        # Compute the orbital energy as expectation value of the averaged Fock
        # matrix (they are not eigenvalues!)
        F_alpha = self.scf_tensors['F_alpha']
        F_beta = self.scf_tensors['F_beta']
        F_avg = 0.5 * (F_alpha + F_beta)

        orbital_energies = np.diag(
            np.linalg.multi_dot([C_natural.T, F_avg, C_natural]))

        # Sort by orbital energies or by occupation numbers?
        # idx = orbital_energies.argsort() # Sort by orbital energies
        idx = occupations.argsort()[::-1]  # Sort by occupation numbers
        orbital_energies = orbital_energies[idx]
        occupations = occupations[idx]
        # eigenvectors = eigenvectors[:,idx]
        C_natural = C_natural[:, idx]

        # Create the MolecularOrbitals object and return
        natural_orbitals = MolecularOrbitals([C_natural], [orbital_energies],
                                             [occupations], molorb.rest)

        return natural_orbitals
