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
from .veloxchemlib import matmul_gpu, eigh_gpu, dot_product_gpu
from .veloxchemlib import compute_error_vector_gpu
from .veloxchemlib import transform_matrix_gpu
from .veloxchemlib import mpi_master, boltzmann_in_hartreeperkelvin
from .molecularorbitals import MolecularOrbitals, molorb
from .outputstream import OutputStream
from .scfdriver import ScfDriver


class ScfRestrictedOpenDriver(ScfDriver):
    """
    Implements spin restricted open shell SCF method with DIIS and
    two-level DIIS convergence accelerators.

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
            dmat_a, dmat_b = den_mat
            fmat_a, fmat_b = fock_mat

            e_mat_a = compute_error_vector_gpu(oao_mat, fmat_a, dmat_a, ovl_mat)
            e_mat_b = compute_error_vector_gpu(oao_mat, fmat_b, dmat_b, ovl_mat)

            e_mat = e_mat_a + e_mat_b
            e_mat *= np.sqrt(2)

            e_mat_shape = e_mat.shape

            e_grad = np.sqrt(dot_product_gpu(e_mat, e_mat))
            max_grad = np.max(np.abs(e_mat))
        else:
            e_mat_shape = None
            e_grad = None
            max_grad = None

        e_mat_shape = self.comm.bcast(e_mat_shape, root=mpi_master())
        e_grad, max_grad = self.comm.bcast(
            (e_grad, max_grad), root=mpi_master())

        if self.rank != mpi_master():
            e_mat = np.zeros(e_mat_shape)
        self.comm.Bcast(e_mat, root=mpi_master())

        return e_mat, e_grad, max_grad

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

            diff_den_a = np.sqrt(dot_product_gpu(ddmat_a, ddmat_a))
            diff_den_b = np.sqrt(dot_product_gpu(ddmat_b, ddmat_b))

            diff_den = max(diff_den_a, diff_den_b)
        else:
            diff_den = 0.0

        diff_den = self.comm.bcast(diff_den, root=mpi_master())

        return diff_den

    def _gen_molecular_orbitals(self, molecule, eff_fock_mat, oao_mat,
                                num_gpus_per_node):
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
            eff_fmo = transform_matrix_gpu(oao_mat, eff_fock_mat[0])

            eigs, evecs_T = eigh_gpu(eff_fmo, num_gpus_per_node)
            evecs = evecs_T.T

            orb_coefs = matmul_gpu(oao_mat, evecs)
            orb_coefs, eigs = self._delete_mos(orb_coefs, eigs)

            occa = molecule.get_aufbau_alpha_occupation(eigs.size)
            occb = molecule.get_aufbau_beta_occupation(eigs.size)

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

                nocc_b = molecule.number_of_beta_electrons()
                e_fermi_b = 0.5 * (eigs[nocc_b - 1] + eigs[nocc_b])
                idx_start_b = max(0, nocc_b - self.pfon_nocc)
                idx_end_b = min(eigs.size, nocc_b + self.pfon_nvir)
                pfon_b = {}
                sum_pfon_b = 0.0
                for idx in range(idx_start_b, idx_end_b):
                    try:
                        exp_ene_kT = math.exp((eigs[idx] - e_fermi_b) * inv_kT)
                    except OverflowError:
                        exp_ene_kT = float('inf')
                    pfon_b[idx] = 1.0 / (1.0 + exp_ene_kT)
                    sum_pfon_b += pfon_b[idx]
                pfon_scale_b = self.pfon_nocc / sum_pfon_b
                for idx in range(idx_start_b, idx_end_b):
                    pfon_b[idx] *= pfon_scale_b
                    occb[idx] = pfon_b[idx]

            return MolecularOrbitals([orb_coefs], [eigs], [occa, occb],
                                     molorb.restopen)

        return MolecularOrbitals()

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
