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
import numpy as np
import time
import sys

from .veloxchemlib import mpi_master
from .veloxchemlib import _RIFockDriver
from .veloxchemlib import TwoCenterElectronRepulsionDriver
from .veloxchemlib import SubMatrix
from .outputstream import OutputStream
from .molecularbasis import MolecularBasis
from .errorhandler import assert_msg_critical

try:
    from scipy.linalg import lu_factor, lu_solve, sqrtm
except ImportError:
    pass


class RIFockDriver:
    """
    Implements RI Fock driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes RI Fock driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream

        self._ri_drv = None

    def prepare_buffers(self,
                        molecule,
                        basis,
                        ri_auxiliary_basis,
                        k_metric=False,
                        verbose=True):
        """
        Prepare buffers for the RI Fock driver.
        """

        #assert_msg_critical(basis.get_label().lower().startswith('def2-'),
        #                    'RI Fock driver: Invalid basis set for RI-J')

        if verbose:
            self.ostream.print_info(
                'Using the resolution of the identity (RI) approximation.')
            self.ostream.print_blank()
            self.ostream.flush()

        if isinstance(ri_auxiliary_basis, str):
            if self.rank == mpi_master():
                basis_ri_j = MolecularBasis.read(molecule, ri_auxiliary_basis)
            else:
                basis_ri_j = None
            basis_ri_j = self.comm.bcast(basis_ri_j, root=mpi_master())
        else:
            basis_ri_j = MolecularBasis(ri_auxiliary_basis)

        if verbose:
            self.ostream.print_info('Dimension of RI auxiliary basis set ' +
                                    f'({basis_ri_j.get_label().upper()}): ' +
                                    f'{basis_ri_j.get_dimensions_of_basis()}')
            self.ostream.print_blank()
            self.ostream.flush()

        ri_prep_t0 = time.time()

        t2c_drv = TwoCenterElectronRepulsionDriver()
        mat_j = t2c_drv.compute(molecule, basis_ri_j)
        mat_j_np = mat_j.to_numpy()

        if verbose:
            self.ostream.print_info('Two-center integrals for RI done in ' +
                                    f'{time.time() - ri_prep_t0:.2f} sec.')
            self.ostream.print_blank()

        ri_prep_t0 = time.time()

        if 'scipy' in sys.modules:
            lu, piv = lu_factor(mat_j_np)
            inv_mat_j_np = lu_solve((lu, piv), np.eye(mat_j_np.shape[0]))
            if k_metric:
                inv_mat_k_np = sqrtm(inv_mat_j_np)
        else:
            inv_mat_j_np = np.linalg.inv(mat_j_np)
            if k_metric:
                assert_msg_critical(
                    False,
                    'RI Fock driver: K metric computation requires scipy module.'
                )

        if verbose:
            self.ostream.print_info('Metric(s) for RI done in ' +
                                    f'{time.time() - ri_prep_t0:.2f} sec.')
            self.ostream.print_blank()

        ri_prep_t0 = time.time()

        inv_mat_j = SubMatrix(
            [0, 0, inv_mat_j_np.shape[0], inv_mat_j_np.shape[1]])
        inv_mat_j.set_values(inv_mat_j_np)

        if k_metric:
            inv_mat_k = SubMatrix(
                [0, 0, inv_mat_k_np.shape[0], inv_mat_k_np.shape[1]])
            inv_mat_k.set_values(inv_mat_k_np)

        # note _RIFockDriver from C++
        if k_metric:
            self._ri_drv = _RIFockDriver(inv_mat_j, inv_mat_k)
        else:
            self._ri_drv = _RIFockDriver(inv_mat_j)

        local_atoms = molecule.partition_atoms(self.comm)
        self._ri_drv.prepare_buffers(molecule, basis, basis_ri_j, local_atoms)

        if verbose:
            self.ostream.print_info('Buffer preparation for RI done in ' +
                                    f'{time.time() - ri_prep_t0:.2f} sec.')
            self.ostream.print_blank()
            self.ostream.flush()

    def compute_bq_vector(self, den_mat_for_fock):
        """
        Compute bq vector using the RI Fock driver.
        """

        local_gvec = np.array(
            self._ri_drv.compute_local_bq_vector(den_mat_for_fock))

        gvec = np.zeros(local_gvec.shape)
        self.comm.Allreduce(local_gvec, gvec, op=MPI.SUM)

        return gvec

    def compute(self, den_mat_for_fock, fock_type):
        """
        Compute Fock matrix using the RI Fock driver.
        """

        assert_msg_critical(fock_type.lower() == 'j',
                            'RI Fock driver: Expecting "j" fock_type')

        gvec = self.compute_bq_vector(den_mat_for_fock)

        return self._ri_drv.local_compute(den_mat_for_fock, gvec, fock_type)

    def compute_mo_bq_vectors(self, lambda_p, lambda_h):
        """
        Compute bq vectors in MO basis using the RI Fock driver.
        """

        if self.comm.Get_size() == 1:
            return self._ri_drv.compute_bq_vector(lambda_p, lambda_h)
        else:
            return None
