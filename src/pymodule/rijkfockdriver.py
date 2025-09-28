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
from .veloxchemlib import _RIJKFockDriver
from .veloxchemlib import TwoCenterElectronRepulsionDriver
from .veloxchemlib import SubMatrix
from .outputstream import OutputStream
from .molecularbasis import MolecularBasis
from .errorhandler import assert_msg_critical

try:
    from scipy.linalg import lu_factor, lu_solve, sqrtm
except ImportError:
    pass


class RIJKFockDriver:
    """
    Implements RI-JK Fock driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes RI-JK Fock driver.
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
        
        self.metric = None

        self._ri_drv = _RIJKFockDriver()

    def compute_metric(self,
                       molecule,
                       basis,
                       verbose=True):
        """
        Computes Cholesky decomposed J metric for the RI-JK Fock driver.
        
        :param molecule:
            The molecule to compute J metric.
        :param basis:
            The auxilary basis to compute J metric.
        :param verbose:
            The information printout level.
        """
        
        # NOTE: For numerical stability J metric is computed on master node
        # and distributed to worker nodes
        if self.rank == mpi_master():

            if verbose:
                self.ostream.print_info(
                    'Using the resolution of the identity (RI) approximation.')
                self.ostream.print_blank()
                self.ostream.flush()

            if isinstance(basis, str):
                basis_ri = MolecularBasis.read(molecule, basis)
            else:
                basis_ri = MolecularBasis(basis)

            if verbose:
                self.ostream.print_info('Dimension of RI auxiliary basis set ' +
                                        f'({basis_ri.get_label().upper()}): ' +
                                        f'{basis_ri.get_dimensions_of_basis()}')
                self.ostream.print_blank()
                self.ostream.flush()

            ri_prep_t0 = time.time()

            t2c_drv = TwoCenterElectronRepulsionDriver()
            mat_j = t2c_drv.compute(molecule, basis_ri)
            mat_j_np = mat_j.to_numpy()

            if verbose:
                self.ostream.print_info('Two-center integrals for RI done in ' +
                                        f'{time.time() - ri_prep_t0:.2f} sec.')
                self.ostream.print_blank()

            ri_prep_t0 = time.time()

            if 'scipy' in sys.modules:
                lu, piv = lu_factor(mat_j_np)
                inv_mat_j_np = lu_solve((lu, piv), np.eye(mat_j_np.shape[0]))
            else:
                inv_mat_j_np = np.linalg.inv(mat_j_np)
            metric_np = np.linalg.cholesky(inv_mat_j_np)
            
            self.metric = SubMatrix(
                [0, 0, metric_np.shape[0], metric_np.shape[1]])
            self.metric.set_values(metric_np)

            if verbose:
                self.ostream.print_info('Metric for RI done in ' +
                                        f'{time.time() - ri_prep_t0:.2f} sec.')
                self.ostream.print_blank()
                self.ostream.flush()

        # broadcast Cholesky decomposed J metric
        self.metric = self.comm.bcast(self.metric)
        
    def compute_bq_vectors(self,
                           molecule,
                           basis,
                           auxiliary_basis,
                           verbose=True):
        """
        Computes cholesky decomposed J metric for the RI-JK Fock driver.
        
        :param molecule:
            The molecule to compute three-center integrals.
        :param basis:
            The basis to compute three-center integrals.
        :param auxiliary_basis:
            The auxiliary basis to compute three-center integrals.
        :param verbose:
            The information printout level.
        """
        
        if verbose:
            self.ostream.print_info(
                'Using the resolution of the identity (RI) approximation.')
            self.ostream.print_blank()
            self.ostream.flush()

        if isinstance(auxiliary_basis, str):
            if self.rank == mpi_master():
                basis_ri = MolecularBasis.read(molecule, auxiliary_basis)
            else:
                basis_ri = None
            basis_ri = self.comm.bcast(basis_ri, root=mpi_master())
        else:
            basis_ri = MolecularBasis(auxiliary_basis)

        if verbose:
            self.ostream.print_info('Dimension of RI auxiliary basis set ' +
                                    f'({basis_ri.get_label().upper()}): ' +
                                    f'{basis_ri.get_dimensions_of_basis()}')
            self.ostream.print_blank()
            self.ostream.flush()
        
        ri_prep_t0 = time.time()
        
        self._ri_drv.compute_bq_vectors(molecule, basis, basis_ri, self.metric, self.rank, self.nodes)
        
        if verbose:
            self.ostream.print_info('B^Q vectors for RI done in ' +
                                    f'{time.time() - ri_prep_t0:.2f} sec.')
            self.ostream.print_blank()
            self.ostream.flush()
        
