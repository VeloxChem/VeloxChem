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

from .veloxchemlib import RIFockDriver
from .veloxchemlib import partition_atoms
from .veloxchemlib import mpi_master
from .veloxchemlib import Matrix

def _RIFockDriver_mpi_prepare_buffers(self, comm, molecule, basis, aux_basis):
    """
    Computes distributed three-center electron integrals tensor.

    :param comm:
        The MPI communicator.
    :param molecule:
        The molecule to compute gradient.
    :param basis:
        The basis set to compute gradient.
    :param aux_basis:
        The auxilary basis set to compute gradient.
    """

    # set up rank, node data
    
    rank = comm.Get_rank()
    nodes = comm.Get_size()

    # select three-center integrals computation method
    
    natoms = molecule.number_of_atoms()
    atoms = partition_atoms(natoms, rank, nodes)
    self.prepare_buffers(molecule, basis, aux_basis, atoms)
        
def _RIFockDriver_mpi_compute_bq_vector(self, comm, density):
    """
    Computes distributed Bq vector.

    :param comm:
        The MPI communicator.
    :param density:
        The density to compute Bq vector.
        
    :return:
        The Bq vector.
    """

    # set up rank, node data
    gv = np.array(self.compute_local_bq_vector(density))
    tv = np.zeros(gv.shape)
    comm.Allreduce([gv, MPI.DOUBLE], [tv, MPI.DOUBLE], op=MPI.SUM)
    return list(tv)

def _RIFockDriver_mpi_compute(self, comm, gamma, density, label):
    """
    Computes Fock matrix.

    :param comm:
        The MPI communicator.
    :param gamma:
        The Bq vector.    
    :param density:
        The density to compute Fock matrix.
    :param label:
        The type of Fock matrix.
        
    :return:
        The Fock matrix.
    """

    # set up rank, node data
    fmat = self.local_compute(density, gamma, label)
    rfmat = Matrix.reduce(fmat, comm, mpi_master())
    return rfmat

RIFockDriver.mpi_prepare_buffers = _RIFockDriver_mpi_prepare_buffers
RIFockDriver.mpi_compute_bq_vector = _RIFockDriver_mpi_compute_bq_vector
RIFockDriver.mpi_compute = _RIFockDriver_mpi_compute
