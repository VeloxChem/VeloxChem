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

from .veloxchemlib import RIFockGradDriver
from .veloxchemlib import RIFockDriver
from .veloxchemlib import TwoCenterElectronRepulsionDriver
from .veloxchemlib import SubMatrix
from .veloxchemlib import partition_atoms
from .veloxchemlib import mpi_master


def _RIFockGradDriver_mpi_compute(self, comm, molecule, basis, aux_basis,
                                  density):
    """
    Computes RI-J part of molecular gradient.

    :param comm:
        The MPI communicator.
    :param molecule:
        The molecule to compute gradient.
    :param basis:
        The basis set to compute gradient.
    :param aux_basis:
        The auxilary basis set to compute gradient.
    :param density:
        The density matrix to compute gradient.

    :return:
        The RI-J part of molecular gradient.
    """

    # set up rank, node data
    rank = comm.Get_rank()
    nodes = comm.Get_size()

    natoms = molecule.number_of_atoms()
    grad_mat = np.zeros((natoms, 3))
    tot_grad = np.zeros((natoms, 3))

    # compute inversion of B_q metric on master node

    gvec = None
    if comm.Get_rank() == mpi_master():
        # inverted J metrix
        t2c_drv = TwoCenterElectronRepulsionDriver()
        matj = t2c_drv.compute(molecule, aux_basis)
        rmatj = np.linalg.inv(matj.full_matrix().to_numpy())
        invmatj = SubMatrix([0, 0, rmatj.shape[0], rmatj.shape[0]])
        invmatj.set_values(rmatj)
        # compute B_q vector
        ri_fock_drv = RIFockDriver(invmatj)
        ri_fock_drv.prepare_buffers(molecule, basis, aux_basis)
        gvec = ri_fock_drv.compute_bq_vector(density)
    gvec = comm.bcast(gvec, root=mpi_master())

    # compute local gradient on MPI rank
    atoms = partition_atoms(natoms, rank, nodes)
    ri_grad_drv = RIFockGradDriver()
    loc_grad = ri_grad_drv.compute(basis, aux_basis, molecule, gvec, density,
                                   atoms)
    for i in range(len(atoms)):
        gxyz = loc_grad[i].coordinates()
        grad_mat[atoms[i], 0] = gxyz[0]
        grad_mat[atoms[i], 1] = gxyz[1]
        grad_mat[atoms[i], 2] = gxyz[2]
    comm.Reduce([grad_mat, MPI.DOUBLE], [tot_grad, MPI.DOUBLE],
                op=MPI.SUM,
                root=mpi_master())

    return tot_grad


RIFockGradDriver.mpi_compute = _RIFockGradDriver_mpi_compute
