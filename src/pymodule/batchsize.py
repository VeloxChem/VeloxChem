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
from os import environ
import psutil
import ctypes

from .veloxchemlib import mpi_master
from .environment import get_ri_free_memory
from .errorhandler import assert_msg_critical


def get_batch_size(input_batch_size, n_total, n_ao, comm):
    """
    Gets batch size for Fock matrix computation.

    :param input_batch_size:
        The batch size from input.
    :param n_total:
        The total number of Fock matrices.
    :param n_ao:
        The number of atomic orbitals in one Fock matrix.
    :param comm:
        The communicator.

    :return:
        The batch size for Fock matrix computation..
    """

    batch_size = input_batch_size

    # check available memories on the nodes
    total_mem = psutil.virtual_memory().total
    total_mem_list = comm.gather(total_mem, root=mpi_master())

    # NOTE: the ranks on one node share its free memory. SLURM_NTASKS_PER_NODE
    # says how many there are; where it is not set, each rank counts the ranks
    # with its own host name, and the master's count is used below. The
    # allgather is a collective, and it relies on every rank of a job seeing
    # the same environment: either every rank finds SLURM_NTASKS_PER_NODE set
    # and none gathers, or every rank gathers.
    if 'SLURM_NTASKS_PER_NODE' in environ:
        ntasks_per_node = int(environ['SLURM_NTASKS_PER_NODE'])
    else:
        import platform
        host = platform.node()
        ntasks_per_node = comm.allgather(host).count(host)

    if comm.Get_rank() == mpi_master():

        # check if master node has larger memory
        mem_adjust = 0.0
        if total_mem > min(total_mem_list):
            mem_adjust = total_mem - min(total_mem_list)

        # computes maximum batch size from available memory
        avail_mem = psutil.virtual_memory().available - mem_adjust
        avail_mem //= ntasks_per_node
        mem_per_mat = n_ao**2 * ctypes.sizeof(ctypes.c_double)
        nthreads = int(environ['OMP_NUM_THREADS'])
        # TODO: double check the estimation of max_batch_size
        max_batch_size = int(avail_mem / (mem_per_mat * nthreads))
        max_batch_size = max(1, max_batch_size)

        # note: batch_size will be zero if n_total is zero
        if batch_size is None:
            batch_size = min(n_total, max_batch_size)

    batch_size = comm.bcast(batch_size, root=mpi_master())

    return batch_size


def get_ri_batch_size(input_batch_size, n_total, item_cost, fixed_cost, comm):
    """
    Gets the batch size of a Fock build with the resolution of the identity,
    from the memory free to each rank when the build starts. Every rank must
    call this. Builds without the resolution of the identity take
    get_batch_size.

    :param input_batch_size:
        The batch size from input, or None.
    :param n_total:
        The number of items to build: densities, or groups of densities which
        the caller does not split.
    :param item_cost:
        The memory each item takes on this rank during the build, in bytes.
    :param fixed_cost:
        The memory the build takes on this rank whatever the number of items,
        beyond what the rank already holds, in bytes.
    :param comm:
        The communicator.

    :return:
        The batch size, the same on every rank: the input batch size where one
        is given, and otherwise the most items every rank has room for, at
        least one and at most n_total (zero where n_total is zero).
    """

    if input_batch_size is not None:
        return comm.bcast(input_batch_size, root=mpi_master())

    assert_msg_critical(item_cost > 0,
                        'get_ri_batch_size: an item has to cost memory')

    # NOTE: read when the build starts, so what the RI driver and the solver
    # already hold is not counted as free. get_ri_free_memory keeps back the
    # node's reserve and divides the rest among the ranks on the node.
    room = get_ri_free_memory(comm, 0) - fixed_cost

    # NOTE: the costs differ between the ranks (the master holds more for each
    # density than the others), and every rank builds the same batch, so each
    # rank works out how many items it has room for and the smallest count is
    # taken.
    count = comm.allreduce(max(room, 0) // item_cost, op=MPI.MIN)

    batch_size = None

    if comm.Get_rank() == mpi_master():
        batch_size = min(max(count, 1), n_total)

    return comm.bcast(batch_size, root=mpi_master())


def get_number_of_batches(n_total, batch_size, comm):
    """
    Gets number of batches for Fock matrix computation.

    :param n_total:
        The total number of Fock matrices.
    :param batch_size:
        The batch size for Fock matrix computation.
    :param comm:
        The communicator.

    :return:
        The number of batches for Fock matrix computation.
    """

    num_batches = None

    # note: num_batches will be zero if batch_size is zero
    if comm.Get_rank() == mpi_master():
        if batch_size > 0:
            num_batches = n_total // batch_size
            if n_total % batch_size != 0:
                num_batches += 1
        else:
            num_batches = 0

    num_batches = comm.bcast(num_batches, root=mpi_master())

    return num_batches
