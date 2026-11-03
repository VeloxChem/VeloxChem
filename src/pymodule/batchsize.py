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

from os import environ
import psutil
import ctypes

from .veloxchemlib import mpi_master
from .veloxchemlib import get_number_of_threads
from .environment import get_available_memory
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


def get_ri_batch_size(input_batch_size, n_total, n_ao, density_memory, comm):
    """
    Gets the batch size of a Fock build with the resolution of the identity,
    from the memory free to each rank (the smaller of the node's and what its
    control groups leave, see environment.get_available_memory). Every rank
    must call this. Builds from the four-center integrals take
    get_batch_size.

    :param input_batch_size:
        The batch size from input, or None.
    :param n_total:
        The total number of Fock matrices.
    :param n_ao:
        The number of atomic orbitals in one Fock matrix; read on the master
        rank only.
    :param density_memory:
        The memory the solver holds on the master rank for each density of a
        batch, in bytes; read on the master rank only.
    :param comm:
        The communicator.

    :return:
        The batch size, the same on every rank: the input batch size where one
        is given, and otherwise as many densities as the smallest room of a
        rank holds, at least one and at most n_total (zero where n_total is
        zero).

    .. note::
        The RI build's own matrices are not part of density_memory: the build
        cuts its densities into chunks that fit the RI driver's budget (see
        rijkresponse._max_densities and rijresponse._max_densities).
    """

    if input_batch_size is not None:
        return comm.bcast(input_batch_size, root=mpi_master())

    # NOTE: the ranks on one node share its free memory, and are counted by
    # their host name, as get_ri_free_memory counts them for the RI budgets.
    # Each rank's room is its node's free memory over the ranks on it, and the
    # smallest room of any rank is taken.
    import platform
    hosts_and_memory = comm.allgather((platform.node(), get_available_memory()))
    hosts = [host for host, _ in hosts_and_memory]
    room = min(available // hosts.count(host)
               for host, available in hosts_and_memory)

    batch_size = None

    if comm.Get_rank() == mpi_master():
        assert_msg_critical(density_memory is not None,
                            'get_ri_batch_size: no memory of a density given')

        # NOTE: the memory of a density is no less than 8 bytes of n_ao^2
        # for each thread, as get_batch_size takes it. With many threads this
        # is above density_memory, and what lies between the two is room for
        # what neither counts, such as the matrices the RI build allocates
        # beside its densities.
        estimate = n_ao * n_ao * 8 * get_number_of_threads()
        per_density = max(estimate, density_memory)

        # NOTE: batch_size is zero where n_total is zero.
        batch_size = min(n_total, max(room // per_density, 1))

    return comm.bcast(batch_size, root=mpi_master())


def get_ri_density_memory(n_ao, n_matrices, n_factor_columns, n_sliced,
                          xc_density_memory):
    """
    Gets the memory the solver holds on the master rank for each density of a
    Fock build with the resolution of the identity, in bytes, as
    get_ri_batch_size takes it.

    :param n_ao:
        The number of atomic orbitals.
    :param n_matrices:
        The matrices of n_ao by n_ao the solver holds for each density.
    :param n_factor_columns:
        The columns of n_ao rows the solver holds for each density: the right
        factors of an RI-JK build.
    :param n_sliced:
        The densities the exchange-correlation integration slices for each
        density.
    :param xc_density_memory:
        The memory one sliced density takes, from
        dftutils.get_response_fock_xc_density_memory.

    :return:
        The memory in bytes.
    """

    return int(8 * (n_matrices * n_ao * n_ao + n_factor_columns * n_ao) +
               n_sliced * xc_density_memory)


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
