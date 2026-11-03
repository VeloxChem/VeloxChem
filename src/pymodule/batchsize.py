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

    return _batch_size_from_memory(input_batch_size, n_total, n_ao, None, comm)


def get_ri_batch_size(input_batch_size, n_total, n_ao, density_memory, comm):
    """
    Gets the batch size of a Fock build with the resolution of the identity:
    as get_batch_size, but with the memory of a density no less than what the
    caller says the solver holds for it. Builds without the resolution of the
    identity take get_batch_size.

    :param input_batch_size:
        The batch size from input, or None.
    :param n_total:
        The total number of Fock matrices.
    :param n_ao:
        The number of atomic orbitals in one Fock matrix.
    :param density_memory:
        The memory the solver holds on the master rank for each density of a
        batch, in bytes. Only the master's value is read.
    :param comm:
        The communicator.

    :return:
        The batch size, which is never more than get_batch_size gives.

    .. note::
        The RI build's own matrices are not part of density_memory: the build
        cuts its densities into chunks that fit the RI driver's budget (see
        rijkresponse._max_densities and rijresponse._max_densities).
    """

    return _batch_size_from_memory(input_batch_size, n_total, n_ao,
                                   density_memory, comm)


def _batch_size_from_memory(input_batch_size, n_total, n_ao, density_memory,
                            comm):
    """
    Gets the batch size for get_batch_size and get_ri_batch_size.

    :param density_memory:
        The memory the solver holds on the master rank for each density, in
        bytes, or None for the estimate of get_batch_size alone.

    The other parameters are those of get_ri_batch_size.
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
        mem_per_density = mem_per_mat * nthreads
        # NOTE: the larger of this estimate and the caller's density_memory,
        # so that density_memory can make a batch smaller than get_batch_size
        # gives and never larger.
        if density_memory is not None:
            mem_per_density = max(mem_per_density, density_memory)
        # TODO: double check the estimation of max_batch_size
        max_batch_size = int(avail_mem / mem_per_density)
        max_batch_size = max(1, max_batch_size)

        # note: batch_size will be zero if n_total is zero
        if batch_size is None:
            batch_size = min(n_total, max_batch_size)

    batch_size = comm.bcast(batch_size, root=mpi_master())

    return batch_size


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
