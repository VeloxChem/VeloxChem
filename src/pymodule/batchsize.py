import os
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

    batch_size = input_batch_size

    # check available memories on the nodes
    total_mem = psutil.virtual_memory().total
    total_mem_list = comm.gather(total_mem, root=mpi_master())

    if comm.Get_rank() == mpi_master():

        # check if master node has larger memory
        mem_adjust = 0.0
        if total_mem > min(total_mem_list):
            mem_adjust = total_mem - min(total_mem_list)

        # computes maximum batch size from available memory
        avail_mem = psutil.virtual_memory().available - mem_adjust
        mem_per_mat = n_ao**2 * ctypes.sizeof(ctypes.c_double)
        nthreads = int(os.environ['OMP_NUM_THREADS'])
        max_batch_size = int(avail_mem / mem_per_mat / (0.625 * nthreads))
        max_batch_size = max(1, max_batch_size)

        if batch_size is None:
            batch_size = min(n_total, max_batch_size)

    batch_size = comm.bcast(batch_size, root=mpi_master())

    return batch_size


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

    if comm.Get_rank() == mpi_master():
        # get number of batches
        num_batches = n_total // batch_size
        if n_total % batch_size != 0:
            num_batches += 1

    num_batches = comm.bcast(num_batches, root=mpi_master())

    return num_batches
