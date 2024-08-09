from mpi4py import MPI

from .veloxchemlib import Matrices
from .matrix import Matrix


@staticmethod
def _Matrices_bcast(matrices, comm, mpi_id):
    """
    Broadcasts matrices.

    :param matrices:
        The matrices to be broadcasted.
    :param comm:
        The MPI communicator.
    :param mpi_id:
        The identifier of root process.

    :return:
        Broadcasted matrices.
    """

    keys = None
    if comm.Get_rank() == mpi_id:
        keys = matrices.keys()
    keys = comm.bcast(keys)
    comm.Barrier()

    new_matrices = Matrices()
    for key in keys:
        matrix = None
        if comm.Get_rank() == mpi_id:
            matrix = matrices.matrix(key)
        matrix = comm.bcast(matrix, mpi_id)
        comm.Barrier()
        new_matrices.add(matrix, key)
        comm.Barrier()

    return new_matrices


@staticmethod
def _Matrices_reduce(matrices, comm, mpi_id):
    """
    Reduces matrices over MPI communicator to specific root process.
    
    :param matrices:
        The matrices to reduce.
    :param comm:
        The MPI communicator.
    :param mpi_id:
        The identifier of root process.

    :return:
        The reduced matrices.
    """

    red_matrices = None
    if comm.Get_rank() == mpi_id:
        red_matrices = Matrices()

    for key in matrices.keys():
        matrix = matrices.matrix(key)
        red_matrix = Matrix.reduce(matrix, comm, mpi_id)
        comm.Barrier()
        if red_matrices is not None:
            red_matrices.add(red_matrix, key)
        comm.Barrier()

    return red_matrices


Matrices.bcast = _Matrices_bcast
Matrices.reduce = _Matrices_reduce
