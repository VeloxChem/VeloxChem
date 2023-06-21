import numpy as np

from mpi4py import MPI

from .veloxchemlib import SubMatrix
from .mpitools import is_master


@staticmethod
def _SubMatrix_bcast(mat, rank, comm):
    """
    Broadcast submatrix over givem MPI communicator.

    :param mat:
        The submatrix.
    :param rank:
        The rank of MPI process.
    :param comm:
        MPI communicator.
    """

    if comm.Get_size() == 1:
        return mat

    # communicate submatrix dimensions
    if is_master(rank):
        dims = mat.get_dimensions()
    else:
        dims = None
    dims = comm.bcast(dims)

    # broadcast submatrix data
    if is_master(rank):
        mdata = mat.to_numpy()
    else:
        mat = SubMatrix(dims)
        mdata = np.zeros([mat.number_of_rows(), mat.number_of_columns()])
    comm.Barrier()
    comm.Bcast(mdata)

    # update submatrix on workers
    if not is_master(rank):
        mat.set_values(mdata)

    return mat


@staticmethod
def _SubMatrix_reduce(mat, rank, comm):
    """
    Reduces submatrix over givem MPI communicator.

    :param mat:
        The submatrix.
    :param rank:
        The rank of MPI process.
    :param comm:
        MPI communicator.
    """

    if comm.Get_size() == 1:
        return mat

    if is_master(rank):
        mdata = mat.to_numpy()
        rdata = np.zeros(mdata.shape)
    else:
        mdata = mat.to_numpy()
        rdata = None

    comm.Barrier()
    comm.Reduce(mdata, rdata, op=MPI.SUM)

    if is_master(rank):
        mat.set_values(rdata)

    return mat


def _SubMatrix_deepcopy(self, memo):
    """
    Implements deepcopy.

    :param memo:
        The memo dictionary for deepcopy.

    :return:
        A deepcopy of self.
    """

    return SubMatrix(self)


SubMatrix.bcast = _SubMatrix_bcast
SubMatrix.reduce = _SubMatrix_reduce
SubMatrix.__deepcopy__ = _SubMatrix_deepcopy
