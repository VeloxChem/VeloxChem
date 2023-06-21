from .veloxchemlib import Matrix
from .mpitools import is_master
from .submatrix import SubMatrix


@staticmethod
def _Matrix_bcast(mat, rank, comm):
    """
    Broadcast matrix over givem MPI communicator.

    :param mat:
        The matrix.
    :param rank:
        The rank of MPI process.
    :param comm:
        MPI communicator.
    """

    if comm.Get_size() == 1:
        return mat

    # communicate angular pairs
    if is_master(rank):
        angpairs = mat.get_angular_pairs()
    else:
        angpairs = None
    angpairs = comm.bcast(angpairs)

    # allocate matrix
    if not is_master(rank):
        mat = Matrix()

    # broadcast submatrices
    for angpair in angpairs:
        if is_master(rank):
            submat = mat.get_submatrix(angpair)
        else:
            submat = None
        submat = SubMatrix.bcast(submat, rank, comm)
        if not is_master(rank):
            mat.add(submat, angpair)
        comm.Barrier()

    # broadcast matrix type
    if is_master(rank):
        mtype = mat.get_type()
    else:
        mtype = None
    mtype = comm.bcast(mtype)
    if not is_master(rank):
        mat.set_type(mtype)

    return mat


@staticmethod
def _Matrix_reduce(mat, rank, comm):
    """
    Reduces matrix over givem MPI communicator.

    :param mat:
        The matrix.
    :param rank:
        The rank of MPI process.
    :param comm:
        MPI communicator.
    """

    if comm.Get_size() == 1:
        return mat

    # set up reduced matrix type
    if is_master(rank):
        rmat = Matrix()
        rmat.set_type(mat.get_type())

    # reduce submatrices
    for angpair in mat.get_angular_pairs():
        submat = mat.get_submatrix(angpair)
        submat = SubMatrix.reduce(submat, rank, comm)
        if is_master(rank):
            rmat.add(submat, angpair)
        print('rank : ', rank, ' ', angpair, ' ', submat.to_numpy())
        comm.Barrier()

    if is_master(rank):
        return rmat
    else:
        return mat


def _Matrix_deepcopy(self, memo):
    """
    Implements deepcopy.

    :param memo:
        The memo dictionary for deepcopy.

    :return:
        A deepcopy of self.
    """

    return Matrix(self)


Matrix.bcast = _Matrix_bcast
Matrix.reduce = _Matrix_reduce
Matrix.__deepcopy__ = _Matrix_deepcopy
