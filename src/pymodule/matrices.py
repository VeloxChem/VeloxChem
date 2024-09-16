from .veloxchemlib import Matrices


@staticmethod
def _Matrices_bcast(matrices, comm, root_rank):
    """
    Broadcasts matrices.

    :param matrices:
        The matrices to be broadcasted.
    :param comm:
        The MPI communicator.
    :param root_rank:
        The rank of root process.

    :return:
        Broadcasted matrices.
    """

    if comm.Get_rank() == root_rank:
        keys = matrices.keys()
    else:
        keys = None
    keys = comm.bcast(keys, root=root_rank)

    new_matrices = Matrices()
    for key in keys:
        if comm.Get_rank() == root_rank:
            matrix = matrices.matrix(key)
        else:
            matrix = None
        matrix = comm.bcast(matrix, root_rank)
        new_matrices.add(matrix, key)

    return new_matrices


def _Matrices_reduce(self, comm, root_rank):
    """
    Reduces matrices over MPI communicator to specific root process.

    :param matrices:
        The matrices to reduce.
    :param comm:
        The MPI communicator.
    :param root_rank:
        The rank of root process.

    :return:
        The reduced matrices.
    """

    reduced_matrices = None

    if comm.Get_rank() == root_rank:
        reduced_matrices = Matrices()

    for key in self.keys():
        matrix = self.matrix(key)
        reduced_matrix = matrix.reduce(comm, root_rank)
        if comm.Get_rank() == root_rank:
            reduced_matrices.add(reduced_matrix, key)

    return reduced_matrices


Matrices.bcast = _Matrices_bcast
Matrices.reduce = _Matrices_reduce
