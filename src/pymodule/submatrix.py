from mpi4py import MPI

from .veloxchemlib import SubMatrix


def _SubMatrix_sum(sub_matrix_a, sub_matrix_b, datatype):
    """
    Helper function to enable submatrices addition.

    :param sub_matrix_a:
        The first submatrix to be added.
    :param sub_matrix_b:
        The second submatrix to be added.
    :param datatype:
        The datatype information (currently, not used).

    :return:
        The sum of two submatrices.
    """

    return sub_matrix_a + sub_matrix_b


@staticmethod
def _SubMatrix_reduce(sub_matrix, comm, mpi_id):
    """
    Reduces submatrix over MPI communicator to specific root process.

    :param sub_matrix:
        The submatrix to reduce.
    :param comm:
        The MPI communicator.
    :param mpi_id:
        The identifier of root process.

    :return:
        The reduced submatrix.
    """

    sum_op = MPI.Op.Create(_SubMatrix_sum, commute=True)
    return comm.reduce(sub_matrix, op=sum_op, root=mpi_id)


SubMatrix.reduce = _SubMatrix_reduce
