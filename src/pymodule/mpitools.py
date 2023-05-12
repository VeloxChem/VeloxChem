from mpi4py import MPI

def is_master(comm=MPI.COMM_WORLD):
    """
    Checks if this rank of MPI communicator is master rank.

    :param comm:
        The MPI communicator.
    :return:
        True if this rank of MPI communicator is master rank, False otherwise.
    """

    if comm.Get_rank() == 0:
        return True
    else:
        return False


