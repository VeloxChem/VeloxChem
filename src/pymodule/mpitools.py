def is_master(rank):
    """
    Checks if this rank of MPI communicator is master rank.

    :param rank:
        The rank of MPI communicator.
    :return:
        True if this rank of MPI communicator is master rank, False otherwise.
    """

    if rank == 0:
        return True
    else:
        return False
