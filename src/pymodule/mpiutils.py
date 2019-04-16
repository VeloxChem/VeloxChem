from .veloxchemlib import mpi_master
from .errorhandler import assert_msg_critical


def split_comm(comm, grps):

    global_comm = comm
    global_rank = global_comm.Get_rank()

    assert_msg_critical(global_comm.Get_size() == len(grps),
                        'split_comm: inconsistent size')

    local_group = grps[global_rank]
    local_comm = global_comm.Split(local_group, global_rank)
    local_rank = local_comm.Get_rank()
    local_master = (local_rank == mpi_master())

    cross_group = 0 if local_master else 1
    cross_comm = global_comm.Split(cross_group, global_rank)

    return local_comm, cross_comm
