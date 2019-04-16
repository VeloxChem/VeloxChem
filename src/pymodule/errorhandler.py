from mpi4py import MPI

from .veloxchemlib import assert_msg_critical as vlx_assert


def assert_msg_critical(condition, msg=''):
    if __debug__ and MPI.COMM_WORLD.Get_size() == 1:
        assert condition, msg
    else:
        vlx_assert(condition, msg)
