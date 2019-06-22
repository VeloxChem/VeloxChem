from mpi4py import MPI

from .veloxchemlib import assert_msg_critical as vlx_assert


def assert_msg_critical(condition, msg=''):
    """Asserts that the condition is true.

    Asserts that the condition is true. Otherwise terminates the program with
    an error message.

    Parameters
    ----------
    condition
        The condition.
    msg
        The error message.
    """
    if __debug__ and MPI.COMM_WORLD.Get_size() == 1:
        assert condition, msg
    else:
        vlx_assert(condition, msg)
