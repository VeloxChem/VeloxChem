from mpi4py import MPI
from sys import stderr

def assert_msg_critical(condition, msg=''):
    """
    Asserts that the condition is true. Otherwise terminates the program with
    an error message.

    :param condition:
        The condition.
    :param msg:
        The error message.
    """
    if __debug__ and MPI.COMM_WORLD.Get_size() == 1:
        assert condition, msg
    else:
        if not condition:
            stderr.write(' **** Critical Error (process {}) **** {}\n'
                         .format(MPI.COMM_WORLD.Get_rank(), msg))
            stderr.flush()
            MPI.COMM_WORLD.Abort()
