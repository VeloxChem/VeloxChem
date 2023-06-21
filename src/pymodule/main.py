from mpi4py import MPI

from .errorhandler import assert_msg_critical


def main():
    """
    Runs VeloxChem with command line arguments.
    """

    assert_msg_critical(False, "I am dummy!")
