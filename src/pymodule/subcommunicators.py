from .veloxchemlib import mpi_master
from .errorhandler import assert_msg_critical


class SubCommunicators:
    """
    Implements the MPI subcommunicator.

    :param global_comm:
        The global communicator.
    :param grps:
        The color group for creating MPI subcommunicators.

    Instance variables
        - local_comm: The local subcommunicator.
        - cross_comm: The cross subcommunicator consisting of the local master
          nodes.
    """

    def __init__(self, global_comm, grps):
        """
        Initializes the MPI subcommunicator.
        """

        global_rank = global_comm.Get_rank()
        assert_msg_critical(global_comm.Get_size() == len(grps),
                            'split_comm: inconsistent size')

        local_group = grps[global_rank]
        self.local_comm = global_comm.Split(local_group, global_rank)

        local_rank = self.local_comm.Get_rank()
        local_master = (local_rank == mpi_master())

        cross_group = 0 if local_master else 1
        self.cross_comm = global_comm.Split(cross_group, global_rank)

    def __del__(self):
        """
        Deletes the MPI subcommunicator.
        """

        self.local_comm.Disconnect()
        self.cross_comm.Disconnect()
