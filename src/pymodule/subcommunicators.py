#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

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
                            'SubCommunicators: inconsistent size')

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
