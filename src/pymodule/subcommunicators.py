#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
