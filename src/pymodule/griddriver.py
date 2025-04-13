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

from mpi4py import MPI
import numpy as np
import sys

from .veloxchemlib import mpi_master
from .veloxchemlib import _GridDriver
from .veloxchemlib import MolecularGrid, DenseMatrix
from .outputstream import OutputStream


class GridDriver:
    """
    Implements grid driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes grid driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream

        self._grid_drv = _GridDriver()

    def set_level(self, grid_level):
        """
        Sets grid level for grid driver.

        :param grid_level:
            The grid level.
        """

        grid_level = self.comm.bcast(grid_level, root=mpi_master())

        self._grid_drv.set_level(grid_level)

    def generate(self, molecule, num_gpus_per_node, max_points_per_box=None, verbose=False):
        """
        Generates molecular grid for molecule.

        :param molecule:
            The molecule.

        :return:
            The molecular grid.
        """

        local_grid = self._grid_drv._generate_local_grid(
            molecule, self.rank, self.nodes, num_gpus_per_node)

        grid_np_arrays = self.comm.allgather(local_grid.grid_to_numpy())
        grid_np_arrays = [arr for arr in grid_np_arrays if arr.size > 0]

        # Note: use hstack since local_grid is of shape (4,N)
        if max_points_per_box is None:
            mol_grid = MolecularGrid(DenseMatrix(np.hstack(grid_np_arrays)))
        else:
            mol_grid = MolecularGrid(DenseMatrix(np.hstack(grid_np_arrays)),
                                     max_points_per_box)

        grid_summary = mol_grid.partition_grid_points()

        if (self.rank == mpi_master()) and verbose:
            for line in grid_summary.splitlines():
                self.ostream.print_info(line)
            self.ostream.print_blank()
            self.ostream.flush()

        mol_grid.distribute_counts_and_displacements(self.rank, self.nodes)

        return mol_grid
