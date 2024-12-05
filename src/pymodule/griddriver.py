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

    def generate(self, molecule, max_points_per_box=None, verbose=False):
        """
        Generates molecular grid for molecule.

        :param molecule:
            The molecule.

        :return:
            The molecular grid.
        """

        local_grid = self._grid_drv._generate_local_grid(
            molecule, self.rank, self.nodes)

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
