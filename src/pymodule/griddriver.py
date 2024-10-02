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

from .veloxchemlib import GridDriver
from .veloxchemlib import MolecularGrid, DenseMatrix


def _GridDriver_generate(self, molecule, comm=None):
    """
    Generates molecular grid for molecule.

    :param molecule:
        The molecule.
    :param comm:
        The MPI communicator.

    :return:
        The molecular grid.
    """

    if comm is None:
        comm = MPI.COMM_WORLD

    local_grid = self.generate_local_grid(molecule, comm.Get_rank(),
                                          comm.Get_size())

    grid_np_arrays = comm.allgather(local_grid.grid_to_numpy())
    grid_np_arrays = [arr for arr in grid_np_arrays if arr.size > 0]

    # Note: use hstack since local_grid is of shape (4,N)
    mol_grid = MolecularGrid(DenseMatrix(np.hstack(grid_np_arrays)))

    mol_grid.partition_grid_points()
    mol_grid.distribute_counts_and_displacements(comm.Get_rank(),
                                                 comm.Get_size())

    return mol_grid


GridDriver.generate = _GridDriver_generate
