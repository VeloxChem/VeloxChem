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

from .veloxchemlib import AODensityMatrix
from .veloxchemlib import denmat
from .veloxchemlib import mpi_master


def _AODensityMatrix_broadcast(self, comm, root=mpi_master()):
    """
    Broadcasts AODensityMatrix object.

    :param comm:
        The MPI communicator.
    :param root:
        The root rank to broadcast from.

    :return:
        The AODensityMatrix object.
    """

    densities = []

    if comm.Get_rank() == root:
        density_type = ('rest'
                        if self.get_density_type() == denmat.rest else 'unrest')
        num_dens = self.number_of_density_matrices()
    else:
        density_type = None
        num_dens = None

    density_type = comm.bcast(density_type, root=root)
    num_dens = comm.bcast(num_dens, root=root)

    if comm.Get_rank() == root:
        for i in range(num_dens):
            densities.append(self.alpha_to_numpy(i))
            if density_type == 'unrest':
                densities.append(self.beta_to_numpy(i))

    else:
        for i in range(num_dens):
            densities.append(None)
            if density_type == 'unrest':
                densities.append(None)

    if density_type == 'rest':
        for i in range(num_dens):
            densities[i] = comm.bcast(densities[i], root=root)
    elif density_type == 'unrest':
        for i in range(num_dens * 2):
            densities[i] = comm.bcast(densities[i], root=root)

    density_type = (denmat.rest if density_type == 'rest' else denmat.unrest)

    return AODensityMatrix(densities, density_type)


AODensityMatrix.broadcast = _AODensityMatrix_broadcast
