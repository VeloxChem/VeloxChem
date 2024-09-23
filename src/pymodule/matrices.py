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

from .veloxchemlib import Matrices


def _Matrices_bcast(self, comm, root_rank):
    """
    Broadcasts matrices.

    :param comm:
        The MPI communicator.
    :param root_rank:
        The rank of root process.

    :return:
        Broadcasted matrices.
    """

    if comm.Get_rank() == root_rank:
        keys = self.keys()
    else:
        keys = None
    keys = comm.bcast(keys, root=root_rank)

    new_matrices = Matrices()
    for key in keys:
        if comm.Get_rank() == root_rank:
            matrix = self.matrix(key)
        else:
            matrix = None
        matrix = comm.bcast(matrix, root_rank)
        new_matrices.add(matrix, key)

    return new_matrices


def _Matrices_reduce(self, comm, root_rank):
    """
    Reduces matrices over MPI communicator to specific root process.

    :param matrices:
        The matrices to reduce.
    :param comm:
        The MPI communicator.
    :param root_rank:
        The rank of root process.

    :return:
        The reduced matrices.
    """

    reduced_matrices = None

    if comm.Get_rank() == root_rank:
        reduced_matrices = Matrices()

    for key in self.keys():
        matrix = self.matrix(key)
        reduced_matrix = matrix.reduce(comm, root_rank)
        if comm.Get_rank() == root_rank:
            reduced_matrices.add(reduced_matrix, key)

    return reduced_matrices


Matrices.bcast = _Matrices_bcast
Matrices.reduce = _Matrices_reduce
