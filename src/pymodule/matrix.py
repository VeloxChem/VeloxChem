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

from .veloxchemlib import Matrix


def _Matrix_sum(matrix_a, matrix_b, datatype):
    """
    Helper function to enable matrices addition.

    :param matrix_a:
        The first matrix to be added.
    :param sub_matrix_b:
        The second matrix to be added.
    :param datatype:
        The datatype information (currently, not used).

    :return:
        The sum of two matrices.
    """
    return matrix_a + matrix_b


def _Matrix_reduce(self, comm, root_rank):
    """
    Reduces matrix over MPI communicator to specific root process.

    TODO: Replace this implementation with direct submatrices reduction
    code to overcome 2Gb limit in MPI4PY communications.

    :param comm:
        The MPI communicator.
    :param root_rank:
        The rank of root process.

    :return:
        The reduced matrix.
    """

    sum_op = MPI.Op.Create(_Matrix_sum, commute=True)
    mat = comm.reduce(self, op=sum_op, root=root_rank)
    sum_op.Free()
    return mat


Matrix.reduce = _Matrix_reduce
