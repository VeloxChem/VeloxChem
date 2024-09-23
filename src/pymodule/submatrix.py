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

from .veloxchemlib import SubMatrix


def _SubMatrix_sum(sub_matrix_a, sub_matrix_b, datatype):
    """
    Helper function to enable submatrices addition.

    :param sub_matrix_a:
        The first submatrix to be added.
    :param sub_matrix_b:
        The second submatrix to be added.
    :param datatype:
        The datatype information (currently, not used).

    :return:
        The sum of two submatrices.
    """

    return sub_matrix_a + sub_matrix_b


@staticmethod
def _SubMatrix_reduce(sub_matrix, comm, mpi_id):
    """
    Reduces submatrix over MPI communicator to specific root process.

    :param sub_matrix:
        The submatrix to reduce.
    :param comm:
        The MPI communicator.
    :param mpi_id:
        The identifier of root process.

    :return:
        The reduced submatrix.
    """

    sum_op = MPI.Op.Create(_SubMatrix_sum, commute=True)
    return comm.reduce(sub_matrix, op=sum_op, root=mpi_id)


SubMatrix.reduce = _SubMatrix_reduce
