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
