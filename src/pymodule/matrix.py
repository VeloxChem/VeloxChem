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
