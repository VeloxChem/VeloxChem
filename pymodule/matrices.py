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
