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

from .veloxchemlib import PackedMatrix
from .veloxchemlib import mpi_master
from .errorhandler import assert_msg_critical


def _PackedMatrix_broadcast(self, comm, root=mpi_master()):
    """
    Broadcasts a packed matrix.

    :param comm:
        The MPI communicator.
    :param root:
        The rank to broadcast from.

    :return:
        The packed matrix, which is this one on the root and a new one
        elsewhere.
    """

    # NOTE: the values are sent as a buffer rather than pickled. The metric of a
    # fitting basis of sixteen thousand functions is a gigabyte of them, which is
    # slow to pickle and near the limit of what mpi4py will pickle at all. The
    # shape has to travel first, as the ranks which receive have nothing to make
    # a matrix of the right size from.

    if comm.Get_size() == 1:
        return self

    if comm.Get_rank() == root:
        shape = (self.number_of_rows(), self.number_of_columns(),
                 self.get_type())
    else:
        shape = None

    shape = comm.bcast(shape, root=root)

    matrix = self if comm.Get_rank() == root else PackedMatrix(*shape)

    comm.Bcast(matrix.values_view(), root=root)

    return matrix


def _PackedMatrix_reduce(self, comm, root=mpi_master(), op=MPI.SUM):
    """
    Reduces packed matrices of the same shape onto the root rank.

    :param comm:
        The MPI communicator.
    :param root:
        The rank to reduce onto.
    :param op:
        The reduction operation.

    :return:
        The packed matrix on the root rank and None elsewhere.
    """

    if comm.Get_size() == 1:
        return self if comm.Get_rank() == root else None

    _PackedMatrix_check_shapes(self, comm)

    if comm.Get_rank() == root:
        total = PackedMatrix(self.number_of_rows(), self.number_of_columns(),
                             self.get_type())

        comm.Reduce(self.values_view(), total.values_view(), op=op, root=root)

        return total

    comm.Reduce(self.values_view(), None, op=op, root=root)

    return None


def _PackedMatrix_allreduce(self, comm, op=MPI.SUM):
    """
    Reduces packed matrices of the same shape onto every rank.

    :param comm:
        The MPI communicator.
    :param op:
        The reduction operation.

    :return:
        The packed matrix.
    """

    if comm.Get_size() == 1:
        return self

    _PackedMatrix_check_shapes(self, comm)

    total = PackedMatrix(self.number_of_rows(), self.number_of_columns(),
                         self.get_type())

    comm.Allreduce(self.values_view(), total.values_view(), op=op)

    return total


def _PackedMatrix_check_shapes(self, comm):
    """
    Checks that every rank holds a matrix of the same shape.

    :param comm:
        The MPI communicator.
    """

    # NOTE: a reduction of buffers of different lengths is not an error a
    # communicator reports usefully, so it is refused here instead.

    shape = (self.number_of_rows(), self.number_of_columns(), self.get_type())

    shapes = comm.allgather(shape)

    assert_msg_critical(
        all(other == shapes[0] for other in shapes),
        'PackedMatrix.reduce: every rank must hold a matrix of the same shape')


PackedMatrix.broadcast = _PackedMatrix_broadcast
PackedMatrix.reduce = _PackedMatrix_reduce
PackedMatrix.allreduce = _PackedMatrix_allreduce
