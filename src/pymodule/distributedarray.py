#
#                           VELOXCHEM 1.0-RC
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
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
import time as tm
import h5py

from .veloxchemlib import mpi_master
from .veloxchemlib import mpi_size_limit


class DistributedArray:
    """
    Impements distributed array.

    :param array:
        The numpy array stored on the master node.
    :param comm:
        The communicator.
    :param distribute:
        The flag for distributing the array via scatter.

    Instance variable
        - comm: The communicator.
        - rank: The MPI rank.
        - nodes: Number of MPI processes.
        - data: The numpy array stored in the distributed array.
    """

    def __init__(self, array, comm, distribute=True):
        """
        Initializes distributed array.
        """

        self.comm = comm
        self.rank = comm.Get_rank()
        self.nodes = comm.Get_size()

        self.data = None

        if not distribute:
            self.data = array.copy()
            return

        if self.rank == mpi_master():
            # determine batch size for scatter
            batch_size = mpi_size_limit() // array.itemsize
            if array.ndim == 2 and array.shape[1] != 0:
                batch_size //= array.shape[1]
            batch_size //= self.nodes
            batch_size = max(1, batch_size)

            # determine counts and displacements for scatter
            ave, res = divmod(array.shape[0], self.nodes)
            counts = [ave + 1 if p < res else ave for p in range(self.nodes)]
            displacements = [sum(counts[:p]) for p in range(self.nodes)]

            # determine number of batches for scatter
            num_batches = max(counts) // batch_size
            if max(counts) % batch_size != 0:
                num_batches += 1
        else:
            num_batches = None

        num_batches = comm.bcast(num_batches, root=mpi_master())

        for batch_ind in range(num_batches):
            if self.rank == mpi_master():
                array_list = []
                for p in range(self.nodes):
                    row_start = batch_size * batch_ind + displacements[p]
                    row_end = min(row_start + batch_size,
                                  counts[p] + displacements[p])
                    array_list.append(array[row_start:row_end])
            else:
                array_list = None

            recv_data = comm.scatter(array_list, root=mpi_master())

            if self.data is None:
                self.data = recv_data.copy()
            elif self.data.ndim == 1:
                self.data = np.hstack((self.data, recv_data))
            elif self.data.ndim == 2:
                self.data = np.vstack((self.data, recv_data))

    def __sizeof__(self):
        """
        Estimates the size of the distributed array, regardless whether it owns
        the data or not.

        :return:
            The esitmated size of the distributed array.
        """

        return self.data.size * self.data.itemsize

    def array(self):
        """
        Returns the numpy array stored in the object.

        :return:
            The numpy array.
        """

        return self.data

    def shape(self, axis):
        """
        Returns the shape of the distributed array along axis.

        :param axis:
            The axis.
        :return:
            The shape along axis.
        """

        return self.data.shape[axis]

    def squared_norm(self, axis=None):
        """
        Returns the squared norm of the distributed array along axis.

        :param axis:
            The axis.
        :return:
            The squared norm along axis.
        """

        n2 = np.linalg.norm(self.data, axis=axis)**2
        n2 = self.comm.allreduce(n2, op=MPI.SUM)
        return n2

    def norm(self, axis=None):
        """
        Returns the norm of the distributed array along axis.

        :param axis:
            The axis.
        :return:
            The norm along axis.
        """

        n2 = np.linalg.norm(self.data, axis=axis)**2
        n2 = self.comm.allreduce(n2, op=MPI.SUM)
        return np.sqrt(n2)

    def dot(self, i, dist_array, j):
        """
        Returns the dot product between a column vector in self and a column
        vector in a distributed array.

        :param i:
            The column index in self.
        :param dist_array:
            The distributed array.
        :param j:
            The column index in dist_array.
        :return:
            The dot product.
        """

        if self.data.ndim == 2 and dist_array.data.ndim == 2:
            dot_prod = np.dot(self.data[:, i], dist_array.data[:, j])
            dot_prod = self.comm.allreduce(dot_prod, op=MPI.SUM)
            return dot_prod
        else:
            return None

    def get_full_vector(self, col=None):
        """
        Gets a full column vector from a distributed array.

        :param: col:
            The column index (used only when self.data.ndim is 2).
        :return:
            The full vector on the master node, None on other nodes.
        """

        data = None

        if self.data.ndim == 1:
            data = self.comm.gather(self.data, root=mpi_master())
        elif self.data.ndim == 2 and col is not None:
            data = self.comm.gather(self.data[:, col], root=mpi_master())

        if self.rank == mpi_master():
            full_shape_0 = sum([m.shape[0] for m in data])
            return np.hstack(data).reshape(full_shape_0)
        else:
            return None

    def matmul_AtB(self, dist_array, factor=None):
        """
        Computes matrix-matrix multiplication between self.T and a distributed
        array.

        :param dist_array:
            The distributed array.
        :param factor:
            The factor to be multiplied to the result.
        :return:
            A numpy array on the master node, None on other nodes.
        """

        mat = np.matmul(self.data.T, dist_array.data)

        if factor is not None:
            mat *= factor

        mat = self.comm.reduce(mat, op=MPI.SUM, root=mpi_master())

        return mat

    def matmul_AtB_allreduce(self, dist_array, factor=None):
        """
        Computes matrix-matrix multiplication between self.T and a distributed
        array, and makes the result available on all nodes.

        :param dist_array:
            The distributed array.
        :param factor:
            The factor to be multiplied to the result.
        :return:
            A numpy array that is available on all nodes.
        """

        mat = np.matmul(self.data.T, dist_array.data)

        if factor is not None:
            mat *= factor

        mat = self.comm.allreduce(mat, op=MPI.SUM)

        return mat

    def matmul_AB_no_gather(self, array, factor=None):
        """
        Computes matrix-matrix multiplication between self and a numpy array
        that is available on all nodes.

        :param array:
            The numpy array.
        :param factor:
            The factor to be multiplied to the result.
        :return:
            A distributed array.
        """

        seg_mat = np.matmul(self.data, array)

        if factor is not None:
            seg_mat *= factor

        return DistributedArray(seg_mat, self.comm, distribute=False)

    def append(self, dist_array, axis=None):
        """
        Appends a distributed array to self.

        :param dist_array:
            The distributed array.
        :param axis:
            The axis parameter as in numpy.append.
        """

        self.data = np.append(self.data, dist_array.data, axis=axis)

    def get_column(self, col_index):
        """
        Gets a column as a distributed 1D array.

        :param col_index:
            The index of the column.

        :return:
            A distributed 1D array.
        """

        return DistributedArray(self.data[:, col_index],
                                self.comm,
                                distribute=False)

    @classmethod
    def read_from_hdf5_file(cls, fname, label, comm):
        """
        Reads an array from checkpoint file and returns a distributed array.

        :param fname:
            The name of the checkpoint file.
        :param label:
            The label for the array.
        :param comm:
            The communicator.

        :return:
            The distributed array.
        """

        rank = comm.Get_rank()
        nodes = comm.Get_size()

        if rank == mpi_master():
            hf = h5py.File(fname, 'r')

            dset = hf[label]

            ave, res = divmod(dset.shape[0], nodes)
            counts = [ave + 1 if p < res else ave for p in range(nodes)]
            displacements = [sum(counts[:p]) for p in range(nodes)]

            data = np.array(dset[0:1])
            batch_size = mpi_size_limit() // data.itemsize
            if data.ndim == 2 and data.shape[1] != 0:
                batch_size //= data.shape[1]
            batch_size = max(1, batch_size)

            pack_list = [batch_size] + counts
        else:
            pack_list = None

        pack_list = comm.bcast(pack_list, root=mpi_master())
        batch_size = pack_list[0]
        counts = pack_list[1:]

        if rank == mpi_master():
            data = np.array(dset[displacements[0]:displacements[0] + counts[0]])

            for i in range(1, nodes):
                for batch_start in range(0, counts[i], batch_size):
                    row_start = batch_start + displacements[i]
                    row_end = min(row_start + batch_size,
                                  counts[i] + displacements[i])
                    comm.send(np.array(dset[row_start:row_end]),
                              dest=i,
                              tag=i + batch_start)

            hf.close()
        else:
            data = None

            for batch_start in range(0, counts[rank], batch_size):
                recv_data = comm.recv(source=mpi_master(),
                                      tag=rank + batch_start)
                if data is None:
                    data = recv_data.copy()
                else:
                    data = np.vstack((data, recv_data))

        return cls(data, comm, distribute=False)

    def append_to_hdf5_file(self, fname, label):
        """
        Appends an array to checkpoint file.

        :param fname:
            The name of the checkpoint file.
        :param label:
            The label for the array.

        :return:
            The time spent in writing to checkpoint file.
        """

        counts = self.comm.allgather(self.shape(0))
        displacements = [sum(counts[:p]) for p in range(self.nodes)]

        batch_size = mpi_size_limit() // self.data.itemsize
        if self.data.ndim == 2 and self.shape(1) != 0:
            batch_size //= self.shape(1)
        batch_size = max(1, batch_size)

        t0 = tm.time()

        if self.rank == mpi_master():
            hf = h5py.File(fname, 'a')

            dset = hf.create_dataset(label, (sum(counts), self.shape(1)),
                                     dtype=self.data.dtype,
                                     compression='gzip')

            dset[displacements[0]:displacements[0] + counts[0]] = self.data[:]

            for i in range(1, self.nodes):
                for batch_start in range(0, counts[i], batch_size):
                    row_start = batch_start + displacements[i]
                    row_end = min(row_start + batch_size,
                                  counts[i] + displacements[i])
                    recv_data = self.comm.recv(source=i, tag=i + batch_start)
                    dset[row_start:row_end] = recv_data[:]

            hf.close()
        else:
            for batch_start in range(0, counts[self.rank], batch_size):
                row_start = batch_start
                row_end = min(row_start + batch_size, counts[self.rank])
                self.comm.send(self.data[row_start:row_end],
                               dest=mpi_master(),
                               tag=self.rank + batch_start)

        return tm.time() - t0
