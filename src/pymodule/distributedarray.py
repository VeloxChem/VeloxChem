#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
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
            # determine counts and displacements for scatter
            ave, res = divmod(array.shape[0], self.nodes)
            counts = [ave + 1 if p < res else ave for p in range(self.nodes)]
            displacements = [sum(counts[:p]) for p in range(self.nodes)]

            # determine batch size for number of columns
            batch_size = mpi_size_limit() // (array.itemsize * array.shape[0])
            batch_size = max(1, batch_size)
            if array.ndim == 1:
                n_total = 1
            elif array.ndim == 2:
                n_total = array.shape[1]
        else:
            counts = None
            batch_size = None
            n_total = None
        batch_size = comm.bcast(batch_size, root=mpi_master())
        n_total = comm.bcast(n_total, root=mpi_master())

        if n_total == 0:
            counts = comm.bcast(counts, root=mpi_master())
            self.data = np.zeros((counts[self.rank], 0))
            return

        for batch_start in range(0, n_total, batch_size):
            batch_end = min(batch_start + batch_size, n_total)

            if self.rank == mpi_master():
                if array.ndim == 1:
                    array_list = [
                        array[displacements[p]:displacements[p] + counts[p]]
                        for p in range(self.nodes)
                    ]
                elif array.ndim == 2:
                    array_list = [
                        array[displacements[p]:displacements[p] + counts[p],
                              batch_start:batch_end] for p in range(self.nodes)
                    ]
            else:
                array_list = None

            recvbuf = comm.scatter(array_list, root=mpi_master())
            if self.data is None:
                self.data = recvbuf
            else:
                self.data = np.hstack((self.data, recvbuf))

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

            # determine batch size for number of columns
            itemsize = np.array(dset[0]).itemsize
            batch_size = mpi_size_limit() // (itemsize * dset.shape[0])
            batch_size = max(1, batch_size)
            n_total = dset.shape[1]
        else:
            counts = None
            batch_size = None
            n_total = None
        batch_size = comm.bcast(batch_size, root=mpi_master())
        n_total = comm.bcast(n_total, root=mpi_master())

        data = None

        if n_total == 0:
            counts = comm.bcast(counts, root=mpi_master())
            data = np.zeros((counts[rank], 0))

        for batch_start in range(0, n_total, batch_size):
            batch_end = min(batch_start + batch_size, n_total)

            if rank == mpi_master():
                array_list = [
                    np.array(dset[displacements[i]:displacements[i] + counts[i],
                                  batch_start:batch_end]) for i in range(nodes)
                ]

            else:
                array_list = None

            recvbuf = comm.scatter(array_list, root=mpi_master())
            if data is None:
                data = recvbuf
            else:
                data = np.hstack((data, recvbuf))

        if rank == mpi_master():
            hf.close()

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

        t0 = tm.time()

        counts = self.comm.gather(self.shape(0))

        if self.rank == mpi_master():
            displacements = [sum(counts[:p]) for p in range(self.nodes)]

            batch_size = mpi_size_limit() // (self.data.itemsize * sum(counts))
            batch_size = max(1, batch_size)
            n_total = self.shape(1)
        else:
            batch_size = None
            n_total = None
        batch_size = self.comm.bcast(batch_size, root=mpi_master())
        n_total = self.comm.bcast(n_total, root=mpi_master())

        if self.rank == mpi_master():
            hf = h5py.File(fname, 'a')

            dset = hf.create_dataset(label, (sum(counts), self.shape(1)),
                                     dtype=self.data.dtype,
                                     compression='gzip')

        for batch_start in range(0, n_total, batch_size):
            batch_end = min(batch_start + batch_size, n_total)

            array_list = self.comm.gather(self.data[:, batch_start:batch_end],
                                          root=mpi_master())

            if self.rank == mpi_master():
                for i in range(self.nodes):
                    dset[displacements[i]:displacements[i] + counts[i],
                         batch_start:batch_end] = array_list[i][:, :]

        if self.rank == mpi_master():
            hf.close()

        return tm.time() - t0
