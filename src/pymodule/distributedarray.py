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
import time as tm
import h5py

from .veloxchemlib import mpi_master
from .veloxchemlib import mpi_size_limit
from .veloxchemlib import matmul_gpu


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

    def __init__(self, array, comm, distribute=True, root=None):
        """
        Initializes distributed array.
        """

        if root is None:
            root = mpi_master()

        self.comm = comm
        self.rank = comm.Get_rank()
        self.nodes = comm.Get_size()

        self.data = None

        if not distribute:
            self.data = array.copy()
            return

        if self.rank == root:
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
        batch_size = comm.bcast(batch_size, root=root)
        n_total = comm.bcast(n_total, root=root)

        if n_total == 0:
            counts = comm.bcast(counts, root=root)
            self.data = np.zeros((counts[self.rank], 0))
            return

        for batch_start in range(0, n_total, batch_size):
            batch_end = min(batch_start + batch_size, n_total)

            if self.rank == root:
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

            recvbuf = comm.scatter(array_list, root=root)
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

    def get_full_vector(self, col=None, root=None):
        """
        Gets a full column vector from a distributed array.

        :param: col:
            The column index (used only when self.data.ndim is 2).
        :param: root:
            The root rank.
        :return:
            The full vector on the master node, None on other nodes.
        """

        if root is None:
            root = mpi_master()

        data = None

        if self.data.ndim == 1:
            data = self.comm.gather(self.data, root=root)
        elif self.data.ndim == 2 and col is not None:
            data = self.comm.gather(self.data[:, col], root=root)

        if self.rank == root:
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

        if (self.data.ndim == 2 and dist_array.data.ndim == 2 and
                self.data.dtype == np.float64 and
                dist_array.data.dtype == np.float64):
            mat = matmul_gpu(self.data.T, dist_array.data)
        else:
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

        if (self.data.ndim == 2 and dist_array.data.ndim == 2 and
                self.data.dtype == np.float64 and
                dist_array.data.dtype == np.float64):
            mat = matmul_gpu(self.data.T, dist_array.data)
        else:
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

        if (self.data.ndim == 2 and array.ndim == 2 and
                self.data.dtype == np.float64 and array.dtype == np.float64):
            seg_mat = matmul_gpu(self.data, array)
        else:
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

        h5cfg = h5py.get_config()

        if hasattr(h5cfg, 'mpi') and h5cfg.mpi:
            return cls.read_from_hdf5_file_parallel(fname, label, comm)
        else:
            return cls.read_from_hdf5_file_serial(fname, label, comm)

    @classmethod
    def read_from_hdf5_file_parallel(cls, fname, label, comm):
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

        hf = h5py.File(fname, 'r', driver='mpio', comm=comm)
        dset = hf[label]

        ave, res = divmod(dset.shape[0], nodes)
        counts = [ave + 1 if p < res else ave for p in range(nodes)]
        displacements = [sum(counts[:p]) for p in range(nodes)]

        row_start = displacements[rank]
        row_end = row_start + counts[rank]
        data = np.array(dset[row_start:row_end, :])

        hf.close()

        return cls(data, comm, distribute=False)

    @classmethod
    def read_from_hdf5_file_serial(cls, fname, label, comm):
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
            n_total = dset.shape[1]
            hf.close()
        else:
            counts = None
            displacements = None
            n_total = None
        counts, displacements, n_total = comm.bcast(
            (counts, displacements, n_total), root=mpi_master())

        data = None

        if n_total == 0:
            data = np.zeros((counts[rank], 0))
        else:
            hf = h5py.File(fname, 'r')
            dset = hf[label]
            row_start = displacements[rank]
            row_end = row_start + counts[rank]
            data = np.array(dset[row_start:row_end, :])
            hf.close()
            comm.barrier()

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

        h5cfg = h5py.get_config()

        if hasattr(h5cfg, 'mpi') and h5cfg.mpi:
            return self.append_to_hdf5_file_parallel(fname, label)
        else:
            return self.append_to_hdf5_file_serial(fname, label)

    def append_to_hdf5_file_parallel(self, fname, label):
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

        counts = self.comm.allgather(self.shape(0))
        displacements = [sum(counts[:p]) for p in range(self.nodes)]

        hf = h5py.File(fname, 'a', driver='mpio', comm=self.comm)
        shape = (sum(counts), self.shape(1))
        dset = hf.create_dataset(label, shape=shape, dtype=self.data.dtype)

        row_start = displacements[self.rank]
        row_end = row_start + counts[self.rank]
        dset[row_start:row_end, :] = self.data[:, :]

        hf.close()

        return tm.time() - t0

    def append_to_hdf5_file_serial(self, fname, label):
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

        counts = self.comm.allgather(self.shape(0))

        if self.rank == mpi_master():
            displacements = [sum(counts[:p]) for p in range(self.nodes)]
            n_total = self.shape(1)
        else:
            displacements = None
            n_total = None
        displacements, n_total = self.comm.bcast((displacements, n_total),
                                                 root=mpi_master())

        if self.rank == mpi_master():
            hf = h5py.File(fname, 'a')
            shape = (sum(counts), self.shape(1))
            dset = hf.create_dataset(label, shape=shape, dtype=self.data.dtype)
            hf.close()
        self.comm.barrier()

        if n_total > 0:
            for i in range(self.nodes):
                if self.rank == i:
                    hf = h5py.File(fname, 'r+')
                    dset = hf[label]
                    row_start = displacements[self.rank]
                    row_end = row_start + counts[self.rank]
                    dset[row_start:row_end, :] = self.data[:, :]
                    hf.close()
                self.comm.barrier()

        return tm.time() - t0
