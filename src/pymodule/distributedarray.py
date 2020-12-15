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

        if not distribute:
            self.data = array.copy()
            return

        if self.rank == mpi_master():
            dtype = array.dtype
            if array.ndim == 1:
                ncols = 1
            elif array.ndim == 2:
                ncols = array.shape[1]

            ave, res = divmod(array.shape[0], self.nodes)
            counts = [ave + 1 if p < res else ave for p in range(self.nodes)]
            if ncols != 0:  # note that ncols can be 0
                counts = [x * ncols for x in counts]
            displacements = [sum(counts[:p]) for p in range(self.nodes)]
            pack_list = [dtype, ncols, array.ndim] + counts
        else:
            counts = None
            displacements = None
            pack_list = None

        pack_list = self.comm.bcast(pack_list, root=mpi_master())
        dtype, ncols, ndim = pack_list[:3]
        counts = pack_list[3:]

        if ndim == 1:
            self.data = np.zeros(counts[self.rank], dtype=dtype)
        elif ndim == 2:
            local_rows = counts[self.rank]
            if ncols != 0:  # note that ncols can be 0
                local_rows = counts[self.rank] // ncols
            self.data = np.zeros((local_rows, ncols), dtype=dtype)

        if ncols == 0:
            return

        if dtype == np.int:
            mpi_data_type = MPI.INT
        elif dtype == np.float:
            mpi_data_type = MPI.DOUBLE
        elif dtype == np.complex:
            mpi_data_type = MPI.C_DOUBLE_COMPLEX

        self.comm.Scatterv([array, counts, displacements, mpi_data_type],
                           self.data,
                           root=mpi_master())

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

    def matmul_AB(self, array, factor=None):
        """
        Computes matrix-matrix multiplication between self and a numpy array
        that is available on all nodes.

        :param array:
            The numpy array.
        :param factor:
            The factor to be multiplied to the result.
        :return:
            A numpy array on the master node, None on other nodes.
        """

        seg_mat = np.matmul(self.data, array)

        if factor is not None:
            seg_mat *= factor

        seg_mat = self.comm.gather(seg_mat, root=mpi_master())

        if self.rank == mpi_master():
            full_shape_0 = sum([m.shape[0] for m in seg_mat])

            if array.ndim == 1:
                return np.hstack(seg_mat).reshape(full_shape_0)

            elif array.ndim == 2:
                return np.vstack(seg_mat).reshape(full_shape_0, array.shape[1])

        else:
            return None

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
        else:
            counts = None
        counts = comm.bcast(counts, root=mpi_master())

        if rank == mpi_master():
            data = np.array(dset[displacements[0]:displacements[0] + counts[0]])

            for i in range(1, nodes):
                batch_size = mpi_size_limit() // data.itemsize
                if data.ndim == 2 and data.shape[1] != 0:
                    batch_size //= data.shape[1]
                batch_size = min(max(1, batch_size), counts[i])
                for row_start in range(displacements[i],
                                       displacements[i] + counts[i],
                                       batch_size):
                    row_end = min(row_start + batch_size,
                                  displacements[i] + counts[i])
                    comm.send(np.array(dset[row_start:row_end]), dest=i, tag=i)

            hf.close()
        else:
            data = None

            row_count = 0
            while row_count < counts[rank]:
                recv_data = comm.recv(source=mpi_master(), tag=rank)
                if data is None:
                    data = recv_data.copy()
                else:
                    data = np.vstack((data, recv_data))
                row_count += recv_data.shape[0]

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

        counts = self.comm.gather(self.shape(0), root=mpi_master())
        if self.rank == mpi_master():
            displacements = [sum(counts[:p]) for p in range(self.nodes)]

        t0 = tm.time()

        if self.rank == mpi_master():
            hf = h5py.File(fname, 'a')

            dset = hf.create_dataset(label, (sum(counts), self.shape(1)),
                                     dtype=self.data.dtype,
                                     compression='gzip')

            dset[displacements[0]:displacements[0] + counts[0]] = self.data[:]

            for i in range(1, self.nodes):
                row_count = 0
                while row_count < counts[i]:
                    recv_data = self.comm.recv(source=i, tag=i)
                    row_start = displacements[i] + row_count
                    row_end = row_start + recv_data.shape[0]
                    dset[row_start:row_end] = recv_data[:]
                    row_count += recv_data.shape[0]

            hf.close()
        else:
            batch_size = mpi_size_limit() // self.data.itemsize
            if self.data.ndim == 2 and self.shape(1) != 0:
                batch_size //= self.shape(1)
            batch_size = min(max(1, batch_size), self.shape(0))
            for row_start in range(0, self.shape(0), batch_size):
                row_end = min(row_start + batch_size, self.shape(0))
                self.comm.send(self.data[row_start:row_end],
                               dest=mpi_master(),
                               tag=self.rank)

        return tm.time() - t0
