from mpi4py import MPI
import numpy as np
import time as tm
import h5py

from .veloxchemlib import mpi_master


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

        if self.rank == mpi_master() and distribute:
            ave, res = divmod(array.shape[0], self.nodes)
            counts = [ave + 1 if p < res else ave for p in range(self.nodes)]
            displacements = [sum(counts[:p]) for p in range(self.nodes)]

            if array.ndim == 1:
                array_list = [
                    array[displacements[p]:displacements[p] + counts[p]]
                    for p in range(self.nodes)
                ]

            elif array.ndim == 2:
                array_list = [
                    array[displacements[p]:displacements[p] + counts[p], :]
                    for p in range(self.nodes)
                ]
        else:
            array_list = None

        if distribute:
            self.data = comm.scatter(array_list, root=mpi_master())
        else:
            self.data = array.copy()

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

            data = np.array(dset[displacements[0]:displacements[0] + counts[0]])

            for i in range(1, nodes):
                send_data = np.array(dset[displacements[i]:displacements[i] +
                                          counts[i]])
                comm.send(send_data, dest=i, tag=i)

            hf.close()
        else:
            data = comm.recv(source=mpi_master(), tag=rank)

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

        counts = self.comm.gather(self.data.shape[0], root=mpi_master())
        if self.rank == mpi_master():
            displacements = [sum(counts[:p]) for p in range(self.nodes)]

        t0 = tm.time()

        if self.rank == mpi_master():
            hf = h5py.File(fname, 'a')

            dset = hf.create_dataset(label, (sum(counts), self.data.shape[1]),
                                     dtype=self.data.dtype,
                                     compression='gzip')

            dset[displacements[0]:displacements[0] + counts[0]] = self.data[:]

            for i in range(1, self.nodes):
                data = self.comm.recv(source=i, tag=i)
                dset[displacements[i]:displacements[i] + counts[i]] = data[:]

            hf.close()
        else:
            self.comm.send(self.data, dest=mpi_master(), tag=self.rank)

        return tm.time() - t0
