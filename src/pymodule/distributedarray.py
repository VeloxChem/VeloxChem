from mpi4py import MPI
import numpy as np

from .veloxchemlib import mpi_master


class DistributedArray:
    """
    Impements distributed array.

    :param array:
        The numpy array stored on the master node.
    :param comm:
        The communicator.
    :param counts:
        The number of elements on each node.
    :param displacements:
        The displacement of element on each node.

    Instance variable
        - comm: The communicator.
        - data: The distributed array segment.
    """

    def __init__(self, array, comm, counts=None, displacements=None):
        """
        Initializes distributed array.
        """

        if comm.Get_rank() == mpi_master():
            n_elems = array.shape[0]
            n_ranks = comm.Get_size()

            if counts is None or displacements is None:
                ave, res = divmod(n_elems, n_ranks)
                counts = [ave + 1 if p < res else ave for p in range(n_ranks)]
                displacements = [sum(counts[:p]) for p in range(n_ranks)]

            if array.ndim == 1:
                array_list = [
                    array[displacements[p]:displacements[p] + counts[p]]
                    for p in range(n_ranks)
                ]

            elif array.ndim == 2:
                array_list = [
                    array[displacements[p]:displacements[p] + counts[p], :]
                    for p in range(n_ranks)
                ]
        else:
            array_list = None

        self.comm = comm
        self.data = comm.scatter(array_list, root=mpi_master())

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

        if self.comm.Get_rank() == mpi_master():
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
