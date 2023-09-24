from mpi4py import MPI
import numpy as np

from veloxchem.veloxchemlib import DenseMatrix
from veloxchem.veloxchemlib import (mpi_master, bcast_scalar,
                                    bcast_dense_matrix, scatter_vector,
                                    gather_dense_matrices_by_columns)


class TestMPI:

    def test_bcast_scalar(self):

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

        ref_val, other_val = 99, 1
        a = ref_val if rank == 0 else other_val
        a = bcast_scalar(a, comm)
        assert a == ref_val

        ref_val, other_val = 0.99, 0.01
        a = ref_val if rank == 0 else other_val
        a = bcast_scalar(a, comm)
        assert a == ref_val

    def test_bcast_dense_matrix(self):

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

        ref_mat = np.arange(12.).reshape(4, 3) + 99.99
        other_mat = np.arange(20.).reshape(4, 5) + 0.01
        a = DenseMatrix(ref_mat) if rank == 0 else DenseMatrix(other_mat)
        a = bcast_dense_matrix(a, comm)
        assert np.max(np.abs(a.to_numpy() - ref_mat)) < 1.0e-13

    def test_scatter_vector(self):

        comm = MPI.COMM_WORLD
        nodes = comm.Get_size()

        ref_vec = np.arange(97)

        scatter_vec = scatter_vector(ref_vec, comm)

        ave, rem = divmod(ref_vec.size, nodes)
        counts = [(ave + 1 if p < rem else ave) for p in range(nodes)]
        displs = [sum(counts[:p]) for p in range(nodes)]
        data_list = [
            ref_vec[displs[p]:displs[p] + counts[p]] for p in range(nodes)
        ]
        scatter_data = comm.scatter(data_list, root=mpi_master())

        assert scatter_vec == list(scatter_data)

    def test_gather_dense_matrix(self):

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nodes = comm.Get_size()

        ref_size = 97
        ave, rem = divmod(ref_size, nodes)
        counts = [(ave + 1 if p < rem else ave) for p in range(nodes)]
        mat = np.arange(5.0 * counts[rank]).reshape(5, -1) + 0.1 * rank

        gathered_mat = gather_dense_matrices_by_columns(DenseMatrix(mat), comm)

        gathered_data = comm.gather(mat, root=mpi_master())

        if rank == mpi_master():
            assert np.max(
                np.abs(gathered_mat.to_numpy() -
                       np.hstack(gathered_data))) < 1.0e-13
