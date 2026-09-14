import numpy as np
from mpi4py import MPI

from veloxchem.veloxchemlib import PackedMatrix, mat_t

# NOTE: the shapes cover a triangular matrix, whose storage is not its shape, and
# a rectangular one, whose rows and columns differ. Run under mpirun to exercise
# the communication; on one rank the methods take their shortcuts and only those
# are tested.

SHAPES = [(6, 6, mat_t.symmetric), (6, 6, mat_t.general), (7, 4, mat_t.general)]


class TestPackedMatrixMpi:

    def test_values_view_writes_into_the_matrix(self):

        matrix = PackedMatrix(5, 5, mat_t.symmetric)
        matrix.zero()

        view = matrix.values_view()
        view[0] = -1.5

        assert matrix.at(0, 0) == -1.5
        assert not view.flags.owndata
        assert view.size == matrix.number_of_elements()

    def test_broadcast_reaches_every_rank(self):

        comm = MPI.COMM_WORLD

        for nrows, ncols, mtype in SHAPES:
            matrix = PackedMatrix(nrows, ncols, mtype)

            if comm.Get_rank() == 0:
                view = matrix.values_view()
                view[:] = np.arange(view.size, dtype=float) + 1.0
            else:
                matrix.zero()

            matrix = matrix.broadcast(comm)

            expected = np.arange(matrix.number_of_elements(), dtype=float) + 1.0

            assert np.array_equal(matrix.values_view(), expected)

    def test_reduce_sums_onto_the_root(self):

        comm = MPI.COMM_WORLD
        rank, nodes = comm.Get_rank(), comm.Get_size()

        for nrows, ncols, mtype in SHAPES:
            matrix = PackedMatrix(nrows, ncols, mtype)
            matrix.values_view()[:] = float(rank + 1)

            total = matrix.reduce(comm)

            if rank == 0:
                assert np.allclose(total.values_view(),
                                   float(nodes * (nodes + 1) // 2))
            else:
                assert total is None

            # a reduction leaves what it reduced alone
            assert np.allclose(matrix.values_view(), float(rank + 1))

    def test_allreduce_sums_onto_every_rank(self):

        comm = MPI.COMM_WORLD
        rank, nodes = comm.Get_rank(), comm.Get_size()

        for nrows, ncols, mtype in SHAPES:
            matrix = PackedMatrix(nrows, ncols, mtype)
            matrix.values_view()[:] = float(rank + 1)

            every = matrix.allreduce(comm)

            assert np.allclose(every.values_view(),
                               float(nodes * (nodes + 1) // 2))
