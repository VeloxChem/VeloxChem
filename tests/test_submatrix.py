import numpy as np
import math as mt

from mpi4py import MPI
from veloxchem.submatrix import SubMatrix
from veloxchem.mpitools import is_master
from tester import Tester


class TestSubMatrix:

    def get_values(self):

        return [
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 1.2,
            2.2, 3.2, 4.2, 5.2, 6.2
        ]

    def get_np_values(self):

        return np.reshape(np.array(self.get_values()), (3, 6))

    def test_at_mutable(self):

        mat_a = SubMatrix(self.get_values(), [1, 5, 3, 6])

        mat_a.at(3, 9, 9.0)
        mat_a.at(2, 7, 6.0)

        mat_b = SubMatrix([
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.1, 2.1, 6.0, 4.1, 5.1, 6.1, 1.2,
            2.2, 3.2, 4.2, 9.0, 6.2
        ], [1, 5, 3, 6])

        Tester.compare_submatrices(mat_a, mat_b)

    def test_at_imutable(self):

        tol = 1.0e-12

        mat = SubMatrix(self.get_values(), [1, 5, 3, 6])

        assert mt.isclose(mat.at(1, 5), 1.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at(1, 6), 2.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at(1, 7), 3.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at(1, 8), 4.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at(1, 9), 5.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at(1, 10), 6.0, rel_tol=tol, abs_tol=tol)

        assert mt.isclose(mat.at(2, 5), 1.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at(2, 6), 2.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at(2, 7), 3.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at(2, 8), 4.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at(2, 9), 5.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at(2, 10), 6.1, rel_tol=tol, abs_tol=tol)

        assert mt.isclose(mat.at(3, 5), 1.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at(3, 6), 2.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at(3, 7), 3.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at(3, 8), 4.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at(3, 9), 5.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at(3, 10), 6.2, rel_tol=tol, abs_tol=tol)

    def test_set_offsets(self):

        mat_a = SubMatrix(self.get_values(), [0, 0, 3, 6])

        mat_b = SubMatrix(self.get_values(), [2, 7, 3, 6])

        mat_a.set_offsets(2, 7)

        Tester.compare_submatrices(mat_a, mat_b)

    def test_set_value(self):

        mat_a = SubMatrix(self.get_values(), [1, 5, 3, 6])

        mat_a.set_value(2, 4, 9.0)
        mat_a.set_value(1, 2, 6.0)

        mat_b = SubMatrix([
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.1, 2.1, 6.0, 4.1, 5.1, 6.1, 1.2,
            2.2, 3.2, 4.2, 9.0, 6.2
        ], [1, 5, 3, 6])

        Tester.compare_submatrices(mat_a, mat_b)

    def test_get_value(self):

        tol = 1.0e-12

        mat = SubMatrix(self.get_values(), [1, 5, 3, 6])

        assert mt.isclose(mat.get_value(0, 0), 1.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.get_value(0, 1), 2.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.get_value(0, 2), 3.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.get_value(0, 3), 4.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.get_value(0, 4), 5.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.get_value(0, 5), 6.0, rel_tol=tol, abs_tol=tol)

        assert mt.isclose(mat.get_value(1, 0), 1.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.get_value(1, 1), 2.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.get_value(1, 2), 3.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.get_value(1, 3), 4.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.get_value(1, 4), 5.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.get_value(1, 5), 6.1, rel_tol=tol, abs_tol=tol)

        assert mt.isclose(mat.get_value(2, 0), 1.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.get_value(2, 1), 2.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.get_value(2, 2), 3.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.get_value(2, 3), 4.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.get_value(2, 4), 5.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.get_value(2, 5), 6.2, rel_tol=tol, abs_tol=tol)

    def test_set_values(self):

        mat_a = SubMatrix([1, 5, 3, 6])

        mat_a.set_values(self.get_np_values())

        mat_b = SubMatrix(self.get_values(), [1, 5, 3, 6])

        Tester.compare_submatrices(mat_a, mat_b)

    def test_zero(self):

        mat_a = SubMatrix(self.get_values(), [1, 5, 3, 6])

        mat_b = SubMatrix([1, 5, 3, 6])

        mat_b.set_values(np.zeros((3, 6)))

        mat_a.zero()

        Tester.compare_submatrices(mat_a, mat_b)

    def test_to_numpy(self):

        tol = 1.0e-12

        mat = SubMatrix(self.get_values(), [1, 5, 3, 6])

        assert np.allclose(self.get_np_values(), mat.to_numpy(), tol, tol,
                           False)

    def test_get_dimensions(self):

        mat = SubMatrix(self.get_values(), [1, 5, 3, 6])

        assert mat.get_dimensions() == [1, 5, 3, 6]

    def test_offset_of_rows(self):

        mat = SubMatrix(self.get_values(), [1, 5, 3, 6])

        assert mat.offset_of_rows() == 1

    def test_offset_of_columns(self):

        mat = SubMatrix(self.get_values(), [1, 5, 3, 6])

        assert mat.offset_of_columns() == 5

    def test_number_of_rows(self):

        mat = SubMatrix(self.get_values(), [1, 5, 3, 6])

        assert mat.number_of_rows() == 3

    def test_number_of_columns(self):

        mat = SubMatrix(self.get_values(), [1, 5, 3, 6])

        assert mat.number_of_columns() == 6

    def test_mpi_bcast(self):

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

        if is_master(rank):
            mat_a = SubMatrix(self.get_values(), [1, 5, 3, 6])
        else:
            mat_a = None

        mat_a = SubMatrix.bcast(mat_a, rank, comm)

        mat_b = SubMatrix(self.get_values(), [1, 5, 3, 6])

        Tester.compare_submatrices(mat_a, mat_b)

    def test_mpi_reduce(self):

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nodes = comm.Get_size()

        if is_master(rank):
            mat_a = SubMatrix(self.get_values(), [1, 5, 3, 6])
        elif rank == 1:
            mat_a = SubMatrix([
                1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 1.4, 2.4, 3.4, 4.4, 5.4, 6.4, 1.5,
                2.5, 3.5, 4.5, 5.5, 6.5
            ], [1, 5, 3, 6])
        else:
            mat_a = SubMatrix([1, 5, 3, 6])

        mat_a = SubMatrix.reduce(mat_a, rank, comm)

        if is_master(rank):
            if nodes == 1:
                mat_b = SubMatrix(self.get_values(), [1, 5, 3, 6])
            else:
                mat_b = SubMatrix([
                    2.2, 4.2, 6.2, 8.2, 10.2, 12.2, 2.5, 4.5, 6.5, 8.5, 10.5,
                    12.5, 2.7, 4.7, 6.7, 8.7, 10.7, 12.7
                ], [1, 5, 3, 6])
            Tester.compare_submatrices(mat_a, mat_b)
