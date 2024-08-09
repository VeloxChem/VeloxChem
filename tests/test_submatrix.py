import pickle
import numpy as np
import math as mt

from mpi4py import MPI

from veloxchem import SubMatrix


class TestSubMatrix:

    def get_values(self):

        return [
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 1.2,
            2.2, 3.2, 4.2, 5.2, 6.2
        ]

    def get_np_values(self):

        return np.reshape(np.array(self.get_values()), (3, 6))

    def test_pickle(self):

        mat_a = SubMatrix(self.get_values(), [1, 5, 3, 6])
        bobj = pickle.dumps(mat_a)
        mat_b = pickle.loads(bobj)
        assert mat_a == mat_b

    def test_sum_op(self):

        mat_a = SubMatrix(self.get_values(), [1, 5, 3, 6])
        mat_a.scale(0.9)
        mat_b = SubMatrix(self.get_values(), [1, 5, 3, 6])
        mat_c = SubMatrix([1, 5, 3, 6])
        mat_c.set_values(1.9 * self.get_np_values())
        assert mat_c == mat_a + mat_b

    def test_subscript_mutable(self):

        mat_a = SubMatrix(self.get_values(), [1, 5, 3, 6])
        mat_a[[3, 9]] = 9.0
        mat_a[[2, 7]] = 6.0
        mat_b = SubMatrix([
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.1, 2.1, 6.0, 4.1, 5.1, 6.1, 1.2,
            2.2, 3.2, 4.2, 9.0, 6.2
        ], [1, 5, 3, 6])
        assert mat_a == mat_b

    def test_subsript_imutable(self):

        tol = 1.0e-12

        mat = SubMatrix(self.get_values(), [1, 5, 3, 6])

        assert mt.isclose(mat[[1, 5]], 1.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat[[1, 6]], 2.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat[[1, 7]], 3.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat[[1, 8]], 4.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat[[1, 9]], 5.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat[[1, 10]], 6.0, rel_tol=tol, abs_tol=tol)

        assert mt.isclose(mat[[2, 5]], 1.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat[[2, 6]], 2.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat[[2, 7]], 3.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat[[2, 8]], 4.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat[[2, 9]], 5.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat[[2, 10]], 6.1, rel_tol=tol, abs_tol=tol)

        assert mt.isclose(mat[[3, 5]], 1.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat[[3, 6]], 2.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat[[3, 7]], 3.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat[[3, 8]], 4.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat[[3, 9]], 5.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat[[3, 10]], 6.2, rel_tol=tol, abs_tol=tol)

    def test_set_offsets(self):

        mat_a = SubMatrix(self.get_values(), [0, 0, 3, 6])
        mat_b = SubMatrix(self.get_values(), [2, 7, 3, 6])
        mat_a.set_offsets([2, 7])
        assert mat_a == mat_b

    def test_at_mutable(self):

        mat_a = SubMatrix(self.get_values(), [1, 5, 3, 6])
        mat_a.at([2, 4], 9.0)
        mat_a.at([1, 2], 6.0)
        mat_b = SubMatrix([
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.1, 2.1, 6.0, 4.1, 5.1, 6.1, 1.2,
            2.2, 3.2, 4.2, 9.0, 6.2
        ], [1, 5, 3, 6])
        assert mat_a == mat_b

    def test_at_imutable(self):

        tol = 1.0e-12

        mat = SubMatrix(self.get_values(), [1, 5, 3, 6])

        assert mt.isclose(mat.at([0, 0]), 1.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at([0, 1]), 2.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at([0, 2]), 3.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at([0, 3]), 4.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at([0, 4]), 5.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at([0, 5]), 6.0, rel_tol=tol, abs_tol=tol)

        assert mt.isclose(mat.at([1, 0]), 1.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at([1, 1]), 2.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at([1, 2]), 3.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at([1, 3]), 4.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at([1, 4]), 5.1, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at([1, 5]), 6.1, rel_tol=tol, abs_tol=tol)

        assert mt.isclose(mat.at([2, 0]), 1.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at([2, 1]), 2.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at([2, 2]), 3.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at([2, 3]), 4.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at([2, 4]), 5.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(mat.at([2, 5]), 6.2, rel_tol=tol, abs_tol=tol)

    def test_set_values(self):

        mat_a = SubMatrix([1, 5, 3, 6])
        mat_a.set_values(self.get_np_values())
        mat_b = SubMatrix(self.get_values(), [1, 5, 3, 6])
        assert mat_a == mat_b

    def test_zero(self):

        mat_a = SubMatrix(self.get_values(), [1, 5, 3, 6])
        mat_b = SubMatrix([1, 5, 3, 6])
        mat_b.set_values(np.zeros((3, 6)))
        mat_a.zero()
        assert mat_a == mat_b

    def test_scale(self):

        mat_a = SubMatrix(self.get_values(), [1, 5, 3, 6])
        mat_a.scale(0.9)
        mat_b = SubMatrix([1, 5, 3, 6])
        mat_b.set_values(0.9 * self.get_np_values())
        assert mat_a == mat_b

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

    def test_number_of_elements(self):

        mat = SubMatrix(self.get_values(), [1, 5, 3, 6])
        assert mat.number_of_elements() == 18

    def test_is_square(self):

        mat = SubMatrix([1, 5, 3, 6])
        assert mat.is_square() == False
        mat = SubMatrix([1, 5, 3, 3])
        assert mat.is_square() == True

    def test_mpi_bcast(self):

        comm = MPI.COMM_WORLD

        mat_a = None
        if comm.Get_rank() == 0:
            mat_a = SubMatrix(self.get_values(), [1, 5, 3, 6])
        mat_a = comm.bcast(mat_a)
        mat_b = SubMatrix(self.get_values(), [1, 5, 3, 6])
        assert mat_a == mat_b

    def test_mpi_reduce(self):

        comm = MPI.COMM_WORLD

        mat_a = SubMatrix(self.get_values(), [1, 5, 3, 6])
        if comm.Get_rank() == 0:
            mat_a.scale(0.82)
        mat_r = SubMatrix.reduce(mat_a, comm, 0)
        if comm.Get_rank() == 0:
            mat_b = SubMatrix(self.get_values(), [1, 5, 3, 6])
            mat_b.scale(0.82 + comm.Get_size() - 1)
            assert mat_r == mat_b
