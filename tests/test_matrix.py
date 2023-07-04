from mpi4py import MPI

from veloxchem.matrix import Matrix
from veloxchem.submatrix import SubMatrix
from veloxchem.mpitools import is_master
from veloxchem.veloxchemlib import mat_t
from tester import Tester


class TestMatrix:

    def get_mat_ss(self):

        return SubMatrix([1.0, 2.0, 3.0, 1.1, 2.1, 3.1, 1.2, 2.2, 3.2],
                         [0, 0, 3, 3])

    def get_mat_ss_r1(self):

        return SubMatrix([1.2, 2.2, 3.2, 1.5, 2.5, 3.5, 1.3, 2.3, 3.3],
                         [0, 0, 3, 3])

    def get_mat_ss_f(self):

        return SubMatrix([2.2, 4.2, 6.2, 2.6, 4.6, 6.6, 2.5, 4.5, 6.5],
                         [0, 0, 3, 3])

    def get_mat_sp(self):

        return SubMatrix([
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 1.2,
            2.2, 3.2, 4.2, 5.2, 6.2
        ], [0, 3, 3, 6])

    def get_mat_sp_r1(self):

        return SubMatrix([
            1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 1.2, 2.2, 3.2, 4.2, 5.2, 6.2, 1.3,
            2.3, 3.3, 4.3, 5.3, 6.3
        ], [0, 3, 3, 6])

    def get_mat_sp_f(self):

        return SubMatrix([
            2.1, 4.1, 6.1, 8.1, 10.1, 12.1, 2.3, 4.3, 6.3, 8.3, 10.3, 12.3, 2.5,
            4.5, 6.5, 8.5, 10.5, 12.5
        ], [0, 3, 3, 6])

    def get_mat_pp(self):

        return SubMatrix([
            2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 2.2,
            3.2, 4.2, 5.2, 6.2, 7.2, 2.3, 3.3, 4.3, 5.3, 6.3, 7.3, 2.4, 3.4,
            4.4, 5.4, 6.4, 7.4, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5
        ], [3, 3, 6, 6])

    def get_mat_pp_r1(self):

        return SubMatrix([
            2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 2.1,
            3.1, 4.1, 5.1, 6.1, 7.1, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 2.1, 3.1,
            4.1, 5.1, 6.1, 7.1, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0
        ], [3, 3, 6, 6])

    def get_mat_pp_f(self):

        return SubMatrix([
            4.1, 6.1, 8.1, 10.1, 12.1, 14.1, 4.1, 6.1, 8.1, 10.1, 12.1, 14.1,
            4.3, 6.3, 8.3, 10.3, 12.3, 14.3, 4.3, 6.3, 8.3, 10.3, 12.3, 14.3,
            4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5
        ], [3, 3, 6, 6])

    def get_mat_ps_symm(self):

        return SubMatrix([
            1.0, 1.1, 1.2, 2.0, 2.1, 2.2, 3.0, 3.1, 3.2, 4.0, 4.1, 4.2, 5.0,
            5.1, 5.2, 6.0, 6.1, 6.2
        ], [3, 0, 6, 3])

    def get_mat_ps_antisymm(self):

        return SubMatrix([
            -1.0, -1.1, -1.2, -2.0, -2.1, -2.2, -3.0, -3.1, -3.2, -4.0, -4.1,
            -4.2, -5.0, -5.1, -5.2, -6.0, -6.1, -6.2
        ], [3, 0, 6, 3])

    def get_mat_ps(self):

        return SubMatrix([
            1.5, 1.6, 1.7, 2.5, 2.6, 2.7, 3.5, 3.6, 3.7, 4.5, 4.6, 4.7, 5.5,
            5.6, 5.7, 6.5, 6.6, 6.7
        ], [3, 0, 6, 3])

    def get_mat_full_symm(self):

        return SubMatrix([
            1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.1, 2.1, 3.1, 1.1,
            2.1, 3.1, 4.1, 5.1, 6.1, 1.2, 2.2, 3.2, 1.2, 2.2, 3.2, 4.2, 5.2,
            6.2, 1.0, 1.1, 1.2, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 2.0, 2.1, 2.2,
            2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 3.0, 3.1, 3.2, 2.2, 3.2, 4.2, 5.2,
            6.2, 7.2, 4.0, 4.1, 4.2, 2.3, 3.3, 4.3, 5.3, 6.3, 7.3, 5.0, 5.1,
            5.2, 2.4, 3.4, 4.4, 5.4, 6.4, 7.4, 6.0, 6.1, 6.2, 2.5, 3.5, 4.5,
            5.5, 6.5, 7.5
        ], [0, 0, 9, 9])

    def get_mat_full_antisymm(self):

        return SubMatrix([
            1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.1, 2.1, 3.1, 1.1,
            2.1, 3.1, 4.1, 5.1, 6.1, 1.2, 2.2, 3.2, 1.2, 2.2, 3.2, 4.2, 5.2,
            6.2, -1.0, -1.1, -1.2, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, -2.0, -2.1,
            -2.2, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, -3.0, -3.1, -3.2, 2.2, 3.2, 4.2,
            5.2, 6.2, 7.2, -4.0, -4.1, -4.2, 2.3, 3.3, 4.3, 5.3, 6.3, 7.3, -5.0,
            -5.1, -5.2, 2.4, 3.4, 4.4, 5.4, 6.4, 7.4, -6.0, -6.1, -6.2, 2.5,
            3.5, 4.5, 5.5, 6.5, 7.5
        ], [0, 0, 9, 9])

    def get_mat_full_gen(self):

        return SubMatrix([
            1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.1, 2.1, 3.1, 1.1,
            2.1, 3.1, 4.1, 5.1, 6.1, 1.2, 2.2, 3.2, 1.2, 2.2, 3.2, 4.2, 5.2,
            6.2, 1.5, 1.6, 1.7, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 2.5, 2.6, 2.7,
            2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 3.5, 3.6, 3.7, 2.2, 3.2, 4.2, 5.2,
            6.2, 7.2, 4.5, 4.6, 4.7, 2.3, 3.3, 4.3, 5.3, 6.3, 7.3, 5.5, 5.6,
            5.7, 2.4, 3.4, 4.4, 5.4, 6.4, 7.4, 6.5, 6.6, 6.7, 2.5, 3.5, 4.5,
            5.5, 6.5, 7.5
        ], [0, 0, 9, 9])

    def test_add(self):

        mat = Matrix()
        mat.set_type(mat_t.symm)
        mat.add(self.get_mat_ss(), (0, 0))
        mat.add(self.get_mat_sp(), (0, 1))
        mat.add(self.get_mat_pp(), (1, 1))

        Tester.compare_submatrices(mat.get_submatrix((0, 0)), self.get_mat_ss())
        Tester.compare_submatrices(mat.get_submatrix((0, 1)), self.get_mat_sp())
        Tester.compare_submatrices(mat.get_submatrix((1, 1)), self.get_mat_pp())

    def test_set_and_get_type(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (1, 0): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symm)

        assert mat.get_type() == mat_t.symm

        mat.set_type(mat_t.antisymm)

        assert mat.get_type() == mat_t.antisymm

    def test_zero(self):

        smat_ss = self.get_mat_ss()
        smat_sp = self.get_mat_sp()
        smat_pp = self.get_mat_pp()

        mat_a = Matrix()
        mat_a.set_type(mat_t.symm)
        mat_a.add(smat_ss, (0, 0))
        mat_a.add(smat_sp, (0, 1))
        mat_a.add(smat_pp, (1, 1))

        smat_ss.zero()
        smat_sp.zero()
        smat_pp.zero()

        mat_b = Matrix()
        mat_b.set_type(mat_t.symm)
        mat_b.add(smat_ss, (0, 0))
        mat_b.add(smat_sp, (0, 1))
        mat_b.add(smat_pp, (1, 1))

        mat_a.zero()

        Tester.compare_matrices(mat_a, mat_b)

    def test_get_angular_pairs(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symm)

        assert mat.get_angular_pairs() == [(0, 0), (0, 1), (1, 1)]

    def test_get_submatrix(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symm)

        Tester.compare_submatrices(mat.get_submatrix((0, 0)), self.get_mat_ss())
        Tester.compare_submatrices(mat.get_submatrix((0, 1)), self.get_mat_sp())
        Tester.compare_submatrices(mat.get_submatrix((1, 1)), self.get_mat_pp())

    def test_is_angular_order(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symm)

        assert mat.is_angular_order((0, 0))
        assert mat.is_angular_order((0, 1))
        assert mat.is_angular_order((1, 1))

        assert mat.is_angular_order((1, 0)) is False
        assert mat.is_angular_order((1, 2)) is False
        assert mat.is_angular_order((2, 1)) is False
        assert mat.is_angular_order((2, 2)) is False

    def test_is_number_of_rows(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symm)

        assert mat.number_of_rows() == 9

        mat = Matrix({
            (0, 0): self.get_mat_ss(),
            (0, 1): self.get_mat_sp()
        }, mat_t.gen)

        assert mat.number_of_rows() == 3

    def test_is_number_of_columns(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symm)

        assert mat.number_of_columns() == 9

        mat = Matrix({
            (0, 0): self.get_mat_ss(),
            (0, 1): self.get_mat_sp()
        }, mat_t.gen)

        assert mat.number_of_columns() == 9

    def test_get_full_matrix_symm(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symm)

        Tester.compare_submatrices(mat.get_full_matrix(),
                                   self.get_mat_full_symm())

    def test_get_full_matrix_antisymm(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.antisymm)

        Tester.compare_submatrices(mat.get_full_matrix(),
                                   self.get_mat_full_antisymm())

    def test_get_full_matrix_gen(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 0): self.get_mat_ps(),
                (1, 1): self.get_mat_pp()
            }, mat_t.gen)

        Tester.compare_submatrices(mat.get_full_matrix(),
                                   self.get_mat_full_gen())

    def test_mpi_bcast(self):

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()

        if is_master(rank):
            mat_a = Matrix(
                {
                    (0, 0): self.get_mat_ss(),
                    (0, 1): self.get_mat_sp(),
                    (1, 1): self.get_mat_pp()
                }, mat_t.symm)
        else:
            mat_a = None

        mat_a = Matrix.bcast(mat_a, rank, comm)

        mat_b = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symm)

        # TODO: fix "'Matrix' object has no attribute 'get_submatrices'"
        # Tester.compare_matrices(mat_a, mat_b)
        assert mat_a.get_type() == mat_b.get_type()
        for key in [(0, 0), (0, 1), (1, 1)]:
            Tester.compare_submatrices(mat_a.get_submatrix(key),
                                       mat_b.get_submatrix(key))

    def test_mpi_reduce(self):

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nodes = comm.Get_size()

        if is_master(rank):
            mat_a = Matrix(
                {
                    (0, 0): self.get_mat_ss(),
                    (0, 1): self.get_mat_sp(),
                    (1, 1): self.get_mat_pp()
                }, mat_t.symm)
        elif rank == 1:
            mat_a = Matrix(
                {
                    (0, 0): self.get_mat_ss_r1(),
                    (0, 1): self.get_mat_sp_r1(),
                    (1, 1): self.get_mat_pp_r1()
                }, mat_t.symm)
        else:
            mat_a = Matrix(
                {
                    (0, 0): SubMatrix([0, 0, 3, 3]),
                    (0, 1): SubMatrix([0, 3, 3, 6]),
                    (1, 1): SubMatrix([3, 3, 6, 6])
                }, mat_t.symm)

        mat_a = Matrix.reduce(mat_a, rank, comm)

        if is_master(rank):
            if nodes == 1:
                mat_b = Matrix(
                    {
                        (0, 0): self.get_mat_ss(),
                        (0, 1): self.get_mat_sp(),
                        (1, 1): self.get_mat_pp()
                    }, mat_t.symm)
            else:
                mat_b = Matrix(
                    {
                        (0, 0): self.get_mat_ss_f(),
                        (0, 1): self.get_mat_sp_f(),
                        (1, 1): self.get_mat_pp_f()
                    }, mat_t.symm)

            # TODO: fix "'Matrix' object has no attribute 'get_submatrices'"
            # Tester.compare_matrices(mat_a, mat_b)
            assert mat_a.get_type() == mat_b.get_type()
            for key in [(0, 0), (0, 1), (1, 1)]:
                Tester.compare_submatrices(mat_a.get_submatrix(key),
                                           mat_b.get_submatrix(key))
