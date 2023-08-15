import math as mt

from veloxchem.matrix import Matrix
from veloxchem.submatrix import SubMatrix
from veloxchem.veloxchemlib import FockMatrix
from veloxchem.veloxchemlib import fock_t
from veloxchem.veloxchemlib import mat_t
from tester import Tester


class TestFockMatrix:

    def get_mat_ss(self):

        return SubMatrix([1.0, 2.0, 3.0, 1.1, 2.1, 3.1, 1.2, 2.2, 3.2],
                         [0, 0, 3, 3])

    def get_mat_sp(self):

        return SubMatrix([
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 1.2,
            2.2, 3.2, 4.2, 5.2, 6.2
        ], [0, 3, 3, 6])

    def get_mat_pp(self):

        return SubMatrix([
            2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 2.2,
            3.2, 4.2, 5.2, 6.2, 7.2, 2.3, 3.3, 4.3, 5.3, 6.3, 7.3, 2.4, 3.4,
            4.4, 5.4, 6.4, 7.4, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5
        ], [3, 3, 6, 6])

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

    def get_symm_mat(self):

        return Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symm)

    def test_set_and_get_fock_type(self):

        fmat = FockMatrix(self.get_symm_mat(), 1.0, fock_t.restk)

        assert fmat.get_fock_type() == fock_t.restk

        fmat.set_fock_type(fock_t.restjkx)

        assert fmat.get_fock_type() == fock_t.restjkx

    def test_set_and_get_exchange_scale(self):

        tol = 1.0e-12

        fmat = FockMatrix(self.get_symm_mat(), 1.0, fock_t.restk)

        assert mt.isclose(fmat.get_exchange_scale(),
                          1.0,
                          rel_tol=tol,
                          abs_tol=tol)

        fmat.set_exchange_scale(0.5)

        assert mt.isclose(fmat.get_exchange_scale(),
                          0.5,
                          rel_tol=tol,
                          abs_tol=tol)

    def test_zero(self):

        fmat = FockMatrix(self.get_symm_mat(), 1.0, fock_t.restk)

        smat_ss = self.get_mat_ss()
        smat_sp = self.get_mat_sp()
        smat_pp = self.get_mat_pp()

        Tester.compare_submatrices(smat_ss, fmat.get_submatrix((0, 0)))
        Tester.compare_submatrices(smat_sp, fmat.get_submatrix((0, 1)))
        Tester.compare_submatrices(smat_pp, fmat.get_submatrix((1, 1)))

        smat_ss.zero()
        smat_sp.zero()
        smat_pp.zero()

        fmat.zero()

        Tester.compare_submatrices(smat_ss, fmat.get_submatrix((0, 0)))
        Tester.compare_submatrices(smat_sp, fmat.get_submatrix((0, 1)))
        Tester.compare_submatrices(smat_pp, fmat.get_submatrix((1, 1)))

    def test_get_angular_pairs(self):

        fmat = FockMatrix(self.get_symm_mat(), 1.0, fock_t.restk)

        assert fmat.get_angular_pairs() == [(0, 0), (0, 1), (1, 1)]

    def test_get_matrix(self):

        fmat = FockMatrix(self.get_symm_mat(), 1.0, fock_t.restk)

        Tester.compare_matrices(fmat.get_matrix(), self.get_symm_mat())

    def test_get_submatrix(self):

        fmat = FockMatrix(self.get_symm_mat(), 1.0, fock_t.restk)

        Tester.compare_submatrices(fmat.get_submatrix((0, 0)),
                                   self.get_mat_ss())
        Tester.compare_submatrices(fmat.get_submatrix((0, 1)),
                                   self.get_mat_sp())
        Tester.compare_submatrices(fmat.get_submatrix((1, 1)),
                                   self.get_mat_pp())

    def test_is_angular_order(self):

        fmat = FockMatrix(self.get_symm_mat(), 1.0, fock_t.restk)

        assert fmat.is_angular_order((0, 0))
        assert fmat.is_angular_order((0, 1))
        assert fmat.is_angular_order((1, 1))

        assert fmat.is_angular_order((1, 0)) is False
        assert fmat.is_angular_order((1, 2)) is False
        assert fmat.is_angular_order((2, 1)) is False
        assert fmat.is_angular_order((2, 2)) is False

    def test_is_number_of_rows(self):

        fmat = FockMatrix(self.get_symm_mat(), 1.0, fock_t.restk)

        assert fmat.number_of_rows() == 9

    def test_is_number_of_columns(self):

        fmat = FockMatrix(self.get_symm_mat(), 1.0, fock_t.restk)

        assert fmat.number_of_columns() == 9

    def test_get_full_matrix_symm(self):

        fmat = FockMatrix(self.get_symm_mat(), 1.0, fock_t.restk)

        Tester.compare_submatrices(fmat.get_full_matrix(),
                                   self.get_mat_full_symm())
