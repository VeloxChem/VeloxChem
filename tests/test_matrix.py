import pickle
import numpy as np
import math as mt

from mpi4py import MPI

from veloxchem.veloxchemlib import mat_t
from veloxchem.submatrix import SubMatrix
from veloxchem.matrix import Matrix


class TestMatrix:

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
            -2.2, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, -3.0, -3.1, -3.2, 2.2, 3.2,
            4.2, 5.2, 6.2, 7.2, -4.0, -4.1, -4.2, 2.3, 3.3, 4.3, 5.3, 6.3, 7.3,
            -5.0, -5.1, -5.2, 2.4, 3.4, 4.4, 5.4, 6.4, 7.4, -6.0, -6.1, -6.2,
            2.5, 3.5, 4.5, 5.5, 6.5, 7.5
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

    def test_pickle(self):

        mat_a = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symmetric)
        bobj = pickle.dumps(mat_a)
        mat_b = pickle.loads(bobj)

        #print(mat_a.full_matrix().to_numpy())
        #print(mat_b.full_matrix().to_numpy())
        assert mat_a == mat_b

    def test_sum_op(self):

        mat_a = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symmetric)
        mat_a.scale(0.9)
        mat_b = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symmetric)
        mat_c = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symmetric)
        mat_c.scale(1.9)
        assert mat_c == mat_a + mat_b

    def test_add(self):

        mat_a = Matrix()
        mat_a.set_type(mat_t.symmetric)
        mat_a.add(self.get_mat_ss(), (0, 0))
        mat_a.add(self.get_mat_sp(), (0, 1))
        mat_a.add(self.get_mat_pp(), (1, 1))
        mat_b = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symmetric)
        assert mat_a == mat_b

    def test_set_and_get_type(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (1, 0): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symmetric)
        assert mat.get_type() == mat_t.symmetric
        mat.set_type(mat_t.antisymmetric)
        assert mat.get_type() == mat_t.antisymmetric

    def test_zero(self):

        smat_ss = self.get_mat_ss()
        smat_sp = self.get_mat_sp()
        smat_pp = self.get_mat_pp()

        mat_a = Matrix()
        mat_a.set_type(mat_t.symmetric)
        mat_a.add(smat_ss, (0, 0))
        mat_a.add(smat_sp, (0, 1))
        mat_a.add(smat_pp, (1, 1))
        mat_a.zero()
        smat_ss.zero()
        smat_sp.zero()
        smat_pp.zero()
        mat_b = Matrix()
        mat_b.set_type(mat_t.symmetric)
        mat_b.add(smat_ss, (0, 0))
        mat_b.add(smat_sp, (0, 1))
        mat_b.add(smat_pp, (1, 1))
        assert mat_a == mat_b

    def test_scale(self):

        smat_ss = self.get_mat_ss()
        smat_sp = self.get_mat_sp()
        smat_pp = self.get_mat_pp()

        mat_a = Matrix()
        mat_a.set_type(mat_t.symmetric)
        mat_a.add(smat_ss, (0, 0))
        mat_a.add(smat_sp, (0, 1))
        mat_a.add(smat_pp, (1, 1))
        mat_a.scale(0.71)
        smat_ss.scale(0.71)
        smat_sp.scale(0.71)
        smat_pp.scale(0.71)
        mat_b = Matrix()
        mat_b.set_type(mat_t.symmetric)
        mat_b.add(smat_ss, (0, 0))
        mat_b.add(smat_sp, (0, 1))
        mat_b.add(smat_pp, (1, 1))
        assert mat_a == mat_b

    def test_symmetrize(self):

        smat_ss = self.get_mat_ss()
        smat_sp = self.get_mat_sp()
        smat_pp = self.get_mat_pp()

        mat_a = Matrix()
        mat_a.set_type(mat_t.symmetric)
        mat_a.add(smat_ss, (0, 0))
        mat_a.add(smat_sp, (0, 1))
        mat_a.add(smat_pp, (1, 1))
        mat_a.symmetrize()
        smat_ss.symmetrize()
        smat_pp.symmetrize()
        mat_b = Matrix()
        mat_b.set_type(mat_t.symmetric)
        mat_b.add(smat_ss, (0, 0))
        mat_b.add(smat_sp, (0, 1))
        mat_b.add(smat_pp, (1, 1))
        assert mat_a == mat_b

    def test_angular_pairs(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symmetric)
        assert mat.angular_pairs() == [(0, 0), (0, 1), (1, 1)]

    def test_submatrix(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symmetric)
        assert mat.submatrix((0, 0)) == self.get_mat_ss()
        assert mat.submatrix((0, 1)) == self.get_mat_sp()
        assert mat.submatrix((1, 1)) == self.get_mat_pp()

    def test_is_angular_order(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symmetric)
        assert mat.is_angular_order((0, 0))
        assert mat.is_angular_order((0, 1))
        assert mat.is_angular_order((1, 1))
        assert mat.is_angular_order((1, 0)) is False
        assert mat.is_angular_order((1, 2)) is False
        assert mat.is_angular_order((2, 1)) is False
        assert mat.is_angular_order((2, 2)) is False

    def test_number_of_rows(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symmetric)
        assert mat.number_of_rows() == 9
        mat = Matrix({
            (0, 0): self.get_mat_ss(),
            (0, 1): self.get_mat_sp()
        }, mat_t.general)
        assert mat.number_of_rows() == 3

    def test_number_of_columns(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symmetric)
        assert mat.number_of_columns() == 9
        mat = Matrix({
            (0, 0): self.get_mat_ss(),
            (0, 1): self.get_mat_sp()
        }, mat_t.general)
        assert mat.number_of_columns() == 9

    def test_full_matrix_symm(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symmetric)
        assert mat.full_matrix() == self.get_mat_full_symm()

    def test_full_matrix_antisymm(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.antisymmetric)
        assert mat.full_matrix() == self.get_mat_full_antisymm()

    def test_full_matrix_gen(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 0): self.get_mat_ps(),
                (1, 1): self.get_mat_pp()
            }, mat_t.general)
        assert mat.full_matrix() == self.get_mat_full_gen()

    def test_set_values_symm(self):

        tol = 1.0e-12

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symmetric)
        mvals = np.random.rand(9, 9)
        mat.set_values(mvals)

        ssmat = mat.submatrix((0, 0))
        spmat = mat.submatrix((0, 1))
        ppmat = mat.submatrix((1, 1))
        assert np.allclose(ssmat.to_numpy(), mvals[:3, :3], tol, tol, False)
        assert np.allclose(spmat.to_numpy(), mvals[:3, 3:9], tol, tol, False)
        assert np.allclose(ppmat.to_numpy(), mvals[3:9, 3:9], tol, tol, False)

    def test_set_values_full(self):

        tol = 1.0e-12

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 0): self.get_mat_ps(),
                (1, 1): self.get_mat_pp()
            }, mat_t.general)
        mvals = np.random.rand(9, 9)
        mat.set_values(mvals)

        ssmat = mat.submatrix((0, 0))
        spmat = mat.submatrix((0, 1))
        psmat = mat.submatrix((1, 0))
        ppmat = mat.submatrix((1, 1))
        assert np.allclose(ssmat.to_numpy(), mvals[:3, :3], tol, tol, False)
        assert np.allclose(spmat.to_numpy(), mvals[:3, 3:9], tol, tol, False)
        assert np.allclose(psmat.to_numpy(), mvals[3:9, :3], tol, tol, False)
        assert np.allclose(ppmat.to_numpy(), mvals[3:9, 3:9], tol, tol, False)

        fmat = mat.full_matrix()
        assert np.allclose(fmat.to_numpy(), mvals, tol, tol, False)

    def test_mpi_bcast(self):

        comm = MPI.COMM_WORLD

        mat_a = None
        if comm.Get_rank() == 0:
            mat_a = Matrix(
                {
                    (0, 0): self.get_mat_ss(),
                    (0, 1): self.get_mat_sp(),
                    (1, 0): self.get_mat_ps(),
                    (1, 1): self.get_mat_pp()
                }, mat_t.general)
        mat_a = comm.bcast(mat_a, 0)
        mat_b = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 0): self.get_mat_ps(),
                (1, 1): self.get_mat_pp()
            }, mat_t.general)
        assert mat_a == mat_b

    def test_flat_values_symm(self):

        mat = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symmetric)

        mat_vals = mat.flat_values()

        ref_vals = [
            1.0, 4.0, 6.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 2.1, 6.2, 2.2, 4.2,
            6.2, 8.2, 10.2, 12.2, 3.2, 2.4, 4.4, 6.4, 8.4, 10.4, 12.4, 2.0,
            6.0, 8.0, 10.0, 12.0, 14.0, 3.1, 8.2, 10.2, 12.2, 14.2, 4.2, 10.4,
            12.4, 14.4, 5.3, 12.6, 14.6, 6.4, 14.8, 7.5
        ]

        for i in range(45):
            assert mt.isclose(ref_vals[i],
                              mat_vals[i],
                              rel_tol=1.0e-12,
                              abs_tol=1.0e-12)

    def test_mpi_reduce(self):

        comm = MPI.COMM_WORLD

        mat_a = Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 0): self.get_mat_ps(),
                (1, 1): self.get_mat_pp()
            }, mat_t.general)
        if comm.Get_rank() == 0:
            mat_a.scale(0.82)
        mat_r = Matrix.reduce(mat_a, comm, 0)
        if comm.Get_rank() == 0:
            mat_b = Matrix(
                {
                    (0, 0): self.get_mat_ss(),
                    (0, 1): self.get_mat_sp(),
                    (1, 0): self.get_mat_ps(),
                    (1, 1): self.get_mat_pp()
                }, mat_t.general)
            mat_b.scale(0.82 + comm.Get_size() - 1)
            assert mat_r == mat_b
