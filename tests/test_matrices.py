from mpi4py import MPI

from veloxchem.veloxchemlib import mat_t
from veloxchem.submatrix import SubMatrix
from veloxchem.matrix import Matrix
from veloxchem.matrices import Matrices


class TestMatrices:

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

    def get_mat_ps(self):

        return SubMatrix([
            1.5, 1.6, 1.7, 2.5, 2.6, 2.7, 3.5, 3.6, 3.7, 4.5, 4.6, 4.7, 5.5,
            5.6, 5.7, 6.5, 6.6, 6.7
        ], [3, 0, 6, 3])

    def get_symm_matrix(self):

        return Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symmetric)

    def get_gen_matrix(self):

        return Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 0): self.get_mat_ps(),
                (1, 1): self.get_mat_pp()
            }, mat_t.general)

    def test_add_with_string(self):

        mats_a = Matrices()
        mats_a.add(self.get_symm_matrix(), "0")
        mats_a.add(self.get_gen_matrix(), "2")
        mats_b = Matrices({
            "0": self.get_symm_matrix(),
            "2": self.get_gen_matrix()
        })
        assert mats_a == mats_b

    def test_add_with_integer(self):

        mats_a = Matrices()
        mats_a.add(self.get_symm_matrix(), 0)
        mats_a.add(self.get_gen_matrix(), 2)
        mats_b = Matrices({
            "0": self.get_symm_matrix(),
            "2": self.get_gen_matrix()
        })
        assert mats_a == mats_b

    def test_zero(self):

        mat_sym = self.get_symm_matrix()
        mat_gen = self.get_gen_matrix()
        mats_a = Matrices({"0": mat_sym, "2": mat_gen})
        mats_a.zero()
        mats_b = Matrices()
        mat_sym.zero()
        mat_gen.zero()
        mats_b.add(mat_sym, "0")
        mats_b.add(mat_gen, "2")
        assert mats_a == mats_b

    def test_scale(self):

        mat_sym = self.get_symm_matrix()
        mat_gen = self.get_gen_matrix()
        mats_a = Matrices({"0": mat_sym, "2": mat_gen})
        mats_a.scale(0.71)
        mats_b = Matrices()
        mat_sym.scale(0.71)
        mat_gen.scale(0.71)
        mats_b.add(mat_sym, "0")
        mats_b.add(mat_gen, "2")
        assert mats_a == mats_b

    def test_symmetrize(self):

        mat_sym = self.get_symm_matrix()
        mat_gen = self.get_gen_matrix()
        mats_a = Matrices({"0": mat_sym, "2": mat_gen})
        mats_a.symmetrize()
        mats_b = Matrices()
        mat_sym.symmetrize()
        mat_gen.symmetrize()
        mats_b.add(mat_sym, "0")
        mats_b.add(mat_gen, "2")
        assert mats_a == mats_b

    def test_get_keys(self):

        mats = Matrices({
            "0": self.get_symm_matrix(),
            "2": self.get_gen_matrix()
        })
        assert mats.keys() == ["0", "2"]

    def test_get_matrix_with_string(self):

        mat_sym = self.get_symm_matrix()
        mat_gen = self.get_gen_matrix()
        mats = Matrices({"0": mat_sym, "2": mat_gen})
        assert mats.matrix("0") == mat_sym
        assert mats.matrix("2") == mat_gen

    def test_get_matrix_with_integer(self):

        mat_sym = self.get_symm_matrix()
        mat_gen = self.get_gen_matrix()
        mats = Matrices({"0": mat_sym, "2": mat_gen})
        assert mats.matrix(0) == mat_sym
        assert mats.matrix(2) == mat_gen

    def test_mpi_bcast(self):

        comm = MPI.COMM_WORLD

        mats_a = None
        if comm.Get_rank() == 0:
            mats_a = Matrices({
                "0": self.get_symm_matrix(),
                "2": self.get_gen_matrix()
            })
        mats_a = Matrices.bcast(mats_a, comm, 0)
        mats_b = Matrices({
            "0": self.get_symm_matrix(),
            "2": self.get_gen_matrix()
        })
        assert mats_a == mats_b

    def test_mpi_reduce(self):

        comm = MPI.COMM_WORLD

        mats_a = Matrices({
            "0": self.get_symm_matrix(),
            "2": self.get_gen_matrix()
        })
        if comm.Get_rank() == 0:
            mats_a.scale(0.82)
        mats_r = Matrices.reduce(mats_a, comm, 0)
        if comm.Get_rank() == 0:
            mats_b = Matrices({
                "0": self.get_symm_matrix(),
                "2": self.get_gen_matrix()
            })
            mats_b.scale(0.82 + comm.Get_size() - 1)
            assert mats_r == mats_b
