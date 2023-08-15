from veloxchem.matrix import Matrix
from veloxchem.submatrix import SubMatrix
from veloxchem.veloxchemlib import FockMatrix
from veloxchem.veloxchemlib import FockMatrices
from veloxchem.veloxchemlib import fock_t
from veloxchem.veloxchemlib import mat_t
from tester import Tester


class TestFockMatrices:

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

    def get_symm_mat(self):

        return Matrix(
            {
                (0, 0): self.get_mat_ss(),
                (0, 1): self.get_mat_sp(),
                (1, 1): self.get_mat_pp()
            }, mat_t.symm)

    def get_fock_matrix(self, scale, ftype):

        return FockMatrix(self.get_symm_mat(), scale, ftype)

    def test_add(self):

        fmats_a = FockMatrices([
            self.get_fock_matrix(0.2, fock_t.restjk),
            self.get_fock_matrix(0.8, fock_t.restjkx)
        ])

        fmats_b = FockMatrices()
        fmats_b.add(self.get_fock_matrix(0.2, fock_t.restjk))
        fmats_b.add(self.get_fock_matrix(0.8, fock_t.restjkx))

        Tester.compare_list_of_fockmatrices(fmats_a, fmats_b)

    def test_get_matrix(self):

        fmats = FockMatrices([
            self.get_fock_matrix(0.2, fock_t.restjk),
            self.get_fock_matrix(0.8, fock_t.restjkx)
        ])

        Tester.compare_fockmatrices(fmats.get_matrix(0),
                                    self.get_fock_matrix(0.2, fock_t.restjk))
        Tester.compare_fockmatrices(fmats.get_matrix(1),
                                    self.get_fock_matrix(0.8, fock_t.restjkx))

    def test_number_of_matrices(self):

        fmats = FockMatrices([
            self.get_fock_matrix(0.2, fock_t.restjk),
            self.get_fock_matrix(0.8, fock_t.restjkx)
        ])

        assert fmats.number_of_matrices() == 2
