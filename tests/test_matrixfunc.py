from veloxchem.veloxchemlib import SubMatrix
from veloxchem.veloxchemlib import make_matrix
from veloxchem.veloxchemlib import mat_t
from veloxchem.veloxchemlib import MolecularBasis
from veloxchem.veloxchemlib import Molecule
from tester import Tester


class TestMatrixFunc:

    def get_data(self):

        h2ostr = """O   0.000   0.000  -1.000
                    H   0.000   1.400  -2.100
                    H   0.000  -1.400  -2.100"""

        mol = Molecule.read_str(h2ostr, 'au')

        bas = MolecularBasis.read(mol, 'DEF2-SVP', 'basis')

        return bas

    def test_make_matrix(self):

        bas_svp = self.get_data()

        mat = make_matrix(bas_svp, mat_t.antisymm)
        Tester.compare_submatrices(mat.get_submatrix((0, 0)),
                                   SubMatrix([0, 0, 7, 7]))
        Tester.compare_submatrices(mat.get_submatrix((0, 1)),
                                   SubMatrix([0, 7, 7, 12]))
        Tester.compare_submatrices(mat.get_submatrix((0, 2)),
                                   SubMatrix([0, 19, 7, 5]))
        Tester.compare_submatrices(mat.get_submatrix((1, 1)),
                                   SubMatrix([7, 7, 12, 12]))
        Tester.compare_submatrices(mat.get_submatrix((1, 2)),
                                   SubMatrix([7, 19, 12, 5]))
        Tester.compare_submatrices(mat.get_submatrix((2, 2)),
                                   SubMatrix([19, 19, 5, 5]))
        assert mat_t.antisymm == mat.get_type()

        mat = make_matrix(bas_svp, mat_t.symm)
        Tester.compare_submatrices(mat.get_submatrix((0, 0)),
                                   SubMatrix([0, 0, 7, 7]))
        Tester.compare_submatrices(mat.get_submatrix((0, 1)),
                                   SubMatrix([0, 7, 7, 12]))
        Tester.compare_submatrices(mat.get_submatrix((0, 2)),
                                   SubMatrix([0, 19, 7, 5]))
        Tester.compare_submatrices(mat.get_submatrix((1, 1)),
                                   SubMatrix([7, 7, 12, 12]))
        Tester.compare_submatrices(mat.get_submatrix((1, 2)),
                                   SubMatrix([7, 19, 12, 5]))
        Tester.compare_submatrices(mat.get_submatrix((2, 2)),
                                   SubMatrix([19, 19, 5, 5]))
        assert mat_t.symm == mat.get_type()

    def test_make_matrix_gen(self):

        bas_svp = self.get_data()

        mat = make_matrix(bas_svp, mat_t.gen)
        Tester.compare_submatrices(mat.get_submatrix((0, 0)),
                                   SubMatrix([0, 0, 7, 7]))
        Tester.compare_submatrices(mat.get_submatrix((0, 1)),
                                   SubMatrix([0, 7, 7, 12]))
        Tester.compare_submatrices(mat.get_submatrix((0, 2)),
                                   SubMatrix([0, 19, 7, 5]))
        Tester.compare_submatrices(mat.get_submatrix((1, 0)),
                                   SubMatrix([7, 0, 12, 7]))
        Tester.compare_submatrices(mat.get_submatrix((1, 1)),
                                   SubMatrix([7, 7, 12, 12]))
        Tester.compare_submatrices(mat.get_submatrix((1, 2)),
                                   SubMatrix([7, 19, 12, 5]))
        Tester.compare_submatrices(mat.get_submatrix((2, 0)),
                                   SubMatrix([19, 0, 5, 7]))
        Tester.compare_submatrices(mat.get_submatrix((2, 1)),
                                   SubMatrix([19, 7, 5, 12]))
        Tester.compare_submatrices(mat.get_submatrix((2, 2)),
                                   SubMatrix([19, 19, 5, 5]))
        assert mat_t.gen == mat.get_type()

        mat = make_matrix(bas_svp, bas_svp)
        Tester.compare_submatrices(mat.get_submatrix((0, 0)),
                                   SubMatrix([0, 0, 7, 7]))
        Tester.compare_submatrices(mat.get_submatrix((0, 1)),
                                   SubMatrix([0, 7, 7, 12]))
        Tester.compare_submatrices(mat.get_submatrix((0, 2)),
                                   SubMatrix([0, 19, 7, 5]))
        Tester.compare_submatrices(mat.get_submatrix((1, 0)),
                                   SubMatrix([7, 0, 12, 7]))
        Tester.compare_submatrices(mat.get_submatrix((1, 1)),
                                   SubMatrix([7, 7, 12, 12]))
        Tester.compare_submatrices(mat.get_submatrix((1, 2)),
                                   SubMatrix([7, 19, 12, 5]))
        Tester.compare_submatrices(mat.get_submatrix((2, 0)),
                                   SubMatrix([19, 0, 5, 7]))
        Tester.compare_submatrices(mat.get_submatrix((2, 1)),
                                   SubMatrix([19, 7, 5, 12]))
        Tester.compare_submatrices(mat.get_submatrix((2, 2)),
                                   SubMatrix([19, 19, 5, 5]))
        assert mat_t.gen == mat.get_type()
