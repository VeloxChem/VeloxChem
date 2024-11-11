from veloxchem.submatrix import SubMatrix
from veloxchem.matrix import Matrix
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.veloxchemlib import make_matrix
from veloxchem.veloxchemlib import mat_t


class TestMatrixFunc:

    def get_data(self):

        h2ostr = """O   0.000   0.000  -1.000
                    H   0.000   1.400  -2.100
                    H   0.000  -1.400  -2.100"""

        mol = Molecule.read_str(h2ostr, 'au')
        bas = MolecularBasis.read(mol, 'DEF2-SVP', ostream=None)

        return bas

    def test_make_matrix(self):

        bas_svp = self.get_data()

        mat_a = make_matrix(bas_svp, mat_t.antisymmetric)
        mat_a.zero()
        mat_b = Matrix()
        mat_b.set_type(mat_t.antisymmetric)
        mat_b.add(SubMatrix([0, 0, 7, 7]), (0, 0))
        mat_b.add(SubMatrix([0, 7, 7, 12]), (0, 1))
        mat_b.add(SubMatrix([0, 19, 7, 5]), (0, 2))
        mat_b.add(SubMatrix([7, 7, 12, 12]), (1, 1))
        mat_b.add(SubMatrix([7, 19, 12, 5]), (1, 2))
        mat_b.add(SubMatrix([19, 19, 5, 5]), (2, 2))
        mat_b.zero()
        assert mat_a == mat_b
        mat_a = make_matrix(bas_svp, mat_t.symmetric)
        mat_a.zero()
        mat_b.set_type(mat_t.symmetric)
        assert mat_a == mat_b

    def test_make_matrix_gen(self):

        bas_svp = self.get_data()
        mat_a = make_matrix(bas_svp, mat_t.general)
        mat_a.zero()
        mat_b = Matrix()
        mat_b.set_type(mat_t.general)
        mat_b.add(SubMatrix([0, 0, 7, 7]), (0, 0))
        mat_b.add(SubMatrix([0, 7, 7, 12]), (0, 1))
        mat_b.add(SubMatrix([0, 19, 7, 5]), (0, 2))
        mat_b.add(SubMatrix([7, 0, 12, 7]), (1, 0))
        mat_b.add(SubMatrix([7, 7, 12, 12]), (1, 1))
        mat_b.add(SubMatrix([7, 19, 12, 5]), (1, 2))
        mat_b.add(SubMatrix([19, 0, 5, 7]), (2, 0))
        mat_b.add(SubMatrix([19, 7, 5, 12]), (2, 1))
        mat_b.add(SubMatrix([19, 19, 5, 5]), (2, 2))
        mat_b.zero()
        assert mat_a == mat_b
        mat_a = make_matrix(bas_svp, bas_svp)
        mat_a.zero()
        assert mat_a == mat_b
