from veloxchem.matrix import Matrix
from veloxchem.matrices import Matrices
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.veloxchemlib import make_matrix
from veloxchem.veloxchemlib import make_matrices
from veloxchem.veloxchemlib import mat_t


class TestMatricesFunc:

    def get_data(self):

        h2ostr = """O   0.000   0.000  -1.000
                    H   0.000   1.400  -2.100
                    H   0.000  -1.400  -2.100"""

        mol = Molecule.read_str(h2ostr, 'au')
        bas = MolecularBasis.read(mol, 'DEF2-SVP', ostream=None)

        return bas

    def test_make_matrices_asym(self):

        bas_svp = self.get_data()
        mats_a = make_matrices(2, bas_svp, mat_t.antisymmetric)
        mats_a.zero()
        matrix = make_matrix(bas_svp, mat_t.antisymmetric)
        mats_b = Matrices()
        mats_b.add(matrix, "XX")
        mats_b.add(matrix, "XY")
        mats_b.add(matrix, "XZ")
        mats_b.add(matrix, "YY")
        mats_b.add(matrix, "YZ")
        mats_b.add(matrix, "ZZ")
        mats_b.zero()
        assert mats_a == mats_b

    def test_make_matrices_gen(self):

        bas_svp = self.get_data()
        mats_a = make_matrices(
            2,
            bas_svp,
            bas_svp,
        )
        mats_a.zero()
        matrix = make_matrix(bas_svp, mat_t.general)
        mats_b = Matrices()
        mats_b.add(matrix, "XX")
        mats_b.add(matrix, "XY")
        mats_b.add(matrix, "XZ")
        mats_b.add(matrix, "YY")
        mats_b.add(matrix, "YZ")
        mats_b.add(matrix, "ZZ")
        mats_b.zero()
        assert mats_a == mats_b

    def test_make_matrices_tdp(self):

        bas_svp = self.get_data()
        mats_a = make_matrices([2, 1], bas_svp, mat_t.antisymmetric)
        mats_a.zero()
        matrix = make_matrix(bas_svp, mat_t.antisymmetric)
        mats_b = Matrices()
        mats_b.add(matrix, "XX_X")
        mats_b.add(matrix, "XX_Y")
        mats_b.add(matrix, "XX_Z")
        mats_b.add(matrix, "XY_X")
        mats_b.add(matrix, "XY_Y")
        mats_b.add(matrix, "XY_Z")
        mats_b.add(matrix, "XZ_X")
        mats_b.add(matrix, "XZ_Y")
        mats_b.add(matrix, "XZ_Z")
        mats_b.add(matrix, "YY_X")
        mats_b.add(matrix, "YY_Y")
        mats_b.add(matrix, "YY_Z")
        mats_b.add(matrix, "YZ_X")
        mats_b.add(matrix, "YZ_Y")
        mats_b.add(matrix, "YZ_Z")
        mats_b.add(matrix, "ZZ_X")
        mats_b.add(matrix, "ZZ_Y")
        mats_b.add(matrix, "ZZ_Z")
        mats_b.zero()
        assert mats_a == mats_b

    def test_make_matrices_tpd(self):

        bas_svp = self.get_data()
        mats_a = make_matrices([1, 2], bas_svp, mat_t.antisymmetric)
        mats_a.zero()
        matrix = make_matrix(bas_svp, mat_t.antisymmetric)
        mats_b = Matrices()
        mats_b.add(matrix, "X_XX")
        mats_b.add(matrix, "X_XY")
        mats_b.add(matrix, "X_XZ")
        mats_b.add(matrix, "X_YY")
        mats_b.add(matrix, "X_YZ")
        mats_b.add(matrix, "X_ZZ")
        mats_b.add(matrix, "Y_XX")
        mats_b.add(matrix, "Y_XY")
        mats_b.add(matrix, "Y_XZ")
        mats_b.add(matrix, "Y_YY")
        mats_b.add(matrix, "Y_YZ")
        mats_b.add(matrix, "Y_ZZ")
        mats_b.add(matrix, "Z_XX")
        mats_b.add(matrix, "Z_XY")
        mats_b.add(matrix, "Z_XZ")
        mats_b.add(matrix, "Z_YY")
        mats_b.add(matrix, "Z_YZ")
        mats_b.add(matrix, "Z_ZZ")

        mats_b.zero()
        assert mats_a == mats_b

    def test_make_matrices_tppp(self):

        bas_svp = self.get_data()
        mats_a = make_matrices(
            [1, 1, 1],
            bas_svp,
            bas_svp,
        )
        mats_a.zero()
        matrix = make_matrix(bas_svp, mat_t.general)
        mats_b = Matrices()
        mats_b.add(matrix, "X_X_X")
        mats_b.add(matrix, "X_X_Y")
        mats_b.add(matrix, "X_X_Z")
        mats_b.add(matrix, "X_Y_X")
        mats_b.add(matrix, "X_Y_Y")
        mats_b.add(matrix, "X_Y_Z")
        mats_b.add(matrix, "X_Z_X")
        mats_b.add(matrix, "X_Z_Y")
        mats_b.add(matrix, "X_Z_Z")
        mats_b.add(matrix, "Y_X_X")
        mats_b.add(matrix, "Y_X_Y")
        mats_b.add(matrix, "Y_X_Z")
        mats_b.add(matrix, "Y_Y_X")
        mats_b.add(matrix, "Y_Y_Y")
        mats_b.add(matrix, "Y_Y_Z")
        mats_b.add(matrix, "Y_Z_X")
        mats_b.add(matrix, "Y_Z_Y")
        mats_b.add(matrix, "Y_Z_Z")
        mats_b.add(matrix, "Z_X_X")
        mats_b.add(matrix, "Z_X_Y")
        mats_b.add(matrix, "Z_X_Z")
        mats_b.add(matrix, "Z_Y_X")
        mats_b.add(matrix, "Z_Y_Y")
        mats_b.add(matrix, "Z_Y_Z")
        mats_b.add(matrix, "Z_Z_X")
        mats_b.add(matrix, "Z_Z_Y")
        mats_b.add(matrix, "Z_Z_Z")
        mats_b.zero()
        assert mats_a == mats_b

    def test_make_matrices_tpppp(self):

        bas_svp = self.get_data()
        mats_a = make_matrices(
            [1, 1, 1, 1],
            bas_svp,
            bas_svp,
        )
        mats_a.zero()
        matrix = make_matrix(bas_svp, mat_t.general)
        mats_b = Matrices()
        mats_b.add(matrix, "X_X_X_X")
        mats_b.add(matrix, "X_X_X_Y")
        mats_b.add(matrix, "X_X_X_Z")
        mats_b.add(matrix, "X_X_Y_X")
        mats_b.add(matrix, "X_X_Y_Y")
        mats_b.add(matrix, "X_X_Y_Z")
        mats_b.add(matrix, "X_X_Z_X")
        mats_b.add(matrix, "X_X_Z_Y")
        mats_b.add(matrix, "X_X_Z_Z")
        mats_b.add(matrix, "X_Y_X_X")
        mats_b.add(matrix, "X_Y_X_Y")
        mats_b.add(matrix, "X_Y_X_Z")
        mats_b.add(matrix, "X_Y_Y_X")
        mats_b.add(matrix, "X_Y_Y_Y")
        mats_b.add(matrix, "X_Y_Y_Z")
        mats_b.add(matrix, "X_Y_Z_X")
        mats_b.add(matrix, "X_Y_Z_Y")
        mats_b.add(matrix, "X_Y_Z_Z")
        mats_b.add(matrix, "X_Z_X_X")
        mats_b.add(matrix, "X_Z_X_Y")
        mats_b.add(matrix, "X_Z_X_Z")
        mats_b.add(matrix, "X_Z_Y_X")
        mats_b.add(matrix, "X_Z_Y_Y")
        mats_b.add(matrix, "X_Z_Y_Z")
        mats_b.add(matrix, "X_Z_Z_X")
        mats_b.add(matrix, "X_Z_Z_Y")
        mats_b.add(matrix, "X_Z_Z_Z")
        mats_b.add(matrix, "Y_X_X_X")
        mats_b.add(matrix, "Y_X_X_Y")
        mats_b.add(matrix, "Y_X_X_Z")
        mats_b.add(matrix, "Y_X_Y_X")
        mats_b.add(matrix, "Y_X_Y_Y")
        mats_b.add(matrix, "Y_X_Y_Z")
        mats_b.add(matrix, "Y_X_Z_X")
        mats_b.add(matrix, "Y_X_Z_Y")
        mats_b.add(matrix, "Y_X_Z_Z")
        mats_b.add(matrix, "Y_Y_X_X")
        mats_b.add(matrix, "Y_Y_X_Y")
        mats_b.add(matrix, "Y_Y_X_Z")
        mats_b.add(matrix, "Y_Y_Y_X")
        mats_b.add(matrix, "Y_Y_Y_Y")
        mats_b.add(matrix, "Y_Y_Y_Z")
        mats_b.add(matrix, "Y_Y_Z_X")
        mats_b.add(matrix, "Y_Y_Z_Y")
        mats_b.add(matrix, "Y_Y_Z_Z")
        mats_b.add(matrix, "Y_Z_X_X")
        mats_b.add(matrix, "Y_Z_X_Y")
        mats_b.add(matrix, "Y_Z_X_Z")
        mats_b.add(matrix, "Y_Z_Y_X")
        mats_b.add(matrix, "Y_Z_Y_Y")
        mats_b.add(matrix, "Y_Z_Y_Z")
        mats_b.add(matrix, "Y_Z_Z_X")
        mats_b.add(matrix, "Y_Z_Z_Y")
        mats_b.add(matrix, "Y_Z_Z_Z")
        mats_b.add(matrix, "Z_X_X_X")
        mats_b.add(matrix, "Z_X_X_Y")
        mats_b.add(matrix, "Z_X_X_Z")
        mats_b.add(matrix, "Z_X_Y_X")
        mats_b.add(matrix, "Z_X_Y_Y")
        mats_b.add(matrix, "Z_X_Y_Z")
        mats_b.add(matrix, "Z_X_Z_X")
        mats_b.add(matrix, "Z_X_Z_Y")
        mats_b.add(matrix, "Z_X_Z_Z")
        mats_b.add(matrix, "Z_Y_X_X")
        mats_b.add(matrix, "Z_Y_X_Y")
        mats_b.add(matrix, "Z_Y_X_Z")
        mats_b.add(matrix, "Z_Y_Y_X")
        mats_b.add(matrix, "Z_Y_Y_Y")
        mats_b.add(matrix, "Z_Y_Y_Z")
        mats_b.add(matrix, "Z_Y_Z_X")
        mats_b.add(matrix, "Z_Y_Z_Y")
        mats_b.add(matrix, "Z_Y_Z_Z")
        mats_b.add(matrix, "Z_Z_X_X")
        mats_b.add(matrix, "Z_Z_X_Y")
        mats_b.add(matrix, "Z_Z_X_Z")
        mats_b.add(matrix, "Z_Z_Y_X")
        mats_b.add(matrix, "Z_Z_Y_Y")
        mats_b.add(matrix, "Z_Z_Y_Z")
        mats_b.add(matrix, "Z_Z_Z_X")
        mats_b.add(matrix, "Z_Z_Z_Y")
        mats_b.add(matrix, "Z_Z_Z_Z")
        mats_b.zero()
        assert mats_a == mats_b
