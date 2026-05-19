import numpy as np

from veloxchem.tabulalib import TabulaDenseMatrix, TabulaSymmetry


class TestTabulaDenseMatrix:
    """Tests for tabula::DenseMatrix — Tabula's own dense integral-result
    storage, exposed from the tabula module of veloxchemlib."""

    def test_default_construction(self):

        mat = TabulaDenseMatrix()
        assert mat.rows() == 0
        assert mat.columns() == 0
        assert mat.symmetry() == TabulaSymmetry.general

    def test_shape_and_symmetry(self):

        gen = TabulaDenseMatrix(2, 4)
        assert gen.rows() == 2
        assert gen.columns() == 4
        assert gen.symmetry() == TabulaSymmetry.general

        sym = TabulaDenseMatrix(3, 3, TabulaSymmetry.symmetric)
        assert sym.symmetry() == TabulaSymmetry.symmetric

        asym = TabulaDenseMatrix(3, 3, TabulaSymmetry.antisymmetric)
        assert asym.symmetry() == TabulaSymmetry.antisymmetric

    def test_element_access(self):

        mat = TabulaDenseMatrix(2, 3)
        mat[0, 0] = 1.5
        mat[1, 2] = -4.0
        assert mat[0, 0] == 1.5
        assert mat[1, 2] == -4.0
        # an untouched element stays zero-initialized
        assert mat[0, 2] == 0.0

    def test_zero(self):

        mat = TabulaDenseMatrix(2, 2)
        mat[0, 1] = 7.0
        mat[1, 0] = 3.0
        mat.zero()
        assert mat[0, 1] == 0.0 and mat[1, 0] == 0.0

    def test_scale(self):

        mat = TabulaDenseMatrix(2, 2)
        mat[0, 1] = 2.0
        mat[1, 1] = -5.0
        mat.scale(3.0)
        assert mat[0, 1] == 6.0 and mat[1, 1] == -15.0

    def test_symmetrize_symmetric(self):

        mat = TabulaDenseMatrix(3, 3, TabulaSymmetry.symmetric)
        # fill the strict upper triangle
        mat[0, 1] = 2.0
        mat[0, 2] = -1.0
        mat[1, 2] = 5.0
        mat.symmetrize()
        # the lower triangle now mirrors it
        assert mat[1, 0] == 2.0
        assert mat[2, 0] == -1.0
        assert mat[2, 1] == 5.0

    def test_symmetrize_antisymmetric(self):

        mat = TabulaDenseMatrix(3, 3, TabulaSymmetry.antisymmetric)
        mat[0, 1] = 3.0
        mat[0, 2] = -2.0
        mat[1, 2] = 4.0
        # the diagonal is independent storage — symmetrize leaves it untouched
        mat[1, 1] = 9.0
        mat.symmetrize()
        assert mat[1, 0] == -3.0
        assert mat[2, 0] == 2.0
        assert mat[2, 1] == -4.0
        assert mat[1, 1] == 9.0

    def test_to_numpy(self):

        mat = TabulaDenseMatrix(2, 3)
        mat[0, 0] = 1.0
        mat[1, 2] = 6.0
        arr = mat.to_numpy()
        assert arr.shape == (2, 3)
        expected = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 6.0]])
        assert np.allclose(arr, expected, 1.0e-13, 1.0e-13, False)

    def test_to_numpy_is_zero_copy(self):

        mat = TabulaDenseMatrix(2, 2, TabulaSymmetry.symmetric)
        mat[0, 1] = 2.0
        arr = mat.to_numpy()
        # a C++-side mutation is visible through the existing NumPy view
        mat.scale(2.0)
        assert arr[0, 1] == 4.0
        # and a NumPy-side write is visible C++-side
        arr[1, 0] = 7.0
        assert mat[1, 0] == 7.0
