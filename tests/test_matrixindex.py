from veloxchem.veloxchemlib import uplo_index


class TestMatrixIndex:
    """
    Implements tests for src/math/MatrixIndex.hpp
    """

    def test_uplo_index(self):

        assert uplo_index(0, 0) == 0
        assert uplo_index(0, 1) == 1
        assert uplo_index(0, 2) == 3
        assert uplo_index(1, 1) == 2
        assert uplo_index(1, 2) == 4
        assert uplo_index(2, 2) == 5
