from veloxchem.veloxchemlib import cantor_index
from veloxchem.veloxchemlib import cantor_pair


class TestMathConst:
    """
    Implements tests for src/math/CantorFunc.hpp
    """

    def test_cantor_index(self):

        assert cantor_index((0, 0)) == 0
        assert cantor_index((1, 0)) == 1
        assert cantor_index((0, 1)) == 2
        assert cantor_index((1, 1)) == 4
        assert cantor_index((0, 2)) == 5
        assert cantor_index((2, 0)) == 3
        assert cantor_index((1, 2)) == 8
        assert cantor_index((2, 1)) == 7
        assert cantor_index((2, 2)) == 12

    def test_cantor_apir(self):

        assert cantor_pair(0) == (0, 0)
        assert cantor_pair(1) == (1, 0)
        assert cantor_pair(2) == (0, 1)
        assert cantor_pair(4) == (1, 1)
        assert cantor_pair(5) == (0, 2)
        assert cantor_pair(3) == (2, 0)
        assert cantor_pair(8) == (1, 2)
        assert cantor_pair(7) == (2, 1)
        assert cantor_pair(12) == (2, 2)
