import math as mt

from veloxchem.veloxchemlib import get_pi


class TestMathConst:
    """
    Implements tests for src/math/MathConst.hpp
    """

    def test_get_pi(self):

        tol = 1.0e-12

        assert mt.isclose(get_pi(),
                          3.14159265358979323846,
                          rel_tol=tol,
                          abs_tol=tol)
