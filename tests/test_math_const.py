import math as mt

from veloxchem.veloxchemlib import pi_value


class TestMathConst:
    """
    Implements tests for src/math/MathConst.hpp
    """

    def test_pi_value(self):

        tol = 1.0e-12
        assert mt.isclose(pi_value(),
                          3.14159265358979323846,
                          rel_tol=tol,
                          abs_tol=tol)
