import math as mt

from veloxchem import spherical_momentum_s_factors
from veloxchem import spherical_momentum_p_factors
from veloxchem import spherical_momentum_d_factors
from veloxchem import spherical_momentum_f_factors
from veloxchem import spherical_momentum_g_factors


class TestSphericalMomentum:
    """
    Implements tests for src/orbdata/SphericalMomentum.hpp
    """

    def test_spherical_momentum_s_factors(self):

        tol = 1.0e-12

        facts = spherical_momentum_s_factors(0)
        (idx, f) = facts[0]
        assert idx == 0
        assert mt.isclose(f, 1.0, rel_tol=tol, abs_tol=tol)

    def test_spherical_momentum_p_factors(self):

        tol = 1.0e-12

        facts = spherical_momentum_p_factors(0)
        (idx, f) = facts[0]
        assert idx == 1
        assert mt.isclose(f, 1.0, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_p_factors(1)
        (idx, f) = facts[0]
        assert idx == 2
        assert mt.isclose(f, 1.0, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_p_factors(2)
        (idx, f) = facts[0]
        assert idx == 0
        assert mt.isclose(f, 1.0, rel_tol=tol, abs_tol=tol)

    def test_spherical_momentum_d_factors(self):

        tol = 1.0e-12

        f3 = mt.sqrt(3.0)

        facts = spherical_momentum_d_factors(0)
        (idx, f) = facts[0]
        assert idx == 1
        assert mt.isclose(f, f3, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_d_factors(1)
        (idx, f) = facts[0]
        assert idx == 4
        assert mt.isclose(f, f3, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_d_factors(2)
        (idx, f) = facts[0]
        assert idx == 0
        assert mt.isclose(f, -0.5, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 3
        assert mt.isclose(f, -0.5, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[2]
        assert idx == 5
        assert mt.isclose(f, 1.0, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_d_factors(3)
        (idx, f) = facts[0]
        assert idx == 2
        assert mt.isclose(f, f3, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_d_factors(4)
        (idx, f) = facts[0]
        assert idx == 0
        assert mt.isclose(f, 0.5 * f3, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 3
        assert mt.isclose(f, -0.5 * f3, rel_tol=tol, abs_tol=tol)

    def test_spherical_momentum_f_factors(self):

        tol = 1.0e-12

        f5 = 0.25 * mt.sqrt(10)
        f15 = mt.sqrt(15.0)
        f3 = 0.25 * mt.sqrt(6.0)

        facts = spherical_momentum_f_factors(0)
        (idx, f) = facts[0]
        assert idx == 1
        assert mt.isclose(f, 3.0 * f5, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 6
        assert mt.isclose(f, -f5, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_f_factors(1)
        (idx, f) = facts[0]
        assert idx == 4
        assert mt.isclose(f, f15, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_f_factors(2)
        (idx, f) = facts[0]
        assert idx == 8
        assert mt.isclose(f, 4.0 * f3, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 1
        assert mt.isclose(f, -f3, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[2]
        assert idx == 6
        assert mt.isclose(f, -f3, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_f_factors(3)
        (idx, f) = facts[0]
        assert idx == 9
        assert mt.isclose(f, 1.0, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 2
        assert mt.isclose(f, -1.5, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[2]
        assert idx == 7
        assert mt.isclose(f, -1.5, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_f_factors(4)
        (idx, f) = facts[0]
        assert idx == 5
        assert mt.isclose(f, 4.0 * f3, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 0
        assert mt.isclose(f, -f3, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[2]
        assert idx == 3
        assert mt.isclose(f, -f3, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_f_factors(5)
        (idx, f) = facts[0]
        assert idx == 2
        assert mt.isclose(f, 0.5 * f15, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 7
        assert mt.isclose(f, -0.5 * f15, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_f_factors(6)
        (idx, f) = facts[0]
        assert idx == 0
        assert mt.isclose(f, f5, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 3
        assert mt.isclose(f, -3.0 * f5, rel_tol=tol, abs_tol=tol)

    def test_spherical_momentum_g_factors(self):

        tol = 1.0e-12

        f35 = 0.50 * mt.sqrt(35)
        f17 = 0.25 * mt.sqrt(70)
        f5 = 0.50 * mt.sqrt(5.0)
        f2 = 0.25 * mt.sqrt(10)

        facts = spherical_momentum_g_factors(0)
        (idx, f) = facts[0]
        assert idx == 1
        assert mt.isclose(f, f35, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 6
        assert mt.isclose(f, -f35, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_g_factors(1)
        (idx, f) = facts[0]
        assert idx == 4
        assert mt.isclose(f, 3.0 * f17, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 11
        assert mt.isclose(f, -f17, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_g_factors(2)
        (idx, f) = facts[0]
        assert idx == 8
        assert mt.isclose(f, 6.0 * f5, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 1
        assert mt.isclose(f, -f5, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[2]
        assert idx == 6
        assert mt.isclose(f, -f5, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_g_factors(3)
        (idx, f) = facts[0]
        assert idx == 13
        assert mt.isclose(f, 4.0 * f2, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 4
        assert mt.isclose(f, -3.0 * f2, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[2]
        assert idx == 11
        assert mt.isclose(f, -3.0 * f2, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_g_factors(4)
        (idx, f) = facts[0]
        assert idx == 14
        assert mt.isclose(f, 1.0, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 0
        assert mt.isclose(f, 0.375, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[2]
        assert idx == 10
        assert mt.isclose(f, 0.375, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[3]
        assert idx == 3
        assert mt.isclose(f, 0.75, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[4]
        assert idx == 5
        assert mt.isclose(f, -3.0, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[5]
        assert idx == 12
        assert mt.isclose(f, -3.0, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_g_factors(5)
        (idx, f) = facts[0]
        assert idx == 9
        assert mt.isclose(f, 4.0 * f2, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 2
        assert mt.isclose(f, -3.0 * f2, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[2]
        assert idx == 7
        assert mt.isclose(f, -3.0 * f2, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_g_factors(6)
        (idx, f) = facts[0]
        assert idx == 5
        assert mt.isclose(f, 3.0 * f5, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 12
        assert mt.isclose(f, -3.0 * f5, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[2]
        assert idx == 0
        assert mt.isclose(f, -0.5 * f5, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[3]
        assert idx == 10
        assert mt.isclose(f, 0.5 * f5, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_g_factors(7)
        (idx, f) = facts[0]
        assert idx == 2
        assert mt.isclose(f, f17, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 7
        assert mt.isclose(f, -3.0 * f17, rel_tol=tol, abs_tol=tol)

        facts = spherical_momentum_g_factors(8)
        (idx, f) = facts[0]
        assert idx == 0
        assert mt.isclose(f, 0.25 * f35, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[1]
        assert idx == 10
        assert mt.isclose(f, 0.25 * f35, rel_tol=tol, abs_tol=tol)
        (idx, f) = facts[2]
        assert idx == 3
        assert mt.isclose(f, -1.5 * f35, rel_tol=tol, abs_tol=tol)
