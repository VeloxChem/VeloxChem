import math as mt

from veloxchem.veloxchemlib import Point


class TestPoint:
    """
    Implements tests for src/general/Point.hpp
    """

    def test_scale(self):

        r_a = Point([1.0, 2.2, 3.8])
        r_a.scale(0.5)
        r_b = Point([0.5, 1.1, 1.9])

        assert r_a == r_b

    def test_coordinates(self):

        r_a = Point([1.0, 2.2, 3.8])
        coords = r_a.coordinates()

        tol = 1.0e-12
        assert mt.isclose(coords[0], 1.0, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(coords[1], 2.2, rel_tol=tol, abs_tol=tol)
        assert mt.isclose(coords[2], 3.8, rel_tol=tol, abs_tol=tol)

    def test_length_square(self):

        r_a = Point([1.0, 2.2, 3.8])
        l2 = r_a.length_square()

        tol = 1.0e-12
        assert mt.isclose(l2, 20.28, rel_tol=tol, abs_tol=tol)

    def test_length(self):

        r_a = Point([1.0, 2.2, 3.8])
        l = r_a.length()

        tol = 1.0e-12
        assert mt.isclose(l,
                          4.5033320996790809631713604879153,
                          rel_tol=tol,
                          abs_tol=tol)

    def test_distance_square(self):

        r_a = Point([1.0, 2.2, 3.8])
        r_b = Point([1.0, 1.2, 1.8])
        d2 = r_a.distance_square(r_b)

        tol = 1.0e-12
        assert mt.isclose(d2, 5.0, rel_tol=tol, abs_tol=tol)

    def test_distance(self):

        r_a = Point([1.0, 2.2, 3.8])
        r_b = Point([1.0, 1.2, 1.8])
        d = r_a.distance(r_b)

        tol = 1.0e-12
        assert mt.isclose(d,
                          2.2360679774997896964091736687313,
                          rel_tol=tol,
                          abs_tol=tol)
