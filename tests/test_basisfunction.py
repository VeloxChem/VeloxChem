import numpy as np

from veloxchem.veloxchemlib import BasisFunction
from tester import Tester


class TestBasisFunction:

    def basis_funtion(self, angmom):

        return BasisFunction([0.8, 1.5, 2.7], [2.4, 0.7, -0.5], angmom)

    def test_add(self):

        # test adding primitves
        bf_a = self.basis_funtion(2)
        bf_b = BasisFunction([], [], 2)
        bf_b.add(0.8, 2.4)
        bf_b.add(1.5, 0.7)
        bf_b.add(2.7, -0.5)
        Tester.compare_basis_functions(bf_a, bf_b)

    def test_set_exponents(self):

        # test setting exponents
        bf_a = self.basis_funtion(2)
        bf_a.set_exponents([0.2, 0.9, 1.7])
        bf_b = BasisFunction([0.2, 0.9, 1.7], [2.4, 0.7, -0.5], 2)
        Tester.compare_basis_functions(bf_a, bf_b)

    def test_set_normalization_factors(self):

        # test setting normalization factors
        bf_a = self.basis_funtion(2)
        bf_a.set_normalization_factors([0.2, 0.9, 1.7])
        bf_b = BasisFunction([0.8, 1.5, 2.7], [0.2, 0.9, 1.7], 2)
        Tester.compare_basis_functions(bf_a, bf_b)

    def test_set_angular_momentum(self):

        # test setting angular momentum
        bf_a = self.basis_funtion(2)
        bf_a.set_angular_momentum(1)
        bf_b = self.basis_funtion(1)
        Tester.compare_basis_functions(bf_a, bf_b)

    def test_normalize(self):

        # testing normalization for s functions
        bf_a = self.basis_funtion(0)
        bf_a.normalize()
        bf_b = BasisFunction(
            [0.800000000000000, 1.50000000000000, 2.700000000000000],
            [0.542250667462678, 0.253418226705534, -0.281296410762365], 0)
        Tester.compare_basis_functions(bf_a, bf_b)

        # testing normalization for p functions
        bf_a = self.basis_funtion(1)
        bf_a.normalize()
        bf_b = BasisFunction(
            [0.800000000000000, 1.50000000000000, 2.700000000000000],
            [0.958298883615608, 0.613252561079874, -0.913275834565931], 1)
        Tester.compare_basis_functions(bf_a, bf_b)

        # testing normalization for d functions
        bf_a = self.basis_funtion(2)
        bf_a.normalize()
        bf_b = BasisFunction(
            [0.800000000000000, 1.50000000000000, 2.700000000000000],
            [0.490395036432180, 0.429719528168190, -0.858586268491259], 2)
        Tester.compare_basis_functions(bf_a, bf_b)

        # testing normalization for f functions
        bf_a = self.basis_funtion(3)
        bf_a.normalize()
        bf_b = BasisFunction(
            [0.800000000000000, 1.50000000000000, 2.70000000000000],
            [0.389718603429863, 0.467617545540924, -1.25350450414293], 3)
        Tester.compare_basis_functions(bf_a, bf_b)

        # testing normalization for g functions
        bf_a = self.basis_funtion(4)
        bf_a.normalize()
        bf_b = BasisFunction(
            [0.8000000000000000, 1.50000000000000, 2.700000000000000],
            [0.0655648873703003, 0.107723787859937, -0.387420831892870], 4)
        Tester.compare_basis_functions(bf_a, bf_b)

        # testing normalization for h functions
        bf_a = self.basis_funtion(5)
        bf_a.normalize()
        bf_b = BasisFunction(
            [0.8000000000000000, 1.5000000000000000, 2.700000000000000],
            [0.0389717954232193, 0.0876781437334342, -0.423057065413809], 5)
        Tester.compare_basis_functions(bf_a, bf_b)

        # testing normalization for i functions
        bf_a = self.basis_funtion(6)
        bf_a.normalize()
        bf_b = BasisFunction(
            [0.8000000000000000, 1.5000000000000000, 2.700000000000000],
            [0.0104896922997734, 0.0323150117007786, -0.209193495104709], 6)
        Tester.compare_basis_functions(bf_a, bf_b)

    def test_get_exponents(self):

        tol = 1.0e-12

        bf = self.basis_funtion(2)
        fe_a = np.array(bf.get_exponents())
        fe_b = np.array([0.8, 1.5, 2.7])
        assert np.allclose(fe_a, fe_b, tol, tol, False)

    def test_get_normalization_factors(self):

        tol = 1.0e-12

        bf = self.basis_funtion(2)
        fn_a = np.array(bf.get_normalization_factors())
        fn_b = np.array([2.4, 0.7, -0.5])
        assert np.allclose(fn_a, fn_b, tol, tol, False)

    def test_get_angular_momentum(self):

        bf = self.basis_funtion(2)
        assert bf.get_angular_momentum() == 2

    def test_number_of_primitives(self):

        bf = self.basis_funtion(2)
        assert bf.number_of_primitives() == 3
