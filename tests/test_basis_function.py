import numpy as np
import pickle

from mpi4py import MPI

from veloxchem.veloxchemlib import BasisFunction


class TestBasisFunction:

    def basis_funtion(self, angmom):

        return BasisFunction([0.8, 1.5, 2.7], [2.4, 0.7, -0.5], angmom)

    def test_pickle(self):

        bf_a = self.basis_funtion(2)
        bobj = pickle.dumps(bf_a)
        bf_b = pickle.loads(bobj)
        assert bf_a == bf_b

    def test_add(self):

        bf_a = self.basis_funtion(2)
        bf_b = BasisFunction([], [], 2)
        bf_b.add(0.8, 2.4)
        bf_b.add(1.5, 0.7)
        bf_b.add(2.7, -0.5)
        assert bf_a == bf_b

    def test_set_exponents(self):

        bf_a = self.basis_funtion(2)
        bf_a.set_exponents([0.2, 0.9, 1.7])
        bf_b = BasisFunction([0.2, 0.9, 1.7], [2.4, 0.7, -0.5], 2)
        assert bf_a == bf_b

    def test_set_normalization_factors(self):

        bf_a = self.basis_funtion(2)
        bf_a.set_normalization_factors([0.2, 0.9, 1.7])
        bf_b = BasisFunction([0.8, 1.5, 2.7], [0.2, 0.9, 1.7], 2)
        assert bf_a == bf_b

    def test_set_angular_momentum(self):

        bf_a = self.basis_funtion(2)
        bf_a.set_angular_momentum(1)
        bf_b = self.basis_funtion(1)
        assert bf_a == bf_b

    def test_normalize(self):

        # testing normalization for s functions
        bf_a = self.basis_funtion(0)
        bf_a.normalize()
        bf_b = BasisFunction(
            [0.800000000000000, 1.50000000000000, 2.700000000000000],
            [0.542250667462678, 0.253418226705534, -0.281296410762365], 0)
        assert bf_a == bf_b

        # testing normalization for p functions
        bf_a = self.basis_funtion(1)
        bf_a.normalize()
        bf_b = BasisFunction(
            [0.800000000000000, 1.50000000000000, 2.700000000000000],
            [0.958298883615608, 0.613252561079874, -0.913275834565931], 1)
        assert bf_a == bf_b

        # testing normalization for d functions
        bf_a = self.basis_funtion(2)
        bf_a.normalize()
        bf_b = BasisFunction(
            [0.800000000000000, 1.50000000000000, 2.700000000000000],
            [0.980790072864359, 0.8594390563363751, -1.7171725369825064], 2)
        assert bf_a == bf_b

        # testing normalization for f functions
        bf_a = self.basis_funtion(3)
        bf_a.normalize()
        bf_b = BasisFunction(
            [0.800000000000000, 1.50000000000000, 2.70000000000000],
            [0.779437206859724, 0.9352350910818505, -2.5070090082858605], 3)
        assert bf_a == bf_b

        # testing normalization for g functions
        bf_a = self.basis_funtion(4)
        bf_a.normalize()
        bf_b = BasisFunction(
            [0.8000000000000000, 1.50000000000000, 2.700000000000000],
            [0.5245190989624026, 0.8617903028794913, -3.099366655142966], 4)
        assert bf_a == bf_b

        # testing normalization for h functions
        bf_a = self.basis_funtion(5)
        bf_a.normalize()
        bf_b = BasisFunction(
            [0.8000000000000000, 1.5000000000000000, 2.700000000000000],
            [0.3117743633857545, 0.7014251498674728, -3.3844565233104653], 5)
        assert bf_a == bf_b

        # testing normalization for i functions
        bf_a = self.basis_funtion(6)
        bf_a.normalize()
        bf_b = BasisFunction(
            [0.8000000000000000, 1.5000000000000000, 2.700000000000000],
            [0.1678350767963739, 0.5170401872124587, -3.3470959216753484], 6)
        assert bf_a == bf_b

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

    def test_mpi_bcast(self):

        comm = MPI.COMM_WORLD

        bf_a = None
        if comm.Get_rank() == 0:
            bf_a = self.basis_funtion(4)
        bf_a = comm.bcast(bf_a)
        bf_b = self.basis_funtion(4)
        assert bf_a == bf_b
