import math as mt

from veloxchem.veloxchemlib import bohr_in_angstroms
from veloxchem.veloxchemlib import hartree_in_ev


class TestCodata:
    """
    Implements tests for src/general/Codata.hpp
    """

    def test_bohr_in_angstroms(self):

        tol = 1.0e-12

        assert mt.isclose(bohr_in_angstroms(),
                          0.529177210903,
                          rel_tol=tol,
                          abs_tol=tol)

    def test_hartree_in_ev(self):

        tol = 1.0e-12

        assert mt.isclose(hartree_in_ev(),
                          27.211386245988,
                          rel_tol=tol,
                          abs_tol=tol)
