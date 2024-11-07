import math as mt

from veloxchem.veloxchemlib import is_chemical_element
from veloxchem.veloxchemlib import chemical_element_name
from veloxchem.veloxchemlib import chemical_element_label
from veloxchem.veloxchemlib import chemical_element_identifier
from veloxchem.veloxchemlib import chemical_element_mass
from veloxchem.veloxchemlib import chemical_element_max_angular_momentum
from veloxchem.veloxchemlib import chemical_element_max_identifier


class TestChemicalElement:
    """
    Implements tests for src/moldata/ChemicalElement.hpp
    """

    def test_is_chemical_element(self):

        assert is_chemical_element(0) == True
        assert is_chemical_element(6) == True
        assert is_chemical_element(86) == True

        assert is_chemical_element(-1) == False
        assert is_chemical_element(87) == False

    def test_chemical_element_name(self):

        assert chemical_element_name(0) == "BQ"
        assert chemical_element_name(7) == "N"
        assert chemical_element_name(29) == "CU"

    def test_chemical_element_label(self):

        assert chemical_element_label(0) == "Bq"
        assert chemical_element_label(7) == "N"
        assert chemical_element_label(29) == "Cu"

    def test_chemical_element_identifier(self):

        assert chemical_element_identifier("BQ") == 0
        assert chemical_element_identifier("N") == 7
        assert chemical_element_identifier("CU") == 29

    def test_chemical_element_mass(self):

        tol = 1.0e-12
        assert mt.isclose(chemical_element_mass(0),
                          0.0,
                          rel_tol=tol,
                          abs_tol=tol)
        assert mt.isclose(chemical_element_mass(1),
                          1.007825,
                          rel_tol=tol,
                          abs_tol=tol)

    def test_chemical_element_max_angular_momentum(self):

        assert chemical_element_max_angular_momentum(0) == -1
        assert chemical_element_max_angular_momentum(1) == 0
        assert chemical_element_max_angular_momentum(7) == 1
        assert chemical_element_max_angular_momentum(29) == 2
        assert chemical_element_max_angular_momentum(79) == 3

    def test_chemical_element_max_identifier(self):

        assert chemical_element_max_identifier() == 86
