import math as mt

from veloxchem.veloxchemlib import ChemicalElement


class TestChemicalElement:
    """
    Implements tests for src/moldata/ChemicalElement.hpp
    """

    def test_set_atom_type_with_label(self):

        tol = 1.0e-12

        elem = ChemicalElement()

        elem.set_atom_type("CU")

        assert elem.get_name() == "Cu"

        assert elem.get_identifier() == 29

        assert mt.isclose(elem.get_mass(), 62.929598, rel_tol=tol, abs_tol=tol)

        assert mt.isclose(elem.get_charge(), 29.0, rel_tol=tol, abs_tol=tol)

        assert elem.get_max_angular_momentum() == 2

    def test_set_atom_type_with_identifier(self):

        tol = 1.0e-12

        elem = ChemicalElement()

        elem.set_atom_type(29)

        assert elem.get_name() == "Cu"

        assert elem.get_identifier() == 29

        assert mt.isclose(elem.get_mass(), 62.929598, rel_tol=tol, abs_tol=tol)

        assert mt.isclose(elem.get_charge(), 29.0, rel_tol=tol, abs_tol=tol)

        assert elem.get_max_angular_momentum() == 2

    def test_set_isotope(self):

        tol = 1.0e-12

        elem = ChemicalElement()

        elem.set_atom_type(29)

        assert elem.get_name() == "Cu"

        assert elem.get_identifier() == 29

        assert mt.isclose(elem.get_mass(), 62.929598, rel_tol=tol, abs_tol=tol)

        assert mt.isclose(elem.get_charge(), 29.0, rel_tol=tol, abs_tol=tol)

        assert elem.get_max_angular_momentum() == 2

        elem.set_isotope(65)

        assert elem.get_name() == "Cu"

        assert elem.get_identifier() == 29

        assert mt.isclose(elem.get_mass(), 64.927790, rel_tol=tol, abs_tol=tol)

        assert mt.isclose(elem.get_charge(), 29.0, rel_tol=tol, abs_tol=tol)

        assert elem.get_max_angular_momentum() == 2

    def test_get_name(self):

        elem = ChemicalElement()

        elem.set_atom_type(6)

        assert elem.get_name() == "C"

        elem.set_atom_type(0)

        assert elem.get_name() == "Bq"

    def test_get_identifier(self):

        elem = ChemicalElement()

        elem.set_atom_type("SC")

        assert elem.get_identifier() == 21

        elem.set_atom_type("BQ")

        assert elem.get_identifier() == 0

    def test_get_mass(self):

        tol = 1.0e-12

        elem = ChemicalElement()

        elem.set_atom_type("H")

        assert mt.isclose(elem.get_mass(), 1.007825, rel_tol=tol, abs_tol=tol)

        elem.set_atom_type("BQ")

        assert mt.isclose(elem.get_mass(), 0.0, rel_tol=tol, abs_tol=tol)

    def test_get_charge(self):

        tol = 1.0e-12

        elem = ChemicalElement()

        elem.set_atom_type("N")

        assert mt.isclose(elem.get_charge(), 7.0, rel_tol=tol, abs_tol=tol)

        elem.set_atom_type("BQ")

        assert mt.isclose(elem.get_charge(), 0.0, rel_tol=tol, abs_tol=tol)

    def test_get_max_angular_momentum(self):

        elem = ChemicalElement()

        elem.set_atom_type("H")

        assert elem.get_max_angular_momentum() == 0

        elem.set_atom_type("N")

        assert elem.get_max_angular_momentum() == 1

        elem.set_atom_type("CU")

        assert elem.get_max_angular_momentum() == 2

        elem.set_atom_type("AU")

        assert elem.get_max_angular_momentum() == 3

        elem.set_atom_type("BQ")

        assert elem.get_max_angular_momentum() == -1

    def test_get_max_identifier(self):

        elem = ChemicalElement()

        assert elem.get_max_identifier() == 86
