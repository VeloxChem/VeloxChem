from veloxchem.veloxchemlib import upper_case
from veloxchem.veloxchemlib import lower_case


class TestStringFormat:
    """
    Implements tests for src/general/StringFormat.hpp
    """

    def test_upcase(self):

        assert upper_case("VeloxChem-1.0") == "VELOXCHEM-1.0"

    def test_lowcase(self):

        assert lower_case("VeloxChem-1.0") == "veloxchem-1.0"
