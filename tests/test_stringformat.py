from veloxchem.veloxchemlib import fmt
from veloxchem.veloxchemlib import upcase
from veloxchem.veloxchemlib import format
from veloxchem.veloxchemlib import to_string
from veloxchem.veloxchemlib import to_angular_momentum

class TestStringFormat:
    """
    Implements tests for src/general/StringFormat.hpp
    """

    def test_upcase(self):

        assert upcase("VeloxChem-1.0") == "VELOXCHEM-1.0"

    def test_format(self):
        
        assert format("Slate", 9, fmt.center) == "  Slate  "

        assert format("Slate", 9, fmt.left)   == "Slate    "
        
        assert format("Slate", 9, fmt.right)  == "    Slate"
        
        assert format("Slate", 2, fmt.center) == "Sl"
        
        assert format("Slate", 3, fmt.left)   == "Sla"
        
        assert format("Slate", 4, fmt.right)  == "Slat"
        
    def test_to_string_float(self):

        assert to_string( 0.52917721067, 5, 12, fmt.center) == "   0.52918  "

        assert to_string(-0.52917721067, 5, 12, fmt.center) == "  -0.52918  "
        
        assert to_string( 0.52917721067, 5, 12, fmt.left)   == " 0.52918    "

        assert to_string(-0.52917721067, 5, 12, fmt.left)   == "-0.52918    "
        
        assert to_string( 0.52917721067, 5, 12, fmt.right)  == "     0.52918"

        assert to_string(-0.52917721067, 5, 12, fmt.right)  == "    -0.52918"
        
    def test_to_string_float_with_precision(self):

        assert to_string( 0.52917721067, 3) == "0.529"

        assert to_string(-0.52917721067, 3) == "-0.529"
        
        assert to_string( 0.52917721067, 2) == "0.53"

        assert to_string(-0.52917721067, 2) == "-0.53"
        
    def test_to_string_int(self):
        
        assert to_string( 124, 6, fmt.center) == "  124 "
        
        assert to_string(-124, 6, fmt.center) == " -124 "
        
        assert to_string( 124, 6, fmt.left)   == " 124  "
        
        assert to_string(-124, 6, fmt.left)   == "-124  "
        
        assert to_string( 124, 6, fmt.right)  == "   124"
        
        assert to_string(-124, 6, fmt.right)  == "  -124"
        
    def test_to_string_bool(self):
    
        assert to_string(True)  == "True"

        assert to_string(False) == "False"
        
    def test_to_angular_momentum_from_string(self):
    
        assert to_angular_momentum("X") == -1
        
        assert to_angular_momentum("XZ") == -1
        
        assert to_angular_momentum("S") == 0
                
        assert to_angular_momentum("s") == 0
        
        assert to_angular_momentum("P") == 1
                
        assert to_angular_momentum("p") == 1
        
        assert to_angular_momentum("D") == 2
                
        assert to_angular_momentum("d") == 2
        
        assert to_angular_momentum("F") == 3
                
        assert to_angular_momentum("f") == 3
        
        assert to_angular_momentum("G") == 4
                
        assert to_angular_momentum("g") == 4
        
        assert to_angular_momentum("H") == 5
                
        assert to_angular_momentum("h") == 5
        
        assert to_angular_momentum("I") == 6
                
        assert to_angular_momentum("i") == 6
        
    def test_to_angular_momentum_from_int(self):
    
        assert to_angular_momentum(0) == "S"
        
        assert to_angular_momentum(1) == "P"
                
        assert to_angular_momentum(2) == "D"
                        
        assert to_angular_momentum(3) == "F"
        
        assert to_angular_momentum(4) == "G"
                
        assert to_angular_momentum(5) == "H"
                        
        assert to_angular_momentum(6) == "I"
        
        assert to_angular_momentum(8) == ""

