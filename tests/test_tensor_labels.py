from veloxchem import tensor_cartesian_labels
from veloxchem import tensor_spherical_labels
from veloxchem import tensor_cartesian_index
from veloxchem import tensor_label
from veloxchem import tensor_order


class TestTensorLabels:
    """
    Implements tests for src/general/TensorLabels.hpp
    """

    def test_tensor_cartesian_labels(self):

        assert tensor_cartesian_labels(0) == []
        assert tensor_cartesian_labels(1) == ["X", "Y", "Z"]
        assert tensor_cartesian_labels(2) == [
            "XX", "XY", "XZ", "YY", "YZ", "ZZ"
        ]
        assert tensor_cartesian_labels(3) == [
            "XXX", "XXY", "XXZ", "XYY", "XYZ", "XZZ", "YYY", "YYZ", "YZZ",
            "ZZZ"
        ]
        assert tensor_cartesian_labels(4) == [
            "XXXX", "XXXY", "XXXZ", "XXYY", "XXYZ", "XXZZ", "XYYY", "XYYZ",
            "XYZZ", "XZZZ", "YYYY", "YYYZ", "YYZZ", "YZZZ", "ZZZZ"
        ]
        assert tensor_cartesian_labels(5) == [
            "XXXXX", "XXXXY", "XXXXZ", "XXXYY", "XXXYZ", "XXXZZ", "XXYYY",
            "XXYYZ", "XXYZZ", "XXZZZ", "XYYYY", "XYYYZ", "XYYZZ", "XYZZZ",
            "XZZZZ", "YYYYY", "YYYYZ", "YYYZZ", "YYZZZ", "YZZZZ", "ZZZZZ"
        ]
        assert tensor_cartesian_labels(6) == [
            "XXXXXX", "XXXXXY", "XXXXXZ", "XXXXYY", "XXXXYZ", "XXXXZZ",
            "XXXYYY", "XXXYYZ", "XXXYZZ", "XXXZZZ", "XXYYYY", "XXYYYZ",
            "XXYYZZ", "XXYZZZ", "XXZZZZ", "XYYYYY", "XYYYYZ", "XYYYZZ",
            "XYYZZZ", "XYZZZZ", "XZZZZZ", "YYYYYY", "YYYYYZ", "YYYYZZ",
            "YYYZZZ", "YYZZZZ", "YZZZZZ", "ZZZZZZ"
        ]

    def test_tensor_spherical_labels(self):

        assert tensor_spherical_labels(0) == [
            "S",
        ]
        assert tensor_spherical_labels(1) == ["P-1", "P0", "P+1"]
        assert tensor_spherical_labels(2) == ["D-2", "D-1", "D0", "D+1", "D+2"]
        assert tensor_spherical_labels(3) == [
            "F-3", "F-2", "F-1", "F0", "F+1", "F+2", "F+3"
        ]
        assert tensor_spherical_labels(4) == [
            "G-4", "G-3", "G-2", "G-1", "G0", "G+1", "G+2", "G+3", "G+4"
        ]
        assert tensor_spherical_labels(5) == [
            "H-5", "H-4", "H-3", "H-2", "H-1", "H0", "H+1", "H+2", "H+3",
            "H+4", "H+5"
        ]
        assert tensor_spherical_labels(6) == [
            "I-6", "I-5", "I-4", "I-3", "I-2", "I-1", "I0", "I+1", "I+2",
            "I+3", "I+4", "I+5", "I+6"
        ]

    def test_tensor_cartesian_index(self):

        # tensor order: 1
        assert tensor_cartesian_index("X") == 0
        assert tensor_cartesian_index("Y") == 1
        assert tensor_cartesian_index("Z") == 2

        # tensor order: 2
        assert tensor_cartesian_index("XX") == 0
        assert tensor_cartesian_index("XY") == 1
        assert tensor_cartesian_index("XZ") == 2
        assert tensor_cartesian_index("YY") == 3
        assert tensor_cartesian_index("YZ") == 4
        assert tensor_cartesian_index("ZZ") == 5

        # tensor order: 3
        assert tensor_cartesian_index("XXX") == 0
        assert tensor_cartesian_index("XXY") == 1
        assert tensor_cartesian_index("XXZ") == 2
        assert tensor_cartesian_index("XYY") == 3
        assert tensor_cartesian_index("XYZ") == 4
        assert tensor_cartesian_index("XZZ") == 5
        assert tensor_cartesian_index("YYY") == 6
        assert tensor_cartesian_index("YYZ") == 7
        assert tensor_cartesian_index("YZZ") == 8
        assert tensor_cartesian_index("ZZZ") == 9

        # tensor order: 4
        assert tensor_cartesian_index("XXXX") == 0
        assert tensor_cartesian_index("XXXY") == 1
        assert tensor_cartesian_index("XXXZ") == 2
        assert tensor_cartesian_index("XXYY") == 3
        assert tensor_cartesian_index("XXYZ") == 4
        assert tensor_cartesian_index("XXZZ") == 5
        assert tensor_cartesian_index("XYYY") == 6
        assert tensor_cartesian_index("XYYZ") == 7
        assert tensor_cartesian_index("XYZZ") == 8
        assert tensor_cartesian_index("XZZZ") == 9
        assert tensor_cartesian_index("YYYY") == 10
        assert tensor_cartesian_index("YYYZ") == 11
        assert tensor_cartesian_index("YYZZ") == 12
        assert tensor_cartesian_index("YZZZ") == 13
        assert tensor_cartesian_index("ZZZZ") == 14

        # tensor order: 5
        assert tensor_cartesian_index("XXXXX") == 0
        assert tensor_cartesian_index("XXXXY") == 1
        assert tensor_cartesian_index("XXXXZ") == 2
        assert tensor_cartesian_index("XXXYY") == 3
        assert tensor_cartesian_index("XXXYZ") == 4
        assert tensor_cartesian_index("XXXZZ") == 5
        assert tensor_cartesian_index("XXYYY") == 6
        assert tensor_cartesian_index("XXYYZ") == 7
        assert tensor_cartesian_index("XXYZZ") == 8
        assert tensor_cartesian_index("XXZZZ") == 9
        assert tensor_cartesian_index("XYYYY") == 10
        assert tensor_cartesian_index("XYYYZ") == 11
        assert tensor_cartesian_index("XYYZZ") == 12
        assert tensor_cartesian_index("XYZZZ") == 13
        assert tensor_cartesian_index("XZZZZ") == 14
        assert tensor_cartesian_index("YYYYY") == 15
        assert tensor_cartesian_index("YYYYZ") == 16
        assert tensor_cartesian_index("YYYZZ") == 17
        assert tensor_cartesian_index("YYZZZ") == 18
        assert tensor_cartesian_index("YZZZZ") == 19
        assert tensor_cartesian_index("ZZZZZ") == 20

        # tensor order: 6
        assert tensor_cartesian_index("XXXXXX") == 0
        assert tensor_cartesian_index("XXXXXY") == 1
        assert tensor_cartesian_index("XXXXXZ") == 2
        assert tensor_cartesian_index("XXXXYY") == 3
        assert tensor_cartesian_index("XXXXYZ") == 4
        assert tensor_cartesian_index("XXXXZZ") == 5
        assert tensor_cartesian_index("XXXYYY") == 6
        assert tensor_cartesian_index("XXXYYZ") == 7
        assert tensor_cartesian_index("XXXYZZ") == 8
        assert tensor_cartesian_index("XXXZZZ") == 9
        assert tensor_cartesian_index("XXYYYY") == 10
        assert tensor_cartesian_index("XXYYYZ") == 11
        assert tensor_cartesian_index("XXYYZZ") == 12
        assert tensor_cartesian_index("XXYZZZ") == 13
        assert tensor_cartesian_index("XXZZZZ") == 14
        assert tensor_cartesian_index("XYYYYY") == 15
        assert tensor_cartesian_index("XYYYYZ") == 16
        assert tensor_cartesian_index("XYYYZZ") == 17
        assert tensor_cartesian_index("XYYZZZ") == 18
        assert tensor_cartesian_index("XYZZZZ") == 19
        assert tensor_cartesian_index("XZZZZZ") == 20
        assert tensor_cartesian_index("YYYYYY") == 21
        assert tensor_cartesian_index("YYYYYZ") == 22
        assert tensor_cartesian_index("YYYYZZ") == 23
        assert tensor_cartesian_index("YYYZZZ") == 24
        assert tensor_cartesian_index("YYZZZZ") == 25
        assert tensor_cartesian_index("YZZZZZ") == 26
        assert tensor_cartesian_index("ZZZZZZ") == 27

    def test_tensor_label(self):

        assert tensor_label(0) == "S"
        assert tensor_label(1) == "P"
        assert tensor_label(2) == "D"
        assert tensor_label(3) == "F"
        assert tensor_label(4) == "G"
        assert tensor_label(5) == "H"
        assert tensor_label(6) == "I"
        assert tensor_label(7) == "K"
        assert tensor_label(8) == "L"
        assert tensor_label(9) == "M"
        assert tensor_label(10) == "N"
        assert tensor_label(11) == "O"
        assert tensor_label(12) == "Q"
        assert tensor_label(13) == "R"
        assert tensor_label(14) == "T"
        assert tensor_label(15) == "U"
        assert tensor_label(16) == "V"

    def test_tensor_order(self):

        assert tensor_order("X") == -1
        assert tensor_order("S") == 0
        assert tensor_order("P") == 1
        assert tensor_order("D") == 2
        assert tensor_order("F") == 3
        assert tensor_order("G") == 4
        assert tensor_order("H") == 5
        assert tensor_order("I") == 6
        assert tensor_order("K") == 7
        assert tensor_order("L") == 8
        assert tensor_order("M") == 9
        assert tensor_order("N") == 10
        assert tensor_order("O") == 11
        assert tensor_order("Q") == 12
        assert tensor_order("R") == 13
        assert tensor_order("T") == 14
        assert tensor_order("U") == 15
        assert tensor_order("V") == 16
