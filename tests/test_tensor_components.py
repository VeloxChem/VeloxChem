from veloxchem import number_of_cartesian_components
from veloxchem import number_of_spherical_components


class TestTensorComponents:
    """
    Implements tests for src/general/TensorComponents.hpp
    """

    def test_number_of_cartesian_components(self):

        assert number_of_cartesian_components(0) == 1
        assert number_of_cartesian_components(1) == 3
        assert number_of_cartesian_components(2) == 6
        assert number_of_cartesian_components(3) == 10
        assert number_of_cartesian_components(4) == 15
        assert number_of_cartesian_components(5) == 21
        assert number_of_cartesian_components(6) == 28

        assert number_of_cartesian_components([1, 3]) == 30
        assert number_of_cartesian_components([3, 1]) == 30

        assert number_of_cartesian_components([1, 2, 3]) == 180
        assert number_of_cartesian_components([3, 2, 1]) == 180
        assert number_of_cartesian_components([2, 1, 3]) == 180
        assert number_of_cartesian_components([2, 3, 1]) == 180

        assert number_of_cartesian_components([1, 2, 3, 0]) == 180
        assert number_of_cartesian_components([3, 2, 1, 0]) == 180
        assert number_of_cartesian_components([2, 1, 3, 0]) == 180
        assert number_of_cartesian_components([2, 3, 1, 0]) == 180
        assert number_of_cartesian_components([0, 1, 2, 3]) == 180
        assert number_of_cartesian_components([0, 3, 2, 1]) == 180
        assert number_of_cartesian_components([0, 2, 1, 3]) == 180
        assert number_of_cartesian_components([0, 2, 3, 1]) == 180
        assert number_of_cartesian_components([1, 0, 2, 3]) == 180
        assert number_of_cartesian_components([3, 0, 2, 1]) == 180
        assert number_of_cartesian_components([2, 0, 1, 3]) == 180
        assert number_of_cartesian_components([2, 0, 3, 1]) == 180
        assert number_of_cartesian_components([1, 2, 0, 3]) == 180
        assert number_of_cartesian_components([3, 2, 0, 1]) == 180
        assert number_of_cartesian_components([2, 1, 0, 3]) == 180
        assert number_of_cartesian_components([2, 3, 0, 1]) == 180

    def test_number_of_spherical_components(self):

        assert number_of_spherical_components(0) == 1
        assert number_of_spherical_components(1) == 3
        assert number_of_spherical_components(2) == 5
        assert number_of_spherical_components(3) == 7
        assert number_of_spherical_components(4) == 9
        assert number_of_spherical_components(5) == 11
        assert number_of_spherical_components(6) == 13

        assert number_of_spherical_components([1, 3]) == 21
        assert number_of_spherical_components([3, 1]) == 21

        assert number_of_spherical_components([1, 2, 3]) == 105
        assert number_of_spherical_components([3, 2, 1]) == 105
        assert number_of_spherical_components([2, 1, 3]) == 105
        assert number_of_spherical_components([2, 3, 1]) == 105

        assert number_of_spherical_components([1, 2, 3, 0]) == 105
        assert number_of_spherical_components([3, 2, 1, 0]) == 105
        assert number_of_spherical_components([2, 1, 3, 0]) == 105
        assert number_of_spherical_components([2, 3, 1, 0]) == 105
        assert number_of_spherical_components([0, 1, 2, 3]) == 105
        assert number_of_spherical_components([0, 3, 2, 1]) == 105
        assert number_of_spherical_components([0, 2, 1, 3]) == 105
        assert number_of_spherical_components([0, 2, 3, 1]) == 105
        assert number_of_spherical_components([1, 0, 2, 3]) == 105
        assert number_of_spherical_components([3, 0, 2, 1]) == 105
        assert number_of_spherical_components([2, 0, 1, 3]) == 105
        assert number_of_spherical_components([2, 0, 3, 1]) == 105
        assert number_of_spherical_components([1, 2, 0, 3]) == 105
        assert number_of_spherical_components([3, 2, 0, 1]) == 105
        assert number_of_spherical_components([2, 1, 0, 3]) == 105
        assert number_of_spherical_components([2, 3, 0, 1]) == 105
