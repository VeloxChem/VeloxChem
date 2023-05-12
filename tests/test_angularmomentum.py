from veloxchem.veloxchemlib import to_spherical_components
from veloxchem.veloxchemlib import to_cartesian_components
from veloxchem.veloxchemlib import angular_component_to_str


class TestAngularMomentum:
    """
    Implements tests for src/math/AngularMomentum.hpp
    """

    def test_to_spherical_components(self):
    
        assert to_spherical_components(0) == 1
        assert to_spherical_components(1) == 3
        assert to_spherical_components(2) == 5
        assert to_spherical_components(3) == 7
        assert to_spherical_components(4) == 9
        assert to_spherical_components(5) == 11
        assert to_spherical_components(6) == 13
        
        assert to_spherical_components(0, 0) == 1
        assert to_spherical_components(0, 1) == 3
        assert to_spherical_components(1, 0) == 3
        assert to_spherical_components(0, 2) == 5
        assert to_spherical_components(2, 0) == 5
        assert to_spherical_components(1, 1) == 9
        assert to_spherical_components(1, 2) == 15
        assert to_spherical_components(2, 1) == 15
        assert to_spherical_components(2, 2) == 25
        
    def test_to_cartesian_components(self):
    
        assert to_cartesian_components(0) == 1
        assert to_cartesian_components(1) == 3
        assert to_cartesian_components(2) == 6
        assert to_cartesian_components(3) == 10
        assert to_cartesian_components(4) == 15
        assert to_cartesian_components(5) == 21
        assert to_cartesian_components(6) == 28
        
        assert to_cartesian_components(0, 0) == 1
        assert to_cartesian_components(0, 1) == 3
        assert to_cartesian_components(1, 0) == 3
        assert to_cartesian_components(0, 2) == 6
        assert to_cartesian_components(2, 0) == 6
        assert to_cartesian_components(1, 1) == 9
        assert to_cartesian_components(1, 2) == 18
        assert to_cartesian_components(2, 1) == 18
        assert to_cartesian_components(2, 2) == 36
        
    def test_angular_component_to_str(self):
        
        assert angular_component_to_str(0, 0) == 's  '
        
        assert angular_component_to_str(1, 0) == 'p-1'
        assert angular_component_to_str(1, 1) == 'p0 '
        assert angular_component_to_str(1, 2) == 'p+1'
   
        assert angular_component_to_str(2, 0) == 'd-2'
        assert angular_component_to_str(2, 1) == 'd-1'
        assert angular_component_to_str(2, 2) == 'd0 '
        assert angular_component_to_str(2, 3) == 'd+1'
        assert angular_component_to_str(2, 4) == 'd+2'
        
        assert angular_component_to_str(3, 0) == 'f-3'
        assert angular_component_to_str(3, 1) == 'f-2'
        assert angular_component_to_str(3, 2) == 'f-1'
        assert angular_component_to_str(3, 3) == 'f0 '
        assert angular_component_to_str(3, 4) == 'f+1'
        assert angular_component_to_str(3, 5) == 'f+2'
        assert angular_component_to_str(3, 6) == 'f+3'
        
        assert angular_component_to_str(4, 0) == 'g-4'
        assert angular_component_to_str(4, 1) == 'g-3'
        assert angular_component_to_str(4, 2) == 'g-2'
        assert angular_component_to_str(4, 3) == 'g-1'
        assert angular_component_to_str(4, 4) == 'g0 '
        assert angular_component_to_str(4, 5) == 'g+1'
        assert angular_component_to_str(4, 6) == 'g+2'
        assert angular_component_to_str(4, 7) == 'g+3'
        assert angular_component_to_str(4, 8) == 'g+4'
        
        assert angular_component_to_str(5, 0) == 'h-5'
        assert angular_component_to_str(5, 1) == 'h-4'
        assert angular_component_to_str(5, 2) == 'h-3'
        assert angular_component_to_str(5, 3) == 'h-2'
        assert angular_component_to_str(5, 4) == 'h-1'
        assert angular_component_to_str(5, 5) == 'h0 '
        assert angular_component_to_str(5, 6) == 'h+1'
        assert angular_component_to_str(5, 7) == 'h+2'
        assert angular_component_to_str(5, 8) == 'h+3'
        assert angular_component_to_str(5, 9) == 'h+4'
        assert angular_component_to_str(5, 10) == 'h+5'
        
        assert angular_component_to_str(6, 0) == 'i-6'
        assert angular_component_to_str(6, 1) == 'i-5'
        assert angular_component_to_str(6, 2) == 'i-4'
        assert angular_component_to_str(6, 3) == 'i-3'
        assert angular_component_to_str(6, 4) == 'i-2'
        assert angular_component_to_str(6, 5) == 'i-1'
        assert angular_component_to_str(6, 6) == 'i0 '
        assert angular_component_to_str(6, 7) == 'i+1'
        assert angular_component_to_str(6, 8) == 'i+2'
        assert angular_component_to_str(6, 9) == 'i+3'
        assert angular_component_to_str(6, 10) == 'i+4'
        assert angular_component_to_str(6, 11) == 'i+5'
        assert angular_component_to_str(6, 12) == 'i+6'
