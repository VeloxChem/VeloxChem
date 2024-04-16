import numpy as np
from itertools import permutations

from .veloxchemlib import bohr_in_angstrom
from .molecule import Molecule

class SymmetryAnalyzer:
    """
    A class to determine the point group of nearly symmetric molecules within selected
    tolerance tresholds (in Schoenflies convention) and then to symmetrize it perfectly
    in the detected point group or in one od its subgroups.

    1. The coordinates are reexpressed in the center of mass frame.
    2. The principal axis and moment are obtained by computing the eigenvalues and 
       eigenvectors of the inertia tensor.
        - Linear molecule have one vanishing eigenvalue.
        - Asymmetric top molecules have non-degenerated eigenvalues.
        - Symmetric molecules have doubly degenerated eigenvalues.
        - Spherical molecules have triply degenerated eigenvalues.
    3. Each category is then treated following a typical detection tree.

    The available tolerance parameters are: loose, tight, very tight.

    Instance variables
        - ref_orientation: the reference frame.
        - schoenflies_symbol: the pointgroup schoenflies symbol.
        - max_order: the maximum proper rotational axis order. 
        - list_of_elements: a list of the detected symmetry elements.

    Note: based on the open source pointgroup 0.4.1 package.
    """

    def __init__(self): 
        """
        Initializes the SymmetryAnalyzer instances.
        """

        # Initializes the orientation of the reference frame
        self._ref_orientation = np.identity(3)

        # Initializes variables
        self._schoenflies_symbol = ''
        self._max_order = 1
        self._list_of_elements = ["E"]
        
    def identify_pointgroup(self,
                            molecule,
                            tolerance = 'tight'):
        """
        Analyze the symmetry of a given molecule based on its nuclear framework.

        :param molecule:
            A VeloxChem molecule object.

        :param tolerance:
            A tolerance threshold string available in tolerance_keywords.

        :return:
            A dictionary containing:
                - The pointgroup symbol.
                - (The list of detected elements in spectific cases.)
                - The list of all expected elements.
                - The list of subgroups available for symmetrization.
                - An array with the reoriented geometry in bohr.
        """

        # Define the tolerance parameters
        tolerance = tolerance.lower()

        if tolerance == 'loose':
            self._tolerance_eig = 0.8
            self._tolerance_ang = 0.085
        elif tolerance == 'tight':
            self._tolerance_eig = 0.4
            self._tolerance_ang = 0.070    # about 4 degrees   
        elif tolerance == 'very tight':
            self._tolerance_eig = 0.002
            self._tolerance_ang = 0.035     # about 2 degrees

        # Read and express geometry in center of mass (COM) frame
        # Geom and COM in bohr because moments of inertia are defined in bohr in Molecule module
        coordinates = Molecule.get_coordinates_in_bohr(molecule)
        self._symbols = Molecule.get_labels(molecule)
        self._natoms = len(self._symbols)
        center_of_mass = Molecule.center_of_mass_in_bohr(molecule)
        self._cent_coord = coordinates - center_of_mass

        # Get the principal momemts amd axes of inertia
        Ivals, Ivecs = Molecule.moments_of_inertia(molecule)
        self._Ivals = Ivals
        self._Ivecs = Ivecs

        # Initializes ensembles of results
        symmetry_analysis = {}
        symmetry_elements_list = []
        subgroups_list = []
        self._reoriented_coordinates = np.zeros((self._natoms,3))

        # Handle case of isolated atom
        if self._natoms == 1:
            self._schoenflies_symbol = "O(3)"
            symmetry_analysis["degeneracy"] = "Input structure is an isolated atom."

        else:
            # Get the degeneracy of the eigenvalues of the inertia tensor
            eig_degeneracy = get_degeneracy(self._Ivals, self._tolerance_eig)

            # Linear groups
            if np.min(abs(self._Ivals)) < self._tolerance_eig:
                self._linear()
                symmetry_analysis["degeneracy"] = "Molecule is linear."

            # Asymmetric group
            elif eig_degeneracy == 1:
                self._asymmetric()
                symmetry_analysis["degeneracy"] = "Principal moments of inertia are not degenerated."

            # Symmetric group
            elif eig_degeneracy == 2:
                self._symmetric()
                symmetry_analysis["degeneracy"] = "Principal moments of inertia are doubly degenerated."

            # Spherical group
            elif eig_degeneracy == 3:
                self._spherical()
                symmetry_analysis["degeneracy"] = "Principal moments of inertia are triply degenerated."

        # Collect the results
        if self._max_order <= 2 and self._schoenflies_symbol != "O(3)":
            symmetry_analysis["Elements_found"] = self._list_of_elements

        symmetry_analysis["Point_group"] = self._schoenflies_symbol

        if self._schoenflies_symbol in all_symmetry_elements:
            symmetry_elements_list = all_symmetry_elements[self._schoenflies_symbol]
            symmetry_analysis["Expected_elements"] = symmetry_elements_list

        if self._schoenflies_symbol in subgroups:
            subgroups_list = subgroups[self._schoenflies_symbol]
            symmetry_analysis["Subgroups"] = subgroups_list

        symmetry_analysis["reoriented_geom"] = self._reoriented_coordinates

        return symmetry_analysis
    
    def print_symmetry_results(self, results_dict):
        """
        Build a small output with the results of the symmetry analysis.

        :param results_dict:
            The dictionary containing the different results from the pointgroup_identify function.
        
        :return:
            The results in a visual friendly output.
        """
        
        print(results_dict["degeneracy"],"\n")

        if "Elements_found" in results_dict:
            print("Detected elements: {}\n".format(results_dict["Elements_found"]))
        
        print("Point group: {}\n".format(results_dict["Point_group"]))

        if "Expected_elements" in results_dict:
            print("All expected symmetry elements: {}".format(results_dict["Expected_elements"]))

        if "Subgroups" in results_dict:
            print("Available Abelian subgroups for symmetrization: {}".format(results_dict["Subgroups"]))
        else:
            print("No available subgroups for symmetrization.")

    def print_tolerance_keywords(self):
        """
        :return:
            The available tolerance keywords.
        """

        tolerance_keywords = ["loose", "tight", "very tight"]
        return tolerance_keywords

    def symmetrize_pointgroup(self, results_dict, pointgoup_to_symmetrize=None):
        """
        Symmetrize and reorient molecules. Only reorient if more than 50 atoms.

        :param results_dict:
            The dictionary containing the different results from the pointgroup_identify function.

        !!!Under construction!!!
        """

        # Initializes
        if pointgoup_to_symmetrize == None:
            pointgoup_to_symmetrize = results_dict["Point_group"]

        symmetrized_data = {}

        # Store the geometry reoriented in the desired reference frame for symmetrization
        symmetrized_data["reoriented_geom_in_angstrom"] = results_dict["reoriented_geom"] * bohr_in_angstrom()

        symmetrized_data["new_pointgroup"] = pointgoup_to_symmetrize

        # if self._natoms <= 50:
            

        return symmetrized_data

    def print_symmetrized_molecule(self, results_dict, xyz_filename=None):
        """
        Build the output with the results of the symmetrized molecule.
        The user can choose to print the new coordinates or create a new xyz file.

        :param results_dict:
            The dictionary containing the different results from the pointgroup_symmetrize function.
        
        :param xyz_filename:
            The name of the xyz file (without the .xyz extension).
        
        :return:
            The symmetrized coordinates (in a xyz file if a name is specified).
        """

        if "symmetrized_coord" in results_dict:
            coord = results_dict["symmetrized_coord"]
        else:
            coord = results_dict["reoriented_geom_in_angstrom"]

        pointgroup = results_dict["new_pointgroup"]

        if xyz_filename:
            write_xyz_file(self._symbols,coord,pointgroup,xyz_filename)
        else:
            print(coord)

    def _linear(self):
        """
        Handle linear molecules.

        :return:
            The Schoenflies symbol, list of detected elements, and reoriented coordinates.
        """

        # Set orientation
        idx = np.argmin(self._Ivals)
        main_axis = self._Ivecs[idx]
        p_axis = get_perpendicular(main_axis)
        self._set_orientation(main_axis, p_axis)

        # Check for inversion center at center of mass
        if self._check_op(Inversion()):
            self._schoenflies_symbol = 'Dinfh'
            self._list_of_elements.append("i")
        else:
            self._schoenflies_symbol = 'Cinfv'
        
        # Save reoriented center of mass frame coordinates for symmetrization
        main_axis = [0.,0.,1.]
        p_axis = get_perpendicular(main_axis)
        self._set_orientation(main_axis, p_axis)
        self._reoriented_coordinates = self._cent_coord

    def _asymmetric(self):
        """
        Handle asymmetric top molecules.
        """

        # Set orientation
        self._set_orientation(self._Ivecs[0], self._Ivecs[1])

        # Check for any C2 axis
        n_axis_c2 = 0
        main_axis = [1, 0, 0]
        for axis in np.identity(3):
            c2 = Rotation(axis, order=2)
            if self._check_op(c2):
                n_axis_c2 += 1
                main_axis = axis
        
        self._max_order = 2

        # Check and save in list the rotation axis of order 2 along the main axis
        if self._check_op(Rotation(main_axis, self._max_order)):
            self._list_of_elements.append("C{}".format(self._max_order))

        if n_axis_c2 == 0:
            self._max_order = 0
            self._no_rot_axis()
        elif n_axis_c2 == 1:
            self._cyclic(main_axis)
        else:
            self._dihedral(main_axis)

    def _symmetric(self):
            """
            Handle symmetric molecules.
            """

            # Get the only non-degenareted principal moment fo inertia and set the main axis
            # along the associated eigenvector
            idx = get_non_degenerated(self._Ivals, self._tolerance_eig)
            main_axis = self._Ivecs[idx]

            # Determine the highest possible rotation axis order along the main axis
            self._max_order = self._get_axis_rot_order(main_axis, n_max=9)

            # Check for rotation of maximum order
            if self._check_op(Rotation(main_axis, order=self._max_order)):
                self._list_of_elements.append("C{}".format(self._max_order))

            # Get the perpendicualar axis to main axis and check for C2 rotation axis 
            # along p_axis by rotating p_axis along the main axis
            p_axis = get_perpendicular(main_axis)
            for angle in np.arange(0, np.pi, 0.1* np.pi / self._max_order):
                axis = np.dot(p_axis, rotation_matrix(main_axis, angle))
                c2 = Rotation(axis, order=2)
                if self._check_op(c2):
                    self._dihedral(main_axis)
                    return
                        
            self._cyclic(main_axis)

    def _spherical(self):
        """
        Handle spherical groups (I, O, T) in iterative way by increasing tolerance if no axis if found.
        """

        main_axis = None
        while main_axis is None:
            for axis in get_cubed_sphere_grid_points(self._tolerance_ang):
                c5 = Rotation(axis, order=5)
                c4 = Rotation(axis, order=4)
                c3 = Rotation(axis, order=3)

                # Check for C5 axis
                if self._check_op(c5, tol_factor=1.0):
                    self._schoenflies_symbol = "I"
                    main_axis = axis
                    self._max_order = 5
                    break
                # Check for C4 axis
                elif self._check_op(c4, tol_factor=1.5):
                    self._schoenflies_symbol = "O"
                    main_axis = axis
                    self._max_order = 4
                    break
                # Check for C3 axis
                elif self._check_op(c3, tol_factor=1.5):
                    self._schoenflies_symbol = "T"
                    main_axis = axis
                    self._max_order = 3

            if main_axis is None:
                print('Increase of the angular tolerance.')
                self._tolerance_ang *= 1.05
                print("New angular tolerance: {} degrees\n".format(self._tolerance_ang))

        p_axis_base = get_perpendicular(main_axis)

        # Save reoriented center of mass frame coordinates for symmetrization
        # Main axis is z by convention
        orientation = np.array([p_axis_base, np.cross(main_axis, p_axis_base), main_axis])
        self._reoriented_coordinates = np.dot(self._cent_coord, orientation.T)

        # I or Ih
        if self._schoenflies_symbol == 'I':

            def _determine_orientation_I(main_axis):
                """
                Determine the orientation of the molecule within the icosahedral symmetry group.
                The angle is derived from geometric considerations to align the orientation
                axis with the symmetry axis.

                :param main_axis:
                    The main axis obtained with _set_orientation. 

                :return:
                    The current orientation in the reference frame.
                """
                
                r_matrix = rotation_matrix(p_axis_base, np.arcsin((np.sqrt(5)+1)/(2*np.sqrt(3))))
                axis = np.dot(main_axis, r_matrix.T)

                # set molecule orientation in I
                for angle in np.arange(0, 2*np.pi+self._tolerance_ang, self._tolerance_ang):
                    rot_matrix = rotation_matrix(main_axis, angle)

                    c5_axis = np.dot(axis, rot_matrix.T)
                    c5 = Rotation(c5_axis, order=5)

                    if self._check_op(c5):
                        t_axis = np.dot(main_axis, rotation_matrix(p_axis_base, np.pi/2).T)
                        return np.dot(t_axis, rot_matrix.T)

            p_axis = _determine_orientation_I(main_axis)
            # Set orientation
            self._set_orientation(main_axis, p_axis)

            # Check for inverison center at center of mass
            if self._check_op(Inversion()):
                self._schoenflies_symbol += 'h'

        # O or Oh
        if self._schoenflies_symbol == 'O':
            
            def _determine_orientation_O(main_axis):
                """
                Determine the orientation of the molecule within the Octahedron symmetry group. 

                :param main_axis:
                    The main axis obtained with _set_orientation. 

                :return:
                    The current orientation in the reference frame.
                """

                r_matrix = rotation_matrix(p_axis_base, np.pi/2)
                axis = np.dot(main_axis, r_matrix.T)

                for angle in np.arange(0, 2*np.pi / self._max_order+self._tolerance_ang, self._tolerance_ang):
                    rot_matrix = rotation_matrix(main_axis, angle)

                    c4_axis = np.dot(axis, rot_matrix.T)
                    c4 = Rotation(c4_axis, order=4)

                    if self._check_op(c4):
                        return axis

            p_axis = _determine_orientation_O(main_axis)
            # Set orientation
            self._set_orientation(main_axis, p_axis)

            # Check for inverison center at center of mass
            if self._check_op(Inversion()):
                self._schoenflies_symbol += 'h'

        # T or Td, Th
        if self._schoenflies_symbol == 'T':

            def _determine_orientation_T(main_axis):
                """
                Determine the orientation of the molecule within the Tetrahedron symmetry group. 

                :param main_axis:
                    The main axis obtained with _set_orientation. 

                :return:
                    The current orientation in the reference frame.
                """

                r_matrix = rotation_matrix(p_axis_base, -np.arccos(-1/3))
                axis = np.dot(main_axis, r_matrix.T)

                for angle in np.arange(0, 2*np.pi / self._max_order + self._tolerance_ang, self._tolerance_ang):
                    rot_matrix = rotation_matrix(main_axis, angle)

                    c3_axis = np.dot(axis, rot_matrix.T)
                    c3 = Rotation(c3_axis, order=3)

                    if self._check_op(c3):
                        t_axis = np.dot(main_axis, rotation_matrix(p_axis_base, np.pi/2).T)
                        return np.dot(t_axis, rot_matrix.T)

            p_axis = _determine_orientation_T(main_axis)
            # Set orientation
            self._set_orientation(main_axis, p_axis)

            # Check for inverison center at center of mass
            if self._check_op(Inversion()):
                self._schoenflies_symbol += 'h'
                return
            
            # Check for any reflexion plane 
            if self._check_op(Reflection([0, 0, 1])):
                self._schoenflies_symbol += 'd'
                return

    def _no_rot_axis(self):
        """
        Detect point group of molecule with no rotation axis.

        :return:
            The Schoenflies symbol and list of detected elements.
        """

        for i, vector in enumerate(np.identity(3)):
            if self._check_op(Reflection(vector)):
                self._list_of_elements.append("sigma")
                self._schoenflies_symbol = 'Cs'
                p_axis = get_perpendicular(vector)
                self._set_orientation(vector, p_axis)
                break
            else:
                if self._check_op(Inversion()):
                    self._schoenflies_symbol = 'Ci'
                    self._list_of_elements.append("i")
                    break
                else:
                    self._schoenflies_symbol = 'C1'
        
        # Save reoriented center of mass frame coordinates for symmetrization
        self._reoriented_coordinates = self._cent_coord

    def _cyclic(self, main_axis):
        """
        Detect point group of cylic group molecules.

        :param main_axis:
            The main axis obtained with _set_orientation.

        :return:
            The Schoenflies symbol and list of detected elements.
        """

        self._schoenflies_symbol = "C{}".format(self._max_order)
        
        # Check for reflexion planes perpenducular to the main axis
        if self._check_op(Reflection(main_axis)):
            self._schoenflies_symbol += 'h'
            self._list_of_elements.append("sigma_h")

        # Check for reflexion planes containing the main axis
        v_symbol = set()
        p_axis = get_perpendicular(main_axis)
        for angle in np.arange(0, np.pi + self._tolerance_ang, 0.5*np.pi / self._max_order + self._tolerance_ang):
            axis = np.dot(p_axis, rotation_matrix(main_axis, angle))
            if self._check_op(Reflection(axis)):
                self._list_of_elements.append("sigma_v")
                v_symbol.add('v')

        self._schoenflies_symbol += ''.join(v_symbol)

        # Check for inversion center at center of mass
        if self._check_op(Inversion()):
            self._list_of_elements.append("i")

        # Check for improper rotational axis along the main axis.
        if self._check_op(ImproperRotation(main_axis, order=2*self._max_order)):
            self._list_of_elements.append("S{}".format(2*self._max_order))
            self._schoenflies_symbol = "S{}".format(2 * self._max_order)

        # Save reoriented center of mass frame coordinates for symmetrization
        # Main axis is z by convention
        orientation = np.array([p_axis, np.cross(main_axis, p_axis), main_axis])
        self._cent_coord = np.dot(self._cent_coord, orientation.T)
        self._reoriented_coordinates = self._cent_coord
            
    def _dihedral(self, main_axis):
        """
        Detect point group of dihedral group molecules.

        :param main_axis:
            The main axis obtained with _set_orientation.

        :return:
            The Schoenflies symbol and list of detected elements.
        """

        # Determine perpendicular axis to main axis
        p_axis = get_perpendicular(main_axis)

        if self._max_order == 1:
            # D1 is equivalent to C2
            self._schoenflies_symbol = "C2"
        else:
            self._schoenflies_symbol = "D{}".format(self._max_order)

        # Check for C2 axis perpendicular to the main axis 
        if self._max_order == 2:
            for angle in np.arange(0, np.pi + self._tolerance_ang, np.pi / self._max_order + self._tolerance_ang):
                axis = np.dot(p_axis, rotation_matrix(main_axis, angle))
                if self._check_op(Rotation(axis, order=self._max_order)):
                    self._list_of_elements.append("C{}".format(self._max_order))

        # Check for inversion center at center of mass
        if self._check_op(Inversion()):
            self._list_of_elements.append("i")

        # Check for reflexion planes perpenducular to the main axis
        h_symbols = False
        if self._check_op(Reflection(main_axis)):
            self._list_of_elements.append("sigma_h")
            h_symbols = True
            if h_symbols:
                self._schoenflies_symbol += 'h'
        
        # Check for reflexion planes containing the main axis
        d_symbol = False
        for angle in np.arange(0, np.pi, 0.5*np.pi / self._max_order):
            axis = np.dot(p_axis, rotation_matrix(main_axis, angle))
            if self._check_op(Reflection(axis)):
                self._list_of_elements.append("sigma")
                if not h_symbols:
                    if not d_symbol:
                        self._schoenflies_symbol += 'd'
                        d_symbol = True
                        
        # Check for improper rotational axis along the main axis.
        if self._check_op(ImproperRotation(main_axis, 2*self._max_order)):
            self._list_of_elements.append("S{}".format(2*self._max_order))

        # Save reoriented center of mass frame coordinates for symmetrization
        # Main axis is z by convention
        orientation = np.array([p_axis, np.cross(main_axis, p_axis), main_axis])
        self._cent_coord = np.dot(self._cent_coord, orientation.T)
        self._reoriented_coordinates = self._cent_coord
    
    def _get_axis_rot_order(self, axis, n_max):
            """
            Get rotation order for a given axis.

            :param axis:
                The axis

            :param n_max:
                Maximum order to scan

            :return:
                The order
            """

            def max_rotation_order(tolerance):
                """
                Set the range of maximum order possible.

                :param tolerance:
                    The tolerance parameter.
                """

                for i in range(2, 15):
                    if 2*np.pi / (i * (i - 1)) <= tolerance:
                        return i-1

            n_max = np.min([max_rotation_order(self._tolerance_ang), n_max])

            for i in range(n_max, 1, -1):
                Cn = Rotation(axis, order=i)
                if self._check_op(Cn):
                    return i
            return 1

    def _check_op(self, operation, tol_factor=1.0):
        """
        Check if the given symmetry operation exists in the point group of the molecule.

        :param operation:
            The matrix representative of the symmetry operation.

        :param tol_factor:
            A factor to scale the tolerance value by.

        :return:
            True if the symmetry operation exists in the point group, False otherwise.
        """

        sym_matrix = operation.get_matrix()
        error_abs_rad = abs_to_rad(self._tolerance_eig, coord=self._cent_coord)

        op_coordinates = np.matmul(self._cent_coord, sym_matrix)
        for idx, op_coord in enumerate(op_coordinates):

            difference_rad = radius_diff_in_radiants(op_coord, self._cent_coord, self._tolerance_eig)
            difference_ang = angle_between_vector_matrix(op_coord, self._cent_coord, self._tolerance_eig)

            def check_diff(diff, diff2):
                for idx_2, (d1, d2) in enumerate(zip(diff, diff2)):
                    if self._symbols[idx_2] != self._symbols[idx]:
                        continue
                    tolerance_total = self._tolerance_ang * tol_factor + error_abs_rad[idx_2]
                    if d1 < tolerance_total and d2 < tolerance_total:
                        return True
                return False

            if not check_diff(difference_ang, difference_rad):
                return False

        return True

    def _set_orientation(self, main_axis, p_axis):
        """
        Set molecular orientation along main_axis (x) and p_axis (y).

        :param main_axis:
            Principal orientation axis (must be unitary)

        :param p_axis:
            Secondary axis perpendicular to principal (must be unitary)
        """

        assert np.linalg.norm(main_axis) > 1e-1
        assert np.linalg.norm(p_axis) > 1e-1

        orientation = np.array([main_axis, p_axis, np.cross(main_axis, p_axis)])
        self._cent_coord = np.dot(self._cent_coord, orientation.T)
        self._ref_orientation = np.dot(self._ref_orientation, orientation.T)



"""
The followig classes determine the representatives of the symmetry operations
(to the n^th order, if necessary).
"""
class Inversion:
    def get_matrix(self):
        return -np.identity(3)
    
class Rotation:
    def __init__(self, axis, order=1):
        self._order = order

        self._axis = np.array(axis)

    def get_matrix(self):

        return rotation_matrix(self._axis, 2*np.pi / self._order)
    
class Reflection:
    def __init__(self, axis):

        norm = np.linalg.norm(axis)
        assert abs(norm) > 1e-8
        self._axis = np.array(axis) / norm   # normalize axis

    def get_matrix(self):
        uax = np.dot(self._axis, self._axis)

        return np.identity(3) - 2*np.outer(self._axis, self._axis)/uax
    
class ImproperRotation:
    def __init__(self, axis, order=1):
        self._order = order

        self._axis = np.array(axis)

    def get_matrix(self):

        rot_matrix = rotation_matrix(self._axis, 2*np.pi / self._order)

        uax = np.dot(self._axis, self._axis)
        refl_matrix = np.identity(3) - 2 * np.outer(self._axis, self._axis) / uax

        return np.dot(rot_matrix, refl_matrix.T)



@staticmethod
def get_degeneracy(Ivals, tolerance):
        """
        Get the degeneracy of the principal inertia moments.

        :param Ivals:
            The array of eigenvalues.

        :param tolerance:
            The tolerance parameter on the eigenvalues.

        :return:
            The degree of degeneracy.
        """

        for ev1 in Ivals:
            single_deg = 0
            for ev2 in Ivals:
                if abs(ev1 - ev2) < tolerance:
                    single_deg += 1
            if single_deg > 1:
                return single_deg
        return 1

@staticmethod
def get_perpendicular(vector, tol=1e-8):
    """
    Generate a vector perpendicular to another vector or axis. 

    :param vector:
        The vector or axis with respect to which the perpendicular axis is determined.

    :param tol:
        An additional tolerance parameter to condisder the axis as perpendicular.

    :return:
        An array of coordinates of the perpendicular and normalized vector.
    """

    index = np.argmin(np.abs(vector))
    p_vector = np.identity(3)[index]
    pp_vector = np.cross(vector, p_vector)
    pp_vector = pp_vector / np.linalg.norm(pp_vector)

    assert np.dot(pp_vector, vector) < tol  # check perpendicular
    assert abs(np.linalg.norm(pp_vector) - 1) < tol  # check normalized

    return pp_vector

@staticmethod
def get_non_degenerated(Ivals, tolerance):
    """
    Get the index of the non-degenerate eigenvalue from the array of eigenvalues.

    :param Ivals:
        The array of eigenvalues.

    :param tolerance:
        The tolerance parameter on the eigenvalues.

    :return:
        The index of the non-degenerate eigenvalue.
    """

    for i, ev1 in enumerate(Ivals):
        single_deg = 0
        index = 0
        for ev2 in Ivals:
            if not abs(ev1 - ev2) < tolerance:
                single_deg += 1
                index = i
        if single_deg == 2:
            return index

    raise Exception("Non-degenerate not found.")

@staticmethod
def get_cubed_sphere_grid_points(tolerance):
    """
    Generate a cubed-grid points grid on the surface of an unitary sphere.

    :param tolerance: 
        Maximum angle between points (radians).

    :return:
        List of points.
    """

    num_points = int(1.0 / tolerance)

    if num_points < 1:
        return [(1, 0, 0)]

    for i in range(-num_points, num_points+1):
        x = i * tolerance
        for j in range(-num_points, num_points+1):
            y = j * tolerance
            for p in permutations([x, y, 1]):
                norm = np.linalg.norm([x, y, 1])
                yield np.array(p)/norm

@staticmethod
def abs_to_rad(tolerance, coord):
    """
    Converts the tolerance from absolute units to radians for an array of coordinates.

    :param tolerance:
        The tolerance defined in absolute units (e.g. the eigenvalue tolerance)

    :param coord:
        The array of coordinates (e.g. the molecule coordinates in center of mass frame).
    
    :return:
        The equivalent of the tolerance in radians.
    """

    coord = np.array(coord)
    return tolerance / np.clip(np.linalg.norm(coord, axis=1), tolerance, None)

@staticmethod
def angle_between_vector_matrix(vector, coord, tolerance=1e-5):
    """
    Calculates the angles between position vectors in the center of mass frame
    and position vectors from another array. 

    :param vector:
        The array of coordinates to compares with center of mass frame coordinates
        (e.g. the coordinates after an operation).
    
    :param coord:
        The reference coordinates (e.g. the center of mass frame coordinates).
        
    :return:
        An array of angles.
    """

    norm_coor = np.linalg.norm(coord, axis=1)
    norm_op_coor = np.linalg.norm(vector)

    angles = []
    for v, n in zip(np.dot(vector, coord.T), norm_coor*norm_op_coor):
        if n < tolerance:
            angles.append(0)
        else:
            angles.append(np.arccos(np.clip(v/n, -1.0, 1.0)))
    return np.array(angles)

@staticmethod
def radius_diff_in_radiants(vector, coord, tolerance=1e-5):
    """
    Calculates the difference between the radii of the vectors in the coord matrix and another vector.

    :param vector:
        The array to evaluate the differences in position vectors.

    :param coord:
        The reference array to evaluate the differences.

    :return:
        An array of differences.
    """

    norm_coor = np.linalg.norm(coord, axis=1)
    norm_op_coor = np.linalg.norm(vector)

    average_radii = np.clip((norm_coor + norm_op_coor) / 2, tolerance, None)
    return np.abs(norm_coor - norm_op_coor) / average_radii

@staticmethod
def rotation_matrix(axis, angle):
    """
    Build a rotation matrix to rotate of a given angle around a vector.

    :param axis:
        The normalized axis or vector around which the rotation is effectuated.

    :param angle:
        The angle to be rotated in radians.

    :return:
        The rotation matrix.
    """

    norm = np.linalg.norm(axis)
    assert norm > 1e-8
    axis = np.array(axis) / norm  # normalize axis

    cos_term = 1 - np.cos(angle)
    rot_matrix = [[axis[0]**2*cos_term + np.cos(angle),              axis[0]*axis[1]*cos_term - axis[2]*np.sin(angle), axis[0]*axis[2]*cos_term + axis[1]*np.sin(angle)],
                  [axis[1]*axis[0]*cos_term + axis[2]*np.sin(angle), axis[1]**2*cos_term + np.cos(angle),              axis[1]*axis[2]*cos_term - axis[0]*np.sin(angle)],
                  [axis[2]*axis[0]*cos_term - axis[1]*np.sin(angle), axis[1]*axis[2]*cos_term + axis[0]*np.sin(angle), axis[2]**2*cos_term + np.cos(angle)]]

    return np.array(rot_matrix)

@staticmethod
def write_xyz_file(symbols,
                   coordinates_in_angstrom,
                   pointgroup,
                   xyz_filename):
        """
        Writes molecular geometry to xyz file.

        :param symbols:
            The list of atomic labels.

        :param coordinates_in_angstrom:
            An array of coordinates.

        :param pointgroup:
            The point group of symmetry.

        :param xyz_filename:
            The name of the xyz file.
        
        :return:
            An xyz file with molecular geometry.
        """

        natoms = len(symbols)
        xyz = f''

        for a in range(natoms):
            xa, ya, za = coordinates_in_angstrom[a]
            xyz += f'{symbols[a]:<6s} {xa:22.12f} {ya:22.12f} {za:22.12f}\n'

        with open(str(xyz_filename) + ".xyz", 'w') as new_geom:
            new_geom.write(str(natoms)+'\n')
            new_geom.write(pointgroup+'\n')
            new_geom.write(str(xyz))



all_symmetry_elements = {
    "C1": ["E"],
    "Cs": ["E", "sigma"],
    "Ci": ["E", "i"],
    "C2": ["E", "C2"],
    "C3": ["E", "2C3"],
    "C4": ["E", "2C4", "C2"],
    "C5": ["E", "4C5"],
    "C6": ["E", "2C6", "2C3", "C2"],
    "C7": ["E", "6C7"],
    "C8": ["E", "4C8", "2C4", "C2"],
    "D2": ["E", "3C2"],
    "D3": ["E", "2C3", "3C2"],
    "D4": ["E", "2C4", "5C2"],
    "D5": ["E", "4C5", "5C2"],
    "D6": ["E", "2C6", "2C3", "7C2"],
    "D7": ["E", "6C7", "7C2"],
    "D8": ["E", "4C8", "2C4", "9C2"],
    "C2v": ["E", "C2", "2sigma_v"],
    "C3v": ["E", "2C3", "3sigma_v"],
    "C4v": ["E", "2C4", "C2", "2sigma_v", "2sigma_d"],
    "C5v": ["E", "4C5", "5sigma_v"],
    "C6v": ["E", "2C6", "2C3", "C2", "3sigma_v", "3sigma_d"],
    "C7v": ["E", "6C7", "7sigma_v"],
    "C8v": ["E", "4C8", "2C4", "C2", "4sigma_v", "4sigma_d"],
    "D2d": ["E", "3C2", "2S4", "2sigma_d"],
    "D3d": ["E", "2C3", "3C2", "i", "2S6", "3sigma_d"],
    "D4d": ["E", "2C4", "5C2", "4S8", "4sigma_d"],
    "D5d": ["E", "4C5", "5C2", "i", "4S10", "5sigma_d"],
    "D6d": ["E", "2C6", "2C3", "7C2", "4S12", "2S4", "6sigma_d"],
    "D7d": ["E", "6C7", "7C2", "i", "6S14", "7sigma_d"],
    "D8d": ["E", "4C8", "2C4", "9C2", "8S16", "8sigma_d"],
    "C2h": ["E", "C2", "i", "sigma_h"],
    "C3h": ["E", "2C3", "sigma_h", "2S3"],
    "C4h": ["E", "2C4", "C2", "i", "2S4", "sigma_h"],
    "C5h": ["E", "4C5", "4S5", "sigma_h"],
    "C6h": ["E", "2C6", "2C3", "C2", "i", "2S3", "2S6", "sigma_h"],
    "C7h": ["E", "6C5", "6S7", "sigma_h"],
    "C8h": ["E", "4C8", "2C4", "C2", "i", "2S4", "4S8", "sigma_h"],
    "D2h": ["E", "3C2", "i", "sigma_h", "sigma_v", "sigma_d"],
    "D3h": ["E", "2C3", "3C2", "2S3", "sigma_h", "3sigma_v"],
    "D4h": ["E", "2C4", "5C2", "i", "S4", "sigma_h", "2sigma_v", "2sigma_d"],
    "D5h": ["E", "4C5", "5C2", "4S5", "sigma_h", "5sigma_v"],
    "D6h": ["E", "2C6", "2C3", "7C2", "i", "S3", "2S6", "sigma_h", "3sigma_v", "3sigma_d"],
    "D7h": ["E", "6C7", "7C2", "6S7", "sigma_h", "7sigma_v"],
    "D8h": ["E", "4C8", "2C4", "9C2", "i", "2S4", "4S8", "sigma_h", "4sigma_v", "4sigma_d"],
    "S4": ["E", "C2", "2S4"],
    "S6": ["E", "2C3", "i", "2S6"],
    "S8": ["E", "2C4", "C2", "4S8"],
    "S10": ["E", "4C5", "i", "4S10"],
    "T": ["E", "8C3", "3C2"],
    "Th": ["E", "8C3", "3C2", "i", "8S6", "3sigma_h"],
    "Td": ["E", "8C3", "3C2", "6S4", "6sigma_d"],
    "O": ["E", "6C4", "8C3", "9C2"],
    "Oh": ["E", "6C4", "8C3", "9C2", "i", "8S6", "6S4", "3sigma_h", "6sigma_d"],
    "I": ["E", "24C5", "20C3", "15C2"],
    "Ih": ["E", "24C5", "20C3", "15C2", "i", "24S10", "20S6", "15sigma"],
    "Cinfv": ["E", "Cinf"],
    "Dinfh": ["E", "Cinf", "i"]
}

subgroups = {
    "C1": ["C1"],
    "Cs": ["Cs", "C1"],
    "Ci": ["C1"],
    "C2": ["C1"],
    "C3": ["C1"],
    "C4": ["C2", "C1"],
    "C5": ["C1"],
    "C6": ["C2", "C1"],
    "C7": ["C1"],
    "C8": ["C4", "C2", "C1"],
    "D2": ["C2", "C1"],
    "D3": ["C2", "C1"],
    "D4": ["D2", "C2", "C1"],
    "D5": ["C2", "C1"],
    "D6": ["D2", "C2", "C1"],
    "D7": ["C2", "C1"],
    "D8": ["C2", "C1"],
    "C2v": ["C2", "Cs", "C1"],
    "C3v": ["Cs", "C1"],
    "C4v": ["C2v", "C2", "Cs", "C1"],
    "C5v": ["Cs", "C1"],
    "C6v": ["C2v", "C2", "Cs", "C1"],
    "C7v": ["Cs", "C1"],
    "C8v": ["C2v", "C2", "Cs", "C1"],
    "D2d": ["C2v", "D2", "C2", "Cs", "C1"],
    "D3d": ["C2h", "C2", "Cs", "Ci", "C1"],
    "D4d": ["C2v", "D2", "C2", "Cs", "C1"],
    "D5d": ["C2h", "C2", "Cs", "Ci", "C1"],
    "D6d": ["C2v", "D2", "C2", "Cs", "C1"],
    "D7d": ["C2h", "C2", "Cs", "Ci", "C1"],
    "D8d": ["C2v", "D2", "C2", "Cs", "C1"],
    "C2h": ["C2", "Cs", "Ci", "C1"],
    "C3h": ["Cs", "C1"],
    "C4h": ["C2h", "C2", "Cs", "Ci", "C1"],
    "C5h": ["Cs", "C1"],
    "C6h": ["C2h", "C2", "Cs", "Ci", "C1"], 
    "C7h": ["Cs", "C1"],
    "C8h": ["C2h", "C2", "Cs", "Ci", "C1"],
    "D2h": ["D2", "C2h", "C2v", "C2", "Cs", "Ci", "C1"],
    "D3h": ["C2v", "C2", "Cs", "C1"],
    "D4h": ["D2h", "D2", "C2h", "C2v", "C2", "Cs", "Ci", "C1"],
    "D5h": ["C2v", "C2", "Cs", "C1"],
    "D6h": ["D2h", "D2", "C2h", "C2", "Cs", "Ci", "C1"],
    "D7h": ["C2v", "C2", "Cs", "C1"],
    "D8h": ["D2h", "D2", "C2h", "C2v", "C2", "Cs", "Ci", "C1"],
    "S4": ["C2", "C1"],
    "S6": ["C3", "Ci", "C1"],
    "S8": ["C2", "C1"],
    "S10": ["Ci", "C1"],
    "T": ["D2", "C3", "C2", "C1"],
    "Th": ["D2h", "D2", "C2h", "C2v", "C2", "Cs", "Ci", "C1"],
    "Td": ["D2", "C2v", "C2", "Cs", "C1"],
    "O": ["D2", "C2", "C1"],
    "Oh": ["D2h", "D2", "C2h", "C2v", "C2", "Cs", "Ci", "C1"],
    "I": ["D2", "C2", "C1"],
    "Ih": ["D2h", "D2", "C2h", "C2v", "C2", "Cs", "Ci", "C1"],
    "Cinfv": ["C2v", "C2", "Cs", "C1"],
    "Dinfh": ["D2h", "D2", "C2v", "C2h", "C2v", "C2", "Ci", "Cs", "C1"]
}