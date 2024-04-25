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
        - inequivalent_atoms: a dictionary containing the operations names and the associated inequivalent atoms.
        - primary_axis: axis defined as z cartesian axis for reorientation in point group convention.
        - secondary_axis: axis perpendicular to primary axis to build conventional cartesian frame for reorientation.
        - molecule_type: the type of molecule on a detection tree (e.g., linear, cyclic, etc.).

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
        self._inequivalent_atoms = {'operations': [], 'ineq_atoms': []}
        self._primary_axis = [1., 0., 0.]
        self._secondary_axis = [0., 1., 0.]
        self._molecule_type = ''
        
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
                - The list of available Abelian groups for symmetrization.
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
        self._reoriented_coordinates = np.zeros((self._natoms,3))

        # Handle case of isolated atom
        if self._natoms == 1:
            self._schoenflies_symbol = "O(3)"
            self._molecule_type = "isolated_atom"
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
        if self._max_order <= 2 and self._schoenflies_symbol not in ["O(3)", "D2d"]:
            symmetry_analysis["Elements_found"] = self._inequivalent_atoms["operations"]

        symmetry_analysis["Point_group"] = self._schoenflies_symbol

        if self._schoenflies_symbol in all_symmetry_elements:
            symmetry_elements_list = all_symmetry_elements[self._schoenflies_symbol]
            symmetry_analysis["Expected_elements"] = symmetry_elements_list

        if self._schoenflies_symbol in groups_for_symmetrization:
            Groups_for_symm_list = groups_for_symmetrization[self._schoenflies_symbol]
            symmetry_analysis["Groups_for_symm"] = Groups_for_symm_list

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

        if "Groups_for_symm" in results_dict:
            print("Available Abelian groups for symmetrization: {}".format(results_dict["Groups_for_symm"]))
        else:
            print("No available subgroups for symmetrization.")

    def print_tolerance_keywords(self):
        """
        :return:
            The available tolerance keywords.
        """

        tolerance_keywords = ["loose", "tight", "very tight"]
        return tolerance_keywords

    def symmetrize_pointgroup(self, symmetry_data, pointgoup_to_symmetrize=None):
        """
        Symmetrize in chosen available Abelian group with real representation
        and reorient molecules. Only reorient if more than 60 atoms.

        :param pointgoup_to_symmetrize:
            The chosen point group in which themolecule will be symmetrized. 
            Default: the detected point group.

        !!!Under construction!!!
        """

        symmetrized_data = {}

        # Initializes
        if pointgoup_to_symmetrize == None:
            pointgoup_to_symmetrize = self._schoenflies_symbol

        # Check that the chosen point group is available for symmetrization
        if pointgoup_to_symmetrize not in symmetry_data["Groups_for_symm"]:
            if pointgoup_to_symmetrize == self._schoenflies_symbol:
               # Store the geometry reoriented in the conventional point group orientation
               symmetrized_data["new_pointgroup"] = "Initial geometry reoriented"  
               self._conventional_orientation()
               symmetrized_data["reoriented_geom_in_angstrom"] = self._reoriented_coordinates * bohr_in_angstrom()

            else:
                raise KeyError("Point group not available for symmetrization.")
        
        else:
            # Temporary
            symmetrized_data["new_pointgroup"] = "Initial geometry reoriented"
            self._conventional_orientation()
            symmetrized_data["reoriented_geom_in_angstrom"] = self._reoriented_coordinates * bohr_in_angstrom()
            
            # if self._natoms <= 60:
            #     symmetrized_coords = self._symmetrize_molecule(pointgoup_to_symmetrize)
            # else:
            #     self._conventional_orientation()
            #     symmetrized_data["reoriented_geom_in_angstrom"] = self._reoriented_coordinates * bohr_in_angstrom()

        # Temporary
        symmetrized_data["inequiv_atoms_by_op"] = self._inequivalent_atoms

        return symmetrized_data

    def print_symmetrized_molecule(self, results_symmetrization, xyz_filename=None):
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

        if "symmetrized_coord" in results_symmetrization:
            coord = results_symmetrization["symmetrized_coord"]
        else:
            coord = results_symmetrization["reoriented_geom_in_angstrom"]

        pointgroup = results_symmetrization["new_pointgroup"]

        if xyz_filename:
            write_xyz_file(self._symbols,coord,pointgroup,xyz_filename)
        else:
            print(coord)

        # Temporary
        print(results_symmetrization["inequiv_atoms_by_op"])

    def _linear(self):
        """
        Handle linear molecules.

        :return:
            The Schoenflies symbol, list of detected elements, and reoriented coordinates.
        """

        # Set type of molecule for reorientation
        self._molecule_type = "linear"

        # Set orientation
        idx = np.argmin(self._Ivals)
        main_axis = self._Ivecs[idx]
        p_axis = get_perpendicular(main_axis)
        self._set_orientation(main_axis, p_axis)

        # Check for inversion center at center of mass
        if self._check_op(Inversion(), "i"):
            self._schoenflies_symbol = 'Dinfh'
        else:
            self._schoenflies_symbol = 'Cinfv'
        
        # Set primary axis as z cartesian axis for conventional reorientation
        self._primary_axis = [0.,0.,1.]

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
            if self._check_op(c2, "C2"):
                n_axis_c2 += 1
                main_axis = axis
        
        self._max_order = 2

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

            # Check for C2 axis along main axis
            # Imperative for symmetrization in abelian with real characters subgroups
            self._check_op(Rotation(main_axis, order=2), "C2")

            # Get the perpendicualar axis to main axis and check for C2 rotation axis 
            # along p_axis by rotating p_axis along the main axis
            p_axis = get_perpendicular(main_axis)
            for angle in np.arange(0, np.pi, 0.1* np.pi / self._max_order):
                axis = np.dot(p_axis, rotation_matrix(main_axis, angle))
                c2 = Rotation(axis, order=2)
                if self._check_op(c2, "C2"):
                    self._dihedral(main_axis)
                    return
                        
            self._cyclic(main_axis)

    def _spherical(self):
        """
        Handle spherical groups (I, O, T).
        """

        # Set type of molecule for reorientation
        self._molecule_type = "spherical"

        for axis in get_cubed_sphere_grid_points(self._tolerance_ang):
            c5 = Rotation(axis, order=5)
            c4 = Rotation(axis, order=4)
            c3 = Rotation(axis, order=3)

            # Check for C5 axis
            if self._check_op(c5, "C5", tol_factor=1.0):
                self._schoenflies_symbol = "I"
                main_axis = axis
                self._max_order = 5
                break
            # Check for C4 axis
            elif self._check_op(c4, "C4", tol_factor=1.5):
                self._schoenflies_symbol = "O"
                main_axis = axis
                self._max_order = 4
                break
            # Check for C3 axis
            elif self._check_op(c3, "C3", tol_factor=1.5):
                self._schoenflies_symbol = "T"
                main_axis = axis
                self._max_order = 3

        p_axis_base = get_perpendicular(main_axis)

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

                    if self._check_op(c5, "C5"):
                        t_axis = np.dot(main_axis, rotation_matrix(p_axis_base, np.pi/2).T)
                        return np.dot(t_axis, rot_matrix.T)

            # Set orientation
            p_axis = _determine_orientation_I(main_axis)
            self._set_orientation(main_axis, p_axis)

            # Check for inverison center at center of mass
            if self._check_op(Inversion(), "i"):
                self._schoenflies_symbol += 'h'

                # Check for one reflexion plane containing a C5 axis as it is a generator for some subgroups
                for angle in np.arange(0, np.pi, 0.1*np.pi / self._max_order + self._tolerance_ang):
                    axis = np.dot(p_axis, rotation_matrix(main_axis, angle))
                    if self._check_op(Reflection(axis), "sigma"):
                        break

            # Check for C2 axes perpendicular to main axis
            # Useful for reorienting and symmetrizing in D subgroups
            for angle in np.arange(0, np.pi/2 + self._tolerance_ang, self._tolerance_ang):
                axis = np.dot(p_axis, rotation_matrix(main_axis, angle))
                self._check_op(Rotation(axis, order=2), "C2")

        # O or Oh
        if self._schoenflies_symbol == 'O':
            
            def _determine_orientation_O(main_axis):
                """
                Determine the orientation of the molecule within the octahedron symmetry group. 

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

                    if self._check_op(c4, "C4"):
                        return axis

            # Set orientation
            p_axis = _determine_orientation_O(main_axis)
            self._set_orientation(main_axis, p_axis)

            # Check for inverison center at center of mass
            if self._check_op(Inversion(), "i"):
                self._schoenflies_symbol += 'h'

                # Check for one reflexion plane containing a C4 axis as it is a generator for some subgroups
                for angle in np.arange(0, np.pi, 0.1*np.pi / self._max_order + self._tolerance_ang):
                    axis = np.dot(p_axis, rotation_matrix(main_axis, angle))
                    if self._check_op(Reflection(axis), "sigma"):
                        break

            # Check for C2 axis perpendicular to main axis
            # Useful for reorienting and symmetrizing in D subgroups
            for angle in np.arange(0, np.pi/2 + self._tolerance_ang, 0.2*np.pi / self._max_order + self._tolerance_ang):
                axis = np.dot(p_axis, rotation_matrix(main_axis, angle))
                self._check_op(Rotation(axis, order=2), "C2")

        # T or Td, Th
        if self._schoenflies_symbol == 'T':

            def _determine_orientation_T(main_axis):
                """
                Determine the orientation of the molecule within the tetrahedron symmetry group. 

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

                    if self._check_op(c3, "C3"):
                        t_axis = np.dot(main_axis, rotation_matrix(p_axis_base, np.pi/2).T)
                        return np.dot(t_axis, rot_matrix.T)

            p_axis = _determine_orientation_T(main_axis)
            # Set orientation
            self._set_orientation(main_axis, p_axis)

            # Check for C2 axis perpendicular to main axis
            # Useful for reorienting and symmetrizing in D subgroups
            for angle in np.arange(0, np.pi + self._tolerance_ang, 0.05*np.pi / self._max_order + self._tolerance_ang):
                axis = np.dot(main_axis, rotation_matrix(p_axis, angle))
                if self._check_op(Rotation(axis, order=2), "C2"):
                    temporary_main_axis = axis
                    
                    # Check for another C2 axis perpendicular to the first one
                    # (Collect the generators of D2x subgroups)
                    temporary_p_axis = get_perpendicular(temporary_main_axis)
                    for angle in np.arange(0, np.pi/2 + self._tolerance_ang, 0.1*np.pi / self._max_order + self._tolerance_ang):
                        axis = np.dot(temporary_p_axis, rotation_matrix(temporary_main_axis, angle))
                        if self._check_op(Rotation(axis, order=2), "C2"):
                            break
                    # Check for one reflexion plane containing a C2 axis as it is a generator for some subgroups
                    self._check_op(Reflection(temporary_main_axis), "sigma")

            # Check for inverison center at center of mass
            if self._check_op(Inversion(), "i"):
                self._schoenflies_symbol += 'h'
                return
            
            # Check for any reflexion plane 
            if self._check_op(Reflection([0, 0, 1]), "sigma"):
                self._schoenflies_symbol += 'd'
                return

    def _no_rot_axis(self):
        """
        Detect point group of molecule with no rotation axis.

        :return:
            The Schoenflies symbol and list of detected elements.
        """

        # Set type of molecule for reorientation
        self._molecule_type = "asym_top"

        for i, vector in enumerate(np.identity(3)):
            if self._check_op(Reflection(vector), "sigma"):
                self._schoenflies_symbol = 'Cs'
                self._primary_axis = vector
                self._secondary_axis = get_perpendicular(vector)
                self._set_orientation(vector, self._secondary_axis)
                break
            else:
                if self._check_op(Inversion(), "i"):
                    self._schoenflies_symbol = 'Ci'
                    break
                else:
                    self._schoenflies_symbol = 'C1'

    def _cyclic(self, main_axis):
        """
        Detect point group of cylic group molecules.

        :param main_axis:
            The main axis obtained with _set_orientation.

        :return:
            The Schoenflies symbol and list of detected elements.
        """

        # Set type of molecule for reorientation
        self._molecule_type = "cyclic"
        self._primary_axis = main_axis

        self._schoenflies_symbol = "C{}".format(self._max_order)
        
        # Check for reflexion planes perpenducular to the main axis
        if self._check_op(Reflection(main_axis), "sigma_h"):
            self._schoenflies_symbol += 'h'

        # Check for reflexion planes containing the main axis
        v_symbol = set()
        self._secondary_axis = get_perpendicular(main_axis)
        for angle in np.arange(0, np.pi + self._tolerance_ang, 0.5*np.pi / self._max_order + self._tolerance_ang):
            axis = np.dot(self._secondary_axis, rotation_matrix(main_axis, angle))
            if self._check_op(Reflection(axis), "sigma_v"):
                v_symbol.add('v')

        self._schoenflies_symbol += ''.join(v_symbol)

        # Check for inversion center at center of mass
        self._check_op(Inversion(), "i")

        # Check for improper rotational axis along the main axis.
        if self._check_op(ImproperRotation(main_axis, order=2*self._max_order), "S{}".format(2 * self._max_order)):
            self._schoenflies_symbol = "S{}".format(2 * self._max_order)
            
    def _dihedral(self, main_axis):
        """
        Detect point group of dihedral group molecules.

        :param main_axis:
            The main axis obtained with _set_orientation.

        :return:
            The Schoenflies symbol and list of detected elements.
        """

        # Set type of molecule for reorientation
        self._molecule_type = "cyclic"
        self._primary_axis = main_axis

        # Determine perpendicular axis to main axis
        self._secondary_axis = get_perpendicular(main_axis)

        if self._max_order == 1:
            # D1 is equivalent to C2
            self._schoenflies_symbol = "C2"
        else:
            self._schoenflies_symbol = "D{}".format(self._max_order)

        # Check for inversion center at center of mass
        self._check_op(Inversion(), "i")

        # Check for reflexion planes perpenducular to the main axis
        h_symbols = False
        if self._check_op(Reflection(main_axis), "sigma_h") and not h_symbols:
            self._schoenflies_symbol += 'h'
            h_symbols = True
        
        # Check for reflexion planes containing the main axis
        d_symbol = False
        for angle in np.arange(0, np.pi, 0.5*np.pi / self._max_order):
            axis = np.dot(self._secondary_axis, rotation_matrix(main_axis, angle))
            if self._check_op(Reflection(axis), "sigma") and not h_symbols and not d_symbol:
                self._schoenflies_symbol += 'd'
                d_symbol = True
    
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
                if self._check_op(Cn, "C{}".format(i)):
                    return i
            return 1

    def _check_op(self, operation, element_string, tol_factor=1.0):
        """
        Check if the given symmetry operation exists in the point group of the molecule.

        :param operation:
            The matrix representative of the symmetry operation.

        :param tol_factor:
            A factor to scale the tolerance value.

        :return:
            True if the symmetry operation exists in the point group, False otherwise.
        """

        # Get representative of the operation
        sym_matrix = operation.get_matrix()

        # Define absolte tolerance from the eigenvalue tolerance
        error_abs_rad = abs_to_rad(self._tolerance_eig, coord=self._cent_coord)

        # Get COM frame coordinates after the operation
        op_coordinates = np.matmul(self._cent_coord, sym_matrix)

        # Initialize objects to obtain inequivalent atoms
        mapping = []
        inequivalent_atoms_by_operation = set()

        # Check if operation exists
        for idx, op_coord in enumerate(op_coordinates):
            # Calculate the differences in radii and angles
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
            
            # Save the original atom index, the manipulated atom index, and their associated original coordinates
            manipulated_atom_idx = np.argmin(np.linalg.norm(self._cent_coord - op_coord, axis=1))
            mapping.append((idx, manipulated_atom_idx))

            # Check if the original atom is already in the list of inequivalent atoms
            if idx == manipulated_atom_idx:
                inequivalent_atoms_by_operation.add(idx)
            elif idx != manipulated_atom_idx:
                for mapping_tuple in mapping:
                    if mapping_tuple[0] not in [indice_tuple[1] for indice_tuple in mapping]:
                        inequivalent_atoms_by_operation.add(mapping_tuple[0])

        # Save string and inequivalent atom indices for detected operations
        self._inequivalent_atoms["operations"].append(element_string)
        self._inequivalent_atoms["ineq_atoms"].append(inequivalent_atoms_by_operation)

        return True

    def _set_orientation(self, main_axis, p_axis):
        """
        Set molecular orientation along main_axis (x) and p_axis (y).

        :param main_axis:
            Principal orientation axis (must be unitary).

        :param p_axis:
            Secondary axis perpendicular to principal (must be unitary).
        """

        assert np.linalg.norm(main_axis) > 1e-1
        assert np.linalg.norm(p_axis) > 1e-1

        orientation = np.array([main_axis, p_axis, np.cross(main_axis, p_axis)])
        self._cent_coord = np.dot(self._cent_coord, orientation.T)
        self._ref_orientation = np.dot(self._ref_orientation, orientation.T)
    
    def _conventional_orientation(self):
        """
        Set molecular orientation in point group convention:
            - a right-handed cartesian,
            - origin at center of mass,
            - highest rotational axis asz cartesian axis,
            - if z in along a sigma plane, x is chosen as perpendicular to the plane,
            - if z is perpendicular to a sigma plane, x and y are in the plane,
            - if no rotational axis but a sigma plane, z is set along the plane,
            - if no rotational axis nor sigma plane, x and y are defined along the two principal inertia axes
              associated with the two smallest principal moments.
            
        :return:
            A list of arrays with the cartesian coordinates.
        """

        if self._molecule_type not in ["isolated_atom", "linear", "asym_top", "cyclic", "dihedral", "spherical"]:
            raise KeyError("Molecule type not available.")
        
        if self._molecule_type =="isolated_atom":
            return
        elif self._molecule_type == "linear":
            self._set_orientation(self._primary_axis, self._secondary_axis)
            self._reoriented_coordinates = self._cent_coord
        elif self._molecule_type == "asym_top":
            self._reoriented_coordinates = self._cent_coord
        else:
            orientation = np.array([self._secondary_axis, np.cross(self._primary_axis, self._secondary_axis), self._primary_axis])
            self._reoriented_coordinates = np.dot(self._cent_coord, orientation.T)



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
    # Define all the expected symmetry elements for each point group
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

groups_for_symmetrization = {
    # Define all the Abelian group available for symmetrization for each point group
    "C1": ["C1"],
    "Cs": ["Cs", "C1"],
    "Ci": ["Ci", "C1"],
    "C2": ["C2", "C1"],
    "C3": ["C1"],
    "C4": ["C2", "C1"],
    "C5": ["C1"],
    "C6": ["C2", "C1"],
    "C7": ["C1"],
    "C8": ["C2", "C1"],
    "D2": ["D2", "C2", "C1"],
    "D3": ["C2", "C1"],
    "D4": ["D2", "C2", "C1"],
    "D5": ["C2", "C1"],
    "D6": ["D2", "C2", "C1"],
    "D7": ["C2", "C1"],
    "D8": ["D2", "C2", "C1"],
    "C2v": ["C2v", "C2", "Cs", "C1"],
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
    "C2h": ["C2h", "C2", "Cs", "Ci", "C1"],
    "C3h": ["Cs", "C1"],
    "C4h": ["C2h", "C2", "Cs", "Ci", "C1"],
    "C5h": ["Cs", "C1"],
    "C6h": ["C2h", "C2", "Cs", "Ci", "C1"], 
    "C7h": ["Cs", "C1"],
    "C8h": ["C2h", "C2", "Cs", "Ci", "C1"],
    "D2h": ["D2h", "D2", "C2h", "C2v", "C2", "Cs", "Ci", "C1"],
    "D3h": ["C2v", "C2", "Cs", "C1"],
    "D4h": ["D2h", "D2", "C2h", "C2v", "C2", "Cs", "Ci", "C1"],
    "D5h": ["C2v", "C2", "Cs", "C1"],
    "D6h": ["D2h", "D2", "C2h", "C2", "Cs", "Ci", "C1"],
    "D7h": ["C2v", "C2", "Cs", "C1"],
    "D8h": ["D2h", "D2", "C2h", "C2v", "C2", "Cs", "Ci", "C1"],
    "S4": ["C2", "C1"],
    "S6": ["Ci", "C1"],
    "S8": ["C2", "C1"],
    "S10": ["Ci", "C1"],
    "T": ["D2", "C2", "C1"],
    "Th": ["D2h", "D2", "C2h", "C2v", "C2", "Cs", "Ci", "C1"],
    "Td": ["D2", "C2v", "C2", "Cs", "C1"],
    "O": ["D2", "C2", "C1"],
    "Oh": ["D2h", "D2", "C2h", "C2v", "C2", "Cs", "Ci", "C1"],
    "I": ["D2", "C2", "C1"],
    "Ih": ["D2h", "D2", "C2h", "C2v", "C2", "Cs", "Ci", "C1"],
    "Cinfv": ["C2v", "C2", "Cs", "C1"],
    "Dinfh": ["D2h", "D2", "C2h", "C2v", "C2", "Ci", "Cs", "C1"]
}

generators = {
    # C2_p means a C2 axis perpendicular to the first C2 axis
    "Cs": ["sigma"],
    "Ci": ["i"],
    "C2": ["C2"],
    "D2": ["C2", "C2_p"],
    "C2v": ["C2", "sigma_v"],
    "C2h": ["C2", "sigma_h"],
    "D2h": ["C2", "C2_p", "sigma_h"],
}