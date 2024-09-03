import numpy as np
import networkx as nx
import re

from .molecule import Molecule
from .symmetryoperations import (Inversion, Rotation, Reflection,
                                 ImproperRotation)
from .symmetryoperations import rotation_matrix
from .errorhandler import assert_msg_critical, safe_arccos


class SymmetryAnalyzer:
    """
    A class to determine the point group of nearly symmetric molecules within selected
    tolerance tresholds (in Schoenflies convention) and then to symmetrize it perfectly
    in the detected point group or in one od its subgroups.

    1. The coordinates are reexpressed in the center of mass frame.
    2. The principal axis and moment are obtained by computing the eigenvalues and
       eigenvectors of the inertia tensor.
        - Linear molecule have one vanishing eigenvalue.
        - Asymmetric top molecules have nondegenerate eigenvalues.
        - Symmetric molecules have doubly degenerate eigenvalues.
        - Spherical molecules have triply degenerate eigenvalues.
    3. Each category is then treated following a typical detection tree.

    Instance variables
        - schoenflies_symbol: the pointgroup schoenflies symbol.
        - max_order: the maximum proper rotational axis order.
        - primary_axis: axis defined as z cartesian axis for reorientation in
          point group convention.
        - secondary_axis: axis perpendicular to primary axis to build
          conventional cartesian frame for reorientation.
        - molecule_type: the type of molecule on a detection tree (e.g.,
          linear, cyclic, etc.).
    """

    def __init__(self):
        """
        Initializes the SymmetryAnalyzer instance.
        """

        # Initializes variables
        self._schoenflies_symbol = ''
        self._max_order = 1
        self._primary_axis = [1., 0., 0.]
        self._secondary_axis = [0., 1., 0.]
        self._molecule_type = ''

        # Define all the expected symmetry elements for each point group
        self._all_symmetry_elements = {
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
            "D4h": [
                "E", "2C4", "5C2", "i", "S4", "sigma_h", "2sigma_v", "2sigma_d"
            ],
            "D5h": ["E", "4C5", "5C2", "4S5", "sigma_h", "5sigma_v"],
            "D6h": [
                "E", "2C6", "2C3", "7C2", "i", "S3", "2S6", "sigma_h",
                "3sigma_v", "3sigma_d"
            ],
            "D7h": ["E", "6C7", "7C2", "6S7", "sigma_h", "7sigma_v"],
            "D8h": [
                "E", "4C8", "2C4", "9C2", "i", "2S4", "4S8", "sigma_h",
                "4sigma_v", "4sigma_d"
            ],
            "S4": ["E", "C2", "2S4"],
            "S6": ["E", "2C3", "i", "2S6"],
            "S8": ["E", "2C4", "C2", "4S8"],
            "S10": ["E", "4C5", "i", "4S10"],
            "T": ["E", "8C3", "3C2"],
            "Th": ["E", "8C3", "3C2", "i", "8S6", "3sigma_h"],
            "Td": ["E", "8C3", "3C2", "6S4", "6sigma_d"],
            "O": ["E", "6C4", "8C3", "9C2"],
            "Oh": [
                "E", "6C4", "8C3", "9C2", "i", "8S6", "6S4", "3sigma_h",
                "6sigma_d"
            ],
            "I": ["E", "24C5", "20C3", "15C2"],
            "Ih": [
                "E", "24C5", "20C3", "15C2", "i", "24S10", "20S6", "15sigma"
            ],
            "Cinfv": ["E", "Cinf"],
            "Dinfh": ["E", "Cinf", "i"],
        }

    def identify_pointgroup(self, molecule, tolerance='tight'):
        """
        Analyze the symmetry of a given molecule based on its nuclear framework.

        :param molecule:
            A VeloxChem molecule object.
        :param tolerance:
            A tolerance threshold string available in tolerance_keywords.
            Default: tight

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

        if tolerance == 'very loose':
            self._tolerance_eig = 0.8
            self._tolerance_ang = np.radians(5.0)
        elif tolerance == 'loose':
            self._tolerance_eig = 0.4
            self._tolerance_ang = np.radians(4.0)
        elif tolerance == 'tight':
            self._tolerance_eig = 0.1
            self._tolerance_ang = np.radians(3.0)
        elif tolerance == 'very tight':
            self._tolerance_eig = 0.002
            self._tolerance_ang = np.radians(2.0)
        else:
            raise KeyError(
                "SymmetryAnalyzer: Tolerance criterion not available.")

        # Read and express geometry in center of mass (COM) frame
        # Geom and COM in bohr because moments of inertia are defined in bohr
        # in Molecule module
        coordinates = molecule.get_coordinates_in_bohr()
        self._symbols = molecule.get_labels()
        self._natoms = molecule.number_of_atoms()
        center_of_mass = molecule.center_of_mass_in_bohr()
        self._centered_coords = coordinates - center_of_mass

        # Get the principal momemts amd axes of inertia
        Ivals, Ivecs = molecule.moments_of_inertia(principal_axes=True)
        self._Ivals = Ivals
        self._Ivecs = Ivecs

        # Initializes ensembles of results
        symmetry_analysis = {}

        # Handle case of isolated atom
        if self._natoms == 1:
            self._schoenflies_symbol = "O(3)"
            self._molecule_type = "isolated_atom"
            symmetry_analysis[
                "degeneracy"] = "Input structure is an isolated atom."

        else:
            # Get the degeneracy of the eigenvalues of the inertia tensor
            eig_degeneracy = self._get_degeneracy(self._Ivals,
                                                  self._tolerance_eig)

            # Linear groups
            if np.min(abs(self._Ivals)) < self._tolerance_eig:
                self._handle_linear()
                symmetry_analysis["degeneracy"] = "Molecule is linear."

            # Asymmetric group
            elif eig_degeneracy == 1:
                self._handle_asymmetric()
                symmetry_analysis["degeneracy"] = (
                    "Principal moments of inertia: Nondegenerate")

            # Symmetric group
            elif eig_degeneracy == 2:
                self._handle_symmetric()
                symmetry_analysis["degeneracy"] = (
                    "Principal moments of inertia: Doubly degenerate")

            # Spherical group
            elif eig_degeneracy == 3:
                self._handle_spherical()
                symmetry_analysis["degeneracy"] = (
                    "Principal moments of inertia: Triply degenerate")

            else:
                assert_msg_critical(False,
                                    'SymmetryAnalyzer: Invalid eig_degeneracy')

        # Collect the results
        symmetry_analysis["point_group"] = self._schoenflies_symbol

        if self._schoenflies_symbol in self._all_symmetry_elements:
            symmetry_elements_list = self._all_symmetry_elements[
                self._schoenflies_symbol]
            symmetry_analysis[
                "expected_symmetry_elements"] = symmetry_elements_list

        return symmetry_analysis

    def print_symmetry_results(self, results_dict, symmetry_info=True):
        """
        Build a small output with the results of the symmetry analysis.

        :param results_dict:
            The dictionary containing the different results from the
            pointgroup_identify function.
        :param symmetry_info:
            The flag for printing more information about symmetry.
        """

        if symmetry_info:
            print(results_dict['degeneracy'])

        print('Point group: {}'.format(results_dict['point_group']))

        if symmetry_info:
            if 'expected_symmetry_elements' in results_dict:
                print('Expected symmetry elements: ' +
                      ', '.join(results_dict["expected_symmetry_elements"]))

    @staticmethod
    def _reorder_symmetry_elements(symmetry_elements_list):
        """
        Reorder a list of symmetry elements based on order and type of symmetry
        element. ImproperRotation and Rotation are prioritized.

        :param symmetry_elements_list:
            The list of symmetry elements. Example: ['3C2', '2S4', 'E']

        :return:
            The reordered list. Example: ['2S4', '3C2', 'E']
        """

        def get_op_name(sym_elem):
            m = re.search(r'(\d*)(\D+.*)$', sym_elem)
            return m.group(2)

        def is_improper_rotation(sym_op_name):
            return (sym_op_name.startswith('S') and sym_op_name[1:].isdigit())

        def is_rotation(sym_op_name):
            return (sym_op_name.startswith('C') and sym_op_name[1:].isdigit())

        reordered_symmetry_elements = []

        # find maximum order
        max_op_order = 0
        for sym_elem in symmetry_elements_list:
            sym_op_name = get_op_name(sym_elem)
            if is_improper_rotation(sym_op_name) or is_rotation(sym_op_name):
                op_order = int(sym_op_name[1:])
                if max_op_order < op_order:
                    max_op_order = op_order

        # prioritize ImproperRotation and Rotation
        # also prioritize higher order elements
        for op_order in range(max_op_order, 0, -1):

            for sym_elem in symmetry_elements_list:
                sym_op_name = get_op_name(sym_elem)
                if (is_improper_rotation(sym_op_name) and
                        int(sym_op_name[1:]) == op_order):
                    reordered_symmetry_elements.append(sym_elem)

            for sym_elem in symmetry_elements_list:
                sym_op_name = get_op_name(sym_elem)
                if (is_rotation(sym_op_name) and
                        int(sym_op_name[1:]) == op_order):
                    reordered_symmetry_elements.append(sym_elem)

        # add the remaining elements
        for sym_elem in symmetry_elements_list:
            sym_op_name = get_op_name(sym_elem)
            if not (is_improper_rotation(sym_op_name) or
                    is_rotation(sym_op_name)):
                reordered_symmetry_elements.append(sym_elem)

        return reordered_symmetry_elements

    def symmetrize_pointgroup(self, symmetry_data, point_group=None):
        """
        Symmetrize and reorient molecule.

        :param symmetry_data:
            The dictionary containing the different results from the
            pointgroup_symmetrize function.
        :param point_group:
            The chosen point group in which the molecule will be symmetrized.

        :return:
            The symmetrized molecule.
        """

        if point_group is None:
            point_group = symmetry_data['point_group']

        centered_mol = Molecule(self._symbols, self._centered_coords, 'au')

        if self._natoms == 1:
            return centered_mol

        # go through symmetry elements

        symmetry_elements_list = self._reorder_symmetry_elements(
            self._all_symmetry_elements[point_group])

        symmetry_mapping = set()
        symmetry_operations = []

        mol_grid_points = self._get_mol_grid_points()

        for sym_elem in symmetry_elements_list:
            if sym_elem == 'E':
                continue

            m = re.search(r'(\d*)(\D+.*)$', sym_elem)
            # n_sym_ops = 1
            # if m.group(1):
            #     n_sym_ops = int(m.group(1))
            sym_op_name = m.group(2)

            # find rotation and improper rotation axes
            if (sym_op_name.startswith('C') or
                    sym_op_name.startswith('S')) and sym_op_name[1:].isdigit():
                if sym_op_name.startswith('C'):
                    m = re.search(r'^C(\d+)$', sym_op_name)
                elif sym_op_name.startswith('S'):
                    m = re.search(r'^S(\d+)$', sym_op_name)
                order = int(m.group(1))
                assert_msg_critical(
                    order > 1, 'SymmetryAnalyzer.symmetrize_pointgroup: ' +
                    'Rotation order must be greater than 1')

                for axis in mol_grid_points:
                    for power in range(1, order):
                        if sym_op_name.startswith('C'):
                            sym_op = Rotation(axis, order=order, power=power)
                        elif sym_op_name.startswith('S'):
                            sym_op = ImproperRotation(axis,
                                                      order=order,
                                                      power=power)
                        sym_op_exists = self._check_symmetry_operation(
                            sym_op, sym_op_name, mapping=True)
                        if power == 1 and not sym_op_exists:
                            break
                        for pair in self._mapping:
                            symmetry_mapping.add(pair)
                        symmetry_operations.append(sym_op)

            # find reflection planes
            elif sym_op_name.startswith('sigma') or sym_op_name == 'Cs':
                rotation_op_axes = [
                    np.array(sym_op._axis)
                    for sym_op in symmetry_operations
                    if isinstance(sym_op, (Rotation, ImproperRotation))
                ]

                for axis in mol_grid_points + rotation_op_axes:
                    sym_op = Reflection(axis)
                    sym_op_exists = self._check_symmetry_operation(sym_op,
                                                                   sym_op_name,
                                                                   mapping=True)
                    if sym_op_exists:
                        for pair in self._mapping:
                            symmetry_mapping.add(pair)
                        symmetry_operations.append(sym_op)

            # find inversion
            elif sym_op_name == 'Ci':
                sym_op = Inversion()
                sym_op_exists = self._check_symmetry_operation(sym_op,
                                                               sym_op_name,
                                                               mapping=True)
                if sym_op_exists:
                    for pair in self._mapping:
                        symmetry_mapping.add(pair)
                    symmetry_operations.append(sym_op)

        # average rotation axes

        symmetry_axes = []
        symmetry_indices = []

        for sym_idx, sym_op in enumerate(symmetry_operations):
            if isinstance(sym_op, (Rotation, ImproperRotation)):
                symmetry_axes.append(np.array(sym_op._axis))
                symmetry_indices.append(sym_idx)
        symmetry_axes = np.array(symmetry_axes)

        for i in range(symmetry_axes.shape[0]):
            sum_axis = np.array(symmetry_axes[i])
            count_axis = 1

            for sym_op in symmetry_operations:
                if isinstance(sym_op, (Rotation, ImproperRotation)):
                    op_axes = np.matmul(symmetry_axes, sym_op.get_matrix())

                    for j in range(op_axes.shape[0]):
                        dot = np.dot(symmetry_axes[i], op_axes[j])
                        factor = 1.0 if dot > 0.0 else -1.0
                        angle_radian = safe_arccos(abs(dot))

                        if abs(angle_radian) < self._tolerance_ang:
                            sum_axis += factor * op_axes[j]
                            count_axis += 1

            sym_idx = symmetry_indices[i]
            symmetry_operations[sym_idx]._axis = sum_axis / count_axis

        # update parallel and perpendicular symmetry elements

        # TODO: idealize symmetry elements for spherical case

        for i in range(len(symmetry_operations)):
            if isinstance(symmetry_operations[i], Inversion):
                continue

            for j in range(len(symmetry_operations)):
                if isinstance(symmetry_operations[j], Inversion):
                    continue

                ax_i = np.array(symmetry_operations[i]._axis)
                ax_j = np.array(symmetry_operations[j]._axis)

                dot = np.dot(ax_i, ax_j)
                factor = 1.0 if dot > 0.0 else -1.0
                angle_radian = safe_arccos(abs(dot))

                # perpendicular
                if abs(angle_radian - np.pi / 2.0) < self._tolerance_ang:
                    y_axis = np.cross(ax_j, ax_i)
                    y_axis /= np.linalg.norm(y_axis)
                    symmetry_operations[i]._axis = np.cross(y_axis, ax_j)
                    break

                # parallel
                elif abs(angle_radian) < self._tolerance_ang:
                    symmetry_operations[i]._axis = np.array(factor * ax_j)
                    break

        # find unique atoms

        redundant_atoms = set()
        for i, j in symmetry_mapping:
            if i < j:
                redundant_atoms.add(j)
            else:
                redundant_atoms.add(i)

        unique_atoms = set()
        for i in range(self._natoms):
            if i not in redundant_atoms:
                unique_atoms.add(i)

        # align unique atoms to symmetry elements

        tol_sq = self._tolerance_eig**2

        atoms = []
        atom_orig_indices = []
        for i in unique_atoms:
            atoms.append(np.array(self._centered_coords[i]))
            atom_orig_indices.append(i)

        # in case central atom is very close to origin
        for a in range(len(atoms)):
            vec_a = np.array(atoms[a])
            if np.sum(vec_a**2) < tol_sq:
                for b in range(len(atoms)):
                    atoms[b] -= vec_a

        for a in range(len(atoms)):
            vec_a = np.array(atoms[a])
            if np.sum(vec_a**2) < tol_sq:
                continue

            vec_a_norm = np.linalg.norm(vec_a)
            u_vec_a = vec_a / vec_a_norm

            for sym_op in symmetry_operations:
                if isinstance(sym_op, (Rotation, ImproperRotation)):
                    axis = np.array(sym_op._axis)

                    dot = np.dot(u_vec_a, axis)
                    factor = 1.0 if dot > 0.0 else -1.0
                    angle_radian = safe_arccos(abs(dot))

                    if abs(angle_radian) < self._tolerance_ang:
                        atoms[a] = np.array(factor * axis * vec_a_norm)
                        break

                elif isinstance(sym_op, (Reflection, ImproperRotation)):
                    axis = np.array(sym_op._axis)

                    dot = np.dot(u_vec_a, axis)
                    angle_radian = safe_arccos(abs(dot))

                    if abs(angle_radian - np.pi / 2.0) < self._tolerance_ang:
                        y_axis = np.cross(axis, u_vec_a)
                        y_axis /= np.linalg.norm(y_axis)
                        atoms[a] = np.array(np.cross(y_axis, axis) * vec_a_norm)
                        break

        # generate molecule from unique atoms

        for sym_op in symmetry_operations:
            sym_mat = sym_op.get_matrix()
            op_coords = np.matmul(np.array(atoms), sym_mat)

            for idx_i, vec_i in enumerate(op_coords):
                idx_a = None
                for a in range(self._natoms):
                    r2_ia = np.sum((vec_i - self._centered_coords[a])**2)
                    if r2_ia < tol_sq:
                        idx_a = a
                        break

                if idx_a is not None and idx_a not in atom_orig_indices:
                    atoms.append(np.array(vec_i))
                    atom_orig_indices.append(idx_a)

        assert_msg_critical(
            len(atoms) == self._natoms and
            len(atom_orig_indices) == self._natoms,
            'SymmetryAnalyzer.symmetrize_pointgroup: ' +
            'Inconsistent number of atoms')

        for idx_i, coord_i in zip(atom_orig_indices, atoms):
            centered_mol.set_atom_coordinates(idx_i, coord_i)

        # reorient molecule

        z_axis = None
        for sym_op in symmetry_operations:
            if isinstance(sym_op, (Rotation, ImproperRotation)):
                z_axis = np.array(sym_op._axis)
                break
        if z_axis is not None:
            min_dot, min_idx = None, None
            for i in range(self._natoms):
                if np.linalg.norm(atoms[i]) < 1e-8:
                    continue
                u_vec_i = atoms[i] / np.linalg.norm(atoms[i])
                dot = abs(np.dot(z_axis, u_vec_i))
                if min_dot is None or min_dot > dot:
                    min_dot = dot
                    min_idx = i
            y_axis = np.cross(z_axis, atoms[min_idx])
            y_axis /= np.linalg.norm(y_axis)
            x_axis = np.cross(y_axis, z_axis)
            rot_mat = np.array([x_axis, y_axis, z_axis])
            rot_coords = np.matmul(centered_mol.get_coordinates_in_bohr(),
                                   rot_mat.T)
            for i in range(self._natoms):
                centered_mol.set_atom_coordinates(i, rot_coords[i])

        return centered_mol

    def _handle_linear(self):
        """
        Handle linear molecules.

        :return:
            The Schoenflies symbol, list of detected elements, and reoriented
            coordinates.
        """

        # Set type of molecule for reorientation
        self._molecule_type = "linear"

        idx = np.argmin(self._Ivals)
        principal_axis = self._Ivecs[idx]
        p_axis = self._get_perpendicular(principal_axis)

        # Detect C2 axis along principal axis (for symmetrization)
        self._check_symmetry_operation(Rotation(principal_axis, order=2), "C2")

        # Detect sigma plane along principal axis (for symmetrization)
        self._check_symmetry_operation(Reflection(p_axis), "sigma_v")

        # Check for inversion center at center of mass
        if self._check_symmetry_operation(Inversion(), "i"):
            self._schoenflies_symbol = 'Dinfh'

            # Detect sigma plane and C2 axis perpendicular to principal axis
            # (for symmetrization)
            self._check_symmetry_operation(Rotation(p_axis, order=2), "C2")
            self._check_symmetry_operation(Reflection(principal_axis),
                                           "sigma_h")

        else:
            self._schoenflies_symbol = 'Cinfv'

        # Set orientation
        self._set_orientation(principal_axis, p_axis)

        # Set primary axis as z cartesian axis for conventional reorientation
        self._primary_axis = [0., 0., 1.]

    def _handle_asymmetric(self):
        """
        Handle asymmetric top molecules.
        """

        # Set orientation
        self._set_orientation(self._Ivecs[0], self._Ivecs[1])

        # Check for any C2 axis
        n_axis_c2 = 0
        principal_axis = [1, 0, 0]
        for axis in np.identity(3):
            c2 = Rotation(axis, order=2)
            if self._check_symmetry_operation(c2, "C2"):
                n_axis_c2 += 1
                principal_axis = axis
                break

        p_axis = self._get_perpendicular(principal_axis)
        for angle in np.arange(0, np.pi + self._tolerance_ang,
                               0.01 * np.pi / 2):
            axis = np.dot(p_axis, rotation_matrix(principal_axis, angle))
            c2 = Rotation(axis, order=2)
            if self._check_symmetry_operation(c2, "C2_p"):
                n_axis_c2 += 1
                self._secondary_axis = p_axis
                break

        self._max_order = 2

        if n_axis_c2 == 0:
            self._max_order = 0
            self._handle_no_rotation_axis()
        elif n_axis_c2 == 1:
            self._handle_cyclic(principal_axis)
        else:
            self._handle_dihedral(principal_axis)

    def _handle_symmetric(self):
        """
        Handle symmetric molecules.
        """

        # Get the only non-degenareted principal moment fo inertia and set the
        # principal axis along the associated eigenvector
        idx = self._get_nondegenerate(self._Ivals, self._tolerance_eig)
        principal_axis = self._Ivecs[idx]

        # Determine the highest possible rotation axis order along the principal axis
        self._max_order = self._get_axis_rot_order(principal_axis, n_max=9)

        # Check for C2 axis along principal axis
        # Imperative for symmetrization in abelian with real characters subgroups
        self._check_symmetry_operation(Rotation(principal_axis, order=2), "C2")

        # Get the perpendicualar axis to principal axis and check for C2 rotation axis
        # along p_axis by rotating p_axis along the principal axis
        p_axis = self._get_perpendicular(principal_axis)
        for angle in np.arange(0, np.pi + self._tolerance_ang,
                               0.01 * np.pi / self._max_order):
            axis = np.dot(p_axis, rotation_matrix(principal_axis, angle))
            c2 = Rotation(axis, order=2)
            if self._check_symmetry_operation(c2, "C2_p"):
                self._handle_dihedral(principal_axis)
                return

        self._handle_cyclic(principal_axis)

    def _handle_spherical(self):
        """
        Handle spherical groups (I, O, T) in iterative way by increasing
        tolerance if no axis is found.
        """

        # Set type of molecule for reorientation
        self._molecule_type = "spherical"
        # Set additional axis for reorientation in D2x group
        self._main_axis = [1., 0., 0.]
        self._p_axis = [0., 1., 0.]

        principal_axis = None

        mol_grid_points = self._get_mol_grid_points()

        # Check for C5 axis
        for axis in mol_grid_points:
            c5 = Rotation(axis, order=5)
            if self._check_symmetry_operation(c5, "C5"):
                self._schoenflies_symbol = "I"
                principal_axis = axis
                self._max_order = 5
                break

        # Check for C4 axis
        if principal_axis is None:
            for axis in mol_grid_points:
                c4 = Rotation(axis, order=4)
                if self._check_symmetry_operation(c4, "C4"):
                    self._schoenflies_symbol = "O"
                    principal_axis = axis
                    self._max_order = 4
                    break

        # Check for C3 axis
        if principal_axis is None:
            for axis in mol_grid_points:
                c3 = Rotation(axis, order=3)
                if self._check_symmetry_operation(c3, "C3"):
                    self._schoenflies_symbol = "T"
                    principal_axis = axis
                    self._max_order = 3
                    break

        assert_msg_critical(
            principal_axis is not None,
            'SymmetryAnalyzer: Could not find principal axis for spherical group'
        )

        p_axis_base = self._get_perpendicular(principal_axis)

        # I or Ih
        if self._schoenflies_symbol == 'I':

            def determine_orientation_I(principal_axis):
                """
                Determine the orientation of the molecule within the
                icosahedral symmetry group.  The angle is derived from
                geometric considerations to align the orientation axis with the
                symmetry axis.

                :param principal_axis:
                    The principal axis obtained with _set_orientation.

                :return:
                    The current orientation in the reference frame.
                """

                r_matrix = rotation_matrix(
                    p_axis_base, np.arcsin((np.sqrt(5) + 1) / (2 * np.sqrt(3))))
                axis = np.dot(principal_axis, r_matrix.T)

                # set molecule orientation in I
                for angle in np.arange(0, 2 * np.pi + self._tolerance_ang,
                                       self._tolerance_ang):
                    rot_matrix = rotation_matrix(principal_axis, angle)

                    c5_axis = np.dot(axis, rot_matrix.T)
                    c5 = Rotation(c5_axis, order=5)

                    if self._check_symmetry_operation(c5, "C5"):
                        t_axis = np.dot(
                            principal_axis,
                            rotation_matrix(p_axis_base, np.pi / 2).T)
                        return np.dot(t_axis, rot_matrix.T)

            # Set orientation
            p_axis = determine_orientation_I(principal_axis)
            self._main_axis = principal_axis
            self._p_axis = p_axis

            # Check for inverison center at center of mass
            if self._check_symmetry_operation(Inversion(), "i"):
                self._schoenflies_symbol += 'h'

            # Check for any C2 axis (for reorientation)
            for another_axis in self._get_mol_grid_points():
                c2 = Rotation(another_axis, order=2)
                if self._check_symmetry_operation(c2, "C2"):
                    self._primary_axis = another_axis

                    # Check for one reflexion plane (for symmetrization)
                    self._check_symmetry_operation(Reflection(another_axis),
                                                   "sigma_h")
                    break

            # Check for a C2 axis perpendicular to the first one (for reorientation)
            C2_p_axis = self._get_perpendicular(another_axis)
            for angle in np.arange(0, np.pi + self._tolerance_ang,
                                   0.01 * np.pi / 2):
                h_axis = np.dot(C2_p_axis, rotation_matrix(another_axis, angle))
                if self._check_symmetry_operation(Rotation(h_axis, order=2),
                                                  "C2_p"):
                    self._secondary_axis = h_axis
                    break

        # O or Oh
        if self._schoenflies_symbol == 'O':

            def determine_orientation_O(principal_axis):
                """
                Determine the orientation of the molecule within the octahedron
                symmetry group.

                :param principal_axis:
                    The principal axis obtained with _set_orientation.

                :return:
                    The current orientation in the reference frame.
                """

                r_matrix = rotation_matrix(p_axis_base, np.pi / 2)
                axis = np.dot(principal_axis, r_matrix.T)

                for angle in np.arange(0, 2 * np.pi + self._tolerance_ang,
                                       self._tolerance_ang):
                    rot_matrix = rotation_matrix(principal_axis, angle)

                    c4_axis = np.dot(axis, rot_matrix.T)
                    c4 = Rotation(c4_axis, order=4)

                    if self._check_symmetry_operation(c4, "C4"):
                        t_axis = np.dot(
                            principal_axis,
                            rotation_matrix(p_axis_base, np.pi / 2).T)
                        return np.dot(t_axis, rot_matrix.T)

            # Set orientation
            p_axis = determine_orientation_O(principal_axis)
            self._main_axis = principal_axis
            self._p_axis = p_axis

            # Check for inverison center at center of mass
            if self._check_symmetry_operation(Inversion(), "i"):
                self._schoenflies_symbol += 'h'

            # Check for any C2 axis (for reorientation)
            if self._check_symmetry_operation(
                    Rotation(self._main_axis, order=2), "C2"):

                # Check for one reflexion plane (for symmetrization)
                self._check_symmetry_operation(Reflection(self._main_axis),
                                               "sigma_h")

            # Check for a C2 axis perpendicular to the first one (for reorientation)
            C2_p_axis = self._get_perpendicular(principal_axis)
            for angle in np.arange(0, np.pi + self._tolerance_ang,
                                   0.01 * np.pi / 2):
                h_axis = np.dot(C2_p_axis,
                                rotation_matrix(principal_axis, angle))
                if self._check_symmetry_operation(Rotation(h_axis, order=2),
                                                  "C2_p"):
                    break

        # T or Td, Th
        if self._schoenflies_symbol == 'T':

            def determine_orientation_T(principal_axis):
                """
                Determine the orientation of the molecule within the
                tetrahedron symmetry group.

                :param principal_axis:
                    The principal axis obtained with _set_orientation.

                :return:
                    The current orientation in the reference frame.
                """

                r_matrix = rotation_matrix(p_axis_base, -np.arccos(-1 / 3))
                axis = np.dot(principal_axis, r_matrix.T)

                for angle in np.arange(0, 2 * np.pi + self._tolerance_ang,
                                       self._tolerance_ang):
                    rot_matrix = rotation_matrix(principal_axis, angle)

                    c3_axis = np.dot(axis, rot_matrix.T)
                    c3 = Rotation(c3_axis, order=3)

                    if self._check_symmetry_operation(c3, "C3"):
                        t_axis = np.dot(
                            principal_axis,
                            rotation_matrix(p_axis_base, np.pi / 2).T)
                        return np.dot(t_axis, rot_matrix.T)

            p_axis = determine_orientation_T(principal_axis)
            self._main_axis = principal_axis
            self._p_axis = p_axis

            # Check for any C2 axis (for reorientation)
            for another_axis in self._get_mol_grid_points():
                c2 = Rotation(another_axis, order=2)
                if self._check_symmetry_operation(c2, "C2"):
                    self._primary_axis = another_axis
                    break

            # Check for a C2 axis perpendicular to the first one (for reorientation)
            C2_p_axis = self._get_perpendicular(principal_axis)
            for angle in np.arange(0, np.pi + self._tolerance_ang,
                                   0.01 * np.pi / 2):
                h_axis = np.dot(C2_p_axis,
                                rotation_matrix(principal_axis, angle))
                if self._check_symmetry_operation(Rotation(h_axis, order=2),
                                                  "C2_p"):
                    self._p_axis = h_axis
                    break

            # Check for inverison center at center of mass
            if self._check_symmetry_operation(Inversion(), "i"):
                self._schoenflies_symbol += 'h'
                return

            # Check for any reflexion plane
            # C2_p_axis = self._get_perpendicular(principal_axis)
            for angle in np.arange(0, np.pi + self._tolerance_ang,
                                   0.01 * np.pi / 2):
                h_axis = np.dot(C2_p_axis,
                                rotation_matrix(principal_axis, angle))
                if self._check_symmetry_operation(Reflection(h_axis),
                                                  "sigma_v"):
                    self._schoenflies_symbol += 'd'
                    return

    def _handle_no_rotation_axis(self):
        """
        Detect point group of molecule with no rotation axis.

        :return:
            The Schoenflies symbol and list of detected elements.
        """

        # Set type of molecule for reorientation
        self._molecule_type = "asym_top"

        for i, vector in enumerate(np.identity(3)):
            if self._check_symmetry_operation(Reflection(vector), "sigma"):
                self._schoenflies_symbol = 'Cs'
                self._primary_axis = vector
                self._secondary_axis = self._get_perpendicular(vector)
                self._set_orientation(vector, self._secondary_axis)
                break
            else:
                if self._check_symmetry_operation(Inversion(), "i"):
                    self._schoenflies_symbol = 'Ci'
                    break
                else:
                    self._schoenflies_symbol = 'C1'

    def _handle_cyclic(self, principal_axis):
        """
        Detect point group of cylic group molecules.

        :param principal_axis:
            The principal axis obtained with _set_orientation.

        :return:
            The Schoenflies symbol and list of detected elements.
        """

        # Set type of molecule for reorientation
        self._molecule_type = "cyclic"
        self._primary_axis = principal_axis

        self._schoenflies_symbol = "C{}".format(self._max_order)

        # Check for reflexion planes perpenducular to the principal axis
        h_symbols = False
        if self._check_symmetry_operation(Reflection(principal_axis),
                                          "sigma_h"):
            self._schoenflies_symbol += 'h'
            h_symbols = True

        p_axis = self._get_perpendicular(principal_axis)
        self._secondary_axis = p_axis

        # Check for reflexion planes containing the principal axis
        v_symbol = set()
        for angle in np.arange(
                0, np.pi, 0.01 * np.pi / self._max_order + self._tolerance_ang):
            axis = np.dot(p_axis, rotation_matrix(principal_axis, angle))
            if self._check_symmetry_operation(Reflection(axis),
                                              "sigma_v") and not h_symbols:
                v_symbol.add('v')

        self._schoenflies_symbol += ''.join(v_symbol)

        # Check for inversion center at center of mass
        self._check_symmetry_operation(Inversion(), "i")

        # Check for improper rotational axis along the principal axis.
        if self._check_symmetry_operation(
                ImproperRotation(principal_axis, order=2 * self._max_order),
                "S{}".format(2 * self._max_order)):
            self._schoenflies_symbol = "S{}".format(2 * self._max_order)

    def _handle_dihedral(self, principal_axis):
        """
        Detect point group of dihedral group molecules.

        :param principal_axis:
            The principal axis obtained with _set_orientation.

        :return:
            The Schoenflies symbol and list of detected elements.
        """

        # Set type of molecule for reorientation
        self._molecule_type = "cyclic"
        self._primary_axis = principal_axis

        # Determine perpendicular axis to principal axis
        p_axis = self._get_perpendicular(principal_axis)

        if self._max_order == 1:
            # D1 is equivalent to C2
            self._schoenflies_symbol = "C2"
        else:
            self._schoenflies_symbol = "D{}".format(self._max_order)

        # Check for inversion center at center of mass
        self._check_symmetry_operation(Inversion(), "i")

        # Check for reflexion planes perpenducular to the principal axis
        h_symbols = False
        if self._check_symmetry_operation(Reflection(principal_axis),
                                          "sigma_h") and not h_symbols:
            self._schoenflies_symbol += 'h'
            h_symbols = True

        # Check for reflexion planes containing the principal axis
        d_symbol = False
        for angle in np.arange(0, np.pi + self._tolerance_ang,
                               0.01 * np.pi / self._max_order):
            axis = np.dot(p_axis, rotation_matrix(principal_axis, angle))
            if self._check_symmetry_operation(
                    Reflection(axis),
                    "sigma_v") and not h_symbols and not d_symbol:
                self._schoenflies_symbol += 'd'
                d_symbol = True

    def _get_axis_rot_order(self, axis, n_max):
        """
            Get rotation order for a given axis.

            :param axis:
                The axis.

            :param n_max:
                Maximum order to scan.

            :return:
                The order.
            """

        def max_rotation_order(tolerance):
            """
                Set the range of maximum order possible.

                :param tolerance:
                    The tolerance parameter.
                """

            for i in range(2, 15):
                if 2 * np.pi / (i * (i - 1)) <= tolerance:
                    return i - 1

        n_max = np.min([max_rotation_order(self._tolerance_ang), n_max])

        for i in range(n_max, 1, -1):
            Cn = Rotation(axis, order=i)
            if self._check_symmetry_operation(Cn, "C{}".format(i)):
                return i
        return 1

    def _check_symmetry_operation(self,
                                  operation,
                                  element_string,
                                  mapping=False):
        """
        Check if the given symmetry operation exists in the point group of the molecule.

        :param operation:
            The matrix representative of the symmetry operation.

        :return:
            True if the symmetry operation exists in the point group, False otherwise.
        """

        # Get representative of the operation
        sym_matrix = operation.get_matrix()

        assert_msg_critical(
            np.max(np.abs(sym_matrix - np.identity(3))) >= 1e-8,
            'SymmetryAnalyzer: Identity matrix should not be checked')

        # Get COM frame coordinates after the operation
        op_coordinates = np.matmul(self._centered_coords, sym_matrix)

        self._mapping = set()

        sym_op_exists = True
        tol_sq = self._tolerance_eig**2

        for i in range(self._natoms):
            min_r2, min_idx = None, None
            for j in range(self._natoms):
                if self._symbols[i] == self._symbols[j]:
                    rvec = op_coordinates[i] - self._centered_coords[j]
                    r2 = np.sum(rvec**2)
                    if min_r2 is None or min_r2 > r2:
                        min_r2 = r2
                        min_idx = j

            if min_r2 < tol_sq:
                if self._symbols[i] == self._symbols[min_idx] and i != min_idx:
                    self._mapping.add(tuple(sorted([i, min_idx])))
            else:
                sym_op_exists = False
                if not mapping:
                    self._mapping = None
                    break

        return sym_op_exists

    def _set_orientation(self, principal_axis, p_axis):
        """
        Set molecular orientation along principal_axis (x) and p_axis (y).

        :param principal_axis:
            Principal orientation axis (must be unitary).
        :param p_axis:
            Secondary axis perpendicular to principal (must be unitary).
        """

        assert_msg_critical(
            abs(np.linalg.norm(principal_axis) - 1.0) < 1e-8,
            "SymmetryAnalyzer: Principal axis is not unitary.")

        assert_msg_critical(
            abs(np.linalg.norm(p_axis) - 1.0) < 1e-8,
            "SymmetryAnalyzer: p_axis is not unitary.")

        orientation = np.array([
            principal_axis,
            p_axis,
            np.cross(principal_axis, p_axis),
        ])

        self._centered_coords = np.matmul(self._centered_coords, orientation.T)

    @staticmethod
    def _get_degeneracy(Ivals, tolerance):
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
    def _get_nondegenerate(Ivals, tolerance):
        """
        Get the index of the nondegenerate eigenvalue from the array of
        eigenvalues.

        :param Ivals:
            The array of eigenvalues.

        :param tolerance:
            The tolerance parameter on the eigenvalues.

        :return:
            The index of the nondegenerate eigenvalue.
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

        assert_msg_critical(False, "SymmetryAnalyzer: Nondegenerate not found.")

    def _get_mol_grid_points(self):
        """
        Generate points surrounding the molecule.

        :return:
            List of points.
        """

        centered_mol = Molecule(self._symbols, self._centered_coords, 'au')
        connectivity_matrix = centered_mol.get_connectivity_matrix()
        Ivals, Ivecs = centered_mol.moments_of_inertia(principal_axes=True)

        points = []

        # principal axes
        for i in range(3):
            vec = np.array(Ivecs[i])
            norm = np.linalg.norm(vec)
            if norm > 1e-8:
                points.append(vec / norm)

        # atoms
        for i in range(self._natoms):
            vec = np.array(self._centered_coords[i])
            norm = np.linalg.norm(vec)
            if norm > 1e-8:
                points.append(vec / norm)

        # centers of rings
        graph = nx.Graph()
        for i in range(self._natoms):
            graph.add_node(i)
            for j in range(i + 1, self._natoms):
                if connectivity_matrix[i][j] == 1:
                    graph.add_edge(i, j)
        cycles = list(nx.simple_cycles(graph, length_bound=8))
        for cycle in cycles:
            vec = np.zeros(3)
            for i in cycle:
                vec += self._centered_coords[i]
            vec /= len(cycle)
            norm = np.linalg.norm(vec)
            if norm > 1e-8:
                points.append(vec / norm)

        # centers of bonds
        for i in range(self._natoms):
            for j in range(i + 1, self._natoms):
                if (connectivity_matrix[i][j] == 1 and
                        self._symbols[i] == self._symbols[j]):
                    vec = 0.5 * (self._centered_coords[i] +
                                 self._centered_coords[j])
                    norm = np.linalg.norm(vec)
                    if norm > 1e-8:
                        points.append(vec / norm)

        # find duplicate points
        duplicate_indices = []
        tol_sq = self._tolerance_eig**2
        for i in range(len(points)):
            p_i = points[i]

            for j in range(i + 1, len(points)):
                if j in duplicate_indices:
                    continue
                p_j = points[j]

                # check if p_i == p_j or p_i == -p_j
                # note: scale axis to approx. C-H bond length
                if np.dot(p_i, p_j) > 0.0:
                    rvec = (p_i - p_j) * 2.06
                else:
                    rvec = (p_i + p_j) * 2.06
                r2 = np.sum(rvec**2)
                if r2 < tol_sq:
                    duplicate_indices.append(j)

        # remove duplicate points
        unique_points = []
        for i, p in enumerate(points):
            if i not in duplicate_indices:
                unique_points.append(p)

        return unique_points

    @staticmethod
    def _get_perpendicular(vector, tol=1e-8):
        """
        Generate a vector perpendicular to another vector or axis.

        :param vector:
            The vector or axis with respect to which the perpendicular axis is
            determined.
        :param tol:
            An additional tolerance parameter to condisder the axis as perpendicular.

        :return:
            An array of coordinates of the perpendicular and normalized vector.
        """

        index = np.argmin(np.abs(vector))
        p_vector = np.identity(3)[index]
        pp_vector = np.cross(vector, p_vector)
        pp_vector = pp_vector / np.linalg.norm(pp_vector)

        # check perpendicular
        assert_msg_critical(
            abs(np.dot(pp_vector, vector)) < tol,
            "SymmetryAnalyzer: perpendicular check failed")
        # check normalization
        assert_msg_critical(
            abs(np.linalg.norm(pp_vector) - 1.0) < tol,
            "SymmetryAnalyzer: normalization check failed")

        return pp_vector
