from pathlib import Path
import numpy as np

from .veloxchemlib import Point
from .veloxchemlib import Molecule
from .veloxchemlib import bohr_in_angstrom
from .veloxchemlib import chemical_element_identifier

from .outputstream import OutputStream
from .inputparser import print_keywords
from .errorhandler import assert_msg_critical


@staticmethod
def _Molecule_smiles_to_xyz(smiles_str, optimize=True, hydrogen=True):
    """
    Converts SMILES string to xyz string.

    :param smiles_str:
        The SMILES string.
    :param optimize:
        Boolean indicating whether to perform geometry optimization.
    :param hydrogen:
        Boolean indicating whether to remove hydrogens.

    :return:
        An xyz string (including number of atoms).
    """

    try:
        from openbabel import pybel as pb

        mol = pb.readstring('smiles', smiles_str)
        mol.make3D()

        if optimize:
            # TODO: Double check if UFF is needed
            mol.localopt(forcefield="mmff94", steps=300)

        if not hydrogen:
            # remove hydrogens
            mol.removeh()
            return mol.write(format="xyz")

        else:
            return mol.write(format="xyz")

    except ImportError:
        raise ImportError('Unable to import openbabel')


@staticmethod
def _Molecule_read_smiles(smiles_str):
    """
    Reads molecule from SMILES string.

    :param smiles_str:
        The SMILES string.

    :return:
        The molecule.
    """

    xyz = Molecule.smiles_to_xyz(smiles_str, optimize=True)

    return Molecule.read_xyz_string(xyz)


@staticmethod
def _Molecule_read_molecule_string(mol_str, units='angstrom'):
    """
    Reads molecule from a string containing Cartesian coordinates.

    :param mol_str:
        The string containing Cartesian coordinates.
    :param units:
        The unit of coordinates.

    :return:
        The molecule.
    """

    labels = []
    coords = []

    for line in mol_str.strip().splitlines():
        if line:
            content = line.split()
            labels.append(content[0].upper())
            coords.append(Point([float(x) for x in content[1:4]]))

    elements = []

    for label in labels:
        elem_id = chemical_element_identifier(label)
        assert_msg_critical(
            elem_id != -1,
            f'Molecule: Unsupported chemical element {label} in XYZ string')
        elements.append(elem_id)

    return Molecule(elements, coords, units)


@staticmethod
def _Molecule_read_xyz_file(xyzfile):
    """
    Reads molecule from file in XYZ format.

    :param xyzfile:
        File with molecular structure in XYZ format.

    :return:
        The molecule.
    """

    with Path(xyzfile).open('r') as fh:
        xyzstr = fh.read()

    return Molecule.read_xyz_string(xyzstr)


@staticmethod
def _Molecule_read_xyz_string(xyz):
    """
    Generate molecule from string in XYZ format.

    :param xyz:
        String with XYZ structure.

    :return:
        The molecule.
    """

    lines = xyz.strip().splitlines()

    try:
        natoms = int(lines[0].strip())
    except (ValueError, TypeError):
        assert_msg_critical(False,
                            'Molecule: Invalid number of atoms in XYZ input')

    assert_msg_critical(natoms == len(lines[2:]),
                        'Molecule: Inconsistent number of atoms in XYZ input')

    mol_str = '\n'.join(lines[2:])
    return Molecule.read_molecule_string(mol_str, 'angstrom')


@staticmethod
def _Molecule_from_dict(mol_dict):
    """
    Reads molecule from a dictionary.

    :param mol_dict:
        The molecule dictionary.

    :return:
        The molecule.
    """

    assert_msg_critical('xyz' in mol_dict or 'xyzfile' in mol_dict,
                        'Molecule: Expecting either "xyz" or "xyzfile" input')

    assert_msg_critical(not ('xyz' in mol_dict and 'xyzfile' in mol_dict),
                        'Molecule: Cannot have both "xyz" and "xyzfile" input')

    if 'xyz' in mol_dict:
        mol_str = '\n'.join(mol_dict['xyz'])
        units = 'angstrom'
        if 'units' in mol_dict:
            units = mol_dict['units'].lower()
        mol = Molecule.read_molecule_string(mol_str, units)

    elif 'xyzfile' in mol_dict:
        assert_msg_critical(
            'units' not in mol_dict,
            'Molecule: Cannot have both "units" and "xyzfile" input')
        mol = Molecule.read_xyz_file(mol_dict['xyzfile'])

    charge = 0.0
    if 'charge' in mol_dict:
        charge = float(mol_dict['charge'])

    multiplicity = 1
    if 'multiplicity' in mol_dict:
        multiplicity = int(mol_dict['multiplicity'])

    mol.set_charge(charge)
    mol.set_multiplicity(multiplicity)

    assert_msg_critical(
        mol.check_multiplicity(),
        'Molecule: Incompatible multiplicity and number of electrons')
    assert_msg_critical(
        mol.check_proximity(0.1),
        'Molecule: Corrupted geometry with closely located atoms')

    return mol


def _Molecule_center_of_mass_in_bohr(self):
    """
    Computes center of mass of a molecule in Bohr.

    :return:
        The center of mass in Bohr.
    """

    masses = np.array(self.get_masses())
    coords = self.get_coordinates_in_bohr()

    x_center = np.sum(coords[:, 0] * masses) / np.sum(masses)
    y_center = np.sum(coords[:, 1] * masses) / np.sum(masses)
    z_center = np.sum(coords[:, 2] * masses) / np.sum(masses)

    return np.array([x_center, y_center, z_center])


def _Molecule_center_of_mass_in_angstrom(self):
    """
    Computes center of mass of a molecule in Angstrom.

    :return:
        The center of mass in Angstrom.
    """

    return self.center_of_mass_in_bohr() * bohr_in_angstrom()


def _Molecule_get_string(self):
    """
    Returns string representation of molecule.

    :return:
        A string with representation of molecule.
    """

    labels = self.get_labels()
    coords_in_angstrom = self.get_coordinates_in_angstrom()

    mol_str = 'Molecular Geometry (Angstroms)\n'
    mol_str += '================================\n\n'
    mol_str += '  Atom         Coordinate X          Coordinate Y          Coordinate Z  \n\n'
    for label, coords in zip(labels, coords_in_angstrom):
        mol_str += f'  {label:<4s}{coords[0]:>22.12f}{coords[1]:>22.12f}{coords[2]:>22.12f}\n'
    mol_str += '\n'

    return mol_str


def _Molecule_more_info(self):
    """
    Returns more information about the molecule.

    :return:
        Molecular information in plain text.
    """

    width = 70
    mol_info = []

    mol_info.append(
        f'Molecular charge            : {self.get_charge():.0f}'.ljust(width))
    mol_info.append(
        f'Spin multiplicity           : {self.get_multiplicity():d}'.ljust(
            width))
    mol_info.append(
        f'Number of atoms             : {self.number_of_atoms():d}'.ljust(
            width))
    mol_info.append(
        f'Number of alpha electrons   : {self.number_of_alpha_electrons():d}'.
        ljust(width))
    mol_info.append(
        f'Number of beta  electrons   : {self.number_of_beta_electrons():d}'.
        ljust(width))

    return '\n'.join(mol_info)


def _Molecule_get_coordinates_in_bohr(self):
    """
    Returns atom coordinates in Bohr.

    :return:
        A numpy array of atom coordinates (nx3) in Bohr.
    """

    coords = []
    for r in self.get_coordinates():
        coords.append(r.coordinates())

    return np.array(coords)


def _Molecule_get_coordinates_in_angstrom(self):
    """
    Returns atom coordinates in Angstrom.

    :return:
        A numpy array of atom coordinates (nx3) in Angstrom.
    """

    coords = []
    for r in self.get_coordinates('angstrom'):
        coords.append(r.coordinates())

    return np.array(coords)


def _Molecule_get_xyz_string(self):
    """
    Returns xyz string of molecule.

    :return:
        An xyz string (including number of atoms).
    """

    labels = self.get_labels()
    coords_in_angstrom = self.get_coordinates_in_angstrom()

    natoms = len(labels)
    xyz = f'{natoms}\n\n'

    for a in range(natoms):
        xa, ya, za = coords_in_angstrom[a]
        xyz += f'{labels[a]:<6s} {xa:22.12f} {ya:22.12f} {za:22.12f}\n'
    return xyz


def _Molecule_write_xyz_file(self, xyz_filename):
    """
    Writes molecular geometry to xyz file.

    :param xyz_filename:
        The name of the xyz file.
    """

    with open(str(xyz_filename), 'w') as fh:
        fh.write(self.get_xyz_string())


def _Molecule_show(self,
                   width=400,
                   height=300,
                   atom_indices=False,
                   atom_labels=False):
    """
    Creates a 3D view with py3dmol.

    :param width:
        The width.
    :param height:
        The height.
    :param atom_indices:
        The flag for showing atom indices (1-based).
    :param atom_labels:
        The flag for showing atom labels.
    """

    try:
        import py3Dmol
        viewer = py3Dmol.view(width=width, height=height)
        viewer.addModel(self.get_xyz_string())
        viewer.setViewStyle({"style": "outline", "width": 0.05})
        viewer.setStyle({"stick": {}, "sphere": {"scale": 0.25}})
        if atom_indices or atom_labels:
            coords = self.get_coordinates_in_angstrom()
            labels = self.get_labels()
            for i in range(coords.shape[0]):
                text = ''
                if atom_labels:
                    text += f'{labels[i]}'
                if atom_indices:
                    text += f'{i + 1}'
                viewer.addLabel(
                    text, {
                        'position': {
                            'x': coords[i, 0],
                            'y': coords[i, 1],
                            'z': coords[i, 2],
                        },
                        'alignment': 'center',
                        'fontColor': 0x000000,
                        'backgroundColor': 0xffffff,
                        'backgroundOpacity': 0.0,
                    })
        viewer.zoomTo()
        viewer.show()

    except ImportError:
        raise ImportError('Unable to import py3Dmol')


def _Molecule_draw_2d(self, width=400, height=300):
    """
    Generates 2D representation of the molecule.

    :param width:
        The width.
    :param height:
        The height.
    """

    try:
        from openbabel import pybel as pb
        from IPython.display import SVG, display

        molecule = self.get_xyz_string()

        mol = pb.readstring('xyz', molecule)

        mol.make2D()
        mol.removeh()

        # Convert to SVG using pybel's drawing method
        svg_string = mol.write(format='svg', opt={'w': width, 'h': height})

        # Display SVG
        display(SVG(svg_string))

    except ImportError:
        raise ImportError('Unable to import openbabel and/or IPython.display.')


def _Molecule_moments_of_inertia(self):
    """
    Calculates the moment of inertia tensor and principle axes.

    :return:
        The principle moments of inertia.
    """

    masses = np.array(self.get_masses())
    coordinates = self.get_coordinates_in_bohr()
    center_of_mass = np.array(self.center_of_mass_in_bohr())
    natm = self.number_of_atoms()

    # Coordinates in the center-of-mass frame
    coords_com = coordinates - center_of_mass[np.newaxis, :]

    # Moment of inertia tensor
    Imat_atoms = [
        masses[i] * (np.eye(3) * (np.dot(coords_com[i], coords_com[i])) -
                     np.outer(coords_com[i], coords_com[i]))
        for i in range(natm)
    ]
    Imom = np.sum(Imat_atoms, axis=0)

    # Principal moments
    Ivals, Ivecs = np.linalg.eigh(Imom)
    # Eigenvectors are in the rows after transpose
    # Ivecs = Ivecs.T

    return Ivals


def _Molecule_is_linear(self):
    """
    Checks if a molecule is linear or not.

    :return:
        True if linear, False otherwise.
    """

    assert_msg_critical(self.number_of_atoms() >= 2,
                        'Molecule.is_linear: Need at least two atoms')

    # Get principle moments of inertia
    Ivals = self.moments_of_inertia()

    # Obtain the number of rotational degrees of freedom (DoF)
    Rotational_DoF = 0
    for i in range(3):
        if abs(Ivals[i]) > 1.0e-10:
            Rotational_DoF += 1

    assert_msg_critical(
        Rotational_DoF in [2, 3],
        'Molecule.is_linear: Unexpected rotational degrees of freedom')

    if Rotational_DoF == 2:
        return True
    elif Rotational_DoF == 3:
        return False


def _Molecule_get_aufbau_alpha_occupation(self, n_mo):
    """
    Gets occupation numbers for alpha spin based on the aufbau principle.

    :param n_mo:
        The number of molecular orbitals.

    :return:
        The occupation numbers for alpha spin.
    """

    nalpha = self.number_of_alpha_electrons()

    return np.hstack((np.ones(nalpha), np.zeros(n_mo - nalpha)))


def _Molecule_get_aufbau_beta_occupation(self, n_mo):
    """
    Gets occupation numbers for beta spin based on the aufbau principle.

    :param n_mo:
        The number of molecular orbitals.

    :return:
        The occupation numbers for beta spin.
    """

    nbeta = self.number_of_beta_electrons()

    return np.hstack((np.ones(nbeta), np.zeros(n_mo - nbeta)))


def _Molecule_get_aufbau_occupation(self, n_mo, flag='restricted'):
    """
    Gets occupation vector(s) based on the aufbau principle.

    :param n_mo:
        The number of molecular orbitals.
    :param flag:
        The flag (restricted or unrestricted).

    :return:
        The occupation vector(s).
    """

    occ_a = self.get_aufbau_alpha_occupation(n_mo)
    occ_b = self.get_aufbau_beta_occupation(n_mo)

    if flag == 'restricted':
        return 0.5 * (occ_a + occ_b)

    elif flag == 'unrestricted':
        return occ_a, occ_b

    return None


@staticmethod
def _Molecule_get_input_keywords():
    """
    Returns input keywords for Molecule.
    """

    return {
        'molecule': {
            'charge': ('int', 'net charge'),
            'multiplicity': ('int', 'spin multiplicity'),
            'units': ('str_lower', 'unit of coordinates, default is Angstrom'),
            'xyz': ('list', 'atom and Cartesian coordinates'),
            'xyzfile': ('str', 'XYZ file name (conflicts with units/xyz)'),
        },
    }


@staticmethod
def _Molecule_print_keywords():
    """
    Prints keywords for Molecule.
    """

    input_keywords = Molecule._get_input_keywords()
    ostream = OutputStream()

    print_keywords(input_keywords, ostream)


def _Molecule_check_multiplicity(self):
    """
    Returns True if multiplicity and charge of molecule consistent, False otherwise.
    """

    multip = self.get_multiplicity() % 2
    nelec = self.number_of_electrons() % 2

    if (multip == 0) and (nelec != 1):
        return False
    if (multip == 1) and (nelec != 0):
        return False
    return True


def _Molecule_number_of_alpha_electrons(self):
    """
    Returns number of alpha electrons in molecule.
    """

    return (self.number_of_electrons() + self.get_multiplicity() - 1) // 2


def _Molecule_number_of_beta_electrons(self):
    """
    Returns number of beta electrons in molecule.
    """

    return (self.number_of_electrons() - self.get_multiplicity() + 1) // 2


Molecule._get_input_keywords = _Molecule_get_input_keywords

Molecule.smiles_to_xyz = _Molecule_smiles_to_xyz
Molecule.show = _Molecule_show
Molecule.draw_2d = _Molecule_draw_2d
Molecule.read_smiles = _Molecule_read_smiles
Molecule.read_molecule_string = _Molecule_read_molecule_string
Molecule.read_xyz_file = _Molecule_read_xyz_file
Molecule.read_xyz_string = _Molecule_read_xyz_string
Molecule.from_dict = _Molecule_from_dict
Molecule.center_of_mass_in_bohr = _Molecule_center_of_mass_in_bohr
Molecule.center_of_mass_in_angstrom = _Molecule_center_of_mass_in_angstrom
Molecule.more_info = _Molecule_more_info
Molecule.get_coordinates_in_bohr = _Molecule_get_coordinates_in_bohr
Molecule.get_coordinates_in_angstrom = _Molecule_get_coordinates_in_angstrom
Molecule.get_xyz_string = _Molecule_get_xyz_string
Molecule.write_xyz_file = _Molecule_write_xyz_file
Molecule.moments_of_inertia = _Molecule_moments_of_inertia
Molecule.is_linear = _Molecule_is_linear
Molecule.get_aufbau_alpha_occupation = _Molecule_get_aufbau_alpha_occupation
Molecule.get_aufbau_beta_occupation = _Molecule_get_aufbau_beta_occupation
Molecule.get_aufbau_occupation = _Molecule_get_aufbau_occupation
Molecule.print_keywords = _Molecule_print_keywords
Molecule.check_multiplicity = _Molecule_check_multiplicity
Molecule.number_of_alpha_electrons = _Molecule_number_of_alpha_electrons
Molecule.number_of_beta_electrons = _Molecule_number_of_beta_electrons
Molecule.get_string = _Molecule_get_string

# aliases for backward compatibility
Molecule.read_xyz = _Molecule_read_xyz_file
Molecule.from_xyz_string = _Molecule_read_xyz_string
Molecule.write_xyz = _Molecule_write_xyz_file
Molecule.read_str = _Molecule_read_molecule_string
