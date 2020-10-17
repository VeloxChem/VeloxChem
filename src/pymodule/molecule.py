import numpy as np
import os

from .veloxchemlib import Molecule
from .veloxchemlib import ChemicalElement


@staticmethod
def _Molecule_read_str(xyzstr, units='angs'):
    """
    Reads molecule from xyz string.

    :param xyzstr:
        The xyz string.
    :param units:
        The unit of coordinates.

    :return:
        The molecule.
    """

    labels = []
    coords = []

    for line in xyzstr.strip().split(os.linesep):
        if line:
            content = line.split()
            labels.append(content[0])
            coords.append([float(x) for x in content[1:4]])

    return Molecule(labels, coords, units)


@staticmethod
def _Molecule_read_xyz(xyzfile):
    """
    Reads molecule from xyz file.

    :param xyzfile:
        The name of the xyz file.

    :return:
        The molecule.
    """

    xyzstr = ''

    with open(xyzfile, 'r') as f_xyz:
        natoms = int(f_xyz.readline().split()[0])
        f_xyz.readline()
        for a in range(natoms):
            xyzstr += f_xyz.readline().strip() + os.linesep

    return Molecule.read_str(xyzstr, 'angs')


@staticmethod
def _Molecule_from_dict(mol_dict):
    """
    Reads molecule from a dictionary.

    :param mol_dict:
        The molecule dictionary.

    :return:
        The molecule.
    """

    xyzstr = mol_dict['xyzstr']

    units = 'angs'
    if 'units' in mol_dict:
        units = mol_dict['units'].lower()

    charge = 0.0
    if 'charge' in mol_dict:
        charge = float(mol_dict['charge'])

    multiplicity = 1
    if 'multiplicity' in mol_dict:
        multiplicity = int(mol_dict['multiplicity'])

    mol = Molecule.read_str(xyzstr, units)
    mol.set_charge(charge)
    mol.set_multiplicity(multiplicity)
    mol.check_multiplicity()
    mol.check_proximity(0.1)

    return mol


def _Molecule_center_of_mass(self):
    """
    Computes center of mass of a molecule.

    :return:
        The center of mass.
    """

    masses = self.masses_to_numpy()
    x_coords = self.x_to_numpy()
    y_coords = self.y_to_numpy()
    z_coords = self.z_to_numpy()

    x_center = np.sum(x_coords * masses) / np.sum(masses)
    y_center = np.sum(y_coords * masses) / np.sum(masses)
    z_center = np.sum(z_coords * masses) / np.sum(masses)

    return x_center, y_center, z_center


def _Molecule_more_info(self):
    """
    Returns more information about the molecule.

    :return:
        Molecular information in plain text.
    """

    width = 70
    mol_info = []

    mol_info.append('Molecular charge            : {:.0f}'.format(
        self.get_charge()).ljust(width))
    mol_info.append('Spin multiplicity           : {:d}'.format(
        self.get_multiplicity()).ljust(width))
    mol_info.append('Number of atoms             : {:d}'.format(
        self.number_of_atoms()).ljust(width))
    mol_info.append('Number of alpha electrons   : {:d}'.format(
        self.number_of_alpha_electrons()).ljust(width))
    mol_info.append('Number of beta  electrons   : {:d}'.format(
        self.number_of_beta_electrons()).ljust(width))

    return os.linesep.join(mol_info)


def _Molecule_get_labels(self):
    """
    Returns atom labels.

    :return:
        A list of atom labels.
    """

    labels = []

    for elem_id in self.elem_ids_to_numpy():
        elem = ChemicalElement()
        elem.set_atom_type(elem_id)
        labels.append(elem.get_name())

    return labels


def _Molecule_get_coordinates(self):
    """
    Returns atom coordinates.

    :return:
        A numpy array of atom coordinates (nx3).
    """

    return np.array([
        self.x_to_numpy(),
        self.y_to_numpy(),
        self.z_to_numpy(),
    ]).T.copy()


Molecule.read_str = _Molecule_read_str
Molecule.read_xyz = _Molecule_read_xyz
Molecule.from_dict = _Molecule_from_dict
Molecule.center_of_mass = _Molecule_center_of_mass
Molecule.more_info = _Molecule_more_info
Molecule.get_labels = _Molecule_get_labels
Molecule.get_coordinates = _Molecule_get_coordinates
