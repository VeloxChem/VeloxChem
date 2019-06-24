import numpy as np
import os

from .veloxchemlib import Molecule


@staticmethod
def _Molecule_read_str(xyzstr, units='angs'):
    """Reads molecule from xyz string.

    Reads molecule from xyz string.

    Parameters
    ----------
    xyzstr
        The xyz string.
    units
        The unit of coordinates.

    Returns
    -------
    Molecule
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
    """Reads molecule from xyz file.

    Reads molecule from xyz file.

    Parameters
    ----------
    xyzfile
        The name of the xyz file.

    Returns
    -------
    Molecule
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
    """Reads molecule from a dictionary.

    Reads molecule from a dictionary.

    Parameters
    ----------
    mol_dict
        The molecule dictionary.

    Returns
    -------
    Molecule
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
    """Computes center of mass of a molecule.

    Computes center of mass of a molecule.

    Returns
    -------
    tuple
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


Molecule.read_str = _Molecule_read_str
Molecule.read_xyz = _Molecule_read_xyz
Molecule.from_dict = _Molecule_from_dict
Molecule.center_of_mass = _Molecule_center_of_mass
