from pathlib import Path
import numpy as np
import math as mt
import sys

from .veloxchemlib import Molecule
from .veloxchemlib import ChemicalElement
from .veloxchemlib import bohr_in_angstroms

@staticmethod
def _Molecule_read_str(xyzstr, units='angstrom'):
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
    
    for line in xyzstr.strip().splitlines():
        if line:
            content = line.split()
            labels.append(content[0])
            coords.append([float(content[1]),
                           float(content[2]),
                           float(content[3])])

    return Molecule(labels, coords, units)

@staticmethod
def _Molecule_read_xyz(xyzfile):
    """
    Reads molecule from file in XYZ format.

    :param xyzfile:
        File with molecular structure in XYZ format.

    :return:
        The molecule.
    """

    with Path(xyzfile).open('r') as fh:
        xyzstr = '\n'.join(fh.readlines()[2:])

    return Molecule.read_str(xyzstr)

@staticmethod
def _Molecule_from_xyz_string(xyz):
    """
    Generate molecule from string in XYZ format.

    :param xyz:
        String with XYZ structure.

    :return:
        The molecule.
    """

    xyzstr = '\n'.join(xyz.strip().splitlines()[2:])

    return Molecule.read_str(xyzstr)

@staticmethod
def _Molecule_from_dict(mol_dict):
    """
    Reads molecule from a dictionary.

    :param mol_dict:
        The molecule dictionary.

    :return:
        The molecule.
    """

    xyzstr = '\n'.join(mol_dict['xyz'])

    units = 'angstrom'
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
    if not mol.check_multiplicity():
        sys.exit('Molecule: Incompatible multiplicity and number of electrons.')
    if not mol.check_proximity(0.1):
        sys.exit('Molecule: Atoms are too close.')

    return mol

def _Molecule_center_of_mass(self):
    """
    Computes center of mass of a molecule.

    :return:
        The center of mass.
    """

    coords = self.get_coordinates()
    x_coords = np.array([r[0] for r in coords])
    y_coords = np.array([r[1] for r in coords])
    z_coords = np.array([r[2] for r in coords])

    masses = np.array(self.get_masses())
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

def _Molecule_write_xyz(self, xyz_filename):
    """
    Writes molecular geometry to xyz file.

    :param xyz_filename:
        The name of the xyz file.
    """

    with open(str(xyz_filename), 'w') as fh:

        fh.write(f"{self.number_of_atoms():d}\n\n")

        for name, r in zip(self.get_labels(), self.get_coordinates("angstrom")):
            fh.write(f'{name:<6s} {r[0]:22.12f} {r[1]:22.12f} {r[2]:22.12f}\n')
            
def _Molecule_check_multiplicity(self):
    """
    Checks if a molecule has consistent multiplicity.

    :return:
        True if multiplicity is consistent, False otherwise.
    """
    
    multip = self.get_multiplicity() % 2
    nelecs = self.number_of_electrons() % 2
    
    if (multip == 0) and (nelecs != 1):
        return False
    if (multip == 1) and (nelecs != 0):
        return False
    return True
    
def _Molecule_number_of_alpha_electrons(self):
    """
    Gets number of alpha electrons in molecule.

    :return:
        The number of alpha electrons.
    """
    
    multip = self.get_multiplicity() - 1
    nelecs = self.number_of_electrons()
    return (nelecs + multip) // 2;

def _Molecule_number_of_beta_electrons(self):
    """
    Gets number of beta electrons in molecule.

    :return:
        The number of beta electrons.
    """
    
    multip = self.get_multiplicity() - 1
    nelecs = self.number_of_electrons()
    return (nelecs - multip) // 2;
        
def _Molecule_moments_of_inertia(self):
    """
    Calculates the moment of inertia tensor and principle axes

    :return:
        The principle moments of inertia.
    """

    masses = np.array(self.get_masses())
    coordinates = np.array(self.get_coordinates())
    center_of_mass = np.array(self.center_of_mass())
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

    # Get principle moments of inertia
    Ivals = self.moments_of_inertia()

    # Obtain the number of rotational degrees of freedom (DoF)
    Rotational_DoF = 0
    for i in range(3):
        if abs(Ivals[i]) > 1.0e-10:
            Rotational_DoF += 1

    if Rotational_DoF == 2:
        return True
    elif Rotational_DoF == 3:
        return False
    # TODO: Raise an error if rotational DoFs are not 2 or 3
    else:
        pass

def _Molecule_get_aufbau_occupation(self, norb, flag='restricted'):
    """
    Creates an occupation vector based on the aufbau principle.

    :param norb:
        The number of molecular orbitals
    :param flag:
        The flag (restricted or unrestricted).

    :return:
        flag=='restricted': single vector assuming doubly occupied orbitals.
        flag=='unrestricted': two vectors assuming singly occupied spin-orbitals.
    """

    nalpha = self.number_of_alpha_electrons()
    nbeta = self.number_of_beta_electrons()

    if flag == 'restricted':
        occ = [
            2.0 if x < nbeta else 1.0 if x < nalpha else 0.0
            for x in range(norb)
        ]
        return occ
    else:
        occa = [1.0 if x < nalpha else 0.0 for x in range(norb)]
        occb = [1.0 if x < nbeta else 0.0 for x in range(norb)]
        return occa, occb

def _Molecule_deepcopy(self, memo):
    """
    Implements deepcopy.

    :param memo:
        The memo dictionary for deepcopy.

    :return:
        A deepcopy of self.
    """

    return Molecule(self)
    
Molecule.read_str = _Molecule_read_str
Molecule.read_xyz = _Molecule_read_xyz
Molecule.from_xyz_string = _Molecule_from_xyz_string
Molecule.from_dict = _Molecule_from_dict
Molecule.center_of_mass = _Molecule_center_of_mass
Molecule.number_of_alpha_electrons = _Molecule_number_of_alpha_electrons
Molecule.number_of_beta_electrons = _Molecule_number_of_beta_electrons
Molecule.more_info = _Molecule_more_info
Molecule.write_xyz = _Molecule_write_xyz
Molecule.check_multiplicity = _Molecule_check_multiplicity
Molecule.moments_of_inertia = _Molecule_moments_of_inertia
Molecule.is_linear = _Molecule_is_linear
Molecule.get_aufbau_occupation = _Molecule_get_aufbau_occupation
Molecule.__deepcopy__ = _Molecule_deepcopy
