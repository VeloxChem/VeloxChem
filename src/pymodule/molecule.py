#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

from pathlib import Path
import numpy as np

from .veloxchemlib import Molecule
from .veloxchemlib import ChemicalElement
from .veloxchemlib import bohr_in_angstrom

from .outputstream import OutputStream
from .inputparser import print_keywords
from .errorhandler import assert_msg_critical


@staticmethod
def _Molecule_smiles_to_xyz(smiles_str, optimize=True, no_hydrogen=False):
    """
    Converts SMILES string to xyz string.

    :param smiles_str:
        The SMILES string.

    :return:
        An xyz string (including number of atoms).
    """

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol_bare = Chem.MolFromSmiles(smiles_str)
        mol_full = Chem.AddHs(mol_bare)
        AllChem.EmbedMolecule(mol_full)
        if optimize:
            AllChem.UFFOptimizeMolecule(mol_full)

        if no_hydrogen:
            return Chem.RemoveHs(mol_full)
        else:
            return Chem.MolToXYZBlock(mol_full)

    except ImportError:
        raise ImportError('Unable to import rdkit.')


@staticmethod
def _Molecule_draw_2d_svg(smiles_str, width=300, height=300):
    """
    Draw 2D representation for SMILES string.

    :param smiles_str:
        The SMILES string.
    :param width:
        The width of the drawing area.
    :param height:
        The height of the drawing area.
    """

    try:
        from rdkit import Chem
        from IPython.display import SVG
        from IPython.display import display

        mol_no_hydrogen = Molecule.smiles_to_xyz(smiles_str,
                                                 optimize=True,
                                                 no_hydrogen=True)

        drawer = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol_no_hydrogen)
        drawer.FinishDrawing()

        display(SVG(drawer.GetDrawingText()))

    except ImportError:
        raise ImportError('Unable to import rdkit.Chem and/or IPython.display.')


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
            labels.append(content[0])
            coords.append([float(x) for x in content[1:4]])

    return Molecule(labels, coords, units)


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
    mol.check_multiplicity()
    mol.check_proximity(0.1)

    return mol


def _Molecule_center_of_mass(self):
    """
    Computes center of mass of a molecule in Bohr (for backward compatibility).

    :return:
        The center of mass in Bohr.
    """

    return self.center_of_mass_in_bohr()


def _Molecule_center_of_mass_in_bohr(self):
    """
    Computes center of mass of a molecule in Bohr.

    :return:
        The center of mass in Bohr.
    """

    masses = self.masses_to_numpy()
    x_coords = self.x_to_numpy()
    y_coords = self.y_to_numpy()
    z_coords = self.z_to_numpy()

    x_center = np.sum(x_coords * masses) / np.sum(masses)
    y_center = np.sum(y_coords * masses) / np.sum(masses)
    z_center = np.sum(z_coords * masses) / np.sum(masses)

    return x_center, y_center, z_center


def _Molecule_center_of_mass_in_angstrom(self):
    """
    Computes center of mass of a molecule in Angstrom.

    :return:
        The center of mass in Angstrom.
    """

    return self.center_of_mass_in_bohr() * bohr_in_angstrom()


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
    Returns atom coordinates in Bohr (for backward compatibility).

    :return:
        A numpy array of atom coordinates (nx3) in Bohr.
    """

    return self.get_coordinates_in_bohr()


def _Molecule_get_coordinates_in_bohr(self):
    """
    Returns atom coordinates in Bohr.

    :return:
        A numpy array of atom coordinates (nx3) in Bohr.
    """

    return np.array([
        self.x_to_numpy(),
        self.y_to_numpy(),
        self.z_to_numpy(),
    ]).T.copy()


def _Molecule_get_coordinates_in_angstrom(self):
    """
    Returns atom coordinates in Angstrom.

    :return:
        A numpy array of atom coordinates (nx3) in Angstrom.
    """

    return self.get_coordinates_in_bohr() * bohr_in_angstrom()


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


def _Molecule_moments_of_inertia(self):
    """
    Calculates the moment of inertia tensor and principle axes.

    :return:
        The principle moments of inertia.
    """

    masses = self.masses_to_numpy()
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


def _Molecule_deepcopy(self, memo):
    """
    Implements deepcopy.

    :param memo:
        The memo dictionary for deepcopy.

    :return:
        A deepcopy of self.
    """

    return Molecule(self)


Molecule._get_input_keywords = _Molecule_get_input_keywords

Molecule.smiles_to_xyz = _Molecule_smiles_to_xyz
Molecule.draw_2d_svg = _Molecule_draw_2d_svg
Molecule.read_smiles = _Molecule_read_smiles
Molecule.read_molecule_string = _Molecule_read_molecule_string
Molecule.read_xyz_file = _Molecule_read_xyz_file
Molecule.read_xyz_string = _Molecule_read_xyz_string
Molecule.from_dict = _Molecule_from_dict
Molecule.center_of_mass = _Molecule_center_of_mass
Molecule.center_of_mass_in_bohr = _Molecule_center_of_mass_in_bohr
Molecule.center_of_mass_in_angstrom = _Molecule_center_of_mass_in_angstrom
Molecule.more_info = _Molecule_more_info
Molecule.get_labels = _Molecule_get_labels
Molecule.get_coordinates = _Molecule_get_coordinates
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
Molecule.__deepcopy__ = _Molecule_deepcopy

# aliases for backward compatibility
Molecule.read_xyz = _Molecule_read_xyz_file
Molecule.from_xyz_string = _Molecule_read_xyz_string
Molecule.write_xyz = _Molecule_write_xyz_file
Molecule.read_str = _Molecule_read_molecule_string
