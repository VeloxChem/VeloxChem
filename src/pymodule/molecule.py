#
#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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
import geometric

from .veloxchemlib import Molecule
from .veloxchemlib import ChemicalElement
from .veloxchemlib import bohr_in_angstroms
from .veloxchemlib import mathconst_pi


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
            coords.append([float(x) for x in content[1:4]])

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

    return Molecule.read_str(xyzstr, units='angstrom')


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

    return Molecule.read_str(xyzstr, 'angstrom')


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
    Returns atom coordinates in atomic units.

    :return:
        A numpy array of atom coordinates (nx3) in atomic units.
    """

    return np.array([
        self.x_to_numpy(),
        self.y_to_numpy(),
        self.z_to_numpy(),
    ]).T.copy()


def _Molecule_write_xyz(self, xyz_filename):
    """
    Writes molecular geometry to xyz file.

    :param xyz_filename:
        The name of the xyz file.
    """

    elem_ids = self.elem_ids_to_numpy()

    xs = self.x_to_numpy() * bohr_in_angstroms()
    ys = self.y_to_numpy() * bohr_in_angstroms()
    zs = self.z_to_numpy() * bohr_in_angstroms()

    with open(str(xyz_filename), 'w') as fh:

        fh.write(f"{self.number_of_atoms():d}\n\n")

        for elem_id, x, y, z in zip(elem_ids, xs, ys, zs):
            elem = ChemicalElement()
            elem.set_atom_type(elem_id)
            fh.write(
                f'{elem.get_name():<6s} {x:22.12f} {y:22.12f} {z:22.12f}\n')


def _Molecule_get_ic_rmsd(self, ref_mol):
    """
    Gets statistical deviation of bonds, angles and dihedral angles between
    self and reference geometry.

    :param ref_mol:
        The reference molecule (or xyz filename).
    """

    if isinstance(ref_mol, str):
        errmsg = '*** Note: invalid reference xyz file!'
    else:
        errmsg = '*** Note: invalid reference molecule!'

    if isinstance(ref_mol, str):
        if Path(ref_mol).is_file():
            ref_mol = Molecule.read_xyz(ref_mol)
        else:
            return errmsg

    if ref_mol.get_labels() != self.get_labels():
        return errmsg

    g_mol = geometric.molecule.Molecule()
    g_mol.elem = self.get_labels()
    g_mol.xyzs = [self.get_coordinates() * geometric.nifty.bohr2ang]

    ic = geometric.internal.DelocalizedInternalCoordinates(g_mol, build=True)

    ref_geom = ref_mol.get_coordinates() * geometric.nifty.bohr2ang
    opt_geom = self.get_coordinates() * geometric.nifty.bohr2ang

    bonds = []
    angles = []
    dihedrals = []

    for internal in ic.Prims.Internals:
        if isinstance(internal, geometric.internal.Distance):
            v1 = internal.value(ref_geom)
            v2 = internal.value(opt_geom)
            bonds.append(abs(v1 - v2))
        elif isinstance(internal, geometric.internal.Angle):
            v1 = internal.value(ref_geom)
            v2 = internal.value(opt_geom)
            angles.append(abs(v1 - v2) * 180.0 / mathconst_pi())
        elif isinstance(internal, geometric.internal.Dihedral):
            v1 = internal.value(ref_geom)
            v2 = internal.value(opt_geom)
            diff_in_deg = (v1 - v2) * 180.0 / mathconst_pi()
            if diff_in_deg > 180.0:
                diff_in_deg -= 360.0
            elif diff_in_deg < -180.0:
                diff_in_deg += 360.0
            dihedrals.append(abs(diff_in_deg))

    ic_rmsd = {'bonds': None, 'angles': None, 'dihedrals': None}

    if bonds:
        np_bonds = np.array(bonds)
        rms_bonds = np.sqrt(np.mean(np_bonds**2))
        max_bonds = np.max(np_bonds)
        ic_rmsd['bonds'] = {
            'rms': rms_bonds,
            'max': max_bonds,
            'unit': 'Angstrom'
        }

    if angles:
        np_angles = np.array(angles)
        rms_angles = np.sqrt(np.mean(np_angles**2))
        max_angles = np.max(np_angles)
        ic_rmsd['angles'] = {
            'rms': rms_angles,
            'max': max_angles,
            'unit': 'degree'
        }

    if dihedrals:
        np_dihedrals = np.array(dihedrals)
        rms_dihedrals = np.sqrt(np.mean(np_dihedrals**2))
        max_dihedrals = np.max(np_dihedrals)
        ic_rmsd['dihedrals'] = {
            'rms': rms_dihedrals,
            'max': max_dihedrals,
            'unit': 'degree'
        }

    return ic_rmsd


def _Molecule_moments_of_inertia(self):
    """
    Calculates the moment of inertia tensor and principle axes

    :return:
        The principle moments of inertia.
    """

    masses = self.masses_to_numpy()
    coordinates = self.get_coordinates()
    center_of_mass = np.array(self.center_of_mass())
    natm = self.number_of_atoms()

    # Coordinates in the center-of-mass frame
    coords_com = coordinates - center_of_mass[np.newaxis, :]

    # Moment of inertia tensor
    Imom = np.sum([
        masses[i] * (np.eye(3) * (np.dot(coords_com[i], coords_com[i])) -
                     np.outer(coords_com[i], coords_com[i]))
        for i in range(natm)], axis=0)

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


Molecule.read_str = _Molecule_read_str
Molecule.read_xyz = _Molecule_read_xyz
Molecule.from_xyz_string = _Molecule_from_xyz_string
Molecule.from_dict = _Molecule_from_dict
Molecule.center_of_mass = _Molecule_center_of_mass
Molecule.more_info = _Molecule_more_info
Molecule.get_labels = _Molecule_get_labels
Molecule.get_coordinates = _Molecule_get_coordinates
Molecule.get_ic_rmsd = _Molecule_get_ic_rmsd
Molecule.write_xyz = _Molecule_write_xyz
Molecule.moments_of_inertia = _Molecule_moments_of_inertia
Molecule.is_linear = _Molecule_is_linear
Molecule.get_aufbau_occupation = _Molecule_get_aufbau_occupation
