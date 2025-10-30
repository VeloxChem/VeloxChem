#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from pathlib import Path
import numpy as np
import math

from .veloxchemlib import Molecule
from .veloxchemlib import bohr_in_angstrom, mpi_master
from .outputstream import OutputStream
from .inputparser import print_keywords
from .errorhandler import assert_msg_critical, safe_arccos
from .pubchemfetcher import (get_data_from_name, get_all_conformer_IDs,
                             get_conformer_data, get_pubchem_sketcher_reference)


@staticmethod
def _Molecule_smiles_formal_charge(smiles_str):
    """
    Gets formal charge from SMILES string.

    :param smiles_str:
        The SMILES string.

    :return:
        The formal charge.
    """

    try:
        from rdkit import Chem

        mol_bare = Chem.MolFromSmiles(smiles_str)
        mol_chg = Chem.GetFormalCharge(mol_bare)

        return mol_chg

    except ImportError:
        raise ImportError('Unable to import rdkit.')


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
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol_bare = Chem.MolFromSmiles(smiles_str)
        mol_full = Chem.AddHs(mol_bare)

        use_random_coords = (mol_full.GetNumConformers() == 0)
        AllChem.EmbedMolecule(mol_full, useRandomCoords=use_random_coords)

        if optimize:
            AllChem.UFFOptimizeMolecule(mol_full)

        if hydrogen:
            return Chem.MolToXYZBlock(mol_full)
        else:
            return Chem.RemoveHs(mol_full)

    except ImportError:
        raise ImportError('Unable to import rdkit.')


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
    mol = Molecule.read_xyz_string(xyz)

    mol_chg = Molecule.smiles_formal_charge(smiles_str)
    mol.set_charge(mol_chg)

    return mol


def _element_guesser(atom_name, residue_name):
    """
    Guesses the chemical element of an atom based on its name. Needed by GRO
    and PDB file format.

    :param atom_name:
        The atom name.
    :param residue_name:
        The residue name.
    """

    periodic_table = [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm',
        'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
        'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
        'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
        'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
        'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
    ]

    two_letter_elements = [elem for elem in periodic_table if len(elem) == 2]

    ion_residues = {
        'IB+': 'I',
        'CA': 'Ca',
        'CL': 'Cl',
        'NA': 'Na',
        'MG': 'Mg',
        'K': 'K',
        'RB': 'Rb',
        'CS': 'Cs',
        'LI': 'Li',
        'ZN': 'Zn'
    }

    protein_residues = [
        'URE', 'ACE', 'NME', 'NHE', 'NH2', 'ALA', 'GLY', 'SER', 'THR', 'LEU',
        'ILE', 'VAL', 'ASN', 'GLN', 'ARG', 'HID', 'HIE', 'HIP', 'TRP', 'PHE',
        'TYR', 'GLU', 'ASP', 'LYS', 'ORN', 'DAB', 'LYN', 'PRO', 'HYP', 'CYS',
        'CYM', 'CYX', 'MET', 'ASH', 'GLH', 'CALA', 'CGLY', 'CSER', 'CTHR',
        'CLEU', 'CILE', 'CVAL', 'CASN', 'CGLN', 'CARG', 'CHID', 'CHIE', 'CHIP',
        'CTRP', 'CPHE', 'CTYR', 'CGLU', 'CASP', 'CLYS', 'CPRO', 'CCYS', 'CCYX',
        'CMET', 'NALA', 'NGLY', 'NSER', 'NTHR', 'NLEU', 'NILE', 'NVAL', 'NASN',
        'NGLN', 'NARG', 'NHID', 'NHIE', 'NHIP', 'NTRP', 'NPHE', 'NTYR', 'NGLU',
        'NASP', 'NLYS', 'NORN', 'NDAB', 'NPRO', 'NCYS', 'NCYX', 'NMET'
    ]

    dna_residues = [
        "DA5", "DA", "DA3", "DAN", "DT5", "DT", "DT3", "DTN", "DG5", "DG",
        "DG3", "DGN", "DC5", "DC", "DC3", "DCN"
    ]

    rna_residues = [
        "RA5", "RA", "RA3", "RAN", "RU5", "RU", "RU3", "RUN", "RG5", "RG",
        "RG3", "RGN", "RC5", "RC", "RC3", "RCN"
    ]

    standard_residues = (protein_residues + dna_residues + rna_residues)

    if residue_name in ion_residues:
        element = ion_residues[atom_name]
        if element not in periodic_table:
            raise NameError(f"Element {element} not in periodic table")
        return element

    elif residue_name in standard_residues:
        # For standard residues, take first letter as element name
        element = atom_name[0]
        if element not in periodic_table:
            raise NameError(f"Element {element} not in periodic table")
        return element

    else:
        # Take the first characters not being a digit
        name = ''
        for c in atom_name:
            if not c.isdigit():
                name += c
            else:
                break

        # Check if the guessed name is a valid element
        name = name.capitalize()
        # Special cases where two character elements are capitalized
        # Sometimes one can find OP instead of O in HETAM residues
        if str(name) not in two_letter_elements:
            name = name[0]
        if name not in periodic_table:
            raise NameError(f"Element {name} not in periodic table")
        return name


@staticmethod
def _Molecule_read_gro_file(grofile):
    """
    Reads molecule from file in GRO format.

    :param grofile:
        File with molecular structure in GRO format.

    :return:
        The molecule.
    """

    with Path(grofile).open('r') as fh:
        grostr = fh.read()

    coordinates = []
    labels = []

    lines = grostr.strip().splitlines()

    for line in lines[2:-1]:
        if line:
            # To access the content we will stick to the C format:
            # %5d%-5s%5s%5d%8.3f%8.3f%8.3f

            residue_name = line[5:10].strip()
            atom_name = line[10:15].strip()

            # coordinates in angstroms
            x = float(line[20:28].strip()) * 10.0
            y = float(line[28:36].strip()) * 10.0
            z = float(line[36:44].strip()) * 10.0
            coordinates.append([x, y, z])

            # Assign label
            name = _element_guesser(atom_name, residue_name)
            labels.append(name)

    return Molecule(labels, coordinates, 'angstrom')


@staticmethod
def _Molecule_read_pdb_file(pdbfile):
    """
    Reads molecule from file in PDB format.

    :param pdbfile:
        File with molecular structure in PDB format.

    :return:
        The molecule.
    """

    with Path(pdbfile).open('r') as fh:
        pdbstr = fh.read()

    coordinates = []
    labels = []

    lines = pdbstr.strip().splitlines()
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):

            if line[76:78].strip() == '':
                atom_name = line[12:15].strip()
                residue_name = line[17:19].strip()
                # Guess element
                name = _element_guesser(atom_name, residue_name)
            else:
                name = str(line[76:78]).strip()

            labels.append(name)
            coordinates.append(
                [float(line[30:38]),
                 float(line[38:46]),
                 float(line[46:54])])

    return Molecule(labels, coordinates, 'angstrom')


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
    basis_set_labels = []

    for line in mol_str.strip().splitlines():
        if line:
            content = line.split()

            elem_name = content[0].upper()
            basis_elem_name = elem_name

            if elem_name.startswith('BQ_'):
                elem_name, basis_elem_name = elem_name.split('_')
            elif elem_name.endswith('_BQ'):
                basis_elem_name, elem_name = elem_name.split('_')

            labels.append(elem_name)
            coords.append([float(x) for x in content[1:4]])

            if len(content) > 4:
                basis_set_labels.append([content[4].upper(), basis_elem_name])
            else:
                basis_set_labels.append(['', basis_elem_name])

    coords = np.array(coords)

    return Molecule(labels, coords, units, basis_set_labels)


@staticmethod
def _Molecule_read_xyz_file(xyzfile):
    """
    Reads molecule from file in XYZ format.

    :param xyzfile:
        File with molecular structure in XYZ format.

    :return:
        The molecule.
    """

    assert_msg_critical(
        Path(xyzfile).is_file(), f'Molecule: xyzfile {xyzfile} does not exist')

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
    assert_msg_critical(mol.check_proximity(0.1), 'Molecule: Atoms too close')

    return mol


def _Molecule_get_connectivity_matrix(self, factor=1.3, H2_factor=1.8):
    """
    Gets connectivity matrix.

    :param factor:
        Scaling factor for the covalent radii to account for the bond
        threshold.

    :return:
        The connectivity matrix as a numpy array of integers.
    """

    coords_in_au = self.get_coordinates_in_bohr()
    covalent_radii_in_au = self.covalent_radii_to_numpy()

    natoms = coords_in_au.shape[0]
    connectivity_matrix = np.zeros((natoms, natoms), dtype='int32')
    labels = self.get_labels()

    for i in range(natoms):
        for j in range(i + 1, natoms):
            distance = np.linalg.norm(coords_in_au[j] - coords_in_au[i])
            threshold = (covalent_radii_in_au[i] + covalent_radii_in_au[j])

            # Special case for H2 molecule
            if labels[i] == 'H' and labels[j] == 'H':
                threshold *= H2_factor
            else:
                threshold *= factor

            if distance <= threshold:
                connectivity_matrix[i, j] = 1
                connectivity_matrix[j, i] = 1

    return connectivity_matrix


def _Molecule_find_connected_atoms(self, atom_idx, connectivity_matrix=None):
    """
    Gets all atoms indices that are (directly or indirectly) connected to a
    given atom.

    :param atom_idx:
        The index of the give atom.
    :param connectivity_matrix:
        The connectivity matrix.

    :return:
        A set containing all the atom indices that are connected to the given
        atom.
    """

    if connectivity_matrix is None:
        connectivity_matrix = self.get_connectivity_matrix()

    connected_atoms = set()
    connected_atoms.add(atom_idx)

    while True:
        more_connected_atoms = set()
        for a in connected_atoms:
            for b in range(connectivity_matrix.shape[0]):
                if (b not in connected_atoms
                        and connectivity_matrix[a, b] == 1):
                    more_connected_atoms.add(b)
        if more_connected_atoms:
            connected_atoms.update(more_connected_atoms)
        else:
            break

    return connected_atoms


def _Molecule_rotate_around_vector(self, coords, origin, vector, rotation_angle,
                                   angle_unit):
    """
    Returns coordinates after rotation around a given vector.

    :param coords:
        The coordinates.
    :param origin:
        The origin of the vector.
    :param vector:
        The vector.
    :param rotation_angle:
        The rotation angle.
    :param angle_unit:
        The unit of rotation angle.

    :return:
        The coordinates after rotation.
    """

    assert_msg_critical(angle_unit.lower() in ['degree', 'radian'],
                        'Molecule: Invalid angle unit for rotation')

    if angle_unit.lower() == 'degree':
        rotation_angle_in_radian = math.pi * rotation_angle / 180.0
    else:
        rotation_angle_in_radian = rotation_angle

    uvec = vector / np.linalg.norm(vector)

    cos_theta = math.cos(rotation_angle_in_radian)
    sin_theta = math.sin(rotation_angle_in_radian)
    m_cos_theta = 1.0 - cos_theta

    rotation_mat = np.zeros((3, 3))

    rotation_mat[0, 0] = cos_theta + m_cos_theta * uvec[0]**2
    rotation_mat[1, 1] = cos_theta + m_cos_theta * uvec[1]**2
    rotation_mat[2, 2] = cos_theta + m_cos_theta * uvec[2]**2

    rotation_mat[0, 1] = m_cos_theta * uvec[0] * uvec[1] - sin_theta * uvec[2]
    rotation_mat[1, 0] = m_cos_theta * uvec[1] * uvec[0] + sin_theta * uvec[2]

    rotation_mat[1, 2] = m_cos_theta * uvec[1] * uvec[2] - sin_theta * uvec[0]
    rotation_mat[2, 1] = m_cos_theta * uvec[2] * uvec[1] + sin_theta * uvec[0]

    rotation_mat[2, 0] = m_cos_theta * uvec[2] * uvec[0] - sin_theta * uvec[1]
    rotation_mat[0, 2] = m_cos_theta * uvec[0] * uvec[2] + sin_theta * uvec[1]

    return np.matmul(coords - origin, rotation_mat.T) + origin


def _Molecule_get_dihedral_in_degrees(self, dihedral_indices_one_based):
    """
    Gets dihedral angle.

    :param dihedral_indices_one_based:
        The dihedral indices (1-based).

    :return:
        The dihedral angle.
    """

    return self.get_dihedral(dihedral_indices_one_based, 'degree')


def _Molecule_get_dihedral(self, dihedral_indices_one_based, angle_unit):
    """
    Gets dihedral angle.

    :param dihedral_indices_one_based:
        The dihedral indices (1-based).
    :param angle_unit:
        The unit of angle (degree or radian).

    :return:
        The dihedral angle.
    """

    assert_msg_critical(
        len(dihedral_indices_one_based) == 4,
        'Molecule.get_dihedral: Expecting four atom indices (1-based)')

    a = dihedral_indices_one_based[0] - 1
    b = dihedral_indices_one_based[1] - 1
    c = dihedral_indices_one_based[2] - 1
    d = dihedral_indices_one_based[3] - 1

    coords_in_au = self.get_coordinates_in_bohr()

    # J. Comput. Chem. 2000, 21, 553-561

    v21 = coords_in_au[a] - coords_in_au[b]
    v32 = coords_in_au[b] - coords_in_au[c]
    v43 = coords_in_au[c] - coords_in_au[d]

    u21 = v21 / np.linalg.norm(v21)
    u32 = v32 / np.linalg.norm(v32)
    u43 = v43 / np.linalg.norm(v43)

    cos_theta_123 = -np.vdot(u21, u32)
    cos_theta_234 = -np.vdot(u32, u43)

    sin_theta_123 = math.sqrt(1.0 - cos_theta_123**2)
    sin_theta_234 = math.sqrt(1.0 - cos_theta_234**2)

    cos_phi = ((cos_theta_123 * cos_theta_234 - np.vdot(u21, u43)) /
               (sin_theta_123 * sin_theta_234))
    sin_phi = -(np.vdot(u43, np.cross(u21, u32)) /
                (sin_theta_123 * sin_theta_234))

    phi_in_radian = safe_arccos(cos_phi)
    if sin_phi < 0.0:
        phi_in_radian *= -1.0

    assert_msg_critical(angle_unit.lower() in ['degree', 'radian'],
                        'Molecule.get_dihedral: Invalid angle unit')

    if angle_unit.lower() == 'degree':
        return 180.0 * phi_in_radian / math.pi
    else:
        return phi_in_radian


def _Molecule_set_dihedral_in_degrees(self, dihedral_indices_one_based,
                                      target_angle):
    """
    Sets dihedral angle.

    :param dihedral_indices_one_based:
        The dihedral indices (1-based).
    :param target_angle:
        The target value of dihedral angle.
    """

    self.set_dihedral(dihedral_indices_one_based, target_angle, 'degree')


def _Molecule_set_dihedral(self, dihedral_indices_one_based, target_angle,
                           angle_unit):
    """
    Sets dihedral angle.

    :param dihedral_indices_one_based:
        The dihedral indices (1-based).
    :param target_angle:
        The target value of dihedral angle.
    :param angle_unit:
        The unit of angle (degree or radian).
    """

    assert_msg_critical(
        len(dihedral_indices_one_based) == 4,
        'Molecule.set_dihedral: Expecting four atom indices (1-based)')

    # get the 0-based atom indices for central bond
    i = dihedral_indices_one_based[1] - 1
    j = dihedral_indices_one_based[2] - 1

    # disconnect i-j and find all atoms that at connected to j
    connectivity_matrix = self.get_connectivity_matrix()
    connectivity_matrix[i, j] = 0
    connectivity_matrix[j, i] = 0

    atoms_connected_to_j = self._find_connected_atoms(j, connectivity_matrix)

    assert_msg_critical(
        i not in atoms_connected_to_j,
        'Molecule.set_dihedral: Cannot rotate dihedral ' +
        '(Maybe it is part of a ring?)')

    # rotate whole molecule around vector i->j
    coords_in_au = self.get_coordinates_in_bohr()

    vij = coords_in_au[j] - coords_in_au[i]

    current_angle = self.get_dihedral(dihedral_indices_one_based, angle_unit)

    updated_target_angle = target_angle
    half_period = 180.0 if angle_unit.lower() == 'degree' else math.pi
    while updated_target_angle - current_angle > half_period:
        updated_target_angle -= 2.0 * half_period
    while updated_target_angle - current_angle < -half_period:
        updated_target_angle += 2.0 * half_period

    # make several attempts to rotate the dihedral angle, with the constraint
    # that the connectivity matrix should not change
    for attempt in range(10, -1, -1):

        rotation_angle = (updated_target_angle - current_angle) * (0.1 *
                                                                   attempt)

        new_coords_in_au = self._rotate_around_vector(coords_in_au,
                                                      coords_in_au[j], vij,
                                                      rotation_angle,
                                                      angle_unit)

        new_mol = Molecule(self)
        for idx in atoms_connected_to_j:
            new_mol.set_atom_coordinates(idx, new_coords_in_au[idx])

        new_conn_mat = new_mol.get_connectivity_matrix()
        conn_mat = self.get_connectivity_matrix()
        if np.max(np.abs(new_conn_mat - conn_mat)) < 1.0e-10:
            for idx in atoms_connected_to_j:
                self.set_atom_coordinates(idx, new_coords_in_au[idx])
            return

    assert_msg_critical(
        False, 'Molecule.set_dihedral: Cannot set dihedral angle due to ' +
        'overlapping atoms')


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

    mol_str = 'Molecular Geometry (Angstroms)\n'
    mol_str += '================================\n\n'
    mol_str += '  Atom'
    mol_str += '         Coordinate X '
    mol_str += '         Coordinate Y '
    mol_str += '         Coordinate Z  \n\n'

    xyz_lines = self.get_xyz_string().splitlines()

    for line in xyz_lines[2:]:
        content = line.split()

        elem_name = content[0]
        if elem_name.startswith('Bq_') or elem_name.endswith('_Bq'):
            elem_name = ' ' + elem_name
        else:
            elem_name = '  ' + elem_name

        x_coord = float(content[1])
        y_coord = float(content[2])
        z_coord = float(content[3])

        mol_str += f'{elem_name:<6s}'
        mol_str += f'{x_coord:>22.12f}'
        mol_str += f'{y_coord:>22.12f}'
        mol_str += f'{z_coord:>22.12f}'
        mol_str += '\n'

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
    for r in self._get_coordinates():
        coords.append(r.coordinates())

    return np.array(coords)


def _Molecule_get_coordinates_in_angstrom(self):
    """
    Returns atom coordinates in Angstrom.

    :return:
        A numpy array of atom coordinates (nx3) in Angstrom.
    """

    coords = []
    for r in self._get_coordinates('angstrom'):
        coords.append(r.coordinates())

    return np.array(coords)


def _Molecule_get_distance_matrix_in_angstrom(self):
    """
    Returns distance matrix in Angstrom.

    :return:
        A numpy array of distance matrix (nxn) in Angstrom.
    """

    coords = self.get_coordinates_in_angstrom()
    natoms = coords.shape[0]
    distance_matrix = np.zeros((natoms, natoms))

    for i in range(natoms):
        for j in range(i, natoms):
            rij = np.linalg.norm(coords[i, :] - coords[j, :])
            distance_matrix[i, j] = rij
            if i != j:
                distance_matrix[j, i] = rij

    return distance_matrix


def _Molecule_get_xyz_string(self, precision=12, comment=''):
    """
    Returns xyz string of molecule.

    :return:
        An xyz string (including number of atoms).
    """

    labels = self.get_labels()
    coords_in_angstrom = self.get_coordinates_in_angstrom()

    elem_ids = self.get_identifiers()
    atom_basis_labels = self.get_atom_basis_labels()

    if comment and isinstance(comment, str):
        trimmed_comment = comment.strip().splitlines()[0][:80]
    else:
        trimmed_comment = ''

    natoms = len(labels)
    xyz = f'{natoms}\n{trimmed_comment}\n'

    for a in range(natoms):
        xa, ya, za = coords_in_angstrom[a]

        if elem_ids[a] == 0 and atom_basis_labels[a][1]:
            if atom_basis_labels[a][1].capitalize() == 'Bq':
                elem_name = 'Bq'
            else:
                elem_name = 'Bq_' + atom_basis_labels[a][1].capitalize()
        else:
            elem_name = labels[a]

        xyz += f'{elem_name:<6s}'
        xyz += f' {xa:{precision + 10}.{precision}f}'
        xyz += f' {ya:{precision + 10}.{precision}f}'
        xyz += f' {za:{precision + 10}.{precision}f}\n'

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
                   atom_labels=False,
                   starting_index=1,
                   bonds=None,
                   dashed_bonds=None):
    """
    Creates a 3D view with py3dmol.

    :param width:
        The width.
    :param height:
        The height.
    :param atom_indices:
        The flag for showing atom indices (1-based).
    :param atom_labels:
        The flag for showing atom labels. If provided with a list, will use
        that list as labels.
    :starting_index:
        The starting index for atom indices.
    :bonds:
        A list of tuples with bonds to draw. If None, connectivity is based on
        proximity.
    :dashed_bonds:
        A list of tuples with bonds to draw as dashed lines.
    """

    try:
        import py3Dmol
        viewer = py3Dmol.view(width=width, height=height)

        if bonds is None:
            viewer.addModel(self.get_xyz_string())
        else:
            from rdkit import Chem
            import re

            rdmol = Chem.MolFromXYZBlock(self.get_xyz_string())
            edit_mol = Chem.EditableMol(rdmol)

            for bond in bonds:
                if dashed_bonds is not None:
                    if bond in dashed_bonds or (bond[1],
                                                bond[0]) in dashed_bonds:
                        continue
                edit_mol.AddBond(bond[0], bond[1], Chem.BondType.SINGLE)

            sdf = Chem.MolToMolBlock(edit_mol.GetMol())

            if dashed_bonds is not None:
                lines = sdf.split('\n')
                last_line = lines[-2]
                lines = lines[:-2]
                for key in sorted(dashed_bonds):
                    line = f" {key[0] + 1:<2} {key[1] + 1:<2} 0.5  0"
                    lines.append(line)
                lines.append(last_line)

                splitline = re.split(r'(\s+)', lines[3])
                splitline[4] = str(int(splitline[4]) + len(dashed_bonds))
                splitline[3] = " " * (2 - len(splitline[4]))
                lines[3] = ''.join(splitline)
                sdf = '\n'.join(lines)

            viewer.addModel(sdf, 'sdf')

        if atom_indices or atom_labels:
            coords = self.get_coordinates_in_angstrom()
            if isinstance(atom_labels, list):
                labels = atom_labels
                assert_msg_critical(
                    len(labels) == coords.shape[0],
                    'Molecule.show: Inconsistent number of atom_labels')
            else:
                labels = self.get_labels()
            for i in range(coords.shape[0]):
                text = ''
                if atom_labels:
                    text += f'{labels[i]}'
                if atom_indices:
                    text += f'{i + starting_index}'
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
        viewer.setViewStyle({"style": "outline", "width": 0.05})
        viewer.setStyle({"stick": {}, "sphere": {"scale": 0.25}})
        viewer.zoomTo()
        viewer.show()

    except ImportError:
        raise ImportError('Unable to import py3Dmol or rdkit')


@staticmethod
def _Molecule_draw_2d(smiles_str, width=400, height=300):
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
                                                 hydrogen=False)

        drawer = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol_no_hydrogen)
        drawer.FinishDrawing()

        display(SVG(drawer.GetDrawingText()))

    except ImportError:
        raise ImportError('Unable to import rdkit.Chem and/or IPython.display.')


def _Molecule_moments_of_inertia(self, principal_axes=False):
    """
    Calculates the moment of inertia tensor and principal axes.

    :return:
        The principal moments of inertia and principal axes.
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

    if principal_axes:
        # Note: Eigenvectors are in the rows after transpose
        return Ivals, Ivecs.T
    else:
        return Ivals


def _Molecule_is_linear(self):
    """
    Checks if a molecule is linear or not.

    :return:
        True if linear, False otherwise.
    """

    assert_msg_critical(self.number_of_atoms() >= 2,
                        'Molecule.is_linear: Need at least two atoms')

    # Get principal moments of inertia
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


def _Molecule_partition_atoms(self, comm):
    """
    Partition atoms for a given MPI communicator.

    :param comm:
        The MPI communicator.

    :return:
        The list of atom indices for the current MPI rank.
    """

    rank = comm.Get_rank()
    nnodes = comm.Get_size()

    if rank == mpi_master():
        elem_ids = self.get_identifiers()
        coords = self.get_coordinates_in_bohr()
        mol_com = self.center_of_mass_in_bohr()

        r2_array = np.sum((coords - mol_com)**2, axis=1)
        sorted_r2_list = sorted([
            (r2, nchg, i)
            for i, (r2, nchg) in enumerate(zip(r2_array, elem_ids))
        ])

        dict_atoms = {}
        for r2, nchg, i in sorted_r2_list:
            if nchg not in dict_atoms:
                dict_atoms[nchg] = []
            dict_atoms[nchg].append(i)

        list_atoms = []
        for nchg in sorted(dict_atoms.keys(), reverse=True):
            list_atoms += dict_atoms[nchg]

    else:
        list_atoms = None

    list_atoms = comm.bcast(list_atoms, root=mpi_master())

    return list(list_atoms[rank::nnodes])


def _Molecule_is_water_molecule(self):
    """
    Checks if a molecule is a water molecule.

    :return:
        True if the molecule is a water molecule, False otherwise.
    """

    natoms = self.number_of_atoms()
    if natoms != 3:
        return False

    labels = self.get_labels()
    if sorted(labels) != ['H', 'H', 'O']:
        return False

    conn = self.get_connectivity_matrix()

    bond_labels = [
        sorted([labels[i], labels[j]]) for i in range(natoms)
        for j in range(i, natoms) if conn[i, j] == 1
    ]

    if bond_labels != [['H', 'O'], ['H', 'O']]:
        return False

    return True


def _Molecule_contains_water_molecule(self):
    """
    Checks if a molecule contains a water molecule.

    :return:
        True if the molecule contains a water molecule, False otherwise.
    """

    conn = self.get_connectivity_matrix()
    natoms = self.number_of_atoms()
    visited = set()
    labels = self.get_labels()

    for i in range(natoms):
        if i in visited:
            continue
        if labels[i] != 'O':
            continue

        connected_atoms = self._find_connected_atoms(i, conn)
        if len(connected_atoms) != 3:
            continue
        
        connected_atoms.remove(i)
        connected_atoms = list(connected_atoms)
        if labels[connected_atoms[0]] == 'H' and labels[
                connected_atoms[1]] == 'H':
            return True

        visited.update(connected_atoms)

    return False


@staticmethod
def _Molecule_read_name(mol_name):
    """
    Reads molecule from its name as a string using PubChem data.

    :param mol_name:
        The molecule name string.

    :return:
        The molecule.
    """

    smiles_str = Molecule.name_to_smiles(mol_name)
    mol = Molecule.read_smiles(smiles_str)

    return mol


@staticmethod
def _Molecule_name_to_smiles(mol_name):
    """
    Returns SMILES-string for a given molecule name using PubChem data.

    :param mol_name:
        The molecule name string.

    :return smiles_str:
        The SMILES-string of the molecule
    """

    smiles_str, title, cid = get_data_from_name(mol_name)

    return smiles_str


@staticmethod
def _Molecule_name_to_xyz(mol_name):
    """
    Returns xyz-string for the first conformer
    of a given molecule name using Pubchem data.

    :param mol_name:
        The molecule name string.

    :return xyz:
        xyz-string of molecule
    """

    conformerID_list = get_all_conformer_IDs(mol_name)
    conformer_ID = conformerID_list[0]
    xyz = get_conformer_data(conformer_ID)

    return xyz


@staticmethod
def _Molecule_get_all_conformer_data(mol_name):
    """
    Gets coordinates for all the conformers of a given compound
    Data from PubChem.

    :param mol_name:
        The molecule name string.

    :return all_conformer_info:
        The coordinates for all of the conformers of the molecule.
    """

    conformerID_list = get_all_conformer_IDs(mol_name)
    all_conformer_info = {}
    for key in conformerID_list:
        all_conformer_info.update({key: get_conformer_data(key)})
    return all_conformer_info


@staticmethod
def _Molecule_builder():
    """
    Displays the PubChem molecule builder.
    """

    try:
        from IPython.display import display, IFrame

        print("The molecule builder uses the PubChem Sketcher.\n")
        print("Direct link to the PubChem Sketcher: " +
              "https://pubchem.ncbi.nlm.nih.gov/edit/ \n")
        print("Documentation for the PubChem Sketcher: " +
              "https://pubchem.ncbi.nlm.nih.gov/docs/sketcher-help \n")
        print("Reference: " + get_pubchem_sketcher_reference() + "\n")

        display(
            IFrame('https://pubchem.ncbi.nlm.nih.gov/edit/',
                   width=900,
                   height=450))

    except ImportError:
        raise ImportError('Unable to import IPython.')


Molecule._get_input_keywords = _Molecule_get_input_keywords
Molecule._find_connected_atoms = _Molecule_find_connected_atoms
Molecule._rotate_around_vector = _Molecule_rotate_around_vector

Molecule.smiles_formal_charge = _Molecule_smiles_formal_charge
Molecule.smiles_to_xyz = _Molecule_smiles_to_xyz
Molecule.read_gro_file = _Molecule_read_gro_file
Molecule.read_pdb_file = _Molecule_read_pdb_file
Molecule.show = _Molecule_show
Molecule.draw_2d = _Molecule_draw_2d
Molecule.read_smiles = _Molecule_read_smiles
Molecule.read_molecule_string = _Molecule_read_molecule_string
Molecule.read_xyz_file = _Molecule_read_xyz_file
Molecule.read_xyz_string = _Molecule_read_xyz_string
Molecule.from_dict = _Molecule_from_dict
Molecule.get_connectivity_matrix = _Molecule_get_connectivity_matrix
Molecule.get_dihedral = _Molecule_get_dihedral
Molecule.set_dihedral = _Molecule_set_dihedral
Molecule.get_dihedral_in_degrees = _Molecule_get_dihedral_in_degrees
Molecule.set_dihedral_in_degrees = _Molecule_set_dihedral_in_degrees
Molecule.center_of_mass_in_bohr = _Molecule_center_of_mass_in_bohr
Molecule.center_of_mass_in_angstrom = _Molecule_center_of_mass_in_angstrom
Molecule.get_string = _Molecule_get_string
Molecule.more_info = _Molecule_more_info
Molecule.get_coordinates_in_bohr = _Molecule_get_coordinates_in_bohr
Molecule.get_coordinates_in_angstrom = _Molecule_get_coordinates_in_angstrom
Molecule.get_distance_matrix_in_angstrom = _Molecule_get_distance_matrix_in_angstrom
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
Molecule.partition_atoms = _Molecule_partition_atoms
Molecule.is_water_molecule = _Molecule_is_water_molecule
Molecule.contains_water_molecule = _Molecule_contains_water_molecule

Molecule.read_name = _Molecule_read_name
Molecule.name_to_smiles = _Molecule_name_to_smiles
Molecule.name_to_xyz = _Molecule_name_to_xyz
Molecule.get_all_conformer_data = _Molecule_get_all_conformer_data
Molecule.builder = _Molecule_builder

# aliases for backward compatibility
Molecule.read_xyz = _Molecule_read_xyz_file
Molecule.from_xyz_string = _Molecule_read_xyz_string
Molecule.write_xyz = _Molecule_write_xyz_file
Molecule.read_str = _Molecule_read_molecule_string
