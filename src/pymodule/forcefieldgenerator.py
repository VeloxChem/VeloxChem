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
from scipy.optimize import curve_fit
import numpy as np
import shutil
import re

from .molecule import Molecule
from .veloxchemlib import bohr_in_angstroms
from .veloxchemlib import hartree_in_kcalpermol
from .inputparser import parse_input
from .errorhandler import assert_msg_critical


class ForceFieldGenerator:
    """
    Parameterizes Amber force fields and creates GROMACS topologies.

    Instance variables:
        - molecule: The molecule.
        - molecule_name: The name of the molecule.
        - eq_param: If equilibrium bond lengths and angles should be used.
        - r_thresh: The threshold for warning if bond lenghts deviate (nm).
        - add_impropers: If impropers avoiding inversion should be added.
        - gen_pairs: If 1-4 parameters should be generated.
        - fudgeLJ: The factor by which to multiply Lennard-Jones 1-4 interactions.
        - fudgeQQ: The factor by which to multiply electrostatic 1-4 interactions.
        - comb_rule: The number of the combination rule.
        - nbfunc: The non-bonded function type (1: Lennard-Jones, 2: Buckingham).
        - force_field_data: The filename of the data file with force field parameters.
        - data_extension: The filename of the data file with user specified parameters.
        - box_length: The box length for energy minimizations with OpenMM (Angstroms).
        - temperature: The temperature (Kelvin).
        - path_to_required_gromacs_files: The path to share/gromacs/top of the
          GROMACS directory.
        - kfac: The dihedral force constant to restrain dihedrals.
        - remove_restraint_energy: If the restraint function energy should be removed.
        - scan_angles: The angles for MM dihedral scans.
        - dft_energies: The energies of relaxed dihedral DFT scans.
        - scan_geometries: The DFT optimized geometries for the angles.
        - dihedrals: The dihedrals.
    """

    def __init__(self, molecule, dft_scans=None):
        """
        Initializes force field generator.
        """

        # topology settings
        self.molecule = molecule
        self.molecule_name = 'molecule'
        self.eq_param = True
        self.r_thresh = 0.005
        self.theta_thresh = 10.
        self.add_impropers = False
        self.gen_pairs = True
        self.fudgeLJ = 0.5
        self.fudgeQQ = 0.8333
        self.comb_rule = 2
        self.nbfunc = 1
        self.nrexcl = 3
        self.force_field_data = None
        self.data_extension = None

        # mm settings
        self.box_length = 100.0
        self.temperature = 293.15
        self.path_to_required_gromacs_files = None
        self.kfac = 10000.0
        self.remove_restraint_energy = False

        # scan settings
        self.scan_angles = []
        self.dft_energies = []
        self.scan_geometries = []
        self.dihedrals = []

        # reading DFT data

        if dft_scans is not None:
            try:
                import MDAnalysis as mda
            except ImportError:
                raise ImportError(
                    'Unable to import MDAnalysis. Please install ' +
                    'MDAnalysis via \'conda install MDAnalysis\'')

            for xyz_file in dft_scans:
                geometries = []
                u = mda.Universe(xyz_file)
                for ts in u.trajectory:
                    geometries.append(
                        Molecule(u.atoms.names, u.atoms.positions, 'angstrom'))
                self.scan_geometries.append(geometries)
                with open(xyz_file, 'rt') as xyz:
                    e = []
                    phi = []
                    pattern = re.compile('Scan')
                    for line in xyz:
                        if re.search(pattern, line):
                            strings = line.split()
                            e.append(float(strings[-1]))
                            phi.append(float(strings[7]))
                            dih = [int(i) for i in strings[5].split('-')]
                    self.dft_energies.append(e)
                    self.scan_angles.append(phi)
                    self.dihedrals.append(dih)

    def update_settings(self, ffg_dict):
        """
        Updates settings in force field generator.

        :param ffg_dict:
            The input dictionary of force field group.
        """

        ffg_keywords = {
            'molecule_name': 'str',
            'eq_param': 'bool',
            'r_thresh': 'float',
            'theta_thresh': 'float',
            'add_impropers': 'bool',
            'gen_pairs': 'bool',
            'fudgeLJ': 'float',
            'fudgeQQ': 'float',
            'comb_rule': 'int',
            'nbfunc': 'int',
            'box_length': 'float',
            'temperature': 'float',
            'path_to_required_gromacs_files': 'str',
            'scan_angles': 'seq_fixed',
            'dft_energies': 'seq_fixed',
            'scan_geometries': 'seq_fixed',
            'dihedrals': 'seq_fixed',
            'kfac': 'float',
            'remove_restraint_energy': 'bool',
            'force_field_data': 'str',
            'data_extension': 'str',
        }

        parse_input(self, ffg_keywords, ffg_dict)

        if self.data_extension is None:
            self.data_extension = '{}_extension.{}'.format(
                self.force_field_data.split('.')[0],
                self.force_field_data.split('.')[1])

    def write_top(self, filename, itp_file):
        """
        Writes a topology file.

        :param filename:
            The filename of the topology file.
        :param itp_file:
            The filename of the included itp file.
        """

        with open(filename, 'w') as top_file:

            # header

            top_file.write('; Generated by VeloxChem\n')

            # defaults

            top_file.write('\n[ defaults ]\n')
            cur_str = '; nbfunc        comb-rule       gen-pairs'
            cur_str += '       fudgeLJ fudgeQQ\n'
            top_file.write(cur_str)
            gen_pairs = 'yes' if self.gen_pairs else 'no'
            top_file.write('{}{:16}{:>18}{:19.4f}{:8.4f}\n'.format(
                self.nbfunc, self.comb_rule, gen_pairs, self.fudgeLJ,
                self.fudgeQQ))

            # include itp

            top_file.write('\n#include "' + itp_file + '"\n')

            # system

            top_file.write('\n[ system ]\n')
            top_file.write(' {}\n'.format(self.molecule_name))

            # molecules

            top_file.write('\n[ molecules ]\n')
            top_file.write('; Compound        nmols\n')
            top_file.write('{:>10}{:9}\n'.format(self.molecule_name, 1))

            top_file.close()

    def write_original_itp(self, filename, atom_types, charges):
        """
        Writes an itp file with the original parameters.

        :param filename:
            The filename.
        :param atom_types:
            The atom types.
        :param charges:
            The charges.
        """

        # molecular information

        coords = self.molecule.get_coordinates()
        n_atoms = self.molecule.number_of_atoms()
        con = self.get_connectivity()

        # preparing atom types and atom names

        assert_msg_critical(
            len(atom_types) == n_atoms,
            'ForceFieldGenerator: inconsistent atom_types')

        for i in range(n_atoms):
            atom_types[i] = f'{atom_types[i].strip():<2s}'
        unique_atom_types = list(set(atom_types))

        atom_names = self.get_atom_names()

        # preparing 1-4 pairs

        bond_indices = set()
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                if con[i, j]:
                    bond_indices.add((i, j))
        bond_indices = sorted(list(bond_indices))

        angle_indices = set()
        for i, j in bond_indices:
            for k in range(n_atoms):
                if k in [i, j]:
                    continue
                if con[j, k]:
                    inds = (i, j, k) if i < k else (k, j, i)
                    angle_indices.add(inds)
                if con[k, i]:
                    inds = (k, i, j) if k < j else (j, i, k)
                    angle_indices.add(inds)

        dihedral_indices = set()
        for i, j, k in angle_indices:
            for l in range(n_atoms):
                if l in [i, j, k]:
                    continue
                if con[k, l]:
                    inds = (i, j, k, l) if i < l else (l, k, j, i)
                    dihedral_indices.add(inds)
                if con[l, i]:
                    inds = (l, i, j, k) if l < k else (k, j, i, l)
                    dihedral_indices.add(inds)

        exclusion_indices = []
        if self.nrexcl >= 2:
            for i, j in bond_indices:
                exclusion_indices.append((i, j))
        if self.nrexcl >= 3:
            for i, j, k in angle_indices:
                exclusion_indices.append((i, k))

        pairs_14 = []
        for i, j, k, l in dihedral_indices:
            if (i, l) not in exclusion_indices:
                pairs_14.append((i, l))
        pairs_14.sort()

        with open(self.force_field_data, 'r') as ff_data:
            ff_data_lines = ff_data.readlines()

        with open(filename, 'w') as itp_file:

            # header

            itp_file.write('; Generated by VeloxChem\n')

            # atom types

            itp_file.write('\n[ atomtypes ]\n')
            cur_str = ';name   bond_type     mass     charge'
            cur_str += '   ptype   sigma         epsilon\n'
            itp_file.write(cur_str)

            for at in unique_atom_types:
                atom_type_found = False

                for line in ff_data_lines:
                    if line.startswith(f'  {at}  '):
                        cur_str = '{:>3}{:>9}{:17.5f}{:9.5f}{:>4}'.format(
                            at, at, 0., 0., 'A')
                        cur_str += '{:16.5e}{:14.5e}\n'.format(
                            float(line.split()[1]) * 2**(-1 / 6) * 2 / 10,
                            float(line.split()[2]) * 4.184)
                        itp_file.write(cur_str)
                        atom_type_found = True
                        break

                assert_msg_critical(
                    atom_type_found,
                    f'ForceFieldGenerator: Unknown atom type {at}')

            # molecule type

            itp_file.write('\n[ moleculetype ]\n')
            itp_file.write(';name            nrexcl\n')
            itp_file.write('{:>10}{:10}\n'.format(self.molecule_name,
                                                  self.nrexcl))

            # atoms

            itp_file.write('\n[ atoms ]\n')
            cur_str = ';   nr  type  resi  res  atom  cgnr'
            cur_str += '     charge      mass\n'
            itp_file.write(cur_str)

            for i, at in enumerate(atom_types):
                for line in ff_data_lines:
                    if line.startswith(f'{at} '):
                        mass_i = float(line.split()[1])
                        cur_str = '{:6}{:>5}{:6}{:>6}{:>6}'.format(
                            i + 1, atom_types[i], 1, 'RES', atom_names[i])
                        cur_str += '{:5}{:13.6f}{:13.5f}\n'.format(
                            i + 1, charges[i], mass_i)
                        itp_file.write(cur_str)

            # bonds

            itp_file.write('\n[ bonds ]\n')
            itp_file.write(';    i      j    funct       r           k_r\n')

            for i, j in bond_indices:
                r_eq = np.linalg.norm(coords[i] - coords[j])
                r_eq *= bohr_in_angstroms() * 0.1

                at_1 = atom_types[i]
                at_2 = atom_types[j]

                str_1 = r'\A' + f'{at_1}-{at_2}  '
                str_2 = r'\A' + f'{at_2}-{at_1}  '
                pattern_1 = re.compile(str_1)
                pattern_2 = re.compile(str_2)

                bond_found = False

                for line in ff_data_lines:
                    if (re.search(pattern_1, line) or
                            re.search(pattern_2, line)):
                        bond_ff = line[5:].strip().split()
                        r = float(bond_ff[1]) * 0.1
                        k_r = float(bond_ff[0]) * 4.184 * 2 * 100
                        bond_found = True
                        break

                errmsg = f'ForceFieldGenerator: bond {at_1}-{at_2}'
                errmsg += ' is not available'
                assert_msg_critical(bond_found, errmsg)

                if abs(r - r_eq) > self.r_thresh:
                    msg = 'Warning: Length of bond {}-{}'.format(at_1, at_2)
                    msg += ' in data does not match length in XYZ file'
                    msg += ' (atoms no. {} and {}):'.format(i + 1, j + 1)
                    msg += ' {:.3f} vs. {:.3f}!'.format(r[-1], r_eq[-1])
                    msg += ' Check atom types!'
                    print(msg)

                if self.eq_param:
                    r = r_eq

                cur_str = '{:6}{:7}{:7}{:14.4e}{:14.4e} ;{:>7} -{:>3}\n'.format(
                    i + 1, j + 1, 1, r, k_r, atom_names[i], atom_names[j])
                itp_file.write(cur_str)

            # pairs

            if self.gen_pairs:
                itp_file.write('\n[ pairs ]\n')
                itp_file.write(';    i      j    funct\n')

                for pair in pairs_14:
                    itp_file.write('{:6}{:7}{:7}\n'.format(
                        pair[0] + 1, pair[1] + 1, 1))

            # angles

            itp_file.write('\n[ angles ]\n')
            cur_str = ';    i      j      k    funct'
            cur_str += '     theta       k_theta\n'
            itp_file.write(cur_str)

            theta = []
            theta_eq = []
            k_theta = []
            for i, at1 in enumerate(atom_types):
                for j, at2 in enumerate(atom_types):
                    if con[i, j]:
                        for k, at3 in enumerate(atom_types[i + 1:],
                                                start=i + 1):
                            if con[j, k]:
                                a = coords[i] - coords[j]
                                b = coords[k] - coords[j]
                                theta_eq.append(
                                    np.arccos(
                                        np.dot(a, b) / np.linalg.norm(a) /
                                        np.linalg.norm(b)) * 180 / np.pi)
                                str1 = r'\A' + f'{at1}-{at2}-{at3} '
                                str2 = r'\A' + f'{at3}-{at2}-{at1} '
                                pattern1 = re.compile(str1)
                                pattern2 = re.compile(str2)

                                theta_available = False
                                ff_data = open(self.force_field_data, 'rt')

                                for line in ff_data:
                                    if re.search(pattern1, line):
                                        n_spaces = at1.count(' ') + at2.count(
                                            ' ')
                                        theta.append(
                                            float(line.split()[2 + n_spaces]))
                                        k_theta.append(
                                            float(line.split()[1 + n_spaces]) *
                                            4.184 * 2)
                                        theta_available = True
                                    elif re.search(pattern2, line):
                                        n_spaces = at3.count(' ') + at2.count(
                                            ' ')
                                        theta.append(
                                            float(line.split()[2 + n_spaces]))
                                        k_theta.append(
                                            float(line.split()[1 + n_spaces]) *
                                            4.184 * 2)
                                        theta_available = True

                                errmsg = 'ForceFieldGenerator: Angle {}-{}-{}'.format(
                                    at1, at2, at3)
                                errmsg += ' is not available in data! Check atom types!'
                                assert_msg_critical(theta_available, errmsg)

                                if abs(theta[-1] -
                                       theta_eq[-1]) > self.theta_thresh:
                                    msg = 'Warning: Angle {}-{}-{}'.format(
                                        at1, at2, at3)
                                    msg += ' in data does not match angle in XYZ file'
                                    msg += ' (atoms no. {}, {} and {}):'.format(
                                        i + 1, j + 1, k + 1)
                                    msg += ' {:.1f} vs. {:.1f}!'.format(
                                        theta[-1], theta_eq[-1])
                                    msg += ' Check atom types!'
                                    print(msg)

            if self.eq_param:
                theta = theta_eq

            n = 0
            for i in range(n_atoms):
                for j in range(n_atoms):
                    if con[i, j]:
                        for k in range(i + 1, n_atoms):
                            if con[j, k]:
                                itp_file.write(
                                    '{:6}{:7}{:7}{:7}{:14.4e}{:14.4e} ;{:>7} -{:>3} -{:>3}\n'
                                    .format(i + 1, j + 1, k + 1, 1, theta[n],
                                            k_theta[n], atom_names[i],
                                            atom_names[j], atom_names[k]))
                                n += 1

            # proper dihedrals

            itp_file.write('\n[ dihedrals ] ; propers\n')
            cur_str = ';    i      j      k      l    funct'
            cur_str += '    phase     k_d      n\n'
            itp_file.write(cur_str)
            cur_str = ';                                        '
            cur_str += 'C0         C1         C2         C3         C4         C5\n'
            itp_file.write(cur_str)

            dihedrals = []
            ats = []
            for i, at1 in enumerate(atom_types):
                for j, at2 in enumerate(atom_types):
                    if con[i, j]:
                        for k, at3 in enumerate(atom_types):
                            if con[j, k] and i != k:
                                for l, at4 in enumerate(atom_types[i + 1:],
                                                        start=i + 1):
                                    if con[k, l] and j != l:
                                        dihedrals.append([i, j, k, l])
                                        ats.append([at1, at2, at3, at4])

            for ([i, j, k, l], [at1, at2, at3, at4]) in zip(dihedrals, ats):
                dihedral_available = False

                str1 = '{}-{}-{}-{}'.format(at1, at2, at3, at4)
                str2 = '{}-{}-{}-{}'.format(at4, at3, at2, at1)
                pattern1 = re.compile(str1)
                pattern2 = re.compile(str2)
                ff_data = open(self.force_field_data, 'rt')

                for line in ff_data:
                    if re.search(pattern1, line):
                        n_spaces = (at1 + at2 + at3).count(' ')
                        cur_str = '{:6}{:7}{:7}{:7}'.format(
                            i + 1, j + 1, k + 1, l + 1)
                        cur_str += '{:7}{:11.2f}{:11.5f}{:4}'.format(
                            1, float(line.split()[3 + n_spaces]),
                            float(line.split()[2 + n_spaces]) /
                            float(line.split()[1 + n_spaces]) * 4.184,
                            abs(int(float(line.split()[4 + n_spaces]))))
                        cur_str += ' ;{:>7} -{:>3} -{:>3}-{:>3}\n'.format(
                            atom_names[i], atom_names[j], atom_names[k],
                            atom_names[l])
                        itp_file.write(cur_str)
                        dihedral_available = True
                    elif re.search(pattern2, line):
                        n_spaces = (at2 + at3 + at4).count(' ')
                        cur_str = '{:6}{:7}{:7}{:7}'.format(
                            i + 1, j + 1, k + 1, l + 1)
                        cur_str += '{:7}{:11.2f}{:11.5f}{:4}'.format(
                            1, float(line.split()[3 + n_spaces]),
                            float(line.split()[2 + n_spaces]) /
                            float(line.split()[1 + n_spaces]) * 4.184,
                            abs(int(float(line.split()[4 + n_spaces]))))
                        cur_str += ' ;{:>7} -{:>3} -{:>3}-{:>3}\n'.format(
                            atom_names[i], atom_names[j], atom_names[k],
                            atom_names[l])
                        itp_file.write(cur_str)
                        dihedral_available = True

                if not dihedral_available:
                    str1 = 'X -{}-{}-X'.format(at2, at3)
                    str2 = 'X -{}-{}-X'.format(at3, at2)
                    pattern1 = re.compile(str1)
                    pattern2 = re.compile(str2)
                    n_spaces = 1 + (at2 + at3).count(' ')

                    ff_data = open(self.force_field_data, 'rt')
                    for line in ff_data:
                        if re.search(pattern1, line) or re.search(
                                pattern2, line):
                            cur_str = '{:6}{:7}{:7}{:7}'.format(
                                i + 1, j + 1, k + 1, l + 1)
                            cur_str += '{:7}{:11.2f}{:11.5f}{:4}'.format(
                                1, float(line.split()[3 + n_spaces]),
                                float(line.split()[2 + n_spaces]) /
                                float(line.split()[1 + n_spaces]) * 4.184,
                                int(float(line.split()[4 + n_spaces])))
                            cur_str += ' ;{:>7} -{:>3} -{:>3}-{:>3}\n'.format(
                                atom_names[i], atom_names[j], atom_names[k],
                                atom_names[l])
                            itp_file.write(cur_str)
                            dihedral_available = True

                if not dihedral_available:
                    if Path(self.data_extension).is_file():
                        ff_data = open(self.data_extension, 'rt')
                        for line in ff_data:
                            if re.search(pattern1, line) or re.search(
                                    pattern2, line):
                                cur_str = '{:6}{:7}{:7}{:7}'.format(
                                    i + 1, j + 1, k + 1, l + 1)
                                cur_str += '{:7}{:11.2f}{:11.5f}{:4}'.format(
                                    1, float(line.split()[3 + n_spaces]),
                                    float(line.split()[2 + n_spaces]) /
                                    float(line.split()[1 + n_spaces]) * 4.184,
                                    int(float(line.split()[4 + n_spaces])))
                                cur_str += ' ;{:>7} -{:>3} -{:>3}-{:>3}\n'.format(
                                    atom_names[i], atom_names[j], atom_names[k],
                                    atom_names[l])
                                itp_file.write(cur_str)
                                dihedral_available = True

                if not dihedral_available:
                    cur_str = 'Dihedral {} is not available in data!'.format(
                        str1)
                    cur_str += 'Add parameters to extension file:'
                    print(cur_str)
                    ff_data = open(self.data_extension, 'a')
                    multiplicity = int(input('Multiplicity: '))
                    barrier = float(input('Rotational barrier in kcal/mol: '))
                    phase = float(input('Phase angle: '))
                    periodicity = int(input('Periodicity: '))

                    ff_data.write('{}{:5}{:9.3f}{:14.3f}{:16.3f}\n'.format(
                        str1, multiplicity, barrier, phase, periodicity))
                    cur_str = '{:6}{:7}{:7}{:7}'.format(i + 1, j + 1, k + 1,
                                                        l + 1)
                    cur_str += '{:7}{:11.2f}{:11.5f}{:4}'.format(
                        1, phase, barrier / multiplicity * 4.184,
                        int(periodicity))
                    cur_str += ' ;{:>7} -{:>3} -{:>3}-{:>3}\n'.format(
                        atom_names[i], atom_names[j], atom_names[k],
                        atom_names[l])
                    itp_file.write(cur_str)
                    ff_data.close()

            # improper dihedrals

            itp_file.write('\n[ dihedrals ] ; impropers\n')
            cur_str = ';    i      j      k      l    funct'
            cur_str += '    phase     k_d      n\n'
            itp_file.write(cur_str)

            # impropers for sp2 atom types in GAFF

            sp2_atom_types = [
                'c2', 'n ', 'na', 'nh', 'no', 'p2', 'ce', 'cf', 'cp', 'cq',
                'cu', 'cv', 'px'
            ]
            impropers = []
            for i, at1 in enumerate(atom_types):
                for j, at2 in enumerate(atom_types):
                    if con[i, j]:
                        for k, at3 in enumerate(atom_types[i + 1:],
                                                start=i + 1):
                            if con[j, k]:
                                for l, at4 in enumerate(atom_types[k + 1:],
                                                        start=k + 1):
                                    k_d = 0.
                                    sp2 = True
                                    if con[j, l] and at2 in sp2_atom_types:
                                        k_d = 1.1
                                    elif con[j, l] and at2 == 'c ':
                                        a = [at1, at3, at4]
                                        if a.count('o ') == 1:
                                            k_d = 10.5
                                        else:
                                            k_d = 1.1
                                    elif con[j,
                                             l] and at2 in ['ca', 'cc', 'cd']:
                                        a = [at1, at3, at4]
                                        if a.count('n2') == 2:
                                            k_d = 10.5
                                        else:
                                            k_d = 1.1
                                    elif con[j, l]:
                                        impropers.append([i, j, l, k])
                                        sp2 = False
                                    else:
                                        sp2 = False
                                    if sp2:
                                        cur_str = '{:6}{:7}{:7}{:7}'.format(
                                            i + 1, j + 1, k + 1, l + 1)
                                        cur_str += '{:7}{:11.2f}{:11.5f}{:4}'.format(
                                            1, 180.0, k_d * 4.184, 2)
                                        cur_str += ' ;{:>7} -{:>3} -{:>3}-{:>3}\n'.format(
                                            atom_names[i], atom_names[j],
                                            atom_names[k], atom_names[l])
                                        itp_file.write(cur_str)

            # possible impropers that are left

            if self.add_impropers and len(impropers) > 0:
                print('Add impropers to topology to avoid inversions:')
                print('No.   i   j   k   l')
                for m, imp in enumerate(impropers):
                    [i, j, k, l] = imp
                    print('{:3}{:4}{:4}{:4}{:4}'.format(m + 1, i + 1, j + 1,
                                                        k + 1, l + 1))
                n = int(input('How many impropers do you want to add? '))

                for m in range(n):
                    [i, j, k, l] = impropers[int(input('No.: ')) - 1]
                    barrier = float(input('Rotational barrier in kcal/mol: '))
                    phase = float(input('Phase angle: '))
                    periodicity = int(input('Periodicity: '))
                    cur_str = '{:6}{:7}{:7}{:7}'.format(i + 1, j + 1, k + 1,
                                                        l + 1)
                    cur_str += '{:7}{:11.2f}{:11.5f}{:4}'.format(
                        1, phase, barrier * 4.184, periodicity)
                    cur_str += ' ;{:>7} -{:>3} -{:>3}-{:>3}\n'.format(
                        atom_names[i], atom_names[j], atom_names[k],
                        atom_names[l])
                    itp_file.write(cur_str)

            itp_file.close()

    def dihedral_correction(self, top_file, dihedrals=None):
        """
        Corrects dihedral parameters.

        :param top_file:
            The topology.
        :param dihedrals:
            The dihedrals.
        """

        # Ryckaert-Bellemans function

        def rbpot(phi, c0, c1, c2, c3, c4, c5):
            v = c0 + c1 * np.cos((180 - phi) * 2 * np.pi / 360)
            v += c2 * np.cos((180 - phi) * 2 * np.pi / 360)**2
            v += c3 * np.cos((180 - phi) * 2 * np.pi / 360)**3
            v += c4 * np.cos((180 - phi) * 2 * np.pi / 360)**4
            v += c5 * np.cos((180 - phi) * 2 * np.pi / 360)**5
            return v

        if dihedrals is None:
            dihedrals = self.dihedrals

        itp_file = self.get_included_file(top_file)

        output_dir = Path('{}_dih_corr'.format(self.molecule_name))
        output_dir.mkdir(parents=True, exist_ok=True)

        if itp_file != '{}_new.itp'.format(self.molecule_name):
            shutil.copy(itp_file, '{}_new.itp'.format(self.molecule_name))
            self.write_top('{}_new.top'.format(self.molecule_name),
                           '{}_new.itp'.format(self.molecule_name))

        for i, dih in enumerate(dihedrals):
            dft_scan = self.dft_energies[self.dihedrals.index(dih)]
            angles = self.scan_angles[self.dihedrals.index(dih)]

            output_dir = Path('{}_dih_corr/{}'.format(self.molecule_name,
                                                      i + 1))
            output_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy(
                '{}_new.itp'.format(self.molecule_name),
                '{}_dih_corr/{}/{}_zero.itp'.format(self.molecule_name, i + 1,
                                                    self.molecule_name))
            self.write_top(
                '{}_dih_corr/{}/{}_zero.top'.format(self.molecule_name, i + 1,
                                                    self.molecule_name),
                '{}_zero.itp'.format(self.molecule_name))

            # MM scan with dihedral parameters set to zero

            self.set_dihedral_parameters(
                '{}_dih_corr/{}/{}_zero.itp'.format(self.molecule_name, i + 1,
                                                    self.molecule_name), dih,
                [0.] * 6)
            mm_zero_scan = self.perform_mm_scan(
                '{}_dih_corr/{}/{}_zero.top'.format(self.molecule_name, i + 1,
                                                    self.molecule_name), dih)

            # fitting Ryckaert-Bellemans function to difference between MM zero
            # scan and DFT

            difference = hartree_in_kcalpermol() * 4.184 * (
                np.array(dft_scan) - min(dft_scan)) - (np.array(mm_zero_scan) -
                                                       min(mm_zero_scan))
            initial_coef = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            coef, cv = curve_fit(rbpot, angles, difference, initial_coef)

            # updating parameters

            self.set_dihedral_parameters(
                '{}_new.itp'.format(self.molecule_name), dih, coef.tolist())

    def validate_force_field(self, top_file, dihedrals=None):
        """
        Validates force field by visualization of dihedral potentials.

        :param top_file:
            The topology.
        :param dihedrals:
            The dihedrals.
        """

        if dihedrals is None:
            dihedrals = self.dihedrals

        for i, dih in enumerate(dihedrals):
            mm_scan = self.perform_mm_scan(top_file, dih)
            self.visualize(mm_scan, dih)

    def perform_mm_scan(self, top_file, dihedral):
        """
        Performs MM scan of a specific dihedral.

        :param top_file:
            The topology.
        :param dihedral:
            The dihedral.
        """

        # select scan angles and geometries from DFT data

        angles = self.scan_angles[self.dihedrals.index(dihedral)]
        geometries = self.scan_geometries[self.dihedrals.index(dihedral)]

        itp_file = self.get_included_file(top_file)

        energies = []
        scan_dir = Path(top_file.split('.')[0] + '_scan')
        scan_dir.mkdir(parents=True, exist_ok=True)

        for i, (angle, geom) in enumerate(zip(angles, geometries)):
            output_dir = Path('{}/{}'.format(scan_dir, i + 1))
            output_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy(itp_file, '{}/{}/{}.itp'.format(scan_dir, i + 1, i + 1))
            self.write_top('{}/{}/{}.top'.format(scan_dir, i + 1, i + 1),
                           '{}.itp'.format(i + 1))
            self.write_pdb_file('{}/{}/{}.pdb'.format(scan_dir, i + 1, i + 1),
                                geom)

            # fix dihedral with harmonic restraint function and calculate
            # potential energy

            self.fix_dihedral('{}/{}/{}.itp'.format(scan_dir, i + 1, i + 1),
                              dihedral, angle)
            state = self.minimize_mm_energy(
                '{}/{}/{}.pdb'.format(scan_dir, i + 1, i + 1),
                '{}/{}/{}.top'.format(scan_dir, i + 1, i + 1))
            pot_energy = float(str(state.getPotentialEnergy()).split()[0])

            # remove energy of harmonic restraint function

            if self.remove_restraint_energy:
                coords = np.array(
                    state.getPositions(asNumpy=True).tolist()) * 10
                final_angle = self.get_dihedral_angle(
                    Molecule(self.molecule.get_labels(), coords), dihedral)
                if abs(final_angle - angle) > 180:
                    delta_phi = abs(abs(final_angle - angle) - 360)
                else:
                    delta_phi = abs(final_angle - angle)
                pot_energy -= self.kfac / 2 * (np.deg2rad(delta_phi))**2

            energies.append(pot_energy)

        return energies

    def minimize_mm_energy(self, pdb_file, top_file):
        """
        Minimizes MM energy of a topology using OpenMM.

        :param pdb_file:
            The pdb file containing the initial geometry.
        :param top_file:
            The topology.
        """

        try:
            import openmm
            import openmm.app
            import openmm.unit
        except ImportError:
            raise ImportError(
                'Unable to import OpenMM. Please install ' +
                'OpenMM via \'conda install -c conda forge openmm\'')

        # setup of the system and simulation parameters

        pdb = openmm.app.PDBFile(pdb_file)
        box_length_in_nm = self.box_length / 10
        box_vectors = (openmm.Vec3(box_length_in_nm, 0.,
                                   0.), openmm.Vec3(0., box_length_in_nm, 0.),
                       openmm.Vec3(0., 0., box_length_in_nm))
        top = openmm.app.GromacsTopFile(
            top_file,
            periodicBoxVectors=box_vectors,
            includeDir=self.path_to_required_gromacs_files)
        system = top.createSystem(nonbondedMethod=openmm.app.PME,
                                  nonbondedCutoff=1 * openmm.unit.nanometer,
                                  constraints=openmm.app.HBonds)
        integrator = openmm.LangevinIntegrator(
            self.temperature * openmm.unit.kelvin, 1 / openmm.unit.picosecond,
            0.004 * openmm.unit.picoseconds)

        # minimize energy for given simulation parameters

        simulation = openmm.app.Simulation(top.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy(maxIterations=50000)

        return simulation.context.getState(getPositions=True, getEnergy=True)

    def set_dihedral_parameters(self, itp_filename, dihedral, coefficients):
        """
        Sets dihedral parameters of topology.

        :param itp_filename:
            The internal topology.
        :param dihedral:
            The dihedral.
        :param coefficients:
            The coefficients for Ryckaert-Bellemans funciton.
        """

        atom_names = self.get_atom_names()

        # updating itp file

        with open(itp_filename, 'r+') as itp_file:

            # removing existing dihedral parameters for chosen dihedral

            text = []
            for line in itp_file:
                try:
                    i = int(line.split()[0])
                    j = int(line.split()[1])
                    k = int(line.split()[2])
                    l = int(line.split()[3])
                    if int(line.split()[4]) == 1 or int(line.split()[4]) == 3:
                        updated_dihedral = False
                        if [i, j, k, l] == dihedral or [l, k, j, i] == dihedral:
                            updated_dihedral = True
                        elif {j, k} == {dihedral[1], dihedral[2]}:
                            updated_dihedral = True
                        if not updated_dihedral:
                            text.append(line)
                    else:
                        text.append(line)
                except (ValueError, IndexError):
                    text.append(line)

            # writing new dihedral parameters at the end of the file

            header_present = False
            pattern = re.compile('[ dihedrals ] ; edited')
            itp_file.seek(0)
            for line in text:
                itp_file.write(line)
                if re.search(pattern, line):
                    header_present = True

            if not header_present:
                itp_file.write('\n[ dihedrals ] ; edited\n')
                cur_str = ';    i      j      k      l    type'
                cur_str += '      C0         C1         C2         C3         C4         C5\n'
                itp_file.write(cur_str)

            cur_str = '{:6}{:7}{:7}{:7}{:7}'.format(dihedral[0], dihedral[1],
                                                    dihedral[2], dihedral[3], 3)
            cur_str += '{:11.5f}{:11.5f}{:11.5f}{:11.5f}{:11.5f}{:11.5f}'.format(
                coefficients[0], coefficients[1], coefficients[2],
                coefficients[3], coefficients[4], coefficients[5])
            cur_str += ' ;{:>7} -{:>3} -{:>3}-{:>3}\n'.format(
                atom_names[dihedral[0] - 1], atom_names[dihedral[1] - 1],
                atom_names[dihedral[2] - 1], atom_names[dihedral[3] - 1])
            itp_file.write(cur_str)
            itp_file.truncate()
            itp_file.close()

    def fix_dihedral(self, itp_filename, dihedral, angle):
        """
        Fixing dihedral with harmonic restraint function.

        :param itp_file:
            The internal topology.
        :param dihedral:
            The dihedral.
        :param angle:
            The angle.
        """

        atom_names = self.get_atom_names()

        with open(itp_filename, 'a') as itp_file:
            itp_file.write('\n[ dihedrals ] ; restrained\n')
            itp_file.write(
                ';    i      j      k      l    type      phi      kfac \n')
            cur_str = '{:6}{:7}{:7}{:7}{:7}{:11.2f}{:11.2f}'.format(
                dihedral[0], dihedral[1], dihedral[2], dihedral[3], 2, angle,
                self.kfac)
            cur_str += ' ;{:>7} -{:>3} -{:>3}-{:>3}\n'.format(
                atom_names[dihedral[0] - 1], atom_names[dihedral[1] - 1],
                atom_names[dihedral[2] - 1], atom_names[dihedral[3] - 1])
            itp_file.write(cur_str)
            itp_file.close()

    def write_pdb_file(self, filename, molecule):
        """
        Writes data in PDB file.

        :param filename:
            The name of the file.
        :param molecule:
            The molecule.
        """

        atom_names = self.get_atom_names()
        coords = molecule.get_coordinates() * bohr_in_angstroms()

        with open(filename, 'w') as pdb_file:
            for i in range(molecule.number_of_atoms()):
                cur_str = '{:6s}{:5d} {:^4s} {:3s}  {:4d}    '.format(
                    'ATOM', i + 1, atom_names[i], 'RES', 1)
                cur_str += '{:8.3f}{:8.3f}{:8.3f}{:6.2f}\n'.format(
                    coords[i][0], coords[i][1], coords[i][2], 1.0)
                pdb_file.write(cur_str)
            pdb_file.close()

    def visualize(self, mm_scan, dihedral):
        """
        Visualizes dihedral potential.

        :param mm_scan:
            The MM scan energies.
        :param dihedral:
            The dihedral.
        """

        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError('Unable to import Matplotlib. Please install ' +
                              'Matplotlib via \'conda install matplotlib\'')

        angles = self.scan_angles[self.dihedrals.index(dihedral)]
        dft_scan = self.dft_energies[self.dihedrals.index(dihedral)]
        dft_scan = hartree_in_kcalpermol() * (np.array(dft_scan) -
                                              min(dft_scan))
        mm_scan = (np.array(mm_scan) - min(mm_scan)) / 4.184
        [i, j, k, l] = dihedral

        plt.plot(angles, dft_scan, '-o', label="DFT")
        plt.plot(angles, mm_scan, '-o', label="MM")

        plt.grid()
        plt.legend(loc='upper right')
        plt.xlabel('dihedral angle {}-{}-{}-{}'.format(i, j, k, l))
        plt.ylabel('E in kcal/mol')
        plt.title('dihedral potential of {}'.format(self.molecule_name))
        plt.show()

    def get_dihedral_angle(self, molecule, dihedral):
        """
        Gets dihedral angle of a molecule in range of -180 and 180 degrees.

        :param molecule:
            The molecule.
        :param dihedral:
            The dihedral.
        """

        # this function is from https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python

        coords = molecule.get_coordinates()
        p = []
        for i in dihedral:
            p.append(coords[i - 1])
        p = np.array(p)
        b = p[:-1] - p[1:]
        b[0] *= -1
        v = np.array(
            [v - (v.dot(b[1]) / b[1].dot(b[1])) * b[1] for v in [b[0], b[2]]])
        v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1, 1)
        b1 = b[1] / np.linalg.norm(b[1])
        x = np.dot(v[0], v[1])
        y = np.dot(np.cross(v[0], b1), v[1])

        return np.degrees(np.arctan2(y, x))

    def get_connectivity(self):
        """
        Gets connectivity.
        """

        coords = self.molecule.get_coordinates()
        n_atoms = self.molecule.number_of_atoms()
        covalent_radii = self.molecule.covalent_radii_to_numpy()

        connectivity = np.full((n_atoms, n_atoms), False, dtype=bool)
        tolerance = 0.4 / bohr_in_angstroms()
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                r_ij = np.linalg.norm(coords[i] - coords[j])
                if abs(r_ij - covalent_radii[i] -
                       covalent_radii[j]) < tolerance:
                    connectivity[i, j] = True
                    connectivity[j, i] = True

        return connectivity

    def get_included_file(self, top_file):
        """
        Gets the name of the included itp file.

        :param top_file:
            The topology file.
        """

        with open(top_file, 'rt') as top:
            pattern = re.compile('include')
            for line in top:
                if re.search(pattern, line):
                    return '{}/{}'.format(str(Path(top_file).parent),
                                          line.split('\"')[1])

    def get_atom_names(self):
        """
        Gets the atom names.
        """

        atom_names = []
        counter = {}
        for label in self.molecule.get_labels():
            if label not in counter:
                counter[label] = 1
            else:
                counter[label] += 1
            atom_names.append(label + str(counter[label]))

        for i, label in enumerate(self.molecule.get_labels()):
            if counter[label] == 1:
                atom_names[i] = label

        return atom_names
