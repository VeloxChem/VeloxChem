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

from pathlib import Path, PurePath
from scipy.optimize import curve_fit
import numpy as np
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
        - force_field_data_extension: The filename of the data file with user-specified
          parameters.
        - box_length: The box length for energy minimizations with OpenMM (Angstroms).
        - temperature: The temperature (Kelvin).
        - gromacs_include_path: The path to share/gromacs/top of the GROMACS directory.
        - kfac: The dihedral force constant to restrain dihedrals.
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
        self.force_field_data = 'gaff2.dat'
        self.force_field_data_extension = None

        # mm settings
        self.box_length = 100.0
        self.temperature = 293.15
        self.gromacs_include_path = None
        self.kfac = 10000.0

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
                with open(xyz_file, 'r') as xyz:
                    e = []
                    phi = []
                    pattern = re.compile(r'\AScan')
                    for line in xyz:
                        if re.search(pattern, line):
                            e.append(float(line.split('Energy')[1].split()[0]))
                            phi.append(float(line.split('=')[1].split()[0]))
                            dih = [
                                int(i) for i in line.split('Dihedral')
                                [1].split()[0].split('-')
                            ]
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
            'gromacs_include_path': 'str',
            'kfac': 'float',
            'force_field_data': 'str',
            'force_field_data_extension': 'str',
        }

        parse_input(self, ffg_keywords, ffg_dict)

        if self.force_field_data_extension is None:
            ff_file = PurePath(self.force_field_data)
            self.force_field_data_extension = str(
                ff_file.parent / (ff_file.stem + '_extension.dat'))

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
        angle_indices = sorted(list(angle_indices))

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
        dihedral_indices = sorted(list(dihedral_indices))

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

        if Path(self.force_field_data_extension).is_file():
            with open(self.force_field_data, 'r') as ff_extension:
                ff_data_lines += ff_extension.readlines()

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
                    msg = f'Warning: bond length {at_1}-{at_2}'
                    msg += ' does not match length in XYZ file'
                    msg += f' (atoms no. {i+1} and {j+1}):'
                    msg += f' {r:.3f} vs. {r_eq:.3f}'
                    print(msg)

                if self.eq_param:
                    r = r_eq

                cur_str = '{:6}{:7}{:7}{:14.4e}{:14.4e} ; {}-{}\n'.format(
                    i + 1, j + 1, 1, r, k_r, atom_names[i], atom_names[j])
                itp_file.write(cur_str)

            # pairs

            if self.gen_pairs:
                itp_file.write('\n[ pairs ]\n')
                itp_file.write(';    i      j    funct\n')
                for i, j in pairs_14:
                    itp_file.write('{:6}{:7}{:7}\n'.format(i + 1, j + 1, 1))

            # angles

            itp_file.write('\n[ angles ]\n')
            cur_str = ';    i      j      k    funct'
            cur_str += '     theta       k_theta\n'
            itp_file.write(cur_str)

            for i, j, k in angle_indices:
                a = coords[i] - coords[j]
                b = coords[k] - coords[j]
                theta_eq = np.arccos(
                    np.dot(a, b) / np.linalg.norm(a) /
                    np.linalg.norm(b)) * 180 / np.pi

                at_1 = atom_types[i]
                at_2 = atom_types[j]
                at_3 = atom_types[k]
                str_1 = r'\A' + f'{at_1}-{at_2}-{at_3} '
                str_2 = r'\A' + f'{at_3}-{at_2}-{at_1} '
                pattern_1 = re.compile(str_1)
                pattern_2 = re.compile(str_2)

                angle_found = False

                for line in ff_data_lines:
                    if (re.search(pattern_1, line) or
                            re.search(pattern_2, line)):
                        angle_ff = line[8:].strip().split()
                        theta = float(angle_ff[1])
                        k_theta = float(angle_ff[0]) * 4.184 * 2
                        angle_found = True

                errmsg = f'ForceFieldGenerator: angle {at_1}-{at_2}-{at_3}'
                errmsg += ' is not available.'
                assert_msg_critical(angle_found, errmsg)

                if abs(theta - theta_eq) > self.theta_thresh:
                    msg = f'Warning: angle {at_1}-{at_2}-{at_3}'
                    msg += ' does not match angle in XYZ file'
                    msg += f' (atoms no. {i+1}, {j+1} and {k+1}):'
                    msg += f' {theta:.1f} vs. {theta_eq:.1f}!'
                    print(msg)

                if self.eq_param:
                    theta = theta_eq

                cur_str = '{:6}{:7}{:7}{:7}{:14.4e}{:14.4e}'.format(
                    i + 1, j + 1, k + 1, 1, theta, k_theta)
                cur_str += ' ; {}-{}-{}\n'.format(atom_names[i], atom_names[j],
                                                  atom_names[k])
                itp_file.write(cur_str)

            # proper dihedrals

            itp_file.write('\n[ dihedrals ] ; propers\n')
            cur_str = ';    i      j      k      l    funct'
            cur_str += '    phase     k_d      n\n'
            itp_file.write(cur_str)
            cur_str = ';                                        '
            cur_str += 'C0         C1         C2         C3         C4         C5\n'
            itp_file.write(cur_str)

            for i, j, k, l in dihedral_indices:
                at_1 = atom_types[i]
                at_2 = atom_types[j]
                at_3 = atom_types[k]
                at_4 = atom_types[l]
                str_1 = '{}-{}-{}-{}'.format(at_1, at_2, at_3, at_4)
                str_2 = '{}-{}-{}-{}'.format(at_4, at_3, at_2, at_1)
                pattern_1 = re.compile(str_1)
                pattern_2 = re.compile(str_2)

                dihedral_found = False

                dihedral_ff_lines = []
                for line in ff_data_lines:
                    if (re.search(pattern_1, line) or
                            re.search(pattern_2, line)):
                        if '.' not in line[11:].strip().split()[0]:
                            dihedral_ff_lines.append(line)
                            dihedral_found = True

                if not dihedral_found:
                    str_1 = 'X -{}-{}-X '.format(at_2, at_3)
                    str_2 = 'X -{}-{}-X '.format(at_3, at_2)
                    pattern_1 = re.compile(str_1)
                    pattern_2 = re.compile(str_2)

                    dihedral_ff_lines = []
                    for line in ff_data_lines:
                        if (re.search(pattern_1, line) or
                                re.search(pattern_2, line)):
                            if '.' not in line[11:].strip().split()[0]:
                                dihedral_ff_lines.append(line)
                                dihedral_found = True

                errmsg = 'ForceFieldGenerator: dihedral'
                errmsg += f' {at_1}-{at_2}-{at_3}-{at_4} is not available.'
                assert_msg_critical(dihedral_found, errmsg)

                for line in dihedral_ff_lines:
                    dihedral_ff = line[11:].strip().split()
                    multiplicity = int(dihedral_ff[0])
                    barrier = float(dihedral_ff[1]) * 4.184 / multiplicity
                    phase = float(dihedral_ff[2])
                    try:
                        periodicity = abs(int(dihedral_ff[3]))
                    except ValueError:
                        periodicity = abs(int(float(dihedral_ff[3])))

                    cur_str = '{:6}{:7}{:7}{:7}'.format(i + 1, j + 1, k + 1,
                                                        l + 1)
                    cur_str += '{:7}{:11.2f}{:11.5f}{:4}'.format(
                        1, phase, barrier, periodicity)
                    cur_str += ' ; {}-{}-{}-{}\n'.format(
                        atom_names[i], atom_names[j], atom_names[k],
                        atom_names[l])
                    itp_file.write(cur_str)

            # improper dihedrals

            itp_file.write('\n[ dihedrals ] ; impropers\n')
            cur_str = ';    i      j      k      l    funct'
            cur_str += '    phase     k_d      n\n'
            itp_file.write(cur_str)

            for i, j, k in angle_indices:
                at_1 = atom_types[i]
                at_2 = atom_types[j]
                at_3 = atom_types[k]

                for l in range(n_atoms):
                    if (l in [i, j, k]) or (not con[l, j]):
                        continue
                    at_4 = atom_types[l]

                    str_1 = '{}-{}-{}-{}'.format(at_4, at_1, at_2, at_3)
                    str_2 = '{}-{}-{}-{}'.format(at_3, at_4, at_2, at_1)
                    str_3 = '{}-{}-{}-{}'.format(at_1, at_3, at_2, at_4)
                    pattern_1 = re.compile(str_1)
                    pattern_2 = re.compile(str_2)
                    pattern_3 = re.compile(str_3)

                    dihedral_found = False

                    for line in ff_data_lines:
                        if (re.search(pattern_1, line) or
                                re.search(pattern_2, line) or
                                re.search(pattern_3, line)):
                            if '.' in line[11:].strip().split()[0]:
                                dihedral_ff = line[11:].strip().split()
                                dihedral_found = True
                                break

                    if not dihedral_found:
                        str_1 = 'X -{}-{}-{}'.format(at_1, at_2, at_3)
                        str_2 = 'X -{}-{}-{}'.format(at_4, at_2, at_1)
                        str_3 = 'X -{}-{}-{}'.format(at_3, at_2, at_4)
                        pattern_1 = re.compile(str_1)
                        pattern_2 = re.compile(str_2)
                        pattern_3 = re.compile(str_3)

                        for line in ff_data_lines:
                            if (re.search(pattern_1, line) or
                                    re.search(pattern_2, line) or
                                    re.search(pattern_3, line)):
                                if '.' in line[11:].strip().split()[0]:
                                    dihedral_ff = line[11:].strip().split()
                                    dihedral_found = True
                                    break

                    if not dihedral_found:
                        str_1 = 'X -X -{}-{}'.format(at_2, at_3)
                        str_2 = 'X -X -{}-{}'.format(at_2, at_1)
                        str_3 = 'X -X -{}-{}'.format(at_2, at_4)
                        pattern_1 = re.compile(str_1)
                        pattern_2 = re.compile(str_2)
                        pattern_3 = re.compile(str_3)

                        for line in ff_data_lines:
                            if (re.search(pattern_1, line) or
                                    re.search(pattern_2, line) or
                                    re.search(pattern_3, line)):
                                if '.' in line[11:].strip().split()[0]:
                                    dihedral_ff = line[11:].strip().split()
                                    dihedral_found = True
                                    break

                    if dihedral_found:
                        barrier = float(dihedral_ff[0]) * 4.184
                        phase = float(dihedral_ff[2])
                        periodicity = abs(int(dihedral_ff[3]))

                        assert_msg_critical(
                            phase == 180.0,
                            'ForceFieldGenerator: invalid improper dihedral phase'
                        )
                        assert_msg_critical(
                            periodicity == 2,
                            'ForceFieldGenerator: invalid improper dihedral periodicity'
                        )

                        cur_str = '{:6}{:7}{:7}{:7}'.format(
                            l + 1, i + 1, j + 1, k + 1)
                        cur_str += '{:7}{:11.2f}{:11.5f}{:4}'.format(
                            1, phase, barrier, periodicity)
                        cur_str += ' ; {}-{}-{}-{}\n'.format(
                            atom_names[l], atom_names[i], atom_names[j],
                            atom_names[k])
                        itp_file.write(cur_str)

    @staticmethod
    def copy_file(src, dest):
        """
        Copies file (from src to dest).

        :param src:
            The source of copy.
        :param dest:
            The destination of copy.
        """

        dest.write_text(src.read_text())

    def dihedral_correction(self, top_filename):
        """
        Corrects dihedral parameters.

        :param top_filename:
            The name of topology file.
        """

        # Ryckaert-Bellemans function

        def rbpot(phi, c0, c1, c2, c3, c4, c5):
            v = c0 + c1 * np.cos((180 - phi) * 2 * np.pi / 360)
            v += c2 * np.cos((180 - phi) * 2 * np.pi / 360)**2
            v += c3 * np.cos((180 - phi) * 2 * np.pi / 360)**3
            v += c4 * np.cos((180 - phi) * 2 * np.pi / 360)**4
            v += c5 * np.cos((180 - phi) * 2 * np.pi / 360)**5
            return v

        itp_fname = self.get_included_file(top_filename)

        new_itp_fname = f'{self.molecule_name}_new.itp'
        new_top_fname = f'{self.molecule_name}_new.top'

        if itp_fname != new_itp_fname:
            self.copy_file(Path(itp_fname), Path(new_itp_fname))
            self.write_top(new_top_fname, new_itp_fname)

        for i, (dih, geom, dft_scan, angles) in enumerate(
                zip(self.dihedrals, self.scan_geometries, self.dft_energies,
                    self.scan_angles)):

            output_dir = Path(f'{self.molecule_name}_dih_corr', f'{i+1}')
            output_dir.mkdir(parents=True, exist_ok=True)

            zero_itp_file = output_dir / f'{self.molecule_name}_zero.itp'
            zero_top_file = output_dir / f'{self.molecule_name}_zero.top'

            self.copy_file(Path(new_itp_fname), zero_itp_file)
            self.write_top(str(zero_top_file), zero_itp_file.name)

            # MM scan with dihedral parameters set to zero

            self.set_dihedral_parameters(str(zero_itp_file), dih, [0.] * 6)
            mm_zero_scan = self.perform_mm_scan(str(zero_top_file), dih, geom,
                                                angles)

            # fitting Ryckaert-Bellemans function to difference between MM zero
            # scan and DFT

            rel_e_dft = np.array(dft_scan) - min(dft_scan)
            rel_e_dft *= hartree_in_kcalpermol() * 4.184
            rel_e_mm = np.array(mm_zero_scan) - min(mm_zero_scan)

            difference = rel_e_dft - rel_e_mm
            initial_coef = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            coef, cv = curve_fit(rbpot, angles, difference, initial_coef)

            # updating parameters

            self.set_dihedral_parameters(
                '{}_new.itp'.format(self.molecule_name), dih, coef.tolist())

    def validate_force_field(self, top_file):
        """
        Validates force field by RMSD of dihedral potentials.

        :param top_file:
            The topology.
        :param dihedrals:
            The dihedrals.
        """

        for i, dih in enumerate(self.dihedrals):
            geom = self.scan_geometries[i]
            angles = self.scan_angles[i]
            mm_scan = self.perform_mm_scan(top_file, dih, geom, angles)

            dft_scan = np.array(self.dft_energies[i]) - min(
                self.dft_energies[i])
            dft_scan *= hartree_in_kcalpermol() * 4.184
            dft_scan += min(mm_scan)

            rmsd = 0.0
            for e_dft, e_mm in zip(dft_scan, mm_scan):
                rmsd += (e_dft - e_mm)**2
            rmsd = np.sqrt(rmsd / len(mm_scan))

            print(f'RMSD: {rmsd:.3f} kJ/mol')

    def perform_mm_scan(self, top_filename, dihedral, geometries, angles):
        """
        Performs MM scan of a specific dihedral.

        :param top_filename:
            The name of topology file.
        :param dihedral:
            The dihedral.
        """

        # select scan angles and geometries from DFT data

        itp_fname = self.get_included_file(top_filename)

        top_file = Path(top_filename)

        scan_dir = top_file.parent / (top_file.stem + '_scan')
        scan_dir.mkdir(parents=True, exist_ok=True)

        print('== MM scan ==')

        energies = []
        for i, (geom, angle) in enumerate(zip(geometries, angles)):

            output_dir = scan_dir / f'{i+1}'
            output_dir.mkdir(parents=True, exist_ok=True)

            local_itp_fname = str(output_dir / f'{i+1}.itp')
            local_top_fname = str(output_dir / f'{i+1}.top')
            local_pdb_fname = str(output_dir / f'{i+1}.pdb')

            self.copy_file(Path(itp_fname), Path(local_itp_fname))
            self.write_top(local_top_fname, Path(local_itp_fname).name)
            self.write_pdb_file(local_pdb_fname, geom)

            # fix dihedral with harmonic restraint function and calculate
            # potential energy
            self.fix_dihedral(local_itp_fname, dihedral, angle)
            state = self.minimize_mm_energy(local_pdb_fname, local_top_fname)

            pot_energy_str = str(state.getPotentialEnergy())
            pot_energy_unit = pot_energy_str.split()[1]
            assert_msg_critical(
                pot_energy_unit == 'kJ/mol',
                'ForceFieldGenerator.perform_mm_scan: ' +
                'unexpected unit for potential energy')
            pot_energy = float(pot_energy_str.split()[0])
            energies.append(pot_energy)

            print(f'  {local_top_fname}: {pot_energy:.3f} kJ/mol')

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
        except ImportError:
            raise ImportError(
                'Unable to import OpenMM. Please install ' +
                'OpenMM via \'conda install -c conda-forge openmm\'')

        # setup of the system and simulation parameters

        pdb = openmm.app.PDBFile(pdb_file)
        box_length_in_nm = self.box_length / 10
        box_vectors = (
            openmm.Vec3(box_length_in_nm, 0., 0.),
            openmm.Vec3(0., box_length_in_nm, 0.),
            openmm.Vec3(0., 0., box_length_in_nm),
        )

        top = openmm.app.GromacsTopFile(top_file,
                                        periodicBoxVectors=box_vectors,
                                        includeDir=self.gromacs_include_path)
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

        # read itp file and remove existing parameters for chosen dihedral

        dihedral_flag = False
        dihedral_pattern = re.compile(r'\[\s*dihedrals\s*\]')
        saved_itp_lines = []

        with open(itp_filename, 'r') as f_itp:
            for line in f_itp:
                title = line.split(';')[0].strip()

                if re.search(dihedral_pattern, title):
                    dihedral_flag = True
                elif title.startswith('['):
                    dihedral_flag = False

                if dihedral_flag:
                    content = line.split()
                    try:
                        i, j, k, l, funct = tuple(
                            [int(content[n]) for n in range(5)])
                        condition_1 = (funct in [1, 3])
                        condition_2 = ([i, j, k, l] == dihedral or
                                       [l, k, j, i] == dihedral or
                                       [j, k] == dihedral[1:3] or
                                       [k, j] == dihedral[1:3])
                        if not (condition_1 and condition_2):
                            saved_itp_lines.append(line)
                    except (ValueError, IndexError):
                        saved_itp_lines.append(line)
                else:
                    saved_itp_lines.append(line)

        # update itp file with constraints for chosen dihedral

        with open(itp_filename, 'w') as f_itp:
            for line in saved_itp_lines:
                f_itp.write(line)

            f_itp.write('\n[ dihedrals ] ; edited\n')
            cur_str = ';    i      j      k      l    type'
            cur_str += '      C0         C1         C2         C3'
            cur_str += '         C4         C5\n'
            f_itp.write(cur_str)

            cur_str = '{:6}{:7}{:7}{:7}{:7}'.format(*dihedral, 3)
            for coef in coefficients:
                cur_str += '{:11.5f}'.format(coef)
            cur_str += ' ; {}-{}-{}-{}\n'.format(atom_names[dihedral[0] - 1],
                                                 atom_names[dihedral[1] - 1],
                                                 atom_names[dihedral[2] - 1],
                                                 atom_names[dihedral[3] - 1])
            f_itp.write(cur_str)

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

        with open(itp_filename, 'a') as f_itp:
            f_itp.write('\n[ dihedrals ] ; restrained\n')
            f_itp.write(
                ';    i      j      k      l    type      phi      kfac \n')
            cur_str = '{:6}{:7}{:7}{:7}{:7}{:11.2f}{:11.2f}'.format(
                *dihedral, 2, angle, self.kfac)
            cur_str += ' ; {}-{}-{}-{}\n'.format(atom_names[dihedral[0] - 1],
                                                 atom_names[dihedral[1] - 1],
                                                 atom_names[dihedral[2] - 1],
                                                 atom_names[dihedral[3] - 1])
            f_itp.write(cur_str)

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

    def visualize(self, dft_scan, mm_scan, dihedral, angles):
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

        dft_scan_kcal = (np.array(dft_scan) -
                         min(dft_scan)) * hartree_in_kcalpermol()
        mm_scan_kcal = (np.array(mm_scan) - min(mm_scan)) / 4.184

        plt.plot(angles, dft_scan_kcal, '-o', label="DFT")
        plt.plot(angles, mm_scan_kcal, '-o', label="MM")

        plt.grid()
        plt.legend(loc='upper right')
        plt.xlabel('dihedral angle {}-{}-{}-{}'.format(*dihedral))
        plt.ylabel('E in kcal/mol')
        plt.title('dihedral potential of {}'.format(self.molecule_name))
        plt.show()

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

        with open(top_file, 'r') as top:
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
