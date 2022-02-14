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

from mpi4py import MPI
from pathlib import Path
from scipy.optimize import curve_fit
import numpy as np
import sys
import re

from .veloxchemlib import mpi_master, bohr_in_angstroms, hartree_in_kcalpermol
from .molecule import Molecule
from .outputstream import OutputStream
from .respchargesdriver import RespChargesDriver
from .inputparser import parse_input, get_datetime_string
from .errorhandler import assert_msg_critical


class ForceFieldGenerator:
    """
    Parameterizes general Amber force field and creates Gromacs topologies.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - molecule_name: The name of the molecule.
        - eq_param: If equilibrium bond lengths and angles should be used.
        - r_thresh: The threshold for warning if bond lenghts deviate (nm).
        - theta_thresh: The threshold for warning if bond angle deviate (deg).
        - gen_pairs: If 1-4 parameters should be generated.
        - fudgeLJ: The factor by which to multiply Lennard-Jones 1-4 interactions.
        - fudgeQQ: The factor by which to multiply electrostatic 1-4 interactions.
        - comb_rule: The number of the combination rule in Gromacs.
        - nbfunc: The non-bonded function type (1: Lennard-Jones, 2: Buckingham).
        - nrexcl: Number of neighboring bonds to be excluded for non-bonded interaction.
        - force_field_data: The filename of the data file with force field parameters.
        - force_field_data_extension: The filename of the data file with user-specified
          parameters.
        - box_length: The box length for energy minimization.
        - temperature: The temperature.
        - gromacs_include_path: The path to topology folder of Gromacs.
        - kfac: The force constant to restrain dihedrals.
        - scan_dih_angles: The dihedral angles for MM scans (list of list).
        - scan_energies: The energies from QM scans (list of list).
        - scan_geometries: The optimized geometries from QM scan (list of list).
        - target_dihedrals: The target dihedral angles for parameterization.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes force field generator.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        self.molecule_name = f'veloxchem_ff_{get_datetime_string()}'
        self.scan_xyz_files = None
        self.atom_types = None

        # topology settings
        self.eq_param = True
        self.r_thresh = 0.005
        self.theta_thresh = 10.
        self.gen_pairs = True
        self.fudgeLJ = 0.5
        self.fudgeQQ = (1.0 / 1.2)
        self.comb_rule = 2
        self.nbfunc = 1
        self.nrexcl = 3
        self.force_field_data = 'gaff2.dat'
        self.force_field_data_extension = None

        # MM settings
        self.box_length = 100.0  # Angstrom
        self.temperature = 293.15  # Kelvin
        self.gromacs_include_path = None
        self.kfac = 10000.0  # force constant in kJ/mol/rad^2

        # scan settings
        self.scan_dih_angles = None
        self.scan_energies = None
        self.scan_geometries = None
        self.target_dihedrals = None

        # resp settings
        self.resp_dict = None

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

    def update_settings(self, ffg_dict, resp_dict=None):
        """
        Updates settings in force field generator.

        :param ffg_dict:
            The input dictionary of force field group.
        :param resp_dict:
            The input dictionary of resp charges group.
        """

        ffg_keywords = {
            'molecule_name': 'str',
            'scan_xyz_files': 'seq_fixed_str',
            'atom_types': 'seq_fixed_str',
            'eq_param': 'bool',
            'r_thresh': 'float',
            'theta_thresh': 'float',
            'gen_pairs': 'bool',
            'fudgeLJ': 'float',
            'fudgeQQ': 'float',
            'comb_rule': 'int',
            'nbfunc': 'int',
            'nrexcl': 'int',
            'box_length': 'float',
            'temperature': 'float',
            'gromacs_include_path': 'str',
            'kfac': 'float',
            'force_field_data': 'str',
            'force_field_data_extension': 'str',
        }

        parse_input(self, ffg_keywords, ffg_dict)

        if self.force_field_data_extension is None:
            ff_file = Path(self.force_field_data)
            self.force_field_data_extension = str(
                ff_file.parent / (ff_file.stem + '_extension.dat'))

        if 'filename' in ffg_dict and 'molecule_name' not in ffg_dict:
            self.molecule_name = ffg_dict['filename']

        if resp_dict is None:
            resp_dict = {}
        self.resp_dict = dict(resp_dict)

    def compute(self, molecule, basis):
        """
        Runs force field optimization.
        Note: Only runs on the master MPI process.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        """

        # sanity check

        assert_msg_critical(
            self.scan_xyz_files is not None,
            'ForceFieldGenerator.compute: scan_xyz_files not defined ')

        assert_msg_critical(
            self.atom_types is not None,
            'ForceFieldGenerator.compute: atom_types not defined ')

        assert_msg_critical(
            len(self.atom_types) == molecule.number_of_atoms(),
            'ForceFieldGenerator.compute: inconsistent number of atom_types')

        self.molecule = molecule

        # RESP charges

        resp_drv = RespChargesDriver(self.comm, self.ostream)
        resp_drv.update_settings(self.resp_dict)
        resp_chg = resp_drv.compute(self.molecule, basis, 'resp')

        # read QM scan

        title = 'Force Field Generator'
        self.ostream.print_header(title)
        self.ostream.print_header('=' * (len(title) + 2))
        self.ostream.print_blank()

        self.ostream.print_info('Reading QM scan from file...')
        for xyz_fname in self.scan_xyz_files:
            self.ostream.print_info(f'  {xyz_fname}')
        self.ostream.print_blank()
        self.ostream.flush()

        self.scan_dih_angles = []
        self.scan_energies = []
        self.scan_geometries = []
        self.target_dihedrals = []

        if self.rank == mpi_master():
            self.read_qm_scan_xyz_files(self.scan_xyz_files)

            inp_dir = Path(self.molecule_name).parent
            mol_name = Path(self.molecule_name).stem

            original_itp_file = inp_dir / (mol_name + '_original.itp')
            original_top_file = original_itp_file.with_suffix('.top')

            self.write_original_itp(str(original_itp_file),
                                    list(self.atom_types), resp_chg)
            self.write_top(str(original_top_file), str(original_itp_file))

            self.validate_force_field(str(original_top_file))

            self.dihedral_correction(str(original_top_file))

            new_top_file = inp_dir / (mol_name + '_new.top')
            self.validate_force_field(str(new_top_file))

    def read_qm_scan_xyz_files(self, scan_xyz_files):
        """
        Reads QM scan xyz files.

        :param scan_xyz_files:
            The list of xyz files from QM scan.
        """

        # reading QM data

        try:
            import MDAnalysis as mda
        except ImportError:
            raise ImportError('Unable to import MDAnalysis. Please install ' +
                              'MDAnalysis via \'conda install MDAnalysis\'')

        for xyz_fname in scan_xyz_files:
            geometries = []
            energies = []
            dih_angles = []

            u = mda.Universe(xyz_fname)
            for ts in u.trajectory:
                geometries.append(
                    Molecule(u.atoms.names, u.atoms.positions, 'angstrom'))
            self.scan_geometries.append(geometries)

            pattern = re.compile(r'\AScan')
            with open(xyz_fname, 'r') as f_xyz:
                for line in f_xyz:
                    if re.search(pattern, line):
                        energies.append(
                            float(line.split('Energy')[1].split()[0]))
                        dih_angles.append(float(line.split('=')[1].split()[0]))
                        dih_inds = [
                            int(i) for i in line.split('Dihedral')[1].split()
                            [0].split('-')
                        ]
                self.scan_energies.append(energies)
                self.scan_dih_angles.append(dih_angles)
                self.target_dihedrals.append(dih_inds)

    def write_top(self, top_filename, itp_filename):
        """
        Writes a topology file.

        :param top_filename:
            The name of the topology file.
        :param itp_filename:
            The name of the included itp file.
        """

        mol_name = Path(self.molecule_name).stem
        itp_name = Path(itp_filename).name

        with open(top_filename, 'w') as top_file:

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

            top_file.write('\n#include "' + itp_name + '"\n')

            # system

            top_file.write('\n[ system ]\n')
            top_file.write(' {}\n'.format(mol_name))

            # molecules

            top_file.write('\n[ molecules ]\n')
            top_file.write('; Compound        nmols\n')
            top_file.write('{:>10}{:9}\n'.format(mol_name, 1))

    def write_original_itp(self, itp_fname, atom_types, charges):
        """
        Writes an itp file with the original parameters.

        :param itp_fname:
            The name of itp file.
        :param atom_types:
            The atom types.
        :param charges:
            The charges.
        """

        mol_name = Path(self.molecule_name).stem

        # molecular information

        coords = self.molecule.get_coordinates()
        n_atoms = self.molecule.number_of_atoms()
        connected = self.get_connectivity()

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
                if connected[i, j]:
                    bond_indices.add((i, j))
        bond_indices = sorted(list(bond_indices))

        angle_indices = set()
        for i, j in bond_indices:
            for k in range(n_atoms):
                if k in [i, j]:
                    continue
                if connected[j, k]:
                    inds = (i, j, k) if i < k else (k, j, i)
                    angle_indices.add(inds)
                if connected[k, i]:
                    inds = (k, i, j) if k < j else (j, i, k)
                    angle_indices.add(inds)
        angle_indices = sorted(list(angle_indices))

        dihedral_indices = set()
        for i, j, k in angle_indices:
            for l in range(n_atoms):
                if l in [i, j, k]:
                    continue
                if connected[k, l]:
                    inds = (i, j, k, l) if i < l else (l, k, j, i)
                    dihedral_indices.add(inds)
                if connected[l, i]:
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

        with open(itp_fname, 'w') as f_itp:

            # header

            f_itp.write('; Generated by VeloxChem\n')

            # atom types

            f_itp.write('\n[ atomtypes ]\n')
            cur_str = ';name   bond_type     mass     charge'
            cur_str += '   ptype   sigma         epsilon\n'
            f_itp.write(cur_str)

            for at in unique_atom_types:
                atom_type_found = False

                for line in ff_data_lines:
                    if line.startswith(f'  {at}  '):
                        cur_str = '{:>3}{:>9}{:17.5f}{:9.5f}{:>4}'.format(
                            at, at, 0., 0., 'A')
                        cur_str += '{:16.5e}{:14.5e}\n'.format(
                            float(line.split()[1]) * 2**(-1 / 6) * 2 / 10,
                            float(line.split()[2]) * 4.184)
                        f_itp.write(cur_str)
                        atom_type_found = True
                        break

                assert_msg_critical(
                    atom_type_found,
                    f'ForceFieldGenerator: Unknown atom type {at}')

            # molecule type

            f_itp.write('\n[ moleculetype ]\n')
            f_itp.write(';name            nrexcl\n')
            f_itp.write('{:>10}{:10}\n'.format(mol_name, self.nrexcl))

            # atoms

            f_itp.write('\n[ atoms ]\n')
            cur_str = ';   nr  type  resi  res  atom  cgnr'
            cur_str += '     charge      mass\n'
            f_itp.write(cur_str)

            for i, at in enumerate(atom_types):
                for line in ff_data_lines:
                    if line.startswith(f'{at} '):
                        mass_i = float(line.split()[1])
                        cur_str = '{:6}{:>5}{:6}{:>6}{:>6}'.format(
                            i + 1, atom_types[i], 1, 'RES', atom_names[i])
                        cur_str += '{:5}{:13.6f}{:13.5f}\n'.format(
                            i + 1, charges[i], mass_i)
                        f_itp.write(cur_str)

            # bonds

            f_itp.write('\n[ bonds ]\n')
            f_itp.write(';    i      j    funct       r           k_r\n')

            for i, j in bond_indices:
                r_eq = np.linalg.norm(coords[i] - coords[j])
                r_eq *= bohr_in_angstroms() * 0.1

                at_1 = atom_types[i]
                at_2 = atom_types[j]
                patterns = [
                    re.compile(r'\A' + f'{at_1}-{at_2}  '),
                    re.compile(r'\A' + f'{at_2}-{at_1}  '),
                ]

                bond_found = False

                for line in ff_data_lines:
                    matches = [re.search(p, line) for p in patterns]
                    if any(matches):
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
                f_itp.write(cur_str)

            # pairs

            f_itp.write('\n[ pairs ]\n')
            f_itp.write(';    i      j    funct\n')
            for i, j in pairs_14:
                f_itp.write('{:6}{:7}{:7}\n'.format(i + 1, j + 1, 1))

            # angles

            f_itp.write('\n[ angles ]\n')
            cur_str = ';    i      j      k    funct'
            cur_str += '     theta       k_theta\n'
            f_itp.write(cur_str)

            for i, j, k in angle_indices:
                a = coords[i] - coords[j]
                b = coords[k] - coords[j]
                theta_eq = np.arccos(
                    np.dot(a, b) / np.linalg.norm(a) /
                    np.linalg.norm(b)) * 180 / np.pi

                at_1 = atom_types[i]
                at_2 = atom_types[j]
                at_3 = atom_types[k]
                patterns = [
                    re.compile(r'\A' + f'{at_1}-{at_2}-{at_3} '),
                    re.compile(r'\A' + f'{at_3}-{at_2}-{at_1} '),
                ]

                angle_found = False

                for line in ff_data_lines:
                    matches = [re.search(p, line) for p in patterns]
                    if any(matches):
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
                f_itp.write(cur_str)

            # proper dihedrals

            f_itp.write('\n[ dihedrals ] ; propers\n')
            cur_str = ';    i      j      k      l    funct'
            cur_str += '    phase     k_d      n\n'
            f_itp.write(cur_str)
            cur_str = ';                                        '
            cur_str += 'C0         C1         C2         C3         C4         C5\n'
            f_itp.write(cur_str)

            for i, j, k, l in dihedral_indices:
                at_1 = atom_types[i]
                at_2 = atom_types[j]
                at_3 = atom_types[k]
                at_4 = atom_types[l]
                patterns = [
                    re.compile(r'\A' + f'{at_1}-{at_2}-{at_3}-{at_4} '),
                    re.compile(r'\A' + f'{at_4}-{at_3}-{at_2}-{at_1} '),
                ]

                dihedral_found = False

                dihedral_ff_lines = []
                for line in ff_data_lines:
                    matches = [re.search(p, line) for p in patterns]
                    if any(matches):
                        if '.' not in line[11:].strip().split()[0]:
                            dihedral_ff_lines.append(line)
                            dihedral_found = True

                if not dihedral_found:
                    patterns = [
                        re.compile(r'\A' + f'X -{at_2}-{at_3}-X  '),
                        re.compile(r'\A' + f'X -{at_3}-{at_2}-X  '),
                    ]

                    dihedral_ff_lines = []
                    for line in ff_data_lines:
                        matches = [re.search(p, line) for p in patterns]
                        if any(matches):
                            if '.' not in line[11:].strip().split()[0]:
                                dihedral_ff_lines.append(line)
                                dihedral_found = True

                errmsg = 'ForceFieldGenerator: proper dihedral'
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
                    f_itp.write(cur_str)

            # improper dihedrals

            sp2_atom_types = [
                'c ', 'cs', 'c2', 'ca', 'cp', 'cq', 'cc', 'cd', 'ce', 'cf',
                'cu', 'cv', 'cz', 'n ', 'n2', 'na', 'nb', 'nc', 'nd', 'ne',
                'nf', 'pb', 'pc', 'pd', 'pe', 'pf'
            ]

            f_itp.write('\n[ dihedrals ] ; impropers\n')
            cur_str = ';    i      j      k      l    funct'
            cur_str += '    phase     k_d      n\n'
            f_itp.write(cur_str)

            for i, j, k in angle_indices:
                at_1 = atom_types[i]
                at_2 = atom_types[j]
                at_3 = atom_types[k]

                if at_2 not in sp2_atom_types:
                    continue

                for l in range(n_atoms):
                    if (l in [i, j, k]) or (not connected[l, j]):
                        continue
                    at_4 = atom_types[l]

                    patterns = [
                        re.compile(r'\A' + f'{at_4}-{at_1}-{at_2}-{at_3} '),
                        re.compile(r'\A' + f'{at_1}-{at_4}-{at_2}-{at_3} '),
                        re.compile(r'\A' + f'{at_3}-{at_4}-{at_2}-{at_1} '),
                        re.compile(r'\A' + f'{at_4}-{at_3}-{at_2}-{at_1} '),
                        re.compile(r'\A' + f'{at_1}-{at_3}-{at_2}-{at_4} '),
                        re.compile(r'\A' + f'{at_3}-{at_1}-{at_2}-{at_4} '),
                    ]

                    dihedral_found = False

                    for line in ff_data_lines:
                        matches = [re.search(p, line) for p in patterns]
                        if any(matches):
                            if '.' in line[11:].strip().split()[0]:
                                dihedral_ff = line[11:].strip().split()
                                dihedral_found = True
                                break

                    if not dihedral_found:
                        patterns = [
                            re.compile(r'\A' + f'X -{at_1}-{at_2}-{at_3} '),
                            re.compile(r'\A' + f'X -{at_4}-{at_2}-{at_3} '),
                            re.compile(r'\A' + f'X -{at_4}-{at_2}-{at_1} '),
                            re.compile(r'\A' + f'X -{at_3}-{at_2}-{at_1} '),
                            re.compile(r'\A' + f'X -{at_3}-{at_2}-{at_4} '),
                            re.compile(r'\A' + f'X -{at_1}-{at_2}-{at_4} '),
                        ]

                        for line in ff_data_lines:
                            matches = [re.search(p, line) for p in patterns]
                            if any(matches):
                                if '.' in line[11:].strip().split()[0]:
                                    dihedral_ff = line[11:].strip().split()
                                    dihedral_found = True
                                    break

                    if not dihedral_found:
                        patterns = [
                            re.compile(r'\A' + f'X -X -{at_2}-{at_3} '),
                            re.compile(r'\A' + f'X -X -{at_2}-{at_1} '),
                            re.compile(r'\A' + f'X -X -{at_2}-{at_4} '),
                        ]

                        for line in ff_data_lines:
                            matches = [re.search(p, line) for p in patterns]
                            if any(matches):
                                if '.' in line[11:].strip().split()[0]:
                                    dihedral_ff = line[11:].strip().split()
                                    dihedral_found = True
                                    break

                    errmsg = 'ForceFieldGenerator: improper dihedral'
                    errmsg += f' {at_4}-{at_3}-{at_2}-{at_1} is not available.'
                    assert_msg_critical(dihedral_found, errmsg)

                    barrier = float(dihedral_ff[0]) * 4.184
                    phase = float(dihedral_ff[2])
                    periodicity = abs(int(dihedral_ff[3]))

                    assert_msg_critical(
                        phase == 180.0,
                        'ForceFieldGenerator: invalid improper dihedral phase')
                    assert_msg_critical(
                        periodicity == 2,
                        'ForceFieldGenerator: invalid improper dihedral periodicity'
                    )

                    cur_str = '{:6}{:7}{:7}{:7}'.format(l + 1, i + 1, j + 1,
                                                        k + 1)
                    cur_str += '{:7}{:11.2f}{:11.5f}{:4}'.format(
                        1, phase, barrier, periodicity)
                    cur_str += ' ; {}-{}-{}-{}\n'.format(
                        atom_names[l], atom_names[i], atom_names[j],
                        atom_names[k])
                    f_itp.write(cur_str)

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

        inp_dir = Path(self.molecule_name).parent
        mol_name = Path(self.molecule_name).stem

        # Ryckaert-Bellemans function

        def rbpot(phi, c0, c1, c2, c3, c4, c5):
            v = c0 + c1 * np.cos((180 - phi) * 2 * np.pi / 360)
            v += c2 * np.cos((180 - phi) * 2 * np.pi / 360)**2
            v += c3 * np.cos((180 - phi) * 2 * np.pi / 360)**3
            v += c4 * np.cos((180 - phi) * 2 * np.pi / 360)**4
            v += c5 * np.cos((180 - phi) * 2 * np.pi / 360)**5
            return v

        itp_fname = self.get_included_file(top_filename)

        new_itp_fname = str(inp_dir / f'{mol_name}_new.itp')
        new_top_fname = str(inp_dir / f'{mol_name}_new.top')

        if itp_fname != new_itp_fname:
            self.copy_file(Path(itp_fname), Path(new_itp_fname))
            self.write_top(new_top_fname, new_itp_fname)

        for i, (dih, geom, dft_scan, angles) in enumerate(
                zip(self.target_dihedrals, self.scan_geometries,
                    self.scan_energies, self.scan_dih_angles)):

            output_dir = inp_dir / f'{mol_name}_dih_corr' / f'{i+1}'
            output_dir.mkdir(parents=True, exist_ok=True)

            zero_itp_file = output_dir / f'{mol_name}_zero.itp'
            zero_top_file = output_dir / f'{mol_name}_zero.top'

            self.copy_file(Path(new_itp_fname), zero_itp_file)
            self.write_top(str(zero_top_file), str(zero_itp_file))

            # MM scan with dihedral parameters set to zero

            self.set_dihedral_parameters(str(zero_itp_file), dih, [0.] * 6)
            mm_zero_scan = self.perform_mm_scan(str(zero_top_file), dih, geom,
                                                angles)

            # fitting Ryckaert-Bellemans function to difference
            # between MM zero scan and QM scan

            rel_e_dft = np.array(dft_scan) - min(dft_scan)
            rel_e_dft *= hartree_in_kcalpermol() * 4.184
            rel_e_mm = np.array(mm_zero_scan) - min(mm_zero_scan)

            difference = rel_e_dft - rel_e_mm
            initial_coef = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            coef, cv = curve_fit(rbpot, angles, difference, initial_coef)

            # updating parameters

            self.set_dihedral_parameters(new_itp_fname, dih, coef.tolist())

    def validate_force_field(self, top_file):
        """
        Validates force field by RMSD of dihedral potentials.

        :param top_file:
            The topology.
        :param dihedrals:
            The dihedrals.
        """

        for i, dih in enumerate(self.target_dihedrals):
            geom = self.scan_geometries[i]
            angles = self.scan_dih_angles[i]
            mm_scan = self.perform_mm_scan(top_file, dih, geom, angles)

            dft_scan = np.array(self.scan_energies[i]) - min(
                self.scan_energies[i])
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

        # select scan angles and geometries from QM data

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
            self.write_top(local_top_fname, local_itp_fname)
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
        box_length_in_nm = self.box_length * 0.1
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

        plt.plot(angles, dft_scan_kcal, '-o', label="QM")
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

        :return:
            A 2d array containing the connectivity information of the molecule.
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

        itp_file = None

        with open(top_file, 'r') as top:
            pattern = re.compile(r'\A#include')
            for line in top:
                if re.search(pattern, line):
                    itp_file = Path(top_file).parent / line.split('"')[1]

        assert_msg_critical(
            itp_file is not None,
            'ForceFieldGenerator.get_included_file: could not find ' +
            'included file')

        return str(itp_file)

    def get_atom_names(self):
        """
        Gets the atom names.

        :return:
            A list of atom names.
        """

        atom_names = []
        counter = {}

        for label in self.molecule.get_labels():
            if label not in counter:
                counter[label] = 0
            counter[label] += 1
            atom_names.append(label + str(counter[label]))

        for i, label in enumerate(self.molecule.get_labels()):
            if counter[label] == 1:
                atom_names[i] = label

        return atom_names
