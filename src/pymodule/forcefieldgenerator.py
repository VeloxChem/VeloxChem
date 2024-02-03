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

from mpi4py import MPI
from pathlib import Path
import numpy as np
import tempfile
import sys
import re

from .veloxchemlib import mpi_master, bohr_in_angstrom, hartree_in_kcalpermol
from .atomtypeidentifier import AtomTypeIdentifier
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .outputstream import OutputStream
from .respchargesdriver import RespChargesDriver
from .openmmdriver import OpenMMDriver
from .openmmgradientdriver import OpenMMGradientDriver
from .optimizationdriver import OptimizationDriver
from .inputparser import parse_input, get_random_string_parallel
from .errorhandler import assert_msg_critical
from .seminario import Seminario
from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver


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
        - gromacs_include_path: The path to topology folder of Gromacs.
        - scan_dih_angles: The dihedral angles for MM scans (list of list).
        - scan_energies: The energies from QM scans (list of list).
        - scan_geometries: The optimized geometries from QM scan (list of list).
        - target_dihedrals: The target dihedral angles for parameterization.
        - ffversion: The version of the force field.
        - workdir: The working directory.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes force field generator.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # molecule
        self.molecule_name = 'veloxchem_ff_' + get_random_string_parallel(
            self.comm)
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
        self.force_field_data = None
        self.force_field_data_extension = None

        # number of rounds for fitting dihedral potentials
        self.n_rounds = 3

        self.partial_charges = None
        self.original_top_file = None

        # MM settings
        self.gromacs_include_path = None

        # scan settings
        self.scan_dih_angles = None
        self.scan_energies = None
        self.scan_geometries = None
        self.target_dihedrals = None
        self.ffversion = 0

        self._workdir = None

        # resp settings
        self.resp_dict = None

        self.keep_files = True

    @property
    def workdir(self):
        """
        Getter function for protected workdir attribute.
        """

        return self._workdir

    @workdir.setter
    def workdir(self, value):
        """
        Setter function for protected workdir attribute.
        """

        self._workdir = value

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
            'workdir': 'str',
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
            'gromacs_include_path': 'str',
            'force_field_data': 'str',
            'force_field_data_extension': 'str',
            'n_rounds': 'int',
            'partial_charges': 'seq_fixed',
            'original_top_file': 'str',
            'keep_files': 'bool',
        }

        parse_input(self, ffg_keywords, ffg_dict)

        if 'filename' in ffg_dict and 'molecule_name' not in ffg_dict:
            self.molecule_name = ffg_dict['filename']

        if self.force_field_data is not None:
            force_field_file = Path(
                self.molecule_name).parent / self.force_field_data
            if force_field_file.is_file():
                self.force_field_data = str(force_field_file)
                assert_msg_critical(
                    'gaff' in Path(self.force_field_data).name.lower(),
                    'ForceFieldGenerator: unrecognized force field ' +
                    f'{self.force_field_data}. Only GAFF is supported.')
                if self.force_field_data_extension is None:
                    ff_file = Path(self.force_field_data)
                    self.force_field_data_extension = str(
                        ff_file.parent / (ff_file.stem + '_extension.dat'))

        if resp_dict is None:
            resp_dict = {}
        self.resp_dict = dict(resp_dict)

    def compute(self, molecule, basis):
        """
        Runs force field optimization.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        """

        # atom type identification

        self.create_topology(molecule, basis)

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

        # read QM scan

        self.ostream.print_blank()
        title = 'Force Field Generator'
        self.ostream.print_header(title)
        self.ostream.print_header('=' * (len(title) + 2))
        self.ostream.print_blank()

        self.scan_dih_angles = []
        self.scan_energies = []
        self.scan_geometries = []
        self.target_dihedrals = []

        inp_dir = Path(self.molecule_name).parent
        self.read_qm_scan_xyz_files(self.scan_xyz_files, inp_dir)

        mol_name = Path(self.molecule_name).stem

        self.ffversion = 0

        use_temp_dir = (self._workdir is None)

        if use_temp_dir:
            try:
                temp_dir = tempfile.TemporaryDirectory(
                    ignore_cleanup_errors=True)
            except TypeError:
                temp_dir = tempfile.TemporaryDirectory()
            self._workdir = Path(temp_dir.name)

        if self.original_top_file is None:
            original_itp_file = self._workdir / (mol_name +
                                                 f'_{self.ffversion:02d}.itp')
            original_top_file = original_itp_file.with_suffix('.top')

            if self.rank == mpi_master():
                self.write_itp(original_itp_file)
                self.write_top(original_top_file, original_itp_file)
            self.comm.barrier()

        else:
            original_top_file = Path(self.original_top_file)

        self.ostream.print_info(
            f'Original topology file: {original_top_file.name}')
        self.ostream.print_blank()
        self.ostream.flush()

        # validate original force field

        for i, dih in enumerate(self.target_dihedrals):
            self.validate_force_field(original_top_file, i)

        # fit dihedral potentials

        n_rounds = self.n_rounds if len(self.scan_dih_angles) > 1 else 1
        target_top_file = original_top_file
        for i_round in range(n_rounds):
            for i, dih in enumerate(self.target_dihedrals):
                target_top_file = self.dihedral_correction(target_top_file, i)

        # validate final force field

        for i, dih in enumerate(self.target_dihedrals):
            self.validate_force_field(target_top_file, i)

        # save output files

        if self.rank == mpi_master() and self.keep_files:
            out_dir = Path(self.molecule_name + '_files')
            for ffver in range(self.ffversion + 1):
                for ftype in ['itp', 'top']:
                    fname = mol_name + f'_{ffver:02d}.{ftype}'
                    self.copy_file(self._workdir / fname, out_dir / fname)
                    valstr = f'Saving file: {str(out_dir / fname)}'
                    self.ostream.print_info(valstr)

            self.ostream.print_blank()
            self.ostream.flush()

        if use_temp_dir:
            try:
                temp_dir.cleanup()
            except (NotADirectoryError, PermissionError):
                pass
            self._workdir = None

    def read_qm_scan_xyz_files(self, scan_xyz_files, inp_dir=None):
        """
        Reads QM scan xyz files.

        :param scan_xyz_files:
            The list of xyz files from QM scan.
        """

        if self.scan_dih_angles is None:
            self.scan_dih_angles = []

        if self.scan_energies is None:
            self.scan_energies = []

        if self.scan_geometries is None:
            self.scan_geometries = []

        if self.target_dihedrals is None:
            self.target_dihedrals = []

        # reading QM data

        if inp_dir is None:
            inp_dir = Path('.')
        elif isinstance(inp_dir, str):
            inp_dir = Path(inp_dir)

        self.ostream.print_info('Reading QM scan from file...')

        for xyz in scan_xyz_files:
            xyz_fname = str(inp_dir / xyz)

            self.ostream.print_info(f'  {xyz_fname}')
            self.ostream.flush()

            geometries = []
            energies = []
            dih_angles = []

            # read geometries

            xyz_lines = None
            if self.rank == mpi_master():
                with open(xyz_fname, 'r') as f_xyz:
                    xyz_lines = f_xyz.readlines()
            xyz_lines = self.comm.bcast(xyz_lines, root=mpi_master())

            n_atoms = int(xyz_lines[0].split()[0])
            n_geoms = len(xyz_lines) // (n_atoms + 2)

            for i_geom in range(n_geoms):
                i_start = i_geom * (n_atoms + 2)
                i_end = i_start + (n_atoms + 2)

                assert_msg_critical(
                    int(xyz_lines[i_start].split()[0]) == n_atoms,
                    'ForceFieldGenerator.read_qm_scan_xyz_files: ' +
                    'inconsistent number of atoms')

                mol_str = ''.join(xyz_lines[i_start + 2:i_end])
                geometries.append(
                    Molecule.read_molecule_string(mol_str, units='angstrom'))

            self.scan_geometries.append(geometries)

            # read energies and dihedral angles

            pattern = re.compile(r'\AScan')

            for line in xyz_lines:
                if re.search(pattern, line):
                    energies.append(float(line.split('Energy')[1].split()[0]))
                    dih_angles.append(float(line.split('=')[1].split()[0]))
                    dih_inds = [
                        int(i)
                        for i in line.split('Dihedral')[1].split()[0].split('-')
                    ]
            self.scan_energies.append(energies)
            self.scan_dih_angles.append(dih_angles)
            self.target_dihedrals.append(dih_inds)

        self.ostream.print_blank()

    def update_dihedral_range(self, dih_start_end, i):
        """
        Updates the range of dihedral angles for i-th dihedral (0-based index).

        :param dih_start_end:
            A tuple containing the starting and ending values of dihedral angle.
        :param i:
            The index of the target dihedral.
        """

        new_scan_dih_angles = []
        new_scan_energies = []
        new_scan_geometries = []

        for angle, ene, geom in zip(
                self.scan_dih_angles[i],
                self.scan_energies[i],
                self.scan_geometries[i],
        ):
            if angle >= dih_start_end[0] and angle <= dih_start_end[1]:
                new_scan_dih_angles.append(angle)
                new_scan_energies.append(ene)
                new_scan_geometries.append(geom)

        self.scan_dih_angles[i] = new_scan_dih_angles
        self.scan_energies[i] = new_scan_energies
        self.scan_geometries[i] = new_scan_geometries

    @staticmethod
    def get_gaff_data_lines():
        """
        Reads GAFF data lines into a list.
        """

        from urllib.request import urlopen

        openmmff_commit = 'b3e92a373c80bfb8fd791e4a72beafc035fcc722'
        gaff_url = (
            'https://raw.githubusercontent.com/openmm/openmmforcefields/' +
            openmmff_commit + '/openmmforcefields/ffxml/amber/gaff/dat/' +
            'gaff-2.11.dat')

        with urlopen(gaff_url) as f_gaff:
            content = f_gaff.read().decode('utf-8')

        return content.splitlines()

    def create_topology(self, molecule, basis=None):
        """
        Analizes the topology of the molecule and create dictionaries
        for the atoms, bonds, angles, dihedrals, impropers and pairs.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        """

        # Read the force field data lines.

        if self.force_field_data is None:
            ff_data_lines = self.get_gaff_data_lines()
        else:
            with open(self.force_field_data, 'r') as ff_data:
                ff_data_lines = ff_data.readlines()

        if ff_data_lines[0].startswith('AMBER General Force Field'):
            ff_data_version = ff_data_lines[0].split('Version')[1].split()[0]
            ff_data_version = ff_data_version.replace(',', '')
        else:
            ff_data_version = None

        # Molecular information

        self.molecule = molecule

        coords = self.molecule.get_coordinates_in_angstrom()
        n_atoms = self.molecule.number_of_atoms()

        atomtypeidentifier = AtomTypeIdentifier(self.comm)
        atomtypeidentifier.ostream.mute()

        self.atom_types = atomtypeidentifier.generate_gaff_atomtypes(
            self.molecule)
        atomtypeidentifier.identify_equivalences()

        self.connectivity_matrix = np.copy(
            atomtypeidentifier.connectivity_matrix)

        if self.partial_charges is None:
            if basis is None:
                if self.rank == mpi_master():
                    basis = MolecularBasis.read(self.molecule,
                                                '6-31G*',
                                                ostream=None)
                else:
                    basis = MolecularBasis()
                basis.broadcast(self.rank, self.comm)
                msg = 'Using 6-31G* basis set for RESP charges...'
                self.ostream.print_info(msg)
                self.ostream.flush()

            resp_drv = RespChargesDriver(self.comm)
            resp_drv.filename = self.molecule_name
            if self.resp_dict is not None:
                resp_drv.update_settings(self.resp_dict)
            if resp_drv.equal_charges is None:
                resp_drv.equal_charges = atomtypeidentifier.equivalent_charges

            self.partial_charges = resp_drv.compute(self.molecule, basis,
                                                    'resp')
            self.partial_charges = self.comm.bcast(self.partial_charges,
                                                   root=mpi_master())

        # preparing atomtypes and atoms

        assert_msg_critical(
            len(self.atom_types) == n_atoms,
            'ForceFieldGenerator: inconsistent atom_types')

        for i in range(n_atoms):
            self.atom_types[i] = f'{self.atom_types[i].strip():<2s}'
        self.unique_atom_types = sorted(list(set(self.atom_types)))

        # Bonds

        bond_indices = set()
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                if self.connectivity_matrix[i, j] == 1:
                    bond_indices.add((i, j))
        bond_indices = sorted(list(bond_indices))

        # Angles

        angle_indices = set()

        for i, j in bond_indices:
            for k in range(n_atoms):
                if k in [i, j]:
                    continue
                if self.connectivity_matrix[j, k] == 1:
                    inds = (i, j, k) if i < k else (k, j, i)
                    angle_indices.add(inds)
                if self.connectivity_matrix[k, i] == 1:
                    inds = (k, i, j) if k < j else (j, i, k)
                    angle_indices.add(inds)
        angle_indices = sorted(list(angle_indices))

        # Dihedrals

        dihedral_indices = set()

        for i, j, k in angle_indices:
            for l in range(n_atoms):
                if l in [i, j, k]:
                    continue
                if self.connectivity_matrix[k, l] == 1:
                    inds = (i, j, k, l) if i < l else (l, k, j, i)
                    dihedral_indices.add(inds)
                if self.connectivity_matrix[l, i] == 1:
                    inds = (l, i, j, k) if l < k else (k, j, i, l)
                    dihedral_indices.add(inds)
        dihedral_indices = sorted(list(dihedral_indices))

        # Exclusions

        exclusion_indices = []
        if self.nrexcl >= 2:
            for i, j in bond_indices:
                exclusion_indices.append((i, j))
        if self.nrexcl >= 3:
            for i, j, k in angle_indices:
                exclusion_indices.append((i, k))

        # 1-4 pairs

        pairs_14 = set()
        for i, j, k, l in dihedral_indices:
            if (i, l) not in exclusion_indices:
                pairs_14.add((i, l))
        pairs_14 = sorted(list(pairs_14))

        # Read the force field and include the data in the topology dictionary.

        # Atomtypes analysis

        atom_type_params = {}

        for at in self.unique_atom_types:
            atom_type_found = False

            for line in ff_data_lines:
                if line.startswith(f'  {at}     '):
                    atom_ff = line[5:].strip().split()
                    sigma = float(atom_ff[0]) * 2**(-1 / 6) * 2 / 10
                    epsilon = float(atom_ff[1]) * 4.184
                    comment = 'GAFF'
                    atom_type_found = True
                    break

            if not atom_type_found:
                warnmsg = f'ForceFieldGenerator: atom type {at} is not in GAFF.'
                self.ostream.print_warning(warnmsg)
                sigma, epsilon, comment = 0.0, 0.0, 'Unknown'

            atom_type_params[at] = {
                'sigma': sigma,
                'epsilon': epsilon,
                'comment': comment
            }

        # Atoms analysis

        self.atoms = {}

        atom_names = self.get_atom_names()
        atom_masses = self.molecule.masses_to_numpy()
        equivalent_atoms = list(atomtypeidentifier.equivalent_atoms)

        for i in range(n_atoms):
            at = self.atom_types[i]
            self.atoms[i] = {
                'type': at,
                'name': atom_names[i],
                'mass': atom_masses[i],
                'charge': self.partial_charges[i],
                'sigma': atom_type_params[at]['sigma'],
                'epsilon': atom_type_params[at]['epsilon'],
                'equivalent_atom': equivalent_atoms[i],
            }

        # Bonds analysis

        self.bonds = {}

        for i, j in bond_indices:

            r_eq = np.linalg.norm(coords[i] - coords[j]) * 0.1

            at_1 = self.atom_types[i]
            at_2 = self.atom_types[j]
            patterns = [
                re.compile(r'\A' + f'{at_1}-{at_2}  '),
                re.compile(r'\A' + f'{at_2}-{at_1}  '),
            ]

            bond_found = False
            r, k_r, comment = None, None, None

            for line in ff_data_lines:
                for p in patterns:
                    m = re.search(p, line)
                    if m is not None:
                        bond_ff = line[5:].strip().split()
                        r = float(bond_ff[1]) * 0.1
                        k_r = float(bond_ff[0]) * 4.184 * 2 * 100
                        comment = m.group(0)
                        bond_found = True
                        break

            if not bond_found:
                warnmsg = f'ForceFieldGenerator: bond {at_1}-{at_2}'
                warnmsg += ' is not available'
                self.ostream.print_warning(warnmsg)
                # Default value for bonds
                r, k_r, comment = r_eq, 2.5e+5, 'Unknown'

            if self.eq_param:
                if abs(r - r_eq) > self.r_thresh:
                    msg = f'Updated bond length {i+1}-{j+1} '
                    msg += f'({at_1}-{at_2}) to {r_eq:.3f} nm'
                    self.ostream.print_info(msg)
                    self.ostream.flush()
                r = r_eq

            self.bonds[(i, j)] = {
                'force_constant': k_r,
                'equilibrium': r,
                'comment': comment
            }

        # Pairs writing

        self.pairs = {}

        for i, j in pairs_14:

            self.pairs[(i, j)] = {'comment': None}

        # Angles analysis

        self.angles = {}

        for i, j, k in angle_indices:

            a = coords[i] - coords[j]
            b = coords[k] - coords[j]
            theta_eq = np.arccos(
                np.dot(a, b) / np.linalg.norm(a) /
                np.linalg.norm(b)) * 180 / np.pi

            at_1 = self.atom_types[i]
            at_2 = self.atom_types[j]
            at_3 = self.atom_types[k]
            patterns = [
                re.compile(r'\A' + f'{at_1}-{at_2}-{at_3} '),
                re.compile(r'\A' + f'{at_3}-{at_2}-{at_1} '),
            ]

            angle_found = False
            theta, k_theta, comment = None, None, None

            for line in ff_data_lines:
                for p in patterns:
                    m = re.search(p, line)
                    if m is not None:
                        angle_ff = line[8:].strip().split()
                        theta = float(angle_ff[1])
                        k_theta = float(angle_ff[0]) * 4.184 * 2
                        comment = m.group(0)
                        angle_found = True
                        break

            if not angle_found:
                warnmsg = f'ForceFieldGenerator: angle {at_1}-{at_2}-{at_3}'
                warnmsg += ' is not available.'
                self.ostream.print_warning(warnmsg)
                # Default value for angles
                theta, k_theta, comment = theta_eq, 1000, 'Unknown'

            if self.eq_param:
                if abs(theta - theta_eq) > self.theta_thresh:
                    msg = f'Updated bond angle {i+1}-{j+1}-{k+1} '
                    msg += f'({at_1}-{at_2}-{at_3}) to {theta_eq:.3f} deg'
                    self.ostream.print_info(msg)
                    self.ostream.flush()
                theta = theta_eq

            self.angles[(i, j, k)] = {
                'force_constant': k_theta,
                'equilibrium': theta,
                'comment': comment
            }

        # Dihedrals analysis

        self.dihedrals = {}

        for i, j, k, l in dihedral_indices:

            at_1 = self.atom_types[i]
            at_2 = self.atom_types[j]
            at_3 = self.atom_types[k]
            at_4 = self.atom_types[l]

            patterns = [
                re.compile(r'\A' + f'{at_1}-{at_2}-{at_3}-{at_4} '),
                re.compile(r'\A' + f'{at_4}-{at_3}-{at_2}-{at_1} '),
            ]

            dihedral_found = False

            dihedral_ff_lines = []
            dihedral_matches = []
            for line in ff_data_lines:
                for p in patterns:
                    m = re.search(p, line)
                    if m is not None:
                        dihedral_ff = line[11:60].strip().split()
                        if len(dihedral_ff) == 4:
                            dihedral_ff_lines.append(line)
                            dihedral_matches.append(m.group(0))
                            dihedral_found = True
                            break

            if not dihedral_found:
                patterns = [
                    re.compile(r'\A' + f'X -{at_2}-{at_3}-X  '),
                    re.compile(r'\A' + f'X -{at_3}-{at_2}-X  '),
                ]

                dihedral_ff_lines = []
                dihedral_matches = []
                for line in ff_data_lines:
                    for p in patterns:
                        m = re.search(p, line)
                        if m is not None:
                            dihedral_ff = line[11:60].strip().split()
                            if len(dihedral_ff) == 4:
                                dihedral_ff_lines.append(line)
                                dihedral_matches.append(m.group(0))
                                dihedral_found = True
                                break

            if not dihedral_found:
                warnmsg = f'ForceFieldGenerator: dihedral {at_1}-{at_2}-{at_3}-{at_4}'
                warnmsg += ' is not available.'
                self.ostream.print_warning(warnmsg)
                # Default value for dihedrals
                self.dihedrals[(i, j, k, l)] = [{
                    'barrier': 0.0,
                    'phase': 0.0,
                    'periodicity': 1,
                    'comment': 'Unknown'
                }]
            else:
                self.dihedrals[(i, j, k, l)] = []

            for line, comment in zip(dihedral_ff_lines, dihedral_matches):
                dihedral_ff = line[11:60].strip().split()

                multiplicity = int(dihedral_ff[0])
                barrier = float(dihedral_ff[1]) * 4.184 / multiplicity
                phase = float(dihedral_ff[2])
                # Note: negative periodicity implies multitermed dihedral
                # See https://ambermd.org/FileFormats.php
                try:
                    periodicity = int(dihedral_ff[3])
                except ValueError:
                    periodicity = int(float(dihedral_ff[3]))

                # TODO: use RB type for multitermed dihedral

                self.dihedrals[(i, j, k, l)].append({
                    'type': 'Fourier',
                    'barrier': barrier,
                    'phase': phase,
                    'periodicity': periodicity,
                    'comment': comment
                })

                if periodicity > 0:
                    break

        # Impropers

        self.impropers = {}

        sp2_atom_types = [
            'c ', 'cs', 'c2', 'ca', 'cp', 'cq', 'cc', 'cd', 'ce', 'cf', 'cu',
            'cv', 'cz', 'n ', 'n2', 'na', 'nb', 'nc', 'nd', 'ne', 'nf', 'pb',
            'pc', 'pd', 'pe', 'pf'
        ]

        improper_atom_inds = []

        for i, j, k in angle_indices:
            at_1 = self.atom_types[i]
            at_2 = self.atom_types[j]
            at_3 = self.atom_types[k]

            if at_2 not in sp2_atom_types:
                continue

            if j not in improper_atom_inds:
                improper_atom_inds.append(j)
            else:
                continue

            for l in range(n_atoms):
                if (l in [i, j, k]) or (self.connectivity_matrix[l, j] != 1):
                    continue
                at_4 = self.atom_types[l]

                patterns = [
                    re.compile(r'\A' + f'{at_4}-{at_1}-{at_2}-{at_3} '),
                    re.compile(r'\A' + f'{at_4}-{at_3}-{at_2}-{at_1} '),
                    re.compile(r'\A' + f'{at_1}-{at_3}-{at_2}-{at_4} '),
                    re.compile(r'\A' + f'{at_1}-{at_4}-{at_2}-{at_3} '),
                    re.compile(r'\A' + f'{at_3}-{at_1}-{at_2}-{at_4} '),
                    re.compile(r'\A' + f'{at_3}-{at_4}-{at_2}-{at_1} '),
                ]

                dihedral_found = False
                barrier, phase, periodicity, comment = None, None, None, None

                for line in ff_data_lines:
                    for p in patterns:
                        m = re.search(p, line)
                        if m is not None:
                            dihedral_ff = line[11:60].strip().split()
                            if len(dihedral_ff) == 3:
                                barrier = float(dihedral_ff[0]) * 4.184
                                phase = float(dihedral_ff[1])
                                periodicity = int(float(dihedral_ff[2]))
                                comment = m.group(0)
                                dihedral_found = True
                                break

                if not dihedral_found:
                    patterns = [
                        re.compile(r'\A' + f'X -{at_1}-{at_2}-{at_3} '),
                        re.compile(r'\A' + f'X -{at_3}-{at_2}-{at_1} '),
                        re.compile(r'\A' + f'X -{at_3}-{at_2}-{at_4} '),
                        re.compile(r'\A' + f'X -{at_4}-{at_2}-{at_3} '),
                        re.compile(r'\A' + f'X -{at_1}-{at_2}-{at_4} '),
                        re.compile(r'\A' + f'X -{at_4}-{at_2}-{at_1} '),
                    ]

                    for line in ff_data_lines:
                        for p in patterns:
                            m = re.search(p, line)
                            if m is not None:
                                dihedral_ff = line[11:60].strip().split()
                                if len(dihedral_ff) == 3:
                                    barrier = float(dihedral_ff[0]) * 4.184
                                    phase = float(dihedral_ff[1])
                                    periodicity = int(float(dihedral_ff[2]))
                                    comment = m.group(0)
                                    dihedral_found = True
                                    break

                if not dihedral_found:
                    patterns = [
                        re.compile(r'\A' + f'X -X -{at_2}-{at_3} '),
                        re.compile(r'\A' + f'X -X -{at_2}-{at_1} '),
                        re.compile(r'\A' + f'X -X -{at_2}-{at_4} '),
                    ]

                    for line in ff_data_lines:
                        for p in patterns:
                            m = re.search(p, line)
                            if m is not None:
                                dihedral_ff = line[11:60].strip().split()
                                if len(dihedral_ff) == 3:
                                    barrier = float(dihedral_ff[0]) * 4.184
                                    phase = float(dihedral_ff[1])
                                    periodicity = int(float(dihedral_ff[2]))
                                    comment = m.group(0)
                                    dihedral_found = True
                                    break

                if not dihedral_found:
                    warnmsg = 'ForceFieldGenerator: '
                    warnmsg += f'improper {at_1}-{at_2}-{at_3}-{at_4} '
                    warnmsg += 'is not available.'
                    self.ostream.print_warning(warnmsg)
                    # Default values for impropers
                    barrier, phase, periodicity = 1.1 * 4.184, 180.0, 2
                    comment = 'Unknown'

                assert_msg_critical(
                    phase == 180.0,
                    'ForceFieldGenerator: invalid improper dihedral phase')

                assert_msg_critical(
                    periodicity == 2,
                    'ForceFieldGenerator: invalid improper dihedral periodicity'
                )

            self.impropers[(i, j, k, l)] = {
                'barrier': barrier,
                'phase': phase,
                'periodicity': periodicity,
                'comment': comment
            }

    def reparameterize(self,
                       hessian=None,
                       reparameterize_all=False,
                       reparameterize_keys=None):
        """
        Reparameterizes all unknown parameters with the Seminario method using
        the given Hessian matrix.

        :param hessian:
            The Hessian matrix, or the method to generate Hessian.
        :param reparameterize_all:
            If True, all parameters are reparameterized. If False, only unknown
            parameters are reparameterized.
        :param reparameterize_keys:
            List of specific keys to reparameterize, can be bonds and angles.
        """

        # Hessian matrix

        if hessian is None:
            # TODO: generate Hessian using VeloxChem
            assert_msg_critical(
                False, 'ForceFieldGenerator.reparameterize: expecting Hessian')

        elif isinstance(hessian, str):
            assert_msg_critical(
                hessian.lower() == 'xtb',
                'ForceFieldGenerator.reparameterize: invalid Hessian option')

            # XTB optimization
            self.ostream.print_info('Optimizing molecule using XTB...')
            self.ostream.flush()
            xtb_drv = XtbDriver(self.comm)
            xtb_drv.mute()
            xtb_grad_drv = XtbGradientDriver(self.comm)
            xtb_grad_drv.ostream.state = False
            xtb_opt_drv = OptimizationDriver(xtb_grad_drv)
            xtb_opt_drv.filename = self.molecule_name
            self.molecule = xtb_opt_drv.compute(self.molecule, xtb_drv)

            # XTB Hessian
            self.ostream.print_info('Computing Hessian using XTB...')
            self.ostream.flush()
            xtb_hessian_drv = XtbHessianDriver(self.comm)
            xtb_hessian_drv.ostream.state = False
            xtb_hessian_drv.compute(self.molecule, xtb_drv)
            hessian = np.copy(xtb_hessian_drv.hessian)

            self.ostream.print_blank()
            self.ostream.print_reference('Reference:')
            self.ostream.print_reference(xtb_drv.get_reference())
            self.ostream.print_blank()
            self.ostream.flush()

        elif isinstance(hessian, np.ndarray):
            natoms = self.molecule.number_of_atoms()
            assert_msg_critical(
                hessian.shape == (natoms * 3, natoms * 3),
                'ForceFieldGenerator.reparameterize: invalid Hessian matrix')

        else:
            assert_msg_critical(
                False,
                'ForceFieldGenerator.reparameterize: invalid Hessian option')

        angstrom_to_nm = 0.1  # 1 angstrom is 0.1 nm
        bohr_to_nm = bohr_in_angstrom() * angstrom_to_nm
        cal_to_joule = 4.184  # 1 calorie is 4.184 joule
        hartree_to_kJmol = hartree_in_kcalpermol() * cal_to_joule

        coords_in_au = self.molecule.get_coordinates_in_bohr()

        seminario = Seminario(hessian, coords_in_au)

        self.ostream.print_info(
            'Force-field reparameterization based on the Seminario method')
        self.ostream.print_blank()
        self.ostream.print_reference('Reference:')
        self.ostream.print_reference(seminario.get_reference())
        self.ostream.print_blank()
        self.ostream.flush()

        # Reparameterize bonds

        for i, j in self.bonds:

            if not reparameterize_all:
                if reparameterize_keys is None:
                    if self.bonds[(i, j)]['comment'].capitalize() != 'Unknown':
                        continue
                elif (i, j) not in reparameterize_keys:
                    continue

            new_equilibrium = np.linalg.norm(coords_in_au[i] -
                                             coords_in_au[j]) * bohr_to_nm

            new_force_constant = seminario.calculate_bond_force_constant(i, j)
            new_force_constant *= hartree_to_kJmol / (bohr_to_nm**2)

            self.bonds[(i, j)]['equilibrium'] = new_equilibrium
            self.bonds[(i, j)]['force_constant'] = new_force_constant
            self.bonds[(i, j)]['comment'] += ' from Hessian'

        # Average over equivalent bonds

        uniq_bonds_data = {}

        for (i, j), bond in self.bonds.items():
            eq_atoms_ij = tuple(
                sorted([
                    self.atoms[i]['equivalent_atom'],
                    self.atoms[j]['equivalent_atom']
                ]))

            if eq_atoms_ij not in uniq_bonds_data:
                uniq_bonds_data[eq_atoms_ij] = {
                    'indices': [],
                    'r_sum': 0.0,
                    'k_r_sum': 0.0,
                    'count': 0
                }

            uniq_bonds_data[eq_atoms_ij]['indices'].append((i, j))
            uniq_bonds_data[eq_atoms_ij]['r_sum'] += bond['equilibrium']
            uniq_bonds_data[eq_atoms_ij]['k_r_sum'] += bond['force_constant']
            uniq_bonds_data[eq_atoms_ij]['count'] += 1

        for atom_pair, bond_data in uniq_bonds_data.items():
            for i, j in bond_data['indices']:
                aver_r = bond_data['r_sum'] / bond_data['count']
                aver_k_r = bond_data['k_r_sum'] / bond_data['count']
                self.bonds[(i, j)]['equilibrium'] = aver_r
                self.bonds[(i, j)]['force_constant'] = aver_k_r

        # Reparameterize angles

        for i, j, k in self.angles:

            if not reparameterize_all:
                if reparameterize_keys is None:
                    if (self.angles[(i, j, k)]['comment'].capitalize()
                            != 'Unknown'):
                        continue
                elif (i, j, k) not in reparameterize_keys:
                    continue

            a = coords_in_au[i] - coords_in_au[j]
            b = coords_in_au[k] - coords_in_au[j]
            new_equilibrium = np.arccos(
                np.dot(a, b) / np.linalg.norm(a) /
                np.linalg.norm(b)) * 180 / np.pi

            new_force_constant = seminario.calculate_angle_force_constant(
                i, j, k)
            new_force_constant *= hartree_to_kJmol

            self.angles[(i, j, k)]['equilibrium'] = new_equilibrium
            self.angles[(i, j, k)]['force_constant'] = new_force_constant
            self.angles[(i, j, k)]['comment'] += ' from Hessian'

        # Average over equivalent angles

        uniq_angles_data = {}

        for (i, j, k), angle in self.angles.items():
            eq_atoms_ijk = sorted([
                self.atoms[i]['equivalent_atom'],
                self.atoms[k]['equivalent_atom']
            ])
            eq_atoms_ijk.insert(1, self.atoms[j]['equivalent_atom'])
            eq_atoms_ijk = tuple(eq_atoms_ijk)

            if eq_atoms_ijk not in uniq_angles_data:
                uniq_angles_data[eq_atoms_ijk] = {
                    'indices': [],
                    'theta_sum': 0.0,
                    'k_theta_sum': 0.0,
                    'count': 0
                }

            uniq_angles_data[eq_atoms_ijk]['indices'].append((i, j, k))
            uniq_angles_data[eq_atoms_ijk]['theta_sum'] += angle['equilibrium']
            uniq_angles_data[eq_atoms_ijk]['k_theta_sum'] += angle[
                'force_constant']
            uniq_angles_data[eq_atoms_ijk]['count'] += 1

        for atom_triple, angle_data in uniq_angles_data.items():
            for i, j, k in angle_data['indices']:
                aver_theta = angle_data['theta_sum'] / angle_data['count']
                aver_k_theta = angle_data['k_theta_sum'] / angle_data['count']
                self.angles[(i, j, k)]['equilibrium'] = aver_theta
                self.angles[(i, j, k)]['force_constant'] = aver_k_theta

    def write_top(self, top_file, itp_file):
        """
        Writes a topology file.

        :param top_file:
            The topology file.
        :param itp_file:
            The included itp file.
        """

        top_fname = str(top_file)
        itp_fname = str(itp_file)

        mol_name = Path(self.molecule_name).stem

        with open(top_fname, 'w') as f_top:

            # header

            f_top.write('; Generated by VeloxChem\n')

            # defaults

            f_top.write('\n[ defaults ]\n')
            cur_str = '; nbfunc        comb-rule       gen-pairs'
            cur_str += '       fudgeLJ fudgeQQ\n'
            f_top.write(cur_str)
            gen_pairs = 'yes' if self.gen_pairs else 'no'
            f_top.write('{}{:16}{:>18}{:19.4f}{:8.4f}\n'.format(
                self.nbfunc, self.comb_rule, gen_pairs, self.fudgeLJ,
                self.fudgeQQ))

            # include itp

            f_top.write('\n#include "' + Path(itp_fname).name + '"\n')

            # system

            f_top.write('\n[ system ]\n')
            f_top.write(' {}\n'.format(mol_name))

            # molecules

            f_top.write('\n[ molecules ]\n')
            f_top.write('; Compound        nmols\n')
            f_top.write('{:>10}{:9}\n'.format(mol_name, 1))

    def write_itp(self, itp_file, mol_name='MOL'):
        """
        Writes an ITP file with the original parameters.

        :param itp_file:
            The ITP file path.
        """

        itp_filename = str(itp_file)
        mol_name = Path(self.molecule_name).stem

        with open(itp_filename, 'w') as f_itp:
            # Header
            f_itp.write('; Generated by VeloxChem\n')

            # Atom types
            f_itp.write('\n[ atomtypes ]\n')
            line_str = ';name   bond_type     mass     charge'
            line_str += '   ptype   sigma         epsilon\n'
            f_itp.write(line_str)

            for at in self.unique_atom_types:
                for i, atom in self.atoms.items():
                    if atom['type'] == at:
                        line_str = '{:>3}{:>9}{:17.5f}{:9.5f}{:>4}'.format(
                            at, at, 0., 0., 'A')
                        line_str += '{:16.5e}{:14.5e}\n'.format(
                            atom['sigma'], atom['epsilon'])
                        f_itp.write(line_str)
                        break

            # Molecule type
            f_itp.write('\n[ moleculetype ]\n')
            f_itp.write(';name            nrexcl\n')
            f_itp.write(f'{mol_name:>10}{self.nrexcl:10}\n')

            # Atoms
            f_itp.write('\n[ atoms ]\n')
            line_str = ';   nr  type  resi  res  atom  cgnr'
            line_str += '     charge       mass\n'
            f_itp.write(line_str)

            total_charge = 0.0
            for i, atom in self.atoms.items():
                total_charge += atom['charge']
                line_str = '{:6}{:>5}{:6}{:>6}{:>6}'.format(
                    i + 1, atom['type'], 1, mol_name, atom['name'])
                line_str += '{:5}{:13.6f}{:13.5f}'.format(
                    i + 1, atom['charge'], atom['mass'])
                line_str += ' ; qtot{:7.3f}  equiv. {}\n'.format(
                    total_charge, atom['equivalent_atom'])
                f_itp.write(line_str)

            # Bonds
            f_itp.write('\n[ bonds ]\n')
            f_itp.write(';   ai     aj    funct       r           k_r\n')
            for (i, j), bond in self.bonds.items():
                line_str = '{:6}{:7}{:7}{:14.4e}{:14.4e} ; {}\n'.format(
                    i + 1, j + 1, 1, bond['equilibrium'],
                    bond['force_constant'], bond['comment'])
                f_itp.write(line_str)

            # Pairs
            f_itp.write('\n[ pairs ]\n')
            f_itp.write(';   ai     aj    funct\n')
            for i, j in self.pairs:
                f_itp.write('{:6}{:7}{:7}\n'.format(i + 1, j + 1, 1))

            # Angles
            f_itp.write('\n[ angles ]\n')
            f_itp.write(
                ';   ai     aj     ak    funct     theta       k_theta\n')
            for (i, j, k), angle in self.angles.items():
                line_str = '{:6}{:7}{:7}{:7}{:14.4e}{:14.4e} ; {}\n'.format(
                    i + 1, j + 1, k + 1, 1, angle['equilibrium'],
                    angle['force_constant'], angle['comment'])
                f_itp.write(line_str)

            # Proper dihedrals
            f_itp.write('\n[ dihedrals ]\n')
            f_itp.write('; propers\n')
            f_itp.write(
                ';   ai     aj     ak     al    funct    phase     k_d      n\n'
            )
            for (i, j, k, l), dihedral_params in self.dihedrals.items():
                for dih in dihedral_params:
                    line_str = '{:6}{:7}{:7}{:7}'.format(
                        i + 1, j + 1, k + 1, l + 1)
                    line_str += '{:7}{:11.2f}{:11.5f}{:4} ; {}\n'.format(
                        9, dih['phase'], dih['barrier'],
                        abs(dih['periodicity']), dih['comment'])
                    f_itp.write(line_str)

            # Improper dihedrals
            f_itp.write('\n[ dihedrals ]\n')
            f_itp.write('; impropers\n')
            f_itp.write(
                ';   ai     aj     ak     al    funct    phase     k_d      n\n'
            )
            for (i, j, k, l), dih in self.impropers.items():
                line_str = '{:6}{:7}{:7}{:7}'.format(l + 1, i + 1, j + 1, k + 1)
                line_str += '{:7}{:11.2f}{:11.5f}{:4} ; {}\n'.format(
                    4, dih['phase'], dih['barrier'], abs(dih['periodicity']),
                    dih['comment'])
                f_itp.write(line_str)

    def write_gro(self, gro_file, mol_name='MOL'):
        """
        Writes a GRO file with the original coordinates.

        :param gro_file:
            The GRO file path.
        """

        gro_filename = str(gro_file)
        mol_name = Path(self.molecule_name).stem

        coords_in_nm = self.molecule.get_coordinates_in_angstrom() * 0.1

        with open(gro_filename, 'w') as f_gro:
            # Header
            f_gro.write(f'GRO file of {mol_name}, generated by VeloxChem\n')
            f_gro.write(f'{len(self.atoms):>5d}\n')

            # Atoms
            for i, atom in self.atoms.items():
                atom_name = atom['name']
                line_str = f'{1:>5d}{mol_name:<5s}{atom_name:<5s}{i + 1:>5d}'
                for d in range(3):
                    line_str += f'{coords_in_nm[i][d]:12.7f}'
                line_str += '\n'
                f_gro.write(line_str)

            # Box
            box_dimension = 10.0
            line_str = f'{box_dimension:10.5f}' * 3
            f_gro.write(line_str)

    def write_gromacs_files(self, filename, mol_name=None):
        """
        Writes all the needed files for a MD simulation with GROMACS.

        :param filename:
            The name of the molecule.
        :param mol_name:
            Str. The name of the molecule.
        """

        if mol_name is None:
            self.molecule_name = 'MOL'
        else:
            self.molecule_name = mol_name

        itp_file = Path(filename).with_suffix('.itp')
        top_file = Path(filename).with_suffix('.top')
        gro_file = Path(filename).with_suffix('.gro')

        self.write_itp(itp_file)
        self.write_top(top_file, itp_file)
        self.write_gro(gro_file)

    @staticmethod
    def copy_file(src, dest):
        """
        Copies file (from src to dest).

        :param src:
            The source of copy.
        :param dest:
            The destination of copy.
        """

        if (not dest.is_file()) or (not src.samefile(dest)):
            if not dest.parent.is_dir():
                dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_text(src.read_text())

    def dihedral_correction(self, top_file, i, kT=None):
        """
        Corrects dihedral parameters.

        :param top_filename:
            The topology file.
        :param i:
            The index of the target dihedral.
        :param kT:
            kT for Boltzmann factor (used in weighted fitting).
        """

        try:
            from scipy.optimize import curve_fit
        except ImportError:
            raise ImportError('Unable to import scipy. Please install scipy ' +
                              'via pip or conda.')

        # Ryckaert-Bellemans function

        def rbpot(phi, c0, c1, c2, c3, c4, c5):
            v = c0 + c1 * np.cos((180 - phi) * 2 * np.pi / 360)
            v += c2 * np.cos((180 - phi) * 2 * np.pi / 360)**2
            v += c3 * np.cos((180 - phi) * 2 * np.pi / 360)**3
            v += c4 * np.cos((180 - phi) * 2 * np.pi / 360)**4
            v += c5 * np.cos((180 - phi) * 2 * np.pi / 360)**5
            return v

        top_fname = str(top_file)

        if self.rank == mpi_master():
            itp_fname = self.get_included_file(top_fname)
        else:
            itp_fname = None
        itp_fname = self.comm.bcast(itp_fname, root=mpi_master())

        mol_name = Path(self.molecule_name).stem

        workdir = Path('.') if self._workdir is None else self._workdir

        self.ffversion += 1
        new_itp_fname = str(workdir / f'{mol_name}_{self.ffversion:02d}.itp')
        new_top_fname = str(workdir / f'{mol_name}_{self.ffversion:02d}.top')

        if itp_fname != new_itp_fname:
            if self.rank == mpi_master():
                self.copy_file(Path(itp_fname), Path(new_itp_fname))
                self.write_top(new_top_fname, new_itp_fname)
            self.comm.barrier()

        dih = self.target_dihedrals[i]
        geom = self.scan_geometries[i]
        qm_scan = self.scan_energies[i]
        angles = self.scan_dih_angles[i]

        self.ostream.print_info('Fitting dihedral angle ' +
                                '{}-{}-{}-{}'.format(*dih) + '...')
        self.ostream.flush()

        # MM scan with dihedral parameters set to zero

        output_dir = workdir / f'{mol_name}_dih_corr' / f'{i+1}'
        zero_itp_file = output_dir / f'{mol_name}_zero.itp'
        zero_top_file = output_dir / f'{mol_name}_zero.top'

        if self.rank == mpi_master():
            output_dir.mkdir(parents=True, exist_ok=True)
            self.copy_file(Path(new_itp_fname), zero_itp_file)
            self.write_top(zero_top_file, zero_itp_file)
            self.set_dihedral_parameters(zero_itp_file, dih, [0.] * 6)
        self.comm.barrier()

        mm_zero_scan = self.perform_mm_scan(zero_top_file,
                                            dih,
                                            geom,
                                            angles,
                                            print_energies=False)

        # fitting Ryckaert-Bellemans function

        if self.rank == mpi_master():
            rel_e_qm = np.array(qm_scan) - min(qm_scan)
            rel_e_qm *= hartree_in_kcalpermol() * 4.184
            rel_e_mm = np.array(mm_zero_scan) - min(mm_zero_scan)

            difference = rel_e_qm - rel_e_mm
            initial_coef = tuple([0.] * 6)

            if kT is not None:
                sigma = 1.0 / np.exp(-rel_e_qm / kT)
                coef, cv = curve_fit(rbpot,
                                     angles,
                                     difference,
                                     initial_coef,
                                     sigma,
                                     absolute_sigma=False)
            else:
                coef, cv = curve_fit(rbpot, angles, difference, initial_coef)

            self.set_dihedral_parameters(new_itp_fname, dih, coef.tolist())
        self.comm.barrier()

        self.ostream.print_info('...done.')
        self.ostream.print_blank()
        self.ostream.print_info(
            f'Generated new topology file: {Path(new_top_fname).name}')
        self.ostream.print_blank()
        self.ostream.flush()

        return Path(new_top_fname)

    def validate_force_field(self, top_file, i):
        """
        Validates force field by RMSD of dihedral potentials.

        :param top_file:
            The topology file.
        :param i:
            The index of the target dihedral.

        :return:
            A dictionary containing the results of validation.
        """

        self.ostream.print_info(f'Validating {top_file.name} ...')
        self.ostream.print_blank()
        self.ostream.flush()

        dih = self.target_dihedrals[i]

        self.ostream.print_info(
            '  Target dihedral angle: {}-{}-{}-{}'.format(*dih))
        self.ostream.print_blank()

        geom = self.scan_geometries[i]
        angles = self.scan_dih_angles[i]
        mm_scan = self.perform_mm_scan(top_file, dih, geom, angles)
        mm_scan = np.array(mm_scan) - min(mm_scan)

        qm_scan = np.array(self.scan_energies[i]) - min(self.scan_energies[i])
        qm_scan *= hartree_in_kcalpermol() * 4.184

        self.ostream.print_blank()
        self.ostream.print_info(
            '      Dihedral      MM energy(rel)      QM energy(rel)       diff')
        self.ostream.print_info(
            '  ---------------------------------------------------------------')
        for angle, e_mm, e_qm in zip(angles, mm_scan, qm_scan):
            self.ostream.print_info(
                f'  {angle:8.1f} deg {e_mm:12.3f} kJ/mol {e_qm:12.3f} kJ/mol ' +
                f'{(e_mm - e_qm):10.3f}')
        self.ostream.print_blank()
        self.ostream.flush()

        return {
            'dihedral_indices': list(dih),
            'dihedral_angles': list(angles),
            'mm_scan_kJpermol': mm_scan.copy(),
            'qm_scan_kJpermol': qm_scan.copy(),
        }

    def perform_mm_scan(self,
                        top_file,
                        dihedral,
                        geometries,
                        angles,
                        print_energies=True):
        """
        Performs MM scan of a specific dihedral.

        :param top_file:
            The topology file.
        :param dihedral:
            The dihedral.
        :param geometries:
            The scanned geometris for this dihedral.
        :param angles:
            The scanned angles for this dihedral.
        """

        # select scan angles and geometries from QM data

        top_fname = str(top_file)

        if self.rank == mpi_master():
            itp_fname = self.get_included_file(top_fname)
        else:
            itp_fname = None
        itp_fname = self.comm.bcast(itp_fname, root=mpi_master())

        top_file = Path(top_fname)

        scan_dir = top_file.parent / (top_file.stem + '_scan')
        if self.rank == mpi_master():
            scan_dir.mkdir(parents=True, exist_ok=True)
        self.comm.barrier()

        if print_energies:
            self.ostream.print_info('      Dihedral           MM energy')
            self.ostream.print_info('  --------------------------------')
            self.ostream.flush()

        energies = []
        for i, (geom, angle) in enumerate(zip(geometries, angles)):

            output_dir = scan_dir / f'{i+1}'
            local_itp_fname = str(output_dir / f'{i+1}.itp')
            local_top_fname = str(output_dir / f'{i+1}.top')

            if self.rank == mpi_master():
                output_dir.mkdir(parents=True, exist_ok=True)
                self.copy_file(Path(itp_fname), Path(local_itp_fname))
                self.write_top(local_top_fname, local_itp_fname)
            self.comm.barrier()

            # energy minimization with dihedral constraint
            constraints = ['$set']
            constraints += ['dihedral {} {} {} {} {}'.format(*dihedral, angle)]
            pot_energy = self.minimize_mm_energy(geom, local_top_fname,
                                                 constraints)
            # convert to kJ/mol
            pot_energy *= 4.184 * hartree_in_kcalpermol()
            energies.append(pot_energy)

            if print_energies:
                self.ostream.print_info(
                    f'  {angle:8.1f} deg {pot_energy:12.3f} kJ/mol')
                self.ostream.flush()

        return energies

    def minimize_mm_energy(self, molecule, top_file, constraints):
        """
        Minimizes MM energy of a topology using OpenMM.

        :param molecule:
            The molecule.
        :param top_file:
            The topology.
        :param constraints:
            The constraints.
        """

        openmm_drv = OpenMMDriver(self.comm)
        openmm_drv.add_topology(top_file, self.gromacs_include_path)

        grad_drv = OpenMMGradientDriver(self.comm)
        grad_drv.ostream.mute()

        opt_drv = OptimizationDriver(grad_drv)
        opt_drv.update_settings({
            'constraints': constraints,
            'filename': str(Path(top_file).parent / Path(top_file).stem),
            'keep_files': self.keep_files,
        })
        final_mol = opt_drv.compute(molecule, openmm_drv)

        openmm_drv.compute(final_mol)

        return openmm_drv.get_energy()

    def set_dihedral_parameters(self, itp_file, dihedral, coefficients):
        """
        Sets dihedral parameters of topology.

        :param itp_file:
            The internal topology file.
        :param dihedral:
            The dihedral.
        :param coefficients:
            The coefficients for Ryckaert-Bellemans funciton.
        """

        itp_fname = str(itp_file)

        atom_names = self.get_atom_names()

        # read itp file and remove existing parameters for chosen dihedral

        dihedral_flag = False
        dihedral_pattern = re.compile(r'\[\s*dihedrals\s*\]')
        saved_itp_lines = []

        added_dih_str = False
        dih_str = '{:6}{:7}{:7}{:7}{:7}'.format(*dihedral, 3)
        for coef in coefficients:
            dih_str += '{:11.5f}'.format(coef)
        dih_str += ' ; {}-{}-{}-{}\n'.format(atom_names[dihedral[0]],
                                             atom_names[dihedral[1]],
                                             atom_names[dihedral[2]],
                                             atom_names[dihedral[3]])

        with open(itp_fname, 'r') as f_itp:
            for line in f_itp:
                title = line.split(';')[0].strip()

                if title.startswith('['):
                    if re.search(dihedral_pattern, title):
                        dihedral_flag = True
                    else:
                        dihedral_flag = False

                if dihedral_flag:
                    content = line.split()
                    try:
                        i, j, k, l, funct = tuple(
                            [int(content[n]) for n in range(5)])
                        condition_1 = (funct in [1, 3, 9])
                        condition_2 = ([i, j, k, l] == dihedral or
                                       [l, k, j, i] == dihedral or
                                       [j, k] == dihedral[1:3] or
                                       [k, j] == dihedral[1:3])
                        if condition_1 and condition_2:
                            if not added_dih_str:
                                saved_itp_lines.append(dih_str)
                                added_dih_str = True
                        else:
                            saved_itp_lines.append(line)
                    except (ValueError, IndexError):
                        saved_itp_lines.append(line)
                else:
                    saved_itp_lines.append(line)

        # update itp file with constraints for chosen dihedral

        with open(itp_fname, 'w') as f_itp:
            for line in saved_itp_lines:
                f_itp.write(line)

    def visualize(self, validation_result):
        """
        Visualizes dihedral potential.

        :param validation_result:
            The dictionary containing the result of validation.
        """

        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError('Unable to import Matplotlib. Please install ' +
                              'Matplotlib via \'conda install matplotlib\'')

        qm_scan_kJpermol = validation_result['qm_scan_kJpermol']
        mm_scan_kJpermol = validation_result['mm_scan_kJpermol']
        dihedral_angles = validation_result['dihedral_angles']
        dihedral_indices = validation_result['dihedral_indices']

        plt.plot(dihedral_angles, qm_scan_kJpermol, '-o', label="QM")
        plt.plot(dihedral_angles, mm_scan_kJpermol, '-o', label="MM")

        plt.grid()
        plt.legend(loc='upper right')
        plt.xlabel('dihedral angle {}-{}-{}-{}'.format(*dihedral_indices))
        plt.ylabel('E in kJ/mol')
        plt.title('dihedral potential of {}'.format(self.molecule_name))
        plt.show()

    def get_included_file(self, top_fname):
        """
        Gets the name of the included itp file.

        :param top_fname:
            The topology file.
        """

        itp_file = None

        with open(top_fname, 'r') as top:
            pattern = re.compile(r'\A#include')
            for line in top:
                if re.search(pattern, line):
                    itp_file = Path(top_fname).parent / line.split('"')[1]

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
