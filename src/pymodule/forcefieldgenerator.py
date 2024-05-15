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
from pathlib import Path, PurePath
import numpy as np
import tempfile
import sys
import re
import xml.etree.ElementTree as ET
from xml.dom import minidom

from .veloxchemlib import mpi_master, bohr_in_angstrom, hartree_in_kcalpermol
from .atomtypeidentifier import AtomTypeIdentifier
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .outputstream import OutputStream
from .respchargesdriver import RespChargesDriver
from .mmdriver import MMDriver
from .mmgradientdriver import MMGradientDriver
from .optimizationdriver import OptimizationDriver
from .inputparser import parse_input, get_random_string_parallel
from .errorhandler import assert_msg_critical, safe_arccos
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
        self.molecule_name = 'vlx_' + get_random_string_parallel(self.comm)
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

        self.workdir = None

        # resp settings
        self.resp_dict = None

        self.keep_files = True
        
        # UFF LJ parameters for all elements.
        # JACS, 1992, 114.25: 10024-10035.
        self.uff_param = {
        'hx': {'sigma': 0.2886, 'epsilon': 0.1839},
        'cx': {'sigma': 0.3851, 'epsilon': 0.4389},
        'nx': {'sigma': 0.3662, 'epsilon': 0.3219},
        'ox': {'sigma': 0.3405, 'epsilon': 0.4013},
        'px': {'sigma': 0.415, 'epsilon': 0.1338},
        'sx': {'sigma': 0.403, 'epsilon': 0.1437},
        'He': {'sigma': 0.2104, 'epsilon': 0.2343},
        'Li': {'sigma': 0.2184, 'epsilon': 0.1046},
        'Be': {'sigma': 0.2446, 'epsilon': 0.3556},
        'B': {'sigma': 0.3638, 'epsilon': 0.7531},
        'F': {'sigma': 0.2997, 'epsilon': 0.2092},
        'Na': {'sigma': 0.2658, 'epsilon': 0.1255},
        'Mg': {'sigma': 0.2691, 'epsilon': 0.4644},
        'Al': {'sigma': 0.4008, 'epsilon': 2.1128},
        'Si': {'sigma': 0.3826, 'epsilon': 1.6819},
        'Cl': {'sigma': 0.3516, 'epsilon': 0.9497},
        'Ar': {'sigma': 0.3446, 'epsilon': 0.774},
        'K': {'sigma': 0.3396, 'epsilon': 0.1464},
        'Ca': {'sigma': 0.3028, 'epsilon': 0.9957},
        'Sc': {'sigma': 0.2936, 'epsilon': 0.0795},
        'Ti': {'sigma': 0.2829, 'epsilon': 0.0711},
        'V': {'sigma': 0.2801, 'epsilon': 0.0669},
        'Cr': {'sigma': 0.2693, 'epsilon': 0.0628},
        'Mn': {'sigma': 0.2638, 'epsilon': 0.0544},
        'Fe': {'sigma': 0.2594, 'epsilon': 0.0544},
        'Co': {'sigma': 0.2559, 'epsilon': 0.0586},
        'Ni': {'sigma': 0.2525, 'epsilon': 0.0628},
        'Cu': {'sigma': 0.3114, 'epsilon': 0.0209},
        'Zn': {'sigma': 0.2462, 'epsilon': 0.5188},
        'Ga': {'sigma': 0.3905, 'epsilon': 1.7363},
        'Ge': {'sigma': 0.3813, 'epsilon': 1.5856},
        'As': {'sigma': 0.3769, 'epsilon': 1.2928},
        'Se': {'sigma': 0.3746, 'epsilon': 1.2175},
        'Br': {'sigma': 0.3732, 'epsilon': 1.0501},
        'Kr': {'sigma': 0.3689, 'epsilon': 0.9204},
        'Rb': {'sigma': 0.3665, 'epsilon': 0.1674},
        'Sr': {'sigma': 0.3244, 'epsilon': 0.9832},
        'Y': {'sigma': 0.298, 'epsilon': 0.3012},
        'Zr': {'sigma': 0.2783, 'epsilon': 0.2887},
        'Nb': {'sigma': 0.282, 'epsilon': 0.2468},
        'Mo': {'sigma': 0.2719, 'epsilon': 0.2343},
        'Tc': {'sigma': 0.2671, 'epsilon': 0.2008},
        'Ru': {'sigma': 0.264, 'epsilon': 0.2343},
        'Rh': {'sigma': 0.2609, 'epsilon': 0.2217},
        'Pd': {'sigma': 0.2583, 'epsilon': 0.2008},
        'Ag': {'sigma': 0.2805, 'epsilon': 0.1506},
        'Cd': {'sigma': 0.2537, 'epsilon': 0.9539},
        'In': {'sigma': 0.3976, 'epsilon': 2.5061},
        'Sn': {'sigma': 0.3913, 'epsilon': 2.3722},
        'Sb': {'sigma': 0.3938, 'epsilon': 1.8785},
        'Te': {'sigma': 0.3982, 'epsilon': 1.6651},
        'I': {'sigma': 0.4009, 'epsilon': 1.4183},
        'Xe': {'sigma': 0.3924, 'epsilon': 1.389},
        'Cs': {'sigma': 0.4024, 'epsilon': 0.1883},
        'Ba': {'sigma': 0.3299, 'epsilon': 1.5229},
        'La': {'sigma': 0.3138, 'epsilon': 0.0711},
        'Ce': {'sigma': 0.3168, 'epsilon': 0.0544},
        'Pr': {'sigma': 0.3213, 'epsilon': 0.0418},
        'Nd': {'sigma': 0.3185, 'epsilon': 0.0418},
        'Pm': {'sigma': 0.316, 'epsilon': 0.0377},
        'Sm': {'sigma': 0.3136, 'epsilon': 0.0335},
        'Eu': {'sigma': 0.3112, 'epsilon': 0.0335},
        'Gd': {'sigma': 0.3001, 'epsilon': 0.0377},
        'Tb': {'sigma': 0.3074, 'epsilon': 0.0293},
        'Dy': {'sigma': 0.3054, 'epsilon': 0.0293},
        'Ho': {'sigma': 0.3037, 'epsilon': 0.0293},
        'Er': {'sigma': 0.3021, 'epsilon': 0.0293},
        'Tm': {'sigma': 0.3006, 'epsilon': 0.0251},
        'Yb': {'sigma': 0.2989, 'epsilon': 0.9539},
        'Lu': {'sigma': 0.3243, 'epsilon': 0.1715},
        'Hf': {'sigma': 0.2798, 'epsilon': 0.3012},
        'Ta': {'sigma': 0.2824, 'epsilon': 0.3389},
        'W': {'sigma': 0.2734, 'epsilon': 0.2803},
        'Re': {'sigma': 0.2632, 'epsilon': 0.2761},
        'Os': {'sigma': 0.278, 'epsilon': 0.1548},
        'Ir': {'sigma': 0.253, 'epsilon': 0.3054},
        'Pt': {'sigma': 0.2454, 'epsilon': 0.3347},
        'Au': {'sigma': 0.2934, 'epsilon': 0.1632},
        'Hg': {'sigma': 0.241, 'epsilon': 1.6108},
        'Tl': {'sigma': 0.3873, 'epsilon': 2.845},
        'Pb': {'sigma': 0.3828, 'epsilon': 2.7738},
        'Bi': {'sigma': 0.3893, 'epsilon': 2.1672},
        'Po': {'sigma': 0.4195, 'epsilon': 1.3597},
        'At': {'sigma': 0.4232, 'epsilon': 1.1882},
        'Rn': {'sigma': 0.4245, 'epsilon': 1.0376},
        'Fr': {'sigma': 0.4365, 'epsilon': 0.2092},
        'Ra': {'sigma': 0.3276, 'epsilon': 1.6902},
        'Ac': {'sigma': 0.3099, 'epsilon': 0.1381},
        'Th': {'sigma': 0.3025, 'epsilon': 0.1088},
        'Pa': {'sigma': 0.305, 'epsilon': 0.092},
        'U': {'sigma': 0.3025, 'epsilon': 0.092},
        'Np': {'sigma': 0.305, 'epsilon': 0.0795},
        'Pu': {'sigma': 0.305, 'epsilon': 0.0669},
        'Am': {'sigma': 0.3012, 'epsilon': 0.0586},
        'Cm': {'sigma': 0.2963, 'epsilon': 0.0544},
        'Bk': {'sigma': 0.2975, 'epsilon': 0.0544},
        'Cf': {'sigma': 0.2952, 'epsilon': 0.0544},
        'Es': {'sigma': 0.2939, 'epsilon': 0.0502},
        'Fm': {'sigma': 0.2927, 'epsilon': 0.0502},
        'Md': {'sigma': 0.2917, 'epsilon': 0.046},
        'No': {'sigma': 0.2894, 'epsilon': 0.046},
        'Lw': {'sigma': 0.2883, 'epsilon': 0.046},
        }
        
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

        use_temp_dir = (self.workdir is None)

        if use_temp_dir:
            try:
                temp_dir = tempfile.TemporaryDirectory(
                    ignore_cleanup_errors=True)
            except TypeError:
                temp_dir = tempfile.TemporaryDirectory()
            self.workdir = Path(temp_dir.name)

        if self.original_top_file is None:
            original_itp_file = self.workdir / (mol_name +
                                                f'_{self.ffversion:02d}.itp')
            original_top_file = original_itp_file.with_suffix('.top')

            if self.rank == mpi_master():
                self.write_gromacs_files(original_itp_file)
            self.comm.barrier()

        else:
            # TODO: read topology file into forcefieldgenerator
            original_top_file = Path(self.original_top_file)

        self.ostream.print_info(
            f'Original topology file: {original_top_file.name}')
        self.ostream.print_blank()
        self.ostream.flush()

        # validate original force field

        for i, dih in enumerate(self.target_dihedrals):
            self.validate_force_field(i)

        # fit dihedral potentials

        n_rounds = self.n_rounds if len(self.scan_dih_angles) > 1 else 1
        for i_round in range(n_rounds):
            for i, dih in enumerate(self.target_dihedrals):
                self.dihedral_correction(i)

        # validate final force field

        for i, dih in enumerate(self.target_dihedrals):
            self.validate_force_field(i)

        # save output files

        if self.rank == mpi_master() and self.keep_files:
            out_dir = Path(self.molecule_name + '_files')
            for ffver in range(self.ffversion + 1):
                for ftype in ['itp', 'top']:
                    fname = mol_name + f'_{ffver:02d}.{ftype}'
                    self.copy_file(self.workdir / fname, out_dir / fname)
                    valstr = f'Saving file: {str(out_dir / fname)}'
                    self.ostream.print_info(valstr)

            self.ostream.print_blank()
            self.ostream.flush()

        if use_temp_dir:
            try:
                temp_dir.cleanup()
            except (NotADirectoryError, PermissionError):
                pass
            self.workdir = None

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
                    i, j, k, l = [
                        int(i) - 1
                        for i in line.split('Dihedral')[1].split()[0].split('-')
                    ]
                    dih_inds = [i, j, k, l] if i < l else [l, k, j, i]
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

    def create_topology(self, molecule, basis=None, scf_result=None, no_resp=False):
        """
        Analyzes the topology of the molecule and create dictionaries
        for the atoms, bonds, angles, dihedrals, impropers and pairs.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_result:
            A converged SCF result.
        :param no_resp:
            If RESP charges should not be computed.
            Partial charges will be set to zero.
        """

        # Read the force field data lines.

        if self.force_field_data is None:
            ff_data_lines = self.get_gaff_data_lines()
        else:
            with open(self.force_field_data, 'r') as ff_data:
                ff_data_lines = ff_data.readlines()

        # check GAFF version
        gaff_version = None
        if ff_data_lines[0].startswith('AMBER General Force Field'):
            ff_data_version = ff_data_lines[0].split('Version')[1].split()[0]
            ff_data_version = ff_data_version.replace(',', '')
            if '.' in ff_data_version:
                version_major = ff_data_version.split('.')[0]
                version_minor = ff_data_version.split('.')[1]
                if version_major.isdigit() and version_minor.isdigit():
                    gaff_version = f'{version_major}.{version_minor}'

        # Molecular information

        self.molecule = molecule

        coords = self.molecule.get_coordinates_in_angstrom()
        n_atoms = self.molecule.number_of_atoms()

        atomtypeidentifier = AtomTypeIdentifier(self.comm)
        atomtypeidentifier.ostream.mute()
        # set GAFF version
        atomtypeidentifier.gaff_version = gaff_version

        self.atom_types = atomtypeidentifier.generate_gaff_atomtypes(
            self.molecule)
        atomtypeidentifier.identify_equivalences()

        self.connectivity_matrix = np.copy(
            atomtypeidentifier.connectivity_matrix)
        
        if no_resp:
            # Skip calculations
            self.partial_charges = np.zeros(self.molecule.number_of_atoms())
            msg = 'RESP calculation disabled: All partial charges are set to zero.'
            self.ostream.print_info(msg)
            self.ostream.flush()

        # Default behavior: compute RESP charges
        if self.partial_charges is None:
            if scf_result is None:
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

                resp_drv = RespChargesDriver(self.comm, self.ostream)
                resp_drv.filename = self.molecule_name
                if self.resp_dict is not None:
                    resp_drv.update_settings(self.resp_dict)
                if resp_drv.equal_charges is None:
                    resp_drv.equal_charges = atomtypeidentifier.equivalent_charges

                self.partial_charges = resp_drv.compute(self.molecule, basis,
                                                        'resp')
                self.partial_charges = self.comm.bcast(self.partial_charges,
                                                    root=mpi_master())

            # Else use the provided SCF result   
            else:
                if basis is None:
                    error_msg = 'Basis is required for RESP charges.'
                    assert_msg_critical(False, error_msg)

                resp_drv = RespChargesDriver(self.comm, self.ostream)
                resp_drv.filename = self.molecule_name
                msg = 'Using provided SCF result for RESP charges'
                self.ostream.print_info(msg)
                self.ostream.flush()
                if self.resp_dict is not None:
                    resp_drv.update_settings(self.resp_dict)
                if resp_drv.equal_charges is None:
                    resp_drv.equal_charges = atomtypeidentifier.equivalent_charges

                self.partial_charges = resp_drv.compute(self.molecule, basis,
                                                        scf_result, 'resp')
                self.partial_charges = self.comm.bcast(self.partial_charges,
                                                    root=mpi_master())
                
        # self.partial charges should be summing to a natural number
        # it should be rounded to double precision

        self.partial_charges = np.round(self.partial_charges, 16)


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
            
            # Auxilary variable for finding parameters in UFF
            element = ''.join([i for i in at if not i.isdigit()])

            for line in ff_data_lines:
                if line.startswith(f'  {at}     '):
                    atom_ff = line[5:].strip().split()
                    sigma = float(atom_ff[0]) * 2**(-1 / 6) * 2 / 10
                    epsilon = float(atom_ff[1]) * 4.184
                    comment = 'GAFF'
                    atom_type_found = True
                    break

            if not atom_type_found:
                if at == 'ow':
                    sigma, epsilon, comment = 3.15061e-01, 6.36386e-01, 'OW'
                elif at == 'hw':
                    sigma, epsilon, comment = 0.0, 0.0, 'HW'
                # Case for atoms in UFF but not in GAFF
                elif element in self.uff_param.keys():
                    warnmsg = f'ForceFieldGenerator: atom type {at} is not in GAFF. Sigma and Epsilon from UFF.'
                    self.ostream.print_warning(warnmsg)
                    sigma = self.uff_param[element]['sigma']
                    epsilon = self.uff_param[element]['epsilon']
                    comment = 'UFF'
                else:
                    warnmsg = f'ForceFieldGenerator: atom type {at} is ill defined. Default values are used.'
                    self.ostream.print_warning(warnmsg)
                    sigma = 0.35
                    epsilon = 0.1
                    comment = 'Default'

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
                # Default value for bonds
                r, k_r, comment = r_eq, 2.5e+5, 'Guessed'

            if self.eq_param:
                if abs(r - r_eq) > self.r_thresh:
                    msg = f'Updated bond length {i+1}-{j+1} '
                    msg += f'({at_1}-{at_2}) to {r_eq:.3f} nm'
                    self.ostream.print_info(msg)
                    self.ostream.flush()
                r = r_eq

            self.bonds[(i, j)] = {
                'type': 'harmonic',
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
            theta_eq = safe_arccos(
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
                # Default value for angles
                theta, k_theta, comment = theta_eq, 1000, 'Guessed'

            if self.eq_param:
                if abs(theta - theta_eq) > self.theta_thresh:
                    msg = f'Updated bond angle {i+1}-{j+1}-{k+1} '
                    msg += f'({at_1}-{at_2}-{at_3}) to {theta_eq:.3f} deg'
                    self.ostream.print_info(msg)
                    self.ostream.flush()
                theta = theta_eq

            self.angles[(i, j, k)] = {
                'type': 'harmonic',
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

            # special treatment for rotatable bonds, e.g. cc-cc, cd-cd, cc-na
            # and cd-na bonds between non-pure aromatic rings
            special_comment = ''

            if [at_2, at_3] in [['cc', 'cc'], ['cd', 'cd']]:
                if not atomtypeidentifier.get_common_cycles(
                        j, k, 'non_pure_aromatic'):
                    patterns = [
                        re.compile(r'\A' + 'X -cp-cp-X  '),
                    ]
                    special_comment = ('(Guessed for rotatable ' +
                                       f'{at_1}-{at_2}-{at_3}-{at_4})')

            elif [at_2, at_3] in [['cc', 'na'], ['na', 'cc'], ['cd', 'na'],
                                  ['na', 'cd']]:
                if not atomtypeidentifier.get_common_cycles(
                        j, k, 'non_pure_aromatic'):
                    patterns = [
                        re.compile(r'\A' + 'X -ca-na-X  '),
                        re.compile(r'\A' + 'X -na-ca-X  '),
                    ]
                    special_comment = ('(Guessed for rotatable ' +
                                       f'{at_1}-{at_2}-{at_3}-{at_4})')

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
                            dihedral_matches.append(
                                m.group(0) + special_comment)
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
                # guesses for proper dihedrals
                patterns = self.get_dihedral_guess_patterns(at_2, at_3)

                dihedral_ff_lines = []
                dihedral_matches = []
                for line in ff_data_lines:
                    for p in patterns:
                        m = re.search(p, line)
                        if m is not None:
                            dihedral_ff = line[11:60].strip().split()
                            if len(dihedral_ff) == 4:
                                dihedral_ff_lines.append(line)
                                dihedral_matches.append(
                                    m.group(0) + '(Guessed for ' +
                                    f'{at_1}-{at_2}-{at_3}-{at_4})')
                                dihedral_found = True
                                break

            if not dihedral_found:
                warnmsg = f'ForceFieldGenerator: dihedral {at_1}-{at_2}-{at_3}-{at_4}'
                warnmsg += ' is not available.'
                self.ostream.print_warning(warnmsg)
                # Default value for dihedrals
                self.dihedrals[(i, j, k, l)] = {
                    'type': 'Fourier',
                    'multiple': False,
                    'barrier': 0.0,
                    'phase': 0.0,
                    'periodicity': 1,
                    'comment': f'Unknown {at_1}-{at_2}-{at_3}-{at_4}'
                }

            dihedral_barriers = []
            dihedral_phases = []
            dihedral_periodicities = []
            dihedral_comments = []

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

                dihedral_barriers.append(barrier)
                dihedral_phases.append(phase)
                dihedral_periodicities.append(periodicity)
                dihedral_comments.append(comment)

                if periodicity > 0:
                    break

            if len(dihedral_barriers) == 1:

                self.dihedrals[(i, j, k, l)] = {
                    'type': 'Fourier',
                    'multiple': False,
                    'barrier': dihedral_barriers[0],
                    'phase': dihedral_phases[0],
                    'periodicity': dihedral_periodicities[0],
                    'comment': dihedral_comments[0]
                }

            elif len(dihedral_barriers) > 1:

                self.dihedrals[(i, j, k, l)] = {
                    'type': 'Fourier',
                    'multiple': True,
                    'barrier': dihedral_barriers,
                    'phase': dihedral_phases,
                    'periodicity': dihedral_periodicities,
                    'comment': dihedral_comments,
                }

        # convert multi-termed dihedral to RB type
        for (i, j, k, l), dih in self.dihedrals.items():

            if dih['multiple']:

                valid_phases = True
                for phase in dih['phase']:
                    if not (abs(phase) < 1.0e-6 or abs(phase - 180.0) < 1.0e-6):
                        valid_phases = False
                        break

                valid_periodicity = True
                for periodicity in dih['periodicity']:
                    if abs(periodicity) not in [1, 2, 3, 4]:
                        valid_periodicity = False
                        break

                if not (valid_phases and valid_periodicity):
                    continue

                F_coefs = {x: 0.0 for x in [1, 2, 3, 4]}
                E_shift = 0.0

                for barrier, phase, periodicity in zip(dih['barrier'],
                                                       dih['phase'],
                                                       dih['periodicity']):
                    if abs(phase) < 1.0e-6:
                        # phase == 0 degree
                        F_coefs[abs(periodicity)] += barrier
                    else:
                        # phase == 180 degree
                        F_coefs[abs(periodicity)] -= barrier
                        E_shift += 2.0 * barrier

                C_coefs = [0.0 for x in range(6)]

                # JPCA 2021, 125, 2673-2681
                # Note that we also take into account the phases in Fourier series
                C_coefs[0] = (F_coefs[1] + F_coefs[3] + 2.0 * F_coefs[4] +
                              E_shift)
                C_coefs[1] = -1.0 * F_coefs[1] + 3.0 * F_coefs[3]
                C_coefs[2] = 2.0 * F_coefs[2] - 8.0 * F_coefs[4]
                C_coefs[3] = -4.0 * F_coefs[3]
                C_coefs[4] = 8.0 * F_coefs[4]
                C_coefs[5] = 0.0

                self.dihedrals[(i, j, k, l)] = {
                    'type': 'RB',
                    'RB_coefficients': C_coefs,
                    'comment': dih['comment'][0] + ' RB',
                }

        # Impropers

        self.impropers = {}

        sp2_atom_types = [
            'c ', 'cs', 'c2', 'ca', 'cp', 'cq', 'cc', 'cd', 'ce', 'cf', 'cu',
            'cv', 'cz', 'n ', 'ns', 'nt', 'n2', 'na', 'nb', 'nc', 'nd', 'ne',
            'nf', 'pb', 'pc', 'pd', 'pe', 'pf'
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
                    # Default values for impropers
                    barrier, phase, periodicity = 1.1 * 4.184, 180.0, 2
                    comment = 'Guessed'

                assert_msg_critical(
                    phase == 180.0,
                    'ForceFieldGenerator: invalid improper dihedral phase')

                assert_msg_critical(
                    periodicity == 2,
                    'ForceFieldGenerator: invalid improper dihedral periodicity'
                )

                self.impropers[(i, j, k, l)] = {
                    'type': 'Fourier',
                    'barrier': barrier,
                    'phase': phase,
                    'periodicity': periodicity,
                    'comment': comment
                }

    def add_bond(self, bond, force_constant=250000.00, equilibrium=None):
        """
        Adds a bond to the topology.

        :param bond:
            The bond to be added. As a tuple of atom indices.
        :param force_constant:
            The force constant of the bond. Default is 250000.00 kJ/mol/nm^2.
        :param equilibrium:
            The equilibrium distance of the bond. If none it will be calculated.
        """

        i, j = bond

        # Convert to zero-based indices
        i = i - 1
        j = j - 1

        if equilibrium is None:
            coords = self.molecule.get_coordinates_in_angstrom()
            equilibrium = np.linalg.norm(coords[i] - coords[j]) * 0.1

        self.bonds[bond] = {
            'type': 'harmonic',
            'force_constant': force_constant,
            'equilibrium': equilibrium,
            'comment': 'User-defined'
        }

    def add_angle(self, angle, force_constant=1000.00, equilibrium=None):
        """
        Adds an angle to the topology.

        :param angle:
            The angle to be added. As a tuple of atom indices.
        :param force_constant:
            The force constant of the angle. Default is 1000.00 kJ/mol/rad^2.
        :param equilibrium:
            The equilibrium angle of the angle. If none it will be calculated.
        """

        i, j, k = angle

        # Convert to zero-based indices
        i = i - 1
        j = j - 1
        k = k - 1

        if equilibrium is None:
            coords = self.molecule.get_coordinates_in_angstrom()
            a = coords[i] - coords[j]
            b = coords[k] - coords[j]
            equilibrium = safe_arccos(
                np.dot(a, b) / np.linalg.norm(a) /
                np.linalg.norm(b)) * 180 / np.pi

        self.angles[angle] = {
            'type': 'harmonic',
            'force_constant': force_constant,
            'equilibrium': equilibrium,
            'comment': 'User-defined'
        }

    def get_dihedral_guess_patterns(self, at_2, at_3):
        """
        Gets guesses for dihedral parameters.

        :param at_2:
            The index of the second atom in the dihedral.
        :param at_3:
            The index of the third atom in the dihedral.

        :return:
            A list of patterns.
        """

        atomtype_pairs_mapping = {
            ('cp', 'cq'): ('ca', 'ca'),
            # ---
            ('nb', 'nb'): ('ca', 'nb'),
            ('nb', 'cp'): ('ca', 'cp'),
            # ---
            ('cc', 'no'): ('ca', 'no'),
            ('cd', 'no'): ('ca', 'no'),
            # ---
            ('ce', 'c3'): ('c2', 'c3'),
            ('ce', 'c5'): ('c2', 'c3'),
            ('ce', 'c6'): ('c2', 'c3'),
            ('cf', 'c3'): ('c2', 'c3'),
            ('cf', 'c5'): ('c2', 'c3'),
            ('cf', 'c6'): ('c2', 'c3'),
            # ---
            ('ce', 'cc'): ('ce', 'ce'),
            ('cf', 'cd'): ('cf', 'cf'),
            ('ce', 'cd'): ('ce', 'cf'),
            ('cc', 'cf'): ('ce', 'cf'),
            # ---
            ('ne', 'cc'): ('ne', 'ce'),
            ('nf', 'cd'): ('nf', 'cf'),
            # ---
            ('cc', 'n2'): ('cc', 'nc'),
            ('cd', 'n2'): ('cd', 'nd'),
            # ---
            ('ce', 'nf'): ('c2', 'n2'),
            ('cf', 'ne'): ('c2', 'n2'),
            ('ce', 'n2'): ('c2', 'n2'),
            ('cf', 'n2'): ('c2', 'n2'),
            # ---
            ('ce', 'nu'): ('c2', 'nh'),
            ('ce', 'nv'): ('c2', 'nh'),
            ('cf', 'nu'): ('c2', 'nh'),
            ('cf', 'nv'): ('c2', 'nh'),
        }

        for at_pair, new_at_pair in atomtype_pairs_mapping.items():
            if (at_2, at_3) == at_pair or (at_3, at_2) == at_pair:
                new_at_2, new_at_3 = new_at_pair
                if new_at_2 == new_at_3:
                    return [
                        re.compile(r'\A' + f'X -{new_at_2}-{new_at_3}-X  '),
                    ]
                else:
                    return [
                        re.compile(r'\A' + f'X -{new_at_2}-{new_at_3}-X  '),
                        re.compile(r'\A' + f'X -{new_at_3}-{new_at_2}-X  '),
                    ]

        for at_val in ['ca', 'os', 'ss', 'oh', 'sh']:
            condition_1 = (at_2 == at_val and at_3 in ['cc', 'cd', 'ce', 'cf'])
            condition_2 = (at_2 in ['cc', 'cd', 'ce', 'cf'] and at_3 == at_val)
            if condition_1 or condition_2:
                return [
                    re.compile(r'\A' + f'X -c2-{at_val}-X  '),
                    re.compile(r'\A' + f'X -{at_val}-c2-X  '),
                ]

        atom_types_mapping = {
            ('cu', 'cv'): 'c2',
            ('cx', 'cy', 'c5', 'c6'): 'c3',
            ('nt', 'ns'): 'n ',
            ('nu', 'nv'): 'nh',
            ('n7', 'n8', 'n5', 'n6'): 'n3',
            ('cs',): 'c ',
        }
        new_at_2, new_at_3 = at_2, at_3
        for key, val in atom_types_mapping.items():
            if at_2 in key:
                new_at_2 = val
            if at_3 in key:
                new_at_3 = val
        if new_at_2 != at_2 or new_at_3 != at_3:
            if new_at_2 == new_at_3:
                return [
                    re.compile(r'\A' + f'X -{new_at_2}-{new_at_3}-X  '),
                ]
            else:
                return [
                    re.compile(r'\A' + f'X -{new_at_2}-{new_at_3}-X  '),
                    re.compile(r'\A' + f'X -{new_at_3}-{new_at_2}-X  '),
                ]

        return []
    
    def add_bond(self, bond, force_constant=250000.00, equilibrium=None):
        """
        Adds a bond to the topology.

        :param bond:
            The bond to be added. As a tuple of atom indices.
        :param force_constant:
            The force constant of the bond. Default is 250000.00 kJ/mol/nm^2.
        :param equilibrium:
            The equilibrium distance of the bond. If none it will be calculated.
        """

        i, j = bond

        # Convert to zero-based indices
        i = i - 1
        j = j - 1

        if equilibrium is None:
            coords = self.molecule.get_coordinates_in_angstrom()
            equilibrium = np.linalg.norm(coords[i] - coords[j]) * 0.1

        self.bonds[bond] = {
            'type': 'harmonic',
            'force_constant': force_constant,
            'equilibrium': equilibrium,
            'comment': 'User-defined'
        }

    def add_angle(self, angle, force_constant=1000.00, equilibrium=None):
        """
        Adds an angle to the topology.

        :param angle:
            The angle to be added. As a tuple of atom indices.
        :param force_constant:
            The force constant of the angle. Default is 1000.00 kJ/mol/rad^2.
        :param equilibrium:
            The equilibrium angle of the angle. If none it will be calculated.
        """

        i, j, k = angle

        # Convert to zero-based indices
        i = i - 1
        j = j - 1
        k = k - 1

        if equilibrium is None:
            coords = self.molecule.get_coordinates_in_angstrom()
            a = coords[i] - coords[j]
            b = coords[k] - coords[j]
            equilibrium = safe_arccos(
                np.dot(a, b) / np.linalg.norm(a) /
                np.linalg.norm(b)) * 180 / np.pi

        self.angles[angle] = {
            'type': 'harmonic',
            'force_constant': force_constant,
            'equilibrium': equilibrium,
            'comment': 'User-defined'
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

            xtb_drv = XtbDriver(self.comm, self.ostream)
            xtb_grad_drv = XtbGradientDriver(xtb_drv)
            xtb_opt_drv = OptimizationDriver(xtb_grad_drv)
            xtb_opt_drv.filename = self.molecule_name
            xtb_opt_results = xtb_opt_drv.compute(self.molecule)
            self.molecule = Molecule.read_xyz_string(
                xtb_opt_results['final_geometry'])

            # XTB Hessian
            self.ostream.print_info('Computing Hessian using XTB...')
            self.ostream.flush()

            xtb_hessian_drv = XtbHessianDriver(xtb_drv)
            xtb_hessian_drv.compute(self.molecule)
            hessian = np.copy(xtb_hessian_drv.hessian)

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
                    if self.bonds[(i, j)]['comment'].capitalize() != 'Guessed':
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
                            != 'Guessed'):
                        continue
                elif (i, j, k) not in reparameterize_keys:
                    continue

            a = coords_in_au[i] - coords_in_au[j]
            b = coords_in_au[k] - coords_in_au[j]
            new_equilibrium = safe_arccos(
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

    def write_top(self,
                  top_file,
                  itp_file,
                  mol_name=None,
                  amber_ff=None,
                  water_model=None):
        """
        Writes a topology file.

        :param top_file:
            The topology file.
        :param itp_file:
            The included itp file.
        :param mol_name:
            The name of the molecule.
        :param amber_ff:
            The name of the Amber force field.
        :param water_model:
            The name of the water model.
        """

        top_fname = str(top_file)
        itp_fname = str(itp_file)

        if mol_name is None:
            mol_name = Path(self.molecule_name).stem

        with open(top_fname, 'w') as f_top:

            # header

            f_top.write('; Generated by VeloxChem\n')

            # defaults

            if amber_ff is not None:
                assert_msg_critical(
                    amber_ff.startswith('amber'),
                    'ForceFieldGenerator.write_top: Invalid amber force field name'
                )
                ff_include = str(PurePath(f'{amber_ff}.ff') / 'forcefield.itp')
                f_top.write(f'\n#include "{ff_include}"\n')
            else:
                f_top.write('\n[ defaults ]\n')
                cur_str = '; nbfunc        comb-rule       gen-pairs'
                cur_str += '        fudgeLJ   fudgeQQ\n'
                f_top.write(cur_str)
                gen_pairs = 'yes' if self.gen_pairs else 'no'
                f_top.write('{}{:16}{:>18}{:21.6f}{:10.6f}\n'.format(
                    self.nbfunc, self.comb_rule, gen_pairs, self.fudgeLJ,
                    self.fudgeQQ))

            # include itp

            f_top.write('\n#include "' + Path(itp_fname).name + '"\n')

            if water_model is not None:
                # very rudimentary check for water model names
                assert_msg_critical(
                    water_model.startswith('tip') or
                    water_model.startswith('spc'),
                    'ForceFieldGenerator.write_top: Invalid water model name')
                assert_msg_critical(
                    amber_ff is not None, 'ForceFieldGenerator.write_top: ' +
                    'amber_ff is required for water_model')
                water_include = str(
                    PurePath(f'{amber_ff}.ff') / f'{water_model}.itp')
                f_top.write(f'\n#include "{water_include}"\n')

            # system

            f_top.write('\n[ system ]\n')
            f_top.write('{}\n'.format(mol_name))

            # molecules

            f_top.write('\n[ molecules ]\n')
            f_top.write('; Compound        nmols\n')
            f_top.write('{:<10}{:9}\n'.format(mol_name, 1))

    def write_itp(self, itp_file, mol_name=None):
        """
        Writes an ITP file with the original parameters.

        :param itp_file:
            The ITP file path.
        """

        itp_filename = str(itp_file)
        if mol_name is None:
            mol_name = Path(self.molecule_name).stem
        res_name = 'MOL'

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
            f_itp.write(f'{mol_name:<10}  {self.nrexcl:10}\n')

            # Atoms
            f_itp.write('\n[ atoms ]\n')
            line_str = ';   nr  type  resi  res  atom  cgnr'
            line_str += '     charge       mass\n'
            f_itp.write(line_str)

            total_charge = 0.0
            for i, atom in self.atoms.items():
                total_charge += atom['charge']
                line_str = '{:6}{:1}{:>5}{:6}{:>6}{:>6}'.format(
                    i + 1,'', atom['type'], 1, res_name, atom['name'])
                line_str += '{:5}{:13.6f}{:13.5f}'.format(
                    i + 1, atom['charge'], atom['mass'])
                line_str += ' ; qtot{:7.3f}  equiv. {}\n'.format(
                    total_charge, atom['equivalent_atom'])
                f_itp.write(line_str)

            # Bonds
            if self.bonds:
                f_itp.write('\n[ bonds ]\n')
                f_itp.write(';   ai     aj    funct       r           k_r\n')

            for (i, j), bond in self.bonds.items():
                line_str = '{:6}{:7}{:7}{:14.4e}{:14.4e} ; {}\n'.format(
                    i + 1, j + 1, 1, bond['equilibrium'],
                    bond['force_constant'], bond['comment'])
                f_itp.write(line_str)

            # Pairs
            if self.pairs:
                f_itp.write('\n[ pairs ]\n')
                f_itp.write(';   ai     aj    funct\n')

            for i, j in self.pairs:
                f_itp.write('{:6}{:7}{:7}\n'.format(i + 1, j + 1, 1))

            # Angles
            if self.angles:
                f_itp.write('\n[ angles ]\n')
                f_itp.write(
                    ';   ai     aj     ak    funct     theta       k_theta\n')

            for (i, j, k), angle in self.angles.items():
                line_str = '{:6}{:7}{:7}{:7}{:14.4e}{:14.4e} ; {}\n'.format(
                    i + 1, j + 1, k + 1, 1, angle['equilibrium'],
                    angle['force_constant'], angle['comment'])
                f_itp.write(line_str)

            # Proper dihedrals
            if self.dihedrals:
                f_itp.write('\n[ dihedrals ]\n')
                f_itp.write('; propers\n')

            dih_RB_lines = []
            dih_fourier_lines = []

            for (i, j, k, l), dih in self.dihedrals.items():
                if dih['type'] == 'RB':
                    line_str = '{:6}{:7}{:7}{:7}{:7}'.format(
                        i + 1, j + 1, k + 1, l + 1, 3)
                    for coef in dih['RB_coefficients']:
                        line_str += '{:11.5f}'.format(coef)
                    line_str += ' ; {}\n'.format(dih['comment'])
                    dih_RB_lines.append(line_str)

                elif dih['type'] == 'Fourier':
                    if dih['multiple']:
                        for barrier, phase, periodicity, comment in zip(
                                dih['barrier'], dih['phase'],
                                dih['periodicity'], dih['comment']):
                            line_str = '{:6}{:7}{:7}{:7}'.format(
                                i + 1, j + 1, k + 1, l + 1)
                            line_str += '{:7}{:11.2f}{:11.5f}{:4} ; {}\n'.format(
                                9, phase, barrier, abs(periodicity), comment)
                            dih_fourier_lines.append(line_str)
                    else:
                        line_str = '{:6}{:7}{:7}{:7}'.format(
                            i + 1, j + 1, k + 1, l + 1)
                        line_str += '{:7}{:11.2f}{:11.5f}{:4} ; {}\n'.format(
                            9, dih['phase'], dih['barrier'],
                            abs(dih['periodicity']), dih['comment'])
                        dih_fourier_lines.append(line_str)

                else:
                    errmsg = 'ForceFieldGenerator.create_topology:'
                    errmsg += ' Invalid dihedral type ' + dih['type']
                    assert_msg_critical(False, errmsg)

            if dih_RB_lines:
                line_str = ';   ai     aj     ak     al    funct'
                line_str += '     C0         C1         C2         C3'
                line_str += '         C4         C5\n'
                f_itp.write(line_str)

            for line_str in dih_RB_lines:
                f_itp.write(line_str)

            if dih_fourier_lines:
                f_itp.write(
                    ';   ai     aj     ak     al    funct    phase     k_d      n\n'
                )

            for line_str in dih_fourier_lines:
                f_itp.write(line_str)

            # Improper dihedrals
            if self.impropers:
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
    
    def generate_residue_xml(self, xml_file, mol_name='MOL'):
        """
        Generates an XML force field file for a single residue.
        
        :param mol_name:
            The name of the molecule.
        :param filename:
            The name of the XML file.
        """

        filename = str(xml_file)

        atoms = self.atoms
        bonds = self.bonds

        # Create the root element of the XML file
        ForceField = ET.Element("ForceField")
        
        # AtomTypes section
        AtomTypes = ET.SubElement(ForceField, "AtomTypes")

        for i, atom in self.atoms.items():
            element = ''.join([i for i in atom['name'] if not i.isdigit()])  
            attributes = {
                # Name is the atom type_molname
                "name": atom['name'] + '_' + mol_name,
                "class": str(i + 1),
                "element": element,
                "mass": str(atom['mass']) 
            }
            ET.SubElement(AtomTypes, "Type", **attributes)

        # Residues section
        Residues = ET.SubElement(ForceField, "Residues")
        Residue = ET.SubElement(Residues, "Residue", name=mol_name)
        for atom_id, atom_data in atoms.items():
            ET.SubElement(Residue, "Atom", name=atom_data['name'], type=atom_data['name'] + '_' + mol_name, charge=str(atom_data['charge']))
        for bond_id, bond_data in bonds.items():
            ET.SubElement(Residue, "Bond", atomName1=atoms[bond_id[0]]['name'], atomName2=atoms[bond_id[1]]['name'])

        # Bonds section
        Bonds = ET.SubElement(ForceField, "HarmonicBondForce")
        for bond_id, bond_data in bonds.items():
            attributes = {
                "class1": str(bond_id[0] + 1),
                "class2": str(bond_id[1] + 1),
                "length": str(bond_data['equilibrium']),
                "k": str(bond_data['force_constant'])
            }
            ET.SubElement(Bonds, "Bond", **attributes)

        # Angles section
        Angles = ET.SubElement(ForceField, "HarmonicAngleForce")
        for angle_id, angle_data in self.angles.items():
            attributes = {
                "class1": str(angle_id[0] + 1),
                "class2": str(angle_id[1] + 1),
                "class3": str(angle_id[2] + 1),
                "angle": str(angle_data['equilibrium'] * np.pi / 180),
                "k": str(angle_data['force_constant'])
            }
            ET.SubElement(Angles, "Angle", **attributes)

        # Periodic Dihedrals section
        Dihedrals = ET.SubElement(ForceField, "PeriodicTorsionForce")
        for dihedral_id, dihedral_data in self.dihedrals.items():
            if dihedral_data['type'] == 'RB':
                continue
            attributes = {
                "class1": str(dihedral_id[0] + 1),
                "class2": str(dihedral_id[1] + 1),
                "class3": str(dihedral_id[2] + 1),
                "class4": str(dihedral_id[3] + 1),
                "periodicity1": str(dihedral_data['periodicity']),
                "phase1": str(dihedral_data['phase'] * np.pi / 180),
                "k1": str(dihedral_data['barrier'])
            }
            ET.SubElement(Dihedrals, "Proper", **attributes)

        # RB Dihedrals section
        RB_Dihedrals = ET.SubElement(ForceField, "RBTorsionForce")
        for dihedral_id, dihedral_data in self.dihedrals.items():
            if dihedral_data['type'] == 'Fourier':
                continue
            attributes = {
                "class1": str(dihedral_id[0] + 1),
                "class2": str(dihedral_id[1] + 1),
                "class3": str(dihedral_id[2] + 1),
                "class4": str(dihedral_id[3] + 1),
                "c0": str(dihedral_data['RB_coefficients'][0]),
                "c1": str(dihedral_data['RB_coefficients'][1]),
                "c2": str(dihedral_data['RB_coefficients'][2]),
                "c3": str(dihedral_data['RB_coefficients'][3]),
                "c4": str(dihedral_data['RB_coefficients'][4]),
                "c5": str(dihedral_data['RB_coefficients'][5])
            }
            ET.SubElement(RB_Dihedrals, "Proper", **attributes)

        # Improper Dihedrals section
        Impropers = ET.SubElement(ForceField, "PeriodicTorsionForce")
        for improper_id, improper_data in self.impropers.items():

            # The order of the atoms is defined in the OpenMM documentation
            # http://docs.openmm.org/latest/userguide/application/06_creating_ffs.html

            attributes = {
                "class1": str(improper_id[1] + 1),
                "class2": str(improper_id[0] + 1),
                "class3": str(improper_id[2] + 1),
                "class4": str(improper_id[3] + 1),
                "periodicity1": str(improper_data['periodicity']),
                "phase1": str(improper_data['phase'] * np.pi / 180),
                "k1": str(improper_data['barrier'])
            }
            ET.SubElement(Impropers, "Improper", **attributes)

        # NonbondedForce section
        NonbondedForce = ET.SubElement(ForceField, "NonbondedForce", coulomb14scale=str(self.fudgeQQ), lj14scale=str(self.fudgeLJ))
        for atom_id, atom_data in atoms.items():
            attributes = {
                "type": atom_data['name'] + '_' + mol_name,
                "charge": str(atom_data['charge']),
                "sigma": str(atom_data['sigma']),
                "epsilon": str(atom_data['epsilon'])
            }
            ET.SubElement(NonbondedForce, "Atom", **attributes)

        # Generate the tree and write to file
        tree = ET.ElementTree(ForceField)
        rough_string = ET.tostring(ForceField, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        indented_string = reparsed.toprettyxml(indent="    ")  

        with open(filename, 'w') as output_file:
            output_file.write(indented_string)

    def write_gro(self, gro_file, mol_name=None, gro_precision=3):
        """
        Writes a GRO file with the original coordinates.

        :param gro_file:
            The GRO file path.
        :param mol_name:
            The name of the molecule.
        :param gro_precision:
            The number of decimal places in gro file.
        """

        gro_filename = str(gro_file)
        if mol_name is None:
            mol_name = Path(self.molecule_name).stem
        res_name = 'MOL'

        coords_in_nm = self.molecule.get_coordinates_in_angstrom() * 0.1

        with open(gro_filename, 'w') as f_gro:
            # Header
            f_gro.write(f'GRO file of {mol_name}, generated by VeloxChem\n')
            f_gro.write(f'{len(self.atoms):>5d}\n')

            ndec = gro_precision

            # Atoms
            for i, atom in self.atoms.items():
                atom_name = atom['name']
                line_str = f'{1:>5d}{res_name:<5s}{atom_name:<5s}{i + 1:>5d}'
                for d in range(3):
                    line_str += f'{coords_in_nm[i][d]:{ndec+5}.{ndec}f}'
                line_str += '\n'
                f_gro.write(line_str)

            # Box
            box_dimension = 10.0
            line_str = f'{box_dimension:10.5f}' * 3 + '\n'
            f_gro.write(line_str)

    def write_pdb(self, pdb_file, mol_name=None):
        """
        Writes a PDB file with the original coordinates.

        :param pdb_file:
            The PDB file path.
        :param mol_name:
            The name of the molecule.
        """

        #PDB format from http://deposit.rcsb.org/adit/docs/pdb_atom_format.html:

        # COLUMNS        DATA TYPE       CONTENTS
        # --------------------------------------------------------------------------------
        #  1 -  6        Record name     "HETATM" or "ATOM  "
        #  7 - 11        Integer         Atom serial number.
        # 13 - 16        Atom            Atom name.
        # 17             Character       Alternate location indicator.
        # 18 - 20        Residue name    Residue name.
        # 22             Character       Chain identifier.
        # 23 - 26        Integer         Residue sequence number.
        # 27             AChar           Code for insertion of residues.
        # 31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
        # 39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
        # 47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
        # 55 - 60        Real(6.2)       Occupancy (Default = 1.0).
        # 61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
        # 73 - 76        LString(4)      Segment identifier, left-justified.
        # 77 - 78        LString(2)      Element symbol, right-justified.
        # 79 - 80        LString(2)      Charge on the atom.
        
        pdb_filename = str(pdb_file)
        if mol_name is None:
            mol_name = Path(self.molecule_name).stem

        coords_in_angstrom = self.molecule.get_coordinates_in_angstrom()
        molecule_elements = self.molecule.get_labels()

        with open(pdb_filename, 'w') as f_pdb:
            # Header
            f_pdb.write(f'TITLE     PDB file of {mol_name}, generated by VeloxChem\n')
            f_pdb.write('MODEL        1\n')

            # Atoms
            for i, (atom, element) in enumerate(zip(self.atoms.values(), molecule_elements), 1):
                atom_name = atom['name']
                occupancy = 1.00
                temp_factor = 0.00
                element_symbol = element[:2].rjust(2)

                # Format string from https://cupnet.net/pdb-format/ 
                line_str = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}".format(
                    'HETATM', i, atom_name[:4],'',mol_name[:3],'A', 1,'', 
                    coords_in_angstrom[i-1][0], coords_in_angstrom[i-1][1], coords_in_angstrom[i-1][2], 
                    occupancy, temp_factor, 
                    element_symbol)

                f_pdb.write(line_str + '\n')
                
            # CONECT section in the PDB file stating the connectivity.
            # Required by OpenMM to correctly assign topology.bonds
            for (i, j) in self.bonds:
                f_pdb.write(f'CONECT{i+1:>5}{j+1:>5}\n')

            f_pdb.write('TER\n')
            f_pdb.write('ENDMDL\n')
            f_pdb.write('END\n')

    def write_gromacs_files(self,
                            filename,
                            mol_name=None,
                            amber_ff=None,
                            water_model=None,
                            gro_precision=3):
        """
        Writes all the needed files for a MD simulation with GROMACS.

        :param filename:
            The name of the molecule.
        :param mol_name:
            The name of the molecule.
        :param amber_ff:
            The name of the Amber force field.
        :param water_model:
            The name of the water model.
        :param gro_precision:
            The number of decimal places in gro file.
        """

        if mol_name is None:
            mol_name = Path(self.molecule_name).stem

        itp_file = Path(filename).with_suffix('.itp')
        top_file = Path(filename).with_suffix('.top')
        gro_file = Path(filename).with_suffix('.gro')

        self.write_itp(itp_file, mol_name)
        self.write_top(top_file, itp_file, mol_name, amber_ff, water_model)
        self.write_gro(gro_file, mol_name, gro_precision)

    def write_openmm_files(self, filename, mol_name=None):
        """
        Writes all the needed files for a MD simulation with OpenMM.

        :param filename:
            The name of the molecule.
        :param mol_name:
            The name of the molecule.
        """

        if mol_name is None:
            mol_name = Path(self.molecule_name).stem

        xml_file = Path(filename).with_suffix('.xml')
        pdb_file = Path(filename).with_suffix('.pdb')

        self.generate_residue_xml(xml_file, mol_name)
        self.write_pdb(pdb_file, mol_name)

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

    def dihedral_correction(self, i, kT=None):
        """
        Corrects dihedral parameters.

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

        dih = self.target_dihedrals[i]
        geom = self.scan_geometries[i]
        qm_scan = self.scan_energies[i]
        angles = self.scan_dih_angles[i]

        dih_str = f'{dih[0]+1}-{dih[1]+1}-{dih[2]+1}-{dih[3]+1}'
        self.ostream.print_info(f'Fitting dihedral angle {dih_str}...')
        self.ostream.flush()

        # MM scan with dihedral parameters set to zero

        self.set_dihedral_parameters(dih, [0.0 for x in range(6)])

        mm_zero_scan = self.perform_mm_scan(dih,
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

            coef_list = coef.tolist()
        else:
            coef_list = None

        coef_list = self.comm.bcast(coef_list, root=mpi_master())

        self.set_dihedral_parameters(dih, coef_list)

        # write itp and gro files

        mol_name = Path(self.molecule_name).stem

        workdir = Path('.') if self.workdir is None else self.workdir

        self.ffversion += 1
        new_itp_fname = str(workdir / f'{mol_name}_{self.ffversion:02d}.itp')
        new_top_fname = str(workdir / f'{mol_name}_{self.ffversion:02d}.top')

        if self.rank == mpi_master():
            self.write_itp(new_itp_fname, mol_name)
            self.write_top(new_top_fname, new_itp_fname, mol_name)
        self.comm.barrier()

        self.ostream.print_info('...done.')
        self.ostream.print_blank()
        self.ostream.print_info(
            f'Generated new topology file: {Path(new_top_fname).name}')
        self.ostream.print_blank()
        self.ostream.flush()

    def validate_force_field(self, i):
        """
        Validates force field by RMSD of dihedral potentials.

        :param i:
            The index of the target dihedral.

        :return:
            A dictionary containing the results of validation.
        """

        self.ostream.print_info('Validating force field ...')
        self.ostream.print_blank()
        self.ostream.flush()

        dih = self.target_dihedrals[i]

        dih_str = f'{dih[0]+1}-{dih[1]+1}-{dih[2]+1}-{dih[3]+1}'
        self.ostream.print_info(f'  Target dihedral angle: {dih_str}')
        self.ostream.print_blank()

        geom = self.scan_geometries[i]
        angles = self.scan_dih_angles[i]
        mm_scan = self.perform_mm_scan(dih, geom, angles)
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
                        dihedral,
                        geometries,
                        angles,
                        print_energies=True):
        """
        Performs MM scan of a specific dihedral.

        :param dihedral:
            The dihedral (list of four atom ids).
        :param geometries:
            The scanned geometris for this dihedral.
        :param angles:
            The scanned angles for this dihedral.
        """

        # select scan angles and geometries from QM data

        if print_energies:
            self.ostream.print_info('      Dihedral           MM energy')
            self.ostream.print_info('  --------------------------------')
            self.ostream.flush()

        energies = []

        for i, (geom, angle) in enumerate(zip(geometries, angles)):
            # energy minimization with dihedral constraint
            constraints = [
                'set dihedral {} {} {} {} {}'.format(dihedral[0] + 1,
                                                     dihedral[1] + 1,
                                                     dihedral[2] + 1,
                                                     dihedral[3] + 1, angle)
            ]
            pot_energy = self.minimize_mm_energy(geom, constraints)
            # convert to kJ/mol
            pot_energy *= 4.184 * hartree_in_kcalpermol()
            energies.append(pot_energy)

            if print_energies:
                self.ostream.print_info(
                    f'  {angle:8.1f} deg {pot_energy:12.3f} kJ/mol')
                self.ostream.flush()

        return energies

    def minimize_mm_energy(self, molecule, constraints):
        """
        Minimizes MM energy of a topology using MM driver.

        :param molecule:
            The molecule.
        :param constraints:
            The constraints.
        """

        mm_drv = MMDriver(self.comm, self.ostream)
        mm_drv.load_force_field(self)

        self.ostream.mute()
        grad_drv = MMGradientDriver(mm_drv)
        opt_drv = OptimizationDriver(grad_drv)
        opt_drv.constraints = constraints
        opt_results = opt_drv.compute(molecule)
        final_mol = Molecule.read_xyz_string(opt_results['final_geometry'])
        self.ostream.unmute()

        mm_drv.compute(final_mol)

        return mm_drv.get_energy()

    def set_dihedral_parameters(self, dihedral_key, coefficients):
        """
        Sets dihedral parameters of topology.

        :param dihedral_key:
            The dihedral (list of four atom ids).
        :param coefficients:
            The coefficients for Ryckaert-Bellemans funciton.
        """

        dih = self.dihedrals[tuple(dihedral_key)]

        if dih['type'] == 'Fourier' and dih['multiple']:
            comment = dih['comment'][0]
        else:
            comment = dih['comment']

        self.dihedrals[tuple(dihedral_key)] = {
            'type': 'RB',
            'RB_coefficients': list(coefficients),
            'comment': comment
        }

        dih_keys_to_remove = []
        for i, j, k, l in self.dihedrals:
            if [i, j, k, l] == dihedral_key:
                continue
            if [j, k] == dihedral_key[1:3] or [k, j] == dihedral_key[1:3]:
                dih_keys_to_remove.append((i, j, k, l))
        for dih_key in dih_keys_to_remove:
            del self.dihedrals[dih_key]

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
        plt.xlabel('dihedral angle {}-{}-{}-{}'.format(dihedral_indices[0] + 1,
                                                       dihedral_indices[1] + 1,
                                                       dihedral_indices[2] + 1,
                                                       dihedral_indices[3] + 1))
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
