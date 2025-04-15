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

from mpi4py import MPI
from pathlib import Path, PurePath
import numpy as np
import tempfile
import sys
import re
import xml.etree.ElementTree as ET
from xml.dom import minidom
from collections import defaultdict

from .veloxchemlib import mpi_master, bohr_in_angstrom, hartree_in_kjpermol
from .atomtypeidentifier import AtomTypeIdentifier
from .molecule import Molecule
from .molecularbasis import MolecularBasis
from .outputstream import OutputStream
from .respchargesdriver import RespChargesDriver
from .mmdriver import MMDriver
from .mmgradientdriver import MMGradientDriver
from .scfrestdriver import ScfRestrictedDriver
from .optimizationdriver import OptimizationDriver
from .inputparser import parse_input
from .errorhandler import assert_msg_critical, safe_arccos
from .seminario import Seminario
from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .uffparameters import get_uff_parameters
from .tmparameters import get_tm_parameters
from .waterparameters import get_water_parameters
from .environment import get_data_path


class MMForceFieldGenerator:
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
        self.molecule_name = 'MOL'
        self.scan_xyz_files = None
        self.atom_types = None
        self.rotatable_bonds = []

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
        self.topology_update_flag = False

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

        self.workdir = None

        # resp settings
        self.resp_dict = None

        self.keep_files = True

        # UFF parameters
        self.uff_parameters = get_uff_parameters()

        # TM parameters
        self.tm_parameters = get_tm_parameters()

        # Water parameters
        self.water_parameters = get_water_parameters()

        # Summary of fitting
        self.fitting_summary = None

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
                    'MMForceFieldGenerator: unrecognized force field ' +
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
            'MMForceFieldGenerator.compute: scan_xyz_files not defined ')

        assert_msg_critical(
            self.atom_types is not None,
            'MMForceFieldGenerator.compute: atom_types not defined ')

        assert_msg_critical(
            len(self.atom_types) == molecule.number_of_atoms(),
            'MMForceFieldGenerator.compute: inconsistent number of atom_types')

        # read QM scan

        # TODO: enable use of multiple scan files
        assert_msg_critical(
            len(self.scan_xyz_files) == 1,
            'MMForceFieldGenerator.compute: only single scan file is supported for now')

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

        use_temp_dir = (self.workdir is None)

        if use_temp_dir:
            try:
                temp_dir = tempfile.TemporaryDirectory(
                    ignore_cleanup_errors=True)
            except TypeError:
                temp_dir = tempfile.TemporaryDirectory()
            self.workdir = Path(temp_dir.name)

        if self.original_top_file is None:
            original_itp_file = self.workdir / (mol_name + '.itp')
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

        # fit dihedral potentials

        n_rounds = self.n_rounds if len(self.scan_dih_angles) > 1 else 1
        for i_round in range(n_rounds):
            for i, dih in enumerate(self.target_dihedrals):
                dih_central_bond = [dih[1] + 1, dih[2] + 1]
                scan_file = str(Path(inp_dir) / self.scan_xyz_files[i])
                self.reparameterize_dihedrals(
                    dih_central_bond, scan_file=scan_file, fit_extrema=False)

        # save output files

        if self.rank == mpi_master():
            self.write_gromacs_files(original_top_file)
        self.comm.barrier()

        if self.rank == mpi_master() and self.keep_files:
            out_dir = Path(self.molecule_name + '_files')
            for ftype in ['itp', 'top']:
                fname = mol_name + f'.{ftype}'
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

    def reparameterize_dihedrals(self,
                                 rotatable_bond,
                                 scan_file=None,
                                 scf_drv=None,
                                 basis=None,
                                 scf_results=None,
                                 scan_range=[0, 360],
                                 n_points=19,
                                 scan_verbose=False,
                                 visualize=False,
                                 fit_extrema=False,
                                 initial_validation=True,
                                 show_diff=False,
                                 verbose=True):
        """
        Changes the dihedral constants for a specific rotatable bond in order to
        fit the QM scan.

        :param rotatable_bond:
            The list of indices of the rotatable bond. (1-indexed)
        :param scan_file:
            The file with the QM scan. If None is provided it a QM scan will be performed.
        :param scf_drv:
            The SCF driver. If None is provided it will use HF.
        :param basis:
            The AO basis set. If None is provided it will use 6-31G*.
        :param scf_results:
            The dictionary containing converged SCF results.
        :param scan_range:
            List with the range of dihedral angles. Default is [0, 360].
        :param n_points:
            The number of points to be calculated. Default is 19.
        :param scan_verbose:
            Whether the QM scan should print all the information.
        :param visualize:
            Whether the dihedral scans should be visualized.
        :param fit_extrema:
            Whether the dihedral parameters should be fitted to the QM scan extrema.
        """

        try:
            from scipy.optimize import least_squares
            from scipy.signal import argrelextrema
        except ImportError:
            error_msg = 'Scipy is required for reparameterize_dihedrals.'
            assert_msg_critical(False, error_msg)

        # If scan file is provided, read it
        if scan_file is not None:
            with open(scan_file, 'r') as f:
                lines = f.readlines()

            for line in lines:
                if line.startswith('Scan'):
                    # Scan Cycle 1/19 ; Dihedral 3-4-6-7 = 0.00 ; Iteration 17 Energy -1089.05773546
                    # Extract the dihedral indices
                    scanned_dih = [int(i) for i in line.split('Dihedral')[1].split()[0].split('-')]
                    central_atom_1 = scanned_dih[1] - 1
                    central_atom_2 = scanned_dih[2] - 1

                    # Check if the rotatable bond matches the scan file
                    if sorted([central_atom_1, central_atom_2]) != sorted(
                            [rotatable_bond[0] - 1, rotatable_bond[1] - 1]):
                        raise ValueError('The rotatable bond does not match the scan file.')

                    break

        else:
            central_atom_1 = rotatable_bond[0] - 1
            central_atom_2 = rotatable_bond[1] - 1

        # Identify the dihedral indices for the rotatable bond
        dihedral_indices = []
        dihedral_types = []

        for (i,j,k,l), dihedral in self.dihedrals.items():
            if (j == central_atom_1 and k == central_atom_2) or (j == central_atom_2 and k == central_atom_1):
                dihedral_indices.append([i,j,k,l])
                dihedral_types.append(dihedral['comment'])

        # Transform the dihedral indices to 1-indexed for printing
        dihedral_indices_print = [[i+1,j+1,k+1,l+1] for i,j,k,l in dihedral_indices]
        
        # Print a header
        header = 'VeloxChem Dihedral Reparameterization'
        self.ostream.print_header(header)
        self.ostream.print_header('=' * len(header))
        self.ostream.print_blank()

        self.ostream.print_info(f'Rotatable bond selected: {rotatable_bond[0]}-{rotatable_bond[1]}')

        # Print the dihedral indices *in 1 index* and types
        self.ostream.print_info(f'Dihedrals involved:{dihedral_indices_print}')

        # If the scan file is provided, read it
        if scan_file is not None:
            # Check if the name of the file is correct 1-2-3-4.xyz
            file_path = Path(scan_file)
            # Note: Sometimes geomeTRIC writes the dihedral atoms in reverse order
            if file_path.stem not in [
                    f"{scanned_dih[0]}-{scanned_dih[1]}-{scanned_dih[2]}-{scanned_dih[3]}",
                    f"{scanned_dih[3]}-{scanned_dih[2]}-{scanned_dih[1]}-{scanned_dih[0]}"]:
                raise ValueError('The scan file name does not match the dihedral indices. Format should be 1-2-3-4.xyz')
            self.read_qm_scan_xyz_files([scan_file])
        else:
            # Perform a QM scan
            self.ostream.print_info('No scan file provided. Performing QM scan...')
        
            if scf_drv is None:
                scf_drv = ScfRestrictedDriver()
                if not scan_verbose:
                    scf_drv.ostream.mute()
                self.ostream.print_info('SCF driver not provided. Using default: RHF')

            if basis is None:
                basis = MolecularBasis.read(self.molecule,'6-31G*')
                self.ostream.print_info('Basis set not provided. Using dafault: 6-31G*')

            # Perform the reference SCF calculation
            if scf_results is None:
                self.ostream.print_info('Performing SCF calculation...')
                self.ostream.flush()
                scf_results = scf_drv.compute(self.molecule, basis)
            
            # Take one of the dihedrals to perform the scan
            reference_dih = dihedral_indices[0]
            reference_dih_name = f"{reference_dih[0] + 1}-{reference_dih[1] + 1}-{reference_dih[2] + 1}-{reference_dih[3] + 1}"

            opt_drv = OptimizationDriver(scf_drv)
            
            if not scan_verbose:
                opt_drv.ostream.mute()

            constraint = f"scan dihedral {reference_dih[0]+1} {reference_dih[1]+1} {reference_dih[2]+1} {reference_dih[3]+1} {scan_range[0]} {scan_range[1]} {n_points}"
            opt_drv.constraints = [constraint]

            # Scan the dihedral
            self.ostream.print_info(f'Dihedral: {reference_dih_name} taken for scan...')
            self.ostream.print_info(f'Angle range: {scan_range[0]}-{scan_range[1]}')
            self.ostream.print_info(f'Number of points: {n_points}')
            self.ostream.print_info(f'Scanning dihedral {reference_dih_name}...')
            self.ostream.flush()

            opt_drv.compute(self.molecule, basis, scf_results)

            self.ostream.print_info('Scan completed.')

            # Change the default file name to the dihedral indices
            file_path = Path("scan-final.xyz")
            file_path.rename(f"{reference_dih_name}.xyz")

            # Read the QM scan
            self.read_qm_scan_xyz_files([f"{reference_dih_name}.xyz"])

        # Group dihedrals by their types
        dihedral_groups = defaultdict(list)
        for i, j, k, l in dihedral_indices:
            dihedral = self.dihedrals[(i, j, k, l)]
            # If the dihedral is multiple add all the types
            if dihedral['multiple']:
                # Handle multiple dihedrals
                for idx, dihedral_type in enumerate(dihedral['comment']):
                    dihedral_groups[dihedral_type].append((i, j, k, l))
            else:
                dihedral_type = dihedral['comment']
                dihedral_groups[dihedral_type].append((i, j, k, l))

        # Extract initial parameters for each dihedral instance
        barriers = []
        phases = []
        periodicities = []
        for i, j, k, l in dihedral_indices:
            dihedral = self.dihedrals[(i, j, k, l)]
            if dihedral['multiple']:
                # Handle multiple dihedrals
                for idx, dihedral_type in enumerate(dihedral['comment']):
                    barrier = dihedral['barrier'][idx]
                    phase = dihedral['phase'][idx]
                    periodicity = abs(dihedral['periodicity'][idx])
                    barriers.append(barrier)
                    phases.append(phase)
                    periodicities.append(periodicity)

            else:
                barrier = dihedral['barrier']
                phase = dihedral['phase']
                periodicity = dihedral['periodicity']
                barriers.append(barrier)
                phases.append(phase)
                periodicities.append(periodicity)

        # Make a copy of the original dihedrals dictionary
        # original_dihedrals_dict = self.dihedrals.copy()

        # Convert to NumPy arrays
        barriers = np.array(barriers, dtype=float)
        phases_array = np.array(phases, dtype=float)
        periodicities_array = np.array(periodicities, dtype=float)
        phases_array_rad = np.deg2rad(phases_array)
        dihedral_angles_rad = np.deg2rad(self.scan_dih_angles[0])

        # Set the dihedral barriers to zero for the scan
        for i, j, k, l in dihedral_indices:
            if self.dihedrals[(i, j, k, l)]['multiple']:
                for idx, dihedral_type in enumerate(self.dihedrals[(i, j, k, l)]['comment']):
                    self.dihedrals[(i, j, k, l)]['barrier'][idx] = 0.0
            else:
                self.dihedrals[(i, j, k, l)]['barrier'] = 0.0

        # Scan the dihedral
        self.ostream.print_info(
            'Performing dihedral scan for MM baseline ' +
            'by excluding the involved dihedral barriers...')
        self.ostream.print_blank()
        self.ostream.flush()

        initial_data = self.validate_force_field(0, verbose=verbose)

        qm_energies = np.array(initial_data['qm_scan_kJpermol'])
        mm_baseline = np.array(initial_data['mm_scan_kJpermol'])

        def dihedral_potential(phi, barriers, phases_rad, periodicities):

            total_potential = np.zeros_like(phi)

            for barrier, phase_rad, periodicity in zip(barriers, phases_rad, periodicities):
                total_potential += barrier * (1 + np.cos(periodicity * phi - phase_rad))

            return total_potential
        
        def extract_maxima(barriers, qm_energies):

            # Identify maxima indices for QM energies
            qm_maxima_indices = argrelextrema(qm_energies, np.greater)[0]   

            # Identify the minima for the QM energies
            qm_minima_indices = argrelextrema(qm_energies, np.less)[0]

            # Use first and last points as maxima (usual behavior for cosine dihedrals)
            qm_extrema = [0, len(qm_energies) - 1]

            # Build the arrays with the maxima and minima
            qm_maxima_indices = (list(qm_maxima_indices) +
                                 list(qm_minima_indices) +
                                 qm_extrema)

            qm_maxima_indices = sorted(list(set(qm_maxima_indices)))

            return qm_maxima_indices

        def objective_function(barriers_to_fit, qm_maxima_indices):

            # Update the dihedral parameters in the force field
            dihedral_energies = dihedral_potential(
                dihedral_angles_rad,
                barriers_to_fit[:-1],
                phases_array_rad,
                periodicities_array
            )
            mm_energies_fit = dihedral_energies + mm_baseline

            # Relative energies
            mm_energies_fit_rel = np.array(mm_energies_fit) - np.min(mm_energies_fit)
            qm_energies_rel = np.array(qm_energies) - np.min(qm_energies)

            if fit_extrema:

                # Match QM and MM maxima using nearest x-axis values
                matched_qm_maxima = []
                matched_mm_maxima = []
                for qm_idx in qm_maxima_indices:
                    matched_qm_maxima.append(qm_energies_rel[qm_idx])
                    matched_mm_maxima.append(mm_energies_fit_rel[qm_idx])

                # Convert to arrays
                qm_energies_rel = np.array(matched_qm_maxima) - min(matched_qm_maxima)
                mm_energies_fit_rel = np.array(matched_mm_maxima) - min(matched_mm_maxima)

            # Residuals
            # TODO: consider using weights
            residuals = (qm_energies_rel - mm_energies_fit_rel)
            residuals += barriers_to_fit[-1]

            # Return the squared residuals for optimization
            return residuals**2

        # Print initial barriers
        self.ostream.print_info(f"Dihedral barriers {barriers} will be used as initial guess.")

        # Store the original barriers and perform the initial validation
        original_barriers = barriers.copy()

        if initial_validation:
            self.ostream.print_info('Validating the initial force field...')
            self.ostream.print_blank()
            self.ostream.flush()

            mm_energies = dihedral_potential(
                dihedral_angles_rad,
                original_barriers,
                phases_array_rad,
                periodicities_array
            )
            mm_energies += mm_baseline

            fitted_dihedral_results = {
                'dihedral_indices': list(initial_data['dihedral_indices']),
                'dihedral_angles': list(initial_data['dihedral_angles']),
                'mm_scan_kJpermol': np.copy(mm_energies) - np.min(mm_energies),
                'qm_scan_kJpermol': np.copy(qm_energies) - np.min(qm_energies),
            }

            self.print_validation_summary(fitted_dihedral_results)

            if visualize:
                self.visualize(fitted_dihedral_results, show_diff=show_diff)

        self.ostream.print_info('Fitting the dihedral parameters...')
        if fit_extrema:
            self.ostream.print_info('Only minimum/maximum points are used for fitting.')

        # Use the original barriers as the initial guess
        # Note that we add an additional parameter for shifting QM and MM
        # energies
        initial_guess = np.zeros(original_barriers.size + 1)
        initial_guess[:-1] = original_barriers[:]
        
        # Extract maxima for QM and MM energies
        if fit_extrema:
            qm_maxima_indices = extract_maxima(barriers, qm_energies)
            args = (qm_maxima_indices,)
        else:
            args = (None,)

        # Set as bound that the barriers should be positive
        bounds = (0, np.inf)

        # Perform the optimization
        self.ostream.print_info('Optimizing dihedral via least squares fitting...')

        result = least_squares(
            fun=objective_function,
            x0=initial_guess,
            jac='2-point',
            method='trf',
            bounds=bounds,
            args=args,
            verbose=0,
            xtol=1e-8,
            ftol=1e-8,
            gtol=1e-8,
            max_nfev=2000,
        )

        # Update the force field parameters with the optimized barriers
        # The additional parameter for shifting energy is no longer needed
        fitted_barriers = np.copy(result.x[:-1])

        # Get the final MM energy profile
        fitted_dihedral_energies = dihedral_potential(
            dihedral_angles_rad,
            fitted_barriers,
            phases_array_rad,
            periodicities_array
        ) + mm_baseline

        fitted_dihedral_energies = (np.array(fitted_dihedral_energies) -
                                    np.min(fitted_dihedral_energies))

        # Adjust the barriers to the dihedral types
        param_tuples = list(zip(
            original_barriers, phases_array_rad, periodicities_array))

        if len(set(param_tuples)) != len(param_tuples):
            for unique_param_tuple in set(param_tuples):
                barrier_indices = [
                    idx for idx in range(len(original_barriers))
                    if unique_param_tuple == (original_barriers[idx],
                                              phases_array_rad[idx],
                                              periodicities_array[idx])
                ]
                mean_barrier = np.mean(fitted_barriers[barrier_indices])
                for idx in barrier_indices:
                    fitted_barriers[idx] = mean_barrier

        # List of the fitted barriers
        fit_barrier_to_print = fitted_barriers.copy()
        self.ostream.print_info(f'New fitted barriers: {fit_barrier_to_print}')

        # If there are multiple dihedrals, group them in list of lists
        fitted_barriers_grouped = []

        for i, j, k, l in dihedral_indices:
            dihedral = self.dihedrals[(i, j, k, l)]
            if dihedral['multiple']:
                number_dihedrals = len(dihedral['comment'])
                fitted_barriers_grouped.append(fitted_barriers[:number_dihedrals])
                fitted_barriers = np.array(fitted_barriers[number_dihedrals:])
            else:
                fitted_barriers_grouped.append([fitted_barriers[0]])
                fitted_barriers = np.array(fitted_barriers[1:])

        fitted_barriers = fitted_barriers_grouped

        # Assign the fitted barriers to the dihedrals
        for idx, (i, j, k, l) in enumerate(dihedral_indices):
            dihedral = self.dihedrals[(i, j, k, l)]
            if dihedral['multiple']:
                # make sure that the barrier is a list of float and not a list of numpy.float
                dihedral['barrier'] = [float(x) for x in fitted_barriers[idx]]
                for comment_idx in range(len(dihedral['comment'])):
                    if not dihedral['comment'][comment_idx].endswith(' (fitted)'):
                        dihedral['comment'][comment_idx] += ' (fitted)'
            else:
                # make sure that the barrier is float and not numpy.float
                dihedral['barrier'] = float(fitted_barriers[idx][0])
                if not dihedral['comment'].endswith(' (fitted)'):
                    dihedral['comment'] += ' (fitted)'

        # Validate the fitted parameters
        self.ostream.print_info('Validating the fitted force field...')
        self.ostream.print_blank()
        self.ostream.flush()

        fitted_dihedral_results = {
            'dihedral_indices': list(initial_data['dihedral_indices']),
            'dihedral_angles': list(initial_data['dihedral_angles']),
            'mm_scan_kJpermol': np.copy(fitted_dihedral_energies),
            'qm_scan_kJpermol': np.copy(qm_energies) - np.min(qm_energies),
        }

        self.print_validation_summary(fitted_dihedral_results)

        if visualize:
            self.visualize(fitted_dihedral_results, show_diff=show_diff)

        self.ostream.print_info('Dihedral MM parameters have been reparameterized and updated in the topology.')
        self.ostream.flush()
   
    def read_qm_scan_xyz_files(self, scan_xyz_files, inp_dir=None):
        """
        Reads QM scan xyz files.

        :param scan_xyz_files:
            The list of xyz files from QM scan.
        """

        self.scan_dih_angles = []
        self.scan_energies = []
        self.scan_geometries = []
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
                    'MMForceFieldGenerator.read_qm_scan_xyz_files: ' +
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

    @staticmethod
    def get_gaff_data_dict():
        """
        Gets GAFF data into a dictionary.
        """

        data = {
            'version': None,
            'refs': [],
            'atom_types': [],
            'bonds': [],
            'angles': [],
            'dihedrals': [],
            'impropers': [],
            'fudgeQQ': None,
            'fudgeLJ': None,
        }

        data_path = get_data_path()
        gaff_path = Path(data_path) / 'gaff-2.11.xml'
        tree = ET.parse(str(gaff_path))
        root = tree.getroot()

        for node in root:

            if node.tag == 'Info':
                for subnode in node:
                    if subnode.tag == 'Source':
                        info = subnode.attrib
                        data['version'] = info['Source']
                    elif subnode.tag == 'Reference':
                        data['refs'].append(subnode.text)

            elif node.tag == 'AtomTypes':
                for subnode in node:
                    if subnode.tag == 'Type':
                        info = subnode.attrib
                        data['atom_types'].append(info)

            elif node.tag == 'HarmonicBondForce':
                for subnode in node:
                    if subnode.tag == 'Bond':
                        info = subnode.attrib
                        data['bonds'].append(info)

            elif node.tag == 'HarmonicAngleForce':
                for subnode in node:
                    if subnode.tag == 'Angle':
                        info = subnode.attrib
                        data['angles'].append(info)

            elif node.tag == 'PeriodicTorsionForce':
                for subnode in node:
                    if subnode.tag == 'Proper':
                        info = subnode.attrib
                        data['dihedrals'].append(info)
                    elif subnode.tag == 'Improper':
                        info = subnode.attrib
                        data['impropers'].append(info)

            elif node.tag == 'NonbondedForce':
                info = node.attrib
                data['fudgeQQ'] = info['coulomb14scale']
                data['fudgeLJ'] = info['lj14scale']

                for subnode in node:
                    if subnode.tag == 'Atom':
                        info = subnode.attrib
                        for idx in range(len(data['atom_types'])):
                            if data['atom_types'][idx]['class'] == info['class']:
                                data['atom_types'][idx]['sigma'] = info['sigma']
                                data['atom_types'][idx]['epsilon'] = info['epsilon']

        return data

    def create_topology(self, molecule, basis=None, scf_results=None, resp=True, water_model=None, use_xml=True):
        """
        Analyzes the topology of the molecule and create dictionaries
        for the atoms, bonds, angles, dihedrals, impropers and pairs.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_results:
            A converged SCF result.
        :param resp:
            If RESP charges should be computed.
            If False partial charges will be set to zero.
        """

        ff_data_dict = None
        ff_data_lines = None

        # Read the force field data

        if use_xml:
            ff_data_dict = self.get_gaff_data_dict()
        elif self.force_field_data is None:
            ff_data_lines = self.get_gaff_data_lines()
        else:
            with open(self.force_field_data, 'r') as ff_data:
                ff_data_lines = ff_data.readlines()

        # check GAFF version

        gaff_version = None

        if use_xml:
            ff_data_version = ff_data_dict['version'].replace('gaff-', '').replace('.dat', '')
            if '.' in ff_data_version:
                version_major = ff_data_version.split('.')[0]
                version_minor = ff_data_version.split('.')[1]
                if version_major.isdigit() and version_minor.isdigit():
                    gaff_version = f'{version_major}.{version_minor}'

        elif ff_data_lines[0].startswith('AMBER General Force Field'):
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

        if self.topology_update_flag:
            self.atom_types = atomtypeidentifier.generate_gaff_atomtypes(
                self.molecule, self.connectivity_matrix)
            # The partial charges have to be recalculated
            self.partial_charges = None
        else:
            self.atom_types = atomtypeidentifier.generate_gaff_atomtypes(
                self.molecule)
            self.connectivity_matrix = np.copy(
                atomtypeidentifier.connectivity_matrix)
            
        atomtypeidentifier.identify_equivalences()
        
        self.atom_info_dict = atomtypeidentifier.atom_info_dict
        ## TODO: change this to skip when water is used
        if not resp:
            # skip RESP charges calculation
            self.partial_charges = np.zeros(self.molecule.number_of_atoms())
            msg = 'RESP calculation disabled: All partial charges are set to zero.'
            self.ostream.print_info(msg)

        if self.partial_charges is None:
            if scf_results is None:
                # compute RESP charges from scratch
                if basis is None:
                    if self.rank == mpi_master():
                        basis = MolecularBasis.read(self.molecule,
                                                    '6-31G*',
                                                    ostream=None)
                    else:
                        basis = MolecularBasis()
                    basis = self.comm.bcast(basis, root=mpi_master())
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

            else:
                # compute RESP charges using the provided SCF result
                if basis is None:
                    error_msg = 'Basis is required for RESP charges.'
                    assert_msg_critical(False, error_msg)

                resp_drv = RespChargesDriver(self.comm, self.ostream)
                resp_drv.filename = self.molecule_name
                msg = 'Using provided SCF result for RESP charges.'
                self.ostream.print_info(msg)

                if self.resp_dict is not None:
                    resp_drv.update_settings(self.resp_dict)
                if resp_drv.equal_charges is None:
                    resp_drv.equal_charges = atomtypeidentifier.equivalent_charges

                self.partial_charges = resp_drv.compute(self.molecule, basis,
                                                        scf_results, 'resp')
                self.partial_charges = self.comm.bcast(self.partial_charges,
                                                       root=mpi_master())

        # self.partial charges should be summing to a whole number
        # Round to the itp writting precision (6 decimal places)
        self.partial_charges = np.round(self.partial_charges, 6)

        # Removing the tail of the partial charges to ensure the sum is a whole number
        excess_charge = (sum(self.partial_charges) -
                         round(sum(self.partial_charges)))

        max_charge_index = np.argmax(np.abs(self.partial_charges))
        self.partial_charges[max_charge_index] -= excess_charge

        if abs(excess_charge) > 1.0e-8:
            msg = 'Sum of partial charges is not a whole number.'
            self.ostream.print_info(msg)
            msg = f'Compensating by removing {excess_charge:.3e}'
            msg += ' from the largest charge.'
            self.ostream.print_info(msg)
            self.ostream.print_blank()

        # preparing atomtypes and atoms

        assert_msg_critical(
            len(self.atom_types) == n_atoms,
            'MMForceFieldGenerator: inconsistent atom_types')

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

        use_gaff = False
        use_uff = False
        use_tm = False
        use_water_model = False

        for at in self.unique_atom_types:
            atom_type_found = False
            
            # Auxilary variable for finding parameters in UFF
            element = ''
            for c in at:
                if not c.isdigit():
                    element += c
                else:
                    break

            if use_xml:
                for atom_type_data in ff_data_dict['atom_types']:
                    # Note: need strip() for converting e.g. 'c ' to 'c'
                    if atom_type_data['class'] == at.strip():
                        sigma = float(atom_type_data['sigma'])
                        epsilon = float(atom_type_data['epsilon'])
                        comment = 'GAFF'
                        atom_type_found = True
                        use_gaff = True
                        break
            
            else:
                for line in ff_data_lines:
                    if line.startswith(f'  {at}     '):
                        atom_ff = line[5:].strip().split()
                        sigma = float(atom_ff[0]) * 2**(-1 / 6) * 2 / 10
                        epsilon = float(atom_ff[1]) * 4.184
                        comment = 'GAFF'
                        atom_type_found = True
                        use_gaff = True
                        break

            if not atom_type_found:
                if at in ['ow','hw']:
                    assert_msg_critical(water_model is not None, 'MMForceFieldGenerator: water model not specified.')
                    assert_msg_critical(water_model in self.water_parameters, 
                        f"Error: '{water_model}' is not available. Available models are: {list(self.water_parameters.keys())}")
                    
                    sigma = self.water_parameters[water_model][at]['sigma']
                    epsilon = self.water_parameters[water_model][at]['epsilon']

                    water_bonds = self.water_parameters[water_model]['bonds']
                    water_angles = self.water_parameters[water_model]['angles']
                    self.partial_charges = [self.water_parameters[water_model][a]['charge'] for a in self.atom_types]
                    atom_type_found = True
                    use_water_model = True
                    self.eq_param = False
                    comment = water_model

                elif element in self.tm_parameters:
                    tmmsg = f'MMForceFieldGenerator: atom type {at} is not in GAFF.'
                    tmmsg += ' Taking TM parameters sigma and epsilon from vlx library.' ##TODO: rephrase
                    self.ostream.print_info(tmmsg)
                    sigma = self.tm_parameters[element]['sigma']
                    epsilon = self.tm_parameters[element]['epsilon']
                    comment = 'TM'
                    use_tm = True
                
                # Case for atoms in UFF but not in GAFF
                elif element in self.uff_parameters:
                    uffmsg = f'MMForceFieldGenerator: atom type {at} is not in GAFF.'
                    uffmsg += ' Taking sigma and epsilon from UFF.'
                    self.ostream.print_info(uffmsg)
                    sigma = self.uff_parameters[element]['sigma']
                    epsilon = self.uff_parameters[element]['epsilon']
                    comment = 'UFF'
                    use_uff = True

                else: 
                    assert_msg_critical(
                        False,
                        f'MMForceFieldGenerator: atom type {at} not found in GAFF or UFF.'
                    )

            atom_type_params[at] = {
                'sigma': sigma,
                'epsilon': epsilon,
                'comment': comment
            }

        if use_gaff:
            if gaff_version is not None:
                self.ostream.print_info(
                    f'Using GAFF (v{gaff_version}) parameters.')
            else:
                self.ostream.print_info('Using GAFF parameters.')
            gaff_ref = 'J. Wang, R. M. Wolf, J. W. Caldwell, P. A. Kollman,'
            gaff_ref += ' D. A. Case, J. Comput. Chem. 2004, 25, 1157-1174.'
            self.ostream.print_reference('Reference: ' + gaff_ref)
            self.ostream.print_blank()
            self.ostream.flush()

        if use_uff:
            self.ostream.print_info('Using UFF parameters.')
            uff_ref = 'A. K. Rappé, C. J. Casewit, K. S.  Colwell, W. A. Goddard III,'
            uff_ref += ' W. M. Skiff, J. Am. Chem. Soc. 1992, 114, 10024-10035.'
            self.ostream.print_reference('Reference: ' + uff_ref)
            self.ostream.print_blank()
            self.ostream.flush()

        if use_tm:
            self.ostream.print_info('Using TM parameters.')
            tm_ref = 'F. Šebesta, V. Sláma, J. Melcr, Z. Futera, and J. V. Burda.'
            tm_ref += 'J. Chem. Theory Comput. 2016 12 (8), 3681-3688.'
            self.ostream.print_reference('Reference: ' + tm_ref)
            self.ostream.print_blank()
            self.ostream.flush()

        if use_water_model:
            self.ostream.print_info(f'Using modified water model parameters for {water_model}.')
            wff_ref = 'T. Luchko, S. Gusarov, D. R. Roe, C. Simmerling, D. A. Case, J. Tuszynski,'
            wff_ref += 'A. Kovalenko. J. Chem. Theory Comput. 2010 6 (3), 607-624.'
            self.ostream.print_reference('Reference: ' + wff_ref)
            self.ostream.print_blank()
            self.ostream.flush()

        # Atoms analysis

        self.atoms = {}

        atom_names = self.get_atom_names()
        atom_masses = self.molecule.get_masses()
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
            target_bond_types = [
                # Note: need strip() for converting e.g. 'c ' to 'c'
                (at_1.strip(), at_2.strip()),
                (at_2.strip(), at_1.strip()),
            ]

            bond_found = False
            r, k_r, comment = None, None, None

            if use_xml:
                for bond_data in ff_data_dict['bonds']:
                    for target_bond in target_bond_types:
                        if target_bond == (bond_data['class1'], bond_data['class2']):
                            r = float(bond_data['length'])
                            k_r = float(bond_data['k'])
                            comment = '-'.join(target_bond)
                            bond_found = True
                            break
            
            elif use_water_model: ##TODO: Double check the UNITS
                r = water_bonds['equilibrium']
                k_r = water_bonds['force_constant']
                comment = 'ow-hw'
                bond_found = True

            else:
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
                    msg = f'Updated bond length {i + 1}-{j + 1} '
                    msg += f'({at_1}-{at_2}) to {r_eq:.3f} nm'
                    self.ostream.print_info(msg)
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
            target_angle_types = [
                # Note: need strip() for converting e.g. 'c ' to 'c'
                (at_1.strip(), at_2.strip(), at_3.strip()),
                (at_3.strip(), at_2.strip(), at_1.strip()),
            ]

            angle_found = False
            theta, k_theta, comment = None, None, None

            if use_xml:
                for angle_data in ff_data_dict['angles']:
                    for target_angle in target_angle_types:
                        if target_angle == (angle_data['class1'], angle_data['class2'], angle_data['class3']):
                            theta = float(angle_data['angle']) / np.pi * 180.0
                            k_theta = float(angle_data['k'])
                            comment = '-'.join(target_angle)
                            angle_found = True
                            break
            
            elif use_water_model:
                k_theta = water_angles['force_constant']
                theta = water_angles['equilibrium']
                comment = water_angles['comment']
                angle_found = True

            else:
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
                    msg = f'Updated bond angle {i + 1}-{j + 1}-{k + 1} '
                    msg += f'({at_1}-{at_2}-{at_3}) to {theta_eq:.3f} deg'
                    self.ostream.print_info(msg)
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
            target_dihedral_types = [
                (at_1.strip(), at_2.strip(), at_3.strip(), at_4.strip()),
                (at_4.strip(), at_3.strip(), at_2.strip(), at_1.strip()),
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
                    target_dihedral_types = [
                        ('', 'cp', 'cp', ''),
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
                    target_dihedral_types = [
                        ('', 'ca', 'na', ''),
                        ('', 'na', 'ca', ''),
                    ]
                    special_comment = ('(Guessed for rotatable ' +
                                       f'{at_1}-{at_2}-{at_3}-{at_4})')

            dihedral_found = False

            dihedral_ff_lines = []
            dihedral_matches = []

            dihedral_ff_data = []

            if use_xml:
                for dihedral_data in ff_data_dict['dihedrals']:
                    for target_dihedral in target_dihedral_types:
                        if target_dihedral == (dihedral_data['class1'],
                                               dihedral_data['class2'],
                                               dihedral_data['class3'],
                                               dihedral_data['class4']):
                            dihedral_ff_data.append(dict(dihedral_data))
                            dihedral_matches.append(self.get_dihedral_type_string(target_dihedral) + special_comment)
                            dihedral_found = True
                            break
            else:
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
                target_dihedral_types = [
                    ('', at_2.strip(), at_3.strip(), ''),
                    ('', at_3.strip(), at_2.strip(), ''),
                ]

                dihedral_ff_lines = []
                dihedral_matches = []

                dihedral_ff_data = []

                if use_xml:
                    for dihedral_data in ff_data_dict['dihedrals']:
                        for target_dihedral in target_dihedral_types:
                            if target_dihedral == (dihedral_data['class1'],
                                                   dihedral_data['class2'],
                                                   dihedral_data['class3'],
                                                   dihedral_data['class4']):
                                dihedral_ff_data.append(dict(dihedral_data))
                                dihedral_matches.append(self.get_dihedral_type_string(target_dihedral))
                                dihedral_found = True
                                break
                else:
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
                target_dihedral_types = []
                for p in patterns:
                    target_dih_str = p.pattern.replace(r'\A', '')
                    target_dih_tuple = tuple([
                        x.strip().replace('X', '')
                        for x in target_dih_str.split('-')
                    ])
                    target_dihedral_types.append(target_dih_tuple)

                dihedral_ff_lines = []
                dihedral_matches = []

                dihedral_ff_data = []

                if use_xml:
                    for dihedral_data in ff_data_dict['dihedrals']:
                        for target_dihedral in target_dihedral_types:
                            if target_dihedral == (dihedral_data['class1'],
                                                   dihedral_data['class2'],
                                                   dihedral_data['class3'],
                                                   dihedral_data['class4']):
                                dihedral_ff_data.append(dict(dihedral_data))
                                dihedral_matches.append(self.get_dihedral_type_string(target_dihedral) +
                                                        '(Guessed for ' +
                                                        f'{at_1}-{at_2}-{at_3}-{at_4})')
                                dihedral_found = True
                                break
                else:
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
                warnmsg = f'MMForceFieldGenerator: dihedral {at_1}-{at_2}-{at_3}-{at_4}'
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

            if use_xml:
                for dihedral_data, comment in zip(dihedral_ff_data, dihedral_matches):
                    dih_param_idx = 1
                    while f'periodicity{dih_param_idx}' in dihedral_data:
                        periodicity = int(dihedral_data[f'periodicity{dih_param_idx}'])
                        barrier = float(dihedral_data[f'k{dih_param_idx}'])
                        phase = float(dihedral_data[f'phase{dih_param_idx}']) / np.pi * 180.0

                        dihedral_barriers.append(barrier)
                        dihedral_phases.append(phase)
                        dihedral_periodicities.append(periodicity)
                        dihedral_comments.append(comment)

                        dih_param_idx += 1
            else:
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

        # Fetch the rotatable bonds from the molecule and the atom types involved
        rotatable_bonds_types = {}
        for i, j, k, l in self.dihedrals:
            # Ensure consistent ordering of bond indices
            bond_indices = (min(j, k), max(j, k))
            bond_types = (self.atom_types[bond_indices[0]], self.atom_types[bond_indices[1]])
            rotatable_bonds_types[bond_indices] = bond_types

        # Check if the rotatable bonds are indeed rotatable or not
        updated_rotatable_bonds = self.check_rotatable_bonds(rotatable_bonds_types)

        # Create a 1-indexed list of rotatable bonds without duplicates
        self.rotatable_bonds = [[bond[0] + 1, bond[1] + 1] for bond in updated_rotatable_bonds.keys()]

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
                target_dihedral_types = [
                    (at_2.strip(), at_3.strip(), at_4.strip(), at_1.strip()),
                    (at_2.strip(), at_1.strip(), at_4.strip(), at_3.strip()),
                    (at_2.strip(), at_4.strip(), at_1.strip(), at_3.strip()),
                    (at_2.strip(), at_3.strip(), at_1.strip(), at_4.strip()),
                    (at_2.strip(), at_4.strip(), at_3.strip(), at_1.strip()),
                    (at_2.strip(), at_1.strip(), at_3.strip(), at_4.strip()),
                ]

                dihedral_found = False
                barrier, phase, periodicity, comment = None, None, None, None

                if use_xml:
                    for dihedral_data in ff_data_dict['impropers']:
                        for target_dihedral in target_dihedral_types:
                            if target_dihedral == (dihedral_data['class1'],
                                                   dihedral_data['class2'],
                                                   dihedral_data['class3'],
                                                   dihedral_data['class4']):
                                periodicity = int(dihedral_data[f'periodicity1'])
                                barrier = float(dihedral_data[f'k1'])
                                phase = float(dihedral_data[f'phase1']) / np.pi * 180.0
                                comment = self.get_dihedral_type_string(target_dihedral)
                                dihedral_found = True
                                break
                else:
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
                    target_dihedral_types = [
                        (at_2.strip(), '', at_3.strip(), at_1.strip()),
                        (at_2.strip(), '', at_1.strip(), at_3.strip()),
                        (at_2.strip(), '', at_4.strip(), at_3.strip()),
                        (at_2.strip(), '', at_3.strip(), at_4.strip()),
                        (at_2.strip(), '', at_4.strip(), at_1.strip()),
                        (at_2.strip(), '', at_1.strip(), at_4.strip()),
                    ]

                    if use_xml:
                        for dihedral_data in ff_data_dict['impropers']:
                            for target_dihedral in target_dihedral_types:
                                if target_dihedral == (dihedral_data['class1'],
                                                       dihedral_data['class2'],
                                                       dihedral_data['class3'],
                                                       dihedral_data['class4']):
                                    periodicity = int(dihedral_data[f'periodicity1'])
                                    barrier = float(dihedral_data[f'k1'])
                                    phase = float(dihedral_data[f'phase1']) / np.pi * 180.0
                                    comment = self.get_dihedral_type_string(target_dihedral)
                                    dihedral_found = True
                                    break
                    else:
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
                    target_dihedral_types = [
                        (at_2.strip(), '', '', at_3.strip()),
                        (at_2.strip(), '', '', at_1.strip()),
                        (at_2.strip(), '', '', at_4.strip()),
                    ]

                    if use_xml:
                        for dihedral_data in ff_data_dict['impropers']:
                            for target_dihedral in target_dihedral_types:
                                if target_dihedral == (dihedral_data['class1'],
                                                       dihedral_data['class2'],
                                                       dihedral_data['class3'],
                                                       dihedral_data['class4']):
                                    periodicity = int(dihedral_data[f'periodicity1'])
                                    barrier = float(dihedral_data[f'k1'])
                                    phase = float(dihedral_data[f'phase1']) / np.pi * 180.0
                                    comment = self.get_dihedral_type_string(target_dihedral)
                                    dihedral_found = True
                                    break
                    else:
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
                    'MMForceFieldGenerator: invalid improper dihedral phase')

                assert_msg_critical(
                    periodicity == 2,
                    'MMForceFieldGenerator: invalid improper dihedral periodicity'
                )

                self.impropers[(i, j, k, l)] = {
                    'type': 'Fourier',
                    'barrier': barrier,
                    'phase': phase,
                    'periodicity': periodicity,
                    'comment': comment
                }

        self.ostream.flush()

    @staticmethod
    def get_dihedral_type_string(target_dihedral):
        """
        Gets type string of a dihedral.

        :param target_dihedral:
            A tuple of atom types.

        :return:
            The type string of a dihedral.
        """

        new_dih_types = []

        for at in target_dihedral:
            if at:
                new_dih_types.append(at)
            else:
                new_dih_types.append('X')

        return '-'.join(new_dih_types)

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

    def add_bond(self, bond):
        """
        Adds a bond to the topology.

        :param bond:
            The bond to be added. As a list of 1-based atom indices.
        :param force_constant:
            The force constant of the bond. Default is 250000.00 kJ/mol/nm^2.
        :param equilibrium:
            The equilibrium distance of the bond. If none it will be calculated.
        """

        # Extract indices from the list
        i, j = bond

        # Convert to zero-based indices
        i = i - 1
        j = j - 1

        # Update the connectivity matrix
        self.connectivity_matrix[i, j] = 1

        # Change the topology update flag to True
        self.topology_update_flag = True

        # Print the information
        msg = f'Added bond {i + 1}-{j + 1}'
        self.ostream.print_info(msg)
        self.ostream.flush()

        msg = "Re-run create_topology() to update the topology."
        self.ostream.print_info(msg)
        self.ostream.flush()
    
    def add_dihedral(self, dihedral, barrier=1, phase=0, periodicity=1):
        """
        Adds a dihedral to the an existing dihedral in the topology
        converting it in a multiple dihedral.

        :param dihedral:
            The dihedral to be added. As a list of 1-based atom indices.
        :param barrier:
            The barrier of the dihedral. Default is 1.00 kJ/mol.
        :param phase:
            The phase of the dihedral. Default is 0.00 degrees.
        :param periodicity:
            The periodicity of the dihedral. Default is 1.
        """

        # Extract indices from the list
        i, j, k, l = dihedral

        # Convert to zero-based indices
        i = i - 1
        j = j - 1
        k = k - 1
        l = l - 1

        # Exctract the original dihedral parameters
        original_dihedral = self.dihedrals[(i, j, k, l)]
        if original_dihedral['multiple']:
            # Values are already in lists
            barriers_list = original_dihedral['barrier']
            phases_list = original_dihedral['phase']
            periodicities_list = original_dihedral['periodicity']
            comments_list = original_dihedral['comment']
        else:
            barriers_list = [original_dihedral['barrier']]
            phases_list = [original_dihedral['phase']]
            periodicities_list = [original_dihedral['periodicity']]
            comments_list = [original_dihedral['comment']]

        # Update the dihedral's data lists
        barriers_list.append(barrier)
        phases_list.append(phase)
        # Negative periodicity implies multitermed dihedral
        periodicities_list.append(-periodicity)
        comments_list.append(f'Added {i + 1}-{j + 1}-{k + 1}-{l + 1} dihedral')

        # Update the dihedral in the topology
        self.dihedrals[(i, j, k, l)] = {
            'type': 'Fourier',
            'multiple': True,
            'barrier': barriers_list,
            'phase': phases_list,
            'periodicity': periodicities_list,
            'comment': comments_list
        }

        # Print the information
        msg = f'Added dihedral {i + 1}-{j + 1}-{k + 1}-{l + 1}'
        self.ostream.print_info(msg)
        self.ostream.flush()

    def get_dihedral_params(self, atom_indices_for_dihedral):
        """
        Gets dihedral parameters.

        :param atom_indices_for_dihedral:
            One-based atom indices for the dihedral.

        :return:
            The dihedral parameters in a dictionary.
        """

        assert_msg_critical(
            len(atom_indices_for_dihedral) == 4,
            'MMForceFieldGenerator.get_dihedral_params: ' +
            'Expecting a tuple of four atom indices')

        # convert 1-based indices to 0-based indices
        key = tuple([x - 1 for x in atom_indices_for_dihedral])

        return dict(self.dihedrals[key])

    def set_dihedral_params(self, atom_indices_for_dihedral, dihedral_params):
        """
        Sets dihedral parameters.

        :param atom_indices_for_dihedral:
            One-based atom indices for the dihedral.
        :param dihedral_params:
            The dihedral parameters in a dictionary.
        """

        assert_msg_critical(
            len(atom_indices_for_dihedral) == 4,
            'MMForceFieldGenerator.set_dihedral_params: ' +
            'Expecting a tuple of four atom indices')

        assert_msg_critical(
            isinstance(dihedral_params, dict),
            'MMForceFieldGenerator.set_dihedral_params: ' +
            'Expecting a dictionary of dihedral parameters')

        # convert 1-based indices to 0-based indices
        key = tuple([x - 1 for x in atom_indices_for_dihedral])

        self.dihedrals[key] = dict(dihedral_params)

    def check_rotatable_bonds(self, rotatable_bonds_types):
        """
        Checks the rotatable bonds in the molecule and the atom types involved.

        :param rotatable_bonds_types:
            The dictionary of possible rotatable bonds and the atom types involved.
            key: tuple of atom indices
            value: tuple of atom types
        """

        non_rotatable_bonds = set()
        for atom1 in ['c2', 'n2', 'cc', 'ce', 'nc', 'ne']:
            for atom2 in ['c2', 'n2', 'cd', 'cf', 'nd', 'nf']:
                non_rotatable_bonds.add((atom1, atom2))
                non_rotatable_bonds.add((atom2, atom1))
        for atom1 in ['c ']:
            for atom2 in ['n ', 'ns']:
                non_rotatable_bonds.add((atom1, atom2))
                non_rotatable_bonds.add((atom2, atom1))
        non_rotatable_bonds = list(non_rotatable_bonds)

        # Identify bonds to delete based on criteria
        bonds_to_delete = []

        for (i, j), bond in rotatable_bonds_types.items():
            # Check if the bond is non-rotatable due to atom types
            if bond in non_rotatable_bonds:
                bonds_to_delete.append((i, j))
                continue

            # Check if the bond is part of a ring
            if self.is_bond_in_ring(i, j):
                bonds_to_delete.append((i, j))
                continue

        # Remove identified non-rotatable bonds
        for key in bonds_to_delete:
            rotatable_bonds_types.pop(key, None)

        # Exclude bonds involving terminal atoms
        bonds_to_delete = []
        for (i, j), bond in rotatable_bonds_types.items():
            atom_i_connections = np.where(self.connectivity_matrix[i] == 1)[0]
            atom_j_connections = np.where(self.connectivity_matrix[j] == 1)[0]

            if len(atom_i_connections) == 1 or len(atom_j_connections) == 1:
                bonds_to_delete.append((i, j))

        # Remove terminal atom bonds
        for key in bonds_to_delete:
            rotatable_bonds_types.pop(key, None)

        return rotatable_bonds_types

    def is_bond_in_ring(self, atom_i, atom_j):
        """
        Determines if the bond between atom_i and atom_j is part of a ring.

        :param atom_i: 
            Index of the first atom in the bond.
        :param atom_j: 
            Index of the second atom in the bond.
        :return: 
            True if the bond is part of a ring, False otherwise.
        """

        atom_i_info = self.atom_info_dict[atom_i + 1]  # +1 because atom_info_dict keys start from 1
        atom_j_info = self.atom_info_dict[atom_j + 1]

        # Check if both atoms are in cycles
        if atom_i_info["CyclicStructure"] == "cycle" and atom_j_info["CyclicStructure"] == "cycle":
            # Find common cycles
            common_cycles = set(atom_i_info["CycleNumber"]).intersection(atom_j_info["CycleNumber"])
            if common_cycles:
                return True
        return False

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
                False, 'MMForceFieldGenerator.reparameterize: expecting Hessian')

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
            # Mute all drivers
            xtb_drv.ostream.mute()
            xtb_grad_drv.ostream.mute()
            xtb_opt_drv.ostream.mute()
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
                'MMForceFieldGenerator.reparameterize: invalid Hessian matrix')

        else:
            assert_msg_critical(
                False,
                'MMForceFieldGenerator.reparameterize: invalid Hessian option')

        bohr_to_nm = bohr_in_angstrom() * 0.1

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
            new_force_constant *= hartree_in_kjpermol() / (bohr_to_nm**2)

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
            new_force_constant *= hartree_in_kjpermol()

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
            mol_name = 'MOL'

        with open(top_fname, 'w') as f_top:

            # header

            f_top.write('; Generated by VeloxChem\n')

            # defaults

            if amber_ff is not None:
                assert_msg_critical(
                    amber_ff.startswith('amber'),
                    'MMForceFieldGenerator.write_top: Invalid amber force field name'
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

            ##TODO: change this to fit with the self.water_parameters -- e.g., if water_model in self.water_parameters
            if water_model is not None:
                # very rudimentary check for water model names
                assert_msg_critical(
                    water_model.startswith('tip') or
                    water_model.startswith('spc'),
                    'MMForceFieldGenerator.write_top: Invalid water model name')
                assert_msg_critical(
                    amber_ff is not None, 'MMForceFieldGenerator.write_top: ' +
                    'amber_ff is required for water_model')
                water_include = str(
                    PurePath(f'{amber_ff}.ff') / f'{water_model}.itp') ##TODO: add maybe changed water model or similar
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
            mol_name = 'MOL'

        res_name = mol_name

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
                line_str = '{:6} {:>5} {:6} {:>6} {:>6}'.format(
                    i + 1, atom['type'], 1, res_name, atom['name'])
                line_str += ' {:5} {:13.6f} {:13.5f}'.format(
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

            dih_fourier_lines = []

            for (i, j, k, l), dih in self.dihedrals.items():

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
                line_str = '{:6}{:7}{:7}{:7}'.format(i + 1, j + 1, k + 1, l + 1)
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

        # Create the root element of the XML file
        ForceField = ET.Element("ForceField")

        # AtomTypes section
        AtomTypes = ET.SubElement(ForceField, "AtomTypes")

        for i, atom in self.atoms.items():

            # Get the element of the atom
            element = ''
            for c in atom['name']:
                if not c.isdigit():
                    element += c
                else:
                    break

            attributes = {
                # Name is the atom type_molname
                "name": atom['name'] + f'_{mol_name}',
                "class": str(i + 1) + f'_{mol_name}',
                "element": element,
                "mass": str(atom['mass'])
            }
            ET.SubElement(AtomTypes, "Type", **attributes)

        # Residues section
        Residues = ET.SubElement(ForceField, "Residues")
        Residue = ET.SubElement(Residues, "Residue", name=mol_name)
        for atom_id, atom_data in self.atoms.items():
            ET.SubElement(Residue,
                          "Atom",
                          name=atom_data['name'],
                          type=atom_data['name'] + f'_{mol_name}',
                          charge=str(atom_data['charge']))
        for bond_id, bond_data in self.bonds.items():
            ET.SubElement(Residue,
                          "Bond",
                          atomName1=self.atoms[bond_id[0]]['name'],
                          atomName2=self.atoms[bond_id[1]]['name'])

        # Bonds section
        Bonds = ET.SubElement(ForceField, "HarmonicBondForce")
        for bond_id, bond_data in self.bonds.items():
            attributes = {
                "class1": str(bond_id[0] + 1) + f'_{mol_name}',
                "class2": str(bond_id[1] + 1) + f'_{mol_name}',
                "length": str(bond_data['equilibrium']),
                "k": str(bond_data['force_constant'])
            }
            ET.SubElement(Bonds, "Bond", **attributes)

        # Angles section
        Angles = ET.SubElement(ForceField, "HarmonicAngleForce")
        for angle_id, angle_data in self.angles.items():
            attributes = {
                "class1": str(angle_id[0] + 1) + f'_{mol_name}',
                "class2": str(angle_id[1] + 1) + f'_{mol_name}',
                "class3": str(angle_id[2] + 1) + f'_{mol_name}',
                "angle": str(angle_data['equilibrium'] * np.pi / 180),
                "k": str(angle_data['force_constant'])
            }
            ET.SubElement(Angles, "Angle", **attributes)

        # Periodic Dihedrals section
        # Proper dihedrals
        Dihedrals = ET.SubElement(ForceField, "PeriodicTorsionForce")
        for dihedral_id, dihedral_data in self.dihedrals.items():
            # Not multiple dihedrals has periodicity1, phase1, k1
            if not dihedral_data['multiple']:
                attributes = {
                    "class1": str(dihedral_id[0] + 1) + f'_{mol_name}',
                    "class2": str(dihedral_id[1] + 1) + f'_{mol_name}',
                    "class3": str(dihedral_id[2] + 1) + f'_{mol_name}',
                    "class4": str(dihedral_id[3] + 1) + f'_{mol_name}',
                    "periodicity1": str(dihedral_data['periodicity']),
                    "phase1": str(dihedral_data['phase'] * np.pi / 180),
                    "k1": str(dihedral_data['barrier'])
                }
                ET.SubElement(Dihedrals, "Proper", **attributes)

            # Multiple dihedrals have periodicityi, phasei, ki for i in range(len(periodicity))
            # Format: http://docs.openmm.org/7.6.0/userguide/application/05_creating_ffs.html#periodictorsionforce
            else:
                
                # One set of classes
                attributes = {
                    "class1": str(dihedral_id[0] + 1) + f'_{mol_name}',
                    "class2": str(dihedral_id[1] + 1) + f'_{mol_name}',
                    "class3": str(dihedral_id[2] + 1) + f'_{mol_name}',
                    "class4": str(dihedral_id[3] + 1) + f'_{mol_name}',
                }
                # Multiple sets of periodicity, phase, k
                for i in range(len(dihedral_data['periodicity'])):
                    attributes.update({
                        f"periodicity{i+1}": str(abs(dihedral_data['periodicity'][i])),
                        f"phase{i+1}": str(dihedral_data['phase'][i] * np.pi / 180),
                        f"k{i+1}": str(dihedral_data['barrier'][i])
                    })

                ET.SubElement(Dihedrals, "Proper", **attributes)

        # Improper dihedrals
        for improper_id, improper_data in self.impropers.items():

            # The order of the atoms is defined in the OpenMM documentation
            # http://docs.openmm.org/latest/userguide/application/06_creating_ffs.html

            attributes = {
                "class1": str(improper_id[1] + 1) + f'_{mol_name}',
                "class2": str(improper_id[0] + 1) + f'_{mol_name}',
                "class3": str(improper_id[2] + 1) + f'_{mol_name}',
                "class4": str(improper_id[3] + 1) + f'_{mol_name}',
                "periodicity1": str(improper_data['periodicity']),
                "phase1": str(improper_data['phase'] * np.pi / 180),
                "k1": str(improper_data['barrier'])
            }
            ET.SubElement(Dihedrals, "Improper", **attributes)

        # NonbondedForce section
        NonbondedForce = ET.SubElement(ForceField,
                                       "NonbondedForce",
                                       coulomb14scale=str(self.fudgeQQ),
                                       lj14scale=str(self.fudgeLJ))
        for atom_id, atom_data in self.atoms.items():
            attributes = {
                "type": atom_data['name'] + f'_{mol_name}',
                "charge": str(atom_data['charge']),
                "sigma": str(atom_data['sigma']),
                "epsilon": str(atom_data['epsilon'])
            }
            ET.SubElement(NonbondedForce, "Atom", **attributes)

        # Generate the tree and write to file
        # tree = ET.ElementTree(ForceField)
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
            mol_name = 'MOL'

        res_name = mol_name

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
                    line_str += f'{coords_in_nm[i][d]:{ndec + 5}.{ndec}f}'
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

        # PDB format from http://deposit.rcsb.org/adit/docs/pdb_atom_format.html

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
            f_pdb.write(
                f'TITLE     PDB file of {mol_name}, generated by VeloxChem\n')
            f_pdb.write('MODEL        1\n')

            # Atoms
            for i, (atom, element) in enumerate(
                    zip(self.atoms.values(), molecule_elements), 1):
                atom_name = atom['name']
                occupancy = 1.00
                temp_factor = 0.00
                element_symbol = element[:2].rjust(2)

                # Format string from https://cupnet.net/pdb-format/

                line_str = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   ".format(
                    'HETATM', i, atom_name[:4], '', mol_name[:3], 'A', 1, '')

                line_str += "{:8.3f}{:8.3f}{:8.3f}".format(
                    coords_in_angstrom[i - 1][0], coords_in_angstrom[i - 1][1],
                    coords_in_angstrom[i - 1][2])

                line_str += "{:6.2f}{:6.2f}          {:>2s}".format(
                    occupancy, temp_factor, element_symbol)

                f_pdb.write(line_str + '\n')

            # CONECT section in the PDB file stating the connectivity.
            # Required by OpenMM to correctly assign topology.bonds
            for (i, j) in self.bonds:
                f_pdb.write(f'CONECT{i + 1:>5}{j + 1:>5}\n')

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
            mol_name = 'MOL'

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

    def validate_force_field(self, i, verbose=True):
        """
        Validates force field by RMSD of dihedral potentials.

        :param i:
            The index of the target dihedral.

        :return:
            A dictionary containing the results of validation.
        """

        dih = self.target_dihedrals[i]

        dih_str = f'{dih[0] + 1}-{dih[1] + 1}-{dih[2] + 1}-{dih[3] + 1}'
        self.ostream.print_info(f'  Target dihedral angle: {dih_str}')
        self.ostream.print_blank()

        geom = self.scan_geometries[i]
        angles = self.scan_dih_angles[i]
        mm_scan = self.perform_mm_scan(dih, geom, angles, verbose=verbose)
        mm_scan = np.array(mm_scan) - min(mm_scan)

        qm_scan = np.array(self.scan_energies[i]) - min(self.scan_energies[i])
        qm_scan *= hartree_in_kjpermol()

        if verbose:
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

        self.fitting_summary = {
            'maximum_difference': np.max(np.abs(mm_scan - qm_scan)),
            'standard_deviation': np.std(mm_scan - qm_scan),
        }

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
                        verbose=True):
        """
        Performs MM scan of a specific dihedral.

        :param dihedral:
            The dihedral (list of four atom ids).
        :param geometries:
            The scanned geometries for this dihedral.
        :param angles:
            The scanned angles for this dihedral.
        """

        # select scan angles and geometries from QM data

        if verbose:
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
            pot_energy *= hartree_in_kjpermol()
            energies.append(pot_energy)

            if verbose:
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

    def print_validation_summary(self, fitted_dihedral_results):
        """
        Prints validation summary.

        :param validation_result:
            The dictionary containing the result of validation.
        """

        self.ostream.print_info('Summary of validation')
        self.ostream.print_info('---------------------')

        scan_diff = (fitted_dihedral_results['mm_scan_kJpermol'] -
                     fitted_dihedral_results['qm_scan_kJpermol'])
        max_diff = np.max(np.abs(scan_diff))
        std_diff = np.std(scan_diff)
        self.fitting_summary = {
            'maximum_difference': max_diff,
            'standard_deviation': std_diff,
        }

        self.ostream.print_info(f'Maximum difference: {max_diff:.3f} kJ/mol')
        self.ostream.print_info(f'Standard deviation: {std_diff:.3f} kJ/mol')
        self.ostream.print_blank()
        self.ostream.flush()

    def visualize(self, validation_result, show_diff=False):
        """
        Visualizes dihedral potential.

        :param validation_result:
            The dictionary containing the result of validation.
        """

        try:
            import matplotlib.pyplot as plt
            from scipy.interpolate import make_interp_spline
        except ImportError:
            raise ImportError('Unable to import Matplotlib. Please install ' +
                              'Matplotlib via \'conda install matplotlib\'')

        qm_scan_kJpermol = validation_result['qm_scan_kJpermol']
        mm_scan_kJpermol = validation_result['mm_scan_kJpermol']
        dihedral_angles = validation_result['dihedral_angles']
        dihedral_indices = validation_result['dihedral_indices']

        # Fit spline
        dihedrals_dense = np.linspace(min(dihedral_angles), max(dihedral_angles), 300)
        spl = make_interp_spline(dihedral_angles, qm_scan_kJpermol, k=3)
        qm_scan_kJpermol_spl = spl(dihedrals_dense)
        spl = make_interp_spline(dihedral_angles, mm_scan_kJpermol, k=3)
        mm_scan_kJpermol_spl = spl(dihedrals_dense)

        # Plot spline
        plt.plot(dihedrals_dense, qm_scan_kJpermol_spl, color='black', linewidth=4,  label='QM (spline)', alpha=0.7)
        plt.plot(dihedrals_dense, mm_scan_kJpermol_spl, color='darkcyan', linewidth=4, label='MM (spline)' , alpha=0.7)

        if show_diff:
            plt.plot(dihedrals_dense, qm_scan_kJpermol_spl - mm_scan_kJpermol_spl, color='orange',
                     linewidth=2, linestyle=':', label='diff (spline)', alpha=0.7)

        # Print the original points 
        plt.scatter(dihedral_angles, qm_scan_kJpermol, color='black', s=25, label='QM (points)')
        plt.scatter(dihedral_angles, mm_scan_kJpermol, color='darkcyan', s=25, label='MM (points)')

        # Legend center right outside the plot
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.xlabel('Dihedral angle {}-{}-{}-{} (deg)'.format(dihedral_indices[0] + 1,
                                                       dihedral_indices[1] + 1,
                                                       dihedral_indices[2] + 1,
                                                       dihedral_indices[3] + 1))
        plt.ylabel('Energy (kJ/mol)')
        plt.title('Dihedral potential for rotatable bond {}-{}'.format(dihedral_indices[1] + 1,
                                                                       dihedral_indices[2] + 1))
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
            'MMForceFieldGenerator.get_included_file: could not find ' +
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
