#
#                              VELOXCHEM
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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
from datetime import datetime
from collections import deque
import numpy as np
import time as tm
import math
import sys

from .oneeints import compute_nuclear_potential_integrals
from .oneeints import compute_electric_dipole_integrals
from .veloxchemlib import OverlapDriver, KineticEnergyDriver
from .veloxchemlib import FockDriver, T4CScreener
from .veloxchemlib import GridDriver, XCIntegrator
from .veloxchemlib import AODensityMatrix, Matrices
from .veloxchemlib import make_matrix
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat, mat_t
from .veloxchemlib import xcfun
from .profiler import Profiler
from .molecularbasis import MolecularBasis
from .molecularorbitals import MolecularOrbitals, molorb
from .sadguessdriver import SadGuessDriver
from .subcommunicators import SubCommunicators
from .firstorderprop import FirstOrderProperties
from .inputparser import (parse_input, print_keywords, print_attributes,
                          get_random_string_parallel)
from .dftutils import get_default_grid_level, print_libxc_reference
from .sanitychecks import molecule_sanity_check, dft_sanity_check
from .errorhandler import assert_msg_critical
from .checkpoint import create_hdf5, write_scf_results_to_hdf5


class ScfDriver:
    """
    Implements SCF method with C2-DIIS and two-level C2-DIIS convergence
    accelerators.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - den_guess: The initial density guess driver.
        - acc_type: The type of SCF convergence accelerator.
        - max_err_vecs: The maximum number of error vectors.
        - max_iter: The maximum number of SCF iterations.
        - first_step: The flag for first step in two-level C2-DIIS convergence
          acceleration.
        - conv_thresh: The SCF convergence threshold.
        - eri_thresh: The electron repulsion integrals screening threshold.
        - ovl_thresh: The atomic orbitals linear dependency threshold.
        - diis_thresh: The C2-DIIS switch on threshold.
        - iter_data: The dictionary of SCF iteration data (scf energy, scf
          energy change, gradient, density change, etc.).
        - is_converged: The flag for SCF convergence.
        - scf_energy: The SCF energy.
        - num_iter: The current number of SCF iterations.
        - fock_matrices: The list of stored Fock/Kohn-Sham matrices.
        - den_matrices: The list of stored density matrices.
        - density: The current density matrix.
        - mol_orbs: The current molecular orbitals.
        - nuc_energy: The nuclear repulsion energy of molecule.
        - comm: The MPI communicator.
        - rank: The rank of MPI process.
        - nodes: The number of MPI processes.
        - restart: The flag for restarting from checkpoint file.
        - checkpoint_file: The name of checkpoint file.
        - ref_mol_orbs: The reference molecular orbitals read from checkpoint
          file.
        - scf_type: The type of SCF calculation (restricted, unrestricted, or
          restricted_openshell).
        - dft: The flag for running DFT.
        - grid_level: The accuracy level of DFT grid.
        - xcfun: The XC functional.
        - molgrid: The molecular grid.
        - pe: The flag for running polarizable embedding calculation.
        - V_es: The polarizable embedding matrix.
        - pe_options: The dictionary with options for polarizable embedding.
        - pe_summary: The summary string for polarizable embedding.
        - use_split_comm: The flag for using split communicators.
        - split_comm_ratio: The list of ratios for split communicators.
        - dispersion: The flag for calculating D4 dispersion correction.
        - d4_energy: The D4 dispersion correction to energy.
        - electric_field: The static electric field.
        - ef_nuc_energy: The electric potential energy of the nuclei in the
          static electric field.
        - dipole_origin: The origin of the dipole operator.
        - timing: The flag for printing timing information.
        - profiling: The flag for printing profiling information.
        - memory_profiling: The flag for printing memory usage.
        - memory_tracing: The flag for tracing memory allocation using
        - program_end_time: The end time of the program.
        - filename: The filename.
    """

    def __init__(self, comm, ostream):
        """
        Initializes SCF driver to default setup (convergence threshold, initial
        guess, etc).
        """

        # scf accelerator
        self.acc_type = 'L2_DIIS'
        self.max_err_vecs = 10
        self.max_iter = 50
        self._first_step = False

        # thresholds
        self.conv_thresh = 1.0e-6
        self.ovl_thresh = 1.0e-6
        self.diis_thresh = 1000.0
        self.eri_thresh = 1.0e-12
        self.eri_thresh_tight = 1.0e-15

        # iterations data
        self._history = None
        self._iter_data = None
        self._is_converged = False
        self._scf_energy = 0.0
        self._num_iter = 0

        # DIIS data
        self._fock_matrices_alpha = deque()
        self._fock_matrices_beta = deque()
        self._fock_matrices_proj = deque()

        self._density_matrices_alpha = deque()
        self._density_matrices_beta = deque()

        # density matrix and molecular orbitals
        self._density = None
        self._molecular_orbitals = MolecularOrbitals()

        # nuclear repulsion energy
        self._nuc_energy = 0.0

        # mpi information
        self._comm = comm
        self._rank = self._comm.Get_rank()
        self._nodes = self._comm.Get_size()

        # output stream
        self._ostream = ostream

        # restart information
        self.restart = True
        self.checkpoint_file = None
        self._ref_mol_orbs = None

        # Maximum overlap constraint
        self._mom = None

        # closed shell?
        self._scf_type = 'restricted'

        # D4 dispersion correction
        self.dispersion = False
        self._d4_energy = 0.0

        # for open-shell system: unpaired electrons for initial guess
        self.guess_unpaired_electrons = ''

        # dft
        self.xcfun = None
        self.grid_level = None
        self._dft = False
        self._mol_grid = None

        # polarizable embedding
        self.potfile = None
        self.pe_options = {}
        self._pe = False
        self._V_es = None
        self._pe_summary = ''

        # split communicators
        self.use_split_comm = False
        self._split_comm_ratio = None

        # static electric field
        self.electric_field = None
        self._ef_nuc_energy = 0.0
        self._dipole_origin = None

        # scf tensors
        self._scf_tensors = None

        # scf properties
        self._scf_prop = None

        # timing and profiling
        self.timing = False
        self.profiling = False
        self.memory_profiling = False
        self.memory_tracing = False

        # verbosity of output (1-3)
        self.print_level = 2

        # program end time for graceful exit
        self.program_end_time = None

        # filename
        self.filename = None

        # input keywords
        self._input_keywords = {
            'scf': {
                'acc_type':
                    ('str_upper', 'type of SCF convergence accelerator'),
                'max_iter': ('int', 'maximum number of SCF iterations'),
                'conv_thresh': ('float', 'SCF convergence threshold'),
                'eri_thresh': ('float', 'ERI screening threshold'),
                'restart': ('bool', 'restart from checkpoint file'),
                'filename': ('str', 'base name of output files'),
                'checkpoint_file': ('str', 'name of checkpoint file'),
                'timing': ('bool', 'print timing information'),
                'profiling': ('bool', 'print profiling information'),
                'memory_profiling': ('bool', 'print memory usage'),
                'memory_tracing': ('bool', 'trace memory allocation'),
                'print_level': ('int', 'verbosity of output (1-3)'),
                'guess_unpaired_electrons':
                    ('str', 'unpaired electrons for initila guess'),
            },
            'method_settings': {
                'dispersion': ('bool', 'use D4 dispersion correction'),
                'xcfun': ('str_upper', 'exchange-correlation functional'),
                'grid_level': ('int', 'accuracy level of DFT grid (1-8)'),
                'potfile': ('str', 'potential file for polarizable embedding'),
                'electric_field': ('seq_fixed', 'static electric field'),
                # 'use_split_comm': ('bool', 'use split communicators'),
            },
        }

    @property
    def scf_type(self):
        """
        Returns the SCF type.
        """

        return self._scf_type

    @property
    def comm(self):
        """
        Returns the MPI communicator.
        """

        return self._comm

    @property
    def rank(self):
        """
        Returns the MPI rank.
        """

        return self._rank

    @property
    def nodes(self):
        """
        Returns the number of MPI processes.
        """

        return self._nodes

    @property
    def nnodes(self):
        """
        Returns the number of MPI processes.
        """

        return self._nodes

    @property
    def ostream(self):
        """
        Returns the output stream.
        """

        return self._ostream

    @property
    def num_iter(self):
        """
        Returns the current number of SCF iterations.
        """

        return self._num_iter

    @property
    def is_converged(self):
        """
        Returns whether SCF is converged.
        """

        return self._is_converged

    @property
    def scf_energy(self):
        """
        Returns SCF energy.
        """

        return self._scf_energy

    @property
    def density(self):
        """
        Returns the density matrix.
        """

        return self._density

    @property
    def molecular_orbitals(self):
        """
        Returns the molecular orbitals.
        """

        return self._molecular_orbitals

    @property
    def mol_orbs(self):
        """
        Returns the molecular orbitals (for backward compatibility).
        """

        return self._molecular_orbitals

    @property
    def scf_tensors(self):
        """
        Returns the SCF tensors.
        """

        return self._scf_tensors

    @property
    def history(self):
        """
        Returns the SCF history.
        """

        return self._history

    def print_keywords(self):
        """
        Prints input keywords in SCF driver.
        """

        print_keywords(self._input_keywords, self.ostream)

    def print_attributes(self):
        """
        Prints attributes in SCF driver.
        """

        print_attributes(self._input_keywords, self.ostream)

    def update_settings(self, scf_dict, method_dict=None):
        """
        Updates settings in SCF driver.

        :param scf_dict:
            The dictionary of scf input.
        :param method_dict:
            The dicitonary of method settings.
        """

        if method_dict is None:
            method_dict = {}

        scf_keywords = {
            key: val[0] for key, val in self._input_keywords['scf'].items()
        }

        parse_input(self, scf_keywords, scf_dict)

        if 'program_end_time' in scf_dict:
            self.program_end_time = scf_dict['program_end_time']
        if 'filename' in scf_dict:
            self.filename = scf_dict['filename']
            if 'checkpoint_file' not in scf_dict:
                self.checkpoint_file = f'{self.filename}.scf.h5'

        method_keywords = {
            key: val[0]
            for key, val in self._input_keywords['method_settings'].items()
        }

        parse_input(self, method_keywords, method_dict)

        dft_sanity_check(self, 'update_settings')

        self._pe_sanity_check(method_dict)

        if self.electric_field is not None:
            assert_msg_critical(
                len(self.electric_field) == 3,
                'SCF driver: Expecting 3 values in \'electric field\' input')
            assert_msg_critical(
                not self._pe,
                'SCF driver: \'electric field\' input is incompatible with ' +
                'polarizable embedding')
            # disable restart of calculation with static electric field since
            # checkpoint file does not contain information about the electric
            # field
            self.restart = False

    def _pe_sanity_check(self, method_dict=None):
        """
        Checks PE settings and updates relevant attributes.

        :param method_dict:
            The dicitonary of method settings.
        """

        if method_dict:
            if 'pe_options' in method_dict:
                self.pe_options = dict(method_dict['pe_options'])
            else:
                self.pe_options = {}

        if self.potfile:
            self.pe_options['potfile'] = self.potfile

        self._pe = ('potfile' in self.pe_options)

        if self._pe:
            potfile = None
            if self.rank == mpi_master():
                potfile = self.pe_options['potfile']
                if not Path(potfile).is_file():
                    potfile = str(
                        Path(self.filename).parent / Path(potfile).name)
            potfile = self.comm.bcast(potfile, root=mpi_master())
            self.pe_options['potfile'] = potfile

    def compute(self, molecule, ao_basis, min_basis=None):
        """
        Performs SCF calculation using molecular data.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        """

        profiler = Profiler()

        if min_basis is None:
            if self.rank == mpi_master():
                min_basis = MolecularBasis.read(molecule,
                                                'AO-START-GUESS',
                                                basis_path='.',
                                                ostream=None)
            else:
                min_basis = None
            min_basis = self.comm.bcast(min_basis, root=mpi_master())

        # check molecule
        molecule_sanity_check(molecule)

        # check dft setup
        dft_sanity_check(self, 'compute')

        # check pe setup
        self._pe_sanity_check()

        # check print level (verbosity of output)
        if self.print_level < 2:
            self.print_level = 1
        if self.print_level > 2:
            self.print_level = 3

        if self.restart:
            self.restart = self.validate_checkpoint(molecule.get_element_ids(),
                                                    ao_basis.get_label(),
                                                    self.scf_type)

        if self.restart:
            self.acc_type = 'DIIS'
            if self.rank == mpi_master():
                self._ref_mol_orbs = MolecularOrbitals.read_hdf5(
                    self.checkpoint_file)

        # nuclear repulsion energy
        self._nuc_energy = molecule.nuclear_repulsion_energy()

        if self.rank == mpi_master():
            self._print_header()
            valstr = 'Nuclear repulsion energy: {:.10f} a.u.'.format(
                self._nuc_energy)
            self.ostream.print_info(valstr)
            self.ostream.print_blank()

        # D4 dispersion correction
        if self.dispersion:
            if self.rank == mpi_master():
                disp = DispersionModel()
                xc_label = self.xcfun.get_func_label() if self._dft else 'HF'
                disp.compute(molecule, xc_label)
                self._d4_energy = disp.get_energy()
            else:
                self._d4_energy = 0.0
            self._d4_energy = self.comm.bcast(self._d4_energy,
                                              root=mpi_master())
        else:
            self._d4_energy = 0.0

        # generate integration grid
        if self._dft:
            print_libxc_reference(self.xcfun, self.ostream)

            grid_drv = GridDriver(self.comm)
            grid_level = (get_default_grid_level(self.xcfun)
                          if self.grid_level is None else self.grid_level)
            grid_drv.set_level(grid_level)

            grid_t0 = tm.time()
            self._mol_grid = grid_drv.generate(molecule)
            n_grid_points = self._mol_grid.number_of_points()
            self.ostream.print_info(
                'Molecular grid with {0:d} points generated in {1:.2f} sec.'.
                format(n_grid_points,
                       tm.time() - grid_t0))
            self.ostream.print_blank()

        # set up polarizable embedding
        if self._pe:
            # TODO: add PE info
            pot_info = 'Reading polarizable embedding potential: {}'.format(
                self.pe_options['potfile'])
            self.ostream.print_info(pot_info)
            self.ostream.print_blank()

        # C2-DIIS method
        if self.acc_type.upper() == 'DIIS':
            if self.rank == mpi_master():
                if self.restart:
                    den_mat = self.gen_initial_density_restart(molecule)
                else:
                    den_mat = self.gen_initial_density_sad(
                        molecule, ao_basis, min_basis)
            else:
                den_mat = None
            den_mat = self.comm.bcast(den_mat, root=mpi_master())
            self._comp_diis(molecule, ao_basis, min_basis, den_mat, profiler)

        # two level C2-DIIS method
        if self.acc_type.upper() == 'L2_DIIS':

            # first step
            self._first_step = True

            old_thresh = self.conv_thresh
            self.conv_thresh = 1.0e-3

            old_max_iter = self.max_iter
            self.max_iter = 5

            val_basis = ao_basis.reduce_to_valence_basis()
            if self.rank == mpi_master():
                den_mat = self.gen_initial_density_sad(molecule, val_basis,
                                                       min_basis)
            else:
                den_mat = None
            den_mat = self.comm.bcast(den_mat, root=mpi_master())
            self._comp_diis(molecule, val_basis, min_basis, den_mat, profiler)

            # second step
            self._first_step = False

            self.diis_thresh = 1000.0
            self.conv_thresh = old_thresh
            self.max_iter = old_max_iter

            if self.rank == mpi_master():
                den_mat = self.gen_initial_density_proj(
                    molecule, ao_basis, val_basis, self._molecular_orbitals)
            else:
                den_mat = None
            den_mat = self.comm.bcast(den_mat, root=mpi_master())
            self._comp_diis(molecule, ao_basis, val_basis, den_mat, profiler)

        self._fock_matrices_alpha.clear()
        self._fock_matrices_beta.clear()
        self._fock_matrices_proj.clear()

        self._density_matrices_alpha.clear()
        self._density_matrices_beta.clear()

        profiler.end(self.ostream, scf_flag=True)

        if not self.is_converged:
            self.ostream.print_header(
                '*** Warning: SCF is not converged!'.ljust(92))
            self.ostream.print_blank()
            self.ostream.flush()
            return

        if self.rank == mpi_master():
            self._print_scf_energy()

            s2 = self.compute_s2(molecule, self.scf_tensors)
            self._print_ground_state(molecule, s2)

            if self.print_level == 2:
                self.molecular_orbitals.print_orbitals(molecule, ao_basis, None,
                                                       self.ostream)
            if self.print_level == 3:
                self.molecular_orbitals.print_orbitals(
                    molecule, ao_basis,
                    (0, self.molecular_orbitals.number_of_mos()), self.ostream)

            self._scf_prop.print_properties(molecule)

            self.ostream.flush()

        return self.scf_tensors

    def gen_initial_density_sad(self, molecule, ao_basis, min_basis):
        """
        Computes initial AO density using superposition of atomic densities
        scheme.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis.
        :param min_basis:
            The minimal AO basis for generation of atomic densities.

        :return:
            The AO density matrix.
        """

        sad_drv = SadGuessDriver()

        return sad_drv.compute(molecule, min_basis, ao_basis, self.scf_type)

    def gen_initial_density_proj(self, molecule, ao_basis, valence_basis,
                                 valence_mo):
        """
        Computes initial AO density from molecular orbitals obtained by
        inserting molecular orbitals from valence basis into molecular
        orbitals in full AO basis.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis.
        :param valence_basis:
            The valence AO basis for generation of molecular orbitals.
        :param valence_mo:
            The molecular orbitals in valence AO basis.

        :return:
            The AO density matrix.
        """

        # nao_1 = valence_mo.number_of_aos()
        nmo_1 = valence_mo.number_of_mos()

        nao_2 = 0
        for ang in range(ao_basis.max_angular_momentum() + 1):
            nao_2 += ao_basis.number_of_basis_functions(ang) * (2 * ang + 1)

        mo_1a = valence_mo.alpha_to_numpy()
        mo_2a = np.zeros((nao_2, nmo_1))

        if self.scf_type == 'unrestricted':
            mo_1b = valence_mo.beta_to_numpy()
            mo_2b = np.zeros((nao_2, nmo_1))

        row_1 = 0
        row_2 = 0
        for ang in range(valence_basis.max_angular_momentum() + 1):
            for s in range(-ang, ang + 1):
                for a in range(molecule.number_of_atoms()):
                    nbf_1 = valence_basis.number_of_basis_functions([a], ang)
                    nbf_2 = ao_basis.number_of_basis_functions([a], ang)
                    if nbf_1 > 0:
                        mo_2a[row_2:row_2 + nbf_2, :] = mo_1a[row_1:row_1 +
                                                              nbf_1, :]
                        if self.scf_type == 'unrestricted':
                            mo_2b[row_2:row_2 + nbf_2, :] = mo_1b[row_1:row_1 +
                                                                  nbf_1, :]
                    row_1 += nbf_1
                    row_2 += nbf_2

        if self.scf_type == 'restricted':
            proj_mo = MolecularOrbitals(
                [mo_2a], [np.zeros(nmo_1)],
                [molecule.get_aufbau_alpha_occupation(nmo_1)],
                valence_mo.get_orbitals_type())

        elif self.scf_type == 'unrestricted':
            proj_mo = MolecularOrbitals(
                [mo_2a, mo_2b],
                [np.zeros(nmo_1), np.zeros(nmo_1)], [
                    molecule.get_aufbau_alpha_occupation(nmo_1),
                    molecule.get_aufbau_beta_occupation(nmo_1)
                ], valence_mo.get_orbitals_type())

        elif self.scf_type == 'restricted_openshell':
            proj_mo = MolecularOrbitals([mo_2a], [np.zeros(nmo_1)], [
                molecule.get_aufbau_alpha_occupation(nmo_1),
                molecule.get_aufbau_beta_occupation(nmo_1)
            ], valence_mo.get_orbitals_type())

        else:
            proj_mo = None
            assert_msg_critical(
                False, 'ScfDriver.gen_initial_density_proj: ' +
                'Invalid molecular orbitals type')

        return proj_mo.get_density(molecule, self.scf_type)

    def gen_initial_density_restart(self, molecule):
        """
        Reads initial molecular orbitals and AO density from checkpoint file.

        :param molecule:
            The molecule.

        :return:
            The AO density matrix.
        """

        self._molecular_orbitals = MolecularOrbitals.read_hdf5(
            self.checkpoint_file)
        den_mat = self._molecular_orbitals.get_density(molecule, self.scf_type)

        restart_text = 'Restarting from checkpoint file: '
        restart_text += self.checkpoint_file
        self.ostream.print_info(restart_text)
        self.ostream.print_blank()

        return den_mat

    def validate_checkpoint(self, nuclear_charges, basis_set, scf_type):
        """
        Validates the checkpoint file by checking nuclear charges and basis set.

        :param nuclear_charges:
            Numpy array of the nuclear charges.
        :param basis_set:
            Name of the AO basis.
        :param scf_type:
            The type of SCF calculation (restricted, unrestricted, or
            restricted_openshell).

        :return:
            Validity of the checkpoint file.
        """

        valid = False

        if self.rank == mpi_master():
            if (isinstance(self.checkpoint_file, str) and
                    Path(self.checkpoint_file).is_file()):
                valid = MolecularOrbitals.match_hdf5(self.checkpoint_file,
                                                     nuclear_charges, basis_set,
                                                     scf_type)

        valid = self.comm.bcast(valid, root=mpi_master())

        return valid

    def maximum_overlap(self, molecule, basis, orbitals, alpha_list, beta_list):
        """
        Constraint the SCF calculation to find orbitals that maximize overlap
        with a reference set.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param orbitals:
            The reference MolecularOrbital object.
        :param alpha_list:
            The list of alpha occupied orbitals.
        :param beta_list:
            The list of beta occupied orbitals.
        """

        C_start_a, C_start_b = None, None

        if self.rank == mpi_master():
            n_alpha = molecule.number_of_alpha_electrons()
            n_beta = molecule.number_of_beta_electrons()

            # Reorder alpha to match beta
            if self.scf_type == 'restricted_openshell':
                alpha_list = beta_list + [
                    x for x in alpha_list if x not in beta_list
                ]

            err_excitations = 'ScfDriver.maximum_overlap: '
            err_excitations += 'incorrect definition of occupation lists'
            assert_msg_critical(len(alpha_list) == n_alpha, err_excitations)
            assert_msg_critical(len(beta_list) == n_beta, err_excitations)
            if self.scf_type == 'restricted':
                assert_msg_critical(alpha_list == beta_list, err_excitations)

            n_mo = orbitals.number_of_mos()
            mo_a = orbitals.alpha_to_numpy()
            mo_b = orbitals.beta_to_numpy()

            C_a = mo_a[:, alpha_list]
            C_b = mo_b[:, beta_list]
            self._mom = (C_a, C_b)

            # Create guess orbitals
            virtual_alpha = [x for x in range(n_mo) if x not in alpha_list]
            C_start_a = mo_a[:, alpha_list + virtual_alpha]
            if self.scf_type == 'unrestricted':
                virtual_beta = [x for x in range(n_mo) if x not in beta_list]
                C_start_b = mo_b[:, beta_list + virtual_beta]

        if self.scf_type == 'unrestricted':
            self.set_start_orbitals(molecule, basis, (C_start_a, C_start_b))
        else:
            self.set_start_orbitals(molecule, basis, C_start_a)

    def set_start_orbitals(self, molecule, basis, array):
        """
        Creates checkpoint file from numpy array containing starting orbitals.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param array:
            The numpy array (or list/tuple of numpy arrays).
        """

        # create MolecularOrbitals object from numpy array

        if self.rank == mpi_master():
            assert_msg_critical(
                isinstance(array, (np.ndarray, tuple, list)),
                'ScfDriver.set_start_orbitals: invalid input for alpha ' +
                'orbitals')

            C_alpha, C_beta = None, None

            if isinstance(array, np.ndarray):
                C_alpha, C_beta = array, None

            elif isinstance(array, (tuple, list)):
                if len(array) == 1:
                    C_alpha, C_beta = array[0], None
                elif len(array) == 2:
                    C_alpha, C_beta = array[0], array[1]
                else:
                    assert_msg_critical(
                        False,
                        'ScfDriver.set_start_orbitals: expecting one or two ' +
                        'input orbitals')

            err_array = 'ScfDriver.set_start_orbitals: expecting numpy array'
            err_mo = 'ScfDriver.set_start_orbitals: inconsistent number of MOs'
            err_ao = 'ScfDriver.set_start_orbitals: inconsistent number of AOs'

            n_ao = basis.get_dimension_of_basis(molecule)

            assert_msg_critical(isinstance(C_alpha, np.ndarray), err_array)
            assert_msg_critical(n_ao == C_alpha.shape[0], err_ao)
            n_mo = C_alpha.shape[1]

            ene_a = np.zeros(n_mo)
            occ_a = molecule.get_aufbau_alpha_occupation(n_mo)

            if self.scf_type == 'restricted':
                self._molecular_orbitals = MolecularOrbitals([C_alpha], [ene_a],
                                                             [occ_a],
                                                             molorb.rest)

            elif self.scf_type == 'unrestricted':
                assert_msg_critical(isinstance(C_beta, np.ndarray), err_array)
                assert_msg_critical(n_ao == C_beta.shape[0], err_ao)
                assert_msg_critical(n_mo == C_beta.shape[1], err_mo)
                ene_b = np.zeros(n_mo)
                occ_b = molecule.get_aufbau_beta_occupation(n_mo)
                self._molecular_orbitals = MolecularOrbitals([C_alpha, C_beta],
                                                             [ene_a, ene_b],
                                                             [occ_a, occ_b],
                                                             molorb.unrest)

            elif self.scf_type == 'restricted_openshell':
                occ_b = molecule.get_aufbau_beta_occupation(n_mo)
                self._molecular_orbitals = MolecularOrbitals([C_alpha], [ene_a],
                                                             [occ_a, occ_b],
                                                             molorb.restopen)

        else:
            self._molecular_orbitals = MolecularOrbitals()

        # write checkpoint file and sychronize MPI processes

        self.restart = True
        if self.checkpoint_file is None:
            if self.filename is not None:
                base_fname = self.filename
            else:
                name_string = get_random_string_parallel(self.comm)
                base_fname = 'vlx_' + name_string
            self.checkpoint_file = f'{base_fname}.scf.h5'
        self.write_checkpoint(molecule.elem_ids_to_numpy(), basis.get_label())
        self.comm.barrier()

    def write_checkpoint(self, nuclear_charges, basis_set):
        """
        Writes molecular orbitals to checkpoint file.

        :param nuclear_charges:
            The nuclear charges.
        :param basis_set:
            Name of the basis set.
        """

        if self.rank == mpi_master():
            if self.checkpoint_file and isinstance(self.checkpoint_file, str):
                self.molecular_orbitals.write_hdf5(self.checkpoint_file,
                                                   nuclear_charges, basis_set)
                self.ostream.print_blank()
                checkpoint_text = 'Checkpoint written to file: '
                checkpoint_text += self.checkpoint_file
                self.ostream.print_info(checkpoint_text)

    def _comp_diis(self, molecule, ao_basis, min_basis, den_mat, profiler):
        """
        Performs SCF calculation with C2-DIIS acceleration.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param min_basis:
            The minimal AO basis set.
        :param profiler:
            The profiler.
        """

        self._scf_tensors = None

        self._scf_prop = FirstOrderProperties(self.comm, self.ostream)

        self._history = []

        if not self._first_step:
            profiler.begin({
                'timing': self.timing,
                'profiling': self.profiling,
                'memory_profiling': self.memory_profiling,
                'memory_tracing': self.memory_tracing,
            })

        diis_start_time = tm.time()

        self._fock_matrices_alpha.clear()
        self._fock_matrices_beta.clear()
        self._fock_matrices_proj.clear()

        self._density_matrices_alpha.clear()
        self._density_matrices_beta.clear()

        ovl_mat, kin_mat, npot_mat, dipole_mats = self._comp_one_ints(
            molecule, ao_basis)

        if self.rank == mpi_master() and self.electric_field is not None:
            dipole_ints = dipole_mats

        linear_dependency = False

        if self.rank == mpi_master():
            t0 = tm.time()

            S = ovl_mat
            eigvals, eigvecs = np.linalg.eigh(S)
            num_eigs = sum(eigvals > self.ovl_thresh)
            if num_eigs < eigvals.size:
                eigvals = eigvals[-num_eigs:]
                eigvecs = eigvecs[:, -num_eigs:]
            oao_mat = eigvecs * (1.0 / np.sqrt(eigvals))

            self.ostream.print_info('Orthogonalization matrix computed in' +
                                    ' {:.2f} sec.'.format(tm.time() - t0))
            self.ostream.print_blank()

            nrow = oao_mat.shape[0]
            ncol = oao_mat.shape[1]
            linear_dependency = (nrow != ncol)

            if linear_dependency:
                ndim = nrow - ncol
                self.ostream.print_info(
                    'Removed ' + str(ndim) + ' linearly dependent' +
                    ' vector{:s}.'.format('' if ndim == 1 else 's'))
                self.ostream.print_blank()
            self.ostream.flush()

        else:
            oao_mat = None

        linear_dependency = self.comm.bcast(linear_dependency,
                                            root=mpi_master())

        if (linear_dependency and self.eri_thresh > self.eri_thresh_tight):
            self.eri_thresh = self.eri_thresh_tight

            if self.rank == mpi_master():
                self.ostream.print_info('ERI screening threshold tightened to' +
                                        ' {:.1e}.'.format(self.eri_thresh))
                self.ostream.print_blank()

        self._density = tuple([x.copy() for x in den_mat])

        if self.use_split_comm:
            self.use_split_comm = ((self._dft or self._pe) and self.nodes >= 8)

        if self.use_split_comm and not self._first_step:
            screener = None
            if not self._first_step:
                valstr = 'ERI'
                if self._dft:
                    valstr += '/DFT'
                if self._pe:
                    valstr += '/PE'
                self.ostream.print_info(
                    'Using sub-communicators for {}.'.format(valstr))
        else:
            if self.rank == mpi_master():
                screener = T4CScreener()
                screener.partition(ao_basis, molecule, 'eri')
            else:
                screener = None
            screener = self.comm.bcast(screener, root=mpi_master())

        profiler.check_memory_usage('Initial guess')

        self._split_comm_ratio = None

        e_grad = None

        if self.rank == mpi_master():
            self._print_scf_title()

        for i in self._get_scf_range():

            # set the current number of SCF iterations
            # (note the extra SCF cycle when starting from scratch)
            if self.restart:
                self._num_iter = i + 1
            else:
                self._num_iter = i

            profiler.set_timing_key(f'Iteration {self._num_iter:d}')

            iter_start_time = tm.time()

            fock_mat, vxc_mat, e_pe, V_pe = self._comp_2e_fock(
                den_mat, molecule, ao_basis, screener, e_grad, profiler)

            profiler.start_timer('ErrVec')

            e_el = self._comp_energy(fock_mat, vxc_mat, e_pe, kin_mat, npot_mat,
                                     den_mat)

            self._comp_full_fock(fock_mat, vxc_mat, V_pe, kin_mat, npot_mat)

            if self.rank == mpi_master() and self.electric_field is not None:
                efpot = sum([
                    -1.0 * ef * mat
                    for ef, mat in zip(self.electric_field, dipole_ints)
                ])

                if self.scf_type == 'restricted':
                    e_el += 2.0 * np.trace(np.matmul(efpot, den_mat[0]))
                    fock_mat[0] += efpot
                else:
                    e_el += np.trace(np.matmul(efpot,
                                               (den_mat[0] + den_mat[1])))
                    fock_mat[0] += efpot
                    fock_mat[1] += efpot

                self._ef_nuc_energy = 0.0
                coords = molecule.get_coordinates_in_bohr()
                elem_ids = molecule.get_element_ids()
                for i in range(molecule.number_of_atoms()):
                    self._ef_nuc_energy -= np.dot(
                        elem_ids[i] * (coords[i] - self._dipole_origin),
                        self.electric_field)

            e_grad, max_grad = self._comp_gradient(fock_mat, ovl_mat, den_mat,
                                                   oao_mat)

            # compute density change and energy change

            diff_den = self._comp_density_change(den_mat, self._density)

            e_scf = (e_el + self._nuc_energy + self._d4_energy +
                     self._ef_nuc_energy)

            diff_e_scf = e_scf - self.scf_energy

            self._iter_data = {
                'energy': e_scf,
                'gradient_norm': e_grad,
                'max_gradient': max_grad,
                'diff_density': diff_den,
                'diff_energy': diff_e_scf,
            }

            self._history.append(self._iter_data)

            # update density and energy

            self._density = tuple([x.copy() for x in den_mat])

            self._scf_energy = e_scf

            profiler.stop_timer('ErrVec')
            profiler.check_memory_usage('Iteration {:d} Fock build'.format(
                self._num_iter))

            # print iteration and check convergence

            self._print_iter_data(i)

            self._check_convergence(molecule, ovl_mat)

            if self.is_converged:
                break

            # compute new Fock matrix, molecular orbitals and density

            profiler.start_timer('FockDiag')

            self._store_diis_data(fock_mat, den_mat, ovl_mat, e_grad)

            eff_fock_mat = self._get_effective_fock(fock_mat, ovl_mat, oao_mat)

            self._molecular_orbitals = self._gen_molecular_orbitals(
                molecule, eff_fock_mat, oao_mat)

            if self._mom is not None:
                self._apply_mom(molecule, ovl_mat)

            self._update_mol_orbs_phase()

            if self.rank == mpi_master():
                den_mat = self.molecular_orbitals.get_density(molecule)
            else:
                den_mat = None
            den_mat = self.comm.bcast(den_mat, root=mpi_master())

            profiler.stop_timer('FockDiag')

            profiler.check_memory_usage('Iteration {:d} Fock diag.'.format(
                self._num_iter))

            if not self._first_step:
                iter_in_hours = (tm.time() - iter_start_time) / 3600
                if self._need_graceful_exit(iter_in_hours):
                    self._graceful_exit(molecule, ao_basis)

        if not self._first_step:
            self.write_checkpoint(molecule.get_element_ids(),
                                  ao_basis.get_label())

        if (not self._first_step) and self.is_converged:

            if self.rank == mpi_master():
                S = ovl_mat

                C_alpha = self.molecular_orbitals.alpha_to_numpy()
                C_beta = self.molecular_orbitals.beta_to_numpy()

                E_alpha = self.molecular_orbitals.ea_to_numpy()
                E_beta = self.molecular_orbitals.eb_to_numpy()

                n_mo = C_alpha.shape[1]
                occ_alpha = molecule.get_aufbau_alpha_occupation(n_mo)
                occ_beta = molecule.get_aufbau_beta_occupation(n_mo)

                if self.scf_type == 'restricted':
                    D_alpha = self._density[0]
                    D_beta = self._density[0]
                    F_alpha = fock_mat[0]
                    F_beta = fock_mat[0]
                else:
                    D_alpha = self._density[0]
                    D_beta = self._density[1]
                    F_alpha = fock_mat[0]
                    F_beta = fock_mat[1]

                den_type = (denmat.rest
                            if self.scf_type == 'restricted' else denmat.unrest)
                self._density = AODensityMatrix(self._density, den_type)

                self._scf_tensors = {
                    # eri info
                    'eri_thresh': self.eri_thresh,
                    # scf info
                    'scf_type': self.scf_type,
                    'scf_energy': self.scf_energy,
                    'restart': self.restart,
                    # scf tensors
                    'S': S,
                    'C_alpha': C_alpha,
                    'C_beta': C_beta,
                    'E_alpha': E_alpha,
                    'E_beta': E_beta,
                    'occ_alpha': occ_alpha,
                    'occ_beta': occ_beta,
                    'D_alpha': D_alpha,
                    'D_beta': D_beta,
                    'F_alpha': F_alpha,
                    'F_beta': F_beta,
                }

                if self._dft:
                    # dft info
                    self._scf_tensors['xcfun'] = self.xcfun.get_func_label()
                    if self.grid_level is not None:
                        self._scf_tensors['grid_level'] = self.grid_level

                if self._pe:
                    # pe info
                    self._scf_tensors['potfile'] = self.potfile

            else:
                self._scf_tensors = None

            self._scf_prop.compute_scf_prop(molecule, ao_basis,
                                            self.scf_tensors)

            if self.rank == mpi_master():
                self._scf_tensors['dipole_moment'] = np.array(
                    self._scf_prop.get_property('dipole_moment'))

                self._write_final_hdf5(molecule, ao_basis)

        if self.rank == mpi_master():
            self._print_scf_finish(diis_start_time)

        profiler.check_memory_usage('End of SCF')

    def _need_graceful_exit(self, iter_in_hours):
        """
        Checks if a graceful exit is needed.

        :param iter_in_hours:
            The time spent in one iteration (in hours).

        :return:
            True if a graceful exit is needed, False otherwise.
        """

        if self.program_end_time is not None:
            remaining_hours = (self.program_end_time -
                               datetime.now()).total_seconds() / 3600
            # exit gracefully when the remaining time is not sufficient to
            # complete the next iteration (plus 25% to be on the safe side).
            if remaining_hours < iter_in_hours * 1.25:
                return True
        return False

    def _graceful_exit(self, molecule, basis):
        """
        Gracefully exits the program.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.

        :return:
            The return code.
        """

        self.ostream.print_blank()
        self.ostream.print_info('Preparing for a graceful termination...')
        self.ostream.flush()

        self.write_checkpoint(molecule.elem_ids_to_numpy(), basis.get_label())

        self.ostream.print_blank()
        self.ostream.print_info('...done.')
        self.ostream.print_blank()
        self.ostream.print_info('Exiting program.')
        self.ostream.print_blank()
        self.ostream.flush()

        sys.exit(0)

    def _comp_one_ints(self, molecule, basis):
        """
        Computes one-electron integrals (overlap, kinetic energy and nuclear
        potential) using molecular data.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.

        :return:
            The one-electron integrals.
        """

        if self.rank == mpi_master():

            t0 = tm.time()

            ovl_drv = OverlapDriver()
            ovl_mat = ovl_drv.compute(molecule, basis)
            ovl_mat = ovl_mat.full_matrix().to_numpy()

            t1 = tm.time()

            kin_drv = KineticEnergyDriver()
            kin_mat = kin_drv.compute(molecule, basis)
            kin_mat = kin_mat.full_matrix().to_numpy()

            t2 = tm.time()

            # TODO: parallelize npot_mat

            npot_mat = compute_nuclear_potential_integrals(molecule, basis)

            t3 = tm.time()

        else:

            ovl_mat = None
            kin_mat = None
            npot_mat = None

        ovl_mat = self.comm.bcast(ovl_mat, root=mpi_master())
        kin_mat = self.comm.bcast(kin_mat, root=mpi_master())
        npot_mat = self.comm.bcast(npot_mat, root=mpi_master())

        if self.electric_field is not None:
            if molecule.get_charge() != 0:
                coords = molecule.get_coordinates_in_bohr()
                nuclear_charges = molecule.elem_ids_to_numpy()
                self._dipole_origin = np.sum(coords.T * nuclear_charges,
                                             axis=1) / np.sum(nuclear_charges)
            else:
                self._dipole_origin = np.zeros(3)

            dipole_mats = compute_electric_dipole_integrals(
                molecule, basis, list(self._dipole_origin))
        else:
            dipole_mats = None

        t4 = tm.time()

        if self.rank == mpi_master() and self.print_level > 1:

            self.ostream.print_info('Overlap matrix computed in' +
                                    ' {:.2f} sec.'.format(t1 - t0))
            self.ostream.print_blank()

            self.ostream.print_info('Kinetic energy matrix computed in' +
                                    ' {:.2f} sec.'.format(t2 - t1))
            self.ostream.print_blank()

            self.ostream.print_info('Nuclear potential matrix computed in' +
                                    ' {:.2f} sec.'.format(t3 - t2))
            self.ostream.print_blank()

            if self.electric_field is not None:
                self.ostream.print_info('Electric dipole matrices computed in' +
                                        ' {:.2f} sec.'.format(t4 - t3))
                self.ostream.print_blank()

            self.ostream.flush()

        return ovl_mat, kin_mat, npot_mat, dipole_mats

    def _comp_npot_mat_split_comm(self, molecule, basis):
        """
        Computes one-electron nuclear potential integral on split
        communicators.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.

        :return:
            The one-electron nuclear potential matrix.
        """

        node_grps = [p for p in range(self.nodes)]
        subcomm = SubCommunicators(self.comm, node_grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm

        ave, res = divmod(molecule.number_of_atoms(), self.nodes)
        counts = [ave + 1 if p < res else ave for p in range(self.nodes)]

        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        charges = molecule.elem_ids_to_numpy()[start:end].astype(float)
        coords = np.vstack(
            (molecule.x_to_numpy()[start:end], molecule.y_to_numpy()[start:end],
             molecule.z_to_numpy()[start:end])).T

        npot_drv = NuclearPotentialIntegralsDriver(local_comm)
        npot_mat = npot_drv.compute(molecule, basis, charges, coords)

        if local_comm.Get_rank() == mpi_master():
            npot_mat.reduce_sum(cross_comm.Get_rank(), cross_comm.Get_size(),
                                cross_comm)

        return npot_mat

    def _comp_2e_fock(self,
                      den_mat,
                      molecule,
                      basis,
                      screener,
                      e_grad=None,
                      profiler=None):
        """
        Computes Fock/Kohn-Sham matrix (only 2e part).

        :param den_mat:
            The AO density matrix.
        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param screener:
            The screening container object.
        :param e_grad:
            The electronic gradient.
        :param profiler:
            The profiler.

        :return:
            The Fock matrix, AO Kohn-Sham (Vxc) matrix, etc.
        """

        if self.use_split_comm and not self._first_step:
            fock_mat, vxc_mat, e_pe, V_pe = self._comp_2e_fock_split_comm(
                den_mat, molecule, basis, screener, e_grad, profiler)

        else:
            fock_mat, vxc_mat, e_pe, V_pe = self._comp_2e_fock_single_comm(
                den_mat, molecule, basis, screener, e_grad, profiler)

        return fock_mat, vxc_mat, e_pe, V_pe

    def _comp_2e_fock_single_comm(self,
                                  den_mat,
                                  molecule,
                                  basis,
                                  screener,
                                  e_grad=None,
                                  profiler=None):
        """
        Computes Fock/Kohn-Sham matrix on single communicator.

        :param den_mat:
            The AO density matrix.
        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param screener:
            The screening container object.
        :param e_grad:
            The electronic gradient.
        :param profiler:
            The profiler.

        :return:
            The Fock matrix, AO Kohn-Sham (Vxc) matrix, etc.
        """

        if self.scf_type == 'restricted':
            if self.rank == mpi_master():
                den_mat_for_fock = make_matrix(basis, mat_t.symmetric)
                den_mat_for_fock.set_values(den_mat[0])
            else:
                den_mat_for_fock = None

            den_mat_for_fock = self.comm.bcast(den_mat_for_fock,
                                               root=mpi_master())

        else:
            if self.rank == mpi_master():
                # for now we calculate Ka, Kb and Jab separately for open-shell
                den_mat_for_Ka = make_matrix(basis, mat_t.symmetric)
                den_mat_for_Ka.set_values(den_mat[0])

                den_mat_for_Kb = make_matrix(basis, mat_t.symmetric)
                den_mat_for_Kb.set_values(den_mat[1])

                den_mat_for_Jab = make_matrix(basis, mat_t.symmetric)
                den_mat_for_Jab.set_values(den_mat[0] + den_mat[1])
            else:
                den_mat_for_Ka = None
                den_mat_for_Kb = None
                den_mat_for_Jab = None

            den_mat_for_fock = Matrices()

            dm_fock = self.comm.bcast(den_mat_for_Ka, root=mpi_master())
            den_mat_for_fock.add(dm_fock, '0')

            dm_fock = self.comm.bcast(den_mat_for_Kb, root=mpi_master())
            den_mat_for_fock.add(dm_fock, '1')

            dm_fock = self.comm.bcast(den_mat_for_Jab, root=mpi_master())
            den_mat_for_fock.add(dm_fock, '2')

            dm_fock = None

        if e_grad is None:
            thresh_int = int(-math.log10(self.eri_thresh))
        else:
            thresh_int = int(-math.log10(self._get_dyn_threshold(e_grad)))

        eri_t0 = tm.time()

        fock_drv = FockDriver()

        # determine fock_type and exchange_scaling_factor
        fock_type = '2jk'
        exchange_scaling_factor = 1.0
        if self._dft and not self._first_step:
            if self.xcfun.is_hybrid():
                fock_type = '2jkx'
                exchange_scaling_factor = self.xcfun.get_frac_exact_exchange()
            else:
                fock_type = 'j'
                exchange_scaling_factor = 0.0

        # further determine exchange_scaling_factor, erf_k_coef and omega
        need_omega = (self._dft and (not self._first_step) and
                      self.xcfun.is_range_separated())
        if need_omega:
            exchange_scaling_factor = (self.xcfun.get_rs_alpha() +
                                       self.xcfun.get_rs_beta())
            erf_k_coef = -self.xcfun.get_rs_beta()
            omega = self.xcfun.get_rs_omega()
        else:
            erf_k_coef, omega = None, None

        fock_mat = None

        if self.scf_type == 'restricted':
            # restricted SCF
            fock_mat = fock_drv.compute(screener, self.rank, self.nodes,
                                        den_mat_for_fock, fock_type,
                                        exchange_scaling_factor, 0.0,
                                        thresh_int)

            fock_mat_np = fock_mat.full_matrix().to_numpy()

            if fock_type == 'j':
                # for pure functional
                fock_mat_np *= 2.0

            if need_omega:
                # for range-separated functional
                fock_mat = fock_drv.compute(screener, self.rank, self.nodes,
                                            den_mat_for_fock, 'kx_rs',
                                            erf_k_coef, omega, thresh_int)

                fock_mat_np -= fock_mat.full_matrix().to_numpy()

            fock_mat_np = self.comm.reduce(fock_mat_np, root=mpi_master())

            if self.rank == mpi_master():
                # Note: make fock_mat a list
                fock_mat = [fock_mat_np]
            else:
                fock_mat = None

        else:
            # unrestricted SCF or restricted open-shell SCF
            if fock_type == 'j':
                # for pure functional
                # den_mat_for_fock.matrix('2') is D_total
                fock_mat = fock_drv.compute(screener, self.rank, self.nodes,
                                            den_mat_for_fock.matrix('2'), 'j',
                                            0.0, 0.0, thresh_int)

                J_ab_np = fock_mat.full_matrix().to_numpy()

                fock_mat_a_np = J_ab_np
                fock_mat_b_np = J_ab_np.copy()

            else:
                fock_mat = fock_drv.compute(screener, self.rank, self.nodes,
                                            den_mat_for_fock, ['kx', 'kx', 'j'],
                                            exchange_scaling_factor, 0.0,
                                            thresh_int)

                K_a_np = fock_mat.matrix('0').full_matrix().to_numpy()
                K_b_np = fock_mat.matrix('1').full_matrix().to_numpy()
                J_ab_np = fock_mat.matrix('2').full_matrix().to_numpy()

                fock_mat_a_np = J_ab_np - K_a_np
                fock_mat_b_np = J_ab_np - K_b_np

            if need_omega:
                # for range-separated functional
                den_mat_for_erf_k = Matrices()
                den_mat_for_erf_k.add(den_mat_for_fock.matrix('0'), '0')
                den_mat_for_erf_k.add(den_mat_for_fock.matrix('1'), '1')

                fock_mat = fock_drv.compute(screener, self.rank, self.nodes,
                                            den_mat_for_erf_k,
                                            ['kx_rs', 'kx_rs'], erf_k_coef,
                                            omega, thresh_int)

                fock_mat_a_np -= fock_mat.matrix('0').full_matrix().to_numpy()
                fock_mat_b_np -= fock_mat.matrix('1').full_matrix().to_numpy()

            fock_mat_a_np = self.comm.reduce(fock_mat_a_np, root=mpi_master())
            fock_mat_b_np = self.comm.reduce(fock_mat_b_np, root=mpi_master())

            if self.rank == mpi_master():
                # Note: make fock_mat a list
                fock_mat = [fock_mat_a_np, fock_mat_b_np]
            else:
                fock_mat = None

        if self.timing:
            profiler.add_timing_info('FockERI', tm.time() - eri_t0)
        vxc_t0 = tm.time()

        if self._dft and not self._first_step:
            if self.xcfun.get_func_type() in [xcfun.lda, xcfun.gga, xcfun.mgga]:
                xc_drv = XCIntegrator(self.comm)
                # Note: vxc_mat will remain distributed across MPI processes.
                # XC energy and Vxc matrix will be reduced in _comp_energy
                # and _comp_full_fock
                vxc_mat = xc_drv.integrate_vxc_fock(molecule, basis, den_mat,
                                                    self._mol_grid, self.xcfun)
            else:
                assert_msg_critical(
                    False, 'SCF driver: Unsupported XC functional type')
        else:
            vxc_mat = None

        if self.timing and self._dft:
            profiler.add_timing_info('FockXC', tm.time() - vxc_t0)
        pe_t0 = tm.time()

        if self._pe and not self._first_step:
            # TODO: add PE contribution and update pe_summary
            pass
        else:
            e_pe, V_pe = 0.0, None

        if self.timing and self._pe:
            profiler.add_timing_info('FockPE', tm.time() - pe_t0)

        return fock_mat, vxc_mat, e_pe, V_pe

    def _comp_2e_fock_split_comm(self,
                                 fock_mat,
                                 den_mat,
                                 molecule,
                                 basis,
                                 screening,
                                 e_grad=None,
                                 profiler=None):
        """
        Computes Fock/Kohn-Sham matrix on split communicators.

        :param fock_mat:
            The AO Fock matrix (only 2e-part).
        :param den_mat:
            The AO density matrix.
        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param screening:
            The screening container object.
        :param e_grad:
            The electronic gradient.
        :param profiler:
            The profiler.

        :return:
            The AO Kohn-Sham (Vxc) matrix.
        """

        return None

    def _comp_energy(self, fock_mat, vxc_mat, e_pe, kin_mat, npot_mat, den_mat):
        """
        Computes the sum of SCF energy components: electronic energy, kinetic
        energy, and nuclear potential energy.

        :param fock_mat:
            The Fock/Kohn-Sham matrix (only 2e-part).
        :param vxc_mat:
            The Vxc matrix.
        :param e_pe:
            The polarizable embedding energy.
        :param kin_mat:
            The kinetic energy matrix.
        :param npot_mat:
            The nuclear potential matrix.
        :param den_mat:
            The density matrix.

        :return:
            The sum of electronic energy, kinetic energy and nuclear potential
            energy.
        """

        xc_ene = 0.0
        if self._dft and not self._first_step:
            xc_ene = self.comm.reduce(vxc_mat.get_energy(), root=mpi_master())

        if self.rank == mpi_master():
            # electronic, kinetic, nuclear energy
            D = den_mat
            F = fock_mat
            T = kin_mat
            V = npot_mat
            if self.scf_type == 'restricted':
                e_ee = np.sum(D[0] * F[0])
                e_kin = 2.0 * np.sum(D[0] * T)
                e_en = 2.0 * np.sum(D[0] * V)
            else:
                e_ee = 0.5 * (np.sum(D[0] * F[0]) + np.sum(D[1] * F[1]))
                e_kin = np.sum((D[0] + D[1]) * T)
                e_en = np.sum((D[0] + D[1]) * V)
            if self._dft and not self._first_step:
                e_ee += xc_ene
            if self._pe and not self._first_step:
                e_ee += e_pe
            e_sum = e_ee + e_kin + e_en
        else:
            e_sum = 0.0
        e_sum = self.comm.bcast(e_sum, root=mpi_master())

        return e_sum

    def _comp_full_fock(self, fock_mat, vxc_mat, pe_mat, kin_mat, npot_mat):
        """
        Computes full Fock/Kohn-Sham matrix by adding to 2e-part of
        Fock/Kohn-Sham matrix the kinetic energy and nuclear potential
        matrices.

        :param fock_mat:
            The Fock/Kohn-Sham matrix (2e-part).
        :param vxc_mat:
            The Vxc matrix.
        :param pe_mat:
            The polarizable embedding matrix.
        :param kin_mat:
            The kinetic energy matrix.
        :param npot_mat:
            The nuclear potential matrix.
        """

        np_xcmat_a, np_xcmat_b = None, None
        if self._dft and not self._first_step:
            np_xcmat_a = self.comm.reduce(vxc_mat.alpha_to_numpy(),
                                          root=mpi_master())
            if self.scf_type != 'restricted':
                np_xcmat_b = self.comm.reduce(vxc_mat.beta_to_numpy(),
                                              root=mpi_master())

        if self.rank == mpi_master():
            T = kin_mat
            V = npot_mat
            fock_mat[0] += (T + V)
            if self.scf_type != 'restricted':
                fock_mat[1] += (T + V)

            if self._dft and not self._first_step:
                fock_mat[0] += np_xcmat_a
                if self.scf_type != 'restricted':
                    fock_mat[1] += np_xcmat_b

            # TODO: add PE contribution to Fock
            if self._pe and not self._first_step:
                pass

    def _comp_gradient(self, fock_mat, ovl_mat, den_mat, oao_mat):
        """
        Computes electronic gradient using Fock/Kohn-Sham matrix.

        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param ovl_mat:
            The overlap matrix.
        :param den_mat:
            The density matrix.
        :param oao_mat:
            The orthogonalization matrix.

        :return:
            The electronic gradient.
        """

        return 0.0

    def _comp_density_change(self, den_mat, old_den_mat):
        """
        Computes norm of density change between two density matrices.

        :param den_mat:
            The current density matrix.
        :param old_den_mat:
            The previous density matrix.

        :return:
            The norm of change between two density matrices.
        """

        return 0.0

    def _store_diis_data(self, fock_mat, den_mat, ovl_mat, e_grad):
        """
        Stores Fock/Kohn-Sham and density matrices for current iteration.

        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param den_mat:
            The density matrix.
        :param ovl_mat:
            The overlap matrix (used in ROSCF).
        :param e_grad:
            The electronic gradient.
        """

        return

    def _get_effective_fock(self, fock_mat, ovl_mat, oao_mat):
        """
        Computes effective Fock/Kohn-Sham matrix in OAO basis by applying
        Lowdin or canonical orthogonalization to AO Fock/Kohn-Sham matrix.

        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param ovl_mat:
            The overlap matrix.
        :param oao_mat:
            The orthogonalization matrix.

        :return:
            The effective Fock/Kohn-Sham matrix.
        """

        return None

    def _gen_molecular_orbitals(self, molecule, fock_mat, oao_mat):
        """
        Generates molecular orbital by diagonalizing Fock/Kohn-Sham matrix.

        :param molecule:
            The molecule.
        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param oao_mat:
            The orthogonalization matrix.

        :return:
            The molecular orbitals.
        """

        return MolecularOrbitals()

    def _apply_mom(self, molecule, ovl_mat):
        """
        Apply the maximum overlap constraint.

        :param molecule:
            The molecule.
        :param ovl_mat:
            The overlap matrix..
        """

        if self.rank == mpi_master():
            smat = ovl_mat

            mo_a = self.molecular_orbitals.alpha_to_numpy()
            ea = self.molecular_orbitals.ea_to_numpy()
            occ_a = self.molecular_orbitals.occa_to_numpy()
            n_alpha = molecule.number_of_alpha_electrons()

            ovl = np.linalg.multi_dot([self._mom[0].T, smat, mo_a])
            argsort = np.argsort(np.sum(np.abs(ovl), 0))[::-1]
            # restore energy ordering
            argsort[:n_alpha] = np.sort(argsort[:n_alpha])
            argsort[n_alpha:] = np.sort(argsort[n_alpha:])
            mo_a = mo_a[:, argsort]
            ea = ea[argsort]

            if self.scf_type == 'restricted':
                self._molecular_orbitals = MolecularOrbitals([mo_a], [ea],
                                                             [occ_a],
                                                             molorb.rest)

            else:
                n_beta = molecule.number_of_beta_electrons()
                occ_b = self.molecular_orbitals.occb_to_numpy()

                if self.scf_type == 'unrestricted':
                    mo_b = self.molecular_orbitals.beta_to_numpy()
                    eb = self.molecular_orbitals.eb_to_numpy()
                elif self.scf_type == 'restricted_openshell':
                    # For ROHF, the beta orbitals have to be a subset of the alpha
                    mo_b = mo_a[:, :n_alpha]

                ovl = np.linalg.multi_dot([self._mom[1].T, smat, mo_b])
                argsort_b = np.argsort(np.sum(np.abs(ovl), 0))[::-1]
                # restore energy ordering
                argsort_b[:n_beta] = np.sort(argsort_b[:n_beta])
                argsort_b[n_beta:] = np.sort(argsort_b[n_beta:])

                if self.scf_type == 'unrestricted':
                    mo_b = mo_b[:, argsort_b]
                    eb = eb[argsort_b]
                    self._molecular_orbitals = MolecularOrbitals([mo_a, mo_b],
                                                                 [ea, eb],
                                                                 [occ_a, occ_b],
                                                                 molorb.unrest)
                elif self.scf_type == 'restricted_openshell':
                    mo_a[:, :n_alpha] = mo_a[:, argsort_b]
                    ea[:n_alpha] = ea[argsort_b]
                    self._molecular_orbitals = MolecularOrbitals(
                        [mo_a], [ea], [occ_a, occ_b], molorb.restopen)

    def _update_mol_orbs_phase(self):
        """
        Updates phase of molecular orbitals.
        """

        if self.rank == mpi_master():
            if self._ref_mol_orbs is None:
                return

            ref_mo_a = self._ref_mol_orbs.alpha_to_numpy()
            mo_a = self.molecular_orbitals.alpha_to_numpy()
            e_a = self.molecular_orbitals.ea_to_numpy()
            occ_a = self.molecular_orbitals.occa_to_numpy()

            for col in range(mo_a.shape[1]):
                if np.dot(mo_a[:, col], ref_mo_a[:, col]) < 0.0:
                    mo_a[:, col] *= -1.0

            if self.molecular_orbitals.get_orbitals_type() == molorb.rest:
                self._molecular_orbitals = MolecularOrbitals([mo_a], [e_a],
                                                             [occ_a],
                                                             molorb.rest)

            elif self.molecular_orbitals.get_orbitals_type() == molorb.unrest:
                ref_mo_b = self._ref_mol_orbs.beta_to_numpy()
                mo_b = self.molecular_orbitals.beta_to_numpy()
                e_b = self.molecular_orbitals.eb_to_numpy()
                occ_b = self.molecular_orbitals.occb_to_numpy()

                for col in range(mo_b.shape[1]):
                    if np.dot(mo_b[:, col], ref_mo_b[:, col]) < 0.0:
                        mo_b[:, col] *= -1.0

                self._molecular_orbitals = MolecularOrbitals([mo_a, mo_b],
                                                             [e_a, e_b],
                                                             [occ_a, occ_b],
                                                             molorb.unrest)

            elif self.molecular_orbitals.get_orbitals_type() == molorb.restopen:
                occ_b = self.molecular_orbitals.occb_to_numpy()
                self._molecular_orbitals = MolecularOrbitals([mo_a], [e_a],
                                                             [occ_a, occ_b],
                                                             molorb.restopen)

    def _get_dyn_threshold(self, e_grad):
        """
        Computes screening threshold for electron repulsion integrals based on
        value of electronic gradient.

        :param e_grad:
            The electronic gradient.

        :return:
            The screening threshold.
        """

        if e_grad < 1.0e-6:
            return self.eri_thresh

        nteri = math.pow(10, math.floor(math.log10(e_grad)))

        nteri = 1.0e-10 * nteri

        if nteri > 1.0e-10:
            return 1.0e-10

        if nteri < self.eri_thresh:
            return self.eri_thresh

        return nteri

    def _check_convergence(self, molecule, ovl_mat):
        """
        Sets SCF convergence flag by checking if convergence condition for
        electronic gradient is fullfiled.

        :param molecule:
            The molecule.
        :param ovl_mat:
            The overlap matrix.
        """

        self._is_converged = False

        if self._num_iter > 0:

            e_grad = self._iter_data['gradient_norm']

            if e_grad < self.conv_thresh:
                if self.restart:
                    # Note: when restarting from checkpoint, double check that the
                    # number of electrons are reasonable
                    nalpha = molecule.number_of_alpha_electrons()
                    nbeta = molecule.number_of_beta_electrons()
                    calc_nelec = self._comp_number_of_electrons(ovl_mat)
                    if (abs(calc_nelec[0] - nalpha) < 1.0e-3 and
                            abs(calc_nelec[1] - nbeta) < 1.0e-3):
                        self._is_converged = True
                else:
                    self._is_converged = True

    def _comp_number_of_electrons(self, ovl_mat):
        """
        Computes number of alpha and beta electrons from density matrices and
        overlap matrix.

        :param ovl_mat:
            The overlap matrix.

        :return:
            The number of alpha and beta electrons.
        """

        if self.rank == mpi_master():
            if self.scf_type == 'restricted':
                D_alpha = self._density[0]
                D_beta = self._density[0]
            else:
                D_alpha = self._density[0]
                D_beta = self._density[1]
            S = ovl_mat
            calc_nelec = (np.sum(D_alpha * S), np.sum(D_beta * S))
        else:
            calc_nelec = None
        calc_nelec = self.comm.bcast(calc_nelec, root=mpi_master())

        return calc_nelec

    def _get_scf_range(self):
        """
        Creates range of SCF iterations from maximum number of SCF iterations.

        :return:
            The range of SCF iterations.
        """

        # set the maximum number of SCF iterations
        # (note the extra SCF cycle when starting from scratch)
        if self.restart:
            return range(self.max_iter)
        else:
            return range(self.max_iter + 1)

    def _print_scf_energy(self):
        """
        Prints SCF energy information to output stream.
        """

        valstr = self.get_scf_type_str() + ':'
        self.ostream.print_header(valstr.ljust(92))
        self.ostream.print_header(('-' * len(valstr)).ljust(92))
        self._print_energy_components()

        if self._pe:
            self.ostream.print_blank()
            for line in self._pe_summary.splitlines():
                self.ostream.print_header(line.ljust(92))
            self.ostream.flush()

    def _print_header(self):
        """
        Prints SCF calculation setup details to output stream,
        """

        self.ostream.print_blank()
        self.ostream.print_header('Self Consistent Field Driver Setup')
        self.ostream.print_header(36 * '=')
        self.ostream.print_blank()

        str_width = 84
        cur_str = 'Wave Function Model             : ' + self.get_scf_type_str()
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Initial Guess Model             : ' + self._get_guess_type()
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'Convergence Accelerator         : ' + self._get_acc_type()
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Max. Number of Iterations       : ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Max. Number of Error Vectors    : ' + str(self.max_err_vecs)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold           : {:.1e}'.format(
            self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'ERI Screening Threshold         : {:.1e}'.format(
            self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Linear Dependence Threshold     : {:.1e}'.format(
            self.ovl_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        if self._dft:
            cur_str = 'Exchange-Correlation Functional : '
            cur_str += self.xcfun.get_func_label().upper()
            self.ostream.print_header(cur_str.ljust(str_width))
            grid_level = (get_default_grid_level(self.xcfun)
                          if self.grid_level is None else self.grid_level)
            cur_str = 'Molecular Grid Level            : ' + str(grid_level)
            self.ostream.print_header(cur_str.ljust(str_width))

        if self.electric_field is not None:
            cur_str = 'Static Electric Field           : '
            cur_str += str(self.electric_field)
            self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()

    def _print_scf_title(self):
        """
        Prints SCF cycles header to output stream.
        """

        if self._first_step:
            self.ostream.print_info('Starting Reduced Basis SCF calculation...')

        else:
            self.ostream.print_blank()
            if self._dft:
                valstr = '{} | {} | {} | {} | {} | {}'.format(
                    'Iter.', '   Kohn-Sham Energy', 'Energy Change',
                    'Gradient Norm', 'Max. Gradient', 'Density Change')
                self.ostream.print_header(valstr)
            else:
                valstr = '{} | {} | {} | {} | {} | {}'.format(
                    'Iter.', 'Hartree-Fock Energy', 'Energy Change',
                    'Gradient Norm', 'Max. Gradient', 'Density Change')
                self.ostream.print_header(valstr)
            self.ostream.print_header(92 * '-')

    def _print_scf_finish(self, start_time):
        """
        Prints SCF calculation finish message to output stream,

        :param start_time:
            The start time of SCF calculation.
        """

        if self._first_step:
            valstr = '...done. SCF energy in reduced basis set: '
            valstr += '{:.12f}'.format(self._scf_energy)
            valstr += ' a.u. Time: '
            valstr += '{:.2f}'.format(tm.time() - start_time) + ' sec.'
            self.ostream.print_info(valstr)
            self.ostream.print_blank()

        else:
            valstr = '*** SCF '
            if self.is_converged:
                valstr += 'converged in '
            else:
                valstr += 'NOT converged in '
            valstr += str(self._num_iter)
            valstr += ' iterations. Time: '
            valstr += '{:.2f}'.format(tm.time() - start_time) + ' sec.'
            self.ostream.print_blank()
            self.ostream.print_header(valstr.ljust(92))
            self.ostream.print_blank()

        self.ostream.flush()

    def _print_iter_data(self, i):
        """
        Prints SCF iteration data to output stream,

        :param i:
            The current SCF iteration.
        """

        if self.rank == mpi_master():
            # no output for first step in two level DIIS
            if self._first_step:
                return

            # DIIS or second step in two level DIIS
            if self._num_iter > 0:

                if self._iter_data:
                    te = self._iter_data['energy']
                    diff_te = self._iter_data['diff_energy']
                    e_grad = self._iter_data['gradient_norm']
                    max_grad = self._iter_data['max_gradient']
                    diff_den = self._iter_data['diff_density']

                if self._num_iter == 1:
                    diff_te = 0.0
                    diff_den = 0.0

                valstr = ' {:3d}   {:20.12f} {:15.10f} '.format(
                    self._num_iter, te, diff_te)
                valstr += '{:15.8f} {:15.8f} {:15.8f} '.format(
                    e_grad, max_grad, diff_den)

                self.ostream.print_header(valstr)
                self.ostream.flush()

    def get_scf_energy(self):
        """
        Gets SCF energy from previous SCF iteration.

        :return:
            The SCF energy.
        """

        return self._scf_energy

    def get_scf_type_str(self):
        """
        Gets string with type of SCF calculation (defined in derrived classes).

        :return:
            The string with type of SCF calculation.
        """

        return 'Undefined'

    def _get_guess_type(self):
        """
        Gets string with type of initial guess (superposition of atomic
        densities or projection of molecular orbitals).

        :return:
            The string with type of initial guess.
        """

        if self.restart:
            return 'Restart from Checkpoint'
        else:
            return 'Superposition of Atomic Densities'

    def _get_acc_type(self):
        """
        Gets string with type of SCF convergence accelerator (DIIS or two level
        DIIS).

        :return:
            The string with type of SCF convergence accelerator.
        """

        if self.acc_type.upper() == 'DIIS':
            return 'Direct Inversion of Iterative Subspace'

        if self.acc_type.upper() == 'L2_DIIS':
            return 'Two Level Direct Inversion of Iterative Subspace'

        return 'Undefined'

    def _delete_mos(self, mol_orbs, mol_eigs):
        """
        Generates trimmed molecular orbital by deleting MOs with coeficients
        exceeding 1.0 / sqrt(ovl_thresh).

        :param mol_orbs:
            The molecular orbitals.
        :param mol_eigs:
            The eigenvalues of molecular orbitals.

        :return:
            The tuple (trimmed molecular orbitals, eigenvalues).
        """

        fmax = 1.0 / math.sqrt(self.ovl_thresh)

        mvec = np.amax(np.abs(mol_orbs), axis=0)

        molist = []
        for i in range(mvec.shape[0]):
            if mvec[i] < fmax:
                molist.append(i)

        return (mol_orbs[:, molist], mol_eigs[molist])

    def compute_s2(self, molecule, scf_tensors):
        """
        Computes expectation value of the S**2 operator.

        :param molecule:
            The molecule.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            Expectation value <S**2>.
        """

        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()

        smat = scf_tensors['S']
        Cocc_a = scf_tensors['C_alpha'][:, :nalpha].copy()
        Cocc_b = scf_tensors['C_beta'][:, :nbeta].copy()

        a_b = float(nalpha - nbeta) / 2.0
        s2_exact = a_b * (a_b + 1.0)

        ovl_a_b = np.matmul(Cocc_a.T, np.matmul(smat, Cocc_b))
        s2 = s2_exact + nbeta - np.sum(ovl_a_b**2)

        return s2

    def _print_ground_state(self, molecule, s2):
        """
        Prints ground state information to output stream.

        :param molecule:
            The molecule.
        :param s2:
            The expectation value of S**2.
        """

        self.ostream.print_blank()

        self.ostream.print_header('Ground State Information'.ljust(92))
        self.ostream.print_header('------------------------'.ljust(92))

        chg = molecule.get_charge()
        valstr = 'Charge of Molecule            :{:5.1f}'.format(chg)
        self.ostream.print_header(valstr.ljust(92))

        mult = molecule.get_multiplicity()
        if self.scf_type == 'restricted':
            valstr = 'Multiplicity (2S+1)           :{:5.1f}'.format(mult)
            self.ostream.print_header(valstr.ljust(92))

        sz = 0.5 * (mult - 1.0)
        valstr = 'Magnetic Quantum Number (M_S) :{:5.1f}'.format(sz)
        self.ostream.print_header(valstr.ljust(92))

        if self.scf_type in ['unrestricted', 'restricted_openshell']:
            valstr = 'Expectation value of S**2     :{:8.4f}'.format(s2)
            self.ostream.print_header(valstr.ljust(92))

        self.ostream.print_blank()

    def _print_energy_components(self):
        """
        Prints SCF energy components to output stream.
        """

        enuc = self._nuc_energy

        e_d4 = self._d4_energy

        e_ef_nuc = self._ef_nuc_energy

        etot = self._iter_data['energy']

        e_el = etot - enuc - e_d4 - e_ef_nuc

        valstr = f'Total Energy                       :{etot:20.10f} a.u.'
        self.ostream.print_header(valstr.ljust(92))

        valstr = f'Electronic Energy                  :{e_el:20.10f} a.u.'
        self.ostream.print_header(valstr.ljust(92))

        valstr = f'Nuclear Repulsion Energy           :{enuc:20.10f} a.u.'
        self.ostream.print_header(valstr.ljust(92))

        if self.dispersion:
            valstr = f'D4 Dispersion Correction           :{e_d4:20.10f} a.u.'
            self.ostream.print_header(valstr.ljust(92))

        if self.electric_field is not None:
            valstr = f'Nuclei in Static Electric Field    :{e_ef_nuc:20.10f} a.u.'
            self.ostream.print_header(valstr.ljust(92))

        self.ostream.print_header(
            '------------------------------------'.ljust(92))

        grad = self._iter_data['gradient_norm']
        valstr = 'Gradient Norm                      :{:20.10f} a.u.'.format(
            grad)
        self.ostream.print_header(valstr.ljust(92))

        self.ostream.print_blank()

        if self.dispersion:
            valstr = '*** Reference for D4 dispersion correction: '
            self.ostream.print_header(valstr.ljust(92))
            valstr = 'E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, '
            valstr += 'S. Spicher, C. Bannwarth'
            self.ostream.print_header(valstr.ljust(92))
            valstr = 'and S. Grimme, J. Chem Phys, 2019, 150, 154122.'
            self.ostream.print_header(valstr.ljust(92))

    def _write_final_hdf5(self, molecule, ao_basis):
        """
        Writes final HDF5 that contains SCF results.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        """

        if self.checkpoint_file is None:
            return

        final_h5_fname = str(
            Path(self.checkpoint_file).with_suffix('.results.h5'))

        if self._dft:
            xc_label = self.xcfun.get_func_label()
        else:
            xc_label = 'HF'

        if self._pe:
            with open(str(self.pe_options['potfile']), 'r') as f_pot:
                potfile_text = '\n'.join(f_pot.readlines())
        else:
            potfile_text = ''

        create_hdf5(final_h5_fname, molecule, ao_basis, xc_label, potfile_text)
        write_scf_results_to_hdf5(final_h5_fname, self.scf_tensors,
                                  self.history)

        self.ostream.print_blank()
        checkpoint_text = 'SCF results written to file: '
        checkpoint_text += final_h5_fname
        self.ostream.print_info(checkpoint_text)
