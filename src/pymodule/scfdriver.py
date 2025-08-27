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
import re

from .oneeints import compute_nuclear_potential_integrals
from .oneeints import compute_electric_dipole_integrals
from .veloxchemlib import TwoCenterElectronRepulsionDriver
from .veloxchemlib import RIFockDriver, SubMatrix
from .veloxchemlib import OverlapDriver, KineticEnergyDriver
from .veloxchemlib import T4CScreener
from .veloxchemlib import XCIntegrator
from .veloxchemlib import mpi_master
from .veloxchemlib import bohr_in_angstrom, hartree_in_kjpermol
from .veloxchemlib import xcfun
from .veloxchemlib import denmat, mat_t
from .veloxchemlib import make_matrix
from .matrix import Matrix
from .aodensitymatrix import AODensityMatrix
from .fockdriver import FockDriver
from .profiler import Profiler
from .griddriver import GridDriver
from .molecularbasis import MolecularBasis
from .molecularorbitals import MolecularOrbitals, molorb
from .sadguessdriver import SadGuessDriver
from .firstorderprop import FirstOrderProperties
from .cpcmdriver import CpcmDriver
from .dispersionmodel import DispersionModel
from .inputparser import (parse_input, print_keywords, print_attributes,
                          get_random_string_parallel)
from .dftutils import get_default_grid_level, print_xc_reference
from .sanitychecks import (molecule_sanity_check, dft_sanity_check,
                           pe_sanity_check, gostshyp_sanity_check, solvation_model_sanity_check)
from .errorhandler import assert_msg_critical
from .checkpoint import create_hdf5, write_scf_results_to_hdf5

try:
    from scipy.linalg import lu_factor, lu_solve
except ImportError:
    pass


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
        - pe_options: The dictionary with options for polarizable embedding.
        - pe_summary: The summary string for polarizable embedding.
        - gostshyp: The flag for running a gostshyp pressure calculation.
        - pressure: The applied hydrostatic pressure.
        - pressure_units: The units of the applied pressure.
        - num_leb_points: The number of Lebedev points per van der Waals sphere.
        - tssf: The tessellation sphere scaling factor.
        - discretization: The surface discretization method.
        - homemade: The flag for running with a homemade cavity for gostshyp pressure calculation
        - tess_file: The data file containing tessellation data for homemade cavity
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

        # pseudo fractional occupation number (pFON)
        # J. Chem. Phys. 110, 695-700 (1999)
        self.pfon = False
        self.pfon_temperature = 1250
        self.pfon_delta_temperature = 50
        self.pfon_nocc = 5
        self.pfon_nvir = 5

        # level-shifting
        # Int. J. Quantum Chem. 7, 699-705 (1973).
        self.level_shifting = 0.0
        self.level_shifting_delta = 0.01

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

        # RI-J
        self.ri_coulomb = False
        self.ri_auxiliary_basis = 'def2-universal-jfit'
        self._ri_drv = None

        # dft
        self.xcfun = None
        self.grid_level = None
        self._dft = False
        self._mol_grid = None
        self.semi_num_exchange = False
        self._semi_num_grid = None

        # polarizable embedding
        self.potfile = None
        self.pe_options = {}
        self._pe = False
        self.embedding = None
        self._embedding_drv = None

        # gostshyp setup
        self._gostshyp = False
        self.pressure = 0.0
        self._pressure_in_input_units = 0.0
        self.pressure_units = 'MPa'
        self.num_leb_points = 110
        self.tssf = 1.2
        self.discretization = 'fixed'
        self.homemade = False #TODO: remove (added for testing of gradient with fixed cavity)
        self.tess_file = None #TODO: remove (added for testing of gradient with fixed cavity)

        # solvation model
        self.solvation_model = None
        self._cpcm = False
        self.cpcm_drv = None
        self.cpcm_epsilon = 78.39
        self.cpcm_grid_per_sphere = 194
        self.cpcm_x = 0

        # point charges (in case we want a simple MM environment without PE)
        self.point_charges = None
        self.qm_vdw_params = None
        self._nuc_mm_energy = 0.0
        self._V_es = None

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

        self._debug = False
        self._block_size_factor = 8
        self._xcfun_ldstaging = 1024
        self._semi_num_ldstaging = 128

        # input keywords
        self._input_keywords = {
            'scf': {
                'ri_coulomb': ('bool', 'use RI-J approximation'),
                'ri_auxiliary_basis': ('str', 'RI-J auxiliary basis set'),
                'acc_type':
                    ('str_upper', 'type of SCF convergence accelerator'),
                'max_iter': ('int', 'maximum number of SCF iterations'),
                'max_err_vecs': ('int', 'maximum number of DIIS error vectors'),
                'pfon': ('bool', 'use pFON to accelerate convergence'),
                'pfon_temperature': ('float', 'pFON temperature'),
                'pfon_delta_temperature': ('float', 'pFON delta temperature'),
                'pfon_nocc':
                    ('int', 'number of occupied orbitals used in pFON'),
                'pfon_nvir': ('int', 'number of virtual orbitals used in pFON'),
                'level_shifting': ('float', 'level shifting parameter'),
                'level_shifting_delta': ('float', 'level shifting delta'),
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
                'point_charges': ('str', 'potential file for point charges'),
                'qm_vdw_params': ('str', 'vdw parameter file for QM atoms'),
                '_debug': ('bool', 'print debug info'),
                '_block_size_factor': ('int', 'block size factor for ERI'),
                '_xcfun_ldstaging': ('int', 'max batch size for DFT grid'),
            },
            'method_settings': {
                'dispersion': ('bool', 'use D4 dispersion correction'),
                'xcfun': ('str_upper', 'exchange-correlation functional'),
                'grid_level': ('int', 'accuracy level of DFT grid (1-8)'),
                'potfile': ('str', 'potential file for polarizable embedding'),
                'pressure': ('float', 'applied hydrostatic pressure'),
                'pressure_units': ('str', 'units of the applied pressure'),
                'num_leb_points': ('int', 'number of grid points per sphere'),
                'tssf': ('float', 'tessellation sphere scaling factor'),
                'discretization': ('str', 'surface discretization method'),
                'homemade': ('bool', 'read tessellation data from file'), #TODO: remove (added for testing of gradient with fixed cavity)
                'tess_file': ('str', 'tessellation data file name'),#TODO: remove (added for testing of gradient with fixed cavity)
                'solvation_model': ('str', 'solvation model'),
                'cpcm_grid_per_sphere':
                    ('int', 'number of grid points per sphere (C-PCM)'),
                'cpcm_epsilon':
                    ('float', 'dielectric constant of solvent (C-PCM)'),
                'cpcm_x': ('float', 'parameter for scaling function (C-PCM)'),
                'electric_field': ('seq_fixed', 'static electric field'),
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

    # @property
    # def filename(self):
    #     """
    #     Getter function for protected filename attribute.
    #     """

    #     return self.filename

    # @filename.setter
    # def filename(self, value):
    #     """
    #     Setter function for protected filename attribute.
    #     """

    #     self.filename = value

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

        pe_sanity_check(self, method_dict)

        #print(self.tssf)
        #print(self.num_leb_points)

        gostshyp_sanity_check(self, method_dict)

        #print(self.tssf)
        #print(self.num_leb_points)

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

        if self.point_charges is not None:
            assert_msg_critical(
                not self._pe,
                'SCF driver: The \'point_charges\' option is incompatible ' +
                'with polarizable embedding')
            # Note: we allow restarting SCF with point charges

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

        # check RI-J
        # for now, force DIIS for RI-J
        # TODO: double check
        if self.ri_coulomb:
            self.acc_type = 'DIIS'

        # check molecule
        molecule_sanity_check(molecule)

        # check dft setup
        dft_sanity_check(self, 'compute')

        # check pe setup
        pe_sanity_check(self)

        # check gostshyp setup
        #gostshyp_sanity_check(self)

        # check solvation model setup
        solvation_model_sanity_check(self)
        if self._cpcm:
            self.cpcm_drv = CpcmDriver(self.comm, self.ostream)
            self.cpcm_drv.grid_per_sphere = self.cpcm_grid_per_sphere
            self.cpcm_drv.epsilon = self.cpcm_epsilon
            self.cpcm_drv.x = self.cpcm_x
        else:
            self.cpcm_drv = None

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

        # generate integration grid
        if self._dft:
            print_xc_reference(self.xcfun, self.ostream)

            grid_drv = GridDriver(self.comm)
            grid_level = (get_default_grid_level(self.xcfun)
                          if self.grid_level is None else self.grid_level)
            grid_drv.set_level(grid_level)

            grid_t0 = tm.time()
            self._mol_grid = grid_drv.generate(molecule, self._xcfun_ldstaging)
            n_grid_points = self._mol_grid.number_of_points()
            self.ostream.print_info(
                'Molecular grid with {0:d} points generated in {1:.2f} sec.'.
                format(n_grid_points,
                       tm.time() - grid_t0))
            self.ostream.print_blank()
            
        # generate integration grid for seminumerical exchange
        if self.semi_num_exchange:
            grid_drv = GridDriver(self.comm)
            grid_drv.set_level(4)

            grid_t0 = tm.time()
            self._semi_num_grid = grid_drv.generate(molecule, self._semi_num_ldstaging)
            n_grid_points = self._semi_num_grid.number_of_points()
            self.ostream.print_info(
                'Molecular grid for seminumerical exchange with {0:d} points generated in {1:.2f} sec.'.
                format(n_grid_points,
                       tm.time() - grid_t0))
            self.ostream.print_blank()

        # D4 dispersion correction
        if self.dispersion:
            if self.rank == mpi_master():
                disp = DispersionModel()
                xc_label = self.xcfun.get_func_label() if self._dft else 'HF'
                disp.compute(molecule, xc_label)
                self._d4_energy = disp.get_energy()

                dftd4_info = 'Using the D4 dispersion correction.'
                self.ostream.print_info(dftd4_info)
                self.ostream.print_blank()
                for dftd4_ref in disp.get_references():
                    self.ostream.print_reference(dftd4_ref)
                self.ostream.print_blank()
                self.ostream.flush()

            else:
                self._d4_energy = 0.0

            self._d4_energy = self.comm.bcast(self._d4_energy,
                                              root=mpi_master())
        else:
            self._d4_energy = 0.0

        # set up polarizable continuum model
        if self._cpcm:
            cpcm_info = 'Using C-PCM with the ISWIG discretization method.'
            self.ostream.print_info(cpcm_info)
            self.ostream.print_blank()
            iswig_ref = 'A. W. Lange, J. M. Herbert,'
            iswig_ref += ' J. Chem. Phys. 2010, 133, 244111.'
            self.ostream.print_reference(iswig_ref)
            self.ostream.print_blank()
            self.ostream.flush()

            (self._cpcm_grid,
             self._cpcm_sw_func) = self.cpcm_drv.generate_cpcm_grid(molecule)
            self._cpcm_Amat = self.cpcm_drv.form_matrix_A(
                self._cpcm_grid, self._cpcm_sw_func)
            self._cpcm_Bmat = self.cpcm_drv.form_matrix_B(
                self._cpcm_grid, molecule)
            self._cpcm_Bzvec = np.dot(self._cpcm_Bmat,
                                      molecule.get_element_ids())

        # set up polarizable embedding
        if self._pe:

            # TODO: print PyFraME info

            pot_info = 'Reading polarizable embedding: {}'.format(
                self.pe_options['potfile'])
            self.ostream.print_info(pot_info)
            self.ostream.print_blank()

            assert_msg_critical(
                self.embedding['settings']['embedding_method'] == 'PE',
                'PolarizableEmbedding: Invalid embedding_method. Only PE is supported.'
            )

            settings = self.embedding['settings']

            from .embedding import PolarizableEmbeddingSCF

            self._embedding_drv = PolarizableEmbeddingSCF(
                molecule=molecule,
                ao_basis=ao_basis,
                options=self.embedding,
                comm=self.comm)

            emb_info = 'Embedding settings:'
            self.ostream.print_info(emb_info)

            for key in ['embedding_method', 'environment_energy']:
                if key in settings:
                    self.ostream.print_info(f'- {key:<15s} : {settings[key]}')

            if 'vdw' in settings:
                self.ostream.print_info(f'- {"vdw":<15s}')
                for key in ['method', 'combination_rule']:
                    self.ostream.print_info(
                        f'  - {key:<15s} : {settings["vdw"][key]}')

            if 'induced_dipoles' in settings:
                self.ostream.print_info(f'- {"induced_dipoles":<15s}')
                default_values = {
                    'solver': 'jacobi',
                    'threshold': 1e-8,
                    'max_iterations': 100,
                    'mic': False,
                }
                for key, default in default_values.items():
                    value = settings['induced_dipoles'].get(key, default)
                    self.ostream.print_info(f'  - {key:<15s} : {value}')

            self.ostream.print_blank()

        # set up point charges without PE
        elif self.point_charges is not None:

            pot_info = 'Preparing for QM/MM calculation.'
            self.ostream.print_info(pot_info)
            self.ostream.print_blank()

            if isinstance(self.point_charges, str):
                potfile = self.point_charges

                pot_info = 'Reading point charges: {}'.format(potfile)
                self.ostream.print_info(pot_info)
                self.ostream.print_blank()

                if self.rank == mpi_master():
                    with Path(potfile).open('r') as fh:
                        lines = fh.read().strip().splitlines()
                        try:
                            npoints = int(lines[0].strip())
                        except (ValueError, TypeError):
                            assert_msg_critical(
                                False, 'potfile: Invalid number of points')
                        assert_msg_critical(
                            npoints == len(lines[2:]),
                            'potfile: Inconsistent number of points')
                        self.point_charges = np.zeros((6, npoints))
                        for idx, line in enumerate(lines[2:]):
                            content = line.split()
                            label, x, y, z, q = content[:5]
                            self.point_charges[0, idx] = (float(x) /
                                                          bohr_in_angstrom())
                            self.point_charges[1, idx] = (float(y) /
                                                          bohr_in_angstrom())
                            self.point_charges[2, idx] = (float(z) /
                                                          bohr_in_angstrom())
                            self.point_charges[3, idx] = float(q)
                            if self.qm_vdw_params is not None:
                                # Note: read MM vdw parameters only when QM vdw
                                # parameters are present
                                sigma, epsilon = content[5:7]
                                self.point_charges[4, idx] = float(sigma)
                                self.point_charges[5, idx] = float(epsilon)
                else:
                    self.point_charges = None
                self.point_charges = self.comm.bcast(self.point_charges,
                                                     root=mpi_master())

            if isinstance(self.qm_vdw_params, str):
                vdw_param_file = self.qm_vdw_params

                vdw_info = 'Reading vdw parameters for QM atoms: {}'.format(
                    vdw_param_file)
                self.ostream.print_info(vdw_info)
                self.ostream.print_blank()

                if self.rank == mpi_master():
                    natoms = molecule.number_of_atoms()
                    self.qm_vdw_params = np.zeros((natoms, 2))
                    with Path(vdw_param_file).open('r') as fh:
                        for a in range(natoms):
                            sigma, epsilon = fh.readline().split()
                            self.qm_vdw_params[a, 0] = float(sigma)
                            self.qm_vdw_params[a, 1] = float(epsilon)
                else:
                    self.qm_vdw_params = None
                self.qm_vdw_params = self.comm.bcast(self.qm_vdw_params,
                                                     root=mpi_master())

            # nuclei - point charges interaction
            self._nuc_mm_energy = 0.0

            natoms = molecule.number_of_atoms()
            coords = molecule.get_coordinates_in_bohr()
            nuclear_charges = molecule.get_element_ids()
            npoints = self.point_charges.shape[1]

            for a in range(self.rank, natoms, self.nodes):
                xyz_a = coords[a]
                chg_a = nuclear_charges[a]

                for p in range(npoints):
                    xyz_p = self.point_charges[:3, p]
                    chg_p = self.point_charges[3, p]
                    r_ap = np.linalg.norm(xyz_a - xyz_p)

                    self._nuc_mm_energy += chg_a * chg_p / r_ap

            if self.qm_vdw_params is not None:
                vdw_ene = 0.0

                for a in range(self.rank, natoms, self.nodes):
                    xyz_i = coords[a]
                    sigma_i = self.qm_vdw_params[a, 0]
                    epsilon_i = self.qm_vdw_params[a, 1]

                    for p in range(npoints):
                        xyz_j = self.point_charges[:3, p]
                        sigma_j = self.point_charges[4, p]
                        epsilon_j = self.point_charges[5, p]

                        distance_ij = np.linalg.norm(xyz_i - xyz_j)
                        # bohr to nm
                        distance_ij *= bohr_in_angstrom() * 0.1

                        epsilon_ij = np.sqrt(epsilon_i * epsilon_j)
                        sigma_ij = 0.5 * (sigma_i + sigma_j)

                        sigma_r_6 = (sigma_ij / distance_ij)**6
                        sigma_r_12 = sigma_r_6**2

                        vdw_ene += 4.0 * epsilon_ij * (sigma_r_12 - sigma_r_6)

                # kJ/mol to Hartree
                vdw_ene /= hartree_in_kjpermol()

                self._nuc_mm_energy += vdw_ene

            self._nuc_mm_energy = self.comm.allreduce(self._nuc_mm_energy)

        # set up gostshyp method by creating a surface tessellation
        if self._gostshyp:
            from .gostshyp import GostshypDriver
            self._gostshyp_drv = GostshypDriver(molecule, ao_basis,
                    self.pressure, self.pressure_units, self.comm, self.ostream)

            tessellation_settings = {
                'num_leb_points': self.num_leb_points,
                'tssf': self.tssf,
                'discretization': self.discretization,
                'filename': self.filename,
                'homemade': self.homemade, #TODO: remove (added for testing of gradient with fixed cavity)
                'tess_file': self.tess_file #TODO: remove (added for testing of gradient with fixed cavity)
            }

            tess_t0 = tm.time()
            tessellation = self._gostshyp_drv.generate_tessellation(
                                                        tessellation_settings)

            tess_info = 'Surface tessellation with '
            tess_info += '{0:d} grid points generated in {1:.2f} sec.'.format(
                tessellation.shape[1], tm.time() - tess_t0)
            self.ostream.print_info(tess_info)
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

        if self.scf_type == 'restricted':
            density_type = 'restricted'
        else:
            density_type = 'unrestricted'

        natoms = molecule.number_of_atoms()
        unpaired_electrons_on_atoms = [0 for a in range(natoms)]

        if (self.guess_unpaired_electrons and density_type == 'restricted'):
            warn_msg = 'Ignoring "guess_unpaired_electrons" in '
            warn_msg += 'spin-restricted SCF calculation.'
            self.ostream.print_warning(warn_msg)
            self.ostream.print_blank()

        if (self.guess_unpaired_electrons and density_type == 'unrestricted'):
            for entry in self.guess_unpaired_electrons.split(','):
                m = re.search(r'^(.*)\((.*)\)$', entry.strip())
                assert_msg_critical(
                    m is not None,
                    'Initial Guess: Invalid input for unpaired electrons')
                atom_index = int(m.group(1).strip()) - 1
                num_unpaired_elec = float(m.group(2).strip())
                unpaired_electrons_on_atoms[atom_index] = num_unpaired_elec

            sad_drv.set_number_of_unpaired_electrons_on_atoms(
                unpaired_electrons_on_atoms)

            guess_msg = 'Generating initial guess with '
            guess_msg += 'user-provided information...'
            self.ostream.print_info(guess_msg)

            labels = molecule.get_labels()
            for a, (num_unpaired_elec, atom_name) in enumerate(
                    zip(unpaired_electrons_on_atoms, labels)):
                if num_unpaired_elec != 0:
                    spin = 'alpha' if num_unpaired_elec > 0 else 'beta '
                    abs_num_unpaired_elec = abs(num_unpaired_elec)
                    self.ostream.print_info(
                        f'  {abs_num_unpaired_elec} unpaired {spin} ' +
                        f'electrons on atom {a + 1} ({atom_name})')
            self.ostream.print_blank()

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

            n_ao = basis.get_dimensions_of_basis()

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
        self.write_checkpoint(molecule.get_element_ids(), basis.get_label())

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

        if self.point_charges is not None and not self._first_step:
            t0point_charges = tm.time()

            ave, res = divmod(self.point_charges.shape[1], self.nodes)
            counts = [ave + 1 if p < res else ave for p in range(self.nodes)]

            start = sum(counts[:self.rank])
            end = sum(counts[:self.rank + 1])

            charges = self.point_charges[3, :].copy()[start:end]
            coords = self.point_charges[:3, :].T.copy()[start:end, :]

            V_es = compute_nuclear_potential_integrals(molecule, ao_basis,
                                                       charges, coords)

            self._V_es = self.comm.reduce(V_es, root=mpi_master())
            if self.rank != mpi_master():
                self._V_es = np.zeros(V_es.shape)

            self.ostream.print_info(
                'Point charges one-electron integral computed in' +
                ' {:.2f} sec.'.format(tm.time() - t0point_charges))
            self.ostream.print_blank()

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

            # TODO: double check if it is necessary to tighten threshold

            self.eri_thresh = self.eri_thresh_tight

            if self.rank == mpi_master():
                self.ostream.print_info('ERI screening threshold tightened to' +
                                        ' {:.1e}.'.format(self.eri_thresh))
                self.ostream.print_blank()

        self._density = tuple([x.copy() for x in den_mat])

        self._print_debug_info('before screener')
        if self.rank == mpi_master():
            screener = T4CScreener()
            screener.partition(ao_basis, molecule, 'eri')
        else:
            screener = None
        self._print_debug_info('after  screener')
        screener = self.comm.bcast(screener, root=mpi_master())
        self._print_debug_info('after  bcast screener')

        profiler.check_memory_usage('Initial guess')

        if self.ri_coulomb and self.rank == mpi_master():
            assert_msg_critical(
                ao_basis.get_label().lower().startswith('def2-'),
                'SCF Driver: Invalid basis set for RI-J')

            self.ostream.print_info(
                'Using the resolution of the identity (RI) approximation.')
            self.ostream.print_blank()
            self.ostream.flush()

            basis_ri_j = MolecularBasis.read(molecule, self.ri_auxiliary_basis)

            self.ostream.print_info('Dimension of RI auxiliary basis set ' +
                                    f'({self.ri_auxiliary_basis.upper()}): ' +
                                    f'{basis_ri_j.get_dimensions_of_basis()}')
            self.ostream.print_blank()
            self.ostream.flush()

            ri_prep_t0 = tm.time()

            t2c_drv = TwoCenterElectronRepulsionDriver()
            mat_j = t2c_drv.compute(molecule, basis_ri_j)
            mat_j_np = mat_j.to_numpy()

            self.ostream.print_info('Two-center integrals for RI done in ' +
                                    f'{tm.time() - ri_prep_t0:.2f} sec.')
            self.ostream.print_blank()

            ri_prep_t0 = tm.time()

            if 'scipy' in sys.modules:
                lu, piv = lu_factor(mat_j_np)
                inv_mat_j_np = lu_solve((lu, piv), np.eye(mat_j_np.shape[0]))
            else:
                inv_mat_j_np = np.linalg.inv(mat_j_np)

            self.ostream.print_info(
                f'Matrix inversion for RI done in {tm.time() - ri_prep_t0:.2f} sec.'
            )
            self.ostream.print_blank()

            ri_prep_t0 = tm.time()

            inv_mat_j = SubMatrix(
                [0, 0, inv_mat_j_np.shape[0], inv_mat_j_np.shape[1]])
            inv_mat_j.set_values(inv_mat_j_np)

            self._ri_drv = RIFockDriver(inv_mat_j)
            self._ri_drv.prepare_buffers(molecule, ao_basis, basis_ri_j)

            self.ostream.print_info(
                f'Buffer preparation for RI done in {tm.time() - ri_prep_t0:.2f} sec.'
            )
            self.ostream.print_blank()
            self.ostream.flush()

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

            fock_mat, vxc_mat, vkx_mat, e_emb, V_emb, e_pr, V_pr = self._comp_2e_fock(
                den_mat, molecule, ao_basis, screener, e_grad, profiler)
            
            self._e_gost = e_pr

            profiler.start_timer('ErrVec')

            e_el = self._comp_energy(fock_mat, vxc_mat, vkx_mat, e_emb, e_pr, kin_mat,
                                     npot_mat, den_mat)

            self._comp_full_fock(fock_mat, vxc_mat, vkx_mat, V_emb, V_pr, kin_mat, npot_mat)

            if self._cpcm:
                if self.scf_type == 'restricted':
                    Cvec = self.cpcm_drv.form_vector_C(molecule, ao_basis,
                                                       self._cpcm_grid,
                                                       den_mat[0] * 2.0)
                else:
                    Cvec = self.cpcm_drv.form_vector_C(molecule, ao_basis,
                                                       self._cpcm_grid,
                                                       den_mat[0] + den_mat[1])

                if self.rank == mpi_master():
                    scale_f = -(self.cpcm_drv.epsilon - 1) / (
                        self.cpcm_drv.epsilon + self.cpcm_drv.x)
                    rhs = scale_f * (self._cpcm_Bzvec + Cvec)
                    self._cpcm_q = np.linalg.solve(self._cpcm_Amat, rhs)

                    e_sol = self.cpcm_drv.compute_solv_energy(
                        self._cpcm_Bzvec, Cvec, self._cpcm_q)
                    e_el += e_sol
                    self.cpcm_epol = e_sol
                else:
                    self._cpcm_q = None
                    self.cpcm_epol = None
                self._cpcm_q = self.comm.bcast(self._cpcm_q, root=mpi_master())

                Fock_sol = self.cpcm_drv.get_contribution_to_Fock(
                    molecule, ao_basis, self._cpcm_grid, self._cpcm_q)

                if self.rank == mpi_master():
                    fock_mat[0] += Fock_sol
                    if self.scf_type != 'restricted':
                        fock_mat[1] += Fock_sol

            if (self.rank == mpi_master() and i > 0 and
                    self.level_shifting > 0.0):

                self.ostream.print_info(
                    f'Applying level-shifting ({self.level_shifting:.2f}au)')

                C_alpha = self.molecular_orbitals.alpha_to_numpy()
                nocc_a = molecule.number_of_alpha_electrons()
                fmo_a = np.linalg.multi_dot([C_alpha.T, fock_mat[0], C_alpha])
                for idx in range(nocc_a, fmo_a.shape[0]):
                    fmo_a[idx, idx] += self.level_shifting
                fock_mat[0] = np.linalg.multi_dot(
                    [S, C_alpha, fmo_a, C_alpha.T, S])

                if self.scf_type != 'restricted':

                    C_beta = self.molecular_orbitals.beta_to_numpy()
                    nocc_b = molecule.number_of_beta_electrons()
                    fmo_b = np.linalg.multi_dot([C_beta.T, fock_mat[1], C_beta])
                    for idx in range(nocc_b, fmo_b.shape[0]):
                        fmo_b[idx, idx] += self.level_shifting
                    fock_mat[1] = np.linalg.multi_dot(
                        [S, C_beta, fmo_b, C_beta.T, S])

                self.level_shifting -= self.level_shifting_delta
                if self.level_shifting < 0.0:
                    self.level_shifting = 0.0

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

            # threshold for deactivating pseudo-FON and level-shifting
            if e_grad < 1.0e-4:
                if self.pfon:
                    self.pfon_temperature = 0
                if self.level_shifting > 0.0:
                    self.level_shifting = 0.0

            # compute density change and energy change

            diff_den = self._comp_density_change(den_mat, self._density)

            e_scf = (e_el + self._nuc_energy + self._nuc_mm_energy +
                     self._d4_energy + self._ef_nuc_energy)

            e_scf = self.comm.bcast(e_scf, root=mpi_master())

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

            profiler.start_timer('EffFock')

            self._store_diis_data(fock_mat, den_mat, ovl_mat, e_grad)

            eff_fock_mat = self._get_effective_fock(fock_mat, ovl_mat, oao_mat)

            profiler.stop_timer('EffFock')

            profiler.start_timer('NewMO')

            self._molecular_orbitals = self._gen_molecular_orbitals(
                molecule, eff_fock_mat, oao_mat)

            if self.pfon:
                self.pfon_temperature -= self.pfon_delta_temperature
                if self.pfon_temperature < 0:
                    self.pfon_temperature = 0

            if self._mom is not None:
                self._apply_mom(molecule, ovl_mat)

            self._update_mol_orbs_phase()

            profiler.stop_timer('NewMO')

            profiler.start_timer('NewDens')

            if self.rank == mpi_master():
                den_mat = self.molecular_orbitals.get_density(molecule)
            else:
                den_mat = None
            den_mat = self.comm.bcast(den_mat, root=mpi_master())

            profiler.stop_timer('NewDens')

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

                # for backward compatibility only
                self._scf_tensors['F'] = (F_alpha, F_beta)

                if self._dft:
                    # dft info
                    self._scf_tensors['xcfun'] = self.xcfun.get_func_label()
                    if self.grid_level is not None:
                        self._scf_tensors['grid_level'] = self.grid_level

                if self._pe:
                    # pe info, energy and potential matrix
                    self._scf_tensors['potfile'] = self.potfile
                    self._scf_tensors['E_emb'] = e_emb
                    self._scf_tensors['F_emb'] = V_emb

                if self.point_charges is not None:
                    self._scf_tensors['point_charges'] = self.point_charges
                if self.qm_vdw_params is not None:
                    self._scf_tensors['qm_vdw_params'] = self.qm_vdw_params

                if self.solvation_model is not None:
                    self._scf_tensors['solvation_model'] = self.solvation_model
                    self._scf_tensors[
                        'dielectric_constant'] = self.cpcm_drv.epsilon

            else:
                self._scf_tensors = None
                self._density = AODensityMatrix()

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

        need_exit = False

        if self.program_end_time is not None:
            remaining_hours = (self.program_end_time -
                               datetime.now()).total_seconds() / 3600
            # exit gracefully when the remaining time is not sufficient to
            # complete the next iteration (plus 25% to be on the safe side).
            if remaining_hours < iter_in_hours * 1.25:
                need_exit = True

        need_exit = self.comm.bcast(need_exit, root=mpi_master())

        return need_exit

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

        self.write_checkpoint(molecule.get_element_ids(), basis.get_label())

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
            ovl_mat = ovl_mat.to_numpy()

            ovl_dt = tm.time() - t0
            t0 = tm.time()

            kin_drv = KineticEnergyDriver()
            kin_mat = kin_drv.compute(molecule, basis)
            kin_mat = kin_mat.to_numpy()

            kin_dt = tm.time() - t0
        else:
            ovl_mat = None
            kin_mat = None

        ovl_mat = self.comm.bcast(ovl_mat, root=mpi_master())
        kin_mat = self.comm.bcast(kin_mat, root=mpi_master())

        t0 = tm.time()

        if molecule.number_of_atoms() >= self.nodes and self.nodes > 1:
            npot_mat = self._comp_npot_mat_parallel(molecule, basis)
        else:
            npot_mat = compute_nuclear_potential_integrals(molecule, basis)

        npot_dt = tm.time() - t0

        npot_mat = self.comm.bcast(npot_mat, root=mpi_master())

        t0 = tm.time()

        if self.electric_field is not None:
            if molecule.get_charge() != 0:
                coords = molecule.get_coordinates_in_bohr()
                nuclear_charges = molecule.get_element_ids()
                self._dipole_origin = np.sum(coords.T * nuclear_charges,
                                             axis=1) / np.sum(nuclear_charges)
            else:
                self._dipole_origin = np.zeros(3)

            dipole_mats = compute_electric_dipole_integrals(
                molecule, basis, list(self._dipole_origin))
        else:
            dipole_mats = None

        dipole_dt = tm.time() - t0

        if self.rank == mpi_master() and self.print_level > 1:

            self.ostream.print_info('Overlap matrix computed in' +
                                    ' {:.2f} sec.'.format(ovl_dt))
            self.ostream.print_blank()

            self.ostream.print_info('Kinetic energy matrix computed in' +
                                    ' {:.2f} sec.'.format(kin_dt))
            self.ostream.print_blank()

            self.ostream.print_info('Nuclear potential matrix computed in' +
                                    ' {:.2f} sec.'.format(npot_dt))
            self.ostream.print_blank()

            if self.electric_field is not None:
                self.ostream.print_info('Electric dipole matrices computed in' +
                                        ' {:.2f} sec.'.format(dipole_dt))
                self.ostream.print_blank()

            self.ostream.flush()

        return ovl_mat, kin_mat, npot_mat, dipole_mats

    def _comp_npot_mat_parallel(self, molecule, basis):
        """
        Computes one-electron nuclear potential integral in parallel.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.

        :return:
            The one-electron nuclear potential matrix.
        """

        ave, res = divmod(molecule.number_of_atoms(), self.nodes)
        counts = [ave + 1 if p < res else ave for p in range(self.nodes)]

        start = sum(counts[:self.rank])
        end = sum(counts[:self.rank + 1])

        charges = molecule.get_element_ids()[start:end]
        coords = molecule.get_coordinates_in_bohr()[start:end, :]

        npot_mat = compute_nuclear_potential_integrals(molecule, basis, charges,
                                                       coords)

        npot_mat = self.comm.reduce(npot_mat, root=mpi_master())

        return npot_mat

    def _print_debug_info(self, label):
        """
        Prints debug information.

        :param label:
            The label of debug information.
        """

        if self._debug:
            profiler = Profiler()
            self.ostream.print_info(f'==DEBUG==   available memory {label}: ' +
                                    profiler.get_available_memory())
            self.ostream.flush()

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

        fock_mat, vxc_mat, vkx_mat, e_emb, V_emb, e_pr, V_pr = self._comp_2e_fock_single_comm(
            den_mat, molecule, basis, screener, e_grad, profiler)

        return fock_mat, vxc_mat, vkx_mat, e_emb, V_emb, e_pr, V_pr

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

            den_mat_for_Ka = self.comm.bcast(den_mat_for_Ka, root=mpi_master())
            den_mat_for_Kb = self.comm.bcast(den_mat_for_Kb, root=mpi_master())
            den_mat_for_Jab = self.comm.bcast(den_mat_for_Jab,
                                              root=mpi_master())

        if e_grad is None:
            thresh_int = int(-math.log10(self.eri_thresh))
        else:
            thresh_int = int(-math.log10(self._get_dyn_threshold(e_grad)))

        eri_t0 = tm.time()

        fock_drv = FockDriver(self.comm)
        fock_drv._set_block_size_factor(self._block_size_factor)

        # determine fock_type and exchange_scaling_factor
        fock_type = '2jk'
        exchange_scaling_factor = 1.0
        if self._dft and not self._first_step:
            if self.xcfun.is_hybrid():
                if self.semi_num_exchange:
                    fock_type = 'j'
                    exchange_scaling_factor = self.xcfun.get_frac_exact_exchange()
                else:
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

        self._print_debug_info('before Fock build')

        if self.scf_type == 'restricted':
            # restricted SCF
            self._print_debug_info('before rest Fock build')

            if self.ri_coulomb:
                assert_msg_critical(
                    fock_type == 'j',
                    'SCF driver: RI is only applicable to pure DFT functional')

            if self.ri_coulomb and fock_type == 'j':
                if self.rank == mpi_master():
                    fock_mat = self._ri_drv.compute(den_mat_for_fock, 'j')
                    fock_mat_np = fock_mat.to_numpy()
                    fock_mat = Matrix()
                else:
                    fock_mat_np = np.zeros(den_mat[0].shape)
            else:
                fock_mat = fock_drv.compute(screener, den_mat_for_fock,
                                            fock_type, exchange_scaling_factor,
                                            0.0, thresh_int)
                fock_mat_np = fock_mat.to_numpy()
                fock_mat = Matrix()

            self._print_debug_info('after  rest Fock build')

            if fock_type == 'j':
                # for pure functional
                fock_mat_np *= 2.0

            if need_omega:
                # for range-separated functional
                self._print_debug_info('before rest erf Fock build')
                fock_mat = fock_drv.compute(screener, den_mat_for_fock, 'kx_rs',
                                            erf_k_coef, omega, thresh_int)
                self._print_debug_info('after  rest erf Fock build')

                fock_mat_np -= fock_mat.to_numpy()
                fock_mat = Matrix()

            fock_mat_np = self.comm.reduce(fock_mat_np, root=mpi_master())

            den_mat_for_fock = Matrix()

            if self.rank == mpi_master():
                # Note: make fock_mat a list
                fock_mat = [fock_mat_np]
            else:
                fock_mat = None

        else:
            # unrestricted SCF or restricted open-shell SCF

            if self.ri_coulomb:
                assert_msg_critical(
                    fock_type == 'j',
                    'SCF driver: RI is only applicable to pure DFT functional')

            if fock_type == 'j':
                # for pure functional
                # den_mat_for_Jab is D_total
                self._print_debug_info('before unrest Fock J build')

                if self.ri_coulomb:
                    if self.rank == mpi_master():
                        fock_mat = self._ri_drv.compute(den_mat_for_Jab, 'j')
                        J_ab_np = fock_mat.to_numpy()
                        fock_mat = Matrix()
                    else:
                        J_ab_np = np.zeros(den_mat[0].shape)
                else:
                    fock_mat = fock_drv.compute(screener, den_mat_for_Jab, 'j',
                                                0.0, 0.0, thresh_int)
                    J_ab_np = fock_mat.to_numpy()
                    fock_mat = Matrix()

                self._print_debug_info('after  unrest Fock J build')

                fock_mat_a_np = J_ab_np
                fock_mat_b_np = J_ab_np.copy()

            else:
                self._print_debug_info('before unrest Fock Ka build')
                fock_mat = fock_drv.compute(screener, den_mat_for_Ka, 'kx',
                                            exchange_scaling_factor, 0.0,
                                            thresh_int)
                self._print_debug_info('after  unrest Fock Kb build')

                K_a_np = fock_mat.to_numpy()
                fock_mat = Matrix()

                self._print_debug_info('before unrest Fock Kb build')
                fock_mat = fock_drv.compute(screener, den_mat_for_Kb, 'kx',
                                            exchange_scaling_factor, 0.0,
                                            thresh_int)
                self._print_debug_info('after  unrest Fock Kb build')

                K_b_np = fock_mat.to_numpy()
                fock_mat = Matrix()

                self._print_debug_info('before unrest Fock J build')
                fock_mat = fock_drv.compute(screener, den_mat_for_Jab, 'j',
                                            exchange_scaling_factor, 0.0,
                                            thresh_int)
                self._print_debug_info('after  unrest Fock J build')

                J_ab_np = fock_mat.to_numpy()
                fock_mat = Matrix()

                fock_mat_a_np = J_ab_np - K_a_np
                fock_mat_b_np = J_ab_np - K_b_np

            if need_omega:
                # for range-separated functional
                self._print_debug_info('before unrest erf Fock Ka build')
                fock_mat = fock_drv.compute(screener, den_mat_for_Ka, 'kx_rs',
                                            erf_k_coef, omega, thresh_int)
                self._print_debug_info('after  unrest erf Fock Ka build')

                fock_mat_a_np -= fock_mat.to_numpy()
                fock_mat = Matrix()

                self._print_debug_info('before unrest erf Fock Kb build')
                fock_mat = fock_drv.compute(screener, den_mat_for_Kb, 'kx_rs',
                                            erf_k_coef, omega, thresh_int)
                self._print_debug_info('after  unrest erf Fock Kb build')

                fock_mat_b_np -= fock_mat.to_numpy()
                fock_mat = Matrix()

            fock_mat_a_np = self.comm.reduce(fock_mat_a_np, root=mpi_master())
            fock_mat_b_np = self.comm.reduce(fock_mat_b_np, root=mpi_master())

            den_mat_for_Ka = Matrix()
            den_mat_for_Kb = Matrix()
            den_mat_for_Jab = Matrix()

            if self.rank == mpi_master():
                # Note: make fock_mat a list
                fock_mat = [fock_mat_a_np, fock_mat_b_np]
            else:
                fock_mat = None

        self._print_debug_info('after  Fock build')

        if self.timing:
            profiler.add_timing_info('FockERI', tm.time() - eri_t0)
        vxc_t0 = tm.time()

        if self._dft and not self._first_step:
            if self.xcfun.get_func_type() in [xcfun.lda, xcfun.gga, xcfun.mgga]:
                xc_drv = XCIntegrator()
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
                            
            
        if self.semi_num_exchange:
            kx_drv = XCIntegrator()
            vkx_mat = kx_drv.integrate_kx_fock(molecule, basis, den_mat, self._semi_num_grid, exchange_scaling_factor)
        else:
            vkx_mat = None
                                        
        if self.timing and self._dft:
            profiler.add_timing_info('FockXC', tm.time() - vxc_t0)
        pe_t0 = tm.time()

        if self._pe and not self._first_step:
            if self.scf_type == 'restricted':
                density_matrix = 2.0 * den_mat[0]
            else:
                density_matrix = den_mat[0] + den_mat[1]
            from .embedding import PolarizableEmbeddingSCF
            assert_msg_critical(
                isinstance(self._embedding_drv, PolarizableEmbeddingSCF),
                'ScfDriver: Inconsistent embedding driver for SCF')
            e_emb, V_emb = self._embedding_drv.compute_pe_contributions(
                density_matrix=density_matrix)
        elif self.point_charges is not None and not self._first_step:
            if self.scf_type == 'restricted':
                density_matrix = 2.0 * den_mat[0]
            else:
                density_matrix = den_mat[0] + den_mat[1]
            e_emb = np.sum(density_matrix * self._V_es)
            V_emb = self._V_es
        else:
            e_emb, V_emb = 0.0, None

        if self.timing and self._pe:
            profiler.add_timing_info('FockPE', tm.time() - pe_t0)
        
        gostshyp_t0 = tm.time()
        if self._gostshyp and not self._first_step:
            # which density matrix is needed?
            if self.scf_type == 'restricted':
                density_matrix = 2.0 * den_mat[0]
            else:
                density_matrix = den_mat[0] + den_mat[1]
            e_pr, V_pr = self._gostshyp_drv.get_gostshyp_contribution(density_matrix)
            #e_pr, V_pr, g_time, f_time = self._gostshyp_drv.get_gostshyp_contribution(density_matrix)
            #print('Energy contribution from  GOSTSHYP: ', e_pr) #should be printed in ostream?
            #print(g_time, f_time)
        else:
            e_pr, V_pr = 0.0, None

        if self.timing and self._gostshyp and not self._first_step:
            profiler.add_timing_info('FockGOSTSHYP', tm.time() - gostshyp_t0)
            # profiler.add_timing_info('g_mat', g_time)
            # profiler.add_timing_info('f_mat', f_time)
            # profiler.add_timing_info('g_mat+f_mat', g_time + f_time)

        return fock_mat, vxc_mat, vkx_mat, e_emb, V_emb, e_pr, V_pr

    def _comp_energy(self, fock_mat, vxc_mat, vkx_mat, e_emb, e_pr, kin_mat, npot_mat,
                     den_mat):
        """
        Computes the sum of SCF energy components: electronic energy, kinetic
        energy, and nuclear potential energy.

        :param fock_mat:
            The Fock/Kohn-Sham matrix (only 2e-part).
        :param vxc_mat:
            The Vxc matrix.
        :param vkx_mat:
            The seminumerical exchange matrix.
        :param e_emb:
            The embedding energy.
        :param e_pr:
            The pressure (GOSTSHYP) energy.
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
                if self.semi_num_exchange:
                    e_ee -= np.sum(D[0] * vkx_mat.alpha_to_numpy())
            else:
                e_ee = 0.5 * (np.sum(D[0] * F[0]) + np.sum(D[1] * F[1]))
                e_kin = np.sum((D[0] + D[1]) * T)
                e_en = np.sum((D[0] + D[1]) * V)

            if self._dft and not self._first_step:
                e_ee += xc_ene

            if self._pe and not self._first_step:
                e_ee += e_emb

            if self._gostshyp and not self._first_step:
                e_ee += e_pr
            
            elif self.point_charges is not None and not self._first_step:
                e_ee += e_emb

            e_sum = e_ee + e_kin + e_en
        else:
            e_sum = 0.0
        e_sum = self.comm.bcast(e_sum, root=mpi_master())

        return e_sum

    def _comp_full_fock(self, fock_mat, vxc_mat, vkx_mat, V_emb, pr_mat, kin_mat, npot_mat):
        """
        Computes full Fock/Kohn-Sham matrix by adding to 2e-part of
        Fock/Kohn-Sham matrix the kinetic energy and nuclear potential
        matrices.

        :param fock_mat:
            The Fock/Kohn-Sham matrix (2e-part).
        :param vxc_mat:
            The Vxc matrix.
        :param vkx_mat:
            The seminumerical exchamge matrix.
        :param V_emb:
            The embedding Fock matrix contributions.        
        :param pr_mat:
            The pressure (GOSTSHYP) matrix.
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
                                              
        np_kxmat_a, np_kxmat_b = None, None
        if self.semi_num_exchange:
            np_kxmat_a = self.comm.reduce(vkx_mat.alpha_to_numpy(),
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
                    
            if self.semi_num_exchange:
                fock_mat[0] -= np_kxmat_a
                                
            if self._pe and not self._first_step:
                fock_mat[0] += V_emb
                if self.scf_type != 'restricted':
                    fock_mat[1] += V_emb

            elif self._gostshyp and not self._first_step:
                fock_mat[0] += pr_mat
                if self.scf_type != 'restricted':
                    fock_mat[1] += pr_mat
            
            elif self.point_charges is not None and not self._first_step:
                fock_mat[0] += V_emb
                if self.scf_type != 'restricted':
                    fock_mat[1] += V_emb

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
            for line in self._embedding_drv.get_pe_summary():
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

        if self._cpcm:
            cur_str = 'Solvation Model                 : '
            cur_str += 'C-PCM with ISWIG Discretization'
            self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'C-PCM Dielectric Constant       : '
            cur_str += f'{self.cpcm_drv.epsilon}'
            self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'C-PCM Points per Atomic Sphere  : '
            cur_str += f'{self.cpcm_drv.grid_per_sphere}'
            self.ostream.print_header(cur_str.ljust(str_width))

        if self._gostshyp:
            cur_str = 'Pressure Model                  : '
            cur_str += 'GOSTSHYP'
            self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'Input Pressure                  : '
            cur_str += f'{self._pressure_in_input_units} '
            cur_str += self.pressure_units
            self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'Discretization Choice for Cavity: '
            cur_str += self.discretization
            self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'Tessellation Sphere Scaling     : '
            cur_str += f'{self.tssf}'
            self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'Grid Points per Atomic Sphere   : '
            cur_str += f'{self.num_leb_points}'
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

                if self._gostshyp:
                    valstr = 'Total of {:d} Negative Gaussian Amplitudes'.format(self._gostshyp_drv._neg_p_amp)
                    self.ostream.print_header(valstr.ljust(92))

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
        valstr = 'Multiplicity (2S+1)           :{:3.0f}'.format(mult)
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

        # TODO: extract nuclei-MM energy also for polarizable embedding

        enuc_mm = self._nuc_mm_energy

        e_d4 = self._d4_energy

        e_ef_nuc = self._ef_nuc_energy

        etot = self._iter_data['energy']

        e_el = etot - enuc - enuc_mm - e_d4 - e_ef_nuc
        if self._cpcm:
            e_el -= self.cpcm_epol

        if self._gostshyp:
            e_el -= self._e_gost

        valstr = f'Total Energy                       :{etot:20.10f} a.u.'
        self.ostream.print_header(valstr.ljust(92))

        valstr = f'Electronic Energy                  :{e_el:20.10f} a.u.'
        self.ostream.print_header(valstr.ljust(92))

        if self._cpcm:
            valstr = 'Electrostatic Solvation Energy     :'
            valstr += f'{self.cpcm_epol:20.10f} a.u.'
            self.ostream.print_header(valstr.ljust(92))

        if self._gostshyp:
            valstr = 'GOSTSHYP Pressure Energy           :'
            valstr += f'{self._e_gost:20.10f} a.u.'
            self.ostream.print_header(valstr.ljust(92))

        valstr = f'Nuclear Repulsion Energy           :{enuc:20.10f} a.u.'
        self.ostream.print_header(valstr.ljust(92))

        if self.point_charges is not None:
            valstr = f'Nuclei-Point Charges Energy        :{enuc_mm:20.10f} a.u.'
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
