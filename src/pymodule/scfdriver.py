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

from mpi4py import MPI
from pathlib import Path
from datetime import datetime
import numpy as np
import time as tm
import math
import sys
import os

from .veloxchemlib import AODensityMatrix, denmat
from .veloxchemlib import ScreeningData, GpuDevices
from .veloxchemlib import mpi_master
from .veloxchemlib import (
    compute_fock_gpu,
    matmul_gpu,
    eigh_gpu,
    dot_product_gpu,
    compute_overlap_and_kinetic_energy_integrals_gpu,
    compute_nuclear_potential_integrals_gpu,
)
from .profiler import Profiler
from .molecularbasis import MolecularBasis
from .molecularorbitals import MolecularOrbitals, molorb
from .sadguessdriver import SadGuessDriver
from .diis import Diis
from .inputparser import (parse_input, print_keywords, print_attributes,
                          get_random_string_parallel)
from .sanitychecks import molecule_sanity_check
from .errorhandler import assert_msg_critical


class ScfDriver:
    """
    Implements SCF method with DIIS and two-level DIIS convergence
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
        - first_step: The flag for first step in two-level DIIS convergence
          acceleration.
        - conv_thresh: The SCF convergence threshold.
        - eri_thresh: The electron repulsion integrals screening threshold.
        - ovl_thresh: The atomic orbitals linear dependency threshold.
        - diis_thresh: The DIIS switch on threshold.
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
        - dispersion: The flag for calculating D4 dispersion correction.
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

        self.pair_thresh = 1.0e-10
        self.density_thresh = 1.0e-10

        self.prelink_thresh = 5.0e-6

        # iterations data
        self._history = None
        self._iter_data = None
        self._is_converged = False
        self._scf_energy = 0.0
        self._num_iter = 0

        # density matrix and molecular orbitals
        self._density = None
        self._molecular_orbitals = MolecularOrbitals()

        # nuclear repulsion energy
        self._nuc_energy = 0.0
        self._nuc_mm_energy = 0.0

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

        # for open-shell system: unpaired electrons for initial guess
        # TODO: enable guess_unpaired_electrons input for open-shell
        self.guess_unpaired_electrons = ''

        # scf tensors
        self._scf_tensors = None

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
        self._filename = 'veloxchem_scf_' + get_random_string_parallel(
            self.comm)

        # input keywords
        self._input_keywords = {
            'scf': {
                'acc_type':
                    ('str_upper', 'type of SCF convergence accelerator'),
                'max_iter': ('int', 'maximum number of SCF iterations'),
                'conv_thresh': ('float', 'SCF convergence threshold'),
                'eri_thresh': ('float', 'ERI screening threshold'),
                'pair_thresh': ('float', 'GTO pair screening threshold'),
                'density_thresh': ('float', 'density screening threshold'),
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

    @property
    def filename(self):
        """
        Getter function for protected filename attribute.
        """

        return self._filename

    @filename.setter
    def filename(self, value):
        """
        Setter function for protected filename attribute.
        """

        self._filename = value

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
            self._filename = scf_dict['filename']
            if 'checkpoint_file' not in scf_dict:
                self.checkpoint_file = f'{self._filename}.scf.h5'

        method_keywords = {
            key: val[0]
            for key, val in self._input_keywords['method_settings'].items()
        }

        parse_input(self, method_keywords, method_dict)

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

        if min_basis is None:
            if self.rank == mpi_master():
                min_basis = MolecularBasis.read(molecule,
                                                'AO-START-GUESS',
                                                basis_path='.',
                                                ostream=None)
            else:
                min_basis = MolecularBasis()

        # check molecule
        molecule_sanity_check(molecule)

        # check print level (verbosity of output)
        if self.print_level < 2:
            self.print_level = 1
        if self.print_level > 2:
            self.print_level = 3

        if self.restart:
            self.acc_type = 'DIIS'

        # nuclear repulsion energy
        self._nuc_energy = molecule.nuclear_repulsion_energy()

        if self.rank == mpi_master():
            self._print_header()
            valstr = 'Nuclear repulsion energy: {:.10f} a.u.'.format(
                self._nuc_energy)
            self.ostream.print_info(valstr)
            self.ostream.print_blank()

        # ERI references
        valstr = 'Using the following algorithms for ERI evaluation '
        valstr += 'and Fock matrix formation.'
        self.ostream.print_info(valstr)
        self.ostream.print_blank()
        valstr = 'S. Obara and A. Saika, J. Chem. Phys. 84, 3963-3974 (1986)'
        self.ostream.print_reference(valstr)
        self.ostream.print_blank()
        valstr = 'I. S. Ufimtsev and T. J. Martinez, '
        valstr += 'J. Chem. Theory Comput. 5, 1004-1015 (2009)'
        self.ostream.print_reference(valstr)
        self.ostream.print_blank()
        valstr = 'J. Kussmann and C. Ochsenfeld, '
        valstr += 'J. Chem. Phys. 138, 134114 (2013)'
        self.ostream.print_reference(valstr)
        self.ostream.print_blank()

        # DIIS method
        if self.acc_type.upper() == 'DIIS':

            profiler = Profiler()

            if self.rank == mpi_master():
                den_mat = self.gen_initial_density_sad(
                    molecule, ao_basis, min_basis)
                naos = den_mat.shape[0]
            else:
                naos = None
            naos = self.comm.bcast(naos, root=mpi_master())

            if self.rank != mpi_master():
                den_mat = np.zeros((naos, naos))

            self.comm.Bcast(den_mat, root=mpi_master())

            self._comp_diis(molecule, ao_basis, den_mat, profiler)

            profiler.end(self.ostream, scf_flag=True)

        # two level DIIS method
        elif self.acc_type.upper() == 'L2_DIIS':

            profiler = Profiler()

            # first step
            self._first_step = True

            old_thresh = self.conv_thresh
            self.conv_thresh = 1.0e-3

            old_max_iter = self.max_iter
            self.max_iter = 5

            val_basis = ao_basis.reduce_to_valence_basis()

            if self.rank == mpi_master():
                den_mat = self.gen_initial_density_sad(
                    molecule, val_basis, min_basis)
                naos = den_mat.shape[0]
            else:
                naos = None
            naos = self.comm.bcast(naos, root=mpi_master())

            if self.rank != mpi_master():
                den_mat = np.zeros((naos, naos))

            self.comm.Bcast(den_mat, root=mpi_master())

            self._comp_diis(molecule, val_basis, den_mat, profiler)

            profiler.end(self.ostream, scf_flag=True)

            profiler = Profiler()

            # second step
            self._first_step = False

            self.diis_thresh = 1000.0
            self.conv_thresh = old_thresh
            self.max_iter = old_max_iter

            if self.rank == mpi_master():
                den_mat = self.gen_initial_density_proj(
                    molecule, ao_basis, val_basis, self._molecular_orbitals)
                naos = den_mat.shape[0]
            else:
                naos = None
            naos = self.comm.bcast(naos, root=mpi_master())

            if self.rank != mpi_master():
                den_mat = np.zeros((naos, naos))

            self.comm.Bcast(den_mat, root=mpi_master())

            self._comp_diis(molecule, ao_basis, den_mat, profiler)

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

        mo_1 = valence_mo.alpha_to_numpy()
        mo_2 = np.zeros((nao_2, nmo_1))

        # TODO: unrestricted MO

        row_1 = 0
        row_2 = 0
        for ang in range(valence_basis.max_angular_momentum() + 1):
            for s in range(-ang, ang + 1):
                for a in range(molecule.number_of_atoms()):
                    nbf_1 = valence_basis.number_of_basis_functions([a], ang)
                    nbf_2 = ao_basis.number_of_basis_functions([a], ang)
                    if nbf_1 > 0:
                        mo_2[row_2:row_2 + nbf_2, :] = mo_1[row_1:row_1 +
                                                            nbf_1, :]
                    row_1 += nbf_1
                    row_2 += nbf_2

        proj_mo = MolecularOrbitals(
            [mo_2], [np.zeros(nmo_1)],
            [molecule.get_aufbau_alpha_occupation(nmo_1)],
            valence_mo.get_orbitals_type())

        return proj_mo.get_density(molecule, self.scf_type)

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

    def _get_num_gpus_per_node(self):
        """
        Gets number of GPUs per MPI process.

        :return:
            The number of GPUs per MPI process.
        """

        if 'VLX_NUM_GPUS_PER_NODE' in os.environ:
            num_gpus_per_node = int(os.environ['VLX_NUM_GPUS_PER_NODE'])
        else:
            devices = GpuDevices()
            num_gpus_per_node = devices.get_number_devices()
            if 'SLURM_NTASKS_PER_NODE' in os.environ:
                num_gpus_per_node //= int(os.environ['SLURM_NTASKS_PER_NODE'])

        assert_msg_critical(
            num_gpus_per_node > 0,
            'SCF driver: Invalid number of GPUs per MPI process')

        return num_gpus_per_node

    def _comp_diis(self, molecule, ao_basis, den_mat, profiler):
        """
        Performs SCF calculation with DIIS acceleration.

        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param profiler:
            The profiler.
        """

        self._scf_tensors = None

        self._history = []

        profiler.begin({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        diis_start_time = tm.time()

        if self.rank == mpi_master():
            acc_diis = Diis(self.max_err_vecs, self.diis_thresh)

        t0 = tm.time()

        num_gpus_per_node = self._get_num_gpus_per_node()

        screener = ScreeningData(molecule, ao_basis, num_gpus_per_node,
                                 self.pair_thresh, self.density_thresh)

        self.ostream.print_info(
            f'Using {num_gpus_per_node} GPU devices per MPI rank.')
        self.ostream.print_blank()

        self.ostream.print_info(
            'Screening data computed in {:.2f} sec.'.format(tm.time() - t0))
        self.ostream.print_blank()

        ovl_mat, kin_mat, npot_mat = self._comp_one_ints(
            molecule, ao_basis, screener)

        linear_dependency = False

        if self.rank == mpi_master():
            t0 = tm.time()

            S = ovl_mat

            eigvals, eigvecs_T = eigh_gpu(S, screener.get_num_gpus_per_node())
            eigvecs = eigvecs_T.T

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

        # TODO: double check self.eri_thresh_tight

        if (linear_dependency and self.eri_thresh > self.eri_thresh_tight):
            self.eri_thresh = self.eri_thresh_tight

            if self.rank == mpi_master():
                self.ostream.print_info('ERI screening threshold tightened to' +
                                        ' {:.1e}.'.format(self.eri_thresh))
                self.ostream.print_blank()

        self._density = den_mat.copy()

        profiler.check_memory_usage('Initial guess')

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

            dmat = AODensityMatrix([den_mat], denmat.rest)

            fock_mat = self._comp_2e_fock(dmat, molecule, ao_basis, screener,
                                          e_grad, profiler)

            profiler.start_timer('ErrVec')

            e_el = self._comp_energy(fock_mat, kin_mat, npot_mat, den_mat)

            fock_mat = self._comp_full_fock(fock_mat, kin_mat, npot_mat)

            e_mat, e_grad, max_grad = self._comp_gradient(
                fock_mat, ovl_mat, den_mat, oao_mat)

            # compute density change and energy change

            diff_den = self._comp_density_change(den_mat, self._density)

            e_scf = (e_el + self._nuc_energy)

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

            self._density = den_mat.copy()

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

            if self.rank == mpi_master():
                acc_diis.store_diis_data(fock_mat, e_mat, e_grad)

            if self.rank == mpi_master():
                eff_fock_mat = acc_diis.get_effective_fock(
                    fock_mat, self.ostream)
            else:
                eff_fock_mat = None

            self._molecular_orbitals = self._gen_molecular_orbitals(
                molecule, eff_fock_mat, oao_mat,
                screener.get_num_gpus_per_node())

            if self._mom is not None:
                self._apply_mom(molecule, ovl_mat)

            self._update_mol_orbs_phase()

            if self.rank == mpi_master():
                # TODO unrestricted and restricted open-shell
                mo = self.molecular_orbitals.alpha_to_numpy()
                occ = self.molecular_orbitals.occa_to_numpy()
                occ_mo = occ * mo
                den_mat = matmul_gpu(occ_mo, occ_mo.T)
                naos = den_mat.shape[0]
            else:
                naos = None
            naos = self.comm.bcast(naos, root=mpi_master())

            if self.rank != mpi_master():
                den_mat = np.zeros((naos, naos))

            self.comm.Bcast(den_mat, root=mpi_master())

            profiler.stop_timer('FockDiag')

            profiler.check_memory_usage('Iteration {:d} Fock diag.'.format(
                self._num_iter))

        if self.rank == mpi_master():
            acc_diis.clear()
            self._V_es = None

        if (not self._first_step) and self.is_converged:

            if self.rank == mpi_master():
                S = ovl_mat

                C_alpha = self.molecular_orbitals.alpha_to_numpy()
                C_beta = self.molecular_orbitals.beta_to_numpy()

                E_alpha = self.molecular_orbitals.ea_to_numpy()
                E_beta = self.molecular_orbitals.eb_to_numpy()

                occ_alpha = self.molecular_orbitals.occa_to_numpy()
                occ_beta = self.molecular_orbitals.occb_to_numpy()

                if self.scf_type == 'restricted':
                    D_alpha = self._density
                    D_beta = self._density
                else:
                    D_alpha = self._density[0]
                    D_beta = self._density[1]

                F_alpha = fock_mat
                # F_beta = fock_mat.beta_to_numpy(0)

                # TODO: add F_beta

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
                }

            else:
                self._scf_tensors = None

        if self.rank == mpi_master():
            self._print_scf_finish(diis_start_time)

        profiler.check_memory_usage('End of SCF')

    def _comp_one_ints(self, molecule, basis, screener):
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

        t0 = tm.time()

        ovl_mat, kin_mat = compute_overlap_and_kinetic_energy_integrals_gpu(
            molecule, basis, screener)

        npot_mat = compute_nuclear_potential_integrals_gpu(
            molecule, basis, screener)

        naos = ovl_mat.number_of_rows()

        if self.rank == mpi_master():
            ovl_mat_np = np.zeros((naos, naos))
            kin_mat_np = np.zeros((naos, naos))
            npot_mat_np = np.zeros((naos, naos))
        else:
            ovl_mat_np = None
            kin_mat_np = None
            npot_mat_np = None

        self.comm.Reduce(ovl_mat.to_numpy(),
                         ovl_mat_np,
                         op=MPI.SUM,
                         root=mpi_master())
        self.comm.Reduce(kin_mat.to_numpy(),
                         kin_mat_np,
                         op=MPI.SUM,
                         root=mpi_master())
        self.comm.Reduce(npot_mat.to_numpy(),
                         npot_mat_np,
                         op=MPI.SUM,
                         root=mpi_master())

        ovl_mat = ovl_mat_np
        kin_mat = kin_mat_np
        npot_mat = npot_mat_np

        t1 = tm.time()

        if self.rank == mpi_master() and self.print_level > 1:
            self.ostream.print_info(
                'One-electron integral matrices computed in' +
                ' {:.2f} sec.'.format(t1 - t0))
            self.ostream.print_blank()
            self.ostream.flush()

        return ovl_mat, kin_mat, npot_mat

    def _comp_2e_fock(self,
                      dmat,
                      molecule,
                      basis,
                      screener,
                      e_grad=None,
                      profiler=None):
        """
        Computes Fock/Kohn-Sham matrix (only 2e part).

        :param dmat:
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
            The AO Kohn-Sham (Vxc) matrix.
        """

        # TODO: use dynamic threshold

        # if e_grad is None:
        #     thresh_int = 12
        # else:
        #     thresh_int = int(-math.log10(self._get_dyn_threshold(e_grad)))

        eri_t0 = tm.time()
        coulomb_timing = 0.0
        exchange_timing = 0.0

        # Hartree-Fock
        fock_mat = compute_fock_gpu(molecule, basis, dmat, 2.0, 1.0, 0.0,
                                    'symm', self.eri_thresh,
                                    self.prelink_thresh, screener)
        fock_mat_local = fock_mat.to_numpy()

        coulomb_timing += np.array(
            [float(dt.split()[0]) for dt in screener.get_coulomb_time()])

        exchange_timing += np.array(
            [float(dt.split()[0]) for dt in screener.get_exchange_time()])

        all_coulomb_timing = self.comm.allgather(coulomb_timing)
        all_exchange_timing = self.comm.allgather(exchange_timing)

        all_coulomb_timing = np.array(all_coulomb_timing).reshape(-1)
        all_exchange_timing = np.array(all_exchange_timing).reshape(-1)

        max_coulomb_timing = np.max(all_coulomb_timing)
        max_exchange_timing = np.max(all_exchange_timing)

        coulomb_load_imb = 1.0 - np.sum(all_coulomb_timing) / (
            all_coulomb_timing.size * max_coulomb_timing)
        exchange_load_imb = 1.0 - np.sum(all_exchange_timing) / (
            all_exchange_timing.size * max_exchange_timing)

        if self.rank == mpi_master():
            fock_mat = np.zeros(fock_mat_local.shape)
        else:
            fock_mat = None

        self.comm.Reduce(fock_mat_local,
                         fock_mat,
                         op=MPI.SUM,
                         root=mpi_master())

        # TODO: add beta density

        if self.timing:
            eri_dt = tm.time() - eri_t0
            profiler.add_timing_info(
                'FockPrep',
                max(0.0, eri_dt - max_coulomb_timing - max_exchange_timing))
            profiler.add_timing_info('FockJ', max_coulomb_timing)
            profiler.add_timing_info('LoadImbJ', coulomb_load_imb)
            profiler.add_timing_info('FockK', max_exchange_timing)
            profiler.add_timing_info('LoadImbK', exchange_load_imb)

        return fock_mat

    def _comp_energy(self, fock_mat, kin_mat, npot_mat, den_mat):
        """
        Computes the sum of SCF energy components: electronic energy, kinetic
        energy, and nuclear potential energy.

        :param fock_mat:
            The Fock/Kohn-Sham matrix (only 2e-part).
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

        if self.rank == mpi_master():
            hcore_mat = kin_mat - npot_mat
            hcore_plus_full_Fock = 2.0 * hcore_mat + fock_mat
            e_sum = dot_product_gpu(den_mat, hcore_plus_full_Fock)
        else:
            e_sum = 0.0

        e_sum = self.comm.allreduce(e_sum, op=MPI.SUM)

        return e_sum

    def _comp_full_fock(self, fock_mat, kin_mat, npot_mat):
        """
        Computes full Fock/Kohn-Sham matrix by adding to 2e-part of
        Fock/Kohn-Sham matrix the kinetic energy and nuclear potential
        matrices.

        :param fock_mat:
            The Fock/Kohn-Sham matrix (2e-part).
        :param kin_mat:
            The kinetic energy matrix.
        :param npot_mat:
            The nuclear potential matrix.
        """

        if self.rank == mpi_master():
            fock_mat += (kin_mat - npot_mat)

        return fock_mat

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

    def _gen_molecular_orbitals(self, molecule, fock_mat, oao_mat,
                                num_gpus_per_node):
        """
        Generates molecular orbital by diagonalizing Fock/Kohn-Sham matrix.

        :param molecule:
            The molecule.
        :param fock_mat:
            The Fock/Kohn-Sham matrix.
        :param oao_mat:
            The orthogonalization matrix.
        :param num_gpus_per_node:
            The number of GPUs per MPI process.

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
                D_alpha = self._density
                D_beta = self._density
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

        self.ostream.print_blank()

    def _print_scf_title(self):
        """
        Prints SCF cycles header to output stream.
        """

        if self._first_step:
            self.ostream.print_info('Starting Reduced Basis SCF calculation...')

        else:
            self.ostream.print_blank()
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

        ovl_a_b = matmul_gpu(Cocc_a.T, matmul_gpu(smat, Cocc_b))
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

        etot = self._iter_data['energy']

        e_el = etot - enuc

        valstr = f'Total Energy                       :{etot:20.10f} a.u.'
        self.ostream.print_header(valstr.ljust(92))

        valstr = f'Electronic Energy                  :{e_el:20.10f} a.u.'
        self.ostream.print_header(valstr.ljust(92))

        valstr = f'Nuclear Repulsion Energy           :{enuc:20.10f} a.u.'
        self.ostream.print_header(valstr.ljust(92))

        self.ostream.print_header(
            '------------------------------------'.ljust(92))

        grad = self._iter_data['gradient_norm']
        valstr = 'Gradient Norm                      :{:20.10f} a.u.'.format(
            grad)
        self.ostream.print_header(valstr.ljust(92))

        self.ostream.print_blank()
