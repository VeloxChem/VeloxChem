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
from pathlib import Path
from copy import deepcopy
import numpy as np
import time as tm
import math
import sys

from .veloxchemlib import XCFunctional, MolecularGrid
from .veloxchemlib import mpi_master, rotatory_strength_in_cgs
from .veloxchemlib import denmat
from .aodensitymatrix import AODensityMatrix
from .outputstream import OutputStream
from .profiler import Profiler
from .linearsolver import LinearSolver
from .blockdavidson import BlockDavidsonSolver
from .molecularorbitals import MolecularOrbitals
from .visualizationdriver import VisualizationDriver
from .cubicgrid import CubicGrid
from .oneeints import (compute_electric_dipole_integrals,
                       compute_linear_momentum_integrals,
                       compute_angular_momentum_integrals)
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           dft_sanity_check, pe_sanity_check,
                           solvation_model_sanity_check)
from .errorhandler import assert_msg_critical
from .checkpoint import (read_rsp_hdf5, write_rsp_hdf5, write_rsp_solution,
                         write_lr_rsp_results_to_hdf5,
                         write_detach_attach_to_hdf5)


class TdaEigenSolver(LinearSolver):
    """
    Implements TDA excited states computation schheme for Hartree-Fock/Kohn-Sham
    level of theory.

    # vlxtag: RHF, Absorption, ECD, CIS
    # vlxtag: RKS, Absorption, ECD, TDA

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - nstates: The number of excited states to be determined.
        - core_excitation: The flag for computing core-excited states.
        - num_core_orbitals: The number of core-orbitals involved in
          the core-excited states.
        - solver: The eigenvalues solver.
        - nto: The flag for natural transition orbital analysis.
        - nto_pairs: The number of NTO pairs in NTO analysis.
        - detach_attach: The flag for detachment/attachment density analysis.
        - cube_origin: The origin of cubic grid points.
        - cube_stepsize: The step size of cubic grid points in X, Y and Z
          directions.
        - cube_points: The number of cubic grid points in X, Y and Z directions.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes TDA excited states computation drived to default setup.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        super().__init__(comm, ostream)

        # excited states information
        self.nstates = 3

        self.core_excitation = False
        self.num_core_orbitals = 0

        # restricted subspace
        self.rsa = False
        self.num_vir_orbitals = 0
        self.num_val_orbitals = 0

        # solver setup
        self.solver = None

        # NTO and detachment/attachment density
        self.nto = False
        self.nto_pairs = None
        self.nto_cubes = False
        self.detach_attach = False
        self.detach_attach_cubes = False
        self.cube_origin = None
        self.cube_stepsize = None
        self.cube_points = [80, 80, 80]

        self._input_keywords['response'].update({
            'nstates': ('int', 'number of excited states'),
            'core_excitation': ('bool', 'compute core-excited states'),
            'num_core_orbitals': ('int', 'number of involved core-orbitals'),
            'nto': ('bool', 'analyze natural transition orbitals'),
            'nto_pairs': ('int', 'number of NTO pairs in NTO analysis'),
            'nto_cubes': ('bool', 'write NTO cube files'),
            'detach_attach': ('bool', 'analyze detachment/attachment density'),
            'detach_attach_cubes':
                ('bool', 'write detachment/attachment density cube files'),
            'cube_origin': ('seq_fixed', 'origin of cubic grid points'),
            'cube_stepsize': ('seq_fixed', 'step size of cubic grid points'),
            'cube_points': ('seq_fixed_int', 'number of cubic grid points'),
        })

        self._input_keywords['response'].pop('lindep_thresh', None)

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in TDA excited states computation.

        :param rsp_dict:
            The dictionary of response input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}

        super().update_settings(rsp_dict, method_dict)

        if self.cube_origin is not None:
            assert_msg_critical(
                len(self.cube_origin) == 3,
                'TdaEigenSolver: cube origin needs 3 numbers')

        if self.cube_stepsize is not None:
            assert_msg_critical(
                len(self.cube_stepsize) == 3,
                'TdaEigenSolver: cube stepsize needs 3 numbers')

        if self.cube_points is not None:
            assert_msg_critical(
                len(self.cube_points) == 3,
                'TdaEigenSolver: cube points needs 3 integers')

        # If the detachemnt and attachment cube files are requested,
        # set the detach_attach flag to True to get the detachment and
        # attachment densities.
        if self.detach_attach_cubes:
            self.detach_attach = True

    def compute(self, molecule, basis, scf_tensors):
        """
        Performs TDA excited states calculation using molecular data.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            A dictionary containing eigenvalues, eigenvectors, transition
            dipole moments, oscillator strengths and rotatory strengths.
        """

        # check molecule
        molecule_sanity_check(molecule)

        # check SCF results
        scf_results_sanity_check(self, scf_tensors)

        # update checkpoint_file after scf_results_sanity_check
        if self.filename is not None and self.checkpoint_file is None:
            self.checkpoint_file = f'{self.filename}_rsp.h5'

        # check dft setup
        dft_sanity_check(self, 'compute')

        # check pe setup
        pe_sanity_check(self, molecule=molecule)

        # check solvation setup
        solvation_model_sanity_check(self)

        # check print level (verbosity of output)
        if self.print_level < 2:
            self.print_level = 1
        if self.print_level > 2:
            self.print_level = 3

        # initialize profiler
        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        if self.rank == mpi_master():
            self._print_header('TDA Driver', nstates=self.nstates)

        # set start time

        self.start_time = tm.time()

        # sanity check

        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'TdaEigenSolver: not implemented for unrestricted case')

        # prepare molecular orbitals

        if self.rank == mpi_master():
            orb_ene = scf_tensors['E_alpha']
            norb = orb_ene.shape[0]
            nocc = molecule.number_of_alpha_electrons()
            if self.core_excitation:
                assert_msg_critical(
                    self.nstates <= self.num_core_orbitals * (norb - nocc),
                    'TdaEigenSolver: too many excited states')
            elif self.rsa:
                assert_msg_critical(
                    self.nstates <= (self.num_core_orbitals + self.num_val_orbitals) * self.num_vir_orbitals,
                    'TdaEigenSolver: too many excited states'
                )
            else:
                assert_msg_critical(self.nstates <= nocc * (norb - nocc),
                                    'TdaEigenSolver: too many excited states')
        else:
            orb_ene = None
            nocc = None

        # ERI information
        eri_dict = self._init_eri(molecule, basis)

        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors)

        # PE information
        pe_dict = self._init_pe(molecule, basis)

        # CPCM_information
        self._init_cpcm(molecule)

        # set up trial excitation vectors on master node

        diag_mat, trial_mat = self._gen_trial_vectors(molecule, orb_ene, nocc)

        # block Davidson algorithm setup

        self.solver = BlockDavidsonSolver()

        # read initial guess from restart file

        n_restart_vectors = 0
        n_restart_iterations = 0

        if self.restart:
            if self.rank == mpi_master():
                rst_trial_mat, rst_sig_mat = read_rsp_hdf5(
                    self.checkpoint_file, ['TDA_trials', 'TDA_sigmas'],
                    molecule, basis, dft_dict, pe_dict, self.ostream)
                self.restart = (rst_trial_mat is not None and
                                rst_sig_mat is not None)
                if rst_trial_mat is not None:
                    n_restart_vectors = rst_trial_mat.shape[1]
            self.restart = self.comm.bcast(self.restart, root=mpi_master())
            n_restart_vectors = self.comm.bcast(n_restart_vectors,
                                                root=mpi_master())
            n_restart_iterations = n_restart_vectors // self.nstates
            if n_restart_vectors % self.nstates != 0:
                n_restart_iterations += 1

        profiler.check_memory_usage('Initial guess')

        # start TDA iteration

        for i in range(self.max_iter):

            profiler.set_timing_key(f'Iteration {i + 1}')

            # perform linear transformation of trial vectors

            if i >= n_restart_iterations:
                tdens = self._get_trans_densities(trial_mat, scf_tensors,
                                                  molecule)
                fock = self._comp_lr_fock(tdens, molecule, basis, eri_dict,
                                          dft_dict, pe_dict, profiler)

            profiler.start_timer('ReducedSpace')

            # solve eigenvalues problem on master node

            if self.rank == mpi_master():

                if i >= n_restart_iterations:
                    sig_mat = self._get_sigmas(fock, scf_tensors, molecule,
                                               trial_mat)
                else:
                    istart = i * self.nstates
                    iend = (i + 1) * self.nstates
                    if iend > n_restart_vectors:
                        iend = n_restart_vectors
                    sig_mat = np.copy(rst_sig_mat[:, istart:iend])
                    trial_mat = np.copy(rst_trial_mat[:, istart:iend])

                self.solver.add_iteration_data(sig_mat, trial_mat, i)

                trial_mat = self.solver.compute(diag_mat)

                self._print_iter_data(i)

            profiler.stop_timer('ReducedSpace')

            profiler.check_memory_usage(f'Iteration {i + 1}')

            profiler.print_memory_tracing(self.ostream)

            self._cur_iter = i

            # check convergence

            self._check_convergence()

            # write checkpoint file

            if (self.rank == mpi_master() and i >= n_restart_iterations):
                trials = self.solver.trial_matrices
                sigmas = self.solver.sigma_matrices
                write_rsp_hdf5(self.checkpoint_file, [trials, sigmas],
                               ['TDA_trials', 'TDA_sigmas'], molecule, basis,
                               dft_dict, pe_dict, self.ostream)

            # finish TDA after convergence

            if self._is_converged:
                break

        # converged?
        if self.rank == mpi_master():
            self._print_convergence('{:d} excited states'.format(self.nstates))

            # final hdf5 file to save response results
            if self.filename is not None:
                final_h5_fname = f'{self.filename}.h5'
            else:
                final_h5_fname = None

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of TDA eigensolver')
        profiler.print_memory_usage(self.ostream)

        # compute 1e dipole integrals

        integrals = self._comp_onee_integrals(molecule, basis)

        # print converged excited states

        if self.rank == mpi_master() and self._is_converged:
            if self.core_excitation:
                mo_occ = scf_tensors['C_alpha'][:, :self.
                                                num_core_orbitals].copy()
                mo_vir = scf_tensors['C_alpha'][:, nocc:].copy()
            elif self.rsa:
                mo_occ = np.hstack((scf_tensors['C_alpha'][:, :self.num_core_orbitals].copy(),
                                scf_tensors['C_alpha'][:, nocc - self.num_val_orbitals:nocc].copy()))
                mo_vir = scf_tensors['C_alpha'][:, nocc:nocc+self.num_vir_orbitals].copy()
            else:
                mo_occ = scf_tensors['C_alpha'][:, :nocc].copy()
                mo_vir = scf_tensors['C_alpha'][:, nocc:].copy()

            eigvals, rnorms = self.solver.get_eigenvalues()
            eigvecs = self.solver.ritz_vectors

            trans_dipoles = self._comp_trans_dipoles(integrals, eigvals,
                                                     eigvecs, mo_occ, mo_vir)

            oscillator_strengths = (2.0 / 3.0) * np.sum(
                trans_dipoles['electric']**2, axis=1) * eigvals

            rotatory_strengths = np.sum(
                trans_dipoles['velocity'] * trans_dipoles['magnetic'],
                axis=1) * rotatory_strength_in_cgs()

        # natural transition orbitals and detachment/attachment densities

        nto_lambdas = []
        nto_cube_files = []
        dens_cube_files = []

        excitation_details = []

        for s in range(self.nstates):
            if self.rank == mpi_master():
                t_mat = eigvecs[:, s].reshape(mo_occ.shape[1], mo_vir.shape[1])

                # save excitation details
                excitation_details.append(
                    self.get_excitation_details(eigvecs[:, s], mo_occ.shape[1],
                                                mo_vir.shape[1]))

            if self.nto or self.detach_attach:
                vis_drv = VisualizationDriver(self.comm)
                if self.cube_origin is None or self.cube_stepsize is None:
                    cubic_grid = vis_drv.gen_cubic_grid(molecule,
                                                        self.cube_points)
                else:
                    cubic_grid = CubicGrid(self.cube_origin, self.cube_stepsize,
                                           self.cube_points)

            if self.nto and self._is_converged:
                self.ostream.print_info(
                    'Running NTO analysis for S{:d}...'.format(s + 1))
                self.ostream.flush()

                if self.rank == mpi_master():
                    nto_mo = self.get_nto(t_mat, mo_occ, mo_vir)

                    nto_lam = nto_mo.occa_to_numpy()
                    lam_start = mo_occ.shape[1]
                    lam_end = lam_start + min(mo_occ.shape[1], mo_vir.shape[1])
                    nto_lambdas.append(nto_lam[lam_start:lam_end])

                    # Add the NTO to the final hdf5 file.
                    nto_label = f'NTO_S{s + 1}'
                    if final_h5_fname is not None:
                        nto_mo.write_hdf5(final_h5_fname,
                                          label=f'rsp/{nto_label}')
                else:
                    nto_mo = MolecularOrbitals()
                nto_mo = nto_mo.broadcast(self.comm, root=mpi_master())

                if self.nto_cubes:
                    lam_diag, nto_cube_fnames = self.write_nto_cubes(
                        cubic_grid, molecule, basis, s, nto_mo, self.nto_pairs)

                    if self.rank == mpi_master():
                        nto_cube_files.append(nto_cube_fnames)

            if self.detach_attach and self._is_converged:
                self.ostream.print_info(
                    'Running detachment/attachment analysis for S{:d}...'.
                    format(s + 1))
                self.ostream.flush()

                if self.rank == mpi_master():
                    dens_D, dens_A = self.get_detach_attach_densities(
                        t_mat, None, mo_occ, mo_vir)

                    # Add the detachment and attachment density matrices
                    # to the final hdf5 file
                    state_label = f'S{s + 1}'
                    if final_h5_fname is not None:
                        write_detach_attach_to_hdf5(final_h5_fname, state_label,
                                                    dens_D, dens_A)

                # Generate and save cube files
                if self.detach_attach_cubes:
                    if self.rank == mpi_master():
                        dens_DA = AODensityMatrix([dens_D, dens_A], denmat.rest)
                    else:
                        dens_DA = AODensityMatrix()
                    dens_DA = dens_DA.broadcast(self.comm, root=mpi_master())

                    dens_cube_fnames = self.write_detach_attach_cubes(
                        cubic_grid, molecule, basis, s, dens_DA)

                    if self.rank == mpi_master():
                        dens_cube_files.append(dens_cube_fnames)

        if (self.nto or self.detach_attach) and self._is_converged:
            self.ostream.print_blank()
            self.ostream.flush()

        # results

        if self.rank == mpi_master() and self._is_converged:

            if self.rsa:
                orbital_details = {
                    'nstates': self.nstates,
                    'num_core': self.num_core_orbitals,
                    'num_val': self.num_val_orbitals,
                    'num_vir': self.num_vir_orbitals
                }
            else:
                orbital_details = {
                    'nstates': self.nstates,
                    'num_core': self.num_core_orbitals,
                    'num_val': nocc,
                    'num_vir': norb - nocc
                }

            ret_dict = {
                'eigenvalues': eigvals,
                'eigenvectors': eigvecs,
                'electric_transition_dipoles': trans_dipoles['electric'],
                'velocity_transition_dipoles': trans_dipoles['velocity'],
                'magnetic_transition_dipoles': trans_dipoles['magnetic'],
                'oscillator_strengths': oscillator_strengths,
                'rotatory_strengths': rotatory_strengths,
                'excitation_details': excitation_details,
                'num_core': orbital_details['num_core'],
                'num_val': orbital_details['num_val'],
                'num_vir': orbital_details['num_vir'],
            }

            if self.nto:
                ret_dict['nto_lambdas'] = nto_lambdas
                if self.nto_cubes:
                    ret_dict['nto_cubes'] = nto_cube_files

            if self.detach_attach_cubes:
                ret_dict['density_cubes'] = dens_cube_files

            self._write_final_hdf5(final_h5_fname, molecule, basis,
                                   dft_dict['dft_func_label'],
                                   pe_dict['potfile_text'], eigvecs)

            if (self.save_solutions and final_h5_fname is not None):
                # Write response results to final checkpoint file.
                write_lr_rsp_results_to_hdf5(final_h5_fname, ret_dict)

            self._print_results(ret_dict)

            return ret_dict
        else:
            return None

    def _gen_trial_vectors(self, molecule, orb_ene, nocc):
        """
        Generates set of TDA trial vectors for given number of excited states
        by selecting primitive excitations wirh lowest approximate energies
        E_ai = e_a-e_i.

        :param molecule:
            The molecule.
        :param orb_ene:
            The orbital energies.
        :param nocc:
            The number of occupied orbitals.

        :return:
            tuple (approximate diagonal of symmetric A, set of trial vectors).
        """

        if self.rank == mpi_master():
            norb = orb_ene.shape[0]
            ea = orb_ene

            if self.core_excitation:
                excitations = [(i, a)
                               for i in range(self.num_core_orbitals)
                               for a in range(nocc, norb)]
                n_exc = self.num_core_orbitals * (norb - nocc)
            elif self.rsa:
                core_and_val_indices = list(range(self.num_core_orbitals)) + list(range(nocc - self.num_val_orbitals, nocc, 1)) #list(range(nocc+1, nocc+1 - self.num_val_orbitals, -1)) 
                excitations = [(i, a)
                               for i in core_and_val_indices
                               for a in range(nocc, nocc+self.num_vir_orbitals)]
                n_exc = (self.num_core_orbitals + self.num_val_orbitals) * self.num_vir_orbitals
            else:
                excitations = [
                    (i, a) for i in range(nocc) for a in range(nocc, norb)
                ]
                n_exc = nocc * (norb - nocc)

            excitation_energies = [ea[a] - ea[i] for i, a in excitations]

            w = {ia: w for ia, w in zip(excitations, excitation_energies)}

            diag_mat = np.array(excitation_energies)
            trial_mat = np.zeros((n_exc, self.nstates))

            for s, (i, a) in enumerate(sorted(w, key=w.get)[:self.nstates]):
                if self.rank == mpi_master():
                    ia = excitations.index((i, a))
                    trial_mat[ia, s] = 1.0

            return diag_mat, trial_mat

        return None, None

    def _check_convergence(self):
        """
        Checks convergence of excitation energies and set convergence flag on
        all processes within MPI communicator.

        :param iteration:
            The current excited states solver iteration.
        """

        self._is_converged = False

        if self.rank == mpi_master():
            self._is_converged = self.solver.check_convergence(self.conv_thresh)

        self._is_converged = self.comm.bcast(self._is_converged,
                                             root=mpi_master())

    def _get_trans_densities(self, trial_mat, tensors, molecule):
        """
        Computes the transition densities.

        :param trial_mat:
            The matrix containing the Z vectors as columns.
        :param tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param molecule:
            The molecule.

        :return:
            The transition density matrix.
        """

        # form transition densities

        if self.rank == mpi_master():
            nocc = molecule.number_of_alpha_electrons()
            norb = tensors['C_alpha'].shape[1]
            nvir = norb - nocc

            if self.core_excitation:
                mo_occ = tensors['C_alpha'][:, :self.num_core_orbitals].copy()
                mo_vir = tensors['C_alpha'][:, nocc:].copy()
            elif self.rsa:
                #mo_occ = tensors['C_alpha'][:, :self.num_core_orbitals].copy()
                mo_occ = np.hstack((tensors['C_alpha'][:, :self.num_core_orbitals].copy(),
                                    tensors['C_alpha'][:, nocc - self.num_val_orbitals:nocc].copy()))
                mo_vir = tensors['C_alpha'][:, nocc:nocc + self.num_vir_orbitals].copy() #mo_vir = tensors['C_alpha'][:, nocc:norb - self.num_vir_orbitals].copy()
            else:
                mo_occ = tensors['C_alpha'][:, :nocc].copy()
                mo_vir = tensors['C_alpha'][:, nocc:].copy()

            tdens = []
            for k in range(trial_mat.shape[1]):
                if self.core_excitation:
                    mat = trial_mat[:, k].reshape(self.num_core_orbitals, nvir)
                    mat = np.matmul(mo_occ, np.matmul(mat, mo_vir.T))
                elif self.rsa:
                    mat = trial_mat[:, k].reshape(self.num_core_orbitals + self.num_val_orbitals, 
                                                  self.num_vir_orbitals)
                    mat = np.matmul(mo_occ, np.matmul(mat, mo_vir.T))
                else:
                    mat = trial_mat[:, k].reshape(nocc, nvir)
                    mat = np.matmul(mo_occ, np.matmul(mat, mo_vir.T))
                tdens.append(mat)
        else:
            tdens = None

        return tdens

    def _get_sigmas(self, fock, tensors, molecule, trial_mat):
        """
        Computes the sigma vectors.

        :param fock:
            The Fock matrix.
        :param tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param molecule:
            The molecule.
        :param trial_mat:
            The trial vectors as 2D Numpy array.

        :return:
            The sigma vectors as 2D Numpy array.
        """

        nocc = molecule.number_of_alpha_electrons()
        norb = tensors['C_alpha'].shape[1]
        nvir = norb - nocc

        if self.core_excitation:
            mo_occ = tensors['C_alpha'][:, :self.num_core_orbitals].copy()
            mo_vir = tensors['C_alpha'][:, nocc:].copy()
        elif self.rsa:
            mo_occ = np.hstack((tensors['C_alpha'][:, :self.num_core_orbitals].copy(),
                                tensors['C_alpha'][:, nocc - self.num_val_orbitals:nocc].copy()))
            mo_vir = tensors['C_alpha'][:, nocc:nocc+self.num_vir_orbitals].copy()
        else:
            mo_occ = tensors['C_alpha'][:, :nocc].copy()
            mo_vir = tensors['C_alpha'][:, nocc:].copy()
        orb_ene = tensors['E_alpha']

        sigma_vecs = []
        for fockind in range(len(fock)):
            # 2e contribution
            mat = fock[fockind].copy()
            mat = np.matmul(mo_occ.T, np.matmul(mat, mo_vir))
            # 1e contribution
            if self.core_excitation:
                cjb = trial_mat[:, fockind].reshape(self.num_core_orbitals,
                                                    nvir)
                mat += np.matmul(cjb, np.diag(orb_ene[nocc:]).T)
                mat -= np.matmul(np.diag(orb_ene[:self.num_core_orbitals]), cjb)
                sigma_vecs.append(mat.reshape(self.num_core_orbitals * nvir, 1))
            elif self.rsa:
                cjb = trial_mat[:, fockind].reshape(self.num_core_orbitals + self.num_val_orbitals,
                                                    self.num_vir_orbitals)
                mat += np.matmul(cjb, np.diag(orb_ene[nocc:nocc+self.num_vir_orbitals]).T)
                mat -= np.matmul(np.diag(np.concatenate((orb_ene[:self.num_core_orbitals],
                                                        orb_ene[nocc - self.num_val_orbitals:nocc]))), cjb)
                sigma_vecs.append(mat.reshape((self.num_core_orbitals + self.num_val_orbitals) * self.num_vir_orbitals, 1))
            else:
                cjb = trial_mat[:, fockind].reshape(nocc, nvir)
                mat += np.matmul(cjb, np.diag(orb_ene[nocc:]).T)
                mat -= np.matmul(np.diag(orb_ene[:nocc]), cjb)
                sigma_vecs.append(mat.reshape(nocc * nvir, 1))

        sigma_mat = sigma_vecs[0]
        for vec in sigma_vecs[1:]:
            sigma_mat = np.hstack((sigma_mat, vec))

        return sigma_mat

    def _comp_onee_integrals(self, molecule, basis):
        """
        Computes one-electron integrals.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.

        :return:
            The one-electron integrals.
        """

        dipole_mats = compute_electric_dipole_integrals(molecule, basis,
                                                        [0.0, 0.0, 0.0])
        linmom_mats = compute_linear_momentum_integrals(molecule, basis)
        angmom_mats = compute_angular_momentum_integrals(
            molecule, basis, [0.0, 0.0, 0.0])

        integrals = {}

        if self.rank == mpi_master():
            integrals['electric_dipole'] = (
                dipole_mats[0],
                dipole_mats[1],
                dipole_mats[2],
            )

            integrals['linear_momentum'] = (
                -1.0 * linmom_mats[0],
                -1.0 * linmom_mats[1],
                -1.0 * linmom_mats[2],
            )

            integrals['angular_momentum'] = (
                -1.0 * angmom_mats[0],
                -1.0 * angmom_mats[1],
                -1.0 * angmom_mats[2],
            )

            integrals['magnetic_dipole'] = (
                0.5 * angmom_mats[0],
                0.5 * angmom_mats[1],
                0.5 * angmom_mats[2],
            )

        return integrals

    def _comp_trans_dipoles(self, integrals, eigvals, eigvecs, mo_occ, mo_vir):
        """
        Computes transition dipole moments.

        :param integrals:
            The one-electron integrals.
        :param eigvals:
            The eigenvalues.
        :param eigvecs:
            The CI vectors.
        :param mo_occ:
            The occupied MO coefficients.
        :param mo_vir:
            The virtual MO coefficients.

        :return:
            The transition dipole moments.
        """

        transition_dipoles = {
            'electric': np.zeros((self.nstates, 3)),
            'velocity': np.zeros((self.nstates, 3)),
            'magnetic': np.zeros((self.nstates, 3))
        }

        sqrt_2 = math.sqrt(2.0)

        for s in range(self.nstates):
            exc_vec = eigvecs[:, s].reshape(mo_occ.shape[1], mo_vir.shape[1])
            trans_dens = sqrt_2 * np.linalg.multi_dot(
                [mo_occ, exc_vec, mo_vir.T])

            transition_dipoles['electric'][s, :] = np.array([
                np.vdot(trans_dens, integrals['electric_dipole'][d].T)
                for d in range(3)
            ])

            transition_dipoles['velocity'][s, :] = np.array([
                np.vdot(trans_dens, integrals['linear_momentum'][d].T)
                for d in range(3)
            ]) / (-eigvals[s])

            transition_dipoles['magnetic'][s, :] = np.array([
                np.vdot(trans_dens, integrals['magnetic_dipole'][d].T)
                for d in range(3)
            ])

        return transition_dipoles

    def _print_iter_data(self, iteration):
        """
        Prints excited states solver iteration data to output stream.

        :param iteration:
            The current excited states solver iteration.
        """

        # iteration header

        exec_str = ' *** Iteration: ' + (str(iteration + 1)).rjust(3)
        exec_str += ' * Reduced Space: '
        exec_str += (str(self.solver.reduced_space_size())).rjust(4)
        rmax, rmin = self.solver.max_min_residual_norms()
        exec_str += ' * Residues (Max,Min): {:.2e} and {:.2e}'.format(
            rmax, rmin)
        self.ostream.print_header(exec_str)
        self.ostream.print_blank()

        # excited states information

        reigs, rnorms = self.solver.get_eigenvalues()
        for i in range(reigs.shape[0]):
            exec_str = 'State {:2d}: {:5.8f} '.format(i + 1, reigs[i])
            exec_str += 'a.u. Residual Norm: {:3.8f}'.format(rnorms[i])
            self.ostream.print_header(exec_str.ljust(84))

        # flush output stream
        self.ostream.print_blank()
        self.ostream.flush()

    def _write_final_hdf5(self, final_h5_fname, molecule, basis, dft_func_label,
                          potfile_text, eigvecs):
        """
        Writes final HDF5 that contains TDA solution vectors.

        :param final_h5_fname:
            The name of the final hdf5 file.
        :param molecule:
            The molecule.
        :param ao_basis:
            The AO basis set.
        :param dft_func_label:
            The name of DFT functional.
        :param potfile_text:
            The content of potential file for polarizable embedding.
        :param eigvecs:
            The TDA eigenvectors (in columns).
        """

        if (not self.save_solutions) or (final_h5_fname is None):
            return

        for s in range(eigvecs.shape[1]):
            write_rsp_solution(final_h5_fname, 'S{:d}'.format(s + 1),
                               eigvecs[:, s])

        self.ostream.print_info('Response solution vectors written to file: ' +
                                final_h5_fname)
        self.ostream.print_blank()

    def _print_results(self, results):
        """
        Prints results to output stream.

        :param results:
            The dictionary containing response results.
        """

        self._print_transition_dipoles(
            'Electric Transition Dipole Moments (dipole length, a.u.)',
            results['electric_transition_dipoles'])

        self._print_transition_dipoles(
            'Electric Transition Dipole Moments (dipole velocity, a.u.)',
            results['velocity_transition_dipoles'])

        self._print_transition_dipoles(
            'Magnetic Transition Dipole Moments (a.u.)',
            results['magnetic_transition_dipoles'])

        self._print_absorption('One-Photon Absorption', results)
        self._print_ecd('Electronic Circular Dichroism', results)
        self._print_excitation_details('Character of excitations:', results)

    def __deepcopy__(self, memo):
        """
        Implements deepcopy.

        :param memo:
            The memo dictionary for deepcopy.

        :return:
            A deepcopy of self.
        """

        new_rsp_drv = TdaEigenSolver(self.comm, self.ostream)

        for key, val in vars(self).items():
            if isinstance(val, (MPI.Intracomm, OutputStream)):
                pass
            elif isinstance(val, XCFunctional):
                new_rsp_drv.key = XCFunctional(val)
            elif isinstance(val, MolecularGrid):
                new_rsp_drv.key = MolecularGrid(val)
            else:
                new_rsp_drv.key = deepcopy(val)

        return new_rsp_drv