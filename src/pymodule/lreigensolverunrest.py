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
from copy import deepcopy
import numpy as np
import time as tm
import h5py
import sys

from .oneeints import compute_electric_dipole_integrals
from .veloxchemlib import XCFunctional, MolecularGrid
from .veloxchemlib import (mpi_master, rotatory_strength_in_cgs, hartree_in_ev,
                           hartree_in_inverse_nm, fine_structure_constant,
                           extinction_coefficient_from_beta, avogadro_constant,
                           bohr_in_angstrom, hartree_in_wavenumber)
from .veloxchemlib import denmat
from .aodensitymatrix import AODensityMatrix
from .outputstream import OutputStream
from .profiler import Profiler
from .distributedarray import DistributedArray
from .linearsolver import LinearSolver
from .molecularorbitals import MolecularOrbitals
from .visualizationdriver import VisualizationDriver
from .cubicgrid import CubicGrid
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           dft_sanity_check, pe_sanity_check,
                           solvation_model_sanity_check)
from .errorhandler import assert_msg_critical
from .checkpoint import (check_rsp_hdf5, write_rsp_solution,
                         write_lr_rsp_results_to_hdf5,
                         write_detach_attach_to_hdf5)

try:
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
except ImportError:
    pass


class LinearResponseUnrestrictedEigenSolver(LinearSolver):
    """
    Implements linear response unrestricted eigensolver.

    # vlxtag: UHF, Absorption, ECD, TDHF
    # vlxtag: UKS, Absorption, ECD, TDDFT

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - nstates: Number of excited states.
        - core_excitation: The flag for computing core-excited states.
        - num_core_orbitals: The number of core-orbitals involved in
          the core-excited states.
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
        Initializes linear response eigensolver to default setup.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        super().__init__(comm, ostream)

        self.nstates = 3

        self.core_excitation = False
        self.num_core_orbitals = 0

        self.nto = False
        self.nto_pairs = None
        self.nto_cubes = False
        self.detach_attach = False
        self.detach_attach_cubes = False
        self.cube_origin = None
        self.cube_stepsize = None
        self.cube_points = [80, 80, 80]

        self.esa = False
        self.esa_from_state = None

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
            'esa': ('bool', 'compute excited state absorption'),
            'esa_from_state':
                ('int', 'the state to excite from (e.g. 1 for S1)'),
            'cube_origin': ('seq_fixed', 'origin of cubic grid points'),
            'cube_stepsize': ('seq_fixed', 'step size of cubic grid points'),
            'cube_points': ('seq_fixed_int', 'number of cubic grid points'),
        })

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in linear response eigensolver.

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
                'LinearResponseUnrestrictedEigenSolver: cube origin needs 3 numbers')

        if self.cube_stepsize is not None:
            assert_msg_critical(
                len(self.cube_stepsize) == 3,
                'LinearResponseUnrestrictedEigenSolver: cube stepsize needs 3 numbers')

        if self.cube_points is not None:
            assert_msg_critical(
                len(self.cube_points) == 3,
                'LinearResponseUnrestrictedEigenSolver: cube points needs 3 integers')

        # If the detachemnt and attachment cube files are requested
        # set the detach_attach flag to True to get the detachment and
        # attachment densities.
        if self.detach_attach_cubes:
            self.detach_attach = True

    def compute(self, molecule, basis, scf_tensors):
        """
        Performs linear response calculation for a molecule and a basis set.

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

        if self.norm_thresh is None:
            self.norm_thresh = self.conv_thresh * 1.0e-6
        if self.lindep_thresh is None:
            self.lindep_thresh = self.conv_thresh * 1.0e-2

        self._dist_bger = None
        self._dist_bung = None
        self._dist_e2bger = None
        self._dist_e2bung = None

        self._dist_fock_ger = None
        self._dist_fock_ung = None

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

        # check solvation
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
            self._print_header('Linear Response EigenSolver',
                               nstates=self.nstates)

        self.start_time = tm.time()

        if self.rank == mpi_master():
            orb_ene_a = scf_tensors['E_alpha']
            orb_ene_b = scf_tensors['E_beta']
        else:
            orb_ene_a = None
            orb_ene_b = None
        orb_ene_a = self.comm.bcast(orb_ene_a, root=mpi_master())
        orb_ene_b = self.comm.bcast(orb_ene_b, root=mpi_master())

        norb = orb_ene_a.shape[0]
        nocc_a = molecule.number_of_alpha_electrons()
        nocc_b = molecule.number_of_beta_electrons()

        if self.rank == mpi_master():
            assert_msg_critical(
                self.nstates <= (nocc_a * (norb - nocc_a) +
                                 nocc_b * (norb - nocc_b)),
                'LinearResponseUnrestrictedEigenSolver: too many excited states')

            if self.core_excitation:
                assert_msg_critical(
                    self.num_core_orbitals > 0,
                    'LinearResponseUnrestrictedEigenSolver: num_core_orbitals not set or invalid')
                assert_msg_critical(
                    (self.num_core_orbitals < nocc_a and
                     self.num_core_orbitals < nocc_b),
                    'LinearResponseUnrestrictedEigenSolver: num_core_orbitals too large')

        # ERI information
        eri_dict = self._init_eri(molecule, basis)

        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors)

        # PE information
        pe_dict = self._init_pe(molecule, basis)

        # CPCM_information
        self._init_cpcm(molecule)

        # For now, 'nonlinear' is not supported for unrestricted case.
        assert_msg_critical(
            not self.nonlinear,
            'LinearResponseUnrestrictedEigenSolver: ' +
            'not implemented for nonlinear')

        # TODO: enable PE
        assert_msg_critical(
            not self._pe,
            'LinearResponseUnrestrictedEigenSolver: ' +
            'not yet implemented for polarizable embedding')

        if self.nonlinear:
            # TODO: unrestricted
            pass
        else:
            rsp_vector_labels = [
                'LR_eigen_bger_half_size',
                'LR_eigen_bung_half_size',
                'LR_eigen_e2bger_half_size',
                'LR_eigen_e2bung_half_size',
            ]

        # check validity of restart file
        if self.restart:
            if self.rank == mpi_master():
                self.restart = check_rsp_hdf5(self.checkpoint_file,
                                              rsp_vector_labels, molecule,
                                              basis, dft_dict, pe_dict)
            self.restart = self.comm.bcast(self.restart, root=mpi_master())

        # read initial guess from restart file
        if self.restart:
            self._read_checkpoint(rsp_vector_labels)

            checkpoint_nstates = self._read_nstates_from_checkpoint()

            # print warning if nstates is not present in the restart file
            if checkpoint_nstates is None:
                self.ostream.print_warning(
                    'Could not find the nstates key in the checkpoint file.')
                self.ostream.print_blank()
                self.ostream.print_info(
                    'Assuming that nstates is not changed before and after ' +
                    'the restart.')
                self.ostream.print_blank()
                self.ostream.flush()

            # generate necessary initial guesses if more states are requested
            # in a restart calculation
            elif checkpoint_nstates < self.nstates:
                self.ostream.print_info(
                    'Generating initial guesses for ' +
                    f'{self.nstates - checkpoint_nstates} more states...')
                self.ostream.print_blank()

                igs = self._initial_excitations(self.nstates, orb_ene, nocc,
                                                norb, checkpoint_nstates)
                bger, bung = self._setup_trials(igs, None, self._dist_bger,
                                                self._dist_bung)

                profiler.set_timing_key('Preparation')

                self._e2n_half_size(bger, bung, molecule, basis, scf_tensors,
                                    eri_dict, dft_dict, pe_dict, profiler)

        # generate initial guess from scratch
        else:
            igs = self._initial_excitations(
                self.nstates, (orb_ene_a, orb_ene_b), (nocc_a, nocc_b), norb)
            bger, bung = self._setup_trials(igs, None)

            profiler.set_timing_key('Preparation')

            self._e2n_half_size(bger, bung, molecule, basis, scf_tensors,
                                eri_dict, dft_dict, pe_dict, profiler,
                                method_type='unrestricted')

        profiler.check_memory_usage('Initial guess')

        exc_energies = {}
        exc_focks = {}
        exc_solutions = {}
        exc_residuals = {}
        relative_residual_norm = {}

        iter_per_trial_in_hours = None

        # start iterations
        for iteration in range(self.max_iter):

            iter_start_time = tm.time()

            profiler.set_timing_key(f'Iteration {iteration + 1}')

            profiler.start_timer('ReducedSpace')

            self._cur_iter = iteration

            e2gg = self._dist_bger.matmul_AtB(self._dist_e2bger)
            e2uu = self._dist_bung.matmul_AtB(self._dist_e2bung)
            s2ug = self._dist_bung.matmul_AtB(self._dist_bger)

            if self.rank == mpi_master():

                # Equations:
                # E[2] X_g - w S[2] X_u = 0
                # E[2] X_u - w S[2] X_g = 0

                # Solutions:
                # (S_gu (E_uu)^-1 S_ug) X_g = 1/w^2 E_gg X_g
                # X_u = w (E_uu)^-1 S_ug X_g

                evals, evecs = np.linalg.eigh(e2uu)
                e2uu_inv = np.linalg.multi_dot(
                    [evecs, np.diag(1.0 / evals), evecs.T])
                ses = np.linalg.multi_dot([s2ug.T, e2uu_inv, s2ug])

                evals, evecs = np.linalg.eigh(e2gg)
                tmat = np.linalg.multi_dot(
                    [evecs, np.diag(1.0 / np.sqrt(evals)), evecs.T])
                ses_tilde = np.linalg.multi_dot([tmat.T, ses, tmat])

                evals, evecs = np.linalg.eigh(ses_tilde)
                p = list(reversed(evals.argsort()))
                evals = evals[p]
                evecs = evecs[:, p]

                wn = 1.0 / np.sqrt(evals[:self.nstates])
                c_ger = np.matmul(tmat, evecs[:, :self.nstates].copy())
                c_ung = wn * np.linalg.multi_dot([e2uu_inv, s2ug, c_ger])

                for k in range(self.nstates):
                    c_ger_k = c_ger[:, k].copy()
                    c_ung_k = c_ung[:, k].copy()
                    norm = np.sqrt(
                        np.linalg.multi_dot([c_ung_k.T, s2ug, c_ger_k]) +
                        np.linalg.multi_dot([c_ger_k.T, s2ug.T, c_ung_k]))
                    c_ger[:, k] /= norm
                    c_ung[:, k] /= norm
            else:
                wn = None
                c_ger, c_ung = None, None
            wn = self.comm.bcast(wn, root=mpi_master())
            c_ger = self.comm.bcast(c_ger, root=mpi_master())
            c_ung = self.comm.bcast(c_ung, root=mpi_master())

            for k in range(self.nstates):
                w = wn[k]

                x_ger = self._dist_bger.matmul_AB_no_gather(c_ger[:, k])
                x_ung = self._dist_bung.matmul_AB_no_gather(c_ung[:, k])

                e2x_ger = self._dist_e2bger.matmul_AB_no_gather(c_ger[:, k])
                e2x_ung = self._dist_e2bung.matmul_AB_no_gather(c_ung[:, k])

                if self.nonlinear:
                    # TODO: unrestricted
                    pass

                s2x_ger = x_ger.data
                s2x_ung = x_ung.data

                r_ger = e2x_ger.data - w * s2x_ung
                r_ung = e2x_ung.data - w * s2x_ger

                r_data = np.hstack((
                    r_ger.reshape(-1, 1),
                    r_ung.reshape(-1, 1),
                ))

                r = DistributedArray(r_data, self.comm, distribute=False)

                x_data = np.hstack((
                    x_ger.data.reshape(-1, 1),
                    x_ung.data.reshape(-1, 1),
                ))

                x = DistributedArray(x_data, self.comm, distribute=False)

                r_norms_2 = r.squared_norm(axis=0)
                x_norms_2 = x.squared_norm(axis=0)

                rn = np.sqrt(np.sum(r_norms_2))
                xn = np.sqrt(np.sum(x_norms_2))

                if xn != 0:
                    relative_residual_norm[k] = 2.0 * rn / xn
                else:
                    relative_residual_norm[k] = 2.0 * rn

                if relative_residual_norm[k] < self.conv_thresh:
                    exc_energies[k] = w
                    exc_solutions[k] = x
                else:
                    exc_residuals[k] = r

            # write to output
            if self.rank == mpi_master():
                self.ostream.print_info(
                    '{:d} gerade trial vectors in reduced space'.format(
                        self._dist_bger.shape(1)))
                self.ostream.print_info(
                    '{:d} ungerade trial vectors in reduced space'.format(
                        self._dist_bung.shape(1)))
                self.ostream.print_blank()

                profiler.print_memory_subspace(
                    {
                        'dist_bger': self._dist_bger,
                        'dist_bung': self._dist_bung,
                        'dist_e2bger': self._dist_e2bger,
                        'dist_e2bung': self._dist_e2bung,
                        'exc_solutions': exc_solutions,
                        'exc_residuals': exc_residuals,
                    }, self.ostream)

                profiler.check_memory_usage(
                    'Iteration {:d} subspace'.format(iteration + 1))

                profiler.print_memory_tracing(self.ostream)

                self._print_iteration(relative_residual_norm, wn)

            profiler.stop_timer('ReducedSpace')

            # check convergence
            self._check_convergence(relative_residual_norm)

            if self.is_converged:
                break

            profiler.start_timer('Orthonorm.')

            # update trial vectors
            precond = {
                k: self._get_precond((orb_ene_a, orb_ene_b), (nocc_a, nocc_b), norb, w)
                for k, w in enumerate(list(wn))
            }

            new_trials_ger, new_trials_ung = self._setup_trials(
                exc_residuals, precond, self._dist_bger, self._dist_bung)

            exc_residuals.clear()

            profiler.stop_timer('Orthonorm.')

            if self.rank == mpi_master():
                n_new_trials = new_trials_ger.shape(1) + new_trials_ung.shape(1)
            else:
                n_new_trials = None
            n_new_trials = self.comm.bcast(n_new_trials, root=mpi_master())

            if iter_per_trial_in_hours is not None:
                next_iter_in_hours = iter_per_trial_in_hours * n_new_trials
                if self._need_graceful_exit(next_iter_in_hours):
                    self._graceful_exit(molecule, basis, dft_dict, pe_dict,
                                        rsp_vector_labels)

            if self.force_checkpoint:
                self._write_checkpoint(molecule, basis, dft_dict, pe_dict,
                                       rsp_vector_labels)
                self._add_nstates_to_checkpoint()

            self._e2n_half_size(new_trials_ger, new_trials_ung,
                                molecule, basis, scf_tensors, eri_dict,
                                dft_dict, pe_dict, profiler,
                                method_type='unrestricted')

            iter_in_hours = (tm.time() - iter_start_time) / 3600
            iter_per_trial_in_hours = iter_in_hours / n_new_trials

            profiler.check_memory_usage(
                'Iteration {:d} sigma build'.format(iteration + 1))

        self._write_checkpoint(molecule, basis, dft_dict, pe_dict,
                               rsp_vector_labels)
        self._add_nstates_to_checkpoint()

        # converged?
        if self.rank == mpi_master():
            self._print_convergence('Linear response')

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of LR eigensolver')
        profiler.print_memory_usage(self.ostream)

        self._dist_bger = None
        self._dist_bung = None
        self._dist_e2bger = None
        self._dist_e2bung = None

        self._dist_fock_ger = None
        self._dist_fock_ung = None

        # calculate properties
        if self.is_converged:
            edip_grad_a = self.get_prop_grad('electric dipole', 'xyz', molecule,
                                           basis, scf_tensors, spin='alpha')
            lmom_grad_a = self.get_prop_grad('linear momentum', 'xyz', molecule,
                                           basis, scf_tensors, spin='alpha')
            mdip_grad_a = self.get_prop_grad('magnetic dipole', 'xyz', molecule,
                                           basis, scf_tensors, spin='alpha')

            edip_grad_b = self.get_prop_grad('electric dipole', 'xyz', molecule,
                                           basis, scf_tensors, spin='beta')
            lmom_grad_b = self.get_prop_grad('linear momentum', 'xyz', molecule,
                                           basis, scf_tensors, spin='beta')
            mdip_grad_b = self.get_prop_grad('magnetic dipole', 'xyz', molecule,
                                           basis, scf_tensors, spin='beta')

            eigvals = np.array([exc_energies[s] for s in range(self.nstates)])

            elec_trans_dipoles = np.zeros((self.nstates, 3))
            velo_trans_dipoles = np.zeros((self.nstates, 3))
            magn_trans_dipoles = np.zeros((self.nstates, 3))

            if self.rank == mpi_master():
                # final h5 file for response solutions
                if self.filename is not None:
                    final_h5_fname = f'{self.filename}.h5'
                else:
                    final_h5_fname = None

            nto_lambdas_a = []
            nto_lambdas_b = []
            nto_cube_files_a = []
            nto_cube_files_b = []
            dens_cube_files = []

            excitation_details = []

            for s in range(self.nstates):
                eigvec_full = self.get_full_solution_vector(exc_solutions[s])

                if self.rank == mpi_master():
                    if self.core_excitation:
                        n_ov_a = self.num_core_orbitals * (norb - nocc_a)
                        n_ov_b = self.num_core_orbitals * (norb - nocc_b)
                    else:
                        n_ov_a = nocc_a * (norb - nocc_a)
                        n_ov_b = nocc_b * (norb - nocc_b)

                    eigvec_a = np.hstack((
                        eigvec_full[:n_ov_a],
                        eigvec_full[n_ov_a + n_ov_b:n_ov_a + n_ov_b + n_ov_a],
                        ))
                    eigvec_b = np.hstack((
                        eigvec_full[n_ov_a:n_ov_a + n_ov_b],
                        eigvec_full[n_ov_a + n_ov_b + n_ov_a:],
                        ))

                    if self.core_excitation:
                        mo_occ_a = scf_tensors['C_alpha'][:, :self.num_core_orbitals]
                        mo_vir_a = scf_tensors['C_alpha'][:, nocc_a:]
                        z_mat_a = eigvec_a[:eigvec_a.size // 2].reshape(self.num_core_orbitals, -1)
                        y_mat_a = eigvec_a[eigvec_a.size // 2:].reshape(self.num_core_orbitals, -1)

                        mo_occ_b = scf_tensors['C_beta'][:, :self.num_core_orbitals]
                        mo_vir_b = scf_tensors['C_beta'][:, nocc_b:]
                        z_mat_b = eigvec_b[:eigvec_b.size // 2].reshape(self.num_core_orbitals, -1)
                        y_mat_b = eigvec_b[eigvec_b.size // 2:].reshape(self.num_core_orbitals, -1)
                    else:
                        mo_occ_a = scf_tensors['C_alpha'][:, :nocc_a]
                        mo_vir_a = scf_tensors['C_alpha'][:, nocc_a:]
                        z_mat_a = eigvec_a[:eigvec_a.size // 2].reshape(nocc_a, -1)
                        y_mat_a = eigvec_a[eigvec_a.size // 2:].reshape(nocc_a, -1)

                        mo_occ_b = scf_tensors['C_beta'][:, :nocc_b]
                        mo_vir_b = scf_tensors['C_beta'][:, nocc_b:]
                        z_mat_b = eigvec_b[:eigvec_b.size // 2].reshape(nocc_b, -1)
                        y_mat_b = eigvec_b[eigvec_b.size // 2:].reshape(nocc_b, -1)

                # TODO: enable nto and detach_attach
                #assert_msg_critical(
                #    not (self.nto or self.detach_attach),
                #    'LinearResponseUnrestrictedEigenSolver: ' +
                #    'not yet implemented for nto or detach_attach')

                if self.nto or self.detach_attach:
                    vis_drv = VisualizationDriver(self.comm)
                    if self.cube_origin is None or self.cube_stepsize is None:
                        cubic_grid = vis_drv.gen_cubic_grid(
                            molecule, self.cube_points)
                    else:
                        cubic_grid = CubicGrid(self.cube_origin,
                                               self.cube_stepsize,
                                               self.cube_points)

                if self.nto:
                    self.ostream.print_info(
                        'Running NTO analysis for S{:d}...'.format(s + 1))
                    self.ostream.flush()

                    if self.rank == mpi_master():
                        nto_mo = self.get_nto_unrestricted(
                            (z_mat_a - y_mat_a, z_mat_b - y_mat_b),
                            (mo_occ_a, mo_occ_b), (mo_vir_a, mo_vir_b))

                        nto_lam_a = nto_mo.occa_to_numpy()
                        lam_start_a = mo_occ_a.shape[1]
                        lam_end_a = lam_start_a + min(mo_occ_a.shape[1],
                                                      mo_vir_a.shape[1])
                        nto_lambdas_a.append(nto_lam_a[lam_start_a:lam_end_a])

                        nto_lam_b = nto_mo.occb_to_numpy()
                        lam_start_b = mo_occ_b.shape[1]
                        lam_end_b = lam_start_b + min(mo_occ_b.shape[1],
                                                      mo_vir_b.shape[1])
                        nto_lambdas_b.append(nto_lam_b[lam_start_b:lam_end_b])

                        # Add the NTO to the final checkpoint file.
                        nto_label = f'NTO_S{s + 1}'
                        if final_h5_fname is not None:
                            nto_mo.write_hdf5(final_h5_fname,
                                              label=f'rsp/{nto_label}')
                    else:
                        nto_mo = MolecularOrbitals()
                    nto_mo = nto_mo.broadcast(self.comm, root=mpi_master())

                    if self.nto_cubes:
                        lam_diag_a, nto_cube_fnames_a = self.write_nto_cubes(
                            cubic_grid, molecule, basis, s, nto_mo,
                            self.nto_pairs, nto_spin='alpha')
                        lam_diag_b, nto_cube_fnames_b = self.write_nto_cubes(
                            cubic_grid, molecule, basis, s, nto_mo,
                            self.nto_pairs, nto_spin='beta')

                        if self.rank == mpi_master():
                            nto_cube_files_a.append(nto_cube_fnames_a)
                            nto_cube_files_b.append(nto_cube_fnames_b)

                if self.detach_attach:
                    self.ostream.print_info(
                        'Running detachment/attachment analysis for S{:d}...'.
                        format(s + 1))
                    self.ostream.flush()

                    if self.rank == mpi_master():
                        dens_D, dens_A = self.get_detach_attach_densities(
                            z_mat, y_mat, mo_occ, mo_vir)

                        # Add the detachment and attachment density matrices
                        # to the checkpoint file
                        state_label = f'S{s + 1}'
                        if final_h5_fname is not None:
                            write_detach_attach_to_hdf5(final_h5_fname,
                                                        state_label, dens_D,
                                                        dens_A)

                    if self.detach_attach_cubes:
                        if self.rank == mpi_master():
                            dens_DA = AODensityMatrix([dens_D, dens_A],
                                                      denmat.rest)
                        else:
                            dens_DA = AODensityMatrix()
                        dens_DA = dens_DA.broadcast(self.comm,
                                                    root=mpi_master())

                        dens_cube_fnames = self.write_detach_attach_cubes(
                            cubic_grid, molecule, basis, s, dens_DA)

                        if self.rank == mpi_master():
                            dens_cube_files.append(dens_cube_fnames)

                if self.esa:
                    if self.esa_from_state is None:
                        source_states = list(range(self.nstates))
                    else:
                        source_states = [self.esa_from_state - 1]

                    esa_pairs = [(s_1, s_2)
                                 for s_1 in source_states
                                 for s_2 in range(s_1 + 1, self.nstates)]

                    if self.rank == mpi_master():
                        esa_results = []
                        dipole_integrals = compute_electric_dipole_integrals(
                            molecule, basis, [0.0, 0.0, 0.0])

                    for s_1, s_2 in esa_pairs:
                        eigvec_1 = self.get_full_solution_vector(
                            exc_solutions[s_1])
                        eigvec_2 = self.get_full_solution_vector(
                            exc_solutions[s_2])

                        if self.rank == mpi_master():
                            half_size = eigvec_1.shape[0] // 2

                            z_mat_1 = eigvec_1[:half_size].reshape(nocc, -1)
                            y_mat_1 = eigvec_1[half_size:].reshape(nocc, -1)

                            z_mat_2 = eigvec_2[:half_size].reshape(nocc, -1)
                            y_mat_2 = eigvec_2[half_size:].reshape(nocc, -1)

                            esa_trans_dens = (
                                np.linalg.multi_dot(
                                    [mo_vir, z_mat_1.T, z_mat_2, mo_vir.T]) -
                                np.linalg.multi_dot(
                                    [mo_occ, z_mat_1, z_mat_2.T, mo_occ.T]))

                            esa_trans_dens += (
                                np.linalg.multi_dot(
                                    [mo_occ, y_mat_1, y_mat_2.T, mo_occ.T]) -
                                np.linalg.multi_dot(
                                    [mo_vir, y_mat_1.T, y_mat_2, mo_vir.T]))

                            esa_trans_dipole = np.array([
                                np.sum(esa_trans_dens * dipole_integrals[i])
                                for i in range(3)
                            ])

                            esa_exc_ene = exc_energies[s_2] - exc_energies[s_1]
                            esa_osc_str = (2.0 / 3.0) * esa_exc_ene * np.sum(
                                esa_trans_dipole**2)

                            esa_results.append({
                                'from_state': f'S{s_1 + 1}',
                                'to_state': f'S{s_2 + 1}',
                                'excitation_energy': esa_exc_ene,
                                'oscillator_strength': esa_osc_str,
                                'transition_dipole': esa_trans_dipole,
                            })

                if self.rank == mpi_master():
                    for ind, comp in enumerate('xyz'):
                        elec_trans_dipoles[s, ind] = np.vdot(
                            edip_grad_a[ind], eigvec_a)
                        velo_trans_dipoles[s, ind] = np.vdot(
                            lmom_grad_a[ind], eigvec_a) / (-eigvals[s])
                        magn_trans_dipoles[s, ind] = np.vdot(
                            mdip_grad_a[ind], eigvec_a)

                        elec_trans_dipoles[s, ind] += np.vdot(
                            edip_grad_b[ind], eigvec_b)
                        velo_trans_dipoles[s, ind] += np.vdot(
                            lmom_grad_b[ind], eigvec_b) / (-eigvals[s])
                        magn_trans_dipoles[s, ind] += np.vdot(
                            mdip_grad_b[ind], eigvec_b)

                        # TODO: move factor elsewhere
                        elec_trans_dipoles[s, ind] /= 2.0
                        velo_trans_dipoles[s, ind] /= 2.0
                        magn_trans_dipoles[s, ind] /= 2.0

                    # write to h5 file for response solutions
                    if (self.save_solutions and final_h5_fname is not None):
                        write_rsp_solution(final_h5_fname,
                                           'S{:d}(a)'.format(s + 1), eigvec_a)
                        write_rsp_solution(final_h5_fname,
                                           'S{:d}(a)'.format(s + 1), eigvec_b)

                    # save excitation details
                    excitation_details.append(
                        self.get_excitation_details_unrestricted(
                            (eigvec_a, eigvec_b),
                            (mo_occ_a.shape[1], mo_occ_b.shape[1]),
                            (mo_vir_a.shape[1], mo_vir_b.shape[1])))

            if self.nto or self.detach_attach:
                self.ostream.print_blank()
                self.ostream.flush()

            if not self.nonlinear:

                if self.rank == mpi_master():
                    osc = (2.0 / 3.0) * np.sum(elec_trans_dipoles**2,
                                               axis=1) * eigvals
                    rot_vel = np.sum(velo_trans_dipoles * magn_trans_dipoles,
                                     axis=1) * rotatory_strength_in_cgs()

                    ret_dict = {
                        'eigenvalues': eigvals,
                        'eigenvectors_distributed': exc_solutions,
                        'electric_transition_dipoles': elec_trans_dipoles,
                        'velocity_transition_dipoles': velo_trans_dipoles,
                        'magnetic_transition_dipoles': magn_trans_dipoles,
                        'oscillator_strengths': osc,
                        'rotatory_strengths': rot_vel,
                        'excitation_details': excitation_details,
                        'number_of_states': self.nstates,
                    }

                    if self.nto:
                        ret_dict['nto_lambdas_a'] = nto_lambdas_a
                        ret_dict['nto_lambdas_b'] = nto_lambdas_b
                        if self.nto_cubes:
                            ret_dict['nto_cubes_a'] = nto_cube_files_a
                            ret_dict['nto_cubes_b'] = nto_cube_files_b

                    if self.detach_attach_cubes:
                        ret_dict['density_cubes'] = dens_cube_files

                    if self.esa:
                        ret_dict['esa_results'] = esa_results

                    if (self.save_solutions and final_h5_fname is not None):
                        self.ostream.print_info(
                            'Response solution vectors written to file: ' +
                            final_h5_fname)
                        self.ostream.print_blank()

                        # Write the response results to the final checkpoint file
                        write_lr_rsp_results_to_hdf5(final_h5_fname, ret_dict)

                    self._print_results(ret_dict)

                    return ret_dict
                else:
                    return {
                        'eigenvalues': eigvals,
                        'eigenvectors_distributed': exc_solutions,
                    }

            else:

                # TODO: unrestricted for nonlinear
                pass

        return None

    def _add_nstates_to_checkpoint(self):
        """
        Add nstates to checkpoint file.
        """

        if self.checkpoint_file is None:
            return

        if self.rank == mpi_master():
            hf = h5py.File(self.checkpoint_file, 'a')
            key = 'nstates'
            if key in hf:
                del hf[key]
            hf.create_dataset(key, data=np.array([self.nstates]))
            hf.close()

        self.comm.barrier()

    def _read_nstates_from_checkpoint(self):
        """
        Read nstates from checkpoint file.

        :return:
            The number of states.
        """

        if self.checkpoint_file is None:
            return None

        nstates = None

        if self.rank == mpi_master():
            hf = h5py.File(self.checkpoint_file, 'r')
            key = 'nstates'
            if key in hf:
                nstates = np.array(hf.get(key))[0]
            hf.close()

        nstates = self.comm.bcast(nstates, root=mpi_master())

        return nstates

    @staticmethod
    def get_full_solution_vector(solution):
        """
        Gets a full solution vector from the distributed solution.

        :param solution:
            The distributed solution as a tuple.

        :return:
            The full solution vector.
        """

        x_ger = solution.get_full_vector(0)
        x_ung = solution.get_full_vector(1)

        if solution.rank == mpi_master():
            x_ger_full = np.hstack((x_ger, x_ger))
            x_ung_full = np.hstack((x_ung, -x_ung))
            return x_ger_full + x_ung_full
        else:
            return None

    def _print_iteration(self, relative_residual_norm, ws):
        """
        Prints information of the iteration.

        :param relative_residual_norm:
            Relative residual norms.
        :param ws:
            Excitation energies.
        """

        width = 92
        output_header = '*** Iteration:   {} '.format(self._cur_iter + 1)
        output_header += '* Residuals (Max,Min): '
        output_header += '{:.2e} and {:.2e}'.format(
            max(relative_residual_norm.values()),
            min(relative_residual_norm.values()))
        self.ostream.print_header(output_header.ljust(width))
        self.ostream.print_blank()
        for k, w in enumerate(ws):
            state_label = 'Excitation {}'.format(k + 1)
            rel_res = relative_residual_norm[k]
            output_iter = '{:<15s}: {:15.8f} '.format(state_label, w)
            output_iter += 'Residual Norm: {:.8f}'.format(rel_res)
            if relative_residual_norm[k] < self.conv_thresh:
                output_iter += '   converged'
            self.ostream.print_header(output_iter.ljust(width))
        self.ostream.print_blank()
        self.ostream.flush()

    def _initial_excitations(self, nstates, orb_ene, nocc, norb, n_excl_states=0):
        """
        Gets initial guess for excitations.

        :param nstates:
            Number of excited states.
        :param orb_ene:
            Orbital energies (alpha and beta).
        :param nocc:
            Number of occupied orbitals (alpha and beta).
        :param norb:
            Number of orbitals.
        :param n_excl_states:
            Number of states to exclude. Useful for generating initial guess
            for a restarting calculation that requests more states.

        :return:
            A list of initial excitations (excitation energy and distributed
            vector).
        """

        nocc_a, nocc_b = nocc
        ea, eb = orb_ene

        if self.core_excitation:
            excitations_a = [(i, a)
                             for i in range(self.num_core_orbitals)
                             for a in range(nocc_a, norb)]
            excitations_b = [(j, b)
                             for j in range(self.num_core_orbitals)
                             for b in range(nocc_b, norb)]
        else:
            excitations_a = [
                (i, a) for i in range(nocc_a) for a in range(nocc_a, norb)
            ]
            excitations_b = [
                (j, b) for j in range(nocc_b) for b in range(nocc_b, norb)
            ]

        excitation_energies_a = [ea[a] - ea[i] for i, a in excitations_a]
        excitation_energies_b = [eb[b] - eb[j] for j, b in excitations_b]

        w_a = {ia: w for ia, w in zip(excitations_a, excitation_energies_a)}
        w_b = {jb: w for jb, w in zip(excitations_b, excitation_energies_b)}
        n_exc_a = len(excitations_a)
        n_exc_b = len(excitations_b)

        final = {}

        for k, (i, a) in enumerate(sorted(w_a, key=w_a.get)[n_excl_states:nstates]):
            if self.rank == mpi_master():
                ia = excitations_a.index((i, a))

                Xn = np.zeros(2 * n_exc_a)
                Xn[ia] = 1.0

                Xn_T = np.zeros(2 * n_exc_a)
                Xn_T[:n_exc_a] = Xn[n_exc_a:]
                Xn_T[n_exc_a:] = Xn[:n_exc_a]

                Xn_ger = 0.5 * (Xn + Xn_T)[:n_exc_a]
                Xn_ung = 0.5 * (Xn - Xn_T)[:n_exc_a]

                X_alpha = np.hstack((
                    Xn_ger.reshape(-1, 1),
                    Xn_ung.reshape(-1, 1),
                ))
                X_beta = np.hstack((
                    np.zeros(n_exc_b).reshape(-1, 1),
                    np.zeros(n_exc_b).reshape(-1, 1),
                ))

                # put alpha and beta together
                X = np.vstack((X_alpha, X_beta))
            else:
                X = None

            final[k] = DistributedArray(X, self.comm)

        k_offset = len(final)

        for k, (j, b) in enumerate(sorted(w_b, key=w_b.get)[n_excl_states:nstates]):
            if self.rank == mpi_master():
                jb = excitations_b.index((j, b))

                Xn = np.zeros(2 * n_exc_b)
                Xn[jb] = 1.0

                Xn_T = np.zeros(2 * n_exc_b)
                Xn_T[:n_exc_b] = Xn[n_exc_b:]
                Xn_T[n_exc_b:] = Xn[:n_exc_b]

                Xn_ger = 0.5 * (Xn + Xn_T)[:n_exc_b]
                Xn_ung = 0.5 * (Xn - Xn_T)[:n_exc_b]

                X_alpha = np.hstack((
                    np.zeros(n_exc_a).reshape(-1, 1),
                    np.zeros(n_exc_a).reshape(-1, 1),
                ))
                X_beta = np.hstack((
                    Xn_ger.reshape(-1, 1),
                    Xn_ung.reshape(-1, 1),
                ))

                # put alpha and beta together
                X = np.vstack((X_alpha, X_beta))
            else:
                X = None

            final[k + k_offset] = DistributedArray(X, self.comm)

        return final

    def _precond_trials(self, vectors, precond):
        """
        Applies preconditioner to distributed trial vectors.

        :param vectors:
            The set of vectors.
        :param precond:
            The preconditioner.

        :return:
            The preconditioned gerade and ungerade trial vectors.
        """

        trials_ger = []
        trials_ung = []

        for k, X in vectors.items():
            if precond is not None:
                v = self._preconditioning(precond[k], X)
            else:
                v = X
            norms_2 = 2.0 * v.squared_norm(axis=0)
            vn = np.sqrt(np.sum(norms_2))

            if vn > self.norm_thresh:
                norms = np.sqrt(norms_2)
                # gerade
                if norms[0] > self.norm_thresh:
                    trials_ger.append(v.data[:, 0])
                # ungerade
                if norms[1] > self.norm_thresh:
                    trials_ung.append(v.data[:, 1])

        new_ger = np.array(trials_ger).T
        new_ung = np.array(trials_ung).T

        dist_new_ger = DistributedArray(new_ger, self.comm, distribute=False)
        dist_new_ung = DistributedArray(new_ung, self.comm, distribute=False)

        return dist_new_ger, dist_new_ung

    def _get_precond(self, orb_ene, nocc, norb, w):
        """
        Constructs the preconditioners.

        :param orb_ene:
            The orbital energies.
        :param nocc:
            The number of doubly occupied orbitals.
        :param norb:
            The number of orbitals.
        :param w:
            The frequency.

        :return:
            The distributed preconditioners.
        """

        orb_ene_a, orb_ene_b = orb_ene
        nocc_a, nocc_b = nocc

        # spawning needed components

        if self.core_excitation:
            ediag_a, sdiag_a = self.construct_ediag_sdiag_half(
                orb_ene_a, nocc_a, norb, self.num_core_orbitals)
            ediag_b, sdiag_b = self.construct_ediag_sdiag_half(
                orb_ene_b, nocc_b, norb, self.num_core_orbitals)
        else:
            ediag_a, sdiag_a = self.construct_ediag_sdiag_half(orb_ene_a, nocc_a, norb)
            ediag_b, sdiag_b = self.construct_ediag_sdiag_half(orb_ene_b, nocc_b, norb)

        ediag_a_sq = ediag_a**2
        sdiag_a_sq = sdiag_a**2

        ediag_b_sq = ediag_b**2
        sdiag_b_sq = sdiag_b**2

        w_sq = w**2

        # constructing matrix block diagonals

        pa_alpha_diag = ediag_a / (ediag_a_sq - w_sq * sdiag_a_sq)
        pb_alpha_diag = (w * sdiag_a) / (ediag_a_sq - w_sq * sdiag_a_sq)

        pa_beta_diag = ediag_b / (ediag_b_sq - w_sq * sdiag_b_sq)
        pb_beta_diag = (w * sdiag_b) / (ediag_b_sq - w_sq * sdiag_b_sq)

        p_mat_alpha = np.hstack((
            pa_alpha_diag.reshape(-1, 1),
            pb_alpha_diag.reshape(-1, 1),
        ))
        p_mat_beta = np.hstack((
            pa_beta_diag.reshape(-1, 1),
            pb_beta_diag.reshape(-1, 1),
        ))

        # put alpha and beta together
        p_mat = np.vstack((p_mat_alpha, p_mat_beta))

        return DistributedArray(p_mat, self.comm)

    def _preconditioning(self, precond, v_in):
        """
        Applies preconditioner to a tuple of distributed trial vectors.

        :param precond:
            The preconditioner.
        :param v_in:
            The input trial vectors.

        :return:
            A tuple of distributed trial vectors after preconditioning.
        """

        pa = precond.data[:, 0]
        pb = precond.data[:, 1]

        v_in_rg = v_in.data[:, 0]
        v_in_ru = v_in.data[:, 1]

        v_out_rg = pa * v_in_rg + pb * v_in_ru
        v_out_ru = pb * v_in_rg + pa * v_in_ru

        v_mat = np.hstack((
            v_out_rg.reshape(-1, 1),
            v_out_ru.reshape(-1, 1),
        ))

        return DistributedArray(v_mat, self.comm, distribute=False)

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

        if self.esa:
            self._print_excited_state_absorption('Excited state absorption:',
                                                 results)

    def __deepcopy__(self, memo):
        """
        Implements deepcopy.

        :param memo:
            The memo dictionary for deepcopy.

        :return:
            A deepcopy of self.
        """

        new_rsp_drv = LinearResponseUnrestrictedEigenSolver(self.comm, self.ostream)

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

    def plot_uv_vis(self,
                    rpa_results,
                    broadening_type="lorentzian",
                    broadening_value=(1000.0 / hartree_in_wavenumber() *
                                      hartree_in_ev()),
                    ax=None):
        """
        Plot the UV spectrum from the response calculation.

        :param rpa_results:
            The dictionary containing RPA results.
        :param broadening_type:
            The type of broadening to use. Either 'lorentzian' or 'gaussian'.
        :param broadening_value:
            The broadening value in eV.
        :param ax:
            The matplotlib axis to plot on.
        """

        assert_msg_critical('matplotlib' in sys.modules,
                            'matplotlib is required.')

        ev_x_nm = hartree_in_ev() / hartree_in_inverse_nm()
        au2ev = hartree_in_ev()
        ev2au = 1.0 / au2ev

        # initialize the plot
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 5))

        ax.set_xlabel('Wavelength [nm]')
        ax.set_ylabel(r'$\epsilon$ [L mol$^{-1}$ cm$^{-1}$]')

        ax.set_title("Absorption Spectrum")

        x = (rpa_results['eigenvalues'])
        y = rpa_results['oscillator_strengths']
        xmin = min(x) - 0.03
        xmax = max(x) + 0.03
        xstep = 0.0001

        ax2 = ax.twinx()

        for i in np.arange(len(rpa_results['eigenvalues'])):
            ax2.plot(
                [
                    ev_x_nm / (rpa_results['eigenvalues'][i] * au2ev),
                    ev_x_nm / (rpa_results['eigenvalues'][i] * au2ev),
                ],
                [0.0, rpa_results['oscillator_strengths'][i]],
                alpha=0.7,
                linewidth=2,
                color="darkcyan",
            )

        c = 1.0 / fine_structure_constant()
        NA = avogadro_constant()
        a_0 = bohr_in_angstrom() * 1.0e-10

        if broadening_type.lower() == "lorentzian":
            xi, yi = self.lorentzian_uv_vis(x, y, xmin, xmax, xstep,
                                            broadening_value * ev2au)

        elif broadening_type.lower() == "gaussian":
            xi, yi = self.gaussian_uv_vis(x, y, xmin, xmax, xstep,
                                          broadening_value * ev2au)

        sigma = (2 * np.pi * np.pi * xi * yi) / c
        sigma_m2 = sigma * a_0**2
        sigma_cm2 = sigma_m2 * 10**4
        epsilon = sigma_cm2 * NA / (np.log(10) * 10**3)
        ax.plot(ev_x_nm / (xi * au2ev),
                epsilon,
                color="black",
                alpha=0.9,
                linewidth=2.5)

        legend_bars = mlines.Line2D([], [],
                                    color='darkcyan',
                                    alpha=0.7,
                                    linewidth=2,
                                    label='Oscillator strength')
        label_spectrum = f'{broadening_type.capitalize()} '
        label_spectrum += f'broadening ({broadening_value:.3f} eV)'
        legend_spectrum = mlines.Line2D([], [],
                                        color='black',
                                        linestyle='-',
                                        linewidth=2.5,
                                        label=label_spectrum)
        ax2.legend(handles=[legend_bars, legend_spectrum],
                   frameon=False,
                   borderaxespad=0.,
                   loc='center left',
                   bbox_to_anchor=(1.15, 0.5))
        ax2.set_ylim(0, max(abs(rpa_results['oscillator_strengths'])) * 1.1)
        ax.set_ylim(0, max(epsilon) * 1.1)
        ax.set_ylim(bottom=0)
        ax2.set_ylim(bottom=0)
        ax2.set_ylabel("Oscillator strength")
        ax.set_xlim(ev_x_nm / (xmax * au2ev), ev_x_nm / (xmin * au2ev))

    def plot_ecd(self,
                 rpa_results,
                 broadening_type="lorentzian",
                 broadening_value=(1000.0 / hartree_in_wavenumber() *
                                   hartree_in_ev()),
                 ax=None):
        """
        Plot the ECD spectrum from the response calculation.

        :param rpa_results:
            The dictionary containing RPA results.
        :param broadening_type:
            The type of broadening to use. Either 'lorentzian' or 'gaussian'.
        :param broadening_value:
            The broadening value in eV.
        :param ax:
            The matplotlib axis to plot on.
        """

        assert_msg_critical('matplotlib' in sys.modules,
                            'matplotlib is required.')

        ev_x_nm = hartree_in_ev() / hartree_in_inverse_nm()
        au2ev = hartree_in_ev()

        # initialize the plot
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 5))

        ax.set_xlabel("Wavelength [nm]")
        ax.set_title("ECD Spectrum")
        ax.set_ylabel(r'$\Delta \epsilon$ [L mol$^{-1}$ cm$^{-1}$]')

        ax2 = ax.twinx()
        ax2.set_ylabel('Rotatory strength [10$^{-40}$ cgs]')

        for i in np.arange(len(rpa_results["eigenvalues"])):
            ax2.plot(
                [
                    ev_x_nm / (rpa_results["eigenvalues"][i] * au2ev),
                    ev_x_nm / (rpa_results["eigenvalues"][i] * au2ev),
                ],
                [0.0, rpa_results["rotatory_strengths"][i]],
                alpha=0.7,
                linewidth=2,
                color="darkcyan",
            )
        ax2.set_ylim(-max(abs(rpa_results["rotatory_strengths"])) * 1.1,
                     max(abs(rpa_results["rotatory_strengths"])) * 1.1)

        ax.axhline(y=0,
                   marker=',',
                   color='k',
                   linestyle='-.',
                   markersize=0,
                   linewidth=0.2)

        x = (rpa_results["eigenvalues"]) * au2ev
        y = rpa_results["rotatory_strengths"]
        xmin = min(x) - 0.8
        xmax = max(x) + 0.8
        xstep = 0.003

        if broadening_type.lower() == "lorentzian":
            xi, yi = self.lorentzian_ecd(x, y, xmin, xmax, xstep,
                                         broadening_value)

        elif broadening_type.lower() == "gaussian":
            xi, yi = self.gaussian_ecd(x, y, xmin, xmax, xstep,
                                       broadening_value)

        # denorm_factor is roughly 22.96 * PI
        denorm_factor = (rotatory_strength_in_cgs() /
                         (extinction_coefficient_from_beta() / 3.0))
        yi = (yi * xi) / denorm_factor

        ax.set_ylim(-max(abs(yi)) * 1.1, max(abs(yi)) * 1.1)

        ax.plot(ev_x_nm / xi, yi, color="black", alpha=0.9, linewidth=2.5)
        ax.set_xlim(ev_x_nm / xmax, ev_x_nm / xmin)

        # include a legend for the bar and for the broadened spectrum
        legend_bars = mlines.Line2D([], [],
                                    color='darkcyan',
                                    alpha=0.7,
                                    linewidth=2,
                                    label='Rotatory strength')
        label_spectrum = f'{broadening_type.capitalize()} '
        label_spectrum += f'broadening ({broadening_value:.3f} eV)'
        legend_spectrum = mlines.Line2D([], [],
                                        color='black',
                                        linestyle='-',
                                        linewidth=2.5,
                                        label=label_spectrum)
        ax.legend(handles=[legend_bars, legend_spectrum],
                  frameon=False,
                  borderaxespad=0.,
                  loc='center left',
                  bbox_to_anchor=(1.15, 0.5))

    @staticmethod
    def lorentzian_uv_vis(x, y, xmin, xmax, xstep, gamma):
        xi = np.arange(xmin, xmax, xstep)
        yi = np.zeros(len(xi))
        for i in range(len(xi)):
            for k in range(len(x)):
                yi[i] = yi[i] + y[k] / x[k] * gamma / (
                    (xi[i] - x[k])**2 + gamma**2)
        yi = yi / np.pi
        return xi, yi

    @staticmethod
    def gaussian_uv_vis(x, y, xmin, xmax, xstep, sigma):
        xi = np.arange(xmin, xmax, xstep)
        yi = np.zeros(len(xi))
        for i in range(len(xi)):
            for k in range(len(x)):
                yi[i] = yi[i] + y[k] / x[k] * np.exp(-((xi[i] - x[k])**2) /
                                                     (2 * sigma**2))
        yi = yi / (sigma * np.sqrt(2 * np.pi))
        return xi, yi

    @staticmethod
    def lorentzian_ecd(x, y, xmin, xmax, xstep, gamma):
        xi = np.arange(xmin, xmax, xstep)
        yi = np.zeros(len(xi))
        for i in range(len(xi)):
            for k in range(len(x)):
                yi[i] = yi[i] + y[k] * (gamma) / ((xi[i] - x[k])**2 +
                                                  (gamma)**2)
        return xi, yi

    @staticmethod
    def gaussian_ecd(x, y, xmin, xmax, xstep, sigma):
        xi = np.arange(xmin, xmax, xstep)
        yi = np.zeros(len(xi))
        for i in range(len(xi)):
            for k in range(len(x)):
                yi[i] = yi[i] + y[k] * np.exp(-((xi[i] - x[k])**2) /
                                              (2 * sigma**2))
        yi = np.pi * yi / (sigma * np.sqrt(2 * np.pi))
        return xi, yi

    def plot(self,
             rpa_results,
             broadening_type="lorentzian",
             broadening_value=(1000.0 / hartree_in_wavenumber() *
                               hartree_in_ev()),
             plot_type="electronic"):
        """
        Plot the UV or ECD spectrum from the response calculation.

        :param rpa_results:
            The dictionary containing RPA results.
        :param broadening_type:
            The type of broadening to use. Either 'lorentzian' or 'gaussian'.
        :param broadening_value:
            The broadening value in eV.
        :param plot_type:
            The type of plot to generate. Either 'uv', 'ecd' or 'electronic'.
        """

        assert_msg_critical('matplotlib' in sys.modules,
                            'matplotlib is required.')

        if plot_type.lower() in ["uv", "uv-vis", "uv_vis"]:
            self.plot_uv_vis(rpa_results,
                             broadening_type=broadening_type,
                             broadening_value=broadening_value)

        elif plot_type.lower() == "ecd":
            self.plot_ecd(rpa_results,
                          broadening_type=broadening_type,
                          broadening_value=broadening_value)

        elif plot_type.lower() == "electronic":
            fig, axs = plt.subplots(2, 1, figsize=(8, 10))
            # Increase the height space between subplots
            fig.subplots_adjust(hspace=0.3)
            self.plot_uv_vis(rpa_results,
                             broadening_type=broadening_type,
                             broadening_value=broadening_value,
                             ax=axs[0])
            self.plot_ecd(rpa_results,
                          broadening_type=broadening_type,
                          broadening_value=broadening_value,
                          ax=axs[1])

        else:
            assert_msg_critical(False, 'Invalid plot type')

        plt.show()

    @staticmethod
    def get_absorption_spectrum(rsp_results, x_data, x_unit, b_value, b_unit):
        """
        Gets absorption spectrum.

        :param rsp_results:
            A dictonary containing the result of response calculation.
        :param x_data:
            The list or array of x values.
        :param x_unit:
            The unit of x values.
        :param b_value:
            The value of the broadening parameter.
        :param b_unit:
            The unit of the broadening parameter.

        :return:
            A dictionary containing the spectrum.
        """

        assert_msg_critical(
            x_unit.lower() in ['au', 'ev', 'nm'],
            'LinearResponseUnrestrictedEigenSolver.get_absorption_spectrum: ' +
            'x_data should be au, ev or nm')

        assert_msg_critical(
            b_unit.lower() in ['au', 'ev'],
            'LinearResponseUnrestrictedEigenSolver.get_absorption_spectrum: ' +
            'broadening parameter should be au or ev')

        au2ev = hartree_in_ev()
        auxnm = 1.0 / hartree_in_inverse_nm()

        exc_ene_au = rsp_results['eigenvalues']
        osc_str = rsp_results['oscillator_strengths']

        spectrum = {}

        if x_unit.lower() == 'au':
            spectrum['x_label'] = 'Photon energy [a.u.]'
        elif x_unit.lower() == 'ev':
            spectrum['x_label'] = 'Photon energy [eV]'
        elif x_unit.lower() == 'nm':
            spectrum['x_label'] = 'Wavelength [nm]'

        spectrum['y_label'] = 'Absorption cross-section [a.u.]'

        if x_unit.lower() == 'au':
            x_data_au = list(x_data)
        elif x_unit.lower() == 'ev':
            x_data_au = [x / au2ev for x in x_data]
        elif x_unit.lower() == 'nm':
            x_data_au = [auxnm / x for x in x_data]

        if b_unit.lower() == 'au':
            b_au = b_value
        elif b_unit.lower() == 'ev':
            b_au = b_value / au2ev

        y_data = []

        sigma_factor = 2.0 * np.pi * fine_structure_constant()

        for x_au in x_data_au:
            y = 0.0
            for e, f in zip(exc_ene_au, osc_str):
                b_factor = b_au / ((e - x_au)**2 + b_au**2)
                y += sigma_factor * b_factor * f
            y_data.append(y)

        spectrum['x_data'] = list(x_data)
        spectrum['y_data'] = y_data

        return spectrum

    @staticmethod
    def get_ecd_spectrum(rsp_results, x_data, x_unit, b_value, b_unit):
        """
        Gets ECD spectrum.

        :param rsp_results:
            A dictonary containing the result of response calculation.
        :param x_data:
            The list or array of x values.
        :param x_unit:
            The unit of x values.
        :param b_value:
            The value of the broadening parameter.
        :param b_unit:
            The unit of the broadening parameter.

        :return:
            A dictionary containing the spectrum.
        """

        assert_msg_critical(
            x_unit.lower() in ['au', 'ev', 'nm'],
            'LinearResponseUnrestrictedEigenSolver.get_ecd_spectrum: ' +
            'x_data should be au, ev or nm')

        assert_msg_critical(
            b_unit.lower() in ['au', 'ev'],
            'LinearResponseUnrestrictedEigenSolver.get_ecd_spectrum: ' +
            'broadening parameter should be au or ev')

        au2ev = hartree_in_ev()
        auxnm = 1.0 / hartree_in_inverse_nm()

        exc_ene_au = rsp_results['eigenvalues']
        rot_str_au = rsp_results[
            'rotatory_strengths'] / rotatory_strength_in_cgs()

        spectrum = {}

        if x_unit.lower() == 'au':
            spectrum['x_label'] = 'Photon energy [a.u.]'
        elif x_unit.lower() == 'ev':
            spectrum['x_label'] = 'Photon energy [eV]'
        elif x_unit.lower() == 'nm':
            spectrum['x_label'] = 'Wavelength [nm]'

        spectrum['y_label'] = 'Molar circular dichroism '
        spectrum['y_label'] += '[L mol$^{-1}$ cm$^{-1}$]'

        if x_unit.lower() == 'au':
            x_data_au = list(x_data)
        elif x_unit.lower() == 'ev':
            x_data_au = [x / au2ev for x in x_data]
        elif x_unit.lower() == 'nm':
            x_data_au = [auxnm / x for x in x_data]

        if b_unit.lower() == 'au':
            b_au = b_value
        elif b_unit.lower() == 'ev':
            b_au = b_value / au2ev

        y_data = []

        delta_eps_factor = extinction_coefficient_from_beta() / 3.0

        for x_au in x_data_au:
            y = 0.0
            for e, r in zip(exc_ene_au, rot_str_au):
                b_factor = b_au / ((e - x_au)**2 + b_au**2)
                y += delta_eps_factor * b_factor * e * r
            y_data.append(y)

        spectrum['x_data'] = list(x_data)
        spectrum['y_data'] = y_data

        return spectrum
