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

import numpy as np
import time as tm

from .veloxchemlib import (mpi_master, fine_structure_constant,
                           extinction_coefficient_from_beta)
from .profiler import Profiler
from .distributedarray import DistributedArray
from .cppsolverbase import ComplexResponseSolverBase
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           ri_sanity_check, dft_sanity_check, pe_sanity_check,
                           solvation_model_sanity_check)
from .errorhandler import assert_msg_critical
from .mathutils import safe_solve
from .checkpoint import (check_rsp_hdf5, write_rsp_solution_with_multiple_keys)
from .inputparser import parse_seq_fixed


class ComplexResponseUnrestrictedSolver(ComplexResponseSolverBase):
    """
    Implements the complex linear response solver.

    # vlxtag: UHF, Absorption, CPP
    # vlxtag: UKS, Absorption, CPP
    # vlxtag: UHF, ECD, CPP
    # vlxtag: UKS, ECD, CPP
    """

    def _get_precond(self, orb_ene, nocc, norb, w, d):
        """
        Constructs the preconditioners.
        """

        orb_ene_a, orb_ene_b = orb_ene
        nocc_a, nocc_b = nocc

        # spawning needed components
        ediag_a, sdiag_a = self.construct_ediag_sdiag_half(
            orb_ene_a, nocc_a, norb)
        ediag_b, sdiag_b = self.construct_ediag_sdiag_half(
            orb_ene_b, nocc_b, norb)

        ediag_a_sq = ediag_a**2
        sdiag_a_sq = sdiag_a**2
        sdiag_a_fp = sdiag_a**4

        ediag_b_sq = ediag_b**2
        sdiag_b_sq = sdiag_b**2
        sdiag_b_fp = sdiag_b**4

        w_sq = w**2
        d_sq = d**2

        # constructing matrix block diagonals
        a_alpha_diag = ediag_a * (ediag_a_sq - (w_sq - d_sq) * sdiag_a_sq)
        b_alpha_diag = (w * sdiag_a) * (ediag_a_sq - (w_sq + d_sq) * sdiag_a_sq)
        c_alpha_diag = (d * sdiag_a) * (ediag_a_sq + (w_sq + d_sq) * sdiag_a_sq)
        d_alpha_diag = (2 * w * d * ediag_a) * sdiag_a_sq
        p_alpha_diag = 1.0 / ((ediag_a_sq - (w_sq - d_sq) * sdiag_a_sq)**2 +
                              (4 * w_sq * d_sq * sdiag_a_fp))

        a_beta_diag = ediag_b * (ediag_b_sq - (w_sq - d_sq) * sdiag_b_sq)
        b_beta_diag = (w * sdiag_b) * (ediag_b_sq - (w_sq + d_sq) * sdiag_b_sq)
        c_beta_diag = (d * sdiag_b) * (ediag_b_sq + (w_sq + d_sq) * sdiag_b_sq)
        d_beta_diag = (2 * w * d * ediag_b) * sdiag_b_sq
        p_beta_diag = 1.0 / ((ediag_b_sq - (w_sq - d_sq) * sdiag_b_sq)**2 +
                             (4 * w_sq * d_sq * sdiag_b_fp))

        pa_alpha_diag = p_alpha_diag * a_alpha_diag
        pb_alpha_diag = p_alpha_diag * b_alpha_diag
        pc_alpha_diag = p_alpha_diag * c_alpha_diag
        pd_alpha_diag = p_alpha_diag * d_alpha_diag

        pa_beta_diag = p_beta_diag * a_beta_diag
        pb_beta_diag = p_beta_diag * b_beta_diag
        pc_beta_diag = p_beta_diag * c_beta_diag
        pd_beta_diag = p_beta_diag * d_beta_diag

        p_mat_alpha = np.hstack((
            pa_alpha_diag.reshape(-1, 1),
            pb_alpha_diag.reshape(-1, 1),
            pc_alpha_diag.reshape(-1, 1),
            pd_alpha_diag.reshape(-1, 1),
        ))

        p_mat_beta = np.hstack((
            pa_beta_diag.reshape(-1, 1),
            pb_beta_diag.reshape(-1, 1),
            pc_beta_diag.reshape(-1, 1),
            pd_beta_diag.reshape(-1, 1),
        ))

        # put alpha and beta together
        p_mat = np.vstack((p_mat_alpha, p_mat_beta))

        return DistributedArray(p_mat, self.comm)

    def compute(self, molecule, basis, scf_results, v_grad=None):
        """
        Solves for the response vector iteratively while checking the residuals
        for convergence.
        """

        # take care of quadrupole components
        if self.is_quadrupole(self.a_operator):
            if isinstance(self.a_components, str):
                self.a_components = parse_seq_fixed(self.a_components, 'str')
        if self.is_quadrupole(self.b_operator):
            if isinstance(self.b_components, str):
                self.b_components = parse_seq_fixed(self.b_components, 'str')

        # check operator components
        for comp in self.a_components:
            assert_msg_critical(
                self.is_valid_component(comp, self.a_operator),
                f'{type(self).__name__}: Undefined or invalid a_component')
        for comp in self.b_components:
            assert_msg_critical(
                self.is_valid_component(comp, self.b_operator),
                f'{type(self).__name__}: Undefined or invalid b_component')

        if self.norm_thresh is None:
            self.norm_thresh = self.conv_thresh * 1.0e-6
        if self.lindep_thresh is None:
            self.lindep_thresh = self.conv_thresh * 1.0e-2

        self._dist_bger = None
        self._dist_bung = None
        self._dist_e2bger = None
        self._dist_e2bung = None

        self.nonlinear = False
        self._dist_fock_ger = None
        self._dist_fock_ung = None

        # make sure that cpp_property is properly set
        if self.property is not None:
            self.set_cpp_property(self.property)

        # check molecule
        molecule_sanity_check(molecule)
        # check SCF results
        scf_results_sanity_check(self, scf_results)

        # update checkpoint_file after scf_results_sanity_check
        if self.filename is not None and self.checkpoint_file is None:
            self.checkpoint_file = f'{self.filename}_rsp.h5'

        # check RI setup
        ri_sanity_check(self)
        # check dft setup
        dft_sanity_check(self, 'compute')
        # check pe setup
        pe_sanity_check(self, molecule=molecule)
        # check solvation setup
        solvation_model_sanity_check(self)

        # check print level (verbosity of output)
        self.print_level = max(1, min(self.print_level, 3))

        # initialize profiler
        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        if self.rank == mpi_master():
            self._print_header('Complex Response Solver',
                               n_freqs=len(self.frequencies))

        self.start_time = tm.time()

        if self.rank == mpi_master():
            orb_ene_a = scf_results['E_alpha']
            orb_ene_b = scf_results['E_beta']
        else:
            orb_ene_a = None
            orb_ene_b = None

        orb_ene_a = self.comm.bcast(orb_ene_a, root=mpi_master())
        orb_ene_b = self.comm.bcast(orb_ene_b, root=mpi_master())

        norb = orb_ene_a.shape[0]

        nocc_a = molecule.number_of_alpha_occupied_orbitals(basis)
        nocc_b = molecule.number_of_beta_occupied_orbitals(basis)

        self._check_mpi_oversubscription(
            self._get_excitation_space_dimension_unrestricted(
                nocc_a, nocc_b, norb), 'response space')

        # ERI information
        eri_dict = self._init_eri(molecule, basis)
        # DFT information
        dft_dict = self._init_dft(molecule, scf_results)
        # PE information
        pe_dict = self._init_pe(molecule, basis)
        # CPCM information
        self._init_cpcm(molecule, basis)

        # TODO: enable PE
        assert_msg_critical(
            not self._pe, f'{type(self).__name__}: ' +
            'not yet implemented for polarizable embedding')

        # right-hand side (gradient)
        if self.rank == mpi_master():
            self.nonlinear = (v_grad is not None)
        self.nonlinear = self.comm.bcast(self.nonlinear, root=mpi_master())

        # For now, 'nonlinear' is not supported for unrestricted case.
        assert_msg_critical(
            not self.nonlinear,
            f'{type(self).__name__}: ' + 'not implemented for nonlinear')

        sqrt_2 = np.sqrt(2.0)

        if not self.nonlinear:
            b_grad_alpha = self.get_complex_prop_grad(self.b_operator,
                                                      self.b_components,
                                                      molecule,
                                                      basis,
                                                      scf_results,
                                                      spin='alpha')
            b_grad_beta = self.get_complex_prop_grad(self.b_operator,
                                                     self.b_components,
                                                     molecule,
                                                     basis,
                                                     scf_results,
                                                     spin='beta')

            # for unrestricted
            b_grad_alpha /= sqrt_2
            b_grad_beta /= sqrt_2

            if self.rank == mpi_master():
                v_grad = {
                    (op, w): (va, vb) for op, va, vb in zip(
                        self.b_components, b_grad_alpha, b_grad_beta)
                    for w in self.frequencies
                }

        # operators, frequencies and preconditioners
        if self.rank == mpi_master():
            op_freq_keys = list(v_grad.keys())
        else:
            op_freq_keys = None
        op_freq_keys = self.comm.bcast(op_freq_keys, root=mpi_master())

        d = self.damping

        self.frequencies = []
        for (op, w) in op_freq_keys:
            if w not in self.frequencies:
                self.frequencies.append(w)

        precond = {
            w: self._get_precond((orb_ene_a, orb_ene_b), (nocc_a, nocc_b), norb,
                                 w, d) for w in self.frequencies
        }

        # distribute the gradient and right-hand side:
        # dist_grad will be used for calculating the subspace matrix
        # equation and residuals, dist_rhs for the initial guess
        dist_grad = {}
        dist_rhs = {}
        for key in op_freq_keys:
            if self.rank == mpi_master():
                gradger_a, gradung_a = self._decomp_grad(v_grad[key][0])
                gradger_b, gradung_b = self._decomp_grad(v_grad[key][1])

                grad_mat_a = np.hstack((
                    gradger_a.real.reshape(-1, 1),
                    gradung_a.real.reshape(-1, 1),
                    gradung_a.imag.reshape(-1, 1),
                    gradger_a.imag.reshape(-1, 1),
                ))
                rhs_mat_a = np.hstack((
                    gradger_a.real.reshape(-1, 1),
                    gradung_a.real.reshape(-1, 1),
                    -gradung_a.imag.reshape(-1, 1),
                    -gradger_a.imag.reshape(-1, 1),
                ))

                grad_mat_b = np.hstack((
                    gradger_b.real.reshape(-1, 1),
                    gradung_b.real.reshape(-1, 1),
                    gradung_b.imag.reshape(-1, 1),
                    gradger_b.imag.reshape(-1, 1),
                ))
                rhs_mat_b = np.hstack((
                    gradger_b.real.reshape(-1, 1),
                    gradung_b.real.reshape(-1, 1),
                    -gradung_b.imag.reshape(-1, 1),
                    -gradger_b.imag.reshape(-1, 1),
                ))

                # put alpha and beta together
                grad_mat = np.vstack((grad_mat_a, grad_mat_b))
                rhs_mat = np.vstack((rhs_mat_a, rhs_mat_b))
            else:
                grad_mat = None
                rhs_mat = None
            dist_grad[key] = DistributedArray(grad_mat, self.comm)
            dist_rhs[key] = DistributedArray(rhs_mat, self.comm)

        if self.nonlinear:
            rsp_vector_labels = [
                'CLR_bger_half_size', 'CLR_bung_half_size',
                'CLR_e2bger_half_size', 'CLR_e2bung_half_size', 'CLR_Fock_ger',
                'CLR_Fock_ung'
            ]
        else:
            rsp_vector_labels = [
                'CLR_bger_half_size', 'CLR_bung_half_size',
                'CLR_e2bger_half_size', 'CLR_e2bung_half_size'
            ]

        # check validity of checkpoint file
        if self.restart:
            if self.rank == mpi_master():
                self.restart = check_rsp_hdf5(self.checkpoint_file,
                                              rsp_vector_labels, molecule,
                                              basis, dft_dict, pe_dict)
                if self.restart:
                    self.restart = self.match_settings(self.checkpoint_file)
            self.restart = self.comm.bcast(self.restart, root=mpi_master())

        # read initial guess from restart file
        if self.restart:
            self._read_checkpoint(rsp_vector_labels)

            checkpoint_frequencies = self._read_frequencies_from_checkpoint()

            # Note that we only handle the restart with additional frequencies when
            # `self.nonlinear` is inactive

            # print warning if frequencies is not present in the restart file
            if checkpoint_frequencies is None and not self.nonlinear:
                self.ostream.print_warning(
                    'Could not find the frequencies key in the checkpoint file.'
                )
                self.ostream.print_blank()
                self.ostream.print_info(
                    'Assuming that frequencies is not changed before and after '
                    + 'the restart.')
                self.ostream.print_blank()
                self.ostream.flush()

            # generate necessary initial guesses if more frequencies are requested
            # in a restart calculation
            elif not self.nonlinear:
                extra_freqs = []
                for w in self.frequencies:
                    if w not in checkpoint_frequencies:
                        extra_freqs.append(w)

                if extra_freqs:
                    self.ostream.print_info(
                        f'Generating initial guesses for {len(extra_freqs)} ' +
                        'more frequencies...')
                    self.ostream.print_blank()

                    if self.rank == mpi_master():
                        extra_v_grad = {
                            (op, w): (va, vb) for op, va, vb in zip(
                                self.b_components, b_grad_alpha, b_grad_beta)
                            for w in extra_freqs
                        }
                        extra_op_freq_keys = list(extra_v_grad.keys())
                    else:
                        extra_op_freq_keys = None
                    extra_op_freq_keys = self.comm.bcast(extra_op_freq_keys,
                                                         root=mpi_master())

                    extra_precond = {
                        w: self._get_precond((orb_ene_a, orb_ene_b),
                                             (nocc_a, nocc_b), norb, w,
                                             d) for w in extra_freqs
                    }

                    extra_dist_rhs = {}
                    for key in extra_op_freq_keys:
                        if self.rank == mpi_master():
                            gradger_a, gradung_a = self._decomp_grad(
                                extra_v_grad[key][0])
                            gradger_b, gradung_b = self._decomp_grad(
                                extra_v_grad[key][1])

                            grad_mat_a = np.hstack((
                                gradger_a.real.reshape(-1, 1),
                                gradung_a.real.reshape(-1, 1),
                                gradung_a.imag.reshape(-1, 1),
                                gradger_a.imag.reshape(-1, 1),
                            ))
                            rhs_mat_a = np.hstack((
                                gradger_a.real.reshape(-1, 1),
                                gradung_a.real.reshape(-1, 1),
                                -gradung_a.imag.reshape(-1, 1),
                                -gradger_a.imag.reshape(-1, 1),
                            ))

                            grad_mat_b = np.hstack((
                                gradger_b.real.reshape(-1, 1),
                                gradung_b.real.reshape(-1, 1),
                                gradung_b.imag.reshape(-1, 1),
                                gradger_b.imag.reshape(-1, 1),
                            ))
                            rhs_mat_b = np.hstack((
                                gradger_b.real.reshape(-1, 1),
                                gradung_b.real.reshape(-1, 1),
                                -gradung_b.imag.reshape(-1, 1),
                                -gradger_b.imag.reshape(-1, 1),
                            ))

                            # put alpha and beta together
                            grad_mat = np.vstack((grad_mat_a, grad_mat_b))
                            rhs_mat = np.vstack((rhs_mat_a, rhs_mat_b))
                        else:
                            grad_mat = None
                            rhs_mat = None
                        dist_grad[key] = DistributedArray(grad_mat, self.comm)
                        extra_dist_rhs[key] = DistributedArray(
                            rhs_mat, self.comm)

                    bger, bung = self._setup_trials(extra_dist_rhs,
                                                    extra_precond)

                    profiler.set_timing_key('Preparation')

                    self._e2n_half_size(bger,
                                        bung,
                                        molecule,
                                        basis,
                                        scf_results,
                                        eri_dict,
                                        dft_dict,
                                        pe_dict,
                                        profiler,
                                        method_type='unrestricted')

        # generate initial guess from scratch
        else:
            bger, bung = self._setup_trials(dist_rhs, precond)

            profiler.set_timing_key('Preparation')

            self._e2n_half_size(bger,
                                bung,
                                molecule,
                                basis,
                                scf_results,
                                eri_dict,
                                dft_dict,
                                pe_dict,
                                profiler,
                                method_type='unrestricted')

        profiler.check_memory_usage('Initial guess')

        focks = {}
        solutions = {}
        residuals = {}
        relative_residual_norm = {}

        iter_per_trial_in_hours = None

        # start iterations
        for iteration in range(self.max_iter):

            iter_start_time = tm.time()

            profiler.set_timing_key(f'Iteration {iteration + 1}')
            profiler.start_timer('ReducedSpace')

            xvs = []
            current_solutions = {}
            self._cur_iter = iteration
            self.collapsed_subspace = False
            self.collapsed_from_dim = None
            self.collapsed_to_dim = None

            n_ger = self._dist_bger.shape(1)
            n_ung = self._dist_bung.shape(1)

            e2gg = self._dist_bger.matmul_AtB(self._dist_e2bger)
            e2uu = self._dist_bung.matmul_AtB(self._dist_e2bung)
            s2ug = self._dist_bung.matmul_AtB(self._dist_bger)

            for op, w in op_freq_keys:
                if (iteration == 0 or
                        relative_residual_norm[(op, w)] > self.conv_thresh):
                    solve_data = self._solve_reduced_problem(
                        op, w, n_ger, n_ung, e2gg, e2uu, s2ug, dist_grad,
                        v_grad, nocc_a, nocc_b, norb)
                    x = solve_data['solution']
                    r = solve_data['residual']
                    current_solutions[(op, w)] = x

                    if self.nonlinear:
                        focks[(op, w)] = solve_data['fock']

                    if self.rank == mpi_master():
                        xvs.append((op, w, solve_data['property']))

                    r_norms_2 = r.squared_norm(axis=0)
                    x_norms_2 = x.squared_norm(axis=0)

                    rn = np.sqrt(np.sum(r_norms_2))
                    xn = np.sqrt(np.sum(x_norms_2))

                    if xn != 0:
                        relative_residual_norm[(op, w)] = 2.0 * rn / xn
                    else:
                        relative_residual_norm[(op, w)] = 2.0 * rn

                    if relative_residual_norm[(op, w)] < self.conv_thresh:
                        solutions[(op, w)] = x
                    else:
                        residuals[(op, w)] = r

            # write to output
            if self.rank == mpi_master():

                self.ostream.print_info(
                    '{:d} gerade trial vectors in reduced space'.format(n_ger))
                self.ostream.print_info(
                    '{:d} ungerade trial vectors in reduced space'.format(
                        n_ung))
                self.ostream.print_blank()

                if self.print_level > 1:
                    profiler.print_memory_subspace(
                        {
                            'dist_bger': self._dist_bger,
                            'dist_bung': self._dist_bung,
                            'dist_e2bger': self._dist_e2bger,
                            'dist_e2bung': self._dist_e2bung,
                            'precond': precond,
                            'solutions': solutions,
                            'residuals': residuals,
                        }, self.ostream)

                profiler.check_memory_usage(
                    'Iteration {:d} subspace'.format(iteration + 1))

                profiler.print_memory_tracing(self.ostream)
                self._print_iteration(relative_residual_norm, xvs)

            profiler.stop_timer('ReducedSpace')

            # check convergence
            self._check_convergence(relative_residual_norm)

            if self.is_converged:
                break

            active_keys = list(residuals.keys())
            if active_keys and self._should_collapse_subspace():
                self.ostream.print_info('Collapsing reduced space...')
                self.ostream.print_blank()

                self._collapse_current_subspace(active_keys, current_solutions,
                                                relative_residual_norm,
                                                molecule, basis, scf_results,
                                                eri_dict, dft_dict, pe_dict,
                                                profiler)

                collapse_str = 'Collapsed reduced space: {:d}->{:d}'.format(
                    self.collapsed_from_dim, self.collapsed_to_dim)
                self.ostream.print_info(collapse_str)
                self.ostream.print_blank()
                self.ostream.flush()

            profiler.start_timer('Orthonorm.')

            # spawning new trial vectors from residuals
            new_trials_ger, new_trials_ung = self._setup_trials(
                residuals, precond, self._dist_bger, self._dist_bung)

            residuals.clear()

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
                self._add_frequencies_to_checkpoint()

            # creating new sigma and rho linear transformations
            self._e2n_half_size(new_trials_ger,
                                new_trials_ung,
                                molecule,
                                basis,
                                scf_results,
                                eri_dict,
                                dft_dict,
                                pe_dict,
                                profiler,
                                method_type='unrestricted')

            iter_in_hours = (tm.time() - iter_start_time) / 3600
            iter_per_trial_in_hours = iter_in_hours / n_new_trials

            profiler.check_memory_usage(
                'Iteration {:d} sigma build'.format(iteration + 1))

        self._write_checkpoint(molecule, basis, dft_dict, pe_dict,
                               rsp_vector_labels)
        self._add_frequencies_to_checkpoint()

        # converged?
        if self.rank == mpi_master():
            self._print_convergence('Complex response')

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of CPP solver')
        profiler.print_memory_usage(self.ostream)

        self._dist_bger = None
        self._dist_bung = None
        self._dist_e2bger = None
        self._dist_e2bung = None

        self._dist_fock_ger = None
        self._dist_fock_ung = None

        # calculate response functions
        if not self.nonlinear:
            a_grad_alpha = self.get_complex_prop_grad(self.a_operator,
                                                      self.a_components,
                                                      molecule,
                                                      basis,
                                                      scf_results,
                                                      spin='alpha')
            a_grad_beta = self.get_complex_prop_grad(self.a_operator,
                                                     self.a_components,
                                                     molecule,
                                                     basis,
                                                     scf_results,
                                                     spin='beta')

            # for unrestricted
            a_grad_alpha /= sqrt_2
            a_grad_beta /= sqrt_2

            if self.is_converged:
                if self.rank == mpi_master():
                    va = {
                        op: (va, vb) for op, va, vb in zip(
                            self.a_components, a_grad_alpha, a_grad_beta)
                    }
                    rsp_funcs = {}

                    # final h5 file for response solutions
                    if self.filename is not None:
                        final_h5_fname = f'{self.filename}.h5'
                    else:
                        final_h5_fname = None

                for bop, w in solutions:
                    x = self.get_full_solution_vector(solutions[(bop, w)])

                    if self.rank == mpi_master():
                        n_ov_a = nocc_a * (norb - nocc_a)
                        n_ov_b = nocc_b * (norb - nocc_b)

                        x_alpha = np.hstack((
                            x[:n_ov_a],
                            x[n_ov_a + n_ov_b:n_ov_a + n_ov_b + n_ov_a],
                        ))
                        x_beta = np.hstack((
                            x[n_ov_a:n_ov_a + n_ov_b],
                            x[n_ov_a + n_ov_b + n_ov_a:],
                        ))

                        for aop in self.a_components:
                            rsp_funcs[(aop, bop,
                                       w)] = -(np.dot(va[aop][0], x_alpha) +
                                               np.dot(va[aop][1], x_beta))

                            # Note: flip sign for imaginary a_operator
                            if self.is_imag(self.a_operator):
                                rsp_funcs[(aop, bop, w)] *= -1.0

                        # write to h5 file for response solutions
                        if (self.save_solutions and final_h5_fname is not None):
                            solution_keys = [
                                '{:s}_{:s}_{:.8f}'.format(aop, bop, w)
                                for aop in self.a_components
                            ]
                            write_rsp_solution_with_multiple_keys(
                                final_h5_fname, solution_keys, x,
                                self.group_label)

                if self.rank == mpi_master():
                    # print information about h5 file for response solutions
                    if (self.save_solutions and final_h5_fname is not None):
                        self.ostream.print_info(
                            'Response solution vectors written to file: ' +
                            final_h5_fname)
                        self.ostream.print_blank()

                    ret_dict = {
                        'a_operator': self.a_operator,
                        'a_components': self.a_components,
                        'b_operator': self.b_operator,
                        'b_components': self.b_components,
                        'frequencies': list(self.frequencies),
                        'response_functions': rsp_funcs,
                        'solutions': solutions,
                    }

                    self._print_results(ret_dict)

                    # write spectrum to h5 file
                    if final_h5_fname is not None:
                        h5_ret_dict = {
                            key: value
                            for key, value in ret_dict.items()
                            if key != 'solutions'
                        }
                        self.write_cpp_rsp_results_to_hdf5(
                            final_h5_fname, h5_ret_dict)

                    return ret_dict
                else:
                    # non-master rank
                    return {'solutions': solutions}
            else:
                # not converged
                return {}

        else:
            # nonlinear
            if self.is_converged:
                return {'focks': focks, 'solutions': solutions}
            else:
                return {}

    def _solve_reduced_problem(self, op, w, n_ger, n_ung, e2gg, e2uu, s2ug,
                               dist_grad, v_grad, nocc_a, nocc_b, norb):
        """
        Solves one reduced unrestricted CPP linear system and forms the
        corresponding solution, residual, and projected property.
        """

        grad_rg = dist_grad[(op, w)].get_column(0)
        grad_ru = dist_grad[(op, w)].get_column(1)
        grad_iu = dist_grad[(op, w)].get_column(2)
        grad_ig = dist_grad[(op, w)].get_column(3)

        g_realger = self._dist_bger.matmul_AtB(grad_rg)
        g_imagger = self._dist_bger.matmul_AtB(grad_ig)
        g_realung = self._dist_bung.matmul_AtB(grad_ru)
        g_imagung = self._dist_bung.matmul_AtB(grad_iu)

        size = 2 * (n_ger + n_ung)

        if self.rank == mpi_master():
            g = np.zeros(size)

            g[:n_ger] = g_realger[:]
            g[n_ger:n_ger + n_ung] = g_realung[:]
            g[n_ger + n_ung:size - n_ger] = -g_imagung[:]
            g[size - n_ger:] = -g_imagger[:]

            mat = np.zeros((size, size))

            mat[:n_ger, :n_ger] = e2gg[:, :]
            mat[size - n_ger:, size - n_ger:] = -e2gg[:, :]

            mat[n_ger:n_ger + n_ung, n_ger:n_ger + n_ung] = e2uu[:, :]
            mat[n_ger + n_ung:size - n_ger,
                n_ger + n_ung:size - n_ger] = -e2uu[:, :]

            mat[n_ger:n_ger + n_ung, :n_ger] = -w * s2ug[:, :]
            mat[n_ger + n_ung:size - n_ger, :n_ger] = self.damping * s2ug[:, :]
            mat[n_ger:n_ger + n_ung, size - n_ger:] = self.damping * s2ug[:, :]
            mat[n_ger + n_ung:size - n_ger, size - n_ger:] = w * s2ug[:, :]

            mat[:n_ger, n_ger:n_ger + n_ung] = -w * s2ug.T[:, :]
            mat[:n_ger,
                n_ger + n_ung:size - n_ger] = self.damping * s2ug.T[:, :]
            mat[size - n_ger:,
                n_ger:n_ger + n_ung] = self.damping * s2ug.T[:, :]
            mat[size - n_ger:, n_ger + n_ung:size - n_ger] = w * s2ug.T[:, :]

            c = safe_solve(mat, g)
        else:
            c = None
        c = self.comm.bcast(c, root=mpi_master())

        c_realger = c[:n_ger]
        c_realung = c[n_ger:n_ger + n_ung]
        c_imagung = c[n_ger + n_ung:size - n_ger]
        c_imagger = c[size - n_ger:]

        x_realger = self._dist_bger.matmul_AB_no_gather(c_realger)
        x_realung = self._dist_bung.matmul_AB_no_gather(c_realung)
        x_imagung = self._dist_bung.matmul_AB_no_gather(c_imagung)
        x_imagger = self._dist_bger.matmul_AB_no_gather(c_imagger)

        e2realger = self._dist_e2bger.matmul_AB_no_gather(c_realger)
        e2imagger = self._dist_e2bger.matmul_AB_no_gather(c_imagger)
        e2realung = self._dist_e2bung.matmul_AB_no_gather(c_realung)
        e2imagung = self._dist_e2bung.matmul_AB_no_gather(c_imagung)

        fock = None
        if self.nonlinear:
            fock_realger = self._dist_fock_ger.matmul_AB_no_gather(c_realger)
            fock_imagger = self._dist_fock_ger.matmul_AB_no_gather(c_imagger)
            fock_realung = self._dist_fock_ung.matmul_AB_no_gather(c_realung)
            fock_imagung = self._dist_fock_ung.matmul_AB_no_gather(c_imagung)

            fock_full_data = (fock_realger.data + fock_realung.data - 1j *
                              (fock_imagger.data + fock_imagung.data))

            fock = DistributedArray(fock_full_data, self.comm, distribute=False)

        s2realger = x_realger.data
        s2imagger = x_imagger.data
        s2realung = x_realung.data
        s2imagung = x_imagung.data

        r_realger = (e2realger.data - w * s2realung + self.damping * s2imagung -
                     grad_rg.data)
        r_realung = (e2realung.data - w * s2realger + self.damping * s2imagger -
                     grad_ru.data)
        r_imagung = (-e2imagung.data + w * s2imagger +
                     self.damping * s2realger + grad_iu.data)
        r_imagger = (-e2imagger.data + w * s2imagung +
                     self.damping * s2realung + grad_ig.data)

        r_data = np.hstack((
            r_realger.reshape(-1, 1),
            r_realung.reshape(-1, 1),
            r_imagung.reshape(-1, 1),
            r_imagger.reshape(-1, 1),
        ))

        r = DistributedArray(r_data, self.comm, distribute=False)

        x_data = np.hstack((
            x_realger.data.reshape(-1, 1),
            x_realung.data.reshape(-1, 1),
            x_imagung.data.reshape(-1, 1),
            x_imagger.data.reshape(-1, 1),
        ))

        x = DistributedArray(x_data, self.comm, distribute=False)

        x_full = self.get_full_solution_vector(x)
        if self.rank == mpi_master():
            n_ov_a = nocc_a * (norb - nocc_a)
            n_ov_b = nocc_b * (norb - nocc_b)

            x_full_a = np.hstack((
                x_full[:n_ov_a],
                x_full[n_ov_a + n_ov_b:n_ov_a + n_ov_b + n_ov_a],
            ))
            x_full_b = np.hstack((
                x_full[n_ov_a:n_ov_a + n_ov_b],
                x_full[n_ov_a + n_ov_b + n_ov_a:],
            ))

            xv = (np.dot(x_full_a, v_grad[(op, w)][0]) +
                  np.dot(x_full_b, v_grad[(op, w)][1]))
        else:
            xv = None

        return {
            'solution': x,
            'residual': r,
            'property': xv,
            'fock': fock,
        }

    def _collapse_current_subspace(self, active_keys, solutions, residual_norms,
                                   molecule, basis, scf_results, eri_dict,
                                   dft_dict, pe_dict, profiler):
        """
        Collapses the reduced space to a basis built from the largest
        unconverged solution vectors and rebuilds associated sigma data.
        """

        ranked_keys = sorted(active_keys,
                             key=lambda key: residual_norms[key],
                             reverse=True)
        keep_keys = ranked_keys[:min(
            self._get_collapse_nvec(),
            len(ranked_keys),
        )]

        prev_dim = self._reduced_space_size()

        ger_cols = []
        ung_cols = []
        for key in keep_keys:
            x = solutions[key]
            ger_cols.append(x.data[:, 0].copy())
            ger_cols.append(x.data[:, 3].copy())
            ung_cols.append(x.data[:, 1].copy())
            ung_cols.append(x.data[:, 2].copy())

        if ger_cols:
            ger_mat = np.array(ger_cols).T
        else:
            ger_mat = np.zeros((self._dist_bger.shape(0), 0))

        if ung_cols:
            ung_mat = np.array(ung_cols).T
        else:
            ung_mat = np.zeros((self._dist_bung.shape(0), 0))

        new_bger = DistributedArray(ger_mat, self.comm, distribute=False)
        new_bung = DistributedArray(ung_mat, self.comm, distribute=False)

        new_bger = self._orthonormalize_collapsed_space(new_bger)
        new_bung = self._orthonormalize_collapsed_space(new_bung)

        self._clear_subspace_data()
        self._e2n_half_size(new_bger,
                            new_bung,
                            molecule,
                            basis,
                            scf_results,
                            eri_dict,
                            dft_dict,
                            pe_dict,
                            profiler,
                            method_type='unrestricted')

        self.collapsed_subspace = True
        self.collapsed_from_dim = prev_dim
        self.collapsed_to_dim = self._reduced_space_size()

    def get_cpp_property_densities(self,
                                   molecule,
                                   basis,
                                   scf_results,
                                   cpp_results,
                                   w,
                                   normalize_densities=True):
        """
        Gets CPP property densities.
        """

        assert_msg_critical(self.property in ['absorption', 'ecd'],
                            'get_cpp_property_densities: Invalid CPP property')

        assert_msg_critical(
            ('x', w) in cpp_results['solutions'] and
            ('y', w) in cpp_results['solutions'] and
            ('z', w) in cpp_results['solutions'],
            f'get_cpp_property_densities: Could not find frequency {w} in ' +
            'CPP results')

        # solution vectors
        cpp_solution_vector_x = self.get_full_solution_vector(
            cpp_results['solutions'][('x', w)])
        cpp_solution_vector_y = self.get_full_solution_vector(
            cpp_results['solutions'][('y', w)])
        cpp_solution_vector_z = self.get_full_solution_vector(
            cpp_results['solutions'][('z', w)])

        # property gradient for a operator
        a_grad_alpha = self.get_complex_prop_grad(self.a_operator,
                                                  self.a_components,
                                                  molecule,
                                                  basis,
                                                  scf_results,
                                                  spin='alpha')
        a_grad_beta = self.get_complex_prop_grad(self.a_operator,
                                                 self.a_components,
                                                 molecule,
                                                 basis,
                                                 scf_results,
                                                 spin='beta')

        # for unrestricted
        sqrt_2 = np.sqrt(2.0)
        a_grad_alpha /= sqrt_2
        a_grad_beta /= sqrt_2

        if self.rank == mpi_master():
            nocc_a = molecule.number_of_alpha_occupied_orbitals(basis)
            nocc_b = molecule.number_of_beta_occupied_orbitals(basis)
            norb = scf_results['E_alpha'].shape[0]

            nvir_a = norb - nocc_a
            n_ov_a = nocc_a * nvir_a

            nvir_b = norb - nocc_b
            n_ov_b = nocc_b * nvir_b

            mo_occ_a = scf_results['C_alpha'][:, :nocc_a].copy()
            mo_vir_a = scf_results['C_alpha'][:, nocc_a:].copy()

            mo_occ_b = scf_results['C_beta'][:, :nocc_b].copy()
            mo_vir_b = scf_results['C_beta'][:, nocc_b:].copy()

            x_alpha = np.hstack((
                cpp_solution_vector_x[:n_ov_a],
                cpp_solution_vector_x[n_ov_a + n_ov_b:n_ov_a + n_ov_b + n_ov_a],
            ))
            x_beta = np.hstack((
                cpp_solution_vector_x[n_ov_a:n_ov_a + n_ov_b],
                cpp_solution_vector_x[n_ov_a + n_ov_b + n_ov_a:],
            ))

            y_alpha = np.hstack((
                cpp_solution_vector_y[:n_ov_a],
                cpp_solution_vector_y[n_ov_a + n_ov_b:n_ov_a + n_ov_b + n_ov_a],
            ))
            y_beta = np.hstack((
                cpp_solution_vector_y[n_ov_a:n_ov_a + n_ov_b],
                cpp_solution_vector_y[n_ov_a + n_ov_b + n_ov_a:],
            ))

            z_alpha = np.hstack((
                cpp_solution_vector_z[:n_ov_a],
                cpp_solution_vector_z[n_ov_a + n_ov_b:n_ov_a + n_ov_b + n_ov_a],
            ))
            z_beta = np.hstack((
                cpp_solution_vector_z[n_ov_a:n_ov_a + n_ov_b],
                cpp_solution_vector_z[n_ov_a + n_ov_b + n_ov_a:],
            ))

            if self.property == 'absorption':
                # vector representation of absorption cross-section
                vec_a = (a_grad_alpha[0] * x_alpha + a_grad_alpha[1] * y_alpha +
                         a_grad_alpha[2] * z_alpha).imag / 3.0
                vec_b = (a_grad_beta[0] * x_beta + a_grad_beta[1] * y_beta +
                         a_grad_beta[2] * z_beta).imag / 3.0
                vec_a *= 4.0 * np.pi * w * fine_structure_constant()
                vec_b *= 4.0 * np.pi * w * fine_structure_constant()
            elif self.property == 'ecd':
                # vector representation of Delta epsilon
                vec_a = (a_grad_alpha[0] * x_alpha / w +
                         a_grad_alpha[1] * y_alpha / w +
                         a_grad_alpha[2] * z_alpha / w).imag
                vec_b = (a_grad_beta[0] * x_beta / w +
                         a_grad_beta[1] * y_beta / w +
                         a_grad_beta[2] * z_beta / w).imag
                vec_a /= (3.0 * w)
                vec_b /= (3.0 * w)
                vec_a *= w**2 * extinction_coefficient_from_beta()
                vec_b *= w**2 * extinction_coefficient_from_beta()

            # excitation and de-excitation
            z_mat_ov_a = vec_a[:n_ov_a].reshape(nocc_a, nvir_a)
            y_mat_ov_a = vec_a[n_ov_a:].reshape(nocc_a, nvir_a)

            prop_diag_D_a = np.sum(z_mat_ov_a, axis=1) + np.sum(y_mat_ov_a,
                                                                axis=1)
            prop_diag_A_a = np.sum(z_mat_ov_a, axis=0) + np.sum(y_mat_ov_a,
                                                                axis=0)

            prop_dens_D_a = -np.linalg.multi_dot(
                [mo_occ_a, np.diag(prop_diag_D_a), mo_occ_a.T])
            prop_dens_A_a = np.linalg.multi_dot(
                [mo_vir_a, np.diag(prop_diag_A_a), mo_vir_a.T])

            z_mat_ov_b = vec_b[:n_ov_b].reshape(nocc_b, nvir_b)
            y_mat_ov_b = vec_b[n_ov_b:].reshape(nocc_b, nvir_b)

            prop_diag_D_b = np.sum(z_mat_ov_b, axis=1) + np.sum(y_mat_ov_b,
                                                                axis=1)
            prop_diag_A_b = np.sum(z_mat_ov_b, axis=0) + np.sum(y_mat_ov_b,
                                                                axis=0)

            prop_dens_D_b = -np.linalg.multi_dot(
                [mo_occ_b, np.diag(prop_diag_D_b), mo_occ_b.T])
            prop_dens_A_b = np.linalg.multi_dot(
                [mo_vir_b, np.diag(prop_diag_A_b), mo_vir_b.T])

            prop_dens_D = prop_dens_D_a + prop_dens_D_b
            prop_dens_A = prop_dens_A_a + prop_dens_A_b

            if normalize_densities:
                abs_sum_val = abs(np.sum(prop_dens_D * scf_results['S']))
                prop_dens_D /= abs_sum_val
                prop_dens_A /= abs_sum_val

            return {
                'property_density_detachment': prop_dens_D,
                'property_density_attachment': prop_dens_A,
            }
        else:
            return None
