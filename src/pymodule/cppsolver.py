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


class ComplexResponseSolver(ComplexResponseSolverBase):
    """
    Implements the complex linear response solver.

    # vlxtag: RHF, Absorption, CPP
    # vlxtag: RKS, Absorption, CPP
    # vlxtag: RHF, ECD, CPP
    # vlxtag: RKS, ECD, CPP
    """

    def _get_precond(self, orb_ene, nocc, norb, w, d):
        """
        Constructs the preconditioners.
        """

        # spawning needed components
        ediag, sdiag = self.construct_ediag_sdiag_half(orb_ene, nocc, norb)
        ediag_sq = ediag**2
        sdiag_sq = sdiag**2
        sdiag_fp = sdiag**4
        w_sq = w**2
        d_sq = d**2

        # constructing matrix block diagonals
        a_diag = ediag * (ediag_sq - (w_sq - d_sq) * sdiag_sq)
        b_diag = (w * sdiag) * (ediag_sq - (w_sq + d_sq) * sdiag_sq)
        c_diag = (d * sdiag) * (ediag_sq + (w_sq + d_sq) * sdiag_sq)
        d_diag = (2 * w * d * ediag) * sdiag_sq
        p_diag = 1.0 / ((ediag_sq - (w_sq - d_sq) * sdiag_sq)**2 +
                        (4 * w_sq * d_sq * sdiag_fp))

        pa_diag = p_diag * a_diag
        pb_diag = p_diag * b_diag
        pc_diag = p_diag * c_diag
        pd_diag = p_diag * d_diag

        p_mat = np.hstack((
            pa_diag.reshape(-1, 1),
            pb_diag.reshape(-1, 1),
            pc_diag.reshape(-1, 1),
            pd_diag.reshape(-1, 1),
        ))

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
        molecule_sanity_check(molecule, 'restricted')
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
            orb_ene = scf_results['E_alpha']
        else:
            orb_ene = None
        orb_ene = self.comm.bcast(orb_ene, root=mpi_master())
        norb = orb_ene.shape[0]
        nocc = molecule.number_of_alpha_occupied_orbitals(basis)

        # ERI information
        eri_dict = self._init_eri(molecule, basis)
        # DFT information
        dft_dict = self._init_dft(molecule, scf_results)
        # PE information
        pe_dict = self._init_pe(molecule, basis)
        # CPCM information
        self._init_cpcm(molecule, basis)

        # right-hand side (gradient)
        if self.rank == mpi_master():
            self.nonlinear = (v_grad is not None)
        self.nonlinear = self.comm.bcast(self.nonlinear, root=mpi_master())

        if not self.nonlinear:
            b_grad = self.get_complex_prop_grad(self.b_operator,
                                                self.b_components, molecule,
                                                basis, scf_results)
            if self.rank == mpi_master():
                v_grad = {
                    (op, w): v for op, v in zip(self.b_components, b_grad)
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
            w: self._get_precond(orb_ene, nocc, norb, w, d)
            for w in self.frequencies
        }

        # distribute the gradient and right-hand side:
        # dist_grad will be used for calculating the subspace matrix
        # equation and residuals, dist_rhs for the initial guess
        dist_grad = {}
        dist_rhs = {}
        for key in op_freq_keys:
            if self.rank == mpi_master():
                gradger, gradung = self._decomp_grad(v_grad[key])
                grad_mat = np.hstack((
                    gradger.real.reshape(-1, 1),
                    gradung.real.reshape(-1, 1),
                    gradung.imag.reshape(-1, 1),
                    gradger.imag.reshape(-1, 1),
                ))
                rhs_mat = np.hstack((
                    gradger.real.reshape(-1, 1),
                    gradung.real.reshape(-1, 1),
                    -gradung.imag.reshape(-1, 1),
                    -gradger.imag.reshape(-1, 1),
                ))
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
                            (op, w): v
                            for op, v in zip(self.b_components, b_grad)
                            for w in extra_freqs
                        }
                        extra_op_freq_keys = list(extra_v_grad.keys())
                    else:
                        extra_op_freq_keys = None
                    extra_op_freq_keys = self.comm.bcast(extra_op_freq_keys,
                                                         root=mpi_master())

                    extra_precond = {
                        w: self._get_precond(orb_ene, nocc, norb, w, d)
                        for w in extra_freqs
                    }

                    extra_dist_rhs = {}
                    for key in extra_op_freq_keys:
                        if self.rank == mpi_master():
                            gradger, gradung = self._decomp_grad(
                                extra_v_grad[key])
                            grad_mat = np.hstack((
                                gradger.real.reshape(-1, 1),
                                gradung.real.reshape(-1, 1),
                                gradung.imag.reshape(-1, 1),
                                gradger.imag.reshape(-1, 1),
                            ))
                            rhs_mat = np.hstack((
                                gradger.real.reshape(-1, 1),
                                gradung.real.reshape(-1, 1),
                                -gradung.imag.reshape(-1, 1),
                                -gradger.imag.reshape(-1, 1),
                            ))
                        else:
                            grad_mat = None
                            rhs_mat = None
                        dist_grad[key] = DistributedArray(grad_mat, self.comm)
                        extra_dist_rhs[key] = DistributedArray(
                            rhs_mat, self.comm)

                    bger, bung = self._setup_trials(extra_dist_rhs,
                                                    extra_precond)

                    profiler.set_timing_key('Preparation')

                    self._e2n_half_size(bger, bung, molecule, basis,
                                        scf_results, eri_dict, dft_dict,
                                        pe_dict, profiler)

        # generate initial guess from scratch
        else:
            bger, bung = self._setup_trials(dist_rhs, precond)

            profiler.set_timing_key('Preparation')

            self._e2n_half_size(bger, bung, molecule, basis, scf_results,
                                eri_dict, dft_dict, pe_dict, profiler)

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

            e2gg = self._dist_bger.matmul_AtB(self._dist_e2bger, 2.0)
            e2uu = self._dist_bung.matmul_AtB(self._dist_e2bung, 2.0)
            s2ug = self._dist_bung.matmul_AtB(self._dist_bger, 2.0)

            for op, w in op_freq_keys:
                if (iteration == 0 or
                        relative_residual_norm[(op, w)] > self.conv_thresh):
                    solve_data = self._solve_reduced_problem(
                        op, w, n_ger, n_ung, e2gg, e2uu, s2ug, dist_grad,
                        v_grad)
                    x = solve_data['solution']
                    r = solve_data['residual']
                    current_solutions[(op, w)] = x

                    if self.nonlinear:
                        focks[(op, w)] = solve_data['fock']

                    if self.rank == mpi_master():
                        xvs.append((op, w, solve_data['property']))

                    r_norms_2 = 2.0 * r.squared_norm(axis=0)
                    x_norms_2 = 2.0 * x.squared_norm(axis=0)

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
            self._e2n_half_size(new_trials_ger, new_trials_ung, molecule, basis,
                                scf_results, eri_dict, dft_dict, pe_dict,
                                profiler)

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
            a_grad = self.get_complex_prop_grad(self.a_operator,
                                                self.a_components, molecule,
                                                basis, scf_results)

            if self.is_converged:
                if self.rank == mpi_master():
                    va = {op: v for op, v in zip(self.a_components, a_grad)}
                    rsp_funcs = {}

                    # final h5 file for response solutions
                    if self.filename is not None:
                        final_h5_fname = f'{self.filename}.h5'
                    else:
                        final_h5_fname = None

                for bop, w in solutions:
                    x = self.get_full_solution_vector(solutions[(bop, w)])

                    if self.rank == mpi_master():
                        for aop in self.a_components:
                            rsp_funcs[(aop, bop, w)] = -np.dot(va[aop], x)

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

                    ret_dict = self._attach_molecular_metadata(
                        ret_dict, molecule)

                    self._print_results(ret_dict)

                    # write spectrum to h5 file
                    if final_h5_fname is not None:
                        self.write_cpp_rsp_results_to_hdf5(
                            final_h5_fname, ret_dict)

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
                               dist_grad, v_grad):
        """
        Solves one reduced CPP linear system and forms the corresponding
        solution, residual, and projected property.
        """

        grad_rg = dist_grad[(op, w)].get_column(0)
        grad_ru = dist_grad[(op, w)].get_column(1)
        grad_iu = dist_grad[(op, w)].get_column(2)
        grad_ig = dist_grad[(op, w)].get_column(3)

        g_realger = self._dist_bger.matmul_AtB(grad_rg, 2.0)
        g_imagger = self._dist_bger.matmul_AtB(grad_ig, 2.0)
        g_realung = self._dist_bung.matmul_AtB(grad_ru, 2.0)
        g_imagung = self._dist_bung.matmul_AtB(grad_iu, 2.0)

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
            xv = np.dot(x_full, v_grad[(op, w)])
        else:
            xv = None

        return {
            'coeffs': c,
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
        self._e2n_half_size(new_bger, new_bung, molecule, basis, scf_results,
                            eri_dict, dft_dict, pe_dict, profiler)

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
        a_prop_grad = self.get_complex_prop_grad(self.a_operator,
                                                 self.a_components, molecule,
                                                 basis, scf_results)

        if self.rank == mpi_master():
            nocc = molecule.number_of_alpha_occupied_orbitals(basis)
            norb = scf_results['E_alpha'].shape[0]
            nvir = norb - nocc
            n_ov = nocc * nvir

            mo_occ = scf_results['C_alpha'][:, :nocc].copy()
            mo_vir = scf_results['C_alpha'][:, nocc:].copy()

            if self.property == 'absorption':
                # vector representation of absorption cross-section
                vec = (a_prop_grad[0] * cpp_solution_vector_x +
                       a_prop_grad[1] * cpp_solution_vector_y +
                       a_prop_grad[2] * cpp_solution_vector_z).imag / 3.0
                vec *= 4.0 * np.pi * w * fine_structure_constant()
            elif self.property == 'ecd':
                # vector representation of Delta epsilon
                vec = (a_prop_grad[0] * cpp_solution_vector_x / w +
                       a_prop_grad[1] * cpp_solution_vector_y / w +
                       a_prop_grad[2] * cpp_solution_vector_z / w).imag
                vec /= (3.0 * w)
                vec *= w**2 * extinction_coefficient_from_beta()

            # excitation and de-excitation
            z_mat_ov = vec[:n_ov].reshape(nocc, nvir)
            y_mat_ov = vec[n_ov:].reshape(nocc, nvir)

            prop_diag_D = np.sum(z_mat_ov, axis=1) + np.sum(y_mat_ov, axis=1)
            prop_diag_A = np.sum(z_mat_ov, axis=0) + np.sum(y_mat_ov, axis=0)

            prop_dens_D = -np.linalg.multi_dot(
                [mo_occ, np.diag(prop_diag_D), mo_occ.T])
            prop_dens_A = np.linalg.multi_dot(
                [mo_vir, np.diag(prop_diag_A), mo_vir.T])

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
