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

from .veloxchemlib import mpi_master
from .profiler import Profiler
from .distributedarray import DistributedArray
from .lrsolverbase import LinearResponseSolverBase
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           ri_sanity_check, dft_sanity_check, pe_sanity_check,
                           solvation_model_sanity_check)
from .errorhandler import assert_msg_critical
from .mathutils import safe_solve
from .checkpoint import (check_rsp_hdf5, write_rsp_solution_with_multiple_keys)


class LinearResponseSolver(LinearResponseSolverBase):
    """
    Implements linear response solver.

    # vlxtag: RHF, Polarizability, LR
    # vlxtag: RKS, Polarizability, LR

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - a_operator: The A operator.
        - a_components: Cartesian components of the A operator.
        - b_operator: The B operator.
        - b_components: Cartesian components of the B operator.
        - frequencies: The frequencies.
    """

    def compute(self, molecule, basis, scf_results, v_grad=None):
        """
        Performs linear response calculation for a molecule and a basis set.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_results:
            The dictionary of tensors from converged SCF wavefunction.
        :param v_grad:
            The gradients on the right-hand side. If not provided, v_grad will
            be computed for the B operator.

        :return:
            A dictionary containing response functions, solutions and a
            dictionary containing solutions and kappa values when called from
            a non-linear response module.
        """

        if self.norm_thresh is None:
            self.norm_thresh = self.conv_thresh * 1.0e-6
        if self.lindep_thresh is None:
            self.lindep_thresh = self.conv_thresh * 1.0e-2

        self._dist_bger = None
        self._dist_bung = None
        self._dist_e2bger = None
        self._dist_e2bung = None

        self.has_external_rhs = False

        # make sure that lr_property is properly set
        if self.property is not None:
            self.set_lr_property(self.property)

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

        # check solvation model setup
        if self.rank == mpi_master():
            assert_msg_critical(
                'solvation_model' not in scf_results,
                type(self).__name__ + ': Solvation model not implemented')

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
            self._print_header('Linear Response Solver',
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
        self._init_cpcm(molecule)

        # right-hand side (gradient)
        if self.rank == mpi_master():
            self.has_external_rhs = (v_grad is not None)
        self.has_external_rhs = self.comm.bcast(self.has_external_rhs,
                                                root=mpi_master())

        if not self.has_external_rhs:
            b_grad = self.get_prop_grad(self.b_operator, self.b_components,
                                        molecule, basis, scf_results)
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

        self.frequencies = []
        for (op, w) in op_freq_keys:
            if w not in self.frequencies:
                self.frequencies.append(w)

        precond = {
            w: self._get_precond(orb_ene, nocc, norb, w)
            for w in self.frequencies
        }

        # distribute the right-hand side
        # dist_grad will also serve as initial guess

        dist_grad = {}
        for key in op_freq_keys:
            if self.rank == mpi_master():
                gradger, gradung = self._decomp_grad(v_grad[key])
                grad_mat = np.hstack((
                    gradger.reshape(-1, 1),
                    gradung.reshape(-1, 1),
                ))
            else:
                grad_mat = None
            dist_grad[key] = DistributedArray(grad_mat, self.comm)

        rsp_vector_labels = [
            'LR_bger_half_size',
            'LR_bung_half_size',
            'LR_e2bger_half_size',
            'LR_e2bung_half_size',
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
            # `self.nonlinear` is inactive (and also `self.has_external_rhs` for
            # LR solver)

            # print warning if frequencies is not present in the restart file
            if checkpoint_frequencies is None and not (self.nonlinear or
                                                       self.has_external_rhs):
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
            elif not (self.nonlinear or self.has_external_rhs):
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
                        w: self._get_precond(orb_ene, nocc, norb, w)
                        for w in extra_freqs
                    }

                    extra_dist_grad = {}
                    for key in extra_op_freq_keys:
                        if self.rank == mpi_master():
                            gradger, gradung = self._decomp_grad(
                                extra_v_grad[key])
                            grad_mat = np.hstack((
                                gradger.reshape(-1, 1),
                                gradung.reshape(-1, 1),
                            ))
                        else:
                            grad_mat = None
                        extra_dist_grad[key] = DistributedArray(
                            grad_mat, self.comm)

                    dist_grad.update(extra_dist_grad)

                    bger, bung = self._setup_trials(extra_dist_grad,
                                                    extra_precond)

                    profiler.set_timing_key('Preparation')

                    self._e2n_half_size(bger, bung, molecule, basis,
                                        scf_results, eri_dict, dft_dict,
                                        pe_dict, profiler)

        # generate initial guess from scratch
        else:
            bger, bung = self._setup_trials(dist_grad, precond)

            profiler.set_timing_key('Preparation')

            self._e2n_half_size(bger, bung, molecule, basis, scf_results,
                                eri_dict, dft_dict, pe_dict, profiler)

        profiler.check_memory_usage('Initial guess')

        solutions = {}
        residuals = {}
        relative_residual_norm = {}

        iter_per_trial_in_hours = None

        # start iterations
        for iteration in range(self.max_iter):

            iter_start_time = tm.time()

            profiler.set_timing_key(f'Iteration {iteration + 1}')

            profiler.start_timer('ReducedSpace')

            n_ger = self._dist_bger.shape(1)
            n_ung = self._dist_bung.shape(1)

            e2gg = self._dist_bger.matmul_AtB(self._dist_e2bger, 2.0)
            e2uu = self._dist_bung.matmul_AtB(self._dist_e2bung, 2.0)
            s2ug = self._dist_bung.matmul_AtB(self._dist_bger, 2.0)

            xvs = []
            self._cur_iter = iteration

            for op, freq in op_freq_keys:
                if (iteration > 0 and
                        relative_residual_norm[(op, freq)] < self.conv_thresh):
                    continue

                gradger = dist_grad[(op, freq)].get_column(0)
                gradung = dist_grad[(op, freq)].get_column(1)

                g_ger = self._dist_bger.matmul_AtB(gradger, 2.0)
                g_ung = self._dist_bung.matmul_AtB(gradung, 2.0)

                if self.rank == mpi_master():
                    mat = np.zeros((n_ger + n_ung, n_ger + n_ung))
                    mat[:n_ger, :n_ger] = e2gg[:, :]
                    mat[:n_ger, n_ger:] = -freq * s2ug.T[:, :]
                    mat[n_ger:, :n_ger] = -freq * s2ug[:, :]
                    mat[n_ger:, n_ger:] = e2uu[:, :]

                    g = np.zeros(n_ger + n_ung)
                    g[:n_ger] = g_ger[:]
                    g[n_ger:] = g_ung[:]

                    c = safe_solve(mat, g)
                else:
                    c = None
                c = self.comm.bcast(c, root=mpi_master())

                c_ger = c[:n_ger]
                c_ung = c[n_ger:]

                x_ger = self._dist_bger.matmul_AB_no_gather(c_ger)
                x_ung = self._dist_bung.matmul_AB_no_gather(c_ung)

                e2x_ger = self._dist_e2bger.matmul_AB_no_gather(c_ger)
                e2x_ung = self._dist_e2bung.matmul_AB_no_gather(c_ung)

                s2x_ger = x_ger.data
                s2x_ung = x_ung.data

                r_ger = e2x_ger.data - freq * s2x_ung - gradger.data
                r_ung = e2x_ung.data - freq * s2x_ger - gradung.data

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

                x_full = self.get_full_solution_vector(x)

                if self.rank == mpi_master():
                    xv = np.dot(x_full, v_grad[(op, freq)])
                    xvs.append((op, freq, xv))

                r_norms_2 = 2.0 * r.squared_norm(axis=0)
                x_norms_2 = 2.0 * x.squared_norm(axis=0)

                rn = np.sqrt(np.sum(r_norms_2))
                xn = np.sqrt(np.sum(x_norms_2))

                if xn != 0:
                    relative_residual_norm[(op, freq)] = 2.0 * rn / xn
                else:
                    relative_residual_norm[(op, freq)] = 2.0 * rn

                if relative_residual_norm[(op, freq)] < self.conv_thresh:
                    solutions[(op, freq)] = x
                else:
                    residuals[(op, freq)] = r

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

            profiler.start_timer('Orthonorm.')

            # update trial vectors
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
            self._print_convergence('Linear response')

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of LR solver')
        profiler.print_memory_usage(self.ostream)

        self._dist_bger = None
        self._dist_bung = None
        self._dist_e2bger = None
        self._dist_e2bung = None

        # calculate response functions
        if not self.has_external_rhs:
            a_grad = self.get_prop_grad(self.a_operator, self.a_components,
                                        molecule, basis, scf_results)

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

                    self._print_results(rsp_funcs, self.ostream)

                    return {
                        'a_operator': self.a_operator,
                        'a_components': self.a_components,
                        'b_operator': self.b_operator,
                        'b_components': self.b_components,
                        'response_functions': rsp_funcs,
                        'solutions': solutions,
                    }
                else:
                    # non-master rank
                    return {'solutions': solutions}
            else:
                # not converged
                return {}

        else:
            # has_external_rhs
            if self.is_converged:
                return {'solutions': solutions}
            else:
                return {}

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

        # spawning needed components

        ediag, sdiag = self.construct_ediag_sdiag_half(orb_ene, nocc, norb)

        ediag_sq = ediag**2
        sdiag_sq = sdiag**2
        w_sq = w**2

        # constructing matrix block diagonals

        pa_diag = ediag / (ediag_sq - w_sq * sdiag_sq)
        pb_diag = (w * sdiag) / (ediag_sq - w_sq * sdiag_sq)

        p_mat = np.hstack((
            pa_diag.reshape(-1, 1),
            pb_diag.reshape(-1, 1),
        ))

        return DistributedArray(p_mat, self.comm)
