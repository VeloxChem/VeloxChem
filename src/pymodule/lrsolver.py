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
import numpy as np
import time as tm
import sys

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .profiler import Profiler
from .distributedarray import DistributedArray
from .linearsolver import LinearSolver
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           dft_sanity_check, pe_sanity_check)
from .errorhandler import assert_msg_critical
from .checkpoint import (check_rsp_hdf5, create_hdf5,
                         write_rsp_solution_with_multiple_keys)


class LinearResponseSolver(LinearSolver):
    """
    Implements linear response solver.

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

    def __init__(self, comm=None, ostream=None):
        """
        Initializes linear response solver to default setup.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        super().__init__(comm, ostream)

        # operators and frequencies
        self.a_operator = 'electric dipole'
        self.a_components = 'xyz'
        self.b_operator = 'electric dipole'
        self.b_components = 'xyz'
        self.frequencies = (0,)

        self._input_keywords['response'].update({
            'a_operator': ('str_lower', 'A operator'),
            'a_components': ('str_lower', 'Cartesian components of A operator'),
            'b_operator': ('str_lower', 'B operator'),
            'b_components': ('str_lower', 'Cartesian components of B operator'),
            'frequencies': ('seq_range', 'frequencies'),
        })

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in linear response solver.

        :param rsp_dict:
            The dictionary of response input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}

        super().update_settings(rsp_dict, method_dict)

    def compute(self, molecule, basis, scf_tensors, v_grad=None):
        """
        Performs linear response calculation for a molecule and a basis set.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param v_grad:
            The gradients on the right-hand side. If not provided, v_grad will
            be computed for the B operator.

        :return:
            A dictionary containing response functions, solutions and a
            dictionarry containing solutions and kappa values when called from
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

        # check molecule
        molecule_sanity_check(molecule)

        # check SCF results
        scf_results_sanity_check(self, scf_tensors)

        # check dft setup
        dft_sanity_check(self, 'compute')

        # check pe setup
        pe_sanity_check(self)

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
            self._print_header('Linear Response Solver',
                               n_freqs=len(self.frequencies))

        self.start_time = tm.time()

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'LinearResponseSolver: not implemented for unrestricted case')

        if self.rank == mpi_master():
            orb_ene = scf_tensors['E_alpha']
        else:
            orb_ene = None
        orb_ene = self.comm.bcast(orb_ene, root=mpi_master())
        norb = orb_ene.shape[0]
        nocc = molecule.number_of_alpha_electrons()

        # ERI information
        eri_dict = self._init_eri(molecule, basis)

        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors)

        # PE information
        pe_dict = self._init_pe(molecule, basis)

        # right-hand side (gradient)
        if self.rank == mpi_master():
            self.has_external_rhs = (v_grad is not None)
        self.has_external_rhs = self.comm.bcast(self.has_external_rhs,
                                                root=mpi_master())

        if not self.has_external_rhs:
            b_grad = self.get_prop_grad(self.b_operator, self.b_components,
                                        molecule, basis, scf_tensors)
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

        # generate initial guess from scratch
        else:
            bger, bung = self._setup_trials(dist_grad, precond)

            self._e2n_half_size(bger, bung, molecule, basis, scf_tensors,
                                eri_dict, dft_dict, pe_dict)

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

                    c = np.linalg.solve(mat, g)
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

            self._e2n_half_size(new_trials_ger, new_trials_ung, molecule, basis,
                                scf_tensors, eri_dict, dft_dict, pe_dict,
                                profiler)

            iter_in_hours = (tm.time() - iter_start_time) / 3600
            iter_per_trial_in_hours = iter_in_hours / n_new_trials

            profiler.check_memory_usage(
                'Iteration {:d} sigma build'.format(iteration + 1))

        self._write_checkpoint(molecule, basis, dft_dict, pe_dict,
                               rsp_vector_labels)

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
                                        molecule, basis, scf_tensors)

            if self.is_converged:
                if self.rank == mpi_master():
                    va = {op: v for op, v in zip(self.a_components, a_grad)}
                    rsp_funcs = {}

                    # create h5 file for response solutions
                    if (self.save_solutions and
                            self.checkpoint_file is not None):
                        final_h5_fname = str(
                            Path(self.checkpoint_file).with_suffix(
                                '.solutions.h5'))
                        create_hdf5(final_h5_fname, molecule, basis,
                                    dft_dict['dft_func_label'],
                                    pe_dict['potfile_text'])

                for bop, w in solutions:
                    x = self.get_full_solution_vector(solutions[(bop, w)])

                    if self.rank == mpi_master():
                        for aop in self.a_components:
                            rsp_funcs[(aop, bop, w)] = -np.dot(va[aop], x)

                        # write to h5 file for response solutions
                        if (self.save_solutions and
                                self.checkpoint_file is not None):
                            solution_keys = [
                                '{:s}_{:s}_{:.8f}'.format(aop, bop, w)
                                for aop in self.a_components
                            ]
                            write_rsp_solution_with_multiple_keys(
                                final_h5_fname, solution_keys, x)

                if self.rank == mpi_master():
                    # print information about h5 file for response solutions
                    if (self.save_solutions and
                            self.checkpoint_file is not None):
                        checkpoint_text = 'Response solution vectors written to file: '
                        checkpoint_text += final_h5_fname
                        self.ostream.print_info(checkpoint_text)
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
                    return {'solutions': solutions}

        else:
            if self.is_converged:
                return {'solutions': solutions}

        return None

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

    def _print_iteration(self, relative_residual_norm, xvs):
        """
        Prints information of the iteration.

        :param relative_residual_norm:
            Relative residual norms.
        :param xvs:
            A list of tuples containing operator component, frequency, and
            property.
        """

        width = 92
        output_header = '*** Iteration:   {} '.format(self._cur_iter + 1)
        output_header += '* Residuals (Max,Min): '
        output_header += '{:.2e} and {:.2e}'.format(
            max(relative_residual_norm.values()),
            min(relative_residual_norm.values()))
        self.ostream.print_header(output_header.ljust(width))
        self.ostream.print_blank()

        if self.print_level > 1:
            for op, freq, xv in xvs:
                ops_label = '<<{};{}>>_{:.4f}'.format(op, op, freq)
                rel_res = relative_residual_norm[(op, freq)]
                output_iter = '{:<15s}: {:15.8f} '.format(ops_label, -xv)
                output_iter += 'Residual Norm: {:.8f}'.format(rel_res)
                self.ostream.print_header(output_iter.ljust(width))
            self.ostream.print_blank()
            self.ostream.flush()

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

        for (op, w), vec in vectors.items():
            v = self._preconditioning(precond[w], vec)
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

    def _print_results(self, rsp_funcs, ostream):
        """
        Prints polarizability to output stream.

        :param rsp_funcs:
            The response functions.
        :param ostream:
            The output stream.
        """

        width = 92

        dipole_ops = ['dipole', 'electric dipole', 'electric_dipole']

        if self.a_operator in dipole_ops and self.b_operator in dipole_ops:

            for w in self.frequencies:
                w_str = 'Polarizability (w={:.4f})'.format(w)
                ostream.print_header(w_str.ljust(width))
                ostream.print_header(('-' * len(w_str)).ljust(width))

                valstr = '{:<5s}'.format('')
                for b in self.b_components:
                    valstr += '{:>15s}'.format(b.upper())
                ostream.print_header(valstr.ljust(width))

                for a in self.a_components:
                    valstr = '{:<5s}'.format(a.upper())
                    for b in self.b_components:
                        prop = -rsp_funcs[(a, b, w)]
                        valstr += '{:15.8f}'.format(prop)
                    ostream.print_header(valstr.ljust(width))
                ostream.print_blank()

        ostream.flush()
