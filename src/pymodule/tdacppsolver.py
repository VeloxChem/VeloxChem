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
import numpy as np
import time as tm
import math
import sys

from .veloxchemlib import (mpi_master, hartree_in_wavenumber, hartree_in_ev,
                           hartree_in_inverse_nm, fine_structure_constant,
                           extinction_coefficient_from_beta)
from .outputstream import OutputStream
from .profiler import Profiler
from .distributedarray import DistributedArray
from .linearsolver import LinearSolver
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           dft_sanity_check, pe_sanity_check,
                           solvation_model_sanity_check)
from .errorhandler import assert_msg_critical
from .checkpoint import (check_rsp_hdf5, write_rsp_hdf5,
                         write_rsp_solution_with_multiple_keys)
from .inputparser import parse_seq_fixed


class ComplexResponseTDA(LinearSolver):
    """
    Implements the complex linear response solver using the Tamm-Dancoff approximation.

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
        - damping: The damping parameter.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes complex linear response solver to default setup.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        super().__init__(comm, ostream)

        self.a_operator = 'electric dipole'
        self.a_components = 'xyz'
        self.b_operator = 'electric dipole'
        self.b_components = 'xyz'

        self.cpp_flag = None

        self.frequencies = (0,)
        self.damping = 1000.0 / hartree_in_wavenumber()

        self._input_keywords['response'].update({
            'a_operator': ('str_lower', 'A operator'),
            'a_components': ('str_lower', 'Cartesian components of A operator'),
            'b_operator': ('str_lower', 'B operator'),
            'b_components': ('str_lower', 'Cartesian components of B operator'),
            'frequencies': ('seq_range', 'frequencies'),
            'damping': ('float', 'damping parameter'),
        })

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in complex liner response solver.

        :param rsp_dict:
            The dictionary of response input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}

        super().update_settings(rsp_dict, method_dict)

    def set_cpp_flag(self, flag):
        """
        Sets CPP flag (absorption or ecd).

        :param flag:
            The flag (absorption or ecd).
        """

        assert_msg_critical(flag.lower() in ['absorption', 'ecd'],
                            'ComplexResponse: invalid CPP flag')

        self.cpp_flag = flag.lower()

        if self.cpp_flag == 'absorption':
            self.a_operator = 'electric dipole'
            self.a_components = 'xyz'
            self.b_operator = 'electric dipole'
            self.b_components = 'xyz'

        elif self.cpp_flag == 'ecd':
            self.a_operator = 'magnetic dipole'
            self.a_components = 'xyz'
            self.b_operator = 'linear momentum'
            self.b_components = 'xyz'

    def _get_precond(self, orb_ene, nocc, norb, w, d):
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
        :param d:
            The damping parameter.

        :return:
            The distributed preconditioners.
        """

        # spawning needed components

        ediag, sdiag_unused = self.construct_ediag_sdiag_half(
            orb_ene, nocc, norb)

        # constructing matrix block diagonals

        A0mg = ediag - w
        invdiag = 1.0 / (A0mg * A0mg + d * d)

        pa_diag = invdiag * A0mg
        pb_diag = invdiag * d

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
        v_in_ig = v_in.data[:, 1]

        v_out_rg = (pa * v_in_rg + pb * v_in_ig)
        v_out_ig = (pb * v_in_rg - pa * v_in_ig)

        v_mat = np.hstack((
            v_out_rg.reshape(-1, 1),
            v_out_ig.reshape(-1, 1),
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

        for (op, w), vec in vectors.items():
            v = self._preconditioning(precond[w], vec)
            norms_2 = v.squared_norm(axis=0)
            vn = np.sqrt(np.sum(norms_2))

            if vn > self.norm_thresh:
                norms = np.sqrt(norms_2)
                # real gerade
                if norms[0] > self.norm_thresh:
                    trials_ger.append(v.data[:, 0])
                # imaginary gerade
                if norms[1] > self.norm_thresh:
                    trials_ger.append(v.data[:, 1])

        new_ger = np.array(trials_ger).T

        dist_new_ger = DistributedArray(new_ger, self.comm, distribute=False)

        return dist_new_ger

    def compute(self, molecule, basis, scf_tensors, v_grad=None):
        """
        Solves for the response vector iteratively while checking the residuals
        for convergence.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param v_grad:
            The gradients on the right-hand side. If not provided, v_grad will
            be computed for the B operator.

        :return:
            A dictionary containing response functions, solutions and a
            dictionary containing solutions and kappa values when called from
            a non-linear response module.
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
                'ComplexResponse: Undefined or invalid a_component')
        for comp in self.b_components:
            assert_msg_critical(
                self.is_valid_component(comp, self.b_operator),
                'ComplexResponse: Undefined or invalid b_component')

        if self.norm_thresh is None:
            self.norm_thresh = self.conv_thresh * 1.0e-6
        if self.lindep_thresh is None:
            self.lindep_thresh = self.conv_thresh * 1.0e-2

        self._dist_bger = None
        self._dist_e2bger = None

        self.nonlinear = False
        self._dist_fock_ger = None

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
            self._print_header('TDA Complex Response Solver',
                               n_freqs=len(self.frequencies))

        self.start_time = tm.time()

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'ComplexResponse: not implemented for unrestricted case')

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

        # CPCM information
        self._init_cpcm(molecule)

        # right-hand side (gradient)
        if self.rank == mpi_master():
            self.nonlinear = (v_grad is not None)
        self.nonlinear = self.comm.bcast(self.nonlinear, root=mpi_master())

        if not self.nonlinear:
            b_grad = self.get_complex_prop_grad(self.b_operator,
                                                self.b_components, molecule,
                                                basis, scf_tensors)
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
                # No decomposing into ger ung
                gradger = v_grad[key]
                grad_mat = np.hstack((
                    gradger.real.reshape(-1, 1),
                    gradger.imag.reshape(-1, 1),
                ))
                rhs_mat = np.hstack((
                    gradger.real.reshape(-1, 1),
                    -gradger.imag.reshape(-1, 1),
                ))
            else:
                grad_mat = None
                rhs_mat = None

            dist_grad[key] = DistributedArray(grad_mat, self.comm)
            dist_rhs[key] = DistributedArray(rhs_mat, self.comm)

        if self.nonlinear:
            rsp_vector_labels = ['CLR_bger', 'CLR_e2bger', 'CLR_Fock_ger']
        else:
            rsp_vector_labels = ['CLR_bger', 'CLR_e2bger']

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
            bger = self.setup_trials(dist_rhs, precond)

            profiler.set_timing_key('Preparation')

            # We set up sigma in another way since we do not have ger ung

            bger_loc = bger.get_full_matrix(root=mpi_master())

            tdens = self._get_trans_densities(bger_loc, scf_tensors, molecule)

            fock = self._comp_lr_fock(tdens, molecule, basis, eri_dict,
                                      dft_dict, pe_dict)

            if self.rank == mpi_master():
                sig_mat = self._get_sigmas(fock, scf_tensors, molecule,
                                           bger_loc)
            else:
                sig_mat = None

            sig_mat = DistributedArray(sig_mat, self.comm)

            self._append_trial_sigma_vectors(bger, sig_mat)

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
            self._cur_iter = iteration

            n_ger = self._dist_bger.shape(1)

            e2gg = self._dist_bger.matmul_AtB(self._dist_e2bger)
            s2gg = self._dist_bger.matmul_AtB(self._dist_bger)

            for op, w in op_freq_keys:
                if (iteration == 0 or
                        relative_residual_norm[(op, w)] > self.conv_thresh):

                    grad_rg = dist_grad[(op, w)].get_column(0)
                    grad_ig = dist_grad[(op, w)].get_column(1)

                    g_realger = self._dist_bger.matmul_AtB(grad_rg)
                    g_imagger = self._dist_bger.matmul_AtB(grad_ig)

                    # creating gradient and matrix for linear equation

                    size = 2 * n_ger

                    if self.rank == mpi_master():

                        g = np.zeros(size)

                        g[:n_ger] = g_realger[:]
                        g[size - n_ger:] = -g_imagger[:]

                        mat = np.zeros((size, size))

                        mat[:n_ger, :n_ger] += e2gg[:, :]
                        mat[size - n_ger:, size - n_ger:] += -e2gg[:, :]

                        mat[:n_ger, :n_ger] += -w * s2gg[:, :]

                        mat[:n_ger, size - n_ger:] += d * s2gg[:, :]

                        mat[size - n_ger:, :n_ger] += d * s2gg[:, :]

                        mat[size - n_ger:, size - n_ger:] += w * s2gg[:, :]

                        # solving matrix equation

                        c = np.linalg.solve(mat, g)

                    else:
                        c = None
                    c = self.comm.bcast(c, root=mpi_master())

                    # extracting the 4 components of c...

                    c_realger = c[:n_ger]
                    c_imagger = c[size - n_ger:]

                    # ...and projecting them onto respective subspace

                    x_realger = self._dist_bger.matmul_AB_no_gather(c_realger)
                    x_imagger = self._dist_bger.matmul_AB_no_gather(c_imagger)

                    # composing E2 matrices projected onto solution subspace

                    e2realger = self._dist_e2bger.matmul_AB_no_gather(c_realger)
                    e2imagger = self._dist_e2bger.matmul_AB_no_gather(c_imagger)

                    if self.nonlinear:
                        fock_realger = self._dist_fock_ger.matmul_AB_no_gather(
                            c_realger)
                        fock_imagger = self._dist_fock_ger.matmul_AB_no_gather(
                            c_imagger)

                        fock_full_data = (fock_realger.data -
                                          1j * fock_imagger.data)

                        focks[(op, w)] = DistributedArray(fock_full_data,
                                                          self.comm,
                                                          distribute=False)

                    # calculating the residual components

                    s2realger = x_realger.data
                    s2imagger = x_imagger.data

                    r_realger = (e2realger.data - w * s2realger +
                                 d * s2imagger - grad_rg.data)
                    r_imagger = (-e2imagger.data + w * s2imagger +
                                 d * s2realger + grad_ig.data)

                    r_data = np.hstack((
                        r_realger.reshape(-1, 1),
                        r_imagger.reshape(-1, 1),
                    ))

                    r = DistributedArray(r_data, self.comm, distribute=False)

                    # calculating relative residual norm
                    # for convergence check

                    x_data = np.hstack((
                        x_realger.data.reshape(-1, 1),
                        x_imagger.data.reshape(-1, 1),
                    ))

                    x = DistributedArray(x_data, self.comm, distribute=False)

                    x_full = self.get_full_solution_vector(x)
                    if self.rank == mpi_master():
                        xv = np.dot(x_full, v_grad[(op, w)])
                        xvs.append((op, w, xv))

                    r_norms_2 = r.squared_norm(axis=0)
                    x_norms_2 = x.squared_norm(axis=0)

                    rn = np.sqrt(np.sum(r_norms_2))
                    xn = np.sqrt(np.sum(x_norms_2))

                    if xn != 0:
                        relative_residual_norm[(op, w)] = rn / xn
                    else:
                        relative_residual_norm[(op, w)] = rn

                    if relative_residual_norm[(op, w)] < self.conv_thresh:
                        solutions[(op, w)] = x
                    else:
                        residuals[(op, w)] = r

            # write to output
            if self.rank == mpi_master():

                self.ostream.print_info(
                    '{:d} gerade trial vectors in reduced space'.format(n_ger))
                self.ostream.print_blank()

                if self.print_level > 1:
                    profiler.print_memory_subspace(
                        {
                            'dist_bger': self._dist_bger,
                            'dist_e2bger': self._dist_e2bger,
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

            # spawning new trial vectors from residuals

            new_trials_ger = self.setup_trials(residuals, precond,
                                               self._dist_bger)

            residuals.clear()

            profiler.stop_timer('Orthonorm.')

            if self.rank == mpi_master():
                n_new_trials = new_trials_ger.shape(1)
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

            # creating new sigma and rho linear transformations

            # Once again using the other way of computing sigma

            new_trials_ger_loc = new_trials_ger.get_full_matrix(
                root=mpi_master())

            tdens = self._get_trans_densities(new_trials_ger_loc, scf_tensors,
                                              molecule)

            fock = self._comp_lr_fock(tdens, molecule, basis, eri_dict,
                                      dft_dict, pe_dict)

            if self.rank == mpi_master():
                sig_mat = self._get_sigmas(fock, scf_tensors, molecule,
                                           new_trials_ger_loc)
            else:
                sig_mat = None

            sig_mat = DistributedArray(sig_mat, self.comm)

            self._append_trial_sigma_vectors(new_trials_ger, sig_mat)

            iter_in_hours = (tm.time() - iter_start_time) / 3600
            iter_per_trial_in_hours = iter_in_hours / n_new_trials

            profiler.check_memory_usage(
                'Iteration {:d} sigma build'.format(iteration + 1))

        self._write_checkpoint(molecule, basis, dft_dict, pe_dict,
                               rsp_vector_labels)

        # converged?
        if self.rank == mpi_master():
            self._print_convergence('Complex response')

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of CPP solver')
        profiler.print_memory_usage(self.ostream)

        self._dist_bger = None
        self._dist_e2bger = None

        self._dist_fock_ger = None

        # calculate response functions
        if not self.nonlinear:
            a_grad = self.get_complex_prop_grad(self.a_operator,
                                                self.a_components, molecule,
                                                basis, scf_tensors)

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
                                final_h5_fname, solution_keys, x)

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

                    return ret_dict
                else:
                    return {'solutions': solutions}

        else:
            if self.is_converged:
                return {'focks': focks, 'solutions': solutions}

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

        x_realger = solution.get_full_vector(0)
        x_imagger = solution.get_full_vector(1)

        if solution.rank == mpi_master():
            return x_realger + 1j * x_imagger
        else:
            return None

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

            mo_occ = tensors['C_alpha'][:, :nocc].copy()
            mo_vir = tensors['C_alpha'][:, nocc:].copy()

            tdens = []
            for k in range(trial_mat.shape[1]):
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

        mo_occ = tensors['C_alpha'][:, :nocc].copy()
        mo_vir = tensors['C_alpha'][:, nocc:].copy()
        orb_ene = tensors['E_alpha']

        sigma_vecs = []

        for fockind in range(len(fock)):
            # 2e contribution
            mat = fock[fockind].copy()
            mat = np.matmul(mo_occ.T, np.matmul(mat, mo_vir))
            # 1e contribution
            cjb = trial_mat[:, fockind].reshape(nocc, nvir)
            mat += np.matmul(cjb, np.diag(orb_ene[nocc:]).T)
            mat -= np.matmul(np.diag(orb_ene[:nocc]), cjb)
            sigma_vecs.append(mat.reshape(nocc * nvir, 1))

        sigma_mat = sigma_vecs[0]
        for vec in sigma_vecs[1:]:
            sigma_mat = np.hstack((sigma_mat, vec))

        return sigma_mat

    def _append_trial_sigma_vectors(self, b, e2b):
        """
        Appends distributed trial vectors and sigma vectors.

        :param bger:
            The distributed gerade trial vectors.
        :param bung:
            The distributed ungerade trial vectors.
        """

        if self._dist_bger is None:
            self._dist_bger = DistributedArray(b.data,
                                               self.comm,
                                               distribute=False)
        else:
            self._dist_bger.append(b, axis=1)

        if self._dist_e2bger is None:
            self._dist_e2bger = DistributedArray(e2b.data,
                                                 self.comm,
                                                 distribute=False)
        else:
            self._dist_e2bger.append(e2b, axis=1)

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

        if (not self.nonlinear) and (self.print_level > 1):
            output_header = 'Operator:  {} ({})'.format(self.b_operator,
                                                        self.b_components)
            self.ostream.print_header(output_header.ljust(width))
            self.ostream.print_blank()

            for op, freq, xv in xvs:
                ops_label = '<<{};{}>>_{:.4f}'.format(op, op, freq)
                rel_res = relative_residual_norm[(op, freq)]
                output_iter = '{:<15s}: {:15.8f} {:15.8f}j   '.format(
                    ops_label, -xv.real, -xv.imag)
                output_iter += 'Residual Norm: {:.8f}'.format(rel_res)
                self.ostream.print_header(output_iter.ljust(width))
            self.ostream.print_blank()

        self.ostream.flush()

    def get_spectrum(self, rsp_results, x_unit):
        """
        Gets spectrum.

        :param rsp_results:
            The dictionary containing response results.
        :param x_unit:
            The unit of x-axis.

        :return:
            A dictionary containing the spectrum.
        """

        if self.cpp_flag == 'absorption':
            return self._get_absorption_spectrum(rsp_results, x_unit)

        elif self.cpp_flag == 'ecd':
            return self._get_ecd_spectrum(rsp_results, x_unit)

        return None

    def _get_absorption_spectrum(self, rsp_results, x_unit):
        """
        Gets absorption spectrum.

        :param rsp_results:
            The dictionary containing response results.
        :param x_unit:
            The unit of x-axis.

        :return:
            A dictionary containing the absorption spectrum.
        """

        assert_msg_critical(
            x_unit.lower() in ['au', 'ev', 'nm'],
            'ComplexResponse.get_spectrum: x_unit should be au, ev or nm')

        au2ev = hartree_in_ev()
        auxnm = 1.0 / hartree_in_inverse_nm()

        spectrum = {'x_data': [], 'y_data': []}

        if x_unit.lower() == 'au':
            spectrum['x_label'] = 'Photon energy [a.u.]'
        elif x_unit.lower() == 'ev':
            spectrum['x_label'] = 'Photon energy [eV]'
        elif x_unit.lower() == 'nm':
            spectrum['x_label'] = 'Wavelength [nm]'

        spectrum['y_label'] = 'Absorption cross-section [a.u.]'

        freqs = rsp_results['frequencies']
        rsp_funcs = rsp_results['response_functions']

        for w in freqs:
            if w == 0.0:
                continue

            if x_unit.lower() == 'au':
                spectrum['x_data'].append(w)
            elif x_unit.lower() == 'ev':
                spectrum['x_data'].append(au2ev * w)
            elif x_unit.lower() == 'nm':
                spectrum['x_data'].append(auxnm / w)

            axx = -rsp_funcs[('x', 'x', w)].imag
            ayy = -rsp_funcs[('y', 'y', w)].imag
            azz = -rsp_funcs[('z', 'z', w)].imag

            alpha_bar = (axx + ayy + azz) / 3.0
            sigma = 4.0 * math.pi * w * alpha_bar * fine_structure_constant()

            spectrum['y_data'].append(sigma)

        return spectrum

    def _get_ecd_spectrum(self, rsp_results, x_unit):
        """
        Gets circular dichroism spectrum.

        :param rsp_results:
            The dictionary containing response results.
        :param x_unit:
            The unit of x-axis.

        :return:
            A dictionary containing the circular dichroism spectrum.
        """

        assert_msg_critical(
            x_unit.lower() in ['au', 'ev', 'nm'],
            'ComplexResponse.get_spectrum: x_unit should be au, ev or nm')

        au2ev = hartree_in_ev()
        auxnm = 1.0 / hartree_in_inverse_nm()

        spectrum = {'x_data': [], 'y_data': []}

        if x_unit.lower() == 'au':
            spectrum['x_label'] = 'Photon energy [a.u.]'
        elif x_unit.lower() == 'ev':
            spectrum['x_label'] = 'Photon energy [eV]'
        elif x_unit.lower() == 'nm':
            spectrum['x_label'] = 'Wavelength [nm]'

        spectrum['y_label'] = 'Molar circular dichroism '
        spectrum['y_label'] += '[L mol$^{-1}$ cm$^{-1}$]'

        freqs = rsp_results['frequencies']
        rsp_funcs = rsp_results['response_functions']

        for w in freqs:
            if w == 0.0:
                continue

            if x_unit.lower() == 'au':
                spectrum['x_data'].append(w)
            elif x_unit.lower() == 'ev':
                spectrum['x_data'].append(au2ev * w)
            elif x_unit.lower() == 'nm':
                spectrum['x_data'].append(auxnm / w)

            Gxx = -rsp_funcs[('x', 'x', w)].imag / w
            Gyy = -rsp_funcs[('y', 'y', w)].imag / w
            Gzz = -rsp_funcs[('z', 'z', w)].imag / w

            beta = -(Gxx + Gyy + Gzz) / (3.0 * w)
            Delta_epsilon = beta * w**2 * extinction_coefficient_from_beta()

            spectrum['y_data'].append(Delta_epsilon)

        return spectrum

    def setup_trials(self, vectors, precond, dist_b=None, renormalize=True):
        """
        Computes orthonormalized trial vectors.

        :param vectors:
            The set of vectors.
        :param precond:
            The preconditioner.
        :param dist_bger:
            The distributed gerade subspace.
        :param dist_bung:
            The distributed ungerade subspace.
        :param renormalize:
            The flag for normalization.

        :return:
            The orthonormalized gerade and ungerade trial vectors.
        """
        dist_new_b = self._precond_trials(vectors, precond)

        if dist_new_b.data.size == 0:
            dist_new_b.data = np.zeros((dist_new_b.shape(0), 0))

        if dist_b is not None:
            # t = t - (b (b.T t))
            bT_new = dist_b.matmul_AtB_allreduce(dist_new_b)
            dist_new_proj = dist_b.matmul_AB_no_gather(bT_new)
            dist_new_b.data -= dist_new_proj.data

        if renormalize:
            if dist_new_b.data.ndim > 0 and dist_new_b.shape(0) > 0:
                dist_new_b = self.remove_linear_dependence(
                    dist_new_b, self.lindep_thresh)

                dist_new_b = self.orthogonalize_gram_schmidt(dist_new_b)

                dist_new_b = self.normalize(dist_new_b)

        if self.rank == mpi_master():
            assert_msg_critical(dist_new_b.data.size > 0,
                                'LinearSolver: trial vectors are empty')

        return dist_new_b

    @staticmethod
    def lrmat2vec(mat, nocc, norb):
        """
        Converts matrices to vectors.

        :param mat:
            The matrices.
        :param nocc:
            Number of occupied orbitals.
        :param norb:
            Number of orbitals.

        :return:
            The vectors.
        """

        # changed the linearsolver function such that we only include the
        # excitation part

        # TODO: Should be changed at some point to still call linearsolver.py
        # function

        nvir = norb - nocc

        n_ov = nocc * nvir
        vec = np.zeros(n_ov, dtype=mat.dtype)
        # excitation only
        vec[:n_ov] = mat[:nocc, nocc:].reshape(n_ov)

        return vec

    def _print_results(self, rsp_results, ostream=None):
        """
        Prints response results to output stream.

        :param rsp_results:
            The dictionary containing response results.
        :param ostream:
            The output stream.
        """

        self._print_response_functions(rsp_results, ostream)

        if self.cpp_flag == 'absorption':
            self._print_absorption_results(rsp_results, ostream)

        elif self.cpp_flag == 'ecd':
            self._print_ecd_results(rsp_results, ostream)

    def _print_response_functions(self, rsp_results, ostream=None):
        """
        Prints response functions to output stream.

        :param rsp_results:
            The dictionary containing response results.
        :param ostream:
            The output stream.
        """

        if ostream is None:
            ostream = self.ostream

        width = 92

        freqs = rsp_results['frequencies']
        rsp_funcs = rsp_results['response_functions']

        title = 'Response Functions at Given Frequencies'
        ostream.print_header(title.ljust(width))
        ostream.print_header(('=' * len(title)).ljust(width))
        ostream.print_blank()

        operator_to_name = {
            'dipole': 'Dipole',
            'electric dipole': 'Dipole',
            'electric_dipole': 'Dipole',
            'quadrupole': 'Quadru',
            'electric quadrupole': 'Quadru',
            'electric_quadrupole': 'Quadru',
            'linear_momentum': 'LinMom',
            'linear momentum': 'LinMom',
            'angular_momentum': 'AngMom',
            'angular momentum': 'AngMom',
            'magnetic dipole': 'MagDip',
            'magnetic_dipole': 'MagDip',
        }
        a_name = operator_to_name[self.a_operator]
        b_name = operator_to_name[self.b_operator]

        for w in freqs:
            title = '{:<7s} {:<7s} {:>10s} {:>15s} {:>16s}'.format(
                a_name, b_name, 'Frequency', 'Real', 'Imaginary')
            ostream.print_header(title.ljust(width))
            ostream.print_header(('-' * len(title)).ljust(width))

            for a in self.a_components:
                for b in self.b_components:
                    rsp_func_val = rsp_funcs[(a, b, w)]
                    ops_label = '<<{:>3s}  ;  {:<3s}>> {:10.4f}'.format(
                        a.lower(), b.lower(), w)
                    output = '{:<15s} {:15.8f} {:15.8f}j'.format(
                        ops_label, rsp_func_val.real, rsp_func_val.imag)
                    ostream.print_header(output.ljust(width))
            ostream.print_blank()
        ostream.flush()

    def _print_absorption_results(self, rsp_results, ostream=None):
        """
        Prints absorption results to output stream.

        :param rsp_results:
            The dictionary containing response results.
        :param ostream:
            The output stream.
        """

        if ostream is None:
            ostream = self.ostream

        width = 92

        spectrum = self.get_spectrum(rsp_results, 'au')

        title = 'Linear Absorption Cross-Section'
        ostream.print_header(title.ljust(width))
        ostream.print_header(('=' * len(title)).ljust(width))
        ostream.print_blank()

        freqs = rsp_results['frequencies']

        if len(freqs) == 1 and freqs[0] == 0.0:
            text = '*** No linear absorption spectrum at zero frequency.'
            ostream.print_header(text.ljust(width))
            ostream.print_blank()
            return

        title = 'Reference: '
        title += 'J. Kauczor and P. Norman, '
        title += 'J. Chem. Theory Comput. 2014, 10, 2449-2455.'
        ostream.print_header(title.ljust(width))
        ostream.print_blank()

        assert_msg_critical(
            '[a.u.]' in spectrum['x_label'],
            'ComplexResponse._print_absorption_results: In valid unit in x_label'
        )
        assert_msg_critical(
            '[a.u.]' in spectrum['y_label'],
            'ComplexResponse._print_absorption_results: In valid unit in y_label'
        )

        title = '{:<20s}{:<20s}{:>15s}'.format('Frequency[a.u.]',
                                               'Frequency[eV]',
                                               'sigma(w)[a.u.]')
        ostream.print_header(title.ljust(width))
        ostream.print_header(('-' * len(title)).ljust(width))

        for w, sigma in zip(spectrum['x_data'], spectrum['y_data']):
            output = '{:<20.4f}{:<20.5f}{:>13.8f}'.format(
                w, w * hartree_in_ev(), sigma)
            ostream.print_header(output.ljust(width))

        ostream.print_blank()
        ostream.flush()

    def _print_ecd_results(self, rsp_results, ostream=None):
        """
        Prints ECD results to output stream.

        :param rsp_results:
            The dictionary containing response results.
        :param ostream:
            The output stream.
        """

        if ostream is None:
            ostream = self.ostream

        width = 92

        spectrum = self.get_spectrum(rsp_results, 'au')

        title = 'Circular Dichroism Spectrum'
        ostream.print_header(title.ljust(width))
        ostream.print_header(('=' * len(title)).ljust(width))
        ostream.print_blank()

        freqs = rsp_results['frequencies']

        if len(freqs) == 1 and freqs[0] == 0.0:
            text = '*** No circular dichroism spectrum at zero frequency.'
            ostream.print_header(text.ljust(width))
            ostream.print_blank()
            return

        title = 'Reference: '
        title += 'A. Jiemchooroj and P. Norman, '
        title += 'J. Chem. Phys. 126, 134102 (2007).'
        ostream.print_header(title.ljust(width))
        ostream.print_blank()

        assert_msg_critical(
            '[a.u.]' in spectrum['x_label'],
            'ComplexResponse._print_ecd_results: In valid unit in x_label')
        assert_msg_critical(
            r'[L mol$^{-1}$ cm$^{-1}$]' in spectrum['y_label'],
            'ComplexResponse._print_ecd_results: In valid unit in y_label')

        title = '{:<20s}{:<20s}{:>28s}'.format('Frequency[a.u.]',
                                               'Frequency[eV]',
                                               'Delta_epsilon[L mol^-1 cm^-1]')
        ostream.print_header(title.ljust(width))
        ostream.print_header(('-' * len(title)).ljust(width))

        for w, Delta_epsilon in zip(spectrum['x_data'], spectrum['y_data']):
            output = '{:<20.4f}{:<20.5f}{:>18.8f}'.format(
                w, w * hartree_in_ev(), Delta_epsilon)
            ostream.print_header(output.ljust(width))

        ostream.print_blank()
        ostream.flush()

    def _write_checkpoint(self, molecule, basis, dft_dict, pe_dict, labels):
        """
        Writes checkpoint file. Copied from linearsolver. Changed to work
        without ungerade.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param dft_dict:
            The dictionary containing DFT information.
        :param pe_dict:
            The dictionary containing PE information.
        :param labels:
            The list of labels.
        """

        if self.checkpoint_file is None and self.filename is None:
            return

        if self.checkpoint_file is None and self.filename is not None:
            self.checkpoint_file = f'{self.filename}.rsp.h5'

        t0 = tm.time()

        if self.rank == mpi_master():
            success = write_rsp_hdf5(self.checkpoint_file, [], [], molecule,
                                     basis, dft_dict, pe_dict, self.ostream)
        else:
            success = False
        success = self.comm.bcast(success, root=mpi_master())

        if success:
            if self.nonlinear:
                dist_arrays = [
                    self._dist_bger, self._dist_e2bger, self._dist_fock_ger
                ]
            else:
                dist_arrays = [self._dist_bger, self._dist_e2bger]

            for dist_array, label in zip(dist_arrays, labels):
                dist_array.append_to_hdf5_file(self.checkpoint_file, label)

            checkpoint_text = 'Time spent in writing checkpoint file: '
            checkpoint_text += f'{(tm.time() - t0):.2f} sec'
            self.ostream.print_info(checkpoint_text)
            self.ostream.print_blank()

    def _read_checkpoint(self, rsp_vector_labels):
        """
        Reads distributed arrays from checkpoint file. Copied from
        linearsolver.py and adjusted for TDA

        :param rsp_vector_labels:
            The list of labels of vectors.
        """

        dist_arrays = [
            DistributedArray.read_from_hdf5_file(self.checkpoint_file, label,
                                                 self.comm)
            for label in rsp_vector_labels
        ]

        if self.nonlinear:
            (self._dist_bger, self._dist_e2bger,
             self._dist_fock_ger) = dist_arrays
        else:
            (self._dist_bger, self._dist_e2bger) = dist_arrays

        checkpoint_text = 'Restarting from checkpoint file: '
        checkpoint_text += self.checkpoint_file
        self.ostream.print_info(checkpoint_text)
        self.ostream.print_blank()

    def remove_linear_dependence(self, basis, threshold):
        """
        Removes linear dependence in a set of vectors.
        Based on the function in linearsolver.py, modified to work with full size
        distributed arrays.

        :param basis:
            The set of vectors.
        :param threshold:
            The threshold for removing linear dependence.

        :return:
            The new set of vectors.
        """

        Sb = basis.matmul_AtB(basis)
        if self.rank == mpi_master():
            l, T = np.linalg.eigh(Sb)
            b_norm = np.sqrt(Sb.diagonal())
            mask = l > b_norm * threshold
            Tmask = T[:, mask].copy()
        else:
            Tmask = None
        Tmask = self.comm.bcast(Tmask, root=mpi_master())

        return basis.matmul_AB_no_gather(Tmask)

    @staticmethod
    def orthogonalize_gram_schmidt(tvecs):
        """
        Applies modified Gram Schmidt orthogonalization to trial vectors.
        Based on the function in linearsolver.py, modified to work with full size
        distributed arrays.

        :param tvecs:
            The trial vectors.

        :return:
            The orthogonalized trial vectors.
        """

        if tvecs.shape(1) > 0:

            n2 = tvecs.dot(0, tvecs, 0)
            f = 1.0 / n2
            tvecs.data[:, 0] *= f

            for i in range(1, tvecs.shape(1)):
                for j in range(i):
                    dot_ij = tvecs.dot(i, tvecs, j)
                    dot_jj = tvecs.dot(j, tvecs, j)
                    f = dot_ij / dot_jj
                    tvecs.data[:, i] -= f * tvecs.data[:, j]

                n2 = tvecs.dot(i, tvecs, i)
                f = 1.0 / n2
                tvecs.data[:, i] *= f

        return tvecs

    @staticmethod
    def normalize(vecs):
        """
        Normalizes vectors by dividing by vector norm.
        Based on the function in linearsolver.py, modified to work with full size
        distributed arrays.

        :param vecs:
            The vectors.

        :param Retruns:
            The normalized vectors.
        """

        invnorm = 1 / vecs.norm(axis=0)

        vecs.data *= invnorm

        return vecs
