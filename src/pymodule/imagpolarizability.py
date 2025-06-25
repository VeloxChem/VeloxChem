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

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .profiler import Profiler
from .distributedarray import DistributedArray
from .linearsolver import LinearSolver
from .c6driver import C6Driver
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           dft_sanity_check, pe_sanity_check,
                           solvation_model_sanity_check)
from .errorhandler import assert_msg_critical, safe_solve
from .checkpoint import (check_rsp_hdf5, write_rsp_solution_with_multiple_keys)


class ImaginaryPolarizability(C6Driver):
    """
    Implements the imaginary polarizability from the C6 value response solver.

    # vlxtag: RHF, C6, LR
    # vlxtag: RKS, C6, LR

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - a_operator: The A operator.
        - a_components: Cartesian components of the A operator.
        - b_operator: The B operator.
        - b_components: Cartesian components of the B operator.
        - n_points: The number of integration points.
        - w0: The transformation function prefactor.
    """

    def compute(self, molecule, basis, scf_tensors, frequencies):
        """
        Solves for the response vector iteratively while checking the residuals
        for convergence.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param frequencies:
            The imaginary frequencies.

        :return:
            A dictionary containing response functions and solutions.
        """

        if self.norm_thresh is None:
            self.norm_thresh = self.conv_thresh * 1.0e-6
        if self.lindep_thresh is None:
            self.lindep_thresh = self.conv_thresh * 1.0e-6

        self._dist_bger = None
        self._dist_bung = None
        self._dist_e2bger = None
        self._dist_e2bung = None

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

        # check solvation model setup
        if self.rank == mpi_master():
            assert_msg_critical(
                'solvation_model' not in scf_tensors,
                type(self).__name__ + ': Solvation model not implemented')

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
            self._print_header('C6 Value Response Solver',
                               n_points=self.n_points)

        self.start_time = tm.time()

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(nalpha == nbeta,
                            'C6Driver: not implemented for unrestricted case')

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
        b_grad = self.get_complex_prop_grad(self.b_operator, self.b_components,
                                            molecule, basis, scf_tensors)

        points, weights = np.polynomial.legendre.leggauss(self.n_points)
        imagfreqs = frequencies
        imagfreqs = np.append(imagfreqs, 0.0)

        v_grad = None
        if self.rank == mpi_master():
            v_grad = {
                (op, iw): v for op, v in zip(self.b_components, b_grad)
                for iw in imagfreqs
            }

        # operators, frequencies and preconditioners
        if self.rank == mpi_master():
            op_imagfreq_keys = list(v_grad.keys())
        else:
            op_imagfreq_keys = None
        op_imagfreq_keys = self.comm.bcast(op_imagfreq_keys, root=mpi_master())

        precond = {
            iw: self._get_precond(orb_ene, nocc, norb, iw) for iw in imagfreqs
        }

        # distribute the gradient and right-hand side:
        # dist_grad will be used for calculating the subspace matrix
        # equation and residuals, dist_rhs for the initial guess

        dist_grad = {}
        dist_rhs = {}
        for key in op_imagfreq_keys:
            if self.rank == mpi_master():
                gradger, gradung = self._decomp_grad(v_grad[key])
                grad_mat = np.hstack((
                    gradung.real.reshape(-1, 1),
                    gradger.imag.reshape(-1, 1),
                ))
                rhs_mat = np.hstack((
                    gradung.real.reshape(-1, 1),
                    -gradger.imag.reshape(-1, 1),
                ))
            else:
                grad_mat = None
                rhs_mat = None
            dist_grad[key] = DistributedArray(grad_mat, self.comm)
            dist_rhs[key] = DistributedArray(rhs_mat, self.comm)

        rsp_vector_labels = [
            'C6_bger_half_size',
            'C6_bung_half_size',
            'C6_e2bger_half_size',
            'C6_e2bung_half_size',
        ]

        # read initial guess from restart file

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
            bger, bung = self._setup_trials(dist_rhs, precond)

            profiler.set_timing_key('Preparation')

            self._e2n_half_size(bger, bung, molecule, basis, scf_tensors,
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

            xvs = []
            self._cur_iter = iteration

            n_ger = self._dist_bger.shape(1)
            n_ung = self._dist_bung.shape(1)

            e2gg = self._dist_bger.matmul_AtB(self._dist_e2bger, 2.0)
            e2uu = self._dist_bung.matmul_AtB(self._dist_e2bung, 2.0)
            s2ug = self._dist_bung.matmul_AtB(self._dist_bger, 2.0)

            for op, iw in op_imagfreq_keys:
                if (iteration == 0 or
                        relative_residual_norm[(op, iw)] > self.conv_thresh):

                    grad_ru = dist_grad[(op, iw)].get_column(0)
                    grad_ig = dist_grad[(op, iw)].get_column(1)

                    # projections onto gerade and ungerade subspaces:

                    g_realung = self._dist_bung.matmul_AtB(grad_ru, 2.0)

                    # creating gradient and matrix for linear equation

                    size = n_ger + n_ung

                    if self.rank == mpi_master():

                        # gradient

                        g = np.zeros(size)

                        g[:n_ung] = g_realung[:]

                        # matrix

                        mat = np.zeros((size, size))

                        mat[n_ung:n_ung + n_ger,
                            n_ung:n_ung + n_ger] = -e2gg[:, :]
                        mat[:n_ung, :n_ung] = e2uu[:, :]
                        mat[:n_ung, n_ung:n_ung + n_ger] = iw * s2ug[:, :]
                        mat[n_ung:n_ung + n_ger, :n_ung] = iw * s2ug.T[:, :]

                        # solving matrix equation

                        c = safe_solve(mat, g)
                    else:
                        c = None
                    c = self.comm.bcast(c, root=mpi_master())

                    # extracting the 2 components of c...

                    c_realung = c[:n_ung]
                    c_imagger = c[n_ung:n_ung + n_ger]

                    # ...and projecting them onto respective subspace

                    x_realung = self._dist_bung.matmul_AB_no_gather(c_realung)
                    x_imagger = self._dist_bger.matmul_AB_no_gather(c_imagger)

                    # composing E2 matrices projected onto solution subspace

                    e2imagger = self._dist_e2bger.matmul_AB_no_gather(c_imagger)
                    e2realung = self._dist_e2bung.matmul_AB_no_gather(c_realung)

                    # calculating the residual components

                    s2imagger = x_imagger.data
                    s2realung = x_realung.data

                    r_realung = (e2realung.data + iw * s2imagger - grad_ru.data)
                    r_imagger = (-e2imagger.data + iw * s2realung +
                                 grad_ig.data)

                    r_data = np.hstack((
                        r_realung.reshape(-1, 1),
                        r_imagger.reshape(-1, 1),
                    ))

                    r = DistributedArray(r_data, self.comm, distribute=False)

                    # calculating relative residual norm
                    # for convergence check

                    x_data = np.hstack((
                        x_realung.data.reshape(-1, 1),
                        x_imagger.data.reshape(-1, 1),
                    ))

                    x = DistributedArray(x_data, self.comm, distribute=False)

                    x_full = self.get_full_solution_vector(x)
                    if self.rank == mpi_master():
                        xv = np.dot(x_full, v_grad[(op, iw)])
                        xvs.append((op, iw, xv))

                    r_norms_2 = 2.0 * r.squared_norm(axis=0)
                    x_norms_2 = 2.0 * x.squared_norm(axis=0)

                    rn = np.sqrt(np.sum(r_norms_2))
                    xn = np.sqrt(np.sum(x_norms_2))

                    if xn != 0:
                        relative_residual_norm[(op, iw)] = 2.0 * rn / xn
                    else:
                        relative_residual_norm[(op, iw)] = 2.0 * rn

                    if relative_residual_norm[(op, iw)] < self.conv_thresh:
                        solutions[(op, iw)] = x
                    else:
                        residuals[(op, iw)] = r

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

            # creating new sigma and rho linear transformations

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
            self._print_convergence('Complex response')

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of C6 solver')
        profiler.print_memory_usage(self.ostream)

        self._dist_bger = None
        self._dist_bung = None
        self._dist_e2bger = None
        self._dist_e2bung = None

        # calculate response functions
        a_grad = self.get_complex_prop_grad(self.a_operator, self.a_components,
                                            molecule, basis, scf_tensors)

        if self.is_converged:
            if self.rank == mpi_master():
                va = {op: v for op, v in zip(self.a_components, a_grad)}
                rsp_funcs = {}

                # final h5 file for response solutions
                if self.filename is not None:
                    final_h5_fname = f'{self.filename}.h5'
                else:
                    final_h5_fname = None

            for bop, iw in solutions:
                x = self.get_full_solution_vector(solutions[(bop, iw)])

                if self.rank == mpi_master():
                    for aop in self.a_components:
                        rsp_funcs[(aop, bop, iw)] = -np.dot(va[aop], x)

                    # write to h5 file for response solutions
                    if (self.save_solutions and final_h5_fname is not None):
                        solution_keys = [
                            '{:s}_{:s}_{:.8f}'.format(aop, bop, iw)
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

                #c6 = self._integrate_c6(self.w0, points, weights, imagfreqs,
                #                        rsp_funcs)
                c6 = None

                # calculate the isotropic averages
                alphabar = self._calculate_isotropic_average(imagfreqs, rsp_funcs)

                #self._print_results(c6, rsp_funcs, self.ostream)
                self._print_results(c6, alphabar, imagfreqs, rsp_funcs, self.ostream)

                return {
                    #'c6': c6,
                    'response_functions': rsp_funcs,
                    'solutions': solutions,
                    'isotropic_averages': alphabar
                }
            else:
                return {'solutions': solutions}

        return None

    def _integrate_c6(self, w0, points, weights, imagfreqs, rsp_funcs):
        """
        Calculates the C6 value with a Gauss-Legendre quadrature for the
        integral in the Casimir-Polder relation using integration by
        substitution.

        :param w0:
            A constant conversion factor.
        :param points:
            The list of integration points.
        :param weights:
            The list of weights for integration points.
        :param imagfreqs:
            The list of imaginary frequencies.
        :param rsp_funcs:
            The response functions.

        :return:
            The C6 value.
        """

        integral = 0.0

        # note: excluding the last element in imagfreqs (which is 0)
        for ind, iw in enumerate(imagfreqs[:-1]):

            Gxx = rsp_funcs[('x', 'x', iw)].real
            Gyy = rsp_funcs[('y', 'y', iw)].real
            Gzz = rsp_funcs[('z', 'z', iw)].real

            alpha = -(Gxx + Gyy + Gzz) / 3.0
            point = points[ind]
            weight = weights[ind]
            derivative = w0 * 2 / (1 + point)**2
            integral += alpha * alpha * weight * derivative

        # Casimir-Polder relation

        c6 = 3 * integral / math.pi

        return c6

        self.ostream.print_header(output_header.ljust(width))
        self.ostream.print_blank()

    def _calculate_isotropic_average(self, imagfreqs, rsp_funcs):
        """
        Calculates the isotropic average of the polarizability tensor

        :param imagfreqs:
            The list of imaginary frequencies.
        :param rsp_funcs:
            The response functions.

        :return iso_aver:
            The isotropic averages at imaginary frequencies.
        """

        iso_aver = {}

        for iw in imagfreqs:

            Gxx = rsp_funcs[('x', 'x', iw)].real
            Gyy = rsp_funcs[('y', 'y', iw)].real
            Gzz = rsp_funcs[('z', 'z', iw)].real

            alpha = -(Gxx + Gyy + Gzz) / 3.0

            iso_aver[iw] = alpha

        return iso_aver

    def _print_results(self, c6, iso_aver, frequencies, rsp_funcs, ostream):
        """
        Prints response property to output stream.

        :param c6:
            The C6 value.
        :param frequencies:
            The imaginary frequencies.
        :param rsp_funcs:
            The response functions.
        :param iso_aver:
            The isotropic averages of the polarizability tensor.
        :param ostream:
            The output stream.
        """

        width = 92

        title = 'Response Functions at Given Imaginary Frequencies'
        ostream.print_header(title.ljust(width))
        ostream.print_header(('=' * len(title)).ljust(width))
        ostream.print_blank()

        points, weights = np.polynomial.legendre.leggauss(self.n_points)
        imagfreqs = frequencies
        printfreqs = imagfreqs #np.append(imagfreqs, 0.0)

        for iw in printfreqs:
            title = '{:<7s} {:<7s} {:>10s} {:>15s} {:>16s}'.format(
                'Dipole', 'Dipole', 'Frequency', 'Real', 'Imaginary')
            ostream.print_header(title.ljust(width))
            ostream.print_header(('-' * len(title)).ljust(width))

            for a in self.a_components:
                for b in self.b_components:
                    rsp_func_val = rsp_funcs[(a, b, iw)]
                    ops_label = '<<{:>3s}  ;  {:<3s}>> {:10.4f}'.format(
                        a.lower(), b.lower(), iw)
                    output = '{:<15s} {:15.8f} {:15.8f}j'.format(
                        ops_label, rsp_func_val.real, rsp_func_val.imag)
                    ostream.print_header(output.ljust(width))
            ostream.print_blank()

        title = 'Isotropic Averages at Given Imaginary Frequencies'
        ostream.print_header(title.ljust(width))
        ostream.print_header(('=' * len(title)).ljust(width))
        ostream.print_blank()

        for iw in imagfreqs:
            title = '{:>10s} {:>15s}'.format(
                'Frequency', 'Alpha')
            ostream.print_header(title.ljust(width))
            ostream.print_header(('-' * len(title)).ljust(width))

            alpha_val = iso_aver[iw]
            output = '{:10.4f} {:15.8f}'.format(
                iw, alpha_val)
            ostream.print_header(output.ljust(width))
            ostream.print_blank()

        #title = 'C6 Dispersion Coefficient'
        #ostream.print_header(title.ljust(width))
        #ostream.print_header(('=' * len(title)).ljust(width))
        #ostream.print_blank()

        #title = 'Reference: '
        #title += 'Amos et al., '
        #title += 'J. Chem. Phys. 89, 2186 (1985).'
        #ostream.print_header(title.ljust(width))
        #ostream.print_blank()

        #Gxx_i0 = rsp_funcs[('x', 'x', 0.0)].real
        #Gyy_i0 = rsp_funcs[('y', 'y', 0.0)].real
        #Gzz_i0 = rsp_funcs[('z', 'z', 0.0)].real

        #alpha_i0 = -(Gxx_i0 + Gyy_i0 + Gzz_i0) / 3.0

        #output = 'Homomolecular C_6 value        :    {:10.6f} a.u.'.format(c6)
        #ostream.print_header(output.ljust(width))
        #ostream.print_blank()
        #output = 'Static polarizability alpha(0) :    {:10.6f} a.u.'.format(
        #    alpha_i0)
        ostream.print_header(output.ljust(width))

        ostream.print_blank()
        ostream.flush()
