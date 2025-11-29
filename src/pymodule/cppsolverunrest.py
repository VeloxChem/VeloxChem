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
import h5py

from .veloxchemlib import (mpi_master, hartree_in_wavenumber, hartree_in_ev,
                           hartree_in_inverse_nm, fine_structure_constant,
                           extinction_coefficient_from_beta, avogadro_constant,
                           bohr_in_angstrom)
from .outputstream import OutputStream
from .profiler import Profiler
from .distributedarray import DistributedArray
from .linearsolver import LinearSolver
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           dft_sanity_check, pe_sanity_check,
                           solvation_model_sanity_check)
from .errorhandler import assert_msg_critical, safe_solve
from .checkpoint import (check_rsp_hdf5, write_rsp_solution_with_multiple_keys)
from .inputparser import parse_seq_fixed

try:
    import matplotlib.pyplot as plt
    from scipy.interpolate import make_interp_spline
except ImportError:
    pass


class ComplexResponseUnrestricted(LinearSolver):
    """
    Implements the complex linear response solver.

    # vlxtag: UHF, Absorption, CPP
    # vlxtag: UKS, Absorption, CPP
    # vlxtag: UHF, ECD, CPP
    # vlxtag: UKS, ECD, CPP

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
        pc = precond.data[:, 2]
        pd = precond.data[:, 3]

        v_in_rg = v_in.data[:, 0]
        v_in_ru = v_in.data[:, 1]
        v_in_iu = v_in.data[:, 2]
        v_in_ig = v_in.data[:, 3]

        v_out_rg = pa * v_in_rg + pb * v_in_ru + pc * v_in_iu + pd * v_in_ig
        v_out_ru = pb * v_in_rg + pa * v_in_ru + pd * v_in_iu + pc * v_in_ig
        v_out_iu = pc * v_in_rg + pd * v_in_ru - pa * v_in_iu - pb * v_in_ig
        v_out_ig = pd * v_in_rg + pc * v_in_ru - pb * v_in_iu - pa * v_in_ig

        v_mat = np.hstack((
            v_out_rg.reshape(-1, 1),
            v_out_ru.reshape(-1, 1),
            v_out_iu.reshape(-1, 1),
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
        trials_ung = []

        for (op, w), vec in vectors.items():
            v = self._preconditioning(precond[w], vec)
            norms_2 = 2.0 * v.squared_norm(axis=0)
            vn = np.sqrt(np.sum(norms_2))

            if vn > self.norm_thresh:
                norms = np.sqrt(norms_2)
                # real gerade
                if norms[0] > self.norm_thresh:
                    trials_ger.append(v.data[:, 0])
                # real ungerade
                if norms[1] > self.norm_thresh:
                    trials_ung.append(v.data[:, 1])
                # imaginary ungerade
                if norms[2] > self.norm_thresh:
                    trials_ung.append(v.data[:, 2])
                # imaginary gerade
                if norms[3] > self.norm_thresh:
                    trials_ger.append(v.data[:, 3])

        new_ger = np.array(trials_ger).T
        new_ung = np.array(trials_ung).T

        dist_new_ger = DistributedArray(new_ger, self.comm, distribute=False)
        dist_new_ung = DistributedArray(new_ung, self.comm, distribute=False)

        return dist_new_ger, dist_new_ung

    def compute(self, molecule, basis, scf_results, v_grad=None):
        """
        Solves for the response vector iteratively while checking the residuals
        for convergence.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
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
        self._dist_bung = None
        self._dist_e2bger = None
        self._dist_e2bung = None

        self.nonlinear = False
        self._dist_fock_ger = None
        self._dist_fock_ung = None

        # make sure that cpp_flag is properly set
        if self.cpp_flag is not None:
            self.set_cpp_flag(self.cpp_flag)

        # check molecule
        molecule_sanity_check(molecule)

        # check SCF results
        scf_results_sanity_check(self, scf_results)

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

        nocc_a = molecule.number_of_alpha_electrons()
        nocc_b = molecule.number_of_beta_electrons()

        # ERI information
        eri_dict = self._init_eri(molecule, basis)

        # DFT information
        dft_dict = self._init_dft(molecule, scf_results)

        # PE information
        pe_dict = self._init_pe(molecule, basis)

        # CPCM information
        self._init_cpcm(molecule)

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
            self.restart = self.comm.bcast(self.restart, root=mpi_master())

        # read initial guess from restart file
        if self.restart:
            self._read_checkpoint(rsp_vector_labels)

            # TODO: handle restarting with different frequencies

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
            self._cur_iter = iteration

            n_ger = self._dist_bger.shape(1)
            n_ung = self._dist_bung.shape(1)

            e2gg = self._dist_bger.matmul_AtB(self._dist_e2bger)
            e2uu = self._dist_bung.matmul_AtB(self._dist_e2bung)
            s2ug = self._dist_bung.matmul_AtB(self._dist_bger)

            for op, w in op_freq_keys:
                if (iteration == 0 or
                        relative_residual_norm[(op, w)] > self.conv_thresh):

                    grad_rg = dist_grad[(op, w)].get_column(0)
                    grad_ru = dist_grad[(op, w)].get_column(1)
                    grad_iu = dist_grad[(op, w)].get_column(2)
                    grad_ig = dist_grad[(op, w)].get_column(3)

                    # projections onto gerade and ungerade subspaces:

                    g_realger = self._dist_bger.matmul_AtB(grad_rg)
                    g_imagger = self._dist_bger.matmul_AtB(grad_ig)
                    g_realung = self._dist_bung.matmul_AtB(grad_ru)
                    g_imagung = self._dist_bung.matmul_AtB(grad_iu)

                    # creating gradient and matrix for linear equation

                    size = 2 * (n_ger + n_ung)

                    if self.rank == mpi_master():

                        # gradient

                        g = np.zeros(size)

                        g[:n_ger] = g_realger[:]
                        g[n_ger:n_ger + n_ung] = g_realung[:]
                        g[n_ger + n_ung:size - n_ger] = -g_imagung[:]
                        g[size - n_ger:] = -g_imagger[:]

                        # matrix

                        mat = np.zeros((size, size))

                        # filling E2gg

                        mat[:n_ger, :n_ger] = e2gg[:, :]
                        mat[size - n_ger:, size - n_ger:] = -e2gg[:, :]

                        # filling E2uu

                        mat[n_ger:n_ger + n_ung,
                            n_ger:n_ger + n_ung] = e2uu[:, :]

                        mat[n_ger + n_ung:size - n_ger,
                            n_ger + n_ung:size - n_ger] = -e2uu[:, :]

                        # filling S2ug

                        mat[n_ger:n_ger + n_ung, :n_ger] = -w * s2ug[:, :]

                        mat[n_ger + n_ung:size - n_ger, :n_ger] = d * s2ug[:, :]

                        mat[n_ger:n_ger + n_ung, size - n_ger:] = d * s2ug[:, :]

                        mat[n_ger + n_ung:size - n_ger,
                            size - n_ger:] = w * s2ug[:, :]

                        # filling S2ug.T (interchanging of row and col)

                        mat[:n_ger, n_ger:n_ger + n_ung] = -w * s2ug.T[:, :]

                        mat[:n_ger,
                            n_ger + n_ung:size - n_ger] = d * s2ug.T[:, :]

                        mat[size - n_ger:,
                            n_ger:n_ger + n_ung] = d * s2ug.T[:, :]

                        mat[size - n_ger:,
                            n_ger + n_ung:size - n_ger] = w * s2ug.T[:, :]

                        # solving matrix equation

                        c = safe_solve(mat, g)
                    else:
                        c = None
                    c = self.comm.bcast(c, root=mpi_master())

                    # extracting the 4 components of c...

                    c_realger = c[:n_ger]
                    c_realung = c[n_ger:n_ger + n_ung]
                    c_imagung = c[n_ger + n_ung:size - n_ger]
                    c_imagger = c[size - n_ger:]

                    # ...and projecting them onto respective subspace

                    x_realger = self._dist_bger.matmul_AB_no_gather(c_realger)
                    x_realung = self._dist_bung.matmul_AB_no_gather(c_realung)
                    x_imagung = self._dist_bung.matmul_AB_no_gather(c_imagung)
                    x_imagger = self._dist_bger.matmul_AB_no_gather(c_imagger)

                    # composing E2 matrices projected onto solution subspace

                    e2realger = self._dist_e2bger.matmul_AB_no_gather(c_realger)
                    e2imagger = self._dist_e2bger.matmul_AB_no_gather(c_imagger)
                    e2realung = self._dist_e2bung.matmul_AB_no_gather(c_realung)
                    e2imagung = self._dist_e2bung.matmul_AB_no_gather(c_imagung)

                    if self.nonlinear:
                        fock_realger = self._dist_fock_ger.matmul_AB_no_gather(
                            c_realger)
                        fock_imagger = self._dist_fock_ger.matmul_AB_no_gather(
                            c_imagger)
                        fock_realung = self._dist_fock_ung.matmul_AB_no_gather(
                            c_realung)
                        fock_imagung = self._dist_fock_ung.matmul_AB_no_gather(
                            c_imagung)

                        fock_full_data = (
                            fock_realger.data + fock_realung.data - 1j *
                            (fock_imagger.data + fock_imagung.data))

                        focks[(op, w)] = DistributedArray(fock_full_data,
                                                          self.comm,
                                                          distribute=False)

                    # calculating the residual components

                    s2realger = x_realger.data
                    s2imagger = x_imagger.data
                    s2realung = x_realung.data
                    s2imagung = x_imagung.data

                    r_realger = (e2realger.data - w * s2realung +
                                 d * s2imagung - grad_rg.data)
                    r_realung = (e2realung.data - w * s2realger +
                                 d * s2imagger - grad_ru.data)
                    r_imagung = (-e2imagung.data + w * s2imagger +
                                 d * s2realger + grad_iu.data)
                    r_imagger = (-e2imagger.data + w * s2imagung +
                                 d * s2realung + grad_ig.data)

                    r_data = np.hstack((
                        r_realger.reshape(-1, 1),
                        r_realung.reshape(-1, 1),
                        r_imagung.reshape(-1, 1),
                        r_imagger.reshape(-1, 1),
                    ))

                    r = DistributedArray(r_data, self.comm, distribute=False)

                    # calculating relative residual norm
                    # for convergence check

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
                        xvs.append((op, w, xv))

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
        x_realung = solution.get_full_vector(1)
        x_imagung = solution.get_full_vector(2)
        x_imagger = solution.get_full_vector(3)

        if solution.rank == mpi_master():
            x_real = np.hstack((x_realger, x_realger)) + np.hstack(
                (x_realung, -x_realung))
            x_imag = np.hstack((x_imagung, -x_imagung)) + np.hstack(
                (x_imagger, x_imagger))
            return x_real + 1j * x_imag
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

    def plot(self, cpp_results, x_unit='nm', plot_scatter=True):
        """
        Plot absorption or ECD spectrum from the CPP calculation.

        :param cpp_results:
            The dictionary containing CPP results.
        :param x_unit:
            The dictionary containing CPP results.
        """

        assert_msg_critical(
            x_unit.lower() in ['au', 'ev', 'nm'],
            'ComplexResponse.plot: x_unit should be au, ev or nm')

        assert_msg_critical('matplotlib' in sys.modules,
                            'matplotlib is required.')

        assert_msg_critical('scipy' in sys.modules, 'scipy is required.')

        cpp_spec = self.get_spectrum(cpp_results, x_unit)

        if cpp_spec is None:
            print('Nothing to plot for this complex response calculation.')
            return

        if x_unit.lower() == 'nm':
            # make sure the values in x_data are strictly increasing
            x_data = np.array(cpp_spec['x_data'][::-1])
            y_data = np.array(cpp_spec['y_data'][::-1])
        else:
            x_data = np.array(cpp_spec['x_data'])
            y_data = np.array(cpp_spec['y_data'])

        if self.cpp_flag == 'absorption':
            assert_msg_critical(
                '[a.u.]' in cpp_spec['y_label'],
                'ComplexResponse.plot: In valid unit in y_label')
            # Note: use epsilon for absorption
            NA = avogadro_constant()
            a_0 = bohr_in_angstrom() * 1.0e-10
            sigma_to_epsilon = a_0**2 * 10**4 * NA / (np.log(10) * 10**3)
            y_data *= sigma_to_epsilon

        spl = make_interp_spline(x_data, y_data, k=3)
        x_spl = np.linspace(x_data[0], x_data[-1], x_data.size * 10)
        y_spl = spl(x_spl)

        fig, ax = plt.subplots(figsize=(8, 5))

        if x_unit.lower() == 'nm':
            ax.set_xlabel('Wavelength [nm]')
        else:
            ax.set_xlabel(f'Excitation Energy [{x_unit.lower()}]')

        y_max = np.max(np.abs(y_data))

        if self.cpp_flag == 'absorption':
            # Note: use epsilon for absorption
            ax.set_ylabel(r'$\epsilon$ [L mol$^{-1}$ cm$^{-1}$]')
            ax.set_title("Absorption Spectrum")

            ax.set_ylim(0.0, y_max * 1.1)

        elif self.cpp_flag == 'ecd':
            ax.set_ylabel(r'$\Delta \epsilon$ [L mol$^{-1}$ cm$^{-1}$]')
            ax.set_title("ECD Spectrum")

            ax.set_ylim(-y_max * 1.1, y_max * 1.1)

            ax.axhline(y=0,
                       marker=',',
                       color='k',
                       linestyle='-.',
                       markersize=0,
                       linewidth=0.2)

        plt.plot(x_spl, y_spl, color='black', alpha=0.8, linewidth=2.0)

        if plot_scatter:
            plt.scatter(x_data, y_data, color='darkcyan', alpha=0.9, s=15)

        plt.show()

    def _print_results(self, rsp_results, ostream=None):
        """
        Prints response results to output stream.

        :param rsp_results:
            The dictionary containing response results.
        :param ostream:
            The output stream.
        """

        if self.print_level > 1:
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

        # Note: flush is needed at the end of every print method
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

        # Note: flush is needed at the end of every print method
        ostream.flush()

    def write_cpp_rsp_results_to_hdf5(self, fname, rsp_results):
        """
        Writes the results of a linear response calculation to HDF5 file.

        :param fname:
            Name of the HDF5 file.
        :param rsp_results:
            The dictionary containing the linear response results.
        """

        if fname and isinstance(fname, str):

            hf = h5py.File(fname, 'a')

            # Write frequencies
            xlabel = self.group_label + '/frequencies'
            if xlabel in hf:
                del hf[xlabel]
            hf.create_dataset(xlabel, data=rsp_results['frequencies'])

            spectrum = self.get_spectrum(rsp_results, 'au')

            # Write spectrum if an absorption or ecd calculation
            # has been performed. Otherwise there is nothing to write.
            if spectrum is not None:
                y_data = np.array(spectrum['y_data'])

                if self.cpp_flag == 'absorption':
                    assert_msg_critical(
                        '[a.u.]' in spectrum['y_label'],
                        'ComplexResponse.write_cpp_rsp_results_to_hdf5: ' +
                        'In valid unit in y_label')
                    ylabel = self.group_label + '/sigma'
                elif self.cpp_flag == 'ecd':
                    ylabel = self.group_label + '/delta-epsilon'
                if ylabel in hf:
                    del hf[ylabel]
                hf.create_dataset(ylabel, data=y_data)

            hf.close()

    def get_cpp_property_densities(self,
                                   molecule,
                                   basis,
                                   scf_results,
                                   cpp_results,
                                   w,
                                   normalize_densities=True):
        """
        Gets CPP property densities.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_results:
            The SCF results dictionary.
        :param cpp_results:
            The CPP results dictionary.
        :param w:
            The given frequency.
        :param normalize_densities:
            Whether or not to normalize the property densities.

        :return:
            A dictionary containing property densities at given frequency.
        """

        assert_msg_critical(self.cpp_flag in ['absorption', 'ecd'],
                            'get_cpp_property_densities: Invalid cpp_flag')

        assert_msg_critical(
            ('x', w) in cpp_results['solutions'] and
            ('y', w) in cpp_results['solutions'] and
            ('z', w) in cpp_results['solutions'],
            'get_cpp_property_densities: Could not find frequency ' +
            f'{w} in CPP results')

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
            nocc_a = molecule.number_of_alpha_electrons()
            nocc_b = molecule.number_of_beta_electrons()
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

            if self.cpp_flag == 'absorption':
                # vector representation of absorption cross-section
                vec_a = (a_grad_alpha[0] * x_alpha + a_grad_alpha[1] * y_alpha +
                         a_grad_alpha[2] * z_alpha).imag / 3.0
                vec_b = (a_grad_beta[0] * x_beta + a_grad_beta[1] * y_beta +
                         a_grad_beta[2] * z_beta).imag / 3.0
                vec_a *= 4.0 * np.pi * w * fine_structure_constant()
                vec_b *= 4.0 * np.pi * w * fine_structure_constant()
            elif self.cpp_flag == 'ecd':
                # vector representation of Delta epsilon
                vec_a = (a_grad_alpha[0] * x_alpha / w +
                         a_grad_alpha[1] * y_alpha / w +
                         a_grad_alpha[2] * z_alpha / w).imag
                vec_b = (a_grad_beta[0] * x_beta / w +
                         a_grad_beta[1] * x_beta / w +
                         a_grad_beta[2] * x_beta / w).imag
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
