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
import sys

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .profiler import Profiler
from .distributedarray import DistributedArray
from .subcommunicators import SubCommunicators
from .linearsolver import LinearSolver
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           dft_sanity_check, pe_sanity_check)
from .errorhandler import assert_msg_critical, safe_solve
from .inputparser import parse_input
from .checkpoint import write_rsp_hdf5, check_rsp_hdf5
from .batchsize import get_batch_size
from .batchsize import get_number_of_batches


class CphfSolver(LinearSolver):
    """
    Implements solver for the coupled-perturbed Hartree-Fock (CPHF) equations.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes CPHF solver to default setup.

        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        super().__init__(comm, ostream)

        self.print_residuals = False
        self.max_iter = 150

        self.orbrsp_type = 'default'
        self._method_type = 'restricted'

        self._input_keywords['orbitalresponse'] = {
            'print_residuals': ('bool', 'print iteration to output'),
            'max_iter': ('int', 'maximum number of iterations'),
            'force_checkpoint':
                ('bool', 'flag for writing checkpoint every iteration'),
        }

    def update_settings(self, cphf_dict, method_dict=None):
        """
        Updates response and method settings in CPHF solver.

        :param cphf_dict:
            The dictionary of CPHF (orbital response) settings.
        :param method_dict:
            The dictionary of method settings.
        """

        super().update_settings(cphf_dict, method_dict)

        cphf_keywords = {
            key: val[0]
            for key, val in self._input_keywords['orbitalresponse'].items()
        }

        parse_input(self, cphf_keywords, cphf_dict)

    def compute_solution_vectors(self, molecule, basis, scf_tensors, dist_rhs,
                                 eri_dict, dft_dict, pe_dict):
        """
        Performs CPHF calculation for a specific atom of a molecule and a basis set.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param dist_rhs:
            The right-hand side as a list of distributed arrays.

        :return:
            A dictionary containing the RHS and solution (ov block)
            of the CPHF equations.
        """

        if self.norm_thresh is None:
            self.norm_thresh = self.conv_thresh * 1.0e-6
        if self.lindep_thresh is None:
            self.lindep_thresh = self.conv_thresh * 1.0e-2

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'CphfSolver: not implemented for unrestricted case')

        self.start_time = tm.time()

        # check molecule
        molecule_sanity_check(molecule)

        # check SCF results
        scf_results_sanity_check(self, scf_tensors)

        # update checkpoint_file after scf_results_sanity_check
        if self.filename is not None and self.checkpoint_file is None:
            self.checkpoint_file = f'{self.filename}_orbrsp.h5'
        elif (self.checkpoint_file is not None and
              self.checkpoint_file.endswith('_rsp.h5')):
            self.checkpoint_file = (self.checkpoint_file[:-len('_rsp.h5')] +
                                    '_orbrsp.h5')

        # check dft setup
        dft_sanity_check(self, 'compute')

        # check pe setup
        pe_sanity_check(self, molecule=molecule)

        if self.rank == mpi_master():
            if self._dft:
                self.print_cphf_header('Coupled-Perturbed Kohn-Sham Solver')
            else:
                self.print_cphf_header('Coupled-Perturbed Hartree-Fock Solver')

        if self.rank == mpi_master():
            mo_energies = scf_tensors['E_alpha']
            # nmo is sometimes different than nao (because of linear
            # dependencies which get removed during SCF)
            nmo = mo_energies.shape[0]
            nocc = molecule.number_of_alpha_electrons()
            nvir = nmo - nocc
            nao = scf_tensors['C_alpha'].shape[0]
            eocc = mo_energies[:nocc]
            evir = mo_energies[nocc:]
            eov = eocc.reshape(-1, 1) - evir
        else:
            mo_energies = None
            eov = None
            nao = None

        mo_energies = self.comm.bcast(mo_energies, root=mpi_master())
        eov = self.comm.bcast(eov, root=mpi_master())
        nao = self.comm.bcast(nao, root=mpi_master())

        nmo = mo_energies.shape[0]
        nocc = molecule.number_of_alpha_electrons()
        nvir = nmo - nocc

        dof = len(dist_rhs)

        # Initialize trial and sigma vectors
        self.dist_trials = None
        self.dist_sigmas = None

        # the preconditioner: 1 / (eocc - evir)
        precond = (1.0 / eov).reshape(nocc * nvir)
        dist_precond = DistributedArray(precond, self.comm)

        # no restart for atom calculation

        # setup and precondition trial vectors
        dist_trials = self.setup_trials(molecule, dist_precond, dist_rhs)

        # construct the sigma (E*t) vectors
        self.build_sigmas(molecule, basis, scf_tensors, dist_trials, eri_dict,
                          dft_dict, pe_dict)

        # lists that will hold the solutions and residuals
        # TODO: double check residuals in setup_trials
        solutions = [None for x in range(dof)]
        residuals = list(np.zeros((dof)))
        relative_residual_norm = [None for x in range(dof)]

        # start iterations
        for iteration in range(self.max_iter):

            # Orbital Hessian in reduced subspace
            orbhess_red = self.dist_trials.matmul_AtB(self.dist_sigmas)

            self._cur_iter = iteration
            num_vecs = self.dist_trials.shape(1)

            for x in range(dof):
                if (iteration > 0 and
                        relative_residual_norm[x] < self.conv_thresh):
                    continue

                # CPHF RHS in reduced space
                cphf_rhs_red = self.dist_trials.matmul_AtB(dist_rhs[x])

                if self.rank == mpi_master():
                    # solve the equations exactly in the subspace
                    u_red = safe_solve(orbhess_red, cphf_rhs_red)
                else:
                    u_red = None
                u_red = self.comm.bcast(u_red, root=mpi_master())

                # solution vector in full space
                u = self.dist_trials.matmul_AB_no_gather(u_red)

                # calculate residuals
                sigmas_u_red = self.dist_sigmas.matmul_AB_no_gather(u_red)
                residual = sigmas_u_red.data - dist_rhs[x].data

                # make distributed array out of residuals and current vectors
                dist_residual = DistributedArray(residual,
                                                 self.comm,
                                                 distribute=False)

                # calculate norms
                r_norm = np.sqrt(dist_residual.squared_norm(axis=0))
                u_norm = np.sqrt(u.squared_norm(axis=0))

                if u_norm != 0:
                    relative_residual_norm[x] = r_norm / u_norm
                else:
                    relative_residual_norm[x] = r_norm

                if relative_residual_norm[x] < self.conv_thresh:
                    solutions[x] = u
                else:
                    residuals[x] = dist_residual

            # write to output
            if self.rank == mpi_master():
                self.ostream.print_info(
                    '{:d} trial vectors in reduced space'.format(num_vecs))
                self.ostream.print_blank()

                self.print_iteration(relative_residual_norm, molecule)

            # check convergence
            self.check_convergence(relative_residual_norm)

            if self.is_converged:
                break

            # update trial vectors
            new_trials = self.setup_trials(molecule, dist_precond, residuals,
                                           self.dist_trials)

            if self.rank == mpi_master():
                n_new_trials = new_trials.shape(1)
            else:
                n_new_trials = None
            n_new_trials = self.comm.bcast(n_new_trials, root=mpi_master())

            # update sigma vectors
            self.build_sigmas(molecule, basis, scf_tensors, new_trials,
                              eri_dict, dft_dict, pe_dict)

        # converged?
        if self.rank == mpi_master():
            if self._dft:
                self._print_convergence('Coupled-Perturbed Kohn-Sham')
            else:
                self._print_convergence('Coupled-Perturbed Hartree-Fock')

        # only return the solution
        return solutions

    def compute(self, molecule, basis, scf_tensors, *args):
        """
        Performs CPHF calculation for a molecule and a basis set.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param *args:
            Additional arguments, same as in compute_rhs.

        :return:
            A dictionary containing the RHS and solution (ov block)
            of the CPHF equations, the derivative of the AO Fock matrix,
            partial derivative of the overlap matrix (oo block),
            and an auxiliary Fock matrix (oo block of the CPHF coefficients
            contracted with the two-electron integrals).
        """

        if self.norm_thresh is None:
            self.norm_thresh = self.conv_thresh * 1.0e-6
        if self.lindep_thresh is None:
            self.lindep_thresh = self.conv_thresh * 1.0e-2

        # check print level (verbosity of output)
        self.print_level = max(1, min(self.print_level, 3))

        self.cphf_results = self.compute_subspace_solver(
            molecule, basis, scf_tensors, *args)

    def compute_subspace_solver(self, molecule, basis, scf_tensors, *args):
        """
        Performs CPHF calculation for a molecule and a basis set.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            To be determined. (Dict with RHS, CPHF coefficients, ...)
        """

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        self.start_time = tm.time()

        # check molecule
        molecule_sanity_check(molecule)

        # check SCF results
        scf_results_sanity_check(self, scf_tensors)

        # update checkpoint_file after scf_results_sanity_check
        if self.filename is not None and self.checkpoint_file is None:
            self.checkpoint_file = f'{self.filename}_orbrsp.h5'
        elif (self.checkpoint_file is not None and
              self.checkpoint_file.endswith('_rsp.h5')):
            self.checkpoint_file = (self.checkpoint_file[:-len('_rsp.h5')] +
                                    '_orbrsp.h5')

        # check dft setup
        dft_sanity_check(self, 'compute')

        # check pe setup
        pe_sanity_check(self, molecule=molecule)

        if self.rank == mpi_master():
            if self._dft:
                self.print_cphf_header('Coupled-Perturbed Kohn-Sham Solver')
            else:
                self.print_cphf_header('Coupled-Perturbed Hartree-Fock Solver')

        # ERI information
        eri_dict = self._init_eri(molecule, basis)

        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors)

        # PE information
        pe_dict = self._init_pe(molecule, basis)

        cphf_rhs_dict = self.compute_rhs(molecule, basis, scf_tensors, eri_dict,
                                         dft_dict, pe_dict, *args)

        dist_rhs = cphf_rhs_dict['dist_cphf_rhs']
        dof = len(dist_rhs)

        # Initialize trial and sigma vectors
        self.dist_trials = None
        self.dist_sigmas = None

        dist_precond = self._get_precond(molecule, basis, scf_tensors,
                                         self._method_type)

        orbrsp_vector_labels = [
            self.orbrsp_type.upper() + '_orbrsp_trials',
            self.orbrsp_type.upper() + '_orbrsp_sigmas',
        ]

        # check validity of checkpoint file
        if self.restart:
            if self.rank == mpi_master():
                self.restart = check_rsp_hdf5(self.checkpoint_file,
                                              orbrsp_vector_labels, molecule,
                                              basis, dft_dict, pe_dict)
            self.restart = self.comm.bcast(self.restart, root=mpi_master())

        # read initial guess from restart file
        if self.restart:
            self._read_checkpoint_for_cphf(orbrsp_vector_labels)
        else:
            # setup and precondition trial vectors
            dist_trials = self.setup_trials(molecule, dist_precond, dist_rhs)

            profiler.set_timing_key('Preparation')

            # construct the sigma (E*t) vectors
            self.build_sigmas(molecule, basis, scf_tensors, dist_trials,
                              eri_dict, dft_dict, pe_dict, profiler,
                              self._method_type)

        # lists that will hold the solutions and residuals
        # TODO: double check residuals in setup_trials
        solutions = [None for x in range(dof)]
        residuals = list(np.zeros((dof)))
        relative_residual_norm = [None for x in range(dof)]

        iter_per_trial_in_hours = None

        # start iterations
        for iteration in range(self.max_iter):

            iter_start_time = tm.time()

            profiler.set_timing_key(f'Iteration {iteration + 1}')
            profiler.start_timer('ReducedSpace')

            # Orbital Hessian in reduced subspace
            orbhess_red = self.dist_trials.matmul_AtB(self.dist_sigmas)

            self._cur_iter = iteration
            num_vecs = self.dist_trials.shape(1)

            for x in range(dof):
                if (iteration > 0 and
                        relative_residual_norm[x] < self.conv_thresh):
                    continue

                # CPHF RHS in reduced space
                cphf_rhs_red = self.dist_trials.matmul_AtB(dist_rhs[x])

                if self.rank == mpi_master():
                    # solve the equations exactly in the subspace
                    u_red = safe_solve(orbhess_red, cphf_rhs_red)
                else:
                    u_red = None
                u_red = self.comm.bcast(u_red, root=mpi_master())

                # solution vector in full space
                u = self.dist_trials.matmul_AB_no_gather(u_red)

                # calculate residuals
                sigmas_u_red = self.dist_sigmas.matmul_AB_no_gather(u_red)
                residual = sigmas_u_red.data - dist_rhs[x].data

                # make distributed array out of residuals and current vectors
                dist_residual = DistributedArray(residual,
                                                 self.comm,
                                                 distribute=False)

                # calculate norms
                r_norm = np.sqrt(dist_residual.squared_norm(axis=0))
                u_norm = np.sqrt(u.squared_norm(axis=0))

                if u_norm != 0:
                    relative_residual_norm[x] = r_norm / u_norm
                else:
                    relative_residual_norm[x] = r_norm

                if relative_residual_norm[x] < self.conv_thresh:
                    solutions[x] = u
                else:
                    residuals[x] = dist_residual

            # write to output
            if self.rank == mpi_master():
                self.ostream.print_info(
                    '{:d} trial vectors in reduced space'.format(num_vecs))
                self.ostream.print_blank()

                if self.print_level > 1:
                    profiler.print_memory_subspace(
                        {
                            'dist_trials': self.dist_trials,
                            'dist_sigmas': self.dist_sigmas,
                            'dist_precond': dist_precond,
                            'solutions': solutions,
                            'residuals': residuals,
                        }, self.ostream)

                profiler.check_memory_usage(
                    'Iteration {:d} subspace'.format(iteration + 1))

                profiler.print_memory_tracing(self.ostream)

                self.print_iteration(relative_residual_norm, molecule)

            profiler.stop_timer('ReducedSpace')

            # check convergence
            self.check_convergence(relative_residual_norm)

            if self.is_converged:
                break

            profiler.start_timer('Orthonorm.')

            # update trial vectors
            new_trials = self.setup_trials(molecule, dist_precond, residuals,
                                           self.dist_trials)

            profiler.stop_timer('Orthonorm.')

            if self.rank == mpi_master():
                n_new_trials = new_trials.shape(1)
            else:
                n_new_trials = None
            n_new_trials = self.comm.bcast(n_new_trials, root=mpi_master())

            if iter_per_trial_in_hours is not None:
                next_iter_in_hours = iter_per_trial_in_hours * n_new_trials
                if self._need_graceful_exit(next_iter_in_hours):
                    self._graceful_exit_for_cphf(molecule, basis, dft_dict,
                                                 pe_dict, orbrsp_vector_labels)

            if self.force_checkpoint:
                self._write_checkpoint_for_cphf(molecule, basis, dft_dict,
                                                pe_dict, orbrsp_vector_labels)

            # update sigma vectors
            self.build_sigmas(molecule, basis, scf_tensors, new_trials,
                              eri_dict, dft_dict, pe_dict, profiler,
                              self._method_type)

            iter_in_hours = (tm.time() - iter_start_time) / 3600
            iter_per_trial_in_hours = iter_in_hours / n_new_trials

            profiler.check_memory_usage(
                'Iteration {:d} sigma build'.format(iteration + 1))

        self._write_checkpoint_for_cphf(molecule, basis, dft_dict, pe_dict,
                                        orbrsp_vector_labels)

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage(
            'End of ComputeCPHF iterative subspace algorithm')
        profiler.print_memory_usage(self.ostream)

        # converged?
        if self.rank == mpi_master():
            if self._dft:
                self._print_convergence('Coupled-Perturbed Kohn-Sham')
            else:
                self._print_convergence('Coupled-Perturbed Hartree-Fock')

        # merge the rhs dict with the solution
        return {
            **cphf_rhs_dict,
            'dist_cphf_ov': solutions,
        }

    def build_sigmas(self,
                     molecule,
                     basis,
                     scf_tensors,
                     dist_trials,
                     eri_dict,
                     dft_dict,
                     pe_dict,
                     profiler=None,
                     method_type='restricted'):
        """
        Compute the matrix-vector product corresponding to the
        left-hand-side of the CPHF/CPKS equations.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param scf_tensors:
            The dictionary of tensors from a converged SCF calculation.
        :param dist_trials:
            The trial vectors as a list of distributed arrays.
        :param eri_dict:
            The dictionary containing ERI information.
        :param dft_dict:
            The dictionary containing DFT information.
        :param pe_dict:
            The dictionary containing PE information.
        :param profiler:
            The profiler.
        :param method_type:
            The method type (restricted, unrestricted).
        """

        if self.use_subcomms:
            if method_type == 'restricted':
                self.build_sigmas_subcomms(molecule, basis, scf_tensors,
                                           dist_trials, eri_dict, dft_dict,
                                           pe_dict, profiler)
            else:
                # TODO: enable subcomms for unrestricted
                assert_msg_critical(
                    False, 'CphfSolver.build_sigmas: ' +
                    'Cannot use subcomms for unrestricted case.')
        else:
            if method_type == 'restricted':
                self.build_sigmas_single_comm(molecule, basis, scf_tensors,
                                              dist_trials, eri_dict, dft_dict,
                                              pe_dict, profiler)
            else:
                self.build_sigmas_single_comm_unrestricted(
                    molecule, basis, scf_tensors, dist_trials, eri_dict,
                    dft_dict, pe_dict, profiler)

    def build_sigmas_subcomms(self,
                              molecule,
                              basis,
                              scf_tensors,
                              dist_trials,
                              eri_dict,
                              dft_dict,
                              pe_dict,
                              profiler=None):
        """
        Compute the matrix-vector product corresponding to the
        left-hand-side of the CPHF/CPKS equations.
        Appends sigma and trial vectors to member variable.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param dist_trials:
            Distributed array of trial vectors
        """

        if self.rank == mpi_master():
            mo = scf_tensors['C_alpha']
            mo_energies = scf_tensors['E_alpha']
            # nmo is sometimes different than nao (because of linear
            # dependencies which get removed during SCF)
            nmo = mo_energies.shape[0]
            nocc = molecule.number_of_alpha_electrons()
            nvir = nmo - nocc
            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            eocc = mo_energies[:nocc]
            evir = mo_energies[nocc:]
            eov = eocc.reshape(-1, 1) - evir
        else:
            mo_occ = None
            mo_vir = None
            nocc = None
            nvir = None
            eov = None

        num_vecs = dist_trials.shape(1)

        n_total = num_vecs

        prep_t0 = tm.time()

        # determine subcomm size

        if self.rank == mpi_master():
            dt_and_subcomm_size = []

            for subcomm_size in range(1, self.nodes + 1):
                if self.nodes % subcomm_size != 0:
                    continue

                n_subcomms = self.nodes // subcomm_size

                # make sure that number of subcomms does not exceed number
                # of trial vectors
                if n_subcomms > n_total:
                    continue

                ave, res = divmod(n_total, n_subcomms)
                counts = [
                    ave + 1 if p < res else ave for p in range(n_subcomms)
                ]

                time_per_fock = self.serial_ratio + (
                    1 - self.serial_ratio) / subcomm_size
                dt = max(counts) * time_per_fock

                dt_and_subcomm_size.append((dt, subcomm_size))

            # find subcomm_size with smallest dt
            dt_and_subcomm_size.sort()
            subcomm_size = dt_and_subcomm_size[0][1]
        else:
            subcomm_size = None
        subcomm_size = self.comm.bcast(subcomm_size, root=mpi_master())

        # create subcomms

        grps = [p // subcomm_size for p in range(self.nodes)]
        subcomm = SubCommunicators(self.comm, grps)
        local_comm = subcomm.local_comm
        cross_comm = subcomm.cross_comm

        is_local_master = (local_comm.Get_rank() == mpi_master())

        # make data available on local master processes
        # and determine subcomm indices

        subcomm_index = None
        local_master_ranks = None

        if is_local_master:
            mo_occ = cross_comm.bcast(mo_occ, root=mpi_master())
            mo_vir = cross_comm.bcast(mo_vir, root=mpi_master())
            nocc, nvir = cross_comm.bcast((nocc, nvir), root=mpi_master())
            eov = cross_comm.bcast(eov, root=mpi_master())

            subcomm_index = cross_comm.Get_rank()
            local_master_ranks = cross_comm.allgather(self.rank)

        subcomm_index = local_comm.bcast(subcomm_index, root=mpi_master())
        local_master_ranks = local_comm.bcast(local_master_ranks,
                                              root=mpi_master())

        # batch_size corresponds to the number of subcomms
        batch_size = self.nodes // subcomm_size

        num_batches = n_total // batch_size
        if n_total % batch_size != 0:
            num_batches += 1

        if profiler is not None:
            profiler.add_timing_info('PreProc', tm.time() - prep_t0)

        vecs_sigma_data = None

        # go through batches

        if self.rank == mpi_master() and self.print_level > 1:
            batch_str = f'Processing {n_total} Fock builds'
            if batch_size > 1:
                batch_str += f' on {batch_size} subcommunicators'
            batch_str += '...'
            self.ostream.print_info(batch_str)
            self.ostream.print_blank()
            self.ostream.flush()

        for batch_ind in range(num_batches):

            batch_start = batch_size * batch_ind
            batch_end = min(batch_start + batch_size, num_vecs)

            vec_list = [None for idx in range(len(local_master_ranks))]

            prep_t0 = tm.time()

            # Here we need to collect the trial vectors (stored as distributed
            # arrays) to local master ranks of the subcommunicators, using
            # Alltoallv "transpose"

            # sendbuf data
            sendbuf = dist_trials.data[:, batch_start:batch_end].copy()

            # sendbuf counts and displacements
            sendbuf_counts = []
            for p in range(self.nodes):
                if p in local_master_ranks:
                    idx = local_master_ranks.index(p)
                    if idx + batch_start < batch_end:
                        sendbuf_counts.append(sendbuf.shape[0])
                    else:
                        sendbuf_counts.append(0)
                else:
                    sendbuf_counts.append(0)
            sendbuf_displs = [
                sum(sendbuf_counts[:p]) for p in range(self.nodes)
            ]

            # recvbuf counts and displacements
            vec_counts = self.comm.allgather(sendbuf.shape[0])
            if is_local_master:
                idx = subcomm_index
                if idx + batch_start < batch_end:
                    recvbuf_counts = list(vec_counts)
                else:
                    recvbuf_counts = [0 for x in vec_counts]
            else:
                recvbuf_counts = [0 for x in vec_counts]
            recvbuf_displs = [
                sum(recvbuf_counts[:p]) for p in range(self.nodes)
            ]
            recvbuf = np.zeros(sum(recvbuf_counts))

            # Alltoallv
            sendbuf = sendbuf.T.copy()
            self.comm.Alltoallv(
                [sendbuf, sendbuf_counts, sendbuf_displs, MPI.DOUBLE],
                [recvbuf, recvbuf_counts, recvbuf_displs, MPI.DOUBLE])

            if subcomm_index + batch_start < batch_end:
                if is_local_master:
                    vec_list[subcomm_index] = recvbuf

            self.comm.barrier()

            if profiler is not None:
                profiler.add_timing_info('PreProc', tm.time() - prep_t0)

            if subcomm_index + batch_start < batch_end:

                prep_t0 = tm.time()

                if is_local_master:
                    vec = vec_list[subcomm_index]
                    vec_ao = np.linalg.multi_dot(
                        [mo_occ, vec.reshape(nocc, nvir), mo_vir.T])
                    dens = [vec_ao]
                else:
                    dens = None

                if profiler is not None:
                    profiler.add_timing_info('PreProc', tm.time() - prep_t0)

                # create Fock matrices and contract with two-electron integrals
                fock = self._comp_lr_fock(dens, molecule, basis, eri_dict,
                                          dft_dict, pe_dict, profiler,
                                          local_comm)

                if profiler is not None:
                    # only increment FockCount on local master
                    if is_local_master:
                        profiler.add_timing_info('_FockCount_', 1)

                prep_t0 = tm.time()

                if is_local_master:
                    sigmas = np.zeros((nocc * nvir, 1))
                    # -np.linalg.multi_dot([mo_occ.T, fock[0], mo_vir])
                    # -np.linalg.multi_dot([mo_vir.T, fock[0], mo_occ]).T
                    cphf_mo = -np.linalg.multi_dot(
                        [mo_occ.T, (fock[0] + fock[0].T), mo_vir])
                    cphf_mo += eov * vec.reshape(nocc, nvir)
                    sigmas[:, 0] = cphf_mo.reshape(nocc * nvir)
                else:
                    sigmas = None

                if profiler is not None:
                    profiler.add_timing_info('PostProc', tm.time() - prep_t0)

            self.comm.barrier()

            prep_t0 = tm.time()

            # Here we need to distribute the sigma vectors (on local master
            # ranks), using Alltoallv "transpose"

            # sendbuf data
            if is_local_master:
                idx = subcomm_index
                if idx + batch_start < batch_end:
                    sendbuf_2 = sigmas
                else:
                    sendbuf_2 = np.zeros(0)
            else:
                sendbuf_2 = np.zeros(0)

            # sendbuf counts and displacements
            sendbuf_counts_2 = list(recvbuf_counts)
            sendbuf_displs_2 = list(recvbuf_displs)

            # revvdbuf, counts and displacements
            recvbuf_2 = np.zeros(sendbuf.shape)
            recvbuf_counts_2 = list(sendbuf_counts)
            recvbuf_displs_2 = list(sendbuf_displs)

            # Alltoallv
            self.comm.Alltoallv(
                [sendbuf_2, sendbuf_counts_2, sendbuf_displs_2, MPI.DOUBLE],
                [recvbuf_2, recvbuf_counts_2, recvbuf_displs_2, MPI.DOUBLE])
            recvbuf_2 = recvbuf_2.T.copy()

            local_sigma_data = recvbuf_2

            # Note: accumulate per-batch data to avoid incremental appending

            if vecs_sigma_data is None:
                vecs_sigma_data = local_sigma_data.copy()
            else:
                vecs_sigma_data = np.hstack((vecs_sigma_data, local_sigma_data))

            if profiler is not None:
                profiler.add_timing_info('PostProc', tm.time() - prep_t0)

        prep_t0 = tm.time()

        # Note: append sigma and trial vectors only once

        vecs_sigma = DistributedArray(vecs_sigma_data,
                                      self.comm,
                                      distribute=False)

        if self.dist_sigmas is None:
            self.dist_sigmas = DistributedArray(vecs_sigma.data,
                                                self.comm,
                                                distribute=False)
        else:
            self.dist_sigmas.append(vecs_sigma, axis=1)

        if self.dist_trials is None:
            self.dist_trials = DistributedArray(dist_trials.data,
                                                self.comm,
                                                distribute=False)
        else:
            self.dist_trials.append(dist_trials, axis=1)

        if profiler is not None:
            profiler.add_timing_info('AppendVec', tm.time() - prep_t0)

    def build_sigmas_single_comm(self,
                                 molecule,
                                 basis,
                                 scf_tensors,
                                 dist_trials,
                                 eri_dict,
                                 dft_dict,
                                 pe_dict,
                                 profiler=None):
        """
        Compute the matrix-vector product corresponding to the
        left-hand-side of the CPHF/CPKS equations.
        Appends sigma and trial vectors to member variable.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param dist_trials:
            Distributed array of trial vectors
        """

        prep_t0 = tm.time()

        if self.rank == mpi_master():
            mo = scf_tensors['C_alpha']
            nao = basis.get_dimension_of_basis()
            mo_energies = scf_tensors['E_alpha']
            # nmo is sometimes different than nao (because of linear
            # dependencies which get removed during SCF)
            nmo = mo_energies.shape[0]
            nocc = molecule.number_of_alpha_electrons()
            nvir = nmo - nocc
            mo_occ = mo[:, :nocc]
            mo_vir = mo[:, nocc:]
            eocc = mo_energies[:nocc]
            evir = mo_energies[nocc:]
            eov = eocc.reshape(-1, 1) - evir
        else:
            nao = None

        num_vecs = dist_trials.shape(1)

        batch_size = get_batch_size(self.batch_size, num_vecs, nao, self.comm)
        num_batches = get_number_of_batches(num_vecs, batch_size, self.comm)

        if profiler is not None:
            profiler.add_timing_info('PreProc', tm.time() - prep_t0)

        vecs_sigma_data = None

        if self.rank == mpi_master() and self.print_level > 1:
            batch_str = f'Processing {num_vecs} Fock build'
            if num_vecs > 1:
                batch_str += 's'
            batch_str += '...'
            self.ostream.print_info(batch_str)
            self.ostream.print_blank()
            self.ostream.flush()

        for batch_ind in range(num_batches):

            batch_start = batch_size * batch_ind
            batch_end = min(batch_start + batch_size, num_vecs)

            if self.rank == mpi_master():
                vec_list = []
            else:
                vec_list = None

            prep_t0 = tm.time()

            # loop over columns / trial vectors
            for col in range(batch_start, batch_end):
                vec = dist_trials.get_full_vector(col)

                if self.rank == mpi_master():
                    vec = vec.reshape(nocc, nvir)
                    vec_ao = np.linalg.multi_dot([mo_occ, vec, mo_vir.T])
                    vec_list.append(vec_ao)

            if profiler is not None:
                profiler.add_timing_info('PreProc', tm.time() - prep_t0)

            # create Fock matrices and contract with two-electron integrals
            fock = self._comp_lr_fock(vec_list, molecule, basis, eri_dict,
                                      dft_dict, pe_dict, profiler)

            if profiler is not None:
                # only increment FockCount on master rank
                if self.rank == mpi_master():
                    profiler.add_timing_info('_FockCount_', len(vec_list))

            prep_t0 = tm.time()

            # create sigma vectors
            if self.rank == mpi_master():
                sigmas = np.zeros((nocc * nvir, batch_end - batch_start))
            else:
                sigmas = None

            for col in range(batch_start, batch_end):
                vec = dist_trials.get_full_vector(col)
                ifock = col - batch_start

                if self.rank == mpi_master():
                    fock_vec = fock[ifock]
                    cphf_mo = (
                        -np.linalg.multi_dot([mo_occ.T, fock_vec, mo_vir]) -
                        np.linalg.multi_dot([mo_vir.T, fock_vec, mo_occ]).T +
                        vec.reshape(nocc, nvir) * eov)
                    sigmas[:, ifock] = cphf_mo.reshape(nocc * nvir)

            dist_sigma = DistributedArray(sigmas, self.comm)

            if vecs_sigma_data is None:
                vecs_sigma_data = dist_sigma.data.copy()
            else:
                vecs_sigma_data = np.hstack((vecs_sigma_data, dist_sigma.data))

            if profiler is not None:
                profiler.add_timing_info('PostProc', tm.time() - prep_t0)

        prep_t0 = tm.time()

        # Note: append sigma and trial vectors only once

        vecs_sigma = DistributedArray(vecs_sigma_data,
                                      self.comm,
                                      distribute=False)

        if self.dist_sigmas is None:
            self.dist_sigmas = DistributedArray(vecs_sigma.data,
                                                self.comm,
                                                distribute=False)
        else:
            self.dist_sigmas.append(vecs_sigma, axis=1)

        if self.dist_trials is None:
            self.dist_trials = DistributedArray(dist_trials.data,
                                                self.comm,
                                                distribute=False)
        else:
            self.dist_trials.append(dist_trials, axis=1)

        if profiler is not None:
            profiler.add_timing_info('AppendVec', tm.time() - prep_t0)

    def build_sigmas_single_comm_unrestricted(self,
                                              molecule,
                                              basis,
                                              scf_tensors,
                                              dist_trials,
                                              eri_dict,
                                              dft_dict,
                                              pe_dict,
                                              profiler=None):
        """
        Compute the matrix-vector product corresponding to the
        left-hand-side of the CPHF/CPKS equations.
        Appends sigma and trial vectors to member variable.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param dist_trials:
            Distributed array of trial vectors
        """

        prep_t0 = tm.time()

        if self.rank == mpi_master():
            mo_a = scf_tensors['C_alpha']
            mo_b = scf_tensors['C_beta']
            nao = basis.get_dimension_of_basis()
            orb_ene_a = scf_tensors['E_alpha']
            orb_ene_b = scf_tensors['E_beta']
            nocc_a = molecule.number_of_alpha_electrons()
            nocc_b = molecule.number_of_beta_electrons()
            mo_occ_a = mo_a[:, :nocc_a]
            mo_occ_b = mo_b[:, :nocc_b]
            mo_vir_a = mo_a[:, nocc_a:]
            mo_vir_b = mo_b[:, nocc_b:]
            eocc_a = orb_ene_a[:nocc_a]
            evir_a = orb_ene_a[nocc_a:]
            eocc_b = orb_ene_b[:nocc_b]
            evir_b = orb_ene_b[nocc_b:]
            nvir_a = evir_a.shape[0]
            nvir_b = evir_b.shape[0]
            nov_a = nocc_a * nvir_a
            eov_a = eocc_a.reshape(-1, 1) - evir_a
            eov_b = eocc_b.reshape(-1, 1) - evir_b
        else:
            nao = None

        num_vecs = dist_trials.shape(1)

        batch_size = get_batch_size(self.batch_size, num_vecs, nao, self.comm)
        num_batches = get_number_of_batches(num_vecs, batch_size, self.comm)

        if profiler is not None:
            profiler.add_timing_info('PreProc', tm.time() - prep_t0)

        vecs_sigma_data = None

        if self.rank == mpi_master() and self.print_level > 1:
            batch_str = f'Processing {num_vecs} Fock build'
            if num_vecs > 1:
                batch_str += 's'
            batch_str += '...'
            self.ostream.print_info(batch_str)
            self.ostream.print_blank()
            self.ostream.flush()

        for batch_ind in range(num_batches):

            batch_start = batch_size * batch_ind
            batch_end = min(batch_start + batch_size, num_vecs)

            if self.rank == mpi_master():
                vec_list_a = []
                vec_list_b = []
            else:
                vec_list_a = None
                vec_list_b = None

            prep_t0 = tm.time()

            # loop over columns / trial vectors
            for col in range(batch_start, batch_end):
                vec = dist_trials.get_full_vector(col)

                if self.rank == mpi_master():
                    vec_a = vec[:nov_a].reshape(nocc_a, nvir_a)
                    vec_b = vec[nov_a:].reshape(nocc_b, nvir_b)
                    vec_ao_a = np.linalg.multi_dot(
                        [mo_occ_a, vec_a, mo_vir_a.T])
                    vec_ao_b = np.linalg.multi_dot(
                        [mo_occ_b, vec_b, mo_vir_b.T])
                    vec_list_a.append(vec_ao_a)
                    vec_list_b.append(vec_ao_b)

            if profiler is not None:
                profiler.add_timing_info('PreProc', tm.time() - prep_t0)

            # create Fock matrices and contract with two-electron integrals
            fock = self._comp_lr_fock_unrestricted([vec_list_a, vec_list_b],
                                                   molecule, basis, eri_dict,
                                                   dft_dict, pe_dict, profiler)

            if profiler is not None:
                # only increment FockCount on master rank
                if self.rank == mpi_master():
                    # Note: Jab/Ka/Kb -> 3 Fock builds for openshell
                    profiler.add_timing_info('_FockCount_', len(vec_list_a) * 3)

            prep_t0 = tm.time()

            # create sigma vectors
            if self.rank == mpi_master():
                sigmas_a = np.zeros((nocc_a * nvir_a, batch_end - batch_start))
                sigmas_b = np.zeros((nocc_b * nvir_b, batch_end - batch_start))
            else:
                sigmas = None

            for col in range(batch_start, batch_end):
                vec = dist_trials.get_full_vector(col)
                ifock = col - batch_start

                if self.rank == mpi_master():
                    vec_a = vec[:nov_a].reshape(nocc_a, nvir_a)
                    vec_b = vec[nov_a:].reshape(nocc_b, nvir_b)
                    fock_vec_a = fock[ifock * 2 + 0]
                    fock_vec_b = fock[ifock * 2 + 1]
                    cphf_mo_a = (
                        -np.linalg.multi_dot([mo_occ_a.T, fock_vec_a, mo_vir_a])
                        - np.linalg.multi_dot(
                            [mo_vir_a.T, fock_vec_a, mo_occ_a]).T +
                        vec_a.reshape(nocc_a, nvir_a) * eov_a)
                    cphf_mo_b = (
                        -np.linalg.multi_dot([mo_occ_b.T, fock_vec_b, mo_vir_b])
                        - np.linalg.multi_dot(
                            [mo_vir_b.T, fock_vec_b, mo_occ_b]).T +
                        vec_b.reshape(nocc_b, nvir_b) * eov_b)
                    sigmas_a[:, ifock] = cphf_mo_a.reshape(nocc_a * nvir_a)
                    sigmas_b[:, ifock] = cphf_mo_b.reshape(nocc_b * nvir_b)
                    sigmas = np.vstack((sigmas_a, sigmas_b))

            dist_sigma = DistributedArray(sigmas, self.comm)

            if vecs_sigma_data is None:
                vecs_sigma_data = dist_sigma.data.copy()
            else:
                vecs_sigma_data = np.hstack((vecs_sigma_data, dist_sigma.data))

            if profiler is not None:
                profiler.add_timing_info('PostProc', tm.time() - prep_t0)

        prep_t0 = tm.time()

        # Note: append sigma and trial vectors only once

        vecs_sigma = DistributedArray(vecs_sigma_data,
                                      self.comm,
                                      distribute=False)

        if self.dist_sigmas is None:
            self.dist_sigmas = DistributedArray(vecs_sigma.data,
                                                self.comm,
                                                distribute=False)
        else:
            self.dist_sigmas.append(vecs_sigma, axis=1)

        if self.dist_trials is None:
            self.dist_trials = DistributedArray(dist_trials.data,
                                                self.comm,
                                                distribute=False)
        else:
            self.dist_trials.append(dist_trials, axis=1)

        if profiler is not None:
            profiler.add_timing_info('AppendVec', tm.time() - prep_t0)

    def _get_precond(self,
                     molecule,
                     basis,
                     scf_tensors,
                     method_type="restricted"):
        """
        Construct the preconditioners.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param scf_tensors:
            The SCF tensors.
        :param method_type:
            The type of method (restricted, unresticted).

        :return:
            The distributed preconditioners.
        """
        if self.rank == mpi_master():
            orb_ene_a = scf_tensors['E_alpha']
            orb_ene_b = scf_tensors['E_beta']
            nocc_a = molecule.number_of_alpha_electrons()
            nocc_b = molecule.number_of_beta_electrons()
            eocc_a = orb_ene_a[:nocc_a]
            eocc_b = orb_ene_b[:nocc_b]
            evir_a = orb_ene_a[nocc_a:]
            evir_b = orb_ene_b[nocc_b:]
            eov_a = eocc_a.reshape(-1, 1) - evir_a
            eov_b = eocc_b.reshape(-1, 1) - evir_b
        else:
            eov_a = None
            eov_b = None

        eov_a = self.comm.bcast(eov_a, root=mpi_master())
        eov_b = self.comm.bcast(eov_b, root=mpi_master())

        if method_type == "restricted":
            nocc_a = eov_a.shape[0]
            nvir_a = eov_a.shape[1]
            precond = (1.0 / eov_a).reshape(nocc_a * nvir_a)
        else:
            nocc_a = eov_a.shape[0]
            nvir_a = eov_a.shape[1]
            nocc_b = eov_b.shape[0]
            nvir_b = eov_b.shape[1]
            precond_a = (1.0 / eov_a).reshape(nocc_a * nvir_a)
            precond_b = (1.0 / eov_b).reshape(nocc_b * nvir_b)
            # stack alpha and beta
            precond = np.hstack((precond_a, precond_b))

        dist_precond = DistributedArray(precond, self.comm)

        return dist_precond

    def setup_trials(self,
                     molecule,
                     dist_precond,
                     dist_rhs,
                     old_trials=None,
                     renormalize=True):
        """
        Set up trial vectors and apply preconditioner to them.

        :param molecule:
            The molecule.
        :param dist_precond:
            The preconditioner as a distributed array.
        :param dist_rhs:
            The CPHF RHS as distributed array.
        :param old_trials:
            The previous trial vectors.
        :param renormalize:
            The flag for normalization.

        :return:
            The new set of orthonormalized trial vectors.
        """

        # apply preconditioner to all right-hand sides and create
        # distributed array for all trial vectors
        trials = []

        for k in range(len(dist_rhs)):
            v = DistributedArray(dist_precond.data * dist_rhs[k].data,
                                 self.comm,
                                 distribute=False)
            norm = np.sqrt(v.squared_norm())

            if norm > self.norm_thresh:
                trials.append(v.data[:])

        dist_trials = DistributedArray(np.array(trials).T,
                                       self.comm,
                                       distribute=False)

        if dist_trials.data.size == 0:
            dist_trials.data = np.zeros((dist_trials.shape(0), 0))

        if old_trials is not None:
            # t = t - (b (b.T t))
            bT_new_trials = old_trials.matmul_AtB_allreduce(dist_trials)
            dist_trials_proj = old_trials.matmul_AB_no_gather(bT_new_trials)
            dist_trials.data -= dist_trials_proj.data

        # remove linear dependencies and orthonormalize trial vectors
        if renormalize:
            if dist_trials.data.ndim > 0 and dist_trials.shape(0) > 0:
                dist_trials = self._remove_linear_dependence_half_size(
                    dist_trials, self.lindep_thresh)
                dist_trials = (
                    self._orthogonalize_gram_schmidt_half_size(dist_trials))
                dist_trials = self._normalize_half_size(dist_trials)

        if self.rank == mpi_master():
            assert_msg_critical(dist_trials.data.size > 0,
                                'CphfSolver: trial vectors are empty')

        return dist_trials

    def check_convergence(self, relative_residual_norm):
        """
        Checks convergence.

        :param relative_residual_norm:
            Relative residual norms.
        """

        if self.rank == mpi_master():
            max_residual = max(relative_residual_norm)
            if max_residual < self.conv_thresh:
                self._is_converged = True

        self._is_converged = self.comm.bcast(self.is_converged,
                                             root=mpi_master())

    def compute_rhs(self, molecule, basis, scf_tensors, eri_dict, dft_dict,
                    pe_dict, *args):
        """
        Computes the right hand side for the CPHF equations for
        the analytical Hessian, all atomic coordinates.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The tensors from the converged SCF calculation.

        :returns:
            The RHS of the CPHF equations.
        """

        return None

    def print_iteration(self, relative_residual_norm, molecule):
        """
        Prints information of the iteration.

        :param relative_residual_norm:
            Relative residual norms.
        :param molecule:
            The molecule.
        """
        natm = molecule.number_of_atoms()
        atom_symbols = molecule.get_labels()

        width = 92
        output_header = '*** Iteration:   {} '.format(self._cur_iter + 1)
        output_header += '* Residuals (Max,Min): '
        output_header += '{:.2e} and {:.2e}'.format(max(relative_residual_norm),
                                                    min(relative_residual_norm))
        self.ostream.print_header(output_header.ljust(width))
        self.ostream.print_blank()

        if self.print_residuals:
            coord = 'xyz'

            # TODO: change power based on the convergence threshold
            if self.conv_thresh < 1e-4:
                power = int(str(self.conv_thresh).split("-")[1]) + 1
            else:
                power = 5
            residual_norm_string = 'Residual Norm: {:.%df}' % (power)

            for i in range(natm):
                for x in range(3):
                    coord_label = '{:16s}'.format('   ' + str(i + 1) + ' ' +
                                                  atom_symbols[i] + '(' +
                                                  coord[x] + ')')
                    rel_res = relative_residual_norm[3 * i + x]
                    output_iter = residual_norm_string.format(rel_res)
                    self.ostream.print_line(coord_label + output_iter)
            self.ostream.print_blank()

        self.ostream.flush()

    def print_cphf_header(self, title):
        """
        Prints information on the solver setup
        """

        return None

    def _write_checkpoint_for_cphf(self, molecule, basis, dft_dict, pe_dict,
                                   labels):
        """
        Writes checkpoint file.

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
            self.checkpoint_file = f'{self.filename}_orbrsp.h5'
        elif (self.checkpoint_file is not None and
              self.checkpoint_file.endswith('_rsp.h5')):
            self.checkpoint_file = (self.checkpoint_file[:-len('_rsp.h5')] +
                                    '_orbrsp.h5')

        t0 = tm.time()

        if self.rank == mpi_master():
            success = write_rsp_hdf5(self.checkpoint_file, [], [], molecule,
                                     basis, dft_dict, pe_dict, self.ostream)
        else:
            success = False
        success = self.comm.bcast(success, root=mpi_master())

        if success:
            dist_arrays = [self.dist_trials, self.dist_sigmas]

            for dist_array, label in zip(dist_arrays, labels):
                dist_array.append_to_hdf5_file(self.checkpoint_file, label)

            checkpoint_text = 'Time spent in writing checkpoint file: '
            checkpoint_text += f'{(tm.time() - t0):.2f} sec'
            self.ostream.print_info(checkpoint_text)
            self.ostream.print_blank()

    def _read_checkpoint_for_cphf(self, labels):
        """
        Reads distributed arrays from checkpoint file.

        :param labels:
            The list of labels.
        """

        dist_arrays = [
            DistributedArray.read_from_hdf5_file(self.checkpoint_file, label,
                                                 self.comm) for label in labels
        ]

        self.dist_trials, self.dist_sigmas = dist_arrays

        checkpoint_text = 'Restarting from checkpoint file: '
        checkpoint_text += self.checkpoint_file
        self.ostream.print_info(checkpoint_text)
        self.ostream.print_blank()

    def _graceful_exit_for_cphf(self, molecule, basis, dft_dict, pe_dict,
                                labels):
        """
        Gracefully exits the program.

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

        self.ostream.print_blank()
        self.ostream.print_info('Preparing for a graceful termination...')
        self.ostream.print_blank()
        self.ostream.flush()

        self._write_checkpoint_for_cphf(molecule, basis, dft_dict, pe_dict,
                                        labels)

        self.ostream.print_info('...done.')
        self.ostream.print_blank()
        self.ostream.print_info('Exiting program.')
        self.ostream.print_blank()
        self.ostream.flush()

        sys.exit(0)
