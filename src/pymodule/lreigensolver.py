#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
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
from copy import deepcopy
import numpy as np
import time as tm
import sys

# TODO import electric dipole driver
# TODO import rotatory_strength_in_cgs, hartree_in_inverse_nm, fine_structure_constant, extinction_coefficient_from_beta
# TODO import VisualizationDriver and CubicGrid
from .veloxchemlib import XCFunctional, MolecularGrid
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import denmat
from .veloxchemlib import mpi_master, hartree_in_ev
from .outputstream import OutputStream
from .profiler import Profiler
from .distributedarray import DistributedArray
from .signalhandler import SignalHandler
from .linearsolver import LinearSolver
from .molecularorbitals import MolecularOrbitals
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           dft_sanity_check)
from .errorhandler import assert_msg_critical
from .inputparser import get_random_string_parallel
from .checkpoint import check_rsp_hdf5, create_hdf5, write_rsp_solution


class LinearResponseEigenSolver(LinearSolver):
    """
    Implements linear response eigensolver.

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

        self.nto = True
        self.nto_pairs = None
        self.nto_cubes = False
        self.detach_attach = False
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
                'LinearResponseEigenSolver: cube origin needs 3 numbers')

        if self.cube_stepsize is not None:
            assert_msg_critical(
                len(self.cube_stepsize) == 3,
                'LinearResponseEigenSolver: cube stepsize needs 3 numbers')

        if self.cube_points is not None:
            assert_msg_critical(
                len(self.cube_points) == 3,
                'LinearResponseEigenSolver: cube points needs 3 integers')

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

        # check dft setup
        dft_sanity_check(self, 'compute')

        # check pe setup
        self._pe_sanity_check()

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

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'LinearResponseEigenSolver: not implemented for unrestricted case')

        if self.rank == mpi_master():
            orb_ene = scf_tensors['E_alpha']
        else:
            orb_ene = None
        orb_ene = self.comm.bcast(orb_ene, root=mpi_master())
        norb = orb_ene.shape[0]
        nocc = molecule.number_of_alpha_electrons()

        if self.rank == mpi_master():
            assert_msg_critical(
                self.nstates <= nocc * (norb - nocc),
                'LinearResponseEigenSolver: too many excited states')

        # ERI information
        eri_dict = self._init_eri(molecule, basis)

        # DFT information
        dft_dict = self._init_dft(molecule, scf_tensors)

        # PE information
        pe_dict = self._init_pe(molecule, basis)

        if self.nonlinear:
            rsp_vector_labels = [
                'LR_eigen_bger_half_size',
                'LR_eigen_bung_half_size',
                'LR_eigen_e2bger_half_size',
                'LR_eigen_e2bung_half_size',
                'LR_eigen_fock_ger',
                'LR_eigen_fock_ung',
            ]
        else:
            rsp_vector_labels = [
                'LR_eigen_bger_half_size',
                'LR_eigen_bung_half_size',
                'LR_eigen_e2bger_half_size',
                'LR_eigen_e2bung_half_size',
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
            igs = self._initial_excitations(self.nstates, orb_ene, nocc, norb)
            bger, bung = self._setup_trials(igs, None)

            self._e2n_half_size(bger, bung, molecule, basis, scf_tensors,
                                eri_dict, dft_dict, pe_dict)

        profiler.check_memory_usage('Initial guess')

        exc_energies = {}
        exc_focks = {}
        exc_solutions = {}
        exc_residuals = {}
        relative_residual_norm = {}

        signal_handler = SignalHandler()
        signal_handler.add_sigterm_function(self._graceful_exit, molecule,
                                            basis, dft_dict, pe_dict,
                                            rsp_vector_labels)

        iter_per_trial_in_hours = None

        # start iterations
        for iteration in range(self.max_iter):

            iter_start_time = tm.time()

            profiler.set_timing_key(f'Iteration {iteration+1}')

            profiler.start_timer('ReducedSpace')

            self._cur_iter = iteration

            e2gg = self._dist_bger.matmul_AtB(self._dist_e2bger, 2.0)
            e2uu = self._dist_bung.matmul_AtB(self._dist_e2bung, 2.0)
            s2ug = self._dist_bung.matmul_AtB(self._dist_bger, 2.0)

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
                    fock_realger = self._dist_fock_ger.matmul_AB_no_gather(
                        c_ger[:, k])
                    fock_realung = self._dist_fock_ung.matmul_AB_no_gather(
                        c_ung[:, k])

                    fock_full_data = (fock_realger.data + fock_realung.data)

                    exc_focks[k] = DistributedArray(fock_full_data,
                                                    self.comm,
                                                    distribute=False)

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

                r_norms_2 = 2.0 * r.squared_norm(axis=0)
                x_norms_2 = 2.0 * x.squared_norm(axis=0)

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
                k: self._get_precond(orb_ene, nocc, norb, w)
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

            self._e2n_half_size(new_trials_ger, new_trials_ung, molecule, basis,
                                scf_tensors, eri_dict, dft_dict, pe_dict,
                                profiler)

            iter_in_hours = (tm.time() - iter_start_time) / 3600
            iter_per_trial_in_hours = iter_in_hours / n_new_trials

            profiler.check_memory_usage(
                'Iteration {:d} sigma build'.format(iteration + 1))

        signal_handler.remove_sigterm_function()

        self._write_checkpoint(molecule, basis, dft_dict, pe_dict,
                               rsp_vector_labels)

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

            # TODO: enable property calculation
            edip_grad = self.get_prop_grad('electric dipole', 'xyz', molecule,
                                           basis, scf_tensors, eri_dict['screening'])
            lmom_grad = self.get_prop_grad('linear momentum', 'xyz', molecule,
                                           basis, scf_tensors, eri_dict['screening'])
            mdip_grad = self.get_prop_grad('magnetic dipole', 'xyz', molecule,
                                           basis, scf_tensors, eri_dict['screening'])

            eigvals = np.array([exc_energies[s] for s in range(self.nstates)])

            elec_trans_dipoles = np.zeros((self.nstates, 3))
            velo_trans_dipoles = np.zeros((self.nstates, 3))
            magn_trans_dipoles = np.zeros((self.nstates, 3))

            if self.rank == mpi_master():
                # create h5 file for response solutions
                if (self.save_solutions and self.checkpoint_file is not None):
                    final_h5_fname = str(
                        Path(self.checkpoint_file).with_suffix('.solutions.h5'))
                    create_hdf5(final_h5_fname, molecule, basis,
                                dft_dict['dft_func_label'],
                                pe_dict['potfile_text'])

            nto_lambdas = []
            nto_h5_files = []
            nto_cube_files = []
            dens_cube_files = []

            excitation_details = []

            for s in range(self.nstates):
                eigvec = self.get_full_solution_vector(exc_solutions[s])

                if self.rank == mpi_master():
                    if self.core_excitation:
                        mo_occ = scf_tensors['C_alpha'][:, :self.
                                                        num_core_orbitals]
                        mo_vir = scf_tensors['C_alpha'][:, nocc:]
                        z_mat = eigvec[:eigvec.size // 2].reshape(
                            self.num_core_orbitals, -1)
                        y_mat = eigvec[eigvec.size // 2:].reshape(
                            self.num_core_orbitals, -1)
                    else:
                        mo_occ = scf_tensors['C_alpha'][:, :nocc]
                        mo_vir = scf_tensors['C_alpha'][:, nocc:]
                        z_mat = eigvec[:eigvec.size // 2].reshape(nocc, -1)
                        y_mat = eigvec[eigvec.size // 2:].reshape(nocc, -1)

                """
                if self.nto or self.detach_attach:
                    vis_drv = VisualizationDriver(self.comm)
                    if self.cube_origin is None or self.cube_stepsize is None:
                        cubic_grid = vis_drv.gen_cubic_grid(
                            molecule, self.cube_points)
                    else:
                        cubic_grid = CubicGrid(self.cube_origin,
                                               self.cube_stepsize,
                                               self.cube_points)
                """

                if self.nto:
                    self.ostream.print_info(
                        'Running NTO analysis for S{:d}...'.format(s + 1))
                    self.ostream.flush()

                    if self.filename is not None:
                        base_fname = self.filename
                    else:
                        name_string = get_random_string_parallel(self.comm)
                        base_fname = 'vlx_' + name_string

                    if self.rank == mpi_master():
                        nto_mo = self.get_nto(z_mat - y_mat, mo_occ, mo_vir)

                        nto_lam = nto_mo.occa_to_numpy()
                        lam_start = mo_occ.shape[1]
                        lam_end = lam_start + min(mo_occ.shape[1],
                                                  mo_vir.shape[1])
                        nto_lambdas.append(nto_lam[lam_start:lam_end])

                        nto_h5_fname = f'{base_fname}_S{s+1}_NTO.h5'
                        nto_mo.write_hdf5(nto_h5_fname)
                        nto_h5_files.append(nto_h5_fname)
                    """
                    else:
                        nto_mo = MolecularOrbitals()
                    nto_mo.broadcast(self.rank, self.comm)

                    if self.nto_cubes:
                        lam_diag, nto_cube_fnames = self.write_nto_cubes(
                            cubic_grid, molecule, basis, s, nto_mo,
                            self.nto_pairs)

                        if self.rank == mpi_master():
                            nto_cube_files.append(nto_cube_fnames)
                    """

                """
                if self.detach_attach:
                    self.ostream.print_info(
                        'Running detachment/attachment analysis for S{:d}...'.
                        format(s + 1))
                    self.ostream.flush()

                    if self.rank == mpi_master():
                        dens_D, dens_A = self.get_detach_attach_densities(
                            z_mat, y_mat, mo_occ, mo_vir)
                        dens_DA = AODensityMatrix([dens_D, dens_A], denmat.rest)
                    else:
                        dens_DA = AODensityMatrix()
                    dens_DA.broadcast(self.rank, self.comm)

                    dens_cube_fnames = self.write_detach_attach_cubes(
                        cubic_grid, molecule, basis, s, dens_DA)

                    if self.rank == mpi_master():
                        dens_cube_files.append(dens_cube_fnames)

                if self.esa:
                    dipole_drv = ElectricDipoleIntegralsDriver(self.comm)
                    dipole_matrices = dipole_drv.compute(molecule, basis)

                    if self.rank == mpi_master():
                        dipole_integrals = (dipole_matrices.x_to_numpy(),
                                            dipole_matrices.y_to_numpy(),
                                            dipole_matrices.z_to_numpy())
                    else:
                        dipole_integrals = None

                    if self.esa_from_state is None:
                        source_states = list(range(self.nstates))
                    else:
                        source_states = [self.esa_from_state - 1]

                    esa_pairs = [(s_1, s_2)
                                 for s_1 in source_states
                                 for s_2 in range(s_1 + 1, self.nstates)]

                    if self.rank == mpi_master():
                        esa_results = []
                    else:
                        esa_results = None

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
                                -1.0 *
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
                """

                if self.rank == mpi_master():
                    for ind, comp in enumerate('xyz'):
                        elec_trans_dipoles[s, ind] = np.vdot(
                            edip_grad[ind], eigvec) * (-1.0)
                        velo_trans_dipoles[s, ind] = np.vdot(
                            lmom_grad[ind], eigvec) / eigvals[s]
                        magn_trans_dipoles[s, ind] = np.vdot(
                            mdip_grad[ind], eigvec)

                    # write to h5 file for response solutions
                    if (self.save_solutions and
                            self.checkpoint_file is not None):
                        write_rsp_solution(final_h5_fname,
                                           'S{:d}'.format(s + 1), eigvec)

                    # save excitation details
                    excitation_details.append(
                        self.get_excitation_details(eigvec, mo_occ.shape[1],
                                                    mo_vir.shape[1]))

            if self.nto or self.detach_attach:
                self.ostream.print_blank()
                self.ostream.flush()

            if not self.nonlinear:

                if self.rank == mpi_master():

                    # TODO: enable property calculation
                    osc = (2.0 / 3.0) * np.sum(elec_trans_dipoles**2,
                                               axis=1) * eigvals
                    # TODO: fix rotatory_strength_in_cgs constant
                    rotatory_strength_in_cgs = 471.443648175
                    rot_vel = np.sum(velo_trans_dipoles * magn_trans_dipoles,
                                     axis=1) * rotatory_strength_in_cgs

                    ret_dict = {
                        'eigenvalues': eigvals,
                        'eigenvectors_distributed': exc_solutions,
                        'electric_transition_dipoles': elec_trans_dipoles,
                        'velocity_transition_dipoles': velo_trans_dipoles,
                        'magnetic_transition_dipoles': magn_trans_dipoles,
                        'oscillator_strengths': osc,
                        'rotatory_strengths': rot_vel,
                        'excitation_details': excitation_details,
                    }

                    if self.nto:
                        ret_dict['nto_lambdas'] = nto_lambdas
                        ret_dict['nto_h5_files'] = nto_h5_files
                        """
                        if self.nto_cubes:
                            ret_dict['nto_cubes'] = nto_cube_files
                        """

                    """
                    if self.detach_attach:
                        ret_dict['density_cubes'] = dens_cube_files

                    if self.esa:
                        ret_dict['esa_results'] = esa_results
                    """

                    if (self.save_solutions and
                            self.checkpoint_file is not None):
                        checkpoint_text = 'Response solution vectors written to file: '
                        checkpoint_text += final_h5_fname
                        self.ostream.print_info(checkpoint_text)
                        self.ostream.print_blank()

                    self._print_results(ret_dict)

                    return ret_dict
                else:
                    return {
                        'eigenvalues': eigvals,
                        'eigenvectors_distributed': exc_solutions,
                    }

            else:
                if self.rank != mpi_master():
                    elec_trans_dipoles = None
                    excitation_details = None

                elec_trans_dipoles = self.comm.bcast(elec_trans_dipoles,
                                                     root=mpi_master())
                excitation_details = self.comm.bcast(excitation_details,
                                                     root=mpi_master())

                osc = (2.0 / 3.0) * np.sum(elec_trans_dipoles**2,
                                           axis=1) * eigvals

                ret_dict = {
                    'eigenvalues': eigvals,
                    'eigenvectors_distributed': exc_solutions,
                    'focks': exc_focks,
                    'oscillator_strengths': osc,
                    'excitation_details': excitation_details,
                    'electric_transition_dipoles': elec_trans_dipoles,
                }

                return ret_dict

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

    def _initial_excitations(self, nstates, ea, nocc, norb):
        """
        Gets initial guess for excitations.

        :param nstates:
            Number of excited states.
        :param ea:
            Orbital energies.
        :param nocc:
            Number of occupied orbitals.
        :param norb:
            Number of orbitals.

        :return:
            A list of initial excitations (excitation energy and distributed
            vector).
        """

        if self.core_excitation:
            excitations = [(i, a)
                           for i in range(self.num_core_orbitals)
                           for a in range(nocc, norb)]
        else:
            excitations = [
                (i, a) for i in range(nocc) for a in range(nocc, norb)
            ]

        excitation_energies = [ea[a] - ea[i] for i, a in excitations]

        w = {ia: w for ia, w in zip(excitations, excitation_energies)}
        n_exc = len(excitations)

        final = {}
        for k, (i, a) in enumerate(sorted(w, key=w.get)[:nstates]):
            if self.rank == mpi_master():
                ia = excitations.index((i, a))

                Xn = np.zeros(2 * n_exc)
                Xn[ia] = 1.0

                Xn_T = np.zeros(2 * n_exc)
                Xn_T[:n_exc] = Xn[n_exc:]
                Xn_T[n_exc:] = Xn[:n_exc]

                Xn_ger = 0.5 * (Xn + Xn_T)[:n_exc]
                Xn_ung = 0.5 * (Xn - Xn_T)[:n_exc]

                X = np.hstack((
                    Xn_ger.reshape(-1, 1),
                    Xn_ung.reshape(-1, 1),
                ))
            else:
                X = None

            final[k] = DistributedArray(X, self.comm)

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

        # spawning needed components

        if self.core_excitation:
            ediag, sdiag = self.construct_ediag_sdiag_half(
                orb_ene, nocc, norb, self.num_core_orbitals)
        else:
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

    def get_e2(self, molecule, basis, scf_tensors):
        """
        Calculates the E[2] matrix.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            The E[2] matrix as numpy array.
        """

        if self.norm_thresh is None:
            self.norm_thresh = self.conv_thresh * 1.0e-6
        if self.lindep_thresh is None:
            self.lindep_thresh = self.conv_thresh * 1.0e-2

        self._dist_bger = None
        self._dist_bung = None
        self._dist_e2bger = None
        self._dist_e2bung = None

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'LinearResponseEigenSolver: not implemented for unrestricted case')

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

        # generate initial guess from scratch

        igs = {}
        n_exc = nocc * (norb - nocc)

        for i in range(2 * n_exc):
            Xn = np.zeros(2 * n_exc)
            Xn[i] = 1.0

            Xn_T = np.zeros(2 * n_exc)
            Xn_T[:n_exc] = Xn[n_exc:]
            Xn_T[n_exc:] = Xn[:n_exc]

            Xn_ger = 0.5 * (Xn + Xn_T)[:n_exc]
            Xn_ung = 0.5 * (Xn - Xn_T)[:n_exc]

            X = np.hstack((
                Xn_ger.reshape(-1, 1),
                Xn_ung.reshape(-1, 1),
            ))

            igs[i] = DistributedArray(X, self.comm)

        bger, bung = self._setup_trials(igs, precond=None, renormalize=False)

        self._e2n_half_size(bger, bung, molecule, basis, scf_tensors, eri_dict,
                            dft_dict, pe_dict)

        if self.rank == mpi_master():
            E2 = np.zeros((2 * n_exc, 2 * n_exc))

        for i in range(2 * n_exc):
            e2b_data = np.hstack((
                self._dist_e2bger.data[:, i:i + 1],
                self._dist_e2bung.data[:, i:i + 1],
            ))

            e2b = DistributedArray(e2b_data, self.comm, distribute=False)

            sigma = self.get_full_solution_vector(e2b)

            if self.rank == mpi_master():
                E2[:, i] = sigma[:]

        if self.rank == mpi_master():
            return E2
        else:
            return None

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

        # TODO: print more results
        """
        if self.esa:
            self._print_excited_state_absorption('Excited state absorption:',
                                                 results)
        """

    def __deepcopy__(self, memo):
        """
        Implements deepcopy.

        :param memo:
            The memo dictionary for deepcopy.

        :return:
            A deepcopy of self.
        """

        new_rsp_drv = LinearResponseEigenSolver(self.comm, self.ostream)

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
            'LinearResponseEigenSolver.get_absorption_spectrum: ' +
            'x_data should be au, ev or nm')

        assert_msg_critical(
            b_unit.lower() in ['au', 'ev'],
            'LinearResponseEigenSolver.get_absorption_spectrum: ' +
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
            'LinearResponseEigenSolver.get_ecd_spectrum: ' +
            'x_data should be au, ev or nm')

        assert_msg_critical(
            b_unit.lower() in ['au', 'ev'],
            'LinearResponseEigenSolver.get_ecd_spectrum: ' +
            'broadening parameter should be au or ev')

        au2ev = hartree_in_ev()
        auxnm = 1.0 / hartree_in_inverse_nm()

        exc_ene_au = rsp_results['eigenvalues']
        # TODO: fix rotatory_strength_in_cgs constant
        rotatory_strength_in_cgs = 471.443648175
        rot_str_au = rsp_results[
            'rotatory_strengths'] / rotatory_strength_in_cgs

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
