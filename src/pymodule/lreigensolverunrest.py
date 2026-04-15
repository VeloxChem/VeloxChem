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

from .oneeints import compute_electric_dipole_integrals
from .veloxchemlib import mpi_master, rotatory_strength_in_cgs
from .veloxchemlib import denmat
from .aodensitymatrix import AODensityMatrix
from .profiler import Profiler
from .distributedarray import DistributedArray
from .lreigensolverbase import LinearResponseEigenSolverBase
from .molecularorbitals import MolecularOrbitals
from .visualizationdriver import VisualizationDriver
from .cubicgrid import CubicGrid
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           ri_sanity_check, dft_sanity_check, pe_sanity_check,
                           solvation_model_sanity_check)
from .errorhandler import assert_msg_critical
from .mathutils import screened_eigh, symmetric_matrix_function
from .checkpoint import check_rsp_hdf5
from .resultsio import (write_lr_rsp_results_to_hdf5,
                        write_detach_attach_to_hdf5, write_rsp_solution)


class LinearResponseUnrestrictedEigenSolver(LinearResponseEigenSolverBase):
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

        super().__init__(comm, ostream)

    def compute(self, molecule, basis, scf_results):
        """
        Performs linear response calculation for a molecule and a basis set.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_results:
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

        # check solvation
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
            self._print_header('Linear Response EigenSolver',
                               nstates=self.nstates)

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

        # check number of excited states, core excitation, restricted subspace
        assert_msg_critical(
            self.nstates
            <= (nocc_a * (norb - nocc_a) + nocc_b * (norb - nocc_b)),
            f'{type(self).__name__}: too many excited states')

        if self.core_excitation:
            assert_msg_critical(
                self.num_core_orbitals > 0,
                f'{type(self).__name__}: num_core_orbitals not set or invalid')
            assert_msg_critical(
                (self.num_core_orbitals < nocc_a and
                 self.num_core_orbitals < nocc_b),
                f'{type(self).__name__}: num_core_orbitals too large')

        elif getattr(self, 'restricted_subspace', False):
            assert_msg_critical(
                False,
                f'{type(self).__name__}: restricted_subspace not implemented')

        self._check_mpi_oversubscription(
            self._get_excitation_space_dimension_unrestricted(
                nocc_a, nocc_b, norb), 'excitation space')

        # ERI information
        eri_dict = self._init_eri(molecule, basis)

        # DFT information
        dft_dict = self._init_dft(molecule, scf_results)

        # PE information
        pe_dict = self._init_pe(molecule, basis)

        # CPCM_information
        self._init_cpcm(molecule, basis)

        # For now, 'nonlinear' is not supported for unrestricted case.
        assert_msg_critical(
            not self.nonlinear,
            f'{type(self).__name__}: ' + 'not implemented for nonlinear')

        # TODO: enable PE
        assert_msg_critical(
            not self._pe, f'{type(self).__name__}: ' +
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
                if self.restart:
                    self.restart = self.match_settings(self.checkpoint_file)
            self.restart = self.comm.bcast(self.restart, root=mpi_master())

        # read initial guess from restart file
        if self.restart:
            self._read_checkpoint(rsp_vector_labels)

            checkpoint_nstates = self._read_nstates_from_checkpoint()

            # Note that we only handle the restart with additional nstates when
            # `self.nonlinear` is inactive

            # print warning if nstates is not present in the restart file
            if checkpoint_nstates is None and not self.nonlinear:
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
            elif checkpoint_nstates < self.nstates and not self.nonlinear:
                self.ostream.print_info(
                    'Generating initial guesses for ' +
                    f'{self.nstates - checkpoint_nstates} more states...')
                self.ostream.print_blank()

                igs = self._initial_excitations(self.nstates,
                                                (orb_ene_a, orb_ene_b),
                                                (nocc_a, nocc_b), norb,
                                                checkpoint_nstates)

                if igs:
                    bger, bung = self._setup_trials(igs, None, self._dist_bger,
                                                    self._dist_bung)

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
            igs = self._initial_excitations(self.nstates,
                                            (orb_ene_a, orb_ene_b),
                                            (nocc_a, nocc_b), norb)
            bger, bung = self._setup_trials(igs, None)

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

        exc_energies = {}
        # exc_focks = {}
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
            self.collapsed_subspace = False
            self.collapsed_from_dim = None
            self.collapsed_to_dim = None

            nroots = max(self.nstates, self._get_collapse_nvec())
            wn_all, c_ger_all, c_ung_all = self._solve_reduced_space(nroots)
            wn = wn_all[:self.nstates].copy()
            c_ger = c_ger_all[:, :self.nstates].copy()
            c_ung = c_ung_all[:, :self.nstates].copy()

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

                if self.print_level > 1:
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

            if self._should_collapse_subspace():
                self.ostream.print_info('Collapsing reduced space...')
                self.ostream.print_blank()

                self._collapse_current_subspace(c_ger_all, c_ung_all, molecule,
                                                basis, scf_results, eri_dict,
                                                dft_dict, pe_dict, profiler)

                collapse_str = 'Collapsed reduced space: {:d}->{:d}'.format(
                    self.collapsed_from_dim, self.collapsed_to_dim)
                self.ostream.print_info(collapse_str)
                self.ostream.print_blank()
                self.ostream.flush()

            profiler.start_timer('Orthonorm.')

            # update trial vectors
            precond = {
                k: self._get_precond((orb_ene_a, orb_ene_b), (nocc_a, nocc_b),
                                     norb, w) for k, w in enumerate(list(wn))
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
        self._add_nstates_to_checkpoint()

        # converged?
        if self.rank == mpi_master():
            self._print_convergence('Linear response')

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of LR eigensolver')
        profiler.print_memory_usage(self.ostream)

        # for unrestricted
        sqrt_2 = np.sqrt(2.0)
        for k in exc_solutions:
            exc_solutions[k].data /= sqrt_2

        self._dist_bger = None
        self._dist_bung = None
        self._dist_e2bger = None
        self._dist_e2bung = None

        self._dist_fock_ger = None
        self._dist_fock_ung = None

        # calculate properties
        if self.is_converged:
            edip_grad_a = self.get_prop_grad('electric dipole',
                                             'xyz',
                                             molecule,
                                             basis,
                                             scf_results,
                                             spin='alpha')
            lmom_grad_a = self.get_prop_grad('linear momentum',
                                             'xyz',
                                             molecule,
                                             basis,
                                             scf_results,
                                             spin='alpha')
            mdip_grad_a = self.get_prop_grad('magnetic dipole',
                                             'xyz',
                                             molecule,
                                             basis,
                                             scf_results,
                                             spin='alpha')

            edip_grad_b = self.get_prop_grad('electric dipole',
                                             'xyz',
                                             molecule,
                                             basis,
                                             scf_results,
                                             spin='beta')
            lmom_grad_b = self.get_prop_grad('linear momentum',
                                             'xyz',
                                             molecule,
                                             basis,
                                             scf_results,
                                             spin='beta')
            mdip_grad_b = self.get_prop_grad('magnetic dipole',
                                             'xyz',
                                             molecule,
                                             basis,
                                             scf_results,
                                             spin='beta')

            # for unrestricted
            edip_grad_a /= sqrt_2
            lmom_grad_a /= sqrt_2
            mdip_grad_a /= sqrt_2
            edip_grad_b /= sqrt_2
            lmom_grad_b /= sqrt_2
            mdip_grad_b /= sqrt_2

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
                        mo_occ_a = scf_results[
                            'C_alpha'][:, :self.num_core_orbitals].copy()
                        mo_vir_a = scf_results['C_alpha'][:, nocc_a:].copy()
                        z_mat_a = eigvec_a[:eigvec_a.size // 2].reshape(
                            self.num_core_orbitals, -1)
                        y_mat_a = eigvec_a[eigvec_a.size // 2:].reshape(
                            self.num_core_orbitals, -1)

                        mo_occ_b = scf_results[
                            'C_beta'][:, :self.num_core_orbitals].copy()
                        mo_vir_b = scf_results['C_beta'][:, nocc_b:].copy()
                        z_mat_b = eigvec_b[:eigvec_b.size // 2].reshape(
                            self.num_core_orbitals, -1)
                        y_mat_b = eigvec_b[eigvec_b.size // 2:].reshape(
                            self.num_core_orbitals, -1)
                    else:
                        mo_occ_a = scf_results['C_alpha'][:, :nocc_a].copy()
                        mo_vir_a = scf_results['C_alpha'][:, nocc_a:].copy()
                        z_mat_a = eigvec_a[:eigvec_a.size // 2].reshape(
                            nocc_a, -1)
                        y_mat_a = eigvec_a[eigvec_a.size // 2:].reshape(
                            nocc_a, -1)

                        mo_occ_b = scf_results['C_beta'][:, :nocc_b].copy()
                        mo_vir_b = scf_results['C_beta'][:, nocc_b:].copy()
                        z_mat_b = eigvec_b[:eigvec_b.size // 2].reshape(
                            nocc_b, -1)
                        y_mat_b = eigvec_b[eigvec_b.size // 2:].reshape(
                            nocc_b, -1)

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
                                              label=f'rsp/nto/{nto_label}_')
                    else:
                        nto_mo = MolecularOrbitals()
                    nto_mo = nto_mo.broadcast(self.comm, root=mpi_master())

                    if self.nto_cubes:
                        lam_diag_a, nto_cube_fnames_a = self.write_nto_cubes(
                            cubic_grid,
                            molecule,
                            basis,
                            s,
                            nto_mo,
                            self.nto_pairs,
                            nto_spin='alpha')
                        lam_diag_b, nto_cube_fnames_b = self.write_nto_cubes(
                            cubic_grid,
                            molecule,
                            basis,
                            s,
                            nto_mo,
                            self.nto_pairs,
                            nto_spin='beta')

                        if self.rank == mpi_master():
                            nto_cube_files_a.append(nto_cube_fnames_a)
                            nto_cube_files_b.append(nto_cube_fnames_b)

                if self.detach_attach:
                    self.ostream.print_info(
                        'Running detachment/attachment analysis for S{:d}...'.
                        format(s + 1))
                    self.ostream.flush()

                    if self.rank == mpi_master():
                        dens_D_alpha, dens_A_alpha = self.get_detach_attach_densities(
                            z_mat_a, y_mat_a, mo_occ_a, mo_vir_a)
                        dens_D_beta, dens_A_beta = self.get_detach_attach_densities(
                            z_mat_b, y_mat_b, mo_occ_b, mo_vir_b)

                        dens_D = dens_D_alpha + dens_D_beta
                        dens_A = dens_A_alpha + dens_A_beta

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

                # TODO: enable esa
                assert_msg_critical(
                    not self.esa, f'{type(self).__name__}: ' +
                    'not yet implemented for excited state absorption')

                if self.esa:
                    if self.rank == mpi_master():
                        esa_results = []
                    """
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
                    """

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

                        # Keep the legacy rsp HDF5 layout for compatibility.
                        # Solution vectors are written separately as S1/S2/...
                        # datasets, so the distributed in-memory vectors do not
                        # belong in this HDF5-facing payload.
                        h5_ret_dict = {
                            key: value
                            for key, value in ret_dict.items()
                            if key != 'eigenvectors_distributed'
                        }
                        write_lr_rsp_results_to_hdf5(final_h5_fname,
                                                     h5_ret_dict)

                    self._print_results(ret_dict)

                    return ret_dict
                else:
                    # non-master rank
                    return {
                        'eigenvalues': eigvals,
                        'eigenvectors_distributed': exc_solutions,
                    }

            else:
                # nonlinear
                # TODO: unrestricted for nonlinear
                return {}

        else:
            # not converged
            return {}

    def _solve_reduced_space(self, nroots):
        """
        Solves the reduced-space eigenvalue problem.

        :param nroots:
            The number of Ritz pairs to retain from the reduced problem.

        :return:
            Tuple of excitation energies and gerade/ungerade coefficients.
        """

        e2gg = self._dist_bger.matmul_AtB(self._dist_e2bger)
        e2uu = self._dist_bung.matmul_AtB(self._dist_e2bung)
        s2ug = self._dist_bung.matmul_AtB(self._dist_bger)

        if self.rank == mpi_master():

            e2uu_inv = symmetric_matrix_function(e2uu,
                                                 lambda x: 1.0 / x,
                                                 thresh=1.0e-12)
            ses = np.linalg.multi_dot([s2ug.T, e2uu_inv, s2ug])

            tmat = symmetric_matrix_function(e2gg,
                                             lambda x: 1.0 / np.sqrt(x),
                                             thresh=1.0e-12)
            ses_tilde = np.linalg.multi_dot([tmat.T, ses, tmat])

            evals, evecs = screened_eigh(ses_tilde,
                                         thresh=None,
                                         descending=True)

            nroots = min(nroots, evals.size)
            wn = 1.0 / np.sqrt(evals[:nroots])
            c_ger = np.matmul(tmat, evecs[:, :nroots].copy())
            c_ung = wn * np.linalg.multi_dot([e2uu_inv, s2ug, c_ger])

            for k in range(nroots):
                c_ger_k = c_ger[:, k].copy()
                c_ung_k = c_ung[:, k].copy()
                norm = np.sqrt(
                    np.linalg.multi_dot([c_ung_k.T, s2ug, c_ger_k]) +
                    np.linalg.multi_dot([c_ger_k.T, s2ug.T, c_ung_k]))
                c_ger[:, k] /= norm
                c_ung[:, k] /= norm
        else:
            wn = None
            c_ger = None
            c_ung = None

        wn = self.comm.bcast(wn, root=mpi_master())
        c_ger = self.comm.bcast(c_ger, root=mpi_master())
        c_ung = self.comm.bcast(c_ung, root=mpi_master())

        return wn, c_ger, c_ung

    def _collapse_current_subspace(self, c_ger, c_ung, molecule, basis,
                                   scf_results, eri_dict, dft_dict, pe_dict,
                                   profiler):
        """
        Collapses the reduced space to retained Ritz vectors and rebuilds
        associated sigma data.
        """

        keep_nvec = min(self._get_collapse_nvec(), c_ger.shape[1])
        prev_dim = self._reduced_space_size()

        new_bger = self._dist_bger.matmul_AB_no_gather(c_ger[:, :keep_nvec])
        new_bung = self._dist_bung.matmul_AB_no_gather(c_ung[:, :keep_nvec])

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

    def _initial_excitations(self,
                             nstates,
                             orb_ene,
                             nocc,
                             norb,
                             n_excl_states=0):
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

        # number of excitations to be excluded from initial guess
        guess_excl_nstates = self._get_initial_guess_size_for_excitations(
            n_excl_states)

        # total number of excitations in initial guess
        guess_nstates = self._get_initial_guess_size_for_excitations(nstates)

        for k, (i, a) in enumerate(
                sorted(w_a, key=w_a.get)[guess_excl_nstates:guess_nstates]):
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

        for k, (j, b) in enumerate(
                sorted(w_b, key=w_b.get)[guess_excl_nstates:guess_nstates]):
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
            ediag_a, sdiag_a = self.construct_ediag_sdiag_half(
                orb_ene_a, nocc_a, norb)
            ediag_b, sdiag_b = self.construct_ediag_sdiag_half(
                orb_ene_b, nocc_b, norb)

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
