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

from .veloxchemlib import mpi_master, rotatory_strength_in_cgs
from .veloxchemlib import denmat
from .aodensitymatrix import AODensityMatrix
from .profiler import Profiler
from .tdaeigensolverbase import TdaEigenSolverBase
from .blockdavidson import BlockDavidsonSolver
from .molecularorbitals import MolecularOrbitals
from .visualizationdriver import VisualizationDriver
from .cubicgrid import CubicGrid
from .sanitychecks import (molecule_sanity_check, scf_results_sanity_check,
                           ri_sanity_check, dft_sanity_check, pe_sanity_check,
                           solvation_model_sanity_check)
from .errorhandler import assert_msg_critical
from .checkpoint import read_rsp_hdf5, write_rsp_hdf5
from .resultsio import (write_lr_rsp_results_to_hdf5,
                        write_detach_attach_to_hdf5)


class TdaUnrestrictedEigenSolver(TdaEigenSolverBase):
    """
    Implements TDA excited states computation schheme for Hartree-Fock/Kohn-Sham
    level of theory.

    # vlxtag: UHF, Absorption, ECD, CIS
    # vlxtag: UKS, Absorption, ECD, TDA

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - nstates: The number of excited states to be determined.
        - core_excitation: The flag for computing core-excited states.
        - num_core_orbitals: The number of core-orbitals involved in
          the core-excited states.
        - solver: The eigenvalues solver.
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
        Initializes TDA excited states computation drived to default setup.
        """

        super().__init__(comm, ostream)

    def compute(self, molecule, basis, scf_results):
        """
        Performs TDA excited states calculation using molecular data.

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

        if self.lindep_thresh is None:
            self.lindep_thresh = self.conv_thresh * 1.0e-2

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
            self._print_header('TDA Eigensolver', nstates=self.nstates)

        # set start time

        self.start_time = tm.time()

        # prepare molecular orbitals

        if self.rank == mpi_master():
            orb_ene_a = scf_results['E_alpha']
            orb_ene_b = scf_results['E_beta']

            norb = orb_ene_a.shape[0]

            nocc_a = molecule.number_of_alpha_occupied_orbitals(basis)
            nocc_b = molecule.number_of_beta_occupied_orbitals(basis)

            # check nstates, core excitation, restricted subspace
            assert_msg_critical(
                self.nstates
                <= (nocc_a * (norb - nocc_a) + nocc_b * (norb - nocc_b)),
                f'{type(self).__name__}: too many excited states')

            if self.core_excitation:
                assert_msg_critical(
                    self.num_core_orbitals > 0,
                    f'{type(self).__name__}: num_core_orbitals not set or invalid'
                )
                assert_msg_critical(
                    (self.num_core_orbitals < nocc_a and
                     self.num_core_orbitals < nocc_b),
                    f'{type(self).__name__}: num_core_orbitals too large')

            elif getattr(self, 'restricted_subspace', False):
                assert_msg_critical(
                    False,
                    f'{type(self).__name__}: restricted_subspace not implemented'
                )

        else:
            orb_ene_a = None
            orb_ene_b = None
            norb = None
            nocc_a = None
            nocc_b = None

        norb, nocc_a, nocc_b = self.comm.bcast((norb, nocc_a, nocc_b),
                                               root=mpi_master())

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

        # TODO: enable PE
        assert_msg_critical(
            not self._pe, f'{type(self).__name__}: ' +
            'not yet implemented for polarizable embedding')

        # set up trial excitation vectors on master node

        diag_mat, trial_mat = self._gen_trial_vectors(molecule,
                                                      (orb_ene_a, orb_ene_b),
                                                      (nocc_a, nocc_b))

        # block Davidson algorithm setup

        self.solver = BlockDavidsonSolver(self.max_subspace_dim,
                                          self.collapse_nvec,
                                          self.lindep_thresh)

        # read initial guess from restart file

        if self.restart:
            if self.rank == mpi_master():
                rst_trial_mat, rst_sig_mat = read_rsp_hdf5(
                    self.checkpoint_file, ['TDA_trials', 'TDA_sigmas'],
                    molecule, basis, dft_dict, pe_dict, self.ostream)
                self.restart = (rst_trial_mat is not None and
                                rst_sig_mat is not None)
                if self.restart:
                    self.restart = self.match_settings(self.checkpoint_file)
            self.restart = self.comm.bcast(self.restart, root=mpi_master())

        if self.restart:
            if self.rank == mpi_master():
                self.solver.add_iteration_data(rst_sig_mat, rst_trial_mat,
                                               self.nstates)

            checkpoint_nstates = self._read_nstates_from_checkpoint()

            # print warning if nstates is not present in the restart file
            if checkpoint_nstates is None:
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
            elif checkpoint_nstates < self.nstates:
                self.ostream.print_info(
                    'Generating initial guesses for ' +
                    f'{self.nstates - checkpoint_nstates} more states...')
                self.ostream.print_blank()

                diag_mat_not_used, trial_mat = self._gen_trial_vectors(
                    molecule, (orb_ene_a, orb_ene_b), (nocc_a, nocc_b),
                    checkpoint_nstates)

                # project out trial vector components already in reduced space
                if self.rank == mpi_master():
                    trial_mat = self.solver.project_trial_vectors(trial_mat)
                    n_trials = trial_mat.shape[1]
                else:
                    n_trials = None
                n_trials = self.comm.bcast(n_trials, root=mpi_master())

                if n_trials > 0:
                    tdens_a = self._get_trans_densities(trial_mat,
                                                        scf_results,
                                                        molecule,
                                                        basis,
                                                        spin='alpha')
                    tdens_b = self._get_trans_densities(trial_mat,
                                                        scf_results,
                                                        molecule,
                                                        basis,
                                                        spin='beta')
                    fock = self._comp_lr_fock_unrestricted(
                        (tdens_a, tdens_b), molecule, basis, eri_dict, dft_dict,
                        pe_dict, profiler)
                    if self.rank == mpi_master():
                        sig_mat = self._get_sigmas(fock, scf_results, molecule,
                                                   basis, trial_mat)
                        self.solver.add_iteration_data(sig_mat, trial_mat,
                                                       self.nstates)

            if self.rank == mpi_master():
                trial_mat = self.solver.compute(diag_mat)

        profiler.check_memory_usage('Initial guess')

        # start TDA iteration

        for i in range(self.max_iter):

            # in case of restart, check convergence at first iteration
            if self.restart and i == 0:
                self._check_convergence()
                if self._is_converged:
                    if self.rank == mpi_master():
                        self._print_iter_data(i)
                    break

            profiler.set_timing_key(f'Iteration {i + 1}')

            # perform linear transformation of trial vectors

            tdens_a = self._get_trans_densities(trial_mat,
                                                scf_results,
                                                molecule,
                                                basis,
                                                spin='alpha')
            tdens_b = self._get_trans_densities(trial_mat,
                                                scf_results,
                                                molecule,
                                                basis,
                                                spin='beta')
            fock = self._comp_lr_fock_unrestricted((tdens_a, tdens_b), molecule,
                                                   basis, eri_dict, dft_dict,
                                                   pe_dict, profiler)

            profiler.start_timer('ReducedSpace')

            # solve eigenvalues problem on master node

            if self.rank == mpi_master():

                sig_mat = self._get_sigmas(fock, scf_results, molecule, basis,
                                           trial_mat)

                self.solver.add_iteration_data(sig_mat, trial_mat, self.nstates)

                trial_mat = self.solver.compute(diag_mat)

                self._print_iter_data(i)

            profiler.stop_timer('ReducedSpace')

            profiler.check_memory_usage(f'Iteration {i + 1}')

            profiler.print_memory_tracing(self.ostream)

            self._cur_iter = i

            # check convergence

            self._check_convergence()

            # write checkpoint file

            if self.rank == mpi_master():
                trials = self.solver.trial_matrices
                sigmas = self.solver.sigma_matrices
                write_rsp_hdf5(self.checkpoint_file, [trials, sigmas],
                               ['TDA_trials', 'TDA_sigmas'], molecule, basis,
                               dft_dict, pe_dict, self.ostream)
            self._write_settings_to_checkpoint()
            self._add_nstates_to_checkpoint()

            # finish TDA after convergence

            if self._is_converged:
                break

        # converged?
        if self.rank == mpi_master():
            self._print_convergence('{:d} excited states'.format(self.nstates))

            # final hdf5 file to save response results
            if self.filename is not None:
                final_h5_fname = f'{self.filename}.h5'
            else:
                final_h5_fname = None

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of TDA eigensolver')
        profiler.print_memory_usage(self.ostream)

        # compute 1e dipole integrals

        integrals = self._comp_onee_integrals(molecule, basis)

        # print converged excited states

        if self.rank == mpi_master() and self._is_converged:
            if self.core_excitation:
                mo_occ_a = scf_results['C_alpha'][:, :self.
                                                  num_core_orbitals].copy()
                mo_occ_b = scf_results['C_beta'][:, :self.
                                                 num_core_orbitals].copy()
            else:
                mo_occ_a = scf_results['C_alpha'][:, :nocc_a].copy()
                mo_occ_b = scf_results['C_beta'][:, :nocc_b].copy()

            mo_vir_a = scf_results['C_alpha'][:, nocc_a:].copy()
            mo_vir_b = scf_results['C_beta'][:, nocc_b:].copy()

            if self.core_excitation:
                n_ov_a = self.num_core_orbitals * mo_vir_a.shape[1]
            else:
                n_ov_a = mo_occ_a.shape[1] * mo_vir_a.shape[1]

            eigvals, rnorms = self.solver.get_eigenvalues()
            eigvecs = self.solver.ritz_vectors.copy()

            eigvecs_a = eigvecs[:n_ov_a, :]
            eigvecs_b = eigvecs[n_ov_a:, :]

            trans_dipoles_a = self._comp_trans_dipoles(integrals, eigvals,
                                                       eigvecs_a, mo_occ_a,
                                                       mo_vir_a)
            trans_dipoles_b = self._comp_trans_dipoles(integrals, eigvals,
                                                       eigvecs_b, mo_occ_b,
                                                       mo_vir_b)

            trans_dipoles = {}
            for trans_dipole_key in trans_dipoles_a:
                trans_dipoles[trans_dipole_key] = (
                    trans_dipoles_a[trans_dipole_key] +
                    trans_dipoles_b[trans_dipole_key])

            oscillator_strengths = (2.0 / 3.0) * np.sum(
                trans_dipoles['electric']**2, axis=1) * eigvals

            rotatory_strengths = np.sum(
                trans_dipoles['velocity'] * trans_dipoles['magnetic'],
                axis=1) * rotatory_strength_in_cgs()

        # natural transition orbitals and detachment/attachment densities

        nto_lambdas_a = []
        nto_lambdas_b = []
        nto_cube_files_a = []
        nto_cube_files_b = []
        dens_cube_files = []

        excitation_details = []

        for s in range(self.nstates):
            if self.rank == mpi_master() and self._is_converged:
                t_mat_a = eigvecs[:n_ov_a, s].reshape(mo_occ_a.shape[1],
                                                      mo_vir_a.shape[1])
                t_mat_b = eigvecs[n_ov_a:, s].reshape(mo_occ_b.shape[1],
                                                      mo_vir_b.shape[1])

                # save excitation details
                excitation_details.append(
                    self.get_excitation_details_unrestricted(
                        (eigvecs[:n_ov_a, s], eigvecs[n_ov_a:, s]),
                        (mo_occ_a.shape[1], mo_occ_b.shape[1]),
                        (mo_vir_a.shape[1], mo_vir_b.shape[1])))

            if self.nto or self.detach_attach:
                vis_drv = VisualizationDriver(self.comm)
                if self.cube_origin is None or self.cube_stepsize is None:
                    cubic_grid = vis_drv.gen_cubic_grid(molecule,
                                                        self.cube_points)
                else:
                    cubic_grid = CubicGrid(self.cube_origin, self.cube_stepsize,
                                           self.cube_points)

            if self.nto and self._is_converged:
                self.ostream.print_info(
                    'Running NTO analysis for S{:d}...'.format(s + 1))
                self.ostream.flush()

                if self.rank == mpi_master():
                    nto_mo = self._compute_nto_unrestricted((t_mat_a, t_mat_b),
                                                            (mo_occ_a, mo_occ_b),
                                                            (mo_vir_a, mo_vir_b))

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

                    # Add the NTO to the final hdf5 file.
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

            if self.detach_attach and self._is_converged:
                self.ostream.print_info(
                    'Running detachment/attachment analysis for S{:d}...'.
                    format(s + 1))
                self.ostream.flush()

                if self.rank == mpi_master():
                    dens_D_alpha, dens_A_alpha = self.get_detach_attach_densities(
                        t_mat_a, None, mo_occ_a, mo_vir_a)
                    dens_D_beta, dens_A_beta = self.get_detach_attach_densities(
                        t_mat_b, None, mo_occ_b, mo_vir_b)

                    dens_D = dens_D_alpha + dens_D_beta
                    dens_A = dens_A_alpha + dens_A_beta

                    # Add the detachment and attachment density matrices
                    # to the final hdf5 file
                    state_label = f'S{s + 1}'
                    if final_h5_fname is not None:
                        write_detach_attach_to_hdf5(final_h5_fname, state_label,
                                                    dens_D, dens_A)

                # Generate and save cube files
                if self.detach_attach_cubes:
                    if self.rank == mpi_master():
                        dens_DA = AODensityMatrix([dens_D, dens_A], denmat.rest)
                    else:
                        dens_DA = AODensityMatrix()
                    dens_DA = dens_DA.broadcast(self.comm, root=mpi_master())

                    dens_cube_fnames = self.write_detach_attach_cubes(
                        cubic_grid, molecule, basis, s, dens_DA)

                    if self.rank == mpi_master():
                        dens_cube_files.append(dens_cube_fnames)

        if (self.nto or self.detach_attach) and self._is_converged:
            self.ostream.print_blank()
            self.ostream.flush()

        # results

        if self.rank == mpi_master() and self._is_converged:

            ret_dict = {
                'eigenvalues': eigvals,
                'eigenvectors': eigvecs,
                'electric_transition_dipoles': trans_dipoles['electric'],
                'velocity_transition_dipoles': trans_dipoles['velocity'],
                'magnetic_transition_dipoles': trans_dipoles['magnetic'],
                'oscillator_strengths': oscillator_strengths,
                'rotatory_strengths': rotatory_strengths,
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

            self._write_final_hdf5(final_h5_fname, molecule, basis,
                                   dft_dict['dft_func_label'],
                                   pe_dict['potfile_text'], eigvecs)

            if (self.save_solutions and final_h5_fname is not None):
                # Write response results to final checkpoint file.
                # Keep the legacy rsp HDF5 layout for compatibility.
                # Eigenvectors are written separately as S1/S2/... datasets, so
                # they do not belong in this HDF5-facing payload.
                h5_ret_dict = {
                    key: value
                    for key, value in ret_dict.items()
                    if key != 'eigenvectors'
                }
                write_lr_rsp_results_to_hdf5(final_h5_fname, h5_ret_dict)

            self._print_results(ret_dict)

            return ret_dict

        elif self.rank != mpi_master() and self._is_converged:
            # non-master rank
            # TODO: return eigenvalues on non-master ranks
            # TODO: return distributed eigenvectors
            return {
                'eigenvalues': None,
                'eigenvectors': None,
            }

        else:
            # not converged
            return {}

    def _gen_trial_vectors(self, molecule, orb_ene, nocc, n_excl_states=0):
        """
        Generates set of TDA trial vectors for given number of excited states
        by selecting primitive excitations wirh lowest approximate energies
        E_ai = e_a-e_i.

        :param molecule:
            The molecule.
        :param orb_ene:
            The orbital energies.
        :param nocc:
            The number of occupied orbitals.
        :param n_excl_states:
            Number of states to exclude. Useful for generating initial guess
            for a restarting calculation that requests more states.

        :return:
            tuple (approximate diagonal of symmetric A, set of trial vectors).
        """

        if self.rank == mpi_master():
            orb_ene_a, orb_ene_b = orb_ene
            nocc_a, nocc_b = nocc

            norb = orb_ene_a.shape[0]

            ea = orb_ene_a
            eb = orb_ene_b

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

            diag_mat = np.hstack((np.array(excitation_energies_a),
                                  np.array(excitation_energies_b)))

            trial_mat = np.zeros((n_exc_a + n_exc_b, 0))

            # number of excitations to be excluded from initial guess
            guess_excl_nstates = self._get_initial_guess_size_for_excitations(
                n_excl_states)

            # total number of excitations in initial guess
            guess_nstates = self._get_initial_guess_size_for_excitations(
                self.nstates)

            for i, a in sorted(w_a,
                               key=w_a.get)[guess_excl_nstates:guess_nstates]:
                ia = excitations_a.index((i, a))
                trial_mat_a = np.zeros((n_exc_a + n_exc_b, 1))
                trial_mat_a[ia, 0] = 1.0
                trial_mat = np.hstack((trial_mat, trial_mat_a))

            for j, b in sorted(w_b,
                               key=w_b.get)[guess_excl_nstates:guess_nstates]:
                jb = excitations_b.index((j, b))
                trial_mat_b = np.zeros((n_exc_a + n_exc_b, 1))
                trial_mat_b[jb + n_exc_a, 0] = 1.0
                trial_mat = np.hstack((trial_mat, trial_mat_b))

            return diag_mat, trial_mat

        return None, None

    def _get_trans_densities(self,
                             trial_mat,
                             tensors,
                             molecule,
                             basis,
                             spin='alpha'):
        """
        Computes the transition densities.

        :param trial_mat:
            The matrix containing the Z vectors as columns.
        :param tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param spin:
            The spin.

        :return:
            The transition density matrix.
        """

        # form transition densities

        if self.rank == mpi_master():
            if spin == 'alpha':
                nocc = molecule.number_of_alpha_occupied_orbitals(basis)
                mo_key = 'C_alpha'
            elif spin == 'beta':
                nocc = molecule.number_of_beta_occupied_orbitals(basis)
                mo_key = 'C_beta'

            norb = tensors[mo_key].shape[1]
            nvir = norb - nocc

            if self.core_excitation:
                mo_occ = tensors[mo_key][:, :self.num_core_orbitals].copy()
            else:
                mo_occ = tensors[mo_key][:, :nocc].copy()
            mo_vir = tensors[mo_key][:, nocc:].copy()

            nocc_a = molecule.number_of_alpha_occupied_orbitals(basis)
            if self.core_excitation:
                n_ov_a = self.num_core_orbitals * (norb - nocc_a)
            else:
                n_ov_a = nocc_a * (norb - nocc_a)
            nocc_a = None

            tdens = []
            for k in range(trial_mat.shape[1]):
                if self.core_excitation:
                    if spin == 'alpha':
                        mat = trial_mat[:n_ov_a,
                                        k].reshape(self.num_core_orbitals, nvir)
                    elif spin == 'beta':
                        mat = trial_mat[n_ov_a:,
                                        k].reshape(self.num_core_orbitals, nvir)
                else:
                    if spin == 'alpha':
                        mat = trial_mat[:n_ov_a, k].reshape(nocc, nvir)
                    elif spin == 'beta':
                        mat = trial_mat[n_ov_a:, k].reshape(nocc, nvir)
                mat = np.matmul(mo_occ, np.matmul(mat, mo_vir.T))
                tdens.append(mat)
        else:
            tdens = None

        return tdens

    def _get_sigmas(self, fock, tensors, molecule, basis, trial_mat):
        """
        Computes the sigma vectors.

        :param fock:
            The Fock matrix.
        :param tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param trial_mat:
            The trial vectors as 2D Numpy array.

        :return:
            The sigma vectors as 2D Numpy array.
        """

        nocc_a = molecule.number_of_alpha_occupied_orbitals(basis)
        nocc_b = molecule.number_of_beta_occupied_orbitals(basis)

        norb = tensors['C_alpha'].shape[1]

        nvir_a = norb - nocc_a
        nvir_b = norb - nocc_b

        if self.core_excitation:
            mo_occ_a = tensors['C_alpha'][:, :self.num_core_orbitals].copy()
            mo_occ_b = tensors['C_beta'][:, :self.num_core_orbitals].copy()
        else:
            mo_occ_a = tensors['C_alpha'][:, :nocc_a].copy()
            mo_occ_b = tensors['C_beta'][:, :nocc_b].copy()

        mo_vir_a = tensors['C_alpha'][:, nocc_a:].copy()
        mo_vir_b = tensors['C_beta'][:, nocc_b:].copy()

        orb_ene_a = tensors['E_alpha']
        orb_ene_b = tensors['E_beta']

        sigma_vecs = []
        for idx in range(len(fock) // 2):
            # 2e contribution
            mat_a = fock[idx * 2 + 0]
            mat_b = fock[idx * 2 + 1]
            mat_a = np.matmul(mo_occ_a.T, np.matmul(mat_a, mo_vir_a))
            mat_b = np.matmul(mo_occ_b.T, np.matmul(mat_b, mo_vir_b))
            # 1e contribution
            if self.core_excitation:
                n_ov_a = self.num_core_orbitals * nvir_a

                cjb_alpha = trial_mat[:n_ov_a,
                                      idx].reshape(self.num_core_orbitals,
                                                   nvir_a)
                mat_a += np.matmul(cjb_alpha, np.diag(orb_ene_a[nocc_a:]).T)
                mat_a -= np.matmul(np.diag(orb_ene_a[:self.num_core_orbitals]),
                                   cjb_alpha)

                cjb_beta = trial_mat[n_ov_a:,
                                     idx].reshape(self.num_core_orbitals,
                                                  nvir_b)
                mat_b += np.matmul(cjb_beta, np.diag(orb_ene_b[nocc_b:]).T)
                mat_b -= np.matmul(np.diag(orb_ene_b[:self.num_core_orbitals]),
                                   cjb_beta)

                mat = np.vstack((mat_a.reshape(-1, 1), mat_b.reshape(-1, 1)))
                sigma_vecs.append(mat)
            else:
                n_ov_a = nocc_a * nvir_a

                cjb_alpha = trial_mat[:n_ov_a, idx].reshape(nocc_a, nvir_a)
                mat_a += np.matmul(cjb_alpha, np.diag(orb_ene_a[nocc_a:]).T)
                mat_a -= np.matmul(np.diag(orb_ene_a[:nocc_a]), cjb_alpha)

                cjb_beta = trial_mat[n_ov_a:, idx].reshape(nocc_b, nvir_b)
                mat_b += np.matmul(cjb_beta, np.diag(orb_ene_b[nocc_b:]).T)
                mat_b -= np.matmul(np.diag(orb_ene_b[:nocc_b]), cjb_beta)

                mat = np.vstack((mat_a.reshape(-1, 1), mat_b.reshape(-1, 1)))
                sigma_vecs.append(mat)

        sigma_mat = sigma_vecs[0]
        for vec in sigma_vecs[1:]:
            sigma_mat = np.hstack((sigma_mat, vec))

        return sigma_mat

    def _comp_trans_dipoles(self, integrals, eigvals, eigvecs, mo_occ, mo_vir):
        """
        Computes transition dipole moments.

        :param integrals:
            The one-electron integrals.
        :param eigvals:
            The eigenvalues.
        :param eigvecs:
            The CI vectors.
        :param mo_occ:
            The occupied MO coefficients.
        :param mo_vir:
            The virtual MO coefficients.

        :return:
            The transition dipole moments.
        """

        transition_dipoles = {
            'electric': np.zeros((self.nstates, 3)),
            'velocity': np.zeros((self.nstates, 3)),
            'magnetic': np.zeros((self.nstates, 3))
        }

        for s in range(self.nstates):
            exc_vec = eigvecs[:, s].reshape(mo_occ.shape[1],
                                            mo_vir.shape[1]).copy()

            # for unrestricted
            trans_dens = np.linalg.multi_dot([mo_occ, exc_vec, mo_vir.T])

            transition_dipoles['electric'][s, :] = np.array([
                np.vdot(trans_dens, integrals['electric_dipole'][d].T)
                for d in range(3)
            ])

            transition_dipoles['velocity'][s, :] = np.array([
                np.vdot(trans_dens, integrals['linear_momentum'][d].T)
                for d in range(3)
            ]) / (-eigvals[s])

            transition_dipoles['magnetic'][s, :] = np.array([
                np.vdot(trans_dens, integrals['magnetic_dipole'][d].T)
                for d in range(3)
            ])

        return transition_dipoles

    def get_nto(self, molecule, basis, scf_results, rsp_results, state_label):
        """
        Computes natural transition orbitals for a given excited state.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_results:
            The dictionary of tensors from converged SCF wavefunction.
        :param rsp_results:
            The dictionary of results from a linear response calculation.
        :param state_label:
            A string label identifying the state, e.g. 'S1', 'S2'.

        :return:
            The NTOs as a MolecularOrbitals object.
        """

        # Extract 1-based state index from label (e.g. 'S1' -> 1)
        state_index = int(state_label[1:])

        # Extract eigenvector from rsp_results
        if 'eigenvectors' in rsp_results:
            if self.rank == mpi_master():
                eigvec = rsp_results['eigenvectors'][:, state_index - 1].copy()
            else:
                eigvec = None
        else:
            if self.rank == mpi_master():
                # for rsp_results read from h5 file
                assert_msg_critical(
                    state_label in rsp_results,
                    f'{type(self).__name__}: No eigenvector found for {state_label}')
                eigvec = rsp_results[state_label].copy()
            else:
                eigvec = None

        if self.rank == mpi_master():
            # Get MO coefficients and orbital information
            mo_a = scf_results['C_alpha']
            nocc_a = molecule.number_of_alpha_occupied_orbitals(basis)
            norb = mo_a.shape[1]
            nvir_a = norb - nocc_a
            mo_occ_a = mo_a[:, :nocc_a].copy()
            mo_vir_a = mo_a[:, nocc_a:].copy()

            mo_b = scf_results['C_beta']
            nocc_b = molecule.number_of_beta_occupied_orbitals(basis)
            nvir_b = norb - nocc_b
            mo_occ_b = mo_b[:, :nocc_b].copy()
            mo_vir_b = mo_b[:, nocc_b:].copy()

            # Split eigenvector into alpha and beta components
            n_ov_a = nocc_a * nvir_a
            eigvec_a = eigvec[:n_ov_a]
            eigvec_b = eigvec[n_ov_a:]

            # Build transition density matrices in MO basis (TDA)
            t_mat_a = eigvec_a.reshape(nocc_a, nvir_a)
            t_mat_b = eigvec_b.reshape(nocc_b, nvir_b)

            nto_mo = self._compute_nto_unrestricted(
                (t_mat_a, t_mat_b),
                (mo_occ_a, mo_occ_b),
                (mo_vir_a, mo_vir_b))
        else:
            nto_mo = MolecularOrbitals()
        nto_mo = nto_mo.broadcast(self.comm, root=mpi_master())

        return nto_mo
