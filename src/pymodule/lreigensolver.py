import numpy as np
import time as tm
import psutil
import sys

from .veloxchemlib import AODensityMatrix
from .veloxchemlib import mpi_master
from .veloxchemlib import rotatory_strength_in_cgs
from .veloxchemlib import denmat
from .profiler import Profiler
from .distributedarray import DistributedArray
from .signalhandler import SignalHandler
from .linearsolver import LinearSolver
from .molecularorbitals import MolecularOrbitals
from .errorhandler import assert_msg_critical
from .checkpoint import check_rsp_hdf5
from .checkpoint import append_rsp_solution_hdf5


class LinearResponseEigenSolver(LinearSolver):
    """
    Implements linear response eigensolver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - nstates: Number of excited states.
        - nto: The flag for natural transition orbital analysis.
        - nto_pairs: The number of NTO pairs in NTO analysis.
        - detach_attach: The flag for detachment/attachment density analysis.
        - cube_points: The number of cubic grid points in X, Y and Z directions.
    """

    def __init__(self, comm, ostream):
        """
        Initializes linear response eigensolver to default setup.
        """

        super().__init__(comm, ostream)

        self.nstates = 3

        self.nto = False
        self.nto_pairs = None
        self.detach_attach = False
        self.cube_points = [80, 80, 80]

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in linear response eigensolver.

        :param rsp_dict:
            The dictionary of response dict.
        :param method_dict:
            The dictionary of method rsp_dict.
        """

        if method_dict is None:
            method_dict = {}

        super().update_settings(rsp_dict, method_dict)

        if 'nstates' in rsp_dict:
            self.nstates = int(rsp_dict['nstates'])

        if 'nto' in rsp_dict:
            key = rsp_dict['nto'].lower()
            self.nto = True if key == 'yes' else False
        if 'nto_pairs' in rsp_dict:
            self.nto_pairs = int(rsp_dict['nto_pairs'])

        if 'detach_attach' in rsp_dict:
            key = rsp_dict['detach_attach'].lower()
            self.detach_attach = True if key == 'yes' else False

        if 'cube_points' in rsp_dict:
            self.cube_points = [
                int(x)
                for x in rsp_dict['cube_points'].replace(',', ' ').split()
            ]
            assert_msg_critical(
                len(self.cube_points) == 3, 'cube points: Need 3 integers')

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

        self.dist_bger = None
        self.dist_bung = None
        self.dist_e2bger = None
        self.dist_e2bung = None

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        if self.rank == mpi_master():
            self.print_header('Linear Response EigenSolver',
                              nstates=self.nstates)

        self.start_time = tm.time()

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'LinearResponseEigenSolver: not implemented for unrestricted case')

        if self.rank == mpi_master():
            orb_ene = scf_tensors['E']
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
        eri_dict = self.init_eri(molecule, basis)

        # DFT information
        dft_dict = self.init_dft(molecule, scf_tensors)

        # PE information
        pe_dict = self.init_pe(molecule, basis)

        timing_dict = {}

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
            self.read_vectors(rsp_vector_labels)

        # generate initial guess from scratch
        else:
            igs = self.initial_excitations(self.nstates, orb_ene, nocc, norb)
            bger, bung = self.setup_trials(igs, None)

            self.e2n_half_size(bger, bung, molecule, basis, scf_tensors,
                               eri_dict, dft_dict, pe_dict, timing_dict)

        profiler.check_memory_usage('Initial guess')

        excitations = {}
        exresiduals = {}
        relative_residual_norm = {}

        signal_handler = SignalHandler()
        signal_handler.add_sigterm_function(self.graceful_exit, molecule, basis,
                                            dft_dict, pe_dict,
                                            rsp_vector_labels)

        iter_per_trail_in_hours = None

        # start iterations
        for iteration in range(self.max_iter):

            iter_start_time = tm.time()

            profiler.start_timer(iteration, 'ReducedSpace')

            self.cur_iter = iteration

            e2gg = self.dist_bger.matmul_AtB(self.dist_e2bger, 2.0)
            e2uu = self.dist_bung.matmul_AtB(self.dist_e2bung, 2.0)
            s2ug = self.dist_bung.matmul_AtB(self.dist_bger, 4.0)

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

                x_ger = self.dist_bger.matmul_AB_no_gather(c_ger[:, k])
                x_ung = self.dist_bung.matmul_AB_no_gather(c_ung[:, k])

                e2x_ger = self.dist_e2bger.matmul_AB_no_gather(c_ger[:, k])
                e2x_ung = self.dist_e2bung.matmul_AB_no_gather(c_ung[:, k])

                s2x_ger = 2.0 * x_ger.data
                s2x_ung = 2.0 * x_ung.data

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
                    relative_residual_norm[k] = rn / xn
                else:
                    relative_residual_norm[k] = rn

                if relative_residual_norm[k] < self.conv_thresh:
                    excitations[k] = (w, x)
                else:
                    exresiduals[k] = (w, r)

            # write to output
            if self.rank == mpi_master():
                self.ostream.print_info(
                    '{:d} gerade trial vectors in reduced space'.format(
                        self.dist_bger.shape(1)))
                self.ostream.print_info(
                    '{:d} ungerade trial vectors in reduced space'.format(
                        self.dist_bung.shape(1)))
                self.ostream.print_blank()

                profiler.print_memory_subspace(
                    {
                        'dist_bger': self.dist_bger,
                        'dist_bung': self.dist_bung,
                        'dist_e2bger': self.dist_e2bger,
                        'dist_e2bung': self.dist_e2bung,
                        'exsolutions': excitations,
                        'exresiduals': exresiduals,
                    }, self.ostream)

                profiler.check_memory_usage(
                    'Iteration {:d} subspace'.format(iteration + 1))

                profiler.print_memory_tracing(self.ostream)

                self.print_iteration(relative_residual_norm, wn)

            profiler.stop_timer(iteration, 'ReducedSpace')

            # check convergence
            self.check_convergence(relative_residual_norm)

            if self.is_converged:
                break

            profiler.start_timer(iteration, 'Orthonorm.')

            # update trial vectors
            precond = {
                k: self.get_precond(orb_ene, nocc, norb, w)
                for k, w in enumerate(list(wn))
            }

            new_trials_ger, new_trials_ung = self.setup_trials(
                exresiduals, precond, self.dist_bger, self.dist_bung)

            exresiduals.clear()

            profiler.stop_timer(iteration, 'Orthonorm.')

            if self.rank == mpi_master():
                n_new_trials = new_trials_ger.shape(1) + new_trials_ung.shape(1)
            else:
                n_new_trials = None
            n_new_trials = self.comm.bcast(n_new_trials, root=mpi_master())

            if iter_per_trail_in_hours is not None:
                next_iter_in_hours = iter_per_trail_in_hours * n_new_trials
                if self.need_graceful_exit(next_iter_in_hours):
                    self.graceful_exit(molecule, basis, dft_dict, pe_dict,
                                       rsp_vector_labels)

            profiler.start_timer(iteration, 'FockBuild')

            self.e2n_half_size(new_trials_ger, new_trials_ung, molecule, basis,
                               scf_tensors, eri_dict, dft_dict, pe_dict,
                               timing_dict)

            iter_in_hours = (tm.time() - iter_start_time) / 3600
            iter_per_trail_in_hours = iter_in_hours / n_new_trials

            profiler.stop_timer(iteration, 'FockBuild')
            if self.dft or self.pe:
                profiler.update_timer(iteration, timing_dict)

            profiler.check_memory_usage(
                'Iteration {:d} sigma build'.format(iteration + 1))

        signal_handler.remove_sigterm_function()

        self.write_checkpoint(molecule, basis, dft_dict, pe_dict,
                              rsp_vector_labels)

        # converged?
        if self.rank == mpi_master():
            self.print_convergence('Linear response')

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of LR eigensolver')
        profiler.print_memory_usage(self.ostream)

        # calculate properties
        if self.is_converged:
            edip_rhs = self.get_prop_grad('electric dipole', 'xyz', molecule,
                                          basis, scf_tensors)
            lmom_rhs = self.get_prop_grad('linear momentum', 'xyz', molecule,
                                          basis, scf_tensors)
            mdip_rhs = self.get_prop_grad('magnetic dipole', 'xyz', molecule,
                                          basis, scf_tensors)

            eigvals = np.array([excitations[s][0] for s in range(self.nstates)])

            elec_trans_dipoles = np.zeros((self.nstates, 3))
            velo_trans_dipoles = np.zeros((self.nstates, 3))
            magn_trans_dipoles = np.zeros((self.nstates, 3))

            key_0 = list(excitations.keys())[0]
            x_0 = self.get_full_solution_vector(excitations[key_0][1])

            if self.rank == mpi_master():
                avail_mem = psutil.virtual_memory().available
                write_solution_to_file = (avail_mem <
                                          sys.getsizeof(x_0) * self.nstates * 2)
                if not write_solution_to_file:
                    eigvecs = np.zeros((x_0.size, self.nstates))

            for s in range(self.nstates):
                eigvec = self.get_full_solution_vector(excitations[s][1])

                if self.rank == mpi_master():
                    mo_occ = scf_tensors['C'][:, :nocc]
                    mo_vir = scf_tensors['C'][:, nocc:]
                    z_mat = eigvec[:eigvec.shape[0] // 2].reshape(
                        mo_occ.shape[1], mo_vir.shape[1]) * np.sqrt(2.0)
                    y_mat = eigvec[eigvec.shape[0] // 2:].reshape(
                        mo_occ.shape[1], mo_vir.shape[1]) * np.sqrt(2.0)

                if self.nto:
                    self.ostream.print_info(
                        'Running NTO analysis for S{:d}...'.format(s + 1))
                    self.ostream.flush()

                    if self.rank == mpi_master():
                        lam_diag, nto_mo = self.get_nto(z_mat, mo_occ, mo_vir)
                    else:
                        lam_diag = None
                        nto_mo = MolecularOrbitals()
                    lam_diag = self.comm.bcast(lam_diag, root=mpi_master())
                    nto_mo.broadcast(self.rank, self.comm)

                    self.write_nto_cubes(self.cube_points, molecule, basis, s,
                                         lam_diag, nto_mo, self.nto_pairs)

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

                    self.write_detach_attach_cubes(self.cube_points, molecule,
                                                   basis, s, dens_DA)

                if self.rank == mpi_master():
                    for ind, comp in enumerate('xyz'):
                        elec_trans_dipoles[s, ind] = np.vdot(
                            edip_rhs[ind], eigvec)
                        velo_trans_dipoles[s, ind] = np.vdot(
                            lmom_rhs[ind], eigvec) / (-eigvals[s])
                        magn_trans_dipoles[s, ind] = np.vdot(
                            mdip_rhs[ind], eigvec)

                    if write_solution_to_file:
                        append_rsp_solution_hdf5(self.checkpoint_file,
                                                 'S{:d}'.format(s + 1), eigvec)
                    else:
                        eigvecs[:, s] = eigvec[:]

            if self.nto or self.detach_attach:
                self.ostream.print_blank()
                self.ostream.flush()

            if self.rank == mpi_master():
                osc = (2.0 / 3.0) * np.sum(elec_trans_dipoles**2,
                                           axis=1) * eigvals
                rot_vel = (-1.0) * np.sum(
                    velo_trans_dipoles * magn_trans_dipoles,
                    axis=1) * rotatory_strength_in_cgs()

                ret_dict = {
                    'eigenvalues': eigvals,
                    'electric_transition_dipoles': elec_trans_dipoles,
                    'velocity_transition_dipoles': velo_trans_dipoles,
                    'magnetic_transition_dipoles': magn_trans_dipoles,
                    'oscillator_strengths': osc,
                    'rotatory_strengths': rot_vel,
                }

                if not write_solution_to_file:
                    ret_dict['eigenvectors'] = eigvecs

                return ret_dict

        return {}

    def get_full_solution_vector(self, solution):
        """
        Gets a full solution vector from the distributed solution.

        :param solution:
            The distributed solution as a tuple.

        :return:
            The full solution vector.
        """

        x_ger = solution.get_full_vector(0)
        x_ung = solution.get_full_vector(1)

        if self.rank == mpi_master():
            x_ger_full = np.hstack((x_ger, x_ger))
            x_ung_full = np.hstack((x_ung, -x_ung))
            return x_ger_full + x_ung_full
        else:
            return None

    def print_iteration(self, relative_residual_norm, ws):
        """
        Prints information of the iteration.

        :param relative_residual_norm:
            Relative residual norms.
        :param ws:
            Excitation energies.
        """

        width = 92
        output_header = '*** Iteration:   {} '.format(self.cur_iter + 1)
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

    def initial_excitations(self, nstates, ea, nocc, norb):
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

        excitations = [(i, a) for i in range(nocc) for a in range(nocc, norb)]
        excitation_energies = [ea[a] - ea[i] for i, a in excitations]

        w = {ia: w for ia, w in zip(excitations, excitation_energies)}
        n_exc = nocc * (norb - nocc)

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

            final[k] = (w[(i, a)], DistributedArray(X, self.comm))

        return final

    def precond_trials(self, excitations, precond):
        """
        Applies preconditioner to distributed trial vectors.

        :param excitations:
            The set of excitations.
        :param precond:
            The preconditioner.

        :return:
            The preconditioned gerade and ungerade trial vectors.
        """

        trials_ger = []
        trials_ung = []

        for k, (w, X) in excitations.items():
            if precond is not None:
                v = self.preconditioning(precond[k], X)
            else:
                v = X
            norms_2 = 2.0 * v.squared_norm(axis=0)
            vn = np.sqrt(np.sum(norms_2))

            if vn > self.small_thresh:
                norms = np.sqrt(norms_2)
                # gerade
                if norms[0] > self.small_thresh:
                    trials_ger.append(v.data[:, 0])
                # ungerade
                if norms[1] > self.small_thresh:
                    trials_ung.append(v.data[:, 1])

        new_ger = np.array(trials_ger).T
        new_ung = np.array(trials_ung).T

        dist_new_ger = DistributedArray(new_ger, self.comm, distribute=False)
        dist_new_ung = DistributedArray(new_ung, self.comm, distribute=False)

        return dist_new_ger, dist_new_ung

    def get_precond(self, orb_ene, nocc, norb, w):
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

        ediag, sdiag = self.construct_ed_sd_half(orb_ene, nocc, norb)

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

    def preconditioning(self, precond, v_in):
        """
        Applies preconditioner to a tuple of distributed trial vectors.

        :param precond:
            The preconditioner.
        :param v_in:
            The input trial vectors.

        :return:
            A tuple of distributed trail vectors after preconditioning.
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

        self.dist_bger = None
        self.dist_bung = None
        self.dist_e2bger = None
        self.dist_e2bung = None

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'LinearResponseEigenSolver: not implemented for unrestricted case')

        if self.rank == mpi_master():
            orb_ene = scf_tensors['E']
        else:
            orb_ene = None
        orb_ene = self.comm.bcast(orb_ene, root=mpi_master())
        norb = orb_ene.shape[0]
        nocc = molecule.number_of_alpha_electrons()

        # ERI information
        eri_dict = self.init_eri(molecule, basis)

        # DFT information
        dft_dict = self.init_dft(molecule, scf_tensors)

        # PE information
        pe_dict = self.init_pe(molecule, basis)

        timing_dict = {}

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

            igs[i] = (1.0, DistributedArray(X, self.comm))

        bger, bung = self.setup_trials(igs, precond=None, renormalize=False)

        self.e2n_half_size(bger, bung, molecule, basis, scf_tensors, eri_dict,
                           dft_dict, pe_dict, timing_dict)

        if self.rank == mpi_master():
            E2 = np.zeros((2 * n_exc, 2 * n_exc))

        for i in range(2 * n_exc):
            e2b_data = np.hstack((
                self.dist_e2bger.data[:, i:i + 1],
                self.dist_e2bung.data[:, i:i + 1],
            ))

            e2b = DistributedArray(e2b_data, self.comm, distribute=False)

            sigma = self.get_full_solution_vector(e2b)

            if self.rank == mpi_master():
                E2[:, i] = sigma[:]

        if self.rank == mpi_master():
            return E2
        else:
            return None
