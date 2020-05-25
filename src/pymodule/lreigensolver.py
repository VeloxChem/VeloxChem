import itertools
import numpy as np
import time as tm
import math

from .veloxchemlib import ExcitationVector
from .veloxchemlib import mpi_master
from .veloxchemlib import rotatory_strength_in_cgs
from .veloxchemlib import szblock
from .profiler import Profiler
from .distributedarray import DistributedArray
from .signalhandler import SignalHandler
from .linearsolver import LinearSolver
from .errorhandler import assert_msg_critical
from .checkpoint import check_rsp_hdf5


class LinearResponseEigenSolver(LinearSolver):
    """
    Implements linear response eigensolver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - nstates: Number of excited states.
    """

    def __init__(self, comm, ostream):
        """
        Initializes linear response eigensolver to default setup.
        """

        super().__init__(comm, ostream)

        self.nstates = 3

    def update_settings(self, rsp_dict, method_dict={}):
        """
        Updates response and method settings in linear response eigensolver.

        :param rsp_dict:
            The dictionary of response dict.
        :param method_dict:
            The dictionary of method rsp_dict.
        """

        super().update_settings(rsp_dict, method_dict)

        if 'nstates' in rsp_dict:
            self.nstates = int(rsp_dict['nstates'])

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
        sqrt_2 = math.sqrt(2.0)

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
                c_ger = np.matmul(tmat, evecs[:, :self.nstates])
                c_ung = wn * np.linalg.multi_dot([e2uu_inv, s2ug, c_ger])

                for k in range(self.nstates):
                    c_ger_k = c_ger[:, k]
                    c_ung_k = c_ung[:, k]
                    norm = np.sqrt(
                        np.linalg.multi_dot([c_ung_k.T, s2ug, c_ger_k]) +
                        np.linalg.multi_dot([c_ger_k.T, s2ug.T, c_ung_k]))
                    c_ger[:, k] /= norm
                    c_ung[:, k] /= norm
            else:
                c_ger = None
                c_ung = None
                wn = None
            c_ger = self.comm.bcast(c_ger, root=mpi_master())
            c_ung = self.comm.bcast(c_ung, root=mpi_master())
            wn = self.comm.bcast(wn, root=mpi_master())

            for k in range(self.nstates):
                w = wn[k]

                x_ger = self.dist_bger.matmul_AB(c_ger[:, k])
                x_ung = self.dist_bung.matmul_AB(c_ung[:, k])

                e2x_ger = self.dist_e2bger.matmul_AB(c_ger[:, k])
                e2x_ung = self.dist_e2bung.matmul_AB(c_ung[:, k])

                if self.rank == mpi_master():
                    s2x_ger = 2.0 * x_ger
                    s2x_ung = 2.0 * x_ung

                    r_ger = e2x_ger - w * s2x_ung
                    r_ung = e2x_ung - w * s2x_ger
                else:
                    r_ger, r_ung = None, None

                if self.rank == mpi_master():
                    r = (r_ger, r_ung)
                    x = (x_ger, x_ung)

                    rn = sqrt_2 * math.sqrt(sum([np.dot(p, p) for p in r]))
                    xn = sqrt_2 * math.sqrt(sum([np.dot(p, p) for p in x]))

                    if xn != 0:
                        relative_residual_norm[k] = rn / xn
                    else:
                        relative_residual_norm[k] = rn

                    state_k_converged = (relative_residual_norm[k] <
                                         self.conv_thresh)
                else:
                    state_k_converged = False
                state_k_converged = self.comm.bcast(state_k_converged,
                                                    root=mpi_master())

                if state_k_converged:
                    x = (DistributedArray(x_ger, self.comm),
                         DistributedArray(x_ung, self.comm))
                    excitations[k] = (w, x)
                else:
                    r = (DistributedArray(r_ger, self.comm),
                         DistributedArray(r_ung, self.comm))
                    exresiduals[k] = (w, r)

            relative_residual_norm = self.comm.bcast(relative_residual_norm,
                                                     root=mpi_master())

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
        dipole_rhs = self.get_rhs('dipole', 'xyz', molecule, basis, scf_tensors)
        linmom_rhs = self.get_rhs('linear_momentum', 'xyz', molecule, basis,
                                  scf_tensors)
        angmom_rhs = self.get_rhs('angular_momentum', 'xyz', molecule, basis,
                                  scf_tensors)

        eigvals = np.array([excitations[s][0] for s in range(self.nstates)])
        # TODO: check system memory
        eigvecs = [
            self.get_full_solution_vector(excitations[s][1])
            for s in range(self.nstates)
        ]

        if self.rank == mpi_master() and self.is_converged:
            V_dipole = {op: V for op, V in zip('xyz', dipole_rhs)}
            V_linmom = {op: V for op, V in zip('xyz', linmom_rhs)}
            V_angmom = {op: V for op, V in zip('xyz', angmom_rhs)}

            elec_tms = {}
            velo_tms = {}
            magn_tms = {}

            for comp in 'xyz':
                elec_tms[comp] = np.array(
                    [np.dot(V_dipole[comp], vec) for vec in eigvecs])
                velo_tms[comp] = -1.0 / eigvals * np.array(
                    [np.dot(V_linmom[comp], vec) for vec in eigvecs])
                magn_tms[comp] = 0.5 * np.array(
                    [np.dot(V_angmom[comp], vec) for vec in eigvecs])

            elec_trans_dipoles = [
                np.array([elec_tms['x'][s], elec_tms['y'][s], elec_tms['z'][s]])
                for s in range(self.nstates)
            ]

            velo_trans_dipoles = [
                np.array([velo_tms['x'][s], velo_tms['y'][s], velo_tms['z'][s]])
                for s in range(self.nstates)
            ]

            magn_trans_dipoles = [
                np.array([magn_tms['x'][s], magn_tms['y'][s], magn_tms['z'][s]])
                for s in range(self.nstates)
            ]

            osc = 2.0 / 3.0 * eigvals * (elec_tms['x']**2 + elec_tms['y']**2 +
                                         elec_tms['z']**2)

            rot_vel = (velo_tms['x'] * magn_tms['x'] +
                       velo_tms['y'] * magn_tms['y'] +
                       velo_tms['z'] * magn_tms['z'])

            rot_vel *= rotatory_strength_in_cgs()

            return {
                'eigenvalues': eigvals,
                'eigenvectors': np.array(eigvecs).T,
                'electric_transition_dipoles': elec_trans_dipoles,
                'velocity_transition_dipoles': velo_trans_dipoles,
                'magnetic_transition_dipoles': magn_trans_dipoles,
                'oscillator_strengths': osc,
                'rotatory_strengths': rot_vel,
            }
        else:
            return {}

    def get_full_solution_vector(self, solution):
        """
        Gets a full solution vector from the distributed solution.

        :param solution:
            The distributed solution as a tuple.

        :return:
            The full solution vector.
        """

        (dist_x_ger, dist_x_ung) = solution

        x_ger = dist_x_ger.get_full_vector()
        x_ung = dist_x_ung.get_full_vector()

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

        xv = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
        excitations = list(
            itertools.product(xv.bra_unique_indexes(), xv.ket_unique_indexes()))

        excitation_energies = [ea[a] - ea[i] for i, a in excitations]

        w = {ia: w for ia, w in zip(excitations, excitation_energies)}

        final = {}
        for k, (i, a) in enumerate(sorted(w, key=w.get)[:nstates]):
            if self.rank == mpi_master():
                ia = excitations.index((i, a))
                n_exc = len(excitations)

                Xn = np.zeros(2 * n_exc)
                Xn[ia] = 1.0

                Xn_T = np.zeros(2 * n_exc)
                Xn_T[:n_exc] = Xn[n_exc:]
                Xn_T[n_exc:] = Xn[:n_exc]

                Xn_ger = 0.5 * (Xn + Xn_T)[:n_exc]
                Xn_ung = 0.5 * (Xn - Xn_T)[:n_exc]

                X = (Xn_ger, Xn_ung)
            else:
                X = (None, None)

            final[k] = ((w[(i, a)],
                         tuple([DistributedArray(v, self.comm) for v in X])))

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
        sqrt_2 = math.sqrt(2.0)

        for k, (w, X) in excitations.items():
            if precond is not None:
                v = self.preconditioning(precond[k], X)
            else:
                v = X
            norms = [sqrt_2 * p.norm() for p in v]
            vn = math.sqrt(sum([n**2 for n in norms]))

            # TODO: check gerade and ungerade norms

            if vn > self.small_thresh:
                # gerade
                trials_ger.append(v[0].data)
                # ungerade
                trials_ung.append(v[1].data)

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

        return tuple(
            [DistributedArray(p, self.comm) for p in [pa_diag, pb_diag]])

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

        pa, pb = [p.data for p in precond]

        v_in_rg, v_in_ru = [v.data for v in v_in]

        v_out_rg = pa * v_in_rg + pb * v_in_ru
        v_out_ru = pb * v_in_rg + pa * v_in_ru

        return tuple([
            DistributedArray(v, self.comm, distribute=False)
            for v in [v_out_rg, v_out_ru]
        ])
