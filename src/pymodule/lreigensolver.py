import itertools
import numpy as np
import time as tm

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import ExcitationVector
from .veloxchemlib import mpi_master
from .veloxchemlib import szblock
from .veloxchemlib import rotatory_strength_in_cgs
from .lrmatvecdriver import LinearResponseMatrixVectorDriver
from .lrmatvecdriver import remove_linear_dependence
from .lrmatvecdriver import orthogonalize_gram_schmidt
from .lrmatvecdriver import normalize
from .lrmatvecdriver import construct_ed_sd
from .lrmatvecdriver import get_rhs
from .lrmatvecdriver import swap_xy
from .qqscheme import get_qq_scheme
from .qqscheme import get_qq_type
from .errorhandler import assert_msg_critical


class LinearResponseEigenSolver:
    """
    Implements linear response eigensolver.

    :param nstates:
        Number of excited states.
    :param eri_thresh:
        The electron repulsion integrals screening threshold.
    :param qq_type:
        The electron repulsion integrals screening scheme.
    :param conv_thresh:
        The convergence threshold for the solver.
    :param max_iter:
        The maximum number of solver iterations.
    :param cur_iter:
        Index of the current iteration.
    :param small_thresh:
        The norm threshold for a vector to be considered a zero vector.
    :param lindep_thresh:
        The threshold for removing linear dependence in the trial vectors.
    :param is_converged:
        The flag for convergence.
    :param comm:
        The MPI communicator.
    :param rank:
        The MPI rank.
    :param nodes:
        Number of MPI processes.
    :param ostream:
        The output stream.
    :param profiling:
        The flag for printing profiling information.
    """

    def __init__(self, comm, ostream):
        """
        Initializes linear response eigensolver to default setup.

        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
        """

        # number of states
        self.nstates = 3

        # ERI settings
        self.eri_thresh = 1.0e-15
        self.qq_type = 'QQ_DEN'

        # solver setup
        self.conv_thresh = 1.0e-5
        self.max_iter = 50
        self.cur_iter = 0
        self.small_thresh = 1.0e-10
        self.lindep_thresh = 1.0e-6
        self.is_converged = False

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        self.profiling = False

    def update_settings(self, settings):
        """
        Updates settings in linear response eigensolver.

        :param settings:
            The settings dictionary.
        """

        if 'nstates' in settings:
            self.nstates = int(settings['nstates'])

        if 'eri_thresh' in settings:
            self.eri_thresh = float(settings['eri_thresh'])
        if 'qq_type' in settings:
            self.qq_type = str(settings['qq_type'])

        if 'conv_thresh' in settings:
            self.conv_thresh = float(settings['conv_thresh'])
        if 'max_iter' in settings:
            self.max_iter = int(settings['max_iter'])
        if 'lindep_thresh' in settings:
            self.lindep_thresh = float(settings['lindep_thresh'])

        if 'profiling' in settings:
            key = settings['profiling'].lower()
            self.profiling = True if key in ['yes', 'y'] else False

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

        if self.profiling:
            import cProfile
            import pstats
            import io
            import os
            pr = cProfile.Profile()
            pr.enable()

        if self.rank == mpi_master():
            self.print_header()

        self.start_time = tm.time()

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'LinearResponseSolver: not implemented for unrestricted case')

        # make preparations
        if self.rank == mpi_master():
            mo = scf_tensors['C']
            ea = scf_tensors['E']
            nocc = nalpha
            norb = mo.shape[1]
            od, sd = construct_ed_sd(ea, nocc, norb)
        else:
            nocc = None
            norb = None

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        e2x_drv = LinearResponseMatrixVectorDriver(self.comm)

        # start linear response calculation
        if self.rank == mpi_master():
            igs = self.initial_excitations(self.nstates, ea, nocc, norb)
            b = self.setup_trials(igs)

            assert_msg_critical(
                np.any(b),
                'LinearResponseSolver.compute: trial vector is empty')
        else:
            b = None

        e2b = e2x_drv.e2n(b, scf_tensors, screening, molecule, basis)
        if self.rank == mpi_master():
            s2b = e2x_drv.s2n(b, scf_tensors, nocc)

        excitations = [None] * self.nstates
        exresiduals = [None] * self.nstates
        relative_residual_norm = {}

        # start iterations
        for i in range(self.max_iter):

            if self.rank == mpi_master():
                self.cur_iter = i
                ws = []

                # next solution

                E2 = b.T @ e2b
                S2 = b.T @ s2b

                wn, Xn = np.linalg.eig((np.linalg.solve(E2, S2)))
                p = list(reversed(wn.argsort()))
                wn = wn[p]
                Xn = Xn[:, p]
                for i in range(self.nstates):
                    norm = np.sqrt(Xn[:, i].T @ S2 @ Xn[:, i])
                    Xn[:, i] /= norm
                reduced_ev = zip(1.0 / wn[:self.nstates],
                                 Xn[:, :self.nstates].T)

                for k, (w, reduced_X) in enumerate(reduced_ev):
                    r = (e2b - w * s2b) @ reduced_X
                    X = b @ reduced_X
                    rn = np.linalg.norm(r)
                    xn = np.linalg.norm(X)
                    exresiduals[k] = (w, r)
                    excitations[k] = (w, X)
                    relative_residual_norm[k] = rn / xn
                    ws.append(w)

                # write to output
                self.print_iteration(relative_residual_norm, ws)

            # check convergence
            self.check_convergence(relative_residual_norm)

            if self.is_converged:
                break

            # update trial vectors
            if self.rank == mpi_master():
                tdx = {w: od - w * sd for w, x in excitations}
                new_trials = self.setup_trials(exresiduals, tdx=tdx, b=b)
                b = np.append(b, new_trials, axis=1)
            else:
                new_trials = None

            new_e2b = e2x_drv.e2n(new_trials, scf_tensors, screening, molecule,
                                  basis)
            if self.rank == mpi_master():
                new_s2b = e2x_drv.s2n(new_trials, scf_tensors, nocc)
                e2b = np.append(e2b, new_e2b, axis=1)
                s2b = np.append(s2b, new_s2b, axis=1)

        if self.profiling:
            pr.disable()
            s = io.StringIO()
            sortby = 'cumulative'
            ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
            ps.print_stats(20)
            if self.rank == mpi_master():
                for line in s.getvalue().split(os.linesep):
                    self.ostream.print_info(line)

        # converged?
        if self.rank == mpi_master():
            self.print_convergence()

            assert_msg_critical(
                self.is_converged,
                'LinearResponseEigenSolver.compute: failed to converge')

        dipole_rhs = get_rhs('dipole', 'xyz', molecule, basis, scf_tensors,
                             self.rank, self.comm)
        linmom_rhs = get_rhs('linear_momentum', 'xyz', molecule, basis,
                             scf_tensors, self.rank, self.comm)
        angmom_rhs = get_rhs('angular_momentum', 'xyz', molecule, basis,
                             scf_tensors, self.rank, self.comm)

        if self.rank == mpi_master():
            V_dipole = {op: V for op, V in zip('xyz', dipole_rhs)}
            V_linmom = {op: V for op, V in zip('xyz', linmom_rhs)}
            V_angmom = {op: V for op, V in zip('xyz', angmom_rhs)}

            elec_tms = {}
            velo_tms = {}
            magn_tms = {}

            eigvals = np.array([s[0] for s in excitations])
            eigvecs = [s[1] for s in excitations]

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

    def print_header(self):
        """
        Prints linear response eigensolver setup header to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header("Linear Response EigenSolver Setup")
        self.ostream.print_header(35 * "=")
        self.ostream.print_blank()

        str_width = 60

        cur_str = "Number of States          : " + str(self.nstates)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = "Max. Number of Iterations : " + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "Convergence Threshold     : " + \
            "{:.1e}".format(self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = "ERI Screening Scheme      : " + get_qq_type(self.qq_type)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "ERI Screening Threshold   : " + \
            "{:.1e}".format(self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_blank()

        self.ostream.flush()

    def print_iteration(self, relative_residual_norm, ws):
        """
        Prints information of the iteration.

        :param relative_residual_norm:
            Relative residual norms.
        :param ws:
            Excitation energies.
        """

        output_header = '*** Iteration:   {} '.format(self.cur_iter + 1)
        output_header += '* Residuals (Max,Min): '
        output_header += '{:.2e} and {:.2e}'.format(
            max(relative_residual_norm.values()),
            min(relative_residual_norm.values()))
        self.ostream.print_header(output_header.ljust(68))
        self.ostream.print_blank()
        for k, w in enumerate(ws):
            state_label = 'Excitation {}'.format(k + 1)
            rel_res = relative_residual_norm[k]
            output_iter = '{:<15s}: {:15.8f} '.format(state_label, w)
            output_iter += 'Residual Norm: {:.8f}'.format(rel_res)
            self.ostream.print_header(output_iter.ljust(68))
        self.ostream.print_blank()
        self.ostream.flush()

    def print_convergence(self):
        """
        Prints information after convergence.
        """

        output_conv = '*** '
        if self.is_converged:
            output_conv += 'Linear response converged'
        else:
            output_conv += 'Linear response NOT converged'
        output_conv += ' in {:d} iterations. '.format(self.cur_iter + 1)
        output_conv += 'Time: {:.2f} sec'.format(tm.time() - self.start_time)
        self.ostream.print_header(output_conv.ljust(68))
        self.ostream.print_blank()

    def check_convergence(self, relative_residual_norm):
        """
        Checks convergence.

        :param relative_residual_norm:
            Relative residual norms.
        """

        if self.rank == mpi_master():
            max_residual = max(relative_residual_norm.values())
            if max_residual < self.conv_thresh:
                self.is_converged = True

        self.is_converged = self.comm.bcast(self.is_converged,
                                            root=mpi_master())

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
            A list of initial excitations (excitation energy and vector).
        """

        xv = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
        excitations = list(
            itertools.product(xv.bra_unique_indexes(), xv.ket_unique_indexes()))

        ediag, sdiag = construct_ed_sd(ea, nocc, norb)
        excitation_energies = 0.5 * ediag

        w = {ia: w for ia, w in zip(excitations, excitation_energies)}

        final = []
        for (i, a) in sorted(w, key=w.get)[:nstates]:
            ia = excitations.index((i, a))
            Xn = np.zeros(2 * len(excitations))
            Xn[ia] = 1.0
            final.append((w[(i, a)], Xn))
        return final

    def setup_trials(self, excitations, tdx=None, b=None, renormalize=True):
        """
        Computes orthonormalized trial vectors.

        :param excitations:
            The set of excitations.
        :param tdx:
            The preconditioner.
        :param b:
            The subspace.
        :param renormalize:
            The flag for normalization.

        :return:
            The orthonormalized trial vectors.
        """

        trials = []

        for w, X in excitations:
            if tdx:
                trials.append(X / tdx[w])
            else:
                trials.append(X)
            trials.append(swap_xy(X))

        new_trials = np.array(trials).T

        if b is not None:
            new_trials = new_trials - np.matmul(b, np.matmul(b.T, new_trials))

        if trials and renormalize:
            new_trials = remove_linear_dependence(new_trials,
                                                  self.lindep_thresh)
            new_trials = orthogonalize_gram_schmidt(new_trials)
            new_trials = normalize(new_trials)

        return new_trials
