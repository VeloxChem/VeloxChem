import itertools
import numpy as np
import time as tm

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import ExcitationVector
from .veloxchemlib import mpi_master
from .veloxchemlib import szblock
from .lrmatvecdriver import LinearResponseMatrixVectorDriver
from .lrmatvecdriver import truncate_and_normalize
from .lrmatvecdriver import construct_ed_sd
from .lrmatvecdriver import lrmat2vec
from .lrmatvecdriver import swap_xy
from .qqscheme import get_qq_scheme
from .qqscheme import get_qq_type
from .errorhandler import assert_msg_critical


class LinearResponseEigenSolver:
    """Implements linear response solver"""

    def __init__(self, comm, ostream):
        """Initializes linear response solver.

        Initializes linear response solver to default setup.

        Parameters
        ----------
        comm
            The MPI communicator.
        ostream
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
        self.is_converged = False

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

    def update_settings(self, settings):
        """Updates settings in linear response solver.

        Updates settings in linear response solver.

        Parameters
        ----------
        settings
            The settings for the driver.
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

    def compute(self, molecule, basis, scf_tensors):
        """Performs linear response calculation.

        Performs linear response calculation for a molecule and a basis set.

        Parameters
        ----------
        molecule
            The molecule.
        basis
            The basis set.
        scf_tensors
            The tensors from converged SCF wavefunction.
        """

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

        # converged?
        if self.rank == mpi_master():
            self.print_convergence()

            assert_msg_critical(
                self.is_converged,
                'LinearResponseEigenSolver.compute: failed to converge')

        if self.rank == mpi_master():
            ops = 'xyz'
            rhs = self.get_rhs(ops, scf_tensors, nocc)
            V1 = {op: V for op, V in zip(ops, rhs)}

            tms = {}
            tms['w'] = np.array([s[0] for s in excitations])

            eigvecs = [s[1] for s in excitations]
            for op in ops:
                tms[op] = np.array([np.dot(V1[op], vec) for vec in eigvecs])

            trans_dipoles = [
                np.array([tms['x'][s], tms['y'][s], tms['z'][s]])
                for s in range(self.nstates)
            ]

            osc = 2.0 / 3.0 * tms['w'] * (tms['x']**2 + tms['y']**2 +
                                          tms['z']**2)

            return {
                'eigenvalues': tms['w'],
                'eigenvectors': np.array(eigvecs).T,
                'transition_dipoles': trans_dipoles,
                'oscillator_strengths': osc,
            }
        else:
            return {}

    def print_header(self):
        """Prints response eigen solver setup header to output stream.

        Prints excitation calculation setup details to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header("Linear Response EigenSolver Setup")
        self.ostream.print_header(35 * "=")
        self.ostream.print_blank()

        str_width = 60

        cur_str = "Number Of States          : " + str(self.nstates)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = "Max. Number Of Iterations : " + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "Convergence Threshold     : " + \
            "{:.1e}".format(self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = "ERI screening scheme      : " + get_qq_type(self.qq_type)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "ERI Screening Threshold   : " + \
            "{:.1e}".format(self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_blank()

        self.ostream.flush()

    def print_iteration(self, relative_residual_norm, ws):
        """Prints information of the iteration"""

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
        """Prints information after convergence"""

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
        """Checks convergence"""

        if self.rank == mpi_master():
            max_residual = max(relative_residual_norm.values())
            if max_residual < self.conv_thresh:
                self.is_converged = True

        self.is_converged = self.comm.bcast(self.is_converged,
                                            root=mpi_master())

    def get_rhs(self, ops, scf_tensors, nocc):
        """Create right-hand sides of linear response equations"""

        mo = scf_tensors['C']
        S = scf_tensors['S']
        D = scf_tensors['D'][0] + scf_tensors['D'][1]
        dipoles = scf_tensors['Mu']

        norb = mo.shape[1]

        if 'x' in ops or 'y' in ops or 'z' in ops:
            props = {k: v for k, v in zip('xyz', dipoles)}

        matrices = tuple(
            mo.T @ (S @ D @ props[p].T - props[p].T @ D @ S) @ mo for p in ops)

        gradients = tuple(lrmat2vec(m, nocc, norb) for m in matrices)
        return gradients

    def initial_excitations(self, nstates, ea, nocc, norb):

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
            new_trials = truncate_and_normalize(new_trials, self.small_thresh)

        return new_trials
