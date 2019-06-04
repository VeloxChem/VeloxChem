import numpy as np
import time as tm

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import mpi_master
from .lrmatvecdriver import LinearResponseMatrixVectorDriver
from .lrmatvecdriver import truncate_and_normalize
from .lrmatvecdriver import construct_ed_sd
from .lrmatvecdriver import lrmat2vec
from .lrmatvecdriver import swap_xy
from .qqscheme import get_qq_scheme
from .errorhandler import assert_msg_critical


class LinearResponseSolver:
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

        # operators and frequencies
        self.a_ops = 'xyz'
        self.b_ops = 'xyz'
        self.frequencies = (0,)

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

        if 'a_ops' in settings:
            self.a_ops = settings['a_ops']
        if 'b_ops' in settings:
            self.b_ops = settings['b_ops']
        if 'frequencies' in settings:
            self.frequencies = settings['frequencies']

        if 'eri_thresh' in settings:
            self.eri_thresh = settings['eri_thresh']
        if 'qq_type' in settings:
            self.qq_type = settings['qq_type']

        if 'conv_thresh' in settings:
            self.conv_thresh = settings['conv_thresh']
        if 'max_iter' in settings:
            self.max_iter = settings['max_iter']

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
            td = {w: od - w * sd for w in self.frequencies}
        else:
            nocc = None
            norb = None

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        e2x_drv = LinearResponseMatrixVectorDriver(self.comm)

        # start linear response calculation
        if self.rank == mpi_master():
            V1 = {
                op: v for op, v in zip(
                    self.b_ops, self.get_rhs(self.b_ops, scf_tensors, nocc))
            }
            igs = self.initial_guess(self.frequencies, V1, od, sd, td)
            b = self.setup_trials(igs)

            assert_msg_critical(
                np.any(b),
                'LinearResponseSolver.compute: trial vector is empty')
        else:
            b = None

        e2b = e2x_drv.e2n(b, scf_tensors, screening, molecule, basis)
        if self.rank == mpi_master():
            s2b = e2x_drv.s2n(b, scf_tensors, nocc)

        solutions = {}
        residuals = {}
        relative_residual_norm = {}

        # start iterations
        for i in range(self.max_iter):

            if self.rank == mpi_master():
                self.cur_iter = i
                nvs = []

                # next solution
                for op, freq in igs:
                    v = V1[op]
                    reduced_solution = np.linalg.solve(
                        np.matmul(b.T, e2b - freq * s2b), np.matmul(b.T, v))
                    solutions[(op, freq)] = np.matmul(b, reduced_solution)
                    residuals[(op, freq)] = np.matmul(e2b - freq * s2b,
                                                      reduced_solution) - v

                    r = residuals[(op, freq)]
                    n = solutions[(op, freq)]

                    nv = np.dot(n, v)
                    nvs.append((op, freq, nv))

                    rn = np.linalg.norm(r)
                    nn = np.linalg.norm(n)
                    relative_residual_norm[(op, freq)] = rn / nn

                # write to output
                self.print_iteration(relative_residual_norm, nvs)

            # check convergence
            self.check_convergence(relative_residual_norm)

            if self.is_converged:
                break

            # update trial vectors
            if self.rank == mpi_master():
                new_trials = self.setup_trials(residuals, td=td, b=b)
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
                'LinearResponseSolver.compute: failed to converge')

        # calculate properties
        if self.rank == mpi_master():
            v1 = {
                op: v for op, v in zip(
                    self.a_ops, self.get_rhs(self.a_ops, scf_tensors, nocc))
            }

            lrs = {}
            for aop in self.a_ops:
                for bop, w in solutions:
                    lrs[(aop, bop, w)] = -np.dot(v1[aop], solutions[(bop, w)])
            return lrs
        else:
            return None

    def print_iteration(self, relative_residual_norm, nvs):
        """Prints information of the iteration"""

        output_header = '*** Iteration:   {} '.format(self.cur_iter + 1)
        output_header += '* Residuals (Max,Min): '
        output_header += '{:.2e} and {:.2e}'.format(
            max(relative_residual_norm.values()),
            min(relative_residual_norm.values()))
        self.ostream.print_header(output_header.ljust(68))
        self.ostream.print_blank()
        for op, freq, nv in nvs:
            ops_label = '<<{};{}>>_{}'.format(op, op, freq)
            rel_res = relative_residual_norm[(op, freq)]
            output_iter = '{:<15s}: {:15.8f} '.format(ops_label, -nv)
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

    def initial_guess(self, freqs, V1, od, sd, td):

        dim = od.shape[0]

        ig = {}
        for op, grad in V1.items():
            gn = np.linalg.norm(grad)
            for w in freqs:
                if gn < self.small_thresh:
                    ig[(op, w)] = np.zeros(dim)
                else:
                    ig[(op, w)] = grad / td[w]
        return ig

    def setup_trials(self, vectors, td=None, b=None, renormalize=True):

        trials = []

        for (op, freq) in vectors:
            vec = vectors[(op, freq)]

            if td is not None:
                v = vec / td[freq]
            else:
                v = vec

            if np.linalg.norm(v) > self.small_thresh:
                trials.append(v)
                if freq > self.small_thresh:
                    trials.append(swap_xy(v))

        new_trials = np.array(trials).T

        if b is not None:
            new_trials = new_trials - np.matmul(b, np.matmul(b.T, new_trials))

        if trials and renormalize:
            new_trials = truncate_and_normalize(new_trials, self.small_thresh)

        return new_trials
