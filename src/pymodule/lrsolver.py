import numpy as np
import time as tm
import itertools
import math

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import ExcitationVector
from .veloxchemlib import mpi_master
from .veloxchemlib import szblock
from .lrmatvecdriver import LinearResponseMatrixVectorDriver
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
        """Updates settings in LR solver.

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

    def compute(self, mol_orbs, molecule, basis):
        """Performs linear response calculation"""

        basis = basis

        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'LinearResponseSolver: not yet implemented for unrestricted case')

        if self.rank == mpi_master():
            nocc = nalpha
            norb = mol_orbs.number_mos()
            mo = mol_orbs.alpha_to_numpy()
            ea = mol_orbs.ea_to_numpy()
            density = mol_orbs.get_density(molecule)
            dens = (density.alpha_to_numpy(0), density.beta_to_numpy(0))

            xv = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
            excitations = list(
                itertools.product(xv.bra_unique_indexes(),
                                  xv.ket_unique_indexes()))

            z = [2.0 * (ea[j] - ea[i]) for i, j in excitations]
            self.orb_diag = np.array(z + z)

            lz = len(excitations)
            self.ovl_diag = 2.0 * np.ones(2 * lz)
            self.ovl_diag[lz:] = -2.0

        else:
            nocc = None
            norb = None
            dens = None
            mo = None

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        e2x_drv = LinearResponseMatrixVectorDriver(self.comm)

        # TODO: make use of 1e integrals from scf
        overlap, hcore, dipoles = e2x_drv.comp_1e_ints(molecule, basis)

        # TODO: make use of Fock matrices from scf
        focks = e2x_drv.comp_fock(hcore, dens, screening, molecule, basis)

        tensors = {'C': mo, 'S': overlap, 'D': dens, 'F': focks, 'Mu': dipoles}
        shapes = {'nocc': nocc, 'norb': norb}

        # start linear response calculation

        self.start_time = tm.time()

        ops = self.b_ops
        freqs = self.frequencies

        if self.rank == mpi_master():
            V1 = {
                op: v for op, v in zip(ops, self.get_rhs(ops, tensors, shapes))
            }
            igs = self.initial_guess(freqs, V1)
            b = self.setup_trials(igs)

            assert_msg_critical(
                np.any(b),
                'LinearResponseSolver.compute: trial vector is empty')
        else:
            b = None

        e2b = e2x_drv.e2n(b, tensors, screening, molecule, basis)

        if self.rank == mpi_master():
            s2b = e2x_drv.s2n(b, tensors, shapes)

            od = self.orb_diag
            sd = self.ovl_diag
            td = {w: od - w * sd for w in freqs}

            op_freq_keys = list(igs.keys())
            solutions = np.zeros((b.shape[0], len(igs)))
            e2nn = np.zeros((e2b.shape[0], len(igs)))
            residuals = {}
            relative_residual_norm = {}

        for i in range(self.max_iter):

            if self.rank == mpi_master():
                self.cur_iter = i

                # next solution
                for col, (op, freq) in enumerate(op_freq_keys):
                    v = V1[op]

                    # reduced_solution = (b.T*v)/(b.T*(e2b-freq*s2b))
                    reduced_solution = np.linalg.solve(
                        np.matmul(b.T, e2b - freq * s2b), np.matmul(b.T, v))

                    solutions[:, col] = np.matmul(b, reduced_solution)
                    e2nn[:, col] = np.matmul(e2b, reduced_solution)

                s2nn = e2x_drv.s2n(solutions, tensors, shapes)
                nvs = []

                # next residual
                for col, (op, freq) in enumerate(op_freq_keys):
                    v = V1[op]
                    n = solutions[:, col]
                    r = e2nn[:, col] - freq * s2nn[:, col] - v
                    residuals[(op, freq)] = r
                    nvs.append(np.dot(n, v))
                    rn = np.linalg.norm(r)
                    nn = np.linalg.norm(n)
                    relative_residual_norm[(op, freq)] = rn / nn

                # write to output
                output_header = '*** Iteration:   {} '.format(i + 1)
                output_header += '* Residuals (Max,Min): '
                output_header += '{:.2e} and {:.2e}'.format(
                    max(relative_residual_norm.values()),
                    min(relative_residual_norm.values()))
                self.ostream.print_header(output_header.ljust(68))
                self.ostream.print_blank()
                for (op, freq), nv in zip(op_freq_keys, nvs):
                    ops_label = '<<{};{}>>_{}'.format(op, op, freq)
                    rel_res = relative_residual_norm[(op, freq)]
                    output_iter = '{:<15s}: {:15.8f} '.format(ops_label, -nv)
                    output_iter += 'Residual Norm: {:.8f}'.format(rel_res)
                    self.ostream.print_header(output_iter.ljust(68))
                self.ostream.print_blank()
                self.ostream.flush()

                max_residual = max(relative_residual_norm.values())
                if max_residual < self.conv_thresh:
                    self.is_converged = True

            self.is_converged = self.comm.bcast(self.is_converged,
                                                root=mpi_master())
            if self.is_converged:
                break

            if self.rank == mpi_master():
                new_trials = self.setup_trials(residuals, td=td, b=b)
                b = np.append(b, new_trials, axis=1)
            else:
                new_trials = None

            new_e2b = e2x_drv.e2n(new_trials, tensors, screening, molecule,
                                  basis)

            if self.rank == mpi_master():
                new_s2b = e2x_drv.s2n(new_trials, tensors, shapes)
                e2b = np.append(e2b, new_e2b, axis=1)
                s2b = np.append(s2b, new_s2b, axis=1)

        if self.rank == mpi_master():
            output_conv = '*** '
            if self.is_converged:
                output_conv += 'Linear response converged'
            else:
                output_conv += 'Linear response NOT converged'
            output_conv += ' in {:d} iterations. '.format(self.cur_iter + 1)
            output_conv += 'Time: {:.2f} sec'.format(tm.time() -
                                                     self.start_time)
            self.ostream.print_header(output_conv.ljust(68))
            self.ostream.print_blank()

            assert_msg_critical(
                self.is_converged,
                'LinearResponseSolver.compute: failed to converge')

        if self.rank == mpi_master():
            v1 = {
                op: v for op, v in zip(
                    self.a_ops, self.get_rhs(self.a_ops, tensors, shapes))
            }

            lrs = {}
            for aop in self.a_ops:
                for col, (bop, w) in enumerate(op_freq_keys):
                    lrs[(aop, bop, w)] = -np.dot(v1[aop], solutions[:, col])
            return lrs
        else:
            return None

    def get_rhs(self, ops, tensors, shapes):
        """Create right-hand sides of linear response equations"""

        mo = tensors['C']
        S = tensors['S']
        D = tensors['D'][0] + tensors['D'][1]
        dipoles = tensors['Mu']

        nocc = shapes['nocc']
        norb = shapes['norb']

        if 'x' in ops or 'y' in ops or 'z' in ops:
            props = {k: v for k, v in zip('xyz', dipoles)}

        matrices = tuple(
            mo.T @ (S @ D @ props[p].T - props[p].T @ D @ S) @ mo for p in ops)

        gradients = tuple(self.mat2vec(m, nocc, norb) for m in matrices)
        return gradients

    def mat2vec(self, mat, nocc, norb):

        xv = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
        excitations = list(
            itertools.product(xv.bra_unique_indexes(), xv.ket_unique_indexes()))

        z = [mat[i, j] for i, j in excitations]
        y = [mat[j, i] for i, j in excitations]
        return np.array(z + y)

    def initial_guess(self, freqs, V1):

        od = self.orb_diag
        sd = self.ovl_diag
        dim = od.shape[0]

        ig = {}
        for op, grad in V1.items():
            gn = np.linalg.norm(grad)
            for w in freqs:
                if gn < self.small_thresh:
                    ig[(op, w)] = np.zeros(dim)
                else:
                    td = od - w * sd
                    ig[(op, w)] = grad / td
        return ig

    def setup_trials(self, vectors, td=None, b=None, renormalize=True):

        trials = []

        for (op, freq) in vectors:
            vec = np.array(vectors[(op, freq)])

            if td is not None:
                v = vec / td[freq]
            else:
                v = vec

            if np.linalg.norm(v) > self.small_thresh:
                trials.append(v)
                if freq > self.small_thresh:
                    trials.append(self.swap(v))

        new_trials = np.array(trials).transpose()

        if b is not None:
            new_trials = new_trials - np.matmul(b, np.matmul(b.T, new_trials))

        if trials and renormalize:
            t = self.get_transform(new_trials)
            truncated = np.matmul(new_trials, t)

            # S12 = (truncated.T*truncated).invsqrt()
            tt = np.matmul(truncated.T, truncated)
            evals, evecs = np.linalg.eigh(tt)
            evals_invsqrt = [1.0 / math.sqrt(x) for x in evals]
            S12 = np.matmul(evecs, np.matmul(np.diag(evals_invsqrt), evecs.T))

            new_trials = np.matmul(truncated, S12)

        return new_trials

    def get_transform(self, basis):

        Sb = np.matmul(basis.T, basis)
        l, T = np.linalg.eigh(Sb)
        b_norm = np.sqrt(Sb.diagonal())
        mask = l > b_norm * self.small_thresh
        return T[:, mask]

    @staticmethod
    def swap(xy):
        """Swaps X and Y parts of response vector"""

        assert_msg_critical(
            len(xy.shape) == 1 or len(xy.shape) == 2,
            'LinearResponseSolver.swap: invalid shape of XY')

        half_rows = xy.shape[0] // 2
        yx = xy.copy()

        if len(xy.shape) == 1:
            yx[:half_rows] = xy[half_rows:]
            yx[half_rows:] = xy[:half_rows]

        elif len(xy.shape) == 2:
            yx[:half_rows, :] = xy[half_rows:, :]
            yx[half_rows:, :] = xy[:half_rows, :]

        return yx
