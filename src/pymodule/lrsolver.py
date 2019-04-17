import numpy as np
import time as tm
import itertools
import math
import sys

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import NuclearPotentialIntegralsDriver
from .veloxchemlib import KineticEnergyIntegralsDriver
from .veloxchemlib import OverlapIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import ExcitationVector
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import szblock
from .outputstream import OutputStream
from .qqscheme import get_qq_scheme
from .errorhandler import assert_msg_critical


class LinearResponseSolver:
    """Implements linear response solver"""

    def __init__(self, a_ops=None, b_ops=None, frequencies=None):
        """Initializes linear response solver"""

        # operators and frequencies
        if a_ops:
            self.a_ops = a_ops
        else:
            self.a_ops = 'xyz'

        if b_ops:
            self.b_ops = b_ops
        else:
            self.b_ops = 'xyz'

        if frequencies:
            self.frequencies = frequencies
        else:
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

    def set_eri(self, eri_thresh, qq_type):
        """Sets screening in computation of electron repulsion integrals"""

        self.eri_thresh = eri_thresh
        self.qq_type = qq_type

    def set_solver(self, conv_thresh, max_iter):
        """Sets convergence threshold and maximum number of iterations"""

        self.conv_thresh = conv_thresh
        self.max_iter = max_iter

    def compute(self,
                mol_orbs,
                molecule,
                basis,
                comm,
                ostream=OutputStream(sys.stdout)):
        """Performs linear response calculation"""

        self.comm = comm.Clone()
        self.rank = comm.Get_rank()
        self.nodes = comm.Get_size()

        if self.rank == mpi_master():
            self.ostream = ostream

        self.molecule = molecule
        self.basis = basis

        nalpha = self.molecule.number_of_alpha_electrons()
        nbeta = self.molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'LinearResponseSolver: not yet implemented for unrestricted case')

        if self.rank == mpi_master():
            self.nocc = nalpha
            self.norb = mol_orbs.number_mos()
            self.mo = mol_orbs.alpha_to_numpy()
            self.ea = mol_orbs.ea_to_numpy()
            dens = mol_orbs.get_density(self.molecule)
            self.dens = (dens.alpha_to_numpy(0), dens.beta_to_numpy(0))

        # TODO: make use of 1e integrals from scf
        self.comp_1e_ints()

        self.eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        self.screening = self.eri_drv.compute(
            get_qq_scheme(self.qq_type), self.eri_thresh, self.molecule,
            self.basis)

        # TODO: make use of Fock matrices from scf
        self.comp_fock()

        # start linear response calculation

        self.start_time = tm.time()

        op_freq_keys, solutions = self.lr_solve(self.b_ops, self.frequencies)

        self.comm.Free()

        if self.rank == mpi_master():
            v1 = {op: v for op, v in zip(self.a_ops, self.get_rhs(self.a_ops))}

            lrs = {}
            for aop in self.a_ops:
                for col, (bop, w) in enumerate(op_freq_keys):
                    lrs[(aop, bop, w)] = -np.dot(v1[aop], solutions[:, col])
            return lrs
        else:
            return None

    def comp_1e_ints(self):
        """Computes 1e integrals"""

        overlap_drv = OverlapIntegralsDriver(self.comm)
        kinetic_drv = KineticEnergyIntegralsDriver(self.comm)
        potential_drv = NuclearPotentialIntegralsDriver(self.comm)
        dipole_drv = ElectricDipoleIntegralsDriver(self.comm)

        S = overlap_drv.compute(self.molecule, self.basis)
        T = kinetic_drv.compute(self.molecule, self.basis)
        V = potential_drv.compute(self.molecule, self.basis)
        Dpl = dipole_drv.compute(self.molecule, self.basis)

        if self.rank == mpi_master():
            self.overlap = S.to_numpy()
            self.hcore = T.to_numpy() - V.to_numpy()
            self.dipoles = (Dpl.x_to_numpy(), Dpl.y_to_numpy(),
                            Dpl.z_to_numpy())

    def comp_fock(self):
        """Computes Fock matrices"""

        if self.rank == mpi_master():
            dabs = (self.dens,)
        else:
            dabs = None

        fabs = self.get_two_el_fock(dabs)

        if self.rank == mpi_master():
            fa, fb = fabs[0]
            fa += self.hcore
            fb += self.hcore
            self.focks = (fa, fb)

    def get_rhs(self, ops):
        """Create right-hand sides of linear response equations"""

        if 'x' in ops or 'y' in ops or 'z' in ops:
            props = {k: v for k, v in zip('xyz', self.dipoles)}

        D = self.dens[0] + self.dens[1]
        S = self.overlap
        mo = self.mo

        matrices = tuple(
            mo.T @ (S @ D @ props[p].T - props[p].T @ D @ S) @ mo for p in ops)

        gradients = tuple(self.mat2vec(m) for m in matrices)
        return gradients

    def mat2vec(self, mat):

        excitations = list(self.get_excitations())

        z = [mat[i, j] for i, j in excitations]
        y = [mat[j, i] for i, j in excitations]

        return np.array(z + y)

    def lr_solve(self, ops, freqs):

        if self.rank == mpi_master():
            V1 = {op: v for op, v in zip(ops, self.get_rhs(ops))}
            igs = self.initial_guess(ops, freqs)
            b = self.setup_trials(igs)

            assert_msg_critical(
                np.any(b),
                'LinearResponseSolver.lr_solve: trial vector is empty')
        else:
            b = None

        e2b = self.e2n(b)

        if self.rank == mpi_master():
            s2b = self.s2n(b)

            od = self.get_orbital_diagonal()
            sd = self.get_overlap_diagonal()
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

                s2nn = self.s2n(solutions)
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

                conv_info = {'is_converged': self.is_converged}
            else:
                conv_info = None

            conv_info = self.comm.bcast(conv_info, root=mpi_master())
            self.is_converged = conv_info['is_converged']
            if self.is_converged:
                break

            if self.rank == mpi_master():
                new_trials = self.setup_trials(residuals, td=td, b=b)
                b = np.append(b, new_trials, axis=1)
            else:
                new_trials = None

            new_e2b = self.e2n(new_trials)

            if self.rank == mpi_master():
                new_s2b = self.s2n(new_trials)
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

            return op_freq_keys, solutions
        else:
            return None, None

    def initial_guess(self, ops, freqs):

        od = self.get_orbital_diagonal()
        sd = self.get_overlap_diagonal()
        dim = od.shape[0]

        ig = {}
        for op, grad in zip(ops, self.get_rhs(ops)):
            gn = np.linalg.norm(grad)
            for w in freqs:
                if gn < self.small_thresh:
                    ig[(op, w)] = np.zeros(dim)
                else:
                    td = od - w * sd
                    ig[(op, w)] = grad / td
        return ig

    def get_orbital_diagonal(self):

        orben = self.ea
        z = [2 * (orben[j] - orben[i]) for i, j in self.get_excitations()]
        e2c = np.array(z + z)
        return e2c

    def get_overlap_diagonal(self):

        lz = len(list(self.get_excitations()))
        s2d = 2.0 * np.ones(2 * lz)
        s2d[lz:] = -2.0
        return s2d

    def get_excitations(self):

        nocc = self.nocc
        norb = self.norb
        xv = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
        cre = xv.bra_unique_indexes()
        ann = xv.ket_unique_indexes()
        excitations = itertools.product(cre, ann)
        return excitations

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

    def e2n(self, vecs):

        if self.rank == mpi_master():
            assert_msg_critical(
                len(vecs.shape) == 2,
                'LinearResponseSolver.e2n: invalid shape of vecs')

            da, db = self.dens
            S = self.overlap
            mo = self.mo

            dks = []
            kns = []

            for col in range(vecs.shape[1]):
                vec = vecs[:, col]

                kN = self.vec2mat(vec).T
                kn = mo @ kN @ mo.T

                dak = kn.T @ S @ da - da @ S @ kn.T
                dbk = kn.T @ S @ db - db @ S @ kn.T

                dks.append((dak, dbk))
                kns.append(kn)

            dks = tuple(dks)
        else:
            dks = None

        fks = self.get_two_el_fock(dks)

        if self.rank == mpi_master():
            gv = np.zeros(vecs.shape)

            fa, fb = self.focks

            for col, (kn, (fak, fbk)) in enumerate(zip(kns, fks)):

                kfa = S @ kn @ fa - fa @ kn @ S
                kfb = S @ kn @ fb - fb @ kn @ S

                fat = fak + kfa
                fbt = fbk + kfb

                gao = S @ (da @ fat.T + db @ fbt.T) - (
                    fat.T @ da + fbt.T @ db) @ S
                gmo = mo.T @ gao @ mo

                gv[:, col] = -self.mat2vec(gmo)
            return gv
        else:
            return None

    def get_two_el_fock(self, dabs):

        # TODO: make this routine more general (for both rest and unrest)

        if self.rank == mpi_master():
            dts = []
            for dab in dabs:
                da, db = dab
                dt = da + db
                dts.append(dt)

                # Note: skip spin density for restricted case
                # ds = da - db
                # dts.append(ds)

            dens = AODensityMatrix(dts, denmat.rest)
        else:
            dens = AODensityMatrix()
        dens.broadcast(self.rank, self.comm)

        fock = AOFockMatrix(dens)

        # for i in range(0, 2 * len(dabs), 2):
        #    fock.set_fock_type(fockmat.rgenjk, i)
        #    fock.set_fock_type(fockmat.rgenk, i + 1)

        # Note: skip spin density for restricted case
        # for i in range(len(dabs)):
        #    fock.set_fock_type(fockmat.rgenjk, i)
        for i in range(fock.number_of_fock_matrices()):
            fock.set_fock_type(fockmat.rgenjk, i)

        self.eri_drv.compute(fock, dens, self.molecule, self.basis,
                             self.screening)
        fock.reduce_sum(self.rank, self.nodes, self.comm)

        fabs = []
        if self.rank == mpi_master():

            # for i in range(0, 2 * len(dabs), 2):
            #    ft = fock.to_numpy(i).T
            #    fs = -fock.to_numpy(i + 1).T
            #    fa = (ft + fs) / 2
            #    fb = (ft - fs) / 2
            #    fabs.append((fa, fb))

            # Note: skip spin density for restricted case
            for i in range(len(dabs)):
                ft = fock.to_numpy(i).T
                fa = ft * 0.5
                fb = ft * 0.5
                fabs.append((fa, fb))

            return tuple(fabs)
        else:
            return None

    def np2vlx(self, vec):

        nocc = self.nocc
        norb = self.norb
        xv = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
        zlen = len(vec) // 2
        z, y = vec[:zlen], vec[zlen:]
        xv.set_yzcoefficients(z, y)
        return xv

    def vec2mat(self, vec):

        xv = self.np2vlx(vec)
        kz = xv.get_zmatrix()
        ky = xv.get_ymatrix()

        rows = kz.number_of_rows() + ky.number_of_rows()
        cols = kz.number_of_columns() + ky.number_of_columns()

        kzy = np.zeros((rows, cols))
        kzy[:kz.number_of_rows(), ky.number_of_columns():] = kz.to_numpy()
        kzy[kz.number_of_rows():, :ky.number_of_columns()] = ky.to_numpy()

        return kzy

    def s2n(self, vecs):

        assert_msg_critical(
            len(vecs.shape) == 2,
            'LinearResponseSolver.s2n: invalid shape of vecs')

        D = self.dens[0] + self.dens[1]
        S = self.overlap
        mo = self.mo

        s2n_vecs = np.ndarray(vecs.shape)
        rows, columns = vecs.shape

        for c in range(columns):
            kappa = self.vec2mat(vecs[:, c]).T
            kappa_ao = mo @ kappa @ mo.T

            s2n_ao = kappa_ao.T @ S @ D - D @ S @ kappa_ao.T
            s2n_mo = mo.T @ S @ s2n_ao @ S @ mo
            s2n_vecs[:, c] = -self.mat2vec(s2n_mo)

        return s2n_vecs
