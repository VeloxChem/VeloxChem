import numpy as np
import time as tm

from .veloxchemlib import ExcitationVector
from .veloxchemlib import ElectricDipoleIntegralsDriver
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import OverlapIntegralsDriver
from .veloxchemlib import KineticEnergyIntegralsDriver
from .veloxchemlib import NuclearPotentialIntegralsDriver
from .veloxchemlib import mpi_master
from .veloxchemlib import szblock
from .veloxchemlib import denmat, fockmat
from .veloxchemlib import ericut
from .lrmatvecdriver import LinearResponseMatrixVectorDriver
from .lrmatvecdriver import lrmat2vec
from .aofockmatrix import AOFockMatrix
from .aodensitymatrix import AODensityMatrix
from .qqscheme import get_qq_scheme


class ComplexResponse:
    """ Provides functionality to solve complex linear
    response equations. """

    def __init__(self, comm):

        self.qq_type = 'QQ_DEN'
        self.eri_thresh = 1.0e-12

        self.is_converged = False

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

    def get_orbital_numbers(self, mol_orbs, task):
        """Requests information about orbitals like the total number,
        number of occupied ones and indices of orbitals involved in
        excitation processes.
        """

        nocc = task.molecule.number_of_electrons() // 2
        norb = mol_orbs.number_mos()
        xv = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
        cre = xv.bra_indexes()
        ann = xv.ket_indexes()

        return nocc, norb, xv, cre, ann

    def paired(self, v_xy):
        """Returns paired trial vector.
        """

        v_yx = v_xy.copy()
        half_rows = v_xy.shape[0] // 2
        v_yx[:half_rows] = v_xy[half_rows:]
        v_yx[half_rows:] = v_xy[:half_rows]

        return v_yx

    def decomp_trials(self, vecs):
        """Decomposes trial vectors into their 4 respective parts.
        """

        quarter_rows = vecs.shape[0] // 4
        half_rows = 2 * quarter_rows

        if len(vecs.shape) != 1:
            realger = []
            realung = []
            imagger = []
            imagung = []
            for vec in range(len(vecs[0, :])):
                realger.append(vecs[:quarter_rows, vec])
                realung.append(vecs[quarter_rows:half_rows, vec])
                imagung.append(vecs[half_rows:-quarter_rows, vec])
                imagger.append(vecs[-quarter_rows:, vec])
        else:
            realger = vecs[:quarter_rows]
            realung = vecs[quarter_rows:half_rows]
            imagung = vecs[half_rows:-quarter_rows]
            imagger = vecs[-quarter_rows:]

        return np.array(realger).T, np.array(realung).T, np.array(
            imagung).T, np.array(imagger).T

    def assemble_subsp(self, realvec, imagvec, threshold=1.0e-10):
        """Assembles subspace out of real and imaginary parts of trials,
        if their norm exceeds a certain threshold (zero vectors shouldn't
        be added).
        """

        space = []
        if len(realvec.shape) != 1:
            for vec in range(len(realvec[0, :])):
                if np.linalg.norm(realvec[:, vec]) > threshold:
                    space.append(realvec[:, vec])
                if np.linalg.norm(imagvec[:, vec]) > threshold:
                    space.append(imagvec[:, vec])

        else:
            if np.linalg.norm(realvec) > threshold:
                space.append(realvec)
            if np.linalg.norm(imagvec) > threshold:
                space.append(imagvec)

        return np.array(space).T

    def decomp_sym(self, vecs):
        """Decomposes gradient into gerade and ungerade parts.
        """

        if len(vecs.shape) != 1:
            ger = []
            ung = []
            for vec in range(len(vecs[0, :])):
                vecp = self.paired(vec)
                ger.append(.5 * (vec + vecp))
                ung.append(.5 * (vec - vecp))
        else:
            vecp = self.paired(vecs)
            ger = .5 * (vecs + vecp)
            ung = .5 * (vecs - vecp)

        return np.array(ger).T, np.array(ung).T

    def rm_lin_depend(self, basis, threshold=1e-15):
        """Removes linear dependencies from input basis vectors.
        """

        sb = np.matmul(basis.T, basis)
        l, T = np.linalg.eigh(sb)
        mask = l > threshold

        return T[:, mask]

    def orthogonalize_gram_schmidt(self, tvecs):
        """Applies modified Gram Schmidt orthogonalization to trial vectors.

        Applies modified Gram Schmidt orthogonalization to trial vectors.

        Parameters
        ----------
        tvecs
            The trial vectors.
        """

        if tvecs.shape[1] > 0:

            f = 1.0 / np.linalg.norm(tvecs[:, 0])
            tvecs[:, 0] *= f

            for i in range(1, tvecs.shape[1]):
                for j in range(i):
                    f = np.dot(tvecs[:, i], tvecs[:, j]) / np.dot(
                        tvecs[:, j], tvecs[:, j])
                    tvecs[:, i] -= f * tvecs[:, j]
                f = 1.0 / np.linalg.norm(tvecs[:, i])
                tvecs[:, i] *= f

        return tvecs

    def normalize(self, vecs):
        """Normalizes vectors by dividing by vector norm.
        """

        if len(vecs.shape) != 1:
            for vec in range(vecs.shape[1]):
                invnorm = 1.0 / np.linalg.norm(vecs[:, vec])
                vecs[:, vec] *= invnorm
        else:
            invnorm = 1.0 / np.linalg.norm(vecs)
            vecs *= invnorm

        return vecs

    def construct_ed(self, mol_orbs, task):
        """Returns the E0 matrix and its diagonal elements as an array.
        """

        nocc, norb, xv, cre, ann = self.get_orbital_numbers(mol_orbs, task)
        ediag = 2 * xv.diagonal_to_numpy(mol_orbs) + 0.0001
        ediag = np.append(ediag, ediag)
        e0 = np.zeros((len(ediag), len(ediag)))
        np.fill_diagonal(e0, ediag)

        return e0, ediag

    def construct_sd(self, mol_orbs, task):
        """Returns the S0 matrix and its diagonal elements as an array.
        """

        nocc, norb, xv, cre, ann = self.get_orbital_numbers(mol_orbs, task)
        sdiag = 2 * np.ones((2, nocc * (norb - nocc)))
        sdiag[1, :] *= -1
        sdiag = sdiag.flatten()
        s0 = np.zeros((len(sdiag), len(sdiag)))
        np.fill_diagonal(s0, sdiag)

        return s0, sdiag

    def get_precond(self, mol_orbs, task, w, d):
        """Constructs the preconditioner matrix.
        """

        start_time = tm.time()

        # spawning needed components

        e0, ediag = self.construct_ed(mol_orbs, task)
        s0, sdiag = self.construct_sd(mol_orbs, task)
        ediag_sq = ediag**2
        sdiag_sq = sdiag**2
        sdiag_fp = sdiag**4
        e0_sq = e0.copy()
        np.fill_diagonal(e0_sq, ediag_sq)
        s0_sq = s0.copy()
        np.fill_diagonal(s0_sq, sdiag_sq)
        s0_fp = s0.copy()
        np.fill_diagonal(s0_fp, sdiag_fp)
        w_sq = w**2
        d_sq = d**2

        # constructing matrix blocks

        a_mat = np.matmul(e0, (e0_sq - (w_sq - d_sq) * s0_sq))
        b_mat = np.matmul((w * s0), (e0_sq - (w_sq + d_sq) * s0_sq))
        c_mat = np.matmul((d * s0), (e0_sq + (w_sq + d_sq) * s0_sq))
        d_mat = np.matmul((2 * w * d * e0), s0_sq)
        p_mat = np.linalg.inv((e0_sq - (w_sq - d_sq) * s0_sq)**2 +
                              (4 * w_sq * d_sq * s0_fp))
        pa_mat = np.matmul(p_mat, a_mat)
        pb_mat = np.matmul(p_mat, b_mat)
        pc_mat = np.matmul(p_mat, c_mat)
        pd_mat = np.matmul(p_mat, d_mat)

        # creating empty matrix with wanted size

        size = e0.shape[0]
        precond = np.zeros((4 * size, 4 * size))

        # filling diagonal elements (de) of matrix blocks:

        for de in range(size):
            precond[de, de] = pa_mat[de, de]
            precond[de + size, de + size] = pa_mat[de, de]
            precond[-de - 1, -de - 1] = -pa_mat[-de - 1, -de - 1]
            precond[-(de + size) - 1, -(de + size) -
                    1] = -pa_mat[-de - 1, -de - 1]

            precond[de, de + size] = pb_mat[de, de]
            precond[de + size, de] = pb_mat[de, de]
            precond[-de - 1, -(de + size) - 1] = -pb_mat[-de - 1, -de - 1]
            precond[-(de + size) - 1, -de - 1] = -pb_mat[-de - 1, -de - 1]

            precond[size - de - 1, -(de + size) - 1] = pc_mat[-de - 1, -de - 1]
            precond[-(de + size) - 1, size - de - 1] = pc_mat[-de - 1, -de - 1]
            precond[de + size, de - size] = pc_mat[de, de]
            precond[de - size, de + size] = pc_mat[de, de]

            precond[de, de - size] = pd_mat[de, de]
            precond[de - size, de] = pd_mat[de, de]
            precond[de + size, de - 2 * size] = pd_mat[de, de]
            precond[de - 2 * size, de + size] = pd_mat[de, de]

        print('Precondition matrix created in',
              '{:.2f}'.format(tm.time() - start_time), 'sec.')

        return precond

    def get_rhs(self, mol_orbs, task, scf_tensors, ops):
        """Creates right-hand sides of complex linear response equations
        enabled gradients: dipole length
        """

        mo = scf_tensors['C']
        s = scf_tensors['S']
        d = scf_tensors['D'][0] + scf_tensors['D'][1]
        dipoles = scf_tensors['Mu']

        nocc = task.molecule.number_of_alpha_electrons()
        norb = mo.shape[1]

        if 'x' in ops or 'y' in ops or 'z' in ops:
            prop = {k: v for k, v in zip('xyz', dipoles)}

        # creating mo gradient matrices and converting them into vectors

        matrices = tuple(
            [mo.T @ (s @ d @ prop[p] - prop[p] @ d @ s) @ mo for p in ops])
        gradients = tuple([lrmat2vec(m, nocc, norb) for m in matrices])

        return gradients

    def initial_guess(self, mol_orbs, task, scf_tensors, d, ops, freqs, precond):
        """Creating initial guess (un-orthonormalized trials) out of gradients.
        """

        ig = {}
        for op, grad in zip(ops, self.get_rhs(mol_orbs, task, scf_tensors, ops)):
            gradger, gradung = self.decomp_sym(grad)
            grad = np.array(
                [gradger.real, gradung.real, -gradung.imag,
                 -gradger.imag]).flatten()
            gn = np.linalg.norm(grad)
            for w in freqs:
                if gn < 1e-10:
                    ig[(op, w)] = np.zeros(grad.shape[0])
                else:
                    ig[(op, w)] = np.matmul(precond[w], grad)

        return ig

    def setup_trials(self,
                     vectors,
                     td=None,
                     bger=np.array([]),
                     bung=np.array([]),
                     threshold=1e-10,
                     normalize=True):
        """Returns orthonormalized trial vectors. Takes set of vectors,
        preconditioner matrix, gerade and ungerade subspaces as input
        arguments.
        """

        trials = []
        for (op, w) in vectors:
            vec = np.array(vectors[(op, w)])

            # preconditioning trials:

            if td is not None:
                v = np.matmul(td[w], vec)
            else:
                v = vec
            if np.linalg.norm(v) > threshold:
                trials.append(v)

        new_trials = np.array(trials).T

        # decomposing the full space trial vectors...

        new_realger, new_realung, new_imagung, new_imagger = self.decomp_trials(
            new_trials)

        # ...and assembling gerade and ungerade subspaces

        new_ger = self.assemble_subsp(new_realger, new_imagger)
        new_ung = self.assemble_subsp(new_realung, new_imagung)

        # orthogonalizing new trial vectors against existing ones

        if bger.any():
            new_ger = new_ger - np.matmul(np.matmul(bger, bger.T), new_ger)
        if bung.any():
            new_ung = new_ung - np.matmul(np.matmul(bung, bung.T), new_ung)

        # normalizing new trial vectors

        if new_ger.any() and normalize:

            # removing linear dependencies in gerade trials:

            t_ger = self.rm_lin_depend(new_ger)
            new_ger = np.matmul(new_ger, t_ger)

            # normalizing gerade trials

            new_ger = self.normalize(new_ger)

        if new_ung.any() and normalize:

            # removing linear dependencies in ungerade trials:

            t_ung = self.rm_lin_depend(new_ung)
            new_ung = np.matmul(new_ung, t_ung)

            # normalizing ungerade trials

            new_ung = self.normalize(new_ung)

        return new_ger, new_ung

    def compute(self,
                mol_orbs,
                task,
                scf_tensors,
                ops='z',
                freqs=(
                    0.0,
                    0.05,
                    0.1,
                    0.15,
                    0.2,
                ),
                d=0.004556335294880438,
                maxit=150,
                threshold=1.0e-6):
        """Solves for the approximate response vector iteratively
        while checking the residuals for convergence.

        Input arguments are the calculation parameters as operators,
        frequencies, damping parameter, maximim number of iterations and
        convergence threshold.
        """

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, task.molecule, task.ao_basis)

        e2x_drv = LinearResponseMatrixVectorDriver(self.comm)

        if self.rank == mpi_master():
            nocc = task.molecule.number_of_alpha_electrons()
        else:
            nocc = None

        trials_info = {'bger': False, 'bung': False}

        if self.rank == mpi_master():

            # calling the gradients

            v1 = {op: v for op, v in zip(ops, self.get_rhs(mol_orbs, task, scf_tensors, ops))}

            # creating the preconditioner matrix

            precond = {w: self.get_precond(mol_orbs, task, w, d) for w in freqs}

            # spawning initial trial vectors

            igs = self.initial_guess(mol_orbs, task, scf_tensors, d, ops, freqs, precond)
            bger, bung = self.setup_trials(igs)

            # creating sigma and rho linear transformations

            if bger.any():
                trials_info['bger'] = True
                bger = self.orthogonalize_gram_schmidt(bger)
                bger = self.normalize(bger)

            if bung.any():
                trials_info['bung'] = True
                bung = self.orthogonalize_gram_schmidt(bung)
                bung = self.normalize(bung)

        else:
            bger = None
            bung = None

        trials_info = self.comm.bcast(trials_info, root=mpi_master())

        if trials_info['bger']:
            e2bger = e2x_drv.e2n(bger, scf_tensors, screening, task.molecule, task.ao_basis)
            if self.rank == mpi_master():
                s2bung = e2x_drv.s2n(bger, scf_tensors, nocc)

        if trials_info['bung']:
            e2bung = e2x_drv.e2n(bung, scf_tensors, screening, task.molecule, task.ao_basis)
            if self.rank == mpi_master():
                s2bger = e2x_drv.s2n(bung, scf_tensors, nocc)

        solutions = {}
        residuals = {}
        relative_residual_norm = {}

        for i in range(maxit):

            if self.rank == mpi_master():

                for op, w in igs:
                    grad = v1[op]

                    gradger, gradung = self.decomp_sym(grad)
                    full_size = gradger.shape[0]

                    # projections onto gerade and ungerade subspaces:

                    if bger.any():
                        g_realger = np.matmul(bger.T, gradger.real)
                        g_imagger = np.matmul(bger.T, gradger.imag)

                        e2gg = np.matmul(bger.T, e2bger)

                        ntrials_ger = bger.shape[1]

                    else:
                        ntrials_ger = 0

                    if bung.any():
                        g_realung = np.matmul(bung.T, gradung.real)
                        g_imagung = np.matmul(bung.T, gradung.imag)

                        e2uu = np.matmul(bung.T, e2bung)

                        if bger.any():
                            s2ug = np.matmul(bung.T, s2bung)

                        ntrials_ung = bung.shape[1]

                    else:
                        ntrials_ung = 0

                    # creating gradient and matrix for linear equation

                    print(ntrials_ger, 'gerade trial vectors')
                    print(ntrials_ung, 'ungerade trial vectors')
                    size = 2 * (ntrials_ger + ntrials_ung)

                    # gradient

                    g = np.zeros(size)

                    for pos in range(ntrials_ger):
                        g[pos] = g_realger[pos]
                        g[-pos - 1] = -g_imagger[-pos - 1]

                    for pos in range(ntrials_ung):
                        g[pos + ntrials_ger] = g_realung[pos]
                        g[-(pos + ntrials_ger) - 1] = -g_imagung[-pos - 1]

                    # matrix

                    mat = np.zeros((size, size))

                    # filling E2gg

                    for row in range(ntrials_ger):
                        for col in range(ntrials_ger):
                            mat[row, col] = e2gg[row, col]
                            mat[-row - 1, -col - 1] = -e2gg[-row - 1, -col - 1]

                    # filling E2uu

                    for row in range(ntrials_ung):
                        for col in range(ntrials_ung):
                            mat[(row + ntrials_ger),
                                (col + ntrials_ger)] = e2uu[row, col]
                            mat[-(row + ntrials_ger) - 1, -(col + ntrials_ger) -
                                1] = -e2uu[-row - 1, -col - 1]

                    for row in range(ntrials_ung):
                        for col in range(ntrials_ger):

                            # filling S2ug

                            mat[(row + ntrials_ger), col] = -w * s2ug[row, col]
                            mat[-(row + ntrials_ger) -
                                1, col] = d * s2ug[-row - 1, col]
                            mat[(row + ntrials_ger), -col -
                                1] = d * s2ug[row, -col - 1]
                            mat[-(row + ntrials_ger) - 1, -col -
                                1] = w * s2ug[-row - 1, -col - 1]

                            # filling S2ug.T (interchanging of row and col)

                            mat[col, (row + ntrials_ger)] = -w * s2ug[row, col]
                            mat[col, -(row + ntrials_ger) -
                                1] = d * s2ug[-row - 1, col]
                            mat[-col - 1,
                                (row + ntrials_ger)] = d * s2ug[row, -col - 1]
                            mat[-col - 1, -(row + ntrials_ger) -
                                1] = w * s2ug[-row - 1, -col - 1]

                    # solving matrix equation

                    c = np.linalg.solve(mat, g)

                    # extracting the 4 components of c...

                    c_realger, c_imagger = np.zeros(ntrials_ger), np.zeros(
                        ntrials_ger)
                    c_realung, c_imagung = np.zeros(ntrials_ung), np.zeros(
                        ntrials_ung)

                    for pos in range(ntrials_ger):
                        c_realger[pos] = c[pos]
                        c_imagger[-pos - 1] = c[-pos - 1]

                    for pos in range(ntrials_ung):
                        c_realung[pos] = c[pos + ntrials_ger]
                        c_imagung[-pos - 1] = c[-(pos + ntrials_ger) - 1]

                    # ...and projecting them onto respective subspace

                    x_realger = np.matmul(bger, c_realger)
                    x_imagger = np.matmul(bger, c_imagger)
                    x_realung = np.matmul(bung, c_realung)
                    x_imagung = np.matmul(bung, c_imagung)

                    # composing response vector

                    x_real = x_realger + x_realung
                    x_imag = x_imagung + x_imagger
                    x = np.zeros(len(x_real), dtype=complex)

                    for pos in range(len(x_real)):
                        x[pos] = complex(x_real[pos], x_imag[pos])

                    solutions[(op, w)] = x

                    # composing E2 and S2 matrices projected onto solution subspace

                    if bger.any():
                        e2realger = np.matmul(e2bger, c_realger)
                        e2imagger = np.matmul(e2bger, c_imagger)
                        s2realger = np.matmul(s2bung, c_realger)
                        s2imagger = np.matmul(s2bung, c_imagger)

                    else:
                        e2realger = np.zeros(full_size)
                        e2imagger = np.zeros(full_size)
                        s2realger = np.zeros(full_size)
                        s2imagger = np.zeros(full_size)

                    if bung.any():
                        e2realung = np.matmul(e2bung, c_realung)
                        e2imagung = np.matmul(e2bung, c_imagung)
                        s2realung = np.matmul(s2bger, c_realung)
                        s2imagung = np.matmul(s2bger, c_imagung)

                    else:
                        e2realung = np.zeros(full_size)
                        e2imagung = np.zeros(full_size)
                        s2realung = np.zeros(full_size)
                        s2imagung = np.zroes(full_size)

                    # calculating the residual components

                    r_realger = (e2realger - w * s2realung + d * s2imagung -
                                 gradger.real)
                    r_realung = (e2realung - w * s2realger + d * s2imagger -
                                 gradung.real)
                    r_imagung = (-e2imagung + w * s2imagger + d * s2realger +
                                 gradung.imag)
                    r_imagger = (-e2imagger + w * s2imagung + d * s2realung +
                                 gradger.imag)

                    # composing total residual

                    r_real = r_realger + r_realung
                    r_imag = r_imagung + r_imagger
                    r = np.zeros(len(r_real), dtype=complex)

                    for pos in range(len(r_real)):
                        r[pos] = complex(r_real[pos], r_imag[pos])

                    residuals[(op, w)] = np.array(
                        [r_realger, r_realung, r_imagung, r_imagger]).flatten()

                    n = solutions[(op, w)]

                    # calculating relative residual norm for convergence check

                    nv = np.matmul(n, grad)
                    rn = np.linalg.norm(r)
                    nn = np.linalg.norm(n)
                    if nn != 0:
                        relative_residual_norm[(op, w)] = rn / nn
                    else:
                        relative_residual_norm[(op, w)] = 0
                    print(f"{i+1} <<{op};{op}>>({w})={-nv:.10f} rn={rn:.5e}")

                print()

                # checking for convergence

                max_residual = max(relative_residual_norm.values())

                if max_residual < threshold:
                    self.is_converged = True
                    #print('Converged')
                    #break

            self.is_converged = self.comm.bcast(self.is_converged,
                                                root=mpi_master())

            if self.is_converged:
                break

            trials_info = {'new_trials_ger': False, 'new_trials_ung': False}

            # spawning new trial vectors from residuals

            if self.rank == mpi_master():

                new_trials_ger, new_trials_ung = self.setup_trials(residuals,
                                                                   td=precond,
                                                                   bger=bger,
                                                                   bung=bung)

                # creating new sigma and rho linear transformations

                if new_trials_ger.any():
                    trials_info['new_trials_ger'] = True

                    bger = np.append(bger, new_trials_ger, axis=1)

                    bger = self.orthogonalize_gram_schmidt(bger)
                    bger = self.normalize(bger)

                if new_trials_ung.any():
                    trials_info['new_trials_ung'] = True

                    bung = np.append(bung, new_trials_ung, axis=1)

                    bung = self.orthogonalize_gram_schmidt(bung)
                    bung = self.normalize(bung)

            else:
                new_trials_ger = None
                new_trials_ung = None

            trials_info = self.comm.bcast(trials_info, root=mpi_master())

            if trials_info['new_trials_ger']:
                new_e2bger = e2x_drv.e2n(new_trials_ger, scf_tensors, screening, task.molecule, task.ao_basis)
                if self.rank == mpi_master():
                    new_s2bung = e2x_drv.s2n(new_trials_ger, scf_tensors, nocc)
                    e2bger = np.append(e2bger, new_e2bger, axis=1)
                    s2bung = np.append(s2bung, new_s2bung, axis=1)

            if trials_info['new_trials_ung']:
                new_e2bung = e2x_drv.e2n(new_trials_ung, scf_tensors, screening, task.molecule, task.ao_basis)
                if self.rank == mpi_master():
                    new_s2bger = e2x_drv.s2n(new_trials_ung, scf_tensors, nocc)
                    e2bung = np.append(e2bung, new_e2bung, axis=1)
                    s2bger = np.append(s2bger, new_s2bger, axis=1)

        return solutions
