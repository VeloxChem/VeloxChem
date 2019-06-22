import numpy as np
import time as tm

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import mpi_master
from .lrmatvecdriver import LinearResponseMatrixVectorDriver
from .lrmatvecdriver import remove_linear_dependence
from .lrmatvecdriver import construct_ed_sd
from .lrmatvecdriver import lrvec2mat
from .lrmatvecdriver import get_rhs
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_scheme
from .qqscheme import get_qq_type
from .inputparser import parse_frequencies


class ComplexResponse:
    """Implements the complex linear response solver.

    Implements the complex linear response solver.

    Attributes
    ----------
    a_operator
        The A operator
    a_components
        Cartesian components of the A operator
    b_operator
        The B operator
    b_components
        Cartesian components of the B operator
    frequencies
        The frequencies.
    damping
        The damping parameter.
    qq_type
        The electron repulsion integrals screening scheme.
    eri_thresh
        The electron repulsion integrals screening threshold.
    max_iter
        The maximum number of solver iterations.
    conv_thresh
        The convergence threshold for the solver.
    cur_iter
        Index of the current iteration.
    is_converged
        The flag for convergence.
    small_thresh
        The norm threshold for a vector to be considered a zero vector.
    lindep_thresh
        The threshold for removing linear dependence in the trial vectors.
    comm
        The MPI communicator.
    rank
        The MPI rank.
    nodes
        Number of MPI processes.
    ostream
        The output stream.
    timing
        The flag for printing timing information.
    profiling
        The flag for printing profiling information.
    """

    def __init__(self, comm, ostream):
        """Initializes complex linear response solver.

        Initializes complex linear response solver to default setup.

        Parameters
        ----------
        comm
            The MPI communicator.
        ostream
            The output stream.
        """

        self.a_operator = 'dipole'
        self.a_components = 'xyz'
        self.b_operator = 'dipole'
        self.b_components = 'xyz'

        self.frequencies = (0.1,)
        self.damping = 0.004556335294880438

        self.qq_type = 'QQ_DEN'
        self.eri_thresh = 1.0e-12

        self.max_iter = 150
        self.conv_thresh = 1.0e-5

        self.cur_iter = 0
        self.is_converged = False
        self.small_thresh = 1.0e-10
        self.lindep_thresh = 1.0e-6

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream

        self.timing = False
        self.profiling = False

    def update_settings(self, settings):
        """Updates settings in complex linear response solver.

        Updates settings in complex liner response solver.

        Parameters
        ----------
        settings
            The settings dictionary.
        """

        if 'frequencies' in settings:
            self.frequencies = parse_frequencies(settings['frequencies'])
        if 'damping' in settings:
            self.damping = float(settings['damping'])

        if 'lindep_thresh' in settings:
            self.lindep_thresh = float(settings['lindep_thresh'])

        if 'conv_thresh' in settings:
            self.conv_thresh = float(settings['conv_thresh'])
        if 'max_iter' in settings:
            self.max_iter = int(settings['max_iter'])

        if 'eri_thresh' in settings:
            self.eri_thresh = float(settings['eri_thresh'])
        if 'qq_type' in settings:
            self.qq_type = settings['qq_type'].upper()

        if 'timing' in settings:
            key = settings['timing'].lower()
            self.timing = True if key in ['yes', 'y'] else False
        if 'profiling' in settings:
            key = settings['profiling'].lower()
            self.profiling = True if key in ['yes', 'y'] else False

    def paired(self, v_xy):
        """Computes paired trial vector.

        Computes paired trial vector.

        Parameters
        ----------
        v_xy
            The trial vector.

        Returns
        -------
            The paired trial vector.
        """

        v_yx = v_xy.copy()
        half_rows = v_xy.shape[0] // 2
        v_yx[:half_rows] = v_xy[half_rows:]
        v_yx[half_rows:] = v_xy[:half_rows]

        return v_yx

    def decomp_trials(self, vecs):
        """Decomposes trial vectors into their 4 respective parts.

        Decomposes trial vectors into their 4 respective parts (real gerade,
        real ungerade, imaginary gerade, and imaginary ungerad).

        Parameters
        ----------
        vecs
            The trial vectors.

        Returns
        -------
            A tuple containing respective parts of the trial vectors.
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

    def assemble_subsp(self, realvec, imagvec):
        """Assembles subspace out of real and imaginary parts of trials.

        Assembles subspace out of real and imaginary parts of trials,
        if their norm exceeds a certain threshold (zero vectors shouldn't
        be added).

        Parameters
        ----------
        realvec
            The real part of trial vectors.
        imagvec
            The imaginary part of trial vectors.

        Returns
        -------
            The assembled trial vectors.
        """

        space = []
        if len(realvec.shape) != 1:
            for vec in range(len(realvec[0, :])):
                if np.linalg.norm(realvec[:, vec]) > self.small_thresh:
                    space.append(realvec[:, vec])
                if np.linalg.norm(imagvec[:, vec]) > self.small_thresh:
                    space.append(imagvec[:, vec])

        else:
            if np.linalg.norm(realvec) > self.small_thresh:
                space.append(realvec)
            if np.linalg.norm(imagvec) > self.small_thresh:
                space.append(imagvec)

        return np.array(space).T

    def decomp_sym(self, vecs):
        """Decomposes gradient into gerade and ungerade parts.

        Decomposes gradient into gerade and ungerade parts.

        Parameters
        ----------
        vecs
            The trial vectors.

        Returns
        -------
            A tuple containing gerade and ungerade parts of vectors.
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

    def orthogonalize_gram_schmidt(self, tvecs):
        """Applies modified Gram Schmidt orthogonalization to trial vectors.

        Applies modified Gram Schmidt orthogonalization to trial vectors.

        Parameters
        ----------
        tvecs
            The trial vectors.

        Returns
        -------
            The orthogonalized trial vectors.
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
        """Normalizes vectors.

        Normalizes vectors by dividing by vector norm.

        Parameters
        ----------
            The vectors.
        Retruns
        -------
            The normalized vectors.
        """

        if len(vecs.shape) != 1:
            for vec in range(vecs.shape[1]):
                invnorm = 1.0 / np.linalg.norm(vecs[:, vec])
                vecs[:, vec] *= invnorm
        else:
            invnorm = 1.0 / np.linalg.norm(vecs)
            vecs *= invnorm

        return vecs

    def get_precond(self, orb_ene, nocc, norb, w, d):
        """Constructs the preconditioner matrix.

        Constructs the preconditioner matrix.

        Parameters
        ----------
        orb_ene
            The orbital energies.
        nocc
            The number of doubly occupied orbitals.
        norb
            The number of orbitals.
        w
            The frequency.
        d
            The damping parameter.

        Returns
        -------
            The preconditioner matrix.
        """

        # spawning needed components

        ediag, sdiag = construct_ed_sd(orb_ene, nocc, norb)
        ediag_sq = ediag**2
        sdiag_sq = sdiag**2
        sdiag_fp = sdiag**4
        w_sq = w**2
        d_sq = d**2

        # constructing matrix block diagonals

        a_diag = ediag * (ediag_sq - (w_sq - d_sq) * sdiag_sq)
        b_diag = (w * sdiag) * (ediag_sq - (w_sq + d_sq) * sdiag_sq)
        c_diag = (d * sdiag) * (ediag_sq + (w_sq + d_sq) * sdiag_sq)
        d_diag = (2 * w * d * ediag) * sdiag_sq
        p_diag = 1.0 / ((ediag_sq - (w_sq - d_sq) * sdiag_sq)**2 +
                        (4 * w_sq * d_sq * sdiag_fp))

        pa_diag = p_diag * a_diag
        pb_diag = p_diag * b_diag
        pc_diag = p_diag * c_diag
        pd_diag = p_diag * d_diag

        precond = np.array([pa_diag, pb_diag, pc_diag, pd_diag])

        return precond

    def preconditioning(self, precond, v_in):
        """Creates trial vectors out of residuals and the preconditioner matrix.

        Creates trial vectors out of residuals and the preconditioner matrix.

        Parameters
        ----------
        precond
            The preconditioner matrix.
        v_in
            The input trial vectors.

        Returns
        -------
            The trail vectors after preconditioning.
        """

        pa, pb, pc, pd = precond[0], precond[1], precond[2], precond[3]

        v_in_rg, v_in_ru, v_in_iu, v_in_ig = self.decomp_trials(v_in)

        v_out_rg = pa * v_in_rg + pb * v_in_ru + pc * v_in_iu + pd * v_in_ig
        v_out_ru = pb * v_in_rg + pa * v_in_ru + pd * v_in_iu + pc * v_in_ig
        v_out_iu = pc * v_in_rg + pd * v_in_ru - pa * v_in_iu - pb * v_in_ig
        v_out_ig = pd * v_in_rg + pc * v_in_ru - pb * v_in_iu - pa * v_in_ig

        v_out = np.array([v_out_rg, v_out_ru, v_out_iu, v_out_ig]).flatten()

        return v_out

    def initial_guess(self, op_grads, d, freqs, precond):
        """Creating initial guess (un-orthonormalized trials) out of gradients.

        Creating initial guess (un-orthonormalized trials) out of gradients.

        Parameters
        ----------
        op_grads
            The dictionary containing operator components (key) and right-hand
            sides (values).
        d
            The damping parameter.
        freqs
            The frequencies.
        precond
            The preconditioner matrices.

        Returns
        -------
            The initial guess.
        """

        ig = {}
        for op, grad in op_grads.items():
            gradger, gradung = self.decomp_sym(grad)
            grad = np.array(
                [gradger.real, gradung.real, -gradung.imag,
                 -gradger.imag]).flatten()
            gn = np.linalg.norm(grad)
            for w in freqs:
                if gn < self.small_thresh:
                    ig[(op, w)] = np.zeros(grad.shape[0])
                else:
                    ig[(op, w)] = self.preconditioning(precond[w], grad)

        return ig

    def setup_trials(self,
                     vectors,
                     pre=None,
                     bger=np.array([]),
                     bung=np.array([]),
                     res_norm=None,
                     normalize=True):
        """Computes orthonormalized trial vectors.

        Computes orthonormalized trial vectors. Takes set of vectors,
        preconditioner matrix, gerade and ungerade subspaces as input
        arguments.

        Parameters
        ----------
        vectors
            The set of vectors.
        pre
            The preconditioner matrix.
        bger
            The gerade subspace.
        bung
            The ungerade subspace.
        res_norm
            The relative residual norm.
        normalize
            The flag for normalization.

        Returns
        -------
            The orthonormalized trial vectors.
        """

        trials = []
        for (op, w) in vectors:
            if res_norm is None or res_norm[(op, w)] > self.conv_thresh:
                vec = np.array(vectors[(op, w)])

                # preconditioning trials:

                if pre is not None:
                    v = self.preconditioning(pre[w], vec)
                else:
                    v = vec
                if np.linalg.norm(v) > self.small_thresh:
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
            new_ger = new_ger - np.matmul(bger, np.matmul(bger.T, new_ger))
        if bung.any():
            new_ung = new_ung - np.matmul(bung, np.matmul(bung.T, new_ung))

        # normalizing new trial vectors

        if new_ger.any() and normalize:

            # removing linear dependencies in gerade trials
            # and normalizing gerade trials

            new_ger = remove_linear_dependence(new_ger, self.lindep_thresh)
            new_ger = self.orthogonalize_gram_schmidt(new_ger)
            new_ger = self.normalize(new_ger)

        if new_ung.any() and normalize:

            # removing linear dependencies in ungerade trials:
            # and normalizing ungerade trials

            new_ung = remove_linear_dependence(new_ung, self.lindep_thresh)
            new_ung = self.orthogonalize_gram_schmidt(new_ung)
            new_ung = self.normalize(new_ung)

        return new_ger, new_ung

    def compute(self, molecule, basis, scf_tensors, b_rhs=None):
        """Solves for the response vector.

        Solves for the response vector iteratively while checking the residuals
        for convergence.

        Parameters
        ----------
        molecule
            The molecule.
        basis
            The AO basis.
        scf_tensors
            The dictionary of tensors from converged SCF wavefunction.
        b_rhs
            The right-hand side. If not provided, b_rhs will be computed for
            the B operator.

        Returns
        -------
            A dictionary containing properties, solutions, and kappas.
        """

        if self.profiling:
            import cProfile
            import pstats
            import io
            import os
            pr = cProfile.Profile()
            pr.enable()

        if self.timing:
            prep_t0 = tm.time()
            self.timing_dict = {}

        if self.rank == mpi_master():
            self.print_header()

        self.start_time = tm.time()

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'ComplexResponseSolver: not implemented for unrestricted case')

        d = self.damping
        freqs = self.frequencies

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        e2x_drv = LinearResponseMatrixVectorDriver(self.comm)

        if self.rank == mpi_master():
            nocc = molecule.number_of_alpha_electrons()
            norb = scf_tensors['C'].shape[1]
            orb_ene = scf_tensors['E']
        else:
            nocc = None

        if b_rhs is None:
            b_rhs = get_rhs(self.b_operator, self.b_components, molecule, basis,
                            scf_tensors, self.rank, self.comm)

        if self.rank == mpi_master():

            trials_info = {'bger': False, 'bung': False}

            # calling the gradients

            v1 = {op: v for op, v in zip(self.b_components, b_rhs)}

            # creating the preconditioner matrix

            precond_start_time = tm.time()
            precond = {
                w: self.get_precond(orb_ene, nocc, norb, w, d) for w in freqs
            }
            self.ostream.print_info(
                'Precondition {} created in {:.2f} sec.'.format(
                    'matrices' if len(freqs) > 1 else 'matrix',
                    tm.time() - precond_start_time))
            self.ostream.print_blank()
            self.ostream.flush()

            # spawning initial trial vectors

            igs = self.initial_guess(v1, d, freqs, precond)
            bger, bung = self.setup_trials(igs)

            assert_msg_critical(
                bger.any() or bung.any(),
                'ComplexResponseSolver: trial vectors are empty')

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
            trials_info = {}
            bger = None
            bung = None

        trials_info = self.comm.bcast(trials_info, root=mpi_master())

        if trials_info['bger']:
            e2bger = e2x_drv.e2n(bger, scf_tensors, screening, molecule, basis)
            if self.rank == mpi_master():
                s2bung = e2x_drv.s2n(bger, scf_tensors, nocc)

        if trials_info['bung']:
            e2bung = e2x_drv.e2n(bung, scf_tensors, screening, molecule, basis)
            if self.rank == mpi_master():
                s2bger = e2x_drv.s2n(bung, scf_tensors, nocc)

        solutions = {}
        residuals = {}
        relative_residual_norm = {}
        kappas = {}

        if self.timing:
            self.timing_dict['initial_guess'] = tm.time() - prep_t0
            self.timing_dict['reduced_space'] = []
            self.timing_dict['new_trials'] = []

        for iteration in range(self.max_iter):
            self.cur_iter = iteration

            if self.timing:
                red_space_t0 = tm.time()

            if self.rank == mpi_master():
                nvs = []

                for op, w in igs:
                    if iteration == 0 or (relative_residual_norm[(op, w)] >
                                          self.conv_thresh):
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

                        size = 2 * (ntrials_ger + ntrials_ung)

                        # gradient

                        g = np.zeros(size)

                        g[:ntrials_ger] = g_realger[:]
                        g[-ntrials_ger:] = -g_imagger[:]

                        g[ntrials_ger:ntrials_ger + ntrials_ung] = g_realung[:]
                        g[-ntrials_ger -
                          ntrials_ung:-ntrials_ger] = g_imagung[:]

                        # matrix

                        mat = np.zeros((size, size))

                        # filling E2gg

                        mat[:ntrials_ger, :ntrials_ger] = e2gg[:, :]

                        mat[-ntrials_ger:, -ntrials_ger:] = -e2gg[:, :]

                        # filling E2uu

                        mat[ntrials_ger:ntrials_ger +
                            ntrials_ung, ntrials_ger:ntrials_ger +
                            ntrials_ung] = e2uu[:, :]

                        mat[-ntrials_ger -
                            ntrials_ung:-ntrials_ger, -ntrials_ger -
                            ntrials_ung:-ntrials_ger] = -e2uu[:, :]

                        # filling S2ug

                        mat[ntrials_ger:ntrials_ger +
                            ntrials_ung, :ntrials_ger] = -w * s2ug[:, :]

                        mat[-ntrials_ger - ntrials_ung:-ntrials_ger, :
                            ntrials_ger] = d * s2ug[:, :]

                        mat[ntrials_ger:ntrials_ger +
                            ntrials_ung, -ntrials_ger:] = d * s2ug[:, :]

                        mat[-ntrials_ger - ntrials_ung:-ntrials_ger,
                            -ntrials_ger:] = w * s2ug[:, :]

                        # filling S2ug.T (interchanging of row and col)

                        mat[:ntrials_ger, ntrials_ger:ntrials_ger +
                            ntrials_ung] = -w * s2ug.T[:, :]

                        mat[:ntrials_ger, -ntrials_ger -
                            ntrials_ung:-ntrials_ger] = d * s2ug.T[:, :]

                        mat[-ntrials_ger:, ntrials_ger:ntrials_ger +
                            ntrials_ung] = d * s2ug.T[:, :]

                        mat[-ntrials_ger:, -ntrials_ger -
                            ntrials_ung:-ntrials_ger] = w * s2ug.T[:, :]

                        # solving matrix equation

                        c = np.linalg.solve(mat, g)

                        # extracting the 4 components of c...

                        c_realger, c_imagger = np.zeros(ntrials_ger), np.zeros(
                            ntrials_ger)
                        c_realung, c_imagung = np.zeros(ntrials_ung), np.zeros(
                            ntrials_ung)

                        c_realger[:] = c[:ntrials_ger]
                        c_imagger[:] = c[-ntrials_ger:]

                        c_realung[:] = c[ntrials_ger:ntrials_ger + ntrials_ung]
                        c_imagung[:] = c[-ntrials_ger -
                                         ntrials_ung:-ntrials_ger]

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

                        kappas[(op, w)] = (lrvec2mat(x.real, nocc, norb) +
                                           1j * lrvec2mat(x.imag, nocc, norb))

                        # composing E2 and S2 matrices projected onto solution
                        # subspace

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
                            s2imagung = np.zeros(full_size)

                        # calculating the residual components

                        r_realger = (e2realger - w * s2realung + d * s2imagung -
                                     gradger.real)
                        r_realung = (e2realung - w * s2realger + d * s2imagger -
                                     gradung.real)
                        r_imagung = (-e2imagung + w * s2imagger +
                                     d * s2realger + gradung.imag)
                        r_imagger = (-e2imagger + w * s2imagung +
                                     d * s2realung + gradger.imag)

                        # composing total residual

                        r_real = r_realger + r_realung
                        r_imag = r_imagung + r_imagger
                        r = np.zeros(len(r_real), dtype=complex)

                        for pos in range(len(r_real)):
                            r[pos] = complex(r_real[pos], r_imag[pos])

                        residuals[(op, w)] = np.array(
                            [r_realger, r_realung, r_imagung,
                             r_imagger]).flatten()

                        n = solutions[(op, w)]

                        # calculating relative residual norm
                        # for convergence check

                        nv = np.matmul(n, grad)
                        nvs.append((op, w, nv))

                        rn = np.linalg.norm(r)
                        nn = np.linalg.norm(n)
                        if nn != 0:
                            relative_residual_norm[(op, w)] = rn / nn
                        else:
                            relative_residual_norm[(op, w)] = 0

                # write to output

                self.ostream.print_info(
                    '{:d} gerade trial vectors'.format(ntrials_ger))
                self.ostream.print_info(
                    '{:d} ungerade trial vectors'.format(ntrials_ung))
                self.ostream.print_blank()

                self.print_iteration(relative_residual_norm, nvs)

            if self.timing:
                self.timing_dict['reduced_space'].append(tm.time() -
                                                         red_space_t0)

            # check convergence

            self.check_convergence(relative_residual_norm)

            if self.is_converged:
                break

            if self.timing:
                new_trials_t0 = tm.time()

            # spawning new trial vectors from residuals

            if self.rank == mpi_master():

                trials_info = {'new_trials_ger': False, 'new_trials_ung': False}

                new_trials_ger, new_trials_ung = self.setup_trials(
                    residuals,
                    pre=precond,
                    bger=bger,
                    bung=bung,
                    res_norm=relative_residual_norm)

                assert_msg_critical(
                    new_trials_ger.any() or new_trials_ung.any(),
                    'ComplexResponseSolver: unable to add new trial vectors')

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
                trials_info = {}
                new_trials_ger = None
                new_trials_ung = None

            trials_info = self.comm.bcast(trials_info, root=mpi_master())

            if trials_info['new_trials_ger']:
                new_e2bger = e2x_drv.e2n(new_trials_ger, scf_tensors, screening,
                                         molecule, basis)
                if self.rank == mpi_master():
                    new_s2bung = e2x_drv.s2n(new_trials_ger, scf_tensors, nocc)
                    e2bger = np.append(e2bger, new_e2bger, axis=1)
                    s2bung = np.append(s2bung, new_s2bung, axis=1)

            if trials_info['new_trials_ung']:
                new_e2bung = e2x_drv.e2n(new_trials_ung, scf_tensors, screening,
                                         molecule, basis)
                if self.rank == mpi_master():
                    new_s2bger = e2x_drv.s2n(new_trials_ung, scf_tensors, nocc)
                    e2bung = np.append(e2bung, new_e2bung, axis=1)
                    s2bger = np.append(s2bger, new_s2bger, axis=1)

            if self.timing:
                self.timing_dict['new_trials'].append(tm.time() - new_trials_t0)

        # converged?
        if self.rank == mpi_master():
            self.print_convergence()

            assert_msg_critical(self.is_converged,
                                'ComplexResponseSolver: failed to converge')

            if self.timing:
                self.print_timing()

        if self.profiling:
            pr.disable()
            s = io.StringIO()
            sortby = 'cumulative'
            ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
            ps.print_stats(20)
            if self.rank == mpi_master():
                for line in s.getvalue().split(os.linesep):
                    self.ostream.print_info(line)

        a_rhs = get_rhs(self.a_operator, self.a_components, molecule, basis,
                        scf_tensors, self.rank, self.comm)

        if self.rank == mpi_master():
            va = {op: v for op, v in zip(self.a_components, a_rhs)}
            props = {}
            for aop in self.a_components:
                for bop, w in solutions:
                    props[(aop, bop, w)] = -np.dot(va[aop], solutions[(bop, w)])
            self.print_properties(props)

            return {
                'properties': props,
                'solutions': solutions,
                'kappas': kappas,
            }
        else:
            return {}

    def check_convergence(self, relative_residual_norm):
        """Checks convergence.

        Checks convergence.

        Parameters
        ----------
        relative_residual_norm
            Relative residual norms.
        """

        if self.rank == mpi_master():
            max_residual = max(relative_residual_norm.values())
            if max_residual < self.conv_thresh:
                self.is_converged = True

        self.is_converged = self.comm.bcast(self.is_converged,
                                            root=mpi_master())

    def print_iteration(self, relative_residual_norm, nvs):
        """Prints information of the iteration.

        Prints information of the iteration.

        Parameters
        ----------
        relative_residual_norm
            Relative residual norms.
        nvs
            A list of tuples containing operator component, frequency, and
            property.
        """

        output_header = '*** Iteration:   {} '.format(self.cur_iter + 1)
        output_header += '* Residuals (Max,Min): '
        output_header += '{:.2e} and {:.2e}'.format(
            max(relative_residual_norm.values()),
            min(relative_residual_norm.values()))
        self.ostream.print_header(output_header.ljust(82))
        self.ostream.print_blank()

        output_header = 'Operator:  {} ({})'.format(self.b_operator,
                                                    self.b_components)
        self.ostream.print_header(output_header.ljust(82))
        self.ostream.print_blank()

        for op, freq, nv in nvs:
            ops_label = '<<{};{}>>_{:.4f}'.format(op, op, freq)
            rel_res = relative_residual_norm[(op, freq)]
            output_iter = '{:<15s}: {:15.8f} {:15.8f}j   '.format(
                ops_label, -nv.real, -nv.imag)
            output_iter += 'Residual Norm: {:.8f}'.format(rel_res)
            self.ostream.print_header(output_iter.ljust(82))
        self.ostream.print_blank()
        self.ostream.flush()

    def print_header(self):
        """Prints complex linear response solver setup header.

        Prints complex linear response solver setup header to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header("Complex Response Driver Setup")
        self.ostream.print_header(31 * "=")
        self.ostream.print_blank()

        str_width = 60

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

    def print_convergence(self):
        """Prints information after convergence.
        """

        output_conv = '*** '
        if self.is_converged:
            output_conv += 'Complex response converged'
        else:
            output_conv += 'Complex response NOT converged'
        output_conv += ' in {:d} iterations. '.format(self.cur_iter + 1)
        output_conv += 'Time: {:.2f} sec'.format(tm.time() - self.start_time)
        self.ostream.print_header(output_conv.ljust(82))
        self.ostream.print_blank()

    def print_properties(self, props):
        """Prints properties.

        Prints properties.

        Parameters
        ----------
        props
            The dictionary of properties.
        """

        for w in self.frequencies:
            w_str = '{}, {}, w={:.4f}'.format(self.a_operator, self.b_operator,
                                              w)
            self.ostream.print_header(w_str.ljust(82))
            self.ostream.print_header(('-' * len(w_str)).ljust(82))
            for a in self.a_components:
                for b in self.b_components:
                    ops_label = '<<{};{}>>_{:.4f}'.format(a, b, w)
                    output_alpha = '{:<15s} {:15.8f} {:15.8f}j'.format(
                        ops_label, props[(a, b, w)].real, props[(a, b, w)].imag)
                    self.ostream.print_header(output_alpha.ljust(82))
            self.ostream.print_blank()

    def print_timing(self):
        """Prints timing.

        Prints timing for the complex linear response solver.
        """

        valstr = 'Timing:'
        self.ostream.print_header(valstr.ljust(82))
        self.ostream.print_header(('-' * len(valstr)).ljust(82))

        valstr = '{:<22s}: {:12.2f} sec'.format(
            'Initial guess', self.timing_dict['initial_guess'])
        self.ostream.print_header(valstr.ljust(82))

        for i, (a, b) in enumerate(
                zip(self.timing_dict['reduced_space'],
                    self.timing_dict['new_trials'])):
            valstr = '{:<22s}: {:12.2f} sec'.format(
                'Iteration {:d}'.format(i + 1), a + b)
            self.ostream.print_header(valstr.ljust(82))

        valstr = '{:<22s}: {:12.2f} sec'.format(
            'Last iteration', self.timing_dict['reduced_space'][-1])
        self.ostream.print_header(valstr.ljust(82))
        self.ostream.print_blank()
