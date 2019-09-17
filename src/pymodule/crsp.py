import numpy as np
import time as tm

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import mpi_master
from .veloxchemlib import GridDriver
from .veloxchemlib import MolecularGrid
from .veloxchemlib import XCFunctional
from .veloxchemlib import denmat
from .veloxchemlib import parse_xc_func
from .aodensitymatrix import AODensityMatrix
from .lrmatvecdriver import LinearResponseMatrixVectorDriver
from .lrmatvecdriver import remove_linear_dependence
from .lrmatvecdriver import orthogonalize_gram_schmidt
from .lrmatvecdriver import normalize
from .lrmatvecdriver import construct_ed_sd
from .lrmatvecdriver import lrvec2mat
from .lrmatvecdriver import get_rhs
from .lrmatvecdriver import read_rsp_hdf5
from .lrmatvecdriver import write_rsp_hdf5
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_scheme
from .qqscheme import get_qq_type
from .inputparser import parse_frequencies


class ComplexResponse:
    """
    Implements the complex linear response solver.

    :param a_operator:
        The A operator
    :param a_components:
        Cartesian components of the A operator
    :param b_operator:
        The B operator
    :param b_components:
        Cartesian components of the B operator
    :param frequencies:
        The frequencies.
    :param damping:
        The damping parameter.
    :param qq_type:
        The electron repulsion integrals screening scheme.
    :param eri_thresh:
        The electron repulsion integrals screening threshold.
    :param dft:
        The flag for running DFT.
    :param grid_level:
        The accuracy level of DFT grid.
    :param xcfun:
        The XC functional.
    :param molgrid:
        The molecular grid.
    :param gs_density:
        The ground state density matrix.
    :param max_iter:
        The maximum number of solver iterations.
    :param conv_thresh:
        The convergence threshold for the solver.
    :param cur_iter:
        Index of the current iteration.
    :param is_converged:
        The flag for convergence.
    :param small_thresh:
        The norm threshold for a vector to be considered a zero vector.
    :param lindep_thresh:
        The threshold for removing linear dependence in the trial vectors.
    :param comm:
        The MPI communicator.
    :param rank:
        The MPI rank.
    :param nodes:
        Number of MPI processes.
    :param ostream:
        The output stream.
    :param restart:
        The flag for restarting from checkpoint file.
    :param checkpoint_file:
        The name of checkpoint file.
    :param timing:
        The flag for printing timing information.
    :param profiling:
        The flag for printing profiling information.
    """

    def __init__(self, comm, ostream):
        """
        Initializes complex linear response solver to default setup.

        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
        """

        self.a_operator = 'dipole'
        self.a_components = 'xyz'
        self.b_operator = 'dipole'
        self.b_components = 'xyz'

        self.frequencies = (0.1,)
        self.damping = 0.004556335294880438

        self.nonlinear = False

        self.qq_type = 'QQ_DEN'
        self.eri_thresh = 1.0e-15

        self.dft = False
        self.grid_level = 4
        self.xcfun = XCFunctional()
        self.molgrid = MolecularGrid()
        self.gs_density = AODensityMatrix()

        self.max_iter = 150
        self.conv_thresh = 1.0e-4

        self.cur_iter = 0
        self.is_converged = False
        self.small_thresh = 1.0e-10
        self.lindep_thresh = 1.0e-6

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream

        self.restart = True
        self.checkpoint_file = None

        self.timing = False
        self.profiling = False

    def update_settings(self, rsp_dict, method_dict={}):
        """
        Updates response and method settings in complex liner response solver.

        :param rsp_dict:
            The dictionary of response dict.
        :param method_dict:
            The dictionary of method rsp_dict.
        """

        if 'a_operator' in rsp_dict:
            self.a_operator = rsp_dict['a_operator'].lower()
        if 'a_components' in rsp_dict:
            self.a_components = rsp_dict['a_components'].lower()
        if 'b_operator' in rsp_dict:
            self.b_operator = rsp_dict['b_operator'].lower()
        if 'b_components' in rsp_dict:
            self.b_components = rsp_dict['b_components'].lower()

        if 'frequencies' in rsp_dict:
            self.frequencies = parse_frequencies(rsp_dict['frequencies'])
        if 'damping' in rsp_dict:
            self.damping = float(rsp_dict['damping'])

        if 'lindep_thresh' in rsp_dict:
            self.lindep_thresh = float(rsp_dict['lindep_thresh'])

        if 'conv_thresh' in rsp_dict:
            self.conv_thresh = float(rsp_dict['conv_thresh'])
        if 'max_iter' in rsp_dict:
            self.max_iter = int(rsp_dict['max_iter'])

        if 'eri_thresh' in rsp_dict:
            self.eri_thresh = float(rsp_dict['eri_thresh'])
        if 'qq_type' in rsp_dict:
            self.qq_type = rsp_dict['qq_type'].upper()

        if 'restart' in rsp_dict:
            key = rsp_dict['restart'].lower()
            self.restart = True if key == 'yes' else False
        if 'checkpoint_file' in rsp_dict:
            self.checkpoint_file = rsp_dict['checkpoint_file']

        if 'timing' in rsp_dict:
            key = rsp_dict['timing'].lower()
            self.timing = True if key in ['yes', 'y'] else False
        if 'profiling' in rsp_dict:
            key = rsp_dict['profiling'].lower()
            self.profiling = True if key in ['yes', 'y'] else False

        if 'dft' in method_dict:
            key = method_dict['dft'].lower()
            self.dft = True if key == 'yes' else False
        if 'grid_level' in method_dict:
            self.grid_level = int(method_dict['grid_level'])
        if 'xcfun' in method_dict:
            if 'dft' not in method_dict:
                self.dft = True
            self.xcfun = parse_xc_func(method_dict['xcfun'].upper())
            assert_msg_critical(not self.xcfun.is_undefined(),
                                'Undefined XC functional')

    def paired(self, v_xy):
        """
        Computes paired trial vector.

        :param v_xy:
            The trial vector.

        :return:
            The paired trial vector.
        """

        v_yx = v_xy.copy()
        half_rows = v_xy.shape[0] // 2
        v_yx[:half_rows] = v_xy[half_rows:]
        v_yx[half_rows:] = v_xy[:half_rows]

        return v_yx

    def decomp_trials(self, vecs):
        """
        Decomposes trial vectors into their 4 respective parts (real gerade,
        real ungerade, imaginary gerade, and imaginary ungerad).

        :param vecs:
            The trial vectors.

        :return:
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
        """
        Assembles subspace out of real and imaginary parts of trials,
        if their norm exceeds a certain threshold (zero vectors shouldn't
        be added).

        :param realvec:
            The real part of trial vectors.
        :param imagvec:
            The imaginary part of trial vectors.

        :return:
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
        """
        Decomposes gradient into gerade and ungerade parts.

        :param vecs:
            The trial vectors.

        :return:
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

    def discard_lower_halfs(self, vecs):
        """
        Discards lower halfs of vectors sent in.

        :param vecs:
            The trial vectors.

        :return:
            The upper halfs of the vectors sent in.
        """

        half = vecs.shape[0] // 2
        if len(vecs.shape) != 1:
            upperhalf = []
            for vec in range(len(vecs[0, :])):
                upperhalf.append(vecs[:half, vec])
        else:
            upperhalf = vecs[:half]

        return np.array(upperhalf)


    def get_precond(self, orb_ene, nocc, norb, w, d):
        """
        Constructs the preconditioner matrix.

        :param orb_ene:
            The orbital energies.
        :param nocc:
            The number of doubly occupied orbitals.
        :param norb:
            The number of orbitals.
        :param w:
            The frequency.
        :param d:
            The damping parameter.

        :return:
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
        """
        Creates trial vectors out of residuals and the preconditioner matrix.

        :param precond:
            The preconditioner matrix.
        :param v_in:
            The input trial vectors.

        :return:
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
        """
        Creating initial guess (un-orthonormalized trials) out of gradients.

        :param op_grads:
            The dictionary containing operator components (key) and right-hand
            sides (values).
        :param d:
            The damping parameter.
        :param freqs:
            The frequencies.
        :param precond:
            The preconditioner matrices.

        :return:
            The initial guess.
        """

        ig = {}
        for op, grad in op_grads.items():
            gradger, gradung = self.decomp_sym(grad)
            grad = np.array(
                [gradger.real, gradung.real, -gradung.imag,
                 -gradger.imag]).flatten()
            gn = np.linalg.norm(grad)
            if self.nonlinear:
                if gn < self.small_thresh:
                    ig[op] = np.zeros(grad.shape[0])
                else:
                    fullsize = self.preconditioning(precond[op[1]], grad)
                    ig[op] = self.discard_lower_halfs(np.array(
                             self.decomp_trials(fullsize)).T).flatten()
            else:
                for w in freqs:
                    if gn < self.small_thresh:
                        ig[(op, w)] = np.zeros(grad.shape[0])
                    else:
                        fullsize = self.preconditioning(precond[w], grad)
                        ig[(op, w)] = self.discard_lower_halfs(np.array(
                            self.decomp_trials(fullsize)).T).flatten()

        return ig

    def setup_trials(self,
                     vectors,
                     pre=None,
                     bger=None,
                     bung=None,
                     res_norm=None,
                     renormalize=True):
        """
        Computes orthonormalized trial vectors. Takes set of vectors,
        preconditioner matrix, gerade and ungerade subspaces as input
        arguments.

        :param vectors:
            The set of vectors.
        :param pre:
            The preconditioner matrix.
        :param bger:
            The gerade subspace.
        :param bung:
            The ungerade subspace.
        :param res_norm:
            The relative residual norm.
        :param renormalize:
            The flag for normalization.

        :return:
            The orthonormalized gerade and ungerade trial vectors.
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

        if bger is not None and bger.any():
            new_ger = new_ger - 2 * np.matmul(bger, np.matmul(bger.T, new_ger))

        if bung is not None and bung.any():
            new_ung = new_ung - 2 * np.matmul(bung, np.matmul(bung.T, new_ung))

        # normalizing new trial vectors

        if new_ger.any() and renormalize:

            # removing linear dependencies in gerade trials
            # and normalizing gerade trials

            new_ger = remove_linear_dependence(new_ger, self.lindep_thresh)
            new_ger = orthogonalize_gram_schmidt(new_ger)
            new_ger = normalize(new_ger)

        if new_ung.any() and renormalize:

            # removing linear dependencies in ungerade trials:
            # and normalizing ungerade trials

            new_ung = remove_linear_dependence(new_ung, self.lindep_thresh)
            new_ung = orthogonalize_gram_schmidt(new_ung)
            new_ung = normalize(new_ung)

        return new_ger, new_ung

    def compute(self, molecule, basis, scf_tensors, b_rhs=None):
        """
        Solves for the response vector iteratively while checking the residuals
        for convergence.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param b_rhs:
            The right-hand side. If not provided, b_rhs will be computed for
            the B operator.

        :return:
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
            self.timing_dict = {
                'reduced_space': [0.0],
                'ortho_norm': [0.0],
                'fock_build': [0.0],
            }
            timing_t0 = tm.time()

        if self.rank == mpi_master() and (not self.nonlinear):
            self.print_header()

        self.start_time = tm.time()

        # generate integration grid
        if self.dft:
            grid_drv = GridDriver(self.comm)
            grid_drv.set_level(self.grid_level)

            grid_t0 = tm.time()
            self.molgrid = grid_drv.generate(molecule)
            n_grid_points = self.molgrid.number_of_points()
            self.ostream.print_info(
                'Molecular grid with {0:d} points generated in {1:.2f} sec.'.
                format(n_grid_points,
                       tm.time() - grid_t0))
            self.ostream.print_blank()

            if self.rank == mpi_master():
                self.gs_density = AODensityMatrix([scf_tensors['D'][0]],
                                                  denmat.rest)
            self.gs_density.broadcast(self.rank, self.comm)

        if self.dft:
            dft_func_label = self.xcfun.get_func_label().upper()
        else:
            dft_func_label = 'HF'

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'ComplexResponseSolver: not implemented for unrestricted case')

        d = self.damping
        if b_rhs:
            freqs = [i[1] for i in b_rhs.keys()]
        else:
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
        else:
            v1 = b_rhs
            self.nonlinear = True

        if self.rank == mpi_master():

            # calling the gradients

            if not self.nonlinear:
                v1 = {op: v for op, v in zip(self.b_components, b_rhs)}
                op_freq_keys = [(op, w) for op in v1 for w in freqs]
            else:
                op_freq_keys = [op for op in v1]

            # creating the preconditioner matrix

            precond = {
                w: self.get_precond(orb_ene, nocc, norb, w, d) for w in freqs
            }

        # read initial guess from restart file

        if self.restart:
            if self.rank == mpi_master():
                bger, bung, e2bger, e2bung = read_rsp_hdf5(
                    self.checkpoint_file,
                    ['CLR_bger', 'CLR_bung', 'CLR_e2bger', 'CLR_e2bung'],
                    molecule.nuclear_repulsion_energy(),
                    molecule.elem_ids_to_numpy(), basis.get_label(),
                    dft_func_label, self.ostream)
                self.restart = (bger is not None and bung is not None and
                                e2bger is not None and e2bung is not None)
            self.restart = self.comm.bcast(self.restart, root=mpi_master())

        # generate initial guess from scratch

        if not self.restart:
            if self.rank == mpi_master():

                # spawning initial trial vectors
                igs = self.initial_guess(v1, d, freqs, precond)
                bger, bung = self.setup_trials(igs)

                assert_msg_critical(
                    bger.any() or bung.any(),
                    'ComplexResponseSolver: trial vectors are empty')

                if not bger.any():
                    bger = np.zeros((bung.shape[0], 0))
                if not bung.any():
                    bung = np.zeros((bger.shape[0], 0))

                # creating sigma and rho linear transformations

            if self.timing:
                self.timing_dict['ortho_norm'][0] += tm.time() - timing_t0
                timing_t0 = tm.time()

            half_bger = None
            half_bung = None
            if self.rank == mpi_master():
                if bger is not None:
                    half_bger = bger#[:bger.shape[0] // 2]
                if bung is not None:
                    half_bung = bung#[:bung.shape[0] // 2]

            half_e2bger, half_e2bung = e2x_drv.e2n_half_size(
                half_bger, half_bung, scf_tensors, screening, molecule, basis,
                self.dft, self.xcfun, self.molgrid, self.gs_density)

            if self.rank == mpi_master():
                e2bger = np.vstack((half_e2bger, half_e2bger))
                e2bung = np.vstack((half_e2bung, -half_e2bung))

        if self.rank == mpi_master():
            half_s2bung, half_s2bger = e2x_drv.s2n_half_size(
                half_bger, half_bung, scf_tensors, nocc)
            s2bung = np.vstack((half_s2bung, -half_s2bung))
            s2bger = np.vstack((half_s2bger, half_s2bger))

        solutions = {}
        residuals = {}
        relative_residual_norm = {}
        kappas = {}

        if self.timing:
            self.timing_dict['fock_build'][0] += tm.time() - timing_t0
            timing_t0 = tm.time()

        for iteration in range(self.max_iter):
            self.cur_iter = iteration

            if self.timing:
                self.timing_dict['reduced_space'].append(0.0)
                self.timing_dict['ortho_norm'].append(0.0)
                self.timing_dict['fock_build'].append(0.0)

            if self.rank == mpi_master():
                nvs = []

                for op, w in op_freq_keys:
                    if iteration == 0 or (relative_residual_norm[(op, w)] >
                                          self.conv_thresh):
                        if self.nonlinear:
                            grad = v1[(op, w)]
                        else:
                            grad = v1[op]

                        gradger, gradung = self.decomp_sym(grad)

                        # projections onto gerade and ungerade subspaces:

                        g_realger = np.matmul(bger.T, gradger.real)
                        g_imagger = np.matmul(bger.T, gradger.imag)
                        g_realung = np.matmul(bung.T, gradung.real)
                        g_imagung = np.matmul(bung.T, gradung.imag)

                        e2gg = np.matmul(bger.T, e2bger)
                        e2uu = np.matmul(bung.T, e2bung)
                        s2ug = np.matmul(bung.T, s2bung)

                        n_ger = bger.shape[1]
                        n_ung = bung.shape[1]

                        # creating gradient and matrix for linear equation

                        size = 2 * (n_ger + n_ung)

                        # gradient

                        g = np.zeros(size)

                        g[:n_ger] = g_realger[:]
                        g[n_ger:n_ger + n_ung] = g_realung[:]
                        g[n_ger + n_ung:size - n_ger] = -g_imagung[:]
                        g[size - n_ger:] = -g_imagger[:]

                        # matrix

                        mat = np.zeros((size, size))

                        # filling E2gg

                        mat[:n_ger, :n_ger] = e2gg[:, :]
                        mat[size - n_ger:, size - n_ger:] = -e2gg[:, :]

                        # filling E2uu

                        mat[n_ger:n_ger + n_ung, n_ger:n_ger +
                            n_ung] = e2uu[:, :]

                        mat[n_ger + n_ung:size - n_ger, n_ger + n_ung:size -
                            n_ger] = -e2uu[:, :]

                        # filling S2ug

                        mat[n_ger:n_ger + n_ung, :n_ger] = -w * s2ug[:, :]

                        mat[n_ger + n_ung:size - n_ger, :n_ger] = d * s2ug[:, :]

                        mat[n_ger:n_ger + n_ung, size - n_ger:] = d * s2ug[:, :]

                        mat[n_ger + n_ung:size - n_ger, size -
                            n_ger:] = w * s2ug[:, :]

                        # filling S2ug.T (interchanging of row and col)

                        mat[:n_ger, n_ger:n_ger + n_ung] = -w * s2ug.T[:, :]

                        mat[:n_ger, n_ger + n_ung:size -
                            n_ger] = d * s2ug.T[:, :]

                        mat[size - n_ger:, n_ger:n_ger +
                            n_ung] = d * s2ug.T[:, :]

                        mat[size - n_ger:, n_ger + n_ung:size -
                            n_ger] = w * s2ug.T[:, :]

                        # solving matrix equation

                        c = np.linalg.solve(mat, g)

                        # extracting the 4 components of c...

                        c_realger = c[:n_ger]
                        c_realung = c[n_ger:n_ger + n_ung]
                        c_imagung = c[n_ger + n_ung:size - n_ger]
                        c_imagger = c[size - n_ger:]

                        # ...and projecting them onto respective subspace

                        x_realger = np.matmul(bger, c_realger)
                        x_realung = np.matmul(bung, c_realung)
                        x_imagung = np.matmul(bung, c_imagung)
                        x_imagger = np.matmul(bger, c_imagger)

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

                        e2realger = np.matmul(e2bger, c_realger)
                        e2imagger = np.matmul(e2bger, c_imagger)
                        s2realger = np.matmul(s2bung, c_realger)
                        s2imagger = np.matmul(s2bung, c_imagger)

                        e2realung = np.matmul(e2bung, c_realung)
                        e2imagung = np.matmul(e2bung, c_imagung)
                        s2realung = np.matmul(s2bger, c_realung)
                        s2imagung = np.matmul(s2bger, c_imagung)

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

                if not self.nonlinear:
                    self.ostream.print_info(
                        '{:d} gerade trial vectors'.format(n_ger))
                    self.ostream.print_info(
                        '{:d} ungerade trial vectors'.format(n_ung))
                    self.ostream.print_blank()

                    self.print_iteration(relative_residual_norm, nvs)

            if self.timing:
                tid = iteration + 1
                self.timing_dict['reduced_space'][tid] += tm.time() - timing_t0
                timing_t0 = tm.time()

            # check convergence

            self.check_convergence(relative_residual_norm)

            if self.is_converged:
                break

            # spawning new trial vectors from residuals

            if self.rank == mpi_master():

                new_trials_ger, new_trials_ung = self.setup_trials(
                    residuals,
                    pre=precond,
                    bger=bger,
                    bung=bung,
                    res_norm=relative_residual_norm)

                assert_msg_critical(
                    new_trials_ger.any() or new_trials_ung.any(),
                    'ComplexResponseSolver: unable to add new trial vectors')

                if not new_trials_ger.any():
                    new_trials_ger = np.zeros((new_trials_ung.shape[0], 0))
                if not new_trials_ung.any():
                    new_trials_ung = np.zeros((new_trials_ger.shape[0], 0))

                # creating new sigma and rho linear transformations

                bger = np.append(bger, new_trials_ger, axis=1)
                bung = np.append(bung, new_trials_ung, axis=1)

            if self.timing:
                tid = iteration + 1
                self.timing_dict['ortho_norm'][tid] += tm.time() - timing_t0
                timing_t0 = tm.time()

            half_new_ger = None
            half_new_ung = None
            if self.rank == mpi_master():
                if new_trials_ger is not None:
                    half_new_ger = new_trials_ger[:new_trials_ger.shape[0] // 2]
                if new_trials_ung is not None:
                    half_new_ung = new_trials_ung[:new_trials_ung.shape[0] // 2]

            half_new_e2bger, half_new_e2bung = e2x_drv.e2n_half_size(
                half_new_ger, half_new_ung, scf_tensors, screening, molecule,
                basis, self.dft, self.xcfun, self.molgrid, self.gs_density)

            if self.rank == mpi_master():
                new_e2bger = np.vstack((half_new_e2bger, half_new_e2bger))
                new_e2bung = np.vstack((half_new_e2bung, -half_new_e2bung))

                half_new_s2bung, half_new_s2bger = e2x_drv.s2n_half_size(
                    half_new_ger, half_new_ung, scf_tensors, nocc)
                new_s2bung = np.vstack((half_new_s2bung, -half_new_s2bung))
                new_s2bger = np.vstack((half_new_s2bger, half_new_s2bger))

                e2bger = np.append(e2bger, new_e2bger, axis=1)
                e2bung = np.append(e2bung, new_e2bung, axis=1)

                s2bung = np.append(s2bung, new_s2bung, axis=1)
                s2bger = np.append(s2bger, new_s2bger, axis=1)

                write_rsp_hdf5(
                    self.checkpoint_file, [bger, bung, e2bger, e2bung],
                    ['CLR_bger', 'CLR_bung', 'CLR_e2bger', 'CLR_e2bung'],
                    molecule.nuclear_repulsion_energy(),
                    molecule.elem_ids_to_numpy(), basis.get_label(),
                    dft_func_label, self.ostream)

            if self.timing:
                tid = iteration + 1
                self.timing_dict['fock_build'][tid] += tm.time() - timing_t0
                timing_t0 = tm.time()

        # converged?
        if self.rank == mpi_master():
            if not self.nonlinear:
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

        if not self.nonlinear:
            a_rhs = get_rhs(self.a_operator, self.a_components, molecule, basis,
                            scf_tensors, self.rank, self.comm)

        if self.rank == mpi_master():
            if not self.nonlinear:
                va = {op: v for op, v in zip(self.a_components, a_rhs)}
                props = {}
                for aop in self.a_components:
                    for bop, w in solutions:
                        props[(aop, bop, w)] = -np.dot(va[aop],
                                                       solutions[(bop, w)])
                self.print_properties(props)
                return {
                    'properties': props,
                    'solutions': solutions,
                    'kappas': kappas,
                }
            else:
                return {
                    'solutions': solutions,
                    'kappas': kappas,
                }
        else:
            return {}

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

    def print_iteration(self, relative_residual_norm, nvs):
        """
        Prints information of the iteration.

        :param relative_residual_norm:
            Relative residual norms.
        :param nvs:
            A list of tuples containing operator component, frequency, and
            property.
        """

        width = 92

        output_header = '*** Iteration:   {} '.format(self.cur_iter + 1)
        output_header += '* Residuals (Max,Min): '
        output_header += '{:.2e} and {:.2e}'.format(
            max(relative_residual_norm.values()),
            min(relative_residual_norm.values()))
        self.ostream.print_header(output_header.ljust(width))
        self.ostream.print_blank()

        output_header = 'Operator:  {} ({})'.format(self.b_operator,
                                                    self.b_components)
        self.ostream.print_header(output_header.ljust(width))
        self.ostream.print_blank()

        for op, freq, nv in nvs:
            ops_label = '<<{};{}>>_{:.4f}'.format(op, op, freq)
            rel_res = relative_residual_norm[(op, freq)]
            output_iter = '{:<15s}: {:15.8f} {:15.8f}j   '.format(
                ops_label, -nv.real, -nv.imag)
            output_iter += 'Residual Norm: {:.8f}'.format(rel_res)
            self.ostream.print_header(output_iter.ljust(width))
        self.ostream.print_blank()
        self.ostream.flush()

    def print_header(self):
        """
        Prints complex linear response solver setup header to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header("Complex Response Driver Setup")
        self.ostream.print_header(31 * "=")
        self.ostream.print_blank()

        width = 60

        cur_str = "Damping Parameter (gamma) : {:.6e}".format(self.damping)
        self.ostream.print_header(cur_str.ljust(width))

        cur_str = "Max. Number of Iterations : " + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = "Convergence Threshold     : " + \
            "{:.1e}".format(self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(width))

        cur_str = "ERI Screening Scheme      : " + get_qq_type(self.qq_type)
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = "ERI Screening Threshold   : " + \
            "{:.1e}".format(self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(width))
        self.ostream.print_blank()

        if self.dft:
            cur_str = "Exchange-Correlation Functional : "
            cur_str += self.xcfun.get_func_label().upper()
            self.ostream.print_header(cur_str.ljust(width))
            cur_str = "Molecular Grid Level            : " + str(
                self.grid_level)
            self.ostream.print_header(cur_str.ljust(width))
        self.ostream.print_blank()

        self.ostream.flush()

    def print_convergence(self):
        """
        Prints information after convergence.
        """

        width = 92
        output_conv = '*** '
        if self.is_converged:
            output_conv += 'Complex response converged'
        else:
            output_conv += 'Complex response NOT converged'
        output_conv += ' in {:d} iterations. '.format(self.cur_iter + 1)
        output_conv += 'Time: {:.2f} sec'.format(tm.time() - self.start_time)
        self.ostream.print_header(output_conv.ljust(width))
        self.ostream.print_blank()

    def print_properties(self, props):
        """
        Prints properties.

        :param props:
            The dictionary of properties.
        """

        width = 92
        for w in self.frequencies:
            w_str = '{}, {}, w={:.4f}'.format(self.a_operator, self.b_operator,
                                              w)
            self.ostream.print_header(w_str.ljust(width))
            self.ostream.print_header(('-' * len(w_str)).ljust(width))
            for a in self.a_components:
                for b in self.b_components:
                    ops_label = '<<{};{}>>_{:.4f}'.format(a, b, w)
                    output_alpha = '{:<15s} {:15.8f} {:15.8f}j'.format(
                        ops_label, props[(a, b, w)].real, props[(a, b, w)].imag)
                    self.ostream.print_header(output_alpha.ljust(width))
            self.ostream.print_blank()

    def print_timing(self):
        """
        Prints timing for the complex response solver.
        """

        width = 92

        valstr = 'Timing (in sec):'
        self.ostream.print_header(valstr.ljust(width))
        self.ostream.print_header(('-' * len(valstr)).ljust(width))

        valstr = '{:<15s} {:>15s} {:>15s} {:>15s}'.format(
            '', 'ReducedSpace', 'Orthonorm.', 'FockBuild')
        self.ostream.print_header(valstr.ljust(width))

        for i, (a, b, c) in enumerate(
                zip(self.timing_dict['reduced_space'],
                    self.timing_dict['ortho_norm'],
                    self.timing_dict['fock_build'])):
            if i == 0:
                title = 'Initial guess'
            else:
                title = 'Iteration {:<5d}'.format(i)
            valstr = '{:<15s} {:15.3f} {:15.3f} {:15.3f}'.format(title, a, b, c)
            self.ostream.print_header(valstr.ljust(width))

        valstr = '---------'
        self.ostream.print_header(valstr.ljust(width))

        valstr = '{:<15s} {:15.3f} {:15.3f} {:15.3f}'.format(
            'Sum', sum(self.timing_dict['reduced_space']),
            sum(self.timing_dict['ortho_norm']),
            sum(self.timing_dict['fock_build']))
        self.ostream.print_header(valstr.ljust(width))

        self.ostream.print_blank()
