import numpy as np
import time as tm
import os

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import mpi_master
from .veloxchemlib import GridDriver
from .veloxchemlib import MolecularGrid
from .veloxchemlib import XCFunctional
from .veloxchemlib import denmat
from .veloxchemlib import parse_xc_func
from .aodensitymatrix import AODensityMatrix
from .lrmatvecdriver import LinearResponseMatrixVectorDriver
from .lrmatvecdriver import remove_linear_dependence_half
from .lrmatvecdriver import orthogonalize_gram_schmidt_half
from .lrmatvecdriver import normalize_half
from .lrmatvecdriver import construct_ed_sd_half
from .lrmatvecdriver import lrvec2mat
from .lrmatvecdriver import get_complex_rhs
from .lrmatvecdriver import read_rsp_hdf5
from .lrmatvecdriver import write_rsp_hdf5
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_scheme
from .qqscheme import get_qq_type
from .inputparser import parse_frequencies


class C6Solver:
    """
    Implements the complex linear response solver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - a_operator: The A operator.
        - a_components: Cartesian components of the A operator.
        - b_operator: The B operator.
        - b_components: Cartesian components of the B operator.
        - imagfrequencies: The imaginary frequencies.
        - qq_type: The electron repulsion integrals screening scheme.
        - eri_thresh: The electron repulsion integrals screening threshold.
        - dft: The flag for running DFT.
        - grid_level: The accuracy level of DFT grid.
        - xcfun: The XC functional.
        - pe: The flag for running polarizable embedding calculation.
        - potfile: The name of the potential file for polarizable embedding.
        - use_split_comm: The flag for using split communicators.
        - max_iter: The maximum number of solver iterations.
        - conv_thresh: The convergence threshold for the solver.
        - cur_iter: Index of the current iteration.
        - is_converged: The flag for convergence.
        - small_thresh: The norm threshold for a vector to be considered a zero
          vector.
        - lindep_thresh: The threshold for removing linear dependence in the
          trial vectors.
        - comm: The MPI communicator.
        - rank: The MPI rank.
        - nodes: Number of MPI processes.
        - ostream: The output stream.
        - restart: The flag for restarting from checkpoint file.
        - checkpoint_file: The name of checkpoint file.
        - checkpoint_time: The timer of checkpoint file.
        - timing: The flag for printing timing information.
        - profiling: The flag for printing profiling information.
    """

    def __init__(self, comm, ostream):
        """
        Initializes complex linear response solver to default setup.
        """

        self.a_operator = 'dipole'
        self.a_components = 'xyz'
        self.b_operator = 'dipole'
        self.b_components = 'xyz'

        self.n_points = 12
        self.imagfrequencies = [0.3*(1-t)/(1+t) for t in np.polynomial.legendre.leggauss(self.n_points)][0]

        self.qq_type = 'QQ_DEN'
        self.eri_thresh = 1.0e-15

        self.dft = False
        self.grid_level = 4
        self.xcfun = XCFunctional()

        self.pe = False
        self.potfile = None

        self.use_split_comm = False

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
        self.checkpoint_time = None

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

        if 'n_points' in rsp_dict:
            self.n_points = int(rsp_dict['n_points'])
            self.imagfrequencies = [0.3*(1-t)/(1+t) for t in np.polynomial.legendre.leggauss(self.n_points)][0]
        #if 'imagfrequencies' in rsp_dict:
        #    self.imagfrequencies = parse_frequencies(rsp_dict['imagfrequencies'])

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

        if 'pe' in method_dict:
            key = method_dict['pe'].lower()
            self.pe = True if key == 'yes' else False
        if 'potfile' in method_dict:
            if 'pe' not in method_dict:
                self.pe = True
            self.potfile = method_dict['potfile']

        if 'use_split_comm' in method_dict:
            key = method_dict['use_split_comm'].lower()
            self.use_split_comm = True if key == 'yes' else False

    def decomp_trials(self, vecs):
        """
        Decomposes trial vectors into their 4 respective parts (real gerade,
        real ungerade, imaginary gerade, and imaginary ungerad).

        :param vecs:
            The trial vectors.

        :return:
            A tuple containing respective parts of the trial vectors.
        """

        half_rows = vecs.shape[0] // 2

        if len(vecs.shape) != 1:
            realung = []
            imagger = []
            for vec in range(len(vecs[0, :])):
                realung.append(vecs[:half_rows, vec])
                imagger.append(vecs[half_rows:, vec])
        else:
            realung = vecs[:half_rows]
            imagger = vecs[half_rows:]

        return np.array(realung).T, np.array(imagger).T

    def decomp_grad(self, grad):
        """
        Decomposes gradient into gerade and ungerade parts.

        :param grad:
            The trial vectors.

        :return:
            A tuple containing gerade and ungerade parts of vectors.
        """

        assert_msg_critical(
            len(grad.shape) == 1, 'decomp_grad: Expecting a 1D array')

        assert_msg_critical(grad.shape[0] % 2 == 0,
                            'decomp_grad: size of array should be even')

        half_size = grad.shape[0] // 2

        grad_T = np.zeros_like(grad)
        grad_T[:half_size] = grad[half_size:]
        grad_T[half_size:] = grad[:half_size]

        ger = 0.5 * (grad + grad_T)[:half_size]
        ung = 0.5 * (grad - grad_T)[:half_size]

        return ger.T, ung.T

    def get_precond(self, orb_ene, nocc, norb, iw):
        """
        Constructs the preconditioner matrix.

        :param orb_ene:
            The orbital energies.
        :param nocc:
            The number of doubly occupied orbitals.
        :param norb:
            The number of orbitals.
        :param iw:
            The imaginary frequency.

        :return:
            The preconditioner matrix.
        """

        # spawning needed components

        ediag, sdiag = construct_ed_sd_half(orb_ene, nocc, norb)
        ediag_sq = ediag**2
        sdiag_sq = sdiag**2
        iw_sq = iw**2

        # constructing matrix block diagonals

        a_diag = ediag * (ediag_sq + iw_sq * sdiag_sq)
        c_diag = (iw * sdiag) * (ediag_sq + iw_sq * sdiag_sq)
        p_diag = 1.0 / (ediag_sq + iw_sq * sdiag_sq)**2

        pa_diag = p_diag * a_diag
        pc_diag = p_diag * c_diag

        precond = np.array([pa_diag, pc_diag])

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

        pa, pc = precond[0], precond[1]

        v_in_ru, v_in_ig = self.decomp_trials(v_in)

        v_out_ru = pa * v_in_ru + pc * v_in_ig
        v_out_ig = pc * v_in_ru - pa * v_in_ig

        v_out = np.array([v_out_ru, v_out_ig]).flatten()

        return v_out

    def initial_guess(self, v1, precond):
        """
        Creating initial guess (un-orthonormalized trials) out of gradients.

        :param v1:
            The dictionary containing (operator, imaginary frequency) as keys
            and right-hand sides as values.
        :param precond:
            The preconditioner matrices.

        :return:
            The initial guess.
        """

        ig = {}
        for (op, iw), grad in v1.items():
            gradger, gradung = self.decomp_grad(grad)

            grad = np.array(
                [gradung.real, -gradger.imag]).flatten()
            gn = np.sqrt(2.0) * np.linalg.norm(grad)

            if gn < self.small_thresh:
                ig[(op, iw)] = np.zeros(grad.shape[0])
            else:
                ig[(op, iw)] = self.preconditioning(precond[iw], grad)

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
                if np.linalg.norm(v) * np.sqrt(2.0) > self.small_thresh:
                    trials.append(v)

        new_trials = np.array(trials).T

        # decomposing the full space trial vectors...

        new_realung, new_imagger = self.decomp_trials(new_trials)

        # ...and assembling gerade and ungerade subspaces

        new_ger = new_imagger
        new_ung = new_realung

        # orthogonalizing new trial vectors against existing ones

        if bger is not None and bger.any():
            new_ger = new_ger - 2.0 * np.matmul(bger, np.matmul(
                bger.T, new_ger))

        if bung is not None and bung.any():
            new_ung = new_ung - 2.0 * np.matmul(bung, np.matmul(
                bung.T, new_ung))

        # orthonormalizing new trial vectors

        if new_ger.any() and renormalize:

            new_ger = remove_linear_dependence_half(new_ger, self.lindep_thresh)
            new_ger = orthogonalize_gram_schmidt_half(new_ger)
            new_ger = normalize_half(new_ger)

        if new_ung.any() and renormalize:

            new_ung = remove_linear_dependence_half(new_ung, self.lindep_thresh)
            new_ung = orthogonalize_gram_schmidt_half(new_ung)
            new_ung = normalize_half(new_ung)

        return new_ger, new_ung

    def compute(self, molecule, basis, scf_tensors, v1=None):
        """
        Solves for the response vector iteratively while checking the residuals
        for convergence.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param v1:
            The gradients on the right-hand side. If not provided, v1 will be
            computed for the B operator.

        :return:
            A dictionary containing response functions, solutions, and kappas.
        """

        if self.profiling:
            import cProfile
            import pstats
            import io
            pr = cProfile.Profile()
            pr.enable()

        if self.timing:
            self.timing_dict = {
                'reduced_space': [0.0],
                'ortho_norm': [0.0],
                'fock_build': [0.0],
            }
            timing_t0 = tm.time()

        if self.rank == mpi_master():
            self.print_header()

        self.start_time = tm.time()
        self.checkpoint_time = self.start_time

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'ComplexResponseSolver: not implemented for unrestricted case')

        # generate integration grid
        if self.dft:
            grid_drv = GridDriver(self.comm)
            grid_drv.set_level(self.grid_level)

            grid_t0 = tm.time()
            molgrid = grid_drv.generate(molecule)
            n_grid_points = molgrid.number_of_points()
            self.ostream.print_info(
                'Molecular grid with {0:d} points generated in {1:.2f} sec.'.
                format(n_grid_points,
                       tm.time() - grid_t0))
            self.ostream.print_blank()

            if self.rank == mpi_master():
                gs_density = AODensityMatrix([scf_tensors['D'][0]], denmat.rest)
            else:
                gs_density = AODensityMatrix()
            gs_density.broadcast(self.rank, self.comm)

            dft_func_label = self.xcfun.get_func_label().upper()
        else:
            molgrid = MolecularGrid()
            gs_density = AODensityMatrix()
            dft_func_label = 'HF'

        # set up polarizable embedding
        if self.pe:
            from .polembed import PolEmbed
            pe_drv = PolEmbed(molecule, basis, self.comm, self.potfile)
            V_es = pe_drv.compute_multipole_potential_integrals()

            pot_info = "Reading polarizable embedding potential: {}".format(
                self.potfile)
            self.ostream.print_info(pot_info)
            self.ostream.print_blank()

            with open(self.potfile, 'r') as f_pot:
                potfile_text = os.linesep.join(f_pot.readlines())
        else:
            pe_drv = None
            V_es = None
            potfile_text = ''

        # generate screening for ERI

        if self.use_split_comm:
            self.use_split_comm = ((self.dft or self.pe) and self.nodes >= 8)

        if self.use_split_comm:
            screening = None
            valstr = 'ERI'
            if self.dft:
                valstr += '/DFT'
            if self.pe:
                valstr += '/PE'
            self.ostream.print_info(
                'Using sub-communicators for {}.'.format(valstr))
            self.ostream.print_blank()
        else:
            eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
            screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                        self.eri_thresh, molecule, basis)

        e2x_drv = LinearResponseMatrixVectorDriver(self.comm,
                                                   self.use_split_comm)
        e2x_drv.update_settings(self.eri_thresh, self.qq_type, self.dft,
                                self.xcfun, self.pe, self.potfile)
        timing_dict = {}

        if self.rank == mpi_master():
            nocc = molecule.number_of_alpha_electrons()
            norb = scf_tensors['C'].shape[1]
            orb_ene = scf_tensors['E']
        else:
            nocc = None

        nonlinear_flag = False

        if not v1:
            b_rhs = get_complex_rhs(self.b_operator, self.b_components,
                                    molecule, basis, scf_tensors, self.rank,
                                    self.comm)
            if self.rank == mpi_master():
                v1 = {(op, iw): v for op, v in zip(self.b_components, b_rhs)
                      for iw in self.imagfrequencies}
        else:
            nonlinear_flag = True

        if self.rank == mpi_master():
            #d = self.damping
            imagfreqs = [op_imagfreq[1] for op_imagfreq in v1]
            op_imagfreq_keys = list(v1.keys())

            # creating the preconditioner matrix

            precond = {
                iw: self.get_precond(orb_ene, nocc, norb, iw) for iw in imagfreqs
            }

        bger = None
        bung = None
        new_trials_ger = None
        new_trials_ung = None

        # read initial guess from restart file

        if self.restart:
            if self.rank == mpi_master():
                bger, bung, e2bger, e2bung = read_rsp_hdf5(
                    self.checkpoint_file,
                    ['CLR_bger', 'CLR_bung', 'CLR_e2bger', 'CLR_e2bung'],
                    molecule.nuclear_repulsion_energy(),
                    molecule.elem_ids_to_numpy(), basis.get_label(),
                    dft_func_label, potfile_text, self.ostream)
                self.restart = (bger is not None and bung is not None and
                                e2bger is not None and e2bung is not None)
            self.restart = self.comm.bcast(self.restart, root=mpi_master())

        # generate initial guess from scratch

        if not self.restart:
            if self.rank == mpi_master():

                # spawning initial trial vectors

                igs = self.initial_guess(v1, precond)
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

            e2bger, e2bung = e2x_drv.e2n_half_size(bger, bung, scf_tensors,
                                                   screening, molecule, basis,
                                                   molgrid, gs_density, V_es,
                                                   pe_drv, timing_dict)

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

                n_ger = bger.shape[1]
                n_ung = bung.shape[1]

                e2gg = 2.0 * np.matmul(bger.T, e2bger)
                e2uu = 2.0 * np.matmul(bung.T, e2bung)
                s2ug = 4.0 * np.matmul(bung.T, bger)

                for op, iw in op_imagfreq_keys:
                    if iteration == 0 or (relative_residual_norm[(op, iw)] >
                                          self.conv_thresh):
                        grad = v1[(op, iw)]
                        gradger, gradung = self.decomp_grad(grad)

                        # projections onto gerade and ungerade subspaces:

                        g_realung = 2.0 * np.matmul(bung.T, gradung.real)

                        # creating gradient and matrix for linear equation

                        size = n_ger + n_ung

                        # gradient

                        g = np.zeros(size)

                        g[:n_ung] = g_realung[:]

                        # matrix

                        mat = np.zeros((size, size))

                        # filling E2gg

                        mat[n_ung:n_ung + n_ger, n_ung:n_ung +
                            n_ger] = -e2gg[:, :]

                        # filling E2uu

                        mat[:n_ung, :n_ung] = e2uu[:, :]

                        # filling S2ug

                        mat[:n_ung, n_ung:n_ung + n_ger] = iw * s2ug[:, :]

                        # filling S2ug.T (interchanging of row and col)

                        mat[n_ung:n_ung + n_ger, :n_ung] = iw * s2ug.T[:, :]

                        # solving matrix equation

                        c = np.linalg.solve(mat, g)

                        # extracting the 2 components of c...

                        c_realung = c[:n_ung]
                        c_imagger = c[n_ung:n_ung + n_ger]

                        # ...and projecting them onto respective subspace

                        x_realung = np.matmul(bung, c_realung)
                        x_imagger = np.matmul(bger, c_imagger)

                        # composing full size response vector

                        x_realung_full = np.hstack((x_realung, -x_realung))
                        x_imagger_full = np.hstack((x_imagger, x_imagger))

                        x_real = x_realung_full
                        x_imag = x_imagger_full
                        x = x_real + 1j * x_imag

                        solutions[(op, iw)] = x

                        kappas[(op, iw)] = (lrvec2mat(x.real, nocc, norb) +
                                           1j * lrvec2mat(x.imag, nocc, norb))

                        # composing E2 and S2 matrices projected onto solution
                        # subspace

                        e2imagger = np.matmul(e2bger, c_imagger)
                        s2imagger = 2.0 * np.matmul(bger, c_imagger)

                        e2realung = np.matmul(e2bung, c_realung)
                        s2realung = 2.0 * np.matmul(bung, c_realung)

                        # calculating the residual components

                        r_realung = (e2realung + iw * s2imagger -
                                     gradung.real)
                        r_imagger = (-e2imagger +
                                     iw * s2realung + gradger.imag)

                        # composing total half-sized residual

                        residuals[(op, iw)] = np.array(
                            [r_realung, r_imagger]).flatten()

                        r = residuals[(op, iw)]
                        n = solutions[(op, iw)]

                        # calculating relative residual norm
                        # for convergence check

                        nv = np.matmul(n, grad)
                        nvs.append((op, iw, nv))

                        rn = np.sqrt(2.0) * np.linalg.norm(r)
                        nn = np.linalg.norm(n)
                        if nn != 0:
                            relative_residual_norm[(op, iw)] = rn / nn
                        else:
                            relative_residual_norm[(op, iw)] = 0

                # write to output

                self.ostream.print_info(
                    '{:d} gerade trial vectors in reduced space'.format(n_ger))
                self.ostream.print_info(
                    '{:d} ungerade trial vectors in reduced space'.format(
                        n_ung))
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

            new_e2bger, new_e2bung = e2x_drv.e2n_half_size(
                new_trials_ger, new_trials_ung, scf_tensors, screening,
                molecule, basis, molgrid, gs_density, V_es, pe_drv, timing_dict)

            if self.rank == mpi_master():

                e2bger = np.append(e2bger, new_e2bger, axis=1)
                e2bung = np.append(e2bung, new_e2bung, axis=1)

                if tm.time() - self.checkpoint_time > 900.0:
                    write_rsp_hdf5(
                        self.checkpoint_file, [bger, bung, e2bger, e2bung],
                        ['CLR_bger', 'CLR_bung', 'CLR_e2bger', 'CLR_e2bung'],
                        molecule.nuclear_repulsion_energy(),
                        molecule.elem_ids_to_numpy(), basis.get_label(),
                        dft_func_label, potfile_text, self.ostream)
                    self.checkpoint_time = tm.time()

            if self.timing:
                tid = iteration + 1
                self.timing_dict['fock_build'][tid] += tm.time() - timing_t0
                timing_t0 = tm.time()

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

        if not nonlinear_flag:
            a_rhs = get_complex_rhs(self.a_operator, self.a_components,
                                    molecule, basis, scf_tensors, self.rank,
                                    self.comm)

            if self.rank == mpi_master():
                va = {op: v for op, v in zip(self.a_components, a_rhs)}
                rsp_funcs = {}
                for aop in self.a_components:
                    for bop, iw in solutions:
                        rsp_funcs[(aop, bop, iw)] = -np.dot(va[aop],
                                                    solutions[(bop, iw)])
                return {
                    'response_functions': rsp_funcs,
                    'solutions': solutions,
                    'kappas': kappas
                }

        else:
            if self.rank == mpi_master():
                return {'solutions': solutions, 'kappas': kappas}

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
            A list of tuples containing operator component, imaginary frequency,
            and property.
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

        for op, imagfreq, nv in nvs:
            ops_label = '<<{};{}>>_i{:.4f}'.format(op, op, imagfreq)
            rel_res = relative_residual_norm[(op, imagfreq)]
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
        self.ostream.print_header("C6 Value Response Driver Setup")
        self.ostream.print_header(31 * "=")
        self.ostream.print_blank()

        width = 60

        cur_str = "Number of integration points    : " + str(self.n_points)
        self.ostream.print_header(cur_str.ljust(width))

        cur_str = "Max. Number of Iterations       : " + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = "Convergence Threshold           : " + \
            "{:.1e}".format(self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(width))

        cur_str = "ERI Screening Scheme            : " + get_qq_type(
            self.qq_type)
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = "ERI Screening Threshold         : " + \
            "{:.1e}".format(self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(width))

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
