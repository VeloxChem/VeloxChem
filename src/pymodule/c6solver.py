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
from .profiler import Profiler
from .distributedarray import DistributedArray
from .aodensitymatrix import AODensityMatrix
from .lrmatvecdriver import LinearResponseMatrixVectorDriver
from .lrmatvecdriver import remove_linear_dependence_half
from .lrmatvecdriver import orthogonalize_gram_schmidt_half
from .lrmatvecdriver import normalize_half
from .lrmatvecdriver import construct_ed_sd_half
from .lrmatvecdriver import get_complex_rhs
from .lrmatvecdriver import check_rsp_hdf5
from .lrmatvecdriver import write_rsp_hdf5
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_scheme
from .qqscheme import get_qq_type


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
        - n_points: The number of integration points.
        - w0: The transformation function prefactor.
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
        - memory_profiling: The flag for printing memory usage.
        - memory_tracing: The flag for tracing memory allocation using
          tracemalloc.
    """

    def __init__(self, comm, ostream):
        """
        Initializes C6 solver to default setup.
        """

        self.a_operator = 'dipole'
        self.a_components = 'xyz'
        self.b_operator = 'dipole'
        self.b_components = 'xyz'

        self.n_points = 9
        self.w0 = 0.3

        self.qq_type = 'QQ_DEN'
        self.eri_thresh = 1.0e-15

        self.dft = False
        self.grid_level = 4
        self.xcfun = XCFunctional()

        self.pe = False
        self.potfile = None

        self.use_split_comm = False

        self.max_iter = 150
        self.conv_thresh = 1.0e-3

        self.cur_iter = 0
        self.is_converged = False
        self.small_thresh = 1.0e-10
        self.lindep_thresh = 1.0e-10

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream

        self.restart = True
        self.checkpoint_file = None
        self.checkpoint_time = None
        self.checkpoint_interval = 4.0 * 3600

        self.timing = False
        self.profiling = False
        self.memory_profiling = False
        self.memory_tracing = False

    def update_settings(self, rsp_dict, method_dict={}):
        """
        Updates response and method settings in C6 solver.

        :param rsp_dict:
            The dictionary of response dict.
        :param method_dict:
            The dictionary of method rsp_dict.
        """

        if 'n_points' in rsp_dict:
            self.n_points = int(rsp_dict['n_points'])
        if 'w0' in rsp_dict:
            self.w0 = float(rsp_dict['w0'])

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
        if 'checkpoint_interval' in rsp_dict:
            if 'h' in rsp_dict['checkpoint_interval'].lower():
                try:
                    n_hours = float(
                        rsp_dict['checkpoint_interval'].lower().split('h')[0])
                    self.checkpoint_interval = max(0, n_hours) * 3600
                except ValueError:
                    pass

        if 'timing' in rsp_dict:
            key = rsp_dict['timing'].lower()
            self.timing = True if key in ['yes', 'y'] else False
        if 'profiling' in rsp_dict:
            key = rsp_dict['profiling'].lower()
            self.profiling = True if key in ['yes', 'y'] else False
        if 'memory_profiling' in rsp_dict:
            key = rsp_dict['memory_profiling'].lower()
            self.memory_profiling = True if key in ['yes', 'y'] else False
        if 'memory_tracing' in rsp_dict:
            key = rsp_dict['memory_tracing'].lower()
            self.memory_tracing = True if key in ['yes', 'y'] else False

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
        Decomposes trial vectors into their 2 respective non-zero parts (real
        ungerade and imaginary gerade).

        :param vecs:
            The trial vectors.

        :return:
            A tuple containing respective non-zero parts of the trial vectors.
        """

        realung, imagger = None, None
        half_rows = vecs.shape[0] // 2

        if vecs.ndim == 1:
            realung = vecs[:half_rows]
            imagger = vecs[half_rows:]

        elif vecs.ndim == 2:
            realung = vecs[:half_rows, :]
            imagger = vecs[half_rows:, :]

        return realung, imagger

    def decomp_grad(self, grad):
        """
        Decomposes gradient into gerade and ungerade parts.

        :param grad:
            The trial vectors.

        :return:
            A tuple containing gerade and ungerade parts of vectors.
        """

        assert_msg_critical(grad.ndim == 1, 'decomp_grad: Expecting a 1D array')

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

        return (pa_diag, pc_diag)

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

        pa, pc = precond

        v_in_ru, v_in_ig = self.decomp_trials(v_in)

        v_out_ru = pa * v_in_ru + pc * v_in_ig
        v_out_ig = pc * v_in_ru - pa * v_in_ig

        return np.hstack((v_out_ru, v_out_ig))

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

            grad = np.hstack((gradung.real, -gradger.imag))
            gn = np.sqrt(2.0) * np.linalg.norm(grad)

            if gn < self.small_thresh:
                ig[(op, iw)] = np.zeros(grad.shape[0])
            else:
                ig[(op, iw)] = self.preconditioning(precond[iw], grad)

        return ig

    def setup_trials(self,
                     vectors,
                     precond=None,
                     dist_bger=None,
                     dist_bung=None,
                     renormalize=True):
        """
        Computes orthonormalized trial vectors. Takes set of vectors,
        preconditioner matrix, gerade and ungerade subspaces as input
        arguments.

        :param vectors:
            The set of vectors.
        :param precond:
            The preconditioner matrix.
        :param dist_bger:
            The distributed gerade subspace.
        :param dist_bung:
            The distributed ungerade subspace.
        :param renormalize:
            The flag for normalization.

        :return:
            The orthonormalized gerade and ungerade trial vectors.
        """

        if self.rank == mpi_master():
            trials = []
            for (op, iw) in vectors:
                vec = np.array(vectors[(op, iw)])

                # preconditioning trials:

                if precond is not None:
                    v = self.preconditioning(precond[iw], vec)
                else:
                    v = vec
                if np.linalg.norm(v) * np.sqrt(2.0) > self.small_thresh:
                    trials.append(v)

            new_trials = np.array(trials).T

            # decomposing the full space trial vectors

            new_ung, new_ger = self.decomp_trials(new_trials)
        else:
            new_ung, new_ger = None, None

        # orthogonalizing new trial vectors against existing ones

        if dist_bger is not None:
            # t = t - (b (b.T t))
            dist_new_ger = DistributedArray(new_ger, self.comm)
            bT_new_ger = dist_bger.matmul_AtB_allreduce(dist_new_ger, 2.0)
            new_ger_proj = dist_bger.matmul_AB(bT_new_ger)
            if self.rank == mpi_master():
                new_ger -= new_ger_proj

        if dist_bung is not None:
            # t = t - (b (b.T t))
            dist_new_ung = DistributedArray(new_ung, self.comm)
            bT_new_ung = dist_bung.matmul_AtB_allreduce(dist_new_ung, 2.0)
            new_ung_proj = dist_bung.matmul_AB(bT_new_ung)
            if self.rank == mpi_master():
                new_ung -= new_ung_proj

        # orthonormalizing new trial vectors

        if self.rank == mpi_master() and renormalize:
            new_ger = remove_linear_dependence_half(new_ger, self.lindep_thresh)
            new_ger = orthogonalize_gram_schmidt_half(new_ger)
            new_ger = normalize_half(new_ger)

            new_ung = remove_linear_dependence_half(new_ung, self.lindep_thresh)
            new_ung = orthogonalize_gram_schmidt_half(new_ung)
            new_ung = normalize_half(new_ung)

        return new_ger, new_ung

    def compute(self, molecule, basis, scf_tensors):
        """
        Solves for the response vector iteratively while checking the residuals
        for convergence.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            A dictionary containing response functions and solutions.
        """

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        if self.rank == mpi_master():
            self.print_header()

        self.start_time = tm.time()
        self.checkpoint_time = self.start_time

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(nalpha == nbeta,
                            'C6Solver: not implemented for unrestricted case')

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

        b_rhs = get_complex_rhs(self.b_operator, self.b_components, molecule,
                                basis, scf_tensors, self.rank, self.comm)
        if self.rank == mpi_master():
            imagfreqs = [
                self.w0 * (1 - t) / (1 + t)
                for t in np.polynomial.legendre.leggauss(self.n_points)
            ][0]
            imagfreqs = np.append(imagfreqs, 0.0)
            v1 = {(op, iw): v for op, v in zip(self.b_components, b_rhs)
                  for iw in imagfreqs}

        if self.rank == mpi_master():
            op_imagfreq_keys = list(v1.keys())
        else:
            op_imagfreq_keys = None
        op_imagfreq_keys = self.comm.bcast(op_imagfreq_keys, root=mpi_master())

        # creating the preconditioner matrix
        if self.rank == mpi_master():
            precond = {
                iw: self.get_precond(orb_ene, nocc, norb, iw)
                for iw in imagfreqs
            }
        else:
            precond = None

        rsp_vector_labels = [
            'C6_bger_half_size',
            'C6_bung_half_size',
            'C6_e2bger_half_size',
            'C6_e2bung_half_size',
        ]

        # read initial guess from restart file

        if self.restart:
            if self.rank == mpi_master():
                self.restart = check_rsp_hdf5(
                    self.checkpoint_file, rsp_vector_labels,
                    molecule.nuclear_repulsion_energy(),
                    molecule.elem_ids_to_numpy(), basis.get_label(),
                    dft_func_label, potfile_text, self.ostream)
            self.restart = self.comm.bcast(self.restart, root=mpi_master())

        # read initial guess from restart file
        if self.restart:
            dist_bger = DistributedArray.read_from_hdf5_file(
                self.checkpoint_file, 'C6_bger_half_size', self.comm)
            dist_bung = DistributedArray.read_from_hdf5_file(
                self.checkpoint_file, 'C6_bung_half_size', self.comm)
            dist_e2bger = DistributedArray.read_from_hdf5_file(
                self.checkpoint_file, 'C6_e2bger_half_size', self.comm)
            dist_e2bung = DistributedArray.read_from_hdf5_file(
                self.checkpoint_file, 'C6_e2bung_half_size', self.comm)

        # generate initial guess from scratch
        else:
            if self.rank == mpi_master():
                igs = self.initial_guess(v1, precond)
                bger, bung = self.setup_trials(igs)

                assert_msg_critical(bger.any() or bung.any(),
                                    'C6Solver: trial vectors are empty')

                if bger is None or not bger.any():
                    bger = np.zeros((bung.shape[0], 0))
                if bung is None or not bung.any():
                    bung = np.zeros((bger.shape[0], 0))
            else:
                bger, bung = None, None

            e2bger, e2bung = e2x_drv.e2n_half_size(bger, bung, scf_tensors,
                                                   screening, molecule, basis,
                                                   molgrid, gs_density, V_es,
                                                   pe_drv, timing_dict)

            dist_bger = DistributedArray(bger, self.comm)
            dist_bung = DistributedArray(bung, self.comm)

            bger, bung = None, None

            dist_e2bger = DistributedArray(e2bger, self.comm)
            dist_e2bung = DistributedArray(e2bung, self.comm)

            e2bger, e2bung = None, None

        profiler.check_memory_usage('Initial guess')

        solutions = {}
        residuals = {}
        relative_residual_norm = {}

        # start iterations
        for iteration in range(self.max_iter):

            profiler.start_timer(iteration, 'ReducedSpace')

            xvs = []
            self.cur_iter = iteration

            n_ger = dist_bger.shape(1)
            n_ung = dist_bung.shape(1)

            e2gg = dist_bger.matmul_AtB(dist_e2bger, 2.0)
            e2uu = dist_bung.matmul_AtB(dist_e2bung, 2.0)
            s2ug = dist_bung.matmul_AtB(dist_bger, 4.0)

            for op, iw in op_imagfreq_keys:
                if (iteration == 0 or
                        relative_residual_norm[(op, iw)] > self.conv_thresh):

                    if self.rank == mpi_master():
                        grad = v1[(op, iw)]
                        gradger, gradung = self.decomp_grad(grad)
                        gradung_real = gradung.real
                    else:
                        gradung_real = None
                    dist_grad_realung = DistributedArray(
                        gradung_real, self.comm)

                    # projections onto gerade and ungerade subspaces:

                    g_realung = dist_bung.matmul_AtB(dist_grad_realung, 2.0)

                    # creating gradient and matrix for linear equation

                    size = n_ger + n_ung

                    if self.rank == mpi_master():

                        # gradient

                        g = np.zeros(size)

                        g[:n_ung] = g_realung[:]

                        # matrix

                        mat = np.zeros((size, size))

                        mat[n_ung:n_ung + n_ger,
                            n_ung:n_ung + n_ger] = -e2gg[:, :]
                        mat[:n_ung, :n_ung] = e2uu[:, :]
                        mat[:n_ung, n_ung:n_ung + n_ger] = iw * s2ug[:, :]
                        mat[n_ung:n_ung + n_ger, :n_ung] = iw * s2ug.T[:, :]

                        # solving matrix equation

                        c = np.linalg.solve(mat, g)
                    else:
                        c = None
                    c = self.comm.bcast(c, root=mpi_master())

                    # extracting the 2 components of c...

                    c_realung = c[:n_ung]
                    c_imagger = c[n_ung:n_ung + n_ger]

                    # ...and projecting them onto respective subspace

                    x_realung = dist_bung.matmul_AB(c_realung)
                    x_imagger = dist_bger.matmul_AB(c_imagger)

                    # composing full size response vector

                    if self.rank == mpi_master():

                        x_realung_full = np.hstack((x_realung, -x_realung))
                        x_imagger_full = np.hstack((x_imagger, x_imagger))

                        x_real = x_realung_full
                        x_imag = x_imagger_full
                        x = x_real + 1j * x_imag

                    # composing E2 matrices projected onto solution subspace

                    e2imagger = dist_e2bger.matmul_AB(c_imagger)
                    e2realung = dist_e2bung.matmul_AB(c_realung)

                    # calculating the residual components

                    if self.rank == mpi_master():

                        s2imagger = 2.0 * x_imagger
                        s2realung = 2.0 * x_realung

                        r_realung = (e2realung + iw * s2imagger - gradung.real)
                        r_imagger = (-e2imagger + iw * s2realung + gradger.imag)

                        # composing total half-sized residual

                        r = np.hstack((r_realung, r_imagger))

                        # calculating relative residual norm
                        # for convergence check

                        xv = np.matmul(x, grad)
                        xvs.append((op, iw, xv))

                        rn = np.sqrt(2.0) * np.linalg.norm(r)
                        xn = np.linalg.norm(x)
                        if xn != 0:
                            relative_residual_norm[(op, iw)] = rn / xn
                        else:
                            relative_residual_norm[(op, iw)] = rn

                        if relative_residual_norm[(op, iw)] < self.conv_thresh:
                            solutions[(op, iw)] = x
                        else:
                            residuals[(op, iw)] = r

            relative_residual_norm = self.comm.bcast(relative_residual_norm,
                                                     root=mpi_master())

            # write to output
            if self.rank == mpi_master():

                self.ostream.print_info(
                    '{:d} gerade trial vectors in reduced space'.format(n_ger))
                self.ostream.print_info(
                    '{:d} ungerade trial vectors in reduced space'.format(
                        n_ung))
                self.ostream.print_blank()

                mem_usage, mem_detail = profiler.get_memory_dictionary({
                    'dist_bger': dist_bger.array(),
                    'dist_bung': dist_bung.array(),
                    'dist_e2bger': dist_e2bger.array(),
                    'dist_e2bung': dist_e2bung.array(),
                    'precond': precond,
                    'solutions': solutions,
                    'residuals': residuals,
                })
                mem_avail = profiler.get_available_memory()

                self.ostream.print_info(
                    '{:s} of memory used for subspace procedure'.format(
                        mem_usage))
                if self.memory_profiling:
                    for m in mem_detail:
                        self.ostream.print_info('  {:<15s} {:s}'.format(*m))
                self.ostream.print_info(
                    '{:s} of memory available for the solver'.format(mem_avail))
                self.ostream.print_blank()

                profiler.check_memory_usage(
                    'Iteration {:d} subspace'.format(iteration + 1))

                profiler.print_memory_tracing(self.ostream)

                self.print_iteration(relative_residual_norm, xvs)

            profiler.stop_timer(iteration, 'ReducedSpace')

            # check convergence

            self.check_convergence(relative_residual_norm)

            if self.is_converged:
                break

            profiler.start_timer(iteration, 'Orthonorm.')

            # spawning new trial vectors from residuals

            new_trials_ger, new_trials_ung = self.setup_trials(
                residuals, precond, dist_bger, dist_bung)

            residuals.clear()

            if self.rank == mpi_master():

                assert_msg_critical(
                    new_trials_ger.any() or new_trials_ung.any(),
                    'C6Solver: unable to add new trial vectors')

                if new_trials_ger is None or not new_trials_ger.any():
                    new_trials_ger = np.zeros((new_trials_ung.shape[0], 0))
                if new_trials_ung is None or not new_trials_ung.any():
                    new_trials_ung = np.zeros((new_trials_ger.shape[0], 0))

            profiler.stop_timer(iteration, 'Orthonorm.')
            profiler.start_timer(iteration, 'FockBuild')

            # creating new sigma and rho linear transformations

            new_e2bger, new_e2bung = e2x_drv.e2n_half_size(
                new_trials_ger, new_trials_ung, scf_tensors, screening,
                molecule, basis, molgrid, gs_density, V_es, pe_drv, timing_dict)

            dist_bger.append(DistributedArray(new_trials_ger, self.comm),
                             axis=1)
            dist_bung.append(DistributedArray(new_trials_ung, self.comm),
                             axis=1)

            new_trials_ger, new_trials_ung = None, None

            dist_e2bger.append(DistributedArray(new_e2bger, self.comm), axis=1)
            dist_e2bung.append(DistributedArray(new_e2bung, self.comm), axis=1)

            new_e2bger, new_e2bung = None, None

            # write to checkpoint file
            if tm.time() - self.checkpoint_time >= self.checkpoint_interval:
                if self.rank == mpi_master():
                    write_rsp_hdf5(self.checkpoint_file, [], [],
                                   molecule.nuclear_repulsion_energy(),
                                   molecule.elem_ids_to_numpy(),
                                   basis.get_label(), dft_func_label,
                                   potfile_text, self.ostream)
                    self.checkpoint_time = tm.time()

                dt = 0.0
                dt += dist_bger.append_to_hdf5_file(self.checkpoint_file,
                                                    'C6_bger_half_size')
                dt += dist_bung.append_to_hdf5_file(self.checkpoint_file,
                                                    'C6_bung_half_size')
                dt += dist_e2bger.append_to_hdf5_file(self.checkpoint_file,
                                                      'C6_e2bger_half_size')
                dt += dist_e2bung.append_to_hdf5_file(self.checkpoint_file,
                                                      'C6_e2bung_half_size')
                dt_text = 'Time spent in writing to checkpoint file: '
                dt_text += '{:.2f} sec'.format(dt)
                self.ostream.print_info(dt_text)
                self.ostream.print_blank()

                self.checkpoint_time = tm.time()

            profiler.stop_timer(iteration, 'FockBuild')
            profiler.update_timer(iteration, timing_dict)

            profiler.check_memory_usage(
                'Iteration {:d} sigma build'.format(iteration + 1))

        # converged?
        if self.rank == mpi_master():
            self.print_convergence()

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of CPP solver')
        profiler.print_memory_usage(self.ostream)

        a_rhs = get_complex_rhs(self.a_operator, self.a_components, molecule,
                                basis, scf_tensors, self.rank, self.comm)

        if self.rank == mpi_master() and self.is_converged:
            va = {op: v for op, v in zip(self.a_components, a_rhs)}
            rsp_funcs = {}
            for aop in self.a_components:
                for bop, iw in solutions:
                    rsp_funcs[(aop, bop,
                               iw)] = -np.dot(va[aop], solutions[(bop, iw)])
            return {
                'response_functions': rsp_funcs,
                'solutions': solutions,
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

    def print_iteration(self, relative_residual_norm, xvs):
        """
        Prints information of the iteration.

        :param relative_residual_norm:
            Relative residual norms.
        :param xvs:
            A list of tuples containing operator component, imaginary
            frequency, and property.
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

        for op, imagfreq, xv in xvs:
            ops_label = '<<{};{}>>_{:.4f}j'.format(op, op, imagfreq)
            rel_res = relative_residual_norm[(op, imagfreq)]
            output_iter = '{:<17s}: {:15.8f} {:15.8f}j   '.format(
                ops_label, -xv.real, -xv.imag)
            output_iter += 'Residual Norm: {:.8f}'.format(rel_res)
            self.ostream.print_header(output_iter.ljust(width))
        self.ostream.print_blank()
        self.ostream.flush()

    def print_header(self):
        """
        Prints C6 solver setup header to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header('C6 Value Response Driver Setup')
        self.ostream.print_header(32 * '=')
        self.ostream.print_blank()

        width = 60

        cur_str = 'Number of integration points    : ' + str(self.n_points)
        self.ostream.print_header(cur_str.ljust(width))

        cur_str = 'Max. Number of Iterations       : ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = 'Convergence Threshold           : ' + \
            '{:.1e}'.format(self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(width))

        cur_str = 'ERI Screening Scheme            : ' + get_qq_type(
            self.qq_type)
        self.ostream.print_header(cur_str.ljust(width))
        cur_str = 'ERI Screening Threshold         : ' + \
            '{:.1e}'.format(self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(width))

        if self.dft:
            cur_str = 'Exchange-Correlation Functional : '
            cur_str += self.xcfun.get_func_label().upper()
            self.ostream.print_header(cur_str.ljust(width))
            cur_str = 'Molecular Grid Level            : ' + str(
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
        Prints timing for the C6 solver.
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
