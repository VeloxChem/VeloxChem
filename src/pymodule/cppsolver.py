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
from .lrmatvecdriver import lrvec2mat
from .lrmatvecdriver import get_complex_rhs
from .lrmatvecdriver import check_rsp_hdf5
from .lrmatvecdriver import write_rsp_hdf5
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_scheme
from .qqscheme import get_qq_type
from .inputparser import parse_frequencies


class ComplexResponse:
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
        - frequencies: The frequencies.
        - damping: The damping parameter.
        - qq_type: The electron repulsion integrals screening scheme.
        - eri_thresh: The electron repulsion integrals screening threshold.
        - batch_size: The batch size for computation of Fock matrices.
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
        - nonlinear: The flag for returning kappa values.
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
        Initializes complex linear response solver to default setup.
        """

        self.a_operator = 'dipole'
        self.a_components = 'xyz'
        self.b_operator = 'dipole'
        self.b_components = 'xyz'

        self.frequencies = (0.1,)
        self.damping = 0.004556335294880438

        self.qq_type = 'QQ_DEN'
        self.eri_thresh = 1.0e-15
        self.batch_size = None

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
        self.nonlinear = False

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
        if 'response' in rsp_dict:
            key = rsp_dict['response'].lower()
            self.nonlinear = True if key == 'nonlinear' else False

        if 'conv_thresh' in rsp_dict:
            self.conv_thresh = float(rsp_dict['conv_thresh'])
        if 'max_iter' in rsp_dict:
            self.max_iter = int(rsp_dict['max_iter'])

        if 'eri_thresh' in rsp_dict:
            self.eri_thresh = float(rsp_dict['eri_thresh'])
        if 'qq_type' in rsp_dict:
            self.qq_type = rsp_dict['qq_type'].upper()
        if 'batch_size' in rsp_dict:
            self.batch_size = int(rsp_dict['batch_size'])

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
        Decomposes trial vectors into their 4 respective parts (real gerade,
        real ungerade, imaginary gerade, and imaginary ungerade).

        :param vecs:
            The trial vectors.

        :return:
            A tuple containing respective parts of the trial vectors.
        """

        quarter = vecs.shape[0] // 4
        quarter_2 = quarter * 2
        quarter_3 = quarter * 3

        if vecs.ndim == 2:
            return (
                vecs[:quarter, :],
                vecs[quarter:quarter_2, :],
                vecs[quarter_2:quarter_3, :],
                vecs[quarter_3:, :],
            )

        elif vecs.ndim == 1:
            return (
                vecs[:quarter],
                vecs[quarter:quarter_2],
                vecs[quarter_2:quarter_3],
                vecs[quarter_3:],
            )

        else:
            return (None, None, None, None)

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

        if realvec.ndim == 2:
            for vec in range(realvec.shape[1]):
                if np.linalg.norm(realvec[:, vec]) > self.small_thresh:
                    space.append(realvec[:, vec])
                if np.linalg.norm(imagvec[:, vec]) > self.small_thresh:
                    space.append(imagvec[:, vec])

        elif realvec.ndim == 1:
            if np.linalg.norm(realvec) > self.small_thresh:
                space.append(realvec)
            if np.linalg.norm(imagvec) > self.small_thresh:
                space.append(imagvec)

        return np.array(space).T

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

        ediag, sdiag = construct_ed_sd_half(orb_ene, nocc, norb)
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

        return (pa_diag, pb_diag, pc_diag, pd_diag)

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

        pa, pb, pc, pd = precond

        v_in_rg, v_in_ru, v_in_iu, v_in_ig = self.decomp_trials(v_in)

        v_out_rg = pa * v_in_rg + pb * v_in_ru + pc * v_in_iu + pd * v_in_ig
        v_out_ru = pb * v_in_rg + pa * v_in_ru + pd * v_in_iu + pc * v_in_ig
        v_out_iu = pc * v_in_rg + pd * v_in_ru - pa * v_in_iu - pb * v_in_ig
        v_out_ig = pd * v_in_rg + pc * v_in_ru - pb * v_in_iu - pa * v_in_ig

        return np.hstack((v_out_rg, v_out_ru, v_out_iu, v_out_ig))

    def initial_guess(self, v1, d, precond):
        """
        Creating initial guess (un-orthonormalized trials) out of gradients.

        :param v1:
            The dictionary containing (operator, frequency) as keys and
            right-hand sides as values.
        :param d:
            The damping parameter.
        :param precond:
            The preconditioner matrices.

        :return:
            The initial guess.
        """

        ig = {}
        for (op, w), grad in v1.items():
            gradger, gradung = self.decomp_grad(grad)

            grad = np.hstack(
                (gradger.real, gradung.real, -gradung.imag, -gradger.imag))
            gn = np.sqrt(2.0) * np.linalg.norm(grad)

            if gn < self.small_thresh:
                ig[(op, w)] = np.zeros(grad.shape[0])
            else:
                ig[(op, w)] = self.preconditioning(precond[w], grad)

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
        :param res_norm:
            The relative residual norm.
        :param renormalize:
            The flag for normalization.

        :return:
            The orthonormalized gerade and ungerade trial vectors.
        """

        if self.rank == mpi_master():
            trials = []
            for (op, w) in vectors:
                vec = np.array(vectors[(op, w)])

                # preconditioning trials:

                if precond is not None:
                    v = self.preconditioning(precond[w], vec)
                else:
                    v = vec
                if np.linalg.norm(v) * np.sqrt(2.0) > self.small_thresh:
                    trials.append(v)

            new_trials = np.array(trials).T

            # decomposing the full space trial vectors...

            (new_realger, new_realung, new_imagung,
             new_imagger) = self.decomp_trials(new_trials)

            # ...and assembling gerade and ungerade subspaces

            new_ger = self.assemble_subsp(new_realger, new_imagger)
            new_ung = self.assemble_subsp(new_realung, new_imagung)

        else:
            new_ger, new_ung = None, None

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
            A dictionary containing response functions, solutions and a
            dictionarry containing solutions and kappa values when called from
            a non-linear response module.
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
                                self.xcfun, self.pe, self.potfile,
                                self.batch_size)
        timing_dict = {}

        if self.rank == mpi_master():
            nocc = molecule.number_of_alpha_electrons()
            norb = scf_tensors['C'].shape[1]
            orb_ene = scf_tensors['E']
        else:
            nocc = None

        valid_v1 = False
        if self.rank == mpi_master():
            valid_v1 = (v1 is not None)
        valid_v1 = self.comm.bcast(valid_v1, root=mpi_master())

        if not valid_v1:
            b_rhs = get_complex_rhs(self.b_operator, self.b_components,
                                    molecule, basis, scf_tensors, self.rank,
                                    self.comm)
            if self.rank == mpi_master():
                v1 = {(op, w): v for op, v in zip(self.b_components, b_rhs)
                      for w in self.frequencies}
        else:
            self.nonlinear = True

        if self.rank == mpi_master():
            op_freq_keys = list(v1.keys())
        else:
            op_freq_keys = None
        op_freq_keys = self.comm.bcast(op_freq_keys, root=mpi_master())

        # creating the preconditioner matrix
        d = self.damping
        if self.rank == mpi_master():
            freqs = set([w for (op, w) in op_freq_keys])
            precond = {
                w: self.get_precond(orb_ene, nocc, norb, w, d) for w in freqs
            }
        else:
            precond = None

        rsp_vector_labels = [
            'CLR_bger_half_size',
            'CLR_bung_half_size',
            'CLR_e2bger_half_size',
            'CLR_e2bung_half_size',
        ]

        # check validity of checkpoint file
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
                self.checkpoint_file, 'CLR_bger_half_size', self.comm)
            dist_bung = DistributedArray.read_from_hdf5_file(
                self.checkpoint_file, 'CLR_bung_half_size', self.comm)
            dist_e2bger = DistributedArray.read_from_hdf5_file(
                self.checkpoint_file, 'CLR_e2bger_half_size', self.comm)
            dist_e2bung = DistributedArray.read_from_hdf5_file(
                self.checkpoint_file, 'CLR_e2bung_half_size', self.comm)

        # generate initial guess from scratch
        else:
            if self.rank == mpi_master():
                igs = self.initial_guess(v1, d, precond)
                bger, bung = self.setup_trials(igs)

                assert_msg_critical(
                    bger.any() or bung.any(),
                    'ComplexResponseSolver: trial vectors are empty')

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

            for op, w in op_freq_keys:
                if (iteration == 0 or
                        relative_residual_norm[(op, w)] > self.conv_thresh):

                    if self.rank == mpi_master():
                        grad = v1[(op, w)]
                        gradger, gradung = self.decomp_grad(grad)
                        gradger_real = gradger.real
                        gradger_imag = gradger.imag
                        gradung_real = gradung.real
                        gradung_imag = gradung.imag
                    else:
                        gradger_real, gradger_imag = None, None
                        gradung_real, gradung_imag = None, None
                    dist_grad_realger = DistributedArray(
                        gradger_real, self.comm)
                    dist_grad_imagger = DistributedArray(
                        gradger_imag, self.comm)
                    dist_grad_realung = DistributedArray(
                        gradung_real, self.comm)
                    dist_grad_imagung = DistributedArray(
                        gradung_imag, self.comm)

                    # projections onto gerade and ungerade subspaces:

                    g_realger = dist_bger.matmul_AtB(dist_grad_realger, 2.0)
                    g_imagger = dist_bger.matmul_AtB(dist_grad_imagger, 2.0)
                    g_realung = dist_bung.matmul_AtB(dist_grad_realung, 2.0)
                    g_imagung = dist_bung.matmul_AtB(dist_grad_imagung, 2.0)

                    # creating gradient and matrix for linear equation

                    size = 2 * (n_ger + n_ung)

                    if self.rank == mpi_master():

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

                        mat[n_ger:n_ger + n_ung,
                            n_ger:n_ger + n_ung] = e2uu[:, :]

                        mat[n_ger + n_ung:size - n_ger,
                            n_ger + n_ung:size - n_ger] = -e2uu[:, :]

                        # filling S2ug

                        mat[n_ger:n_ger + n_ung, :n_ger] = -w * s2ug[:, :]

                        mat[n_ger + n_ung:size - n_ger, :n_ger] = d * s2ug[:, :]

                        mat[n_ger:n_ger + n_ung, size - n_ger:] = d * s2ug[:, :]

                        mat[n_ger + n_ung:size - n_ger,
                            size - n_ger:] = w * s2ug[:, :]

                        # filling S2ug.T (interchanging of row and col)

                        mat[:n_ger, n_ger:n_ger + n_ung] = -w * s2ug.T[:, :]

                        mat[:n_ger,
                            n_ger + n_ung:size - n_ger] = d * s2ug.T[:, :]

                        mat[size - n_ger:,
                            n_ger:n_ger + n_ung] = d * s2ug.T[:, :]

                        mat[size - n_ger:,
                            n_ger + n_ung:size - n_ger] = w * s2ug.T[:, :]

                        # solving matrix equation

                        c = np.linalg.solve(mat, g)
                    else:
                        c = None
                    c = self.comm.bcast(c, root=mpi_master())

                    # extracting the 4 components of c...

                    c_realger = c[:n_ger]
                    c_realung = c[n_ger:n_ger + n_ung]
                    c_imagung = c[n_ger + n_ung:size - n_ger]
                    c_imagger = c[size - n_ger:]

                    # ...and projecting them onto respective subspace

                    x_realger = dist_bger.matmul_AB(c_realger)
                    x_realung = dist_bung.matmul_AB(c_realung)
                    x_imagger = dist_bger.matmul_AB(c_imagger)
                    x_imagung = dist_bung.matmul_AB(c_imagung)

                    # composing full size response vector

                    if self.rank == mpi_master():

                        x_realger_full = np.hstack((x_realger, x_realger))
                        x_realung_full = np.hstack((x_realung, -x_realung))
                        x_imagung_full = np.hstack((x_imagung, -x_imagung))
                        x_imagger_full = np.hstack((x_imagger, x_imagger))

                        x_real = x_realger_full + x_realung_full
                        x_imag = x_imagung_full + x_imagger_full
                        x = x_real + 1j * x_imag

                    # composing E2 and S2 matrices projected onto
                    # solution subspace

                    e2realger = dist_e2bger.matmul_AB(c_realger)
                    e2imagger = dist_e2bger.matmul_AB(c_imagger)
                    e2realung = dist_e2bung.matmul_AB(c_realung)
                    e2imagung = dist_e2bung.matmul_AB(c_imagung)

                    # calculating the residual components

                    if self.rank == mpi_master():

                        s2realger = 2.0 * x_realger
                        s2imagger = 2.0 * x_imagger
                        s2realung = 2.0 * x_realung
                        s2imagung = 2.0 * x_imagung

                        r_realger = (e2realger - w * s2realung + d * s2imagung -
                                     gradger.real)
                        r_realung = (e2realung - w * s2realger + d * s2imagger -
                                     gradung.real)
                        r_imagung = (-e2imagung + w * s2imagger +
                                     d * s2realger + gradung.imag)
                        r_imagger = (-e2imagger + w * s2imagung +
                                     d * s2realung + gradger.imag)

                        # composing total half-sized residual

                        r = np.hstack(
                            (r_realger, r_realung, r_imagung, r_imagger))

                        # calculating relative residual norm
                        # for convergence check

                        xv = np.matmul(x, grad)
                        xvs.append((op, w, xv))

                        rn = np.sqrt(2.0) * np.linalg.norm(r)
                        xn = np.linalg.norm(x)
                        if xn != 0:
                            relative_residual_norm[(op, w)] = rn / xn
                        else:
                            relative_residual_norm[(op, w)] = rn

                        if relative_residual_norm[(op, w)] < self.conv_thresh:
                            solutions[(op, w)] = x
                        else:
                            residuals[(op, w)] = r

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

            # update trial vectors
            new_trials_ger, new_trials_ung = self.setup_trials(
                residuals, precond, dist_bger, dist_bung)

            residuals.clear()

            if self.rank == mpi_master():

                assert_msg_critical(
                    new_trials_ger.any() or new_trials_ung.any(),
                    'ComplexResponseSolver: unable to add new trial vectors')

                if new_trials_ger is None or not new_trials_ger.any():
                    new_trials_ger = np.zeros((new_trials_ung.shape[0], 0))
                if new_trials_ung is None or not new_trials_ung.any():
                    new_trials_ung = np.zeros((new_trials_ger.shape[0], 0))

                # creating new sigma and rho linear transformations

            profiler.stop_timer(iteration, 'Orthonorm.')
            profiler.start_timer(iteration, 'FockBuild')

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

                dt = 0.0
                dt += dist_bger.append_to_hdf5_file(self.checkpoint_file,
                                                    'CLR_bger_half_size')
                dt += dist_bung.append_to_hdf5_file(self.checkpoint_file,
                                                    'CLR_bung_half_size')
                dt += dist_e2bger.append_to_hdf5_file(self.checkpoint_file,
                                                      'CLR_e2bger_half_size')
                dt += dist_e2bung.append_to_hdf5_file(self.checkpoint_file,
                                                      'CLR_e2bung_half_size')
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

        if not self.nonlinear:
            a_rhs = get_complex_rhs(self.a_operator, self.a_components,
                                    molecule, basis, scf_tensors, self.rank,
                                    self.comm)

            if self.rank == mpi_master() and self.is_converged:
                va = {op: v for op, v in zip(self.a_components, a_rhs)}
                rsp_funcs = {}
                for aop in self.a_components:
                    for bop, w in solutions:
                        rsp_funcs[(aop, bop,
                                   w)] = -np.dot(va[aop], solutions[(bop, w)])
                return {'response_functions': rsp_funcs, 'solutions': solutions}
            else:
                return {}

        else:
            if self.rank == mpi_master() and self.is_converged:
                kappas = {}
                for op, w in solutions:
                    kappas[(op, w)] = (
                        lrvec2mat(solutions[(op, w)].real, nocc, norb) +
                        1j * lrvec2mat(solutions[(op, w)].imag, nocc, norb))
                return {'solutions': solutions, 'kappas': kappas}
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

        for op, freq, xv in xvs:
            ops_label = '<<{};{}>>_{:.4f}'.format(op, op, freq)
            rel_res = relative_residual_norm[(op, freq)]
            output_iter = '{:<15s}: {:15.8f} {:15.8f}j   '.format(
                ops_label, -xv.real, -xv.imag)
            output_iter += 'Residual Norm: {:.8f}'.format(rel_res)
            self.ostream.print_header(output_iter.ljust(width))
        self.ostream.print_blank()
        self.ostream.flush()

    def print_header(self):
        """
        Prints complex linear response solver setup header to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Complex Response Driver Setup')
        self.ostream.print_header(31 * '=')
        self.ostream.print_blank()

        width = 60

        cur_str = 'Damping Parameter (gamma)       : {:.6e}'.format(
            self.damping)
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
        if self.batch_size is not None:
            cur_str = 'Batch Size of Fock Matrices     : ' + \
                '{:d}'.format(self.batch_size)
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
        self.ostream.print_blank()
