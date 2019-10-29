import numpy as np
import time as tm
import os

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import GridDriver
from .veloxchemlib import MolecularGrid
from .veloxchemlib import XCFunctional
from .veloxchemlib import denmat
from .veloxchemlib import mpi_master
from .veloxchemlib import parse_xc_func
from .aodensitymatrix import AODensityMatrix
from .lrmatvecdriver import LinearResponseMatrixVectorDriver
from .lrmatvecdriver import remove_linear_dependence_half
from .lrmatvecdriver import orthogonalize_gram_schmidt_half
from .lrmatvecdriver import normalize_half
from .lrmatvecdriver import construct_ed_sd_half
from .lrmatvecdriver import get_rhs
from .lrmatvecdriver import read_rsp_hdf5
from .lrmatvecdriver import write_rsp_hdf5
from .qqscheme import get_qq_scheme
from .qqscheme import get_qq_type
from .errorhandler import assert_msg_critical
from .inputparser import parse_frequencies


class LinearResponseSolver:
    """
    Implements linear response solver.

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
    :param eri_thresh:
        The electron repulsion integrals screening threshold.
    :param qq_type:
        The electron repulsion integrals screening scheme.
    :param dft:
        The flag for running DFT.
    :param grid_level:
        The accuracy level of DFT grid.
    :param xcfun:
        The XC functional.
    :param pe:
        The flag for running polarizable embedding calculation.
    :param potfile:
        The name of the potential file for polarizable embedding.
    :param conv_thresh:
        The convergence threshold for the solver.
    :param max_iter:
        The maximum number of solver iterations.
    :param cur_iter:
        Index of the current iteration.
    :param small_thresh:
        The norm threshold for a vector to be considered a zero vector.
    :param lindep_thresh:
        The threshold for removing linear dependence in the trial vectors.
    :param is_converged:
        The flag for convergence.
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
    :param checkpoint_time:
        The timer of checkpoint file.
    :param timing:
        The flag for printing timing information.
    :param profiling:
        The flag for printing profiling information.
    """

    def __init__(self, comm, ostream):
        """
        Initializes linear response solver to default setup.

        :param comm:
            The MPI communicator.
        :param ostream:
            The output stream.
        """

        # operators and frequencies
        self.a_operator = 'dipole'
        self.a_components = 'xyz'
        self.b_operator = 'dipole'
        self.b_components = 'xyz'
        self.frequencies = (0,)

        # ERI settings
        self.eri_thresh = 1.0e-15
        self.qq_type = 'QQ_DEN'

        # dft
        self.dft = False
        self.grid_level = 4
        self.xcfun = XCFunctional()

        # polarizable embedding
        self.pe = False
        self.potfile = None

        self.use_split_comm = False

        # solver setup
        self.conv_thresh = 1.0e-4
        self.max_iter = 50
        self.cur_iter = 0
        self.small_thresh = 1.0e-10
        self.lindep_thresh = 1.0e-6
        self.is_converged = False

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # output stream
        self.ostream = ostream

        # restart information
        self.restart = True
        self.checkpoint_file = None
        self.checkpoint_time = None

        self.timing = False
        self.profiling = False

    def update_settings(self, rsp_dict, method_dict={}):
        """
        Updates response and method settings in linear response solver.

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

        if 'eri_thresh' in rsp_dict:
            self.eri_thresh = float(rsp_dict['eri_thresh'])
        if 'qq_type' in rsp_dict:
            self.qq_type = rsp_dict['qq_type']

        if 'conv_thresh' in rsp_dict:
            self.conv_thresh = float(rsp_dict['conv_thresh'])
        if 'max_iter' in rsp_dict:
            self.max_iter = int(rsp_dict['max_iter'])
        if 'lindep_thresh' in rsp_dict:
            self.lindep_thresh = float(rsp_dict['lindep_thresh'])

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

    def compute(self, molecule, basis, scf_tensors):
        """
        Performs linear response calculation for a molecule and a basis set.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            A dictionary containing properties.
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
                'new_trials': [0.0],
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
            'LinearResponseSolver: not implemented for unrestricted case')

        # make preparations
        if self.rank == mpi_master():
            mo = scf_tensors['C']
            ea = scf_tensors['E']
            nocc = nalpha
            norb = mo.shape[1]

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

        e2x_drv = LinearResponseMatrixVectorDriver(self.comm, self.qq_type,
                                                   self.eri_thresh,
                                                   self.use_split_comm)

        a_rhs = get_rhs(self.a_operator, self.a_components, molecule, basis,
                        scf_tensors, self.rank, self.comm)
        b_rhs = get_rhs(self.b_operator, self.b_components, molecule, basis,
                        scf_tensors, self.rank, self.comm)

        if self.rank == mpi_master():
            V1 = {op: v for op, v in zip(self.b_components, b_rhs)}
            op_freq_keys = [(op, w) for op in V1 for w in self.frequencies]
            precond = {
                w: self.get_precond(ea, nocc, norb, w) for w in self.frequencies
            }

        rsp_vector_labels = [
            'LR_bger_half_size',
            'LR_bung_half_size',
            'LR_e2bger_half_size',
            'LR_e2bung_half_size',
        ]

        bger = None
        bung = None
        new_trials_ger = None
        new_trials_ung = None

        # read initial guess from restart file
        if self.restart:
            if self.rank == mpi_master():
                bger, bung, e2bger, e2bung = read_rsp_hdf5(
                    self.checkpoint_file, rsp_vector_labels,
                    molecule.nuclear_repulsion_energy(),
                    molecule.elem_ids_to_numpy(), basis.get_label(),
                    dft_func_label, potfile_text, self.ostream)
                self.restart = (bger is not None and bung is not None and
                                e2bger is not None and e2bung is not None)
            self.restart = self.comm.bcast(self.restart, root=mpi_master())

        # generate initial guess from scratch
        if not self.restart:
            if self.rank == mpi_master():

                igs = self.initial_guess(self.frequencies, V1, precond)
                bger, bung = self.setup_trials(igs)

                if self.timing:
                    elapsed_time = tm.time() - timing_t0
                    self.timing_dict['reduced_space'][0] += elapsed_time
                    timing_t0 = tm.time()

                assert_msg_critical(
                    bger.any() or bung.any(),
                    'LinearResponseSolver.compute: trial vector is empty')

                if bger is None or not bger.any():
                    bger = np.zeros((bung.shape[0], 0))
                if bung is None or not bung.any():
                    bung = np.zeros((bger.shape[0], 0))

            e2bger, e2bung = e2x_drv.e2n_half_size(bger, bung, scf_tensors,
                                                   screening, molecule, basis,
                                                   self.dft, self.xcfun,
                                                   molgrid, gs_density, self.pe,
                                                   self.potfile, V_es, pe_drv)

        solutions = {}
        residuals = {}
        relative_residual_norm = {}
        converged = {}

        if self.timing:
            self.timing_dict['new_trials'][0] += tm.time() - timing_t0
            timing_t0 = tm.time()

        # start iterations
        for iteration in range(self.max_iter):

            if self.timing:
                self.timing_dict['reduced_space'].append(0.0)
                self.timing_dict['new_trials'].append(0.0)

            if self.rank == mpi_master():
                self.cur_iter = iteration
                nvs = []

                n_ger = bger.shape[1]
                n_ung = bung.shape[1]

                e2gg = np.matmul(bger.T, e2bger) * 2.0
                e2uu = np.matmul(bung.T, e2bung) * 2.0
                s2ug = np.matmul(bung.T, bger) * 4.0

                # next solution
                for op, freq in op_freq_keys:
                    v = V1[op]

                    gradger, gradung = self.decomp_grad(v)

                    g_ger = np.matmul(bger.T, gradger) * 2.0
                    g_ung = np.matmul(bung.T, gradung) * 2.0

                    mat = np.zeros((n_ger + n_ung, n_ger + n_ung))
                    mat[:n_ger, :n_ger] = e2gg[:, :]
                    mat[:n_ger, n_ger:] = -freq * s2ug.T[:, :]
                    mat[n_ger:, :n_ger] = -freq * s2ug[:, :]
                    mat[n_ger:, n_ger:] = e2uu[:, :]

                    g = np.zeros(n_ger + n_ung)
                    g[:n_ger] = g_ger[:]
                    g[n_ger:] = g_ung[:]

                    c = np.linalg.solve(mat, g)

                    c_ger = c[:n_ger]
                    c_ung = c[n_ger:]

                    x_ger = np.matmul(bger, c_ger)
                    x_ung = np.matmul(bung, c_ung)

                    x_ger_full = np.hstack((x_ger, x_ger))
                    x_ung_full = np.hstack((x_ung, -x_ung))

                    solutions[(op, freq)] = x_ger_full + x_ung_full

                    r_ger = np.matmul(e2bger, c_ger) - freq * 2.0 * np.matmul(
                        bung, c_ung) - gradger
                    r_ung = np.matmul(e2bung, c_ung) - freq * 2.0 * np.matmul(
                        bger, c_ger) - gradung

                    residuals[(op, freq)] = np.array([r_ger, r_ung]).flatten()

                    r = residuals[(op, freq)]
                    n = solutions[(op, freq)]

                    nv = np.dot(n, v)
                    nvs.append((op, freq, nv))

                    rn = np.linalg.norm(r) * np.sqrt(2.0)
                    nn = np.linalg.norm(n)
                    relative_residual_norm[(op, freq)] = rn / nn

                    converged[(op, freq)] = (rn / nn < self.conv_thresh)

                # write to output
                self.ostream.print_info(
                    '{:d} gerade trial vectors in reduced space'.format(n_ger))
                self.ostream.print_info(
                    '{:d} ungerade trial vectors in reduced space'.format(
                        n_ung))
                self.ostream.print_blank()

                self.print_iteration(relative_residual_norm, converged, nvs)

            if self.timing:
                tid = iteration + 1
                self.timing_dict['reduced_space'][tid] += tm.time() - timing_t0
                timing_t0 = tm.time()

            # check convergence
            self.check_convergence(relative_residual_norm)

            if self.is_converged:
                break

            # update trial vectors
            if self.rank == mpi_master():
                new_trials_ger, new_trials_ung = self.setup_trials(
                    residuals, converged, precond, bger, bung)

                assert_msg_critical(
                    new_trials_ger.any() or new_trials_ung.any(),
                    'LinearResponseSolver: unable to add new trial vector')

                if new_trials_ger is None or not new_trials_ger.any():
                    new_trials_ger = np.zeros((new_trials_ung.shape[0], 0))
                if new_trials_ung is None or not new_trials_ung.any():
                    new_trials_ung = np.zeros((new_trials_ger.shape[0], 0))

                bger = np.append(bger, new_trials_ger, axis=1)
                bung = np.append(bung, new_trials_ung, axis=1)

            if self.timing:
                tid = iteration + 1
                self.timing_dict['reduced_space'][tid] += tm.time() - timing_t0
                timing_t0 = tm.time()

            new_e2bger, new_e2bung = e2x_drv.e2n_half_size(
                new_trials_ger, new_trials_ung, scf_tensors, screening,
                molecule, basis, self.dft, self.xcfun, molgrid, gs_density,
                self.pe, self.potfile, V_es, pe_drv)

            if self.rank == mpi_master():
                e2bger = np.append(e2bger, new_e2bger, axis=1)
                e2bung = np.append(e2bung, new_e2bung, axis=1)

                if tm.time() - self.checkpoint_time > 900.0:
                    write_rsp_hdf5(self.checkpoint_file,
                                   [bger, bung, e2bger, e2bung],
                                   rsp_vector_labels,
                                   molecule.nuclear_repulsion_energy(),
                                   molecule.elem_ids_to_numpy(),
                                   basis.get_label(), dft_func_label,
                                   potfile_text, self.ostream)
                    self.checkpoint_time = tm.time()

            if self.timing:
                tid = iteration + 1
                self.timing_dict['new_trials'][tid] += tm.time() - timing_t0
                timing_t0 = tm.time()

        # converged?
        if self.rank == mpi_master():
            self.print_convergence()

            assert_msg_critical(
                self.is_converged,
                'LinearResponseSolver.compute: failed to converge')

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

        # calculate properties
        if self.rank == mpi_master():
            v1 = {op: v for op, v in zip(self.a_components, a_rhs)}
            lrs = {}
            for aop in self.a_components:
                for bop, w in solutions:
                    lrs[(aop, bop, w)] = -np.dot(v1[aop], solutions[(bop, w)])
            return lrs
        else:
            return None

    def print_header(self):
        """
        Prints linear response solver setup header to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header("Linear Response Solver Setup")
        self.ostream.print_header(30 * "=")
        self.ostream.print_blank()

        str_width = 60

        cur_str = "Max. Number of Iterations       : " + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "Convergence Threshold           : " + \
            "{:.1e}".format(self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = "ERI Screening Scheme            : " + get_qq_type(
            self.qq_type)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "ERI Screening Threshold         : " + \
            "{:.1e}".format(self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        if self.dft:
            cur_str = "Exchange-Correlation Functional : "
            cur_str += self.xcfun.get_func_label().upper()
            self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = "Molecular Grid Level            : " + str(
                self.grid_level)
            self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()

    def print_iteration(self, relative_residual_norm, converged, nvs):
        """
        Prints information of the iteration.

        :param relative_residual_norm:
            Relative residual norms.
        :param converge:
            Flags for converged vectors.
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
        for op, freq, nv in nvs:
            ops_label = '<<{};{}>>_{:.4f}'.format(op, op, freq)
            rel_res = relative_residual_norm[(op, freq)]
            output_iter = '{:<15s}: {:15.8f} '.format(ops_label, -nv)
            output_iter += 'Residual Norm: {:.8f}'.format(rel_res)
            if converged[(op, freq)]:
                output_iter += '   converged'
            self.ostream.print_header(output_iter.ljust(width))
        self.ostream.print_blank()
        self.ostream.flush()

    def print_convergence(self):
        """
        Prints information after convergence.
        """

        width = 92
        output_conv = '*** '
        if self.is_converged:
            output_conv += 'Linear response converged'
        else:
            output_conv += 'Linear response NOT converged'
        output_conv += ' in {:d} iterations. '.format(self.cur_iter + 1)
        output_conv += 'Time: {:.2f} sec'.format(tm.time() - self.start_time)
        self.ostream.print_header(output_conv.ljust(width))
        self.ostream.print_blank()

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

    def initial_guess(self, freqs, V1, precond):
        """
        Creating initial guess for the linear response solver.

        :param freq:
            The frequencies.
        :param V1:
            The dictionary containing operator components (key) and right-hand
            sides (values).
        :param precond:
            The preconditioner.

        :return:
            The initial guess.
        """

        ig = {}
        for op, grad in V1.items():
            gradger, gradung = self.decomp_grad(grad)

            grad = np.array([gradger, gradung]).flatten()
            gn = np.linalg.norm(grad) * np.sqrt(2.0)

            for w in freqs:
                if gn < self.small_thresh:
                    ig[(op, w)] = np.zeros(grad.shape[0])
                else:
                    ig[(op, w)] = self.preconditioning(precond[w], grad)
        return ig

    def decomp_grad(self, grad):
        """
        Decomposes gradient into gerade and ungerade parts.

        :param grad:
            The gradient.

        :return:
            A tuple containing gerade and ungerade parts of gradient.
        """

        assert_msg_critical(
            len(grad.shape) == 1, 'decomp_grad: Expecting a 1D array')

        assert_msg_critical(grad.shape[0] % 2 == 0,
                            'decomp_grad: size of array should be even')

        half_size = grad.shape[0] // 2

        grad_T = np.zeros(grad.shape)
        grad_T[:half_size] = grad[half_size:]
        grad_T[half_size:] = grad[:half_size]

        ger = 0.5 * (grad + grad_T)[:half_size]
        ung = 0.5 * (grad - grad_T)[:half_size]

        return ger.T, ung.T

    def get_precond(self, orb_ene, nocc, norb, w):
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

        :return:
            The preconditioner matrix.
        """

        # spawning needed components

        ediag, sdiag = construct_ed_sd_half(orb_ene, nocc, norb)

        ediag_sq = ediag**2
        sdiag_sq = sdiag**2
        w_sq = w**2

        # constructing matrix block diagonals

        pa_diag = ediag / (ediag_sq - w_sq * sdiag_sq)
        pb_diag = (w * sdiag) / (ediag_sq - w_sq * sdiag_sq)

        precond = np.array([pa_diag, pb_diag])

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

        pa, pb = precond[0], precond[1]

        v_in_rg, v_in_ru = self.decomp_trials(v_in)

        v_out_rg = pa * v_in_rg + pb * v_in_ru
        v_out_ru = pb * v_in_rg + pa * v_in_ru

        v_out = np.array([v_out_rg, v_out_ru]).flatten()

        return v_out

    def setup_trials(self,
                     vectors,
                     converged={},
                     precond=None,
                     bger=None,
                     bung=None,
                     renormalize=True):
        """
        Computes orthonormalized trial vectors.

        :param vectors:
            The set of vectors.
        :param converged:
            The flags for converged vectors.
        :param precond:
            The preconditioner.
        :param bger:
            The gerade subspace.
        :param bung:
            The ungerade subspace.
        :param renormalize:
            The flag for normalization.

        :return:
            The orthonormalized gerade and ungerade trial vectors.
        """

        trials = []

        for (op, freq) in vectors:
            if converged and converged[(op, freq)]:
                continue

            vec = vectors[(op, freq)]

            if precond is not None:
                v = self.preconditioning(precond[freq], vec)
            else:
                v = vec

            if np.linalg.norm(v) * np.sqrt(2.0) > self.small_thresh:
                trials.append(v)

        new_trials = np.array(trials).T

        # decomposing the full space trial vectors...

        new_ger, new_ung = self.decomp_trials(new_trials)

        if bger is not None and bger.any():
            new_ger_proj = np.matmul(bger, 2.0 * np.matmul(bger.T, new_ger))
            new_ger = new_ger - new_ger_proj

        if bung is not None and bung.any():
            new_ung_proj = np.matmul(bung, 2.0 * np.matmul(bung.T, new_ung))
            new_ung = new_ung - new_ung_proj

        if new_ger.any() and renormalize:
            new_ger = remove_linear_dependence_half(new_ger, self.lindep_thresh)
            new_ger = orthogonalize_gram_schmidt_half(new_ger)
            new_ger = normalize_half(new_ger)

        if new_ung.any() and renormalize:
            new_ung = remove_linear_dependence_half(new_ung, self.lindep_thresh)
            new_ung = orthogonalize_gram_schmidt_half(new_ung)
            new_ung = normalize_half(new_ung)

        return new_ger, new_ung

    def decomp_trials(self, vecs):
        """
        Decomposes trial vectors into gerade and ungerade parts.

        :param vecs:
            The trial vectors.

        :return:
            A tuple containing gerade and ungerade parts of the trial vectors.
        """

        assert_msg_critical(vecs.shape[0] % 2 == 0,
                            'decomp_trials: shape[0] of array should be even')

        ger, ung = None, None
        half_rows = vecs.shape[0] // 2

        if len(vecs.shape) == 1:
            ger = vecs[:half_rows]
            ung = vecs[half_rows:]

        elif len(vecs.shape) == 2:
            ger = vecs[:half_rows, :]
            ung = vecs[half_rows:, :]

        return ger, ung

    def print_timing(self):
        """
        Prints timing for the linear response eigensolver.
        """

        width = 92

        valstr = 'Timing (in sec):'
        self.ostream.print_header(valstr.ljust(width))
        self.ostream.print_header(('-' * len(valstr)).ljust(width))

        valstr = '{:<15s} {:>15s} {:>18s}'.format('', 'ReducedSpace',
                                                  'NewTrialVectors')
        self.ostream.print_header(valstr.ljust(width))

        for i, (a, b) in enumerate(
                zip(self.timing_dict['reduced_space'],
                    self.timing_dict['new_trials'])):
            if i == 0:
                title = 'Initial guess'
            else:
                title = 'Iteration {:<5d}'.format(i)
            valstr = '{:<15s} {:15.3f} {:18.3f}'.format(title, a, b)
            self.ostream.print_header(valstr.ljust(width))

        valstr = '---------'
        self.ostream.print_header(valstr.ljust(width))

        valstr = '{:<15s} {:15.3f} {:18.3f}'.format(
            'Sum', sum(self.timing_dict['reduced_space']),
            sum(self.timing_dict['new_trials']))
        self.ostream.print_header(valstr.ljust(width))

        self.ostream.print_blank()
