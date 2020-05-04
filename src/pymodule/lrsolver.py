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
from .profiler import Profiler
from .distributedarray import DistributedArray
from .aodensitymatrix import AODensityMatrix
from .signalhandler import SignalHandler
from .lrmatvecdriver import LinearResponseMatrixVectorDriver
from .lrmatvecdriver import remove_linear_dependence_half
from .lrmatvecdriver import orthogonalize_gram_schmidt_half
from .lrmatvecdriver import normalize_half
from .lrmatvecdriver import construct_ed_sd_half
from .lrmatvecdriver import lrvec2mat
from .lrmatvecdriver import get_rhs
from .lrmatvecdriver import check_rsp_hdf5
from .lrmatvecdriver import write_rsp_hdf5
from .qqscheme import get_qq_scheme
from .qqscheme import get_qq_type
from .errorhandler import assert_msg_critical
from .inputparser import parse_frequencies


class LinearResponseSolver:
    """
    Implements linear response solver.

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
        - eri_thresh: The electron repulsion integrals screening threshold.
        - qq_type: The electron repulsion integrals screening scheme.
        - dft: The flag for running DFT.
        - grid_level: The accuracy level of DFT grid.
        - xcfun: The XC functional.
        - pe: The flag for running polarizable embedding calculation.
        - potfile: The name of the potential file for polarizable embedding.
        - use_split_comm: The flag for using split communicators.
        - conv_thresh: The convergence threshold for the solver.
        - max_iter: The maximum number of solver iterations.
        - cur_iter: Index of the current iteration.
        - small_thresh: The norm threshold for a vector to be considered a zero
          vector.
        - lindep_thresh: The threshold for removing linear dependence in the
          trial vectors.
        - is_converged: The flag for convergence.
        - comm: The MPI communicator.
        - rank: The MPI rank.
        - nodes: Number of MPI processes.
        - ostream: The output stream.
        - restart: The flag for restarting from checkpoint file.
        - checkpoint_file: The name of checkpoint file.
        - timing: The flag for printing timing information.
        - profiling: The flag for printing profiling information.
        - memory_profiling: The flag for printing memory usage.
        - memory_tracing: The flag for tracing memory allocation using
    """

    def __init__(self, comm, ostream):
        """
        Initializes linear response solver to default setup.
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

        # split communicators
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

        self.program_start_time = None
        self.maximum_hours = None

        self.timing = False
        self.profiling = False
        self.memory_profiling = False
        self.memory_tracing = False

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

        if 'program_start_time' in rsp_dict:
            self.program_start_time = rsp_dict['program_start_time']
        if 'maximum_hours' in rsp_dict:
            self.maximum_hours = rsp_dict['maximum_hours']

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

    def compute(self, molecule, basis, scf_tensors, v1=None):
        """
        Performs linear response calculation for a molecule and a basis set.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
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

        e2x_drv = LinearResponseMatrixVectorDriver(self.comm,
                                                   self.use_split_comm,
                                                   self.ostream)
        e2x_drv.update_settings(self.eri_thresh, self.qq_type, self.dft,
                                self.xcfun, self.pe, self.potfile)
        timing_dict = {}

        valid_v1 = False
        if self.rank == mpi_master():
            valid_v1 = (v1 is not None)
        valid_v1 = self.comm.bcast(valid_v1, root=mpi_master())

        if not valid_v1:
            nonlinear_flag = False
            b_rhs = get_rhs(self.b_operator, self.b_components, molecule, basis,
                            scf_tensors, self.rank, self.comm)
            if self.rank == mpi_master():
                v1 = {(op, w): v for op, v in zip(self.b_components, b_rhs)
                      for w in self.frequencies}
        else:
            nonlinear_flag = True

        if self.rank == mpi_master():
            op_freq_keys = list(v1.keys())
        else:
            op_freq_keys = None
        op_freq_keys = self.comm.bcast(op_freq_keys, root=mpi_master())

        if self.rank == mpi_master():
            freqs = set([w for (op, w) in op_freq_keys])
            precond = {w: self.get_precond(ea, nocc, norb, w) for w in freqs}
        else:
            precond = None

        rsp_vector_labels = [
            'LR_bger_half_size',
            'LR_bung_half_size',
            'LR_e2bger_half_size',
            'LR_e2bung_half_size',
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
            dist_bger, dist_bung, dist_e2bger, dist_e2bung = [
                DistributedArray.read_from_hdf5_file(self.checkpoint_file,
                                                     label, self.comm)
                for label in rsp_vector_labels
            ]
            checkpoint_text = 'Restarting from checkpoint file: '
            checkpoint_text += self.checkpoint_file
            self.ostream.print_info(checkpoint_text)
            self.ostream.print_blank()

        # generate initial guess from scratch
        else:
            if self.rank == mpi_master():
                igs = self.initial_guess(v1, precond)
                bger, bung = self.setup_trials(igs)

                assert_msg_critical(
                    bger.any() or bung.any(),
                    'LinearResponseSolver.compute: trial vector is empty')

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

        signal_handler = SignalHandler()
        signal_handler.add_sigterm_function(
            self.graceful_exit, molecule, basis, dft_func_label, potfile_text,
            [dist_bger, dist_bung, dist_e2bger, dist_e2bung], rsp_vector_labels)

        # start iterations
        for iteration in range(self.max_iter):

            profiler.start_timer(iteration, 'ReducedSpace')

            n_ger = dist_bger.shape(1)
            n_ung = dist_bung.shape(1)

            e2gg = dist_bger.matmul_AtB(dist_e2bger, 2.0)
            e2uu = dist_bung.matmul_AtB(dist_e2bung, 2.0)
            s2ug = dist_bung.matmul_AtB(dist_bger, 4.0)

            xvs = []
            self.cur_iter = iteration

            for op, freq in op_freq_keys:
                if (iteration > 0 and
                        relative_residual_norm[(op, freq)] < self.conv_thresh):
                    continue

                if self.rank == mpi_master():
                    v = v1[(op, freq)]
                    gradger, gradung = self.decomp_grad(v)
                else:
                    gradger, gradung = None, None
                dist_gradger = DistributedArray(gradger, self.comm)
                dist_gradung = DistributedArray(gradung, self.comm)

                g_ger = dist_bger.matmul_AtB(dist_gradger, 2.0)
                g_ung = dist_bung.matmul_AtB(dist_gradung, 2.0)

                if self.rank == mpi_master():
                    mat = np.zeros((n_ger + n_ung, n_ger + n_ung))
                    mat[:n_ger, :n_ger] = e2gg[:, :]
                    mat[:n_ger, n_ger:] = -freq * s2ug.T[:, :]
                    mat[n_ger:, :n_ger] = -freq * s2ug[:, :]
                    mat[n_ger:, n_ger:] = e2uu[:, :]

                    g = np.zeros(n_ger + n_ung)
                    g[:n_ger] = g_ger[:]
                    g[n_ger:] = g_ung[:]

                    c = np.linalg.solve(mat, g)
                else:
                    c = None
                c = self.comm.bcast(c, root=mpi_master())

                c_ger = c[:n_ger]
                c_ung = c[n_ger:]

                x_ger = dist_bger.matmul_AB(c_ger)
                x_ung = dist_bung.matmul_AB(c_ung)

                e2x_ger = dist_e2bger.matmul_AB(c_ger)
                e2x_ung = dist_e2bung.matmul_AB(c_ung)

                if self.rank == mpi_master():
                    s2x_ger = 2.0 * x_ger
                    s2x_ung = 2.0 * x_ung

                    x_ger_full = np.hstack((x_ger, x_ger))
                    x_ung_full = np.hstack((x_ung, -x_ung))
                    x = x_ger_full + x_ung_full

                    r_ger = e2x_ger - freq * s2x_ung - gradger
                    r_ung = e2x_ung - freq * s2x_ger - gradung
                    r = np.hstack((r_ger, r_ung))

                    xv = np.dot(x, v)
                    xvs.append((op, freq, xv))

                    rn = np.linalg.norm(r) * np.sqrt(2.0)
                    xn = np.linalg.norm(x)
                    if xn != 0:
                        relative_residual_norm[(op, freq)] = rn / xn
                    else:
                        relative_residual_norm[(op, freq)] = rn

                    if relative_residual_norm[(op, freq)] < self.conv_thresh:
                        solutions[(op, freq)] = x
                    else:
                        residuals[(op, freq)] = r

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
                    'LinearResponseSolver: unable to add new trial vector')

                if new_trials_ger is None or not new_trials_ger.any():
                    new_trials_ger = np.zeros((new_trials_ung.shape[0], 0))
                if new_trials_ung is None or not new_trials_ung.any():
                    new_trials_ung = np.zeros((new_trials_ger.shape[0], 0))

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

            if self.maximum_hours is not None:
                remaining_hours = (self.maximum_hours -
                                   (tm.time() - self.program_start_time) / 3600)
                if self.maximum_hours < 10.0:
                    if remaining_hours < 0.1 * self.maximum_hours:
                        signal_handler.raise_signal('SIGTERM')
                else:
                    if remaining_hours < 1.0:
                        signal_handler.raise_signal('SIGTERM')

            profiler.stop_timer(iteration, 'FockBuild')
            if self.dft or self.pe:
                profiler.update_timer(iteration, timing_dict)

            profiler.check_memory_usage(
                'Iteration {:d} sigma build'.format(iteration + 1))

        signal_handler.remove_sigterm_function()

        self.write_checkpoint(molecule, basis, dft_func_label, potfile_text,
                              [dist_bger, dist_bung, dist_e2bger, dist_e2bung],
                              rsp_vector_labels)

        # converged?
        if self.rank == mpi_master():
            self.print_convergence()

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of LR solver')
        profiler.print_memory_usage(self.ostream)

        # calculate response functions
        if not nonlinear_flag:
            a_rhs = get_rhs(self.a_operator, self.a_components, molecule, basis,
                            scf_tensors, self.rank, self.comm)

            if self.rank == mpi_master() and self.is_converged:
                va = {op: v for op, v in zip(self.a_components, a_rhs)}
                rsp_funcs = {}
                for aop in self.a_components:
                    for bop, w in solutions:
                        rsp_funcs[(aop, bop,
                                   w)] = -np.dot(va[aop], solutions[(bop, w)])
                return {
                    'response_functions': rsp_funcs,
                    'solutions': solutions,
                }
            else:
                return {}

        else:
            if self.rank == mpi_master() and self.is_converged:
                kappas = {}
                for (op, freq), x in solutions.items():
                    kappas[(op, freq)] = lrvec2mat(x, nocc, norb)
                return {
                    'solutions': solutions,
                    'kappas': kappas,
                }
            else:
                return {}

    def write_checkpoint(self, molecule, basis, dft_func_label, potfile_text,
                         dist_arrays, labels):
        """
        Writes checkpoint file.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param dft_func_label:
            The name of the density functional.
        :param potfile_text:
            The text in the potential file.
        :param dist_arrays:
            The list of distributed arrays.
        :param labels:
            The list of labels.
        """

        if self.checkpoint_file is None:
            return

        if self.rank == mpi_master():
            success = write_rsp_hdf5(self.checkpoint_file, [], [],
                                     molecule.nuclear_repulsion_energy(),
                                     molecule.elem_ids_to_numpy(),
                                     basis.get_label(), dft_func_label,
                                     potfile_text, self.ostream)
        else:
            success = False
        success = self.comm.bcast(success, root=mpi_master())

        if success:
            for dist_array, label in zip(dist_arrays, labels):
                dist_array.append_to_hdf5_file(self.checkpoint_file, label)

    def graceful_exit(self, molecule, basis, dft_func_label, potfile_text,
                      dist_arrays, labels):
        """
        Gracefully exits SCF driver.

        :param molecule:
            The molecule.
        :param basis:
            The basis set.
        :param dft_func_label:
            The name of the density functional.
        :param potfile_text:
            The text in the potential file.
        :param dist_arrays:
            The list of distributed arrays.
        :param labels:
            The list of labels.

        :return:
            The return code.
        """

        self.ostream.print_blank()
        self.ostream.print_info('Preparing for a graceful termination...')
        self.ostream.print_blank()
        self.ostream.flush()

        self.write_checkpoint(molecule, basis, dft_func_label, potfile_text,
                              dist_arrays, labels)

        self.ostream.print_info('...done.')
        self.ostream.print_blank()
        self.ostream.print_info('Exiting program.')
        self.ostream.print_blank()
        self.ostream.flush()

        return 0

    def print_header(self):
        """
        Prints linear response solver setup header to output stream.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Linear Response Solver Setup')
        self.ostream.print_header(30 * '=')
        self.ostream.print_blank()

        str_width = 60

        cur_str = 'Max. Number of Iterations       : ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold           : ' + \
            '{:.1e}'.format(self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = 'ERI Screening Scheme            : ' + get_qq_type(
            self.qq_type)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'ERI Screening Threshold         : ' + \
            '{:.1e}'.format(self.eri_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        if self.dft:
            cur_str = 'Exchange-Correlation Functional : '
            cur_str += self.xcfun.get_func_label().upper()
            self.ostream.print_header(cur_str.ljust(str_width))
            cur_str = 'Molecular Grid Level            : ' + str(
                self.grid_level)
            self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()

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
        for op, freq, xv in xvs:
            ops_label = '<<{};{}>>_{:.4f}'.format(op, op, freq)
            rel_res = relative_residual_norm[(op, freq)]
            output_iter = '{:<15s}: {:15.8f} '.format(ops_label, -xv)
            output_iter += 'Residual Norm: {:.8f}'.format(rel_res)
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

    def initial_guess(self, v1, precond):
        """
        Creating initial guess for the linear response solver.

        :param v1:
            The dictionary containing (operator, frequency) as keys and
            right-hand sides as values.
        :param precond:
            The preconditioner.

        :return:
            The initial guess.
        """

        ig = {}
        for (op, w), grad in v1.items():
            gradger, gradung = self.decomp_grad(grad)

            grad = np.hstack((gradger, gradung))
            gn = np.linalg.norm(grad) * np.sqrt(2.0)

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

        assert_msg_critical(grad.ndim == 1, 'decomp_grad: Expecting a 1D array')

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

        return (pa_diag, pb_diag)

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

        pa, pb = precond

        v_in_rg, v_in_ru = self.decomp_trials(v_in)

        v_out_rg = pa * v_in_rg + pb * v_in_ru
        v_out_ru = pb * v_in_rg + pa * v_in_ru

        return np.hstack((v_out_rg, v_out_ru))

    def setup_trials(self,
                     vectors,
                     precond=None,
                     dist_bger=None,
                     dist_bung=None,
                     renormalize=True):
        """
        Computes orthonormalized trial vectors.

        :param vectors:
            The set of vectors.
        :param precond:
            The preconditioner.
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

            for (op, freq) in vectors:
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
        else:
            new_ger, new_ung = None, None

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

        if self.rank == mpi_master() and renormalize:
            new_ger = remove_linear_dependence_half(new_ger, self.lindep_thresh)
            new_ger = orthogonalize_gram_schmidt_half(new_ger)
            new_ger = normalize_half(new_ger)

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

        if vecs.ndim == 1:
            ger = vecs[:half_rows]
            ung = vecs[half_rows:]

        elif vecs.ndim == 2:
            ger = vecs[:half_rows, :]
            ung = vecs[half_rows:, :]

        return ger, ung
