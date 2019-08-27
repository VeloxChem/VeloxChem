import numpy as np
import time as tm

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import GridDriver
from .veloxchemlib import MolecularGrid
from .veloxchemlib import XCFunctional
from .veloxchemlib import mpi_master
from .veloxchemlib import parse_xc_func
from .lrmatvecdriver import LinearResponseMatrixVectorDriver
from .lrmatvecdriver import remove_linear_dependence
from .lrmatvecdriver import orthogonalize_gram_schmidt
from .lrmatvecdriver import normalize
from .lrmatvecdriver import construct_ed_sd
from .lrmatvecdriver import get_rhs
from .lrmatvecdriver import swap_xy
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
        self.molgrid = MolecularGrid()

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
            self.xcfun = parse_xc_func(method_dict['xcfun'].upper())

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
            import os
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
            od, sd = construct_ed_sd(ea, nocc, norb)
            td = {w: od - w * sd for w in self.frequencies}
        else:
            nocc = None
            norb = None

        # generate integration grid
        if self.dft:
            grid_drv = GridDriver(self.comm)
            grid_drv.set_level(self.grid_level)

            grid_t0 = tm.time()
            self.molgrid = grid_drv.generate(molecule)
            n_grid_points = self.molgrid.number_of_points()
            self.molgrid.distribute(self.rank, self.nodes, self.comm)
            self.ostream.print_info(
                'Molecular grid with {0:d} points generated in {1:.2f} sec.'.
                format(n_grid_points,
                       tm.time() - grid_t0))
            self.ostream.print_blank()

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        e2x_drv = LinearResponseMatrixVectorDriver(self.comm)

        a_rhs = get_rhs(self.a_operator, self.a_components, molecule, basis,
                        scf_tensors, self.rank, self.comm)
        b_rhs = get_rhs(self.b_operator, self.b_components, molecule, basis,
                        scf_tensors, self.rank, self.comm)

        # start linear response calculation
        if self.rank == mpi_master():
            V1 = {op: v for op, v in zip(self.b_components, b_rhs)}
            igs = self.initial_guess(self.frequencies, V1, od, sd, td)
            b = self.setup_trials(igs)

            assert_msg_critical(
                np.any(b),
                'LinearResponseSolver.compute: trial vector is empty')
        else:
            b = None

        if self.timing:
            self.timing_dict['reduced_space'][0] += tm.time() - timing_t0
            timing_t0 = tm.time()

        e2b = e2x_drv.e2n(b, scf_tensors, screening, molecule, basis)
        if self.rank == mpi_master():
            s2b = e2x_drv.s2n(b, scf_tensors, nocc)

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

                # next solution
                for op, freq in igs:
                    v = V1[op]
                    reduced_solution = np.linalg.solve(
                        np.matmul(b.T, e2b - freq * s2b), np.matmul(b.T, v))
                    solutions[(op, freq)] = np.matmul(b, reduced_solution)
                    residuals[(op, freq)] = np.matmul(e2b - freq * s2b,
                                                      reduced_solution) - v

                    r = residuals[(op, freq)]
                    n = solutions[(op, freq)]

                    nv = np.dot(n, v)
                    nvs.append((op, freq, nv))

                    rn = np.linalg.norm(r)
                    nn = np.linalg.norm(n)
                    relative_residual_norm[(op, freq)] = rn / nn
                    converged[(op, freq)] = (rn / nn < self.conv_thresh)

                # write to output
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
                new_trials = self.setup_trials(residuals, converged, td=td, b=b)
                b = np.append(b, new_trials, axis=1)
            else:
                new_trials = None

            if self.timing:
                tid = iteration + 1
                self.timing_dict['reduced_space'][tid] += tm.time() - timing_t0
                timing_t0 = tm.time()

            new_e2b = e2x_drv.e2n(new_trials, scf_tensors, screening, molecule,
                                  basis)
            if self.rank == mpi_master():
                new_s2b = e2x_drv.s2n(new_trials, scf_tensors, nocc)
                e2b = np.append(e2b, new_e2b, axis=1)
                s2b = np.append(s2b, new_s2b, axis=1)

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

        cur_str = "Max. Number of Iterations : " + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "Convergence Threshold     : " + \
            "{:.1e}".format(self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        cur_str = "ERI Screening Scheme      : " + get_qq_type(self.qq_type)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = "ERI Screening Threshold   : " + \
            "{:.1e}".format(self.eri_thresh)
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

    def initial_guess(self, freqs, V1, od, sd, td):
        """
        Creating initial guess for the linear response solver.

        :param freq:
            The frequencies.
        :param V1:
            The dictionary containing operator components (key) and right-hand
            sides (values).
        :param od:
            The array of E0 diagonal elements.
        :param sd:
            The array of S0 diagonal elements.
        :param td:
            The preconditioner.

        :return:
            The initial guess.
        """

        dim = od.shape[0]

        ig = {}
        for op, grad in V1.items():
            gn = np.linalg.norm(grad)
            for w in freqs:
                if gn < self.small_thresh:
                    ig[(op, w)] = np.zeros(dim)
                else:
                    ig[(op, w)] = grad / td[w]
        return ig

    def setup_trials(self,
                     vectors,
                     converged={},
                     td=None,
                     b=None,
                     renormalize=True):
        """
        Computes orthonormalized trial vectors.

        :param vectors:
            The set of vectors.
        :param converged:
            The flags for converged vectors.
        :param td:
            The preconditioner.
        :param b:
            The subspace.
        :param renormalize:
            The flag for normalization.

        :return:
            The orthonormalized trial vectors.
        """

        trials = []

        for (op, freq) in vectors:
            if converged and converged[(op, freq)]:
                continue

            vec = vectors[(op, freq)]

            if td is not None:
                v = vec / td[freq]
            else:
                v = vec

            if np.linalg.norm(v) > self.small_thresh:
                trials.append(v)
                if freq > self.small_thresh:
                    trials.append(swap_xy(v))

        new_trials = np.array(trials).T

        if b is not None:
            new_trials = new_trials - np.matmul(b, np.matmul(b.T, new_trials))

        if trials and renormalize:
            new_trials = remove_linear_dependence(new_trials,
                                                  self.lindep_thresh)
            new_trials = orthogonalize_gram_schmidt(new_trials)
            new_trials = normalize(new_trials)

        return new_trials

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
