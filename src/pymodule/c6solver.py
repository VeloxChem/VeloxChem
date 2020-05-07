import numpy as np
import time as tm

from .veloxchemlib import mpi_master
from .profiler import Profiler
from .distributedarray import DistributedArray
from .signalhandler import SignalHandler
from .linearsolver import LinearSolver
from .errorhandler import assert_msg_critical
from .checkpoint import check_rsp_hdf5


class C6Solver(LinearSolver):
    """
    Implements the C6 value response solver.

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
    """

    def __init__(self, comm, ostream):
        """
        Initializes C6 solver to default setup.
        """

        super().__init__(comm, ostream)

        self.a_operator = 'dipole'
        self.a_components = 'xyz'
        self.b_operator = 'dipole'
        self.b_components = 'xyz'

        self.n_points = 9
        self.w0 = 0.3

        self.conv_thresh = 1.0e-3
        self.lindep_thresh = 1.0e-10

    def update_settings(self, rsp_dict, method_dict={}):
        """
        Updates response and method settings in C6 solver.

        :param rsp_dict:
            The dictionary of response dict.
        :param method_dict:
            The dictionary of method rsp_dict.
        """

        super().update_settings(rsp_dict, method_dict)

        if 'a_operator' in rsp_dict:
            self.a_operator = rsp_dict['a_operator'].lower()
        if 'a_components' in rsp_dict:
            self.a_components = rsp_dict['a_components'].lower()
        if 'b_operator' in rsp_dict:
            self.b_operator = rsp_dict['b_operator'].lower()
        if 'b_components' in rsp_dict:
            self.b_components = rsp_dict['b_components'].lower()

        if 'n_points' in rsp_dict:
            self.n_points = int(rsp_dict['n_points'])
        if 'w0' in rsp_dict:
            self.w0 = float(rsp_dict['w0'])

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

        ediag, sdiag = self.construct_ed_sd_half(orb_ene, nocc, norb)
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

    def precond_trials(self, vectors, precond=None):
        """
        Applies preconditioner to trial vectors.

        :param vectors:
            The set of vectors.
        :param precond:
            The preconditioner.

        :return:
            The preconditioned gerade and ungerade trial vectors.
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

        # Note: return the gerade and ungerade parts in order.
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
            self.print_header('C6 Value Response Solver',
                              n_points=self.n_points)

        self.start_time = tm.time()

        # sanity check
        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(nalpha == nbeta,
                            'C6Solver: not implemented for unrestricted case')

        if self.rank == mpi_master():
            nocc = molecule.number_of_alpha_electrons()
            norb = scf_tensors['C'].shape[1]
            orb_ene = scf_tensors['E']

        # ERI information
        eri_dict = self.init_eri(molecule, basis)

        # DFT information
        dft_dict = self.init_dft(molecule, scf_tensors)

        # PE information
        pe_dict = self.init_pe(molecule, basis)

        timing_dict = {}

        # right-hand side (gradient)
        b_rhs = self.get_complex_rhs(self.b_operator, self.b_components,
                                     molecule, basis, scf_tensors)
        if self.rank == mpi_master():
            imagfreqs = [
                self.w0 * (1 - t) / (1 + t)
                for t in np.polynomial.legendre.leggauss(self.n_points)
            ][0]
            imagfreqs = np.append(imagfreqs, 0.0)
            v1 = {(op, iw): v for op, v in zip(self.b_components, b_rhs)
                  for iw in imagfreqs}

        # operators, frequencies and preconditioners
        if self.rank == mpi_master():
            op_imagfreq_keys = list(v1.keys())
        else:
            op_imagfreq_keys = None
        op_imagfreq_keys = self.comm.bcast(op_imagfreq_keys, root=mpi_master())

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
                self.restart = check_rsp_hdf5(self.checkpoint_file,
                                              rsp_vector_labels, molecule,
                                              basis, dft_dict, pe_dict)
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

                assert_msg_critical(bger.any() or bung.any(),
                                    'C6Solver: trial vectors are empty')

                if bger is None or not bger.any():
                    bger = np.zeros((bung.shape[0], 0))
                if bung is None or not bung.any():
                    bung = np.zeros((bger.shape[0], 0))
            else:
                bger, bung = None, None

            e2bger, e2bung = self.e2n_half_size(bger, bung, molecule, basis,
                                                scf_tensors, eri_dict, dft_dict,
                                                pe_dict, timing_dict)

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
            self.graceful_exit, molecule, basis, dft_dict, pe_dict,
            [dist_bger, dist_bung, dist_e2bger, dist_e2bung], rsp_vector_labels)

        iter_per_trail_in_hours = None

        # start iterations
        for iteration in range(self.max_iter):

            iter_start_time = tm.time()

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

            profiler.stop_timer(iteration, 'Orthonorm.')

            if self.rank == mpi_master():

                assert_msg_critical(
                    new_trials_ger.any() or new_trials_ung.any(),
                    'C6Solver: unable to add new trial vectors')

                if new_trials_ger is None or not new_trials_ger.any():
                    new_trials_ger = np.zeros((new_trials_ung.shape[0], 0))
                if new_trials_ung is None or not new_trials_ung.any():
                    new_trials_ung = np.zeros((new_trials_ger.shape[0], 0))

                n_new_trials = new_trials_ger.shape[1] + new_trials_ung.shape[1]
            else:
                n_new_trials = None
            n_new_trials = self.comm.bcast(n_new_trials, root=mpi_master())

            if iter_per_trail_in_hours is not None:
                next_iter_in_hours = iter_per_trail_in_hours * n_new_trials
                if self.need_graceful_exit(next_iter_in_hours):
                    self.graceful_exit(
                        molecule, basis, dft_dict, pe_dict,
                        [dist_bger, dist_bung, dist_e2bger, dist_e2bung],
                        rsp_vector_labels)

            profiler.start_timer(iteration, 'FockBuild')

            # creating new sigma and rho linear transformations

            new_e2bger, new_e2bung = self.e2n_half_size(
                new_trials_ger, new_trials_ung, molecule, basis, scf_tensors,
                eri_dict, dft_dict, pe_dict, timing_dict)

            dist_bger.append(DistributedArray(new_trials_ger, self.comm),
                             axis=1)
            dist_bung.append(DistributedArray(new_trials_ung, self.comm),
                             axis=1)

            new_trials_ger, new_trials_ung = None, None

            dist_e2bger.append(DistributedArray(new_e2bger, self.comm), axis=1)
            dist_e2bung.append(DistributedArray(new_e2bung, self.comm), axis=1)

            new_e2bger, new_e2bung = None, None

            iter_in_hours = (tm.time() - iter_start_time) / 3600
            iter_per_trail_in_hours = iter_in_hours / n_new_trials

            profiler.stop_timer(iteration, 'FockBuild')
            if self.dft or self.pe:
                profiler.update_timer(iteration, timing_dict)

            profiler.check_memory_usage(
                'Iteration {:d} sigma build'.format(iteration + 1))

        signal_handler.remove_sigterm_function()

        self.write_checkpoint(molecule, basis, dft_dict, pe_dict,
                              [dist_bger, dist_bung, dist_e2bger, dist_e2bung],
                              rsp_vector_labels)

        # converged?
        if self.rank == mpi_master():
            self.print_convergence('Complex response')

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of C6 solver')
        profiler.print_memory_usage(self.ostream)

        # calculate response functions
        a_rhs = self.get_complex_rhs(self.a_operator, self.a_components,
                                     molecule, basis, scf_tensors)

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
