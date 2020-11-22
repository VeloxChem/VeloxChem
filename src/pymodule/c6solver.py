import numpy as np
import time as tm
import psutil
import math
import sys

from .veloxchemlib import mpi_master
from .profiler import Profiler
from .distributedarray import DistributedArray
from .signalhandler import SignalHandler
from .linearsolver import LinearSolver
from .errorhandler import assert_msg_critical
from .checkpoint import check_rsp_hdf5
from .checkpoint import append_rsp_solution_hdf5


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

        self.a_operator = 'electric dipole'
        self.a_components = 'xyz'
        self.b_operator = 'electric dipole'
        self.b_components = 'xyz'

        self.n_points = 9
        self.w0 = 0.3

        self.conv_thresh = 1.0e-3
        self.lindep_thresh = 1.0e-10

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in C6 solver.

        :param rsp_dict:
            The dictionary of response dict.
        :param method_dict:
            The dictionary of method rsp_dict.
        """

        if method_dict is None:
            method_dict = {}

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

    def get_precond(self, orb_ene, nocc, norb, iw):
        """
        Constructs the preconditioners.

        :param orb_ene:
            The orbital energies.
        :param nocc:
            The number of doubly occupied orbitals.
        :param norb:
            The number of orbitals.
        :param iw:
            The imaginary frequency.

        :return:
            The distributed preconditioners.
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

        return tuple(
            [DistributedArray(p, self.comm) for p in [pa_diag, pc_diag]])

    def preconditioning(self, precond, v_in):
        """
        Applies preconditioner to a tuple of distributed trial vectors.

        :param precond:
            The preconditioner.
        :param v_in:
            The input trial vectors.

        :return:
            A tuple of distributed trail vectors after preconditioning.
        """

        pa, pc = [p.data for p in precond]

        v_in_ru, v_in_ig = [v.data for v in v_in]

        v_out_ru = pa * v_in_ru + pc * v_in_ig
        v_out_ig = pc * v_in_ru - pa * v_in_ig

        return tuple([
            DistributedArray(v, self.comm, distribute=False)
            for v in [v_out_ru, v_out_ig]
        ])

    def precond_trials(self, vectors, precond):
        """
        Applies preconditioner to distributed trial vectors.

        :param vectors:
            The set of vectors.
        :param precond:
            The preconditioner.

        :return:
            The preconditioned gerade and ungerade trial vectors.
        """

        trials_ger = []
        trials_ung = []
        sqrt_2 = math.sqrt(2.0)

        for (op, w), vec in vectors.items():
            v = self.preconditioning(precond[w], vec)
            norms = [sqrt_2 * p.norm() for p in v]
            vn = math.sqrt(sum([n**2 for n in norms]))

            if vn > self.small_thresh:
                # real ungerade
                if norms[0] > self.small_thresh:
                    trials_ung.append(v[0].data)
                # imaginary gerade
                if norms[1] > self.small_thresh:
                    trials_ger.append(v[1].data)

        new_ger = np.array(trials_ger).T
        new_ung = np.array(trials_ung).T

        dist_new_ger = DistributedArray(new_ger, self.comm, distribute=False)
        dist_new_ung = DistributedArray(new_ung, self.comm, distribute=False)

        return dist_new_ger, dist_new_ung

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

        self.dist_bger = None
        self.dist_bung = None
        self.dist_e2bger = None
        self.dist_e2bung = None

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
            orb_ene = scf_tensors['E']
        else:
            orb_ene = None
        orb_ene = self.comm.bcast(orb_ene, root=mpi_master())
        norb = orb_ene.shape[0]
        nocc = molecule.number_of_alpha_electrons()

        # ERI information
        eri_dict = self.init_eri(molecule, basis)

        # DFT information
        dft_dict = self.init_dft(molecule, scf_tensors)

        # PE information
        pe_dict = self.init_pe(molecule, basis)

        timing_dict = {}

        # right-hand side (gradient)
        b_rhs = self.get_complex_prop_grad(self.b_operator, self.b_components,
                                           molecule, basis, scf_tensors)

        imagfreqs = [
            self.w0 * (1 - t) / (1 + t)
            for t in np.polynomial.legendre.leggauss(self.n_points)
        ][0]
        imagfreqs = np.append(imagfreqs, 0.0)

        v1 = None
        if self.rank == mpi_master():
            v1 = {(op, iw): v for op, v in zip(self.b_components, b_rhs)
                  for iw in imagfreqs}

        # operators, frequencies and preconditioners
        if self.rank == mpi_master():
            op_imagfreq_keys = list(v1.keys())
        else:
            op_imagfreq_keys = None
        op_imagfreq_keys = self.comm.bcast(op_imagfreq_keys, root=mpi_master())

        precond = {
            iw: self.get_precond(orb_ene, nocc, norb, iw) for iw in imagfreqs
        }

        # distribute the right-hand side
        # dist_v1 will also serve as initial guess
        dist_v1 = {}
        for key in op_imagfreq_keys:
            if self.rank == mpi_master():
                gradger, gradung = self.decomp_grad(v1[key])
                grad_ru = gradung.real
                grad_ig = gradger.imag
            else:
                grad_ru, grad_ig = None, None
            dist_v1[key] = (DistributedArray(grad_ru, self.comm),
                            DistributedArray(grad_ig, self.comm))

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
            (self.dist_bger, self.dist_bung, self.dist_e2bger,
             self.dist_e2bung) = self.read_vectors(rsp_vector_labels)

        # generate initial guess from scratch
        else:
            bger, bung = self.setup_trials(dist_v1, precond)

            self.e2n_half_size(bger, bung, molecule, basis, scf_tensors,
                               eri_dict, dft_dict, pe_dict, timing_dict)

        profiler.check_memory_usage('Initial guess')

        solutions = {}
        residuals = {}
        relative_residual_norm = {}

        signal_handler = SignalHandler()
        signal_handler.add_sigterm_function(self.graceful_exit, molecule, basis,
                                            dft_dict, pe_dict,
                                            rsp_vector_labels)

        iter_per_trail_in_hours = None
        sqrt_2 = math.sqrt(2.0)

        # start iterations
        for iteration in range(self.max_iter):

            iter_start_time = tm.time()

            profiler.start_timer(iteration, 'ReducedSpace')

            xvs = []
            self.cur_iter = iteration

            n_ger = self.dist_bger.shape(1)
            n_ung = self.dist_bung.shape(1)

            e2gg = self.dist_bger.matmul_AtB(self.dist_e2bger, 2.0)
            e2uu = self.dist_bung.matmul_AtB(self.dist_e2bung, 2.0)
            s2ug = self.dist_bung.matmul_AtB(self.dist_bger, 4.0)

            for op, iw in op_imagfreq_keys:
                if (iteration == 0 or
                        relative_residual_norm[(op, iw)] > self.conv_thresh):

                    grad_ru, grad_ig = dist_v1[(op, iw)]

                    # projections onto gerade and ungerade subspaces:

                    g_realung = self.dist_bung.matmul_AtB(grad_ru, 2.0)

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

                    x_realung = self.dist_bung.matmul_AB_no_gather(c_realung)
                    x_imagger = self.dist_bger.matmul_AB_no_gather(c_imagger)

                    # composing E2 matrices projected onto solution subspace

                    e2imagger = self.dist_e2bger.matmul_AB_no_gather(c_imagger)
                    e2realung = self.dist_e2bung.matmul_AB_no_gather(c_realung)

                    # calculating the residual components

                    s2imagger = 2.0 * x_imagger.data
                    s2realung = 2.0 * x_realung.data

                    r_realung = (e2realung.data + iw * s2imagger - grad_ru.data)
                    r_imagger = (-e2imagger.data + iw * s2realung +
                                 grad_ig.data)

                    r_realung = DistributedArray(r_realung,
                                                 self.comm,
                                                 distribute=False)
                    r_imagger = DistributedArray(r_imagger,
                                                 self.comm,
                                                 distribute=False)

                    # calculating relative residual norm
                    # for convergence check

                    r = (r_realung, r_imagger)
                    x = (x_realung, x_imagger)

                    x_full = self.get_full_solution_vector(x)
                    if self.rank == mpi_master():
                        xv = np.dot(x_full, v1[(op, iw)])
                        xvs.append((op, iw, xv))

                    r_norms = [sqrt_2 * p.norm() for p in r]
                    x_norms = [sqrt_2 * p.norm() for p in x]

                    rn = math.sqrt(sum([n**2 for n in r_norms]))
                    xn = math.sqrt(sum([n**2 for n in x_norms]))

                    if xn != 0:
                        relative_residual_norm[(op, iw)] = rn / xn
                    else:
                        relative_residual_norm[(op, iw)] = rn

                    if relative_residual_norm[(op, iw)] < self.conv_thresh:
                        solutions[(op, iw)] = x
                    else:
                        residuals[(op, iw)] = r

            # write to output
            if self.rank == mpi_master():

                self.ostream.print_info(
                    '{:d} gerade trial vectors in reduced space'.format(n_ger))
                self.ostream.print_info(
                    '{:d} ungerade trial vectors in reduced space'.format(
                        n_ung))
                self.ostream.print_blank()

                profiler.print_memory_subspace(
                    {
                        'dist_bger': self.dist_bger,
                        'dist_bung': self.dist_bung,
                        'dist_e2bger': self.dist_e2bger,
                        'dist_e2bung': self.dist_e2bung,
                        'precond': precond,
                        'solutions': solutions,
                        'residuals': residuals,
                    }, self.ostream)

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
                residuals, precond, self.dist_bger, self.dist_bung)

            residuals.clear()

            profiler.stop_timer(iteration, 'Orthonorm.')

            if self.rank == mpi_master():
                n_new_trials = new_trials_ger.shape(1) + new_trials_ung.shape(1)
            else:
                n_new_trials = None
            n_new_trials = self.comm.bcast(n_new_trials, root=mpi_master())

            if iter_per_trail_in_hours is not None:
                next_iter_in_hours = iter_per_trail_in_hours * n_new_trials
                if self.need_graceful_exit(next_iter_in_hours):
                    self.graceful_exit(molecule, basis, dft_dict, pe_dict,
                                       rsp_vector_labels)

            profiler.start_timer(iteration, 'FockBuild')

            # creating new sigma and rho linear transformations

            self.e2n_half_size(new_trials_ger, new_trials_ung, molecule, basis,
                               scf_tensors, eri_dict, dft_dict, pe_dict,
                               timing_dict)

            iter_in_hours = (tm.time() - iter_start_time) / 3600
            iter_per_trail_in_hours = iter_in_hours / n_new_trials

            profiler.stop_timer(iteration, 'FockBuild')
            if self.dft or self.pe:
                profiler.update_timer(iteration, timing_dict)

            profiler.check_memory_usage(
                'Iteration {:d} sigma build'.format(iteration + 1))

        signal_handler.remove_sigterm_function()

        dist_arrays = [
            self.dist_bger,
            self.dist_bung,
            self.dist_e2bger,
            self.dist_e2bung,
        ]

        self.write_checkpoint(molecule, basis, dft_dict, pe_dict, dist_arrays,
                              rsp_vector_labels)

        # converged?
        if self.rank == mpi_master():
            self.print_convergence('Complex response')

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of C6 solver')
        profiler.print_memory_usage(self.ostream)

        # calculate response functions
        a_rhs = self.get_complex_prop_grad(self.a_operator, self.a_components,
                                           molecule, basis, scf_tensors)

        if self.is_converged:
            key_0 = list(solutions.keys())[0]
            x_0 = self.get_full_solution_vector(solutions[key_0])

            if self.rank == mpi_master():
                avail_mem = psutil.virtual_memory().available
                write_solution_to_file = (
                    avail_mem < sys.getsizeof(x_0) * len(solutions) * 2)
                full_solutions = {}

                va = {op: v for op, v in zip(self.a_components, a_rhs)}
                rsp_funcs = {}

            for bop, iw in solutions:
                x = self.get_full_solution_vector(solutions[(bop, iw)])

                if self.rank == mpi_master():
                    for aop in self.a_components:
                        rsp_funcs[(aop, bop, iw)] = -np.dot(va[aop], x)

                        if write_solution_to_file:
                            append_rsp_solution_hdf5(
                                self.checkpoint_file,
                                '{:s}_{:s}_{:.8f}'.format(aop, bop, iw), x)
                        else:
                            full_solutions[(bop, iw)] = x

            if self.rank == mpi_master():
                if write_solution_to_file:
                    return {'response_functions': rsp_funcs}
                else:
                    return {
                        'response_functions': rsp_funcs,
                        'solutions': full_solutions
                    }

        return {}

    def get_full_solution_vector(self, solution):
        """
        Gets a full solution vector from the distributed solution.

        :param solution:
            The distributed solution as a tuple.

        :return:
            The full solution vector.
        """

        (dist_x_realung, dist_x_imagger) = solution

        x_realung = dist_x_realung.get_full_vector()
        x_imagger = dist_x_imagger.get_full_vector()

        if self.rank == mpi_master():
            x_real = np.hstack((x_realung, -x_realung))
            x_imag = np.hstack((x_imagger, x_imagger))
            return x_real + 1j * x_imag
        else:
            return None

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
