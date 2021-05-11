import numpy as np
import time as tm

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .veloxchemlib import GridDriver
from .veloxchemlib import XCFunctional
from .veloxchemlib import XCIntegrator
from .veloxchemlib import MolecularGrid
from .veloxchemlib import parse_xc_func
from .profiler import Profiler
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_scheme
from scipy.sparse import linalg


class OrbitalResponse:
    """
    Implements orbital response Lagrange multipliers computation using a
    conjugate gradient scheme for the time-dependent Hartree-Fock or DFT
    level of theory.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - is_tda: Flag if Tamm-Dancoff approximation is employed.
        - n_state_deriv: The number of the excited state of interest.
        - is_converged: The flag for convergence.
        - comm: The MPI communicator.
        - rank: The MPI rank.
        - dft: The flag for running DFT.
        - grid_level: The accuracy level of DFT grid.
        - xcfun: The XC functional.
        - eri_thresh: The electron repulsion integrals screening threshold.
        - qq_type: The electron repulsion integrals screening scheme.
        - nodes: Number of MPI processes.
        - ostream: The output stream.
        - conv_thresh: The convergence threshold for the solver.
        - max_iter: The maximum number of solver iterations.
        - iter_count: Index of the current iteration.
        - timing: The flag for printing timing information.
        - start_time: The start time of the calculation.
        - profiling: The flag for printing profiling information.
        - memory_profiling: The flag for printing memory usage.
        - memory_tracing: The flag for tracing memory allocation.
    """

    def __init__(self, comm, ostream):
        """
        Initializes orbital response computation driver to default setup.
        """

        # ERI settings
        self.eri_thresh = 1.0e-15
        self.qq_type = 'QQ_DEN'

        # Solver setup
        self.conv_thresh = 1.0e-4
        self.max_iter = 50
        self.iter_count = 0
        self.is_converged = False

        # MPI information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        # DFT information
        self.dft = False
        self.grid_level = 4
        self.xcfun = XCFunctional()

        # Output stream
        self.ostream = ostream

        # Flag on whether RPA or TDA is calculated
        self.is_tda = False

        # Excited state information, default to first excited state
        self.n_state_deriv = 0

        # Timing and profiling
        self.timing = False
        self.profiling = False
        self.memory_profiling = False
        self.memory_tracing = False

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in orbital response computation
        driver.

        :param rsp_dict:
            The dictionary of response settings.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}

        # ERI settings
        if 'eri_thresh' in rsp_dict:
            self.eri_thresh = float(rsp_dict['eri_thresh'])
        if 'qq_type' in rsp_dict:
            self.qq_type = rsp_dict['qq_type']

        # Solver setup
        # TODO: use specific orbital response keywords here?
        if 'conv_thresh' in rsp_dict:
            self.conv_thresh = float(rsp_dict['conv_thresh'])
        if 'max_iter' in rsp_dict:
            self.max_iter = int(rsp_dict['max_iter'])

        # Use TDA or not
        if 'tamm_dancoff' in rsp_dict:
            key = rsp_dict['tamm_dancoff'].lower()
            self.is_tda = True if key in ['yes', 'y'] else False

        # Excited state of interest
        if 'n_state_deriv' in rsp_dict:
            # user gives '1' for first excited state, but internal index is 0
            self.n_state_deriv = int(rsp_dict['n_state_deriv']) - 1

        # DFT
        if 'dft' in method_dict:
            key = method_dict['dft'].lower()
            self.dft = True if key in ['yes', 'y'] else False
        if 'grid_level' in method_dict:
            self.grid_level = int(method_dict['grid_level'])
        if 'xcfun' in method_dict:
            if 'dft' not in method_dict:
                self.dft = True
            self.xcfun = parse_xc_func(method_dict['xcfun'].upper())
            assert_msg_critical(not self.xcfun.is_undefined(),
                                'Response solver: Undefined XC functional')

        # Timing and profiling
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

    def init_dft(self, molecule, scf_tensors):
        """
        Initializes DFT.

        :param molecule:
            The molecule.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            The dictionary of DFT information.
        """

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
                gs_density = AODensityMatrix([scf_tensors['D_alpha']],
                                             denmat.rest)
            else:
                gs_density = AODensityMatrix()
            gs_density.broadcast(self.rank, self.comm)

            dft_func_label = self.xcfun.get_func_label().upper()
        else:
            molgrid = MolecularGrid()
            gs_density = AODensityMatrix()
            dft_func_label = 'HF'

        return {
            'molgrid': molgrid,
            'gs_density': gs_density,
            'dft_func_label': dft_func_label,
        }

    def compute(self, molecule, basis, scf_tensors, rsp_results):
        """
        Computes orbital response Lagrange multipliers and relaxed density
        for the calculation of energy gradients.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from the converged SCF wavefunction.
        :param rsp_results:
            The results from the RPA or TDA excited states calculation.
        :param profiler:
            The profiler.
        """

        profiler = Profiler({
            'timing': self.timing,
            'profiling': self.profiling,
            'memory_profiling': self.memory_profiling,
            'memory_tracing': self.memory_tracing,
        })

        if self.rank == mpi_master():
            self.print_orbrsp_header('Orbital Response Driver',
                                     self.n_state_deriv)

        # DFT information

        dft_dict = self.init_dft(molecule, scf_tensors)

        # set start time

        self.start_time = tm.time()

        # sanity check

        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'OrbitalResponse: not implemented for unrestricted case')

        # Workflow:
        # 1) Construct the necessary density matrices => in child classes
        # 2) Construct the RHS => in child classes
        # 3) Construct the initial guess
        # 4) Write the linear operator for matrix-vector product
        # 5) Run the conjugate gradient

        rhs_results = self.compute_rhs(molecule, basis, scf_tensors,
                                       rsp_results, dft_dict, profiler)

        if self.rank == mpi_master():
            rhs_mo = rhs_results['rhs_mo']
            dm_oo = rhs_results['dm_oo']
            dm_vv = rhs_results['dm_vv']
            unrel_dm_ao = rhs_results['unrel_dm_ao']
            fock_ao_rhs = rhs_results['fock_ao_rhs']
        else:
            rhs_mo = None

        rhs_mo = self.comm.bcast(rhs_mo, root=mpi_master())

        # Calculate the lambda multipliers in the parent class
        lambda_multipliers = self.compute_lambda(molecule, basis, scf_tensors,
                                                 rhs_mo, dft_dict, profiler)

        profiler.start_timer(0, 'omega')

        # Prerequesites for the overlap matrix multipliers
        if self.rank == mpi_master():

            # 1) Compute an energy-weighted density matrix
            nocc = molecule.number_of_alpha_electrons()
            mo = scf_tensors['C']
            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            mo_energies = scf_tensors['E']
            eocc = mo_energies[:nocc]
            evir = mo_energies[nocc:]
            eo_diag = np.diag(eocc)
            ev_diag = np.diag(evir)

            epsilon_dm_ao = np.linalg.multi_dot(
                [mo_occ,
                 np.matmul(eo_diag, 0.5 * dm_oo) + eo_diag, mo_occ.T])
            epsilon_dm_ao += np.linalg.multi_dot(
                [mo_vir, np.matmul(ev_diag, 0.5 * dm_vv), mo_vir.T])
            epsilon_dm_ao += np.linalg.multi_dot(
                [mo_occ,
                 np.matmul(eo_diag, lambda_multipliers), mo_vir.T])
            epsilon_dm_ao += np.linalg.multi_dot(
                [mo_occ,
                 np.matmul(eo_diag, lambda_multipliers), mo_vir.T]).T

            # 2) Transform the lambda multipliers to AO basis:
            lambda_ao = np.linalg.multi_dot(
                [mo_occ, lambda_multipliers, mo_vir.T])
            ao_density_lambda = AODensityMatrix([lambda_ao], denmat.rest)
        else:
            ao_density_lambda = AODensityMatrix()
        ao_density_lambda.broadcast(self.rank, self.comm)

        molgrid = dft_dict['molgrid']
        gs_density = dft_dict['gs_density']

        fock_lambda = AOFockMatrix(ao_density_lambda)
        fock_flag = fockmat.rgenjk
        if self.dft:
            if self.xcfun.is_hybrid():
                fock_flag = fockmat.rgenjkx
                fact_xc = self.xcfun.get_frac_exact_exchange()
                fock_lambda.set_scale_factor(fact_xc, 0)
            else:
                fock_flag = fockmat.rgenj

        fock_lambda.set_fock_type(fock_flag, 0)

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)
        eri_drv.compute(fock_lambda, ao_density_lambda, molecule, basis,
                        screening)
        if self.dft:
            if not self.xcfun.is_hybrid():
                fock_lambda.scale(2.0, 0)
            xc_drv = XCIntegrator(self.comm)
            molgrid.distribute(self.rank, self.nodes, self.comm)
            xc_drv.integrate(fock_lambda, ao_density_lambda, gs_density,
                             molecule, basis, molgrid, self.xcfun.get_func_label())

        fock_lambda.reduce_sum(self.rank, self.nodes, self.comm)

        # Compute the omega multipliers
        if self.rank == mpi_master():
            ovlp = scf_tensors['S']
            omega_ao = self.compute_omega(ovlp, mo_occ, mo_vir, epsilon_dm_ao,
                                          rsp_results, fock_ao_rhs, fock_lambda)

        self.ostream.print_blank()
        self.ostream.flush()

        profiler.stop_timer(0, 'omega')

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of CG')
        profiler.print_memory_usage(self.ostream)

        if self.rank == mpi_master():
            # print warning if Lambda did not converge
            if not self.is_converged:
                warn_msg = '*** Warning: Orbital response did not converge.'
                self.ostream.print_header(warn_msg.ljust(56))
            else:
                orbrsp_time = tm.time() - self.start_time
                self.ostream.print_info(
                    'Orbital response converged after {:d} iterations.'.format(
                        self.iter_count))
                self.ostream.print_info(
                    'Total time needed for orbital response: {:.2f} s.'.format(
                        orbrsp_time))

            self.ostream.flush()

        if self.rank == mpi_master():
            if self.is_converged:
                # Calculate the relaxed one-particle density matrix
                # Factor 4: (ov + vo)*(alpha + beta)
                rel_dm_ao = unrel_dm_ao + 4 * lambda_ao
                return {
                    'lambda_ao': lambda_ao,
                    'omega_ao': omega_ao,
                    'unrel_dm_ao': unrel_dm_ao,
                    'rel_dm_ao': rel_dm_ao,
                }
            else:
                # return only unrelaxed density matrix if not converged
                return {'unrel_dm_ao': unrel_dm_ao}
        else:
            return {}

    def compute_lambda(self, molecule, basis, scf_tensors, rhs_mo, dft_dict, profiler):
        """
        Performs orbital response Lagrange multipliers calculation for the
        occupied-virtual alpha block using molecular data.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param rhs_mo:
            The right-hand side of orbital response equation in MO basis.
        :param dft_dict:
            The dictionary containing DFT information.
        :param profiler:
            The profiler.

        :return:
            A dictionary containing the Lagrange multipliers and relaxed
            one-particle density.
        """

        # count variable for conjugate gradient iterations
        self.iter_count = 0

        nocc = molecule.number_of_alpha_electrons()

        if self.rank == mpi_master():
            mo = scf_tensors['C']
            mo_energies = scf_tensors['E']

            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nvir = mo_vir.shape[1]

            eocc = mo_energies[:nocc]
            evir = mo_energies[nocc:]
            eov = eocc.reshape(-1, 1) - evir
        else:
            nvir = None
            eov = None
        nvir = self.comm.bcast(nvir, root=mpi_master())
        eov = self.comm.bcast(eov, root=mpi_master())

        # Calculate the initial guess for the Lagrange multipliers given by
        # the RHS divided by orbital-energy differences
        lambda_guess = rhs_mo / eov

        if self.rank == mpi_master():
            # Create AODensityMatrix object from lambda in AO
            lambda_ao = np.linalg.multi_dot([mo_occ, lambda_guess, mo_vir.T])
            ao_density_lambda = AODensityMatrix([lambda_ao], denmat.rest)
        else:
            ao_density_lambda = AODensityMatrix()
        ao_density_lambda.broadcast(self.rank, self.comm)

        # TODO: make sure that ERI settings are consistent with response solver

        # ERI driver
        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        molgrid = dft_dict['molgrid']
        gs_density = dft_dict['gs_density']

        # Create a Fock Matrix Object (initialized with zeros)
        fock_lambda = AOFockMatrix(ao_density_lambda)
        fock_flag = fockmat.rgenjk
        if self.dft:
            if self.xcfun.is_hybrid():
                fock_flag = fockmat.rgenjkx
                fact_xc = self.xcfun.get_frac_exact_exchange()
                fock_lambda.set_scale_factor(fact_xc, 0)
            else:
                fock_flag = fockmat.rgenj

        fock_lambda.set_fock_type(fock_flag, 0)

        # Matrix-vector product of orbital Hessian with trial vector
        def orb_rsp_matvec(v):
            """
            Function to carry out matrix multiplication of Lagrange multiplier
            vector with orbital Hessian matrix.
            """

            profiler.start_timer(self.iter_count, 'CG')

            # Create AODensityMatrix object from lambda in AO
            if self.rank == mpi_master():
                lambda_ao = np.linalg.multi_dot(
                    [mo_occ, v.reshape(nocc, nvir), mo_vir.T])
                ao_density_lambda = AODensityMatrix([lambda_ao], denmat.rest)
            else:
                ao_density_lambda = AODensityMatrix()
            ao_density_lambda.broadcast(self.rank, self.comm)

            eri_drv.compute(fock_lambda, ao_density_lambda, molecule, basis,
                            screening)
            if self.dft:
                #t0 = tm.time()
                if not self.xcfun.is_hybrid():
                    fock_lambda.scale(2.0, 0)
                xc_drv = XCIntegrator(self.comm)
                molgrid.distribute(self.rank, self.nodes, self.comm)
                xc_drv.integrate(fock_lambda, ao_density_lambda, gs_density,
                                 molecule, basis, molgrid, self.xcfun.get_func_label())
                #if timing_dict is not None:
                #    timing_dict['DFT'] = tm.time() - t0

            fock_lambda.reduce_sum(self.rank, self.nodes, self.comm)

            # Transform to MO basis (symmetrized w.r.t. occ. and virt.)
            # and add diagonal part
            if self.rank == mpi_master():
                fock_lambda_0 = fock_lambda.alpha_to_numpy(0)
                lambda_mo = (
                    -(np.linalg.multi_dot([mo_occ.T, fock_lambda_0, mo_vir]) +
                      np.linalg.multi_dot([mo_vir.T, fock_lambda_0, mo_occ]).T)
                    + v.reshape(nocc, nvir) * eov)
            else:
                lambda_mo = None

            lambda_mo = self.comm.bcast(lambda_mo, root=mpi_master())

            profiler.stop_timer(self.iter_count, 'CG')

            profiler.check_memory_usage(
                'CG Iteration {:d}'.format(self.iter_count + 1))

            profiler.print_memory_tracing(self.ostream)

            # increase iteration counter every time this function is called
            self.iter_count += 1

            return lambda_mo.reshape(nocc * nvir)

        # Matrix-vector product for preconditioner using the
        # inverse of the diagonal (i.e. eocc - evir)
        def precond_matvec(v):
            """
            Function that defines the matrix-vector product
            required by the pre-conditioner for the conjugate gradient.
            It is an approximation for the inverse of matrix A in Ax = b.
            """
            current_v = v.reshape(nocc, nvir)
            M_dot_v = current_v / eov

            return M_dot_v.reshape(nocc * nvir)

        # 5) Define the linear operators and run conjugate gradient
        LinOp = linalg.LinearOperator((nocc * nvir, nocc * nvir),
                                      matvec=orb_rsp_matvec)
        PrecondOp = linalg.LinearOperator((nocc * nvir, nocc * nvir),
                                          matvec=precond_matvec)

        b = rhs_mo.reshape(nocc * nvir)
        x0 = lambda_guess.reshape(nocc * nvir)

        lambda_multipliers, cg_conv = linalg.cg(A=LinOp,
                                                b=b,
                                                x0=x0,
                                                M=PrecondOp,
                                                tol=self.conv_thresh,
                                                atol=0,
                                                maxiter=self.max_iter)

        self.is_converged = (cg_conv == 0)

        return lambda_multipliers.reshape(nocc, nvir)

    def print_orbrsp_header(self, title, n_state_deriv):
        self.ostream.print_blank()
        self.ostream.print_header('{:s} Setup'.format(title))
        self.ostream.print_header('=' * (len(title) + 8))
        self.ostream.print_blank()

        str_width = 60

        # print solver-specific info

        if n_state_deriv is not None:
            cur_str = 'Excited State of Interest       : ' + str(n_state_deriv +
                                                                 1)
            self.ostream.print_header(cur_str.ljust(str_width))

        # print general info

        cur_str = 'Max. Number of Iterations       : ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold           : {:.1e}'.format(
            self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()
