import numpy as np
import time as tm

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .linearsolver import LinearSolver
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_scheme
from scipy.sparse import linalg


class OrbitalResponse(LinearSolver):
    """
    Implements orbital response Lagrange multipliers computation using a
    conjugate gradient scheme for the time-dependent Hartree-Fock or DFT
    level of theory.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - n_state_deriv: The number of the excited state of interest.
        - solver: The linear equations solver.
    """

    def __init__(self, comm, ostream):
        """
        Initializes orbital response computation driver to default setup.
        """

        super().__init__(comm, ostream)

        # excited state information, default to first excited state
        self.n_state_deriv = 0

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

        # Many settings updated in LinearSolver
        super().update_settings(rsp_dict, method_dict)

        if 'n_state_deriv' in rsp_dict:
            # user gives '1' for first excited state, but internal index is 0
            self.n_state_deriv = int(rsp_dict['n_state_deriv']) - 1

    def compute_lambda(self, molecule, basis, scf_tensors, rhs_mo, profiler):
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
        :param profiler:
            The profiler.

        :return:
            A dictionary containing the Lagrange multipliers and relaxed
            one-particle density.
        """

        if self.rank == mpi_master():
            self.print_orbrsp_header('Orbital Response Driver',
                                     self.n_state_deriv)

        # set start time

        self.start_time = tm.time()

        # sanity check

        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'OrbitalResponse: not implemented for unrestricted case')

        # count variable for conjugate gradient iterations
        self.iter_count = 0

        # Workflow:
        # 1) Construct the necessary density matrices => in child classes
        # 2) Construct the RH => in child classes
        # 3) Construct the initial gues
        # 4) Write the linear operator for matrix-vector product
        # 5) Run the conjugate gradient

        nocc = molecule.number_of_alpha_electrons()

        if self.rank == mpi_master():
            mo = scf_tensors['C']
            ea = scf_tensors['E']

            mo_occ = mo[:, :nocc].copy()
            mo_vir = mo[:, nocc:].copy()
            nvir = mo_vir.shape[1]

            eocc = ea[:nocc]
            evir = ea[nocc:]
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

        # Create a Fock Matrix Object (initialized with zeros)
        fock_lambda = AOFockMatrix(ao_density_lambda)
        fock_lambda.set_fock_type(fockmat.rgenjk, 0)

        # matrix-vector product
        def orb_rsp_matvec(v):
            """
            Function to carry out matrix multiplication of Lagrange multipier
            vector with orbital Hessian matrix
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

        # 5) Define the linear operator and run conjugate gradient
        LinOp = linalg.LinearOperator((nocc * nvir, nocc * nvir),
                                      matvec=orb_rsp_matvec)

        b = rhs_mo.reshape(nocc * nvir)
        x0 = lambda_guess.reshape(nocc * nvir)

        lambda_multipliers, cg_conv = linalg.cg(A=LinOp,
                                                b=b,
                                                x0=x0,
                                                tol=self.conv_thresh,
                                                atol=0,
                                                maxiter=self.max_iter)

        self.is_converged = (cg_conv == 0)

        # TODO: print warning or something if not converged

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
