import numpy as np
import time as tm

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .profiler import Profiler
from .linearsolver import LinearSolver
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_scheme
from .checkpoint import read_rsp_hdf5
from .checkpoint import write_rsp_hdf5
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
        - unrel_dm_ao: The unrelaxed density matrix in AO basis.
        - rel_dm_ao: The relaxed density matrix in AO basis.
        - rhs_mo: The right-hand side of orbital response equation in MO basis.
        - lambda_ao: The converged lambda multipliers in AO basis.
        - omega_ao: The omega multipliers in AO basis.
    """

    # Give number of state of interest and all vectors in compute function?
    # Or only the specific vector in compute function?
    def __init__(self, comm, ostream):
        """
        Initializes orbital response computation driver to default setup.
        """

        # many settings from LinearSolver (ERI, DFT, MPI, ...)
        super().__init__(comm, ostream)

        # excited state information, default to first excited state
        self.n_state_deriv = 0

        # solver setup
        self.solver = None

        # More instance variables
        self.unrel_dm_ao = None
        self.rel_dm_ao = None
        self.rhs_mo = None
        self.lambda_ao = None
        self.omega_ao = None

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

    def compute_lambda(self, molecule, basis, scf_tensors):
        """
        Performs orbital response Lagrange multipliers calculation using molecular data.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.
        :param excitation_vecs:
            The excitation vectors from converged RPA or TDA calculation.

        :return:
            A dictionary containing the Lagrange multipliers and relaxed
            one-particle density.
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

        # set start time

        self.start_time = tm.time()

        # sanity check

        nalpha = molecule.number_of_alpha_electrons()
        nbeta = molecule.number_of_beta_electrons()
        assert_msg_critical(
            nalpha == nbeta,
            'OrbitalResponse: not implemented for unrestricted case')

        # ERI information
        eri_dict = self.init_eri(molecule, basis)

        # DFT information
        dft_dict = self.init_dft(molecule, scf_tensors)

        # PE information
        pe_dict = self.init_pe(molecule, basis)

        timing_dict = {}

        # block Davidson algorithm setup

        #self.solver = BlockDavidsonSolver() # something with linalg.cg?

        # read initial guess from restart file

        # It might be possible to restart by using the current lambda
        # and the right-hand side
        # Restart iterations number
        #n_restart_iterations = 4

        # count variable for conjugate gradient iterations
        self.iter_count = 0

		# => Reading from checkpoint files does not improve convergence
		#    when using scipy's conjugate gradient, it would probably
		#    need information from two previous iterations...
        #if self.restart:
        #    print("Restarting orbital response in compute_lambda:")
        #    print("checkpoint file =", self.checkpoint_file)
        #    if self.rank == mpi_master():
        #        #rhs_mo_from_chk = read_rsp_hdf5(
        #        #    self.checkpoint_file, ['OrbRsp_RHS'],
        #        #    molecule, basis, dft_dict, pe_dict, self.ostream)
        #        self.rhs_mo, rst_orbrsp_lambda = read_rsp_hdf5(
        #            self.checkpoint_file, ['OrbRsp_RHS', 'OrbRsp_lambda'],
        #            molecule, basis, dft_dict, pe_dict, self.ostream)

        #        print("This is what I'm reading from checkpoint:\n")
        #        print("RHS:\n", self.rhs_mo)
        #        print("\nLambda\n", rst_orbrsp_lambda)
        #        # What does this do??
        #        #self.restart = (self.rhs_mo is not None and
        #        #                rst_orbrsp_lambda is not None)
        #else:
        #    # Write RHS to checkpoint file
        #    if self.rank == mpi_master():
        #        print("Writing checkpoint file, RHS")
        #        write_rsp_hdf5(self.checkpoint_file, [self.rhs_mo],
        #                       ['OrbRsp_RHS'], molecule, basis,
        #                       dft_dict, pe_dict, self.ostream)
        #        rhs_mo_from_chk = read_rsp_hdf5(
        #            self.checkpoint_file, ['OrbRsp_RHS'],
        #            molecule, basis, dft_dict, pe_dict, self.ostream)
        #        print(rhs_mo_from_chk)


        # TODO
        # 1) Construct the necessary density matrices => in child classes
        # 2) Construct the RH => in child classesS
        # 3) Construct the initial gues
        # 4) Write the linear operator for matrix-vector product
        # 5) Run the conjugate gradient

        # 3) Calculate the initial guess for the Lagrange multipliers
        #   given by the RHS divided by orbital-energy differences
        nocc = molecule.number_of_alpha_electrons()
        mo_occ = scf_tensors['C'][:, :nocc]  # occupied MO coefficients
        mo_vir = scf_tensors['C'][:, nocc:]  # virtual MO coefficients
        nvir = mo_vir.shape[1]
        eocc = scf_tensors['E'][:nocc]
        evir = scf_tensors['E'][nocc:]
        eov = eocc.reshape(-1, 1) - evir

        #if self.restart:
        #    lambda_guess = rst_orbrsp_lambda
        #else:
        lambda_guess = self.rhs_mo / eov

        # Transform to AO
        lambda_ao = np.matmul(mo_occ,
                              np.matmul(lambda_guess, mo_vir.T))

        # Create AODensityMatrix object from lambda in AO
        ao_density_lambda = AODensityMatrix([lambda_ao], denmat.rest)

        # ERI driver
        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        # Create a Fock Matrix Object (initialized with zeros)
        fock_flag = fockmat.rgenjk
        fock_lambda = AOFockMatrix(ao_density_lambda)
        fock_lambda.set_fock_type(fock_flag, 0)

        # 4) Write a function for matrix-vector product
        #   (inside compute since it needs to know many variables from outside)
        def OrbRsp_MatVec(v):
             """
             Function to carry out matrix multiplication
             of Lagrange multipier vector with orbital Hessian
             matrix, using AODensityMatrix, AOFockMatrix

             v: lambda at current iteration
             """
             # Transform to AO
             lambda_ao = np.matmul(mo_occ,
                                   np.matmul(v.reshape(nocc, nvir), mo_vir.T))

             # Create AODensityMatrix object from lambda in AO
             ao_density_lambda = AODensityMatrix([lambda_ao], denmat.rest)

             eri_drv.compute(fock_lambda, ao_density_lambda, molecule, basis,
                             screening)

             # Transform to MO basis (symmetrized w.r.t. occ. and virt.) and add diagonal part
             lambda_mo = (
                 -(np.matmul(mo_occ.T,
                             np.matmul(fock_lambda.alpha_to_numpy(0), mo_vir)) +
                   np.matmul(mo_vir.T,
                             np.matmul(fock_lambda.alpha_to_numpy(0), mo_occ)).T)
                 + v.reshape(nocc, nvir) * eov)

             # increase iteration counter every time this function is called
             self.iter_count += 1

             # Write current lambda to checkpoint file
             #if (self.rank == mpi_master() and self.iter_count % n_restart_iterations == 0):
             #    print("What we write to the checkpoint file:")
             #    print("self.iter_count =", self.iter_count)
             #    print("self.rhs_mo =\n", self.rhs_mo)
             #    print("lambda_mo =\n", lambda_mo)
             #    write_rsp_hdf5(self.checkpoint_file, [self.rhs_mo, lambda_mo],
             #                   ['OrbRsp_RHS', 'OrbRsp_lambda'], molecule, basis,
             #                   dft_dict, pe_dict, self.ostream)

             return lambda_mo.reshape(nocc * nvir)

        # 5) Define the linear operator and run conjugate gradient
        LinOp = linalg.LinearOperator((nocc * nvir, nocc * nvir),
                                      matvec=OrbRsp_MatVec)

        profiler.start_timer(0, 'Conjugate Gradient')

        lambda_multipliers, cg_conv = linalg.cg(
            A=LinOp,
            b=self.rhs_mo.reshape(nocc * nvir),
            x0=lambda_guess.reshape(nocc * nvir),
            tol=self.conv_thresh, atol=0,
            maxiter=self.max_iter)

        if cg_conv == 0:
            self.is_converged = True
        # TODO: print warning or something if not converged

        profiler.check_memory_usage('Conjugate Gradient')
        profiler.stop_timer(0, 'Conjugate Gradient')

        # Transform to AO
        self.lambda_ao = np.matmul(
            mo_occ, np.matmul(lambda_multipliers.reshape(nocc, nvir), mo_vir.T))

        # Calculate the relaxed one-particle density matrix
        # Factor 4: (ov + vo)*(alpha + beta)
        self.rel_dm_ao = self.unrel_dm_ao + 4 * self.lambda_ao

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

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


