import numpy as np
import time as tm

from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import mpi_master
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .profiler import Profiler
from .orbitalresponse import OrbitalResponse
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_scheme
from .checkpoint import read_rsp_hdf5
from .checkpoint import write_rsp_hdf5


class TdaOrbitalResponse(OrbitalResponse):
    """
    Implements orbital response Lagrange multipliers computation using a
    conjugate gradient scheme for the Tamm-Dancoff Approximation (TDA)
    level of theory.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - n_state_deriv: The number of the excited state of interest.
        - solver: The linear equations solver.
    """

    # Give number of state of interest and all vectors in compute function?
    # Or only the specific vector in compute function?
    def __init__(self, comm, ostream):
        """
        Initializes orbital response computation driver to default setup.
        """

        # Settings from OrbitalResponse parent class 
        super().__init__(comm, ostream)


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

		# Update setting in parent class
        super().update_settings(rsp_dict, method_dict)


    def compute(self, molecule, basis, scf_tensors, excitation_vecs):
        """
        Performs orbital response Lagrange multipliers calculation using molecular data.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.
        :param scf_tensors:
            The dictionary of tensors from converged SCF calculation.
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
        # eri_dict = self.init_eri(molecule, basis)

        # DFT information
        # dft_dict = self.init_dft(molecule, scf_tensors)

        # PE information
        # pe_dict = self.init_pe(molecule, basis)

        # timing_dict = {}

        # block Davidson algorithm setup

        #self.solver = BlockDavidsonSolver() # something with linalg.cg?

        # read initial guess from restart file

        # It might be possible to restart by using the current lambda
        # and the right-hand side
        #n_restart_vectors = 0
        #n_restart_iterations = 0

        #if self.restart:
        #    if self.rank == mpi_master():
        #        rst_trial_mat, rst_sig_mat = read_rsp_hdf5(
        #            self.checkpoint_file, ['TDA_trials', 'TDA_sigmas'],
        #            molecule, basis, dft_dict, pe_dict, self.ostream)
        #        self.restart = (rst_trial_mat is not None and
        #                        rst_sig_mat is not None)
        #        if rst_trial_mat is not None:
        #            n_restart_vectors = rst_trial_mat.shape[1]
        #    self.restart = self.comm.bcast(self.restart, root=mpi_master())
        #    n_restart_vectors = self.comm.bcast(n_restart_vectors,
        #                                        root=mpi_master())
        #    n_restart_iterations = n_restart_vectors // self.nstates
        #    if n_restart_vectors % self.nstates != 0:
        #        n_restart_iterations += 1

        # TODO
        # 1) Construct the necessary density matrices
        # 2) Construct the RHS
        # 3) Construct the initial guess
        # 4) Write the linear operator for matrix-vector product
        # 5) Run the conjugate gradient

        profiler.start_timer(0, 'Orbital Response')

        # 1) Calculate unrelaxed one-particle and transition density matrix
        nocc = molecule.number_of_alpha_electrons()
        mo_occ = scf_tensors['C'][:, :nocc]  # occupied MO coefficients
        mo_vir = scf_tensors['C'][:, nocc:]  # virtual MO coefficients
        nvir = mo_vir.shape[1]
        ovlp = scf_tensors['S']  # overlap matrix

        # Take vector of interest and convert to matrix form
        exc_vec = excitation_vecs[:, self.n_state_deriv].copy().reshape(
            nocc, nvir)

        # Calcuate the unrelaxed one-particle density matrix in MO basis
        dm_oo = -np.einsum('ia,ja->ij', exc_vec, exc_vec)
        dm_vv = np.einsum('ia,ib->ab', exc_vec, exc_vec)

        # Transform unrelaxed one-particle density matrix to the AO basis
        self.dm_unrel_ao = (np.matmul(mo_occ, np.matmul(dm_oo, mo_occ.T)) +
                       np.matmul(mo_vir, np.matmul(dm_vv, mo_vir.T)))

        # Transform the excitation vectors to the AO basis
        exc_vec_ao = np.matmul(mo_occ, np.matmul(exc_vec, mo_vir.T))

        # 2) Construct the right-hand side
        dm_ao_rhs = AODensityMatrix([dm_unrel_ao, exc_vec_ao], denmat.rest)
        fock_ao_rhs = AOFockMatrix(dm_ao_rhs)
        fock_flag = fockmat.rgenjk
        fock_ao_rhs.set_fock_type(fock_flag, 1)

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        eri_drv.compute(fock_ao_rhs, dm_ao_rhs, molecule, basis, screening)

        # Calculate the RHS and transform it to the MO basis
        self.rhs_mo = (np.einsum(
            'pi,pa->ia', mo_occ,
            np.einsum('ta,pt->pa', mo_vir, 0.5 * fock_ao_rhs.alpha_to_numpy(0)))
                  + np.einsum(
                      'mi,ma->ia', mo_occ,
                      np.einsum(
                          'za,mz->ma', mo_vir,
                          np.einsum(
                              'mr,rz->mz', ovlp,
                              np.einsum('rp,zp->rz', exc_vec_ao,
                                        0.5 * fock_ao_rhs.alpha_to_numpy(1))) -
                          np.einsum(
                              'rz,mr->mz', ovlp,
                              np.einsum('pr,mp->mr', exc_vec_ao, 0.5 *
                                        fock_ao_rhs.alpha_to_numpy(1).T)))))


        # Calculate the lambda multipliers and the relaxed one-particle density
		# in the parent class
        lambda_multipliers = self.compute_lambda(molecule, basis, scf_tensors)

        profiler.stop_timer(0, 'Orbital Response')
        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        # Calculate the overlap matrix multipliers

		# TODO: As soon as we figured out how to do everything in AO basis,
		# this should be put in a separate routine
        omega_VV = -(np.einsum(
            'ta,tb->ab', mo_vir,
            np.einsum('rb,tr->tb', mo_vir, (np.einsum(
                'zr,tz->tr', ovlp,
                np.einsum('nz,tn->tz', exc_vec_ao,
                          0.5 * fock_ao_rhs.alpha_to_numpy(1).T))))) +
                     np.einsum('a,ab->ab', evir, 0.5 * dm_vv))

        omega_OV = -(
            np.einsum(
                'mi,ma->ia', mo_occ,
                np.einsum('za,mz->ma', mo_vir, (np.einsum(
                    'rz,mr->mz', ovlp,
                    np.einsum('tr,mt->mr', exc_vec_ao,
                              0.5 * fock_ao_rhs.alpha_to_numpy(1).T))))) +
            np.einsum('i,ia->ia', eocc, lambda_multipliers))

        #TODO: We need the object from inside the matrix-vector product.
        # We recompute it here. Should it be defined outside the function instead?

        #Create AODensityMatrix object from lambda in AO
        ao_density_lambda = AODensityMatrix([self.lambda_ao], denmat.rest)

        #Create a Fock Matrix Object (initialized with zeros)
        fock_lambda = AOFockMatrix(ao_density_lambda)
        fock_lambda.set_fock_type(fock_flag, 0)
        eri_drv.compute(fock_lambda, ao_density_lambda, molecule, basis,
                        screening)

        omega_OO = -(
            np.einsum(
                'mi,mj->ij', mo_occ,
                np.einsum('zj,mz->mj', mo_occ, (np.einsum(
                    'zr,mr->mz', ovlp,
                    np.einsum('rp,mp->mr', exc_vec_ao, 0.5 *
                              fock_ao_rhs.alpha_to_numpy(1)))))) + np.einsum(
                                  'pi,pj->ij', mo_occ,
                                  np.einsum('tj,pt->pj', mo_occ, 0.5 *
                                            fock_ao_rhs.alpha_to_numpy(0))) +
            np.matmul(mo_occ.T, np.matmul(fock_lambda.alpha_to_numpy(0),
                                          mo_occ)) +
            np.matmul(mo_occ.T, np.matmul(fock_lambda.alpha_to_numpy(0),
                                          mo_occ)).T + np.diag(eocc) +
            np.einsum('i,ij->ij', eocc, 0.5 * dm_oo))

        self.omega_ao = (np.matmul(mo_occ, np.matmul(omega_OO, mo_occ.T)) +
                    np.matmul(mo_vir, np.matmul(omega_VV, mo_vir.T)) +
                    np.matmul(mo_occ, np.matmul(omega_OV, mo_vir.T)) +
                    np.matmul(mo_vir, np.matmul(omega_OV.T, mo_occ.T)))

        profiler.check_memory_usage('End of Orbital Response Driver')
        profiler.print_memory_usage(self.ostream)

        if self.rank == mpi_master() and self.is_converged:
            self.ostream.print_blank()
            self.ostream.flush()



##    def compute_omega_mult(self, molecule, scf_tensors, exc_vectors):
##        """
##        Calculates the Lagrange multipliers for the overlap matrix.
##
##        :param molecule:
##            The molecule.
##        :param scf_tensors:
##            The tensors from the converged SCF calculation.
##        :param exc_vectors:
##            The excitation vectors from the TDA calculation.
##
##        :return:
##            A dictionary containing the Lagrange multipliers.
##        """
##        # Orbital information
##        nocc = molecule.number_of_alpha_electrons()
##        mo_occ = scf_tensors['C'][:, :nocc] # occupied MO coefficients
##        mo_vir = scf_tensors['C'][:, nocc:] # virtual MO coefficients
##        nvir = mo_vir.shape[1]
##        ovlp = scf_tensors['S'] # overlap matrix
##        eocc = scf_tensors['E'][:nocc]
##        evir = scf_tensors['E'][nocc:]
##
##        # Take vector of interest and convert to matrix form
##        exc_vec = exc_vectors[:, self.n_state_deriv].copy().reshape(nocc, nvir)
##
##        #TODO: to be continued...
