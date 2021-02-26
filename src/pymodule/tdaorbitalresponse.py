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

        print("self.restart before update_settings in linearsolver:\n", self.restart)
        # Update setting in parent class
        super().update_settings(rsp_dict, method_dict)
        print("self.restart after update_settings in linearsolver:\n", self.restart)


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
            The excitation vectors from converged TDA calculation.

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
        eocc = scf_tensors['E'][:nocc]
        evir = scf_tensors['E'][nocc:]

        # Take vector of interest and convert to matrix form
        exc_vec = excitation_vecs[:, self.n_state_deriv].copy().reshape(
            nocc, nvir)

        # Transform the excitation vectors to the AO basis
        exc_vec_ao = np.matmul(mo_occ, np.matmul(exc_vec, mo_vir.T))

        # Calcuate the unrelaxed one-particle density matrix in MO basis
        dm_oo = -np.einsum('ia,ja->ij', exc_vec, exc_vec)
        dm_vv = np.einsum('ia,ib->ab', exc_vec, exc_vec)

        # Transform unrelaxed one-particle density matrix to the AO basis
        self.unrel_dm_ao = (np.matmul(mo_occ, np.matmul(dm_oo, mo_occ.T)) +
                       np.matmul(mo_vir, np.matmul(dm_vv, mo_vir.T)))

        # 2) Construct the right-hand side
        dm_ao_rhs = AODensityMatrix([self.unrel_dm_ao, exc_vec_ao], denmat.rest)
        fock_ao_rhs = AOFockMatrix(dm_ao_rhs)
        fock_flag = fockmat.rgenjk
        fock_ao_rhs.set_fock_type(fock_flag, 1)

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        eri_drv.compute(fock_ao_rhs, dm_ao_rhs, molecule, basis, screening)

        print("tdaorbitalresponse.py:")
        #TODO: save RHS in checkpoint file or not?
        if not self.restart:
            print("self.restart was false, calculating the RHS")
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
        else:
            print("self.restart was true, so RHS not calculated")


        # Calculate the lambda multipliers and the relaxed one-particle density
        # in the parent class
        lambda_multipliers = self.compute_lambda(molecule, basis, scf_tensors)


        # Calculate the overlap matrix multipliers
        # 1. compute an energy-weighted density matrix
        # (needed to compute omega)
        epsilon_dm_ao = (np.matmul(mo_occ,
                                   np.matmul((np.einsum("i,ij->ij",eocc,0.5*dm_oo)
                                              +np.diag(eocc)),
                                              mo_occ.T))
                        +np.matmul(mo_vir,
                                   np.matmul(np.einsum("a,ab->ab",evir,0.5*dm_vv),
                                             mo_vir.T)
                                    )
                        +np.matmul(mo_occ,
                                   np.matmul(np.einsum("i,ia->ia",eocc,
                                                        lambda_multipliers),
                                             mo_vir.T))
                        +np.matmul(mo_occ,
                                   np.matmul(np.einsum("i,ia->ia",eocc,
                                                        lambda_multipliers),
                                             mo_vir.T)).T

                        )

        # 2. compute the omega multipliers in AO basis:
        self.omega_ao = self.compute_omega(molecule, basis, scf_tensors,
                                           epsilon_dm_ao,
                                           exc_vec_ao, fock_ao_rhs) 

        profiler.stop_timer(0, 'Orbital Response')
        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)

        profiler.check_memory_usage('End of Orbital Response Driver')
        profiler.print_memory_usage(self.ostream)
        profiler.print_memory_tracing(self.ostream)

        if self.rank == mpi_master() and self.is_converged:
            self.ostream.print_blank()
            self.ostream.flush()



    def compute_omega(self, molecule, basis, scf_tensors,
                      epsilon_dm_ao, exc_vec_ao, fock_ao_rhs):
        """
        Calculates the Lagrange multipliers for the overlap matrix.

        :param molecule:
            The molecule
        :param basis:
            The basis set
        :param epsilon_dm_ao:
            The energy-weighted relaxed density matrix
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param exc_vec_ao:
            The excitation vector of interest in AO basis
        :param fock_ao_rhs:
            The AOFockMatrix from the right-hand side of the orbital response eq. 

        :return:
            a numpy array containing the Lagrange multipliers in AO basis.
        """

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                                    self.eri_thresh, molecule, basis)

        ao_density_lambda = AODensityMatrix([self.lambda_ao], denmat.rest)

        #Create a Fock Matrix Object (initialized with zeros)
        fock_flag = fockmat.rgenjk
        fock_lambda = AOFockMatrix(ao_density_lambda)
        fock_lambda.set_fock_type(fock_flag, 0)
        eri_drv.compute(fock_lambda, ao_density_lambda, molecule, basis,
                        screening)

        nocc = molecule.number_of_alpha_electrons()
        mo_vir = scf_tensors['C'][:, nocc:]  # virtual MO coefficients
        ovlp = scf_tensors['S']

        # The density matrix; only alpha block;
        # Only works for the restricted case
        D_occ = scf_tensors['D'][0]
        D_vir = np.matmul(mo_vir, mo_vir.T) 

        # Because the excitation vector is not symmetric,
        # we need both the matrix (OO block in omega, and probably VO)
        # and its transpose (VV, OV blocks)
        # this comes from the transformation of the 2PDM contribution
        # from MO to AO basis
        Ft = ( np.einsum('zr,tz->tr', ovlp,
                  np.einsum('nz,tn->tz', exc_vec_ao, 
                               0.5*fock_ao_rhs.alpha_to_numpy(1).T
                            )
                 )
            )

        F = ( np.einsum('zr,mr->mz', ovlp,
                   np.einsum('rp,mp->mr', exc_vec_ao, 
                                 0.5*fock_ao_rhs.alpha_to_numpy(1)
                                )
                    )
            )
                        
 
        # Compute the contributions from the 2PDM and the relaxed 1PDM
        # to the omega Lagrange multipliers:
        O_1pdm_2pdm_contribs = (np.matmul(D_occ, np.matmul(F, D_occ))
                                +np.matmul(D_occ, np.matmul(Ft, D_vir))
                                +np.matmul(D_occ, np.matmul(Ft, D_vir)).T
                                +np.matmul(D_vir, np.matmul(Ft, D_vir))
                                +np.matmul(D_occ,
                                           np.matmul((fock_lambda.alpha_to_numpy(0)
                                             +fock_lambda.alpha_to_numpy(0).T
                                             +0.5*fock_ao_rhs.alpha_to_numpy(0)),
                                             D_occ)
                                    )
                                )

        omega = -epsilon_dm_ao - O_1pdm_2pdm_contribs

        return omega

