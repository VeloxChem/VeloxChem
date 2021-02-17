import numpy as np
import time as tm
import math

#from .veloxchemlib import ElectricDipoleIntegralsDriver
#from .veloxchemlib import LinearMomentumIntegralsDriver
#from .veloxchemlib import AngularMomentumIntegralsDriver
from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import AODensityMatrix
from .veloxchemlib import AOFockMatrix
from .veloxchemlib import mpi_master
#from .veloxchemlib import rotatory_strength_in_cgs
#from .veloxchemlib import molorb
from .veloxchemlib import denmat
from .veloxchemlib import fockmat
from .profiler import Profiler
from .linearsolver import LinearSolver
#from .blockdavidson import BlockDavidsonSolver
#from .molecularorbitals import MolecularOrbitals
from .errorhandler import assert_msg_critical
from .qqscheme import get_qq_scheme
from .checkpoint import read_rsp_hdf5
from .checkpoint import write_rsp_hdf5
from scipy.sparse import linalg


class OrbitalResponse(LinearSolver):
    """
    Implements orbital response Lagrange multipliers computation using a
    conjugate gradient scheme for the time-dependent Hartree-Fock (DFT)
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

        # many settings from LinearSolver (ERI, DFT, MPI, ...)
        super().__init__(comm, ostream)

        # excited state information, default to first excited state
        self.n_state_deriv = 0

        # solver setup
        self.solver = None

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

        # many settings updated in LinearSolver, ok for now!?
        super().update_settings(rsp_dict, method_dict)

        if 'n_state_deriv' in rsp_dict:
            # user gives '1' for first excited state, but internal index is 0
            self.n_state_deriv = int(rsp_dict['n_state_deriv']) - 1


    def compute(self, molecule, basis, scf_tensors, excitation_vecs):
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
            self.print_orbrsp_header('Orbital Response Driver', self.n_state_deriv)

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

        # 1) Calculate unrelaxed one-particle and transition density matrix
        nocc = molecule.number_of_alpha_electrons()
        mo_occ = scf_tensors['C'][:, :nocc] # occupied MO coefficients
        mo_vir = scf_tensors['C'][:, nocc:] # virtual MO coefficients
        nvir = mo_vir.shape[1]
        ovlp = scf_tensors['S']

        # Take vector of interest and convert to matrix form
        exc_vec = excitation_vecs[:, self.n_state_deriv].copy().reshape(nocc, nvir)

        # TODO: These densities are correct only for TDA! Needs to be adapted for RPA.
        # Calcuate the unrelaxed one-particle density matrix in MO basis
        dm_oo = -np.einsum('ia,ja->ij', exc_vec, exc_vec)
        dm_vv =  np.einsum('ia,ib->ab', exc_vec, exc_vec)

        # Transform unrelaxed one-particle density matrix to the AO basis
        dm_unrel_ao = (np.matmul(mo_occ, np.matmul(dm_oo, mo_occ.T))
                        + np.matmul(mo_vir, np.matmul(dm_vv, mo_vir.T)))

        # Transform the excitation vectors to the AO basis
        exc_vec_ao = np.matmul(mo_occ, np.matmul(exc_vec, mo_vir.T))

        # 2) Construct the right-hand side
        dm_ao_rhs = AODensityMatrix([dm_unrel_ao, exc_vec_ao], denmat.rest)
        fock_ao_rhs = AOFockMatrix(dm_ao_rhs)
        fock_flag = fockmat.rgenjk
        
        for i in range(fock_ao_rhs.number_of_fock_matrices()):
            fock_ao_rhs.set_fock_type(fock_flag, i)

        eri_drv = ElectronRepulsionIntegralsDriver(self.comm)
        screening = eri_drv.compute(get_qq_scheme(self.qq_type),
                            self.eri_thresh, molecule, basis)

        eri_drv.compute(fock_ao_rhs, dm_ao_rhs, molecule, basis, screening)

        # Calculate the RHS and transform it to the MO basis
        rhs_mo = (np.einsum('pi,pa->ia', mo_occ,
                    np.einsum('ta,pt->pa', mo_vir, 0.5*fock_ao_rhs.alpha_to_numpy(0)))
                + np.einsum('mi,ma->ia', mo_occ,
                    np.einsum('za,mz->ma', mo_vir,
                        np.einsum('mr,rz->mz', ovlp,
                            np.einsum('rp,zp->rz',
                                exc_vec_ao, 0.5*fock_ao_rhs.alpha_to_numpy(1)
                                    )
                                )
                - np.einsum('rz,mr->mz', ovlp,
                    np.einsum('pr,mp->mr',
                        exc_vec_ao, 0.5*fock_ao_rhs.alpha_to_numpy(1).T
                            )
                        )
                    )
                ))

        # 3) Calculate the initial guess for the Lagrange multipliers
        #   given by the RHS divided by orbital-energy differences
        eocc = scf_tensors['E'][:nocc]
        evir = scf_tensors['E'][nocc:]
        eov = eocc.reshape(-1, 1) - evir

        lambda_guess = rhs_mo / eov


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
            lambda_ao = np.matmul(mo_occ, np.matmul(v.reshape(nocc,nvir), mo_vir.T))
            
            #Create AODensityMatrix object from lambda in AO
            ao_density_lambda = AODensityMatrix([lambda_ao], denmat.rest)
            
            #Create a Fock Matrix Object (initialized with zeros)    
            fock_lambda = AOFockMatrix(ao_density_lambda)
            fock_lambda.set_fock_type(fock_flag, 0)
            eri_drv.compute(fock_lambda, ao_density_lambda, molecule, basis, screening) 
            
            #Transform to MO basis (symmetrized w.r.t. occ. and virt.) and add diagonal part
            lambda_mo = (-( np.matmul(mo_occ.T,
                            np.matmul(fock_lambda.alpha_to_numpy(0), mo_vir)) 
                        + np.matmul(mo_vir.T,
                            np.matmul(fock_lambda.alpha_to_numpy(0), mo_occ)).T )
                        +v.reshape(nocc,nvir) * eov) 
            
            return lambda_mo.reshape(nocc*nvir)

        # 5) Define the linear operator and run conjugate gradient
        LinOp = linalg.LinearOperator((nocc*nvir,nocc*nvir), matvec=OrbRsp_MatVec)

        profiler.start_timer(0, 'Conjugate Gradient')
        lambda_multipliers, cg_conv = linalg.cg(A=LinOp, b=rhs_mo.reshape(nocc*nvir),
                                        x0=lambda_guess.reshape(nocc*nvir),
                                        tol=1e-8, maxiter=50) # take from some variables

        profiler.check_memory_usage('Conjugate Gradient')
        profiler.stop_timer(0, 'Conjugate Gradient')

        if cg_conv == 0:
            self.is_converged = True

        # Factor 4: (ov + vo)*(alpha + beta)
        dm_rel_ao = (dm_unrel_ao + 4*np.matmul(mo_occ,
                            np.matmul(lambda_multipliers.reshape(nocc,nvir), mo_vir.T)))

        profiler.print_timing(self.ostream)
        profiler.print_profiling_summary(self.ostream)
        
        profiler.check_memory_usage('End of Orbital Response Driver')
        profiler.print_memory_usage(self.ostream) 

        if self.rank == mpi_master() and self.is_converged:
            self.ostream.print_blank()
            self.ostream.flush()
            
            return {
                'lambda_multipliers': lambda_multipliers.reshape(nocc,nvir),
                'relaxed_density': dm_rel_ao,
                'unrelaxed_density': dm_unrel_ao,
            }
                



    # do this in the parent class function? (add the extra variable)
    def print_orbrsp_header(self, title, n_state_deriv):
        self.ostream.print_blank()
        self.ostream.print_header('{:s} Setup'.format(title))
        self.ostream.print_header('=' * (len(title) + 8))
        self.ostream.print_blank()

        str_width = 60

        # print solver-specific info

        if n_state_deriv is not None:
            cur_str = 'State of interest           : ' + str(n_state_deriv + 1)
            self.ostream.print_header(cur_str.ljust(str_width))

        # print general info

        cur_str = 'Max. Number of Iterations       : ' + str(self.max_iter)
        self.ostream.print_header(cur_str.ljust(str_width))
        cur_str = 'Convergence Threshold           : {:.1e}'.format(
            self.conv_thresh)
        self.ostream.print_header(cur_str.ljust(str_width))

        self.ostream.print_blank()
        self.ostream.flush()

    def get_ao_densities(molecule, scf_tensors, excitation_vecs, n_state_deriv):
        """
        Calculate the unrelaxed one-particle density matrix
        in AO basis and transform the excitation vector of interest
        (the transition density) from MO to AO basis.
        
        :param molecule:
            The molecule object.
        :param scf_tensors:
            The tensors from the converged SCF calculation.
        :param excitation_vecs:
            The set of eigenvectors from converged RPA or TDA calculation.
        :param n_state_deriv:
            The number of the excited state of interest.

        :return:
            The unrelaxed one-particle density matrix and the excitation
            vector in AO basis.
        """
        # the compute function which calls this one knows "molecule"
        # and "scf_tensors"
        # should these be more general somewhere outside?
        nocc = molecule.number_of_alpha_electrons()
        mo_occ = scf_tensors['C'][:, :nocc] # occupied MO coefficients
        mo_vir = scf_tensors['C'][:, nocc:] # virtual MO coefficients
        nvir = mo_vir.shape[1]

        # Take vector of interest and convert to matrix form
        exc_vec = excitation_vecs[:,n_state_deriv].copy().reshape(nocc, nvir)

        # TODO: These densities are correct only for TDA! Needs to be adapted for RPA.
        # Calcuate the unrelaxed one-particle density matrix in MO basis
        dm_oo = -np.einsum('ia,ja->ij', exc_vec, exc_vec)
        dm_vv =  np.einsum('ia,ib->ab', exc_vec, exc_vec)

        # Transform unrelaxed one-particle density matrix to the AO basis
        dm_unrel_ao = (np.matmul(mo_occ, np.matmul(dm_oo, mo_occ.T))
                        + np.matmul(mo_vir, np.matmul(dm_vv, mo_vir.T)))

        # Transform the excitation vectors to the AO basis
        exc_vec_ao = np.matmul(mo_occ, np.matmul(exc_vec, mo_vir.T))

        return dm_unrel_ao, exc_vec_ao

