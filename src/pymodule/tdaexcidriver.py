from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import MolecularOrbitals
from .veloxchemlib import ScreeningContainer
from .veloxchemlib import ExcitationVector
from .veloxchemlib import TDASigmaVectorDriver
from .veloxchemlib import mpi_master
from .veloxchemlib import ericut
from .veloxchemlib import molorb
from .veloxchemlib import szblock

from .outputstream import OutputStream
from .blockdavidson import BlockDavidsonSolver

import numpy as np
import time as tm

class TDAExciDriver:
    """Implements TDA excited states computation driver.
        
    Implements TDA excited states computation schheme for Hatree-Fock/Kohn-Sham
    level of theory.
    
    Attributes
    ----------
    nstates
        The number of excited states determined by driver.
    triplet
        The triplet excited states flag.
    eri_thresh
        The electron repulsion integrals screening threshold.
    conv_thresh
        The excited states convergence threshold.
    max_iter
        The maximum number of excited states driver iterations.
    cur_iter
        The current number of excited states driver iterations.
    solver
        The eigenvalues solver.
    is_converged
        The flag for excited states convergence.
    rank
        The rank of MPI process.
    nodes
        The number of MPI processes.
    """

    def __init__(self, rank, nodes):
        """Initializes TDA excited states computation driver.
            
        Initializes TDA excited states computation drived to default setup.
        """
    
        # excited states information
        self.nstates = 0
        self.triplet = False
        
        # thresholds
        self.eri_thresh = 1.0e-15
        self.conv_thesh = 1.0-4
        
        # solver setup
        self.max_iter = 50
        self.cur_iter = 0
        self.solver = None
        self.is_converged = None
        
        # mpi information
        self.rank = rank
        self.nodes = nodes
    
<<<<<<< HEAD
    def set_number_states(self, nstates):
        """Sets number of excited states determined by solver.
            
        Sets number of excited states determined by solver.
        
        Parameters
        ----------
        nstates
            The number of excited states.
        """
        
        self.nstates = nstates

    def set_eri_threshold(self, eri_thresh):
        """Sets threshold in computation of electron repulsion integrals.
            
        Sets threshold in computation of electron repulsion integrals.
        
        Parameters
        ----------
        eri_thresh
            The threshold for computation of electron repulsion integrals.
        """
        self.eri_thresh = eri_thresh
    
    def set_solver(self, conv_thresh, max_iter):
        """Sets convergence threshold and maximum number of iterations.
            
        Sets convergence threshold and maximum number of iterations in excited
        states solver.
        
        Parameters
        ----------
        conv_thresh
            The convergence threshold of excited states.
        max_iter
            The maximum number of iterations.
        """
        
        self.conv_thresh = conv_thresh
        self.max_iter = max_iter
 
    def compute(self, qq_data, mol_orbs, molecule, ao_basis, comm,
                ostream=OutputStream("")):
        """Performs TDA excited states calculation.
            
        Performs TDA excited states calculation using molecular data, MPI
        communicator and output stream.
        
        Parameters
        ----------
        mol_orbs
            The molecular orbitals.
        molecule
            The molecule.
        ao_basis
            The AO basis set.
        min_basis
            The minimal AO basis set.
        comm
            The MPI communicator.
        ostream
            The output stream.
        """
       
        # set start time
       
        start_time = tm.time()
       
        # set up trial excitation vectors on master node
        
        diag_mat, trial_vecs = self.gen_trial_vectors(mol_orbs, molecule)
        
        # initalize sigma vectors driver
        
        a2x_drv = TDASigmaVectorDriver(self.rank, self.nodes, comm)

        # block Davidson algorithm setup
        
        self.solver = BlockDavidsonSolver()
        
        for i in range(self.max_iter):
            
            # perform linear transformation of trial vectors
            
            sig_vecs = a2x_drv.compute(trial_vecs, qq_data, mol_orbs, molecule,
                                       ao_basis, comm)

            # solve eigenvalues problem on master node
            
            if self.rank == mpi_master():
            
                sig_mat = self.convert_to_sigma_matrix(sig_vecs)
                trial_mat = self.convert_to_trial_matrix(trial_vecs)

                self.solver.add_iteration_data(sig_mat, trial_mat, i)
                                       
                zvecs = self.solver.compute(diag_mat)
            
                self.print_iter_data(i, ostream)
                    
            # check convergence
            
            self.check_convergence(i, comm)
            
            if self.is_converged:
                break
            
            # update trial vectors on master node
            
            if self.rank == mpi_master():
                self.update_trial_vectors(trial_vecs, zvecs)

        # print converged excited states

        if self.rank == mpi_master():
            self.print_excited_states(start_time, ostream)
        
    def gen_trial_vectors(self, mol_orbs, molecule):
        """Generates set of TDA trial vectors.
            
        Generates set of TDA trial vectors for given number of excited states by
        selecting primitive excitations wirh lowest approximate energies
        E_ai = e_a-e_i.
            
        Parameters
        ----------
        mol_orbs
            The molecular orbitals.
        molecule
            The molecule.
        Returns
        -------
            The set of trial vectors.
        """
        
        if self.rank == mpi_master():
        
            nocc = molecule.number_of_electrons() // 2
            norb = mol_orbs.number_mos()
        
            zvec = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
            exci_list = zvec.small_energy_identifiers(mol_orbs, self.nstates)
            
            diag_mat = zvec.diagonal_to_numpy(mol_orbs)

            trial_vecs = []
            for i in exci_list:
                trial_vecs.append(ExcitationVector(szblock.aa, 0, nocc, nocc,
                                                   norb, True))
                trial_vecs[-1].set_zcoefficient(1.0, i)

            return (diag_mat, trial_vecs)
    
        return (None, [])

    def update_trial_vectors(self, trial_vecs, zvecs):
        """Updates set of TDA trial vectors.
            
        Updates set of TDA trial vectors by replacing Z vector coefficients with
        values generated by reduced space solver.
            
        Parameters
        ----------
        trial_vecs
            The set of trial vectors.
        zvecs
            The set of Z vector coefficients.
        """
    
        for i in range(zvecs.shape[1]):
            for j in range(zvecs.shape[0]):
                trial_vecs[i].set_zcoefficient(zvecs[j,i], j)

    def check_convergence(self, iter, comm):
        """Checks convergence of excitation energies and set convergence flag.
            
        Checks convergence of excitation energies and set convergence flag on
        all processes within MPI communicator.
            
        Parameters
        ----------
        iter
            The current excited states solver iteration.
        comm
            The MPI communicator.
        """
        
        self.cur_iter = iter
        
        if self.rank == mpi_master():
            self.is_converged = self.solver.check_convergence(self.conv_thresh)
        else:
            self.is_converged = False
    
        self.is_converged = comm.bcast(self.is_converged, root=mpi_master())

    def convert_to_sigma_matrix(self, sig_vecs):
        """Converts set of sigma vectors from C++ data format to numpy data
        format.
            
        Converts set of sigma vectors from std::vector<CDenseMatrix> to numpy
        2D array.
            
        Parameters
        ----------
        sig_vecs
            The sigma vectors as std::vector<CDenseMatrix>.
        Returns
        -------
            The 2D numpy array.
        """
    
        nvecs = len(sig_vecs)
        
        if nvecs > 0:
            
            sig_mat = sig_vecs[0].to_numpy()
            for i in range(1, nvecs):
                sig_mat = np.hstack((sig_mat, sig_vecs[i].to_numpy()))

            return sig_mat
    
        return None
    
    def convert_to_trial_matrix(self, trial_vecs):
        """Converts set of Z vectors from C++ data format to numpy data
        format.
        
        Converts set of Z vectors from std::vector<CExcitationVector> to numpy
        2D array.
        
        Parameters
        ----------
        trial_vecs
            The Z vectors as std::vector<CExcitationVector>.
        Returns
        -------
            The 2D numpy array.
        """
    
        nvecs = len(trial_vecs)
    
        if nvecs > 0:
    
            trial_mat = trial_vecs[0].zvector_to_numpy()
            for i in range(1, nvecs):
                trial_mat = np.hstack((trial_mat,
                                       trial_vecs[i].zvector_to_numpy()))

            return trial_mat

        return None

    def print_iter_data(self, iter, ostream):
        """Prints excited states solver iteration data to output stream.
            
        Prints excited states solver iteration data to output stream.
        
        Parameters
        ----------
        iter
            The current excited states solver iteration.
        ostream
            The output stream.
        """
    
        # iteration header
        
        exec_str  = " *** Iteration: " + (str(iter + 1)).rjust(3)
        exec_str += " * Reduced Space: "
        exec_str += (str(self.solver.reduced_space_size())).rjust(4)
        rmin, rmax = self.solver.max_min_residual_norms()
        exec_str += " * Residues (Max,Min): {:.2e} and {:.2e}".format(rmax,rmin)
        ostream.print_header(exec_str)
        ostream.print_blank()
        
        # excited states information
            
        reigs, rnorms = self.solver.get_eigenvalues()
        for i in range(reigs.shape[0]):
            exec_str  = "State {:2d}: {:5.8f} ".format(i + 1, reigs[i])
            exec_str += "au Residual Norm: {:3.8f}".format(rnorms[i])
            ostream.print_header(exec_str.ljust(84))
        
        # flush output stream
        ostream.print_blank()
        ostream.flush()

    def print_excited_states(self, start_time, ostream):
        """Prints excited states information to output stream.
        
        Prints excited states information to output stream.
        
        Parameters
        ----------
        start_time
            The start time of SCF calculation.
        ostream
            The output stream.
        """
        
        valstr = "*** {:d} excited states ".format(self.nstates)
        if self.is_converged:
            valstr += "converged"
        else:
            valstr += "not converged"
        valstr += " in {:d} iterations. ".format(self.cur_iter + 1)
        valstr += "Time: {:.2f}".format(tm.time() - start_time) + " sec."
        ostream.print_blank()
        ostream.print_header(valstr.ljust(92))
        ostream.print_blank()







