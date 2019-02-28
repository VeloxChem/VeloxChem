from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import MolecularOrbitals
from .veloxchemlib import ScreeningContainer
from .veloxchemlib import ExcitationVector
from .veloxchemlib import TDASigmaVectorDriver
from .veloxchemlib import mpi_master
from .veloxchemlib import ericut
from .veloxchemlib import molorb
from .veloxchemlib import szblock

from .blockdavidson import BlockDavidsonSolver

import numpy as np

class TDAExciDriver:
    """Implements TDA excited states computation driver.
        
    Implements TDA excited states computation schheme for Hatree-Fock/Kohn-Sham
    level of theory.
    
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
        self.solver = None
        self.is_converged = None
        
        # mpi information
        self.rank = rank
        self.nodes = nodes
    
    def compute(self, qq_data, mol_orbs, molecule, ao_basis, comm, ostream):
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
            
            self.check_convergence(comm)
            
            if self.is_converged:
                break
            
            # update trial vectors on master node
            
            if self.rank == mpi_master():
                self.update_trial_vectors(trial_vecs, zvecs)

        # print converged excited states
        
    def gen_trial_vectors(self, mol_orbs, molecule):
        
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
    
        for i in range(zvecs.shape[1]):
            for j in range(zvecs.shape[0]):
                trial_vecs[i].set_zcoefficient(zvecs[j,i], j)

    def check_convergence(self, comm):
        
        if self.rank == mpi_master():
            self.is_converged = self.solver.check_convergence(self.conv_thresh)
        else:
            self.is_converged = False
    
        self.is_converged = comm.bcast(self.is_converged, root=mpi_master())

    def set_number_states(self, nstates):
        self.nstates = nstates

    def set_eri_threshold(self, eri_thresh):
        self.eri_thresh = eri_thresh

    def set_solver(self, conv_thresh, max_iter):
        self.conv_thresh = conv_thresh
        self.max_iter = max_iter

    def convert_to_sigma_matrix(self, sig_vecs):
    
        nvecs = len(sig_vecs)
        
        if nvecs > 0:

            sig_mat = sig_vecs[0].to_numpy()
            
            for i in range(1, nvecs):
                sig_mat = np.hstack((sig_mat, sig_vecs[i].to_numpy()))

            return sig_mat
    
        return None
    
    def convert_to_trial_matrix(self, trial_vecs):
    
        nvecs = len(trial_vecs)
    
        if nvecs > 0:
    
            trial_mat = trial_vecs[0].zvector_to_numpy()
            
            for i in range(1, nvecs):
                trial_mat = np.hstack((trial_mat, trial_vecs[i].zvector_to_numpy()))

            return trial_mat

        return None

    def print_iter_data(self, iter, ostream):
    
        exec_str  = " *** Iteration: " + (str(iter + 1)).rjust(3)
        exec_str += " * Reduced Space: "
        exec_str += (str(self.solver.reduced_space_size())).rjust(4)
        rmin, rmax = self.solver.max_min_residual_norms()
        exec_str += " * Residues (Max,Min): {:.2e} and {:.2e}".format(rmax,rmin)
        ostream.print_header(exec_str)
        ostream.print_blank()
            
        reigs, rnorms = self.solver.get_eigenvalues()
        for i in range(reigs.shape[0]):
            exec_str  = "State {:2d}: {:5.8f} ".format(i + 1, reigs[i])
            exec_str += "au Residual Norm: {:3.8f}".format(rnorms[i])
            ostream.print_header(exec_str.ljust(84))
            
        ostream.print_blank()
        ostream.flush()







