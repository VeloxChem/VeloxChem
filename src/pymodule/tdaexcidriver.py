from .veloxchemlib import ElectronRepulsionIntegralsDriver
from .veloxchemlib import MolecularOrbitals
from .veloxchemlib import ScreeningContainer
from .veloxchemlib import ExcitationVector
from .veloxchemlib import TDASigmaVectorDriver
from .veloxchemlib import mpi_master
from .veloxchemlib import ericut
from .veloxchemlib import molorb
from .veloxchemlib import szblock

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
        self.max_iter = 2
        
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
       
        # set up excitation vectors on master node
        
        trial_vecs = self.gen_zvectors(mol_orbs, molecule)
    
        # initalize sigma vectors driver
        
        a2x_drv = TDASigmaVectorDriver(self.rank, self.nodes, comm)

        # block Davidson algorithm
        
        sig_vecs = []
        evecs = []
        evals = []
        
        for i in range(self.max_iter):
            
            a2x_vecs = a2x_drv.compute(trial_vecs, qq_data, mol_orbs, molecule,
                                       ao_basis, comm)
            
            sig_vecs.extend(a2x_vecs)
            evecs.extend(trial_vecs)
            
            rl_mat = self.comp_rayleigh_matrix(sig_vecs, evecs)
            print("Rayleight matrix:")
            print(rl_mat)
            self.comp_ritz_vector(rl_mat)


        
    def gen_zvectors(self, mol_orbs, molecule):
        
        trial_vecs = []
        
        if self.rank == mpi_master():
        
            nocc = molecule.number_of_electrons() // 2
            norb = mol_orbs.number_mos()
        
            zvec = ExcitationVector(szblock.aa, 0, nocc, nocc, norb, True)
            exci_list = zvec.small_energy_identifiers(mol_orbs, self.nstates)


            for i in exci_list:
                trial_vecs.append(ExcitationVector(szblock.aa, 0, nocc, nocc,
                                               norb, True))
                trial_vecs[-1].set_zcoefficient(1.0, i)

        return trial_vecs

    def set_number_states(self, nstates):
        self.nstates = nstates

    def set_eri_threshold(self, eri_thresh):
        self.eri_thresh = eri_thresh

    def set_solver(self, conv_thresh, max_iter):
        self.conv_thresh = conv_thresh
        self.max_iter = max_iter

    def comp_rayleigh_matrix(self, sig_vecs, evecs):

        rdim = len(sig_vecs)
    
        rmat = np.zeros(shape=(rdim, rdim), dtype=float)
        ridx = np.triu_indices(rdim)
    
        for i, j in zip(ridx[0], ridx[1]):
            
            fij = evecs[i].dot_z_matrix(sig_vecs[j])
        
            rmat[i][j] = fij
            rmat[j][i] = fij
        
        return rmat

    def comp_ritz_vector(self, rl_mat):

        eigs, evecs = np.linalg.eigh(rl_mat)

        print("Solve Rayleigh matrix:")
        print(evecs)
        print("EigenValues: ", eigs[0 : self.nstates])

        for i in range(self.nstates):
            print("EigenVector(", i, "):")
            print(evecs[:, i])







