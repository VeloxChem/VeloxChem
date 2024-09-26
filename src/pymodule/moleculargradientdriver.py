import numpy as np

from .veloxchemlib import OverlapGeom100Driver
from .veloxchemlib import KineticEnergyGeom100Driver
from .veloxchemlib import NuclearPotentialGeom100Driver
from .veloxchemlib import NuclearPotentialGeom010Driver
from .veloxchemlib import FockGeom1000Driver
from .veloxchemlib import make_matrix
from .veloxchemlib import mat_t
from .veloxchemlib import partition_atoms

class MolecularGradientDriver:
    """
    Implements the molecular gradient driver.
    """
    
    def __init__(self):
        """
        Initializes the molecular gradient driver.
        """
        
        self.comm = None
        
    def mpi_comp_electronic_grad(self, comm, molecule, basis, density, wdensity, xcfactor, omega):
        """
        Computes electronic part of molecular gradient.

        :param comm:
            The MPI communicator.
        :param molecule:
            The molecule to compute gradient.
        :param basis:
            The basis set to compute gradient.
        :param density:
            The density matrix to compute gradient.
        :param wdensity:
            The weighted density matrix to compute gradient.
        :param xcfactor: 
            The exchange matrix scaling factor.
        :param omega:
            The range separation factor.

        :return:
            The electronic part of molecular gradient.
        """

        natoms = molecule.number_of_atoms()
        grad_mat = np.zeros((natoms, 3))
        
        for i in partition_atoms(natoms, comm.Get_rank(), comm.Get_size()):
            self.comp_kin_grad(grad_mat, molecule, basis, density, i)
            self.comp_npot_grad(grad_mat, molecule, basis, density, i)
            self.comp_orb_grad(grad_mat, molecule, basis, wdensity, i)
            self.comp_fock_grad(grad_mat, molecule, basis, density, xcfactor, omega, i)
        
        return grad_mat

    def comp_electronic_grad(self, molecule, basis, density, wdensity, xcfactor, omega):
        """
        Computes electronic part of molecular gradient.

        :param molecule:
            The molecule to compute gradient.
        :param basis:
            The basis set to compute gradient.
        :param density:
            The density matrix to compute gradient.
        :param wdensity:
            The weighted density matrix to compute gradient.
        :param xcfactor: 
            The exchange matrix scaling factor.
        :param omega:
            The range separation factor.

        :return:
            The electronic part of molecular gradient.
        """

        natoms = molecule.number_of_atoms()
        grad_mat = np.zeros((natoms, 3))
                
        for i in range(natoms):
            #self.comp_kin_grad(grad_mat, molecule, basis, density, i)
            #self.comp_npot_grad(grad_mat, molecule, basis, density, i)
            #self.comp_orb_grad(grad_mat, molecule, basis, wdensity, i)
            self.comp_fock_grad(grad_mat, molecule, basis, density, xcfactor, omega, i)
                
        return grad_mat
        
    def comp_kin_grad(self, grad_mat, molecule, basis, density, iatom):
        """
        Computes kinetic energy contribution to electronic gradient for specific atom.

        :param molecule:
            The molecule to compute gradient.
        :param basis:
            The basis set to compute gradient.
        :param density:
            The density matrix to compute gradient.
        :param iatom:
            The selected atom to compute gradient.
        """

        grad_drv = KineticEnergyGeom100Driver()
        gmats = grad_drv.compute(molecule, basis, iatom)
        for i, label in enumerate(['X', 'Y', 'Z']):
            gmat = gmats.matrix(label).full_matrix().to_numpy()
            grad_mat[iatom, i] += 2.0 * np.trace(np.matmul(gmat + gmat.T, density))
        
    def comp_npot_grad(self, grad_mat, molecule, basis, density, iatom):
        """
        Computes nuclear potential contribution to electronic gradient for specific atom.

        :param molecule:
            The molecule to compute gradient.
        :param basis:
            The basis set to compute gradient.
        :param density:
            The density matrix to compute gradient.
        :param iatom:
            The selected atom to compute gradient.
        """

        grad_drv = NuclearPotentialGeom100Driver()
        gmats = grad_drv.compute(molecule, basis, iatom)
        for i, label in enumerate(['X', 'Y', 'Z']):
            gmat = gmats.matrix(label).full_matrix().to_numpy()
            grad_mat[iatom, i] -= 2.0 * np.trace(np.matmul(gmat + gmat.T, density))
        
        grad_drv = NuclearPotentialGeom010Driver()
        gmats = grad_drv.compute(molecule, basis, iatom)
        for i, label in enumerate(['X', 'Y', 'Z']):
            gmat = gmats.matrix(label).full_matrix().to_numpy()
            grad_mat[iatom, i] -= 2.0 * np.trace(np.matmul(gmat, density))
            
    def comp_orb_grad(self, grad_mat, molecule, basis, density, iatom):
        """
        Computes molecular orbitals contribution to electronic gradient for specific atom.

        :param molecule:
            The molecule to compute gradient.
        :param basis:
            The basis set to compute gradient.
        :param density:
            The weighted density matrix to compute gradient.
        :param iatom:
            The selected atom to compute gradient.
        """
        
        grad_drv = OverlapGeom100Driver()
        gmats = grad_drv.compute(molecule, basis, iatom)
        for i, label in enumerate(['X', 'Y', 'Z']):
            gmat = gmats.matrix(label).full_matrix().to_numpy()
            grad_mat[iatom, i] += 2.0 * np.trace(np.matmul(gmat + gmat.T, density))
            
    def comp_fock_grad(self, grad_mat, molecule, basis, density, xcfactor, omega, iatom):
        """
        Computes Fock matrix contribution to electronic gradient for specific atom.

        :param molecule:
            The molecule to compute gradient.
        :param basis:
            The basis set to compute gradient.
        :param density:
            The density matrix to compute gradient. 
        :param xcfactor: 
            The exchange matrix scaling factor.
        :param omega:
            The range separation factor.
        :param iatom:
            The selected atom to compute gradient.
        """
        
        den_mat = make_matrix(basis, mat_t.symmetric)
        den_mat.set_values(density)
        
        fock_drv = FockGeom1000Driver()
        if xcfactor > 1.0e-10:
            fmats = fock_drv.compute(basis, molecule, den_mat, iatom, "2jkx", xcfactor, 0.0)
        else:
            fmats = fock_drv.compute(basis, molecule, den_mat, iatom, "2jk", 0.0, 0.0)
            
        for i, label in enumerate(['X', 'Y', 'Z']):
            fmat = fmats.matrix(label).full_matrix().to_numpy()
            grad_mat[iatom, i] += np.trace(np.matmul(fmat, density))
