import numpy as np

from .veloxchemlib import OverlapGeom100Driver
from .veloxchemlib import KineticEnergyGeom100Driver
from .veloxchemlib import NuclearPotentialGeom100Driver
from .veloxchemlib import NuclearPotentialGeom010Driver
from .veloxchemlib import FockGeom1000Driver

class MolecularGradientDriver:
    """
    Implements the molecular gradient driver.
    """
    
    def __init__(self):
        """
        Initializes the molecular gradient driver.
        """
        
        self.comm = None

    def compute(self, molecule, basis, density, wdensity):
        """
        Computes molecular gradient for Hatree-Fock.

        :param molecule:
            The molecule to compute gradient.
        :param basis:
            The basis set to compute gradient.
        :param density:
            The density matrix to compute gradient.
        :param wdensity:
            The weighted density matrix to compute gradient.

        :return:
            The molecular gradient.
        """

        # set up density matrix
        den_mat = density.full_matrix().to_numpy()

        natoms = molecule.number_of_atoms()
        grad_mat = np.zeros((natoms, 3))
        
        ovl_grad_drv = OverlapGeom100Driver()
        kin_grad_drv = KineticEnergyGeom100Driver()
        npot_grad_drv = NuclearPotentialGeom100Driver()
        efield_drv = NuclearPotentialGeom010Driver()
        fock_drv = FockGeom1000Driver()
        
        for i in range(natoms):
            # kinetic energy contribution
            kin_mats = kin_grad_drv.compute(molecule, basis, i)
            for j, label in enumerate(['X', 'Y', 'Z']):
                kmat = kin_mats.matrix(label).full_matrix().to_numpy()
                grad_mat[i, j] = 2.0 * np.trace(np.matmul(kmat + kmat.T, den_mat))
            # nuclear potential contributions
            npot_mats = npot_grad_drv.compute(molecule, basis, i)
            for j, label in enumerate(['X', 'Y', 'Z']):
                npmat = npot_mats.matrix(label).full_matrix().to_numpy()
                grad_mat[i, j] -= 2.0 * np.trace(np.matmul(npmat + npmat.T, den_mat))
            # nuclear potential contributions
            efield_mats = efield_drv.compute(molecule, basis, i)
            for j, label in enumerate(['X', 'Y', 'Z']):
                efmat = efield_mats.matrix(label).full_matrix().to_numpy()
                grad_mat[i, j] -= 2.0 * np.trace(np.matmul(efmat, den_mat))
            # fock matrix contribution
            fock_mats = fock_drv.compute(basis, molecule, density, i, "2jk", 0.0, 0.0)
            for j, label in enumerate(['X', 'Y', 'Z']):
                fmat = fock_mats.matrix(label).full_matrix().to_numpy()
                grad_mat[i, j] += np.trace(np.matmul(fmat, den_mat))
            # weighted density contributions
            ovl_mats = ovl_grad_drv.compute(molecule, basis, i)
            for j, label in enumerate(['X', 'Y', 'Z']):
                omat = ovl_mats.matrix(label).full_matrix().to_numpy()
                grad_mat[i, j] += 2.0 * np.trace(np.matmul(omat + omat.T, wdensity))
                
        return grad_mat
        
