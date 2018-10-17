from .VeloxChemLib import mpi_master
from .VeloxChemLib import MolecularOrbitals

from .VeloxChemLib import molorb

from .scfdriver import ScfDriver

import numpy as np

class ScfRestrictedDriver(ScfDriver):

    def __init__(self):
        super().__init__()
    
    def comp_energy(self, fock_mat, kin_mat, npot_mat, den_mat):
        
        dmat = den_mat.total_to_numpy(0)
        
        # electronic energy
        
        gmat = fock_mat.to_numpy(0)
        gd = np.matmul(gmat, dmat)
        e_ee = gd.trace()
        
        # kinetic energy
        
        kmat = kin_mat.to_numpy()
        kd = np.matmul(kmat, dmat)
        e_kin = 2.0 * kd.trace()
        
        # nuclear potential energy
        
        npmat = npot_mat.to_numpy()
        npd = np.matmul(npmat, dmat)
        e_en = -2.0 * npd.trace()
        
        return (e_ee, e_kin, e_en)
    
    def comp_full_fock(self, fock_mat, kin_mat, npot_mat):
        fock_mat.add_hcore(kin_mat, npot_mat, 0)

    def comp_gradient(self, fock_mat, ovl_mat, den_mat):
        
        dmat = den_mat.total_to_numpy(0)
        smat = ovl_mat.to_numpy()
        fmat = fock_mat.to_numpy(0)
        
        # compute SDF and FDS matrices
        
        fa = np.matmul(smat, np.matmul(dmat, fmat))
        fb = np.matmul(fmat, np.matmul(dmat, smat))
        
        return 2.0 * np.linalg.norm(fa - fb)

    def comp_density_change(self, den_mat, old_den_mat):
        
        ndmat = den_mat.total_to_numpy(0)
        odmat = old_den_mat.total_to_numpy(0)
        
        return np.linalg.norm(ndmat - odmat)
    
    def gen_molecular_orbitals(self, fock_mat, oao_mat, ostream):
    
        fmat = fock_mat.to_numpy(0)
        smat = oao_mat.to_numpy()
    
        fmo = np.matmul(smat.transpose(), np.matmul(fmat, smat))
    
        eigs, evecs = np.linalg.eigh(fmo)
    
        orb_coefs = np.matmul(smat, evecs)
        
        print(orb_coefs);
        
        mol_orbs = MolecularOrbitals()
        
        return mol_orbs.from_numpy_list([orb_coefs], [eigs], molorb.rest)
    
    def gen_new_density(self, mol_orbs, molecule):
        return mol_orbs.get_rest_density(molecule)
    
    def get_scf_type(self):
        return "Spin-Restricted Hatree-Fock"
