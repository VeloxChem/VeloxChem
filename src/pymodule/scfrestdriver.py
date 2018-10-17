from .VeloxChemLib import mpi_master

from .scfdriver import ScfDriver

import numpy as np

class ScfRestrictedDriver(ScfDriver):

    def __init__(self):
        super().__init__()
    
    def comp_energy(self, fockmat, kinmat, npotmat, denmat):
        
        dmat = denmat.total_to_numpy(0)
        
        # electronic energy
        
        gmat = fockmat.to_numpy(0)
        gd = np.matmul(gmat, dmat)
        eelec = gd.trace()
        
        # kinetic energy
        
        kmat = kinmat.to_numpy()
        kd = np.matmul(kmat, dmat)
        ekin = 2.0 * kd.trace()
        
        # nuclear potential energy
        
        npmat = npotmat.to_numpy()
        npd = np.matmul(npmat, dmat)
        enpot = -2.0 * npd.trace()
        
        return (eelec, ekin, enpot)
    
    def comp_full_fock(self, fockmat, kinmat, npotmat):
        fockmat.add_hcore(kinmat, npotmat, 0)

    def comp_gradient(self, fockmat, ovlmat, denmat):
        
        dmat = denmat.total_to_numpy(0)
        smat = ovlmat.to_numpy()
        fmat = fockmat.to_numpy(0)
        
        # compute SDF and FDS matrices
        
        fa = np.matmul(smat, np.matmul(dmat, fmat))
        fb = np.matmul(fmat, np.matmul(dmat, smat))
        
        return 2.0 * np.linalg.norm(fa - fb)

    def comp_density_change(self, newden, oldden):
        
        ndmat = newden.total_to_numpy(0)
        odmat = oldden.total_to_numpy(0)
        
        return np.linalg.norm(ndmat - odmat)
    
    def gen_molecular_orbitals(self, fockmat, oaomat, ostream):
    
        fmat = fockmat.to_numpy(0)
        smat = oaomat.to_numpy()
    
        fmat_mo = np.matmul(smat.transpose(), np.matmul(fmat, smat))
    
        evals, evecs = np.linalg.eigh(fmat_mo)
    
        morbs = np.matmul(smat, evecs)
    
        #print(evals)
    
        #print(morbs)
    
    def get_scf_type(self):
        return "Spin-Restricted Hatree-Fock"
