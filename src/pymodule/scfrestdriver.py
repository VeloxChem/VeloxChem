from .veloxchemlib import mpi_master
from .veloxchemlib import MolecularOrbitals

from .veloxchemlib import molorb

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
        
        fa = np.matmul(fmat, np.matmul(dmat, smat))
        fb = np.matmul(smat, np.matmul(dmat, fmat))
        
        return 2.0 * np.linalg.norm(fa - fb)

    def comp_density_change(self, den_mat, old_den_mat):
        
        ndmat = den_mat.total_to_numpy(0)
        odmat = old_den_mat.total_to_numpy(0)
        return np.linalg.norm(ndmat - odmat)
    
    def store_diis_data(self, i, fock_mat, den_mat):
        
        if not self.skip_iter:
            if len(self.fock_matrices) == self.max_err_vecs:
                self.fock_matrices.popleft()
                self.den_matrices.popleft()
            
            self.fock_matrices.append(np.copy(fock_mat.to_numpy(0)))
            self.den_matrices.append(np.copy(den_mat.total_to_numpy(0)))

    def get_effective_fock(self, fock_mat, ovl_mat, oao_mat):
        
        if len(self.fock_matrices) == 1:
            return np.copy(self.fock_matrices[0])
        
        if len(self.fock_matrices) > 1:
            weights = self.comp_c2diis_weigths(ovl_mat, oao_mat)
            return self.get_scaled_fock(weights)
        
        return np.copy(fock_mat.to_numpy(0))
            
    def comp_c2diis_weigths(self, ovl_mat, oao_mat):
        
        evecs = self.comp_error_vectors(ovl_mat, oao_mat)
        
        bmat = self.comp_bmatrix(evecs)
        beigs, bvecs = np.linalg.eigh(bmat);
        
        weights = self.pick_weights(evecs, self.norm_bvectors(bvecs))

        return weights
    
    def get_scaled_fock(self, weights):
        
        effmat = np.zeros(self.fock_matrices[0].shape, dtype=float)
      
        for w, fmat in zip(weights, self.fock_matrices):
            effmat = effmat + w * fmat
        
        return effmat
    
    def comp_error_vectors(self, ovl_mat, oao_mat):
        
        smat = ovl_mat.to_numpy()
        tmat = oao_mat.to_numpy()
        
        evecs = []
        for fmat, dmat in zip(self.fock_matrices, self.den_matrices):
            fa = np.matmul(fmat, np.matmul(dmat, smat))
            fb = np.matmul(smat, np.matmul(dmat, fmat))
            evecs.append(np.matmul(tmat.transpose(), np.matmul(fa - fb, tmat)))
        
        return evecs
    
    def comp_bmatrix(self, evecs):
        
        bdim = len(evecs)
        bmat = np.zeros(shape=(bdim, bdim), dtype=float)
    
        for i in range(bdim):
            for j in range(i, bdim):
                fij = np.vdot(evecs[i], evecs[j])
                bmat[i][j] = fij
                bmat[j][i] = fij
        
        return bmat
    
    def norm_bvectors(self, bvectors):
    
        bsums = np.sum(bvectors, axis=0)
        
        normvecs = []
        for i in range(len(bsums)):
            if abs(bsums[i]) > 1.0e-6:
                normvecs.append(bvectors[:,i] / bsums[i]);

        return normvecs
    
    def pick_weights(self, evecs, weights):
        
        fmin = 1.0e8
        wmin = weights[0]
        
        for w in weights:
            ev = np.zeros(evecs[0].shape, dtype=float)
            for f, v in zip(w, evecs):
                ev = ev + f * v
            fact = np.vdot(ev, ev)
            if fmin > fact:
                fmin = fact
                wmin = w
    
        return wmin
    
    def solve_diis_weigts(self, bmat):
        
        bdim = bmat.shape[0]
        vmat = np.zeros(shape=(bdim), dtype=float)
        vmat[bdim-1] = 1.0
        
        weights = np.linalg.solve(bmat, vmat)
        return weights
    
    def check_diis_weights(self, weights):
        
        wdim = len(weights)
        for i in range(wdim-1):
            if abs(weights[i]) > 1.5:
                self.fock_matrices.popleft()
                self.den_matrices.popleft()
                if len(self.fock_matrices) == 1:
                    return False
                else:
                    return True
    
        return False
    
    def apply_level_shift(self, fock_mat, oao_mat, molecule):
        
        if self.use_level_shift:
            tmat = oao_mat.to_numpy()
            ndomo = molecule.number_of_electrons() // 2
            ndim = fock_mat.shape[0]
            fdisp = np.zeros(shape=(ndim, ndim), dtype=float)
            
            for i in range(ndomo):
                fdisp[i][i] = -self.level_shift;
            
            fock_mat = fock_mat + np.matmul(tmat.transpose(),
                                            np.matmul(fdisp, tmat))

    def gen_molecular_orbitals(self, fock_mat, oao_mat, ostream):
        
        smat = oao_mat.to_numpy()
        
        fmo = np.matmul(smat.transpose(), np.matmul(fock_mat, smat))
        
        eigs, evecs = np.linalg.eigh(fmo)
        orb_coefs = np.matmul(smat, evecs)
        
        return MolecularOrbitals.from_numpy_list([orb_coefs], [eigs],
                                                 molorb.rest)
    
    def gen_new_density(self, mol_orbs, molecule):
        
        return mol_orbs.get_rest_density(molecule)
    
    def get_scf_type(self):
        
        return "Spin-Restricted Hatree-Fock"
