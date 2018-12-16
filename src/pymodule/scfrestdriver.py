from .veloxchemlib import mpi_master

from .veloxchemlib import MolecularOrbitals

from .veloxchemlib import molorb

from .aodensitymatrix import AODensityMatrix

from .scfdriver import ScfDriver

from .c2diis import CTwoDiis

import numpy as np

class ScfRestrictedDriver(ScfDriver):

    def __init__(self):
        
        super().__init__()
    
    def comp_energy(self, fock_mat, kin_mat, npot_mat, den_mat, comm):
        
        if self.rank == mpi_master():
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
        else:
            e_ee = 0.0
            e_kin = 0.0
            e_en = 0.0
        
        e_ee = comm.bcast(e_ee, mpi_master())
        e_kin = comm.bcast(e_kin, mpi_master())
        e_en = comm.bcast(e_en, mpi_master())
        
        return (e_ee, e_kin, e_en)
    
    def comp_full_fock(self, fock_mat, kin_mat, npot_mat):

        if self.rank == mpi_master():
            fock_mat.add_hcore(kin_mat, npot_mat, 0)

    def comp_gradient(self, fock_mat, ovl_mat, den_mat, comm):
        
        if self.rank == mpi_master():
            smat = ovl_mat.to_numpy()
            dmat = den_mat.total_to_numpy(0)
            fmat = fock_mat.to_numpy(0)
        
            fa = np.matmul(fmat, np.matmul(dmat, smat))
            fb = np.matmul(smat, np.matmul(dmat, fmat))
        
            e_grad = 2.0 * np.linalg.norm(fa - fb)
        else:
            e_grad = 0.0

        e_grad = comm.bcast(e_grad, mpi_master())

        return e_grad

    def comp_density_change(self, den_mat, old_den_mat, comm):
        
        if self.rank == mpi_master():
            diff_mat = den_mat.sub(old_den_mat)
            ddmat = diff_mat.total_to_numpy(0)
        
            diff_den = np.linalg.norm(ddmat)
        else:
            diff_den = 0.0
        
        diff_den = comm.bcast(diff_den, mpi_master())
        
        return diff_den
    
    def store_diis_data(self, i, fock_mat, den_mat):
        
        if self.rank == mpi_master():
            
            if not self.skip_iter:
            
                if len(self.fock_matrices) == self.max_err_vecs:
                
                    self.fock_matrices.popleft()
                    self.den_matrices.popleft()
            
                self.fock_matrices.append(np.copy(fock_mat.to_numpy(0)))
                self.den_matrices.append(np.copy(den_mat.total_to_numpy(0)))

    def get_effective_fock(self, fock_mat, ovl_mat, oao_mat):
        
        if self.rank == mpi_master():

            if len(self.fock_matrices) == 1:
            
                return np.copy(self.fock_matrices[0])
        
            if len(self.fock_matrices) > 1:
            
                acc_diis = CTwoDiis()
                
                acc_diis.compute_error_vectors(self.fock_matrices,
                                               self.den_matrices,
                                               ovl_mat, oao_mat)
            
                weights = acc_diis.compute_weights()
            
                return self.get_scaled_fock(weights)
        
            return np.copy(fock_mat.to_numpy(0))

        return None
    
    def get_scaled_fock(self, weights):
        
        effmat = np.zeros(self.fock_matrices[0].shape, dtype=float)
      
        for w, fmat in zip(weights, self.fock_matrices):
            
            effmat = effmat + w * fmat
        
        return effmat
    
    def gen_molecular_orbitals(self, fock_mat, oao_mat, ostream):
        
        if self.rank == mpi_master():
            
            tmat = oao_mat.to_numpy()
        
            fmo = np.matmul(tmat.transpose(), np.matmul(fock_mat, tmat))
        
            eigs, evecs = np.linalg.eigh(fmo)
            orb_coefs = np.matmul(tmat, evecs)
        
            return MolecularOrbitals.from_numpy_list([orb_coefs], [eigs],
                                                     molorb.rest)
        
        return MolecularOrbitals()
    
    def gen_new_density(self, molecule):
        
        if self.rank == mpi_master():
            
            return self.mol_orbs.get_rest_density(molecule)

        return AODensityMatrix()
    
    def get_scf_type(self):
        
        return "Spin-Restricted Hatree-Fock"
