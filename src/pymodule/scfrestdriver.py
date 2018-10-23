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
            weights = self.comp_diis_weigths(ovl_mat, oao_mat)
            if len(self.fock_matrices) > 1:
                return self.get_scaled_fock(weights)
        return np.copy(fock_mat.to_numpy(0))
            
    def comp_diis_weigths(self, ovl_mat, oao_mat):
        need_weights = True
        while need_weights:
            bmat = self.comp_diis_matrix(ovl_mat, oao_mat)
            weights = self.solve_diis_weigts(bmat)
            need_weights = self.check_diis_weights(weights)
        return weights
    
    def comp_diis_matrix(self, ovl_mat, oao_mat):
        bdim = len(self.fock_matrices)
        bmat = np.zeros(shape=(bdim+1, bdim+1), dtype=float)
        smat = ovl_mat.to_numpy()
        tmat = oao_mat.to_numpy()
        for i in range(bdim):
            fmat = self.fock_matrices[i]
            dmat = self.den_matrices[i]
            fa = np.matmul(fmat, np.matmul(dmat, smat))
            fb = np.matmul(smat, np.matmul(dmat, fmat))
            evec_i = np.matmul(tmat.transpose(), np.matmul(fa - fb, tmat))
            bmat[i][i] = np.vdot(evec_i, evec_i)
            for j in range(i + 1, bdim):
                fmat = self.fock_matrices[j]
                dmat = self.den_matrices[j]
                fa = np.matmul(fmat, np.matmul(dmat, smat))
                fb = np.matmul(smat, np.matmul(dmat, fmat))
                evec_j = np.matmul(tmat.transpose(), np.matmul(fa - fb, tmat))
                feij = np.vdot(evec_i, evec_j)
                bmat[i][j] = feij
                bmat[j][i] = feij
            bmat[bdim, i] = 1.0
            bmat[i, bdim] = 1.0
        return bmat
    
    def solve_diis_weigts(self, bmat):
        bdim = bmat.shape[0]
        vmat = np.zeros(shape=(bdim), dtype=float)
        vmat[bdim-1] = 1.0
        weights = np.linalg.solve(bmat, vmat)
        return weights
    
    def check_diis_weights(self, weights):
        wdim = len(weights)
        for i in range(len(weights) - 1):
            if abs(weights[i]) > 1.3:
                self.fock_matrices.popleft()
                self.den_matrices.popleft()
                if len(self.fock_matrices) == 1:
                    return False
                else:
                    return True
        return False
    
    def get_scaled_fock(self, weights):
        wdim = len(weights) - 1
        fdim = len(self.fock_matrices)
        if len(self.fock_matrices) != wdim:
            return np.copy(self.fock_matrices[-1])
        fmat = np.copy(self.fock_matrices[0])
        fmat = weights[0] * fmat
        for i in range(1, wdim):
            fmat = fmat + weights[i] * self.fock_matrices[i]
        return fmat
    
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
