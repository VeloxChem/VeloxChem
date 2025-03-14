from mpi4py import MPI
from pathlib import Path
import math as mt
import numpy as np

from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import FockGeom1000Driver
from veloxchem import RIFockDriver
from veloxchem import TwoCenterElectronRepulsionDriver
from veloxchem import T4CScreener
from veloxchem import SubMatrix
from veloxchem import Matrix
from veloxchem import Matrices
from veloxchem import make_matrix
from veloxchem import mat_t
from veloxchem.rigradientdriver import RIFockGradDriver

class TestRIJFockGeomExGradDriver:

    def get_data_h2o(self):

        h2ostr = """
            O    0.000000000000        0.000000000000        0.000000000000
            H    0.000000000000        0.740848095288        0.582094932012
            H    0.000000000000       -0.740848095288        0.582094932012
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'sto-3g')
        aux_bas = MolecularBasis.read(mol, 'def2-universal-jkfit')

        return mol, bas, aux_bas
    
    def test_h2o_fock_2j_grad_h2_sto3g_screened(self):

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        density = np.load(npyfile)
        gs_den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        gs_den_mat.set_values(density)
        dims = density.shape
        density = np.random.rand(dims[0], dims[1])
        density = density + density.T
        rw_den_sym_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        rw_den_sym_mat.set_values(0.68 * density)
        rw_den_gen_mat = make_matrix(bas_sto3g, mat_t.general)
        rw_den_gen_mat.set_values(0.68 * density)
        
        ket_t4c = T4CScreener()
        ket_t4c.partition(bas_sto3g, mol_h2o, 'eri')
        
        fock_grad_drv = FockGeom1000Driver()
        
        ref_grads = []
        for i in range(3):
            bra_t4c = T4CScreener()
            bra_t4c.partition_atom(bas_sto3g, mol_h2o, 'eri', i)
            grad = fock_grad_drv.compute(bas_sto3g, bra_t4c, ket_t4c,
                                         gs_den_mat, rw_den_gen_mat,
                                         i, 'j', 0.0, 0.0, 15)
            ref_grads.append([grad[0], grad[1], grad[2]])
        
        print(ref_grads)
                
        # compute J metric
        t2c_drv = TwoCenterElectronRepulsionDriver()
        matj = t2c_drv.compute(mol_h2o, bas_aux)
        rmatj = np.linalg.inv(matj.full_matrix().to_numpy())
        invmatj = SubMatrix([0, 0, rmatj.shape[0], rmatj.shape[0]])
        invmatj.set_values(rmatj)
        
        ri_fock_drv = RIFockDriver(invmatj)
        ri_fock_drv.prepare_buffers(mol_h2o, bas_sto3g, bas_aux)
        gs_bq = ri_fock_drv.compute_bq_vector(gs_den_mat)
        rw_bq = ri_fock_drv.compute_bq_vector(rw_den_sym_mat)
        
        ri_grad_drv = RIFockGradDriver()
        
        ri_grads = []
        for i in range(3):
            grad = ri_grad_drv.direct_compute(ket_t4c, bas_sto3g, bas_aux, mol_h2o, rw_bq, gs_bq, rw_den_sym_mat, gs_den_mat, i, 15)
            ri_grads.append(grad.coordinates())
        
        print(ri_grads)
    
        for i in range(3):
            for j in range(3):
                assert mt.isclose(ref_grads[i][j], ri_grads[i][j], rel_tol=1.0e-4, abs_tol=1.0e-4)
        
        assert False
        

       
            
    
