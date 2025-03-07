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

class TestRIJFockGeomGradDriver:

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
    
    def get_data_h2o_trimer(self):

        h2ostr = """
            O    1.206898105612   1.227024188072  -0.113458607566
            H    0.490554304233   1.138628097468   0.508036765414
            H    1.992568601678   1.153431301285   0.406410355109
            O   -0.995720200000   0.016041500000   1.242255600000
            H   -1.454270300000  -0.566974100000   1.847281700000
            H   -0.937795000000  -0.481791200000   0.426756200000
            O   -0.243234300000  -1.019856600000  -1.195380800000
            H    0.436753600000  -0.375943300000  -0.997329700000
            H   -0.503183500000  -0.825149200000  -2.095795900000
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'def2-svp')
        aux_bas = MolecularBasis.read(mol, 'def2-universal-jkfit')

        return mol, bas, aux_bas

    def test_h2o_fock_2j_grad_h2_sto3g(self):

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        density = np.load(npyfile)
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(density)


        fock_drv = FockGeom1000Driver()
        fmats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, 2, "j", 0.0, 0.0)
        fmatx = fmats.matrix('X').full_matrix().to_numpy()
        gradx = 2.0 * np.trace(np.matmul(fmatx, density))
        fmaty = fmats.matrix('Y').full_matrix().to_numpy()
        grady = 2.0 * np.trace(np.matmul(fmaty, density))
        fmatz = fmats.matrix('Z').full_matrix().to_numpy()
        gradz = 2.0 * np.trace(np.matmul(fmatz, density))
        ref_grad = [gradx, grady, gradz]
        
        # compute J metric
        t2c_drv = TwoCenterElectronRepulsionDriver()
        matj = t2c_drv.compute(mol_h2o, bas_aux)
        rmatj = np.linalg.inv(matj.full_matrix().to_numpy())
        invmatj = SubMatrix([0, 0, rmatj.shape[0], rmatj.shape[0]])
        invmatj.set_values(rmatj)
        
        ri_fock_drv = RIFockDriver(invmatj)
        ri_fock_drv.prepare_buffers(mol_h2o, bas_sto3g, bas_aux)
        gv = ri_fock_drv.compute_bq_vector(den_mat)
        
        ri_grad_drv = RIFockGradDriver()
        g = ri_grad_drv.compute(bas_sto3g, bas_aux, mol_h2o, gv, den_mat, 2)
        grad = g.coordinates()
        
        for i in range(3):
            assert mt.isclose(ref_grad[i], grad[i], rel_tol=1.0e-5, abs_tol=1.0e-5)
            
    def test_h2o_fock_2j_grad_h2_sto3g_screened(self):

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        density = np.load(npyfile)
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(density)

        fock_drv = FockGeom1000Driver()
        fmats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, 2, "2jkx", 0.0, 0.0)
        fmatx = fmats.matrix('X').full_matrix().to_numpy()
        gradx = np.trace(np.matmul(fmatx, density))
        fmaty = fmats.matrix('Y').full_matrix().to_numpy()
        grady = np.trace(np.matmul(fmaty, density))
        fmatz = fmats.matrix('Z').full_matrix().to_numpy()
        gradz = np.trace(np.matmul(fmatz, density))
        ref_grad = [gradx, grady, gradz]
        
        # compute J metric
        t2c_drv = TwoCenterElectronRepulsionDriver()
        matj = t2c_drv.compute(mol_h2o, bas_aux)
        rmatj = np.linalg.inv(matj.full_matrix().to_numpy())
        invmatj = SubMatrix([0, 0, rmatj.shape[0], rmatj.shape[0]])
        invmatj.set_values(rmatj)
        
        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_sto3g, mol_h2o, "eri")
        
        ri_fock_drv = RIFockDriver(invmatj)
        ri_fock_drv.prepare_buffers(mol_h2o, bas_sto3g, bas_aux)
        gv = ri_fock_drv.compute_bq_vector(den_mat)
        
        ri_grad_drv = RIFockGradDriver()
        g = ri_grad_drv.compute(t4c_drv, bas_sto3g, bas_aux, mol_h2o, gv, den_mat, 2, 12)
        grad = g.coordinates()
        
        for i in range(3):
            assert mt.isclose(ref_grad[i], grad[i], rel_tol=1.0e-5, abs_tol=1.0e-5)
            
    def test_h2o_fock_2j_grad_h2_sto3g_screened_direct(self):

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        density = np.load(npyfile)
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(density)

        fock_drv = FockGeom1000Driver()
        fmats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, 2, "2jkx", 0.0, 0.0)
        fmatx = fmats.matrix('X').full_matrix().to_numpy()
        gradx = np.trace(np.matmul(fmatx, density))
        fmaty = fmats.matrix('Y').full_matrix().to_numpy()
        grady = np.trace(np.matmul(fmaty, density))
        fmatz = fmats.matrix('Z').full_matrix().to_numpy()
        gradz = np.trace(np.matmul(fmatz, density))
        ref_grad = [gradx, grady, gradz]
        
        # compute J metric
        t2c_drv = TwoCenterElectronRepulsionDriver()
        matj = t2c_drv.compute(mol_h2o, bas_aux)
        rmatj = np.linalg.inv(matj.full_matrix().to_numpy())
        invmatj = SubMatrix([0, 0, rmatj.shape[0], rmatj.shape[0]])
        invmatj.set_values(rmatj)
        
        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_sto3g, mol_h2o, "eri")
        
        ri_fock_drv = RIFockDriver(invmatj)
        ri_fock_drv.prepare_buffers(mol_h2o, bas_sto3g, bas_aux)
        gv = ri_fock_drv.compute_bq_vector(den_mat)
        
        ri_grad_drv = RIFockGradDriver()
        g = ri_grad_drv.direct_compute(t4c_drv, bas_sto3g, bas_aux, mol_h2o, gv, den_mat, 2, 12)
        grad = g.coordinates()
        
        for i in range(3):
            assert mt.isclose(ref_grad[i], grad[i], rel_tol=1.0e-5, abs_tol=1.0e-5)
            
    def test_h2o_fock_2j_grad_h2o_sto3g(self):

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        density = np.load(npyfile)
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(density)

        # convectional gradient
        fock_drv = FockGeom1000Driver()
        ref_grad = []
        
        for i in range(3):
            fmats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, i, "2jkx", 0.0, 0.0)
            fmatx = fmats.matrix('X').full_matrix().to_numpy()
            gradx = np.trace(np.matmul(fmatx, density))
            fmaty = fmats.matrix('Y').full_matrix().to_numpy()
            grady = np.trace(np.matmul(fmaty, density))
            fmatz = fmats.matrix('Z').full_matrix().to_numpy()
            gradz = np.trace(np.matmul(fmatz, density))
            ref_grad.append(gradx)
            ref_grad.append(grady)
            ref_grad.append(gradz)
        
        # compute J metric
        t2c_drv = TwoCenterElectronRepulsionDriver()
        matj = t2c_drv.compute(mol_h2o, bas_aux)
        rmatj = np.linalg.inv(matj.full_matrix().to_numpy())
        invmatj = SubMatrix([0, 0, rmatj.shape[0], rmatj.shape[0]])
        invmatj.set_values(rmatj)
        
        # compute B_q vector
        ri_fock_drv = RIFockDriver(invmatj)
        ri_fock_drv.prepare_buffers(mol_h2o, bas_sto3g, bas_aux)
        gv = ri_fock_drv.compute_bq_vector(den_mat)
        
        # compute RI gradients for all atoms
        ri_grad_drv = RIFockGradDriver()
        g = ri_grad_drv.compute(bas_sto3g, bas_aux, mol_h2o, gv, den_mat, [0, 1, 2])
       
        # check results
        for i in range(3):
            grad = g[i].coordinates()
            for j in range(3):
                assert mt.isclose(ref_grad[3 * i + j], grad[j], rel_tol=1.0e-5, abs_tol=1.0e-5)
                
    def test_h2o_fock_2j_grad_h2o_sto3g_screened(self):

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        density = np.load(npyfile)
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(density)

        # convectional gradient
        fock_drv = FockGeom1000Driver()
        ref_grad = []
        
        for i in range(3):
            fmats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, i, "2jkx", 0.0, 0.0)
            fmatx = fmats.matrix('X').full_matrix().to_numpy()
            gradx = np.trace(np.matmul(fmatx, density))
            fmaty = fmats.matrix('Y').full_matrix().to_numpy()
            grady = np.trace(np.matmul(fmaty, density))
            fmatz = fmats.matrix('Z').full_matrix().to_numpy()
            gradz = np.trace(np.matmul(fmatz, density))
            ref_grad.append(gradx)
            ref_grad.append(grady)
            ref_grad.append(gradz)
        
        # compute J metric
        t2c_drv = TwoCenterElectronRepulsionDriver()
        matj = t2c_drv.compute(mol_h2o, bas_aux)
        rmatj = np.linalg.inv(matj.full_matrix().to_numpy())
        invmatj = SubMatrix([0, 0, rmatj.shape[0], rmatj.shape[0]])
        invmatj.set_values(rmatj)
        
        # compute B_q vector
        ri_fock_drv = RIFockDriver(invmatj)
        ri_fock_drv.prepare_buffers(mol_h2o, bas_sto3g, bas_aux)
        gv = ri_fock_drv.compute_bq_vector(den_mat)
        
        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_sto3g, mol_h2o, "eri")
        
        # compute RI gradients for all atoms
        ri_grad_drv = RIFockGradDriver()
        g = ri_grad_drv.compute(t4c_drv, bas_sto3g, bas_aux, mol_h2o, gv, den_mat, [0, 1, 2], 12)
       
        # check results
        for i in range(3):
            grad = g[i].coordinates()
            for j in range(3):
                assert mt.isclose(ref_grad[3 * i + j], grad[j], rel_tol=1.0e-5, abs_tol=1.0e-5)
                
    def test_h2o_trimer_fock_2j_grad_h2o_svp_screened(self):

        mol_h2o_trimer, bas_svp, bas_aux = self.get_data_h2o_trimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.trimer.svp.density.npy')
        density = np.load(npyfile)
        den_mat = make_matrix(bas_svp, mat_t.symmetric)
        den_mat.set_values(density)

        # convectional gradient
        fock_drv = FockGeom1000Driver()
        ref_grad = []
        
        for i in range(9):
            fmats = fock_drv.compute(bas_svp, mol_h2o_trimer, den_mat, i, "2jkx", 0.0, 0.0)
            fmatx = fmats.matrix('X').full_matrix().to_numpy()
            gradx = np.trace(np.matmul(fmatx, density))
            fmaty = fmats.matrix('Y').full_matrix().to_numpy()
            grady = np.trace(np.matmul(fmaty, density))
            fmatz = fmats.matrix('Z').full_matrix().to_numpy()
            gradz = np.trace(np.matmul(fmatz, density))
            ref_grad.append(gradx)
            ref_grad.append(grady)
            ref_grad.append(gradz)
        
        # compute J metric
        t2c_drv = TwoCenterElectronRepulsionDriver()
        matj = t2c_drv.compute(mol_h2o_trimer, bas_aux)
        rmatj = np.linalg.inv(matj.full_matrix().to_numpy())
        invmatj = SubMatrix([0, 0, rmatj.shape[0], rmatj.shape[0]])
        invmatj.set_values(rmatj)
        
        # compute B_q vector
        ri_fock_drv = RIFockDriver(invmatj)
        ri_fock_drv.prepare_buffers(mol_h2o_trimer, bas_svp, bas_aux)
        gv = ri_fock_drv.compute_bq_vector(den_mat)
        
        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_svp, mol_h2o_trimer, "eri")
        
        # compute RI gradients for all atoms
        ri_grad_drv = RIFockGradDriver()
        for i in range(9):
            g = ri_grad_drv.compute(t4c_drv, bas_svp, bas_aux, mol_h2o_trimer, gv, den_mat, i, 12)
            grad = g.coordinates()
            for j in range(3):
                assert mt.isclose(ref_grad[3 * i + j], grad[j], rel_tol=1.0e-5, abs_tol=1.0e-5)
       
    def test_h2o_mpi_fock_2j_grad_h2o_sto3g(self):

        comm = MPI.COMM_WORLD

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()
        
        # distribute density matrix
                    
        den_mat = None
        if comm.Get_rank() == 0:
            here = Path(__file__).parent
            npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
            density = np.load(npyfile)
            den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
            den_mat.set_values(density)
        den_mat = comm.bcast(den_mat, 0)
        
        # compute RI gradients for all atoms
        ri_grad_drv = RIFockGradDriver()
        g = ri_grad_drv.mpi_compute(comm, mol_h2o, bas_sto3g, bas_aux, den_mat)
        
        if comm.Get_rank() == 0:
            fock_drv = FockGeom1000Driver()
            ref_grad = []
        
            for i in range(3):
                fmats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, i, "2jkx", 0.0, 0.0)
                fmatx = fmats.matrix('X').full_matrix().to_numpy()
                gradx = np.trace(np.matmul(fmatx, density))
                fmaty = fmats.matrix('Y').full_matrix().to_numpy()
                grady = np.trace(np.matmul(fmaty, density))
                fmatz = fmats.matrix('Z').full_matrix().to_numpy()
                gradz = np.trace(np.matmul(fmatz, density))
                ref_grad.append(gradx)
                ref_grad.append(grady)
                ref_grad.append(gradz)
                
            # check results
            for i in range(3):
                for j in range(3):
                    assert mt.isclose(ref_grad[3 * i + j], g[i, j], rel_tol=1.0e-5, abs_tol=1.0e-5)
                    
    def test_h2o_mpi_fock_2j_grad_h2o_sto3g(self):

        comm = MPI.COMM_WORLD

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()
        
        # distribute density matrix
                    
        den_mat = None
        if comm.Get_rank() == 0:
            here = Path(__file__).parent
            npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
            density = np.load(npyfile)
            den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
            den_mat.set_values(density)
        den_mat = comm.bcast(den_mat, 0)
        
        # compute RI gradients for all atoms
        ri_grad_drv = RIFockGradDriver()
        g = ri_grad_drv.mpi_compute(comm, mol_h2o, bas_sto3g, bas_aux, den_mat)
        
        if comm.Get_rank() == 0:
            fock_drv = FockGeom1000Driver()
            ref_grad = []
        
            for i in range(3):
                fmats = fock_drv.compute(bas_sto3g, mol_h2o, den_mat, i, "2jkx", 0.0, 0.0)
                fmatx = fmats.matrix('X').full_matrix().to_numpy()
                gradx = np.trace(np.matmul(fmatx, density))
                fmaty = fmats.matrix('Y').full_matrix().to_numpy()
                grady = np.trace(np.matmul(fmaty, density))
                fmatz = fmats.matrix('Z').full_matrix().to_numpy()
                gradz = np.trace(np.matmul(fmatz, density))
                ref_grad.append(gradx)
                ref_grad.append(grady)
                ref_grad.append(gradz)
                
            # check results
            for i in range(3):
                for j in range(3):
                    assert mt.isclose(ref_grad[3 * i + j], g[i, j], rel_tol=1.0e-5, abs_tol=1.0e-5)
