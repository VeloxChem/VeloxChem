from pathlib import Path
import math as mt
import numpy as np

from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import FockGeom1000Driver
from veloxchem import RIFockDriver
from veloxchem import RIFockGradDriver
from veloxchem import TwoCenterElectronRepulsionDriver
from veloxchem import SubMatrix
from veloxchem import Matrix
from veloxchem import Matrices
from veloxchem import make_matrix
from veloxchem import mat_t


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

    def test_h2o_fock_2j_grad_o1_sto3g(self):

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
