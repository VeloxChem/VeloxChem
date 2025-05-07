from pathlib import Path
import math as mt
import numpy as np
import pytest

from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import FockGeom1000Driver
from veloxchem import T4CScreener
from veloxchem import make_matrix, mat_t
from veloxchem.rifockdriver import RIFockDriver
from veloxchem.veloxchemlib import RIFockGradDriver


@pytest.mark.solvers
class TestRIJFockGeomGradDriver:

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
            fmats = fock_drv.compute(bas_svp, mol_h2o_trimer, den_mat, i,
                                     "2jkx", 0.0, 0.0)
            fmatx = fmats.matrix('X').full_matrix().to_numpy()
            gradx = np.trace(np.matmul(fmatx, density))
            fmaty = fmats.matrix('Y').full_matrix().to_numpy()
            grady = np.trace(np.matmul(fmaty, density))
            fmatz = fmats.matrix('Z').full_matrix().to_numpy()
            gradz = np.trace(np.matmul(fmatz, density))
            ref_grad.append(gradx)
            ref_grad.append(grady)
            ref_grad.append(gradz)

        # compute B_q vector
        ri_fock_drv = RIFockDriver()
        ri_fock_drv.prepare_buffers(mol_h2o_trimer,
                                    bas_svp,
                                    bas_aux,
                                    verbose=False)
        gv = ri_fock_drv.compute_bq_vector(den_mat)

        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_svp, mol_h2o_trimer, "eri")

        # compute RI gradients for all atoms
        ri_grad_drv = RIFockGradDriver()
        for i in range(9):
            g = ri_grad_drv.compute(t4c_drv, bas_svp, bas_aux, mol_h2o_trimer,
                                    gv, den_mat, i, 12)
            grad = g.coordinates()
            for j in range(3):
                assert mt.isclose(ref_grad[3 * i + j],
                                  grad[j],
                                  rel_tol=1.0e-5,
                                  abs_tol=1.0e-5)
