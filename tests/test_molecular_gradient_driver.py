from pathlib import Path
import numpy as np

from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import MolecularGradientDriver
from veloxchem import SubMatrix
from veloxchem import Matrix
from veloxchem import Matrices
from veloxchem import make_matrix
from veloxchem import mat_t


class TestMolecularGradientDriver:

    def get_data_h2o(self):

        h2ostr = """
            O    0.000000000000        0.000000000000        0.000000000000
            H    0.000000000000        0.740848095288        0.582094932012
            H    0.000000000000       -0.740848095288        0.582094932012
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'sto-3g')

        return mol, bas

    def test_h2o_hf(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))
        
        grad_drv = MolecularGradientDriver()
        mol_grad = grad_drv.compute(mol_h2o, bas_sto3g, den_mat)
        
        print(mol_grad)
        assert False

        
        
