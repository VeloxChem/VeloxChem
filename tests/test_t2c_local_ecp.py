from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import LocalECPDriver
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.submatrix import SubMatrix
from veloxchem.veloxchemlib import BaseCorePotential


class TestLocalECPDriver:

    def get_data(self):

        costr = """
            Au   0.000   0.000   0.000
             H   0.200   1.800  -1.100
             H   0.100  -0.900  -1.000
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-svp')

        return mol, bas

    def test_local_ecp_auh2_svp(self):

        mol_auh2, bas_svp = self.get_data()
        
        lpot = BaseCorePotential([4.78982000, 2.39491000],
                                 [30.49008890, 5.17107381],
                                 [2, 2])
                                 
        ecp_drv = LocalECPDriver()
        ecp_mat = ecp_drv.compute(mol_auh2, bas_svp, lpot, 0)
        
        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'auh2.def2svp.au.ecp.only.ul.npy')
        ref_mat = np.load(npyfile)
        
        print(ecp_mat.full_matrix().to_numpy()[5,5])
        
        print(ref_mat[5, 5])
        
        assert False
