from pathlib import Path
import numpy as np
import math as mt

from veloxchem import ThreeCenterElectronRepulsionDriver
from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import SubMatrix
from veloxchem import T3FlatBuffer

class TestThreeCenterElectronRepulsionDriver:

    def get_data_svp(self):

        h2ostr = """
            O    0.000000000000        0.000000000000        0.000000000000
            H    0.000000000000        0.740848095288        0.582094932012
            H    0.000000000000       -0.740848095288        0.582094932012
        """
        
        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'def2-svp')

        return mol, bas
    
    def get_data_qzvp(self):

        h2ostr = """
            O    0.000000000000        0.000000000000        0.000000000000
            H    0.000000000000        0.740848095288        0.582094932012
            H    0.000000000000       -0.740848095288        0.582094932012
        """
        
        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'def2-qzvp')

        return mol, bas

    def test_electron_repulsion_h2o_svp(self):

        mol, bas = self.get_data_svp()

        # compute electron repulsion matrix
        eri_drv = ThreeCenterElectronRepulsionDriver()
        eri_buf = eri_drv.compute(mol, bas, bas)

        # load reference kinetic energy data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.svp.int3c2e.npy')
        ref_buf = np.load(npyfile)
        
        indexes = np.triu_indices(24)
        
        for i in range(24):
            for k, l in zip(indexes[0], indexes[1]):
                print('(',i,"|",k,",",l,") = ", eri_buf.value(i,k,l)," vs. ", ref_buf[k,l,i])
                assert mt.isclose(eri_buf.value(i,k,l), ref_buf[k,l,i], rel_tol=1.0e-12, abs_tol=1.0e-12)
