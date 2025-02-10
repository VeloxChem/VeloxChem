from pathlib import Path
import numpy as np
import math as mt

from veloxchem import ThreeCenterElectronRepulsionGeom100Driver
from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import SubMatrix
from veloxchem import T3FlatBuffer


class TestThreeCenterElectronRepulsionGeom100Driver:

    def get_data_svp(self):

        h2ostr = """
            O    0.000000000000        0.000000000000        0.000000000000
            H    0.000000000000        0.740848095288        0.582094932012
            H    0.000000000000       -0.740848095288        0.582094932012
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'def2-svp')

        return mol, bas

    def test_electron_repulsion_geom_100_h2o_svp(self):

        mol, bas = self.get_data_svp()

        # compute electron repulsion matrix gradient on auxilary basis
        eri_drv = ThreeCenterElectronRepulsionGeom100Driver()
        eri_buf = eri_drv.compute(mol, bas, bas, 0)
    
        assert False
