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
        
        # load reference kinetic energy data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.svp.int3c2e.geom.001.o1.npy')
        ref_buf = -np.load(npyfile)
        
        indexes = np.triu_indices(24)
        bra_ids = eri_buf.indices()
        
        for i in range(3):
            for j in range(3):
                for k, l in zip(indexes[0], indexes[1]):
                    assert mt.isclose(eri_buf.value(14 * i + j, k, l),
                                      ref_buf[i, k, l, bra_ids[j]],
                                      rel_tol=1.0e-12,
                                      abs_tol=1.0e-12)

