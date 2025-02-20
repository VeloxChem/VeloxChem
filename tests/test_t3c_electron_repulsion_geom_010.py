from pathlib import Path
import numpy as np
import math as mt

from veloxchem import ThreeCenterElectronRepulsionGeom010Driver
from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import SubMatrix
from veloxchem import T3FlatBuffer


class TestThreeCenterElectronRepulsionGeom010Driver:

    def get_data_svp(self):

        h2ostr = """
            O    0.000000000000        0.000000000000        0.000000000000
            H    0.000000000000        0.740848095288        0.582094932012
            H    0.000000000000       -0.740848095288        0.582094932012
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'def2-svp')

        return mol, bas

    def test_electron_repulsion_geom_010_h2o_svp(self):

        mol, bas = self.get_data_svp()

        # compute electron repulsion matrix gradient on auxilary basis
        eri_drv = ThreeCenterElectronRepulsionGeom010Driver()
        eri_buf = eri_drv.compute(mol, bas, bas, 0)
        
        # load reference kinetic energy data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.svp.int3c2e.geom.010.o1.npy')
        ref_buf = -np.load(npyfile)
        
        mask_ids = eri_buf.mask_indices()
        
        print(mask_ids)
        
        for i in range(3):
            for j in range(24):
                for k in [0, 1, 2, 7, 8, 11, 12, 15, 16, 19, 20, 21, 22, 23]:
                    for l in range(24):
                        print("grad ", i, " (", j, "|", k, ",", l, ") = ",
                              eri_buf.value(24 * i + j, k, l), " vs. ",
                              ref_buf[i, k, l, j])
                        assert mt.isclose(eri_buf.value(24 * i + j, k, l),
                                          ref_buf[i, k, l, j],
                                          rel_tol=1.0e-12,
                                          abs_tol=1.0e-12)

