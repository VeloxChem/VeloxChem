from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import T4CScreener
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.fockdriver import FockDriver


class TestFockERI:

    def get_data_h2o(self):

        xyzstr = """
        3
        xyz
        O    0.000000000000        0.000000000000        0.000000000000
        H    0.000000000000        0.740848095288        0.582094932012
        H    0.000000000000       -0.740848095288        0.582094932012
        """
        mol = Molecule.read_xyz_string(xyzstr)
        bas = MolecularBasis.read(mol, 'sto-3g')

        return mol, bas

    def test_h2o_fock_eri_sto3g(self):

        mol, bas = self.get_data_h2o()

        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.eri.npy')
        ref_eri = np.load(npyfile)

        naos = ref_eri.shape[0]
        ithresh = 12

        t4c_drv = T4CScreener()
        t4c_drv.partition(bas, mol, "eri")

        fock_drv = FockDriver()
        eri = fock_drv.compute_eri(t4c_drv, naos, ithresh)

        assert np.max(np.abs(eri - ref_eri)) < 1e-12
