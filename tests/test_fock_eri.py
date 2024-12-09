from pathlib import Path
import numpy as np

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

        fock_drv = FockDriver()
        eri = fock_drv.compute_eri(mol, bas)

        assert np.max(np.abs(eri - ref_eri)) < 1e-12
