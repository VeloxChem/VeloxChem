from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import GridDriver
from veloxchem.molecule import Molecule


class TestDftGrid:

    def test_dft_grid(self):

        molstr = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """
        mol = Molecule.read_str(molstr, 'au')

        grid_drv = GridDriver()
        grid_drv.set_level(1)
        mol_grid = grid_drv.generate(mol)

        gx = mol_grid.x_to_numpy()
        gy = mol_grid.y_to_numpy()
        gz = mol_grid.z_to_numpy()
        gw = mol_grid.w_to_numpy()

        g_data = np.vstack((gx, gy, gz, gw))

        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dftgrid.npy')
        ref_g_data = np.load(npyfile)

        assert np.max(np.abs(g_data - ref_g_data)) < 1e-10
