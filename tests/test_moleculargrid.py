from pathlib import Path
import numpy as np

from veloxchem.griddriver import GridDriver
from veloxchem.molecule import Molecule


class TestMolecularGrid:

    def test_molecular_grid(self):

        molstr = """N  -3.710   3.019  -0.037
                    H  -3.702   4.942   0.059
                    H  -4.704   2.415   1.497
                    H  -4.780   2.569  -1.573"""

        mol = Molecule.read_molecule_string(molstr, units='au')

        grid_drv = GridDriver()
        grid_drv.set_level(1)

        num_gpus_per_node = 1
        mol_grid = grid_drv.generate(mol, num_gpus_per_node)

        grid_data = mol_grid.grid_to_numpy()

        here = Path(__file__).parent
        npy_file = str(here / 'data' / 'nh3_grid.npy')
        ref_grid_data = np.load(npy_file)

        assert np.max(np.abs(grid_data - ref_grid_data)) < 1.0e-10
