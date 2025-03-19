from pathlib import Path
import numpy as np
import h5py

from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.veloxchemlib import DenseMatrix
from veloxchem.veloxchemlib import MolecularGrid
from veloxchem.veloxchemlib import XCIntegrator


class TestGtoValues:

    def test_gto_values_and_derivatives(self):

        mol_string = """
        O   0.0   0.0   0.0
        H   0.0   1.4   1.1
        H   0.0  -1.4   1.1
        """
        basis_label = 'def2-svp'

        mol = Molecule.read_molecule_string(mol_string, units='au')
        bas = MolecularBasis.read(mol, basis_label)

        here = Path(__file__).parent
        hf = h5py.File(str(here / 'data' / 'water_gto_values.h5'), 'r')
        grid_data = np.array(hf.get('grid_data'))
        ref_gto_values = np.array(hf.get('gto_values'))
        ref_gto_values_x = np.array(hf.get('gto_values_x'))
        ref_gto_values_y = np.array(hf.get('gto_values_y'))
        ref_gto_values_z = np.array(hf.get('gto_values_z'))
        hf.close()

        mol_grid_copy = MolecularGrid(DenseMatrix(grid_data))
        mol_grid_copy.partition_grid_points()

        xc_drv = XCIntegrator()
        gto_values = xc_drv.compute_gto_values(mol, bas, mol_grid_copy)
        (gto_values_2, gto_values_x, gto_values_y,
         gto_values_z) = xc_drv.compute_gto_values_and_derivatives(
             mol, bas, mol_grid_copy, 1)

        assert np.max(np.abs(gto_values - ref_gto_values)) < 1.0e-12
        assert np.max(np.abs(gto_values_2 - ref_gto_values)) < 1.0e-12
        assert np.max(np.abs(gto_values_x - ref_gto_values_x)) < 1.0e-12
        assert np.max(np.abs(gto_values_y - ref_gto_values_y)) < 1.0e-12
        assert np.max(np.abs(gto_values_z - ref_gto_values_z)) < 1.0e-12
