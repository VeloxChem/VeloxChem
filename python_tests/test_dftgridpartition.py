import numpy as np
import pytest

from veloxchem.veloxchemlib import (GridDriver, DensityGridDriver, XCIntegrator,
                                    XCNewIntegrator)
from veloxchem.veloxchemlib import is_single_node
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestDftGridPartition:

    def run_dft_grid_partition(self, mol_str, basis_label, xcfun_label,
                               grid_level):

        molecule = Molecule.read_str(mol_str, units='angstrom')
        basis = MolecularBasis.read(molecule, basis_label)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.grid_level = grid_level
        scf_drv.compute(molecule, basis)
        density = scf_drv.density

        grid_drv = GridDriver()
        grid_drv.set_level(grid_level)
        mol_grid = grid_drv.generate(molecule)

        dengrid_drv = DensityGridDriver()
        vxc = mol_grid.test_partition(molecule, basis, density, xcfun_label,
                                      dengrid_drv)

        xc_drv = XCIntegrator()
        vxc_ref = xc_drv.integrate(density, molecule, basis, mol_grid,
                                   xcfun_label)

        max_diff = np.max(
            np.abs(vxc.to_numpy() - vxc_ref.get_matrix().to_numpy()))
        assert max_diff < 1.0e-11

        xc_drv = XCNewIntegrator()
        xc_drv.initialize_grid(mol_grid)
        xc_drv.partition_grid()
        vxc2 = xc_drv.integrate_vxc(molecule, basis, density, xcfun_label)

        max_diff2 = np.max(
            np.abs(vxc2.to_numpy() - vxc_ref.get_matrix().to_numpy()))
        assert max_diff2 < 1.0e-11

    @pytest.mark.skipif(not is_single_node(), reason='single node only')
    def test_slater(self):

        mol_str = """
            O  0.0000000000   0.0000000000  -0.0254395383
            H  0.0000000000   0.7695699584   0.5948147012
            H  0.0000000000  -0.7695699584   0.5948147012
            O  4.0000000000   0.0000000000  -0.0254395383
            H  4.0000000000   0.7695699584   0.5948147012
            H  4.0000000000  -0.7695699584   0.5948147012
        """
        basis_label = 'def2-svp'
        xcfun_label = 'slater'
        grid_level = 4

        self.run_dft_grid_partition(mol_str, basis_label, xcfun_label,
                                    grid_level)
