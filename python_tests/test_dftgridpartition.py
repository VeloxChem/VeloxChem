import numpy as np
import pytest

from veloxchem.veloxchemlib import (GridDriver, XCIntegrator, XCNewIntegrator)
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

        xc_drv = XCIntegrator()
        vxc_ref = xc_drv.integrate(density, molecule, basis, mol_grid,
                                   xcfun_label)
        vxc_ref_mat = vxc_ref.get_matrix()

        xc_drv = XCNewIntegrator()
        mol_grid.partition_grid_points()
        vxc = xc_drv.integrate_vxc_fock(molecule, basis, density, mol_grid,
                                        xcfun_label)
        vxc_mat = vxc.get_matrix()

        max_diff2 = np.max(np.abs(vxc_mat.to_numpy() - vxc_ref_mat.to_numpy()))
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
