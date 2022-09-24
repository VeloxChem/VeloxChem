import numpy as np

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import (GridDriver, XCMolecularGradient,
                                    XCNewMolecularGradient)
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestXCMolGrad:

    def run_xc_mol_grad(self, molecule, basis, xcfun, ref_grad):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.state = False

        scf_drv.dft = True
        scf_drv.xcfun = xcfun
        scf_drv.conv_thresh = 1.0e-8

        scf_drv.compute(molecule, basis)
        density = scf_drv.density

        grid_drv = GridDriver(scf_drv.comm)
        grid_drv.set_level(scf_drv.grid_level)
        mol_grid = grid_drv.generate(molecule)
        mol_grid.distribute(scf_drv.rank, scf_drv.nodes, scf_drv.comm)

        grad_drv = XCMolecularGradient(scf_drv.comm)
        mol_grad = grad_drv.integrate_vxc_gradient(density, molecule, basis,
                                                   mol_grid, xcfun)
        mol_grad = scf_drv.comm.reduce(mol_grad, root=mpi_master())

        mol_grid_new = grid_drv.generate(molecule)
        mol_grid_new.partition_grid_points()
        mol_grid_new.distribute_counts_and_displacements(
            scf_drv.rank, scf_drv.nodes, scf_drv.comm)

        grad_drv_new = XCNewMolecularGradient(scf_drv.comm)
        mol_grad_new = grad_drv_new.integrate_vxc_gradient(
            molecule, basis, density, mol_grid_new, xcfun)
        mol_grad_new = scf_drv.comm.reduce(mol_grad_new, root=mpi_master())

        if scf_drv.rank == mpi_master():
            assert np.max(np.abs(mol_grad - ref_grad)) < 1.0e-4
            assert np.max(np.abs(mol_grad_new - ref_grad)) < 1.0e-4

    def test_xc_mol_grad_slater(self):

        mol_str = """
            O  0.0000000000   0.0000000000  -0.0254395383
            H  0.0000000000   0.7695699584   0.5948147012
            H  0.0000000000  -0.7695699584   0.5948147012
        """

        molecule = Molecule.read_str(mol_str, units='angstrom')
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        xcfun = 'slater'

        ref_grad = np.array([
            [0., 0., -0.4075911],
            [-0., 0.2540621, 0.2037956],
            [0., -0.2540621, 0.2037956],
        ])

        self.run_xc_mol_grad(molecule, basis, xcfun, ref_grad)

    def test_xc_mol_grad_blyp(self):

        mol_str = """
            O  0.0000000000   0.0000000000  -0.0254395383
            H  0.0000000000   0.7695699584   0.5948147012
            H  0.0000000000  -0.7695699584   0.5948147012
        """

        molecule = Molecule.read_str(mol_str, units='angstrom')
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        xcfun = 'blyp'

        ref_grad = np.array([
            [-0., 0.0000000, -0.4453751],
            [0., 0.2801306, 0.2226875],
            [0., -0.2801306, 0.2226875],
        ])

        self.run_xc_mol_grad(molecule, basis, xcfun, ref_grad)
