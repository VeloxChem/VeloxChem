import numpy as np

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import GridDriver, XCMolecularGradient
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestXCMolGrad:

    def run_xc_mol_grad(self, molecule, basis, xcfun, ref_grad):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.state = False

        scf_drv.dft = True
        scf_drv.xcfun = xcfun
        scf_drv.grid_level = 5
        scf_drv.conv_thresh = 1.0e-8

        scf_drv.compute(molecule, basis)
        density = scf_drv.density

        grid_drv = GridDriver(scf_drv.comm)
        grid_drv.set_level(scf_drv.grid_level)
        mol_grid = grid_drv.generate(molecule)
        mol_grid.distribute(scf_drv.rank, scf_drv.nodes, scf_drv.comm)

        grad_drv = XCMolecularGradient(scf_drv.comm)
        atom_ids = list(range(molecule.number_of_atoms()))
        mol_grad = grad_drv.integrate(atom_ids, density, molecule, basis,
                                      mol_grid, xcfun)
        mol_grad = scf_drv.comm.reduce(mol_grad, root=mpi_master())

        if scf_drv.rank == mpi_master():
            assert np.max(np.abs(mol_grad - ref_grad)) < 1.0e-4

    def test_xc_mol_grad(self):

        mol_str = """
            O  0.0000000000   0.0000000000  -0.0254395383
            H  0.0000000000   0.7695699584   0.5948147012
            H  0.0000000000  -0.7695699584   0.5948147012
        """

        molecule = Molecule.read_str(mol_str, units='angstrom')
        basis = MolecularBasis.read(molecule, 'def2-svp')
        xcfun = 'blyp'

        ref_grad = np.array([
            [-0., 0.0000000, -0.4453755],
            [0., 0.2801308, 0.2226875],
            [0., -0.2801308, 0.2226875],
        ])

        self.run_xc_mol_grad(molecule, basis, xcfun, ref_grad)
