import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import XCMolecularGradient
from veloxchem.griddriver import GridDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.dftutils import get_default_grid_level


@pytest.mark.solvers
class TestXCMolGrad:

    def run_xc_mol_grad(self, molecule, basis, xcfun, ri, ref_grad):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()

        scf_drv.xcfun = xcfun
        scf_drv.conv_thresh = 1.0e-8

        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        density = scf_drv.density
        density = density.broadcast(scf_drv.comm, root=mpi_master())

        grid_drv = GridDriver(scf_drv.comm)
        grid_level = (get_default_grid_level(scf_drv.xcfun)
                      if scf_drv.grid_level is None else scf_drv.grid_level)
        grid_drv.set_level(grid_level)
        mol_grid = grid_drv.generate(molecule)

        grad_drv = XCMolecularGradient()
        mol_grad = grad_drv.integrate_vxc_gradient(molecule, basis,
                                                   [density.alpha_to_numpy(0)],
                                                   mol_grid, xcfun)
        mol_grad = scf_drv.comm.reduce(mol_grad, root=mpi_master())

        if scf_drv.rank == mpi_master():
            assert np.max(np.abs(mol_grad - ref_grad)) < 1.0e-4

    @pytest.mark.parametrize('ri', [False, True])
    def test_xc_mol_grad_slater(self, ri):

        mol_str = """
            O  0.0000000000   0.0000000000  -0.0254395383
            H  0.0000000000   0.7695699584   0.5948147012
            H  0.0000000000  -0.7695699584   0.5948147012
        """

        molecule = Molecule.read_molecule_string(mol_str, units='angstrom')
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        xcfun = 'slater'

        ref_grad = {
            False: np.array([
                [0., 0., -0.4075911],
                [-0., 0.2540621, 0.2037956],
                [0., -0.2540621, 0.2037956],
            ]),
            True: np.array([
                [0.0000000, 0.0000000, -0.4076932],
                [0.0000000, 0.2540673, 0.2038515],
                [0.0000000, -0.2540673, 0.2038515],
            ]),
        }[ri]

        self.run_xc_mol_grad(molecule, basis, xcfun, ri, ref_grad)

    @pytest.mark.parametrize('ri', [False, True])
    def test_xc_mol_grad_blyp(self, ri):

        mol_str = """
            O  0.0000000000   0.0000000000  -0.0254395383
            H  0.0000000000   0.7695699584   0.5948147012
            H  0.0000000000  -0.7695699584   0.5948147012
        """

        molecule = Molecule.read_molecule_string(mol_str, units='angstrom')

        xcfun = 'blyp'

        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        ref_grad = {
            False: np.array([
                [-0., 0.0000000, -0.4453751],
                [0., 0.2801306, 0.2226875],
                [0., -0.2801306, 0.2226875],
            ]),
            True: np.array([
                [0.0000000, 0.0000000, -0.4454790],
                [0.0000000, 0.2801448, 0.2227503],
                [0.0000000, -0.2801448, 0.2227503],
            ]),
        }[ri]
        self.run_xc_mol_grad(molecule, basis, xcfun, ri, ref_grad)

        basis = MolecularBasis.read(molecule, 'def2-tzvp', ostream=None)
        ref_grad = {
            False: np.array([
                [0., 0.0000000, -0.3753406],
                [0., 0.2442551, 0.1876845],
                [0., -0.2442551, 0.1876845],
            ]),
            True: np.array([
                [0.0000000, 0.0000000, -0.3748129],
                [0.0000000, 0.2441086, 0.1874206],
                [0.0000000, -0.2441086, 0.1874206],
            ]),
        }[ri]
        self.run_xc_mol_grad(molecule, basis, xcfun, ri, ref_grad)

    @pytest.mark.parametrize('ri', [False, True])
    def test_xc_mol_grad_tpss(self, ri):

        mol_str = """
            O  0.0000000000   0.0000000000  -0.0254395383
            H  0.0000000000   0.7695699584   0.5948147012
            H  0.0000000000  -0.7695699584   0.5948147012
        """

        molecule = Molecule.read_molecule_string(mol_str, units='angstrom')

        xcfun = 'tpss'

        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        ref_grad = {
            False: np.array([
                [-0., 0., -0.4584369],
                [0., 0.2891046, 0.2292038],
                [0., -0.2891046, 0.2292038],
            ]),
            True: np.array([
                [0.0000000, 0.0000000, -0.4585718],
                [0.0000000, 0.2891194, 0.2292717],
                [0.0000000, -0.2891194, 0.2292717],
            ]),
        }[ri]
        self.run_xc_mol_grad(molecule, basis, xcfun, ri, ref_grad)
