import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver


@pytest.mark.solvers
class TestScfUnrestrictedGradientDriverWithRI:

    def run_scf_grad(self, mol_xyz, xcfun_label, ref_scf_energy, ref_scf_grad,
                     tol_energy, tol_grad):

        basis_label = 'def2-svp'

        mol = Molecule.read_xyz_string(mol_xyz)
        mol.set_multiplicity(2)
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.ri_coulomb = True
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        scf_grad_drv = ScfGradientDriver(scf_drv)
        scf_grad_drv.compute(mol, bas, scf_results)
        grad = scf_grad_drv.get_gradient()

        if scf_drv.rank == mpi_master():
            assert abs(ref_scf_energy - scf_results['scf_energy']) < tol_energy
            assert np.max(np.abs(ref_scf_grad - grad)) < tol_grad

    def test_ri_blyp(self):

        mol_xyz = """4
            xyz
            C          -1.85334300       -0.63945100        1.29623300
            H          -2.40884500       -1.56570200        1.04276400
            H          -2.24160900       -0.22442700        2.24900500
            H          -1.98830700        0.10613700        0.48589200
        """

        ref_scf_energy = -39.7656435146
        tol_energy = 1.0e-4

        ref_scf_grad = np.array(
            [[0.037856611800, -0.008221343122, 0.003901538936],
             [-0.012024423766, 0.005781000981, -0.000643701820],
             [-0.012530159571, 0.001721679671, -0.004294268094],
             [-0.013297084764, 0.000721229411, 0.001041950109]])
        tol_grad = 5.0e-5

        self.run_scf_grad(mol_xyz, 'blyp', ref_scf_energy, ref_scf_grad,
                          tol_energy, tol_grad)
