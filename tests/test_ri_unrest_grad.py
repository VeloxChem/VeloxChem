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
        #scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_drv.ri_coulomb = True
        scf_drv.conv_thresh = 1.0e-5
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
            C          -1.853343000000       -0.639451000000        1.296233000000                         
            H          -2.408845000000       -1.565702000000        1.042764000000                         
            H          -2.241609000000       -0.224427000000        2.249005000000                         
            H          -1.988307000000        0.106137000000        0.485892000000
            """
        
        ref_scf_energy = -39.7656435146
        tol_energy = 1.0e-4

        ref_scf_grad = np.array([
           [ 0.037856611800, -0.008221343122,   0.003901538936],
           [-0.012024423766,  0.005781000981,  -0.000643701820],
           [-0.012530159571,  0.001721679671,  -0.004294268094],
           [-0.013297084764,  0.000721229411,   0.001041950109]
        ])
        tol_grad = 5.0e-5

        self.run_scf_grad(mol_xyz, 'blyp', ref_scf_energy, ref_scf_grad, tol_energy,
                          tol_grad)


