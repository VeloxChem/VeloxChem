import numpy as np
import pytest

from veloxchem import Molecule, MolecularBasis
from veloxchem import ScfRestrictedDriver
from veloxchem import ScfGradientDriver


@pytest.mark.solvers
class TestGrad:

    def run_grad(self, molecule, basis):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.delta_h = 0.0005

        grad_drv.compute_numerical_gradient(molecule, basis, scf_results)
        num_grad = grad_drv.get_gradient()

        grad_drv.compute(molecule, basis, scf_results)
        ana_grad = grad_drv.get_gradient()

        assert np.max(np.abs(num_grad - ana_grad)) < 1.0e-5

    def test_nh3_sto3g(self):

        molstr = """
        N         -1.96309        1.59755       -0.01963
        H         -1.95876        2.61528        0.03109
        H         -2.48929        1.27814        0.79244
        H         -2.52930        1.35928       -0.83265
        """
        mol = Molecule.read_molecule_string(molstr, units='bohr')

        bas = MolecularBasis.read(mol, 'sto-3g')

        self.run_grad(mol, bas)

    def test_c2h4_sto3g(self):

        xyzstr = """
        6

        C      -0.667     -0.000     -0.000
        C       0.667      0.000     -0.000
        H      -1.227      0.930      0.000
        H      -1.227     -0.930     -0.000
        H       1.227      0.930     -0.000
        H       1.227     -0.930      0.000
        """
        mol = Molecule.read_xyz_string(xyzstr)

        bas = MolecularBasis.read(mol, 'sto-3g')

        self.run_grad(mol, bas)
