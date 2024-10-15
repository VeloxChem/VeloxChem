import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver


@pytest.mark.solvers
class TestGradUnrestrictedSCF:

    def run_grad(self,
                 molecule,
                 xcfun_label,
                 basis_label,
                 ref_grad,
                 tol,
                 d4_flag='none'):

        basis = MolecularBasis.read(molecule, basis_label)

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_drv.dispersion = (d4_flag.lower() == 'd4')
        scf_results = scf_drv.compute(molecule, basis)

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.compute(molecule, basis, scf_results)
        grad = grad_drv.get_gradient()

        assert np.max(np.abs(grad - ref_grad)) < tol

    def test_nh3_tzvp(self):

        molstr = """
        N         -1.96309        1.59755       -0.01963
        H         -1.95876        2.61528        0.03109
        H         -2.48929        1.27814        0.79244
        H         -2.52930        1.35928       -0.83265
        """
        mol = Molecule.read_molecule_string(molstr, units='bohr')
        mol.set_charge(1)
        mol.set_multiplicity(2)

        ref_grad = np.array([
            [-2.406223382, 1.017390065, 0.110065654],
            [-0.111995746, -2.491875636, -0.121685384],
            [1.209288923, 0.838285105, -2.017829831],
            [1.308930132, 0.636200679, 2.029449566],
        ])

        self.run_grad(mol, 'hf', 'def2-tzvp', ref_grad, 1.0e-4)

    def test_c2h4_svp(self):

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
        mol.set_charge(0)
        mol.set_multiplicity(3)

        ref_grad = np.array([
            [0.206766783, -0., 0.],
            [-0.206766783, 0., 0.],
            [0.001566781, 0.004679599, -0.],
            [0.001566781, -0.004679599, 0.],
            [-0.001566781, 0.004679599, -0.],
            [-0.001566781, -0.004679599, -0.],
        ])

        self.run_grad(mol, 'hf', 'def2-svp', ref_grad, 1.0e-5)
