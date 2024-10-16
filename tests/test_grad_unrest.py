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

        """
        grad_drv.delta_h = 0.0005
        grad_drv.compute_numerical_gradient(molecule, basis, scf_results)
        numgrad = grad_drv.get_gradient()
        np.set_printoptions(suppress=True, precision=9)
        print()
        print('    ref_grad = np.array(',
              np.array2string(numgrad, separator=', '), ')')
        ref_grad = numgrad
        """

        print(np.max(np.abs(grad - ref_grad)))
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
            [-2.406225127, 1.017383924, 0.110065855],
            [-0.111995782, -2.491871473, -0.12168558],
            [1.209289899, 0.8382861, -2.017828853],
            [1.308930992, 0.636201502, 2.029448579],
        ])

        self.run_grad(mol, 'hf', 'def2-tzvp', ref_grad, 1.0e-4)

        ref_grad = np.array([
            [-2.429691597, 1.027438566, 0.111298759],
            [-0.10204324, -2.490299784, -0.121855711],
            [1.216181028, 0.832324072, -2.013620701],
            [1.315553795, 0.6305372, 2.024177654],
        ])

        self.run_grad(mol, 'slda', 'def2-tzvp', ref_grad, 1.0e-3)

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
            [0.206766526, -0., 0.],
            [-0.206766526, -0., 0.],
            [0.001566777, 0.00467968, 0.],
            [0.001566777, -0.00467968, 0.],
            [-0.001566777, 0.00467968, 0.],
            [-0.001566777, -0.00467968, -0.],
        ])

        self.run_grad(mol, 'hf', 'def2-svp', ref_grad, 1.0e-5)

        ref_grad = np.array([
            [0.173586245, 0., 0.],
            [-0.173586245, -0., -0.],
            [0.008790969, -0.006423456, -0.],
            [0.008790969, 0.006423456, -0.],
            [-0.008790969, -0.006423456, -0.],
            [-0.008790969, 0.006423456, 0.],
        ])

        self.run_grad(mol, 'slda', 'def2-svp', ref_grad, 1.0e-4)
