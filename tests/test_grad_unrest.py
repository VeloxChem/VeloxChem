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

        ref_grad = np.array([[-2.407691473, 1.018158314, 0.11031899],
                             [-0.107273535, -2.482263778, -0.121327734],
                             [1.207918349, 0.83272141, -2.008731851],
                             [1.307046644, 0.631384109, 2.019740596]])

        self.run_grad(mol, 'blyp', 'def2-tzvp', ref_grad, 1.0e-3)

        ref_grad = np.array([[-2.404693831, 1.016859978, 0.110144742],
                             [-0.108078426, -2.481352867, -0.121261112],
                             [1.206827125, 0.832875139, -2.008254371],
                             [1.305945115, 0.631617805, 2.019370742]])

        self.run_grad(mol, 'b3lyp', 'def2-tzvp', ref_grad, 1.0e-3)

        ref_grad = np.array([
            [-2.405765356322, 1.017221798072, 0.110096775899],
            [-0.109028651398, -2.484488849564, -0.121379513686],
            [1.207761315992, 0.834344125110, -2.011111344412],
            [1.307032618595, 0.632923133068, 2.022394084388],
        ])

        self.run_grad(mol, 'tpssh', 'def2-tzvp', ref_grad, 1.0e-3)

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

        ref_grad = np.array([[0.189529949, 0., 0.], [-0.189529949, -0., 0.],
                             [0.008659375, -0.005825038, -0.],
                             [0.008659375, 0.005825038, -0.],
                             [-0.008659375, -0.005825038, 0.],
                             [-0.008659375, 0.005825038, -0.]])

        self.run_grad(mol, 'blyp', 'def2-svp', ref_grad, 1.0e-4)

        ref_grad = np.array([[0.190927027, 0., -0.], [-0.190927026, 0., -0.],
                             [0.00593267, -0.00172427, -0.],
                             [0.00593267, 0.001724269, 0.],
                             [-0.00593267, -0.00172427, -0.],
                             [-0.00593267, 0.00172427, -0.]])

        self.run_grad(mol, 'b3lyp', 'def2-svp', ref_grad, 1.0e-4)

        ref_grad = np.array([
            [0.182943548, 0.000000000, -0.000000000],
            [-0.182943548, 0.000000000, 0.000000000],
            [0.005993681, -0.000833749, -0.000000000],
            [0.005993681, 0.000833749, -0.000000000],
            [-0.005993681, -0.000833749, 0.000000000],
            [-0.005993681, 0.000833749, 0.000000000],
        ])

        self.run_grad(mol, 'm06-l', 'def2-svp', ref_grad, 1.0e-4)

        ref_grad = np.array([
            [0.190663659, -0.000000000, 0.000000000],
            [-0.190663659, 0.000000000, 0.000000000],
            [0.005775999, -0.001793691, 0.000000000],
            [0.005775999, 0.001793691, 0.000000000],
            [-0.005775999, -0.001793691, -0.000000000],
            [-0.005775999, 0.001793691, -0.000000000],
        ])

        self.run_grad(mol, 'lrc-wpbeh', 'def2-svp', ref_grad, 1.0e-4)
