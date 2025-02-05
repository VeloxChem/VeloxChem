import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver


@pytest.mark.solvers
class TestGradCpcm:

    def run_grad(self, molecule, xcfun_label, basis_label, ref_grad, tol):

        basis = MolecularBasis.read(molecule, basis_label)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.solvation_model = 'cpcm'
        scf_drv.cpcm_grid_per_sphere = 110
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.compute(molecule, basis, scf_results)
        grad = grad_drv.get_gradient()

        assert np.max(np.abs(grad - ref_grad)) < tol

    def test_nh3_def2svp(self):

        molstr = """
        N         -1.96309        1.59755       -0.01963
        H         -1.95876        2.61528        0.03109
        H         -2.48929        1.27814        0.79244
        H         -2.52930        1.35928       -0.83265
        """
        mol = Molecule.read_molecule_string(molstr, units='bohr')

        ref_grad = np.array([[-2.57769517, 1.08987991, 0.11790935],
                             [-0.04727107, -2.49820434, -0.12359626],
                             [1.263075, 0.80436803, -2.00404575],
                             [1.36189124, 0.6039564, 2.00973266]])

        self.run_grad(mol, 'hf', 'def2-svp', ref_grad, 1.0e-5)

        ref_grad = np.array([[-2.59710917, 1.09805904, 0.1187649],
                             [-0.04439474, -2.5096928, -0.12424109],
                             [1.27129737, 0.80635098, -2.01236183],
                             [1.37051574, 0.60512271, 2.01779083]])

        self.run_grad(mol, 'slater', 'def2-svp', ref_grad, 1.0e-3)

        ref_grad = np.array([[-2.56212, 1.0832703, 0.11717353],
                             [-0.04534075, -2.47949799, -0.12270971],
                             [1.25484747, 0.7974705, -1.98858115],
                             [1.35289684, 0.59861351, 1.99408099]])

        self.run_grad(mol, 'b3lyp', 'def2-svp', ref_grad, 1.0e-3)
