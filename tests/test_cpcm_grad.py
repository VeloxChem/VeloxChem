import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver


@pytest.mark.solvers
class TestGradCpcm:

    def run_grad(self, molecule, xcfun_label, basis_label, ref_grad,
                 cpcm_custom_vdw_radii, tol, *, ri):

        basis = MolecularBasis.read(molecule, basis_label)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.xcfun = xcfun_label

        scf_drv.solvation_model = 'cpcm'
        scf_drv.cpcm_grid_per_sphere = (110, 110)
        scf_drv.cpcm_custom_vdw_radii = cpcm_custom_vdw_radii

        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.compute(molecule, basis, scf_results)
        grad = grad_drv.get_gradient()

        assert np.max(np.abs(grad - ref_grad)) < tol

    @pytest.mark.parametrize('ri', [False, True])
    def test_nh3_def2svp(self, ri):

        molstr = """
        N         -1.96309        1.59755       -0.01963
        H         -1.95876        2.61528        0.03109
        H         -2.48929        1.27814        0.79244
        H         -2.52930        1.35928       -0.83265
        """
        mol = Molecule.read_molecule_string(molstr, units='bohr')

        ref_grad = {
            False: np.array([[-2.57769517, 1.08987991, 0.11790935],
                             [-0.04727107, -2.49820434, -0.12359626],
                             [1.263075, 0.80436803, -2.00404575],
                             [1.36189124, 0.6039564, 2.00973266]]),
            True: np.array([[-2.57766268, 1.08986617, 0.11790786],
                            [-0.04729726, -2.49823592, -0.12359720],
                            [1.26307101, 0.80439239, -2.00407853],
                            [1.36188892, 0.60397737, 2.00976787]])
        }[ri]

        self.run_grad(mol, 'hf', 'def2-svp', ref_grad, None, 1.0e-5, ri=ri)

        ref_grad = {
            False: np.array([[-2.59710917, 1.09805904, 0.1187649],
                             [-0.04439474, -2.5096928, -0.12424109],
                             [1.27129737, 0.80635098, -2.01236183],
                             [1.37051574, 0.60512271, 2.01779083]]),
            True: np.array([[-2.59837472, 1.09859410, 0.11882278],
                            [-0.04398870, -2.50990841, -0.12426185],
                            [1.27172627, 0.80619301, -2.01241544],
                            [1.37094637, 0.60496123, 2.01780731]])
        }[ri]

        self.run_grad(mol, 'slater', 'def2-svp', ref_grad, None, 1.0e-3, ri=ri)

        ref_grad = np.array([[-2.56212, 1.0832703, 0.11717353],
                             [-0.04534075, -2.47949799, -0.12270971],
                             [1.25484747, 0.7974705, -1.98858115],
                             [1.35289684, 0.59861351, 1.99408099]])

        self.run_grad(mol, 'b3lyp', 'def2-svp', ref_grad, None, 1.0e-3, ri=ri)

        ref_grad = np.array([[-2.55730976, 1.0813965, 0.11713555],
                             [-0.04737833, -2.47972691, -0.12267538],
                             [1.25332429, 0.79870663, -1.98929246],
                             [1.3513638, 0.59962406, 1.99483228]])

        self.run_grad(mol,
                      'b3lyp',
                      'def2-svp',
                      ref_grad, ['N', 1.8],
                      1.0e-3,
                      ri=ri)

    def run_grad_with_ecp(self, molecule, basis_label, ref_grad, tol, *, ri):

        basis = MolecularBasis.read(molecule, basis_label)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.solvation_model = 'cpcm'
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.compute(molecule, basis, scf_results)
        grad = grad_drv.get_gradient()

        assert np.max(np.abs(grad - ref_grad)) < tol

    @pytest.mark.parametrize('ri', [False, True])
    def test_auh(self, ri):

        xyzstr = """2
        xyz
        Au 0 0 0
        H  0 0 1.55
        """
        mol = Molecule.read_xyz_string(xyzstr)

        ref_grad = np.array([[0, 0, 0.0128007496], [0, 0, -0.0128007496]])

        self.run_grad_with_ecp(mol, 'def2-svp', ref_grad, 1.0e-5, ri=ri)
