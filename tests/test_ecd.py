import numpy as np
import pytest

from veloxchem import mpi_master
from veloxchem import Molecule, MolecularBasis
from veloxchem import ScfRestrictedDriver
from veloxchem import LinearResponseEigenSolver
#from veloxchem import LinearResponseEigenSolver, TdaEigenSolver


@pytest.mark.solvers
class TestECD:

    def run_ecd(self, xcfun_label, flag, ref_edip, ref_vdip, ref_mdip,
                ref_exc_ene, ref_osc_str, ref_rot_str):

        xyz_string = """6
        xyz
        O -3.42904  1.55532  0.01546
        C -1.99249  1.74379  0.02665
        H -1.74709  2.74160  0.44749
        H -1.59636  1.67836 -1.00861
        H -1.51398  0.95881  0.64937
        H -3.84726  2.33620 -0.34927
        """
        mol = Molecule.read_xyz_string(xyz_string)
        mol.check_multiplicity()

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        if flag == 'rpa':
            lr_drv = LinearResponseEigenSolver()
        elif flag == 'tda':
            lr_drv = TdaEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.nstates = 5
        lr_results = lr_drv.compute(mol, bas, scf_results)

        if lr_drv.rank == mpi_master():
            assert np.max(np.abs(ref_exc_ene -
                                 lr_results['eigenvalues'])) < 1.0e-6
            assert np.max(
                np.abs(ref_osc_str -
                       lr_results['oscillator_strengths'])) < 1.0e-4
            assert np.max(np.abs(ref_rot_str -
                                 lr_results['rotatory_strengths'])) < 1.0e-4

            edip = lr_results['electric_transition_dipoles']
            vdip = lr_results['velocity_transition_dipoles']
            mdip = lr_results['magnetic_transition_dipoles']

            for s in range(ref_exc_ene.size):
                factor = -1.0 if np.sum(edip[s] * ref_edip[s]) < 0.0 else 1.0
                assert np.max(np.abs(edip[s] * factor - ref_edip[s])) < 1.0e-4
                assert np.max(np.abs(vdip[s] * factor - ref_vdip[s])) < 1.0e-4
                assert np.max(np.abs(mdip[s] * factor - ref_mdip[s])) < 1.0e-4

    def test_hf_rpa(self):

        ref_edip = np.array([[0.006676, 0.036815, -0.003599],
                             [0.008225, 0.053100, 0.086421],
                             [0.046176, -0.299265, -0.636825],
                             [-0.069385, 0.069203, 0.522331],
                             [-0.201731, -1.007463, 0.327814]])

        ref_vdip = np.array([[0.024821, -0.048399, -0.207146],
                             [-0.023091, -0.050053, 0.148208],
                             [0.049042, -0.351892, -0.722428],
                             [-0.080638, 0.065384, 0.538794],
                             [-0.228599, -1.089253, 0.354952]])

        ref_mdip = np.array([[-0.167257, -0.500029, 0.193794],
                             [0.048102, 0.357994, 0.598159],
                             [-0.502564, -0.552159, 0.290751],
                             [0.443203, 0.209185, 0.093352],
                             [0.214914, 0.319150, 1.093141]])

        ref_exc_ene = np.array(
            [0.32938946, 0.39989948, 0.40921479, 0.43429957, 0.45311866])
        ref_osc_str = np.array([0.0003, 0.0028, 0.1357, 0.0818, 0.3514])
        ref_rot_str = np.array([-9.4733, 32.8231, -19.0429, 13.3116, -4.1259])

        self.run_ecd('hf', 'rpa', ref_edip, ref_vdip, ref_mdip, ref_exc_ene,
                     ref_osc_str, ref_rot_str)

    """
    def test_hf_tda(self):

        ref_edip = np.array([[-0.005942, -0.045678, -0.011238],
                             [-0.011930, -0.063188, -0.084152],
                             [-0.050277, 0.312144, 0.676748],
                             [0.068417, -0.041465, -0.504707],
                             [0.190288, 1.040164, -0.334470]])

        ref_vdip = np.array([[-0.028588, 0.094991, 0.298493],
                             [0.022311, 0.030726, -0.135916],
                             [-0.052879, 0.367115, 0.778772],
                             [0.065827, -0.013180, -0.408485],
                             [0.194158, 1.017040, -0.341502]])

        ref_mdip = np.array([[0.212632, 0.601168, -0.240445],
                             [-0.040315, -0.347648, -0.573949],
                             [0.543002, 0.618808, -0.307563],
                             [-0.350373, -0.083502, -0.134883],
                             [-0.202904, -0.314138, -1.003555]])

        ref_exc_ene = np.array(
            [0.33148263, 0.40168416, 0.41121722, 0.43694079, 0.45467313])
        ref_osc_str = np.array([0.0005, 0.0030, 0.1530, 0.0761, 0.3728])
        ref_rot_str = np.array([-9.7798, 31.3168, -19.3580, 15.6209, -7.6236])

        self.run_ecd('hf', 'tda', ref_edip, ref_vdip, ref_mdip, ref_exc_ene,
                     ref_osc_str, ref_rot_str)
    """
