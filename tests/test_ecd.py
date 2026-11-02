import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tdaeigensolver import TdaEigenSolver


@pytest.mark.solvers
class TestECD:

    def run_ecd(self, xcfun_label, flag, ri, ref_edip, ref_vdip, ref_mdip,
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

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.xcfun = xcfun_label
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

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

    @pytest.mark.parametrize('ri', [False, True])
    def test_hf_rpa(self, ri):

        # vlxtag: RHF, Absorption, ECD, CIS

        ref_edip = {
            False: np.array([[0.006676, 0.036815, -0.003599],
                             [0.008225, 0.053100, 0.086421],
                             [0.046176, -0.299265, -0.636825],
                             [-0.069385, 0.069203, 0.522331],
                             [-0.201731, -1.007463, 0.327814]]),
            True: np.array([[0.006683, 0.036817, -0.003593],
                            [0.008224, 0.053047, 0.086406],
                            [0.046180, -0.299280, -0.636832],
                            [-0.069409, 0.069334, 0.522480],
                            [-0.201775, -1.007533, 0.327830]])
        }[ri]

        ref_vdip = {
            False: np.array([[0.024821, -0.048399, -0.207146],
                             [-0.023091, -0.050053, 0.148208],
                             [0.049042, -0.351892, -0.722428],
                             [-0.080638, 0.065384, 0.538794],
                             [-0.228599, -1.089253, 0.354952]]),
            True: np.array([[0.024822, -0.048399, -0.207143],
                            [-0.023087, -0.050130, 0.148187],
                            [0.049050, -0.351896, -0.722404],
                            [-0.080662, 0.065547, 0.538927],
                            [-0.228669, -1.089292, 0.354957]])
        }[ri]

        ref_mdip = {
            False: np.array([[-0.167257, -0.500029, 0.193794],
                             [0.048102, 0.357994, 0.598159],
                             [-0.502564, -0.552159, 0.290751],
                             [0.443203, 0.209185, 0.093352],
                             [0.214914, 0.319150, 1.093141]]),
            True: np.array([[-0.167267, -0.500016, 0.193787],
                            [0.048080, 0.357991, 0.598245],
                            [-0.502547, -0.552152, 0.290748],
                            [0.443335, 0.209268, 0.093242],
                            [0.214914, 0.319146, 1.093239]])
        }[ri]

        ref_exc_ene = {
            False: np.array(
                [0.32938946, 0.39989948, 0.40921479, 0.43429957, 0.45311866]),
            True: np.array(
                [0.32939424, 0.39989954, 0.40921942, 0.43430115, 0.45311913])
        }[ri]
        ref_osc_str = np.array([0.0003, 0.0028, 0.1357, 0.0818, 0.3514])
        ref_rot_str = {
            False: np.array([-9.4755, 32.8233, -19.0445, 13.3045, -4.1229]),
            True: np.array(
                [-9.472829, 32.810620, -19.040397, 13.298018, -4.117345])
        }[ri]

        self.run_ecd('hf', 'rpa', ri, ref_edip, ref_vdip, ref_mdip, ref_exc_ene,
                     ref_osc_str, ref_rot_str)

    @pytest.mark.parametrize('ri', [False, True])
    def test_hf_tda(self, ri):

        # vlxtag: RHF, Absorption, ECD, TDHF

        ref_edip = {
            False: np.array([[-0.005942, -0.045678, -0.011238],
                             [-0.011930, -0.063188, -0.084152],
                             [-0.050277, 0.312144, 0.676748],
                             [0.068417, -0.041465, -0.504707],
                             [0.190288, 1.040164, -0.334470]]),
            True: np.array([[-0.005944, -0.045670, -0.011234],
                            [-0.011925, -0.063136, -0.084142],
                            [-0.050279, 0.312161, 0.676766],
                            [0.068444, -0.041575, -0.504856],
                            [0.190326, 1.040241, -0.334490]])
        }[ri]

        ref_vdip = {
            False: np.array([[-0.028588, 0.094991, 0.298493],
                             [0.022311, 0.030726, -0.135916],
                             [-0.052879, 0.367115, 0.778772],
                             [0.065827, -0.013180, -0.408485],
                             [0.194158, 1.017040, -0.341502]]),
            True: np.array([[-0.028580, 0.094973, 0.298430],
                            [0.022301, 0.030786, -0.135894],
                            [-0.052876, 0.367123, 0.778739],
                            [0.065868, -0.013321, -0.408677],
                            [0.194199, 1.017090, -0.341494]])
        }[ri]

        ref_mdip = {
            False: np.array([[0.212632, 0.601168, -0.240445],
                             [-0.040315, -0.347648, -0.573949],
                             [0.543002, 0.618808, -0.307563],
                             [-0.350373, -0.083502, -0.134883],
                             [-0.202904, -0.314138, -1.003555]]),
            True: np.array([[0.212610, 0.601110, -0.240426],
                            [-0.040291, -0.347647, -0.574024],
                            [0.542983, 0.618796, -0.307574],
                            [-0.350548, -0.083628, -0.134813],
                            [-0.202890, -0.314132, -1.003635]])
        }[ri]

        ref_exc_ene = {
            False: np.array(
                [0.33148263, 0.40168416, 0.41121722, 0.43694079, 0.45467313]),
            True: np.array(
                [0.33148729, 0.40168417, 0.41122175, 0.43694226, 0.45467346])
        }[ri]
        ref_osc_str = {
            False: np.array([0.0005, 0.0030, 0.1530, 0.0761, 0.3728]),
            True: np.array([0.000497, 0.003001, 0.152970, 0.076113, 0.372895])
        }[ri]
        ref_rot_str = {
            False: np.array([-9.7797, 31.3169, -19.3585, 15.6197, -7.6241]),
            True: np.array(
                [-9.776502, 31.306423, -19.355995, 15.613821, -7.621020])
        }[ri]

        self.run_ecd('hf', 'tda', ri, ref_edip, ref_vdip, ref_mdip, ref_exc_ene,
                     ref_osc_str, ref_rot_str)
