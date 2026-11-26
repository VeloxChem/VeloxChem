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

    def run_ecd(self,
                xcfun_label,
                flag,
                ri,
                ref_edip,
                ref_vdip,
                ref_mdip,
                ref_exc_ene,
                ref_osc_str,
                ref_rot_str,
                conv_thresh=None):

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
        if conv_thresh is not None:
            lr_drv.conv_thresh = conv_thresh
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
            False: np.array([[0.006679, 0.036818, -0.003596],
                             [0.008224, 0.053103, 0.086419],
                             [0.046175, -0.299262, -0.636826],
                             [-0.069375, 0.069214, 0.522326],
                             [-0.201737, -1.007467, 0.327815]]),
            True: np.array([[0.006680, 0.036814, -0.003593],
                            [0.008226, 0.053056, 0.086403],
                            [0.046175, -0.299273, -0.636830],
                            [-0.069403, 0.069326, 0.522484],
                            [-0.201775, -1.007531, 0.327830]])
        }[ri]

        ref_vdip = {
            False: np.array([[0.024825, -0.048391, -0.207142],
                             [-0.023093, -0.050046, 0.148203],
                             [0.049042, -0.351884, -0.722428],
                             [-0.080618, 0.065405, 0.538779],
                             [-0.228610, -1.089260, 0.354960]]),
            True: np.array([[0.024818, -0.048396, -0.207125],
                            [-0.023082, -0.050118, 0.148182],
                            [0.049044, -0.351891, -0.722403],
                            [-0.080654, 0.065529, 0.538921],
                            [-0.228666, -1.089293, 0.354956]])
        }[ri]

        ref_mdip = {
            False: np.array([[-0.167255, -0.500028, 0.193785],
                             [0.048099, 0.357990, 0.598155],
                             [-0.502562, -0.552158, 0.290740],
                             [0.443193, 0.209180, 0.093308],
                             [0.214919, 0.319153, 1.093161]]),
            True: np.array([[-0.167257, -0.500008, 0.193784],
                            [0.048077, 0.357982, 0.598223],
                            [-0.502544, -0.552159, 0.290744],
                            [0.443326, 0.209263, 0.093249],
                            [0.214912, 0.319144, 1.093236]])
        }[ri]

        ref_exc_ene = {
            False: np.array(
                [0.32938947, 0.39989947, 0.40921479, 0.43429958, 0.45311866]),
            True: np.array(
                [0.32939424, 0.39989954, 0.40921942, 0.43430115, 0.45311913])
        }[ri]
        ref_osc_str = {
            False: np.array([0.000310, 0.002761, 0.135651, 0.081772, 0.351364]),
            True: np.array([0.000310, 0.002759, 0.135656, 0.081826, 0.351411])
        }[ri]
        ref_rot_str = {
            False: np.array(
                [-9.474075, 32.822745, -19.041582, 13.306341, -4.122787]),
            True: np.array(
                [-9.471445, 32.810070, -19.037454, 13.299826, -4.117262])
        }[ri]

        # conv_thresh below the default: the results depend on the path of
        # the Davidson iterations by less than at the default threshold, and
        # by much less than the tolerances of this test
        self.run_ecd('hf',
                     'rpa',
                     ri,
                     ref_edip,
                     ref_vdip,
                     ref_mdip,
                     ref_exc_ene,
                     ref_osc_str,
                     ref_rot_str,
                     conv_thresh=1.0e-7)

    @pytest.mark.parametrize('ri', [False, True])
    def test_hf_tda(self, ri):

        # vlxtag: RHF, Absorption, ECD, TDHF

        ref_edip = {
            False: np.array([[-0.005942, -0.045677, -0.011239],
                             [-0.011929, -0.063188, -0.084152],
                             [-0.050278, 0.312144, 0.676749],
                             [0.068409, -0.041465, -0.504704],
                             [0.190293, 1.040172, -0.334475]]),
            True: np.array([[-0.005944, -0.045670, -0.011234],
                            [-0.011928, -0.063139, -0.084140],
                            [-0.050280, 0.312162, 0.676765],
                            [0.068441, -0.041572, -0.504858],
                            [0.190330, 1.040240, -0.334490]])
        }[ri]

        ref_vdip = {
            False: np.array([[-0.028588, 0.094991, 0.298493],
                             [0.022311, 0.030726, -0.135916],
                             [-0.052879, 0.367115, 0.778773],
                             [0.065825, -0.013182, -0.408480],
                             [0.194162, 1.017049, -0.341504]]),
            True: np.array([[-0.028580, 0.094973, 0.298431],
                            [0.022297, 0.030781, -0.135893],
                            [-0.052878, 0.367124, 0.778739],
                            [0.065862, -0.013317, -0.408677],
                            [0.194201, 1.017097, -0.341496]])
        }[ri]

        ref_mdip = {
            False: np.array([[0.212632, 0.601168, -0.240445],
                             [-0.040316, -0.347648, -0.573949],
                             [0.543003, 0.618809, -0.307562],
                             [-0.350370, -0.083499, -0.134878],
                             [-0.202908, -0.314140, -1.003565]]),
            True: np.array([[0.212610, 0.601110, -0.240426],
                            [-0.040291, -0.347645, -0.574014],
                            [0.542983, 0.618796, -0.307572],
                            [-0.350545, -0.083627, -0.134810],
                            [-0.202890, -0.314132, -1.003642]])
        }[ri]

        ref_exc_ene = {
            False: np.array(
                [0.33148259, 0.40168416, 0.41121722, 0.43694078, 0.45467313]),
            True: np.array(
                [0.33148729, 0.40168417, 0.41122175, 0.43694226, 0.45467346])
        }[ri]
        ref_osc_str = {
            False: np.array([0.000497, 0.003004, 0.152959, 0.076064, 0.372845]),
            True: np.array([0.000497, 0.003002, 0.152970, 0.076114, 0.372895])
        }[ri]
        ref_rot_str = {
            False: np.array(
                [-9.779799, 31.316844, -19.357930, 15.620101, -7.623613]),
            True: np.array(
                [-9.776558, 31.306340, -19.355404, 15.614216, -7.620553])
        }[ri]

        # conv_thresh below the default: the results depend on the path of
        # the Davidson iterations by less than at the default threshold, and
        # by much less than the tolerances of this test
        self.run_ecd('hf',
                     'tda',
                     ri,
                     ref_edip,
                     ref_vdip,
                     ref_mdip,
                     ref_exc_ene,
                     ref_osc_str,
                     ref_rot_str,
                     conv_thresh=1.0e-7)
