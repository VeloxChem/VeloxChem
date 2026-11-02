from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.errorhandler import VeloxChemError
from veloxchem.tdaeigensolver import TdaEigenSolver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tddftgradientdriver import TddftGradientDriver


@pytest.mark.solvers
class TestGrad:

    def run_tddft_grad(self, xcfun_label, tamm_dancoff, ri, ref_grad):

        molecule_string = """
        O 0.0000000000    0.0000000000   -0.0254395383
        H 0.0000000000    0.7695699584    0.5948147012
        H 0.0000000000   -0.7695699584    0.5948147012
        """

        basis_set_label = "def2-svp"

        molecule = Molecule.read_str(molecule_string, units='angstrom')
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri

        scf_dict = {'conv_thresh': 1e-8}
        method_dict = {'xcfun': xcfun_label, 'grid_level': 5}

        scf_drv.update_settings(scf_dict, method_dict)

        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        if tamm_dancoff:
            rsp_solver = TdaEigenSolver()
        else:
            rsp_solver = LinearResponseEigenSolver()

        if tamm_dancoff:
            tda = 'yes'
        else:
            tda = 'no'
        rsp_dict = {"tamm_dancoff": tda}

        rsp_solver.update_settings(rsp_dict, method_dict)
        rsp_solver.conv_thresh = 1e-5
        rsp_solver.nstates = 3

        rsp_solver.ostream.mute()
        rsp_results = rsp_solver.compute(molecule, basis, scf_drv.scf_results)

        tddft_grad = TddftGradientDriver(scf_drv)

        orbrsp_dict = {"conv_thresh": 1e-7}
        tddft_grad.update_settings({}, rsp_dict, orbrsp_dict, method_dict)

        tddft_grad.ostream.mute()
        tddft_grad.compute(molecule, basis, scf_drv, rsp_solver, rsp_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            grad = tddft_grad.get_gradient()[0]
            assert np.max(np.abs(grad - ref_grad)) < 1.0e-6

    @pytest.mark.parametrize('ri', [False, True])
    def test_tda_slater(self, ri):
        xcfun_label = 'slater'
        tamm_dancoff = True
        ref_grad = {
            False: np.array([[-0., -0.0, 0.096828839746124],
                             [0., -0.059512179812975, -0.048405146528462],
                             [-0., 0.059512179812977, -0.048405146528465]]),
            True: np.array([[0.0, 0.0, 0.096477705739322],
                            [0.0, -0.059333908732665, -0.048229578679697],
                            [0.0, 0.05933390873266, -0.048229578679703]]),
        }[ri]

        self.run_tddft_grad(xcfun_label, tamm_dancoff, ri, ref_grad)

    @pytest.mark.parametrize('ri', [False, True])
    def test_tda_blyp(self, ri):
        xcfun_label = 'blyp'
        tamm_dancoff = True
        ref_grad = {
            False: np.array([[0., 0.0, 0.081702810507772],
                             [-0., -0.052036797551521, -0.040837361620094],
                             [0., 0.052036797551518, -0.040837361620095]]),
            True: np.array([[0., 0.0, 0.081353914537983],
                            [0., -0.051842808410082, -0.040662912704928],
                            [0., 0.051842808410084, -0.040662912704952]]),
        }[ri]

        self.run_tddft_grad(xcfun_label, tamm_dancoff, ri, ref_grad)

    @pytest.mark.parametrize('ri', [False, True])
    def test_tda_lrcwpbeh(self, ri):
        xcfun_label = "lrc-wpbeh"
        tamm_dancoff = True
        ref_grad = {
            False: np.array([[0., -0.0, 0.076205604655348],
                             [-0., -0.050007199772849, -0.038080925818616],
                             [0., 0.050007199772858, -0.038080925818619]]),
            True: np.array([[0.0, 0.0, 0.076174713794153],
                            [0.0, -0.049994387376423, -0.038065459396769],
                            [0.0, 0.04999438738611, -0.038065459391165]]),
        }[ri]

        self.run_tddft_grad(xcfun_label, tamm_dancoff, ri, ref_grad)

    @pytest.mark.parametrize('ri', [False, True])
    def test_rpa_slater(self, ri):
        xcfun_label = 'slater'
        tamm_dancoff = False
        ref_grad = {
            False: np.array([[0., -0.0, 0.097148642876961],
                             [0., -0.060393461796934, -0.048565041680259],
                             [0., 0.060393461796936, -0.048565041680257]]),
            True: np.array([[0.0, 0.0, 0.096785558076883],
                            [0.0, -0.060211106217539, -0.048383498400645],
                            [0.0, 0.060211106217503, -0.048383498400653]]),
        }[ri]

        self.run_tddft_grad(xcfun_label, tamm_dancoff, ri, ref_grad)

    @pytest.mark.parametrize('ri', [False, True])
    def test_rpa_blyp(self, ri):
        xcfun_label = 'blyp'
        tamm_dancoff = False
        ref_grad = {
            False: np.array([[-0., -0.0, 0.0818110121047],
                             [-0., -0.052800208718226, -0.040891498953149],
                             [0., 0.052800208718228, -0.040891498953151]]),
            True: np.array([[0., 0.0, 0.081452905413363],
                            [-0., -0.052602723644523, -0.040712444660867],
                            [0., 0.052602723644504, -0.040712444660883]]),
        }[ri]

        self.run_tddft_grad(xcfun_label, tamm_dancoff, ri, ref_grad)

    @pytest.mark.parametrize('ri', [False, True])
    def test_rpa_lrcwpbeh(self, ri):
        xcfun_label = "lrc-wpbeh"
        tamm_dancoff = False
        ref_grad = {
            False: np.array([[-0., 0.0, 0.076746538972703],
                             [-0., -0.050926793834548, -0.038346546437516],
                             [0., 0.050926793834542, -0.03834654643751]]),
            True: np.array([[0.0, 0.0, 0.076715657764122],
                            [0.0, -0.050914200877127, -0.03833108824137],
                            [0.0, 0.050914200886843, -0.0383310882358]]),
        }[ri]
        self.run_tddft_grad(xcfun_label, tamm_dancoff, ri, ref_grad)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    @pytest.mark.parametrize('ri', [False, True])
    def test_analytical_gradient_rejects_pe(self, ri):

        molecule_string = """
        O 0.0000000000    0.0000000000   -0.0254395383
        H 0.0000000000    0.7695699584    0.5948147012
        H 0.0000000000   -0.7695699584    0.5948147012
        """
        molecule = Molecule.read_str(molecule_string, units='angstrom')
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        # normally set by pe_sanity_check when embedding is configured
        scf_drv._pe = True

        tddft_grad = TddftGradientDriver(scf_drv)

        with pytest.raises(
                VeloxChemError,
                match='Analytical TDDFT gradient is not supported with '
                      'polarizable embedding'):
            tddft_grad.compute_analytical(molecule, basis, {})

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    @pytest.mark.parametrize('ri', [False, True])
    def test_excited_state_gradient_rejects_unrestricted_reference(self, ri):

        scf_drv = ScfUnrestrictedDriver()
        scf_drv.ri = ri

        with pytest.raises(
                VeloxChemError,
                match='excited-state gradient not implemented for open-shell'):
            TddftGradientDriver(scf_drv)

    @pytest.mark.parametrize('ri', [False, True])
    def test_tda_point_charges(self, ri):
        """
        Analytical TDDFT gradient of the two lowest excited states with
        external point charges. The four-center references were verified
        against finite differences; the RI references are the analytical RI
        gradients, which differ from them by 1.2e-4.
        """

        molecule_string = """
        O 0.0000000000    0.0000000000   -0.0254395383
        H 0.0000000000    0.7695699584    0.5948147012
        H 0.0000000000   -0.7695699584    0.5948147012
        """
        molecule = Molecule.read_str(molecule_string, units='angstrom')
        basis = MolecularBasis.read(molecule, 'sto-3g', ostream=None)

        # one point charge fixed at (2.0, 1.5, 3.0) bohr with charge +1.0
        point_charges = np.zeros((6, 1))
        point_charges[0, 0] = 2.0
        point_charges[1, 0] = 1.5
        point_charges[2, 0] = 3.0
        point_charges[3, 0] = 1.0

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.point_charges = point_charges
        scf_drv.conv_thresh = 1.0e-10
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        rsp_drv = TdaEigenSolver()
        rsp_drv.nstates = 2
        rsp_drv.conv_thresh = 1.0e-8
        rsp_drv.ostream.mute()
        rsp_results = rsp_drv.compute(molecule, basis, scf_drv.scf_results)

        tddft_grad = TddftGradientDriver(scf_drv)
        tddft_grad.update_settings({}, {'tamm_dancoff': 'yes'}, {},
                                   {'xcfun': 'hf'})
        tddft_grad.ostream.mute()
        tddft_grad.compute(molecule, basis, scf_drv, rsp_drv, rsp_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            grads = tddft_grad.get_gradient()
            ref_grads = {
                False: [
                    np.array([[0.01226097, 0.26217545, 0.29824809],
                              [-0.03585263, -0.34321184, -0.27535775],
                              [0.00466140, 0.09629800, -0.03346097]]),
                    np.array([[0.00496588, -0.26089465, 0.43299797],
                              [0.00921538, -0.06707517, -0.09743285],
                              [-0.00279742, 0.32777274, -0.32197701]]),
                ],
                True: [
                    np.array([[0.01226556, 0.26227162, 0.29834682],
                              [-0.03586291, -0.34332940, -0.27545052],
                              [0.00466077, 0.09632161, -0.03347017]]),
                    np.array([[0.00497025, -0.26096852, 0.43311757],
                              [0.00921090, -0.06708937, -0.09746243],
                              [-0.00280000, 0.32786016, -0.32206699]]),
                ],
            }[ri]
            for s in range(2):
                assert np.max(np.abs(grads[s] - ref_grads[s])) < 1.0e-5
