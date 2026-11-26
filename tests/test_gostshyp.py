import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem import GostshypCavityDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.outputstream import OutputStream
from veloxchem.errorhandler import VeloxChemError


def test_write_grid_to_file(tmp_path):
    """Checks that the grid writer uses its mandatory XYZ filename."""

    tessellation_drv = GostshypCavityDriver()
    vdw_surface = np.zeros((15, 2))
    vdw_surface[:3, 0] = [1.0, 2.0, 3.0]
    vdw_surface[:3, 1] = [4.0, 5.0, 6.0]
    xyz_filename = tmp_path / 'tessellation.xyz'

    tessellation_drv.write_grid_to_file(vdw_surface, xyz_filename)

    lines = xyz_filename.read_text().splitlines()
    assert lines[0].strip() == '2'
    assert lines[2].split() == ['x', '1.0', '2.0', '3.0']
    assert lines[3].split() == ['x', '4.0', '5.0', '6.0']


@pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                    reason='skip pytest.raises for multiple MPI processes')
def test_parse_pressure_units_rejects_non_string_units():
    """Checks that parse_pressure_units rejects non-string units."""

    from veloxchem.gostshypdriver import GostshypDriver

    with pytest.raises(VeloxChemError,
                       match='GOSTSHYP: Invalid unit for pressure'):
        GostshypDriver.parse_pressure_units(1.0, None)


@pytest.mark.solvers
class TestGostshyp:
    """GOSTSHYP hydrostatic pressure tests: SCF, analytical gradient and
    linear response for small molecules.
    """

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    @pytest.mark.parametrize('discretization', ['iwig', 'becke'])
    @pytest.mark.parametrize('ri', [False, True])
    def test_rejects_invalid_discretization(self, discretization, ri):

        xyz_string = """3

        O   -3.3278470    3.1951799   -0.0000000
        H   -4.2057717    2.7370843   -0.0000000
        H   -2.6643996    2.4600330   -0.0000000
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.pressure = 1.0
        scf_drv.gostshyp_discretization = discretization

        with pytest.raises(VeloxChemError,
                           match='GOSTSHYP: Invalid discretization'):
            scf_results_not_used = scf_drv.compute(mol, bas)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    @pytest.mark.parametrize('tco_tol', [0.0, -1.0e-14])
    @pytest.mark.parametrize('ri', [False, True])
    def test_rejects_nonpositive_tco_tolerance(self, tco_tol, ri):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()

        with pytest.raises(
                VeloxChemError,
                match='Three-center overlap integral screening threshold '
                      'must be positive'):
            scf_drv.update_settings({}, {
                'pressure': 1.0,
                'gostshyp_tco_tol': tco_tol,
            })

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    @pytest.mark.parametrize('tssf', [0.0, -1.0, np.nan, np.inf, -np.inf])
    @pytest.mark.parametrize('ri', [False, True])
    def test_rejects_invalid_tessellation_scaling_factor(self, tssf, ri):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()

        with pytest.raises(
                VeloxChemError,
                match='Tessellation sphere scaling factor must be finite and '
                      'positive'):
            scf_drv.update_settings({}, {
                'pressure': 1.0,
                'gostshyp_tssf': tssf,
            })

    @pytest.mark.parametrize('ri', [False, True])
    def test_custom_tco_tolerance_is_inherited_by_response(self, ri):

        xyz_string = """3

        O   -3.3278470    3.1951799   -0.0000000
        H   -4.2057717    2.7370843   -0.0000000
        H   -2.6643996    2.4600330   -0.0000000
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)
        custom_tco_tol = 2.5e-13

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.update_settings({}, {
            'pressure': 20000.0,
            'gostshyp_tco_tol': custom_tco_tol,
        })
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        if scf_drv.rank == mpi_master():
            assert scf_results['gostshyp_tco_tol'] == custom_tco_tol

        rsp_drv = LinearResponseEigenSolver()
        rsp_drv.ostream.mute()
        rsp_drv.nstates = 1
        rsp_results_not_used = rsp_drv.compute(mol, bas, scf_results)

        assert rsp_drv.gostshyp_tco_tol == custom_tco_tol
        assert rsp_drv._gostshyp_drv.tco_tol == custom_tco_tol

    @pytest.mark.parametrize('ri', [False, True])
    def test_scf_prints_negative_amplitude_exclusion_info(self, tmp_path, ri):

        xyz_string = """14

        C         -0.87931       -4.77039        0.03554
        C         -0.68902       -4.75793        1.56474
        C         -0.13597       -3.42950        2.03558
        O          0.03224       -3.23762        3.23105
        O          0.16603       -2.44864        1.15545
        N          0.21623       -5.83383        1.98004
        S         -1.69196       -6.29641       -0.54799
        H          0.09753       -4.64794       -0.48155
        H         -1.52147       -3.92199       -0.28390
        H         -1.68351       -4.90147        2.04299
        H          0.51499       -1.59560        1.45037
        H          0.30340       -5.83393        3.02290
        H         -0.18435       -6.76203        1.71757
        H         -0.50079       -6.84368       -1.06167
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        # write the SCF output to a file on the master rank only
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            ostream = OutputStream(tmp_path / 'scf_neg_amp.out')
        else:
            ostream = OutputStream(None)

        scf_drv = ScfRestrictedDriver(ostream=ostream)
        scf_drv.ri = ri
        scf_drv.pressure = 50
        scf_drv.pressure_units = 'GPa'
        scf_drv.gostshyp_discretization = 'swig'
        scf_drv.gostshyp_num_lebedev_points = 110
        scf_drv.gostshyp_tssf = 1.2
        scf_drv.gostshyp_r_ext = 0.0
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            ostream.close()
            output_text = (tmp_path / 'scf_neg_amp.out').read_text()

            # the count is reduced to the master rank only
            num_neg_amp = scf_drv._gostshyp_drv.num_neg_amp
            assert num_neg_amp > 0

            expected = (
                '*** GOSTSHYP information: A total number of '
                '{} grid points with negative amplitudes were '
                'excluded ***'.format(num_neg_amp))
            assert expected in output_text

            # the information is printed at the end of SCF, after the
            # convergence message
            assert (output_text.index(
                'grid points with negative amplitudes were excluded') >
                output_text.index('converged in'))

    def run_gostshyp(self,
                     mol,
                     basis_label,
                     pressure,
                     pressure_units,
                     discretization,
                     num_lebedev_points,
                     tssf,
                     r_ext,
                     ref_energy,
                     ref_dipole,
                     ref_grad,
                     ref_exc_energies,
                     ref_osc_strengths,
                     *,
                     ri):

        basis = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ri = ri
        scf_drv.ostream.mute()
        scf_drv.pressure = pressure
        scf_drv.pressure_units = pressure_units
        scf_drv.gostshyp_discretization = discretization
        scf_drv.gostshyp_num_lebedev_points = num_lebedev_points
        scf_drv.gostshyp_tssf = tssf
        scf_drv.gostshyp_r_ext = r_ext
        scf_results = scf_drv.compute(mol, basis)

        if scf_drv.rank == mpi_master():
            if ri:
                assert scf_results['ri_jk'] or scf_results['ri_coulomb']
                assert not scf_results['ri_fallback']

        # SCF energy, dipole moment and GOSTSHYP keywords are only available
        # on the master rank
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert abs(scf_results['scf_energy'] - ref_energy) < 1.0e-8
            assert abs(np.linalg.norm(scf_results['dipole_moment']) -
                       ref_dipole) < 1.0e-5

            # GOSTSHYP keywords are recorded in the SCF results
            assert scf_results['pressure'] == pytest.approx(pressure)
            assert scf_results['pressure_units'] == pressure_units
            assert scf_results['gostshyp_discretization'] == discretization
            assert scf_results['gostshyp_num_lebedev_points'] == num_lebedev_points
            assert scf_results['gostshyp_tssf'] == pytest.approx(tssf)
            assert scf_results['gostshyp_r_ext'] == pytest.approx(r_ext)
            assert scf_results['gostshyp_tco_tol'] == pytest.approx(1.0e-14)

        # Analytical gradient (the full gradient is available on all ranks)
        grad_drv = ScfGradientDriver(scf_drv)
        grad_drv.compute(mol, basis, scf_results)
        grad = grad_drv.get_gradient()
        assert np.max(np.abs(grad - ref_grad)) < 1.0e-6

        # Linear response (10 excited states)
        rsp_drv = LinearResponseEigenSolver()
        rsp_drv.ostream.mute()
        rsp_drv.nstates = 10
        rsp_results = rsp_drv.compute(mol, basis, scf_results)

        # Response results are only available on the master rank
        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert rsp_results['eigenvalues'].size == 10
            assert np.max(np.abs(rsp_results['eigenvalues'] -
                                 ref_exc_energies)) < 1.0e-5
            assert np.max(np.abs(rsp_results['oscillator_strengths'] -
                                 ref_osc_strengths)) < 1.0e-4

    @pytest.mark.parametrize('ri', [False, True])
    def test_water_631g(self, ri):

        xyz_string = """3

        O   -3.3278470    3.1951799   -0.0000000
        H   -4.2057717    2.7370843   -0.0000000
        H   -2.6643996    2.4600330   -0.0000000
        """
        mol = Molecule.read_xyz_string(xyz_string)

        ref_energy = {False: -75.8747887413, True: -75.8747792316}[ri]
        ref_dipole = 1.104423
        ref_grad = {
            False: np.array([
                [0.012209839697, 0.077470676071, 0.000000000000],
                [-0.041425742215, -0.031721548093, -0.000000000000],
                [0.029215902525, -0.045749127985, 0.000000000000],
            ]),
            True: np.array([
                [0.012213680610, 0.077491556479, 0.000000000000],
                [-0.041436630262, -0.031730409932, 0.000000000000],
                [0.029222949660, -0.045761146555, 0.000000000000],
            ]),
        }[ri]
        ref_exc_energies = {
            False: np.array([
                0.41324878, 0.48168018, 0.49322128, 0.57090623, 0.60811004,
                0.74283691, 1.14340276, 1.16001904, 1.18150362, 1.18449630,
            ]),
            True: np.array([
                0.41327169, 0.48169098, 0.49323423, 0.57091496, 0.60812787,
                0.74283140, 1.14340150, 1.15999633, 1.18149588, 1.18452782,
            ]),
        }[ri]
        ref_osc_strengths = {
            False: np.array([
                0.0127, 0.0000, 0.1506, 0.3356, 0.3081, 0.3124, 0.0766, 0.0214,
                0.0038, 0.0000,
            ]),
            True: np.array([
                0.012732, 0.000000, 0.150574, 0.335832, 0.307866, 0.312318,
                0.076645, 0.021417, 0.003772, 0.000000,
            ]),
        }[ri]

        self.run_gostshyp(mol, '6-31G', 50, 'GPa', 'swig', 110, 1.2, 0.0,
                          ref_energy, ref_dipole, ref_grad, ref_exc_energies,
                          ref_osc_strengths, ri=ri)

    @pytest.mark.parametrize('ri', [False, True])
    def test_ammonia_def2svpd(self, ri):

        xyz_string = """4

        N              1.257702000000         0.723638000000        -1.017351000000
        H              1.399549000000         0.070264000000        -0.214983000000
        H              1.950604000000         1.496433000000        -0.901218000000
        H              1.530872000000         0.202796000000        -1.880434000000
        """
        mol = Molecule.read_xyz_string(xyz_string)

        ref_energy = {False: -56.1279246511, True: -56.1279059014}[ri]
        ref_dipole = {False: 0.647947, True: 0.6479404}[ri]
        ref_grad = {
            False: np.array([
                [-0.036647008305, 0.013387125522, -0.000558714708],
                [0.003885093472, -0.022557151775, 0.028133684760],
                [0.023860048317, 0.027407288277, 0.003160775717],
                [0.008901866516, -0.018237262024, -0.030735745769],
            ]),
            True: np.array([
                [-0.036643201381, 0.013385751139, -0.000558491014],
                [0.003883067138, -0.022558378475, 0.028136171609],
                [0.023859842768, 0.027410692747, 0.003161009076],
                [0.008900291475, -0.018238065411, -0.030738689672],
            ]),
        }[ri]
        ref_exc_energies = {
            False: np.array([
                0.30911108, 0.38051211, 0.38054433, 0.45641532, 0.45706488,
                0.47120903, 0.48843886, 0.48848174, 0.51307889, 0.52254160,
            ]),
            True: np.array([
                0.30910397, 0.38051057, 0.38054279, 0.45642388, 0.45707349,
                0.47120209, 0.48844605, 0.48848893, 0.51294881, 0.52255181,
            ]),
        }[ri]
        ref_osc_strengths = {
            False: np.array([
                0.0903, 0.0333, 0.0334, 0.0000, 0.0000, 0.0002, 0.2311, 0.2306,
                0.2364, 0.1192,
            ]),
            True: np.array([
                0.090338, 0.033330, 0.033379, 0.000005, 0.000005, 0.000214,
                0.231140, 0.230660, 0.236300, 0.119098,
            ]),
        }[ri]

        self.run_gostshyp(mol, 'def2-svpd', 10, 'GPa', 'iswig', 110, 1.2,
                          0.25, ref_energy, ref_dipole, ref_grad,
                          ref_exc_energies, ref_osc_strengths, ri=ri)

    @pytest.mark.timeconsuming
    @pytest.mark.parametrize('ri', [False, True])
    def test_ethanol_def2svp(self, ri):

        xyz_string = """9

        C             -0.648457000000         1.139431000000         1.296608000000
        C              0.180279000000        -0.065100000000         0.867461000000
        O              0.908184000000        -0.568364000000         1.954230000000
        H             -1.197467000000         1.546810000000         0.421624000000
        H              0.014721000000        1.930578000000         1.706269000000
        H             -1.383109000000        0.839873000000         2.073787000000
        H              0.892270000000        0.248030000000         0.074931000000
        H             -0.483888000000       -0.852074000000         0.443428000000
        H              0.265522000000       -1.075860000000         2.515754000000
        """
        mol = Molecule.read_xyz_string(xyz_string)

        ref_energy = {False: -153.8555975233, True: -153.8555118170}[ri]
        ref_dipole = {False: 0.748203, True: 0.7482812}[ri]
        ref_grad = {
            False: np.array([
                [0.002151441718, -0.002347638328, 0.001949940303],
                [0.002508416427, -0.007421447257, 0.006257771148],
                [0.040482395173, 0.030054618932, -0.020039447690],
                [-0.006643660709, 0.004739948225, -0.015626276798],
                [0.013620462901, 0.013413389638, 0.005947484205],
                [-0.011697273670, -0.007782177499, 0.012364599460],
                [0.010568107155, 0.009566476962, -0.015006139020],
                [-0.012885406326, -0.012852668609, -0.003680279601],
                [-0.038104482672, -0.027370502067, 0.027832347996],
            ]),
            True: np.array([
                [0.002148285870, -0.002346912037, 0.001955363985],
                [0.002521028113, -0.007431141771, 0.006283884076],
                [0.040468277573, 0.030078544262, -0.020068639914],
                [-0.006639887256, 0.004735444206, -0.015621429201],
                [0.013618705102, 0.013404669123, 0.005947672810],
                [-0.011686931309, -0.007781672988, 0.012359496094],
                [0.010575453533, 0.009570264764, -0.015008409681],
                [-0.012890502430, -0.012860363331, -0.003678438280],
                [-0.038114429197, -0.027368832229, 0.027830500113],
            ]),
        }[ri]
        ref_exc_energies = {
            False: np.array([
                0.36943517, 0.43245564, 0.45910379, 0.46645888, 0.46948806,
                0.48061238, 0.48939530, 0.49717328, 0.50527400, 0.51345538,
            ]),
            True: np.array([
                0.36944277, 0.43246154, 0.45911565, 0.46647090, 0.46949758,
                0.48062605, 0.48940350, 0.49718977, 0.50528955, 0.51346836,
            ]),
        }[ri]
        ref_osc_strengths = {
            False: np.array([
                0.0034, 0.0257, 0.0345, 0.0995, 0.0356, 0.3733, 0.3376, 0.8105,
                0.3279, 0.0485,
            ]),
            True: np.array([
                0.003416, 0.025739, 0.034473, 0.099520, 0.035709, 0.373325,
                0.337779, 0.810478, 0.327873, 0.048488,
            ]),
        }[ri]

        self.run_gostshyp(mol, 'def2-svp', 25, 'GPa', 'iswig', 194, 1.2, 0.25,
                          ref_energy, ref_dipole, ref_grad, ref_exc_energies,
                          ref_osc_strengths, ri=ri)
