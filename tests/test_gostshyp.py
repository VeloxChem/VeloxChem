import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfgradientdriver import ScfGradientDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.tessellation import TessellationDriver
from veloxchem.errorhandler import VeloxChemError


def test_write_grid_to_file(tmp_path):
    """Checks that the grid writer uses its mandatory XYZ filename."""

    tessellation_drv = TessellationDriver()
    vdw_surface = np.zeros((15, 2))
    vdw_surface[:3, 0] = [1.0, 2.0, 3.0]
    vdw_surface[:3, 1] = [4.0, 5.0, 6.0]
    xyz_filename = tmp_path / 'tessellation.xyz'

    tessellation_drv.write_grid_to_file(vdw_surface, xyz_filename)

    lines = xyz_filename.read_text().splitlines()
    assert lines[0].strip() == '2'
    assert lines[2].split() == ['x', '1.0', '2.0', '3.0']
    assert lines[3].split() == ['x', '4.0', '5.0', '6.0']


@pytest.mark.solvers
class TestGostshyp:
    """GOSTSHYP hydrostatic pressure tests: SCF, analytical gradient and
    linear response for small molecules.
    """

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    @pytest.mark.parametrize('discretization', ['iwig', 'becke'])
    def test_rejects_invalid_discretization(self, discretization):

        xyz_string = """3

        O   -3.3278470    3.1951799   -0.0000000
        H   -4.2057717    2.7370843   -0.0000000
        H   -2.6643996    2.4600330   -0.0000000
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.pressure = 1.0
        scf_drv.gostshyp_discretization = discretization

        with pytest.raises(VeloxChemError,
                           match='GOSTSHYP: Invalid discretization'):
            scf_results_not_used = scf_drv.compute(mol, bas)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    @pytest.mark.parametrize('tco_tol', [0.0, -1.0e-14])
    def test_rejects_nonpositive_tco_tolerance(self, tco_tol):

        scf_drv = ScfRestrictedDriver()
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
    def test_rejects_invalid_tessellation_scaling_factor(self, tssf):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()

        with pytest.raises(
                VeloxChemError,
                match='Tessellation sphere scaling factor must be finite and '
                      'positive'):
            scf_drv.update_settings({}, {
                'pressure': 1.0,
                'gostshyp_tssf': tssf,
            })

    def test_custom_tco_tolerance_is_inherited_by_response(self):

        xyz_string = """3

        O   -3.3278470    3.1951799   -0.0000000
        H   -4.2057717    2.7370843   -0.0000000
        H   -2.6643996    2.4600330   -0.0000000
        """
        mol = Molecule.read_xyz_string(xyz_string)
        bas = MolecularBasis.read(mol, 'sto-3g', ostream=None)
        custom_tco_tol = 2.5e-13

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.update_settings({}, {
            'pressure': 20000.0,
            'gostshyp_tco_tol': custom_tco_tol,
        })
        scf_results = scf_drv.compute(mol, bas)

        if scf_drv.rank == mpi_master():
            assert scf_results['gostshyp_tco_tol'] == custom_tco_tol

        rsp_drv = LinearResponseEigenSolver()
        rsp_drv.ostream.mute()
        rsp_drv.nstates = 1
        rsp_results_not_used = rsp_drv.compute(mol, bas, scf_results)

        assert rsp_drv.gostshyp_tco_tol == custom_tco_tol
        assert rsp_drv._gostshyp_drv.tco_tol == custom_tco_tol

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
                     ref_osc_strengths):

        basis = MolecularBasis.read(mol, basis_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.pressure = pressure
        scf_drv.pressure_units = pressure_units
        scf_drv.gostshyp_discretization = discretization
        scf_drv.gostshyp_num_lebedev_points = num_lebedev_points
        scf_drv.gostshyp_tssf = tssf
        scf_drv.gostshyp_r_ext = r_ext
        scf_results = scf_drv.compute(mol, basis)

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

    def test_water_631g(self):

        xyz_string = """3

        O   -3.3278470    3.1951799   -0.0000000
        H   -4.2057717    2.7370843   -0.0000000
        H   -2.6643996    2.4600330   -0.0000000
        """
        mol = Molecule.read_xyz_string(xyz_string)

        ref_energy = -75.8747887413
        ref_dipole = 1.104423
        ref_grad = np.array([
            [0.012209839697, 0.077470676071, 0.000000000000],
            [-0.041425742215, -0.031721548093, -0.000000000000],
            [0.029215902525, -0.045749127985, 0.000000000000],
        ])
        ref_exc_energies = np.array([
            0.41324878, 0.48168018, 0.49322128, 0.57090623, 0.60811004,
            0.74283691, 1.14340276, 1.16001904, 1.18150362, 1.18449630,
        ])
        ref_osc_strengths = np.array([
            0.0127, 0.0000, 0.1506, 0.3356, 0.3081, 0.3124, 0.0766, 0.0214,
            0.0038, 0.0000,
        ])

        self.run_gostshyp(mol, '6-31G', 50, 'GPa', 'swig', 110, 1.2, 0.0,
                          ref_energy, ref_dipole, ref_grad, ref_exc_energies,
                          ref_osc_strengths)

    def test_ammonia_def2svpd(self):

        xyz_string = """4

        N              1.257702000000         0.723638000000        -1.017351000000
        H              1.399549000000         0.070264000000        -0.214983000000
        H              1.950604000000         1.496433000000        -0.901218000000
        H              1.530872000000         0.202796000000        -1.880434000000
        """
        mol = Molecule.read_xyz_string(xyz_string)

        ref_energy = -56.1279246511
        ref_dipole = 0.647947
        ref_grad = np.array([
            [-0.036647008305, 0.013387125522, -0.000558714708],
            [0.003885093472, -0.022557151775, 0.028133684760],
            [0.023860048317, 0.027407288277, 0.003160775717],
            [0.008901866516, -0.018237262024, -0.030735745769],
        ])
        ref_exc_energies = np.array([
            0.30911108, 0.38051211, 0.38054433, 0.45641532, 0.45706488,
            0.47120903, 0.48843886, 0.48848174, 0.51307889, 0.52254160,
        ])
        ref_osc_strengths = np.array([
            0.0903, 0.0333, 0.0334, 0.0000, 0.0000, 0.0002, 0.2311, 0.2306,
            0.2364, 0.1192,
        ])

        self.run_gostshyp(mol, 'def2-svpd', 10, 'GPa', 'iswig', 110, 1.2,
                          0.25, ref_energy, ref_dipole, ref_grad,
                          ref_exc_energies, ref_osc_strengths)

    @pytest.mark.timeconsuming
    def test_ethanol_def2svp(self):

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

        ref_energy = -153.8555975233
        ref_dipole = 0.748203
        ref_grad = np.array([
            [0.002151441718, -0.002347638328, 0.001949940303],
            [0.002508416427, -0.007421447257, 0.006257771148],
            [0.040482395173, 0.030054618932, -0.020039447690],
            [-0.006643660709, 0.004739948225, -0.015626276798],
            [0.013620462901, 0.013413389638, 0.005947484205],
            [-0.011697273670, -0.007782177499, 0.012364599460],
            [0.010568107155, 0.009566476962, -0.015006139020],
            [-0.012885406326, -0.012852668609, -0.003680279601],
            [-0.038104482672, -0.027370502067, 0.027832347996],
        ])
        ref_exc_energies = np.array([
            0.36943517, 0.43245564, 0.45910379, 0.46645888, 0.46948806,
            0.48061238, 0.48939530, 0.49717328, 0.50527400, 0.51345538,
        ])
        ref_osc_strengths = np.array([
            0.0034, 0.0257, 0.0345, 0.0995, 0.0356, 0.3733, 0.3376, 0.8105,
            0.3279, 0.0485,
        ])

        self.run_gostshyp(mol, 'def2-svp', 25, 'GPa', 'iswig', 194, 1.2, 0.25,
                          ref_energy, ref_dipole, ref_grad, ref_exc_energies,
                          ref_osc_strengths)
