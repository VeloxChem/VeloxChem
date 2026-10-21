from mpi4py import MPI
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.outputstream import OutputStream
from veloxchem.errorhandler import VeloxChemError
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.cubicresponsedriver import CubicResponseDriver
from veloxchem.excitedstatemomentdriver import ExcitedStateMomentDriver
from veloxchem.quadraticresponsedriver import QuadraticResponseDriver
from veloxchem.shgdriver import ShgDriver
from veloxchem.thgdriver import ThgDriver
from veloxchem.thgreddriver import ThgReducedDriver
from veloxchem.threepatransitiondriver import ThreePATransitionDriver
from veloxchem.tpafulldriver import TpaFullDriver
from veloxchem.tpareddriver import TpaReducedDriver
from veloxchem.tpatransitiondriver import TpaTransitionDriver


@pytest.mark.solvers
class TestNonlinearEnvironmentGap:

    NONLINEAR_DRIVERS = [
        ShgDriver,
        ThgDriver,
        QuadraticResponseDriver,
        CubicResponseDriver,
        ExcitedStateMomentDriver,
        TpaFullDriver,
        TpaReducedDriver,
        TpaTransitionDriver,
        ThgReducedDriver,
        ThreePATransitionDriver,
    ]

    def setup_method(self):

        xyz_string = """3
        xyz
        O   0.000000    0.000000    0.000000
        H   0.758602    0.000000   -0.504284
        H   0.758602    0.000000    0.504284
        """

        self.molecule = Molecule.read_xyz_string(xyz_string)
        self.basis = MolecularBasis.read(self.molecule, 'sto-3g')

    def run_scf(self, **env_settings):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        for key, value in env_settings.items():
            setattr(scf_drv, key, value)
        return scf_drv.compute(self.molecule, self.basis)

    def assert_compute_rejected(self, driver_cls, scf_results, errmsg):

        drv = driver_cls(MPI.COMM_WORLD, OutputStream())
        drv.ostream.mute()
        # drivers that validate their own settings before the SCF sanity
        # check need valid defaults set explicitly
        for attr in ('a_component', 'b_component', 'c_component',
                     'd_component'):
            if hasattr(drv, attr):
                setattr(drv, attr, 'x')
        with pytest.raises(VeloxChemError, match=errmsg):
            drv.compute(self.molecule, self.basis, scf_results)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='pytest.raises only valid in serial')
    def test_compute_rejects_inherited_solvation_model(self):

        scf_results = self.run_scf(solvation_model='cpcm',
                                   cpcm_epsilon=10.0)

        errmsg = ("NonlinearSolver: The 'solvation_model' keyword is not "
                  "supported in nonlinear response calculation.")
        for driver_cls in self.NONLINEAR_DRIVERS:
            self.assert_compute_rejected(driver_cls, scf_results, errmsg)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='pytest.raises only valid in serial')
    def test_compute_rejects_inherited_pressure(self):

        scf_results = self.run_scf(pressure=20000.0)

        errmsg = ("NonlinearSolver: The 'pressure' keyword is not "
                  "supported in nonlinear response calculation.")
        for driver_cls in self.NONLINEAR_DRIVERS:
            self.assert_compute_rejected(driver_cls, scf_results, errmsg)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='pytest.raises only valid in serial')
    def test_compute_rejects_inherited_potfile(self):

        # the check fires before any tensor is accessed, so an SCF run is
        # not needed
        scf_results = {'potfile': 'fake_pe.pot'}

        errmsg = ("NonlinearSolver: The 'potfile' keyword is not "
                  "supported in nonlinear response calculation.")
        for driver_cls in self.NONLINEAR_DRIVERS:
            self.assert_compute_rejected(driver_cls, scf_results, errmsg)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='pytest.raises only valid in serial')
    def test_update_settings_rejects_environment_settings(self):

        settings = [
            ('solvation_model', 'cpcm',
             "NonlinearSolver: The 'solvation_model' keyword is not "
             "supported in nonlinear response calculation."),
            ('pressure', 20000.0,
             "NonlinearSolver: The 'pressure' keyword is not "
             "supported in nonlinear response calculation."),
            ('potfile', 'fake_pe.pot',
             "NonlinearSolver: The 'potfile' keyword is not "
             "supported in nonlinear response calculation."),
        ]
        # representative drivers to verify the super() chaining
        for driver_cls in [ThgDriver, ShgDriver, TpaFullDriver]:
            for attr, value, errmsg in settings:
                drv = driver_cls(MPI.COMM_WORLD, OutputStream())
                drv.ostream.mute()
                setattr(drv, attr, value)
                with pytest.raises(VeloxChemError, match=errmsg):
                    drv.update_settings({}, {})

    def test_electric_field_is_supported(self):

        scf_results = self.run_scf(electric_field=(0.0, 0.0, 1.0e-4))

        qrf_drv = QuadraticResponseDriver()
        qrf_drv.ostream.mute()
        qrf_drv.update_settings({
            'b_frequencies': [0.0656],
            'c_frequencies': [0.0445],
            'damping': 0.0,
            'a_operator': 'electric dipole',
            'b_operator': 'electric dipole',
            'c_operator': 'electric dipole',
            'a_component': 'x',
            'b_component': 'x',
            'c_component': 'x',
        })
        qrf_result = qrf_drv.compute(self.molecule, self.basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert ('qrf', 0.0656, 0.0445) in qrf_result
