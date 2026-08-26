from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.errorhandler import VeloxChemError
from veloxchem.sanitychecks import (environment_compatibility_sanity_check,
                                    scf_results_sanity_check)


def one_point_charge():
    """
    A minimal point-charge array (6 x 1).
    """

    return np.zeros((6, 1))


def make_recording_output():
    """Creates a minimal output stream that records warnings."""

    warnings = []
    ostream = SimpleNamespace(print_warning=warnings.append,
                              flush=lambda: None)
    return ostream, warnings


@pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                    reason='warning output is checked on the master process')
@pytest.mark.parametrize('environment, value', [
    ('pressure', 100.0),
    ('potfile', 'embedding.json'),
])
def test_warns_about_environment_absent_from_scf(environment, value):

    driver = LinearResponseEigenSolver()
    driver._ostream, warnings = make_recording_output()
    driver.xcfun = 'BLYP'
    setattr(driver, environment, value)

    scf_results_sanity_check(driver, {'xcfun': 'BLYP'})

    assert warnings == [
        'Environment settings active in the current calculation but absent '
        f'from the SCF reference: {environment}. The SCF reference orbitals '
        'and density are not relaxed for these settings.'
    ]


@pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                    reason='warning output is checked on the master process')
@pytest.mark.parametrize('environment, value', [
    ('pressure', 100.0),
    ('potfile', 'embedding.json'),
])
def test_does_not_warn_about_environment_inherited_from_scf(environment,
                                                            value):

    driver = LinearResponseEigenSolver()
    driver._ostream, warnings = make_recording_output()
    driver.xcfun = 'BLYP'
    setattr(driver, environment, value)

    scf_results = {'xcfun': 'BLYP', environment: value}
    if environment == 'pressure':
        scf_results['pressure_units'] = 'MPa'

    scf_results_sanity_check(driver, scf_results)

    assert warnings == []


@pytest.mark.solvers
class TestEnvironmentCompatibility:
    """
    Tests the pairwise compatibility of the five environment settings
    (potfile / PE, solvation_model / CPCM, pressure / GOSTSHYP,
    electric_field, point_charges).

    Policy: each setting can be used individually; the supported
    combinations are 'electric_field' with 'point_charges' and
    'electric_field' with 'pressure' (GOSTSHYP). All other
    combinations are rejected, in 'update_settings' and 'compute'.
    """

    @staticmethod
    def get_water_and_basis():

        xyz_string = """3
        xyz
        O   -0.1858140  -1.1749469   0.7662596
        H   -0.1285513  -0.8984365   1.6808606
        H   -0.0582782  -0.3702550   0.2638279
        """

        molecule = Molecule.read_xyz_string(xyz_string)
        basis = MolecularBasis.read(molecule, "sto-3g", ostream=None)

        return molecule, basis

    @staticmethod
    def make_dummy(**kwargs):
        """
        Builds a minimal driver-like object carrying the five environment
        attributes, for unit-testing the sanity check helper.
        """

        defaults = {
            '_pe': False,
            'potfile': None,
            'solvation_model': None,
            'pressure': 0.0,
            '_gostshyp': False,
            'electric_field': None,
            'point_charges': None,
        }
        defaults.update(kwargs)
        return SimpleNamespace(**defaults)

    # ------------------------------------------------------------------
    # unit tests of the sanity check helper (all nine rejected pairs)
    # ------------------------------------------------------------------

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    @pytest.mark.parametrize(
        'settings_a, settings_b',
        [
            # PE with each of the other four
            ({'_pe': True}, {'solvation_model': 'cpcm'}),
            ({'_pe': True}, {'pressure': 100.0}),
            ({'_pe': True}, {'electric_field': [0.0, 0.0, 1.0e-4]}),
            ({'_pe': True}, {'point_charges': one_point_charge()}),
            # CPCM with each of the other four
            ({'solvation_model': 'cpcm'}, {'pressure': 100.0}),
            ({'solvation_model': 'cpcm'}, {'electric_field': [0.0, 0.0, 1.0e-4]}),
            ({'solvation_model': 'cpcm'}, {'point_charges': one_point_charge()}),
            # GOSTSHYP with the static settings
            ({'pressure': 100.0}, {'point_charges': one_point_charge()}),
        ])
    def test_helper_rejects_incompatible_pairs(self, settings_a, settings_b):

        dummy = self.make_dummy(**settings_a, **settings_b)

        with pytest.raises(VeloxChemError,
                           match='Incompatible environment settings'):
            environment_compatibility_sanity_check(dummy)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_helper_rejects_three_settings(self):

        # only the exact combinations field + point_charges and
        # field + pressure are allowed
        dummy = self.make_dummy(pressure=100.0,
                                electric_field=[0.0, 0.0, 1.0e-4],
                                point_charges=one_point_charge())

        with pytest.raises(VeloxChemError,
                           match='Incompatible environment settings'):
            environment_compatibility_sanity_check(dummy)

    def test_helper_allows_electric_field_with_point_charges(self):

        dummy = self.make_dummy(electric_field=[0.0, 0.0, 1.0e-4],
                                point_charges=one_point_charge())

        environment_compatibility_sanity_check(dummy)

    def test_helper_allows_electric_field_with_pressure(self):

        dummy = self.make_dummy(pressure=100.0,
                                electric_field=[0.0, 0.0, 1.0e-4])

        environment_compatibility_sanity_check(dummy)

    def test_helper_allows_each_setting_individually(self):

        for settings in [
                {'_pe': True},
                {'potfile': 'pot.json'},
                {'solvation_model': 'smd'},
                {'pressure': 100.0},
                {'_gostshyp': True},
                {'electric_field': [0.0, 0.0, 1.0e-4]},
                {'point_charges': one_point_charge()},
        ]:
            environment_compatibility_sanity_check(self.make_dummy(**settings))

    def test_helper_allows_no_environment(self):

        environment_compatibility_sanity_check(self.make_dummy())

    def test_helper_treats_none_pressure_as_inactive(self):

        # gostshyp_sanity_check normalizes pressure None to 0.0
        dummy = self.make_dummy(pressure=None,
                                electric_field=[0.0, 0.0, 1.0e-4])

        environment_compatibility_sanity_check(dummy)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_helper_validates_electric_field_length(self):

        dummy = self.make_dummy(electric_field=[1.0e-4, 2.0e-4])

        with pytest.raises(VeloxChemError,
                           match="Expecting 3 values in 'electric field'"):
            environment_compatibility_sanity_check(dummy)

    # ------------------------------------------------------------------
    # integration tests through ScfRestrictedDriver.compute
    # ------------------------------------------------------------------

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_cpcm_with_electric_field(self):

        molecule, basis = self.get_water_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.solvation_model = 'cpcm'
        scf_drv.electric_field = [0.0, 0.0, 1.0e-4]

        with pytest.raises(VeloxChemError,
                           match='Incompatible environment settings'):
            scf_drv.compute(molecule, basis)

    def test_compute_allows_gostshyp_with_electric_field(self):

        # the supported combination pressure + electric field runs to
        # convergence: the field is a static one-electron term added to
        # the density-dependent GOSTSHYP potential within the SCF
        # iterations
        molecule, basis = self.get_water_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.pressure = 100.0
        scf_drv.electric_field = [0.0, 0.0, 1.0e-4]

        scf_results = scf_drv.compute(molecule, basis)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert scf_results is not None
            assert scf_drv.is_converged
            assert np.isfinite(scf_results['scf_energy'])
            assert 'pressure' in scf_results
            assert 'electric_field' in scf_results

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_gostshyp_with_point_charges(self):

        molecule, basis = self.get_water_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.pressure = 100.0
        scf_drv.point_charges = one_point_charge()

        with pytest.raises(VeloxChemError,
                           match='Incompatible environment settings'):
            scf_drv.compute(molecule, basis)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_pe_with_gostshyp(self, tmp_path):

        pytest.importorskip('pyframe')
        molecule, basis = self.get_water_and_basis()
        potfile = tmp_path / 'dummy_pe.json'
        potfile.write_text('{}')

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.potfile = str(potfile)
        scf_drv.pressure = 100.0

        with pytest.raises(VeloxChemError,
                           match='Incompatible environment settings'):
            scf_drv.compute(molecule, basis)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_pe_with_electric_field(self, tmp_path):

        # the pre-existing rejection only ran in update_settings; compute
        # must also reject the pair when attributes are set directly
        pytest.importorskip('pyframe')
        molecule, basis = self.get_water_and_basis()
        potfile = tmp_path / 'dummy_pe.json'
        potfile.write_text('{}')

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.potfile = str(potfile)
        scf_drv.electric_field = [0.0, 0.0, 1.0e-4]

        with pytest.raises(VeloxChemError,
                           match='Incompatible environment settings'):
            scf_drv.compute(molecule, basis)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_compute_rejects_pe_with_point_charges(self, tmp_path):

        # same as above: the pair must be rejected in compute as well
        pytest.importorskip('pyframe')
        molecule, basis = self.get_water_and_basis()
        potfile = tmp_path / 'dummy_pe.json'
        potfile.write_text('{}')

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.potfile = str(potfile)
        scf_drv.point_charges = one_point_charge()

        with pytest.raises(VeloxChemError,
                           match='Incompatible environment settings'):
            scf_drv.compute(molecule, basis)

    def test_update_settings_allows_gostshyp_with_electric_field(self):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()

        scf_drv.update_settings(
            {}, {'pressure': 100.0,
                 'electric_field': [0.0, 0.0, 1.0e-4]})

        assert scf_drv.pressure == 100.0
        assert scf_drv.electric_field == (0.0, 0.0, 1.0e-4)

    def test_compute_allows_electric_field_with_point_charges(self):

        # the only supported combination runs to convergence
        molecule, basis = self.get_water_and_basis()

        point_charges = np.zeros((6, 1))
        point_charges[1, 0] = 6.0
        point_charges[3, 0] = 1.0

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.electric_field = [0.0, 0.0, 1.0e-4]
        scf_drv.point_charges = point_charges

        scf_results = scf_drv.compute(molecule, basis)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert scf_results is not None
            assert scf_drv.is_converged
            assert np.isfinite(scf_results['scf_energy'])
            assert 'electric_field' in scf_results
            assert 'point_charges' in scf_results

    # ------------------------------------------------------------------
    # integration tests through linear response drivers
    # (explicitly set response environments survive the SCF inheritance
    # when the SCF results do not carry those keys, so the response
    # drivers must enforce the same compatibility rule themselves)
    # ------------------------------------------------------------------

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_lr_update_settings_rejects_cpcm_with_electric_field(self):

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()

        with pytest.raises(VeloxChemError,
                           match='Incompatible environment settings'):
            lr_drv.update_settings(
                {}, {'solvation_model': 'cpcm',
                     'electric_field': [0.0, 0.0, 1.0e-4]})

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_lr_compute_rejects_cpcm_with_electric_field(self):

        molecule, basis = self.get_water_and_basis()

        # empty scf results: no environment is inherited from SCF, so the
        # explicitly set response environment is what is being checked
        scf_results = {}

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_drv.solvation_model = 'cpcm'
        lr_drv.electric_field = [0.0, 0.0, 1.0e-4]

        with pytest.raises(VeloxChemError,
                           match='Incompatible environment settings'):
            lr_drv.compute(molecule, basis, scf_results)

    def test_lr_compute_allows_gostshyp_with_electric_field(self):

        # the response driver inherits both the pressure and the electric
        # field from the SCF results; the combined environment is
        # supported and the response calculation runs to completion
        molecule, basis = self.get_water_and_basis()

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()
        scf_drv.pressure = 100.0
        scf_drv.electric_field = [0.0, 0.0, 1.0e-4]
        scf_results = scf_drv.compute(molecule, basis)

        lr_drv = LinearResponseEigenSolver()
        lr_drv.ostream.mute()
        lr_results = lr_drv.compute(molecule, basis, scf_results)

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            assert 'eigenvalues' in lr_results
            assert np.all(np.isfinite(lr_results['eigenvalues']))
