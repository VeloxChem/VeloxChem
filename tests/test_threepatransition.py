from mpi4py import MPI
from types import SimpleNamespace
import numpy as np
import pytest

from veloxchem.veloxchemlib import (mpi_master, hartree_in_ev,
                                    hartree_in_inverse_nm,
                                    fine_structure_constant)
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.threepatransitiondriver import ThreePATransitionDriver


@pytest.mark.solvers
class TestThreePATransition:

    @staticmethod
    def run_scf(xcfun_label):

        molecule_string = """
            O  0.0           0.0  0.0
            H   .7586020000  0.0  -.5042840000
            H   .7586020000  0.0   .5042840000
        """
        basis_set_label = 'def2-svpd'
        scf_conv_thresh = 1.0e-8

        molecule = Molecule.read_molecule_string(molecule_string)
        basis = MolecularBasis.read(molecule, basis_set_label, ostream=None)

        scf_drv = ScfRestrictedDriver()
        scf_drv.xcfun = xcfun_label
        scf_drv.conv_thresh = scf_conv_thresh
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(molecule, basis)

        return scf_results, molecule, basis

    @classmethod
    def compute_threepatransition(cls, xcfun_label):

        scf_results, molecule, ao_basis = cls.run_scf(xcfun_label)

        threepa_drv = ThreePATransitionDriver()
        threepa_drv.xcfun = xcfun_label
        threepa_drv.nstates = 2
        threepa_drv.conv_thresh = 1.0e-7
        threepa_drv.ostream.mute()

        return threepa_drv.compute(molecule, ao_basis, scf_results)

    @pytest.fixture(scope='class')
    @classmethod
    def threepa_results(cls):

        return cls.compute_threepatransition('hf')

    def get_synthetic_results(self):

        photon_energies = [0.123456789012345, 0.178901234567890]
        return {
            'photon_energies': photon_energies,
            '3pa_strengths': {
                'linear': dict(zip([-e for e in photon_energies], [1.0, 2.0])),
                'circular': dict(zip([-e for e in photon_energies],
                                     [10.0, 20.0])),
            },
        }

    def test_get_spectrum(self, threepa_results):

        if MPI.COMM_WORLD.Get_rank() != mpi_master():
            return

        results = threepa_results
        energies = np.asarray(results['photon_energies'])
        assert len(energies) == 2
        assert np.all(np.isfinite(energies))
        assert np.all(energies > 0.0)
        x_data_au = np.linspace(0.9 * energies.min(),
                                1.1 * energies.max(), 5).tolist()
        x_data_ev = [x * hartree_in_ev() for x in x_data_au]
        x_data_nm = [1.0 / (hartree_in_inverse_nm() * x) for x in x_data_au]

        spectrum_au = ThreePATransitionDriver.get_spectrum(
            results, x_data_au, 'au', 0.01, 'au')
        spectrum_ev = ThreePATransitionDriver.get_spectrum(
            results, x_data_ev, 'ev', 0.01 * hartree_in_ev(), 'ev')
        spectrum_nm = ThreePATransitionDriver.get_spectrum(
            results, x_data_nm, 'nm', 0.01, 'au')

        assert spectrum_au['x_label'] == 'Photon energy [a.u.]'
        assert spectrum_ev['x_label'] == 'Photon energy [eV]'
        assert spectrum_nm['x_label'] == 'Wavelength [nm]'
        assert spectrum_au['y_label'] == (
            r'3PA cross-section [$10^{-80}$ cm$^6$ s$^2$ photon$^{-2}$]')

        assert spectrum_au['x_data'] == pytest.approx(x_data_au)
        assert spectrum_ev['x_data'] == pytest.approx(x_data_ev)
        assert spectrum_nm['x_data'] == pytest.approx(x_data_nm)
        assert spectrum_au['y_data'] == pytest.approx(
            spectrum_ev['y_data'], rel=1.0e-12, abs=0.0)
        assert spectrum_au['y_data'] == pytest.approx(
            spectrum_nm['y_data'], rel=1.0e-12, abs=0.0)
        assert all(value > 0.0 for value in spectrum_au['y_data'])

    def test_get_spectrum_resonance(self):

        energy = 0.123456789012345
        strength = 2.5e6
        gamma = 0.004
        results = {
            'photon_energies': [energy],
            '3pa_strengths': {'linear': {-energy: strength}},
        }

        spectrum = ThreePATransitionDriver.get_spectrum(
            results, [0.0, energy], 'au', gamma, 'au')

        # a0**6 * t0**2 in cm**6 s**2, calculated from NIST CODATA 2018 constants.
        atomic_cross_section = 1.28480435108e-83
        sigma_au = ((4.0 * np.pi**3 / 3.0)
                    * fine_structure_constant()**3
                    * energy**3 * strength / gamma)
        expected = sigma_au * atomic_cross_section / 1.0e-80

        assert spectrum['y_data'][0] == 0.0
        assert spectrum['y_data'][1] == pytest.approx(
            expected, rel=1.0e-10, abs=0.0)

    def test_get_spectrum_linewidth(self):

        energy = 0.123456789012345
        gamma = 0.004
        results = {
            'photon_energies': [energy],
            '3pa_strengths': {'linear': {-energy: 1.0}},
        }
        x_data = np.array([energy - gamma / 3.0, energy,
                           energy + gamma / 3.0])
        spectrum = ThreePATransitionDriver.get_spectrum(
            results, x_data, 'au', gamma, 'au')

        # Check the Lorentzian width after removing the omega**3 factor.
        profile = np.array(spectrum['y_data']) / x_data**3
        assert profile / profile[1] == pytest.approx(
            [0.5, 1.0, 0.5], rel=1.0e-12, abs=0.0)

    def test_get_spectrum_polarization(self):

        results = self.get_synthetic_results()
        x_data = [0.12, 0.15, 0.18]
        spectrum_linear = ThreePATransitionDriver.get_spectrum(
            results, x_data, 'au', 0.01, 'au', polarization='linear')
        spectrum_circular = ThreePATransitionDriver.get_spectrum(
            results, x_data, 'au', 0.01, 'au', polarization='circular')

        assert spectrum_circular['y_data'] == pytest.approx(
            [10.0 * value for value in spectrum_linear['y_data']],
            rel=1.0e-12, abs=0.0)

    def test_get_spectrum_additivity(self):

        results = self.get_synthetic_results()
        x_data = [0.12, 0.15, 0.18]
        spectrum = ThreePATransitionDriver.get_spectrum(
            results, x_data, 'au', 0.01, 'au')
        expected = np.zeros(len(x_data))
        for energy in results['photon_energies']:
            single_state = {
                'photon_energies': [energy],
                '3pa_strengths': {
                    'linear': {-energy: results['3pa_strengths']['linear'][-energy]},
                },
            }
            contribution = ThreePATransitionDriver.get_spectrum(
                single_state, x_data, 'au', 0.01, 'au')
            expected += contribution['y_data']

        assert spectrum['y_data'] == pytest.approx(
            expected, rel=1.0e-12, abs=0.0)

    @pytest.mark.parametrize('x_unit', ['au', 'ev', 'nm'])
    @pytest.mark.parametrize('polarization', ['linear', 'circular'])
    @pytest.mark.parametrize('show_sticks', [True, False])
    def test_plot_spectrum(self, threepa_results, x_unit, polarization,
                           show_sticks):

        pytest.importorskip('matplotlib')
        import matplotlib.pyplot as plt

        if MPI.COMM_WORLD.Get_rank() != mpi_master():
            return

        results = threepa_results
        driver = ThreePATransitionDriver()
        driver.ostream.mute()
        fig, ax = plt.subplots()

        try:
            returned_ax = driver.plot_spectrum(
                results, x_unit=x_unit, broadening_value=0.124,
                broadening_unit='ev', polarization=polarization,
                ax=ax, show_sticks=show_sticks)

            assert returned_ax is ax
            assert len(ax.lines) == 1
            x_data = np.asarray(ax.lines[0].get_xdata())
            assert np.all(np.isfinite(x_data))
            assert np.all(np.diff(x_data) > 0.0)

            expected = ThreePATransitionDriver.get_spectrum(
                results, x_data, x_unit, 0.124, 'ev', polarization)
            assert ax.get_xlabel() == expected['x_label']
            assert ax.get_ylabel() == expected['y_label']
            assert ax.lines[0].get_ydata() == pytest.approx(
                expected['y_data'], rel=1.0e-12, abs=0.0)

            assert len(fig.axes) == (2 if show_sticks else 1)
            if show_sticks:
                stick_ax = fig.axes[1]
                assert stick_ax.get_ylabel() == (
                    f'3PA strengths ({polarization}) [a.u.]')
                energies = np.asarray(results['photon_energies'])
                if x_unit == 'ev':
                    positions = energies * hartree_in_ev()
                elif x_unit == 'nm':
                    positions = 1.0 / (hartree_in_inverse_nm() * energies)
                else:
                    positions = energies

                strengths = [results['3pa_strengths'][polarization][-energy]
                             for energy in energies]
                assert [bar.get_height() for bar in stick_ax.patches] == (
                    pytest.approx(strengths))
                assert [bar.get_x() + bar.get_width() / 2.0
                        for bar in stick_ax.patches] == pytest.approx(positions)
                assert np.all(positions >= x_data.min())
                assert np.all(positions <= x_data.max())
        finally:
            plt.close(fig)

    def test_plot_spectrum_without_ax_returns_none(self, threepa_results):

        pytest.importorskip('matplotlib')
        import matplotlib.pyplot as plt

        if MPI.COMM_WORLD.Get_rank() != mpi_master():
            return

        driver = ThreePATransitionDriver()
        driver.ostream.mute()
        try:
            assert driver.plot_spectrum(threepa_results) is None
        finally:
            plt.close('all')

    def test_plot_spectrum_nonmaster_returns_none(self):

        driver = SimpleNamespace(rank=mpi_master() + 1)
        assert ThreePATransitionDriver.plot_spectrum(driver, None) is None
