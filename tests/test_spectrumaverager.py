import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.cppsolver import ComplexResponseSolver
from veloxchem.ensembledriver import EnsembleDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver
from veloxchem.spectrumaverager import SpectrumAverager
from veloxchem.veloxchemlib import (extinction_coefficient_from_beta,
                                    fine_structure_constant, mpi_master,
                                    rotatory_strength_in_cgs)


def _make_tddft_results():
    return [
        (6, {
            'eigenvalues': np.array([0.21, 0.27]),
            'oscillator_strengths': np.array([0.8, 0.2]),
            'rotatory_strengths': np.array([30.0, -10.0]),
        }),
        (8, {
            'eigenvalues': np.array([0.22, 0.28]),
            'oscillator_strengths': np.array([0.7, 0.3]),
            'rotatory_strengths': np.array([20.0, -5.0]),
        }),
    ]


def _make_cpp_results():
    frequencies = np.array([0.20, 0.22])

    def make_result():
        response_functions = {}
        for w in frequencies:
            for component in ('x', 'y', 'z'):
                response_functions[(component, component, w)] = (
                    complex(0.0, -1.0e-3))

        return {
            'frequencies': frequencies,
            'response_functions': response_functions,
        }

    return [(6, make_result()), (8, make_result())]


def _make_driver():
    driver = EnsembleDriver()
    driver.ostream.mute()
    return driver


@pytest.mark.solvers
class TestSpectrumProperties:

    @pytest.mark.parametrize(
        ('property_name', 'expected'),
        [
            ('absorption', ('absorption', False)),
            ('uv-vis', ('absorption', False)),
            ('absorption (cpp)', ('absorption', True)),
            ('absorption(cpp)', ('absorption', True)),
            ('linear absorption cross-section (cpp)',
             ('absorption', True)),
            ('ecd', ('ecd', False)),
            ('circular dichroism', ('ecd', False)),
            ('ecd (cpp)', ('ecd', True)),
            ('circular dichroism spectrum(cpp)', ('ecd', True)),
        ],
    )
    def test_parse_property_type(self, property_name, expected):

        assert EnsembleDriver._parse_property_type(
            property_name) == expected

    def test_build_property_driver(self):

        driver = _make_driver()

        assert isinstance(
            driver._build_property_driver({
                'property': 'absorption',
                'nstates': 2,
            }), LinearResponseEigenSolver)

        assert isinstance(
            driver._build_property_driver({
                'property': 'ecd',
                'nstates': 2,
            }), LinearResponseEigenSolver)

        assert isinstance(
            driver._build_property_driver({
                'property': 'absorption(cpp)',
                'frequencies': '0.20 - 0.22 (0.01)',
            }), ComplexResponseSolver)

        assert isinstance(
            driver._build_property_driver({
                'property': 'ecd(cpp)',
                'frequencies': '0.20 - 0.22 (0.01)',
            }), ComplexResponseSolver)

        with pytest.raises(ValueError, match='frequencies'):
            driver._build_property_driver({
                'property': 'absorption(cpp)',
                'nstates': 2,
            })

        with pytest.raises(ValueError, match='nstates'):
            driver._build_property_driver({
                'property': 'absorption',
                'frequencies': '0.20 - 0.22 (0.01)',
            })

    def test_get_absorption_spectrum_formula(self):

        driver = LinearResponseEigenSolver()
        driver.ostream.mute()

        rsp_results = {
            'eigenvalues': np.array([0.20]),
            'oscillator_strengths': np.array([1.0]),
        }
        x_data = [0.15, 0.20, 0.25]
        b_value = 0.02

        spectrum = driver.get_absorption_spectrum(
            rsp_results, x_data, 'au', b_value, 'au')

        sigma_factor = 2.0 * np.pi * fine_structure_constant()
        expected = []
        for x in x_data:
            y = 0.0
            for e, f in zip(rsp_results['eigenvalues'],
                            rsp_results['oscillator_strengths']):
                b_factor = b_value / ((e - x)**2 + b_value**2)
                y += sigma_factor * (x / e) * b_factor * f
            expected.append(y)

        np.testing.assert_allclose(spectrum['y_data'], expected)

    def test_get_ecd_spectrum_formula(self):

        driver = LinearResponseEigenSolver()
        driver.ostream.mute()

        rsp_results = {
            'eigenvalues': np.array([0.20, 0.30]),
            'rotatory_strengths': np.array([10.0, -5.0]),
        }
        x_data = [0.15, 0.20, 0.25]
        b_value = 0.02

        spectrum = driver.get_ecd_spectrum(
            rsp_results, x_data, 'au', b_value, 'au')

        delta_eps_factor = extinction_coefficient_from_beta() / 3.0
        rotatory_strengths = (rsp_results['rotatory_strengths'] /
                              rotatory_strength_in_cgs())
        expected = []
        for x in x_data:
            y = 0.0
            for e, r in zip(rsp_results['eigenvalues'],
                            rotatory_strengths):
                b_factor = b_value / ((e - x)**2 + b_value**2)
                y += delta_eps_factor * b_factor * x * r
            expected.append(y)

        np.testing.assert_allclose(spectrum['y_data'], expected)


@pytest.mark.solvers
class TestSpectrumAveragerECD:

    def test_compute_ecd_matches_manual_average(self):

        rsp_all = _make_tddft_results()

        averager = SpectrumAverager()
        averager.broadening_type = 'lorentzian'
        averager.broadening_value_ev = 0.10

        spectrum = averager.compute_ecd(rsp_all)

        lr_driver = LinearResponseEigenSolver()
        lr_driver.ostream.mute()

        manual_spectra = []
        for _, rsp_results in rsp_all:
            result = lr_driver.get_ecd_spectrum(
                rsp_results,
                x_data=spectrum['wavelength_nm'],
                x_unit='nm',
                b_value=averager.broadening_value_ev,
                b_unit='ev',
            )
            manual_spectra.append(result['y_data'])

        manual_average = np.mean(manual_spectra, axis=0)

        np.testing.assert_allclose(spectrum['delta_eps_mean_wl'],
                                   manual_average,
                                   rtol=1.0e-12,
                                   atol=1.0e-12)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip plot for multiple MPI processes')
    def test_plot_ecd_spectra(self, tmp_path):

        plt = pytest.importorskip('matplotlib.pyplot')

        averager = SpectrumAverager()
        averager.broadening_type = 'lorentzian'
        averager.broadening_value_ev = 0.10

        csv_file = tmp_path / 'ecd.csv'
        ax = averager.plot_ecd_spectra(
            _make_tddft_results(),
            show_individual=True,
            show_std=True,
            save_averaged_spectra=True,
            averaged_spectra_filename=csv_file,
        )

        assert ax is not None
        assert csv_file.is_file()

        header = csv_file.read_text().splitlines()[0]
        assert header == (
            'Wavelength [nm],Delta epsilon [L mol^-1 cm^-1]')

        plt.close(ax.figure)


@pytest.mark.solvers
class TestEnsembleDriverSpectrumSaving:

    def test_save_tddft_spectra(self, tmp_path):

        driver = _make_driver()
        rsp_all = _make_tddft_results()

        driver.rsp_dict = {'property': 'absorption'}
        absorption_file = tmp_path / 'absorption.csv'
        absorption = driver.save_uv_vis_spectra(
            {'rsp_all': rsp_all},
            broadening_type='lorentzian',
            broadening_value_ev=0.10,
            averaged_spectra_filename=absorption_file,
        )

        driver.rsp_dict = {'property': 'ecd'}
        ecd_file = tmp_path / 'ecd.csv'
        ecd = driver.save_uv_vis_spectra(
            {'rsp_all': rsp_all},
            broadening_type='lorentzian',
            broadening_value_ev=0.10,
            averaged_spectra_filename=ecd_file,
        )

        if driver.rank == mpi_master():
            absorption_data = np.loadtxt(absorption_file,
                                         delimiter=',',
                                         skiprows=1)
            ecd_data = np.loadtxt(ecd_file, delimiter=',', skiprows=1)

            np.testing.assert_allclose(absorption_data[:, 1],
                                       absorption['eps_mean_wl'])
            np.testing.assert_allclose(ecd_data[:, 1],
                                       ecd['delta_eps_mean_wl'])

            assert absorption_file.read_text().splitlines()[0] == (
                'Wavelength [nm],epsilon [L mol^-1 cm^-1]')
            assert ecd_file.read_text().splitlines()[0] == (
                'Wavelength [nm],Delta epsilon [L mol^-1 cm^-1]')

    @pytest.mark.parametrize(
        ('property_name', 'expected_header'),
        [
            ('absorption(cpp)',
             'Wavelength [nm],epsilon [L mol^-1 cm^-1]'),
            ('ecd(cpp)',
             'Wavelength [nm],delta epsilon [L mol^-1 cm^-1]'),
        ],
    )
    def test_save_cpp_spectra(self, tmp_path, property_name,
                              expected_header):

        driver = _make_driver()
        driver.rsp_dict = {'property': property_name}

        csv_file = tmp_path / f'{property_name}.csv'
        spectrum = driver.save_uv_vis_spectra(
            {'rsp_all': _make_cpp_results()},
            averaged_spectra_filename=csv_file,
        )

        if driver.rank == mpi_master():
            assert csv_file.is_file()
            assert csv_file.read_text().splitlines()[0] == expected_header
            np.testing.assert_allclose(spectrum['y_mean'],
                                       spectrum['y_all'][0])

    def test_save_spectra_requires_property(self, tmp_path):

        driver = _make_driver()

        with pytest.raises(ValueError, match='property'):
            driver.save_uv_vis_spectra(
                {'rsp_all': _make_tddft_results()},
                averaged_spectra_filename=tmp_path / 'spectrum.csv',
            )
