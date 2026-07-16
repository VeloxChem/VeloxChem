from mpi4py import MPI
import numpy as np
import pytest

from veloxchem.veloxchemlib import (fine_structure_constant,
                                    extinction_coefficient_from_beta,
                                    hartree_in_ev, hartree_in_inverse_nm)
from veloxchem.tdacppsolver import ComplexResponseTdaSolver
from veloxchem.errorhandler import VeloxChemError


class RecordingOutput:

    def __init__(self):
        self.headers = []
        self.infos = []
        self.blank_count = 0
        self.flushed = False

    def print_header(self, text):
        self.headers.append(text)

    def print_info(self, text):
        self.infos.append(text)

    def print_blank(self):
        self.blank_count += 1

    def flush(self):
        self.flushed = True


def _build_cpp_results(freqs):

    response_functions = {}
    for w in freqs:
        for i, a in enumerate('xyz', start=1):
            for j, b in enumerate('xyz', start=1):
                response_functions[(a, b, w)] = complex(i + j, -(i * j * w))

    return {
        'frequencies': list(freqs),
        'response_functions': response_functions,
    }


@pytest.mark.solvers
class TestComplexResponseTdaSolverMisc:

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_set_property_and_validate_spectrum_inputs(self):

        driver = ComplexResponseTdaSolver()

        driver.set_cpp_property('Absorption')
        assert driver.property == 'absorption'
        assert driver.a_operator == 'electric dipole'
        assert driver.b_operator == 'electric dipole'
        assert driver.a_components == 'xyz'
        assert driver.b_components == 'xyz'

        with pytest.raises(VeloxChemError,
                           match='ComplexResponseTdaSolver: invalid CPP property'):
            driver.set_cpp_property('raman')

        with pytest.raises(
                VeloxChemError,
                match='ComplexResponseTdaSolver.get_spectrum: x_unit should be au, ev or nm'
        ):
            driver.get_spectrum(_build_cpp_results([0.4]), 'cm-1')

    def test_absorption_and_ecd_spectra(self):

        rsp_results = _build_cpp_results([0.0, 0.4])
        driver = ComplexResponseTdaSolver()

        driver.set_cpp_property('absorption')
        spectrum_au = driver.get_spectrum(rsp_results, 'au')
        spectrum_ev = driver.get_spectrum(rsp_results, 'ev')
        spectrum_nm = driver.get_spectrum(rsp_results, 'nm')

        alpha_bar = (0.4 + 1.6 + 3.6) / 3.0
        sigma_ref = 4.0 * np.pi * 0.4 * alpha_bar * fine_structure_constant()

        assert spectrum_au['x_label'] == 'Photon energy [a.u.]'
        assert spectrum_au['y_label'] == 'Absorption cross-section [a.u.]'
        assert spectrum_au['x_data'] == pytest.approx([0.4])
        assert spectrum_au['y_data'] == pytest.approx([sigma_ref])

        assert spectrum_ev['x_label'] == 'Photon energy [eV]'
        assert spectrum_ev['x_data'] == pytest.approx([0.4 * hartree_in_ev()])

        assert spectrum_nm['x_label'] == 'Wavelength [nm]'
        assert spectrum_nm['x_data'] == pytest.approx(
            [(1.0 / hartree_in_inverse_nm()) / 0.4])

        driver.set_cpp_property('ecd')
        ecd_spectrum = driver.get_spectrum(rsp_results, 'au')
        beta_ref = -((1.0 + 4.0 + 9.0) / (3.0 * 0.4))
        delta_ref = beta_ref * (0.4**2) * extinction_coefficient_from_beta()

        assert ecd_spectrum['x_label'] == 'Photon energy [a.u.]'
        assert ecd_spectrum[
            'y_label'] == 'Molar circular dichroism [L mol$^{-1}$ cm$^{-1}$]'
        assert ecd_spectrum['x_data'] == pytest.approx([0.4])
        assert ecd_spectrum['y_data'] == pytest.approx([delta_ref])

    def test_print_helpers_cover_linear_and_zero_frequency_outputs(self):

        rsp_results = _build_cpp_results([0.0, 0.4])
        ostream = RecordingOutput()
        driver = ComplexResponseTdaSolver(ostream=ostream)
        driver.a_components = 'xy'
        driver.b_components = 'yz'

        driver._cur_iter = 1
        driver.nonlinear = False
        driver.print_level = 2
        driver._print_iteration({('y', 0.4): 2.0e-3}, [('y', 0.4, 1.2 - 0.3j)])

        assert any('Residuals (Max,Min)' in line for line in ostream.headers)
        assert any('Operator:' in line for line in ostream.headers)
        assert any('<<y;y>>_0.4000' in line for line in ostream.headers)
        assert ostream.flushed is True

        driver.property = 'absorption'
        driver._print_results(rsp_results, ostream)
        assert any('Response Functions at Given Frequencies' in line
                   for line in ostream.headers)
        assert any('Linear Absorption Cross-Section' in line
                   for line in ostream.headers)

        zero_freq_results = _build_cpp_results([0.0])

        driver.property = 'absorption'
        driver._print_absorption_results(zero_freq_results, ostream)
        assert any('No linear absorption spectrum at zero frequency.' in line
                   for line in ostream.headers)

        driver.property = 'ecd'
        driver._print_ecd_results(zero_freq_results, ostream)
        assert any('No circular dichroism spectrum at zero frequency.' in line
                   for line in ostream.headers)

    def test_plot_covers_none_absorption_and_ecd_paths(self, monkeypatch,
                                                       capsys):

        plt = pytest.importorskip('matplotlib.pyplot')
        pytest.importorskip('scipy.interpolate')
        monkeypatch.setattr(plt, 'show', lambda: None)

        freqs = [0.38, 0.39, 0.40, 0.41]
        rsp_results = _build_cpp_results(freqs)
        driver = ComplexResponseTdaSolver()

        driver.plot({}, x_unit='au')
        assert 'Nothing to plot for this complex response calculation.' in (
            capsys.readouterr().out)

        driver.set_cpp_property('absorption')
        driver.plot(rsp_results, x_unit='nm', plot_scatter=True)

        driver.set_cpp_property('ecd')
        driver.plot(rsp_results, x_unit='ev', plot_scatter=False)

        plt.close('all')

    def test_checkpoint_helpers_handle_empty_and_nonlinear_restart_states(
            self, tmp_path, monkeypatch):

        class DummyDistArray:

            def __init__(self, name):
                self.name = name
                self.calls = []

            def append_to_hdf5_file(self, checkpoint_file, label):
                self.calls.append((checkpoint_file, label))

        driver = ComplexResponseTdaSolver(ostream=RecordingOutput())
        driver._write_checkpoint(None, None, {}, {}, ['unused'])

        driver.filename = str(tmp_path / 'tdacppsolver_checkpoint')
        driver.nonlinear = True
        driver._dist_bger = DummyDistArray('bger')
        driver._dist_e2bger = DummyDistArray('e2bger')
        driver._dist_fock_ger = DummyDistArray('fock')

        written_settings = []

        monkeypatch.setattr('veloxchem.tdacppsolver.write_rsp_hdf5',
                            lambda *args: True)
        monkeypatch.setattr(driver, '_write_settings_to_checkpoint',
                            lambda filename: written_settings.append(filename))

        driver._write_checkpoint(None, None, {}, {},
                                 ['CLR_bger', 'CLR_e2bger', 'CLR_Fock_ger'])

        expected_checkpoint = f'{driver.filename}.rsp.h5'
        assert driver.checkpoint_file == expected_checkpoint
        assert written_settings == [expected_checkpoint]
        assert driver._dist_bger.calls == [(expected_checkpoint, 'CLR_bger')]
        assert driver._dist_e2bger.calls == [(expected_checkpoint,
                                             'CLR_e2bger')]
        assert driver._dist_fock_ger.calls == [(expected_checkpoint,
                                               'CLR_Fock_ger')]
        assert any('Time spent in writing checkpoint file:' in info
                   for info in driver.ostream.infos)

        restored_labels = []

        def fake_read_from_hdf5_file(checkpoint_file, label, comm):
            restored_labels.append((checkpoint_file, label))
            return f'restored:{label}'

        monkeypatch.setattr(
            'veloxchem.tdacppsolver.DistributedArray.read_from_hdf5_file',
            staticmethod(fake_read_from_hdf5_file))

        driver._read_checkpoint(['CLR_bger', 'CLR_e2bger', 'CLR_Fock_ger'])

        assert driver._dist_bger == 'restored:CLR_bger'
        assert driver._dist_e2bger == 'restored:CLR_e2bger'
        assert driver._dist_fock_ger == 'restored:CLR_Fock_ger'
        assert restored_labels == [
            (expected_checkpoint, 'CLR_bger'),
            (expected_checkpoint, 'CLR_e2bger'),
            (expected_checkpoint, 'CLR_Fock_ger'),
        ]
        assert any('Restarting from checkpoint file:' in info
                   for info in driver.ostream.infos)

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_property_density_validation_errors(self):

        driver = ComplexResponseTdaSolver()

        with pytest.raises(VeloxChemError,
                           match='get_cpp_property_densities: Invalid CPP property'):
            driver.get_cpp_property_densities(None, None, None,
                                              {'solutions': {}}, 0.4)

        driver.property = 'absorption'

        with pytest.raises(
                VeloxChemError,
                match='get_cpp_property_densities: Could not find frequency 0.4 in CPP results'
        ):
            driver.get_cpp_property_densities(None, None, None,
                                              {'solutions': {}}, 0.4)
