from mpi4py import MPI
from pathlib import Path
import h5py
import numpy as np
import pytest

from veloxchem.rixsdriver import RixsDriver


class RecordingOutput:

    def __init__(self):
        self.infos = []
        self.warnings = []
        self.flushed = False

    def print_info(self, text):
        self.infos.append(text)

    def print_warning(self, text):
        self.warnings.append(text)

    def flush(self):
        self.flushed = True


@pytest.mark.solvers
class TestRixsDriverUtilities:

    def test_update_settings_parses_photon_energy_string(self):

        driver = RixsDriver()

        driver.update_settings({'photon_energy': '1.25, 2.50 3.75'})

        assert driver.photon_energy == pytest.approx([1.25, 2.50, 3.75])

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    @pytest.mark.parametrize('rsp_dict', [
        {'photon_energy': ''},
        {'photon_energy': '   '},
        {'photon_energy': []},
        {'photon_energy': ()},
    ])
    def test_update_settings_rejects_empty_photon_energy_input(self, rsp_dict):

        driver = RixsDriver()

        with pytest.raises(ValueError, match='photon_energy'):
            driver.update_settings(rsp_dict)

    def test_print_info_warns_when_no_results_available(self):

        ostream = RecordingOutput()
        driver = RixsDriver(MPI.COMM_WORLD, ostream)

        driver.print_info()

        assert ostream.warnings == ['Compute first!']
        assert ostream.infos == []

    def test_print_info_reports_state_summary_via_output_stream(self, capsys):

        ostream = RecordingOutput()
        driver = RixsDriver(MPI.COMM_WORLD, ostream)
        driver._orb_and_state_dict = {
            'num_intermediate_states': 2,
            'num_final_states': 3,
            'mo_core_indices': [0],
            'mo_valence_indices': [4, 5],
            'mo_virtual_indices': [6, 7],
            'core_states': [0, 2],
            'valence_states': [1, 3, 4],
        }

        driver.print_info()

        captured = capsys.readouterr()
        assert captured.out == ''
        assert ostream.warnings == []
        assert ostream.infos == [
            'Number of states: (intermediate, final): (2, 3)',
            'MO indices (core, valence, virtual): ([0], [4, 5], [6, 7])',
            'State indices (intermediate, final): ([0, 2], [1, 3, 4])',
        ]

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip file I/O for multiple MPI processes')
    def test_write_hdf5_accepts_path_objects_and_appends_extension(self,
                                                                   tmp_path):

        driver = RixsDriver()
        driver.photon_energy = [1.0]
        driver.cross_sections = np.array([[2.0]])
        driver.emission_enes = np.array([[3.0]])
        driver.ene_losses = np.array([[4.0]])
        driver.elastic_cross_sections = np.array([5.0])
        driver.scattering_amplitudes = np.array([[[[6.0 + 0.0j]]]])

        target = tmp_path / 'rixs_results'
        driver._write_hdf5(target)

        assert target.with_suffix('.h5').is_file()

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip file I/O for multiple MPI processes')
    def test_write_hdf5_overwrites_existing_rixs_datasets(self, tmp_path):

        driver = RixsDriver()
        driver.photon_energy = [1.0]
        driver.cross_sections = np.array([[2.0]])
        driver.emission_enes = np.array([[3.0]])
        driver.ene_losses = np.array([[4.0]])
        driver.elastic_cross_sections = np.array([5.0])
        driver.scattering_amplitudes = np.array([[[[6.0 + 0.0j]]]])

        target = tmp_path / 'shared_results.h5'
        driver._write_hdf5(target)

        driver.photon_energy = [7.0]
        driver.cross_sections = np.array([[8.0]])
        driver.emission_enes = np.array([[9.0]])
        driver.ene_losses = np.array([[10.0]])
        driver.elastic_cross_sections = np.array([11.0])
        driver.scattering_amplitudes = np.array([[[[12.0 + 0.0j]]]])
        driver._write_hdf5(target)

        with h5py.File(target, 'r') as h5f:
            assert h5f['rixs/photon_energies'][()] == pytest.approx([7.0])
            assert np.allclose(h5f['rixs/cross_sections'][()], [[8.0]])
            assert np.allclose(h5f['rixs/emission_energies'][()], [[9.0]])
            assert np.allclose(h5f['rixs/energy_losses'][()], [[10.0]])
            assert h5f['rixs/elastic_cross_sections'][()] == pytest.approx([11.0])
            assert np.allclose(h5f['rixs/scattering_amplitudes'][()],
                               [[[[12.0 + 0.0j]]]])
