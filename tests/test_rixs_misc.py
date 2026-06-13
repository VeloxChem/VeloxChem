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

        driver.update_settings({'photon_energy': '1.25, 2.50, 3.75'})

        assert driver.photon_energy == pytest.approx([1.25, 2.50, 3.75])

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
