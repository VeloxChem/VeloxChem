from mpi4py import MPI
from pathlib import Path
from unittest.mock import patch
import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.inputparser import get_random_string_parallel
from veloxchem.errorhandler import VeloxChemError
from veloxchem.main import main, check_checkpoint_files


class TestInputCheckpointFile:
    """
    A checkpoint file set in the input must not be the final h5 file of the
    run, f'{filename}.h5', where the results are written.
    """

    def test_other_names_are_accepted(self):

        check_checkpoint_files({
            'filename': 'data/water',
            'jobs': {
                'task': 'response'
            },
            'scf': {
                'checkpoint_file': 'data/water_chk.h5'
            },
            'response': {
                'checkpoint_file': 'data/water_rsp.h5'
            },
        })

        check_checkpoint_files({
            'filename': 'data/water',
            'scf': {},
        })

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='pytest.raises only valid in serial')
    @pytest.mark.parametrize('group, checkpoint_file', [
        ('scf', 'data/water.h5'),
        ('response', './data/water.h5'),
        ('scf', 'data/../data/water.h5'),
    ])
    def test_the_final_h5_file_is_refused(self, group, checkpoint_file):

        with pytest.raises(VeloxChemError,
                           match=f'@{group}: checkpoint_file .* final h5'):
            check_checkpoint_files({
                'filename': 'data/water',
                group: {
                    'checkpoint_file': checkpoint_file
                },
            })

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='an error in parallel aborts the run')
    def test_input_stops_before_the_scf(self, capsys):

        here = Path(__file__).parent
        random_string = get_random_string_parallel()
        input_file = here / 'data' / f'vlx_chk_{random_string}.inp'
        final_h5 = input_file.with_suffix('.h5')

        input_lines = f"""
            @jobs
            task: scf
            @end

            @method settings
            basis: def2-svp
            @end

            @scf
            checkpoint_file: {final_h5}
            @end

            @molecule
            charge: 0
            multiplicity: 1
            xyz:
            O   0.0   0.0   0.0
            H   0.0   0.757   0.587
            H   0.0  -0.757   0.587
            @end
        """

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            with open(input_file, 'w') as f_inp:
                for line in input_lines.splitlines():
                    if line.strip():
                        f_inp.write(f'{line.strip()}\n')

        try:
            with patch.object(sys, 'argv', ['vlx', str(input_file)]):
                assert main() == 1

            captured = capsys.readouterr()

            assert '@scf: checkpoint_file' in captured.err
            assert 'final h5 file' in captured.err
            assert 'SCF converged' not in captured.out
            assert not final_h5.is_file()

        finally:
            for path in [input_file, final_h5]:
                if path.is_file():
                    path.unlink()
