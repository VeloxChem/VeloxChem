from unittest.mock import patch
import numpy.testing as npt
import numpy as np
import textwrap
import pytest

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.veloxchemlib import is_single_node
from veloxchem.inputparser import InputParser


@patch.object(InputParser, 'parse')
def test_input_dict(mock_parse):

    ip = InputParser('foo.inp')

    assert mock_parse.called
    assert ip.inpname == 'foo.inp'
    assert ip.outname is None
    assert ip.input_dict == {}
    assert ip.get_dict() is ip.input_dict


def test_full_input(tmpdir):

    if not is_mpi_master():
        return

    with open(tmpdir / 'h2o.inp', 'w') as f:
        f.write(
            textwrap.dedent("""
                @jobs
                task: hf
                @end

                @method settings
                basis: cc-pvdz
                @end

                @molecule
                charge: 0
                multiplicity: 1
                units: au
                xyz:
                O   0.0   0.0   0.0
                H   0.0   1.4   1.1
                H   0.0  -1.4   1.1
                @end
                """))

    ip = InputParser(str(tmpdir / 'h2o.inp'))

    expected = {
        'jobs': {
            'task': 'hf',
        },
        'method_settings': {
            'basis': 'cc-pvdz',
        },
        'molecule': {
            'charge': '0',
            'multiplicity': '1',
            'units': 'au',
            'xyz': ['O 0.0 0.0 0.0', 'H 0.0 1.4 1.1', 'H 0.0 -1.4 1.1'],
        },
        'filename': str(tmpdir / 'h2o'),
        'scf': {
            'checkpoint_file': str(tmpdir / 'h2o.scf.h5'),
        },
        'response': {
            'checkpoint_file': str(tmpdir / 'h2o.rsp.h5'),
        },
        'exciton': {
            'checkpoint_file': str(tmpdir / 'h2o.exciton.h5'),
        },
        'loprop': {
            'checkpoint_file': str(tmpdir / 'h2o.loprop.h5'),
        },
    }

    assert ip.input_dict == expected


def test_error_in_input(tmpdir):

    if is_single_node():

        with open(tmpdir / 'h2o.inp', 'w') as f:
            f.write(
                textwrap.dedent("""
                    @jobs
                    task: hf

                    @method settings
                    basis: cc-pvdz
                    @end

                    @molecule
                    charge: 0
                    multiplicity: 1
                    units: au
                    xyz:
                    O   0.0   0.0   0.0
                    H   0.0   1.4   1.1
                    H   0.0  -1.4   1.1
                    @end
                    """))

        with pytest.raises(AssertionError) as info:
            InputParser(str(tmpdir / 'h2o.inp'))

        assert 'bad syntax' in str(info.value)


@pytest.mark.parametrize('input_frequencies, expected', [
    ([], []),
    ([1.0, 2.0], [1.0, 2.0]),
    (np.array([1.0, 2.0]), [1.0, 2.0]),
    (
        "0.0 - 0.1 (0.05), 0.5 - 1.0 (0.1), 2.0",
        [0.0, 0.05, 0.1, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0],
    ),
    (
        "0.0-0.1-0.05, 0.5-1.0-0.1, 2.0",
        [0.0, 0.05, 0.1, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0],
    ),
])
def test_parse_frequencies(input_frequencies, expected):

    if is_mpi_master():
        npt.assert_allclose(InputParser.parse_frequencies(input_frequencies),
                            expected)
