from pathlib import Path
from unittest.mock import patch

import numpy as np
import numpy.testing as npt
import pytest
from veloxchem.inputparser import (InputParser, parse_seq_range,
                                   parse_seq_fixed, parse_bool, parse_str)
from veloxchem.veloxchemlib import is_mpi_master, is_single_node


@patch.object(InputParser, 'parse')
def test_input_dict(mock_parse):

    ip = InputParser('foo.inp')

    assert mock_parse.called
    assert ip.inpname == 'foo.inp'
    assert ip.outname is None
    assert ip.input_dict == {}
    assert ip.get_dict() is ip.input_dict


def test_full_input(tmpdir):

    here = Path(__file__).parent
    inpfile = here / 'inputs' / 'water.inp'

    ip = InputParser(str(inpfile))

    expected = {
        'jobs': {
            'task': 'scf',
        },
        'method_settings': {
            'basis': 'aug-cc-pvdz',
        },
        'molecule': {
            'charge': '0',
            'multiplicity': '1',
            'units': 'au',
            'xyz': ['O 0.0 0.0 0.0', 'H 0.0 1.4 1.1', 'H 0.0 -1.4 1.1'],
        },
        'filename': str(inpfile.with_name(inpfile.stem)),
    }

    assert ip.input_dict == expected


def test_error_in_input(tmpdir):

    if is_single_node():
        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'water_error.inp'

        with pytest.raises(AssertionError) as info:
            InputParser(str(inpfile))
            assert 'bad syntax' in str(info.value)


@pytest.mark.parametrize(
    'input_seq, expected',
    [
        ([], []),
        ([1.0, 2.0], [1.0, 2.0]),
        (np.array([1.0, 2.0]), [1.0, 2.0]),
        (
            '0.0 - 0.1 (0.05), 0.5 - 1.0 (0.1), 2.0',
            [0.0, 0.05, 0.1, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0],
        ),
        (
            '0.0-0.1-0.05, 0.5-1.0-0.1, 2.0',
            [0.0, 0.05, 0.1, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0],
        ),
    ],
)
def test_parse_seq_range(input_seq, expected):

    if is_mpi_master():
        npt.assert_allclose(parse_seq_range(input_seq), expected)


@pytest.mark.parametrize(
    'input_seq, flag, expected',
    [
        ([], 'float', []),
        ([1.0, 2.0], 'float', [1.0, 2.0]),
        (np.array([1.0, 2.0, 3.0]), 'float', [1.0, 2.0, 3.0]),
        ([], 'int', []),
        ([1, 2], 'int', [1, 2]),
        (np.array([1, 2, 3]), 'int', [1, 2, 3]),
    ],
)
def test_parse_seq_fixed(input_seq, flag, expected):

    if is_mpi_master():
        npt.assert_allclose(parse_seq_fixed(input_seq, flag), expected)


@pytest.mark.parametrize(
    'input_bool, expected',
    [
        ('yes', True),
        ('Y', True),
        ('no', False),
        ('N', False),
    ],
)
def test_parse_bool(input_bool, expected):

    if is_mpi_master():
        assert parse_bool(input_bool) == expected


@pytest.mark.parametrize(
    'input_str, flag, expected',
    [
        ('Qq', 'upper', 'QQ'),
        ('Qq', 'lower', 'qq'),
        ('Qq', None, 'Qq'),
    ],
)
def test_parse_str(input_str, flag, expected):

    if is_mpi_master():
        assert parse_str(input_str, flag) == expected
