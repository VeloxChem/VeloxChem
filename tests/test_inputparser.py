from mpi4py import MPI
from pathlib import Path
from unittest.mock import patch
import numpy as np
import numpy.testing as npt
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.inputparser import (InputParser, parse_seq_range,
                                   parse_seq_fixed, parse_bool, parse_str,
                                   parse_input, unparse_input,
                                   write_unparsed_input_to_hdf5,
                                   read_unparsed_input_from_hdf5)


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
    inpfile = here / 'data' / 'water.inp'

    ip = InputParser(str(inpfile))

    expected = {
        'keyline': '',
        'jobs': {
            'task': 'scf',
        },
        'method_settings': {
            'basis': 'aug-cc-pvdz',
        },
        'scf': {},
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

    if MPI.COMM_WORLD.Get_size() == 1:
        here = Path(__file__).parent
        inpfile = here / 'data' / 'water_error.inp'

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

    if MPI.COMM_WORLD.Get_rank() == mpi_master():
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

    if MPI.COMM_WORLD.Get_rank() == mpi_master():
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

    if MPI.COMM_WORLD.Get_rank() == mpi_master():
        assert parse_bool(input_bool) == expected


@pytest.mark.parametrize(
    'input_str, flag, expected',
    [
        ('Diis', 'upper', 'DIIS'),
        ('Diis', 'lower', 'diis'),
        ('Diis', None, 'Diis'),
    ],
)
def test_parse_str(input_str, flag, expected):

    if MPI.COMM_WORLD.Get_rank() == mpi_master():
        assert parse_str(input_str, flag) == expected


def test_unparse_input_roundtrip():

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    class Dummy:
        pass

    keyword_types = {
        'name': 'str',
        'acc_type': 'str_upper',
        'coordsys': 'str_lower',
        'max_iter': 'int',
        'eri_thresh': 'float',
        'restart': 'bool',
        'constraints': 'list',
        'atom_types': 'seq_fixed_str',
        'cube_points': 'seq_fixed_int',
        'cube_origin': 'seq_fixed',
        'frequencies': 'seq_range',
    }

    input_dictionary = {
        'name': 'checkpoint.h5',
        'acc_type': 'Diis',
        'coordsys': 'TrIc',
        'max_iter': '300',
        'eri_thresh': '1.0e-12',
        'restart': 'Y',
        'constraints': ['bond 1 2', 'angle 1 2 3'],
        'atom_types': 'c3,c3,hc',
        'cube_points': '80, 80, 80',
        'cube_origin': '0.0, 0.1, 0.2',
        'frequencies': '0.0 - 0.1 (0.05), 0.5 - 1.0 (0.1), 2.0',
    }

    parsed = Dummy()
    parse_input(parsed, keyword_types, input_dictionary)

    roundtrip_dictionary = unparse_input(parsed, keyword_types)

    reparsed = Dummy()
    parse_input(reparsed, keyword_types, roundtrip_dictionary)

    for key in keyword_types:
        assert getattr(reparsed, key) == getattr(parsed, key)


def test_unparsed_input_hdf5_roundtrip(tmpdir):

    if MPI.COMM_WORLD.Get_rank() != mpi_master():
        return

    input_dictionary = {
        'label': 'checkpoint.h5',
        'constraints': ['bond 1 2', 'angle 1 2 3'],
        'optional': None,
    }

    h5file = str(Path(tmpdir) / 'input.h5')

    write_unparsed_input_to_hdf5(h5file, input_dictionary, group_name='test')
    recovered = read_unparsed_input_from_hdf5(h5file, group_name='test')

    assert recovered == input_dictionary
