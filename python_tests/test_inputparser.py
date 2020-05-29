from mpi4py import MPI
from unittest.mock import patch, mock_open
import numpy.testing as npt
import numpy as np
import textwrap
import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.inputparser import InputParser
from veloxchem.inputparser import parse_frequencies


@patch.object(InputParser, 'parse')
def test_input_dict(mock_parse):

    ip = InputParser('foo.inp')

    assert mock_parse.called
    assert ip.filename == 'foo.inp'
    assert ip.outname is None
    assert ip.input_dict == {}
    assert ip.get_dict() is ip.input_dict


@pytest.mark.skipif(sys.version_info < (3, 7),
                    reason="mock_open does not support iteration in 3.6")
@patch.object(InputParser, 'parse')
def test_read_file(mock_parse):
    with patch(
            'veloxchem.inputparser.open',
            mock_open(read_data=textwrap.dedent("""
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
                """))) as mocked_open:

        ip = InputParser('foo.inp')
        ip.read_file()

    mocked_open.assert_called_once_with('foo.inp', 'r')
    assert ip.content == textwrap.dedent("""\
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
        O 0.0 0.0 0.0
        H 0.0 1.4 1.1
        H 0.0 -1.4 1.1
        @end
        """)


@pytest.mark.parametrize('content', [
    '@end@foo',
    '@foo@bar',
    '@foo@end@bar',
])
@patch.object(InputParser, 'parse')
def test_incomplete_group(mock_parse, content):

    ip = InputParser('foo.inp')
    ip.content = content

    with pytest.raises(SyntaxError):
        ip.incomplete_group_check()


@pytest.mark.parametrize('content', [
    '@foo\n@end',
])
@patch.object(InputParser, 'parse')
def test_empty_group_check(mock_parse, content):

    ip = InputParser('foo.inp')
    ip.content = content

    with pytest.raises(SyntaxError):
        ip.empty_group_check()


@pytest.mark.parametrize('content, cleared', [
    ('@foo@end...@bar@end', '@foo@end\n@bar@end'),
    ('@foo@end\n...@bar@end\n', '@foo@end\n@bar@end\n'),
])
@patch.object(InputParser, 'parse')
def test_interspace(mock_parse, content, cleared):

    ip = InputParser('foo.inp')
    ip.content = content
    ip.clear_interspace()

    assert ip.content == cleared


@pytest.mark.parametrize(
    'content, grouped',
    [
        ('', []),
        (
            '@foo\none foo-line\n@end',
            [['foo', 'one foo-line']],
        ),
        (
            '@foo\none foo-line\n@end\n@bar\none bar-line\n@end',
            [['foo', 'one foo-line'], ['bar', 'one bar-line']],
        ),
    ],
    ids=['empty', 'one-group', 'two-groups'],
)
@patch.object(InputParser, 'parse')
def test_grouped(mock_parse, content, grouped):

    ip = InputParser('foo.inp')
    ip.content = content
    ip.groupsplit()

    assert ip.grouplist == grouped


@pytest.mark.parametrize(
    'grouped, todict',
    [
        ([], {}),
        (
            [['foo', 'one: foo-line']],
            {
                'foo': {
                    'one': 'foo-line'
                }
            },
        ),
        (
            [['foo', 'one: foo-line'], ['bar', 'two: bar-line']],
            {
                'foo': {
                    'one': 'foo-line'
                },
                'bar': {
                    'two': 'bar-line'
                }
            },
        ),
        (
            [['molecule', 'charge: 0', 'units: au']],
            {
                'molecule': {
                    'charge': '0',
                    'units': 'au',
                    'xyzstr': ''
                }
            },
        ),
        (
            [['molecule', 'charge: 0', 'units: au', 'xyz', 'H 0 0 0']],
            {
                'molecule': {
                    'charge': '0',
                    'units': 'au',
                    'xyzstr': 'xyz\nH 0 0 0'
                }
            },
        ),
    ],
    # ids=['empty', 'one-group', 'two-groups'],
)
@patch.object(InputParser, 'parse')
def test_convert_dict(mock_parse, grouped, todict):
    if MPI.COMM_WORLD.Get_rank() == mpi_master():

        ip = InputParser('foo.inp')
        expected = {
            'input_file': 'foo.inp',
            'scf': {
                'checkpoint_file': 'foo.scf.h5'
            },
            'response': {
                'checkpoint_file': 'foo.rsp.h5'
            },
            'exciton': {
                'checkpoint_file': 'foo.exciton.h5'
            },
            'loprop': {
                'checkpoint_file': 'foo.loprop.h5'
            },
            **todict,
        }

        ip.grouplist = grouped
        ip.convert_dict()

        assert ip.input_dict == expected


@pytest.mark.parametrize(
    'grouped, todict',
    [
        ([], {}),
        (
            [['foo', 'one: foo-line']],
            {
                'foo': {
                    'one': 'foo-line'
                }
            },
        ),
        (
            [['foo', 'one: foo-line'], ['bar', 'two: bar-line']],
            {
                'foo': {
                    'one': 'foo-line'
                },
                'bar': {
                    'two': 'bar-line'
                }
            },
        ),
    ],
    ids=['empty', 'one-group', 'two-groups'],
)
@patch.object(InputParser, 'parse')
def test_convert_dict_with_output_file(mock_parse, grouped, todict):

    ip = InputParser('bar.inp', 'foo.out')
    expected = {
        'input_file': 'bar.inp',
        'scf': {
            'checkpoint_file': 'foo.scf.h5'
        },
        'response': {
            'checkpoint_file': 'foo.rsp.h5'
        },
        'exciton': {
            'checkpoint_file': 'foo.exciton.h5'
        },
        'loprop': {
            'checkpoint_file': 'foo.loprop.h5'
        },
        **todict,
    }

    ip.grouplist = grouped
    ip.convert_dict()

    assert ip.input_dict == expected


def test_full_input(tmpdir):

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
            'input_file': str(tmpdir / 'h2o.inp'),
            'scf': {
                'checkpoint_file': str(tmpdir / 'h2o.scf.h5')
            },
            'response': {
                'checkpoint_file': str(tmpdir / 'h2o.rsp.h5')
            },
            'exciton': {
                'checkpoint_file': str(tmpdir / 'h2o.exciton.h5')
            },
            'loprop': {
                'checkpoint_file': str(tmpdir / 'h2o.loprop.h5')
            },
        }

        assert ip.input_dict == expected


def test_missing_file():
    if MPI.COMM_WORLD.Get_size() == 1:

        with pytest.raises(AssertionError) as nofileinfo:
            InputParser('no_file')

        assert str(nofileinfo.value) == 'input parser: cannot open file no_file'


def test_error_in_input(tmpdir):
    if MPI.COMM_WORLD.Get_size() == 1:

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


@pytest.mark.parametrize('input_frequencies, expected',
                         [([], []), ([1.0, 2.0], [1.0, 2.0]),
                          (np.array([1.0, 2.0]), [1.0, 2.0]),
                          (
                              "0.0 - 0.1 (0.05), 0.5 - 1.0 (0.1), 2.0",
                              [0.0, 0.05, 0.5, 0.6, 0.7, 0.8, 0.9, 2.0],
                          ),
                          (
                              "0.0-0.1-0.05, 0.5-1.0-0.1, 2.0",
                              [0.0, 0.05, 0.5, 0.6, 0.7, 0.8, 0.9, 2.0],
                          )])
def test_parse_frequencies(input_frequencies, expected):
    if MPI.COMM_WORLD.Get_rank() == mpi_master():
        npt.assert_allclose(parse_frequencies(input_frequencies), expected)
