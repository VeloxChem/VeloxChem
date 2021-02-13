from mpi4py import MPI
from unittest.mock import patch
from unittest.mock import mock_open
import numpy.testing as npt
import numpy as np
import textwrap
import pytest
import sys

from veloxchem.veloxchemlib import mpi_master
from veloxchem.inputparser import InputParser


@patch.object(InputParser, 'parse')
def test_input_dict(mock_parse):

    ip = InputParser('foo.inp')

    assert mock_parse.called
    assert ip.inpname == 'foo.inp'
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
            'filename': str(tmpdir / 'h2o'),
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
    if MPI.COMM_WORLD.Get_rank() == mpi_master():
        npt.assert_allclose(InputParser.parse_frequencies(input_frequencies),
                            expected)
