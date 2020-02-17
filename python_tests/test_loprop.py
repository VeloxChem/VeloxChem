import io
import sys
import textwrap
from unittest.mock import MagicMock, patch

import pytest

from veloxchem.main import main
from veloxchem.inputparser import InputParser
from veloxchem.loprop import LoPropDriver


@pytest.fixture
def sample():
    inp = textwrap.dedent(
        """
        @jobs
        task: loprop
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
        """
    )
    return inp


@patch('veloxchem.main.MpiTask')
def test_loprop_called_from_main(mock_mpi, sample, tmpdir):
    # out = io.StringIO()

    task = mock_mpi()
    task.input_dict = {'jobs': {'task': 'loprop'}}
    input_file = f'{tmpdir/"water.inp"}'
    with open(input_file, 'w') as f:
        f.write(sample)

    with patch('veloxchem.main.LoPropDriver', autospec=True) as mock_prop:
        sys.argv[1:] = [input_file]
        task = mock_mpi()
        main()

    mock_prop.assert_called_with(task)
    mock_prop(task).compute.assert_called()


@patch('veloxchem.loprop.h5py')
@patch('veloxchem.loprop.MpiTask')
@patch('veloxchem.loprop.OverlapIntegralsDriver')
def test_overlap_called(mock_ovldrv, mock_mpi, mock_h5py):

    task = mock_mpi()
    task.input_dict = {'loprop': {'checkpoint_file': 'water.loprop.h5'}}
    lpd = LoPropDriver(task)
    lpd.save_overlap()

    mock_ovldrv.assert_called_with(task.mpi_comm)
    mock_ovldrv().compute.assert_called_with(task.molecule, task.ao_basis)

    mock_h5py.File.assert_called_with('water.loprop.h5', 'w')


def test_input_dict(sample, tmpdir):
    input_file = f'{tmpdir/"water.inp"}'
    with open(input_file, 'w') as f:
        f.write(sample)

    ip = InputParser(input_file)

    assert ip.input_dict['jobs']['task'] == 'loprop'
    assert 'loprop' in ip.input_dict
