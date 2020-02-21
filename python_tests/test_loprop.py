import sys
import textwrap
from unittest.mock import patch, MagicMock

import pytest
import numpy.testing as npt

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


@patch('veloxchem.main.ScfUnrestrictedDriver')
@patch('veloxchem.main.ScfRestrictedDriver')
@patch('veloxchem.main.MpiTask')
def test_loprop_called_from_main(mock_mpi, mock_scf, mock_uscf, sample, tmpdir):
    """
    Verify that LoPropDriver is called from main, loprop task in input
    """

    # given
    task = mock_mpi()
    task.input_dict = {'jobs': {'task': 'loprop'}, 'loprop': {}}
    input_file = f'{tmpdir/"water.inp"}'
    with open(input_file, 'w') as f:
        f.write(sample)
    sys.argv[1:] = [input_file]

    # when
    with patch('veloxchem.main.LoPropDriver', autospec=True) as mock_prop:
        main()

    # then
    mock_prop.assert_called_with(task)
    mock_prop(task).compute.assert_called()


@patch('veloxchem.main.ScfUnrestrictedDriver')
@patch('veloxchem.main.ScfRestrictedDriver')
@patch('veloxchem.main.MpiTask')
def test_scf_called_from_main(mock_mpi, mock_scf, mock_uscf, sample, tmpdir):
    """
    Verify that SCF is called from main, loprop task in input
    """

    # given
    task = mock_mpi()
    task.input_dict = {'jobs': {'task': 'loprop'}, 'loprop': {}}
    input_file = f'{tmpdir/"water.inp"}'
    with open(input_file, 'w') as f:
        f.write(sample)
    sys.argv[1:] = [input_file]

    # when
    with patch('veloxchem.main.LoPropDriver', autospec=True):
        with patch('veloxchem.main.ScfRestrictedDriver', autospec=True) as mock_scf:
            main()

    # then
    mock_scf.assert_called_with(task.mpi_comm, task.ostream)
    mock_scf().compute.assert_called


@patch('veloxchem.loprop.h5py')
@patch('veloxchem.loprop.OverlapIntegralsDriver')
@patch('veloxchem.loprop.ao_matrix_to_dalton')
def test_overlap_called(mock_to_dalton, mock_ovldrv, mock_h5py):
    """
    Verify that veloxchem overlap driver from loprop module
    """

    # given
    mock_mpi = MagicMock()
    task = mock_mpi()
    task.input_dict = {
        'loprop': {'checkpoint_file': 'water.loprop.h5'},
        'scf': {'checkpoint_file': 'water.scf.h5'},
    }

    # when
    lpd = LoPropDriver(task)
    lpd.save_overlap()

    # then
    mock_ovldrv.assert_called_with(task.mpi_comm)
    mock_ovldrv().compute.assert_called_with(task.molecule, task.ao_basis)
    mock_h5py.File.assert_called_with('water.loprop.h5', 'w')
    mock_to_dalton.assert_called


@patch('veloxchem.main.MpiTask')
def test_save_coordinates(mock_mpi, sample, tmpdir):
    # Given
    task = mock_mpi()
    task.input_dict = {
        'loprop': {'checkpoint_file': f'{tmpdir}/water.loprop.h5'},
        'scf': {'checkpoint_file': f'{tmpdir}/water.scf.h5'},
    }
    task.molecule.x_to_numpy.return_value = [0., 0., 0.]
    task.molecule.y_to_numpy.return_value = [0., 1.4, -1.4]
    task.molecule.z_to_numpy.return_value = [0., 1.1, 1.1]

    # when
    loprop_driver = LoPropDriver(task)
    loprop_driver.save_coordinates()

    # then
    npt.assert_allclose(
        loprop_driver.get_coordinates(),
        [[0.0,  0.0,  0.0],
         [0.0,  1.4,  1.1],
         [0.0, -1.4,  1.1]]
    )


def test_input_dict(sample, tmpdir):
    """
    Verify that input parser sets a loprop key
    """
    # given
    input_file = f'{tmpdir/"water.inp"}'
    with open(input_file, 'w') as f:
        f.write(sample)

    # when
    ip = InputParser(input_file)

    # then
    assert ip.input_dict['jobs']['task'] == 'loprop'
    assert 'loprop' in ip.input_dict


@pytest.mark.skip('not there yet')
class TestIntegrations:

    @patch.object(LoPropDriver, 'get_cpa')
    @patch.object(LoPropDriver, 'get_opa')
    @patch.object(LoPropDriver, 'get_coordinates')
    def test_h2o_only_charges(self, get_coor, get_opa, get_cpa, capsys, tmpdir):
        get_cpa.return_value = (14, 5, 5)
        get_opa.return_value = ((0, 1, 3, 4, 5), (0,), (0,))
        get_coor.return_value = (
            (0.0,  0.0,  0.0),
            (0.0,  1.4,  1.1),
            (0.0, -1.4,  1.1),
            )
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
        input_file = f'{tmpdir/"water.inp"}'
        with open(input_file, 'w') as f:
            f.write(inp)

        sys.argv[1:] = [input_file]
        main()
        std = capsys.readouterr()

        expected = textwrap.dedent(
            """
            AU
            3 0 0 1
            1     0.000     0.000     0.000    -0.666
            1     0.000     1.400     1.100     0.333
            1     0.000    -1.400     1.100     0.333
            """
        )

        assert expected in std.out
