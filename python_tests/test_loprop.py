import os
import sys
import textwrap
from unittest.mock import patch, MagicMock, mock_open

import pytest
import numpy.testing as npt
from mpi4py import MPI

from veloxchem.main import main
from veloxchem.mpitask import MpiTask
from veloxchem.inputparser import InputParser, InputError
from veloxchem.loprop import (
    LoPropDriver,
    count_contracted,
    count_contracted_on_atom,
    get_basis_file,
)


@pytest.fixture
def sample():
    inp = textwrap.dedent(
        """
        @jobs
        task: loprop
        @end

        @method settings
        basis: def2-svp
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
        with patch(
            'veloxchem.main.ScfRestrictedDriver', autospec=True
        ) as mock_scf:
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


def test_input_settings(tmpdir):
    """
    Verify loprop options
    """
    input_file = f'{tmpdir/"water.inp"}'
    input_content = textwrap.dedent(
        """
        @jobs
        task: 'loprop'
        @end

        @loprop
        localize: charges
        @end
        """
    )
    with open(input_file, 'w') as f:
        f.write(input_content)

    # when
    ip = InputParser(input_file)
    assert ip.input_dict['loprop']['localize'] == 'charges'


def test_wrong_input(tmpdir):
    """
    Verify loprop options
    """
    input_file = f'{tmpdir/"water.inp"}'
    input_content = textwrap.dedent(
        """
        @jobs
        task: loprop
        @end

        @loprop
        localize: notimplemented
        @end
        """
    )
    with open(input_file, 'w') as f:
        f.write(input_content)

    # when
    with pytest.raises(InputError) as nie:
        InputParser(input_file)
    assert nie.value.args == ('localize: notimplemented illegal value',)


def test_cpa(sample, tmpdir):
    """
    Verify count of orbitals per atom
    """
    # given
    input_file = f'{tmpdir/"water.inp"}'
    with open(input_file, 'w') as f:
        f.write(sample)

    # when
    task = MpiTask([input_file], MPI.COMM_WORLD)
    loprop_driver = LoPropDriver(task)

    # then

    assert loprop_driver.get_cpa() == [14, 5, 5]


def test_opa(sample, tmpdir):
    """
    Verify list of atomic orbitals considered occupied
    """

    # given
    input_file = f'{tmpdir/"water.inp"}'
    with open(input_file, 'w') as f:
        f.write(sample)

    # when
    task = MpiTask([input_file], MPI.COMM_WORLD)
    loprop_driver = LoPropDriver(task)

    # then
    assert loprop_driver.get_opa() == [[0, 1, 3, 4, 5], [0], [0]]


@pytest.mark.parametrize(
    'input, expected',
    [
        (
            textwrap.dedent(
                """
                @ATOMBASIS H
                S 3  1
                1.301070100000e+01  1.968215800000e-02
                1.962257200000e+00  1.379652400000e-01
                4.445379600000e-01  4.783193500000e-01
                S 1  1
                1.219496200000e-01  1.000000000000e+00
                P 1  1
                8.000000000000e-01  1.000000000000e+00
                @END
                """
            ),
            {'S': 2, 'P': 1}
        ),
        (
            textwrap.dedent(
                """
                @ATOMBASIS O
                S 5  1
                2.266176778500e+03 -5.343180992600e-03
                3.408701019100e+02 -3.989003923000e-02
                7.736313516700e+01 -1.785391198500e-01
                2.147964494000e+01 -4.642768495900e-01
                6.658943312400e+00 -4.430974517200e-01
                S 1  1
                8.097597566800e-01  1.000000000000e+00
                S 1  1
                2.553077223400e-01  1.000000000000e+00
                P 3  1
                1.772150431700e+01  4.339457319300e-02
                3.863550544000e+00  2.309412076500e-01
                1.048092088300e+00  5.137531106400e-01
                P 1  1
                2.764154441100e-01  1.000000000000e+00
                D 1  1
                1.200000000000e+00  1.000000000000e+00
                @END
                """
            ),
            {'S': 3, 'P': 2, 'D': 1}
        ),
    ],
    ids=['H:2S1P', 'O:3S2P1D']
)
def test_count_contracted(input, expected):
    assert count_contracted(input.split('\n')) == expected


@pytest.mark.parametrize(
    'input, expected',
    [
        ({'S': 2, 'P': 1}, 5),
    ]
)
def test_count_contracted_on_atom(input, expected):
    assert count_contracted_on_atom(input) == expected


def test_get_local_basis_file():
    basis = 'STO-3G'

    with patch('veloxchem.loprop.os.path') as mock_path:
        mock_path.exists.return_value = True
        mock_path.abspath.return_value = '/somepath/STO-3G'
        full_path = get_basis_file(basis)

    assert full_path == '/somepath/STO-3G'


def test_get_lib_basis_file():

    with patch('veloxchem.loprop.os.path.exists') as mock_exists:
        with patch.dict(os.environ, {'VLXBASISPATH': '/vlxlib'}):
            mock_exists.side_effect = [False, True]
            full_path = get_basis_file('STO-3G')

    assert full_path == '/vlxlib/STO-3G'


class TestIntegrations:

    def test_h2o_only_charges(self, capsys, tmpdir):
        inp = textwrap.dedent(
            """
            @jobs
            task: loprop
            @end

            @method settings
            basis: def2-svp
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
            1     0.000     0.000     0.000     0.060
            1     0.000     1.400     1.100    -0.030
            1     0.000    -1.400     1.100    -0.030
            """

        )

        assert expected in std.out
