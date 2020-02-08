import io
import textwrap
from unittest.mock import MagicMock

from veloxchem.mpitask import MpiTask
from mpi4py import MPI


def test_loprop_main(tmpdir):
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
    out = io.StringIO()

    with open(tmpdir/'water.inp', 'w'):

    task = MagicMock()  # MpiTask((inp, out), MPI.COMM_WORLD)
    molecule = MagicMock()  # task.molecule
    basis = MagicMock()  # task.ao_basis
    comm = MagicMock()  # task.mpi_comm
    rank = MagicMock()  # task.mpi_rank


    with patch('veloxchem.LoPropDriver') as mock_prop:
        main()

    assert mock_prop.called
