from pathlib import Path

from veloxchem.veloxchemlib import fockmat, is_mpi_master
from veloxchem.aofockmatrix import AOFockMatrix


class TestAOFockMatrix:
    """
    Implements tests for src/pymodule/aofockmatrix.py
    """

    def test_aofockmatrix(self):

        data_a = [[1., .2], [.2, 1.]]
        data_b = [[.9, .5], [.5, .9]]
        data_c = [[.8, .6], [.6, .8]]
        data_d = [[.7, .5], [.5, .7]]

        f_rest = AOFockMatrix(
            [data_a, data_b, data_c, data_d],
            [fockmat.restjk, fockmat.restj, fockmat.rgenj, fockmat.rgenk],
            [1.0, 2.0, 0.5, 1.2],
            [0, 1, 2, 3],
        )

        f_unrest = AOFockMatrix(
            [data_a, data_b, data_c, data_d],
            [
                fockmat.unrestjk, fockmat.unrestj, fockmat.unrestjkx,
                fockmat.unrestj
            ],
            [1.0, 2.0, 0.2, 1.5],
            [0, 0, 1, 1],
        )

        # hdf5 read/write tests

        if is_mpi_master():

            here = Path(__file__).parent
            h5file = here / 'inputs' / 'dummy.h5'

            f_rest.write_hdf5(h5file)
            dummy = AOFockMatrix.read_hdf5(h5file)
            assert f_rest == dummy

            f_unrest.write_hdf5(h5file)
            dummy = AOFockMatrix.read_hdf5(h5file)
            assert f_unrest == dummy
