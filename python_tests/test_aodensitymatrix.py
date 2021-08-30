from pathlib import Path

from veloxchem.veloxchemlib import denmat, is_mpi_master
from veloxchem.aodensitymatrix import AODensityMatrix


class TestAODensityMatrix:
    """
    Implements tests for src/pymodule/aodensitymatrix.py
    """

    def test_aodensitymatrix(self):

        data_a = [[1., .2], [.2, 1.]]
        data_b = [[.9, .5], [.5, .9]]
        data_c = [[.8, .6], [.6, .8]]
        data_d = [[.7, .5], [.5, .7]]

        d_rest = AODensityMatrix([data_a, data_b], denmat.rest)

        d_unrest = AODensityMatrix([data_a, data_b, data_c, data_d],
                                   denmat.unrest)

        d_osrest = AODensityMatrix([data_a, data_b, data_c, data_d],
                                   denmat.osrest)

        if is_mpi_master():

            here = Path(__file__).parent
            h5file = str(here / 'inputs' / 'dummy.h5')

            d_rest.write_hdf5(h5file)
            dummy = AODensityMatrix.read_hdf5(h5file)
            assert d_rest == dummy

            d_unrest.write_hdf5(h5file)
            dummy = AODensityMatrix.read_hdf5(h5file)
            assert d_unrest == dummy

            d_osrest.write_hdf5(h5file)
            dummy = AODensityMatrix.read_hdf5(h5file)
            assert d_osrest == dummy
