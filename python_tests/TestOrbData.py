from mpi4py import MPI
from veloxchem.taskparser import LocalTask
from veloxchem.veloxchemlib import denmat
from veloxchem.veloxchemlib import molorb
from veloxchem.veloxchemlib import mpi_master

from veloxchem.aodensitymatrix import AODensityMatrix
from veloxchem.molecularorbitals import MolecularOrbitals

import numpy as np
import unittest


class TestOrbData(unittest.TestCase):

    def test_get_label(self):

        task = LocalTask("inputs/dimer.inp", "inputs/dimer.out")
        basis = task.ao_basis

        self.assertEqual(basis.get_label(), "DEF2-SVP")

    def test_density_matrix(self):

        data_a = [[ 1., .2, ], [ .2, 1., ]]
        data_b = [[ .9, .5, ], [ .5, .9, ]]

        d_rest = AODensityMatrix.from_numpy_list([data_a], denmat.rest)
        d_unrest = AODensityMatrix.from_numpy_list([data_a, data_b],
                                                   denmat.unrest)

        den_a1 = d_rest.total_to_numpy(0)
        den_a2 = d_unrest.alpha_to_numpy(0)
        den_b1 = d_unrest.beta_to_numpy(0)

        self.assertEqual(0, np.max(np.abs(data_a - den_a1)))
        self.assertEqual(0, np.max(np.abs(data_a - den_a2)))
        self.assertEqual(0, np.max(np.abs(data_b - den_b1)))

        self.assertEqual(denmat.rest, d_rest.get_density_type())
        self.assertEqual(denmat.unrest, d_unrest.get_density_type())

    def test_density_hdf5(self):

        data_a = [[ 1., .2, ], [ .2, 1., ]]
        data_b = [[ .9, .5, ], [ .5, .9, ]]
        data_c = [[ .8, .6, ], [ .6, .8, ]]
        data_d = [[ .7, .5, ], [ .5, .7, ]]

        d_rest = AODensityMatrix.from_numpy_list(
            [data_a, data_b], denmat.rest)

        d_unrest = AODensityMatrix.from_numpy_list(
            [data_a, data_b, data_c, data_d], denmat.unrest)

        # hdf5 read/write tests

        if MPI.COMM_WORLD.Get_rank() == mpi_master():

            d_rest.write_hdf5("inputs/dummy.h5")
            dummy = AODensityMatrix.read_hdf5("inputs/dummy.h5")
            self.assertEqual(d_rest, dummy)

            d_unrest.write_hdf5("inputs/dummy.h5")
            dummy = AODensityMatrix.read_hdf5("inputs/dummy.h5")
            self.assertEqual(d_unrest, dummy)

    def test_orbitals_matrix(self):

        data_a = [[ .9, .2, ], [ .1, .3, ], [ .4, .9, ]]
        data_b = [[ .8, .3, ], [ .2, .4, ], [ .5, .8, ]]

        ener_a = [ 0.5, 0.9 ]
        ener_b = [ 0.3, 0.6 ]

        d_rest = MolecularOrbitals.from_numpy_list(
            [data_a], [ener_a], molorb.rest) 

        d_unrest = MolecularOrbitals.from_numpy_list(
            [data_a, data_b], [ener_a, ener_b], molorb.unrest)

        den_a1 = d_rest.alpha_to_numpy()
        den_a2 = d_unrest.alpha_to_numpy()
        den_b1 = d_unrest.beta_to_numpy()

        self.assertEqual(0, np.max(np.abs(data_a - den_a1)))
        self.assertEqual(0, np.max(np.abs(data_a - den_a2)))
        self.assertEqual(0, np.max(np.abs(data_b - den_b1)))

        ene_a1 = d_rest.ea_to_numpy()
        ene_a2 = d_unrest.ea_to_numpy()
        ene_b1 = d_unrest.eb_to_numpy()

        self.assertEqual(0, np.max(np.abs(ener_a - ene_a1)))
        self.assertEqual(0, np.max(np.abs(ener_a - ene_a2)))
        self.assertEqual(0, np.max(np.abs(ener_b - ene_b1)))

        self.assertEqual(molorb.rest, d_rest.get_orbitals_type())
        self.assertEqual(molorb.unrest, d_unrest.get_orbitals_type())

        # hdf5 read/write tests

        if MPI.COMM_WORLD.Get_rank() == mpi_master():

            d_rest.write_hdf5("inputs/dummy.h5")
            dummy = MolecularOrbitals.read_hdf5("inputs/dummy.h5")
            self.assertEqual(d_rest, dummy)

            d_unrest.write_hdf5("inputs/dummy.h5")
            dummy = MolecularOrbitals.read_hdf5("inputs/dummy.h5")
            self.assertEqual(d_unrest, dummy)


if __name__ == "__main__":
    unittest.main()
