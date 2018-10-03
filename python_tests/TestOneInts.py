from mpi4py import MPI
from HelperClass import Task
from VeloxChemLib import OverlapMatrix
from VeloxChemLib import KineticEnergyMatrix
from VeloxChemLib import NuclearPotentialMatrix

import numpy as np
import unittest


class TestOneInts(unittest.TestCase):

    def test_overlap_matrix(self):

        data = [[ 1., .2, ], [ .2, 1., ]]

        array = np.array(data)
        matrix = OverlapMatrix.from_numpy(array)
        array2 = matrix.to_numpy()
        matrix2 = OverlapMatrix.from_numpy(array2)

        self.assertEqual(0, np.max(np.abs(array - array2)))
        self.assertEqual(matrix, matrix2)

    def test_kinetic_energy_matrix(self):

        data = [[ 1., .2, ], [ .2, 1., ]]

        array = np.array(data)
        matrix = KineticEnergyMatrix.from_numpy(array)
        array2 = matrix.to_numpy()
        matrix2 = KineticEnergyMatrix.from_numpy(array2)

        self.assertEqual(0, np.max(np.abs(array - array2)))
        self.assertEqual(matrix, matrix2)

    def test_nuclear_potential_matrix(self):

        data = [[ 1., .2, ], [ .2, 1., ]]

        array = np.array(data)
        matrix = NuclearPotentialMatrix.from_numpy(array)
        array2 = matrix.to_numpy()
        matrix2 = NuclearPotentialMatrix.from_numpy(array2)

        self.assertEqual(0, np.max(np.abs(array - array2)))
        self.assertEqual(matrix, matrix2)


if __name__ == "__main__":
    unittest.main()
