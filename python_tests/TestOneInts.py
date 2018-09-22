from mpi4py import MPI
from VeloxChemMP import *
from HelperClass import *

import numpy as np
import unittest

class TestOneInts(unittest.TestCase):

    def print_title(self, label):

        print("\n[ Running ] TestOneInts " + label)

    def test_overlap_matrix(self):

        self.print_title("overlap_matrix")

        data = [[1., .2,], [.2, 1.,]]

        array = np.array(data)

        matrix = OverlapMatrix.from_numpy(array)

        array2 = matrix.to_numpy()

        matrix2 = OverlapMatrix.from_numpy(array2)

        self.assertEqual(0, np.max(np.abs(array-array2)))

        self.assertEqual(matrix, matrix2)

    def test_kinetic_energy_matrix(self):

        self.print_title("kinetic_energy_matrix")

        data = [[1., .2,], [.2, 1.,]]

        array = np.array(data)

        matrix = KineticEnergyMatrix.from_numpy(array)

        array2 = matrix.to_numpy()

        matrix2 = KineticEnergyMatrix.from_numpy(array2)

        self.assertEqual(0, np.max(np.abs(array-array2)))

        self.assertEqual(matrix, matrix2)

    def test_nuclear_potential_matrix(self):

        self.print_title("nuclear_potential_matrix")

        data = [[1., .2,], [.2, 1.,]]

        array = np.array(data)

        matrix = NuclearPotentialMatrix.from_numpy(array)

        array2 = matrix.to_numpy()

        matrix2 = NuclearPotentialMatrix.from_numpy(array2)

        self.assertEqual(0, np.max(np.abs(array-array2)))

        self.assertEqual(matrix, matrix2)
