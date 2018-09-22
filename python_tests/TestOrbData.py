from mpi4py import MPI
from VeloxChemMP import *
from HelperClass import *

import numpy as np
import unittest

class TestOrbData(unittest.TestCase):

    def print_title(self, label):

        print("\n[ Running ] TestOrbData " + label)

    def test_get_label(self):

        self.print_title("get_label")

        task = Task("inputs/dimer.inp", "inputs/dimer.out")

        basis = task.ao_basis

        self.assertEqual(basis.get_label(), "DEF2-SVP")

    def test_density_matrix(self):

        self.print_title("density_matrix")

        data_a = [[1., .2,], [.2, 1.,]]

        data_b = [[1., .5,], [.5, 1.,]]

        array_a = np.array(data_a)

        array_b = np.array(data_b)

        d_rest = AODensityMatrix.from_numpy_list([data_a], True)

        d_unrest = AODensityMatrix.from_numpy_list([data_a, data_b], False)

        den_a1 = d_rest.total_to_numpy(0)

        den_a2 = d_unrest.alpha_to_numpy(0)

        den_b1 = d_unrest.beta_to_numpy(0)

        self.assertEqual(0, np.max(np.abs(data_a-den_a1)))

        self.assertEqual(0, np.max(np.abs(data_a-den_a2)))

        self.assertEqual(0, np.max(np.abs(data_b-den_b1)))
