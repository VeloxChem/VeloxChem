from mpi4py import MPI
from VeloxChemMP import *
from HelperClass import *

import numpy as np
import unittest

class TestReaders(unittest.TestCase):

    def print_title(self, label):

        print("\n[ Running ] TestReaders " + label)

    def test_get_state(self):

        self.print_title("get_state")

        task = Task("inputs/water.inp", "inputs/water.out")

        self.assertTrue(task.xyz_reader.get_state())

        self.assertTrue(task.env_reader.get_state())

        self.assertTrue(task.basis_reader.get_state())
