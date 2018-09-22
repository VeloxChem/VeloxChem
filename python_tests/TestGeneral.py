from mpi4py import MPI
from VeloxChemMP import *

import numpy as np
import unittest

class TestGeneral(unittest.TestCase):

    def print_title(self, label):

        print("\n[ Running ] TestGeneral " + label)

    def test_mpi_master(self):

        self.print_title("mpi_master")

        self.assertEqual(0, mpi_master())

    def test_mpi_initialized(self):

        self.print_title("mpi_initialized")

        self.assertTrue(mpi_initialized())

    def test_assert_msg_critical(self):

        self.print_title("assert_msg_critical")

        assert_msg_critical(True, "")

        self.assertTrue(True)
