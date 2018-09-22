from mpi4py import MPI
from VeloxChemMP import *
from HelperClass import *
import numpy as np

import unittest

class TestStreams(unittest.TestCase):

    def print_title(self, label):

        print()
        print("[ Running ] TestStreams " + label)

    def test_get_state(self):

        self.print_title("get_state")

        ostream = OutputStream("")
        istream = InputStream("inputs/water.inp", ostream)

        self.assertTrue(istream.get_state())
        self.assertTrue(ostream.get_state())
