from mpi4py import MPI
from veloxchem.veloxchemlib import OutputStream

import numpy as np
import unittest


class TestStreams(unittest.TestCase):

    def test_get_state(self):

        ostream = OutputStream("")

        self.assertTrue(ostream.get_state())


if __name__ == "__main__":
    unittest.main()
