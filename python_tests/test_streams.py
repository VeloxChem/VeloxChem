from mpi4py import MPI
from veloxchem.outputstream import OutputStream

import numpy as np
import unittest


class TestStreams(unittest.TestCase):

    def test_get_state(self):

        ostream = OutputStream("inputs/dummy.out")
        ostream.print_blank()
        self.assertTrue(ostream.get_state())

        ostream = OutputStream("")
        ostream.print_blank()
        self.assertFalse(ostream.get_state())


if __name__ == "__main__":
    unittest.main()
