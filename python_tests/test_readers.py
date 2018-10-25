from mpi4py import MPI
from veloxchem.taskparser import LocalTask

import numpy as np
import unittest


class TestReaders(unittest.TestCase):

    def test_get_state(self):

        task = LocalTask("inputs/water.inp", "inputs/water.out")
        self.assertTrue(task.xyz_reader.get_state())
        self.assertTrue(task.env_reader.get_state())
        self.assertTrue(task.basis_reader.get_state())


if __name__ == "__main__":
    unittest.main()
