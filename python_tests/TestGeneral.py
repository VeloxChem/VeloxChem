from mpi4py import MPI
from veloxchem.VeloxChemLib import mpi_master
from veloxchem.VeloxChemLib import mpi_initialized
from veloxchem.VeloxChemLib import assert_msg_critical

import numpy as np
import unittest


class TestGeneral(unittest.TestCase):

    def test_mpi_master(self):

        self.assertEqual(0, mpi_master())

    def test_mpi_initialized(self):

        self.assertTrue(mpi_initialized())

    def test_assert_msg_critical(self):

        assert_msg_critical(True, "")

        self.assertTrue(True)


if __name__ == "__main__":
    unittest.main()
