from mpi4py import MPI
from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import mpi_initialized
from veloxchem.veloxchemlib import assert_msg_critical
from veloxchem.veloxchemlib import bohr_in_angstroms
from veloxchem.veloxchemlib import hartree_in_ev

import numpy as np
import unittest


class TestGeneral(unittest.TestCase):

    def test_mpi_master(self):

        self.assertEqual(0, mpi_master())

    def test_mpi_initialized(self):

        self.assertTrue(mpi_initialized())

    def test_assert_msg_critical(self):

        assert_msg_critical(True, "")

    def test_constants(self):

        self.assertAlmostEqual(bohr_in_angstroms(), 0.529177, 6)
        self.assertAlmostEqual(hartree_in_ev(), 27.2114, 4)


if __name__ == "__main__":
    unittest.main()
