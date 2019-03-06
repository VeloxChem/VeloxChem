from mpi4py import MPI
from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import mpi_initialized
from veloxchem.veloxchemlib import bohr_in_angstroms
from veloxchem.veloxchemlib import hartree_in_ev
from veloxchem.veloxchemlib import to_angular_momentum
from veloxchem.errorhandler import assert_msg_critical

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

    def test_angular_momentum(self):

        angmoms = "SPDFGH"

        for ind, ang in enumerate(angmoms):
            self.assertEqual(to_angular_momentum(ind), ang)
            self.assertEqual(to_angular_momentum(ang), ind)


if __name__ == "__main__":
    unittest.main()
