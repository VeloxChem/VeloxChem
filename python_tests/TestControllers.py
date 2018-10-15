from mpi4py import MPI
from veloxchem.VeloxChemLib import AppManager

import numpy as np
import unittest


class TestControllers(unittest.TestCase):

    def test_app_manager_0(self):

        app = AppManager("inputs/dummy.inp", "")
        app.execute()

        self.assertFalse(app.get_state())

    def test_app_manager_1(self):

        app = AppManager("inputs/water.inp", "inputs/water.out")
        app.execute()

        self.assertTrue(app.get_state())


if __name__ == "__main__":
    unittest.main()
