from mpi4py import MPI
from VeloxChemMP import *

import numpy as np
import unittest

class TestControllers(unittest.TestCase):

    def test_app_manager_0(self):

        app = AppManager.create("dummy.inp", "")

        app.execute()

        self.assertFalse(app.get_state())

    def test_app_manager_1(self):

        app = AppManager.create("inputs/test.inp", "inputs/test.out")

        app.execute()

        self.assertTrue(app.get_state())

if __name__ == "__main__":
    unittest.main()
