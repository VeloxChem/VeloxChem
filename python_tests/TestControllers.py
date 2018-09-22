from mpi4py import MPI
from VeloxChemMP import *

import numpy as np
import unittest

class TestControllers(unittest.TestCase):

    def print_title(self, label):

        print("\n[ Running ] TestControllers " + label)

    def test_app_manager(self):

        self.print_title("app_manager")

        app = AppManager.create("dummy.inp", "")

        app.execute()

        self.assertFalse(app.get_state())
