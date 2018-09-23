from mpi4py import MPI
from VeloxChemMP import *
from HelperClass import *

import numpy as np
import unittest

class TestMolData(unittest.TestCase):

    def test_get_sub_molecule(self):

        task = Task("inputs/dimer.inp", "inputs/dimer.out")

        molecule = task.molecule

        mol_1 = molecule.get_sub_molecule(0,4)

        mol_2 = molecule.get_sub_molecule(0,9)

        mol_3 = mol_2.get_sub_molecule(0,4)

        self.assertEqual(mol_1, mol_3)

        self.assertEqual(mol_2, molecule)

if __name__ == "__main__":
    unittest.main()
