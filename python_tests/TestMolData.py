from mpi4py import MPI
from HelperClass import Task
from veloxchem.VeloxChemLib import Molecule

import numpy as np
import unittest


class TestMolData(unittest.TestCase):

    def test_get_sub_molecule(self):

        task = Task("inputs/dimer.inp", "inputs/dimer.out")
        molecule = task.molecule

        mol_1 = molecule.get_sub_molecule(0, 4)
        mol_2 = molecule.get_sub_molecule(0, 9)
        mol_3 = mol_2.get_sub_molecule(0, 4)
        mol_4 = Molecule.from_list(
            [ -3.710, -3.702, -4.704, -4.780,
               3.019,  4.942,  2.415,  2.569,
              -0.037,  0.059,  1.497, -1.573 ],
            [ 7.0, 1.0, 1.0, 1.0 ],
            [ 14.003074, 1.007825, 1.007825, 1.007825 ],
            [ "N", "H", "H", "H" ],
            [ 7, 1, 1, 1 ])

        self.assertEqual(mol_1, mol_3)
        self.assertEqual(mol_1, mol_4)
        self.assertEqual(mol_2, molecule)


if __name__ == "__main__":
    unittest.main()
