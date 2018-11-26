from mpi4py import MPI
from veloxchem.taskparser import LocalTask
from veloxchem.veloxchemlib import Molecule

import numpy as np
import unittest


class TestMolData(unittest.TestCase):

    def test_get_sub_molecule(self):

        task = LocalTask("inputs/dimer.inp", "inputs/dimer.out")
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

    def test_coordinates_to_numpy(self):

        mol = Molecule.from_list(
            [ -3.710, -3.702, -4.704, -4.780,
               3.019,  4.942,  2.415,  2.569,
              -0.037,  0.059,  1.497, -1.573 ],
            [ 7.0, 1.0, 1.0, 1.0 ],
            [ 14.003074, 1.007825, 1.007825, 1.007825 ],
            [ "N", "H", "H", "H" ],
            [ 7, 1, 1, 1 ])

        ref_coords = np.array(
            [[ -3.710, -3.702, -4.704, -4.780 ],
             [  3.019,  4.942,  2.415,  2.569 ],
             [ -0.037,  0.059,  1.497, -1.573 ]])

        mol_coords = mol.coordinates_to_numpy()

        self.assertEqual(mol_coords.shape[0], 3)
        self.assertEqual(mol_coords.shape[1], mol.number_of_atoms())

        self.assertTrue((mol_coords == ref_coords).all())

    def test_number_of_atoms(self):

        mol = Molecule.from_list(
            [ -3.710, -3.702, -4.704, -4.780,
               3.019,  4.942,  2.415,  2.569,
              -0.037,  0.059,  1.497, -1.573 ],
            [ 7.0, 1.0, 1.0, 1.0 ],
            [ 14.003074, 1.007825, 1.007825, 1.007825 ],
            [ "N", "H", "H", "H" ],
            [ 7, 1, 1, 1 ])

        self.assertEqual(mol.number_of_atoms(1), 3)
        self.assertEqual(mol.number_of_atoms(7), 1)

        self.assertEqual(mol.number_of_atoms(0, 1, 1), 0)
        self.assertEqual(mol.number_of_atoms(0, 2, 1), 1)


if __name__ == "__main__":
    unittest.main()
