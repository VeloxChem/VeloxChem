from mpi4py import MPI
from veloxchem.mpitask import MpiTask
from veloxchem.veloxchemlib import Molecule
from veloxchem.veloxchemlib import bohr_in_angstroms

import numpy as np
import unittest


class TestMolData(unittest.TestCase):

    def test_get_sub_molecule(self):

        task = MpiTask("inputs/dimer.inp", "inputs/dimer.out", MPI.COMM_WORLD)
        molecule = task.molecule

        mol_1 = molecule.get_sub_molecule(0, 4)
        mol_2 = molecule.get_sub_molecule(0, 9)
        mol_3 = mol_2.get_sub_molecule(0, 4)

        mol_4 = Molecule.from_xyz(
            [ "N", "H", "H", "H" ],
            [ -3.710, -3.702, -4.704, -4.780 ],
            [  3.019,  4.942,  2.415,  2.569 ],
            [ -0.037,  0.059,  1.497, -1.573 ])

        self.assertEqual(mol_1, mol_3)
        self.assertEqual(mol_1, mol_4)
        self.assertEqual(mol_2, molecule)

    def test_coordinates_to_numpy(self):

        x_list = [ -3.710, -3.702, -4.704, -4.780 ]
        y_list = [  3.019,  4.942,  2.415,  2.569 ]
        z_list = [ -0.037,  0.059,  1.497, -1.573 ]

        mol = Molecule.from_xyz([ "N", "H", "H", "H" ], x_list, y_list, z_list)

        x_arr = np.array(x_list)
        y_arr = np.array(y_list)
        z_arr = np.array(z_list)

        x = mol.x_to_numpy()
        y = mol.y_to_numpy()
        z = mol.z_to_numpy()

        self.assertTrue((x == x_arr).all())
        self.assertTrue((y == y_arr).all())
        self.assertTrue((z == z_arr).all())

    def test_number_of_atoms(self):

        mol = Molecule.from_xyz(
            [ "N", "H", "H", "H" ],
            [ -3.710, -3.702, -4.704, -4.780 ],
            [  3.019,  4.942,  2.415,  2.569 ],
            [ -0.037,  0.059,  1.497, -1.573 ])

        self.assertEqual(mol.number_of_atoms(1), 3)
        self.assertEqual(mol.number_of_atoms(7), 1)

        self.assertEqual(mol.number_of_atoms(0, 1, 1), 0)
        self.assertEqual(mol.number_of_atoms(0, 2, 1), 1)

    def test_vdw_radii_and_elem_ids(self):

        # fake molecule made of H,Li,C,N,O,S,Cu,Zn,Br,Ag,Au,Hg

        mol = Molecule.from_xyz(
            [ "H", "Li", "C", "N", "O", "S", "Cu", "Zn",
              "Br", "Ag", "Au", "Hg", ],
            [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0 ],
            [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0 ],
            [ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0 ])

        atom_radii = mol.vdw_radii_to_numpy()

        ref_radii = np.array([ 1.09, 1.82, 1.70, 1.55, 1.52, 1.80, 1.40, 1.39,
                               1.85, 1.72, 1.66, 1.55 ])

        ref_radii /= bohr_in_angstroms()

        self.assertTrue((atom_radii == ref_radii).all())

        elem_ids = mol.elem_ids_to_numpy()

        ref_ids = np.array([ 1, 3, 6, 7, 8, 16, 29, 30, 35, 47, 79, 80 ])

        self.assertTrue((elem_ids == ref_ids).all())


if __name__ == "__main__":
    unittest.main()
