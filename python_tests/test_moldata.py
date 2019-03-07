from mpi4py import MPI
from veloxchem.mpitask import MpiTask
from veloxchem.molecule import Molecule
from veloxchem.veloxchemlib import ChemicalElement
from veloxchem.veloxchemlib import bohr_in_angstroms

import numpy as np
import unittest


class TestMolData(unittest.TestCase):

    def nh3_labels(self):

        return ["N", "H", "H", "H"]

    def nh3_coords(self):

        return [[-3.710, 3.019, -0.037], [-3.702, 4.942, 0.059],
                [-4.704, 2.415, 1.497], [-4.780, 2.569, -1.573]]

    def nh3_xyzstr(self):

        return """N  -3.710   3.019  -0.037
                  H  -3.702   4.942   0.059
                  H  -4.704   2.415   1.497
                  H  -4.780   2.569  -1.573"""

    def nh3_molecule(self):

        labels = self.nh3_labels()
        coords = self.nh3_coords()

        return Molecule(labels, coords, 'au')

    def test_constructors(self):

        labels = self.nh3_labels()
        coords = self.nh3_coords()
        xyzstr = self.nh3_xyzstr()

        mol_1 = Molecule(labels, coords, 'au')
        mol_2 = Molecule.read_str(xyzstr, 'au')

        array = np.array(coords)
        arrayT = np.zeros((array.shape[1], array.shape[0]))
        for i in range(array.shape[0]):
            for j in range(array.shape[1]):
                arrayT[j][i] = array[i][j]

        mol_3 = Molecule(labels, array, 'au')
        mol_4 = Molecule(labels, arrayT.T, 'au')

        array_angs = array * bohr_in_angstroms()

        mol_5 = Molecule(labels, array_angs)
        mol_6 = Molecule(labels, array_angs, 'angs')

        mol_7 = Molecule([7, 1, 1, 1], array, 'au')
        mol_8 = Molecule([7, 1, 1, 1], array_angs, 'angs')
        mol_9 = Molecule([7, 1, 1, 1], array_angs)

        self.assertEqual(mol_1, mol_2)
        self.assertEqual(mol_1, mol_3)
        self.assertEqual(mol_1, mol_4)
        self.assertEqual(mol_1, mol_5)
        self.assertEqual(mol_1, mol_6)
        self.assertEqual(mol_1, mol_7)
        self.assertEqual(mol_1, mol_8)
        self.assertEqual(mol_1, mol_9)

    def test_get_sub_molecule(self):

        task = MpiTask(["inputs/dimer.inp", "inputs/dimer.out"], MPI.COMM_WORLD)
        molecule = task.molecule

        mol_1 = molecule.get_sub_molecule(0, 4)
        mol_2 = molecule.get_sub_molecule(0, 9)
        mol_3 = mol_2.get_sub_molecule(0, 4)

        mol_4 = self.nh3_molecule()

        self.assertEqual(mol_1, mol_3)
        self.assertEqual(mol_1, mol_4)
        self.assertEqual(mol_2, molecule)

    def test_coordinates_to_numpy(self):

        mol = self.nh3_molecule()

        x = mol.x_to_numpy()
        y = mol.y_to_numpy()
        z = mol.z_to_numpy()

        x_arr = np.array([-3.710, -3.702, -4.704, -4.780])
        y_arr = np.array([3.019, 4.942, 2.415, 2.569])
        z_arr = np.array([-0.037, 0.059, 1.497, -1.573])

        self.assertTrue((x == x_arr).all())
        self.assertTrue((y == y_arr).all())
        self.assertTrue((z == z_arr).all())

    def test_setters_and_getters(self):

        mol = self.nh3_molecule()

        mol.set_charge(-1)
        self.assertEqual(-1.0, mol.get_charge())

        mol.set_multiplicity(2)
        self.assertEqual(2.0, mol.get_multiplicity())

        mol.check_multiplicity()
        mol.check_proximity(0.1)

        elem_comp = mol.get_elemental_composition()
        self.assertTrue(elem_comp == [1, 7])

    def test_number_of_atoms(self):

        mol = self.nh3_molecule()

        self.assertEqual(mol.number_of_atoms(1), 3)
        self.assertEqual(mol.number_of_atoms(7), 1)

        self.assertEqual(mol.number_of_atoms(0, 1, 1), 0)
        self.assertEqual(mol.number_of_atoms(0, 2, 1), 1)

    def test_vdw_radii_and_elem_ids(self):

        # fake molecule made of H,Li,C,N,O,S,Cu,Zn,Br,Ag,Au,Hg

        mol = Molecule.read_str(
            """H    0.0   0.0   0.0
               Li   0.0   0.0   1.0
               C    0.0   0.0   2.0
               N    0.0   0.0   3.0
               O    0.0   0.0   4.0
               S    0.0   0.0   5.0
               Cu   0.0   0.0   6.0
               Zn   0.0   0.0   7.0
               Br   0.0   0.0   8.0
               Ag   0.0   0.0   9.0
               Au   0.0   0.0  10.0
               Hg   0.0   0.0  11.0""", 'angs', 0, 2)

        atom_radii = mol.vdw_radii_to_numpy()

        ref_radii = np.array([
            1.09, 1.82, 1.70, 1.55, 1.52, 1.80, 1.40, 1.39, 1.85, 1.72, 1.66,
            1.55
        ])

        ref_radii /= bohr_in_angstroms()

        self.assertTrue((atom_radii == ref_radii).all())

        elem_ids = mol.elem_ids_to_numpy()

        ref_ids = np.array([1, 3, 6, 7, 8, 16, 29, 30, 35, 47, 79, 80])

        self.assertTrue((elem_ids == ref_ids).all())

    def test_chemical_element(self):

        elem = ChemicalElement()
        self.assertEqual("", elem.get_name())
        elem.set_atom_type("BR")
        self.assertEqual("Br", elem.get_name())

        elem2 = ChemicalElement()
        elem2.set_atom_type(35)
        self.assertEqual("Br", elem2.get_name())

        self.assertEqual(elem, elem2)


if __name__ == "__main__":
    unittest.main()
