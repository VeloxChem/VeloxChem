from mpi4py import MPI
from pathlib import Path
import numpy as np
import unittest
import tempfile

from veloxchem.veloxchemlib import bohr_in_angstroms
from veloxchem.veloxchemlib import ChemicalElement
from veloxchem.mpitask import MpiTask
from veloxchem.molecule import Molecule


class TestMolData(unittest.TestCase):

    def nh3_labels(self):

        return ['N', 'H', 'H', 'H']

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

        array_ang = array * bohr_in_angstroms()

        mol_5 = Molecule(labels, array_ang)
        mol_6 = Molecule(labels, array_ang, 'angstrom')

        mol_7 = Molecule([7, 1, 1, 1], array, 'au')
        mol_8 = Molecule([7, 1, 1, 1], array_ang, 'angstrom')
        mol_9 = Molecule([7, 1, 1, 1], array_ang)

        self.assertEqual(mol_1, mol_2)
        self.assertEqual(mol_1, mol_3)
        self.assertEqual(mol_1, mol_4)
        self.assertEqual(mol_1, mol_5)
        self.assertEqual(mol_1, mol_6)
        self.assertEqual(mol_1, mol_7)
        self.assertEqual(mol_1, mol_8)
        self.assertEqual(mol_1, mol_9)

    def test_get_sub_molecule(self):

        here = Path(__file__).parent
        inpfile = here / 'inputs' / 'dimer.inp'
        outfile = inpfile.with_suffix('.out')

        task = MpiTask([str(inpfile), str(outfile)], MPI.COMM_WORLD)
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

    def test_center_of_mass(self):

        mol = self.nh3_molecule()
        mol_com = mol.center_of_mass()
        ref_com = np.array([-3.831697, 3.070437, -0.031436])
        self.assertTrue(np.max(np.abs(ref_com - mol_com)) < 1.0e-6)

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

        mol = Molecule.read_str("""H    0.0   0.0   0.0
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
                                   Hg   0.0   0.0  11.0""")

        atom_radii = mol.vdw_radii_to_numpy()

        ref_radii = np.array([
            1.09, 1.82, 1.70, 1.55, 1.52, 1.80, 1.40, 1.39, 1.85, 1.72, 1.66,
            1.55
        ])

        ref_radii /= bohr_in_angstroms()

        self.assertTrue((atom_radii == ref_radii).all())

        elem_ids = mol.elem_ids_to_numpy()

        self.assertTrue(elem_ids.dtype.type is np.int32)

        ref_ids = np.array([1, 3, 6, 7, 8, 16, 29, 30, 35, 47, 79, 80])

        self.assertTrue((elem_ids == ref_ids).all())

    def test_chemical_element(self):

        elem = ChemicalElement()
        self.assertEqual('', elem.get_name())
        elem.set_atom_type('BR')
        self.assertEqual('Br', elem.get_name())

        elem2 = ChemicalElement()
        elem2.set_atom_type(35)
        self.assertEqual('Br', elem2.get_name())

        self.assertEqual(elem, elem2)

    def test_get_ic_rmsd(self):

        init_xyz_str = """
            C          -2.723479309037        0.388400197499       -0.032041680121
            C          -1.198591674506        0.292402159666       -0.001957955680
            H          -3.161468702971        0.002148459476        0.916810101463
            H          -3.113340035638       -0.249559972670       -0.853012496915
            O          -3.131199766507        1.709538730510       -0.265228909885
            H          -0.880032285304       -0.767190536848        0.124520689502
            O          -0.670239980035        1.084030100114        1.027238801838
            H          -0.785679988597        0.645580322006       -0.970638007358
            H          -0.821288323120        0.587492539563        1.873631291844
            H          -3.112122928053        2.174373284354        0.611781773545
        """

        final_xyz_str = """
            C          -2.760169377434        0.322719976286       -0.001925407564
            C          -1.228941816065        0.304680448557       -0.002845411826
            H          -3.128135202340       -0.132424338111        0.933568102246
            H          -3.149271541311       -0.251138690089       -0.845341019529
            O          -3.252173881237        1.631288380215       -0.127670008266
            H          -0.863067365090       -0.729220194991        0.037102867132
            O          -0.717621604087        1.071775860470        1.064783463761
            H          -0.854080755853        0.782453035355       -0.910452256526
            H          -0.974559560638        0.659440589631        1.897920584004
            H          -2.669421889714        2.207640216349        0.385962694800
        """

        init_mol = Molecule.read_str(init_xyz_str, units='angstrom')
        final_mol = Molecule.read_str(final_xyz_str, units='angstrom')

        ic_rmsd = final_mol.get_ic_rmsd(init_mol)
        ref_ic_rmsd = {
            'bonds': {
                'rms': 0.017,
                'max': 0.028,
                'unit': 'Angstrom'
            },
            'angles': {
                'rms': 1.463,
                'max': 3.234,
                'unit': 'degree'
            },
            'dihedrals': {
                'rms': 20.573,
                'max': 45.089,
                'unit': 'degree'
            }
        }

        for ic in ['bonds', 'angles', 'dihedrals']:
            for val in ['rms', 'max']:
                self.assertTrue(
                    abs(ic_rmsd[ic][val] - ref_ic_rmsd[ic][val]) < 1.0e-3)

    def test_write_xyz(self):

        with tempfile.TemporaryDirectory() as temp_dir:
            fname = str(Path(temp_dir, 'mol.xyz'))

            mol = self.nh3_molecule()
            mol.write_xyz(fname)

            with open(fname, 'r') as f_xyz:
                lines = f_xyz.readlines()

                labels = [line.split()[0] for line in lines[2:]]
                coords = np.array([
                    [float(x) for x in line.split()[1:]] for line in lines[2:]
                ])

                for a, b in zip(labels, mol.get_labels()):
                    self.assertEqual(a.lower(), b.lower())

                ref_coords = mol.get_coordinates() * bohr_in_angstroms()
                self.assertTrue(np.max(np.abs(coords - ref_coords)) < 1.0e-10)


if __name__ == "__main__":
    unittest.main()
