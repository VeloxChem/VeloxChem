from mpi4py import MPI
import numpy as np
import unittest
import os

from veloxchem.veloxchemlib import ChemicalElement
from veloxchem.veloxchemlib import DispersionModel
from veloxchem.veloxchemlib import bohr_in_angstroms
from veloxchem.mpitask import MpiTask
from veloxchem.molecule import Molecule


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

        inpfile = os.path.join('inputs', 'dimer.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)
        outfile = inpfile.replace('.inp', '.out')

        task = MpiTask([inpfile, outfile], MPI.COMM_WORLD)
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
        self.assertEqual("", elem.get_name())
        elem.set_atom_type("BR")
        self.assertEqual("Br", elem.get_name())

        elem2 = ChemicalElement()
        elem2.set_atom_type(35)
        self.assertEqual("Br", elem2.get_name())

        self.assertEqual(elem, elem2)

    def test_dispersion_model(self):

        ref_energies = {
            'b3lyp': (-0.2421986012033e-02, 0.4584078631740e-06),
            'blyp': (-0.3504378057810e-02, 0.9461341515956e-06),
            'hf': (-0.6671313029358e-02, 0.5473122492978e-05),
        }

        ref_gradient = {}

        ref_gradient['b3lyp'] = np.array([
            [-0.1299360785248e-03, 0.2173451050590e-03, -0.3709704540840e-05],
            [0.3994969804870e-05, -0.4278600323727e-04, -0.3004853785695e-05],
            [0.2248231831000e-04, 0.4826552264307e-04, -0.4026908304668e-04],
            [0.2585427033048e-04, 0.3687322138623e-04, 0.3605346888461e-04],
            [0.3668558637179e-04, -0.1301671081015e-03, 0.5463254511935e-05],
            [-0.3229412701673e-05, 0.4922085484071e-05, 0.5884251321327e-05],
            [0.1936253825266e-04, -0.5468305617267e-04, 0.4693862097277e-05],
            [0.1839250629302e-04, -0.6545014186048e-04, 0.3934710919238e-05],
            [0.6393301863664e-05, -0.1431962520046e-04, -0.9045906361170e-05],
        ])

        ref_gradient['blyp'] = np.array([
            [-0.1920479181714e-03, 0.2729735760896e-03, -0.2511659875033e-05],
            [0.9234017559838e-05, -0.7983770503499e-04, -0.4956308627117e-05],
            [0.4283192296416e-04, 0.6348404201255e-04, -0.6717298399052e-04],
            [0.4829075493265e-04, 0.4755520001250e-04, 0.6235137755103e-04],
            [0.4290323248918e-04, -0.1543006045860e-03, 0.6640266553842e-05],
            [-0.4698918137746e-05, 0.8662903487663e-05, 0.7308871341866e-05],
            [0.2376967810450e-04, -0.6486040126786e-04, 0.5846333490862e-05],
            [0.2249172537731e-04, -0.7949176018520e-04, 0.4807049578003e-05],
            [0.7225504881460e-05, -0.1418525052826e-04, -0.1231294602293e-04],
        ])

        ref_gradient['hf'] = np.array([
            [-0.3686765721713e-03, 0.3397769090112e-03, 0.4534616906695e-05],
            [0.3281830111142e-04, -0.2014976035809e-03, -0.1052292444208e-04],
            [0.1161479269369e-03, 0.9219122684571e-04, -0.1528315355729e-03],
            [0.1280011820621e-03, 0.6531340670539e-04, 0.1479566767179e-03],
            [0.4036674113369e-04, -0.1594768392328e-03, 0.1040317973416e-04],
            [0.4515452238650e-05, -0.5312330663685e-05, -0.6266158009138e-05],
            [0.4492572018223e-05, -0.5927187807585e-04, -0.3729549690795e-05],
            [0.3911077473617e-04, -0.6105721802557e-04, 0.2269617645742e-05],
            [0.3223621934107e-05, -0.1066567298349e-04, 0.8186076710362e-05],
        ])

        inpfile = os.path.join('inputs', 'dimer.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        molecule = task.molecule

        disp = DispersionModel()
        for xc_label in ref_energies:
            disp.compute(molecule, xc_label)

            e_disp = disp.get_energy()
            e_ref = sum(ref_energies[xc_label])
            self.assertTrue(abs(e_disp - e_ref) < 1.0e-13)
            self.assertTrue(abs(e_disp - e_ref) / abs(e_ref) < 1.0e-10)

            g_disp = disp.get_gradient().to_numpy()
            g_ref = ref_gradient[xc_label]
            max_diff = np.max(np.abs(g_disp - g_ref))
            max_rel_diff = np.max(np.abs(g_disp - g_ref) / np.abs(g_ref))
            self.assertTrue(max_diff < 1.0e-13)
            self.assertTrue(max_rel_diff < 1.0e-10)


if __name__ == "__main__":
    unittest.main()
