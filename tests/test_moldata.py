from mpi4py import MPI
from pathlib import Path
import numpy as np
import math
import pytest

from veloxchem.veloxchemlib import bohr_in_angstrom, mpi_master
from veloxchem.molecule import Molecule
from veloxchem.mpitask import MpiTask
from veloxchem.optimizationdriver import OptimizationDriver
from veloxchem.inputparser import get_random_string_serial


@pytest.mark.filterwarnings(
    'ignore:.*tostring.*tobytes:DeprecationWarning:geometric')
class TestMolData:

    def nh3_labels(self):

        return ['N', 'H', 'H', 'H']

    def nh3_coords(self):

        return [
            [-3.710, 3.019, -0.037],
            [-3.702, 4.942, 0.059],
            [-4.704, 2.415, 1.497],
            [-4.780, 2.569, -1.573],
        ]

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
        mol_2 = Molecule.read_molecule_string(xyzstr, 'au')

        array = np.array(coords)
        arrayT = np.zeros((array.shape[1], array.shape[0]))
        for i in range(array.shape[0]):
            for j in range(array.shape[1]):
                arrayT[j][i] = array[i][j]

        mol_3 = Molecule(labels, array, 'au')
        mol_4 = Molecule(labels, arrayT.T.copy(), 'au')

        array_ang = array * bohr_in_angstrom()

        mol_5 = Molecule(labels, array_ang)
        mol_6 = Molecule(labels, array_ang, 'angstrom')

        mol_7 = Molecule([7, 1, 1, 1], array, 'au')
        mol_8 = Molecule([7, 1, 1, 1], array_ang, 'angstrom')
        mol_9 = Molecule([7, 1, 1, 1], array_ang)

        assert mol_1 == mol_2
        assert mol_1 == mol_3
        assert mol_1 == mol_4
        assert mol_1 == mol_5
        assert mol_1 == mol_6
        assert mol_1 == mol_7
        assert mol_1 == mol_8
        assert mol_1 == mol_9

    def disabled_test_get_sub_molecule(self):

        # TODO: enable get_sub_molecule

        here = Path(__file__).parent
        inpfile = here / 'data' / 'dimer.inp'
        outfile = inpfile.with_suffix('.out')

        task = MpiTask([str(inpfile), str(outfile)])
        molecule = task.molecule

        mol_1 = molecule.get_sub_molecule(0, 4)
        mol_2 = molecule.get_sub_molecule(0, 9)
        mol_3 = mol_2.get_sub_molecule(0, 4)

        mol_4 = self.nh3_molecule()

        assert mol_1 == mol_3
        assert mol_1 == mol_4
        assert mol_2 == molecule

    def test_coordinates_to_numpy(self):

        mol = self.nh3_molecule()

        coords = mol.get_coordinates_in_bohr()

        x = coords[:, 0]
        y = coords[:, 1]
        z = coords[:, 2]

        x_arr = np.array([-3.710, -3.702, -4.704, -4.780])
        y_arr = np.array([3.019, 4.942, 2.415, 2.569])
        z_arr = np.array([-0.037, 0.059, 1.497, -1.573])

        assert (x == x_arr).all()
        assert (y == y_arr).all()
        assert (z == z_arr).all()

    def test_center_of_mass(self):

        mol = self.nh3_molecule()
        mol_com = mol.center_of_mass_in_bohr()
        ref_com = np.array([-3.831697, 3.070437, -0.031436])
        assert np.max(np.abs(ref_com - mol_com)) < 1.0e-6

    def test_is_linear(self):

        mol = self.nh3_molecule()
        assert (not mol.is_linear())

        mol_str = """
            H  0.0  0.0  0.0
            H  0.0  0.0  1.5
        """
        mol = Molecule.read_molecule_string(mol_str, units='au')
        assert mol.is_linear()

        mol_str = """
            C  0.0  0.0  0.0
            O  0.0  0.0  2.2
            O  0.0  0.0 -2.2
        """
        mol = Molecule.read_molecule_string(mol_str, units='au')
        assert mol.is_linear()

    def test_setters_and_getters(self):

        mol = self.nh3_molecule()

        mol.set_charge(-1)
        assert mol.get_charge() == -1.0

        mol.set_multiplicity(2)
        assert mol.get_multiplicity() == 2.0

        mol.check_multiplicity()
        mol.check_proximity(0.1)

        elem_comp = mol.get_elemental_composition()
        assert elem_comp == set((7, 1))

    def test_number_of_atoms(self):

        mol = self.nh3_molecule()

        assert mol.number_of_atoms(1) == 3
        assert mol.number_of_atoms(7) == 1

        assert mol.number_of_atoms(0, 1, 1) == 0
        assert mol.number_of_atoms(0, 2, 1) == 1

    def test_vdw_radii_and_elem_ids(self):

        # fake molecule made of H,Li,C,N,O,S,Cu,Zn,Br,Ag,Au,Hg

        mol = Molecule.read_molecule_string("""H    0.0   0.0   0.0
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

        ref_radii /= bohr_in_angstrom()

        assert (atom_radii == ref_radii).all()

        elem_ids = mol.get_identifiers()

        ref_ids = np.array([1, 3, 6, 7, 8, 16, 29, 30, 35, 47, 79, 80])

        assert (elem_ids == ref_ids).all()

    def test_distance_matrix(self):

        xyz_string = """
            3
            water
            O       0.944703950608    0.070615471110   -0.056477643546
            H       1.903732770975    0.056427544902   -0.044142825330
            H       0.642753278417   -0.616673016012    0.540640468876
        """
        mol = Molecule.read_xyz_string(xyz_string)
        distance_matrix = mol.get_distance_matrix_in_angstrom()

        ref_distrance_matrix = np.array([[0., 0.95921308, 0.95921307],
                                         [0.95921308, 0., 1.54437856],
                                         [0.95921307, 1.54437856, 0.]])

        assert np.max(np.abs(distance_matrix - ref_distrance_matrix)) < 1.0e-8

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

        init_mol = Molecule.read_molecule_string(init_xyz_str, units='angstrom')
        final_mol = Molecule.read_molecule_string(final_xyz_str,
                                                  units='angstrom')

        ic_rmsd = OptimizationDriver.get_ic_rmsd(final_mol, init_mol)
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
                assert abs(ic_rmsd[ic][val] - ref_ic_rmsd[ic][val]) < 1.0e-3

    def test_write_xyz(self):

        if MPI.COMM_WORLD.Get_rank() == mpi_master():
            here = Path(__file__).parent
            random_string = get_random_string_serial()
            fpath = here / 'data' / f'vlx_molecule_{random_string}.xyz'
            fname = str(fpath)

            mol = self.nh3_molecule()
            mol.write_xyz(fname)

            ref_labels = mol.get_labels()
            ref_coords_au = mol.get_coordinates_in_bohr()
            ref_coords_ang = mol.get_coordinates_in_angstrom()

            mol_2 = Molecule.read_xyz_file(fname)
            assert ref_labels == mol_2.get_labels()
            assert np.max(
                np.abs(ref_coords_au -
                       mol_2.get_coordinates_in_bohr())) < 1.0e-10
            assert np.max(
                np.abs(ref_coords_ang -
                       mol_2.get_coordinates_in_angstrom())) < 1.0e-10

            with open(fname, 'r') as f_xyz:
                lines = f_xyz.readlines()
                mol_3 = Molecule.read_xyz_string(''.join(lines))
                assert ref_labels == mol_3.get_labels()
                assert np.max(
                    np.abs(ref_coords_au -
                           mol_3.get_coordinates_in_bohr())) < 1.0e-10
                assert np.max(
                    np.abs(ref_coords_ang -
                           mol_3.get_coordinates_in_angstrom())) < 1.0e-10

            if fpath.is_file():
                fpath.unlink()

    def test_dihedral(self):

        xyz_string = """
            9
            xyz
            O    1.086900000000    0.113880000000   -0.060730000000
            C    2.455250000000    0.132120000000   -0.071390000000
            C    3.171673900000   -0.838788100000    0.496389800000
            C    2.491492369403   -1.966504408464    1.155862765438
            O    1.664691816845   -2.650313648401    0.565927537003
            H    0.786520000000   -0.686240000000    0.407170000000
            H    2.871553600000    0.995167700000   -0.576074500000
            H    4.254316047964   -0.842628605717    0.498660431748
            H    2.767678583706   -2.159148582998    2.205812810612
        """
        mol = Molecule.read_xyz_string(xyz_string)
        assert abs(mol.get_dihedral_in_degrees((2, 3, 4, 5)) - 55.0) < 1e-4

        mol.set_dihedral_in_degrees((2, 3, 4, 5), 270.0)
        assert abs(mol.get_dihedral((2, 3, 4, 5), 'radian') +
                   math.pi / 2.0) < 1e-4

        mol.set_dihedral((2, 3, 4, 5), math.pi / 2.0, 'radian')
        assert abs(mol.get_dihedral_in_degrees((2, 3, 4, 5)) - 90.0) < 1e-4
