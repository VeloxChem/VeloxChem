from pathlib import Path
import numpy as np
import pytest
from mpi4py import MPI

from veloxchem.molecule import Molecule
from veloxchem.symmetryanalyzer import SymmetryAnalyzer
from veloxchem.symmetryoperations import ImproperRotation, Rotation
from veloxchem.errorhandler import VeloxChemError


class TestPointGroup:

    def gen_ammonia_with_ghost_nitrogen(self):

        nh3_bq_str = """
        Bq_N   -3.710    3.019   -0.037
        H      -3.702    4.942    0.059
        H      -4.704    2.415    1.497
        H      -4.780    2.569   -1.573
        """

        mol = Molecule.read_molecule_string(nh3_bq_str.strip(), units='bohr')
        mol.set_multiplicity(2)
        return mol

    def gen_ammonia_with_ghost_nitrogen_mixed_basis(self):

        nh3_bq_str = """
        N_Bq   -3.710    3.019   -0.037
        H      -3.702    4.942    0.059  def2-svpd
        H      -4.704    2.415    1.497
        H      -4.780    2.569   -1.573
        """

        mol = Molecule.read_molecule_string(nh3_bq_str.strip(), units='bohr')
        mol.set_multiplicity(2)
        return mol

    def rotation_matrix(self):

        return np.array([[0.53667456, -0.38518553, -0.75074131],
                         [0.75403351, 0.6182462, 0.22182223],
                         [0.37870025, -0.68513046, 0.62223981]])

    def rotation_matrix_2(self):

        return np.array([[0.89061711, -0.44964967, -0.06794363],
                         [0.36751046, 0.79967259, -0.47482609],
                         [0.26783805, 0.39791825, 0.87745304]])

    def rotate_molecule(self, mol, rot_mat):

        labels = mol.get_labels()
        coords = mol.get_coordinates_in_bohr()
        u_mat, _, vh_mat = np.linalg.svd(rot_mat)
        regularized_rot_mat = np.matmul(u_mat, vh_mat)
        if np.linalg.det(regularized_rot_mat) < 0.0:
            u_mat[:, -1] *= -1.0
            regularized_rot_mat = np.matmul(u_mat, vh_mat)
        new_coords = np.matmul(coords, regularized_rot_mat.T)
        return Molecule(labels, new_coords, 'au')

    def distort_molecule(self, mol, scale=1.0e-4):

        coords = mol.get_coordinates_in_bohr().copy()
        mol_seed = int(sum(mol.get_identifiers()) * coords.shape[0])
        rng = np.random.default_rng(seed=mol_seed)
        for idx in range(coords.shape[0]):
            coords[idx] += scale * rng.uniform(-0.6, 0.6, 3)
        return Molecule(mol.get_labels(), coords, 'au')

    def gen_mol_Cinfv(self):

        molxyz = """
        3

        C              1.090000         0.000000         0.000000
        N              2.250000         0.000000         0.000000
        H              0.025000         0.000000         0.000000
        """
        return Molecule.read_xyz_string(molxyz)

    def gen_mol_Dinfh(self):

        molxyz = """
        4

        C              1.057900         0.000000         0.000000
        C              2.258200         0.000000         0.000000
        H             -0.007800         0.000000         0.000000
        H              3.324000         0.000000         0.000000
        """
        return Molecule.read_xyz_string(molxyz)

    def gen_mol_Cs(self):

        molxyz = """
        3

        O              0.000000        -0.064219         1.117721
        Cl            -0.000000         0.005396        -0.554310
        H             -0.000000         0.831957         1.494046
        """
        return Molecule.read_xyz_string(molxyz)

    def gen_mol_C2h(self):

        molxyz = """
        6

        C              0.437626         0.500889        -0.000000
        C             -0.437626        -0.500889         0.000000
        Cl            -0.042779         2.144720         0.000000
        Cl             0.042779        -2.144720        -0.000000
        H              1.513008         0.379042         0.000000
        H             -1.513008        -0.379042        -0.000000
        """
        return Molecule.read_xyz_string(molxyz)

    def gen_mol_C2v(self):

        molxyz = """
        5

        C             -0.000000        -0.000000         0.817262
        Cl             0.000000         1.457662        -0.181912
        Cl             0.000000        -1.457662        -0.181912
        H             -0.893486         0.000000         1.446368
        H              0.893496        -0.000000         1.446368
        """
        return Molecule.read_xyz_string(molxyz)

    def gen_mol_C3v(self):

        molxyz = """
        4

        N              0.069986         0.000001        -0.000000
        H             -0.324140        -0.000000        -0.939686
        H             -0.324140        -0.813800         0.469849
        H             -0.324140         0.813800         0.469849
        """
        return Molecule.read_xyz_string(molxyz)

    def gen_mol_D2h(self):

        molxyz = """
        6

        C              0.000000        -0.000000        -0.667880
        C              0.000000        -0.000000         0.667880
        H             -0.000000         0.929999        -1.227710
        H             -0.000000        -0.929999        -1.227710
        H             -0.000000         0.929999         1.227710
        H             -0.000000        -0.929999         1.227710
        """
        return Molecule.read_xyz_string(molxyz)

    def gen_mol_D2d(self):

        molxyz = """
        7

        C             -1.297220        -0.000000         0.000000
        C             -0.000000        -0.000000         0.000000
        C              1.297220        -0.000000         0.000000
        H             -1.857700        -0.000000         0.930000
        H             -1.857700        -0.000000        -0.930000
        H              1.857700         0.930000        -0.000000
        H              1.857700        -0.930000        -0.000000
        """
        return Molecule.read_xyz_string(molxyz)

    def gen_mol_D3d(self):

        molxyz = """
        8

        C              0.000000         0.000000        -0.756010
        C              0.000000         0.000000         0.756010
        H             -0.887000        -0.512110        -1.140350
        H              0.887000        -0.512110        -1.140350
        H             -0.000000         1.024220        -1.140350
        H              0.000000        -1.024220         1.140350
        H             -0.887000         0.512110         1.140350
        H              0.887000         0.512110         1.140350
        """
        return Molecule.read_xyz_string(molxyz)

    def gen_mol_D6h(self):

        molxyz = """
        12

        C              0.215190        -1.378099         0.000000
        C             -0.215190         1.378099         0.000000
        C              1.085890         0.875416         0.000000
        C             -1.085890        -0.875416         0.000000
        C              1.301060        -0.502702        -0.000000
        C             -1.301060         0.502702        -0.000000
        H              0.382840        -2.451814         0.000000
        H             -0.382840         2.451814         0.000000
        H              1.931900         1.557454         0.000000
        H             -1.931900        -1.557454         0.000000
        H              2.314750        -0.894350        -0.000000
        H             -2.314750         0.894350        -0.000000
        """
        return Molecule.read_xyz_string(molxyz)

    def gen_mol_Td(self):

        molxyz = """
        5

        C              0.000000         0.000000         0.000000
        H              0.771000         0.771000         0.771000
        H             -0.771000        -0.771000         0.771000
        H              0.771000        -0.771000        -0.771000
        H             -0.771000         0.771000        -0.771000
        """
        return Molecule.read_xyz_string(molxyz)

    def gen_mol_Oh(self):

        molxyz = """
        7

        S              0.000000         0.000000         0.000000
        F              1.564000         0.000000         0.000000
        F             -1.564000         0.000000         0.000000
        F              0.000000         1.564000         0.000000
        F              0.000000        -1.564000         0.000000
        F              0.000000         0.000000         1.564000
        F              0.000000         0.000000        -1.564000
        """
        return Molecule.read_xyz_string(molxyz)

    def gen_mol_Ih(self):

        here = Path(__file__).parent
        xyzfile = here / 'data' / 'c60_ih.xyz'
        return Molecule.read_xyz_file(xyzfile)

    def gen_mol_C1(self):

        molxyz = """
        5

        C              0.000000         0.000000         0.000000
        H              0.629000         0.629000         0.629000
        F             -0.691000        -0.691000         0.691000
        Cl             0.890000        -0.890000        -0.890000
        Br            -1.050000         1.050000        -1.050000
        """
        return Molecule.read_xyz_string(molxyz)

    def gen_mol_atom(self):

        molxyz = """
        1

        He             0.000000         0.000000         0.000000
        """
        return Molecule.read_xyz_string(molxyz)

    def test_pointgroup(self):

        sym_analyzer = SymmetryAnalyzer()

        pg_mol = {
            'Cinfv': self.gen_mol_Cinfv(),
            'Dinfh': self.gen_mol_Dinfh(),
            'Cs': self.gen_mol_Cs(),
            'C2h': self.gen_mol_C2h(),
            'C2v': self.gen_mol_C2v(),
            'C3v': self.gen_mol_C3v(),
            'D2h': self.gen_mol_D2h(),
            'D2d': self.gen_mol_D2d(),
            'D3d': self.gen_mol_D3d(),
            'D6h': self.gen_mol_D6h(),
            'Td': self.gen_mol_Td(),
            'Oh': self.gen_mol_Oh(),
        }

        for pg, mol in pg_mol.items():
            # original molecule
            sym_res = sym_analyzer.identify_pointgroup(mol)
            assert sym_res['point_group'] == pg

            # rotated molecule
            rot_mat = self.rotation_matrix()
            mol = self.rotate_molecule(mol, rot_mat)
            sym_res = sym_analyzer.identify_pointgroup(mol)
            assert sym_res['point_group'] == pg

            # further rotated molecule
            rot_mat_2 = self.rotation_matrix_2()
            mol = self.rotate_molecule(mol, rot_mat_2)
            sym_res = sym_analyzer.identify_pointgroup(mol)
            assert sym_res['point_group'] == pg

    def test_symmetrize_pointgroup(self):

        pg_mol = {
            'Cinfv': self.gen_mol_Cinfv(),
            'Dinfh': self.gen_mol_Dinfh(),
            'Cs': self.gen_mol_Cs(),
            'C2h': self.gen_mol_C2h(),
            'C2v': self.gen_mol_C2v(),
            'C3v': self.gen_mol_C3v(),
            'D2h': self.gen_mol_D2h(),
            'D2d': self.gen_mol_D2d(),
            'D3d': self.gen_mol_D3d(),
            'D6h': self.gen_mol_D6h(),
            'Td': self.gen_mol_Td(),
            'Oh': self.gen_mol_Oh(),
        }

        for pg, mol in pg_mol.items():
            for rot_mat in [
                    None,
                    self.rotation_matrix(),
                    self.rotation_matrix_2()
            ]:
                if rot_mat is not None:
                    mol = self.rotate_molecule(mol, rot_mat)

                mol = self.distort_molecule(mol)

                sym_analyzer = SymmetryAnalyzer()
                sym_res = sym_analyzer.identify_pointgroup(mol)
                assert sym_res['point_group'] == pg

                sym_mol = sym_analyzer.symmetrize_pointgroup(sym_res)

                assert sym_mol.number_of_atoms() == mol.number_of_atoms()
                assert sym_mol.get_labels() == mol.get_labels()
                np.testing.assert_allclose(sym_mol.center_of_mass_in_bohr(),
                                           np.zeros(3),
                                           atol=1.0e-12)

                sym_res_after = SymmetryAnalyzer().identify_pointgroup(
                    sym_mol, tolerance='very tight')
                assert sym_res_after['point_group'] == pg

    @pytest.mark.timeconsuming
    def test_symmetrize_pointgroup_ih_from_distorted_geometry(self):

        mol = self.gen_mol_Ih()
        distorted_mol = self.distort_molecule(mol, 5.0e-5)

        sym_analyzer = SymmetryAnalyzer()
        sym_res = sym_analyzer.identify_pointgroup(distorted_mol)
        assert sym_res['point_group'] == 'Ih'

        sym_mol = sym_analyzer.symmetrize_pointgroup(sym_res)

        assert sym_mol.number_of_atoms() == distorted_mol.number_of_atoms()
        assert sym_mol.get_labels() == distorted_mol.get_labels()
        np.testing.assert_allclose(sym_mol.center_of_mass_in_bohr(),
                                   np.zeros(3),
                                   atol=1.0e-12)

        sym_res_after = SymmetryAnalyzer().identify_pointgroup(
            sym_mol, tolerance='very tight')
        assert sym_res_after['point_group'] == 'Ih'

    def test_symmetrize_pointgroup_d6h_from_loose_geometry(self):

        mol = self.gen_mol_D6h()
        distorted_mol = self.distort_molecule(mol, scale=1.0e-3)

        distorted_res = SymmetryAnalyzer().identify_pointgroup(distorted_mol)
        assert distorted_res['point_group'] != 'D6h'

        sym_analyzer = SymmetryAnalyzer()
        sym_res = sym_analyzer.identify_pointgroup(distorted_mol,
                                                   tolerance='loose')
        assert sym_res['point_group'] == 'D6h'

        sym_mol = sym_analyzer.symmetrize_pointgroup(sym_res)

        assert sym_mol.number_of_atoms() == distorted_mol.number_of_atoms()
        assert sym_mol.get_labels() == distorted_mol.get_labels()
        np.testing.assert_allclose(sym_mol.center_of_mass_in_bohr(),
                                   np.zeros(3),
                                   atol=1.0e-12)

        sym_res_after = SymmetryAnalyzer().identify_pointgroup(
            sym_mol, tolerance='very tight')
        assert sym_res_after['point_group'] == 'D6h'

    def test_symmetrize_pointgroup_preserves_ghost_atom_info(self):

        mol = self.gen_ammonia_with_ghost_nitrogen()

        sym_analyzer = SymmetryAnalyzer()
        sym_res = sym_analyzer.identify_pointgroup(mol)
        assert sym_res['point_group'] == 'C3v'

        sym_mol = sym_analyzer.symmetrize_pointgroup(sym_res)

        assert sym_mol.number_of_atoms() == mol.number_of_atoms()
        assert sym_mol.get_labels() == mol.get_labels()
        assert sym_mol.get_identifiers() == mol.get_identifiers()
        assert sym_mol.get_atom_basis_labels() == mol.get_atom_basis_labels()
        assert 'Bq_N' in sym_mol.get_xyz_string()
        assert sym_mol.get_xyz_string().count('\nH') == 3
        np.testing.assert_allclose(sym_mol.center_of_mass_in_bohr(),
                                   np.zeros(3),
                                   atol=1.0e-12)

    def test_symmetrize_pointgroup_preserves_mixed_basis_labels(self):

        mol = self.gen_ammonia_with_ghost_nitrogen_mixed_basis()

        sym_analyzer = SymmetryAnalyzer()
        sym_res = sym_analyzer.identify_pointgroup(mol)
        assert sym_res['point_group'] == 'C3v'

        sym_mol = sym_analyzer.symmetrize_pointgroup(sym_res)

        assert sym_mol.number_of_atoms() == mol.number_of_atoms()
        assert sym_mol.get_labels() == mol.get_labels()
        assert sym_mol.get_identifiers() == mol.get_identifiers()
        assert sym_mol.get_atom_basis_labels() == [
            ('', 'N'),
            ('DEF2-SVPD', 'H'),
            ('', 'H'),
            ('', 'H'),
        ]
        assert sym_mol.get_atom_basis_labels() == mol.get_atom_basis_labels()
        assert 'Bq_N' in sym_mol.get_xyz_string()
        assert sym_mol.get_xyz_string().count('\nH') == 3
        np.testing.assert_allclose(sym_mol.center_of_mass_in_bohr(),
                                   np.zeros(3),
                                   atol=1.0e-12)

    def test_c7h_expected_symmetry_elements(self):

        sym_analyzer = SymmetryAnalyzer()

        assert sym_analyzer._all_symmetry_elements['C7h'] == [
            'E',
            '6C7',
            '6S7',
            'sigma_h',
        ]

    def test_isolated_atom_analysis_and_printing(self, capsys):

        sym_analyzer = SymmetryAnalyzer()
        results = sym_analyzer.identify_pointgroup(self.gen_mol_atom(),
                                                   tolerance='VERY TIGHT')

        assert results == {
            'degeneracy': 'Input structure is an isolated atom.',
            'point_group': 'O(3)',
        }
        assert sym_analyzer._molecule_type == 'isolated_atom'

        sym_analyzer.print_symmetry_results(results)
        sym_analyzer.print_tolerance_keywords()

        out = capsys.readouterr().out
        assert 'Input structure is an isolated atom.' in out
        assert 'Point group: O(3)' in out
        assert 'very tight' in out

        symmetrized_mol = sym_analyzer.symmetrize_pointgroup(results)
        np.testing.assert_allclose(symmetrized_mol.center_of_mass_in_bohr(),
                                   np.zeros(3))
        np.testing.assert_allclose(symmetrized_mol.get_coordinates_in_bohr(),
                                   np.zeros((1, 3)))

    @pytest.mark.skipif(MPI.COMM_WORLD.Get_size() > 1,
                        reason='skip pytest.raises for multiple MPI processes')
    def test_invalid_tolerance_and_identity_operation_raise(self):

        sym_analyzer = SymmetryAnalyzer()

        with pytest.raises(VeloxChemError, match='Invalid tolerance keyword'):
            sym_analyzer.identify_pointgroup(self.gen_mol_C2v(),
                                             tolerance='invalid')

        sym_analyzer.identify_pointgroup(self.gen_mol_C2v())
        with pytest.raises(VeloxChemError,
                           match='Identity matrix should not be checked'):
            sym_analyzer._check_symmetry_operation(Rotation([1.0, 0.0, 0.0]),
                                                   'C1')

    def test_c1_symmetrization_and_operation_mapping(self):

        sym_analyzer = SymmetryAnalyzer()
        results = sym_analyzer.identify_pointgroup(self.gen_mol_C1())

        assert results['point_group'] == 'C1'
        assert results['expected_symmetry_elements'] == ['E']
        assert results['degeneracy'] == (
            'Principal moments of inertia: Nondegenerate')

        symmetrized_mol = sym_analyzer.symmetrize_pointgroup(results)
        np.testing.assert_allclose(symmetrized_mol.center_of_mass_in_bohr(),
                                   np.zeros(3),
                                   atol=1.0e-10)
        assert symmetrized_mol.number_of_atoms() == 5

        assert sym_analyzer._check_symmetry_operation(Rotation([1.0, 0.0, 0.0],
                                                               order=2),
                                                      'C2',
                                                      mapping=True) is False
        assert sym_analyzer._mapping == set()

    def test_helper_methods_for_symmetry_elements_and_axes(self):

        sym_analyzer = SymmetryAnalyzer()

        assert sym_analyzer._get_sym_op_name('4S10') == 'S10'
        assert sym_analyzer._is_rotation('C8') is True
        assert sym_analyzer._is_rotation('sigma_h') is False
        assert sym_analyzer._is_improper_rotation('S6') is True
        assert sym_analyzer._is_improper_rotation('C6') is False

        reordered = sym_analyzer._reorder_symmetry_elements(
            ['E', 'sigma_h', 'C2', '3C2', '2S4', '2C4'])
        assert reordered == ['2S4', '2C4', 'C2', '3C2', 'E', 'sigma_h']

        symmetry_operations = [
            Rotation([1.0, 0.0, 0.0], order=2),
            Rotation([0.0, 1.0, 0.0], order=4),
            Rotation([0.0, 0.0, 1.0], order=4, power=2),
            ImproperRotation([1.0, 1.0, 0.0], order=4),
        ]
        highest_axes = sym_analyzer._get_highest_order_rotation_axes(
            symmetry_operations)
        assert len(highest_axes) == 1
        np.testing.assert_allclose(highest_axes[0], np.array([0.0, 1.0, 0.0]))

        t_axes = sym_analyzer._get_idealized_rotation_axes_T(
            [np.array([1.0, 0.0, 0.0]),
             np.array([0.0, 1.0, 0.0])])
        assert len(t_axes) == 10

        o_axes = sym_analyzer._get_idealized_rotation_axes_O([
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0])
        ])
        assert len(o_axes) == 15

        i_axes = sym_analyzer._get_idealized_rotation_axes_I([
            np.array([0.0, 0.0, 1.0]),
            np.array([0.5257311121, 0.0, 0.8506508083]),
        ])
        assert len(i_axes) == 66

        for axis in t_axes + o_axes + i_axes:
            np.testing.assert_allclose(np.linalg.norm(axis), 1.0)

    def test_helper_methods_for_orientation_grid_and_orders(self):

        sym_analyzer = SymmetryAnalyzer()
        sym_analyzer.identify_pointgroup(self.gen_mol_D6h())

        assert sym_analyzer._get_degeneracy(np.array([1.0, 2.0, 3.0]),
                                            1.0e-3) == 1
        assert sym_analyzer._get_degeneracy(np.array([1.0, 1.0001, 3.0]),
                                            1.0e-3) == 2
        assert sym_analyzer._get_degeneracy(np.array([1.0, 1.0001, 1.0002]),
                                            1.0e-3) == 3

        assert sym_analyzer._get_nondegenerate(np.array([1.0, 1.0001, 3.0]),
                                               1.0e-3) == 2
        assert sym_analyzer._get_nondegenerate(np.array([1.0, 3.0, 3.0001]),
                                               1.0e-3) == 0

        grid_points = sym_analyzer._get_mol_grid_points()
        spherical_grid_points = sym_analyzer._get_mol_grid_points('spherical')
        assert len(grid_points) >= 9
        assert len(grid_points) == len(spherical_grid_points)
        for axis in grid_points + spherical_grid_points:
            np.testing.assert_allclose(np.linalg.norm(axis), 1.0)

        custom_analyzer = SymmetryAnalyzer()
        custom_analyzer._tolerance_ang = 0.2
        custom_analyzer._tolerance_eig = 1.0e-3
        custom_analyzer._natoms = 2
        custom_analyzer._symbols = ['H', 'H']
        custom_analyzer._centered_coords = np.array([[-1.0, 0.0, 0.0],
                                                     [1.0, 0.0, 0.0]])

        assert custom_analyzer._get_axis_rot_order(np.array([0.0, 1.0, 0.0]),
                                                   8) == 2

        custom_analyzer._set_orientation(np.array([0.0, 1.0, 0.0]),
                                         np.array([0.0, 0.0, 1.0]))
        np.testing.assert_allclose(
            custom_analyzer._centered_coords,
            np.array([[0.0, 0.0, -1.0], [0.0, 0.0, 1.0]]))
