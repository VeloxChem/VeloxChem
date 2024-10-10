from pathlib import Path
import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.symmetryanalyzer import SymmetryAnalyzer


class TestPointGroup:

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
        new_coords = np.matmul(coords, rot_mat.T)
        return Molecule(labels, new_coords, 'au')

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
            'Ih': self.gen_mol_Ih(),
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
