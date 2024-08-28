from veloxchem.molecule import Molecule
from veloxchem.symmetryanalyzer import SymmetryAnalyzer


class TestPointGroup:

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

    def test_pointgroup(self):

        sym_analyzer = SymmetryAnalyzer()

        pg_mol = {
            'Cs': self.gen_mol_Cs(),
            'C2h': self.gen_mol_C2h(),
            'C2v': self.gen_mol_C2v(),
            'D2h': self.gen_mol_D2h(),
            'D3d': self.gen_mol_D3d(),
            'D6h': self.gen_mol_D6h(),
            'Td': self.gen_mol_Td(),
            'Oh': self.gen_mol_Oh(),
        }

        for pg, mol in pg_mol.items():
            sym_res = sym_analyzer.identify_pointgroup(mol)
            assert sym_res['Point_group'] == pg
