import numpy as np

from veloxchem.veloxchemlib import BasisFunction
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.veloxchemlib import GtoBlock
from veloxchem.veloxchemlib import Point


class TestGtoBlock:

    def get_data(self):

        h2ostr = """O   0.000   0.000  -1.000
                    H   0.000   1.400  -2.100
                    H   0.000  -1.400  -2.100"""

        mol = Molecule.read_str(h2ostr, 'au')
        bas = MolecularBasis.read(mol, 'DEF2-SVP', ostream=None)

        return (mol, bas)

    def test_constructor(self):

        mol_h2o, bas_svp = self.get_data()
        a_block = GtoBlock(bas_svp, mol_h2o, 0, 3)
        b_block = GtoBlock(bas_svp, mol_h2o, [1, 2], 0, 3)
        assert a_block == b_block

    def test_coordinates(self):

        tol = 1.0e-12

        mol_h2o, bas_svp = self.get_data()
        gblock = GtoBlock(bas_svp, mol_h2o, 0, 1)
        a_coords = gblock.coordinates()
        b_coords = [
            Point([0.000, 0.000, -1.000]),
            Point([0.000, 0.000, -1.000]),
            Point([0.000, 1.400, -2.100]),
            Point([0.000, -1.400, -2.100])
        ]
        assert a_coords == b_coords
        gblock = GtoBlock(bas_svp, mol_h2o, 0, 3)
        a_coords = gblock.coordinates()
        b_coords = [
            Point([0.000, 1.400, -2.100]),
            Point([0.000, -1.400, -2.100])
        ]
        assert a_coords == b_coords
        
    def test_coordinates_with_origin(self):

        tol = 1.0e-12

        mol_h2o, bas_svp = self.get_data()
        gblock = GtoBlock(bas_svp, mol_h2o, 0, 1)
        a_coords = gblock.coordinates(Point([1.0, -0.4, 2.0]))
        b_coords = [
            Point([-1.000, 0.400, -3.000]),
            Point([-1.000, 0.400, -3.000]),
            Point([-1.000, 1.800, -4.100]),
            Point([-1.000, -1.000, -4.100])
        ]
        assert a_coords == b_coords

    def test_exponents(self):

        tol = 1.0e-12

        mol_h2o, bas_svp = self.get_data()
        gblock = GtoBlock(bas_svp, mol_h2o, 0, 3)
        a_exps = np.array(gblock.exponents())
        bf1s = BasisFunction(
            [1.301070100000e+01, 1.962257200000e+00, 4.445379600000e-01],
            [1.968215800000e-02, 1.379652400000e-01, 4.783193500000e-01], 0)
        bf1s.normalize()
        fe1s = bf1s.get_exponents()
        b_exps = np.array(
            [fe1s[0], fe1s[0], fe1s[1], fe1s[1], fe1s[2], fe1s[2]])
        assert np.allclose(a_exps, b_exps, tol, tol, False)

    def test_get_normalization_factors(self):

        tol = 1.0e-12

        mol_h2o, bas_svp = self.get_data()
        gblock = GtoBlock(bas_svp, mol_h2o, 0, 3)
        a_norms = np.array(gblock.normalization_factors())
        bf1s = BasisFunction(
            [1.301070100000e+01, 1.962257200000e+00, 4.445379600000e-01],
            [1.968215800000e-02, 1.379652400000e-01, 4.783193500000e-01], 0)
        bf1s.normalize()
        fn1s = bf1s.get_normalization_factors()
        b_norms = np.array(
            [fn1s[0], fn1s[0], fn1s[1], fn1s[1], fn1s[2], fn1s[2]])
        assert np.allclose(a_norms, b_norms, tol, tol, False)

    def test_get_orbital_indices(self):

        mol_h2o, bas_svp = self.get_data()
        gblock = GtoBlock(bas_svp, mol_h2o, 0, 3)
        a_indexes = gblock.orbital_indices()
        b_indexes = [7, 3, 5]
        assert a_indexes == b_indexes
        gblock = GtoBlock(bas_svp, mol_h2o, 0, 1)
        a_indexes = gblock.orbital_indices()
        b_indexes = [7, 1, 2, 4, 6]
        assert a_indexes == b_indexes

    def test_get_atomic_indices(self):

        mol_h2o, bas_svp = self.get_data()
        gblock = GtoBlock(bas_svp, mol_h2o, 0, 3)
        a_indexes = gblock.atomic_indices()
        b_indexes = [1, 2]
        assert a_indexes == b_indexes
        gblock = GtoBlock(bas_svp, mol_h2o, 0, 1)
        a_indexes = gblock.atomic_indices()
        b_indexes = [0, 0, 1, 2]
        assert a_indexes == b_indexes

    def test_get_angular_momentums(self):

        mol_h2o, bas_svp = self.get_data()
        gblock = GtoBlock(bas_svp, mol_h2o, 0, 3)
        assert gblock.angular_momentum() == 0

    def test_number_of_primitives(self):

        mol_h2o, bas_svp = self.get_data()
        gblock = GtoBlock(bas_svp, mol_h2o, 0, 3)
        assert gblock.number_of_primitives() == 3

    def test_number_of_basis_functions(self):

        mol_h2o, bas_svp = self.get_data()
        gblock = GtoBlock(bas_svp, mol_h2o, 0, 3)
        assert gblock.number_of_basis_functions() == 2
        gblock = GtoBlock(bas_svp, mol_h2o, 0, 1)
        assert gblock.number_of_basis_functions() == 4
