import numpy as np

from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import GtoBlock
from veloxchem import GtoPairBlock


class TestGtoPairBlock:

    def eq_coordinates(self, lhs, rhs):
        if len(lhs) == len(rhs):
            for lpnt, rpnt in zip(lhs, rhs):
                if lpnt != rpnt:
                    return False
            return True
        else:
            return False

    def get_data(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        return mol, bas

    def test_number_of_contracted_pairs(self):

        mol_co, bas_svp = self.get_data()

        # GTO blocks

        p3_gtos = GtoBlock(bas_svp, mol_co, 1, 3)
        p1_gtos = GtoBlock(bas_svp, mol_co, 1, 1)

        # GTO pairs

        p3x3_pairs = GtoPairBlock(p3_gtos)
        assert p3x3_pairs.number_of_contracted_pairs() == 3

        p3x1_pairs = GtoPairBlock(p3_gtos, p1_gtos)
        assert p3x1_pairs.number_of_contracted_pairs() == 4

    def test_number_of_primitive_pairs(self):

        mol_co, bas_svp = self.get_data()

        # GTO blocks

        p3_gtos = GtoBlock(bas_svp, mol_co, 1, 3)
        p1_gtos = GtoBlock(bas_svp, mol_co, 1, 1)

        # GTO pairs

        p3x3_pairs = GtoPairBlock(p3_gtos)
        assert p3x3_pairs.number_of_primitive_pairs() == 9

        p3x1_pairs = GtoPairBlock(p3_gtos, p1_gtos)
        assert p3x1_pairs.number_of_primitive_pairs() == 3

    def test_angular_momentums(self):

        mol_co, bas_svp = self.get_data()

        # GTO blocks

        s5_gtos = GtoBlock(bas_svp, mol_co, 0, 5)
        p1_gtos = GtoBlock(bas_svp, mol_co, 1, 1)

        # GTO pairs

        s5x5_pairs = GtoPairBlock(s5_gtos)
        assert s5x5_pairs.angular_momentums() == (0, 0)

        s5xp1_pairs = GtoPairBlock(s5_gtos, p1_gtos)
        assert s5xp1_pairs.angular_momentums() == (0, 1)

        p1xs5_pairs = GtoPairBlock(p1_gtos, s5_gtos)
        assert p1xs5_pairs.angular_momentums() == (1, 0)

    def test_atomic_indices(self):

        mol_co, bas_svp = self.get_data()

        # GTO blocks

        s5_gtos = GtoBlock(bas_svp, mol_co, 0, 5)
        p1_gtos = GtoBlock(bas_svp, mol_co, 1, 1)

        # GTO pairs

        s5x5_pairs = GtoPairBlock(s5_gtos)
        assert s5x5_pairs.bra_atomic_indices() == [0, 0, 1]
        assert s5x5_pairs.ket_atomic_indices() == [0, 1, 1]

        s5xp1_pairs = GtoPairBlock(s5_gtos, p1_gtos)
        assert s5xp1_pairs.bra_atomic_indices() == [0, 0, 1, 1]
        assert s5xp1_pairs.ket_atomic_indices() == [0, 1, 0, 1]

    def test_orbital_indices(self):

        mol_co, bas_svp = self.get_data()

        # GTO blocks

        s5_gtos = GtoBlock(bas_svp, mol_co, 0, 5)
        p1_gtos = GtoBlock(bas_svp, mol_co, 1, 1)

        # GTO pairs

        s5x5_pairs = GtoPairBlock(s5_gtos)
        assert s5x5_pairs.bra_orbital_indices() == [6, 0, 0, 3]
        assert s5x5_pairs.ket_orbital_indices() == [6, 0, 3, 3]

        s5xp1_pairs = GtoPairBlock(s5_gtos, p1_gtos)
        assert s5xp1_pairs.bra_orbital_indices() == [6, 0, 0, 3, 3]
        assert s5xp1_pairs.ket_orbital_indices() == [4, 7, 9, 7, 9]

    def test_normalization_factors(self):

        tol = 1.0e-12

        mol_co, bas_svp = self.get_data()

        # GTO blocks

        p3_gtos = GtoBlock(bas_svp, mol_co, 1, 3)
        p1_gtos = GtoBlock(bas_svp, mol_co, 1, 1)

        # GTO norms

        p3norms = p3_gtos.normalization_factors()
        p1norms = p1_gtos.normalization_factors()

        # GTO pairs

        p3x3_pairs = GtoPairBlock(p3_gtos)
        fnorms = np.array(p3x3_pairs.normalization_factors())
        rnorms = np.array([
            p3norms[0] * p3norms[0], p3norms[0] * p3norms[1], p3norms[1] *
            p3norms[1], p3norms[0] * p3norms[2], p3norms[0] * p3norms[3],
            p3norms[1] * p3norms[3], p3norms[0] * p3norms[4],
            p3norms[0] * p3norms[5], p3norms[1] * p3norms[5],
            p3norms[2] * p3norms[0], p3norms[2] * p3norms[1],
            p3norms[3] * p3norms[1], p3norms[2] * p3norms[2],
            p3norms[2] * p3norms[3], p3norms[3] * p3norms[3],
            p3norms[2] * p3norms[4], p3norms[2] * p3norms[5],
            p3norms[3] * p3norms[5], p3norms[4] * p3norms[0],
            p3norms[4] * p3norms[1], p3norms[5] * p3norms[1],
            p3norms[4] * p3norms[2], p3norms[4] * p3norms[3],
            p3norms[5] * p3norms[3], p3norms[4] * p3norms[4],
            p3norms[4] * p3norms[5], p3norms[5] * p3norms[5]
        ])
        assert np.allclose(fnorms, rnorms, tol, tol, False)

        p3x1_pairs = GtoPairBlock(p3_gtos, p1_gtos)
        fnorms = np.array(p3x1_pairs.normalization_factors())
        rnorms = np.array([
            p3norms[0] * p1norms[0], p3norms[0] * p1norms[1],
            p3norms[1] * p1norms[0], p3norms[1] * p1norms[1],
            p3norms[2] * p1norms[0], p3norms[2] * p1norms[1],
            p3norms[3] * p1norms[0], p3norms[3] * p1norms[1],
            p3norms[4] * p1norms[0], p3norms[4] * p1norms[1],
            p3norms[5] * p1norms[0], p3norms[5] * p1norms[1]
        ])
        assert np.allclose(fnorms, rnorms, tol, tol, False)

    def test_exponents(self):

        tol = 1.0e-12

        mol_co, bas_svp = self.get_data()

        # GTO blocks

        p3_gtos = GtoBlock(bas_svp, mol_co, 1, 3)
        p1_gtos = GtoBlock(bas_svp, mol_co, 1, 1)

        # GTO exponents

        p3exps = p3_gtos.exponents()
        p1exps = p1_gtos.exponents()

        # GTO pairs

        p3x3_pairs = GtoPairBlock(p3_gtos)

        fexps = np.array(p3x3_pairs.bra_exponents())
        rexps = np.array([
            p3exps[0], p3exps[0], p3exps[1], p3exps[0], p3exps[0], p3exps[1],
            p3exps[0], p3exps[0], p3exps[1], p3exps[2], p3exps[2], p3exps[3],
            p3exps[2], p3exps[2], p3exps[3], p3exps[2], p3exps[2], p3exps[3],
            p3exps[4], p3exps[4], p3exps[5], p3exps[4], p3exps[4], p3exps[5],
            p3exps[4], p3exps[4], p3exps[5]
        ])
        assert np.allclose(fexps, rexps, tol, tol, False)

        fexps = np.array(p3x3_pairs.ket_exponents())
        rexps = np.array([
            p3exps[0], p3exps[1], p3exps[1], p3exps[2], p3exps[3], p3exps[3],
            p3exps[4], p3exps[5], p3exps[5], p3exps[0], p3exps[1], p3exps[1],
            p3exps[2], p3exps[3], p3exps[3], p3exps[4], p3exps[5], p3exps[5],
            p3exps[0], p3exps[1], p3exps[1], p3exps[2], p3exps[3], p3exps[3],
            p3exps[4], p3exps[5], p3exps[5]
        ])
        assert np.allclose(fexps, rexps, tol, tol, False)

        p3x1_pairs = GtoPairBlock(p3_gtos, p1_gtos)

        fexps = np.array(p3x1_pairs.bra_exponents())
        rexps = np.array([
            p3exps[0], p3exps[0], p3exps[1], p3exps[1], p3exps[2], p3exps[2],
            p3exps[3], p3exps[3], p3exps[4], p3exps[4], p3exps[5], p3exps[5]
        ])
        assert np.allclose(fexps, rexps, tol, tol, False)

        fexps = np.array(p3x1_pairs.ket_exponents())
        rexps = np.array([
            p1exps[0], p1exps[1], p1exps[0], p1exps[1], p1exps[0], p1exps[1],
            p1exps[0], p1exps[1], p1exps[0], p1exps[1], p1exps[0], p1exps[1]
        ])
        assert np.allclose(fexps, rexps, tol, tol, False)

    def test_coordinates(self):

        tol = 1.0e-12

        mol_co, bas_svp = self.get_data()

        # GTO blocks

        p3_gtos = GtoBlock(bas_svp, mol_co, 1, 3)
        p1_gtos = GtoBlock(bas_svp, mol_co, 1, 1)

        # GTO coordinates

        p3coords = p3_gtos.coordinates()
        p1coords = p1_gtos.coordinates()

        # GTO pairs

        p3x3_pairs = GtoPairBlock(p3_gtos)
        fcoords = p3x3_pairs.bra_coordinates()
        rcoords = [p3coords[0], p3coords[0], p3coords[1]]
        assert self.eq_coordinates(fcoords, rcoords)

        fcoords = p3x3_pairs.ket_coordinates()
        rcoords = [p3coords[0], p3coords[1], p3coords[1]]
        assert self.eq_coordinates(fcoords, rcoords)

        p3x1_pairs = GtoPairBlock(p3_gtos, p1_gtos)
        fcoords = p3x1_pairs.bra_coordinates()
        rcoords = [p3coords[0], p3coords[0], p3coords[1], p3coords[1]]
        assert self.eq_coordinates(fcoords, rcoords)

        fcoords = p3x1_pairs.ket_coordinates()
        rcoords = [p1coords[0], p1coords[1], p1coords[0], p1coords[1]]
        assert self.eq_coordinates(fcoords, rcoords)
