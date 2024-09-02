import numpy as np

from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import GtoBlock
from veloxchem import GtoPairBlock
from veloxchem import BlockedGtoPairBlock
from veloxchem import Point


class TestBlockedGtoPairBlock:

    def get_data(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)

        return mol, bas

    def test_is_empty_gto_pair_block(self):

        mol_co, bas_svp = self.get_data()

        # GTO pairs
        p3_gtos = GtoBlock(bas_svp, mol_co, 1, 3)
        p3x3_pairs = GtoPairBlock(p3_gtos)

        # Blocked GTO pairs
        bp3x3_pairs = BlockedGtoPairBlock(p3x3_pairs, [1.1, 0.01, 1.2])

        # check partitioning
        assert not bp3x3_pairs.is_empty_gto_pair_block(0)
        assert not bp3x3_pairs.is_empty_gto_pair_block(1)
        assert bp3x3_pairs.is_empty_gto_pair_block(2)
        assert bp3x3_pairs.is_empty_gto_pair_block(3)
        assert bp3x3_pairs.is_empty_gto_pair_block(4)
        assert bp3x3_pairs.is_empty_gto_pair_block(5)
        assert bp3x3_pairs.is_empty_gto_pair_block(6)
        assert bp3x3_pairs.is_empty_gto_pair_block(7)
        assert bp3x3_pairs.is_empty_gto_pair_block(8)
        assert bp3x3_pairs.is_empty_gto_pair_block(9)
        assert bp3x3_pairs.is_empty_gto_pair_block(10)
        assert bp3x3_pairs.is_empty_gto_pair_block(11)
        assert bp3x3_pairs.is_empty_gto_pair_block(12)
        assert bp3x3_pairs.is_empty_gto_pair_block(13)
        assert bp3x3_pairs.is_empty_gto_pair_block(14)
        assert bp3x3_pairs.is_empty_gto_pair_block(15)

    def test_get_gto_pair_block(self):

        tol = 1.0e-12

        mol_co, bas_svp = self.get_data()

        # GTO pairs
        p3_gtos = GtoBlock(bas_svp, mol_co, 1, 3)
        p3x3_pairs = GtoPairBlock(p3_gtos)

        # GTO norms and exponents
        p3norms = p3_gtos.normalization_factors()

        # GTO exponents
        p3exps = p3_gtos.exponents()

        # Blocked GTO pairs
        bp3x3_pairs = BlockedGtoPairBlock(p3x3_pairs, [1.1, 0.01, 1.2])

        # Partitioned GTO pairs
        t0_pairs = bp3x3_pairs.gto_pair_block(0)
        t1_pairs = bp3x3_pairs.gto_pair_block(1)

        # Check first GTO pair
        assert t0_pairs.number_of_contracted_pairs() == 2
        assert t0_pairs.number_of_primitive_pairs() == 9
        assert t0_pairs.angular_momentums() == (1, 1)
        assert t0_pairs.bra_atomic_indices() == [0, 1]
        assert t0_pairs.ket_atomic_indices() == [0, 1]
        assert t0_pairs.bra_orbital_indices() == [4, 6, 8]
        assert t0_pairs.ket_orbital_indices() == [4, 6, 8]
        fnorms = np.array(t0_pairs.normalization_factors())
        rnorms = np.array([
            p3norms[0] * p3norms[0], p3norms[1] * p3norms[1],
            p3norms[0] * p3norms[2], p3norms[1] * p3norms[3],
            p3norms[0] * p3norms[4], p3norms[1] * p3norms[5],
            p3norms[2] * p3norms[0], p3norms[3] * p3norms[1],
            p3norms[2] * p3norms[2], p3norms[3] * p3norms[3],
            p3norms[2] * p3norms[4], p3norms[3] * p3norms[5],
            p3norms[4] * p3norms[0], p3norms[5] * p3norms[1],
            p3norms[4] * p3norms[2], p3norms[5] * p3norms[3],
            p3norms[4] * p3norms[4], p3norms[5] * p3norms[5]
        ])
        assert np.allclose(fnorms, rnorms, tol, tol, False)
        fexps = np.array(t0_pairs.bra_exponents())
        rexps = np.array([
            p3exps[0], p3exps[1], p3exps[0], p3exps[1], p3exps[0], p3exps[1],
            p3exps[2], p3exps[3], p3exps[2], p3exps[3], p3exps[2], p3exps[3],
            p3exps[4], p3exps[5], p3exps[4], p3exps[5], p3exps[4], p3exps[5]
        ])
        assert np.allclose(fexps, rexps, tol, tol, False)
        fexps = np.array(t0_pairs.ket_exponents())
        rexps = np.array([
            p3exps[0], p3exps[1], p3exps[2], p3exps[3], p3exps[4], p3exps[5],
            p3exps[0], p3exps[1], p3exps[2], p3exps[3], p3exps[4], p3exps[5],
            p3exps[0], p3exps[1], p3exps[2], p3exps[3], p3exps[4], p3exps[5]
        ])
        assert np.allclose(fexps, rexps, tol, tol, False)
        fcoords = t0_pairs.bra_coordinates()
        assert len(fcoords) == 2
        assert fcoords[0] == Point([0.100, -0.400, -1.000])
        assert fcoords[1] == Point([0.300,  1.400, -2.100])
        fcoords = t0_pairs.ket_coordinates()
        assert len(fcoords) == 2
        assert fcoords[0] == Point([0.100, -0.400, -1.000])
        assert fcoords[1] == Point([0.300,  1.400, -2.100])
        
        # Check second GTO pair
        assert t1_pairs.number_of_contracted_pairs() == 1
        assert t1_pairs.number_of_primitive_pairs() == 9
        assert t1_pairs.angular_momentums() == (1, 1)
        assert t1_pairs.bra_atomic_indices() == [
            0,
        ]
        assert t1_pairs.ket_atomic_indices() == [
            1,
        ]
        assert t1_pairs.bra_orbital_indices() == [4, 6]
        assert t1_pairs.ket_orbital_indices() == [4, 8]
        fnorms = np.array(t1_pairs.normalization_factors())
        rnorms = np.array([
            p3norms[0] * p3norms[1], p3norms[0] * p3norms[3],
            p3norms[0] * p3norms[5], p3norms[2] * p3norms[1],
            p3norms[2] * p3norms[3], p3norms[2] * p3norms[5],
            p3norms[4] * p3norms[1], p3norms[4] * p3norms[3],
            p3norms[4] * p3norms[5]
        ])
        assert np.allclose(fnorms, rnorms, tol, tol, False)
        fexps = np.array(t1_pairs.bra_exponents())
        rexps = np.array([
            p3exps[0], p3exps[0], p3exps[0], p3exps[2], p3exps[2], p3exps[2],
            p3exps[4], p3exps[4], p3exps[4]
        ])
        assert np.allclose(fexps, rexps, tol, tol, False)
        fexps = np.array(t1_pairs.ket_exponents())
        rexps = np.array([
            p3exps[1], p3exps[3], p3exps[5], p3exps[1], p3exps[3], p3exps[5],
            p3exps[1], p3exps[3], p3exps[5]
        ])
        assert np.allclose(fexps, rexps, tol, tol, False)
        fcoords = t1_pairs.bra_coordinates()
        assert len(fcoords) == 1
        assert fcoords[0] == Point([0.100, -0.400, -1.000])
        fcoords = t1_pairs.ket_coordinates()
        assert len(fcoords) == 1
        assert fcoords[0] == Point([0.300, 1.400, -2.100])
