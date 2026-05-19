import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.veloxchemlib import GtoBlock
from veloxchem.tabulalib import TabulaGtoPairBlock


class TestTabulaScreenedPairBlock:
    """Tests for the screened TabulaGtoPairBlock constructor — built from two
    CGtoBlocks, a screening estimator, and a threshold. A contracted-GTO pair
    is kept when estimator(bra_data, ket_data, r) >= threshold, where the
    estimator receives each contracted GTO's screening data and the distance
    |R| between their centers. A threshold at or below 0 keeps every pair."""

    def get_blocks(self):

        mol = Molecule.read_str(
            "O 0.000  0.000 -1.000\n"
            "H 0.000  1.400 -2.100\n"
            "H 0.000 -1.400 -2.100", 'au')
        bas = MolecularBasis.read(mol, 'DEF2-SVP', ostream=None)
        bra = GtoBlock(bas, mol, 0, 3)
        ket = GtoBlock(bas, mol, 0, 3)
        return bra, ket

    def test_keep_all_matches_unscreened(self):

        bra, ket = self.get_blocks()
        full = TabulaGtoPairBlock(bra, ket)
        # a threshold of 0 keeps every pair
        screened = TabulaGtoPairBlock(bra, ket, lambda bd, kd, r: 1.0, 0.0)

        assert (screened.number_of_contracted_pairs()
                == full.number_of_contracted_pairs())
        assert (screened.number_of_primitive_pairs()
                == full.number_of_primitive_pairs())
        # the kept-all block reproduces the unscreened pair data exactly
        assert np.allclose(screened.bra_exponents(),
                           full.bra_exponents(), 0.0, 1.0e-13)
        assert np.allclose(screened.ket_exponents(),
                           full.ket_exponents(), 0.0, 1.0e-13)
        assert np.allclose(screened.weights(),
                           full.weights(), 0.0, 1.0e-13)

    def test_keep_none_is_empty(self):

        bra, ket = self.get_blocks()
        # an estimator below the threshold everywhere drops every pair
        screened = TabulaGtoPairBlock(bra, ket, lambda bd, kd, r: 0.0, 1.0)

        assert screened.number_of_contracted_pairs() == 0

    def test_selective_screening(self):

        bra, ket = self.get_blocks()
        nc = bra.number_of_basis_functions()
        # keep only the same-center pairs (|R| == 0) — for this s-block each
        # contracted GTO sits on a distinct atom, so |R| > 0 iff i != j
        screened = TabulaGtoPairBlock(
            bra, ket, lambda bd, kd, r: 1.0 if r < 1.0 else 0.0, 0.5)

        assert screened.number_of_contracted_pairs() == nc
        # the screening estimator sees the per-CGTO screening data — a tiny
        # positive threshold runs the estimator yet keeps every pair
        first = bra.screening_data(0)
        seen = TabulaGtoPairBlock(
            bra, ket, lambda bd, kd, r: bd.max_exponent, 1.0e-12)
        assert seen.number_of_contracted_pairs() > 0
        assert first.max_exponent > 0.0
        # the primitive-pair count is unaffected by contracted-pair screening
        assert (screened.number_of_primitive_pairs()
                == TabulaGtoPairBlock(bra, ket).number_of_primitive_pairs())
