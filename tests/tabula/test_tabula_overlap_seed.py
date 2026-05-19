import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.veloxchemlib import GtoBlock
from veloxchem.tabulalib import TabulaGtoPairBlock, tabula_overlap_seed


class TestTabulaOverlapSeed:
    """Tests for step (a) of the late-contraction overlap recursion — the
    seed ladder [0]^m, m = 0 .. l_a + l_c, over every primitive pair."""

    def get_basis(self):

        mol = Molecule.read_str(
            "O 0.000  0.000 -1.000\n"
            "H 0.000  1.400 -2.100\n"
            "H 0.000 -1.400 -2.100", 'au')
        bas = MolecularBasis.read(mol, 'DEF2-SVP', ostream=None)
        return mol, bas

    def reference_seed(self, pair_block):

        # [0]^0 = weight * (1/2a)^l_a * (-1/2g)^l_c, the weight folding the
        # normalization and overlap factors (and c_i*c_j) into one
        # [0]^m = (-2 rho) * [0]^(m-1),  rho = a*g/(a+g)
        a = np.array(pair_block.bra_exponents())
        g = np.array(pair_block.ket_exponents())
        w = np.array(pair_block.weights())
        l_a, l_c = pair_block.angular_momentums()

        rho = a * g / (a + g)
        s0 = w * (0.5 / a)**l_a * (-0.5 / g)**l_c
        return np.array([(-2.0 * rho)**m * s0 for m in range(l_a + l_c + 1)])

    def test_s_block_seed(self):

        mol, bas = self.get_basis()
        block = GtoBlock(bas, mol, 0, 3)          # s functions, 3 primitives
        pair_block = TabulaGtoPairBlock(block, block)   # (0, 0) -> L = 0

        seed = tabula_overlap_seed(pair_block)
        ref = self.reference_seed(pair_block)

        assert ref.shape[0] == 1                  # only [0]^0
        assert ref.shape[1] > 0
        assert seed.shape == ref.shape
        assert np.allclose(seed, ref, 0.0, 1.0e-12)

    def test_p_block_seed_ladder(self):

        mol, bas = self.get_basis()
        block = GtoBlock(bas, mol, 1, 1)          # p functions, 1 primitive
        pair_block = TabulaGtoPairBlock(block, block)   # (1, 1) -> L = 2

        seed = tabula_overlap_seed(pair_block)
        ref = self.reference_seed(pair_block)

        assert ref.shape[0] == 3                  # m = 0, 1, 2
        assert ref.shape[1] > 0
        assert seed.shape == ref.shape
        # exercises the (1/2a)^l / (-1/2g)^l factors and the m-ladder
        assert np.allclose(seed, ref, 0.0, 1.0e-12)
