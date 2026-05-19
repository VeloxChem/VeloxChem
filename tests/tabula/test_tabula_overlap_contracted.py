import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.veloxchemlib import GtoBlock, GtoPairBlock
from veloxchem.tabulalib import tabula_overlap_seed, tabula_overlap_contracted


class TestTabulaOverlapContracted:
    """Tests for step (b) of the late-contraction overlap recursion — the
    contraction of the seed ladder [0]^m over the primitive pairs."""

    def get_basis(self):

        mol = Molecule.read_str(
            "O 0.000  0.000 -1.000\n"
            "H 0.000  1.400 -2.100\n"
            "H 0.000 -1.400 -2.100", 'au')
        bas = MolecularBasis.read(mol, 'DEF2-SVP', ostream=None)
        return mol, bas

    def check_block(self, angmom, npgtos):

        mol, bas = self.get_basis()
        block = GtoBlock(bas, mol, angmom, npgtos)
        pair_block = GtoPairBlock(block, block)

        cdim = pair_block.number_of_contracted_pairs()
        nppairs = pair_block.number_of_primitive_pairs()
        assert cdim > 0 and nppairs > 0

        seed = tabula_overlap_seed(pair_block)              # (rows, cdim*nppairs)
        contracted = tabula_overlap_contracted(pair_block)  # (rows, cdim)

        rows = seed.shape[0]
        assert contracted.shape == (rows, cdim)

        # the seed row is primitive-pair-major (ijoff = pp*cdim + ij), so the
        # contraction is the sum over the primitive-pair axis
        reference = seed.reshape(rows, nppairs, cdim).sum(axis=1)
        assert np.allclose(contracted, reference, 0.0, 1.0e-12)

    def test_s_block(self):

        # s functions, 3 primitives -> L = 0
        self.check_block(0, 3)

    def test_p_block(self):

        # p functions, 1 primitive -> L = 2, exercises the m-ladder rows
        self.check_block(1, 1)
