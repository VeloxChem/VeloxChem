import numpy as np

from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.veloxchemlib import GtoBlock, GtoPairBlock
from veloxchem.tabulalib import tabula_overlap_contracted, tabula_overlap_rterms


def _monomial_index(x, y, d):
    dx = d - x
    return dx * (dx + 1) // 2 + (dx - y)


def _md_reference(contracted, order, cdim, ac):
    """Reference single-centre MD recursion:
    [r]^m = AC_i*[r-1_i]^(m+1) + (r_i-1)*[r-2_i]^(m+1)."""

    level = [np.zeros(((d + 1) * (d + 2) // 2, order - d + 1, cdim))
             for d in range(order + 1)]
    level[0][0, :, :] = contracted                       # [0]^m

    for d in range(1, order + 1):
        r_index = 0
        for x in range(d, -1, -1):
            for y in range(d - x, -1, -1):
                z = d - x - y
                if x >= 1:
                    axis, r_i = 0, x
                    idx1 = _monomial_index(x - 1, y, d - 1)
                elif y >= 1:
                    axis, r_i = 1, y
                    idx1 = _monomial_index(x, y - 1, d - 1)
                else:
                    axis, r_i = 2, z
                    idx1 = _monomial_index(x, y, d - 1)
                for m in range(order - d + 1):
                    val = ac[:, axis] * level[d - 1][idx1, m + 1, :]
                    if r_i >= 2:
                        if axis == 0:
                            idx2 = _monomial_index(x - 2, y, d - 2)
                        elif axis == 1:
                            idx2 = _monomial_index(x, y - 2, d - 2)
                        else:
                            idx2 = _monomial_index(x, y, d - 2)
                        val = val + (r_i - 1) * level[d - 2][idx2, m + 1, :]
                    level[d][r_index, m, :] = val
                r_index += 1
    return level[order][:, 0, :]


class TestTabulaOverlapRterms:
    """Tests for step (c) of the overlap recursion — the single-centre MD
    recursion building [r]^0 for |r| = l_a + l_c."""

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
        l_a, l_c = pair_block.angular_momentums()
        order = l_a + l_c

        rterms = tabula_overlap_rterms(pair_block)
        assert rterms.shape == ((order + 1) * (order + 2) // 2, cdim)

        # reference — the contracted seed driven through the MD recursion
        contracted = tabula_overlap_contracted(pair_block)
        ac = np.array([
            np.array(ra.coordinates()) - np.array(rc.coordinates())
            for ra, rc in zip(pair_block.bra_coordinates(),
                              pair_block.ket_coordinates())
        ])
        reference = _md_reference(contracted, order, cdim, ac)

        assert np.allclose(rterms, reference, 0.0, 1.0e-12)

    def test_s_block(self):

        # s functions -> L = 0, [r]^0 is just [0]^0
        self.check_block(0, 3)

    def test_p_block(self):

        # p functions -> L = 2, exercises the full degree ladder
        self.check_block(1, 1)
