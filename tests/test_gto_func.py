from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import GtoBlock
from veloxchem import GtoPairBlock
from veloxchem import make_gto_blocks
from veloxchem import make_gto_pair_blocks


class TestGtoFunc:

    def get_data(self):

        h2ostr = """O   0.000   0.000  -1.000
                    H   0.000   1.400  -2.100
                    H   0.000  -1.400  -2.100"""

        mol = Molecule.read_str(h2ostr, 'au')
        bas = MolecularBasis.read(mol, 'DEF2-SVP', ostream=None)

        return (mol, bas)

    def test_make_gto_blocks(self):

        mol_h2o, bas_svp = self.get_data()

        # full molecule
        b1s = GtoBlock(bas_svp, mol_h2o, 0, 1)
        b3s = GtoBlock(bas_svp, mol_h2o, 0, 3)
        b5s = GtoBlock(bas_svp, mol_h2o, 0, 5)
        b1p = GtoBlock(bas_svp, mol_h2o, 1, 1)
        b3p = GtoBlock(bas_svp, mol_h2o, 1, 3)
        b1d = GtoBlock(bas_svp, mol_h2o, 2, 1)
        a_blocks = make_gto_blocks(bas_svp, mol_h2o)
        b_blocks = [b1s, b3s, b5s, b1p, b3p, b1d]
        assert len(a_blocks) == len(b_blocks)
        for a_blk, b_blk in zip(a_blocks, b_blocks):
            assert a_blk == b_blk

        # selected atoms in molecule
        b1s = GtoBlock(bas_svp, mol_h2o, [0, 2], 0, 1)
        b3s = GtoBlock(bas_svp, mol_h2o, [0, 2], 0, 3)
        b5s = GtoBlock(bas_svp, mol_h2o, [0, 2], 0, 5)
        b1p = GtoBlock(bas_svp, mol_h2o, [0, 2], 1, 1)
        b3p = GtoBlock(bas_svp, mol_h2o, [0, 2], 1, 3)
        b1d = GtoBlock(bas_svp, mol_h2o, [0, 2], 2, 1)
        a_blocks = make_gto_blocks(bas_svp, mol_h2o, [0, 2])
        b_blocks = [b1s, b3s, b5s, b1p, b3p, b1d]
        assert len(a_blocks) == len(b_blocks)
        for a_blk, b_blk in zip(a_blocks, b_blocks):
            assert a_blk == b_blk

        # selected atoms in molecule
        b1s = GtoBlock(bas_svp, mol_h2o, [1, 2], 0, 1)
        b3s = GtoBlock(bas_svp, mol_h2o, [1, 2], 0, 3)
        b1p = GtoBlock(bas_svp, mol_h2o, [1, 2], 1, 1)
        a_blocks = make_gto_blocks(bas_svp, mol_h2o, [1, 2])
        b_blocks = [b1s, b3s, b1p]
        assert len(a_blocks) == len(b_blocks)
        for a_blk, b_blk in zip(a_blocks, b_blocks):
            assert a_blk == b_blk

    def test_make_gto_pair_blocks(self):

        mol_h2o, bas_svp = self.get_data()

        # full molecule
        b1s = GtoBlock(bas_svp, mol_h2o, 0, 1)
        b3s = GtoBlock(bas_svp, mol_h2o, 0, 3)
        b5s = GtoBlock(bas_svp, mol_h2o, 0, 5)
        b1p = GtoBlock(bas_svp, mol_h2o, 1, 1)
        b3p = GtoBlock(bas_svp, mol_h2o, 1, 3)
        b1d = GtoBlock(bas_svp, mol_h2o, 2, 1)

        # all GTO pair blocks
        gp_b1s1s = GtoPairBlock(b1s)
        gp_b1s3s = GtoPairBlock(b1s, b3s)
        gp_b1s5s = GtoPairBlock(b1s, b5s)
        gp_b1s1p = GtoPairBlock(b1s, b1p)
        gp_b1s3p = GtoPairBlock(b1s, b3p)
        gp_b1s1d = GtoPairBlock(b1s, b1d)
        gp_b3s3s = GtoPairBlock(b3s)
        gp_b3s5s = GtoPairBlock(b3s, b5s)
        gp_b3s1p = GtoPairBlock(b3s, b1p)
        gp_b3s3p = GtoPairBlock(b3s, b3p)
        gp_b3s1d = GtoPairBlock(b3s, b1d)
        gp_b5s5s = GtoPairBlock(b5s)
        gp_b5s1p = GtoPairBlock(b5s, b1p)
        gp_b5s3p = GtoPairBlock(b5s, b3p)
        gp_b5s1d = GtoPairBlock(b5s, b1d)
        gp_b1p1p = GtoPairBlock(b1p)
        gp_b1p3p = GtoPairBlock(b1p, b3p)
        gp_b1p1d = GtoPairBlock(b1p, b1d)
        gp_b3p3p = GtoPairBlock(b3p)
        gp_b3p1d = GtoPairBlock(b3p, b1d)
        gp_b1d1d = GtoPairBlock(b1d)

        # compare GTO pair blocks
        a_pair_blocks = make_gto_pair_blocks(bas_svp, mol_h2o)
        b_pair_blocks = [
            gp_b1s1s, gp_b1s3s, gp_b1s5s, gp_b1s1p, gp_b1s3p, gp_b1s1d,
            gp_b3s3s, gp_b3s5s, gp_b3s1p, gp_b3s3p, gp_b3s1d, gp_b5s5s,
            gp_b5s1p, gp_b5s3p, gp_b5s1d, gp_b1p1p, gp_b1p3p, gp_b1p1d,
            gp_b3p3p, gp_b3p1d, gp_b1d1d
        ]
        assert len(a_pair_blocks) == len(b_pair_blocks)
        for a_blk, b_blk in zip(a_pair_blocks, b_pair_blocks):
            assert a_blk == b_blk
