from veloxchem.veloxchemlib import MolecularBasis
from veloxchem.veloxchemlib import Molecule
from veloxchem.veloxchemlib import GtoBlock
from veloxchem.veloxchemlib import make_gto_blocks
from tester import Tester


class TestGtoBlock:

    def get_data(self):

        h2ostr = """O   0.000   0.000  -1.000
                    H   0.000   1.400  -2.100
                    H   0.000  -1.400  -2.100"""

        mol = Molecule.read_str(h2ostr, 'au')

        bas = MolecularBasis.read(mol, 'DEF2-SVP', 'basis')

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
            Tester.compare_gto_blocks(a_blk, b_blk)

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
            Tester.compare_gto_blocks(a_blk, b_blk)

        # selected atoms in molecule

        b1s = GtoBlock(bas_svp, mol_h2o, [1, 2], 0, 1)
        b3s = GtoBlock(bas_svp, mol_h2o, [1, 2], 0, 3)
        b1p = GtoBlock(bas_svp, mol_h2o, [1, 2], 1, 1)

        a_blocks = make_gto_blocks(bas_svp, mol_h2o, [1, 2])
        b_blocks = [b1s, b3s, b1p]
        assert len(a_blocks) == len(b_blocks)
        for a_blk, b_blk in zip(a_blocks, b_blocks):
            Tester.compare_gto_blocks(a_blk, b_blk)
