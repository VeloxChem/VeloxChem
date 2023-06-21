from veloxchem.veloxchemlib import OverlapDriver
from veloxchem.veloxchemlib import MolecularBasis
from veloxchem.veloxchemlib import Molecule
from tester import Tester


class TestOverlapDriver:

    def get_data(self):

        of2str = """O   0.100  -0.400  -1.000
                    F   0.300   1.400  -2.100
                    F   0.200  -1.300  -2.000"""

        mol = Molecule.read_str(of2str, 'au')

        bas = MolecularBasis.read(mol, 'DEF2-SVPD', 'basis')

        return (mol, bas)

    def test_overlap_of2_svp(self):

        mol_h2o, bas_svp = self.get_data()

        ovl_drv = OverlapDriver()
        ovl_mat = ovl_drv.compute(bas_svp, mol_h2o)

        ssmat = ovl_mat.get_submatrix((0, 0))
        print(ssmat.to_numpy())

        spmat = ovl_mat.get_submatrix((0, 1))
        print(spmat.to_numpy())

        ppmat = ovl_mat.get_submatrix((1, 1))
        print(ppmat.to_numpy())
        assert False
