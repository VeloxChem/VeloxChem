from pathlib import Path
import numpy as np

from veloxchem import OverlapDriver
from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import SubMatrix


class TestOverlapDriver:

    def get_data(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-QZVP', ostream=None)

        return mol, bas

    def get_bra_ket_data(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')
        bra_bas = MolecularBasis.read(mol, 'def2-svp', ostream=None)
        ket_bas = MolecularBasis.read(mol, 'def2-svpd', ostream=None)

        return mol, bra_bas, ket_bas
        
    def get_data_sto3g(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'cc-pV5Z')

        return mol, bas
        
    def test_overlap_co_sto3g(self):

        mol_co, bas_sto3g = self.get_data_sto3g()

        # compute overlap matrix
        ovl_drv = OverlapDriver()
        ovl_mat = ovl_drv.compute(mol_co, bas_sto3g)

        # check full overlap matrix
        fmat = ovl_mat.full_matrix().to_numpy()
        for i in range(11):
            for j in range(11):
                print("(", i, ",", j, ") = ", fmat[160 + 2 * i, 160 + 2 *j])
        assert False

    def test_overlap_co_qzvp(self):

        mol_co, bas_qzvp = self.get_data()

        # compute overlap matrix
        ovl_drv = OverlapDriver()
        ovl_mat = ovl_drv.compute(mol_co, bas_qzvp)

        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.qzvp.overlap.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(5)
        basdims = [0, 14, 38, 68, 96, 114]

        # check individual overlap submatrices
        for i, j in zip(indexes[0], indexes[1]):
            # bra side
            sbra = basdims[i]
            ebra = basdims[i + 1]
            # ket side
            sket = basdims[j]
            eket = basdims[j + 1]
            # load computed submatrix
            cmat = ovl_mat.submatrix((i, j))
            # load reference submatrix
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra,
                                                         sket:eket]))
            # compare submatrices
            assert cmat == rmat

        # check full overlap matrix
        fmat = ovl_mat.full_matrix()
        fref = SubMatrix([0, 0, 114, 114])
        fref.set_values(np.ascontiguousarray(ref_mat))
        assert fmat == fref

    def test_overlap_co_svp_svpd(self):

        mol_co, bas_svp, bas_svpd = self.get_bra_ket_data()

        # compute overlap matrix
        ovl_drv = OverlapDriver()
        ovl_mat = ovl_drv.compute(mol_co, bas_svp, bas_svpd)

        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.svp.svpd.overlap.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        bradims = [0, 6, 18, 28]
        ketdims = [0, 8, 23, 43]

        # check individual overlap submatrices
        for i in range(3):
            # bra side
            sbra = bradims[i]
            ebra = bradims[i + 1]
            for j in range(3):
                # ket side
                sket = ketdims[j]
                eket = ketdims[j + 1]
                # load computed submatrix
                cmat = ovl_mat.submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(
                    np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
                # compare submatrices
                assert cmat == rmat

        # check full overlap matrix
        fmat = ovl_mat.full_matrix()
        fref = SubMatrix([0, 0, 28, 43])
        fref.set_values(np.ascontiguousarray(ref_mat))
        assert fmat == fref
