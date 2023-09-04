from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import OverlapDriver
from veloxchem.veloxchemlib import MolecularBasis
from veloxchem.veloxchemlib import Molecule
from veloxchem.submatrix import SubMatrix
from tester import Tester


class TestOverlapDriver:

    def get_data(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')

        here = Path(__file__).parent
        basis_path = str(here.parent / 'basis')
        bas = MolecularBasis.read(mol, 'def2-tzvpp', basis_path)

        return mol, bas
        
    def get_bra_ket_data(self):

        costr = """
            C   0.100  -0.400  -1.000
            O   0.300   1.400  -2.100
        """
        mol = Molecule.read_str(costr, 'au')

        here = Path(__file__).parent
        basis_path = str(here.parent / 'basis')
        bra_bas = MolecularBasis.read(mol, 'def2-svp', basis_path)
        ket_bas = MolecularBasis.read(mol, 'def2-tzvp', basis_path)

        return mol, bra_bas, ket_bas

    def test_overlap_co_tzvpp(self):

        mol_co, bas_tzvpp = self.get_data()

        # compute overlap matrix
        ovl_drv = OverlapDriver()
        ovl_mat = ovl_drv.compute(mol_co, bas_tzvpp)

        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.tzvpp.overlap.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(4)
        basdims = [0, 10, 28, 48, 62]

        # check individual overlap submatrices
        for i, j in zip(indexes[0], indexes[1]):
            # bra side
            sbra = basdims[i]
            ebra = basdims[i + 1]
            # ket side
            sket = basdims[j]
            eket = basdims[j + 1]
            # load computed submatrix
            cmat = ovl_mat.get_submatrix((i, j))
            # load reference submatrix
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
            # compare submatrices
            Tester.compare_submatrices(cmat, rmat)

        # check full overlap matrix
        fmat = ovl_mat.get_full_matrix()
        fref = SubMatrix([0, 0, 62, 62])
        fref.set_values(np.ascontiguousarray(ref_mat))
        Tester.compare_submatrices(fmat, fref)
        
    def test_overlap_co_svp_tzvp(self):

        mol_co, bas_svp, bas_tzvp = self.get_bra_ket_data()

        # compute overlap matrix
        ovl_drv = OverlapDriver()
        ovl_mat = ovl_drv.compute(mol_co, bas_svp, bas_tzvp)

        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'co.svp.tzvp.overlap.npy')
        ref_mat = np.load(npyfile)
        
        # dimension of molecular basis
        bradims = [0, 6, 18, 28]
        ketdims = [0, 10, 28, 48, 62]
        
        # check individual overlap submatrices
        for i in range(3):
            # bra side
            sbra = bradims[i]
            ebra = bradims[i + 1]
            for j in range(4):
                # ket side
                sket = ketdims[j]
                eket = ketdims[j + 1]
                # load computed submatrix
                cmat = ovl_mat.get_submatrix((i, j))
                # load reference submatrix
                rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra, sket:eket]))
                # compare submatrices
                Tester.compare_submatrices(cmat, rmat)
        
        # check full overlap matrix
        fmat = ovl_mat.get_full_matrix()
        fref = SubMatrix([0, 0, 28, 62])
        fref.set_values(np.ascontiguousarray(ref_mat))
        Tester.compare_submatrices(fmat, fref)
