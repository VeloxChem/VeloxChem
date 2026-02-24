from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import LocalECPDriver
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.submatrix import SubMatrix
from veloxchem.veloxchemlib import BaseCorePotential


class TestLocalECPDriver:

    def get_svp_data(self):

        costr = """
            Au   0.000   0.000   0.000
             H   0.200   1.800  -1.100
             H   0.100  -0.900  -1.000
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-svp')

        return mol, bas
        
    def get_tzvpp_data(self):

        costr = """
            Au   0.000   0.000   0.000
             H   0.200   1.800  -1.100
             H   0.100  -0.900  -1.000
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-tzvpp')

        return mol, bas

    def test_local_ecp_auh2_svp(self):

        mol_auh2, bas_svp = self.get_svp_data()
        
        lpot = BaseCorePotential([4.78982000, 2.39491000],
                                 [30.49008890, 5.17107381],
                                 [2, 2])
                                 
        ecp_drv = LocalECPDriver()
        ecp_mat = ecp_drv.compute(mol_auh2, bas_svp, lpot, 0)
        
        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'auh2.def2svp.au.ecp.only.ul.npy')
        ref_mat = np.load(npyfile)
        
        # dimension of molecular basis
        indexes = np.triu_indices(4)
        basdims = [0, 10, 25, 35, 42]
        
        # check individual overlap submatrices
        for i, j in zip(indexes[0], indexes[1]):
            # bra side
            sbra = basdims[i]
            ebra = basdims[i + 1]
            # ket side
            sket = basdims[j]
            eket = basdims[j + 1]
            # load computed submatrix
            cmat = ecp_mat.submatrix((i, j))
            # load reference submatrix
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra,
                                                         sket:eket]))
            # compare submatrices
            assert cmat == rmat

        # check full overlap matrix
        fmat = ecp_mat.full_matrix()
        fref = SubMatrix([0, 0, 42, 42])
        fref.set_values(np.ascontiguousarray(ref_mat))
        assert fmat == fref
        
    def test_local_ecp_auh2_tzvpp(self):

        mol_auh2, bas_tzvpp = self.get_tzvpp_data()
        
        lpot = BaseCorePotential([4.78982000, 2.39491000],
                                 [30.49008890, 5.17107381],
                                 [2, 2])
                                 
        ecp_drv = LocalECPDriver()
        ecp_mat = ecp_drv.compute(mol_auh2, bas_tzvpp, lpot, 0)
        
        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'auh2.def2tzvpp.au.ecp.only.ul.npy')
        ref_mat = np.load(npyfile)
        
        # dimension of molecular basis
        indexes = np.triu_indices(5)
        basdims = [0, 12, 36, 61, 75, 84]
        
        # check individual overlap submatrices
        for i, j in zip(indexes[0], indexes[1]):
            # bra side
            sbra = basdims[i]
            ebra = basdims[i + 1]
            # ket side
            sket = basdims[j]
            eket = basdims[j + 1]
            # load computed submatrix
            cmat = ecp_mat.submatrix((i, j))
            # load reference submatrix
            rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
            rmat.set_values(np.ascontiguousarray(ref_mat[sbra:ebra,
                                                         sket:eket]))
            # compare submatrices
            assert cmat == rmat

        # check full overlap matrix
        fmat = ecp_mat.full_matrix()
        fref = SubMatrix([0, 0, 84, 84])
        fref.set_values(np.ascontiguousarray(ref_mat))
        assert fmat == fref

