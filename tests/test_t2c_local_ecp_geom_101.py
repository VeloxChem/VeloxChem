from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import LocalECPGeom101Driver
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.submatrix import SubMatrix
from veloxchem.veloxchemlib import BaseCorePotential


class TestProjectedSECPGeom100Driver:

    def get_svp_data(self):

        costr = """
            Au   0.000   0.000   0.000
             H   0.200   1.800  -1.100
             H   0.100  -0.900  -1.000
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-svp')

        return mol, bas
        
    def test_projected_ecp_auh2_svp_for_au1_au1(self):

        mol_auh2, bas_svp = self.get_svp_data()
        
        lpot = BaseCorePotential([4.78982000, 2.39491000],
                                 [30.49008890, 5.17107381],
                                 [2, 2])
                                 
        ecp_drv = LocalECPGeom101Driver()
        ecp_mats = ecp_drv.compute(mol_auh2, bas_svp, lpot, 0, 0)
        
        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'auh2.def2svp.au.ecp.only.ul.geom.101.au1.au1.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 10, 25, 35, 42]
        
        # indices map
        labels = ['X_X', 'X_Y', 'X_Z', 'Y_X', 'Y_Y', 'Y_Z', 'Z_X', 'Z_Y', 'Z_Z']

        for k, label in enumerate(labels):
            fmat = ecp_mats.matrix(label)
            for i in range(0, 4):
                for j in range(0, 4):
                    # bra side
                    sbra = basdims[i]
                    ebra = basdims[i + 1]
                    # ket side
                    sket = basdims[j]
                    eket = basdims[j + 1]
                    # load computed submatrix
                    cmat = fmat.submatrix((i, j))
                    # load reference submatrix
                    rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                    rmat.set_values(
                        np.ascontiguousarray(ref_mat[k, sbra:ebra, sket:eket]))
                    # compare submatrices
                    assert cmat == rmat
            smat = fmat.full_matrix()
            fref = SubMatrix([0, 0, 42, 42])
            fref.set_values(np.ascontiguousarray(ref_mat[k]))
            assert smat == fref
            
    def test_projected_ecp_auh2_svp_for_h2_au1(self):

        mol_auh2, bas_svp = self.get_svp_data()
        
        lpot = BaseCorePotential([4.78982000, 2.39491000],
                                 [30.49008890, 5.17107381],
                                 [2, 2])
                                 
        ecp_drv = LocalECPGeom101Driver()
        ecp_mats = ecp_drv.compute(mol_auh2, bas_svp, lpot, 1, 0)
        
        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'auh2.def2svp.au.ecp.only.ul.geom.101.h2.au1.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 10, 25, 35, 42]
        
        # indices map
        labels = ['X_X', 'X_Y', 'X_Z', 'Y_X', 'Y_Y', 'Y_Z', 'Z_X', 'Z_Y', 'Z_Z']

        for k, label in enumerate(labels):
            fmat = ecp_mats.matrix(label)
            for i in range(0, 4):
                for j in range(0, 4):
                    # bra side
                    sbra = basdims[i]
                    ebra = basdims[i + 1]
                    # ket side
                    sket = basdims[j]
                    eket = basdims[j + 1]
                    # load computed submatrix
                    cmat = fmat.submatrix((i, j))
                    # load reference submatrix
                    rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                    rmat.set_values(
                        np.ascontiguousarray(ref_mat[k, sbra:ebra, sket:eket]))
                    # compare submatrices
                    assert cmat == rmat
            smat = fmat.full_matrix()
            fref = SubMatrix([0, 0, 42, 42])
            fref.set_values(np.ascontiguousarray(ref_mat[k]))
            assert smat == fref
            
    def test_projected_ecp_auh2_svp_for_h2_h3(self):

        mol_auh2, bas_svp = self.get_svp_data()
        
        lpot = BaseCorePotential([4.78982000, 2.39491000],
                                 [30.49008890, 5.17107381],
                                 [2, 2])
                                 
        ecp_drv = LocalECPGeom101Driver()
        ecp_mats = ecp_drv.compute(mol_auh2, bas_svp, lpot, 1, 2)
        
        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'auh2.def2svp.au.ecp.only.ul.geom.101.h2.h3.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 10, 25, 35, 42]
        
        # indices map
        labels = ['X_X', 'X_Y', 'X_Z', 'Y_X', 'Y_Y', 'Y_Z', 'Z_X', 'Z_Y', 'Z_Z']

        for k, label in enumerate(labels):
            fmat = ecp_mats.matrix(label)
            for i in range(0, 4):
                for j in range(0, 4):
                    # bra side
                    sbra = basdims[i]
                    ebra = basdims[i + 1]
                    # ket side
                    sket = basdims[j]
                    eket = basdims[j + 1]
                    # load computed submatrix
                    cmat = fmat.submatrix((i, j))
                    # load reference submatrix
                    rmat = SubMatrix([sbra, sket, ebra - sbra, eket - sket])
                    rmat.set_values(
                        np.ascontiguousarray(ref_mat[k, sbra:ebra, sket:eket]))
                    # compare submatrices
                    assert cmat == rmat
            smat = fmat.full_matrix()
            fref = SubMatrix([0, 0, 42, 42])
            fref.set_values(np.ascontiguousarray(ref_mat[k]))
            assert smat == fref
