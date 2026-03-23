from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import ProjectedECPGeom100Driver
from veloxchem.veloxchemlib import BasisFunction
from veloxchem.veloxchemlib import AtomBasis
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
        
    def get_tzvpp_data(self):

        costr = """
            Au   0.000   0.000   0.000
             H   0.200   1.800  -1.100
             H   0.100  -0.900  -1.000
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-tzvpp')

        return mol, bas
        
    def get_prim_data(self):
    
        costr = """
            Au   0.000   0.000   0.000
             H   0.200   1.800  -1.100
             H   0.100  -0.900  -1.000
        """
        mol = Molecule.read_str(costr, 'au')
        
        b1s = BasisFunction([1.000], [1.000], 0)
        b1s.normalize()
        
        b1p = BasisFunction([1.000], [1.000], 1)
        b1p.normalize()
        
        h_basis = AtomBasis([b1s, b1p], 'PRIMITIVE', 1)
        au_basis = AtomBasis([b1s, b1p], 'PRIMITIVE', 79)
        
        return mol, MolecularBasis([au_basis, h_basis], [0, 0, 1])
        
    def test_projected_ecp_auh2_svp_for_au(self):

        mol_auh2, bas_svp = self.get_svp_data()
        
        lpot = BaseCorePotential([ 13.20510000,  6.60255000,   4.78982000,  2.39491000],
                                 [426.84667920, 37.00708285, -30.49008890, -5.17107381],
                                 [2, 2, 2, 2])
                                 
        ecp_drv = ProjectedECPGeom100Driver()
        ecp_mats = ecp_drv.compute(mol_auh2, bas_svp, lpot, 0, 0)
        
        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'auh2.def2svp.au.ecp.only.s.geom.100.au1.npy')
        ref_mat = -np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 10, 25, 35, 42]
        
        # indices map
        labels = ['X', 'Y', 'Z']

        for k, label in enumerate(labels):
            fmat = ecp_mats.matrix(label)
            for i in range(0, 2):
                for j in range(0, 2):
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
                    print("(k,i,j) = ", k, " ", i, " ", j)
                    print(np.max(np.abs(cmat.to_numpy()-rmat.to_numpy())))
                    assert cmat == rmat
            smat = fmat.full_matrix()
            fref = SubMatrix([0, 0, 42, 42])
            fref.set_values(np.ascontiguousarray(ref_mat[k]))
            assert smat == fref
       
    def test_projected_ecp_auh2_svp_for_h2(self):

        mol_auh2, bas_svp = self.get_svp_data()
        
        lpot = BaseCorePotential([ 13.20510000,  6.60255000,   4.78982000,  2.39491000],
                                 [426.84667920, 37.00708285, -30.49008890, -5.17107381],
                                 [2, 2, 2, 2])
                                 
        ecp_drv = ProjectedECPGeom100Driver()
        ecp_mats = ecp_drv.compute(mol_auh2, bas_svp, lpot, 0, 1)
        
        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'auh2.def2svp.au.ecp.only.s.geom.100.h2.npy')
        ref_mat = -np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 10, 25, 35, 42]
        
        # indices map
        labels = ['X', 'Y', 'Z']

        for k, label in enumerate(labels):
            fmat = ecp_mats.matrix(label)
            for i in range(0, 2):
                for j in range(0, 2):
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
                    print("(k,i,j) = ", k, " ", i, " ", j)
                    print(np.max(np.abs(cmat.to_numpy()-rmat.to_numpy())))
                    #assert cmat == rmat
            smat = fmat.full_matrix()
            fref = SubMatrix([0, 0, 42, 42])
            fref.set_values(np.ascontiguousarray(ref_mat[k]))
            #assert smat == fref
            
    def test_projected_ecp_auh2_svp_for_h3(self):

        mol_auh2, bas_svp = self.get_svp_data()
        
        lpot = BaseCorePotential([ 13.20510000,  6.60255000,   4.78982000,  2.39491000],
                                 [426.84667920, 37.00708285, -30.49008890, -5.17107381],
                                 [2, 2, 2, 2])
                                 
        ecp_drv = ProjectedECPGeom100Driver()
        ecp_mats = ecp_drv.compute(mol_auh2, bas_svp, lpot, 0, 2)
        
        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'auh2.def2svp.au.ecp.only.s.geom.100.h3.npy')
        ref_mat = -np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 10, 25, 35, 42]
        
        # indices map
        labels = ['X', 'Y', 'Z']

        for k, label in enumerate(labels):
            fmat = ecp_mats.matrix(label)
            for i in range(0, 2):
                for j in range(0, 2):
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
                    print("(k,i,j) = ", k, " ", i, " ", j)
                    print(np.max(np.abs(cmat.to_numpy()-rmat.to_numpy())))
                    #assert cmat == rmat
            smat = fmat.full_matrix()
            fref = SubMatrix([0, 0, 42, 42])
            fref.set_values(np.ascontiguousarray(ref_mat[k]))
            #assert smat == fref
    
    def test_projected_ecp_auh2_primitive_for_h3(self):

        mol_auh2, bas_prim = self.get_prim_data()
        
        lpot = BaseCorePotential([ 13.20510000,],
                                 [426.84667920,],
                                 [2, ])
                                 
        ecp_drv = ProjectedECPGeom100Driver()
        ecp_mats = ecp_drv.compute(mol_auh2, bas_prim, lpot, 0, 2)
        
        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'auh2.primitive.s.p.au.ecp.only.s.geom.100.h3.npy')
        ref_mat = -np.load(npyfile)

        print(bas_prim.info_str("Basis Info"))

        # dimension of molecular basis
        basdims = [0, 3, 12]
        
        # indices map
        labels = ['X', 'Y', 'Z']

        for k, label in enumerate(labels):
            fmat = ecp_mats.matrix(label)
            for i in range(0, 2):
                for j in range(0, 2):
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
            fref = SubMatrix([0, 0, 12, 12])
            fref.set_values(np.ascontiguousarray(ref_mat[k]))
            assert smat == fref
