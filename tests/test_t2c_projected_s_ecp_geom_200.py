from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import ProjectedECPGeom200Driver
from veloxchem.veloxchemlib import BasisFunction
from veloxchem.veloxchemlib import AtomBasis
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.submatrix import SubMatrix
from veloxchem.veloxchemlib import BaseCorePotential


class TestProjectedSECPGeom200Driver:

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
        
    def test_projected_ecp_auh2_svp_for_au(self):

        mol_auh2, bas_svp = self.get_svp_data()
        
        lpot = BaseCorePotential([ 13.20510000,  6.60255000,   4.78982000,  2.39491000],
                                 [426.84667920, 37.00708285, -30.49008890, -5.17107381],
                                 [2, 2, 2, 2])
                                 
        ecp_drv = ProjectedECPGeom200Driver()
        ecp_mats = ecp_drv.compute(mol_auh2, bas_svp, lpot, 0, 0)
        
        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'auh2.def2svp.au.ecp.only.s.geom.200.au1.au1.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 10, 25, 35, 42]
        
        # indices map
        labels = ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ']
        rindex = [0, 1, 2, 4, 5, 8]

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
                        np.ascontiguousarray(ref_mat[rindex[k], sbra:ebra, sket:eket]))
                    # compare submatrices
                    print("(k,i,j) = ", k, " ", i, " ", j)
                    print(np.max(np.abs(cmat.to_numpy()-rmat.to_numpy())))
                    assert cmat == rmat
            smat = fmat.full_matrix()
            fref = SubMatrix([0, 0, 42, 42])
            fref.set_values(np.ascontiguousarray(ref_mat[rindex[k]]))
            assert smat == fref
       
    def test_projected_ecp_auh2_svp_for_h(self):

        mol_auh2, bas_svp = self.get_svp_data()
        
        lpot = BaseCorePotential([ 13.20510000,  6.60255000,   4.78982000,  2.39491000],
                                 [426.84667920, 37.00708285, -30.49008890, -5.17107381],
                                 [2, 2, 2, 2])
                                 
        ecp_drv = ProjectedECPGeom200Driver()
        ecp_mats = ecp_drv.compute(mol_auh2, bas_svp, lpot, 0, 1)
        
        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'auh2.def2svp.au.ecp.only.s.geom.200.h2.h2.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        basdims = [0, 10, 25, 35, 42]
        
        # indices map
        labels = ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ']
        rindex = [0, 1, 2, 4, 5, 8]

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
                        np.ascontiguousarray(ref_mat[rindex[k], sbra:ebra, sket:eket]))
                    # compare submatrices
                    print("(k,i,j) = ", k, " ", i, " ", j)
                    print(np.max(np.abs(cmat.to_numpy()-rmat.to_numpy())))
                    #assert cmat == rmat
            smat = fmat.full_matrix()
            fref = SubMatrix([0, 0, 42, 42])
            fref.set_values(np.ascontiguousarray(ref_mat[rindex[k]]))
            #assert smat == fref
        
        assert False
    
