from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import ProjectedECPDriver
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.submatrix import SubMatrix
from veloxchem.veloxchemlib import BaseCorePotential


class TestProjectedPECPDriver:

    def get_data(self):

        costr = """
            Au   0.000   0.000   0.000
             H   0.200   1.800  -1.100
             H   0.100  -0.900  -1.000
        """
        mol = Molecule.read_str(costr, 'au')
        bas = MolecularBasis.read(mol, 'def2-svp')

        return mol, bas

    def test_projected_ecp_auh2_svp(self):

        mol_auh2, bas_svp = self.get_data()


        print(bas_svp.info_str("Basis set info"))
        
        lpot = BaseCorePotential([ 10.45202000, 5.22601000, 4.78982000, 2.39491000],
                                 [261.19958038, 26.96249604, -30.49008890, -5.17107381],
                                 [2, 2, 2, 2])
                                 
        ecp_drv = ProjectedECPDriver()
        ecp_mat = ecp_drv.compute(mol_auh2, bas_svp, lpot, 1, 0)
        
        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'auh2.def2svp.au.ecp.only.p.npy')
        ref_mat = np.load(npyfile)

        print(ref_mat.shape)
        
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
            print("XXX : ", i, " ", j)
            print("Ref. Mat. : ", rmat.to_numpy())
            print("Cal. Mat. : ", cmat.to_numpy())
            print("Max. diff : ", np.max(np.abs(rmat.to_numpy()-cmat.to_numpy())))
            assert cmat == rmat

        # check full overlap matrix
        fmat = ecp_mat.full_matrix()
        fref = SubMatrix([0, 0, 42, 42])
        fref.set_values(np.ascontiguousarray(ref_mat))
        assert fmat == fref

