from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import ECPDriver
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.molecule import Molecule
from veloxchem.submatrix import SubMatrix
from veloxchem.veloxchemlib import BaseCorePotential
from veloxchem.veloxchemlib import AtomCorePotential


class TestECPDriver:

    def get_auh2_svp_data(self):

        auh2str = """
            Au   0.000   0.000   0.000
             H   0.200   1.800  -1.100
             H   0.100  -0.900  -1.000
        """
        mol = Molecule.read_str(auh2str, 'au')
        bas = MolecularBasis.read(mol, 'def2-svp')

        return mol, bas
        
    def get_aucl_svp_data(self):

        auclstr = """
            Au 0.000   0.000   0.000
            Cl 0.000   2.550   0.000
        """
        mol = Molecule.read_str(auclstr)
        bas = MolecularBasis.read(mol, 'def2-svp')

        return mol, bas
        
    def get_gdh3_svp_data(self):

        gdh3str = """
             Gd   0.000   0.000   0.000
              H   0.200   1.800  -1.100
              H   0.100  -0.900  -1.000
              H   0.300   1.600   2.000
        """
        mol = Molecule.read_str(gdh3str, 'au')
        bas = MolecularBasis.read(mol, 'def2-svp')

        return mol, bas
        
    def get_gdh3_tzvpp_data(self):

        gdh3str = """
             Gd   0.000   0.000   0.000
              H   0.200   1.800  -1.100
              H   0.100  -0.900  -1.000
              H   0.300   1.600   2.000
        """
        mol = Molecule.read_str(gdh3str, 'au')
        bas = MolecularBasis.read(mol, 'def2-tzvpp')

        return mol, bas
    
    def test_ecp_auh2_svp(self):

        mol_auh2, bas_svp = self.get_auh2_svp_data()
        
        lpot = BaseCorePotential([4.78982000, 2.39491000],
                                 [30.49008890, 5.17107381],
                                 [2, 2])
        
        spot = BaseCorePotential([ 13.20510000,  6.60255000,   4.78982000,  2.39491000],
                                 [426.84667920, 37.00708285, -30.49008890, -5.17107381],
                                 [2, 2, 2, 2])
                                 
        ppot = BaseCorePotential([ 10.45202000, 5.22601000, 4.78982000, 2.39491000],
                                 [261.19958038, 26.96249604, -30.49008890, -5.17107381],
                                 [2, 2, 2, 2])
                                 
        dpot = BaseCorePotential([7.85110000, 3.92555000, 4.78982000, 2.39491000],
                                 [124.79066561, 16.30072573, -30.49008890, -5.17107381],
                                 [2, 2, 2, 2])
        
        atom_pot = AtomCorePotential(lpot, [spot, ppot, dpot], [0, 1, 2], 60);
        
        ecp_drv = ECPDriver()
        ecp_mat = ecp_drv.compute(mol_auh2, bas_svp, atom_pot)
        
        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'auh2.def2svp.au.ecp.full.npy')
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
        
    def test_ecp_gdh3_svp(self):

        mol_gdh3, bas_svp = self.get_gdh3_svp_data()
        
        lpot = BaseCorePotential([1.000000000000e+00, ], [0.000000000000e+00, ], [2,])
                
        spot = BaseCorePotential([2.460215100000e+01, ], [6.372008690000e+02, ], [2,])
                                 
        ppot = BaseCorePotential([1.688925000000e+01, ], [2.616896010000e+02, ], [2,])
        
        dpot = BaseCorePotential([1.364335800000e+01, ], [1.068565330000e+02, ], [2,])
        
        fpot = BaseCorePotential([2.412691700000e+01, ], [-5.068359000000e+01, ], [2,])
        
        gpot = BaseCorePotential([2.213188700000e+01, ], [-2.757963000000e+01, ], [2,])
        
        atom_pot = AtomCorePotential(lpot, [spot, ppot, dpot, fpot, gpot], [0, 1, 2, 3, 4], 60);
        
        ecp_drv = ECPDriver()
        ecp_mat = ecp_drv.compute(mol_gdh3, bas_svp, atom_pot)
        
        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'gdh3.def2svp.gd.ecp.full.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(5)
        basdims = [0, 16, 43, 63, 84, 93]
        
        #print(ref_mat.shape)
        
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
        fref = SubMatrix([0, 0, 93, 93])
        fref.set_values(np.ascontiguousarray(ref_mat))
        assert fmat == fref
        
    def test_ecp_gdh3_tzvpp(self):

        mol_gdh3, bas_tzvpp = self.get_gdh3_tzvpp_data()
        
        lpot = BaseCorePotential([1.000000000000e+00, ], [0.000000000000e+00, ], [2,])
                
        spot = BaseCorePotential([2.460215100000e+01, ], [6.372008690000e+02, ], [2,])
                                 
        ppot = BaseCorePotential([1.688925000000e+01, ], [2.616896010000e+02, ], [2,])
        
        dpot = BaseCorePotential([1.364335800000e+01, ], [1.068565330000e+02, ], [2,])
        
        fpot = BaseCorePotential([2.412691700000e+01, ], [-5.068359000000e+01, ], [2,])
        
        gpot = BaseCorePotential([2.213188700000e+01, ], [-2.757963000000e+01, ], [2,])
        
        atom_pot = AtomCorePotential(lpot, [spot, ppot, dpot, fpot, gpot], [0, 1, 2, 3, 4], 60);
        
        ecp_drv = ECPDriver()
        ecp_mat = ecp_drv.compute(mol_gdh3, bas_tzvpp, atom_pot)
        
        # load reference overlap data
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'gdh3.def2tzvpp.gd.ecp.full.npy')
        ref_mat = np.load(npyfile)

        # dimension of molecular basis
        indexes = np.triu_indices(5)
        basdims = [0, 19, 58, 98, 126, 144]
        
        #print(ref_mat.shape)
        
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
        fref = SubMatrix([0, 0, 144, 144])
        fref.set_values(np.ascontiguousarray(ref_mat))
        assert fmat == fref

    def test_ecp_aucl_svp(self):

        mol_aucl, bas_svp = self.get_aucl_svp_data()
        
        lpot = BaseCorePotential([4.78982000, 2.39491000],
                                 [30.49008890, 5.17107381],
                                 [2, 2])
        
        spot = BaseCorePotential([ 13.20510000,  6.60255000,   4.78982000,  2.39491000],
                                 [426.84667920, 37.00708285, -30.49008890, -5.17107381],
                                 [2, 2, 2, 2])
                                 
        ppot = BaseCorePotential([ 10.45202000, 5.22601000, 4.78982000, 2.39491000],
                                 [261.19958038, 26.96249604, -30.49008890, -5.17107381],
                                 [2, 2, 2, 2])
                                 
        dpot = BaseCorePotential([7.85110000, 3.92555000, 4.78982000, 2.39491000],
                                 [124.79066561, 16.30072573, -30.49008890, -5.17107381],
                                 [2, 2, 2, 2])
        
        atom_pot = AtomCorePotential(lpot, [spot, ppot, dpot], [0, 1, 2], 60);
        
        #atom_pot = AtomCorePotential(lpot, [], [], 60);
        
        #print(bas_svp.info_str('TEST'))
        
        ecp_drv = ECPDriver()
        ecp_mat = ecp_drv.compute(mol_aucl, bas_svp, atom_pot)
       
        #print(np.max(ecp_mat.full_matrix().to_numpy()))
        
        #print(np.min(ecp_mat.full_matrix().to_numpy()))
       
        #assert False
