from pathlib import Path
import numpy as np
import math as mt

from mpi4py import MPI

from veloxchem.veloxchemlib import T4CScreener, SubMatrix, Matrices
from veloxchem.veloxchemlib import make_matrix, mat_t
from veloxchem.veloxchemlib import mpi_master
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem import TwoCenterElectronRepulsionDriver
from veloxchem.rijkfockdriver import RIJKFockDriver
<<<<<<< HEAD
from veloxchem.rifockdriver import RIFockDriver
=======
from veloxchem.scfrestdriver import ScfRestrictedDriver
>>>>>>> ef6622c93f36830652a4f89a0ac01d467a45b70e
from veloxchem.fockdriver import FockDriver


class TestRIMODriver:

    def get_data_h2o(self):

        h2ostr = """
            O    0.000000000000        0.000000000000        0.000000000000
            H    0.000000000000        0.740848095288        0.582094932012
            H    0.000000000000       -0.740848095288        0.582094932012
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'sto-3g')
        aux_bas = MolecularBasis.read(mol, 'def2-universal-jkfit')

        return mol, bas, aux_bas

    def test_h2o_compute_j_fock(self):

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()
        
        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
        den_mat.set_values(np.load(npyfile))
    
        ri_fock_drv.compute_bq_vectors(mol_h2o, bas_sto3g, bas_sto3g)
        
        
    def test_h2o_compute_k_fock(self):

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()
        
        scf_drv = ScfRestrictedDriver()
        scf_results = scf_drv.compute(mol_h2o, bas_sto3g)
        
        dmat = make_matrix(bas_sto3g, mat_t.symmetric)
        dmat.set_values(scf_results['D_alpha'])
        
        molorbs = scf_drv.molecular_orbitals
        
        fmat = ri_fock_drv.compute_j_fock(den_mat, "j")
        
        #print(fmat.full_matrix().to_numpy())
        fmat = ri_fock_drv.compute_k_fock(dmat, molorbs)
        
        assert False
            
