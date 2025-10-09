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
from veloxchem.rifockdriver import RIFockDriver
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
        
        ri_fock_drv = RIJKFockDriver()
        ri_fock_drv.compute_metric(mol_h2o, bas_aux)
        
        ri_fock_drv.compute_bq_vectors(mol_h2o, bas_sto3g, bas_aux)
        
        fmat = ri_fock_drv.compute_j_fock(den_mat, "j")
        
        #print(fmat.full_matrix().to_numpy())
        
        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(bas_sto3g, mol_h2o, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        rmat = fock_drv._compute_fock_omp(t4c_drv, den_mat, "j", 0.0, 0.0,
                                              15)
        
        print(np.max(np.abs(rmat.full_matrix().to_numpy() - fmat.full_matrix().to_numpy())))
       
        
        ri_fock_drv = RIFockDriver()
        ri_fock_drv.prepare_buffers(mol_h2o, bas_sto3g, bas_aux)
        cfock = ri_fock_drv.compute(den_mat, 'j')
        print(np.max(np.abs(cfock.full_matrix().to_numpy() - fmat.full_matrix().to_numpy())))
       
        assert False
            
