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
from veloxchem.rifockdriver import RIFockDriver


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

    def test_h2o_compute_bq_vector_with_k_metric(self):

        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()
        
        ri_fock_drv = RIFockDriver()
        
        ri_fock_drv.prepare_buffers(mol_h2o, bas_sto3g, bas_aux, k_metric=True, verbose=True)
        
        lambda_p = SubMatrix([0, 0, 7, 5], 0.7)
        
        lambda_h = SubMatrix([0, 0, 7, 5], 0.3)
        
        mints = ri_fock_drv.compute_mo_bq_vectors(lambda_p, lambda_h)
       
        assert True
            
