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
        
        ri_fock_drv = RIJKFockDriver()
        ri_fock_drv.compute_metric(mol_h2o, bas_sto3g)
        
        ri_fock_drv.compute_bq_vectors(mol_h2o, bas_sto3g, bas_sto3g)
       
        assert False
            
