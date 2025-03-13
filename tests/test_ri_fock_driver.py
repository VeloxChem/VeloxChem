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


class TestRIFockDriver:

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

    def test_h2o_compute_bq_vector_mpi(self):

        comm = MPI.COMM_WORLD
        
        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()
        
        # compute J metric
        t2c_drv = TwoCenterElectronRepulsionDriver()
        matj = t2c_drv.compute(mol_h2o, bas_aux)
        rmatj = np.linalg.inv(matj.full_matrix().to_numpy())
        invmatj = SubMatrix([0, 0, rmatj.shape[0], rmatj.shape[0]])
        invmatj.set_values(rmatj)

        # load density matrix
        den_mat = None
        if comm.Get_rank() == 0:
            here = Path(__file__).parent
            npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
            density = np.load(npyfile)
            den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
            den_mat.set_values(density)
        den_mat = comm.bcast(den_mat, 0)
        
        gri_fock_drv = RIFockDriver(invmatj)
        gri_fock_drv.prepare_buffers(mol_h2o, bas_sto3g, bas_aux)
        mgv = gri_fock_drv.compute_bq_vector(den_mat)
        
        ri_fock_drv = RIFockDriver(invmatj)
        ri_fock_drv.mpi_prepare_buffers(comm, mol_h2o, bas_sto3g, bas_aux)
        gv = ri_fock_drv.mpi_compute_bq_vector(comm, den_mat)
        
        for rval, val in zip(mgv, gv):
            assert mt.isclose(rval, val, rel_tol=1.0e-8, abs_tol=1.0e-8)
            
    def test_h2o_compute_mpi(self):
    
        comm = MPI.COMM_WORLD
        
        mol_h2o, bas_sto3g, bas_aux = self.get_data_h2o()
        
        # compute J metric
        t2c_drv = TwoCenterElectronRepulsionDriver()
        matj = t2c_drv.compute(mol_h2o, bas_aux)
        rmatj = np.linalg.inv(matj.full_matrix().to_numpy())
        invmatj = SubMatrix([0, 0, rmatj.shape[0], rmatj.shape[0]])
        invmatj.set_values(rmatj)

        # load density matrix
        den_mat = None
        if comm.Get_rank() == 0:
            here = Path(__file__).parent
            npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
            density = np.load(npyfile)
            den_mat = make_matrix(bas_sto3g, mat_t.symmetric)
            den_mat.set_values(density)
        den_mat = comm.bcast(den_mat, 0)
        
        if comm.Get_rank() == 0:
            gri_fock_drv = RIFockDriver(invmatj)
            gri_fock_drv.prepare_buffers(mol_h2o, bas_sto3g, bas_aux)
            rfock = gri_fock_drv.compute(den_mat, 'j')
            
        ri_fock_drv = RIFockDriver(invmatj)
        ri_fock_drv.mpi_prepare_buffers(comm, mol_h2o, bas_sto3g, bas_aux)
        gv = ri_fock_drv.mpi_compute_bq_vector(comm, den_mat)
        cfock = ri_fock_drv.mpi_compute(comm, gv, den_mat, 'j')
    
        if comm.Get_rank() == 0:
            maxval = np.max(np.abs(cfock.full_matrix().to_numpy() -
                                   rfock.full_matrix().to_numpy()))
            assert maxval < 1.0e-8
