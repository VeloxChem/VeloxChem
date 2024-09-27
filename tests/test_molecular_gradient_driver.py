from pathlib import Path
import numpy as np

from mpi4py import MPI

from veloxchem import MolecularBasis
from veloxchem import Molecule
from veloxchem import MolecularGradientDriver
from veloxchem import SubMatrix
from veloxchem import Matrix
from veloxchem import Matrices
from veloxchem import make_matrix
from veloxchem import mat_t


class TestMolecularGradientDriver:

    def get_data_h2o(self):

        h2ostr = """
            O    0.000000000000        0.000000000000        0.000000000000
            H    0.000000000000        0.740848095288        0.582094932012
            H    0.000000000000       -0.740848095288        0.582094932012
        """

        mol = Molecule.read_str(h2ostr, 'angstrom')
        bas = MolecularBasis.read(mol, 'sto-3g')

        return mol, bas
        
    def get_data_h2o_dimer(self):

        h2o_dimer_str = """
            O -1.464  0.099  0.300
            H -1.956  0.624 -0.340
            H -1.797 -0.799  0.206
            O  1.369  0.146 -0.395
            H  1.894  0.486  0.335
            H  0.451  0.163 -0.083
        """
        mol = Molecule.read_str(h2o_dimer_str, 'angstrom')
        bas = MolecularBasis.read(mol, 'def2-tzvp')

        return mol, bas

    def test_h2o_mpi_comp_electronic_grad(self):
    
        comm = MPI.COMM_WORLD
            
        mol_h2o, bas_sto3g = self.get_data_h2o()
        
        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = np.load(npyfile)
        
        # load weighted density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.energy.weighted.density.npy')
        wden_mat = np.load(npyfile)
        
        grad_drv = MolecularGradientDriver()
        mol_grad = grad_drv.mpi_comp_electronic_grad(comm, mol_h2o, bas_sto3g, den_mat, wden_mat, 0.0, 0.0)
        
        tot_grad = np.zeros((3, 3))
        comm.Reduce(mol_grad, tot_grad, MPI.SUM, 0)
        
        if comm.Get_rank() == 0:
            # electronic gradient contribution
            here = Path(__file__).parent
            npyfile = str(here / 'data' / 'h2o.sto3g.electronic.gradient.npy')
            ref_grad = np.load(npyfile)
        
            assert np.allclose(ref_grad, tot_grad, 1.0e-12, 1.0e-12, False)
        
        
    def test_h2o_comp_electronic_grad(self):

        mol_h2o, bas_sto3g = self.get_data_h2o()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.density.npy')
        den_mat = np.load(npyfile)
        
        # load weighted density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.energy.weighted.density.npy')
        wden_mat = np.load(npyfile)
        
        grad_drv = MolecularGradientDriver()
        mol_grad = grad_drv.comp_electronic_grad(mol_h2o, bas_sto3g, den_mat, wden_mat, 0.0, 0.0)
        
        # electronic gradient contribution
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.sto3g.electronic.gradient.npy')
        ref_grad = np.load(npyfile)
        
        assert np.allclose(ref_grad, mol_grad, 1.0e-12, 1.0e-12, False)
        
    def test_h2o_dimer_tzvp_comp_electronic_grad(self):

        mol_h2o_dimer, bas_tzvp = self.get_data_h2o_dimer()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.tzvp.density.npy')
        den_mat = np.load(npyfile)
        
        # load weighted density matrix (need replacement)
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'h2o.dimer.tzvp.density.npy')
        wden_mat = np.load(npyfile)
        
        grad_drv = MolecularGradientDriver()
        #mol_grad = grad_drv.comp_electronic_grad(mol_h2o_dimer, bas_tzvp, den_mat, wden_mat, 0.0, 0.0)
        mol_grad = np.zeros((6, 3))
        grad_drv.comp_fock_grad(mol_grad, mol_h2o_dimer, bas_tzvp, den_mat, 0.0, 0.0, 0)
        
        assert False
        #assert np.allclose(ref_grad, mol_grad, 1.0e-12, 1.0e-12, False)

        
        
