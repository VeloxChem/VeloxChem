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
        
    def get_data_c2h4(self):

        c2h4str = """
            C             -0.667        -0.000        -0.000
            C              0.667         0.000        -0.000
            H             -1.227         0.930         0.000
            H             -1.227        -0.930        -0.000
            H              1.227         0.930        -0.000
            H              1.227        -0.930         0.000
        """

        mol = Molecule.read_str(c2h4str, 'angstrom')
        bas = MolecularBasis.read(mol, 'sto-3g')

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
        
    def test_c2h4_comp_electronic_grad(self):

        mol_c2h4, bas_sto3g = self.get_data_c2h4()

        # load density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'c2h4.sto3g.density.npy')
        den_mat = np.load(npyfile)
        
        # load weighted density matrix
        here = Path(__file__).parent
        npyfile = str(here / 'data' / 'c2h4.sto3g.energy.weighted.density.npy')
        wden_mat = np.load(npyfile)
        
        grad_drv = MolecularGradientDriver()
        mol_grad = grad_drv.comp_electronic_grad(mol_c2h4, bas_sto3g, den_mat, wden_mat, 0.0, 0.0)
        
        #assert np.allclose(ref_grad, mol_grad, 1.0e-12, 1.0e-12, False)

        
        
