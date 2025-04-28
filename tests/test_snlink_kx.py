import numpy as np
import time as tm

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import XCIntegrator
from veloxchem.griddriver import GridDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.dftutils import get_default_grid_level
from veloxchem.fockdriver import FockDriver
from veloxchem.veloxchemlib import T4CScreener
from veloxchem.veloxchemlib import make_matrix, mat_t

class TestSNLinkK:

    def run_sn_link_k(self, molecule, basis):

        scf_drv = ScfRestrictedDriver()
        scf_drv.ostream.mute()

        scf_drv.conv_thresh = 1.0e-8

        scf_drv.compute(molecule, basis)
        density = scf_drv.density
        density = density.broadcast(scf_drv.comm, root=mpi_master())

        grid_drv = GridDriver(scf_drv.comm)
        grid_drv.set_level(3)
        mol_grid = grid_drv.generate(molecule)

        xc_drv = XCIntegrator()
        kx_mat = xc_drv.integrate_kx_fock(molecule, basis,
                                          [density.alpha_to_numpy(0)],
                                          mol_grid, 0.68)
        
        dmat = make_matrix(basis, mat_t.symmetric)
        dmat.set_values(density.alpha_to_numpy(0))
        
        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(basis, molecule, "eri")

        # compute Fock matrix
        fock_drv = FockDriver()
        start = tm.time()
        fock_mat = fock_drv._compute_fock_omp(t4c_drv, dmat, "kx", 0.68, 0.0, 12)
        end = tm.time()
        print("Analytical exchange matrix: ", end - start, " sec.")
        
        #print("*** SEMINUMERICAL FOCK ***")
        #print(kx_mat.alpha_to_numpy())
        
        #print("*** ANALYTICAL FOCK ***")
        #print(fock_mat.full_matrix().to_numpy())
        
        print(" * Max. abs. diff. : ", np.max(np.abs(kx_mat.alpha_to_numpy() - fock_mat.full_matrix().to_numpy())), " For grid points : ", mol_grid.number_of_points())
        return kx_mat

    def test_sn_link_kx_h2o_svp(self):

        mol_str = """
            O  0.0000000000   0.0000000000  -0.0254395383
            H  0.0000000000   0.7695699584   0.5948147012
            H  0.0000000000  -0.7695699584   0.5948147012
        """
       
        molecule = Molecule.read_molecule_string(mol_str, units='angstrom')
        basis = MolecularBasis.read(molecule, 'def2-svp')
    
        self.run_sn_link_k(molecule, basis)
        
        assert False
