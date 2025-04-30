import numpy as np
import time as tm
from scipy.linalg import sqrtm

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import XCIntegrator
from veloxchem.griddriver import GridDriver
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.dftutils import get_default_grid_level
from veloxchem.submatrix import SubMatrix
from veloxchem.veloxchemlib import TwoCenterElectronRepulsionDriver
from veloxchem.veloxchemlib import mat_t
from veloxchem.veloxchemlib import T4CScreener
from veloxchem.fockdriver import FockDriver
from veloxchem.veloxchemlib import make_matrix
from veloxchem.veloxchemlib import RIFockDriver


class TestRIKDriver:

    def test_ri_k(self):
    
        #xyz_string = """6
        #xyz
        #O   -0.1858140  -1.1749469   0.7662596
        #H   -0.1285513  -0.8984365   1.6808606
        #H   -0.0582782  -0.3702550   0.2638279
        #O    0.1747051  11.1050002  -0.7244430
        #H   -0.5650842  11.3134964  -1.2949455
        #H    0.9282185  11.0652990  -1.3134026
        #"""
        xyz_string = """24
            xyz
            O    -0.537797810763    2.548484593229    0.000000000000
            O     3.060053657153   -0.342736262655   -0.000000000000
            N     0.981027030126   -1.286353215374   -0.000000000000
            N    -2.292611960807    0.008027029000    0.000000000000
            N     1.262169862164    1.085372972620    0.000000000000
            N    -1.346854892845   -2.034482811562   -0.000000000000
            C    -0.915797659810    0.195858564021    0.000000000000
            C    -0.379747568515   -1.078485797996   -0.000000000000
            C    -0.121435139475    1.390724912803    0.000000000000
            C     1.846035071709   -0.196767477313   -0.000000000000
            C    -2.478662573702   -1.335323738891   -0.000000000000
            C     1.562698463388   -2.628418219958   -0.000000000000
            C    -3.317085271030    1.051202855269    0.000000000000
            C     2.214919159110    2.200764309972    0.000000000000
            H    -3.468470456205   -1.773774496739   -0.000000000000
            H     0.742628036719   -3.346033930132   -0.000000000000
            H     2.186800646177   -2.766305911953    0.886988671775
            H     2.186800646177   -2.766305911953   -0.886988671775
            H    -2.808345859284    2.015050628143    0.000000000000
            H    -3.944028242129    0.965315311095    0.892919708821
            H    -3.944028242129    0.965315311095   -0.892919708821
            H     1.635263989242    3.121966936999    0.000000000000
            H     2.852229840507    2.147373554317   -0.886283579224
            H     2.852229840507    2.147373554317    0.886283579224
        """
        molecule = Molecule.read_xyz_string(xyz_string)
        basis = MolecularBasis.read(molecule, 'def2-svp', ostream=None)
        aux_basis = MolecularBasis.read(molecule, 'def2-universal-jkfit')

        # carry out reference scf
        scf_drv = ScfRestrictedDriver()
        scf_drv.conv_thresh = 1.0e-6
        scf_results = scf_drv.compute(molecule, basis)
       
        # set up molecular orbitals
        nmos = molecule.number_of_electrons() // 2
        cmos = scf_results['C_alpha'][:,0:nmos]
        naos = cmos.shape[0]
        orbmat = SubMatrix([0, 0, naos, nmos])
        orbmat.set_values(cmos)
        
        # compute J^(-1/2) metric
        t2c_drv = TwoCenterElectronRepulsionDriver()
        matj = t2c_drv.compute(molecule, aux_basis).full_matrix().to_numpy()
        rmatj = sqrtm(np.linalg.inv(matj))
        sqmatj = SubMatrix([0, 0, rmatj.shape[0], rmatj.shape[1]])
        sqmatj.set_values(rmatj)
        
        # set up density matix
        density = scf_results['D_alpha']
        denmat = make_matrix(basis, mat_t.symmetric)
        denmat.set_values(density)
        
        # screen basis function pairs
        t4c_drv = T4CScreener()
        t4c_drv.partition(basis, molecule, "eri")
        
        # compute exchange matrix
        fock_drv = FockDriver()
        start = tm.time()
        fockmat = fock_drv._compute_fock_omp(t4c_drv, denmat, "k", 0.0, 0.0, 12)
        end = tm.time()
        print(" * Compute full exchange matrix : ", end - start, " sec.")
        
        # compute 3C integrals
        start = tm.time()
        ri_fock_drv = RIFockDriver(sqmatj)
        ri_fock_drv.prepare_buffers(molecule, basis, aux_basis)
        end = tm.time()
        print(" * Compute 3C-integrals         : ", end - start, " sec.")
        
        # compute B^Q_it matrices
        start = tm.time()
        bmats = ri_fock_drv.compute_bq_matrices(orbmat)
        end = tm.time()
        print(" * Compute B^Q_it matrices      : ", end - start, " sec.")
        
        assert False
