import numpy as np

from veloxchem.veloxchemlib import GradientScreeningData
from veloxchem.veloxchemlib import DenseMatrix
from veloxchem.veloxchemlib import AODensityMatrix, denmat
from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import compute_fock_gradient_gpu
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestFockGradient:

    def run_fock_grad(self, mol, bas, min_bas, coulomb_coef, exchange_coef,
                      ref_grad, tol):

        scf_drv = ScfRestrictedDriver()
        scf_drv.filename = None
        scf_drv.checkpoint_file = None
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas, min_bas)

        if scf_drv.rank == 0:
            Da = scf_results['D_alpha']
            nocc = mol.number_of_alpha_electrons()
            ene_occ = scf_results['E_alpha'][:nocc]
            mo_occ = scf_results['C_alpha'][:, :nocc].copy()
            W = np.linalg.multi_dot([mo_occ, np.diag(ene_occ), mo_occ.T])
        else:
            Da = None
            W = None
        Da = scf_drv.comm.bcast(Da, root=mpi_master())
        W = scf_drv.comm.bcast(W, root=mpi_master())

        dmat = AODensityMatrix([Da], denmat.rest)
        wmat = DenseMatrix(W)

        omega = 0.0

        num_gpus = scf_drv._get_num_gpus_per_node()

        rank = scf_drv.comm.Get_rank()
        nnodes = scf_drv.comm.Get_size()

        grad_screener = GradientScreeningData(mol, bas, dmat, wmat, num_gpus,
                                              1e-10, 1e-10, rank, nnodes)

        fock_grad = compute_fock_gradient_gpu(mol, bas, dmat, coulomb_coef,
                                              exchange_coef, omega, 'symm',
                                              1e-12, 5e-6, grad_screener,
                                              rank, nnodes)

        fock_grad = fock_grad.to_numpy()

        fock_grad = scf_drv.comm.reduce(fock_grad, root=mpi_master())

        if scf_drv.rank == 0:
            assert np.max(np.abs(fock_grad - ref_grad)) < tol

    def test_fock_grad_h2_jk(self):

        molstr = """
        H   0.7   0.0   0.0
        H  -0.7   0.0   0.0
        H   3.0   0.7   0.0
        H   3.0  -0.7   0.0
        H   0.0   4.0   0.7
        H   0.0   4.0  -0.7
        """
        mol = Molecule.read_molecule_string(molstr, units='bohr')

        basis_label = 'sto-3g'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)
        min_bas = MolecularBasis.read(mol, 'ao-start-guess', ostream=None)

        coulomb_coef = 2.0
        exchange_coef = 1.0

        ref_grad = np.array([
            [-0.88928406, 0.07261110, 0.00000000],
            [0.88844959, 0.08858140, 0.00000000],
            [-0.04663407, -0.76533970, 0.00000000],
            [-0.03951433, 0.88338200, 0.00000000],
            [0.04349144, -0.13961740, -0.77804365],
            [0.04349144, -0.13961740, 0.77804365],
        ])

        self.run_fock_grad(mol, bas, min_bas, coulomb_coef, exchange_coef,
                           ref_grad, 1e-6)

    def test_fock_grad_h2o_jk_sto3g(self):

        xyz_string = """12
        xyz
        O          -2.377339964731        1.679515893320       -0.510254883359
        H          -1.500145191104        1.314112718648       -0.241275722341
        H          -2.488440809368        2.443900504459        0.105495197187
        O           0.020898138147        0.701676611620        0.214138162174
        H           0.628544514094        0.555958451144       -0.549321702475
        H          -0.029634629308       -0.206974125261        0.594937691184
        O          -2.493912410575        4.261162052941        0.144933582271
        H          -1.632235820282        4.680959145241        0.378770123493
        H          -2.436571839272        4.235776227388       -0.848712904704
        O          -2.268528607902        3.452276473874       -2.300196783642
        H          -3.135165145954        3.418043400412       -2.767789746482
        H          -2.279698230112        2.594412630314       -1.784163013358
        """
        mol = Molecule.read_xyz_string(xyz_string)

        basis_label = '6-31g'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)
        min_bas = MolecularBasis.read(mol, 'ao-start-guess', ostream=None)

        coulomb_coef = 2.0
        exchange_coef = 1.0

        ref_grad = np.array([
            [6.02377219, 6.32582750, 2.42293662],
            [-3.04643746, 1.29247803, -1.10946098],
            [0.15041522, -2.50089340, -2.23469489],
            [-2.39152368, -0.67029145, -3.07487562],
            [-2.24808682, 0.83847966, 2.45494051],
            [-0.11735348, 3.19842575, -1.17354842],
            [4.70301782, -4.22696990, -6.55605298],
            [-2.96525293, -1.74779320, -0.76844255],
            [-0.38791367, -0.14948745, 2.89475431],
            [-3.19759647, -4.93520650, 6.71441447],
            [3.07719425, 0.20576762, 1.79325165],
            [0.39976503, 2.36966335, -1.36322213],
        ])

        self.run_fock_grad(mol, bas, min_bas, coulomb_coef, exchange_coef,
                           ref_grad, 1e-5)

    def test_fock_grad_h2o_jk_svp(self):

        xyz_string = """12
        xyz
        O          -2.377339964731        1.679515893320       -0.510254883359
        H          -1.500145191104        1.314112718648       -0.241275722341
        H          -2.488440809368        2.443900504459        0.105495197187
        O           0.020898138147        0.701676611620        0.214138162174
        H           0.628544514094        0.555958451144       -0.549321702475
        H          -0.029634629308       -0.206974125261        0.594937691184
        O          -2.493912410575        4.261162052941        0.144933582271
        H          -1.632235820282        4.680959145241        0.378770123493
        H          -2.436571839272        4.235776227388       -0.848712904704
        O          -2.268528607902        3.452276473874       -2.300196783642
        H          -3.135165145954        3.418043400412       -2.767789746482
        H          -2.279698230112        2.594412630314       -1.784163013358
        """
        mol = Molecule.read_xyz_string(xyz_string)

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)
        min_bas = MolecularBasis.read(mol, 'ao-start-guess', ostream=None)

        coulomb_coef = 2.0
        exchange_coef = 1.0

        ref_grad = np.array([
            [6.79926789, 6.76971270, 3.03231206],
            [-4.17111378, 1.97709047, -1.26659463],
            [0.54693248, -3.53484269, -3.05188471],
            [-1.88913515, -1.71007896, -3.44992636],
            [-3.01611535, 0.87721045, 3.45817872],
            [0.05356851, 4.30234926, -1.78731808],
            [5.61175864, -3.81931046, -7.33819155],
            [-4.06890958, -2.35831063, -1.37059444],
            [-0.28847997, -0.19312297, 4.17755472],
            [-4.08500466, -5.80817544, 6.92765242],
            [4.22020854, -0.01574608, 2.58052342],
            [0.28702242, 3.51322432, -1.91171157],
        ])

        self.run_fock_grad(mol, bas, min_bas, coulomb_coef, exchange_coef,
                           ref_grad, 1e-5)
