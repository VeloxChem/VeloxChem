import numpy as np

from veloxchem.veloxchemlib import GradientScreeningData
from veloxchem.veloxchemlib import DenseMatrix
from veloxchem.veloxchemlib import AODensityMatrix, denmat
from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import compute_kinetic_energy_gradient_gpu
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestKineticEnergyGradient:

    def run_kinetic_energy_grad(self, mol, bas, min_bas, ref_grad, tol):

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

        num_gpus = scf_drv._get_num_gpus_per_node()

        rank = scf_drv.comm.Get_rank()
        nnodes = scf_drv.comm.Get_size()

        grad_screener = GradientScreeningData(mol, bas, dmat, wmat, num_gpus,
                                              1e-10, 1e-10, rank, nnodes)

        T_grad = compute_kinetic_energy_gradient_gpu(mol, bas, grad_screener, rank, nnodes)
        T_grad = T_grad.to_numpy()
        T_grad *= 2.0

        T_grad = scf_drv.comm.reduce(T_grad, root=mpi_master())

        if scf_drv.rank == 0:
            assert np.max(np.abs(T_grad - ref_grad)) < tol

    def test_kinetic_energy_grad_h2(self):

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

        ref_grad = np.array([
            [-0.43719321, -0.00079452, 0.],
            [0.40007615, -0.00035274, 0.],
            [0.01933955, -0.38827443, -0.],
            [0.01755299, 0.38909375, -0.],
            [0.00011227, 0.00016397, -0.38394741],
            [0.00011227, 0.00016397, 0.38394741],
        ])

        self.run_kinetic_energy_grad(mol, bas, min_bas, ref_grad, 1e-6)

    def test_kinetic_energy_grad_h2o(self):

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

        ref_grad = np.array([
            [0.58181038, 0.3320438, 0.60122203],
            [-0.64048474, 0.26216713, -0.18454763],
            [0.07462018, -0.56323495, -0.43026206],
            [0.3836813, -0.76986667, -0.29248991],
            [-0.44525676, 0.10260869, 0.57137658],
            [0.04662023, 0.67284307, -0.28368951],
            [0.67804207, 0.24284973, -0.54035004],
            [-0.63109763, -0.30245299, -0.17506406],
            [-0.04290258, 0.03332512, 0.69087138],
            [-0.65396588, -0.59430456, 0.0640736],
            [0.6377395, 0.01774016, 0.33739195],
            [0.01119393, 0.56628147, -0.35853232],
        ])

        self.run_kinetic_energy_grad(mol, bas, min_bas, ref_grad, 1e-6)
