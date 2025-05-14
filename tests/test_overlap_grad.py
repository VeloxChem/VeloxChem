import numpy as np

from veloxchem.veloxchemlib import GradientScreeningData
from veloxchem.veloxchemlib import DenseMatrix
from veloxchem.veloxchemlib import AODensityMatrix, denmat
from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import compute_overlap_gradient_gpu
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestOverlapGradient:

    def run_overlap_grad(self, mol, bas, min_bas, ref_grad, tol):

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

        S_grad = compute_overlap_gradient_gpu(mol, bas, grad_screener, rank, nnodes)
        S_grad = S_grad.to_numpy()
        S_grad *= -2.0

        S_grad = scf_drv.comm.reduce(S_grad, root=mpi_master())

        if scf_drv.rank == 0:
            assert np.max(np.abs(S_grad - ref_grad)) < tol

    def test_overlap_grad_h2(self):

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
            [-0.20660057, -0.0050195, -0.],
            [0.1992455, -0.00362082, 0.],
            [0.00486902, -0.24699538, -0.],
            [0.0033819, 0.24647427, -0.],
            [-0.00044793, 0.00458072, -0.24058029],
            [-0.00044793, 0.00458072, 0.24058029],
        ])

        self.run_overlap_grad(mol, bas, min_bas, ref_grad, 1e-6)

    def test_overlap_grad_h2o(self):

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
            [0.27027957, 0.16916293, 0.25829908],
            [-0.26907165, 0.10417319, -0.08076347],
            [0.01913165, -0.23830326, -0.19001137],
            [0.17361041, -0.38051475, -0.15391178],
            [-0.21281791, 0.05922033, 0.27436262],
            [0.01883637, 0.32750915, -0.1282277],
            [0.32030438, 0.09853893, -0.245071],
            [-0.28644891, -0.1374864, -0.06885483],
            [-0.02631816, 0.00860775, 0.28682694],
            [-0.31499678, -0.2634244, 0.04609405],
            [0.29190199, 0.01584321, 0.14855116],
            [0.01558905, 0.23667333, -0.1472937],
        ])

        self.run_overlap_grad(mol, bas, min_bas, ref_grad, 1e-6)
