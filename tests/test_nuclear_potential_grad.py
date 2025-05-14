import numpy as np

from veloxchem.veloxchemlib import GradientScreeningData
from veloxchem.veloxchemlib import DenseMatrix
from veloxchem.veloxchemlib import AODensityMatrix, denmat
from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import compute_nuclear_potential_gradient_gpu
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestNuclearPotentialGradient:

    def run_nuclear_potential_grad(self, mol, bas, min_bas, ref_grad, tol):

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

        V_grad = compute_nuclear_potential_gradient_gpu(mol, bas, grad_screener, rank, nnodes)
        V_grad = V_grad.to_numpy()
        V_grad *= 2.0

        V_grad = scf_drv.comm.reduce(V_grad, root=mpi_master())

        if scf_drv.rank == 0:
            assert np.max(np.abs(V_grad - ref_grad)) < tol

    def test_nuclear_potential_grad_h2(self):

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
            [1.85268673, -0.1728567, 0.],
            [-2.17425485, -0.19147729, 0.],
            [0.27068195, 1.93107893, 0.],
            [0.23656806, -2.17503966, 0.],
            [-0.09284094, 0.30414736, 1.97589694],
            [-0.09284094, 0.30414736, -1.97589694],
        ])

        self.run_nuclear_potential_grad(mol, bas, min_bas, ref_grad, 1e-6)

    def test_nuclear_potential_grad_h2o(self):

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
            [-12.07461035, -12.73118511, -4.42976943],
            [6.62566605, -3.43603286, 2.11344344],
            [-1.24839297, 5.46136757, 5.55989071],
            [4.40366432, 2.35968011, 6.41100034],
            [5.60814702, -1.61513024, -6.17753122],
            [0.01831149, -7.86766327, 3.35694751],
            [-9.83691015, 7.90222004, 13.5145974],
            [7.17378817, 4.35817698, 2.59313424],
            [0.27049611, 1.06264529, -6.84888563],
            [6.9020097, 10.35705267, -13.63084712],
            [-7.50161843, 0.09698743, -4.83653088],
            [-0.34055095, -5.94811859, 2.37455064],
        ])

        self.run_nuclear_potential_grad(mol, bas, min_bas, ref_grad, 1e-5)
