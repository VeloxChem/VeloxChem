import numpy as np

from veloxchem.veloxchemlib import GradientScreeningData
from veloxchem.veloxchemlib import DenseMatrix
from veloxchem.veloxchemlib import AODensityMatrix, denmat
from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import compute_fock_hessian_gpu_2000
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestFockHessian:

    def run_fock_hessian(self, mol, bas, coulomb_coef, exchange_coef,
                         ref_hessian, tol):

        scf_drv = ScfRestrictedDriver()
        scf_drv.filename = None
        scf_drv.checkpoint_file = None
        scf_drv.ostream.mute()
        scf_results = scf_drv.compute(mol, bas)

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

        fock_hess_2000 = compute_fock_hessian_gpu_2000(mol, bas, dmat,
                                                       coulomb_coef,
                                                       [exchange_coef], [omega],
                                                       'symm', 1e-12, 5e-6,
                                                       grad_screener)

        fock_hess_2000 = fock_hess_2000.to_numpy()

        fock_hess_2000 = scf_drv.comm.reduce(fock_hess_2000, root=mpi_master())

        if scf_drv.rank == 0:
            assert np.max(np.abs(fock_hess_2000 - ref_hessian)) < tol

    def test_fock_hessian_methanol(self):

        xyzstr = """6
        methanol
        H      1.2001      0.0363      0.8431
        C      0.7031      0.0083     -0.1305
        H      0.9877      0.8943     -0.7114
        H      1.0155     -0.8918     -0.6742
        O     -0.6582     -0.0067      0.1730
        H     -1.1326     -0.0311     -0.6482
        """
        mol = Molecule.read_xyz_string(xyzstr)

        basis_label = 'sto-3g'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        coulomb_coef = 2.0

        # TODO: exchange Hessian
        exchange_coef = 0.0

        ref_hessian_2000 = np.array([
            [-6.204107, 0.0255, 0.882136, -6.662291, 0.044925, -5.164297],
            [
                -830.258383, 0.032366, -1.984103, -834.93212, -0.032633,
                -835.002495
            ],
            [-6.558089, 0.551889, -0.350025, -5.54811, -0.825046, -6.204037],
            [-6.522697, -0.597758, -0.351035, -5.51711, 0.784687, -6.271384],
            [
                -2109.21271, 1.121386, 8.617113, -2169.445436, 0.942879,
                -2130.384044
            ],
            [-7.177555, 0.025357, 0.963311, -7.524177, 0.044137, -6.119375],
        ])

        self.run_fock_hessian(mol, bas, coulomb_coef, exchange_coef,
                              ref_hessian_2000, 1e-5)
