import numpy as np

from veloxchem.veloxchemlib import GradientScreeningData
from veloxchem.veloxchemlib import DenseMatrix
from veloxchem.veloxchemlib import AODensityMatrix, denmat
from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import compute_fock_hessian_gpu_2000
from veloxchem.veloxchemlib import compute_fock_hessian_gpu_1100
from veloxchem.molecule import Molecule
from veloxchem.molecularbasis import MolecularBasis
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestFockHessian:

    def run_fock_hessian(self, mol, bas, coulomb_coef, exchange_coef,
                         hessian_flag, ref_hessian, tol):

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

        if hessian_flag == 'hessian_2000':

            fock_hess_2000 = compute_fock_hessian_gpu_2000(
                mol, bas, dmat, coulomb_coef, [exchange_coef], [omega], 'symm',
                1e-12, 5e-6, grad_screener)

            fock_hess_2000 = fock_hess_2000.to_numpy()

            fock_hess_2000 = scf_drv.comm.reduce(fock_hess_2000,
                                                 root=mpi_master())

            if scf_drv.rank == 0:
                assert np.max(np.abs(fock_hess_2000 - ref_hessian)) < tol

        elif hessian_flag == 'hessian_1100':

            fock_hess_1100 = compute_fock_hessian_gpu_1100(
                mol, bas, dmat, coulomb_coef, [exchange_coef], [omega], 'symm',
                1e-12, 5e-6, grad_screener)

            fock_hess_1100 = fock_hess_1100.to_numpy()

            fock_hess_1100 = scf_drv.comm.reduce(fock_hess_1100,
                                                 root=mpi_master())

            if scf_drv.rank == 0:
                assert np.max(np.abs(fock_hess_1100 - ref_hessian)) < tol

    def test_fock_hessian_2000_methanol(self):

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
                              'hessian_2000', ref_hessian_2000, 1e-5)

    def test_fock_hessian_1100_h2(self):

        xyzstr = """6
        h2
        H   0.35  0.0   0.0
        H  -0.35  0.0   0.0
        H   1.5   0.35  0.0
        H   1.5  -0.35  0.0
        H   0.0   2.0   0.35
        H   0.0   2.0  -0.35
        """
        mol = Molecule.read_xyz_string(xyzstr)

        basis_label = 'sto-3g'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        coulomb_coef = 2.0

        # TODO: exchange Hessian
        exchange_coef = 0.0

        ref_hessian_1100_raw = np.array([
            [
                1.769523, -0.003658, 0., -0.003658, 1.692694, -0., 0., -0.,
                1.66782
            ],
            [
                -0.120988, -0.028446, 0., 0.0284, 0.848696, -0., -0., -0.,
                0.83178
            ],
            [
                0.064885, 0.040068, 0., 0.03936, -0.046006, 0., -0., 0.,
                -0.056305
            ],
            [
                0.056629, -0.038724, 0., -0.028748, -0.040158, 0., -0., 0.,
                -0.049193
            ],
            [
                -0.007254, -0.005549, -0.001145, -0.005828, 0.024341, 0.006653,
                -0.000955, 0.005289, -0.006895
            ],
            [
                -0.007254, -0.005549, 0.001145, -0.005828, 0.024341, -0.006653,
                0.000955, -0.005289, -0.006895
            ],
            [
                1.882644, 0.004264, -0., 0.004264, 1.83475, -0., -0., -0.,
                1.808871
            ],
            [
                0.134437, 0.036395, 0., 0.035816, -0.040527, 0., -0., 0.,
                -0.045454
            ],
            [
                0.125507, -0.040611, 0., -0.0248, -0.037779, 0., -0., 0.,
                -0.042504
            ],
            [
                -0.006295, 0.006808, 0.001355, 0.004116, 0.021542, 0.00592,
                0.000666, 0.004676, -0.006115
            ],
            [
                -0.006295, 0.006808, -0.001355, 0.004116, 0.021542, -0.00592,
                -0.000666, -0.004676, -0.006115
            ],
            [1.600701, 0.000878, 0., 0.000878, 1.60061, -0., 0., -0., 1.550943],
            [
                0.794244, -0.061827, 0., 0.051898, -0.138007, 0., 0., -0.,
                0.761352
            ],
            [
                -0.000188, 0.000425, 0.000103, 0.000351, -0.000246, -0.000102,
                0.000073, -0.000086, 0.000143
            ],
            [
                -0.000188, 0.000425, -0.000103, 0.000351, -0.000246, 0.000102,
                -0.000073, 0.000086, 0.000143
            ],
            [
                1.520871, -0.013897, 0., -0.013897, 1.522843, -0., 0., -0.,
                1.476557
            ],
            [
                0.000615, -0.002168, -0.000375, -0.001942, 0.002526, 0.000554,
                -0.000261, 0.000428, -0.000584
            ],
            [
                0.000615, -0.002168, 0.000375, -0.001942, 0.002526, -0.000554,
                0.000261, -0.000428, -0.000584
            ],
            [
                1.378984, -0.005619, -0.001049, -0.005619, 1.403001, 0.005039,
                -0.001049, 0.005039, 1.410151
            ],
            [
                0.696383, -0.003995, 0.01062, -0.003995, 0.713543, -0.041386,
                -0.01062, 0.041386, -0.13018
            ],
            [
                1.378984, -0.005619, 0.001049, -0.005619, 1.403001, -0.005039,
                0.001049, -0.005039, 1.410151
            ],
        ])

        natoms = mol.number_of_atoms()
        ref_hessian_1100 = np.zeros((natoms, natoms, 3, 3))

        ij_index = 0
        for i in range(natoms):
            for j in range(i, natoms):
                ref_hessian_1100[i, j] = ref_hessian_1100_raw[ij_index].reshape(
                    3, 3)
                ref_hessian_1100[j, i] = ref_hessian_1100_raw[ij_index].reshape(
                    3, 3).T
                ij_index += 1

        ref_hessian_1100 = ref_hessian_1100.reshape(natoms * natoms, 3 * 3)

        self.run_fock_hessian(mol, bas, coulomb_coef, exchange_coef,
                              'hessian_1100', ref_hessian_1100, 1e-5)
