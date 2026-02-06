import numpy as np

from veloxchem.veloxchemlib import GradientScreeningData
from veloxchem.veloxchemlib import DenseMatrix
from veloxchem.veloxchemlib import AODensityMatrix, denmat
from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import compute_fock_hessian_gpu_2000
from veloxchem.veloxchemlib import compute_fock_hessian_gpu_1100
from veloxchem.veloxchemlib import compute_fock_hessian_gpu_1010
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

        elif hessian_flag == 'hessian_1010':

            fock_hess_1010 = compute_fock_hessian_gpu_1010(
                mol, bas, dmat, coulomb_coef, [exchange_coef], [omega], 'symm',
                1e-12, 5e-6, grad_screener)

            fock_hess_1010 = fock_hess_1010.to_numpy()

            fock_hess_1010 = scf_drv.comm.reduce(fock_hess_1010,
                                                 root=mpi_master())

            if scf_drv.rank == 0:
                assert np.max(np.abs(fock_hess_1010 - ref_hessian)) < tol

    def test_hessian_contributions_methanol(self):

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

        natoms = mol.number_of_atoms()

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

        ref_hessian_1100_triu = np.array([
            [
                3.967622, 0.004698, 0.125483, 0.004698, 3.834203, 0.006106,
                0.125483, 0.006106, 4.033421
            ],
            [
                1.885454, -0.007854, -0.316796, -0.005484, 1.967358, -0.01508,
                -0.202384, -0.01328, 1.479717
            ],
            [
                -0.066431, -0.029054, 0.043305, -0.005729, -0.008353, -0.092388,
                0.008836, -0.123857, 0.128773
            ],
            [
                -0.067489, 0.028787, 0.039121, 0.004027, 0.001959, 0.098978,
                0.004688, 0.129373, 0.119516
            ],
            [
                0.313557, 0.007908, 0.065506, 0.009476, -0.105475, 0.001292,
                0.14074, 0.002466, -0.092694
            ],
            [
                -0.050386, -0.001858, -0.042377, -0.001849, 0.012201, -0.001224,
                -0.041897, -0.001216, -0.014808
            ],
            [
                695.926305, 0.010603, 0.884596, 0.010603, 696.415987, 0.025144,
                0.884596, 0.025144, 696.956777
            ],
            [
                1.866403, -0.172974, 0.109325, -0.222112, 1.526739, 0.262969,
                0.135716, 0.24889, 1.751848
            ],
            [
                1.853816, 0.187879, 0.109674, 0.238141, 1.518887, -0.24954,
                0.133578, -0.234644, 1.772748
            ],
            [
                3.096913, -0.002003, -0.303376, 0.001362, 2.819333, -0.006405,
                -0.13859, -0.003843, 2.739568
            ],
            [
                0.377934, 0.013927, 0.243814, 0.011633, -0.190259, 0.004708,
                0.132885, 0.002974, -0.147172
            ],
            [
                4.046691, 0.082074, -0.06274, 0.082074, 4.151042, -0.121507,
                -0.06274, -0.121507, 4.057745
            ],
            [
                -0.08372, 0.015749, 0.002993, -0.025457, 0.232955, 0.017734,
                0.004607, -0.030826, -0.083986
            ],
            [
                0.160512, 0.022918, -0.087387, 0.093171, -0.093527, -0.024663,
                -0.118833, -0.002945, -0.040314
            ],
            [
                0.032297, 0.016123, 0.003246, 0.021187, -0.001507, 0.001777,
                -0.00487, -0.001614, -0.010284
            ],
            [
                4.052133, -0.087931, -0.062899, -0.087931, 4.155744, 0.115363,
                -0.062899, 0.115363, 4.048017
            ],
            [
                0.164096, -0.018556, -0.08733, -0.090065, -0.09596, 0.024195,
                -0.115148, 0.001546, -0.041448
            ],
            [
                0.033412, -0.014877, 0.003949, -0.020273, -0.002657, -0.002005,
                -0.004047, 0.001133, -0.010257
            ],
            [
                1768.66154, -1.326603, -11.738791, -1.326603, 1837.873493,
                -1.186532, -11.738791, -1.186532, 1789.392537
            ],
            [
                2.338379, -0.000802, 0.028691, 0.00083, 2.428036, 0.012696,
                0.10753, 0.013926, 3.019607
            ],
            [
                4.294822, 0.006667, 0.188952, 0.006667, 4.118586, 0.009054,
                0.188952, 0.009054, 4.413266
            ],
        ])

        ref_hessian_1100 = np.zeros((natoms, natoms, 3, 3))

        ij_index = 0
        for i in range(natoms):
            for j in range(i, natoms):
                hess_ij = ref_hessian_1100_triu[ij_index]
                ref_hessian_1100[i, j] = hess_ij.reshape(3, 3)
                ref_hessian_1100[j, i] = hess_ij.reshape(3, 3).T
                ij_index += 1

        ref_hessian_1100 = ref_hessian_1100.reshape(natoms * natoms, 3 * 3)

        self.run_fock_hessian(mol, bas, coulomb_coef, exchange_coef,
                              'hessian_1100', ref_hessian_1100, 1e-5)

        ref_hessian_1010_triu = np.array([
            [
                0.28842, 0.001408, 0.051288, 0.001408, 0.266544, 0.002903,
                0.051288, 0.002903, 0.367631
            ],
            [
                0.209834, -0.019122, -0.658809, -0.019277, 0.556316, -0.037574,
                -0.666351, -0.037694, -0.755336
            ],
            [
                0.022839, 0.014516, -0.022681, 0.000476, -0.005117, 0.04019,
                -0.002773, 0.06714, -0.073597
            ],
            [
                0.023291, -0.014465, -0.020771, 0.000406, -0.010165, -0.043517,
                -0.000643, -0.069847, -0.069001
            ],
            [
                -0.315744, -0.010034, -0.133327, -0.013215, 0.14811, -0.004075,
                -0.286747, -0.006479, 0.052726
            ],
            [
                -0.006858, -0.000431, -0.01208, -0.000438, 0.004711, -0.000537,
                -0.012378, -0.000541, -0.012052
            ],
            [
                129.223131, -0.00015, -0.313449, -0.00015, 128.822018,
                -0.006452, -0.313449, -0.006452, 128.749042
            ],
            [
                0.440202, -0.332423, 0.219072, -0.345822, -0.5153, 0.702255,
                0.224413, 0.692769, 0.104745
            ],
            [
                0.41788, 0.371681, 0.225771, 0.385234, -0.549631, -0.667596,
                0.230228, -0.657972, 0.161416
            ],
            [
                -4.870016, -0.077555, 1.79777, -0.07839, 2.467964, 0.021884,
                1.756562, 0.021244, 2.16763
            ],
            [
                -0.169473, -0.00659, -0.12036, -0.005887, 0.092708, -0.00336,
                -0.08623, -0.002823, 0.021512
            ],
            [
                0.28562, 0.028984, -0.018252, 0.028984, 0.356983, -0.052598,
                -0.018252, -0.052598, 0.31248
            ],
            [
                0.02533, -0.00882, -0.000433, 0.012771, -0.102801, -0.012667,
                -0.001359, 0.017948, 0.024521
            ],
            [
                -0.154014, -0.126145, 0.155544, -0.226431, 0.006015, 0.142115,
                0.1972, 0.10124, 0.02061
            ],
            [
                -0.017641, -0.009585, -0.009137, -0.018253, 0.000982, -0.010086,
                0.00525, 0.002322, 0.011497
            ],
            [
                0.287516, -0.031974, -0.018574, -0.031974, 0.359368, 0.050028,
                -0.018574, 0.050028, 0.308209
            ],
            [
                -0.165081, 0.127394, 0.154618, 0.229345, 0.007249, -0.136496,
                0.190787, -0.094385, 0.030492
            ],
            [
                -0.018485, 0.00858, -0.009829, 0.017846, 0.002163, 0.010197,
                0.00456, -0.001764, 0.011158
            ],
            [
                339.867138, 0.318783, 2.759886, 0.318783, 323.15526, 0.280239,
                2.759886, 0.280239, 334.615163
            ],
            [
                0.11543, -0.038016, -1.29287, -0.038324, 0.834938, -0.068126,
                -1.307756, -0.068358, -1.479833
            ],
            [
                0.248125, 0.001505, 0.05469, 0.001505, 0.224274, 0.002769,
                0.05469, 0.002769, 0.316741
            ],
        ])

        ref_hessian_1010 = np.zeros((natoms, natoms, 3, 3))

        ij_index = 0
        for i in range(natoms):
            for j in range(i, natoms):
                hess_ij = ref_hessian_1010_triu[ij_index]
                ref_hessian_1010[i, j] = hess_ij.reshape(3, 3)
                ref_hessian_1010[j, i] = hess_ij.reshape(3, 3).T
                ij_index += 1

        ref_hessian_1010 = ref_hessian_1010.reshape(natoms * natoms, 3 * 3)

        self.run_fock_hessian(mol, bas, coulomb_coef, exchange_coef,
                              'hessian_1010', ref_hessian_1010, 1e-5)

    def test_Coulomb_hessian_methanol_blyp(self):

        ref_hessian = np.array([[[[-2.219942, 0.033711, 1.183905],
                                  [0.033711, -2.800731, 0.061748],
                                  [1.183905, 0.061748, -0.721549]],
                                 [[2.376752, -0.026581, -0.966563],
                                  [-0.023792, 2.792364, -0.054394],
                                  [-0.831826, -0.052268, 0.90246]],
                                 [[-0.036467, -0.00612, 0.008008],
                                  [-0.00104, -0.027333, -0.014843],
                                  [-0.00234, -0.004467, -0.015261]],
                                 [[-0.036676, 0.006129, 0.007276],
                                  [0.000632, -0.026288, 0.015534],
                                  [-0.002528, 0.00483, -0.016106]],
                                 [[0.017616, -0.003057, -0.136335],
                                  [-0.00533, 0.029933, -0.004912],
                                  [-0.246341, -0.006643, -0.103314]],
                                 [[-0.101283, -0.004082, -0.096291],
                                  [-0.00418, 0.032056, -0.003133],
                                  [-0.100871, -0.003201, -0.046231]]],
                                [[[2.376752, -0.023792, -0.831826],
                                  [-0.026581, 2.792364, -0.052268],
                                  [-0.966563, -0.054394, 0.90246]],
                                 [[-4.305567, 0.078127, -1.689199],
                                  [0.078127, -11.515956, -0.014799],
                                  [-1.689199, -0.014799, -10.954168]],
                                 [[2.532114, -0.476605, 0.30865],
                                  [-0.545019, 1.194878, 1.001125],
                                  [0.344004, 0.97109, 2.073888]],
                                 [[2.498872, 0.531382, 0.317975],
                                  [0.60123, 1.148428, -0.95209],
                                  [0.349513, -0.92098, 2.154202]],
                                 [[-2.933636, -0.104489, 1.992397],
                                  [-0.105176, 6.383006, 0.02007],
                                  [1.961041, 0.019563, 5.850329]],
                                 [[-0.168535, -0.004623, -0.097997],
                                  [-0.00258, -0.00272, -0.002038],
                                  [0.001205, -0.00048, -0.026711]]],
                                [[[-0.036467, -0.00104, -0.00234],
                                  [-0.00612, -0.027333, -0.004467],
                                  [0.008008, -0.014843, -0.015261]],
                                 [[2.532114, -0.545019, 0.344004],
                                  [-0.476605, 1.194878, 0.97109],
                                  [0.30865, 1.001125, 2.073888]],
                                 [[-2.462361, 0.701816, -0.417865],
                                  [0.701816, -1.030319, -1.135722],
                                  [-0.417865, -1.135722, -1.897734]],
                                 [[-0.085126, 0.01418, 0.004671],
                                  [-0.019408, 0.08937, 0.007602],
                                  [0.005719, -0.014627, -0.083156]],
                                 [[0.022025, -0.185445, 0.076678],
                                  [-0.202232, -0.22155, 0.176506],
                                  [0.091281, 0.163881, -0.077057]],
                                 [[0.029814, 0.015508, -0.005149],
                                  [0.002548, -0.005047, -0.01501],
                                  [0.004207, 0.000186, -0.000681]]],
                                [[[-0.036676, 0.000632, -0.002528],
                                  [0.006129, -0.026288, 0.00483],
                                  [0.007276, 0.015534, -0.016106]],
                                 [[2.498872, 0.60123, 0.349513],
                                  [0.531382, 1.148428, -0.92098],
                                  [0.317975, -0.95209, 2.154202]],
                                 [[-0.085126, -0.019408, 0.005719],
                                  [0.01418, 0.08937, -0.014627],
                                  [0.004671, 0.007602, -0.083156]],
                                 [[-2.417022, -0.763835, -0.421819],
                                  [-0.763835, -0.983645, 1.083594],
                                  [-0.421819, 1.083594, -1.99035]],
                                 [[0.009636, 0.195991, 0.074107],
                                  [0.213398, -0.222914, -0.16784],
                                  [0.087605, -0.154751, -0.063299]],
                                 [[0.030315, -0.01461, -0.004991],
                                  [-0.001254, -0.004952, 0.015023],
                                  [0.004292, 0.000111, -0.001292]]],
                                [[[0.017616, -0.00533, -0.246341],
                                  [-0.003057, 0.029933, -0.006643],
                                  [-0.136335, -0.004912, -0.103314]],
                                 [[-2.933636, -0.105176, 1.961041],
                                  [-0.104489, 6.383006, 0.019563],
                                  [1.992397, 0.02007, 5.850329]],
                                 [[0.022025, -0.202232, 0.091281],
                                  [-0.185445, -0.22155, 0.163881],
                                  [0.076678, 0.176506, -0.077057]],
                                 [[0.009636, 0.213398, 0.087605],
                                  [0.195991, -0.222914, -0.154751],
                                  [0.074107, -0.16784, -0.063299]],
                                 [[-0.692689, 0.16294, 0.134498],
                                  [0.16294, -10.927883, 0.080517],
                                  [0.134498, 0.080517, -7.142008]],
                                 [[3.577048, -0.063601, -2.028084],
                                  [-0.065941, 4.959408, -0.102567],
                                  [-2.141344, -0.10434, 1.535348]]],
                                [[[-0.101283, -0.00418, -0.100871],
                                  [-0.004082, 0.032056, -0.003201],
                                  [-0.096291, -0.003133, -0.046231]],
                                 [[-0.168535, -0.00258, 0.001205],
                                  [-0.004623, -0.00272, -0.00048],
                                  [-0.097997, -0.002038, -0.026711]],
                                 [[0.029814, 0.002548, 0.004207],
                                  [0.015508, -0.005047, 0.000186],
                                  [-0.005149, -0.01501, -0.000681]],
                                 [[0.030315, -0.001254, 0.004292],
                                  [-0.01461, -0.004952, 0.000111],
                                  [-0.004991, 0.015023, -0.001292]],
                                 [[3.577048, -0.065941, -2.141344],
                                  [-0.063601, 4.959408, -0.10434],
                                  [-2.028084, -0.102567, 1.535348]],
                                 [[-3.367359, 0.071407, 2.232512],
                                  [0.071407, -4.978746, 0.107724],
                                  [2.232512, 0.107724, -1.460433]]]])

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

        basis_label = 'def2-svp'
        bas = MolecularBasis.read(mol, basis_label, ostream=None)

        xcfun_label = 'blyp'
        coulomb_coef = 2.0
        exchange_coef = 0.0
        omega = 0.0

        natoms = mol.number_of_atoms()

        scf_drv = ScfRestrictedDriver()
        scf_drv.filename = None
        scf_drv.checkpoint_file = None
        scf_drv.xcfun = xcfun_label
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

        fock_hess_1100 = compute_fock_hessian_gpu_1100(mol, bas, dmat,
                                                       coulomb_coef,
                                                       [exchange_coef], [omega],
                                                       'symm', 1e-12, 5e-6,
                                                       grad_screener)
        fock_hess_1100 = fock_hess_1100.to_numpy()
        fock_hess_1100 = scf_drv.comm.reduce(fock_hess_1100, root=mpi_master())

        fock_hess_1010 = compute_fock_hessian_gpu_1010(mol, bas, dmat,
                                                       coulomb_coef,
                                                       [exchange_coef], [omega],
                                                       'symm', 1e-12, 5e-6,
                                                       grad_screener)
        fock_hess_1010 = fock_hess_1010.to_numpy()
        fock_hess_1010 = scf_drv.comm.reduce(fock_hess_1010, root=mpi_master())

        if scf_drv.rank == mpi_master():
            j_hess = np.zeros((natoms, natoms, 3, 3))

            for a in range(natoms):
                atom_hess = np.zeros((3, 3))
                index = 0
                for x in range(3):
                    for y in range(x, 3):
                        atom_hess[x, y] = fock_hess_2000[a][index]
                        if x != y:
                            atom_hess[y, x] = atom_hess[x, y]
                        index += 1
                j_hess[a, a] += atom_hess

            j_hess += fock_hess_1100.reshape(natoms, natoms, 3, 3)
            j_hess += fock_hess_1010.reshape(natoms, natoms, 3, 3)

            assert np.max(np.abs(j_hess - ref_hessian)) < 1.0e-5
