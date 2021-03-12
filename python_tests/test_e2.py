from pathlib import Path
import numpy as np
import unittest

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lreigensolver import LinearResponseEigenSolver


class TestE2(unittest.TestCase):

    def test_E2_Be(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'be.inp')
        outfile = str(here / 'inputs' / 'be.out')

        task = MpiTask([inpfile, outfile])
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        scf_tensors = scf_drv.scf_tensors

        # E[2]

        lreig_solver = LinearResponseEigenSolver(task.mpi_comm, task.ostream)
        E2 = lreig_solver.get_e2(task.molecule, task.ao_basis, scf_tensors)

        if is_mpi_master(task.mpi_comm):

            n = E2.shape[0] // 2

            # 1s->3s, 2s->3s

            s_mat = np.zeros((4, 4))
            for ind_i, i in enumerate([4, 11, 4 + n, 11 + n]):
                for ind_k, k in enumerate([4, 11, 4 + n, 11 + n]):
                    s_mat[ind_i, ind_k] = E2[i - 1, k - 1]

            ref_s_mat = np.array([
                [9.52482128, -0.00864730, -0.04135515, 0.00944951],
                [-0.00864730, 1.06134246, 0.00944951, -0.09537676],
                [-0.04135515, 0.00944951, 9.52482128, -0.00864730],
                [0.00944951, -0.09537676, -0.00864730, 1.06134246],
            ])

            s_diag = np.diag(s_mat)
            ref_s_diag = np.diag(ref_s_mat)
            self.assertTrue(np.max(np.abs(s_diag - ref_s_diag)) < 1e-6)

            s_offdiag = s_mat - np.diag(s_diag)
            ref_s_offdiag = ref_s_mat - np.diag(ref_s_diag)
            self.assertTrue(
                np.max(np.abs(np.abs(s_offdiag) -
                              np.abs(ref_s_offdiag))) < 1e-6)

            # 1s->2p, 1s->3p, 2s->2p, 2s->3p

            ref_p_diag = np.array([
                8.85773484, 9.28555595, 0.40866678, 0.96062385, 8.85773484,
                9.28555595, 0.40866678, 0.96062385
            ])

            p_list = [1, 5, 8, 12, 1 + n, 5 + n, 8 + n, 12 + n]

            for d in range(3):
                p_diag = np.zeros(8)
                for ind_i, i in enumerate(p_list):
                    p_diag[ind_i] = E2[i + d - 1, i + d - 1]
                self.assertTrue(np.max(np.abs(p_diag - ref_p_diag)) < 1e-6)

            # eigenvalues

            ref_evals = np.array([
                0.18956763, 0.18956763, 0.18956763, 0.47989224, 0.47989224,
                0.47989224, 0.52851628, 4.34055925, 4.34055925, 4.34055925,
                4.73163430, 4.73163430, 4.73163430, 4.76236576
            ])

            S2 = np.zeros(E2.shape)
            for i in range(n):
                S2[i, i] = 1.0
                S2[i + n, i + n] = -1.0
            evals, evecs = np.linalg.eig((np.linalg.solve(E2 / 2.0, S2)))
            evals = np.sort(1.0 / evals.real)[n:]
            self.assertTrue(np.max(np.abs(evals - ref_evals)) < 1.0e-6)


if __name__ == "__main__":
    unittest.main()
