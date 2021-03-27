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

            ref_evals = np.array([
                0.18956763, 0.18956763, 0.18956763, 0.47989224, 0.47989224,
                0.47989224, 0.52851628, 4.34055925, 4.34055925, 4.34055925,
                4.73163430, 4.73163430, 4.73163430, 4.76236576
            ])

            half_size = E2.shape[0] // 2
            S2 = np.diag([1.0] * half_size + [-1.0] * half_size)

            evals, evecs = np.linalg.eig(np.matmul(S2, E2))
            evals = np.sort(evals.real)[half_size:]
            self.assertTrue(np.max(np.abs(evals - ref_evals)) < 1.0e-6)


if __name__ == "__main__":
    unittest.main()
