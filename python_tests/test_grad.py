from mpi4py import MPI
import numpy as np
import unittest
import os

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.gradientdriver import GradientDriver


class TestSCF(unittest.TestCase):

    def run_grad(self, inpfile, ref_grad):

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None

        grad_drv = GradientDriver(task.mpi_comm, task.ostream)
        grad_drv.update_settings(task.input_dict['scf'],
                                 task.input_dict['method_settings'])
        grad_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if task.mpi_rank == mpi_master():
            grad = grad_drv.get_gradient()
            self.assertTrue(np.max(np.abs(grad - ref_grad)) < 1.0e-6)

    def test_nh3(self):

        inpfile = os.path.join('inputs', 'nh3.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        ref_grad = np.array([
            [0.0133714, -0.0004248, -0.0062534, -0.0066932],
            [-0.0056538, 0.0113812, -0.0033086, -0.0024188],
            [-0.0006128, 0.0005803, 0.0089449, -0.0089123],
        ]).T

        self.run_grad(inpfile, ref_grad)


if __name__ == "__main__":
    unittest.main()
