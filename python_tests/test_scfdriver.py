from mpi4py import MPI
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.mpitask import MpiTask

import unittest


class TestScfDriver(unittest.TestCase):

    def test_h2se_scf(self):

        scf_drv = ScfRestrictedDriver()
        task = MpiTask(["inputs/h2se.inp", "inputs/h2se.out"], MPI.COMM_WORLD)

        scf_drv.compute_task(task)
        task.finish()

        e_scf = scf_drv.get_scf_energy()

        e_ref = -2400.704613197391

        self.assertAlmostEqual(e_ref, e_scf, 10)


if __name__ == "__main__":
    unittest.main()
