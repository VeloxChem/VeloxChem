from mpi4py import MPI
import unittest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import moints
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.mointsdriver import MOIntegralsDriver
from veloxchem.mp2driver import Mp2Driver
from veloxchem.mpitask import MpiTask


class TestMOIntegralsDriver(unittest.TestCase):

    def test_moints_type(self):

        moints_drv = MOIntegralsDriver()

        self.assertEqual(moints_drv.get_moints_type("OOOO"), moints.oooo)
        self.assertEqual(moints_drv.get_moints_type("OOOV"), moints.ooov)
        self.assertEqual(moints_drv.get_moints_type("OVOV"), moints.ovov)
        self.assertEqual(moints_drv.get_moints_type("OOVV"), moints.oovv)
        self.assertEqual(moints_drv.get_moints_type("OVVV"), moints.ovvv)
        self.assertEqual(moints_drv.get_moints_type("VVVV"), moints.vvvv)

    def test_h2se_mp2(self):

        # scf
        scf_drv = ScfRestrictedDriver()
        task = MpiTask(["inputs/h2se.inp", None], MPI.COMM_WORLD)

        scf_drv.compute_task(task)
        mol_orbs = scf_drv.mol_orbs

        # mp2
        mp2_drv = Mp2Driver()
        mp2_drv.compute_task(task, mol_orbs)

        if task.mpi_rank == mpi_master():
            e_mp2 = mp2_drv.e_mp2
            e_ref = -0.2852908836
            self.assertAlmostEqual(e_ref, e_mp2, 8)

        task.finish()


if __name__ == "__main__":
    unittest.main()
