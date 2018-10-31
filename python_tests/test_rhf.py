from mpi4py import MPI
import veloxchem as vlx
import unittest


class TestRHF(unittest.TestCase):

    def run_rhf(self, molname):

        comm = MPI.COMM_WORLD

        task = vlx.GlobalTask("inputs/" + molname + ".inp",
                              "inputs/" + molname + ".out", comm)

        enuc = task.molecule.nuclear_repulsion_energy()

        scf_drv = vlx.ScfRestrictedDriver()

        scf_drv.compute_task(task, comm)

        return scf_drv.get_scf_energy(enuc)

    def test_caffeine(self):

        e_scf = self.run_rhf("caf")

        self.assertAlmostEqual(e_scf, -675.828421206170, 10)


if __name__ == "__main__":
    unittest.main()
