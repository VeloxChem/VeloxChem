from mpi4py import MPI
import veloxchem as vlx
import unittest


class TestRHF(unittest.TestCase):

    def run_rhf(self, molname):

        comm = MPI.COMM_WORLD

        task = vlx.MpiTask(["inputs/" + molname + ".inp",
                            "inputs/" + molname + ".out"], comm)

        scf_drv = vlx.ScfRestrictedDriver()

        scf_drv.compute_task(task)

        task.finish()

        return scf_drv.get_scf_energy()

    def test_h2o(self):

        e_scf = self.run_rhf("h2o")

        self.assertAlmostEqual(e_scf, -75.922903268112, 10)

    def test_nh3(self):

        e_scf = self.run_rhf("nh3")

        self.assertAlmostEqual(e_scf, -56.195395860545, 10)

    def test_ch4(self):

        e_scf = self.run_rhf("ch4")

        self.assertAlmostEqual(e_scf, -40.155481408646, 10)

    def test_c2h4(self):

        e_scf = self.run_rhf("c2h4")

        self.assertAlmostEqual(e_scf, -78.043152739545, 10)


if __name__ == "__main__":
    unittest.main()
