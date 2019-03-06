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

    def test_small_molecules(self):

        mol_list = ["h2se"]

        scf_ener = [-2400.704613197391]

        for mol_name, e_ref in zip(mol_list, scf_ener):
            e_scf = self.run_rhf(mol_name)
            self.assertAlmostEqual(e_scf, e_ref, 10)


if __name__ == "__main__":
    unittest.main()
