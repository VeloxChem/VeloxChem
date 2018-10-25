from mpi4py import MPI
import veloxchem as vlx
import unittest


class TestRHF(unittest.TestCase):

    def run_rhf(self, molname):

        comm = MPI.COMM_WORLD

        task = vlx.GlobalTask("inputs/" + molname + ".inp",
                              "inputs/" + molname + ".out", comm)

        molecule = task.molecule
        ao_basis = task.ao_basis
        min_basis = task.min_basis
        ostream = task.ostream

        scf_drv = vlx.ScfRestrictedDriver()
        scf_drv.compute(molecule, ao_basis, min_basis, comm, ostream)

        return scf_drv.old_energy + molecule.nuclear_repulsion_energy()

    def test_caffeine(self):

        e_scf = self.run_rhf("caf")
        self.assertAlmostEqual(e_scf, -675.828421206170, 10)


if __name__ == "__main__":
    unittest.main()
