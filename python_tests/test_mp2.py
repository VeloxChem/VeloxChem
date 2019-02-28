from mpi4py import MPI
import veloxchem as vlx
import unittest


class TestMP2(unittest.TestCase):

    def run_mp2(self, molname):

        comm = MPI.COMM_WORLD

        task = vlx.MpiTask(["inputs/" + molname + ".inp",
                            "inputs/" + molname + ".out"], comm)

        # scf

        scf_drv = vlx.ScfRestrictedDriver()

        scf_drv.compute_task(task)

        task.finish()

        # molecular orbitals

        if task.mpi_rank == vlx.mpi_master():
            mol_orbs = scf_drv.mol_orbs
        else:
            mol_orbs = vlx.MolecularOrbitals()

        # MO integrals

        moints_drv = vlx.MOIntegralsDriver()
        oovv = moints_drv.compute_task(task, mol_orbs, "OOVV")

        # MP2 energy

        nocc = task.molecule.number_of_alpha_electrons()
        e_mp2 = moints_drv.compute_mp2_energy(mol_orbs, nocc, oovv)

        return e_mp2

    def test_small_molecules(self):

        mol_list = ["h2se"]

        mp2_ener = [-0.2852908836]

        for mol_name, e_ref in zip(mol_list, mp2_ener):
            e_mp2 = self.run_mp2(mol_name)
            if MPI.COMM_WORLD.Get_rank() == vlx.mpi_master():
                self.assertAlmostEqual(e_mp2, e_ref, 8)


if __name__ == "__main__":
    unittest.main()
