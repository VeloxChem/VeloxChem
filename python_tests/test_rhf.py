from mpi4py import MPI
import veloxchem as vlx
import unittest


class TestRHF(unittest.TestCase):

    def run_rhf(self, molname):

        comm = MPI.COMM_WORLD
        task = vlx.MpiTask(
            ["inputs/" + molname + ".inp", "inputs/" + molname + ".out"], comm)

        scf_drv = vlx.ScfRestrictedDriver()
        scf_drv.compute_task(task)

        # TODO: add real test for VisualizationDriver
        vis_drv = vlx.VisualizationDriver()

        if task.mpi_rank == vlx.mpi_master():
            mol_orbs = scf_drv.mol_orbs
            density = scf_drv.density

            nelec = task.molecule.number_of_electrons()
            homo = nelec // 2 - 1

            vis_grid = vis_drv.gen_grid(task.molecule, 3, 3, 3)
            vis_drv.write_cube("homo.cube", task.molecule, task.ao_basis,
                               mol_orbs, homo, "alpha", vis_grid)
            vis_drv.write_cube("density.cube", task.molecule, task.ao_basis,
                               density, 0, "alpha", vis_grid)

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
