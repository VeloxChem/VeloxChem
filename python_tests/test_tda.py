from mpi4py import MPI
import veloxchem as vlx
import unittest


class TestTDA(unittest.TestCase):

    def run_tda(self, molname):

        # scf

        comm = MPI.COMM_WORLD
        task = vlx.MpiTask(["inputs/" + molname + ".inp",
                            "inputs/" + molname + ".out"], comm)

        scf_drv = vlx.ScfRestrictedDriver()
        scf_drv.compute_task(task)
        task.finish()

        # molecular orbitals

        if task.mpi_rank == vlx.mpi_master():
            mol_orbs = scf_drv.mol_orbs
        else:
            mol_orbs = vlx.MolecularOrbitals()

        # TDA singlet excited states

        tda_exci = vlx.TDAExciDriver(task.mpi_rank, task.mpi_size)

        tda_exci.set_number_states(3)
        tda_exci.set_eri_threshold(scf_drv.eri_thresh)
        tda_exci.set_solver(1.0e-4, 50)

        eri_drv = vlx.ElectronRepulsionIntegralsDriver(
            task.mpi_rank, task.mpi_size, task.mpi_comm)

        qq_data = eri_drv.compute(
            vlx.ericut.qqden, scf_drv.eri_thresh, task.molecule, task.ao_basis)

        tda_exci.compute(qq_data, mol_orbs, task.molecule, task.ao_basis,
                         task.mpi_comm, vlx.OutputStream())

        if task.mpi_rank == vlx.mpi_master():
            reigs, rnorms = tda_exci.solver.get_eigenvalues()
            return reigs
        else:
            return None

    def test_small_molecules(self):

        mol_list = ["h2se"]

        tda_ener = [[0.207436, 0.257474, 0.368358]]

        for mol_name, e_ref in zip(mol_list, tda_ener):
            e_calc = self.run_tda(mol_name)
            if MPI.COMM_WORLD.Get_rank() == vlx.mpi_master():
                for e1, e2 in zip(e_calc, e_ref):
                    self.assertAlmostEqual(e1, e2, 6)


if __name__ == "__main__":
    unittest.main()
