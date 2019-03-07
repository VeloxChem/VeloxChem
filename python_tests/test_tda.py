from mpi4py import MPI
import veloxchem as vlx
import unittest


class TestTDA(unittest.TestCase):

    def run_scf(self, molname):

        comm = MPI.COMM_WORLD
        task = vlx.MpiTask(["inputs/" + molname + ".inp", ""], comm)

        # scf
        scf_drv = vlx.ScfRestrictedDriver()
        scf_drv.compute_task(task)

        # molecular orbitals
        if task.mpi_rank == vlx.mpi_master():
            return task, scf_drv.mol_orbs
        else:
            return task, vlx.MolecularOrbitals()

    def run_rsp(self, molname):

        task, mol_orbs = self.run_scf(molname)

        # response
        rsp_drv = vlx.ResponseDriver()
        rsp_drv.compute_task(mol_orbs, task)
        task.finish()

    def run_tda(self, molname):

        task, mol_orbs = self.run_scf(molname)
        eri_thresh = 1.0e-15

        # TDA singlet excited states
        tda_exci = vlx.TDAExciDriver(task.mpi_rank, task.mpi_size)

        tda_exci.set_number_states(3)
        tda_exci.set_eri_threshold(eri_thresh)
        tda_exci.set_solver(1.0e-4, 50)

        eri_drv = vlx.ElectronRepulsionIntegralsDriver(
            task.mpi_rank, task.mpi_size, task.mpi_comm)

        qq_data = eri_drv.compute(vlx.ericut.qq, eri_thresh, task.molecule,
                                  task.ao_basis)

        tda_exci.compute(qq_data, mol_orbs, task.molecule, task.ao_basis,
                         task.mpi_comm, vlx.OutputStream())
        task.finish()

        if task.mpi_rank == vlx.mpi_master():
            reigs, rnorms = tda_exci.solver.get_eigenvalues()
            return reigs
        else:
            return None

    def test_small_molecules(self):

        mol_list = ["h2se"]
        tda_ener = [[0.207436, 0.257474, 0.368358]]

        #TODO: add real test for ResponseDriver
        for mol_name, e_ref in zip(mol_list, tda_ener):
            self.run_rsp(mol_name)

        for mol_name, e_ref in zip(mol_list, tda_ener):
            e_calc = self.run_tda(mol_name)
            if MPI.COMM_WORLD.Get_rank() == vlx.mpi_master():
                for e1, e2 in zip(e_calc, e_ref):
                    self.assertAlmostEqual(e1, e2, 6)


if __name__ == "__main__":
    unittest.main()
