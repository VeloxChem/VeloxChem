from mpi4py import MPI
import veloxchem as vlx
import numpy as np
import unittest


class TestMP2(unittest.TestCase):

    def run_mp2(self, molname):

        comm = MPI.COMM_WORLD
        task = vlx.MpiTask(["inputs/" + molname + ".inp", None], comm)

        # scf
        scf_drv = vlx.ScfRestrictedDriver()
        scf_drv.compute_task(task)
        self.assertTrue(scf_drv.need_min_basis())

        # molecular orbitals
        if task.mpi_rank == vlx.mpi_master():
            mol_orbs = scf_drv.mol_orbs
        else:
            mol_orbs = vlx.MolecularOrbitals()

        # MO integrals
        moints_drv = vlx.MOIntegralsDriver()
        oovv = moints_drv.compute_task(task, mol_orbs, "OOVV")

        oooo = moints_drv.compute_task(task, mol_orbs, "OOOO")

        self.assertEqual(moints_drv.get_moints_type("OOOO"), vlx.moints.oooo)
        self.assertEqual(moints_drv.get_moints_type("OOOV"), vlx.moints.ooov)
        self.assertEqual(moints_drv.get_moints_type("OOVV"), vlx.moints.oovv)
        self.assertEqual(moints_drv.get_moints_type("OVOV"), vlx.moints.ovov)
        self.assertEqual(moints_drv.get_moints_type("OVVV"), vlx.moints.ovvv)
        self.assertEqual(moints_drv.get_moints_type("VVVV"), vlx.moints.vvvv)

        # MP2 energy
        nocc = task.molecule.number_of_alpha_electrons()
        e_mp2 = moints_drv.compute_mp2_energy(mol_orbs, nocc, oovv)

        # extra test: collect moints batches to master node
        moints_drv.collect_moints_batches(oovv)
        if task.mpi_rank == vlx.mpi_master():
            orb_ene = mol_orbs.ea_to_numpy()
            eocc = orb_ene[:nocc]
            evir = orb_ene[nocc:]
            eab = evir.reshape(-1, 1) + evir
            master_e_mp2 = 0.0
            # loop over generator pairs
            for pair in oovv.get_gen_pairs():
                ij = oovv.to_numpy(pair)
                ij_antisym = ij - ij.T
                denom = eocc[pair.first()] + eocc[pair.second()] - eab
                master_e_mp2 += np.sum(ij * (ij + ij_antisym) / denom)
            self.assertEqual(e_mp2, master_e_mp2)

        task.finish()

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
