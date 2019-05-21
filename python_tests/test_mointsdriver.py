from mpi4py import MPI
import numpy as np
import unittest
import os

from veloxchem.veloxchemlib import mpi_master
from veloxchem.veloxchemlib import moints
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.mointsdriver import MOIntegralsDriver
from veloxchem.mp2driver import Mp2Driver
from veloxchem.mpitask import MpiTask
from veloxchem.outputstream import OutputStream


class TestMOIntegralsDriver(unittest.TestCase):

    def test_moints_type(self):

        moints_drv = MOIntegralsDriver(MPI.COMM_WORLD, OutputStream())

        self.assertEqual(moints_drv.get_moints_type("OOOO"), moints.oooo)
        self.assertEqual(moints_drv.get_moints_type("OOOV"), moints.ooov)
        self.assertEqual(moints_drv.get_moints_type("OVOV"), moints.ovov)
        self.assertEqual(moints_drv.get_moints_type("OOVV"), moints.oovv)
        self.assertEqual(moints_drv.get_moints_type("OVVV"), moints.ovvv)
        self.assertEqual(moints_drv.get_moints_type("VVVV"), moints.vvvv)

    def test_h2se_mp2(self):

        # scf
        inpfile = os.path.join('inputs', 'h2se.inp')
        task = MpiTask([inpfile, None], MPI.COMM_WORLD)

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        mol_orbs = scf_drv.mol_orbs

        # mp2
        mp2_drv = Mp2Driver(task.mpi_comm, task.ostream)
        mp2_drv.compute(task.molecule, task.ao_basis, mol_orbs)

        if task.mpi_rank == mpi_master():
            e_mp2 = mp2_drv.e_mp2
            e_ref = -0.2852908836
            self.assertAlmostEqual(e_ref, e_mp2, 8)

        # extra test: collect moints batches to master node
        moints_drv = MOIntegralsDriver(task.mpi_comm, task.ostream)
        grps = [p for p in range(task.mpi_size)]
        oovv = moints_drv.compute(task.molecule, task.ao_basis, mol_orbs,
                                  "OOVV", grps)
        moints_drv.collect_moints_batches(oovv, grps)

        if task.mpi_rank == mpi_master():
            orb_ene = mol_orbs.ea_to_numpy()
            nocc = task.molecule.number_of_alpha_electrons()
            eocc = orb_ene[:nocc]
            evir = orb_ene[nocc:]
            eab = evir.reshape(-1, 1) + evir

            # loop over generator pairs for mp2 energy
            master_e_mp2 = 0.0
            for pair in oovv.get_gen_pairs():
                ij = oovv.to_numpy(pair)
                ij_antisym = ij - ij.T
                denom = eocc[pair.first()] + eocc[pair.second()] - eab
                master_e_mp2 += np.sum(ij * (ij + ij_antisym) / denom)

            self.assertAlmostEqual(e_mp2, master_e_mp2, 10)

        task.finish()


if __name__ == "__main__":
    unittest.main()
