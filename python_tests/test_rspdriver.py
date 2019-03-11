from mpi4py import MPI
from veloxchem.veloxchemlib import ElectronRepulsionIntegralsDriver
from veloxchem.veloxchemlib import ericut
from veloxchem.veloxchemlib import mpi_master
from veloxchem.outputstream import OutputStream
from veloxchem.molecularorbitals import MolecularOrbitals
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tdaexcidriver import TDAExciDriver
from veloxchem.rspdriver import ResponseDriver
from veloxchem.mpitask import MpiTask

import numpy as np
import unittest


class TestRspDriver(unittest.TestCase):

    def test_h2se_tda(self):

        # scf
        task = MpiTask(["inputs/h2se.inp", None], MPI.COMM_WORLD)
        scf_drv = ScfRestrictedDriver()

        scf_drv.compute_task(task)
        mol_orbs = scf_drv.mol_orbs

        # response
        # TODO: add test for ResponseDriver
        rsp_drv = ResponseDriver()
        rsp_drv.compute_task(mol_orbs, task)

        # tda
        tda_exci = TDAExciDriver(task.mpi_rank, task.mpi_size)

        tda_exci.set_number_states(3)
        tda_exci.set_eri_threshold(1.0e-12)
        tda_exci.set_solver(1.0e-4, 50)

        eri_drv = ElectronRepulsionIntegralsDriver(task.mpi_rank, task.mpi_size,
                                                   task.mpi_comm)

        qq_data = eri_drv.compute(ericut.qq, 1.0e-12, task.molecule,
                                  task.ao_basis)

        tda_exci.compute(qq_data, mol_orbs, task.molecule, task.ao_basis,
                         task.mpi_comm, OutputStream())

        if task.mpi_rank == mpi_master():

            reigs, rnorms = tda_exci.solver.get_eigenvalues()

            ref_eigs = np.array([0.207436, 0.257474, 0.368358])

            self.assertTrue(np.max(np.abs(reigs - ref_eigs)) < 1.0e-6)

        task.finish()


if __name__ == "__main__":
    unittest.main()
