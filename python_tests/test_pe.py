from mpi4py import MPI
import unittest
import os

from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestPolEmbed(unittest.TestCase):

    def test_pe_dft(self):

        inpfile = os.path.join('inputs', 'pe_water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        potfile = os.path.join('inputs', 'pe_water.pot')
        if not os.path.isfile(potfile):
            potfile = os.path.join('python_tests', potfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['method_settings']['potfile'] = potfile

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        #    Final DFT energy:            -76.426499129387
        #    Nuclear repulsion:             9.194963912878
        #    Electronic energy:           -85.593468025478
        #    Embedding energy:             -0.027995016788

        ref_scf_energy = -76.426499129387
        ref_pe_energy = -0.027995016788

        scf_energy = scf_drv.get_scf_energy()
        self.assertTrue(abs(scf_energy - ref_scf_energy) < 1.0e-5)


if __name__ == "__main__":
    unittest.main()
