from mpi4py import MPI
import unittest
import os

from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.mpitask import MpiTask


class TestScfDriver(unittest.TestCase):

    def test_h2se_scf(self):

        inpfile = os.path.join('inputs', 'h2se.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)
        outfile = inpfile.replace('.inp', '.out')

        task = MpiTask([inpfile, outfile], MPI.COMM_WORLD)
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        task.finish()

        e_scf = scf_drv.get_scf_energy()

        e_ref = -2400.704613197391

        self.assertAlmostEqual(e_ref, e_scf, 10)


if __name__ == "__main__":
    unittest.main()
