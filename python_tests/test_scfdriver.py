from mpi4py import MPI
import unittest
import os

from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
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

        e_scf = scf_drv.get_scf_energy()

        self.assertAlmostEqual(-2400.704613197391, e_scf, 10)

        # unrestricted scf

        task.molecule.set_charge(1)
        task.molecule.set_multiplicity(2)
        task.molecule.check_multiplicity()

        scf_unrest_drv = ScfUnrestrictedDriver(task.mpi_comm, task.ostream)
        scf_unrest_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        e_uhf = scf_unrest_drv.get_scf_energy()
        s2 = scf_unrest_drv.compute_s2(task.molecule,
                                       scf_unrest_drv.scf_tensors['S'],
                                       scf_unrest_drv.mol_orbs)

        self.assertAlmostEqual(-2400.38319890, e_uhf, 8)
        self.assertAlmostEqual(0.7619, s2, 4)

        task.finish()


if __name__ == "__main__":
    unittest.main()
