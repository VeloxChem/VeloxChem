from mpi4py import MPI
import numpy as np
import unittest
import os

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestSCF(unittest.TestCase):

    def test_scf_hf(self):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        #    Final HF energy:             -76.041697549811
        ref_e_scf = -76.041697549811

        if task.mpi_rank == mpi_master():
            e_scf = scf_drv.get_scf_energy()
            self.assertTrue(np.max(np.abs(e_scf - ref_e_scf)) < 1.0e-6)

    def test_scf_dft(self):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None
        task.input_dict['method_settings']['xcfun'] = 'b3lyp'

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        #    Final DFT energy:            -76.443545741524
        ref_e_scf = -76.443545741524

        if task.mpi_rank == mpi_master():
            e_scf = scf_drv.get_scf_energy()
            self.assertTrue(np.max(np.abs(e_scf - ref_e_scf)) < 1.0e-5)

    def test_scf_dft_slda(self):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None
        task.input_dict['method_settings']['xcfun'] = 'slda'

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        #    Final DFT energy:            -76.074208234637
        ref_e_scf = -76.074208234637

        if task.mpi_rank == mpi_master():
            e_scf = scf_drv.get_scf_energy()
            self.assertTrue(np.max(np.abs(e_scf - ref_e_scf)) < 1.0e-5)

    def test_scf_hf_pe(self):

        try:
            import cppe
        except ImportError:
            return

        inpfile = os.path.join('inputs', 'pe_water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        potfile = os.path.join('inputs', 'pe_water.pot')
        if not os.path.isfile(potfile):
            potfile = os.path.join('python_tests', potfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None
        task.input_dict['method_settings']['potfile'] = potfile

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        #    Final HF energy:             -76.067159426565
        ref_e_scf = -76.067159426565

        if task.mpi_rank == mpi_master():
            e_scf = scf_drv.get_scf_energy()
            self.assertTrue(np.max(np.abs(e_scf - ref_e_scf)) < 1.0e-6)

    def test_scf_dft_pe(self):

        try:
            import cppe
        except ImportError:
            return

        inpfile = os.path.join('inputs', 'pe_water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        potfile = os.path.join('inputs', 'pe_water.pot')
        if not os.path.isfile(potfile):
            potfile = os.path.join('python_tests', potfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None
        task.input_dict['method_settings']['xcfun'] = 'b3lyp'
        task.input_dict['method_settings']['potfile'] = potfile

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        #    Final DFT energy:            -76.468733754150
        ref_e_scf = -76.468733754150

        if task.mpi_rank == mpi_master():
            e_scf = scf_drv.get_scf_energy()
            self.assertTrue(np.max(np.abs(e_scf - ref_e_scf)) < 1.0e-5)


if __name__ == "__main__":
    unittest.main()
