from mpi4py import MPI
import numpy as np
import unittest
import pytest
import sys
import os
try:
    import cppe
except ImportError:
    pass

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver


class TestSCF(unittest.TestCase):

    def run_scf(self, inpfile, potfile, xcfun_label, ref_e_scf):

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None

        if potfile is not None:
            task.input_dict['method_settings']['potfile'] = potfile

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        if task.mpi_rank == mpi_master():
            e_scf = scf_drv.get_scf_energy()
            tol = 1.0e-5 if xcfun_label is not None else 1.0e-6
            self.assertTrue(np.max(np.abs(e_scf - ref_e_scf)) < tol)

    def test_scf_hf(self):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        potfile = None

        xcfun_label = None

        #    Final HF energy:             -76.041697549811
        ref_e_scf = -76.041697549811

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf)

    def test_scf_dft(self):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        potfile = None

        xcfun_label = 'b3lyp'

        #    Final DFT energy:            -76.443545741524
        ref_e_scf = -76.443545741524

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf)

    def test_scf_dft_slda(self):

        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        potfile = None

        xcfun_label = 'slda'

        #    Final DFT energy:            -76.074208234637
        ref_e_scf = -76.074208234637

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf)

    @pytest.mark.skipif('cppe' not in sys.modules, reason='cppe not available')
    def test_scf_hf_pe(self):

        inpfile = os.path.join('inputs', 'pe_water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        potfile = os.path.join('inputs', 'pe_water.pot')
        if not os.path.isfile(potfile):
            potfile = os.path.join('python_tests', potfile)

        xcfun_label = None

        #    Final HF energy:             -76.067159426565
        ref_e_scf = -76.067159426565

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf)

    @pytest.mark.skipif('cppe' not in sys.modules, reason='cppe not available')
    def test_scf_dft_pe(self):

        inpfile = os.path.join('inputs', 'pe_water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        potfile = os.path.join('inputs', 'pe_water.pot')
        if not os.path.isfile(potfile):
            potfile = os.path.join('python_tests', potfile)

        xcfun_label = 'b3lyp'

        #    Final DFT energy:            -76.468733754150
        ref_e_scf = -76.468733754150

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf)


if __name__ == "__main__":
    unittest.main()
