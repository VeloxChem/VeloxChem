from mpi4py import MPI
import numpy as np
import unittest
import pytest
import sys
from pathlib import Path
try:
    import cppe
except ImportError:
    pass

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scffirstorderprop import ScfFirstOrderProperties


class TestSCF(unittest.TestCase):

    def run_scf(self, inpfile, potfile, xcfun_label, ref_e_scf, ref_dip):

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

        scf_prop = ScfFirstOrderProperties(task.mpi_comm, task.ostream)
        scf_prop.compute(task.molecule, task.ao_basis, scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            e_scf = scf_drv.get_scf_energy()
            tol = 1.0e-5 if xcfun_label is not None else 1.0e-6
            self.assertTrue(np.max(np.abs(e_scf - ref_e_scf)) < tol)

            dip = np.linalg.norm(scf_prop.get_property('dipole moment'))
            if ref_dip is not None:
                self.assertTrue(np.max(np.abs(dip - ref_dip)) < tol)

    def test_scf_hf(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')

        potfile = None

        xcfun_label = None

        #    Final HF energy:             -76.041697549811
        ref_e_scf = -76.041697549811

        ref_dip = 0.7867699  # a.u.

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf, ref_dip)

    def test_scf_dft(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')

        potfile = None

        xcfun_label = 'b3lyp'

        #    Final DFT energy:            -76.443545741524
        ref_e_scf = -76.443545741524

        ref_dip = 0.7312569  # a.u.

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf, ref_dip)

    def test_scf_dft_slda(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')

        potfile = None

        xcfun_label = 'slda'

        #    Final DFT energy:            -76.074208234637
        ref_e_scf = -76.074208234637

        ref_dip = 0.7312887  # a.u.

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf, ref_dip)

    @pytest.mark.skipif('cppe' not in sys.modules, reason='cppe not available')
    def test_scf_hf_pe(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'pe_water.inp')
        potfile = str(here / 'inputs' / 'pe_water.pot')

        xcfun_label = None

        #    Final HF energy:             -76.067159426565
        ref_e_scf = -76.067159426565

        ref_dip = None

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf, ref_dip)

    @pytest.mark.skipif('cppe' not in sys.modules, reason='cppe not available')
    def test_scf_dft_pe(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'pe_water.inp')
        potfile = str(here / 'inputs' / 'pe_water.pot')

        xcfun_label = 'b3lyp'

        #    Final DFT energy:            -76.468733754150
        ref_e_scf = -76.468733754150

        ref_dip = None

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf, ref_dip)


if __name__ == "__main__":
    unittest.main()
