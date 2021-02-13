from mpi4py import MPI
import numpy as np
import unittest
from pathlib import Path
try:
    import cppe
except ImportError:
    pass

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.scffirstorderprop import ScfFirstOrderProperties


class TestSCF(unittest.TestCase):

    def run_scf(self, inpfile, potfile, xcfun_label, ref_e_scf, ref_dip):

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None

        if potfile is not None:
            task.input_dict['method_settings']['potfile'] = potfile

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_drv = ScfUnrestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        scf_prop = ScfFirstOrderProperties(task.mpi_comm, task.ostream)
        scf_prop.compute(task.molecule, task.ao_basis, scf_drv.scf_tensors)

        if task.mpi_rank == mpi_master():
            e_scf = scf_drv.get_scf_energy()
            tol = 1.0e-5 if xcfun_label is not None else 1.0e-6
            self.assertTrue(np.max(np.abs(e_scf - ref_e_scf)) < tol)

    def test_scf_dft_blyp(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_unrest.inp')

        potfile = None

        xcfun_label = 'blyp'

        ref_e_scf = -76.119729

        ref_dip = None  # needs to be added later

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf, ref_dip)

    def test_scf_dft_b3lyp(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_unrest.inp')

        potfile = None

        xcfun_label = 'b3lyp'

        ref_e_scf = -76.140226

        ref_dip = None  # needs to be added later

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf, ref_dip)

    def test_scf_dft_slda(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_unrest.inp')

        potfile = None

        xcfun_label = 'slda'

        ref_e_scf = -75.75613

        ref_dip = None  # needs to be added later

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf, ref_dip)

    def test_scf_dft_b88x(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_unrest.inp')

        potfile = None

        xcfun_label = 'B88X'

        ref_e_scf = -75.80179

        ref_dip = None  # needs to be added later

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf, ref_dip)


if __name__ == "__main__":
    unittest.main()
