from pathlib import Path
import numpy as np
import unittest
try:
    import cppe
except ImportError:
    pass

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.scffirstorderprop import ScfFirstOrderProperties


class TestTripletH2(unittest.TestCase):

    def run_scf(self, inpfile, potfile, xcfun_label, ref_e_scf, ref_dip):

        task = MpiTask([inpfile, None])
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

        if is_mpi_master(task.mpi_comm):
            e_scf = scf_drv.get_scf_energy()
            tol = 1.0e-5 if xcfun_label is not None else 1.0e-6
            self.assertTrue(np.max(np.abs(e_scf - ref_e_scf)) < tol)

            dip = scf_prop.get_property('dipole moment')
            self.assertTrue(np.max(np.abs(dip - ref_dip)) < 1.0e-5)

    def test_uscf_hf(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'h2_t.inp')

        potfile = None

        xcfun_label = None

        ref_e_scf = -0.77331601

        ref_dip = np.array([0.000000, 0.000000, 0.000000])

        self.run_scf(inpfile, potfile, xcfun_label, ref_e_scf, ref_dip)


if __name__ == "__main__":
    unittest.main()
