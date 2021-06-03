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
from veloxchem.firstorderprop import FirstOrderProperties


class TestScfUnrestricted(unittest.TestCase):

    def run_scf(self, inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                ref_dip):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        if potfile is not None:
            task.input_dict['method_settings']['potfile'] = potfile
        elif electric_field is not None:
            task.input_dict['scf']['electric_field'] = electric_field

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_drv = ScfUnrestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        scf_prop = FirstOrderProperties(task.mpi_comm, task.ostream)
        total_density = scf_drv.scf_tensors['D_alpha'] + scf_drv.scf_tensors['D_beta']
        scf_prop.compute(task.molecule, task.ao_basis, total_density)

        if is_mpi_master(task.mpi_comm):
            e_scf = scf_drv.get_scf_energy()
            tol = 1.0e-5 if xcfun_label is not None else 1.0e-6
            self.assertTrue(np.max(np.abs(e_scf - ref_e_scf)) < tol)

            if ref_dip is not None:
                dip = scf_prop.get_property('dipole moment')
                self.assertTrue(np.max(np.abs(dip - ref_dip)) < 1.0e-5)

    def test_scf_dft_blyp(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_unrest.inp')
        potfile = None

        xcfun_label = 'blyp'
        electric_field = None

        ref_e_scf = -76.119737
        ref_dip = np.array([0.000000, 0.000000, -0.238897])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

    def test_scf_dft_blyp_efield(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_unrest.inp')
        potfile = None

        xcfun_label = 'blyp'
        electric_field = '0, -0.002, 0'

        ref_e_scf = -76.119843
        ref_dip = np.array([0.000000, -0.106245, -0.238934])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

    def test_scf_dft_b3lyp(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_unrest.inp')
        potfile = None

        xcfun_label = 'b3lyp'
        electric_field = None

        ref_e_scf = -76.140232
        ref_dip = np.array([0.000000, 0.000000, -0.229996])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

    def test_scf_dft_slda(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_unrest.inp')
        potfile = None

        xcfun_label = 'slda'
        electric_field = None

        ref_e_scf = -75.756138
        ref_dip = np.array([0.000000, 0.000000, -0.226886])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

    def test_scf_dft_b88x(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_unrest.inp')
        potfile = None

        xcfun_label = 'b88x'
        electric_field = None

        ref_e_scf = -75.801799
        ref_dip = np.array([0.000000, 0.000000, -0.249728])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)


if __name__ == "__main__":
    unittest.main()
