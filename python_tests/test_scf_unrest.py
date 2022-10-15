from pathlib import Path
import numpy as np

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfunrestdriver import ScfUnrestrictedDriver
from veloxchem.firstorderprop import FirstOrderProperties


class TestScfUnrestricted:

    def run_scf(self, inpfile, potfile, xcfun_label, efield, ref_e_scf,
                ref_dip):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        if potfile is not None:
            task.input_dict['method_settings']['potfile'] = potfile
        elif efield is not None:
            task.input_dict['method_settings']['electric_field'] = efield

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_drv = ScfUnrestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        scf_prop = FirstOrderProperties(task.mpi_comm, task.ostream)
        scf_prop.compute_scf_prop(task.molecule, task.ao_basis,
                                  scf_drv.scf_tensors)

        if is_mpi_master(task.mpi_comm):
            e_scf = scf_drv.get_scf_energy()
            tol = 1.0e-5 if xcfun_label is not None else 1.0e-6
            assert np.max(np.abs(e_scf - ref_e_scf)) < tol

            if ref_dip is not None:
                dip = scf_prop.get_property('dipole moment')
                assert np.max(np.abs(dip - ref_dip)) < 1.0e-5

    def test_scf_dft_blyp(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_triplet.inp')
        potfile = None

        xcfun_label = 'blyp'
        electric_field = None

        ref_e_scf = -76.119737
        ref_dip = np.array([0.000000, 0.000000, -0.238897])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

    def test_scf_dft_blyp_efield(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_triplet.inp')
        potfile = None

        xcfun_label = 'blyp'
        electric_field = '0, -0.002, 0'

        ref_e_scf = -76.119843
        ref_dip = np.array([0.000000, -0.106245, -0.238934])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

    def test_scf_dft_b3lyp(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_triplet.inp')
        potfile = None

        xcfun_label = 'b3lyp'
        electric_field = None

        ref_e_scf = -76.140232
        ref_dip = np.array([0.000000, 0.000000, -0.229996])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

    def test_scf_dft_slda(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_triplet.inp')
        potfile = None

        xcfun_label = 'slda'
        electric_field = None

        ref_e_scf = -75.7561378
        ref_dip = np.array([0.000000, 0.000000, -0.226850])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)

    def test_scf_dft_b88x(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_triplet.inp')
        potfile = None

        xcfun_label = 'b88x'
        electric_field = None

        ref_e_scf = -75.801799
        ref_dip = np.array([0.000000, 0.000000, -0.249728])

        self.run_scf(inpfile, potfile, xcfun_label, electric_field, ref_e_scf,
                     ref_dip)
