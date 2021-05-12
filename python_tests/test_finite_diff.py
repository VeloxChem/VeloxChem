from pathlib import Path

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scffirstorderprop import ScfFirstOrderProperties
from veloxchem.rsppolarizability import Polarizability


class TestFiniteDifference:

    def run_mu_and_beta_fd(self, inpfile, xcfun_label):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        scf_prop = ScfFirstOrderProperties(task.mpi_comm, task.ostream)
        scf_prop.compute(task.molecule, task.ao_basis, scf_drv.scf_tensors)

        cart = 2  # z-component
        delta_ef = 1.0e-4
        efield_plus = [0.0, 0.0, delta_ef]
        efield_minus = [0.0, 0.0, -delta_ef]

        method_dict_plus = dict(task.input_dict['method_settings'])
        method_dict_minus = dict(task.input_dict['method_settings'])
        method_dict_plus['electric_field'] = efield_plus
        method_dict_minus['electric_field'] = efield_minus

        scf_drv_plus = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv_plus.update_settings(task.input_dict['scf'], method_dict_plus)
        scf_drv_plus.compute(task.molecule, task.ao_basis, task.min_basis)

        lr_prop_plus = Polarizability({'frequencies': '0'}, method_dict_plus)
        lr_prop_plus.init_driver(task.mpi_comm, task.ostream)
        lr_prop_plus.compute(task.molecule, task.ao_basis,
                             scf_drv_plus.scf_tensors)

        scf_drv_minus = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv_minus.update_settings(task.input_dict['scf'], method_dict_minus)
        scf_drv_minus.compute(task.molecule, task.ao_basis, task.min_basis)

        lr_prop_minus = Polarizability({'frequencies': '0'}, method_dict_minus)
        lr_prop_minus.init_driver(task.mpi_comm, task.ostream)
        lr_prop_minus.compute(task.molecule, task.ao_basis,
                              scf_drv_minus.scf_tensors)

        if is_mpi_master(task.mpi_comm):
            e_scf_plus = scf_drv_plus.get_scf_energy()
            e_scf_minus = scf_drv_minus.get_scf_energy()
            mu_z_fd = -(e_scf_plus - e_scf_minus) / (2.0 * delta_ef)
            mu_z = scf_prop.get_property('dipole moment')[cart]
            assert abs(mu_z_fd - mu_z) < 1.0e-6

            alpha_zz_plus = -lr_prop_plus.get_property('response_functions')[
                ('z', 'z', 0)]
            alpha_zz_minus = -lr_prop_minus.get_property('response_functions')[
                ('z', 'z', 0)]
            beta_zzz_fd = (alpha_zz_plus - alpha_zz_minus) / (2.0 * delta_ef)
            beta_zzz = -4.64532629
            assert abs(beta_zzz_fd - beta_zzz) < 1.0e-3

    def test_mu_and_beta_fd(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')

        xcfun_label = None

        self.run_mu_and_beta_fd(inpfile, xcfun_label)
