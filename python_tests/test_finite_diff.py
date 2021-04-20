from pathlib import Path

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.scffirstorderprop import ScfFirstOrderProperties


class TestFiniteDifference:

    def run_dipole_fd(self, inpfile, xcfun_label):

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

        for cart in range(3):
            delta_ef = 1.0e-4

            efield_plus = [0.0] * 3
            efield_plus[cart] = delta_ef
            scf_dict_plus = dict(task.input_dict['scf'])
            scf_dict_plus['electric_field'] = ','.join(
                [str(x) for x in efield_plus])

            scf_drv_plus = ScfRestrictedDriver(task.mpi_comm, task.ostream)
            scf_drv_plus.update_settings(scf_dict_plus,
                                         task.input_dict['method_settings'])
            scf_drv_plus.compute(task.molecule, task.ao_basis, task.min_basis)

            efield_minus = [0.0] * 3
            efield_minus[cart] = -delta_ef
            scf_dict_minus = dict(task.input_dict['scf'])
            scf_dict_minus['electric_field'] = ','.join(
                [str(x) for x in efield_minus])

            scf_drv_minus = ScfRestrictedDriver(task.mpi_comm, task.ostream)
            scf_drv_minus.update_settings(scf_dict_minus,
                                          task.input_dict['method_settings'])
            scf_drv_minus.compute(task.molecule, task.ao_basis, task.min_basis)

            if is_mpi_master(task.mpi_comm):
                e_scf_plus = scf_drv_plus.get_scf_energy()
                e_scf_minus = scf_drv_minus.get_scf_energy()
                dipole_comp_fd = -1.0 * (e_scf_plus - e_scf_minus) / (2.0 *
                                                                      delta_ef)
                dipole_comp = scf_prop.get_property('dipole moment')[cart]
                assert abs(dipole_comp_fd - dipole_comp) < 1.0e-6

    def test_dipole_fd(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water.inp')

        xcfun_label = None

        self.run_dipole_fd(inpfile, xcfun_label)
