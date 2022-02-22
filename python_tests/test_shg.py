from pathlib import Path
import pytest

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.rspshg import SHG


@pytest.mark.solvers
class TestSHG:

    def run_scf(self, task):

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        return scf_drv.scf_tensors

    def run_shg(self, inpfile, w, ref_result):

        task = MpiTask([inpfile, None])

        task.input_dict['scf']['checkpoint_file'] = None

        scf_tensors = self.run_scf(task)

        shg_prop = SHG({
            'damping': task.input_dict['response']['damping'],
            'frequencies': task.input_dict['response']['frequencies'],
            'conv_thresh': '1.0e-8',
        })

        shg_prop.init_driver(task.mpi_comm, task.ostream)
        shg_prop.compute(task.molecule, task.ao_basis, scf_tensors)
        shg_result = shg_prop.rsp_property

        if task.mpi_rank == mpi_master():
            freq = list(shg_result['beta'].keys())[0]

            for ind, comp in enumerate('xyz'):
                assert abs(shg_result['beta'][freq][ind].real -
                           ref_result[comp].real) < 1.0e-5
                assert abs(shg_result['beta'][freq][ind].imag -
                           ref_result[comp].imag) < 1.0e-5

            assert abs(shg_result['beta_bar'][freq].real -
                       ref_result['beta_bar'].real) < 1.0e-4
            assert abs(shg_result['beta_bar'][freq].imag -
                       ref_result['beta_bar'].imag) < 1.0e-4

    def test_shg(self):

        w = 0.1

        ref_result = {
            'x': -35.83087082 - 19.02338541j,
            'y': -0.13248044 - 0.04458826j,
            'z': 20.54249935 + 12.04788926j,
            'beta_bar': -12.94957848 - 7.99384834j,
        }

        here = Path(__file__).parent

        inpfile = str(here / 'inputs' / 'methanol_shg.inp')

        self.run_shg(inpfile, w, ref_result)
