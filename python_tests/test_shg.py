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

            # x-component

            assert abs(shg_result[0.1][0].real - ref_result['x'].real) < 1.0e-5
            assert abs(shg_result[0.1][0].imag - ref_result['x'].imag) < 1.0e-5

            # y-component

            assert abs(shg_result[0.1][1].real - ref_result['y'].real) < 1.0e-5
            assert abs(shg_result[0.1][1].imag - ref_result['y'].imag) < 1.0e-5

            # z-component

            assert abs(shg_result[0.1][2].real - ref_result['z'].real) < 1.0e-5
            assert abs(shg_result[0.1][2].imag - ref_result['z'].imag) < 1.0e-5

    def test_shg(self):

        w = 0.1

        ref_result = {
            'x': -35.83087082066745-19.023385411711395j,
            'y': -0.1324804367386819-0.0445882565520379j,
            'z': 20.54249934914053+12.047889262521506j,
        }

        here = Path(__file__).parent

        inpfile = str(here / 'inputs' / 'methanol_shg.inp')

        self.run_shg(inpfile, w, ref_result)
