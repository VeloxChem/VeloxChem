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

    def run_shg(self, inpfile, xcfun_label, ref_result):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        scf_tensors = self.run_scf(task)

        shg_prop = SHG(task.input_dict['response'],
                       task.input_dict['method_settings'])
        shg_prop.init_driver(task.mpi_comm, task.ostream)
        shg_prop.compute(task.molecule, task.ao_basis, scf_tensors)
        shg_result = shg_prop.rsp_property

        if task.mpi_rank == mpi_master():
            freq = list(shg_result['beta'].keys())[0]

            for ind, comp in enumerate('xyz'):
                assert abs(shg_result['beta'][freq][ind].real -
                           ref_result[comp].real) < 1.0e-6
                assert abs(shg_result['beta'][freq][ind].imag -
                           ref_result[comp].imag) < 1.0e-6

            assert abs(shg_result['beta_bar'][freq].real -
                       ref_result['beta_bar'].real) < 1.0e-6
            assert abs(shg_result['beta_bar'][freq].imag -
                       ref_result['beta_bar'].imag) < 1.0e-6

    def test_shg_hf(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'methanol_shg.inp')

        xcfun_label = None

        ref_result = {
            'x': -40.52766331 - 11.31636048j,
            'y': -0.11777265 - 0.02104016j,
            'z': 24.95268906 + 7.53241250j,
            'beta_bar': -39.79479881 - 11.62336253j,
        }

        self.run_shg(inpfile, xcfun_label, ref_result)

    def test_shg_b3lyp(self):

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'methanol_shg.inp')

        xcfun_label = 'b3lyp'

        ref_result = {
            'x': -48.55306041 - 20.31300835j,
            'y': -0.09594022 - 0.00331360j,
            'z': 32.09729558 + 15.20016566j,
            'beta_bar': -47.13851686 - 21.37968687j,
        }

        self.run_shg(inpfile, xcfun_label, ref_result)
