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

    def run_shg(self, inpfile, xcfun_label, shg_type, ref_result):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        if xcfun_label is not None:
            task.input_dict['method_settings']['xcfun'] = xcfun_label

        if shg_type is not None:
            task.input_dict['response']['shg_type'] = shg_type

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
                           ref_result[comp].real) < 1.0e-5
                assert abs(shg_result['beta'][freq][ind].imag -
                           ref_result[comp].imag) < 1.0e-5

            assert abs(shg_result['beta_bar'][freq].real -
                       ref_result['beta_bar'].real) < 1.0e-5
            assert abs(shg_result['beta_bar'][freq].imag -
                       ref_result['beta_bar'].imag) < 1.0e-5

    def test_shg_hf(self):

        # vlxtag: RHF, SHG, QR

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'methanol_shg.inp')

        xcfun_label = None
        shg_type = 'full'

        ref_result = {
            'x': -40.52766331 - 11.31636048j,
            'y': -0.11777265 - 0.02104016j,
            'z': 24.95268906 + 7.53241250j,
            'beta_bar': -39.79479881 - 11.62336253j,
        }

        self.run_shg(inpfile, xcfun_label, shg_type, ref_result)

    def test_shg_b3lyp(self):

        # vlxtag: RKS, SHG, QR

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'methanol_shg.inp')

        xcfun_label = 'b3lyp'
        shg_type = 'full'

        ref_result = {
            'x': -48.55320751 - 20.31312232j,
            'y': -0.09597950 - 0.00334644j,
            'z': 32.09728459 + 15.20019960j,
            'beta_bar': -47.13862975 - 21.37978710j,
        }

        self.run_shg(inpfile, xcfun_label, shg_type, ref_result)

    def test_shg_reduced_hf(self):

        # vlxtag: RHF, SHG, QR

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'methanol_shg.inp')

        xcfun_label = None
        shg_type = 'reduced'

        ref_result = {
            'x': -40.53256501 - 11.43101669j,
            'y': -0.11772170 - 0.02074979j,
            'z': 24.95878533 + 7.63353776j,
            'beta_bar': -39.80239953 - 11.76359430j,
        }

        self.run_shg(inpfile, xcfun_label, shg_type, ref_result)

    def test_shg_reduced_b3lyp(self):

        # vlxtag: RKS, SHG, QR

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'methanol_shg.inp')

        xcfun_label = 'b3lyp'
        shg_type = 'reduced'

        ref_result = {
            'x': -48.45584560 - 20.92686505j,
            'y': -0.09586086 - 0.00356429j,
            'z': 32.02926487 + 15.65493314j,
            'beta_bar': -47.04068190 - 22.02152032j,
        }

        self.run_shg(inpfile, xcfun_label, shg_type, ref_result)
