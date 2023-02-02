from mpi4py import MPI
from pathlib import Path
import pytest

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.outputstream import OutputStream
from veloxchem.tpadriver import TpaDriver
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.rsptpa import TPA


@pytest.mark.solvers
class TestTPA:

    def run_scf(self, task):
        scf_settings = {'conv_thresh': 1.0e-8}
        method_settings = {'xcfun': 'BP86'}
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(scf_settings, method_settings)
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)

        return scf_drv.scf_tensors

    def run_tpa(self, inpfile, tpa_type, w, ref_result):

        task = MpiTask([inpfile, None])
        task.input_dict['scf']['checkpoint_file'] = None

        scf_tensors = self.run_scf(task)

        tpa_prop = TPA(
            {
                'damping': task.input_dict['response']['damping'],
                'frequencies': task.input_dict['response']['frequencies'],
                'conv_thresh': '1.0e-8',
                'tpa_type': tpa_type,
            }, {'xcfun': 'BP86'})

        tpa_prop.init_driver(task.mpi_comm, task.ostream)
        tpa_prop.compute(task.molecule, task.ao_basis, scf_tensors)

        if is_mpi_master(task.mpi_comm):
            tpa_result = tpa_prop.rsp_property

            for key in [
                    't4_dict',
                    't3_dict',
                    'NaX3NyNz',
                    'NaA3NxNy',
                    'NaX2Nyz',
                    'NxA2Nyz',
                    'gamma',
            ]:
                if key in tpa_result and key in ref_result:
                    if tpa_result[key] is None:
                        continue
                    assert abs(tpa_result[key][
                        (w, -w, w)].real / ref_result[key].real - 1.0) < 5.0e-5
                    assert abs(tpa_result[key][
                        (w, -w, w)].imag / ref_result[key].imag - 1.0) < 5.0e-5

    def test_tpa_full(self):

        w = 0.05

        ref_result = {
            't4_dict': -39.80683343 - 0.37934801j,
            't3_dict': -92.30911828 - 1.15693267j,
            'NaX3NyNz': -154.10953954 - 1.18508642j,
            'NaA3NxNy': -51.38824383 - 0.10346367j,
            'NaX2Nyz': 621.18409570 + 10.77095784j,
            'NxA2Nyz': 621.79842988 + 2.18029474j,
            'gamma': 905.36879049 + 10.12642181j,
        }

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_tpa.inp')

        self.run_tpa(inpfile, 'full', w, ref_result)

    def test_tpa_reduced(self):

        w = 0.05

        ref_result = {
            't3_dict': -35.37821846 - 0.92780165j,
            'NaX2Nyz': 230.17509800 + 7.10707465j,
            'NxA2Nyz': 230.43689707 + 2.16398432j,
            'gamma': 425.23377660 + 8.34325731j,
        }

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_tpa.inp')

        self.run_tpa(inpfile, 'reduced', w, ref_result)

    def test_update_settings(self):

        tpa_dict = {
            'frequencies': (0.1, 0.12, 0.16, 0.20),
            'damping': 0.01,
            'eri_thresh': 1e-13,
            'qq_type': 'QQ',
            'batch_size': 99,
            'max_iter': 199,
            'conv_thresh': 1e-5,
            'lindep_thresh': 1e-11,
            'restart': False,
            'checkpoint_file': 'mycheckpoint.h5',
            'timing': True,
            'profiling': True,
            'memory_profiling': True,
            'memory_tracing': True,
        }

        tpa_drv = TpaDriver(MPI.COMM_WORLD, OutputStream(None))

        for key, val in tpa_dict.items():
            assert getattr(tpa_drv, key) != val

        tpa_drv.update_settings(tpa_dict)

        for key, val in tpa_dict.items():
            assert getattr(tpa_drv, key) == val
