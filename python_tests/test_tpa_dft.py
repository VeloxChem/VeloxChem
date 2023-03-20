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

            for key in ref_result:
                assert abs(tpa_result[key][
                    (w, -w, w)].real / ref_result[key].real - 1.0) < 1.0e-6
                assert abs(tpa_result[key][
                    (w, -w, w)].imag / ref_result[key].imag - 1.0) < 1.0e-6

    def test_tpa_full(self):

        w = 0.05

        ref_result = {
            't4_dict': -39.81102482 - 0.37943318j,
            't3_dict': -92.30945471 - 1.15690967j,
            'NaX3NyNz': -154.11176435 - 1.18509718j,
            'NaA3NxNy': -51.38898554 - 0.10346488j,
            'NaX2Nyz': 621.17712859 + 10.77037598j,
            'NxA2Nyz': 621.79142488 + 2.18017301j,
            'gamma': 905.34732406 + 10.12564409j,
        }

        here = Path(__file__).parent
        inpfile = str(here / 'inputs' / 'water_tpa.inp')

        self.run_tpa(inpfile, 'full', w, ref_result)

    def test_tpa_reduced(self):

        w = 0.05

        ref_result = {
            't3_dict': -35.37827126 - 0.92778160j,
            'NaX2Nyz': 230.17180603 + 7.10669620j,
            'NxA2Nyz': 230.43358578 + 2.16388980j,
            'gamma': 425.22712055 + 8.34280441j,
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
