from mpi4py import MPI
import unittest
import os

from veloxchem.veloxchemlib import mpi_master
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.tpafulldriver import TpaFullDriver
from veloxchem.tpareddriver import TpaReducedDriver


class TestTPA(unittest.TestCase):

    def run_tpa(self, inpfile, tpa_type, w, ref_result):

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        task.input_dict['scf']['checkpoint_file'] = None

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.update_settings(task.input_dict['scf'],
                                task.input_dict['method_settings'])
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        scf_tensors = scf_drv.scf_tensors

        if tpa_type.lower() == 'full':
            tpa_drv = TpaFullDriver(task.mpi_comm, task.ostream)
        elif tpa_type.lower() == 'reduced':
            tpa_drv = TpaReducedDriver(task.mpi_comm, task.ostream)

        tpa_drv.update_settings({
            'damping': task.input_dict['response']['damping'],
            'frequencies': task.input_dict['response']['frequencies'],
            'conv_thresh': '1.0e-8',
        })

        tpa_result = tpa_drv.compute(task.molecule, task.ao_basis, scf_tensors)

        if task.mpi_rank == mpi_master():
            for key in [
                    't4_dict', 't3_dict', 'NaX3NyNz', 'NaA3NxNy', 'NaX2Nyz',
                    'NxA2Nyz', 'gamma'
            ]:
                if key in tpa_result and key in ref_result:
                    if tpa_result[key] is None:
                        continue
                    self.assertTrue(
                        abs(tpa_result[key][(w, -w, w)].real -
                            ref_result[key].real) < 1.0e-4)
                    self.assertTrue(
                        abs(tpa_result[key][(w, -w, w)].imag -
                            ref_result[key].imag) < 1.0e-4)

    def test_tpa_full(self):

        w = 0.05

        ref_result = {
            't4_dict': 11.43071305 + 0.04957732j,
            't3_dict': -42.19841751 - 0.28695214j,
            'NaX3NyNz': -81.62345190 - 0.35812832j,
            'NaA3NxNy': -27.21320341 - 0.03029788j,
            'NaX2Nyz': 270.69041328 + 2.67837597j,
            'NxA2Nyz': 270.83461366 + 0.52758094j,
            'gamma': 401.92066716 + 2.58015589j,
        }

        inpfile = os.path.join('inputs', 'water_tpa.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        self.run_tpa(inpfile, 'full', w, ref_result)

    def test_tpa_reduced(self):

        w = 0.05

        ref_result = {
            't3_dict': -15.12982062 - 0.19793495j,
            'NaX2Nyz': 96.30910639 + 1.72679037j,
            'NxA2Nyz': 96.36431088 + 0.51886895j,
            'gamma': 177.54359664 + 2.04772438j,
        }

        inpfile = os.path.join('inputs', 'water_tpa.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        self.run_tpa(inpfile, 'reduced', w, ref_result)


if __name__ == "__main__":
    unittest.main()
