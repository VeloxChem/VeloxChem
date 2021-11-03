from mpi4py import MPI
from pathlib import Path
import unittest
import pytest

from veloxchem.veloxchemlib import is_mpi_master
from veloxchem.outputstream import OutputStream
from veloxchem.shgdriver import SHGDriver
from veloxchem.mpitask import MpiTask
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.rspshg import SHG


class TestSHG(unittest.TestCase):
    """
    Tests the SHG code 
    """

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
        
        # x-component 


        self.assertTrue(abs(shg_result[0.2][0].real - ref_result['x'].real) < 1.0e-6)

        self.assertTrue(abs(shg_result[0.2][0].imag - ref_result['x'].imag) < 1.0e-6)
        
        # y-component 
        
        self.assertTrue(abs(shg_result[0.2][1].real - ref_result['y'].real) < 1.0e-6)

        self.assertTrue(abs(shg_result[0.2][1].imag - ref_result['y'].imag) < 1.0e-6)

        # z-component 
        
        self.assertTrue(abs(shg_result[0.2][2].real - ref_result['z'].real) < 1.0e-6)

        self.assertTrue(abs(shg_result[0.2][2].imag - ref_result['z'].imag) < 1.0e-6)


    def test_shg(self):

        w = 0.2

        ref_result = { 
            'x': 163.69097360 + 177.80276777j,
            'y': 0 + 0j,
            'z': 0 + 0J,
        }

        here = Path(__file__).parent

        inpfile = str(here / 'inputs' / 'water_shg.inp')

        self.run_shg(inpfile, w, ref_result)


if __name__ == '__main__':
    unittest.main()
