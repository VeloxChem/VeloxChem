from mpi4py import MPI
import unittest
import os

from veloxchem.veloxchemlib import mpi_master
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.cppsolver import ComplexResponse
from veloxchem.mpitask import MpiTask


class TestCppDriver(unittest.TestCase):

    def test_water_cpp(self):

        # scf
        inpfile = os.path.join('inputs', 'water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)

        task = MpiTask([inpfile, None], MPI.COMM_WORLD)
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        scf_tensors = scf_drv.scf_tensors

        # cpp
        cpp_solver = ComplexResponse(task.mpi_comm, task.ostream)

        cpp_solver.update_settings({
            'frequencies': '0, 0.1',
            'eri_thresh': '1.0e-12',
            'conv_thresh': '1.0e-4',
        })

        results = cpp_solver.compute(task.molecule, task.ao_basis, scf_tensors)

        if task.mpi_rank == mpi_master():

            reference = {
                ('x', 'x', 0.0): complex(7.251356, 0.000000),
                ('x', 'y', 0.0): complex(0.000000, 0.000000),
                ('x', 'z', 0.0): complex(0.000000, 0.000000),
                ('y', 'x', 0.0): complex(0.000000, 0.000000),
                ('y', 'y', 0.0): complex(8.724147, 0.000000),
                ('y', 'z', 0.0): complex(0.000000, 0.000000),
                ('z', 'x', 0.0): complex(0.000000, 0.000000),
                ('z', 'y', 0.0): complex(0.000000, 0.000000),
                ('z', 'z', 0.0): complex(7.880163, 0.000000),
                ('x', 'x', 0.1): complex(7.501754, 0.024595),
                ('x', 'y', 0.1): complex(0.000000, 0.000000),
                ('x', 'z', 0.1): complex(0.000000, 0.000000),
                ('y', 'x', 0.1): complex(0.000000, 0.000000),
                ('y', 'y', 0.1): complex(8.912304, 0.017626),
                ('y', 'z', 0.1): complex(0.000000, 0.000000),
                ('z', 'x', 0.1): complex(0.000000, 0.000000),
                ('z', 'y', 0.1): complex(0.000000, 0.000000),
                ('z', 'z', 0.1): complex(8.084276, 0.019433),
            }

            for key in results['properties']:
                real_val = -results['properties'][key].real
                imag_val = -results['properties'][key].imag
                real_ref = reference[key].real
                imag_ref = reference[key].imag
                self.assertTrue(abs(real_val - real_ref) < 1.0e-5)
                self.assertTrue(abs(imag_val - imag_ref) < 1.0e-5)

        task.finish()


if __name__ == "__main__":
    unittest.main()
