from mpi4py import MPI
import numpy as np
import os
import unittest
import pickle

from veloxchem.veloxchemlib import mpi_master
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.crsp import ComplexResponse
from veloxchem.mpitask import MpiTask
from veloxchem.pulsedrsp import PulsedResponse
from veloxchem.rspabsorption import Absorption


class TestComplexResponse(unittest.TestCase):

#    def test_pulsed_response(self):
#
#        print('Running pulsed response test on water')
#
#        # scf
#        task = MpiTask(["inputs/water.inp", None], MPI.COMM_WORLD)
#        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
#
#        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
#        scf_tensors = scf_drv.scf_tensors
#
#        rsp_dict = {'envelope': 'gaussian',
#                    'smallest_field_ratio': 1e-5,
#                    'file_prefix' : 'water',
#                    'frequencies' : "0.2-0.4(0.007),",
#                    }
#
#        crsp_dict = {
#
#                    }
#
#        pulsed_response = PulsedResponse(task.mpi_comm, task.ostream,  rsp_dict, crsp_dict)
#
#        results = pulsed_response.compute(task.molecule, task.ao_basis, scf_tensors)
#
#        for key, item in results['properties'].items():
#            print(key, item)
#
#        print('Now printing results')
#        if 'properties_zeropad' in results:
#            print("Printing zero padded properties")
#            for key, item in results['properties_zeropad'].items():
#                print(key, item)
#
#        
#        with open('prt_water_function.pickle', 'wb') as f:
#            # Pickle the 'data' dictionary using the highest protocol available.
#            pickle.dump(results, f, protocol=2)
#
#        task.finish()

    def test_noradrenalin(self):

        print('Running pulsed response test on noradrenalin')

        inputfile = os.path.expanduser('inputs/noradrenalin.inp')

        # scf
        task = MpiTask([inputfile, None], MPI.COMM_WORLD)
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        scf_tensors = scf_drv.scf_tensors

        absorp = Absorption({'property' : 'absorption'})
        absorp.init_driver(task.mpi_comm, task.ostream)

        absorp.compute(task.molecule, task.ao_basis, scf_tensors)
        print("abs_results")
        print(absorp.rsp_property)

        
        # Prepare the input to the pulsed response calculation
        rsp_dict = {'envelope': 'gaussian',
                    'frequencies' : "0.0-1.0(0.0018),",
                    'file_prefix' : 'noradrenalin_',
                    'carrier_frequency' : absorp.rsp_property['eigenvalues'][0],
                    'smallest_field_ratio': 1e-5}

        crsp_dict = {

                    }

        pulsed_response = PulsedResponse(task.mpi_comm, task.ostream,  rsp_dict, crsp_dict)
        results = pulsed_response.compute(task.molecule, task.ao_basis, scf_tensors)

        print('Printing properties')
        for key, item in results['properties'].items():
            print(key, item)

        if 'properties_zeropad' in results:
            print("Printing zero padded properties")
            for key, item in results['properties_zeropad'].items():
                print(key, item)

        
        results['rsp_results'] = absorp.rsp_property
        with open('prt_noradrenalin_function.pickle', 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(results, f, protocol=2)

        task.finish()

if __name__ == "__main__":
    unittest.main()
