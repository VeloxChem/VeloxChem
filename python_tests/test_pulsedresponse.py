from mpi4py import MPI
import numpy as np
import os
import unittest
import h5py
import pickle

from veloxchem.veloxchemlib import mpi_master
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.lrsolver import LinearResponseSolver
from veloxchem.crsp import ComplexResponse
from veloxchem.mpitask import MpiTask
from veloxchem.pulsedrsp import PulsedResponse
from veloxchem.rspabsorption import Absorption


class TestComplexResponse(unittest.TestCase):


    def test_pulsed_filesave(self):
        print("Testing the h5 output file")
        h5fname = "test_h5pulse"

        inpfile = os.path.join('inputs', 'pulsed_water.inp')
        outfile = inpfile.replace('.inp', '.out')

        task = MpiTask([inpfile, outfile], MPI.COMM_WORLD)

        # scf
        pulse_input = task.input_dict['pulses']
        pulse_input['h5'] = h5fname
        crsp_input = {} 

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        scf_tensors = scf_drv.scf_tensors


        # Run the computaiton 
        pulsed_response = PulsedResponse(task.mpi_comm, task.ostream)
        pulsed_response.update_settings(pulse_input, crsp_input)
        results = pulsed_response.compute(task.molecule, task.ao_basis, scf_tensors)

        directions = ['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz']
        expected_keys = ['amplitudes', 'frequencies','zero_padded_frequencies','zero_padded'] + directions


        # Test the contents of the file
        print("Testing the presence of datasets in the h5 output file")


        with h5py.File('{}.h5'.format(h5fname), 'r') as hf:
            for key in expected_keys:
                if key not in hf.keys():
                    self.fail("Error - expected key '{}' but did not find it".format(key))

            primary_key = "frequencies"
            if hf.get('zero_padded')[()]:
                primary_key = "zero_padded_frequencies"

            for key in directions:
              try:
                assert(len(hf.get(primary_key)[()]) == len(hf.get(key)[()]))
              except:
                self.fail("Len of {}[{}] did not match data length {}!= {}".format(primary_key, key, 
                            len(hf.get('zero_padded_frequencies')[()],
                            len(hf.get(key)[()]))))

        task.finish()

    def test_pulsed_response(self):

        inpfile = os.path.join('inputs', 'pulsed_water.inp')
        outfile = inpfile.replace('.inp', '.out')

        task = MpiTask([inpfile, outfile], MPI.COMM_WORLD)

        # scf
        pulse_input = task.input_dict['pulses']
        crsp_input = {} 

        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)

        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        scf_tensors = scf_drv.scf_tensors


        pulsed_response = PulsedResponse(task.mpi_comm, task.ostream)
        pulsed_response.update_settings(pulse_input, crsp_input)

        results = pulsed_response.compute(task.molecule, task.ao_basis, scf_tensors)


        for key, item in results['properties'].items():
            print(key, item)

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


        task.finish()


if __name__ == "__main__":
    unittest.main()
