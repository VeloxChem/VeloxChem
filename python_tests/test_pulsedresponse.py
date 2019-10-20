from mpi4py import MPI
import numpy as np
import os
import unittest
import h5py

from veloxchem.veloxchemlib import mpi_master
from veloxchem.scfrestdriver import ScfRestrictedDriver
from veloxchem.mpitask import MpiTask
from veloxchem.pulsedrsp import PulsedResponse

# ---------------------------------------------------------------------------
# Test data for comparison of Pulse Response Computed results
# ---------------------------------------------------------------------------
amplitudes = np.array([
    0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j,
    0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j,
    0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j,
    0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j,
    0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j,
    1.56529426e-09 + 2.92866741e-10j, -5.08292880e-09 + 5.86403396e-09j,
    -1.07602222e-08 - 3.16800521e-08j, 1.25027182e-07 + 2.55755615e-08j,
    -2.87496573e-07 + 3.20621981e-07j, -3.92951818e-07 - 1.22413688e-06j,
    3.31493710e-06 + 7.36381572e-07j, -5.39592725e-06 + 5.81780137e-06j,
    -4.75041339e-06 - 1.57010595e-05j, 2.91747820e-05 + 6.99765027e-06j,
    -3.36065744e-05 + 3.50337288e-05j, -1.90041005e-05 - 6.68473558e-05j,
    8.52315566e-05 + 2.19648904e-05j, -6.94570698e-05 + 7.00113338e-05j,
    -2.51482949e-05 - 9.44706968e-05j, 8.26518543e-05 + 2.27887354e-05j,
    -4.76377685e-05 + 4.64298183e-05j, -1.10027375e-05 - 4.43168455e-05j,
    2.66049829e-05 + 7.81916594e-06j, -1.08426945e-05 + 1.02179518e-05j,
    -1.59060922e-06 - 6.90081727e-06j, 2.84269521e-06 + 8.87654593e-07j,
    -8.18995446e-07 + 7.46208068e-07j, -7.59234571e-08 - 3.56691886e-07j
])

frequencies = np.array([
    0., 0.007, 0.014, 0.021, 0.028, 0.035, 0.042, 0.049, 0.056, 0.063, 0.07,
    0.077, 0.084, 0.091, 0.098, 0.105, 0.112, 0.119, 0.126, 0.133, 0.14, 0.147,
    0.154, 0.161, 0.168, 0.175, 0.182, 0.189, 0.196, 0.203, 0.21, 0.217, 0.224,
    0.231, 0.238, 0.245, 0.252, 0.259, 0.266, 0.273, 0.28, 0.287, 0.294, 0.301,
    0.308, 0.315, 0.322, 0.329, 0.336, 0.343, 0.35, 0.357, 0.364, 0.371, 0.378,
    0.385, 0.392
])

xx = np.array([
    0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j,
    0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j,
    0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j,
    0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j,
    0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 9.38450147 + 0.15495907j,
    9.63961388 + 0.17972249j, 9.93719219 + 0.21126062j,
    10.28938819 + 0.25240444j, 10.71365508 + 0.30763993j,
    11.23598443 + 0.38441952j, 11.89682064 + 0.49586375j,
    12.76260525 + 0.66682352j, 13.95020958 + 0.94905214j,
    15.68405007 + 1.46507308j, 18.44600447 + 2.56311645j,
    23.39852404 + 5.55077652j, 32.37389248 + 17.93780658j,
    0.38545281 + 52.27765737j, -16.04673524 + 14.08693823j,
    -7.65070008 + 4.76397713j, -3.19292747 + 2.31623591j,
    -0.6336971 + 1.37180542j, 1.02092542 + 0.91721761j,
    2.19004003 + 0.66663608j, 3.07350937 + 0.51583667j, 3.7778817 + 0.41984064j,
    4.36546063 + 0.35690017j, 4.87587013 + 0.3157651j
])
# ---------------------------------------------------------------------------


class TestComplexResponse(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """ setUpClass - initiates the Test object making a single calculcation for all
            subsequent tests
        """

        cls.h5fname = "pulsed"

        inpfile = os.path.join('inputs', 'pulsed_water.inp')
        if not os.path.isfile(inpfile):
            inpfile = os.path.join('python_tests', inpfile)
        outfile = inpfile.replace('.inp', '.out')

        task = MpiTask([inpfile, outfile], MPI.COMM_WORLD)

        # scf
        pulse_input = task.input_dict['pulses']
        pulse_input['h5'] = cls.h5fname

        # Fix the test parameters to match the test comparison results
        pulse_input['field_max'] = 1.0e-5
        pulse_input['number_pulses'] = 1
        pulse_input['centers'] = 300
        pulse_input['field_cutoff_ratio'] = 1e-5
        pulse_input['frequency_range'] = "0.2-0.4(0.007)"
        pulse_input['carrier_frequencies'] = 0.325
        pulse_input['pulse_widths'] = 50
        pulse_input['pol_dir'] = "xyz"

        cpp_input = {}

        # Setup the Task parameters
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        scf_tensors = scf_drv.scf_tensors

        # Run the pulsed response computation
        pulsed_response = PulsedResponse(task.mpi_comm, task.ostream)
        pulsed_response.update_settings(pulse_input, cpp_input)
        cls.results = pulsed_response.compute(task.molecule, task.ao_basis,
                                              scf_tensors)

        task.finish()

    def test_center_freq_peak(self):

        def find_nearest(array, value):
            err = 9999999999
            best_val = 0
            for val in array:
                if abs(val - value) < err:
                    err = abs(val - value)
                    best_val = val

            return best_val

        if MPI.COMM_WORLD.Get_rank() == mpi_master():

            hf = h5py.File('{}.h5'.format(self.h5fname), 'r')
            calc_amplitudes = hf.get('amplitudes')[()]
            calc_frequencies = hf.get('frequencies')[()]

            max_index = np.argmax(np.abs(calc_amplitudes))
            cloest_freq = find_nearest(calc_frequencies, 0.325)

            if not calc_frequencies[max_index] == cloest_freq:
                self.fail(
                    "Amplitude max did not match center freq: {} != {}".format(
                        cloest_freq, 0.325))

    def test_filesave(self):

        directions = ['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz']
        expected_keys = ['amplitudes', 'frequencies', 'zero_padded'
                         ] + directions

        # Test the contents of the file

        if MPI.COMM_WORLD.Get_rank() == mpi_master():

            hf = h5py.File('{}.h5'.format(self.h5fname), 'r')
            for key in expected_keys:
                if key not in hf.keys():
                    self.fail(
                        "Error - expected key '{}' but did not find it".format(
                            key))

            # Verify that the non-zero padded list of amplitudes and
            # frequencies match
            self.assertTrue(
                len(hf.get('amplitudes')[()]) == len(hf.get('frequencies')[()]))

            primary_key = "frequencies"

            for key in directions:
                try:
                    # Match the length of the frequency list
                    # (zero padded or non zero padded)
                    # and the polarization direction lists of amplitudes
                    self.assertTrue(
                        len(hf.get(primary_key)[()]) == len(hf.get(key)[()]))
                except ValueError:
                    self.fail("Len of {}[{}] did not match data length {}!= {}".
                              format(
                                  primary_key, key,
                                  len(
                                      hf.get('zero_padded_frequencies')[()],
                                      len(hf.get(key)[()]))))

            # Verify that the results stored are the same as expected
            self.assertTrue(np.allclose(hf.get('amplitudes')[()], amplitudes))
            self.assertTrue(np.allclose(hf.get('frequencies')[()], frequencies))
            self.assertTrue(np.allclose(hf.get('xx')[()], xx))

            hf.close()

    def test_pulsed_response(self):

        expected_keys = [
            'properties', 'solutions', 'kappas', 'pulse_settings',
            'properties_zeropad'
        ]

        if MPI.COMM_WORLD.Get_rank() == mpi_master():

            for key in expected_keys:
                if key not in self.results:
                    self.fail(
                        "Error - expected key '{}' but did not find it".format(
                            key))

            # Extract the truncated part of the data
            trunc_amplitudes = amplitudes[33:]
            trunc_frequencies = frequencies[33:]
            # Verify that the results stored are the same as expected
            self.assertTrue(
                np.allclose(self.results['pulse_settings']['freq_amplitude'],
                            trunc_amplitudes))
            self.assertTrue(
                np.allclose(self.results['pulse_settings']['frequencies'],
                            trunc_frequencies))


if __name__ == "__main__":
    unittest.main()
