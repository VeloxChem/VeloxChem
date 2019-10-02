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
    1.59140758e-09 - 5.77836923e-11j, -3.67213109e-09 + 6.83655709e-09j,
    -1.74502523e-08 - 2.85463969e-08j, 1.27592033e-07 - 2.48570048e-09j,
    -2.10126192e-07 + 3.75898547e-07j, -6.52015265e-07 - 1.10806061e-06j,
    3.39573025e-06 - 9.04784541e-09j, -3.98764543e-06 + 6.86014049e-06j,
    -8.08028193e-06 - 1.42758096e-05j, 2.99992466e-05 + 4.24524335e-07j,
    -2.50990214e-05 + 4.15548207e-05j, -3.32107960e-05 - 6.10472590e-05j,
    8.79741389e-05 + 2.72498271e-06j, -5.23995798e-05 + 8.35473246e-05j,
    -4.52672827e-05 - 8.66488453e-05j, 8.56381194e-05 + 4.09487227e-06j,
    -3.62872478e-05 + 5.57523153e-05j, -2.04600590e-05 - 4.08219182e-05j,
    2.76723777e-05 + 1.78996808e-06j, -8.33600161e-06 + 1.23483456e-05j,
    -3.06625102e-06 - 6.38352735e-06j, 2.96819791e-06 + 2.42173236e-07j,
    -6.35271703e-07 + 9.07749902e-07j, -1.52350347e-07 - 3.31334943e-07j
])

frequencies = np.array([
    0.231, 0.238, 0.245, 0.252, 0.259, 0.266, 0.273, 0.28, 0.287, 0.294, 0.301,
    0.308, 0.315, 0.322, 0.329, 0.336, 0.343, 0.35, 0.357, 0.364, 0.371, 0.378,
    0.385, 0.392
])

zero_padded_frequencies = np.array([
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
    0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, 0. + 0.j, -9.38450147 - 0.15495907j,
    -9.63961388 - 0.17972249j, -9.93719219 - 0.21126062j,
    -10.28938819 - 0.25240444j, -10.71365508 - 0.30763993j,
    -11.23598443 - 0.38441952j, -11.89682064 - 0.49586375j,
    -12.76260525 - 0.66682352j, -13.95020958 - 0.94905214j,
    -15.68405007 - 1.46507308j, -18.44600447 - 2.56311645j,
    -23.39852404 - 5.55077652j, -32.37389248 - 17.93780658j,
    -0.3854528 - 52.27765737j, 16.04673524 - 14.08693823j,
    7.65070008 - 4.76397713j, 3.19292747 - 2.31623591j, 0.6336971 - 1.37180542j,
    -1.02092542 - 0.91721761j, -2.19004003 - 0.66663608j,
    -3.07350937 - 0.51583667j, -3.7778817 - 0.41984064j,
    -4.36546063 - 0.35690017j, -4.87587013 - 0.3157651j
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

        crsp_input = {}

        # Setup the Task parameters
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
        scf_drv.compute(task.molecule, task.ao_basis, task.min_basis)
        scf_tensors = scf_drv.scf_tensors

        # Run the pulsed response computation
        pulsed_response = PulsedResponse(task.mpi_comm, task.ostream)
        pulsed_response.update_settings(pulse_input, crsp_input)
        cls.results = pulsed_response.compute(task.molecule, task.ao_basis,
                                              scf_tensors)

        task.finish()

    def test_filesave(self):

        directions = ['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz']
        expected_keys = [
            'amplitudes', 'frequencies', 'zero_padded_frequencies',
            'zero_padded'
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
            if hf.get('zero_padded')[()]:
                primary_key = "zero_padded_frequencies"
                self.assertTrue(
                    len(hf.get('zero_padded_frequencies')[(
                    )]) == len(hf.get('zero_padded_amplitudes')[()]))

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
            self.assertTrue(
                np.allclose(
                    hf.get('zero_padded_frequencies')[()],
                    zero_padded_frequencies))
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

            # Verify that the results stored are the same as expected
            self.assertTrue(
                np.allclose(self.results['pulse_settings']['freq_amplitude'],
                            amplitudes))
            self.assertTrue(
                np.allclose(self.results['pulse_settings']['frequencies'],
                            frequencies))


if __name__ == "__main__":
    unittest.main()
