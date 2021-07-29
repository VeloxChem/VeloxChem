#
#                           VELOXCHEM 1.0-RC2
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

from mpi4py import MPI
import numpy as np
import h5py
import sys

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .cppsolver import ComplexResponse
from .inputparser import parse_seq_range


class PulsedResponse:
    """
    Pulsed Reponse class for computing molecular responses to
    a resonant or non-resonant laser pulse.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variable
        - comm: The MPI communicator.
        - rank: The MPI rank.
        - nodes: The number of MPI processes.
        - ostream: The output stream.
        - field_cutoff_ratio: Float - Ratio between the smallest field
          amplitude to be included in the calculation wrt. the maximum field
          amplitude.  Frequencies with associated field amplitudes smaller than
          this ratio will NOT be included in the calculation.
        - envelope: String - Envelope of the of pulse - available arguments:
          ['gaussian', ]
        - number_pulses: Integer - The number of different pulse settings to
          expect 'N'.  Currently limited to one.
        - pulse_widths: List of floats (len(N)) - pulse widths in [a.u.].
        - carrier_frequencies: List of floats (len(N)) - carrier frequencies in
          [a.u.].
        - field_max: List of floats (len(N)) - pulse amplitudes in [a.u.].
        - centers: List of floats (len(N)) - time centers for the pulses in
          [a.u.].
        - pol_dir: string - polarization directions, arguments given in
          combinations of x, y and z - e.g.: 'x' for [1 0 0], yz for [0 1 1].
        - frequency_range: List of frequencies to map solution to.  Given as
          range: 'start-end(df)' in [a.u.], e.g. 0.2-0.4(0.007) for 0.2 -> 0.4 a.u.
          in steps of 0.007 a.u.  Zero-padding will be applied if range does not
          start at 0.
        - CEP: List of floats (len(N)) - carrier envelope phases in [radians].
        - h5file: String - optional - name of requested h5 formatted result
          file.
        - ascii: String - optional - name of requested ASCII formatted file.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes pulsed response.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            ostream = OutputStream(sys.stdout)

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream

        self.zero_pad = False
        self.pulse_settings = {}
        self.cpp_settings = {}

        # Default Pulse settings
        # Smallest fraction field amplitude to be included
        self.pulse_settings['field_cutoff_ratio'] = 1.e-7
        self.pulse_settings['envelope'] = 'gaussian'
        self.pulse_settings['pulse_widths'] = 50.
        self.pulse_settings['carrier_frequencies'] = 0.325
        self.pulse_settings['centers'] = 300.0
        self.pulse_settings['CEP'] = self.pulse_settings[
            'carrier_frequencies'] * self.pulse_settings['centers']
        self.pulse_settings['field_max'] = 1.0e-5
        self.pulse_settings['pol_dir'] = 'x'
        self.pulse_settings['frequency_range'] = '0.0 - 1.0 (0.0025)'

        self.multi_input_keys = [
            'pulse_widths', 'carrier_frequencies', 'field_max', 'centers', 'CEP'
        ]

        # List of parameters that are floats for parsing
        self.float_input_keys = [
            'field_cutoff_ratio', 'pulse_widths', 'carrier_frequencies',
            'centers', 'CEP', 'field_max'
        ]

    def update_settings(self, settings, cpp_settings, method_settings=None):
        """
        Updates settings in PulsedRespnse

        :param settings:
            The settings dictionary for the driver.
        :param cpp_settings:
            The settings dictionary for complex response driver.
        :param method_settings:
            The dictionary of method settings.
        """

        if method_settings is None:
            method_settings = {}

        # Update the default args with the user provided inputs
        self.pulse_settings.update(settings)
        self.cpp_settings.update(cpp_settings)

        # Construct frequency array based on input
        self.pulse_settings['given_frequencies'] = self.pulse_settings[
            'frequency_range']
        freq_range = parse_seq_range(self.pulse_settings['frequency_range'])

        # Convert frequency range to a zero padded and a truncated list
        # of frequncies
        self.zero_padded_freqs, truncated_freqs = self.verify_freqs(freq_range)
        self.end_w = truncated_freqs[-1]

        # TODO Override of number of pulses - set to 1
        self.pulse_settings['number_pulses'] = 1

        field_w = np.zeros_like(truncated_freqs, dtype=np.complex128)

        # Check that all Pulse list parameters are given correct
        for key in self.multi_input_keys:
            if key not in self.pulse_settings:
                raise KeyError(
                    'Key "{}" not defined for PulsedResponse'.format(key))

            if self.pulse_settings['number_pulses'] == 1:
                self.pulse_settings[key] = [self.pulse_settings[key]]
                continue
            else:
                try:
                    self.pulse_settings[key] = self.pulse_settings[key].split()
                except AttributeError:
                    raise AttributeError('Could not split key "{}" : {}'.format(
                        key, self.pulse_settings[key]))

            length = len(self.pulse_settings[key])
            if length != self.pulse_settings['number_pulses']:
                raise KeyError(
                    'Lengt of key "{}" did not match "number_pulses": {}'.
                    format(key, self.pulse_settings['number_pulses']))

        # Loop over number of pulses
        for p in range(self.pulse_settings['number_pulses']):
            for key in self.float_input_keys:
                if key in self.pulse_settings:
                    try:
                        if key in self.multi_input_keys:
                            self.pulse_settings[key][p] = float(
                                self.pulse_settings[key][p])
                        else:
                            self.pulse_settings[key] = float(
                                self.pulse_settings[key])
                    except ValueError:
                        raise ValueError(
                            'Pulse key: "{}" not parseable as float: {}'.format(
                                key, self.pulse_settings[key]))

            # Compute electric field in frequency domain
            field_w += self.gauss_env_pulse_w(
                truncated_freqs,
                float(self.pulse_settings['field_max'][p]),
                float(self.pulse_settings['centers'][p]),
                float(self.pulse_settings['pulse_widths'][p]),
                float(self.pulse_settings['carrier_frequencies'][p]),
                float(self.pulse_settings['CEP'][p]),
            )

        # Threshold to limit the number of complex polarizability evaluations
        # Based on the field strength and bandwidth of the pulse
        if self.pulse_settings['field_cutoff_ratio'] > 0.:
            self.zero_pad = True
            least_field = self.pulse_settings['field_cutoff_ratio'] * max(
                field_w)
            truncated_freqs = truncated_freqs[np.abs(field_w) > least_field]
            field_w = field_w[np.abs(field_w) > least_field]

        # If the input frequency range is limited - set the 'use zero pad' flag
        # Request zero-padding if start freq > zero
        if (truncated_freqs[0] != 0.0):
            self.zero_pad = True

        self.pulse_settings['zero_pad'] = 'True' if self.zero_pad else 'False'
        self.pulse_settings['frequencies'] = truncated_freqs
        self.pulse_settings['dw'] = truncated_freqs[1] - truncated_freqs[0]
        self.pulse_settings['freq_amplitude'] = field_w

        self.truncated_freqs = truncated_freqs
        self.dw = truncated_freqs[1] - truncated_freqs[0]
        self.amplitudes = field_w

        # Set up complex response solver
        cpp_settings.update({
            'frequencies': truncated_freqs,
            'a_components': self.pulse_settings['pol_dir'],
            'b_components': self.pulse_settings['pol_dir']
        })
        self.rsp_driver = ComplexResponse(self.comm, self.ostream)
        self.rsp_driver.update_settings(cpp_settings, method_settings)

    def compute(self, molecule, ao_basis, scf_tensors):
        """
        Invoke complex response solver to compute complex response
        function at requested frequencies

        :param molecule:
            The molecule.
        :param aobasis:
            The AO basis.
        :param scf_tensors:
            The dictionary of tensors from converged SCF wavefunction.

        :return:
            A dictionary containing properties, solutions and kappas
            and possibly zeropadded properties if requested (self.zero_pad).
        """

        # Start ny printing a pulsed response header
        if self.rank == mpi_master():
            self.print_header()

            self.ostream.print_info(
                'Entering CPP Solver from Pulsed Linear Response')

        # Launch the rsp_driver calculations on the selected frequencies
        results = self.rsp_driver.compute(molecule, ao_basis, scf_tensors)

        if self.rank == mpi_master():
            results['pulse_settings'] = self.pulse_settings
            results['properties'] = {}
            # Loop over all directions
            for xyz1 in ['x', 'y', 'z']:
                for xyz2 in ['x', 'y', 'z']:
                    # Convert response function to polarizability
                    # by multiplication with -1
                    for freq in self.truncated_freqs:
                        results['properties'][(
                            xyz1, xyz2, freq
                        )] = -results['response_functions'][(xyz1, xyz2, freq)]

            self.ostream.print_info(
                'Exiting CPP Solver and returning to Pulsed Linear Response')

            # Store the results internally for saving to file is chosen
            self.results = results
            if self.zero_pad:
                results = self.apply_zero_pad(results)
                self.results.update(results)

            # footer
            self.ostream.print_blank()
            self.ostream.print_info(
                'Post-processing results and saving to files...')

            # Store results in h5 data file if requested
            if 'h5' in self.pulse_settings:
                self.write_hdf5(self.pulse_settings['h5'])

            # Store results in ascii txt file if requested
            if 'ascii' in self.pulse_settings:
                self.write_ascii(self.pulse_settings['ascii'])

            self.ostream.print_info('...done.')
            self.ostream.print_blank()

        return results

    def write_ascii(self, fname):
        """
        Writes the Pulsed response vectors to the specified output file in
        ASCII format.

        :param fname:
            Name of the checkpoint file.

        :return:
            The ASCII file saved contains the the amplitudes for all
            frequencies for each of the following directions:
                - 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz'
                    =>  Amplitudes for all directions
        """

        if not fname:
            raise ValueError('No filename given to write_ascii()')

        # Add the .h5 extension if not given
        if not fname[-3:] == '.txt':
            fname += '.txt'

        with open(str(fname), 'w') as f:
            for xyz1 in ['x', 'y', 'z']:
                for xyz2 in ['x', 'y', 'z']:
                    f.write('Frequency   Amplitude   {}{}\n'.format(xyz1, xyz2))
                    for freq, amp in zip(self.truncated_freqs, self.amplitudes):
                        cur_str = [
                            '{0:12.6f} {1:>12.8f}{2:>+12.8f}j'.format(
                                freq, np.real(amp), np.imag(amp))
                        ]
                        cur_str.append('{0.real:12.8f}{0.imag:+.8f}j'.format(
                            self.results['properties_zeropad'][(xyz1, xyz2,
                                                                freq)]))
                        f.write(' '.join(cur_str) + '\n')

    def write_hdf5(self, fname):
        """
        Writes the Pulsed response vectors to the specified output file in h5
        format. The h5 file saved contains the following datasets:

        - amplitudes
            The pulse amplitudes for the calculated truncated_freqs
        - zero_padded
            Is the dataset zero padded or not
        - 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz'
            =>  Amplitudes for all directions
        - zero_padded_freqs
            The zero padded frequency list
        - zero_padded_amplitudes
            The pulse amplitudes for the calculated frequencies zero
            padded to match th zero padded frequencies.

        :param fname:
            Name of the checkpoint file.
        """

        if not fname:
            raise ValueError('No filename given to write_hdf5()')

        # Add the .h5 extension if not given
        if not fname[-3:] == '.h5':
            fname += '.h5'

        # Convert the boolean to a a numpy array value
        zeropad = np.array([self.zero_pad])

        # Save all the internal data to the h5 datafile named 'fname'
        try:
            with h5py.File(fname, 'w') as hf:
                hf.create_dataset('frequencies',
                                  data=self.zero_padded_freqs,
                                  compression='gzip')
                hf.create_dataset('amplitudes',
                                  data=self.zero_padded_amplitudes,
                                  compression='gzip')
                hf.create_dataset('zero_padded',
                                  data=zeropad,
                                  compression='gzip')

                # Loop over all directions
                for xyz1 in ['x', 'y', 'z']:
                    for xyz2 in ['x', 'y', 'z']:
                        polarizability = []
                        # Add all polarizability for the give direction
                        for freq in self.zero_padded_freqs:
                            polarizability.append(
                                self.results['properties_zeropad'][(xyz1, xyz2,
                                                                    freq)])

                        hf.create_dataset('{}{}'.format(xyz1, xyz2),
                                          data=np.array(polarizability),
                                          compression='gzip')

        except Exception as e:
            print('Pulsed response failed to create h5 data file: {}'.format(e),
                  file=sys.stdout)

    def apply_zero_pad(self, results):
        """
        Adds zeros between 0 and the requested start frequency
        with dw increments.

        The zero padding may be necessary due to:
            1)  the limited number of frequencies for which
                the complex response function is being computed
                as dictated by 'field_cutoff_ratio'
            2)  the need for extending the frequency range to achieve
                higher temporal resolution

        Example of zero padding:
            w: 0~~~zero_pad~~~|===DATA===|~~~zero_pad~~~end_freq

        :param results:
            A dictionary containing properties, solutions and kappas
            (the results dict from .compute()).

        :return:
            A dictionary containing properties, solutions and kappas
            and zeropadded properties.
        """

        # Use the original frequencies to zero pad the results if necessary
        # |~~~~~~~~~~~~|======|~~~~~~~~~~~~|
        # 0  zero pad    DATA   zero pad   end freq
        zero_padded_results = {}
        for xyz1 in ['x', 'y', 'z']:
            for xyz2 in ['x', 'y', 'z']:
                for freq in self.zero_padded_freqs:
                    if (xyz1, xyz2, freq) not in results['properties']:
                        zero_padded_results[(xyz1, xyz2,
                                             freq)] = complex(0.0, 0.0)
                    else:
                        zero_padded_results[(xyz1, xyz2,
                                             freq)] = results['properties'][(
                                                 xyz1, xyz2, freq)]

        # Update the results dictionary with the new zero padded result
        results['properties_zeropad'] = zero_padded_results

        # Zero pad the amplitudes as well
        dw = self.truncated_freqs[1] - self.truncated_freqs[0]

        # Find the number of zeros before the truncated pulse
        prepended_zeros = int(self.truncated_freqs[0] / dw) + 1  # for zero

        # Find the number of zeros truncated at the end
        appended_zeros = (len(self.zero_padded_freqs) - len(self.amplitudes) -
                          prepended_zeros)

        # Assemble the array of amplitudes that matches the frequency array
        pre_zeros = [0 for x in range(prepended_zeros)]
        app_zeros = [0 for x in range(appended_zeros)]
        self.zero_padded_amplitudes = np.array(pre_zeros +
                                               self.amplitudes.tolist() +
                                               app_zeros)

        return results

    def verify_freqs(self, frequencies):
        """
        Takes a list or nd-array of frequencies
        and verifies that it have or can be extended to include
        zero frequency

        :param frequencies:
            List of frequencies.

        :return:
            Numpy array of zero-padded frequencies and
            numpy array of truncated frequencies.
        """
        if isinstance(frequencies, list):
            frequencies = np.array(frequencies)

        # If the frequency already contains 0-frequency, return
        if frequencies[0] == 0.0:
            return frequencies, frequencies

        # Check if the frequencies can be extended to 0
        dw = np.round(frequencies[1] - frequencies[0], 6)
        if frequencies[0] % dw == 0.0:
            zero_padded_frequencies = np.arange(0.0, frequencies[-1], dw)
            return zero_padded_frequencies, frequencies

        # Create a new equidistant frequency list
        if self.rank == mpi_master():
            self.ostream.print_info(
                'Pulsed response module adjusts frequencies to intersect with 0'
            )
        zero_padded_frequencies = np.arange(0.0, frequencies[-1], dw)
        truncated_frequencies = zero_padded_frequencies[
            zero_padded_frequencies >= frequencies[0]]
        return zero_padded_frequencies, truncated_frequencies

    def gauss_env_pulse_w(self, w, F0, t0, delta_t, w_carrier, cep):
        """
        Gaussian pulse from frequency domain input.

        :param w:
            np.array - list of frequencies.
        :param F0:
            float - pulse amplitude in [au].
        :param t0:
            float - time center for pulse in [au].
        :param delta_t:
            float - pulse width in [au].
        :param w_carrier:
            float - carrier frequency in [au].
        :param cep:
            float - carrier envelope phases in [radians].

        :return:
            Numpy array of the pulse amplitude in the frequency domain.
        """

        return ((F0 * delta_t) / (2. * (2. * np.pi)**0.5) *
                (np.exp(-(delta_t**2 * (w - w_carrier)**2) / 2.0) *
                 np.exp(1.j * (w - w_carrier) * t0 + 1.j * cep) +
                 np.exp(-(delta_t**2 * (w + w_carrier)**2) / 2.0) *
                 np.exp(1.j * (w + w_carrier) * t0 - 1.j * cep)))

    def print_header(self):
        """
        Prints Pulsed Lionear Response calculation setup details to output
        stream.
        """

        # Print string width (global norm)
        str_width = 60

        # PRT header
        self.ostream.print_blank()
        title = 'Pulsed Linear Reponse Theory Calculation'
        self.ostream.print_header(title)
        self.ostream.print_header('=' * (len(title) + 2))
        self.ostream.print_blank()

        # Print all settings
        header_fields = {
            'field_cutoff_ratio': 'Field cutoff ratio',
            'envelope': 'Envelope',
            'pulse_widths': 'Pulse Duration',
            'carrier_frequencies': 'Carrier Frequency',
            'centers': 'Pulse Center time',
            'field_max': 'Max Field',
            'pol_dir': 'Polarization Direction',
            'frequency_range': 'Frequency Range',
            'zero_pad': 'Zero padding results',
        }

        # Print the header information fields
        for key, text in header_fields.items():
            cur_str = '{0:30s} : {1:s}'.format(text,
                                               str(self.pulse_settings[key]))
            self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_blank()

        # Print the list of truncated frequencies and their amplitudes
        valstr = '{:<12s}  |  {:>18s}'.format('Frequency', 'Amplitude')
        self.ostream.print_header(valstr.ljust(str_width))
        self.ostream.print_header(('-' * 45).ljust(str_width))

        for freq, amp in zip(self.truncated_freqs, self.amplitudes):
            cur_str = '{:<12.6f}  :  {:>12.8f}   {:>+12.8f}j'.format(
                freq, amp.real, amp.imag)
            self.ostream.print_header(cur_str.ljust(str_width))
        self.ostream.print_blank()
