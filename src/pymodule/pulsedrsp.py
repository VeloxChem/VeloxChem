import numpy as np
import copy
import h5py
from .veloxchemlib import mpi_master
from .rspproperty import ResponseProperty
from .crsp import ComplexResponse
from .inputparser import parse_frequencies


class PulsedResponse():
    """
    Pulsed Reponse class for computing molecular responses to
    a resonant or non-resonant laser pulse.


    :param comm:
        The MPI communicator.
    :param rank:
        The MPI rank.
    :param smallest_field_ratio:
        
    :param envelope:
        Envelope of pulse - available arguments: ['gaussian', ]
    :param number_pulses:
        The number of different pulse settings to expect 'N' 
    :param pulse_widths:
        List(N) of pulse widths give in [au]
    :param carrier_frequencies:
        List(N) of carrier frequencies given in [au]
    :param centers:
        List(N) of time centers for the pulses given in [au]
    :param pol_dir:
        List(N) of polarization directions, arguments given in combinations of x, y and z
        e.g.: 'x' for [1 0 0], yz for [0 1 1].
    :param frequency_range:
        Frequencies to map solution to - given as range 'start-end(df)' in [au]
        e.g. 0.2-0.4(0.007) for 0.2 -> 0.4 au in steps of 0.007 au
        zero-padding will be applied if range does not start at 0
    :param CEP:
        
    :param field_max:

    :param h5file:
        If a name is given - a h5 formatted file is saved using gzip compression
    :param ascii:
        If a name is given - a ASCII formatted file is saved
    """

    # List of parameters that are floats for parsing
    float_input_keys = ['smallest_field_ratio',
                        'pulse_duration',
                        'carrier_frequencies',
                        'center',
                        'CEP',
                        'field_max',
                        ]

    def __init__(self, comm, ostream):
        """
        Initializes pulsed response.
        """

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.comm = comm
        self.ostream = ostream
        self.zero_pad = False
        self.pulse_settings = {}
        self.crsp_settings = {}

        # TODO: make a E-field class that takes pulse settings
        #       as input and generate field_w and frequency window 
        # Default Pulse settings
        self.pulse_settings['smallest_field_ratio']  = 1.e-7 # Smallest fraction field amplitude to be included
        self.pulse_settings['envelope']              = 'gaussian'
        self.pulse_settings['pulse_duration']        = 50. 
        self.pulse_settings['carrier_frequencies']   = 0.325 
        self.pulse_settings['center']                = 300.0 
        self.pulse_settings['CEP']                   = self.pulse_settings['carrier_frequencies'] * self.pulse_settings['center'] 
        self.pulse_settings['field_max']             = 1.0e-5
        self.pulse_settings['pol_dir']               = 'x'
        self.pulse_settings['frequency_range']       = "0.0 - 1.0 (0.0025)"



    def update_settings(self, settings, crsp_settings):
        """
        Updates settings in PulsedRespnse

        :param settings:
            The settings dictionary for the driver.

        :param crsp_settings:
            The settings dictionary for complex response driver.
        """
        # Default CRSP settings (if nothing else given, it will use defaults in crsp.cpp)
#        crsp_settings['rot_averaging']           = False

        # Update the default args with the user provided inputs
        self.pulse_settings.update(settings)
        self.crsp_settings.update(crsp_settings)

        for key in PulsedResponse.float_input_keys:
            if key in self.pulse_settings: 
                try:
                    self.pulse_settings[key] = float(self.pulse_settings[key])
                except:
                    raise TypeError("Pulsed response key: '{}' was not parseable as a float: {}".format(key, self.pulse_settings[key]))
        

        # TODO: Setup response operators based on polarization vector
        #       The crsp.py file needs to be updated with a/b operator input option
         
        # Construct frequency array based on input
        self.pulse_settings['given_frequencies'] = self.pulse_settings['frequency_range']
        freqs = parse_frequencies(self.pulse_settings['frequency_range'])
        
        self.w, freqs = self.verify_frequencies(freqs)
        self.end_w = freqs[-1]

        # Compute electric field in frequency domain
        field_w = self.gauss_env_pulse_w(
                        freqs,
                        self.pulse_settings['field_max'],
                        self.pulse_settings['center'],
                        self.pulse_settings['pulse_duration'],
                        self.pulse_settings['carrier_frequencies'],
                        self.pulse_settings['CEP'],
                        )

        # Save frequencies for later zero padding
        self.freqs = copy.copy(freqs)

        # Threshold to limit the number of complex polarizability evaluations
        # Based on the field strength and bandwidth of the pulse
        if self.pulse_settings['smallest_field_ratio'] > 0.:
            self.zero_pad = True
            least_field = self.pulse_settings['smallest_field_ratio'] * max(field_w)
            freqs = freqs[np.abs(field_w) > least_field]
            field_w = field_w[np.abs(field_w) > least_field]

        # If the input frequency range is limited - set the 'use zero pad' flag
        # Request zero-padding if start freq > zero
        if (self.freqs[0] != 0.0):
            self.zero_pad = True
    
        self.pulse_settings['zero_pad'] = "True" if self.zero_pad else "False"
        self.pulse_settings['frequencies'] = freqs
        self.pulse_settings['dw'] = freqs[1] - freqs[0]
        self.pulse_settings['freq_amplitude'] = field_w

        self.frequencies = freqs
        self.dw = freqs[1] - freqs[0]
        self.amplitudes = field_w 

        # Set up complex response solver
        crsp_settings.update({'frequencies' : freqs})
        self.rsp_driver = ComplexResponse(self.comm, self.ostream)
        self.rsp_driver.update_settings(crsp_settings)


    def compute(self,  molecule, ao_basis, scf_tensors):
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
            and possibly zeropadded properties if requested (self.zero_pad)
        """


        # Start ny printing a pulsed response header
        if self.rank == mpi_master():
            self.print_header()


        self.ostream.print_info("Entering CPP Solver from Pulsed Linear Response")

        # Launch the rsp_driver calculations on the selected frequencies
        results = self.rsp_driver.compute(molecule, ao_basis, scf_tensors)
        self.ostream.print_info("Exiting CPP Solver returning to Pulsed Linear Response")

        results['pulse_settings'] = self.pulse_settings
        
        self.results = results
        if self.zero_pad:
            results = self.apply_zero_pad(results)
            self.results.update(results)

        # Store the results internally for saving to file is chosen
        self.results = results


        # footer
        self.print_footer()


        # Store results in h5 data file if requested
        if 'h5' in self.pulse_settings:
            self.write_hdf5(self.pulse_settings['h5'])

        # Store results in ascii txt file if requested
        if 'ascii' in self.pulse_settings:
            self.write_ascii(self.pulse_settings['ascii'])



        return results


    def write_ascii(self, fname):
        """
        Writes the Pulsed response vectors to the specified 
        output file in ASCII format 

        :param fname:
            Name of the checkpoint file.

        :return:
            The ASCII file saved contains the the amplitudes for all frequencies
            for each of the following directions 
                - 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz'
                   =>  Amplitudes for all directions
        """

        if not fname: raise ValueError("No filename given to write_hdf5()") 


        # Add the .h5 extension if not given
        if not fname[-3:] == ".txt":
            fname += ".txt"

        directions = ['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz']

        with open(fname, 'w') as f:
            for xyz1 in ['x', 'y','z']:
              for xyz2 in ['x', 'y','z']:
                f.write("Frequency   Amplitude   {}{}\n".format(xyz1, xyz2))
                for freq, amp in zip (self.frequencies, self.amplitudes):
                    cur_str = ["{0:12.6f} {1:>12.8f}{2:>+12.8f}j".format(freq, np.real(amp), np.imag(amp))]
                    cur_str.append("{0.real:12.8f}{0.imag:+.8f}j".format(self.results['properties_zeropad'][(xyz1, xyz2, freq)]))
                    f.write(" ".join(cur_str)+"\n")
            


    def write_hdf5(self, fname):
        """
        Writes the Pulsed response vectors to the specified 
        output file in h5 format 

        :param fname:
            Name of the checkpoint file.

        :return:
            The h5 file saved contains the following datasets:
                - frequencies : The calculated frequencies
                - amplitudes : The pulse amplitudes for the calculated frequencies
                - zero_padded_frequencies : The zero padded frequency list
                - zero_padded : Is the dataset zero padded or not
                - 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz'
                   =>  Amplitudes for all directions
        """

        if not fname: raise ValueError("No filename given to write_hdf5()") 


        # Add the .h5 extension if not given
        if not fname[-3:] == ".h5":
            fname += ".h5"

        # Convert the boolean to a a numpy array value
        zeropad = np.array([self.zero_pad])

        """ Save all the internal data to the h5 datafile named 'fname'
        """
        try:
            with h5py.File(fname, 'w') as hf:
                hf.create_dataset('frequencies', data=self.freqs, compression="gzip")
                hf.create_dataset('amplitudes', data=self.amplitudes, compression="gzip")
                hf.create_dataset('zero_padded_frequencies', data=self.w, compression="gzip")
                hf.create_dataset('zero_padded', data=zeropad, compression="gzip")

                # Loop over all directions
                for xyz1 in ['x', 'y','z']:
                    for xyz2 in ['x', 'y','z']:
                        amplitudes = []
                        # Add all amplitudes for the give direction 
                        for freq in self.w:
                            if self.zero_pad:
                                amplitudes.append(self.results['properties_zeropad'][(xyz1, xyz2, freq)])
                            elif (xyz1, xyz2, freq) in self.results['properties']:
                                amplitudes.append(self.results['properties_zeropad'][(xyz1, xyz2, freq)])

                        hf.create_dataset("{}{}".format(xyz1,xyz2), data=np.array(amplitudes), compression="gzip")

            print("Pulsed response data saved to '{}'".format(fname))
        except Exception as e:
            print("Pulsed response failed to create h5 data file: {}".format(e))
            
        return


    def apply_zero_pad(self, results):
        """
        Adds zeros between 0 and the requested start frequency
        with dw increments

        The zero padding may be necessary due to:
            1)  the limited number of frequencies for which 
                the complex response function is being computed
                as dictated by "smallest_field_ratio"
            2)  the need for extending the frequency range to achieve
                higher temporal resolution
        
        Example of zero padding 
         w:     |~~~~~~~~~~~~|======|~~~~~~~~~~~~|
                0  zero pad    DATA   zero pad   end freq

        :param results:
            A dictionary containing properties, solutions and kappas
            (the results dict from .compute()

        :return:
            A dictionary containing properties, solutions and kappas
            and zeropadded properties
        """

        # Use the original frequencies to zero pad the results if necessary
        # |~~~~~~~~~~~~|======|~~~~~~~~~~~~|
        # 0  zero pad    DATA   zero pad   end freq
        zero_padded_results = {}
        for xyz1 in ['x', 'y','z']:
            for xyz2 in ['x', 'y','z']:
                for freq in self.w:
                    if (xyz1, xyz2, freq) not in results['properties']:
                        zero_padded_results[(xyz1, xyz2, freq)] = complex(0.0, 0.0)
                    else:
                        zero_padded_results[(xyz1, xyz2, freq)] = results['properties'][(xyz1, xyz2, freq)]
        
        # Update the results dictionary with the new zero padded result
        results['properties_zeropad'] = zero_padded_results

        return results 
        
    def verify_frequencies(self, frequencies):
        """
        Takes a list or nd-array of frequencies
        and verifies that it have or can be extended to include
        zero frequency

        :param frequencies:
            List of frequencies

        :return:
            numpy array of frequencies
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
        self.ostream.print_info("Pulsed response module adjusts frequencies to intersect with 0")
        zero_padded_frequencies = np.arange(0.0, frequencies[-1], dw)
        truncated_frequencies = zero_padded_frequencies[zero_padded_frequencies >= frequencies[0]]
        return zero_padded_frequencies, truncated_frequencies
        

    def gauss_env_pulse_w(self, w, F0, t0, delta_t, w_carrier, cep):

        return ((F0 * delta_t)/(2.*(2.*np.pi)**0.5) 
                * ( np.exp( -(delta_t**2 * (w - w_carrier)**2 ) / 2.0 ) * np.exp(1.j * (w - w_carrier) * t0 - 1.j * cep) 
                  + np.exp( -(delta_t**2 * (w + w_carrier)**2 ) / 2.0 ) * np.exp(1.j * (w + w_carrier) * t0 + 1.j * cep))
               )

    def print_footer(self):
        """ A footer to be printed after finalizing the calculation  
        """

        # PRT footer
        self.ostream.print_blank()
        self.ostream.print_header("Post processing results and saving to files")
        self.ostream.print_header(44 * "-")
        self.ostream.print_blank()


    def print_header(self):
        """
        Prints Pulsed Lionear Response calculation setup details to output stream,
        """

        # Print string width (global norm)
        str_width = 60

        # PRT header
        self.ostream.print_blank()
        self.ostream.print_header("Pulsed Linear Reponse Theory Calculation")
        self.ostream.print_header(44 * "=")
        self.ostream.print_blank()

        # Print all settings
        header_fields = {
                "smallest_field_ratio" : "Smallest Field ratio",
                "envelope" : "Envelope",
                "pulse_duration" : "Pulse Duration",
                "carrier_frequencies" : "Carrier Frequency",
                "center" : "Pulse Center time",
                "field_max" : "Max Field",
                "pol_dir" : "Polarization Direction",
                "frequency_range" : "Frequency Range",
                "zero_pad" : "Zero padding results",
            }

        # Print the header information fields
        for key, text in header_fields.items():
            cur_str = "{0:30s} : {1:s}".format(text, str(self.pulse_settings[key]))
            self.ostream.print_header(cur_str.ljust(str_width))


        # Print the list of frequencies and their amplitudes
        self.ostream.print_blank()
        self.ostream.print_header("Frequency | Amplitude")
        self.ostream.print_header(21 * "=")

        str_width = 40
        for freq, amp in zip (self.frequencies, self.amplitudes):
            cur_str = "{0:12.6f} : {1:>12.8f}{2:>+12.8f}j".format(freq, np.real(amp), np.imag(amp))
            self.ostream.print_header(cur_str.ljust(str_width))
            

        self.ostream.print_blank()

