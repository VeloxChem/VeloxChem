import numpy as np
import copy

from .rspproperty import ResponseProperty
from .crsp import ComplexResponse
from .inputparser import parse_frequencies


class PulsedResponse():
    """
    Pulsed Reponse class for computing molecular responses to
    a resonant or non-resonant laser pulse.
    """

    def __init__(self, comm, ostream, pulse_input, crsp_input):
        """
        Initializes pulsed response.
        """

        self.comm = comm
        self.ostream = ostream
        self.zero_pad = False
        pulse_settings = {}
        crsp_settings = {}

        # TODO: make a E-field class that takes pulse settings
        #       as input and generate field_w and frequency window 
        # Default Pulse settings
        pulse_settings['smallest_field_ratio']  = 1.e-7 # Smallest fraction field amplitude to be included
        pulse_settings['envelope']              = 'gaussian'
        pulse_settings['pulse_duration']        = 50. 
        pulse_settings['carrier_frequency']     = 0.325 
        pulse_settings['center']                = 300.0 
        pulse_settings['CEP']                   = pulse_settings['carrier_frequency'] * pulse_settings['center'] 
        pulse_settings['field_max']             = 1.0e-5
        pulse_settings['pol_dir']               = (1.0, 0., 0.)
        pulse_settings['frequencies']           = "0.0 - 1.0 (0.0025)"
        pulse_settings['file_prefix']           = ''

        # Default CRSP settings (if nothing else given, it will use defaults in crsp.cpp)
#        crsp_settings['rot_averaging']           = False

        # Update the default args with the user provided inputs
        pulse_settings.update(pulse_input)
        crsp_settings.update(crsp_input)
        
        # TODO: Setup response operators based on polarization vector
        #       The crsp.py file needs to be updated with a/b operator input option
         
        # Construct frequency array based on input
        freqs = parse_frequencies(pulse_settings['frequencies'])
        self.w, freqs = self.verify_frequencies(freqs)
        self.end_w = freqs[-1]

        # Compute electric field in frequency domain
        field_w = self.gauss_env_pulse_w(
                        freqs,
                        pulse_settings['field_max'],
                        pulse_settings['center'],
                        pulse_settings['pulse_duration'],
                        pulse_settings['carrier_frequency'],
                        pulse_settings['CEP'],
                        )

        # Save frequency and field (frequency domain) in numpy arrays
        np.save('{}freq.npy'.format(pulse_input['file_prefix']),freqs)
        np.save('{}field_w.npy'.format(pulse_input['file_prefix']),field_w)

        # Save frequencies for later zero padding
        self.freqs = copy.copy(freqs)

        # Threshold to limit the number of complex polarizability evaluations
        # Based on the field strength and bandwidth of the pulse
        if pulse_settings['smallest_field_ratio'] > 0.:
            self.zero_pad = True
            least_field = pulse_settings['smallest_field_ratio'] * max(field_w)
            freqs = freqs[np.abs(field_w) > least_field]
            field_w = field_w[np.abs(field_w) > least_field]

        # If the input frequency range is limited - set the 'use zero pad' flag
        # Request zero-padding if start freq > zero
        if (self.freqs[0] != 0.0):
            self.zero_pad = True

        pulse_settings['frequencies'] = freqs
        pulse_settings['dw'] = freqs[1] - freqs[0]
        pulse_settings['freq_amplitude'] = field_w

        # Store the input parameters
        self.pulse_settings = pulse_settings

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

        results = self.rsp_driver.compute(molecule, ao_basis, scf_tensors)
        results['pulse_settings'] = self.pulse_settings

        if self.zero_pad:
            results = self.apply_zero_pad(results)
        return results

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

