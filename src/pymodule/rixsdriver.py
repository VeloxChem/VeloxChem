#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
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
import math
import h5py
import sys

from .veloxchemlib import bohr_in_angstrom, mpi_master
from .veloxchemlib import hartree_in_ev, bohr_in_angstroms
from .oneeints import compute_electric_dipole_integrals
from .subcommunicators import SubCommunicators
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical

from .inputparser import (parse_input, print_keywords, print_attributes,
                          get_random_string_parallel)

class RixsDriver:
    """
    Implements RIXS driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - theta: Angle.
        - gamma_n: Broadening term.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes (CVS-)RIXS driver.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # mpi information
        self.comm = comm
        self.rank = comm.Get_rank()
        self.nodes = comm.Get_size()

        # outputstream
        self.ostream = ostream

        # method settings
        self.nr_ve = None # nr_ce

        # CVS
        self.nr_ce = None # nr_ce
        self.nr_CO = None # nr_ce
        self.photon_energy = None

        # SRA
        self.sra = False

        self.theta = 0
        self.gamma_n = .124/hartree_in_ev() # a.u.

        # input keywords
        self.input_keywords = {
            'rixs': {
                'theta': ('float', 'angle between incident polarization vector and propagation vector of outgoing'),
                'gamma': ('float', 'broadening term'),
                'nr_CO': ('int', 'number of involved core-orbitals'),
            },
        }

    def update_settings(self, rixs_dict):
        """
        Updates settings in SCF driver.

        :param rixs_dict:
            The dictionary of rixs input.
        """

        if method_dict is None:
            method_dict = {}

        rixs_keywords = {
            key: val[0] for key, val in self._input_keywords['rixs'].items()
        }

        parse_input(self, rixs_keywords, rixs_dict)

        if 'program_end_time' in rixs_dict:
            self.program_end_time = rixs_dict['program_end_time']
        if 'filename' in rixs_dict:
            self.filename = rixs_dict['filename']

        """
        if self.electric_field is not None:
            assert_msg_critical(
                len(self.electric_field) == 3,
                'SCF driver: Expecting 3 values in \'electric field\' input')
            assert_msg_critical(
                not self._pe,
                'SCF driver: \'electric field\' input is incompatible with ' +
                'polarizable embedding')
            # disable restart of calculation with static electric field since
            # checkpoint file does not contain information about the electric
            # field
            self.restart = False
        """
    
    def sts_tdm(self, scf_results, valence_eigenvectors, core_eigenvectors):
        """
        The state-to-state transition-density-matrix

        :param molecule    : molecule object (could be removed w/ small modifications)
        :param ao_basis    : atomic orbital basis (could be removed)
        :param scf_results : results of the SCF calculation
        :tda_res_val       : linear-response results for the valence calculation
        :tda_res_core      : linear-response results for the core calculation

        returns the state-to-state transition-density-matrix in ao-basis
        """

        if self.sra:
            X_val_mo = np.zeros((self.nocc, self.nvir, self.nr_ve))
            X_val_mo[:self.nr_CO] = valence_eigenvectors[:self.nr_CO]
            X_val_mo[-self.nr_val:] = valence_eigenvectors[self.nr_CO:]
            X_cor_mo = np.zeros((self.nocc, self.nvir, self.nr_ce))
            X_cor_mo[:self.nr_CO] = core_eigenvectors[:self.nr_CO]
            X_cor_mo[-self.nr_val:] = core_eigenvectors[self.nr_CO:]
        elif self.nr_CO > 0:
            X_cor_mo = np.zeros((self.nocc, self.nvir, self.nr_ce))
            X_cor_mo[:self.nr_CO] = core_eigenvectors.reshape(self.nr_CO, self.nvir, self.nr_ce)
            X_val_mo = valence_eigenvectors.reshape(self.nocc, self.nvir, self.nr_ve)
        else:
            X_cor_mo = core_eigenvectors.reshape(self.nocc, self.nvir, self.nr_ce)
            X_val_mo = valence_eigenvectors.reshape(self.nocc, self.nvir, self.nr_ve)

        # Electron- and hole density matrices
        gamma_ab = np.einsum('ipJ, iqI -> pqIJ', X_cor_mo, X_val_mo, optimize=True)
        gamma_ij = np.einsum('qaJ, paI -> pqIJ', X_cor_mo, X_val_mo, optimize=True)

        # Transform
        if self.sra:
            C_ab = scf_results['C_alpha'][:,self.nocc:self.nocc + self.nvir]
        else:
            C_ab = scf_results['C_alpha'][:,self.nocc:]
        C_ij = scf_results['C_alpha'][:,:self.nocc]
        gamma_ab_ao = np.einsum('vp, pqIJ, wq -> vwIJ', C_ab, gamma_ab, C_ab, optimize=True)
        gamma_ij_ao = np.einsum('vp, pqIJ, wq -> vwIJ', C_ij, gamma_ij, C_ij, optimize=True)

        gamma_ao = gamma_ab_ao - gamma_ij_ao
        gamma_ao *= np.sqrt(2) #2 # alpha and beta spin
        return gamma_ao

    def gts_tdm(self, scf_results, core_eigenvectors):
        """
        The GS-to-state transition-density-matrix (GS to core) in ao --
        essentially transformed and reshaped eigenvector

        :param molecule: 
            molecule object (could be removed w/ small modifications)
        :param ao_basis: 
            atomic orbital basis (could be removed)
        :param scf_results: 
            results of the SCF calculation (could be removed w/ small modifications)
        :tda_res_val:
            linear-response results for the valence calculation
        :tda_res_core:
            linear-response results for the core calculation

        :returns: 
            the GS-to-state transition-density-matrix in ao-basis
        """

        if self.sra:
            X_cor_mo = np.zeros((self.nocc, self.nvir, self.nr_ce))
            X_cor_mo[:self.nr_CO] = core_eigenvectors[:self.nr_CO]
            X_cor_mo[-self.nr_val:] = core_eigenvectors[self.nr_CO:]
        elif self.nr_CO > 0:
            X_cor_mo = np.zeros((self.nocc, self.nvir, self.nr_ce))
            X_cor_mo[:self.nr_CO] = core_eigenvectors.reshape(self.nr_CO, self.nvir, self.nr_ce)
        else:
            X_cor_mo = core_eigenvectors.reshape(self.nocc, self.nvir, self.nr_ce)

        # Transform
        if self.sra:
            C_ab = scf_results['C_alpha'][:,self.nocc:self.nocc + self.nvir]
        #    C_ij = scf_results['C_alpha'][:,:self.nr_CO]
        else:
            C_ab = scf_results['C_alpha'][:,self.nocc:]
        C_ij = scf_results['C_alpha'][:,:self.nocc]
        gamma_ng_ao = np.einsum('vp, pqJ, wq -> vwJ', C_ij, X_cor_mo, C_ab, optimize=True)
        gamma_ng_ao *= np.sqrt(2) # 2 # alpha and beta spin
        return gamma_ng_ao

    def F_xy(self, w, f, gamma_ao,
         gamma_ng_ao, core_eigenvalues,
         dipole_ints):
        """
        The RIXS scattering amptlitude (sum-over-states transition strength).
        
        :param w:
            Energy of incoming photon 
        :param f: 
            Which, out of the nr_ve, final valence states to end at
        :param gamma_factor: 
            Broadening factor, or half-width of the core-excited state; 
            in the resonant case corresponding to state with eigenvalue w
    
        :return: 
            The scattering amplitude/transition strength 
            for all molecular directions and outgoing photon polarisation
        """
        
        eigenvalues = core_eigenvalues # excitation energies
        e_n = (1 / (eigenvalues - w - self.gamma_n*1j)) # denominator
        #print(f'e_n:{e_n.shape}, dip:{dipole_ints.shape} gamma:{gamma_ao[:,:,f].shape}, gamma_ng:{gamma_ng_ao.shape}')
        scatt_amp = np.einsum('n, xij, ijn, yab, abn -> xy', e_n, -dipole_ints, gamma_ao[:,:,f], -dipole_ints, gamma_ng_ao, optimize='greedy')
        return scatt_amp
    
    def osc_str(self, tdpm, omega_prime):# tda_res_core, tda_res_val):
        """
        Computes the oscillator strength between the core and valence states

        returns: oscillator strengths as array, shape = (core_exc, valence_exc)
        """
        if tdpm.ndim == 3:
            unweighted_f = (2/3) * np.einsum('fnx->nf', tdpm ** 2)
        elif tdpm.ndim == 2:
            unweighted_f = (2/3) * np.einsum('nx->n', tdpm ** 2)
        f = omega_prime * unweighted_f
        return f

    def rixs_xsection(self, w, f, gamma_ao,
         gamma_ng_ao, core_eigenvalues, valence_eigenvalues,
         dipole_ints, theta, w_p=None):
        """
        Calculate the RIXS cross-section, sigma.
    
        :param w            : Energy of the incident photon
        :param w_p          : Energy of the outgoing/scattered photon, if None: these are matched with valence energies
        :param f            : Index of the final/target excited state
        :param theta        : Scattering angle 
        :param gamma_factor : Broadening factor for Lorentzian profile
        :param tda_res(_val) : The results tensor of the TDDFT/TDA ground to valence excited state
        
        returns the RIXS transition intensity for the given frequencies, w (and) and state f
        """
        if w_p == None:
            w_prime = w - valence_eigenvalues[f] #w - w_f0 
        else:
            w_prime = w_p
            
        F = self.F_xy(w, f, gamma_ao,
         gamma_ng_ao, core_eigenvalues,
         dipole_ints)
        self.F_mat = F
        
        sigma = w_prime/w * 1/15 * ((2 - (1/2) * np.sin(theta) ** 2) * np.sum(np.abs(F)**2) 
                   + ((3/4) * np.sin(theta) ** 2 - 1/2) * (np.sum(F * F.T.conj())
                                                   + np.trace(np.abs(F)**2)))
        return sigma

    def transition_dipole_mom(self, gamma_ao, dipole_integrals):
        if gamma_ao.ndim == 3:
            T_fn = np.einsum('ijn,xij->nx', gamma_ao, -dipole_integrals)
        elif gamma_ao.ndim == 4:
            T_fn = np.einsum('ijfn,xij->fnx', gamma_ao, -dipole_integrals)
        return T_fn

    def omega_p(self, core_eigenvalues, val_eigenvalues):
        d_E = [ce - val_eigenvalues for ce in core_eigenvalues]
        return np.array(d_E)
    
    def find_nr_CO(self, excitation_details):
        largest_core = 0  # Initialize the largest core index as 0

        for detail in excitation_details:
            entry = detail[0]
            if 'core_' in entry:
                core_index = int(entry.split('core_')[1].split()[0])
                largest_core = max(largest_core, core_index)
        return largest_core

    def elastic_rixs_xsection(self, eigenvals_core, dipole_ints, gamma_ng_ao, theta):
        sigma = np.zeros((self.nr_ce), dtype=np.complex128)
        # excitation energies
        e_n = np.array([(1 / (en - eigenvals_core - self.gamma_n*1j)) for en in eigenvals_core]) # denominator
        #print(f'e_n:{e_n.shape}, dip:{dipole_ints.shape} gamma:{gamma_ao[:,:,f].shape}, gamma_ng:{gamma_ng_ao.shape}')

        F = np.einsum('kn, xij, ijn, yab, abn -> kxy', e_n, -dipole_ints, gamma_ng_ao, -dipole_ints, gamma_ng_ao, optimize='greedy')
        #F = e_n * np.einsum('xij, ijf, yab, abf -> nxy', -dipole_ints, gamma_ng_ao, -dipole_ints, gamma_ng_ao, optimize='greedy')

        for n in range(self.nr_ce):
            _F = F[n]
            sigma[n] = 1/15 * ((2 - (1/2) * np.sin(theta) ** 2) * np.sum(np.abs(_F)**2) 
                   + ((3/4) * np.sin(theta) ** 2 - 1/2) * (np.sum(_F * _F.T.conj())
                                                   + np.trace(np.abs(_F)**2)))
        return sigma

    def compute(self, molecule, ao_basis, scf_results, tda_res_val, tda_res_core=None, nr_CO=None, nr_vir=None, nr_val=None):
        self.norb = scf_results['C_alpha'].shape[0]
        self.nocc = int(sum(scf_results['occ_alpha'])) #molecule.number_of_alpha_electrons()
        self.nvir = self.norb - self.nocc

        if not tda_res_core:
            # SRA
            assert_msg_critical(nr_CO is not None,
                                'RixsDriver: need the number of core, valence, and virtual orbitals involved for subspace-restricted approx.')
            self.sra = True
            self.nr_CO = nr_CO
            self.nr_val = nr_val
            self.nvir = nr_vir
            self.nr_ce = nr_CO * nr_vir
            self.nr_ve = nr_val * nr_vir

            eigenvector = tda_res_val['eigenvectors'].reshape(self.nr_CO + self.nr_val, self.nvir, self.nr_ce + self.nr_ve)
            core_eigenvectors     = eigenvector[:, : , self.nr_ve:]
            core_eigenvalues      = tda_res_val['eigenvalues'][self.nr_ve:]
            valence_eigenvectors  = eigenvector[:, : , :self.nr_ve]
            valence_eigenvalues   = tda_res_val['eigenvalues'][:self.nr_ve]
        else: 
            # CVS, i.e., two-shot approach
            self.nr_ve = len(tda_res_val['excitation_details'])
            self.nr_ce = len(tda_res_core['excitation_details'])
            self.nr_CO = self.find_nr_CO(tda_res_core['excitation_details'])

            core_eigenvectors     = tda_res_core['eigenvectors']
            core_eigenvalues      = tda_res_core['eigenvalues']
            valence_eigenvectors  = tda_res_val['eigenvectors']
            valence_eigenvalues   = tda_res_val['eigenvalues']
        
        tdm_ng = self.gts_tdm(scf_results, core_eigenvectors)
        tdm_fn = self.sts_tdm(scf_results, valence_eigenvectors, core_eigenvectors)

        dip_tuple   = compute_electric_dipole_integrals(molecule, ao_basis)
        dipole_ints = -1.0 * np.array([dip_tuple[0],
                                       dip_tuple[1],
                                       dip_tuple[2]])

        omega_f = valence_eigenvalues
        omega_n = core_eigenvalues
        if self.photon_energy is None:
            self.photon_energy = omega_n # resonant
        
        ene_loss              = np.zeros((self.nr_ve, len(self.photon_energy)))
        emiss                 = np.zeros((self.nr_ve, len(self.photon_energy)))
        crossections          = np.zeros((self.nr_ve, len(self.photon_energy)))
        scattering_amp_tensor = np.zeros((3, 3, self.nr_ve, self.nr_ce), dtype=np.complex128)

        for j,vs_i in enumerate(range(self.nr_ve)):
            for k, w_n in enumerate(self.photon_energy):
                emiss[vs_i, k]        = w_n - omega_f[vs_i]
                ene_loss[vs_i, k]     = w_n - (w_n - omega_f[vs_i])
                crossections[vs_i, k] = self.rixs_xsection(w_n, vs_i, tdm_fn, 
                                                            tdm_ng, core_eigenvalues, valence_eigenvalues, 
                                                            dipole_ints, self.theta).real
                scattering_amp_tensor[:,:, j, k] = self.F_mat

        self.emission      = emiss
        self.ene_loss      = ene_loss
        self.crossections  = crossections
        self.scatt_amp_mat = scattering_amp_tensor
        
        emission_ene             = self.omega_p(core_eigenvalues, valence_eigenvalues)
        self.emission_energy_map = emission_ene
        self.ene_loss_map        = np.array([ce - emission_ene[i] for i, ce in enumerate(core_eigenvalues)])

        T_fn = self.transition_dipole_mom(tdm_fn, dipole_ints)
        self.oscillator_strength  = self.osc_str(T_fn, emission_ene) 

        self.elastic_rixs = self.elastic_rixs_xsection(omega_n, dipole_ints, tdm_ng, self.theta).real

    def write_hdf5(self, fname):
        """
        TODO
        Writes the RIXS {output?} to the specified output file in h5
        format. The h5 file saved contains the following datasets:

        - amplitudes
            The scattering amplitudes for the calculated truncated_freqs,
            as 'xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz'
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
                hf.create_dataset('frequencies', data=self.zero_padded_freqs)
                hf.create_dataset('amplitudes',
                                  data=self.zero_padded_amplitudes)
                hf.create_dataset('zero_padded', data=zeropad)

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
                                          data=np.array(polarizability))

        except Exception as e:
            print('Pulsed response failed to create h5 data file: {}'.format(e),
                  file=sys.stdout)
            
    def print_header(self):
        """
        TODO
        Prints RIXS calculation setup details to output
        stream.
        """

        # Print string width (global norm)
        str_width = 60

        # PRT header
        self.ostream.print_blank()
        title = 'RIXS Linear Reponse CVS Calculation'
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


    


