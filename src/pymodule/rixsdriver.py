##
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

from .veloxchemlib import hartree_in_ev, bohr_in_angstrom, mpi_master
from .oneeints import compute_electric_dipole_integrals
from .subcommunicators import SubCommunicators
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical

from .inputparser import (parse_input, print_keywords, print_attributes,
                          get_random_string_parallel)

class RixsDriver:
    """
    Implements the RIXS driver.
    Note: includes only the inelastic, resonant term in
    the 2nd order pert. expression.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - photon_energy: Incoming photon energy; omega (a.u.)
        - theta: Angle (rad)
        - gamma: Life-time broadening (FWHM) (a.u.)
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes the RIXS driver.
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
        self.photon_energy = None
        self.theta = 0
        # a.u., FWHM
        self.gamma = .124/hartree_in_ev() 
        self.gamma_hwhm = self.gamma / 2

        # input keywords
        self.input_keywords = {
            'rixs': {
                'theta': ('float', 'angle between incident polarization vector and propagation vector of outgoing'),
                'gamma': ('float', 'broadening term (FWHM)'),
                'photon_energy': ('float', 'incoming photon energy'),
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

    def print_info(self):
        """
        Prints relevant orbital and state information
        """
        if self.photon_energy is None:
            print('Compute first!')
        else:
            info = self.orb_and_state_dict
            intermed_states = info['num_intermediate_states']
            final_states = info['num_final_states']
            mo_c_ind = info['mo_core_indices']
            mo_val_ind = info['mo_val_indices']
            mo_vir_ind = info['mo_vir_indices']
            ce_states = info['core_states']
            ve_states = info['val_states']
            print(f'\nNumber of states: (intermediate, final): ({intermed_states}, {final_states})')
            print(f'\nMO indices (core, valence, virtual): ({mo_c_ind}, {mo_val_ind}, {mo_vir_ind})')
            print(f'\nState indices (intermediate, final): ({ce_states}, {ve_states})')
    
    def scattering_amplitude_tensor(self, omega, core_eigenvalue, val_eigenvalue,
                                    intermediate_tdens, final_tdens, dipole_integrals, elastic=False):
        """
        The RIXS scattering amptlitude (sum-over-states transition amplitude).
        
        :param omega:
            Energy of incoming photon 
        :param f: 
            The final valence states to end
        :param eigenvalue: 
            Energy of intermediate/core excited state
        :param intermediate_tdens:
            Transition density from ground to
            intermediate/core excited state (AO)
        :param final_tdens:
            Transition density matrix from 
            intermediate to final excited state (AO)
        :param dipole_integrals:
            Electric dipole integrals (length gauge)
            in AO-basis
        
        :return: 
            The scattering amplitude tensor; shape: (3,3)
        """
        
        e_n = 1 / (omega - (core_eigenvalue + 1j * self.gamma_hwhm))

        if elastic:
            core_eigvals2 = core_eigenvalue**2
            scatt_amp = core_eigvals2 * e_n * np.einsum('xij, ij, yab, ab -> xy', dipole_integrals, intermediate_tdens, dipole_integrals, intermediate_tdens, optimize='greedy')

        else:
            omega_product = (val_eigenvalue - core_eigenvalue) * core_eigenvalue
            scatt_amp = e_n * omega_product * np.einsum('xij, ij, yab, ab -> xy', dipole_integrals, final_tdens, dipole_integrals, intermediate_tdens, optimize='greedy')

        return scatt_amp

    def cross_section(self, F, omegaprime_omega=1):
        """
        Computes the RIXS cross-section.
        
        :param F:
            The scattering amplitude tensor.
        :param omegaprime_omega:
            The ratio between outgoing- and incoming frequencies, w'/w.
            
        :return:
            The scattering cross-section.
        """
        
        # F_xy (F^(yx))^*
        FF_T_conj = np.sum(F * F.T.conjugate())
        # F_xx (F^(yy))^*
        trace_F2 = np.abs(np.trace(F))**2
        # F_xy (F^(xy))^*
        F2 = np.sum(np.abs(F)**2)

        sigma_f = omegaprime_omega * (1.0/15.0) * (
                    (2.0 - 0.5*np.sin(self.theta)**2) * F2 +
                    (0.75*np.sin(self.theta)**2 - 0.5) * (FF_T_conj + trace_F2))
        
        return sigma_f
    
    def compute(self, molecule, basis, scf_tensors, 
                rsp_tensors, cvs_rsp_tensors=None,
                fulldiag_thresh=None, tda=False, 
                final_cutoff=None, core_cutoff=None):
        """
        Computes RIXS properties.
        
        :param molecule:
            The molecule object
        :param basis:
            The atomic orbital basis object
        :scf_tensors:
            The converged SCF results
        :rsp_tensors:
            The linear-response results
        :cvs_rsp_tensors:
            The core-valence-separated linear-response results
        :num_core_orbitals:
            Defines the split between intermediate and final states
            (not yet implemented)
        NOTE: To run full diagonalization, do a subspace-restricted calculation
                with the full space, indicating which orbitals define the core.

        :return:
            Dictionary with cross-sections, outgoing photon energy in (energy loss and full energy (emission)),
            and the scattering amplitude tensor.
        """

        norb = scf_tensors['C_alpha'].shape[0]
        nocc = molecule.number_of_alpha_electrons()
        nvir = norb - nocc

        self.twoshot = False
        self.tda = tda
        
        num_vir_orbitals  = rsp_tensors['num_vir']

        if cvs_rsp_tensors is not None:
            self.twoshot = True
            # TODO: improve the handling of these types
            num_core_orbitals = cvs_rsp_tensors['num_core']
            num_val_orbitals  = nocc - num_core_orbitals # self.get_num_val_orbs(rsp_tensors) # = nocc

            num_intermediate_states = len(cvs_rsp_tensors['eigenvalues'])
            self.ostream.print_info(f'Running RIXS in the two-shot approach with {num_intermediate_states}'
                                    ' intermediate states')
            num_final_states = len(rsp_tensors['eigenvalues'])
            num_tot_states = num_final_states + num_intermediate_states

            core_states = list(range(num_intermediate_states))
            occupied_core = num_core_orbitals
            val_states = list(range(num_final_states))

        else:
            # What if too few states to include all valence orbitals?
            #self.ostream.print_info(
            #    'Assuming subspace-restricted approximation.')
            # TODO: add assertion for splitting of response vector
            #assert_msg_critical(
            #   np.any(rsp_tensors['excitation_details'] core) is not None,
            #   '')

            num_core_orbitals = np.squeeze(rsp_tensors['num_core'])
            num_val_orbitals  = np.squeeze(rsp_tensors['num_val'])

            num_tot_states = len(rsp_tensors['eigenvalues'])
            if fulldiag_thresh is not None:
                # full matrix
                if self.photon_energy is None:
                    # assume first "non-zero" resonance, i.e., first with large osc. strength, if not given
                    for k, entries in enumerate(rsp_tensors['excitation_details']):
                        if entries[0].split()[0].startswith("core") and rsp_tensors['oscillator_strengths'][k] > 1e-3:
                            self.photon_energy = [rsp_tensors['eigenvalues'][k]]
                            break

                #mask = (np.abs(rsp_tensors['eigenvalues'] - self.photon_energy) < fulldiag_thresh)
                detuning = rsp_tensors['eigenvalues'] - self.photon_energy
                mask = (detuning <= fulldiag_thresh) & (detuning >= 0)

                init_core_states = np.where(mask)[0]

                core_states = []
                """
                for state in init_core_states:
                    entry = rsp_tensors['excitation_details'][state][0].split()
                    label = entry[0]
                    if label.startswith("core"):
                        core_states.append(state)
                """
                if core_cutoff:
                    init_core_states = np.where(mask)[0]

                    core_states = []
                    for state in init_core_states:
                        entry = rsp_tensors['excitation_details'][state][0].split()
                        label = entry[0]
                        if label.startswith("core"):
                            core_states.append(state)

                    #try:
                    #    core_states = np.squeeze(core_states)[:core_cutoff]
                    #except:
                core_states = np.squeeze(core_states)
                
                #core_states = np.where(mask)[0]
                num_intermediate_states = len(core_states)
                val_states  = np.where(~mask)[0]
                num_final_states = len(val_states)
                if final_cutoff is not None:
                    num_final_states = final_cutoff
                    val_states = val_states[:final_cutoff]

            else:
                # normal SRA
                num_final_states = num_val_orbitals * num_vir_orbitals
                # do something with full_diag based on ene_thresh
                core_states = list(range(num_final_states, num_tot_states))
                val_states = list(range(num_final_states))
                num_intermediate_states = num_tot_states - num_final_states
            assert_msg_critical(num_intermediate_states > 0,
                                 'Too few excited states included in response calculation')
            cvs_rsp_tensors = rsp_tensors
            occupied_core = num_core_orbitals + num_val_orbitals
        
        mo_core_indices = list(range(num_core_orbitals))
        mo_val_indices  = list(range(nocc - num_val_orbitals, nocc))
        mo_vir_indices  = list(range(nocc, nocc + num_vir_orbitals))
        
        mo_occ = scf_tensors['C_alpha'][:, mo_core_indices + mo_val_indices]
        mo_vir = scf_tensors['C_alpha'][:, mo_vir_indices]

        core_eigvals    = cvs_rsp_tensors['eigenvalues'][core_states]
        valence_eigvals = rsp_tensors['eigenvalues'][val_states]

        try:
            core_eigvecs    = np.array([
                self.get_full_solution_vector(cvs_rsp_tensors['eigenvectors_distributed'][state]) for state in core_states])
            valence_eigvecs = np.array([
                self.get_full_solution_vector(rsp_tensors['eigenvectors_distributed'][state]) for state in val_states])
        except KeyError:
            try:
                core_eigvecs    = np.array([
                    cvs_rsp_tensors['eigenvectors'].T[state] for state in core_states])
                valence_eigvecs = np.array([
                    rsp_tensors['eigenvectors'].T[state] for state in val_states])
                self.tda = True
            except KeyError:
                core_eigvecs    = np.array([cvs_rsp_tensors['S' + str(i + 1)] for i in core_states])
                valence_eigvecs = np.array([rsp_tensors['S' + str(i + 1)] for i in val_states])

        # For bookkeeping
        self.orb_and_state_dict = {
            'num_intermediate_states': num_intermediate_states,
            'num_final_states': num_final_states,
            'mo_core_indices': mo_core_indices,
            'mo_val_indices': mo_val_indices,
            'mo_vir_indices': mo_vir_indices,
            'core_states': core_states,
            'val_states': val_states,
        }
        
        if self.photon_energy is None:
            # assume first core resonance
            if cvs_rsp_tensors is None:
                osc_arr = rsp_tensors['oscillator_strengths']
            else:
                osc_arr = cvs_rsp_tensors['oscillator_strengths']

            for k, osc in enumerate(osc_arr[core_states]):
                if osc > 1e-3:
                    self.photon_energy = [core_eigvals[k]]
                    self.ostream.print_info(
                        'Incoming photon energy not set; computing ' \
                        'only for the first core resonance at: ' \
                        f'{self.photon_energy[0]:.4f} a.u. = ' \
                        f'{self.photon_energy[0] * hartree_in_ev():.2f} eV')
                    break
        elif isinstance(self.photon_energy, (float, int, np.floating)):
            self.ostream.print_info(
                f'Incoming photon energy: {self.photon_energy:.2f} a.u. = ' \
                f'{self.photon_energy*hartree_in_ev():.2f} eV')
            self.photon_energy = [self.photon_energy]
        else:
            formatted_au = ', '.join(f'{enes:.2f}' for enes in self.photon_energy)
            formatted_ev = ', '.join(f'{enes * hartree_in_ev():.2f}' for enes in self.photon_energy)
            self.ostream.print_info(f'Incoming photon energies: ({formatted_au}) a.u. = ({formatted_ev}) eV')

        dipole_integrals = compute_electric_dipole_integrals(
                molecule, basis, [0.0,0.0,0.0])

        ene_losses = np.zeros((num_final_states, len(self.photon_energy)))
        emission_enes = np.zeros((num_final_states, len(self.photon_energy)))
        cross_sections = np.zeros((num_final_states, len(self.photon_energy)))
        elastic_cross_sections = np.zeros((len(self.photon_energy)))
        scattering_amplitudes = np.zeros((num_final_states, len(self.photon_energy),
                                           3, 3), dtype=np.complex128)

        for w_ind, omega in enumerate(self.photon_energy):
            F_elastic = np.zeros((3,3), dtype=np.complex128)
            for f in range(num_final_states):
                F_inelastic = np.zeros((3,3), dtype=np.complex128)

                if self.tda:
                    valence_z_mat = valence_eigvecs[f].reshape(num_core_orbitals + num_val_orbitals, num_vir_orbitals)
                    valence_y_mat = np.zeros_like(valence_z_mat)
                else:
                    half_val = valence_eigvecs.shape[1] // 2
                    valence_z_mat = valence_eigvecs[f, :half_val].reshape(
                                                    num_core_orbitals + num_val_orbitals, num_vir_orbitals)
                    valence_y_mat = valence_eigvecs[f, half_val:].reshape(
                                                    num_core_orbitals + num_val_orbitals, num_vir_orbitals)
                    
                for n in range(num_intermediate_states):
                    
                    if self.tda:
                        core_z_mat = core_eigvecs[n].reshape(occupied_core, num_vir_orbitals)
                        core_y_mat = np.zeros_like(core_z_mat)
                    else:
                        half_core = core_eigvecs.shape[1] // 2
                        core_z_mat = core_eigvecs[n, :half_core].reshape(occupied_core, num_vir_orbitals)
                        core_y_mat = core_eigvecs[n, half_core:].reshape(occupied_core, num_vir_orbitals)

                    if self.twoshot:
                        # pad with zeroes
                        padding = ((0, num_val_orbitals), (0, 0))

                        core_z_mat = np.pad(core_z_mat, padding)
                        core_y_mat = np.pad(core_y_mat, padding)
 
                    gs_to_core_tdens = mo_occ @ (core_z_mat - core_y_mat) @ mo_vir.T  
                    gs_to_core_tdens *= np.sqrt(2)

                    core_to_val_tdens = (
                                np.linalg.multi_dot(
                                    [mo_vir, valence_z_mat.T, core_z_mat, mo_vir.T]) -
                                np.linalg.multi_dot(
                                    [mo_occ, valence_z_mat, core_z_mat.T, mo_occ.T]))

                    core_to_val_tdens += (
                                np.linalg.multi_dot(
                                    [mo_occ, valence_y_mat, core_y_mat.T, mo_occ.T]) -
                                np.linalg.multi_dot(
                                    [mo_vir, valence_y_mat.T, core_y_mat, mo_vir.T]))
                    
                    F_inelastic += self.scattering_amplitude_tensor(omega, core_eigvals[n], valence_eigvals[f], gs_to_core_tdens,
                                                            core_to_val_tdens, dipole_integrals)
                    F_elastic  += self.scattering_amplitude_tensor(omega, core_eigvals[n], None, gs_to_core_tdens,
                                                            None, dipole_integrals, elastic=True)
                
                # TODO: improve results dictionary structure
                emission_enes[f, w_ind] = omega - valence_eigvals[f]
                ene_losses[f, w_ind] = valence_eigvals[f]
                prefactor_ratio = emission_enes[f, w_ind] / omega # w'/w

                sigma = self.cross_section(F_inelastic, prefactor_ratio)
                sigma_elastic = self.cross_section(F_elastic)

                elastic_cross_sections[w_ind] = sigma_elastic.real
                scattering_amplitudes[f, w_ind] = F_inelastic
                cross_sections[f, w_ind] = sigma.real
        
        return_dict = {
                    'cross_sections': cross_sections,
                    'elastic_cross_sections': elastic_cross_sections,
                    'elastic_emission': self.photon_energy,
                    'scattering_amplitudes': scattering_amplitudes,
                    'emission_energies': emission_enes,
                    'energy_losses': ene_losses,
                    'excitation_energies': core_eigvals,
                    }

        return return_dict
    
    @staticmethod
    def get_full_solution_vector(solution):
        """
        Gets a full solution vector from the distributed solution.

        :param solution:
            The distributed solution as a tuple.

        :return:
            The full solution vector.
        """

        x_ger = solution.get_full_vector(0)
        x_ung = solution.get_full_vector(1)

        if solution.rank == mpi_master():
            x_ger_full = np.hstack((x_ger, x_ger))
            x_ung_full = np.hstack((x_ung, -x_ung))
            return x_ger_full + x_ung_full
        else:
            return None


