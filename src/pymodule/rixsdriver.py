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
        self.gamma = .124/hartree_in_ev() # a.u.

        # input keywords
        self.input_keywords = {
            'rixs': {
                'theta': ('float', 'angle between incident polarization vector and propagation vector of outgoing'),
                'gamma': ('float', 'broadening term (FWHM)'),
                #'nr_CO': ('int', 'number of involved core-orbitals'),
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
    
    def scattering_amplitude_tensor(self, omega, core_eigenvalues, val_eigenvalues,
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
        gamma_hwhm = self.gamma / 2
        e_n = 1 / (omega - (core_eigenvalues + 1j * gamma_hwhm))
        omega_product = (val_eigenvalues[:, np.newaxis] - core_eigenvalues) * core_eigenvalues

        if elastic:
            core_eigvals_2 = core_eigenvalues**2
            scatt_amp = np.einsum('n, xij, ijn, yab, abn -> nxy', core_eigvals_2 * e_n, dipole_integrals, intermediate_tdens, dipole_integrals, intermediate_tdens, optimize='greedy')

        else:
            #scatt_amp = np.einsum('n, xij, ijfn, yab, abn -> fxy', e_n, dipole_integrals, final_tdens, dipole_integrals, intermediate_tdens, optimize='greedy')
            scatt_amp = np.einsum('n, fn, xij, ijfn, yab, abn -> fxy', e_n, omega_product, dipole_integrals, final_tdens, dipole_integrals, intermediate_tdens, optimize='greedy')

        return scatt_amp

    def cross_section(self, F, prefactor_ratio=1):
        """
        Computes the cross section
        
        :param F:
            The scattering amplitude tensor

        :return:
            The scattering cross section
        """
        FF_T_conj = np.sum(F * F.transpose(0,2,1).conjugate(), axis=(1, 2))
        trace_F_2 = np.sum(np.diagonal(np.abs(F)**2, axis1=1, axis2=2), axis=-1)

        # Combine into cross section for each final state
        sigma_f = prefactor_ratio * (1.0 / 15.0) * (
            (2.0 - 0.5 * np.sin(self.theta)**2) * np.sum(np.abs(F)**2, axis=(1, 2)) +
            (0.75 * np.sin(self.theta)**2 - 0.5) * (FF_T_conj + trace_F_2))

        #sigma = 1/15 * ((2 - (1/2) * np.sin(self.theta) ** 2) * np.sum(np.abs(F)**2) 
        #           + ((3/4) * np.sin(self.theta) ** 2 - 1/2) * (np.sum(F * F.T.conj())
        #                                           + np.trace(np.abs(F)**2)))
        return sigma_f
    
    def compute(self, molecule, basis, scf_tensors, 
                rsp_tensors, cvs_rsp_tensors=None,
                fulldiag_thresh=None, tda=False, final_cutoff=None):
        """
        Computes the relevant RIXS properties
        
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
                  with the full space, indicating which orbitals define the core

        :return:
            Dictionary with outgoing photon energy in energy loss and full energy (emisison)
            as well as scattering amplitudes and cross-sections
        """

        norb = scf_tensors['C_alpha'].shape[0]
        nocc = molecule.number_of_alpha_electrons()
        nvir = norb - nocc

        self.cvs = False
        #self.sra = False
        self.tda = tda
        
        num_vir_orbitals  = np.squeeze(rsp_tensors['num_vir'])
        # TODO: add safeguards for when both num_intermediate_states and cvs is given

        if cvs_rsp_tensors is not None:
            #self.ostream.print_info(
            #    'Running RIXS with CVS approximation.')
            self.cvs = True
            # TODO: improve the handling of these types
            #num_core_orbitals = cvs_rsp_tensors['num_core']
            num_core_orbitals = np.squeeze(cvs_rsp_tensors['num_core'])
            num_val_orbitals  = nocc - num_core_orbitals # self.get_num_val_orbs(rsp_tensors) # = nocc

            num_intermediate_states = len(cvs_rsp_tensors['eigenvalues'])
            num_final_states = len(rsp_tensors['eigenvalues'])
            num_tot_states = num_final_states + num_intermediate_states

            core_states = list(range(num_intermediate_states))
            occupied_core = num_core_orbitals
            val_states = list(range(num_final_states))

        else:
            # What if too few states to include all valence orbitals?
            self.ostream.print_info(
                'Assuming subspace-restricted approximation.')
            # TODO: add assertion for splitting of response vector
            #assert_msg_critical(
            #   np.any(rsp_tensors['excitation_details'] core) is not None,
            #   '')
            #self.sra = True
            num_core_orbitals = np.squeeze(rsp_tensors['num_core'])
            num_val_orbitals  = np.squeeze(rsp_tensors['num_val'])

            num_tot_states = len(rsp_tensors['eigenvalues'])
            if num_core_orbitals == 0 or fulldiag_thresh is not None:
                # full diagonalization
                mask = (np.abs(rsp_tensors['eigenvalues'] - self.photon_energy) < fulldiag_thresh)
                # perhaps remove all states from final_states with eigenvalues LARGER than core/intermediate?
                # otherwise we get negative emission energies
                num_final_states = len(rsp_tensors['eigenvalues'][~mask])
                num_intermediate_states = len(rsp_tensors['eigenvalues'][mask])
                core_states = np.where(mask)[0]
                val_states  = np.where(~mask)[0]
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
            #self.ostream.print_info(
            #    'Incoming photon energy not set; calculating only for the first core resonance.')
            self.photon_energy = [core_eigvals[0]]
        elif type(self.photon_energy) == float:
            self.photon_energy = [self.photon_energy]

        dipole_integrals = compute_electric_dipole_integrals(
                molecule, basis, [0.0,0.0,0.0])

        ene_losses = np.zeros((num_final_states, len(self.photon_energy)))
        emission_enes = np.zeros((num_final_states, len(self.photon_energy)))
        cross_sections = np.zeros((num_final_states, len(self.photon_energy)))
        elastic_cross_sections = np.zeros((num_intermediate_states, len(self.photon_energy)))
        scattering_amplitudes = np.zeros((num_final_states, len(self.photon_energy),
                                           3, 3), dtype=np.complex128)

        for w_ind, omega in enumerate(self.photon_energy):
                    
            if self.tda:
                valence_z_mat = valence_eigvecs.reshape(num_final_states,
                            num_core_orbitals + num_val_orbitals, num_vir_orbitals)
                valence_y_mat = np.zeros_like(valence_z_mat)
    
                core_z_mat = core_eigvecs.reshape(num_intermediate_states, occupied_core, num_vir_orbitals)
                core_y_mat = np.zeros_like(core_z_mat)

            else:
                half_val = valence_eigvecs.shape[1] // 2
                valence_z_mat = valence_eigvecs[:, :half_val].reshape(num_final_states,
                                                num_core_orbitals + num_val_orbitals, num_vir_orbitals)
                valence_y_mat = valence_eigvecs[:, half_val:].reshape(num_final_states,
                                                num_core_orbitals + num_val_orbitals, num_vir_orbitals)
            
                half_core = core_eigvecs.shape[1] // 2
                core_z_mat = core_eigvecs[:, :half_core].reshape(num_intermediate_states,
                                                                 occupied_core, num_vir_orbitals)
                core_y_mat = core_eigvecs[:, half_core:].reshape(num_intermediate_states,
                                                                 occupied_core, num_vir_orbitals)

            if self.cvs:
                new_shape = (
                    num_intermediate_states,
                    occupied_core + num_val_orbitals,
                    core_z_mat.shape[2],
                )
                tmp_z = np.zeros(new_shape, dtype=core_z_mat.dtype)
                tmp_y = np.zeros(new_shape, dtype=core_y_mat.dtype)

                tmp_z[:, :occupied_core, :] = core_z_mat
                tmp_y[:, :occupied_core, :] = core_y_mat

                core_z_mat = tmp_z
                core_y_mat = tmp_y

            mo_occ_expanded  = mo_occ[np.newaxis, :, :]
            mo_vir_expanded = mo_vir[:, :, np.newaxis] 

            # want a + sign here but
            difference_mat  = core_z_mat - core_y_mat
            gs_core = mo_occ_expanded @ difference_mat @ mo_vir_expanded.T  

            # rearrange to (norb, norb, num_intermediate_states)
            gs_to_core_tdens = np.sqrt(2) * gs_core.transpose(1, 2, 0)

            core_to_val_tdens  = np.einsum('ui, ijf, njk, kv -> uvfn', mo_vir, valence_z_mat.T, core_z_mat, mo_vir.T, optimize=True)
            core_to_val_tdens -= np.einsum('ui, fij, jkn, kv -> uvfn', mo_occ, valence_z_mat, core_z_mat.T, mo_occ.T, optimize=True)
            core_to_val_tdens += np.einsum('ui, fij, jkn, kv -> uvfn', mo_occ, valence_y_mat, core_y_mat.T, mo_occ.T, optimize=True)
            core_to_val_tdens -= np.einsum('ui, ijf, njk, kv -> uvfn', mo_vir, valence_y_mat.T, core_y_mat, mo_vir.T, optimize=True)
                
            # TODO: improve results dictionary structure
            emission_ene = omega - valence_eigvals #[f]
            emission_enes[..., w_ind] = emission_ene
            prefactor_ratio = emission_ene / omega #w'/w
            energy_loss = valence_eigvals #omega - emission_ene # independent of intermediate state, but 
            ene_losses[..., w_ind] = energy_loss # keeping it like this for consistency
            
            F = self.scattering_amplitude_tensor(omega, core_eigvals, valence_eigvals, gs_to_core_tdens,
                                                    core_to_val_tdens, dipole_integrals)
            F_elastic = self.scattering_amplitude_tensor(omega, core_eigvals, valence_eigvals, gs_to_core_tdens,
                                                    core_to_val_tdens, dipole_integrals, elastic=True)
            
            sigma = self.cross_section(F, prefactor_ratio)
            sigma_elastic = self.cross_section(F_elastic)

            scattering_amplitudes[:, w_ind, ...] = F
            cross_sections[..., w_ind] = sigma.real
            elastic_cross_sections[..., w_ind] = sigma_elastic.real
        
        return_dict = {
                    'cross_sections': cross_sections,
                    'elastic_cross_sections': elastic_cross_sections,
                    'elastic_emission': self.photon_energy,
                    'scattering_amplitudes': scattering_amplitudes,
                    'emission_energies': emission_enes,
                    'energy_losses': ene_losses}

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


