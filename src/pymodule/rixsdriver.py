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
        - gamma_n: Life-time broadening (a.u.)
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
        self.gamma_n = .124/hartree_in_ev() # a.u.

        # input keywords
        self.input_keywords = {
            'rixs': {
                'theta': ('float', 'angle between incident polarization vector and propagation vector of outgoing'),
                'gamma': ('float', 'broadening term'),
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

    def scattering_amplitude_tensor(self, omega, eigenvalue, intermediate_tdens,
                             final_tdens, dipole_integrals, elastic=False):
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

        e_n = 1 / (omega - (eigenvalue + 1j * self.gamma_n))

        if elastic:
            # TODO: addd missing term (F.T_f) in eq.7 
            # ref. Gel'mukhanov, F., & Ågren, H. (1999). Physics Reports, 312(3-6), 87-330.
            scatt_amp = np.einsum('n, xij, ijn, yab, abn -> xyn', e_n, dipole_integrals, intermediate_tdens, dipole_integrals, intermediate_tdens, optimize='greedy')

        else:
            scatt_amp = np.einsum('n, xij, ijn, yab, abn -> xy', e_n, dipole_integrals, final_tdens, dipole_integrals, intermediate_tdens, optimize='greedy')

        return scatt_amp

    def cross_section(self, F):
        """
        Computes the cross-section
        
        :param F:
            The scattering amplitude tensor

        :return:
            The scattering cross-section
        """
        sigma = 1/15 * ((2 - (1/2) * np.sin(self.theta) ** 2) * np.sum(np.abs(F)**2) 
                   + ((3/4) * np.sin(self.theta) ** 2 - 1/2) * (np.sum(F * F.T.conj())
                                                   + np.trace(np.abs(F)**2)))
        return sigma
    
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
    
    def compute(self, molecule, basis, scf_tensors, rsp_tensors,
                 cvs_rsp_tensors=None, num_core_orbitals=None):
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
        self.sra = False

        # TODO: add safeguards for when both num_intermediate_states and cvs is given
        # TODO: think if there is ever a case where get_num_val_orbs() don't find all valence?

        if num_core_orbitals is not None:
            self.ostream.print_info(
                'Full space eigenvector assumed for RIXS. Number of intermediate/source/core orbitals {:d}'.format(
                    num_core_orbitals))
            # use intermediate states to define the split b/w core and valence
            num_val_orbitals = nocc - num_core_orbitals
            #num_val_orbitals = nocc
            num_vir_orbitals = nvir # =nvir

            num_tot_states = len(rsp_tensors['eigenvalues'])
            num_final_states = num_val_orbitals * num_vir_orbitals
            num_intermediate_states = num_tot_states - num_final_states

        elif cvs_rsp_tensors is not None:
            self.ostream.print_info(
                'Running RIXS with CVS approximation.')
            self.cvs = True
            num_core_orbitals = self.get_num_core_orbs(cvs_rsp_tensors)
            num_val_orbitals  = nocc - num_core_orbitals # self.get_num_val_orbs(rsp_tensors) # = nocc
            num_vir_orbitals  = nvir # self.get_num_vir_orbs(rsp_tensors) # = nvir

            num_final_states = len(rsp_tensors['eigenvalues'])
            num_intermediate_states = len(cvs_rsp_tensors['eigenvalues'])
            num_tot_states = num_final_states + num_intermediate_states
            core_states = list(range(num_intermediate_states))
            occupied_core = num_core_orbitals

        else:
            self.ostream.print_info(
                'Assuming subspace-restricted approximation.')
            # TODO: add assertion for splitting of response vector
            #assert_msg_critical(
            #   np.any(rsp_tensors['excitation_details'] core) is not None,
            #   '')
            self.sra = True
            num_core_orbitals = self.get_num_core_orbs(rsp_tensors)
            num_val_orbitals  = self.get_num_val_orbs(rsp_tensors)
            num_vir_orbitals  = self.get_num_vir_orbs(rsp_tensors)
            
            num_tot_states = len(rsp_tensors['eigenvalues'])
            num_final_states = num_val_orbitals * num_vir_orbitals
            num_intermediate_states = num_tot_states - num_final_states
            assert_msg_critical(num_intermediate_states > 0, 'Too few excited states included in response calculation')
            core_states = list(range(num_final_states, num_tot_states))
            cvs_rsp_tensors = rsp_tensors
            occupied_core = num_core_orbitals + num_val_orbitals

        
        mo_core_indices = list(range(num_core_orbitals))
        mo_val_indices  = list(range(nocc - num_val_orbitals, nocc))
        mo_vir_indices  = list(range(nocc, nocc + num_vir_orbitals))
        mo_occ = scf_tensors['C_alpha'][:, mo_core_indices + mo_val_indices]
        mo_vir = scf_tensors['C_alpha'][:, mo_vir_indices]

        val_states = list(range(num_final_states))
        core_eigvecs_dist = np.array([cvs_rsp_tensors['eigenvectors_distributed'][state] for state in core_states])
        core_eigvals = cvs_rsp_tensors['eigenvalues'][core_states]
        valence_eigvecs_dist = np.array([rsp_tensors['eigenvectors_distributed'][state] for state in val_states])
        valence_eigvals = rsp_tensors['eigenvalues'][val_states]
        
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
            self.ostream.print_info(
                'Incoming photon energy not set; calculating only for the first core resonance.')
            self.photon_energy = [core_eigvals[0]]

        dipole_integrals = compute_electric_dipole_integrals(
                molecule, basis, [0.0,0.0,0.0])

        ene_losses = []
        emission_enes = []
        cross_sections = []
        scattering_amplitudes = []

        for w_ind, omega in enumerate(self.photon_energy):

            for f in range(num_final_states):

                valence_eigvec = self.get_full_solution_vector(valence_eigvecs_dist[f])
                    
                valence_z_mat = valence_eigvec[:len(valence_eigvec) // 2].reshape(
                            num_core_orbitals + num_val_orbitals, num_vir_orbitals)
                valence_y_mat = valence_eigvec[len(valence_eigvec) // 2:].reshape(
                            num_core_orbitals + num_val_orbitals, num_vir_orbitals)
                
                gs_to_core_tdens = np.zeros((norb, norb, num_intermediate_states))
                core_to_val_tdens = np.zeros((norb, norb, num_intermediate_states))

                for n in range(num_intermediate_states):

                    core_eigvec = self.get_full_solution_vector(core_eigvecs_dist[n])
                    
                    core_z_mat = core_eigvec[:len(core_eigvec) // 2].reshape(
                            occupied_core, num_vir_orbitals)
                    core_y_mat = core_eigvec[len(core_eigvec) // 2:].reshape(
                            occupied_core, num_vir_orbitals)

                    if self.cvs:
                        core_z_mat = np.vstack([core_z_mat, np.zeros((num_val_orbitals, core_z_mat.shape[1]))])
                        core_y_mat = np.vstack([core_y_mat, np.zeros((num_val_orbitals, core_y_mat.shape[1]))])
                   
                    # want a + sign here but
                    gs_to_core_tdens[..., n] = np.sqrt(2) * np.linalg.multi_dot([mo_occ, core_z_mat - core_y_mat, mo_vir.T])

                    core_to_val_tdens[..., n] = (
                                np.linalg.multi_dot(
                                    [mo_vir, valence_z_mat.T, core_z_mat, mo_vir.T]) -
                                np.linalg.multi_dot(
                                    [mo_occ, valence_z_mat, core_z_mat.T, mo_occ.T]))

                    core_to_val_tdens[..., n] += (
                                np.linalg.multi_dot(
                                    [mo_occ, valence_y_mat, core_y_mat.T, mo_occ.T]) -
                                np.linalg.multi_dot(
                                    [mo_vir, valence_y_mat.T, core_y_mat, mo_vir.T]))
                    
                emission_ene = omega - valence_eigvals[f]
                emission_enes.append([w_ind, n, emission_ene])
                energy_loss = omega - emission_ene
                ene_losses.append([w_ind, n, energy_loss])
                    
                F = self.scattering_amplitude_tensor(omega, core_eigvals, gs_to_core_tdens,
                                                       core_to_val_tdens, dipole_integrals)
                scattering_amplitudes.append(F)

                sigma = self.cross_section(F)
                cross_sections.append(sigma.real)
        
        return_dict = {
                    'cross_sections': cross_sections,
                    'scattering_amplitudes': scattering_amplitudes,
                    'emission_energies': emission_enes,
                    'energy_losses': ene_losses}

        return return_dict

    # TODO: find better solution for the below functions
    @staticmethod
    def get_num_core_orbs(rsp_results):
        excitation_details = rsp_results['excitation_details']
        largest_core = 0  # Initialize the largest core index as 0

        for detail in excitation_details:
            entry = detail[0]
            if 'core_' in entry:
                core_index = int(entry.split('core_')[1].split()[0])
                largest_core = max(largest_core, core_index)
        return largest_core
    
    @staticmethod
    def get_num_val_orbs(rsp_results):
        excitation_details = rsp_results['excitation_details']
        largest_valence = 0  # Initialize the largest valence index as 0

        for detail in excitation_details:
            entry = detail[0]
            if 'HOMO-' in entry:
                val_index = int(entry.split('HOMO-')[1].split()[0])
                largest_valence = max(largest_valence, val_index)
        return largest_valence + 1

    @staticmethod
    def get_num_vir_orbs(rsp_results):
        excitation_details = rsp_results['excitation_details']
        largest_virtual = 0  # Initialize the largest virtual index as 0

        for detail in excitation_details:
            entry = detail[0]
            if 'LUMO+' in entry:
                vir_index = int(entry.split('LUMO+')[1].split()[0])
                largest_virtual = max(largest_virtual, vir_index)
        return largest_virtual + 1
    
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

    


