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
        - theta: Angle (rad.)
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

        self.tda = False
        self.fulldiag_threshold = np.inf
        self.num_final_states = None
        self.core_cutoff = None

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
        The RIXS scattering amptlitude (sum-over-states) transition amplitude.
        
        :param omega:
            Energy of incoming photon.
        :param f: 
            The final valence state index.
        :param eigenvalue: 
            Energy of intermediate/core excited state.
        :param intermediate_tdens:
            Transition density from ground to
            intermediate/core excited state (AO).
        :param final_tdens:
            Transition density matrix from 
            intermediate to final excited state (AO).
        :param dipole_integrals:
            Electric dipole integrals (length gauge)
            in AO-basis.
        :param elastic:
            Bool to compute the elastic line.

        :return: 
            The scattering amplitude tensor: shape = (3,3)
        """
        
        e_n = 1 / (omega - (core_eigenvalue + 1j * self.gamma_hwhm))

        if elastic:
            core_eigvals2 = core_eigenvalue**2
            gs_exc_tdipole = np.array([
                                np.sum(intermediate_tdens * dipole_integrals[i])
                                for i in range(3)])
            outer_product = np.matmul(gs_exc_tdipole[:, np.newaxis], gs_exc_tdipole[np.newaxis, :])
            scatt_amp = e_n * core_eigvals2 * outer_product

        else:
            omega_product = (val_eigenvalue - core_eigenvalue) * core_eigenvalue
            gs_exc_tdipole = np.array([
                                np.sum(intermediate_tdens * dipole_integrals[i])
                                for i in range(3)])
            exc_exc_tdipole = np.array([
                                np.sum(final_tdens * dipole_integrals[i])
                                for i in range(3)])
            outer_product = np.matmul(exc_exc_tdipole[:, np.newaxis], gs_exc_tdipole[np.newaxis, :])
            scatt_amp = e_n * omega_product * outer_product

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
                rsp_tensors, cvs_rsp_tensors=None):
        """
        Computes RIXS properties.
        
        :param molecule:
            The molecule object.
        :param basis:
            The AO basis set.
        :scf_tensors:
            The dictionary of tensors from converged SCF 
            wavefunction.
        :rsp_tensors:
            The linear-response results dictionary.
        :cvs_rsp_tensors:
            The core-valence-separated linear-response 
            results dictionary
        NOTE: To run full diagonalization, do a restricted
              subspace calculation with the full space, 
              indicating which orbitals define the core.

        :return:
            Dictionary with cross-sections, outgoing photon energy
            in (energy loss and full energy (emission)),
            and the scattering amplitude tensor.
        """

        norb = scf_tensors['C_alpha'].shape[0]
        nocc = molecule.number_of_alpha_electrons()

        self.twoshot = False
        
        num_vir_orbitals  = rsp_tensors['num_vir']

        if cvs_rsp_tensors is not None:
            self.twoshot = True
            num_core_orbitals = cvs_rsp_tensors['num_core']
            num_val_orbitals  = nocc - num_core_orbitals 

            num_intermediate_states = len(cvs_rsp_tensors['eigenvalues'])
            num_final_states = len(rsp_tensors['eigenvalues'])
            self.ostream.print_info('Running RIXS in the two-shot approach '\
                                    f'with {num_intermediate_states} intermediate states')

            occupied_core = num_core_orbitals
            core_states = list(range(num_intermediate_states))
            val_states = list(range(num_final_states))

        else:
            num_core_orbitals = rsp_tensors['num_core']
            assert_msg_critical(num_core_orbitals > 0,
                                 'No core orbitals indicated in the response tensor.')
            num_val_orbitals  = rsp_tensors['num_val']


            for k, entries in enumerate(rsp_tensors['excitation_details']):
                if entries[0].split()[0].startswith("core"):
                    first_core_ene = rsp_tensors['eigenvalues'][k]
                    break

            detuning = rsp_tensors['eigenvalues'] - first_core_ene
            tol = 1e-10
            mask = (detuning >= -tol) & (detuning <= self.fulldiag_threshold)
            #mask = (detuning >= 0) & (detuning <= self.fulldiag_threshold)
            init_core_states = np.where(mask)[0]
            core_states = []
            for state in init_core_states:
                entry = rsp_tensors['excitation_details'][state][0].split()
                label = entry[0]
                if label.startswith("core"):
                    core_states.append(state)

            num_intermediate_states = len(core_states)
            assert_msg_critical(num_intermediate_states > 0,
                                 'Too few excited states included in response calculation.')

            val_states = np.where(detuning < 0)[0]
            num_final_states = len(val_states)
            if self.num_final_states is not None:
                num_final_states = self.num_final_states
                val_states = val_states[:self.num_final_states]

            cvs_rsp_tensors = rsp_tensors
            occupied_core = num_core_orbitals + num_val_orbitals
        
        mo_core_indices = list(range(num_core_orbitals))
        mo_val_indices  = list(range(nocc - num_val_orbitals, nocc))
        mo_vir_indices  = list(range(nocc, nocc + num_vir_orbitals))
        
        mo_occ = scf_tensors['C_alpha'][:, mo_core_indices + mo_val_indices]
        mo_vir = scf_tensors['C_alpha'][:, mo_vir_indices]

        core_eigvals    = cvs_rsp_tensors['eigenvalues'][core_states]
        valence_eigvals = rsp_tensors['eigenvalues'][val_states]
        
        core_eigvecs    = self.get_eigvecs(cvs_rsp_tensors, core_states, "core")
        valence_eigvecs = self.get_eigvecs(rsp_tensors, val_states, "valence")

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
            for eig, osc in zip(core_eigvals, osc_arr[core_states]):
                if osc > 1e-3:
                    self.photon_energy = [eig]
                    self.ostream.print_info(
                        'Incoming photon energy not set; computing ' \
                        'RIXS for the first core resonance at: ' \
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

        ene_losses             = np.zeros((num_final_states, len(self.photon_energy)))
        emission_enes          = np.zeros((num_final_states, len(self.photon_energy)))
        cross_sections         = np.zeros((num_final_states, len(self.photon_energy)))
        elastic_cross_sections = np.zeros((len(self.photon_energy)))
        scattering_amplitudes  = np.zeros((num_final_states, len(self.photon_energy),
                                           3, 3), dtype=complex)

        for w_ind, omega in enumerate(self.photon_energy):
            F_elastic = np.zeros((3,3), dtype=complex)

            for f in range(num_final_states):
                F_inelastic = np.zeros((3,3), dtype=complex)
                z_val, y_val = self.split_eigvec(valence_eigvecs[f], num_core_orbitals + num_val_orbitals,
                                            num_vir_orbitals, self.tda)
                    
                for n in range(num_intermediate_states):
                    z_core, y_core = self.split_eigvec(core_eigvecs[n], occupied_core,
                                            num_vir_orbitals, self.tda)
                    if self.twoshot:
                        z_core, y_core = self.pad_if_twoshot(z_core, y_core, num_val_orbitals)

                    gs2core, core2val = self.build_tdms(mo_occ, mo_vir, z_val, y_val, z_core, y_core)
                    
                    F_inelastic += self.scattering_amplitude_tensor(omega, core_eigvals[n], valence_eigvals[f],
                                                                    gs2core, core2val, dipole_integrals)
                    F_elastic   += self.scattering_amplitude_tensor(omega, core_eigvals[n], None, gs2core,
                                                                    None, dipole_integrals, elastic=True)
                

                emission_enes[f, w_ind] = omega - valence_eigvals[f]
                ene_losses[f, w_ind] = valence_eigvals[f]
                # w'/w
                prefactor_ratio = emission_enes[f, w_ind] / omega

                sigma = self.cross_section(F_inelastic, prefactor_ratio)
                sigma_elastic = self.cross_section(F_elastic)

                elastic_cross_sections[w_ind]   = sigma_elastic.real
                cross_sections[f, w_ind] = sigma.real
                scattering_amplitudes[f, w_ind] = F_inelastic
        
        results_dict = {
                    'cross_sections': cross_sections,
                    'elastic_cross_sections': elastic_cross_sections,
                    'elastic_emission': self.photon_energy,
                    'scattering_amplitudes': scattering_amplitudes,
                    'emission_energies': emission_enes,
                    'energy_losses': ene_losses,
                    'excitation_energies': core_eigvals,
                    }

        return results_dict
    
    def get_eigvecs(self, tensor, states, label):
        if 'eigenvectors_distributed' in tensor:
            return np.array([
                self.get_full_solution_vector(tensor['eigenvectors_distributed'][i]) for i in states])
        elif 'eigenvectors' in tensor:
            self.tda = True
            return np.array([
                tensor['eigenvectors'].T[i] for i in states])
        else:
            # h5-file as input
            try:
                return np.array([tensor[f'S{i + 1}'] for i in states])
            except KeyError as e:
                raise RuntimeError(f"Could not extract {label} eigenvectors using the key: {e}")
            
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
    
    @staticmethod
    def split_eigvec(eigvec, n_occ, n_vir, tda=False):
        if tda:
            z = eigvec.reshape(n_occ, n_vir)
            y = np.zeros_like(z)
        else:
            half = eigvec.shape[0] // 2
            z = eigvec[:half].reshape(n_occ, n_vir)
            y = eigvec[half:].reshape(n_occ, n_vir)
        return z, y
    
    @staticmethod
    def pad_if_twoshot(z_mat, y_mat, pad_occ):
        padding = ((0, pad_occ), (0, 0))
        return np.pad(z_mat, padding), np.pad(y_mat, padding)
    
    @staticmethod
    def build_tdms(mo_occ, mo_vir, z_val, y_val, z_core, y_core):
        gs_to_core = mo_occ @ (z_core - y_core) @ mo_vir.T
        gs_to_core *= np.sqrt(2)

        core_to_val = (
            np.linalg.multi_dot([mo_vir, z_val.T, z_core, mo_vir.T]) -
            np.linalg.multi_dot([mo_occ, z_val, z_core.T, mo_occ.T]) +
            np.linalg.multi_dot([mo_occ, y_val, y_core.T, mo_occ.T]) -
            np.linalg.multi_dot([mo_vir, y_val.T, y_core, mo_vir.T])
        )
        return gs_to_core, core_to_val


