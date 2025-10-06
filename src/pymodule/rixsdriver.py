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

        self.tda = False
        self.orb_and_state_dict = None

        #self.rixs_dict = {}
        self.filename = None

        # TODO?: define both in terms of either nr. of states or energy
        self.num_final_states = None
        self.core_cutoff = np.inf

        # input keywords
        self.input_keywords = {
            'rixs': {
                'theta': ('float', 'angle between incident polarization vector and propagation vector of outgoing'),
                'gamma': ('float', 'broadening term (FWHM)'),
                'photon_energy': ('list', 'list of incoming photon energies'),
                'num_final_states': ('int', 'number of final states to include'),
                'core_cutoff': ('float', 'energy threshold, above first core excited, for which to include core-excited states'),
            },
        }

    def update_settings(self, rixs_dict, method_dict=None):
        """
        Updates settings in RIXS driver.

        :param rixs_dict:
            The dictionary of rixs input.
        """

        if method_dict is None:
            method_dict = {}

        rixs_keywords = {
            key: val[0] for key, val in self.input_keywords['rixs'].items()
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
        if self.orb_and_state_dict is None:
            self.ostream.print_warning('Compute first!')
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
                                    intermediate_tdens, final_tdens, dipole_integrals):
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
        gamma_hwhm = self.gamma / 2
        e_n = 1 / (omega - (core_eigenvalue + 1j * gamma_hwhm))

        if val_eigenvalue is None or intermediate_tdens is None:
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
        
        # F_xy (F^(xy))^*
        F2 = np.sum(np.abs(F)**2)
        # F_xy (F^(yx))^*
        FF_T_conj = np.sum(F * F.T.conjugate())
        # F_xx (F^(yy))^*
        trace_F2 = np.abs(np.trace(F))**2

        sigma_f = omegaprime_omega * (1.0/15.0) * (
                    (2.0 - 0.5*np.sin(self.theta)**2) * F2 +
                    (0.75*np.sin(self.theta)**2 - 0.5) * (FF_T_conj + trace_F2))
        
        return sigma_f
    
    def compute(self, molecule, basis, scf_tensors, 
                rsp_tensors, cvs_rsp_tensors=None,
                cvs_scf_tensors=None):
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
            The core-valence-separated (CVS) linear-response 
            results dictionary.
        :cvs_scf_tensors:
            The dictionary of tensors from converged CVS SCF 
            wavefunction. If given -- assumes that the CVS response
            was obtained from an independent SCF wavefunction.
            
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
        init_photon_set = True
        
        num_vir_orbitals  = rsp_tensors['num_vir']

        if cvs_rsp_tensors is not None:
            self.twoshot = True
            num_core_orbitals = cvs_rsp_tensors['num_core']
            num_val_orbitals  = nocc - num_core_orbitals 

            num_intermediate_states = len(cvs_rsp_tensors['eigenvalues'])
            num_final_states = len(rsp_tensors['eigenvalues'])

            occupied_core = num_core_orbitals
            core_states = list(range(num_intermediate_states))
            val_states = list(range(num_final_states))

            self._approach_string = (f'Running RIXS in the two‑shot approach '
                             f'with {num_intermediate_states} intermediate states.')

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
            mask = (detuning >= -tol) & (detuning <= self.core_cutoff)

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
                self._approach_string = (f'Running RIXS in the restricted‑subspace approach with '
                                 f'{num_intermediate_states} intermediate states; '
                                 f'{self.num_final_states} final states kept.')

                num_final_states = self.num_final_states
                val_states = val_states[:self.num_final_states]
            else:
                self._approach_string = (f'Running RIXS in the restricted‑subspace approach with '
                                 f'{num_intermediate_states} intermediate states.')

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

        # get transformation matrices from valence to core basis
        # returns identity if cvs_scf_tensors are not given (None)
        U_occ, U_vir = None, None
        if cvs_scf_tensors is not None:
            U_occ, U_vir = self.get_transformation_mats(scf_tensors, cvs_scf_tensors, 
                                                    mo_core_indices + mo_val_indices, mo_vir_indices)
        core_mats = self.preprocess_core_eigvecs(core_eigvecs, occupied_core,
                                                 num_val_orbitals, num_vir_orbitals, U_occ, U_vir)
        
        dipole_integrals = compute_electric_dipole_integrals(molecule, basis, [0.0,0.0,0.0])

        # store state and orbital information used in computation
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
            init_photon_set = False
            # assume first core resonance
            if cvs_rsp_tensors is None:
                osc_arr = rsp_tensors['oscillator_strengths']
            else:
                osc_arr = cvs_rsp_tensors['oscillator_strengths']
            for eig, osc in zip(core_eigvals, osc_arr[core_states]):
                if osc > 1e-3:
                    self.photon_energy = [eig]
                    break
        elif isinstance(self.photon_energy, (float, int, np.floating)):
            self.photon_energy = [self.photon_energy]

        self.ene_losses             = np.zeros((num_final_states, len(self.photon_energy)))
        self.emission_enes          = np.zeros((num_final_states, len(self.photon_energy)))
        self.cross_sections         = np.zeros((num_final_states, len(self.photon_energy)))
        self.elastic_cross_sections = np.zeros((len(self.photon_energy)))
        self.scattering_amplitudes  = np.zeros((num_final_states, len(self.photon_energy),
                                           3, 3), dtype=complex)

        self.print_header()

        for w_ind, omega in enumerate(self.photon_energy):
            F_elastic = np.zeros((3,3), dtype=complex)

            for f in range(num_final_states):
                F_inelastic = np.zeros((3,3), dtype=complex)
                z_val, y_val = self.split_eigvec(valence_eigvecs[f], num_core_orbitals + num_val_orbitals,
                                            num_vir_orbitals, self.tda)
                    
                for n in range(num_intermediate_states):
                    z_core, y_core = core_mats[n]

                    gs2core, core2val = self.get_tdms(mo_occ, mo_vir, z_val, y_val, z_core, y_core)
                    
                    F_inelastic += self.scattering_amplitude_tensor(omega, core_eigvals[n], valence_eigvals[f],
                                                                    gs2core, core2val, dipole_integrals)
                    if f == 0:
                        F_elastic += self.scattering_amplitude_tensor(omega, core_eigvals[n], None,
                                                                      gs2core, None, dipole_integrals)
                
                self.emission_enes[f, w_ind] = omega - valence_eigvals[f]
                self.ene_losses[f, w_ind] = valence_eigvals[f]
                # w'/w
                prefactor_ratio = self.emission_enes[f, w_ind] / omega

                sigma = self.cross_section(F_inelastic, prefactor_ratio)

                self.cross_sections[f, w_ind] = sigma.real
                self.scattering_amplitudes[f, w_ind] = F_inelastic

            sigma_elastic = self.cross_section(F_elastic)
            self.elastic_cross_sections[w_ind] = sigma_elastic.real

            self.ostream.print_info(f'Computed RIXS cross-sections for {num_final_states} ' 
                                    f'final states at photon energy: {omega*hartree_in_ev():.2f} eV.')
            self.ostream.print_blank()
        
        results_dict = {
                    'cross_sections': self.cross_sections,
                    'elastic_cross_sections': self.elastic_cross_sections,
                    'elastic_emission': self.photon_energy,
                    'scattering_amplitudes': self.scattering_amplitudes,
                    'emission_energies': self.emission_enes,
                    'energy_losses': self.ene_losses,
                    #'excitation_energies': self.core_eigvals,
                    }

        if self.filename is not None:
            self.ostream.print_info('Writing to files...')
            self.ostream.print_blank()
            self.write_hdf5(self.filename + '_rixs')

        if not init_photon_set:
            self.photon_energy = None

        self.ostream.print_info('...done.')
        self.ostream.print_blank()

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
    def split_eigvec(eigvec, nocc, nvir, tda=False):
        """
        Gets the excitation and deexcitation matrices 
        from the eigenvectors. 
        
        :param eigvec:
            The eigenvector from a response calculation,
            either from TDA or RPA.
        :param nocc:
            Number of occupied orbitals.
        :param nvir:
            Number of virtual orbitals.
        :param tda:
            Flag for if the Tamm-Dancoff approximation is used.
        
        :return:
            The excitation and deexcitation vectors
            as matrices, in shape: (nocc, nvir).
        """
        if tda:
            z = eigvec.reshape(nocc, nvir)
            y = np.zeros_like(z)
        else:
            half = eigvec.shape[0] // 2
            z = eigvec[:half].reshape(nocc, nvir)
            y = eigvec[half:].reshape(nocc, nvir)
        return z, y
    
    @staticmethod
    def pad_matrices(z_mat, y_mat, pad_occ):
        """
        Pad matrices (e.g., from a CVS computation) with zeros
        to get them in full size.

        :param z_mat:
            The excitation matrix.
        :param y_mat:
            The deexcitation matrix.
        :pad_occ:
            The number of missing orbitals.

        :return:
            The padded, full-size matrices z and y.
        """
        padding = ((0, pad_occ), (0, 0))
        return np.pad(z_mat, padding), np.pad(y_mat, padding)
    
    def preprocess_core_eigvecs(self, core_eigvecs, occ_core,
                                num_val_orbs, num_vir_orbs, U_occ=None, U_vir=None):
        """
        Split, pad with zeros (if twoshot) and transform core eigenvectors.
        """
        pp_core_eigvecs = []

        for vec in core_eigvecs:
            z, y = self.split_eigvec(vec, occ_core,
                                    num_vir_orbs, self.tda)
            
            if self.twoshot:
                z, y = self.pad_matrices(z, y, num_val_orbs)

            if U_occ is not None and U_vir is not None:
                z = U_occ.T @ z @ U_vir
                y = U_occ.T @ y @ U_vir

            pp_core_eigvecs.append((z, y))

        return pp_core_eigvecs
    
    @staticmethod
    def get_tdms(mo_occ, mo_vir, z_val, y_val, z_core, y_core):
        """
        Get the transition density matrices, both from ground-state
        to excited state, and between two excited states (here between)
        core-excited and valence-excited states.

        :param mo_occ:
            The occupied molecular orbitals.
        :param mo_vir:
            The unoccupied/virtual molecular orbitals.
        :z_val:
            The excitation matrix (valence-excited state).
        :y_val:
            The dexcitation matrix (valence-excited state).
        :z_core:
            The excitation matrix (core-excited state).
        :y_core:
            The dexcitation matrix (core-excited state).
        """
        gs_to_core = mo_occ @ (z_core - y_core) @ mo_vir.T
        gs_to_core *= np.sqrt(2)

        core_to_val = (
            np.linalg.multi_dot([mo_vir, z_val.T, z_core, mo_vir.T]) -
            np.linalg.multi_dot([mo_occ, z_val, z_core.T, mo_occ.T]) +
            np.linalg.multi_dot([mo_occ, y_val, y_core.T, mo_occ.T]) -
            np.linalg.multi_dot([mo_vir, y_val.T, y_core, mo_vir.T])
        )
        return gs_to_core, core_to_val

    @staticmethod
    def get_transformation_mats(target_scf_tensors, initial_scf_tensors, nocc, nvir):
        """
        Gets the transformation matrices of the excitation spaces between
        two -- up-to a phase -- equal SCF wavefunctions, i.e.,
        from "initial" space to "target" space.

        :param target_scf_tensors:
            SCF tensors of the target space.
        :param initial_scf_tensors:
            SCF tensors of the initial space.
        :param nocc:
            Number of occupied orbitals.
        :param nvir:
            Number of virtual orbitals.
        
        :return:
            Transformation matrices for the occupied, U_occ, and
            for the virtual, U_vir.
        """
        if initial_scf_tensors is None:
            return np.eye(len(nocc)), np.eye(len(nvir))

        S      = target_scf_tensors['S']
        C_val  = target_scf_tensors['C_alpha']
        C_core = initial_scf_tensors['C_alpha']

        C_occ_val  = C_val[:, nocc]
        C_vir_val  = C_val[:, nvir]
        C_occ_core = C_core[:, nocc]
        C_vir_core = C_core[:, nvir]

        U_occ = C_occ_val.T @ S @ C_occ_core
        U_vir = C_vir_val.T @ S @ C_vir_core

        return U_occ, U_vir
    
    def write_hdf5(self, fname):
        """
        Writes the Pulsed response vectors to the specified output file in h5
        format. The h5 file saved contains the following datasets:

        - photon_energies
            The incoming photon energies (w, or omega), at which the RIXS amplitudes
            are computed and so also the cross-sections.
        - cross_sections
            The RIXS cross-sections per photon energy (w, or omega), per final state (f),
            with shape: (f,w).
        - emission_energies
            The outgoing, or scattered, photon energy.
        - energy_losses
            The energy (> 0) which the molecule is left with, i.e., incoming - outgoing.
        - elastic_cross_sections
            The cross-section for elastic scattering, i.e., energy loss = 0.
        - scattering_amplitudes
            The scattering ampltitude tensor, with shape: (f,w,x,y).

        :param fname:
            Name of the h5 file.
        """

        if not fname:
            raise ValueError('No filename given to write_hdf5()')

        # Add the .h5 extension if not given
        if not fname[-3:] == '.h5':
            fname += '.h5'

        # Save all the internal data to the h5 datafile named 'fname'
        try:
            with h5py.File(fname, 'w') as hf:
                hf.create_dataset('photon_energies', data=self.photon_energy)
                hf.create_dataset('cross_sections', data=self.cross_sections)
                hf.create_dataset('emission_energies', data=self.emission_enes)
                hf.create_dataset('energy_losses', data=self.ene_losses)
                hf.create_dataset('elastic_cross_sections', data=self.elastic_cross_sections)
                hf.create_dataset('scattering_amplitudes', data=self.scattering_amplitudes)

        except Exception as e:
            print('Failed to create h5 data file: {}'.format(e),
                  file=sys.stdout)

    def print_header(self):
        """
        Prints RIXS calculation setup details to output stream.
        """
        
        def _fmt_indices(lst, max_show=5):
            """
            Shring list of indices if too long.
            """
            if len(lst) > max_show:
                return f"[{lst[0]} .. {lst[-1]}]"
            return str(list(lst))

        label_width = 42
        str_width   = 90

        self.ostream.print_blank()
        title = 'Resonant Inelastic X‑ray Scattering (RIXS) Calculation'
        self.ostream.print_header(f'{title:^{str_width}}')
        self.ostream.print_header(f'{"=" * len(title):^{str_width}}')
        self.ostream.print_blank()

        gamma_ev = self.gamma * hartree_in_ev()
        basic_fields = {
            'Scattering angle (theta) (rad)'   : f'{self.theta}',
            'Lifetime broadening (gamma) (eV)'  : f'{gamma_ev:.3g}',
            'Core‑excited energy cutoff (a.u.)' : f'{self.core_cutoff}',
        }
        for label, val in basic_fields.items():
            self.ostream.print_header(f'{label:<{label_width}} : {val}'.ljust(str_width))
        self.ostream.print_blank()

        if self.photon_energy:
            au_str = ', '.join(f'{e:.4f}'        for e in self.photon_energy)
            ev_str = ', '.join(f'{e*hartree_in_ev():.2f}' for e in self.photon_energy)
            self.ostream.print_header(f'Incoming photon energies (a.u.)'.ljust(label_width) +
                                    f': {au_str}'.ljust(str_width-label_width-3))
            self.ostream.print_header(f'Incoming photon energies (eV) '.ljust(label_width) +
                                    f': {ev_str}'.ljust(str_width-label_width-3))
            self.ostream.print_blank()

        if getattr(self, "_approach_string", None):
            self.ostream.print_header(self._approach_string.ljust(str_width))
            self.ostream.print_blank()

        if self.orb_and_state_dict:
            od = self.orb_and_state_dict
            orb_fields = {
                'Core orbitals'    : _fmt_indices(od["mo_core_indices"]),
                'Valence orbitals' : _fmt_indices(od["mo_val_indices"]),
                'Virtual orbitals' : _fmt_indices(od["mo_vir_indices"]),
            }
            self.ostream.print_header('Orbital index ranges:'.center(str_width))
            for lab, val in orb_fields.items():
                self.ostream.print_header(f'{lab:<{label_width-2}} : {val}'.ljust(str_width))
            self.ostream.print_blank()

            state_fields = {
                'Intermediate/core states' : _fmt_indices(od["core_states"]),
                'Final/valence states'     : _fmt_indices(od["val_states"]),
            }
            self.ostream.print_header('State index sets:'.center(str_width))
            for lab, val in state_fields.items():
                self.ostream.print_header(f'{lab:<{label_width-2}} : {val}'.ljust(str_width))
            self.ostream.print_blank()

            self.ostream.print_header(
                f'Number of intermediate states : {od["num_intermediate_states"]}'.ljust(str_width))
            self.ostream.print_header(
                f'Number of final states        : {od["num_final_states"]}'.ljust(str_width))
            self.ostream.print_blank()

