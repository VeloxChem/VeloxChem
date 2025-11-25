#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from mpi4py import MPI
import numpy as np
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
    Implements the RIXS driver in a linear-response framework with two 
    approaches: the two-shot and restricted-subspace approximation.

    # vlxtag: RHF, RIXS
    # vlxtag: RKS, RIXS

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - photon_energy: Incoming photon energy; omega (a.u.)
        - theta: Angle between incident polarization-
          and outgoing propagation vectors (rad.)
        - gamma: Life-time broadening (FWHM) (a.u.)
        - final_state_cutoff: Energy window of final states to include (a.u.)'),
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
        self.comm  = comm
        self.rank  = comm.Get_rank()
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

        self.rixs_dict = {}
        self.filename = None

        # TODO?: should this be in terms of energy or number of states?
        #self.num_final_states = None
        self.final_state_cutoff = None

        # input keywords
        self.input_keywords = {
            'rixs': {
                'theta': ('float', 'angle between incident polarization vector and propagation vector of outgoing'),
                'gamma': ('float', 'broadening term (FWHM)'),
                'photon_energy': ('list', 'list of incoming photon energies'),
                'final_state_cutoff': ('float', 'energy window of final states to include'),
            },
        }

    def update_settings(self, rixs_dict=None, method_dict=None):
        """
        Updates settings in RixsDriver.

        :param rixs_dict:
            The dictionary of rixs input.
        """

        if method_dict is None:
            method_dict = {}
        if rixs_dict is None:
            rixs_dict = {}

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

        nocc = molecule.number_of_alpha_electrons()

        self.twoshot = (cvs_rsp_tensors is not None)
        init_photon_set = True

        if self.rank == mpi_master():
            num_vir_orbitals = rsp_tensors['num_vir']
        else:
            num_vir_orbitals = None
        num_vir_orbitals = self.comm.bcast(num_vir_orbitals, root=mpi_master())

        if self.twoshot:
            self._approach_string = (f'Running RIXS calculation in the two‑shot approach')

            num_core_orbitals       = cvs_rsp_tensors['num_core']
            num_valence_orbitals    = nocc - num_core_orbitals 
            num_intermediate_states = len(cvs_rsp_tensors['eigenvalues'])
            num_final_states        = len(rsp_tensors['eigenvalues'])
            occupied_core           = num_core_orbitals
            core_states             = list(range(num_intermediate_states))
            val_states              = list(range(num_final_states))

        else:
            self._approach_string = (f'Running RIXS calculation in the restricted‑subspace approach')

            num_valence_orbitals = rsp_tensors['num_val']
            num_core_orbitals    = rsp_tensors['num_core']
            assert_msg_critical(num_core_orbitals > 0,
                                 'No core orbitals indicated in the response tensor.')

            # identify the energy of the lowest core-excited state
            first_core_ene = self._first_core_energy(rsp_tensors)
            detuning = rsp_tensors['eigenvalues'] - first_core_ene

            # identify (and possibly remove unphysical valence-excited states) the core-excited states
            core_states = self._core_state_indices(rsp_tensors, detuning)
            num_intermediate_states = len(core_states)
            assert_msg_critical(num_intermediate_states > 0,
                                'Too few excited states included in response calculation.')
            # identify the valence-excited states
            val_states = self._valence_state_indices(detuning, rsp_tensors['eigenvalues'])
            num_final_states = len(val_states)

            # for compatibiltiy with the two-shot approach
            cvs_rsp_tensors = rsp_tensors
            # for compatibiltiy with the two-shot approach
            occupied_core   = num_core_orbitals + num_valence_orbitals

        mo_core_indices = list(range(num_core_orbitals))
        mo_val_indices  = list(range(nocc - num_valence_orbitals, nocc))
        mo_vir_indices  = list(range(nocc, nocc + num_vir_orbitals))
        mo_occ          = scf_tensors['C_alpha'][:, mo_core_indices + mo_val_indices]
        mo_vir          = scf_tensors['C_alpha'][:, mo_vir_indices]

        core_eigvals    = cvs_rsp_tensors['eigenvalues'][core_states]
        valence_eigvals = rsp_tensors['eigenvalues'][val_states]
        core_eigvecs    = self._get_eigvecs(cvs_rsp_tensors, core_states, "core")
        valence_eigvecs = self._get_eigvecs(rsp_tensors, val_states, "valence")

        # get transformation matrices from valence to core basis, if the
        # valence- and core-excited states were obtained from independent SCFs.
        # returns identity if cvs_scf_tensors are not given (None).
        U_occ, U_vir = None, None
        if cvs_scf_tensors is not None:
            U_occ, U_vir = self.get_transformation_mats(scf_tensors, cvs_scf_tensors, 
                                                    mo_core_indices + mo_val_indices, mo_vir_indices)
        # TODO parallelise, and broadcast?
        core_mats = self._preprocess_core_eigvecs(core_eigvecs, occupied_core,
                                                 num_valence_orbitals, num_vir_orbitals, U_occ, U_vir)
        
        # TODO parallelise, and broadcast?
        if self.rank == mpi_master():
            dipole_integrals = compute_electric_dipole_integrals(molecule, basis, [0.0,0.0,0.0])
        else:
            dipole_integrals = None
        dipole_integrals = self.comm.bcast(dipole_integrals, root=mpi_master())

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
        
        # if incoming photon energy is not set, set it to match the
        # first core-excited state with osc_strength > 1e-3
        if self.photon_energy is None:
            init_photon_set = False
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

        # Define the return objects
        self.ene_losses             = np.zeros((num_final_states, len(self.photon_energy)))
        self.emission_enes          = np.zeros((num_final_states, len(self.photon_energy)))
        self.cross_sections         = np.zeros((num_final_states, len(self.photon_energy)))
        self.elastic_cross_sections = np.zeros((len(self.photon_energy)))
        self.scattering_amplitudes  = np.zeros((num_final_states, len(self.photon_energy),
                                           3, 3), dtype=complex)

        if self.rank == mpi_master():
            self._print_header(molecule)

        ave, rem = divmod(num_final_states, self.nodes)
        counts = [ave + 1 if p < rem else ave for p in range(self.nodes)]
        f_start = sum(counts[:self.rank])
        f_end   = sum(counts[:self.rank + 1])
        n_local = f_end - f_start

        for w_ind, omega in enumerate(self.photon_energy):
            F_elastic_local = np.zeros((3, 3), dtype=complex)
            local_cross_sections    = np.empty(n_local)
            local_amplitude_tensors = np.empty((n_local, 3, 3), dtype=complex)
            local_emission_energies = np.empty(n_local)
            local_energy_losses     = np.empty(n_local)

            for local_i, f in enumerate(range(f_start, f_end)):
                F_inelastic = np.zeros((3, 3), dtype=complex)
                z_val, y_val = self.split_eigvec(valence_eigvecs[f],
                                                num_core_orbitals + num_valence_orbitals,
                                                num_vir_orbitals, self.tda)
                for n in range(num_intermediate_states):
                    z_core, y_core = core_mats[n]
                    gs2core, core2val = self.get_tdms(mo_occ, mo_vir, z_val, y_val, z_core, y_core)

                    F_inelastic += self.scattering_amplitude_tensor(
                        omega, core_eigvals[n], valence_eigvals[f], gs2core, core2val, dipole_integrals
                    )
                    if f == 0:
                        F_elastic_local += self.scattering_amplitude_tensor(
                            omega, core_eigvals[n], None, gs2core, None, dipole_integrals
                        )

                emission_energy = omega - valence_eigvals[f]
                energy_loss     = valence_eigvals[f]
                prefactor       = emission_energy / omega

                # Fill local result slices
                local_emission_energies[local_i] = emission_energy
                local_energy_losses[local_i]     = energy_loss
                local_cross_sections[local_i]    = self.cross_section(F_inelastic, prefactor).real
                local_amplitude_tensors[local_i] = F_inelastic

            # Reduce per omega
            F_elastic = self.comm.allreduce(F_elastic_local, op=MPI.SUM)
            parts = self.comm.allgather(
                (f_start, f_end,
                local_emission_energies,
                local_energy_losses,
                local_cross_sections,
                local_amplitude_tensors)
            )

            if self.rank == mpi_master():
                self.elastic_cross_sections[w_ind] = self.cross_section(F_elastic).real
                for start, end, em_en_part, loss_part, cs_part, amp_part in parts:
                    self.emission_enes[start:end, w_ind]         = em_en_part
                    self.ene_losses[start:end, w_ind]            = loss_part
                    self.cross_sections[start:end, w_ind]        = cs_part
                    self.scattering_amplitudes[start:end, w_ind] = amp_part

                self.ostream.print_info(
                    f'Computed RIXS cross-sections for {num_final_states} final states '
                    f'at photon energy: {omega*hartree_in_ev():.2f} eV.'
                )
                self.ostream.print_blank()
                self.ostream.flush()

        results_dict = {
            'cross_sections': self.cross_sections,
            'elastic_cross_sections': self.elastic_cross_sections,
            'elastic_emission': self.photon_energy,
            'scattering_amplitudes': self.scattering_amplitudes,
            'emission_energies': self.emission_enes,
            'energy_losses': self.ene_losses,
        }

        if self.filename is not None and self.rank == mpi_master():
            self.ostream.print_info('Writing to files...')
            self.ostream.print_blank()
            self.ostream.flush()
            self._write_hdf5(self.filename + '_rixs')

        if not init_photon_set:
            self.photon_energy = None

        if self.rank == mpi_master():
            self.ostream.print_info('...done.')
            self.ostream.print_blank()
            self.ostream.flush()

        return results_dict
    
    def _first_core_energy(self, rsp_tensors):
        """
        Return energy of the first core-excited state (by label order).

        :param rsp_tensors:
            Dictionary containing excitation information.

        :return:
            The energy (float) corresponding to the first
            state labeled as 'core' in exciitation_details.
        """
        for k, entries in enumerate(rsp_tensors['excitation_details']):
            if entries[0].split()[0].startswith("core"):
                return rsp_tensors['eigenvalues'][k]
        raise ValueError("No core-labeled state found in excitation_details.")

    def _core_state_indices(self, rsp_tensors, detuning, tol=1e-6):
        """
        Filter out (possibly) unphysical non-core states.
        Keep states with detuning >= -tol and labeled 'core'.

        :param rsp_tensors:
            Dictionary containing the response data.
        :param detuning:
            Array of detuning values relative to the
            first core-excited state.
        :param tol:
            Numerical tolerance used for detuning cutoff.

        :return: 
            List of core-excited states and array of detuning.
        """

        # -tol to make sure the first_core_ene-state is included
        init_core_states = np.where(detuning >= -tol)[0]

        core_states = []
        for state in init_core_states:
            entry = rsp_tensors['excitation_details'][state][0].split()
            label = entry[0]
            if label.startswith("core"):
                core_states.append(int(state))

        return core_states

    def _valence_state_indices(self, detuning, eigenvalues, tol=1e-6):
        """
        Final/valence states: detuning < 0, 
        possibly windowed by self.final_state_cutoff.
        
        :param detuning:
            Array of detuning values relative
            to the first core-excited state.
        :param tol:
            Numerical tolerance used for detuning cutoff.
        
        :return: 
            Array containing the indices of 
            final (valence) states.
        """

        val_states = np.where(detuning < 0)[0]
        if self.final_state_cutoff is not None:
            cutoff = float(self.final_state_cutoff) + tol
            keep = eigenvalues[val_states] <= cutoff
            val_states = val_states[keep]
            assert_msg_critical(val_states.size > 0, "No valence states pass the cutoff.")
        return val_states.astype(int)

    def scattering_amplitude_tensor(self, omega, core_eigenvalue, val_eigenvalue,
                                    intermediate_tdens, final_tdens, dipole_integrals):
        """
        The RIXS scattering amplitude.
        
        :param omega:
            Energy of incoming photon.
        :param f: 
            The final (valence) state index.
        :param eigenvalue: 
            Energy of intermediate (core) excited state.
        :param intermediate_tdens:
            Transition density from ground to
            intermediate (core) excited state (AO).
        :param final_tdens:
            Transition density matrix from 
            intermediate to final excited state (AO).
        :param dipole_integrals:
            Electric dipole integrals (length gauge) in AO-basis.

        :return: 
            The scattering amplitude tensor: shape = (3,3)

        """
        gamma_hwhm = self.gamma / 2.0
        e_n = 1.0 / (omega - (core_eigenvalue + 1j * gamma_hwhm))

        if val_eigenvalue is None or intermediate_tdens is None:
            # Elastic line, note that the term prop. to A^2 is not included
            # Phys. Rep. 324, 1–105 (2000), DOI: 10.1016/S0370-1573(99)00003-4
            core_eigvals2 = core_eigenvalue**2
            gs_exc_tdipole = np.array([
                                np.sum(intermediate_tdens * dipole_integrals[i])
                                for i in range(3)])
            outer_product = np.matmul(gs_exc_tdipole[:, np.newaxis], gs_exc_tdipole[np.newaxis, :])
            scatt_amp = e_n * core_eigvals2 * outer_product

        else:
            # Inelastic line
            omega_product = (core_eigenvalue - val_eigenvalue) * core_eigenvalue
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
        
        Journal of Chemical Theory and Computation 2021 17 (5), 3031-3038
        DOI: 10.1021/acs.jctc.1c00144

        Journal of Chemical Theory and Computation 2017 13 (11), 5552-5559
        DOI: 10.1021/acs.jctc.7b00636
        
        :param F:
            The scattering amplitude tensor.
        :param omegaprime_omega:
            The ratio between outgoing- and incoming frequencies, w'/w.
            
        :return:
            The scattering cross-section.
        """
        
        # F_xy (F^(xy))^*
        F2 = np.vdot(F, F).real
        # F_xy (F^(yx))^*
        FF_T_conj = np.vdot(F, F.T).real
        # F_xx (F^(yy))^*
        tr_F2 = np.abs(np.trace(F))**2

        sigma_f = omegaprime_omega * (1.0/15.0) * (
                    (2.0 - 0.5*np.sin(self.theta)**2) * F2 +
                    (0.75*np.sin(self.theta)**2 - 0.5) * (FF_T_conj + tr_F2))
        
        return sigma_f
    
    def _get_eigvecs(self, tensor, states, label):
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
    
    def _preprocess_core_eigvecs(self, core_eigvecs, occ_core,
                                num_val_orbs, num_vir_orbs,
                                U_occ=None, U_vir=None):
        """
        Split, pad with zeros (if twoshot) and 
        possibly transform core eigenvectors if needed.
        """
        pp_core_eigvecs = []

        for vec in core_eigvecs:
            z, y = self.split_eigvec(vec, occ_core,
                                    num_vir_orbs, self.tda)
            
            if self.twoshot:
                z, y = self.pad_matrices(z, y, num_val_orbs)

            if U_occ is not None and U_vir is not None:
                z = np.linalg.multi_dot([U_occ.T, z, U_vir])
                y = np.linalg.multi_dot([U_occ.T, y, U_vir])

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

        gs_to_core = np.linalg.multi_dot([mo_occ, z_core - y_core, mo_vir.T])
        gs_to_core *= np.sqrt(2.0)

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

        C_occ_val  = C_val[:, nocc].copy()
        C_vir_val  = C_val[:, nvir].copy()
        C_occ_core = C_core[:, nocc].copy()
        C_vir_core = C_core[:, nvir].copy()

        U_occ = np.linalg.multi_dot([C_occ_val.T, S, C_occ_core])
        U_vir = np.linalg.multi_dot([C_vir_val.T, S, C_vir_core])

        return U_occ, U_vir
    
    def _write_hdf5(self, fname):
        """
        Writes the RIXS results to the specified output file in h5
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
            raise ValueError('No filename given to _write_hdf5()')

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

    def _print_header(self, molecule):
        """
        Prints RIXS calculation setup details to output stream.
        """
        
        def _format_indices(lst, max_show=3):
            """
            Shrink list of indices if too long.
            """
            arr = np.asarray(lst)
            if arr.size == 0:
                return "[]"
            
            items = [str(x) for x in arr.tolist()]
            if len(items) > max_show:
                return f"[{items[0]} ... {items[-1]}]"
            return "[" + ", ".join(items) + "]"

        
        def _join_indices(parts, n):
            sep = ", " if n <= 3 else " "
            return sep.join(parts)
        
        def _format_orbital_names(indices, orb_type, nocc):
            """
            Convert orbital indices to orbital labels.

            :param indices:
                List or array of orbital indices.
            :param orb_type:
                One of 'core', 'valence', or 'virtual'.
            :param nocc:
                Number of occupied orbitals.

            :return:
                List of orbital labels.
            """
            labels = []
            if orb_type == "core":
                # Core orbitals: Core1, Core2, ...
                labels = [f"Core{idx + 1}" for idx in indices]

            elif orb_type == "valence":
                # Valence orbitals: HOMO-1, HOMO, ...
                max_val = max(indices)
                for idx in indices:
                    diff = max_val - idx
                    labels.append("HOMO" if diff == 0 else f"HOMO-{diff}")

            elif orb_type == "virtual":
                # Virtual orbitals: LUMO, LUMO+1, ...
                for idx in indices:
                    diff = idx - nocc
                    labels.append("LUMO" if diff == 0 else f"LUMO+{diff}")
            
            return labels

        str_width   = 60
        nocc        = molecule.number_of_alpha_electrons()

        self.ostream.print_blank()
        title = 'Resonant Inelastic X‑ray Scattering (RIXS) Setup'
        self.ostream.print_header(f'{title:^{str_width}}')
        self.ostream.print_header(f'{"=" * len(title):^{str_width}}')
        self.ostream.print_blank()

        cur_str = 'Scattering angle (theta) [rad]   : {:.2f}'.format(self.theta)
        self.ostream.print_header(cur_str.ljust(str_width))
        gamma_ev = self.gamma * hartree_in_ev()
        cur_str = 'Lifetime broadening (gamma) [eV] : {:.2f}'.format(gamma_ev)
        self.ostream.print_header(cur_str.ljust(str_width))

        if self.photon_energy:
            if len(self.photon_energy) > 3:
                display_energies = self.photon_energy[:1] + ["..."] + self.photon_energy[-1:]
            else:
                display_energies = self.photon_energy

            au_parts = [f"{e:.4f}" if isinstance(e, (int, float)) else e for e in display_energies]
            ev_parts = [f"{e * hartree_in_ev():.2f}" if isinstance(e, (int, float)) else e for e in display_energies]
            au_str = _join_indices(au_parts, len(self.photon_energy))
            ev_str = _join_indices(ev_parts, len(self.photon_energy))

            self.ostream.print_header(
                f"Incoming photon energies [a.u.]  : {au_str}".ljust(str_width))
            self.ostream.print_header(
                f"Incoming photon energies [eV]    : {ev_str}".ljust(str_width))

        if self.orb_and_state_dict:
            od = self.orb_and_state_dict

            self.ostream.print_header(
                f'Number of intermediate states    : {od["num_intermediate_states"]}'.ljust(str_width))
            self.ostream.print_header(
                f'Number of final states           : {od["num_final_states"]}'.ljust(str_width))
            
            if self.final_state_cutoff is not None:
                self.ostream.print_header(
                f'Final state energy cutoff [eV]   : {hartree_in_ev() * self.final_state_cutoff:.2f}'.ljust(str_width))

            self.ostream.print_blank()

            self.ostream.print_header('State index sets'.center(str_width))
            state_list = _format_indices(od["core_states"])
            self.ostream.print_header(
                f'Intermediate/core states         : {state_list}'.ljust(str_width)
                )

            state_list = _format_indices(od["val_states"])
            self.ostream.print_header(
                f'Final/valence states             : {state_list}'.ljust(str_width)
                )
            self.ostream.print_blank()

            self.ostream.print_header('Orbital sets'.center(str_width))
            orbital_labels = _format_orbital_names(od["mo_core_indices"], "core", nocc)
            orbital_str    = _format_indices(orbital_labels)
            self.ostream.print_header(
                f'Core orbitals                    : {orbital_str}'.ljust(str_width)
                )

            orbital_labels = _format_orbital_names(od["mo_val_indices"], "valence", nocc)
            orbital_str    = _format_indices(orbital_labels)
            self.ostream.print_header(
                f'Valence orbitals                 : {orbital_str}'.ljust(str_width)
                )
            
            orbital_labels = _format_orbital_names(od["mo_vir_indices"], "virtual", nocc)
            orbital_str    = _format_indices(orbital_labels)
            self.ostream.print_header(
                f'Virtual orbitals                 : {orbital_str}'.ljust(str_width)
                )
            self.ostream.print_blank()

        if getattr(self, "_approach_string", None):
            self.ostream.print_info(self._approach_string)
            self.ostream.print_blank()
        
        self.ostream.flush()

    @staticmethod
    def lorentzian_broadening(x, y, xmin, xmax, xstep, br):
        xi = np.arange(xmin, xmax, xstep)
        yi = np.zeros(len(xi))
        for i in range(len(xi)):
            for k in range(len(x)):
                yi[i] = (yi[i] + y[k] * br / ((xi[i] - x[k])**2 +
                                              (br / 2.0)**2) / np.pi)
        return xi, yi

    @staticmethod
    def gaussian_broadening(x, y, xmin, xmax, xstep, br):
        br_g = br / np.sqrt(4.0 * 2.0 * np.log(2))
        xi = np.arange(xmin, xmax, xstep)
        yi = np.zeros(len(xi))
        for i in range(len(xi)):
            for k in range(len(y)):
                yi[i] = yi[i] + y[k] * np.exp(-((xi[i] - x[k])**2) /
                                              (2 * br_g**2))
        return xi, yi
    
