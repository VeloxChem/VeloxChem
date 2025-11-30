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
from pathlib import Path
import numpy as np
import h5py
import sys

from .veloxchemlib import hartree_in_ev, mpi_master
from .oneeints import compute_electric_dipole_integrals
from .outputstream import OutputStream
from .linearsolver import LinearSolver
from .lreigensolver import LinearResponseEigenSolver
from .tdaeigensolver import TdaEigenSolver
from .errorhandler import assert_msg_critical
from .inputparser import parse_input


class RixsDriver(LinearSolver):
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

        super().__init__(comm, ostream)

        # for RIXS, the default convergence threshold is tighter
        self.conv_thresh = 2.0e-5

        self.tamm_dancoff = False

        # method settings
        self.photon_energy = None
        self.theta = 0
        # a.u., FWHM
        # 1000.0 / hartree_in_wavenumber()
        self.gamma = 0.124 / hartree_in_ev()

        self.nstates = 5
        self.num_core_states = 5

        # restricted subspace
        self.restricted_subspace = True

        self.num_core_orbitals = 0
        self.num_virtual_orbitals = 0
        self.num_valence_orbitals = 0

        self._orb_and_state_dict = None

        # TODO?: should this be in terms of energy or number of states?
        #self.num_final_states = None
        self.final_state_cutoff = None

        # input keywords
        self._input_keywords['response'].update({
            'theta':
                ('float', 'angle between incident polarization vector and ' +
                    'outgoing propagation vector'),
            'gamma': ('float', 'broadening term (FWHM)'),
            'photon_energy': ('list', 'list of incoming photon energies'),
            'final_state_cutoff':
                ('float', 'energy window of final states to include'),
            'nstates': ('int', 'number of excited states'),
            'num_core_states': ('int', 'number of core-excited states'),
            'num_core_orbitals': ('int', 'number of involved core-orbitals'),
            'restricted_subspace':
                ('bool', 'restricted subspace approximation'),
            'num_valence_orbitals':
                ('int', 'number of involved valence orbitals'),
            'num_virtual_orbitals':
                ('int', 'number of involved virtual orbitals'),
        })

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates settings in RixsDriver.

        :param rsp_dict:
            The dictionary of response input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}
        
        key = 'photon_energy'
        if key in rsp_dict:
            val = rsp_dict[key]

            if isinstance(val, list):
                rsp_dict[key] = [float(str(x).strip()) for x in val]

            elif isinstance(val, str):
                parts = val.replace(',', ' ').split()
                rsp_dict[key] = [float(x) for x in parts]

        super().update_settings(rsp_dict, method_dict)

    def print_info(self):
        """
        Prints relevant orbital and state information
        """
        if self._orb_and_state_dict is None:
            self.ostream.print_warning('Compute first!')
        else:
            info = self._orb_and_state_dict
            intermed_states = info['num_intermediate_states']
            final_states = info['num_final_states']
            mo_c_ind = info['mo_core_indices']
            mo_val_ind = info['mo_valence_indices']
            mo_vir_ind = info['mo_virtual_indices']
            ce_states = info['core_states']
            ve_states = info['valence_states']
            print(f'\nNumber of states: (intermediate, final): ({intermed_states}, {final_states})')
            print(f'\nMO indices (core, valence, virtual): ({mo_c_ind}, {mo_val_ind}, {mo_vir_ind})')
            print(f'\nState indices (intermediate, final): ({ce_states}, {ve_states})')

    def compute(self, molecule, basis, scf_results,
                rsp_results=None, cvs_rsp_results=None):
        """
        Computes RIXS properties.

        :param molecule:
            The molecule object.
        :param basis:
            The AO basis set.
        :scf_results:
            The results dictionary from converged SCF wavefunction.
        :rsp_results:
            The linear-response results dictionary.
        :cvs_rsp_results:
            The core-valence-separated (CVS) linear-response 
            results dictionary.

        NOTE: To run full diagonalization, do a restricted
              subspace calculation with the full space, 
              indicating which orbitals define the core.

        :return:
            Dictionary with cross-sections, outgoing photon energy
            in (energy loss and full energy (emission)),
            and the scattering amplitude tensor.
        """

        if rsp_results is None:

            rsp_keys = [
                'restricted_subspace',
                'nstates',
                'num_virtual_orbitals',
                'num_valence_orbitals',
                'num_core_orbitals',
                'conv_thresh',
                'restart',
            ]

            if self.tamm_dancoff:
                rsp_drv = TdaEigenSolver(self.comm, self.ostream)
            else:
                rsp_drv = LinearResponseEigenSolver(self.comm, self.ostream)

            for key in rsp_keys:
                if hasattr(self, key):
                    setattr(rsp_drv, key, getattr(self, key))

            if rsp_drv.restricted_subspace:
                assert_msg_critical(
                    (rsp_drv.num_core_orbitals > 0 and
                     rsp_drv.num_valence_orbitals > 0 and
                     rsp_drv.num_virtual_orbitals > 0),
                    'Invalid number of core/valence/virtual orbitals for ' +
                    'restricted_subspace')

            if self.checkpoint_file is not None:
                fpath = Path(self.checkpoint_file)
                fpath = fpath.with_name(fpath.stem)
                rsp_drv.checkpoint_file = str(fpath) + '_rixs.h5'

            rsp_results = rsp_drv.compute(molecule, basis, scf_results)

        if cvs_rsp_results is None and (not self.restricted_subspace):

            cvs_rsp_keys = [
                'num_core_orbitals',
                'conv_thresh',
                'restart',
            ]

            if self.tamm_dancoff:
                cvs_rsp_drv = TdaEigenSolver(self.comm, self.ostream)
            else:
                cvs_rsp_drv = LinearResponseEigenSolver(self.comm, self.ostream)

            for key in cvs_rsp_keys:
                if hasattr(self, key):
                    setattr(cvs_rsp_drv, key, getattr(self, key))

            cvs_rsp_drv.core_excitation = True
            cvs_rsp_drv.nstates = self.num_core_states

            assert_msg_critical(cvs_rsp_drv.num_core_orbitals > 0,
                                'Invalid number of core orbitals')

            if self.checkpoint_file is not None:
                fpath = Path(self.checkpoint_file)
                fpath = fpath.with_name(fpath.stem)
                cvs_rsp_drv.checkpoint_file = str(fpath) + '_rixs_cvs.h5'

            cvs_rsp_results = cvs_rsp_drv.compute(molecule, basis, scf_results)

        nocc = molecule.number_of_alpha_electrons()

        self.twoshot = (not self.restricted_subspace)
        init_photon_set = True

        if self.rank == mpi_master():
            num_vir_orbitals = rsp_results['num_virtual']
        else:
            num_vir_orbitals = None
        num_vir_orbitals = self.comm.bcast(num_vir_orbitals, root=mpi_master())

        if self.twoshot:
            self._approach_string = 'Running RIXS calculation in the two‑shot approach'

            if self.rank == mpi_master():
                num_core_orbitals = cvs_rsp_results['num_core']
                num_valence_orbitals = nocc - num_core_orbitals
                num_intermediate_states = len(cvs_rsp_results['eigenvalues'])
                num_final_states = len(rsp_results['eigenvalues'])
            else:
                num_core_orbitals = None
                num_valence_orbitals = None
                num_intermediate_states = None
                num_final_states = None

            num_valence_orbitals, num_core_orbitals = self.comm.bcast(
                (num_valence_orbitals, num_core_orbitals), root=mpi_master())

            num_intermediate_states, num_final_states = self.comm.bcast(
                (num_intermediate_states, num_final_states), root=mpi_master())

            occupied_core           = num_core_orbitals
            core_states             = list(range(num_intermediate_states))
            val_states              = list(range(num_final_states))

        else:
            self._approach_string = 'Running RIXS calculation in the restricted‑subspace approach'

            if self.rank == mpi_master():
                num_valence_orbitals = rsp_results['num_valence']
                num_core_orbitals    = rsp_results['num_core']
                assert_msg_critical(num_core_orbitals > 0,
                                    'No core orbitals indicated in the response results.')

                # identify the energy of the lowest core-excited state
                first_core_ene = self._first_core_energy(rsp_results)
                detuning = rsp_results['eigenvalues'] - first_core_ene

                # identify (and possibly remove unphysical valence-excited states) the core-excited states
                core_states = self._core_state_indices(rsp_results, detuning)
                num_intermediate_states = len(core_states)
                assert_msg_critical(num_intermediate_states > 0,
                                    'Too few excited states included in response calculation.')
                # identify the valence-excited states
                val_states = self._valence_state_indices(detuning, rsp_results['eigenvalues'])
                num_final_states = len(val_states)
            else:
                num_valence_orbitals = None
                num_core_orbitals    = None
                num_intermediate_states = None
                num_final_states = None
                core_states = None
                val_states = None

            num_valence_orbitals, num_core_orbitals = self.comm.bcast(
                (num_valence_orbitals, num_core_orbitals), root=mpi_master())

            num_intermediate_states, num_final_states = self.comm.bcast(
                (num_intermediate_states, num_final_states), root=mpi_master())

            core_states = self.comm.bcast(core_states, root=mpi_master())
            val_states = self.comm.bcast(val_states, root=mpi_master())

            # for compatibility with the two-shot approach
            cvs_rsp_results = rsp_results
            occupied_core   = num_core_orbitals + num_valence_orbitals

        mo_core_indices = list(range(num_core_orbitals))
        mo_val_indices  = list(range(nocc - num_valence_orbitals, nocc))
        mo_vir_indices  = list(range(nocc, nocc + num_vir_orbitals))

        if self.rank == mpi_master():
            mo_occ = scf_results['C_alpha'][:, mo_core_indices + mo_val_indices].copy()
            mo_vir = scf_results['C_alpha'][:, mo_vir_indices].copy()

            core_eigvals = cvs_rsp_results['eigenvalues'][core_states].copy()
            valence_eigvals = rsp_results['eigenvalues'][val_states].copy()
        else:
            mo_occ = None
            mo_vir = None

            core_eigvals = None
            valence_eigvals = None

        mo_occ = self.comm.bcast(mo_occ, root=mpi_master())
        mo_vir = self.comm.bcast(mo_vir, root=mpi_master())

        core_eigvals = self.comm.bcast(core_eigvals, root=mpi_master())
        valence_eigvals = self.comm.bcast(valence_eigvals, root=mpi_master())

        core_eigvecs    = self._get_eigvecs(cvs_rsp_results, core_states, "core")
        valence_eigvecs = self._get_eigvecs(rsp_results, val_states, "valence")

        core_mats = self._preprocess_core_eigvecs(core_eigvecs, occupied_core,
                                                  num_valence_orbitals, num_vir_orbitals)
        
        if self.rank == mpi_master():
            dipole_integrals = compute_electric_dipole_integrals(molecule, basis, [0.0,0.0,0.0])
        else:
            dipole_integrals = None
        dipole_integrals = self.comm.bcast(dipole_integrals, root=mpi_master())

        # store state and orbital information used in computation
        self._orb_and_state_dict = {
            'num_intermediate_states': num_intermediate_states,
            'num_final_states': num_final_states,
            'mo_core_indices': mo_core_indices,
            'mo_valence_indices': mo_val_indices,
            'mo_virtual_indices': mo_vir_indices,
            'core_states': core_states,
            'valence_states': val_states,
        }
        
        # if incoming photon energy is not set, set it to match the
        # first core-excited state with osc_strength > 1e-3
        if self.photon_energy is None:
            init_photon_set = False

            if self.rank == mpi_master():
                if cvs_rsp_results is None:
                    osc_arr = rsp_results['oscillator_strengths']
                else:
                    osc_arr = cvs_rsp_results['oscillator_strengths']
            else:
                osc_arr = None
            osc_arr = self.comm.bcast(osc_arr, root=mpi_master())

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
                                                 num_vir_orbitals, self.tamm_dancoff)
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
            parts = self.comm.allgather((
                f_start,
                f_end,
                local_emission_energies,
                local_energy_losses,
                local_cross_sections,
                local_amplitude_tensors,
            ))

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
                self._print_rixs_data(f'RIXS cross-sections at incident X-ray energy '
                                      f'{omega*hartree_in_ev():.2f} eV, energy-loss mode',
                                      self.ene_losses[:,w_ind], self.cross_sections[:,w_ind], self.elastic_cross_sections[w_ind])
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
            # TODO: append to final h5 (instead of writting new file)
            self._write_hdf5(self.filename + '_rixs')

        if not init_photon_set:
            self.photon_energy = None

        if self.rank == mpi_master():
            self.ostream.print_info('...done.')
            self.ostream.print_blank()
            self.ostream.flush()

        return results_dict
    
    def _first_core_energy(self, rsp_results):
        """
        Return energy of the first core-excited state (by label order).

        :param rsp_results:
            Dictionary containing excitation information.

        :return:
            The energy (float) corresponding to the first
            state labeled as 'core' in exciitation_details.
        """

        # TODO: implement more robust way of determining core-excited states
        # (other than relying on the excitation_details text)

        for k, entries in enumerate(rsp_results['excitation_details']):
            if entries[0].split()[0].startswith("core"):
                return rsp_results['eigenvalues'][k]
        raise ValueError("No core-labeled state found in excitation_details.")

    def _core_state_indices(self, rsp_results, detuning, tol=1e-6):
        """
        Filter out (possibly) unphysical non-core states.
        Keep states with detuning >= -tol and labeled 'core'.

        :param rsp_results:
            Dictionary containing the response data.
        :param detuning:
            Array of detuning values relative to the
            first core-excited state.
        :param tol:
            Numerical tolerance used for detuning cutoff.

        :return: 
            List of core-excited states and array of detuning.
        """

        # TODO: implement more robust way of determining core-excited states
        # (other than relying on the excitation_details text)

        # -tol to make sure the first_core_ene-state is included
        init_core_states = np.where(detuning >= -tol)[0]

        core_states = []
        for state in init_core_states:
            entry = rsp_results['excitation_details'][state][0].split()
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
            gs_exc_tdipole = np.array([np.sum(intermediate_tdens * dipole_integrals[i])
                                       for i in range(3)])
            outer_product = np.matmul(gs_exc_tdipole[:, np.newaxis], gs_exc_tdipole[np.newaxis, :])
            scatt_amp = e_n * core_eigvals2 * outer_product

        else:
            # Inelastic line
            omega_product = (core_eigenvalue - val_eigenvalue) * core_eigenvalue
            gs_exc_tdipole = np.array([np.sum(intermediate_tdens * dipole_integrals[i])
                                       for i in range(3)])
            exc_exc_tdipole = np.array([np.sum(final_tdens * dipole_integrals[i])
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
    
    def _get_eigvecs(self, rsp_results, states, label):
        """
        Gets eigenvectors.

        :param rsp_results:
            The dictionay containint response results.
        :param states:
            The excited states.
        :param label:
            The label.

        :return:
            The eigenvectors as a list of 1D numpy arrays.
        """

        if not self.tamm_dancoff:
            # RPA
            assert_msg_critical(
                'eigenvectors_distributed' in rsp_results,
                'RixsDriver._get_eigvecs: Incorrect key in rsp_results')

            vecs = []
            for i in states:
                full_v = self.get_full_solution_vector(
                    rsp_results['eigenvectors_distributed'][i])
                full_v = self.comm.bcast(full_v, root=mpi_master())
                vecs.append(full_v)
            return np.array(vecs)

        else:
            # TDA
            assert_msg_critical(
                'eigenvectors' in rsp_results,
                'RixsDriver._get_eigvecs: Incorrect key in rsp_results')

            vecs = []
            for i in states:
                if self.rank == mpi_master():
                    full_v = rsp_results['eigenvectors'][:, i].copy()
                else:
                    full_v = None
                full_v = self.comm.bcast(full_v, root=mpi_master())
                vecs.append(full_v)
            return np.array(vecs)

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
                                 num_val_orbs, num_vir_orbs):
        """
        Split, pad with zeros (if twoshot) and 
        possibly transform core eigenvectors if needed.
        """
        pp_core_eigvecs = []

        for vec in core_eigvecs:
            z, y = self.split_eigvec(vec, occ_core,
                                     num_vir_orbs, self.tamm_dancoff)
            
            if self.twoshot:
                z, y = self.pad_matrices(z, y, num_val_orbs)

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

        # TODO: use pathlib for handling file extension
        # Add the .h5 extension if not given
        if not fname[-3:] == '.h5':
            fname += '.h5'

        # TODO: append to final h5 (instead of writting new file)

        # Save all the internal data to the h5 datafile named 'fname'
        with h5py.File(fname, 'w') as hf:
            hf.create_dataset('photon_energies', data=self.photon_energy)
            hf.create_dataset('cross_sections', data=self.cross_sections)
            hf.create_dataset('emission_energies', data=self.emission_enes)
            hf.create_dataset('energy_losses', data=self.ene_losses)
            hf.create_dataset('elastic_cross_sections', data=self.elastic_cross_sections)
            hf.create_dataset('scattering_amplitudes', data=self.scattering_amplitudes)

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
            """
            Join indices.
            """
            # TODO: double check
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

        if self._orb_and_state_dict:
            od = self._orb_and_state_dict

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

            state_list = _format_indices(od["valence_states"])
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

            orbital_labels = _format_orbital_names(od["mo_valence_indices"], "valence", nocc)
            orbital_str    = _format_indices(orbital_labels)
            self.ostream.print_header(
                f'Valence orbitals                 : {orbital_str}'.ljust(str_width)
            )
            
            orbital_labels = _format_orbital_names(od["mo_virtual_indices"], "virtual", nocc)
            orbital_str    = _format_indices(orbital_labels)
            self.ostream.print_header(
                f'Virtual orbitals                 : {orbital_str}'.ljust(str_width)
            )
            self.ostream.print_blank()

        if getattr(self, "_approach_string", None):
            self.ostream.print_info(self._approach_string)
            self.ostream.print_blank()
        
        self.ostream.flush()

    def _print_rixs_data(self, title, ene_losses,
                         cross_sections, elastic_cross_section):
        """
        Prints rixs-data to output stream.

        :param title:
            The title to be printed to the output stream.
        :param results:
            The dictionary containing response results.
        """

        spin_str = 'S'

        valstr = title
        self.ostream.print_header(valstr.ljust(92))
        self.ostream.print_header(('-' * len(valstr)).ljust(92))
        valstr = 'Ground State  {:>5s}: '.format(spin_str + str(0))
        valstr += '{:15.8f} a.u. '.format(0)
        valstr += '{:12.5f} eV'.format(0)
        valstr += '    Cross-section   {:9.2e}'.format(elastic_cross_section)
        self.ostream.print_header(valstr.ljust(92))
    
        for s, e in enumerate(ene_losses):
            valstr = 'Excited State {:>5s}: '.format(spin_str + str(s + 1))
            valstr += '{:15.8f} a.u. '.format(e)
            valstr += '{:12.5f} eV'.format(e * hartree_in_ev())
            f = cross_sections[s]
            valstr += '    Cross-section   {:9.2e}'.format(f)
            self.ostream.print_header(valstr.ljust(92))
        self.ostream.print_blank()
        self.ostream.flush()