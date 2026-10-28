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
import sys

from .veloxchemlib import hartree_in_ev, mpi_master, fine_structure_constant
from .oneeints import compute_electric_dipole_integrals
from .spectrumplot import plot_rixs_spectrum, plot_rixs_map
from .outputstream import OutputStream
from .linearsolver import LinearSolver
from .lreigensolver import LinearResponseEigenSolver
from .tdaeigensolver import TdaEigenSolver
from .cppsolver import ComplexResponseSolver
from .errorhandler import assert_msg_critical
from .resultsio import clear_group_in_hdf5
from .resultsio import write_rsp_results_to_hdf5


class RixsDriver(LinearSolver):
    """
    Implements the RIXS driver in a linear-response framework with two
    approximations: the two-shot and restricted-subspace approximation.

    # vlxtag: RHF, RIXS
    # vlxtag: RKS, RIXS

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - photon_energy: Incoming photon energies (a.u.).
        - theta: Angle between incident polarization
          and outgoing propagation vectors (rad.).
        - gamma: Life-time broadening (FWHM) (a.u.).
        - final_state_cutoff: Energy window of final states to include (a.u.).
        - tamm_dancoff: Use the Tamm-Dancoff approximation.
        - nstates: Number of valence-excited states.
        - num_core_states: Number of core-excited states.
        - restricted_subspace: Use the restricted-subspace approximation.
        - num_core_orbitals: Number of involved core orbitals.
        - num_valence_orbitals: Number of involved valence orbitals.
        - num_virtual_orbitals: Number of involved virtual orbitals.
        - twoshot: Whether the two-shot approach is used
          (read-only; derived from restricted_subspace).
        - cpp: Use the CPP routes: hybrid CPP-SOS, or 2x-CPP when
          energyloss_frequencies are given.
        - components: Cartesian components of the incoming photon
          (CPP routes).
        - energyloss_frequencies: Energy-loss frequencies (a.u.) of the
          2x-CPP route.
        - final_state_gamma: Final-state broadening (FWHM) (a.u.) of the
          2x-CPP route.
        - core_resonances: Number of core resonances used as incoming
          photon energies (SOS route).
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

        # For RIXS, the default convergence threshold is tighter
        self.conv_thresh = 2.0e-5

        self.tamm_dancoff = False

        # method settings
        self.photon_energy = None
        self.core_resonances = None
        self.theta = 0
        # a.u., FWHM
        # 1000.0 / hartree_in_wavenumber()
        self.gamma = 0.124 / hartree_in_ev()
        self.final_state_gamma = 0.124 / hartree_in_ev()

        self.nstates = 5
        self.num_core_states = 5

        # Restricted subspace
        self.restricted_subspace = False

        self.compute_in_mo = False
        # CPP-approach
        self.cpp = False
        self.components = 'xyz'
        self.energyloss_frequencies = None
        # TODO: enable emission freqs. in a good way
        self.emission_frequencies = None

        self.num_core_orbitals = 0
        self.num_virtual_orbitals = 0
        self.num_valence_orbitals = 0

        self._orb_and_state_dict = None

        # Set the energy range of the final states
        self.final_state_cutoff = None

        # Input keywords
        self._input_keywords['response'].update({
            'theta':
                ('float', 'angle between incident polarization vector and ' +
                    'outgoing propagation vector'),
            'gamma': ('float', 'broadening term (FWHM)'),
            'photon_energy': ('seq_range', 'list of incoming photon energies'),
            'core_resonances': ('int', 'number of core-resonances to take as incoming photon energy'),
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
            'compute_in_mo':
                ('bool', 'evaluate transition moments in the MO basis'),
            'cpp': ('bool', 'use the complex polarization propagator routes'),
            'components':
                ('str', 'Cartesian components of the incoming photon ' +
                    'polarization'),
            'energyloss_frequencies':
                ('seq_range', 'energy-loss frequencies for the double-CPP route'),
            'final_state_gamma':
                ('float', 'final-state broadening term (FWHM)'),
        })

    @property
    def twoshot(self):
        """
        Returns whether the two-shot approach is used.

        :return:
            Whether the two-shot approach is used.
        """

        return not self.restricted_subspace

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
        if key in rsp_dict and 'core_resonances' in rsp_dict:
            raise ValueError("Give either photon energies or core resonances, not both!")

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
            self.ostream.print_info(
                f'Number of states: (intermediate, final): '
                f'({intermed_states}, {final_states})')
            self.ostream.print_info(
                f'MO indices (core, valence, virtual): '
                f'({mo_c_ind}, {mo_val_ind}, {mo_vir_ind})')
            self.ostream.print_info(
                f'State indices (intermediate, final): '
                f'({ce_states}, {ve_states})')

    def compute(self, molecule, basis, scf_results,
                rsp_results=None, core_rsp_results=None):
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
        :core_rsp_results:
            The core-valence-separated (CVS) linear-response
            results dictionary.

        NOTE: To run full diagonalization, do a restricted
              subspace calculation with the full space,
              indicating which orbitals define the core.

        :return:
            Dictionary with cross sections, outgoing photon energy
            in (energy loss and full energy (emission)),
            and the scattering amplitude tensor.
        """

        if isinstance(self.photon_energy, (float, int, np.floating)):
            self.photon_energy = [self.photon_energy]

        if self.cpp and (self.energyloss_frequencies is not None or self.emission_frequencies is not None):
            return self.compute_cpp(molecule, basis, scf_results,
                             rsp_results=rsp_results, core_rsp_results=core_rsp_results)
        elif self.cpp and (self.energyloss_frequencies is None and self.emission_frequencies is None) and self.nstates is not None:
            return self.compute_hybrid(molecule, basis, scf_results,
                            rsp_results=rsp_results, core_rsp_results=core_rsp_results)
        else:
            return self.compute_sos(molecule, basis, scf_results,
                            rsp_results=rsp_results, core_rsp_results=core_rsp_results)

    def compute_sos(self, molecule, basis, scf_results,
                    rsp_results=None, core_rsp_results=None):
        """
        Computes RIXS properties.

        :param molecule:
            The molecule object.
        :param basis:
            The AO basis set.
        :param scf_results:
            The results dictionary from converged SCF wavefunction.
        :param rsp_results:
            The linear-response results dictionary.
        :param core_rsp_results:
            The core-valence-separated (CVS) linear-response
            results dictionary.

        NOTE: To run full diagonalization, do a restricted
              subspace calculation with the full space,
              indicating which orbitals define the core.

        :return:
            Dictionary containing the inelastic and elastic RIXS
            cross sections, incoming photon energies, energy losses,
            emission energies, scattering amplitude tensors,
            core-excitation energies and oscillator strengths, the
            lifetime broadening in eV, and the response type.
        """

        assert_msg_critical(
            not basis.has_ecp(),
            f'{type(self).__name__}.compute: ECPs are not supported for ' +
            'RIXS calculations. Explicit core orbitals are required.')

        self._is_converged = False
        solver_convergence = []

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
            solver_convergence.append(rsp_drv.is_converged)

            if not rsp_drv.is_converged:
                return {}

        if core_rsp_results is None and (not self.restricted_subspace):

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

            core_rsp_results = cvs_rsp_drv.compute(molecule, basis, scf_results)
            solver_convergence.append(cvs_rsp_drv.is_converged)

            if not cvs_rsp_drv.is_converged:
                return {}

        nocc = molecule.number_of_alpha_occupied_orbitals(basis)

        init_photon_set = True

        if self.rank == mpi_master():
            num_vir_orbitals = rsp_results['num_virtual']
        else:
            num_vir_orbitals = None
        num_vir_orbitals = self.comm.bcast(num_vir_orbitals, root=mpi_master())

        if self.twoshot:
            if self.rank == mpi_master():
                num_core_orbitals = core_rsp_results['num_core']
                num_valence_orbitals = nocc - num_core_orbitals
                num_intermediate_states = len(core_rsp_results['eigenvalues'])
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

            occupied_core = num_core_orbitals
            core_states = list(range(num_intermediate_states))
            val_states = list(range(num_final_states))

        else:
            if self.rank == mpi_master():
                num_valence_orbitals = rsp_results['num_valence']
                num_core_orbitals = rsp_results['num_core']
                assert_msg_critical(num_core_orbitals > 0,
                                    'No core orbitals indicated in the response results.')

                # Identify the energy of the lowest core-excited state
                first_core_ene = self._first_core_energy(rsp_results)
                detuning = rsp_results['eigenvalues'] - first_core_ene

                # Identify (and remove unphysical valence-excited states) the core-excited states
                core_states = self._core_state_indices(rsp_results, detuning)
                num_intermediate_states = len(core_states)
                assert_msg_critical(num_intermediate_states > 0,
                                    'Too few excited states included in response calculation.')
                # Identify the valence-excited states
                val_states = self._valence_state_indices(detuning, rsp_results['eigenvalues'])
                num_final_states = len(val_states)
            else:
                num_valence_orbitals = None
                num_core_orbitals = None
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

            # For compatibility with the two-shot approach
            core_rsp_results = rsp_results
            occupied_core = num_core_orbitals + num_valence_orbitals

        mo_core_indices = list(range(num_core_orbitals))
        mo_val_indices = list(range(nocc - num_valence_orbitals, nocc))
        mo_vir_indices = list(range(nocc, nocc + num_vir_orbitals))

        if self.rank == mpi_master():
            mo_occ = scf_results['C_alpha'][:, mo_core_indices + mo_val_indices].copy()
            mo_vir = scf_results['C_alpha'][:, mo_vir_indices].copy()
            mo = np.hstack((mo_occ, mo_vir))

            core_eigvals = core_rsp_results['eigenvalues'][core_states].copy()
            valence_eigvals = rsp_results['eigenvalues'][val_states].copy()
            core_osc_strengths = core_rsp_results['oscillator_strengths'][core_states].copy()
        else:
            mo_occ = None
            mo_vir = None

            core_eigvals = None
            valence_eigvals = None
            core_osc_strengths = None

        mo_occ = self.comm.bcast(mo_occ, root=mpi_master())
        mo_vir = self.comm.bcast(mo_vir, root=mpi_master())

        core_eigvals = self.comm.bcast(core_eigvals, root=mpi_master())
        valence_eigvals = self.comm.bcast(valence_eigvals, root=mpi_master())
        core_osc_strengths = self.comm.bcast(core_osc_strengths, root=mpi_master())

        core_eigvecs = self._get_eigvecs(core_rsp_results, core_states)
        valence_eigvecs = self._get_eigvecs(rsp_results, val_states)

        core_mats = self._preprocess_core_eigvecs(core_eigvecs, occupied_core,
                                                  num_valence_orbitals, num_vir_orbitals)

        if self.rank == mpi_master():
            dipole_integrals_ao = compute_electric_dipole_integrals(molecule, basis, [0.0,0.0,0.0])
            # perhaps this should be parallelized
            dipole_integrals_mo = np.array([mo.T @ dipole_integrals_ao[i] @ mo for i in range(3)])
            d_ov = np.array([mo_occ.T @ dipole_integrals_ao[i] @ mo_vir for i in range(3)])
            d_oo = np.array([mo_occ.T @ dipole_integrals_ao[i] @ mo_occ for i in range(3)])
            d_vv = np.array([mo_vir.T @ dipole_integrals_ao[i] @ mo_vir for i in range(3)])
        else:
            dipole_integrals_ao = None
            dipole_integrals_mo = None
            d_ov = None
            d_oo = None
            d_vv = None
        dipole_integrals_ao = self.comm.bcast(dipole_integrals_ao, root=mpi_master())
        dipole_integrals_mo = self.comm.bcast(dipole_integrals_mo, root=mpi_master())
        # d_ov = self.comm.bcast(d_ov, root=mpi_master())
        # d_oo = self.comm.bcast(d_oo, root=mpi_master())
        # d_vv = self.comm.bcast(d_vv, root=mpi_master())

        # Store state and orbital information used in computation
        self._orb_and_state_dict = {
            'num_intermediate_states': num_intermediate_states,
            'num_final_states': num_final_states,
            'mo_core_indices': mo_core_indices,
            'mo_valence_indices': mo_val_indices,
            'mo_virtual_indices': mo_vir_indices,
            'core_states': core_states,
            'valence_states': val_states,
        }

        if self.core_resonances is not None:
            # Take incoming photon energies as those of the first n core-excited state eigenvalues
            if self.rank == mpi_master():
                assert_msg_critical(self.core_resonances <= len(core_eigvals),
                                    f'Requested {self.core_resonances} core resonances, '
                                    f'but only {len(core_eigvals)} core states are available.')

                photon_list = core_eigvals[:self.core_resonances].tolist()

            else:
                photon_list = None

            photon_list = self.comm.bcast(photon_list, root=mpi_master())
            self.photon_energy = photon_list

        # If incoming photon energy is not set, set it to match the
        # first core-excited state with osc_strength > 1e-3
        if self.photon_energy is None:
            init_photon_set = False

            if self.rank == mpi_master():
                if core_rsp_results is None:
                    osc_arr = rsp_results['oscillator_strengths']
                else:
                    osc_arr = core_rsp_results['oscillator_strengths']
            else:
                osc_arr = None
            osc_arr = self.comm.bcast(osc_arr, root=mpi_master())

            for eig, osc in zip(core_eigvals, osc_arr[core_states]):
                if osc > 1e-3:
                    self.photon_energy = [eig]
                    break
        elif isinstance(self.photon_energy, (float, int, np.floating)):
            self.photon_energy = [self.photon_energy]

        assert_msg_critical(
            self.photon_energy is not None,
            'RixsDriver.compute: Unable to determine an incident photon '
            'energy automatically. Set photon_energy explicitly or include '
            'at least one core-excited state with oscillator strength > 1e-3.')
        assert_msg_critical(
            len(self.photon_energy) > 0,
            'RixsDriver.compute: photon_energy must contain at least one '
            'value.')

        photon_energy_arr = np.asarray(self.photon_energy, dtype=float)
        assert_msg_critical(
            np.all(np.isfinite(photon_energy_arr)),
            'RixsDriver.compute: photon_energy must contain only finite '
            'numeric values.')
        assert_msg_critical(
            np.all(np.abs(photon_energy_arr) > 0.0),
            'RixsDriver.compute: photon_energy must not contain zero, since '
            'the RIXS prefactor divides by the incident photon energy.')
        self.photon_energy = photon_energy_arr.tolist()

        # Define the return objects
        self.ene_losses = np.zeros((num_final_states, len(self.photon_energy)))
        self.emission_enes = np.zeros((num_final_states, len(self.photon_energy)))
        self.cross_sections = np.zeros((num_final_states, len(self.photon_energy)))
        self.elastic_cross_sections = np.zeros((len(self.photon_energy)))
        self.scattering_amplitudes = np.zeros((num_final_states, len(self.photon_energy),
                                               3, 3), dtype=complex)

        if self.rank == mpi_master():
            self._print_header(molecule, basis)

        ave, rem = divmod(num_final_states, self.nodes)
        counts = [ave + 1 if p < rem else ave for p in range(self.nodes)]
        f_start = sum(counts[:self.rank])
        f_end = sum(counts[:self.rank + 1])
        n_local = f_end - f_start

        for w_ind, omega in enumerate(self.photon_energy):
            F_elastic_local = np.zeros((3, 3), dtype=complex)
            local_cross_sections = np.empty(n_local)
            local_amplitude_tensors = np.empty((n_local, 3, 3), dtype=complex)
            local_emission_energies = np.empty(n_local)
            local_energy_losses = np.empty(n_local)

            for local_i, f in enumerate(range(f_start, f_end)):
                F_inelastic = np.zeros((3, 3), dtype=complex)
                z_val, y_val = self.split_eigvec(valence_eigvecs[f],
                                                 num_core_orbitals + num_valence_orbitals,
                                                 num_vir_orbitals, self.tamm_dancoff)
                for n in range(num_intermediate_states):
                    z_core, y_core = core_mats[n]
                    if self.compute_in_mo:
                        gs2core, core2val = self.get_tdms_mo(z_val, y_val, z_core, y_core)
                    else:
                        gs2core, core2val = self.get_tdms(mo_occ, mo_vir, z_val, y_val, z_core, y_core)

                    F_inelastic += self.scattering_amplitude_tensor(
                        omega, core_eigvals[n], valence_eigvals[f], gs2core, core2val, dipole_integrals_mo if self.compute_in_mo else dipole_integrals_ao
                    )
                    if f == 0:
                        F_elastic_local += self.scattering_amplitude_tensor(
                            omega, core_eigvals[n], None, gs2core, None, dipole_integrals_mo if self.compute_in_mo else dipole_integrals_ao
                        )

                emission_energy = omega - valence_eigvals[f]
                energy_loss = valence_eigvals[f]
                prefactor = emission_energy / omega

                # Fill local result slices
                local_emission_energies[local_i] = emission_energy
                local_energy_losses[local_i] = energy_loss
                local_cross_sections[local_i] = self.cross_section(F_inelastic, prefactor).real
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
                    self.emission_enes[start:end, w_ind] = em_en_part
                    self.ene_losses[start:end, w_ind] = loss_part
                    self.cross_sections[start:end, w_ind] = cs_part
                    self.scattering_amplitudes[start:end, w_ind] = amp_part

                self.ostream.print_info(
                    f'Computed RIXS cross sections for {num_final_states} final states '
                    f'at photon energy: {omega*hartree_in_ev():.2f} eV.'
                )
                self.ostream.print_blank()
                self._print_rixs_data(f'RIXS cross sections at incident X-ray energy '
                                      f'{omega*hartree_in_ev():.2f} eV, energy-loss mode',
                                      self.ene_losses[:,w_ind], self.cross_sections[:,w_ind], self.elastic_cross_sections[w_ind])
                self.ostream.print_blank()
                self.ostream.flush()

        results_dict = {
            'cross_sections': self.cross_sections,
            'elastic_cross_sections': self.elastic_cross_sections,
            'photon_energies': np.array(self.photon_energy),
            'scattering_amplitudes': self.scattering_amplitudes,
            'emission_energies': self.emission_enes,
            'energy_losses': self.ene_losses,
            'core_eigenvalues': core_eigvals,
            'core_osc_strengths': core_osc_strengths,
            'gamma_fwhm_ev': self.gamma * hartree_in_ev(),
        }

        # add rsp type
        results_dict.update({'rsp_type': 'rixs'})

        if self.filename is not None and self.rank == mpi_master():
            final_h5_fname = f'{self.filename}.h5'
            if Path(final_h5_fname).is_file():
                clear_group_in_hdf5(final_h5_fname, 'rsp')
                write_rsp_results_to_hdf5(final_h5_fname, results_dict)
                valstr = f'RIXS results written to file: {final_h5_fname}'
                self.ostream.print_info(valstr)
                self.ostream.print_blank()
                self.ostream.flush()

        if not init_photon_set:
            self.photon_energy = None

        if solver_convergence:
            self._is_converged = all(solver_convergence)
        else:
            # Externally supplied response result dictionaries are trusted as
            # converged inputs RIXS.
            self._is_converged = True

        return results_dict

    def compute_cpp(self, molecule, basis, scf_results,
                    rsp_results=None, core_rsp_results=None):
        """
        Implements the 2x-CPP variant of LR-TDDFT RIXS.

        :param rsp_results:
            The CPP results at the energy-loss frequencies.
        :param core_rsp_results:
            The CPP results at the incoming photon energies.
        """

        # TODO: enable ECP
        assert_msg_critical(
            not basis.has_ecp(),
            f'{type(self).__name__}.compute: ECP is not yet supported')

        assert_msg_critical(self.photon_energy is not None,
                            'Must provide incoming photon energies for the CPP approach!')

        assert_msg_critical(self.energyloss_frequencies is not None,
                            'Must provide energy-loss frequencies for the 2x-CPP approach!')

        self.photon_energy = [float(w) for w in self.photon_energy]

        self._is_converged = False
        solver_convergence = []

        # set up the response driver for the final states
        # no lindep_thresh for final/valence manifold
        rsp_drv = self._setup_cpp_solver(self.energyloss_frequencies,
                                         self.final_state_gamma, '_rixs.h5')

        if self.rank == mpi_master():
            self._print_header(molecule, basis)

        if core_rsp_results is None:
            core_rsp_results, converged = self._solve_core_cpp(
                molecule, basis, scf_results, '_rixs_cvs.h5', 1.0e-8)
            solver_convergence.append(converged)

            if not converged:
                return {}

        n_comp = len(self.components)
        n_photons = len(self.photon_energy)
        num_final_freqs = len(self.energyloss_frequencies)

        self.ene_losses             = np.zeros((num_final_freqs, n_photons))
        self.emission_enes          = np.zeros((num_final_freqs, n_photons))
        self.cross_sections         = np.zeros((num_final_freqs, n_photons))
        self.elastic_cross_sections = np.zeros((n_photons))
        self.rsp_tensor = np.zeros((num_final_freqs, n_photons, n_comp, n_comp, n_comp, n_comp), dtype=complex)

        nocc = molecule.number_of_alpha_occupied_orbitals(basis)

        if self.rank == mpi_master():
            C = scf_results['C_alpha']
            dip = compute_electric_dipole_integrals(molecule, basis, [0.0, 0.0, 0.0])
            dip_ints_mo = {comp: C.T @ dip_x_ao @ C for comp, dip_x_ao in zip('xyz', dip)}
        else:
            dip_ints_mo = None
        dip_ints_mo = self.comm.bcast(dip_ints_mo, root=mpi_master())

        norb = dip_ints_mo['x'].shape[0]

        # TODO: make stable for len(self.photon_energy) = 1
        for n_ind, w_ex in enumerate(self.photon_energy):
            self.ostream.print_info(
                    f'Computing cross sections for components: {self.components}, of {num_final_freqs} energy-loss frequencies '
                    f'at photon energy: {w_ex*hartree_in_ev():.2f} eV.'
                )
            self.ostream.print_blank()

            # TODO: enable elastic cross-sections

            # CPP solutions at the incoming photon energy, one per component,
            # at w_ex for incident component a
            A = {
                a: ComplexResponseSolver.get_full_solution_vector(
                    core_rsp_results['solutions'][(a, w_ex)]
                )
                for a in self.components
            }
            A = self.comm.bcast(A, root=mpi_master())

            # pair order (a, b)
            # a = incident, b = outgoing
            pairs = [(a, b) for a in self.components for b in self.components]

            # eta[(a, b)] = -[mu_b, A_a]: RHS of the second leg and bra of M
            eta = {}
            for a in self.components:
                grads = self._eff_excited_grad(
                    A[a], dip_ints_mo, nocc, norb, components=self.components
                )

                for b, grad in zip(self.components, grads):
                    eta[(a, b)] = grad

            B = {}
            for a in self.components:
                vec_grad = {(b, w): eta[(a, b)]
                            for b in self.components 
                            for w in rsp_drv.frequencies
                            }

                B[a] = rsp_drv.compute(
                    molecule, basis, scf_results, v_grad=vec_grad
                )
                solver_convergence.append(rsp_drv.is_converged)
                
            """
            # slower but cleaner alternative to the above
            #vec_grad = {(f'{a}{b}', w): eta[(a, b)] 
            #            for (a, b) in pairs for w in rsp_drv.frequencies}

            B = rsp_drv.compute(
                molecule, basis, scf_results, v_grad=vec_grad
            )
            solver_convergence.append(rsp_drv.is_converged)
            """

            # (9,n)
            eta_mat = np.array([eta[ab] for ab in pairs])
            #eta_mat = np.array([
            #    eta[x][y]
            #    for y in self.components
            #    for x in self.components
            #])

            for f_ind, w in enumerate(rsp_drv.frequencies):
                w_em = w_ex - w

                B_mat = [
                    ComplexResponseSolver.get_full_solution_vector(
                        B[a]['solutions'][(b, w)])
                    for (a, b) in pairs
                ]

                if self.rank == mpi_master():
                    B_mat = np.array(B_mat)

                    # M[(a', b'), (a, b)] = <eta_{a'b'} | B_{ab}(w_loss)>
                    M = (np.conjugate(eta_mat) @ B_mat.T).reshape(n_comp, n_comp, n_comp, n_comp)
                    T = M.imag

                    # TODO: rewrite without einsum
                    # sum_xy |F_xy|^2
                    T1 = np.einsum('xyxy ->', T)
                    # sum_xy F_xy F^*_yx
                    T2 = np.einsum('xyyx ->', T)
                    # sum_xy F_xx F^*_yy
                    T3 = np.einsum('xxyy ->', T)

                    pref = fine_structure_constant()**4 * w_em**3 * w_ex / np.pi
                    sigma = pref * (1.0/15.0) * (
                                (2.0 - 0.5 * np.sin(self.theta)**2) * T1 +
                                (0.75 * np.sin(self.theta)**2 - 0.5) * (T2 + T3)
                                )
                else:
                    M = None
                    sigma = None

                M = self.comm.bcast(M, root=mpi_master())
                sigma = self.comm.bcast(sigma, root=mpi_master())

                self.ene_losses[f_ind, n_ind] = w
                self.emission_enes[f_ind, n_ind] = w_em
                self.cross_sections[f_ind, n_ind] = sigma
                self.rsp_tensor[f_ind, n_ind] = M

            if self.rank == mpi_master():
                self.ostream.print_info(
                    f'Computed RIXS cross sections for {num_final_freqs} energy-loss frequencies '
                    f'at photon energy: {w_ex*hartree_in_ev():.2f} eV.'
                )
                self.ostream.print_blank()
                self._print_rixs_data(f'RIXS cross sections at incident X-ray energy '
                                      f'{w_ex*hartree_in_ev():.2f} eV, energy-loss mode',
                                      self.ene_losses[:,n_ind], self.cross_sections[:,n_ind], self.elastic_cross_sections[n_ind])
                self.ostream.print_blank()
                self.ostream.flush()

        results_dict = {
            'cross_sections': self.cross_sections,
            'elastic_cross_sections': self.elastic_cross_sections,
            'photon_energies': np.array(self.photon_energy),
            'rsp_tensor': self.rsp_tensor,
            'emission_energies': self.emission_enes,
            'energy_losses': self.ene_losses,
            # no explicit core-excited states in the CPP routes
            'core_eigenvalues': np.array([]),
            'core_osc_strengths': np.array([]),
            'gamma_fwhm_ev': self.gamma * hartree_in_ev(),
            'final_state_gamma_fwhm_ev': self.final_state_gamma * hartree_in_ev(),
            'rsp_type': 'rixs',
        }

        self._write_rixs_results(results_dict)

        if self.rank == mpi_master():
            self.ostream.print_info('...done.')
            self.ostream.print_blank()
            self.ostream.flush()

        if solver_convergence:
            self._is_converged = all(solver_convergence)
        else:
            self._is_converged = True

        return results_dict

    def compute_hybrid(self, molecule, basis, scf_results,
                       rsp_results=None, core_rsp_results=None):
        """
        Implements the hybrid CPP-SOS variant of LR-TDDFT RIXS: the sum over
        core-excited states is replaced by the CPP solutions at the incoming
        photon energies, while the final states are explicit eigenpairs.

        :param rsp_results:
            The valence linear-response results dictionary (optional).
        :param core_rsp_results:
            The CPP results at the incoming photon energies (optional).
        """

        assert_msg_critical(
            not basis.has_ecp(),
            f'{type(self).__name__}.compute: ECP is not yet supported')

        assert_msg_critical(not self.tamm_dancoff,
                            'Tamm-Dancoff approximation is not supported in the hybrid approach!')

        assert_msg_critical(self.photon_energy is not None,
                            'Must provide incoming photon energies for the CPP approach!')

        self.photon_energy = [float(w) for w in self.photon_energy]

        self._is_converged = False
        solver_convergence = []

        if core_rsp_results is None:
            core_rsp_results, converged = self._solve_core_cpp(
                molecule, basis, scf_results, '_rixs_cpp.h5', 1.0e-10)
            solver_convergence.append(converged)

            if not converged:
                return {}

        if rsp_results is None:

            rsp_keys = [
                'nstates',
                'conv_thresh',
                'restart',
            ]

            rsp_drv = LinearResponseEigenSolver(self.comm, self.ostream)

            for key in rsp_keys:
                if hasattr(self, key):
                    setattr(rsp_drv, key, getattr(self, key))

            if self.checkpoint_file is not None:
                fpath = Path(self.checkpoint_file)
                fpath = fpath.with_name(fpath.stem)
                rsp_drv.checkpoint_file = str(fpath) + '_rixs.h5'

            rsp_results = rsp_drv.compute(molecule, basis, scf_results)
            solver_convergence.append(rsp_drv.is_converged)

            if not rsp_drv.is_converged:
                return {}

        if self.rank == mpi_master():
            self._print_header(molecule, basis)

        nocc = molecule.number_of_alpha_occupied_orbitals(basis)

        if self.rank == mpi_master():
            C = scf_results['C_alpha']
            valence_eigvals = rsp_results['eigenvalues'].copy()
            dipole_integrals_ao = np.array(compute_electric_dipole_integrals(molecule, basis, [0.0,0.0,0.0]))
        else:
            C = None
            valence_eigvals = None
            dipole_integrals_ao = None
        C = self.comm.bcast(C, root=mpi_master())
        valence_eigvals = self.comm.bcast(valence_eigvals, root=mpi_master())
        dipole_integrals_ao = self.comm.bcast(dipole_integrals_ao, root=mpi_master())

        norb = C.shape[1]
        nvir = norb - nocc

        mo_occ = C[:, :nocc].copy()
        mo_vir = C[:, nocc:].copy()
        dipole_integrals_mo = np.array([C.T @ dipole_integrals_ao[i] @ C for i in range(3)])

        components = self.components

        num_final_states = valence_eigvals.size
        valence_eigvecs = self._get_eigvecs(rsp_results, list(range(num_final_states)))

        # Define the return objects
        self.ene_losses             = np.zeros((num_final_states, len(self.photon_energy)))
        self.emission_enes          = np.zeros((num_final_states, len(self.photon_energy)))
        self.cross_sections         = np.zeros((num_final_states, len(self.photon_energy)))
        self.elastic_cross_sections = np.zeros((len(self.photon_energy)))
        self.scattering_amplitudes  = np.zeros((num_final_states, len(self.photon_energy),
                                               3, 3), dtype=complex)

        # TODO: make stable for len(self.photon_energy) = 1
        for n_ind, w_ex in enumerate(self.photon_energy):

            # CPP solutions (Z, Y) at the incoming photon energy, one per component
            core_mats = [ComplexResponseSolver.get_full_solution_vector(core_rsp_results['solutions'][(comp, w_ex)]) for comp in components]
            core_mats = [self.comm.bcast(eigvec, root=mpi_master()) for eigvec in core_mats]
            core_vecs = [(eigvec[:eigvec.size // 2].reshape(nocc, -1), eigvec[eigvec.size // 2:].reshape(nocc, -1)) for eigvec in core_mats]

            for f_ind, w in enumerate(valence_eigvals):

                val_vec = valence_eigvecs[f_ind]

                z_val = val_vec[:val_vec.size // 2].reshape(nocc, -1)
                y_val = val_vec[val_vec.size // 2:].reshape(nocc, -1)

                omega_product = (w_ex - w) * w_ex

                # F[x, y]: x = emitted, y = incoming component; the minus sign
                # on the transition density matches the SOS amplitude convention
                F = np.zeros((3, 3), dtype=complex)
                for comp, (z_core, y_core) in zip(components, core_vecs):
                    if self.compute_in_mo:
                        core_to_val = self.get_tdms_mo(z_val, y_val, z_core, y_core)[1]
                        dipole_integrals = dipole_integrals_mo
                    else:
                        core_to_val = self._get_esa_transition_density(z_core, y_core, z_val, y_val, mo_occ, mo_vir)
                        dipole_integrals = dipole_integrals_ao

                    A = dipole_integrals.reshape(dipole_integrals.shape[0], -1)
                    B = (-core_to_val).reshape(-1)
                    F[:, 'xyz'.index(comp)] = omega_product * np.matmul(A, B)
                self.scattering_amplitudes[f_ind, n_ind] = F

                omegaprim = (w_ex - w) / w_ex
                xsec = self.cross_section(F, omegaprime_omega=omegaprim)
                self.cross_sections[f_ind, n_ind] = xsec

                self.ene_losses[f_ind, n_ind] = w
                self.emission_enes[f_ind, n_ind] = w_ex - w

            if self.rank == mpi_master():
                self.ostream.print_info(
                    f'Computed RIXS cross sections for {num_final_states} energy-loss frequencies '
                    f'at photon energy: {w_ex*hartree_in_ev():.2f} eV.'
                )
                self.ostream.print_blank()
                self._print_rixs_data(f'RIXS cross sections at incident X-ray energy '
                                      f'{w_ex*hartree_in_ev():.2f} eV, energy-loss mode',
                                      self.ene_losses[:,n_ind], self.cross_sections[:,n_ind], self.elastic_cross_sections[n_ind])
                self.ostream.print_blank()
                self.ostream.flush()

        results_dict = {
            'cross_sections': self.cross_sections,
            'elastic_cross_sections': self.elastic_cross_sections,
            'photon_energies': np.array(self.photon_energy),
            'scattering_amplitudes': self.scattering_amplitudes,
            'emission_energies': self.emission_enes,
            'energy_losses': self.ene_losses,
            # no explicit core-excited states in the CPP routes
            'core_eigenvalues': np.array([]),
            'core_osc_strengths': np.array([]),
            'gamma_fwhm_ev': self.gamma * hartree_in_ev(),
            'rsp_type': 'rixs',
        }

        self._write_rixs_results(results_dict)

        if self.rank == mpi_master():
            self.ostream.print_info('...done.')
            self.ostream.print_blank()
            self.ostream.flush()

        if solver_convergence:
            self._is_converged = all(solver_convergence)
        else:
            self._is_converged = True

        return results_dict

    def _setup_cpp_solver(self, frequencies, gamma_fwhm, checkpoint_suffix,
                          lindep_thresh=None):
        """
        Sets up a CPP solver for the given frequencies and FWHM broadening.
        """

        cpp_drv = ComplexResponseSolver(self.comm, self.ostream)

        for key in ['conv_thresh', 'restart']:
            if hasattr(self, key):
                setattr(cpp_drv, key, getattr(self, key))

        cpp_drv.lindep_thresh = lindep_thresh
        cpp_drv.frequencies = [float(w) for w in frequencies]
        # given broadening is FWHM but CPP-driver takes HWHM -> HWHM
        cpp_drv.damping = gamma_fwhm / 2

        if self.checkpoint_file is not None:
            fpath = Path(self.checkpoint_file)
            fpath = fpath.with_name(fpath.stem)
            cpp_drv.checkpoint_file = str(fpath) + checkpoint_suffix

        return cpp_drv

    def _solve_core_cpp(self, molecule, basis, scf_results, checkpoint_suffix,
                        single_energy_lindep_thresh):
        """
        Solves the CPP equations at the incoming photon energies.
        """

        lindep_thresh = self.lindep_thresh

        if len(self.photon_energy) == 1 and lindep_thresh is None:
            # tighten the linear dependence threshold for single photon energy
            lindep_thresh = single_energy_lindep_thresh
            self.ostream.print_info(f'Single photon energy provided, tightening linear dependence threshold to {lindep_thresh:.1e}')

        core_rsp_drv = self._setup_cpp_solver(self.photon_energy, self.gamma,
                                              checkpoint_suffix, lindep_thresh)

        energies_ev = np.array(self.photon_energy) * hartree_in_ev()
        energies_str = ", ".join(f"{ene:.2f}" for ene in energies_ev)
        self.ostream.print_blank()
        self.ostream.print_info(
            f"Solving CPP-equations at incoming photon energies: {energies_str} eV."
        )
        self.ostream.print_blank()

        core_rsp_results = core_rsp_drv.compute(molecule, basis, scf_results)
        self.ostream.print_info('...done.')

        return core_rsp_results, core_rsp_drv.is_converged

    def _write_rixs_results(self, results_dict):
        """
        Writes the RIXS results to the final HDF5 file, if it exists.
        """

        if self.filename is not None and self.rank == mpi_master():
            final_h5_fname = f'{self.filename}.h5'
            if Path(final_h5_fname).is_file():
                clear_group_in_hdf5(final_h5_fname, 'rsp')
                write_rsp_results_to_hdf5(final_h5_fname, results_dict)
                valstr = f'RIXS results written to file: {final_h5_fname}'
                self.ostream.print_info(valstr)
                self.ostream.print_blank()
                self.ostream.flush()

    def _eff_excited_grad(self, response_vec, dip_ints_mo, nocc, norb, components='xyz'):
        """
        Get effective RHS/gradient vectors of a LR calculation.

        :param response_vec:
            The full response vector (Z, Y).
        :param dip_ints_mo:
            The dipole integrals in the MO basis, keyed by component.
        :param nocc:
            The number of occupied orbitals.
        :param norb:
            The number of orbitals.
        :param components:
            The Cartesian components.

        :return:
            Excited state gradient (tuple)
        """

        Z = self.lrvec2mat(response_vec, nocc, norb)
        Z[:nocc, nocc:] *= -1.0

        grads = []
        for comp in components:
            # dip_x_mo @ Z - Z @ dip_x_mo
            grad = self.commut(dip_ints_mo[comp], Z)

            grad_vec = -1.0 * self.lrmat2vec(grad, nocc, norb)
            grads.append(grad_vec)

        return tuple(grads)

    def _first_core_energy(self, rsp_results):
        """
        Return energy of the first core-excited state (by label order).

        :param rsp_results:
            Dictionary containing excitation information.

        :return:
            The energy (float) corresponding to the first
            state labeled as 'core' in excitation_details.
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
            List of core-excited state indices.
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
        :param eigenvalues:
            Array of excitation energies used
            for the final-state cutoff.
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
        :param core_eigenvalue:
            Energy of the intermediate (core) excited state.
        :param val_eigenvalue:
            Energy of the final (valence) excited state. If None,
            the elastic line is computed.
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
        Computes the RIXS cross section.

        Journal of Chemical Theory and Computation 2021 17 (5), 3031-3038
        DOI: 10.1021/acs.jctc.1c00144

        Journal of Chemical Theory and Computation 2017 13 (11), 5552-5559
        DOI: 10.1021/acs.jctc.7b00636

        :param F:
            The scattering amplitude tensor.
        :param omegaprime_omega:
            The ratio between outgoing- and incoming frequencies, w'/w.

        :return:
            The scattering cross section.
        """

        # F_xy (F^(xy))^*
        F2 = np.vdot(F, F).real
        # F_xy (F^(yx))^*
        FF_T_conj = np.vdot(F, F.T).real
        # F_xx (F^(yy))^*
        tr_F2 = np.abs(np.trace(F))**2

        sigma_f = (fine_structure_constant()**4) * omegaprime_omega * (1.0/15.0) * (
            (2.0 - 0.5*np.sin(self.theta)**2) * F2 +
            (0.75*np.sin(self.theta)**2 - 0.5) * (FF_T_conj + tr_F2))

        return sigma_f

    def _get_eigvecs(self, rsp_results, states):
        """
        Gets eigenvectors.

        :param rsp_results:
            The dictionary containing response results.
        :param states:
            The excited states.

        :return:
            The eigenvectors as a 2D numpy array with one row
            per state.
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
        :param pad_occ:
            The number of missing orbitals.

        :return:
            The padded, full-size matrices z and y.
        """
        padding = ((0, pad_occ), (0, 0))
        return np.pad(z_mat, padding), np.pad(y_mat, padding)

    def _preprocess_core_eigvecs(self, core_eigvecs, occ_core,
                                 num_val_orbs, num_vir_orbs):
        """
        Splits and, for the two-shot approach, zero-pads the core
        eigenvectors.

        :param core_eigvecs:
            The eigenvectors of the core-excited states.
        :param occ_core:
            Number of occupied orbitals in the core-excited
            response calculation.
        :param num_val_orbs:
            Number of valence orbitals.
        :param num_vir_orbs:
            Number of virtual orbitals.

        :return:
            The processed core eigenvectors as a list of
            (excitation, deexcitation) matrix pairs.
        """
        pp_core_eigvecs = []

        for vec in core_eigvecs:
            z, y = self.split_eigvec(vec, occ_core,
                                     num_vir_orbs, self.tamm_dancoff)

            if self.twoshot:
                z, y = self.pad_matrices(z, y, num_val_orbs)

            pp_core_eigvecs.append((z, y))

        return pp_core_eigvecs

    def get_tdms(self, mo_occ, mo_vir, z_val, y_val, z_core, y_core):
        """
        Get the transition density matrices, both from ground-state
        to excited state, and between two excited states (here between)
        core-excited and valence-excited states.

        :param mo_occ:
            The occupied molecular orbitals.
        :param mo_vir:
            The unoccupied/virtual molecular orbitals.
        :param z_val:
            The excitation matrix (valence-excited state).
        :param y_val:
            The dexcitation matrix (valence-excited state).
        :param z_core:
            The excitation matrix (core-excited state).
        :param y_core:
            The dexcitation matrix (core-excited state).

        :return:
            The ground-to-core and core-to-valence transition densities.
        """

        gs_to_core = np.linalg.multi_dot([mo_occ, z_core - y_core, mo_vir.T])
        gs_to_core *= np.sqrt(2.0)

        core_to_val = self._get_esa_transition_density(
            z_core, y_core, z_val, y_val, mo_occ, mo_vir)

        return gs_to_core, core_to_val

    def get_tdms_mo(self, z_val, y_val, z_core, y_core):
        """
        Get the transition density matrices, both from ground-state
        to excited state, and between two excited states (here between)
        core-excited and valence-excited states.

        :param mo_occ:
            The occupied molecular orbitals.
        :param mo_vir:
            The unoccupied/virtual molecular orbitals.
        :param z_val:
            The excitation matrix (valence-excited state).
        :param y_val:
            The dexcitation matrix (valence-excited state).
        :param z_core:
            The excitation matrix (core-excited state).
        :param y_core:
            The dexcitation matrix (core-excited state).

        :return:
            The ground-to-core and core-to-valence transition densities.
        """
        nocc, nvir = z_core.shape
        nmo = nocc + nvir

        # keep the dtype of complex (CPP) response vectors
        dtype = np.result_type(z_val, y_val, z_core, y_core)

        gs_to_core = np.zeros((nmo, nmo), dtype=dtype)
        gs_to_core[:nocc, nocc:] = np.sqrt(2.0) * (z_core - y_core)

        #core_to_val = self._get_esa_transition_density(
        #    z_core, y_core, z_val, y_val, mo_occ, mo_vir)
        core_to_val = np.zeros((nmo, nmo), dtype=dtype)
        core_to_val[nocc:, nocc:] = z_core.T @ z_val + y_core.T @ y_val
        core_to_val[:nocc, :nocc] = -(z_core @ z_val.T + y_core @ y_val.T)

        return gs_to_core, core_to_val

    def _print_header(self, molecule, basis):
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

        str_width = 60
        nocc = molecule.number_of_alpha_occupied_orbitals(basis)

        self.ostream.print_blank()
        title = 'Resonant Inelastic X-ray Scattering (RIXS) Setup'
        self.ostream.print_header(f'{title:^{str_width}}')
        self.ostream.print_header(f'{"=" * len(title):^{str_width}}')
        self.ostream.print_blank()

        cur_str = 'Scattering angle (theta) [rad]   : {:.2f}'.format(self.theta)
        self.ostream.print_header(cur_str.ljust(str_width))
        gamma_ev = self.gamma * hartree_in_ev()
        cur_str = 'Lifetime broadening (gamma) [eV] : {:.2f}'.format(gamma_ev)
        self.ostream.print_header(cur_str.ljust(str_width))

        if self.cpp:
            cur_str = 'Incoming photon components       : {:s}'.format(
                ''.join(self.components))
            self.ostream.print_header(cur_str.ljust(str_width))

        if self.cpp and self.energyloss_frequencies is not None:
            fs_gamma_ev = self.final_state_gamma * hartree_in_ev()
            cur_str = 'Final-state broadening [eV]      : {:.2f}'.format(
                fs_gamma_ev)
            self.ostream.print_header(cur_str.ljust(str_width))
            loss_ev = np.asarray(self.energyloss_frequencies,
                                 dtype=float) * hartree_in_ev()
            cur_str = 'Energy-loss frequencies [eV]     : {:.2f} ... {:.2f} ({:d} points)'.format(
                loss_ev.min(), loss_ev.max(), loss_ev.size)
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
            orbital_str = _format_indices(orbital_labels)
            self.ostream.print_header(
                f'Core orbitals                    : {orbital_str}'.ljust(str_width)
            )

            orbital_labels = _format_orbital_names(od["mo_valence_indices"], "valence", nocc)
            orbital_str = _format_indices(orbital_labels)
            self.ostream.print_header(
                f'Valence orbitals                 : {orbital_str}'.ljust(str_width)
            )

            orbital_labels = _format_orbital_names(od["mo_virtual_indices"], "virtual", nocc)
            orbital_str = _format_indices(orbital_labels)
            self.ostream.print_header(
                f'Virtual orbitals                 : {orbital_str}'.ljust(str_width)
            )
            self.ostream.print_blank()

        if self.cpp and self.energyloss_frequencies is not None:
            self.ostream.print_info(
                'Running RIXS calculation in the double-CPP approach')
            self.ostream.print_blank()
        elif self.cpp:
            self.ostream.print_info(
                'Running RIXS calculation in the two-shot hybrid (CPP-SOS) approach')
            self.ostream.print_blank()
        elif self.twoshot:
            self.ostream.print_info(
                'Running RIXS calculation in the two-shot approach')
            self.ostream.print_blank()
            self.ostream.print_reference(
                'J. Chem. Theory Comput. (2021) 17 (5): 3031-3038, DOI: 10.1021/acs.jctc.1c00144')
            self.ostream.print_reference(
                'J. Phys. Chem. A, 2025, 129, 8783-8797, DOI: 10.1021/acs.jpca.5c04528')
            self.ostream.print_blank()
        else:
            self.ostream.print_info(
                'Running RIXS calculation in the restricted-subspace approach')
            self.ostream.print_blank()
            self.ostream.print_reference(
                'Phys. Chem. Chem. Phys., 2021, 23, 1835-1848, DOI: 10.1039/D0CP04726K')
            self.ostream.print_reference(
                'J. Phys. Chem. A, 2025, 129, 8783-8797, DOI: 10.1021/acs.jpca.5c04528')
            self.ostream.print_blank()

        self.ostream.flush()

    def _print_rixs_data(self, title, ene_losses,
                         cross_sections, elastic_cross_section):
        """
        Prints rixs-data to output stream.

        :param title:
            The title to be printed to the output stream.
        :param ene_losses:
            The list of energy losses.
        :param cross_sections:
            The list of cross sections.
        :param elastic_cross_section:
            The elastic cross section.
        """

        spin_str = 'S'

        valstr = title
        self.ostream.print_header(valstr.ljust(92))
        self.ostream.print_header(('-' * len(valstr)).ljust(92))
        valstr = 'Ground State  {:>5s}: '.format(spin_str + str(0))
        valstr += '{:15.8f} a.u. '.format(0)
        valstr += '{:12.5f} eV'.format(0)
        valstr += '    Cross section   {:9.2e}'.format(elastic_cross_section)
        self.ostream.print_header(valstr.ljust(92))

        for s, e in enumerate(ene_losses):
            valstr = 'Excited State {:>5s}: '.format(spin_str + str(s + 1))
            valstr += '{:15.8f} a.u. '.format(e)
            valstr += '{:12.5f} eV'.format(e * hartree_in_ev())
            f = cross_sections[s]
            valstr += '    Cross section   {:9.2e}'.format(f)
            self.ostream.print_header(valstr.ljust(92))
        self.ostream.print_blank()
        self.ostream.flush()

    def plot(self,
             results,
             broadening_type="lorentzian",
             broadening_value_ev=0.24,
             x_unit='ev',
             xstep=1e-4,
             energy_loss=True,
             photon_energy_ev=None,
             plot_elastic_line=False,
             ax=None):
        """
        Plot the RIXS spectrum.

        :param results:
            The results dictionary from compute method.
        :param broadening_type:
            The type of broadening to use. Either 'lorentzian' or 'gaussian'.
        :param broadening_value_ev:
            The broadening value, FWHM in eV.
        :param x_unit:
            Unit for the x-axis. Currently only 'ev' is supported.
            NOTE: does not affect other inputs, which are expected to be in eV.
        :param xstep:
            Grid spacing in eV for the broadened RIXS spectrum.
        :param energy_loss:
            If True, plot versus energy loss.
            If False, plot versus emission energy.
        :param photon_energy_ev:
            Incoming photon energy in eV. The closest
            available calculated photon energy is used.
        :param plot_elastic_line:
            If True, include the elastic line in the plot.
        :param ax:
            The matplotlib axis to plot on.

        """

        if photon_energy_ev is not None:
            # Find the closest calculated photon energy index
            incoming_photon_energies_ev = np.asarray(results['photon_energies']) * hartree_in_ev()
            photon_index = int(np.argmin(np.abs(incoming_photon_energies_ev - float(photon_energy_ev))))
            closest_energy_ev = incoming_photon_energies_ev[photon_index]

            if abs(closest_energy_ev - photon_energy_ev) > 0.1:
                self.ostream.print_warning(
                    f'{type(self).__name__}: requested photon energy {photon_energy_ev} eV is more '
                    f'than 0.1 eV far from the closest calculated photon energy {closest_energy_ev:.2f} eV.'
                )
            self.ostream.flush()
        else:
            photon_index = 0

        plot_rixs_spectrum(
            results,
            broadening_type=broadening_type,
            broadening_value=broadening_value_ev,
            photon_index=photon_index,
            x_unit=x_unit,
            xstep=xstep,
            energy_loss=energy_loss,
            plot_elastic_line=plot_elastic_line,
            ax=ax
        )

    def plot_map(self,
                 results,
                 broadening_type="lorentzian",
                 broadening_value_ev=0.24,
                 xstep=1e-2,
                 x_unit='ev',
                 energy_loss=True,
                 min_photon_points=10,
                 colormap='viridis',
                 normalize='global_max',
                 plot_elastic_line=False,
                 ax=None):
        """
        Plot the 2D RIXS map together with the corresponding XAS spectrum.

        :param results:
            The RIXS results dictionary from the compute method.
        :param broadening_type:
            The type of broadening to use. Either 'lorentzian' or 'gaussian'.
        :param broadening_value_ev:
            The broadening value, FWHM in eV.
        :param xstep:
            Grid spacing in eV for the broadened RIXS map.
        :param x_unit:
            Unit for the x-axis. Currently only 'ev' is supported.
            NOTE: does not affect other inputs, which are expected to be in eV.
        :param energy_loss:
            If True, plot versus energy loss.
            If False, plot versus emission energy.
        :param min_photon_points:
            Minimum number of incoming photon energies required for a meaningful map.
        :param colormap:
            The colormap to use for the RIXS map.
        :param normalize:
            Normalization method for the RIXS map.
            'global_max': normalize by the global maximum of the RIXS intensities -> range of intensity [0-1] (default).
        :param plot_elastic_line:
            If True, include the elastic line in the RIXS map.
        :param ax:
            If None, new axes are created.
            Otherwise pass a tuple (ax_map, ax_xas).
        """

        # incoming photon energies
        nr_incoming_photons = len(np.asarray(results['photon_energies']))

        assert_msg_critical(
            nr_incoming_photons >= min_photon_points,
            f'{type(self).__name__}: too few incoming photon energies for a meaningful 2D RIXS map.'
            f' (min: {min_photon_points}, got: {nr_incoming_photons})'
        )

        plot_rixs_map(
            results,
            broadening_type=broadening_type,
            broadening_value=broadening_value_ev,
            xstep=xstep,
            x_unit=x_unit,
            energy_loss=energy_loss,
            cmap=colormap,
            normalize=normalize,
            plot_elastic_line=plot_elastic_line,
            ax=ax
        )
