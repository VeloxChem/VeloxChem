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
from copy import deepcopy
import numpy as np
import sys

from .veloxchemlib import XCFunctional, MolecularGrid
from .veloxchemlib import (mpi_master, rotatory_strength_in_cgs, hartree_in_ev,
                           hartree_in_inverse_nm, hartree_in_wavenumber,
                           fine_structure_constant,
                           extinction_coefficient_from_beta)
from .outputstream import OutputStream
from .distributedarray import DistributedArray
from .linearsolver import LinearSolver
from .errorhandler import assert_msg_critical
from .spectrumplot import (plot_uv_vis_spectrum, plot_xas_spectrum,
                           plot_ecd_spectrum, plot_xcd_spectrum)

try:
    import matplotlib.pyplot as plt
except ImportError:
    pass


class LinearResponseEigenSolverBase(LinearSolver):
    """
    Thin shared scaffold for restricted and unrestricted LR eigensolvers.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes linear response eigensolver to default setup.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        super().__init__(comm, ostream)

        self.nstates = 3

        self.initial_guess_multiplier = 3
        self.guess_scaling_threshold = 10
        self.max_subspace_dim = None
        self.collapse_nvec = None

        self.core_excitation = False
        self.num_core_orbitals = 0

        self.nto = False
        self.nto_pairs = None
        self.nto_cubes = False
        self.detach_attach = False
        self.detach_attach_cubes = False
        self.cube_origin = None
        self.cube_stepsize = None
        self.cube_points = [80, 80, 80]

        self.esa = False
        self.esa_from_state = None

        self.collapsed_subspace = False
        self.collapsed_from_dim = None
        self.collapsed_to_dim = None

        self._input_keywords['response'].update({
            'nstates': ('int', 'number of excited states'),
            'initial_guess_multiplier':
                ('int', 'multiplier for initial guess size'),
            'guess_scaling_threshold':
                ('int', 'threshold for guess size to increase linearly'),
            'max_subspace_dim':
                ('int', 'maximum reduced-space dimension before collapse'),
            'collapse_nvec':
                ('int', 'number of Ritz states kept after collapse'),
            'core_excitation': ('bool', 'compute core-excited states'),
            'num_core_orbitals': ('int', 'number of involved core-orbitals'),
            'nto': ('bool', 'analyze natural transition orbitals'),
            'nto_pairs': ('int', 'number of NTO pairs in NTO analysis'),
            'nto_cubes': ('bool', 'write NTO cube files'),
            'detach_attach': ('bool', 'analyze detachment/attachment density'),
            'detach_attach_cubes':
                ('bool', 'write detachment/attachment density cube files'),
            'esa': ('bool', 'compute excited state absorption'),
            'esa_from_state':
                ('int', 'the state to excite from (e.g. 1 for S1)'),
            'cube_origin': ('seq_fixed', 'origin of cubic grid points'),
            'cube_stepsize': ('seq_fixed', 'step size of cubic grid points'),
            'cube_points': ('seq_fixed_int', 'number of cubic grid points'),
        })

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in linear response eigensolver.

        :param rsp_dict:
            The dictionary of response input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}

        super().update_settings(rsp_dict, method_dict)

        if self.cube_origin is not None:
            assert_msg_critical(
                len(self.cube_origin) == 3,
                f'{type(self).__name__}: cube origin needs 3 numbers')

        if self.cube_stepsize is not None:
            assert_msg_critical(
                len(self.cube_stepsize) == 3,
                f'{type(self).__name__}: cube stepsize needs 3 numbers')

        if self.cube_points is not None:
            assert_msg_critical(
                len(self.cube_points) == 3,
                f'{type(self).__name__}: cube points needs 3 integers')

        if self.detach_attach_cubes:
            self.detach_attach = True

        if getattr(self, 'detach_attach_charges', False):
            self.detach_attach = True

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
        return None

    def _reduced_space_size(self):
        """
        Gets the total reduced-space dimension across both parity sectors.
        """

        return self._dist_bger.shape(1) + self._dist_bung.shape(1)

    def _get_max_subspace_dim(self):
        """
        Gets the maximum allowed reduced-space dimension before collapse.
        """

        if self.max_subspace_dim is not None:
            return self.max_subspace_dim

        return 20 * self.nstates * 2

    def _get_collapse_nvec(self):
        """
        Gets the number of Ritz states retained during collapse.
        """

        if self.collapse_nvec is not None:
            return max(self.nstates, self.collapse_nvec)

        return 2 * self.nstates

    def _should_collapse_subspace(self):
        """
        Checks whether the reduced-space basis should be collapsed.
        """

        if self.nonlinear:
            return False

        return self._reduced_space_size() > self._get_max_subspace_dim()

    def _orthonormalize_collapsed_space(self, vectors):
        """
        Removes linear dependence and orthonormalizes a collapsed basis.

        :param vectors:
            The distributed half-sized basis vectors.

        :return:
            The orthonormalized basis.
        """

        if vectors.data.size == 0:
            return vectors

        vectors = self._remove_linear_dependence_half_size(
            vectors, self.lindep_thresh)
        vectors = self._orthogonalize_gram_schmidt_half_size(vectors)
        vectors = self._normalize_half_size(vectors)

        return vectors

    def _print_iteration(self, relative_residual_norm, ws):
        """
        Prints information of the iteration.

        :param relative_residual_norm:
            Relative residual norms.
        :param ws:
            Excitation energies.
        """

        width = 92
        output_header = '*** Iteration:   {} '.format(self._cur_iter + 1)
        output_header += '* Residuals (Max,Min): '
        output_header += '{:.2e} and {:.2e}'.format(
            max(relative_residual_norm.values()),
            min(relative_residual_norm.values()))
        self.ostream.print_header(output_header.ljust(width))
        self.ostream.print_blank()

        if self.print_level > 1:
            for k, w in enumerate(ws):
                state_label = 'Excitation {}'.format(k + 1)
                rel_res = relative_residual_norm[k]
                output_iter = '{:<15s}: {:15.8f} '.format(state_label, w)
                output_iter += 'Residual Norm: {:.8f}'.format(rel_res)
                if relative_residual_norm[k] < self.conv_thresh:
                    output_iter += '   converged'
                self.ostream.print_header(output_iter.ljust(width))
            self.ostream.print_blank()

        self.ostream.flush()

    def _precond_trials(self, vectors, precond):
        """
        Applies preconditioner to distributed trial vectors.

        :param vectors:
            The set of vectors.
        :param precond:
            The preconditioner.

        :return:
            The preconditioned gerade and ungerade trial vectors.
        """

        trials_ger = []
        trials_ung = []

        for k, X in vectors.items():
            if precond is not None:
                v = self._preconditioning(precond[k], X)
            else:
                v = X
            norms_2 = 2.0 * v.squared_norm(axis=0)
            vn = np.sqrt(np.sum(norms_2))

            if vn > self.norm_thresh:
                norms = np.sqrt(norms_2)
                if norms[0] > self.norm_thresh:
                    trials_ger.append(v.data[:, 0])
                if norms[1] > self.norm_thresh:
                    trials_ung.append(v.data[:, 1])

        new_ger = np.array(trials_ger).T
        new_ung = np.array(trials_ung).T

        dist_new_ger = DistributedArray(new_ger, self.comm, distribute=False)
        dist_new_ung = DistributedArray(new_ung, self.comm, distribute=False)

        return dist_new_ger, dist_new_ung

    def _preconditioning(self, precond, v_in):
        """
        Applies preconditioner to a tuple of distributed trial vectors.

        :param precond:
            The preconditioner.
        :param v_in:
            The input trial vectors.

        :return:
            A tuple of distributed trial vectors after preconditioning.
        """

        pa = precond.data[:, 0]
        pb = precond.data[:, 1]

        v_in_rg = v_in.data[:, 0]
        v_in_ru = v_in.data[:, 1]

        v_out_rg = pa * v_in_rg + pb * v_in_ru
        v_out_ru = pb * v_in_rg + pa * v_in_ru

        v_mat = np.hstack((
            v_out_rg.reshape(-1, 1),
            v_out_ru.reshape(-1, 1),
        ))

        return DistributedArray(v_mat, self.comm, distribute=False)

    def _print_results(self, results):
        """
        Prints results to output stream.

        :param results:
            The dictionary containing response results.
        """

        self._print_transition_dipoles(
            'Electric Transition Dipole Moments (dipole length, a.u.)',
            results['electric_transition_dipoles'])

        self._print_transition_dipoles(
            'Electric Transition Dipole Moments (dipole velocity, a.u.)',
            results['velocity_transition_dipoles'])

        self._print_transition_dipoles(
            'Magnetic Transition Dipole Moments (a.u.)',
            results['magnetic_transition_dipoles'])

        self._print_absorption('One-Photon Absorption', results)
        self._print_ecd('Electronic Circular Dichroism', results)
        self._print_excitation_details('Character of excitations:', results)

        if self.esa:
            self._print_excited_state_absorption('Excited state absorption:',
                                                 results)

    def __deepcopy__(self, memo):
        """
        Implements deepcopy.

        :param memo:
            The memo dictionary for deepcopy.

        :return:
            A deepcopy of self.
        """

        new_rsp_drv = type(self)(self.comm, self.ostream)

        for key, val in vars(self).items():
            if isinstance(val, (MPI.Intracomm, OutputStream)):
                continue

            if isinstance(val, XCFunctional):
                new_rsp_drv.key = XCFunctional(val)
            elif isinstance(val, MolecularGrid):
                new_rsp_drv.key = MolecularGrid(val)
            else:
                new_rsp_drv.key = deepcopy(val)

        return new_rsp_drv

    def get_absorption_spectrum(self, rsp_results, x_data, x_unit, b_value,
                                b_unit):
        """
        Gets absorption spectrum.
        """

        assert_msg_critical(
            x_unit.lower() in ['au', 'ev', 'nm'],
            f'{type(self).__name__}.get_absorption_spectrum: ' +
            'x_data should be au, ev or nm')

        assert_msg_critical(
            b_unit.lower() in ['au', 'ev'],
            f'{type(self).__name__}.get_absorption_spectrum: ' +
            'broadening parameter should be au or ev')

        au2ev = hartree_in_ev()
        auxnm = 1.0 / hartree_in_inverse_nm()

        exc_ene_au = rsp_results['eigenvalues']
        osc_str = rsp_results['oscillator_strengths']

        spectrum = {}

        if x_unit.lower() == 'au':
            spectrum['x_label'] = 'Photon energy [a.u.]'
        elif x_unit.lower() == 'ev':
            spectrum['x_label'] = 'Photon energy [eV]'
        elif x_unit.lower() == 'nm':
            spectrum['x_label'] = 'Wavelength [nm]'

        spectrum['y_label'] = 'Absorption cross-section [a.u.]'

        if x_unit.lower() == 'au':
            x_data_au = list(x_data)
        elif x_unit.lower() == 'ev':
            x_data_au = [x / au2ev for x in x_data]
        elif x_unit.lower() == 'nm':
            x_data_au = [auxnm / x for x in x_data]

        if b_unit.lower() == 'au':
            b_au = b_value
        elif b_unit.lower() == 'ev':
            b_au = b_value / au2ev

        y_data = []

        sigma_factor = 2.0 * np.pi * fine_structure_constant()

        for x_au in x_data_au:
            y = 0.0
            for e, f in zip(exc_ene_au, osc_str):
                b_factor = b_au / ((e - x_au)**2 + b_au**2)
                y += sigma_factor * b_factor * f
            y_data.append(y)

        spectrum['x_data'] = list(x_data)
        spectrum['y_data'] = y_data

        return spectrum

    def get_ecd_spectrum(self, rsp_results, x_data, x_unit, b_value, b_unit):
        """
        Gets ECD spectrum.
        """

        assert_msg_critical(
            x_unit.lower() in ['au', 'ev', 'nm'],
            f'{type(self).__name__}.get_ecd_spectrum: ' +
            'x_data should be au, ev or nm')

        assert_msg_critical(
            b_unit.lower() in ['au', 'ev'],
            f'{type(self).__name__}.get_ecd_spectrum: ' +
            'broadening parameter should be au or ev')

        au2ev = hartree_in_ev()
        auxnm = 1.0 / hartree_in_inverse_nm()

        exc_ene_au = rsp_results['eigenvalues']
        rot_str_au = rsp_results[
            'rotatory_strengths'] / rotatory_strength_in_cgs()

        spectrum = {}

        if x_unit.lower() == 'au':
            spectrum['x_label'] = 'Photon energy [a.u.]'
        elif x_unit.lower() == 'ev':
            spectrum['x_label'] = 'Photon energy [eV]'
        elif x_unit.lower() == 'nm':
            spectrum['x_label'] = 'Wavelength [nm]'

        spectrum['y_label'] = 'Molar circular dichroism '
        spectrum['y_label'] += '[L mol$^{-1}$ cm$^{-1}$]'

        if x_unit.lower() == 'au':
            x_data_au = list(x_data)
        elif x_unit.lower() == 'ev':
            x_data_au = [x / au2ev for x in x_data]
        elif x_unit.lower() == 'nm':
            x_data_au = [auxnm / x for x in x_data]

        if b_unit.lower() == 'au':
            b_au = b_value
        elif b_unit.lower() == 'ev':
            b_au = b_value / au2ev

        y_data = []

        delta_eps_factor = extinction_coefficient_from_beta() / 3.0

        for x_au in x_data_au:
            y = 0.0
            for e, r in zip(exc_ene_au, rot_str_au):
                b_factor = b_au / ((e - x_au)**2 + b_au**2)
                y += delta_eps_factor * b_factor * e * r
            y_data.append(y)

        spectrum['x_data'] = list(x_data)
        spectrum['y_data'] = y_data

        return spectrum

    def plot_xas(self,
                 rsp_results,
                 broadening_type="lorentzian",
                 broadening_value=(1000.0 / hartree_in_wavenumber() *
                                   hartree_in_ev()),
                 ax=None):
        """
        Plot the X-ray absorption spectrum from the response calculation.

        :param rsp_results:
            The dictionary containing the linear response results.
        :param broadening_type:
            The type of broadening to use. Either 'lorentzian' or 'gaussian'.
        :param broadening_value:
            The broadening value in eV.
        :param ax:
            The matplotlib axis to plot on.
        """

        assert_msg_critical(
            not getattr(self, 'restricted_subspace', False),
            'Plotting spectrum for restricted_subspace is not implemented.')

        assert_msg_critical(self.core_excitation,
                            'Please use plot_uv_vis for valence excitation.')

        plot_xas_spectrum(rsp_results,
                          broadening_type=broadening_type,
                          broadening_value=broadening_value,
                          ax=ax)

    def plot_uv_vis(self,
                    rsp_results,
                    broadening_type="lorentzian",
                    broadening_value=(1000.0 / hartree_in_wavenumber() *
                                      hartree_in_ev()),
                    ax=None):
        """
        Plot the UV-Vis absorption spectrum from the response calculation.

        :param rsp_results:
            The dictionary containing the linear response results.
        :param broadening_type:
            The type of broadening to use. Either 'lorentzian' or 'gaussian'.
        :param broadening_value:
            The broadening value in eV.
        :param ax:
            The matplotlib axis to plot on.
        """

        assert_msg_critical(
            not getattr(self, 'restricted_subspace', False),
            'Plotting spectrum for restricted_subspace is not implemented.')

        assert_msg_critical(not self.core_excitation,
                            'Please use plot_xas for core excitation.')

        plot_uv_vis_spectrum(rsp_results,
                             broadening_type=broadening_type,
                             broadening_value=broadening_value,
                             ax=ax)

    def plot_xcd(self,
                 rsp_results,
                 broadening_type="lorentzian",
                 broadening_value=(1000.0 / hartree_in_wavenumber() *
                                   hartree_in_ev()),
                 ax=None):
        """
        Plot the X-ray CD spectrum from the response calculation.

        :param rsp_results:
            The dictionary containing linear response results.
        :param broadening_type:
            The type of broadening to use. Either 'lorentzian' or 'gaussian'.
        :param broadening_value:
            The broadening value in eV.
        :param ax:
            The matplotlib axis to plot on.
        """

        assert_msg_critical(
            not getattr(self, 'restricted_subspace', False),
            'Plotting spectrum for restricted_subspace is not implemented.')

        assert_msg_critical(self.core_excitation,
                            'Please use plot_ecd for valence excitation.')

        plot_xcd_spectrum(rsp_results,
                          broadening_type=broadening_type,
                          broadening_value=broadening_value,
                          ax=ax)

    def plot_ecd(self,
                 rsp_results,
                 broadening_type="lorentzian",
                 broadening_value=(1000.0 / hartree_in_wavenumber() *
                                   hartree_in_ev()),
                 ax=None):
        """
        Plot the CD spectrum from the response calculation.

        :param rsp_results:
            The dictionary containing linear response results.
        :param broadening_type:
            The type of broadening to use. Either 'lorentzian' or 'gaussian'.
        :param broadening_value:
            The broadening value in eV.
        :param ax:
            The matplotlib axis to plot on.
        """

        assert_msg_critical(
            not getattr(self, 'restricted_subspace', False),
            'Plotting spectrum for restricted_subspace is not implemented.')

        assert_msg_critical(not self.core_excitation,
                            'Please use plot_xcd for core excitation.')

        plot_ecd_spectrum(rsp_results,
                          broadening_type=broadening_type,
                          broadening_value=broadening_value,
                          ax=ax)

    def plot(self,
             rsp_results,
             broadening_type="lorentzian",
             broadening_value=(1000.0 / hartree_in_wavenumber() *
                               hartree_in_ev()),
             plot_type="electronic"):
        """
        Plot the absorption or ECD spectrum from the response calculation.

        :param rsp_results:
            The dictionary containing linear response results.
        :param broadening_type:
            The type of broadening to use. 'lorentzian' or 'gaussian'.
        :param broadening_value:
            The broadening value in eV.
        :param plot_type:
            The type of plot to generate. 'uv', 'xas', 'ecd', 'xcd', or
            'electronic'.
        """

        assert_msg_critical('matplotlib' in sys.modules,
                            'matplotlib is required.')

        assert_msg_critical(
            not getattr(self, 'restricted_subspace', False),
            'Plotting spectrum for restricted_subspace is not implemented.')

        if plot_type.lower() in ["uv", "uv-vis", "uv_vis"]:
            self.plot_uv_vis(rsp_results,
                             broadening_type=broadening_type,
                             broadening_value=broadening_value)
        elif plot_type.lower() == "xas":
            self.plot_xas(rsp_results,
                          broadening_type=broadening_type,
                          broadening_value=broadening_value)
        elif plot_type.lower() == "ecd":
            self.plot_ecd(rsp_results,
                          broadening_type=broadening_type,
                          broadening_value=broadening_value)
        elif plot_type.lower() == "xcd":
            self.plot_xcd(rsp_results,
                          broadening_type=broadening_type,
                          broadening_value=broadening_value)
        elif plot_type.lower() == "electronic":
            fig, axs = plt.subplots(2, 1, figsize=(8, 10))
            # Increase the height space between subplots
            fig.subplots_adjust(hspace=0.3)

            if self.core_excitation:
                self.plot_xas(rsp_results,
                              broadening_type=broadening_type,
                              broadening_value=broadening_value,
                              ax=axs[0])
                self.plot_xcd(rsp_results,
                              broadening_type=broadening_type,
                              broadening_value=broadening_value,
                              ax=axs[1])
            else:
                self.plot_uv_vis(rsp_results,
                                 broadening_type=broadening_type,
                                 broadening_value=broadening_value,
                                 ax=axs[0])
                self.plot_ecd(rsp_results,
                              broadening_type=broadening_type,
                              broadening_value=broadening_value,
                              ax=axs[1])
        else:
            assert_msg_critical(False, 'Invalid plot type')

        plt.show()
