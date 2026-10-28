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
import sys

from .veloxchemlib import mpi_master, hartree_in_wavenumber, hartree_in_ev
from .outputstream import OutputStream
from .linearsolver import LinearSolver
from .oneeints import (compute_electric_dipole_integrals,
                       compute_linear_momentum_integrals,
                       compute_angular_momentum_integrals)
from .errorhandler import assert_msg_critical
from .spectrumplot import (plot_uv_vis_spectrum, plot_xas_spectrum,
                           plot_ecd_spectrum, plot_xcd_spectrum)

try:
    import matplotlib.pyplot as plt
except ImportError:
    pass


class TdaEigenSolverBase(LinearSolver):
    """
    Thin shared scaffold for restricted and unrestricted TDA eigensolvers.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - nstates: The number of excited states to be determined.
        - core_excitation: The flag for computing core-excited states.
        - num_core_orbitals: The number of core-orbitals involved in
          the core-excited states.
        - solver: The eigenvalues solver.
        - nto: The flag for natural transition orbital analysis.
        - nto_pairs: The number of NTO pairs in NTO analysis.
        - detach_attach: The flag for detachment/attachment density analysis.
        - cube_origin: The origin of cubic grid points.
        - cube_stepsize: The step size of cubic grid points in X, Y and Z
          directions.
        - cube_points: The number of cubic grid points in X, Y and Z directions.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes TDA excited states computation drived to default setup.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        super().__init__(comm, ostream)

        # excited states information
        self.nstates = 3

        self.initial_guess_multiplier = 3
        self.guess_scaling_threshold = 10
        self.max_subspace_dim = None
        self.collapse_nvec = None

        self.core_excitation = False
        self.num_core_orbitals = 0

        # solver setup
        self.solver = None

        # NTO and detachment/attachment density
        self.nto = False
        self.nto_pairs = None
        self.nto_cubes = False
        self.detach_attach = False
        self.detach_attach_cubes = False
        self.cube_origin = None
        self.cube_stepsize = None
        self.cube_points = [80, 80, 80]

        self._input_keywords['response'].update({
            'nstates': ('int', 'number of excited states'),
            'initial_guess_multiplier':
                ('int', 'multiplier for initial guess size'),
            'guess_scaling_threshold':
                ('int', 'threshold for guess size to increase linearly'),
            'max_subspace_dim':
                ('int', 'maximum Davidson reduced-space dimension'),
            'collapse_nvec':
                ('int', 'number of Davidson vectors kept after collapse'),
            'core_excitation': ('bool', 'compute core-excited states'),
            'num_core_orbitals': ('int', 'number of involved core-orbitals'),
            'nto': ('bool', 'analyze natural transition orbitals'),
            'nto_pairs': ('int', 'number of NTO pairs in NTO analysis'),
            'nto_cubes': ('bool', 'write NTO cube files'),
            'detach_attach': ('bool', 'analyze detachment/attachment density'),
            'detach_attach_cubes':
                ('bool', 'write detachment/attachment density cube files'),
            'cube_origin': ('seq_fixed', 'origin of cubic grid points'),
            'cube_stepsize': ('seq_fixed', 'step size of cubic grid points'),
            'cube_points': ('seq_fixed_int', 'number of cubic grid points'),
        })

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in TDA excited states computation.

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

    def _check_convergence(self):
        """
        Checks convergence of excitation energies and set convergence flag on
        all processes within MPI communicator.

        :param iteration:
            The current excited states solver iteration.
        """

        self._is_converged = False

        if self.rank == mpi_master():
            self._is_converged = self.solver.check_convergence(self.conv_thresh)

        self._is_converged = self.comm.bcast(self._is_converged,
                                             root=mpi_master())

    def _comp_onee_integrals(self, molecule, basis):
        """
        Computes one-electron integrals.

        :param molecule:
            The molecule.
        :param basis:
            The AO basis set.

        :return:
            The one-electron integrals.
        """

        dipole_mats = compute_electric_dipole_integrals(molecule, basis,
                                                        [0.0, 0.0, 0.0])
        linmom_mats = compute_linear_momentum_integrals(molecule, basis)
        angmom_mats = compute_angular_momentum_integrals(
            molecule, basis, [0.0, 0.0, 0.0])

        integrals = {}

        if self.rank == mpi_master():
            integrals['electric_dipole'] = (
                dipole_mats[0],
                dipole_mats[1],
                dipole_mats[2],
            )

            integrals['linear_momentum'] = (
                -1.0 * linmom_mats[0],
                -1.0 * linmom_mats[1],
                -1.0 * linmom_mats[2],
            )

            integrals['angular_momentum'] = (
                -1.0 * angmom_mats[0],
                -1.0 * angmom_mats[1],
                -1.0 * angmom_mats[2],
            )

            integrals['magnetic_dipole'] = (
                0.5 * angmom_mats[0],
                0.5 * angmom_mats[1],
                0.5 * angmom_mats[2],
            )

        return integrals

    def _print_iter_data(self, iteration):
        """
        Prints excited states solver iteration data to output stream.

        :param iteration:
            The current excited states solver iteration.
        """

        if self.solver.collapsed_subspace:
            exec_str = 'Collapsed reduced space: {:d}->{:d}'.format(
                self.solver.collapsed_from_dim, self.solver.collapsed_to_dim)
            self.ostream.print_info(exec_str)
            self.ostream.print_blank()

        exec_str = ' *** Iteration: ' + (str(iteration + 1)).rjust(3)
        exec_str += ' * Reduced Space: '
        exec_str += (str(self.solver.reduced_space_size())).rjust(4)
        rmax, rmin = self.solver.max_min_residual_norms()
        exec_str += ' * Residues (Max,Min): {:.2e} and {:.2e}'.format(
            rmax, rmin)
        self.ostream.print_header(exec_str)
        self.ostream.print_blank()

        reigs, rnorms = self.solver.get_eigenvalues()
        for i in range(reigs.shape[0]):
            exec_str = 'State {:2d}: {:5.8f} '.format(i + 1, reigs[i])
            exec_str += 'a.u. Residual Norm: {:3.8f}'.format(rnorms[i])
            self.ostream.print_header(exec_str.ljust(84))

        self.ostream.print_blank()
        self.ostream.flush()

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
            setattr(new_rsp_drv, key, deepcopy(val))

        return new_rsp_drv

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
