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

import numpy as np
import sys

from .veloxchemlib import (rotatory_strength_in_cgs, hartree_in_ev,
                           hartree_in_inverse_nm, fine_structure_constant,
                           extinction_coefficient_from_beta, avogadro_constant,
                           bohr_in_angstrom, hartree_in_wavenumber)
from .errorhandler import assert_msg_critical

try:
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
except ImportError:
    pass


def plot_xas_spectrum(rsp_results,
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

    plot_absorption_spectrum(rsp_results,
                             broadening_type=broadening_type,
                             broadening_value=broadening_value,
                             ax=ax,
                             x_unit="ev")


def plot_uv_vis_spectrum(rsp_results,
                         broadening_type="lorentzian",
                         broadening_value=(1000.0 / hartree_in_wavenumber() *
                                           hartree_in_ev()),
                         ax=None,
                         x_unit="nm"):
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
    :param x_unit:
        The unit of x-axis. Either 'nm' or 'ev'.
    """

    plot_absorption_spectrum(rsp_results,
                             broadening_type=broadening_type,
                             broadening_value=broadening_value,
                             ax=ax,
                             x_unit=x_unit)


def plot_absorption_spectrum(
        rsp_results,
        broadening_type="lorentzian",
        broadening_value=(1000.0 / hartree_in_wavenumber() * hartree_in_ev()),
        ax=None,
        x_unit="nm"):
    """
    Plot the absoprtion spectrum from the response calculation.

    :param rsp_results:
        The dictionary containing linear response results.
    :param broadening_type:
        The type of broadening to use. Either 'lorentzian' or 'gaussian'.
    :param broadening_value:
        The broadening value in eV.
    :param ax:
        The matplotlib axis to plot on.
    :param x_unit:
        The unit of x-axis.
    """

    assert_msg_critical('matplotlib' in sys.modules, 'matplotlib is required.')

    assert_msg_critical(x_unit.lower() in ['nm', 'ev'], 'plot: Invalid x_unit')

    use_ev = (x_unit.lower() == 'ev')

    ev_x_nm = hartree_in_ev() / hartree_in_inverse_nm()
    au2ev = hartree_in_ev()
    ev2au = 1.0 / au2ev

    # initialize the plot
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    if use_ev:
        ax.set_xlabel('Photon energy [eV]')
    else:
        ax.set_xlabel('Wavelength [nm]')
    ax.set_ylabel(r'$\epsilon$ [L mol$^{-1}$ cm$^{-1}$]')

    ax.set_title("Absorption Spectrum")

    x = (rsp_results['eigenvalues'])
    y = rsp_results['oscillator_strengths']
    xmin = max(0.0, min(x) - 0.03)
    xmax = max(x) + 0.03
    xstep = 0.0001

    ax2 = ax.twinx()

    for i in np.arange(len(rsp_results['eigenvalues'])):
        if use_ev:
            x_val = rsp_results['eigenvalues'][i] * au2ev,
        else:
            x_val = ev_x_nm / (rsp_results['eigenvalues'][i] * au2ev)
        ax2.plot(
            [x_val, x_val],
            [0.0, rsp_results['oscillator_strengths'][i]],
            alpha=0.7,
            linewidth=2,
            color="darkcyan",
        )

    c = 1.0 / fine_structure_constant()
    NA = avogadro_constant()
    a_0 = bohr_in_angstrom() * 1.0e-10

    if broadening_type.lower() == "lorentzian":
        xi, yi = lorentzian_absorption(x, y, xmin, xmax, xstep,
                                       broadening_value * ev2au)

    elif broadening_type.lower() == "gaussian":
        xi, yi = gaussian_absorption(x, y, xmin, xmax, xstep,
                                     broadening_value * ev2au)

    sigma = (2 * np.pi * np.pi * xi * yi) / c
    sigma_m2 = sigma * a_0**2
    sigma_cm2 = sigma_m2 * 10**4
    epsilon = sigma_cm2 * NA / (np.log(10) * 10**3)

    if use_ev:
        x_data = xi * au2ev
    else:
        x_data = ev_x_nm / (xi * au2ev)
    ax.plot(x_data, epsilon, color="black", alpha=0.9, linewidth=2.5)

    legend_bars = mlines.Line2D([], [],
                                color='darkcyan',
                                alpha=0.7,
                                linewidth=2,
                                label='Oscillator strength')
    label_spectrum = f'{broadening_type.capitalize()} '
    label_spectrum += f'broadening ({broadening_value:.3f} eV)'
    legend_spectrum = mlines.Line2D([], [],
                                    color='black',
                                    linestyle='-',
                                    linewidth=2.5,
                                    label=label_spectrum)
    ax2.legend(handles=[legend_bars, legend_spectrum],
               frameon=False,
               borderaxespad=0.,
               loc='center left',
               bbox_to_anchor=(1.15, 0.5))
    ax2.set_ylim(0, max(abs(rsp_results['oscillator_strengths'])) * 1.1)
    ax.set_ylim(0, max(epsilon) * 1.1)
    ax.set_ylim(bottom=0)
    ax2.set_ylim(bottom=0)
    ax2.set_ylabel("Oscillator strength")

    if use_ev:
        x_lim = (xmin * au2ev, xmax * au2ev)
    else:
        x_lim = (ev_x_nm / (xmax * au2ev), ev_x_nm / (xmin * au2ev))
    ax.set_xlim(x_lim)


def plot_xcd_spectrum(rsp_results,
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

    plot_cd_spectrum(rsp_results,
                     broadening_type=broadening_type,
                     broadening_value=broadening_value,
                     ax=ax,
                     x_unit="ev",
                     cd_name="xcd")


def plot_ecd_spectrum(rsp_results,
                      broadening_type="lorentzian",
                      broadening_value=(1000.0 / hartree_in_wavenumber() *
                                        hartree_in_ev()),
                      ax=None):
    """
    Plot the ECD spectrum from the response calculation.

    :param rsp_results:
        The dictionary containing linear response results.
    :param broadening_type:
        The type of broadening to use. Either 'lorentzian' or 'gaussian'.
    :param broadening_value:
        The broadening value in eV.
    :param ax:
        The matplotlib axis to plot on.
    """

    plot_cd_spectrum(rsp_results,
                     broadening_type=broadening_type,
                     broadening_value=broadening_value,
                     ax=ax,
                     x_unit="nm",
                     cd_name="ecd")


def plot_cd_spectrum(rsp_results,
                     broadening_type="lorentzian",
                     broadening_value=(1000.0 / hartree_in_wavenumber() *
                                       hartree_in_ev()),
                     ax=None,
                     x_unit="nm",
                     cd_name="ecd"):
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

    assert_msg_critical('matplotlib' in sys.modules, 'matplotlib is required.')

    assert_msg_critical(x_unit.lower() in ['nm', 'ev'], 'plot: Invalid x_unit')

    assert_msg_critical(cd_name.lower() in ['ecd', 'xcd'],
                        'plot: Invalid cd_name')

    use_ev = (x_unit.lower() == 'ev')

    ev_x_nm = hartree_in_ev() / hartree_in_inverse_nm()
    au2ev = hartree_in_ev()

    # initialize the plot
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    if use_ev:
        ax.set_xlabel("Photon energy [eV]")
    else:
        ax.set_xlabel("Wavelength [nm]")
    ax.set_title(f"{cd_name.upper()} Spectrum")
    ax.set_ylabel(r'$\Delta \epsilon$ [L mol$^{-1}$ cm$^{-1}$]')

    ax2 = ax.twinx()
    ax2.set_ylabel('Rotatory strength [10$^{-40}$ cgs]')

    for i in np.arange(len(rsp_results["eigenvalues"])):
        if use_ev:
            x_val = rsp_results["eigenvalues"][i] * au2ev
        else:
            x_val = ev_x_nm / (rsp_results["eigenvalues"][i] * au2ev)
        ax2.plot(
            [x_val, x_val],
            [0.0, rsp_results["rotatory_strengths"][i]],
            alpha=0.7,
            linewidth=2,
            color="darkcyan",
        )
    ax2.set_ylim(-max(abs(rsp_results["rotatory_strengths"])) * 1.1,
                 max(abs(rsp_results["rotatory_strengths"])) * 1.1)

    ax.axhline(y=0,
               marker=',',
               color='k',
               linestyle='-.',
               markersize=0,
               linewidth=0.2)

    x = (rsp_results["eigenvalues"]) * au2ev
    y = rsp_results["rotatory_strengths"]
    xmin = max(0.0, min(x) - 0.8)
    xmax = max(x) + 0.8
    xstep = 0.003

    if broadening_type.lower() == "lorentzian":
        xi, yi = lorentzian_ecd(x, y, xmin, xmax, xstep, broadening_value)

    elif broadening_type.lower() == "gaussian":
        xi, yi = gaussian_ecd(x, y, xmin, xmax, xstep, broadening_value)

    # denorm_factor is roughly 22.96 * PI
    denorm_factor = (rotatory_strength_in_cgs() /
                     (extinction_coefficient_from_beta() / 3.0))
    yi = (yi * xi) / denorm_factor

    ax.set_ylim(-max(abs(yi)) * 1.1, max(abs(yi)) * 1.1)

    if use_ev:
        x_data = xi
        x_lim = (xmin, xmax)
    else:
        x_data = ev_x_nm / xi
        x_lim = (ev_x_nm / xmax, ev_x_nm / xmin)
    ax.plot(x_data, yi, color="black", alpha=0.9, linewidth=2.5)
    ax.set_xlim(x_lim)

    # include a legend for the bar and for the broadened spectrum
    legend_bars = mlines.Line2D([], [],
                                color='darkcyan',
                                alpha=0.7,
                                linewidth=2,
                                label='Rotatory strength')
    label_spectrum = f'{broadening_type.capitalize()} '
    label_spectrum += f'broadening ({broadening_value:.3f} eV)'
    legend_spectrum = mlines.Line2D([], [],
                                    color='black',
                                    linestyle='-',
                                    linewidth=2.5,
                                    label=label_spectrum)
    ax.legend(handles=[legend_bars, legend_spectrum],
              frameon=False,
              borderaxespad=0.,
              loc='center left',
              bbox_to_anchor=(1.15, 0.5))


def lorentzian_absorption(x, y, xmin, xmax, xstep, gamma):
    xi = np.arange(xmin, xmax, xstep)
    yi = np.zeros(len(xi))
    for i in range(len(xi)):
        for k in range(len(x)):
            yi[i] = yi[i] + y[k] / x[k] * gamma / ((xi[i] - x[k])**2 + gamma**2)
    yi = yi / np.pi
    return xi, yi


def gaussian_absorption(x, y, xmin, xmax, xstep, sigma):
    xi = np.arange(xmin, xmax, xstep)
    yi = np.zeros(len(xi))
    for i in range(len(xi)):
        for k in range(len(x)):
            yi[i] = yi[i] + y[k] / x[k] * np.exp(-((xi[i] - x[k])**2) /
                                                 (2 * sigma**2))
    yi = yi / (sigma * np.sqrt(2 * np.pi))
    return xi, yi


def lorentzian_ecd(x, y, xmin, xmax, xstep, gamma):
    xi = np.arange(xmin, xmax, xstep)
    yi = np.zeros(len(xi))
    for i in range(len(xi)):
        for k in range(len(x)):
            yi[i] = yi[i] + y[k] * (gamma) / ((xi[i] - x[k])**2 + (gamma)**2)
    return xi, yi


def gaussian_ecd(x, y, xmin, xmax, xstep, sigma):
    xi = np.arange(xmin, xmax, xstep)
    yi = np.zeros(len(xi))
    for i in range(len(xi)):
        for k in range(len(x)):
            yi[i] = yi[i] + y[k] * np.exp(-((xi[i] - x[k])**2) / (2 * sigma**2))
    yi = np.pi * yi / (sigma * np.sqrt(2 * np.pi))
    return xi, yi


def plot_xps_spectrum(xps_results,
                      element=None,
                      broadening_type="lorentzian",
                      broadening_value=0.5,
                      color='vlx',
                      show_atom_labels=True,
                      color_by_atom=False,
                      ax=None):
    """
    Plot the XPS (X-ray Photoelectron Spectroscopy) spectrum for a single element.

    :param xps_results:
        The dictionary containing XPS results from XPSDriver.compute().
        Format: {element: [{
            'mo_index': int,
            'atom_index': int,
            'ionization_energy_ev': float,
            'contribution': float,
            'is_delocalized': bool,
        }, ...]}
    :param element:
        Element symbol to plot (e.g., 'C', 'O', 'N', 'F', 'S').
        If None and xps_results contains only one element, that element is plotted.
        If None and xps_results contains multiple elements, an error is raised.
    :param broadening_type:
        The type of broadening to use. Either 'lorentzian' or 'gaussian'.
    :param broadening_value:
        The broadening value (FWHM) in eV.
    :param color:
        Color scheme. Either 'vlx' for VeloxChem default (darkcyan) or 'cpk' for CPK coloring.
        If 'cpk', uses element-specific colors. Default is 'vlx'.
    :param show_atom_labels:
        If True, display atom indices as labels above peaks. Default is True.
    :param color_by_atom:
        If True, color each peak according to its atom index instead of using element color.
        Default is False.
    :param ax:
        The matplotlib axis to plot on. If None, a new figure is created.

    :return:
        The matplotlib axis object.
    """

    assert_msg_critical('matplotlib' in sys.modules, 'matplotlib is required.')

    assert_msg_critical(
        broadening_type.lower() in ['lorentzian', 'gaussian'],
        f'plot_xps_spectrum: Invalid broadening_type: {broadening_type}')

    assert_msg_critical(
        color.lower() in ['cpk', 'vlx'],
        f'plot_xps_spectrum: Invalid color: {color}')

    # Determine which element to plot
    if element is None:
        if len(xps_results) == 1:
            element = list(xps_results.keys())[0]
        else:
            assert_msg_critical(
                False,
                f'plot_xps_spectrum: Must specify element when multiple elements are present. '
                f'Available elements: {list(xps_results.keys())}')
    else:
        assert_msg_critical(
            element in xps_results,
            f'plot_xps_spectrum: Element {element} not found in results. '
            f'Available elements: {list(xps_results.keys())}')

    # Get data for the specified element
    ionization_data = xps_results[element]

    # Check if there's any data
    if len(ionization_data) == 0:
        assert_msg_critical(False, f'plot_xps_spectrum: No XPS data for element {element}')

    # CPK color scheme (including Fluorine)
    cpk_colors = {
        'C': '#909090',  # Gray
        'N': '#3050F8',  # Blue
        'O': '#FF0D0D',  # Red
        'F': '#90E050',  # Light green
        'S': '#FFFF30',  # Yellow
        'P': '#FF8000',  # Orange
        'H': '#FFFFFF',  # White
    }
    vlx_color = 'darkcyan'

    # Get color for this element
    if color.lower() == 'cpk':
        elem_color = cpk_colors.get(element, vlx_color)
    else:
        elem_color = vlx_color

    energies = np.array(
        [record['ionization_energy_ev'] for record in ionization_data])
    atom_indices = np.array([record['atom_index'] for record in ionization_data])
    intensities = np.ones(len(energies))

    # Create or use provided axis
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    # Plot the spectrum for this single element
    _plot_single_element_xps(ax, element, energies, intensities,
                             broadening_type, broadening_value,
                             elem_color, atom_indices, show_atom_labels,
                             color_by_atom)

    return ax


def plot_rixs_spectrum(rixs_results,
                       photon_index=0,
                       energy_loss=True,
                       broadening_type="lorentzian",
                       broadening_value=0.24,
                       x_unit="ev",
                       x_step=0.01,
                       ax=None):
    """
    Plot the RIXS spectrum from a RIXS calculation.

    :param rixs_results:
        The dictionary containing the RIXS results.
    :param photon_index:
        The index of the incoming photon to plot the spectrum for. (0-based)
    :param energy_loss:
        If True, plot the spectrum as a function of energy loss.
        If False, plot as a function of emission energy.
    :param broadening_type:
        The type of broadening to use. Either 'lorentzian' or 'gaussian'.
    :param broadening_value:
        The FWHM in eV.
    :param x_step:
        Grid spacing in eV for the broadened RIXS spectrum.
    :param ax:
        The matplotlib axis to plot on.
    """
    assert_msg_critical('matplotlib' in sys.modules, 'matplotlib is required.')
    assert_msg_critical(x_unit.lower() in ['nm', 'ev'], 'plot: Invalid x_unit')

    use_ev = (x_unit.lower() == 'ev')

    ev_x_nm = hartree_in_ev() / hartree_in_inverse_nm()
    au2ev = hartree_in_ev()

    # incoming photon energies in eV
    incoming_ev = np.asarray(rixs_results['photon_energies']) * au2ev

    assert_msg_critical(
        photon_index >= 0 and photon_index < incoming_ev.size,
        'plot: index of incoming photon energy (photon_index) out of range.'
    )

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    if energy_loss:
        if use_ev:
            ax.set_xlabel('Energy loss [eV]')
        else:
            ax.set_xlabel('Energy loss [nm]')
    else:
        if use_ev:
            ax.set_xlabel('Emission energy [eV]')
        else:
            ax.set_xlabel('Emission energy [nm]')

    omega_ev = incoming_ev[photon_index]

    ax.set_ylabel(r'Normalized intensity [arb. units]')
    ax.set_title(rf"RIXS Spectrum @ $\hbar \omega = {omega_ev:.2f}$ eV")

    if energy_loss:
        x = rixs_results['energy_losses'][:, photon_index]
    else:
        x = rixs_results['emission_energies'][:, photon_index]

    y = rixs_results['cross_sections'][:, photon_index]

    xmin = max(0.0, min(x) - 0.03)
    xmax = max(x) + 0.03
    xstep = x_step

    x_ev = x * au2ev
    xmin_ev = xmin * au2ev
    xmax_ev = xmax * au2ev
    xstep_ev = xstep * au2ev

    ax2 = ax.twinx()
    ax2.set_ylabel(r'Cross section, $\sigma$ [a.u.]')

    for i in np.arange(len(y)):
        if energy_loss:
            energy_au = rixs_results['energy_losses'][i, photon_index]
        else:
            energy_au = rixs_results['emission_energies'][i, photon_index]

        energy_ev = energy_au * au2ev

        if use_ev:
            x_val = energy_ev
        else:
            x_val = ev_x_nm / energy_ev

        ax2.plot(
            [x_val, x_val],
            [0.0, rixs_results['cross_sections'][i, photon_index]],
            alpha=0.7,
            linewidth=2,
            color="darkcyan",
        )

    if broadening_type.lower() == "lorentzian":
        xi, yi = lorentzian_xps(
            x_ev, y, xmin_ev, xmax_ev, xstep_ev, broadening_value
        )

    elif broadening_type.lower() == "gaussian":
        xi, yi = gaussian_xps(
            x_ev, y, xmin_ev, xmax_ev, xstep_ev, broadening_value
        )

    else:
        assert_msg_critical(False, 'plot: Invalid broadening_type')

    # Normalize to maximum
    max_intensity = np.max(yi)

    assert_msg_critical(
        max_intensity > 0.0,
        'plot_rixs_spectrum: Cannot normalize RIXS map because maximum intensity is zero.'
    )

    norm_intensity = yi / max_intensity

    if use_ev:
        x_data = xi
    else:
        x_data = ev_x_nm / xi

    ax.plot(x_data, norm_intensity, color="black", alpha=0.9, linewidth=2.5)

    legend_bars = mlines.Line2D([], [],
                                color='darkcyan',
                                alpha=0.7,
                                linewidth=2,
                                label='Cross section')

    label_spectrum = f'{broadening_type.capitalize()} '
    label_spectrum += f'broadening ({broadening_value:.3f} eV)'

    legend_spectrum = mlines.Line2D([], [],
                                    color='black',
                                    linestyle='-',
                                    linewidth=2.5,
                                    label=label_spectrum)

    ax2.legend(handles=[legend_bars, legend_spectrum],
               frameon=False,
               borderaxespad=0.,
               loc='center left',
               bbox_to_anchor=(1.15, 0.5))

    ax2.set_ylim(0, max(abs(y)) * 1.1)
    ax2.set_ylim(bottom=-0.02 * max(abs(y)))

    ax.set_ylim(0, max(norm_intensity) * 1.1)
    ax.set_ylim(bottom=-0.02 * max(abs(norm_intensity)))

    if use_ev:
        x_lim = (xmin_ev, xmax_ev)
    else:
        x_lim = (ev_x_nm / xmax_ev, ev_x_nm / xmin_ev)

    if energy_loss:
        x_lim = x_lim[::-1]

    ax.set_xlim(x_lim)
    ax2.set_xlim(x_lim)

    return ax


def plot_rixs_map(rixs_results,
                  energy_loss=True,
                  broadening_type="lorentzian",
                  broadening_value=0.24,
                  x_step=0.01,
                  x_unit="ev",
                  cmap='viridis',
                  normalize='global_max',
                  ax=None):
    """
    Plot a 2D RIXS map together with the corresponding XAS.

    :param rixs_results:
        The dictionary containing the RIXS results.
    :param energy_loss:
        If True, plot the map as a function of energy loss. (RIXS)
        If False, plot the map as a function of emission energy. (res-XES)
    :param broadening_type:
        The type of broadening to use. Either 'lorentzian' or 'gaussian'.
    :param broadening_value:
        The FWHM in eV.
    :param x_step:
        Grid spacing in eV for the broadened RIXS map.
    :param x_unit:
        Only 'ev' is supported.
    :param ax:
        If None, create new axes.
        Otherwise pass a tuple (ax_map, ax_xas).

    :return:
        (ax_map, ax_xas)
    """

    assert_msg_critical('matplotlib' in sys.modules, 'matplotlib is required.')
    assert_msg_critical(x_unit.lower() == 'ev',
                        'plot_rixs_map: only x_unit="ev" is currently supported.')

    au2ev = hartree_in_ev()

    if ax is None:
        fig, (ax_map, ax_xas) = plt.subplots(
            1, 2,
            sharey=True,
            figsize=(10, 6),
            gridspec_kw={'wspace': 0.00, 'width_ratios': [4, 1.25]}
        )
    else:
        ax_map, ax_xas = ax
        fig = ax_map.figure

    incoming_photon_energies_ev = np.array(rixs_results['photon_energies']) * au2ev
    core_eigvals_ev = np.array(rixs_results['core_eigenvalues']) * au2ev
    gamma_fwhm_ev = float(rixs_results['gamma_fwhm_ev'])
    core_osc_strs = np.array(rixs_results['core_osc_strengths'])
    cross_sections = np.array(rixs_results['cross_sections'])

    if energy_loss:
        x_au = np.array(rixs_results['energy_losses'])
    else:
        x_au = np.array(rixs_results['emission_energies'])

    x_ev = x_au * au2ev

    assert_msg_critical(
        core_eigvals_ev.size == core_osc_strs.size,
        'plot_rixs_map: number of core eigenvalues does not match oscillator strengths.'
    )

    sort_idx = np.argsort(incoming_photon_energies_ev)
    incoming_photon_energies_ev = incoming_photon_energies_ev[sort_idx]
    x_ev = x_ev[:, sort_idx]
    cross_sections = cross_sections[:, sort_idx]

    xas_sort_idx = np.argsort(core_eigvals_ev)
    core_eigvals_ev = core_eigvals_ev[xas_sort_idx]
    core_osc_strs = core_osc_strs[xas_sort_idx]

    xmin_ev = np.min(x_ev)
    xmax_ev = np.max(x_ev)

    pad = 0.5
    if energy_loss:
        val_min_ev = max(0.0, xmin_ev - pad)
        val_max_ev = xmax_ev + pad
    else:
        val_min_ev = xmin_ev - pad
        val_max_ev = xmax_ev + pad

    nr_incoming_photons = incoming_photon_energies_ev.size

    xi = None
    rixs_map = None

    for i in range(nr_incoming_photons):

        if broadening_type.lower() == "lorentzian":
            xi_tmp, yi_tmp = lorentzian_xps(
                x_ev[:, i],
                cross_sections[:, i],
                val_min_ev,
                val_max_ev,
                x_step,
                broadening_value
            )

        elif broadening_type.lower() == "gaussian":
            xi_tmp, yi_tmp = gaussian_xps(
                x_ev[:, i],
                cross_sections[:, i],
                val_min_ev,
                val_max_ev,
                x_step,
                broadening_value
            )

        else:
            assert_msg_critical(False, 'plot_rixs_map: Invalid broadening_type')

        if xi is None:
            xi = xi_tmp
            rixs_map = np.zeros((len(xi), nr_incoming_photons))

        assert_msg_critical(
            len(xi_tmp) == len(xi),
            'plot_rixs_map: inconsistent broadening grids between photon energies.'
        )

        rixs_map[:, i] = yi_tmp

    # rows = incoming photon energies, cols = loss/emission
    rixs_map_T = rixs_map.T

    if normalize is not None and normalize is not False:
        normalize = normalize.lower()

        if normalize == 'global_max':
            max_intensity = np.max(rixs_map_T)

            assert_msg_critical(
                max_intensity > 0.0,
                'plot_rixs_map: Cannot normalize RIXS map because maximum intensity is zero.'
            )

            rixs_map_T /= max_intensity

        else:
            assert_msg_critical(
                False,
                'plot_rixs_map: Invalid normalize option. Use None or "global_max".'
            )

    ymin_map = incoming_photon_energies_ev.min()
    ymax_map = incoming_photon_energies_ev.max()

    limits = [xi.min(), xi.max(), ymin_map, ymax_map]

    im = ax_map.imshow(
        rixs_map_T,
        extent=limits,
        origin='lower',
        aspect='auto',
        cmap=cmap
    )

    xas_pos = ax_xas.get_position()

    cbar_pad = 0.012
    cbar_width = 0.020

    cax = fig.add_axes([
        xas_pos.x1 + cbar_pad,
        xas_pos.y0,
        cbar_width,
        xas_pos.height
    ])

    cbar = fig.colorbar(im, cax=cax)

    cbar.ax.ticklabel_format(
        axis='y',
        style='sci',
        scilimits=(0, 0),
        useMathText=True
    )

    offset_text = cbar.ax.yaxis.get_offset_text()
    offset_text.set_ha('left')

    cbar.set_label(
        r'Normalized intensity [arb. units]'
    )

    if energy_loss:
        ax_map.set_xlabel('Energy loss [eV]')
    else:
        ax_map.set_xlabel('Emission energy [eV]')

    ax_map.set_ylabel('Photon energy [eV]')
    ax_map.set_title('RIXS map')

    # reverse x-axis for energy-loss convention
    if energy_loss:
        ax_map.set_xlim(xi.max(), xi.min())
    else:
        ax_map.set_xlim(xi.min(), xi.max())

    xas_min_ev = core_eigvals_ev.min() - 1.0
    xas_max_ev = core_eigvals_ev.max() + 1.0

    if broadening_type.lower() == "lorentzian":
        xas_x, xas_y = lorentzian_xps(
            core_eigvals_ev,
            core_osc_strs,
            xas_min_ev,
            xas_max_ev,
            0.01,
            gamma_fwhm_ev
        )

    elif broadening_type.lower() == "gaussian":
        xas_x, xas_y = gaussian_xps(
            core_eigvals_ev,
            core_osc_strs,
            xas_min_ev,
            xas_max_ev,
            0.01,
            gamma_fwhm_ev
        )

    ax_xas.plot(xas_y, xas_x, lw=2.5, color='black')
    ax_xas.fill_betweenx(xas_x, 0.0, xas_y, alpha=0.2, color='black')
    ax_xas.barh(core_eigvals_ev, core_osc_strs, height=0.05, color='darkcyan', alpha=0.7)

    ax_xas.set_title('XAS')
    ax_xas.set_xlabel('Intensity')
    ax_xas.yaxis.set_visible(False)
    ax_xas.set_xlim(left=0.0)
    ax_xas.set_xticks([])

    ax_map.set_ylim(ymin_map, ymax_map)
    ax_xas.set_ylim(ymin_map, ymax_map)
    ax_map.tick_params(axis='both', which='both',
                       direction='out',
                       bottom=True, top=False,
                       left=True, right=False,
                       labelbottom=True, labeltop=False,
                       labelleft=True, labelright=False)

    ax_xas.tick_params(axis='both', which='both',
                       direction='out',
                       bottom=True, top=False,
                       left=False, right=False,
                       labelbottom=False, labeltop=False,
                       labelleft=False, labelright=False)

    return ax_map, ax_xas


def plot_tpa_transition_spectrum(rsp_results,
                                 spectrum,
                                 x_unit='ev',
                                 polarization='linear',
                                 ax=None,
                                 show_sticks=True):
    """
    Plot the TPA transition spectrum from broadened cross sections.

    :param rsp_results:
        The dictionary containing the TPA transition results.
    :param spectrum:
        The precomputed broadened spectrum dictionary from
        ``TpaTransitionDriver.get_spectrum``.
    :param x_unit:
        The unit of x-axis. Either 'au', 'ev', or 'nm'.
    :param polarization:
        Polarization channel for stick strengths. Either 'linear' or
        'circular'.
    :param ax:
        The matplotlib axis to plot on. If None, a new figure is created.
    :param show_sticks:
        If True, overlay stick/bar strengths on a secondary axis.

    :return:
        The matplotlib axis object.
    """

    assert_msg_critical('matplotlib' in sys.modules, 'matplotlib is required.')

    assert_msg_critical(
        x_unit.lower() in ['au', 'ev', 'nm'],
        'plot_tpa_transition_spectrum: x_unit should be au, ev or nm')

    assert_msg_critical(
        polarization.lower() in ['linear', 'circular'],
        'plot_tpa_transition_spectrum: polarization should be linear or circular')

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    ax.set_title('TPA Transition Spectrum')
    ax.set_xlabel(spectrum['x_label'])
    ax.set_ylabel(spectrum['y_label'])

    ax.plot(spectrum['x_data'],
            spectrum['y_data'],
            color='black',
            alpha=0.9,
            linewidth=2.5)

    if show_sticks:
        photon_energies = np.array(rsp_results['photon_energies'], dtype=float)
        if x_unit.lower() == 'ev':
            stick_x = photon_energies * hartree_in_ev()
        elif x_unit.lower() == 'nm':
            stick_x = 1.0 / (hartree_in_inverse_nm() * photon_energies)
        else:
            stick_x = photon_energies

        strengths = [
            rsp_results['tpa_strengths'][polarization.lower()][idx]
            for idx in range(len(photon_energies))
        ]

        ax2 = ax.twinx()
        if len(stick_x) > 1:
            width = 0.01 * abs(float(np.max(stick_x) - np.min(stick_x)))
        else:
            width = 0.02 * abs(float(stick_x[0])) if len(stick_x) == 1 else 0.02
        if width == 0.0:
            width = 0.02

        ax2.bar(stick_x,
                strengths,
                width=width,
                alpha=0.45,
                color='darkcyan')
        ax2.set_ylabel(f'TPA strengths ({polarization.lower()}) [a.u.]')
        ax2.set_ylim(0, max(strengths) * 1.1 if strengths else 1.0)
        ax2.set_ylim(bottom=0)

        legend_bars = mlines.Line2D([], [],
                                    color='darkcyan',
                                    alpha=0.7,
                                    linewidth=2,
                                    label=f'{polarization.capitalize()} TPA strength')
        legend_spectrum = mlines.Line2D([], [],
                                        color='black',
                                        linestyle='-',
                                        linewidth=2.5,
                                        label='Broadened TPA cross-section')
        ax2.legend(handles=[legend_bars, legend_spectrum],
                   frameon=False,
                   borderaxespad=0.,
                   loc='center left',
                   bbox_to_anchor=(1.15, 0.5))

    ax.set_ylim(0, max(spectrum['y_data']) * 1.1 if spectrum['y_data'] else 1.0)
    ax.set_ylim(bottom=0)
    if spectrum['x_data']:
        ax.set_xlim(min(spectrum['x_data']), max(spectrum['x_data']))

    return ax


def plot_tpa_spectrum(spectrum, ax=None, interpolate=True, show_points=True):
    """
    Plot a reduced/full TPA cross-section spectrum from discrete values.

    :param spectrum:
        The spectrum dictionary from ``TpaDriverBase.get_spectrum``.
    :param ax:
        The matplotlib axis to plot on. If None, a new figure is created.
    :param interpolate:
        If True, use a smooth cubic interpolation when enough points are
        available.
    :param show_points:
        If True, show the discrete computed cross sections as markers.

    :return:
        The matplotlib axis object.
    """

    assert_msg_critical('matplotlib' in sys.modules, 'matplotlib is required.')

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))

    x_data = np.array(spectrum['x_data'], dtype=float)
    y_data = np.array(spectrum['y_data'], dtype=float)

    ax.set_xlabel(spectrum['x_label'])
    ax.set_ylabel(spectrum['y_label'])
    ax.set_title('TPA Spectrum')

    if interpolate and len(x_data) >= 3:
        try:
            from scipy.interpolate import CubicSpline
        except ImportError:
            CubicSpline = None

        if CubicSpline is not None:
            sort_idx = np.argsort(x_data)
            x_sorted = x_data[sort_idx]
            y_sorted = y_data[sort_idx]
            x_dense = np.linspace(float(np.min(x_data)), float(np.max(x_data)),
                                  500)
            cs = CubicSpline(x_sorted, y_sorted)
            ax.plot(x_dense, cs(x_dense),
                    color='black', alpha=0.9, linewidth=2.5)
        else:
            ax.plot(x_data, y_data, color='black', alpha=0.9, linewidth=2.0)
    else:
        ax.plot(x_data, y_data, color='black', alpha=0.9, linewidth=2.0)

    if show_points:
        ax.plot(x_data,
                y_data,
                linestyle='none',
                marker='o',
                markersize=5,
                markerfacecolor='none',
                markeredgecolor='darkcyan')

    ax.set_ylim(bottom=0)
    if len(x_data) > 0:
        ax.set_xlim(float(np.min(x_data)), float(np.max(x_data)))

    legend_handles = []
    if show_points:
        legend_handles.append(
            mlines.Line2D([], [],
                          linestyle='none',
                          marker='o',
                          markersize=5,
                          markerfacecolor='none',
                          markeredgecolor='darkcyan',
                          label='Computed cross sections'))
    legend_handles.append(
        mlines.Line2D([], [],
                      color='black',
                      linestyle='-',
                      linewidth=2.5 if interpolate and len(x_data) >= 3 else 2.0,
                      label='Interpolated spectrum' if interpolate and len(x_data) >= 3 else 'Spectrum'))

    ax.legend(handles=legend_handles,
              frameon=False,
              borderaxespad=0.,
              loc='center left',
              bbox_to_anchor=(1.02, 0.5))

    return ax


def _plot_single_element_xps(ax, element, energies, intensities,
                             broadening_type, broadening_value, color,
                             atom_indices=None, show_atom_labels=True,
                             color_by_atom=False):
    """
    Helper function to plot XPS spectrum for a single element.

    :param ax:
        The matplotlib axis to plot on.
    :param element:
        The element symbol.
    :param energies:
        Array of binding energies.
    :param intensities:
        Array of intensities.
    :param broadening_type:
        'lorentzian' or 'gaussian'.
    :param broadening_value:
        FWHM in eV.
    :param color:
        Color for the element.
    :param atom_indices:
        Array of atom indices corresponding to each peak. If None, no atom attribution.
    :param show_atom_labels:
        If True and atom_indices provided, show atom labels above peaks.
    :param color_by_atom:
        If True and atom_indices provided, color peaks by atom instead of using single color.
    """

    ax.set_xlabel('Binding Energy [eV]')
    ax.set_ylabel('Intensity [a.u.]')
    ax.set_title(f"{element} 1s XPS Spectrum")

    # Create secondary y-axis for stick spectrum
    ax2 = ax.twinx()
    ax2.set_ylabel('Peak Intensity')

    # Generate colors for atom-based coloring
    if color_by_atom and atom_indices is not None:
        unique_atoms = np.unique(atom_indices)
        # Use a colormap for distinguishing atoms
        cmap = plt.cm.get_cmap('tab10')
        atom_colors = {atom_idx: cmap(i % 10) for i, atom_idx in enumerate(unique_atoms)}

    # Plot stick spectrum
    for i in range(len(energies)):
        if color_by_atom and atom_indices is not None:
            peak_color = atom_colors[atom_indices[i]]
        else:
            peak_color = color

        ax2.plot(
            [energies[i], energies[i]],
            [0.0, intensities[i]],
            alpha=0.7,
            linewidth=2,
            color=peak_color,
        )

        # Add atom labels above peaks (1-based indexing)
        if show_atom_labels and atom_indices is not None:
            ax2.text(energies[i], intensities[i] * 1.05,
                     f'{element}{atom_indices[i] + 1}',
                     ha='center', va='bottom', fontsize=9,
                     fontweight='bold')

    # Set up energy range for broadened spectrum
    xmin = max(0.0, min(energies) - 5.0)
    xmax = max(energies) + 5.0
    xstep = 0.01

    # Apply broadening
    if broadening_type.lower() == "lorentzian":
        xi, yi = lorentzian_xps(energies, intensities, xmin, xmax, xstep,
                                broadening_value)
    elif broadening_type.lower() == "gaussian":
        xi, yi = gaussian_xps(energies, intensities, xmin, xmax, xstep,
                              broadening_value)

    # Plot broadened spectrum
    ax.plot(xi, yi, color="black", alpha=0.9, linewidth=2.5)

    # Set plot limits
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(0, max(yi) * 1.1)
    ax2.set_ylim(0, max(intensities) * 1.1)

    # Add legend
    legend_handles = []

    if color_by_atom and atom_indices is not None:
        # Add legend entry for each atom (1-based indexing)
        for atom_idx in np.unique(atom_indices):
            legend_handles.append(
                mlines.Line2D([], [],
                              color=atom_colors[atom_idx],
                              alpha=0.7,
                              linewidth=2,
                              label=f'{element}{atom_idx + 1}')
            )
    else:
        legend_handles.append(
            mlines.Line2D([], [],
                          color=color,
                          alpha=0.7,
                          linewidth=2,
                          label=f'{element} core orbitals')
        )

    label_spectrum = f'{broadening_type.capitalize()} '
    label_spectrum += f'broadening ({broadening_value:.2f} eV FWHM)'
    legend_spectrum = mlines.Line2D([], [],
                                    color='black',
                                    linestyle='-',
                                    linewidth=2.5,
                                    label=label_spectrum)
    legend_handles.append(legend_spectrum)

    ax2.legend(handles=legend_handles,
               frameon=False,
               borderaxespad=0.,
               loc='center left',
               bbox_to_anchor=(1.15, 0.5))


def lorentzian_xps(x, y, xmin, xmax, xstep, fwhm):
    """
    Apply Lorentzian broadening to XPS data.

    :param x:
        Array of binding energies (in eV).
    :param y:
        Array of intensities.
    :param xmin:
        Minimum energy value for result.
    :param xmax:
        Maximum energy value for result.
    :param xstep:
        Energy step for result.
    :param fwhm:
        Full width at half maximum (in eV).

    :return:
        Tuple of (energy_array, intensity_array).
    """
    xi = np.arange(xmin, xmax, xstep)
    yi = np.zeros(len(xi))
    gamma = fwhm / 2.0  # Half-width at half-maximum

    for i in range(len(xi)):
        for k in range(len(x)):
            yi[i] += y[k] * gamma / ((xi[i] - x[k])**2 + gamma**2)

    # Normalize
    yi = yi / np.pi

    return xi, yi


def gaussian_xps(x, y, xmin, xmax, xstep, fwhm):
    """
    Apply Gaussian broadening to XPS data.

    :param x:
        Array of binding energies (in eV).
    :param y:
        Array of intensities.
    :param xmin:
        Minimum energy value for result.
    :param xmax:
        Maximum energy value for result.
    :param xstep:
        Energy step for result.
    :param fwhm:
        Full width at half maximum (in eV).

    :return:
        Tuple of (energy_array, intensity_array).
    """
    xi = np.arange(xmin, xmax, xstep)
    yi = np.zeros(len(xi))
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))  # Convert FWHM to sigma

    for i in range(len(xi)):
        for k in range(len(x)):
            yi[i] += y[k] * np.exp(-((xi[i] - x[k])**2) / (2 * sigma**2))

    # Normalize
    yi = yi / (sigma * np.sqrt(2 * np.pi))

    return xi, yi
