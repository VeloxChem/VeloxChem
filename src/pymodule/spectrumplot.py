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
        Format: {element: [(mo_idx, atom_idx, ionization_energy, contribution), ...]}
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

    # Extract data for the specified element
    # Handle both old format (mo_idx, ie) and new format (mo_idx, atom_idx, ie, contribution)
    if len(ionization_data[0]) == 2:
        energies = np.array([ie for _, ie in ionization_data])
        atom_indices = None
    else:
        energies = np.array([ie for _, _, ie, _ in ionization_data])
        atom_indices = np.array([atom_idx for _, atom_idx, _, _ in ionization_data])
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
