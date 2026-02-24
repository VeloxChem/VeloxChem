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

    plot_absorption_spectrum(rsp_results,
                             broadening_type=broadening_type,
                             broadening_value=broadening_value,
                             ax=ax,
                             x_unit="nm")


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
