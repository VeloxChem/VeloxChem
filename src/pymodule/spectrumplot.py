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
from .cppsolverbase import ComplexResponseSolverBase

try:
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
    from matplotlib.legend_handler import HandlerTuple
    try:
        from scipy.interpolate import make_interp_spline
    except ImportError:
        make_interp_spline = None
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


def _decode_spectrum_value(value):
    """Decode bytes / singleton arrays into a normal lower-case string."""
    if isinstance(value, np.ndarray):
        if value.size == 0:
            return ''
        return _decode_spectrum_value(value.reshape(-1)[0])
    if isinstance(value, (bytes, np.bytes_)):
        text = value.decode(errors='ignore')
    else:
        text = str(value)
    return text.strip().strip('[]').strip("'").strip('"').lower()


def _detect_spectrum_backend_and_property(rsp_results, property_name=None):
    """Infer whether the response data comes from LR or CPP."""
    backend = None
    if 'eigenvalues' in rsp_results:
        backend = 'lr'
    elif 'response_functions' in rsp_results or 'frequencies' in rsp_results:
        backend = 'cpp'

    assert_msg_critical(
        backend is not None,
        'plot_trajectory_spectrum: unsupported response data.'
    )

    prop = property_name.lower() if isinstance(property_name, str) else None
    if prop is None and 'property' in rsp_results:
        prop = _decode_spectrum_value(rsp_results['property'])
    if prop is None:
        if 'rotatory_strengths' in rsp_results or 'delta-epsilon' in rsp_results:
            prop = 'ecd'
        elif ('oscillator_strengths' in rsp_results or 'sigma' in rsp_results or
              'spectrum_y_data' in rsp_results):
            prop = 'absorption'

    assert_msg_critical(
        prop in ['absorption', 'ecd'],
        "plot_trajectory_spectrum: specify property='absorption' or 'ecd'."
    )

    return backend, prop


def _build_user_x_grid(x_unit, x_range, npts=2000):
    """Create a plotting grid from a user-supplied x-range."""
    assert_msg_critical(
        isinstance(x_range, (tuple, list)) and len(x_range) == 2,
        'plot_trajectory_spectrum: x_range must be a 2-element tuple or list.'
    )

    x0 = float(x_range[0])
    x1 = float(x_range[1])
    assert_msg_critical(
        x0 > 0.0 and x1 > 0.0 and x0 != x1,
        'plot_trajectory_spectrum: x_range values must be positive and distinct.'
    )

    x_lo, x_hi = min(x0, x1), max(x0, x1)
    return np.linspace(x_lo, x_hi, npts)


def _build_lr_x_grid(frame_results, x_unit, property_name, broadening_value,
                     x_range=None, npts=2000):
    """Create a common x-grid for averaging LR spectra across frames."""
    if x_range is not None:
        return _build_user_x_grid(x_unit, x_range, npts)

    au2ev = hartree_in_ev()
    ev_x_nm = hartree_in_ev() / hartree_in_inverse_nm()

    eigvals = []
    for _, rsp_results in frame_results:
        if rsp_results is None or 'eigenvalues' not in rsp_results:
            continue
        vals = np.asarray(rsp_results['eigenvalues'], dtype=float).ravel()
        vals = vals[np.isfinite(vals) & (vals > 0.0)]
        if vals.size:
            eigvals.append(vals)

    assert_msg_critical(
        eigvals,
        'plot_trajectory_spectrum: no LR eigenvalues available to plot.'
    )

    eigvals = np.hstack(eigvals)
    pad_ev = max(
        6.0 * float(broadening_value),
        0.6 if property_name == 'ecd' else 0.3,
    )
    xmin_ev = max(1.0e-5, float(eigvals.min()) * au2ev - pad_ev)
    xmax_ev = float(eigvals.max()) * au2ev + pad_ev

    if x_unit.lower() == 'au':
        return np.linspace(xmin_ev / au2ev, xmax_ev / au2ev, npts)
    if x_unit.lower() == 'ev':
        return np.linspace(xmin_ev, xmax_ev, npts)

    nm_lo = ev_x_nm / xmax_ev
    nm_hi = ev_x_nm / xmin_ev
    return np.linspace(nm_lo, nm_hi, npts)


def _lr_absorption_curve(rsp_results, x_data, x_unit, broadening_type,
                         broadening_value):
    """Return a broadened LR absorption curve on the requested grid."""
    au2ev = hartree_in_ev()
    auxnm = 1.0 / hartree_in_inverse_nm()
    ev2au = 1.0 / au2ev

    exc_ene_au = np.asarray(rsp_results['eigenvalues'], dtype=float).ravel()
    osc_str = np.asarray(rsp_results['oscillator_strengths'], dtype=float).ravel()

    x_grid = np.asarray(x_data, dtype=float)
    if x_unit.lower() == 'au':
        x_grid_au = x_grid
    elif x_unit.lower() == 'ev':
        x_grid_au = x_grid * ev2au
    else:
        x_grid_au = auxnm / x_grid

    gamma = float(broadening_value) * ev2au

    xmin = float(np.min(x_grid_au))
    xmax = float(np.max(x_grid_au))
    xstep = float(np.min(np.diff(x_grid_au))) if x_grid_au.size > 1 else 1.0e-4

    if broadening_type.lower() == 'lorentzian':
        xi, lineshape = lorentzian_absorption(exc_ene_au, osc_str,
                                              xmin, xmax + xstep, xstep,
                                              gamma)
    else:
        xi, lineshape = gaussian_absorption(exc_ene_au, osc_str,
                                            xmin, xmax + xstep, xstep,
                                            gamma)

    if xi.size != x_grid_au.size or not np.allclose(xi, x_grid_au):
        lineshape = np.interp(x_grid_au, xi, lineshape)

    c = 1.0 / fine_structure_constant()
    n_a = avogadro_constant()
    a_0 = bohr_in_angstrom() * 1.0e-10
    sigma = (2.0 * np.pi * np.pi * x_grid_au * lineshape) / c
    sigma_m2 = sigma * a_0**2
    sigma_cm2 = sigma_m2 * 10**4
    epsilon = sigma_cm2 * n_a / (np.log(10.0) * 10**3)

    return epsilon


def _lr_ecd_curve(rsp_results, x_data, x_unit, broadening_type,
                  broadening_value):
    """Return a broadened LR ECD curve on the requested grid."""
    au2ev = hartree_in_ev()
    ev_x_nm = hartree_in_ev() / hartree_in_inverse_nm()

    exc_ene_ev = np.asarray(rsp_results['eigenvalues'], dtype=float).ravel() * au2ev
    rot_str = np.asarray(rsp_results['rotatory_strengths'], dtype=float).ravel()

    x_grid = np.asarray(x_data, dtype=float)
    if x_unit.lower() == 'au':
        x_grid_ev = x_grid * au2ev
    elif x_unit.lower() == 'ev':
        x_grid_ev = x_grid
    else:
        x_grid_ev = ev_x_nm / x_grid

    gamma = float(broadening_value)

    xmin = float(np.min(x_grid_ev))
    xmax = float(np.max(x_grid_ev))
    if x_grid_ev.size > 1:
        xdiff = np.abs(np.diff(x_grid_ev))
        xdiff = xdiff[xdiff > 0.0]
        xstep = float(np.min(xdiff)) if xdiff.size > 0 else 0.003
    else:
        xstep = 0.003

    if broadening_type.lower() == 'lorentzian':
        xi, yi = lorentzian_ecd(exc_ene_ev, rot_str,
                                xmin, xmax + xstep, xstep,
                                gamma)
    else:
        xi, yi = gaussian_ecd(exc_ene_ev, rot_str,
                              xmin, xmax + xstep, xstep,
                              gamma)

    if xi.size != x_grid_ev.size or not np.allclose(xi, x_grid_ev):
        yi = np.interp(x_grid_ev, xi, yi)

    denorm_factor = (rotatory_strength_in_cgs() /
                     (extinction_coefficient_from_beta() / 3.0))
    return (yi * x_grid_ev) / denorm_factor


def _cpp_curve(rsp_results, property_name, x_unit):
    """Return a CPP spectrum curve using the native VeloxChem CPP logic."""
    au2ev = hartree_in_ev()
    auxnm = 1.0 / hartree_in_inverse_nm()

    x_data = None
    y_data = None

    if 'response_functions' in rsp_results:
        cpp_driver = ComplexResponseSolverBase()
        cpp_driver.set_cpp_property(property_name)
        spectrum = cpp_driver.get_spectrum(rsp_results, x_unit)
        if spectrum is not None:
            x_data = np.asarray(spectrum.get('x_data', []), dtype=float).ravel()
            y_data = np.asarray(spectrum.get('y_data', []), dtype=float).ravel()

    if x_data is None or y_data is None or y_data.size == 0:
        y_data = np.asarray(rsp_results.get('spectrum_y_data', []), dtype=float).ravel()
        freqs = np.asarray(rsp_results.get('frequencies', []), dtype=float).ravel()
        freqs = freqs[np.abs(freqs) > 1.0e-12]

        if freqs.size == 0 and y_data.size:
            freqs = np.arange(1, y_data.size + 1, dtype=float)

        if freqs.size != y_data.size:
            n = min(freqs.size, y_data.size)
            freqs = freqs[:n]
            y_data = y_data[:n]

        if x_unit.lower() == 'au':
            x_data = freqs
        elif x_unit.lower() == 'ev':
            x_data = au2ev * freqs
        else:
            x_data = auxnm / freqs

    x_data = np.asarray(x_data, dtype=float)
    y_data = np.asarray(y_data, dtype=float)

    order = np.argsort(x_data)
    x_data = x_data[order]
    y_data = y_data[order]

    if property_name == 'absorption':
        n_a = avogadro_constant()
        a_0 = bohr_in_angstrom() * 1.0e-10
        sigma_to_epsilon = a_0**2 * 10**4 * n_a / (np.log(10.0) * 10**3)
        y_data = y_data * sigma_to_epsilon

    return x_data, y_data


def _build_cpp_x_grid(frame_results, property_name, x_unit, x_range=None,
                      npts=2000):
    """Create a common x-grid for averaging CPP spectra across frames."""
    x_values = []
    for _, rsp_results in frame_results:
        x_frame, _ = _cpp_curve(rsp_results, property_name, x_unit)
        if x_frame.size:
            x_values.append(np.asarray(x_frame, dtype=float))

    assert_msg_critical(
        x_values,
        'plot_trajectory_spectrum: no CPP frequencies available to plot.'
    )

    x_all = np.hstack(x_values)
    data_lo = float(np.min(x_all))
    data_hi = float(np.max(x_all))

    if x_range is not None:
        x0 = float(x_range[0])
        x1 = float(x_range[1])
        req_lo, req_hi = min(x0, x1), max(x0, x1)
        x_lo = max(data_lo, req_lo)
        x_hi = min(data_hi, req_hi)
        assert_msg_critical(
            x_hi > x_lo,
            'plot_trajectory_spectrum: requested CPP x_range lies outside the sampled frequency window.'
        )
        return np.linspace(x_lo, x_hi, npts)

    return np.linspace(data_lo, data_hi, npts)


def _interpolate_curve(x_source, y_source, x_target, smooth=False):
    """Interpolate a curve onto a common target grid."""
    x_source = np.asarray(x_source, dtype=float)
    y_source = np.asarray(y_source, dtype=float)
    x_target = np.asarray(x_target, dtype=float)

    order = np.argsort(x_source)
    x_sorted = x_source[order]
    y_sorted = y_source[order]

    unique_x, unique_idx = np.unique(x_sorted, return_index=True)
    x_sorted = unique_x
    y_sorted = y_sorted[unique_idx]

    y_interp = np.full_like(x_target, np.nan, dtype=float)
    if x_sorted.size == 0:
        return y_interp

    in_range = (x_target >= x_sorted[0]) & (x_target <= x_sorted[-1])
    if not np.any(in_range):
        return y_interp

    x_eval = x_target[in_range]
    if (smooth and x_sorted.size >= 3 and 'make_interp_spline' in globals()
            and make_interp_spline is not None):
        k = min(3, x_sorted.size - 1)
        spline = make_interp_spline(x_sorted, y_sorted, k=k)
        y_interp[in_range] = spline(x_eval)
    else:
        y_interp[in_range] = np.interp(x_eval, x_sorted, y_sorted)

    return y_interp


def plot_trajectory_spectrum(frame_results,
                             property=None,
                             x_unit='nm',
                             broadening_type='lorentzian',
                             broadening_value=(1000.0 / hartree_in_wavenumber() *
                                               hartree_in_ev()),
                             x_range=None,
                             show_frames=True,
                             show_average=True,
                             show_sticks=True,
                             frame_alpha=0.70,
                             ax=None):
    """
    Plot one spectrum per converged frame together with the average spectrum.

    Auto x-range selection:
    - LR: from the min/max excitation energies over all frames with padding
      based on the broadening value.
    - CPP: from the min/max sampled response frequencies over all frames.
    Pass x_range=(xmin, xmax) to override this behavior in the chosen x_unit.
    """
    assert_msg_critical('matplotlib' in sys.modules, 'matplotlib is required.')
    assert_msg_critical(x_unit.lower() in ['au', 'ev', 'nm'],
                        'plot_trajectory_spectrum: invalid x_unit')

    clean_results = [
        (int(frame), rsp) for frame, rsp in frame_results if rsp is not None
    ]
    assert_msg_critical(clean_results,
                        'plot_trajectory_spectrum: no converged property frames found.')

    backend, property_name = _detect_spectrum_backend_and_property(
        clean_results[0][1], property
    )

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    else:
        fig = None

    n_frames = len(clean_results)
    frame_colors = plt.cm.viridis(np.linspace(0.10, 0.90, max(n_frames, 2)))[:n_frames]

    frame_curves = []
    stick_values = []
    x_plot = None
    ax2 = None

    if backend == 'lr':
        x_plot = _build_lr_x_grid(clean_results, x_unit, property_name,
                                  broadening_value, x_range=x_range)

        if show_sticks:
            ax2 = ax.twinx()

        for i, (frame, rsp_results) in enumerate(clean_results):
            color = frame_colors[i]
            if property_name == 'ecd':
                y_curve = _lr_ecd_curve(rsp_results, x_plot, x_unit,
                                        broadening_type, broadening_value)
            else:
                y_curve = _lr_absorption_curve(rsp_results, x_plot, x_unit,
                                               broadening_type,
                                               broadening_value)
            frame_curves.append(np.asarray(y_curve, dtype=float))

            if show_frames:
                ax.plot(
                    x_plot,
                    y_curve,
                    color=color,
                    alpha=frame_alpha,
                    linewidth=1.4,
                    label=None,
                )

            if ax2 is not None:
                eig = np.asarray(rsp_results['eigenvalues'], dtype=float)
                if x_unit.lower() == 'au':
                    stick_x = eig
                elif x_unit.lower() == 'ev':
                    stick_x = eig * hartree_in_ev()
                else:
                    stick_x = ((hartree_in_ev() / hartree_in_inverse_nm()) /
                               (eig * hartree_in_ev()))

                if property_name == 'ecd':
                    stick_y = np.asarray(rsp_results['rotatory_strengths'], dtype=float)
                    ax2.set_ylabel('Rotatory strength [10$^{-40}$ cgs]')
                else:
                    stick_y = np.asarray(rsp_results['oscillator_strengths'], dtype=float)
                    ax2.set_ylabel('Oscillator strength')

                stick_values.extend(np.asarray(stick_y, dtype=float).ravel().tolist())
                for x_val, y_val in zip(stick_x, stick_y):
                    ax2.plot([x_val, x_val], [0.0, y_val], color=color,
                             alpha=min(0.90, frame_alpha + 0.10), linewidth=1.2)
    else:
        x_plot = _build_cpp_x_grid(clean_results, property_name, x_unit,
                                   x_range=x_range)
        for i, (frame, rsp_results) in enumerate(clean_results):
            color = frame_colors[i]
            x_frame, y_frame = _cpp_curve(rsp_results, property_name, x_unit)
            y_curve = _interpolate_curve(x_frame, y_frame, x_plot, smooth=True)
            frame_curves.append(np.asarray(y_curve, dtype=float))
            if show_frames:
                ax.plot(
                    x_plot,
                    y_curve,
                    color=color,
                    alpha=frame_alpha,
                    linewidth=1.4,
                    label=None,
                )

    y_stack = np.vstack(frame_curves)
    y_avg = np.nanmean(y_stack, axis=0)

    if show_average:
        ax.plot(x_plot, y_avg, color='black', alpha=0.95, linewidth=2.5,
                label=None)

    if x_unit.lower() == 'nm':
        ax.set_xlabel('Wavelength [nm]')
        ax.set_xlim(float(np.min(x_plot)), float(np.max(x_plot)))
    elif x_unit.lower() == 'ev':
        ax.set_xlabel('Photon energy [eV]')
    else:
        ax.set_xlabel('Photon energy [a.u.]')

    if property_name == 'ecd':
        ax.set_ylabel(r'$\Delta \epsilon$ [L mol$^{-1}$ cm$^{-1}$]')
        ax.set_title('Trajectory ECD Spectrum')
        ymax = np.nanmax(np.abs(y_stack)) if y_stack.size else 1.0
        ax.set_ylim(-1.1 * ymax, 1.1 * ymax)
        ax.axhline(y=0.0, color='k', linestyle='-.', linewidth=0.4)
        if ax2 is not None:
            smax = np.max(np.abs(stick_values)) if stick_values else 1.0
            ax2.set_ylim(-1.1 * smax, 1.1 * smax if smax > 0 else 1.0)
            ax2.margins(y=0.0)
    else:
        ax.set_ylabel(r'$\epsilon$ [L mol$^{-1}$ cm$^{-1}$]')
        ax.set_title('Trajectory Absorption Spectrum')
        ymax = np.nanmax(y_stack) if y_stack.size else 1.0
        ax.set_ylim(bottom=0.0, top=1.1 * ymax if ymax > 0 else 1.0)
        if ax2 is not None:
            smax = np.max(stick_values) if stick_values else 1.0
            ax2.set_ylim(bottom=0.0, top=1.1 * smax if smax > 0 else 1.0)
            ax2.margins(y=0.0)

    legend_handles = []
    legend_labels = []
    legend_handler_map = None

    if show_frames:
        viridis_handle = tuple(
            mlines.Line2D([], [], color=color, linewidth=2.0,
                          alpha=min(0.95, frame_alpha + 0.05))
            for color in plt.cm.viridis(np.linspace(0.10, 0.90, 4))
        )
        frame_label = 'Frame 1 (viridis)' if n_frames == 1 else f'Frames 1-{n_frames} (viridis)'
        legend_handles.append(viridis_handle)
        legend_labels.append(frame_label)
        legend_handler_map = {tuple: HandlerTuple(ndivide=None, pad=0.3)}

    if show_average:
        legend_handles.append(mlines.Line2D([], [], color='black', linewidth=2.5))
        legend_labels.append('Average')

    if legend_handles:
        ax.legend(handles=legend_handles,
                  labels=legend_labels,
                  handler_map=legend_handler_map,
                  frameon=False)

    if fig is not None:
        fig.tight_layout()
        plt.show()

    return ax


def _set_x_axis_label(ax, x_unit):
    """Set a consistent x-axis label for spectrum plots."""
    unit = x_unit.lower()
    if unit == 'ev':
        ax.set_xlabel('Photon energy [eV]')
    elif unit == 'au':
        ax.set_xlabel('Photon energy [a.u.]')
    else:
        ax.set_xlabel('Wavelength [nm]')


def _set_spectrum_primary_axis(ax, property_name, title):
    """Set common title and y-axis labeling for absorption / ECD plots."""
    if property_name == 'ecd':
        ax.set_ylabel(r'$\Delta \epsilon$ [L mol$^{-1}$ cm$^{-1}$]')
    else:
        ax.set_ylabel(r'$\epsilon$ [L mol$^{-1}$ cm$^{-1}$]')
    ax.set_title(title)

