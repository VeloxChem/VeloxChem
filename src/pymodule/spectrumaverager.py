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

import sys
from mpi4py import MPI
import numpy as np
from .veloxchemlib import mpi_master
from .outputstream import OutputStream

from .spectrumplot import lorentzian_absorption, gaussian_absorption
from .veloxchemlib import (
    hartree_in_ev,
    fine_structure_constant,
    avogadro_constant,
    bohr_in_angstrom,
    hartree_in_wavenumber,
)

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

au2ev = hartree_in_ev()
ev2au = 1.0 / au2ev
au2nm = 1.0e7 / hartree_in_wavenumber()

class SpectrumAverager:
    """
    Build and plot an averaged UV/Vis spectrum from multiple LR results.
    Inputs are expected to be (frame, rsp_results) tuples.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables:
        - broadening_type: The broadening type, 'lorentzian' or 'gaussian'.
        - broadening_value_ev: The broadening parameter in eV.
        - xstep_ev: Photon-energy grid spacing in eV.
        - padding_ev: Padding added to the global min/max excitation energies (eV)
          when constructing the common energy grid.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes spectrum averaging driver to default setup.
        """
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)
   
        self.broadening_type = "lorentzian"
        self.broadening_value_ev = 0.124
        self.xstep_ev = 0.01
        self.padding_ev = 0.8
        self.spectrum_color = "k"
        self.stick_color = "darkcyan"

    def unpack_snapshots(self, rsp_all):
        """
        Unpack snapshots input.

        :param rsp_all:
            A list of response results dictionaries, or list of
            (frame, rsp_results) tuples.
        :return:
            frames: list with frame indices
            results: list with response results dictionaries
       """
        frames = []
        results = []
        for item in rsp_all:
            if isinstance(item, tuple) and len(item) == 2:
                frame, rsp = item
                frames.append(frame)
                results.append(rsp)
            else:
                frames.append(None)
                results.append(item)
        return frames, results

    def compute(
        self,
        rsp_all,
        *,
        energy_min_ev=None,
        energy_max_ev=None,
    ):
        """
        Compute an averaged UV/Vis spectrum from multiple response results.

        A common energy grid (in a.u.) is constrcuted from the global minimum
        and maximum excitation energies across snapshots (plus padding). Each
        snapshot is broadened and converted to molar absorptivity epsilon on the
        common grid, and the mean and standard deviation are computed.
        Followed Eqn (1) and procedure here:
        https://doi.org/10.1021/acs.jctc.5c01719

        :param rsp_all:
            A list of response results dictionaries, or a list of
            (frame, rsp_results) tuples.
        :param energy_min_ev:
            Minimum photon energy in eV for the common grid. If None, it is set to
            (global min excitation energy - padding).
        :param energy_max_ev:
            Maximum photon energy in eV for the common grid. If None, it is set to
            (global max excitation energy + padding).
        :return:
            A dictionary containing the averaged spectrum and auxilary data:
              - frames: frame indices (or None)
              - xgrid_au: energy grid in a.u.
              - xgrid_ev: energy grid in eV
              - wavelength_nm: wavelength axis in nm (reversed for plotting)
              - mask_wavelength_finite: mask applied to drop a non-finite wavelength points
              - eps_all: list of epsilon arrays for each snapshot (energy axis)
              - eps_mean_ev / eps_std_ev: mean/std epsilon vs energy grid
              - eps_mean_wl / eps_std_wl: mean/std epsilon vs wavelength grid (reversed)
              - transitions_au: excitation energies per snapshot (a.u.)
              - oscillator_strengths: oscillator strengths per snapshot
        """
        frames, results = self.unpack_snapshots(rsp_all)

        if len(results) == 0:
            raise ValueError("No snapshots provided")

        # Collect transitions per snapshot
        all_e_au = []
        all_f = []
        for rsp in results:
            e = np.asarray(rsp["eigenvalues"], dtype=float)
            f = np.asarray(rsp["oscillator_strengths"], dtype=float)
            all_e_au.append(e)
            all_f.append(f)

        # Build xmin/xmax (a.u.) for the ensemble
        if energy_min_ev is None:
            emin_au = min(np.min(e) for e in all_e_au) - float(self.padding_ev) * ev2au
        else:
            emin_au = float(energy_min_ev) * ev2au

        if energy_max_ev is None:
            emax_au = max(np.max(e) for e in all_e_au) + float(self.padding_ev) * ev2au
        else:
            emax_au = float(energy_max_ev) * ev2au

        emin_au = max(0.0, float(emin_au))

        xstep_au = float(self.xstep_ev) * ev2au
        broad_au = float(self.broadening_value_ev) * ev2au

        c = 1.0 / fine_structure_constant()
        NA = avogadro_constant()
        a0 = bohr_in_angstrom() * 1.0e-10

        epsilon_all = []
        bt = self.broadening_type.lower()

        xgrid_au = None

        for e, f in zip(all_e_au, all_f):
            if bt == "lorentzian":
                xi, yi = lorentzian_absorption(e, f, emin_au, emax_au, xstep_au, broad_au)
            elif bt == "gaussian":
                xi, yi = gaussian_absorption(e, f, emin_au, emax_au, xstep_au, broad_au)
            else:
                raise ValueError(f"Unknown broadening_type: {self.broadening_type}")

            if xgrid_au is None:
                xgrid_au = xi
            else:
                if xi.shape != xgrid_au.shape or not np.allclose(xi, xgrid_au):
                    raise ValueError("Internal error: snapshot spectrum grid differs from common grid")

            sigma_au = (2.0 * np.pi * np.pi * xi * yi) / c
            sigma_m2 = sigma_au * (a0 ** 2)
            sigma_cm2 = sigma_m2 * 1.0e4
            epsilon = sigma_cm2 * NA / (np.log(10.0) * 1.0e3)  # L mol^-1 cm^-1

            epsilon_all.append(epsilon)

        nsnaps = len(epsilon_all)
        if nsnaps == 0:
            raise ValueError("No snapshots to average")

        eps0 = np.asarray(epsilon_all[0], dtype=float)
        eps_sum = np.zeros_like(eps0)
        eps_sq_sum = np.zeros_like(eps0)

        for eps in epsilon_all:
            eps = np.asarray(eps, dtype=float)
            if eps.shape != eps0.shape:
                raise ValueError(f"All snapshots must have the same shape, got {eps.shape} vs {eps0.shape}")
            eps_sum += eps
            eps_sq_sum += eps * eps

        epsilon_avg = eps_sum / nsnaps
        epsilon_std = np.sqrt(np.maximum(0.0, eps_sq_sum / nsnaps - epsilon_avg * epsilon_avg))

        # Convert grids
        xgrid_ev = xgrid_au * au2ev

        # Wavelength directly from a.u. grid: wl_nm = (nm * a.u.) / a.u.
        with np.errstate(divide="ignore", invalid="ignore"):
            wl_nm = au2nm / xgrid_au

        # Drop non-finite wavelength points (e.g. xgrid_au == 0)
        mask = np.isfinite(wl_nm)
        wl_nm_f = wl_nm[mask]
        eps_avg_f = epsilon_avg[mask]
        eps_std_f = epsilon_std[mask]

        # wavelength axis used by plot_uv_vis (reversed)
        wl_nm_f = wl_nm_f[::-1]
        eps_mean_wl = eps_avg_f[::-1]
        eps_std_wl = eps_std_f[::-1]

        return {
            "frames": frames,
            "xgrid_au": xgrid_au,
            "xgrid_ev": xgrid_ev,
            "wavelength_nm": wl_nm_f,
            "mask_wavelength_finite": mask,
            "eps_all": epsilon_all,          # each is eps vs xgrid_au (energy axis)
            "eps_mean_ev": epsilon_avg,      # eps vs energy grid
            "eps_std_ev": epsilon_std,
            "eps_mean_wl": eps_mean_wl,      # eps vs wavelength grid (reversed)
            "eps_std_wl": eps_std_wl,
            "transitions_au": all_e_au,
            "oscillator_strengths": all_f,
        }

    def plot_uv_vis_spectra(
        self,
        rsp_all,
        energy_min_ev=None,
        energy_max_ev=None,
        show_individual=False,
        show_sticks=True,
        show_std=False,
        title="Absorption Spectrum (Averaged)",
        ax=None,
        xlim_nm=None,
    ):
        """
        Plot an averaged UV/Vis spectrum from multiple response results.

        :param rsp_all:
            A list of response results dictionaries: (frame, rsp_results) tuples 
            (as returned by the compute method).
        :param energy_min_ev:
            Minimum photon energy in eV for the common grid. If None, it is set to
            (min excitation energy - padding).
        :param energy_max_ev:
            Maximum photon energy in eV for the common grid. If None, it is set to
            (max excitation energy + padding).
        :param show_individual:
            If True, plot individual broadened spectra for each snapshot.
        :param show_sticks:
            If True, plots oscillator strength sticks at transition wavelengths.
        :param show_std:
            If True, shows a shaded area corresponding to +/- one standard deviation.
        :param title:
            Title of the plot.
        :param ax:
            Matplotlib Axes object to plot on. If None, a new figure and axes are
            created.
        :param xlim_nm:
            Tuple (xmin, xmax) to set x-axis limits in nm. If None, automatic limits are used.
        :return:
            The Matplotlib Axes object containing the plot.
        :raises ImportError:
            If matplotlib is not available.
        """
        if plt is None:
            raise ImportError("matplotlib is required for plotting")

        spectrum = self.compute(
                rsp_all,
                energy_min_ev=energy_min_ev,
                energy_max_ev=energy_max_ev,
            )

        if ax is None:
            _, ax = plt.subplots()

        wl = spectrum["wavelength_nm"]
        ymean = spectrum["eps_mean_wl"]
        ystd = spectrum["eps_std_wl"]
        mask = spectrum.get("mask_wavelength_finite", None)

        if show_individual:
            xgrid_au = spectrum["xgrid_au"]
            with np.errstate(divide="ignore", invalid="ignore"):
                wl_forward = au2nm / xgrid_au  # forward (same ordering as eps arrays)

            if mask is None:
                mask = np.isfinite(wl_forward)

            wl_forward = wl_forward[mask][::-1]
            for eps in spectrum["eps_all"]:
                eps = np.asarray(eps, dtype=float)
                ax.plot(wl_forward, eps[mask][::-1], linewidth=1, alpha=0.45)

        ax.plot(wl, ymean, linewidth=2, color=self.spectrum_color)
        ax.set_title(title)
        ax.set_xlabel("Wavelength [nm]")
        ax.set_ylabel(r"$\epsilon$ [L mol$^{-1}$ cm$^{-1}$]")

        if show_std:
            ax.fill_between(wl, ymean - ystd, ymean + ystd, alpha=0.2)

        if xlim_nm is not None:
            if len(xlim_nm) != 2:
                raise ValueError("xlim_nm must be a tuple of (xmin, xmax)")
            ax.set_xlim(float(xlim_nm[0]), float(xlim_nm[1]))

        if show_sticks:
            ax2 = ax.twinx()
            ax2.set_ylabel("Oscillator strength")

            for e_au, f in zip(spectrum["transitions_au"], spectrum["oscillator_strengths"]):
                e_au = np.asarray(e_au, dtype=float)
                with np.errstate(divide="ignore", invalid="ignore"):
                    wl_t = au2nm / e_au
                ax2.vlines(wl_t, 0.0, f, alpha=0.35, linewidth=2, color=self.stick_color)

            ax2.set_ylim(0.0, max(0.2, ax2.get_ylim()[1]))
            ax2.set_xlim(ax.get_xlim())

        return ax
