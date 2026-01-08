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

from .spectrumplot import lorentzian_absorption, gaussian_absorption


from .veloxchemlib import (hartree_in_ev,
                           hartree_in_inverse_nm,
                           fine_structure_constant,
                            avogadro_constant,
                           bohr_in_angstrom,
                           hartree_in_wavenumber,
                           )

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


# Unit conversion factors
au2ev = hartree_in_ev()
ev2au = 1.0 / au2ev
ev2nm = 1e7 / hartree_in_wavenumber()



class SpectrumAveragingDriver:
    """
    Build and plot an averaged UV/Vis spectrum from multiple LR results.

    Inputs are expected to be VeloxChem rsp results dicts (as returned by
    LinearResponseEigenSolver.compute), or (frame, rsp_results) tuples.
    """

    def __init__(self):
        # defaults chosen to match the existing plot_uv_vis style
        self.broadening_type = "lorentzian"
        self.broadening_value_ev = 0.124
        self.xstep_ev = 0.01
        self.padding_ev = 1.0

    def _unpack(self, rsp_all):
        """Accept [rsp_results, ...] or [(frame, rsp_results), ...]."""
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
        nstates=None,
        energy_min_ev=None,
        energy_max_ev=None,
    ):
        """
        Returns a dict with:
            xgrid_au, xgrid_ev, wavelength_nm (reversed for plotting)
            eps_mean, eps_std, eps_all
            transitions (energies and osc strengths per snapshot)
        """
        frames, results = self._unpack(rsp_all)

        if len(results) == 0:
            raise ValueError("No snapshots provided")

        # Collect transitions per snapshot
        all_e_au = []
        all_f = []
        for rsp in results:
            e = np.asarray(rsp["eigenvalues"], dtype=float)
            f = np.asarray(rsp["oscillator_strengths"], dtype=float)
            if nstates is not None:
                e = e[:nstates]
                f = f[:nstates]
            all_e_au.append(e)
            all_f.append(f)

        # Build a common energy grid (Hartree) for the ensemble
        if energy_min_ev is None:
            emin_au = min(np.min(e) for e in all_e_au) - self.padding_ev * ev2au
        else:
            emin_au = float(energy_min_ev) * ev2au

        if energy_max_ev is None:
            emax_au = max(np.max(e) for e in all_e_au) + self.padding_ev * ev2au
        else:
            emax_au = float(energy_max_ev) * ev2au

        emin_au = max(0.0, float(emin_au))

        xstep_au = float(self.xstep_ev) * ev2au
        xgrid_au = np.arange(emin_au, emax_au, xstep_au)

        # Broadening in a.u.
        broad_au = float(self.broadening_value_ev) * ev2au

        # Convert absorption cross section (a.u.) -> molar absorptivity epsilon (L mol^-1 cm^-1)
        c = 1.0 / fine_structure_constant()
        NA = avogadro_constant()
        a0 = bohr_in_angstrom() * 1.0e-10  # meters

        # Compute spectrum for each snapshot on the common grid
        epsilon_all = []
        bt = self.broadening_type.lower()

        for e, f in zip(all_e_au, all_f):
            if bt == "lorentzian":
                xi, yi = lorentzian_absorption(e, f, emin_au, emax_au, xstep_au, broad_au)
            elif bt == "gaussian":
                xi, yi = gaussian_absorption(e, f, emin_au, emax_au, xstep_au, broad_au)
            else:
                raise ValueError(f"Unknown broadening_type: {self.broadening_type}")

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
                raise ValueError(
                    f"All snapshots must have the same shape, got {eps.shape} vs {eps0.shape}"
                )
            eps_sum += eps
            eps_sq_sum += eps * eps

        epsilon_avg = eps_sum / nsnaps
        epsilon_std = np.sqrt(np.maximum(0.0, eps_sq_sum / nsnaps - epsilon_avg * epsilon_avg))

        # Convert grids
        xgrid_ev = xgrid_au * au2ev
        with np.errstate(divide="ignore", invalid="ignore"):
            wl_nm = ev2nm / xgrid_ev

        # Drop non-finite wavelength points (e.g. xgrid_ev == 0)
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
            "eps_all": epsilon_all,          # each is eps vs xgrid_au/xgrid_ev (energy axis)
            "eps_mean_ev": epsilon_avg,      # eps vs energy grid
            "eps_std_ev": epsilon_std,
            "eps_mean_wl": eps_mean_wl,      # eps vs wavelength grid (reversed)
            "eps_std_wl": eps_std_wl,
            "transitions_au": all_e_au,
            "oscillator_strengths": all_f,
        }

    def plot_uv_vis(
        self,
        spectrum,
        *,
        show_individual=False,
        show_sticks=True,
        show_std=False,
        title="Absorption Spectrum (Averaged)",
        ax=None,
    ):
        """
        Plot averaged spectrum in the same "wavelength + oscillator sticks" style.
        """
        if plt is None:
            raise ImportError("matplotlib is required for plotting")

        if ax is None:
            _, ax = plt.subplots()

        wl = spectrum["wavelength_nm"]
        ymean = spectrum["eps_mean_wl"]
        ystd = spectrum["eps_std_wl"]
        mask = spectrum.get("mask_wavelength_finite", None)

        if show_individual:
            xgrid_ev = spectrum["xgrid_ev"]
            with np.errstate(divide="ignore", invalid="ignore"):
                wl_forward = ev2nm / xgrid_ev

            if mask is None:
                mask = np.isfinite(wl_forward)

            wl_forward = wl_forward[mask][::-1]
            for eps in spectrum["eps_all"]:
                eps = np.asarray(eps, dtype=float)
                ax.plot(wl_forward, eps[mask][::-1], linewidth=1, alpha=0.25)

        ax.plot(wl, ymean, linewidth=2)
        ax.set_title(title)
        ax.set_xlabel("Wavelength [nm]")
        ax.set_ylabel(r"$\epsilon$ [L mol$^{-1}$ cm$^{-1}$]")

        if show_std:
            ax.fill_between(wl, ymean - ystd, ymean + ystd, alpha=0.2)

        if show_sticks:
            ax2 = ax.twinx()
            ax2.set_ylabel("Oscillator strength")

            for e_au, f in zip(spectrum["transitions_au"], spectrum["oscillator_strengths"]):
                e_ev = np.asarray(e_au, dtype=float) * au2ev
                with np.errstate(divide="ignore", invalid="ignore"):
                    wl_t = ev2nm / e_ev
                ax2.vlines(wl_t, 0.0, f, alpha=0.35, linewidth=2)

            ax2.set_ylim(0.0, max(0.2, ax2.get_ylim()[1]))

        return ax
