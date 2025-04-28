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
from io import StringIO
from contextlib import redirect_stderr
import numpy as np
import sys

from .veloxchemlib import mpi_master
from .veloxchemlib import boltzmann_in_hartreeperkelvin, hartree_in_kjpermol
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical

with redirect_stderr(StringIO()) as fg_err:
    try:
        import pymbar
    except ImportError:
        pass
    try:
        import scipy
    except ImportError:
        pass


class EvbDataProcessing:

    def __init__(self, comm=None, ostream=None):
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        # output stream
        self.ostream = ostream

        # mpi information
        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.results: dict = {}  # Dictionary of dictionaries with all data

        self.barrier: float = None
        self.free_energy: float = None

        self.alpha: float = None
        self.H12: float = None

        self.max_iter: int = 100
        self.fitting_slowness: float = .5
        self.tol: float = 0.01

        self.alpha_guess: float = 0
        self.H12_guess: float = 10

        # Boltzmann constant: ca 8.314e-3 kJ/mol/K
        self.kb = boltzmann_in_hartreeperkelvin() * hartree_in_kjpermol()
        self.verbose: bool = True

        self.calculate_discrete = True
        self.calculate_analytical = True
        self.smooth_window_size = 10
        self.smooth_polynomial_order = 3
        self.coordinate_bins = np.array([])
        self.bin_size = 10
        self.dens_threshold = 0.1
        self.dens_max_window = 50

    def _beta(self, T) -> float:
        return 1 / (self.kb * T)

    def compute(self, results, barrier, free_energy):
        self.ostream.print_info("Starting data processing")
        self.barrier = barrier
        self.free_energy = free_energy
        self.results = results["configuration_results"]

        self.Lambda = results["Lambda"]
        self.Lambda_frame = results["Lambda_frame"]
        self.Lambda_indices = results["Lambda_indices"]

        if self.alpha is None or self.H12 is None:
            self.ostream.print_info("Fitting H12 and alpha")
            if self.alpha is not None:
                self.ostream.print_info(
                    "Overwriting provided alpha. Provide both H12 and alpha to skip fitting"
                )
            if self.H12 is not None:
                self.ostream.print_info(
                    "Overwriting provided H12. Provide both H12 and alpha to skip fitting"
                )
            self.ostream.flush()
            self.alpha, self.H12 = self._fit_EVB_parameters()
        else:
            self.ostream.print_info("Using provided H12 and alpha")
        self.ostream.print_info("Calculating FEP and EVB curves")
        self.ostream.flush()
        self._get_FEP_and_EVB()
        results["coordinate_bins"] = self.coordinate_bins
        return results

    def _fit_EVB_parameters(self):

        assert_msg_critical('scipy' in sys.modules,
                            'scipy is required for EvbDataProcessing.')

        reference_key = list(self.results.keys())[0]

        E1_ref = self.results[reference_key]["E1_pes"]
        E2_ref = self.results[reference_key]["E2_pes"]
        Temp_set = self.results[reference_key]["Temp_set"]

        def get_barrier_and_free_energy_difference(x):
            alpha, H12 = x

            E2_shifted, V, dE, Eg = self._calculate_Eg_V_dE(
                E1_ref, E2_ref, alpha, H12)
            dGfep = self._calculate_dGfep(dE, Temp_set)
            xi = np.linspace(-10000, 10000, 20000)
            dGevb_ana, shiftxi, fepxi = self._dGevb_analytical(
                dGfep, self.Lambda, H12, xi)
            dGevb_smooth, barrier, free_energy = self._get_free_energies(
                dGevb_ana, fitting=True)
            barrier_dif = self.barrier - barrier
            free_energy_dif = self.free_energy - free_energy

            return barrier_dif, free_energy_dif

        alpha, H12 = scipy.optimize.fsolve(
            get_barrier_and_free_energy_difference,
            [self.alpha_guess, self.H12_guess])
        self.ostream.print_info(f"Fitted alpha: {alpha}, H12: {H12}")
        return alpha, H12

    def _calculate_Eg_V_dE(self, E1, E2, alpha, H12):
        E2_shifted = np.copy(E2) + alpha
        V = (1 - self.Lambda_frame) * E1 + self.Lambda_frame * E2_shifted
        dE = E1 - E2_shifted
        Eg = 0.5 * (
            (E1 + E2_shifted) - np.sqrt((E1 - E2_shifted)**2 + 4 * H12**2))
        return E2_shifted, V, dE, Eg

    def _calculate_dGfep(self, dE, Temp_set):

        assert_msg_critical('pymbar' in sys.modules,
                            'pymbar is required for EvbDataProcessing.')

        de_lambda = self._bin(dE)
        dG_bar = [0.0]
        for i, l in enumerate(self.Lambda[:-1]):
            delta_lambda = self.Lambda[i + 1] - l

            forward_energy = self._beta(Temp_set) * delta_lambda * de_lambda[i]
            backward_energy = self._beta(Temp_set) * delta_lambda * de_lambda[
                i + 1]

            try:
                dF = pymbar.other_estimators.bar(forward_energy,
                                                 -backward_energy,
                                                 False)["Delta_f"]
                dg_bar = -1 / self._beta(Temp_set) * dF
            except Exception as e:
                self.ostream.print_warning(
                    f"Error {e} encountered during BAR calculation, setting dG_bar to 0 for lambda {l}"
                )
                dg_bar = 0

            dG_bar.append(dG_bar[-1] + dg_bar)

        return dG_bar

    def _bin(self, data):
        binned_data = [[] for _ in range(np.max(self.Lambda_indices) + 1)]
        for i, li in enumerate(self.Lambda_indices):
            binned_data[li].append(data[i])

        binned_data = np.array(binned_data)
        return binned_data

    @staticmethod
    def _dGevb_analytical(dGfep, Lambda, H12, xi):

        def R(de):
            return np.sqrt(de**2 + 4 * H12**2)

        def shift(xi):
            return -2 * H12**2 / R(xi)

        def arg(xi):
            return 0.5 * (1 + xi / R(xi))

        shiftxi = shift(xi)
        fepxi = np.interp(arg(xi), Lambda, dGfep)
        dGevb = shiftxi + fepxi
        return dGevb, shiftxi, fepxi

    def _dGevb_discretised(self, dGfep, Eg, V, dE, Temp_set):
        V = np.array(V)
        Eg = np.array(Eg)
        dE = np.array(dE)

        dGfep = np.array(dGfep)
        li = self.Lambda_indices
        bins = self.coordinate_bins

        N = len(
            self.Lambda)  # The amount of frames, every lambda value is a frame
        S = (
            len(bins) + 1
        )  # The amount of bins, in between every value of X is a bin, and outside the range are two bins
        hist = [[[] for x in range(S)] for x in range(N)]
        beta_set = self._beta(Temp_set)
        content = np.exp(-beta_set * (Eg - V))

        Xi = np.searchsorted(bins, dE)
        for i in range(len(li)):
            bin = Xi[i]
            Lambda_index = li[i]
            hist[Lambda_index][bin].append(content[i])

        dGcor = np.zeros((N, S))
        pnscount = np.zeros((N, S))
        pns = np.zeros((N, S))
        for n in range(N):
            for s in range(S):
                if len(hist[n][s]) > 0:
                    dGcor[n, s] = -self.kb * Temp_set * np.log(
                        np.mean(
                            hist[n][s]))  #What to do with the temperature here
                    pnscount[n, s] = len(hist[n][s])
                else:
                    dGcor[n, s] = 0

        for n in range(N):
            for s in range(S):
                pnssum = np.sum(pnscount[:, s])
                if pnssum > 0:
                    pns[n, s] = pnscount[n, s] / pnssum

        pns = pns.transpose()

        pnsfep = pns @ dGfep
        pnscor = np.sum(pns * dGcor.transpose(), axis=1)
        dGevb = pnsfep + pnscor

        return dGevb, pns, dGcor

    def _get_free_energies(self, dGevb, fitting=False):

        assert_msg_critical('scipy' in sys.modules,
                            'scipy is required for EvbDataProcessing.')

        if fitting:
            dGevb_smooth = scipy.signal.savgol_filter(
                dGevb, self.smooth_window_size, self.smooth_polynomial_order)
        else:
            dGevb_smooth = dGevb

        if fitting:
            min_arg = scipy.signal.argrelmin(dGevb_smooth)[0]
            max_arg = scipy.signal.argrelmax(dGevb_smooth)[0]
        else:
            scope = len(dGevb_smooth) // 4
            min_arg = scipy.signal.argrelmin(dGevb_smooth, order=scope)[0]
            max_arg = scipy.signal.argrelmax(dGevb_smooth, order=scope)[0]
            if len(min_arg) < 2:
                min_arg = scipy.signal.argrelmin(dGevb_smooth)[0]
            if len(max_arg) < 1:
                max_arg = scipy.signal.argrelmax(dGevb_smooth)[0]

            # if len(min_arg) < 2:
            #     min_arg = [0,-1]

        if len(min_arg) >= 2:
            Erea = dGevb_smooth[min_arg[0]]
            Epro = dGevb_smooth[min_arg[-1]]
        else:
            Erea = dGevb_smooth[0]
            Epro = dGevb_smooth[-1]

        if not fitting and len(min_arg) != 2:
            self.ostream.print_warning(
                f"Found {len(min_arg)} minima in the EVB profile instead of 2. Confirm the calculated extrema with the plot."
            )

        if len(max_arg) == 1:
            Ebar = dGevb_smooth[max_arg[0]]
        elif len(max_arg) > 1:
            mid_arg = max_arg[len(max_arg) // 2]
            Ebar = dGevb_smooth[mid_arg]
        else:
            Ebar = dGevb_smooth[len(dGevb_smooth) // 2]

        if not fitting and len(max_arg) != 1:
            self.ostream.print_warning(
                f"Found {len(max_arg)} maxima in the EVB profile instead of 1. Confirm the calculated extrema with the plot."
            )

        barrier = Ebar - Erea
        free_energy = Epro - Erea
        dGevb_smooth -= Erea
        self.ostream.flush()
        return dGevb_smooth, barrier, free_energy

    def _get_FEP_and_EVB(self):

        for result in self.results.values():
            # Temp = result["Temp_step"]
            E1_ref = result["E1_pes"]
            E2_ref = result["E2_pes"]
            E2_shifted, V, dE, Eg = self._calculate_Eg_V_dE(
                E1_ref, E2_ref, self.alpha, self.H12)

            E1_avg = np.average(self._bin(E1_ref), axis=1)
            E2_avg = np.average(self._bin(E2_ref), axis=1)
            Eg_avg = np.average(self._bin(Eg), axis=1)
            V_avg = np.average(self._bin(V), axis=1)
            dE_avg = np.average(self._bin(dE), axis=1)
            result.update({
                "dE": dE,
                "Eg": Eg,
                "V": V,
                "E1_avg": E1_avg,
                "E2_avg": E2_avg,
                "Eg_avg": Eg_avg,
                "V_avg": V_avg,
                "dE_avg": dE_avg,
            })

        if self.coordinate_bins.size == 0:
            self.coordinate_bins = self._calculate_coordinate_bins(
                self.Lambda_indices, self.results, self.bin_size,
                self.dens_threshold)
        else:
            self.ostream.print_info(
                f"Using provided coordinate bins, min: {self.coordinate_bins[0]}, max: {self.coordinate_bins[-1]}, length: {len(self.coordinate_bins)}"
            )
            self.ostream.flush()

        for result in self.results.values():
            dE = result["dE"]
            Temp_set = result["Temp_set"]

            dGfep = self._calculate_dGfep(dE, Temp_set)

            result.update({"dGfep": dGfep})
            if self.calculate_discrete:
                dGevb_discrete, pns, dGcor = self._dGevb_discretised(
                    dGfep, result["Eg"], result["V"], result["dE"],
                    result["Temp_set"])

                (
                    dGevb_discrete,
                    barrier_discretised,
                    reaction_free_energy_discretised,
                ) = self._get_free_energies(dGevb_discrete)

                result.update({
                    "discrete": {
                        "EVB": dGevb_discrete,
                        "free_energy": reaction_free_energy_discretised,
                        "barrier": barrier_discretised,
                        "pns": pns,
                        "dGcor": dGcor,
                    }
                })

            if self.calculate_analytical:
                dGevb_analytical, shift, fepxi = self._dGevb_analytical(
                    result["dGfep"],
                    self.Lambda,
                    self.H12,
                    self.coordinate_bins,
                )

                (
                    dGevb_analytical,
                    barrier_analytical,
                    reaction_free_energy_analytical,
                ) = self._get_free_energies(dGevb_analytical)

                result.update({
                    "analytical": {
                        "EVB": dGevb_analytical,
                        "shift": shift,
                        "fep": fepxi,
                        "free_energy": reaction_free_energy_analytical,
                        "barrier": barrier_analytical,
                    }
                })

    def _calculate_coordinate_bins(self, Lambda_indices, results, bin_size,
                                   dens_threshold):

        assert_msg_critical('scipy' in sys.modules,
                            'scipy is required for EvbDataProcessing.')

        dE_min = 0
        dE_max = 0
        for result in results.values():
            dens_max = np.array([])
            dE = result['dE']
            xy = np.vstack([dE, Lambda_indices])
            dens = scipy.stats.gaussian_kde(xy)(xy)
            dens = dens / np.max(dens)
            result.update({"dE_dens": dens})

            minde = np.min(dE)
            maxde = np.max(dE)
            steps = int((maxde - minde) // 2)
            dE_bins = np.linspace(np.min(dE), np.max(dE), steps)
            bin_inds = np.digitize(dE, dE_bins)
            for i, bin in enumerate(dE_bins):
                inds = np.where(bin_inds == i)[0]
                dE_hist = []
                for ind in inds:
                    dE_hist.append(dens[ind])
                if len(dE_hist) > 0:
                    dens_max = np.append(dens_max, np.max(dE_hist))
                else:
                    dens_max = np.append(dens_max, 0)

            middle = len(dens_max) // 2
            window_size = self.dens_max_window
            dens_max_smooth = np.zeros(len(dens_max) - window_size)
            for i in range(1, len(dens_max) - window_size):
                dens_max_smooth[i] = np.max(dens_max[i:i + window_size])
            dens_max = dens_max_smooth
            result.update({"dE_dens_max": dens_max})
            result.update({"dE_dens_threshold": dens_threshold})
            # dens_max = scipy.signal.savgol_filter(dens_max, 20, 3)
            min_inds = np.where(dens_max[:middle] < dens_threshold)[0]
            max_inds = np.where(dens_max[middle:] < dens_threshold)[0]

            if len(min_inds) == 0:
                min_inds = [0]
            if len(max_inds) == 0:
                max_inds = [len(dE_bins) - middle - 1]

            dE_min_new = dE_bins[min_inds[-1]].round()
            dE_max_new = dE_bins[max_inds[0] + middle].round()
            if dE_min_new < dE_min:
                dE_min = dE_min_new
            if dE_max_new > dE_max:
                dE_max = dE_max_new

        return np.arange(dE_min, dE_max, bin_size)

    @staticmethod
    def print_results(results, ostream):

        ostream.print_info(
            f"{'Discrete':<30} {'Barrier (kJ/mol)':>20} {'Free Energy (kJ/mol)':>20}"
        )
        for name, result in results["configuration_results"].items():
            if "discrete" in result.keys():
                ostream.print_info(
                    f"{name:<30} {result['discrete']['barrier']:20.2f} {result['discrete']['free_energy']:20.2f}"
                )

        ostream.print_info("\n")
        ostream.print_info(
            f"{'Analytical':<30} {'Barrier (kJ/mol)':>20} {'Free Energy (kJ/mol)':>20}"
        )
        for name, result in results["configuration_results"].items():
            if "analytical" in result.keys():
                ostream.print_info(
                    f"{name:<30} {result['analytical']['barrier']:20.2f} {result['analytical']['free_energy']:20.2f}"
                )

    @staticmethod
    def plot_dE_density(results):

        import matplotlib.pyplot as plt

        assert_msg_critical('scipy' in sys.modules,
                            'scipy is required for EvbDataProcessing.')

        result_count = len(results["configuration_results"])
        coordinate_bins = results["coordinate_bins"]
        indices = results["Lambda_indices"]
        Lambda = results["Lambda"]

        fig, ax = plt.subplots(result_count, 2, figsize=(12, 4 * result_count))

        dE_min = coordinate_bins[0]
        dE_max = coordinate_bins[-1]
        dE_bins = np.linspace(dE_min, dE_max, int(dE_max - dE_min // 2))

        for j, (name,
                result) in enumerate(results["configuration_results"].items()):
            dE = result["dE"]
            L_values = [Lambda[i] for i in indices]
            dens = result["dE_dens"]
            dens_max = result["dE_dens_max"]
            dens_thres = result["dE_dens_threshold"]

            # dens_max = scipy.signal.savgol_filter(dens_max, 20, 3)

            ax[j, 0].scatter(dE, L_values, c=dens, s=5)
            ax[j, 0].plot([dE_min, dE_min], [0, 0.3])
            ax[j, 0].plot([dE_max, dE_max], [0.7, 1])
            ax[j, 0].set_ylabel(r"$\lambda$")
            ax[j, 0].set_xlabel(r"$\Delta \mathcal{E}$ (kJ/mol)")
            ax[j, 0].set_ylim(0, 1)
            ax[j, 0].set_xlim(min(dE) * 1.1, max(dE) * 1.1)

            min_label = min(dE)
            min_label = min_label - min_label % 200 + 200
            max_label = max(dE)
            max_label = max_label - max_label % 200 + 200
            xlabels = np.arange(min_label, max_label, 200)
            ylabels = np.arange(0, 1.1, 0.2)

            min_tick = min(dE)
            min_tick = min_tick - min_tick % 100 + 100
            max_tick = max(dE)
            max_tick = max_tick - max_tick % 100 + 100

            minor_xticks = np.arange(min_tick, max_tick, 100)
            minor_yticks = np.arange(0, 1.1, 0.1)

            ax[j, 0].set_xticks(xlabels)
            ax[j, 0].set_yticks(ylabels)

            ax[j, 0].set_xticks(minor_xticks, minor=True)
            ax[j, 0].set_yticks(minor_yticks, minor=True)
            ax[j, 0].grid(True, linestyle='-', which='major')
            ax[j, 0].grid(True, linestyle=':', which='minor')

            middle = len(dens_max) // 2
            # start = np.where(dens_max[:middle] == 0)[0][-1]
            # end = np.where(dens_max[middle:] == 0)[0][0] + middle
            start = 0
            end = len(dens_max) - 1

            ax[j, 1].scatter(dE, dens, s=1)
            ax[j, 1].plot([dE_min, dE_max], [dens_thres, dens_thres])
            ax[j, 1].plot(dE_bins[start:end], dens_max[start:end])
            ax[j, 1].plot([dE_min, dE_min], [0, 1])
            ax[j, 1].plot([dE_max, dE_max], [0, 1])
            ax[j, 1].set_ylabel("Density")
            ax[j, 1].set_xlabel(r"$\Delta \mathcal{E}$ (kJ/mol)")
            ax[j, 1].tick_params(
                axis='y',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                bottom=False,  # ticks along the bottom edge are off
                top=False,  # ticks along the top edge are off
                labelbottom=False,  # labels along the bottom edge are off
            )

        return fig, ax

    @staticmethod
    def plot_results(results, plot_analytical=True, plot_discrete=False):

        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from matplotlib.lines import Line2D

        coordinate_bins = results["coordinate_bins"]
        Lambda = results["Lambda"]

        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        bin_indicators = (coordinate_bins[:-1] + coordinate_bins[1:]) / 2
        colors = mcolors.TABLEAU_COLORS

        colorkeys = list(colors.keys())
        legend_lines = []
        legend_labels = []
        if plot_analytical and plot_discrete:
            discrete_linestyle = "--"
        else:
            discrete_linestyle = "-"
        for i, (name,
                result) in enumerate(results["configuration_results"].items()):

            #Shift both averages by the same amount so that their relative differences stay the same
            ax[0].plot(Lambda, result["dGfep"], label=name)
            ax[0].set_xlim(0, 1)
            if plot_discrete:
                if "discrete" in result.keys():
                    ax[1].plot(
                        bin_indicators,
                        result["discrete"]["EVB"][1:-1],
                        label=f"{name} discretised",
                        color=colors[colorkeys[i]],
                        linestyle=discrete_linestyle,
                    )

            if plot_analytical:
                if "analytical" in result.keys():
                    ax[1].plot(
                        bin_indicators,
                        result["analytical"]["EVB"][1:],
                        label=f"{name} analytical",
                        color=colors[colorkeys[i]],
                    )

            ax[1].set_xlim(coordinate_bins[0], coordinate_bins[-1])

            legend_lines.append(Line2D([0], [0], color=colors[colorkeys[i]]))
            legend_labels.append(name)

        if plot_analytical and plot_discrete:
            EVB_legend_lines = []
            EVB_legend_labels = []
            EVB_legend_lines.append(
                Line2D([0], [0], linestyle="-", color="grey"))
            EVB_legend_labels.append("analytical")
            EVB_legend_lines.append(
                Line2D([0], [0], linestyle="--", color="grey"))
            EVB_legend_labels.append("discrete")
            ax[1].legend(EVB_legend_lines, EVB_legend_labels)

        # ax[0].legend()
        ax[0].set_xlabel(r"$\lambda$")
        ax[0].set_ylabel(r"$\Delta G_{FEP}$ (kJ/mol)")

        ax[1].set_xlabel(r"$\Delta \mathcal{E}$ (kJ/mol)")
        ax[1].set_ylabel(r"$\Delta G_{EVB}$ (kJ/mol)")
        fig.legend(
            legend_lines,
            legend_labels,
            loc=(0.22, 0.91),
            ncol=len(legend_labels),
        )
        return fig, ax

    @staticmethod
    def plot_evb_details(results):

        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors

        result_count = len(results["configuration_results"])
        fig, ax = plt.subplots(1, result_count, figsize=(5 * result_count, 5))
        colors = mcolors.TABLEAU_COLORS
        coordinate_bins = results["coordinate_bins"]

        for i, (name,
                result) in enumerate(results["configuration_results"].items()):
            dGfep = result["dGfep"]
            #discrete curves
            pns = result['discrete']['pns']
            dGcor = result['discrete']['dGcor']

            pnsfep = pns @ dGfep
            pnscor = np.sum(pns * dGcor.transpose(), axis=1)
            dGevb_disc = pnsfep + pnscor

            ax[i].plot(coordinate_bins,
                       pnsfep[:-1],
                       colors['tab:blue'],
                       linestyle="--")
            ax[i].plot(coordinate_bins,
                       pnscor[:-1],
                       colors['tab:blue'],
                       linestyle=":")
            ax[i].plot(coordinate_bins, dGevb_disc[:-1], colors['tab:blue'])

            #analytical curves
            shift = result['analytical']['shift']
            fepxi = result['analytical']['fep']
            dGevb_ana = shift + fepxi

            ax[i].plot(coordinate_bins,
                       shift,
                       colors['tab:orange'],
                       linestyle=":")
            ax[i].plot(coordinate_bins,
                       fepxi,
                       colors['tab:orange'],
                       linestyle="--")
            ax[i].plot(coordinate_bins, dGevb_ana, colors['tab:orange'])

            ax[i].set_title(name)

        fig.legend([
            r"$p_{n,s} \Delta G_{FEP}$",
            r"$p_{n,s} \Delta G_{cor}$",
            r"$\Delta G_{EVB,disc.}$",
            r"$\mu$",
            r"$\nu$",
            r"$\Delta G_{EVB,ana.}$",
        ])
