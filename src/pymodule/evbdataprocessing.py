import json
import math
import pymbar
import numpy as np
import scipy
import h5py

import sys
from mpi4py import MPI
from .veloxchemlib import mpi_master
from .outputstream import OutputStream


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

        self.results: dict = {} # Dictionary of dictionaries with all data

        self.barrier: float = None
        self.free_energy: float = None

        self.alpha: float = None
        self.H12: float = None

        self.max_iter: int = 100
        self.fitting_slowness: float = .5
        self.tol: float = 0.01

        self.alpha_guess: float = 0
        self.H12_guess: float = 10

        self.kb = 1.987204259e-3  # kcal/molK #todo use vlx
        self.joule_to_cal = 1 / 4.184
        self.verbose: bool = True

        self.calculate_discrete = True
        self.calculate_analytical = True
        self.smooth_window_size = 10
        self.smooth_polynomial_order = 3
        self.coordinate_bins = np.array([])
        self.bin_size = 10
        self.dens_threshold = 0.05

    def beta(self, T) -> float:
        return 1 / (self.kb * T)

    def compute(self, results, barrier, free_energy):
        self.ostream.print_info("Starting data processing")
        # if isinstance(target_folders, str):
        #     target_folders = [target_folders]

        self.barrier = barrier
        self.free_energy = free_energy
        self.ostream.print_info("Loading files")
        self.results = results
        self.ostream.print_info("Fitting H12 and alpha")
        self.alpha, self.H12 = self.fit_EVB_parameters()
        self.ostream.print_info("Calculating FEP and EVB curves")
        self.get_FEP_and_EVB()
        return results

    def fit_EVB_parameters(self):
        reference_key = list(self.results.keys())[0]
        Lambda = self.results[reference_key]["Lambda"]
        Lambda_frame = self.results[reference_key]["Lambda_frame"]
        E1_ref = self.results[reference_key]["E1_ref"]
        E2_ref = self.results[reference_key]["E2_ref"]
        Temp_set = self.results[reference_key]["Temp_set"]
        Lambda_indices = self.results[reference_key]["Lambda_indices"]

        def get_barrier_and_free_energy_difference(x):
            alpha, H12 = x

            E2_shifted, V, dE, Eg = self.calculate_Eg_V_dE(
                E1_ref,
                E2_ref,
                alpha,
                H12,
                Lambda_frame,
            )
            dGfep = self.calculate_dGfep(
                dE,
                Temp_set,
                Lambda,
                Lambda_indices,
            )  
            xi = np.linspace(-10000, 10000, 20000)
            dGevb_ana, shiftxi, fepxi = self.calculate_dGevb_analytical(dGfep, Lambda, H12, xi)
            dGevb_smooth, barrier, free_energy = self.calculate_free_energies(dGevb_ana, smooth=False)
            barrier_dif = self.barrier - barrier
            free_energy_dif = self.free_energy - free_energy
            
            return barrier_dif, free_energy_dif

        alpha, H12 = scipy.optimize.fsolve(get_barrier_and_free_energy_difference, [self.alpha_guess, self.H12_guess])
        self.ostream.print_info(f"Fitted alpha: {alpha}, H12: {H12}")
        return alpha, H12

    def calculate_Eg_V_dE(self, E1, E2, alpha, H12, lambda_frame):
        E2_shifted = np.copy(E2) + alpha
        V = (1 - lambda_frame) * E1 + lambda_frame * E2_shifted
        dE = E1 - E2_shifted
        Eg = 0.5 * ((E1 + E2_shifted) - np.sqrt((E1 - E2_shifted)**2 + 4 * H12**2))
        return E2_shifted, V, dE, Eg

    def calculate_dGfep(self, dE, Temp_set, Lambda, Lambda_indices):
        de_lambda = self.bin(dE, Lambda_indices)
        # dG_middle, dG_forward, dG_backward, dG_bar = [0.0], [0.0], [0.0], [0.0]
        dG_bar = [0.0]
        for i, l in enumerate(Lambda[:-1]):
            delta_lambda = Lambda[i + 1] - l

            forward_energy = self.beta(Temp_set) * delta_lambda * de_lambda[i]
            backward_energy = self.beta(Temp_set) * delta_lambda * de_lambda[i + 1]
            
            try:
                dF = pymbar.other_estimators.bar(forward_energy, -backward_energy, False)["Delta_f"]
                dg_bar = -1 / self.beta(Temp_set) * dF
            except Exception as e:
                self.ostream.print_warning(f"Error {e} encountered during BAR calculation, setting dG_bar to 0 for lambda {l}")
                dg_bar = 0

            dG_bar.append(dG_bar[-1] + dg_bar)
            
        return dG_bar

    def bin(self, data, lam_i):
        binned_data = [[] for _ in range(np.max(lam_i) + 1)]
        for i, li in enumerate(lam_i):
            binned_data[li].append(data[i])
        
        binned_data = np.array(binned_data)
        return binned_data

    @staticmethod
    def calculate_dGevb_analytical(dGfep, Lambda, H12,xi):

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

    def calculate_dGevb_discretised(self, dGfep, Eg, V, dE, Temp_set, Lambda, Lambda_indices, coordinate_bins):
        V = np.array(V)
        Eg = np.array(Eg)
        dE = np.array(dE)
        
        dGfep = np.array(dGfep)
        li = Lambda_indices
        bins = coordinate_bins

        N = len(Lambda)  # The amount of frames, every lambda value is a frame
        S = (len(bins) + 1)  # The amount of bins, in between every value of X is a bin, and outside the range are two bins
        hist = [[[] for x in range(S)] for x in range(N)]
        beta_set = self.beta(Temp_set)
        content = np.exp(-beta_set*(Eg-V))

        Xi = np.searchsorted(bins,dE)
        for i in range(len(li)):
            bin = Xi[i]
            Lambda_index = li[i]
            hist[Lambda_index][bin].append(content[i])

        dGcor = np.zeros((N,S))
        pnscount = np.zeros((N,S))
        pns = np.zeros((N,S))
        for n in range(N):
            for s in range(S):
                if len(hist[n][s]) > 0:
                    dGcor[n,s] = -self.kb*Temp_set*np.log(np.mean(hist[n][s])) #What to do with the temperature here
                    pnscount[n,s] = len(hist[n][s])
                else:
                    dGcor[n,s] = 0

        for n in range(N):
            for s in range(S):
                pnssum = np.sum(pnscount[:,s])
                if pnssum > 0:
                    pns[n,s] = pnscount[n,s]/pnssum
                
        pns = pns.transpose()

        pnsfep = pns@dGfep
        pnscor = np.sum(pns*dGcor.transpose(),axis=1)
        dGevb = pnsfep+pnscor

        return dGevb, pns, dGcor

    def calculate_free_energies(self, dGevb, smooth=True):
        if smooth:
            dGevb_smooth = scipy.signal.savgol_filter(dGevb, self.smooth_window_size, self.smooth_polynomial_order)
        else:
            dGevb_smooth = dGevb

        min_arg = scipy.signal.argrelmin(dGevb_smooth)[0]
        max_arg = scipy.signal.argrelmax(dGevb_smooth)[0]
        if len(min_arg) != 2:
            Erea = dGevb_smooth[0]
            Epro = dGevb_smooth[-1]
        else:
            Erea = dGevb_smooth[min_arg[0]]
            Epro = dGevb_smooth[min_arg[1]]

        if len(max_arg) != 1:
            Ebar = dGevb_smooth[len(dGevb_smooth) // 2]
        else:
            Ebar = dGevb_smooth[max_arg[0]]
        barrier = Ebar - Erea
        free_energy = Epro - Erea
        dGevb_smooth -= Erea

        return dGevb_smooth, barrier, free_energy

    def get_FEP_and_EVB(self):

        Lambda_indices = np.array([])
        Lambda = np.array([])
        for result in self.results.values():
            Temp = result["Temp_step"]
            Lambda = result["Lambda"]
            Lambda_frame = result["Lambda_frame"]
            Lambda_indices = result["Lambda_indices"]
            E1_ref = result["E1_ref"]
            E2_ref = result["E2_ref"]
            E2_shifted, V, dE, Eg = self.calculate_Eg_V_dE(E1_ref, E2_ref, self.alpha, self.H12, Lambda_frame)

            E1_avg = np.average(self.bin(E1_ref,Lambda_indices), axis =1)
            E2_avg = np.average(self.bin(E2_ref,Lambda_indices), axis =1)
            Eg_avg = np.average(self.bin(Eg,Lambda_indices), axis =1)
            V_avg = np.average(self.bin(V,Lambda_indices), axis =1)
            dE_avg = np.average(self.bin(dE,Lambda_indices), axis =1)
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

        dE_min = 0
        dE_max = 0
        for result in self.results.values():
            dE = result['dE']
            L = result['Lambda_indices']
            xy = np.vstack([dE,L])
            dens = scipy.stats.gaussian_kde(xy)(xy)
            dens = dens / np.max(dens)
            result.update({"dE_dens": dens})

            dE_bins = np.linspace(-2000,2000,2000)
            bin_inds = np.digitize(dE,dE_bins)
            dens_max = np.array([])
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
            dens_max = scipy.signal.savgol_filter(dens_max, 5, 3) #todo maybe get rid of this?
            min_inds = np.where(dens_max[:middle] < self.dens_threshold)[0]
            max_inds = np.where(dens_max[middle:] < self.dens_threshold)[0]

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
        
        coordinate_bins = np.arange(dE_min,dE_max,self.bin_size)
        self.coordinate_bins = coordinate_bins

        for result in self.results.values():
            dE = result["dE"]
            Temp_set = result["Temp_set"]
            Lambda = result["Lambda"]
            Lambda_indices = result["Lambda_indices"]

            dGfep = self.calculate_dGfep(dE, Temp_set, Lambda, Lambda_indices) 

            result.update({"dGfep": dGfep})
            if self.calculate_discrete:
                dGevb_discrete, pns, dGcor = self.calculate_dGevb_discretised(
                    dGfep,
                    result["Eg"],
                    result["V"],
                    result["dE"],
                    result["Temp_set"],
                    result["Lambda"],
                    result["Lambda_indices"],
                    coordinate_bins,
                )

                (
                    dGevb_discrete,
                    barrier_discretised,
                    reaction_free_energy_discretised,
                ) = self.calculate_free_energies(dGevb_discrete)

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
                dGevb_analytical, shift, fepxi = self.calculate_dGevb_analytical(
                    result["dGfep"],
                    result["Lambda"],
                    self.H12,
                    coordinate_bins,
)

                (
                    dGevb_analytical,
                    barrier_analytical,
                    reaction_free_energy_analytical,
                ) = self.calculate_free_energies(dGevb_analytical)

                result.update({
                    "analytical": {
                        "EVB": dGevb_analytical,
                        "shift": shift,
                        "fep": fepxi,
                        "free_energy": reaction_free_energy_analytical,
                        "barrier": barrier_analytical,
                    }
                })

    def print_results(self):
        self.ostream.print_info(f"{'Discrete':<30}\t Barrier \t\t Free Energy")
        for name, result in self.results.items():
            if "discrete" in result.keys():
                self.ostream.print_info(f"{name:<30} \t {result['discrete']['barrier']} \t {result['discrete']['free_energy']}")

        self.ostream.print_info("\n")
        self.ostream.print_info("Analytical\t Barrier \t\t Free Energy")
        for name, result in self.results.items():
            if "analytical" in result.keys():
                self.ostream.print_info(f"{name:<30} \t {result['analytical']['barrier']} \t {result['analytical']['free_energy']}")

    @staticmethod
    def import_matplotlib():
        try:
            import matplotlib.pyplot as plt
            import matplotlib.colors as mcolors
            from matplotlib.lines import Line2D
            return plt, mcolors, Line2D
        except ImportError:
            raise ImportError("Matplotlib is not installed, plotting is not possible.")

    def plot_dE_density(self):
        plt, mcolors, Line2D = self.import_matplotlib()

        fig, ax = plt.subplots(len(self.results), 2, figsize=(12, 4*len(self.results)))
        dE_min = self.coordinate_bins[0]
        dE_max = self.coordinate_bins[-1]
        dE_bins = np.linspace(-2000,2000,2000)
        for j, (name, result) in enumerate(self.results.items()):
            dE = result["dE"]
            indices = result["Lambda_indices"]
            Lambda = result["Lambda"]
            L_values = [Lambda[i] for i in indices]
            dens = result["dE_dens"]

            dens_max = np.array([])
            bin_inds = np.digitize(dE,dE_bins)
            for i, bin in enumerate(dE_bins):
                inds = np.where(bin_inds == i)[0]
                dE_hist = []
                for ind in inds:
                    dE_hist.append(dens[ind])
                if len(dE_hist)>0:
                    dens_max = np.append(dens_max, np.max(dE_hist))
                else:
                    dens_max = np.append(dens_max, 0)
            dens_max = scipy.signal.savgol_filter(dens_max, 5, 3) #todo maybe get rid of this?

            ax[j,0].scatter(dE,L_values,c=dens,s=5)
            ax[j,0].plot([dE_min,dE_min],[0,0.3])
            ax[j,0].plot([dE_max,dE_max],[0.7,1])
            ax[j,0].set_ylabel(r"$\lambda$")
            ax[j,0].set_xlabel(r"$\Delta \mathcal{E}$ (kcal/mol)")
            ax[j,0].set_ylim(0,1)
            ax[j,0].set_xlim(min(dE) * 1.1,max(dE) * 1.1)

            min_label = min(dE)
            min_label = min_label - min_label % 200 + 200
            max_label = max(dE)
            max_label = max_label - max_label % 200 + 200
            xlabels = np.arange(min_label,max_label,200)
            ylabels = np.arange(0,1.1,0.2)

            min_tick = min(dE)
            min_tick = min_tick - min_tick % 100 + 100
            max_tick = max(dE)
            max_tick = max_tick - max_tick % 100 + 100

            minor_xticks = np.arange(min_tick,max_tick,100)
            minor_yticks = np.arange(0,1.1,0.1)

            ax[j,0].set_xticks(xlabels)
            ax[j,0].set_yticks(ylabels)

            ax[j,0].set_xticks(minor_xticks, minor=True)
            ax[j,0].set_yticks(minor_yticks, minor=True)
            ax[j,0].grid(True,linestyle='-', which='major')
            ax[j,0].grid(True,linestyle=':', which='minor')

            middle = len(dens_max) // 2
            start = np.where(dens_max[:middle] ==0)[0][-1]
            end = np.where(dens_max[middle:] ==0)[0][0] + middle

            ax[j,1].scatter(dE,dens,s=1)
            ax[j,1].plot([-start,end],[0.05,0.05])
            ax[j,1].plot(dE_bins[start:end],dens_max[start:end])
            ax[j,1].plot([dE_min,dE_min],[0,1])
            ax[j,1].plot([dE_max,dE_max],[0,1])
            ax[j,1].set_ylabel("Density")
            ax[j,1].set_xlabel(r"$\Delta \mathcal{E}$ (kcal/mol)")
            ax[j,1].tick_params(
                axis='y',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                labelbottom=False, # labels along the bottom edge are off
            )

        return fig, ax

    def plot_results(self, plot_analytical=True, plot_discrete=True):
        plt, mcolors, Line2D = self.import_matplotlib()

        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        bin_indicators = (self.coordinate_bins[:-1] + self.coordinate_bins[1:]) / 2
        colors = mcolors.TABLEAU_COLORS

        colorkeys = list(colors.keys())
        legend_lines = []
        legend_labels = []
        if plot_analytical and plot_discrete:
            discrete_linestyle = "--"
        else:
            discrete_linestyle = "-"
        for i, (name, result) in enumerate(self.results.items()):

            #Shift both averages by the same amount so that their relative differences stay the same
            ax[0].plot(result["Lambda"], result["dGfep"], label=name)
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

            legend_lines.append(Line2D([0], [0], color=colors[colorkeys[i]]))
            legend_labels.append(name)

        E_legend_lines = []
        E_legend_labels = []
        E_legend_lines.append(Line2D([0], [0], linestyle="--", color="grey"))
        E_legend_labels.append("E1")
        E_legend_lines.append(Line2D([0], [0], linestyle=":", color="grey"))
        E_legend_labels.append("E2")
        ax[0].legend(E_legend_lines, E_legend_labels)

        if plot_analytical and plot_discrete:
            EVB_legend_lines = []
            EVB_legend_labels = []
            EVB_legend_lines.append(Line2D([0], [0], linestyle="-", color="grey"))
            EVB_legend_labels.append("analytical")
            EVB_legend_lines.append(Line2D([0], [0], linestyle="--", color="grey"))
            EVB_legend_labels.append("discrete")
            ax[1].legend(EVB_legend_lines, EVB_legend_labels)

        # ax[0].legend()
        ax[0].set_xlabel(r"$\lambda$")
        ax[0].set_ylabel(r"$\Delta G_{FEP}$ (kcal/mol)")

        ax[0].set_xlabel(r"$\Delta \mathcal{E}$ (kcal/mol)")
        ax[0].set_ylabel(r"$\Delta G_{EVB}$ (kcal/mol)")
        fig.legend(
            legend_lines,
            legend_labels,
            loc=(0.22, 0.91),
            ncol=len(legend_labels),
        )
        return fig, ax

    def plot_evb_details(self):
        plt, mcolors, Line2D = self.import_matplotlib()
        fig, ax = plt.subplots(1, len(self.results), figsize=(5 * len(self.results), 5))
        colors = mcolors.TABLEAU_COLORS
        for i, (name,result) in enumerate(self.results.items()):
            dGfep = result["dGfep"]
            #discrete curves
            pns = result['discrete']['pns']
            dGcor = result['discrete']['dGcor']

            pnsfep = pns @ dGfep
            pnscor = np.sum(pns * dGcor.transpose(),axis=1)
            dGevb_disc = pnsfep + pnscor

            ax[i].plot(self.coordinate_bins, pnsfep[:-1], colors['tab:blue'], linestyle="--")
            ax[i].plot(self.coordinate_bins, pnscor[:-1], colors['tab:blue'], linestyle=":")
            ax[i].plot(self.coordinate_bins, dGevb_disc[:-1], colors['tab:blue'])

            #analytical curves
            shift = result['analytical']['shift']
            fepxi = result['analytical']['fep']
            dGevb_ana = shift + fepxi

            ax[i].plot(self.coordinate_bins, shift, colors['tab:orange'], linestyle=":")
            ax[i].plot(self.coordinate_bins, fepxi, colors['tab:orange'], linestyle="--")
            ax[i].plot(self.coordinate_bins, dGevb_ana, colors['tab:orange'])

            ax[i].set_title(name)

        fig.legend([r"$p_{n,s} \Delta G_{FEP}$", r"$p_{n,s} \Delta G_{cor}$", r"$\Delta G_{EVB,disc.}$", r"$\mu$", r"$\nu$", r"$\Delta G_{EVB,ana.}$",])
            
