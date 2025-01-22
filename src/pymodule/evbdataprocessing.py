import json
import math
import pymbar
import numpy as np
import scipy
import h5py


class EvbDataProcessing:

    def __init__(self):
        self.results: dict  # List of dictionaries with all data
        self.reference: str  # Folder name with reference data
        self.targets: list[str]  # List of folder names with target data

        self.barrier: float
        self.free_energy: float

        self.alpha: float
        self.H12: float

        self.max_iter: int = 100
        self.fitting_slowness: float = .5
        self.tol: float = 0.01

        self.alpha_guess: float = 0
        self.H12_guess: float = 10

        self.kb: float = 1.987204259e-3  # kcal/molK #todo use vlx
        self.joule_to_cal: float = 1 / 4.184
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

    def compute(self, reference_folder: str, target_folders: str | list[str], barrier, free_energy):
        print("Starting data processing")
        if isinstance(target_folders, str):
            target_folders = [target_folders]

        self.barrier = barrier
        self.free_energy = free_energy
        print("Loading files")
        self.load_files(reference_folder, target_folders)
        print("Fitting H12 and alpha")
        self.alpha, self.H12 = self.fit_EVB_parameters()
        print("Calculating FEP and EVB curves")
        self.get_FEP_and_EVB()
        print("Saving results")
        self.save_results()

    def save_results(self):
        for name, result in self.results.items():
            with h5py.File(f"{name}.h5", "w") as file:
                for key, value in result.items():
                    if isinstance(value, np.ndarray) or isinstance(value, list):
                        file.create_dataset(key, data=value)
                    elif isinstance(value, dict):
                        group = file.create_group(key)
                        for k, v in value.items():
                            group.create_dataset(k, data=v)
                    else:
                        file.attrs[key] = value

    def load_result(self, file, name):
        with h5py.File(file, "r") as f:
            result = {}
            for key, value in f.items():
                if isinstance(value, h5py.Dataset):
                    result[key] = value[()]
                elif isinstance(value, h5py.Group):
                    group = {}
                    for k, v in value.items():
                        group[k] = v[()]
                    result[key] = group
                else:
                    result[key] = value
            self.results.update({name: result})

    def load_files(self, reference_folder: str, target_folders: list[str]):

        self.reference = reference_folder.split("/")[-1]
        self.targets = []
        for target in target_folders:
            self.targets.append(target.split("/")[-1])

        folders = [reference_folder] + target_folders

        for name, folder in zip([self.reference] + self.targets, folders):

            E_file = f"{folder}/Energies.dat"
            data_file = f"{folder}/Data_combined.dat"
            options_file = f"{folder}/options.json"

            E = np.loadtxt(E_file)
            E *= self.joule_to_cal
            E1_ref, E2_ref, E1_run, E2_run, E_m = E.T

            Data = np.loadtxt(data_file)
            step, Ep, Ek, Temp, Vol, Dens, Lambda_frame = Data.T

            with open(options_file, "r") as file:
                options = json.load(file)
            Lambda = options["Lambda"]
            Temp_set = options["temperature"]
            options.pop("Lambda")
            result = {
                "E1_ref": E1_ref,
                "E2_ref": E2_ref,
                "E1_run": E1_run,
                "E2_run": E2_run,
                "E_m": E_m,
                "step": step,
                "Ep": Ep,
                "Ek": Ek,
                "Temp_step": Temp,
                "Temp_set": Temp_set,
                "Vol": Vol,
                "Dens": Dens,
                "Lambda": Lambda,
                "Lambda_frame": Lambda_frame,
                "Lambda_indices": [np.where(np.round(Lambda, 3) == L)[0][0] for L in Lambda_frame],
                "options": options,
            }
            self.results.update({name: result})

    def fit_EVB_parameters(self):

        Lambda = self.results[self.reference]["Lambda"]
        Lambda_frame = self.results[self.reference]["Lambda_frame"]
        E1_ref = self.results[self.reference]["E1_ref"]
        E2_ref = self.results[self.reference]["E2_ref"]
        Temp = self.results[self.reference]["Temp_step"]
        Lambda_indices = self.results[self.reference]["Lambda_indices"]

        def get_barrier_and_free_energy_difference(x):
            alpha, H12 = x

            E2_shifted, V, dE, Eg = self.calculate_Eg_V_dE(
                E1_ref,
                E2_ref,
                alpha,
                H12,
                Lambda_frame,
            )
            dGfep, _, _, _ = self.calculate_dGfep(
                dE,
                Temp,
                Lambda,
                Lambda_indices,
            )  
            xi = np.linspace(-10000, 10000, 20000)
            dGevb_ana,_,_ = self.calculate_dGevb_analytical(dGfep, Lambda, H12, xi)
            _, barrier, free_energy = self.calculate_free_energies(dGevb_ana, smooth=False)
            barrier_dif = self.barrier - barrier
            free_energy_dif = self.free_energy - free_energy
            # print(f"Barrier: {barrier}, difference: {barrier_dif}, Free energy: {free_energy}, difference: {free_energy_dif}")
            return barrier_dif, free_energy_dif

        alpha, H12 = scipy.optimize.fsolve(get_barrier_and_free_energy_difference, [self.alpha_guess, self.H12_guess])
        print(f"Fitted alpha: {alpha}, H12: {H12}")
        return alpha, H12

    def calculate_Eg_V_dE(self, E1, E2, alpha, H12, lambda_frame):
        E2_shifted = np.copy(E2) + alpha
        V = (1 - lambda_frame) * E1 + lambda_frame * E2_shifted
        dE = E1 - E2_shifted
        Eg = 0.5 * ((E1 + E2_shifted) - np.sqrt((E1 - E2_shifted)**2 + 4 * H12**2))
        return E2_shifted, V, dE, Eg

    def calculate_dGfep(self, dE, Temp, Lambda, Lambda_indices):
        de_lambda = self.bin(dE, Lambda_indices)
        T_lambda = self.bin(Temp, Lambda_indices)
        dG_middle, dG_forward, dG_backward, dG_bar = [0.0], [0.0], [0.0], [0.0]
        for i, l in enumerate(Lambda[:-1]):
            delta_lambda = Lambda[i + 1] - l

            # if remove_correlation:
            #     fw, t0fw = timeseries_analysis(de_l[i])
            #     bw, t0bw = timeseries_analysis(de_l[i + 1])
            #     print(
            #         f"{t0fw} {len(fw)-len(de_l[i])}     {t0bw} {len(bw)-len(de_l[i+1])}"
            #     )
            #     fw = beta * dl * fw
            #     bw = beta * dl * bw
            # else:
            forward_energy = self.beta(T_lambda[i]) * delta_lambda * de_lambda[i]
            forward_avg_temp = np.average(T_lambda[i])
            backward_energy = self.beta(T_lambda[i + 1]) * delta_lambda * de_lambda[i + 1]
            backward_avg_temp = np.average(T_lambda[i + 1])
            avg_temp = (forward_avg_temp + backward_avg_temp) / 2
            try:
                average_forward = np.average(np.exp(forward_energy / 2))
                average_backward = np.average(np.exp(-backward_energy / 2))

                dg_middle = -1 / self.beta(avg_temp) * math.log(average_forward / average_backward)
                dg_forward = -1 / self.beta(forward_avg_temp) * math.log(np.average(np.exp(forward_energy)))
                dg_backward = 1 / self.beta(backward_avg_temp) * math.log(np.average(np.exp(-backward_energy)))
            except ValueError:
                print(
                    f"ValueError encountered during FEP calculation, setting all dG for middle forward and backward to 0 for lambda {l}"
                )
                dg_middle = 0
                dg_forward = 0
                dg_backward = 0
            try:
                dg_bar = (-1 / self.beta(avg_temp) *
                            pymbar.other_estimators.bar(forward_energy, -backward_energy, False)["Delta_f"])
            except Exception as e:
                print(f"Error {e} encountered during BAR calculation, setting dG_bar to 0 for lambda {l}")
                dg_bar = 0
            dG_middle.append(dG_middle[-1] + dg_middle)
            dG_forward.append(dG_forward[-1] + dg_forward)
            dG_backward.append(dG_backward[-1] + dg_backward)
            dG_bar.append(dG_bar[-1] + dg_bar)
            # dG = dGfep_all[-1, :] + [middle, manfw, manbw, bar]
            # dGfep_all = np.row_stack((dGfep_all, dG))
        return dG_bar, dG_middle, dG_forward, dG_backward

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
            _, V, dE, Eg = self.calculate_Eg_V_dE(E1_ref, E2_ref, self.alpha, self.H12, Lambda_frame)

            E1_avg = np.average(self.bin(E1_ref,Lambda_indices),axis =1)
            E2_avg = np.average(self.bin(E2_ref,Lambda_indices),axis =1)
            Eg_avg = np.average(self.bin(Eg,Lambda_indices),axis =1)
            V_avg = np.average(self.bin(V,Lambda_indices),axis =1)
            dE_avg = np.average(self.bin(dE,Lambda_indices),axis =1)
            result.update({
                "dE": dE,
                "Eg": Eg,
                "V": V,
                "Lambda": Lambda,
                "E1_avg": E1_avg,
                "E2_avg": E2_avg,
                "Eg_avg": Eg_avg,
                "V_avg": V_avg,
                "dE_avg": dE_avg,
                "Temp": Temp,
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
            Temp = result["Temp"]
            Lambda = result["Lambda"]
            Lambda_indices = result["Lambda_indices"]

            dGfep, _, _, _ = self.calculate_dGfep(dE, Temp, Lambda, Lambda_indices) 

            result.update({"dGfep": dGfep})
            if self.calculate_discrete:
                dGevb_discrete, pns, dGcor= self.calculate_dGevb_discretised(
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

            # print(np.shape(L))
            # print(np.shape(dE))
            # print(np.shape(dens))
            ax[j,0].scatter(dE,L_values,c=dens,s=5)
            ax[j,0].plot([dE_min,dE_min],[0,0.3])
            ax[j,0].plot([dE_max,dE_max],[0.7,1])
            ax[j,0].set_ylabel(r"$\lambda$")
            ax[j,0].set_xlabel(r"$\Delta \mathcal{E}$ (kcal/mol)")
            ax[j,0].set_ylim(0,1)
            ax[j,0].set_xlim(min(dE)*1.1,max(dE)*1.1)


            min_label = min(dE)
            min_label = min_label -min_label%200 + 200
            max_label = max(dE)
            max_label = max_label -max_label%200 + 200
            xlabels = np.arange(min_label,max_label,200)
            ylabels = np.arange(0,1.1,0.2)

            min_tick = min(dE)
            min_tick = min_tick -min_tick%100 + 100
            max_tick = max(dE)
            max_tick = max_tick -max_tick%100 + 100

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


    def print_results(self):
        print(f"{'Discrete':<30}\t Barrier \t\t Free Energy")
        for name, result in self.results.items():
            if "discrete" in result.keys():
                print(f"{name:<30} \t {result['discrete']['barrier']} \t {result['discrete']['free_energy']}")

        print("\n")
        print("Analytical\t Barrier \t\t Free Energy")
        for name, result in self.results.items():
            if "analytical" in result.keys():
                print(f"{name:<30} \t {result['analytical']['barrier']} \t {result['analytical']['free_energy']}")

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
            E1_plot = result["E1_avg"] - np.min(result["E1_avg"])
            E2_plot = result["E2_avg"] - np.min(result["E1_avg"])
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

            pnsfep = pns@dGfep
            pnscor = np.sum(pns*dGcor.transpose(),axis=1)
            dGevb_disc = pnsfep+pnscor

            ax[i].plot(self.coordinate_bins, pnsfep[:-1], colors['tab:blue'], linestyle="--")
            ax[i].plot(self.coordinate_bins, pnscor[:-1], colors['tab:blue'], linestyle=":")
            ax[i].plot(self.coordinate_bins, dGevb_disc[:-1], colors['tab:blue'])

            #analytical curves
            shift = result['analytical']['shift']
            fepxi = result['analytical']['fep']
            dGevb_ana = shift+fepxi

            ax[i].plot(self.coordinate_bins, shift, colors['tab:orange'], linestyle=":")
            ax[i].plot(self.coordinate_bins, fepxi, colors['tab:orange'], linestyle="--")
            ax[i].plot(self.coordinate_bins, dGevb_ana, colors['tab:orange'])

            ax[i].set_title(name)

        fig.legend([r"$p_{n,s} \Delta G_{FEP}$", r"$p_{n,s} \Delta G_{cor}$", r"$\Delta G_{EVB,disc.}$", r"$\mu$", r"$\nu$", r"$\Delta G_{EVB,ana.}$",])
            

    # def show_snapshots(folder):
    #     import py3Dmol
    #     # Read the pdb file and split it into models
    #     with open(f"{folder}/traj_combined.pdb", "r") as file:
    #         models = file.read().split("ENDMDL")

    #     # Extract the first and last model
    #     first_model = models[0] + "ENDMDL"
    #     last_model = models[-2] + "ENDMDL"  # -2 because the last element is an empty string

    #     # Display the first model
    #     view = py3Dmol.view(width=400, height=300)
    #     view.addModel(first_model, "pdb", {"keepH": True})
    #     view.setStyle({}, {"stick": {}, "sphere": {"scale": 0.25}})
    #     view.zoomTo()
    #     view.show()

    #     # Display the last model
    #     view = py3Dmol.view(width=400, height=300)
    #     view.addModel(last_model, "pdb", {"keepH": True})
    #     view.setStyle({}, {"stick": {}, "sphere": {"scale": 0.25}})
    #     view.zoomTo()
    #     view.show()

    # @staticmethod
    # def discard_data(discard, E1, E2, lambda_frame):
    #     del_indices = np.array([])
    #     unique_lambda_frame = np.unique(lambda_frame)
    #     discard_num = int(len(lambda_frame) * discard / len(unique_lambda_frame))

    #     for l in unique_lambda_frame:
    #         # all indices where lambda_frame = l
    #         indices = np.where(lambda_frame == l)

    #         new_del_indices = np.round(np.linspace(indices[0][0], indices[0][-1] - 1, discard_num))

    #         del_indices = np.append(del_indices, new_del_indices)
    #     del_indices = np.array(del_indices, dtype=int)
    #     return (
    #         np.delete(E1, del_indices),
    #         np.delete(E2, del_indices),
    #         np.delete(lambda_frame, del_indices),
    #     )

    # def plot_energies(folder):
    #     E_file = f"{folder}/Energies.dat"
    #     options_file = f"{folder}/options.json"
    #     with open(options_file, "r") as file:
    #         options = json.load(file)
    #     Lambda = options["Lambda"]
    #     ETV_file = f"{folder}/Data_combined.dat"

    #     E1_ref, E2_ref, E1_run, E2_run, E_m = load_energies(E_file)
    #     steps, E, T, V, lambda_frame = load_ETV(ETV_file)

    #     _, V_ref, dE_ref, Eg_ref = calculate_Eg_V_dE(E1_ref, E2_ref, 0, 0, lambda_frame)
    #     _, V_run, dE_run, Eg_run = calculate_Eg_V_dE(E1_run, E2_run, 0, 0, lambda_frame)

    #     E1f_ref_file = f"{folder}/E1f_ref.dat"
    #     E2f_ref_file = f"{folder}/E2f_ref.dat"
    #     E1f_run_file = f"{folder}/E1f_run.dat"
    #     E2f_run_file = f"{folder}/E2f_run.dat"
    #     Vf_ref, ref_headers = get_Vf(E1f_ref_file, E2f_ref_file, lambda_frame)
    #     Vf_run, run_headers = get_Vf(E1f_run_file, E2f_run_file, lambda_frame)

    #     Efm_file = f"{folder}/Efm.dat"
    #     Efm = np.loadtxt(Efm_file)

    #     # fig = plt.figure(layout="constrained", figsize=(10, 10))
    #     # space = 0.1
    #     # subfigs = fig.subfigures(2, 2, wspace=space, hspace=space)

    #     # ETV_ax = plot_ETV(subfigs[0, 0], E, T, V, steps)
    #     # joule_to_cal: float = 0.239001
    #     # V_ax = plot_V(subfigs[1, 0], E * joule_to_cal, E_m, V_ref, V_run, lambda_frame)
    #     # diabats_ax = plot_diabats(
    #     #     subfigs[0, 1],
    #     #     dE_ref,
    #     #     E1_ref,
    #     #     E2_ref,
    #     #     Eg_ref,
    #     #     dE_run,
    #     #     E1_run,
    #     #     E2_run,
    #     #     Eg_run,
    #     # )
    #     # fc_ax = plot_force_contributions(
    #     #     subfigs[1, 1], Efm, Vf_ref, Vf_run, lambda_frame, ref_headers
    #     # )

    #     fig = plt.figure(layout="constrained", figsize=(8, 4))

    #     diabats_ax = plot_diabats(
    #         fig,
    #         dE_ref,
    #         E1_ref,
    #         E2_ref,
    #         Eg_ref,
    #         dE_run,
    #         E1_run,
    #         E2_run,
    #         Eg_run,
    #     )

    #     handles = [
    #         Line2D(
    #             [0],
    #             [0],
    #             marker="o",
    #             color="w",
    #             label="Scatter",
    #             markerfacecolor=mcolors.TABLEAU_COLORS["tab:blue"],
    #             markersize=5,
    #         ),
    #         Line2D(
    #             [0],
    #             [0],
    #             marker="o",
    #             color="w",
    #             label="Scatter",
    #             markerfacecolor=mcolors.TABLEAU_COLORS["tab:orange"],
    #             markersize=5,
    #         ),
    #         Line2D(
    #             [0],
    #             [0],
    #             marker="o",
    #             color="w",
    #             label="Scatter",
    #             markerfacecolor=mcolors.TABLEAU_COLORS["tab:green"],
    #             markersize=5,
    #         ),
    #     ]
    #     labels = [r"$\mathcal{E}_1$", r"$\mathcal{E}_2$", r"$E_g$"]
    #     fig.legend(handles, labels, ncols=3, loc="upper center", bbox_to_anchor=(0.5, 1.1))
    #     fig.tight_layout(pad=2.0)
    #     return fig

    # def plot_ETV(fig, E, T, V, steps):
    #     ax = fig.subplots(3, 1, sharex=True)
    #     # Plot E
    #     ax[0].plot(steps, E)
    #     ax[0].set_ylabel("E")

    #     # Plot T
    #     ax[1].plot(steps, T)
    #     ax[1].set_ylabel("T")

    #     # Plot V
    #     ax[2].plot(steps, V)
    #     ax[2].set_xlabel("Steps")
    #     ax[2].set_ylabel("V")
    #     return ax

    # # def average_per_lambda(lambda_frame, quant):
    # #     average = []
    # #     unique_lambda_frame = np.unique(lambda_frame)
    # #     for lam in unique_lambda_frame:
    # #         average.append(np.mean(quant[lambda_frame == lam]))
    # #     return np.array(average), unique_lambda_frame

    # def plot_V(fig, E, Em, V_ref, V_run, lambda_frame):
    #     ax = fig.subplots(1, 1)
    #     E_average, Lambda = average_per_lambda(lambda_frame, E)
    #     Em_average, _ = average_per_lambda(lambda_frame, Em)
    #     V_ref_average, _ = average_per_lambda(lambda_frame, V_ref)
    #     V_run_average, _ = average_per_lambda(lambda_frame, V_run)
    #     # ax.plot(Lambda, E_average, label="E", linewidth=3)
    #     ax.plot(Lambda, E_average - Em_average, label="Em", linewidth=2)
    #     ax.plot(Lambda, E_average - V_ref_average, label="V_ref", linewidth=1)
    #     ax.plot(Lambda, E_average - V_run_average, label="V_run", linewidth=1)

    #     # ax.plot(steps, E, linewidth=2, label="E")
    #     # opacity = 0.8
    #     # ax.plot(steps, Em, linewidth=1, alpha=opacity, label="Em")
    #     # ax.plot(steps, V_run, linewidth=0.5, alpha=opacity, label="V_run")
    #     # ax.plot(steps, V_ref, linewidth=0.5, alpha=opacity, label="V_ref")
    #     ax.legend()
    #     ax.set_xlabel(r"$\lambda$")
    #     ax.set_ylabel("Energy (kcal/mol)")
    #     return ax

    # def plot_diabats(fig, dE_ref, E1_ref, E2_ref, Eg_ref, dE_run, E1_run, E2_run, Eg_run):
    #     ax = fig.subplots(1, 3)
    #     dotsize = 0.1
    #     ax[0].scatter(dE_run, E1_run, s=dotsize * 10)
    #     ax[0].scatter(dE_run, E2_run, s=dotsize * 10)
    #     ax[0].scatter(dE_run, Eg_run, s=dotsize, alpha=0.8)
    #     ax[0].set_ylabel("E (kcal/mol)")
    #     ax[0].text(
    #         0.5,
    #         0.9,
    #         r"$V_{\mathrm{sample}}$",
    #         horizontalalignment="center",
    #         transform=ax[0].transAxes,
    #     )
    #     ax[0].set_xlabel(r"$\Delta \mathcal{E}$ (kcal/mol)")
    #     ax[0].set_xticks([-1500, 0, 1500])
    #     ax[1].scatter(dE_run, E1_run, s=dotsize * 10)
    #     ax[1].scatter(dE_run, E2_run, s=dotsize * 10)
    #     ax[1].scatter(dE_run, Eg_run, s=dotsize, alpha=0.8)
    #     ax[1].set_ylabel("E (kcal/mol)")
    #     ax[1].text(
    #         0.5,
    #         0.9,
    #         r"$V_{\mathrm{sample}}$",
    #         horizontalalignment="center",
    #         transform=ax[1].transAxes,
    #     )

    #     ax[1].set_xlabel(r"$\Delta \mathcal{E}$ (kcal/mol)")
    #     ax[1].set_xlim(-300, 300)
    #     ax[1].set_ylim(-20, 300)
    #     ax[2].scatter(dE_ref, E1_ref, s=dotsize * 10)
    #     ax[2].scatter(dE_ref, E2_ref, s=dotsize * 10)
    #     ax[2].scatter(dE_ref, Eg_ref, s=dotsize, alpha=0.8)
    #     ax[2].set_xlabel(r"$\Delta \mathcal{E}$ (kcal/mol)")
    #     ax[2].set_ylabel("E (kcal/mol)")
    #     ax[2].text(
    #         0.5,
    #         0.9,
    #         r"$V_{\mathrm{recalc.}}$",
    #         horizontalalignment="center",
    #         transform=ax[2].transAxes,
    #     )

    #     ax[2].set_xlim(-300, 300)
    #     ax[2].set_ylim(-20, 300)
    #     return ax

    # def get_Vf(E1f_file, E2f_file, lambda_frame):
    #     E1f = np.loadtxt(E1f_file, skiprows=1)
    #     E2f = np.loadtxt(E2f_file, skiprows=1)
    #     with open(E1f_file, "r") as file:
    #         headers = file.readline().split(",")
    #     Lambda_frame_tiled = np.tile(lambda_frame, (len(headers), 1)).transpose()
    #     Vf = (1 - Lambda_frame_tiled) * E1f + Lambda_frame_tiled * E2f
    #     return Vf, headers

    # def plot_force_contributions(fig, Efm, Vf_ref, Vf_run, lambda_frame, headers):
    #     ax = fig.subplots(1, 1)
    #     tol = 0.1
    #     opacity = 1
    #     for i, force_name in enumerate(headers):
    #         ref_dif = Efm[:, i] - Vf_run[:, i]
    #         ref_dif_avg, Lambda = average_per_lambda(lambda_frame, ref_dif)
    #         run_dif = Efm[:, i] - Vf_ref[:, i]
    #         run_dif_avg, Lambda = average_per_lambda(lambda_frame, run_dif)

    #         if np.max(np.abs(ref_dif_avg)) > tol:
    #             ax.plot(Lambda, ref_dif_avg, label=f"Ref {force_name}", alpha=opacity)
    #             opacity = 0.7
    #         if np.max(np.abs(run_dif_avg)) > tol:
    #             ax.plot(Lambda, run_dif_avg, label=f"Run {force_name}", alpha=opacity)
    #             opacity = 0.7
    #     ax.set_xlabel(r"$\lambda$")
    #     ax.set_ylabel("Energy (kcal/mol)")
    #     ax.legend()
    #     return ax
