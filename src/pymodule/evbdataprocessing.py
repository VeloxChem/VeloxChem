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
from .reactionsystembuilder import EvbForceGroup
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical, print_exception_if_debug

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

        # "bar": pairwise BAR between adjacent lambda windows (default).
        # "mbar": joint multistate estimate using every window's overlap with
        # every other window, not just its neighbors.
        self.fep_estimator: str = "mbar"

        # Boltzmann constant: ca 8.314e-3 kJ/mol/K
        self.kb = boltzmann_in_hartreeperkelvin() * hartree_in_kjpermol()
        self.verbose: bool = True

        self.calculate_discrete = True
        self.calculate_analytical = True
        # Per-(replica, direction) FEP/EVB curves + uncertainty, used for
        # hysteresis analysis and per-replica plotting. Only computed when
        # the loaded results carry "Replica_frame"/"Direction_frame" (see
        # EvbDriver._load_output_files). Set to False to skip -- this reruns
        # BAR/MBAR per replica/direction, so it can be expensive for runs
        # with many replicas.
        self.calculate_per_replica = True
        # Also compute each replica's "both" (forward+backward pooled)
        # curve, on top of the forward-only/backward-only curves hysteresis
        # needs. Set to False to skip it and cut per-replica work by ~1/3
        # when only the hysteresis check matters, not the per-replica
        # "both" curve used by print_uncertainties/plot_evb_replicas.
        self.calculate_per_replica_pooled = True
        # Estimator used for the per-(replica, direction) subsets in
        # _calculate_per_replica_results. None (default) reuses whatever
        # fep_estimator is set to. These per-replica curves are primarily a
        # diagnostic (uncertainty/hysteresis), so forcing "bar" here can be
        # much faster than "mbar" when fep_estimator="mbar"
        self.per_replica_fep_estimator = 'bar'

        # Warm-start cache for _calculate_dGfep_mbar
        self._mbar_f_k_warm_start = None
        # The Savitzky-Golay filter only locates the extrema of an EVB
        # profile; the reported energies are read from the unsmoothed
        # profile (see _get_free_energies), so a wide window no longer
        # flattens the barrier.
        self.smooth_window_size = 10
        self.smooth_polynomial_order = 3
        # Energy scale of a feature worth calling a well or a barrier, in
        # kJ/mol (about kT at 300 K). It is both the prominence a candidate
        # extremum must have to survive the noise filter and the depth over
        # which the parabola that refines an extremum is fitted.
        self.extremum_tolerance = 2.5
        self.coordinate_bins = np.array([])
        self.bin_size = 10
        self.dens_threshold = 0.1
        self.dens_max_window = 50
        # gaussian_kde in _calculate_coordinate_bins is O(n_frames^2); replica
        # x direction sampling can push n_frames into the tens of thousands,
        # so cap how many frames are fed to it (kept above the ~6k frames/config
        # in the canned test fixtures, so that regression test is unaffected).
        self.dens_kde_max_samples: int = 8000
        self.kde_seed = 0

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
        results.update({'alpha': self.alpha, 'H12': self.H12})
        self.ostream.print_info("Calculating FEP and EVB curves")
        self.ostream.flush()
        self._get_FEP_and_EVB()
        results["coordinate_bins"] = self.coordinate_bins
        return results

    def _fit_EVB_parameters(self):

        assert_msg_critical('scipy' in sys.modules,
                            'scipy is required for EvbDataProcessing.')

        reference_key = list(self.results.keys())[0]
        self.ostream.print_info(
            f"Configuration with key {reference_key} is used as reference for fitting"
        )

        E1_ref = self.results[reference_key]["E1_pes"]
        E2_ref = self.results[reference_key]["E2_pes"]
        Temp_set = self.results[reference_key]["Temp_set"]

        # fsolve does not raise on a NaN residual, it just stops and hands
        # back whatever it has, so broken input would come out as plausible
        # looking alpha and H12. Catch it here instead.
        assert_msg_critical(
            np.all(np.isfinite(E1_ref)) and np.all(np.isfinite(E2_ref)),
            f"EvbDataProcessing: configuration '{reference_key}' holds " +
            "non-finite state energies, so alpha and H12 cannot be fitted.")

        xi = np.linspace(-10000, 10000, 20000)

        def get_barrier_and_free_energy_difference(x):
            alpha, H12 = x

            E2_shifted, V, dE, Eg = self._calculate_Eg_V_dE(
                E1_ref, E2_ref, alpha, H12)
            dGfep, _ = self._calculate_dGfep(E1_ref, E2_shifted, Temp_set)
            dGevb_ana, shiftxi, fepxi = self._dGevb_analytical(
                dGfep, self.Lambda, H12, xi)
            dGevb_smooth, barrier, free_energy, _, _, _, _, _ = self._get_free_energies(
                dGevb_ana, xi, fitting=True)
            assert_msg_critical(
                np.isfinite(barrier) and np.isfinite(free_energy),
                "EvbDataProcessing: the EVB profile has no barrier or no " +
                f"reaction free energy at alpha = {alpha:.3f}, " +
                f"H12 = {H12:.3f}, so there is nothing to fit to. Provide " +
                "both alpha and H12, or check the sampling.")
            barrier_dif = self.barrier - barrier
            free_energy_dif = self.free_energy - free_energy

            return barrier_dif, free_energy_dif

        alpha, H12 = scipy.optimize.fsolve(
            get_barrier_and_free_energy_difference,
            [self.alpha_guess, self.H12_guess],
            xtol=self.tol)
        self.ostream.print_info(f"Fitted alpha: {alpha}, H12: {H12}")
        return alpha, H12

    def _calculate_Eg_V_dE(self, E1, E2, alpha, H12):
        E2_shifted = np.copy(E2) + alpha
        V = (1 - self.Lambda_frame) * E1 + self.Lambda_frame * E2_shifted
        dE = E1 - E2_shifted
        Eg = 0.5 * (
            (E1 + E2_shifted) - np.sqrt((E1 - E2_shifted)**2 + 4 * H12**2))
        return E2_shifted, V, dE, Eg

    def _calculate_dGfep(self, E1, E2_shifted, Temp_set, lambda_indices=None):
        assert_msg_critical(
            self.fep_estimator in ("bar", "mbar"),
            f"Unknown fep_estimator '{self.fep_estimator}'. Expected 'bar' or 'mbar'."
        )
        if self.fep_estimator == "mbar":
            return self._calculate_dGfep_mbar(E1, E2_shifted, Temp_set,
                                              lambda_indices)
        return self._calculate_dGfep_bar(E1 - E2_shifted, Temp_set,
                                         lambda_indices)

    def _calculate_dGfep_bar(self, dE, Temp_set, lambda_indices=None):

        assert_msg_critical('pymbar' in sys.modules,
                            'pymbar is required for EvbDataProcessing.')

        de_lambda = self._bin(dE, lambda_indices)
        dG_bar = [0.0]
        ddG_bar = [0.0]
        for i, l in enumerate(self.Lambda[:-1]):
            delta_lambda = self.Lambda[i + 1] - l

            # Reduced works in pymbar's convention, with dE = E1 - E2' so that
            # u_{i+1} - u_i = -beta * delta_lambda * dE: forward work on
            # window-i frames, reverse work on window-(i+1) frames. Passing
            # the negated works instead gives the same estimate only for
            # equal frame counts, and a wrong uncertainty in any case.
            forward_work = -self._beta(Temp_set) * delta_lambda * de_lambda[i]
            reverse_work = self._beta(Temp_set) * delta_lambda * de_lambda[i +
                                                                           1]

            try:
                bar_result = pymbar.other_estimators.bar(
                    forward_work, reverse_work)
                dF = bar_result["Delta_f"]
                ddF = bar_result["dDelta_f"]
                dg_bar = 1 / self._beta(Temp_set) * dF
                ddg_bar = 1 / self._beta(Temp_set) * ddF
            except Exception as e:
                print_exception_if_debug()
                self.ostream.print_warning(
                    f"Error {e} encountered during BAR calculation, setting dG_bar to 0 for lambda {l}"
                )
                dg_bar = 0
                ddg_bar = 0

            dG_bar.append(dG_bar[-1] + dg_bar)
            # Adjacent-window BAR estimates are independent, so their
            # variances add along the cumulative sum.
            ddG_bar.append(np.sqrt(ddG_bar[-1]**2 + ddg_bar**2))

        return dG_bar, ddG_bar

    def _calculate_dGfep_mbar(self,
                              E1,
                              E2_shifted,
                              Temp_set,
                              lambda_indices=None):
        assert_msg_critical('pymbar' in sys.modules,
                            'pymbar is required for EvbDataProcessing.')

        if lambda_indices is None:
            lambda_indices = self.Lambda_indices

        beta = self._beta(Temp_set)
        Lambda_indices = np.asarray(lambda_indices)
        Lambda_arr = np.asarray(self.Lambda)
        K = len(Lambda_arr)

        # Group frames by originating state, matching pymbar's N_k convention
        # (the first N_k[0] samples are from state 0, the next N_k[1] from
        # state 1, and so on).
        order = np.argsort(Lambda_indices, kind="stable")
        E1_sorted = np.asarray(E1)[order]
        E2_sorted = np.asarray(E2_shifted)[order]
        N_k = np.array([np.sum(Lambda_indices == k) for k in range(K)])

        # u_kn[k, n] = beta * V_k(frame_n). Cheap to get every frame's energy
        # at every state because V is linear in Lambda (Vi = (1-Lambda_i)*E1 +
        # Lambda_i*E2_shifted), so no re-evaluation at other windows is
        # needed, unlike a typical MBAR use case.
        u_kn = beta * ((1 - Lambda_arr)[:, None] * E1_sorted[None, :] +
                       Lambda_arr[:, None] * E2_sorted[None, :])

        # pymbar's default solver protocol tries 'hybr' (scipy's hybrid
        # Powell method) first, which reliably fails to converge once there
        # are this many states (tens of lambda windows) and only then falls
        # back to 'adaptive' -- same convention as solvationfepdriver.py.
        # Going straight to the 'robust' protocol (adaptive, then L-BFGS-B)
        # skips the doomed-to-fail hybr attempt: same converged result,
        # faster, and no spurious "failed to reach a solution" warning.
        initial_f_k = self._mbar_f_k_warm_start
        if initial_f_k is not None and len(initial_f_k) != K:
            # Mismatched state count (different Lambda ladder / caller) --
            # not a valid warm start, fall back to MBAR's own default.
            initial_f_k = None
        mbar = pymbar.MBAR(u_kn,
                           N_k,
                           solver_protocol='robust',
                           initial_f_k=initial_f_k)
        mbar_results = mbar.compute_free_energy_differences()
        # Cache for the next call (e.g. the next alpha tried during fitting,
        # or the next replica/configuration) -- only changes how fast the
        # solver converges, not what it converges to.
        self._mbar_f_k_warm_start = mbar.f_k
        dGfep = (mbar_results['Delta_f'][0, :] / beta).tolist()
        ddGfep = (mbar_results['dDelta_f'][0, :] / beta).tolist()
        return dGfep, ddGfep

    def _bin(self, data, lambda_indices=None):
        if lambda_indices is None:
            lambda_indices = self.Lambda_indices
        lambda_indices = np.asarray(lambda_indices)
        data = np.asarray(data)
        n_bins = np.max(lambda_indices) + 1

        # Vectorized equivalent of "for i, li in enumerate(lambda_indices):
        # binned_data[li].append(data[i])" -- a stable sort groups frames by
        # bin while preserving each bin's original relative frame order, so
        # this produces the exact same per-bin grouping, just without a
        # pure-Python per-frame loop.
        order = np.argsort(lambda_indices, kind="stable")
        sorted_data = data[order]
        sorted_indices = lambda_indices[order]
        counts = np.bincount(sorted_indices, minlength=n_bins)

        if np.all(counts == counts[0]):
            # Equal-sized bins (the common case): reshape directly into a
            # dense 2D array.
            binned_data = sorted_data.reshape(n_bins, counts[0])
        else:
            split_points = np.cumsum(counts)[:-1]
            binned_data = np.array(np.split(sorted_data, split_points),
                                   dtype=object)
        return binned_data

    def _bin_average(self, data, lambda_indices=None):
        """Mean of data per lambda window.

        _bin returns a ragged object array when windows hold different
        numbers of frames (e.g. after time sub-sampling), which
        np.average(..., axis=1) cannot handle.
        """
        binned_data = self._bin(data, lambda_indices)
        if binned_data.dtype == object:
            return np.array([np.average(b) for b in binned_data])
        return np.average(binned_data, axis=1)

    @staticmethod
    def _ground_state_shift_and_lambda(xi, H12):
        """The two xi-dependent pieces of the analytical EVB profile.

        shift = -2 H12^2 / R and lambda* = (1 + xi/R) / 2, with
        R = sqrt(xi^2 + 4 H12^2). R is 0 only for H12 = 0 at xi = 0, where
        both are 0/0: the shift tends to 0, and the two one-sided limits of
        lambda* meet at the midpoint 1/2.
        """

        xi = np.asarray(xi, dtype=float)
        R = np.sqrt(xi**2 + 4 * H12**2)
        positive = R > 0
        R_safe = np.where(positive, R, 1.0)
        shift = np.where(positive, -2 * H12**2 / R_safe, 0.0)
        lambda_star = 0.5 * (1 + np.where(positive, xi / R_safe, 0.0))
        return shift, lambda_star

    @staticmethod
    def _dGevb_analytical(dGfep, Lambda, H12, xi):

        shiftxi, lambda_star = (
            EvbDataProcessing._ground_state_shift_and_lambda(xi, H12))
        fepxi = np.interp(lambda_star, Lambda, dGfep)
        dGevb = shiftxi + fepxi
        return dGevb, shiftxi, fepxi

    @staticmethod
    def _dGevb_analytical_uncertainty(ddGfep, Lambda, H12, xi):
        """Propagate per-window dGfep uncertainty through the same
        interpolation _dGevb_analytical uses to build the EVB curve.

        shift(xi) is a deterministic function of H12 alone, so only the
        interpolated FEP term carries sampling uncertainty. Linear
        interpolation of the uncertainty is an approximation: it ignores
        correlation between neighbouring lambda windows.
        """

        _, lambda_star = (EvbDataProcessing._ground_state_shift_and_lambda(
            xi, H12))
        return np.interp(lambda_star, Lambda, ddGfep)

    @staticmethod
    def _propagate_free_energy_uncertainty(dGevb_unc, erea_ind, ebar_ind,
                                           epro_ind):
        """Combine pointwise EVB-curve uncertainty at the reactant/barrier/
        product indices in quadrature into barrier/free_energy uncertainty.

        Approximate: ignores correlation between curve points (they share
        the same underlying dGfep uncertainty they were interpolated from).
        """
        barrier_unc = np.sqrt(dGevb_unc[ebar_ind]**2 + dGevb_unc[erea_ind]**2)
        free_energy_unc = np.sqrt(dGevb_unc[epro_ind]**2 +
                                  dGevb_unc[erea_ind]**2)
        return barrier_unc, free_energy_unc

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
                            hist[n][s]))  # What to do with the temperature here
                    pnscount[n, s] = len(hist[n][s])
                else:
                    dGcor[n, s] = 0

        for n in range(N):
            for s in range(S):
                pnssum = np.sum(pnscount[:, s])
                if pnssum > 0:
                    pns[n, s] = pnscount[n, s] / pnssum

        # Window-density term -kT ln P_n(s), with P_n(s) = N_ns / N_n, the
        # fraction of window n's frames that fall in bin s. Without it the
        # sum below reduces to dG_n + (Eg - V_n)(s) per window, which is the
        # free energy of window n's whole ensemble rather than of bin s, and
        # the barrier comes out too low (0.5-2.0 kJ/mol on the exact models,
        # WP-A3 Q3).
        dGdens = np.zeros((N, S))
        window_count = np.sum(pnscount, axis=1)
        for n in range(N):
            if window_count[n] == 0:
                continue
            occupied = pnscount[n] > 0
            dGdens[n, occupied] = -self.kb * Temp_set * np.log(
                pnscount[n, occupied] / window_count[n])

        pns = pns.transpose()

        pnsfep = pns @ dGfep
        pnscor = np.sum(pns * (dGcor + dGdens).transpose(), axis=1)
        dGevb = pnsfep + pnscor
        # A bin no window ever visited has no free energy. Leaving it at 0
        # puts a point hundreds of kJ/mol away from the rest of the profile
        # wherever the sampling has a hole, which the extremum search then
        # picks up; NaN says "not sampled" and is skipped instead.
        dGevb[np.sum(pnscount, axis=0) == 0] = np.nan

        return dGevb, pns, dGcor, dGdens

    def _profile_axis(self, n_points, xi):
        """Coordinate of every point of a profile, and the range to search.

        The analytical profile is evaluated at the coordinate bins, one
        point per bin edge. The discretised profile has one point per bin,
        len(bins) + 1 of them: entry s covers (bins[s-1], bins[s]], and the
        first and last entries collect everything outside the grid, so they
        get nominal coordinates and are never searched for extrema.

        Returns the coordinate array and the first and last searchable
        index.
        """

        if xi is None:
            return np.arange(n_points, dtype=float), 0, n_points - 1
        xi = np.asarray(xi, dtype=float)
        if n_points == len(xi):
            return xi, 0, n_points - 1
        assert_msg_critical(
            n_points == len(xi) + 1,
            f"EvbDataProcessing: a profile of {n_points} points matches " +
            f"neither {len(xi)} coordinate bins nor the {len(xi) + 1} bins " +
            "they define.")
        centers = 0.5 * (xi[:-1] + xi[1:])
        axis = np.concatenate(([xi[0] - 0.5 * (xi[1] - xi[0])], centers,
                               [xi[-1] + 0.5 * (xi[-1] - xi[-2])]))
        return axis, 1, n_points - 2

    @staticmethod
    def _prominent_extrema(profile, lo, hi, is_maximum, prominence):
        """Local extrema between lo and hi, noise-sized ones dropped.

        Prominence is how far the profile has to climb away from an extremum
        before it reaches a deeper one, so it separates a real well or
        barrier from a wiggle on a slope. It is measured inside [lo, hi]
        only, so the out-of-range bins of a discretised profile cannot
        inflate it.

        find_peaks is used rather than argrelmax because a profile often has
        a flat top: binned data repeats values, and an even smoothing window
        turns the top of a symmetric barrier into two equal points, which a
        strict point-by-point comparison misses altogether.
        """

        window = np.asarray(profile[lo:hi + 1], dtype=float)
        signed = window if is_maximum else -window
        found, _ = scipy.signal.find_peaks(signed, prominence=prominence)
        return np.asarray(found, dtype=int) + lo

    def _extremum_window(self,
                         profile,
                         index,
                         lo,
                         hi,
                         is_maximum,
                         adaptive=True):
        """Indices to fit a parabola over: the points within
        extremum_tolerance of profile[index], at least one on each side.

        The width adapts to the profile: wide where it is flat, so that noise
        averages out, and narrow at a sharp barrier, where a wide parabola
        would cut the top off. With adaptive=False a fixed five points are
        used instead; see _refine_extremum.
        """

        if not adaptive:
            return np.arange(max(lo, index - 2), min(hi, index + 2) + 1)

        sign = -1.0 if is_maximum else 1.0
        left, right = index, index
        while (left - 1 >= lo and sign *
               (profile[left - 1] - profile[index]) < self.extremum_tolerance):
            left -= 1
        while (right + 1 <= hi and sign *
               (profile[right + 1] - profile[index]) < self.extremum_tolerance):
            right += 1
        left = min(left, max(lo, index - 1))
        right = max(right, min(hi, index + 1))
        return np.arange(left, right + 1)

    def _refine_extremum(self,
                         axis,
                         profile,
                         index,
                         lo,
                         hi,
                         is_maximum,
                         adaptive=True):
        """Energy, position and grid index of one extremum of the unsmoothed
        profile.

        The smoothed profile only says which point is the extremum; its
        energy comes from a parabola fitted to the raw profile around it.
        That removes the flattening of the Savitzky-Golay window and
        recovers most of what is lost when the true extremum falls between
        two grid points.

        The filter also displaces a peak by a bin or so, which matters when
        the barrier is sharp, so the fit is re-centered on the grid point
        nearest its own vertex and repeated. A step is only taken towards a
        more extreme raw point, so the walk cannot wander off into noise.
        The raw value is kept whenever the parabola curves the wrong way or
        puts its vertex outside the fitted window, which is what noise does.

        The index is returned because the walk may have moved it: everything
        that points back at the profile, the markers and the uncertainty,
        has to use the point the energy actually came from.

        adaptive=False fits five points once and does not walk. Both the
        width rule and the walk change in steps as the profile moves under
        them, which is no help on the dense, well-resolved grid the
        alpha/H12 fit evaluates on, and turns that fit's target into a
        staircase: measured over a 20 kJ/mol scan of H12, the largest step
        in the target is 16 times the median with the walk and 1.1 times
        without, for the same barrier to 0.05 kJ/mol.
        """

        sign = -1.0 if is_maximum else 1.0

        def more_extreme(candidate, current):
            return (lo <= candidate <= hi
                    and sign * (profile[candidate] - profile[current]) < 0)

        for _ in range(4 if adaptive else 1):
            window = self._extremum_window(profile, index, lo, hi, is_maximum,
                                           adaptive)
            if len(window) < 3:
                return profile[index], axis[index], index
            x = axis[window] - axis[index]
            a, b, c = np.polyfit(x, profile[window], 2)
            if a == 0 or (a > 0) == is_maximum:
                return profile[index], axis[index], index
            x_top = -0.5 * b / a
            spacing = 0.5 * (axis[min(index + 1, hi)] -
                             axis[max(index - 1, lo)])
            if abs(x_top) > spacing:
                # A fit that puts the extremum more than a bin away is not
                # describing the point it was centered on. Either the filter
                # located the wrong grid point, in which case stepping that
                # way lands on a more extreme one and the next fit is
                # centered properly, or the profile is too flat or too noisy
                # here for a parabola to say anything, and the raw value
                # stands.
                step = index + (1 if x_top > 0 else -1)
                if adaptive and more_extreme(step, index):
                    index = step
                    continue
                return profile[index], axis[index], index
            return c - 0.25 * b**2 / a, axis[index] + x_top, index
        return profile[index], axis[index], index

    def _get_free_energies(self, dGevb, xi=None, fitting=True):
        """Reactant, transition state and product of one EVB profile.

        Returns the profile shifted to the reactant well, the barrier and
        the reaction free energy, the extrema that survived the noise
        filter, and the three indices the energies were taken at.

        Smoothing only locates the extrema and never sets a reported energy:
        each of the three comes from a parabola fitted to the unsmoothed
        profile (_refine_extremum), so the filter window cannot flatten the
        barrier. The transition state is the highest maximum whose
        prominence clears extremum_tolerance, and the reactant and product
        are the lowest points on either side of it. A profile with no such
        maximum has no barrier to report, and says so with NaN and a
        warning, rather than returning a number read off noise.

        xi is the coordinate of the profile, used to fit the parabolas and
        to recognize the out-of-range bins of a discretised profile.

        With fitting=True, which is how the alpha/H12 fit calls this, the
        wells are always taken on either side of the highest point and a
        number is returned even when nothing clears the noise filter, and
        nothing is printed. The parabola is fitted over a fixed five points
        and does not walk, which keeps the target fsolve differentiates
        smooth (see _refine_extremum). A profile that is all NaN still gives
        a NaN barrier; _fit_EVB_parameters checks for that before it starts.
        """

        assert_msg_critical('scipy' in sys.modules,
                            'scipy is required for EvbDataProcessing.')

        dGevb = np.asarray(dGevb, dtype=float)
        full_axis, full_lo, full_hi = self._profile_axis(len(dGevb), xi)

        # Unsampled bins hold NaN (see _dGevb_discretised). Everything below
        # runs on the points that do hold a free energy, and the indices are
        # mapped back at the end, so a hole in the sampling is stepped over
        # rather than treated as a well or a barrier.
        sampled = full_lo + np.flatnonzero(
            np.isfinite(dGevb[full_lo:full_hi + 1]))
        if len(sampled) < 3:
            if not fitting:
                self.ostream.print_warning(
                    f"Only {len(sampled)} points of the EVB profile hold a "
                    "free energy, which is too few to find a barrier. Check "
                    "the sampling and the coordinate range.")
                self.ostream.flush()
            empty = np.array([], dtype=int)
            return (dGevb, np.nan, np.nan, empty, empty, full_lo, full_hi,
                    full_lo)
        profile, axis = dGevb[sampled], full_axis[sampled]
        lo, hi = 0, len(profile) - 1

        # A short profile (a narrow gap range, or coarse bins) cannot carry
        # the default 10-point window: it would span the whole reaction and
        # smear the barrier into the wells. Shrink it to a third of the
        # profile, which leaves the default untouched for the usual grids of
        # 30 points or more.
        window = int(min(self.smooth_window_size, max(3, len(profile) // 3)))
        order = int(min(self.smooth_polynomial_order, window - 1))
        try:
            smooth = scipy.signal.savgol_filter(profile, window, order)
        except ValueError:
            if not fitting:
                self.ostream.print_warning(
                    f"Could not apply Savitzky-Golay filter with window size "
                    f"{window} and polynomial order {order}. Using "
                    "unfiltered data.")
            smooth = profile

        min_arg = self._prominent_extrema(smooth, lo, hi, False,
                                          self.extremum_tolerance)
        max_arg = self._prominent_extrema(smooth, lo, hi, True,
                                          self.extremum_tolerance)

        # The transition state is found first, as the highest maximum that
        # stands clear of the noise, and the reactant and product are then
        # sought on either side of it.
        barrier_found = len(max_arg) > 0
        if barrier_found:
            ebar_ind = int(max_arg[np.argmax(smooth[max_arg])])
        else:
            # The fit still needs a target, so it gets the highest point of
            # the profile even when nothing stands out.
            ebar_ind = lo + int(np.argmax(smooth[lo:hi + 1]))

        def well(first, last):
            """The reactant or product well between two indices.

            The deepest minimum that cleared the noise filter is taken when
            there is one. Otherwise the lowest point of the stretch is used,
            which is what a shallow well, or one that runs into the edge of
            the sampled range, looks like. That fallback can land on a noisy
            end point, so it says so.
            """

            candidates = min_arg[(min_arg >= first) & (min_arg <= last)]
            if len(candidates) > 0:
                return int(candidates[np.argmin(smooth[candidates])])
            return first + int(np.argmin(smooth[first:last + 1]))

        if barrier_found or fitting:
            erea_ind = well(lo, ebar_ind)
            epro_ind = well(ebar_ind, hi)
        else:
            # Without a transition state there is nothing to split the
            # profile at, so the reaction free energy is quoted between the
            # ends of the sampled gap range.
            erea_ind, epro_ind = lo, hi

        if not fitting and (erea_ind in (lo, hi) or epro_ind in (lo, hi)):
            self.ostream.print_warning(
                "The reactant or product well of the EVB profile sits at the "
                "end of the sampled gap range, so it is an end point rather "
                "than a minimum, and a single noisy bin there sets the zero "
                "of the whole profile. Widen the gap range, or check the "
                "sampling at the ends.")

        def energy_at(index, is_maximum):
            # An extremum against the edge of the range is not a stationary
            # point, so there is nothing for a parabola to refine.
            if index <= lo or index >= hi:
                return profile[index], index
            value, _, settled = self._refine_extremum(axis,
                                                      profile,
                                                      index,
                                                      lo,
                                                      hi,
                                                      is_maximum,
                                                      adaptive=not fitting)
            return value, settled

        Erea, erea_ind = energy_at(erea_ind, False)
        Epro, epro_ind = energy_at(epro_ind, False)
        if barrier_found or fitting:
            Ebar, ebar_ind = energy_at(ebar_ind, True)
            barrier = Ebar - Erea
        else:
            barrier = np.nan
            self.ostream.print_warning(
                "No barrier was found in the EVB profile: no maximum stands "
                f"more than {self.extremum_tolerance} kJ/mol above its "
                "surroundings. The barrier is reported as NaN. A barrierless "
                "or inverted-region reaction looks like this, and so does a "
                "gap range that stops short of the transition state or bins "
                "too coarse to resolve it.")

        free_energy = Epro - Erea

        if not fitting and len(max_arg) > 1:
            self.ostream.print_warning(
                f"Found {len(max_arg)} maxima in the EVB profile instead of "
                "1, and took the highest one. Confirm the calculated extrema "
                "with the plot.")

        # The reported curve is the unsmoothed profile, so that it passes
        # through the reported barrier and reaction free energy.
        dGevb_shifted = dGevb - Erea
        self.ostream.flush()
        return (dGevb_shifted, barrier, free_energy,
                sampled[min_arg], sampled[max_arg], int(sampled[erea_ind]),
                int(sampled[epro_ind]), int(sampled[ebar_ind]))

    def _get_FEP_and_EVB(self):

        for result in self.results.values():
            E1_ref = result["E1_pes"]
            E2_ref = result["E2_pes"]
            E2_shifted, V, dE, Eg = self._calculate_Eg_V_dE(
                E1_ref, E2_ref, self.alpha, self.H12)

            E1_avg = self._bin_average(E1_ref)
            E2_avg = self._bin_average(E2_ref)
            Eg_avg = self._bin_average(Eg)
            V_avg = self._bin_average(V)
            dE_avg = self._bin_average(dE)
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

        for name, result in self.results.items():
            Temp_set = result["Temp_set"]
            # E2_shifted isn't persisted on the result dict (would otherwise
            # change its key set and break the h5-reference comparison test),
            # so recompute the cheap scalar shift instead of storing it.
            E1 = result["E1_pes"]
            E2_shifted = result["E2_pes"] + self.alpha

            dGfep, dGfep_unc = self._calculate_dGfep(E1, E2_shifted, Temp_set)

            result.update({"dGfep": dGfep, "dGfep_unc": dGfep_unc})
            if self.calculate_discrete:
                dGevb_discrete, pns, dGcor, dGdens = self._dGevb_discretised(
                    dGfep, result["Eg"], result["V"], result["dE"],
                    result["Temp_set"])

                (
                    dGevb_discrete,
                    barrier_discretised,
                    reaction_free_energy_discretised,
                    min_arg,
                    max_arg,
                    erea_ind,
                    epro_ind,
                    ebar_ind,
                ) = self._get_free_energies(dGevb_discrete,
                                            self.coordinate_bins,
                                            fitting=False)

                result.update({
                    "discrete": {
                        "EVB": dGevb_discrete,
                        "free_energy": reaction_free_energy_discretised,
                        "barrier": barrier_discretised,
                        "pns": pns,
                        "dGcor": dGcor,
                        "dGdens": dGdens,
                        "min_arg": min_arg,
                        "max_arg": max_arg,
                        "reactant_index": erea_ind,
                        "barrier_index": ebar_ind,
                        "product_index": epro_ind,
                    }
                })

            if self.calculate_analytical:
                self.ostream.print_info(
                    f"Calculating analytical EVB curve for configuration {name}"
                )
                self.ostream.flush()
                dGevb_analytical, shift, fepxi = self._dGevb_analytical(
                    result["dGfep"],
                    self.Lambda,
                    self.H12,
                    self.coordinate_bins,
                )
                dGevb_analytical_unc = self._dGevb_analytical_uncertainty(
                    dGfep_unc,
                    self.Lambda,
                    self.H12,
                    self.coordinate_bins,
                )

                (
                    dGevb_analytical,
                    barrier_analytical,
                    reaction_free_energy_analytical,
                    min_arg,
                    max_arg,
                    erea_ind,
                    epro_ind,
                    ebar_ind,
                ) = self._get_free_energies(dGevb_analytical,
                                            self.coordinate_bins,
                                            fitting=False)

                barrier_analytical_unc, free_energy_analytical_unc = (
                    self._propagate_free_energy_uncertainty(
                        dGevb_analytical_unc, erea_ind, ebar_ind, epro_ind))

                result.update({
                    "analytical": {
                        "EVB": dGevb_analytical,
                        "EVB_unc": dGevb_analytical_unc,
                        "shift": shift,
                        "fep": fepxi,
                        "free_energy": reaction_free_energy_analytical,
                        "free_energy_unc": free_energy_analytical_unc,
                        "barrier": barrier_analytical,
                        "barrier_unc": barrier_analytical_unc,
                        "min_arg": min_arg,
                        "max_arg": max_arg,
                        "reactant_index": erea_ind,
                        "barrier_index": ebar_ind,
                        "product_index": epro_ind,
                    }
                })

            if self.calculate_per_replica and result.get(
                    "Replica_frame") is not None:
                self.ostream.print_info(
                    f"Calculating per-replica FEP/EVB curves for configuration {name}"
                )
                self.ostream.flush()
                self._calculate_per_replica_results(name, result, E1,
                                                    E2_shifted, Temp_set)
                self.ostream.print_info(
                    f"Calculating per-replica hysteresis for configuration {name}"
                )
                self.ostream.flush()
                self._calculate_hysteresis(result)
                self._calculate_hysterisis_data(result)
                self._calculate_replica_average(result)

            if 'E1_fg' in result.keys():
                self.ostream.print_info(
                    f"Calculating per-force-group FEP/EVB curves for configuration {name}"
                )
                self.ostream.flush()
                self._calculate_fg_profiles(result)

    def _calculate_replica_average(self, result):
        """Mean and standard error of the mean of barrier/free_energy across replicas.

        Uses each replica's "both" (forward+backward pooled) analytical
        estimate from _calculate_per_replica_results as that replica's single
        independent data point, so the resulting mean/SEM captures
        replica-to-replica variability -- a different (and often more
        conservative) uncertainty source than the internal BAR/MBAR
        statistical error reported per curve, since it also picks up
        replica-to-replica effects like slow correlated fluctuations that a
        single replica's own BAR/MBAR estimate can't see.

        Requires calculate_per_replica_pooled=True (the default): with it off,
        no replica has a "both" curve to average, so this is skipped.
        """
        by_replica = result.get('replicas')
        if not by_replica:
            return

        barriers = []
        free_energies = []
        for directions in by_replica.values():
            if "both" not in directions:
                continue
            barriers.append(directions["both"]["analytical"]["barrier"])
            free_energies.append(
                directions["both"]["analytical"]["free_energy"])

        if not barriers:
            return

        barriers = np.asarray(barriers)
        free_energies = np.asarray(free_energies)
        n = len(barriers)
        # SEM needs at least 2 replicas to estimate a spread; report NaN
        # rather than a misleading 0.0 when there's only one.
        if n > 1:
            barrier_sem = float(np.std(barriers, ddof=1) / np.sqrt(n))
            free_energy_sem = float(np.std(free_energies, ddof=1) / np.sqrt(n))
        else:
            barrier_sem = float("nan")
            free_energy_sem = float("nan")

        result.update({
            "replica_average": {
                "n_replicas": n,
                "barrier_mean": float(np.mean(barriers)),
                "barrier_sem": barrier_sem,
                "free_energy_mean": float(np.mean(free_energies)),
                "free_energy_sem": free_energy_sem,
            }
        })

    def _calculate_hysterisis_data(self, result):
        """Calculate the average, median, min max and std of the hysteresis
        data for each configuration"""

        by_replica = result.get('replicas')
        if not by_replica or not any('hysteresis' in directions
                                     for directions in by_replica.values()):
            return

        metrics = (
            "max_abs_dGfep_delta",
            "mean_abs_dGfep_delta",
            "endpoint_dGfep_delta",
            "barrier_delta",
            "free_energy_delta",
        )

        values = {metric: [] for metric in metrics}
        replica_ids = {metric: [] for metric in metrics}
        for replica, directions in by_replica.items():
            hysteresis = directions.get("hysteresis")
            if not hysteresis:
                continue
            for metric in metrics:
                values[metric].append(hysteresis[metric])
                replica_ids[metric].append(replica)

        hysteresis_summary = {}
        for metric, data in values.items():
            data = np.asarray(data)
            metric_replicas = replica_ids[metric]
            median = float(np.median(data))
            # The median of an even-sized sample isn't necessarily any single
            # replica's value, so report whichever replica sits closest to it.
            median_replica = metric_replicas[int(
                np.argmin(np.abs(data - median)))]
            hysteresis_summary[metric] = {
                "mean":
                float(np.mean(data)),
                "median":
                median,
                "median_replica":
                median_replica,
                "min":
                float(np.min(data)),
                "min_replica":
                metric_replicas[int(np.argmin(data))],
                "max":
                float(np.max(data)),
                "max_replica":
                metric_replicas[int(np.argmax(data))],
                # std needs at least 2 replicas to estimate a spread; report
                # NaN rather than a misleading 0.0 when there's only one.
                "std":
                float(np.std(data, ddof=1)) if len(data) > 1 else float("nan"),
            }

        result.update({"hysteresis_summary": hysteresis_summary})

    def _calculate_hysteresis(self, result):
        """Per-replica forward-vs-backward reversibility check.

        Compares the forward-only and backward-only dGfep/EVB curves from
        _calculate_per_replica_results. Both are independently anchored to 0
        at lambda=0 (same BAR/MBAR convention as the pooled curve), so if
        sampling were perfectly reversible the forward and backward curves
        would coincide at every lambda; the gap between them is the
        hysteresis.
        """
        replicas = result.get('replicas')
        if not replicas:
            self.ostream.print_info(
                "No per-replica forward/backward curves available, skipping hysteresis analysis"
            )
            return

        for directions in replicas.values():
            if "forward" not in directions or "backward" not in directions:
                continue

            dGfep_fwd = np.asarray(directions["forward"]["dGfep"])
            dGfep_bwd = np.asarray(directions["backward"]["dGfep"])
            delta = dGfep_fwd - dGfep_bwd

            barrier_fwd = directions["forward"]["analytical"]["barrier"]
            barrier_bwd = directions["backward"]["analytical"]["barrier"]
            free_energy_fwd = directions["forward"]["analytical"]["free_energy"]
            free_energy_bwd = directions["backward"]["analytical"][
                "free_energy"]

            directions["hysteresis"] = {
                "dGfep_delta": delta,
                "max_abs_dGfep_delta": float(np.max(np.abs(delta))),
                "mean_abs_dGfep_delta": float(np.mean(np.abs(delta))),
                # Disagreement in the overall (lambda=0 -> 1) free energy
                # estimate between the forward-only and backward-only frames.
                "endpoint_dGfep_delta": float(delta[-1]),
                "barrier_delta": barrier_fwd - barrier_bwd,
                "free_energy_delta": free_energy_fwd - free_energy_bwd,
            }

    def _calculate_per_replica_results(self, name, result, E1, E2_shifted,
                                       Temp_set):
        """Per-(replica, direction) FEP/EVB curves + uncertainty.

        Feeds the hysteresis analysis (_calculate_hysteresis) and the
        per-replica plots (plot_evb_replicas / plot_fep_hysteresis). Each
        subset reruns the same BAR/MBAR estimator used for the pooled curve,
        restricted to that subset's frames -- every frame carries both
        end-state energies regardless of which lambda window it was sampled
        from, so this is a frame-subsetting exercise, not a different
        estimator.
        """
        replica_frame = np.asarray(result["Replica_frame"])
        direction_frame = np.asarray(result["Direction_frame"])
        lambda_indices = np.asarray(self.Lambda_indices)

        subsets = [("forward", 0), ("backward", 1)]
        if self.calculate_per_replica_pooled:
            subsets.append(("both", None))

        # These per-replica curves are a diagnostic (uncertainty/hysteresis),
        # not the main reported free energy, so allow forcing the cheaper
        # BAR estimator here even when fep_estimator="mbar" is used for the
        # pooled curve -- MBAR's joint multistate solve is otherwise repeated
        # from scratch for every one of the (up to) 3 * n_replicas subsets.
        saved_estimator = self.fep_estimator
        if self.per_replica_fep_estimator is not None:
            self.fep_estimator = self.per_replica_fep_estimator

        try:
            by_replica = {}
            for replica in sorted(set(replica_frame.tolist())):
                by_replica[replica] = {}
                for label, direction in subsets:
                    if direction is None:
                        mask = replica_frame == replica
                    else:
                        mask = (replica_frame == replica) & (direction_frame
                                                             == direction)

                    if not np.any(mask):
                        continue

                    try:
                        dGfep_sub, dGfep_unc_sub = (self._calculate_dGfep(
                            E1[mask], E2_shifted[mask], Temp_set,
                            lambda_indices[mask]))

                        dGevb_sub, _, _ = self._dGevb_analytical(
                            dGfep_sub, self.Lambda, self.H12,
                            self.coordinate_bins)
                        dGevb_unc_sub = self._dGevb_analytical_uncertainty(
                            dGfep_unc_sub, self.Lambda, self.H12,
                            self.coordinate_bins)

                        (
                            dGevb_sub,
                            barrier_sub,
                            free_energy_sub,
                            min_arg_sub,
                            max_arg_sub,
                            erea_ind_sub,
                            epro_ind_sub,
                            ebar_ind_sub,
                        ) = self._get_free_energies(dGevb_sub,
                                                    self.coordinate_bins,
                                                    fitting=False)

                        barrier_unc_sub, free_energy_unc_sub = (
                            self._propagate_free_energy_uncertainty(
                                dGevb_unc_sub, erea_ind_sub, ebar_ind_sub,
                                epro_ind_sub))
                    except Exception as e:
                        print_exception_if_debug()
                        self.ostream.print_warning(
                            f"Error {e} encountered while processing replica "
                            f"{replica} ({label}) of configuration '{name}', skipping"
                        )
                        continue

                    by_replica[replica][label] = {
                        "dGfep": dGfep_sub,
                        "dGfep_unc": dGfep_unc_sub,
                        "analytical": {
                            "EVB": dGevb_sub,
                            "EVB_unc": dGevb_unc_sub,
                            "barrier": barrier_sub,
                            "barrier_unc": barrier_unc_sub,
                            "free_energy": free_energy_sub,
                            "free_energy_unc": free_energy_unc_sub,
                            "min_arg": min_arg_sub,
                            "max_arg": max_arg_sub,
                        },
                    }
        finally:
            self.fep_estimator = saved_estimator

        result.update({'replicas': by_replica})
        self.ostream.flush()

    def _calculate_coordinate_bins(self, Lambda_indices, results, bin_size,
                                   dens_threshold):

        assert_msg_critical('scipy' in sys.modules,
                            'scipy is required for EvbDataProcessing.')

        dE_min = 0
        dE_max = 0
        for result in results.values():
            dens_max = np.array([])
            dE = result['dE']
            Lambda_indices_arr = np.asarray(Lambda_indices)

            # gaussian_kde(xy)(xy) is O(n^2); with replica x direction
            # sampling n can reach the tens of thousands, so cap the frames
            # fed to it. The density estimate only ever feeds a coarse,
            # smoothed/thresholded histogram used to pick the coordinate_bins
            # axis bounds below, not the reported free energies, so a
            # representative subsample is safe. minde/maxde (and thus the
            # dE_bins grid) still use the full dE array so the true observed
            # range is preserved.
            n_frames = dE.shape[0]
            subsampled = n_frames > self.dens_kde_max_samples
            if subsampled:
                rng = np.random.default_rng(self.kde_seed)
                kde_idx = rng.choice(n_frames,
                                     size=self.dens_kde_max_samples,
                                     replace=False)
                self.ostream.print_info(
                    f"Subsampling {n_frames} frames to {self.dens_kde_max_samples} "
                    "for the coordinate-bin density estimate")
                self.ostream.flush()
            else:
                kde_idx = np.arange(n_frames)
            dE_kde = dE[kde_idx]
            Lambda_indices_kde = Lambda_indices_arr[kde_idx]

            xy = np.vstack([dE_kde, Lambda_indices_kde])
            assert_msg_critical(
                np.ptp(dE_kde) > 0,
                'EvbDataProcessing: the energy gap dE is constant, so no ' +
                'reaction coordinate can be defined. Check that the two ' +
                'states differ and that alpha and H12 are set.')
            dens = scipy.stats.gaussian_kde(xy)(xy)
            dens = dens / np.max(dens)
            result.update({"dE_dens": dens})
            if subsampled:
                # Only stored when a genuine subsample was taken, so results
                # without subsampling keep exactly the same keys as before
                # (plot_dE_density falls back to the full dE/Lambda_indices
                # when these are absent).
                result.update({
                    "dE_dens_sample":
                    dE_kde,
                    "dE_dens_sample_lambda_indices":
                    Lambda_indices_kde,
                })

            minde = np.min(dE)
            maxde = np.max(dE)
            assert_msg_critical(
                maxde - minde >= 2 * bin_size,
                'EvbDataProcessing: the energy gap spans only ' +
                f'{maxde - minde:.2f} kJ/mol, less than the two bins of ' +
                f'{bin_size} kJ/mol a reaction coordinate needs. Check that ' +
                'the two states differ, and lower bin_size.')
            # One entry per 2 kJ/mol of observed gap; the assert above keeps
            # this at two or more.
            steps = int((maxde - minde) // 2)
            dE_bins = np.linspace(np.min(dE), np.max(dE), steps)
            bin_inds = np.digitize(dE_kde, dE_bins)
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
            # dens_max holds one entry per 2 kJ/mol of observed dE range, so
            # a narrow gap distribution can be shorter than dens_max_window.
            # The running maximum below then gets a negative length. Shrink
            # the window to what the range allows instead.
            window_size = min(self.dens_max_window, max(1, len(dens_max) - 2))
            window_size = min(window_size, len(dens_max))
            if window_size < self.dens_max_window:
                self.ostream.print_warning(
                    f"The dE range only holds {len(dens_max)} density bins, "
                    f"so the running-maximum window is reduced from "
                    f"{self.dens_max_window} to {window_size}.")
                self.ostream.flush()
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

    def _calculate_fg_profiles(self, result):
        assert len(self.coordinate_bins) > 0, "Coordinate bins not set"
        assert len(self.Lambda) > 0, "Lambda not set"

        E1_fg = result["E1_fg"]
        E2_fg = result["E2_fg"]

        names = result.get("E1_fg_names")
        names = names.split(',') if names is not None else [
            fg.name for fg in EvbForceGroup
        ]
        dGfep_fg = []
        dGevb_fg = []
        for i, name in enumerate(names):
            E1 = E1_fg[i]
            E2 = E2_fg[i]
            E2_shifted, V, dE, Eg = self._calculate_Eg_V_dE(
                E1, E2, self.alpha, self.H12)
            dGfep, _ = self._calculate_dGfep(E1, E2_shifted, result["Temp_set"])
            dGevb, shift, fepxi = self._dGevb_analytical(
                dGfep, self.Lambda, self.H12, self.coordinate_bins)
            dGfep_fg.append(dGfep)
            dGevb_fg.append(dGevb)

        result.update({"dGfep_fg_names": names})
        result.update({"dGfep_fg": np.array(dGfep_fg)})
        result.update({"dGevb_fg": np.array(dGevb_fg)})

    @staticmethod
    def print_results(results, ostream=None):
        if ostream == None:
            ostream = OutputStream(sys.stdout)

        ostream.print_info(
            f"{'Analytical':<30} {'Barrier (kJ/mol)':>25} {'Free Energy (kJ/mol)':>25}"
        )
        for name, result in results["configuration_results"].items():
            if "analytical" in result.keys():
                analytical = result["analytical"]
                barrier_str = f"{analytical['barrier']:.2f}"
                if "barrier_unc" in analytical:
                    barrier_str += f" +/- {analytical['barrier_unc']:.2f}"
                free_energy_str = f"{analytical['free_energy']:.2f}"
                if "free_energy_unc" in analytical:
                    free_energy_str += f" +/- {analytical['free_energy_unc']:.2f}"
                ostream.print_info(
                    f"{name:<30} {barrier_str:>25} {free_energy_str:>25}")

        ostream.flush()

    @staticmethod
    def print_uncertainties(results, ostream=None):
        """Print statistical uncertainty, per-replica/direction uncertainty, and hysteresis diagnostics.

        The pooled uncertainties are BAR/MBAR statistical errors, which treat
        every frame as independent. Correlated frames make them too small by
        about sqrt(g), with g the statistical inefficiency, so they are not a
        complete error bar. The replica average (mean +/- SEM across
        replicas) also captures correlation within a replica and is the
        better error bar when there are several replicas.

        The statistical-uncertainty section only needs the pooled dGfep_unc /
        analytical uncertainty, always present once EvbDataProcessing.compute
        has run. The per-replica and hysteresis sections additionally
        require EvbDataProcessing.calculate_per_replica = True (the default)
        and Replica_frame/Direction_frame to have been loaded (see
        EvbDriver._load_output_files) -- if neither is available this prints
        a note and returns after the statistical-uncertainty section.

        Args:
            results (dict): EVB results, as produced by EvbDataProcessing.compute.
            ostream (OutputStream): output stream to print to.
        """
        if ostream == None:
            ostream = OutputStream(sys.stdout)

        ostream.print_header(
            "Statistical uncertainty (pooled, assumes uncorrelated frames)")
        ostream.print_info(
            "Too small for correlated frames; see the replica average below.")
        ostream.print_info(
            f"{'Configuration':<20} {'Barrier (kJ/mol)':>25} {'Free Energy (kJ/mol)':>25} {'dGfep(l=1) unc':>18}"
        )
        for name, result in results["configuration_results"].items():
            if "analytical" not in result:
                continue
            analytical = result["analytical"]
            barrier_str = f"{analytical['barrier']:.2f} +/- {analytical.get('barrier_unc', float('nan')):.2f}"
            free_energy_str = f"{analytical['free_energy']:.2f} +/- {analytical.get('free_energy_unc', float('nan')):.2f}"
            dGfep_unc = result.get("dGfep_unc")
            endpoint_str = (f"{dGfep_unc[-1]:.3f}" if dGfep_unc is not None
                            and len(dGfep_unc) > 0 else "n/a")
            ostream.print_info(
                f"{name:<20} {barrier_str:>25} {free_energy_str:>25} {endpoint_str:>18}"
            )
        ostream.print_blank()

        has_per_replica = any(
            'replicas' in result
            for result in results["configuration_results"].values())
        if not has_per_replica:
            ostream.print_info(
                "No per-replica data available (set EvbDataProcessing.calculate_per_replica "
                "= True and ensure Replica_frame/Direction_frame were loaded).")
            ostream.flush()
            return

        has_replica_average = any(
            "replica_average" in result
            for result in results["configuration_results"].values())
        ostream.print_header("Replica average (mean +/- SEM across replicas)")
        if not has_replica_average:
            ostream.print_info(
                "No replica averages available (requires EvbDataProcessing."
                "calculate_per_replica_pooled = True, the default, so each "
                "replica has a \"both\" forward+backward curve to average).")
        else:
            ostream.print_info(
                f"{'Configuration':<20} {'N':>4} {'Barrier (kJ/mol)':>25} {'Free Energy (kJ/mol)':>25}"
            )
            for name, result in results["configuration_results"].items():
                replica_average = result.get("replica_average")
                if not replica_average:
                    continue
                barrier_str = (f"{replica_average['barrier_mean']:.2f} +/- "
                               f"{replica_average['barrier_sem']:.2f}")
                free_energy_str = (
                    f"{replica_average['free_energy_mean']:.2f} +/- "
                    f"{replica_average['free_energy_sem']:.2f}")
                ostream.print_info(
                    f"{name:<20} {replica_average['n_replicas']:>4} {barrier_str:>25} {free_energy_str:>25}"
                )

        has_hysteresis_summary = any(
            "hysteresis_summary" in result
            for result in results["configuration_results"].values())
        if not has_hysteresis_summary:
            ostream.print_info(
                "No hysteresis summary available (requires per-replica "
                "forward/backward curves, see the Hysteresis section above).")
        else:
            ostream.print_blank()
            ostream.print_header("Hysteresis summary (stats across replicas)")
            ostream.print_info(
                f"{'Configuration':<20} {'Metric':<22} {'Mean':>12} {'Std':>12} "
                f"{'Median':>18} {'Min':>18} {'Max':>18}")
            for name, result in results["configuration_results"].items():
                hysteresis_summary = result.get("hysteresis_summary")
                if not hysteresis_summary:
                    continue
                for metric, stats in hysteresis_summary.items():
                    median_str = f"{stats['median']:.3f} ({stats['median_replica']})"
                    min_str = f"{stats['min']:.3f} ({stats['min_replica']})"
                    max_str = f"{stats['max']:.3f} ({stats['max_replica']})"
                    ostream.print_info(
                        f"{name:<20} {metric:<22} {stats['mean']:12.3f} "
                        f"{stats['std']:12.3f} {median_str:>18} "
                        f"{min_str:>18} {max_str:>18}")
                ostream.print_blank()

        ostream.flush()

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

            # dE_dens (and dens) may come from a subsample of dE (see
            # _calculate_coordinate_bins), so scatter against the matching
            # sample rather than assuming 1:1 alignment with the full dE.
            dE_scatter = result.get("dE_dens_sample", dE)
            if "dE_dens_sample_lambda_indices" in result:
                L_scatter = [
                    Lambda[i] for i in result["dE_dens_sample_lambda_indices"]
                ]
            else:
                L_scatter = L_values

            # dens_max = scipy.signal.savgol_filter(dens_max, 20, 3)

            ax[j, 0].scatter(dE_scatter, L_scatter, c=dens, s=5)
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

            # middle = len(dens_max) // 2
            # start = np.where(dens_max[:middle] == 0)[0][-1]
            # end = np.where(dens_max[middle:] == 0)[0][0] + middle
            start = 0
            end = len(dens_max) - 1

            ax[j, 1].scatter(dE_scatter, dens, s=1)
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
    def plot_results(results,
                     plot_analytical=True,
                     plot_discrete=False,
                     order=None,
                     x_axis_publication=True,
                     figsize=(15, 3.5),
                     rows=1,
                     cols=2):

        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from matplotlib.lines import Line2D

        coordinate_bins = results["coordinate_bins"]
        Lambda = results["Lambda"]

        assert_msg_critical(rows * cols == 2, "2 subplots are required")
        fig, ax = plt.subplots(rows, cols, figsize=figsize)
        bin_indicators = (coordinate_bins[:-1] + coordinate_bins[1:]) / 2
        colors = mcolors.TABLEAU_COLORS

        colorkeys = list(colors.keys())
        legend_lines = []
        legend_labels = []
        if plot_analytical and plot_discrete:
            discrete_linestyle = "--"
        else:
            discrete_linestyle = "-"

        to_plot = []
        names = []
        if order is not None:
            names = order
            for name in order:
                assert name in results['configuration_results'].keys(
                ), f"{name} not in configuration results"
                to_plot.append(results['configuration_results'][name])
        else:
            for name, conf in results['configuration_results'].items():
                names.append(name)
                to_plot.append(conf)
        for i, (name, result) in enumerate(zip(names, to_plot)):

            # Shift both averages by the same amount so that their relative differences stay the same
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
                    EVB = result["analytical"]["EVB"]
                    # The analytical curve is evaluated at the bin edges,
                    # not at the bin centres the discretised profile uses,
                    # so it is drawn against coordinate_bins. Drawing EVB[1:]
                    # against the centres shifted it left by half a bin.
                    analytical_x = coordinate_bins
                    ax[1].plot(
                        analytical_x,
                        EVB,
                        label=f"{name} analytical",
                        color=colors[colorkeys[i]],
                    )
                    max_arg = result['analytical']['max_arg']
                    barrier = result['analytical']['barrier']
                    free_energy = result['analytical']['free_energy']

                    # The points the reported numbers were taken at, as
                    # chosen by _get_free_energies. Results written before
                    # those indices were stored fall back to the ends.
                    erea_ind = result['analytical'].get('reactant_index', 0)
                    epro_ind = result['analytical'].get('product_index',
                                                        len(EVB) - 1)
                    ebar_ind = result['analytical'].get(
                        'barrier_index',
                        int(np.argmax(EVB)) if len(max_arg) else len(EVB) // 2)

                    # EVB is plotted against analytical_x one for one, so
                    # the extremum indices need no shift.
                    zero_ind = erea_ind
                    barrier_ind = ebar_ind
                    free_ind = epro_ind

                    # mark the zero-point
                    ax[1].plot(
                        [
                            analytical_x[max(0, zero_ind - 25)],
                            analytical_x[min(
                                len(analytical_x) - 1, zero_ind + 25)]
                        ],
                        [0, 0],
                        color=colors[colorkeys[i]],
                        linewidth=.5,
                    )
                    ax[1].plot(
                        [analytical_x[zero_ind]] * 2,
                        [-10, +10],
                        color=colors[colorkeys[i]],
                        linewidth=.5,
                    )

                    ax[1].plot(
                        [
                            analytical_x[max(0, barrier_ind - 25)],
                            analytical_x[min(
                                len(analytical_x) - 1, barrier_ind + 25)]
                        ],
                        [barrier, barrier],
                        color=colors[colorkeys[i]],
                        linewidth=.5,
                    )
                    ax[1].plot(
                        [analytical_x[barrier_ind]] * 2,
                        [barrier - 10, barrier + 10],
                        color=colors[colorkeys[i]],
                        linewidth=.5,
                    )
                    barrier_unc = result['analytical'].get('barrier_unc')
                    barrier_label = f"{barrier:.0f}"
                    if barrier_unc is not None:
                        barrier_label += f" ± {barrier_unc:.0f}"
                    ax[1].text(analytical_x[barrier_ind],
                               barrier,
                               barrier_label,
                               ha='left',
                               va='bottom')

                    ax[1].plot(
                        [
                            analytical_x[max(0, free_ind - 25)],
                            analytical_x[min(
                                len(analytical_x) - 1, free_ind + 25)]
                        ],
                        [free_energy, free_energy],
                        color=colors[colorkeys[i]],
                        linewidth=.5,
                    )
                    ax[1].plot(
                        [analytical_x[free_ind]] * 2,
                        [free_energy - 10, free_energy + 10],
                        color=colors[colorkeys[i]],
                        linewidth=.5,
                    )
                    free_energy_unc = result['analytical'].get(
                        'free_energy_unc')
                    free_energy_label = f"{free_energy:.0f}"
                    if free_energy_unc is not None:
                        free_energy_label += f" ± {free_energy_unc:.0f}"
                    ax[1].text(analytical_x[free_ind],
                               free_energy,
                               free_energy_label,
                               ha='left',
                               va='bottom')
                    # #Add barrier
                    # ax[1].plot(

                    # )
                    # #add free energy

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
        ax[0].set_title("FEP profiles", fontsize=12)

        ax[1].set_xlabel(r"$\Delta \mathcal{E}$ (kJ/mol)")
        ax[1].set_ylabel(r"$\Delta G_{EVB}$ (kJ/mol)")
        ax[1].set_title("EVB profiles", fontsize=12)
        if x_axis_publication:
            ax[1].tick_params(labelbottom=False)
            ax[1].set_xlabel("Reaction coordinate")
        # fig.tight_layout()  # Adjust layout
        fig.legend(
            legend_lines,
            legend_labels,
            loc='center left',
            bbox_to_anchor=(0.9, 0.5),
            ncol=1,
        )
        # fig.legend(
        #     legend_lines,
        #     legend_labels,
        #     loc=(0.22, 0.91),
        #     ncol=1,
        # )
        return fig, ax

    @staticmethod
    def plot_fep_hysteresis(results,
                            configuration_name,
                            replicas=None,
                            figsize=(8, 6)):
        """Plot the FEP curve as a continuous forward/backward trace to visualise hysteresis.

        For each replica, the forward sweep (lambda 0 -> 1) is drawn rising
        from the current y-offset; the backward sweep (lambda 1 -> 0) is
        drawn retracing the same lambda range, continuing from the forward
        sweep's end point. If sampling were perfectly reversible the
        backward sweep would land back exactly where the forward sweep
        started; the gap it actually lands on is the hysteresis, and it
        becomes the starting point for the next replica's forward sweep
        (mirroring the physical trajectory, which is seeded from the
        previous sweep's final state -- see EvbFepDriver.run_replicas).

        Args:
            results (dict): EVB results, as produced by EvbDataProcessing.compute.
            configuration_name (str): key into results["configuration_results"].
            replicas (list, optional): which replicas to plot, in order. Defaults to all replicas with both a forward and backward curve, in ascending order.
            figsize (tuple, optional): figure size. Defaults to (12, 4).
        """

        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        from matplotlib.lines import Line2D

        result = results['configuration_results'][configuration_name]
        by_replica = result.get('replicas')
        assert_msg_critical(
            by_replica,
            f"No per-replica data for configuration '{configuration_name}'. "
            "Set EvbDataProcessing.calculate_per_replica = True and ensure "
            "Replica_frame/Direction_frame were loaded (see EvbDriver._load_output_files)."
        )

        Lambda = np.asarray(results['Lambda'])
        colors = list(mcolors.TABLEAU_COLORS.values())

        if replicas is None:
            replicas = sorted(r for r, d in by_replica.items()
                              if 'forward' in d and 'backward' in d)

        fig, ax = plt.subplots(figsize=figsize)

        offset = 0.0
        for i, replica in enumerate(replicas):
            directions = by_replica[replica]
            if 'forward' not in directions or 'backward' not in directions:
                continue

            y_fwd_local = np.asarray(directions['forward']['dGfep'])
            y_bwd_local = np.asarray(directions['backward']['dGfep'])

            y_fwd = offset + y_fwd_local
            ax.plot(Lambda, y_fwd, linewidth=.5, color=colors[0])

            # Backward sweep continues from the forward sweep's end point,
            # then retraces lambda 1 -> 0 using the backward-only estimate
            # (itself 0-anchored at lambda=0) shifted by the same amount so
            # its lambda=1 value matches the forward curve's end point.
            y_bwd = y_fwd[-1] - y_bwd_local[-1] + y_bwd_local
            ax.plot(Lambda[::-1],
                    y_bwd[::-1],
                    linestyle='--',
                    linewidth=.5,
                    color=colors[0])

            # Next replica's forward sweep starts where this one's backward
            # sweep landed -- 0 only if this replica showed no hysteresis.
            offset = y_bwd[0]

        legend_lines = [
            Line2D([0], [0], color='grey', linestyle='-'),
            Line2D([0], [0], color='grey', linestyle='--'),
        ]
        ax.legend(legend_lines, ['forward (l: 0->1)', 'backward (l: 1->0)'])
        ax.set_xlabel(r"$\lambda$")
        ax.set_ylabel(r"$\Delta G_{FEP}$ (kJ/mol)")
        ax.set_title(f"FEP hysteresis: {configuration_name}", fontsize=12)
        return fig, ax

    @staticmethod
    def plot_evb_replicas(results,
                          configuration_names=None,
                          replicas=None,
                          direction='both',
                          figsize=(6, 4)):
        """Overlay each replica's independently-computed EVB profile.

        First-pass visualisation: every replica's EVB curve is plotted as-is
        (each already zeroed at its own reactant-well minimum, per
        EvbDataProcessing._get_free_energies), with no averaging or
        uncertainty band across replicas yet.

        Args:
            results (dict): EVB results, as produced by EvbDataProcessing.compute.
            configuration_names (str | list[str], optional): string or list of strings to use as alternative plotting names.
            replicas (list, optional): which replicas to plot. Defaults to all replicas that have the requested direction.
            direction (str, optional): "both" (pool forward+backward frames per replica), "forward", or "backward". Defaults to "both".
            figsize (tuple, optional): figure size. Defaults to (6, 4).
        """

        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors

        names = []
        if configuration_names is None:
            names = list(results['configuration_results'].keys())
        else:
            if isinstance(configuration_names, str):
                configuration_name = [configuration_names]

            elif isinstance(configuration_names, list):
                names = configuration_names
            else:
                raise ValueError(
                    "configuration_names must be a string or a list of strings")
        plot_colours = list(mcolors.TABLEAU_COLORS.values())
        fig, ax = plt.subplots(figsize=figsize)

        for i, configuration_name in enumerate(names):
            result = results['configuration_results'][configuration_name]
            by_replica = result.get('replicas')
            assert_msg_critical(
                by_replica,
                f"No per-replica data for configuration '{configuration_name}'. "
                "Set EvbDataProcessing.calculate_per_replica = True and ensure "
                "Replica_frame/Direction_frame were loaded (see EvbDriver._load_output_files)."
            )

            coordinate_bins = results['coordinate_bins']

            if replicas is None:
                replicas = sorted(r for r, d in by_replica.items()
                                  if direction in d)

            for j, replica in enumerate(replicas):
                directions = by_replica[replica]
                if direction not in directions:
                    continue
                EVB = directions[direction]['analytical']['EVB']
                if len(names) > 1:
                    if j == 0:
                        label = f"{configuration_name}"
                    else:
                        label = None
                else:
                    label = f"Replica {replica}"
                # The analytical curve is evaluated at the bin edges.
                ax.plot(coordinate_bins,
                        EVB,
                        alpha=0.5,
                        linewidth=1,
                        color=plot_colours[i % len(plot_colours)],
                        label=label)

        ax.set_xlabel(r"$\Delta \mathcal{E}$ (kJ/mol)")
        ax.set_ylabel(r"$\Delta G_{EVB}$ (kJ/mol)")
        ax.set_title(f"EVB profiles per replica", fontsize=12)
        ax.legend()
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
            # discrete curves
            pns = result['discrete']['pns']
            # Both correction terms, as _dGevb_discretised sums them, so the
            # curve drawn here is the one in result['discrete']['EVB'].
            dGcor = result['discrete']['dGcor']
            dGdens = result['discrete'].get('dGdens', np.zeros_like(dGcor))

            pnsfep = pns @ dGfep
            pnscor = np.sum(pns * (dGcor + dGdens).transpose(), axis=1)
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

            # analytical curves
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

    @staticmethod
    def plot_force_contributions(results, dif_to=0):
        from IPython.display import clear_output
        from IPython.display import display
        import matplotlib.pyplot as plt
        import ipywidgets as widgets

        missing_fg = [
            name for name, result in results['configuration_results'].items()
            if 'E1_fg' not in result or 'E2_fg' not in result
        ]

        assert_msg_critical(
            not missing_fg,
            "plot_force_decomp requires force-group contributions "
            "(E1_fg/E2_fg), which are not present for configuration(s): "
            f"{missing_fg}. Run EvbDriver.compute_force_groups() and "
            "EvbDriver.compute_energy_profiles() again before plotting.")

        # Get tableau colors
        tableau_colors = plt.cm.tab10.colors[
            1:]  # skip first color (blue) to match colours with plot_force_decomp

        def rgb_to_hex(rgb):
            """Convert RGB tuple to hex string"""
            return '#{:02x}{:02x}{:02x}'.format(int(rgb[0] * 255),
                                                int(rgb[1] * 255),
                                                int(rgb[2] * 255))

        def lighten_color(rgb, amount=0.3):
            """Make color lighter by moving towards white"""
            r, g, b = rgb
            return (r + (1 - r) * amount, g + (1 - g) * amount,
                    b + (1 - b) * amount)

        def darken_color(rgb, amount=0.3):
            """Make color darker by moving towards black"""
            r, g, b = rgb
            return (r * (1 - amount), g * (1 - amount), b * (1 - amount))

        # Generate variations
        colors = []
        for color in tableau_colors:
            original_hex = rgb_to_hex(color)
            lighter_hex = rgb_to_hex(lighten_color(color, 0.3))
            darker_hex = rgb_to_hex(darken_color(color, 0.3))

            colors.append(lighter_hex)
            colors.append(darker_hex)

        conf_keys = results['configuration_results'].keys()
        ref_result = results['configuration_results'][list(conf_keys)[dif_to]]
        fg_names = ref_result['E1_fg_names']
        max_arg = ref_result['analytical']['max_arg']
        min_arg = ref_result['analytical']['min_arg']
        ref_evb = ref_result['dGevb_fg']

        to_plot = {}
        for name, result in results['configuration_results'].items():
            if name == list(conf_keys)[dif_to]:
                continue
            analytical = result['analytical']
            dif = result['dGevb_fg'] - ref_evb

            # The points the reported barrier and reaction free energy were
            # taken at. min_arg/max_arg hold every extremum that cleared the
            # noise filter, which is neither necessarily two minima and one
            # maximum nor the ones that were used.
            n_points = dif.shape[1]
            erea_ind = analytical.get('reactant_index', 0)
            ebar_ind = analytical.get('barrier_index', n_points // 2)
            epro_ind = analytical.get('product_index', n_points - 1)

            barrier = dif[:, ebar_ind]
            e0 = dif[:, erea_ind]
            e1 = dif[:, epro_ind]
            free_energy = e1 - e0
            to_plot.update({
                f"{name}_barrier_dif": barrier,
                f"{name}_free_energy_dif": free_energy
            })

        to_pop = []
        for i in range(ref_evb.shape[0]):
            keep = False
            for name, data in to_plot.items():
                if abs(data[i]) > 0.01:
                    keep = True
                    break
            if not keep:
                to_pop.append(i)
        for name, data in to_plot.items():
            to_plot[name] = np.delete(data, to_pop)
        tick_labels = np.delete(fg_names.split(','), to_pop)

        fig, ax = plt.subplots(figsize=(20, 6))
        ax.grouped_bar(to_plot, colors=colors, tick_labels=tick_labels)
        ax.legend()

    @staticmethod
    def plot_force_decomp(results, dif_to=0):
        from IPython.display import clear_output
        from IPython.display import display
        import matplotlib.pyplot as plt
        import ipywidgets as widgets

        lam = results['Lambda']
        bins = results['coordinate_bins']
        config_results = results['configuration_results']

        missing_fg = [
            name for name, result in config_results.items()
            if 'E1_fg' not in result or 'E2_fg' not in result
        ]
        assert_msg_critical(
            not missing_fg,
            "plot_force_decomp requires force-group contributions "
            "(E1_fg/E2_fg), which are not present for configuration(s): "
            f"{missing_fg}. Run EvbDriver.compute_force_groups() and "
            "EvbDriver.compute_energy_profiles() again before plotting.")

        def fg_names(result):
            names = result.get('E1_fg_names')
            if names is not None:
                return names.split(',')
            return [fg.name for fg in EvbForceGroup]

        # relevant_column = []
        relevant_fgs = []
        relevant_decomps = []

        dp = EvbDataProcessing()
        dp.Lambda_frame = results["Lambda_frame"]
        dp.Lambda = lam
        dp.Lambda_indices = results["Lambda_indices"]
        dp.alpha = results['alpha']
        dp.H12 = results['H12']

        for result in config_results.values():
            for i, name in enumerate(fg_names(result)):
                if name == EvbForceGroup.FROZEN.name:
                    continue
                if not np.all(result['E1_fg'][i] == 0) or not np.all(
                        result['E2_fg'][i] == 0):
                    if name not in relevant_fgs:
                        relevant_fgs.append(name)
        if dif_to == 0:
            print(
                "Differences are taken with respect to the first configuration in the results dictionary, unless specified with the dif_to argument."
            )
        print(
            f"Static forces stay present over the whole reaction, although their parameters might change"
        )
        print(
            f"Dynamic forces either appear or disappear over the course of the reaction",
            flush=True)
        for result in config_results.values():
            if 'decompositions' in result.keys():
                for name in result['decompositions']['names']:
                    if name not in relevant_decomps:
                        relevant_decomps.append(name)
        fg_checkboxes = [
            widgets.Checkbox(
                True,
                description=fg_name,
            ) for fg_name in relevant_fgs
        ]

        decomp_checkboxes = [
            widgets.Checkbox(
                False,
                description=decomp,
            ) for decomp in relevant_decomps
        ]
        checkboxes = fg_checkboxes + decomp_checkboxes
        plot_output = widgets.Output()

        def update_plot(change=None):
            # x = np.linspace(0, 10, 500)
            with plot_output:
                fig1, ax1 = plt.subplots(1, 3, figsize=(18, 4))
                fig2, ax2 = plt.subplots(1, 1, figsize=(18, 4))

                ax1[0].set_title("FEP")
                ax1[0].set_xlabel(r"$\lambda$")
                ax1[0].set_ylabel(r"$\Delta G_{FEP}$ (kJ/mol)")
                ax1[0].set_xlim(0, 1)
                ax1[1].set_title("EVB")
                ax1[1].set_xlabel(r"$\Delta \mathcal{E}$ (kJ/mol)")
                ax1[1].set_ylabel(r"$\Delta G_{EVB}$ (kJ/mol)")
                ax1[1].set_xlim(bins[0], bins[-1])
                ax1[2].set_title("EVB difference")
                ax1[2].set_xlabel(r"$\Delta \mathcal{E}$ (kJ/mol)")
                ax1[2].set_ylabel(r"$\Delta G_{EVB}$ difference (kJ/mol)")
                ax1[2].set_xlim(bins[0], bins[-1])
                ax2.set_title("E1/E2")
                ax2.set_xlabel(r"$\Delta \mathcal{E}$ (kJ/mol)")
                ax2.set_ylabel(r"$\Delta \mathcal{E}$ (kJ/mol)")
                evbs = []

                for name, result in config_results.items():
                    clear_output(wait=True)

                    # Each result gets its own name->index lookup: two
                    # configurations being compared here may have been built
                    # under different EvbForceGroup numberings (or even
                    # different sets of force groups), so a shared/global
                    # index list would silently misalign one of them.
                    name_to_idx = {n: i for i, n in enumerate(fg_names(result))}
                    fg_to_sum = []
                    for fg_name, cb in zip(relevant_fgs, fg_checkboxes):
                        if cb.value and fg_name in name_to_idx:
                            fg_to_sum.append(name_to_idx[fg_name])
                    E1_fg = np.sum(result['E1_fg'][fg_to_sum], axis=0)
                    E2_fg = np.sum(result['E2_fg'][fg_to_sum], axis=0)

                    if 'decompositions' in result.keys():

                        decomp_to_sum = []
                        for i, cb in enumerate(decomp_checkboxes):
                            if cb.value:
                                decomp_to_sum.append(i)
                        E1_decomp = np.sum(
                            result['decompositions']['E1'][decomp_to_sum],
                            axis=0)
                        E2_decomp = np.sum(
                            result['decompositions']['E2'][decomp_to_sum],
                            axis=0)
                        # todo does the minus sign here work?
                        E1 = E1_fg - E1_decomp
                        E2 = E2_fg - E2_decomp
                    else:
                        E1 = E1_fg
                        E2 = E2_fg

                    E2_shifted, V, dE, Eg = dp._calculate_Eg_V_dE(
                        E1,
                        E2,
                        dp.alpha,
                        dp.H12,
                    )
                    dGfep, _ = dp._calculate_dGfep(E1, E2_shifted,
                                                   result['Temp_set'])
                    dGevb, shift, fepxi = dp._dGevb_analytical(
                        dGfep,
                        lam,
                        dp.H12,
                        bins,
                    )
                    dGevb, _, _, min_arg, max_arg, _, _, _ = dp._get_free_energies(
                        dGevb, bins, fitting=False)
                    evbs.append(dGevb)
                    ax1[0].plot(lam, dGfep)

                    ax1[1].plot(bins, dGevb, label=name)

                    ax2.plot(E1 - np.min(E1), label=f'rea {name}', alpha=0.7)
                    ax2.plot(E2 - np.min(E2), label=f'pro {name}', alpha=0.7)
                for i, (name, evb) in enumerate(zip(config_results.keys(),
                                                    evbs)):

                    dif = evb - evbs[dif_to]
                    ax1[2].plot(bins, dif)

                fig1.legend()
                fig2.legend()
                plt.show()

        recalculate_button = widgets.Button(description="Recalculate")
        recalculate_button.on_click(update_plot)

        # Layout in 5x5 grid
        grid = widgets.GridBox(
            checkboxes,
            layout=widgets.Layout(
                grid_template_columns="repeat(5, 220px)",
                grid_gap="5px",
            ),
        )

        # Initial display
        display(grid, recalculate_button, plot_output)
        update_plot()  # initial empty plot
