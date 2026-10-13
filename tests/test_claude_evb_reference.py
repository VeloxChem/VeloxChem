"""Exact references for the EVB analysis checks (validation WP-A1).

This module holds reference code only, no tests. Other test_claude_evb_*
files import it with ``from test_claude_evb_reference import ...``. It does
not import EvbDataProcessing: the references must not be shaped by the code
they test.

Conventions (WP-A1 plan, section 5.1):

- Units kJ/mol, default T = 300 K, kT = kb*T with the kb the code uses.
- Gap xi = E1 - E2' with E2' = E2 + alpha; reactant states sit at negative xi.
- The models' H22 already contains the offset delta. The code's alpha is an
  extra shift on top of it; ``model.shifted(alpha)`` gives the model whose
  exact answers belong to that alpha. All other outputs are for alpha = 0.
- Lambda values are floats rounded to 3 decimals, as the loader matches them.

Profiles (any model with constant H12), with P_l(xi) the gap density in the
mapping window l (l = 0 is state 1, l = 1 is state 2):

- dg1(xi) = -kT ln P_0(xi)
- dg2(xi) = -kT ln P_1(xi) + dG_map(1)   (equals dg1(xi) - xi exactly)
- dg_ground(xi) = dg1(xi) - xi/2 - sqrt(xi^2 + 4 H12^2)/2

The additive constant is the same for all three, so differences between them
are free energies relative to the state-1 partition function.
"""

import json
from pathlib import Path

import numpy as np
from scipy import integrate, optimize, signal, special, stats

from veloxchem.veloxchemlib import (boltzmann_in_hartreeperkelvin,
                                    hartree_in_kjpermol)

KB = boltzmann_in_hartreeperkelvin() * hartree_in_kjpermol()
T_DEFAULT = 300.0

ENERGIES_HEADER = ('Lambda,reactant PES,product PES,reactant integration,'
                   'product integration,Em,replica,direction')
DATA_HEADER = 'step,Ep,Ek,Temp,Vol,Dens'

# Standard parameter sets (WP-A1 plan, section 5.4). Values in kJ/mol.
# 'barrier' and 'dG0' are the planning-session values (scratch numpy, grid of
# 0.005 kJ/mol); None where the set has no such value or it is computed here.
# 'k_ratio' is k2/k1; for k_ratio != 1, lambda_r is 2*k1*a^2 (state 1).
STANDARD_SETS = {
    'SYM': dict(lambda_r=290.0, delta=0.0, H12=50.0, barrier=31.1207,
                dG0=0.0),
    'EXO': dict(lambda_r=290.0, delta=-60.0, H12=30.0, barrier=20.3628,
                dG0=-58.6368),
    'ENDO': dict(lambda_r=290.0, delta=60.0, H12=30.0, barrier=78.9995,
                 dG0=58.6368),
    'LOWH': dict(lambda_r=290.0, delta=0.0, H12=5.0, barrier=67.5862,
                 dG0=0.0),
    'HIGHH': dict(lambda_r=290.0, delta=0.0, H12=100.0, barrier=6.9828,
                  dG0=0.0),
    'H0': dict(lambda_r=290.0, delta=0.0, H12=0.0, barrier=72.5, dG0=0.0),
    'ENZ': dict(lambda_r=92.0, delta=0.0, H12=20.0, barrier=7.3478, dG0=0.0),
    'SMALLR': dict(lambda_r=20.0, delta=0.0, H12=2.0, barrier=3.2, dG0=0.0),
    'INV': dict(lambda_r=290.0, delta=-400.0, H12=5.0, barrier=None,
                dG0=None),
    'IDENT': dict(lambda_r=0.0, delta=0.0, H12=0.0, barrier=None, dG0=0.0),
    'CURV05': dict(lambda_r=290.0, delta=0.0, H12=30.0, k_ratio=0.5,
                   barrier=None, dG0=None),
    'CURV2': dict(lambda_r=290.0, delta=0.0, H12=30.0, k_ratio=2.0,
                  barrier=None, dG0=None),
    'CURV05_EXO': dict(lambda_r=290.0, delta=-60.0, H12=30.0, k_ratio=0.5,
                       barrier=None, dG0=None),
    'CURV2_EXO': dict(lambda_r=290.0, delta=-60.0, H12=30.0, k_ratio=2.0,
                      barrier=None, dG0=None),
}


def default_lambdas():
    """The code's default 31-point ladder (EvbDriver, non-debug)."""

    lam = np.linspace(0, 0.1, 6)
    lam = np.append(lam[:-1], np.linspace(0.1, 0.9, 21))
    lam = np.append(lam[:-1], np.linspace(0.9, 1, 6))
    return np.round(lam, 3)


def _as_lambdas(lambdas):
    return np.round(np.atleast_1d(np.asarray(lambdas, dtype=float)), 3)


def marcus_barrier(lambda_r, dG0):
    """Marcus: (dG0 + lambda_r)^2 / (4 lambda_r) (digest section 2)."""

    return (dG0 + lambda_r)**2 / (4.0 * lambda_r)


def haw_barrier(lambda_r, dG0, H12):
    """Hwang-Aqvist-Warshel adiabatic barrier with constant H12.

    (dG0 + lambda_r)^2 / (4 lambda_r) - H12 + H12^2 / (dG0 + lambda_r),
    digest section 2 (KW2010 eq 10). Exact for dG0 = 0 in the equal-curvature
    model, perturbative otherwise.
    """

    return marcus_barrier(lambda_r, dG0) - H12 + H12**2 / (dG0 + lambda_r)


def ar1_series(n, n_dof, phi, rng):
    """Stationary AR(1) chains of unit variance, shape (n, n_dof).

    z_0 is an exact N(0, 1) draw and z_t = phi z_{t-1} + sqrt(1-phi^2) e_t, so
    every frame is exactly N(0, 1) and tau_int = (1 + phi) / (1 - phi).
    phi = 0 gives independent draws.
    """

    eps = rng.standard_normal((n, n_dof))
    if phi == 0.0:
        return eps
    eps[1:] *= np.sqrt(1.0 - phi**2)
    return signal.lfilter([1.0], [1.0, -phi], eps, axis=0)


class _ExactModel:
    """Sampling and profile logic shared by the exact models.

    Subclasses provide _draw(lam, z) -> (E1, E2), n_dof, and the gap density
    ln P_lam(xi) through log_gap_density, gap_moments and xi_domain.
    """

    def __init__(self, delta, T):
        self.delta = float(delta)
        self.T = float(T)
        self.kT = KB * self.T

    def sample(self,
               lambdas,
               n_per_window,
               rng,
               n_replicas=1,
               phi=0.0,
               bias=None):
        """Exact draws from V_lam = (1-lam) H11 + lam H22 for every window.

        Per window and replica, n_per_window frames are split evenly into
        direction 0 and direction 1 (direction 0 gets the extra frame when
        n_per_window is odd). Each (window, replica, direction) block is an
        independent draw, or an independent stationary AR(1) chain for
        phi > 0. ``bias(lam, replica, direction)`` returns an energy added to
        the recorded E2 of that block; the frames themselves are still drawn
        from the unbiased V_lam.

        Returns a dict of arrays: E1, E2, Em, xi, Lambda_frame,
        Lambda_indices, Replica_frame, Direction_frame.
        """

        lambdas = _as_lambdas(lambdas)
        n_dir = [(n_per_window + 1) // 2, n_per_window // 2]
        out = {key: [] for key in ('E1', 'E2', 'Lambda_frame',
                                   'Lambda_indices', 'Replica_frame',
                                   'Direction_frame')}
        for il, lam in enumerate(lambdas):
            for rep in range(n_replicas):
                for direction in (0, 1):
                    n = n_dir[direction]
                    if n == 0:
                        continue
                    z = ar1_series(n, self.n_dof, phi, rng)
                    E1, E2 = self._draw(lam, z)
                    if bias is not None:
                        E2 = E2 + bias(lam, rep, direction)
                    out['E1'].append(E1)
                    out['E2'].append(E2)
                    out['Lambda_frame'].append(np.full(n, lam))
                    out['Lambda_indices'].append(np.full(n, il))
                    out['Replica_frame'].append(np.full(n, rep))
                    out['Direction_frame'].append(np.full(n, direction))
        samples = {key: np.concatenate(val) for key, val in out.items()}
        for key in ('Lambda_indices', 'Replica_frame', 'Direction_frame'):
            samples[key] = samples[key].astype(int)
        lam_f = samples['Lambda_frame']
        samples['Em'] = (1 - lam_f) * samples['E1'] + lam_f * samples['E2']
        samples['xi'] = samples['E1'] - samples['E2']
        return samples

    # profiles

    @property
    def degenerate(self):
        """True when the gap is a constant (identical states)."""

        return self.lambda_r == 0.0

    def dg1(self, xi):
        """Diabatic free energy of state 1: -kT ln P_0(xi)."""

        return -self.kT * self.log_gap_density(xi, 0.0)

    def dg2(self, xi):
        """Diabatic free energy of state 2, from the state-2 gap density."""

        return (-self.kT * self.log_gap_density(xi, 1.0)
                + self.dG_map(1.0)[0])

    def dg_ground(self, xi, H12):
        """Ground-state (E_g) free energy profile for constant H12."""

        xi = np.asarray(xi, dtype=float)
        R = np.sqrt(xi**2 + 4.0 * H12**2)
        return self.dg1(xi) - 0.5 * xi - 0.5 * R

    def dg_ground_window(self, xi, H12, lam):
        """dg_ground by the FEP/umbrella route from mapping window lam.

        dG_map(lam) - kT ln P_lam(xi) + (E_g - V_lam)(xi), with
        E_g - V_lam = -xi/2 - R/2 + lam xi (Rev eq 8). Exact for every lam.
        """

        xi = np.asarray(xi, dtype=float)
        R = np.sqrt(xi**2 + 4.0 * H12**2)
        return (self.dG_map(lam)[0] - self.kT * self.log_gap_density(xi, lam)
                - 0.5 * xi - 0.5 * R + lam * xi)

    def lambda_r_gap_means(self):
        """Half the difference of the gap means: (<xi>_2 - <xi>_1) / 2."""

        return 0.5 * (self.gap_moments(1.0)[0] - self.gap_moments(0.0)[0])

    def xi_search_range(self, n_sigma=10.0):
        """xi range covering both diabatic wells, inside the gap's domain."""

        m0, s0 = self.gap_moments(0.0)
        m1, s1 = self.gap_moments(1.0)
        lo = min(m0 - n_sigma * s0, m1 - n_sigma * s1)
        hi = max(m0 + n_sigma * s0, m1 + n_sigma * s1)
        dom_lo, dom_hi = self.xi_domain()
        # Stay clear of the integrable 1/sqrt edge of a quadratic gap, where
        # -kT ln P dips logarithmically.
        margin = 0.01 * (hi - lo)
        return max(lo, dom_lo + margin), min(hi, dom_hi - margin)

    def stationary(self, H12, dxi=0.005):
        """Minima and maxima of dg_ground, grid search then Brent refinement.

        Returns a dict: barrier, dG0, xi_ts, xi_r, xi_p, dg_ts, dg_r, dg_p,
        n_min, lambda_r. Barrier and dG0 are None unless there are two or
        more minima (the reactant is the lowest-xi minimum, the product the
        highest-xi one, the TS the highest maximum between them).
        """

        res = dict(barrier=None, dG0=None, xi_ts=None, xi_r=None, xi_p=None,
                   dg_ts=None, dg_r=None, dg_p=None, n_min=0,
                   lambda_r=self.lambda_r)
        if self.degenerate:
            return res
        lo, hi = self.xi_search_range()
        grid = np.arange(lo, hi, dxi)
        dg = self.dg_ground(grid, H12)

        def refine(i, sign):
            f = lambda x: sign * self.dg_ground(np.array([x]), H12)[0]
            opt = optimize.minimize_scalar(f,
                                           bounds=(grid[i - 1], grid[i + 1]),
                                           method='bounded',
                                           options={'xatol': 1e-10})
            return opt.x, sign * opt.fun

        inner = np.arange(1, len(grid) - 1)
        i_min = inner[(dg[inner] < dg[inner - 1])
                      & (dg[inner] <= dg[inner + 1])]
        i_max = inner[(dg[inner] > dg[inner - 1])
                      & (dg[inner] >= dg[inner + 1])]
        minima = [refine(i, 1.0) for i in i_min]
        maxima = [refine(i, -1.0) for i in i_max]
        res['n_min'] = len(minima)
        if len(minima) == 0:
            return res
        res['xi_r'], res['dg_r'] = minima[0]
        if len(minima) < 2:
            return res
        res['xi_p'], res['dg_p'] = minima[-1]
        between = [m for m in maxima if minima[0][0] < m[0] < minima[-1][0]]
        res['xi_ts'], res['dg_ts'] = max(between, key=lambda m: m[1])
        res['barrier'] = res['dg_ts'] - res['dg_r']
        res['dG0'] = res['dg_p'] - res['dg_r']
        return res


class Model1(_ExactModel):
    """One dimension: H11 = k1/2 (x-a)^2, H22 = k2/2 (x+a)^2 + delta.

    V_lam is quadratic in x, so window draws x ~ N(x_lam, kT/k_lam) are exact.
    The gap xi(x) = H11 - H22 is quadratic in x (linear for k1 = k2), and its
    density is a sum over the roots of xi(x) = xi', each weighted by
    1/|dxi/dx|.
    """

    n_dof = 1

    def __init__(self, k1, k2, a, delta, T=T_DEFAULT):
        super().__init__(delta, T)
        self.k1, self.k2, self.a = float(k1), float(k2), float(a)
        # xi(x) = c2 x^2 + c1 x + c0
        self.c2 = 0.5 * (self.k1 - self.k2)
        self.c1 = -self.a * (self.k1 + self.k2)
        self.c0 = 0.5 * (self.k1 - self.k2) * self.a**2 - self.delta

    @property
    def lambda_r(self):
        """2 k1 a^2 (the state-1 reorganization energy)."""

        return 2.0 * self.k1 * self.a**2

    def shifted(self, alpha):
        """The model whose E2 carries an extra alpha (the code's alpha)."""

        return Model1(self.k1, self.k2, self.a, self.delta + alpha, self.T)

    def H11(self, x):
        return 0.5 * self.k1 * (x - self.a)**2

    def H22(self, x):
        return 0.5 * self.k2 * (x + self.a)**2 + self.delta

    def gap(self, x):
        return self.H11(x) - self.H22(x)

    def window(self, lam):
        """k_lam, x_lam, c_lam (minimum of V_lam) and the x std in window."""

        k = (1 - lam) * self.k1 + lam * self.k2
        x0 = self.a * ((1 - lam) * self.k1 - lam * self.k2) / k
        c = ((1 - lam) * self.H11(x0) + lam * self.H22(x0))
        return k, x0, c, np.sqrt(self.kT / k)

    def _draw(self, lam, z):
        _, x0, _, s = self.window(lam)
        x = x0 + s * z[:, 0]
        return self.H11(x), self.H22(x)

    def dG_map(self, lambdas):
        """Closed form c_lam - c_0 + kT/2 ln(k_lam / k1)."""

        lambdas = _as_lambdas(lambdas)
        out = []
        for lam in lambdas:
            k, _, c, _ = self.window(lam)
            k0, _, c0, _ = self.window(0.0)
            out.append(c - c0 + 0.5 * self.kT * np.log(k / k0))
        return np.array(out)

    def _log_Z(self, lam):
        # Shift by the numerically located minimum, not the closed-form c_lam.
        V = lambda x: (1 - lam) * self.H11(x) + lam * self.H22(x)
        k, x0, _, s = self.window(lam)
        opt = optimize.minimize_scalar(V, bracket=(x0 - s, x0 + s),
                                       options={'xtol': 1e-12})
        vmin, xm = opt.fun, opt.x
        f = lambda x: np.exp(-(V(x) - vmin) / self.kT)
        val, _ = integrate.quad(f, xm - 40 * s, xm + 40 * s, points=[xm],
                                epsabs=0, epsrel=1e-13, limit=200)
        return np.log(val) - vmin / self.kT

    def dG_map_quadrature(self, lambdas):
        """-kT ln(Z_lam / Z_0) by numerical quadrature over x."""

        lambdas = _as_lambdas(lambdas)
        lz0 = self._log_Z(0.0)
        return np.array([-self.kT * (self._log_Z(l) - lz0) for l in lambdas])

    def _roots(self, xi):
        """Roots of xi(x) = xi and |dxi/dx| there (nan where none)."""

        xi = np.asarray(xi, dtype=float)
        if self.c2 == 0.0:
            r = (xi - self.c0) / self.c1
            return r, np.full_like(xi, np.nan), np.full_like(xi, abs(self.c1))
        disc = self.c1**2 - 4.0 * self.c2 * (self.c0 - xi)
        sq = np.sqrt(np.where(disc >= 0, disc, np.nan))
        q = -0.5 * (self.c1 + np.copysign(sq, self.c1))
        return q / self.c2, (self.c0 - xi) / q, sq

    def log_gap_density(self, xi, lam):
        """ln P_lam(xi), the exact gap density in mapping window lam."""

        if self.degenerate:
            return np.full_like(np.asarray(xi, dtype=float), np.nan)
        _, x0, _, s = self.window(lam)
        r1, r2, deriv = self._roots(xi)
        lp1 = stats.norm.logpdf(r1, x0, s)
        lp2 = stats.norm.logpdf(r2, x0, s)
        lp = np.logaddexp(np.where(np.isnan(lp1), -np.inf, lp1),
                          np.where(np.isnan(lp2), -np.inf, lp2))
        return np.where(np.isnan(deriv), -np.inf, lp - np.log(deriv))

    def gap_cdf(self, xi, lam):
        """P(gap <= xi) in mapping window lam (exact)."""

        _, x0, _, s = self.window(lam)
        xi = np.asarray(xi, dtype=float)
        if self.c2 == 0.0:
            r = (xi - self.c0) / self.c1
            return (stats.norm.sf(r, x0, s) if self.c1 < 0
                    else stats.norm.cdf(r, x0, s))
        r1, r2, deriv = self._roots(xi)
        rlo, rhi = np.fmin(r1, r2), np.fmax(r1, r2)
        inside = stats.norm.cdf(rhi, x0, s) - stats.norm.cdf(rlo, x0, s)
        inside = np.where(np.isnan(deriv), 0.0, inside)
        return inside if self.c2 > 0 else 1.0 - inside

    def gap_moments(self, lam):
        """Mean and std of the gap in mapping window lam."""

        _, m, _, s = self.window(lam)
        mean = self.c2 * (m**2 + s**2) + self.c1 * m + self.c0
        var = (2 * self.c2 * m + self.c1)**2 * s**2 + 2 * self.c2**2 * s**4
        return mean, np.sqrt(var)

    def xi_domain(self):
        if self.c2 == 0.0:
            return -np.inf, np.inf
        edge = self.c0 - self.c1**2 / (4.0 * self.c2)
        return (edge, np.inf) if self.c2 > 0 else (-np.inf, edge)


class Model2(_ExactModel):
    """Harmonic bath of M modes plus spectator modes and a large offset.

    H11 = sum k_i/2 (x_i - a_i)^2 + S, H22 = sum k_i/2 (x_i + a_i)^2 + delta + S,
    S = sum k_s/2 y_s^2 + E0 is shared by both states. k_i, a_i are drawn once
    from ``seed`` and a_i is scaled so that sum 2 k_i a_i^2 = lambda_r. The gap
    -sum 2 k_i a_i x_i - delta is Gaussian with variance 2 kT lambda_r in every
    window.
    """

    def __init__(self,
                 lambda_r,
                 delta,
                 n_modes=20,
                 seed=0,
                 n_spectator=10,
                 E0=-5.0e4,
                 T=T_DEFAULT):
        super().__init__(delta, T)
        rng = np.random.default_rng(seed)
        self.k = rng.uniform(50.0, 500.0, n_modes)
        a = rng.uniform(0.2, 1.0, n_modes)
        scale = (np.sqrt(lambda_r / np.sum(2 * self.k * a**2))
                 if lambda_r > 0 else 0.0)
        self.a = a * scale
        self.k_s = rng.uniform(50.0, 500.0, n_spectator)
        self.E0 = float(E0)
        self._lambda_r = float(lambda_r)
        self.seed = seed
        self.n_dof = n_modes + n_spectator

    @property
    def lambda_r(self):
        return self._lambda_r

    def shifted(self, alpha):
        return Model2(self._lambda_r, self.delta + alpha, len(self.k),
                      self.seed, len(self.k_s), self.E0, self.T)

    def _draw(self, lam, z):
        M = len(self.k)
        x = (self.a * (1 - 2 * lam)
             + np.sqrt(self.kT / self.k) * z[:, :M])
        y = np.sqrt(self.kT / self.k_s) * z[:, M:]
        S = np.sum(0.5 * self.k_s * y**2, axis=1) + self.E0
        E1 = np.sum(0.5 * self.k * (x - self.a)**2, axis=1) + S
        E2 = np.sum(0.5 * self.k * (x + self.a)**2, axis=1) + self.delta + S
        return E1, E2

    def dG_map(self, lambdas):
        """Closed form lambda_r lam (1-lam) + lam delta."""

        lam = _as_lambdas(lambdas)
        return self._lambda_r * lam * (1 - lam) + lam * self.delta

    def dG_map_quadrature(self, lambdas):
        """Sum over modes of 1-D quadratures, plus lam delta."""

        lam = _as_lambdas(lambdas)
        total = lam * self.delta
        for k, a in zip(self.k, self.a):
            total = total + Model1(k, k, a, 0.0, self.T).dG_map_quadrature(lam)
        return total

    def gap_moments(self, lam):
        mean = -self._lambda_r * (1 - 2 * lam) - self.delta
        return mean, np.sqrt(2 * self.kT * self._lambda_r)

    def log_gap_density(self, xi, lam):
        if self.degenerate:
            return np.full_like(np.asarray(xi, dtype=float), np.nan)
        m, s = self.gap_moments(lam)
        return stats.norm.logpdf(xi, m, s)

    def gap_cdf(self, xi, lam):
        m, s = self.gap_moments(lam)
        return stats.norm.cdf(xi, m, s)

    def xi_domain(self):
        return -np.inf, np.inf


class GenericModel1D(_ExactModel):
    """Any 1-D pair H11(x), H22(x), sampled by inverse CDF on a fine grid.

    Each window's grid spans +-n_width local widths (from the curvature at the
    minimum of V_lam inside [x_lo, x_hi]). It cross-checks the Gaussian
    sampler of Model1, and carries the quartic option for A3. It provides
    sampling and dG_map_quadrature only; exact profiles of a non-Gaussian
    model are added if A3 asks for them.
    """

    n_dof = 1

    def __init__(self, h11, h22, x_lo, x_hi, delta=0.0, T=T_DEFAULT,
                 n_grid=2_000_000, n_width=12.0):
        super().__init__(delta, T)
        self.h11, self.h22 = h11, h22
        self.x_lo, self.x_hi = x_lo, x_hi
        self.n_grid, self.n_width = int(n_grid), n_width
        self._cache = {}

    lambda_r = np.nan

    def _V(self, lam, x):
        return (1 - lam) * self.h11(x) + lam * self.h22(x)

    def _grid(self, lam):
        if lam not in self._cache:
            coarse = np.linspace(self.x_lo, self.x_hi, 20001)
            xm = coarse[np.argmin(self._V(lam, coarse))]
            h = 1e-4 * (self.x_hi - self.x_lo)
            curv = (self._V(lam, xm + h) - 2 * self._V(lam, xm)
                    + self._V(lam, xm - h)) / h**2
            w = np.sqrt(self.kT / curv)
            x = np.linspace(xm - self.n_width * w, xm + self.n_width * w,
                            self.n_grid)
            V = self._V(lam, x)
            p = np.exp(-(V - V.min()) / self.kT)
            cdf = integrate.cumulative_trapezoid(p, x, initial=0.0)
            self._cache[lam] = (x, cdf / cdf[-1])
        return self._cache[lam]

    def _draw(self, lam, z):
        x, cdf = self._grid(lam)
        u = special.ndtr(z[:, 0])
        xs = np.interp(u, cdf, x)
        return self.h11(xs), self.h22(xs)

    def dG_map_quadrature(self, lambdas):
        lam = _as_lambdas(lambdas)

        def log_Z(l):
            x = np.linspace(self.x_lo, self.x_hi, 400001)
            V = self._V(l, x)
            return (np.log(integrate.trapezoid(np.exp(-(V - V.min()) / self.kT),
                                               x)) - V.min() / self.kT)

        lz0 = log_Z(0.0)
        return np.array([-self.kT * (log_Z(l) - lz0) for l in lam])


def quartic_model(k, a, delta, c4, T=T_DEFAULT, **kwargs):
    """H11 = k/2 (x-a)^2 + c4 (x-a)^4, H22 = k/2 (x+a)^2 + delta (for A3)."""

    h11 = lambda x: 0.5 * k * (x - a)**2 + c4 * (x - a)**4
    h22 = lambda x: 0.5 * k * (x + a)**2 + delta
    return GenericModel1D(h11, h22, -4 * a, 4 * a, delta, T, **kwargs)


def make_model(name, kind=1, T=T_DEFAULT, **kwargs):
    """Model1 (a = 1, k1 = lambda_r/2) or Model2 for a STANDARD_SETS name.

    IDENT uses Model1(100, 100, 0, 0). ``kind=2`` requires equal curvature.
    Extra keyword arguments go to Model2.
    """

    p = STANDARD_SETS[name]
    ratio = p.get('k_ratio', 1.0)
    if kind == 2:
        assert ratio == 1.0, 'Model2 has equal curvature only'
        return Model2(p['lambda_r'], p['delta'], T=T, **kwargs)
    if p['lambda_r'] == 0.0:
        return Model1(100.0, 100.0, 0.0, p['delta'], T)
    k1 = 0.5 * p['lambda_r']
    return Model1(k1, ratio * k1, 1.0, p['delta'], T)


def make_results(samples, lambdas, T=T_DEFAULT, name='synthetic'):
    """The EvbDataProcessing.compute input dict (digest section 1 schema)."""

    lambdas = _as_lambdas(lambdas)
    return {
        'Lambda': [float(l) for l in lambdas],
        'Lambda_frame': samples['Lambda_frame'].copy(),
        'Lambda_indices': samples['Lambda_indices'].copy(),
        'configuration_results': {
            name: {
                'E1_pes': samples['E1'].copy(),
                'E2_pes': samples['E2'].copy(),
                'Temp_set': float(T),
                'Replica_frame': samples['Replica_frame'].copy(),
                'Direction_frame': samples['Direction_frame'].copy(),
            }
        },
    }


def write_energy_files(folder, samples, lambdas, T=T_DEFAULT):
    """Write Energies.csv, Data_combined.csv and options.json into folder.

    The integration columns repeat the PES energies and Ep is the sampled
    mapping energy Em; Ek, Temp, Vol and Dens are placeholders. Floats are
    written with 17 significant digits, so they read back exactly.
    Returns the three paths.
    """

    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    lambdas = _as_lambdas(lambdas)
    n = len(samples['E1'])
    E_file = folder / 'Energies.csv'
    data_file = folder / 'Data_combined.csv'
    options_file = folder / 'options.json'

    energies = np.column_stack([
        samples['Lambda_frame'], samples['E1'], samples['E2'], samples['E1'],
        samples['E2'], samples['Em'], samples['Replica_frame'],
        samples['Direction_frame']
    ])
    fmt = ['%.17g'] * 6 + ['%d', '%d']
    np.savetxt(E_file, energies, fmt=fmt, delimiter=',',
               header=ENERGIES_HEADER, comments='')

    data = np.column_stack([
        np.arange(n), samples['Em'],
        np.zeros(n),
        np.full(n, float(T)),
        np.ones(n),
        np.ones(n)
    ])
    np.savetxt(data_file, data, fmt='%.17g', delimiter=',',
               header=DATA_HEADER, comments='')

    with options_file.open('w') as fh:
        json.dump({'Lambda': [float(l) for l in lambdas],
                   'temperature': float(T)}, fh)
    return E_file, data_file, options_file
