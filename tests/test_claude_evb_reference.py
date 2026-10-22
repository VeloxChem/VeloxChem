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
import os
import sys
from pathlib import Path

import numpy as np
from scipy import integrate, optimize, signal, special, stats

from veloxchem.veloxchemlib import (boltzmann_in_hartreeperkelvin,
                                    hartree_in_kjpermol)

# The EVB validation project (plans, reports, tools/) lives outside the
# repository, in a different place on each machine. EVB_VALIDATION_DIR
# overrides the lookup; otherwise the first candidate that has a tools/
# folder is used. A new machine gets a new candidate here.
EVB_VALIDATION_ENV = 'EVB_VALIDATION_DIR'
EVB_VALIDATION_CANDIDATES = (
    # Mac: <code_fall_26>/VeloxChem/tests -> <code_fall_26>/Agent_work
    Path(__file__).resolve().parents[2] / 'Agent_work' /
    'EVB validation plan',
    # Linux box: ~/code/Agent_work
    Path.home() / 'code' / 'Agent_work' / 'EVB validation plan',
)


def evb_validation_dir():
    """The EVB validation project folder on this machine."""

    override = os.environ.get(EVB_VALIDATION_ENV)
    candidates = (Path(override), ) if override else EVB_VALIDATION_CANDIDATES
    for folder in candidates:
        if (folder / 'tools').is_dir():
            return folder
    raise FileNotFoundError(
        'EVB validation project not found; tried ' +
        ', '.join(str(c) for c in candidates) +
        f'. Set {EVB_VALIDATION_ENV} to the folder that holds tools/.')


def add_evb_tools_to_path():
    """Make the project's tools (evb_reorg, evb_reweight) importable."""

    tools = str(evb_validation_dir() / 'tools')
    if tools not in sys.path:
        sys.path.insert(0, tools)

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

    def _V_sampled(self, lam, x):
        # The potential window lam draws from: the mapping potential here,
        # the integration potential in MismatchedModel.
        return self._V(lam, x)

    def _grid(self, lam):
        if lam not in self._cache:
            coarse = np.linspace(self.x_lo, self.x_hi, 20001)
            xm = coarse[np.argmin(self._V_sampled(lam, coarse))]
            h = 1e-4 * (self.x_hi - self.x_lo)
            curv = (self._V_sampled(lam, xm + h) - 2 * self._V_sampled(lam, xm)
                    + self._V_sampled(lam, xm - h)) / h**2
            w = np.sqrt(self.kT / curv)
            x = np.linspace(xm - self.n_width * w, xm + self.n_width * w,
                            self.n_grid)
            V = self._V_sampled(lam, x)
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


class MismatchedModel(GenericModel1D):
    """1-D model whose windows sample v_int(lam, x), not the mapping potential.

    WP-B2. The analysis reads the PES pair h11(x), h22(x) (h22 includes
    delta): E1 = h11, E2 = h22, V_pes(lam) = (1-lam) h11 + lam h22. Window lam
    draws x from exp(-v_int(lam, x)/kT) by inverse CDF (the GenericModel1D
    grid), so the frames come from another ensemble than the one analysed,
    as in the MD code, which samples the per-lambda integration systems
    (digest section 1). v_int is any confining callable that broadcasts over
    lam and x, and may be nonlinear in lam. ep_offset(lam) is an arbitrary
    per-window constant added to Ep only; a correct reweighting cancels it.

    Optional spectator modes y (n_spec of them) stand for the degrees of
    freedom a 1-D gap cannot show: both PES states carry S_pes = k_spec/2
    sum y^2, so they cancel from the gap and from every PES free energy,
    while window lam samples them with force constant k_spec / (1 +
    spec_mismatch(lam)). They add mismatch V_pes - v_int that is
    orthogonal to the gap.

    sample() adds to the base dict: x, Q = sum y^2, E1_int and E2_int (the
    sampled potential at lam = 0 and 1), Ep (at the frame's own lam, plus
    ep_offset) and with cross_energies=True U_int (n_frames x K, the sampled
    potential at every ladder lambda).

    Exact references, by quadrature on one grid of n_quad points over
    [x_lo, x_hi] (spectator parts in closed form):
    - dG_map (= dG_map_quadrature): Delta G of the PES mapping, the answer a
      correct analysis gives.
    - dG_code_limit: the infinite-sample solution of the code's MBAR
      equations (u_k = beta V_pes,k) for frames drawn from v_int.
    - log_gap_density, gap_moments, xi_domain: the PES gap density, summed
      over the roots of xi(x) = xi. dg1, dg2, dg_ground and stationary come
      from _ExactModel. lambda_r is the energy-gap value
      (<xi>_1 - <xi>_0) / 2 of the two PES ends.
    - delta_v_moments and ess_fraction: mean and std of V_pes - v_int in a
      window, and the Kish ESS/N of the weights exp(-beta (V_pes - v_int)).
    """

    def __init__(self, h11, h22, v_int, x_lo, x_hi, delta=0.0, T=T_DEFAULT,
                 ep_offset=None, n_spec=0, k_spec=0.0, spec_mismatch=None,
                 n_grid=200_000, n_width=12.0, n_quad=400001):
        super().__init__(h11, h22, x_lo, x_hi, delta, T, n_grid, n_width)
        self.v_int = v_int
        self.ep_offset = ep_offset
        self.n_spec = int(n_spec)
        self.k_spec = float(k_spec)
        self.spec_mismatch = spec_mismatch
        self.n_dof = 1 + self.n_spec
        self.n_quad = int(n_quad)
        self._blocks = None
        self._pes_cache = {}

    @property
    def lambda_r(self):
        return self.lambda_r_gap_means()

    def shifted(self, alpha):
        """The model whose E2 carries an extra alpha (sampling unchanged)."""

        h22 = self.h22
        return MismatchedModel(self.h11, lambda x: h22(x) + alpha, self.v_int,
                               self.x_lo, self.x_hi, self.delta + alpha,
                               self.T, self.ep_offset, self.n_spec,
                               self.k_spec, self.spec_mismatch, self.n_grid,
                               self.n_width, self.n_quad)

    def _V_sampled(self, lam, x):
        return self.v_int(lam, x)

    def _spec_ratio(self, lam):
        # k_spec / k_spec_int(lam) - 1
        if self.spec_mismatch is None:
            return np.zeros_like(np.asarray(lam, dtype=float))
        return self.spec_mismatch(lam)

    def sampled_energy(self, lam, x, Q):
        """The sampled potential at lam for frames (x, Q), spectators
        included. Broadcasts over lam, x and Q."""

        return (self.v_int(lam, x) + 0.5 * self.k_spec /
                (1.0 + self._spec_ratio(lam)) * Q)

    def _draw(self, lam, z):
        x, cdf = self._grid(lam)
        xs = np.interp(special.ndtr(z[:, 0]), cdf, x)
        Q = (self.kT * (1.0 + self._spec_ratio(lam)) / self.k_spec *
             np.sum(z[:, 1:]**2, axis=1) if self.n_spec else
             np.zeros(len(xs)))
        if self._blocks is not None:
            self._blocks.append((xs, Q))
        if not self.n_spec:
            return self.h11(xs), self.h22(xs)
        S = 0.5 * self.k_spec * Q
        return self.h11(xs) + S, self.h22(xs) + S

    def sample(self,
               lambdas,
               n_per_window,
               rng,
               n_replicas=1,
               phi=0.0,
               bias=None,
               cross_energies=False):
        """The base sample dict plus x, Q, E1_int, E2_int, Ep (and U_int)."""

        # The base loop draws the blocks in the order it concatenates them,
        # so the coordinates collected by _draw line up with the frames.
        self._blocks = []
        try:
            samples = super().sample(lambdas, n_per_window, rng, n_replicas,
                                     phi, bias)
            x = np.concatenate([b[0] for b in self._blocks])
            Q = np.concatenate([b[1] for b in self._blocks])
        finally:
            self._blocks = None
        lam_f = samples['Lambda_frame']
        samples['x'] = x
        samples['Q'] = Q
        samples['E1_int'] = self.sampled_energy(0.0, x, Q)
        samples['E2_int'] = self.sampled_energy(1.0, x, Q)
        Ep = self.sampled_energy(lam_f, x, Q)
        if self.ep_offset is not None:
            Ep = Ep + self.ep_offset(lam_f)
        samples['Ep'] = Ep
        if cross_energies:
            lam = _as_lambdas(lambdas)
            samples['U_int'] = self.sampled_energy(lam[None, :], x[:, None],
                                                   Q[:, None])
        return samples

    # exact references

    def dG_map(self, lambdas):
        return self.dG_map_quadrature(lambdas)

    def _quad_grid(self):
        x = np.linspace(self.x_lo, self.x_hi, self.n_quad)
        wq = np.full(self.n_quad, x[1] - x[0])
        wq[0] *= 0.5
        wq[-1] *= 0.5
        return x, wq

    def _log_density(self, V, wq):
        # ln of exp(-V/kT), normalised with the trapezoid weights wq
        lp = -(V - V.min()) / self.kT
        return lp - special.logsumexp(lp, b=wq)

    def _pes_window(self, lam):
        lam = float(lam)
        if lam not in self._pes_cache:
            x, wq = self._quad_grid()
            self._pes_cache[lam] = self._log_density(self._V(lam, x), wq)
        return self._pes_cache[lam]

    def log_gap_density(self, xi, lam):
        """ln P_lam(xi) of the PES mapping window lam.

        Sum over the monotone pieces of xi(x) of p_lam(x) / |dxi/dx| at the
        root in that piece, each piece interpolated in xi on the grid.
        """

        x, _ = self._quad_grid()
        logp = self._pes_window(lam)
        g = self.h11(x) - self.h22(x)
        dg = np.gradient(g, x)
        with np.errstate(divide='ignore'):
            lq = logp - np.log(np.abs(dg))
        xi = np.asarray(xi, dtype=float)
        out = np.full(xi.shape, -np.inf)
        cuts = np.flatnonzero(np.diff(np.sign(dg)) != 0) + 1
        for seg in np.split(np.arange(len(x)), cuts):
            if len(seg) < 2:
                continue
            gs, ls = g[seg], lq[seg]
            if gs[0] > gs[-1]:
                gs, ls = gs[::-1], ls[::-1]
            inside = (xi >= gs[0]) & (xi <= gs[-1])
            out[inside] = np.logaddexp(out[inside],
                                       np.interp(xi[inside], gs, ls))
        return out

    def gap_moments(self, lam):
        x, wq = self._quad_grid()
        p = np.exp(self._pes_window(lam)) * wq
        g = self.h11(x) - self.h22(x)
        mean = np.sum(p * g)
        return mean, np.sqrt(np.sum(p * (g - mean)**2))

    def xi_domain(self):
        """xi range of the monotone branch of xi(x) holding both PES wells.

        Where xi(x) turns over (an anharmonic gap), the density has an
        integrable edge, and -kT ln P shows spurious extrema next to it;
        stationary() searches inside this branch only.
        """

        x, _ = self._quad_grid()
        g = self.h11(x) - self.h22(x)
        dg = np.gradient(g, x)
        ends = [np.argmax(self._pes_window(l)) for l in (0.0, 1.0)]
        lo, hi = min(ends), max(ends)
        assert np.all(np.sign(dg[lo:hi + 1]) == np.sign(dg[lo])), \
            'xi(x) is not monotone between the two wells'
        sign = np.sign(dg)
        while lo > 0 and sign[lo - 1] == sign[lo]:
            lo -= 1
        while hi < len(x) - 1 and sign[hi + 1] == sign[hi]:
            hi += 1
        return min(g[lo], g[hi]), max(g[lo], g[hi])

    def dG_code_limit(self, lambdas, n_k=None, tol=1e-13, floor=1e-9,
                      max_iter=100):
        """Infinite-sample limit of the code's MBAR dGfep on these frames.

        MBAR with u_k = beta V_pes,k and frames from p_int,k converges to the
        f that solves int p_mix(x) W_k(x) dx = c_k, with p_mix = sum_k c_k
        p_int,k, W_k = c_k exp(f_k - u_k) / sum_j c_j exp(f_j - u_j) and
        c_k = N_k / N (equal unless n_k is given). The spectators cancel
        from W_k, so the x marginals suffice. Solved by Newton, f_0 = 0,
        halving the step until the largest residual decreases. It stops at
        tol, or at the rounding floor of the quadrature sums when that is
        below floor. With v_int = V_pes the solution is dG_map. Returns
        kJ/mol per window, 0 at lambda = 0.
        """

        lam = _as_lambdas(lambdas)
        K = len(lam)
        c = (np.full(K, 1.0 / K) if n_k is None else
             np.asarray(n_k, dtype=float) / np.sum(n_k))
        x, wq = self._quad_grid()
        log_mix = special.logsumexp(
            [np.log(ck) + self._log_density(self._V_sampled(l, x), wq)
             for l, ck in zip(lam, c)], axis=0)
        mw = np.exp(log_mix) * wq
        u = self._V(lam[:, None], x[None, :]) / self.kT

        def weights(f):
            a = np.log(c)[:, None] + f[:, None] - u
            W = np.exp(a - special.logsumexp(a, axis=0))
            return W, W @ mw - c

        f = self.dG_map(lam) / self.kT
        W, r = weights(f)
        for _ in range(max_iter):
            res = np.max(np.abs(r))
            if res < tol:
                break
            J = np.diag(W @ mw) - (W * mw) @ W.T
            step = np.zeros(K)
            step[1:] = -np.linalg.solve(J[1:, 1:], r[1:])
            t = 1.0
            while t > 1e-6:
                W_new, r_new = weights(f + t * step)
                if np.max(np.abs(r_new)) < res:
                    break
                t *= 0.5
            else:
                if res < floor:
                    break
                raise RuntimeError('dG_code_limit: Newton stalled at '
                                   f'residual {res:.3g}')
            f, W, r = f + t * step, W_new, r_new
        else:
            raise RuntimeError('dG_code_limit: Newton did not converge')
        return (f - f[0]) * self.kT

    def _int_window(self, lam):
        x, wq = self._quad_grid()
        lp = self._log_density(self._V_sampled(lam, x), wq)
        dv = self._V(lam, x) - self._V_sampled(lam, x)
        return lp, wq, dv

    def delta_v_moments(self, lam):
        """Mean and std of V_pes,lam - v_int(lam) over the sampled window.

        The spectator part is (kT/2) r chi^2_n with r = spec_mismatch(lam):
        mean n r kT/2, variance n r^2 kT^2 / 2, independent of x.
        """

        lp, wq, dv = self._int_window(lam)
        p = np.exp(lp) * wq
        mean = np.sum(p * dv)
        var = np.sum(p * (dv - mean)**2)
        r = float(self._spec_ratio(lam))
        mean += 0.5 * self.n_spec * r * self.kT
        var += 0.5 * self.n_spec * (r * self.kT)**2
        return mean, np.sqrt(var)

    def ess_fraction(self, lam, cut_kT=30.0):
        """Kish ESS/N of w = exp(-beta dV) in window lam: E[w]^2 / E[w^2].

        The x part is taken over the region where the sampled density is
        within cut_kT of its maximum. Where V_pes is softer than v_int (a
        Morse plateau against a harmonic restraint), E[w^2] over the whole
        line is dominated by regions no finite sample reaches, and can
        diverge. The spectator part is exact, [sqrt(1 + 2r) / (1 + r)]^n,
        finite for r > -1/2.
        """

        lp, wq, dv = self._int_window(lam)
        keep = lp > lp.max() - cut_kT
        lw = -dv[keep] / self.kT
        lb = np.log(wq[keep])
        frac = np.exp(2 * special.logsumexp(lp[keep] + lb + lw)
                      - special.logsumexp(lp[keep] + lb + 2 * lw)
                      - special.logsumexp(lp[keep] + lb))
        r = float(self._spec_ratio(lam))
        if self.n_spec:
            frac *= (np.sqrt(1.0 + 2.0 * r) / (1.0 + r))**self.n_spec
        return float(frac)


def dynamic_bond(lam, eqA, eqB, fcA, fcB, p=3.25, gamma=0.8):
    """r0 and fc of a reactive pair's dynamic harmonic bond at lam.

    Copy of the integration-system formula in reactionsystembuilder.py (the
    dynamic bond block, about line 3018 at 95512a5c3): r0 is interpolated
    linearly and then pulled towards the shorter end, ((r0 - min) / dif)^p
    dif + min, with p = dynamic_bond_tightening; fc is interpolated linearly
    and scaled by 4 (lam - 1/2)^2 (1 - gamma) + gamma, with gamma =
    dynamic_bond_fc_factor.
    """

    lam = np.asarray(lam, dtype=float)
    eq = eqA * (1 - lam) + eqB * lam
    lo, dif = min(eqA, eqB), abs(eqA - eqB)
    if dif > 0:
        eq = ((eq - lo) / dif)**p * dif + lo
    fc = ((fcA * (1 - lam) + fcB * lam)
          * (4 * (lam - 0.5)**2 * (1 - gamma) + gamma))
    return eq, fc


# WP-B2 toy cases, T = 300 K, kJ/mol. 'H12' is the coupling at which the
# exact ground-state profile is read. The sampled potentials are asymmetric
# on purpose: a mirror-symmetric pair gives a code bias b(lam) = b(1 - lam),
# hence none at lam = 1.
# - T0: v_int = V_pes (identity), the EXO set of Model1 (k 145, a 1).
# - T1: PES as T0; the sampled diabats are harmonic, (k1_int, a1_int) for
#   state 1 at +a1_int and (k2_int, a2_int) for state 2 at -a2_int, mixed
#   linearly in lambda.
# - T2: Morse PES diabats in a bond-difference coordinate x, r1 = re + a - x
#   (bond of state 1), r2 = re + a + x (bond of state 2), b = sqrt(k / 2D).
#   The sampled potential is the code's pair of dynamic harmonic bonds
#   (dynamic_bond), with the dissociated ends at r_far1 and r_far2 and
#   force constant fc_far * k, plus lambda * delta. Nonlinear in lambda.
#   D = 2000, not the code's 500: in this coordinate xi(x) turns over about
#   sqrt(D / 32 kT) gap stds from each well (2.5 at D = 500), and the
#   density edge there spoils the exact profile; at 2000 it is ~4 stds out.
# - T3: T2 plus n_spec spectator modes sampled softer than the PES by
#   spec_mid * 4 lam (1 - lam) (k_spec / k_int - 1), so that the mismatch
#   peaks in the middle windows. In 1-D the mismatch std cannot grow much
#   beyond the gap's own fluctuations (stiffer T2 bonds left it at 0.7 kT).
MISMATCH_CASES = {
    'T0': dict(kind='harmonic', k=145.0, a=1.0, delta=-60.0, H12=30.0),
    'T1': dict(kind='harmonic', k=145.0, a=1.0, delta=-60.0, H12=30.0,
               k1_int=120.0, a1_int=0.95, k2_int=150.0, a2_int=1.0),
    'T2': dict(kind='bond', k=145.0, a=1.0, delta=-60.0, H12=15.0,
               D=2000.0, re=1.0, r_far1=2.3, r_far2=3.5, fc_far=0.2, p=3.25,
               gamma=0.8),
    'T3': dict(kind='bond', k=145.0, a=1.0, delta=-60.0, H12=15.0,
               D=2000.0, re=1.0, r_far1=2.3, r_far2=3.5, fc_far=0.2, p=3.25,
               gamma=0.8, n_spec=20, k_spec=100.0, spec_mid=2.5),
}


def make_mismatched(name, T=T_DEFAULT, ep_offset=None, **overrides):
    """MismatchedModel for a MISMATCH_CASES name; overrides replace entries.

    T0 (without overrides) samples with v_int = V_pes exactly, the same
    callable as the mapping.
    """

    p = dict(MISMATCH_CASES[name], **overrides)
    k, a, delta = p['k'], p['a'], p['delta']
    if p['kind'] == 'harmonic':
        h11 = lambda x: 0.5 * k * (x - a)**2
        h22 = lambda x: 0.5 * k * (x + a)**2 + delta
        if 'k1_int' in p:
            k1, a1, k2, a2 = p['k1_int'], p['a1_int'], p['k2_int'], p['a2_int']
            v_int = lambda lam, x: ((1 - lam) * 0.5 * k1 * (x - a1)**2
                                    + lam * (0.5 * k2 * (x + a2)**2 + delta))
        else:
            v_int = lambda lam, x: (1 - lam) * h11(x) + lam * h22(x)
        x_lim = a + 2.5
    else:
        D, re = p['D'], p['re']
        b = np.sqrt(k / (2.0 * D))
        morse = lambda r: D * (1.0 - np.exp(-b * (r - re)))**2
        r1 = lambda x: re + a - x
        r2 = lambda x: re + a + x
        h11 = lambda x: morse(r1(x))
        h22 = lambda x: morse(r2(x)) + delta
        far = p['fc_far'] * k

        def v_int(lam, x):
            eq1, fc1 = dynamic_bond(lam, re, p['r_far1'], k, far, p['p'],
                                    p['gamma'])
            eq2, fc2 = dynamic_bond(lam, p['r_far2'], re, far, k, p['p'],
                                    p['gamma'])
            return (0.5 * fc1 * (r1(x) - eq1)**2 + 0.5 * fc2 *
                    (r2(x) - eq2)**2 + np.asarray(lam) * delta)

        x_lim = a + 2.0
    spec = {}
    if p.get('n_spec'):
        mid = p['spec_mid']
        spec = dict(n_spec=p['n_spec'], k_spec=p['k_spec'],
                    spec_mismatch=lambda lam: mid * 4 * np.asarray(lam) *
                    (1 - np.asarray(lam)))
    return MismatchedModel(h11, h22, v_int, -x_lim, x_lim, delta, T,
                           ep_offset, **spec)


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
    """The EvbDataProcessing.compute input dict (digest section 1 schema).

    Samples from MismatchedModel also carry E1_int, E2_int, Ep and E_m_pes
    (and U_int when drawn), as the driver's loader returns them.
    """

    lambdas = _as_lambdas(lambdas)
    results = {
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
    if 'E1_int' in samples:
        conf = results['configuration_results'][name]
        for key in ('E1_int', 'E2_int', 'Ep', 'U_int'):
            if key in samples:
                conf[key] = samples[key].copy()
        conf['E_m_pes'] = samples['Em'].copy()
    return results


def write_energy_files(folder, samples, lambdas, T=T_DEFAULT):
    """Write Energies.csv, Data_combined.csv and options.json into folder.

    The integration columns hold E1_int and E2_int and Ep is the frame's own
    sampled energy when the samples carry them (MismatchedModel); otherwise
    they repeat the PES energies and Em. Ek, Temp, Vol and Dens are
    placeholders. Floats are
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
        samples['Lambda_frame'], samples['E1'], samples['E2'],
        samples.get('E1_int', samples['E1']),
        samples.get('E2_int', samples['E2']), samples['Em'],
        samples['Replica_frame'],
        samples['Direction_frame']
    ])
    fmt = ['%.17g'] * 6 + ['%d', '%d']
    np.savetxt(E_file, energies, fmt=fmt, delimiter=',',
               header=ENERGIES_HEADER, comments='')

    data = np.column_stack([
        np.arange(n), samples.get('Ep', samples['Em']),
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


# --------------------------------------------------------------------------
# WP-C3: independent energy evaluator for the PES systems.
#
# Written from the textbook definitions in the WP-C3 plan (section 4), not
# from reactionsystembuilder.py. Pure numpy; units nm, kJ/mol, e, degrees in
# the JSON.
# --------------------------------------------------------------------------

_E_CHARGE = 1.602176634e-19  # C (SI, exact)
_AVOGADRO = 6.02214076e23  # 1/mol (SI, exact)
_EPSILON0 = 8.8541878128e-12  # F/m (CODATA 2018)

# Coulomb constant in kJ mol^-1 nm e^-2 and the Faraday constant in C/mol
COULOMB_KE = _E_CHARGE**2 * _AVOGADRO / (4 * np.pi * _EPSILON0) / 1e-9 / 1e3
FARADAY = _E_CHARGE * _AVOGADRO

FF_TERMS = ('bond', 'morse', 'angle', 'torsion', 'improper', 'coul_static',
            'lj_static', 'coul_dynamic', 'lj_dynamic', 'nonbonded_force',
            'efield')


def _ff_key(key):
    """'(0, 1, 2)' -> (0, 1, 2)."""

    return tuple(int(x) for x in key.strip('()').split(','))


def _ff_dihedral(p0, p1, p2, p3):
    """IUPAC dihedral angle (radian) of four points."""

    b1, b2, b3 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    return np.arctan2(np.dot(m1, n2), np.dot(n1, n2))


def _ff_angle(p0, p1, p2):
    """Angle (radian) at p1."""

    u, v = p0 - p1, p2 - p1
    c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return np.arccos(np.clip(c, -1.0, 1.0))


def coulomb_pair(r, qq):
    """k_e qq / r and its first two derivatives."""

    return (COULOMB_KE * qq / r, -COULOMB_KE * qq / r**2,
            2 * COULOMB_KE * qq / r**3)


def lj_pair(r, sigma, epsilon):
    """4 eps [(s/r)^12 - (s/r)^6] and its first two derivatives."""

    a12, a6 = 4 * epsilon * sigma**12, 4 * epsilon * sigma**6
    return (a12 / r**12 - a6 / r**6, -12 * a12 / r**13 + 6 * a6 / r**7,
            156 * a12 / r**14 - 42 * a6 / r**8)


def soft_core(pair, r, rc, *params):
    """The pair potential for r >= rc; below rc its second-order Taylor
    polynomial about rc. rc <= 0 means the plain potential everywhere."""

    if rc <= 0 or r >= rc:
        return pair(r, *params)[0]
    f, f1, f2 = pair(rc, *params)
    return f + f1 * (r - rc) + 0.5 * f2 * (r - rc)**2


class FFEnergyReference:
    """Force-field energy of one state of a reaction path, term by term.

    The PES system of state k should equal the state-k force field with two
    by-design changes: reactive bonds are Morse, and the nonbonded terms of
    reactive pairs may be soft-core. Reactive pairs are the bonds present in
    one state of a step and not the other, over all steps of the path.

    states: the per-state force-field dicts (the content of the
    <state>_ff_data.json files) or paths to them, in path order.
    soft_core_flags: dict with keys coulomb_static, lj_static,
    coulomb_dynamic, lj_dynamic (defaults: dynamic on, static off).
    E_field: field in V/m, or None.
    elements: element symbols per atom; needed when a state has impropers.

    Impropers: the JSON key is (c, a1, a2, a3), centre first, as an OpenMM
    ffxml <Improper>, which is how MMForceFieldGenerator writes them. The
    energy is that of the torsion OpenMM's default (AMBER) ordering builds:
    the dihedral (x, y, c, a3) with {x, y} = {a1, a2}, the lower index first
    for equal elements, otherwise carbon first, otherwise the heavier first.
    """

    def __init__(self,
                 states,
                 coul14_scale=0.833,
                 lj14_scale=0.5,
                 morse_D=500.0,
                 soft_core_flags=None,
                 sc_alpha_q=0.3,
                 sc_sigma_q=1.0,
                 sc_alpha_lj=0.7,
                 sc_power=1 / 6,
                 E_field=None,
                 elements=None):

        self.states = [
            json.loads(Path(s).read_text()) if isinstance(s, (str, Path))
            else s for s in states
        ]
        self.coul14_scale = coul14_scale
        self.lj14_scale = lj14_scale
        self.morse_D = morse_D
        self.soft_core_flags = {
            'coulomb_static': False,
            'lj_static': False,
            'coulomb_dynamic': True,
            'lj_dynamic': True,
        }
        if soft_core_flags:
            self.soft_core_flags.update(soft_core_flags)
        self.sc_alpha_q = sc_alpha_q
        self.sc_sigma_q = sc_sigma_q
        self.sc_alpha_lj = sc_alpha_lj
        self.sc_power = sc_power
        self.E_field = None if E_field is None else np.asarray(E_field, float)
        self.elements = None if elements is None else list(elements)

        self.reactive_pairs = set()
        for left, right in zip(self.states[:-1], self.states[1:]):
            self.reactive_pairs |= (self.bonds_of(left) ^ self.bonds_of(right))
        self._pairs = {}

    @staticmethod
    def bonds_of(state):
        return {tuple(sorted(_ff_key(k))) for k in state['bonds']}

    def n_atoms(self):
        return len(self.states[0]['atoms'])

    def atom(self, k, i):
        return self.states[k]['atoms'][str(i)]

    def graph_distances(self, k):
        """Shortest bond-path length between atoms in state k (inf if none)."""

        n = self.n_atoms()
        dist = np.full((n, n), np.inf)
        neighbors = [set() for _ in range(n)]
        for i, j in self.bonds_of(self.states[k]):
            neighbors[i].add(j)
            neighbors[j].add(i)
        for start in range(n):
            dist[start, start] = 0
            frontier = [start]
            while frontier:
                nxt = []
                for a in frontier:
                    for b in neighbors[a]:
                        if dist[start, b] == np.inf:
                            dist[start, b] = dist[start, a] + 1
                            nxt.append(b)
                frontier = nxt
        return dist

    def pairs(self, k):
        """Every atom pair of state k, classified by state k's bond graph.

        Returns {(i, j): dict(cls, reactive, qq, sigma, epsilon)} with cls
        'excluded' (1-2, 1-3), '1-4' (scaled) or 'full', and the effective
        (scaled) qq and epsilon; excluded pairs carry qq = epsilon = 0.
        """

        if k in self._pairs:
            return self._pairs[k]
        dist = self.graph_distances(k)
        out = {}
        n = self.n_atoms()
        for i in range(n):
            for j in range(i + 1, n):
                ai, aj = self.atom(k, i), self.atom(k, j)
                qq = ai['charge'] * aj['charge']
                sigma = 0.5 * (ai['sigma'] + aj['sigma'])
                epsilon = np.sqrt(ai['epsilon'] * aj['epsilon'])
                if dist[i, j] <= 2:
                    cls, qq, epsilon = 'excluded', 0.0, 0.0
                elif dist[i, j] == 3:
                    cls = '1-4'
                    qq *= self.coul14_scale
                    epsilon *= self.lj14_scale
                else:
                    cls = 'full'
                out[(i, j)] = dict(cls=cls,
                                   reactive=(i, j) in self.reactive_pairs,
                                   qq=qq,
                                   sigma=sigma,
                                   epsilon=epsilon)
        self._pairs[k] = out
        return out

    def improper_atoms(self, k, key):
        """The dihedral atoms of the improper key (c, a1, a2, a3)."""

        assert self.elements is not None, 'impropers need the elements'
        c, x, y, a3 = key
        ex, ey = self.elements[x], self.elements[y]
        if ex == ey:
            swap = x > y
        else:
            swap = ex != 'C' and (ey == 'C' or self.atom(k, x)['mass'] <
                                  self.atom(k, y)['mass'])
        if swap:
            x, y = y, x
        return (x, y, c, a3)

    def coulomb_rc(self, qq):
        return self.sc_alpha_q * (1 + self.sc_sigma_q * qq)

    def lj_rc(self, sigma, epsilon):
        return self.sc_alpha_lj * ((26 / 7) * 4 * epsilon * sigma**6)**self.sc_power

    def nonbonded(self, k, xyz, kind):
        """(Coulomb, LJ) of state k over the static (non-reactive) or
        dynamic (reactive) pairs."""

        xyz = np.asarray(xyz, float)
        coul = lj = 0.0
        for (i, j), p in self.pairs(k).items():
            if p['cls'] == 'excluded' or p['reactive'] != (kind == 'dynamic'):
                continue
            r = np.linalg.norm(xyz[i] - xyz[j])
            rc_q = (self.coulomb_rc(p['qq'])
                    if self.soft_core_flags[f'coulomb_{kind}'] else 0.0)
            rc_lj = (self.lj_rc(p['sigma'], p['epsilon'])
                     if self.soft_core_flags[f'lj_{kind}'] else 0.0)
            coul += soft_core(coulomb_pair, r, rc_q, p['qq'])
            if p['epsilon'] != 0.0:
                lj += soft_core(lj_pair, r, rc_lj, p['sigma'], p['epsilon'])
        return coul, lj

    def energies(self, k, xyz):
        """Energy terms (kJ/mol) of state k at coordinates xyz (nm, n x 3)."""

        xyz = np.asarray(xyz, float)
        state = self.states[k]
        e = dict.fromkeys(FF_TERMS, 0.0)

        for key, b in state['bonds'].items():
            i, j = _ff_key(key)
            r = np.linalg.norm(xyz[i] - xyz[j])
            r0, fc = b['equilibrium'], b['force_constant']
            if tuple(sorted((i, j))) in self.reactive_pairs:
                D = b.get('D', self.morse_D)
                a = np.sqrt(fc / (2 * D))
                e['morse'] += D * (1 - np.exp(-a * (r - r0)))**2
            else:
                e['bond'] += 0.5 * fc * (r - r0)**2

        for key, ang in state['angles'].items():
            i, j, l = _ff_key(key)
            theta = _ff_angle(xyz[i], xyz[j], xyz[l])
            theta0 = np.deg2rad(ang['equilibrium'])
            e['angle'] += 0.5 * ang['force_constant'] * (theta - theta0)**2

        for section, term in (('dihedrals', 'torsion'), ('impropers',
                                                          'improper')):
            for key, d in state[section].items():
                idx = _ff_key(key)
                if section == 'impropers':
                    idx = self.improper_atoms(k, idx)
                phi = _ff_dihedral(*(xyz[a] for a in idx))
                barriers = np.atleast_1d(d['barrier'])
                phases = np.deg2rad(np.atleast_1d(d['phase']))
                # A negative periodicity marks a multi-term entry (Amber)
                ns = np.abs(np.atleast_1d(d['periodicity']))
                e[term] += float(
                    np.sum(barriers * (1 + np.cos(ns * phi - phases))))

        for kind in ('static', 'dynamic'):
            e[f'coul_{kind}'], e[f'lj_{kind}'] = self.nonbonded(k, xyz, kind)

        if self.E_field is not None:
            q = np.array([self.atom(k, i)['charge'] for i in range(len(xyz))])
            # q [e] * E [V/m] * x [nm] -> kJ/mol
            e['efield'] = float(-FARADAY * 1e-12 *
                                np.sum(q * (xyz @ self.E_field)))

        e['total'] = sum(e[t] for t in FF_TERMS)
        return e
