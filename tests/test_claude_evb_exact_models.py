"""Checks of the exact EVB reference generator (validation WP-A1).

The generator, test_claude_evb_reference.py, is the ground truth for the
phase-A checks of EvbDataProcessing, so these tests check it against closed
forms, independent quadrature and its own exact identities. Section numbers
refer to plans/WP-A1-exact-reference-generator.md in the EVB validation
project.
"""

import numpy as np
import pytest
from scipy import integrate, stats

from veloxchem.evbdriver import EvbDriver
from veloxchem.outputstream import OutputStream

from test_claude_evb_reference import (STANDARD_SETS, GenericModel1D,
                                       ar1_series, default_lambdas,
                                       haw_barrier, make_model, make_results,
                                       marcus_barrier, write_energy_files)

EQUAL_CURVATURE = [
    name for name, p in STANDARD_SETS.items() if p.get('k_ratio', 1.0) == 1.0
]


def _sampled_models():
    """Model1 SYM, Model1 CURV2_EXO, Model2 SYM and Model1 CURV05."""

    return {
        'M1_SYM': make_model('SYM'),
        'M1_CURV2_EXO': make_model('CURV2_EXO'),
        'M2_SYM': make_model('SYM', kind=2),
        'M1_CURV05': make_model('CURV05'),
    }


def _generic_curv05():
    """CURV05 as a GenericModel1D, sampled by inverse CDF on a grid."""

    m = make_model('CURV05')
    return GenericModel1D(m.H11, m.H22, -4.0, 4.0, m.delta)


class TestEvbExactModels:

    def test_sampler_moments(self):
        """Per-window gap mean and variance from 1e5 draws (Q1).

        Expected: the exact gap moments in each mapping window of the 31-point
        default ladder (Gaussian x, so the mean and variance of a quadratic
        or linear gap follow in closed form; Model2 gap mean
        -lambda_r (1 - 2 lam) - delta, variance 2 kT lambda_r). Tolerance:
        4 standard errors, using the sample fourth moment for the error of the
        variance, so non-Gaussian gaps are judged correctly.
        """

        lambdas = default_lambdas()
        rng = np.random.default_rng(101)
        n = 100_000
        for label, model in _sampled_models().items():
            s = model.sample(lambdas, n, rng)
            for il, lam in enumerate(lambdas):
                xi = s['xi'][s['Lambda_indices'] == il]
                mean, sd = model.gap_moments(lam)
                se_mean = sd / np.sqrt(n)
                m4 = np.mean((xi - xi.mean())**4)
                se_var = np.sqrt((m4 - xi.var()**2) / n)
                assert abs(xi.mean() - mean) < 4 * se_mean, (label, lam)
                assert abs(xi.var() - sd**2) < 4 * se_var, (label, lam)

    def test_sampler_ks(self):
        """KS test of sampled gaps against the exact gap CDF (Q1).

        Three windows (lam = 0, 0.5, 1) for each of four samplers: Model1 SYM
        (linear gap), Model1 CURV2_EXO (quadratic gap), Model2 SYM (bath with
        E0 = -5e4 kJ/mol), and the generic inverse-CDF sampler for CURV05,
        checked against Model1 CURV05's analytic CDF. 1e5 draws each, fixed
        seed, so the outcome is deterministic. Pass: p > 1e-3 (plan 5.6).
        """

        models = _sampled_models()
        samplers = dict(models)
        samplers['generic_CURV05'] = _generic_curv05()
        exact = dict(models)
        exact['generic_CURV05'] = models['M1_CURV05']
        rng = np.random.default_rng(202)
        for label in ('M1_SYM', 'M1_CURV2_EXO', 'M2_SYM', 'generic_CURV05'):
            for lam in (0.0, 0.5, 1.0):
                s = samplers[label].sample([lam], 100_000, rng)
                cdf = lambda x: exact[label].gap_cdf(x, lam)
                p = stats.kstest(s['xi'], cdf).pvalue
                assert p > 1e-3, (label, lam, p)

    def test_dG_map_closed_vs_quadrature(self):
        """Closed-form dG_map(lam) against independent quadrature (Q2).

        Model1, every standard set: c_lam - c_0 + kT/2 ln(k_lam/k1) against
        -kT ln(Z_lam/Z_0) from scipy quad around a numerically located
        minimum. Model2, equal-curvature sets: lambda_r lam (1-lam) + lam delta
        against the sum of per-mode quadratures. Tolerance 1e-8 kJ/mol: quad
        runs at relative accuracy 1e-13.
        """

        lambdas = default_lambdas()
        for name in STANDARD_SETS:
            m = make_model(name)
            diff = m.dG_map(lambdas) - m.dG_map_quadrature(lambdas)
            assert np.max(np.abs(diff)) < 1e-8, name
        for name in EQUAL_CURVATURE:
            m = make_model(name, kind=2)
            diff = m.dG_map(lambdas) - m.dG_map_quadrature(lambdas)
            assert np.max(np.abs(diff)) < 1e-8, name

    def test_dG_map_equal_curvature_closed_form(self):
        """Model1 dG_map reduces to lambda_r lam (1-lam) + lam delta (5.2).

        For k1 = k2 the general closed form must equal the Marcus mapping
        form exactly; tolerance 1e-10 kJ/mol (rounding only).
        """

        lam = default_lambdas()
        for name in EQUAL_CURVATURE:
            p = STANDARD_SETS[name]
            ref = p['lambda_r'] * lam * (1 - lam) + lam * p['delta']
            assert np.allclose(make_model(name).dG_map(lam), ref, rtol=0,
                               atol=1e-10), name

    def test_profile_identities(self):
        """Exact profile identities and a brute-force check of dg1 (Q3).

        dg2 is computed from the state-2 gap density, independently of dg1,
        so dg2 - dg1 + xi = 0 is a real check (5.3). dg_ground computed by
        the FEP/umbrella route from windows lam = 0.25, 0.5, 0.75 (Rev eq 8)
        must equal dg1 - xi/2 - R/2. Both to 1e-10 kJ/mol on the xi range
        where the densities involved exceed 1e-8 (rounding only).

        dg1 against a Boltzmann-weighted histogram of 1e7 x points: bins of
        5 kJ/mol where the bin-averaged P1 > 1e-8. Tolerance 1e-4 kJ/mol; the
        histogram's bin-edge error is about dx/dx_bin = 2e-5 relative.
        """

        for name, p in STANDARD_SETS.items():
            if name == 'IDENT':
                continue
            m = make_model(name)
            H12 = p['H12']
            lo, hi = m.xi_search_range()
            xi = np.linspace(lo, hi, 4001)
            keep = ((m.log_gap_density(xi, 0.0) > np.log(1e-8))
                    | (m.log_gap_density(xi, 1.0) > np.log(1e-8)))
            d = m.dg2(xi) - m.dg1(xi) + xi
            assert np.max(np.abs(d[keep])) < 1e-10, name
            for lam in (0.25, 0.5, 0.75):
                keep = m.log_gap_density(xi, lam) > np.log(1e-8)
                d = m.dg_ground_window(xi, H12, lam) - m.dg_ground(xi, H12)
                assert np.max(np.abs(d[keep])) < 1e-10, (name, lam)

        for name in ('SYM', 'CURV05', 'CURV2_EXO'):
            m = make_model(name)
            _, x0, _, s = m.window(0.0)
            x, dx = np.linspace(x0 - 12 * s, x0 + 12 * s, 10_000_000,
                                retstep=True)
            w = np.exp(-m.H11(x) / m.kT)
            w *= dx / np.sqrt(2 * np.pi * m.kT / m.k1)
            mean, sd = m.gap_moments(0.0)
            edges = np.arange(mean - 6 * sd, mean + 6 * sd, 5.0)
            mass, _ = np.histogram(m.gap(x), bins=edges, weights=w)
            for i in np.nonzero(mass / 5.0 > 1e-8)[0]:
                exact, _ = integrate.quad(
                    lambda v: np.exp(-m.dg1(np.array([v]))[0] / m.kT),
                    edges[i], edges[i + 1], epsabs=0, epsrel=1e-12)
                diff = -m.kT * (np.log(mass[i]) - np.log(exact))
                assert abs(diff) < 1e-4, (name, edges[i], diff)

    def test_stationary_values(self):
        """Stationary points of dg_ground against closed forms (Q4).

        - Barrier and dG0 against the planning table (plan 5.4), 1e-3 kJ/mol.
        - Symmetric HAW, lambda_r/4 - H12 + H12^2/lambda_r, exact for
          delta = 0 (plan 5.3 derivation): SYM, LOWH, HIGHH, ENZ, SMALLR,
          1e-6 kJ/mol (Brent refinement with xatol 1e-10).
        - Marcus (delta + lambda_r)^2 / (4 lambda_r) for H0, 1e-6.
        - ENDO mirrors EXO: barrier_ENDO = barrier_EXO - dG0_EXO and
          dG0_ENDO = -dG0_EXO, 1e-6.
        - INV has one ground-state minimum (plan 5.4).
        """

        st = {}
        for name, p in STANDARD_SETS.items():
            st[name] = make_model(name).stationary(p['H12'])
            if p['barrier'] is not None:
                assert abs(st[name]['barrier'] - p['barrier']) < 1e-3, name
            if p['dG0'] is not None and name != 'IDENT':
                assert abs(st[name]['dG0'] - p['dG0']) < 1e-3, name
        for name in ('SYM', 'LOWH', 'HIGHH', 'ENZ', 'SMALLR'):
            p = STANDARD_SETS[name]
            ref = haw_barrier(p['lambda_r'], 0.0, p['H12'])
            assert abs(st[name]['barrier'] - ref) < 1e-6, name
        assert abs(st['H0']['barrier'] - marcus_barrier(290.0, 0.0)) < 1e-6
        assert abs(st['ENDO']['barrier']
                   - (st['EXO']['barrier'] - st['EXO']['dG0'])) < 1e-6
        assert abs(st['ENDO']['dG0'] + st['EXO']['dG0']) < 1e-6
        assert st['INV']['n_min'] == 1 and st['INV']['barrier'] is None

    def test_ar1_autocorrelation(self):
        """AR(1) sampler: exact marginal and the requested tau_int (Q5).

        phi = 0.9 gives tau_int = (1 + phi)/(1 - phi) = 19. One window of
        Model1 SYM (gap linear in x), 1e6 frames in the direction-0 chain.
        - Variance within 3 sigma of 2 kT lambda_r, with sigma =
          var sqrt(2 tau2 / N) and tau2 = (1 + phi^2)/(1 - phi^2), the
          integrated autocorrelation time of the squared AR(1) series.
        - tau_int from the FFT autocorrelation summed to lag 200 (truncation
          bias 2 phi^201/(1-phi) < 1e-8) within 10% of 19.
        """

        phi, n = 0.9, 1_000_000
        m = make_model('SYM')
        s = m.sample([0.5], 2 * n, np.random.default_rng(303), phi=phi)
        xi = s['xi'][s['Direction_frame'] == 0]
        assert len(xi) == n
        var = 2 * m.kT * m.lambda_r
        tau2 = (1 + phi**2) / (1 - phi**2)
        assert abs(xi.var() - var) < 3 * var * np.sqrt(2 * tau2 / n)
        d = xi - xi.mean()
        f = np.fft.rfft(d, 2 * n)
        acf = np.fft.irfft(f * np.conj(f))[:201]
        acf /= acf[0]
        tau = 1 + 2 * np.sum(acf[1:])
        assert abs(tau - 19.0) < 1.9, tau

    def test_ar1_series_marginal(self):
        """ar1_series starts in the stationary distribution (5.2).

        Over 2e5 independent chains, frame 0 and frame 5 each have variance 1
        within 4 standard errors (sqrt(2/N)): no equilibration transient.
        """

        z = ar1_series(6, 200_000, 0.9, np.random.default_rng(404))
        for t in (0, 5):
            assert abs(z[t].var() - 1.0) < 4 * np.sqrt(2 / 200_000), t

    def test_model2_matches_model1(self):
        """Model2 and Model1 with equal lambda_r give identical answers (Q6).

        The gap distribution depends only on lambda_r and delta, so gap
        moments in every window, dG_map, and the SYM stationary values must
        agree to 1e-8 kJ/mol. The xi positions of the extrema are compared to
        1e-5 only: rounding of ~1e-14 in dg moves an extremum by
        sqrt(2 delta_f / f'') ~ 1e-6, so a position is ill-conditioned even
        though the value at it is not.
        """

        lambdas = default_lambdas()
        m1, m2 = make_model('SYM'), make_model('SYM', kind=2)
        for lam in lambdas:
            assert np.allclose(m1.gap_moments(lam), m2.gap_moments(lam),
                               rtol=0, atol=1e-8)
        assert np.allclose(m1.dG_map(lambdas), m2.dG_map(lambdas), rtol=0,
                           atol=1e-8)
        s1, s2 = m1.stationary(50.0), m2.stationary(50.0)
        for key in ('barrier', 'dG0'):
            assert abs(s1[key] - s2[key]) < 1e-8, key
        for key in ('xi_ts', 'xi_r', 'xi_p'):
            assert abs(s1[key] - s2[key]) < 1e-5, key

    def test_writer_roundtrip(self, tmp_path):
        """Files read back by EvbDriver._load_output_files unchanged (Q7).

        Model2 SYM (energies near -5e4 kJ/mol), 31 windows, 2 replicas,
        20 frames per window and replica. Energies, Lambda, Lambda_frame,
        Lambda_indices, Replica_frame and Direction_frame must read back
        exactly: the writer uses 17 significant digits.
        """

        m = make_model('SYM', kind=2)
        lambdas = default_lambdas()
        s = m.sample(lambdas, 20, np.random.default_rng(505), n_replicas=2)
        files = write_energy_files(tmp_path, s, lambdas, m.T)
        evb = EvbDriver(ostream=OutputStream(None))
        specific, common = evb._load_output_files(*files)
        ref = make_results(s, lambdas, m.T)
        conf = ref['configuration_results']['synthetic']
        assert np.array_equal(np.asarray(common['Lambda']), ref['Lambda'])
        assert np.array_equal(common['Lambda_frame'], ref['Lambda_frame'])
        assert np.array_equal(np.asarray(common['Lambda_indices']),
                              ref['Lambda_indices'])
        for key in ('E1_pes', 'E2_pes', 'Replica_frame', 'Direction_frame'):
            assert np.array_equal(specific[key], conf[key]), key
        assert specific['Temp_set'] == m.T

    def test_compute_accepts_results(self):
        """EvbDataProcessing.compute runs on make_results output (Q7).

        Model1 SYM, 31 windows, 200 frames each, alpha = 0 and H12 = 50 fixed
        so the fit is skipped. Pass: compute runs and the per-configuration
        dGfep and analytical profile are finite. This checks the interface
        only, not the values (those are A2-A4).
        """

        pytest.importorskip('pymbar')
        from veloxchem.evbdataprocessing import EvbDataProcessing

        m = make_model('SYM')
        lambdas = default_lambdas()
        s = m.sample(lambdas, 200, np.random.default_rng(606))
        dp = EvbDataProcessing(ostream=OutputStream(None))
        dp.alpha = 0.0
        dp.H12 = 50.0
        out = dp.compute(make_results(s, lambdas, m.T), 31.12, 0.0)
        conf = out['configuration_results']['synthetic']
        assert np.all(np.isfinite(conf['dGfep']))
        analytical = conf['analytical']
        finite = [
            np.all(np.isfinite(np.asarray(v, dtype=float)))
            for v in analytical.values()
            if isinstance(v, np.ndarray) and v.dtype.kind == 'f'
        ]
        assert finite and all(finite)
