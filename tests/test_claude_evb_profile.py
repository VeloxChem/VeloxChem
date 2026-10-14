"""EVB profile reconstruction in EvbDataProcessing against exact dg(xi) (WP-A3).

Dev-only (test_claude_*): must not reach master. Exact references come from
test_claude_evb_reference (WP-A1).

Master identity used throughout (derived in the WP-A3 plan section 4, verified
here to 1e-14 by test_analytical_identity): for any model whose gap density is
P_lam(xi), with dGfep fed exactly,

    analytical(xi) - exact(xi) = kT ln P_lam*(xi),   lam*(xi) = (1 + xi/R)/2

because _dGevb_analytical keeps dG_map(lam*) - 2 H12^2 / R and drops the
-kT ln P_lam*(xi) term of the exact FEP/US expression (digest section 2,
Rev eq 8). Everything else in the noise-free error is linear interpolation of
dG_map on the lambda ladder.

Trap: Model.dG_map() rounds lambda to 3 decimals (_as_lambdas, meant for
ladders). lam*(xi) is not a ladder point, so the predictions here use
_dg_map_exact(), the same closed form without the rounding.
"""

import numpy as np
import pytest
import scipy.signal

from veloxchem.evbdataprocessing import EvbDataProcessing
from veloxchem.outputstream import OutputStream

from test_claude_evb_reference import (KB, T_DEFAULT, STANDARD_SETS,
                                       default_lambdas, make_model,
                                       make_results)

KT = KB * T_DEFAULT
FINE_LADDER = np.round(np.linspace(0, 1, 1001), 6)

# Equal-curvature sets with two minima: the width term cancels between
# stationary points, so the analytical route is exact up to interpolation.
EQUAL_CURVATURE = ['SYM', 'EXO', 'ENDO', 'HIGHH', 'ENZ']


def _processor(H12, alpha=0.0):
    dp = EvbDataProcessing(ostream=OutputStream(None))
    dp.alpha = alpha
    dp.H12 = H12
    dp.calculate_per_replica = False
    return dp


def _dg_map_exact(m, lam):
    """Model1 dG_map closed form without the 3-decimal lambda rounding."""

    k, _, c, _ = m.window(lam)
    k0, _, c0, _ = m.window(0.0)
    return c - c0 + 0.5 * m.kT * np.log(k / k0)


def _lam_star(xi, H12):
    return 0.5 * (1.0 + xi / np.sqrt(xi**2 + 4.0 * H12**2))


def _error_parts(m, H12, lambdas, xi):
    """(interpolation error, width term) of the analytical profile at xi."""

    lam = _lam_star(xi, H12)
    exact_ladder = np.array([_dg_map_exact(m, l) for l in lambdas])
    interp = np.interp(lam, lambdas, exact_ladder) - _dg_map_exact(m, lam)
    width = KT * m.log_gap_density(np.array([xi]), lam)[0]
    return interp, width


def _analytical_noise_free(m, H12, lambdas, xi):
    """_dGevb_analytical fed the exact dG_map on this ladder."""

    dG = np.array([_dg_map_exact(m, l) for l in lambdas])
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    return EvbDataProcessing._dGevb_analytical(dG, lambdas, H12, xi)[0]


@pytest.mark.parametrize('name', ['SYM', 'CURV2'])
def test_analytical_identity(name):
    """analytical(xi) - exact(xi) = kT ln P_lam*(xi), pointwise.

    Expected from the WP-A3 plan section 4 derivation; tolerance 1e-9 kJ/mol
    is pure floating point (both sides are closed forms).
    """

    m = make_model(name)
    H12 = STANDARD_SETS[name]['H12']
    for xi in (-300.0, -100.0, 0.0, 50.0, 200.0):
        lam = _lam_star(xi, H12)
        code = _dg_map_exact(m, lam) - 2 * H12**2 / np.sqrt(xi**2 +
                                                            4 * H12**2)
        exact = m.dg_ground(np.array([xi]), H12)[0]
        width = KT * m.log_gap_density(np.array([xi]), lam)[0]
        assert abs(code - exact - width) < 1e-9


@pytest.mark.parametrize('name', EQUAL_CURVATURE)
def test_analytical_noise_free_exact(name):
    """Equal curvature: the noise-free analytical barrier and dG0 are exact
    up to linear interpolation of dG_map on the lambda ladder.

    Expected values: Model.stationary() (WP-A1 exact table). The residual is
    asserted to equal the predicted interpolation error (plan section 4,
    acceptance Q1) to 1e-6 kJ/mol, and the width term is asserted to cancel
    between stationary points to 1e-9. On the 1001-point ladder the residual
    itself must be below 1e-4 kJ/mol (lambda_r h^2 / 4 = 7e-5 for h = 1e-3).
    """

    m = make_model(name)
    H12 = STANDARD_SETS[name]['H12']
    st = m.stationary(H12)
    for lambdas in (default_lambdas(), FINE_LADDER):
        at = {k: _analytical_noise_free(m, H12, lambdas, st[k])[0]
              for k in ('xi_r', 'xi_ts', 'xi_p')}
        parts = {k: _error_parts(m, H12, lambdas, st[k])
                 for k in ('xi_r', 'xi_ts', 'xi_p')}
        for key, exact, hi in (('barrier', st['barrier'], 'xi_ts'),
                               ('dG0', st['dG0'], 'xi_p')):
            residual = at[hi] - at['xi_r'] - exact
            interp = parts[hi][0] - parts['xi_r'][0]
            width = parts[hi][1] - parts['xi_r'][1]
            assert abs(width) < 1e-9, f'{name} {key} width term {width}'
            assert abs(residual - interp) < 1e-6, f'{name} {key}'
            if len(lambdas) > 100:
                assert abs(residual) < 1e-4, f'{name} {key} {residual}'


@pytest.mark.parametrize('name', ['CURV05', 'CURV2'])
def test_analytical_width_bias_curv(name):
    """Unequal curvature: the analytical barrier carries the width-term bias.

    Characterization. The measured noise-free bias must equal the predicted
    kT[ln P_lam*(xi_ts) - ln P_lam*(xi_r)] from the A1 exact gap density
    (plan acceptance Q1, tolerance 1e-3 kJ/mol). Measured on the 1001-point
    ladder, where interpolation contributes < 1e-3.
    """

    m = make_model(name)
    H12 = STANDARD_SETS[name]['H12']
    st = m.stationary(H12)
    at = {k: _analytical_noise_free(m, H12, FINE_LADDER, st[k])[0]
          for k in ('xi_r', 'xi_ts', 'xi_p')}
    parts = {k: _error_parts(m, H12, FINE_LADDER, st[k])
             for k in ('xi_r', 'xi_ts', 'xi_p')}
    bias = at['xi_ts'] - at['xi_r'] - st['barrier']
    predicted = parts['xi_ts'][1] - parts['xi_r'][1]
    assert abs(predicted) > 1.0, f'{name} should carry a visible bias'
    assert abs(bias - predicted) < 1e-3


def test_analytical_h0_has_no_barrier():
    """H12 = 0: lam*(xi) is 0 or 1 for every xi != 0, so the analytical
    profile is a two-level step and carries no barrier at all.

    Characterization of a documented limit, not a defect: the exact Marcus
    barrier of H0 is 72.5 kJ/mol (WP-A1 table), and the noise-free analytical
    route returns 0 within 1e-9.
    """

    m = make_model('H0')
    st = m.stationary(0.0)
    at = {k: _analytical_noise_free(m, 0.0, default_lambdas(), st[k])[0]
          for k in ('xi_r', 'xi_ts')}
    assert abs(st['barrier'] - 72.5) < 1e-3
    assert abs(at['xi_ts'] - at['xi_r']) < 1e-9


def test_identity_smoothing_off():
    """smooth_window_size = 1, smooth_polynomial_order = 0 is the exact
    identity filter, which is how every 'smoothing off' number in WP-A3 is
    produced.

    A degree w-1 polynomial fits w points exactly, so polyorder = window - 1
    is an identity in exact arithmetic, and WP-A2's addendum used it. It is
    not usable at the default window of 10: the Vandermonde solve is
    ill-conditioned and moves individual points by several kJ/mol, which is
    larger than the smoothing bias it is meant to isolate. Tolerances are
    floating point (1e-9) and the measured conditioning error (> 1 kJ/mol).
    """

    rng = np.random.default_rng(0)
    worst_identity = 0.0
    worst_polyorder = 0.0
    for _ in range(20):
        y = rng.normal(size=120) * 10.0 + np.linspace(0, 50, 120)**2 / 100.0
        off = scipy.signal.savgol_filter(y, 1, 0)
        naive = scipy.signal.savgol_filter(y, 10, 9)
        worst_identity = max(worst_identity, np.max(np.abs(off - y)))
        worst_polyorder = max(worst_polyorder, np.max(np.abs(naive - y)))
    assert worst_identity == 0.0
    assert worst_polyorder > 1.0
    dp = _processor(50.0)
    on = scipy.signal.savgol_filter(y, dp.smooth_window_size,
                                    dp.smooth_polynomial_order)
    assert np.max(np.abs(on - y)) > 1.0


def test_middle_maximum_selection():
    """With several maxima the barrier must be the highest one between the
    two minima, not the positional-middle candidate.

    Constructed profile with maxima of 10, 30 and 12 kJ/mol: the barrier is
    the 30 by the definition of a barrier (digest section 2). The old rule
    took max_arg[len//2], which here is the 12.

    The barrier is measured from the well next to the transition state, the
    dip at xi = 26 that sits 0.29 kJ/mol above zero, not from the lower end
    point of the range, so the expected value is 30.08 - 0.29 = 29.8.
    """

    xi = np.arange(0.0, 100.0, 1.0)
    y = np.zeros_like(xi)
    for centre, height in ((15.0, 10.0), (40.0, 30.0), (75.0, 12.0)):
        y += height * np.exp(-0.5 * ((xi - centre) / 4.0)**2)
    # dips at the two ends of the sampled range, three interior maxima
    y[0] = y[-1] = -0.5
    dp = _processor(50.0)
    dp.smooth_polynomial_order = dp.smooth_window_size - 1  # no smoothing
    # fitting=False is the reporting path, which is what the selection rule
    # is about; the fit's own path uses a fixed narrow window (see
    # _refine_extremum) and reads 28.8 on this deliberately peaky profile.
    out = dp._get_free_energies(y, xi, fitting=False)
    barrier, max_arg, ebar_ind = out[1], out[4], out[7]
    assert len(max_arg) == 3
    assert xi[ebar_ind] == pytest.approx(40.0, abs=1.0)
    assert barrier == pytest.approx(29.8, abs=0.2)


def _sampled_profiles(name, n_per_window, seed, H12=None, dp=None):
    """compute() on one exact-model sample.

    Returns the model, H12, the bins, the configuration result, the
    discretised profile as the code builds it, and the same profile with the
    -kT ln P_n(s) density term taken back out (what the code did before the
    WP-A3 fix).
    """

    m = make_model(name)
    H12 = STANDARD_SETS[name]['H12'] if H12 is None else H12
    lambdas = default_lambdas()
    samples = m.sample(lambdas, n_per_window, np.random.default_rng(seed))
    dp = _processor(H12) if dp is None else dp
    out = dp.compute(make_results(samples, lambdas), 0.0, 0.0)
    c = out['configuration_results']['synthetic']
    bins = np.asarray(out['coordinate_bins'], dtype=float)
    code, pns, _, dGdens = dp._dGevb_discretised(c['dGfep'], c['Eg'], c['V'],
                                                 c['dE'], T_DEFAULT)
    without = code - np.sum(pns * dGdens.transpose(), axis=1)
    return m, H12, bins, c, code, without


def test_discretised_vs_standard_form():
    """The discretised profile carries the -kT ln P_n(s) density term of the
    standard FEP/US form (digest section 2, Rev eq 8).

    SYM, N = 1000, seed 0, alpha = 0, H12 = 50, read at the bins holding the
    exact stationary points so that only the formula is tested, not the
    extremum selection. With the term the barrier must match the exact
    31.1207 kJ/mol (WP-A1) within 0.7 kJ/mol = 3 sigma_seed, sigma_seed =
    0.226 kJ/mol over 200 seeds (work/WP-A3 analyze.py Q3). Without it the
    same profile is low by about 0.6 kJ/mol, which is the defect A3-6 fixed
    here, so the test also asserts that dropping the term moves the barrier
    down by at least 0.3 kJ/mol.
    """

    m, H12, bins, c, code, without = _sampled_profiles('SYM', 1000, 0)
    st = m.stationary(H12)
    i_r, i_ts = (int(np.searchsorted(bins, st[k])) for k in ('xi_r', 'xi_ts'))
    with_term = code[i_ts] - code[i_r]
    without_term = without[i_ts] - without[i_r]
    assert abs(with_term - st['barrier']) < 0.7, (
        f'with the term {with_term}, without {without_term}, '
        f'exact {st["barrier"]}')
    assert with_term - without_term > 0.3


def test_discretised_barrier_matches_exact():
    """What the code reports as discrete/barrier, not just the formula.

    Before WP-A3 the discretised profile was read with the unscoped
    extremum rule, so binning noise supplied the extrema and the reported
    barrier was off by tens of kJ/mol (A3-7). SYM and EXO at N = 1000 must
    now land within 1.0 kJ/mol of the exact barrier (WP-A1), which is about
    4 sigma_seed of the reported value.
    """

    for name in ('SYM', 'EXO'):
        m, H12, _, c, _, _ = _sampled_profiles(name, 1000, 0)
        exact = m.stationary(H12)['barrier']
        assert abs(c['discrete']['barrier'] - exact) < 1.0, name
        assert abs(c['analytical']['barrier'] - exact) < 1.0, name


def test_smoothing_does_not_move_the_barrier():
    """Smoothing locates the extrema; it must not set the reported energy.

    The same samples are processed with the default 10-point window and with
    an exact identity filter (window 1, order 0). Measured differences in
    the reported barrier: 0.000 (SYM), 0.058 (ENZ), 0.137 (EXO) and 0.002
    (LOWH at bin_size 2), and 0.23 for LOWH's discretised reaction free
    energy at bin_size 2, where 2 kJ/mol bins leave few frames each and the
    product well wanders. 0.3 kJ/mol bounds all of them. Before WP-A3 the
    default window cost 0.34 kJ/mol on SYM and 28.9 on LOWH (work/WP-A3 Q4),
    because the reported value was read off the filtered curve.

    LOWH at the default bin_size of 10 is deliberately not included: its
    barrier is one bin wide, so which grid point the filter calls the top
    still moves the answer by about 0.5 kJ/mol. That is resolution, not
    smoothing bias, and bin_size 2 removes it.
    """

    for name, bin_size in (('SYM', 10), ('ENZ', 10), ('EXO', 10),
                           ('LOWH', 2)):
        H12 = STANDARD_SETS[name]['H12']
        smoothed = _processor(H12)
        plain = _processor(H12)
        plain.smooth_window_size = 1
        plain.smooth_polynomial_order = 0
        for dp in (smoothed, plain):
            dp.bin_size = bin_size
        _, _, _, c_smooth, _, _ = _sampled_profiles(name, 1000, 0,
                                                    dp=smoothed)
        _, _, _, c_plain, _, _ = _sampled_profiles(name, 1000, 0, dp=plain)
        for route in ('analytical', 'discrete'):
            for key in ('barrier', 'free_energy'):
                assert c_smooth[route][key] == pytest.approx(
                    c_plain[route][key], abs=0.3), f'{name} {route} {key}'


def test_extremum_refinement_recovers_a_sub_grid_peak():
    """_refine_extremum fits a parabola to the raw profile, so a barrier
    whose top falls between two bins is recovered, not truncated.

    A parabola of curvature -0.02 kJ/mol per (kJ/mol)^2 peaking at xi = 5,
    i.e. exactly between two 10 kJ/mol bins, is sampled on the grid: the
    highest sampled point is 0.5 kJ/mol below the true top. The refinement
    must return the true top to 1e-6, since the profile is exactly
    quadratic.
    """

    dp = _processor(50.0)
    xi = np.arange(-200.0, 201.0, 10.0)
    top, curvature = 40.0, -0.02
    profile = top + curvature * (xi - 5.0)**2
    index = int(np.argmax(profile))
    assert profile[index] == pytest.approx(top - 0.5, abs=1e-9)
    value, position, settled = dp._refine_extremum(xi, profile, index, 0,
                                                   len(xi) - 1, True)
    assert value == pytest.approx(top, abs=1e-6)
    assert position == pytest.approx(5.0, abs=1e-6)
    # the index the energy came from, which the caller reports and marks
    assert settled in (index, index + 1)


def test_coordinate_bins_cover_both_minima():
    """coordinate_bins must contain both exact minima with >= 2 bins of
    margin, for every standard set with two minima (plan acceptance Q6).

    Expected positions: Model.stationary() (WP-A1). bin_size is the default
    10 kJ/mol, so the margin is asserted as 20 kJ/mol.
    """

    for name in ['SYM', 'EXO', 'ENDO', 'LOWH', 'HIGHH', 'ENZ']:
        m = make_model(name)
        H12 = STANDARD_SETS[name]['H12']
        lambdas = default_lambdas()
        samples = m.sample(lambdas, 500, np.random.default_rng(0))
        dp = _processor(H12)
        out = dp.compute(make_results(samples, lambdas), 0.0, 0.0)
        bins = np.asarray(out['coordinate_bins'], dtype=float)
        st = m.stationary(H12)
        assert st['xi_r'] - bins[0] >= 2 * dp.bin_size, name
        assert bins[-1] - st['xi_p'] >= 2 * dp.bin_size, name


def test_edge_h0_zero_coupling_no_nan():
    """H12 = 0 with xi = 0 on the coordinate grid must not produce NaN.

    _dGevb_analytical divides by R = sqrt(xi^2 + 4 H12^2), which is 0 there,
    so xi/R is 0/0. The limit is well defined: the -2 H12^2 / R shift is 0
    and lam* is the midpoint 1/2 (the profile is a step between the two
    endpoint values). Savitzky-Golay then spreads the single NaN over its
    whole window, so the reported barrier is NaN (defect A3-2).
    """

    m = make_model('H0')
    lambdas = default_lambdas()
    samples = m.sample(lambdas, 500, np.random.default_rng(0))
    dp = _processor(0.0)
    dp.coordinate_bins = np.arange(-400.0, 401.0, 10.0)
    assert 0.0 in dp.coordinate_bins
    out = dp.compute(make_results(samples, lambdas), 0.0, 0.0)
    analytical = out['configuration_results']['synthetic']['analytical']
    assert not np.any(np.isnan(analytical['EVB']))
    assert not np.isnan(analytical['barrier'])


def test_edge_identical_states_clear_error():
    """Identical states (IDENT: xi == 0 for every frame) must fail with a
    clear VeloxChemError, not a raw scipy LinAlgError.

    gaussian_kde of a dE that is identically 0 is singular
    ('data appears to lie in a lower-dimensional subspace'), which says
    nothing about EVB input (defect A3-3).
    """

    from veloxchem.errorhandler import VeloxChemError

    m = make_model('IDENT')
    lambdas = default_lambdas()
    samples = m.sample(lambdas, 200, np.random.default_rng(0))
    dp = _processor(0.0)
    with pytest.raises(VeloxChemError) as excinfo:
        dp.compute(make_results(samples, lambdas), 0.0, 0.0)
    assert 'dE' in str(excinfo.value) or 'gap' in str(excinfo.value)


@pytest.mark.parametrize('dens_max_window', [50, 100])
def test_edge_short_dE_range(dens_max_window):
    """A gap range shorter than 2 kJ/mol x dens_max_window must still give
    coordinate bins, not ValueError: negative dimensions are not allowed.

    _calculate_coordinate_bins builds np.zeros(len(dens_max) -
    dens_max_window) where len(dens_max) = (max dE - min dE) // 2. SMALLR
    (lambda_r = 20) spans about 110 kJ/mol, i.e. 55 entries, so the default
    window of 50 is already marginal and 100 crashes (defect A3-4).
    """

    m = make_model('SMALLR')
    lambdas = default_lambdas()
    samples = m.sample(lambdas, 500, np.random.default_rng(0))
    dp = _processor(STANDARD_SETS['SMALLR']['H12'])
    dp.dens_max_window = dens_max_window
    out = dp.compute(make_results(samples, lambdas), 0.0, 0.0)
    bins = np.asarray(out['coordinate_bins'], dtype=float)
    assert len(bins) > 2
    st = m.stationary(STANDARD_SETS['SMALLR']['H12'])
    assert bins[0] < st['xi_r'] and bins[-1] > st['xi_p']


def test_edge_inverted_region_reports_no_barrier():
    """In the inverted region the ground state has a single minimum, so the
    code must say there is no barrier instead of reporting noise.

    INV (lambda_r = 290, Delta = -400, H12 = 5) has n_min = 1 (WP-A1) and a
    monotone noise-free analytical profile. Before WP-A3 every seed got a
    barrier (median 10, worst 440 kJ/mol) with no warning (A3-5); now the
    barrier must be NaN and the reaction free energy must still be reported,
    close to Delta = -400 kJ/mol.
    """

    m, H12, _, c, _, _ = _sampled_profiles('INV', 1000, 0)
    assert m.stationary(H12)['n_min'] == 1
    for route in ('analytical', 'discrete'):
        assert np.isnan(c[route]['barrier']), route
        assert c[route]['free_energy'] < -300.0, route


def test_barrier_below_the_noise_is_not_reported():
    """A profile whose only maximum is smaller than extremum_tolerance has
    no barrier worth reporting.

    A flat profile of 60 points with a single 1 kJ/mol bump 60 kJ/mol wide:
    its prominence is the bump height, which is below the default 2.5 kJ/mol
    (about kT at 300 K), so the barrier must come back NaN. Lowering
    extremum_tolerance under the bump makes it reportable again, which is
    what that knob is for.
    """

    dp = _processor(50.0)
    xi = np.arange(-300.0, 300.0, 10.0)
    profile = 1.0 * np.exp(-0.5 * (xi / 60.0)**2)
    barrier = dp._get_free_energies(profile, xi, fitting=False)[1]
    assert np.isnan(barrier)

    dp.extremum_tolerance = 0.5
    barrier = dp._get_free_energies(profile, xi, fitting=False)[1]
    assert np.isfinite(barrier)
