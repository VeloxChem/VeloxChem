"""FEP estimators of EvbDataProcessing against exact Delta G_map (WP-A2).

Dev-only (test_claude_*): must not reach master. Exact references come from
test_claude_evb_reference (WP-A1): for Model1, Delta G_map(lambda) =
lambda_r lambda (1 - lambda) + lambda Delta exactly (grand plan, A1).

The estimators are called directly on a fresh EvbDataProcessing after setting
Lambda, Lambda_frame and Lambda_indices the way compute() does. Energies are
in kJ/mol.
"""

import numpy as np
import pymbar
import pytest

from veloxchem.evbdataprocessing import EvbDataProcessing
from veloxchem.evbdriver import EvbDriver
from veloxchem.outputstream import OutputStream

from test_claude_evb_reference import (KB, T_DEFAULT, STANDARD_SETS, Model1,
                                       default_lambdas, make_model,
                                       make_results, write_energy_files)

BETA = 1.0 / (KB * T_DEFAULT)


def _processor(estimator='mbar'):
    dp = EvbDataProcessing(ostream=OutputStream(None))
    dp.fep_estimator = estimator
    return dp


def _dgfep(estimator, samples, lambdas, alpha=0.0, dp=None):
    """dGfep and its sigma, set up as compute() would."""

    if dp is None:
        dp = _processor(estimator)
    dp.fep_estimator = estimator
    dp.Lambda = [float(l) for l in lambdas]
    dp.Lambda_frame = np.asarray(samples['Lambda_frame'])
    dp.Lambda_indices = np.asarray(samples['Lambda_indices'])
    dG, ddG = dp._calculate_dGfep(np.asarray(samples['E1']),
                                  np.asarray(samples['E2']) + alpha,
                                  T_DEFAULT)
    return np.asarray(dG, dtype=float), np.asarray(ddG, dtype=float)


def _sample_unequal(model, lambdas, counts, rng):
    """One independent block of counts[k] frames per window k."""

    parts = []
    for k, (lam, n) in enumerate(zip(lambdas, counts)):
        part = model.sample([lam], int(n), rng)
        part['Lambda_indices'][:] = k
        parts.append(part)
    return {key: np.concatenate([p[key] for p in parts]) for key in parts[0]}


def _assert_within_sigma(dG, ddG, exact, n_sigma):
    assert dG[0] == 0.0 and ddG[0] == 0.0
    err = dG[1:] - exact[1:]
    assert np.all(np.isfinite(dG)) and np.all(ddG[1:] > 0)
    # PASS only means something if sigma is not huge (grand plan rule)
    assert np.max(ddG) < 1.0
    assert np.all(np.abs(err) <= n_sigma * ddG[1:]), np.max(
        np.abs(err) / ddG[1:])


@pytest.mark.parametrize('estimator', ['mbar', 'bar'])
@pytest.mark.parametrize('name', ['SYM', 'CURV2'])
def test_dgfep_accuracy(estimator, name):
    """dGfep(lambda) agrees with exact Delta G_map(lambda) (Q1).

    31-point default ladder, 1000 independent frames per window, one seed.
    Expected: Model1.dG_map, closed form (A1). Tolerance 4 sigma at every
    lambda: a fixed-seed check against gross errors; calibration over seeds
    is the job of work/WP-A2. Also requires sigma < 1 kJ/mol, so that a
    pass is not bought with a huge error bar.
    """

    model = make_model(name)
    lambdas = default_lambdas()
    samples = model.sample(lambdas, 1000, np.random.default_rng(11))
    dG, ddG = _dgfep(estimator, samples, lambdas)
    _assert_within_sigma(dG, ddG, model.dG_map(lambdas), 4.0)


@pytest.mark.parametrize('estimator', ['mbar', 'bar'])
def test_invariance_common_energy(estimator):
    """A per-frame energy c_n added to E1 and E2 leaves dGfep unchanged (Q7a).

    c_n = E0 + noise, E0 = -5e4 kJ/mol as in Model2. Exact algebra: c_n
    cancels in every u_kn difference (MBAR) and in dE (BAR). Tolerance
    1e-6 kJ/mol (plan; measured ~5e-12).
    """

    lambdas = default_lambdas()
    samples = make_model('SYM').sample(lambdas, 400, np.random.default_rng(2))
    c = -5e4 + 100.0 * np.random.default_rng(3).standard_normal(
        len(samples['E1']))
    shifted = dict(samples, E1=samples['E1'] + c, E2=samples['E2'] + c)
    dG0, _ = _dgfep(estimator, samples, lambdas)
    dG1, _ = _dgfep(estimator, shifted, lambdas)
    assert np.max(np.abs(dG1 - dG0)) <= 1e-6


@pytest.mark.parametrize('estimator', ['mbar', 'bar'])
def test_alpha_shift_exact(estimator):
    """E2 -> E2 + alpha changes dGfep(lambda_k) by exactly lambda_k alpha (Q7b).

    V_k = (1 - lambda_k) E1 + lambda_k (E2 + alpha) = V_k(alpha=0) +
    lambda_k alpha for every frame. Tolerance 1e-6 kJ/mol.
    """

    lambdas = default_lambdas()
    samples = make_model('EXO').sample(lambdas, 400, np.random.default_rng(4))
    alpha = 37.3
    dG0, _ = _dgfep(estimator, samples, lambdas)
    dGa, _ = _dgfep(estimator, samples, lambdas, alpha=alpha)
    assert np.max(np.abs(dGa - dG0 - lambdas * alpha)) <= 1e-6


def test_warm_start_neutral():
    """The MBAR warm-start cache does not change results (Q7c).

    A processor warmed on unrelated data (EXO) must give the same dGfep on
    SYM data as a fresh one. Tolerance 1e-6 kJ/mol.
    """

    lambdas = default_lambdas()
    other = make_model('EXO').sample(lambdas, 400, np.random.default_rng(9))
    samples = make_model('SYM').sample(lambdas, 400, np.random.default_rng(2))
    dp = _processor()
    _dgfep('mbar', other, lambdas, dp=dp)
    assert dp._mbar_f_k_warm_start is not None
    dG_warm, _ = _dgfep('mbar', samples, lambdas, dp=dp)
    dG_cold, _ = _dgfep('mbar', samples, lambdas)
    assert np.max(np.abs(dG_warm - dG_cold)) <= 1e-6


@pytest.mark.parametrize('estimator', ['mbar', 'bar'])
def test_swap_mirror(estimator):
    """Swapping E1 <-> E2 and reversing the ladder mirrors dGfep (Q7d).

    State j of the swapped problem (lambda' = 1 - lambda_{K-1-j}) has the
    same mapping potential as state K-1-j of the original, so
    dG'(j) = dG(K-1-j) - dG(K-1). Tolerance 1e-6 kJ/mol.
    """

    lambdas = default_lambdas()
    samples = make_model('EXO').sample(lambdas, 400, np.random.default_rng(5))
    K = len(lambdas)
    rev = np.round(1.0 - lambdas[::-1], 3)
    idx = K - 1 - samples['Lambda_indices']
    swapped = dict(samples, E1=samples['E2'], E2=samples['E1'],
                   Lambda_indices=idx, Lambda_frame=rev[idx])
    dG, _ = _dgfep(estimator, samples, lambdas)
    dGs, _ = _dgfep(estimator, swapped, rev)
    assert np.max(np.abs(dGs - (dG[::-1] - dG[-1]))) <= 1e-6


@pytest.mark.parametrize('counts', ['equal', 'unequal'])
def test_bar_matches_pymbar_convention(counts):
    """Per-window BAR equals pymbar.bar on correctly signed works.

    pymbar's convention: w_F = u_{k+1} - u_k on window-k frames and
    w_R = u_k - u_{k+1} on window-(k+1) frames, with
    u_{k+1} - u_k = -beta (lambda_{k+1} - lambda_k) dE, dE = E1 - E2'.
    Expected per-window dG and sigma: pymbar.other_estimators.bar(w_F, w_R)
    divided by beta. Tolerance 1e-8 kJ/mol (solver relative tolerance
    1e-12). Equal counts check sigma; unequal counts (50-400 frames) check
    the estimate as well, which depends on the sign convention through the
    ln(N_F/N_R) term.
    """

    model = make_model('SYM')
    lambdas = np.round(np.linspace(0.0, 1.0, 11), 3)
    rng = np.random.default_rng(6)
    if counts == 'equal':
        n = np.full(len(lambdas), 200)
    else:
        n = rng.integers(50, 401, len(lambdas))
    samples = _sample_unequal(model, lambdas, n, rng)
    dG, ddG = _dgfep('bar', samples, lambdas)

    dE = samples['E1'] - samples['E2']
    idx = samples['Lambda_indices']
    for k in range(len(lambdas) - 1):
        dl = lambdas[k + 1] - lambdas[k]
        ref = pymbar.other_estimators.bar(-BETA * dl * dE[idx == k],
                                          BETA * dl * dE[idx == k + 1])
        dg_ref = ref['Delta_f'] / BETA
        ddg_ref = ref['dDelta_f'] / BETA
        assert abs((dG[k + 1] - dG[k]) - dg_ref) <= 1e-8
        assert abs(np.sqrt(ddG[k + 1]**2 - ddG[k]**2) - ddg_ref) <= 1e-8


@pytest.mark.xfail(strict=True,
                   reason='defect A2-5 open (WP-A2); remove when fixed')
def test_bar_cumulative_sigma_calibration():
    """Cumulative BAR sigma at lambda = 1 is calibrated over seeds.

    SYM, uniform 11-point ladder, 200 frames per window, 200 seeds. The
    z = (dG(1) - exact) / sigma(1) must have var z in [0.70, 1.30] (WP-A2
    statistical acceptance, ~99.9% chi^2 interval for 199 d.o.f.). Adjacent
    BAR steps share a window, so they are positively correlated; adding
    their variances in quadrature under-reports. Linear response predicts
    var z = (4K - 6) / (2K - 2) = 1.9 for K = 11; WP-A2 measured 1.6-2.6
    after the sign fix (A2-5). Runs in about 5 s.
    """

    model = make_model('SYM')
    lambdas = np.round(np.linspace(0.0, 1.0, 11), 3)
    exact = model.dG_map(lambdas)[-1]
    dp = _processor('bar')
    z = []
    for seed in range(200):
        samples = model.sample(lambdas, 200, np.random.default_rng(seed))
        dG, ddG = _dgfep('bar', samples, lambdas, dp=dp)
        z.append((dG[-1] - exact) / ddG[-1])
    assert 0.70 <= np.var(z, ddof=1) <= 1.30, np.var(z, ddof=1)


@pytest.mark.parametrize('estimator', ['mbar', 'bar'])
def test_unequal_counts(estimator):
    """Frame counts per window drawn from 50-400 (Q5).

    Same expectation and 4-sigma tolerance as test_dgfep_accuracy, on SYM
    with the default ladder. BAR goes through the object-array path of _bin.
    """

    model = make_model('SYM')
    lambdas = default_lambdas()
    rng = np.random.default_rng(5)
    counts = rng.integers(50, 401, len(lambdas))
    samples = _sample_unequal(model, lambdas, counts, rng)
    dG, ddG = _dgfep(estimator, samples, lambdas)
    _assert_within_sigma(dG, ddG, model.dG_map(lambdas), 4.0)


def test_compute_time_subsample_unequal_counts(tmp_path):
    """compute() runs when time sub-sampling leaves unequal window counts.

    _load_output_files keeps every time_sub_sample-th frame of the whole
    file, so 200 frames per window with time_sub_sample=3 leaves 66 or 67
    per window. compute() must accept that loader output and give a finite
    dGfep within 4 sigma of exact at lambda = 1 (alpha = 0, H12 fixed).
    """

    model = make_model('SYM')
    lambdas = default_lambdas()
    samples = model.sample(lambdas, 200, np.random.default_rng(8))
    files = write_energy_files(tmp_path, samples, lambdas)
    evb = EvbDriver(ostream=OutputStream(None))
    specific, common = evb._load_output_files(*files, time_sub_sample=3)
    counts = np.bincount(common['Lambda_indices'])
    assert counts.min() != counts.max()

    results = dict(common)
    results['configuration_results'] = {'synthetic': specific}
    dp = _processor()
    dp.alpha = 0.0
    dp.H12 = STANDARD_SETS['SYM']['H12']
    out = dp.compute(results, 31.0, 0.0)['configuration_results']['synthetic']
    dG = np.asarray(out['dGfep'])
    ddG = np.asarray(out['dGfep_unc'])
    assert abs(dG[-1] - model.dG_map(lambdas)[-1]) <= 4.0 * ddG[-1]


def _failed_window_ok(dG, ddG, exact):
    """A window is honest if it is non-finite or its sigma covers its error."""

    finite = np.isfinite(dG) & np.isfinite(ddG)
    covered = np.abs(dG - exact) <= 4.0 * ddG
    return ~finite | covered


@pytest.mark.xfail(strict=True,
                   reason='defect A2-3 open (WP-A2); remove when fixed')
def test_bar_disjoint_windows():
    """Disjoint windows must not give a confident finite answer (Q4).

    lambda_r = 4000 kJ/mol, ladder [0, 0.5, 1]: window means a(1 - 2 lambda)
    are 28 widths apart, so adjacent windows share no configurations. The
    method requires that a failed window is not reported as a finite number
    whose sigma does not cover its error (grand plan: silent failure is
    major). Raising is acceptable. Expected exact values: dG_map = [0, 1000,
    0]. Baseline (negated works): BAR returns ~[0, 903, 27] with sigma
    ~2e-8. With pymbar's sign convention (WP-A2 fix): ~[0, 964, -44] with
    sigma ~[0, 3.0, 4.5], still 12 sigma off. pymbar.bar does not raise
    here, so only an overlap check could flag it (escalated, WP-A2 memo).
    """

    model = Model1(2000.0, 2000.0, 1.0, 0.0)
    lambdas = np.array([0.0, 0.5, 1.0])
    samples = model.sample(lambdas, 1000, np.random.default_rng(0))
    try:
        dG, ddG = _dgfep('bar', samples, lambdas)
    except Exception:
        return
    assert np.all(_failed_window_ok(dG, ddG, model.dG_map(lambdas)))


@pytest.mark.xfail(strict=True,
                   reason='defect A2-3 open (WP-A2); remove when fixed')
def test_bar_empty_window():
    """A window with no frames must not become a finite 0 with sigma 0 (Q4).

    Window 5 (lambda = 0.5) of an 11-point ladder has no frames, as in a
    per-replica subset that misses a window. The two BAR steps that touch it
    cannot be estimated, so the result from lambda = 0.5 on must be
    non-finite, raise, or carry a sigma that covers the error. Baseline:
    pymbar.bar raises on the empty array, the except branch in
    _calculate_dGfep_bar sets both steps to 0 with sigma 0 (with a
    print_warning), and dG(0.7) is 30.1 against exact 18.9, sigma 0.2.
    (An empty last window instead raises IndexError outside the try.)
    """

    model = make_model('EXO')
    lambdas = np.round(np.linspace(0.0, 1.0, 11), 3)
    samples = model.sample(lambdas, 200, np.random.default_rng(1))
    keep = samples['Lambda_indices'] != 5
    samples = {key: val[keep] for key, val in samples.items()}
    try:
        dG, ddG = _dgfep('bar', samples, lambdas)
    except Exception:
        return
    assert np.all(_failed_window_ok(dG, ddG, model.dG_map(lambdas)))


@pytest.mark.parametrize('ends', [False, True])
@pytest.mark.parametrize('n_points', [31, 12])
def test_lambda_subsample_loader(tmp_path, ends, n_points):
    """lambda_sub_sample=2 through _load_output_files (Q6).

    31 points is the default ladder; uniform 12 has an even count, so
    sub-sampling with ends=True drops lambda = 1 and the loader appends it
    again (its frames are in the file). Checks: Lambda[Lambda_indices] ==
    Lambda_frame, every window has frames, the kept frames are exactly the
    file's frames at the kept lambdas, and MBAR dG(1) is within 3 sigma of
    exact (plan acceptance Q6).
    """

    model = make_model('SYM')
    if n_points == 31:
        lambdas = default_lambdas()
    else:
        lambdas = np.round(np.linspace(0.0, 1.0, n_points), 3)
    samples = model.sample(lambdas, 200, np.random.default_rng(7))
    files = write_energy_files(tmp_path, samples, lambdas)
    evb = EvbDriver(ostream=OutputStream(None))
    specific, common = evb._load_output_files(*files,
                                              lambda_sub_sample=2,
                                              lambda_sub_sample_ends=ends)
    L = np.asarray(common['Lambda'], dtype=float)
    Lf = np.asarray(common['Lambda_frame'])
    Li = np.asarray(common['Lambda_indices'])
    assert L[0] == 0.0 and L[-1] == 1.0 and np.all(np.diff(L) > 0)
    assert np.array_equal(np.round(L[Li], 3), Lf)
    assert np.all(np.bincount(Li, minlength=len(L)) > 0)
    kept = np.isin(samples['Lambda_frame'], L)
    assert np.array_equal(specific['E1_pes'], samples['E1'][kept])

    sub = {'E1': specific['E1_pes'], 'E2': specific['E2_pes'],
           'Lambda_frame': Lf, 'Lambda_indices': Li}
    dG, ddG = _dgfep('mbar', sub, L)
    assert abs(dG[-1] - model.dG_map(L)[-1]) <= 3.0 * ddG[-1]


def test_hysteresis_injected_bias():
    """A known backward bias gives the predicted endpoint hysteresis (Q8c).

    Four replicas of SYM; each replica's direction-1 frames are copies of
    its direction-0 frames, then b = 5 kJ/mol is added to the backward E2 of
    replica 2. Adding b to E2 shifts Delta G_map(lambda) by lambda b, and the
    code defines delta = dGfep_fwd - dGfep_bwd, so the expected delta is
    -b lambda: endpoint_dGfep_delta = -b (lambda_end - lambda_0) = -5 for
    replica 2 and exactly 0 for the others. Per-replica BAR gives this
    exactly (forward and reverse works shift by opposite constants).
    Tolerance 1e-6 kJ/mol. The signed endpoint delta is most negative for
    replica 2, so hysteresis_summary names it as min_replica of
    endpoint_dGfep_delta and max_replica of max_abs_dGfep_delta.
    """

    model = make_model('SYM')
    lambdas = default_lambdas()
    b = 5.0
    biased = 2
    s = model.sample(lambdas, 200, np.random.default_rng(12), n_replicas=4)
    for k in range(len(lambdas)):
        for rep in range(4):
            block = (s['Lambda_indices'] == k) & (s['Replica_frame'] == rep)
            fwd = block & (s['Direction_frame'] == 0)
            bwd = block & (s['Direction_frame'] == 1)
            s['E1'][bwd] = s['E1'][fwd]
            s['E2'][bwd] = s['E2'][fwd] + (b if rep == biased else 0.0)

    dp = _processor()
    dp.alpha = 0.0
    dp.H12 = STANDARD_SETS['SYM']['H12']
    out = dp.compute(make_results(s, lambdas), 31.0, 0.0)
    res = out['configuration_results']['synthetic']
    for rep, directions in res['replicas'].items():
        delta = np.asarray(directions['hysteresis']['dGfep_delta'])
        expected = -b * lambdas if rep == biased else np.zeros(len(lambdas))
        assert np.max(np.abs(delta - expected)) <= 1e-6
        assert abs(directions['hysteresis']['endpoint_dGfep_delta'] -
                   expected[-1]) <= 1e-6
    summary = res['hysteresis_summary']
    assert summary['endpoint_dGfep_delta']['min_replica'] == biased
    assert summary['max_abs_dGfep_delta']['max_replica'] == biased
