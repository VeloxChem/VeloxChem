"""Smoke tests that stand in for the WP-A2 to WP-A5 sweeps.

Dev-only (test_claude_*): must not reach master. Each test runs a few seeds
of the cells a sweep's verdict rests on and applies that verdict's
acceptance rule with the error bar the few seeds allow. The sweeps under
<validation project>/work/ are run only when a test here
fails; each docstring names its sweep.

Exact references come from test_claude_evb_reference (WP-A1). Energies are in
kJ/mol. The phase-A accuracy bar is 0.209 kJ/mol (0.05 kcal/mol).
"""

import numpy as np

from veloxchem.evbdataprocessing import EvbDataProcessing
from veloxchem.outputstream import OutputStream

from test_claude_evb_reference import (STANDARD_SETS, T_DEFAULT,
                                       default_lambdas, make_model,
                                       make_results)

PHASE_A_BAR = 0.209


def _processor():
    dp = EvbDataProcessing(ostream=OutputStream(None))
    dp.calculate_per_replica = False
    return dp


def _assert_unbiased(errors, offset=0.0, label=''):
    """Mean error within 3 SE of the expected offset, and within the bar.

    The bar is widened by 3 SE: with a handful of seeds the mean is not
    known better than that, and the sweep is the place to resolve it.
    """

    errors = np.asarray(errors, dtype=float)
    mean = np.mean(errors)
    se = np.std(errors, ddof=1) / np.sqrt(len(errors))
    assert abs(mean - offset) <= 3 * se + 1e-6, (label, mean, se)
    assert abs(mean) <= abs(offset) + PHASE_A_BAR + 3 * se, (label, mean, se)


def test_smoke_fep_mbar():
    """MBAR dGfep is unbiased and its sigma calibrated (WP-A2 Q1, Q2).

    SYM and CURV2, default 31-point ladder, 200 frames per window, seeds
    0-9. At lambda = 0.5 and 1 the error against the exact Delta G_map must
    have mean within 3 SE of 0 (WP-A2 measured <= 0.11 kJ/mol at N >= 200),
    and z = error / sigma must have variance in [0.15, 3.1]: the 99.9%
    chi-squared interval for 10 seeds of a calibrated sigma. That catches
    a sigma off by about 2x or more (A2-1 gave var z 1.8-28).

    If this fails, run work/WP-A2/run_sweep.py.
    """

    lambdas = default_lambdas()
    for name in ('SYM', 'CURV2'):
        model = make_model(name)
        exact = np.asarray(model.dG_map(lambdas))
        errors, z = [], []
        for seed in range(10):
            samples = model.sample(lambdas, 200, np.random.default_rng(seed))
            dp = _processor()
            dp.Lambda = [float(lam) for lam in lambdas]
            dp.Lambda_frame = np.asarray(samples['Lambda_frame'])
            dp.Lambda_indices = np.asarray(samples['Lambda_indices'])
            dG, ddG = dp._calculate_dGfep(np.asarray(samples['E1']),
                                          np.asarray(samples['E2']),
                                          T_DEFAULT)
            err = (np.asarray(dG) - exact)[[15, -1]]
            errors.append(err)
            z.append(err / np.asarray(ddG)[[15, -1]])
        errors, z = np.array(errors), np.array(z)
        for k, lam in enumerate((0.5, 1.0)):
            _assert_unbiased(errors[:, k], label=f'{name} lambda {lam}')
            assert 0.15 <= np.var(z[:, k], ddof=1) <= 3.1, (name, lam)


def test_smoke_profile():
    """The analytical barrier and dG0 match the exact model (WP-A3).

    SYM, EXO and ENZ at their own H12 with alpha = 0 (no fit), 500 frames
    per window, seeds 0-4, default settings. Exact values: the stationary
    points of the A1 model. WP-A3 measured errors of -0.03 to -0.14 kJ/mol
    for these sets after the rewrite, so the mean must be within the
    phase-A bar (plus 3 SE).

    If this fails, run work/WP-A3/run_bias_map.py.
    """

    lambdas = default_lambdas()
    for name in ('SYM', 'EXO', 'ENZ'):
        model = make_model(name)
        H12 = STANDARD_SETS[name]['H12']
        exact = model.stationary(H12)
        barrier_err, dG0_err = [], []
        for seed in range(5):
            samples = model.sample(lambdas, 500, np.random.default_rng(seed))
            dp = _processor()
            dp.alpha, dp.H12 = 0.0, H12
            out = dp.compute(make_results(samples, lambdas, name=name), 0.0,
                             0.0)
            ana = out['configuration_results'][name]['analytical']
            barrier_err.append(ana['barrier'] - exact['barrier'])
            dG0_err.append(ana['free_energy'] - exact['dG0'])
        _assert_unbiased(barrier_err, label=f'{name} barrier')
        _assert_unbiased(dG0_err, label=f'{name} dG0')


def test_smoke_fit():
    """The alpha/H12 fit recovers the true parameters (WP-A4 Q1b, Q3).

    Delta = 0 samples of SYM fitted to the exact targets of SYM shifted by
    alpha = -60 with H12 = 30 (the EXO case, which also pins the sign of
    alpha), 200 frames per window, seeds 0-4. WP-A4 measured at N = 200 a
    mean error of -0.10 to +0.06 kJ/mol with an SD of 0.52-0.62 (alpha) and
    0.21-0.60 (H12), and noise-free offsets of at most 0.15 kJ/mol, so the
    mean error must lie within the phase-A bar plus 3 SE.

    If this fails, run work/WP-A4/run_fit_sweep.py.
    """

    lambdas = default_lambdas()
    alpha_star, H12_star = -60.0, 30.0
    targets = make_model('SYM').shifted(alpha_star).stationary(H12_star)
    alpha_err, H12_err = [], []
    for seed in range(5):
        samples = make_model('SYM').sample(lambdas, 200,
                                           np.random.default_rng(seed))
        dp = _processor()
        dp.compute(make_results(samples, lambdas, name='SYM'),
                   targets['barrier'], targets['dG0'])
        alpha_err.append(dp.alpha - alpha_star)
        H12_err.append(dp.H12 - H12_star)
    _assert_unbiased(alpha_err, label='alpha')
    _assert_unbiased(H12_err, label='H12')


def test_smoke_reorg():
    """The lambda_r tool recovers lambda_r and dG0 (WP-A5).

    tools/evb_reorg.py in the validation project, on EXO (lambda_r 290,
    Delta -60) and ENZ (lambda_r 92, Delta 0), 1000 frames per window,
    seeds 0-4, no bootstrap. For equal curvature R1 (symmetric), R2 (both
    sides) and dG0 have exact answers lambda_r and Delta (WP-A1). The mean
    error must lie within the phase-A bar plus 3 SE. WP-A5 measured R2 low
    by about 0.5 kJ/mol at this N for lambda_r 290 (a finite-sample bias),
    inside that bound. Resolution: about 1.5 kJ/mol on EXO (0.5%); a 1%
    error in R1 was caught, 0.5% at N = 200 was not.

    If this fails, run work/WP-A5/run_sweep.py.
    """

    from test_claude_evb_reference import add_evb_tools_to_path
    add_evb_tools_to_path()
    import evb_reorg

    lambdas = default_lambdas()
    for name in ('EXO', 'ENZ'):
        model = make_model(name)
        errors = {key: [] for key in ('R1_sym', 'R2_left', 'R2_right',
                                      'dG0')}
        for seed in range(5):
            samples = model.sample(lambdas, 1000,
                                   np.random.default_rng(seed))
            data = evb_reorg.load(make_results(samples, lambdas, name=name))
            values = evb_reorg.analyze(data, n_boot=0)['values']
            for key in errors:
                exact = model.delta if key == 'dG0' else model.lambda_r
                errors[key].append(values[key] - exact)
        for key, err in errors.items():
            _assert_unbiased(err, label=f'{name} {key}')


def test_smoke_reorg_code():
    """EvbDataProcessing's lambda_r and dG0 against the exact answers (WP-A6).

    EXO (lambda_r 290, Delta -60) and ENZ (lambda_r 92, Delta 0), 1000
    frames per window, seeds 0-4, alpha 0 and H12 fixed, so there is no fit.
    For equal curvature lambda_r (and both one-sided values) are the exact
    lambda_r and dG0 is Delta (WP-A1); the mean error must lie within the
    phase-A bar plus 3 SE, and the one-sided split below the tolerance. An
    injected 1% error in lambda_r was caught (WP-A6).

    If this fails, run work/WP-A5/run_sweep.py on the tool, then compare the
    code against it with tests/test_claude_evb_reorg_code.py.
    """

    lambdas = default_lambdas()
    for name in ('EXO', 'ENZ'):
        model = make_model(name)
        errors = {key: [] for key in ('lambda_r', 'lambda_reactant',
                                      'lambda_product', 'dG0')}
        for seed in range(5):
            samples = model.sample(lambdas, 1000,
                                   np.random.default_rng(seed))
            results = make_results(samples, lambdas, name=name)
            dp = _processor()
            dp.calculate_discrete = False
            dp.calculate_analytical = False
            # Fixed bins skip the KDE; the curves are not checked here.
            dp.coordinate_bins = np.arange(-1000.0, 1000.1, 10.0)
            dp.alpha = 0.0
            dp.H12 = STANDARD_SETS[name]['H12']
            dp.compute(results, 0.0, 0.0)
            result = results['configuration_results'][name]
            reorganization = result['reorganization']
            for key in errors:
                if key == 'dG0':
                    errors[key].append(result['dGfep'][-1] - model.delta)
                else:
                    errors[key].append(reorganization[key] - model.lambda_r)
            assert reorganization['split'] < reorganization['split_tolerance']
        for key, err in errors.items():
            _assert_unbiased(err, label=f'{name} {key}')


def test_smoke_reweight():
    """The reweighting tool recovers the PES dG0 from mismatched frames (WP-B2).

    tools/evb_reweight.py in the validation project, on the T1 toy of
    MISMATCH_CASES: frames sampled from harmonic diabats other than the
    analysed PES pair, where the code's dGfep(1) is off by -28 kJ/mol.
    Route R (importance weights within each window, then BAR) at 1000
    frames per window, seeds 0-9, and route M (MBAR with exact cross
    energies) at 300 frames, seeds 0-4, without bootstrap. The exact answer
    is Delta G_map(1) = -60. The mean error must lie within 3 SE of 0 and
    within the phase-A bar plus 3 SE (WP-B2 measured R -0.09 +- 0.05 and
    M -0.02 +- 0.03 over 50 seeds). Flipping the sign of the weights moves
    R by tens of kJ/mol (work/WP-B2/output/smoke_mutation.out).

    If this fails, run work/WP-B2/run_sweep.py.
    """

    from test_claude_evb_reference import (add_evb_tools_to_path,
                                           make_mismatched)
    add_evb_tools_to_path()
    import evb_reweight

    lambdas = default_lambdas()
    model = make_mismatched('T1')
    errors = {'R': [], 'M': []}
    for seed in range(10):
        samples = model.sample(lambdas, 1000, np.random.default_rng(seed))
        data = evb_reweight.load(make_results(samples, lambdas))
        errors['R'].append(evb_reweight.route_R_point(data, 0.0)['dG'][-1] +
                           60.0)
    for seed in range(5):
        samples = model.sample(lambdas, 300, np.random.default_rng(seed),
                               cross_energies=True)
        data = evb_reweight.load(make_results(samples, lambdas))
        errors['M'].append(
            evb_reweight.route_M(data, 0.0, exact=True)['dG'][-1] + 60.0)
    for key, err in errors.items():
        _assert_unbiased(err, label=f'T1 route {key}')


def test_smoke_native_reweight():
    """EvbDataProcessing with reweighting='mbar' recovers the PES dG0 and
    the E_g barrier from mismatched frames (WP-B5b).

    T1 of MISMATCH_CASES: frames sampled from harmonic diabats other than
    the analysed PES pair, where reweighting='none' puts dGfep(1) 28 kJ/mol
    off. 300 frames per window, seeds 0-4, alpha 0 and H12 30 fixed. The
    errors of dGfep(1) against the exact -60 and of the binless barrier
    against the exact E_g barrier must be unbiased (WP-B5b measured +0.08
    +- 0.09 and -0.02 +- 0.07 over 25 seeds). Taking the sampled states'
    energies from the mapping potentials instead of U_int, or swapping the
    sampled and the target blocks of the MBAR matrix, fails it
    (work/WP-B5/smoke_mutation.out).

    If this fails, run work/WP-B5/toy_sweep.py.
    """

    from test_claude_evb_reference import make_mismatched

    lambdas = default_lambdas()
    model = make_mismatched('T1')
    exact = model.stationary(30.0)
    errors = {'dGfep(1)': [], 'binless barrier': []}
    for seed in range(5):
        samples = model.sample(lambdas, 300, np.random.default_rng(seed),
                               cross_energies=True)
        dp = _processor()
        dp.reweighting = 'mbar'
        dp.calculate_reorganization = False
        dp.bin_size = 2.0
        dp.alpha = 0.0
        dp.H12 = 30.0
        result = dp.compute(make_results(samples, lambdas), np.nan,
                            np.nan)['configuration_results']['synthetic']
        errors['dGfep(1)'].append(result['dGfep'][-1] + 60.0)
        errors['binless barrier'].append(result['binless']['barrier'] -
                                         exact['barrier'])
    for key, err in errors.items():
        _assert_unbiased(err, label=f'T1 reweighted {key}')
