"""WP-A4: the alpha/H12 fit of EvbDataProcessing against exact models.

Dev-only (never to master). Targets come from the WP-A1 exact models in
test_claude_evb_reference.py: for a model whose samples carry Delta = 0 and
a chosen (alpha*, H12*), the exact barrier and dG0 are those of
``model.shifted(alpha*).stationary(H12*)``. The fit is then asked to recover
(alpha*, H12*) from the unshifted samples.

"Noise free" means ``dp._calculate_dGfep`` is replaced, on the instance
only, by the exact Delta G_map(lambda) of the model shifted by the alpha the
code is probing. That isolates the fit and the profile reading from the FEP
estimator (A2).

Plan: ~/code/Agent_work/EVB validation plan/plans/WP-A4-alpha-h12-fit.md.
"""

from pathlib import Path
import os

import numpy as np
import pytest
import scipy.optimize

os.environ.setdefault('PYMBAR_DISABLE_JAX', 'true')

from veloxchem.evbdataprocessing import EvbDataProcessing
from veloxchem.evbdriver import EvbDriver
from veloxchem.outputstream import OutputStream
from veloxchem.errorhandler import VeloxChemError

from test_claude_evb_reference import (make_model, default_lambdas,
                                       make_results, write_energy_files)

LAMBDAS = default_lambdas()
# Phase-A bar: 0.05 kcal/mol (grand plan, phase A acceptance).
PHASE_A_BAR = 0.209


@pytest.fixture
def fsolve_spy(monkeypatch):
    """Record what scipy.optimize.fsolve returned inside the code under test.

    The call is forwarded with full_output=True so the residual and ier are
    recorded, and the result is handed back in whatever shape the caller
    asked for.
    """

    real = scipy.optimize.fsolve
    calls = []

    def spy(func, x0, *args, **kwargs):
        wanted_full = kwargs.pop('full_output', False)
        x, info, ier, msg = real(func, x0, *args, full_output=True, **kwargs)
        calls.append(dict(x=x, fvec=info['fvec'], nfev=info['nfev'], ier=ier,
                          msg=msg))
        if wanted_full:
            return x, info, ier, msg
        return x

    monkeypatch.setattr(scipy.optimize, 'fsolve', spy)
    return calls


def _samples(name, n_per_window, seed=0):
    return make_model(name).sample(LAMBDAS, n_per_window,
                                   np.random.default_rng(seed))


def _make_dp(ostream=None):
    dp = EvbDataProcessing(ostream=ostream or OutputStream(None))
    dp.calculate_per_replica = False
    return dp


def _patch_exact_dgfep(dp, configs):
    """Replace dp._calculate_dGfep by the exact Delta G_map.

    configs maps id(E1 array) -> (model, E2 array). The alpha being probed is
    recovered as E2_shifted - E2, which is a constant.
    """

    def exact(E1, E2_shifted, Temp_set, lambda_indices=None):
        model, E2 = configs[id(E1)]
        alpha = float(np.mean(E2_shifted - E2))
        g = model.shifted(alpha).dG_map(np.asarray(dp.Lambda, dtype=float))
        return g, np.zeros_like(g)

    dp._calculate_dGfep = exact


def _noise_free_dp(name, samples, ostream=None):
    results = make_results(samples, LAMBDAS, name=name)
    conf = results['configuration_results'][name]
    dp = _make_dp(ostream)
    _patch_exact_dgfep(dp, {id(conf['E1_pes']): (make_model(name),
                                                 conf['E2_pes'])})
    return dp, results


def _exact_targets(name, alpha, H12):
    st = make_model(name).shifted(alpha).stationary(H12)
    return np.array([st['barrier'], st['dG0']])


def _fit_targets_of_code(dp, results, name, alpha, H12):
    """The code's fit target function at (alpha, H12): barrier and dG0.

    Mirrors the closure in EvbDataProcessing._fit_EVB_parameters (same xi
    grid, fitting=True), so that the bias of the analytical formula itself
    can be separated from what the solver does.
    """

    dp.Lambda = results['Lambda']
    dp.Lambda_frame = results['Lambda_frame']
    dp.Lambda_indices = results['Lambda_indices']
    conf = results['configuration_results'][name]
    # The fit estimates dGfep once, unshifted, and adds Lambda * alpha.
    dGfep, _ = dp._calculate_dGfep(conf['E1_pes'], conf['E2_pes'],
                                   conf['Temp_set'])
    dGfep = np.asarray(dGfep) + np.asarray(dp.Lambda) * alpha
    xi = np.linspace(-10000, 10000, 20000)
    profile, _, _ = dp._dGevb_analytical(dGfep, dp.Lambda, H12, xi)
    _, barrier, dG0, *_ = dp._get_free_energies(profile, xi, fitting=True)
    return np.array([barrier, dG0])


def _exact_jacobian(name, alpha, H12, h=0.5):
    """d(barrier, dG0)/d(alpha, H12) of the exact model, central differences."""

    J = np.zeros((2, 2))
    for j, (da, dh) in enumerate(((h, 0.0), (0.0, h))):
        J[:, j] = (_exact_targets(name, alpha + da, H12 + dh) -
                   _exact_targets(name, alpha - da, H12 - dh)) / (2 * h)
    return J


NOISE_FREE_CASES = [
    ('SYM', 0.0, 50.0),
    ('SYM', -60.0, 30.0),
    ('SYM', 60.0, 30.0),
    ('SYM', 0.0, 100.0),
    ('SYM', 0.0, 5.0),
    ('ENZ', 0.0, 20.0),
]


@pytest.mark.parametrize('name,alpha_star,H12_star', NOISE_FREE_CASES)
def test_fit_roundtrip_noise_free(fsolve_spy, name, alpha_star, H12_star):
    """Q1a: with the exact Delta G_map, the fit recovers (alpha*, H12*).

    Expected: (alpha*, H12*), whose exact barrier and dG0 (A1 quadrature,
    Model.stationary) are the targets. Two checks:
    - the solver: the residual of the code's own target equations at the
      returned point is <= 0.01 kJ/mol (plan section 6, one order below the
      phase-A bar);
    - the parameters: |alpha - alpha*| and |H12 - H12*| are within
      max(0.1, |predicted| + 0.05) kJ/mol, where "predicted" is the noise-free
      bias of the analytical formula on the fit grid (the code's target
      function at the true parameters minus the exact targets) propagated
      through the inverse exact Jacobian. The 0.05 covers the 0.01 residual
      times |J^-1| (<= 2 here) plus the finite-difference Jacobian.
    """

    samples = _samples(name, 200)
    dp, results = _noise_free_dp(name, samples)
    targets = _exact_targets(name, alpha_star, H12_star)
    dp.compute(results, *targets)

    residual = np.abs(fsolve_spy[-1]['fvec']).max()
    assert residual <= 0.01, (
        f'{name} ({alpha_star}, {H12_star}): fit residual {residual:.3g} '
        'kJ/mol at the returned point')

    bias = (_fit_targets_of_code(dp, results, name, alpha_star, H12_star) -
            targets)
    predicted = -np.linalg.solve(_exact_jacobian(name, alpha_star, H12_star),
                                 bias)
    error = np.array([dp.alpha - alpha_star, dp.H12 - H12_star])
    tolerance = np.maximum(0.1, np.abs(predicted) + 0.05)
    assert np.all(np.abs(error) <= tolerance), (
        f'{name} ({alpha_star}, {H12_star}): (alpha, H12) error {error}, '
        f'predicted from the formula bias {predicted}, tolerance {tolerance}')


def test_fit_roundtrip_sampled():
    """Q1b: the fit recovers (alpha*, H12*) from sampled data (MBAR).

    EXO in the plan's sense: SYM samples (Delta = 0), N = 1000 per window,
    seed 0, targets the exact barrier and dG0 at (alpha*, H12*) = (-60, 30).
    Expected: (-60, 30). Tolerance: 4x the across-seed SD of the fitted
    values, alpha 0.235 and H12 0.136 kJ/mol (100 seeds, WP-A4
    work/WP-A4/output/tables_baseline.md, "Q1b sampled round trip"), around
    the mean offset expected from the formula's grid bias (alpha +0.03,
    H12 -0.15 kJ/mol, same table, Q1a).
    """

    samples = _samples('SYM', 1000)
    results = make_results(samples, LAMBDAS, name='SYM')
    dp = _make_dp()
    dp.compute(results, *_exact_targets('SYM', -60.0, 30.0))
    assert abs(dp.alpha - (-60.0 + 0.03)) <= 4 * 0.235
    assert abs(dp.H12 - (30.0 - 0.15)) <= 4 * 0.136


@pytest.mark.parametrize('name,alpha_star,H12_star', [
    ('SYM', 0.0, 50.0),
    ('SYM', -60.0, 30.0),
])
def test_reported_barrier_equals_target(name, alpha_star, H12_star):
    """Q2: after the fit, compute() reports the barrier and dG0 it was fitted to.

    Expected: the targets (exact barrier and dG0 from A1 quadrature).
    Tolerance: 0.209 kJ/mol, the phase-A bar (0.05 kcal/mol). Noise free,
    so the only differences are the reading grids: the fit reads a 1 kJ/mol
    grid with a fixed five-point parabola, the report the 10 kJ/mol
    coordinate bins with the adaptive one.
    """

    samples = _samples(name, 200)
    dp, results = _noise_free_dp(name, samples)
    targets = _exact_targets(name, alpha_star, H12_star)
    out = dp.compute(results, *targets)
    ana = out['configuration_results'][name]['analytical']
    reported = np.array([ana['barrier'], ana['free_energy']])
    assert np.all(np.abs(reported - targets) <= PHASE_A_BAR), (
        f'reported {reported} against targets {targets}')


def test_alpha_sign_convention():
    """Q3: alpha is added to state 2, so a lower product needs alpha < 0.

    Samples are drawn at Delta = 0 (SYM, N = 1000, MBAR). The targets are
    the exact barrier and dG0 of the same model with Delta = -60, H12 = 30
    (A1 table: 20.3628 and -58.6368). Expected alpha = -60 (Delta_target -
    Delta_sampled); tolerance 1 kJ/mol, about 5x the across-seed spread of
    alpha seen in the WP-A4 probe.

    Also: shifting E2 by +a shifts dGfep(lambda = 1) by exactly +a (MBAR,
    u_kn is linear in E2), to 1e-6 kJ/mol.
    """

    samples = _samples('SYM', 1000)
    results = make_results(samples, LAMBDAS, name='SYM')
    dp = _make_dp()
    dp.compute(results, 20.3628, -58.6368)
    assert dp.alpha < 0
    assert abs(dp.alpha + 60.0) <= 1.0

    conf = results['configuration_results']['SYM']
    shift = 7.3
    g0 = np.asarray(
        dp._calculate_dGfep(conf['E1_pes'], conf['E2_pes'], 300.0)[0])
    g1 = np.asarray(
        dp._calculate_dGfep(conf['E1_pes'], conf['E2_pes'] + shift, 300.0)[0])
    assert abs((g1[-1] - g0[-1]) - shift) <= 1e-6
    assert np.allclose(g1 - g0, shift * np.asarray(LAMBDAS), atol=1e-6)


def _driver_on(folder_name, samples, ostream=None):
    write_energy_files(Path(folder_name), samples, LAMBDAS)
    drv = EvbDriver(ostream=ostream or OutputStream(None))
    drv.dataprocessing.calculate_per_replica = False
    drv.system_confs = [{'name': folder_name, 'data_folder': folder_name}]
    return drv


def test_alpha_sign_via_driver(tmp_path, monkeypatch):
    """Q3 end to end: EvbDriver.compute_energy_profiles on A1-written files.

    Expected: the same alpha and H12 as compute() on the same samples, to
    1e-6 (the writer stores 17 significant digits, so the data are
    identical), and alpha < 0 for the EXO targets on Delta = 0 samples.
    Runs inside tmp_path because the method writes results.h5 to the
    current directory.
    """

    monkeypatch.chdir(tmp_path)
    samples = _samples('SYM', 200)
    drv = _driver_on('SYM', samples)
    drv.compute_energy_profiles(20.3628, -58.6368)

    dp = _make_dp()
    dp.compute(make_results(samples, LAMBDAS, name='SYM'), 20.3628, -58.6368)

    assert drv.dataprocessing.alpha < 0
    assert abs(drv.dataprocessing.alpha - dp.alpha) <= 1e-6
    assert abs(drv.dataprocessing.H12 - dp.H12) <= 1e-6


@pytest.mark.parametrize('name,targets', [
    ('SYM', (90.0, 0.0)),
    ('ENZ', (0.7126, -17.8300)),
])
def test_fit_unattainable_target(fsolve_spy, name, targets):
    """Q5: a fit that does not solve its equations is not returned silently.

    - SYM, barrier 90: above the H12 = 0 Marcus limit of 72.5 kJ/mol, so no
      (alpha, H12) reaches it.
    - ENZ, the exact targets of (alpha, H12) = (-20, 20): a 0.71 kJ/mol
      barrier, below extremum_tolerance (2.5), so the fit target function
      does not see it. fsolve stops with ier = 1 and a residual of 0.8 kJ/mol
      (WP-A4 probe), so checking ier alone would miss it.

    Expected (ground rules, no silent failure): compute() raises
    VeloxChemError, naming the residual. Noise free.
    """

    samples = _samples(name, 200)
    dp, results = _noise_free_dp(name, samples)
    with pytest.raises(VeloxChemError, match='residual'):
        dp.compute(results, *targets)


def test_partial_input_overwritten(tmp_path):
    """Q6: giving alpha alone refits both, and says so with a warning.

    Expected: alpha is refitted (it does not stay at the value given), and a
    '* Warning *' line naming alpha is written to the output stream, since
    the user's input is being discarded. Noise free, SYM (0, 50).
    """

    log = tmp_path / 'out.txt'
    ostream = OutputStream(log)
    samples = _samples('SYM', 200)
    dp, results = _noise_free_dp('SYM', samples, ostream)
    dp.alpha = 5.0
    dp.compute(results, *_exact_targets('SYM', 0.0, 50.0))
    ostream.flush()
    ostream.close()

    assert abs(dp.alpha - 5.0) > 1.0
    warnings = [line for line in log.read_text().splitlines()
                if line.startswith('* Warning *')]
    assert any('alpha' in line for line in warnings), (
        'no warning about the discarded alpha')


def test_driver_second_call_refits(tmp_path, monkeypatch):
    """Q7: a second compute_energy_profiles call on one driver starts afresh.

    First call: SYM samples with the SYM targets. Second call: the same
    driver pointed at ENZ samples, with the ENZ targets. Expected: alpha, H12
    and the coordinate bins of the second call equal those of a fresh driver
    on the ENZ data (to 1e-6), because nothing the first call derived
    belongs to the second data set.
    """

    monkeypatch.chdir(tmp_path)
    sym, enz = _samples('SYM', 200), _samples('ENZ', 200, seed=1)
    enz_targets = (7.3478, 0.0)

    drv = _driver_on('SYM', sym)
    drv.compute_energy_profiles(31.1207, 0.0)
    write_energy_files(Path('ENZ'), enz, LAMBDAS)
    drv.system_confs = [{'name': 'ENZ', 'data_folder': 'ENZ'}]
    drv.compute_energy_profiles(*enz_targets)

    fresh = _driver_on('ENZ', enz)
    fresh.compute_energy_profiles(*enz_targets)

    second, ref = drv.dataprocessing, fresh.dataprocessing
    assert abs(second.alpha - ref.alpha) <= 1e-6
    assert abs(second.H12 - ref.H12) <= 1e-6
    assert len(second.coordinate_bins) == len(ref.coordinate_bins)
    assert np.allclose(second.coordinate_bins, ref.coordinate_bins,
                       atol=1e-6)


def _transfer(alpha_star, H12_star, bin_size=None):
    """SYM reference plus ENZ second configuration, noise free, fitted to the
    exact SYM targets at (alpha*, H12*). Returns dp and the results."""

    sym, enz = _samples('SYM', 200), _samples('ENZ', 200, seed=1)
    results = make_results(sym, LAMBDAS, name='SYM')
    enz_res = make_results(enz, LAMBDAS, name='ENZ')
    results['configuration_results']['ENZ'] = (
        enz_res['configuration_results']['ENZ'])
    confs = results['configuration_results']
    assert np.array_equal(results['Lambda_frame'], enz_res['Lambda_frame'])

    dp = _make_dp()
    if bin_size is not None:
        dp.bin_size = bin_size
    _patch_exact_dgfep(
        dp, {
            id(confs[name]['E1_pes']): (make_model(name),
                                        confs[name]['E2_pes'])
            for name in ('SYM', 'ENZ')
        })
    out = dp.compute(results, *_exact_targets('SYM', alpha_star, H12_star))
    return dp, out


@pytest.mark.parametrize('alpha_star,H12_star', [(0.0, 20.0), (10.0, 20.0),
                                                 (-5.0, 20.0)])
def test_transfer_second_configuration(alpha_star, H12_star):
    """Q8: calibrate on the first configuration, predict the second.

    Reference: SYM (lambda_r 290); second: ENZ (lambda_r 92), both sampled
    at Delta = 0, noise free. The fit uses only the reference, whose targets
    are the exact values at (alpha*, H12*). Expected: the ENZ profile the
    code builds with the fitted (alpha, H12) has the exact ENZ barrier and
    dG0 under those parameters. It is read on a 0.1 kJ/mol grid over the
    shared coordinate range, so that the reading of the default bins (see
    test_transfer_reported_second_configuration) stays out of it.
    Tolerance: the phase-A bar, 0.209 kJ/mol (measured: <= 0.02).
    """

    dp, out = _transfer(alpha_star, H12_star)
    conf = out['configuration_results']['ENZ']
    fine = np.arange(dp.coordinate_bins[0], dp.coordinate_bins[-1], 0.1)
    profile, _, _ = dp._dGevb_analytical(np.asarray(conf['dGfep']),
                                         dp.Lambda, dp.H12, fine)
    _, barrier, dG0, *_ = dp._get_free_energies(profile, fine, fitting=False)

    expected = _exact_targets('ENZ', dp.alpha, dp.H12)
    assert np.all(np.abs(np.array([barrier, dG0]) - expected) <= PHASE_A_BAR)


@pytest.mark.parametrize('alpha_star,H12_star', [(0.0, 20.0), (10.0, 20.0),
                                                 (-5.0, 20.0)])
def test_transfer_reported_second_configuration(alpha_star, H12_star):
    """Q8, reported (A4-5): compute()'s analytical values for the second
    configuration, read on bins shared with a wider one.

    Expected: the exact ENZ barrier and dG0 under the fitted (alpha, H12),
    within the phase-A bar, 0.209 kJ/mol. The bins are computed once from
    both configurations, so they span SYM's much wider gap range. Before the
    fix ENZ was read over all of it: at (10, 20) its product well is 3.5
    kJ/mol deep and the analytical profile beyond it rises only ~3 kJ/mol
    over the rest of SYM's range, and the adaptive parabola that took this
    shoulder in put dG0 0.58 kJ/mol too high. Now each configuration is read
    over its own sampled range.

    bin_size = 2: on the default 10 kJ/mol bins the ENZ barrier (about
    2*H12 = 40 kJ/mol wide) is itself under-resolved, by -0.15 to -0.27
    depending on where the bins fall, which is the bin resolution (skill
    trap T24), not this defect.
    """

    dp, out = _transfer(alpha_star, H12_star, bin_size=2)
    ana = out['configuration_results']['ENZ']['analytical']
    reported = np.array([ana['barrier'], ana['free_energy']])
    expected = _exact_targets('ENZ', dp.alpha, dp.H12)
    assert np.all(np.abs(reported - expected) <= PHASE_A_BAR), (
        f'ENZ reported {reported}, exact under the fitted parameters '
        f'{expected}')


def test_fit_rejects_impossible_targets():
    """Memo 3 (c): a target barrier at or below the target dG0 is refused
    before the fit starts.

    The transition state is the highest point between the wells, so on any
    profile barrier >= dG0, and a barrier is positive. Expected: a
    VeloxChemError saying so, for the (5, 10) the old test_evb_dataproc
    asked for and for a negative barrier.
    """

    samples = _samples('SYM', 200)
    for targets in ((5.0, 10.0), (-1.0, -5.0)):
        dp, results = _noise_free_dp('SYM', samples)
        with pytest.raises(VeloxChemError, match='larger than the target'):
            dp.compute(results, *targets)


def test_fit_spurious_small_H12_root_raises(fsolve_spy):
    """Memo 1 (A4-4): a fit that lands below what its xi grid resolves is
    refused.

    From H12_guess = 1, fsolve finds an exact root of the fit's equations at
    |H12| = 0.25 for the SYM (0, 50) targets (residual 4e-6, WP-A4 guess
    grid). There the fit grid skips over the top of the barrier, so the
    root is an artefact; the report afterwards has no barrier. Expected: a
    VeloxChemError naming the artefact, although the residual check passes.
    """

    samples = _samples('SYM', 200)
    dp, results = _noise_free_dp('SYM', samples)
    dp.H12_guess = 1.0
    with pytest.raises(VeloxChemError, match='artefact of the fit grid'):
        dp.compute(results, *_exact_targets('SYM', 0.0, 50.0))
    assert np.abs(fsolve_spy[-1]['fvec']).max() <= 0.01


def test_fit_report_mismatch_warns(tmp_path):
    """Memo 1 (Q2 LOWH): a correct fit whose barrier the bins cannot resolve
    gets a warning.

    LOWH (H12 = 5): the barrier is about 10 kJ/mol wide, so the default 10
    kJ/mol bins put one or two points on it and the reported barrier is ~10
    kJ/mol below the target (WP-A4, Q2), while the fit itself is right.
    Expected: H12 within 0.2 of 5, no error, and a '* Warning *' line that
    names fit_report_tolerance.
    """

    log = tmp_path / 'out.txt'
    ostream = OutputStream(log)
    samples = _samples('SYM', 200)
    dp, results = _noise_free_dp('SYM', samples, ostream)
    dp.compute(results, *_exact_targets('SYM', 0.0, 5.0))
    ostream.flush()
    ostream.close()

    assert abs(dp.H12 - 5.0) <= 0.2
    warnings = [line for line in log.read_text().splitlines()
                if line.startswith('* Warning *')]
    assert any('fit_report_tolerance' in line for line in warnings)


def test_fit_max_iter_bounds_evaluations(fsolve_spy):
    """Memo 4 (A4-6): max_iter caps the evaluations of the fit's target.

    Expected: with max_iter = 3, fsolve stops early (ier = 5, "reached
    maxfev") and compute() refuses the unfinished fit. fsolve checks maxfev
    between iterations, and its finite-difference Jacobian costs n + 1 = 3
    calls, so it may go up to 3 over (measured: 6); a converged fit here
    takes 9.
    """

    samples = _samples('SYM', 200)
    dp, results = _noise_free_dp('SYM', samples)
    dp.max_iter = 3
    with pytest.raises(VeloxChemError, match='maxfev'):
        dp.compute(results, *_exact_targets('SYM', 0.0, 50.0))
    assert fsolve_spy[-1]['nfev'] <= 3 + 3


def test_fit_returns_positive_H12():
    """Memo 5: only H12**2 enters, and the fit reports the magnitude.

    From H12_guess = -10, fsolve converges to H12 = -50 for the SYM (0, 50)
    targets. Expected: the fitted H12 is +50 (within 0.1), and the results
    dict stores the positive value.
    """

    samples = _samples('SYM', 200)
    dp, results = _noise_free_dp('SYM', samples)
    dp.H12_guess = -10.0
    out = dp.compute(results, *_exact_targets('SYM', 0.0, 50.0))
    assert abs(dp.H12 - 50.0) <= 0.1
    assert out['H12'] > 0
