"""Dev tests (WP-B5b, feature F5): native reweighting in EvbDataProcessing.

reweighting="mbar" treats the per-lambda integration systems as the sampled
states (U_int, every frame in every lambda's integration system) and the PES
mapping states and E_g as unsampled MBAR targets. On the mismatched toys of
WP-B2 (MISMATCH_CASES): the identities against reweighting="none" (T0) and the
tool's route M (T0-T3), recovery of the exact dG0 and profiles, the fit on
both profiles and the continuity of its target, the per-replica results, the
diagnostics and the refusals.

Sweeps: work/WP-B5/toy_sweep.py and fit_continuity.py in the validation
project.
"""

import copy

import numpy as np
import pytest

from veloxchem.evbdataprocessing import EvbDataProcessing
from veloxchem.evbdriver import EvbDriver
from veloxchem.outputstream import OutputStream

from test_claude_evb_reference import (add_evb_tools_to_path,
                                       default_lambdas, make_mismatched,
                                       make_results)

LAMBDAS = default_lambdas()
H12 = {'T0': 30.0, 'T1': 30.0, 'T2': 15.0, 'T3': 15.0}


def _results(case, n=300, seed=0, n_replicas=1):
    model = make_mismatched(case)
    samples = model.sample(LAMBDAS, n, np.random.default_rng(seed),
                           n_replicas=n_replicas, cross_energies=True)
    return model, make_results(samples, LAMBDAS)


def _compute(results,
             reweighting='mbar',
             alpha=0.0,
             H12_value=30.0,
             targets=(np.nan, np.nan),
             warnings=None,
             **settings):
    stream = OutputStream(None)
    if warnings is not None:
        stream.print_warning = warnings.append
    dp = EvbDataProcessing(ostream=stream)
    dp.reweighting = reweighting
    dp.calculate_per_replica = False
    for key, value in settings.items():
        setattr(dp, key, value)
    if alpha is not None:
        dp.alpha = alpha
        dp.H12 = H12_value
    out = dp.compute(copy.deepcopy(results), *targets)
    return out, out['configuration_results']['synthetic']


def _reweighted(results):
    """An EvbDataProcessing set up to reweight results, and its reweighted
    mapping free energies and log-density over all frames (no fit)."""

    dp = EvbDataProcessing(ostream=OutputStream(None))
    dp.reweighting = 'mbar'
    dp.results = results['configuration_results']
    dp.Lambda = results['Lambda']
    dp.Lambda_frame = results['Lambda_frame']
    dp.Lambda_indices = results['Lambda_indices']
    dp._check_reweighting_settings()
    conf = dp.results['synthetic']
    return dp, dp._reweight('synthetic', conf['E1_pes'], conf['E2_pes'],
                            conf['Temp_set'], dp.Lambda_indices)


def _route_M(results, **kwargs):
    add_evb_tools_to_path()
    import evb_reweight
    return evb_reweight.route_M(evb_reweight.load(results), 0.0, exact=True,
                                **kwargs)


def test_t0_mbar_equals_none():
    """Sampled = PES: the reweighting changes nothing (plan identity 1)."""

    _, results = _results('T0', n=150)
    _, none = _compute(results, 'none')
    _, mbar = _compute(results, 'mbar')
    np.testing.assert_allclose(mbar['dGfep'], none['dGfep'], atol=1e-6)
    np.testing.assert_allclose(mbar['dGfep_unc'], none['dGfep_unc'],
                               atol=1e-6)
    np.testing.assert_allclose(mbar['analytical']['EVB'],
                               none['analytical']['EVB'],
                               atol=1e-6)
    np.testing.assert_allclose(mbar['reweighting']['dGfep_unreweighted'],
                               none['dGfep'],
                               atol=1e-12)
    assert 'discrete' not in none and 'discrete' not in mbar
    assert 'binless' in mbar and 'binless' not in none


@pytest.mark.parametrize('case', ['T0', 'T1', 'T2', 'T3'])
def test_mbar_equals_route_M(case):
    """The code's reweighted dGfep, sigma and mapping-state ESS are the
    tool's route M exact (plan identity 2)."""

    _, results = _results(case, n=100)
    _, mbar = _compute(results, H12_value=H12[case])
    M = _route_M(results)
    np.testing.assert_allclose(mbar['dGfep'], M['dG'], atol=1e-6)
    np.testing.assert_allclose(mbar['dGfep_unc'], M['sigma'], atol=1e-6)
    np.testing.assert_allclose(mbar['reweighting']['ess_pes'],
                               M['ess_pes'],
                               rtol=1e-6)


def test_reweighted_lambda_r_matches_route_M_weights():
    """R1 from the weights of the mapping end states (B3's analyze.py
    lambda_r with the route M weights)."""

    _, results = _results('T1', n=150)
    _, mbar = _compute(results)
    M = _route_M(results, return_weights=True)
    conf = results['configuration_results']['synthetic']
    de = conf['E2_pes'] - conf['E1_pes']
    K = len(LAMBDAS)
    m0 = np.sum(M['weights'][:, 0] * de) / np.sum(M['weights'][:, 0])
    m1 = np.sum(M['weights'][:, K - 1] * de) / np.sum(M['weights'][:, K - 1])
    reweighting = mbar['reweighting']
    assert reweighting['lambda_r'] == pytest.approx(0.5 * (m0 - m1),
                                                    abs=1e-6)
    assert reweighting['lambda_reactant'] == pytest.approx(
        m0 - M['dG'][-1], abs=1e-6)
    assert reweighting['lambda_product'] == pytest.approx(
        M['dG'][-1] - m1, abs=1e-6)


def test_t1_recovers_the_pes():
    """T1 (harmonic, sampled diabats off the PES): reweighted dG0 within 3
    SE of the exact -60 over 4 seeds, where none sits about -28 off (the
    code limit); the analytical and binless barriers near the exact
    stationary barrier (toy_sweep: SD 0.25 and 0.35 at N 300)."""

    model = make_mismatched('T1')
    exact = model.stationary(30.0)
    errors = {'dG1': [], 'a': [], 'b': []}
    for seed in range(4):
        _, results = _results('T1', n=200, seed=seed)
        _, mbar = _compute(results, bin_size=2.0, calculate_reorganization=False)
        errors['dG1'].append(mbar['dGfep'][-1] + 60.0)
        errors['a'].append(mbar['analytical']['barrier'] - exact['barrier'])
        errors['b'].append(mbar['binless']['barrier'] - exact['barrier'])
    for key in ('dG1', 'a', 'b'):
        err = np.asarray(errors[key])
        se = np.std(err, ddof=1) / np.sqrt(len(err))
        assert abs(np.mean(err)) <= 3 * se + 0.1, (key, err)
    assert mbar['reweighting']['dGfep_unreweighted'][-1] + 60.0 < -20


def test_binless_matches_dg_ground_on_t2():
    """T2 (Morse PES, dynamic-bond sampling): the binless barrier against
    the exact E_g profile, where the analytical formula carries its own
    bias (-1.29 kJ/mol here). One seed at N 700 (SD about 0.2); the kernel
    (width 5) lowers the barrier by 0.24 (toy_sweep)."""

    model = make_mismatched('T2')
    exact = model.stationary(15.0)
    _, results = _results('T2', n=700, seed=7)
    _, mbar = _compute(results,
                       H12_value=15.0,
                       bin_size=2.0,
                       calculate_reorganization=False)
    assert mbar['binless']['barrier'] - exact['barrier'] == pytest.approx(
        -0.24, abs=0.7)
    assert mbar['binless']['free_energy'] == pytest.approx(exact['dG0'],
                                                           abs=0.8)


@pytest.mark.parametrize('profile', ['analytical', 'binless'])
def test_fit_on_either_profile(profile):
    """The alpha/H12 fit on the reweighted profile it is told to use; the
    profile reports the targets. T1's PES needs alpha 0 and H12 30."""

    model = make_mismatched('T1')
    exact = model.stationary(30.0)
    _, results = _results('T1', seed=3)
    out, mbar = _compute(results,
                         alpha=None,
                         targets=(exact['barrier'], exact['dG0']),
                         reweighting_profile=profile,
                         H12_guess=20.0)
    assert out['alpha'] == pytest.approx(0.0, abs=0.6)
    assert out['H12'] == pytest.approx(30.0, abs=1.0)
    reported = mbar[profile]
    assert reported['barrier'] == pytest.approx(exact['barrier'], abs=1.0)
    assert reported['free_energy'] == pytest.approx(exact['dG0'], abs=1.0)
    assert mbar['reweighting']['profile'] == profile


@pytest.mark.parametrize('case', ['T1', 'T2'])
@pytest.mark.parametrize('scanned', ['H12', 'alpha'])
def test_binless_fit_target_is_continuous(case, scanned):
    """The barrier the binless fit reads over a 20 kJ/mol scan of H12 or
    alpha (201 points): largest step at most 2x the median (plan Q7).
    Without the dense-range mask, isolated tail frames made a 159x step in
    the H12 scan on T1."""

    _, results = _results(case, seed=3)
    conf = results['configuration_results']['synthetic']
    dp, reweighting = _reweighted(results)
    setup = dp._binless_setup(conf['E1_pes'], conf['E2_pes'],
                              reweighting['log_den'], conf['Temp_set'])
    barriers = []
    for shift in np.linspace(-10, 10, 201):
        alpha, H12_value = ((0.0, H12[case] + shift) if scanned == 'H12' else
                            (shift, H12[case]))
        profile = dp._binless_free_energy(setup, alpha, H12_value)
        barriers.append(
            dp._get_free_energies(profile, setup['gap'] - alpha,
                                  fitting=True)[1])
    steps = np.abs(np.diff(barriers))
    assert np.max(steps) <= 2 * np.median(steps)


def test_per_replica_results():
    """Per replica and direction the same reweighting (no BAR, no BAR
    warning), with the binless barrier and the reweighted lambda_r, and
    their replica averages."""

    _, results = _results('T1', n=100, n_replicas=3)
    warnings = []
    _, mbar = _compute(results, warnings=warnings, calculate_per_replica=True)
    assert not any('BAR' in w for w in warnings)
    replicas = mbar['replicas']
    assert len(replicas) == 3
    for directions in replicas.values():
        for label in ('forward', 'backward', 'both'):
            assert np.isfinite(directions[label]['binless']['barrier'])
        assert 'reweighted_lambda_r' in directions['both']
        assert 'hysteresis' in directions
    average = mbar['replica_average']
    for key in ('binless_barrier_mean', 'binless_barrier_sem',
                'binless_free_energy_mean', 'reweighted_lambda_r_mean',
                'reweighted_lambda_r_sem'):
        assert np.isfinite(average[key]), key
    # one replica's "both" is that replica's frames alone, reweighted
    single = copy.deepcopy(results)
    keep = single['configuration_results']['synthetic']['Replica_frame'] == 1
    conf = single['configuration_results']['synthetic']
    for key in ('E1_pes', 'E2_pes', 'U_int', 'Replica_frame',
                'Direction_frame', 'E1_int', 'E2_int', 'Ep', 'E_m_pes'):
        conf[key] = conf[key][keep]
    single['Lambda_frame'] = single['Lambda_frame'][keep]
    single['Lambda_indices'] = single['Lambda_indices'][keep]
    _, direct = _reweighted(single)
    np.testing.assert_allclose(replicas[1]['both']['dGfep'],
                               direct['dGfep'],
                               atol=1e-6)


def test_ess_and_warnings():
    """The ESS diagnostics and the warnings below reweighting_ess_min."""

    _, results = _results('T1', n=150)
    warnings = []
    _, mbar = _compute(results, warnings=warnings)
    reweighting = mbar['reweighting']
    assert reweighting['ess_pes'].shape == (len(LAMBDAS), )
    assert reweighting['ess_eg'].shape == (3, )
    assert np.all(reweighting['ess_eg'] > 50)
    assert not any('Kish' in w for w in warnings)

    warnings = []
    _compute(results, warnings=warnings, reweighting_ess_min=1e9)
    assert any('mapping states have a Kish' in w for w in warnings)
    assert any('binless EVB profile has a Kish' in w for w in warnings)


@pytest.mark.parametrize('setting, match', [
    ({'reweighting': 'route_M'}, "unknown reweighting"),
    ({'reweighting_profile': 'discrete'}, "unknown reweighting_profile"),
    ({'fep_estimator': 'bar'}, "no BAR form"),
    ({'reweighting_kernel_width': 0.0}, "must be positive"),
])
def test_bad_settings_are_refused(setting, match):
    _, results = _results('T0', n=20)
    with pytest.raises(Exception, match=match):
        _compute(results, **setting)


def test_auto_is_the_default_and_falls_back():
    """reweighting="auto" (default): "mbar" with U_int; "none" with a
    warning without it, or with fep_estimator="bar"."""

    assert EvbDataProcessing(ostream=OutputStream(None)).reweighting == 'auto'
    _, results = _results('T1', n=100)
    _, mbar = _compute(results, 'mbar')
    warnings = []
    _, auto = _compute(results, 'auto', warnings=warnings)
    assert auto['dGfep'] == mbar['dGfep'] and 'binless' in auto
    assert auto['reweighting']['method'] == 'mbar'
    assert not any('U_int' in w for w in warnings)

    missing = copy.deepcopy(results)
    del missing['configuration_results']['synthetic']['U_int']
    _, none = _compute(missing, 'none')
    warnings = []
    _, auto = _compute(missing, 'auto', warnings=warnings)
    assert auto['dGfep'] == none['dGfep'] and 'binless' not in auto
    assert any('report_integration_energies=True' in w for w in warnings)

    warnings = []
    _, auto = _compute(results, 'auto', warnings=warnings, fep_estimator='bar')
    assert 'binless' not in auto
    assert any('no reweighted form' in w for w in warnings)


@pytest.mark.parametrize('bad', [np.nan, np.inf, 1e19])
@pytest.mark.parametrize('key', ['U_int', 'E1_pes'])
def test_broken_energies_are_refused(key, bad):
    """Review finding 3: a non-finite or blown-up energy is an error that
    names it, not a NaN profile or a misleading fit failure."""

    _, results = _results('T1', n=50)
    conf = results['configuration_results']['synthetic']
    conf[key] = conf[key].copy()
    conf[key][17] = bad
    with pytest.raises(Exception, match=f"blown-up energies .* in {key}"):
        _compute(results)


def test_outlying_frames_are_warned_about():
    """A frame whose own-window integration energy is far out of line
    gets a warning, and the analysis goes on."""

    _, results = _results('T1', n=100)
    warnings = []
    _compute(results, warnings=warnings)
    assert not any('robust standard' in w for w in warnings)
    conf = results['configuration_results']['synthetic']
    conf['U_int'] = conf['U_int'].copy()
    frame = 250
    conf['U_int'][frame, results['Lambda_indices'][frame]] += 1e4
    warnings = []
    _, mbar = _compute(results, warnings=warnings)
    assert any('1 frame(s)' in w and 'robust standard' in w
               for w in warnings)
    assert np.all(np.isfinite(mbar['dGfep']))


def test_missing_or_misshapen_cross_energies_are_refused():
    _, results = _results('T0', n=20)
    missing = copy.deepcopy(results)
    del missing['configuration_results']['synthetic']['U_int']
    with pytest.raises(Exception, match="report_integration_energies"):
        _compute(missing)
    _compute(missing, 'none')
    misshapen = copy.deepcopy(results)
    conf = misshapen['configuration_results']['synthetic']
    conf['U_int'] = conf['U_int'][:, :-1]
    with pytest.raises(Exception, match="shape"):
        _compute(misshapen)


def test_printing_and_h5_round_trip(tmp_path, monkeypatch, capsys):
    import sys
    _, results = _results('T1', n=100, n_replicas=2)
    out, mbar = _compute(results, calculate_per_replica=True)
    stream = OutputStream(sys.stdout)
    EvbDataProcessing.print_results(out, stream)
    EvbDataProcessing.print_uncertainties(out, stream)
    stream.flush()
    printed = capsys.readouterr().out
    assert 'Binless (reweighted)' in printed
    assert 'Reweighting (MBAR' in printed
    assert 'synthetic (binless)' in printed

    monkeypatch.chdir(tmp_path)
    driver = EvbDriver(ostream=OutputStream(None))
    driver._save_dict_as_h5(out, 'reweighted')
    back = EvbDriver._load_dict_from_h5(
        tmp_path / 'reweighted.h5')['configuration_results']['synthetic']
    assert back['binless']['barrier'] == mbar['binless']['barrier']
    for key in ('ess_pes', 'ess_eg', 'dGfep_unreweighted', 'de_mean'):
        np.testing.assert_array_equal(np.asarray(back['reweighting'][key]),
                                      np.asarray(mbar['reweighting'][key]))
    assert back['reweighting']['method'] == 'mbar'
    assert back['reweighting']['lambda_r'] == mbar['reweighting']['lambda_r']


def test_driver_passes_the_keyword(tmp_path, monkeypatch):
    _, results = _results('T1', n=100)
    monkeypatch.chdir(tmp_path)
    driver = EvbDriver(ostream=OutputStream(None))
    driver.dataprocessing.ostream = OutputStream(None)
    driver.dataprocessing.calculate_per_replica = False
    monkeypatch.setattr(driver, '_load_output_from_folders',
                        lambda *args: copy.deepcopy(results))
    out = driver.compute_energy_profiles(np.nan,
                                         np.nan,
                                         alpha=0.0,
                                         H12=30.0,
                                         reweighting='mbar')
    assert 'binless' in out['configuration_results']['synthetic']
    assert driver.dataprocessing.reweighting == 'mbar'
    out = driver.compute_energy_profiles(np.nan, np.nan, alpha=0.0, H12=30.0)
    assert 'binless' in out['configuration_results']['synthetic']
    out = driver.compute_energy_profiles(np.nan,
                                         np.nan,
                                         alpha=0.0,
                                         H12=30.0,
                                         reweighting='none')
    assert 'binless' not in out['configuration_results']['synthetic']


def test_solver_protocol_survives_a_stalled_adaptive_step():
    """B5-5: with pymbar's JAX backend (on by default), handing a
    non-converged adaptive result on to L-BFGS-B ('continuation') raised a
    TypeError (f_k must be a numpy array). One adaptive step forces the
    hand-over; the solution must match a plain solve. Needs JAX on to
    reproduce (PYMBAR_DISABLE_JAX unset)."""

    import pymbar
    rng = np.random.default_rng(0)
    x = np.concatenate([rng.normal(k, 1.0, 200) for k in range(6)])
    u_kn = 1.5 * (x[None, :] - np.arange(6)[:, None])**2
    N_k = np.full(6, 200)
    protocol = EvbDataProcessing._mbar_solver_protocol()
    protocol[0]['options']['maxiter'] = 1
    stalled = pymbar.MBAR(u_kn,
                          N_k,
                          solver_protocol=protocol,
                          maximum_iterations=1)
    np.testing.assert_allclose(np.asarray(stalled.f_k),
                               np.asarray(pymbar.MBAR(u_kn, N_k).f_k),
                               atol=1e-6)


def test_plot_evb_details_refuses_without_a_discrete_profile():
    """B5-3: plot_evb_details read result['discrete'] unconditionally, a
    bare KeyError without a discretised profile (now the default), and
    reweighted results have no FEP/US decomposition to draw."""

    matplotlib = pytest.importorskip('matplotlib')
    matplotlib.use('Agg')
    _, results = _results('T1', n=50)
    for settings in ({}, {'calculate_discrete': True}):
        out, _ = _compute(results, **settings)
        with pytest.raises(Exception, match="no FEP/US decomposition"):
            EvbDataProcessing.plot_evb_details(out)
    out, _ = _compute(results)
    with pytest.raises(Exception, match="no discretised profile"):
        EvbDataProcessing.plot_results(out, plot_discrete=True)


def test_reweighted_discrete_profile_on_request():
    """calculate_discrete under reweighting: the E_g target weights summed
    per coordinate bin, read like the other profiles. T1 at bin_size 2
    against the exact barrier and dG0 (toy_sweep SDs 0.25-0.47 at N 300)."""

    model = make_mismatched('T1')
    exact = model.stationary(30.0)
    _, results = _results('T1', seed=1)
    out, mbar = _compute(results, calculate_discrete=True, bin_size=2.0)
    discrete = mbar['discrete']
    assert len(discrete['EVB']) == len(out['coordinate_bins']) + 1
    assert discrete['barrier'] == pytest.approx(exact['barrier'], abs=1.5)
    assert discrete['free_energy'] == pytest.approx(exact['dG0'], abs=1.5)
    assert 'pns' not in discrete
