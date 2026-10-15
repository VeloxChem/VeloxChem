"""Reorganization energy and diabatic curves in EvbDataProcessing (WP-A6).

Dev-only (test_claude_*): must not reach master. The reference
implementation is the WP-A5 tool, imported from
~/code/Agent_work/EVB validation plan/tools/evb_reorg.py, so these tests
need that folder. Exact references come from test_claude_evb_reference
(WP-A1). Energies are in kJ/mol.

Conventions (reports/WP-A5.md): xi = E1 - E2', de = -xi, dG0 = G2' - G1 =
dGfep[-1]. lambda_r = (<de>_0 - <de>_1)/2, reactant side <de>_0 - dG0,
product side dG0 - <de>_1. dg1(s) = -kT ln P1(s) and dg2(s) = -kT ln P2(s)
+ dG0 per bin (bin probabilities; the tool's curves are densities).
"""

import sys
from pathlib import Path

import numpy as np
import pytest

from veloxchem.evbdataprocessing import EvbDataProcessing
from veloxchem.evbdriver import EvbDriver
from veloxchem.outputstream import OutputStream

from test_claude_evb_reference import (STANDARD_SETS, default_lambdas,
                                       make_model, make_results)

sys.path.insert(
    0, str(Path.home() / 'code' / 'Agent_work' / 'EVB validation plan' /
           'tools'))
import evb_reorg  # noqa: E402

LAMBDAS = default_lambdas()
NEW_KEYS = ('lambda_r', 'lambda_r_unc', 'lambda_reactant',
            'lambda_reactant_unc', 'lambda_product', 'lambda_product_unc',
            'split', 'split_tolerance', 'de_mean', 'de_var', 'g', 'n_eff',
            'dg1', 'dg2', 'convention')
# Fixed 10 kJ/mol bins over every gap these models sample: deriving the bins
# (a KDE) costs about half of a compute(), and only the curve test needs them.
FIXED_BINS = np.arange(-1000.0, 1000.1, 10.0)
_CACHE = {}


def _samples(name, n=200, seed=0, n_replicas=1, phi=0.0):
    return make_model(name).sample(LAMBDAS, n, np.random.default_rng(seed),
                                   n_replicas=n_replicas, phi=phi)


def _compute(samples, name, alpha=0.0, per_replica=False, ostream=None,
             reorganization=True, lambdas=LAMBDAS, bins=FIXED_BINS):
    """compute() with alpha and H12 fixed (no fit), on the given bins (None:
    derived by compute). Returns (results, the configuration's result, the
    processor)."""

    results = make_results(samples, lambdas, name=name)
    dp = EvbDataProcessing(ostream=ostream or OutputStream(None))
    dp.calculate_per_replica = per_replica
    dp.calculate_reorganization = reorganization
    dp.alpha = alpha
    dp.H12 = STANDARD_SETS[name]['H12']
    if bins is not None:
        dp.coordinate_bins = bins
    dp.compute(results, 0.0, 0.0)
    return results, results['configuration_results'][name], dp


def _cached(name, alpha=0.0, n=200, seed=0, bins=FIXED_BINS):
    key = (name, alpha, n, seed, bins is None)
    if key not in _CACHE:
        samples = _samples(name, n, seed)
        _CACHE[key] = (samples, _compute(samples, name, alpha, bins=bins))
    return _CACHE[key]


@pytest.mark.parametrize('name', ['SYM', 'EXO', 'CURV2'])
@pytest.mark.parametrize('alpha', [0.0, -60.0])
def test_r1_equals_tool(name, alpha):
    """(a), (b): lambda_r and both one-sided values are the tool's R1 to
    1e-9 with the code's dG0, the gap moments are the samples' own, and
    the code's dG0 is the tool's MBAR dG0 to 1e-6."""

    samples, (_, result, dp) = _cached(name, alpha)
    reorg = result['reorganization']
    E2s = samples['E2'] + alpha
    xi = samples['E1'] - E2s
    li = samples['Lambda_indices']
    T = result['Temp_set']
    dG0 = result['dGfep'][-1]
    tool = evb_reorg.lambda_energy_gap(xi, li, LAMBDAS, dG0, T)
    assert reorg['lambda_r'] == pytest.approx(tool['R1_sym'], abs=1e-9)
    assert reorg['lambda_reactant'] == pytest.approx(tool['R1_reactant'],
                                                     abs=1e-9)
    assert reorg['lambda_product'] == pytest.approx(tool['R1_product'],
                                                    abs=1e-9)
    beta = 1.0 / (dp.kb * T)
    last = len(LAMBDAS) - 1
    for i, k in enumerate((0, last)):
        de = -xi[li == k]
        assert reorg['de_mean'][i] == pytest.approx(de.mean(), abs=1e-9)
        assert reorg['de_var'][i] == pytest.approx(de.var(ddof=1), rel=1e-9)
    assert reorg['de_var'][0] == pytest.approx(2 * tool['R1_var_0'] / beta,
                                               rel=1e-9)
    assert reorg['de_var'][1] == pytest.approx(2 * tool['R1_var_1'] / beta,
                                               rel=1e-9)
    fep = evb_reorg.fep(samples['E1'], E2s, LAMBDAS, li, T)
    assert dG0 == pytest.approx(fep['dG'][-1], abs=1e-6)


def _expected_sigma(model, n, tau, T):
    kT = EvbDataProcessing(ostream=OutputStream(None)).kb * T
    return 0.5 * np.sqrt(2 * 2 * model.lambda_r * kT / n) * np.sqrt(tau)


@pytest.mark.parametrize('tau', [1, 19])
def test_sigma_with_time_correlation(tau):
    """(c): sigma(lambda_r) within a factor 1.3 of
    0.5 sqrt(2 var(de) tau / N) on AR(1) frames (SYM, N = 1000)."""

    phi = (tau - 1) / (tau + 1)
    samples = _samples('SYM', 1000, seed=3, phi=phi)
    _, result, _ = _compute(samples, 'SYM')
    reorg = result['reorganization']
    expected = _expected_sigma(make_model('SYM'), 1000, tau,
                               result['Temp_set'])
    ratio = reorg['lambda_r_unc'] / expected
    assert 1 / 1.3 < ratio < 1.3, (ratio, reorg['g'])
    assert reorg['n_eff'][0] == pytest.approx(1000 / reorg['g'][0])


def _replica_run():
    if 'replicas' not in _CACHE:
        samples = _samples('SYM', 1000, seed=1, n_replicas=4)
        _CACHE['replicas'] = _compute(samples, 'SYM', per_replica=True)
    return _CACHE['replicas']


def test_replica_average_of_lambda_r():
    """(c): with 4 replicas each pooled replica carries lambda_r and the
    replica average gives a finite mean and SEM, of the size of the
    analytic sigma."""

    _, result, _ = _replica_run()
    per = [d['both']['lambda_r'] for d in result['replicas'].values()]
    assert len(per) == 4 and np.all(np.isfinite(per))
    average = result['replica_average']
    assert average['lambda_r_mean'] == pytest.approx(np.mean(per))
    assert average['lambda_r_sem'] == pytest.approx(
        np.std(per, ddof=1) / 2)
    assert np.isfinite(average['lambda_r_sem'])
    # Four replicas of the same pooled data: the mean of the per-replica
    # lambda_r is the pooled one when every replica has equal end windows.
    assert average['lambda_r_mean'] == pytest.approx(
        result['reorganization']['lambda_r'], abs=1e-9)


def test_diabatic_curves_against_exact():
    """(d): dg1, dg2 against the exact bin curves on SYM (N = 1000).

    Compared where the exact bin probability exceeds 1e-6 and the bin
    holds at least 20 frames. Tolerance: 4x the tool's bootstrap sigma on
    the same bins (12 samples), floored at 0.1 kJ/mol so that a bin the
    bootstrap happens to see as noise-free still has a finite band.
    """

    samples, (_, result, dp) = _cached('SYM', 0.0, n=1000, bins=None)
    model = make_model('SYM')
    reorg = result['reorganization']
    bins = np.asarray(dp.coordinate_bins)
    kT = dp.kb * result['Temp_set']
    # Interior entries 1..len(bins)-1 are the bins between the edges.
    dg1 = np.asarray(reorg['dg1'])[1:-1]
    dg2 = np.asarray(reorg['dg2'])[1:-1]
    p1 = model.gap_cdf(bins[1:], 0.0) - model.gap_cdf(bins[:-1], 0.0)
    p2 = model.gap_cdf(bins[1:], 1.0) - model.gap_cdf(bins[:-1], 1.0)
    with np.errstate(divide='ignore'):
        exact1 = -kT * np.log(p1)
        exact2 = -kT * np.log(p2) + model.dG_map(1.0)[0]
    counts = np.histogram(samples['E1'] - samples['E2'], bins=bins)[0]

    data = evb_reorg.load(make_results(samples, LAMBDAS, name='SYM'))
    tool = evb_reorg.analyze(data, n_boot=12, edges=bins)
    w = np.diff(bins)
    for code, exact, p, name in ((dg1, exact1, p1, 'dg1'),
                                 (dg2, exact2, p2, 'dg2')):
        sigma = np.fmax(tool['curves'][f'{name}_sigma'], 0.1)
        tool_curve = tool['curves'][name] - kT * np.log(w)
        use = (p > 1e-6) & (counts >= 20)
        assert np.count_nonzero(use) >= 10
        z = (code[use] - exact[use]) / sigma[use]
        assert np.all(np.abs(z) < 4), (name, np.max(np.abs(z)))
        # The two estimators (FEP/US here, MBAR in the tool) agree too.
        both = use & np.isfinite(tool_curve)
        assert np.all(np.abs(code[both] - tool_curve[both]) < 4 *
                      sigma[both])
    centers = 0.5 * (bins[1:] + bins[:-1])
    finite = np.isfinite(dg1) & np.isfinite(dg2)
    diff = (dg2 - dg1)[finite] + centers[finite]
    assert np.all(np.abs(diff) <= 0.5 * w[finite] + 1e-9)


def test_invariances():
    """(e): lambda_r does not depend on alpha; swapping the states with the
    ladder mirrored swaps the one-sided values and negates dG0."""

    _, (_, r0, _) = _cached('EXO', 0.0)
    _, (_, r60, _) = _cached('EXO', -60.0)
    for key in ('lambda_r', 'lambda_reactant', 'lambda_product'):
        assert r60['reorganization'][key] == pytest.approx(
            r0['reorganization'][key], abs=1e-6 if key != 'lambda_r' else
            1e-9)
    assert r60['dGfep'][-1] - r0['dGfep'][-1] == pytest.approx(-60.0,
                                                               abs=1e-6)

    samples = _samples('EXO')
    last = len(LAMBDAS) - 1
    mirrored = dict(samples, E1=samples['E2'], E2=samples['E1'],
                    Lambda_frame=1.0 - samples['Lambda_frame'],
                    Lambda_indices=last - samples['Lambda_indices'])
    lambdas = np.round(1.0 - np.asarray(LAMBDAS)[::-1], 3)
    assert np.allclose(lambdas, LAMBDAS)
    _, rs, _ = _compute(mirrored, 'EXO', lambdas=lambdas)
    a, b = r0['reorganization'], rs['reorganization']
    assert b['lambda_r'] == pytest.approx(a['lambda_r'], abs=1e-9)
    assert b['lambda_reactant'] == pytest.approx(a['lambda_product'],
                                                 abs=1e-6)
    assert b['lambda_product'] == pytest.approx(a['lambda_reactant'],
                                                abs=1e-6)
    assert rs['dGfep'][-1] == pytest.approx(-r0['dGfep'][-1], abs=1e-6)


@pytest.mark.parametrize('name, fires', [('CURV2', True), ('SYM', False)])
def test_split_note(name, fires, tmp_path):
    """(f): the linear-response note fires on CURV2 (split about 0.67) and
    not on SYM."""

    log = tmp_path / 'reorg.log'
    _, result, dp = _compute(_samples(name), name,
                             ostream=OutputStream(str(log)))
    dp.ostream.flush()
    split = result['reorganization']['split']
    if fires:
        assert split == pytest.approx(2 / 3, abs=0.05)
    else:
        assert split < 0.05
    text = log.read_text()
    assert ('linear response' in text) == fires
    out = tmp_path / 'print.log'
    stream = OutputStream(str(out))
    EvbDataProcessing.print_uncertainties(
        {'configuration_results': {name: result}}, stream)
    EvbDataProcessing.print_results({'configuration_results': {name: result}},
                                    stream)
    stream.flush()
    printed = out.read_text()
    assert 'Reorganization energy' in printed
    assert ('linear response' in printed) == fires


def test_h5_round_trip(tmp_path, monkeypatch):
    """(g): the reorganization dict and the replica fields survive
    _save_dict_as_h5 / _load_dict_from_h5."""

    results, result, _ = _replica_run()
    monkeypatch.chdir(tmp_path)
    driver = EvbDriver(ostream=OutputStream(None))
    driver._save_dict_as_h5(results, 'reorg')
    loaded = EvbDriver._load_dict_from_h5(tmp_path / 'reorg.h5')
    back = loaded['configuration_results']['SYM']
    for key in NEW_KEYS:
        a, b = result['reorganization'][key], back['reorganization'][key]
        if key == 'convention':
            assert b == a
        else:
            np.testing.assert_array_equal(np.asarray(b), np.asarray(a))
    for replica, directions in result['replicas'].items():
        assert back['replicas'][str(replica)]['both']['lambda_r'] == (
            directions['both']['lambda_r'])
    for key in ('lambda_r_mean', 'lambda_r_sem'):
        assert back['replica_average'][key] == result['replica_average'][key]


def test_switch_off_leaves_results_unchanged():
    """(h): calculate_reorganization = False gives exactly the old key set
    and values."""

    samples = _samples('SYM', 200, seed=2, n_replicas=3)
    _, on, _ = _compute(samples, 'SYM', per_replica=True)
    _, off, _ = _compute(samples, 'SYM', per_replica=True,
                         reorganization=False)
    assert 'reorganization' in on and 'reorganization' not in off
    assert set(on) - set(off) == {'reorganization'}
    assert all('lambda_r' not in d['both'] for d in off['replicas'].values())
    assert not any(k.startswith('lambda_r') for k in off['replica_average'])
    for key in off:
        if key in ('replicas', 'replica_average'):
            continue
        a, b = off[key], on[key]
        if isinstance(a, dict):
            for sub in a:
                np.testing.assert_array_equal(np.asarray(b[sub]),
                                              np.asarray(a[sub]))
        else:
            np.testing.assert_array_equal(np.asarray(b), np.asarray(a))
    for key, value in off['replica_average'].items():
        assert on['replica_average'][key] == value


def _gaussian_results(lambda_r, n=200, seed=0, T=300.0, lambdas=LAMBDAS):
    """Linear-response gap model with dG0 = 0: de = E2 - E1 in window l is
    normal with mean lambda_r (1 - 2l) and variance 2 kT lambda_r; E1 = 0."""

    rng = np.random.default_rng(seed)
    kT = EvbDataProcessing().kb * T
    de = np.concatenate([
        rng.normal(lambda_r * (1 - 2 * lam), np.sqrt(2 * kT * lambda_r), n)
        for lam in lambdas
    ])
    indices = np.repeat(np.arange(len(lambdas)), n)
    return {
        'Lambda': list(lambdas),
        'Lambda_frame': np.asarray(lambdas)[indices],
        'Lambda_indices': indices,
        'configuration_results': {
            'G': {
                'E1_pes': np.zeros_like(de),
                'E2_pes': de,
                'Temp_set': T
            }
        },
    }


def test_diabatic_curves_finite_at_large_lambda_r():
    """dg1, dg2 stay finite at lambda_r = 2500 kJ/mol. There exp(-beta
    lambda xi) of the end windows underflows; an unshifted log-mean-exp put
    about 80 of 562 bins at +inf (review of 36161b1d4)."""

    results = _gaussian_results(2500.0)
    dp = EvbDataProcessing(ostream=OutputStream(None))
    dp.alpha, dp.H12 = 0.0, 50.0
    bins = np.arange(-3500.0, 3500.1, 10.0)
    dp.coordinate_bins = bins
    dp.compute(results, 0.0, 0.0)
    reorg = results['configuration_results']['G']['reorganization']
    dg1, dg2 = np.asarray(reorg['dg1']), np.asarray(reorg['dg2'])
    assert not np.any(np.isinf(dg1)) and not np.any(np.isinf(dg2))
    centers = 0.5 * (bins[1:] + bins[:-1])
    finite = np.isfinite(dg1[1:-1]) & np.isfinite(dg2[1:-1])
    assert np.count_nonzero(finite) > 300
    diff = (dg2 - dg1)[1:-1][finite] + centers[finite]
    assert np.all(np.abs(diff) <= 5.0 + 1e-9)


def test_ladder_must_run_from_0_to_1():
    """compute() refuses a ladder without lambda = 1, with or without the
    reorganization energy: dG0 is defined there."""

    lambdas = LAMBDAS[:-1]
    for reorganization in (True, False):
        dp = EvbDataProcessing(ostream=OutputStream(None))
        dp.alpha, dp.H12 = 0.0, 50.0
        dp.calculate_reorganization = reorganization
        dp.coordinate_bins = FIXED_BINS
        with pytest.raises(Exception, match='end at 1'):
            dp.compute(_gaussian_results(300.0, lambdas=lambdas), 0.0, 0.0)
