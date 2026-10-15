"""The external reorganization-energy tool against the exact models (WP-A5).

Dev-only (test_claude_*): must not reach master. The tool lives outside
VeloxChem, in the validation project: it is imported from
~/code/Agent_work/EVB validation plan/tools/evb_reorg.py, so these tests
need that folder. Exact references come from test_claude_evb_reference
(WP-A1). Energies are in kJ/mol.

Conventions under test (reports/WP-A5.md): xi = E1 - E2', de = -xi,
dG0 = G2' - G1. For equal curvature, lambda_r = (<de>_0 - <de>_1)/2
= <de>_0 - dG0 = dG0 - <de>_1, and dg1(xi2) - dg1(xi1) = dg2(xi1) - dg2(xi2)
= lambda_r with the dg_i parabolas of curvature 1/(4 lambda_r).
"""

import sys
from pathlib import Path

import numpy as np
import pytest

from veloxchem.evbdataprocessing import EvbDataProcessing
from veloxchem.evbdriver import EvbDriver
from veloxchem.outputstream import OutputStream

from test_claude_evb_reference import (default_lambdas, make_model,
                                       make_results, write_energy_files)

sys.path.insert(
    0, str(Path.home() / 'code' / 'Agent_work' / 'EVB validation plan' /
           'tools'))
import evb_reorg  # noqa: E402

LAMBDAS = default_lambdas()
N_BOOT = 20
EQUAL_CURVATURE = ['SYM', 'EXO', 'ENDO', 'INV', 'ENZ']
ROUTES = ['R1_sym', 'R1_reactant', 'R1_product', 'R2_left', 'R2_right',
          'R2_curv_1', 'R2_curv_2']
_CACHE = {}


def _data(name, n=200, seed=0, n_replicas=1):
    samples = make_model(name).sample(LAMBDAS, n, np.random.default_rng(seed),
                                      n_replicas=n_replicas)
    return evb_reorg.load(make_results(samples, LAMBDAS, name=name))


def _analyzed(name, n=200, seed=0):
    """analyze() with N_BOOT bootstrap samples, shared between tests."""

    if (name, n, seed) not in _CACHE:
        _CACHE[(name, n, seed)] = evb_reorg.analyze(_data(name, n, seed),
                                                    n_boot=N_BOOT)
    return _CACHE[(name, n, seed)]


def _exact_bin_curves(model, edges):
    """dg1, dg2 integrated over the bins exactly (gap CDFs of WP-A1), and
    the exact bin probabilities on the pure states."""

    w = np.diff(edges)
    p1 = model.gap_cdf(edges[1:], 0.0) - model.gap_cdf(edges[:-1], 0.0)
    p2 = model.gap_cdf(edges[1:], 1.0) - model.gap_cdf(edges[:-1], 1.0)
    with np.errstate(divide='ignore'):
        dg1 = -model.kT * np.log(p1 / w)
        dg2 = -model.kT * np.log(p2 / w) + model.dG_map(1.0)[0]
    dg1[~np.isfinite(dg1)] = np.nan
    dg2[~np.isfinite(dg2)] = np.nan
    return dg1, dg2, p1, p2


def _exact_routes(model, edges, window):
    """Every route evaluated on the exact model with the tool's estimators:
    exact gap moments for R1, exact bin curves for R2."""

    beta = 1.0 / model.kT
    dG0 = model.dG_map(1.0)[0]
    (m0, s0), (m1, s1) = model.gap_moments(0.0), model.gap_moments(1.0)
    de0, de1 = -m0, -m1
    out = dict(dG0=dG0, R1_sym=0.5 * (de0 - de1), R1_reactant=de0 - dG0,
               R1_product=dG0 - de1, R1_var_0=0.5 * beta * s0**2,
               R1_var_1=0.5 * beta * s1**2,
               aqvist_residual=(0.5 * (de0 + de1) - 0.25 * beta *
                                (s0**2 - s1**2) - dG0))
    dg1, dg2, _, _ = _exact_bin_curves(model, edges)
    centers = 0.5 * (edges[1:] + edges[:-1])
    out.update(evb_reorg.lambda_parabola(centers, dg1, dg2, window))
    return out


@pytest.mark.parametrize('name', EQUAL_CURVATURE)
def test_conventions_exact(name):
    """Q1: the section-4 relations hold exactly on equal-curvature models.

    From the A1 closed forms: <de>_0 = dG0 + lambda_r, <de>_1 = dG0 -
    lambda_r, and the R2 forms on the exact pointwise dg1, dg2 (quadratic in
    xi, so the tool's vertex fit is exact). Tolerance 1e-8 kJ/mol for the
    moments and 1e-6 for the fitted parabolas (least squares in floating
    point on values up to ~1e3).
    """

    m = make_model(name)
    lam_r, dG0 = m.lambda_r, m.dG_map(1.0)[0]
    de0, de1 = -m.gap_moments(0.0)[0], -m.gap_moments(1.0)[0]
    assert dG0 == pytest.approx(m.delta, abs=1e-8)
    assert 0.5 * (de0 - de1) == pytest.approx(lam_r, abs=1e-8)
    assert de0 - dG0 == pytest.approx(lam_r, abs=1e-8)
    assert dG0 - de1 == pytest.approx(lam_r, abs=1e-8)

    xi = np.linspace(de1 * -1 - 3 * lam_r, -de0 + 3 * lam_r, 4001)
    xi = np.sort(np.concatenate([xi, [-de0, -de1]]))
    r2 = evb_reorg.lambda_parabola(xi, m.dg1(xi), m.dg2(xi), 4 * m.kT)
    for key in ('R2_left', 'R2_right', 'R2_curv_1', 'R2_curv_2'):
        assert r2[key] == pytest.approx(lam_r, abs=1e-6), key
    assert r2['xi1_min'] == pytest.approx(-de0, abs=1e-6)
    assert r2['xi2_min'] == pytest.approx(-de1, abs=1e-6)


@pytest.mark.parametrize('name', ['SYM', 'EXO', 'ENZ'])
def test_routes_recover_lambda_r(name):
    """Q2-Q4: R1 and R2 on sampled data recover lambda_r and dG0.

    200 frames per window, default ladder, seed 0, 20 bootstrap samples.
    Every route and dG0 must lie within 4 bootstrap sigma of the exact value
    (A1 closed forms: lambda_r = 2 k a^2 for equal curvature, dG0 = delta),
    and every sigma must be finite and nonzero. The sweep (work/WP-A5) checks
    bias and calibration over seeds.
    """

    m = make_model(name)
    res = _analyzed(name)
    v, s = res['values'], res['sigma']
    assert v['dG0'] == pytest.approx(m.delta, abs=4 * s['dG0'])
    for key in ROUTES:
        assert np.isfinite(s[key]) and s[key] > 0, key
        assert abs(v[key] - m.lambda_r) <= 4 * s[key], (key, v[key], s[key])


def test_diabatic_curves_vs_exact():
    """Q3: the tool's dg1 and dg2 match the exact bin-integrated curves.

    SYM, 200 frames per window. Compared bin by bin wherever the exact bin
    probability on the pure state exceeds 1e-6 (the plan's region; outside
    it the exact CDF difference underflows) and the bin holds at least 5
    effective frames (with one or two frames the bootstrap sigma is
    meaningless: the only bins beyond 3.5 sigma held a single frame).
    Expected: |z| <= 4 in at least 99% of those bins and an RMS z within
    [0.6, 1.5] (a 20-sample bootstrap sigma is itself uncertain by ~16%).
    """

    m = make_model('SYM')
    res = _analyzed('SYM')
    data = _data('SYM')
    edges = evb_reorg.default_edges(data, 0.0, 2.0)
    ex1, ex2, p1, p2 = _exact_bin_curves(m, edges)
    for name, exact, p in (('dg1', ex1, p1), ('dg2', ex2, p2)):
        g, sg = res['curves'][name], res['curves'][f'{name}_sigma']
        sel = ((p > 1e-6) & np.isfinite(g) & np.isfinite(sg) & (sg > 0) &
               (res['curves'][f'{name}_neff'] >= evb_reorg.MIN_NEFF))
        z = (g[sel] - exact[sel]) / sg[sel]
        assert sel.sum() > 50, name
        assert np.mean(np.abs(z) <= 4) >= 0.99, (name, np.max(np.abs(z)))
        assert 0.6 <= np.sqrt(np.mean(z**2)) <= 1.5, name


def test_tool_fep_matches_code_dgfep():
    """Cross-check: the tool's MBAR equals EvbDataProcessing's MBAR dGfep.

    Same samples (EXO, 200 frames per window), same kB. Both solve the same
    MBAR equations to convergence, so they agree to 1e-6 kJ/mol.
    """

    samples = make_model('EXO').sample(LAMBDAS, 200, np.random.default_rng(3))
    dp = EvbDataProcessing(ostream=OutputStream(None))
    dp.Lambda = [float(lam) for lam in LAMBDAS]
    dp.Lambda_frame = np.asarray(samples['Lambda_frame'])
    dp.Lambda_indices = np.asarray(samples['Lambda_indices'])
    code, _ = dp._calculate_dGfep_mbar(samples['E1'], samples['E2'], 300.0)
    tool = evb_reorg.fep(samples['E1'], samples['E2'], LAMBDAS,
                         samples['Lambda_indices'], 300.0)['dG']
    assert np.max(np.abs(np.asarray(code) - tool)) <= 1e-6


def test_invariances():
    """Q7: alpha, a common per-frame energy, and swapping the states.

    - alpha: every lambda_r route unchanged and dG0 shifted by alpha, to
      1e-8 kJ/mol (the tool shifts dGfep by lambda*alpha exactly).
    - c_n added to E1 and E2: every value unchanged to 1e-6 (MBAR is
      re-solved; its tolerance, not the algebra, sets this).
    - E1 <-> E2 with the ladder mirrored: R1_reactant and R1_product swap,
      R2_left and R2_right swap, dG0 changes sign, to 1e-6.
    """

    data = _data('EXO', 200, 5)
    base = evb_reorg.analyze(data, n_boot=0)['values']

    shifted = evb_reorg.analyze(data, alpha=7.3, n_boot=0)['values']
    assert shifted['dG0'] == pytest.approx(base['dG0'] + 7.3, abs=1e-8)
    for key in ROUTES:
        assert shifted[key] == pytest.approx(base[key], abs=1e-8), key

    c = np.random.default_rng(1).normal(0.0, 50.0, len(data['E1']))
    common = evb_reorg.analyze(dict(data, E1=data['E1'] + c,
                                    E2=data['E2'] + c), n_boot=0)['values']
    for key in ROUTES + ['dG0']:
        assert common[key] == pytest.approx(base[key], abs=1e-6), key

    K = len(LAMBDAS)
    swapped_data = dict(data, E1=data['E2'], E2=data['E1'],
                        Lambda=np.round(1.0 - data['Lambda'][::-1], 3),
                        Lambda_indices=K - 1 - data['Lambda_indices'])
    swapped = evb_reorg.analyze(swapped_data, n_boot=0)['values']
    assert swapped['dG0'] == pytest.approx(-base['dG0'], abs=1e-6)
    for a, b in (('R1_reactant', 'R1_product'), ('R2_left', 'R2_right'),
                 ('R1_sym', 'R1_sym')):
        assert swapped[a] == pytest.approx(base[b], abs=1e-6), a


def test_curvature_cases_characterized():
    """Q5: CURV2 (k2 = 2 k1): the one-sided routes measure different things.

    Each sampled route must equal the same estimator applied to the exact
    model (exact gap moments; exact bin curves with the same bins and fit
    window) within 4 bootstrap sigma. On this model the reactant-side
    routes (R1_reactant, R2_right) are about the surface-2 energetic
    lambda 2 k2 a^2 = 580 kJ/mol and the product-side ones (R1_product,
    R2_left) about the surface-1 lambda 2 k1 a^2 = 290 kJ/mol; that split is
    asserted to 5%. The Aqvist & Hansson residual is compared the same way.
    """

    m = make_model('CURV2')
    data = _data('CURV2')
    res = evb_reorg.analyze(data, n_boot=N_BOOT)
    edges = evb_reorg.default_edges(data, 0.0, 2.0)
    exact = _exact_routes(m, edges, res['meta']['fit_window_kJ'])
    v, s = res['values'], res['sigma']
    for key in ROUTES + ['dG0', 'R1_var_0', 'R1_var_1', 'aqvist_residual']:
        assert abs(v[key] - exact[key]) <= 4 * s[key], (key, v[key],
                                                         exact[key], s[key])
    lam1, lam2 = 2 * m.k1 * m.a**2, 2 * m.k2 * m.a**2
    for key, ref in (('R1_reactant', lam2), ('R2_right', lam2),
                     ('R1_product', lam1), ('R2_left', lam1)):
        assert exact[key] == pytest.approx(ref, rel=0.05), key


def test_interfaces_identical(tmp_path, monkeypatch):
    """Q9: a results dict, a run folder and a results h5 give one answer.

    The A1 writer stores 17 significant digits, so the three inputs hold
    identical energies; values must agree to 1e-9 kJ/mol. The h5 comes
    from EvbDriver.compute_energy_profiles with alpha and H12 fixed, and
    its stored alpha (2.5) is applied unless overridden. Replica and
    direction labels reach the tool from every interface. A multi-step
    result is refused.
    """

    monkeypatch.chdir(tmp_path)
    samples = make_model('ENZ').sample(LAMBDAS, 50, np.random.default_rng(2),
                                       n_replicas=2)
    write_energy_files(tmp_path / 'ENZ', samples, LAMBDAS)
    drv = EvbDriver(ostream=OutputStream(None))
    drv.dataprocessing.calculate_per_replica = False
    drv.system_confs = [{'name': 'ENZ', 'data_folder': 'ENZ'}]
    drv.compute_energy_profiles(7.0, 0.0, alpha=2.5, H12=20.0)

    from_dict = evb_reorg.load(make_results(samples, LAMBDAS, name='ENZ'))
    from_folder = evb_reorg.load(tmp_path / 'ENZ')
    from_h5 = evb_reorg.load(tmp_path / 'results.h5')
    assert from_h5['alpha'] == 2.5 and from_dict['alpha'] == 0.0

    ref = evb_reorg.analyze(from_dict, alpha=2.5, n_boot=0)['values']
    for data in (from_folder, from_h5):
        alpha = None if data is from_h5 else 2.5
        values = evb_reorg.analyze(data, alpha=alpha, n_boot=0)['values']
        for key in ROUTES + ['dG0']:
            assert values[key] == pytest.approx(ref[key], abs=1e-9), key

    for data in (from_dict, from_folder, from_h5):
        assert len(np.unique(data['Replica_frame'])) == 2
        assert len(np.unique(data['Direction_frame'])) == 2

    path = make_results(samples, LAMBDAS, name='ENZ')
    path['configuration_results']['ENZ']['steps'] = []
    with pytest.raises(ValueError, match='multi-step'):
        evb_reorg.load(path)


def test_decorrelated_bootstrap_and_replica_sem():
    """Decision 1 (WP-A5): sigma follows time correlation; replica SEM.

    SYM, 1000 frames per window, frames from AR(1) chains (A1 sampler: phi
    = (tau-1)/(tau+1), each window, replica and direction its own chain).
    For independent frames sigma(R1_sym) = 0.5 sqrt(2 * 2 lambda_r kT / N)
    = 0.85 kJ/mol; correlated frames multiply it by sqrt(g). With tau_int =
    19 the tool's g must be close to 19 (median over blocks in [12, 26]; the
    estimate is noisy on 500-frame blocks) and its sigma within a factor
    1.6 of 0.85 sqrt(19); with tau = 1 and 4 replicas g must be close to 1,
    sigma within a factor 1.6 of 0.85 (10 to 20 bootstrap samples leave
    ~20% noise), and the replica standard error must be reported.
    """

    m = make_model('SYM')
    expected = 0.5 * np.sqrt(2 * 2 * m.lambda_r * m.kT / 1000)

    def run(tau, n_replicas, n_boot):
        phi = (tau - 1) / (tau + 1)
        samples = m.sample(LAMBDAS, 1000 // n_replicas,
                           np.random.default_rng(0), n_replicas=n_replicas,
                           phi=phi)
        data = evb_reorg.load(make_results(samples, LAMBDAS, name='SYM'))
        return evb_reorg.analyze(data, n_boot=n_boot)

    corr = run(19, 1, N_BOOT)
    assert 12 <= corr['meta']['g_median'] <= 26
    ratio = corr['sigma']['R1_sym'] / (expected * np.sqrt(19))
    assert 1 / 1.6 <= ratio <= 1.6, ratio

    indep = run(1, 4, 10)
    assert indep['meta']['g_median'] <= 1.3
    ratio = indep['sigma']['R1_sym'] / expected
    assert 1 / 1.6 <= ratio <= 1.6, ratio
    assert indep['meta']['n_replicas_sem'] == 4
    assert np.isfinite(indep['sigma_replica']['R1_sym'])
    assert indep['sigma_replica']['R1_sym'] > 0
