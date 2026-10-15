"""The integration-vs-PES bias and the reweighting tool on toy models (WP-B2).

Dev-only (test_claude_*): must not reach master. The tool lives outside
VeloxChem, in the validation project: it is imported from
~/code/Agent_work/EVB validation plan/tools/evb_reweight.py, so these tests
need that folder. The toy models are MismatchedModel cases of
test_claude_evb_reference (MISMATCH_CASES): window lambda samples a
potential v_int that differs from the PES mapping potential the code
analyses. Energies are in kJ/mol.

Routes under test (reports/WP-B2.md): R, weighted BAR with importance
weights exp(-beta (V_pes - E_own)) per window; M, MBAR with the sampled
states and the PES states unsampled, from exact cross energies U_int or
linearised (M-lin).
"""

import shutil
import sys
from pathlib import Path

import numpy as np

from veloxchem.evbdataprocessing import EvbDataProcessing
from veloxchem.outputstream import OutputStream

from test_claude_evb_reference import (Model1, default_lambdas,
                                       make_mismatched, make_results,
                                       write_energy_files)

sys.path.insert(
    0, str(Path.home() / 'code' / 'Agent_work' / 'EVB validation plan' /
           'tools'))
import evb_reweight  # noqa: E402

LAMBDAS = default_lambdas()


def _data(case, n, seed, cross=False, **kwargs):
    model = make_mismatched(case, **kwargs)
    samples = model.sample(LAMBDAS, n, np.random.default_rng(seed),
                           cross_energies=cross)
    return model, samples, evb_reweight.load(make_results(samples, LAMBDAS))


def test_identity_without_mismatch():
    """With v_int = V_pes (T0) both routes reduce to the code (Q7).

    Route R has unit weights, so it is BAR between adjacent windows and must
    equal the code's BAR; route M with exact cross energies holds every state
    twice, so its PES states must equal the code's MBAR. Both to 1e-6
    kJ/mol (the solvers converge to 1e-12 in reduced units). The continuum
    limit of the code's MBAR must equal dG_map to 1e-6 on an 11-point
    ladder (quadrature on the same grid).
    """

    model, _, data = _data('T0', 200, 1, cross=True)
    dp = EvbDataProcessing(ostream=OutputStream(None))
    dp.Lambda = [float(v) for v in LAMBDAS]
    dp.Lambda_indices = data['Lambda_indices']
    dp.fep_estimator = 'bar'
    bar, _ = dp._calculate_dGfep(data['E1'], data['E2'], data['T'])
    mbar, _ = evb_reweight.code_dGfep(data)
    R = evb_reweight.route_R_point(data, 0.0)
    M = evb_reweight.route_M(data, 0.0, exact=True)
    assert np.allclose(R['ess'], np.bincount(data['Lambda_indices']))
    assert np.max(np.abs(R['dG'] - np.asarray(bar))) < 1e-6
    assert np.max(np.abs(M['dG'] - mbar)) < 1e-6
    short = np.linspace(0, 1, 11)
    assert np.max(np.abs(model.dG_code_limit(short) -
                         model.dG_map(short))) < 1e-6


def test_exact_profile_machinery():
    """MismatchedModel's exact PES profile equals Model1's closed form.

    T1's PES pair is Model1(145, 145, 1, -60) (the EXO set). The root-sum
    gap density on the grid must match the closed form to 1e-6 in ln P over
    +-4 gap stds of both end windows, and stationary(30) must give the same
    barrier and dG0 to 1e-4 kJ/mol (grid search plus Brent refinement on
    an interpolated density).
    """

    mm = make_mismatched('T1')
    m1 = Model1(145.0, 145.0, 1.0, -60.0)
    for lam in (0.0, 1.0):
        mean, std = m1.gap_moments(lam)
        xi = np.linspace(mean - 4 * std, mean + 4 * std, 101)
        assert np.max(np.abs(mm.log_gap_density(xi, lam) -
                             m1.log_gap_density(xi, lam))) < 1e-6
    a, b = mm.stationary(30.0), m1.stationary(30.0)
    assert abs(a['barrier'] - b['barrier']) < 1e-4
    assert abs(a['dG0'] - b['dG0']) < 1e-4


def test_routes_recover_pes_dG0():
    """Routes R and M (exact U_int) recover the PES dG0 on T1 (Q2, Q3).

    T1 samples harmonic diabats with other force constants and centres
    (mismatch std 0.4-1.1 kT); the code's dGfep is off by -28 kJ/mol there
    (its continuum limit). For seeds 0-2 at 300 frames per window, each
    route's dG(1) must lie within 3 sigma of the exact Delta G_map(1) = -60
    (route R: bootstrap sigma, 50 samples; route M: MBAR's sigma), and the
    code's must lie more than 20 sigma away. M-lin equals M here (T1 is
    linear in lambda) to 1e-6.
    """

    exact = -60.0
    for seed in range(3):
        _, _, data = _data('T1', 300, seed, cross=True)
        R = evb_reweight.route_R(data, 0.0, n_boot=50, seed=seed)
        M = evb_reweight.route_M(data, 0.0, exact=True)
        M_lin = evb_reweight.route_M(data, 0.0, exact=False)
        code, code_sigma = evb_reweight.code_dGfep(data)
        assert abs(R['dG'][-1] - exact) < 3 * R['sigma'][-1], seed
        assert abs(M['dG'][-1] - exact) < 3 * M['sigma'][-1], seed
        assert abs(code[-1] - exact) > 20 * code_sigma[-1]
        assert np.max(np.abs(M_lin['dG'] - M['dG'])) < 1e-6


def test_code_converges_to_biased_limit():
    """The code's dGfep on T2 sits on its predicted biased limit (Q1).

    T2 samples the code's dynamic harmonic bonds (nonlinear in lambda)
    while the analysis reads Morse diabats. The prediction is the
    infinite-sample solution of the MBAR equations with the sampled
    densities (dG_code_limit, by quadrature): +16.2 kJ/mol above the exact
    Delta G_map(1) = -60. At 1000 frames per window (seed 0) the code's
    dG(1) must lie within 3 of its sigmas of the limit and more than 5
    sigma from the exact value.
    """

    model, _, data = _data('T2', 1000, 0)
    code, sigma = evb_reweight.code_dGfep(data)
    limit = model.dG_code_limit(LAMBDAS)[-1]
    assert abs(code[-1] - limit) < 3 * sigma[-1], (code[-1], limit)
    assert abs(code[-1] - (-60.0)) > 5 * sigma[-1]


def test_large_mismatch_is_flagged():
    """T3: Kish ESS collapses in the middle windows and R is unreliable (Q6).

    T3 adds 20 spectator modes sampled up to 3.5x softer than the PES
    (mismatch std up to 7.9 kT at lambda 0.5; exact ESS/N there ~1e-8). At
    200 frames per window every window with 0.3 <= lambda <= 0.7 must be
    flagged (ESS < ess_min = 20), the end windows must not be, and analyze()
    must mark route R unreliable; M-lin, linear in lambda only for the
    diabats, must be marked unreliable by its own-lambda residual.
    """

    _, _, data = _data('T3', 200, 0)
    table = evb_reweight.windows_table(data, 0.0)
    mid = (LAMBDAS >= 0.3) & (LAMBDAS <= 0.7)
    assert np.all(table['flagged'][mid])
    assert not table['flagged'][0] and not table['flagged'][-1]
    res = evb_reweight.analyze(data, 0.0, n_boot=5)
    assert res['reliable']['R'] is False
    assert res['reliable']['M_lin'] is False


def test_own_energy_constant_cancels():
    """A per-window constant in E_own changes nothing (plan 6.4 e).

    Ep of T1 frames gets an arbitrary constant per window (up to 1e4
    kJ/mol). Route R's dG, its ESS and lambda_r, and the M-lin residual
    std must be unchanged to 1e-9 (relative to the numbers' size).
    """

    offsets = np.random.default_rng(7).uniform(-1e4, 1e4, len(LAMBDAS))
    out = []
    for ep_offset in (None, lambda lam: np.interp(lam, LAMBDAS, offsets)):
        _, _, data = _data('T1', 100, 3, ep_offset=ep_offset)
        R = evb_reweight.route_R_point(data, 0.0)
        lr = evb_reweight.reweighted_lambda_r(data, 0.0, R['dG'][-1],
                                              R['log_w'])
        table = evb_reweight.windows_table(data, 0.0)
        out.append((R['dG'], R['ess'], lr['sym'],
                    table['mlin_residual_std']))
    for a, b in zip(*out):
        assert np.allclose(a, b, rtol=1e-9, atol=1e-9)


def test_loader_reads_own_energy(tmp_path):
    """load() finds the frame's own sampled energy in run folders.

    A folder written by write_energy_files from T1 samples must load with
    E_own = Ep (own_source 'Ep') and E1_int/E2_int as written (17 digits:
    exact). The old-format Sn2 vacuum fixture must load E_own from its
    E_m_int column (own_source 'E_m_int'), one value per frame.
    """

    _, samples, _ = _data('T1', 20, 0)
    folder = tmp_path / 'toy'
    write_energy_files(folder, samples, LAMBDAS)
    data = evb_reweight.load(folder)
    assert data['own_source'] == 'Ep'
    assert np.array_equal(data['E_own'], samples['Ep'])
    assert np.array_equal(data['E1_int'], samples['E1_int'])
    assert np.array_equal(data['E2_int'], samples['E2_int'])

    src = Path(__file__).parent / 'data'
    fixture = tmp_path / 'sn2'
    fixture.mkdir()
    for name, target in (('evb_Sn2_vacuum_Energies.csv', 'Energies.csv'),
                         ('evb_Sn2_vacuum_Data_combined.csv',
                          'Data_combined.csv'),
                         ('evb_options.json', 'options.json')):
        shutil.copy(src / name, fixture / target)
    data = evb_reweight.load(fixture)
    assert data['own_source'] == 'E_m_int'
    assert len(data['E_own']) == len(data['E1'])
