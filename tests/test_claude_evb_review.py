"""Code review of the EVB validation changes (8caf5680a..6b9ef6108).

Dev-only (never to master). One reproducer per approved finding:

- R1: alpha/H12 passed back from an earlier result are used, not re-fitted.
- R4: compute_energy_profiles arguments apply to that call only.
- R2: a replica without a barrier does not turn the replica summaries NaN.
- R3: a short gap range with bin_size < 2 gives the clear error.

Review log: ~/code/Agent_work/EVB validation plan/reports/code-review-log.md.
"""

import numpy as np
import pytest

from veloxchem.evbdataprocessing import EvbDataProcessing
from veloxchem.outputstream import OutputStream
from veloxchem.errorhandler import VeloxChemError

from test_claude_evb_reference import STANDARD_SETS, make_results
from test_claude_evb_fit import LAMBDAS, _samples, _make_dp, _driver_on

EXO = (STANDARD_SETS['EXO']['barrier'], STANDARD_SETS['EXO']['dG0'])


def _count_fits(dp):
    """Wrap dp._fit_EVB_parameters on the instance; returns the call list."""

    calls = []
    fit = dp._fit_EVB_parameters

    def counting():
        calls.append(1)
        return fit()

    dp._fit_EVB_parameters = counting
    return calls


def test_passed_back_fit_is_used():
    """R1, EvbDataProcessing: alpha and H12 taken from a previous result
    (the very objects the fit returned) count as given. Expected: no second
    fit, and the second call uses them unchanged although its targets are
    10 kJ/mol higher."""

    samples = _samples('SYM', 200)
    dp = _make_dp()
    calls = _count_fits(dp)
    first = dp.compute(make_results(samples, LAMBDAS, name='SYM'), *EXO)
    alpha, H12 = first['alpha'], first['H12']

    dp.alpha, dp.H12 = first['alpha'], first['H12']
    dp.compute(make_results(samples, LAMBDAS, name='SYM'), EXO[0] + 10,
               EXO[1])

    assert len(calls) == 1
    assert dp.alpha == alpha and dp.H12 == H12


def test_driver_passed_back_fit_is_used(tmp_path, monkeypatch):
    """R1, EvbDriver: compute_energy_profiles(alpha=results['alpha'],
    H12=results['H12']) uses them. Expected: one fit over both calls."""

    monkeypatch.chdir(tmp_path)
    drv = _driver_on('SYM', _samples('SYM', 200))
    calls = _count_fits(drv.dataprocessing)
    drv.compute_energy_profiles(*EXO)
    alpha, H12 = drv.results['alpha'], drv.results['H12']

    drv.compute_energy_profiles(EXO[0] + 10,
                                EXO[1],
                                alpha=drv.results['alpha'],
                                H12=drv.results['H12'])

    assert len(calls) == 1
    assert drv.results['alpha'] == alpha and drv.results['H12'] == H12


def test_driver_arguments_apply_to_one_call(tmp_path, monkeypatch):
    """R4: alpha, H12 and dE_range given to one compute_energy_profiles
    call do not carry over. Expected: the second call, without them, fits
    (alpha and H12 equal to a fresh driver's, to 1e-6) and derives its
    bins (equal to the fresh driver's)."""

    monkeypatch.chdir(tmp_path)
    samples = _samples('SYM', 200)
    drv = _driver_on('SYM', samples)
    calls = _count_fits(drv.dataprocessing)
    drv.compute_energy_profiles(*EXO,
                                alpha=-50.0,
                                H12=40.0,
                                dE_range=(-600.0, 600.0))
    assert len(calls) == 0
    drv.compute_energy_profiles(*EXO)
    assert len(calls) == 1

    fresh = _driver_on('fresh', samples)
    fresh.compute_energy_profiles(*EXO)

    second, ref = drv.dataprocessing, fresh.dataprocessing
    assert abs(second.alpha - ref.alpha) <= 1e-6
    assert abs(second.H12 - ref.H12) <= 1e-6
    assert np.array_equal(second.coordinate_bins, ref.coordinate_bins)


def _replica(barrier, free_energy):
    return {'analytical': {'barrier': barrier, 'free_energy': free_energy}}


def test_replica_average_skips_nan_barrier(tmp_path):
    """R2: one replica out of three has no barrier (NaN, a legitimate A3
    result). Expected: the barrier mean and SEM come from the other two,
    the free energy from all three, and a warning says one was left out."""

    log = tmp_path / 'r2.log'
    dp = EvbDataProcessing(ostream=OutputStream(str(log)))
    # The hand-built replicas carry no lambda_r (that needs compute()).
    dp.calculate_reorganization = False
    result = {
        'replicas': {
            0: {'both': _replica(30.0, -5.0)},
            1: {'both': _replica(np.nan, -6.0)},
            2: {'both': _replica(31.0, -4.0)},
        }
    }
    dp._calculate_replica_average(result)
    dp.ostream.flush()

    average = result['replica_average']
    assert average['barrier_mean'] == pytest.approx(30.5)
    assert average['barrier_sem'] == pytest.approx(
        np.std([30.0, 31.0], ddof=1) / np.sqrt(2))
    assert average['free_energy_mean'] == pytest.approx(-5.0)
    assert average['free_energy_sem'] == pytest.approx(
        np.std([-5.0, -6.0, -4.0], ddof=1) / np.sqrt(3))
    assert '1 of 3' in log.read_text()


def _hysteresis(barrier_delta, other):
    return {
        'hysteresis': {
            'max_abs_dGfep_delta': other,
            'mean_abs_dGfep_delta': other,
            'endpoint_dGfep_delta': other,
            'barrier_delta': barrier_delta,
            'free_energy_delta': np.nan,
        }
    }


def test_hysteresis_summary_skips_nan(tmp_path):
    """R2: barrier_delta is NaN for replica 1, free_energy_delta for all.
    Expected: barrier_delta statistics over replicas 0 and 2, with min and
    max naming them; free_energy_delta, with nothing to summarize, is left
    out; the complete metrics are untouched."""

    log = tmp_path / 'r2h.log'
    dp = EvbDataProcessing(ostream=OutputStream(str(log)))
    result = {
        'replicas': {
            0: _hysteresis(2.0, 1.0),
            1: _hysteresis(np.nan, 2.0),
            2: _hysteresis(-1.0, 3.0),
        }
    }
    dp._calculate_hysterisis_data(result)
    dp.ostream.flush()

    summary = result['hysteresis_summary']
    barrier = summary['barrier_delta']
    assert barrier['mean'] == pytest.approx(0.5)
    assert barrier['std'] == pytest.approx(np.std([2.0, -1.0], ddof=1))
    assert barrier['min'] == -1.0 and barrier['min_replica'] == 2
    assert barrier['max'] == 2.0 and barrier['max_replica'] == 0
    assert barrier['median'] == pytest.approx(0.5)
    assert 'free_energy_delta' not in summary
    assert summary['max_abs_dGfep_delta']['mean'] == pytest.approx(2.0)
    assert summary['max_abs_dGfep_delta']['median_replica'] == 1
    assert 'barrier_delta' in log.read_text()


def _short_gap_results(span, seed=0):
    rng = np.random.default_rng(seed)
    n_lambda, n = len(LAMBDAS), 200
    E1 = rng.normal(0.0, 1.0, n_lambda * n)
    E2 = E1 + rng.uniform(0.0, span, n_lambda * n)
    return {
        'Lambda': [float(lam) for lam in LAMBDAS],
        'Lambda_indices': np.repeat(np.arange(n_lambda), n),
        'Lambda_frame': np.repeat(np.asarray(LAMBDAS, dtype=float), n),
        'configuration_results': {
            'X': {
                'E1_pes': E1,
                'E2_pes': E2,
                'Temp_set': 300.0
            }
        },
    }


def test_short_gap_range_small_bin_size():
    """R3: bin_size 0.5 and a gap span of 1.5 kJ/mol pass the old
    2 * bin_size guard but leave no 2 kJ/mol density bin. Expected: the
    clear VeloxChemError, not an IndexError."""

    dp = _make_dp()
    dp.bin_size = 0.5
    dp.alpha, dp.H12 = 0.0, 1.0
    with pytest.raises(VeloxChemError, match='energy gap spans only'):
        dp.compute(_short_gap_results(1.5), 10.0, 0.0)
