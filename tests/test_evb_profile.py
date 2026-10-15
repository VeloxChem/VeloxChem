import numpy as np
import pytest

from veloxchem.evbdataprocessing import EvbDataProcessing
from veloxchem.evbdriver import EvbDriver
from veloxchem.outputstream import OutputStream

from test_evb_helper import EvbTestHelper


def _default_lambda():
    """The lambda ladder EvbDriver.build_systems uses by default."""

    lam = np.linspace(0, 0.1, 6)
    lam = np.append(lam[:-1], np.linspace(0.1, 0.9, 21))
    lam = np.append(lam[:-1], np.linspace(0.9, 1, 6))
    return np.round(lam, 3)


def _processor():
    dp = EvbDataProcessing(ostream=OutputStream(None))
    dp.alpha = 0.0
    dp.calculate_per_replica = False
    return dp


@pytest.mark.parametrize('lambda_r, H12', [(290.0, 50.0), (92.0, 20.0),
                                           (290.0, 100.0)])
def test_analytical_profile_matches_marcus(lambda_r, H12):
    """The analytical EVB profile reproduces the adiabatic Marcus barrier.

    For two harmonic states of equal curvature and equal energy the mapping
    free energy is exactly dG(lambda) = lambda_r lambda (1 - lambda), and
    the ground-state barrier with a constant coupling is the
    Hwang-Aqvist-Warshel result lambda_r/4 - H12 + H12^2/lambda_r (digest
    section 2; Warshel 1992 eq 13), which is exact for a symmetric reaction.
    Feeding that mapping free energy to _dGevb_analytical must give it back.

    Tolerance 0.15 kJ/mol: the profile is read on a 10 kJ/mol grid and the
    mapping free energy is interpolated linearly between lambda points, each
    worth a few hundredths of a kJ/mol here (measured: 0.03-0.12).
    """

    lam = _default_lambda()
    dG_map = lambda_r * lam * (1 - lam)
    expected = lambda_r / 4 - H12 + H12**2 / lambda_r

    xi = np.arange(-800.0, 801.0, 10.0)
    profile = EvbDataProcessing._dGevb_analytical(dG_map, lam, H12, xi)[0]
    _, barrier, free_energy = _processor()._get_free_energies(
        profile, xi, fitting=False)[:3]

    assert barrier == pytest.approx(expected, abs=0.15)
    assert free_energy == pytest.approx(0.0, abs=0.05)


def test_zero_coupling_profile_is_finite():
    """H12 = 0 must not put NaN in the profile.

    The analytical profile divides by R = sqrt(xi^2 + 4 H12^2), which is 0
    when an uncoupled pair of states is evaluated at xi = 0. The limit is
    well defined (no shift, and the mapping potential at lambda = 1/2), and
    a single NaN would spread over the whole smoothing window.
    """

    lam = _default_lambda()
    dG_map = 290.0 * lam * (1 - lam)
    xi = np.arange(-400.0, 401.0, 10.0)
    assert 0.0 in xi

    profile = EvbDataProcessing._dGevb_analytical(dG_map, lam, 0.0, xi)[0]
    assert np.all(np.isfinite(profile))

    uncertainty = EvbDataProcessing._dGevb_analytical_uncertainty(
        np.ones_like(dG_map), lam, 0.0, xi)
    assert np.all(np.isfinite(uncertainty))


@pytest.mark.timeconsuming
def test_narrow_gap_range_does_not_crash():
    """A gap range too short for dens_max_window must still give bins.

    _calculate_coordinate_bins takes a running maximum over a histogram with
    one entry per 2 kJ/mol of observed energy gap, so a narrow gap
    distribution (or a large dens_max_window) leaves fewer entries than the
    window is wide. The window is reduced to fit instead of raising
    ValueError: negative dimensions are not allowed.
    """

    folder = EvbTestHelper.evb_data_dir()
    driver = EvbDriver(ostream=OutputStream(None))
    specific, common = driver._load_output_files(
        folder / 'evb_Sn2_vacuum_Energies.csv',
        folder / 'evb_Sn2_vacuum_Data_combined.csv',
        folder / 'evb_options.json',
    )
    results = dict(common)
    results['configuration_results'] = {'vacuum': specific}

    dp = _processor()
    dp.H12 = 50.0
    dp.dens_max_window = 100000
    out = dp.compute(results, 5, 10)

    bins = np.asarray(out['coordinate_bins'])
    assert len(bins) > 2
    assert np.all(np.isfinite(bins))
