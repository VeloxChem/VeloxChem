import math

import numpy as np
import pytest

from veloxchem.veloxchemlib import newints

mpmath = pytest.importorskip("mpmath")

# Region boundaries of the three-region minimax scheme (anc/constants.h).
X0 = 11.899848152108484
X1 = 28.989337738820740


def boys_reference(k, x):
    """High-precision Boys function F_k(x) = int_0^1 t^{2k} exp(-x t^2) dt."""
    mpmath.mp.dps = 50
    if x == 0.0:
        return 1.0 / (2 * k + 1)
    a = mpmath.mpf(k) + mpmath.mpf(1) / 2
    xx = mpmath.mpf(x)
    return float(mpmath.gammainc(a, 0, xx) / (2 * mpmath.power(xx, a)))


class TestNewIntsBoysFunction:

    def test_returns_requested_number_of_values(self):

        vals = newints.boys_function(1.5, 5)
        assert len(vals) == 5

    def test_value_at_zero_is_analytic(self):

        # F_k(0) = 1 / (2k + 1), within the scheme's ~5e-14 bound.
        vals = newints.boys_function(0.0, 20)
        for k, v in enumerate(vals):
            assert v == pytest.approx(1.0 / (2 * k + 1), abs=1.0e-13)

    @pytest.mark.parametrize("x", [
        1.0e-8, 0.05, 0.5, 1.0, 3.0, 7.0, 11.0,        # region A
        X0, 12.5, 18.0, 25.0, 28.5,                     # region B
        X1, 30.0, 45.0, 80.0, 200.0,                    # region C
    ])
    def test_matches_high_precision_reference(self, x):

        nvals = 33  # F_0 .. F_32
        vals = newints.boys_function(x, nvals)
        assert len(vals) == nvals
        for k in range(nvals):
            ref = boys_reference(k, x)
            assert abs(vals[k] - ref) < 1.0e-13, (x, k, vals[k], ref)

    def test_all_regions_within_target_tolerance(self):

        # Sweep all three regions and every supported top order; the scheme
        # guarantees a maximum absolute error of ~5e-14.
        xs = (list(np.linspace(1.0e-6, X0, 200))
              + list(np.linspace(X0, X1, 80))
              + list(np.linspace(X1, 80.0, 80)))
        worst = 0.0
        for nvals in range(1, 34):
            for x in xs:
                vals = newints.boys_function(float(x), nvals)
                for k in range(nvals):
                    err = abs(vals[k] - boys_reference(k, float(x)))
                    worst = max(worst, err)
        assert worst < 1.0e-13, f"worst absolute error {worst:.3e}"

    def test_top_order_then_downward_consistency(self):

        # The batch must agree with single high-order calls at the shared order.
        x = 4.0
        full = newints.boys_function(x, 33)
        for nvals in (1, 5, 12, 33):
            sub = newints.boys_function(x, nvals)
            assert sub[nvals - 1] == pytest.approx(full[nvals - 1], abs=1.0e-13)

    def test_region_boundaries_match_reference(self):

        # Each branch must match the true function right up to the boundary it
        # owns; comparing the branches to each other instead is meaningless,
        # since F_0 itself moves by ~|F_1| * dx across the boundary.
        for xb in (X0, X1):
            for x in (xb - 1.0e-9, xb + 1.0e-9):
                vals = newints.boys_function(x, 10)
                for k in range(10):
                    assert abs(vals[k] - boys_reference(k, x)) < 1.0e-13

    @pytest.mark.parametrize("nvals", [0, 34, 100])
    def test_out_of_range_order_raises(self, nvals):

        with pytest.raises(Exception):
            newints.boys_function(1.0, nvals)
