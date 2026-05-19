import math

from veloxchem.tabulalib import tabula_boys


def boys_reference(n, x):
    """F_n(x) from an independent reference — the convergent series for
    moderate x, the asymptotic form for large x. Shares no code with the
    minimax evaluator."""
    if x < 40.0:
        # F_n(x) = e^{-x} · Σ_k (2x)^k / [(2n+1)(2n+3)…(2n+2k+1)]
        term = 1.0 / (2 * n + 1)
        total = term
        k = 0
        while term > 1.0e-19 * total:
            k += 1
            term *= (2.0 * x) / (2 * n + 2 * k + 1)
            total += term
        return math.exp(-x) * total

    # asymptotic — F_0 = (√π/2)/√x, then F_l = (2l-1)·F_{l-1}/(2x)
    f = 0.5 * math.sqrt(math.pi) / math.sqrt(x)
    for m in range(1, n + 1):
        f = (2 * m - 1) * f / (2.0 * x)
    return f


class TestTabulaBoys:
    """Tests for tabula::boys — the three-region rational-minimax Boys
    function. Validated against the convergent-series / asymptotic
    reference, with arguments straddling the region boundaries."""

    # boundaries: region A < 11.8998, region B < 28.9893, region C above
    ARGUMENTS = [0.001, 0.05, 0.5, 1.0, 3.0, 8.0, 11.8, 11.89, 11.95,
                 12.5, 20.0, 28.9, 28.98, 29.05, 60.0, 150.0, 500.0]

    def test_orders_at_zero(self):

        # F_n(0) = 1 / (2n + 1)
        values = tabula_boys(32, 0.0)
        assert len(values) == 33
        for n, value in enumerate(values):
            assert abs(value - 1.0 / (2 * n + 1)) <= 1.0e-13

    def test_matches_reference(self):

        # F_0 … F_32 against the independent reference, every region
        for x in self.ARGUMENTS:
            values = tabula_boys(32, x)
            for n, value in enumerate(values):
                reference = boys_reference(n, x)
                assert abs(value - reference) <= 1.0e-10 * abs(reference) + 1.0e-13

    def test_f0_matches_erf(self):

        # F_0(x) = (√π/2)·erf(√x)/√x
        for x in (0.05, 0.5, 3.0, 11.8, 20.0, 60.0):
            f0 = tabula_boys(0, x)[0]
            reference = 0.5 * math.sqrt(math.pi) * math.erf(math.sqrt(x)) / math.sqrt(x)
            assert abs(f0 - reference) <= 1.0e-12

    def test_single_order(self):

        # a low-order request returns just F_0 … F_order
        values = tabula_boys(2, 5.0)
        assert len(values) == 3
        for n, value in enumerate(values):
            assert abs(value - boys_reference(n, 5.0)) <= 1.0e-12
