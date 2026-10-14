import numpy as np
import pytest

from veloxchem.evbdataprocessing import EvbDataProcessing
from veloxchem.outputstream import OutputStream


class TestEvbFit:
    """The alpha/H12 fit of EvbDataProcessing on an exactly solvable model.

    Two harmonic states along one coordinate x, with equal force constants:
    H11 = k/2 (x - 1)^2 and H22 = k/2 (x + 1)^2 + Delta, so that the
    reorganization energy is lambda_r = 2k. In the mapping window
    V = (1 - lam) H11 + lam H22 the coordinate is Gaussian around 1 - 2 lam,
    which gives the exact window free energy
    dG_map(lam) = lambda_r lam (1 - lam) + Delta lam,
    and the exact ground-state profile along the gap xi = H11 - H22
    dg(xi) = (xi + lambda_r + Delta)^2 / (4 lambda_r) - xi/2 - R/2,
    with R = sqrt(xi^2 + 4 H12^2).

    The samples are drawn at Delta = 0, and the fit has to find the alpha that
    reproduces the barrier and reaction free energy of the model with
    Delta = alpha*. The FEP step is replaced by the exact dG_map, so the test
    checks the fit and the EVB profile alone.
    """

    lambda_r = 290.0
    temperature = 300.0

    def _samples(self, n_per_window=200, seed=0):
        lambdas = np.round(
            np.concatenate([
                np.linspace(0, 0.1, 6)[:-1],
                np.linspace(0.1, 0.9, 21)[:-1],
                np.linspace(0.9, 1, 6)
            ]), 3)
        k = self.lambda_r / 2
        kb = EvbDataProcessing(ostream=OutputStream(None)).kb
        kT = kb * self.temperature
        rng = np.random.default_rng(seed)
        lambda_frame = np.repeat(lambdas, n_per_window)
        x = rng.normal(1 - 2 * lambda_frame, np.sqrt(kT / k))
        return {
            "Lambda": list(lambdas),
            "Lambda_frame": lambda_frame,
            "Lambda_indices": np.repeat(np.arange(len(lambdas)),
                                        n_per_window),
            "configuration_results": {
                "model": {
                    "E1_pes": k / 2 * (x - 1)**2,
                    "E2_pes": k / 2 * (x + 1)**2,
                    "Temp_set": self.temperature,
                }
            },
        }

    def _exact_targets(self, delta, H12):
        """Barrier and reaction free energy of the exact profile."""

        xi = np.arange(-self.lambda_r - abs(delta) - 200,
                       self.lambda_r + abs(delta) + 200, 0.001)
        dg = ((xi + self.lambda_r + delta)**2 / (4 * self.lambda_r) - xi / 2 -
              np.sqrt(xi**2 + 4 * H12**2) / 2)
        inner = np.arange(1, len(xi) - 1)
        minima = inner[(dg[inner] < dg[inner - 1])
                       & (dg[inner] < dg[inner + 1])]
        reactant, product = minima[0], minima[-1]
        barrier = dg[reactant:product + 1].max() - dg[reactant]
        return barrier, dg[product] - dg[reactant]

    @pytest.mark.parametrize("alpha_star, H12_star", [(0.0, 50.0),
                                                      (-60.0, 30.0)])
    def test_fit_recovers_alpha_and_H12(self, alpha_star, H12_star):
        # alpha is added to the product state, so the exothermic case
        # (alpha* = -60) also pins the sign convention. The tolerances cover
        # the bias of the analytical profile on the default lambda ladder
        # (at most 0.15 kJ/mol in H12 for these cases).
        results = self._samples()
        E2 = results["configuration_results"]["model"]["E2_pes"]
        barrier, free_energy = self._exact_targets(alpha_star, H12_star)

        dp = EvbDataProcessing(ostream=OutputStream(None))
        dp.calculate_per_replica = False

        def exact_dgfep(E1, E2_shifted, Temp_set, lambda_indices=None):
            alpha = np.mean(E2_shifted - E2)
            lam = np.asarray(dp.Lambda)
            dg = self.lambda_r * lam * (1 - lam) + alpha * lam
            return dg, np.zeros_like(dg)

        dp._calculate_dGfep = exact_dgfep
        dp.compute(results, barrier, free_energy)

        assert dp.alpha == pytest.approx(alpha_star, abs=0.1)
        assert dp.H12 == pytest.approx(H12_star, abs=0.2)
