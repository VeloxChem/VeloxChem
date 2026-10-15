import numpy as np
import pytest

from veloxchem.evbdataprocessing import EvbDataProcessing
from veloxchem.outputstream import OutputStream


class TestEvbReorganization:
    """The reorganization energy and diabatic curves of EvbDataProcessing.

    Two harmonic states along one coordinate x, with equal force constants:
    H11 = k/2 (x - a)^2 and H22 = k/2 (x + a)^2 + Delta. The reorganization
    energy is then lambda_r = 2 k a^2, the reaction free energy is Delta,
    and linear response holds exactly, so the energy-gap estimates from the
    reactant and the product side agree. In the mapping window
    V = (1 - lam) H11 + lam H22 the coordinate is Gaussian around
    a (1 - 2 lam) with variance kT/k, so the windows are sampled exactly.
    """

    k = 145.0
    a = 1.0
    delta = -60.0
    H12 = 30.0
    temperature = 300.0

    def _results(self, n_per_window=1000, seed=0):
        lambdas = np.round(
            np.concatenate([
                np.linspace(0, 0.1, 6)[:-1],
                np.linspace(0.1, 0.9, 21)[:-1],
                np.linspace(0.9, 1, 6)
            ]), 3)
        kb = EvbDataProcessing(ostream=OutputStream(None)).kb
        kT = kb * self.temperature
        rng = np.random.default_rng(seed)
        lambda_frame = np.repeat(lambdas, n_per_window)
        x = rng.normal(self.a * (1 - 2 * lambda_frame), np.sqrt(kT / self.k))
        return {
            "Lambda": list(lambdas),
            "Lambda_frame": lambda_frame,
            "Lambda_indices": np.repeat(np.arange(len(lambdas)),
                                        n_per_window),
            "configuration_results": {
                "model": {
                    "E1_pes": self.k / 2 * (x - self.a)**2,
                    "E2_pes": self.k / 2 * (x + self.a)**2 + self.delta,
                    "Temp_set": self.temperature,
                }
            },
        }

    def test_reorganization_energy_and_diabatic_curves(self):
        results = self._results()
        dp = EvbDataProcessing(ostream=OutputStream(None))
        dp.alpha = 0.0
        dp.H12 = self.H12
        dp.compute(results, 0.0, 0.0)
        result = results["configuration_results"]["model"]
        reorganization = result["reorganization"]

        lambda_r = 2 * self.k * self.a**2
        for key in ("lambda_r", "lambda_reactant", "lambda_product"):
            assert reorganization[key] == pytest.approx(
                lambda_r, abs=4 * reorganization[key + "_unc"])
        assert result["dGfep"][-1] == pytest.approx(
            self.delta, abs=4 * result["dGfep_unc"][-1])
        assert reorganization["split"] < reorganization["split_tolerance"]

        # dg2 - dg1 = -xi frame by frame, so in every sampled bin it lies
        # within half a bin of minus the bin center. The first and last
        # entries collect the frames outside the bins and are skipped.
        bins = np.asarray(results["coordinate_bins"])
        centers = 0.5 * (bins[1:] + bins[:-1])
        difference = (np.asarray(reorganization["dg2"]) -
                      np.asarray(reorganization["dg1"]))[1:-1]
        sampled = np.isfinite(difference)
        assert np.count_nonzero(sampled) > 10
        assert np.all(
            np.abs(difference[sampled] + centers[sampled]) <=
            0.5 * np.diff(bins)[sampled] + 1e-9)
