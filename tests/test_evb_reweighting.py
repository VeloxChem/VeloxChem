import numpy as np

from veloxchem.evbdataprocessing import EvbDataProcessing
from veloxchem.outputstream import OutputStream


class TestEvbReweighting:
    """EVB free energies from frames that sampled other potentials than the
    analysed ones, as the MD does (it samples the integration systems).

    A 1-D model: the analysed states are harmonic, E1 = k/2 (x - a)^2 and
    E2 = k/2 (x + a)^2 + delta, so dG0 = delta exactly and the E_g profile
    has a closed form. Window lambda samples (1 - lambda) k1/2 (x - a1)^2 +
    lambda (k2/2 (x + a2)^2 + delta) instead, a Gaussian in x drawn exactly.
    reweighting="mbar", which knows every frame's energy in every sampled
    system (U_int), recovers dG0 and the barrier; "none" does not. The
    default ("auto") reweights when U_int is present.
    """

    kT = 0.0083144626 * 300.0
    k, a, delta, H12 = 145.0, 1.0, -60.0, 30.0
    k1, a1, k2, a2 = 120.0, 0.95, 150.0, 1.0

    def sampled(self, lam, x):
        return ((1 - lam) * 0.5 * self.k1 * (x - self.a1)**2 + lam *
                (0.5 * self.k2 * (x + self.a2)**2 + self.delta))

    def results(self, n_per_window=200, seed=0):
        rng = np.random.default_rng(seed)
        lambdas = np.round(np.linspace(0, 1, 31), 3)
        stiffness = (1 - lambdas) * self.k1 + lambdas * self.k2
        centre = ((1 - lambdas) * self.k1 * self.a1 -
                  lambdas * self.k2 * self.a2) / stiffness
        x = np.concatenate([
            rng.normal(c, np.sqrt(self.kT / s), n_per_window)
            for c, s in zip(centre, stiffness)
        ])
        indices = np.repeat(np.arange(len(lambdas)), n_per_window)
        return {
            "Lambda": lambdas.tolist(),
            "Lambda_frame": lambdas[indices],
            "Lambda_indices": indices,
            "configuration_results": {
                "toy": {
                    "E1_pes": 0.5 * self.k * (x - self.a)**2,
                    "E2_pes": 0.5 * self.k * (x + self.a)**2 + self.delta,
                    "U_int": self.sampled(lambdas[None, :], x[:, None]),
                    "Temp_set": 300.0,
                }
            },
        }

    def exact_profile(self):
        # xi = E1 - E2 is linear in x, so it is Gaussian on state 1, and
        # g(xi) = (xi - <xi>_1)^2 / (4 lambda_r) - xi/2 - R/2 (E_g - E1).
        lambda_r = 2 * self.k * self.a**2
        xi = np.linspace(-600, 600, 1200001)
        R = np.sqrt(xi**2 + 4 * self.H12**2)
        g = (xi + lambda_r + self.delta)**2 / (4 * lambda_r) - xi / 2 - R / 2
        top = np.argmax(np.where(np.abs(xi) < 100, g, -np.inf))
        reactant, product = np.min(g[:top]), np.min(g[top:])
        return g[top] - reactant, product - reactant

    def compute(self, results, reweighting=None):
        dp = EvbDataProcessing(ostream=OutputStream(None))
        if reweighting is not None:
            dp.reweighting = reweighting
        dp.alpha = 0.0
        dp.H12 = self.H12
        return dp.compute(results, np.nan,
                          np.nan)["configuration_results"]["toy"]

    def test_reweighting_recovers_the_analysed_states(self):
        barrier, free_energy = self.exact_profile()
        none = self.compute(self.results(), "none")
        mbar = self.compute(self.results(), "mbar")
        default = self.compute(self.results())

        assert abs(mbar["dGfep"][-1] - self.delta) < 4 * mbar["dGfep_unc"][-1]
        assert abs(mbar["analytical"]["barrier"] - barrier) < 1.5
        assert abs(mbar["binless"]["barrier"] - barrier) < 1.5
        assert abs(mbar["binless"]["free_energy"] - free_energy) < 1.5

        assert abs(none["dGfep"][-1] - self.delta) > 10
        assert abs(none["analytical"]["barrier"] - barrier) > 5

        # The default reweights whenever the cross energies are there.
        assert default["dGfep"] == mbar["dGfep"]
