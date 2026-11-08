#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np
import pytest

from veloxchem.molecule import Molecule
from veloxchem.mmforcefieldgenerator import MMForceFieldGenerator
from veloxchem.veloxchemlib import hartree_in_kjpermol
from veloxchem.autocdam import (AutoCDAM, CalculationResult, CDAMCandidate,
                                MetalComplexLJDriver, MetalComplexLJStudy,
                                format_table, plot_fit, samples_table,
                                save_table)

# Zn LJ parameters encoded by the synthetic QM energies (angstrom, kJ/mol).
TRUE_SIGMA, TRUE_EPSILON = 2.6, 1.4
PROBE = {
    "sigma_angstrom": 2.0,
    "epsilon_kjmol": 0.3,
    "source": "illustrative He probe"
}


class ToyCalculator:
    """Analytic IOD (angstrom) and CN: IOD = 1.90 near sigma = 1.95, CN = 6 at q = 0.5."""

    def calculate(self, batch, parameters, targets, **options):
        nonbonded = batch.system.nonbonded_parameters
        sigma = nonbonded["metal_lj_sigma_angstrom"]
        epsilon = nonbonded["metal_lj_epsilon_kcal_per_mol"]
        dummy_charge = batch.system.charges[1]
        iod = 1.30 + 0.30 * sigma + 0.05 * epsilon - 0.20 * (dummy_charge - 0.5)
        cn = 6.0 - 4.0 * (dummy_charge - 0.5)**2 - 2.0 * max(0.0, sigma - 2.3)
        values = {"IOD": iod, "CN": cn}
        metadata = {"IOD": {"unit": "angstrom"}, "CN": {"unit": "count"}}
        return CalculationResult(values=values, metadata=metadata)


def toy_fit(search, output):
    """A three-parameter Al3+ CDAM fit of the toy model."""
    parameters = {
        "sigma": (1.4, 2.6),
        "epsilon": (0.02, 0.5),
        "dummy_charge": (0.2, 0.7),
        "metal_dummy_distance": 0.9
    }
    if search == "staged":
        options = {
            "method": "staged",
            "max_iterations": 4,
            "samples_per_stage": 8
        }
    else:
        options = {
            "mode": "constrained",
            "max_evaluations": 30,
            "initial_samples": 8,
            "seed": 1
        }
    targets = {
        "IOD": {
            "value": 1.90,
            "tolerance": 0.01,
            "unit": "angstrom",
            "role": "required"
        },
        "CN": {
            "value": 6.0,
            "tolerance": 0.1
        }
    }
    fit = AutoCDAM(candidate=CDAMCandidate("Al", 3),
                   parameters=parameters,
                   fit_targets=targets,
                   calculator=ToyCalculator(),
                   cache=False,
                   checkpoint_path=output / "checkpoint.json",
                   search=options)
    result = fit.fit()
    result.save(output / "result.json")
    return result


def zn_complex_and_baseline():
    """A neutral Zn-O fragment and a hand-written baseline force field (nm, kJ/mol)."""
    molecule = Molecule.read_xyz_string(
        "2\nZn-O toy fragment\nZn 0 0 0\nO 2 0 0\n")
    molecule.set_charge(0)
    molecule.set_multiplicity(1)
    ff = MMForceFieldGenerator()
    ff.ostream.mute()
    ff.molecule = molecule
    ff.atoms = {
        0: {
            "type": "Zn",
            "name": "Zn1",
            "mass": 65.38,
            "charge": 0.5,
            "sigma": 0.22,
            "epsilon": 0.8,
            "comment": "baseline"
        },
        1: {
            "type": "o",
            "name": "O1",
            "mass": 16.0,
            "charge": -0.5,
            "sigma": 0.31,
            "epsilon": 0.4,
            "comment": "baseline"
        }
    }
    ff.bonds = {(0, 1): {"equilibrium": 0.2, "force_constant": 1000.0}}
    ff.angles, ff.dihedrals, ff.impropers, ff.pairs = {}, {}, {}, {}
    ff.partial_charges = np.array([0.5, -0.5])
    return molecule, ff


class SyntheticQM:
    """Host and probe energies are constants; the host-probe energy is LJ."""

    def __init__(self, baseline):
        self.baseline = baseline

    def single_point(self, molecule, *, settings):
        labels = molecule.get_labels()
        if labels[-1] != "He":
            return -100.0
        if len(labels) == 1:
            return -2.0
        coords = molecule.get_coordinates_in_angstrom()
        energy = 0.0
        for i, atom in self.baseline.atoms.items():
            sigma = TRUE_SIGMA if i == 0 else atom["sigma"] * 10.0
            epsilon = TRUE_EPSILON if i == 0 else atom["epsilon"]
            s = (sigma + PROBE["sigma_angstrom"]) / 2
            e = np.sqrt(epsilon * PROBE["epsilon_kjmol"])
            r = np.linalg.norm(coords[-1] - coords[i])
            energy += 4 * e * ((s / r)**12 - (s / r)**6)
        return -102.0 + energy / hartree_in_kjpermol()


class TestAutoCdamToyFits:

    @pytest.mark.parametrize("search", ["bayesian", "staged"])
    def test_toy_fit_converges_with_every_search(self, tmp_path, search):
        result = toy_fit(search, tmp_path)
        assert result.converged
        assert abs(result.diagnostics["best_properties"]["IOD"] - 1.90) <= 0.01
        assert (tmp_path / "result.json").exists()

    def test_toy_fit_results_read_back(self, tmp_path):
        toy_fit("bayesian", tmp_path)
        rows = samples_table(tmp_path)
        assert rows and "objective" in rows[0]
        assert save_table(rows, tmp_path / "samples.csv").exists()
        assert "objective" in format_table(rows[:2], ["sample_id", "objective"])
        pytest.importorskip("matplotlib")
        plot_fit(tmp_path / "result.json").savefig(tmp_path / "fit_plot.png")
        assert (tmp_path / "fit_plot.png").exists()

    def test_metal_complex_fit_recovers_the_encoded_parameters(self, tmp_path):
        molecule, baseline = zn_complex_and_baseline()
        study = MetalComplexLJStudy(molecule=molecule,
                                    target_atom_index=0,
                                    probe="He",
                                    baseline_force_field=baseline,
                                    probe_parameters=PROBE,
                                    scan_distances_angstrom=np.linspace(
                                        2.4, 5.0, 10),
                                    scan_directions=[[0, 0, 1], [0, 1, 0]],
                                    fit_settings={
                                        "sigma_bounds_angstrom": [1.0, 4.0],
                                        "epsilon_bounds_kjmol": [0.01, 5.0],
                                        "energy_tolerance_kjmol": 0.05
                                    },
                                    workdir=tmp_path)
        driver = MetalComplexLJDriver(study, calculator=SyntheticQM(baseline))
        driver.compute_qm_dataset()
        result = driver.fit()
        assert result.parameters["metal_lj_sigma"] == pytest.approx(TRUE_SIGMA,
                                                                    abs=1e-3)
        assert result.parameters["metal_lj_epsilon"] == pytest.approx(
            TRUE_EPSILON, abs=1e-3)
        assert driver.export_fitted_model()
