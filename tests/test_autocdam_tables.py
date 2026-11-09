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
"""Flat per-candidate tables from a fit, a driver, or a saved run."""
import csv

from veloxchem.autocdam import AutoCDAM, CalculationResult, CDAMCandidate, samples_table, save_table


class Surface:

    def calculate(self, batch, parameters, targets, **kwargs):
        sigma = batch.system.nonbonded_parameters["metal_lj_sigma_angstrom"]
        values = {"IOD": 0.9 + 0.5 * sigma, "HFE": -4000.0 - 100 * sigma}
        metadata = {"HFE": {"uncertainty": 12.0}, "IOD": {"spread": 0.03}}
        return CalculationResult(
            values={t.name: values[t.name] for t in targets},
            metadata={t.name: metadata[t.name] for t in targets})


def fit(tmp_path):
    targets = {
        "IOD": {
            "value": 1.9,
            "tolerance": 0.002
        },
        "HFE": {
            "value": -4200,
            "tolerance": 50,
            "weight": 0.1
        }
    }
    driver = AutoCDAM()
    driver.calculator = Surface()
    driver.update_settings({
        "search_method": "bayesian",
        "max_evaluations": 6,
        "seed": 1,
        "checkpoint_file": str(tmp_path / "checkpoint.json"),
    })
    return driver, driver.fit(CDAMCandidate("Al", 3),
                              targets,
                              parameters={"sigma": (1.4, 2.6)})


class TestAutoCdamTables:

    def test_rows_from_result_driver_and_saved_run_agree(self, tmp_path):
        driver, result = fit(tmp_path)
        result.save(tmp_path / "result.json")
        from_result = samples_table(result)
        from_dir = samples_table(tmp_path)
        assert len(from_result) == len(driver.fitting_driver.table.samples) == 6
        assert from_result == from_dir
        row = from_result[0]
        assert row["step"] == 1 and "parameter.metal_lj_sigma" in row
        assert row[
            "property.IOD"] == 0.9 + 0.5 * row["parameter.metal_lj_sigma"]
        assert row["uncertainty.HFE"] == 12.0 and row["spread.IOD"] == 0.03
        assert "objective" in row and "feasible" in row
        # Without scores (driver/table) the measured columns are the same.
        plain = samples_table(driver.fitting_driver)[0]
        assert {
            k: v for k, v in row.items() if k not in ("objective", "feasible")
        } == plain

    def test_saved_result_alone_and_csv(self, tmp_path):
        _, result = fit(tmp_path)
        path = result.save(tmp_path / "only" / "result.json")
        rows = samples_table(path)
        assert len(rows) == 6 and rows[0]["objective"] is not None
        with save_table(rows, tmp_path / "samples.csv").open() as handle:
            read = list(csv.DictReader(handle))
        assert len(read) == 6 and float(
            read[0]["property.HFE"]) == rows[0]["property.HFE"]
