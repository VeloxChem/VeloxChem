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
"""Plot source alignment and public API checks; synthetic inputs are fixtures."""
from copy import deepcopy
from dataclasses import asdict
import json

import matplotlib.pyplot as plt
import numpy as np

from veloxchem.autocdam import FitResult, ParameterSet, plot_fit


def metal_result():
    parameters = ParameterSet({"metal_lj_sigma": 2.5, "metal_lj_epsilon": 0.3})
    history = [{
        "iteration": 1,
        "objective": 3.0
    }, {
        "iteration": 2,
        "objective": 4.0
    }]
    diagnostics = {
        "record_ids": ["b", "a", "d", "c"],
        "qm_interaction_kjmol": [-0.2, 0.4, -0.3, 0.5],
        "mm_interaction_kjmol": [-0.1, 0.6, -0.25, 0.45],
        "scan_distances_angstrom": [3.5, 2.5, 3.5, 2.5],
        "scan_directions": [[1, 0, 0], [1, 0, 0], [0, 0, 1], [0, 0, 1]],
        "energy_tolerance_kjmol": 0.1,
        "stop_reason": "optimizer_terminated"
    }
    return FitResult(parameters, history=history, diagnostics=diagnostics)


class TestAutoCdamPlotting:

    def test_cdam_stage_plot_preserves_missing_values_and_uses_normalized_errors(
            self, tmp_path):
        history = [
            {
                "iteration": 1,
                "objective": 4.0,
                "errors": {
                    "IOD": 2.0,
                    "CN": 0.0
                }
            },
            {
                "iteration": 2,
                "objective": None,
                "errors": {}
            },
            {
                "iteration": 3,
                "objective": 0.25,
                "errors": {
                    "IOD": 0.5,
                    "CN": 0.0
                }
            },
        ]
        result = FitResult(ParameterSet({"metal_lj_sigma": 2.0}),
                           history=history,
                           diagnostics={"stop_reason": "max_iterations"})
        original = deepcopy(result.history)
        figure = result.plot(output_path=tmp_path / "cdam.png")
        assert len(figure.axes) == 2
        np.testing.assert_allclose(figure.axes[0].lines[0].get_ydata(),
                                   [4.0, np.nan, 0.25])
        iod = next(
            line for line in figure.axes[1].lines if line.get_label() == "IOD")
        np.testing.assert_allclose(iod.get_ydata(), [2.0, np.nan, 0.5])
        assert figure.axes[1].get_ylabel() == "Error / target tolerance"
        assert result.history == original
        assert (tmp_path / "cdam.png").read_bytes().startswith(b"\x89PNG")
        plt.close(figure)
        source = tmp_path / "summary.json"
        source.write_text(json.dumps({"result": asdict(result)}))
        figure = plot_fit(source)
        assert len(figure.axes) == 2
        plt.close(figure)

    def test_metal_plot_matches_ids_and_keeps_directions_separate(
            self, tmp_path, monkeypatch):
        from veloxchem.autocdam import Calculator
        import veloxchem as vlx

        def unexpected_calculation(*args, **kwargs):
            raise AssertionError("plotting must not run physical calculations")

        monkeypatch.setattr(Calculator, "single_point", unexpected_calculation)
        monkeypatch.setattr(vlx.MMDriver, "compute_nonbonded",
                            unexpected_calculation)
        result = metal_result()
        dataset = {
            "metadata": {
                "study": {
                    "probe": "He"
                }
            },
            "records": [{
                "record_id": "rejected",
                "distance_angstrom": 2.0,
                "direction": [1, 0, 0],
                "accepted": False
            },]
        }
        figure = result.plot(dataset=dataset,
                             zoom_range=(2.5, 3.5),
                             show_residuals=True,
                             show_rejected=True,
                             output_path=tmp_path / "metal.png")
        assert len(figure.axes) == 3
        overview, expanded, residual = figure.axes
        qm_line = next(line for line in overview.lines
                       if line.get_label() == "QM, direction 1")
        np.testing.assert_allclose(qm_line.get_xdata(), [2.5, 3.5])
        np.testing.assert_allclose(qm_line.get_ydata(), [0.4, -0.2])
        errors = next(line for line in residual.lines
                      if line.get_label() == "MM − QM, direction 1")
        np.testing.assert_allclose(errors.get_ydata(), [0.2, 0.1])
        assert expanded.get_xlim() == (2.5, 3.5)
        assert overview.get_ylabel() == "Interaction energy (kJ/mol)"
        figure.savefig(tmp_path / "metal.pdf")
        assert (tmp_path / "metal.pdf").read_bytes().startswith(b"%PDF")
        plt.close(figure)

    def test_saved_result_file_plots_with_its_probe(self, tmp_path):
        result = metal_result()
        result.diagnostics["probe"] = "He"
        source = tmp_path / "fit_result.json"
        source.write_text(json.dumps(asdict(result)))
        figure = plot_fit(source)
        np.testing.assert_allclose(figure.axes[0].lines[0].get_xdata(),
                                   [2.5, 3.5])
        assert len(figure.axes) == 1
        assert "He" in figure._suptitle.get_text()
        plt.close(figure)
