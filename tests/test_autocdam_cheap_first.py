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

import pytest

from veloxchem.autocdam import (AutoCDAM, CalculationResult, CDAMCandidate,
                                FitTarget)
from veloxchem.autocdam.evaluator import Evaluator, cheap_first_reason
from veloxchem.errorhandler import VeloxChemError


def _iod(role="required"):
    return FitTarget.from_mapping("IOD", {"value": 1.90, "tolerance": 0.05,
                                          "unit": "angstrom", "role": role})


class TestTargetValueAndWeight:

    def test_a_target_without_value_has_weight_zero(self):
        hfe = FitTarget.from_mapping("HFE", {"system": "water"})
        assert hfe.reference is None and hfe.weight == 0.0
        assert not hfe.is_fitted
        fitted = FitTarget.from_mapping("HFE", {"value": -4525.0,
                                                "tolerance": 50.0})
        assert fitted.weight == 1.0 and fitted.is_fitted
        weightless = FitTarget.from_mapping("HFE", {"value": -4525.0,
                                                    "tolerance": 50.0,
                                                    "weight": 0.0})
        assert not weightless.is_fitted
        with pytest.raises(VeloxChemError, match="needs a 'value'"):
            FitTarget.from_mapping("CN", {"role": "required"})

    def test_unfitted_targets_are_not_calculated_unless_monitored(self):
        from veloxchem import OutputStream

        cdamfit = AutoCDAM(ostream=OutputStream(None))
        cdamfit.setup(CDAMCandidate("Al", 3), {
            "IOD": {"value": 1.90, "tolerance": 0.05},
            "HFE": {"tolerance": 50.0},
            "CN": {"value": 6.0, "tolerance": 0.3, "weight": 0.0},
            "SHELL_DIST_RMSD": {"role": "monitor"},
        })
        assert [t.name for t in cdamfit.search_targets] == ["IOD"]
        assert [t.name for t in cdamfit.monitor_targets] == ["SHELL_DIST_RMSD"]

    def test_a_monitor_without_value_is_reported_not_scored(self):
        shell = FitTarget.from_mapping("SHELL_DIST_RMSD", {"role": "monitor"})
        evaluation = Evaluator().evaluate(
            CalculationResult(values={"IOD": 1.92, "SHELL_DIST_RMSD": 0.03}),
            [_iod(), shell])
        assert evaluation.target_results["SHELL_DIST_RMSD"][
            "prediction"] == 0.03
        assert "SHELL_DIST_RMSD" not in evaluation.metrics
        assert evaluation.feasible


class TestCheapFirst:

    def test_reason_only_beyond_the_margin(self):
        assert cheap_first_reason({"IOD": 2.00}, [_iod()], 3.0) is None
        reason = cheap_first_reason({"IOD": 2.10}, [_iod()], 3.0)
        assert reason.startswith("cheap first: IOD misses its value by 4")

    def test_calculator_check_converts_md_units(self):
        from veloxchem.autocdam.calculator import Calculator

        # IOD is analysed in angstrom; a nanometer target is converted first
        gate = [FitTarget("IOD", 0.190, 0.005, unit="nm", role="required")]
        assert Calculator._cheap_first_check({"IOD": 1.95}, gate, 3.0) is None
        assert Calculator._cheap_first_check({"IOD": 2.30}, gate, 3.0)
        assert Calculator._cheap_first_check({"IOD": 2.30}, [], 3.0) is None

    def test_skipped_target_is_not_a_failure(self):
        hfe = FitTarget("HFE", -4525.0, 50.0, role="weighted")
        calculation = CalculationResult(
            values={"IOD": 2.30},
            metadata={"HFE": {"skipped": "cheap first: IOD misses"}})
        evaluation = Evaluator().evaluate(calculation, [_iod(), hfe])
        assert not evaluation.failed
        assert not evaluation.feasible
        assert evaluation.skipped == {"HFE": "cheap first: IOD misses"}
        assert "HFE" not in evaluation.metrics

    def test_router_skips_frameworks_after_a_failed_md_check(self, monkeypatch):
        from veloxchem.autocdam import multisystem

        class Aqueous:

            def __init__(self, iod):
                self.iod = iod

            def calculate(self, batch, parameters, targets, **kwargs):
                return CalculationResult(values={"IOD": self.iod})

        class Framework:
            calls = 0

            def bond_check(self, candidate):
                return None

            def calculate(self, candidate, targets, protocol=None):
                Framework.calls += 1
                return CalculationResult(values={"Al-O": 1.9})

        class Candidate:

            @staticmethod
            def from_generated_batch(batch):
                return None

        monkeypatch.setattr(multisystem, "CDAMCandidate", Candidate)
        distance = FitTarget("Al-O", 1.87, 0.05, role="weighted",
                             metadata={"system": "mof"})
        iod = FitTarget("IOD", 1.90, 0.05, unit="angstrom", role="required",
                        metadata={"system": "water"})
        for value, calls, skipped in ((2.30, 0, True), (1.95, 1, False)):
            Framework.calls = 0
            # the framework is listed first; the aqueous MD still runs first
            router = multisystem.MultiSystemCalculator(
                {"mof": Framework(), "water": None}, Aqueous(value))
            router.cheap_first = ([iod], 3.0)
            result = router.calculate(None, None, [distance, iod])
            assert Framework.calls == calls
            assert ("skipped" in result.metadata.get("Al-O", {})) is skipped
            assert ("Al-O" in result.values) is not skipped

    def test_driver_gates_hfe_and_frameworks_on_required_md_targets(self):
        from veloxchem import OutputStream

        cdamfit = AutoCDAM(ostream=OutputStream(None))
        cdamfit.setup(CDAMCandidate("Al", 3), {
            "IOD": {"value": 1.90, "tolerance": 0.05, "role": "required"},
            "CN": {"value": 6.0, "tolerance": 0.3},
            "HFE": {"value": -4525.0, "tolerance": 50.0},
        })
        aqueous = cdamfit.calculator.inner
        gate, margin = aqueous.cheap_first
        assert [t.name for t in gate] == ["IOD"] and margin == 3.0
        assert cdamfit._study_metadata()["cheap_first"] == {
            "targets": ["IOD"], "margin": 3.0}

        cdamfit.cheap_first = False
        cdamfit.setup(CDAMCandidate("Al", 3), {
            "IOD": {"value": 1.90, "tolerance": 0.05, "role": "required"},
            "HFE": {"value": -4525.0, "tolerance": 50.0},
        })
        assert cdamfit.calculator.inner.cheap_first is None
        assert cdamfit._study_metadata()["cheap_first"] is None


class TestQuickDirtyHfe:

    def test_preset_is_short_and_explicit_options_win(self):
        from veloxchem.autocdam.md import FEP_DEFAULTS, fep_settings

        quick = fep_settings({"quick_dirty": True})
        assert quick["num_steps"] == 250 and quick["num_snapshots"] == 10
        # the production lambda windows, only shorter sampling
        assert quick["lambdas_stage1"] == FEP_DEFAULTS["lambdas_stage1"]
        assert quick["lambdas_stage2"] == FEP_DEFAULTS["lambdas_stage2"]
        assert fep_settings({"quick_dirty": True,
                             "num_steps": 200})["num_steps"] == 200
        assert fep_settings()["num_steps"] == FEP_DEFAULTS["num_steps"]

    def test_pilot_only_and_runs_the_fep_for_every_candidate(self):
        from veloxchem import OutputStream

        targets = {
            "IOD": {"value": 1.90, "tolerance": 0.05, "role": "required"},
            "HFE": {"value": -4525.0, "tolerance": 50.0},
        }
        cdamfit = AutoCDAM(ostream=OutputStream(None))
        cdamfit.quick_dirty_hfe = True
        with pytest.raises(VeloxChemError, match="set pilot = True"):
            cdamfit.setup(CDAMCandidate("Al", 3), targets)

        cdamfit.pilot = True
        cdamfit.setup(CDAMCandidate("Al", 3), targets)
        assert cdamfit.calculator_options["fep_options"] == {
            "quick_dirty": True}
        assert cdamfit.calculator.inner.cheap_first is None
        assert cdamfit._study_metadata()["cheap_first"] is None


class TestProgressAndMbarOutput:

    def test_bar_shows_the_running_calculation_and_keeps_lines_apart(
            self, monkeypatch, capsys):
        import sys
        from veloxchem import OutputStream
        from veloxchem.autocdam import search
        from veloxchem.autocdam.models import PROGRESS_STATUS, ParameterSet

        monkeypatch.setattr(search, "_interactive", lambda: True)
        ostream = OutputStream(sys.stdout)
        progress = search.FitProgress(ostream, total=2)
        progress.start_candidate()
        PROGRESS_STATUS.set("FEP stage 1/2: lambda 4/11")
        assert progress.bar().endswith("(FEP stage 1/2: lambda 4/11)")
        ostream.print_info("a message")
        ostream.flush()
        out = capsys.readouterr().out
        # the bar is erased before the message, which starts its own line
        assert "\r" + " " * ostream.width + "\r* Info * a message" in out
        progress.finish_candidate(
            ParameterSet({"metal_lj_sigma": 2.0}),
            CalculationResult(values={"HFE": -4500.0},
                              metadata={"HFE": {"mbar": "MBAR converged"}}),
            1.0)
        out = capsys.readouterr().out
        assert "HFE: MBAR converged" in out
        assert "flush" not in ostream.__dict__
        assert PROGRESS_STATUS.text is None
        progress.close()

    def test_no_redraw_into_a_file(self, tmp_path):
        from veloxchem import OutputStream
        from veloxchem.autocdam.search import FitProgress

        ostream = OutputStream(str(tmp_path / "log"))
        assert not FitProgress(ostream, total=2)._live

    def test_pymbar_messages_become_one_conclusion(self):
        import logging
        from veloxchem.autocdam.md import _PymbarLog

        logger = logging.getLogger("pymbar.mbar_solvers")
        with _PymbarLog() as log:
            for _ in range(2):
                logger.warning("max_delta = 0.000000e+00, tol = 1.000000e-12,"
                               " maximum_iterations = 10000, iterations "
                               "completed = 9999")
                logger.warning("No solution found to within tolerance.")
                logger.warning("The solution with the smallest gradient "
                               "4.104087e-12 norm is adaptive")
                logger.info("Final gradient norm: 4.1e-12")
        usable, text = log.summary()
        assert usable
        assert text.startswith("MBAR usable: 2 of 2 solve(s) ended at "
                               "gradient 4.1e-12")
        assert logging.getLogger("pymbar").handlers == []

        with _PymbarLog() as log:
            logger.warning("max_delta = 3.000000e-04, tol = 1.000000e-12,")
            logger.warning("No solution found to within tolerance.")
            logger.info("Final gradient norm: 1e-3")
        usable, text = log.summary()
        assert not usable and "still changing by 3.0e-04" in text
