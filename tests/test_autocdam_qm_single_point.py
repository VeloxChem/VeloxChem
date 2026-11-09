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
"""Exercise real VeloxChem single points, spin selection, caching, and units."""
import math

import pytest
from veloxchem import OutputStream, hartree_in_kjpermol

from veloxchem.autocdam import (
    AutoCDAM,
    Calculator,
    CDAMCandidate,
    FitTarget,
    GeneratedBatch,
    ParameterSet,
    SystemData,
)


class TestAutoCdamQmSinglePoint:

    def test_qm_plan_groups_scf_separately_from_shared_md(self, tmp_path):
        ostream = OutputStream(str(tmp_path / "plan.log"))
        try:
            plan = Calculator(ostream=ostream).plan([
                FitTarget("QM_ENERGY", 0, 1),
                FitTarget("IOD", 2, 0.1),
                FitTarget("CN", 6, 0.5),
            ])
            assert plan.jobs == {"qm": ["QM_ENERGY"], "md": ["IOD", "CN"]}
            assert plan.methods[
                "QM_ENERGY"] == "VeloxChem single-point SCF energy"
        finally:
            ostream.close()

    def test_real_al3_single_point_excludes_dummies_and_skips_openmm(
            self, tmp_path, monkeypatch):
        import veloxchem.autocdam.calculator as module

        def unexpected_openmm(*args, **kwargs):
            raise AssertionError(
                "A QM-only request must not initialize OpenMM or build solvent")

        monkeypatch.setattr(Calculator, "_initialize_platform",
                            unexpected_openmm)
        monkeypatch.setattr(module.CoordinateGenerator,
                            "create_openmm_template", unexpected_openmm)
        calls = []
        original_compute = module.ScfRestrictedDriver.compute

        def compute(self, molecule, basis, *args, **kwargs):
            calls.append((molecule.number_of_atoms(), molecule.get_charge(),
                          molecule.get_multiplicity()))
            return original_compute(self, molecule, basis, *args, **kwargs)

        monkeypatch.setattr(module.ScfRestrictedDriver, "compute", compute)
        ostream = OutputStream(str(tmp_path / "al3_qm.log"))
        try:
            cdam = AutoCDAM(ostream=ostream)
            cdam.setup(
                CDAMCandidate("Al", 3, metal_lj_sigma=2, metal_lj_epsilon=0.2),
                {"QM_ENERGY": {
                    "value": 0,
                    "tolerance": 1
                }},
                fit_parameters={
                    "metal_lj_sigma": {
                        "initial": 2,
                        "bounds": [1, 3]
                    }
                })
            result = cdam.calculate()
            energy = result.values["QM_ENERGY"]
            assert math.isfinite(energy) and energy < 0
            assert calls == [(1, 3, 1)]
            metadata = result.metadata["QM_ENERGY"]
            assert metadata["qm_basis"] == "def2-SVP"
            assert metadata["unit"] == "hartree"
            assert metadata["scf_type"] == "restricted"
            assert metadata["scf_converged"]
            assert not metadata["cache_hit"]
            assert cdam.openmm_template is None

            batch = cdam.generated_batch
            converted = cdam.calculator.calculate(
                batch, batch.parameters,
                [FitTarget("QM_ENERGY", 0, 1, unit="kJ/mol")])
            assert converted.values["QM_ENERGY"] == pytest.approx(
                energy * hartree_in_kjpermol())
            assert converted.metadata["QM_ENERGY"]["cache_hit"]
            assert calls == [(1, 3, 1)]
        finally:
            ostream.close()

    def test_real_open_shell_and_cache_invalidation_for_basis_and_charge(
            self, tmp_path):
        # Neutral hydrogen has one electron; its doublet needs unrestricted SCF.
        batch = GeneratedBatch(
            SystemData(atoms=[{
                "element": "H"
            }, {
                "element": "EP",
                "role": "dummy"
            }]),
            [[[0, 0, 0], [0.9, 0, 0]]],
        )
        target = FitTarget("QM_ENERGY", 0, 1)
        ostream = OutputStream(str(tmp_path / "hydrogen_qm.log"))
        try:
            calculator = Calculator(ostream=ostream)
            options = {
                "qm_basis": "STO-3G",
                "formal_charge": 0,
                "multiplicity": 2
            }
            first = calculator.calculate(batch,
                                         ParameterSet(), [target],
                                         qm_options=options)
            assert -0.5 <= first.values["QM_ENERGY"] < -0.4
            assert first.metadata["QM_ENERGY"]["scf_type"] == "unrestricted"
            repeated = calculator.calculate(batch,
                                            ParameterSet(), [target],
                                            qm_options=options)
            assert repeated.metadata["QM_ENERGY"]["cache_hit"]
            assert repeated.values == first.values

            other_basis = calculator.calculate(
                batch,
                ParameterSet(), [target],
                qm_options=dict(options, qm_basis="def2-SVP"))
            assert not other_basis.metadata["QM_ENERGY"]["cache_hit"]
            assert other_basis.values["QM_ENERGY"] != pytest.approx(
                first.values["QM_ENERGY"])
            other_charge = calculator.calculate(
                batch,
                ParameterSet(),
                [target],
                qm_options=dict(options, formal_charge=-1, multiplicity=1),
            )
            assert not other_charge.metadata["QM_ENERGY"]["cache_hit"]
            assert other_charge.metadata["QM_ENERGY"][
                "scf_type"] == "restricted"
            assert math.isfinite(other_charge.values["QM_ENERGY"])
            assert other_charge.values["QM_ENERGY"] != pytest.approx(
                first.values["QM_ENERGY"])
        finally:
            ostream.close()

    def test_qm_rejects_a_frame_containing_only_dummy_sites(self, tmp_path):
        batch = GeneratedBatch(SystemData(atoms=[{
            "element": "EP"
        }]), [[[0, 0, 0]]])
        ostream = OutputStream(str(tmp_path / "no_atoms.log"))
        try:
            with pytest.raises(ValueError, match="at least one real atom"):
                Calculator(ostream=ostream).calculate(
                    batch, ParameterSet(), [FitTarget("QM_ENERGY", 0, 1)])
        finally:
            ostream.close()
