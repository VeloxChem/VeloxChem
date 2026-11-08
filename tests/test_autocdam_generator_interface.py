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
"""Check that parameter proposals produce independent CDAM geometry and charges."""
import pytest

from veloxchem.autocdam import CDAMCandidate, CoordinateGenerator, ParameterSet


class TestAutoCdamGeneratorInterface:

    def test_generator_keeps_candidate_template_unchanged(self):
        candidate = CDAMCandidate("Zn",
                                  2,
                                  metal_lj_sigma=2.0,
                                  metal_lj_epsilon=0.2)
        proposals = [
            ParameterSet({
                "dummy_charge": 0.4,
                "metal_dummy_distance": 0.8
            }),
            ParameterSet({
                "dummy_charge": 0.7,
                "metal_dummy_distance": 1.0
            })
        ]
        batches = CoordinateGenerator().generate_parameter_batch(
            candidate, proposals)
        assert candidate.dummy_charge == 0.5
        assert candidate.metal_dummy_distance == 0.9
        assert batches[0].system is not batches[1].system
        assert batches[0].configurations[0][1][0] == pytest.approx(0.8)
        assert batches[1].configurations[0][1][0] == pytest.approx(1.0)
        assert all(
            sum(batch.system.charges) == pytest.approx(2.0)
            for batch in batches)
        batch_values = [batch.parameters.values for batch in batches]
        assert batch_values == [point.values for point in proposals]
