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

from veloxchem.autocdam import GeneratedBatch, ParameterSpec, SystemData
from veloxchem.errorhandler import VeloxChemError


class TestAutoCdamModels:

    def test_system_data_can_be_instantiated(self) -> None:
        system = SystemData(
            atoms=[{
                "element": "Zn"
            }, {
                "element": "O"
            }],
            atom_types=["Zn", "O"],
            charges=[2.0, -0.5],
        )
        assert len(system.atoms) == 2

    def test_generated_batch_validates_shared_frame_shape(self) -> None:
        system = SystemData(atoms=[{"element": "A"}, {"element": "B"}])
        batch = GeneratedBatch(system,
                               [[[0, 0, 0], [1, 0, 0]], [[0, 0, 0], [2, 0, 0]]])
        assert batch.n_configurations == 2

        with pytest.raises(ValueError, match="same frame shape"):
            GeneratedBatch(system, [[[0, 0, 0], [1, 0, 0]], [[0, 0, 0]]])

    def test_parameter_spec_rejects_initial_outside_bounds(self) -> None:
        with pytest.raises(VeloxChemError, match="inside bounds"):
            ParameterSpec("sigma", initial=5.0, bounds=(0.5, 4.0))
