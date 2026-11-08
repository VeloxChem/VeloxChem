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
"""Regression for the shell-cutoff error found in the Al pilot diagnostics."""
import numpy as np
import openmm.app as app
import pytest

from veloxchem.autocdam.geometry import analyze_trajectory


def analyze(radii):
    topology = app.Topology()
    chain = topology.addChain()
    topology.addAtom("MG", app.element.magnesium,
                     topology.addResidue("CDM", chain))
    for _ in radii:
        topology.addAtom("O", app.element.oxygen,
                         topology.addResidue("HOH", chain))
    # Separated shells with noisy internal extrema, independent of fit targets.
    rng = np.random.default_rng(5)
    frames = np.zeros((50, len(radii) + 1, 3))
    frames[:, 1:, 0] = radii + rng.normal(0, 0.0004, (50, len(radii)))
    boxes = np.repeat((np.eye(3) * 3)[None, :, :], len(frames), axis=0)
    return analyze_trajectory(frames, boxes, topology, ["IOD", "CN"])


class TestAutoCdamTrajectoryShell:

    @pytest.mark.parametrize("n", [4, 6, 8])
    def test_counts_whole_shell_without_assuming_six(self, n):
        first = np.linspace(0.166, 0.175, n)
        result = analyze(np.r_[first, [0.337, 0.345, 0.353]])
        assert result["CN"] == n
        cutoff = result["_analysis"]["cn_cutoff_angstrom"]
        assert first.max() * 10 < cutoff < 3.37
        assert 1.65 < result["IOD"] < 1.76

    def test_missing_outer_shell_is_rejected_instead_of_inventing_cutoff(self):
        with pytest.raises(ValueError,
                           match="No resolved first-shell RDF minimum"):
            analyze(np.linspace(0.166, 0.175, 6))

    @pytest.mark.parametrize("metal", ["Dy", "Mg"])
    def test_metals_starting_with_d_are_not_dummy_sites(self, metal):
        from veloxchem.autocdam.geometry import cdam_indices
        topology = app.Topology()
        residue = topology.addResidue("CDM", topology.addChain())
        topology.addAtom(metal.upper(), app.element.get_by_symbol(metal),
                         residue)
        for index in range(1, 4):
            topology.addAtom(f"D{index}", None, residue)
        assert cdam_indices(topology) == (0, [1, 2, 3])
