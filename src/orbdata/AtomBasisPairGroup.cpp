//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "AtomBasisPairGroup.hpp"

#include <algorithm>
#include <numeric>
#include <ranges>

#include "ErrorHandler.hpp"
#include "Molecule.hpp"

auto
CAtomBasisPairGroup::sort_by_distance(const CMolecule &molecule) -> void
{
    const auto coords = molecule.coordinates();

    errors::assertMsgCritical(
        std::ranges::all_of(_bra_atoms, [&](const int i) { return i < static_cast<int>(coords.size()); }) &&
            std::ranges::all_of(_ket_atoms, [&](const int i) { return i < static_cast<int>(coords.size()); }),
        std::string("AtomBasisPairGroup.sort_by_distance: Atomic index out of range of molecule"));

    _distances.resize(_bra_atoms.size());

    std::ranges::for_each(std::views::iota(size_t{0}, _bra_atoms.size()),
                          [&](const auto i) { _distances[i] = coords[_bra_atoms[i]].distance(coords[_ket_atoms[i]]); });

    // NOTE: parallel vectors are reordered through a permutation of their common
    // indices. Sorting is stable, so atom pairs at equal interatomic distances
    // retain their construction order.

    std::vector<size_t> perm(_bra_atoms.size());

    std::iota(perm.begin(), perm.end(), size_t{0});

    std::ranges::stable_sort(perm, std::ranges::less{}, [&](const auto i) { return _distances[i]; });

    std::vector<int> bra_atoms(_bra_atoms.size());

    std::vector<int> ket_atoms(_ket_atoms.size());

    std::vector<double> distances(_distances.size());

    std::ranges::for_each(std::views::iota(size_t{0}, perm.size()), [&](const auto i) {
        bra_atoms[i] = _bra_atoms[perm[i]];

        ket_atoms[i] = _ket_atoms[perm[i]];

        distances[i] = _distances[perm[i]];
    });

    _bra_atoms = std::move(bra_atoms);

    _ket_atoms = std::move(ket_atoms);

    _distances = std::move(distances);

    _order = pairord::distance;
}
