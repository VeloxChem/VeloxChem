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
#include <cstdint>
#include <cstring>
#include <ranges>
#include <vector>

#include "ErrorHandler.hpp"
#include "Molecule.hpp"

auto
CAtomBasisPairGroup::sort_by_distance(const CMolecule &molecule) -> void
{
    const auto &coords = molecule.coordinates();

    errors::assertMsgCritical(
        std::ranges::all_of(_bra_atoms, [&](const int i) { return i < static_cast<int>(coords.size()); }) &&
            std::ranges::all_of(_ket_atoms, [&](const int i) { return i < static_cast<int>(coords.size()); }),
        std::string("AtomBasisPairGroup.sort_by_distance: Atomic index out of range of molecule"));

    _distances.resize(_bra_atoms.size());

    std::ranges::for_each(std::views::iota(size_t{0}, _bra_atoms.size()),
                          [&](const auto i) { _distances[i] = coords[_bra_atoms[i]].distance(coords[_ket_atoms[i]]); });

    // NOTE: parallel vectors are reordered through a permutation of their common
    // indices, which is formed by sorting the interatomic distances.

    // NOTE: the distances are positive, and a positive double compares as its bit
    // pattern does when the pattern is read as an unsigned integer, so the order
    // is found by sorting the patterns by radix instead of by comparison. The
    // order is the same one a comparison sort gives and not an approximation of
    // it. The sort is stable, as the counting passes preserve the order of equal
    // keys, so atom pairs at equal interatomic distances retain their
    // construction order.

    const auto npairs = _distances.size();

    std::vector<std::pair<uint64_t, uint32_t>> keys(npairs), work(npairs);

    for (size_t i = 0; i < npairs; i++)
    {
        uint64_t bits = 0;

        std::memcpy(&bits, &_distances[i], sizeof(double));

        keys[i] = {bits, static_cast<uint32_t>(i)};
    }

    // NOTE: the sixty four bits of a pattern are consumed in four passes of
    // sixteen bits, so the counters of one pass fit in the cache.

    constexpr auto nbits = 16;

    constexpr auto nbins = size_t{1} << nbits;

    std::vector<uint32_t> counts(nbins);

    for (int pass = 0; pass < 4; pass++)
    {
        const auto shift = pass * nbits;

        std::ranges::fill(counts, 0);

        for (size_t i = 0; i < npairs; i++) counts[(keys[i].first >> shift) & (nbins - 1)]++;

        uint32_t total = 0;

        for (size_t i = 0; i < nbins; i++)
        {
            const auto count = counts[i];

            counts[i] = total;

            total += count;
        }

        for (size_t i = 0; i < npairs; i++) work[counts[(keys[i].first >> shift) & (nbins - 1)]++] = keys[i];

        keys.swap(work);
    }

    std::vector<int> bra_atoms(npairs);

    std::vector<int> ket_atoms(npairs);

    std::vector<double> distances(npairs);

    for (size_t i = 0; i < npairs; i++)
    {
        const auto index = static_cast<size_t>(keys[i].second);

        bra_atoms[i] = _bra_atoms[index];

        ket_atoms[i] = _ket_atoms[index];

        distances[i] = _distances[index];
    }

    _bra_atoms = std::move(bra_atoms);

    _ket_atoms = std::move(ket_atoms);

    _distances = std::move(distances);

    _order = pairord::distance;
}
