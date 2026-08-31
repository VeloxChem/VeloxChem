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
#include <utility>
#include <vector>

#include "ErrorHandler.hpp"
#include "Molecule.hpp"

namespace {  // anonymous namespace

/// @brief The key of an atom pair, as the bit pattern of the interatomic
/// distance and the index of the atom pair in the group.
/// @note The distances are positive, and a positive double compares as its bit
/// pattern does when the pattern is read as an unsigned integer, so the atom
/// pairs are ordered by ordering the patterns. The index is the second half of
/// the key, so that the atom pairs at equal interatomic distances keep the order
/// in which they were created.
using TAtomPairKey = std::pair<uint64_t, uint32_t>;

/// @brief Orders the keys of a group by radix.
/// @param keys The keys of the atom pairs of the group.
/// @param work The scratch space of the counting passes, of the size of the keys.
/// @note The sort is stable, as the counting passes preserve the order of equal
/// keys, and gives the order the ordering of the keys as pairs gives.
static auto
_order_keys(std::vector<TAtomPairKey> &keys, std::vector<TAtomPairKey> &work) -> void
{
    const auto npairs = keys.size();

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
}

/// @brief Forms the key of each atom pair of a group and fills in its
/// interatomic distance.
/// @param bra_atoms The atoms on bra side of the atom pairs.
/// @param ket_atoms The atoms on ket side of the atom pairs.
/// @param coordinates The coordinates of the atoms of the molecule.
/// @param distances The interatomic distances of the atom pairs to fill in.
/// @param keys The keys of the atom pairs to form.
static auto
_make_keys(const std::vector<int>            &bra_atoms,
           const std::vector<int>            &ket_atoms,
           const std::vector<TPoint<double>> &coordinates,
           std::vector<double>               &distances,
           std::vector<TAtomPairKey>         &keys) -> void
{
    for (size_t i = 0; i < bra_atoms.size(); i++)
    {
        distances[i] = coordinates[bra_atoms[i]].distance(coordinates[ket_atoms[i]]);

        uint64_t bits = 0;

        std::memcpy(&bits, &distances[i], sizeof(double));

        keys[i] = {bits, static_cast<uint32_t>(i)};
    }
}

/// @brief Reorders the atom pairs of a group through the permutation carried by
/// the ordered keys.
/// @param keys The ordered keys of the atom pairs.
/// @param bra_atoms The atoms on bra side of the atom pairs.
/// @param ket_atoms The atoms on ket side of the atom pairs.
/// @param distances The interatomic distances of the atom pairs.
/// @param sorted_bra_atoms The reordered atoms on bra side.
/// @param sorted_ket_atoms The reordered atoms on ket side.
/// @param sorted_distances The reordered interatomic distances.
static auto
_move_pairs(const std::vector<TAtomPairKey> &keys,
            const std::vector<int>          &bra_atoms,
            const std::vector<int>          &ket_atoms,
            const std::vector<double>       &distances,
            std::vector<int>                &sorted_bra_atoms,
            std::vector<int>                &sorted_ket_atoms,
            std::vector<double>             &sorted_distances) -> void
{
    for (size_t i = 0; i < keys.size(); i++)
    {
        const auto index = static_cast<size_t>(keys[i].second);

        sorted_bra_atoms[i] = bra_atoms[index];

        sorted_ket_atoms[i] = ket_atoms[index];

        sorted_distances[i] = distances[index];
    }
}

}  // namespace

auto
CAtomBasisPairGroup::sort_by_distance(std::vector<CAtomBasisPairGroup> &groups, const CMolecule &molecule) -> void
{
    if (groups.empty()) return;

    const auto &coords = molecule.coordinates();

    const auto natoms = static_cast<int>(coords.size());

    bool valid = true;

    for (const auto &group : groups)
    {
        valid = valid && std::ranges::all_of(group._bra_atoms, [&](const int i) { return i < natoms; }) &&
                std::ranges::all_of(group._ket_atoms, [&](const int i) { return i < natoms; });
    }

    errors::assertMsgCritical(valid, std::string("AtomBasisPairGroup.sort_by_distance: Atomic index out of range of molecule"));

    // NOTE: the keys and the scratch space of the counting passes are carried
    // across the groups, so that the largest group alone decides how much space
    // is allocated and the groups ordered after it allocate none.

    std::vector<TAtomPairKey> keys, work;

    for (auto &group : groups)
    {
        const auto npairs = group._bra_atoms.size();

        group._distances.resize(npairs);

        group._order = pairord::distance;

        if (npairs == 0) continue;

        keys.resize(npairs);

        work.resize(npairs);

        _make_keys(group._bra_atoms, group._ket_atoms, coords, group._distances, keys);

        _order_keys(keys, work);

        std::vector<int> bra_atoms(npairs), ket_atoms(npairs);

        std::vector<double> distances(npairs);

        _move_pairs(keys, group._bra_atoms, group._ket_atoms, group._distances, bra_atoms, ket_atoms, distances);

        group._bra_atoms = std::move(bra_atoms);

        group._ket_atoms = std::move(ket_atoms);

        group._distances = std::move(distances);
    }
}
