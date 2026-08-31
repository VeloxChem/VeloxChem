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

/// @brief Gets the index of the first item of a block of a partitioned range.
/// @param nitems The number of items to partition.
/// @param nblocks The number of blocks to partition the items into.
/// @param index The index of the block.
/// @return The index of the first item of the block.
/// @note The blocks hold the same number of items but for the remainder of the
/// division, which is spread one over the leading blocks, so the blocks differ
/// in size by at most one item. Passing the number of the blocks as the index
/// gives the number of the items, so the bounds of a block are the offsets of
/// the block and of the one after it.
static auto
_block_offset(const size_t nitems, const size_t nblocks, const size_t index) -> size_t
{
    const auto base = nitems / nblocks;

    const auto rest = nitems % nblocks;

    return index * base + std::min(index, rest);
}

}  // namespace

auto
CAtomBasisPairGroup::sort_by_distance(std::vector<CAtomBasisPairGroup> &groups, const CMolecule &molecule) -> void
{
    if (groups.empty()) return;

    const auto &coords = molecule.coordinates();

    const auto natoms = static_cast<int>(coords.size());

    const auto ngroups = static_cast<int>(groups.size());

    bool valid = true;

#pragma omp parallel for schedule(dynamic) reduction(&& : valid)
    for (int i = 0; i < ngroups; i++)
    {
        valid = valid && std::ranges::all_of(groups[i]._bra_atoms, [&](const int j) { return j < natoms; }) &&
                std::ranges::all_of(groups[i]._ket_atoms, [&](const int j) { return j < natoms; });
    }

    errors::assertMsgCritical(valid, std::string("AtomBasisPairGroup.sort_by_distance: Atomic index out of range of molecule"));

    // NOTE: the groups are ordered independently of each other, so one group is
    // ordered by one thread. Dynamic scheduling is used as the groups need not
    // hold the same number of atom pairs. The keys and the scratch space of the
    // counting passes are held by the thread ordering the group, so that the
    // threads share nothing while a group is ordered.

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < ngroups; i++)
    {
        auto &group = groups[i];

        const auto npairs = group._bra_atoms.size();

        group._distances.resize(npairs);

        group._order = pairord::distance;

        if (npairs == 0) continue;

        std::vector<TAtomPairKey> keys(npairs), work(npairs);

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

auto
CAtomBasisPairGroup::partition(const size_t nblocks) const -> std::vector<CAtomBasisPairGroup>
{
    errors::assertMsgCritical(nblocks > 0, std::string("AtomBasisPairGroup.partition: Number of groups is zero"));

    std::vector<CAtomBasisPairGroup> groups;

    groups.reserve(nblocks);

    const auto npairs = _bra_atoms.size();

    const auto natoms = _diagonal_atoms.size();

    // NOTE: the interatomic distances are sliced along with the atom pairs they
    // belong to, and are empty when the atom pairs are not ordered by them.

    const auto ordered = !_distances.empty();

    for (size_t i = 0; i < nblocks; i++)
    {
        const auto pfirst = _block_offset(npairs, nblocks, i);

        const auto plast = _block_offset(npairs, nblocks, i + 1);

        const auto afirst = _block_offset(natoms, nblocks, i);

        const auto alast = _block_offset(natoms, nblocks, i + 1);

        // NOTE: a group left with neither an atom pair nor an atom of a diagonal
        // atom pair describes nothing and is dropped, which happens when the
        // group holds fewer atom pairs and atoms than the groups requested.

        if ((pfirst == plast) && (afirst == alast)) continue;

        groups.push_back(CAtomBasisPairGroup(_bra_basis.get(),
                                             _ket_basis.get(),
                                             std::vector<int>(_bra_atoms.begin() + pfirst, _bra_atoms.begin() + plast),
                                             std::vector<int>(_ket_atoms.begin() + pfirst, _ket_atoms.begin() + plast),
                                             std::vector<int>(_diagonal_atoms.begin() + afirst, _diagonal_atoms.begin() + alast),
                                             ordered ? std::vector<double>(_distances.begin() + pfirst, _distances.begin() + plast)
                                                     : std::vector<double>{},
                                             _bra_index,
                                             _ket_index,
                                             _symmetric,
                                             _order));
    }

    return groups;
}
