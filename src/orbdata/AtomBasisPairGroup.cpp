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
#include <iterator>
#include <utility>
#include <vector>

#include "ErrorHandler.hpp"
#include "Molecule.hpp"
#include "OpenMPFunc.hpp"

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

/// @brief The atom basis pair group a chunk of atom pairs belongs to, and the
/// atom pairs it forms.
struct TPairChunk
{
    /// @brief The index of the group among the atom basis pair groups.
    size_t group;

    /// @brief The index of the atom basis group on bra side.
    size_t bra;

    /// @brief The index of the atom basis group on ket side.
    size_t ket;

    /// @brief The index of the first atom pair of the chunk within the group.
    size_t first;

    /// @brief The index past the last atom pair of the chunk within the group.
    size_t last;
};

/// @brief Gets the index of the first atom pair of a symmetric group whose atom
/// on bra side is the atom with a given index among the atoms of the group.
/// @param row The index of the atom on bra side among the atoms of the group.
/// @param natoms The number of atoms of the group.
/// @return The index of the first atom pair.
/// @note The atom pairs of a symmetric group are the strict upper triangle of its
/// atoms, formed with the atom on bra side running slowest, so the atom pairs of
/// the atoms before a given one are as many as this returns.
static auto
_triangle_offset(const size_t row, const size_t natoms) -> size_t
{
    return row * (natoms - 1) - row * (row - 1) / 2;
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

auto
CAtomBasisPairGroup::make_pair_groups(const std::vector<CAtomBasisGroup> &groups) -> std::vector<CAtomBasisPairGroup>
{
    const auto ngroups = groups.size();

    // NOTE: the groups are created with their atom pairs sized but not formed,
    // as the number of the atom pairs of a group follows from the atoms of the
    // atom basis groups it pairs and needs no enumeration.

    std::vector<CAtomBasisPairGroup> pair_groups;

    pair_groups.reserve(ngroups * (ngroups + 1) / 2);

    std::vector<std::pair<size_t, size_t>> sides;

    sides.reserve(ngroups * (ngroups + 1) / 2);

    for (size_t i = 0; i < ngroups; i++)
    {
        const auto natoms = groups[i].number_of_atoms();

        pair_groups.push_back(CAtomBasisPairGroup(groups[i], groups[i], (natoms > 1) ? natoms * (natoms - 1) / 2 : 0, true));

        sides.emplace_back(i, i);

        for (size_t j = i + 1; j < ngroups; j++)
        {
            pair_groups.push_back(CAtomBasisPairGroup(groups[i], groups[j], natoms * groups[j].number_of_atoms(), false));

            sides.emplace_back(i, j);
        }
    }

    // NOTE: the atom pairs of the groups are formed by the threads, one chunk of
    // atom pairs at a time, so that the work divides whatever the groups hold.
    // The chunks are drawn on by one pool spanning all the groups, as the groups
    // differ widely in the number of their atom pairs.

    std::vector<TPairChunk> chunks;

    for (size_t i = 0; i < pair_groups.size(); i++)
    {
        const auto npairs = pair_groups[i]._bra_atoms.size();

        for (size_t first = 0; first < npairs; first += _pairs_per_chunk)
        {
            chunks.push_back({i, sides[i].first, sides[i].second, first, std::min(first + _pairs_per_chunk, npairs)});
        }
    }

    const auto nchunks = static_cast<int>(chunks.size());

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nchunks; i++)
    {
        auto &group = pair_groups[chunks[i].group];

        const auto &bra_atoms = groups[chunks[i].bra].atoms();

        const auto &ket_atoms = groups[chunks[i].ket].atoms();

        if (group._symmetric)
        {
            const auto natoms = bra_atoms.size();

            // NOTE: the atom on bra side of the first atom pair of the chunk is
            // found once and is then advanced with the atom pairs, as the atom
            // pairs of a chunk follow one another.

            size_t row = 0;

            while (_triangle_offset(row + 1, natoms) <= chunks[i].first) row++;

            for (size_t k = chunks[i].first; k < chunks[i].last; k++)
            {
                if (k >= _triangle_offset(row + 1, natoms)) row++;

                group._bra_atoms[k] = bra_atoms[row];

                group._ket_atoms[k] = bra_atoms[row + 1 + k - _triangle_offset(row, natoms)];
            }
        }
        else
        {
            const auto nket = ket_atoms.size();

            for (size_t k = chunks[i].first; k < chunks[i].last; k++)
            {
                group._bra_atoms[k] = bra_atoms[k / nket];

                group._ket_atoms[k] = ket_atoms[k % nket];
            }
        }
    }

    return pair_groups;
}

auto
CAtomBasisPairGroup::divide(const std::vector<CAtomBasisPairGroup> &groups, const size_t block_size) -> std::vector<CAtomBasisPairGroup>
{
    std::vector<CAtomBasisPairGroup> blocks;

    blocks.reserve(groups.size());

    for (const auto &group : groups)
    {
        const auto npairs = group.number_of_pairs();

        // NOTE: the blocks are counted by a division and a remainder rather than
        // by rounding the division up, as adding the target size to the atom
        // pairs overflows for a target size near the largest one.

        const auto nblocks = std::max(size_t{1}, npairs / block_size + ((npairs % block_size == 0) ? size_t{0} : size_t{1}));

        auto parts = group.partition(nblocks);

        blocks.insert(blocks.end(), std::make_move_iterator(parts.begin()), std::make_move_iterator(parts.end()));
    }

    return blocks;
}

auto
CAtomBasisPairGroup::make_block_size(const std::vector<CAtomBasisPairGroup> &groups,
                                     const size_t                           blocks_per_thread,
                                     const size_t                           min_block_size,
                                     const size_t                           max_block_size) -> size_t
{
    const auto nthreads = omp_in_parallel() ? 1 : omp::get_number_of_threads();

    if (nthreads < 2) return 0;

    size_t npairs = 0;

    for (const auto &group : groups) npairs += group.number_of_pairs();

    const auto size = npairs / (blocks_per_thread * static_cast<size_t>(nthreads));

    return std::min(std::max(size, min_block_size), max_block_size);
}
