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

/// @brief The atom pairs of one group which are counted and moved by one thread.
struct TAtomPairChunk
{
    /// @brief The index of the group holding the atom pairs.
    size_t group;

    /// @brief The index of the chunk among the chunks of the group.
    size_t index;

    /// @brief The index of the first atom pair of the chunk.
    size_t first;

    /// @brief The index past the last atom pair of the chunk.
    size_t last;
};

/// @brief The atom pairs of one group which fall between two splitters and are
/// ordered by one thread.
struct TAtomPairBucket
{
    /// @brief The index of the group holding the atom pairs.
    size_t group;

    /// @brief The index of the first atom pair of the bucket.
    size_t first;

    /// @brief The index past the last atom pair of the bucket.
    size_t last;
};

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

/// @brief Orders the keys of one bucket by radix.
/// @param keys The keys of the atom pairs of the bucket.
/// @param work The scratch space of the counting passes, of the size of the keys.
/// @param npairs The number of atom pairs of the bucket.
/// @note The keys of a bucket lie between two splitters, so the leading bytes of
/// their patterns are shared by all of them and the passes consuming those bytes
/// leave the order untouched and are skipped. The bytes which differ are found
/// in one pass over the keys.
static auto
_order_bucket(TAtomPairKey *keys, TAtomPairKey *work, const size_t npairs) -> void
{
    if (npairs < 2) return;

    uint64_t diff = 0;

    for (size_t i = 1; i < npairs; i++) diff |= keys[i].first ^ keys[0].first;

    // NOTE: the counters of a pass are as many as the patterns its bits take, so
    // the bits of a pass are as many as the atom pairs of the bucket support and
    // the passes are as few as that allows.

    const auto nbits = (npairs >= (size_t{1} << 16)) ? 16 : ((npairs >= (size_t{1} << 11)) ? 11 : 8);

    const auto nbins = size_t{1} << nbits;

    const auto npasses = (64 + nbits - 1) / nbits;

    std::vector<uint32_t> counts(nbins);

    auto *source = keys;

    auto *target = work;

    for (int pass = 0; pass < npasses; pass++)
    {
        const auto shift = pass * nbits;

        if (((diff >> shift) & (nbins - 1)) == 0) continue;

        std::ranges::fill(counts, uint32_t{0});

        for (size_t i = 0; i < npairs; i++) counts[(source[i].first >> shift) & (nbins - 1)]++;

        uint32_t total = 0;

        for (size_t i = 0; i < nbins; i++)
        {
            const auto count = counts[i];

            counts[i] = total;

            total += count;
        }

        for (size_t i = 0; i < npairs; i++) target[counts[(source[i].first >> shift) & (nbins - 1)]++] = source[i];

        std::swap(source, target);
    }

    if (source != keys) std::copy(source, source + npairs, keys);
}

/// @brief Forms the key of each atom pair of a range and fills in its
/// interatomic distance.
static auto
_make_keys(const std::vector<int>       &bra_atoms,
           const std::vector<int>       &ket_atoms,
           const std::vector<TPoint<double>> &coordinates,
           const size_t                  first,
           const size_t                  last,
           std::vector<double>          &distances,
           std::vector<TAtomPairKey>    &keys) -> void
{
    for (size_t i = first; i < last; i++)
    {
        distances[i] = coordinates[bra_atoms[i]].distance(coordinates[ket_atoms[i]]);

        uint64_t bits = 0;

        std::memcpy(&bits, &distances[i], sizeof(double));

        keys[i] = {bits, static_cast<uint32_t>(i)};
    }
}

/// @brief Reorders a range of the atom pairs of a group through the permutation
/// carried by the ordered keys.
static auto
_move_pairs(const std::vector<TAtomPairKey> &keys,
            const size_t                     first,
            const size_t                     last,
            const std::vector<int>          &bra_atoms,
            const std::vector<int>          &ket_atoms,
            const std::vector<double>       &distances,
            std::vector<int>                &sorted_bra_atoms,
            std::vector<int>                &sorted_ket_atoms,
            std::vector<double>             &sorted_distances) -> void
{
    for (size_t i = first; i < last; i++)
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
    const auto ngroups = static_cast<int>(groups.size());

    if (ngroups == 0) return;

    const auto &coords = molecule.coordinates();

    const auto natoms = static_cast<int>(coords.size());

    bool valid = true;

#pragma omp parallel for schedule(dynamic) reduction(&& : valid)
    for (int i = 0; i < ngroups; i++)
    {
        valid = valid && std::ranges::all_of(groups[i]._bra_atoms, [&](const int j) { return j < natoms; }) &&
                std::ranges::all_of(groups[i]._ket_atoms, [&](const int j) { return j < natoms; });
    }

    errors::assertMsgCritical(valid, std::string("AtomBasisPairGroup.sort_by_distance: Atomic index out of range of molecule"));

    size_t npairs = 0;

    for (int i = 0; i < ngroups; i++)
    {
        groups[i]._distances.resize(groups[i]._bra_atoms.size());

        npairs += groups[i]._bra_atoms.size();
    }

    std::vector<std::vector<TAtomPairKey>> keys(groups.size()), work(groups.size());

    // NOTE: the atom pairs are ordered with the threads only when the caller is
    // not already inside a parallel region, as the nested regions do not compose.

    const auto nthreads = omp_in_parallel() ? 1 : omp::get_number_of_threads();

    size_t largest = 0;

    for (int i = 0; i < ngroups; i++) largest = std::max(largest, groups[i]._bra_atoms.size());

    // NOTE: a group is ordered by one thread unless that leaves the threads idle
    // and the groups hold enough atom pairs to repay dividing them. The threads
    // are left idle when the groups are fewer than the threads, and also when the
    // largest of them holds so many of the atom pairs that it bounds the time
    // whatever the other threads do, which is why the bound and not the number of
    // the groups is compared against the threads. The atom pairs repaying the
    // division are measured on fourteen threads.

    const auto divide = (nthreads > 1) && (largest * static_cast<size_t>(nthreads) > npairs) &&
                        (npairs > size_t{20000} * static_cast<size_t>(nthreads));

    if (!divide)
    {
#pragma omp parallel for schedule(dynamic) if (nthreads > 1)
        for (int i = 0; i < ngroups; i++)
        {
            const auto n = groups[i]._bra_atoms.size();

            if (n == 0)
            {
                groups[i]._order = pairord::distance;

                continue;
            }

            keys[i].resize(n);

            work[i].resize(n);

            _make_keys(groups[i]._bra_atoms, groups[i]._ket_atoms, coords, 0, n, groups[i]._distances, keys[i]);

            _order_keys(keys[i], work[i]);

            std::vector<int> bra_atoms(n), ket_atoms(n);

            std::vector<double> distances(n);

            _move_pairs(keys[i], 0, n, groups[i]._bra_atoms, groups[i]._ket_atoms, groups[i]._distances, bra_atoms,
                        ket_atoms, distances);

            groups[i]._bra_atoms = std::move(bra_atoms);

            groups[i]._ket_atoms = std::move(ket_atoms);

            groups[i]._distances = std::move(distances);

            groups[i]._order = pairord::distance;
        }
    }
    else
    {
        // NOTE: the groups differ widely in the number of atom pairs, and the
        // largest of them holds a third of the atom pairs, so ordering a group
        // per thread is bounded well below the number of the threads. The atom
        // pairs of a group are therefore divided into chunks, which are counted
        // and moved by the threads, and into buckets of ascending interatomic
        // distance, which are ordered by the threads. The unit of work is a chunk
        // or a bucket of a group and never a group, so the threads draw on one
        // pool of work spanning all the groups.

        // NOTE: the buckets are as many as a few per thread and no more. They are
        // drawn on by one pool of work spanning all the groups, so a few per
        // thread balance them, while more of them would deepen the search of the
        // splitters, which is done for every atom pair twice, without dividing
        // the work of ordering a bucket, which does not depend on their number.

        const auto nbuckets_total = static_cast<size_t>(8 * nthreads);

        // NOTE: the counters of a group are as many as its chunks times its
        // buckets, so the chunks are capped at a few per thread. Were they to
        // grow with the atom pairs as the buckets do, the counters would grow as
        // the square of the atom pairs of the group and outweigh the ordering.

        const auto max_chunks = static_cast<size_t>(4 * nthreads);

        std::vector<size_t> nchunks(groups.size(), 0), nbuckets(groups.size(), 0), offsets(groups.size(), 0);

        std::vector<TAtomPairChunk> chunks;

        size_t ncounts = 0;

        for (int i = 0; i < ngroups; i++)
        {
            const auto n = groups[i]._bra_atoms.size();

            offsets[i] = ncounts;

            if (n == 0) continue;

            nchunks[i] = std::min((n + size_t{4095}) / size_t{4096}, max_chunks);

            nbuckets[i] = std::max(size_t{1}, (n * nbuckets_total) / npairs);

            ncounts += nchunks[i] * nbuckets[i];

            const auto chunk_size = (n + nchunks[i] - 1) / nchunks[i];

            for (size_t j = 0; j < nchunks[i]; j++)
            {
                chunks.push_back({static_cast<size_t>(i), j, j * chunk_size, std::min(n, (j + 1) * chunk_size)});
            }
        }

        const auto nchunk_tasks = static_cast<int>(chunks.size());

        // NOTE: the keys of the groups are allocated by the threads, as the
        // groups hold enough atom pairs for the allocation to be seen in the
        // time of the ordering.

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < ngroups; i++)
        {
            keys[i].resize(groups[i]._bra_atoms.size());

            work[i].resize(groups[i]._bra_atoms.size());
        }

        // form the key of each atom pair

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < nchunk_tasks; i++)
        {
            const auto group = chunks[i].group;

            _make_keys(groups[group]._bra_atoms, groups[group]._ket_atoms, coords, chunks[i].first, chunks[i].last,
                       groups[group]._distances, keys[group]);
        }

        // NOTE: the splitters of a group are drawn from a sample of its keys, so
        // that the buckets hold a comparable number of atom pairs whatever the
        // distribution of the interatomic distances is.

        std::vector<std::vector<uint64_t>> splitters(groups.size());

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < ngroups; i++)
        {
            if (nbuckets[i] < 2) continue;

            const auto n = keys[i].size();

            const auto nsample = std::min(n, 32 * nbuckets[i]);

            std::vector<uint64_t> sample(nsample);

            for (size_t j = 0; j < nsample; j++) sample[j] = keys[i][(j * n) / nsample].first;

            std::sort(sample.begin(), sample.end());

            splitters[i].resize(nbuckets[i] - 1);

            for (size_t j = 0; j + 1 < nbuckets[i]; j++)
            {
                splitters[i][j] = sample[((j + 1) * nsample) / nbuckets[i]];
            }
        }

        // count the atom pairs of each bucket of each chunk

        std::vector<uint32_t> counts(ncounts, 0);

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < nchunk_tasks; i++)
        {
            const auto group = chunks[i].group;

            const auto &gkeys = keys[group];

            const auto &gsplitters = splitters[group];

            auto *gcounts = counts.data() + offsets[group];

            const auto nc = nchunks[group];

            for (size_t j = chunks[i].first; j < chunks[i].last; j++)
            {
                const auto bucket = static_cast<size_t>(
                    std::upper_bound(gsplitters.begin(), gsplitters.end(), gkeys[j].first) - gsplitters.begin());

                gcounts[bucket * nc + chunks[i].index]++;
            }
        }

        // NOTE: the counters of one bucket are accumulated in the order of the
        // chunks, so that the atom pairs of a bucket keep the order they have and
        // the bucket is ordered by the ordering of its keys alone.

        std::vector<std::vector<size_t>> bounds(groups.size());

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < ngroups; i++)
        {
            if (nchunks[i] == 0) continue;

            auto *gcounts = counts.data() + offsets[i];

            bounds[i].resize(nbuckets[i] + 1);

            uint32_t total = 0;

            for (size_t j = 0; j < nbuckets[i]; j++)
            {
                bounds[i][j] = total;

                for (size_t k = 0; k < nchunks[i]; k++)
                {
                    const auto count = gcounts[j * nchunks[i] + k];

                    gcounts[j * nchunks[i] + k] = total;

                    total += count;
                }
            }

            bounds[i][nbuckets[i]] = total;
        }

        // move the atom pairs of each chunk to their buckets

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < nchunk_tasks; i++)
        {
            const auto group = chunks[i].group;

            const auto &gkeys = keys[group];

            const auto &gsplitters = splitters[group];

            auto *gcounts = counts.data() + offsets[group];

            auto &gwork = work[group];

            const auto nc = nchunks[group];

            for (size_t j = chunks[i].first; j < chunks[i].last; j++)
            {
                const auto bucket = static_cast<size_t>(
                    std::upper_bound(gsplitters.begin(), gsplitters.end(), gkeys[j].first) - gsplitters.begin());

                gwork[gcounts[bucket * nc + chunks[i].index]++] = gkeys[j];
            }
        }

        // NOTE: the buckets are ordered by ascending interatomic distance and are
        // ordered independently of each other, so the keys of the group are in
        // order once the keys of every bucket are. The keys of the group are no
        // longer read once they are moved to their buckets, so they carry the
        // scratch space of the counting passes and no space is allocated for it.

        std::vector<TAtomPairBucket> buckets;

        for (int i = 0; i < ngroups; i++)
        {
            for (size_t j = 0; j < nbuckets[i]; j++)
            {
                if (bounds[i][j + 1] > bounds[i][j]) buckets.push_back({static_cast<size_t>(i), bounds[i][j], bounds[i][j + 1]});
            }
        }

        const auto nbucket_tasks = static_cast<int>(buckets.size());

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < nbucket_tasks; i++)
        {
            const auto group = buckets[i].group;

            _order_bucket(work[group].data() + buckets[i].first, keys[group].data() + buckets[i].first,
                          buckets[i].last - buckets[i].first);
        }

        for (int i = 0; i < ngroups; i++) keys[i].swap(work[i]);

        // reorder the atom pairs of each chunk

        std::vector<std::vector<int>> bra_atoms(groups.size()), ket_atoms(groups.size());

        std::vector<std::vector<double>> distances(groups.size());

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < ngroups; i++)
        {
            bra_atoms[i].resize(keys[i].size());

            ket_atoms[i].resize(keys[i].size());

            distances[i].resize(keys[i].size());
        }

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < nchunk_tasks; i++)
        {
            const auto group = chunks[i].group;

            _move_pairs(keys[group], chunks[i].first, chunks[i].last, groups[group]._bra_atoms, groups[group]._ket_atoms,
                        groups[group]._distances, bra_atoms[group], ket_atoms[group], distances[group]);
        }

        for (int i = 0; i < ngroups; i++)
        {
            groups[i]._bra_atoms = std::move(bra_atoms[i]);

            groups[i]._ket_atoms = std::move(ket_atoms[i]);

            groups[i]._distances = std::move(distances[i]);

            groups[i]._order = pairord::distance;
        }

    }
}
