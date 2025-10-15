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

#ifndef OpenMPFunc_hpp
#define OpenMPFunc_hpp

#include <array>
#include <cstddef>
#include <vector>
#include <algorithm>
#include <ranges>

#include "GtoBlock.hpp"
#include "GtoPairBlock.hpp"
#include "BlockedGtoPairBlock.hpp"
#include "omp.h"

namespace omp {

/// @brief Sets number of OMP threads available.
/// @param nthreads The number of OMP threads.
inline auto
set_number_of_threads(const int nthreads) -> void
{
    omp_set_num_threads(nthreads);
}

/// @brief Gets number of OMP threads available.
/// @return The number of OMP threads.
inline auto
get_number_of_threads() -> int
{
    return omp_get_max_threads();
}

/// @brief Sets static scheduling for parallel region.
inline auto
set_static_scheduler() -> void
{
    omp_set_dynamic(0);
}

/// @brief Gets maximum size of task work units.
/// @return The maximum size of task work units.
inline auto
max_block_size() -> size_t
{
    return 256;
}

/// @brief Generates work groups for OMP tasks manager.
/// @param gto_blocks The vector of basis functions blocks.
/// @return The vector of work tasks.
auto make_work_tasks(const std::vector<CGtoBlock>& gto_blocks) -> std::vector<std::array<size_t, 6>>;

/// @brief Generates work groups for OMP tasks manager.
/// @param bra_gto_blocks The vector of basis functions blocks on bra side.
/// @param ket_gto_blocks The vector of basis functions blocks on ket side.
/// @return The vector of work tasks.
auto make_work_tasks(const std::vector<CGtoBlock>& bra_gto_blocks,
                     const std::vector<CGtoBlock>& ket_gto_blocks) -> std::vector<std::array<size_t, 6>>;

/// @brief Generates work groups for OMP tasks manager.
/// @param gto_pair_blocks The vector of basis functions pair blocks.
/// @return The vector of work tasks.
auto make_diag_work_group(const std::vector<CGtoPairBlock>& gto_pair_blocks) -> std::vector<std::array<size_t, 3>>;

/// @brief Generates work groups for OMP tasks manager.
/// @param gto_pair_blocks The vector of basis functions pair blocks.
/// @param ithreshold The screening threshold of integrals.
/// @return The vector of work tasks.
auto make_work_group(const std::vector<CBlockedGtoPairBlock>& gto_pair_blocks,
                     const int                                ithreshold,
                     const int                                block_size_factor) -> std::vector<std::array<size_t, 8>>;

/// @brief Generates work groups for OMP tasks manager.
/// @param gto_pair_blocks The vector of basis functions pair blocks.
/// @param min_threshold The minimal screening threshold of integrals.
/// @param max_threshold The maximum screening threshold of integrals.
/// @return The vector of work tasks.
auto make_work_group(const std::vector<CBlockedGtoPairBlock>& gto_pair_blocks,
                     const int                                min_threshold,
                     const int                                max_threshold,
                     const int                                block_size_factor) -> std::vector<std::array<size_t, 8>>;

/// @brief Generates auxilary work groups for OMP tasks manager.
/// @param gto_blocks The vector of basis functions blocks.
/// @return The vector of work tasks.
auto make_aux_work_tasks(const std::vector<CGtoBlock>& gto_blocks) -> std::vector<std::array<size_t, 3>>;

/// @brief Generates work groups for OMP tasks manager.
/// @param bra_gto_pair_blocks The vector of bra basis functions pair blocks.
/// @param ket_gto_pair_blocks The vector of ket basis functions pair blocks.
/// @return The vector of work tasks.
auto make_bra_ket_work_group(const std::vector<CBlockedGtoPairBlock>& bra_gto_pair_blocks,
                             const std::vector<CBlockedGtoPairBlock>& ket_gto_pair_blocks,
                             const int                                ithreshold,
                             const int                                block_size_factor) -> std::vector<std::array<size_t, 8>>;

/// @brief Generates auxilary work groups for OMP tasks manager.
/// @param gto_blocks The vector of basis functions blocks.
/// @param gto_pair_blocks The vector of basis functions pair blocks.
/// @param ithreshold The screening threshold of integrals.
/// @return The vector of work tasks.
auto make_aux_work_tasks(const std::vector<CGtoBlock>&            gto_blocks,
                         const std::vector<CBlockedGtoPairBlock>& gto_pair_blocks,
                         const int                                ithreshold) -> std::vector<std::array<size_t, 3>>;

/// @brief Gets angular momentum scaling factor for SIMD width.
/// @param ang_pair The angular momentum pair.
/// @return The scaling factor for SIMD width
auto angular_momentum_scale(const std::pair<int, int>& ang_pair) -> size_t;

/// @brief Partitions vector of tasks into reduced vection using round robin scheme.
/// @param tasks The vector of task to partition.
/// @param rank The rank of requested node.
/// @param nodes The number of nodes used in partitioning.
/// @return The reduced vector of work tasks.
template<class T>
inline
auto partition_tasks(const std::vector<T>& tasks, const int rank, const int nodes) -> std::vector<T>
{
    std::vector<T> rtasks;
    
    const auto nblocks = tasks.size() / nodes;
    
    rtasks.reserve(nblocks + 1);
    
    if (nblocks > 0)
    {
        std::ranges::for_each(std::views::iota(size_t{0}, nblocks), [&](const auto i) {
            rtasks.push_back(tasks[i * nodes + rank]);
        });
    }
    
    if (const auto rblocks = tasks.size() % nodes; rank < rblocks)
    {
        rtasks.push_back(tasks[nblocks * nodes + rank]);
    }
    
    return rtasks;
}

/// @brief Gets angular momentum scaling factor for SIMD width.
/// @param natoms The number of atoms.
/// @param rank The rank of requested node.
/// @param nodes The number of nodes used in partitioning.
/// @return The reduced vector of atomic indices.
auto partition_atoms(const int natoms, const int rank, const int nodes) -> std::vector<int>;

/// @brief Generates partitioning of flat buffer.
/// @param gto_pair_blocks The vector of basis functions pair blocks.
/// @param ithreshold The screening threshold of integrals.
/// @return The vector of partitioning indices.
auto partition_flat_buffer(const std::vector<CBlockedGtoPairBlock>& gto_pair_blocks,
                           const int                                ithreshold) -> std::vector<std::array<size_t, 4>>;

/// @brief Generates masked indices for reduced flat buffer.
/// @param gto_pair_blocks The vector of basis functions pair blocks.
/// @param indices The vector of partitioned indices.
/// @return The vector of masked indices for reduced flat buffer.
auto generate_flat_buffer_mask(const std::vector<CBlockedGtoPairBlock>&  gto_pair_blocks,
                               const std::vector<std::array<size_t, 4>>& indices) -> std::vector<std::pair<size_t, size_t>>;

}  // namespace omp

#endif /* OpenMPFunc_hpp */
