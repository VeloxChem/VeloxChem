//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center overlap driver.
//

#include "TabulaOverlapDriver.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <utility>
#include <vector>

#include "omp.h"

#include "GtoBlock.hpp"
#include "GtoFunc.hpp"
#include "Point.hpp"
#include "TabulaContraction.hpp"
#include "TabulaGtoPairBlock.hpp"
#include "TabulaMDRecursion.hpp"
#include "TabulaOverlapRecursion.hpp"
#include "TabulaOverlapScreener.hpp"
#include "TabulaOverlapTransform.hpp"

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

/// @brief Wall seconds spanned by a steady-clock interval.
inline auto
seconds(const std::chrono::steady_clock::duration &interval) -> double
{
    return std::chrono::duration<double>(interval).count();
}

// per-thread load balance of the most recent compute's block-pair loop
ThreadBalance g_balance;

// accumulated wall time of the Cartesian-to-spherical transform, keyed by
// the angular-momentum pair as l_a * 5 + l_c — for profiling
std::array<std::atomic<double>, 25> g_transform{};

/// @brief Evaluates a screened `CGtoPairBlock` into the overlap matrix — the
/// late-contraction recursion end to end: the seed ladder (a), the
/// primitive-pair contraction (b), the single-centre MD recursion (c), the
/// Cartesian-to-spherical assembly (d), and the scatter into the matrix (e).
///
/// Each element is scattered once into the matrix's upper triangle; the
/// lower triangle is filled by a single `symmetrize` once every block pair
/// is done. For a `diagonal` pair block (one belonging to a block paired
/// with itself) a contracted pair and its mirror both occur, so only the
/// `r ≤ c` half is written — the mirror entry fills the rest. When `profile`
/// is non-null, the per-phase wall times are accumulated into it.
auto
evaluate_pair_block(DenseMatrix &matrix, const GtoPairBlock &pair_block, const bool diagonal, OverlapProfile *profile) -> void
{
    const auto angular_momentums = pair_block.angular_momentums();
    const auto l_a               = angular_momentums.first;
    const auto l_c               = angular_momentums.second;
    const auto order             = static_cast<std::size_t>(l_a + l_c);

    const auto cdim = pair_block.number_of_contracted_pairs();
    if (cdim == 0) return;

    const auto nppairs = static_cast<std::size_t>(pair_block.number_of_primitive_pairs());

    const auto t_start = std::chrono::steady_clock::now();

    // (a) + (b) — the contracted seed ladder [0]^m
    const auto &seed = compute_overlap_seed(pair_block);

    const auto t_seed = std::chrono::steady_clock::now();

    const auto contracted = contract_primitive_pairs(seed, order + 1, cdim, nppairs);

    const auto t_contract = std::chrono::steady_clock::now();

    // AC = A − C, per contracted pair
    const auto &bra_coords = pair_block.bra_coordinates();
    const auto &ket_coords = pair_block.ket_coordinates();

    std::vector<double> ac_x(cdim), ac_y(cdim), ac_z(cdim);
    for (std::size_t ij = 0; ij < cdim; ij++)
    {
        const auto a = bra_coords[ij].coordinates();
        const auto c = ket_coords[ij].coordinates();
        ac_x[ij]     = a[0] - c[0];
        ac_y[ij]     = a[1] - c[1];
        ac_z[ij]     = a[2] - c[2];
    }

    // (c) — the single-centre MD recursion [r]^0
    const auto rterms = compute_one_center_md(contracted, order, cdim, ac_x, ac_y, ac_z);

    const auto t_md = std::chrono::steady_clock::now();

    // (d) — the Cartesian-to-spherical assembly
    const auto stride         = ((cdim + 7) / 8) * 8;
    const auto bra_components = 2 * l_a + 1;
    const auto ket_components = 2 * l_c + 1;

    std::vector<double> spherical(static_cast<std::size_t>(bra_components * ket_components) * stride, 0.0);
    overlap_transform(l_a, l_c, rterms.data(), cdim, spherical.data());

    const auto t_transform = std::chrono::steady_clock::now();

    g_transform[static_cast<std::size_t>(l_a * 5 + l_c)].fetch_add(seconds(t_transform - t_md), std::memory_order_relaxed);

    // (e) — scatter the spherical block into the matrix; the orbital indices
    // carry the AO component stride ([0]) and the per-pair offset ([ij+1]).
    // Each element is written once into the upper triangle — at (min(r,c),
    // max(r,c)) — and the lower triangle is filled later by `symmetrize`
    const auto &bra_orbitals = pair_block.bra_orbital_indices();
    const auto &ket_orbitals = pair_block.ket_orbital_indices();

    for (int ca = 0; ca < bra_components; ca++)
    {
        for (int cc = 0; cc < ket_components; cc++)
        {
            const auto *row = spherical.data() + static_cast<std::size_t>(ca * ket_components + cc) * stride;

            for (std::size_t ij = 0; ij < cdim; ij++)
            {
                const auto r = static_cast<std::size_t>(ca) * bra_orbitals[0] + bra_orbitals[ij + 1];
                const auto c = static_cast<std::size_t>(cc) * ket_orbitals[0] + ket_orbitals[ij + 1];

                // a diagonal pair block also produces every mirror (j,i), so
                // the r > c half is skipped — its mirror writes that slot
                if (diagonal && r > c) continue;

                matrix(std::min(r, c), std::max(r, c)) = row[ij];
            }
        }
    }

    if (profile != nullptr)
    {
        const auto t_scatter = std::chrono::steady_clock::now();
        profile->seed      += seconds(t_seed - t_start);
        profile->contract  += seconds(t_contract - t_seed);
        profile->md        += seconds(t_md - t_contract);
        profile->transform += seconds(t_transform - t_md);
        profile->scatter   += seconds(t_scatter - t_transform);
    }
}

}  // namespace

auto
overlap_thread_balance() -> ThreadBalance
{
    return g_balance;
}

auto
transform_profile() -> std::array<double, 25>
{
    std::array<double, 25> profile{};
    for (std::size_t k = 0; k < 25; k++)
    {
        profile[k] = g_transform[k].load(std::memory_order_relaxed);
    }
    return profile;
}

auto
reset_transform_profile() -> void
{
    for (auto &accumulator : g_transform)
    {
        accumulator.store(0.0, std::memory_order_relaxed);
    }
}

auto
OverlapDriver::compute(const CMolecule       &molecule,
                       const CMolecularBasis &basis,
                       const double           threshold,
                       OverlapProfile        *profile) const -> DenseMatrix
{
    const auto t_blocks = std::chrono::steady_clock::now();

    const auto gto_blocks = gtofunc::make_gto_blocks(basis, molecule);

    if (profile != nullptr)
    {
        profile->make_blocks += seconds(std::chrono::steady_clock::now() - t_blocks);
    }

    const auto dimension = static_cast<std::size_t>(gtofunc::getNumberOfAtomicOrbitals(gto_blocks));

    DenseMatrix matrix(dimension, dimension, Symmetry::symmetric);

    // the unit of parallel work is a task — one block pair restricted to a
    // sub-range of the bra block's contracted GTOs. The triangular block
    // pairs span orders of magnitude in cost and are far too few and uneven
    // to balance; splitting the bra dimension so each task holds roughly
    // `target_pairs` contracted pairs turns dozens of coarse units into
    // thousands of even ones. Each task writes a disjoint AO region, so the
    // concurrent writes do not race.
    constexpr int target_pairs = 1 << 16;

    struct Task
    {
        std::size_t i;
        std::size_t j;
        int         bra_begin;
        int         bra_end;
        bool        diagonal;
        double      cost;
    };

    std::vector<Task> tasks;
    for (std::size_t i = 0; i < gto_blocks.size(); i++)
    {
        for (std::size_t j = 0; j <= i; j++)
        {
            const int bcgtos = gto_blocks[i].number_of_basis_functions();
            const int kcgtos = gto_blocks[j].number_of_basis_functions();
            const int order  = gto_blocks[i].angular_momentum() + gto_blocks[j].angular_momentum();

            // screening keeps a block pair whole; otherwise split the bra
            // dimension into chunks of about `target_pairs` contracted pairs
            const int chunk = (threshold > 0.0) ? bcgtos : std::max(1, target_pairs / std::max(1, kcgtos));

            for (int a = 0; a < bcgtos; a += chunk)
            {
                const int    b    = std::min(a + chunk, bcgtos);
                const double cost = static_cast<double>(b - a) * kcgtos * (order + 1) * (order + 2);
                tasks.push_back({i, j, a, b, i == j, cost});
            }
        }
    }

    // longest task first — dynamic scheduling then packs the small tasks into
    // the tail rather than stalling on a late big one
    std::sort(tasks.begin(), tasks.end(), [](const Task &lhs, const Task &rhs) { return lhs.cost > rhs.cost; });

    const auto nthreads = static_cast<std::size_t>(omp_get_max_threads());

    // per-thread profile accumulators, summed after the parallel region
    std::vector<OverlapProfile> thread_profiles;
    if (profile != nullptr)
    {
        thread_profiles.resize(nthreads);
    }

    // per-thread load-balance capture of the task loop
    std::vector<double> thread_busy(nthreads, 0.0);
    std::vector<long>   thread_pairs(nthreads, 0);

    const auto t_region = std::chrono::steady_clock::now();

#pragma omp parallel for schedule(dynamic)
    for (int p = 0; p < static_cast<int>(tasks.size()); p++)
    {
        const auto tid    = static_cast<std::size_t>(omp_get_thread_num());
        const auto t_body = std::chrono::steady_clock::now();

        OverlapProfile *slot = (profile != nullptr) ? &thread_profiles[tid] : nullptr;

        const auto &task = tasks[static_cast<std::size_t>(p)];

        const auto t_setup = std::chrono::steady_clock::now();

        // threshold <= 0 disables screening — the task's bra range is built
        // directly; with screening on the whole block pair is built and the
        // estimator drops the negligible contracted-GTO pairs
        const auto pair_block =
            (threshold > 0.0)
                ? GtoPairBlock(gto_blocks[task.i],
                               gto_blocks[task.j],
                               make_overlap_screening_estimator(gto_blocks[task.i].angular_momentum(),
                                                                gto_blocks[task.j].angular_momentum()),
                               threshold)
                : GtoPairBlock(gto_blocks[task.i], gto_blocks[task.j], task.bra_begin, task.bra_end);

        if (slot != nullptr)
        {
            slot->pair_setup += seconds(std::chrono::steady_clock::now() - t_setup);
        }

        evaluate_pair_block(matrix, pair_block, task.diagonal, slot);

        thread_busy[tid] += seconds(std::chrono::steady_clock::now() - t_body);
        thread_pairs[tid] += 1;
    }

    g_balance.wall  = seconds(std::chrono::steady_clock::now() - t_region);
    g_balance.busy  = thread_busy;
    g_balance.pairs = thread_pairs;

    // the scatter filled the upper triangle only — mirror it into the lower
    const auto t_symmetrize = std::chrono::steady_clock::now();
    matrix.symmetrize();

    if (profile != nullptr)
    {
        profile->symmetrize += seconds(std::chrono::steady_clock::now() - t_symmetrize);

        for (const auto &thread_profile : thread_profiles)
        {
            profile->pair_setup += thread_profile.pair_setup;
            profile->seed       += thread_profile.seed;
            profile->contract   += thread_profile.contract;
            profile->md         += thread_profile.md;
            profile->transform  += thread_profile.transform;
            profile->scatter    += thread_profile.scatter;
        }
    }

    return matrix;
}

}  // namespace tabula
