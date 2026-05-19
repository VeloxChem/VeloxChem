//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center overlap driver.
//

#include "TabulaOverlapDriver.hpp"

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

// accumulated wall time of the overlap scatter, split into the direct
// matrix(r,c) writes and the transposed matrix(c,r) writes — for profiling
std::atomic<double> g_scatter_direct{0.0};
std::atomic<double> g_scatter_transpose{0.0};

/// @brief Evaluates a screened `CGtoPairBlock` into the overlap matrix — the
/// late-contraction recursion end to end: the seed ladder (a), the
/// primitive-pair contraction (b), the single-centre MD recursion (c), the
/// Cartesian-to-spherical assembly (d), and the scatter into the matrix (e).
///
/// An off-diagonal block pair also fills its transpose. When `profile` is
/// non-null, the per-phase wall times are accumulated into it.
auto
evaluate_pair_block(DenseMatrix &matrix, const GtoPairBlock &pair_block, const bool diagonal_pair, OverlapProfile *profile) -> void
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

    // (e) — scatter the spherical block into the matrix; the orbital indices
    // carry the AO component stride ([0]) and the per-pair offset ([ij+1])
    const auto &bra_orbitals = pair_block.bra_orbital_indices();
    const auto &ket_orbitals = pair_block.ket_orbital_indices();

    // direct pass — matrix(r, c)
    for (int ca = 0; ca < bra_components; ca++)
    {
        for (int cc = 0; cc < ket_components; cc++)
        {
            const auto *row = spherical.data() + static_cast<std::size_t>(ca * ket_components + cc) * stride;

            for (std::size_t ij = 0; ij < cdim; ij++)
            {
                const auto r = static_cast<std::size_t>(ca) * bra_orbitals[0] + bra_orbitals[ij + 1];
                const auto c = static_cast<std::size_t>(cc) * ket_orbitals[0] + ket_orbitals[ij + 1];

                matrix(r, c) = row[ij];
            }
        }
    }

    const auto t_scatter_direct = std::chrono::steady_clock::now();

    // transposed pass — matrix(c, r); an off-diagonal block pair fills the
    // other triangle too
    if (!diagonal_pair)
    {
        for (int ca = 0; ca < bra_components; ca++)
        {
            for (int cc = 0; cc < ket_components; cc++)
            {
                const auto *row = spherical.data() + static_cast<std::size_t>(ca * ket_components + cc) * stride;

                for (std::size_t ij = 0; ij < cdim; ij++)
                {
                    const auto r = static_cast<std::size_t>(ca) * bra_orbitals[0] + bra_orbitals[ij + 1];
                    const auto c = static_cast<std::size_t>(cc) * ket_orbitals[0] + ket_orbitals[ij + 1];

                    matrix(c, r) = row[ij];
                }
            }
        }
    }

    const auto t_scatter_done = std::chrono::steady_clock::now();

    g_scatter_direct.fetch_add(seconds(t_scatter_direct - t_transform), std::memory_order_relaxed);
    g_scatter_transpose.fetch_add(seconds(t_scatter_done - t_scatter_direct), std::memory_order_relaxed);

    if (profile != nullptr)
    {
        profile->seed      += seconds(t_seed - t_start);
        profile->contract  += seconds(t_contract - t_seed);
        profile->md        += seconds(t_md - t_contract);
        profile->transform += seconds(t_transform - t_md);
        profile->scatter   += seconds(t_scatter_done - t_transform);
    }
}

}  // namespace

auto
scatter_profile() -> ScatterProfile
{
    return {g_scatter_direct.load(std::memory_order_relaxed), g_scatter_transpose.load(std::memory_order_relaxed)};
}

auto
reset_scatter_profile() -> void
{
    g_scatter_direct.store(0.0, std::memory_order_relaxed);
    g_scatter_transpose.store(0.0, std::memory_order_relaxed);
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

    // the triangular basis-function-block pairs are the unit of parallel work
    // — each pair writes a disjoint AO region of the matrix, so the concurrent
    // writes do not race
    std::vector<std::pair<std::size_t, std::size_t>> block_pairs;
    for (std::size_t i = 0; i < gto_blocks.size(); i++)
    {
        for (std::size_t j = 0; j <= i; j++)
        {
            block_pairs.push_back({i, j});
        }
    }

    // per-thread profile accumulators, summed after the parallel region
    std::vector<OverlapProfile> thread_profiles;
    if (profile != nullptr)
    {
        thread_profiles.resize(static_cast<std::size_t>(omp_get_max_threads()));
    }

#pragma omp parallel for schedule(dynamic)
    for (int p = 0; p < static_cast<int>(block_pairs.size()); p++)
    {
        OverlapProfile *slot =
            (profile != nullptr) ? &thread_profiles[static_cast<std::size_t>(omp_get_thread_num())] : nullptr;

        const auto i = block_pairs[static_cast<std::size_t>(p)].first;
        const auto j = block_pairs[static_cast<std::size_t>(p)].second;

        const auto t_setup = std::chrono::steady_clock::now();

        // threshold <= 0 disables screening — every pair is kept, so the
        // plain pair block is built directly and the screening pre-pass (and
        // its per-pair estimator evaluation) is skipped; with screening on a
        // screened pair block drops the negligible pairs at construction
        const auto pair_block =
            (threshold > 0.0)
                ? GtoPairBlock(gto_blocks[i],
                               gto_blocks[j],
                               make_overlap_screening_estimator(gto_blocks[i].angular_momentum(),
                                                                gto_blocks[j].angular_momentum()),
                               threshold)
                : GtoPairBlock(gto_blocks[i], gto_blocks[j]);

        if (slot != nullptr)
        {
            slot->pair_setup += seconds(std::chrono::steady_clock::now() - t_setup);
        }

        evaluate_pair_block(matrix, pair_block, i == j, slot);
    }

    if (profile != nullptr)
    {
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
