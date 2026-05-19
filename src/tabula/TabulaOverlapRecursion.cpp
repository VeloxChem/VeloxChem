//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center overlap recursion.
//

#include "TabulaOverlapRecursion.hpp"

#include <atomic>
#include <chrono>
#include <cstddef>

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

// accumulated per-section wall time of compute_overlap_seed — for profiling
std::atomic<double> g_seed_allocate{0.0};
std::atomic<double> g_seed_row0{0.0};
std::atomic<double> g_seed_ladder{0.0};

/// @brief The current steady-clock instant.
inline auto
seed_now() -> std::chrono::steady_clock::time_point
{
    return std::chrono::steady_clock::now();
}

/// @brief Adds the elapsed wall seconds of an interval to an accumulator.
inline auto
seed_add(std::atomic<double>                         &accumulator,
         const std::chrono::steady_clock::time_point &start,
         const std::chrono::steady_clock::time_point &end) -> void
{
    accumulator.fetch_add(std::chrono::duration<double>(end - start).count(), std::memory_order_relaxed);
}

}  // namespace

auto
seed_profile() -> SeedProfile
{
    return {g_seed_allocate.load(std::memory_order_relaxed),
            g_seed_row0.load(std::memory_order_relaxed),
            g_seed_ladder.load(std::memory_order_relaxed)};
}

auto
reset_seed_profile() -> void
{
    g_seed_allocate.store(0.0, std::memory_order_relaxed);
    g_seed_row0.store(0.0, std::memory_order_relaxed);
    g_seed_ladder.store(0.0, std::memory_order_relaxed);
}

auto
compute_overlap_seed(const GtoPairBlock& pair_block) -> const std::vector<double>&
{
    // angular momenta — plain ints (a structured binding cannot be captured
    // by an OpenMP simd region)
    const auto angular_momentums = pair_block.angular_momentums();
    const int  l_a               = angular_momentums.first;
    const int  l_c               = angular_momentums.second;

    const auto order = l_a + l_c;

    // primitive-pair count — cdim contracted pairs × nppairs primitive pairs;
    // no primitive-pair screening, so the whole vector is processed
    const auto cdim    = pair_block.number_of_contracted_pairs();
    const auto nppairs = static_cast<std::size_t>(pair_block.number_of_primitive_pairs());
    const auto pdim    = cdim * nppairs;

    // row stride — padded up to a multiple of 8 so each seed-order row is
    // aligned for SIMD
    const auto stride = ((pdim + 7) / 8) * 8;

    const auto t_seed_start = seed_now();

    // the seed buffer is a grow-only thread-local scratch — every value the
    // recursion reads is written before use (row 0 and each ladder row write
    // all pdim primitive pairs of their row), so it needs no zero-fill
    thread_local std::vector<double> seed;
    if (const auto need = static_cast<std::size_t>(order + 1) * stride; seed.size() < need)
    {
        seed.resize(need);
    }

    const auto t_allocated = seed_now();
    seed_add(g_seed_allocate, t_seed_start, t_allocated);

    if (pdim == 0) return seed;

    const auto* bra_exps = pair_block.bra_exponents();
    const auto* ket_exps = pair_block.ket_exponents();
    const auto* weights  = pair_block.weights();

    // row 0 — [0]^0 = weight · (1/2α)^l_a · (−1/2γ)^l_c, the weight being the
    // normalization and overlap factors folded into one with the contraction
    // weight c_i·c_j premultiplied in, so the contraction step is a plain
    // sum. The powers are applied as l_a + l_c separate SIMD passes; a
    // runtime-trip inner power loop would block vectorization of the
    // primitive-pair sweep.
    {
        auto* seed0 = seed.data();

#pragma omp simd
        for (std::size_t k = 0; k < pdim; k++)
        {
            seed0[k] = weights[k];
        }

        for (int i = 0; i < l_a; i++)
        {
#pragma omp simd
            for (std::size_t k = 0; k < pdim; k++)
            {
                seed0[k] *= 0.5 / bra_exps[k];
            }
        }

        for (int i = 0; i < l_c; i++)
        {
#pragma omp simd
            for (std::size_t k = 0; k < pdim; k++)
            {
                seed0[k] *= -0.5 / ket_exps[k];
            }
        }
    }

    const auto t_row0 = seed_now();
    seed_add(g_seed_row0, t_allocated, t_row0);

    // rows 1 … L — the seed ladder [0]^m = (−2ρ)·[0]^(m−1). ρ = αγ/(α+γ)
    // depends only on the exponents, so the row multiplier −2ρ is computed
    // once and reused for every ladder row instead of being rebuilt per row.
    if (order > 0)
    {
        thread_local std::vector<double> ladder_factor;
        if (ladder_factor.size() < pdim) ladder_factor.resize(pdim);
        auto* factor = ladder_factor.data();

#pragma omp simd
        for (std::size_t k = 0; k < pdim; k++)
        {
            const auto a = bra_exps[k];
            const auto g = ket_exps[k];

            factor[k] = -2.0 * a * g / (a + g);
        }

        for (int m = 1; m <= order; m++)
        {
            const auto* prev = seed.data() + static_cast<std::size_t>(m - 1) * stride;
            auto*       curr = seed.data() + static_cast<std::size_t>(m) * stride;

#pragma omp simd
            for (std::size_t k = 0; k < pdim; k++)
            {
                curr[k] = factor[k] * prev[k];
            }
        }
    }

    seed_add(g_seed_ladder, t_row0, seed_now());

    return seed;
}

}  // namespace tabula
