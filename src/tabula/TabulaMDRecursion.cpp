//
//  Tabula — custom-recursion molecular-integral machinery.
//  Single-centre McMurchie–Davidson recursion.
//

#include "TabulaMDRecursion.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

// accumulated per-section wall time of compute_one_center_md — for profiling
std::atomic<double>                g_md_allocate{0.0};
std::atomic<double>                g_md_seed_copy{0.0};
std::array<std::atomic<double>, 9> g_md_build{};

/// @brief The current steady-clock instant.
inline auto
md_now() -> std::chrono::steady_clock::time_point
{
    return std::chrono::steady_clock::now();
}

/// @brief Adds the elapsed wall seconds of an interval to an accumulator.
inline auto
md_add(std::atomic<double>                         &accumulator,
       const std::chrono::steady_clock::time_point &start,
       const std::chrono::steady_clock::time_point &end) -> void
{
    accumulator.fetch_add(std::chrono::duration<double>(end - start).count(), std::memory_order_relaxed);
}

/// @brief The number of Cartesian monomials of degree d.
inline auto
monomial_count(const std::size_t d) -> std::size_t
{
    return (d + 1) * (d + 2) / 2;
}

/// @brief The canonical index of the degree-d monomial (x, y, ·) — the order
/// is x descending, then y descending (CartesianMonomial's ordering).
inline auto
monomial_index(const std::size_t x, const std::size_t y, const std::size_t d) -> std::size_t
{
    const auto dx = d - x;

    return dx * (dx + 1) / 2 + (dx - y);
}

}  // namespace

auto
md_recursion_profile() -> MDRecursionProfile
{
    MDRecursionProfile profile;
    profile.allocate  = g_md_allocate.load(std::memory_order_relaxed);
    profile.seed_copy = g_md_seed_copy.load(std::memory_order_relaxed);
    for (std::size_t d = 0; d < profile.build_by_degree.size(); d++)
    {
        profile.build_by_degree[d] = g_md_build[d].load(std::memory_order_relaxed);
    }
    return profile;
}

auto
reset_md_recursion_profile() -> void
{
    g_md_allocate.store(0.0, std::memory_order_relaxed);
    g_md_seed_copy.store(0.0, std::memory_order_relaxed);
    for (auto &accumulator : g_md_build)
    {
        accumulator.store(0.0, std::memory_order_relaxed);
    }
}

auto
compute_one_center_md(const std::vector<double> &contracted_seed,
                      const std::size_t          order,
                      const std::size_t          cdim,
                      const std::vector<double> &ac_x,
                      const std::vector<double> &ac_y,
                      const std::vector<double> &ac_z) -> std::vector<double>
{
    const auto stride = ((cdim + 7) / 8) * 8;
    const auto L      = order;

    const auto t_start = md_now();

    // per-degree buffers — level[d] holds [r]^m for the degree-d monomials and
    // m = 0 … L-d: monomial_count(d) × (L-d+1) rows of `stride`
    std::vector<std::vector<double>> level(L + 1);

    for (std::size_t d = 0; d <= L; d++)
    {
        level[d].assign(monomial_count(d) * (L - d + 1) * stride, 0.0);
    }

    const auto t_allocated = md_now();
    md_add(g_md_allocate, t_start, t_allocated);

    // level 0 — the contracted seed ladder [0]^m, m = 0 … L
    if (contracted_seed.size() >= level[0].size())
    {
        std::copy(contracted_seed.begin(), contracted_seed.begin() + static_cast<long>(level[0].size()), level[0].begin());
    }

    md_add(g_md_seed_copy, t_allocated, md_now());

    // build levels 1 … L
    for (std::size_t d = 1; d <= L; d++)
    {
        const auto t_level = md_now();

        const auto md_count   = L - d + 1;   // m rows at level d
        const auto md_count_1 = L - d + 2;   // m rows at level d-1

        std::size_t r_index = 0;

        // canonical monomial order — x descending, then y descending
        for (std::size_t x = d + 1; x-- > 0;)
        {
            for (std::size_t y = (d - x) + 1; y-- > 0;)
            {
                // first Cartesian axis with exponent ≥ 1, its exponent r_i,
                // the AC component, and the index of [r − 1_i] at level d-1
                int                 axis;
                std::size_t         r_i;
                std::size_t         idx1;
                const double       *ac;

                if (x >= 1)
                {
                    axis = 0;
                    r_i  = x;
                    ac   = ac_x.data();
                    idx1 = monomial_index(x - 1, y, d - 1);
                }
                else if (y >= 1)
                {
                    axis = 1;
                    r_i  = y;
                    ac   = ac_y.data();
                    idx1 = monomial_index(x, y - 1, d - 1);
                }
                else
                {
                    axis = 2;
                    r_i  = d - x - y;
                    ac   = ac_z.data();
                    idx1 = monomial_index(x, y, d - 1);
                }

                // index of [r − 2_i] at level d-2, when r_i ≥ 2
                std::size_t idx2 = 0;
                if (r_i >= 2)
                {
                    if (axis == 0)
                        idx2 = monomial_index(x - 2, y, d - 2);
                    else if (axis == 1)
                        idx2 = monomial_index(x, y - 2, d - 2);
                    else
                        idx2 = monomial_index(x, y, d - 2);
                }

                const auto md_count_2 = (d >= 2) ? (L - d + 3) : std::size_t{0};

                for (std::size_t m = 0; m < md_count; m++)
                {
                    auto       *out   = level[d].data() + (r_index * md_count + m) * stride;
                    const auto *prev1 = level[d - 1].data() + (idx1 * md_count_1 + (m + 1)) * stride;

                    if (r_i >= 2)
                    {
                        const auto *prev2 = level[d - 2].data() + (idx2 * md_count_2 + (m + 1)) * stride;
                        const auto  scale = static_cast<double>(r_i - 1);

#pragma omp simd
                        for (std::size_t ij = 0; ij < cdim; ij++)
                        {
                            out[ij] = ac[ij] * prev1[ij] + scale * prev2[ij];
                        }
                    }
                    else
                    {
#pragma omp simd
                        for (std::size_t ij = 0; ij < cdim; ij++)
                        {
                            out[ij] = ac[ij] * prev1[ij];
                        }
                    }
                }

                r_index++;
            }
        }

        md_add(g_md_build[d], t_level, md_now());
    }

    // level L is [r]^0 for |r| = L — monomial_count(L) rows of `stride`
    return level[L];
}

}  // namespace tabula
