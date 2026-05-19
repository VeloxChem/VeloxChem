//
//  Tabula — custom-recursion molecular-integral machinery.
//  Single-centre McMurchie–Davidson recursion.
//

#include "TabulaMDRecursion.hpp"

#include <algorithm>

namespace tabula {  // tabula namespace

namespace {  // unnamed namespace

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
compute_one_center_md(const std::vector<double> &contracted_seed,
                      const std::size_t          order,
                      const std::size_t          cdim,
                      const std::vector<double> &ac_x,
                      const std::vector<double> &ac_y,
                      const std::vector<double> &ac_z) -> std::vector<double>
{
    const auto stride = ((cdim + 7) / 8) * 8;
    const auto L      = order;

    // per-degree buffers — level[d] holds [r]^m for the degree-d monomials and
    // m = 0 … L-d: monomial_count(d) × (L-d+1) rows of `stride`
    std::vector<std::vector<double>> level(L + 1);

    for (std::size_t d = 0; d <= L; d++)
    {
        level[d].assign(monomial_count(d) * (L - d + 1) * stride, 0.0);
    }

    // level 0 — the contracted seed ladder [0]^m, m = 0 … L
    if (contracted_seed.size() >= level[0].size())
    {
        std::copy(contracted_seed.begin(), contracted_seed.begin() + static_cast<long>(level[0].size()), level[0].begin());
    }

    // build levels 1 … L
    for (std::size_t d = 1; d <= L; d++)
    {
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
    }

    // level L is [r]^0 for |r| = L — monomial_count(L) rows of `stride`
    return level[L];
}

}  // namespace tabula
