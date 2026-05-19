//
//  Tabula — custom-recursion molecular-integral machinery.
//  Single-centre McMurchie–Davidson recursion.
//

#include "TabulaMDRecursion.hpp"

#include <algorithm>
#include <array>

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

    // level L — [r]^0 for the degree-L monomials — is the returned result
    std::vector<double> result(monomial_count(L) * stride, 0.0);

    if (cdim == 0) return result;

    // The recursion is independent per contracted pair, so it is run in
    // tiles of `tile` pairs: the per-tile intermediate buffers stay
    // cache-resident, instead of the tens-to-hundreds of MB a full sweep
    // streams through DRAM as it writes and re-reads every [r]^m.
    constexpr std::size_t tile = 512;

    // the intermediate levels 0 … L-1 of a tile share a grow-only
    // thread-local scratch; every value the recursion reads is written
    // before use, so the scratch needs no zero-fill
    thread_local std::vector<double> md_scratch;

    // level d holds [r]^m for the degree-d monomials and m = 0 … L-d:
    // monomial_count(d) × (L-d+1) rows of `tile`
    std::array<std::size_t, 9> level_offset{};
    std::size_t                scratch_size = 0;
    for (std::size_t d = 0; d < L; d++)
    {
        level_offset[d] = scratch_size;
        scratch_size += monomial_count(d) * (L - d + 1) * tile;
    }
    if (md_scratch.size() < scratch_size) md_scratch.resize(scratch_size);

    auto *scratch = md_scratch.data();

    for (std::size_t t0 = 0; t0 < cdim; t0 += tile)
    {
        const auto width = std::min(tile, cdim - t0);

        // the row of degree-d monomial `idx` at seed order `m`: the
        // intermediate levels 0 … L-1 live tile-local in the scratch, level L
        // directly in the result (which keeps the full `stride`)
        const auto row = [&](const std::size_t d, const std::size_t idx, const std::size_t m) -> double * {
            if (d == L) return result.data() + idx * stride + t0;
            return scratch + level_offset[d] + (idx * (L - d + 1) + m) * tile;
        };

        // level 0 — the contracted seed ladder [0]^m, m = 0 … L, sliced to
        // this tile
        for (std::size_t m = 0; m <= L; m++)
        {
            const auto *src = contracted_seed.data() + m * stride + t0;
            std::copy(src, src + width, row(0, 0, m));
        }

        // build levels 1 … L for this tile
        for (std::size_t d = 1; d <= L; d++)
        {
            const auto  md_count = L - d + 1;   // m rows at level d
            std::size_t r_index  = 0;

            // canonical monomial order — x descending, then y descending
            for (std::size_t x = d + 1; x-- > 0;)
            {
                for (std::size_t y = (d - x) + 1; y-- > 0;)
                {
                    // first Cartesian axis with exponent ≥ 1, its exponent
                    // r_i, the tile-sliced AC component, and the index of
                    // [r − 1_i] at level d-1
                    int           axis;
                    std::size_t   r_i;
                    std::size_t   idx1;
                    const double *ac;

                    if (x >= 1)
                    {
                        axis = 0;
                        r_i  = x;
                        ac   = ac_x.data() + t0;
                        idx1 = monomial_index(x - 1, y, d - 1);
                    }
                    else if (y >= 1)
                    {
                        axis = 1;
                        r_i  = y;
                        ac   = ac_y.data() + t0;
                        idx1 = monomial_index(x, y - 1, d - 1);
                    }
                    else
                    {
                        axis = 2;
                        r_i  = d - x - y;
                        ac   = ac_z.data() + t0;
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

                    for (std::size_t m = 0; m < md_count; m++)
                    {
                        auto       *out   = row(d, r_index, m);
                        const auto *prev1 = row(d - 1, idx1, m + 1);

                        if (r_i >= 2)
                        {
                            const auto *prev2 = row(d - 2, idx2, m + 1);
                            const auto  scale = static_cast<double>(r_i - 1);

#pragma omp simd
                            for (std::size_t k = 0; k < width; k++)
                            {
                                out[k] = ac[k] * prev1[k] + scale * prev2[k];
                            }
                        }
                        else
                        {
#pragma omp simd
                            for (std::size_t k = 0; k < width; k++)
                            {
                                out[k] = ac[k] * prev1[k];
                            }
                        }
                    }

                    r_index++;
                }
            }
        }
    }

    return result;
}

}  // namespace tabula
