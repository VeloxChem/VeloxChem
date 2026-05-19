//
//  Tabula — custom-recursion molecular-integral machinery.
//  Two-center overlap recursion.
//

#include "TabulaOverlapRecursion.hpp"

#include <cstddef>

namespace tabula {  // tabula namespace

auto
compute_overlap_seed(const CGtoPairBlock& pair_block) -> std::vector<double>
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

    std::vector<double> seed(static_cast<std::size_t>(order + 1) * stride, 0.0);

    if (pdim == 0) return seed;

    const auto bra_exps = pair_block.bra_exponents();
    const auto ket_exps = pair_block.ket_exponents();
    const auto overlaps = pair_block.overlap_factors();
    const auto norms    = pair_block.normalization_factors();

    // row 0 — [0]^0 = c_i·c_j · overlap_factor · (1/2α)^l_a · (−1/2γ)^l_c,
    // the primitive contraction weight (normalization_factors) premultiplied
    // in so the contraction step is a plain sum. The powers are applied as
    // l_a + l_c separate SIMD passes; a runtime-trip inner power loop would
    // block vectorization of the primitive-pair sweep.
    {
        auto* seed0 = seed.data();

#pragma omp simd
        for (std::size_t k = 0; k < pdim; k++)
        {
            seed0[k] = norms[k] * overlaps[k];
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

    // rows 1 … L — the seed ladder [0]^m = (−2ρ)·[0]^(m−1), ρ = αγ/(α+γ)
    for (int m = 1; m <= order; m++)
    {
        const auto* prev = seed.data() + static_cast<std::size_t>(m - 1) * stride;
        auto*       curr = seed.data() + static_cast<std::size_t>(m) * stride;

#pragma omp simd
        for (std::size_t k = 0; k < pdim; k++)
        {
            const auto a   = bra_exps[k];
            const auto g   = ket_exps[k];
            const auto rho = a * g / (a + g);

            curr[k] = -2.0 * rho * prev[k];
        }
    }

    return seed;
}

}  // namespace tabula
