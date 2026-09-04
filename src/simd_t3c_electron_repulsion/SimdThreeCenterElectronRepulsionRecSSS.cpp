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



#include "SimdThreeCenterElectronRepulsionRecSSS.hpp"

#include <algorithm>
#include <cmath>
#include <ranges>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdAlign.hpp"
#include "SimdBoysFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdVariableMatrix.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_sss_electron_repulsion(double               *values,
                               const size_t          npairs,
                               const size_t          natoms,
                               const size_t          iatom,
                               const CBasisFunction &a_function,
                               const CBasisFunction &b_function,
                               const CBasisFunction &c_function,
                               const CSimdMatrix    &ab_coordinates,
                               const CSimdMatrix    &bc_coordinates,
                               const double          threshold) -> void
{
    if ((a_function.get_angular_momentum() != 0) || (b_function.get_angular_momentum() != 0) ||
        (c_function.get_angular_momentum() != 0))
    {
        errors::assertMsgCritical(
            false,
            std::string("SimdThreeCenterElectronRepulsionRecSSS.compute_sss_electron_repulsion: Basis functions must be of zero angular momentum"));
    }

    if (npairs > ab_coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdThreeCenterElectronRepulsionRecSSS.compute_sss_electron_repulsion: Number of atom pairs exceeds coordinates"));
    }

    if (iatom >= natoms)
    {
        errors::assertMsgCritical(
            false, std::string("SimdThreeCenterElectronRepulsionRecSSS.compute_sss_electron_repulsion: Index of atom on c side is out of range"));
    }

    // NOTE: the values of one atom on c side are contiguous over the atom pairs,
    // and the atoms are npairs apart, as the sparsity pattern lays them out.

    auto *slice = values + iatom * npairs;

    if (npairs == 0) return;

    const auto &a_exps = a_function.exponents();

    const auto &b_exps = b_function.exponents();

    const auto &c_exps = c_function.exponents();

    const auto &a_norms = a_function.normalization_factors();

    const auto &b_norms = b_function.normalization_factors();

    const auto &c_norms = c_function.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprim_c = c_exps.size();

    const auto nprims = nprim_a * nprim_b * nprim_c;

    // NOTE: the triples of primitives are screened with the threshold of the
    // integrals divided by their number, as their contributions accumulate into
    // a single value and the error of the sum is bounded by the number of terms.

    // NOTE: only the distance between the atoms on a and b sides enters the
    // bound, the dependence on the position of the atom on c side being
    // neglected, so the surviving atom pairs of a triple are a leading subrange
    // of the atom pairs of the block and no gathering is needed.

    const auto dimensions = simdfunc::make_column_dimensions(a_function,
                                                             b_function,
                                                             c_function,
                                                             npairs,
                                                             ab_coordinates,
                                                             screenfunc::three_center_electron_repulsion_primitive_bound,
                                                             threshold / static_cast<double>(nprims));

    // NOTE: the buffers span the atom pairs reached by the triple of primitives
    // reaching furthest, which is searched for rather than assumed. The
    // primitives are sorted by descending exponent, but the bound of a triple
    // carries its prefactor as well as its decay, so a tighter triple with a
    // larger prefactor reaches further than a more diffuse one with a smaller
    // prefactor, and the last triple is not always the furthest reaching.

    const auto nmax = *std::ranges::max_element(dimensions);

    if (nmax == 0)
    {
        std::fill(slice, slice + npairs, 0.0);

        return;
    }

    const auto *a_x = ab_coordinates.data(0);
    const auto *a_y = ab_coordinates.data(1);
    const auto *a_z = ab_coordinates.data(2);

    const auto *ab_2 = ab_coordinates.data(6);

    const auto *b_x = bc_coordinates.data(0);
    const auto *b_y = bc_coordinates.data(1);
    const auto *b_z = bc_coordinates.data(2);

    const auto *c_x = bc_coordinates.data(3);
    const auto *c_y = bc_coordinates.data(4);
    const auto *c_z = bc_coordinates.data(5);

    // NOTE: the integrals of all triples of primitives are accumulated in a
    // single row, which starts at a cache line boundary and spans only the atom
    // pairs reached by the furthest reaching triple.

    auto buffer = CSimdMatrix(1, nmax);

    buffer.zero();

    auto *prim = buffer.data(0);

    // NOTE: the exponential of the pair on a and b sides and the squared distance
    // from their product center to the atom on c side depend on the exponents of
    // that pair alone, the ket exponent entering neither, so both are formed once
    // for a pair of primitives and read by every primitive on c side.

    auto factors = CSimdMatrix(2, nmax);

    auto *e_ab = factors.data(0);

    auto *pc_2 = factors.data(1);

    constexpr auto fpi = mathconst::pi_value();

    // NOTE: the three-center repulsion of three S type primitives is two pi to
    // the five halves over the exponent of the pair times the exponent on c side
    // times the square root of their sum, times the exponential of the pair and
    // the Boys function of order zero.

    const auto fcoul = 2.0 * fpi * fpi * std::sqrt(fpi);

    for (size_t i = 0; i < nprim_a; i++)
    {
        const auto aexp = a_exps[i];

        const auto anorm = a_norms[i];

        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto bexp = b_exps[j];

            // NOTE: the widest of the triples of this pair of primitives is
            // searched for rather than assumed to be the last. The primitives on
            // c side are sorted by descending exponent, but the bound of a triple
            // carries its prefactor as well as its decay, so a tighter primitive
            // with a larger prefactor reaches further than a more diffuse one
            // with a smaller prefactor. The factors are formed over the atom
            // pairs of the widest triple and the narrower ones read a leading
            // part of them.

            const auto first = dimensions.begin() + static_cast<long>((i * nprim_b + j) * nprim_c);

            const auto npair_max = *std::ranges::max_element(first, first + static_cast<long>(nprim_c));

            if (npair_max == 0) continue;

            const auto pexp = aexp + bexp;

            const auto fmu = aexp * bexp / pexp;

            const auto frp = 1.0 / pexp;

            // NOTE: the product center of the pair is the exponent weighted mean
            // of the atoms on a and b sides, so the vector from it to the atom on
            // c side is the same mean of the vectors from those atoms to it.

#pragma omp simd aligned(e_ab, pc_2, ab_2, a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z : simd::cache_line_size())
            for (size_t k = 0; k < npair_max; k++)
            {
                e_ab[k] = std::exp(-fmu * ab_2[k]);

                const auto p_x = frp * (aexp * (a_x[k] - c_x[k]) + bexp * (b_x[k] - c_x[k]));

                const auto p_y = frp * (aexp * (a_y[k] - c_y[k]) + bexp * (b_y[k] - c_y[k]));

                const auto p_z = frp * (aexp * (a_z[k] - c_z[k]) + bexp * (b_z[k] - c_z[k]));

                pc_2[k] = p_x * p_x + p_y * p_y + p_z * p_z;
            }

            // NOTE: the Boys function of every primitive on c side of this pair
            // is computed by one call, which fills the order zero of every row.
            // The rows are the atom pairs the triples reach, so a tighter
            // primitive on c side carries a shorter row.

            auto boys = CSimdVariableMatrix(std::vector<size_t>(first, first + static_cast<long>(nprim_c)), 2);

            for (size_t k = 0; k < nprim_c; k++)
            {
                const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                if (ncols == 0) continue;

                // NOTE: the pair is weighed against the primitive on c side by
                // the reduced exponent of the two, which is what scales the
                // squared distance the Boys function is evaluated at.

                const auto frho = pexp * c_exps[k] / (pexp + c_exps[k]);

                auto *bargs = boys.data(0, k);

#pragma omp simd aligned(bargs, pc_2 : simd::cache_line_size())
                for (size_t l = 0; l < ncols; l++)
                {
                    bargs[l] = frho * pc_2[l];
                }
            }

            simdfunc::compute_boys_function(boys);

            for (size_t k = 0; k < nprim_c; k++)
            {
                const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                if (ncols == 0) continue;

                const auto cexp = c_exps[k];

                const auto qexp = pexp + cexp;

                const auto ffact = fcoul * anorm * b_norms[j] * c_norms[k] / (pexp * cexp * std::sqrt(qexp));

                const auto *bvals = boys.data(1, k);

                // NOTE: the rows of the buffers and of the Boys function start at
                // a cache line boundary, so the loop is vectorized with aligned
                // loads and stores. A triple of primitives contributes only to the
                // atom pairs it reaches, so the loop shortens as they get tighter.

#pragma omp simd aligned(prim, e_ab, bvals : simd::cache_line_size())
                for (size_t l = 0; l < ncols; l++)
                {
                    prim[l] += ffact * e_ab[l] * bvals[l];
                }
            }
        }
    }

    // NOTE: the atom pairs beyond the reach of every triple of primitives have no
    // contribution and are set to zero.

    std::copy(prim, prim + nmax, slice);

    std::fill(slice + nmax, slice + npairs, 0.0);
}

}  // namespace simdt3ceri
