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


#include "SimdOverlapFunc.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_ss_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const CSimdMatrix    &coordinates,
                   const double          threshold) -> void
{
    if ((bra.get_angular_momentum() != 0) || (ket.get_angular_momentum() != 0))
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapFunc.compute_ss_overlap: Basis functions must be of zero angular momentum"));
    }

    if (nvalues > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("SimdOverlapFunc.compute_ss_overlap: Number of values exceeds number of atom pairs"));
    }

    if (nvalues == 0) return;

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_a = a_exps.size();

    const auto nprim_b = b_exps.size();

    const auto nprims = nprim_a * nprim_b;

    // NOTE: the pairs of primitives are screened with a tighter threshold than
    // the integrals, as their contributions accumulate into a single value.

    const auto dimensions = simdfunc::make_column_dimensions(
        bra, ket, nvalues, coordinates, screenfunc::two_center_overlap_primitive_bound, 1.0e-2 * threshold);

    // NOTE: the primitives are sorted by descending exponent, so the numbers of
    // surviving atom pairs do not decrease along either index and the last pair
    // of primitives is the one reaching furthest.

    const auto nmax = dimensions.back();

    if (nmax == 0)
    {
        std::fill(values, values + nvalues, 0.0);

        return;
    }

    // NOTE: each row holds the primitive integrals of one pair of primitives, and
    // the rows span the atom pairs surviving the screening of that pair and not
    // all atom pairs of the sparsity pattern. The atom pairs are ordered by
    // ascending interatomic distance, so the survivors are the leading ones and
    // the integrals of the remaining pairs are below threshold.

    auto buffer = CSimdMatrix(nprims, nmax);

    // NOTE: the squared distances of the atom pairs are carried by the
    // coordinates, so that they are formed once for the whole block instead of
    // once for every combination of basis functions.

    const auto *ab_2 = coordinates.data(6);

    constexpr auto fpi = mathconst::pi_value();

    // compute the primitive integrals of each pair of primitives

    for (size_t i = 0; i < nprim_a; i++)
    {
        const auto aexp = a_exps[i];

        const auto anorm = a_norms[i];

        for (size_t j = 0; j < nprim_b; j++)
        {
            const auto fexp = aexp + b_exps[j];

            const auto fmu = aexp * b_exps[j] / fexp;

            const auto fovl = fpi / fexp;

            const auto ffact = anorm * b_norms[j] * fovl * std::sqrt(fovl);

            const auto ncols = dimensions[i * nprim_b + j];

            if (ncols == 0) continue;

            auto *prim = buffer.data(i * nprim_b + j);

            // NOTE: the rows of the buffer and of the coordinates start at a cache
            // line boundary, so the loop is vectorized with aligned loads and
            // stores.

            // NOTE: the exponential is issued as a scalar call, as no vector math
            // library is linked, and it dominates the cost of the loop. Vectorizing
            // around it still pays, as the squared distance leaves only a load, two
            // multiplications and a store to be packed into the lanes.

#pragma omp simd aligned(prim, ab_2 : simd::cache_line_size())
            for (size_t k = 0; k < ncols; k++)
            {
                prim[k] = ffact * std::exp(-fmu * ab_2[k]);
            }
        }
    }

    // contract the primitive integrals into the values of the sparsity block

    // NOTE: the integrals of the pair of primitives reaching furthest are written
    // to the values and those of the remaining pairs are added to them, so that
    // the contraction needs neither a sweep to initialize the values nor a sweep
    // to copy them out of the buffer. Only the rows of the buffer are declared
    // aligned, as the values start at the offset of this combination of basis
    // functions in the values block.

    const auto *head = buffer.data(nprims - 1);

#pragma omp simd aligned(head : simd::cache_line_size())
    for (size_t k = 0; k < nmax; k++)
    {
        values[k] = head[k];
    }

    // NOTE: the atom pairs beyond the reach of every pair of primitives have no
    // contribution and are set to zero.

    std::fill(values + nmax, values + nvalues, 0.0);

    for (size_t irow = 0; irow + 1 < nprims; irow++)
    {
        const auto ncols = dimensions[irow];

        if (ncols == 0) continue;

        const auto *prim = buffer.data(irow);

#pragma omp simd aligned(prim : simd::cache_line_size())
        for (size_t k = 0; k < ncols; k++)
        {
            values[k] += prim[k];
        }
    }
}

auto
compute_overlap(double               *values,
                const size_t          nvalues,
                const CBasisFunction &bra,
                const CBasisFunction &ket,
                const CSimdMatrix    &coordinates,
                const double          threshold) -> void
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    if ((lbra == 0) && (lket == 0))
    {
        compute_ss_overlap(values, nvalues, bra, ket, coordinates, threshold);

        return;
    }

    errors::assertMsgCritical(
        false, std::string("SimdOverlapFunc.compute_overlap: Overlap integrals of the requested angular momenta are not implemented"));
}

}  // namespace simdovl
