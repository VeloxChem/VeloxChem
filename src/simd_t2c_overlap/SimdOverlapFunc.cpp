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
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_ss_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const CSimdMatrix    &coordinates) -> void
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

    // NOTE: the primitive integrals of a pair of primitives occupy one row, and
    // the last row is reserved for the contracted integrals, so that a single
    // allocation holds both and every row starts on a cache line boundary.

    // NOTE: the rows span the atom pairs surviving the screening of this
    // combination of basis functions and not all atom pairs of the sparsity
    // pattern. The atom pairs are ordered by ascending interatomic distance, so
    // the survivors are the leading ones, and the integrals of the remaining
    // pairs would be discarded when the values are written out.

    auto buffer = CSimdMatrix(nprims + 1, nvalues);

    const auto *a_x = coordinates.data(0);
    const auto *a_y = coordinates.data(1);
    const auto *a_z = coordinates.data(2);
    const auto *b_x = coordinates.data(3);
    const auto *b_y = coordinates.data(4);
    const auto *b_z = coordinates.data(5);

    constexpr auto fpi = mathconst::pi_value();

    // compute the primitive integrals of each pair of primitives

    // NOTE: the loops run over the atom pairs and not over the padded row, as
    // the padding of the coordinates is left undefined by their construction.

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

            auto *prim = buffer.data(i * nprim_b + j);

            // NOTE: the rows of the coordinates and of the buffer start at a cache
            // line boundary, so the loop is vectorized with aligned loads and stores.
            // The vectorization is conditional on the exponential being vectorized as
            // well, as it is otherwise slower than the scalar loop.

#pragma omp simd aligned(prim, a_x, a_y, a_z, b_x, b_y, b_z : simd::cache_line_size()) if (simd::has_vector_exp())
            for (size_t k = 0; k < nvalues; k++)
            {
                const auto rx = a_x[k] - b_x[k];

                const auto ry = a_y[k] - b_y[k];

                const auto rz = a_z[k] - b_z[k];

                prim[k] = ffact * std::exp(-fmu * (rx * rx + ry * ry + rz * rz));
            }
        }
    }

    // contract the primitive integrals into the reserved row

    // NOTE: the reserved row is seeded with the integrals of the first pair of
    // primitives instead of being zeroed, which saves a sweep over the row.

    auto *cvals = buffer.data(nprims);

    std::copy(buffer.data(0), buffer.data(0) + nvalues, cvals);

    for (size_t irow = 1; irow < nprims; irow++)
    {
        const auto *prim = buffer.data(irow);

#pragma omp simd aligned(cvals, prim : simd::cache_line_size())
        for (size_t k = 0; k < nvalues; k++)
        {
            cvals[k] += prim[k];
        }
    }

    std::copy(cvals, cvals + nvalues, values);
}

auto
compute_overlap(double               *values,
                const size_t          nvalues,
                const CBasisFunction &bra,
                const CBasisFunction &ket,
                const CSimdMatrix    &coordinates) -> void
{
    const auto lbra = bra.get_angular_momentum();

    const auto lket = ket.get_angular_momentum();

    if ((lbra == 0) && (lket == 0))
    {
        compute_ss_overlap(values, nvalues, bra, ket, coordinates);

        return;
    }

    errors::assertMsgCritical(
        false, std::string("SimdOverlapFunc.compute_overlap: Overlap integrals of the requested angular momenta are not implemented"));
}

}  // namespace simdovl
