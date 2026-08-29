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
#include <ranges>
#include <string>
#include <vector>

#include "ErrorHandler.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdVariableMatrix.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_ss_overlap(double               *values,
                   const size_t          nvalues,
                   const CBasisFunction &bra,
                   const CBasisFunction &ket,
                   const size_t          npairs,
                   const CSimdMatrix    &coordinates,
                   const double          threshold) -> void
{
    errors::assertMsgCritical((bra.get_angular_momentum() == 0) && (ket.get_angular_momentum() == 0),
                              std::string("SimdOverlapFunc.compute_ss_overlap: Basis functions must be of zero angular momentum"));

    errors::assertMsgCritical(nvalues <= npairs,
                              std::string("SimdOverlapFunc.compute_ss_overlap: Number of values exceeds number of atom pairs"));

    if (nvalues == 0) return;

    // NOTE: the pairs of primitives are screened with a tighter threshold than
    // the integrals, as their contributions accumulate into a single value.

    const auto dimensions =
        simdfunc::make_column_dimensions(bra, ket, npairs, coordinates, screenfunc::two_center_overlap_primitive_bound,
                                         1.0e-2 * threshold);

    auto buffer = CSimdVariableMatrix(dimensions, 1);

    const auto *a_x = coordinates.data(0);
    const auto *a_y = coordinates.data(1);
    const auto *a_z = coordinates.data(2);
    const auto *b_x = coordinates.data(3);
    const auto *b_y = coordinates.data(4);
    const auto *b_z = coordinates.data(5);

    const auto &a_exps = bra.exponents();

    const auto &b_exps = ket.exponents();

    const auto &a_norms = bra.normalization_factors();

    const auto &b_norms = ket.normalization_factors();

    const auto nprim_b = b_exps.size();

    constexpr auto fpi = mathconst::pi_value();

    // compute the primitive integrals of each pair of primitives

    std::ranges::for_each(std::views::iota(size_t{0}, dimensions.size()), [&](const auto irow) {
        const auto ncols = dimensions[irow];

        if (ncols == 0) return;

        const auto aexp = a_exps[irow / nprim_b];

        const auto bexp = b_exps[irow % nprim_b];

        const auto fnorm = a_norms[irow / nprim_b] * b_norms[irow % nprim_b];

        const auto fexp = aexp + bexp;

        const auto fmu = aexp * bexp / fexp;

        const auto fovl = fpi / fexp;

        const auto ffact = fnorm * fovl * std::sqrt(fovl);

        auto *prim = buffer.data(0, irow);

        for (size_t i = 0; i < ncols; i++)
        {
            const auto rx = a_x[i] - b_x[i];

            const auto ry = a_y[i] - b_y[i];

            const auto rz = a_z[i] - b_z[i];

            prim[i] = ffact * std::exp(-fmu * (rx * rx + ry * ry + rz * rz));
        }
    });

    // contract the primitive integrals into the values of the atom pairs

    auto contracted = CSimdMatrix(1, npairs);

    contracted.zero();

    auto *cvals = contracted.data(0);

    std::ranges::for_each(std::views::iota(size_t{0}, dimensions.size()), [&](const auto irow) {
        const auto ncols = std::min(dimensions[irow], npairs);

        const auto *prim = buffer.data(0, irow);

        for (size_t i = 0; i < ncols; i++)
        {
            cvals[i] += prim[i];
        }
    });

    std::copy(cvals, cvals + nvalues, values);
}

}  // namespace simdovl
