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



#include "SimdOverlapPrimRecSLP.hpp"

#include <cmath>

#include "MathConst.hpp"
#include "SimdAlign.hpp"

namespace simdovl {  // simdovl namespace

auto
compute_prim_slp_overlap(double *prim, const double *ab_2, const simdfunc::CPrimitivePair &pair, const bool on_ket) -> void
{
    const auto ncols = pair.ncols;

    const auto fexp = pair.aexp + pair.bexp;

    const auto fmu = pair.aexp * pair.bexp / fexp;

    constexpr auto fpi = mathconst::pi_value();

    const auto fovl = fpi / fexp;

    const auto fbase = pair.anorm * pair.bnorm * fovl * std::sqrt(fovl);

    const auto fr = on_ket ? (pair.aexp / fexp) : (-pair.bexp / fexp);

    // NOTE: the power of the displacement is written out, as the angular momentum
    // is small and a call to pow would be issued for every pair of primitives.

    const auto ffact = fbase * fr;

    // NOTE: the row of the buffer and the row of the coordinates start at a cache
    // line boundary, so the loop is vectorized with aligned loads and stores. A
    // pair of primitives contributes only to the atom pairs it reaches, so the loop
    // shortens as the primitives get tighter.

    // NOTE: the exponential is issued as a call to the vector math library of the
    // platform, so it does not break the vectorization of the loop.

#pragma omp simd aligned(prim, ab_2 : simd::cache_line_size())
    for (size_t k = 0; k < ncols; k++)
    {
        prim[k] += ffact * std::exp(-fmu * ab_2[k]);
    }
}

}  // namespace simdovl
