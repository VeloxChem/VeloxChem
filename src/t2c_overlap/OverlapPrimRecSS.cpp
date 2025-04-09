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

#include "OverlapPrimRecSS.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace ovlrec {  // ovlrec namespace

auto
comp_prim_overlap_ss(CSimdArray<double>& pbuffer, const size_t idx_ovl_ss, CSimdArray<double>& factors, const double a_exp, const double a_norm)
    -> void
{
    const double fpi = mathconst::pi_value();

    // Set up exponents, normalization factors

    auto b_exps = factors.data(0);

    auto b_norms = factors.data(1);

    // Set up R(AB) distances

    auto ab_x = factors.data(5);

    auto ab_y = factors.data(6);

    auto ab_z = factors.data(7);

    // Set up components of targeted buffer : SS

    auto ts_0_0 = pbuffer.data(idx_ovl_ss);

    /// compute primitive integrals

    const auto nelems = pbuffer.number_of_active_elements();

#pragma omp simd aligned(ts_0_0, ab_x, ab_y, ab_z, b_exps, b_norms : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double fe_0 = 1.0 / (a_exp + b_exps[i]);

        double fz_0 = a_exp * b_exps[i] * fe_0 * (ab_x[i] * ab_x[i] + ab_y[i] * ab_y[i] + ab_z[i] * ab_z[i]);

        fe_0 *= fpi;

        ts_0_0[i] = b_norms[i] * a_norm * fe_0 * std::sqrt(fe_0) * std::exp(-fz_0);
    }
}

}  // namespace ovlrec
