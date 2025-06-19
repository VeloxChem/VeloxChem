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

#include "LocalCorePotentialPrimRecSS.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace ecprec {  // ovlrec namespace

auto comp_prim_local_core_potential_ss(      CSimdArray<double>& pbuffer,
                                       const size_t              idx_lpot_ss,
                                             CSimdArray<double>& factors,
                                       const size_t              idx_r,
                                       const TPoint<double>&     r_a,
                                       const double              a_exp,
                                       const double              a_norm,
                                       const double              c_exp,
                                       const double              c_fact) -> void
{
    const double fpi = mathconst::pi_value();

    // Set up exponents, normalization factors

    auto b_exps = factors.data(0);

    auto b_norms = factors.data(1);
    
    // set up Cartesian B coordinates

    auto b_x = factors.data(2);

    auto b_y = factors.data(3);

    auto b_z = factors.data(4);

    // set up Cartesian R coordinates

    auto r_x = factors.data(idx_r);

    auto r_y = factors.data(idx_r + 1);

    auto r_z = factors.data(idx_r + 2);
    
    // set up Cartesian A coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    const auto fza = -a_exp * (a_x * a_x + a_y * a_y + a_z * a_z);
    
    // Set up components of targeted buffer : SS

    auto ts_0_0 = pbuffer.data(idx_lpot_ss);

    /// compute primitive integrals

    const auto nelems = pbuffer.number_of_active_elements();

#pragma omp simd aligned(ts_0_0, b_x, b_y, b_z, b_exps, b_norms, r_x, r_y, r_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double fg_0 = a_exp + b_exps[i] + c_exp;
        
        double fe_0 = 1.0 / fg_0;

        double fz_0 = fza - b_exps[i] * (b_x[i] * b_x[i] + b_y[i] * b_y[i] + b_z[i] * b_z[i])
                   
                    + fg_0 * (r_x[i] * r_x[i] + r_y[i] * r_y[i] + r_z[i] * r_z[i]);

        fe_0 *= fpi;

        ts_0_0[i] = c_fact * b_norms[i] * a_norm * fe_0 * std::sqrt(fe_0) * std::exp(-fz_0);
    }
    
}

}  // namespace ecprec
