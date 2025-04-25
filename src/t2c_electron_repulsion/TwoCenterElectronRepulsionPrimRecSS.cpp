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

#include "TwoCenterElectronRepulsionPrimRecSS.hpp"

#include <cmath>

#include "MathConst.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_ss(CSimdArray<double>& pbuffer,
                                const size_t idx_eri_ss,
                                const CSimdArray<double>& bf_data,
                                const size_t              idx_vals,
                                CSimdArray<double>&       factors,
                                const double              a_exp,
                                const double              a_norm) -> void
{
    const auto fpi = mathconst::pi_value();
    
    const auto fact = 2.0 * fpi * fpi * std::sqrt(fpi);

    // Set up exponents, normalization factors

    auto b_exps = factors.data(0);

    auto b_norms = factors.data(1);

    /// Boys function values

    auto bvals = bf_data.data(idx_vals);

    /// Set up components of auxiliary buffer : SS

    auto tg_0_0 = pbuffer.data(idx_eri_ss);

    /// compute primitive integrals

    const auto nelems = pbuffer.number_of_active_elements();

#pragma omp simd aligned(tg_0_0, bvals, b_exps, b_norms : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        tg_0_0[i] = fact * a_norm * b_norms[i] * bvals[i] / (a_exp * b_exps[i] * std::sqrt(a_exp + b_exps[i])) ;
    }
}

} // t2ceri namespace
