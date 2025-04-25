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

#include "TwoCenterElectronRepulsionPrimRecHS.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_hs(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_hs,
                                const size_t idx_eri_0_fs,
                                const size_t idx_eri_1_fs,
                                const size_t idx_eri_1_gs,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpa,
                                const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : FS

    auto g_xxx_0_0 = pbuffer.data(idx_eri_0_fs);

    auto g_xyy_0_0 = pbuffer.data(idx_eri_0_fs + 3);

    auto g_xzz_0_0 = pbuffer.data(idx_eri_0_fs + 5);

    auto g_yyy_0_0 = pbuffer.data(idx_eri_0_fs + 6);

    auto g_yzz_0_0 = pbuffer.data(idx_eri_0_fs + 8);

    auto g_zzz_0_0 = pbuffer.data(idx_eri_0_fs + 9);

    // Set up components of auxiliary buffer : FS

    auto g_xxx_0_1 = pbuffer.data(idx_eri_1_fs);

    auto g_xyy_0_1 = pbuffer.data(idx_eri_1_fs + 3);

    auto g_xzz_0_1 = pbuffer.data(idx_eri_1_fs + 5);

    auto g_yyy_0_1 = pbuffer.data(idx_eri_1_fs + 6);

    auto g_yzz_0_1 = pbuffer.data(idx_eri_1_fs + 8);

    auto g_zzz_0_1 = pbuffer.data(idx_eri_1_fs + 9);

    // Set up components of auxiliary buffer : GS

    auto g_xxxx_0_1 = pbuffer.data(idx_eri_1_gs);

    auto g_xxxz_0_1 = pbuffer.data(idx_eri_1_gs + 2);

    auto g_xxyy_0_1 = pbuffer.data(idx_eri_1_gs + 3);

    auto g_xxzz_0_1 = pbuffer.data(idx_eri_1_gs + 5);

    auto g_xyyy_0_1 = pbuffer.data(idx_eri_1_gs + 6);

    auto g_xzzz_0_1 = pbuffer.data(idx_eri_1_gs + 9);

    auto g_yyyy_0_1 = pbuffer.data(idx_eri_1_gs + 10);

    auto g_yyyz_0_1 = pbuffer.data(idx_eri_1_gs + 11);

    auto g_yyzz_0_1 = pbuffer.data(idx_eri_1_gs + 12);

    auto g_yzzz_0_1 = pbuffer.data(idx_eri_1_gs + 13);

    auto g_zzzz_0_1 = pbuffer.data(idx_eri_1_gs + 14);

    // Set up components of targeted buffer : HS

    auto g_xxxxx_0_0 = pbuffer.data(idx_eri_0_hs);

    auto g_xxxxy_0_0 = pbuffer.data(idx_eri_0_hs + 1);

    auto g_xxxxz_0_0 = pbuffer.data(idx_eri_0_hs + 2);

    auto g_xxxyy_0_0 = pbuffer.data(idx_eri_0_hs + 3);

    auto g_xxxyz_0_0 = pbuffer.data(idx_eri_0_hs + 4);

    auto g_xxxzz_0_0 = pbuffer.data(idx_eri_0_hs + 5);

    auto g_xxyyy_0_0 = pbuffer.data(idx_eri_0_hs + 6);

    auto g_xxyyz_0_0 = pbuffer.data(idx_eri_0_hs + 7);

    auto g_xxyzz_0_0 = pbuffer.data(idx_eri_0_hs + 8);

    auto g_xxzzz_0_0 = pbuffer.data(idx_eri_0_hs + 9);

    auto g_xyyyy_0_0 = pbuffer.data(idx_eri_0_hs + 10);

    auto g_xyyyz_0_0 = pbuffer.data(idx_eri_0_hs + 11);

    auto g_xyyzz_0_0 = pbuffer.data(idx_eri_0_hs + 12);

    auto g_xyzzz_0_0 = pbuffer.data(idx_eri_0_hs + 13);

    auto g_xzzzz_0_0 = pbuffer.data(idx_eri_0_hs + 14);

    auto g_yyyyy_0_0 = pbuffer.data(idx_eri_0_hs + 15);

    auto g_yyyyz_0_0 = pbuffer.data(idx_eri_0_hs + 16);

    auto g_yyyzz_0_0 = pbuffer.data(idx_eri_0_hs + 17);

    auto g_yyzzz_0_0 = pbuffer.data(idx_eri_0_hs + 18);

    auto g_yzzzz_0_0 = pbuffer.data(idx_eri_0_hs + 19);

    auto g_zzzzz_0_0 = pbuffer.data(idx_eri_0_hs + 20);

    #pragma omp simd aligned(g_xxx_0_0, g_xxx_0_1, g_xxxx_0_1, g_xxxxx_0_0, g_xxxxy_0_0, g_xxxxz_0_0, g_xxxyy_0_0, g_xxxyz_0_0, g_xxxz_0_1, g_xxxzz_0_0, g_xxyy_0_1, g_xxyyy_0_0, g_xxyyz_0_0, g_xxyzz_0_0, g_xxzz_0_1, g_xxzzz_0_0, g_xyy_0_0, g_xyy_0_1, g_xyyy_0_1, g_xyyyy_0_0, g_xyyyz_0_0, g_xyyzz_0_0, g_xyzzz_0_0, g_xzz_0_0, g_xzz_0_1, g_xzzz_0_1, g_xzzzz_0_0, g_yyy_0_0, g_yyy_0_1, g_yyyy_0_1, g_yyyyy_0_0, g_yyyyz_0_0, g_yyyz_0_1, g_yyyzz_0_0, g_yyzz_0_1, g_yyzzz_0_0, g_yzz_0_0, g_yzz_0_1, g_yzzz_0_1, g_yzzzz_0_0, g_zzz_0_0, g_zzz_0_1, g_zzzz_0_1, g_zzzzz_0_0, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxxx_0_0[i] = 4.0 * g_xxx_0_0[i] * fbe_0 - 4.0 * g_xxx_0_1[i] * fz_be_0 + g_xxxx_0_1[i] * pa_x[i];

        g_xxxxy_0_0[i] = g_xxxx_0_1[i] * pa_y[i];

        g_xxxxz_0_0[i] = g_xxxx_0_1[i] * pa_z[i];

        g_xxxyy_0_0[i] = 2.0 * g_xyy_0_0[i] * fbe_0 - 2.0 * g_xyy_0_1[i] * fz_be_0 + g_xxyy_0_1[i] * pa_x[i];

        g_xxxyz_0_0[i] = g_xxxz_0_1[i] * pa_y[i];

        g_xxxzz_0_0[i] = 2.0 * g_xzz_0_0[i] * fbe_0 - 2.0 * g_xzz_0_1[i] * fz_be_0 + g_xxzz_0_1[i] * pa_x[i];

        g_xxyyy_0_0[i] = g_yyy_0_0[i] * fbe_0 - g_yyy_0_1[i] * fz_be_0 + g_xyyy_0_1[i] * pa_x[i];

        g_xxyyz_0_0[i] = g_xxyy_0_1[i] * pa_z[i];

        g_xxyzz_0_0[i] = g_xxzz_0_1[i] * pa_y[i];

        g_xxzzz_0_0[i] = g_zzz_0_0[i] * fbe_0 - g_zzz_0_1[i] * fz_be_0 + g_xzzz_0_1[i] * pa_x[i];

        g_xyyyy_0_0[i] = g_yyyy_0_1[i] * pa_x[i];

        g_xyyyz_0_0[i] = g_yyyz_0_1[i] * pa_x[i];

        g_xyyzz_0_0[i] = g_yyzz_0_1[i] * pa_x[i];

        g_xyzzz_0_0[i] = g_yzzz_0_1[i] * pa_x[i];

        g_xzzzz_0_0[i] = g_zzzz_0_1[i] * pa_x[i];

        g_yyyyy_0_0[i] = 4.0 * g_yyy_0_0[i] * fbe_0 - 4.0 * g_yyy_0_1[i] * fz_be_0 + g_yyyy_0_1[i] * pa_y[i];

        g_yyyyz_0_0[i] = g_yyyy_0_1[i] * pa_z[i];

        g_yyyzz_0_0[i] = 2.0 * g_yzz_0_0[i] * fbe_0 - 2.0 * g_yzz_0_1[i] * fz_be_0 + g_yyzz_0_1[i] * pa_y[i];

        g_yyzzz_0_0[i] = g_zzz_0_0[i] * fbe_0 - g_zzz_0_1[i] * fz_be_0 + g_yzzz_0_1[i] * pa_y[i];

        g_yzzzz_0_0[i] = g_zzzz_0_1[i] * pa_y[i];

        g_zzzzz_0_0[i] = 4.0 * g_zzz_0_0[i] * fbe_0 - 4.0 * g_zzz_0_1[i] * fz_be_0 + g_zzzz_0_1[i] * pa_z[i];
    }
}

} // t2ceri namespace

