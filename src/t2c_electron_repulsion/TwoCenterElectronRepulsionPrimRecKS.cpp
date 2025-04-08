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

#include "TwoCenterElectronRepulsionPrimRecKS.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_ks(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_ks,
                                const size_t idx_eri_0_hs,
                                const size_t idx_eri_1_hs,
                                const size_t idx_eri_1_is,
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

    // Set up components of auxiliary buffer : HS

    auto g_xxxxx_0_0 = pbuffer.data(idx_eri_0_hs);

    auto g_xxxyy_0_0 = pbuffer.data(idx_eri_0_hs + 3);

    auto g_xxxzz_0_0 = pbuffer.data(idx_eri_0_hs + 5);

    auto g_xxyyy_0_0 = pbuffer.data(idx_eri_0_hs + 6);

    auto g_xxzzz_0_0 = pbuffer.data(idx_eri_0_hs + 9);

    auto g_xyyyy_0_0 = pbuffer.data(idx_eri_0_hs + 10);

    auto g_xyyzz_0_0 = pbuffer.data(idx_eri_0_hs + 12);

    auto g_xzzzz_0_0 = pbuffer.data(idx_eri_0_hs + 14);

    auto g_yyyyy_0_0 = pbuffer.data(idx_eri_0_hs + 15);

    auto g_yyyzz_0_0 = pbuffer.data(idx_eri_0_hs + 17);

    auto g_yyzzz_0_0 = pbuffer.data(idx_eri_0_hs + 18);

    auto g_yzzzz_0_0 = pbuffer.data(idx_eri_0_hs + 19);

    auto g_zzzzz_0_0 = pbuffer.data(idx_eri_0_hs + 20);

    // Set up components of auxiliary buffer : HS

    auto g_xxxxx_0_1 = pbuffer.data(idx_eri_1_hs);

    auto g_xxxyy_0_1 = pbuffer.data(idx_eri_1_hs + 3);

    auto g_xxxzz_0_1 = pbuffer.data(idx_eri_1_hs + 5);

    auto g_xxyyy_0_1 = pbuffer.data(idx_eri_1_hs + 6);

    auto g_xxzzz_0_1 = pbuffer.data(idx_eri_1_hs + 9);

    auto g_xyyyy_0_1 = pbuffer.data(idx_eri_1_hs + 10);

    auto g_xyyzz_0_1 = pbuffer.data(idx_eri_1_hs + 12);

    auto g_xzzzz_0_1 = pbuffer.data(idx_eri_1_hs + 14);

    auto g_yyyyy_0_1 = pbuffer.data(idx_eri_1_hs + 15);

    auto g_yyyzz_0_1 = pbuffer.data(idx_eri_1_hs + 17);

    auto g_yyzzz_0_1 = pbuffer.data(idx_eri_1_hs + 18);

    auto g_yzzzz_0_1 = pbuffer.data(idx_eri_1_hs + 19);

    auto g_zzzzz_0_1 = pbuffer.data(idx_eri_1_hs + 20);

    // Set up components of auxiliary buffer : IS

    auto g_xxxxxx_0_1 = pbuffer.data(idx_eri_1_is);

    auto g_xxxxxz_0_1 = pbuffer.data(idx_eri_1_is + 2);

    auto g_xxxxyy_0_1 = pbuffer.data(idx_eri_1_is + 3);

    auto g_xxxxzz_0_1 = pbuffer.data(idx_eri_1_is + 5);

    auto g_xxxyyy_0_1 = pbuffer.data(idx_eri_1_is + 6);

    auto g_xxxzzz_0_1 = pbuffer.data(idx_eri_1_is + 9);

    auto g_xxyyyy_0_1 = pbuffer.data(idx_eri_1_is + 10);

    auto g_xxyyzz_0_1 = pbuffer.data(idx_eri_1_is + 12);

    auto g_xxzzzz_0_1 = pbuffer.data(idx_eri_1_is + 14);

    auto g_xyyyyy_0_1 = pbuffer.data(idx_eri_1_is + 15);

    auto g_xyyyzz_0_1 = pbuffer.data(idx_eri_1_is + 17);

    auto g_xyyzzz_0_1 = pbuffer.data(idx_eri_1_is + 18);

    auto g_xzzzzz_0_1 = pbuffer.data(idx_eri_1_is + 20);

    auto g_yyyyyy_0_1 = pbuffer.data(idx_eri_1_is + 21);

    auto g_yyyyyz_0_1 = pbuffer.data(idx_eri_1_is + 22);

    auto g_yyyyzz_0_1 = pbuffer.data(idx_eri_1_is + 23);

    auto g_yyyzzz_0_1 = pbuffer.data(idx_eri_1_is + 24);

    auto g_yyzzzz_0_1 = pbuffer.data(idx_eri_1_is + 25);

    auto g_yzzzzz_0_1 = pbuffer.data(idx_eri_1_is + 26);

    auto g_zzzzzz_0_1 = pbuffer.data(idx_eri_1_is + 27);

    // Set up components of targeted buffer : KS

    auto g_xxxxxxx_0_0 = pbuffer.data(idx_eri_0_ks);

    auto g_xxxxxxy_0_0 = pbuffer.data(idx_eri_0_ks + 1);

    auto g_xxxxxxz_0_0 = pbuffer.data(idx_eri_0_ks + 2);

    auto g_xxxxxyy_0_0 = pbuffer.data(idx_eri_0_ks + 3);

    auto g_xxxxxyz_0_0 = pbuffer.data(idx_eri_0_ks + 4);

    auto g_xxxxxzz_0_0 = pbuffer.data(idx_eri_0_ks + 5);

    auto g_xxxxyyy_0_0 = pbuffer.data(idx_eri_0_ks + 6);

    auto g_xxxxyyz_0_0 = pbuffer.data(idx_eri_0_ks + 7);

    auto g_xxxxyzz_0_0 = pbuffer.data(idx_eri_0_ks + 8);

    auto g_xxxxzzz_0_0 = pbuffer.data(idx_eri_0_ks + 9);

    auto g_xxxyyyy_0_0 = pbuffer.data(idx_eri_0_ks + 10);

    auto g_xxxyyyz_0_0 = pbuffer.data(idx_eri_0_ks + 11);

    auto g_xxxyyzz_0_0 = pbuffer.data(idx_eri_0_ks + 12);

    auto g_xxxyzzz_0_0 = pbuffer.data(idx_eri_0_ks + 13);

    auto g_xxxzzzz_0_0 = pbuffer.data(idx_eri_0_ks + 14);

    auto g_xxyyyyy_0_0 = pbuffer.data(idx_eri_0_ks + 15);

    auto g_xxyyyyz_0_0 = pbuffer.data(idx_eri_0_ks + 16);

    auto g_xxyyyzz_0_0 = pbuffer.data(idx_eri_0_ks + 17);

    auto g_xxyyzzz_0_0 = pbuffer.data(idx_eri_0_ks + 18);

    auto g_xxyzzzz_0_0 = pbuffer.data(idx_eri_0_ks + 19);

    auto g_xxzzzzz_0_0 = pbuffer.data(idx_eri_0_ks + 20);

    auto g_xyyyyyy_0_0 = pbuffer.data(idx_eri_0_ks + 21);

    auto g_xyyyyyz_0_0 = pbuffer.data(idx_eri_0_ks + 22);

    auto g_xyyyyzz_0_0 = pbuffer.data(idx_eri_0_ks + 23);

    auto g_xyyyzzz_0_0 = pbuffer.data(idx_eri_0_ks + 24);

    auto g_xyyzzzz_0_0 = pbuffer.data(idx_eri_0_ks + 25);

    auto g_xyzzzzz_0_0 = pbuffer.data(idx_eri_0_ks + 26);

    auto g_xzzzzzz_0_0 = pbuffer.data(idx_eri_0_ks + 27);

    auto g_yyyyyyy_0_0 = pbuffer.data(idx_eri_0_ks + 28);

    auto g_yyyyyyz_0_0 = pbuffer.data(idx_eri_0_ks + 29);

    auto g_yyyyyzz_0_0 = pbuffer.data(idx_eri_0_ks + 30);

    auto g_yyyyzzz_0_0 = pbuffer.data(idx_eri_0_ks + 31);

    auto g_yyyzzzz_0_0 = pbuffer.data(idx_eri_0_ks + 32);

    auto g_yyzzzzz_0_0 = pbuffer.data(idx_eri_0_ks + 33);

    auto g_yzzzzzz_0_0 = pbuffer.data(idx_eri_0_ks + 34);

    auto g_zzzzzzz_0_0 = pbuffer.data(idx_eri_0_ks + 35);

    #pragma omp simd aligned(g_xxxxx_0_0, g_xxxxx_0_1, g_xxxxxx_0_1, g_xxxxxxx_0_0, g_xxxxxxy_0_0, g_xxxxxxz_0_0, g_xxxxxyy_0_0, g_xxxxxyz_0_0, g_xxxxxz_0_1, g_xxxxxzz_0_0, g_xxxxyy_0_1, g_xxxxyyy_0_0, g_xxxxyyz_0_0, g_xxxxyzz_0_0, g_xxxxzz_0_1, g_xxxxzzz_0_0, g_xxxyy_0_0, g_xxxyy_0_1, g_xxxyyy_0_1, g_xxxyyyy_0_0, g_xxxyyyz_0_0, g_xxxyyzz_0_0, g_xxxyzzz_0_0, g_xxxzz_0_0, g_xxxzz_0_1, g_xxxzzz_0_1, g_xxxzzzz_0_0, g_xxyyy_0_0, g_xxyyy_0_1, g_xxyyyy_0_1, g_xxyyyyy_0_0, g_xxyyyyz_0_0, g_xxyyyzz_0_0, g_xxyyzz_0_1, g_xxyyzzz_0_0, g_xxyzzzz_0_0, g_xxzzz_0_0, g_xxzzz_0_1, g_xxzzzz_0_1, g_xxzzzzz_0_0, g_xyyyy_0_0, g_xyyyy_0_1, g_xyyyyy_0_1, g_xyyyyyy_0_0, g_xyyyyyz_0_0, g_xyyyyzz_0_0, g_xyyyzz_0_1, g_xyyyzzz_0_0, g_xyyzz_0_0, g_xyyzz_0_1, g_xyyzzz_0_1, g_xyyzzzz_0_0, g_xyzzzzz_0_0, g_xzzzz_0_0, g_xzzzz_0_1, g_xzzzzz_0_1, g_xzzzzzz_0_0, g_yyyyy_0_0, g_yyyyy_0_1, g_yyyyyy_0_1, g_yyyyyyy_0_0, g_yyyyyyz_0_0, g_yyyyyz_0_1, g_yyyyyzz_0_0, g_yyyyzz_0_1, g_yyyyzzz_0_0, g_yyyzz_0_0, g_yyyzz_0_1, g_yyyzzz_0_1, g_yyyzzzz_0_0, g_yyzzz_0_0, g_yyzzz_0_1, g_yyzzzz_0_1, g_yyzzzzz_0_0, g_yzzzz_0_0, g_yzzzz_0_1, g_yzzzzz_0_1, g_yzzzzzz_0_0, g_zzzzz_0_0, g_zzzzz_0_1, g_zzzzzz_0_1, g_zzzzzzz_0_0, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);

        g_xxxxxxx_0_0[i] = 6.0 * g_xxxxx_0_0[i] * fbe_0 - 6.0 * g_xxxxx_0_1[i] * fz_be_0 + g_xxxxxx_0_1[i] * pa_x[i];

        g_xxxxxxy_0_0[i] = g_xxxxxx_0_1[i] * pa_y[i];

        g_xxxxxxz_0_0[i] = g_xxxxxx_0_1[i] * pa_z[i];

        g_xxxxxyy_0_0[i] = 4.0 * g_xxxyy_0_0[i] * fbe_0 - 4.0 * g_xxxyy_0_1[i] * fz_be_0 + g_xxxxyy_0_1[i] * pa_x[i];

        g_xxxxxyz_0_0[i] = g_xxxxxz_0_1[i] * pa_y[i];

        g_xxxxxzz_0_0[i] = 4.0 * g_xxxzz_0_0[i] * fbe_0 - 4.0 * g_xxxzz_0_1[i] * fz_be_0 + g_xxxxzz_0_1[i] * pa_x[i];

        g_xxxxyyy_0_0[i] = 3.0 * g_xxyyy_0_0[i] * fbe_0 - 3.0 * g_xxyyy_0_1[i] * fz_be_0 + g_xxxyyy_0_1[i] * pa_x[i];

        g_xxxxyyz_0_0[i] = g_xxxxyy_0_1[i] * pa_z[i];

        g_xxxxyzz_0_0[i] = g_xxxxzz_0_1[i] * pa_y[i];

        g_xxxxzzz_0_0[i] = 3.0 * g_xxzzz_0_0[i] * fbe_0 - 3.0 * g_xxzzz_0_1[i] * fz_be_0 + g_xxxzzz_0_1[i] * pa_x[i];

        g_xxxyyyy_0_0[i] = 2.0 * g_xyyyy_0_0[i] * fbe_0 - 2.0 * g_xyyyy_0_1[i] * fz_be_0 + g_xxyyyy_0_1[i] * pa_x[i];

        g_xxxyyyz_0_0[i] = g_xxxyyy_0_1[i] * pa_z[i];

        g_xxxyyzz_0_0[i] = 2.0 * g_xyyzz_0_0[i] * fbe_0 - 2.0 * g_xyyzz_0_1[i] * fz_be_0 + g_xxyyzz_0_1[i] * pa_x[i];

        g_xxxyzzz_0_0[i] = g_xxxzzz_0_1[i] * pa_y[i];

        g_xxxzzzz_0_0[i] = 2.0 * g_xzzzz_0_0[i] * fbe_0 - 2.0 * g_xzzzz_0_1[i] * fz_be_0 + g_xxzzzz_0_1[i] * pa_x[i];

        g_xxyyyyy_0_0[i] = g_yyyyy_0_0[i] * fbe_0 - g_yyyyy_0_1[i] * fz_be_0 + g_xyyyyy_0_1[i] * pa_x[i];

        g_xxyyyyz_0_0[i] = g_xxyyyy_0_1[i] * pa_z[i];

        g_xxyyyzz_0_0[i] = g_yyyzz_0_0[i] * fbe_0 - g_yyyzz_0_1[i] * fz_be_0 + g_xyyyzz_0_1[i] * pa_x[i];

        g_xxyyzzz_0_0[i] = g_yyzzz_0_0[i] * fbe_0 - g_yyzzz_0_1[i] * fz_be_0 + g_xyyzzz_0_1[i] * pa_x[i];

        g_xxyzzzz_0_0[i] = g_xxzzzz_0_1[i] * pa_y[i];

        g_xxzzzzz_0_0[i] = g_zzzzz_0_0[i] * fbe_0 - g_zzzzz_0_1[i] * fz_be_0 + g_xzzzzz_0_1[i] * pa_x[i];

        g_xyyyyyy_0_0[i] = g_yyyyyy_0_1[i] * pa_x[i];

        g_xyyyyyz_0_0[i] = g_yyyyyz_0_1[i] * pa_x[i];

        g_xyyyyzz_0_0[i] = g_yyyyzz_0_1[i] * pa_x[i];

        g_xyyyzzz_0_0[i] = g_yyyzzz_0_1[i] * pa_x[i];

        g_xyyzzzz_0_0[i] = g_yyzzzz_0_1[i] * pa_x[i];

        g_xyzzzzz_0_0[i] = g_yzzzzz_0_1[i] * pa_x[i];

        g_xzzzzzz_0_0[i] = g_zzzzzz_0_1[i] * pa_x[i];

        g_yyyyyyy_0_0[i] = 6.0 * g_yyyyy_0_0[i] * fbe_0 - 6.0 * g_yyyyy_0_1[i] * fz_be_0 + g_yyyyyy_0_1[i] * pa_y[i];

        g_yyyyyyz_0_0[i] = g_yyyyyy_0_1[i] * pa_z[i];

        g_yyyyyzz_0_0[i] = 4.0 * g_yyyzz_0_0[i] * fbe_0 - 4.0 * g_yyyzz_0_1[i] * fz_be_0 + g_yyyyzz_0_1[i] * pa_y[i];

        g_yyyyzzz_0_0[i] = 3.0 * g_yyzzz_0_0[i] * fbe_0 - 3.0 * g_yyzzz_0_1[i] * fz_be_0 + g_yyyzzz_0_1[i] * pa_y[i];

        g_yyyzzzz_0_0[i] = 2.0 * g_yzzzz_0_0[i] * fbe_0 - 2.0 * g_yzzzz_0_1[i] * fz_be_0 + g_yyzzzz_0_1[i] * pa_y[i];

        g_yyzzzzz_0_0[i] = g_zzzzz_0_0[i] * fbe_0 - g_zzzzz_0_1[i] * fz_be_0 + g_yzzzzz_0_1[i] * pa_y[i];

        g_yzzzzzz_0_0[i] = g_zzzzzz_0_1[i] * pa_y[i];

        g_zzzzzzz_0_0[i] = 6.0 * g_zzzzz_0_0[i] * fbe_0 - 6.0 * g_zzzzz_0_1[i] * fz_be_0 + g_zzzzzz_0_1[i] * pa_z[i];
    }
}

} // t2ceri namespace

