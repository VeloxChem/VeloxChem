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

#include "ThreeCenterElectronRepulsionPrimRecKSS.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_kss(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_kss,
                                 size_t idx_eri_0_hss,
                                 size_t idx_eri_1_hss,
                                 size_t idx_eri_1_iss,
                                 CSimdArray<double>& factors,
                                 const size_t idx_wa,
                                 const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WA) distances

    auto wa_x = factors.data(idx_wa);

    auto wa_y = factors.data(idx_wa + 1);

    auto wa_z = factors.data(idx_wa + 2);

    /// Set up components of auxilary buffer : HSS

    auto g_xxxxx_0_0_0 = pbuffer.data(idx_eri_0_hss);

    auto g_xxxyy_0_0_0 = pbuffer.data(idx_eri_0_hss + 3);

    auto g_xxxzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 5);

    auto g_xxyyy_0_0_0 = pbuffer.data(idx_eri_0_hss + 6);

    auto g_xxzzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 9);

    auto g_xyyyy_0_0_0 = pbuffer.data(idx_eri_0_hss + 10);

    auto g_xyyzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 12);

    auto g_xzzzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 14);

    auto g_yyyyy_0_0_0 = pbuffer.data(idx_eri_0_hss + 15);

    auto g_yyyzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 17);

    auto g_yyzzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 18);

    auto g_yzzzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 19);

    auto g_zzzzz_0_0_0 = pbuffer.data(idx_eri_0_hss + 20);

    /// Set up components of auxilary buffer : HSS

    auto g_xxxxx_0_0_1 = pbuffer.data(idx_eri_1_hss);

    auto g_xxxyy_0_0_1 = pbuffer.data(idx_eri_1_hss + 3);

    auto g_xxxzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 5);

    auto g_xxyyy_0_0_1 = pbuffer.data(idx_eri_1_hss + 6);

    auto g_xxzzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 9);

    auto g_xyyyy_0_0_1 = pbuffer.data(idx_eri_1_hss + 10);

    auto g_xyyzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 12);

    auto g_xzzzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 14);

    auto g_yyyyy_0_0_1 = pbuffer.data(idx_eri_1_hss + 15);

    auto g_yyyzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 17);

    auto g_yyzzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 18);

    auto g_yzzzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 19);

    auto g_zzzzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 20);

    /// Set up components of auxilary buffer : ISS

    auto g_xxxxxx_0_0_1 = pbuffer.data(idx_eri_1_iss);

    auto g_xxxxxz_0_0_1 = pbuffer.data(idx_eri_1_iss + 2);

    auto g_xxxxyy_0_0_1 = pbuffer.data(idx_eri_1_iss + 3);

    auto g_xxxxzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 5);

    auto g_xxxyyy_0_0_1 = pbuffer.data(idx_eri_1_iss + 6);

    auto g_xxxzzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 9);

    auto g_xxyyyy_0_0_1 = pbuffer.data(idx_eri_1_iss + 10);

    auto g_xxyyzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 12);

    auto g_xxzzzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 14);

    auto g_xyyyyy_0_0_1 = pbuffer.data(idx_eri_1_iss + 15);

    auto g_xyyyzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 17);

    auto g_xyyzzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 18);

    auto g_xzzzzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 20);

    auto g_yyyyyy_0_0_1 = pbuffer.data(idx_eri_1_iss + 21);

    auto g_yyyyyz_0_0_1 = pbuffer.data(idx_eri_1_iss + 22);

    auto g_yyyyzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 23);

    auto g_yyyzzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 24);

    auto g_yyzzzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 25);

    auto g_yzzzzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 26);

    auto g_zzzzzz_0_0_1 = pbuffer.data(idx_eri_1_iss + 27);

    /// Set up components of targeted buffer : KSS

    auto g_xxxxxxx_0_0_0 = pbuffer.data(idx_eri_0_kss);

    auto g_xxxxxxy_0_0_0 = pbuffer.data(idx_eri_0_kss + 1);

    auto g_xxxxxxz_0_0_0 = pbuffer.data(idx_eri_0_kss + 2);

    auto g_xxxxxyy_0_0_0 = pbuffer.data(idx_eri_0_kss + 3);

    auto g_xxxxxyz_0_0_0 = pbuffer.data(idx_eri_0_kss + 4);

    auto g_xxxxxzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 5);

    auto g_xxxxyyy_0_0_0 = pbuffer.data(idx_eri_0_kss + 6);

    auto g_xxxxyyz_0_0_0 = pbuffer.data(idx_eri_0_kss + 7);

    auto g_xxxxyzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 8);

    auto g_xxxxzzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 9);

    auto g_xxxyyyy_0_0_0 = pbuffer.data(idx_eri_0_kss + 10);

    auto g_xxxyyyz_0_0_0 = pbuffer.data(idx_eri_0_kss + 11);

    auto g_xxxyyzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 12);

    auto g_xxxyzzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 13);

    auto g_xxxzzzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 14);

    auto g_xxyyyyy_0_0_0 = pbuffer.data(idx_eri_0_kss + 15);

    auto g_xxyyyyz_0_0_0 = pbuffer.data(idx_eri_0_kss + 16);

    auto g_xxyyyzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 17);

    auto g_xxyyzzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 18);

    auto g_xxyzzzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 19);

    auto g_xxzzzzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 20);

    auto g_xyyyyyy_0_0_0 = pbuffer.data(idx_eri_0_kss + 21);

    auto g_xyyyyyz_0_0_0 = pbuffer.data(idx_eri_0_kss + 22);

    auto g_xyyyyzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 23);

    auto g_xyyyzzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 24);

    auto g_xyyzzzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 25);

    auto g_xyzzzzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 26);

    auto g_xzzzzzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 27);

    auto g_yyyyyyy_0_0_0 = pbuffer.data(idx_eri_0_kss + 28);

    auto g_yyyyyyz_0_0_0 = pbuffer.data(idx_eri_0_kss + 29);

    auto g_yyyyyzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 30);

    auto g_yyyyzzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 31);

    auto g_yyyzzzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 32);

    auto g_yyzzzzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 33);

    auto g_yzzzzzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 34);

    auto g_zzzzzzz_0_0_0 = pbuffer.data(idx_eri_0_kss + 35);

    #pragma omp simd aligned(g_xxxxx_0_0_0, g_xxxxx_0_0_1, g_xxxxxx_0_0_1, g_xxxxxxx_0_0_0, g_xxxxxxy_0_0_0, g_xxxxxxz_0_0_0, g_xxxxxyy_0_0_0, g_xxxxxyz_0_0_0, g_xxxxxz_0_0_1, g_xxxxxzz_0_0_0, g_xxxxyy_0_0_1, g_xxxxyyy_0_0_0, g_xxxxyyz_0_0_0, g_xxxxyzz_0_0_0, g_xxxxzz_0_0_1, g_xxxxzzz_0_0_0, g_xxxyy_0_0_0, g_xxxyy_0_0_1, g_xxxyyy_0_0_1, g_xxxyyyy_0_0_0, g_xxxyyyz_0_0_0, g_xxxyyzz_0_0_0, g_xxxyzzz_0_0_0, g_xxxzz_0_0_0, g_xxxzz_0_0_1, g_xxxzzz_0_0_1, g_xxxzzzz_0_0_0, g_xxyyy_0_0_0, g_xxyyy_0_0_1, g_xxyyyy_0_0_1, g_xxyyyyy_0_0_0, g_xxyyyyz_0_0_0, g_xxyyyzz_0_0_0, g_xxyyzz_0_0_1, g_xxyyzzz_0_0_0, g_xxyzzzz_0_0_0, g_xxzzz_0_0_0, g_xxzzz_0_0_1, g_xxzzzz_0_0_1, g_xxzzzzz_0_0_0, g_xyyyy_0_0_0, g_xyyyy_0_0_1, g_xyyyyy_0_0_1, g_xyyyyyy_0_0_0, g_xyyyyyz_0_0_0, g_xyyyyzz_0_0_0, g_xyyyzz_0_0_1, g_xyyyzzz_0_0_0, g_xyyzz_0_0_0, g_xyyzz_0_0_1, g_xyyzzz_0_0_1, g_xyyzzzz_0_0_0, g_xyzzzzz_0_0_0, g_xzzzz_0_0_0, g_xzzzz_0_0_1, g_xzzzzz_0_0_1, g_xzzzzzz_0_0_0, g_yyyyy_0_0_0, g_yyyyy_0_0_1, g_yyyyyy_0_0_1, g_yyyyyyy_0_0_0, g_yyyyyyz_0_0_0, g_yyyyyz_0_0_1, g_yyyyyzz_0_0_0, g_yyyyzz_0_0_1, g_yyyyzzz_0_0_0, g_yyyzz_0_0_0, g_yyyzz_0_0_1, g_yyyzzz_0_0_1, g_yyyzzzz_0_0_0, g_yyzzz_0_0_0, g_yyzzz_0_0_1, g_yyzzzz_0_0_1, g_yyzzzzz_0_0_0, g_yzzzz_0_0_0, g_yzzzz_0_0_1, g_yzzzzz_0_0_1, g_yzzzzzz_0_0_0, g_zzzzz_0_0_0, g_zzzzz_0_0_1, g_zzzzzz_0_0_1, g_zzzzzzz_0_0_0, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxxx_0_0_0[i] = 6.0 * g_xxxxx_0_0_0[i] * fbe_0 - 6.0 * g_xxxxx_0_0_1[i] * fz_be_0 + g_xxxxxx_0_0_1[i] * wa_x[i];

        g_xxxxxxy_0_0_0[i] = g_xxxxxx_0_0_1[i] * wa_y[i];

        g_xxxxxxz_0_0_0[i] = g_xxxxxx_0_0_1[i] * wa_z[i];

        g_xxxxxyy_0_0_0[i] = 4.0 * g_xxxyy_0_0_0[i] * fbe_0 - 4.0 * g_xxxyy_0_0_1[i] * fz_be_0 + g_xxxxyy_0_0_1[i] * wa_x[i];

        g_xxxxxyz_0_0_0[i] = g_xxxxxz_0_0_1[i] * wa_y[i];

        g_xxxxxzz_0_0_0[i] = 4.0 * g_xxxzz_0_0_0[i] * fbe_0 - 4.0 * g_xxxzz_0_0_1[i] * fz_be_0 + g_xxxxzz_0_0_1[i] * wa_x[i];

        g_xxxxyyy_0_0_0[i] = 3.0 * g_xxyyy_0_0_0[i] * fbe_0 - 3.0 * g_xxyyy_0_0_1[i] * fz_be_0 + g_xxxyyy_0_0_1[i] * wa_x[i];

        g_xxxxyyz_0_0_0[i] = g_xxxxyy_0_0_1[i] * wa_z[i];

        g_xxxxyzz_0_0_0[i] = g_xxxxzz_0_0_1[i] * wa_y[i];

        g_xxxxzzz_0_0_0[i] = 3.0 * g_xxzzz_0_0_0[i] * fbe_0 - 3.0 * g_xxzzz_0_0_1[i] * fz_be_0 + g_xxxzzz_0_0_1[i] * wa_x[i];

        g_xxxyyyy_0_0_0[i] = 2.0 * g_xyyyy_0_0_0[i] * fbe_0 - 2.0 * g_xyyyy_0_0_1[i] * fz_be_0 + g_xxyyyy_0_0_1[i] * wa_x[i];

        g_xxxyyyz_0_0_0[i] = g_xxxyyy_0_0_1[i] * wa_z[i];

        g_xxxyyzz_0_0_0[i] = 2.0 * g_xyyzz_0_0_0[i] * fbe_0 - 2.0 * g_xyyzz_0_0_1[i] * fz_be_0 + g_xxyyzz_0_0_1[i] * wa_x[i];

        g_xxxyzzz_0_0_0[i] = g_xxxzzz_0_0_1[i] * wa_y[i];

        g_xxxzzzz_0_0_0[i] = 2.0 * g_xzzzz_0_0_0[i] * fbe_0 - 2.0 * g_xzzzz_0_0_1[i] * fz_be_0 + g_xxzzzz_0_0_1[i] * wa_x[i];

        g_xxyyyyy_0_0_0[i] = g_yyyyy_0_0_0[i] * fbe_0 - g_yyyyy_0_0_1[i] * fz_be_0 + g_xyyyyy_0_0_1[i] * wa_x[i];

        g_xxyyyyz_0_0_0[i] = g_xxyyyy_0_0_1[i] * wa_z[i];

        g_xxyyyzz_0_0_0[i] = g_yyyzz_0_0_0[i] * fbe_0 - g_yyyzz_0_0_1[i] * fz_be_0 + g_xyyyzz_0_0_1[i] * wa_x[i];

        g_xxyyzzz_0_0_0[i] = g_yyzzz_0_0_0[i] * fbe_0 - g_yyzzz_0_0_1[i] * fz_be_0 + g_xyyzzz_0_0_1[i] * wa_x[i];

        g_xxyzzzz_0_0_0[i] = g_xxzzzz_0_0_1[i] * wa_y[i];

        g_xxzzzzz_0_0_0[i] = g_zzzzz_0_0_0[i] * fbe_0 - g_zzzzz_0_0_1[i] * fz_be_0 + g_xzzzzz_0_0_1[i] * wa_x[i];

        g_xyyyyyy_0_0_0[i] = g_yyyyyy_0_0_1[i] * wa_x[i];

        g_xyyyyyz_0_0_0[i] = g_yyyyyz_0_0_1[i] * wa_x[i];

        g_xyyyyzz_0_0_0[i] = g_yyyyzz_0_0_1[i] * wa_x[i];

        g_xyyyzzz_0_0_0[i] = g_yyyzzz_0_0_1[i] * wa_x[i];

        g_xyyzzzz_0_0_0[i] = g_yyzzzz_0_0_1[i] * wa_x[i];

        g_xyzzzzz_0_0_0[i] = g_yzzzzz_0_0_1[i] * wa_x[i];

        g_xzzzzzz_0_0_0[i] = g_zzzzzz_0_0_1[i] * wa_x[i];

        g_yyyyyyy_0_0_0[i] = 6.0 * g_yyyyy_0_0_0[i] * fbe_0 - 6.0 * g_yyyyy_0_0_1[i] * fz_be_0 + g_yyyyyy_0_0_1[i] * wa_y[i];

        g_yyyyyyz_0_0_0[i] = g_yyyyyy_0_0_1[i] * wa_z[i];

        g_yyyyyzz_0_0_0[i] = 4.0 * g_yyyzz_0_0_0[i] * fbe_0 - 4.0 * g_yyyzz_0_0_1[i] * fz_be_0 + g_yyyyzz_0_0_1[i] * wa_y[i];

        g_yyyyzzz_0_0_0[i] = 3.0 * g_yyzzz_0_0_0[i] * fbe_0 - 3.0 * g_yyzzz_0_0_1[i] * fz_be_0 + g_yyyzzz_0_0_1[i] * wa_y[i];

        g_yyyzzzz_0_0_0[i] = 2.0 * g_yzzzz_0_0_0[i] * fbe_0 - 2.0 * g_yzzzz_0_0_1[i] * fz_be_0 + g_yyzzzz_0_0_1[i] * wa_y[i];

        g_yyzzzzz_0_0_0[i] = g_zzzzz_0_0_0[i] * fbe_0 - g_zzzzz_0_0_1[i] * fz_be_0 + g_yzzzzz_0_0_1[i] * wa_y[i];

        g_yzzzzzz_0_0_0[i] = g_zzzzzz_0_0_1[i] * wa_y[i];

        g_zzzzzzz_0_0_0[i] = 6.0 * g_zzzzz_0_0_0[i] * fbe_0 - 6.0 * g_zzzzz_0_0_1[i] * fz_be_0 + g_zzzzzz_0_0_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

