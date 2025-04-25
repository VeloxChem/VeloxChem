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

#include "ThreeCenterElectronRepulsionPrimRecISS.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_iss(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_iss,
                                 size_t idx_eri_0_gss,
                                 size_t idx_eri_1_gss,
                                 size_t idx_eri_1_hss,
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

    /// Set up components of auxilary buffer : GSS

    auto g_xxxx_0_0_0 = pbuffer.data(idx_eri_0_gss);

    auto g_xxyy_0_0_0 = pbuffer.data(idx_eri_0_gss + 3);

    auto g_xxzz_0_0_0 = pbuffer.data(idx_eri_0_gss + 5);

    auto g_xyyy_0_0_0 = pbuffer.data(idx_eri_0_gss + 6);

    auto g_xzzz_0_0_0 = pbuffer.data(idx_eri_0_gss + 9);

    auto g_yyyy_0_0_0 = pbuffer.data(idx_eri_0_gss + 10);

    auto g_yyzz_0_0_0 = pbuffer.data(idx_eri_0_gss + 12);

    auto g_yzzz_0_0_0 = pbuffer.data(idx_eri_0_gss + 13);

    auto g_zzzz_0_0_0 = pbuffer.data(idx_eri_0_gss + 14);

    /// Set up components of auxilary buffer : GSS

    auto g_xxxx_0_0_1 = pbuffer.data(idx_eri_1_gss);

    auto g_xxyy_0_0_1 = pbuffer.data(idx_eri_1_gss + 3);

    auto g_xxzz_0_0_1 = pbuffer.data(idx_eri_1_gss + 5);

    auto g_xyyy_0_0_1 = pbuffer.data(idx_eri_1_gss + 6);

    auto g_xzzz_0_0_1 = pbuffer.data(idx_eri_1_gss + 9);

    auto g_yyyy_0_0_1 = pbuffer.data(idx_eri_1_gss + 10);

    auto g_yyzz_0_0_1 = pbuffer.data(idx_eri_1_gss + 12);

    auto g_yzzz_0_0_1 = pbuffer.data(idx_eri_1_gss + 13);

    auto g_zzzz_0_0_1 = pbuffer.data(idx_eri_1_gss + 14);

    /// Set up components of auxilary buffer : HSS

    auto g_xxxxx_0_0_1 = pbuffer.data(idx_eri_1_hss);

    auto g_xxxxz_0_0_1 = pbuffer.data(idx_eri_1_hss + 2);

    auto g_xxxyy_0_0_1 = pbuffer.data(idx_eri_1_hss + 3);

    auto g_xxxzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 5);

    auto g_xxyyy_0_0_1 = pbuffer.data(idx_eri_1_hss + 6);

    auto g_xxzzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 9);

    auto g_xyyyy_0_0_1 = pbuffer.data(idx_eri_1_hss + 10);

    auto g_xyyzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 12);

    auto g_xzzzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 14);

    auto g_yyyyy_0_0_1 = pbuffer.data(idx_eri_1_hss + 15);

    auto g_yyyyz_0_0_1 = pbuffer.data(idx_eri_1_hss + 16);

    auto g_yyyzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 17);

    auto g_yyzzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 18);

    auto g_yzzzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 19);

    auto g_zzzzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 20);

    /// Set up components of targeted buffer : ISS

    auto g_xxxxxx_0_0_0 = pbuffer.data(idx_eri_0_iss);

    auto g_xxxxxy_0_0_0 = pbuffer.data(idx_eri_0_iss + 1);

    auto g_xxxxxz_0_0_0 = pbuffer.data(idx_eri_0_iss + 2);

    auto g_xxxxyy_0_0_0 = pbuffer.data(idx_eri_0_iss + 3);

    auto g_xxxxyz_0_0_0 = pbuffer.data(idx_eri_0_iss + 4);

    auto g_xxxxzz_0_0_0 = pbuffer.data(idx_eri_0_iss + 5);

    auto g_xxxyyy_0_0_0 = pbuffer.data(idx_eri_0_iss + 6);

    auto g_xxxyyz_0_0_0 = pbuffer.data(idx_eri_0_iss + 7);

    auto g_xxxyzz_0_0_0 = pbuffer.data(idx_eri_0_iss + 8);

    auto g_xxxzzz_0_0_0 = pbuffer.data(idx_eri_0_iss + 9);

    auto g_xxyyyy_0_0_0 = pbuffer.data(idx_eri_0_iss + 10);

    auto g_xxyyyz_0_0_0 = pbuffer.data(idx_eri_0_iss + 11);

    auto g_xxyyzz_0_0_0 = pbuffer.data(idx_eri_0_iss + 12);

    auto g_xxyzzz_0_0_0 = pbuffer.data(idx_eri_0_iss + 13);

    auto g_xxzzzz_0_0_0 = pbuffer.data(idx_eri_0_iss + 14);

    auto g_xyyyyy_0_0_0 = pbuffer.data(idx_eri_0_iss + 15);

    auto g_xyyyyz_0_0_0 = pbuffer.data(idx_eri_0_iss + 16);

    auto g_xyyyzz_0_0_0 = pbuffer.data(idx_eri_0_iss + 17);

    auto g_xyyzzz_0_0_0 = pbuffer.data(idx_eri_0_iss + 18);

    auto g_xyzzzz_0_0_0 = pbuffer.data(idx_eri_0_iss + 19);

    auto g_xzzzzz_0_0_0 = pbuffer.data(idx_eri_0_iss + 20);

    auto g_yyyyyy_0_0_0 = pbuffer.data(idx_eri_0_iss + 21);

    auto g_yyyyyz_0_0_0 = pbuffer.data(idx_eri_0_iss + 22);

    auto g_yyyyzz_0_0_0 = pbuffer.data(idx_eri_0_iss + 23);

    auto g_yyyzzz_0_0_0 = pbuffer.data(idx_eri_0_iss + 24);

    auto g_yyzzzz_0_0_0 = pbuffer.data(idx_eri_0_iss + 25);

    auto g_yzzzzz_0_0_0 = pbuffer.data(idx_eri_0_iss + 26);

    auto g_zzzzzz_0_0_0 = pbuffer.data(idx_eri_0_iss + 27);

    #pragma omp simd aligned(g_xxxx_0_0_0, g_xxxx_0_0_1, g_xxxxx_0_0_1, g_xxxxxx_0_0_0, g_xxxxxy_0_0_0, g_xxxxxz_0_0_0, g_xxxxyy_0_0_0, g_xxxxyz_0_0_0, g_xxxxz_0_0_1, g_xxxxzz_0_0_0, g_xxxyy_0_0_1, g_xxxyyy_0_0_0, g_xxxyyz_0_0_0, g_xxxyzz_0_0_0, g_xxxzz_0_0_1, g_xxxzzz_0_0_0, g_xxyy_0_0_0, g_xxyy_0_0_1, g_xxyyy_0_0_1, g_xxyyyy_0_0_0, g_xxyyyz_0_0_0, g_xxyyzz_0_0_0, g_xxyzzz_0_0_0, g_xxzz_0_0_0, g_xxzz_0_0_1, g_xxzzz_0_0_1, g_xxzzzz_0_0_0, g_xyyy_0_0_0, g_xyyy_0_0_1, g_xyyyy_0_0_1, g_xyyyyy_0_0_0, g_xyyyyz_0_0_0, g_xyyyzz_0_0_0, g_xyyzz_0_0_1, g_xyyzzz_0_0_0, g_xyzzzz_0_0_0, g_xzzz_0_0_0, g_xzzz_0_0_1, g_xzzzz_0_0_1, g_xzzzzz_0_0_0, g_yyyy_0_0_0, g_yyyy_0_0_1, g_yyyyy_0_0_1, g_yyyyyy_0_0_0, g_yyyyyz_0_0_0, g_yyyyz_0_0_1, g_yyyyzz_0_0_0, g_yyyzz_0_0_1, g_yyyzzz_0_0_0, g_yyzz_0_0_0, g_yyzz_0_0_1, g_yyzzz_0_0_1, g_yyzzzz_0_0_0, g_yzzz_0_0_0, g_yzzz_0_0_1, g_yzzzz_0_0_1, g_yzzzzz_0_0_0, g_zzzz_0_0_0, g_zzzz_0_0_1, g_zzzzz_0_0_1, g_zzzzzz_0_0_0, wa_x, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxx_0_0_0[i] = 5.0 * g_xxxx_0_0_0[i] * fbe_0 - 5.0 * g_xxxx_0_0_1[i] * fz_be_0 + g_xxxxx_0_0_1[i] * wa_x[i];

        g_xxxxxy_0_0_0[i] = g_xxxxx_0_0_1[i] * wa_y[i];

        g_xxxxxz_0_0_0[i] = g_xxxxx_0_0_1[i] * wa_z[i];

        g_xxxxyy_0_0_0[i] = 3.0 * g_xxyy_0_0_0[i] * fbe_0 - 3.0 * g_xxyy_0_0_1[i] * fz_be_0 + g_xxxyy_0_0_1[i] * wa_x[i];

        g_xxxxyz_0_0_0[i] = g_xxxxz_0_0_1[i] * wa_y[i];

        g_xxxxzz_0_0_0[i] = 3.0 * g_xxzz_0_0_0[i] * fbe_0 - 3.0 * g_xxzz_0_0_1[i] * fz_be_0 + g_xxxzz_0_0_1[i] * wa_x[i];

        g_xxxyyy_0_0_0[i] = 2.0 * g_xyyy_0_0_0[i] * fbe_0 - 2.0 * g_xyyy_0_0_1[i] * fz_be_0 + g_xxyyy_0_0_1[i] * wa_x[i];

        g_xxxyyz_0_0_0[i] = g_xxxyy_0_0_1[i] * wa_z[i];

        g_xxxyzz_0_0_0[i] = g_xxxzz_0_0_1[i] * wa_y[i];

        g_xxxzzz_0_0_0[i] = 2.0 * g_xzzz_0_0_0[i] * fbe_0 - 2.0 * g_xzzz_0_0_1[i] * fz_be_0 + g_xxzzz_0_0_1[i] * wa_x[i];

        g_xxyyyy_0_0_0[i] = g_yyyy_0_0_0[i] * fbe_0 - g_yyyy_0_0_1[i] * fz_be_0 + g_xyyyy_0_0_1[i] * wa_x[i];

        g_xxyyyz_0_0_0[i] = g_xxyyy_0_0_1[i] * wa_z[i];

        g_xxyyzz_0_0_0[i] = g_yyzz_0_0_0[i] * fbe_0 - g_yyzz_0_0_1[i] * fz_be_0 + g_xyyzz_0_0_1[i] * wa_x[i];

        g_xxyzzz_0_0_0[i] = g_xxzzz_0_0_1[i] * wa_y[i];

        g_xxzzzz_0_0_0[i] = g_zzzz_0_0_0[i] * fbe_0 - g_zzzz_0_0_1[i] * fz_be_0 + g_xzzzz_0_0_1[i] * wa_x[i];

        g_xyyyyy_0_0_0[i] = g_yyyyy_0_0_1[i] * wa_x[i];

        g_xyyyyz_0_0_0[i] = g_yyyyz_0_0_1[i] * wa_x[i];

        g_xyyyzz_0_0_0[i] = g_yyyzz_0_0_1[i] * wa_x[i];

        g_xyyzzz_0_0_0[i] = g_yyzzz_0_0_1[i] * wa_x[i];

        g_xyzzzz_0_0_0[i] = g_yzzzz_0_0_1[i] * wa_x[i];

        g_xzzzzz_0_0_0[i] = g_zzzzz_0_0_1[i] * wa_x[i];

        g_yyyyyy_0_0_0[i] = 5.0 * g_yyyy_0_0_0[i] * fbe_0 - 5.0 * g_yyyy_0_0_1[i] * fz_be_0 + g_yyyyy_0_0_1[i] * wa_y[i];

        g_yyyyyz_0_0_0[i] = g_yyyyy_0_0_1[i] * wa_z[i];

        g_yyyyzz_0_0_0[i] = 3.0 * g_yyzz_0_0_0[i] * fbe_0 - 3.0 * g_yyzz_0_0_1[i] * fz_be_0 + g_yyyzz_0_0_1[i] * wa_y[i];

        g_yyyzzz_0_0_0[i] = 2.0 * g_yzzz_0_0_0[i] * fbe_0 - 2.0 * g_yzzz_0_0_1[i] * fz_be_0 + g_yyzzz_0_0_1[i] * wa_y[i];

        g_yyzzzz_0_0_0[i] = g_zzzz_0_0_0[i] * fbe_0 - g_zzzz_0_0_1[i] * fz_be_0 + g_yzzzz_0_0_1[i] * wa_y[i];

        g_yzzzzz_0_0_0[i] = g_zzzzz_0_0_1[i] * wa_y[i];

        g_zzzzzz_0_0_0[i] = 5.0 * g_zzzz_0_0_0[i] * fbe_0 - 5.0 * g_zzzz_0_0_1[i] * fz_be_0 + g_zzzzz_0_0_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

