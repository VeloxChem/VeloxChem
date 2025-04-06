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

#include "TwoCenterElectronRepulsionPrimRecSI.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_si(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_si,
                                const size_t idx_eri_0_sg,
                                const size_t idx_eri_1_sg,
                                const size_t idx_eri_1_sh,
                                const CSimdArray<double>& factors,
                                const size_t idx_rpb,
                                const double a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PB) distances

    auto pb_x = factors.data(idx_rpb);

    auto pb_y = factors.data(idx_rpb + 1);

    auto pb_z = factors.data(idx_rpb + 2);

    // Set up components of auxiliary buffer : SG

    auto g_0_xxxx_0 = pbuffer.data(idx_eri_0_sg);

    auto g_0_xxyy_0 = pbuffer.data(idx_eri_0_sg + 3);

    auto g_0_xxzz_0 = pbuffer.data(idx_eri_0_sg + 5);

    auto g_0_xyyy_0 = pbuffer.data(idx_eri_0_sg + 6);

    auto g_0_xzzz_0 = pbuffer.data(idx_eri_0_sg + 9);

    auto g_0_yyyy_0 = pbuffer.data(idx_eri_0_sg + 10);

    auto g_0_yyzz_0 = pbuffer.data(idx_eri_0_sg + 12);

    auto g_0_yzzz_0 = pbuffer.data(idx_eri_0_sg + 13);

    auto g_0_zzzz_0 = pbuffer.data(idx_eri_0_sg + 14);

    // Set up components of auxiliary buffer : SG

    auto g_0_xxxx_1 = pbuffer.data(idx_eri_1_sg);

    auto g_0_xxyy_1 = pbuffer.data(idx_eri_1_sg + 3);

    auto g_0_xxzz_1 = pbuffer.data(idx_eri_1_sg + 5);

    auto g_0_xyyy_1 = pbuffer.data(idx_eri_1_sg + 6);

    auto g_0_xzzz_1 = pbuffer.data(idx_eri_1_sg + 9);

    auto g_0_yyyy_1 = pbuffer.data(idx_eri_1_sg + 10);

    auto g_0_yyzz_1 = pbuffer.data(idx_eri_1_sg + 12);

    auto g_0_yzzz_1 = pbuffer.data(idx_eri_1_sg + 13);

    auto g_0_zzzz_1 = pbuffer.data(idx_eri_1_sg + 14);

    // Set up components of auxiliary buffer : SH

    auto g_0_xxxxx_1 = pbuffer.data(idx_eri_1_sh);

    auto g_0_xxxxz_1 = pbuffer.data(idx_eri_1_sh + 2);

    auto g_0_xxxyy_1 = pbuffer.data(idx_eri_1_sh + 3);

    auto g_0_xxxzz_1 = pbuffer.data(idx_eri_1_sh + 5);

    auto g_0_xxyyy_1 = pbuffer.data(idx_eri_1_sh + 6);

    auto g_0_xxzzz_1 = pbuffer.data(idx_eri_1_sh + 9);

    auto g_0_xyyyy_1 = pbuffer.data(idx_eri_1_sh + 10);

    auto g_0_xyyzz_1 = pbuffer.data(idx_eri_1_sh + 12);

    auto g_0_xzzzz_1 = pbuffer.data(idx_eri_1_sh + 14);

    auto g_0_yyyyy_1 = pbuffer.data(idx_eri_1_sh + 15);

    auto g_0_yyyyz_1 = pbuffer.data(idx_eri_1_sh + 16);

    auto g_0_yyyzz_1 = pbuffer.data(idx_eri_1_sh + 17);

    auto g_0_yyzzz_1 = pbuffer.data(idx_eri_1_sh + 18);

    auto g_0_yzzzz_1 = pbuffer.data(idx_eri_1_sh + 19);

    auto g_0_zzzzz_1 = pbuffer.data(idx_eri_1_sh + 20);

    // Set up components of targeted buffer : SI

    auto g_0_xxxxxx_0 = pbuffer.data(idx_eri_0_si);

    auto g_0_xxxxxy_0 = pbuffer.data(idx_eri_0_si + 1);

    auto g_0_xxxxxz_0 = pbuffer.data(idx_eri_0_si + 2);

    auto g_0_xxxxyy_0 = pbuffer.data(idx_eri_0_si + 3);

    auto g_0_xxxxyz_0 = pbuffer.data(idx_eri_0_si + 4);

    auto g_0_xxxxzz_0 = pbuffer.data(idx_eri_0_si + 5);

    auto g_0_xxxyyy_0 = pbuffer.data(idx_eri_0_si + 6);

    auto g_0_xxxyyz_0 = pbuffer.data(idx_eri_0_si + 7);

    auto g_0_xxxyzz_0 = pbuffer.data(idx_eri_0_si + 8);

    auto g_0_xxxzzz_0 = pbuffer.data(idx_eri_0_si + 9);

    auto g_0_xxyyyy_0 = pbuffer.data(idx_eri_0_si + 10);

    auto g_0_xxyyyz_0 = pbuffer.data(idx_eri_0_si + 11);

    auto g_0_xxyyzz_0 = pbuffer.data(idx_eri_0_si + 12);

    auto g_0_xxyzzz_0 = pbuffer.data(idx_eri_0_si + 13);

    auto g_0_xxzzzz_0 = pbuffer.data(idx_eri_0_si + 14);

    auto g_0_xyyyyy_0 = pbuffer.data(idx_eri_0_si + 15);

    auto g_0_xyyyyz_0 = pbuffer.data(idx_eri_0_si + 16);

    auto g_0_xyyyzz_0 = pbuffer.data(idx_eri_0_si + 17);

    auto g_0_xyyzzz_0 = pbuffer.data(idx_eri_0_si + 18);

    auto g_0_xyzzzz_0 = pbuffer.data(idx_eri_0_si + 19);

    auto g_0_xzzzzz_0 = pbuffer.data(idx_eri_0_si + 20);

    auto g_0_yyyyyy_0 = pbuffer.data(idx_eri_0_si + 21);

    auto g_0_yyyyyz_0 = pbuffer.data(idx_eri_0_si + 22);

    auto g_0_yyyyzz_0 = pbuffer.data(idx_eri_0_si + 23);

    auto g_0_yyyzzz_0 = pbuffer.data(idx_eri_0_si + 24);

    auto g_0_yyzzzz_0 = pbuffer.data(idx_eri_0_si + 25);

    auto g_0_yzzzzz_0 = pbuffer.data(idx_eri_0_si + 26);

    auto g_0_zzzzzz_0 = pbuffer.data(idx_eri_0_si + 27);

    #pragma omp simd aligned(g_0_xxxx_0, g_0_xxxx_1, g_0_xxxxx_1, g_0_xxxxxx_0, g_0_xxxxxy_0, g_0_xxxxxz_0, g_0_xxxxyy_0, g_0_xxxxyz_0, g_0_xxxxz_1, g_0_xxxxzz_0, g_0_xxxyy_1, g_0_xxxyyy_0, g_0_xxxyyz_0, g_0_xxxyzz_0, g_0_xxxzz_1, g_0_xxxzzz_0, g_0_xxyy_0, g_0_xxyy_1, g_0_xxyyy_1, g_0_xxyyyy_0, g_0_xxyyyz_0, g_0_xxyyzz_0, g_0_xxyzzz_0, g_0_xxzz_0, g_0_xxzz_1, g_0_xxzzz_1, g_0_xxzzzz_0, g_0_xyyy_0, g_0_xyyy_1, g_0_xyyyy_1, g_0_xyyyyy_0, g_0_xyyyyz_0, g_0_xyyyzz_0, g_0_xyyzz_1, g_0_xyyzzz_0, g_0_xyzzzz_0, g_0_xzzz_0, g_0_xzzz_1, g_0_xzzzz_1, g_0_xzzzzz_0, g_0_yyyy_0, g_0_yyyy_1, g_0_yyyyy_1, g_0_yyyyyy_0, g_0_yyyyyz_0, g_0_yyyyz_1, g_0_yyyyzz_0, g_0_yyyzz_1, g_0_yyyzzz_0, g_0_yyzz_0, g_0_yyzz_1, g_0_yyzzz_1, g_0_yyzzzz_0, g_0_yzzz_0, g_0_yzzz_1, g_0_yzzzz_1, g_0_yzzzzz_0, g_0_zzzz_0, g_0_zzzz_1, g_0_zzzzz_1, g_0_zzzzzz_0, pb_x, pb_y, pb_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fke_0 = 0.5 / b_exps[i];

        const double fz_ke_0 = a_exp * fke_0 / (a_exp + b_exps[i]);

        g_0_xxxxxx_0[i] = 5.0 * g_0_xxxx_0[i] * fke_0 - 5.0 * g_0_xxxx_1[i] * fz_ke_0 + g_0_xxxxx_1[i] * pb_x[i];

        g_0_xxxxxy_0[i] = g_0_xxxxx_1[i] * pb_y[i];

        g_0_xxxxxz_0[i] = g_0_xxxxx_1[i] * pb_z[i];

        g_0_xxxxyy_0[i] = 3.0 * g_0_xxyy_0[i] * fke_0 - 3.0 * g_0_xxyy_1[i] * fz_ke_0 + g_0_xxxyy_1[i] * pb_x[i];

        g_0_xxxxyz_0[i] = g_0_xxxxz_1[i] * pb_y[i];

        g_0_xxxxzz_0[i] = 3.0 * g_0_xxzz_0[i] * fke_0 - 3.0 * g_0_xxzz_1[i] * fz_ke_0 + g_0_xxxzz_1[i] * pb_x[i];

        g_0_xxxyyy_0[i] = 2.0 * g_0_xyyy_0[i] * fke_0 - 2.0 * g_0_xyyy_1[i] * fz_ke_0 + g_0_xxyyy_1[i] * pb_x[i];

        g_0_xxxyyz_0[i] = g_0_xxxyy_1[i] * pb_z[i];

        g_0_xxxyzz_0[i] = g_0_xxxzz_1[i] * pb_y[i];

        g_0_xxxzzz_0[i] = 2.0 * g_0_xzzz_0[i] * fke_0 - 2.0 * g_0_xzzz_1[i] * fz_ke_0 + g_0_xxzzz_1[i] * pb_x[i];

        g_0_xxyyyy_0[i] = g_0_yyyy_0[i] * fke_0 - g_0_yyyy_1[i] * fz_ke_0 + g_0_xyyyy_1[i] * pb_x[i];

        g_0_xxyyyz_0[i] = g_0_xxyyy_1[i] * pb_z[i];

        g_0_xxyyzz_0[i] = g_0_yyzz_0[i] * fke_0 - g_0_yyzz_1[i] * fz_ke_0 + g_0_xyyzz_1[i] * pb_x[i];

        g_0_xxyzzz_0[i] = g_0_xxzzz_1[i] * pb_y[i];

        g_0_xxzzzz_0[i] = g_0_zzzz_0[i] * fke_0 - g_0_zzzz_1[i] * fz_ke_0 + g_0_xzzzz_1[i] * pb_x[i];

        g_0_xyyyyy_0[i] = g_0_yyyyy_1[i] * pb_x[i];

        g_0_xyyyyz_0[i] = g_0_yyyyz_1[i] * pb_x[i];

        g_0_xyyyzz_0[i] = g_0_yyyzz_1[i] * pb_x[i];

        g_0_xyyzzz_0[i] = g_0_yyzzz_1[i] * pb_x[i];

        g_0_xyzzzz_0[i] = g_0_yzzzz_1[i] * pb_x[i];

        g_0_xzzzzz_0[i] = g_0_zzzzz_1[i] * pb_x[i];

        g_0_yyyyyy_0[i] = 5.0 * g_0_yyyy_0[i] * fke_0 - 5.0 * g_0_yyyy_1[i] * fz_ke_0 + g_0_yyyyy_1[i] * pb_y[i];

        g_0_yyyyyz_0[i] = g_0_yyyyy_1[i] * pb_z[i];

        g_0_yyyyzz_0[i] = 3.0 * g_0_yyzz_0[i] * fke_0 - 3.0 * g_0_yyzz_1[i] * fz_ke_0 + g_0_yyyzz_1[i] * pb_y[i];

        g_0_yyyzzz_0[i] = 2.0 * g_0_yzzz_0[i] * fke_0 - 2.0 * g_0_yzzz_1[i] * fz_ke_0 + g_0_yyzzz_1[i] * pb_y[i];

        g_0_yyzzzz_0[i] = g_0_zzzz_0[i] * fke_0 - g_0_zzzz_1[i] * fz_ke_0 + g_0_yzzzz_1[i] * pb_y[i];

        g_0_yzzzzz_0[i] = g_0_zzzzz_1[i] * pb_y[i];

        g_0_zzzzzz_0[i] = 5.0 * g_0_zzzz_0[i] * fke_0 - 5.0 * g_0_zzzz_1[i] * fz_ke_0 + g_0_zzzzz_1[i] * pb_z[i];
    }
}

} // t2ceri namespace

