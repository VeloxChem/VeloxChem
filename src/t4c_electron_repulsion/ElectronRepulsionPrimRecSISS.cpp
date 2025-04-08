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

#include "ElectronRepulsionPrimRecSISS.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_siss(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_siss,
                                  size_t                idx_eri_0_sgss,
                                  size_t                idx_eri_1_sgss,
                                  size_t                idx_eri_0_shss,
                                  size_t                idx_eri_1_shss,
                                  CSimdArray<double>&   factors,
                                  const size_t          idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double          a_exp,
                                  const double          b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SGSS

    auto g_0_xxxx_0_0_0 = pbuffer.data(idx_eri_0_sgss);

    auto g_0_xxyy_0_0_0 = pbuffer.data(idx_eri_0_sgss + 3);

    auto g_0_xxzz_0_0_0 = pbuffer.data(idx_eri_0_sgss + 5);

    auto g_0_xyyy_0_0_0 = pbuffer.data(idx_eri_0_sgss + 6);

    auto g_0_xzzz_0_0_0 = pbuffer.data(idx_eri_0_sgss + 9);

    auto g_0_yyyy_0_0_0 = pbuffer.data(idx_eri_0_sgss + 10);

    auto g_0_yyzz_0_0_0 = pbuffer.data(idx_eri_0_sgss + 12);

    auto g_0_yzzz_0_0_0 = pbuffer.data(idx_eri_0_sgss + 13);

    auto g_0_zzzz_0_0_0 = pbuffer.data(idx_eri_0_sgss + 14);

    /// Set up components of auxilary buffer : SGSS

    auto g_0_xxxx_0_0_1 = pbuffer.data(idx_eri_1_sgss);

    auto g_0_xxyy_0_0_1 = pbuffer.data(idx_eri_1_sgss + 3);

    auto g_0_xxzz_0_0_1 = pbuffer.data(idx_eri_1_sgss + 5);

    auto g_0_xyyy_0_0_1 = pbuffer.data(idx_eri_1_sgss + 6);

    auto g_0_xzzz_0_0_1 = pbuffer.data(idx_eri_1_sgss + 9);

    auto g_0_yyyy_0_0_1 = pbuffer.data(idx_eri_1_sgss + 10);

    auto g_0_yyzz_0_0_1 = pbuffer.data(idx_eri_1_sgss + 12);

    auto g_0_yzzz_0_0_1 = pbuffer.data(idx_eri_1_sgss + 13);

    auto g_0_zzzz_0_0_1 = pbuffer.data(idx_eri_1_sgss + 14);

    /// Set up components of auxilary buffer : SHSS

    auto g_0_xxxxx_0_0_0 = pbuffer.data(idx_eri_0_shss);

    auto g_0_xxxxz_0_0_0 = pbuffer.data(idx_eri_0_shss + 2);

    auto g_0_xxxyy_0_0_0 = pbuffer.data(idx_eri_0_shss + 3);

    auto g_0_xxxzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 5);

    auto g_0_xxyyy_0_0_0 = pbuffer.data(idx_eri_0_shss + 6);

    auto g_0_xxzzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 9);

    auto g_0_xyyyy_0_0_0 = pbuffer.data(idx_eri_0_shss + 10);

    auto g_0_xyyzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 12);

    auto g_0_xzzzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 14);

    auto g_0_yyyyy_0_0_0 = pbuffer.data(idx_eri_0_shss + 15);

    auto g_0_yyyyz_0_0_0 = pbuffer.data(idx_eri_0_shss + 16);

    auto g_0_yyyzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 17);

    auto g_0_yyzzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 18);

    auto g_0_yzzzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 19);

    auto g_0_zzzzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 20);

    /// Set up components of auxilary buffer : SHSS

    auto g_0_xxxxx_0_0_1 = pbuffer.data(idx_eri_1_shss);

    auto g_0_xxxxz_0_0_1 = pbuffer.data(idx_eri_1_shss + 2);

    auto g_0_xxxyy_0_0_1 = pbuffer.data(idx_eri_1_shss + 3);

    auto g_0_xxxzz_0_0_1 = pbuffer.data(idx_eri_1_shss + 5);

    auto g_0_xxyyy_0_0_1 = pbuffer.data(idx_eri_1_shss + 6);

    auto g_0_xxzzz_0_0_1 = pbuffer.data(idx_eri_1_shss + 9);

    auto g_0_xyyyy_0_0_1 = pbuffer.data(idx_eri_1_shss + 10);

    auto g_0_xyyzz_0_0_1 = pbuffer.data(idx_eri_1_shss + 12);

    auto g_0_xzzzz_0_0_1 = pbuffer.data(idx_eri_1_shss + 14);

    auto g_0_yyyyy_0_0_1 = pbuffer.data(idx_eri_1_shss + 15);

    auto g_0_yyyyz_0_0_1 = pbuffer.data(idx_eri_1_shss + 16);

    auto g_0_yyyzz_0_0_1 = pbuffer.data(idx_eri_1_shss + 17);

    auto g_0_yyzzz_0_0_1 = pbuffer.data(idx_eri_1_shss + 18);

    auto g_0_yzzzz_0_0_1 = pbuffer.data(idx_eri_1_shss + 19);

    auto g_0_zzzzz_0_0_1 = pbuffer.data(idx_eri_1_shss + 20);

    /// Set up components of targeted buffer : SISS

    auto g_0_xxxxxx_0_0_0 = pbuffer.data(idx_eri_0_siss);

    auto g_0_xxxxxy_0_0_0 = pbuffer.data(idx_eri_0_siss + 1);

    auto g_0_xxxxxz_0_0_0 = pbuffer.data(idx_eri_0_siss + 2);

    auto g_0_xxxxyy_0_0_0 = pbuffer.data(idx_eri_0_siss + 3);

    auto g_0_xxxxyz_0_0_0 = pbuffer.data(idx_eri_0_siss + 4);

    auto g_0_xxxxzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 5);

    auto g_0_xxxyyy_0_0_0 = pbuffer.data(idx_eri_0_siss + 6);

    auto g_0_xxxyyz_0_0_0 = pbuffer.data(idx_eri_0_siss + 7);

    auto g_0_xxxyzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 8);

    auto g_0_xxxzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 9);

    auto g_0_xxyyyy_0_0_0 = pbuffer.data(idx_eri_0_siss + 10);

    auto g_0_xxyyyz_0_0_0 = pbuffer.data(idx_eri_0_siss + 11);

    auto g_0_xxyyzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 12);

    auto g_0_xxyzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 13);

    auto g_0_xxzzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 14);

    auto g_0_xyyyyy_0_0_0 = pbuffer.data(idx_eri_0_siss + 15);

    auto g_0_xyyyyz_0_0_0 = pbuffer.data(idx_eri_0_siss + 16);

    auto g_0_xyyyzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 17);

    auto g_0_xyyzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 18);

    auto g_0_xyzzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 19);

    auto g_0_xzzzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 20);

    auto g_0_yyyyyy_0_0_0 = pbuffer.data(idx_eri_0_siss + 21);

    auto g_0_yyyyyz_0_0_0 = pbuffer.data(idx_eri_0_siss + 22);

    auto g_0_yyyyzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 23);

    auto g_0_yyyzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 24);

    auto g_0_yyzzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 25);

    auto g_0_yzzzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 26);

    auto g_0_zzzzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 27);

#pragma omp simd aligned(g_0_xxxx_0_0_0,       \
                             g_0_xxxx_0_0_1,   \
                             g_0_xxxxx_0_0_0,  \
                             g_0_xxxxx_0_0_1,  \
                             g_0_xxxxxx_0_0_0, \
                             g_0_xxxxxy_0_0_0, \
                             g_0_xxxxxz_0_0_0, \
                             g_0_xxxxyy_0_0_0, \
                             g_0_xxxxyz_0_0_0, \
                             g_0_xxxxz_0_0_0,  \
                             g_0_xxxxz_0_0_1,  \
                             g_0_xxxxzz_0_0_0, \
                             g_0_xxxyy_0_0_0,  \
                             g_0_xxxyy_0_0_1,  \
                             g_0_xxxyyy_0_0_0, \
                             g_0_xxxyyz_0_0_0, \
                             g_0_xxxyzz_0_0_0, \
                             g_0_xxxzz_0_0_0,  \
                             g_0_xxxzz_0_0_1,  \
                             g_0_xxxzzz_0_0_0, \
                             g_0_xxyy_0_0_0,   \
                             g_0_xxyy_0_0_1,   \
                             g_0_xxyyy_0_0_0,  \
                             g_0_xxyyy_0_0_1,  \
                             g_0_xxyyyy_0_0_0, \
                             g_0_xxyyyz_0_0_0, \
                             g_0_xxyyzz_0_0_0, \
                             g_0_xxyzzz_0_0_0, \
                             g_0_xxzz_0_0_0,   \
                             g_0_xxzz_0_0_1,   \
                             g_0_xxzzz_0_0_0,  \
                             g_0_xxzzz_0_0_1,  \
                             g_0_xxzzzz_0_0_0, \
                             g_0_xyyy_0_0_0,   \
                             g_0_xyyy_0_0_1,   \
                             g_0_xyyyy_0_0_0,  \
                             g_0_xyyyy_0_0_1,  \
                             g_0_xyyyyy_0_0_0, \
                             g_0_xyyyyz_0_0_0, \
                             g_0_xyyyzz_0_0_0, \
                             g_0_xyyzz_0_0_0,  \
                             g_0_xyyzz_0_0_1,  \
                             g_0_xyyzzz_0_0_0, \
                             g_0_xyzzzz_0_0_0, \
                             g_0_xzzz_0_0_0,   \
                             g_0_xzzz_0_0_1,   \
                             g_0_xzzzz_0_0_0,  \
                             g_0_xzzzz_0_0_1,  \
                             g_0_xzzzzz_0_0_0, \
                             g_0_yyyy_0_0_0,   \
                             g_0_yyyy_0_0_1,   \
                             g_0_yyyyy_0_0_0,  \
                             g_0_yyyyy_0_0_1,  \
                             g_0_yyyyyy_0_0_0, \
                             g_0_yyyyyz_0_0_0, \
                             g_0_yyyyz_0_0_0,  \
                             g_0_yyyyz_0_0_1,  \
                             g_0_yyyyzz_0_0_0, \
                             g_0_yyyzz_0_0_0,  \
                             g_0_yyyzz_0_0_1,  \
                             g_0_yyyzzz_0_0_0, \
                             g_0_yyzz_0_0_0,   \
                             g_0_yyzz_0_0_1,   \
                             g_0_yyzzz_0_0_0,  \
                             g_0_yyzzz_0_0_1,  \
                             g_0_yyzzzz_0_0_0, \
                             g_0_yzzz_0_0_0,   \
                             g_0_yzzz_0_0_1,   \
                             g_0_yzzzz_0_0_0,  \
                             g_0_yzzzz_0_0_1,  \
                             g_0_yzzzzz_0_0_0, \
                             g_0_zzzz_0_0_0,   \
                             g_0_zzzz_0_0_1,   \
                             g_0_zzzzz_0_0_0,  \
                             g_0_zzzzz_0_0_1,  \
                             g_0_zzzzzz_0_0_0, \
                             wp_x,             \
                             wp_y,             \
                             wp_z,             \
                             c_exps,           \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxx_0_0_0[i] =
            5.0 * g_0_xxxx_0_0_0[i] * fi_ab_0 - 5.0 * g_0_xxxx_0_0_1[i] * fti_ab_0 + g_0_xxxxx_0_0_0[i] * pb_x + g_0_xxxxx_0_0_1[i] * wp_x[i];

        g_0_xxxxxy_0_0_0[i] = g_0_xxxxx_0_0_0[i] * pb_y + g_0_xxxxx_0_0_1[i] * wp_y[i];

        g_0_xxxxxz_0_0_0[i] = g_0_xxxxx_0_0_0[i] * pb_z + g_0_xxxxx_0_0_1[i] * wp_z[i];

        g_0_xxxxyy_0_0_0[i] =
            3.0 * g_0_xxyy_0_0_0[i] * fi_ab_0 - 3.0 * g_0_xxyy_0_0_1[i] * fti_ab_0 + g_0_xxxyy_0_0_0[i] * pb_x + g_0_xxxyy_0_0_1[i] * wp_x[i];

        g_0_xxxxyz_0_0_0[i] = g_0_xxxxz_0_0_0[i] * pb_y + g_0_xxxxz_0_0_1[i] * wp_y[i];

        g_0_xxxxzz_0_0_0[i] =
            3.0 * g_0_xxzz_0_0_0[i] * fi_ab_0 - 3.0 * g_0_xxzz_0_0_1[i] * fti_ab_0 + g_0_xxxzz_0_0_0[i] * pb_x + g_0_xxxzz_0_0_1[i] * wp_x[i];

        g_0_xxxyyy_0_0_0[i] =
            2.0 * g_0_xyyy_0_0_0[i] * fi_ab_0 - 2.0 * g_0_xyyy_0_0_1[i] * fti_ab_0 + g_0_xxyyy_0_0_0[i] * pb_x + g_0_xxyyy_0_0_1[i] * wp_x[i];

        g_0_xxxyyz_0_0_0[i] = g_0_xxxyy_0_0_0[i] * pb_z + g_0_xxxyy_0_0_1[i] * wp_z[i];

        g_0_xxxyzz_0_0_0[i] = g_0_xxxzz_0_0_0[i] * pb_y + g_0_xxxzz_0_0_1[i] * wp_y[i];

        g_0_xxxzzz_0_0_0[i] =
            2.0 * g_0_xzzz_0_0_0[i] * fi_ab_0 - 2.0 * g_0_xzzz_0_0_1[i] * fti_ab_0 + g_0_xxzzz_0_0_0[i] * pb_x + g_0_xxzzz_0_0_1[i] * wp_x[i];

        g_0_xxyyyy_0_0_0[i] = g_0_yyyy_0_0_0[i] * fi_ab_0 - g_0_yyyy_0_0_1[i] * fti_ab_0 + g_0_xyyyy_0_0_0[i] * pb_x + g_0_xyyyy_0_0_1[i] * wp_x[i];

        g_0_xxyyyz_0_0_0[i] = g_0_xxyyy_0_0_0[i] * pb_z + g_0_xxyyy_0_0_1[i] * wp_z[i];

        g_0_xxyyzz_0_0_0[i] = g_0_yyzz_0_0_0[i] * fi_ab_0 - g_0_yyzz_0_0_1[i] * fti_ab_0 + g_0_xyyzz_0_0_0[i] * pb_x + g_0_xyyzz_0_0_1[i] * wp_x[i];

        g_0_xxyzzz_0_0_0[i] = g_0_xxzzz_0_0_0[i] * pb_y + g_0_xxzzz_0_0_1[i] * wp_y[i];

        g_0_xxzzzz_0_0_0[i] = g_0_zzzz_0_0_0[i] * fi_ab_0 - g_0_zzzz_0_0_1[i] * fti_ab_0 + g_0_xzzzz_0_0_0[i] * pb_x + g_0_xzzzz_0_0_1[i] * wp_x[i];

        g_0_xyyyyy_0_0_0[i] = g_0_yyyyy_0_0_0[i] * pb_x + g_0_yyyyy_0_0_1[i] * wp_x[i];

        g_0_xyyyyz_0_0_0[i] = g_0_yyyyz_0_0_0[i] * pb_x + g_0_yyyyz_0_0_1[i] * wp_x[i];

        g_0_xyyyzz_0_0_0[i] = g_0_yyyzz_0_0_0[i] * pb_x + g_0_yyyzz_0_0_1[i] * wp_x[i];

        g_0_xyyzzz_0_0_0[i] = g_0_yyzzz_0_0_0[i] * pb_x + g_0_yyzzz_0_0_1[i] * wp_x[i];

        g_0_xyzzzz_0_0_0[i] = g_0_yzzzz_0_0_0[i] * pb_x + g_0_yzzzz_0_0_1[i] * wp_x[i];

        g_0_xzzzzz_0_0_0[i] = g_0_zzzzz_0_0_0[i] * pb_x + g_0_zzzzz_0_0_1[i] * wp_x[i];

        g_0_yyyyyy_0_0_0[i] =
            5.0 * g_0_yyyy_0_0_0[i] * fi_ab_0 - 5.0 * g_0_yyyy_0_0_1[i] * fti_ab_0 + g_0_yyyyy_0_0_0[i] * pb_y + g_0_yyyyy_0_0_1[i] * wp_y[i];

        g_0_yyyyyz_0_0_0[i] = g_0_yyyyy_0_0_0[i] * pb_z + g_0_yyyyy_0_0_1[i] * wp_z[i];

        g_0_yyyyzz_0_0_0[i] =
            3.0 * g_0_yyzz_0_0_0[i] * fi_ab_0 - 3.0 * g_0_yyzz_0_0_1[i] * fti_ab_0 + g_0_yyyzz_0_0_0[i] * pb_y + g_0_yyyzz_0_0_1[i] * wp_y[i];

        g_0_yyyzzz_0_0_0[i] =
            2.0 * g_0_yzzz_0_0_0[i] * fi_ab_0 - 2.0 * g_0_yzzz_0_0_1[i] * fti_ab_0 + g_0_yyzzz_0_0_0[i] * pb_y + g_0_yyzzz_0_0_1[i] * wp_y[i];

        g_0_yyzzzz_0_0_0[i] = g_0_zzzz_0_0_0[i] * fi_ab_0 - g_0_zzzz_0_0_1[i] * fti_ab_0 + g_0_yzzzz_0_0_0[i] * pb_y + g_0_yzzzz_0_0_1[i] * wp_y[i];

        g_0_yzzzzz_0_0_0[i] = g_0_zzzzz_0_0_0[i] * pb_y + g_0_zzzzz_0_0_1[i] * wp_y[i];

        g_0_zzzzzz_0_0_0[i] =
            5.0 * g_0_zzzz_0_0_0[i] * fi_ab_0 - 5.0 * g_0_zzzz_0_0_1[i] * fti_ab_0 + g_0_zzzzz_0_0_0[i] * pb_z + g_0_zzzzz_0_0_1[i] * wp_z[i];
    }
}

}  // namespace erirec
