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

#include "ElectronRepulsionPrimRecSKSS.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_skss(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_skss,
                                  size_t                idx_eri_0_shss,
                                  size_t                idx_eri_1_shss,
                                  size_t                idx_eri_0_siss,
                                  size_t                idx_eri_1_siss,
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

    /// Set up components of auxilary buffer : SHSS

    auto g_0_xxxxx_0_0_0 = pbuffer.data(idx_eri_0_shss);

    auto g_0_xxxyy_0_0_0 = pbuffer.data(idx_eri_0_shss + 3);

    auto g_0_xxxzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 5);

    auto g_0_xxyyy_0_0_0 = pbuffer.data(idx_eri_0_shss + 6);

    auto g_0_xxzzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 9);

    auto g_0_xyyyy_0_0_0 = pbuffer.data(idx_eri_0_shss + 10);

    auto g_0_xyyzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 12);

    auto g_0_xzzzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 14);

    auto g_0_yyyyy_0_0_0 = pbuffer.data(idx_eri_0_shss + 15);

    auto g_0_yyyzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 17);

    auto g_0_yyzzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 18);

    auto g_0_yzzzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 19);

    auto g_0_zzzzz_0_0_0 = pbuffer.data(idx_eri_0_shss + 20);

    /// Set up components of auxilary buffer : SHSS

    auto g_0_xxxxx_0_0_1 = pbuffer.data(idx_eri_1_shss);

    auto g_0_xxxyy_0_0_1 = pbuffer.data(idx_eri_1_shss + 3);

    auto g_0_xxxzz_0_0_1 = pbuffer.data(idx_eri_1_shss + 5);

    auto g_0_xxyyy_0_0_1 = pbuffer.data(idx_eri_1_shss + 6);

    auto g_0_xxzzz_0_0_1 = pbuffer.data(idx_eri_1_shss + 9);

    auto g_0_xyyyy_0_0_1 = pbuffer.data(idx_eri_1_shss + 10);

    auto g_0_xyyzz_0_0_1 = pbuffer.data(idx_eri_1_shss + 12);

    auto g_0_xzzzz_0_0_1 = pbuffer.data(idx_eri_1_shss + 14);

    auto g_0_yyyyy_0_0_1 = pbuffer.data(idx_eri_1_shss + 15);

    auto g_0_yyyzz_0_0_1 = pbuffer.data(idx_eri_1_shss + 17);

    auto g_0_yyzzz_0_0_1 = pbuffer.data(idx_eri_1_shss + 18);

    auto g_0_yzzzz_0_0_1 = pbuffer.data(idx_eri_1_shss + 19);

    auto g_0_zzzzz_0_0_1 = pbuffer.data(idx_eri_1_shss + 20);

    /// Set up components of auxilary buffer : SISS

    auto g_0_xxxxxx_0_0_0 = pbuffer.data(idx_eri_0_siss);

    auto g_0_xxxxxz_0_0_0 = pbuffer.data(idx_eri_0_siss + 2);

    auto g_0_xxxxyy_0_0_0 = pbuffer.data(idx_eri_0_siss + 3);

    auto g_0_xxxxzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 5);

    auto g_0_xxxyyy_0_0_0 = pbuffer.data(idx_eri_0_siss + 6);

    auto g_0_xxxzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 9);

    auto g_0_xxyyyy_0_0_0 = pbuffer.data(idx_eri_0_siss + 10);

    auto g_0_xxyyzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 12);

    auto g_0_xxzzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 14);

    auto g_0_xyyyyy_0_0_0 = pbuffer.data(idx_eri_0_siss + 15);

    auto g_0_xyyyzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 17);

    auto g_0_xyyzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 18);

    auto g_0_xzzzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 20);

    auto g_0_yyyyyy_0_0_0 = pbuffer.data(idx_eri_0_siss + 21);

    auto g_0_yyyyyz_0_0_0 = pbuffer.data(idx_eri_0_siss + 22);

    auto g_0_yyyyzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 23);

    auto g_0_yyyzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 24);

    auto g_0_yyzzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 25);

    auto g_0_yzzzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 26);

    auto g_0_zzzzzz_0_0_0 = pbuffer.data(idx_eri_0_siss + 27);

    /// Set up components of auxilary buffer : SISS

    auto g_0_xxxxxx_0_0_1 = pbuffer.data(idx_eri_1_siss);

    auto g_0_xxxxxz_0_0_1 = pbuffer.data(idx_eri_1_siss + 2);

    auto g_0_xxxxyy_0_0_1 = pbuffer.data(idx_eri_1_siss + 3);

    auto g_0_xxxxzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 5);

    auto g_0_xxxyyy_0_0_1 = pbuffer.data(idx_eri_1_siss + 6);

    auto g_0_xxxzzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 9);

    auto g_0_xxyyyy_0_0_1 = pbuffer.data(idx_eri_1_siss + 10);

    auto g_0_xxyyzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 12);

    auto g_0_xxzzzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 14);

    auto g_0_xyyyyy_0_0_1 = pbuffer.data(idx_eri_1_siss + 15);

    auto g_0_xyyyzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 17);

    auto g_0_xyyzzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 18);

    auto g_0_xzzzzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 20);

    auto g_0_yyyyyy_0_0_1 = pbuffer.data(idx_eri_1_siss + 21);

    auto g_0_yyyyyz_0_0_1 = pbuffer.data(idx_eri_1_siss + 22);

    auto g_0_yyyyzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 23);

    auto g_0_yyyzzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 24);

    auto g_0_yyzzzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 25);

    auto g_0_yzzzzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 26);

    auto g_0_zzzzzz_0_0_1 = pbuffer.data(idx_eri_1_siss + 27);

    /// Set up components of targeted buffer : SKSS

    auto g_0_xxxxxxx_0_0_0 = pbuffer.data(idx_eri_0_skss);

    auto g_0_xxxxxxy_0_0_0 = pbuffer.data(idx_eri_0_skss + 1);

    auto g_0_xxxxxxz_0_0_0 = pbuffer.data(idx_eri_0_skss + 2);

    auto g_0_xxxxxyy_0_0_0 = pbuffer.data(idx_eri_0_skss + 3);

    auto g_0_xxxxxyz_0_0_0 = pbuffer.data(idx_eri_0_skss + 4);

    auto g_0_xxxxxzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 5);

    auto g_0_xxxxyyy_0_0_0 = pbuffer.data(idx_eri_0_skss + 6);

    auto g_0_xxxxyyz_0_0_0 = pbuffer.data(idx_eri_0_skss + 7);

    auto g_0_xxxxyzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 8);

    auto g_0_xxxxzzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 9);

    auto g_0_xxxyyyy_0_0_0 = pbuffer.data(idx_eri_0_skss + 10);

    auto g_0_xxxyyyz_0_0_0 = pbuffer.data(idx_eri_0_skss + 11);

    auto g_0_xxxyyzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 12);

    auto g_0_xxxyzzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 13);

    auto g_0_xxxzzzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 14);

    auto g_0_xxyyyyy_0_0_0 = pbuffer.data(idx_eri_0_skss + 15);

    auto g_0_xxyyyyz_0_0_0 = pbuffer.data(idx_eri_0_skss + 16);

    auto g_0_xxyyyzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 17);

    auto g_0_xxyyzzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 18);

    auto g_0_xxyzzzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 19);

    auto g_0_xxzzzzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 20);

    auto g_0_xyyyyyy_0_0_0 = pbuffer.data(idx_eri_0_skss + 21);

    auto g_0_xyyyyyz_0_0_0 = pbuffer.data(idx_eri_0_skss + 22);

    auto g_0_xyyyyzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 23);

    auto g_0_xyyyzzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 24);

    auto g_0_xyyzzzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 25);

    auto g_0_xyzzzzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 26);

    auto g_0_xzzzzzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 27);

    auto g_0_yyyyyyy_0_0_0 = pbuffer.data(idx_eri_0_skss + 28);

    auto g_0_yyyyyyz_0_0_0 = pbuffer.data(idx_eri_0_skss + 29);

    auto g_0_yyyyyzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 30);

    auto g_0_yyyyzzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 31);

    auto g_0_yyyzzzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 32);

    auto g_0_yyzzzzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 33);

    auto g_0_yzzzzzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 34);

    auto g_0_zzzzzzz_0_0_0 = pbuffer.data(idx_eri_0_skss + 35);

#pragma omp simd aligned(g_0_xxxxx_0_0_0,       \
                             g_0_xxxxx_0_0_1,   \
                             g_0_xxxxxx_0_0_0,  \
                             g_0_xxxxxx_0_0_1,  \
                             g_0_xxxxxxx_0_0_0, \
                             g_0_xxxxxxy_0_0_0, \
                             g_0_xxxxxxz_0_0_0, \
                             g_0_xxxxxyy_0_0_0, \
                             g_0_xxxxxyz_0_0_0, \
                             g_0_xxxxxz_0_0_0,  \
                             g_0_xxxxxz_0_0_1,  \
                             g_0_xxxxxzz_0_0_0, \
                             g_0_xxxxyy_0_0_0,  \
                             g_0_xxxxyy_0_0_1,  \
                             g_0_xxxxyyy_0_0_0, \
                             g_0_xxxxyyz_0_0_0, \
                             g_0_xxxxyzz_0_0_0, \
                             g_0_xxxxzz_0_0_0,  \
                             g_0_xxxxzz_0_0_1,  \
                             g_0_xxxxzzz_0_0_0, \
                             g_0_xxxyy_0_0_0,   \
                             g_0_xxxyy_0_0_1,   \
                             g_0_xxxyyy_0_0_0,  \
                             g_0_xxxyyy_0_0_1,  \
                             g_0_xxxyyyy_0_0_0, \
                             g_0_xxxyyyz_0_0_0, \
                             g_0_xxxyyzz_0_0_0, \
                             g_0_xxxyzzz_0_0_0, \
                             g_0_xxxzz_0_0_0,   \
                             g_0_xxxzz_0_0_1,   \
                             g_0_xxxzzz_0_0_0,  \
                             g_0_xxxzzz_0_0_1,  \
                             g_0_xxxzzzz_0_0_0, \
                             g_0_xxyyy_0_0_0,   \
                             g_0_xxyyy_0_0_1,   \
                             g_0_xxyyyy_0_0_0,  \
                             g_0_xxyyyy_0_0_1,  \
                             g_0_xxyyyyy_0_0_0, \
                             g_0_xxyyyyz_0_0_0, \
                             g_0_xxyyyzz_0_0_0, \
                             g_0_xxyyzz_0_0_0,  \
                             g_0_xxyyzz_0_0_1,  \
                             g_0_xxyyzzz_0_0_0, \
                             g_0_xxyzzzz_0_0_0, \
                             g_0_xxzzz_0_0_0,   \
                             g_0_xxzzz_0_0_1,   \
                             g_0_xxzzzz_0_0_0,  \
                             g_0_xxzzzz_0_0_1,  \
                             g_0_xxzzzzz_0_0_0, \
                             g_0_xyyyy_0_0_0,   \
                             g_0_xyyyy_0_0_1,   \
                             g_0_xyyyyy_0_0_0,  \
                             g_0_xyyyyy_0_0_1,  \
                             g_0_xyyyyyy_0_0_0, \
                             g_0_xyyyyyz_0_0_0, \
                             g_0_xyyyyzz_0_0_0, \
                             g_0_xyyyzz_0_0_0,  \
                             g_0_xyyyzz_0_0_1,  \
                             g_0_xyyyzzz_0_0_0, \
                             g_0_xyyzz_0_0_0,   \
                             g_0_xyyzz_0_0_1,   \
                             g_0_xyyzzz_0_0_0,  \
                             g_0_xyyzzz_0_0_1,  \
                             g_0_xyyzzzz_0_0_0, \
                             g_0_xyzzzzz_0_0_0, \
                             g_0_xzzzz_0_0_0,   \
                             g_0_xzzzz_0_0_1,   \
                             g_0_xzzzzz_0_0_0,  \
                             g_0_xzzzzz_0_0_1,  \
                             g_0_xzzzzzz_0_0_0, \
                             g_0_yyyyy_0_0_0,   \
                             g_0_yyyyy_0_0_1,   \
                             g_0_yyyyyy_0_0_0,  \
                             g_0_yyyyyy_0_0_1,  \
                             g_0_yyyyyyy_0_0_0, \
                             g_0_yyyyyyz_0_0_0, \
                             g_0_yyyyyz_0_0_0,  \
                             g_0_yyyyyz_0_0_1,  \
                             g_0_yyyyyzz_0_0_0, \
                             g_0_yyyyzz_0_0_0,  \
                             g_0_yyyyzz_0_0_1,  \
                             g_0_yyyyzzz_0_0_0, \
                             g_0_yyyzz_0_0_0,   \
                             g_0_yyyzz_0_0_1,   \
                             g_0_yyyzzz_0_0_0,  \
                             g_0_yyyzzz_0_0_1,  \
                             g_0_yyyzzzz_0_0_0, \
                             g_0_yyzzz_0_0_0,   \
                             g_0_yyzzz_0_0_1,   \
                             g_0_yyzzzz_0_0_0,  \
                             g_0_yyzzzz_0_0_1,  \
                             g_0_yyzzzzz_0_0_0, \
                             g_0_yzzzz_0_0_0,   \
                             g_0_yzzzz_0_0_1,   \
                             g_0_yzzzzz_0_0_0,  \
                             g_0_yzzzzz_0_0_1,  \
                             g_0_yzzzzzz_0_0_0, \
                             g_0_zzzzz_0_0_0,   \
                             g_0_zzzzz_0_0_1,   \
                             g_0_zzzzzz_0_0_0,  \
                             g_0_zzzzzz_0_0_1,  \
                             g_0_zzzzzzz_0_0_0, \
                             wp_x,              \
                             wp_y,              \
                             wp_z,              \
                             c_exps,            \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fti_ab_0 = fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxxxxxx_0_0_0[i] =
            6.0 * g_0_xxxxx_0_0_0[i] * fi_ab_0 - 6.0 * g_0_xxxxx_0_0_1[i] * fti_ab_0 + g_0_xxxxxx_0_0_0[i] * pb_x + g_0_xxxxxx_0_0_1[i] * wp_x[i];

        g_0_xxxxxxy_0_0_0[i] = g_0_xxxxxx_0_0_0[i] * pb_y + g_0_xxxxxx_0_0_1[i] * wp_y[i];

        g_0_xxxxxxz_0_0_0[i] = g_0_xxxxxx_0_0_0[i] * pb_z + g_0_xxxxxx_0_0_1[i] * wp_z[i];

        g_0_xxxxxyy_0_0_0[i] =
            4.0 * g_0_xxxyy_0_0_0[i] * fi_ab_0 - 4.0 * g_0_xxxyy_0_0_1[i] * fti_ab_0 + g_0_xxxxyy_0_0_0[i] * pb_x + g_0_xxxxyy_0_0_1[i] * wp_x[i];

        g_0_xxxxxyz_0_0_0[i] = g_0_xxxxxz_0_0_0[i] * pb_y + g_0_xxxxxz_0_0_1[i] * wp_y[i];

        g_0_xxxxxzz_0_0_0[i] =
            4.0 * g_0_xxxzz_0_0_0[i] * fi_ab_0 - 4.0 * g_0_xxxzz_0_0_1[i] * fti_ab_0 + g_0_xxxxzz_0_0_0[i] * pb_x + g_0_xxxxzz_0_0_1[i] * wp_x[i];

        g_0_xxxxyyy_0_0_0[i] =
            3.0 * g_0_xxyyy_0_0_0[i] * fi_ab_0 - 3.0 * g_0_xxyyy_0_0_1[i] * fti_ab_0 + g_0_xxxyyy_0_0_0[i] * pb_x + g_0_xxxyyy_0_0_1[i] * wp_x[i];

        g_0_xxxxyyz_0_0_0[i] = g_0_xxxxyy_0_0_0[i] * pb_z + g_0_xxxxyy_0_0_1[i] * wp_z[i];

        g_0_xxxxyzz_0_0_0[i] = g_0_xxxxzz_0_0_0[i] * pb_y + g_0_xxxxzz_0_0_1[i] * wp_y[i];

        g_0_xxxxzzz_0_0_0[i] =
            3.0 * g_0_xxzzz_0_0_0[i] * fi_ab_0 - 3.0 * g_0_xxzzz_0_0_1[i] * fti_ab_0 + g_0_xxxzzz_0_0_0[i] * pb_x + g_0_xxxzzz_0_0_1[i] * wp_x[i];

        g_0_xxxyyyy_0_0_0[i] =
            2.0 * g_0_xyyyy_0_0_0[i] * fi_ab_0 - 2.0 * g_0_xyyyy_0_0_1[i] * fti_ab_0 + g_0_xxyyyy_0_0_0[i] * pb_x + g_0_xxyyyy_0_0_1[i] * wp_x[i];

        g_0_xxxyyyz_0_0_0[i] = g_0_xxxyyy_0_0_0[i] * pb_z + g_0_xxxyyy_0_0_1[i] * wp_z[i];

        g_0_xxxyyzz_0_0_0[i] =
            2.0 * g_0_xyyzz_0_0_0[i] * fi_ab_0 - 2.0 * g_0_xyyzz_0_0_1[i] * fti_ab_0 + g_0_xxyyzz_0_0_0[i] * pb_x + g_0_xxyyzz_0_0_1[i] * wp_x[i];

        g_0_xxxyzzz_0_0_0[i] = g_0_xxxzzz_0_0_0[i] * pb_y + g_0_xxxzzz_0_0_1[i] * wp_y[i];

        g_0_xxxzzzz_0_0_0[i] =
            2.0 * g_0_xzzzz_0_0_0[i] * fi_ab_0 - 2.0 * g_0_xzzzz_0_0_1[i] * fti_ab_0 + g_0_xxzzzz_0_0_0[i] * pb_x + g_0_xxzzzz_0_0_1[i] * wp_x[i];

        g_0_xxyyyyy_0_0_0[i] =
            g_0_yyyyy_0_0_0[i] * fi_ab_0 - g_0_yyyyy_0_0_1[i] * fti_ab_0 + g_0_xyyyyy_0_0_0[i] * pb_x + g_0_xyyyyy_0_0_1[i] * wp_x[i];

        g_0_xxyyyyz_0_0_0[i] = g_0_xxyyyy_0_0_0[i] * pb_z + g_0_xxyyyy_0_0_1[i] * wp_z[i];

        g_0_xxyyyzz_0_0_0[i] =
            g_0_yyyzz_0_0_0[i] * fi_ab_0 - g_0_yyyzz_0_0_1[i] * fti_ab_0 + g_0_xyyyzz_0_0_0[i] * pb_x + g_0_xyyyzz_0_0_1[i] * wp_x[i];

        g_0_xxyyzzz_0_0_0[i] =
            g_0_yyzzz_0_0_0[i] * fi_ab_0 - g_0_yyzzz_0_0_1[i] * fti_ab_0 + g_0_xyyzzz_0_0_0[i] * pb_x + g_0_xyyzzz_0_0_1[i] * wp_x[i];

        g_0_xxyzzzz_0_0_0[i] = g_0_xxzzzz_0_0_0[i] * pb_y + g_0_xxzzzz_0_0_1[i] * wp_y[i];

        g_0_xxzzzzz_0_0_0[i] =
            g_0_zzzzz_0_0_0[i] * fi_ab_0 - g_0_zzzzz_0_0_1[i] * fti_ab_0 + g_0_xzzzzz_0_0_0[i] * pb_x + g_0_xzzzzz_0_0_1[i] * wp_x[i];

        g_0_xyyyyyy_0_0_0[i] = g_0_yyyyyy_0_0_0[i] * pb_x + g_0_yyyyyy_0_0_1[i] * wp_x[i];

        g_0_xyyyyyz_0_0_0[i] = g_0_yyyyyz_0_0_0[i] * pb_x + g_0_yyyyyz_0_0_1[i] * wp_x[i];

        g_0_xyyyyzz_0_0_0[i] = g_0_yyyyzz_0_0_0[i] * pb_x + g_0_yyyyzz_0_0_1[i] * wp_x[i];

        g_0_xyyyzzz_0_0_0[i] = g_0_yyyzzz_0_0_0[i] * pb_x + g_0_yyyzzz_0_0_1[i] * wp_x[i];

        g_0_xyyzzzz_0_0_0[i] = g_0_yyzzzz_0_0_0[i] * pb_x + g_0_yyzzzz_0_0_1[i] * wp_x[i];

        g_0_xyzzzzz_0_0_0[i] = g_0_yzzzzz_0_0_0[i] * pb_x + g_0_yzzzzz_0_0_1[i] * wp_x[i];

        g_0_xzzzzzz_0_0_0[i] = g_0_zzzzzz_0_0_0[i] * pb_x + g_0_zzzzzz_0_0_1[i] * wp_x[i];

        g_0_yyyyyyy_0_0_0[i] =
            6.0 * g_0_yyyyy_0_0_0[i] * fi_ab_0 - 6.0 * g_0_yyyyy_0_0_1[i] * fti_ab_0 + g_0_yyyyyy_0_0_0[i] * pb_y + g_0_yyyyyy_0_0_1[i] * wp_y[i];

        g_0_yyyyyyz_0_0_0[i] = g_0_yyyyyy_0_0_0[i] * pb_z + g_0_yyyyyy_0_0_1[i] * wp_z[i];

        g_0_yyyyyzz_0_0_0[i] =
            4.0 * g_0_yyyzz_0_0_0[i] * fi_ab_0 - 4.0 * g_0_yyyzz_0_0_1[i] * fti_ab_0 + g_0_yyyyzz_0_0_0[i] * pb_y + g_0_yyyyzz_0_0_1[i] * wp_y[i];

        g_0_yyyyzzz_0_0_0[i] =
            3.0 * g_0_yyzzz_0_0_0[i] * fi_ab_0 - 3.0 * g_0_yyzzz_0_0_1[i] * fti_ab_0 + g_0_yyyzzz_0_0_0[i] * pb_y + g_0_yyyzzz_0_0_1[i] * wp_y[i];

        g_0_yyyzzzz_0_0_0[i] =
            2.0 * g_0_yzzzz_0_0_0[i] * fi_ab_0 - 2.0 * g_0_yzzzz_0_0_1[i] * fti_ab_0 + g_0_yyzzzz_0_0_0[i] * pb_y + g_0_yyzzzz_0_0_1[i] * wp_y[i];

        g_0_yyzzzzz_0_0_0[i] =
            g_0_zzzzz_0_0_0[i] * fi_ab_0 - g_0_zzzzz_0_0_1[i] * fti_ab_0 + g_0_yzzzzz_0_0_0[i] * pb_y + g_0_yzzzzz_0_0_1[i] * wp_y[i];

        g_0_yzzzzzz_0_0_0[i] = g_0_zzzzzz_0_0_0[i] * pb_y + g_0_zzzzzz_0_0_1[i] * wp_y[i];

        g_0_zzzzzzz_0_0_0[i] =
            6.0 * g_0_zzzzz_0_0_0[i] * fi_ab_0 - 6.0 * g_0_zzzzz_0_0_1[i] * fti_ab_0 + g_0_zzzzzz_0_0_0[i] * pb_z + g_0_zzzzzz_0_0_1[i] * wp_z[i];
    }
}

}  // namespace erirec
