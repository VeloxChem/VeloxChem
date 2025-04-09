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

#include "ElectronRepulsionPrimRecSSSK.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sssk(CSimdArray<double>& pbuffer,
                                  const size_t        idx_eri_0_sssk,
                                  size_t              idx_eri_0_sssh,
                                  size_t              idx_eri_1_sssh,
                                  size_t              idx_eri_0_sssi,
                                  size_t              idx_eri_1_sssi,
                                  CSimdArray<double>& factors,
                                  const size_t        idx_qd,
                                  const size_t        idx_wq,
                                  const double        a_exp,
                                  const double        b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(QD) distances

    auto qd_x = factors.data(idx_qd);

    auto qd_y = factors.data(idx_qd + 1);

    auto qd_z = factors.data(idx_qd + 2);

    // Set up R(WQ) distances

    auto wq_x = factors.data(idx_wq);

    auto wq_y = factors.data(idx_wq + 1);

    auto wq_z = factors.data(idx_wq + 2);

    /// Set up components of auxilary buffer : SSSH

    auto g_0_0_0_xxxxx_0 = pbuffer.data(idx_eri_0_sssh);

    auto g_0_0_0_xxxyy_0 = pbuffer.data(idx_eri_0_sssh + 3);

    auto g_0_0_0_xxxzz_0 = pbuffer.data(idx_eri_0_sssh + 5);

    auto g_0_0_0_xxyyy_0 = pbuffer.data(idx_eri_0_sssh + 6);

    auto g_0_0_0_xxzzz_0 = pbuffer.data(idx_eri_0_sssh + 9);

    auto g_0_0_0_xyyyy_0 = pbuffer.data(idx_eri_0_sssh + 10);

    auto g_0_0_0_xyyzz_0 = pbuffer.data(idx_eri_0_sssh + 12);

    auto g_0_0_0_xzzzz_0 = pbuffer.data(idx_eri_0_sssh + 14);

    auto g_0_0_0_yyyyy_0 = pbuffer.data(idx_eri_0_sssh + 15);

    auto g_0_0_0_yyyzz_0 = pbuffer.data(idx_eri_0_sssh + 17);

    auto g_0_0_0_yyzzz_0 = pbuffer.data(idx_eri_0_sssh + 18);

    auto g_0_0_0_yzzzz_0 = pbuffer.data(idx_eri_0_sssh + 19);

    auto g_0_0_0_zzzzz_0 = pbuffer.data(idx_eri_0_sssh + 20);

    /// Set up components of auxilary buffer : SSSH

    auto g_0_0_0_xxxxx_1 = pbuffer.data(idx_eri_1_sssh);

    auto g_0_0_0_xxxyy_1 = pbuffer.data(idx_eri_1_sssh + 3);

    auto g_0_0_0_xxxzz_1 = pbuffer.data(idx_eri_1_sssh + 5);

    auto g_0_0_0_xxyyy_1 = pbuffer.data(idx_eri_1_sssh + 6);

    auto g_0_0_0_xxzzz_1 = pbuffer.data(idx_eri_1_sssh + 9);

    auto g_0_0_0_xyyyy_1 = pbuffer.data(idx_eri_1_sssh + 10);

    auto g_0_0_0_xyyzz_1 = pbuffer.data(idx_eri_1_sssh + 12);

    auto g_0_0_0_xzzzz_1 = pbuffer.data(idx_eri_1_sssh + 14);

    auto g_0_0_0_yyyyy_1 = pbuffer.data(idx_eri_1_sssh + 15);

    auto g_0_0_0_yyyzz_1 = pbuffer.data(idx_eri_1_sssh + 17);

    auto g_0_0_0_yyzzz_1 = pbuffer.data(idx_eri_1_sssh + 18);

    auto g_0_0_0_yzzzz_1 = pbuffer.data(idx_eri_1_sssh + 19);

    auto g_0_0_0_zzzzz_1 = pbuffer.data(idx_eri_1_sssh + 20);

    /// Set up components of auxilary buffer : SSSI

    auto g_0_0_0_xxxxxx_0 = pbuffer.data(idx_eri_0_sssi);

    auto g_0_0_0_xxxxxz_0 = pbuffer.data(idx_eri_0_sssi + 2);

    auto g_0_0_0_xxxxyy_0 = pbuffer.data(idx_eri_0_sssi + 3);

    auto g_0_0_0_xxxxzz_0 = pbuffer.data(idx_eri_0_sssi + 5);

    auto g_0_0_0_xxxyyy_0 = pbuffer.data(idx_eri_0_sssi + 6);

    auto g_0_0_0_xxxzzz_0 = pbuffer.data(idx_eri_0_sssi + 9);

    auto g_0_0_0_xxyyyy_0 = pbuffer.data(idx_eri_0_sssi + 10);

    auto g_0_0_0_xxyyzz_0 = pbuffer.data(idx_eri_0_sssi + 12);

    auto g_0_0_0_xxzzzz_0 = pbuffer.data(idx_eri_0_sssi + 14);

    auto g_0_0_0_xyyyyy_0 = pbuffer.data(idx_eri_0_sssi + 15);

    auto g_0_0_0_xyyyzz_0 = pbuffer.data(idx_eri_0_sssi + 17);

    auto g_0_0_0_xyyzzz_0 = pbuffer.data(idx_eri_0_sssi + 18);

    auto g_0_0_0_xzzzzz_0 = pbuffer.data(idx_eri_0_sssi + 20);

    auto g_0_0_0_yyyyyy_0 = pbuffer.data(idx_eri_0_sssi + 21);

    auto g_0_0_0_yyyyyz_0 = pbuffer.data(idx_eri_0_sssi + 22);

    auto g_0_0_0_yyyyzz_0 = pbuffer.data(idx_eri_0_sssi + 23);

    auto g_0_0_0_yyyzzz_0 = pbuffer.data(idx_eri_0_sssi + 24);

    auto g_0_0_0_yyzzzz_0 = pbuffer.data(idx_eri_0_sssi + 25);

    auto g_0_0_0_yzzzzz_0 = pbuffer.data(idx_eri_0_sssi + 26);

    auto g_0_0_0_zzzzzz_0 = pbuffer.data(idx_eri_0_sssi + 27);

    /// Set up components of auxilary buffer : SSSI

    auto g_0_0_0_xxxxxx_1 = pbuffer.data(idx_eri_1_sssi);

    auto g_0_0_0_xxxxxz_1 = pbuffer.data(idx_eri_1_sssi + 2);

    auto g_0_0_0_xxxxyy_1 = pbuffer.data(idx_eri_1_sssi + 3);

    auto g_0_0_0_xxxxzz_1 = pbuffer.data(idx_eri_1_sssi + 5);

    auto g_0_0_0_xxxyyy_1 = pbuffer.data(idx_eri_1_sssi + 6);

    auto g_0_0_0_xxxzzz_1 = pbuffer.data(idx_eri_1_sssi + 9);

    auto g_0_0_0_xxyyyy_1 = pbuffer.data(idx_eri_1_sssi + 10);

    auto g_0_0_0_xxyyzz_1 = pbuffer.data(idx_eri_1_sssi + 12);

    auto g_0_0_0_xxzzzz_1 = pbuffer.data(idx_eri_1_sssi + 14);

    auto g_0_0_0_xyyyyy_1 = pbuffer.data(idx_eri_1_sssi + 15);

    auto g_0_0_0_xyyyzz_1 = pbuffer.data(idx_eri_1_sssi + 17);

    auto g_0_0_0_xyyzzz_1 = pbuffer.data(idx_eri_1_sssi + 18);

    auto g_0_0_0_xzzzzz_1 = pbuffer.data(idx_eri_1_sssi + 20);

    auto g_0_0_0_yyyyyy_1 = pbuffer.data(idx_eri_1_sssi + 21);

    auto g_0_0_0_yyyyyz_1 = pbuffer.data(idx_eri_1_sssi + 22);

    auto g_0_0_0_yyyyzz_1 = pbuffer.data(idx_eri_1_sssi + 23);

    auto g_0_0_0_yyyzzz_1 = pbuffer.data(idx_eri_1_sssi + 24);

    auto g_0_0_0_yyzzzz_1 = pbuffer.data(idx_eri_1_sssi + 25);

    auto g_0_0_0_yzzzzz_1 = pbuffer.data(idx_eri_1_sssi + 26);

    auto g_0_0_0_zzzzzz_1 = pbuffer.data(idx_eri_1_sssi + 27);

    /// Set up components of targeted buffer : SSSK

    auto g_0_0_0_xxxxxxx_0 = pbuffer.data(idx_eri_0_sssk);

    auto g_0_0_0_xxxxxxy_0 = pbuffer.data(idx_eri_0_sssk + 1);

    auto g_0_0_0_xxxxxxz_0 = pbuffer.data(idx_eri_0_sssk + 2);

    auto g_0_0_0_xxxxxyy_0 = pbuffer.data(idx_eri_0_sssk + 3);

    auto g_0_0_0_xxxxxyz_0 = pbuffer.data(idx_eri_0_sssk + 4);

    auto g_0_0_0_xxxxxzz_0 = pbuffer.data(idx_eri_0_sssk + 5);

    auto g_0_0_0_xxxxyyy_0 = pbuffer.data(idx_eri_0_sssk + 6);

    auto g_0_0_0_xxxxyyz_0 = pbuffer.data(idx_eri_0_sssk + 7);

    auto g_0_0_0_xxxxyzz_0 = pbuffer.data(idx_eri_0_sssk + 8);

    auto g_0_0_0_xxxxzzz_0 = pbuffer.data(idx_eri_0_sssk + 9);

    auto g_0_0_0_xxxyyyy_0 = pbuffer.data(idx_eri_0_sssk + 10);

    auto g_0_0_0_xxxyyyz_0 = pbuffer.data(idx_eri_0_sssk + 11);

    auto g_0_0_0_xxxyyzz_0 = pbuffer.data(idx_eri_0_sssk + 12);

    auto g_0_0_0_xxxyzzz_0 = pbuffer.data(idx_eri_0_sssk + 13);

    auto g_0_0_0_xxxzzzz_0 = pbuffer.data(idx_eri_0_sssk + 14);

    auto g_0_0_0_xxyyyyy_0 = pbuffer.data(idx_eri_0_sssk + 15);

    auto g_0_0_0_xxyyyyz_0 = pbuffer.data(idx_eri_0_sssk + 16);

    auto g_0_0_0_xxyyyzz_0 = pbuffer.data(idx_eri_0_sssk + 17);

    auto g_0_0_0_xxyyzzz_0 = pbuffer.data(idx_eri_0_sssk + 18);

    auto g_0_0_0_xxyzzzz_0 = pbuffer.data(idx_eri_0_sssk + 19);

    auto g_0_0_0_xxzzzzz_0 = pbuffer.data(idx_eri_0_sssk + 20);

    auto g_0_0_0_xyyyyyy_0 = pbuffer.data(idx_eri_0_sssk + 21);

    auto g_0_0_0_xyyyyyz_0 = pbuffer.data(idx_eri_0_sssk + 22);

    auto g_0_0_0_xyyyyzz_0 = pbuffer.data(idx_eri_0_sssk + 23);

    auto g_0_0_0_xyyyzzz_0 = pbuffer.data(idx_eri_0_sssk + 24);

    auto g_0_0_0_xyyzzzz_0 = pbuffer.data(idx_eri_0_sssk + 25);

    auto g_0_0_0_xyzzzzz_0 = pbuffer.data(idx_eri_0_sssk + 26);

    auto g_0_0_0_xzzzzzz_0 = pbuffer.data(idx_eri_0_sssk + 27);

    auto g_0_0_0_yyyyyyy_0 = pbuffer.data(idx_eri_0_sssk + 28);

    auto g_0_0_0_yyyyyyz_0 = pbuffer.data(idx_eri_0_sssk + 29);

    auto g_0_0_0_yyyyyzz_0 = pbuffer.data(idx_eri_0_sssk + 30);

    auto g_0_0_0_yyyyzzz_0 = pbuffer.data(idx_eri_0_sssk + 31);

    auto g_0_0_0_yyyzzzz_0 = pbuffer.data(idx_eri_0_sssk + 32);

    auto g_0_0_0_yyzzzzz_0 = pbuffer.data(idx_eri_0_sssk + 33);

    auto g_0_0_0_yzzzzzz_0 = pbuffer.data(idx_eri_0_sssk + 34);

    auto g_0_0_0_zzzzzzz_0 = pbuffer.data(idx_eri_0_sssk + 35);

#pragma omp simd aligned(g_0_0_0_xxxxx_0,       \
                             g_0_0_0_xxxxx_1,   \
                             g_0_0_0_xxxxxx_0,  \
                             g_0_0_0_xxxxxx_1,  \
                             g_0_0_0_xxxxxxx_0, \
                             g_0_0_0_xxxxxxy_0, \
                             g_0_0_0_xxxxxxz_0, \
                             g_0_0_0_xxxxxyy_0, \
                             g_0_0_0_xxxxxyz_0, \
                             g_0_0_0_xxxxxz_0,  \
                             g_0_0_0_xxxxxz_1,  \
                             g_0_0_0_xxxxxzz_0, \
                             g_0_0_0_xxxxyy_0,  \
                             g_0_0_0_xxxxyy_1,  \
                             g_0_0_0_xxxxyyy_0, \
                             g_0_0_0_xxxxyyz_0, \
                             g_0_0_0_xxxxyzz_0, \
                             g_0_0_0_xxxxzz_0,  \
                             g_0_0_0_xxxxzz_1,  \
                             g_0_0_0_xxxxzzz_0, \
                             g_0_0_0_xxxyy_0,   \
                             g_0_0_0_xxxyy_1,   \
                             g_0_0_0_xxxyyy_0,  \
                             g_0_0_0_xxxyyy_1,  \
                             g_0_0_0_xxxyyyy_0, \
                             g_0_0_0_xxxyyyz_0, \
                             g_0_0_0_xxxyyzz_0, \
                             g_0_0_0_xxxyzzz_0, \
                             g_0_0_0_xxxzz_0,   \
                             g_0_0_0_xxxzz_1,   \
                             g_0_0_0_xxxzzz_0,  \
                             g_0_0_0_xxxzzz_1,  \
                             g_0_0_0_xxxzzzz_0, \
                             g_0_0_0_xxyyy_0,   \
                             g_0_0_0_xxyyy_1,   \
                             g_0_0_0_xxyyyy_0,  \
                             g_0_0_0_xxyyyy_1,  \
                             g_0_0_0_xxyyyyy_0, \
                             g_0_0_0_xxyyyyz_0, \
                             g_0_0_0_xxyyyzz_0, \
                             g_0_0_0_xxyyzz_0,  \
                             g_0_0_0_xxyyzz_1,  \
                             g_0_0_0_xxyyzzz_0, \
                             g_0_0_0_xxyzzzz_0, \
                             g_0_0_0_xxzzz_0,   \
                             g_0_0_0_xxzzz_1,   \
                             g_0_0_0_xxzzzz_0,  \
                             g_0_0_0_xxzzzz_1,  \
                             g_0_0_0_xxzzzzz_0, \
                             g_0_0_0_xyyyy_0,   \
                             g_0_0_0_xyyyy_1,   \
                             g_0_0_0_xyyyyy_0,  \
                             g_0_0_0_xyyyyy_1,  \
                             g_0_0_0_xyyyyyy_0, \
                             g_0_0_0_xyyyyyz_0, \
                             g_0_0_0_xyyyyzz_0, \
                             g_0_0_0_xyyyzz_0,  \
                             g_0_0_0_xyyyzz_1,  \
                             g_0_0_0_xyyyzzz_0, \
                             g_0_0_0_xyyzz_0,   \
                             g_0_0_0_xyyzz_1,   \
                             g_0_0_0_xyyzzz_0,  \
                             g_0_0_0_xyyzzz_1,  \
                             g_0_0_0_xyyzzzz_0, \
                             g_0_0_0_xyzzzzz_0, \
                             g_0_0_0_xzzzz_0,   \
                             g_0_0_0_xzzzz_1,   \
                             g_0_0_0_xzzzzz_0,  \
                             g_0_0_0_xzzzzz_1,  \
                             g_0_0_0_xzzzzzz_0, \
                             g_0_0_0_yyyyy_0,   \
                             g_0_0_0_yyyyy_1,   \
                             g_0_0_0_yyyyyy_0,  \
                             g_0_0_0_yyyyyy_1,  \
                             g_0_0_0_yyyyyyy_0, \
                             g_0_0_0_yyyyyyz_0, \
                             g_0_0_0_yyyyyz_0,  \
                             g_0_0_0_yyyyyz_1,  \
                             g_0_0_0_yyyyyzz_0, \
                             g_0_0_0_yyyyzz_0,  \
                             g_0_0_0_yyyyzz_1,  \
                             g_0_0_0_yyyyzzz_0, \
                             g_0_0_0_yyyzz_0,   \
                             g_0_0_0_yyyzz_1,   \
                             g_0_0_0_yyyzzz_0,  \
                             g_0_0_0_yyyzzz_1,  \
                             g_0_0_0_yyyzzzz_0, \
                             g_0_0_0_yyzzz_0,   \
                             g_0_0_0_yyzzz_1,   \
                             g_0_0_0_yyzzzz_0,  \
                             g_0_0_0_yyzzzz_1,  \
                             g_0_0_0_yyzzzzz_0, \
                             g_0_0_0_yzzzz_0,   \
                             g_0_0_0_yzzzz_1,   \
                             g_0_0_0_yzzzzz_0,  \
                             g_0_0_0_yzzzzz_1,  \
                             g_0_0_0_yzzzzzz_0, \
                             g_0_0_0_zzzzz_0,   \
                             g_0_0_0_zzzzz_1,   \
                             g_0_0_0_zzzzzz_0,  \
                             g_0_0_0_zzzzzz_1,  \
                             g_0_0_0_zzzzzzz_0, \
                             qd_x,              \
                             qd_y,              \
                             qd_z,              \
                             wq_x,              \
                             wq_y,              \
                             wq_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_cd_0 = 0.5 / (c_exps[i] + d_exps[i]);

        const double fti_cd_0 = fi_cd_0 * (a_exp + b_exp) / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_0_0_xxxxxxx_0[i] =
            6.0 * g_0_0_0_xxxxx_0[i] * fi_cd_0 - 6.0 * g_0_0_0_xxxxx_1[i] * fti_cd_0 + g_0_0_0_xxxxxx_0[i] * qd_x[i] + g_0_0_0_xxxxxx_1[i] * wq_x[i];

        g_0_0_0_xxxxxxy_0[i] = g_0_0_0_xxxxxx_0[i] * qd_y[i] + g_0_0_0_xxxxxx_1[i] * wq_y[i];

        g_0_0_0_xxxxxxz_0[i] = g_0_0_0_xxxxxx_0[i] * qd_z[i] + g_0_0_0_xxxxxx_1[i] * wq_z[i];

        g_0_0_0_xxxxxyy_0[i] =
            4.0 * g_0_0_0_xxxyy_0[i] * fi_cd_0 - 4.0 * g_0_0_0_xxxyy_1[i] * fti_cd_0 + g_0_0_0_xxxxyy_0[i] * qd_x[i] + g_0_0_0_xxxxyy_1[i] * wq_x[i];

        g_0_0_0_xxxxxyz_0[i] = g_0_0_0_xxxxxz_0[i] * qd_y[i] + g_0_0_0_xxxxxz_1[i] * wq_y[i];

        g_0_0_0_xxxxxzz_0[i] =
            4.0 * g_0_0_0_xxxzz_0[i] * fi_cd_0 - 4.0 * g_0_0_0_xxxzz_1[i] * fti_cd_0 + g_0_0_0_xxxxzz_0[i] * qd_x[i] + g_0_0_0_xxxxzz_1[i] * wq_x[i];

        g_0_0_0_xxxxyyy_0[i] =
            3.0 * g_0_0_0_xxyyy_0[i] * fi_cd_0 - 3.0 * g_0_0_0_xxyyy_1[i] * fti_cd_0 + g_0_0_0_xxxyyy_0[i] * qd_x[i] + g_0_0_0_xxxyyy_1[i] * wq_x[i];

        g_0_0_0_xxxxyyz_0[i] = g_0_0_0_xxxxyy_0[i] * qd_z[i] + g_0_0_0_xxxxyy_1[i] * wq_z[i];

        g_0_0_0_xxxxyzz_0[i] = g_0_0_0_xxxxzz_0[i] * qd_y[i] + g_0_0_0_xxxxzz_1[i] * wq_y[i];

        g_0_0_0_xxxxzzz_0[i] =
            3.0 * g_0_0_0_xxzzz_0[i] * fi_cd_0 - 3.0 * g_0_0_0_xxzzz_1[i] * fti_cd_0 + g_0_0_0_xxxzzz_0[i] * qd_x[i] + g_0_0_0_xxxzzz_1[i] * wq_x[i];

        g_0_0_0_xxxyyyy_0[i] =
            2.0 * g_0_0_0_xyyyy_0[i] * fi_cd_0 - 2.0 * g_0_0_0_xyyyy_1[i] * fti_cd_0 + g_0_0_0_xxyyyy_0[i] * qd_x[i] + g_0_0_0_xxyyyy_1[i] * wq_x[i];

        g_0_0_0_xxxyyyz_0[i] = g_0_0_0_xxxyyy_0[i] * qd_z[i] + g_0_0_0_xxxyyy_1[i] * wq_z[i];

        g_0_0_0_xxxyyzz_0[i] =
            2.0 * g_0_0_0_xyyzz_0[i] * fi_cd_0 - 2.0 * g_0_0_0_xyyzz_1[i] * fti_cd_0 + g_0_0_0_xxyyzz_0[i] * qd_x[i] + g_0_0_0_xxyyzz_1[i] * wq_x[i];

        g_0_0_0_xxxyzzz_0[i] = g_0_0_0_xxxzzz_0[i] * qd_y[i] + g_0_0_0_xxxzzz_1[i] * wq_y[i];

        g_0_0_0_xxxzzzz_0[i] =
            2.0 * g_0_0_0_xzzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_0_xzzzz_1[i] * fti_cd_0 + g_0_0_0_xxzzzz_0[i] * qd_x[i] + g_0_0_0_xxzzzz_1[i] * wq_x[i];

        g_0_0_0_xxyyyyy_0[i] =
            g_0_0_0_yyyyy_0[i] * fi_cd_0 - g_0_0_0_yyyyy_1[i] * fti_cd_0 + g_0_0_0_xyyyyy_0[i] * qd_x[i] + g_0_0_0_xyyyyy_1[i] * wq_x[i];

        g_0_0_0_xxyyyyz_0[i] = g_0_0_0_xxyyyy_0[i] * qd_z[i] + g_0_0_0_xxyyyy_1[i] * wq_z[i];

        g_0_0_0_xxyyyzz_0[i] =
            g_0_0_0_yyyzz_0[i] * fi_cd_0 - g_0_0_0_yyyzz_1[i] * fti_cd_0 + g_0_0_0_xyyyzz_0[i] * qd_x[i] + g_0_0_0_xyyyzz_1[i] * wq_x[i];

        g_0_0_0_xxyyzzz_0[i] =
            g_0_0_0_yyzzz_0[i] * fi_cd_0 - g_0_0_0_yyzzz_1[i] * fti_cd_0 + g_0_0_0_xyyzzz_0[i] * qd_x[i] + g_0_0_0_xyyzzz_1[i] * wq_x[i];

        g_0_0_0_xxyzzzz_0[i] = g_0_0_0_xxzzzz_0[i] * qd_y[i] + g_0_0_0_xxzzzz_1[i] * wq_y[i];

        g_0_0_0_xxzzzzz_0[i] =
            g_0_0_0_zzzzz_0[i] * fi_cd_0 - g_0_0_0_zzzzz_1[i] * fti_cd_0 + g_0_0_0_xzzzzz_0[i] * qd_x[i] + g_0_0_0_xzzzzz_1[i] * wq_x[i];

        g_0_0_0_xyyyyyy_0[i] = g_0_0_0_yyyyyy_0[i] * qd_x[i] + g_0_0_0_yyyyyy_1[i] * wq_x[i];

        g_0_0_0_xyyyyyz_0[i] = g_0_0_0_yyyyyz_0[i] * qd_x[i] + g_0_0_0_yyyyyz_1[i] * wq_x[i];

        g_0_0_0_xyyyyzz_0[i] = g_0_0_0_yyyyzz_0[i] * qd_x[i] + g_0_0_0_yyyyzz_1[i] * wq_x[i];

        g_0_0_0_xyyyzzz_0[i] = g_0_0_0_yyyzzz_0[i] * qd_x[i] + g_0_0_0_yyyzzz_1[i] * wq_x[i];

        g_0_0_0_xyyzzzz_0[i] = g_0_0_0_yyzzzz_0[i] * qd_x[i] + g_0_0_0_yyzzzz_1[i] * wq_x[i];

        g_0_0_0_xyzzzzz_0[i] = g_0_0_0_yzzzzz_0[i] * qd_x[i] + g_0_0_0_yzzzzz_1[i] * wq_x[i];

        g_0_0_0_xzzzzzz_0[i] = g_0_0_0_zzzzzz_0[i] * qd_x[i] + g_0_0_0_zzzzzz_1[i] * wq_x[i];

        g_0_0_0_yyyyyyy_0[i] =
            6.0 * g_0_0_0_yyyyy_0[i] * fi_cd_0 - 6.0 * g_0_0_0_yyyyy_1[i] * fti_cd_0 + g_0_0_0_yyyyyy_0[i] * qd_y[i] + g_0_0_0_yyyyyy_1[i] * wq_y[i];

        g_0_0_0_yyyyyyz_0[i] = g_0_0_0_yyyyyy_0[i] * qd_z[i] + g_0_0_0_yyyyyy_1[i] * wq_z[i];

        g_0_0_0_yyyyyzz_0[i] =
            4.0 * g_0_0_0_yyyzz_0[i] * fi_cd_0 - 4.0 * g_0_0_0_yyyzz_1[i] * fti_cd_0 + g_0_0_0_yyyyzz_0[i] * qd_y[i] + g_0_0_0_yyyyzz_1[i] * wq_y[i];

        g_0_0_0_yyyyzzz_0[i] =
            3.0 * g_0_0_0_yyzzz_0[i] * fi_cd_0 - 3.0 * g_0_0_0_yyzzz_1[i] * fti_cd_0 + g_0_0_0_yyyzzz_0[i] * qd_y[i] + g_0_0_0_yyyzzz_1[i] * wq_y[i];

        g_0_0_0_yyyzzzz_0[i] =
            2.0 * g_0_0_0_yzzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_0_yzzzz_1[i] * fti_cd_0 + g_0_0_0_yyzzzz_0[i] * qd_y[i] + g_0_0_0_yyzzzz_1[i] * wq_y[i];

        g_0_0_0_yyzzzzz_0[i] =
            g_0_0_0_zzzzz_0[i] * fi_cd_0 - g_0_0_0_zzzzz_1[i] * fti_cd_0 + g_0_0_0_yzzzzz_0[i] * qd_y[i] + g_0_0_0_yzzzzz_1[i] * wq_y[i];

        g_0_0_0_yzzzzzz_0[i] = g_0_0_0_zzzzzz_0[i] * qd_y[i] + g_0_0_0_zzzzzz_1[i] * wq_y[i];

        g_0_0_0_zzzzzzz_0[i] =
            6.0 * g_0_0_0_zzzzz_0[i] * fi_cd_0 - 6.0 * g_0_0_0_zzzzz_1[i] * fti_cd_0 + g_0_0_0_zzzzzz_0[i] * qd_z[i] + g_0_0_0_zzzzzz_1[i] * wq_z[i];
    }
}

}  // namespace erirec
