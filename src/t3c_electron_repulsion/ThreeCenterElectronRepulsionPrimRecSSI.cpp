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

#include "ThreeCenterElectronRepulsionPrimRecSSI.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_ssi(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_ssi,
                                 size_t idx_eri_0_ssg,
                                 size_t idx_eri_1_ssg,
                                 size_t idx_eri_0_ssh,
                                 size_t idx_eri_1_ssh,
                                 CSimdArray<double>& factors,
                                 const size_t idx_qd,
                                 const size_t idx_wq,
                                 const double a_exp) -> void
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

    /// Set up components of auxilary buffer : SSG

    auto g_0_0_xxxx_0 = pbuffer.data(idx_eri_0_ssg);

    auto g_0_0_xxyy_0 = pbuffer.data(idx_eri_0_ssg + 3);

    auto g_0_0_xxzz_0 = pbuffer.data(idx_eri_0_ssg + 5);

    auto g_0_0_xyyy_0 = pbuffer.data(idx_eri_0_ssg + 6);

    auto g_0_0_xzzz_0 = pbuffer.data(idx_eri_0_ssg + 9);

    auto g_0_0_yyyy_0 = pbuffer.data(idx_eri_0_ssg + 10);

    auto g_0_0_yyzz_0 = pbuffer.data(idx_eri_0_ssg + 12);

    auto g_0_0_yzzz_0 = pbuffer.data(idx_eri_0_ssg + 13);

    auto g_0_0_zzzz_0 = pbuffer.data(idx_eri_0_ssg + 14);

    /// Set up components of auxilary buffer : SSG

    auto g_0_0_xxxx_1 = pbuffer.data(idx_eri_1_ssg);

    auto g_0_0_xxyy_1 = pbuffer.data(idx_eri_1_ssg + 3);

    auto g_0_0_xxzz_1 = pbuffer.data(idx_eri_1_ssg + 5);

    auto g_0_0_xyyy_1 = pbuffer.data(idx_eri_1_ssg + 6);

    auto g_0_0_xzzz_1 = pbuffer.data(idx_eri_1_ssg + 9);

    auto g_0_0_yyyy_1 = pbuffer.data(idx_eri_1_ssg + 10);

    auto g_0_0_yyzz_1 = pbuffer.data(idx_eri_1_ssg + 12);

    auto g_0_0_yzzz_1 = pbuffer.data(idx_eri_1_ssg + 13);

    auto g_0_0_zzzz_1 = pbuffer.data(idx_eri_1_ssg + 14);

    /// Set up components of auxilary buffer : SSH

    auto g_0_0_xxxxx_0 = pbuffer.data(idx_eri_0_ssh);

    auto g_0_0_xxxxz_0 = pbuffer.data(idx_eri_0_ssh + 2);

    auto g_0_0_xxxyy_0 = pbuffer.data(idx_eri_0_ssh + 3);

    auto g_0_0_xxxzz_0 = pbuffer.data(idx_eri_0_ssh + 5);

    auto g_0_0_xxyyy_0 = pbuffer.data(idx_eri_0_ssh + 6);

    auto g_0_0_xxzzz_0 = pbuffer.data(idx_eri_0_ssh + 9);

    auto g_0_0_xyyyy_0 = pbuffer.data(idx_eri_0_ssh + 10);

    auto g_0_0_xyyzz_0 = pbuffer.data(idx_eri_0_ssh + 12);

    auto g_0_0_xzzzz_0 = pbuffer.data(idx_eri_0_ssh + 14);

    auto g_0_0_yyyyy_0 = pbuffer.data(idx_eri_0_ssh + 15);

    auto g_0_0_yyyyz_0 = pbuffer.data(idx_eri_0_ssh + 16);

    auto g_0_0_yyyzz_0 = pbuffer.data(idx_eri_0_ssh + 17);

    auto g_0_0_yyzzz_0 = pbuffer.data(idx_eri_0_ssh + 18);

    auto g_0_0_yzzzz_0 = pbuffer.data(idx_eri_0_ssh + 19);

    auto g_0_0_zzzzz_0 = pbuffer.data(idx_eri_0_ssh + 20);

    /// Set up components of auxilary buffer : SSH

    auto g_0_0_xxxxx_1 = pbuffer.data(idx_eri_1_ssh);

    auto g_0_0_xxxxz_1 = pbuffer.data(idx_eri_1_ssh + 2);

    auto g_0_0_xxxyy_1 = pbuffer.data(idx_eri_1_ssh + 3);

    auto g_0_0_xxxzz_1 = pbuffer.data(idx_eri_1_ssh + 5);

    auto g_0_0_xxyyy_1 = pbuffer.data(idx_eri_1_ssh + 6);

    auto g_0_0_xxzzz_1 = pbuffer.data(idx_eri_1_ssh + 9);

    auto g_0_0_xyyyy_1 = pbuffer.data(idx_eri_1_ssh + 10);

    auto g_0_0_xyyzz_1 = pbuffer.data(idx_eri_1_ssh + 12);

    auto g_0_0_xzzzz_1 = pbuffer.data(idx_eri_1_ssh + 14);

    auto g_0_0_yyyyy_1 = pbuffer.data(idx_eri_1_ssh + 15);

    auto g_0_0_yyyyz_1 = pbuffer.data(idx_eri_1_ssh + 16);

    auto g_0_0_yyyzz_1 = pbuffer.data(idx_eri_1_ssh + 17);

    auto g_0_0_yyzzz_1 = pbuffer.data(idx_eri_1_ssh + 18);

    auto g_0_0_yzzzz_1 = pbuffer.data(idx_eri_1_ssh + 19);

    auto g_0_0_zzzzz_1 = pbuffer.data(idx_eri_1_ssh + 20);

    /// Set up components of targeted buffer : SSI

    auto g_0_0_xxxxxx_0 = pbuffer.data(idx_eri_0_ssi);

    auto g_0_0_xxxxxy_0 = pbuffer.data(idx_eri_0_ssi + 1);

    auto g_0_0_xxxxxz_0 = pbuffer.data(idx_eri_0_ssi + 2);

    auto g_0_0_xxxxyy_0 = pbuffer.data(idx_eri_0_ssi + 3);

    auto g_0_0_xxxxyz_0 = pbuffer.data(idx_eri_0_ssi + 4);

    auto g_0_0_xxxxzz_0 = pbuffer.data(idx_eri_0_ssi + 5);

    auto g_0_0_xxxyyy_0 = pbuffer.data(idx_eri_0_ssi + 6);

    auto g_0_0_xxxyyz_0 = pbuffer.data(idx_eri_0_ssi + 7);

    auto g_0_0_xxxyzz_0 = pbuffer.data(idx_eri_0_ssi + 8);

    auto g_0_0_xxxzzz_0 = pbuffer.data(idx_eri_0_ssi + 9);

    auto g_0_0_xxyyyy_0 = pbuffer.data(idx_eri_0_ssi + 10);

    auto g_0_0_xxyyyz_0 = pbuffer.data(idx_eri_0_ssi + 11);

    auto g_0_0_xxyyzz_0 = pbuffer.data(idx_eri_0_ssi + 12);

    auto g_0_0_xxyzzz_0 = pbuffer.data(idx_eri_0_ssi + 13);

    auto g_0_0_xxzzzz_0 = pbuffer.data(idx_eri_0_ssi + 14);

    auto g_0_0_xyyyyy_0 = pbuffer.data(idx_eri_0_ssi + 15);

    auto g_0_0_xyyyyz_0 = pbuffer.data(idx_eri_0_ssi + 16);

    auto g_0_0_xyyyzz_0 = pbuffer.data(idx_eri_0_ssi + 17);

    auto g_0_0_xyyzzz_0 = pbuffer.data(idx_eri_0_ssi + 18);

    auto g_0_0_xyzzzz_0 = pbuffer.data(idx_eri_0_ssi + 19);

    auto g_0_0_xzzzzz_0 = pbuffer.data(idx_eri_0_ssi + 20);

    auto g_0_0_yyyyyy_0 = pbuffer.data(idx_eri_0_ssi + 21);

    auto g_0_0_yyyyyz_0 = pbuffer.data(idx_eri_0_ssi + 22);

    auto g_0_0_yyyyzz_0 = pbuffer.data(idx_eri_0_ssi + 23);

    auto g_0_0_yyyzzz_0 = pbuffer.data(idx_eri_0_ssi + 24);

    auto g_0_0_yyzzzz_0 = pbuffer.data(idx_eri_0_ssi + 25);

    auto g_0_0_yzzzzz_0 = pbuffer.data(idx_eri_0_ssi + 26);

    auto g_0_0_zzzzzz_0 = pbuffer.data(idx_eri_0_ssi + 27);

    #pragma omp simd aligned(g_0_0_xxxx_0, g_0_0_xxxx_1, g_0_0_xxxxx_0, g_0_0_xxxxx_1, g_0_0_xxxxxx_0, g_0_0_xxxxxy_0, g_0_0_xxxxxz_0, g_0_0_xxxxyy_0, g_0_0_xxxxyz_0, g_0_0_xxxxz_0, g_0_0_xxxxz_1, g_0_0_xxxxzz_0, g_0_0_xxxyy_0, g_0_0_xxxyy_1, g_0_0_xxxyyy_0, g_0_0_xxxyyz_0, g_0_0_xxxyzz_0, g_0_0_xxxzz_0, g_0_0_xxxzz_1, g_0_0_xxxzzz_0, g_0_0_xxyy_0, g_0_0_xxyy_1, g_0_0_xxyyy_0, g_0_0_xxyyy_1, g_0_0_xxyyyy_0, g_0_0_xxyyyz_0, g_0_0_xxyyzz_0, g_0_0_xxyzzz_0, g_0_0_xxzz_0, g_0_0_xxzz_1, g_0_0_xxzzz_0, g_0_0_xxzzz_1, g_0_0_xxzzzz_0, g_0_0_xyyy_0, g_0_0_xyyy_1, g_0_0_xyyyy_0, g_0_0_xyyyy_1, g_0_0_xyyyyy_0, g_0_0_xyyyyz_0, g_0_0_xyyyzz_0, g_0_0_xyyzz_0, g_0_0_xyyzz_1, g_0_0_xyyzzz_0, g_0_0_xyzzzz_0, g_0_0_xzzz_0, g_0_0_xzzz_1, g_0_0_xzzzz_0, g_0_0_xzzzz_1, g_0_0_xzzzzz_0, g_0_0_yyyy_0, g_0_0_yyyy_1, g_0_0_yyyyy_0, g_0_0_yyyyy_1, g_0_0_yyyyyy_0, g_0_0_yyyyyz_0, g_0_0_yyyyz_0, g_0_0_yyyyz_1, g_0_0_yyyyzz_0, g_0_0_yyyzz_0, g_0_0_yyyzz_1, g_0_0_yyyzzz_0, g_0_0_yyzz_0, g_0_0_yyzz_1, g_0_0_yyzzz_0, g_0_0_yyzzz_1, g_0_0_yyzzzz_0, g_0_0_yzzz_0, g_0_0_yzzz_1, g_0_0_yzzzz_0, g_0_0_yzzzz_1, g_0_0_yzzzzz_0, g_0_0_zzzz_0, g_0_0_zzzz_1, g_0_0_zzzzz_0, g_0_0_zzzzz_1, g_0_0_zzzzzz_0, qd_x, qd_y, qd_z, wq_x, wq_y, wq_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_cd_0 = 0.5 / (c_exps[i] + d_exps[i]);

        const double fzi_cd_0 = fi_cd_0 * a_exp / (a_exp + c_exps[i] + d_exps[i]) ;

        g_0_0_xxxxxx_0[i] = 5.0 * g_0_0_xxxx_0[i] * fi_cd_0 - 5.0 * g_0_0_xxxx_1[i] * fzi_cd_0 + g_0_0_xxxxx_0[i] * qd_x[i] + g_0_0_xxxxx_1[i] * wq_x[i];

        g_0_0_xxxxxy_0[i] = g_0_0_xxxxx_0[i] * qd_y[i] + g_0_0_xxxxx_1[i] * wq_y[i];

        g_0_0_xxxxxz_0[i] = g_0_0_xxxxx_0[i] * qd_z[i] + g_0_0_xxxxx_1[i] * wq_z[i];

        g_0_0_xxxxyy_0[i] = 3.0 * g_0_0_xxyy_0[i] * fi_cd_0 - 3.0 * g_0_0_xxyy_1[i] * fzi_cd_0 + g_0_0_xxxyy_0[i] * qd_x[i] + g_0_0_xxxyy_1[i] * wq_x[i];

        g_0_0_xxxxyz_0[i] = g_0_0_xxxxz_0[i] * qd_y[i] + g_0_0_xxxxz_1[i] * wq_y[i];

        g_0_0_xxxxzz_0[i] = 3.0 * g_0_0_xxzz_0[i] * fi_cd_0 - 3.0 * g_0_0_xxzz_1[i] * fzi_cd_0 + g_0_0_xxxzz_0[i] * qd_x[i] + g_0_0_xxxzz_1[i] * wq_x[i];

        g_0_0_xxxyyy_0[i] = 2.0 * g_0_0_xyyy_0[i] * fi_cd_0 - 2.0 * g_0_0_xyyy_1[i] * fzi_cd_0 + g_0_0_xxyyy_0[i] * qd_x[i] + g_0_0_xxyyy_1[i] * wq_x[i];

        g_0_0_xxxyyz_0[i] = g_0_0_xxxyy_0[i] * qd_z[i] + g_0_0_xxxyy_1[i] * wq_z[i];

        g_0_0_xxxyzz_0[i] = g_0_0_xxxzz_0[i] * qd_y[i] + g_0_0_xxxzz_1[i] * wq_y[i];

        g_0_0_xxxzzz_0[i] = 2.0 * g_0_0_xzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_xzzz_1[i] * fzi_cd_0 + g_0_0_xxzzz_0[i] * qd_x[i] + g_0_0_xxzzz_1[i] * wq_x[i];

        g_0_0_xxyyyy_0[i] = g_0_0_yyyy_0[i] * fi_cd_0 - g_0_0_yyyy_1[i] * fzi_cd_0 + g_0_0_xyyyy_0[i] * qd_x[i] + g_0_0_xyyyy_1[i] * wq_x[i];

        g_0_0_xxyyyz_0[i] = g_0_0_xxyyy_0[i] * qd_z[i] + g_0_0_xxyyy_1[i] * wq_z[i];

        g_0_0_xxyyzz_0[i] = g_0_0_yyzz_0[i] * fi_cd_0 - g_0_0_yyzz_1[i] * fzi_cd_0 + g_0_0_xyyzz_0[i] * qd_x[i] + g_0_0_xyyzz_1[i] * wq_x[i];

        g_0_0_xxyzzz_0[i] = g_0_0_xxzzz_0[i] * qd_y[i] + g_0_0_xxzzz_1[i] * wq_y[i];

        g_0_0_xxzzzz_0[i] = g_0_0_zzzz_0[i] * fi_cd_0 - g_0_0_zzzz_1[i] * fzi_cd_0 + g_0_0_xzzzz_0[i] * qd_x[i] + g_0_0_xzzzz_1[i] * wq_x[i];

        g_0_0_xyyyyy_0[i] = g_0_0_yyyyy_0[i] * qd_x[i] + g_0_0_yyyyy_1[i] * wq_x[i];

        g_0_0_xyyyyz_0[i] = g_0_0_yyyyz_0[i] * qd_x[i] + g_0_0_yyyyz_1[i] * wq_x[i];

        g_0_0_xyyyzz_0[i] = g_0_0_yyyzz_0[i] * qd_x[i] + g_0_0_yyyzz_1[i] * wq_x[i];

        g_0_0_xyyzzz_0[i] = g_0_0_yyzzz_0[i] * qd_x[i] + g_0_0_yyzzz_1[i] * wq_x[i];

        g_0_0_xyzzzz_0[i] = g_0_0_yzzzz_0[i] * qd_x[i] + g_0_0_yzzzz_1[i] * wq_x[i];

        g_0_0_xzzzzz_0[i] = g_0_0_zzzzz_0[i] * qd_x[i] + g_0_0_zzzzz_1[i] * wq_x[i];

        g_0_0_yyyyyy_0[i] = 5.0 * g_0_0_yyyy_0[i] * fi_cd_0 - 5.0 * g_0_0_yyyy_1[i] * fzi_cd_0 + g_0_0_yyyyy_0[i] * qd_y[i] + g_0_0_yyyyy_1[i] * wq_y[i];

        g_0_0_yyyyyz_0[i] = g_0_0_yyyyy_0[i] * qd_z[i] + g_0_0_yyyyy_1[i] * wq_z[i];

        g_0_0_yyyyzz_0[i] = 3.0 * g_0_0_yyzz_0[i] * fi_cd_0 - 3.0 * g_0_0_yyzz_1[i] * fzi_cd_0 + g_0_0_yyyzz_0[i] * qd_y[i] + g_0_0_yyyzz_1[i] * wq_y[i];

        g_0_0_yyyzzz_0[i] = 2.0 * g_0_0_yzzz_0[i] * fi_cd_0 - 2.0 * g_0_0_yzzz_1[i] * fzi_cd_0 + g_0_0_yyzzz_0[i] * qd_y[i] + g_0_0_yyzzz_1[i] * wq_y[i];

        g_0_0_yyzzzz_0[i] = g_0_0_zzzz_0[i] * fi_cd_0 - g_0_0_zzzz_1[i] * fzi_cd_0 + g_0_0_yzzzz_0[i] * qd_y[i] + g_0_0_yzzzz_1[i] * wq_y[i];

        g_0_0_yzzzzz_0[i] = g_0_0_zzzzz_0[i] * qd_y[i] + g_0_0_zzzzz_1[i] * wq_y[i];

        g_0_0_zzzzzz_0[i] = 5.0 * g_0_0_zzzz_0[i] * fi_cd_0 - 5.0 * g_0_0_zzzz_1[i] * fzi_cd_0 + g_0_0_zzzzz_0[i] * qd_z[i] + g_0_0_zzzzz_1[i] * wq_z[i];
    }
}

} // t3ceri namespace

