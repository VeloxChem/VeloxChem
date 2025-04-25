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

#include "ThreeCenterElectronRepulsionPrimRecPSH.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_psh(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_psh,
                                 size_t idx_eri_1_ssg,
                                 size_t idx_eri_1_ssh,
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

    /// Set up components of auxilary buffer : SSG

    auto g_0_0_xxxx_1 = pbuffer.data(idx_eri_1_ssg);

    auto g_0_0_xxxy_1 = pbuffer.data(idx_eri_1_ssg + 1);

    auto g_0_0_xxxz_1 = pbuffer.data(idx_eri_1_ssg + 2);

    auto g_0_0_xxyy_1 = pbuffer.data(idx_eri_1_ssg + 3);

    auto g_0_0_xxyz_1 = pbuffer.data(idx_eri_1_ssg + 4);

    auto g_0_0_xxzz_1 = pbuffer.data(idx_eri_1_ssg + 5);

    auto g_0_0_xyyy_1 = pbuffer.data(idx_eri_1_ssg + 6);

    auto g_0_0_xyyz_1 = pbuffer.data(idx_eri_1_ssg + 7);

    auto g_0_0_xyzz_1 = pbuffer.data(idx_eri_1_ssg + 8);

    auto g_0_0_xzzz_1 = pbuffer.data(idx_eri_1_ssg + 9);

    auto g_0_0_yyyy_1 = pbuffer.data(idx_eri_1_ssg + 10);

    auto g_0_0_yyyz_1 = pbuffer.data(idx_eri_1_ssg + 11);

    auto g_0_0_yyzz_1 = pbuffer.data(idx_eri_1_ssg + 12);

    auto g_0_0_yzzz_1 = pbuffer.data(idx_eri_1_ssg + 13);

    auto g_0_0_zzzz_1 = pbuffer.data(idx_eri_1_ssg + 14);

    /// Set up components of auxilary buffer : SSH

    auto g_0_0_xxxxx_1 = pbuffer.data(idx_eri_1_ssh);

    auto g_0_0_xxxxy_1 = pbuffer.data(idx_eri_1_ssh + 1);

    auto g_0_0_xxxxz_1 = pbuffer.data(idx_eri_1_ssh + 2);

    auto g_0_0_xxxyy_1 = pbuffer.data(idx_eri_1_ssh + 3);

    auto g_0_0_xxxyz_1 = pbuffer.data(idx_eri_1_ssh + 4);

    auto g_0_0_xxxzz_1 = pbuffer.data(idx_eri_1_ssh + 5);

    auto g_0_0_xxyyy_1 = pbuffer.data(idx_eri_1_ssh + 6);

    auto g_0_0_xxyyz_1 = pbuffer.data(idx_eri_1_ssh + 7);

    auto g_0_0_xxyzz_1 = pbuffer.data(idx_eri_1_ssh + 8);

    auto g_0_0_xxzzz_1 = pbuffer.data(idx_eri_1_ssh + 9);

    auto g_0_0_xyyyy_1 = pbuffer.data(idx_eri_1_ssh + 10);

    auto g_0_0_xyyyz_1 = pbuffer.data(idx_eri_1_ssh + 11);

    auto g_0_0_xyyzz_1 = pbuffer.data(idx_eri_1_ssh + 12);

    auto g_0_0_xyzzz_1 = pbuffer.data(idx_eri_1_ssh + 13);

    auto g_0_0_xzzzz_1 = pbuffer.data(idx_eri_1_ssh + 14);

    auto g_0_0_yyyyy_1 = pbuffer.data(idx_eri_1_ssh + 15);

    auto g_0_0_yyyyz_1 = pbuffer.data(idx_eri_1_ssh + 16);

    auto g_0_0_yyyzz_1 = pbuffer.data(idx_eri_1_ssh + 17);

    auto g_0_0_yyzzz_1 = pbuffer.data(idx_eri_1_ssh + 18);

    auto g_0_0_yzzzz_1 = pbuffer.data(idx_eri_1_ssh + 19);

    auto g_0_0_zzzzz_1 = pbuffer.data(idx_eri_1_ssh + 20);

    /// Set up 0-21 components of targeted buffer : PSH

    auto g_x_0_xxxxx_0 = pbuffer.data(idx_eri_0_psh);

    auto g_x_0_xxxxy_0 = pbuffer.data(idx_eri_0_psh + 1);

    auto g_x_0_xxxxz_0 = pbuffer.data(idx_eri_0_psh + 2);

    auto g_x_0_xxxyy_0 = pbuffer.data(idx_eri_0_psh + 3);

    auto g_x_0_xxxyz_0 = pbuffer.data(idx_eri_0_psh + 4);

    auto g_x_0_xxxzz_0 = pbuffer.data(idx_eri_0_psh + 5);

    auto g_x_0_xxyyy_0 = pbuffer.data(idx_eri_0_psh + 6);

    auto g_x_0_xxyyz_0 = pbuffer.data(idx_eri_0_psh + 7);

    auto g_x_0_xxyzz_0 = pbuffer.data(idx_eri_0_psh + 8);

    auto g_x_0_xxzzz_0 = pbuffer.data(idx_eri_0_psh + 9);

    auto g_x_0_xyyyy_0 = pbuffer.data(idx_eri_0_psh + 10);

    auto g_x_0_xyyyz_0 = pbuffer.data(idx_eri_0_psh + 11);

    auto g_x_0_xyyzz_0 = pbuffer.data(idx_eri_0_psh + 12);

    auto g_x_0_xyzzz_0 = pbuffer.data(idx_eri_0_psh + 13);

    auto g_x_0_xzzzz_0 = pbuffer.data(idx_eri_0_psh + 14);

    auto g_x_0_yyyyy_0 = pbuffer.data(idx_eri_0_psh + 15);

    auto g_x_0_yyyyz_0 = pbuffer.data(idx_eri_0_psh + 16);

    auto g_x_0_yyyzz_0 = pbuffer.data(idx_eri_0_psh + 17);

    auto g_x_0_yyzzz_0 = pbuffer.data(idx_eri_0_psh + 18);

    auto g_x_0_yzzzz_0 = pbuffer.data(idx_eri_0_psh + 19);

    auto g_x_0_zzzzz_0 = pbuffer.data(idx_eri_0_psh + 20);

    #pragma omp simd aligned(g_0_0_xxxx_1, g_0_0_xxxxx_1, g_0_0_xxxxy_1, g_0_0_xxxxz_1, g_0_0_xxxy_1, g_0_0_xxxyy_1, g_0_0_xxxyz_1, g_0_0_xxxz_1, g_0_0_xxxzz_1, g_0_0_xxyy_1, g_0_0_xxyyy_1, g_0_0_xxyyz_1, g_0_0_xxyz_1, g_0_0_xxyzz_1, g_0_0_xxzz_1, g_0_0_xxzzz_1, g_0_0_xyyy_1, g_0_0_xyyyy_1, g_0_0_xyyyz_1, g_0_0_xyyz_1, g_0_0_xyyzz_1, g_0_0_xyzz_1, g_0_0_xyzzz_1, g_0_0_xzzz_1, g_0_0_xzzzz_1, g_0_0_yyyy_1, g_0_0_yyyyy_1, g_0_0_yyyyz_1, g_0_0_yyyz_1, g_0_0_yyyzz_1, g_0_0_yyzz_1, g_0_0_yyzzz_1, g_0_0_yzzz_1, g_0_0_yzzzz_1, g_0_0_zzzz_1, g_0_0_zzzzz_1, g_x_0_xxxxx_0, g_x_0_xxxxy_0, g_x_0_xxxxz_0, g_x_0_xxxyy_0, g_x_0_xxxyz_0, g_x_0_xxxzz_0, g_x_0_xxyyy_0, g_x_0_xxyyz_0, g_x_0_xxyzz_0, g_x_0_xxzzz_0, g_x_0_xyyyy_0, g_x_0_xyyyz_0, g_x_0_xyyzz_0, g_x_0_xyzzz_0, g_x_0_xzzzz_0, g_x_0_yyyyy_0, g_x_0_yyyyz_0, g_x_0_yyyzz_0, g_x_0_yyzzz_0, g_x_0_yzzzz_0, g_x_0_zzzzz_0, wa_x  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_x_0_xxxxx_0[i] = 5.0 * g_0_0_xxxx_1[i] * fi_acd_0 + g_0_0_xxxxx_1[i] * wa_x[i];

        g_x_0_xxxxy_0[i] = 4.0 * g_0_0_xxxy_1[i] * fi_acd_0 + g_0_0_xxxxy_1[i] * wa_x[i];

        g_x_0_xxxxz_0[i] = 4.0 * g_0_0_xxxz_1[i] * fi_acd_0 + g_0_0_xxxxz_1[i] * wa_x[i];

        g_x_0_xxxyy_0[i] = 3.0 * g_0_0_xxyy_1[i] * fi_acd_0 + g_0_0_xxxyy_1[i] * wa_x[i];

        g_x_0_xxxyz_0[i] = 3.0 * g_0_0_xxyz_1[i] * fi_acd_0 + g_0_0_xxxyz_1[i] * wa_x[i];

        g_x_0_xxxzz_0[i] = 3.0 * g_0_0_xxzz_1[i] * fi_acd_0 + g_0_0_xxxzz_1[i] * wa_x[i];

        g_x_0_xxyyy_0[i] = 2.0 * g_0_0_xyyy_1[i] * fi_acd_0 + g_0_0_xxyyy_1[i] * wa_x[i];

        g_x_0_xxyyz_0[i] = 2.0 * g_0_0_xyyz_1[i] * fi_acd_0 + g_0_0_xxyyz_1[i] * wa_x[i];

        g_x_0_xxyzz_0[i] = 2.0 * g_0_0_xyzz_1[i] * fi_acd_0 + g_0_0_xxyzz_1[i] * wa_x[i];

        g_x_0_xxzzz_0[i] = 2.0 * g_0_0_xzzz_1[i] * fi_acd_0 + g_0_0_xxzzz_1[i] * wa_x[i];

        g_x_0_xyyyy_0[i] = g_0_0_yyyy_1[i] * fi_acd_0 + g_0_0_xyyyy_1[i] * wa_x[i];

        g_x_0_xyyyz_0[i] = g_0_0_yyyz_1[i] * fi_acd_0 + g_0_0_xyyyz_1[i] * wa_x[i];

        g_x_0_xyyzz_0[i] = g_0_0_yyzz_1[i] * fi_acd_0 + g_0_0_xyyzz_1[i] * wa_x[i];

        g_x_0_xyzzz_0[i] = g_0_0_yzzz_1[i] * fi_acd_0 + g_0_0_xyzzz_1[i] * wa_x[i];

        g_x_0_xzzzz_0[i] = g_0_0_zzzz_1[i] * fi_acd_0 + g_0_0_xzzzz_1[i] * wa_x[i];

        g_x_0_yyyyy_0[i] = g_0_0_yyyyy_1[i] * wa_x[i];

        g_x_0_yyyyz_0[i] = g_0_0_yyyyz_1[i] * wa_x[i];

        g_x_0_yyyzz_0[i] = g_0_0_yyyzz_1[i] * wa_x[i];

        g_x_0_yyzzz_0[i] = g_0_0_yyzzz_1[i] * wa_x[i];

        g_x_0_yzzzz_0[i] = g_0_0_yzzzz_1[i] * wa_x[i];

        g_x_0_zzzzz_0[i] = g_0_0_zzzzz_1[i] * wa_x[i];
    }

    /// Set up 21-42 components of targeted buffer : PSH

    auto g_y_0_xxxxx_0 = pbuffer.data(idx_eri_0_psh + 21);

    auto g_y_0_xxxxy_0 = pbuffer.data(idx_eri_0_psh + 22);

    auto g_y_0_xxxxz_0 = pbuffer.data(idx_eri_0_psh + 23);

    auto g_y_0_xxxyy_0 = pbuffer.data(idx_eri_0_psh + 24);

    auto g_y_0_xxxyz_0 = pbuffer.data(idx_eri_0_psh + 25);

    auto g_y_0_xxxzz_0 = pbuffer.data(idx_eri_0_psh + 26);

    auto g_y_0_xxyyy_0 = pbuffer.data(idx_eri_0_psh + 27);

    auto g_y_0_xxyyz_0 = pbuffer.data(idx_eri_0_psh + 28);

    auto g_y_0_xxyzz_0 = pbuffer.data(idx_eri_0_psh + 29);

    auto g_y_0_xxzzz_0 = pbuffer.data(idx_eri_0_psh + 30);

    auto g_y_0_xyyyy_0 = pbuffer.data(idx_eri_0_psh + 31);

    auto g_y_0_xyyyz_0 = pbuffer.data(idx_eri_0_psh + 32);

    auto g_y_0_xyyzz_0 = pbuffer.data(idx_eri_0_psh + 33);

    auto g_y_0_xyzzz_0 = pbuffer.data(idx_eri_0_psh + 34);

    auto g_y_0_xzzzz_0 = pbuffer.data(idx_eri_0_psh + 35);

    auto g_y_0_yyyyy_0 = pbuffer.data(idx_eri_0_psh + 36);

    auto g_y_0_yyyyz_0 = pbuffer.data(idx_eri_0_psh + 37);

    auto g_y_0_yyyzz_0 = pbuffer.data(idx_eri_0_psh + 38);

    auto g_y_0_yyzzz_0 = pbuffer.data(idx_eri_0_psh + 39);

    auto g_y_0_yzzzz_0 = pbuffer.data(idx_eri_0_psh + 40);

    auto g_y_0_zzzzz_0 = pbuffer.data(idx_eri_0_psh + 41);

    #pragma omp simd aligned(g_0_0_xxxx_1, g_0_0_xxxxx_1, g_0_0_xxxxy_1, g_0_0_xxxxz_1, g_0_0_xxxy_1, g_0_0_xxxyy_1, g_0_0_xxxyz_1, g_0_0_xxxz_1, g_0_0_xxxzz_1, g_0_0_xxyy_1, g_0_0_xxyyy_1, g_0_0_xxyyz_1, g_0_0_xxyz_1, g_0_0_xxyzz_1, g_0_0_xxzz_1, g_0_0_xxzzz_1, g_0_0_xyyy_1, g_0_0_xyyyy_1, g_0_0_xyyyz_1, g_0_0_xyyz_1, g_0_0_xyyzz_1, g_0_0_xyzz_1, g_0_0_xyzzz_1, g_0_0_xzzz_1, g_0_0_xzzzz_1, g_0_0_yyyy_1, g_0_0_yyyyy_1, g_0_0_yyyyz_1, g_0_0_yyyz_1, g_0_0_yyyzz_1, g_0_0_yyzz_1, g_0_0_yyzzz_1, g_0_0_yzzz_1, g_0_0_yzzzz_1, g_0_0_zzzz_1, g_0_0_zzzzz_1, g_y_0_xxxxx_0, g_y_0_xxxxy_0, g_y_0_xxxxz_0, g_y_0_xxxyy_0, g_y_0_xxxyz_0, g_y_0_xxxzz_0, g_y_0_xxyyy_0, g_y_0_xxyyz_0, g_y_0_xxyzz_0, g_y_0_xxzzz_0, g_y_0_xyyyy_0, g_y_0_xyyyz_0, g_y_0_xyyzz_0, g_y_0_xyzzz_0, g_y_0_xzzzz_0, g_y_0_yyyyy_0, g_y_0_yyyyz_0, g_y_0_yyyzz_0, g_y_0_yyzzz_0, g_y_0_yzzzz_0, g_y_0_zzzzz_0, wa_y  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_y_0_xxxxx_0[i] = g_0_0_xxxxx_1[i] * wa_y[i];

        g_y_0_xxxxy_0[i] = g_0_0_xxxx_1[i] * fi_acd_0 + g_0_0_xxxxy_1[i] * wa_y[i];

        g_y_0_xxxxz_0[i] = g_0_0_xxxxz_1[i] * wa_y[i];

        g_y_0_xxxyy_0[i] = 2.0 * g_0_0_xxxy_1[i] * fi_acd_0 + g_0_0_xxxyy_1[i] * wa_y[i];

        g_y_0_xxxyz_0[i] = g_0_0_xxxz_1[i] * fi_acd_0 + g_0_0_xxxyz_1[i] * wa_y[i];

        g_y_0_xxxzz_0[i] = g_0_0_xxxzz_1[i] * wa_y[i];

        g_y_0_xxyyy_0[i] = 3.0 * g_0_0_xxyy_1[i] * fi_acd_0 + g_0_0_xxyyy_1[i] * wa_y[i];

        g_y_0_xxyyz_0[i] = 2.0 * g_0_0_xxyz_1[i] * fi_acd_0 + g_0_0_xxyyz_1[i] * wa_y[i];

        g_y_0_xxyzz_0[i] = g_0_0_xxzz_1[i] * fi_acd_0 + g_0_0_xxyzz_1[i] * wa_y[i];

        g_y_0_xxzzz_0[i] = g_0_0_xxzzz_1[i] * wa_y[i];

        g_y_0_xyyyy_0[i] = 4.0 * g_0_0_xyyy_1[i] * fi_acd_0 + g_0_0_xyyyy_1[i] * wa_y[i];

        g_y_0_xyyyz_0[i] = 3.0 * g_0_0_xyyz_1[i] * fi_acd_0 + g_0_0_xyyyz_1[i] * wa_y[i];

        g_y_0_xyyzz_0[i] = 2.0 * g_0_0_xyzz_1[i] * fi_acd_0 + g_0_0_xyyzz_1[i] * wa_y[i];

        g_y_0_xyzzz_0[i] = g_0_0_xzzz_1[i] * fi_acd_0 + g_0_0_xyzzz_1[i] * wa_y[i];

        g_y_0_xzzzz_0[i] = g_0_0_xzzzz_1[i] * wa_y[i];

        g_y_0_yyyyy_0[i] = 5.0 * g_0_0_yyyy_1[i] * fi_acd_0 + g_0_0_yyyyy_1[i] * wa_y[i];

        g_y_0_yyyyz_0[i] = 4.0 * g_0_0_yyyz_1[i] * fi_acd_0 + g_0_0_yyyyz_1[i] * wa_y[i];

        g_y_0_yyyzz_0[i] = 3.0 * g_0_0_yyzz_1[i] * fi_acd_0 + g_0_0_yyyzz_1[i] * wa_y[i];

        g_y_0_yyzzz_0[i] = 2.0 * g_0_0_yzzz_1[i] * fi_acd_0 + g_0_0_yyzzz_1[i] * wa_y[i];

        g_y_0_yzzzz_0[i] = g_0_0_zzzz_1[i] * fi_acd_0 + g_0_0_yzzzz_1[i] * wa_y[i];

        g_y_0_zzzzz_0[i] = g_0_0_zzzzz_1[i] * wa_y[i];
    }

    /// Set up 42-63 components of targeted buffer : PSH

    auto g_z_0_xxxxx_0 = pbuffer.data(idx_eri_0_psh + 42);

    auto g_z_0_xxxxy_0 = pbuffer.data(idx_eri_0_psh + 43);

    auto g_z_0_xxxxz_0 = pbuffer.data(idx_eri_0_psh + 44);

    auto g_z_0_xxxyy_0 = pbuffer.data(idx_eri_0_psh + 45);

    auto g_z_0_xxxyz_0 = pbuffer.data(idx_eri_0_psh + 46);

    auto g_z_0_xxxzz_0 = pbuffer.data(idx_eri_0_psh + 47);

    auto g_z_0_xxyyy_0 = pbuffer.data(idx_eri_0_psh + 48);

    auto g_z_0_xxyyz_0 = pbuffer.data(idx_eri_0_psh + 49);

    auto g_z_0_xxyzz_0 = pbuffer.data(idx_eri_0_psh + 50);

    auto g_z_0_xxzzz_0 = pbuffer.data(idx_eri_0_psh + 51);

    auto g_z_0_xyyyy_0 = pbuffer.data(idx_eri_0_psh + 52);

    auto g_z_0_xyyyz_0 = pbuffer.data(idx_eri_0_psh + 53);

    auto g_z_0_xyyzz_0 = pbuffer.data(idx_eri_0_psh + 54);

    auto g_z_0_xyzzz_0 = pbuffer.data(idx_eri_0_psh + 55);

    auto g_z_0_xzzzz_0 = pbuffer.data(idx_eri_0_psh + 56);

    auto g_z_0_yyyyy_0 = pbuffer.data(idx_eri_0_psh + 57);

    auto g_z_0_yyyyz_0 = pbuffer.data(idx_eri_0_psh + 58);

    auto g_z_0_yyyzz_0 = pbuffer.data(idx_eri_0_psh + 59);

    auto g_z_0_yyzzz_0 = pbuffer.data(idx_eri_0_psh + 60);

    auto g_z_0_yzzzz_0 = pbuffer.data(idx_eri_0_psh + 61);

    auto g_z_0_zzzzz_0 = pbuffer.data(idx_eri_0_psh + 62);

    #pragma omp simd aligned(g_0_0_xxxx_1, g_0_0_xxxxx_1, g_0_0_xxxxy_1, g_0_0_xxxxz_1, g_0_0_xxxy_1, g_0_0_xxxyy_1, g_0_0_xxxyz_1, g_0_0_xxxz_1, g_0_0_xxxzz_1, g_0_0_xxyy_1, g_0_0_xxyyy_1, g_0_0_xxyyz_1, g_0_0_xxyz_1, g_0_0_xxyzz_1, g_0_0_xxzz_1, g_0_0_xxzzz_1, g_0_0_xyyy_1, g_0_0_xyyyy_1, g_0_0_xyyyz_1, g_0_0_xyyz_1, g_0_0_xyyzz_1, g_0_0_xyzz_1, g_0_0_xyzzz_1, g_0_0_xzzz_1, g_0_0_xzzzz_1, g_0_0_yyyy_1, g_0_0_yyyyy_1, g_0_0_yyyyz_1, g_0_0_yyyz_1, g_0_0_yyyzz_1, g_0_0_yyzz_1, g_0_0_yyzzz_1, g_0_0_yzzz_1, g_0_0_yzzzz_1, g_0_0_zzzz_1, g_0_0_zzzzz_1, g_z_0_xxxxx_0, g_z_0_xxxxy_0, g_z_0_xxxxz_0, g_z_0_xxxyy_0, g_z_0_xxxyz_0, g_z_0_xxxzz_0, g_z_0_xxyyy_0, g_z_0_xxyyz_0, g_z_0_xxyzz_0, g_z_0_xxzzz_0, g_z_0_xyyyy_0, g_z_0_xyyyz_0, g_z_0_xyyzz_0, g_z_0_xyzzz_0, g_z_0_xzzzz_0, g_z_0_yyyyy_0, g_z_0_yyyyz_0, g_z_0_yyyzz_0, g_z_0_yyzzz_0, g_z_0_yzzzz_0, g_z_0_zzzzz_0, wa_z  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_z_0_xxxxx_0[i] = g_0_0_xxxxx_1[i] * wa_z[i];

        g_z_0_xxxxy_0[i] = g_0_0_xxxxy_1[i] * wa_z[i];

        g_z_0_xxxxz_0[i] = g_0_0_xxxx_1[i] * fi_acd_0 + g_0_0_xxxxz_1[i] * wa_z[i];

        g_z_0_xxxyy_0[i] = g_0_0_xxxyy_1[i] * wa_z[i];

        g_z_0_xxxyz_0[i] = g_0_0_xxxy_1[i] * fi_acd_0 + g_0_0_xxxyz_1[i] * wa_z[i];

        g_z_0_xxxzz_0[i] = 2.0 * g_0_0_xxxz_1[i] * fi_acd_0 + g_0_0_xxxzz_1[i] * wa_z[i];

        g_z_0_xxyyy_0[i] = g_0_0_xxyyy_1[i] * wa_z[i];

        g_z_0_xxyyz_0[i] = g_0_0_xxyy_1[i] * fi_acd_0 + g_0_0_xxyyz_1[i] * wa_z[i];

        g_z_0_xxyzz_0[i] = 2.0 * g_0_0_xxyz_1[i] * fi_acd_0 + g_0_0_xxyzz_1[i] * wa_z[i];

        g_z_0_xxzzz_0[i] = 3.0 * g_0_0_xxzz_1[i] * fi_acd_0 + g_0_0_xxzzz_1[i] * wa_z[i];

        g_z_0_xyyyy_0[i] = g_0_0_xyyyy_1[i] * wa_z[i];

        g_z_0_xyyyz_0[i] = g_0_0_xyyy_1[i] * fi_acd_0 + g_0_0_xyyyz_1[i] * wa_z[i];

        g_z_0_xyyzz_0[i] = 2.0 * g_0_0_xyyz_1[i] * fi_acd_0 + g_0_0_xyyzz_1[i] * wa_z[i];

        g_z_0_xyzzz_0[i] = 3.0 * g_0_0_xyzz_1[i] * fi_acd_0 + g_0_0_xyzzz_1[i] * wa_z[i];

        g_z_0_xzzzz_0[i] = 4.0 * g_0_0_xzzz_1[i] * fi_acd_0 + g_0_0_xzzzz_1[i] * wa_z[i];

        g_z_0_yyyyy_0[i] = g_0_0_yyyyy_1[i] * wa_z[i];

        g_z_0_yyyyz_0[i] = g_0_0_yyyy_1[i] * fi_acd_0 + g_0_0_yyyyz_1[i] * wa_z[i];

        g_z_0_yyyzz_0[i] = 2.0 * g_0_0_yyyz_1[i] * fi_acd_0 + g_0_0_yyyzz_1[i] * wa_z[i];

        g_z_0_yyzzz_0[i] = 3.0 * g_0_0_yyzz_1[i] * fi_acd_0 + g_0_0_yyzzz_1[i] * wa_z[i];

        g_z_0_yzzzz_0[i] = 4.0 * g_0_0_yzzz_1[i] * fi_acd_0 + g_0_0_yzzzz_1[i] * wa_z[i];

        g_z_0_zzzzz_0[i] = 5.0 * g_0_0_zzzz_1[i] * fi_acd_0 + g_0_0_zzzzz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

