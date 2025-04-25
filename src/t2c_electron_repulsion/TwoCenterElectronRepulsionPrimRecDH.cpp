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

#include "TwoCenterElectronRepulsionPrimRecDH.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_dh(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_dh,
                                const size_t idx_eri_0_sh,
                                const size_t idx_eri_1_sh,
                                const size_t idx_eri_1_pg,
                                const size_t idx_eri_1_ph,
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

    // Set up components of auxiliary buffer : SH

    auto g_0_xxxxx_0 = pbuffer.data(idx_eri_0_sh);

    auto g_0_xxxxy_0 = pbuffer.data(idx_eri_0_sh + 1);

    auto g_0_xxxxz_0 = pbuffer.data(idx_eri_0_sh + 2);

    auto g_0_xxxyy_0 = pbuffer.data(idx_eri_0_sh + 3);

    auto g_0_xxxyz_0 = pbuffer.data(idx_eri_0_sh + 4);

    auto g_0_xxxzz_0 = pbuffer.data(idx_eri_0_sh + 5);

    auto g_0_xxyyy_0 = pbuffer.data(idx_eri_0_sh + 6);

    auto g_0_xxyyz_0 = pbuffer.data(idx_eri_0_sh + 7);

    auto g_0_xxyzz_0 = pbuffer.data(idx_eri_0_sh + 8);

    auto g_0_xxzzz_0 = pbuffer.data(idx_eri_0_sh + 9);

    auto g_0_xyyyy_0 = pbuffer.data(idx_eri_0_sh + 10);

    auto g_0_xyyyz_0 = pbuffer.data(idx_eri_0_sh + 11);

    auto g_0_xyyzz_0 = pbuffer.data(idx_eri_0_sh + 12);

    auto g_0_xyzzz_0 = pbuffer.data(idx_eri_0_sh + 13);

    auto g_0_xzzzz_0 = pbuffer.data(idx_eri_0_sh + 14);

    auto g_0_yyyyy_0 = pbuffer.data(idx_eri_0_sh + 15);

    auto g_0_yyyyz_0 = pbuffer.data(idx_eri_0_sh + 16);

    auto g_0_yyyzz_0 = pbuffer.data(idx_eri_0_sh + 17);

    auto g_0_yyzzz_0 = pbuffer.data(idx_eri_0_sh + 18);

    auto g_0_yzzzz_0 = pbuffer.data(idx_eri_0_sh + 19);

    auto g_0_zzzzz_0 = pbuffer.data(idx_eri_0_sh + 20);

    // Set up components of auxiliary buffer : SH

    auto g_0_xxxxx_1 = pbuffer.data(idx_eri_1_sh);

    auto g_0_xxxxy_1 = pbuffer.data(idx_eri_1_sh + 1);

    auto g_0_xxxxz_1 = pbuffer.data(idx_eri_1_sh + 2);

    auto g_0_xxxyy_1 = pbuffer.data(idx_eri_1_sh + 3);

    auto g_0_xxxyz_1 = pbuffer.data(idx_eri_1_sh + 4);

    auto g_0_xxxzz_1 = pbuffer.data(idx_eri_1_sh + 5);

    auto g_0_xxyyy_1 = pbuffer.data(idx_eri_1_sh + 6);

    auto g_0_xxyyz_1 = pbuffer.data(idx_eri_1_sh + 7);

    auto g_0_xxyzz_1 = pbuffer.data(idx_eri_1_sh + 8);

    auto g_0_xxzzz_1 = pbuffer.data(idx_eri_1_sh + 9);

    auto g_0_xyyyy_1 = pbuffer.data(idx_eri_1_sh + 10);

    auto g_0_xyyyz_1 = pbuffer.data(idx_eri_1_sh + 11);

    auto g_0_xyyzz_1 = pbuffer.data(idx_eri_1_sh + 12);

    auto g_0_xyzzz_1 = pbuffer.data(idx_eri_1_sh + 13);

    auto g_0_xzzzz_1 = pbuffer.data(idx_eri_1_sh + 14);

    auto g_0_yyyyy_1 = pbuffer.data(idx_eri_1_sh + 15);

    auto g_0_yyyyz_1 = pbuffer.data(idx_eri_1_sh + 16);

    auto g_0_yyyzz_1 = pbuffer.data(idx_eri_1_sh + 17);

    auto g_0_yyzzz_1 = pbuffer.data(idx_eri_1_sh + 18);

    auto g_0_yzzzz_1 = pbuffer.data(idx_eri_1_sh + 19);

    auto g_0_zzzzz_1 = pbuffer.data(idx_eri_1_sh + 20);

    // Set up components of auxiliary buffer : PG

    auto g_x_xxxx_1 = pbuffer.data(idx_eri_1_pg);

    auto g_x_xxxy_1 = pbuffer.data(idx_eri_1_pg + 1);

    auto g_x_xxxz_1 = pbuffer.data(idx_eri_1_pg + 2);

    auto g_x_xxyy_1 = pbuffer.data(idx_eri_1_pg + 3);

    auto g_x_xxyz_1 = pbuffer.data(idx_eri_1_pg + 4);

    auto g_x_xxzz_1 = pbuffer.data(idx_eri_1_pg + 5);

    auto g_x_xyyy_1 = pbuffer.data(idx_eri_1_pg + 6);

    auto g_x_xyyz_1 = pbuffer.data(idx_eri_1_pg + 7);

    auto g_x_xyzz_1 = pbuffer.data(idx_eri_1_pg + 8);

    auto g_x_xzzz_1 = pbuffer.data(idx_eri_1_pg + 9);

    auto g_x_yyyy_1 = pbuffer.data(idx_eri_1_pg + 10);

    auto g_x_yyyz_1 = pbuffer.data(idx_eri_1_pg + 11);

    auto g_x_yyzz_1 = pbuffer.data(idx_eri_1_pg + 12);

    auto g_x_yzzz_1 = pbuffer.data(idx_eri_1_pg + 13);

    auto g_x_zzzz_1 = pbuffer.data(idx_eri_1_pg + 14);

    auto g_y_xxxx_1 = pbuffer.data(idx_eri_1_pg + 15);

    auto g_y_xxxy_1 = pbuffer.data(idx_eri_1_pg + 16);

    auto g_y_xxxz_1 = pbuffer.data(idx_eri_1_pg + 17);

    auto g_y_xxyy_1 = pbuffer.data(idx_eri_1_pg + 18);

    auto g_y_xxyz_1 = pbuffer.data(idx_eri_1_pg + 19);

    auto g_y_xxzz_1 = pbuffer.data(idx_eri_1_pg + 20);

    auto g_y_xyyy_1 = pbuffer.data(idx_eri_1_pg + 21);

    auto g_y_xyyz_1 = pbuffer.data(idx_eri_1_pg + 22);

    auto g_y_xyzz_1 = pbuffer.data(idx_eri_1_pg + 23);

    auto g_y_xzzz_1 = pbuffer.data(idx_eri_1_pg + 24);

    auto g_y_yyyy_1 = pbuffer.data(idx_eri_1_pg + 25);

    auto g_y_yyyz_1 = pbuffer.data(idx_eri_1_pg + 26);

    auto g_y_yyzz_1 = pbuffer.data(idx_eri_1_pg + 27);

    auto g_y_yzzz_1 = pbuffer.data(idx_eri_1_pg + 28);

    auto g_y_zzzz_1 = pbuffer.data(idx_eri_1_pg + 29);

    auto g_z_xxxx_1 = pbuffer.data(idx_eri_1_pg + 30);

    auto g_z_xxxy_1 = pbuffer.data(idx_eri_1_pg + 31);

    auto g_z_xxxz_1 = pbuffer.data(idx_eri_1_pg + 32);

    auto g_z_xxyy_1 = pbuffer.data(idx_eri_1_pg + 33);

    auto g_z_xxyz_1 = pbuffer.data(idx_eri_1_pg + 34);

    auto g_z_xxzz_1 = pbuffer.data(idx_eri_1_pg + 35);

    auto g_z_xyyy_1 = pbuffer.data(idx_eri_1_pg + 36);

    auto g_z_xyyz_1 = pbuffer.data(idx_eri_1_pg + 37);

    auto g_z_xyzz_1 = pbuffer.data(idx_eri_1_pg + 38);

    auto g_z_xzzz_1 = pbuffer.data(idx_eri_1_pg + 39);

    auto g_z_yyyy_1 = pbuffer.data(idx_eri_1_pg + 40);

    auto g_z_yyyz_1 = pbuffer.data(idx_eri_1_pg + 41);

    auto g_z_yyzz_1 = pbuffer.data(idx_eri_1_pg + 42);

    auto g_z_yzzz_1 = pbuffer.data(idx_eri_1_pg + 43);

    auto g_z_zzzz_1 = pbuffer.data(idx_eri_1_pg + 44);

    // Set up components of auxiliary buffer : PH

    auto g_x_xxxxx_1 = pbuffer.data(idx_eri_1_ph);

    auto g_x_xxxxy_1 = pbuffer.data(idx_eri_1_ph + 1);

    auto g_x_xxxxz_1 = pbuffer.data(idx_eri_1_ph + 2);

    auto g_x_xxxyy_1 = pbuffer.data(idx_eri_1_ph + 3);

    auto g_x_xxxyz_1 = pbuffer.data(idx_eri_1_ph + 4);

    auto g_x_xxxzz_1 = pbuffer.data(idx_eri_1_ph + 5);

    auto g_x_xxyyy_1 = pbuffer.data(idx_eri_1_ph + 6);

    auto g_x_xxyyz_1 = pbuffer.data(idx_eri_1_ph + 7);

    auto g_x_xxyzz_1 = pbuffer.data(idx_eri_1_ph + 8);

    auto g_x_xxzzz_1 = pbuffer.data(idx_eri_1_ph + 9);

    auto g_x_xyyyy_1 = pbuffer.data(idx_eri_1_ph + 10);

    auto g_x_xyyyz_1 = pbuffer.data(idx_eri_1_ph + 11);

    auto g_x_xyyzz_1 = pbuffer.data(idx_eri_1_ph + 12);

    auto g_x_xyzzz_1 = pbuffer.data(idx_eri_1_ph + 13);

    auto g_x_xzzzz_1 = pbuffer.data(idx_eri_1_ph + 14);

    auto g_x_yyyyy_1 = pbuffer.data(idx_eri_1_ph + 15);

    auto g_x_yyyyz_1 = pbuffer.data(idx_eri_1_ph + 16);

    auto g_x_yyyzz_1 = pbuffer.data(idx_eri_1_ph + 17);

    auto g_x_yyzzz_1 = pbuffer.data(idx_eri_1_ph + 18);

    auto g_x_yzzzz_1 = pbuffer.data(idx_eri_1_ph + 19);

    auto g_x_zzzzz_1 = pbuffer.data(idx_eri_1_ph + 20);

    auto g_y_xxxxx_1 = pbuffer.data(idx_eri_1_ph + 21);

    auto g_y_xxxxy_1 = pbuffer.data(idx_eri_1_ph + 22);

    auto g_y_xxxxz_1 = pbuffer.data(idx_eri_1_ph + 23);

    auto g_y_xxxyy_1 = pbuffer.data(idx_eri_1_ph + 24);

    auto g_y_xxxyz_1 = pbuffer.data(idx_eri_1_ph + 25);

    auto g_y_xxxzz_1 = pbuffer.data(idx_eri_1_ph + 26);

    auto g_y_xxyyy_1 = pbuffer.data(idx_eri_1_ph + 27);

    auto g_y_xxyyz_1 = pbuffer.data(idx_eri_1_ph + 28);

    auto g_y_xxyzz_1 = pbuffer.data(idx_eri_1_ph + 29);

    auto g_y_xxzzz_1 = pbuffer.data(idx_eri_1_ph + 30);

    auto g_y_xyyyy_1 = pbuffer.data(idx_eri_1_ph + 31);

    auto g_y_xyyyz_1 = pbuffer.data(idx_eri_1_ph + 32);

    auto g_y_xyyzz_1 = pbuffer.data(idx_eri_1_ph + 33);

    auto g_y_xyzzz_1 = pbuffer.data(idx_eri_1_ph + 34);

    auto g_y_xzzzz_1 = pbuffer.data(idx_eri_1_ph + 35);

    auto g_y_yyyyy_1 = pbuffer.data(idx_eri_1_ph + 36);

    auto g_y_yyyyz_1 = pbuffer.data(idx_eri_1_ph + 37);

    auto g_y_yyyzz_1 = pbuffer.data(idx_eri_1_ph + 38);

    auto g_y_yyzzz_1 = pbuffer.data(idx_eri_1_ph + 39);

    auto g_y_yzzzz_1 = pbuffer.data(idx_eri_1_ph + 40);

    auto g_y_zzzzz_1 = pbuffer.data(idx_eri_1_ph + 41);

    auto g_z_xxxxx_1 = pbuffer.data(idx_eri_1_ph + 42);

    auto g_z_xxxxy_1 = pbuffer.data(idx_eri_1_ph + 43);

    auto g_z_xxxxz_1 = pbuffer.data(idx_eri_1_ph + 44);

    auto g_z_xxxyy_1 = pbuffer.data(idx_eri_1_ph + 45);

    auto g_z_xxxyz_1 = pbuffer.data(idx_eri_1_ph + 46);

    auto g_z_xxxzz_1 = pbuffer.data(idx_eri_1_ph + 47);

    auto g_z_xxyyy_1 = pbuffer.data(idx_eri_1_ph + 48);

    auto g_z_xxyyz_1 = pbuffer.data(idx_eri_1_ph + 49);

    auto g_z_xxyzz_1 = pbuffer.data(idx_eri_1_ph + 50);

    auto g_z_xxzzz_1 = pbuffer.data(idx_eri_1_ph + 51);

    auto g_z_xyyyy_1 = pbuffer.data(idx_eri_1_ph + 52);

    auto g_z_xyyyz_1 = pbuffer.data(idx_eri_1_ph + 53);

    auto g_z_xyyzz_1 = pbuffer.data(idx_eri_1_ph + 54);

    auto g_z_xyzzz_1 = pbuffer.data(idx_eri_1_ph + 55);

    auto g_z_xzzzz_1 = pbuffer.data(idx_eri_1_ph + 56);

    auto g_z_yyyyy_1 = pbuffer.data(idx_eri_1_ph + 57);

    auto g_z_yyyyz_1 = pbuffer.data(idx_eri_1_ph + 58);

    auto g_z_yyyzz_1 = pbuffer.data(idx_eri_1_ph + 59);

    auto g_z_yyzzz_1 = pbuffer.data(idx_eri_1_ph + 60);

    auto g_z_yzzzz_1 = pbuffer.data(idx_eri_1_ph + 61);

    auto g_z_zzzzz_1 = pbuffer.data(idx_eri_1_ph + 62);

    // Set up 0-21 components of targeted buffer : DH

    auto g_xx_xxxxx_0 = pbuffer.data(idx_eri_0_dh);

    auto g_xx_xxxxy_0 = pbuffer.data(idx_eri_0_dh + 1);

    auto g_xx_xxxxz_0 = pbuffer.data(idx_eri_0_dh + 2);

    auto g_xx_xxxyy_0 = pbuffer.data(idx_eri_0_dh + 3);

    auto g_xx_xxxyz_0 = pbuffer.data(idx_eri_0_dh + 4);

    auto g_xx_xxxzz_0 = pbuffer.data(idx_eri_0_dh + 5);

    auto g_xx_xxyyy_0 = pbuffer.data(idx_eri_0_dh + 6);

    auto g_xx_xxyyz_0 = pbuffer.data(idx_eri_0_dh + 7);

    auto g_xx_xxyzz_0 = pbuffer.data(idx_eri_0_dh + 8);

    auto g_xx_xxzzz_0 = pbuffer.data(idx_eri_0_dh + 9);

    auto g_xx_xyyyy_0 = pbuffer.data(idx_eri_0_dh + 10);

    auto g_xx_xyyyz_0 = pbuffer.data(idx_eri_0_dh + 11);

    auto g_xx_xyyzz_0 = pbuffer.data(idx_eri_0_dh + 12);

    auto g_xx_xyzzz_0 = pbuffer.data(idx_eri_0_dh + 13);

    auto g_xx_xzzzz_0 = pbuffer.data(idx_eri_0_dh + 14);

    auto g_xx_yyyyy_0 = pbuffer.data(idx_eri_0_dh + 15);

    auto g_xx_yyyyz_0 = pbuffer.data(idx_eri_0_dh + 16);

    auto g_xx_yyyzz_0 = pbuffer.data(idx_eri_0_dh + 17);

    auto g_xx_yyzzz_0 = pbuffer.data(idx_eri_0_dh + 18);

    auto g_xx_yzzzz_0 = pbuffer.data(idx_eri_0_dh + 19);

    auto g_xx_zzzzz_0 = pbuffer.data(idx_eri_0_dh + 20);

    #pragma omp simd aligned(g_0_xxxxx_0, g_0_xxxxx_1, g_0_xxxxy_0, g_0_xxxxy_1, g_0_xxxxz_0, g_0_xxxxz_1, g_0_xxxyy_0, g_0_xxxyy_1, g_0_xxxyz_0, g_0_xxxyz_1, g_0_xxxzz_0, g_0_xxxzz_1, g_0_xxyyy_0, g_0_xxyyy_1, g_0_xxyyz_0, g_0_xxyyz_1, g_0_xxyzz_0, g_0_xxyzz_1, g_0_xxzzz_0, g_0_xxzzz_1, g_0_xyyyy_0, g_0_xyyyy_1, g_0_xyyyz_0, g_0_xyyyz_1, g_0_xyyzz_0, g_0_xyyzz_1, g_0_xyzzz_0, g_0_xyzzz_1, g_0_xzzzz_0, g_0_xzzzz_1, g_0_yyyyy_0, g_0_yyyyy_1, g_0_yyyyz_0, g_0_yyyyz_1, g_0_yyyzz_0, g_0_yyyzz_1, g_0_yyzzz_0, g_0_yyzzz_1, g_0_yzzzz_0, g_0_yzzzz_1, g_0_zzzzz_0, g_0_zzzzz_1, g_x_xxxx_1, g_x_xxxxx_1, g_x_xxxxy_1, g_x_xxxxz_1, g_x_xxxy_1, g_x_xxxyy_1, g_x_xxxyz_1, g_x_xxxz_1, g_x_xxxzz_1, g_x_xxyy_1, g_x_xxyyy_1, g_x_xxyyz_1, g_x_xxyz_1, g_x_xxyzz_1, g_x_xxzz_1, g_x_xxzzz_1, g_x_xyyy_1, g_x_xyyyy_1, g_x_xyyyz_1, g_x_xyyz_1, g_x_xyyzz_1, g_x_xyzz_1, g_x_xyzzz_1, g_x_xzzz_1, g_x_xzzzz_1, g_x_yyyy_1, g_x_yyyyy_1, g_x_yyyyz_1, g_x_yyyz_1, g_x_yyyzz_1, g_x_yyzz_1, g_x_yyzzz_1, g_x_yzzz_1, g_x_yzzzz_1, g_x_zzzz_1, g_x_zzzzz_1, g_xx_xxxxx_0, g_xx_xxxxy_0, g_xx_xxxxz_0, g_xx_xxxyy_0, g_xx_xxxyz_0, g_xx_xxxzz_0, g_xx_xxyyy_0, g_xx_xxyyz_0, g_xx_xxyzz_0, g_xx_xxzzz_0, g_xx_xyyyy_0, g_xx_xyyyz_0, g_xx_xyyzz_0, g_xx_xyzzz_0, g_xx_xzzzz_0, g_xx_yyyyy_0, g_xx_yyyyz_0, g_xx_yyyzz_0, g_xx_yyzzz_0, g_xx_yzzzz_0, g_xx_zzzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xx_xxxxx_0[i] = g_0_xxxxx_0[i] * fbe_0 - g_0_xxxxx_1[i] * fz_be_0 + 5.0 * g_x_xxxx_1[i] * fe_0 + g_x_xxxxx_1[i] * pa_x[i];

        g_xx_xxxxy_0[i] = g_0_xxxxy_0[i] * fbe_0 - g_0_xxxxy_1[i] * fz_be_0 + 4.0 * g_x_xxxy_1[i] * fe_0 + g_x_xxxxy_1[i] * pa_x[i];

        g_xx_xxxxz_0[i] = g_0_xxxxz_0[i] * fbe_0 - g_0_xxxxz_1[i] * fz_be_0 + 4.0 * g_x_xxxz_1[i] * fe_0 + g_x_xxxxz_1[i] * pa_x[i];

        g_xx_xxxyy_0[i] = g_0_xxxyy_0[i] * fbe_0 - g_0_xxxyy_1[i] * fz_be_0 + 3.0 * g_x_xxyy_1[i] * fe_0 + g_x_xxxyy_1[i] * pa_x[i];

        g_xx_xxxyz_0[i] = g_0_xxxyz_0[i] * fbe_0 - g_0_xxxyz_1[i] * fz_be_0 + 3.0 * g_x_xxyz_1[i] * fe_0 + g_x_xxxyz_1[i] * pa_x[i];

        g_xx_xxxzz_0[i] = g_0_xxxzz_0[i] * fbe_0 - g_0_xxxzz_1[i] * fz_be_0 + 3.0 * g_x_xxzz_1[i] * fe_0 + g_x_xxxzz_1[i] * pa_x[i];

        g_xx_xxyyy_0[i] = g_0_xxyyy_0[i] * fbe_0 - g_0_xxyyy_1[i] * fz_be_0 + 2.0 * g_x_xyyy_1[i] * fe_0 + g_x_xxyyy_1[i] * pa_x[i];

        g_xx_xxyyz_0[i] = g_0_xxyyz_0[i] * fbe_0 - g_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_x_xyyz_1[i] * fe_0 + g_x_xxyyz_1[i] * pa_x[i];

        g_xx_xxyzz_0[i] = g_0_xxyzz_0[i] * fbe_0 - g_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_x_xyzz_1[i] * fe_0 + g_x_xxyzz_1[i] * pa_x[i];

        g_xx_xxzzz_0[i] = g_0_xxzzz_0[i] * fbe_0 - g_0_xxzzz_1[i] * fz_be_0 + 2.0 * g_x_xzzz_1[i] * fe_0 + g_x_xxzzz_1[i] * pa_x[i];

        g_xx_xyyyy_0[i] = g_0_xyyyy_0[i] * fbe_0 - g_0_xyyyy_1[i] * fz_be_0 + g_x_yyyy_1[i] * fe_0 + g_x_xyyyy_1[i] * pa_x[i];

        g_xx_xyyyz_0[i] = g_0_xyyyz_0[i] * fbe_0 - g_0_xyyyz_1[i] * fz_be_0 + g_x_yyyz_1[i] * fe_0 + g_x_xyyyz_1[i] * pa_x[i];

        g_xx_xyyzz_0[i] = g_0_xyyzz_0[i] * fbe_0 - g_0_xyyzz_1[i] * fz_be_0 + g_x_yyzz_1[i] * fe_0 + g_x_xyyzz_1[i] * pa_x[i];

        g_xx_xyzzz_0[i] = g_0_xyzzz_0[i] * fbe_0 - g_0_xyzzz_1[i] * fz_be_0 + g_x_yzzz_1[i] * fe_0 + g_x_xyzzz_1[i] * pa_x[i];

        g_xx_xzzzz_0[i] = g_0_xzzzz_0[i] * fbe_0 - g_0_xzzzz_1[i] * fz_be_0 + g_x_zzzz_1[i] * fe_0 + g_x_xzzzz_1[i] * pa_x[i];

        g_xx_yyyyy_0[i] = g_0_yyyyy_0[i] * fbe_0 - g_0_yyyyy_1[i] * fz_be_0 + g_x_yyyyy_1[i] * pa_x[i];

        g_xx_yyyyz_0[i] = g_0_yyyyz_0[i] * fbe_0 - g_0_yyyyz_1[i] * fz_be_0 + g_x_yyyyz_1[i] * pa_x[i];

        g_xx_yyyzz_0[i] = g_0_yyyzz_0[i] * fbe_0 - g_0_yyyzz_1[i] * fz_be_0 + g_x_yyyzz_1[i] * pa_x[i];

        g_xx_yyzzz_0[i] = g_0_yyzzz_0[i] * fbe_0 - g_0_yyzzz_1[i] * fz_be_0 + g_x_yyzzz_1[i] * pa_x[i];

        g_xx_yzzzz_0[i] = g_0_yzzzz_0[i] * fbe_0 - g_0_yzzzz_1[i] * fz_be_0 + g_x_yzzzz_1[i] * pa_x[i];

        g_xx_zzzzz_0[i] = g_0_zzzzz_0[i] * fbe_0 - g_0_zzzzz_1[i] * fz_be_0 + g_x_zzzzz_1[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : DH

    auto g_xy_xxxxx_0 = pbuffer.data(idx_eri_0_dh + 21);

    auto g_xy_xxxxy_0 = pbuffer.data(idx_eri_0_dh + 22);

    auto g_xy_xxxxz_0 = pbuffer.data(idx_eri_0_dh + 23);

    auto g_xy_xxxyy_0 = pbuffer.data(idx_eri_0_dh + 24);

    auto g_xy_xxxyz_0 = pbuffer.data(idx_eri_0_dh + 25);

    auto g_xy_xxxzz_0 = pbuffer.data(idx_eri_0_dh + 26);

    auto g_xy_xxyyy_0 = pbuffer.data(idx_eri_0_dh + 27);

    auto g_xy_xxyyz_0 = pbuffer.data(idx_eri_0_dh + 28);

    auto g_xy_xxyzz_0 = pbuffer.data(idx_eri_0_dh + 29);

    auto g_xy_xxzzz_0 = pbuffer.data(idx_eri_0_dh + 30);

    auto g_xy_xyyyy_0 = pbuffer.data(idx_eri_0_dh + 31);

    auto g_xy_xyyyz_0 = pbuffer.data(idx_eri_0_dh + 32);

    auto g_xy_xyyzz_0 = pbuffer.data(idx_eri_0_dh + 33);

    auto g_xy_xyzzz_0 = pbuffer.data(idx_eri_0_dh + 34);

    auto g_xy_xzzzz_0 = pbuffer.data(idx_eri_0_dh + 35);

    auto g_xy_yyyyy_0 = pbuffer.data(idx_eri_0_dh + 36);

    auto g_xy_yyyyz_0 = pbuffer.data(idx_eri_0_dh + 37);

    auto g_xy_yyyzz_0 = pbuffer.data(idx_eri_0_dh + 38);

    auto g_xy_yyzzz_0 = pbuffer.data(idx_eri_0_dh + 39);

    auto g_xy_yzzzz_0 = pbuffer.data(idx_eri_0_dh + 40);

    auto g_xy_zzzzz_0 = pbuffer.data(idx_eri_0_dh + 41);

    #pragma omp simd aligned(g_x_xxxxx_1, g_x_xxxxz_1, g_x_xxxzz_1, g_x_xxzzz_1, g_x_xzzzz_1, g_xy_xxxxx_0, g_xy_xxxxy_0, g_xy_xxxxz_0, g_xy_xxxyy_0, g_xy_xxxyz_0, g_xy_xxxzz_0, g_xy_xxyyy_0, g_xy_xxyyz_0, g_xy_xxyzz_0, g_xy_xxzzz_0, g_xy_xyyyy_0, g_xy_xyyyz_0, g_xy_xyyzz_0, g_xy_xyzzz_0, g_xy_xzzzz_0, g_xy_yyyyy_0, g_xy_yyyyz_0, g_xy_yyyzz_0, g_xy_yyzzz_0, g_xy_yzzzz_0, g_xy_zzzzz_0, g_y_xxxxy_1, g_y_xxxy_1, g_y_xxxyy_1, g_y_xxxyz_1, g_y_xxyy_1, g_y_xxyyy_1, g_y_xxyyz_1, g_y_xxyz_1, g_y_xxyzz_1, g_y_xyyy_1, g_y_xyyyy_1, g_y_xyyyz_1, g_y_xyyz_1, g_y_xyyzz_1, g_y_xyzz_1, g_y_xyzzz_1, g_y_yyyy_1, g_y_yyyyy_1, g_y_yyyyz_1, g_y_yyyz_1, g_y_yyyzz_1, g_y_yyzz_1, g_y_yyzzz_1, g_y_yzzz_1, g_y_yzzzz_1, g_y_zzzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xy_xxxxx_0[i] = g_x_xxxxx_1[i] * pa_y[i];

        g_xy_xxxxy_0[i] = 4.0 * g_y_xxxy_1[i] * fe_0 + g_y_xxxxy_1[i] * pa_x[i];

        g_xy_xxxxz_0[i] = g_x_xxxxz_1[i] * pa_y[i];

        g_xy_xxxyy_0[i] = 3.0 * g_y_xxyy_1[i] * fe_0 + g_y_xxxyy_1[i] * pa_x[i];

        g_xy_xxxyz_0[i] = 3.0 * g_y_xxyz_1[i] * fe_0 + g_y_xxxyz_1[i] * pa_x[i];

        g_xy_xxxzz_0[i] = g_x_xxxzz_1[i] * pa_y[i];

        g_xy_xxyyy_0[i] = 2.0 * g_y_xyyy_1[i] * fe_0 + g_y_xxyyy_1[i] * pa_x[i];

        g_xy_xxyyz_0[i] = 2.0 * g_y_xyyz_1[i] * fe_0 + g_y_xxyyz_1[i] * pa_x[i];

        g_xy_xxyzz_0[i] = 2.0 * g_y_xyzz_1[i] * fe_0 + g_y_xxyzz_1[i] * pa_x[i];

        g_xy_xxzzz_0[i] = g_x_xxzzz_1[i] * pa_y[i];

        g_xy_xyyyy_0[i] = g_y_yyyy_1[i] * fe_0 + g_y_xyyyy_1[i] * pa_x[i];

        g_xy_xyyyz_0[i] = g_y_yyyz_1[i] * fe_0 + g_y_xyyyz_1[i] * pa_x[i];

        g_xy_xyyzz_0[i] = g_y_yyzz_1[i] * fe_0 + g_y_xyyzz_1[i] * pa_x[i];

        g_xy_xyzzz_0[i] = g_y_yzzz_1[i] * fe_0 + g_y_xyzzz_1[i] * pa_x[i];

        g_xy_xzzzz_0[i] = g_x_xzzzz_1[i] * pa_y[i];

        g_xy_yyyyy_0[i] = g_y_yyyyy_1[i] * pa_x[i];

        g_xy_yyyyz_0[i] = g_y_yyyyz_1[i] * pa_x[i];

        g_xy_yyyzz_0[i] = g_y_yyyzz_1[i] * pa_x[i];

        g_xy_yyzzz_0[i] = g_y_yyzzz_1[i] * pa_x[i];

        g_xy_yzzzz_0[i] = g_y_yzzzz_1[i] * pa_x[i];

        g_xy_zzzzz_0[i] = g_y_zzzzz_1[i] * pa_x[i];
    }

    // Set up 42-63 components of targeted buffer : DH

    auto g_xz_xxxxx_0 = pbuffer.data(idx_eri_0_dh + 42);

    auto g_xz_xxxxy_0 = pbuffer.data(idx_eri_0_dh + 43);

    auto g_xz_xxxxz_0 = pbuffer.data(idx_eri_0_dh + 44);

    auto g_xz_xxxyy_0 = pbuffer.data(idx_eri_0_dh + 45);

    auto g_xz_xxxyz_0 = pbuffer.data(idx_eri_0_dh + 46);

    auto g_xz_xxxzz_0 = pbuffer.data(idx_eri_0_dh + 47);

    auto g_xz_xxyyy_0 = pbuffer.data(idx_eri_0_dh + 48);

    auto g_xz_xxyyz_0 = pbuffer.data(idx_eri_0_dh + 49);

    auto g_xz_xxyzz_0 = pbuffer.data(idx_eri_0_dh + 50);

    auto g_xz_xxzzz_0 = pbuffer.data(idx_eri_0_dh + 51);

    auto g_xz_xyyyy_0 = pbuffer.data(idx_eri_0_dh + 52);

    auto g_xz_xyyyz_0 = pbuffer.data(idx_eri_0_dh + 53);

    auto g_xz_xyyzz_0 = pbuffer.data(idx_eri_0_dh + 54);

    auto g_xz_xyzzz_0 = pbuffer.data(idx_eri_0_dh + 55);

    auto g_xz_xzzzz_0 = pbuffer.data(idx_eri_0_dh + 56);

    auto g_xz_yyyyy_0 = pbuffer.data(idx_eri_0_dh + 57);

    auto g_xz_yyyyz_0 = pbuffer.data(idx_eri_0_dh + 58);

    auto g_xz_yyyzz_0 = pbuffer.data(idx_eri_0_dh + 59);

    auto g_xz_yyzzz_0 = pbuffer.data(idx_eri_0_dh + 60);

    auto g_xz_yzzzz_0 = pbuffer.data(idx_eri_0_dh + 61);

    auto g_xz_zzzzz_0 = pbuffer.data(idx_eri_0_dh + 62);

    #pragma omp simd aligned(g_x_xxxxx_1, g_x_xxxxy_1, g_x_xxxyy_1, g_x_xxyyy_1, g_x_xyyyy_1, g_xz_xxxxx_0, g_xz_xxxxy_0, g_xz_xxxxz_0, g_xz_xxxyy_0, g_xz_xxxyz_0, g_xz_xxxzz_0, g_xz_xxyyy_0, g_xz_xxyyz_0, g_xz_xxyzz_0, g_xz_xxzzz_0, g_xz_xyyyy_0, g_xz_xyyyz_0, g_xz_xyyzz_0, g_xz_xyzzz_0, g_xz_xzzzz_0, g_xz_yyyyy_0, g_xz_yyyyz_0, g_xz_yyyzz_0, g_xz_yyzzz_0, g_xz_yzzzz_0, g_xz_zzzzz_0, g_z_xxxxz_1, g_z_xxxyz_1, g_z_xxxz_1, g_z_xxxzz_1, g_z_xxyyz_1, g_z_xxyz_1, g_z_xxyzz_1, g_z_xxzz_1, g_z_xxzzz_1, g_z_xyyyz_1, g_z_xyyz_1, g_z_xyyzz_1, g_z_xyzz_1, g_z_xyzzz_1, g_z_xzzz_1, g_z_xzzzz_1, g_z_yyyyy_1, g_z_yyyyz_1, g_z_yyyz_1, g_z_yyyzz_1, g_z_yyzz_1, g_z_yyzzz_1, g_z_yzzz_1, g_z_yzzzz_1, g_z_zzzz_1, g_z_zzzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xz_xxxxx_0[i] = g_x_xxxxx_1[i] * pa_z[i];

        g_xz_xxxxy_0[i] = g_x_xxxxy_1[i] * pa_z[i];

        g_xz_xxxxz_0[i] = 4.0 * g_z_xxxz_1[i] * fe_0 + g_z_xxxxz_1[i] * pa_x[i];

        g_xz_xxxyy_0[i] = g_x_xxxyy_1[i] * pa_z[i];

        g_xz_xxxyz_0[i] = 3.0 * g_z_xxyz_1[i] * fe_0 + g_z_xxxyz_1[i] * pa_x[i];

        g_xz_xxxzz_0[i] = 3.0 * g_z_xxzz_1[i] * fe_0 + g_z_xxxzz_1[i] * pa_x[i];

        g_xz_xxyyy_0[i] = g_x_xxyyy_1[i] * pa_z[i];

        g_xz_xxyyz_0[i] = 2.0 * g_z_xyyz_1[i] * fe_0 + g_z_xxyyz_1[i] * pa_x[i];

        g_xz_xxyzz_0[i] = 2.0 * g_z_xyzz_1[i] * fe_0 + g_z_xxyzz_1[i] * pa_x[i];

        g_xz_xxzzz_0[i] = 2.0 * g_z_xzzz_1[i] * fe_0 + g_z_xxzzz_1[i] * pa_x[i];

        g_xz_xyyyy_0[i] = g_x_xyyyy_1[i] * pa_z[i];

        g_xz_xyyyz_0[i] = g_z_yyyz_1[i] * fe_0 + g_z_xyyyz_1[i] * pa_x[i];

        g_xz_xyyzz_0[i] = g_z_yyzz_1[i] * fe_0 + g_z_xyyzz_1[i] * pa_x[i];

        g_xz_xyzzz_0[i] = g_z_yzzz_1[i] * fe_0 + g_z_xyzzz_1[i] * pa_x[i];

        g_xz_xzzzz_0[i] = g_z_zzzz_1[i] * fe_0 + g_z_xzzzz_1[i] * pa_x[i];

        g_xz_yyyyy_0[i] = g_z_yyyyy_1[i] * pa_x[i];

        g_xz_yyyyz_0[i] = g_z_yyyyz_1[i] * pa_x[i];

        g_xz_yyyzz_0[i] = g_z_yyyzz_1[i] * pa_x[i];

        g_xz_yyzzz_0[i] = g_z_yyzzz_1[i] * pa_x[i];

        g_xz_yzzzz_0[i] = g_z_yzzzz_1[i] * pa_x[i];

        g_xz_zzzzz_0[i] = g_z_zzzzz_1[i] * pa_x[i];
    }

    // Set up 63-84 components of targeted buffer : DH

    auto g_yy_xxxxx_0 = pbuffer.data(idx_eri_0_dh + 63);

    auto g_yy_xxxxy_0 = pbuffer.data(idx_eri_0_dh + 64);

    auto g_yy_xxxxz_0 = pbuffer.data(idx_eri_0_dh + 65);

    auto g_yy_xxxyy_0 = pbuffer.data(idx_eri_0_dh + 66);

    auto g_yy_xxxyz_0 = pbuffer.data(idx_eri_0_dh + 67);

    auto g_yy_xxxzz_0 = pbuffer.data(idx_eri_0_dh + 68);

    auto g_yy_xxyyy_0 = pbuffer.data(idx_eri_0_dh + 69);

    auto g_yy_xxyyz_0 = pbuffer.data(idx_eri_0_dh + 70);

    auto g_yy_xxyzz_0 = pbuffer.data(idx_eri_0_dh + 71);

    auto g_yy_xxzzz_0 = pbuffer.data(idx_eri_0_dh + 72);

    auto g_yy_xyyyy_0 = pbuffer.data(idx_eri_0_dh + 73);

    auto g_yy_xyyyz_0 = pbuffer.data(idx_eri_0_dh + 74);

    auto g_yy_xyyzz_0 = pbuffer.data(idx_eri_0_dh + 75);

    auto g_yy_xyzzz_0 = pbuffer.data(idx_eri_0_dh + 76);

    auto g_yy_xzzzz_0 = pbuffer.data(idx_eri_0_dh + 77);

    auto g_yy_yyyyy_0 = pbuffer.data(idx_eri_0_dh + 78);

    auto g_yy_yyyyz_0 = pbuffer.data(idx_eri_0_dh + 79);

    auto g_yy_yyyzz_0 = pbuffer.data(idx_eri_0_dh + 80);

    auto g_yy_yyzzz_0 = pbuffer.data(idx_eri_0_dh + 81);

    auto g_yy_yzzzz_0 = pbuffer.data(idx_eri_0_dh + 82);

    auto g_yy_zzzzz_0 = pbuffer.data(idx_eri_0_dh + 83);

    #pragma omp simd aligned(g_0_xxxxx_0, g_0_xxxxx_1, g_0_xxxxy_0, g_0_xxxxy_1, g_0_xxxxz_0, g_0_xxxxz_1, g_0_xxxyy_0, g_0_xxxyy_1, g_0_xxxyz_0, g_0_xxxyz_1, g_0_xxxzz_0, g_0_xxxzz_1, g_0_xxyyy_0, g_0_xxyyy_1, g_0_xxyyz_0, g_0_xxyyz_1, g_0_xxyzz_0, g_0_xxyzz_1, g_0_xxzzz_0, g_0_xxzzz_1, g_0_xyyyy_0, g_0_xyyyy_1, g_0_xyyyz_0, g_0_xyyyz_1, g_0_xyyzz_0, g_0_xyyzz_1, g_0_xyzzz_0, g_0_xyzzz_1, g_0_xzzzz_0, g_0_xzzzz_1, g_0_yyyyy_0, g_0_yyyyy_1, g_0_yyyyz_0, g_0_yyyyz_1, g_0_yyyzz_0, g_0_yyyzz_1, g_0_yyzzz_0, g_0_yyzzz_1, g_0_yzzzz_0, g_0_yzzzz_1, g_0_zzzzz_0, g_0_zzzzz_1, g_y_xxxx_1, g_y_xxxxx_1, g_y_xxxxy_1, g_y_xxxxz_1, g_y_xxxy_1, g_y_xxxyy_1, g_y_xxxyz_1, g_y_xxxz_1, g_y_xxxzz_1, g_y_xxyy_1, g_y_xxyyy_1, g_y_xxyyz_1, g_y_xxyz_1, g_y_xxyzz_1, g_y_xxzz_1, g_y_xxzzz_1, g_y_xyyy_1, g_y_xyyyy_1, g_y_xyyyz_1, g_y_xyyz_1, g_y_xyyzz_1, g_y_xyzz_1, g_y_xyzzz_1, g_y_xzzz_1, g_y_xzzzz_1, g_y_yyyy_1, g_y_yyyyy_1, g_y_yyyyz_1, g_y_yyyz_1, g_y_yyyzz_1, g_y_yyzz_1, g_y_yyzzz_1, g_y_yzzz_1, g_y_yzzzz_1, g_y_zzzz_1, g_y_zzzzz_1, g_yy_xxxxx_0, g_yy_xxxxy_0, g_yy_xxxxz_0, g_yy_xxxyy_0, g_yy_xxxyz_0, g_yy_xxxzz_0, g_yy_xxyyy_0, g_yy_xxyyz_0, g_yy_xxyzz_0, g_yy_xxzzz_0, g_yy_xyyyy_0, g_yy_xyyyz_0, g_yy_xyyzz_0, g_yy_xyzzz_0, g_yy_xzzzz_0, g_yy_yyyyy_0, g_yy_yyyyz_0, g_yy_yyyzz_0, g_yy_yyzzz_0, g_yy_yzzzz_0, g_yy_zzzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yy_xxxxx_0[i] = g_0_xxxxx_0[i] * fbe_0 - g_0_xxxxx_1[i] * fz_be_0 + g_y_xxxxx_1[i] * pa_y[i];

        g_yy_xxxxy_0[i] = g_0_xxxxy_0[i] * fbe_0 - g_0_xxxxy_1[i] * fz_be_0 + g_y_xxxx_1[i] * fe_0 + g_y_xxxxy_1[i] * pa_y[i];

        g_yy_xxxxz_0[i] = g_0_xxxxz_0[i] * fbe_0 - g_0_xxxxz_1[i] * fz_be_0 + g_y_xxxxz_1[i] * pa_y[i];

        g_yy_xxxyy_0[i] = g_0_xxxyy_0[i] * fbe_0 - g_0_xxxyy_1[i] * fz_be_0 + 2.0 * g_y_xxxy_1[i] * fe_0 + g_y_xxxyy_1[i] * pa_y[i];

        g_yy_xxxyz_0[i] = g_0_xxxyz_0[i] * fbe_0 - g_0_xxxyz_1[i] * fz_be_0 + g_y_xxxz_1[i] * fe_0 + g_y_xxxyz_1[i] * pa_y[i];

        g_yy_xxxzz_0[i] = g_0_xxxzz_0[i] * fbe_0 - g_0_xxxzz_1[i] * fz_be_0 + g_y_xxxzz_1[i] * pa_y[i];

        g_yy_xxyyy_0[i] = g_0_xxyyy_0[i] * fbe_0 - g_0_xxyyy_1[i] * fz_be_0 + 3.0 * g_y_xxyy_1[i] * fe_0 + g_y_xxyyy_1[i] * pa_y[i];

        g_yy_xxyyz_0[i] = g_0_xxyyz_0[i] * fbe_0 - g_0_xxyyz_1[i] * fz_be_0 + 2.0 * g_y_xxyz_1[i] * fe_0 + g_y_xxyyz_1[i] * pa_y[i];

        g_yy_xxyzz_0[i] = g_0_xxyzz_0[i] * fbe_0 - g_0_xxyzz_1[i] * fz_be_0 + g_y_xxzz_1[i] * fe_0 + g_y_xxyzz_1[i] * pa_y[i];

        g_yy_xxzzz_0[i] = g_0_xxzzz_0[i] * fbe_0 - g_0_xxzzz_1[i] * fz_be_0 + g_y_xxzzz_1[i] * pa_y[i];

        g_yy_xyyyy_0[i] = g_0_xyyyy_0[i] * fbe_0 - g_0_xyyyy_1[i] * fz_be_0 + 4.0 * g_y_xyyy_1[i] * fe_0 + g_y_xyyyy_1[i] * pa_y[i];

        g_yy_xyyyz_0[i] = g_0_xyyyz_0[i] * fbe_0 - g_0_xyyyz_1[i] * fz_be_0 + 3.0 * g_y_xyyz_1[i] * fe_0 + g_y_xyyyz_1[i] * pa_y[i];

        g_yy_xyyzz_0[i] = g_0_xyyzz_0[i] * fbe_0 - g_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_y_xyzz_1[i] * fe_0 + g_y_xyyzz_1[i] * pa_y[i];

        g_yy_xyzzz_0[i] = g_0_xyzzz_0[i] * fbe_0 - g_0_xyzzz_1[i] * fz_be_0 + g_y_xzzz_1[i] * fe_0 + g_y_xyzzz_1[i] * pa_y[i];

        g_yy_xzzzz_0[i] = g_0_xzzzz_0[i] * fbe_0 - g_0_xzzzz_1[i] * fz_be_0 + g_y_xzzzz_1[i] * pa_y[i];

        g_yy_yyyyy_0[i] = g_0_yyyyy_0[i] * fbe_0 - g_0_yyyyy_1[i] * fz_be_0 + 5.0 * g_y_yyyy_1[i] * fe_0 + g_y_yyyyy_1[i] * pa_y[i];

        g_yy_yyyyz_0[i] = g_0_yyyyz_0[i] * fbe_0 - g_0_yyyyz_1[i] * fz_be_0 + 4.0 * g_y_yyyz_1[i] * fe_0 + g_y_yyyyz_1[i] * pa_y[i];

        g_yy_yyyzz_0[i] = g_0_yyyzz_0[i] * fbe_0 - g_0_yyyzz_1[i] * fz_be_0 + 3.0 * g_y_yyzz_1[i] * fe_0 + g_y_yyyzz_1[i] * pa_y[i];

        g_yy_yyzzz_0[i] = g_0_yyzzz_0[i] * fbe_0 - g_0_yyzzz_1[i] * fz_be_0 + 2.0 * g_y_yzzz_1[i] * fe_0 + g_y_yyzzz_1[i] * pa_y[i];

        g_yy_yzzzz_0[i] = g_0_yzzzz_0[i] * fbe_0 - g_0_yzzzz_1[i] * fz_be_0 + g_y_zzzz_1[i] * fe_0 + g_y_yzzzz_1[i] * pa_y[i];

        g_yy_zzzzz_0[i] = g_0_zzzzz_0[i] * fbe_0 - g_0_zzzzz_1[i] * fz_be_0 + g_y_zzzzz_1[i] * pa_y[i];
    }

    // Set up 84-105 components of targeted buffer : DH

    auto g_yz_xxxxx_0 = pbuffer.data(idx_eri_0_dh + 84);

    auto g_yz_xxxxy_0 = pbuffer.data(idx_eri_0_dh + 85);

    auto g_yz_xxxxz_0 = pbuffer.data(idx_eri_0_dh + 86);

    auto g_yz_xxxyy_0 = pbuffer.data(idx_eri_0_dh + 87);

    auto g_yz_xxxyz_0 = pbuffer.data(idx_eri_0_dh + 88);

    auto g_yz_xxxzz_0 = pbuffer.data(idx_eri_0_dh + 89);

    auto g_yz_xxyyy_0 = pbuffer.data(idx_eri_0_dh + 90);

    auto g_yz_xxyyz_0 = pbuffer.data(idx_eri_0_dh + 91);

    auto g_yz_xxyzz_0 = pbuffer.data(idx_eri_0_dh + 92);

    auto g_yz_xxzzz_0 = pbuffer.data(idx_eri_0_dh + 93);

    auto g_yz_xyyyy_0 = pbuffer.data(idx_eri_0_dh + 94);

    auto g_yz_xyyyz_0 = pbuffer.data(idx_eri_0_dh + 95);

    auto g_yz_xyyzz_0 = pbuffer.data(idx_eri_0_dh + 96);

    auto g_yz_xyzzz_0 = pbuffer.data(idx_eri_0_dh + 97);

    auto g_yz_xzzzz_0 = pbuffer.data(idx_eri_0_dh + 98);

    auto g_yz_yyyyy_0 = pbuffer.data(idx_eri_0_dh + 99);

    auto g_yz_yyyyz_0 = pbuffer.data(idx_eri_0_dh + 100);

    auto g_yz_yyyzz_0 = pbuffer.data(idx_eri_0_dh + 101);

    auto g_yz_yyzzz_0 = pbuffer.data(idx_eri_0_dh + 102);

    auto g_yz_yzzzz_0 = pbuffer.data(idx_eri_0_dh + 103);

    auto g_yz_zzzzz_0 = pbuffer.data(idx_eri_0_dh + 104);

    #pragma omp simd aligned(g_y_xxxxy_1, g_y_xxxyy_1, g_y_xxyyy_1, g_y_xyyyy_1, g_y_yyyyy_1, g_yz_xxxxx_0, g_yz_xxxxy_0, g_yz_xxxxz_0, g_yz_xxxyy_0, g_yz_xxxyz_0, g_yz_xxxzz_0, g_yz_xxyyy_0, g_yz_xxyyz_0, g_yz_xxyzz_0, g_yz_xxzzz_0, g_yz_xyyyy_0, g_yz_xyyyz_0, g_yz_xyyzz_0, g_yz_xyzzz_0, g_yz_xzzzz_0, g_yz_yyyyy_0, g_yz_yyyyz_0, g_yz_yyyzz_0, g_yz_yyzzz_0, g_yz_yzzzz_0, g_yz_zzzzz_0, g_z_xxxxx_1, g_z_xxxxz_1, g_z_xxxyz_1, g_z_xxxz_1, g_z_xxxzz_1, g_z_xxyyz_1, g_z_xxyz_1, g_z_xxyzz_1, g_z_xxzz_1, g_z_xxzzz_1, g_z_xyyyz_1, g_z_xyyz_1, g_z_xyyzz_1, g_z_xyzz_1, g_z_xyzzz_1, g_z_xzzz_1, g_z_xzzzz_1, g_z_yyyyz_1, g_z_yyyz_1, g_z_yyyzz_1, g_z_yyzz_1, g_z_yyzzz_1, g_z_yzzz_1, g_z_yzzzz_1, g_z_zzzz_1, g_z_zzzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yz_xxxxx_0[i] = g_z_xxxxx_1[i] * pa_y[i];

        g_yz_xxxxy_0[i] = g_y_xxxxy_1[i] * pa_z[i];

        g_yz_xxxxz_0[i] = g_z_xxxxz_1[i] * pa_y[i];

        g_yz_xxxyy_0[i] = g_y_xxxyy_1[i] * pa_z[i];

        g_yz_xxxyz_0[i] = g_z_xxxz_1[i] * fe_0 + g_z_xxxyz_1[i] * pa_y[i];

        g_yz_xxxzz_0[i] = g_z_xxxzz_1[i] * pa_y[i];

        g_yz_xxyyy_0[i] = g_y_xxyyy_1[i] * pa_z[i];

        g_yz_xxyyz_0[i] = 2.0 * g_z_xxyz_1[i] * fe_0 + g_z_xxyyz_1[i] * pa_y[i];

        g_yz_xxyzz_0[i] = g_z_xxzz_1[i] * fe_0 + g_z_xxyzz_1[i] * pa_y[i];

        g_yz_xxzzz_0[i] = g_z_xxzzz_1[i] * pa_y[i];

        g_yz_xyyyy_0[i] = g_y_xyyyy_1[i] * pa_z[i];

        g_yz_xyyyz_0[i] = 3.0 * g_z_xyyz_1[i] * fe_0 + g_z_xyyyz_1[i] * pa_y[i];

        g_yz_xyyzz_0[i] = 2.0 * g_z_xyzz_1[i] * fe_0 + g_z_xyyzz_1[i] * pa_y[i];

        g_yz_xyzzz_0[i] = g_z_xzzz_1[i] * fe_0 + g_z_xyzzz_1[i] * pa_y[i];

        g_yz_xzzzz_0[i] = g_z_xzzzz_1[i] * pa_y[i];

        g_yz_yyyyy_0[i] = g_y_yyyyy_1[i] * pa_z[i];

        g_yz_yyyyz_0[i] = 4.0 * g_z_yyyz_1[i] * fe_0 + g_z_yyyyz_1[i] * pa_y[i];

        g_yz_yyyzz_0[i] = 3.0 * g_z_yyzz_1[i] * fe_0 + g_z_yyyzz_1[i] * pa_y[i];

        g_yz_yyzzz_0[i] = 2.0 * g_z_yzzz_1[i] * fe_0 + g_z_yyzzz_1[i] * pa_y[i];

        g_yz_yzzzz_0[i] = g_z_zzzz_1[i] * fe_0 + g_z_yzzzz_1[i] * pa_y[i];

        g_yz_zzzzz_0[i] = g_z_zzzzz_1[i] * pa_y[i];
    }

    // Set up 105-126 components of targeted buffer : DH

    auto g_zz_xxxxx_0 = pbuffer.data(idx_eri_0_dh + 105);

    auto g_zz_xxxxy_0 = pbuffer.data(idx_eri_0_dh + 106);

    auto g_zz_xxxxz_0 = pbuffer.data(idx_eri_0_dh + 107);

    auto g_zz_xxxyy_0 = pbuffer.data(idx_eri_0_dh + 108);

    auto g_zz_xxxyz_0 = pbuffer.data(idx_eri_0_dh + 109);

    auto g_zz_xxxzz_0 = pbuffer.data(idx_eri_0_dh + 110);

    auto g_zz_xxyyy_0 = pbuffer.data(idx_eri_0_dh + 111);

    auto g_zz_xxyyz_0 = pbuffer.data(idx_eri_0_dh + 112);

    auto g_zz_xxyzz_0 = pbuffer.data(idx_eri_0_dh + 113);

    auto g_zz_xxzzz_0 = pbuffer.data(idx_eri_0_dh + 114);

    auto g_zz_xyyyy_0 = pbuffer.data(idx_eri_0_dh + 115);

    auto g_zz_xyyyz_0 = pbuffer.data(idx_eri_0_dh + 116);

    auto g_zz_xyyzz_0 = pbuffer.data(idx_eri_0_dh + 117);

    auto g_zz_xyzzz_0 = pbuffer.data(idx_eri_0_dh + 118);

    auto g_zz_xzzzz_0 = pbuffer.data(idx_eri_0_dh + 119);

    auto g_zz_yyyyy_0 = pbuffer.data(idx_eri_0_dh + 120);

    auto g_zz_yyyyz_0 = pbuffer.data(idx_eri_0_dh + 121);

    auto g_zz_yyyzz_0 = pbuffer.data(idx_eri_0_dh + 122);

    auto g_zz_yyzzz_0 = pbuffer.data(idx_eri_0_dh + 123);

    auto g_zz_yzzzz_0 = pbuffer.data(idx_eri_0_dh + 124);

    auto g_zz_zzzzz_0 = pbuffer.data(idx_eri_0_dh + 125);

    #pragma omp simd aligned(g_0_xxxxx_0, g_0_xxxxx_1, g_0_xxxxy_0, g_0_xxxxy_1, g_0_xxxxz_0, g_0_xxxxz_1, g_0_xxxyy_0, g_0_xxxyy_1, g_0_xxxyz_0, g_0_xxxyz_1, g_0_xxxzz_0, g_0_xxxzz_1, g_0_xxyyy_0, g_0_xxyyy_1, g_0_xxyyz_0, g_0_xxyyz_1, g_0_xxyzz_0, g_0_xxyzz_1, g_0_xxzzz_0, g_0_xxzzz_1, g_0_xyyyy_0, g_0_xyyyy_1, g_0_xyyyz_0, g_0_xyyyz_1, g_0_xyyzz_0, g_0_xyyzz_1, g_0_xyzzz_0, g_0_xyzzz_1, g_0_xzzzz_0, g_0_xzzzz_1, g_0_yyyyy_0, g_0_yyyyy_1, g_0_yyyyz_0, g_0_yyyyz_1, g_0_yyyzz_0, g_0_yyyzz_1, g_0_yyzzz_0, g_0_yyzzz_1, g_0_yzzzz_0, g_0_yzzzz_1, g_0_zzzzz_0, g_0_zzzzz_1, g_z_xxxx_1, g_z_xxxxx_1, g_z_xxxxy_1, g_z_xxxxz_1, g_z_xxxy_1, g_z_xxxyy_1, g_z_xxxyz_1, g_z_xxxz_1, g_z_xxxzz_1, g_z_xxyy_1, g_z_xxyyy_1, g_z_xxyyz_1, g_z_xxyz_1, g_z_xxyzz_1, g_z_xxzz_1, g_z_xxzzz_1, g_z_xyyy_1, g_z_xyyyy_1, g_z_xyyyz_1, g_z_xyyz_1, g_z_xyyzz_1, g_z_xyzz_1, g_z_xyzzz_1, g_z_xzzz_1, g_z_xzzzz_1, g_z_yyyy_1, g_z_yyyyy_1, g_z_yyyyz_1, g_z_yyyz_1, g_z_yyyzz_1, g_z_yyzz_1, g_z_yyzzz_1, g_z_yzzz_1, g_z_yzzzz_1, g_z_zzzz_1, g_z_zzzzz_1, g_zz_xxxxx_0, g_zz_xxxxy_0, g_zz_xxxxz_0, g_zz_xxxyy_0, g_zz_xxxyz_0, g_zz_xxxzz_0, g_zz_xxyyy_0, g_zz_xxyyz_0, g_zz_xxyzz_0, g_zz_xxzzz_0, g_zz_xyyyy_0, g_zz_xyyyz_0, g_zz_xyyzz_0, g_zz_xyzzz_0, g_zz_xzzzz_0, g_zz_yyyyy_0, g_zz_yyyyz_0, g_zz_yyyzz_0, g_zz_yyzzz_0, g_zz_yzzzz_0, g_zz_zzzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zz_xxxxx_0[i] = g_0_xxxxx_0[i] * fbe_0 - g_0_xxxxx_1[i] * fz_be_0 + g_z_xxxxx_1[i] * pa_z[i];

        g_zz_xxxxy_0[i] = g_0_xxxxy_0[i] * fbe_0 - g_0_xxxxy_1[i] * fz_be_0 + g_z_xxxxy_1[i] * pa_z[i];

        g_zz_xxxxz_0[i] = g_0_xxxxz_0[i] * fbe_0 - g_0_xxxxz_1[i] * fz_be_0 + g_z_xxxx_1[i] * fe_0 + g_z_xxxxz_1[i] * pa_z[i];

        g_zz_xxxyy_0[i] = g_0_xxxyy_0[i] * fbe_0 - g_0_xxxyy_1[i] * fz_be_0 + g_z_xxxyy_1[i] * pa_z[i];

        g_zz_xxxyz_0[i] = g_0_xxxyz_0[i] * fbe_0 - g_0_xxxyz_1[i] * fz_be_0 + g_z_xxxy_1[i] * fe_0 + g_z_xxxyz_1[i] * pa_z[i];

        g_zz_xxxzz_0[i] = g_0_xxxzz_0[i] * fbe_0 - g_0_xxxzz_1[i] * fz_be_0 + 2.0 * g_z_xxxz_1[i] * fe_0 + g_z_xxxzz_1[i] * pa_z[i];

        g_zz_xxyyy_0[i] = g_0_xxyyy_0[i] * fbe_0 - g_0_xxyyy_1[i] * fz_be_0 + g_z_xxyyy_1[i] * pa_z[i];

        g_zz_xxyyz_0[i] = g_0_xxyyz_0[i] * fbe_0 - g_0_xxyyz_1[i] * fz_be_0 + g_z_xxyy_1[i] * fe_0 + g_z_xxyyz_1[i] * pa_z[i];

        g_zz_xxyzz_0[i] = g_0_xxyzz_0[i] * fbe_0 - g_0_xxyzz_1[i] * fz_be_0 + 2.0 * g_z_xxyz_1[i] * fe_0 + g_z_xxyzz_1[i] * pa_z[i];

        g_zz_xxzzz_0[i] = g_0_xxzzz_0[i] * fbe_0 - g_0_xxzzz_1[i] * fz_be_0 + 3.0 * g_z_xxzz_1[i] * fe_0 + g_z_xxzzz_1[i] * pa_z[i];

        g_zz_xyyyy_0[i] = g_0_xyyyy_0[i] * fbe_0 - g_0_xyyyy_1[i] * fz_be_0 + g_z_xyyyy_1[i] * pa_z[i];

        g_zz_xyyyz_0[i] = g_0_xyyyz_0[i] * fbe_0 - g_0_xyyyz_1[i] * fz_be_0 + g_z_xyyy_1[i] * fe_0 + g_z_xyyyz_1[i] * pa_z[i];

        g_zz_xyyzz_0[i] = g_0_xyyzz_0[i] * fbe_0 - g_0_xyyzz_1[i] * fz_be_0 + 2.0 * g_z_xyyz_1[i] * fe_0 + g_z_xyyzz_1[i] * pa_z[i];

        g_zz_xyzzz_0[i] = g_0_xyzzz_0[i] * fbe_0 - g_0_xyzzz_1[i] * fz_be_0 + 3.0 * g_z_xyzz_1[i] * fe_0 + g_z_xyzzz_1[i] * pa_z[i];

        g_zz_xzzzz_0[i] = g_0_xzzzz_0[i] * fbe_0 - g_0_xzzzz_1[i] * fz_be_0 + 4.0 * g_z_xzzz_1[i] * fe_0 + g_z_xzzzz_1[i] * pa_z[i];

        g_zz_yyyyy_0[i] = g_0_yyyyy_0[i] * fbe_0 - g_0_yyyyy_1[i] * fz_be_0 + g_z_yyyyy_1[i] * pa_z[i];

        g_zz_yyyyz_0[i] = g_0_yyyyz_0[i] * fbe_0 - g_0_yyyyz_1[i] * fz_be_0 + g_z_yyyy_1[i] * fe_0 + g_z_yyyyz_1[i] * pa_z[i];

        g_zz_yyyzz_0[i] = g_0_yyyzz_0[i] * fbe_0 - g_0_yyyzz_1[i] * fz_be_0 + 2.0 * g_z_yyyz_1[i] * fe_0 + g_z_yyyzz_1[i] * pa_z[i];

        g_zz_yyzzz_0[i] = g_0_yyzzz_0[i] * fbe_0 - g_0_yyzzz_1[i] * fz_be_0 + 3.0 * g_z_yyzz_1[i] * fe_0 + g_z_yyzzz_1[i] * pa_z[i];

        g_zz_yzzzz_0[i] = g_0_yzzzz_0[i] * fbe_0 - g_0_yzzzz_1[i] * fz_be_0 + 4.0 * g_z_yzzz_1[i] * fe_0 + g_z_yzzzz_1[i] * pa_z[i];

        g_zz_zzzzz_0[i] = g_0_zzzzz_0[i] * fbe_0 - g_0_zzzzz_1[i] * fz_be_0 + 5.0 * g_z_zzzz_1[i] * fe_0 + g_z_zzzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

