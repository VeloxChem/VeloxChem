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

#include "ThreeCenterElectronRepulsionPrimRecISP.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_isp(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_isp,
                                 size_t idx_eri_0_gsp,
                                 size_t idx_eri_1_gsp,
                                 size_t idx_eri_1_hss,
                                 size_t idx_eri_1_hsp,
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

    /// Set up components of auxilary buffer : GSP

    auto g_xxxx_0_x_0 = pbuffer.data(idx_eri_0_gsp);

    auto g_xxxx_0_y_0 = pbuffer.data(idx_eri_0_gsp + 1);

    auto g_xxxx_0_z_0 = pbuffer.data(idx_eri_0_gsp + 2);

    auto g_xxxy_0_x_0 = pbuffer.data(idx_eri_0_gsp + 3);

    auto g_xxxz_0_x_0 = pbuffer.data(idx_eri_0_gsp + 6);

    auto g_xxyy_0_x_0 = pbuffer.data(idx_eri_0_gsp + 9);

    auto g_xxyy_0_y_0 = pbuffer.data(idx_eri_0_gsp + 10);

    auto g_xxyy_0_z_0 = pbuffer.data(idx_eri_0_gsp + 11);

    auto g_xxzz_0_x_0 = pbuffer.data(idx_eri_0_gsp + 15);

    auto g_xxzz_0_y_0 = pbuffer.data(idx_eri_0_gsp + 16);

    auto g_xxzz_0_z_0 = pbuffer.data(idx_eri_0_gsp + 17);

    auto g_xyyy_0_y_0 = pbuffer.data(idx_eri_0_gsp + 19);

    auto g_xyyy_0_z_0 = pbuffer.data(idx_eri_0_gsp + 20);

    auto g_xzzz_0_y_0 = pbuffer.data(idx_eri_0_gsp + 28);

    auto g_xzzz_0_z_0 = pbuffer.data(idx_eri_0_gsp + 29);

    auto g_yyyy_0_x_0 = pbuffer.data(idx_eri_0_gsp + 30);

    auto g_yyyy_0_y_0 = pbuffer.data(idx_eri_0_gsp + 31);

    auto g_yyyy_0_z_0 = pbuffer.data(idx_eri_0_gsp + 32);

    auto g_yyyz_0_y_0 = pbuffer.data(idx_eri_0_gsp + 34);

    auto g_yyzz_0_x_0 = pbuffer.data(idx_eri_0_gsp + 36);

    auto g_yyzz_0_y_0 = pbuffer.data(idx_eri_0_gsp + 37);

    auto g_yyzz_0_z_0 = pbuffer.data(idx_eri_0_gsp + 38);

    auto g_yzzz_0_x_0 = pbuffer.data(idx_eri_0_gsp + 39);

    auto g_yzzz_0_z_0 = pbuffer.data(idx_eri_0_gsp + 41);

    auto g_zzzz_0_x_0 = pbuffer.data(idx_eri_0_gsp + 42);

    auto g_zzzz_0_y_0 = pbuffer.data(idx_eri_0_gsp + 43);

    auto g_zzzz_0_z_0 = pbuffer.data(idx_eri_0_gsp + 44);

    /// Set up components of auxilary buffer : GSP

    auto g_xxxx_0_x_1 = pbuffer.data(idx_eri_1_gsp);

    auto g_xxxx_0_y_1 = pbuffer.data(idx_eri_1_gsp + 1);

    auto g_xxxx_0_z_1 = pbuffer.data(idx_eri_1_gsp + 2);

    auto g_xxxy_0_x_1 = pbuffer.data(idx_eri_1_gsp + 3);

    auto g_xxxz_0_x_1 = pbuffer.data(idx_eri_1_gsp + 6);

    auto g_xxyy_0_x_1 = pbuffer.data(idx_eri_1_gsp + 9);

    auto g_xxyy_0_y_1 = pbuffer.data(idx_eri_1_gsp + 10);

    auto g_xxyy_0_z_1 = pbuffer.data(idx_eri_1_gsp + 11);

    auto g_xxzz_0_x_1 = pbuffer.data(idx_eri_1_gsp + 15);

    auto g_xxzz_0_y_1 = pbuffer.data(idx_eri_1_gsp + 16);

    auto g_xxzz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 17);

    auto g_xyyy_0_y_1 = pbuffer.data(idx_eri_1_gsp + 19);

    auto g_xyyy_0_z_1 = pbuffer.data(idx_eri_1_gsp + 20);

    auto g_xzzz_0_y_1 = pbuffer.data(idx_eri_1_gsp + 28);

    auto g_xzzz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 29);

    auto g_yyyy_0_x_1 = pbuffer.data(idx_eri_1_gsp + 30);

    auto g_yyyy_0_y_1 = pbuffer.data(idx_eri_1_gsp + 31);

    auto g_yyyy_0_z_1 = pbuffer.data(idx_eri_1_gsp + 32);

    auto g_yyyz_0_y_1 = pbuffer.data(idx_eri_1_gsp + 34);

    auto g_yyzz_0_x_1 = pbuffer.data(idx_eri_1_gsp + 36);

    auto g_yyzz_0_y_1 = pbuffer.data(idx_eri_1_gsp + 37);

    auto g_yyzz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 38);

    auto g_yzzz_0_x_1 = pbuffer.data(idx_eri_1_gsp + 39);

    auto g_yzzz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 41);

    auto g_zzzz_0_x_1 = pbuffer.data(idx_eri_1_gsp + 42);

    auto g_zzzz_0_y_1 = pbuffer.data(idx_eri_1_gsp + 43);

    auto g_zzzz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 44);

    /// Set up components of auxilary buffer : HSS

    auto g_xxxxx_0_0_1 = pbuffer.data(idx_eri_1_hss);

    auto g_xxxyy_0_0_1 = pbuffer.data(idx_eri_1_hss + 3);

    auto g_xxxzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 5);

    auto g_xxyyy_0_0_1 = pbuffer.data(idx_eri_1_hss + 6);

    auto g_xxzzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 9);

    auto g_yyyyy_0_0_1 = pbuffer.data(idx_eri_1_hss + 15);

    auto g_yyyzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 17);

    auto g_yyzzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 18);

    auto g_zzzzz_0_0_1 = pbuffer.data(idx_eri_1_hss + 20);

    /// Set up components of auxilary buffer : HSP

    auto g_xxxxx_0_x_1 = pbuffer.data(idx_eri_1_hsp);

    auto g_xxxxx_0_y_1 = pbuffer.data(idx_eri_1_hsp + 1);

    auto g_xxxxx_0_z_1 = pbuffer.data(idx_eri_1_hsp + 2);

    auto g_xxxxy_0_x_1 = pbuffer.data(idx_eri_1_hsp + 3);

    auto g_xxxxy_0_y_1 = pbuffer.data(idx_eri_1_hsp + 4);

    auto g_xxxxz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 6);

    auto g_xxxxz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 8);

    auto g_xxxyy_0_x_1 = pbuffer.data(idx_eri_1_hsp + 9);

    auto g_xxxyy_0_y_1 = pbuffer.data(idx_eri_1_hsp + 10);

    auto g_xxxyy_0_z_1 = pbuffer.data(idx_eri_1_hsp + 11);

    auto g_xxxzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 15);

    auto g_xxxzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 16);

    auto g_xxxzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 17);

    auto g_xxyyy_0_x_1 = pbuffer.data(idx_eri_1_hsp + 18);

    auto g_xxyyy_0_y_1 = pbuffer.data(idx_eri_1_hsp + 19);

    auto g_xxyyy_0_z_1 = pbuffer.data(idx_eri_1_hsp + 20);

    auto g_xxyzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 24);

    auto g_xxzzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 27);

    auto g_xxzzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 28);

    auto g_xxzzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 29);

    auto g_xyyyy_0_x_1 = pbuffer.data(idx_eri_1_hsp + 30);

    auto g_xyyyy_0_y_1 = pbuffer.data(idx_eri_1_hsp + 31);

    auto g_xyyyy_0_z_1 = pbuffer.data(idx_eri_1_hsp + 32);

    auto g_xyyzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 37);

    auto g_xyyzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 38);

    auto g_xzzzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 42);

    auto g_xzzzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 43);

    auto g_xzzzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 44);

    auto g_yyyyy_0_x_1 = pbuffer.data(idx_eri_1_hsp + 45);

    auto g_yyyyy_0_y_1 = pbuffer.data(idx_eri_1_hsp + 46);

    auto g_yyyyy_0_z_1 = pbuffer.data(idx_eri_1_hsp + 47);

    auto g_yyyyz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 49);

    auto g_yyyyz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 50);

    auto g_yyyzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 51);

    auto g_yyyzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 52);

    auto g_yyyzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 53);

    auto g_yyzzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 54);

    auto g_yyzzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 55);

    auto g_yyzzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 56);

    auto g_yzzzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 57);

    auto g_yzzzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 58);

    auto g_yzzzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 59);

    auto g_zzzzz_0_x_1 = pbuffer.data(idx_eri_1_hsp + 60);

    auto g_zzzzz_0_y_1 = pbuffer.data(idx_eri_1_hsp + 61);

    auto g_zzzzz_0_z_1 = pbuffer.data(idx_eri_1_hsp + 62);

    /// Set up 0-3 components of targeted buffer : ISP

    auto g_xxxxxx_0_x_0 = pbuffer.data(idx_eri_0_isp);

    auto g_xxxxxx_0_y_0 = pbuffer.data(idx_eri_0_isp + 1);

    auto g_xxxxxx_0_z_0 = pbuffer.data(idx_eri_0_isp + 2);

    #pragma omp simd aligned(g_xxxx_0_x_0, g_xxxx_0_x_1, g_xxxx_0_y_0, g_xxxx_0_y_1, g_xxxx_0_z_0, g_xxxx_0_z_1, g_xxxxx_0_0_1, g_xxxxx_0_x_1, g_xxxxx_0_y_1, g_xxxxx_0_z_1, g_xxxxxx_0_x_0, g_xxxxxx_0_y_0, g_xxxxxx_0_z_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxxx_0_x_0[i] = 5.0 * g_xxxx_0_x_0[i] * fbe_0 - 5.0 * g_xxxx_0_x_1[i] * fz_be_0 + g_xxxxx_0_0_1[i] * fi_acd_0 + g_xxxxx_0_x_1[i] * wa_x[i];

        g_xxxxxx_0_y_0[i] = 5.0 * g_xxxx_0_y_0[i] * fbe_0 - 5.0 * g_xxxx_0_y_1[i] * fz_be_0 + g_xxxxx_0_y_1[i] * wa_x[i];

        g_xxxxxx_0_z_0[i] = 5.0 * g_xxxx_0_z_0[i] * fbe_0 - 5.0 * g_xxxx_0_z_1[i] * fz_be_0 + g_xxxxx_0_z_1[i] * wa_x[i];
    }

    /// Set up 3-6 components of targeted buffer : ISP

    auto g_xxxxxy_0_x_0 = pbuffer.data(idx_eri_0_isp + 3);

    auto g_xxxxxy_0_y_0 = pbuffer.data(idx_eri_0_isp + 4);

    auto g_xxxxxy_0_z_0 = pbuffer.data(idx_eri_0_isp + 5);

    #pragma omp simd aligned(g_xxxxx_0_0_1, g_xxxxx_0_x_1, g_xxxxx_0_y_1, g_xxxxx_0_z_1, g_xxxxxy_0_x_0, g_xxxxxy_0_y_0, g_xxxxxy_0_z_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxy_0_x_0[i] = g_xxxxx_0_x_1[i] * wa_y[i];

        g_xxxxxy_0_y_0[i] = g_xxxxx_0_0_1[i] * fi_acd_0 + g_xxxxx_0_y_1[i] * wa_y[i];

        g_xxxxxy_0_z_0[i] = g_xxxxx_0_z_1[i] * wa_y[i];
    }

    /// Set up 6-9 components of targeted buffer : ISP

    auto g_xxxxxz_0_x_0 = pbuffer.data(idx_eri_0_isp + 6);

    auto g_xxxxxz_0_y_0 = pbuffer.data(idx_eri_0_isp + 7);

    auto g_xxxxxz_0_z_0 = pbuffer.data(idx_eri_0_isp + 8);

    #pragma omp simd aligned(g_xxxxx_0_0_1, g_xxxxx_0_x_1, g_xxxxx_0_y_1, g_xxxxx_0_z_1, g_xxxxxz_0_x_0, g_xxxxxz_0_y_0, g_xxxxxz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxxz_0_x_0[i] = g_xxxxx_0_x_1[i] * wa_z[i];

        g_xxxxxz_0_y_0[i] = g_xxxxx_0_y_1[i] * wa_z[i];

        g_xxxxxz_0_z_0[i] = g_xxxxx_0_0_1[i] * fi_acd_0 + g_xxxxx_0_z_1[i] * wa_z[i];
    }

    /// Set up 9-12 components of targeted buffer : ISP

    auto g_xxxxyy_0_x_0 = pbuffer.data(idx_eri_0_isp + 9);

    auto g_xxxxyy_0_y_0 = pbuffer.data(idx_eri_0_isp + 10);

    auto g_xxxxyy_0_z_0 = pbuffer.data(idx_eri_0_isp + 11);

    #pragma omp simd aligned(g_xxxx_0_x_0, g_xxxx_0_x_1, g_xxxxy_0_x_1, g_xxxxyy_0_x_0, g_xxxxyy_0_y_0, g_xxxxyy_0_z_0, g_xxxyy_0_y_1, g_xxxyy_0_z_1, g_xxyy_0_y_0, g_xxyy_0_y_1, g_xxyy_0_z_0, g_xxyy_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxyy_0_x_0[i] = g_xxxx_0_x_0[i] * fbe_0 - g_xxxx_0_x_1[i] * fz_be_0 + g_xxxxy_0_x_1[i] * wa_y[i];

        g_xxxxyy_0_y_0[i] = 3.0 * g_xxyy_0_y_0[i] * fbe_0 - 3.0 * g_xxyy_0_y_1[i] * fz_be_0 + g_xxxyy_0_y_1[i] * wa_x[i];

        g_xxxxyy_0_z_0[i] = 3.0 * g_xxyy_0_z_0[i] * fbe_0 - 3.0 * g_xxyy_0_z_1[i] * fz_be_0 + g_xxxyy_0_z_1[i] * wa_x[i];
    }

    /// Set up 12-15 components of targeted buffer : ISP

    auto g_xxxxyz_0_x_0 = pbuffer.data(idx_eri_0_isp + 12);

    auto g_xxxxyz_0_y_0 = pbuffer.data(idx_eri_0_isp + 13);

    auto g_xxxxyz_0_z_0 = pbuffer.data(idx_eri_0_isp + 14);

    #pragma omp simd aligned(g_xxxxy_0_y_1, g_xxxxyz_0_x_0, g_xxxxyz_0_y_0, g_xxxxyz_0_z_0, g_xxxxz_0_x_1, g_xxxxz_0_z_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xxxxyz_0_x_0[i] = g_xxxxz_0_x_1[i] * wa_y[i];

        g_xxxxyz_0_y_0[i] = g_xxxxy_0_y_1[i] * wa_z[i];

        g_xxxxyz_0_z_0[i] = g_xxxxz_0_z_1[i] * wa_y[i];
    }

    /// Set up 15-18 components of targeted buffer : ISP

    auto g_xxxxzz_0_x_0 = pbuffer.data(idx_eri_0_isp + 15);

    auto g_xxxxzz_0_y_0 = pbuffer.data(idx_eri_0_isp + 16);

    auto g_xxxxzz_0_z_0 = pbuffer.data(idx_eri_0_isp + 17);

    #pragma omp simd aligned(g_xxxx_0_x_0, g_xxxx_0_x_1, g_xxxxz_0_x_1, g_xxxxzz_0_x_0, g_xxxxzz_0_y_0, g_xxxxzz_0_z_0, g_xxxzz_0_y_1, g_xxxzz_0_z_1, g_xxzz_0_y_0, g_xxzz_0_y_1, g_xxzz_0_z_0, g_xxzz_0_z_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxzz_0_x_0[i] = g_xxxx_0_x_0[i] * fbe_0 - g_xxxx_0_x_1[i] * fz_be_0 + g_xxxxz_0_x_1[i] * wa_z[i];

        g_xxxxzz_0_y_0[i] = 3.0 * g_xxzz_0_y_0[i] * fbe_0 - 3.0 * g_xxzz_0_y_1[i] * fz_be_0 + g_xxxzz_0_y_1[i] * wa_x[i];

        g_xxxxzz_0_z_0[i] = 3.0 * g_xxzz_0_z_0[i] * fbe_0 - 3.0 * g_xxzz_0_z_1[i] * fz_be_0 + g_xxxzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 18-21 components of targeted buffer : ISP

    auto g_xxxyyy_0_x_0 = pbuffer.data(idx_eri_0_isp + 18);

    auto g_xxxyyy_0_y_0 = pbuffer.data(idx_eri_0_isp + 19);

    auto g_xxxyyy_0_z_0 = pbuffer.data(idx_eri_0_isp + 20);

    #pragma omp simd aligned(g_xxxy_0_x_0, g_xxxy_0_x_1, g_xxxyy_0_x_1, g_xxxyyy_0_x_0, g_xxxyyy_0_y_0, g_xxxyyy_0_z_0, g_xxyyy_0_y_1, g_xxyyy_0_z_1, g_xyyy_0_y_0, g_xyyy_0_y_1, g_xyyy_0_z_0, g_xyyy_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyy_0_x_0[i] = 2.0 * g_xxxy_0_x_0[i] * fbe_0 - 2.0 * g_xxxy_0_x_1[i] * fz_be_0 + g_xxxyy_0_x_1[i] * wa_y[i];

        g_xxxyyy_0_y_0[i] = 2.0 * g_xyyy_0_y_0[i] * fbe_0 - 2.0 * g_xyyy_0_y_1[i] * fz_be_0 + g_xxyyy_0_y_1[i] * wa_x[i];

        g_xxxyyy_0_z_0[i] = 2.0 * g_xyyy_0_z_0[i] * fbe_0 - 2.0 * g_xyyy_0_z_1[i] * fz_be_0 + g_xxyyy_0_z_1[i] * wa_x[i];
    }

    /// Set up 21-24 components of targeted buffer : ISP

    auto g_xxxyyz_0_x_0 = pbuffer.data(idx_eri_0_isp + 21);

    auto g_xxxyyz_0_y_0 = pbuffer.data(idx_eri_0_isp + 22);

    auto g_xxxyyz_0_z_0 = pbuffer.data(idx_eri_0_isp + 23);

    #pragma omp simd aligned(g_xxxyy_0_0_1, g_xxxyy_0_x_1, g_xxxyy_0_y_1, g_xxxyy_0_z_1, g_xxxyyz_0_x_0, g_xxxyyz_0_y_0, g_xxxyyz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyyz_0_x_0[i] = g_xxxyy_0_x_1[i] * wa_z[i];

        g_xxxyyz_0_y_0[i] = g_xxxyy_0_y_1[i] * wa_z[i];

        g_xxxyyz_0_z_0[i] = g_xxxyy_0_0_1[i] * fi_acd_0 + g_xxxyy_0_z_1[i] * wa_z[i];
    }

    /// Set up 24-27 components of targeted buffer : ISP

    auto g_xxxyzz_0_x_0 = pbuffer.data(idx_eri_0_isp + 24);

    auto g_xxxyzz_0_y_0 = pbuffer.data(idx_eri_0_isp + 25);

    auto g_xxxyzz_0_z_0 = pbuffer.data(idx_eri_0_isp + 26);

    #pragma omp simd aligned(g_xxxyzz_0_x_0, g_xxxyzz_0_y_0, g_xxxyzz_0_z_0, g_xxxzz_0_0_1, g_xxxzz_0_x_1, g_xxxzz_0_y_1, g_xxxzz_0_z_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyzz_0_x_0[i] = g_xxxzz_0_x_1[i] * wa_y[i];

        g_xxxyzz_0_y_0[i] = g_xxxzz_0_0_1[i] * fi_acd_0 + g_xxxzz_0_y_1[i] * wa_y[i];

        g_xxxyzz_0_z_0[i] = g_xxxzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 27-30 components of targeted buffer : ISP

    auto g_xxxzzz_0_x_0 = pbuffer.data(idx_eri_0_isp + 27);

    auto g_xxxzzz_0_y_0 = pbuffer.data(idx_eri_0_isp + 28);

    auto g_xxxzzz_0_z_0 = pbuffer.data(idx_eri_0_isp + 29);

    #pragma omp simd aligned(g_xxxz_0_x_0, g_xxxz_0_x_1, g_xxxzz_0_x_1, g_xxxzzz_0_x_0, g_xxxzzz_0_y_0, g_xxxzzz_0_z_0, g_xxzzz_0_y_1, g_xxzzz_0_z_1, g_xzzz_0_y_0, g_xzzz_0_y_1, g_xzzz_0_z_0, g_xzzz_0_z_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxzzz_0_x_0[i] = 2.0 * g_xxxz_0_x_0[i] * fbe_0 - 2.0 * g_xxxz_0_x_1[i] * fz_be_0 + g_xxxzz_0_x_1[i] * wa_z[i];

        g_xxxzzz_0_y_0[i] = 2.0 * g_xzzz_0_y_0[i] * fbe_0 - 2.0 * g_xzzz_0_y_1[i] * fz_be_0 + g_xxzzz_0_y_1[i] * wa_x[i];

        g_xxxzzz_0_z_0[i] = 2.0 * g_xzzz_0_z_0[i] * fbe_0 - 2.0 * g_xzzz_0_z_1[i] * fz_be_0 + g_xxzzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 30-33 components of targeted buffer : ISP

    auto g_xxyyyy_0_x_0 = pbuffer.data(idx_eri_0_isp + 30);

    auto g_xxyyyy_0_y_0 = pbuffer.data(idx_eri_0_isp + 31);

    auto g_xxyyyy_0_z_0 = pbuffer.data(idx_eri_0_isp + 32);

    #pragma omp simd aligned(g_xxyy_0_x_0, g_xxyy_0_x_1, g_xxyyy_0_x_1, g_xxyyyy_0_x_0, g_xxyyyy_0_y_0, g_xxyyyy_0_z_0, g_xyyyy_0_y_1, g_xyyyy_0_z_1, g_yyyy_0_y_0, g_yyyy_0_y_1, g_yyyy_0_z_0, g_yyyy_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyy_0_x_0[i] = 3.0 * g_xxyy_0_x_0[i] * fbe_0 - 3.0 * g_xxyy_0_x_1[i] * fz_be_0 + g_xxyyy_0_x_1[i] * wa_y[i];

        g_xxyyyy_0_y_0[i] = g_yyyy_0_y_0[i] * fbe_0 - g_yyyy_0_y_1[i] * fz_be_0 + g_xyyyy_0_y_1[i] * wa_x[i];

        g_xxyyyy_0_z_0[i] = g_yyyy_0_z_0[i] * fbe_0 - g_yyyy_0_z_1[i] * fz_be_0 + g_xyyyy_0_z_1[i] * wa_x[i];
    }

    /// Set up 33-36 components of targeted buffer : ISP

    auto g_xxyyyz_0_x_0 = pbuffer.data(idx_eri_0_isp + 33);

    auto g_xxyyyz_0_y_0 = pbuffer.data(idx_eri_0_isp + 34);

    auto g_xxyyyz_0_z_0 = pbuffer.data(idx_eri_0_isp + 35);

    #pragma omp simd aligned(g_xxyyy_0_0_1, g_xxyyy_0_x_1, g_xxyyy_0_y_1, g_xxyyy_0_z_1, g_xxyyyz_0_x_0, g_xxyyyz_0_y_0, g_xxyyyz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyyz_0_x_0[i] = g_xxyyy_0_x_1[i] * wa_z[i];

        g_xxyyyz_0_y_0[i] = g_xxyyy_0_y_1[i] * wa_z[i];

        g_xxyyyz_0_z_0[i] = g_xxyyy_0_0_1[i] * fi_acd_0 + g_xxyyy_0_z_1[i] * wa_z[i];
    }

    /// Set up 36-39 components of targeted buffer : ISP

    auto g_xxyyzz_0_x_0 = pbuffer.data(idx_eri_0_isp + 36);

    auto g_xxyyzz_0_y_0 = pbuffer.data(idx_eri_0_isp + 37);

    auto g_xxyyzz_0_z_0 = pbuffer.data(idx_eri_0_isp + 38);

    #pragma omp simd aligned(g_xxyyzz_0_x_0, g_xxyyzz_0_y_0, g_xxyyzz_0_z_0, g_xxyzz_0_x_1, g_xxzz_0_x_0, g_xxzz_0_x_1, g_xyyzz_0_y_1, g_xyyzz_0_z_1, g_yyzz_0_y_0, g_yyzz_0_y_1, g_yyzz_0_z_0, g_yyzz_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyzz_0_x_0[i] = g_xxzz_0_x_0[i] * fbe_0 - g_xxzz_0_x_1[i] * fz_be_0 + g_xxyzz_0_x_1[i] * wa_y[i];

        g_xxyyzz_0_y_0[i] = g_yyzz_0_y_0[i] * fbe_0 - g_yyzz_0_y_1[i] * fz_be_0 + g_xyyzz_0_y_1[i] * wa_x[i];

        g_xxyyzz_0_z_0[i] = g_yyzz_0_z_0[i] * fbe_0 - g_yyzz_0_z_1[i] * fz_be_0 + g_xyyzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 39-42 components of targeted buffer : ISP

    auto g_xxyzzz_0_x_0 = pbuffer.data(idx_eri_0_isp + 39);

    auto g_xxyzzz_0_y_0 = pbuffer.data(idx_eri_0_isp + 40);

    auto g_xxyzzz_0_z_0 = pbuffer.data(idx_eri_0_isp + 41);

    #pragma omp simd aligned(g_xxyzzz_0_x_0, g_xxyzzz_0_y_0, g_xxyzzz_0_z_0, g_xxzzz_0_0_1, g_xxzzz_0_x_1, g_xxzzz_0_y_1, g_xxzzz_0_z_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzzz_0_x_0[i] = g_xxzzz_0_x_1[i] * wa_y[i];

        g_xxyzzz_0_y_0[i] = g_xxzzz_0_0_1[i] * fi_acd_0 + g_xxzzz_0_y_1[i] * wa_y[i];

        g_xxyzzz_0_z_0[i] = g_xxzzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 42-45 components of targeted buffer : ISP

    auto g_xxzzzz_0_x_0 = pbuffer.data(idx_eri_0_isp + 42);

    auto g_xxzzzz_0_y_0 = pbuffer.data(idx_eri_0_isp + 43);

    auto g_xxzzzz_0_z_0 = pbuffer.data(idx_eri_0_isp + 44);

    #pragma omp simd aligned(g_xxzz_0_x_0, g_xxzz_0_x_1, g_xxzzz_0_x_1, g_xxzzzz_0_x_0, g_xxzzzz_0_y_0, g_xxzzzz_0_z_0, g_xzzzz_0_y_1, g_xzzzz_0_z_1, g_zzzz_0_y_0, g_zzzz_0_y_1, g_zzzz_0_z_0, g_zzzz_0_z_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxzzzz_0_x_0[i] = 3.0 * g_xxzz_0_x_0[i] * fbe_0 - 3.0 * g_xxzz_0_x_1[i] * fz_be_0 + g_xxzzz_0_x_1[i] * wa_z[i];

        g_xxzzzz_0_y_0[i] = g_zzzz_0_y_0[i] * fbe_0 - g_zzzz_0_y_1[i] * fz_be_0 + g_xzzzz_0_y_1[i] * wa_x[i];

        g_xxzzzz_0_z_0[i] = g_zzzz_0_z_0[i] * fbe_0 - g_zzzz_0_z_1[i] * fz_be_0 + g_xzzzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 45-48 components of targeted buffer : ISP

    auto g_xyyyyy_0_x_0 = pbuffer.data(idx_eri_0_isp + 45);

    auto g_xyyyyy_0_y_0 = pbuffer.data(idx_eri_0_isp + 46);

    auto g_xyyyyy_0_z_0 = pbuffer.data(idx_eri_0_isp + 47);

    #pragma omp simd aligned(g_xyyyyy_0_x_0, g_xyyyyy_0_y_0, g_xyyyyy_0_z_0, g_yyyyy_0_0_1, g_yyyyy_0_x_1, g_yyyyy_0_y_1, g_yyyyy_0_z_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyyy_0_x_0[i] = g_yyyyy_0_0_1[i] * fi_acd_0 + g_yyyyy_0_x_1[i] * wa_x[i];

        g_xyyyyy_0_y_0[i] = g_yyyyy_0_y_1[i] * wa_x[i];

        g_xyyyyy_0_z_0[i] = g_yyyyy_0_z_1[i] * wa_x[i];
    }

    /// Set up 48-51 components of targeted buffer : ISP

    auto g_xyyyyz_0_x_0 = pbuffer.data(idx_eri_0_isp + 48);

    auto g_xyyyyz_0_y_0 = pbuffer.data(idx_eri_0_isp + 49);

    auto g_xyyyyz_0_z_0 = pbuffer.data(idx_eri_0_isp + 50);

    #pragma omp simd aligned(g_xyyyy_0_x_1, g_xyyyyz_0_x_0, g_xyyyyz_0_y_0, g_xyyyyz_0_z_0, g_yyyyz_0_y_1, g_yyyyz_0_z_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyyyyz_0_x_0[i] = g_xyyyy_0_x_1[i] * wa_z[i];

        g_xyyyyz_0_y_0[i] = g_yyyyz_0_y_1[i] * wa_x[i];

        g_xyyyyz_0_z_0[i] = g_yyyyz_0_z_1[i] * wa_x[i];
    }

    /// Set up 51-54 components of targeted buffer : ISP

    auto g_xyyyzz_0_x_0 = pbuffer.data(idx_eri_0_isp + 51);

    auto g_xyyyzz_0_y_0 = pbuffer.data(idx_eri_0_isp + 52);

    auto g_xyyyzz_0_z_0 = pbuffer.data(idx_eri_0_isp + 53);

    #pragma omp simd aligned(g_xyyyzz_0_x_0, g_xyyyzz_0_y_0, g_xyyyzz_0_z_0, g_yyyzz_0_0_1, g_yyyzz_0_x_1, g_yyyzz_0_y_1, g_yyyzz_0_z_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyzz_0_x_0[i] = g_yyyzz_0_0_1[i] * fi_acd_0 + g_yyyzz_0_x_1[i] * wa_x[i];

        g_xyyyzz_0_y_0[i] = g_yyyzz_0_y_1[i] * wa_x[i];

        g_xyyyzz_0_z_0[i] = g_yyyzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 54-57 components of targeted buffer : ISP

    auto g_xyyzzz_0_x_0 = pbuffer.data(idx_eri_0_isp + 54);

    auto g_xyyzzz_0_y_0 = pbuffer.data(idx_eri_0_isp + 55);

    auto g_xyyzzz_0_z_0 = pbuffer.data(idx_eri_0_isp + 56);

    #pragma omp simd aligned(g_xyyzzz_0_x_0, g_xyyzzz_0_y_0, g_xyyzzz_0_z_0, g_yyzzz_0_0_1, g_yyzzz_0_x_1, g_yyzzz_0_y_1, g_yyzzz_0_z_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzzz_0_x_0[i] = g_yyzzz_0_0_1[i] * fi_acd_0 + g_yyzzz_0_x_1[i] * wa_x[i];

        g_xyyzzz_0_y_0[i] = g_yyzzz_0_y_1[i] * wa_x[i];

        g_xyyzzz_0_z_0[i] = g_yyzzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 57-60 components of targeted buffer : ISP

    auto g_xyzzzz_0_x_0 = pbuffer.data(idx_eri_0_isp + 57);

    auto g_xyzzzz_0_y_0 = pbuffer.data(idx_eri_0_isp + 58);

    auto g_xyzzzz_0_z_0 = pbuffer.data(idx_eri_0_isp + 59);

    #pragma omp simd aligned(g_xyzzzz_0_x_0, g_xyzzzz_0_y_0, g_xyzzzz_0_z_0, g_xzzzz_0_x_1, g_yzzzz_0_y_1, g_yzzzz_0_z_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_xyzzzz_0_x_0[i] = g_xzzzz_0_x_1[i] * wa_y[i];

        g_xyzzzz_0_y_0[i] = g_yzzzz_0_y_1[i] * wa_x[i];

        g_xyzzzz_0_z_0[i] = g_yzzzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 60-63 components of targeted buffer : ISP

    auto g_xzzzzz_0_x_0 = pbuffer.data(idx_eri_0_isp + 60);

    auto g_xzzzzz_0_y_0 = pbuffer.data(idx_eri_0_isp + 61);

    auto g_xzzzzz_0_z_0 = pbuffer.data(idx_eri_0_isp + 62);

    #pragma omp simd aligned(g_xzzzzz_0_x_0, g_xzzzzz_0_y_0, g_xzzzzz_0_z_0, g_zzzzz_0_0_1, g_zzzzz_0_x_1, g_zzzzz_0_y_1, g_zzzzz_0_z_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzzz_0_x_0[i] = g_zzzzz_0_0_1[i] * fi_acd_0 + g_zzzzz_0_x_1[i] * wa_x[i];

        g_xzzzzz_0_y_0[i] = g_zzzzz_0_y_1[i] * wa_x[i];

        g_xzzzzz_0_z_0[i] = g_zzzzz_0_z_1[i] * wa_x[i];
    }

    /// Set up 63-66 components of targeted buffer : ISP

    auto g_yyyyyy_0_x_0 = pbuffer.data(idx_eri_0_isp + 63);

    auto g_yyyyyy_0_y_0 = pbuffer.data(idx_eri_0_isp + 64);

    auto g_yyyyyy_0_z_0 = pbuffer.data(idx_eri_0_isp + 65);

    #pragma omp simd aligned(g_yyyy_0_x_0, g_yyyy_0_x_1, g_yyyy_0_y_0, g_yyyy_0_y_1, g_yyyy_0_z_0, g_yyyy_0_z_1, g_yyyyy_0_0_1, g_yyyyy_0_x_1, g_yyyyy_0_y_1, g_yyyyy_0_z_1, g_yyyyyy_0_x_0, g_yyyyyy_0_y_0, g_yyyyyy_0_z_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyyy_0_x_0[i] = 5.0 * g_yyyy_0_x_0[i] * fbe_0 - 5.0 * g_yyyy_0_x_1[i] * fz_be_0 + g_yyyyy_0_x_1[i] * wa_y[i];

        g_yyyyyy_0_y_0[i] = 5.0 * g_yyyy_0_y_0[i] * fbe_0 - 5.0 * g_yyyy_0_y_1[i] * fz_be_0 + g_yyyyy_0_0_1[i] * fi_acd_0 + g_yyyyy_0_y_1[i] * wa_y[i];

        g_yyyyyy_0_z_0[i] = 5.0 * g_yyyy_0_z_0[i] * fbe_0 - 5.0 * g_yyyy_0_z_1[i] * fz_be_0 + g_yyyyy_0_z_1[i] * wa_y[i];
    }

    /// Set up 66-69 components of targeted buffer : ISP

    auto g_yyyyyz_0_x_0 = pbuffer.data(idx_eri_0_isp + 66);

    auto g_yyyyyz_0_y_0 = pbuffer.data(idx_eri_0_isp + 67);

    auto g_yyyyyz_0_z_0 = pbuffer.data(idx_eri_0_isp + 68);

    #pragma omp simd aligned(g_yyyyy_0_0_1, g_yyyyy_0_x_1, g_yyyyy_0_y_1, g_yyyyy_0_z_1, g_yyyyyz_0_x_0, g_yyyyyz_0_y_0, g_yyyyyz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyyz_0_x_0[i] = g_yyyyy_0_x_1[i] * wa_z[i];

        g_yyyyyz_0_y_0[i] = g_yyyyy_0_y_1[i] * wa_z[i];

        g_yyyyyz_0_z_0[i] = g_yyyyy_0_0_1[i] * fi_acd_0 + g_yyyyy_0_z_1[i] * wa_z[i];
    }

    /// Set up 69-72 components of targeted buffer : ISP

    auto g_yyyyzz_0_x_0 = pbuffer.data(idx_eri_0_isp + 69);

    auto g_yyyyzz_0_y_0 = pbuffer.data(idx_eri_0_isp + 70);

    auto g_yyyyzz_0_z_0 = pbuffer.data(idx_eri_0_isp + 71);

    #pragma omp simd aligned(g_yyyy_0_y_0, g_yyyy_0_y_1, g_yyyyz_0_y_1, g_yyyyzz_0_x_0, g_yyyyzz_0_y_0, g_yyyyzz_0_z_0, g_yyyzz_0_x_1, g_yyyzz_0_z_1, g_yyzz_0_x_0, g_yyzz_0_x_1, g_yyzz_0_z_0, g_yyzz_0_z_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyzz_0_x_0[i] = 3.0 * g_yyzz_0_x_0[i] * fbe_0 - 3.0 * g_yyzz_0_x_1[i] * fz_be_0 + g_yyyzz_0_x_1[i] * wa_y[i];

        g_yyyyzz_0_y_0[i] = g_yyyy_0_y_0[i] * fbe_0 - g_yyyy_0_y_1[i] * fz_be_0 + g_yyyyz_0_y_1[i] * wa_z[i];

        g_yyyyzz_0_z_0[i] = 3.0 * g_yyzz_0_z_0[i] * fbe_0 - 3.0 * g_yyzz_0_z_1[i] * fz_be_0 + g_yyyzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 72-75 components of targeted buffer : ISP

    auto g_yyyzzz_0_x_0 = pbuffer.data(idx_eri_0_isp + 72);

    auto g_yyyzzz_0_y_0 = pbuffer.data(idx_eri_0_isp + 73);

    auto g_yyyzzz_0_z_0 = pbuffer.data(idx_eri_0_isp + 74);

    #pragma omp simd aligned(g_yyyz_0_y_0, g_yyyz_0_y_1, g_yyyzz_0_y_1, g_yyyzzz_0_x_0, g_yyyzzz_0_y_0, g_yyyzzz_0_z_0, g_yyzzz_0_x_1, g_yyzzz_0_z_1, g_yzzz_0_x_0, g_yzzz_0_x_1, g_yzzz_0_z_0, g_yzzz_0_z_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyzzz_0_x_0[i] = 2.0 * g_yzzz_0_x_0[i] * fbe_0 - 2.0 * g_yzzz_0_x_1[i] * fz_be_0 + g_yyzzz_0_x_1[i] * wa_y[i];

        g_yyyzzz_0_y_0[i] = 2.0 * g_yyyz_0_y_0[i] * fbe_0 - 2.0 * g_yyyz_0_y_1[i] * fz_be_0 + g_yyyzz_0_y_1[i] * wa_z[i];

        g_yyyzzz_0_z_0[i] = 2.0 * g_yzzz_0_z_0[i] * fbe_0 - 2.0 * g_yzzz_0_z_1[i] * fz_be_0 + g_yyzzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 75-78 components of targeted buffer : ISP

    auto g_yyzzzz_0_x_0 = pbuffer.data(idx_eri_0_isp + 75);

    auto g_yyzzzz_0_y_0 = pbuffer.data(idx_eri_0_isp + 76);

    auto g_yyzzzz_0_z_0 = pbuffer.data(idx_eri_0_isp + 77);

    #pragma omp simd aligned(g_yyzz_0_y_0, g_yyzz_0_y_1, g_yyzzz_0_y_1, g_yyzzzz_0_x_0, g_yyzzzz_0_y_0, g_yyzzzz_0_z_0, g_yzzzz_0_x_1, g_yzzzz_0_z_1, g_zzzz_0_x_0, g_zzzz_0_x_1, g_zzzz_0_z_0, g_zzzz_0_z_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 = (c_exps[i] + d_exps[i]) * fbe_0 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyzzzz_0_x_0[i] = g_zzzz_0_x_0[i] * fbe_0 - g_zzzz_0_x_1[i] * fz_be_0 + g_yzzzz_0_x_1[i] * wa_y[i];

        g_yyzzzz_0_y_0[i] = 3.0 * g_yyzz_0_y_0[i] * fbe_0 - 3.0 * g_yyzz_0_y_1[i] * fz_be_0 + g_yyzzz_0_y_1[i] * wa_z[i];

        g_yyzzzz_0_z_0[i] = g_zzzz_0_z_0[i] * fbe_0 - g_zzzz_0_z_1[i] * fz_be_0 + g_yzzzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 78-81 components of targeted buffer : ISP

    auto g_yzzzzz_0_x_0 = pbuffer.data(idx_eri_0_isp + 78);

    auto g_yzzzzz_0_y_0 = pbuffer.data(idx_eri_0_isp + 79);

    auto g_yzzzzz_0_z_0 = pbuffer.data(idx_eri_0_isp + 80);

    #pragma omp simd aligned(g_yzzzzz_0_x_0, g_yzzzzz_0_y_0, g_yzzzzz_0_z_0, g_zzzzz_0_0_1, g_zzzzz_0_x_1, g_zzzzz_0_y_1, g_zzzzz_0_z_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzzz_0_x_0[i] = g_zzzzz_0_x_1[i] * wa_y[i];

        g_yzzzzz_0_y_0[i] = g_zzzzz_0_0_1[i] * fi_acd_0 + g_zzzzz_0_y_1[i] * wa_y[i];

        g_yzzzzz_0_z_0[i] = g_zzzzz_0_z_1[i] * wa_y[i];
    }

    /// Set up 81-84 components of targeted buffer : ISP

    auto g_zzzzzz_0_x_0 = pbuffer.data(idx_eri_0_isp + 81);

    auto g_zzzzzz_0_y_0 = pbuffer.data(idx_eri_0_isp + 82);

    auto g_zzzzzz_0_z_0 = pbuffer.data(idx_eri_0_isp + 83);

    #pragma omp simd aligned(g_zzzz_0_x_0, g_zzzz_0_x_1, g_zzzz_0_y_0, g_zzzz_0_y_1, g_zzzz_0_z_0, g_zzzz_0_z_1, g_zzzzz_0_0_1, g_zzzzz_0_x_1, g_zzzzz_0_y_1, g_zzzzz_0_z_1, g_zzzzzz_0_x_0, g_zzzzzz_0_y_0, g_zzzzzz_0_z_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzzz_0_x_0[i] = 5.0 * g_zzzz_0_x_0[i] * fbe_0 - 5.0 * g_zzzz_0_x_1[i] * fz_be_0 + g_zzzzz_0_x_1[i] * wa_z[i];

        g_zzzzzz_0_y_0[i] = 5.0 * g_zzzz_0_y_0[i] * fbe_0 - 5.0 * g_zzzz_0_y_1[i] * fz_be_0 + g_zzzzz_0_y_1[i] * wa_z[i];

        g_zzzzzz_0_z_0[i] = 5.0 * g_zzzz_0_z_0[i] * fbe_0 - 5.0 * g_zzzz_0_z_1[i] * fz_be_0 + g_zzzzz_0_0_1[i] * fi_acd_0 + g_zzzzz_0_z_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

