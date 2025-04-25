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

#include "ThreeCenterElectronRepulsionPrimRecHSD.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_hsd(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_hsd,
                                 size_t idx_eri_0_fsd,
                                 size_t idx_eri_1_fsd,
                                 size_t idx_eri_1_gsp,
                                 size_t idx_eri_1_gsd,
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

    /// Set up components of auxilary buffer : FSD

    auto g_xxx_0_xx_0 = pbuffer.data(idx_eri_0_fsd);

    auto g_xxx_0_xy_0 = pbuffer.data(idx_eri_0_fsd + 1);

    auto g_xxx_0_xz_0 = pbuffer.data(idx_eri_0_fsd + 2);

    auto g_xxx_0_yy_0 = pbuffer.data(idx_eri_0_fsd + 3);

    auto g_xxx_0_yz_0 = pbuffer.data(idx_eri_0_fsd + 4);

    auto g_xxx_0_zz_0 = pbuffer.data(idx_eri_0_fsd + 5);

    auto g_xxy_0_xx_0 = pbuffer.data(idx_eri_0_fsd + 6);

    auto g_xxy_0_xz_0 = pbuffer.data(idx_eri_0_fsd + 8);

    auto g_xxz_0_xx_0 = pbuffer.data(idx_eri_0_fsd + 12);

    auto g_xxz_0_xy_0 = pbuffer.data(idx_eri_0_fsd + 13);

    auto g_xyy_0_xy_0 = pbuffer.data(idx_eri_0_fsd + 19);

    auto g_xyy_0_yy_0 = pbuffer.data(idx_eri_0_fsd + 21);

    auto g_xyy_0_yz_0 = pbuffer.data(idx_eri_0_fsd + 22);

    auto g_xyy_0_zz_0 = pbuffer.data(idx_eri_0_fsd + 23);

    auto g_xzz_0_xz_0 = pbuffer.data(idx_eri_0_fsd + 32);

    auto g_xzz_0_yy_0 = pbuffer.data(idx_eri_0_fsd + 33);

    auto g_xzz_0_yz_0 = pbuffer.data(idx_eri_0_fsd + 34);

    auto g_xzz_0_zz_0 = pbuffer.data(idx_eri_0_fsd + 35);

    auto g_yyy_0_xx_0 = pbuffer.data(idx_eri_0_fsd + 36);

    auto g_yyy_0_xy_0 = pbuffer.data(idx_eri_0_fsd + 37);

    auto g_yyy_0_xz_0 = pbuffer.data(idx_eri_0_fsd + 38);

    auto g_yyy_0_yy_0 = pbuffer.data(idx_eri_0_fsd + 39);

    auto g_yyy_0_yz_0 = pbuffer.data(idx_eri_0_fsd + 40);

    auto g_yyy_0_zz_0 = pbuffer.data(idx_eri_0_fsd + 41);

    auto g_yyz_0_xy_0 = pbuffer.data(idx_eri_0_fsd + 43);

    auto g_yyz_0_yy_0 = pbuffer.data(idx_eri_0_fsd + 45);

    auto g_yzz_0_xx_0 = pbuffer.data(idx_eri_0_fsd + 48);

    auto g_yzz_0_xz_0 = pbuffer.data(idx_eri_0_fsd + 50);

    auto g_yzz_0_yz_0 = pbuffer.data(idx_eri_0_fsd + 52);

    auto g_yzz_0_zz_0 = pbuffer.data(idx_eri_0_fsd + 53);

    auto g_zzz_0_xx_0 = pbuffer.data(idx_eri_0_fsd + 54);

    auto g_zzz_0_xy_0 = pbuffer.data(idx_eri_0_fsd + 55);

    auto g_zzz_0_xz_0 = pbuffer.data(idx_eri_0_fsd + 56);

    auto g_zzz_0_yy_0 = pbuffer.data(idx_eri_0_fsd + 57);

    auto g_zzz_0_yz_0 = pbuffer.data(idx_eri_0_fsd + 58);

    auto g_zzz_0_zz_0 = pbuffer.data(idx_eri_0_fsd + 59);

    /// Set up components of auxilary buffer : FSD

    auto g_xxx_0_xx_1 = pbuffer.data(idx_eri_1_fsd);

    auto g_xxx_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 1);

    auto g_xxx_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 2);

    auto g_xxx_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 3);

    auto g_xxx_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 4);

    auto g_xxx_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 5);

    auto g_xxy_0_xx_1 = pbuffer.data(idx_eri_1_fsd + 6);

    auto g_xxy_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 8);

    auto g_xxz_0_xx_1 = pbuffer.data(idx_eri_1_fsd + 12);

    auto g_xxz_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 13);

    auto g_xyy_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 19);

    auto g_xyy_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 21);

    auto g_xyy_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 22);

    auto g_xyy_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 23);

    auto g_xzz_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 32);

    auto g_xzz_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 33);

    auto g_xzz_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 34);

    auto g_xzz_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 35);

    auto g_yyy_0_xx_1 = pbuffer.data(idx_eri_1_fsd + 36);

    auto g_yyy_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 37);

    auto g_yyy_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 38);

    auto g_yyy_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 39);

    auto g_yyy_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 40);

    auto g_yyy_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 41);

    auto g_yyz_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 43);

    auto g_yyz_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 45);

    auto g_yzz_0_xx_1 = pbuffer.data(idx_eri_1_fsd + 48);

    auto g_yzz_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 50);

    auto g_yzz_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 52);

    auto g_yzz_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 53);

    auto g_zzz_0_xx_1 = pbuffer.data(idx_eri_1_fsd + 54);

    auto g_zzz_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 55);

    auto g_zzz_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 56);

    auto g_zzz_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 57);

    auto g_zzz_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 58);

    auto g_zzz_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 59);

    /// Set up components of auxilary buffer : GSP

    auto g_xxxx_0_x_1 = pbuffer.data(idx_eri_1_gsp);

    auto g_xxxx_0_y_1 = pbuffer.data(idx_eri_1_gsp + 1);

    auto g_xxxx_0_z_1 = pbuffer.data(idx_eri_1_gsp + 2);

    auto g_xxxz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 8);

    auto g_xxyy_0_x_1 = pbuffer.data(idx_eri_1_gsp + 9);

    auto g_xxyy_0_y_1 = pbuffer.data(idx_eri_1_gsp + 10);

    auto g_xxyy_0_z_1 = pbuffer.data(idx_eri_1_gsp + 11);

    auto g_xxzz_0_x_1 = pbuffer.data(idx_eri_1_gsp + 15);

    auto g_xxzz_0_y_1 = pbuffer.data(idx_eri_1_gsp + 16);

    auto g_xxzz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 17);

    auto g_xyyy_0_y_1 = pbuffer.data(idx_eri_1_gsp + 19);

    auto g_xzzz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 29);

    auto g_yyyy_0_x_1 = pbuffer.data(idx_eri_1_gsp + 30);

    auto g_yyyy_0_y_1 = pbuffer.data(idx_eri_1_gsp + 31);

    auto g_yyyy_0_z_1 = pbuffer.data(idx_eri_1_gsp + 32);

    auto g_yyyz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 35);

    auto g_yyzz_0_x_1 = pbuffer.data(idx_eri_1_gsp + 36);

    auto g_yyzz_0_y_1 = pbuffer.data(idx_eri_1_gsp + 37);

    auto g_yyzz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 38);

    auto g_yzzz_0_y_1 = pbuffer.data(idx_eri_1_gsp + 40);

    auto g_yzzz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 41);

    auto g_zzzz_0_x_1 = pbuffer.data(idx_eri_1_gsp + 42);

    auto g_zzzz_0_y_1 = pbuffer.data(idx_eri_1_gsp + 43);

    auto g_zzzz_0_z_1 = pbuffer.data(idx_eri_1_gsp + 44);

    /// Set up components of auxilary buffer : GSD

    auto g_xxxx_0_xx_1 = pbuffer.data(idx_eri_1_gsd);

    auto g_xxxx_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 1);

    auto g_xxxx_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 2);

    auto g_xxxx_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 3);

    auto g_xxxx_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 4);

    auto g_xxxx_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 5);

    auto g_xxxy_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 6);

    auto g_xxxy_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 7);

    auto g_xxxy_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 8);

    auto g_xxxy_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 9);

    auto g_xxxz_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 12);

    auto g_xxxz_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 13);

    auto g_xxxz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 14);

    auto g_xxxz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 16);

    auto g_xxxz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 17);

    auto g_xxyy_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 18);

    auto g_xxyy_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 19);

    auto g_xxyy_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 20);

    auto g_xxyy_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 21);

    auto g_xxyy_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 22);

    auto g_xxyy_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 23);

    auto g_xxzz_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 30);

    auto g_xxzz_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 31);

    auto g_xxzz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 32);

    auto g_xxzz_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 33);

    auto g_xxzz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 34);

    auto g_xxzz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 35);

    auto g_xyyy_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 36);

    auto g_xyyy_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 37);

    auto g_xyyy_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 39);

    auto g_xyyy_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 40);

    auto g_xyyy_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 41);

    auto g_xzzz_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 54);

    auto g_xzzz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 56);

    auto g_xzzz_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 57);

    auto g_xzzz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 58);

    auto g_xzzz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 59);

    auto g_yyyy_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 60);

    auto g_yyyy_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 61);

    auto g_yyyy_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 62);

    auto g_yyyy_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 63);

    auto g_yyyy_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 64);

    auto g_yyyy_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 65);

    auto g_yyyz_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 67);

    auto g_yyyz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 68);

    auto g_yyyz_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 69);

    auto g_yyyz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 70);

    auto g_yyyz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 71);

    auto g_yyzz_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 72);

    auto g_yyzz_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 73);

    auto g_yyzz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 74);

    auto g_yyzz_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 75);

    auto g_yyzz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 76);

    auto g_yyzz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 77);

    auto g_yzzz_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 78);

    auto g_yzzz_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 79);

    auto g_yzzz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 80);

    auto g_yzzz_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 81);

    auto g_yzzz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 82);

    auto g_yzzz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 83);

    auto g_zzzz_0_xx_1 = pbuffer.data(idx_eri_1_gsd + 84);

    auto g_zzzz_0_xy_1 = pbuffer.data(idx_eri_1_gsd + 85);

    auto g_zzzz_0_xz_1 = pbuffer.data(idx_eri_1_gsd + 86);

    auto g_zzzz_0_yy_1 = pbuffer.data(idx_eri_1_gsd + 87);

    auto g_zzzz_0_yz_1 = pbuffer.data(idx_eri_1_gsd + 88);

    auto g_zzzz_0_zz_1 = pbuffer.data(idx_eri_1_gsd + 89);

    /// Set up 0-6 components of targeted buffer : HSD

    auto g_xxxxx_0_xx_0 = pbuffer.data(idx_eri_0_hsd);

    auto g_xxxxx_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 1);

    auto g_xxxxx_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 2);

    auto g_xxxxx_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 3);

    auto g_xxxxx_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 4);

    auto g_xxxxx_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 5);

    #pragma omp simd aligned(g_xxx_0_xx_0, g_xxx_0_xx_1, g_xxx_0_xy_0, g_xxx_0_xy_1, g_xxx_0_xz_0, g_xxx_0_xz_1, g_xxx_0_yy_0, g_xxx_0_yy_1, g_xxx_0_yz_0, g_xxx_0_yz_1, g_xxx_0_zz_0, g_xxx_0_zz_1, g_xxxx_0_x_1, g_xxxx_0_xx_1, g_xxxx_0_xy_1, g_xxxx_0_xz_1, g_xxxx_0_y_1, g_xxxx_0_yy_1, g_xxxx_0_yz_1, g_xxxx_0_z_1, g_xxxx_0_zz_1, g_xxxxx_0_xx_0, g_xxxxx_0_xy_0, g_xxxxx_0_xz_0, g_xxxxx_0_yy_0, g_xxxxx_0_yz_0, g_xxxxx_0_zz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxxx_0_xx_0[i] = 4.0 * g_xxx_0_xx_0[i] * fbe_0 - 4.0 * g_xxx_0_xx_1[i] * fz_be_0 + 2.0 * g_xxxx_0_x_1[i] * fi_acd_0 + g_xxxx_0_xx_1[i] * wa_x[i];

        g_xxxxx_0_xy_0[i] = 4.0 * g_xxx_0_xy_0[i] * fbe_0 - 4.0 * g_xxx_0_xy_1[i] * fz_be_0 + g_xxxx_0_y_1[i] * fi_acd_0 + g_xxxx_0_xy_1[i] * wa_x[i];

        g_xxxxx_0_xz_0[i] = 4.0 * g_xxx_0_xz_0[i] * fbe_0 - 4.0 * g_xxx_0_xz_1[i] * fz_be_0 + g_xxxx_0_z_1[i] * fi_acd_0 + g_xxxx_0_xz_1[i] * wa_x[i];

        g_xxxxx_0_yy_0[i] = 4.0 * g_xxx_0_yy_0[i] * fbe_0 - 4.0 * g_xxx_0_yy_1[i] * fz_be_0 + g_xxxx_0_yy_1[i] * wa_x[i];

        g_xxxxx_0_yz_0[i] = 4.0 * g_xxx_0_yz_0[i] * fbe_0 - 4.0 * g_xxx_0_yz_1[i] * fz_be_0 + g_xxxx_0_yz_1[i] * wa_x[i];

        g_xxxxx_0_zz_0[i] = 4.0 * g_xxx_0_zz_0[i] * fbe_0 - 4.0 * g_xxx_0_zz_1[i] * fz_be_0 + g_xxxx_0_zz_1[i] * wa_x[i];
    }

    /// Set up 6-12 components of targeted buffer : HSD

    auto g_xxxxy_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 6);

    auto g_xxxxy_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 7);

    auto g_xxxxy_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 8);

    auto g_xxxxy_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 9);

    auto g_xxxxy_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 10);

    auto g_xxxxy_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 11);

    #pragma omp simd aligned(g_xxxx_0_x_1, g_xxxx_0_xx_1, g_xxxx_0_xy_1, g_xxxx_0_xz_1, g_xxxx_0_y_1, g_xxxx_0_yy_1, g_xxxx_0_yz_1, g_xxxx_0_z_1, g_xxxx_0_zz_1, g_xxxxy_0_xx_0, g_xxxxy_0_xy_0, g_xxxxy_0_xz_0, g_xxxxy_0_yy_0, g_xxxxy_0_yz_0, g_xxxxy_0_zz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxy_0_xx_0[i] = g_xxxx_0_xx_1[i] * wa_y[i];

        g_xxxxy_0_xy_0[i] = g_xxxx_0_x_1[i] * fi_acd_0 + g_xxxx_0_xy_1[i] * wa_y[i];

        g_xxxxy_0_xz_0[i] = g_xxxx_0_xz_1[i] * wa_y[i];

        g_xxxxy_0_yy_0[i] = 2.0 * g_xxxx_0_y_1[i] * fi_acd_0 + g_xxxx_0_yy_1[i] * wa_y[i];

        g_xxxxy_0_yz_0[i] = g_xxxx_0_z_1[i] * fi_acd_0 + g_xxxx_0_yz_1[i] * wa_y[i];

        g_xxxxy_0_zz_0[i] = g_xxxx_0_zz_1[i] * wa_y[i];
    }

    /// Set up 12-18 components of targeted buffer : HSD

    auto g_xxxxz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 12);

    auto g_xxxxz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 13);

    auto g_xxxxz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 14);

    auto g_xxxxz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 15);

    auto g_xxxxz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 16);

    auto g_xxxxz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 17);

    #pragma omp simd aligned(g_xxxx_0_x_1, g_xxxx_0_xx_1, g_xxxx_0_xy_1, g_xxxx_0_xz_1, g_xxxx_0_y_1, g_xxxx_0_yy_1, g_xxxx_0_yz_1, g_xxxx_0_z_1, g_xxxx_0_zz_1, g_xxxxz_0_xx_0, g_xxxxz_0_xy_0, g_xxxxz_0_xz_0, g_xxxxz_0_yy_0, g_xxxxz_0_yz_0, g_xxxxz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxxz_0_xx_0[i] = g_xxxx_0_xx_1[i] * wa_z[i];

        g_xxxxz_0_xy_0[i] = g_xxxx_0_xy_1[i] * wa_z[i];

        g_xxxxz_0_xz_0[i] = g_xxxx_0_x_1[i] * fi_acd_0 + g_xxxx_0_xz_1[i] * wa_z[i];

        g_xxxxz_0_yy_0[i] = g_xxxx_0_yy_1[i] * wa_z[i];

        g_xxxxz_0_yz_0[i] = g_xxxx_0_y_1[i] * fi_acd_0 + g_xxxx_0_yz_1[i] * wa_z[i];

        g_xxxxz_0_zz_0[i] = 2.0 * g_xxxx_0_z_1[i] * fi_acd_0 + g_xxxx_0_zz_1[i] * wa_z[i];
    }

    /// Set up 18-24 components of targeted buffer : HSD

    auto g_xxxyy_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 18);

    auto g_xxxyy_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 19);

    auto g_xxxyy_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 20);

    auto g_xxxyy_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 21);

    auto g_xxxyy_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 22);

    auto g_xxxyy_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 23);

    #pragma omp simd aligned(g_xxx_0_xx_0, g_xxx_0_xx_1, g_xxx_0_xz_0, g_xxx_0_xz_1, g_xxxy_0_xx_1, g_xxxy_0_xz_1, g_xxxyy_0_xx_0, g_xxxyy_0_xy_0, g_xxxyy_0_xz_0, g_xxxyy_0_yy_0, g_xxxyy_0_yz_0, g_xxxyy_0_zz_0, g_xxyy_0_xy_1, g_xxyy_0_y_1, g_xxyy_0_yy_1, g_xxyy_0_yz_1, g_xxyy_0_zz_1, g_xyy_0_xy_0, g_xyy_0_xy_1, g_xyy_0_yy_0, g_xyy_0_yy_1, g_xyy_0_yz_0, g_xyy_0_yz_1, g_xyy_0_zz_0, g_xyy_0_zz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxyy_0_xx_0[i] = g_xxx_0_xx_0[i] * fbe_0 - g_xxx_0_xx_1[i] * fz_be_0 + g_xxxy_0_xx_1[i] * wa_y[i];

        g_xxxyy_0_xy_0[i] = 2.0 * g_xyy_0_xy_0[i] * fbe_0 - 2.0 * g_xyy_0_xy_1[i] * fz_be_0 + g_xxyy_0_y_1[i] * fi_acd_0 + g_xxyy_0_xy_1[i] * wa_x[i];

        g_xxxyy_0_xz_0[i] = g_xxx_0_xz_0[i] * fbe_0 - g_xxx_0_xz_1[i] * fz_be_0 + g_xxxy_0_xz_1[i] * wa_y[i];

        g_xxxyy_0_yy_0[i] = 2.0 * g_xyy_0_yy_0[i] * fbe_0 - 2.0 * g_xyy_0_yy_1[i] * fz_be_0 + g_xxyy_0_yy_1[i] * wa_x[i];

        g_xxxyy_0_yz_0[i] = 2.0 * g_xyy_0_yz_0[i] * fbe_0 - 2.0 * g_xyy_0_yz_1[i] * fz_be_0 + g_xxyy_0_yz_1[i] * wa_x[i];

        g_xxxyy_0_zz_0[i] = 2.0 * g_xyy_0_zz_0[i] * fbe_0 - 2.0 * g_xyy_0_zz_1[i] * fz_be_0 + g_xxyy_0_zz_1[i] * wa_x[i];
    }

    /// Set up 24-30 components of targeted buffer : HSD

    auto g_xxxyz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 24);

    auto g_xxxyz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 25);

    auto g_xxxyz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 26);

    auto g_xxxyz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 27);

    auto g_xxxyz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 28);

    auto g_xxxyz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 29);

    #pragma omp simd aligned(g_xxxy_0_xy_1, g_xxxy_0_yy_1, g_xxxyz_0_xx_0, g_xxxyz_0_xy_0, g_xxxyz_0_xz_0, g_xxxyz_0_yy_0, g_xxxyz_0_yz_0, g_xxxyz_0_zz_0, g_xxxz_0_xx_1, g_xxxz_0_xz_1, g_xxxz_0_yz_1, g_xxxz_0_z_1, g_xxxz_0_zz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxyz_0_xx_0[i] = g_xxxz_0_xx_1[i] * wa_y[i];

        g_xxxyz_0_xy_0[i] = g_xxxy_0_xy_1[i] * wa_z[i];

        g_xxxyz_0_xz_0[i] = g_xxxz_0_xz_1[i] * wa_y[i];

        g_xxxyz_0_yy_0[i] = g_xxxy_0_yy_1[i] * wa_z[i];

        g_xxxyz_0_yz_0[i] = g_xxxz_0_z_1[i] * fi_acd_0 + g_xxxz_0_yz_1[i] * wa_y[i];

        g_xxxyz_0_zz_0[i] = g_xxxz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 30-36 components of targeted buffer : HSD

    auto g_xxxzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 30);

    auto g_xxxzz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 31);

    auto g_xxxzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 32);

    auto g_xxxzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 33);

    auto g_xxxzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 34);

    auto g_xxxzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 35);

    #pragma omp simd aligned(g_xxx_0_xx_0, g_xxx_0_xx_1, g_xxx_0_xy_0, g_xxx_0_xy_1, g_xxxz_0_xx_1, g_xxxz_0_xy_1, g_xxxzz_0_xx_0, g_xxxzz_0_xy_0, g_xxxzz_0_xz_0, g_xxxzz_0_yy_0, g_xxxzz_0_yz_0, g_xxxzz_0_zz_0, g_xxzz_0_xz_1, g_xxzz_0_yy_1, g_xxzz_0_yz_1, g_xxzz_0_z_1, g_xxzz_0_zz_1, g_xzz_0_xz_0, g_xzz_0_xz_1, g_xzz_0_yy_0, g_xzz_0_yy_1, g_xzz_0_yz_0, g_xzz_0_yz_1, g_xzz_0_zz_0, g_xzz_0_zz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxzz_0_xx_0[i] = g_xxx_0_xx_0[i] * fbe_0 - g_xxx_0_xx_1[i] * fz_be_0 + g_xxxz_0_xx_1[i] * wa_z[i];

        g_xxxzz_0_xy_0[i] = g_xxx_0_xy_0[i] * fbe_0 - g_xxx_0_xy_1[i] * fz_be_0 + g_xxxz_0_xy_1[i] * wa_z[i];

        g_xxxzz_0_xz_0[i] = 2.0 * g_xzz_0_xz_0[i] * fbe_0 - 2.0 * g_xzz_0_xz_1[i] * fz_be_0 + g_xxzz_0_z_1[i] * fi_acd_0 + g_xxzz_0_xz_1[i] * wa_x[i];

        g_xxxzz_0_yy_0[i] = 2.0 * g_xzz_0_yy_0[i] * fbe_0 - 2.0 * g_xzz_0_yy_1[i] * fz_be_0 + g_xxzz_0_yy_1[i] * wa_x[i];

        g_xxxzz_0_yz_0[i] = 2.0 * g_xzz_0_yz_0[i] * fbe_0 - 2.0 * g_xzz_0_yz_1[i] * fz_be_0 + g_xxzz_0_yz_1[i] * wa_x[i];

        g_xxxzz_0_zz_0[i] = 2.0 * g_xzz_0_zz_0[i] * fbe_0 - 2.0 * g_xzz_0_zz_1[i] * fz_be_0 + g_xxzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 36-42 components of targeted buffer : HSD

    auto g_xxyyy_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 36);

    auto g_xxyyy_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 37);

    auto g_xxyyy_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 38);

    auto g_xxyyy_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 39);

    auto g_xxyyy_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 40);

    auto g_xxyyy_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 41);

    #pragma omp simd aligned(g_xxy_0_xx_0, g_xxy_0_xx_1, g_xxy_0_xz_0, g_xxy_0_xz_1, g_xxyy_0_xx_1, g_xxyy_0_xz_1, g_xxyyy_0_xx_0, g_xxyyy_0_xy_0, g_xxyyy_0_xz_0, g_xxyyy_0_yy_0, g_xxyyy_0_yz_0, g_xxyyy_0_zz_0, g_xyyy_0_xy_1, g_xyyy_0_y_1, g_xyyy_0_yy_1, g_xyyy_0_yz_1, g_xyyy_0_zz_1, g_yyy_0_xy_0, g_yyy_0_xy_1, g_yyy_0_yy_0, g_yyy_0_yy_1, g_yyy_0_yz_0, g_yyy_0_yz_1, g_yyy_0_zz_0, g_yyy_0_zz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyyy_0_xx_0[i] = 2.0 * g_xxy_0_xx_0[i] * fbe_0 - 2.0 * g_xxy_0_xx_1[i] * fz_be_0 + g_xxyy_0_xx_1[i] * wa_y[i];

        g_xxyyy_0_xy_0[i] = g_yyy_0_xy_0[i] * fbe_0 - g_yyy_0_xy_1[i] * fz_be_0 + g_xyyy_0_y_1[i] * fi_acd_0 + g_xyyy_0_xy_1[i] * wa_x[i];

        g_xxyyy_0_xz_0[i] = 2.0 * g_xxy_0_xz_0[i] * fbe_0 - 2.0 * g_xxy_0_xz_1[i] * fz_be_0 + g_xxyy_0_xz_1[i] * wa_y[i];

        g_xxyyy_0_yy_0[i] = g_yyy_0_yy_0[i] * fbe_0 - g_yyy_0_yy_1[i] * fz_be_0 + g_xyyy_0_yy_1[i] * wa_x[i];

        g_xxyyy_0_yz_0[i] = g_yyy_0_yz_0[i] * fbe_0 - g_yyy_0_yz_1[i] * fz_be_0 + g_xyyy_0_yz_1[i] * wa_x[i];

        g_xxyyy_0_zz_0[i] = g_yyy_0_zz_0[i] * fbe_0 - g_yyy_0_zz_1[i] * fz_be_0 + g_xyyy_0_zz_1[i] * wa_x[i];
    }

    /// Set up 42-48 components of targeted buffer : HSD

    auto g_xxyyz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 42);

    auto g_xxyyz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 43);

    auto g_xxyyz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 44);

    auto g_xxyyz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 45);

    auto g_xxyyz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 46);

    auto g_xxyyz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 47);

    #pragma omp simd aligned(g_xxyy_0_x_1, g_xxyy_0_xx_1, g_xxyy_0_xy_1, g_xxyy_0_xz_1, g_xxyy_0_y_1, g_xxyy_0_yy_1, g_xxyy_0_yz_1, g_xxyy_0_z_1, g_xxyy_0_zz_1, g_xxyyz_0_xx_0, g_xxyyz_0_xy_0, g_xxyyz_0_xz_0, g_xxyyz_0_yy_0, g_xxyyz_0_yz_0, g_xxyyz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyyz_0_xx_0[i] = g_xxyy_0_xx_1[i] * wa_z[i];

        g_xxyyz_0_xy_0[i] = g_xxyy_0_xy_1[i] * wa_z[i];

        g_xxyyz_0_xz_0[i] = g_xxyy_0_x_1[i] * fi_acd_0 + g_xxyy_0_xz_1[i] * wa_z[i];

        g_xxyyz_0_yy_0[i] = g_xxyy_0_yy_1[i] * wa_z[i];

        g_xxyyz_0_yz_0[i] = g_xxyy_0_y_1[i] * fi_acd_0 + g_xxyy_0_yz_1[i] * wa_z[i];

        g_xxyyz_0_zz_0[i] = 2.0 * g_xxyy_0_z_1[i] * fi_acd_0 + g_xxyy_0_zz_1[i] * wa_z[i];
    }

    /// Set up 48-54 components of targeted buffer : HSD

    auto g_xxyzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 48);

    auto g_xxyzz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 49);

    auto g_xxyzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 50);

    auto g_xxyzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 51);

    auto g_xxyzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 52);

    auto g_xxyzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 53);

    #pragma omp simd aligned(g_xxyzz_0_xx_0, g_xxyzz_0_xy_0, g_xxyzz_0_xz_0, g_xxyzz_0_yy_0, g_xxyzz_0_yz_0, g_xxyzz_0_zz_0, g_xxzz_0_x_1, g_xxzz_0_xx_1, g_xxzz_0_xy_1, g_xxzz_0_xz_1, g_xxzz_0_y_1, g_xxzz_0_yy_1, g_xxzz_0_yz_1, g_xxzz_0_z_1, g_xxzz_0_zz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyzz_0_xx_0[i] = g_xxzz_0_xx_1[i] * wa_y[i];

        g_xxyzz_0_xy_0[i] = g_xxzz_0_x_1[i] * fi_acd_0 + g_xxzz_0_xy_1[i] * wa_y[i];

        g_xxyzz_0_xz_0[i] = g_xxzz_0_xz_1[i] * wa_y[i];

        g_xxyzz_0_yy_0[i] = 2.0 * g_xxzz_0_y_1[i] * fi_acd_0 + g_xxzz_0_yy_1[i] * wa_y[i];

        g_xxyzz_0_yz_0[i] = g_xxzz_0_z_1[i] * fi_acd_0 + g_xxzz_0_yz_1[i] * wa_y[i];

        g_xxyzz_0_zz_0[i] = g_xxzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 54-60 components of targeted buffer : HSD

    auto g_xxzzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 54);

    auto g_xxzzz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 55);

    auto g_xxzzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 56);

    auto g_xxzzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 57);

    auto g_xxzzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 58);

    auto g_xxzzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 59);

    #pragma omp simd aligned(g_xxz_0_xx_0, g_xxz_0_xx_1, g_xxz_0_xy_0, g_xxz_0_xy_1, g_xxzz_0_xx_1, g_xxzz_0_xy_1, g_xxzzz_0_xx_0, g_xxzzz_0_xy_0, g_xxzzz_0_xz_0, g_xxzzz_0_yy_0, g_xxzzz_0_yz_0, g_xxzzz_0_zz_0, g_xzzz_0_xz_1, g_xzzz_0_yy_1, g_xzzz_0_yz_1, g_xzzz_0_z_1, g_xzzz_0_zz_1, g_zzz_0_xz_0, g_zzz_0_xz_1, g_zzz_0_yy_0, g_zzz_0_yy_1, g_zzz_0_yz_0, g_zzz_0_yz_1, g_zzz_0_zz_0, g_zzz_0_zz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzzz_0_xx_0[i] = 2.0 * g_xxz_0_xx_0[i] * fbe_0 - 2.0 * g_xxz_0_xx_1[i] * fz_be_0 + g_xxzz_0_xx_1[i] * wa_z[i];

        g_xxzzz_0_xy_0[i] = 2.0 * g_xxz_0_xy_0[i] * fbe_0 - 2.0 * g_xxz_0_xy_1[i] * fz_be_0 + g_xxzz_0_xy_1[i] * wa_z[i];

        g_xxzzz_0_xz_0[i] = g_zzz_0_xz_0[i] * fbe_0 - g_zzz_0_xz_1[i] * fz_be_0 + g_xzzz_0_z_1[i] * fi_acd_0 + g_xzzz_0_xz_1[i] * wa_x[i];

        g_xxzzz_0_yy_0[i] = g_zzz_0_yy_0[i] * fbe_0 - g_zzz_0_yy_1[i] * fz_be_0 + g_xzzz_0_yy_1[i] * wa_x[i];

        g_xxzzz_0_yz_0[i] = g_zzz_0_yz_0[i] * fbe_0 - g_zzz_0_yz_1[i] * fz_be_0 + g_xzzz_0_yz_1[i] * wa_x[i];

        g_xxzzz_0_zz_0[i] = g_zzz_0_zz_0[i] * fbe_0 - g_zzz_0_zz_1[i] * fz_be_0 + g_xzzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 60-66 components of targeted buffer : HSD

    auto g_xyyyy_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 60);

    auto g_xyyyy_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 61);

    auto g_xyyyy_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 62);

    auto g_xyyyy_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 63);

    auto g_xyyyy_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 64);

    auto g_xyyyy_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 65);

    #pragma omp simd aligned(g_xyyyy_0_xx_0, g_xyyyy_0_xy_0, g_xyyyy_0_xz_0, g_xyyyy_0_yy_0, g_xyyyy_0_yz_0, g_xyyyy_0_zz_0, g_yyyy_0_x_1, g_yyyy_0_xx_1, g_yyyy_0_xy_1, g_yyyy_0_xz_1, g_yyyy_0_y_1, g_yyyy_0_yy_1, g_yyyy_0_yz_1, g_yyyy_0_z_1, g_yyyy_0_zz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyy_0_xx_0[i] = 2.0 * g_yyyy_0_x_1[i] * fi_acd_0 + g_yyyy_0_xx_1[i] * wa_x[i];

        g_xyyyy_0_xy_0[i] = g_yyyy_0_y_1[i] * fi_acd_0 + g_yyyy_0_xy_1[i] * wa_x[i];

        g_xyyyy_0_xz_0[i] = g_yyyy_0_z_1[i] * fi_acd_0 + g_yyyy_0_xz_1[i] * wa_x[i];

        g_xyyyy_0_yy_0[i] = g_yyyy_0_yy_1[i] * wa_x[i];

        g_xyyyy_0_yz_0[i] = g_yyyy_0_yz_1[i] * wa_x[i];

        g_xyyyy_0_zz_0[i] = g_yyyy_0_zz_1[i] * wa_x[i];
    }

    /// Set up 66-72 components of targeted buffer : HSD

    auto g_xyyyz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 66);

    auto g_xyyyz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 67);

    auto g_xyyyz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 68);

    auto g_xyyyz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 69);

    auto g_xyyyz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 70);

    auto g_xyyyz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 71);

    #pragma omp simd aligned(g_xyyy_0_xx_1, g_xyyy_0_xy_1, g_xyyyz_0_xx_0, g_xyyyz_0_xy_0, g_xyyyz_0_xz_0, g_xyyyz_0_yy_0, g_xyyyz_0_yz_0, g_xyyyz_0_zz_0, g_yyyz_0_xz_1, g_yyyz_0_yy_1, g_yyyz_0_yz_1, g_yyyz_0_z_1, g_yyyz_0_zz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyyz_0_xx_0[i] = g_xyyy_0_xx_1[i] * wa_z[i];

        g_xyyyz_0_xy_0[i] = g_xyyy_0_xy_1[i] * wa_z[i];

        g_xyyyz_0_xz_0[i] = g_yyyz_0_z_1[i] * fi_acd_0 + g_yyyz_0_xz_1[i] * wa_x[i];

        g_xyyyz_0_yy_0[i] = g_yyyz_0_yy_1[i] * wa_x[i];

        g_xyyyz_0_yz_0[i] = g_yyyz_0_yz_1[i] * wa_x[i];

        g_xyyyz_0_zz_0[i] = g_yyyz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 72-78 components of targeted buffer : HSD

    auto g_xyyzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 72);

    auto g_xyyzz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 73);

    auto g_xyyzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 74);

    auto g_xyyzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 75);

    auto g_xyyzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 76);

    auto g_xyyzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 77);

    #pragma omp simd aligned(g_xyyzz_0_xx_0, g_xyyzz_0_xy_0, g_xyyzz_0_xz_0, g_xyyzz_0_yy_0, g_xyyzz_0_yz_0, g_xyyzz_0_zz_0, g_yyzz_0_x_1, g_yyzz_0_xx_1, g_yyzz_0_xy_1, g_yyzz_0_xz_1, g_yyzz_0_y_1, g_yyzz_0_yy_1, g_yyzz_0_yz_1, g_yyzz_0_z_1, g_yyzz_0_zz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyzz_0_xx_0[i] = 2.0 * g_yyzz_0_x_1[i] * fi_acd_0 + g_yyzz_0_xx_1[i] * wa_x[i];

        g_xyyzz_0_xy_0[i] = g_yyzz_0_y_1[i] * fi_acd_0 + g_yyzz_0_xy_1[i] * wa_x[i];

        g_xyyzz_0_xz_0[i] = g_yyzz_0_z_1[i] * fi_acd_0 + g_yyzz_0_xz_1[i] * wa_x[i];

        g_xyyzz_0_yy_0[i] = g_yyzz_0_yy_1[i] * wa_x[i];

        g_xyyzz_0_yz_0[i] = g_yyzz_0_yz_1[i] * wa_x[i];

        g_xyyzz_0_zz_0[i] = g_yyzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 78-84 components of targeted buffer : HSD

    auto g_xyzzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 78);

    auto g_xyzzz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 79);

    auto g_xyzzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 80);

    auto g_xyzzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 81);

    auto g_xyzzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 82);

    auto g_xyzzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 83);

    #pragma omp simd aligned(g_xyzzz_0_xx_0, g_xyzzz_0_xy_0, g_xyzzz_0_xz_0, g_xyzzz_0_yy_0, g_xyzzz_0_yz_0, g_xyzzz_0_zz_0, g_xzzz_0_xx_1, g_xzzz_0_xz_1, g_yzzz_0_xy_1, g_yzzz_0_y_1, g_yzzz_0_yy_1, g_yzzz_0_yz_1, g_yzzz_0_zz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzzz_0_xx_0[i] = g_xzzz_0_xx_1[i] * wa_y[i];

        g_xyzzz_0_xy_0[i] = g_yzzz_0_y_1[i] * fi_acd_0 + g_yzzz_0_xy_1[i] * wa_x[i];

        g_xyzzz_0_xz_0[i] = g_xzzz_0_xz_1[i] * wa_y[i];

        g_xyzzz_0_yy_0[i] = g_yzzz_0_yy_1[i] * wa_x[i];

        g_xyzzz_0_yz_0[i] = g_yzzz_0_yz_1[i] * wa_x[i];

        g_xyzzz_0_zz_0[i] = g_yzzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 84-90 components of targeted buffer : HSD

    auto g_xzzzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 84);

    auto g_xzzzz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 85);

    auto g_xzzzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 86);

    auto g_xzzzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 87);

    auto g_xzzzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 88);

    auto g_xzzzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 89);

    #pragma omp simd aligned(g_xzzzz_0_xx_0, g_xzzzz_0_xy_0, g_xzzzz_0_xz_0, g_xzzzz_0_yy_0, g_xzzzz_0_yz_0, g_xzzzz_0_zz_0, g_zzzz_0_x_1, g_zzzz_0_xx_1, g_zzzz_0_xy_1, g_zzzz_0_xz_1, g_zzzz_0_y_1, g_zzzz_0_yy_1, g_zzzz_0_yz_1, g_zzzz_0_z_1, g_zzzz_0_zz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzzz_0_xx_0[i] = 2.0 * g_zzzz_0_x_1[i] * fi_acd_0 + g_zzzz_0_xx_1[i] * wa_x[i];

        g_xzzzz_0_xy_0[i] = g_zzzz_0_y_1[i] * fi_acd_0 + g_zzzz_0_xy_1[i] * wa_x[i];

        g_xzzzz_0_xz_0[i] = g_zzzz_0_z_1[i] * fi_acd_0 + g_zzzz_0_xz_1[i] * wa_x[i];

        g_xzzzz_0_yy_0[i] = g_zzzz_0_yy_1[i] * wa_x[i];

        g_xzzzz_0_yz_0[i] = g_zzzz_0_yz_1[i] * wa_x[i];

        g_xzzzz_0_zz_0[i] = g_zzzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 90-96 components of targeted buffer : HSD

    auto g_yyyyy_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 90);

    auto g_yyyyy_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 91);

    auto g_yyyyy_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 92);

    auto g_yyyyy_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 93);

    auto g_yyyyy_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 94);

    auto g_yyyyy_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 95);

    #pragma omp simd aligned(g_yyy_0_xx_0, g_yyy_0_xx_1, g_yyy_0_xy_0, g_yyy_0_xy_1, g_yyy_0_xz_0, g_yyy_0_xz_1, g_yyy_0_yy_0, g_yyy_0_yy_1, g_yyy_0_yz_0, g_yyy_0_yz_1, g_yyy_0_zz_0, g_yyy_0_zz_1, g_yyyy_0_x_1, g_yyyy_0_xx_1, g_yyyy_0_xy_1, g_yyyy_0_xz_1, g_yyyy_0_y_1, g_yyyy_0_yy_1, g_yyyy_0_yz_1, g_yyyy_0_z_1, g_yyyy_0_zz_1, g_yyyyy_0_xx_0, g_yyyyy_0_xy_0, g_yyyyy_0_xz_0, g_yyyyy_0_yy_0, g_yyyyy_0_yz_0, g_yyyyy_0_zz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyyy_0_xx_0[i] = 4.0 * g_yyy_0_xx_0[i] * fbe_0 - 4.0 * g_yyy_0_xx_1[i] * fz_be_0 + g_yyyy_0_xx_1[i] * wa_y[i];

        g_yyyyy_0_xy_0[i] = 4.0 * g_yyy_0_xy_0[i] * fbe_0 - 4.0 * g_yyy_0_xy_1[i] * fz_be_0 + g_yyyy_0_x_1[i] * fi_acd_0 + g_yyyy_0_xy_1[i] * wa_y[i];

        g_yyyyy_0_xz_0[i] = 4.0 * g_yyy_0_xz_0[i] * fbe_0 - 4.0 * g_yyy_0_xz_1[i] * fz_be_0 + g_yyyy_0_xz_1[i] * wa_y[i];

        g_yyyyy_0_yy_0[i] = 4.0 * g_yyy_0_yy_0[i] * fbe_0 - 4.0 * g_yyy_0_yy_1[i] * fz_be_0 + 2.0 * g_yyyy_0_y_1[i] * fi_acd_0 + g_yyyy_0_yy_1[i] * wa_y[i];

        g_yyyyy_0_yz_0[i] = 4.0 * g_yyy_0_yz_0[i] * fbe_0 - 4.0 * g_yyy_0_yz_1[i] * fz_be_0 + g_yyyy_0_z_1[i] * fi_acd_0 + g_yyyy_0_yz_1[i] * wa_y[i];

        g_yyyyy_0_zz_0[i] = 4.0 * g_yyy_0_zz_0[i] * fbe_0 - 4.0 * g_yyy_0_zz_1[i] * fz_be_0 + g_yyyy_0_zz_1[i] * wa_y[i];
    }

    /// Set up 96-102 components of targeted buffer : HSD

    auto g_yyyyz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 96);

    auto g_yyyyz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 97);

    auto g_yyyyz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 98);

    auto g_yyyyz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 99);

    auto g_yyyyz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 100);

    auto g_yyyyz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 101);

    #pragma omp simd aligned(g_yyyy_0_x_1, g_yyyy_0_xx_1, g_yyyy_0_xy_1, g_yyyy_0_xz_1, g_yyyy_0_y_1, g_yyyy_0_yy_1, g_yyyy_0_yz_1, g_yyyy_0_z_1, g_yyyy_0_zz_1, g_yyyyz_0_xx_0, g_yyyyz_0_xy_0, g_yyyyz_0_xz_0, g_yyyyz_0_yy_0, g_yyyyz_0_yz_0, g_yyyyz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyyz_0_xx_0[i] = g_yyyy_0_xx_1[i] * wa_z[i];

        g_yyyyz_0_xy_0[i] = g_yyyy_0_xy_1[i] * wa_z[i];

        g_yyyyz_0_xz_0[i] = g_yyyy_0_x_1[i] * fi_acd_0 + g_yyyy_0_xz_1[i] * wa_z[i];

        g_yyyyz_0_yy_0[i] = g_yyyy_0_yy_1[i] * wa_z[i];

        g_yyyyz_0_yz_0[i] = g_yyyy_0_y_1[i] * fi_acd_0 + g_yyyy_0_yz_1[i] * wa_z[i];

        g_yyyyz_0_zz_0[i] = 2.0 * g_yyyy_0_z_1[i] * fi_acd_0 + g_yyyy_0_zz_1[i] * wa_z[i];
    }

    /// Set up 102-108 components of targeted buffer : HSD

    auto g_yyyzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 102);

    auto g_yyyzz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 103);

    auto g_yyyzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 104);

    auto g_yyyzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 105);

    auto g_yyyzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 106);

    auto g_yyyzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 107);

    #pragma omp simd aligned(g_yyy_0_xy_0, g_yyy_0_xy_1, g_yyy_0_yy_0, g_yyy_0_yy_1, g_yyyz_0_xy_1, g_yyyz_0_yy_1, g_yyyzz_0_xx_0, g_yyyzz_0_xy_0, g_yyyzz_0_xz_0, g_yyyzz_0_yy_0, g_yyyzz_0_yz_0, g_yyyzz_0_zz_0, g_yyzz_0_xx_1, g_yyzz_0_xz_1, g_yyzz_0_yz_1, g_yyzz_0_z_1, g_yyzz_0_zz_1, g_yzz_0_xx_0, g_yzz_0_xx_1, g_yzz_0_xz_0, g_yzz_0_xz_1, g_yzz_0_yz_0, g_yzz_0_yz_1, g_yzz_0_zz_0, g_yzz_0_zz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyzz_0_xx_0[i] = 2.0 * g_yzz_0_xx_0[i] * fbe_0 - 2.0 * g_yzz_0_xx_1[i] * fz_be_0 + g_yyzz_0_xx_1[i] * wa_y[i];

        g_yyyzz_0_xy_0[i] = g_yyy_0_xy_0[i] * fbe_0 - g_yyy_0_xy_1[i] * fz_be_0 + g_yyyz_0_xy_1[i] * wa_z[i];

        g_yyyzz_0_xz_0[i] = 2.0 * g_yzz_0_xz_0[i] * fbe_0 - 2.0 * g_yzz_0_xz_1[i] * fz_be_0 + g_yyzz_0_xz_1[i] * wa_y[i];

        g_yyyzz_0_yy_0[i] = g_yyy_0_yy_0[i] * fbe_0 - g_yyy_0_yy_1[i] * fz_be_0 + g_yyyz_0_yy_1[i] * wa_z[i];

        g_yyyzz_0_yz_0[i] = 2.0 * g_yzz_0_yz_0[i] * fbe_0 - 2.0 * g_yzz_0_yz_1[i] * fz_be_0 + g_yyzz_0_z_1[i] * fi_acd_0 + g_yyzz_0_yz_1[i] * wa_y[i];

        g_yyyzz_0_zz_0[i] = 2.0 * g_yzz_0_zz_0[i] * fbe_0 - 2.0 * g_yzz_0_zz_1[i] * fz_be_0 + g_yyzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 108-114 components of targeted buffer : HSD

    auto g_yyzzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 108);

    auto g_yyzzz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 109);

    auto g_yyzzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 110);

    auto g_yyzzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 111);

    auto g_yyzzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 112);

    auto g_yyzzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 113);

    #pragma omp simd aligned(g_yyz_0_xy_0, g_yyz_0_xy_1, g_yyz_0_yy_0, g_yyz_0_yy_1, g_yyzz_0_xy_1, g_yyzz_0_yy_1, g_yyzzz_0_xx_0, g_yyzzz_0_xy_0, g_yyzzz_0_xz_0, g_yyzzz_0_yy_0, g_yyzzz_0_yz_0, g_yyzzz_0_zz_0, g_yzzz_0_xx_1, g_yzzz_0_xz_1, g_yzzz_0_yz_1, g_yzzz_0_z_1, g_yzzz_0_zz_1, g_zzz_0_xx_0, g_zzz_0_xx_1, g_zzz_0_xz_0, g_zzz_0_xz_1, g_zzz_0_yz_0, g_zzz_0_yz_1, g_zzz_0_zz_0, g_zzz_0_zz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzzz_0_xx_0[i] = g_zzz_0_xx_0[i] * fbe_0 - g_zzz_0_xx_1[i] * fz_be_0 + g_yzzz_0_xx_1[i] * wa_y[i];

        g_yyzzz_0_xy_0[i] = 2.0 * g_yyz_0_xy_0[i] * fbe_0 - 2.0 * g_yyz_0_xy_1[i] * fz_be_0 + g_yyzz_0_xy_1[i] * wa_z[i];

        g_yyzzz_0_xz_0[i] = g_zzz_0_xz_0[i] * fbe_0 - g_zzz_0_xz_1[i] * fz_be_0 + g_yzzz_0_xz_1[i] * wa_y[i];

        g_yyzzz_0_yy_0[i] = 2.0 * g_yyz_0_yy_0[i] * fbe_0 - 2.0 * g_yyz_0_yy_1[i] * fz_be_0 + g_yyzz_0_yy_1[i] * wa_z[i];

        g_yyzzz_0_yz_0[i] = g_zzz_0_yz_0[i] * fbe_0 - g_zzz_0_yz_1[i] * fz_be_0 + g_yzzz_0_z_1[i] * fi_acd_0 + g_yzzz_0_yz_1[i] * wa_y[i];

        g_yyzzz_0_zz_0[i] = g_zzz_0_zz_0[i] * fbe_0 - g_zzz_0_zz_1[i] * fz_be_0 + g_yzzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 114-120 components of targeted buffer : HSD

    auto g_yzzzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 114);

    auto g_yzzzz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 115);

    auto g_yzzzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 116);

    auto g_yzzzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 117);

    auto g_yzzzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 118);

    auto g_yzzzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 119);

    #pragma omp simd aligned(g_yzzzz_0_xx_0, g_yzzzz_0_xy_0, g_yzzzz_0_xz_0, g_yzzzz_0_yy_0, g_yzzzz_0_yz_0, g_yzzzz_0_zz_0, g_zzzz_0_x_1, g_zzzz_0_xx_1, g_zzzz_0_xy_1, g_zzzz_0_xz_1, g_zzzz_0_y_1, g_zzzz_0_yy_1, g_zzzz_0_yz_1, g_zzzz_0_z_1, g_zzzz_0_zz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzzz_0_xx_0[i] = g_zzzz_0_xx_1[i] * wa_y[i];

        g_yzzzz_0_xy_0[i] = g_zzzz_0_x_1[i] * fi_acd_0 + g_zzzz_0_xy_1[i] * wa_y[i];

        g_yzzzz_0_xz_0[i] = g_zzzz_0_xz_1[i] * wa_y[i];

        g_yzzzz_0_yy_0[i] = 2.0 * g_zzzz_0_y_1[i] * fi_acd_0 + g_zzzz_0_yy_1[i] * wa_y[i];

        g_yzzzz_0_yz_0[i] = g_zzzz_0_z_1[i] * fi_acd_0 + g_zzzz_0_yz_1[i] * wa_y[i];

        g_yzzzz_0_zz_0[i] = g_zzzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 120-126 components of targeted buffer : HSD

    auto g_zzzzz_0_xx_0 = pbuffer.data(idx_eri_0_hsd + 120);

    auto g_zzzzz_0_xy_0 = pbuffer.data(idx_eri_0_hsd + 121);

    auto g_zzzzz_0_xz_0 = pbuffer.data(idx_eri_0_hsd + 122);

    auto g_zzzzz_0_yy_0 = pbuffer.data(idx_eri_0_hsd + 123);

    auto g_zzzzz_0_yz_0 = pbuffer.data(idx_eri_0_hsd + 124);

    auto g_zzzzz_0_zz_0 = pbuffer.data(idx_eri_0_hsd + 125);

    #pragma omp simd aligned(g_zzz_0_xx_0, g_zzz_0_xx_1, g_zzz_0_xy_0, g_zzz_0_xy_1, g_zzz_0_xz_0, g_zzz_0_xz_1, g_zzz_0_yy_0, g_zzz_0_yy_1, g_zzz_0_yz_0, g_zzz_0_yz_1, g_zzz_0_zz_0, g_zzz_0_zz_1, g_zzzz_0_x_1, g_zzzz_0_xx_1, g_zzzz_0_xy_1, g_zzzz_0_xz_1, g_zzzz_0_y_1, g_zzzz_0_yy_1, g_zzzz_0_yz_1, g_zzzz_0_z_1, g_zzzz_0_zz_1, g_zzzzz_0_xx_0, g_zzzzz_0_xy_0, g_zzzzz_0_xz_0, g_zzzzz_0_yy_0, g_zzzzz_0_yz_0, g_zzzzz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzzz_0_xx_0[i] = 4.0 * g_zzz_0_xx_0[i] * fbe_0 - 4.0 * g_zzz_0_xx_1[i] * fz_be_0 + g_zzzz_0_xx_1[i] * wa_z[i];

        g_zzzzz_0_xy_0[i] = 4.0 * g_zzz_0_xy_0[i] * fbe_0 - 4.0 * g_zzz_0_xy_1[i] * fz_be_0 + g_zzzz_0_xy_1[i] * wa_z[i];

        g_zzzzz_0_xz_0[i] = 4.0 * g_zzz_0_xz_0[i] * fbe_0 - 4.0 * g_zzz_0_xz_1[i] * fz_be_0 + g_zzzz_0_x_1[i] * fi_acd_0 + g_zzzz_0_xz_1[i] * wa_z[i];

        g_zzzzz_0_yy_0[i] = 4.0 * g_zzz_0_yy_0[i] * fbe_0 - 4.0 * g_zzz_0_yy_1[i] * fz_be_0 + g_zzzz_0_yy_1[i] * wa_z[i];

        g_zzzzz_0_yz_0[i] = 4.0 * g_zzz_0_yz_0[i] * fbe_0 - 4.0 * g_zzz_0_yz_1[i] * fz_be_0 + g_zzzz_0_y_1[i] * fi_acd_0 + g_zzzz_0_yz_1[i] * wa_z[i];

        g_zzzzz_0_zz_0[i] = 4.0 * g_zzz_0_zz_0[i] * fbe_0 - 4.0 * g_zzz_0_zz_1[i] * fz_be_0 + 2.0 * g_zzzz_0_z_1[i] * fi_acd_0 + g_zzzz_0_zz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

