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

#include "TwoCenterElectronRepulsionPrimRecFG.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_fg(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_fg,
                                const size_t idx_eri_0_pg,
                                const size_t idx_eri_1_pg,
                                const size_t idx_eri_1_df,
                                const size_t idx_eri_1_dg,
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

    // Set up components of auxiliary buffer : PG

    auto g_x_xxxx_0 = pbuffer.data(idx_eri_0_pg);

    auto g_x_xxxy_0 = pbuffer.data(idx_eri_0_pg + 1);

    auto g_x_xxxz_0 = pbuffer.data(idx_eri_0_pg + 2);

    auto g_x_xxyy_0 = pbuffer.data(idx_eri_0_pg + 3);

    auto g_x_xxyz_0 = pbuffer.data(idx_eri_0_pg + 4);

    auto g_x_xxzz_0 = pbuffer.data(idx_eri_0_pg + 5);

    auto g_x_xyyy_0 = pbuffer.data(idx_eri_0_pg + 6);

    auto g_x_xyyz_0 = pbuffer.data(idx_eri_0_pg + 7);

    auto g_x_xyzz_0 = pbuffer.data(idx_eri_0_pg + 8);

    auto g_x_xzzz_0 = pbuffer.data(idx_eri_0_pg + 9);

    auto g_x_yyyy_0 = pbuffer.data(idx_eri_0_pg + 10);

    auto g_x_yyyz_0 = pbuffer.data(idx_eri_0_pg + 11);

    auto g_x_yyzz_0 = pbuffer.data(idx_eri_0_pg + 12);

    auto g_x_yzzz_0 = pbuffer.data(idx_eri_0_pg + 13);

    auto g_x_zzzz_0 = pbuffer.data(idx_eri_0_pg + 14);

    auto g_y_xxxx_0 = pbuffer.data(idx_eri_0_pg + 15);

    auto g_y_xxxy_0 = pbuffer.data(idx_eri_0_pg + 16);

    auto g_y_xxxz_0 = pbuffer.data(idx_eri_0_pg + 17);

    auto g_y_xxyy_0 = pbuffer.data(idx_eri_0_pg + 18);

    auto g_y_xxyz_0 = pbuffer.data(idx_eri_0_pg + 19);

    auto g_y_xxzz_0 = pbuffer.data(idx_eri_0_pg + 20);

    auto g_y_xyyy_0 = pbuffer.data(idx_eri_0_pg + 21);

    auto g_y_xyyz_0 = pbuffer.data(idx_eri_0_pg + 22);

    auto g_y_xyzz_0 = pbuffer.data(idx_eri_0_pg + 23);

    auto g_y_xzzz_0 = pbuffer.data(idx_eri_0_pg + 24);

    auto g_y_yyyy_0 = pbuffer.data(idx_eri_0_pg + 25);

    auto g_y_yyyz_0 = pbuffer.data(idx_eri_0_pg + 26);

    auto g_y_yyzz_0 = pbuffer.data(idx_eri_0_pg + 27);

    auto g_y_yzzz_0 = pbuffer.data(idx_eri_0_pg + 28);

    auto g_y_zzzz_0 = pbuffer.data(idx_eri_0_pg + 29);

    auto g_z_xxxx_0 = pbuffer.data(idx_eri_0_pg + 30);

    auto g_z_xxxy_0 = pbuffer.data(idx_eri_0_pg + 31);

    auto g_z_xxxz_0 = pbuffer.data(idx_eri_0_pg + 32);

    auto g_z_xxyy_0 = pbuffer.data(idx_eri_0_pg + 33);

    auto g_z_xxyz_0 = pbuffer.data(idx_eri_0_pg + 34);

    auto g_z_xxzz_0 = pbuffer.data(idx_eri_0_pg + 35);

    auto g_z_xyyy_0 = pbuffer.data(idx_eri_0_pg + 36);

    auto g_z_xyyz_0 = pbuffer.data(idx_eri_0_pg + 37);

    auto g_z_xyzz_0 = pbuffer.data(idx_eri_0_pg + 38);

    auto g_z_xzzz_0 = pbuffer.data(idx_eri_0_pg + 39);

    auto g_z_yyyy_0 = pbuffer.data(idx_eri_0_pg + 40);

    auto g_z_yyyz_0 = pbuffer.data(idx_eri_0_pg + 41);

    auto g_z_yyzz_0 = pbuffer.data(idx_eri_0_pg + 42);

    auto g_z_yzzz_0 = pbuffer.data(idx_eri_0_pg + 43);

    auto g_z_zzzz_0 = pbuffer.data(idx_eri_0_pg + 44);

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

    // Set up components of auxiliary buffer : DF

    auto g_xx_xxx_1 = pbuffer.data(idx_eri_1_df);

    auto g_xx_xxy_1 = pbuffer.data(idx_eri_1_df + 1);

    auto g_xx_xxz_1 = pbuffer.data(idx_eri_1_df + 2);

    auto g_xx_xyy_1 = pbuffer.data(idx_eri_1_df + 3);

    auto g_xx_xyz_1 = pbuffer.data(idx_eri_1_df + 4);

    auto g_xx_xzz_1 = pbuffer.data(idx_eri_1_df + 5);

    auto g_xx_yyy_1 = pbuffer.data(idx_eri_1_df + 6);

    auto g_xx_yyz_1 = pbuffer.data(idx_eri_1_df + 7);

    auto g_xx_yzz_1 = pbuffer.data(idx_eri_1_df + 8);

    auto g_xx_zzz_1 = pbuffer.data(idx_eri_1_df + 9);

    auto g_yy_xxx_1 = pbuffer.data(idx_eri_1_df + 30);

    auto g_yy_xxy_1 = pbuffer.data(idx_eri_1_df + 31);

    auto g_yy_xxz_1 = pbuffer.data(idx_eri_1_df + 32);

    auto g_yy_xyy_1 = pbuffer.data(idx_eri_1_df + 33);

    auto g_yy_xyz_1 = pbuffer.data(idx_eri_1_df + 34);

    auto g_yy_xzz_1 = pbuffer.data(idx_eri_1_df + 35);

    auto g_yy_yyy_1 = pbuffer.data(idx_eri_1_df + 36);

    auto g_yy_yyz_1 = pbuffer.data(idx_eri_1_df + 37);

    auto g_yy_yzz_1 = pbuffer.data(idx_eri_1_df + 38);

    auto g_yy_zzz_1 = pbuffer.data(idx_eri_1_df + 39);

    auto g_yz_xyz_1 = pbuffer.data(idx_eri_1_df + 44);

    auto g_yz_yyz_1 = pbuffer.data(idx_eri_1_df + 47);

    auto g_yz_yzz_1 = pbuffer.data(idx_eri_1_df + 48);

    auto g_zz_xxx_1 = pbuffer.data(idx_eri_1_df + 50);

    auto g_zz_xxy_1 = pbuffer.data(idx_eri_1_df + 51);

    auto g_zz_xxz_1 = pbuffer.data(idx_eri_1_df + 52);

    auto g_zz_xyy_1 = pbuffer.data(idx_eri_1_df + 53);

    auto g_zz_xyz_1 = pbuffer.data(idx_eri_1_df + 54);

    auto g_zz_xzz_1 = pbuffer.data(idx_eri_1_df + 55);

    auto g_zz_yyy_1 = pbuffer.data(idx_eri_1_df + 56);

    auto g_zz_yyz_1 = pbuffer.data(idx_eri_1_df + 57);

    auto g_zz_yzz_1 = pbuffer.data(idx_eri_1_df + 58);

    auto g_zz_zzz_1 = pbuffer.data(idx_eri_1_df + 59);

    // Set up components of auxiliary buffer : DG

    auto g_xx_xxxx_1 = pbuffer.data(idx_eri_1_dg);

    auto g_xx_xxxy_1 = pbuffer.data(idx_eri_1_dg + 1);

    auto g_xx_xxxz_1 = pbuffer.data(idx_eri_1_dg + 2);

    auto g_xx_xxyy_1 = pbuffer.data(idx_eri_1_dg + 3);

    auto g_xx_xxyz_1 = pbuffer.data(idx_eri_1_dg + 4);

    auto g_xx_xxzz_1 = pbuffer.data(idx_eri_1_dg + 5);

    auto g_xx_xyyy_1 = pbuffer.data(idx_eri_1_dg + 6);

    auto g_xx_xyyz_1 = pbuffer.data(idx_eri_1_dg + 7);

    auto g_xx_xyzz_1 = pbuffer.data(idx_eri_1_dg + 8);

    auto g_xx_xzzz_1 = pbuffer.data(idx_eri_1_dg + 9);

    auto g_xx_yyyy_1 = pbuffer.data(idx_eri_1_dg + 10);

    auto g_xx_yyyz_1 = pbuffer.data(idx_eri_1_dg + 11);

    auto g_xx_yyzz_1 = pbuffer.data(idx_eri_1_dg + 12);

    auto g_xx_yzzz_1 = pbuffer.data(idx_eri_1_dg + 13);

    auto g_xx_zzzz_1 = pbuffer.data(idx_eri_1_dg + 14);

    auto g_xy_xxxy_1 = pbuffer.data(idx_eri_1_dg + 16);

    auto g_xy_xxyy_1 = pbuffer.data(idx_eri_1_dg + 18);

    auto g_xy_xyyy_1 = pbuffer.data(idx_eri_1_dg + 21);

    auto g_xz_xxxx_1 = pbuffer.data(idx_eri_1_dg + 30);

    auto g_xz_xxxz_1 = pbuffer.data(idx_eri_1_dg + 32);

    auto g_xz_xxzz_1 = pbuffer.data(idx_eri_1_dg + 35);

    auto g_xz_xzzz_1 = pbuffer.data(idx_eri_1_dg + 39);

    auto g_yy_xxxx_1 = pbuffer.data(idx_eri_1_dg + 45);

    auto g_yy_xxxy_1 = pbuffer.data(idx_eri_1_dg + 46);

    auto g_yy_xxxz_1 = pbuffer.data(idx_eri_1_dg + 47);

    auto g_yy_xxyy_1 = pbuffer.data(idx_eri_1_dg + 48);

    auto g_yy_xxyz_1 = pbuffer.data(idx_eri_1_dg + 49);

    auto g_yy_xxzz_1 = pbuffer.data(idx_eri_1_dg + 50);

    auto g_yy_xyyy_1 = pbuffer.data(idx_eri_1_dg + 51);

    auto g_yy_xyyz_1 = pbuffer.data(idx_eri_1_dg + 52);

    auto g_yy_xyzz_1 = pbuffer.data(idx_eri_1_dg + 53);

    auto g_yy_xzzz_1 = pbuffer.data(idx_eri_1_dg + 54);

    auto g_yy_yyyy_1 = pbuffer.data(idx_eri_1_dg + 55);

    auto g_yy_yyyz_1 = pbuffer.data(idx_eri_1_dg + 56);

    auto g_yy_yyzz_1 = pbuffer.data(idx_eri_1_dg + 57);

    auto g_yy_yzzz_1 = pbuffer.data(idx_eri_1_dg + 58);

    auto g_yy_zzzz_1 = pbuffer.data(idx_eri_1_dg + 59);

    auto g_yz_xxyz_1 = pbuffer.data(idx_eri_1_dg + 64);

    auto g_yz_xyyz_1 = pbuffer.data(idx_eri_1_dg + 67);

    auto g_yz_xyzz_1 = pbuffer.data(idx_eri_1_dg + 68);

    auto g_yz_yyyy_1 = pbuffer.data(idx_eri_1_dg + 70);

    auto g_yz_yyyz_1 = pbuffer.data(idx_eri_1_dg + 71);

    auto g_yz_yyzz_1 = pbuffer.data(idx_eri_1_dg + 72);

    auto g_yz_yzzz_1 = pbuffer.data(idx_eri_1_dg + 73);

    auto g_yz_zzzz_1 = pbuffer.data(idx_eri_1_dg + 74);

    auto g_zz_xxxx_1 = pbuffer.data(idx_eri_1_dg + 75);

    auto g_zz_xxxy_1 = pbuffer.data(idx_eri_1_dg + 76);

    auto g_zz_xxxz_1 = pbuffer.data(idx_eri_1_dg + 77);

    auto g_zz_xxyy_1 = pbuffer.data(idx_eri_1_dg + 78);

    auto g_zz_xxyz_1 = pbuffer.data(idx_eri_1_dg + 79);

    auto g_zz_xxzz_1 = pbuffer.data(idx_eri_1_dg + 80);

    auto g_zz_xyyy_1 = pbuffer.data(idx_eri_1_dg + 81);

    auto g_zz_xyyz_1 = pbuffer.data(idx_eri_1_dg + 82);

    auto g_zz_xyzz_1 = pbuffer.data(idx_eri_1_dg + 83);

    auto g_zz_xzzz_1 = pbuffer.data(idx_eri_1_dg + 84);

    auto g_zz_yyyy_1 = pbuffer.data(idx_eri_1_dg + 85);

    auto g_zz_yyyz_1 = pbuffer.data(idx_eri_1_dg + 86);

    auto g_zz_yyzz_1 = pbuffer.data(idx_eri_1_dg + 87);

    auto g_zz_yzzz_1 = pbuffer.data(idx_eri_1_dg + 88);

    auto g_zz_zzzz_1 = pbuffer.data(idx_eri_1_dg + 89);

    // Set up 0-15 components of targeted buffer : FG

    auto g_xxx_xxxx_0 = pbuffer.data(idx_eri_0_fg);

    auto g_xxx_xxxy_0 = pbuffer.data(idx_eri_0_fg + 1);

    auto g_xxx_xxxz_0 = pbuffer.data(idx_eri_0_fg + 2);

    auto g_xxx_xxyy_0 = pbuffer.data(idx_eri_0_fg + 3);

    auto g_xxx_xxyz_0 = pbuffer.data(idx_eri_0_fg + 4);

    auto g_xxx_xxzz_0 = pbuffer.data(idx_eri_0_fg + 5);

    auto g_xxx_xyyy_0 = pbuffer.data(idx_eri_0_fg + 6);

    auto g_xxx_xyyz_0 = pbuffer.data(idx_eri_0_fg + 7);

    auto g_xxx_xyzz_0 = pbuffer.data(idx_eri_0_fg + 8);

    auto g_xxx_xzzz_0 = pbuffer.data(idx_eri_0_fg + 9);

    auto g_xxx_yyyy_0 = pbuffer.data(idx_eri_0_fg + 10);

    auto g_xxx_yyyz_0 = pbuffer.data(idx_eri_0_fg + 11);

    auto g_xxx_yyzz_0 = pbuffer.data(idx_eri_0_fg + 12);

    auto g_xxx_yzzz_0 = pbuffer.data(idx_eri_0_fg + 13);

    auto g_xxx_zzzz_0 = pbuffer.data(idx_eri_0_fg + 14);

    #pragma omp simd aligned(g_x_xxxx_0, g_x_xxxx_1, g_x_xxxy_0, g_x_xxxy_1, g_x_xxxz_0, g_x_xxxz_1, g_x_xxyy_0, g_x_xxyy_1, g_x_xxyz_0, g_x_xxyz_1, g_x_xxzz_0, g_x_xxzz_1, g_x_xyyy_0, g_x_xyyy_1, g_x_xyyz_0, g_x_xyyz_1, g_x_xyzz_0, g_x_xyzz_1, g_x_xzzz_0, g_x_xzzz_1, g_x_yyyy_0, g_x_yyyy_1, g_x_yyyz_0, g_x_yyyz_1, g_x_yyzz_0, g_x_yyzz_1, g_x_yzzz_0, g_x_yzzz_1, g_x_zzzz_0, g_x_zzzz_1, g_xx_xxx_1, g_xx_xxxx_1, g_xx_xxxy_1, g_xx_xxxz_1, g_xx_xxy_1, g_xx_xxyy_1, g_xx_xxyz_1, g_xx_xxz_1, g_xx_xxzz_1, g_xx_xyy_1, g_xx_xyyy_1, g_xx_xyyz_1, g_xx_xyz_1, g_xx_xyzz_1, g_xx_xzz_1, g_xx_xzzz_1, g_xx_yyy_1, g_xx_yyyy_1, g_xx_yyyz_1, g_xx_yyz_1, g_xx_yyzz_1, g_xx_yzz_1, g_xx_yzzz_1, g_xx_zzz_1, g_xx_zzzz_1, g_xxx_xxxx_0, g_xxx_xxxy_0, g_xxx_xxxz_0, g_xxx_xxyy_0, g_xxx_xxyz_0, g_xxx_xxzz_0, g_xxx_xyyy_0, g_xxx_xyyz_0, g_xxx_xyzz_0, g_xxx_xzzz_0, g_xxx_yyyy_0, g_xxx_yyyz_0, g_xxx_yyzz_0, g_xxx_yzzz_0, g_xxx_zzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxx_xxxx_0[i] = 2.0 * g_x_xxxx_0[i] * fbe_0 - 2.0 * g_x_xxxx_1[i] * fz_be_0 + 4.0 * g_xx_xxx_1[i] * fe_0 + g_xx_xxxx_1[i] * pa_x[i];

        g_xxx_xxxy_0[i] = 2.0 * g_x_xxxy_0[i] * fbe_0 - 2.0 * g_x_xxxy_1[i] * fz_be_0 + 3.0 * g_xx_xxy_1[i] * fe_0 + g_xx_xxxy_1[i] * pa_x[i];

        g_xxx_xxxz_0[i] = 2.0 * g_x_xxxz_0[i] * fbe_0 - 2.0 * g_x_xxxz_1[i] * fz_be_0 + 3.0 * g_xx_xxz_1[i] * fe_0 + g_xx_xxxz_1[i] * pa_x[i];

        g_xxx_xxyy_0[i] = 2.0 * g_x_xxyy_0[i] * fbe_0 - 2.0 * g_x_xxyy_1[i] * fz_be_0 + 2.0 * g_xx_xyy_1[i] * fe_0 + g_xx_xxyy_1[i] * pa_x[i];

        g_xxx_xxyz_0[i] = 2.0 * g_x_xxyz_0[i] * fbe_0 - 2.0 * g_x_xxyz_1[i] * fz_be_0 + 2.0 * g_xx_xyz_1[i] * fe_0 + g_xx_xxyz_1[i] * pa_x[i];

        g_xxx_xxzz_0[i] = 2.0 * g_x_xxzz_0[i] * fbe_0 - 2.0 * g_x_xxzz_1[i] * fz_be_0 + 2.0 * g_xx_xzz_1[i] * fe_0 + g_xx_xxzz_1[i] * pa_x[i];

        g_xxx_xyyy_0[i] = 2.0 * g_x_xyyy_0[i] * fbe_0 - 2.0 * g_x_xyyy_1[i] * fz_be_0 + g_xx_yyy_1[i] * fe_0 + g_xx_xyyy_1[i] * pa_x[i];

        g_xxx_xyyz_0[i] = 2.0 * g_x_xyyz_0[i] * fbe_0 - 2.0 * g_x_xyyz_1[i] * fz_be_0 + g_xx_yyz_1[i] * fe_0 + g_xx_xyyz_1[i] * pa_x[i];

        g_xxx_xyzz_0[i] = 2.0 * g_x_xyzz_0[i] * fbe_0 - 2.0 * g_x_xyzz_1[i] * fz_be_0 + g_xx_yzz_1[i] * fe_0 + g_xx_xyzz_1[i] * pa_x[i];

        g_xxx_xzzz_0[i] = 2.0 * g_x_xzzz_0[i] * fbe_0 - 2.0 * g_x_xzzz_1[i] * fz_be_0 + g_xx_zzz_1[i] * fe_0 + g_xx_xzzz_1[i] * pa_x[i];

        g_xxx_yyyy_0[i] = 2.0 * g_x_yyyy_0[i] * fbe_0 - 2.0 * g_x_yyyy_1[i] * fz_be_0 + g_xx_yyyy_1[i] * pa_x[i];

        g_xxx_yyyz_0[i] = 2.0 * g_x_yyyz_0[i] * fbe_0 - 2.0 * g_x_yyyz_1[i] * fz_be_0 + g_xx_yyyz_1[i] * pa_x[i];

        g_xxx_yyzz_0[i] = 2.0 * g_x_yyzz_0[i] * fbe_0 - 2.0 * g_x_yyzz_1[i] * fz_be_0 + g_xx_yyzz_1[i] * pa_x[i];

        g_xxx_yzzz_0[i] = 2.0 * g_x_yzzz_0[i] * fbe_0 - 2.0 * g_x_yzzz_1[i] * fz_be_0 + g_xx_yzzz_1[i] * pa_x[i];

        g_xxx_zzzz_0[i] = 2.0 * g_x_zzzz_0[i] * fbe_0 - 2.0 * g_x_zzzz_1[i] * fz_be_0 + g_xx_zzzz_1[i] * pa_x[i];
    }

    // Set up 15-30 components of targeted buffer : FG

    auto g_xxy_xxxx_0 = pbuffer.data(idx_eri_0_fg + 15);

    auto g_xxy_xxxy_0 = pbuffer.data(idx_eri_0_fg + 16);

    auto g_xxy_xxxz_0 = pbuffer.data(idx_eri_0_fg + 17);

    auto g_xxy_xxyy_0 = pbuffer.data(idx_eri_0_fg + 18);

    auto g_xxy_xxyz_0 = pbuffer.data(idx_eri_0_fg + 19);

    auto g_xxy_xxzz_0 = pbuffer.data(idx_eri_0_fg + 20);

    auto g_xxy_xyyy_0 = pbuffer.data(idx_eri_0_fg + 21);

    auto g_xxy_xyyz_0 = pbuffer.data(idx_eri_0_fg + 22);

    auto g_xxy_xyzz_0 = pbuffer.data(idx_eri_0_fg + 23);

    auto g_xxy_xzzz_0 = pbuffer.data(idx_eri_0_fg + 24);

    auto g_xxy_yyyy_0 = pbuffer.data(idx_eri_0_fg + 25);

    auto g_xxy_yyyz_0 = pbuffer.data(idx_eri_0_fg + 26);

    auto g_xxy_yyzz_0 = pbuffer.data(idx_eri_0_fg + 27);

    auto g_xxy_yzzz_0 = pbuffer.data(idx_eri_0_fg + 28);

    auto g_xxy_zzzz_0 = pbuffer.data(idx_eri_0_fg + 29);

    #pragma omp simd aligned(g_xx_xxx_1, g_xx_xxxx_1, g_xx_xxxy_1, g_xx_xxxz_1, g_xx_xxy_1, g_xx_xxyy_1, g_xx_xxyz_1, g_xx_xxz_1, g_xx_xxzz_1, g_xx_xyy_1, g_xx_xyyy_1, g_xx_xyyz_1, g_xx_xyz_1, g_xx_xyzz_1, g_xx_xzz_1, g_xx_xzzz_1, g_xx_yyy_1, g_xx_yyyy_1, g_xx_yyyz_1, g_xx_yyz_1, g_xx_yyzz_1, g_xx_yzz_1, g_xx_yzzz_1, g_xx_zzz_1, g_xx_zzzz_1, g_xxy_xxxx_0, g_xxy_xxxy_0, g_xxy_xxxz_0, g_xxy_xxyy_0, g_xxy_xxyz_0, g_xxy_xxzz_0, g_xxy_xyyy_0, g_xxy_xyyz_0, g_xxy_xyzz_0, g_xxy_xzzz_0, g_xxy_yyyy_0, g_xxy_yyyz_0, g_xxy_yyzz_0, g_xxy_yzzz_0, g_xxy_zzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxy_xxxx_0[i] = g_xx_xxxx_1[i] * pa_y[i];

        g_xxy_xxxy_0[i] = g_xx_xxx_1[i] * fe_0 + g_xx_xxxy_1[i] * pa_y[i];

        g_xxy_xxxz_0[i] = g_xx_xxxz_1[i] * pa_y[i];

        g_xxy_xxyy_0[i] = 2.0 * g_xx_xxy_1[i] * fe_0 + g_xx_xxyy_1[i] * pa_y[i];

        g_xxy_xxyz_0[i] = g_xx_xxz_1[i] * fe_0 + g_xx_xxyz_1[i] * pa_y[i];

        g_xxy_xxzz_0[i] = g_xx_xxzz_1[i] * pa_y[i];

        g_xxy_xyyy_0[i] = 3.0 * g_xx_xyy_1[i] * fe_0 + g_xx_xyyy_1[i] * pa_y[i];

        g_xxy_xyyz_0[i] = 2.0 * g_xx_xyz_1[i] * fe_0 + g_xx_xyyz_1[i] * pa_y[i];

        g_xxy_xyzz_0[i] = g_xx_xzz_1[i] * fe_0 + g_xx_xyzz_1[i] * pa_y[i];

        g_xxy_xzzz_0[i] = g_xx_xzzz_1[i] * pa_y[i];

        g_xxy_yyyy_0[i] = 4.0 * g_xx_yyy_1[i] * fe_0 + g_xx_yyyy_1[i] * pa_y[i];

        g_xxy_yyyz_0[i] = 3.0 * g_xx_yyz_1[i] * fe_0 + g_xx_yyyz_1[i] * pa_y[i];

        g_xxy_yyzz_0[i] = 2.0 * g_xx_yzz_1[i] * fe_0 + g_xx_yyzz_1[i] * pa_y[i];

        g_xxy_yzzz_0[i] = g_xx_zzz_1[i] * fe_0 + g_xx_yzzz_1[i] * pa_y[i];

        g_xxy_zzzz_0[i] = g_xx_zzzz_1[i] * pa_y[i];
    }

    // Set up 30-45 components of targeted buffer : FG

    auto g_xxz_xxxx_0 = pbuffer.data(idx_eri_0_fg + 30);

    auto g_xxz_xxxy_0 = pbuffer.data(idx_eri_0_fg + 31);

    auto g_xxz_xxxz_0 = pbuffer.data(idx_eri_0_fg + 32);

    auto g_xxz_xxyy_0 = pbuffer.data(idx_eri_0_fg + 33);

    auto g_xxz_xxyz_0 = pbuffer.data(idx_eri_0_fg + 34);

    auto g_xxz_xxzz_0 = pbuffer.data(idx_eri_0_fg + 35);

    auto g_xxz_xyyy_0 = pbuffer.data(idx_eri_0_fg + 36);

    auto g_xxz_xyyz_0 = pbuffer.data(idx_eri_0_fg + 37);

    auto g_xxz_xyzz_0 = pbuffer.data(idx_eri_0_fg + 38);

    auto g_xxz_xzzz_0 = pbuffer.data(idx_eri_0_fg + 39);

    auto g_xxz_yyyy_0 = pbuffer.data(idx_eri_0_fg + 40);

    auto g_xxz_yyyz_0 = pbuffer.data(idx_eri_0_fg + 41);

    auto g_xxz_yyzz_0 = pbuffer.data(idx_eri_0_fg + 42);

    auto g_xxz_yzzz_0 = pbuffer.data(idx_eri_0_fg + 43);

    auto g_xxz_zzzz_0 = pbuffer.data(idx_eri_0_fg + 44);

    #pragma omp simd aligned(g_xx_xxx_1, g_xx_xxxx_1, g_xx_xxxy_1, g_xx_xxxz_1, g_xx_xxy_1, g_xx_xxyy_1, g_xx_xxyz_1, g_xx_xxz_1, g_xx_xxzz_1, g_xx_xyy_1, g_xx_xyyy_1, g_xx_xyyz_1, g_xx_xyz_1, g_xx_xyzz_1, g_xx_xzz_1, g_xx_xzzz_1, g_xx_yyy_1, g_xx_yyyy_1, g_xx_yyyz_1, g_xx_yyz_1, g_xx_yyzz_1, g_xx_yzz_1, g_xx_yzzz_1, g_xx_zzz_1, g_xx_zzzz_1, g_xxz_xxxx_0, g_xxz_xxxy_0, g_xxz_xxxz_0, g_xxz_xxyy_0, g_xxz_xxyz_0, g_xxz_xxzz_0, g_xxz_xyyy_0, g_xxz_xyyz_0, g_xxz_xyzz_0, g_xxz_xzzz_0, g_xxz_yyyy_0, g_xxz_yyyz_0, g_xxz_yyzz_0, g_xxz_yzzz_0, g_xxz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxz_xxxx_0[i] = g_xx_xxxx_1[i] * pa_z[i];

        g_xxz_xxxy_0[i] = g_xx_xxxy_1[i] * pa_z[i];

        g_xxz_xxxz_0[i] = g_xx_xxx_1[i] * fe_0 + g_xx_xxxz_1[i] * pa_z[i];

        g_xxz_xxyy_0[i] = g_xx_xxyy_1[i] * pa_z[i];

        g_xxz_xxyz_0[i] = g_xx_xxy_1[i] * fe_0 + g_xx_xxyz_1[i] * pa_z[i];

        g_xxz_xxzz_0[i] = 2.0 * g_xx_xxz_1[i] * fe_0 + g_xx_xxzz_1[i] * pa_z[i];

        g_xxz_xyyy_0[i] = g_xx_xyyy_1[i] * pa_z[i];

        g_xxz_xyyz_0[i] = g_xx_xyy_1[i] * fe_0 + g_xx_xyyz_1[i] * pa_z[i];

        g_xxz_xyzz_0[i] = 2.0 * g_xx_xyz_1[i] * fe_0 + g_xx_xyzz_1[i] * pa_z[i];

        g_xxz_xzzz_0[i] = 3.0 * g_xx_xzz_1[i] * fe_0 + g_xx_xzzz_1[i] * pa_z[i];

        g_xxz_yyyy_0[i] = g_xx_yyyy_1[i] * pa_z[i];

        g_xxz_yyyz_0[i] = g_xx_yyy_1[i] * fe_0 + g_xx_yyyz_1[i] * pa_z[i];

        g_xxz_yyzz_0[i] = 2.0 * g_xx_yyz_1[i] * fe_0 + g_xx_yyzz_1[i] * pa_z[i];

        g_xxz_yzzz_0[i] = 3.0 * g_xx_yzz_1[i] * fe_0 + g_xx_yzzz_1[i] * pa_z[i];

        g_xxz_zzzz_0[i] = 4.0 * g_xx_zzz_1[i] * fe_0 + g_xx_zzzz_1[i] * pa_z[i];
    }

    // Set up 45-60 components of targeted buffer : FG

    auto g_xyy_xxxx_0 = pbuffer.data(idx_eri_0_fg + 45);

    auto g_xyy_xxxy_0 = pbuffer.data(idx_eri_0_fg + 46);

    auto g_xyy_xxxz_0 = pbuffer.data(idx_eri_0_fg + 47);

    auto g_xyy_xxyy_0 = pbuffer.data(idx_eri_0_fg + 48);

    auto g_xyy_xxyz_0 = pbuffer.data(idx_eri_0_fg + 49);

    auto g_xyy_xxzz_0 = pbuffer.data(idx_eri_0_fg + 50);

    auto g_xyy_xyyy_0 = pbuffer.data(idx_eri_0_fg + 51);

    auto g_xyy_xyyz_0 = pbuffer.data(idx_eri_0_fg + 52);

    auto g_xyy_xyzz_0 = pbuffer.data(idx_eri_0_fg + 53);

    auto g_xyy_xzzz_0 = pbuffer.data(idx_eri_0_fg + 54);

    auto g_xyy_yyyy_0 = pbuffer.data(idx_eri_0_fg + 55);

    auto g_xyy_yyyz_0 = pbuffer.data(idx_eri_0_fg + 56);

    auto g_xyy_yyzz_0 = pbuffer.data(idx_eri_0_fg + 57);

    auto g_xyy_yzzz_0 = pbuffer.data(idx_eri_0_fg + 58);

    auto g_xyy_zzzz_0 = pbuffer.data(idx_eri_0_fg + 59);

    #pragma omp simd aligned(g_xyy_xxxx_0, g_xyy_xxxy_0, g_xyy_xxxz_0, g_xyy_xxyy_0, g_xyy_xxyz_0, g_xyy_xxzz_0, g_xyy_xyyy_0, g_xyy_xyyz_0, g_xyy_xyzz_0, g_xyy_xzzz_0, g_xyy_yyyy_0, g_xyy_yyyz_0, g_xyy_yyzz_0, g_xyy_yzzz_0, g_xyy_zzzz_0, g_yy_xxx_1, g_yy_xxxx_1, g_yy_xxxy_1, g_yy_xxxz_1, g_yy_xxy_1, g_yy_xxyy_1, g_yy_xxyz_1, g_yy_xxz_1, g_yy_xxzz_1, g_yy_xyy_1, g_yy_xyyy_1, g_yy_xyyz_1, g_yy_xyz_1, g_yy_xyzz_1, g_yy_xzz_1, g_yy_xzzz_1, g_yy_yyy_1, g_yy_yyyy_1, g_yy_yyyz_1, g_yy_yyz_1, g_yy_yyzz_1, g_yy_yzz_1, g_yy_yzzz_1, g_yy_zzz_1, g_yy_zzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyy_xxxx_0[i] = 4.0 * g_yy_xxx_1[i] * fe_0 + g_yy_xxxx_1[i] * pa_x[i];

        g_xyy_xxxy_0[i] = 3.0 * g_yy_xxy_1[i] * fe_0 + g_yy_xxxy_1[i] * pa_x[i];

        g_xyy_xxxz_0[i] = 3.0 * g_yy_xxz_1[i] * fe_0 + g_yy_xxxz_1[i] * pa_x[i];

        g_xyy_xxyy_0[i] = 2.0 * g_yy_xyy_1[i] * fe_0 + g_yy_xxyy_1[i] * pa_x[i];

        g_xyy_xxyz_0[i] = 2.0 * g_yy_xyz_1[i] * fe_0 + g_yy_xxyz_1[i] * pa_x[i];

        g_xyy_xxzz_0[i] = 2.0 * g_yy_xzz_1[i] * fe_0 + g_yy_xxzz_1[i] * pa_x[i];

        g_xyy_xyyy_0[i] = g_yy_yyy_1[i] * fe_0 + g_yy_xyyy_1[i] * pa_x[i];

        g_xyy_xyyz_0[i] = g_yy_yyz_1[i] * fe_0 + g_yy_xyyz_1[i] * pa_x[i];

        g_xyy_xyzz_0[i] = g_yy_yzz_1[i] * fe_0 + g_yy_xyzz_1[i] * pa_x[i];

        g_xyy_xzzz_0[i] = g_yy_zzz_1[i] * fe_0 + g_yy_xzzz_1[i] * pa_x[i];

        g_xyy_yyyy_0[i] = g_yy_yyyy_1[i] * pa_x[i];

        g_xyy_yyyz_0[i] = g_yy_yyyz_1[i] * pa_x[i];

        g_xyy_yyzz_0[i] = g_yy_yyzz_1[i] * pa_x[i];

        g_xyy_yzzz_0[i] = g_yy_yzzz_1[i] * pa_x[i];

        g_xyy_zzzz_0[i] = g_yy_zzzz_1[i] * pa_x[i];
    }

    // Set up 60-75 components of targeted buffer : FG

    auto g_xyz_xxxx_0 = pbuffer.data(idx_eri_0_fg + 60);

    auto g_xyz_xxxy_0 = pbuffer.data(idx_eri_0_fg + 61);

    auto g_xyz_xxxz_0 = pbuffer.data(idx_eri_0_fg + 62);

    auto g_xyz_xxyy_0 = pbuffer.data(idx_eri_0_fg + 63);

    auto g_xyz_xxyz_0 = pbuffer.data(idx_eri_0_fg + 64);

    auto g_xyz_xxzz_0 = pbuffer.data(idx_eri_0_fg + 65);

    auto g_xyz_xyyy_0 = pbuffer.data(idx_eri_0_fg + 66);

    auto g_xyz_xyyz_0 = pbuffer.data(idx_eri_0_fg + 67);

    auto g_xyz_xyzz_0 = pbuffer.data(idx_eri_0_fg + 68);

    auto g_xyz_xzzz_0 = pbuffer.data(idx_eri_0_fg + 69);

    auto g_xyz_yyyy_0 = pbuffer.data(idx_eri_0_fg + 70);

    auto g_xyz_yyyz_0 = pbuffer.data(idx_eri_0_fg + 71);

    auto g_xyz_yyzz_0 = pbuffer.data(idx_eri_0_fg + 72);

    auto g_xyz_yzzz_0 = pbuffer.data(idx_eri_0_fg + 73);

    auto g_xyz_zzzz_0 = pbuffer.data(idx_eri_0_fg + 74);

    #pragma omp simd aligned(g_xy_xxxy_1, g_xy_xxyy_1, g_xy_xyyy_1, g_xyz_xxxx_0, g_xyz_xxxy_0, g_xyz_xxxz_0, g_xyz_xxyy_0, g_xyz_xxyz_0, g_xyz_xxzz_0, g_xyz_xyyy_0, g_xyz_xyyz_0, g_xyz_xyzz_0, g_xyz_xzzz_0, g_xyz_yyyy_0, g_xyz_yyyz_0, g_xyz_yyzz_0, g_xyz_yzzz_0, g_xyz_zzzz_0, g_xz_xxxx_1, g_xz_xxxz_1, g_xz_xxzz_1, g_xz_xzzz_1, g_yz_xxyz_1, g_yz_xyyz_1, g_yz_xyz_1, g_yz_xyzz_1, g_yz_yyyy_1, g_yz_yyyz_1, g_yz_yyz_1, g_yz_yyzz_1, g_yz_yzz_1, g_yz_yzzz_1, g_yz_zzzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyz_xxxx_0[i] = g_xz_xxxx_1[i] * pa_y[i];

        g_xyz_xxxy_0[i] = g_xy_xxxy_1[i] * pa_z[i];

        g_xyz_xxxz_0[i] = g_xz_xxxz_1[i] * pa_y[i];

        g_xyz_xxyy_0[i] = g_xy_xxyy_1[i] * pa_z[i];

        g_xyz_xxyz_0[i] = 2.0 * g_yz_xyz_1[i] * fe_0 + g_yz_xxyz_1[i] * pa_x[i];

        g_xyz_xxzz_0[i] = g_xz_xxzz_1[i] * pa_y[i];

        g_xyz_xyyy_0[i] = g_xy_xyyy_1[i] * pa_z[i];

        g_xyz_xyyz_0[i] = g_yz_yyz_1[i] * fe_0 + g_yz_xyyz_1[i] * pa_x[i];

        g_xyz_xyzz_0[i] = g_yz_yzz_1[i] * fe_0 + g_yz_xyzz_1[i] * pa_x[i];

        g_xyz_xzzz_0[i] = g_xz_xzzz_1[i] * pa_y[i];

        g_xyz_yyyy_0[i] = g_yz_yyyy_1[i] * pa_x[i];

        g_xyz_yyyz_0[i] = g_yz_yyyz_1[i] * pa_x[i];

        g_xyz_yyzz_0[i] = g_yz_yyzz_1[i] * pa_x[i];

        g_xyz_yzzz_0[i] = g_yz_yzzz_1[i] * pa_x[i];

        g_xyz_zzzz_0[i] = g_yz_zzzz_1[i] * pa_x[i];
    }

    // Set up 75-90 components of targeted buffer : FG

    auto g_xzz_xxxx_0 = pbuffer.data(idx_eri_0_fg + 75);

    auto g_xzz_xxxy_0 = pbuffer.data(idx_eri_0_fg + 76);

    auto g_xzz_xxxz_0 = pbuffer.data(idx_eri_0_fg + 77);

    auto g_xzz_xxyy_0 = pbuffer.data(idx_eri_0_fg + 78);

    auto g_xzz_xxyz_0 = pbuffer.data(idx_eri_0_fg + 79);

    auto g_xzz_xxzz_0 = pbuffer.data(idx_eri_0_fg + 80);

    auto g_xzz_xyyy_0 = pbuffer.data(idx_eri_0_fg + 81);

    auto g_xzz_xyyz_0 = pbuffer.data(idx_eri_0_fg + 82);

    auto g_xzz_xyzz_0 = pbuffer.data(idx_eri_0_fg + 83);

    auto g_xzz_xzzz_0 = pbuffer.data(idx_eri_0_fg + 84);

    auto g_xzz_yyyy_0 = pbuffer.data(idx_eri_0_fg + 85);

    auto g_xzz_yyyz_0 = pbuffer.data(idx_eri_0_fg + 86);

    auto g_xzz_yyzz_0 = pbuffer.data(idx_eri_0_fg + 87);

    auto g_xzz_yzzz_0 = pbuffer.data(idx_eri_0_fg + 88);

    auto g_xzz_zzzz_0 = pbuffer.data(idx_eri_0_fg + 89);

    #pragma omp simd aligned(g_xzz_xxxx_0, g_xzz_xxxy_0, g_xzz_xxxz_0, g_xzz_xxyy_0, g_xzz_xxyz_0, g_xzz_xxzz_0, g_xzz_xyyy_0, g_xzz_xyyz_0, g_xzz_xyzz_0, g_xzz_xzzz_0, g_xzz_yyyy_0, g_xzz_yyyz_0, g_xzz_yyzz_0, g_xzz_yzzz_0, g_xzz_zzzz_0, g_zz_xxx_1, g_zz_xxxx_1, g_zz_xxxy_1, g_zz_xxxz_1, g_zz_xxy_1, g_zz_xxyy_1, g_zz_xxyz_1, g_zz_xxz_1, g_zz_xxzz_1, g_zz_xyy_1, g_zz_xyyy_1, g_zz_xyyz_1, g_zz_xyz_1, g_zz_xyzz_1, g_zz_xzz_1, g_zz_xzzz_1, g_zz_yyy_1, g_zz_yyyy_1, g_zz_yyyz_1, g_zz_yyz_1, g_zz_yyzz_1, g_zz_yzz_1, g_zz_yzzz_1, g_zz_zzz_1, g_zz_zzzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzz_xxxx_0[i] = 4.0 * g_zz_xxx_1[i] * fe_0 + g_zz_xxxx_1[i] * pa_x[i];

        g_xzz_xxxy_0[i] = 3.0 * g_zz_xxy_1[i] * fe_0 + g_zz_xxxy_1[i] * pa_x[i];

        g_xzz_xxxz_0[i] = 3.0 * g_zz_xxz_1[i] * fe_0 + g_zz_xxxz_1[i] * pa_x[i];

        g_xzz_xxyy_0[i] = 2.0 * g_zz_xyy_1[i] * fe_0 + g_zz_xxyy_1[i] * pa_x[i];

        g_xzz_xxyz_0[i] = 2.0 * g_zz_xyz_1[i] * fe_0 + g_zz_xxyz_1[i] * pa_x[i];

        g_xzz_xxzz_0[i] = 2.0 * g_zz_xzz_1[i] * fe_0 + g_zz_xxzz_1[i] * pa_x[i];

        g_xzz_xyyy_0[i] = g_zz_yyy_1[i] * fe_0 + g_zz_xyyy_1[i] * pa_x[i];

        g_xzz_xyyz_0[i] = g_zz_yyz_1[i] * fe_0 + g_zz_xyyz_1[i] * pa_x[i];

        g_xzz_xyzz_0[i] = g_zz_yzz_1[i] * fe_0 + g_zz_xyzz_1[i] * pa_x[i];

        g_xzz_xzzz_0[i] = g_zz_zzz_1[i] * fe_0 + g_zz_xzzz_1[i] * pa_x[i];

        g_xzz_yyyy_0[i] = g_zz_yyyy_1[i] * pa_x[i];

        g_xzz_yyyz_0[i] = g_zz_yyyz_1[i] * pa_x[i];

        g_xzz_yyzz_0[i] = g_zz_yyzz_1[i] * pa_x[i];

        g_xzz_yzzz_0[i] = g_zz_yzzz_1[i] * pa_x[i];

        g_xzz_zzzz_0[i] = g_zz_zzzz_1[i] * pa_x[i];
    }

    // Set up 90-105 components of targeted buffer : FG

    auto g_yyy_xxxx_0 = pbuffer.data(idx_eri_0_fg + 90);

    auto g_yyy_xxxy_0 = pbuffer.data(idx_eri_0_fg + 91);

    auto g_yyy_xxxz_0 = pbuffer.data(idx_eri_0_fg + 92);

    auto g_yyy_xxyy_0 = pbuffer.data(idx_eri_0_fg + 93);

    auto g_yyy_xxyz_0 = pbuffer.data(idx_eri_0_fg + 94);

    auto g_yyy_xxzz_0 = pbuffer.data(idx_eri_0_fg + 95);

    auto g_yyy_xyyy_0 = pbuffer.data(idx_eri_0_fg + 96);

    auto g_yyy_xyyz_0 = pbuffer.data(idx_eri_0_fg + 97);

    auto g_yyy_xyzz_0 = pbuffer.data(idx_eri_0_fg + 98);

    auto g_yyy_xzzz_0 = pbuffer.data(idx_eri_0_fg + 99);

    auto g_yyy_yyyy_0 = pbuffer.data(idx_eri_0_fg + 100);

    auto g_yyy_yyyz_0 = pbuffer.data(idx_eri_0_fg + 101);

    auto g_yyy_yyzz_0 = pbuffer.data(idx_eri_0_fg + 102);

    auto g_yyy_yzzz_0 = pbuffer.data(idx_eri_0_fg + 103);

    auto g_yyy_zzzz_0 = pbuffer.data(idx_eri_0_fg + 104);

    #pragma omp simd aligned(g_y_xxxx_0, g_y_xxxx_1, g_y_xxxy_0, g_y_xxxy_1, g_y_xxxz_0, g_y_xxxz_1, g_y_xxyy_0, g_y_xxyy_1, g_y_xxyz_0, g_y_xxyz_1, g_y_xxzz_0, g_y_xxzz_1, g_y_xyyy_0, g_y_xyyy_1, g_y_xyyz_0, g_y_xyyz_1, g_y_xyzz_0, g_y_xyzz_1, g_y_xzzz_0, g_y_xzzz_1, g_y_yyyy_0, g_y_yyyy_1, g_y_yyyz_0, g_y_yyyz_1, g_y_yyzz_0, g_y_yyzz_1, g_y_yzzz_0, g_y_yzzz_1, g_y_zzzz_0, g_y_zzzz_1, g_yy_xxx_1, g_yy_xxxx_1, g_yy_xxxy_1, g_yy_xxxz_1, g_yy_xxy_1, g_yy_xxyy_1, g_yy_xxyz_1, g_yy_xxz_1, g_yy_xxzz_1, g_yy_xyy_1, g_yy_xyyy_1, g_yy_xyyz_1, g_yy_xyz_1, g_yy_xyzz_1, g_yy_xzz_1, g_yy_xzzz_1, g_yy_yyy_1, g_yy_yyyy_1, g_yy_yyyz_1, g_yy_yyz_1, g_yy_yyzz_1, g_yy_yzz_1, g_yy_yzzz_1, g_yy_zzz_1, g_yy_zzzz_1, g_yyy_xxxx_0, g_yyy_xxxy_0, g_yyy_xxxz_0, g_yyy_xxyy_0, g_yyy_xxyz_0, g_yyy_xxzz_0, g_yyy_xyyy_0, g_yyy_xyyz_0, g_yyy_xyzz_0, g_yyy_xzzz_0, g_yyy_yyyy_0, g_yyy_yyyz_0, g_yyy_yyzz_0, g_yyy_yzzz_0, g_yyy_zzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyy_xxxx_0[i] = 2.0 * g_y_xxxx_0[i] * fbe_0 - 2.0 * g_y_xxxx_1[i] * fz_be_0 + g_yy_xxxx_1[i] * pa_y[i];

        g_yyy_xxxy_0[i] = 2.0 * g_y_xxxy_0[i] * fbe_0 - 2.0 * g_y_xxxy_1[i] * fz_be_0 + g_yy_xxx_1[i] * fe_0 + g_yy_xxxy_1[i] * pa_y[i];

        g_yyy_xxxz_0[i] = 2.0 * g_y_xxxz_0[i] * fbe_0 - 2.0 * g_y_xxxz_1[i] * fz_be_0 + g_yy_xxxz_1[i] * pa_y[i];

        g_yyy_xxyy_0[i] = 2.0 * g_y_xxyy_0[i] * fbe_0 - 2.0 * g_y_xxyy_1[i] * fz_be_0 + 2.0 * g_yy_xxy_1[i] * fe_0 + g_yy_xxyy_1[i] * pa_y[i];

        g_yyy_xxyz_0[i] = 2.0 * g_y_xxyz_0[i] * fbe_0 - 2.0 * g_y_xxyz_1[i] * fz_be_0 + g_yy_xxz_1[i] * fe_0 + g_yy_xxyz_1[i] * pa_y[i];

        g_yyy_xxzz_0[i] = 2.0 * g_y_xxzz_0[i] * fbe_0 - 2.0 * g_y_xxzz_1[i] * fz_be_0 + g_yy_xxzz_1[i] * pa_y[i];

        g_yyy_xyyy_0[i] = 2.0 * g_y_xyyy_0[i] * fbe_0 - 2.0 * g_y_xyyy_1[i] * fz_be_0 + 3.0 * g_yy_xyy_1[i] * fe_0 + g_yy_xyyy_1[i] * pa_y[i];

        g_yyy_xyyz_0[i] = 2.0 * g_y_xyyz_0[i] * fbe_0 - 2.0 * g_y_xyyz_1[i] * fz_be_0 + 2.0 * g_yy_xyz_1[i] * fe_0 + g_yy_xyyz_1[i] * pa_y[i];

        g_yyy_xyzz_0[i] = 2.0 * g_y_xyzz_0[i] * fbe_0 - 2.0 * g_y_xyzz_1[i] * fz_be_0 + g_yy_xzz_1[i] * fe_0 + g_yy_xyzz_1[i] * pa_y[i];

        g_yyy_xzzz_0[i] = 2.0 * g_y_xzzz_0[i] * fbe_0 - 2.0 * g_y_xzzz_1[i] * fz_be_0 + g_yy_xzzz_1[i] * pa_y[i];

        g_yyy_yyyy_0[i] = 2.0 * g_y_yyyy_0[i] * fbe_0 - 2.0 * g_y_yyyy_1[i] * fz_be_0 + 4.0 * g_yy_yyy_1[i] * fe_0 + g_yy_yyyy_1[i] * pa_y[i];

        g_yyy_yyyz_0[i] = 2.0 * g_y_yyyz_0[i] * fbe_0 - 2.0 * g_y_yyyz_1[i] * fz_be_0 + 3.0 * g_yy_yyz_1[i] * fe_0 + g_yy_yyyz_1[i] * pa_y[i];

        g_yyy_yyzz_0[i] = 2.0 * g_y_yyzz_0[i] * fbe_0 - 2.0 * g_y_yyzz_1[i] * fz_be_0 + 2.0 * g_yy_yzz_1[i] * fe_0 + g_yy_yyzz_1[i] * pa_y[i];

        g_yyy_yzzz_0[i] = 2.0 * g_y_yzzz_0[i] * fbe_0 - 2.0 * g_y_yzzz_1[i] * fz_be_0 + g_yy_zzz_1[i] * fe_0 + g_yy_yzzz_1[i] * pa_y[i];

        g_yyy_zzzz_0[i] = 2.0 * g_y_zzzz_0[i] * fbe_0 - 2.0 * g_y_zzzz_1[i] * fz_be_0 + g_yy_zzzz_1[i] * pa_y[i];
    }

    // Set up 105-120 components of targeted buffer : FG

    auto g_yyz_xxxx_0 = pbuffer.data(idx_eri_0_fg + 105);

    auto g_yyz_xxxy_0 = pbuffer.data(idx_eri_0_fg + 106);

    auto g_yyz_xxxz_0 = pbuffer.data(idx_eri_0_fg + 107);

    auto g_yyz_xxyy_0 = pbuffer.data(idx_eri_0_fg + 108);

    auto g_yyz_xxyz_0 = pbuffer.data(idx_eri_0_fg + 109);

    auto g_yyz_xxzz_0 = pbuffer.data(idx_eri_0_fg + 110);

    auto g_yyz_xyyy_0 = pbuffer.data(idx_eri_0_fg + 111);

    auto g_yyz_xyyz_0 = pbuffer.data(idx_eri_0_fg + 112);

    auto g_yyz_xyzz_0 = pbuffer.data(idx_eri_0_fg + 113);

    auto g_yyz_xzzz_0 = pbuffer.data(idx_eri_0_fg + 114);

    auto g_yyz_yyyy_0 = pbuffer.data(idx_eri_0_fg + 115);

    auto g_yyz_yyyz_0 = pbuffer.data(idx_eri_0_fg + 116);

    auto g_yyz_yyzz_0 = pbuffer.data(idx_eri_0_fg + 117);

    auto g_yyz_yzzz_0 = pbuffer.data(idx_eri_0_fg + 118);

    auto g_yyz_zzzz_0 = pbuffer.data(idx_eri_0_fg + 119);

    #pragma omp simd aligned(g_yy_xxx_1, g_yy_xxxx_1, g_yy_xxxy_1, g_yy_xxxz_1, g_yy_xxy_1, g_yy_xxyy_1, g_yy_xxyz_1, g_yy_xxz_1, g_yy_xxzz_1, g_yy_xyy_1, g_yy_xyyy_1, g_yy_xyyz_1, g_yy_xyz_1, g_yy_xyzz_1, g_yy_xzz_1, g_yy_xzzz_1, g_yy_yyy_1, g_yy_yyyy_1, g_yy_yyyz_1, g_yy_yyz_1, g_yy_yyzz_1, g_yy_yzz_1, g_yy_yzzz_1, g_yy_zzz_1, g_yy_zzzz_1, g_yyz_xxxx_0, g_yyz_xxxy_0, g_yyz_xxxz_0, g_yyz_xxyy_0, g_yyz_xxyz_0, g_yyz_xxzz_0, g_yyz_xyyy_0, g_yyz_xyyz_0, g_yyz_xyzz_0, g_yyz_xzzz_0, g_yyz_yyyy_0, g_yyz_yyyz_0, g_yyz_yyzz_0, g_yyz_yzzz_0, g_yyz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyz_xxxx_0[i] = g_yy_xxxx_1[i] * pa_z[i];

        g_yyz_xxxy_0[i] = g_yy_xxxy_1[i] * pa_z[i];

        g_yyz_xxxz_0[i] = g_yy_xxx_1[i] * fe_0 + g_yy_xxxz_1[i] * pa_z[i];

        g_yyz_xxyy_0[i] = g_yy_xxyy_1[i] * pa_z[i];

        g_yyz_xxyz_0[i] = g_yy_xxy_1[i] * fe_0 + g_yy_xxyz_1[i] * pa_z[i];

        g_yyz_xxzz_0[i] = 2.0 * g_yy_xxz_1[i] * fe_0 + g_yy_xxzz_1[i] * pa_z[i];

        g_yyz_xyyy_0[i] = g_yy_xyyy_1[i] * pa_z[i];

        g_yyz_xyyz_0[i] = g_yy_xyy_1[i] * fe_0 + g_yy_xyyz_1[i] * pa_z[i];

        g_yyz_xyzz_0[i] = 2.0 * g_yy_xyz_1[i] * fe_0 + g_yy_xyzz_1[i] * pa_z[i];

        g_yyz_xzzz_0[i] = 3.0 * g_yy_xzz_1[i] * fe_0 + g_yy_xzzz_1[i] * pa_z[i];

        g_yyz_yyyy_0[i] = g_yy_yyyy_1[i] * pa_z[i];

        g_yyz_yyyz_0[i] = g_yy_yyy_1[i] * fe_0 + g_yy_yyyz_1[i] * pa_z[i];

        g_yyz_yyzz_0[i] = 2.0 * g_yy_yyz_1[i] * fe_0 + g_yy_yyzz_1[i] * pa_z[i];

        g_yyz_yzzz_0[i] = 3.0 * g_yy_yzz_1[i] * fe_0 + g_yy_yzzz_1[i] * pa_z[i];

        g_yyz_zzzz_0[i] = 4.0 * g_yy_zzz_1[i] * fe_0 + g_yy_zzzz_1[i] * pa_z[i];
    }

    // Set up 120-135 components of targeted buffer : FG

    auto g_yzz_xxxx_0 = pbuffer.data(idx_eri_0_fg + 120);

    auto g_yzz_xxxy_0 = pbuffer.data(idx_eri_0_fg + 121);

    auto g_yzz_xxxz_0 = pbuffer.data(idx_eri_0_fg + 122);

    auto g_yzz_xxyy_0 = pbuffer.data(idx_eri_0_fg + 123);

    auto g_yzz_xxyz_0 = pbuffer.data(idx_eri_0_fg + 124);

    auto g_yzz_xxzz_0 = pbuffer.data(idx_eri_0_fg + 125);

    auto g_yzz_xyyy_0 = pbuffer.data(idx_eri_0_fg + 126);

    auto g_yzz_xyyz_0 = pbuffer.data(idx_eri_0_fg + 127);

    auto g_yzz_xyzz_0 = pbuffer.data(idx_eri_0_fg + 128);

    auto g_yzz_xzzz_0 = pbuffer.data(idx_eri_0_fg + 129);

    auto g_yzz_yyyy_0 = pbuffer.data(idx_eri_0_fg + 130);

    auto g_yzz_yyyz_0 = pbuffer.data(idx_eri_0_fg + 131);

    auto g_yzz_yyzz_0 = pbuffer.data(idx_eri_0_fg + 132);

    auto g_yzz_yzzz_0 = pbuffer.data(idx_eri_0_fg + 133);

    auto g_yzz_zzzz_0 = pbuffer.data(idx_eri_0_fg + 134);

    #pragma omp simd aligned(g_yzz_xxxx_0, g_yzz_xxxy_0, g_yzz_xxxz_0, g_yzz_xxyy_0, g_yzz_xxyz_0, g_yzz_xxzz_0, g_yzz_xyyy_0, g_yzz_xyyz_0, g_yzz_xyzz_0, g_yzz_xzzz_0, g_yzz_yyyy_0, g_yzz_yyyz_0, g_yzz_yyzz_0, g_yzz_yzzz_0, g_yzz_zzzz_0, g_zz_xxx_1, g_zz_xxxx_1, g_zz_xxxy_1, g_zz_xxxz_1, g_zz_xxy_1, g_zz_xxyy_1, g_zz_xxyz_1, g_zz_xxz_1, g_zz_xxzz_1, g_zz_xyy_1, g_zz_xyyy_1, g_zz_xyyz_1, g_zz_xyz_1, g_zz_xyzz_1, g_zz_xzz_1, g_zz_xzzz_1, g_zz_yyy_1, g_zz_yyyy_1, g_zz_yyyz_1, g_zz_yyz_1, g_zz_yyzz_1, g_zz_yzz_1, g_zz_yzzz_1, g_zz_zzz_1, g_zz_zzzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzz_xxxx_0[i] = g_zz_xxxx_1[i] * pa_y[i];

        g_yzz_xxxy_0[i] = g_zz_xxx_1[i] * fe_0 + g_zz_xxxy_1[i] * pa_y[i];

        g_yzz_xxxz_0[i] = g_zz_xxxz_1[i] * pa_y[i];

        g_yzz_xxyy_0[i] = 2.0 * g_zz_xxy_1[i] * fe_0 + g_zz_xxyy_1[i] * pa_y[i];

        g_yzz_xxyz_0[i] = g_zz_xxz_1[i] * fe_0 + g_zz_xxyz_1[i] * pa_y[i];

        g_yzz_xxzz_0[i] = g_zz_xxzz_1[i] * pa_y[i];

        g_yzz_xyyy_0[i] = 3.0 * g_zz_xyy_1[i] * fe_0 + g_zz_xyyy_1[i] * pa_y[i];

        g_yzz_xyyz_0[i] = 2.0 * g_zz_xyz_1[i] * fe_0 + g_zz_xyyz_1[i] * pa_y[i];

        g_yzz_xyzz_0[i] = g_zz_xzz_1[i] * fe_0 + g_zz_xyzz_1[i] * pa_y[i];

        g_yzz_xzzz_0[i] = g_zz_xzzz_1[i] * pa_y[i];

        g_yzz_yyyy_0[i] = 4.0 * g_zz_yyy_1[i] * fe_0 + g_zz_yyyy_1[i] * pa_y[i];

        g_yzz_yyyz_0[i] = 3.0 * g_zz_yyz_1[i] * fe_0 + g_zz_yyyz_1[i] * pa_y[i];

        g_yzz_yyzz_0[i] = 2.0 * g_zz_yzz_1[i] * fe_0 + g_zz_yyzz_1[i] * pa_y[i];

        g_yzz_yzzz_0[i] = g_zz_zzz_1[i] * fe_0 + g_zz_yzzz_1[i] * pa_y[i];

        g_yzz_zzzz_0[i] = g_zz_zzzz_1[i] * pa_y[i];
    }

    // Set up 135-150 components of targeted buffer : FG

    auto g_zzz_xxxx_0 = pbuffer.data(idx_eri_0_fg + 135);

    auto g_zzz_xxxy_0 = pbuffer.data(idx_eri_0_fg + 136);

    auto g_zzz_xxxz_0 = pbuffer.data(idx_eri_0_fg + 137);

    auto g_zzz_xxyy_0 = pbuffer.data(idx_eri_0_fg + 138);

    auto g_zzz_xxyz_0 = pbuffer.data(idx_eri_0_fg + 139);

    auto g_zzz_xxzz_0 = pbuffer.data(idx_eri_0_fg + 140);

    auto g_zzz_xyyy_0 = pbuffer.data(idx_eri_0_fg + 141);

    auto g_zzz_xyyz_0 = pbuffer.data(idx_eri_0_fg + 142);

    auto g_zzz_xyzz_0 = pbuffer.data(idx_eri_0_fg + 143);

    auto g_zzz_xzzz_0 = pbuffer.data(idx_eri_0_fg + 144);

    auto g_zzz_yyyy_0 = pbuffer.data(idx_eri_0_fg + 145);

    auto g_zzz_yyyz_0 = pbuffer.data(idx_eri_0_fg + 146);

    auto g_zzz_yyzz_0 = pbuffer.data(idx_eri_0_fg + 147);

    auto g_zzz_yzzz_0 = pbuffer.data(idx_eri_0_fg + 148);

    auto g_zzz_zzzz_0 = pbuffer.data(idx_eri_0_fg + 149);

    #pragma omp simd aligned(g_z_xxxx_0, g_z_xxxx_1, g_z_xxxy_0, g_z_xxxy_1, g_z_xxxz_0, g_z_xxxz_1, g_z_xxyy_0, g_z_xxyy_1, g_z_xxyz_0, g_z_xxyz_1, g_z_xxzz_0, g_z_xxzz_1, g_z_xyyy_0, g_z_xyyy_1, g_z_xyyz_0, g_z_xyyz_1, g_z_xyzz_0, g_z_xyzz_1, g_z_xzzz_0, g_z_xzzz_1, g_z_yyyy_0, g_z_yyyy_1, g_z_yyyz_0, g_z_yyyz_1, g_z_yyzz_0, g_z_yyzz_1, g_z_yzzz_0, g_z_yzzz_1, g_z_zzzz_0, g_z_zzzz_1, g_zz_xxx_1, g_zz_xxxx_1, g_zz_xxxy_1, g_zz_xxxz_1, g_zz_xxy_1, g_zz_xxyy_1, g_zz_xxyz_1, g_zz_xxz_1, g_zz_xxzz_1, g_zz_xyy_1, g_zz_xyyy_1, g_zz_xyyz_1, g_zz_xyz_1, g_zz_xyzz_1, g_zz_xzz_1, g_zz_xzzz_1, g_zz_yyy_1, g_zz_yyyy_1, g_zz_yyyz_1, g_zz_yyz_1, g_zz_yyzz_1, g_zz_yzz_1, g_zz_yzzz_1, g_zz_zzz_1, g_zz_zzzz_1, g_zzz_xxxx_0, g_zzz_xxxy_0, g_zzz_xxxz_0, g_zzz_xxyy_0, g_zzz_xxyz_0, g_zzz_xxzz_0, g_zzz_xyyy_0, g_zzz_xyyz_0, g_zzz_xyzz_0, g_zzz_xzzz_0, g_zzz_yyyy_0, g_zzz_yyyz_0, g_zzz_yyzz_0, g_zzz_yzzz_0, g_zzz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzz_xxxx_0[i] = 2.0 * g_z_xxxx_0[i] * fbe_0 - 2.0 * g_z_xxxx_1[i] * fz_be_0 + g_zz_xxxx_1[i] * pa_z[i];

        g_zzz_xxxy_0[i] = 2.0 * g_z_xxxy_0[i] * fbe_0 - 2.0 * g_z_xxxy_1[i] * fz_be_0 + g_zz_xxxy_1[i] * pa_z[i];

        g_zzz_xxxz_0[i] = 2.0 * g_z_xxxz_0[i] * fbe_0 - 2.0 * g_z_xxxz_1[i] * fz_be_0 + g_zz_xxx_1[i] * fe_0 + g_zz_xxxz_1[i] * pa_z[i];

        g_zzz_xxyy_0[i] = 2.0 * g_z_xxyy_0[i] * fbe_0 - 2.0 * g_z_xxyy_1[i] * fz_be_0 + g_zz_xxyy_1[i] * pa_z[i];

        g_zzz_xxyz_0[i] = 2.0 * g_z_xxyz_0[i] * fbe_0 - 2.0 * g_z_xxyz_1[i] * fz_be_0 + g_zz_xxy_1[i] * fe_0 + g_zz_xxyz_1[i] * pa_z[i];

        g_zzz_xxzz_0[i] = 2.0 * g_z_xxzz_0[i] * fbe_0 - 2.0 * g_z_xxzz_1[i] * fz_be_0 + 2.0 * g_zz_xxz_1[i] * fe_0 + g_zz_xxzz_1[i] * pa_z[i];

        g_zzz_xyyy_0[i] = 2.0 * g_z_xyyy_0[i] * fbe_0 - 2.0 * g_z_xyyy_1[i] * fz_be_0 + g_zz_xyyy_1[i] * pa_z[i];

        g_zzz_xyyz_0[i] = 2.0 * g_z_xyyz_0[i] * fbe_0 - 2.0 * g_z_xyyz_1[i] * fz_be_0 + g_zz_xyy_1[i] * fe_0 + g_zz_xyyz_1[i] * pa_z[i];

        g_zzz_xyzz_0[i] = 2.0 * g_z_xyzz_0[i] * fbe_0 - 2.0 * g_z_xyzz_1[i] * fz_be_0 + 2.0 * g_zz_xyz_1[i] * fe_0 + g_zz_xyzz_1[i] * pa_z[i];

        g_zzz_xzzz_0[i] = 2.0 * g_z_xzzz_0[i] * fbe_0 - 2.0 * g_z_xzzz_1[i] * fz_be_0 + 3.0 * g_zz_xzz_1[i] * fe_0 + g_zz_xzzz_1[i] * pa_z[i];

        g_zzz_yyyy_0[i] = 2.0 * g_z_yyyy_0[i] * fbe_0 - 2.0 * g_z_yyyy_1[i] * fz_be_0 + g_zz_yyyy_1[i] * pa_z[i];

        g_zzz_yyyz_0[i] = 2.0 * g_z_yyyz_0[i] * fbe_0 - 2.0 * g_z_yyyz_1[i] * fz_be_0 + g_zz_yyy_1[i] * fe_0 + g_zz_yyyz_1[i] * pa_z[i];

        g_zzz_yyzz_0[i] = 2.0 * g_z_yyzz_0[i] * fbe_0 - 2.0 * g_z_yyzz_1[i] * fz_be_0 + 2.0 * g_zz_yyz_1[i] * fe_0 + g_zz_yyzz_1[i] * pa_z[i];

        g_zzz_yzzz_0[i] = 2.0 * g_z_yzzz_0[i] * fbe_0 - 2.0 * g_z_yzzz_1[i] * fz_be_0 + 3.0 * g_zz_yzz_1[i] * fe_0 + g_zz_yzzz_1[i] * pa_z[i];

        g_zzz_zzzz_0[i] = 2.0 * g_z_zzzz_0[i] * fbe_0 - 2.0 * g_z_zzzz_1[i] * fz_be_0 + 4.0 * g_zz_zzz_1[i] * fe_0 + g_zz_zzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

