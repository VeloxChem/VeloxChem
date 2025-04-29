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

#include "TwoCenterElectronRepulsionPrimRecDG.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_dg(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_dg,
                                const size_t idx_eri_0_sg,
                                const size_t idx_eri_1_sg,
                                const size_t idx_eri_1_pf,
                                const size_t idx_eri_1_pg,
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

    // Set up components of auxiliary buffer : SG

    auto g_0_xxxx_0 = pbuffer.data(idx_eri_0_sg);

    auto g_0_xxxy_0 = pbuffer.data(idx_eri_0_sg + 1);

    auto g_0_xxxz_0 = pbuffer.data(idx_eri_0_sg + 2);

    auto g_0_xxyy_0 = pbuffer.data(idx_eri_0_sg + 3);

    auto g_0_xxyz_0 = pbuffer.data(idx_eri_0_sg + 4);

    auto g_0_xxzz_0 = pbuffer.data(idx_eri_0_sg + 5);

    auto g_0_xyyy_0 = pbuffer.data(idx_eri_0_sg + 6);

    auto g_0_xyyz_0 = pbuffer.data(idx_eri_0_sg + 7);

    auto g_0_xyzz_0 = pbuffer.data(idx_eri_0_sg + 8);

    auto g_0_xzzz_0 = pbuffer.data(idx_eri_0_sg + 9);

    auto g_0_yyyy_0 = pbuffer.data(idx_eri_0_sg + 10);

    auto g_0_yyyz_0 = pbuffer.data(idx_eri_0_sg + 11);

    auto g_0_yyzz_0 = pbuffer.data(idx_eri_0_sg + 12);

    auto g_0_yzzz_0 = pbuffer.data(idx_eri_0_sg + 13);

    auto g_0_zzzz_0 = pbuffer.data(idx_eri_0_sg + 14);

    // Set up components of auxiliary buffer : SG

    auto g_0_xxxx_1 = pbuffer.data(idx_eri_1_sg);

    auto g_0_xxxy_1 = pbuffer.data(idx_eri_1_sg + 1);

    auto g_0_xxxz_1 = pbuffer.data(idx_eri_1_sg + 2);

    auto g_0_xxyy_1 = pbuffer.data(idx_eri_1_sg + 3);

    auto g_0_xxyz_1 = pbuffer.data(idx_eri_1_sg + 4);

    auto g_0_xxzz_1 = pbuffer.data(idx_eri_1_sg + 5);

    auto g_0_xyyy_1 = pbuffer.data(idx_eri_1_sg + 6);

    auto g_0_xyyz_1 = pbuffer.data(idx_eri_1_sg + 7);

    auto g_0_xyzz_1 = pbuffer.data(idx_eri_1_sg + 8);

    auto g_0_xzzz_1 = pbuffer.data(idx_eri_1_sg + 9);

    auto g_0_yyyy_1 = pbuffer.data(idx_eri_1_sg + 10);

    auto g_0_yyyz_1 = pbuffer.data(idx_eri_1_sg + 11);

    auto g_0_yyzz_1 = pbuffer.data(idx_eri_1_sg + 12);

    auto g_0_yzzz_1 = pbuffer.data(idx_eri_1_sg + 13);

    auto g_0_zzzz_1 = pbuffer.data(idx_eri_1_sg + 14);

    // Set up components of auxiliary buffer : PF

    auto g_x_xxx_1 = pbuffer.data(idx_eri_1_pf);

    auto g_x_xxy_1 = pbuffer.data(idx_eri_1_pf + 1);

    auto g_x_xxz_1 = pbuffer.data(idx_eri_1_pf + 2);

    auto g_x_xyy_1 = pbuffer.data(idx_eri_1_pf + 3);

    auto g_x_xyz_1 = pbuffer.data(idx_eri_1_pf + 4);

    auto g_x_xzz_1 = pbuffer.data(idx_eri_1_pf + 5);

    auto g_x_yyy_1 = pbuffer.data(idx_eri_1_pf + 6);

    auto g_x_yyz_1 = pbuffer.data(idx_eri_1_pf + 7);

    auto g_x_yzz_1 = pbuffer.data(idx_eri_1_pf + 8);

    auto g_x_zzz_1 = pbuffer.data(idx_eri_1_pf + 9);

    auto g_y_xxx_1 = pbuffer.data(idx_eri_1_pf + 10);

    auto g_y_xxy_1 = pbuffer.data(idx_eri_1_pf + 11);

    auto g_y_xxz_1 = pbuffer.data(idx_eri_1_pf + 12);

    auto g_y_xyy_1 = pbuffer.data(idx_eri_1_pf + 13);

    auto g_y_xyz_1 = pbuffer.data(idx_eri_1_pf + 14);

    auto g_y_xzz_1 = pbuffer.data(idx_eri_1_pf + 15);

    auto g_y_yyy_1 = pbuffer.data(idx_eri_1_pf + 16);

    auto g_y_yyz_1 = pbuffer.data(idx_eri_1_pf + 17);

    auto g_y_yzz_1 = pbuffer.data(idx_eri_1_pf + 18);

    auto g_y_zzz_1 = pbuffer.data(idx_eri_1_pf + 19);

    auto g_z_xxx_1 = pbuffer.data(idx_eri_1_pf + 20);

    auto g_z_xxy_1 = pbuffer.data(idx_eri_1_pf + 21);

    auto g_z_xxz_1 = pbuffer.data(idx_eri_1_pf + 22);

    auto g_z_xyy_1 = pbuffer.data(idx_eri_1_pf + 23);

    auto g_z_xyz_1 = pbuffer.data(idx_eri_1_pf + 24);

    auto g_z_xzz_1 = pbuffer.data(idx_eri_1_pf + 25);

    auto g_z_yyy_1 = pbuffer.data(idx_eri_1_pf + 26);

    auto g_z_yyz_1 = pbuffer.data(idx_eri_1_pf + 27);

    auto g_z_yzz_1 = pbuffer.data(idx_eri_1_pf + 28);

    auto g_z_zzz_1 = pbuffer.data(idx_eri_1_pf + 29);

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

    // Set up 0-15 components of targeted buffer : DG

    auto g_xx_xxxx_0 = pbuffer.data(idx_eri_0_dg);

    auto g_xx_xxxy_0 = pbuffer.data(idx_eri_0_dg + 1);

    auto g_xx_xxxz_0 = pbuffer.data(idx_eri_0_dg + 2);

    auto g_xx_xxyy_0 = pbuffer.data(idx_eri_0_dg + 3);

    auto g_xx_xxyz_0 = pbuffer.data(idx_eri_0_dg + 4);

    auto g_xx_xxzz_0 = pbuffer.data(idx_eri_0_dg + 5);

    auto g_xx_xyyy_0 = pbuffer.data(idx_eri_0_dg + 6);

    auto g_xx_xyyz_0 = pbuffer.data(idx_eri_0_dg + 7);

    auto g_xx_xyzz_0 = pbuffer.data(idx_eri_0_dg + 8);

    auto g_xx_xzzz_0 = pbuffer.data(idx_eri_0_dg + 9);

    auto g_xx_yyyy_0 = pbuffer.data(idx_eri_0_dg + 10);

    auto g_xx_yyyz_0 = pbuffer.data(idx_eri_0_dg + 11);

    auto g_xx_yyzz_0 = pbuffer.data(idx_eri_0_dg + 12);

    auto g_xx_yzzz_0 = pbuffer.data(idx_eri_0_dg + 13);

    auto g_xx_zzzz_0 = pbuffer.data(idx_eri_0_dg + 14);

    #pragma omp simd aligned(g_0_xxxx_0, g_0_xxxx_1, g_0_xxxy_0, g_0_xxxy_1, g_0_xxxz_0, g_0_xxxz_1, g_0_xxyy_0, g_0_xxyy_1, g_0_xxyz_0, g_0_xxyz_1, g_0_xxzz_0, g_0_xxzz_1, g_0_xyyy_0, g_0_xyyy_1, g_0_xyyz_0, g_0_xyyz_1, g_0_xyzz_0, g_0_xyzz_1, g_0_xzzz_0, g_0_xzzz_1, g_0_yyyy_0, g_0_yyyy_1, g_0_yyyz_0, g_0_yyyz_1, g_0_yyzz_0, g_0_yyzz_1, g_0_yzzz_0, g_0_yzzz_1, g_0_zzzz_0, g_0_zzzz_1, g_x_xxx_1, g_x_xxxx_1, g_x_xxxy_1, g_x_xxxz_1, g_x_xxy_1, g_x_xxyy_1, g_x_xxyz_1, g_x_xxz_1, g_x_xxzz_1, g_x_xyy_1, g_x_xyyy_1, g_x_xyyz_1, g_x_xyz_1, g_x_xyzz_1, g_x_xzz_1, g_x_xzzz_1, g_x_yyy_1, g_x_yyyy_1, g_x_yyyz_1, g_x_yyz_1, g_x_yyzz_1, g_x_yzz_1, g_x_yzzz_1, g_x_zzz_1, g_x_zzzz_1, g_xx_xxxx_0, g_xx_xxxy_0, g_xx_xxxz_0, g_xx_xxyy_0, g_xx_xxyz_0, g_xx_xxzz_0, g_xx_xyyy_0, g_xx_xyyz_0, g_xx_xyzz_0, g_xx_xzzz_0, g_xx_yyyy_0, g_xx_yyyz_0, g_xx_yyzz_0, g_xx_yzzz_0, g_xx_zzzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xx_xxxx_0[i] = g_0_xxxx_0[i] * fbe_0 - g_0_xxxx_1[i] * fz_be_0 + 4.0 * g_x_xxx_1[i] * fe_0 + g_x_xxxx_1[i] * pa_x[i];

        g_xx_xxxy_0[i] = g_0_xxxy_0[i] * fbe_0 - g_0_xxxy_1[i] * fz_be_0 + 3.0 * g_x_xxy_1[i] * fe_0 + g_x_xxxy_1[i] * pa_x[i];

        g_xx_xxxz_0[i] = g_0_xxxz_0[i] * fbe_0 - g_0_xxxz_1[i] * fz_be_0 + 3.0 * g_x_xxz_1[i] * fe_0 + g_x_xxxz_1[i] * pa_x[i];

        g_xx_xxyy_0[i] = g_0_xxyy_0[i] * fbe_0 - g_0_xxyy_1[i] * fz_be_0 + 2.0 * g_x_xyy_1[i] * fe_0 + g_x_xxyy_1[i] * pa_x[i];

        g_xx_xxyz_0[i] = g_0_xxyz_0[i] * fbe_0 - g_0_xxyz_1[i] * fz_be_0 + 2.0 * g_x_xyz_1[i] * fe_0 + g_x_xxyz_1[i] * pa_x[i];

        g_xx_xxzz_0[i] = g_0_xxzz_0[i] * fbe_0 - g_0_xxzz_1[i] * fz_be_0 + 2.0 * g_x_xzz_1[i] * fe_0 + g_x_xxzz_1[i] * pa_x[i];

        g_xx_xyyy_0[i] = g_0_xyyy_0[i] * fbe_0 - g_0_xyyy_1[i] * fz_be_0 + g_x_yyy_1[i] * fe_0 + g_x_xyyy_1[i] * pa_x[i];

        g_xx_xyyz_0[i] = g_0_xyyz_0[i] * fbe_0 - g_0_xyyz_1[i] * fz_be_0 + g_x_yyz_1[i] * fe_0 + g_x_xyyz_1[i] * pa_x[i];

        g_xx_xyzz_0[i] = g_0_xyzz_0[i] * fbe_0 - g_0_xyzz_1[i] * fz_be_0 + g_x_yzz_1[i] * fe_0 + g_x_xyzz_1[i] * pa_x[i];

        g_xx_xzzz_0[i] = g_0_xzzz_0[i] * fbe_0 - g_0_xzzz_1[i] * fz_be_0 + g_x_zzz_1[i] * fe_0 + g_x_xzzz_1[i] * pa_x[i];

        g_xx_yyyy_0[i] = g_0_yyyy_0[i] * fbe_0 - g_0_yyyy_1[i] * fz_be_0 + g_x_yyyy_1[i] * pa_x[i];

        g_xx_yyyz_0[i] = g_0_yyyz_0[i] * fbe_0 - g_0_yyyz_1[i] * fz_be_0 + g_x_yyyz_1[i] * pa_x[i];

        g_xx_yyzz_0[i] = g_0_yyzz_0[i] * fbe_0 - g_0_yyzz_1[i] * fz_be_0 + g_x_yyzz_1[i] * pa_x[i];

        g_xx_yzzz_0[i] = g_0_yzzz_0[i] * fbe_0 - g_0_yzzz_1[i] * fz_be_0 + g_x_yzzz_1[i] * pa_x[i];

        g_xx_zzzz_0[i] = g_0_zzzz_0[i] * fbe_0 - g_0_zzzz_1[i] * fz_be_0 + g_x_zzzz_1[i] * pa_x[i];
    }

    // Set up 15-30 components of targeted buffer : DG

    auto g_xy_xxxx_0 = pbuffer.data(idx_eri_0_dg + 15);

    auto g_xy_xxxy_0 = pbuffer.data(idx_eri_0_dg + 16);

    auto g_xy_xxxz_0 = pbuffer.data(idx_eri_0_dg + 17);

    auto g_xy_xxyy_0 = pbuffer.data(idx_eri_0_dg + 18);

    auto g_xy_xxyz_0 = pbuffer.data(idx_eri_0_dg + 19);

    auto g_xy_xxzz_0 = pbuffer.data(idx_eri_0_dg + 20);

    auto g_xy_xyyy_0 = pbuffer.data(idx_eri_0_dg + 21);

    auto g_xy_xyyz_0 = pbuffer.data(idx_eri_0_dg + 22);

    auto g_xy_xyzz_0 = pbuffer.data(idx_eri_0_dg + 23);

    auto g_xy_xzzz_0 = pbuffer.data(idx_eri_0_dg + 24);

    auto g_xy_yyyy_0 = pbuffer.data(idx_eri_0_dg + 25);

    auto g_xy_yyyz_0 = pbuffer.data(idx_eri_0_dg + 26);

    auto g_xy_yyzz_0 = pbuffer.data(idx_eri_0_dg + 27);

    auto g_xy_yzzz_0 = pbuffer.data(idx_eri_0_dg + 28);

    auto g_xy_zzzz_0 = pbuffer.data(idx_eri_0_dg + 29);

    #pragma omp simd aligned(g_x_xxxx_1, g_x_xxxz_1, g_x_xxzz_1, g_x_xzzz_1, g_xy_xxxx_0, g_xy_xxxy_0, g_xy_xxxz_0, g_xy_xxyy_0, g_xy_xxyz_0, g_xy_xxzz_0, g_xy_xyyy_0, g_xy_xyyz_0, g_xy_xyzz_0, g_xy_xzzz_0, g_xy_yyyy_0, g_xy_yyyz_0, g_xy_yyzz_0, g_xy_yzzz_0, g_xy_zzzz_0, g_y_xxxy_1, g_y_xxy_1, g_y_xxyy_1, g_y_xxyz_1, g_y_xyy_1, g_y_xyyy_1, g_y_xyyz_1, g_y_xyz_1, g_y_xyzz_1, g_y_yyy_1, g_y_yyyy_1, g_y_yyyz_1, g_y_yyz_1, g_y_yyzz_1, g_y_yzz_1, g_y_yzzz_1, g_y_zzzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xy_xxxx_0[i] = g_x_xxxx_1[i] * pa_y[i];

        g_xy_xxxy_0[i] = 3.0 * g_y_xxy_1[i] * fe_0 + g_y_xxxy_1[i] * pa_x[i];

        g_xy_xxxz_0[i] = g_x_xxxz_1[i] * pa_y[i];

        g_xy_xxyy_0[i] = 2.0 * g_y_xyy_1[i] * fe_0 + g_y_xxyy_1[i] * pa_x[i];

        g_xy_xxyz_0[i] = 2.0 * g_y_xyz_1[i] * fe_0 + g_y_xxyz_1[i] * pa_x[i];

        g_xy_xxzz_0[i] = g_x_xxzz_1[i] * pa_y[i];

        g_xy_xyyy_0[i] = g_y_yyy_1[i] * fe_0 + g_y_xyyy_1[i] * pa_x[i];

        g_xy_xyyz_0[i] = g_y_yyz_1[i] * fe_0 + g_y_xyyz_1[i] * pa_x[i];

        g_xy_xyzz_0[i] = g_y_yzz_1[i] * fe_0 + g_y_xyzz_1[i] * pa_x[i];

        g_xy_xzzz_0[i] = g_x_xzzz_1[i] * pa_y[i];

        g_xy_yyyy_0[i] = g_y_yyyy_1[i] * pa_x[i];

        g_xy_yyyz_0[i] = g_y_yyyz_1[i] * pa_x[i];

        g_xy_yyzz_0[i] = g_y_yyzz_1[i] * pa_x[i];

        g_xy_yzzz_0[i] = g_y_yzzz_1[i] * pa_x[i];

        g_xy_zzzz_0[i] = g_y_zzzz_1[i] * pa_x[i];
    }

    // Set up 30-45 components of targeted buffer : DG

    auto g_xz_xxxx_0 = pbuffer.data(idx_eri_0_dg + 30);

    auto g_xz_xxxy_0 = pbuffer.data(idx_eri_0_dg + 31);

    auto g_xz_xxxz_0 = pbuffer.data(idx_eri_0_dg + 32);

    auto g_xz_xxyy_0 = pbuffer.data(idx_eri_0_dg + 33);

    auto g_xz_xxyz_0 = pbuffer.data(idx_eri_0_dg + 34);

    auto g_xz_xxzz_0 = pbuffer.data(idx_eri_0_dg + 35);

    auto g_xz_xyyy_0 = pbuffer.data(idx_eri_0_dg + 36);

    auto g_xz_xyyz_0 = pbuffer.data(idx_eri_0_dg + 37);

    auto g_xz_xyzz_0 = pbuffer.data(idx_eri_0_dg + 38);

    auto g_xz_xzzz_0 = pbuffer.data(idx_eri_0_dg + 39);

    auto g_xz_yyyy_0 = pbuffer.data(idx_eri_0_dg + 40);

    auto g_xz_yyyz_0 = pbuffer.data(idx_eri_0_dg + 41);

    auto g_xz_yyzz_0 = pbuffer.data(idx_eri_0_dg + 42);

    auto g_xz_yzzz_0 = pbuffer.data(idx_eri_0_dg + 43);

    auto g_xz_zzzz_0 = pbuffer.data(idx_eri_0_dg + 44);

    #pragma omp simd aligned(g_x_xxxx_1, g_x_xxxy_1, g_x_xxyy_1, g_x_xyyy_1, g_xz_xxxx_0, g_xz_xxxy_0, g_xz_xxxz_0, g_xz_xxyy_0, g_xz_xxyz_0, g_xz_xxzz_0, g_xz_xyyy_0, g_xz_xyyz_0, g_xz_xyzz_0, g_xz_xzzz_0, g_xz_yyyy_0, g_xz_yyyz_0, g_xz_yyzz_0, g_xz_yzzz_0, g_xz_zzzz_0, g_z_xxxz_1, g_z_xxyz_1, g_z_xxz_1, g_z_xxzz_1, g_z_xyyz_1, g_z_xyz_1, g_z_xyzz_1, g_z_xzz_1, g_z_xzzz_1, g_z_yyyy_1, g_z_yyyz_1, g_z_yyz_1, g_z_yyzz_1, g_z_yzz_1, g_z_yzzz_1, g_z_zzz_1, g_z_zzzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xz_xxxx_0[i] = g_x_xxxx_1[i] * pa_z[i];

        g_xz_xxxy_0[i] = g_x_xxxy_1[i] * pa_z[i];

        g_xz_xxxz_0[i] = 3.0 * g_z_xxz_1[i] * fe_0 + g_z_xxxz_1[i] * pa_x[i];

        g_xz_xxyy_0[i] = g_x_xxyy_1[i] * pa_z[i];

        g_xz_xxyz_0[i] = 2.0 * g_z_xyz_1[i] * fe_0 + g_z_xxyz_1[i] * pa_x[i];

        g_xz_xxzz_0[i] = 2.0 * g_z_xzz_1[i] * fe_0 + g_z_xxzz_1[i] * pa_x[i];

        g_xz_xyyy_0[i] = g_x_xyyy_1[i] * pa_z[i];

        g_xz_xyyz_0[i] = g_z_yyz_1[i] * fe_0 + g_z_xyyz_1[i] * pa_x[i];

        g_xz_xyzz_0[i] = g_z_yzz_1[i] * fe_0 + g_z_xyzz_1[i] * pa_x[i];

        g_xz_xzzz_0[i] = g_z_zzz_1[i] * fe_0 + g_z_xzzz_1[i] * pa_x[i];

        g_xz_yyyy_0[i] = g_z_yyyy_1[i] * pa_x[i];

        g_xz_yyyz_0[i] = g_z_yyyz_1[i] * pa_x[i];

        g_xz_yyzz_0[i] = g_z_yyzz_1[i] * pa_x[i];

        g_xz_yzzz_0[i] = g_z_yzzz_1[i] * pa_x[i];

        g_xz_zzzz_0[i] = g_z_zzzz_1[i] * pa_x[i];
    }

    // Set up 45-60 components of targeted buffer : DG

    auto g_yy_xxxx_0 = pbuffer.data(idx_eri_0_dg + 45);

    auto g_yy_xxxy_0 = pbuffer.data(idx_eri_0_dg + 46);

    auto g_yy_xxxz_0 = pbuffer.data(idx_eri_0_dg + 47);

    auto g_yy_xxyy_0 = pbuffer.data(idx_eri_0_dg + 48);

    auto g_yy_xxyz_0 = pbuffer.data(idx_eri_0_dg + 49);

    auto g_yy_xxzz_0 = pbuffer.data(idx_eri_0_dg + 50);

    auto g_yy_xyyy_0 = pbuffer.data(idx_eri_0_dg + 51);

    auto g_yy_xyyz_0 = pbuffer.data(idx_eri_0_dg + 52);

    auto g_yy_xyzz_0 = pbuffer.data(idx_eri_0_dg + 53);

    auto g_yy_xzzz_0 = pbuffer.data(idx_eri_0_dg + 54);

    auto g_yy_yyyy_0 = pbuffer.data(idx_eri_0_dg + 55);

    auto g_yy_yyyz_0 = pbuffer.data(idx_eri_0_dg + 56);

    auto g_yy_yyzz_0 = pbuffer.data(idx_eri_0_dg + 57);

    auto g_yy_yzzz_0 = pbuffer.data(idx_eri_0_dg + 58);

    auto g_yy_zzzz_0 = pbuffer.data(idx_eri_0_dg + 59);

    #pragma omp simd aligned(g_0_xxxx_0, g_0_xxxx_1, g_0_xxxy_0, g_0_xxxy_1, g_0_xxxz_0, g_0_xxxz_1, g_0_xxyy_0, g_0_xxyy_1, g_0_xxyz_0, g_0_xxyz_1, g_0_xxzz_0, g_0_xxzz_1, g_0_xyyy_0, g_0_xyyy_1, g_0_xyyz_0, g_0_xyyz_1, g_0_xyzz_0, g_0_xyzz_1, g_0_xzzz_0, g_0_xzzz_1, g_0_yyyy_0, g_0_yyyy_1, g_0_yyyz_0, g_0_yyyz_1, g_0_yyzz_0, g_0_yyzz_1, g_0_yzzz_0, g_0_yzzz_1, g_0_zzzz_0, g_0_zzzz_1, g_y_xxx_1, g_y_xxxx_1, g_y_xxxy_1, g_y_xxxz_1, g_y_xxy_1, g_y_xxyy_1, g_y_xxyz_1, g_y_xxz_1, g_y_xxzz_1, g_y_xyy_1, g_y_xyyy_1, g_y_xyyz_1, g_y_xyz_1, g_y_xyzz_1, g_y_xzz_1, g_y_xzzz_1, g_y_yyy_1, g_y_yyyy_1, g_y_yyyz_1, g_y_yyz_1, g_y_yyzz_1, g_y_yzz_1, g_y_yzzz_1, g_y_zzz_1, g_y_zzzz_1, g_yy_xxxx_0, g_yy_xxxy_0, g_yy_xxxz_0, g_yy_xxyy_0, g_yy_xxyz_0, g_yy_xxzz_0, g_yy_xyyy_0, g_yy_xyyz_0, g_yy_xyzz_0, g_yy_xzzz_0, g_yy_yyyy_0, g_yy_yyyz_0, g_yy_yyzz_0, g_yy_yzzz_0, g_yy_zzzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yy_xxxx_0[i] = g_0_xxxx_0[i] * fbe_0 - g_0_xxxx_1[i] * fz_be_0 + g_y_xxxx_1[i] * pa_y[i];

        g_yy_xxxy_0[i] = g_0_xxxy_0[i] * fbe_0 - g_0_xxxy_1[i] * fz_be_0 + g_y_xxx_1[i] * fe_0 + g_y_xxxy_1[i] * pa_y[i];

        g_yy_xxxz_0[i] = g_0_xxxz_0[i] * fbe_0 - g_0_xxxz_1[i] * fz_be_0 + g_y_xxxz_1[i] * pa_y[i];

        g_yy_xxyy_0[i] = g_0_xxyy_0[i] * fbe_0 - g_0_xxyy_1[i] * fz_be_0 + 2.0 * g_y_xxy_1[i] * fe_0 + g_y_xxyy_1[i] * pa_y[i];

        g_yy_xxyz_0[i] = g_0_xxyz_0[i] * fbe_0 - g_0_xxyz_1[i] * fz_be_0 + g_y_xxz_1[i] * fe_0 + g_y_xxyz_1[i] * pa_y[i];

        g_yy_xxzz_0[i] = g_0_xxzz_0[i] * fbe_0 - g_0_xxzz_1[i] * fz_be_0 + g_y_xxzz_1[i] * pa_y[i];

        g_yy_xyyy_0[i] = g_0_xyyy_0[i] * fbe_0 - g_0_xyyy_1[i] * fz_be_0 + 3.0 * g_y_xyy_1[i] * fe_0 + g_y_xyyy_1[i] * pa_y[i];

        g_yy_xyyz_0[i] = g_0_xyyz_0[i] * fbe_0 - g_0_xyyz_1[i] * fz_be_0 + 2.0 * g_y_xyz_1[i] * fe_0 + g_y_xyyz_1[i] * pa_y[i];

        g_yy_xyzz_0[i] = g_0_xyzz_0[i] * fbe_0 - g_0_xyzz_1[i] * fz_be_0 + g_y_xzz_1[i] * fe_0 + g_y_xyzz_1[i] * pa_y[i];

        g_yy_xzzz_0[i] = g_0_xzzz_0[i] * fbe_0 - g_0_xzzz_1[i] * fz_be_0 + g_y_xzzz_1[i] * pa_y[i];

        g_yy_yyyy_0[i] = g_0_yyyy_0[i] * fbe_0 - g_0_yyyy_1[i] * fz_be_0 + 4.0 * g_y_yyy_1[i] * fe_0 + g_y_yyyy_1[i] * pa_y[i];

        g_yy_yyyz_0[i] = g_0_yyyz_0[i] * fbe_0 - g_0_yyyz_1[i] * fz_be_0 + 3.0 * g_y_yyz_1[i] * fe_0 + g_y_yyyz_1[i] * pa_y[i];

        g_yy_yyzz_0[i] = g_0_yyzz_0[i] * fbe_0 - g_0_yyzz_1[i] * fz_be_0 + 2.0 * g_y_yzz_1[i] * fe_0 + g_y_yyzz_1[i] * pa_y[i];

        g_yy_yzzz_0[i] = g_0_yzzz_0[i] * fbe_0 - g_0_yzzz_1[i] * fz_be_0 + g_y_zzz_1[i] * fe_0 + g_y_yzzz_1[i] * pa_y[i];

        g_yy_zzzz_0[i] = g_0_zzzz_0[i] * fbe_0 - g_0_zzzz_1[i] * fz_be_0 + g_y_zzzz_1[i] * pa_y[i];
    }

    // Set up 60-75 components of targeted buffer : DG

    auto g_yz_xxxx_0 = pbuffer.data(idx_eri_0_dg + 60);

    auto g_yz_xxxy_0 = pbuffer.data(idx_eri_0_dg + 61);

    auto g_yz_xxxz_0 = pbuffer.data(idx_eri_0_dg + 62);

    auto g_yz_xxyy_0 = pbuffer.data(idx_eri_0_dg + 63);

    auto g_yz_xxyz_0 = pbuffer.data(idx_eri_0_dg + 64);

    auto g_yz_xxzz_0 = pbuffer.data(idx_eri_0_dg + 65);

    auto g_yz_xyyy_0 = pbuffer.data(idx_eri_0_dg + 66);

    auto g_yz_xyyz_0 = pbuffer.data(idx_eri_0_dg + 67);

    auto g_yz_xyzz_0 = pbuffer.data(idx_eri_0_dg + 68);

    auto g_yz_xzzz_0 = pbuffer.data(idx_eri_0_dg + 69);

    auto g_yz_yyyy_0 = pbuffer.data(idx_eri_0_dg + 70);

    auto g_yz_yyyz_0 = pbuffer.data(idx_eri_0_dg + 71);

    auto g_yz_yyzz_0 = pbuffer.data(idx_eri_0_dg + 72);

    auto g_yz_yzzz_0 = pbuffer.data(idx_eri_0_dg + 73);

    auto g_yz_zzzz_0 = pbuffer.data(idx_eri_0_dg + 74);

    #pragma omp simd aligned(g_y_xxxy_1, g_y_xxyy_1, g_y_xyyy_1, g_y_yyyy_1, g_yz_xxxx_0, g_yz_xxxy_0, g_yz_xxxz_0, g_yz_xxyy_0, g_yz_xxyz_0, g_yz_xxzz_0, g_yz_xyyy_0, g_yz_xyyz_0, g_yz_xyzz_0, g_yz_xzzz_0, g_yz_yyyy_0, g_yz_yyyz_0, g_yz_yyzz_0, g_yz_yzzz_0, g_yz_zzzz_0, g_z_xxxx_1, g_z_xxxz_1, g_z_xxyz_1, g_z_xxz_1, g_z_xxzz_1, g_z_xyyz_1, g_z_xyz_1, g_z_xyzz_1, g_z_xzz_1, g_z_xzzz_1, g_z_yyyz_1, g_z_yyz_1, g_z_yyzz_1, g_z_yzz_1, g_z_yzzz_1, g_z_zzz_1, g_z_zzzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yz_xxxx_0[i] = g_z_xxxx_1[i] * pa_y[i];

        g_yz_xxxy_0[i] = g_y_xxxy_1[i] * pa_z[i];

        g_yz_xxxz_0[i] = g_z_xxxz_1[i] * pa_y[i];

        g_yz_xxyy_0[i] = g_y_xxyy_1[i] * pa_z[i];

        g_yz_xxyz_0[i] = g_z_xxz_1[i] * fe_0 + g_z_xxyz_1[i] * pa_y[i];

        g_yz_xxzz_0[i] = g_z_xxzz_1[i] * pa_y[i];

        g_yz_xyyy_0[i] = g_y_xyyy_1[i] * pa_z[i];

        g_yz_xyyz_0[i] = 2.0 * g_z_xyz_1[i] * fe_0 + g_z_xyyz_1[i] * pa_y[i];

        g_yz_xyzz_0[i] = g_z_xzz_1[i] * fe_0 + g_z_xyzz_1[i] * pa_y[i];

        g_yz_xzzz_0[i] = g_z_xzzz_1[i] * pa_y[i];

        g_yz_yyyy_0[i] = g_y_yyyy_1[i] * pa_z[i];

        g_yz_yyyz_0[i] = 3.0 * g_z_yyz_1[i] * fe_0 + g_z_yyyz_1[i] * pa_y[i];

        g_yz_yyzz_0[i] = 2.0 * g_z_yzz_1[i] * fe_0 + g_z_yyzz_1[i] * pa_y[i];

        g_yz_yzzz_0[i] = g_z_zzz_1[i] * fe_0 + g_z_yzzz_1[i] * pa_y[i];

        g_yz_zzzz_0[i] = g_z_zzzz_1[i] * pa_y[i];
    }

    // Set up 75-90 components of targeted buffer : DG

    auto g_zz_xxxx_0 = pbuffer.data(idx_eri_0_dg + 75);

    auto g_zz_xxxy_0 = pbuffer.data(idx_eri_0_dg + 76);

    auto g_zz_xxxz_0 = pbuffer.data(idx_eri_0_dg + 77);

    auto g_zz_xxyy_0 = pbuffer.data(idx_eri_0_dg + 78);

    auto g_zz_xxyz_0 = pbuffer.data(idx_eri_0_dg + 79);

    auto g_zz_xxzz_0 = pbuffer.data(idx_eri_0_dg + 80);

    auto g_zz_xyyy_0 = pbuffer.data(idx_eri_0_dg + 81);

    auto g_zz_xyyz_0 = pbuffer.data(idx_eri_0_dg + 82);

    auto g_zz_xyzz_0 = pbuffer.data(idx_eri_0_dg + 83);

    auto g_zz_xzzz_0 = pbuffer.data(idx_eri_0_dg + 84);

    auto g_zz_yyyy_0 = pbuffer.data(idx_eri_0_dg + 85);

    auto g_zz_yyyz_0 = pbuffer.data(idx_eri_0_dg + 86);

    auto g_zz_yyzz_0 = pbuffer.data(idx_eri_0_dg + 87);

    auto g_zz_yzzz_0 = pbuffer.data(idx_eri_0_dg + 88);

    auto g_zz_zzzz_0 = pbuffer.data(idx_eri_0_dg + 89);

    #pragma omp simd aligned(g_0_xxxx_0, g_0_xxxx_1, g_0_xxxy_0, g_0_xxxy_1, g_0_xxxz_0, g_0_xxxz_1, g_0_xxyy_0, g_0_xxyy_1, g_0_xxyz_0, g_0_xxyz_1, g_0_xxzz_0, g_0_xxzz_1, g_0_xyyy_0, g_0_xyyy_1, g_0_xyyz_0, g_0_xyyz_1, g_0_xyzz_0, g_0_xyzz_1, g_0_xzzz_0, g_0_xzzz_1, g_0_yyyy_0, g_0_yyyy_1, g_0_yyyz_0, g_0_yyyz_1, g_0_yyzz_0, g_0_yyzz_1, g_0_yzzz_0, g_0_yzzz_1, g_0_zzzz_0, g_0_zzzz_1, g_z_xxx_1, g_z_xxxx_1, g_z_xxxy_1, g_z_xxxz_1, g_z_xxy_1, g_z_xxyy_1, g_z_xxyz_1, g_z_xxz_1, g_z_xxzz_1, g_z_xyy_1, g_z_xyyy_1, g_z_xyyz_1, g_z_xyz_1, g_z_xyzz_1, g_z_xzz_1, g_z_xzzz_1, g_z_yyy_1, g_z_yyyy_1, g_z_yyyz_1, g_z_yyz_1, g_z_yyzz_1, g_z_yzz_1, g_z_yzzz_1, g_z_zzz_1, g_z_zzzz_1, g_zz_xxxx_0, g_zz_xxxy_0, g_zz_xxxz_0, g_zz_xxyy_0, g_zz_xxyz_0, g_zz_xxzz_0, g_zz_xyyy_0, g_zz_xyyz_0, g_zz_xyzz_0, g_zz_xzzz_0, g_zz_yyyy_0, g_zz_yyyz_0, g_zz_yyzz_0, g_zz_yzzz_0, g_zz_zzzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zz_xxxx_0[i] = g_0_xxxx_0[i] * fbe_0 - g_0_xxxx_1[i] * fz_be_0 + g_z_xxxx_1[i] * pa_z[i];

        g_zz_xxxy_0[i] = g_0_xxxy_0[i] * fbe_0 - g_0_xxxy_1[i] * fz_be_0 + g_z_xxxy_1[i] * pa_z[i];

        g_zz_xxxz_0[i] = g_0_xxxz_0[i] * fbe_0 - g_0_xxxz_1[i] * fz_be_0 + g_z_xxx_1[i] * fe_0 + g_z_xxxz_1[i] * pa_z[i];

        g_zz_xxyy_0[i] = g_0_xxyy_0[i] * fbe_0 - g_0_xxyy_1[i] * fz_be_0 + g_z_xxyy_1[i] * pa_z[i];

        g_zz_xxyz_0[i] = g_0_xxyz_0[i] * fbe_0 - g_0_xxyz_1[i] * fz_be_0 + g_z_xxy_1[i] * fe_0 + g_z_xxyz_1[i] * pa_z[i];

        g_zz_xxzz_0[i] = g_0_xxzz_0[i] * fbe_0 - g_0_xxzz_1[i] * fz_be_0 + 2.0 * g_z_xxz_1[i] * fe_0 + g_z_xxzz_1[i] * pa_z[i];

        g_zz_xyyy_0[i] = g_0_xyyy_0[i] * fbe_0 - g_0_xyyy_1[i] * fz_be_0 + g_z_xyyy_1[i] * pa_z[i];

        g_zz_xyyz_0[i] = g_0_xyyz_0[i] * fbe_0 - g_0_xyyz_1[i] * fz_be_0 + g_z_xyy_1[i] * fe_0 + g_z_xyyz_1[i] * pa_z[i];

        g_zz_xyzz_0[i] = g_0_xyzz_0[i] * fbe_0 - g_0_xyzz_1[i] * fz_be_0 + 2.0 * g_z_xyz_1[i] * fe_0 + g_z_xyzz_1[i] * pa_z[i];

        g_zz_xzzz_0[i] = g_0_xzzz_0[i] * fbe_0 - g_0_xzzz_1[i] * fz_be_0 + 3.0 * g_z_xzz_1[i] * fe_0 + g_z_xzzz_1[i] * pa_z[i];

        g_zz_yyyy_0[i] = g_0_yyyy_0[i] * fbe_0 - g_0_yyyy_1[i] * fz_be_0 + g_z_yyyy_1[i] * pa_z[i];

        g_zz_yyyz_0[i] = g_0_yyyz_0[i] * fbe_0 - g_0_yyyz_1[i] * fz_be_0 + g_z_yyy_1[i] * fe_0 + g_z_yyyz_1[i] * pa_z[i];

        g_zz_yyzz_0[i] = g_0_yyzz_0[i] * fbe_0 - g_0_yyzz_1[i] * fz_be_0 + 2.0 * g_z_yyz_1[i] * fe_0 + g_z_yyzz_1[i] * pa_z[i];

        g_zz_yzzz_0[i] = g_0_yzzz_0[i] * fbe_0 - g_0_yzzz_1[i] * fz_be_0 + 3.0 * g_z_yzz_1[i] * fe_0 + g_z_yzzz_1[i] * pa_z[i];

        g_zz_zzzz_0[i] = g_0_zzzz_0[i] * fbe_0 - g_0_zzzz_1[i] * fz_be_0 + 4.0 * g_z_zzz_1[i] * fe_0 + g_z_zzzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

