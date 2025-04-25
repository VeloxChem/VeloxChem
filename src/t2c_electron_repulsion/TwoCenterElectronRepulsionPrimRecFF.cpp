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

#include "TwoCenterElectronRepulsionPrimRecFF.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_ff(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_ff,
                                const size_t idx_eri_0_pf,
                                const size_t idx_eri_1_pf,
                                const size_t idx_eri_1_dd,
                                const size_t idx_eri_1_df,
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

    // Set up components of auxiliary buffer : PF

    auto g_x_xxx_0 = pbuffer.data(idx_eri_0_pf);

    auto g_x_xxy_0 = pbuffer.data(idx_eri_0_pf + 1);

    auto g_x_xxz_0 = pbuffer.data(idx_eri_0_pf + 2);

    auto g_x_xyy_0 = pbuffer.data(idx_eri_0_pf + 3);

    auto g_x_xyz_0 = pbuffer.data(idx_eri_0_pf + 4);

    auto g_x_xzz_0 = pbuffer.data(idx_eri_0_pf + 5);

    auto g_x_yyy_0 = pbuffer.data(idx_eri_0_pf + 6);

    auto g_x_yyz_0 = pbuffer.data(idx_eri_0_pf + 7);

    auto g_x_yzz_0 = pbuffer.data(idx_eri_0_pf + 8);

    auto g_x_zzz_0 = pbuffer.data(idx_eri_0_pf + 9);

    auto g_y_xxx_0 = pbuffer.data(idx_eri_0_pf + 10);

    auto g_y_xxy_0 = pbuffer.data(idx_eri_0_pf + 11);

    auto g_y_xxz_0 = pbuffer.data(idx_eri_0_pf + 12);

    auto g_y_xyy_0 = pbuffer.data(idx_eri_0_pf + 13);

    auto g_y_xyz_0 = pbuffer.data(idx_eri_0_pf + 14);

    auto g_y_xzz_0 = pbuffer.data(idx_eri_0_pf + 15);

    auto g_y_yyy_0 = pbuffer.data(idx_eri_0_pf + 16);

    auto g_y_yyz_0 = pbuffer.data(idx_eri_0_pf + 17);

    auto g_y_yzz_0 = pbuffer.data(idx_eri_0_pf + 18);

    auto g_y_zzz_0 = pbuffer.data(idx_eri_0_pf + 19);

    auto g_z_xxx_0 = pbuffer.data(idx_eri_0_pf + 20);

    auto g_z_xxy_0 = pbuffer.data(idx_eri_0_pf + 21);

    auto g_z_xxz_0 = pbuffer.data(idx_eri_0_pf + 22);

    auto g_z_xyy_0 = pbuffer.data(idx_eri_0_pf + 23);

    auto g_z_xyz_0 = pbuffer.data(idx_eri_0_pf + 24);

    auto g_z_xzz_0 = pbuffer.data(idx_eri_0_pf + 25);

    auto g_z_yyy_0 = pbuffer.data(idx_eri_0_pf + 26);

    auto g_z_yyz_0 = pbuffer.data(idx_eri_0_pf + 27);

    auto g_z_yzz_0 = pbuffer.data(idx_eri_0_pf + 28);

    auto g_z_zzz_0 = pbuffer.data(idx_eri_0_pf + 29);

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

    // Set up components of auxiliary buffer : DD

    auto g_xx_xx_1 = pbuffer.data(idx_eri_1_dd);

    auto g_xx_xy_1 = pbuffer.data(idx_eri_1_dd + 1);

    auto g_xx_xz_1 = pbuffer.data(idx_eri_1_dd + 2);

    auto g_xx_yy_1 = pbuffer.data(idx_eri_1_dd + 3);

    auto g_xx_yz_1 = pbuffer.data(idx_eri_1_dd + 4);

    auto g_xx_zz_1 = pbuffer.data(idx_eri_1_dd + 5);

    auto g_yy_xx_1 = pbuffer.data(idx_eri_1_dd + 18);

    auto g_yy_xy_1 = pbuffer.data(idx_eri_1_dd + 19);

    auto g_yy_xz_1 = pbuffer.data(idx_eri_1_dd + 20);

    auto g_yy_yy_1 = pbuffer.data(idx_eri_1_dd + 21);

    auto g_yy_yz_1 = pbuffer.data(idx_eri_1_dd + 22);

    auto g_yy_zz_1 = pbuffer.data(idx_eri_1_dd + 23);

    auto g_yz_yz_1 = pbuffer.data(idx_eri_1_dd + 28);

    auto g_zz_xx_1 = pbuffer.data(idx_eri_1_dd + 30);

    auto g_zz_xy_1 = pbuffer.data(idx_eri_1_dd + 31);

    auto g_zz_xz_1 = pbuffer.data(idx_eri_1_dd + 32);

    auto g_zz_yy_1 = pbuffer.data(idx_eri_1_dd + 33);

    auto g_zz_yz_1 = pbuffer.data(idx_eri_1_dd + 34);

    auto g_zz_zz_1 = pbuffer.data(idx_eri_1_dd + 35);

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

    auto g_xy_xxy_1 = pbuffer.data(idx_eri_1_df + 11);

    auto g_xy_xyy_1 = pbuffer.data(idx_eri_1_df + 13);

    auto g_xz_xxx_1 = pbuffer.data(idx_eri_1_df + 20);

    auto g_xz_xxz_1 = pbuffer.data(idx_eri_1_df + 22);

    auto g_xz_xzz_1 = pbuffer.data(idx_eri_1_df + 25);

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

    auto g_yz_yyy_1 = pbuffer.data(idx_eri_1_df + 46);

    auto g_yz_yyz_1 = pbuffer.data(idx_eri_1_df + 47);

    auto g_yz_yzz_1 = pbuffer.data(idx_eri_1_df + 48);

    auto g_yz_zzz_1 = pbuffer.data(idx_eri_1_df + 49);

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

    // Set up 0-10 components of targeted buffer : FF

    auto g_xxx_xxx_0 = pbuffer.data(idx_eri_0_ff);

    auto g_xxx_xxy_0 = pbuffer.data(idx_eri_0_ff + 1);

    auto g_xxx_xxz_0 = pbuffer.data(idx_eri_0_ff + 2);

    auto g_xxx_xyy_0 = pbuffer.data(idx_eri_0_ff + 3);

    auto g_xxx_xyz_0 = pbuffer.data(idx_eri_0_ff + 4);

    auto g_xxx_xzz_0 = pbuffer.data(idx_eri_0_ff + 5);

    auto g_xxx_yyy_0 = pbuffer.data(idx_eri_0_ff + 6);

    auto g_xxx_yyz_0 = pbuffer.data(idx_eri_0_ff + 7);

    auto g_xxx_yzz_0 = pbuffer.data(idx_eri_0_ff + 8);

    auto g_xxx_zzz_0 = pbuffer.data(idx_eri_0_ff + 9);

    #pragma omp simd aligned(g_x_xxx_0, g_x_xxx_1, g_x_xxy_0, g_x_xxy_1, g_x_xxz_0, g_x_xxz_1, g_x_xyy_0, g_x_xyy_1, g_x_xyz_0, g_x_xyz_1, g_x_xzz_0, g_x_xzz_1, g_x_yyy_0, g_x_yyy_1, g_x_yyz_0, g_x_yyz_1, g_x_yzz_0, g_x_yzz_1, g_x_zzz_0, g_x_zzz_1, g_xx_xx_1, g_xx_xxx_1, g_xx_xxy_1, g_xx_xxz_1, g_xx_xy_1, g_xx_xyy_1, g_xx_xyz_1, g_xx_xz_1, g_xx_xzz_1, g_xx_yy_1, g_xx_yyy_1, g_xx_yyz_1, g_xx_yz_1, g_xx_yzz_1, g_xx_zz_1, g_xx_zzz_1, g_xxx_xxx_0, g_xxx_xxy_0, g_xxx_xxz_0, g_xxx_xyy_0, g_xxx_xyz_0, g_xxx_xzz_0, g_xxx_yyy_0, g_xxx_yyz_0, g_xxx_yzz_0, g_xxx_zzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxx_xxx_0[i] = 2.0 * g_x_xxx_0[i] * fbe_0 - 2.0 * g_x_xxx_1[i] * fz_be_0 + 3.0 * g_xx_xx_1[i] * fe_0 + g_xx_xxx_1[i] * pa_x[i];

        g_xxx_xxy_0[i] = 2.0 * g_x_xxy_0[i] * fbe_0 - 2.0 * g_x_xxy_1[i] * fz_be_0 + 2.0 * g_xx_xy_1[i] * fe_0 + g_xx_xxy_1[i] * pa_x[i];

        g_xxx_xxz_0[i] = 2.0 * g_x_xxz_0[i] * fbe_0 - 2.0 * g_x_xxz_1[i] * fz_be_0 + 2.0 * g_xx_xz_1[i] * fe_0 + g_xx_xxz_1[i] * pa_x[i];

        g_xxx_xyy_0[i] = 2.0 * g_x_xyy_0[i] * fbe_0 - 2.0 * g_x_xyy_1[i] * fz_be_0 + g_xx_yy_1[i] * fe_0 + g_xx_xyy_1[i] * pa_x[i];

        g_xxx_xyz_0[i] = 2.0 * g_x_xyz_0[i] * fbe_0 - 2.0 * g_x_xyz_1[i] * fz_be_0 + g_xx_yz_1[i] * fe_0 + g_xx_xyz_1[i] * pa_x[i];

        g_xxx_xzz_0[i] = 2.0 * g_x_xzz_0[i] * fbe_0 - 2.0 * g_x_xzz_1[i] * fz_be_0 + g_xx_zz_1[i] * fe_0 + g_xx_xzz_1[i] * pa_x[i];

        g_xxx_yyy_0[i] = 2.0 * g_x_yyy_0[i] * fbe_0 - 2.0 * g_x_yyy_1[i] * fz_be_0 + g_xx_yyy_1[i] * pa_x[i];

        g_xxx_yyz_0[i] = 2.0 * g_x_yyz_0[i] * fbe_0 - 2.0 * g_x_yyz_1[i] * fz_be_0 + g_xx_yyz_1[i] * pa_x[i];

        g_xxx_yzz_0[i] = 2.0 * g_x_yzz_0[i] * fbe_0 - 2.0 * g_x_yzz_1[i] * fz_be_0 + g_xx_yzz_1[i] * pa_x[i];

        g_xxx_zzz_0[i] = 2.0 * g_x_zzz_0[i] * fbe_0 - 2.0 * g_x_zzz_1[i] * fz_be_0 + g_xx_zzz_1[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : FF

    auto g_xxy_xxx_0 = pbuffer.data(idx_eri_0_ff + 10);

    auto g_xxy_xxy_0 = pbuffer.data(idx_eri_0_ff + 11);

    auto g_xxy_xxz_0 = pbuffer.data(idx_eri_0_ff + 12);

    auto g_xxy_xyy_0 = pbuffer.data(idx_eri_0_ff + 13);

    auto g_xxy_xyz_0 = pbuffer.data(idx_eri_0_ff + 14);

    auto g_xxy_xzz_0 = pbuffer.data(idx_eri_0_ff + 15);

    auto g_xxy_yyy_0 = pbuffer.data(idx_eri_0_ff + 16);

    auto g_xxy_yyz_0 = pbuffer.data(idx_eri_0_ff + 17);

    auto g_xxy_yzz_0 = pbuffer.data(idx_eri_0_ff + 18);

    auto g_xxy_zzz_0 = pbuffer.data(idx_eri_0_ff + 19);

    #pragma omp simd aligned(g_xx_xx_1, g_xx_xxx_1, g_xx_xxy_1, g_xx_xxz_1, g_xx_xy_1, g_xx_xyy_1, g_xx_xyz_1, g_xx_xz_1, g_xx_xzz_1, g_xx_yy_1, g_xx_yyy_1, g_xx_yyz_1, g_xx_yz_1, g_xx_yzz_1, g_xx_zz_1, g_xx_zzz_1, g_xxy_xxx_0, g_xxy_xxy_0, g_xxy_xxz_0, g_xxy_xyy_0, g_xxy_xyz_0, g_xxy_xzz_0, g_xxy_yyy_0, g_xxy_yyz_0, g_xxy_yzz_0, g_xxy_zzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxy_xxx_0[i] = g_xx_xxx_1[i] * pa_y[i];

        g_xxy_xxy_0[i] = g_xx_xx_1[i] * fe_0 + g_xx_xxy_1[i] * pa_y[i];

        g_xxy_xxz_0[i] = g_xx_xxz_1[i] * pa_y[i];

        g_xxy_xyy_0[i] = 2.0 * g_xx_xy_1[i] * fe_0 + g_xx_xyy_1[i] * pa_y[i];

        g_xxy_xyz_0[i] = g_xx_xz_1[i] * fe_0 + g_xx_xyz_1[i] * pa_y[i];

        g_xxy_xzz_0[i] = g_xx_xzz_1[i] * pa_y[i];

        g_xxy_yyy_0[i] = 3.0 * g_xx_yy_1[i] * fe_0 + g_xx_yyy_1[i] * pa_y[i];

        g_xxy_yyz_0[i] = 2.0 * g_xx_yz_1[i] * fe_0 + g_xx_yyz_1[i] * pa_y[i];

        g_xxy_yzz_0[i] = g_xx_zz_1[i] * fe_0 + g_xx_yzz_1[i] * pa_y[i];

        g_xxy_zzz_0[i] = g_xx_zzz_1[i] * pa_y[i];
    }

    // Set up 20-30 components of targeted buffer : FF

    auto g_xxz_xxx_0 = pbuffer.data(idx_eri_0_ff + 20);

    auto g_xxz_xxy_0 = pbuffer.data(idx_eri_0_ff + 21);

    auto g_xxz_xxz_0 = pbuffer.data(idx_eri_0_ff + 22);

    auto g_xxz_xyy_0 = pbuffer.data(idx_eri_0_ff + 23);

    auto g_xxz_xyz_0 = pbuffer.data(idx_eri_0_ff + 24);

    auto g_xxz_xzz_0 = pbuffer.data(idx_eri_0_ff + 25);

    auto g_xxz_yyy_0 = pbuffer.data(idx_eri_0_ff + 26);

    auto g_xxz_yyz_0 = pbuffer.data(idx_eri_0_ff + 27);

    auto g_xxz_yzz_0 = pbuffer.data(idx_eri_0_ff + 28);

    auto g_xxz_zzz_0 = pbuffer.data(idx_eri_0_ff + 29);

    #pragma omp simd aligned(g_xx_xx_1, g_xx_xxx_1, g_xx_xxy_1, g_xx_xxz_1, g_xx_xy_1, g_xx_xyy_1, g_xx_xyz_1, g_xx_xz_1, g_xx_xzz_1, g_xx_yy_1, g_xx_yyy_1, g_xx_yyz_1, g_xx_yz_1, g_xx_yzz_1, g_xx_zz_1, g_xx_zzz_1, g_xxz_xxx_0, g_xxz_xxy_0, g_xxz_xxz_0, g_xxz_xyy_0, g_xxz_xyz_0, g_xxz_xzz_0, g_xxz_yyy_0, g_xxz_yyz_0, g_xxz_yzz_0, g_xxz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxz_xxx_0[i] = g_xx_xxx_1[i] * pa_z[i];

        g_xxz_xxy_0[i] = g_xx_xxy_1[i] * pa_z[i];

        g_xxz_xxz_0[i] = g_xx_xx_1[i] * fe_0 + g_xx_xxz_1[i] * pa_z[i];

        g_xxz_xyy_0[i] = g_xx_xyy_1[i] * pa_z[i];

        g_xxz_xyz_0[i] = g_xx_xy_1[i] * fe_0 + g_xx_xyz_1[i] * pa_z[i];

        g_xxz_xzz_0[i] = 2.0 * g_xx_xz_1[i] * fe_0 + g_xx_xzz_1[i] * pa_z[i];

        g_xxz_yyy_0[i] = g_xx_yyy_1[i] * pa_z[i];

        g_xxz_yyz_0[i] = g_xx_yy_1[i] * fe_0 + g_xx_yyz_1[i] * pa_z[i];

        g_xxz_yzz_0[i] = 2.0 * g_xx_yz_1[i] * fe_0 + g_xx_yzz_1[i] * pa_z[i];

        g_xxz_zzz_0[i] = 3.0 * g_xx_zz_1[i] * fe_0 + g_xx_zzz_1[i] * pa_z[i];
    }

    // Set up 30-40 components of targeted buffer : FF

    auto g_xyy_xxx_0 = pbuffer.data(idx_eri_0_ff + 30);

    auto g_xyy_xxy_0 = pbuffer.data(idx_eri_0_ff + 31);

    auto g_xyy_xxz_0 = pbuffer.data(idx_eri_0_ff + 32);

    auto g_xyy_xyy_0 = pbuffer.data(idx_eri_0_ff + 33);

    auto g_xyy_xyz_0 = pbuffer.data(idx_eri_0_ff + 34);

    auto g_xyy_xzz_0 = pbuffer.data(idx_eri_0_ff + 35);

    auto g_xyy_yyy_0 = pbuffer.data(idx_eri_0_ff + 36);

    auto g_xyy_yyz_0 = pbuffer.data(idx_eri_0_ff + 37);

    auto g_xyy_yzz_0 = pbuffer.data(idx_eri_0_ff + 38);

    auto g_xyy_zzz_0 = pbuffer.data(idx_eri_0_ff + 39);

    #pragma omp simd aligned(g_xyy_xxx_0, g_xyy_xxy_0, g_xyy_xxz_0, g_xyy_xyy_0, g_xyy_xyz_0, g_xyy_xzz_0, g_xyy_yyy_0, g_xyy_yyz_0, g_xyy_yzz_0, g_xyy_zzz_0, g_yy_xx_1, g_yy_xxx_1, g_yy_xxy_1, g_yy_xxz_1, g_yy_xy_1, g_yy_xyy_1, g_yy_xyz_1, g_yy_xz_1, g_yy_xzz_1, g_yy_yy_1, g_yy_yyy_1, g_yy_yyz_1, g_yy_yz_1, g_yy_yzz_1, g_yy_zz_1, g_yy_zzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyy_xxx_0[i] = 3.0 * g_yy_xx_1[i] * fe_0 + g_yy_xxx_1[i] * pa_x[i];

        g_xyy_xxy_0[i] = 2.0 * g_yy_xy_1[i] * fe_0 + g_yy_xxy_1[i] * pa_x[i];

        g_xyy_xxz_0[i] = 2.0 * g_yy_xz_1[i] * fe_0 + g_yy_xxz_1[i] * pa_x[i];

        g_xyy_xyy_0[i] = g_yy_yy_1[i] * fe_0 + g_yy_xyy_1[i] * pa_x[i];

        g_xyy_xyz_0[i] = g_yy_yz_1[i] * fe_0 + g_yy_xyz_1[i] * pa_x[i];

        g_xyy_xzz_0[i] = g_yy_zz_1[i] * fe_0 + g_yy_xzz_1[i] * pa_x[i];

        g_xyy_yyy_0[i] = g_yy_yyy_1[i] * pa_x[i];

        g_xyy_yyz_0[i] = g_yy_yyz_1[i] * pa_x[i];

        g_xyy_yzz_0[i] = g_yy_yzz_1[i] * pa_x[i];

        g_xyy_zzz_0[i] = g_yy_zzz_1[i] * pa_x[i];
    }

    // Set up 40-50 components of targeted buffer : FF

    auto g_xyz_xxx_0 = pbuffer.data(idx_eri_0_ff + 40);

    auto g_xyz_xxy_0 = pbuffer.data(idx_eri_0_ff + 41);

    auto g_xyz_xxz_0 = pbuffer.data(idx_eri_0_ff + 42);

    auto g_xyz_xyy_0 = pbuffer.data(idx_eri_0_ff + 43);

    auto g_xyz_xyz_0 = pbuffer.data(idx_eri_0_ff + 44);

    auto g_xyz_xzz_0 = pbuffer.data(idx_eri_0_ff + 45);

    auto g_xyz_yyy_0 = pbuffer.data(idx_eri_0_ff + 46);

    auto g_xyz_yyz_0 = pbuffer.data(idx_eri_0_ff + 47);

    auto g_xyz_yzz_0 = pbuffer.data(idx_eri_0_ff + 48);

    auto g_xyz_zzz_0 = pbuffer.data(idx_eri_0_ff + 49);

    #pragma omp simd aligned(g_xy_xxy_1, g_xy_xyy_1, g_xyz_xxx_0, g_xyz_xxy_0, g_xyz_xxz_0, g_xyz_xyy_0, g_xyz_xyz_0, g_xyz_xzz_0, g_xyz_yyy_0, g_xyz_yyz_0, g_xyz_yzz_0, g_xyz_zzz_0, g_xz_xxx_1, g_xz_xxz_1, g_xz_xzz_1, g_yz_xyz_1, g_yz_yyy_1, g_yz_yyz_1, g_yz_yz_1, g_yz_yzz_1, g_yz_zzz_1, pa_x, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyz_xxx_0[i] = g_xz_xxx_1[i] * pa_y[i];

        g_xyz_xxy_0[i] = g_xy_xxy_1[i] * pa_z[i];

        g_xyz_xxz_0[i] = g_xz_xxz_1[i] * pa_y[i];

        g_xyz_xyy_0[i] = g_xy_xyy_1[i] * pa_z[i];

        g_xyz_xyz_0[i] = g_yz_yz_1[i] * fe_0 + g_yz_xyz_1[i] * pa_x[i];

        g_xyz_xzz_0[i] = g_xz_xzz_1[i] * pa_y[i];

        g_xyz_yyy_0[i] = g_yz_yyy_1[i] * pa_x[i];

        g_xyz_yyz_0[i] = g_yz_yyz_1[i] * pa_x[i];

        g_xyz_yzz_0[i] = g_yz_yzz_1[i] * pa_x[i];

        g_xyz_zzz_0[i] = g_yz_zzz_1[i] * pa_x[i];
    }

    // Set up 50-60 components of targeted buffer : FF

    auto g_xzz_xxx_0 = pbuffer.data(idx_eri_0_ff + 50);

    auto g_xzz_xxy_0 = pbuffer.data(idx_eri_0_ff + 51);

    auto g_xzz_xxz_0 = pbuffer.data(idx_eri_0_ff + 52);

    auto g_xzz_xyy_0 = pbuffer.data(idx_eri_0_ff + 53);

    auto g_xzz_xyz_0 = pbuffer.data(idx_eri_0_ff + 54);

    auto g_xzz_xzz_0 = pbuffer.data(idx_eri_0_ff + 55);

    auto g_xzz_yyy_0 = pbuffer.data(idx_eri_0_ff + 56);

    auto g_xzz_yyz_0 = pbuffer.data(idx_eri_0_ff + 57);

    auto g_xzz_yzz_0 = pbuffer.data(idx_eri_0_ff + 58);

    auto g_xzz_zzz_0 = pbuffer.data(idx_eri_0_ff + 59);

    #pragma omp simd aligned(g_xzz_xxx_0, g_xzz_xxy_0, g_xzz_xxz_0, g_xzz_xyy_0, g_xzz_xyz_0, g_xzz_xzz_0, g_xzz_yyy_0, g_xzz_yyz_0, g_xzz_yzz_0, g_xzz_zzz_0, g_zz_xx_1, g_zz_xxx_1, g_zz_xxy_1, g_zz_xxz_1, g_zz_xy_1, g_zz_xyy_1, g_zz_xyz_1, g_zz_xz_1, g_zz_xzz_1, g_zz_yy_1, g_zz_yyy_1, g_zz_yyz_1, g_zz_yz_1, g_zz_yzz_1, g_zz_zz_1, g_zz_zzz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzz_xxx_0[i] = 3.0 * g_zz_xx_1[i] * fe_0 + g_zz_xxx_1[i] * pa_x[i];

        g_xzz_xxy_0[i] = 2.0 * g_zz_xy_1[i] * fe_0 + g_zz_xxy_1[i] * pa_x[i];

        g_xzz_xxz_0[i] = 2.0 * g_zz_xz_1[i] * fe_0 + g_zz_xxz_1[i] * pa_x[i];

        g_xzz_xyy_0[i] = g_zz_yy_1[i] * fe_0 + g_zz_xyy_1[i] * pa_x[i];

        g_xzz_xyz_0[i] = g_zz_yz_1[i] * fe_0 + g_zz_xyz_1[i] * pa_x[i];

        g_xzz_xzz_0[i] = g_zz_zz_1[i] * fe_0 + g_zz_xzz_1[i] * pa_x[i];

        g_xzz_yyy_0[i] = g_zz_yyy_1[i] * pa_x[i];

        g_xzz_yyz_0[i] = g_zz_yyz_1[i] * pa_x[i];

        g_xzz_yzz_0[i] = g_zz_yzz_1[i] * pa_x[i];

        g_xzz_zzz_0[i] = g_zz_zzz_1[i] * pa_x[i];
    }

    // Set up 60-70 components of targeted buffer : FF

    auto g_yyy_xxx_0 = pbuffer.data(idx_eri_0_ff + 60);

    auto g_yyy_xxy_0 = pbuffer.data(idx_eri_0_ff + 61);

    auto g_yyy_xxz_0 = pbuffer.data(idx_eri_0_ff + 62);

    auto g_yyy_xyy_0 = pbuffer.data(idx_eri_0_ff + 63);

    auto g_yyy_xyz_0 = pbuffer.data(idx_eri_0_ff + 64);

    auto g_yyy_xzz_0 = pbuffer.data(idx_eri_0_ff + 65);

    auto g_yyy_yyy_0 = pbuffer.data(idx_eri_0_ff + 66);

    auto g_yyy_yyz_0 = pbuffer.data(idx_eri_0_ff + 67);

    auto g_yyy_yzz_0 = pbuffer.data(idx_eri_0_ff + 68);

    auto g_yyy_zzz_0 = pbuffer.data(idx_eri_0_ff + 69);

    #pragma omp simd aligned(g_y_xxx_0, g_y_xxx_1, g_y_xxy_0, g_y_xxy_1, g_y_xxz_0, g_y_xxz_1, g_y_xyy_0, g_y_xyy_1, g_y_xyz_0, g_y_xyz_1, g_y_xzz_0, g_y_xzz_1, g_y_yyy_0, g_y_yyy_1, g_y_yyz_0, g_y_yyz_1, g_y_yzz_0, g_y_yzz_1, g_y_zzz_0, g_y_zzz_1, g_yy_xx_1, g_yy_xxx_1, g_yy_xxy_1, g_yy_xxz_1, g_yy_xy_1, g_yy_xyy_1, g_yy_xyz_1, g_yy_xz_1, g_yy_xzz_1, g_yy_yy_1, g_yy_yyy_1, g_yy_yyz_1, g_yy_yz_1, g_yy_yzz_1, g_yy_zz_1, g_yy_zzz_1, g_yyy_xxx_0, g_yyy_xxy_0, g_yyy_xxz_0, g_yyy_xyy_0, g_yyy_xyz_0, g_yyy_xzz_0, g_yyy_yyy_0, g_yyy_yyz_0, g_yyy_yzz_0, g_yyy_zzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyy_xxx_0[i] = 2.0 * g_y_xxx_0[i] * fbe_0 - 2.0 * g_y_xxx_1[i] * fz_be_0 + g_yy_xxx_1[i] * pa_y[i];

        g_yyy_xxy_0[i] = 2.0 * g_y_xxy_0[i] * fbe_0 - 2.0 * g_y_xxy_1[i] * fz_be_0 + g_yy_xx_1[i] * fe_0 + g_yy_xxy_1[i] * pa_y[i];

        g_yyy_xxz_0[i] = 2.0 * g_y_xxz_0[i] * fbe_0 - 2.0 * g_y_xxz_1[i] * fz_be_0 + g_yy_xxz_1[i] * pa_y[i];

        g_yyy_xyy_0[i] = 2.0 * g_y_xyy_0[i] * fbe_0 - 2.0 * g_y_xyy_1[i] * fz_be_0 + 2.0 * g_yy_xy_1[i] * fe_0 + g_yy_xyy_1[i] * pa_y[i];

        g_yyy_xyz_0[i] = 2.0 * g_y_xyz_0[i] * fbe_0 - 2.0 * g_y_xyz_1[i] * fz_be_0 + g_yy_xz_1[i] * fe_0 + g_yy_xyz_1[i] * pa_y[i];

        g_yyy_xzz_0[i] = 2.0 * g_y_xzz_0[i] * fbe_0 - 2.0 * g_y_xzz_1[i] * fz_be_0 + g_yy_xzz_1[i] * pa_y[i];

        g_yyy_yyy_0[i] = 2.0 * g_y_yyy_0[i] * fbe_0 - 2.0 * g_y_yyy_1[i] * fz_be_0 + 3.0 * g_yy_yy_1[i] * fe_0 + g_yy_yyy_1[i] * pa_y[i];

        g_yyy_yyz_0[i] = 2.0 * g_y_yyz_0[i] * fbe_0 - 2.0 * g_y_yyz_1[i] * fz_be_0 + 2.0 * g_yy_yz_1[i] * fe_0 + g_yy_yyz_1[i] * pa_y[i];

        g_yyy_yzz_0[i] = 2.0 * g_y_yzz_0[i] * fbe_0 - 2.0 * g_y_yzz_1[i] * fz_be_0 + g_yy_zz_1[i] * fe_0 + g_yy_yzz_1[i] * pa_y[i];

        g_yyy_zzz_0[i] = 2.0 * g_y_zzz_0[i] * fbe_0 - 2.0 * g_y_zzz_1[i] * fz_be_0 + g_yy_zzz_1[i] * pa_y[i];
    }

    // Set up 70-80 components of targeted buffer : FF

    auto g_yyz_xxx_0 = pbuffer.data(idx_eri_0_ff + 70);

    auto g_yyz_xxy_0 = pbuffer.data(idx_eri_0_ff + 71);

    auto g_yyz_xxz_0 = pbuffer.data(idx_eri_0_ff + 72);

    auto g_yyz_xyy_0 = pbuffer.data(idx_eri_0_ff + 73);

    auto g_yyz_xyz_0 = pbuffer.data(idx_eri_0_ff + 74);

    auto g_yyz_xzz_0 = pbuffer.data(idx_eri_0_ff + 75);

    auto g_yyz_yyy_0 = pbuffer.data(idx_eri_0_ff + 76);

    auto g_yyz_yyz_0 = pbuffer.data(idx_eri_0_ff + 77);

    auto g_yyz_yzz_0 = pbuffer.data(idx_eri_0_ff + 78);

    auto g_yyz_zzz_0 = pbuffer.data(idx_eri_0_ff + 79);

    #pragma omp simd aligned(g_yy_xx_1, g_yy_xxx_1, g_yy_xxy_1, g_yy_xxz_1, g_yy_xy_1, g_yy_xyy_1, g_yy_xyz_1, g_yy_xz_1, g_yy_xzz_1, g_yy_yy_1, g_yy_yyy_1, g_yy_yyz_1, g_yy_yz_1, g_yy_yzz_1, g_yy_zz_1, g_yy_zzz_1, g_yyz_xxx_0, g_yyz_xxy_0, g_yyz_xxz_0, g_yyz_xyy_0, g_yyz_xyz_0, g_yyz_xzz_0, g_yyz_yyy_0, g_yyz_yyz_0, g_yyz_yzz_0, g_yyz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyz_xxx_0[i] = g_yy_xxx_1[i] * pa_z[i];

        g_yyz_xxy_0[i] = g_yy_xxy_1[i] * pa_z[i];

        g_yyz_xxz_0[i] = g_yy_xx_1[i] * fe_0 + g_yy_xxz_1[i] * pa_z[i];

        g_yyz_xyy_0[i] = g_yy_xyy_1[i] * pa_z[i];

        g_yyz_xyz_0[i] = g_yy_xy_1[i] * fe_0 + g_yy_xyz_1[i] * pa_z[i];

        g_yyz_xzz_0[i] = 2.0 * g_yy_xz_1[i] * fe_0 + g_yy_xzz_1[i] * pa_z[i];

        g_yyz_yyy_0[i] = g_yy_yyy_1[i] * pa_z[i];

        g_yyz_yyz_0[i] = g_yy_yy_1[i] * fe_0 + g_yy_yyz_1[i] * pa_z[i];

        g_yyz_yzz_0[i] = 2.0 * g_yy_yz_1[i] * fe_0 + g_yy_yzz_1[i] * pa_z[i];

        g_yyz_zzz_0[i] = 3.0 * g_yy_zz_1[i] * fe_0 + g_yy_zzz_1[i] * pa_z[i];
    }

    // Set up 80-90 components of targeted buffer : FF

    auto g_yzz_xxx_0 = pbuffer.data(idx_eri_0_ff + 80);

    auto g_yzz_xxy_0 = pbuffer.data(idx_eri_0_ff + 81);

    auto g_yzz_xxz_0 = pbuffer.data(idx_eri_0_ff + 82);

    auto g_yzz_xyy_0 = pbuffer.data(idx_eri_0_ff + 83);

    auto g_yzz_xyz_0 = pbuffer.data(idx_eri_0_ff + 84);

    auto g_yzz_xzz_0 = pbuffer.data(idx_eri_0_ff + 85);

    auto g_yzz_yyy_0 = pbuffer.data(idx_eri_0_ff + 86);

    auto g_yzz_yyz_0 = pbuffer.data(idx_eri_0_ff + 87);

    auto g_yzz_yzz_0 = pbuffer.data(idx_eri_0_ff + 88);

    auto g_yzz_zzz_0 = pbuffer.data(idx_eri_0_ff + 89);

    #pragma omp simd aligned(g_yzz_xxx_0, g_yzz_xxy_0, g_yzz_xxz_0, g_yzz_xyy_0, g_yzz_xyz_0, g_yzz_xzz_0, g_yzz_yyy_0, g_yzz_yyz_0, g_yzz_yzz_0, g_yzz_zzz_0, g_zz_xx_1, g_zz_xxx_1, g_zz_xxy_1, g_zz_xxz_1, g_zz_xy_1, g_zz_xyy_1, g_zz_xyz_1, g_zz_xz_1, g_zz_xzz_1, g_zz_yy_1, g_zz_yyy_1, g_zz_yyz_1, g_zz_yz_1, g_zz_yzz_1, g_zz_zz_1, g_zz_zzz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzz_xxx_0[i] = g_zz_xxx_1[i] * pa_y[i];

        g_yzz_xxy_0[i] = g_zz_xx_1[i] * fe_0 + g_zz_xxy_1[i] * pa_y[i];

        g_yzz_xxz_0[i] = g_zz_xxz_1[i] * pa_y[i];

        g_yzz_xyy_0[i] = 2.0 * g_zz_xy_1[i] * fe_0 + g_zz_xyy_1[i] * pa_y[i];

        g_yzz_xyz_0[i] = g_zz_xz_1[i] * fe_0 + g_zz_xyz_1[i] * pa_y[i];

        g_yzz_xzz_0[i] = g_zz_xzz_1[i] * pa_y[i];

        g_yzz_yyy_0[i] = 3.0 * g_zz_yy_1[i] * fe_0 + g_zz_yyy_1[i] * pa_y[i];

        g_yzz_yyz_0[i] = 2.0 * g_zz_yz_1[i] * fe_0 + g_zz_yyz_1[i] * pa_y[i];

        g_yzz_yzz_0[i] = g_zz_zz_1[i] * fe_0 + g_zz_yzz_1[i] * pa_y[i];

        g_yzz_zzz_0[i] = g_zz_zzz_1[i] * pa_y[i];
    }

    // Set up 90-100 components of targeted buffer : FF

    auto g_zzz_xxx_0 = pbuffer.data(idx_eri_0_ff + 90);

    auto g_zzz_xxy_0 = pbuffer.data(idx_eri_0_ff + 91);

    auto g_zzz_xxz_0 = pbuffer.data(idx_eri_0_ff + 92);

    auto g_zzz_xyy_0 = pbuffer.data(idx_eri_0_ff + 93);

    auto g_zzz_xyz_0 = pbuffer.data(idx_eri_0_ff + 94);

    auto g_zzz_xzz_0 = pbuffer.data(idx_eri_0_ff + 95);

    auto g_zzz_yyy_0 = pbuffer.data(idx_eri_0_ff + 96);

    auto g_zzz_yyz_0 = pbuffer.data(idx_eri_0_ff + 97);

    auto g_zzz_yzz_0 = pbuffer.data(idx_eri_0_ff + 98);

    auto g_zzz_zzz_0 = pbuffer.data(idx_eri_0_ff + 99);

    #pragma omp simd aligned(g_z_xxx_0, g_z_xxx_1, g_z_xxy_0, g_z_xxy_1, g_z_xxz_0, g_z_xxz_1, g_z_xyy_0, g_z_xyy_1, g_z_xyz_0, g_z_xyz_1, g_z_xzz_0, g_z_xzz_1, g_z_yyy_0, g_z_yyy_1, g_z_yyz_0, g_z_yyz_1, g_z_yzz_0, g_z_yzz_1, g_z_zzz_0, g_z_zzz_1, g_zz_xx_1, g_zz_xxx_1, g_zz_xxy_1, g_zz_xxz_1, g_zz_xy_1, g_zz_xyy_1, g_zz_xyz_1, g_zz_xz_1, g_zz_xzz_1, g_zz_yy_1, g_zz_yyy_1, g_zz_yyz_1, g_zz_yz_1, g_zz_yzz_1, g_zz_zz_1, g_zz_zzz_1, g_zzz_xxx_0, g_zzz_xxy_0, g_zzz_xxz_0, g_zzz_xyy_0, g_zzz_xyz_0, g_zzz_xzz_0, g_zzz_yyy_0, g_zzz_yyz_0, g_zzz_yzz_0, g_zzz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzz_xxx_0[i] = 2.0 * g_z_xxx_0[i] * fbe_0 - 2.0 * g_z_xxx_1[i] * fz_be_0 + g_zz_xxx_1[i] * pa_z[i];

        g_zzz_xxy_0[i] = 2.0 * g_z_xxy_0[i] * fbe_0 - 2.0 * g_z_xxy_1[i] * fz_be_0 + g_zz_xxy_1[i] * pa_z[i];

        g_zzz_xxz_0[i] = 2.0 * g_z_xxz_0[i] * fbe_0 - 2.0 * g_z_xxz_1[i] * fz_be_0 + g_zz_xx_1[i] * fe_0 + g_zz_xxz_1[i] * pa_z[i];

        g_zzz_xyy_0[i] = 2.0 * g_z_xyy_0[i] * fbe_0 - 2.0 * g_z_xyy_1[i] * fz_be_0 + g_zz_xyy_1[i] * pa_z[i];

        g_zzz_xyz_0[i] = 2.0 * g_z_xyz_0[i] * fbe_0 - 2.0 * g_z_xyz_1[i] * fz_be_0 + g_zz_xy_1[i] * fe_0 + g_zz_xyz_1[i] * pa_z[i];

        g_zzz_xzz_0[i] = 2.0 * g_z_xzz_0[i] * fbe_0 - 2.0 * g_z_xzz_1[i] * fz_be_0 + 2.0 * g_zz_xz_1[i] * fe_0 + g_zz_xzz_1[i] * pa_z[i];

        g_zzz_yyy_0[i] = 2.0 * g_z_yyy_0[i] * fbe_0 - 2.0 * g_z_yyy_1[i] * fz_be_0 + g_zz_yyy_1[i] * pa_z[i];

        g_zzz_yyz_0[i] = 2.0 * g_z_yyz_0[i] * fbe_0 - 2.0 * g_z_yyz_1[i] * fz_be_0 + g_zz_yy_1[i] * fe_0 + g_zz_yyz_1[i] * pa_z[i];

        g_zzz_yzz_0[i] = 2.0 * g_z_yzz_0[i] * fbe_0 - 2.0 * g_z_yzz_1[i] * fz_be_0 + 2.0 * g_zz_yz_1[i] * fe_0 + g_zz_yzz_1[i] * pa_z[i];

        g_zzz_zzz_0[i] = 2.0 * g_z_zzz_0[i] * fbe_0 - 2.0 * g_z_zzz_1[i] * fz_be_0 + 3.0 * g_zz_zz_1[i] * fe_0 + g_zz_zzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

