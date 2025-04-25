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

#include "TwoCenterElectronRepulsionPrimRecDF.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_df(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_df,
                                const size_t idx_eri_0_sf,
                                const size_t idx_eri_1_sf,
                                const size_t idx_eri_1_pd,
                                const size_t idx_eri_1_pf,
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

    // Set up components of auxiliary buffer : SF

    auto g_0_xxx_0 = pbuffer.data(idx_eri_0_sf);

    auto g_0_xxy_0 = pbuffer.data(idx_eri_0_sf + 1);

    auto g_0_xxz_0 = pbuffer.data(idx_eri_0_sf + 2);

    auto g_0_xyy_0 = pbuffer.data(idx_eri_0_sf + 3);

    auto g_0_xyz_0 = pbuffer.data(idx_eri_0_sf + 4);

    auto g_0_xzz_0 = pbuffer.data(idx_eri_0_sf + 5);

    auto g_0_yyy_0 = pbuffer.data(idx_eri_0_sf + 6);

    auto g_0_yyz_0 = pbuffer.data(idx_eri_0_sf + 7);

    auto g_0_yzz_0 = pbuffer.data(idx_eri_0_sf + 8);

    auto g_0_zzz_0 = pbuffer.data(idx_eri_0_sf + 9);

    // Set up components of auxiliary buffer : SF

    auto g_0_xxx_1 = pbuffer.data(idx_eri_1_sf);

    auto g_0_xxy_1 = pbuffer.data(idx_eri_1_sf + 1);

    auto g_0_xxz_1 = pbuffer.data(idx_eri_1_sf + 2);

    auto g_0_xyy_1 = pbuffer.data(idx_eri_1_sf + 3);

    auto g_0_xyz_1 = pbuffer.data(idx_eri_1_sf + 4);

    auto g_0_xzz_1 = pbuffer.data(idx_eri_1_sf + 5);

    auto g_0_yyy_1 = pbuffer.data(idx_eri_1_sf + 6);

    auto g_0_yyz_1 = pbuffer.data(idx_eri_1_sf + 7);

    auto g_0_yzz_1 = pbuffer.data(idx_eri_1_sf + 8);

    auto g_0_zzz_1 = pbuffer.data(idx_eri_1_sf + 9);

    // Set up components of auxiliary buffer : PD

    auto g_x_xx_1 = pbuffer.data(idx_eri_1_pd);

    auto g_x_xy_1 = pbuffer.data(idx_eri_1_pd + 1);

    auto g_x_xz_1 = pbuffer.data(idx_eri_1_pd + 2);

    auto g_x_yy_1 = pbuffer.data(idx_eri_1_pd + 3);

    auto g_x_yz_1 = pbuffer.data(idx_eri_1_pd + 4);

    auto g_x_zz_1 = pbuffer.data(idx_eri_1_pd + 5);

    auto g_y_xx_1 = pbuffer.data(idx_eri_1_pd + 6);

    auto g_y_xy_1 = pbuffer.data(idx_eri_1_pd + 7);

    auto g_y_xz_1 = pbuffer.data(idx_eri_1_pd + 8);

    auto g_y_yy_1 = pbuffer.data(idx_eri_1_pd + 9);

    auto g_y_yz_1 = pbuffer.data(idx_eri_1_pd + 10);

    auto g_y_zz_1 = pbuffer.data(idx_eri_1_pd + 11);

    auto g_z_xx_1 = pbuffer.data(idx_eri_1_pd + 12);

    auto g_z_xy_1 = pbuffer.data(idx_eri_1_pd + 13);

    auto g_z_xz_1 = pbuffer.data(idx_eri_1_pd + 14);

    auto g_z_yy_1 = pbuffer.data(idx_eri_1_pd + 15);

    auto g_z_yz_1 = pbuffer.data(idx_eri_1_pd + 16);

    auto g_z_zz_1 = pbuffer.data(idx_eri_1_pd + 17);

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

    // Set up 0-10 components of targeted buffer : DF

    auto g_xx_xxx_0 = pbuffer.data(idx_eri_0_df);

    auto g_xx_xxy_0 = pbuffer.data(idx_eri_0_df + 1);

    auto g_xx_xxz_0 = pbuffer.data(idx_eri_0_df + 2);

    auto g_xx_xyy_0 = pbuffer.data(idx_eri_0_df + 3);

    auto g_xx_xyz_0 = pbuffer.data(idx_eri_0_df + 4);

    auto g_xx_xzz_0 = pbuffer.data(idx_eri_0_df + 5);

    auto g_xx_yyy_0 = pbuffer.data(idx_eri_0_df + 6);

    auto g_xx_yyz_0 = pbuffer.data(idx_eri_0_df + 7);

    auto g_xx_yzz_0 = pbuffer.data(idx_eri_0_df + 8);

    auto g_xx_zzz_0 = pbuffer.data(idx_eri_0_df + 9);

    #pragma omp simd aligned(g_0_xxx_0, g_0_xxx_1, g_0_xxy_0, g_0_xxy_1, g_0_xxz_0, g_0_xxz_1, g_0_xyy_0, g_0_xyy_1, g_0_xyz_0, g_0_xyz_1, g_0_xzz_0, g_0_xzz_1, g_0_yyy_0, g_0_yyy_1, g_0_yyz_0, g_0_yyz_1, g_0_yzz_0, g_0_yzz_1, g_0_zzz_0, g_0_zzz_1, g_x_xx_1, g_x_xxx_1, g_x_xxy_1, g_x_xxz_1, g_x_xy_1, g_x_xyy_1, g_x_xyz_1, g_x_xz_1, g_x_xzz_1, g_x_yy_1, g_x_yyy_1, g_x_yyz_1, g_x_yz_1, g_x_yzz_1, g_x_zz_1, g_x_zzz_1, g_xx_xxx_0, g_xx_xxy_0, g_xx_xxz_0, g_xx_xyy_0, g_xx_xyz_0, g_xx_xzz_0, g_xx_yyy_0, g_xx_yyz_0, g_xx_yzz_0, g_xx_zzz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xx_xxx_0[i] = g_0_xxx_0[i] * fbe_0 - g_0_xxx_1[i] * fz_be_0 + 3.0 * g_x_xx_1[i] * fe_0 + g_x_xxx_1[i] * pa_x[i];

        g_xx_xxy_0[i] = g_0_xxy_0[i] * fbe_0 - g_0_xxy_1[i] * fz_be_0 + 2.0 * g_x_xy_1[i] * fe_0 + g_x_xxy_1[i] * pa_x[i];

        g_xx_xxz_0[i] = g_0_xxz_0[i] * fbe_0 - g_0_xxz_1[i] * fz_be_0 + 2.0 * g_x_xz_1[i] * fe_0 + g_x_xxz_1[i] * pa_x[i];

        g_xx_xyy_0[i] = g_0_xyy_0[i] * fbe_0 - g_0_xyy_1[i] * fz_be_0 + g_x_yy_1[i] * fe_0 + g_x_xyy_1[i] * pa_x[i];

        g_xx_xyz_0[i] = g_0_xyz_0[i] * fbe_0 - g_0_xyz_1[i] * fz_be_0 + g_x_yz_1[i] * fe_0 + g_x_xyz_1[i] * pa_x[i];

        g_xx_xzz_0[i] = g_0_xzz_0[i] * fbe_0 - g_0_xzz_1[i] * fz_be_0 + g_x_zz_1[i] * fe_0 + g_x_xzz_1[i] * pa_x[i];

        g_xx_yyy_0[i] = g_0_yyy_0[i] * fbe_0 - g_0_yyy_1[i] * fz_be_0 + g_x_yyy_1[i] * pa_x[i];

        g_xx_yyz_0[i] = g_0_yyz_0[i] * fbe_0 - g_0_yyz_1[i] * fz_be_0 + g_x_yyz_1[i] * pa_x[i];

        g_xx_yzz_0[i] = g_0_yzz_0[i] * fbe_0 - g_0_yzz_1[i] * fz_be_0 + g_x_yzz_1[i] * pa_x[i];

        g_xx_zzz_0[i] = g_0_zzz_0[i] * fbe_0 - g_0_zzz_1[i] * fz_be_0 + g_x_zzz_1[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : DF

    auto g_xy_xxx_0 = pbuffer.data(idx_eri_0_df + 10);

    auto g_xy_xxy_0 = pbuffer.data(idx_eri_0_df + 11);

    auto g_xy_xxz_0 = pbuffer.data(idx_eri_0_df + 12);

    auto g_xy_xyy_0 = pbuffer.data(idx_eri_0_df + 13);

    auto g_xy_xyz_0 = pbuffer.data(idx_eri_0_df + 14);

    auto g_xy_xzz_0 = pbuffer.data(idx_eri_0_df + 15);

    auto g_xy_yyy_0 = pbuffer.data(idx_eri_0_df + 16);

    auto g_xy_yyz_0 = pbuffer.data(idx_eri_0_df + 17);

    auto g_xy_yzz_0 = pbuffer.data(idx_eri_0_df + 18);

    auto g_xy_zzz_0 = pbuffer.data(idx_eri_0_df + 19);

    #pragma omp simd aligned(g_x_xxx_1, g_x_xxz_1, g_x_xzz_1, g_xy_xxx_0, g_xy_xxy_0, g_xy_xxz_0, g_xy_xyy_0, g_xy_xyz_0, g_xy_xzz_0, g_xy_yyy_0, g_xy_yyz_0, g_xy_yzz_0, g_xy_zzz_0, g_y_xxy_1, g_y_xy_1, g_y_xyy_1, g_y_xyz_1, g_y_yy_1, g_y_yyy_1, g_y_yyz_1, g_y_yz_1, g_y_yzz_1, g_y_zzz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xy_xxx_0[i] = g_x_xxx_1[i] * pa_y[i];

        g_xy_xxy_0[i] = 2.0 * g_y_xy_1[i] * fe_0 + g_y_xxy_1[i] * pa_x[i];

        g_xy_xxz_0[i] = g_x_xxz_1[i] * pa_y[i];

        g_xy_xyy_0[i] = g_y_yy_1[i] * fe_0 + g_y_xyy_1[i] * pa_x[i];

        g_xy_xyz_0[i] = g_y_yz_1[i] * fe_0 + g_y_xyz_1[i] * pa_x[i];

        g_xy_xzz_0[i] = g_x_xzz_1[i] * pa_y[i];

        g_xy_yyy_0[i] = g_y_yyy_1[i] * pa_x[i];

        g_xy_yyz_0[i] = g_y_yyz_1[i] * pa_x[i];

        g_xy_yzz_0[i] = g_y_yzz_1[i] * pa_x[i];

        g_xy_zzz_0[i] = g_y_zzz_1[i] * pa_x[i];
    }

    // Set up 20-30 components of targeted buffer : DF

    auto g_xz_xxx_0 = pbuffer.data(idx_eri_0_df + 20);

    auto g_xz_xxy_0 = pbuffer.data(idx_eri_0_df + 21);

    auto g_xz_xxz_0 = pbuffer.data(idx_eri_0_df + 22);

    auto g_xz_xyy_0 = pbuffer.data(idx_eri_0_df + 23);

    auto g_xz_xyz_0 = pbuffer.data(idx_eri_0_df + 24);

    auto g_xz_xzz_0 = pbuffer.data(idx_eri_0_df + 25);

    auto g_xz_yyy_0 = pbuffer.data(idx_eri_0_df + 26);

    auto g_xz_yyz_0 = pbuffer.data(idx_eri_0_df + 27);

    auto g_xz_yzz_0 = pbuffer.data(idx_eri_0_df + 28);

    auto g_xz_zzz_0 = pbuffer.data(idx_eri_0_df + 29);

    #pragma omp simd aligned(g_x_xxx_1, g_x_xxy_1, g_x_xyy_1, g_xz_xxx_0, g_xz_xxy_0, g_xz_xxz_0, g_xz_xyy_0, g_xz_xyz_0, g_xz_xzz_0, g_xz_yyy_0, g_xz_yyz_0, g_xz_yzz_0, g_xz_zzz_0, g_z_xxz_1, g_z_xyz_1, g_z_xz_1, g_z_xzz_1, g_z_yyy_1, g_z_yyz_1, g_z_yz_1, g_z_yzz_1, g_z_zz_1, g_z_zzz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xz_xxx_0[i] = g_x_xxx_1[i] * pa_z[i];

        g_xz_xxy_0[i] = g_x_xxy_1[i] * pa_z[i];

        g_xz_xxz_0[i] = 2.0 * g_z_xz_1[i] * fe_0 + g_z_xxz_1[i] * pa_x[i];

        g_xz_xyy_0[i] = g_x_xyy_1[i] * pa_z[i];

        g_xz_xyz_0[i] = g_z_yz_1[i] * fe_0 + g_z_xyz_1[i] * pa_x[i];

        g_xz_xzz_0[i] = g_z_zz_1[i] * fe_0 + g_z_xzz_1[i] * pa_x[i];

        g_xz_yyy_0[i] = g_z_yyy_1[i] * pa_x[i];

        g_xz_yyz_0[i] = g_z_yyz_1[i] * pa_x[i];

        g_xz_yzz_0[i] = g_z_yzz_1[i] * pa_x[i];

        g_xz_zzz_0[i] = g_z_zzz_1[i] * pa_x[i];
    }

    // Set up 30-40 components of targeted buffer : DF

    auto g_yy_xxx_0 = pbuffer.data(idx_eri_0_df + 30);

    auto g_yy_xxy_0 = pbuffer.data(idx_eri_0_df + 31);

    auto g_yy_xxz_0 = pbuffer.data(idx_eri_0_df + 32);

    auto g_yy_xyy_0 = pbuffer.data(idx_eri_0_df + 33);

    auto g_yy_xyz_0 = pbuffer.data(idx_eri_0_df + 34);

    auto g_yy_xzz_0 = pbuffer.data(idx_eri_0_df + 35);

    auto g_yy_yyy_0 = pbuffer.data(idx_eri_0_df + 36);

    auto g_yy_yyz_0 = pbuffer.data(idx_eri_0_df + 37);

    auto g_yy_yzz_0 = pbuffer.data(idx_eri_0_df + 38);

    auto g_yy_zzz_0 = pbuffer.data(idx_eri_0_df + 39);

    #pragma omp simd aligned(g_0_xxx_0, g_0_xxx_1, g_0_xxy_0, g_0_xxy_1, g_0_xxz_0, g_0_xxz_1, g_0_xyy_0, g_0_xyy_1, g_0_xyz_0, g_0_xyz_1, g_0_xzz_0, g_0_xzz_1, g_0_yyy_0, g_0_yyy_1, g_0_yyz_0, g_0_yyz_1, g_0_yzz_0, g_0_yzz_1, g_0_zzz_0, g_0_zzz_1, g_y_xx_1, g_y_xxx_1, g_y_xxy_1, g_y_xxz_1, g_y_xy_1, g_y_xyy_1, g_y_xyz_1, g_y_xz_1, g_y_xzz_1, g_y_yy_1, g_y_yyy_1, g_y_yyz_1, g_y_yz_1, g_y_yzz_1, g_y_zz_1, g_y_zzz_1, g_yy_xxx_0, g_yy_xxy_0, g_yy_xxz_0, g_yy_xyy_0, g_yy_xyz_0, g_yy_xzz_0, g_yy_yyy_0, g_yy_yyz_0, g_yy_yzz_0, g_yy_zzz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yy_xxx_0[i] = g_0_xxx_0[i] * fbe_0 - g_0_xxx_1[i] * fz_be_0 + g_y_xxx_1[i] * pa_y[i];

        g_yy_xxy_0[i] = g_0_xxy_0[i] * fbe_0 - g_0_xxy_1[i] * fz_be_0 + g_y_xx_1[i] * fe_0 + g_y_xxy_1[i] * pa_y[i];

        g_yy_xxz_0[i] = g_0_xxz_0[i] * fbe_0 - g_0_xxz_1[i] * fz_be_0 + g_y_xxz_1[i] * pa_y[i];

        g_yy_xyy_0[i] = g_0_xyy_0[i] * fbe_0 - g_0_xyy_1[i] * fz_be_0 + 2.0 * g_y_xy_1[i] * fe_0 + g_y_xyy_1[i] * pa_y[i];

        g_yy_xyz_0[i] = g_0_xyz_0[i] * fbe_0 - g_0_xyz_1[i] * fz_be_0 + g_y_xz_1[i] * fe_0 + g_y_xyz_1[i] * pa_y[i];

        g_yy_xzz_0[i] = g_0_xzz_0[i] * fbe_0 - g_0_xzz_1[i] * fz_be_0 + g_y_xzz_1[i] * pa_y[i];

        g_yy_yyy_0[i] = g_0_yyy_0[i] * fbe_0 - g_0_yyy_1[i] * fz_be_0 + 3.0 * g_y_yy_1[i] * fe_0 + g_y_yyy_1[i] * pa_y[i];

        g_yy_yyz_0[i] = g_0_yyz_0[i] * fbe_0 - g_0_yyz_1[i] * fz_be_0 + 2.0 * g_y_yz_1[i] * fe_0 + g_y_yyz_1[i] * pa_y[i];

        g_yy_yzz_0[i] = g_0_yzz_0[i] * fbe_0 - g_0_yzz_1[i] * fz_be_0 + g_y_zz_1[i] * fe_0 + g_y_yzz_1[i] * pa_y[i];

        g_yy_zzz_0[i] = g_0_zzz_0[i] * fbe_0 - g_0_zzz_1[i] * fz_be_0 + g_y_zzz_1[i] * pa_y[i];
    }

    // Set up 40-50 components of targeted buffer : DF

    auto g_yz_xxx_0 = pbuffer.data(idx_eri_0_df + 40);

    auto g_yz_xxy_0 = pbuffer.data(idx_eri_0_df + 41);

    auto g_yz_xxz_0 = pbuffer.data(idx_eri_0_df + 42);

    auto g_yz_xyy_0 = pbuffer.data(idx_eri_0_df + 43);

    auto g_yz_xyz_0 = pbuffer.data(idx_eri_0_df + 44);

    auto g_yz_xzz_0 = pbuffer.data(idx_eri_0_df + 45);

    auto g_yz_yyy_0 = pbuffer.data(idx_eri_0_df + 46);

    auto g_yz_yyz_0 = pbuffer.data(idx_eri_0_df + 47);

    auto g_yz_yzz_0 = pbuffer.data(idx_eri_0_df + 48);

    auto g_yz_zzz_0 = pbuffer.data(idx_eri_0_df + 49);

    #pragma omp simd aligned(g_y_xxy_1, g_y_xyy_1, g_y_yyy_1, g_yz_xxx_0, g_yz_xxy_0, g_yz_xxz_0, g_yz_xyy_0, g_yz_xyz_0, g_yz_xzz_0, g_yz_yyy_0, g_yz_yyz_0, g_yz_yzz_0, g_yz_zzz_0, g_z_xxx_1, g_z_xxz_1, g_z_xyz_1, g_z_xz_1, g_z_xzz_1, g_z_yyz_1, g_z_yz_1, g_z_yzz_1, g_z_zz_1, g_z_zzz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yz_xxx_0[i] = g_z_xxx_1[i] * pa_y[i];

        g_yz_xxy_0[i] = g_y_xxy_1[i] * pa_z[i];

        g_yz_xxz_0[i] = g_z_xxz_1[i] * pa_y[i];

        g_yz_xyy_0[i] = g_y_xyy_1[i] * pa_z[i];

        g_yz_xyz_0[i] = g_z_xz_1[i] * fe_0 + g_z_xyz_1[i] * pa_y[i];

        g_yz_xzz_0[i] = g_z_xzz_1[i] * pa_y[i];

        g_yz_yyy_0[i] = g_y_yyy_1[i] * pa_z[i];

        g_yz_yyz_0[i] = 2.0 * g_z_yz_1[i] * fe_0 + g_z_yyz_1[i] * pa_y[i];

        g_yz_yzz_0[i] = g_z_zz_1[i] * fe_0 + g_z_yzz_1[i] * pa_y[i];

        g_yz_zzz_0[i] = g_z_zzz_1[i] * pa_y[i];
    }

    // Set up 50-60 components of targeted buffer : DF

    auto g_zz_xxx_0 = pbuffer.data(idx_eri_0_df + 50);

    auto g_zz_xxy_0 = pbuffer.data(idx_eri_0_df + 51);

    auto g_zz_xxz_0 = pbuffer.data(idx_eri_0_df + 52);

    auto g_zz_xyy_0 = pbuffer.data(idx_eri_0_df + 53);

    auto g_zz_xyz_0 = pbuffer.data(idx_eri_0_df + 54);

    auto g_zz_xzz_0 = pbuffer.data(idx_eri_0_df + 55);

    auto g_zz_yyy_0 = pbuffer.data(idx_eri_0_df + 56);

    auto g_zz_yyz_0 = pbuffer.data(idx_eri_0_df + 57);

    auto g_zz_yzz_0 = pbuffer.data(idx_eri_0_df + 58);

    auto g_zz_zzz_0 = pbuffer.data(idx_eri_0_df + 59);

    #pragma omp simd aligned(g_0_xxx_0, g_0_xxx_1, g_0_xxy_0, g_0_xxy_1, g_0_xxz_0, g_0_xxz_1, g_0_xyy_0, g_0_xyy_1, g_0_xyz_0, g_0_xyz_1, g_0_xzz_0, g_0_xzz_1, g_0_yyy_0, g_0_yyy_1, g_0_yyz_0, g_0_yyz_1, g_0_yzz_0, g_0_yzz_1, g_0_zzz_0, g_0_zzz_1, g_z_xx_1, g_z_xxx_1, g_z_xxy_1, g_z_xxz_1, g_z_xy_1, g_z_xyy_1, g_z_xyz_1, g_z_xz_1, g_z_xzz_1, g_z_yy_1, g_z_yyy_1, g_z_yyz_1, g_z_yz_1, g_z_yzz_1, g_z_zz_1, g_z_zzz_1, g_zz_xxx_0, g_zz_xxy_0, g_zz_xxz_0, g_zz_xyy_0, g_zz_xyz_0, g_zz_xzz_0, g_zz_yyy_0, g_zz_yyz_0, g_zz_yzz_0, g_zz_zzz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zz_xxx_0[i] = g_0_xxx_0[i] * fbe_0 - g_0_xxx_1[i] * fz_be_0 + g_z_xxx_1[i] * pa_z[i];

        g_zz_xxy_0[i] = g_0_xxy_0[i] * fbe_0 - g_0_xxy_1[i] * fz_be_0 + g_z_xxy_1[i] * pa_z[i];

        g_zz_xxz_0[i] = g_0_xxz_0[i] * fbe_0 - g_0_xxz_1[i] * fz_be_0 + g_z_xx_1[i] * fe_0 + g_z_xxz_1[i] * pa_z[i];

        g_zz_xyy_0[i] = g_0_xyy_0[i] * fbe_0 - g_0_xyy_1[i] * fz_be_0 + g_z_xyy_1[i] * pa_z[i];

        g_zz_xyz_0[i] = g_0_xyz_0[i] * fbe_0 - g_0_xyz_1[i] * fz_be_0 + g_z_xy_1[i] * fe_0 + g_z_xyz_1[i] * pa_z[i];

        g_zz_xzz_0[i] = g_0_xzz_0[i] * fbe_0 - g_0_xzz_1[i] * fz_be_0 + 2.0 * g_z_xz_1[i] * fe_0 + g_z_xzz_1[i] * pa_z[i];

        g_zz_yyy_0[i] = g_0_yyy_0[i] * fbe_0 - g_0_yyy_1[i] * fz_be_0 + g_z_yyy_1[i] * pa_z[i];

        g_zz_yyz_0[i] = g_0_yyz_0[i] * fbe_0 - g_0_yyz_1[i] * fz_be_0 + g_z_yy_1[i] * fe_0 + g_z_yyz_1[i] * pa_z[i];

        g_zz_yzz_0[i] = g_0_yzz_0[i] * fbe_0 - g_0_yzz_1[i] * fz_be_0 + 2.0 * g_z_yz_1[i] * fe_0 + g_z_yzz_1[i] * pa_z[i];

        g_zz_zzz_0[i] = g_0_zzz_0[i] * fbe_0 - g_0_zzz_1[i] * fz_be_0 + 3.0 * g_z_zz_1[i] * fe_0 + g_z_zzz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

