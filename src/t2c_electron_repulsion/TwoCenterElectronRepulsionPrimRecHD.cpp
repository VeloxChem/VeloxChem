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

#include "TwoCenterElectronRepulsionPrimRecHD.hpp"

namespace t2ceri { // t2ceri namespace

auto
comp_prim_electron_repulsion_hd(CSimdArray<double>& pbuffer, 
                                const size_t idx_eri_0_hd,
                                const size_t idx_eri_0_fd,
                                const size_t idx_eri_1_fd,
                                const size_t idx_eri_1_gp,
                                const size_t idx_eri_1_gd,
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

    // Set up components of auxiliary buffer : FD

    auto g_xxx_xx_0 = pbuffer.data(idx_eri_0_fd);

    auto g_xxx_xy_0 = pbuffer.data(idx_eri_0_fd + 1);

    auto g_xxx_xz_0 = pbuffer.data(idx_eri_0_fd + 2);

    auto g_xxx_yy_0 = pbuffer.data(idx_eri_0_fd + 3);

    auto g_xxx_yz_0 = pbuffer.data(idx_eri_0_fd + 4);

    auto g_xxx_zz_0 = pbuffer.data(idx_eri_0_fd + 5);

    auto g_xxy_xx_0 = pbuffer.data(idx_eri_0_fd + 6);

    auto g_xxy_xz_0 = pbuffer.data(idx_eri_0_fd + 8);

    auto g_xxz_xx_0 = pbuffer.data(idx_eri_0_fd + 12);

    auto g_xxz_xy_0 = pbuffer.data(idx_eri_0_fd + 13);

    auto g_xyy_xy_0 = pbuffer.data(idx_eri_0_fd + 19);

    auto g_xyy_yy_0 = pbuffer.data(idx_eri_0_fd + 21);

    auto g_xyy_yz_0 = pbuffer.data(idx_eri_0_fd + 22);

    auto g_xyy_zz_0 = pbuffer.data(idx_eri_0_fd + 23);

    auto g_xzz_xz_0 = pbuffer.data(idx_eri_0_fd + 32);

    auto g_xzz_yy_0 = pbuffer.data(idx_eri_0_fd + 33);

    auto g_xzz_yz_0 = pbuffer.data(idx_eri_0_fd + 34);

    auto g_xzz_zz_0 = pbuffer.data(idx_eri_0_fd + 35);

    auto g_yyy_xx_0 = pbuffer.data(idx_eri_0_fd + 36);

    auto g_yyy_xy_0 = pbuffer.data(idx_eri_0_fd + 37);

    auto g_yyy_xz_0 = pbuffer.data(idx_eri_0_fd + 38);

    auto g_yyy_yy_0 = pbuffer.data(idx_eri_0_fd + 39);

    auto g_yyy_yz_0 = pbuffer.data(idx_eri_0_fd + 40);

    auto g_yyy_zz_0 = pbuffer.data(idx_eri_0_fd + 41);

    auto g_yyz_xy_0 = pbuffer.data(idx_eri_0_fd + 43);

    auto g_yyz_yy_0 = pbuffer.data(idx_eri_0_fd + 45);

    auto g_yzz_xx_0 = pbuffer.data(idx_eri_0_fd + 48);

    auto g_yzz_xz_0 = pbuffer.data(idx_eri_0_fd + 50);

    auto g_yzz_yz_0 = pbuffer.data(idx_eri_0_fd + 52);

    auto g_yzz_zz_0 = pbuffer.data(idx_eri_0_fd + 53);

    auto g_zzz_xx_0 = pbuffer.data(idx_eri_0_fd + 54);

    auto g_zzz_xy_0 = pbuffer.data(idx_eri_0_fd + 55);

    auto g_zzz_xz_0 = pbuffer.data(idx_eri_0_fd + 56);

    auto g_zzz_yy_0 = pbuffer.data(idx_eri_0_fd + 57);

    auto g_zzz_yz_0 = pbuffer.data(idx_eri_0_fd + 58);

    auto g_zzz_zz_0 = pbuffer.data(idx_eri_0_fd + 59);

    // Set up components of auxiliary buffer : FD

    auto g_xxx_xx_1 = pbuffer.data(idx_eri_1_fd);

    auto g_xxx_xy_1 = pbuffer.data(idx_eri_1_fd + 1);

    auto g_xxx_xz_1 = pbuffer.data(idx_eri_1_fd + 2);

    auto g_xxx_yy_1 = pbuffer.data(idx_eri_1_fd + 3);

    auto g_xxx_yz_1 = pbuffer.data(idx_eri_1_fd + 4);

    auto g_xxx_zz_1 = pbuffer.data(idx_eri_1_fd + 5);

    auto g_xxy_xx_1 = pbuffer.data(idx_eri_1_fd + 6);

    auto g_xxy_xz_1 = pbuffer.data(idx_eri_1_fd + 8);

    auto g_xxz_xx_1 = pbuffer.data(idx_eri_1_fd + 12);

    auto g_xxz_xy_1 = pbuffer.data(idx_eri_1_fd + 13);

    auto g_xyy_xy_1 = pbuffer.data(idx_eri_1_fd + 19);

    auto g_xyy_yy_1 = pbuffer.data(idx_eri_1_fd + 21);

    auto g_xyy_yz_1 = pbuffer.data(idx_eri_1_fd + 22);

    auto g_xyy_zz_1 = pbuffer.data(idx_eri_1_fd + 23);

    auto g_xzz_xz_1 = pbuffer.data(idx_eri_1_fd + 32);

    auto g_xzz_yy_1 = pbuffer.data(idx_eri_1_fd + 33);

    auto g_xzz_yz_1 = pbuffer.data(idx_eri_1_fd + 34);

    auto g_xzz_zz_1 = pbuffer.data(idx_eri_1_fd + 35);

    auto g_yyy_xx_1 = pbuffer.data(idx_eri_1_fd + 36);

    auto g_yyy_xy_1 = pbuffer.data(idx_eri_1_fd + 37);

    auto g_yyy_xz_1 = pbuffer.data(idx_eri_1_fd + 38);

    auto g_yyy_yy_1 = pbuffer.data(idx_eri_1_fd + 39);

    auto g_yyy_yz_1 = pbuffer.data(idx_eri_1_fd + 40);

    auto g_yyy_zz_1 = pbuffer.data(idx_eri_1_fd + 41);

    auto g_yyz_xy_1 = pbuffer.data(idx_eri_1_fd + 43);

    auto g_yyz_yy_1 = pbuffer.data(idx_eri_1_fd + 45);

    auto g_yzz_xx_1 = pbuffer.data(idx_eri_1_fd + 48);

    auto g_yzz_xz_1 = pbuffer.data(idx_eri_1_fd + 50);

    auto g_yzz_yz_1 = pbuffer.data(idx_eri_1_fd + 52);

    auto g_yzz_zz_1 = pbuffer.data(idx_eri_1_fd + 53);

    auto g_zzz_xx_1 = pbuffer.data(idx_eri_1_fd + 54);

    auto g_zzz_xy_1 = pbuffer.data(idx_eri_1_fd + 55);

    auto g_zzz_xz_1 = pbuffer.data(idx_eri_1_fd + 56);

    auto g_zzz_yy_1 = pbuffer.data(idx_eri_1_fd + 57);

    auto g_zzz_yz_1 = pbuffer.data(idx_eri_1_fd + 58);

    auto g_zzz_zz_1 = pbuffer.data(idx_eri_1_fd + 59);

    // Set up components of auxiliary buffer : GP

    auto g_xxxx_x_1 = pbuffer.data(idx_eri_1_gp);

    auto g_xxxx_y_1 = pbuffer.data(idx_eri_1_gp + 1);

    auto g_xxxx_z_1 = pbuffer.data(idx_eri_1_gp + 2);

    auto g_xxxz_z_1 = pbuffer.data(idx_eri_1_gp + 8);

    auto g_xxyy_x_1 = pbuffer.data(idx_eri_1_gp + 9);

    auto g_xxyy_y_1 = pbuffer.data(idx_eri_1_gp + 10);

    auto g_xxyy_z_1 = pbuffer.data(idx_eri_1_gp + 11);

    auto g_xxzz_x_1 = pbuffer.data(idx_eri_1_gp + 15);

    auto g_xxzz_y_1 = pbuffer.data(idx_eri_1_gp + 16);

    auto g_xxzz_z_1 = pbuffer.data(idx_eri_1_gp + 17);

    auto g_xyyy_y_1 = pbuffer.data(idx_eri_1_gp + 19);

    auto g_xzzz_z_1 = pbuffer.data(idx_eri_1_gp + 29);

    auto g_yyyy_x_1 = pbuffer.data(idx_eri_1_gp + 30);

    auto g_yyyy_y_1 = pbuffer.data(idx_eri_1_gp + 31);

    auto g_yyyy_z_1 = pbuffer.data(idx_eri_1_gp + 32);

    auto g_yyyz_z_1 = pbuffer.data(idx_eri_1_gp + 35);

    auto g_yyzz_x_1 = pbuffer.data(idx_eri_1_gp + 36);

    auto g_yyzz_y_1 = pbuffer.data(idx_eri_1_gp + 37);

    auto g_yyzz_z_1 = pbuffer.data(idx_eri_1_gp + 38);

    auto g_yzzz_y_1 = pbuffer.data(idx_eri_1_gp + 40);

    auto g_yzzz_z_1 = pbuffer.data(idx_eri_1_gp + 41);

    auto g_zzzz_x_1 = pbuffer.data(idx_eri_1_gp + 42);

    auto g_zzzz_y_1 = pbuffer.data(idx_eri_1_gp + 43);

    auto g_zzzz_z_1 = pbuffer.data(idx_eri_1_gp + 44);

    // Set up components of auxiliary buffer : GD

    auto g_xxxx_xx_1 = pbuffer.data(idx_eri_1_gd);

    auto g_xxxx_xy_1 = pbuffer.data(idx_eri_1_gd + 1);

    auto g_xxxx_xz_1 = pbuffer.data(idx_eri_1_gd + 2);

    auto g_xxxx_yy_1 = pbuffer.data(idx_eri_1_gd + 3);

    auto g_xxxx_yz_1 = pbuffer.data(idx_eri_1_gd + 4);

    auto g_xxxx_zz_1 = pbuffer.data(idx_eri_1_gd + 5);

    auto g_xxxy_xx_1 = pbuffer.data(idx_eri_1_gd + 6);

    auto g_xxxy_xy_1 = pbuffer.data(idx_eri_1_gd + 7);

    auto g_xxxy_xz_1 = pbuffer.data(idx_eri_1_gd + 8);

    auto g_xxxy_yy_1 = pbuffer.data(idx_eri_1_gd + 9);

    auto g_xxxz_xx_1 = pbuffer.data(idx_eri_1_gd + 12);

    auto g_xxxz_xy_1 = pbuffer.data(idx_eri_1_gd + 13);

    auto g_xxxz_xz_1 = pbuffer.data(idx_eri_1_gd + 14);

    auto g_xxxz_yz_1 = pbuffer.data(idx_eri_1_gd + 16);

    auto g_xxxz_zz_1 = pbuffer.data(idx_eri_1_gd + 17);

    auto g_xxyy_xx_1 = pbuffer.data(idx_eri_1_gd + 18);

    auto g_xxyy_xy_1 = pbuffer.data(idx_eri_1_gd + 19);

    auto g_xxyy_xz_1 = pbuffer.data(idx_eri_1_gd + 20);

    auto g_xxyy_yy_1 = pbuffer.data(idx_eri_1_gd + 21);

    auto g_xxyy_yz_1 = pbuffer.data(idx_eri_1_gd + 22);

    auto g_xxyy_zz_1 = pbuffer.data(idx_eri_1_gd + 23);

    auto g_xxzz_xx_1 = pbuffer.data(idx_eri_1_gd + 30);

    auto g_xxzz_xy_1 = pbuffer.data(idx_eri_1_gd + 31);

    auto g_xxzz_xz_1 = pbuffer.data(idx_eri_1_gd + 32);

    auto g_xxzz_yy_1 = pbuffer.data(idx_eri_1_gd + 33);

    auto g_xxzz_yz_1 = pbuffer.data(idx_eri_1_gd + 34);

    auto g_xxzz_zz_1 = pbuffer.data(idx_eri_1_gd + 35);

    auto g_xyyy_xx_1 = pbuffer.data(idx_eri_1_gd + 36);

    auto g_xyyy_xy_1 = pbuffer.data(idx_eri_1_gd + 37);

    auto g_xyyy_yy_1 = pbuffer.data(idx_eri_1_gd + 39);

    auto g_xyyy_yz_1 = pbuffer.data(idx_eri_1_gd + 40);

    auto g_xyyy_zz_1 = pbuffer.data(idx_eri_1_gd + 41);

    auto g_xzzz_xx_1 = pbuffer.data(idx_eri_1_gd + 54);

    auto g_xzzz_xz_1 = pbuffer.data(idx_eri_1_gd + 56);

    auto g_xzzz_yy_1 = pbuffer.data(idx_eri_1_gd + 57);

    auto g_xzzz_yz_1 = pbuffer.data(idx_eri_1_gd + 58);

    auto g_xzzz_zz_1 = pbuffer.data(idx_eri_1_gd + 59);

    auto g_yyyy_xx_1 = pbuffer.data(idx_eri_1_gd + 60);

    auto g_yyyy_xy_1 = pbuffer.data(idx_eri_1_gd + 61);

    auto g_yyyy_xz_1 = pbuffer.data(idx_eri_1_gd + 62);

    auto g_yyyy_yy_1 = pbuffer.data(idx_eri_1_gd + 63);

    auto g_yyyy_yz_1 = pbuffer.data(idx_eri_1_gd + 64);

    auto g_yyyy_zz_1 = pbuffer.data(idx_eri_1_gd + 65);

    auto g_yyyz_xy_1 = pbuffer.data(idx_eri_1_gd + 67);

    auto g_yyyz_xz_1 = pbuffer.data(idx_eri_1_gd + 68);

    auto g_yyyz_yy_1 = pbuffer.data(idx_eri_1_gd + 69);

    auto g_yyyz_yz_1 = pbuffer.data(idx_eri_1_gd + 70);

    auto g_yyyz_zz_1 = pbuffer.data(idx_eri_1_gd + 71);

    auto g_yyzz_xx_1 = pbuffer.data(idx_eri_1_gd + 72);

    auto g_yyzz_xy_1 = pbuffer.data(idx_eri_1_gd + 73);

    auto g_yyzz_xz_1 = pbuffer.data(idx_eri_1_gd + 74);

    auto g_yyzz_yy_1 = pbuffer.data(idx_eri_1_gd + 75);

    auto g_yyzz_yz_1 = pbuffer.data(idx_eri_1_gd + 76);

    auto g_yyzz_zz_1 = pbuffer.data(idx_eri_1_gd + 77);

    auto g_yzzz_xx_1 = pbuffer.data(idx_eri_1_gd + 78);

    auto g_yzzz_xy_1 = pbuffer.data(idx_eri_1_gd + 79);

    auto g_yzzz_xz_1 = pbuffer.data(idx_eri_1_gd + 80);

    auto g_yzzz_yy_1 = pbuffer.data(idx_eri_1_gd + 81);

    auto g_yzzz_yz_1 = pbuffer.data(idx_eri_1_gd + 82);

    auto g_yzzz_zz_1 = pbuffer.data(idx_eri_1_gd + 83);

    auto g_zzzz_xx_1 = pbuffer.data(idx_eri_1_gd + 84);

    auto g_zzzz_xy_1 = pbuffer.data(idx_eri_1_gd + 85);

    auto g_zzzz_xz_1 = pbuffer.data(idx_eri_1_gd + 86);

    auto g_zzzz_yy_1 = pbuffer.data(idx_eri_1_gd + 87);

    auto g_zzzz_yz_1 = pbuffer.data(idx_eri_1_gd + 88);

    auto g_zzzz_zz_1 = pbuffer.data(idx_eri_1_gd + 89);

    // Set up 0-6 components of targeted buffer : HD

    auto g_xxxxx_xx_0 = pbuffer.data(idx_eri_0_hd);

    auto g_xxxxx_xy_0 = pbuffer.data(idx_eri_0_hd + 1);

    auto g_xxxxx_xz_0 = pbuffer.data(idx_eri_0_hd + 2);

    auto g_xxxxx_yy_0 = pbuffer.data(idx_eri_0_hd + 3);

    auto g_xxxxx_yz_0 = pbuffer.data(idx_eri_0_hd + 4);

    auto g_xxxxx_zz_0 = pbuffer.data(idx_eri_0_hd + 5);

    #pragma omp simd aligned(g_xxx_xx_0, g_xxx_xx_1, g_xxx_xy_0, g_xxx_xy_1, g_xxx_xz_0, g_xxx_xz_1, g_xxx_yy_0, g_xxx_yy_1, g_xxx_yz_0, g_xxx_yz_1, g_xxx_zz_0, g_xxx_zz_1, g_xxxx_x_1, g_xxxx_xx_1, g_xxxx_xy_1, g_xxxx_xz_1, g_xxxx_y_1, g_xxxx_yy_1, g_xxxx_yz_1, g_xxxx_z_1, g_xxxx_zz_1, g_xxxxx_xx_0, g_xxxxx_xy_0, g_xxxxx_xz_0, g_xxxxx_yy_0, g_xxxxx_yz_0, g_xxxxx_zz_0, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxxx_xx_0[i] = 4.0 * g_xxx_xx_0[i] * fbe_0 - 4.0 * g_xxx_xx_1[i] * fz_be_0 + 2.0 * g_xxxx_x_1[i] * fe_0 + g_xxxx_xx_1[i] * pa_x[i];

        g_xxxxx_xy_0[i] = 4.0 * g_xxx_xy_0[i] * fbe_0 - 4.0 * g_xxx_xy_1[i] * fz_be_0 + g_xxxx_y_1[i] * fe_0 + g_xxxx_xy_1[i] * pa_x[i];

        g_xxxxx_xz_0[i] = 4.0 * g_xxx_xz_0[i] * fbe_0 - 4.0 * g_xxx_xz_1[i] * fz_be_0 + g_xxxx_z_1[i] * fe_0 + g_xxxx_xz_1[i] * pa_x[i];

        g_xxxxx_yy_0[i] = 4.0 * g_xxx_yy_0[i] * fbe_0 - 4.0 * g_xxx_yy_1[i] * fz_be_0 + g_xxxx_yy_1[i] * pa_x[i];

        g_xxxxx_yz_0[i] = 4.0 * g_xxx_yz_0[i] * fbe_0 - 4.0 * g_xxx_yz_1[i] * fz_be_0 + g_xxxx_yz_1[i] * pa_x[i];

        g_xxxxx_zz_0[i] = 4.0 * g_xxx_zz_0[i] * fbe_0 - 4.0 * g_xxx_zz_1[i] * fz_be_0 + g_xxxx_zz_1[i] * pa_x[i];
    }

    // Set up 6-12 components of targeted buffer : HD

    auto g_xxxxy_xx_0 = pbuffer.data(idx_eri_0_hd + 6);

    auto g_xxxxy_xy_0 = pbuffer.data(idx_eri_0_hd + 7);

    auto g_xxxxy_xz_0 = pbuffer.data(idx_eri_0_hd + 8);

    auto g_xxxxy_yy_0 = pbuffer.data(idx_eri_0_hd + 9);

    auto g_xxxxy_yz_0 = pbuffer.data(idx_eri_0_hd + 10);

    auto g_xxxxy_zz_0 = pbuffer.data(idx_eri_0_hd + 11);

    #pragma omp simd aligned(g_xxxx_x_1, g_xxxx_xx_1, g_xxxx_xy_1, g_xxxx_xz_1, g_xxxx_y_1, g_xxxx_yy_1, g_xxxx_yz_1, g_xxxx_z_1, g_xxxx_zz_1, g_xxxxy_xx_0, g_xxxxy_xy_0, g_xxxxy_xz_0, g_xxxxy_yy_0, g_xxxxy_yz_0, g_xxxxy_zz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxy_xx_0[i] = g_xxxx_xx_1[i] * pa_y[i];

        g_xxxxy_xy_0[i] = g_xxxx_x_1[i] * fe_0 + g_xxxx_xy_1[i] * pa_y[i];

        g_xxxxy_xz_0[i] = g_xxxx_xz_1[i] * pa_y[i];

        g_xxxxy_yy_0[i] = 2.0 * g_xxxx_y_1[i] * fe_0 + g_xxxx_yy_1[i] * pa_y[i];

        g_xxxxy_yz_0[i] = g_xxxx_z_1[i] * fe_0 + g_xxxx_yz_1[i] * pa_y[i];

        g_xxxxy_zz_0[i] = g_xxxx_zz_1[i] * pa_y[i];
    }

    // Set up 12-18 components of targeted buffer : HD

    auto g_xxxxz_xx_0 = pbuffer.data(idx_eri_0_hd + 12);

    auto g_xxxxz_xy_0 = pbuffer.data(idx_eri_0_hd + 13);

    auto g_xxxxz_xz_0 = pbuffer.data(idx_eri_0_hd + 14);

    auto g_xxxxz_yy_0 = pbuffer.data(idx_eri_0_hd + 15);

    auto g_xxxxz_yz_0 = pbuffer.data(idx_eri_0_hd + 16);

    auto g_xxxxz_zz_0 = pbuffer.data(idx_eri_0_hd + 17);

    #pragma omp simd aligned(g_xxxx_x_1, g_xxxx_xx_1, g_xxxx_xy_1, g_xxxx_xz_1, g_xxxx_y_1, g_xxxx_yy_1, g_xxxx_yz_1, g_xxxx_z_1, g_xxxx_zz_1, g_xxxxz_xx_0, g_xxxxz_xy_0, g_xxxxz_xz_0, g_xxxxz_yy_0, g_xxxxz_yz_0, g_xxxxz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxxz_xx_0[i] = g_xxxx_xx_1[i] * pa_z[i];

        g_xxxxz_xy_0[i] = g_xxxx_xy_1[i] * pa_z[i];

        g_xxxxz_xz_0[i] = g_xxxx_x_1[i] * fe_0 + g_xxxx_xz_1[i] * pa_z[i];

        g_xxxxz_yy_0[i] = g_xxxx_yy_1[i] * pa_z[i];

        g_xxxxz_yz_0[i] = g_xxxx_y_1[i] * fe_0 + g_xxxx_yz_1[i] * pa_z[i];

        g_xxxxz_zz_0[i] = 2.0 * g_xxxx_z_1[i] * fe_0 + g_xxxx_zz_1[i] * pa_z[i];
    }

    // Set up 18-24 components of targeted buffer : HD

    auto g_xxxyy_xx_0 = pbuffer.data(idx_eri_0_hd + 18);

    auto g_xxxyy_xy_0 = pbuffer.data(idx_eri_0_hd + 19);

    auto g_xxxyy_xz_0 = pbuffer.data(idx_eri_0_hd + 20);

    auto g_xxxyy_yy_0 = pbuffer.data(idx_eri_0_hd + 21);

    auto g_xxxyy_yz_0 = pbuffer.data(idx_eri_0_hd + 22);

    auto g_xxxyy_zz_0 = pbuffer.data(idx_eri_0_hd + 23);

    #pragma omp simd aligned(g_xxx_xx_0, g_xxx_xx_1, g_xxx_xz_0, g_xxx_xz_1, g_xxxy_xx_1, g_xxxy_xz_1, g_xxxyy_xx_0, g_xxxyy_xy_0, g_xxxyy_xz_0, g_xxxyy_yy_0, g_xxxyy_yz_0, g_xxxyy_zz_0, g_xxyy_xy_1, g_xxyy_y_1, g_xxyy_yy_1, g_xxyy_yz_1, g_xxyy_zz_1, g_xyy_xy_0, g_xyy_xy_1, g_xyy_yy_0, g_xyy_yy_1, g_xyy_yz_0, g_xyy_yz_1, g_xyy_zz_0, g_xyy_zz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxyy_xx_0[i] = g_xxx_xx_0[i] * fbe_0 - g_xxx_xx_1[i] * fz_be_0 + g_xxxy_xx_1[i] * pa_y[i];

        g_xxxyy_xy_0[i] = 2.0 * g_xyy_xy_0[i] * fbe_0 - 2.0 * g_xyy_xy_1[i] * fz_be_0 + g_xxyy_y_1[i] * fe_0 + g_xxyy_xy_1[i] * pa_x[i];

        g_xxxyy_xz_0[i] = g_xxx_xz_0[i] * fbe_0 - g_xxx_xz_1[i] * fz_be_0 + g_xxxy_xz_1[i] * pa_y[i];

        g_xxxyy_yy_0[i] = 2.0 * g_xyy_yy_0[i] * fbe_0 - 2.0 * g_xyy_yy_1[i] * fz_be_0 + g_xxyy_yy_1[i] * pa_x[i];

        g_xxxyy_yz_0[i] = 2.0 * g_xyy_yz_0[i] * fbe_0 - 2.0 * g_xyy_yz_1[i] * fz_be_0 + g_xxyy_yz_1[i] * pa_x[i];

        g_xxxyy_zz_0[i] = 2.0 * g_xyy_zz_0[i] * fbe_0 - 2.0 * g_xyy_zz_1[i] * fz_be_0 + g_xxyy_zz_1[i] * pa_x[i];
    }

    // Set up 24-30 components of targeted buffer : HD

    auto g_xxxyz_xx_0 = pbuffer.data(idx_eri_0_hd + 24);

    auto g_xxxyz_xy_0 = pbuffer.data(idx_eri_0_hd + 25);

    auto g_xxxyz_xz_0 = pbuffer.data(idx_eri_0_hd + 26);

    auto g_xxxyz_yy_0 = pbuffer.data(idx_eri_0_hd + 27);

    auto g_xxxyz_yz_0 = pbuffer.data(idx_eri_0_hd + 28);

    auto g_xxxyz_zz_0 = pbuffer.data(idx_eri_0_hd + 29);

    #pragma omp simd aligned(g_xxxy_xy_1, g_xxxy_yy_1, g_xxxyz_xx_0, g_xxxyz_xy_0, g_xxxyz_xz_0, g_xxxyz_yy_0, g_xxxyz_yz_0, g_xxxyz_zz_0, g_xxxz_xx_1, g_xxxz_xz_1, g_xxxz_yz_1, g_xxxz_z_1, g_xxxz_zz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxxyz_xx_0[i] = g_xxxz_xx_1[i] * pa_y[i];

        g_xxxyz_xy_0[i] = g_xxxy_xy_1[i] * pa_z[i];

        g_xxxyz_xz_0[i] = g_xxxz_xz_1[i] * pa_y[i];

        g_xxxyz_yy_0[i] = g_xxxy_yy_1[i] * pa_z[i];

        g_xxxyz_yz_0[i] = g_xxxz_z_1[i] * fe_0 + g_xxxz_yz_1[i] * pa_y[i];

        g_xxxyz_zz_0[i] = g_xxxz_zz_1[i] * pa_y[i];
    }

    // Set up 30-36 components of targeted buffer : HD

    auto g_xxxzz_xx_0 = pbuffer.data(idx_eri_0_hd + 30);

    auto g_xxxzz_xy_0 = pbuffer.data(idx_eri_0_hd + 31);

    auto g_xxxzz_xz_0 = pbuffer.data(idx_eri_0_hd + 32);

    auto g_xxxzz_yy_0 = pbuffer.data(idx_eri_0_hd + 33);

    auto g_xxxzz_yz_0 = pbuffer.data(idx_eri_0_hd + 34);

    auto g_xxxzz_zz_0 = pbuffer.data(idx_eri_0_hd + 35);

    #pragma omp simd aligned(g_xxx_xx_0, g_xxx_xx_1, g_xxx_xy_0, g_xxx_xy_1, g_xxxz_xx_1, g_xxxz_xy_1, g_xxxzz_xx_0, g_xxxzz_xy_0, g_xxxzz_xz_0, g_xxxzz_yy_0, g_xxxzz_yz_0, g_xxxzz_zz_0, g_xxzz_xz_1, g_xxzz_yy_1, g_xxzz_yz_1, g_xxzz_z_1, g_xxzz_zz_1, g_xzz_xz_0, g_xzz_xz_1, g_xzz_yy_0, g_xzz_yy_1, g_xzz_yz_0, g_xzz_yz_1, g_xzz_zz_0, g_xzz_zz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxxzz_xx_0[i] = g_xxx_xx_0[i] * fbe_0 - g_xxx_xx_1[i] * fz_be_0 + g_xxxz_xx_1[i] * pa_z[i];

        g_xxxzz_xy_0[i] = g_xxx_xy_0[i] * fbe_0 - g_xxx_xy_1[i] * fz_be_0 + g_xxxz_xy_1[i] * pa_z[i];

        g_xxxzz_xz_0[i] = 2.0 * g_xzz_xz_0[i] * fbe_0 - 2.0 * g_xzz_xz_1[i] * fz_be_0 + g_xxzz_z_1[i] * fe_0 + g_xxzz_xz_1[i] * pa_x[i];

        g_xxxzz_yy_0[i] = 2.0 * g_xzz_yy_0[i] * fbe_0 - 2.0 * g_xzz_yy_1[i] * fz_be_0 + g_xxzz_yy_1[i] * pa_x[i];

        g_xxxzz_yz_0[i] = 2.0 * g_xzz_yz_0[i] * fbe_0 - 2.0 * g_xzz_yz_1[i] * fz_be_0 + g_xxzz_yz_1[i] * pa_x[i];

        g_xxxzz_zz_0[i] = 2.0 * g_xzz_zz_0[i] * fbe_0 - 2.0 * g_xzz_zz_1[i] * fz_be_0 + g_xxzz_zz_1[i] * pa_x[i];
    }

    // Set up 36-42 components of targeted buffer : HD

    auto g_xxyyy_xx_0 = pbuffer.data(idx_eri_0_hd + 36);

    auto g_xxyyy_xy_0 = pbuffer.data(idx_eri_0_hd + 37);

    auto g_xxyyy_xz_0 = pbuffer.data(idx_eri_0_hd + 38);

    auto g_xxyyy_yy_0 = pbuffer.data(idx_eri_0_hd + 39);

    auto g_xxyyy_yz_0 = pbuffer.data(idx_eri_0_hd + 40);

    auto g_xxyyy_zz_0 = pbuffer.data(idx_eri_0_hd + 41);

    #pragma omp simd aligned(g_xxy_xx_0, g_xxy_xx_1, g_xxy_xz_0, g_xxy_xz_1, g_xxyy_xx_1, g_xxyy_xz_1, g_xxyyy_xx_0, g_xxyyy_xy_0, g_xxyyy_xz_0, g_xxyyy_yy_0, g_xxyyy_yz_0, g_xxyyy_zz_0, g_xyyy_xy_1, g_xyyy_y_1, g_xyyy_yy_1, g_xyyy_yz_1, g_xyyy_zz_1, g_yyy_xy_0, g_yyy_xy_1, g_yyy_yy_0, g_yyy_yy_1, g_yyy_yz_0, g_yyy_yz_1, g_yyy_zz_0, g_yyy_zz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxyyy_xx_0[i] = 2.0 * g_xxy_xx_0[i] * fbe_0 - 2.0 * g_xxy_xx_1[i] * fz_be_0 + g_xxyy_xx_1[i] * pa_y[i];

        g_xxyyy_xy_0[i] = g_yyy_xy_0[i] * fbe_0 - g_yyy_xy_1[i] * fz_be_0 + g_xyyy_y_1[i] * fe_0 + g_xyyy_xy_1[i] * pa_x[i];

        g_xxyyy_xz_0[i] = 2.0 * g_xxy_xz_0[i] * fbe_0 - 2.0 * g_xxy_xz_1[i] * fz_be_0 + g_xxyy_xz_1[i] * pa_y[i];

        g_xxyyy_yy_0[i] = g_yyy_yy_0[i] * fbe_0 - g_yyy_yy_1[i] * fz_be_0 + g_xyyy_yy_1[i] * pa_x[i];

        g_xxyyy_yz_0[i] = g_yyy_yz_0[i] * fbe_0 - g_yyy_yz_1[i] * fz_be_0 + g_xyyy_yz_1[i] * pa_x[i];

        g_xxyyy_zz_0[i] = g_yyy_zz_0[i] * fbe_0 - g_yyy_zz_1[i] * fz_be_0 + g_xyyy_zz_1[i] * pa_x[i];
    }

    // Set up 42-48 components of targeted buffer : HD

    auto g_xxyyz_xx_0 = pbuffer.data(idx_eri_0_hd + 42);

    auto g_xxyyz_xy_0 = pbuffer.data(idx_eri_0_hd + 43);

    auto g_xxyyz_xz_0 = pbuffer.data(idx_eri_0_hd + 44);

    auto g_xxyyz_yy_0 = pbuffer.data(idx_eri_0_hd + 45);

    auto g_xxyyz_yz_0 = pbuffer.data(idx_eri_0_hd + 46);

    auto g_xxyyz_zz_0 = pbuffer.data(idx_eri_0_hd + 47);

    #pragma omp simd aligned(g_xxyy_x_1, g_xxyy_xx_1, g_xxyy_xy_1, g_xxyy_xz_1, g_xxyy_y_1, g_xxyy_yy_1, g_xxyy_yz_1, g_xxyy_z_1, g_xxyy_zz_1, g_xxyyz_xx_0, g_xxyyz_xy_0, g_xxyyz_xz_0, g_xxyyz_yy_0, g_xxyyz_yz_0, g_xxyyz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyyz_xx_0[i] = g_xxyy_xx_1[i] * pa_z[i];

        g_xxyyz_xy_0[i] = g_xxyy_xy_1[i] * pa_z[i];

        g_xxyyz_xz_0[i] = g_xxyy_x_1[i] * fe_0 + g_xxyy_xz_1[i] * pa_z[i];

        g_xxyyz_yy_0[i] = g_xxyy_yy_1[i] * pa_z[i];

        g_xxyyz_yz_0[i] = g_xxyy_y_1[i] * fe_0 + g_xxyy_yz_1[i] * pa_z[i];

        g_xxyyz_zz_0[i] = 2.0 * g_xxyy_z_1[i] * fe_0 + g_xxyy_zz_1[i] * pa_z[i];
    }

    // Set up 48-54 components of targeted buffer : HD

    auto g_xxyzz_xx_0 = pbuffer.data(idx_eri_0_hd + 48);

    auto g_xxyzz_xy_0 = pbuffer.data(idx_eri_0_hd + 49);

    auto g_xxyzz_xz_0 = pbuffer.data(idx_eri_0_hd + 50);

    auto g_xxyzz_yy_0 = pbuffer.data(idx_eri_0_hd + 51);

    auto g_xxyzz_yz_0 = pbuffer.data(idx_eri_0_hd + 52);

    auto g_xxyzz_zz_0 = pbuffer.data(idx_eri_0_hd + 53);

    #pragma omp simd aligned(g_xxyzz_xx_0, g_xxyzz_xy_0, g_xxyzz_xz_0, g_xxyzz_yy_0, g_xxyzz_yz_0, g_xxyzz_zz_0, g_xxzz_x_1, g_xxzz_xx_1, g_xxzz_xy_1, g_xxzz_xz_1, g_xxzz_y_1, g_xxzz_yy_1, g_xxzz_yz_1, g_xxzz_z_1, g_xxzz_zz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xxyzz_xx_0[i] = g_xxzz_xx_1[i] * pa_y[i];

        g_xxyzz_xy_0[i] = g_xxzz_x_1[i] * fe_0 + g_xxzz_xy_1[i] * pa_y[i];

        g_xxyzz_xz_0[i] = g_xxzz_xz_1[i] * pa_y[i];

        g_xxyzz_yy_0[i] = 2.0 * g_xxzz_y_1[i] * fe_0 + g_xxzz_yy_1[i] * pa_y[i];

        g_xxyzz_yz_0[i] = g_xxzz_z_1[i] * fe_0 + g_xxzz_yz_1[i] * pa_y[i];

        g_xxyzz_zz_0[i] = g_xxzz_zz_1[i] * pa_y[i];
    }

    // Set up 54-60 components of targeted buffer : HD

    auto g_xxzzz_xx_0 = pbuffer.data(idx_eri_0_hd + 54);

    auto g_xxzzz_xy_0 = pbuffer.data(idx_eri_0_hd + 55);

    auto g_xxzzz_xz_0 = pbuffer.data(idx_eri_0_hd + 56);

    auto g_xxzzz_yy_0 = pbuffer.data(idx_eri_0_hd + 57);

    auto g_xxzzz_yz_0 = pbuffer.data(idx_eri_0_hd + 58);

    auto g_xxzzz_zz_0 = pbuffer.data(idx_eri_0_hd + 59);

    #pragma omp simd aligned(g_xxz_xx_0, g_xxz_xx_1, g_xxz_xy_0, g_xxz_xy_1, g_xxzz_xx_1, g_xxzz_xy_1, g_xxzzz_xx_0, g_xxzzz_xy_0, g_xxzzz_xz_0, g_xxzzz_yy_0, g_xxzzz_yz_0, g_xxzzz_zz_0, g_xzzz_xz_1, g_xzzz_yy_1, g_xzzz_yz_1, g_xzzz_z_1, g_xzzz_zz_1, g_zzz_xz_0, g_zzz_xz_1, g_zzz_yy_0, g_zzz_yy_1, g_zzz_yz_0, g_zzz_yz_1, g_zzz_zz_0, g_zzz_zz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_xxzzz_xx_0[i] = 2.0 * g_xxz_xx_0[i] * fbe_0 - 2.0 * g_xxz_xx_1[i] * fz_be_0 + g_xxzz_xx_1[i] * pa_z[i];

        g_xxzzz_xy_0[i] = 2.0 * g_xxz_xy_0[i] * fbe_0 - 2.0 * g_xxz_xy_1[i] * fz_be_0 + g_xxzz_xy_1[i] * pa_z[i];

        g_xxzzz_xz_0[i] = g_zzz_xz_0[i] * fbe_0 - g_zzz_xz_1[i] * fz_be_0 + g_xzzz_z_1[i] * fe_0 + g_xzzz_xz_1[i] * pa_x[i];

        g_xxzzz_yy_0[i] = g_zzz_yy_0[i] * fbe_0 - g_zzz_yy_1[i] * fz_be_0 + g_xzzz_yy_1[i] * pa_x[i];

        g_xxzzz_yz_0[i] = g_zzz_yz_0[i] * fbe_0 - g_zzz_yz_1[i] * fz_be_0 + g_xzzz_yz_1[i] * pa_x[i];

        g_xxzzz_zz_0[i] = g_zzz_zz_0[i] * fbe_0 - g_zzz_zz_1[i] * fz_be_0 + g_xzzz_zz_1[i] * pa_x[i];
    }

    // Set up 60-66 components of targeted buffer : HD

    auto g_xyyyy_xx_0 = pbuffer.data(idx_eri_0_hd + 60);

    auto g_xyyyy_xy_0 = pbuffer.data(idx_eri_0_hd + 61);

    auto g_xyyyy_xz_0 = pbuffer.data(idx_eri_0_hd + 62);

    auto g_xyyyy_yy_0 = pbuffer.data(idx_eri_0_hd + 63);

    auto g_xyyyy_yz_0 = pbuffer.data(idx_eri_0_hd + 64);

    auto g_xyyyy_zz_0 = pbuffer.data(idx_eri_0_hd + 65);

    #pragma omp simd aligned(g_xyyyy_xx_0, g_xyyyy_xy_0, g_xyyyy_xz_0, g_xyyyy_yy_0, g_xyyyy_yz_0, g_xyyyy_zz_0, g_yyyy_x_1, g_yyyy_xx_1, g_yyyy_xy_1, g_yyyy_xz_1, g_yyyy_y_1, g_yyyy_yy_1, g_yyyy_yz_1, g_yyyy_z_1, g_yyyy_zz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyy_xx_0[i] = 2.0 * g_yyyy_x_1[i] * fe_0 + g_yyyy_xx_1[i] * pa_x[i];

        g_xyyyy_xy_0[i] = g_yyyy_y_1[i] * fe_0 + g_yyyy_xy_1[i] * pa_x[i];

        g_xyyyy_xz_0[i] = g_yyyy_z_1[i] * fe_0 + g_yyyy_xz_1[i] * pa_x[i];

        g_xyyyy_yy_0[i] = g_yyyy_yy_1[i] * pa_x[i];

        g_xyyyy_yz_0[i] = g_yyyy_yz_1[i] * pa_x[i];

        g_xyyyy_zz_0[i] = g_yyyy_zz_1[i] * pa_x[i];
    }

    // Set up 66-72 components of targeted buffer : HD

    auto g_xyyyz_xx_0 = pbuffer.data(idx_eri_0_hd + 66);

    auto g_xyyyz_xy_0 = pbuffer.data(idx_eri_0_hd + 67);

    auto g_xyyyz_xz_0 = pbuffer.data(idx_eri_0_hd + 68);

    auto g_xyyyz_yy_0 = pbuffer.data(idx_eri_0_hd + 69);

    auto g_xyyyz_yz_0 = pbuffer.data(idx_eri_0_hd + 70);

    auto g_xyyyz_zz_0 = pbuffer.data(idx_eri_0_hd + 71);

    #pragma omp simd aligned(g_xyyy_xx_1, g_xyyy_xy_1, g_xyyyz_xx_0, g_xyyyz_xy_0, g_xyyyz_xz_0, g_xyyyz_yy_0, g_xyyyz_yz_0, g_xyyyz_zz_0, g_yyyz_xz_1, g_yyyz_yy_1, g_yyyz_yz_1, g_yyyz_z_1, g_yyyz_zz_1, pa_x, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyyz_xx_0[i] = g_xyyy_xx_1[i] * pa_z[i];

        g_xyyyz_xy_0[i] = g_xyyy_xy_1[i] * pa_z[i];

        g_xyyyz_xz_0[i] = g_yyyz_z_1[i] * fe_0 + g_yyyz_xz_1[i] * pa_x[i];

        g_xyyyz_yy_0[i] = g_yyyz_yy_1[i] * pa_x[i];

        g_xyyyz_yz_0[i] = g_yyyz_yz_1[i] * pa_x[i];

        g_xyyyz_zz_0[i] = g_yyyz_zz_1[i] * pa_x[i];
    }

    // Set up 72-78 components of targeted buffer : HD

    auto g_xyyzz_xx_0 = pbuffer.data(idx_eri_0_hd + 72);

    auto g_xyyzz_xy_0 = pbuffer.data(idx_eri_0_hd + 73);

    auto g_xyyzz_xz_0 = pbuffer.data(idx_eri_0_hd + 74);

    auto g_xyyzz_yy_0 = pbuffer.data(idx_eri_0_hd + 75);

    auto g_xyyzz_yz_0 = pbuffer.data(idx_eri_0_hd + 76);

    auto g_xyyzz_zz_0 = pbuffer.data(idx_eri_0_hd + 77);

    #pragma omp simd aligned(g_xyyzz_xx_0, g_xyyzz_xy_0, g_xyyzz_xz_0, g_xyyzz_yy_0, g_xyyzz_yz_0, g_xyyzz_zz_0, g_yyzz_x_1, g_yyzz_xx_1, g_yyzz_xy_1, g_yyzz_xz_1, g_yyzz_y_1, g_yyzz_yy_1, g_yyzz_yz_1, g_yyzz_z_1, g_yyzz_zz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyyzz_xx_0[i] = 2.0 * g_yyzz_x_1[i] * fe_0 + g_yyzz_xx_1[i] * pa_x[i];

        g_xyyzz_xy_0[i] = g_yyzz_y_1[i] * fe_0 + g_yyzz_xy_1[i] * pa_x[i];

        g_xyyzz_xz_0[i] = g_yyzz_z_1[i] * fe_0 + g_yyzz_xz_1[i] * pa_x[i];

        g_xyyzz_yy_0[i] = g_yyzz_yy_1[i] * pa_x[i];

        g_xyyzz_yz_0[i] = g_yyzz_yz_1[i] * pa_x[i];

        g_xyyzz_zz_0[i] = g_yyzz_zz_1[i] * pa_x[i];
    }

    // Set up 78-84 components of targeted buffer : HD

    auto g_xyzzz_xx_0 = pbuffer.data(idx_eri_0_hd + 78);

    auto g_xyzzz_xy_0 = pbuffer.data(idx_eri_0_hd + 79);

    auto g_xyzzz_xz_0 = pbuffer.data(idx_eri_0_hd + 80);

    auto g_xyzzz_yy_0 = pbuffer.data(idx_eri_0_hd + 81);

    auto g_xyzzz_yz_0 = pbuffer.data(idx_eri_0_hd + 82);

    auto g_xyzzz_zz_0 = pbuffer.data(idx_eri_0_hd + 83);

    #pragma omp simd aligned(g_xyzzz_xx_0, g_xyzzz_xy_0, g_xyzzz_xz_0, g_xyzzz_yy_0, g_xyzzz_yz_0, g_xyzzz_zz_0, g_xzzz_xx_1, g_xzzz_xz_1, g_yzzz_xy_1, g_yzzz_y_1, g_yzzz_yy_1, g_yzzz_yz_1, g_yzzz_zz_1, pa_x, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xyzzz_xx_0[i] = g_xzzz_xx_1[i] * pa_y[i];

        g_xyzzz_xy_0[i] = g_yzzz_y_1[i] * fe_0 + g_yzzz_xy_1[i] * pa_x[i];

        g_xyzzz_xz_0[i] = g_xzzz_xz_1[i] * pa_y[i];

        g_xyzzz_yy_0[i] = g_yzzz_yy_1[i] * pa_x[i];

        g_xyzzz_yz_0[i] = g_yzzz_yz_1[i] * pa_x[i];

        g_xyzzz_zz_0[i] = g_yzzz_zz_1[i] * pa_x[i];
    }

    // Set up 84-90 components of targeted buffer : HD

    auto g_xzzzz_xx_0 = pbuffer.data(idx_eri_0_hd + 84);

    auto g_xzzzz_xy_0 = pbuffer.data(idx_eri_0_hd + 85);

    auto g_xzzzz_xz_0 = pbuffer.data(idx_eri_0_hd + 86);

    auto g_xzzzz_yy_0 = pbuffer.data(idx_eri_0_hd + 87);

    auto g_xzzzz_yz_0 = pbuffer.data(idx_eri_0_hd + 88);

    auto g_xzzzz_zz_0 = pbuffer.data(idx_eri_0_hd + 89);

    #pragma omp simd aligned(g_xzzzz_xx_0, g_xzzzz_xy_0, g_xzzzz_xz_0, g_xzzzz_yy_0, g_xzzzz_yz_0, g_xzzzz_zz_0, g_zzzz_x_1, g_zzzz_xx_1, g_zzzz_xy_1, g_zzzz_xz_1, g_zzzz_y_1, g_zzzz_yy_1, g_zzzz_yz_1, g_zzzz_z_1, g_zzzz_zz_1, pa_x, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_xzzzz_xx_0[i] = 2.0 * g_zzzz_x_1[i] * fe_0 + g_zzzz_xx_1[i] * pa_x[i];

        g_xzzzz_xy_0[i] = g_zzzz_y_1[i] * fe_0 + g_zzzz_xy_1[i] * pa_x[i];

        g_xzzzz_xz_0[i] = g_zzzz_z_1[i] * fe_0 + g_zzzz_xz_1[i] * pa_x[i];

        g_xzzzz_yy_0[i] = g_zzzz_yy_1[i] * pa_x[i];

        g_xzzzz_yz_0[i] = g_zzzz_yz_1[i] * pa_x[i];

        g_xzzzz_zz_0[i] = g_zzzz_zz_1[i] * pa_x[i];
    }

    // Set up 90-96 components of targeted buffer : HD

    auto g_yyyyy_xx_0 = pbuffer.data(idx_eri_0_hd + 90);

    auto g_yyyyy_xy_0 = pbuffer.data(idx_eri_0_hd + 91);

    auto g_yyyyy_xz_0 = pbuffer.data(idx_eri_0_hd + 92);

    auto g_yyyyy_yy_0 = pbuffer.data(idx_eri_0_hd + 93);

    auto g_yyyyy_yz_0 = pbuffer.data(idx_eri_0_hd + 94);

    auto g_yyyyy_zz_0 = pbuffer.data(idx_eri_0_hd + 95);

    #pragma omp simd aligned(g_yyy_xx_0, g_yyy_xx_1, g_yyy_xy_0, g_yyy_xy_1, g_yyy_xz_0, g_yyy_xz_1, g_yyy_yy_0, g_yyy_yy_1, g_yyy_yz_0, g_yyy_yz_1, g_yyy_zz_0, g_yyy_zz_1, g_yyyy_x_1, g_yyyy_xx_1, g_yyyy_xy_1, g_yyyy_xz_1, g_yyyy_y_1, g_yyyy_yy_1, g_yyyy_yz_1, g_yyyy_z_1, g_yyyy_zz_1, g_yyyyy_xx_0, g_yyyyy_xy_0, g_yyyyy_xz_0, g_yyyyy_yy_0, g_yyyyy_yz_0, g_yyyyy_zz_0, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyyy_xx_0[i] = 4.0 * g_yyy_xx_0[i] * fbe_0 - 4.0 * g_yyy_xx_1[i] * fz_be_0 + g_yyyy_xx_1[i] * pa_y[i];

        g_yyyyy_xy_0[i] = 4.0 * g_yyy_xy_0[i] * fbe_0 - 4.0 * g_yyy_xy_1[i] * fz_be_0 + g_yyyy_x_1[i] * fe_0 + g_yyyy_xy_1[i] * pa_y[i];

        g_yyyyy_xz_0[i] = 4.0 * g_yyy_xz_0[i] * fbe_0 - 4.0 * g_yyy_xz_1[i] * fz_be_0 + g_yyyy_xz_1[i] * pa_y[i];

        g_yyyyy_yy_0[i] = 4.0 * g_yyy_yy_0[i] * fbe_0 - 4.0 * g_yyy_yy_1[i] * fz_be_0 + 2.0 * g_yyyy_y_1[i] * fe_0 + g_yyyy_yy_1[i] * pa_y[i];

        g_yyyyy_yz_0[i] = 4.0 * g_yyy_yz_0[i] * fbe_0 - 4.0 * g_yyy_yz_1[i] * fz_be_0 + g_yyyy_z_1[i] * fe_0 + g_yyyy_yz_1[i] * pa_y[i];

        g_yyyyy_zz_0[i] = 4.0 * g_yyy_zz_0[i] * fbe_0 - 4.0 * g_yyy_zz_1[i] * fz_be_0 + g_yyyy_zz_1[i] * pa_y[i];
    }

    // Set up 96-102 components of targeted buffer : HD

    auto g_yyyyz_xx_0 = pbuffer.data(idx_eri_0_hd + 96);

    auto g_yyyyz_xy_0 = pbuffer.data(idx_eri_0_hd + 97);

    auto g_yyyyz_xz_0 = pbuffer.data(idx_eri_0_hd + 98);

    auto g_yyyyz_yy_0 = pbuffer.data(idx_eri_0_hd + 99);

    auto g_yyyyz_yz_0 = pbuffer.data(idx_eri_0_hd + 100);

    auto g_yyyyz_zz_0 = pbuffer.data(idx_eri_0_hd + 101);

    #pragma omp simd aligned(g_yyyy_x_1, g_yyyy_xx_1, g_yyyy_xy_1, g_yyyy_xz_1, g_yyyy_y_1, g_yyyy_yy_1, g_yyyy_yz_1, g_yyyy_z_1, g_yyyy_zz_1, g_yyyyz_xx_0, g_yyyyz_xy_0, g_yyyyz_xz_0, g_yyyyz_yy_0, g_yyyyz_yz_0, g_yyyyz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yyyyz_xx_0[i] = g_yyyy_xx_1[i] * pa_z[i];

        g_yyyyz_xy_0[i] = g_yyyy_xy_1[i] * pa_z[i];

        g_yyyyz_xz_0[i] = g_yyyy_x_1[i] * fe_0 + g_yyyy_xz_1[i] * pa_z[i];

        g_yyyyz_yy_0[i] = g_yyyy_yy_1[i] * pa_z[i];

        g_yyyyz_yz_0[i] = g_yyyy_y_1[i] * fe_0 + g_yyyy_yz_1[i] * pa_z[i];

        g_yyyyz_zz_0[i] = 2.0 * g_yyyy_z_1[i] * fe_0 + g_yyyy_zz_1[i] * pa_z[i];
    }

    // Set up 102-108 components of targeted buffer : HD

    auto g_yyyzz_xx_0 = pbuffer.data(idx_eri_0_hd + 102);

    auto g_yyyzz_xy_0 = pbuffer.data(idx_eri_0_hd + 103);

    auto g_yyyzz_xz_0 = pbuffer.data(idx_eri_0_hd + 104);

    auto g_yyyzz_yy_0 = pbuffer.data(idx_eri_0_hd + 105);

    auto g_yyyzz_yz_0 = pbuffer.data(idx_eri_0_hd + 106);

    auto g_yyyzz_zz_0 = pbuffer.data(idx_eri_0_hd + 107);

    #pragma omp simd aligned(g_yyy_xy_0, g_yyy_xy_1, g_yyy_yy_0, g_yyy_yy_1, g_yyyz_xy_1, g_yyyz_yy_1, g_yyyzz_xx_0, g_yyyzz_xy_0, g_yyyzz_xz_0, g_yyyzz_yy_0, g_yyyzz_yz_0, g_yyyzz_zz_0, g_yyzz_xx_1, g_yyzz_xz_1, g_yyzz_yz_1, g_yyzz_z_1, g_yyzz_zz_1, g_yzz_xx_0, g_yzz_xx_1, g_yzz_xz_0, g_yzz_xz_1, g_yzz_yz_0, g_yzz_yz_1, g_yzz_zz_0, g_yzz_zz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyyzz_xx_0[i] = 2.0 * g_yzz_xx_0[i] * fbe_0 - 2.0 * g_yzz_xx_1[i] * fz_be_0 + g_yyzz_xx_1[i] * pa_y[i];

        g_yyyzz_xy_0[i] = g_yyy_xy_0[i] * fbe_0 - g_yyy_xy_1[i] * fz_be_0 + g_yyyz_xy_1[i] * pa_z[i];

        g_yyyzz_xz_0[i] = 2.0 * g_yzz_xz_0[i] * fbe_0 - 2.0 * g_yzz_xz_1[i] * fz_be_0 + g_yyzz_xz_1[i] * pa_y[i];

        g_yyyzz_yy_0[i] = g_yyy_yy_0[i] * fbe_0 - g_yyy_yy_1[i] * fz_be_0 + g_yyyz_yy_1[i] * pa_z[i];

        g_yyyzz_yz_0[i] = 2.0 * g_yzz_yz_0[i] * fbe_0 - 2.0 * g_yzz_yz_1[i] * fz_be_0 + g_yyzz_z_1[i] * fe_0 + g_yyzz_yz_1[i] * pa_y[i];

        g_yyyzz_zz_0[i] = 2.0 * g_yzz_zz_0[i] * fbe_0 - 2.0 * g_yzz_zz_1[i] * fz_be_0 + g_yyzz_zz_1[i] * pa_y[i];
    }

    // Set up 108-114 components of targeted buffer : HD

    auto g_yyzzz_xx_0 = pbuffer.data(idx_eri_0_hd + 108);

    auto g_yyzzz_xy_0 = pbuffer.data(idx_eri_0_hd + 109);

    auto g_yyzzz_xz_0 = pbuffer.data(idx_eri_0_hd + 110);

    auto g_yyzzz_yy_0 = pbuffer.data(idx_eri_0_hd + 111);

    auto g_yyzzz_yz_0 = pbuffer.data(idx_eri_0_hd + 112);

    auto g_yyzzz_zz_0 = pbuffer.data(idx_eri_0_hd + 113);

    #pragma omp simd aligned(g_yyz_xy_0, g_yyz_xy_1, g_yyz_yy_0, g_yyz_yy_1, g_yyzz_xy_1, g_yyzz_yy_1, g_yyzzz_xx_0, g_yyzzz_xy_0, g_yyzzz_xz_0, g_yyzzz_yy_0, g_yyzzz_yz_0, g_yyzzz_zz_0, g_yzzz_xx_1, g_yzzz_xz_1, g_yzzz_yz_1, g_yzzz_z_1, g_yzzz_zz_1, g_zzz_xx_0, g_zzz_xx_1, g_zzz_xz_0, g_zzz_xz_1, g_zzz_yz_0, g_zzz_yz_1, g_zzz_zz_0, g_zzz_zz_1, pa_y, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_yyzzz_xx_0[i] = g_zzz_xx_0[i] * fbe_0 - g_zzz_xx_1[i] * fz_be_0 + g_yzzz_xx_1[i] * pa_y[i];

        g_yyzzz_xy_0[i] = 2.0 * g_yyz_xy_0[i] * fbe_0 - 2.0 * g_yyz_xy_1[i] * fz_be_0 + g_yyzz_xy_1[i] * pa_z[i];

        g_yyzzz_xz_0[i] = g_zzz_xz_0[i] * fbe_0 - g_zzz_xz_1[i] * fz_be_0 + g_yzzz_xz_1[i] * pa_y[i];

        g_yyzzz_yy_0[i] = 2.0 * g_yyz_yy_0[i] * fbe_0 - 2.0 * g_yyz_yy_1[i] * fz_be_0 + g_yyzz_yy_1[i] * pa_z[i];

        g_yyzzz_yz_0[i] = g_zzz_yz_0[i] * fbe_0 - g_zzz_yz_1[i] * fz_be_0 + g_yzzz_z_1[i] * fe_0 + g_yzzz_yz_1[i] * pa_y[i];

        g_yyzzz_zz_0[i] = g_zzz_zz_0[i] * fbe_0 - g_zzz_zz_1[i] * fz_be_0 + g_yzzz_zz_1[i] * pa_y[i];
    }

    // Set up 114-120 components of targeted buffer : HD

    auto g_yzzzz_xx_0 = pbuffer.data(idx_eri_0_hd + 114);

    auto g_yzzzz_xy_0 = pbuffer.data(idx_eri_0_hd + 115);

    auto g_yzzzz_xz_0 = pbuffer.data(idx_eri_0_hd + 116);

    auto g_yzzzz_yy_0 = pbuffer.data(idx_eri_0_hd + 117);

    auto g_yzzzz_yz_0 = pbuffer.data(idx_eri_0_hd + 118);

    auto g_yzzzz_zz_0 = pbuffer.data(idx_eri_0_hd + 119);

    #pragma omp simd aligned(g_yzzzz_xx_0, g_yzzzz_xy_0, g_yzzzz_xz_0, g_yzzzz_yy_0, g_yzzzz_yz_0, g_yzzzz_zz_0, g_zzzz_x_1, g_zzzz_xx_1, g_zzzz_xy_1, g_zzzz_xz_1, g_zzzz_y_1, g_zzzz_yy_1, g_zzzz_yz_1, g_zzzz_z_1, g_zzzz_zz_1, pa_y, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        g_yzzzz_xx_0[i] = g_zzzz_xx_1[i] * pa_y[i];

        g_yzzzz_xy_0[i] = g_zzzz_x_1[i] * fe_0 + g_zzzz_xy_1[i] * pa_y[i];

        g_yzzzz_xz_0[i] = g_zzzz_xz_1[i] * pa_y[i];

        g_yzzzz_yy_0[i] = 2.0 * g_zzzz_y_1[i] * fe_0 + g_zzzz_yy_1[i] * pa_y[i];

        g_yzzzz_yz_0[i] = g_zzzz_z_1[i] * fe_0 + g_zzzz_yz_1[i] * pa_y[i];

        g_yzzzz_zz_0[i] = g_zzzz_zz_1[i] * pa_y[i];
    }

    // Set up 120-126 components of targeted buffer : HD

    auto g_zzzzz_xx_0 = pbuffer.data(idx_eri_0_hd + 120);

    auto g_zzzzz_xy_0 = pbuffer.data(idx_eri_0_hd + 121);

    auto g_zzzzz_xz_0 = pbuffer.data(idx_eri_0_hd + 122);

    auto g_zzzzz_yy_0 = pbuffer.data(idx_eri_0_hd + 123);

    auto g_zzzzz_yz_0 = pbuffer.data(idx_eri_0_hd + 124);

    auto g_zzzzz_zz_0 = pbuffer.data(idx_eri_0_hd + 125);

    #pragma omp simd aligned(g_zzz_xx_0, g_zzz_xx_1, g_zzz_xy_0, g_zzz_xy_1, g_zzz_xz_0, g_zzz_xz_1, g_zzz_yy_0, g_zzz_yy_1, g_zzz_yz_0, g_zzz_yz_1, g_zzz_zz_0, g_zzz_zz_1, g_zzzz_x_1, g_zzzz_xx_1, g_zzzz_xy_1, g_zzzz_xz_1, g_zzzz_y_1, g_zzzz_yy_1, g_zzzz_yz_1, g_zzzz_z_1, g_zzzz_zz_1, g_zzzzz_xx_0, g_zzzzz_xy_0, g_zzzzz_xz_0, g_zzzzz_yy_0, g_zzzzz_yz_0, g_zzzzz_zz_0, pa_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        const double fbe_0 = 0.5 / a_exp;

        const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;

        g_zzzzz_xx_0[i] = 4.0 * g_zzz_xx_0[i] * fbe_0 - 4.0 * g_zzz_xx_1[i] * fz_be_0 + g_zzzz_xx_1[i] * pa_z[i];

        g_zzzzz_xy_0[i] = 4.0 * g_zzz_xy_0[i] * fbe_0 - 4.0 * g_zzz_xy_1[i] * fz_be_0 + g_zzzz_xy_1[i] * pa_z[i];

        g_zzzzz_xz_0[i] = 4.0 * g_zzz_xz_0[i] * fbe_0 - 4.0 * g_zzz_xz_1[i] * fz_be_0 + g_zzzz_x_1[i] * fe_0 + g_zzzz_xz_1[i] * pa_z[i];

        g_zzzzz_yy_0[i] = 4.0 * g_zzz_yy_0[i] * fbe_0 - 4.0 * g_zzz_yy_1[i] * fz_be_0 + g_zzzz_yy_1[i] * pa_z[i];

        g_zzzzz_yz_0[i] = 4.0 * g_zzz_yz_0[i] * fbe_0 - 4.0 * g_zzz_yz_1[i] * fz_be_0 + g_zzzz_y_1[i] * fe_0 + g_zzzz_yz_1[i] * pa_z[i];

        g_zzzzz_zz_0[i] = 4.0 * g_zzz_zz_0[i] * fbe_0 - 4.0 * g_zzz_zz_1[i] * fz_be_0 + 2.0 * g_zzzz_z_1[i] * fe_0 + g_zzzz_zz_1[i] * pa_z[i];
    }

}

} // t2ceri namespace

