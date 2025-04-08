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

#include "ThreeCenterElectronRepulsionPrimRecGSD.hpp"

namespace t3ceri { // t3ceri namespace

auto
comp_prim_electron_repulsion_gsd(CSimdArray<double>& pbuffer,
                                 const size_t idx_eri_0_gsd,
                                 size_t idx_eri_0_dsd,
                                 size_t idx_eri_1_dsd,
                                 size_t idx_eri_1_fsp,
                                 size_t idx_eri_1_fsd,
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

    /// Set up components of auxilary buffer : DSD

    auto g_xx_0_xx_0 = pbuffer.data(idx_eri_0_dsd);

    auto g_xx_0_xy_0 = pbuffer.data(idx_eri_0_dsd + 1);

    auto g_xx_0_xz_0 = pbuffer.data(idx_eri_0_dsd + 2);

    auto g_xx_0_yy_0 = pbuffer.data(idx_eri_0_dsd + 3);

    auto g_xx_0_yz_0 = pbuffer.data(idx_eri_0_dsd + 4);

    auto g_xx_0_zz_0 = pbuffer.data(idx_eri_0_dsd + 5);

    auto g_yy_0_xx_0 = pbuffer.data(idx_eri_0_dsd + 18);

    auto g_yy_0_xy_0 = pbuffer.data(idx_eri_0_dsd + 19);

    auto g_yy_0_xz_0 = pbuffer.data(idx_eri_0_dsd + 20);

    auto g_yy_0_yy_0 = pbuffer.data(idx_eri_0_dsd + 21);

    auto g_yy_0_yz_0 = pbuffer.data(idx_eri_0_dsd + 22);

    auto g_yy_0_zz_0 = pbuffer.data(idx_eri_0_dsd + 23);

    auto g_zz_0_xx_0 = pbuffer.data(idx_eri_0_dsd + 30);

    auto g_zz_0_xy_0 = pbuffer.data(idx_eri_0_dsd + 31);

    auto g_zz_0_xz_0 = pbuffer.data(idx_eri_0_dsd + 32);

    auto g_zz_0_yy_0 = pbuffer.data(idx_eri_0_dsd + 33);

    auto g_zz_0_yz_0 = pbuffer.data(idx_eri_0_dsd + 34);

    auto g_zz_0_zz_0 = pbuffer.data(idx_eri_0_dsd + 35);

    /// Set up components of auxilary buffer : DSD

    auto g_xx_0_xx_1 = pbuffer.data(idx_eri_1_dsd);

    auto g_xx_0_xy_1 = pbuffer.data(idx_eri_1_dsd + 1);

    auto g_xx_0_xz_1 = pbuffer.data(idx_eri_1_dsd + 2);

    auto g_xx_0_yy_1 = pbuffer.data(idx_eri_1_dsd + 3);

    auto g_xx_0_yz_1 = pbuffer.data(idx_eri_1_dsd + 4);

    auto g_xx_0_zz_1 = pbuffer.data(idx_eri_1_dsd + 5);

    auto g_yy_0_xx_1 = pbuffer.data(idx_eri_1_dsd + 18);

    auto g_yy_0_xy_1 = pbuffer.data(idx_eri_1_dsd + 19);

    auto g_yy_0_xz_1 = pbuffer.data(idx_eri_1_dsd + 20);

    auto g_yy_0_yy_1 = pbuffer.data(idx_eri_1_dsd + 21);

    auto g_yy_0_yz_1 = pbuffer.data(idx_eri_1_dsd + 22);

    auto g_yy_0_zz_1 = pbuffer.data(idx_eri_1_dsd + 23);

    auto g_zz_0_xx_1 = pbuffer.data(idx_eri_1_dsd + 30);

    auto g_zz_0_xy_1 = pbuffer.data(idx_eri_1_dsd + 31);

    auto g_zz_0_xz_1 = pbuffer.data(idx_eri_1_dsd + 32);

    auto g_zz_0_yy_1 = pbuffer.data(idx_eri_1_dsd + 33);

    auto g_zz_0_yz_1 = pbuffer.data(idx_eri_1_dsd + 34);

    auto g_zz_0_zz_1 = pbuffer.data(idx_eri_1_dsd + 35);

    /// Set up components of auxilary buffer : FSP

    auto g_xxx_0_x_1 = pbuffer.data(idx_eri_1_fsp);

    auto g_xxx_0_y_1 = pbuffer.data(idx_eri_1_fsp + 1);

    auto g_xxx_0_z_1 = pbuffer.data(idx_eri_1_fsp + 2);

    auto g_xxz_0_z_1 = pbuffer.data(idx_eri_1_fsp + 8);

    auto g_xyy_0_y_1 = pbuffer.data(idx_eri_1_fsp + 10);

    auto g_xzz_0_z_1 = pbuffer.data(idx_eri_1_fsp + 17);

    auto g_yyy_0_x_1 = pbuffer.data(idx_eri_1_fsp + 18);

    auto g_yyy_0_y_1 = pbuffer.data(idx_eri_1_fsp + 19);

    auto g_yyy_0_z_1 = pbuffer.data(idx_eri_1_fsp + 20);

    auto g_yyz_0_z_1 = pbuffer.data(idx_eri_1_fsp + 23);

    auto g_yzz_0_y_1 = pbuffer.data(idx_eri_1_fsp + 25);

    auto g_yzz_0_z_1 = pbuffer.data(idx_eri_1_fsp + 26);

    auto g_zzz_0_x_1 = pbuffer.data(idx_eri_1_fsp + 27);

    auto g_zzz_0_y_1 = pbuffer.data(idx_eri_1_fsp + 28);

    auto g_zzz_0_z_1 = pbuffer.data(idx_eri_1_fsp + 29);

    /// Set up components of auxilary buffer : FSD

    auto g_xxx_0_xx_1 = pbuffer.data(idx_eri_1_fsd);

    auto g_xxx_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 1);

    auto g_xxx_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 2);

    auto g_xxx_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 3);

    auto g_xxx_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 4);

    auto g_xxx_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 5);

    auto g_xxy_0_xx_1 = pbuffer.data(idx_eri_1_fsd + 6);

    auto g_xxy_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 7);

    auto g_xxy_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 8);

    auto g_xxy_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 9);

    auto g_xxz_0_xx_1 = pbuffer.data(idx_eri_1_fsd + 12);

    auto g_xxz_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 13);

    auto g_xxz_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 14);

    auto g_xxz_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 16);

    auto g_xxz_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 17);

    auto g_xyy_0_xx_1 = pbuffer.data(idx_eri_1_fsd + 18);

    auto g_xyy_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 19);

    auto g_xyy_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 21);

    auto g_xyy_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 22);

    auto g_xyy_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 23);

    auto g_xzz_0_xx_1 = pbuffer.data(idx_eri_1_fsd + 30);

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

    auto g_yyz_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 44);

    auto g_yyz_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 45);

    auto g_yyz_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 46);

    auto g_yyz_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 47);

    auto g_yzz_0_xx_1 = pbuffer.data(idx_eri_1_fsd + 48);

    auto g_yzz_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 49);

    auto g_yzz_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 50);

    auto g_yzz_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 51);

    auto g_yzz_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 52);

    auto g_yzz_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 53);

    auto g_zzz_0_xx_1 = pbuffer.data(idx_eri_1_fsd + 54);

    auto g_zzz_0_xy_1 = pbuffer.data(idx_eri_1_fsd + 55);

    auto g_zzz_0_xz_1 = pbuffer.data(idx_eri_1_fsd + 56);

    auto g_zzz_0_yy_1 = pbuffer.data(idx_eri_1_fsd + 57);

    auto g_zzz_0_yz_1 = pbuffer.data(idx_eri_1_fsd + 58);

    auto g_zzz_0_zz_1 = pbuffer.data(idx_eri_1_fsd + 59);

    /// Set up 0-6 components of targeted buffer : GSD

    auto g_xxxx_0_xx_0 = pbuffer.data(idx_eri_0_gsd);

    auto g_xxxx_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 1);

    auto g_xxxx_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 2);

    auto g_xxxx_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 3);

    auto g_xxxx_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 4);

    auto g_xxxx_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 5);

    #pragma omp simd aligned(g_xx_0_xx_0, g_xx_0_xx_1, g_xx_0_xy_0, g_xx_0_xy_1, g_xx_0_xz_0, g_xx_0_xz_1, g_xx_0_yy_0, g_xx_0_yy_1, g_xx_0_yz_0, g_xx_0_yz_1, g_xx_0_zz_0, g_xx_0_zz_1, g_xxx_0_x_1, g_xxx_0_xx_1, g_xxx_0_xy_1, g_xxx_0_xz_1, g_xxx_0_y_1, g_xxx_0_yy_1, g_xxx_0_yz_1, g_xxx_0_z_1, g_xxx_0_zz_1, g_xxxx_0_xx_0, g_xxxx_0_xy_0, g_xxxx_0_xz_0, g_xxxx_0_yy_0, g_xxxx_0_yz_0, g_xxxx_0_zz_0, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxxx_0_xx_0[i] = 3.0 * g_xx_0_xx_0[i] * fbe_0 - 3.0 * g_xx_0_xx_1[i] * fz_be_0 + 2.0 * g_xxx_0_x_1[i] * fi_acd_0 + g_xxx_0_xx_1[i] * wa_x[i];

        g_xxxx_0_xy_0[i] = 3.0 * g_xx_0_xy_0[i] * fbe_0 - 3.0 * g_xx_0_xy_1[i] * fz_be_0 + g_xxx_0_y_1[i] * fi_acd_0 + g_xxx_0_xy_1[i] * wa_x[i];

        g_xxxx_0_xz_0[i] = 3.0 * g_xx_0_xz_0[i] * fbe_0 - 3.0 * g_xx_0_xz_1[i] * fz_be_0 + g_xxx_0_z_1[i] * fi_acd_0 + g_xxx_0_xz_1[i] * wa_x[i];

        g_xxxx_0_yy_0[i] = 3.0 * g_xx_0_yy_0[i] * fbe_0 - 3.0 * g_xx_0_yy_1[i] * fz_be_0 + g_xxx_0_yy_1[i] * wa_x[i];

        g_xxxx_0_yz_0[i] = 3.0 * g_xx_0_yz_0[i] * fbe_0 - 3.0 * g_xx_0_yz_1[i] * fz_be_0 + g_xxx_0_yz_1[i] * wa_x[i];

        g_xxxx_0_zz_0[i] = 3.0 * g_xx_0_zz_0[i] * fbe_0 - 3.0 * g_xx_0_zz_1[i] * fz_be_0 + g_xxx_0_zz_1[i] * wa_x[i];
    }

    /// Set up 6-12 components of targeted buffer : GSD

    auto g_xxxy_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 6);

    auto g_xxxy_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 7);

    auto g_xxxy_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 8);

    auto g_xxxy_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 9);

    auto g_xxxy_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 10);

    auto g_xxxy_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 11);

    #pragma omp simd aligned(g_xxx_0_x_1, g_xxx_0_xx_1, g_xxx_0_xy_1, g_xxx_0_xz_1, g_xxx_0_y_1, g_xxx_0_yy_1, g_xxx_0_yz_1, g_xxx_0_z_1, g_xxx_0_zz_1, g_xxxy_0_xx_0, g_xxxy_0_xy_0, g_xxxy_0_xz_0, g_xxxy_0_yy_0, g_xxxy_0_yz_0, g_xxxy_0_zz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxy_0_xx_0[i] = g_xxx_0_xx_1[i] * wa_y[i];

        g_xxxy_0_xy_0[i] = g_xxx_0_x_1[i] * fi_acd_0 + g_xxx_0_xy_1[i] * wa_y[i];

        g_xxxy_0_xz_0[i] = g_xxx_0_xz_1[i] * wa_y[i];

        g_xxxy_0_yy_0[i] = 2.0 * g_xxx_0_y_1[i] * fi_acd_0 + g_xxx_0_yy_1[i] * wa_y[i];

        g_xxxy_0_yz_0[i] = g_xxx_0_z_1[i] * fi_acd_0 + g_xxx_0_yz_1[i] * wa_y[i];

        g_xxxy_0_zz_0[i] = g_xxx_0_zz_1[i] * wa_y[i];
    }

    /// Set up 12-18 components of targeted buffer : GSD

    auto g_xxxz_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 12);

    auto g_xxxz_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 13);

    auto g_xxxz_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 14);

    auto g_xxxz_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 15);

    auto g_xxxz_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 16);

    auto g_xxxz_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 17);

    #pragma omp simd aligned(g_xxx_0_x_1, g_xxx_0_xx_1, g_xxx_0_xy_1, g_xxx_0_xz_1, g_xxx_0_y_1, g_xxx_0_yy_1, g_xxx_0_yz_1, g_xxx_0_z_1, g_xxx_0_zz_1, g_xxxz_0_xx_0, g_xxxz_0_xy_0, g_xxxz_0_xz_0, g_xxxz_0_yy_0, g_xxxz_0_yz_0, g_xxxz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxxz_0_xx_0[i] = g_xxx_0_xx_1[i] * wa_z[i];

        g_xxxz_0_xy_0[i] = g_xxx_0_xy_1[i] * wa_z[i];

        g_xxxz_0_xz_0[i] = g_xxx_0_x_1[i] * fi_acd_0 + g_xxx_0_xz_1[i] * wa_z[i];

        g_xxxz_0_yy_0[i] = g_xxx_0_yy_1[i] * wa_z[i];

        g_xxxz_0_yz_0[i] = g_xxx_0_y_1[i] * fi_acd_0 + g_xxx_0_yz_1[i] * wa_z[i];

        g_xxxz_0_zz_0[i] = 2.0 * g_xxx_0_z_1[i] * fi_acd_0 + g_xxx_0_zz_1[i] * wa_z[i];
    }

    /// Set up 18-24 components of targeted buffer : GSD

    auto g_xxyy_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 18);

    auto g_xxyy_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 19);

    auto g_xxyy_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 20);

    auto g_xxyy_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 21);

    auto g_xxyy_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 22);

    auto g_xxyy_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 23);

    #pragma omp simd aligned(g_xx_0_xx_0, g_xx_0_xx_1, g_xx_0_xz_0, g_xx_0_xz_1, g_xxy_0_xx_1, g_xxy_0_xz_1, g_xxyy_0_xx_0, g_xxyy_0_xy_0, g_xxyy_0_xz_0, g_xxyy_0_yy_0, g_xxyy_0_yz_0, g_xxyy_0_zz_0, g_xyy_0_xy_1, g_xyy_0_y_1, g_xyy_0_yy_1, g_xyy_0_yz_1, g_xyy_0_zz_1, g_yy_0_xy_0, g_yy_0_xy_1, g_yy_0_yy_0, g_yy_0_yy_1, g_yy_0_yz_0, g_yy_0_yz_1, g_yy_0_zz_0, g_yy_0_zz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxyy_0_xx_0[i] = g_xx_0_xx_0[i] * fbe_0 - g_xx_0_xx_1[i] * fz_be_0 + g_xxy_0_xx_1[i] * wa_y[i];

        g_xxyy_0_xy_0[i] = g_yy_0_xy_0[i] * fbe_0 - g_yy_0_xy_1[i] * fz_be_0 + g_xyy_0_y_1[i] * fi_acd_0 + g_xyy_0_xy_1[i] * wa_x[i];

        g_xxyy_0_xz_0[i] = g_xx_0_xz_0[i] * fbe_0 - g_xx_0_xz_1[i] * fz_be_0 + g_xxy_0_xz_1[i] * wa_y[i];

        g_xxyy_0_yy_0[i] = g_yy_0_yy_0[i] * fbe_0 - g_yy_0_yy_1[i] * fz_be_0 + g_xyy_0_yy_1[i] * wa_x[i];

        g_xxyy_0_yz_0[i] = g_yy_0_yz_0[i] * fbe_0 - g_yy_0_yz_1[i] * fz_be_0 + g_xyy_0_yz_1[i] * wa_x[i];

        g_xxyy_0_zz_0[i] = g_yy_0_zz_0[i] * fbe_0 - g_yy_0_zz_1[i] * fz_be_0 + g_xyy_0_zz_1[i] * wa_x[i];
    }

    /// Set up 24-30 components of targeted buffer : GSD

    auto g_xxyz_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 24);

    auto g_xxyz_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 25);

    auto g_xxyz_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 26);

    auto g_xxyz_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 27);

    auto g_xxyz_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 28);

    auto g_xxyz_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 29);

    #pragma omp simd aligned(g_xxy_0_xy_1, g_xxy_0_yy_1, g_xxyz_0_xx_0, g_xxyz_0_xy_0, g_xxyz_0_xz_0, g_xxyz_0_yy_0, g_xxyz_0_yz_0, g_xxyz_0_zz_0, g_xxz_0_xx_1, g_xxz_0_xz_1, g_xxz_0_yz_1, g_xxz_0_z_1, g_xxz_0_zz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xxyz_0_xx_0[i] = g_xxz_0_xx_1[i] * wa_y[i];

        g_xxyz_0_xy_0[i] = g_xxy_0_xy_1[i] * wa_z[i];

        g_xxyz_0_xz_0[i] = g_xxz_0_xz_1[i] * wa_y[i];

        g_xxyz_0_yy_0[i] = g_xxy_0_yy_1[i] * wa_z[i];

        g_xxyz_0_yz_0[i] = g_xxz_0_z_1[i] * fi_acd_0 + g_xxz_0_yz_1[i] * wa_y[i];

        g_xxyz_0_zz_0[i] = g_xxz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 30-36 components of targeted buffer : GSD

    auto g_xxzz_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 30);

    auto g_xxzz_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 31);

    auto g_xxzz_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 32);

    auto g_xxzz_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 33);

    auto g_xxzz_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 34);

    auto g_xxzz_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 35);

    #pragma omp simd aligned(g_xx_0_xx_0, g_xx_0_xx_1, g_xx_0_xy_0, g_xx_0_xy_1, g_xxz_0_xx_1, g_xxz_0_xy_1, g_xxzz_0_xx_0, g_xxzz_0_xy_0, g_xxzz_0_xz_0, g_xxzz_0_yy_0, g_xxzz_0_yz_0, g_xxzz_0_zz_0, g_xzz_0_xz_1, g_xzz_0_yy_1, g_xzz_0_yz_1, g_xzz_0_z_1, g_xzz_0_zz_1, g_zz_0_xz_0, g_zz_0_xz_1, g_zz_0_yy_0, g_zz_0_yy_1, g_zz_0_yz_0, g_zz_0_yz_1, g_zz_0_zz_0, g_zz_0_zz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_xxzz_0_xx_0[i] = g_xx_0_xx_0[i] * fbe_0 - g_xx_0_xx_1[i] * fz_be_0 + g_xxz_0_xx_1[i] * wa_z[i];

        g_xxzz_0_xy_0[i] = g_xx_0_xy_0[i] * fbe_0 - g_xx_0_xy_1[i] * fz_be_0 + g_xxz_0_xy_1[i] * wa_z[i];

        g_xxzz_0_xz_0[i] = g_zz_0_xz_0[i] * fbe_0 - g_zz_0_xz_1[i] * fz_be_0 + g_xzz_0_z_1[i] * fi_acd_0 + g_xzz_0_xz_1[i] * wa_x[i];

        g_xxzz_0_yy_0[i] = g_zz_0_yy_0[i] * fbe_0 - g_zz_0_yy_1[i] * fz_be_0 + g_xzz_0_yy_1[i] * wa_x[i];

        g_xxzz_0_yz_0[i] = g_zz_0_yz_0[i] * fbe_0 - g_zz_0_yz_1[i] * fz_be_0 + g_xzz_0_yz_1[i] * wa_x[i];

        g_xxzz_0_zz_0[i] = g_zz_0_zz_0[i] * fbe_0 - g_zz_0_zz_1[i] * fz_be_0 + g_xzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 36-42 components of targeted buffer : GSD

    auto g_xyyy_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 36);

    auto g_xyyy_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 37);

    auto g_xyyy_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 38);

    auto g_xyyy_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 39);

    auto g_xyyy_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 40);

    auto g_xyyy_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 41);

    #pragma omp simd aligned(g_xyyy_0_xx_0, g_xyyy_0_xy_0, g_xyyy_0_xz_0, g_xyyy_0_yy_0, g_xyyy_0_yz_0, g_xyyy_0_zz_0, g_yyy_0_x_1, g_yyy_0_xx_1, g_yyy_0_xy_1, g_yyy_0_xz_1, g_yyy_0_y_1, g_yyy_0_yy_1, g_yyy_0_yz_1, g_yyy_0_z_1, g_yyy_0_zz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyy_0_xx_0[i] = 2.0 * g_yyy_0_x_1[i] * fi_acd_0 + g_yyy_0_xx_1[i] * wa_x[i];

        g_xyyy_0_xy_0[i] = g_yyy_0_y_1[i] * fi_acd_0 + g_yyy_0_xy_1[i] * wa_x[i];

        g_xyyy_0_xz_0[i] = g_yyy_0_z_1[i] * fi_acd_0 + g_yyy_0_xz_1[i] * wa_x[i];

        g_xyyy_0_yy_0[i] = g_yyy_0_yy_1[i] * wa_x[i];

        g_xyyy_0_yz_0[i] = g_yyy_0_yz_1[i] * wa_x[i];

        g_xyyy_0_zz_0[i] = g_yyy_0_zz_1[i] * wa_x[i];
    }

    /// Set up 42-48 components of targeted buffer : GSD

    auto g_xyyz_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 42);

    auto g_xyyz_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 43);

    auto g_xyyz_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 44);

    auto g_xyyz_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 45);

    auto g_xyyz_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 46);

    auto g_xyyz_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 47);

    #pragma omp simd aligned(g_xyy_0_xx_1, g_xyy_0_xy_1, g_xyyz_0_xx_0, g_xyyz_0_xy_0, g_xyyz_0_xz_0, g_xyyz_0_yy_0, g_xyyz_0_yz_0, g_xyyz_0_zz_0, g_yyz_0_xz_1, g_yyz_0_yy_1, g_yyz_0_yz_1, g_yyz_0_z_1, g_yyz_0_zz_1, wa_x, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyyz_0_xx_0[i] = g_xyy_0_xx_1[i] * wa_z[i];

        g_xyyz_0_xy_0[i] = g_xyy_0_xy_1[i] * wa_z[i];

        g_xyyz_0_xz_0[i] = g_yyz_0_z_1[i] * fi_acd_0 + g_yyz_0_xz_1[i] * wa_x[i];

        g_xyyz_0_yy_0[i] = g_yyz_0_yy_1[i] * wa_x[i];

        g_xyyz_0_yz_0[i] = g_yyz_0_yz_1[i] * wa_x[i];

        g_xyyz_0_zz_0[i] = g_yyz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 48-54 components of targeted buffer : GSD

    auto g_xyzz_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 48);

    auto g_xyzz_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 49);

    auto g_xyzz_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 50);

    auto g_xyzz_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 51);

    auto g_xyzz_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 52);

    auto g_xyzz_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 53);

    #pragma omp simd aligned(g_xyzz_0_xx_0, g_xyzz_0_xy_0, g_xyzz_0_xz_0, g_xyzz_0_yy_0, g_xyzz_0_yz_0, g_xyzz_0_zz_0, g_xzz_0_xx_1, g_xzz_0_xz_1, g_yzz_0_xy_1, g_yzz_0_y_1, g_yzz_0_yy_1, g_yzz_0_yz_1, g_yzz_0_zz_1, wa_x, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xyzz_0_xx_0[i] = g_xzz_0_xx_1[i] * wa_y[i];

        g_xyzz_0_xy_0[i] = g_yzz_0_y_1[i] * fi_acd_0 + g_yzz_0_xy_1[i] * wa_x[i];

        g_xyzz_0_xz_0[i] = g_xzz_0_xz_1[i] * wa_y[i];

        g_xyzz_0_yy_0[i] = g_yzz_0_yy_1[i] * wa_x[i];

        g_xyzz_0_yz_0[i] = g_yzz_0_yz_1[i] * wa_x[i];

        g_xyzz_0_zz_0[i] = g_yzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 54-60 components of targeted buffer : GSD

    auto g_xzzz_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 54);

    auto g_xzzz_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 55);

    auto g_xzzz_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 56);

    auto g_xzzz_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 57);

    auto g_xzzz_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 58);

    auto g_xzzz_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 59);

    #pragma omp simd aligned(g_xzzz_0_xx_0, g_xzzz_0_xy_0, g_xzzz_0_xz_0, g_xzzz_0_yy_0, g_xzzz_0_yz_0, g_xzzz_0_zz_0, g_zzz_0_x_1, g_zzz_0_xx_1, g_zzz_0_xy_1, g_zzz_0_xz_1, g_zzz_0_y_1, g_zzz_0_yy_1, g_zzz_0_yz_1, g_zzz_0_z_1, g_zzz_0_zz_1, wa_x, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_xzzz_0_xx_0[i] = 2.0 * g_zzz_0_x_1[i] * fi_acd_0 + g_zzz_0_xx_1[i] * wa_x[i];

        g_xzzz_0_xy_0[i] = g_zzz_0_y_1[i] * fi_acd_0 + g_zzz_0_xy_1[i] * wa_x[i];

        g_xzzz_0_xz_0[i] = g_zzz_0_z_1[i] * fi_acd_0 + g_zzz_0_xz_1[i] * wa_x[i];

        g_xzzz_0_yy_0[i] = g_zzz_0_yy_1[i] * wa_x[i];

        g_xzzz_0_yz_0[i] = g_zzz_0_yz_1[i] * wa_x[i];

        g_xzzz_0_zz_0[i] = g_zzz_0_zz_1[i] * wa_x[i];
    }

    /// Set up 60-66 components of targeted buffer : GSD

    auto g_yyyy_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 60);

    auto g_yyyy_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 61);

    auto g_yyyy_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 62);

    auto g_yyyy_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 63);

    auto g_yyyy_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 64);

    auto g_yyyy_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 65);

    #pragma omp simd aligned(g_yy_0_xx_0, g_yy_0_xx_1, g_yy_0_xy_0, g_yy_0_xy_1, g_yy_0_xz_0, g_yy_0_xz_1, g_yy_0_yy_0, g_yy_0_yy_1, g_yy_0_yz_0, g_yy_0_yz_1, g_yy_0_zz_0, g_yy_0_zz_1, g_yyy_0_x_1, g_yyy_0_xx_1, g_yyy_0_xy_1, g_yyy_0_xz_1, g_yyy_0_y_1, g_yyy_0_yy_1, g_yyy_0_yz_1, g_yyy_0_z_1, g_yyy_0_zz_1, g_yyyy_0_xx_0, g_yyyy_0_xy_0, g_yyyy_0_xz_0, g_yyyy_0_yy_0, g_yyyy_0_yz_0, g_yyyy_0_zz_0, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyyy_0_xx_0[i] = 3.0 * g_yy_0_xx_0[i] * fbe_0 - 3.0 * g_yy_0_xx_1[i] * fz_be_0 + g_yyy_0_xx_1[i] * wa_y[i];

        g_yyyy_0_xy_0[i] = 3.0 * g_yy_0_xy_0[i] * fbe_0 - 3.0 * g_yy_0_xy_1[i] * fz_be_0 + g_yyy_0_x_1[i] * fi_acd_0 + g_yyy_0_xy_1[i] * wa_y[i];

        g_yyyy_0_xz_0[i] = 3.0 * g_yy_0_xz_0[i] * fbe_0 - 3.0 * g_yy_0_xz_1[i] * fz_be_0 + g_yyy_0_xz_1[i] * wa_y[i];

        g_yyyy_0_yy_0[i] = 3.0 * g_yy_0_yy_0[i] * fbe_0 - 3.0 * g_yy_0_yy_1[i] * fz_be_0 + 2.0 * g_yyy_0_y_1[i] * fi_acd_0 + g_yyy_0_yy_1[i] * wa_y[i];

        g_yyyy_0_yz_0[i] = 3.0 * g_yy_0_yz_0[i] * fbe_0 - 3.0 * g_yy_0_yz_1[i] * fz_be_0 + g_yyy_0_z_1[i] * fi_acd_0 + g_yyy_0_yz_1[i] * wa_y[i];

        g_yyyy_0_zz_0[i] = 3.0 * g_yy_0_zz_0[i] * fbe_0 - 3.0 * g_yy_0_zz_1[i] * fz_be_0 + g_yyy_0_zz_1[i] * wa_y[i];
    }

    /// Set up 66-72 components of targeted buffer : GSD

    auto g_yyyz_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 66);

    auto g_yyyz_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 67);

    auto g_yyyz_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 68);

    auto g_yyyz_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 69);

    auto g_yyyz_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 70);

    auto g_yyyz_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 71);

    #pragma omp simd aligned(g_yyy_0_x_1, g_yyy_0_xx_1, g_yyy_0_xy_1, g_yyy_0_xz_1, g_yyy_0_y_1, g_yyy_0_yy_1, g_yyy_0_yz_1, g_yyy_0_z_1, g_yyy_0_zz_1, g_yyyz_0_xx_0, g_yyyz_0_xy_0, g_yyyz_0_xz_0, g_yyyz_0_yy_0, g_yyyz_0_yz_0, g_yyyz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yyyz_0_xx_0[i] = g_yyy_0_xx_1[i] * wa_z[i];

        g_yyyz_0_xy_0[i] = g_yyy_0_xy_1[i] * wa_z[i];

        g_yyyz_0_xz_0[i] = g_yyy_0_x_1[i] * fi_acd_0 + g_yyy_0_xz_1[i] * wa_z[i];

        g_yyyz_0_yy_0[i] = g_yyy_0_yy_1[i] * wa_z[i];

        g_yyyz_0_yz_0[i] = g_yyy_0_y_1[i] * fi_acd_0 + g_yyy_0_yz_1[i] * wa_z[i];

        g_yyyz_0_zz_0[i] = 2.0 * g_yyy_0_z_1[i] * fi_acd_0 + g_yyy_0_zz_1[i] * wa_z[i];
    }

    /// Set up 72-78 components of targeted buffer : GSD

    auto g_yyzz_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 72);

    auto g_yyzz_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 73);

    auto g_yyzz_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 74);

    auto g_yyzz_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 75);

    auto g_yyzz_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 76);

    auto g_yyzz_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 77);

    #pragma omp simd aligned(g_yy_0_xy_0, g_yy_0_xy_1, g_yy_0_yy_0, g_yy_0_yy_1, g_yyz_0_xy_1, g_yyz_0_yy_1, g_yyzz_0_xx_0, g_yyzz_0_xy_0, g_yyzz_0_xz_0, g_yyzz_0_yy_0, g_yyzz_0_yz_0, g_yyzz_0_zz_0, g_yzz_0_xx_1, g_yzz_0_xz_1, g_yzz_0_yz_1, g_yzz_0_z_1, g_yzz_0_zz_1, g_zz_0_xx_0, g_zz_0_xx_1, g_zz_0_xz_0, g_zz_0_xz_1, g_zz_0_yz_0, g_zz_0_yz_1, g_zz_0_zz_0, g_zz_0_zz_1, wa_y, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_yyzz_0_xx_0[i] = g_zz_0_xx_0[i] * fbe_0 - g_zz_0_xx_1[i] * fz_be_0 + g_yzz_0_xx_1[i] * wa_y[i];

        g_yyzz_0_xy_0[i] = g_yy_0_xy_0[i] * fbe_0 - g_yy_0_xy_1[i] * fz_be_0 + g_yyz_0_xy_1[i] * wa_z[i];

        g_yyzz_0_xz_0[i] = g_zz_0_xz_0[i] * fbe_0 - g_zz_0_xz_1[i] * fz_be_0 + g_yzz_0_xz_1[i] * wa_y[i];

        g_yyzz_0_yy_0[i] = g_yy_0_yy_0[i] * fbe_0 - g_yy_0_yy_1[i] * fz_be_0 + g_yyz_0_yy_1[i] * wa_z[i];

        g_yyzz_0_yz_0[i] = g_zz_0_yz_0[i] * fbe_0 - g_zz_0_yz_1[i] * fz_be_0 + g_yzz_0_z_1[i] * fi_acd_0 + g_yzz_0_yz_1[i] * wa_y[i];

        g_yyzz_0_zz_0[i] = g_zz_0_zz_0[i] * fbe_0 - g_zz_0_zz_1[i] * fz_be_0 + g_yzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 78-84 components of targeted buffer : GSD

    auto g_yzzz_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 78);

    auto g_yzzz_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 79);

    auto g_yzzz_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 80);

    auto g_yzzz_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 81);

    auto g_yzzz_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 82);

    auto g_yzzz_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 83);

    #pragma omp simd aligned(g_yzzz_0_xx_0, g_yzzz_0_xy_0, g_yzzz_0_xz_0, g_yzzz_0_yy_0, g_yzzz_0_yz_0, g_yzzz_0_zz_0, g_zzz_0_x_1, g_zzz_0_xx_1, g_zzz_0_xy_1, g_zzz_0_xz_1, g_zzz_0_y_1, g_zzz_0_yy_1, g_zzz_0_yz_1, g_zzz_0_z_1, g_zzz_0_zz_1, wa_y, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        g_yzzz_0_xx_0[i] = g_zzz_0_xx_1[i] * wa_y[i];

        g_yzzz_0_xy_0[i] = g_zzz_0_x_1[i] * fi_acd_0 + g_zzz_0_xy_1[i] * wa_y[i];

        g_yzzz_0_xz_0[i] = g_zzz_0_xz_1[i] * wa_y[i];

        g_yzzz_0_yy_0[i] = 2.0 * g_zzz_0_y_1[i] * fi_acd_0 + g_zzz_0_yy_1[i] * wa_y[i];

        g_yzzz_0_yz_0[i] = g_zzz_0_z_1[i] * fi_acd_0 + g_zzz_0_yz_1[i] * wa_y[i];

        g_yzzz_0_zz_0[i] = g_zzz_0_zz_1[i] * wa_y[i];
    }

    /// Set up 84-90 components of targeted buffer : GSD

    auto g_zzzz_0_xx_0 = pbuffer.data(idx_eri_0_gsd + 84);

    auto g_zzzz_0_xy_0 = pbuffer.data(idx_eri_0_gsd + 85);

    auto g_zzzz_0_xz_0 = pbuffer.data(idx_eri_0_gsd + 86);

    auto g_zzzz_0_yy_0 = pbuffer.data(idx_eri_0_gsd + 87);

    auto g_zzzz_0_yz_0 = pbuffer.data(idx_eri_0_gsd + 88);

    auto g_zzzz_0_zz_0 = pbuffer.data(idx_eri_0_gsd + 89);

    #pragma omp simd aligned(g_zz_0_xx_0, g_zz_0_xx_1, g_zz_0_xy_0, g_zz_0_xy_1, g_zz_0_xz_0, g_zz_0_xz_1, g_zz_0_yy_0, g_zz_0_yy_1, g_zz_0_yz_0, g_zz_0_yz_1, g_zz_0_zz_0, g_zz_0_zz_1, g_zzz_0_x_1, g_zzz_0_xx_1, g_zzz_0_xy_1, g_zzz_0_xz_1, g_zzz_0_y_1, g_zzz_0_yy_1, g_zzz_0_yz_1, g_zzz_0_z_1, g_zzz_0_zz_1, g_zzzz_0_xx_0, g_zzzz_0_xy_0, g_zzzz_0_xz_0, g_zzzz_0_yy_0, g_zzzz_0_yz_0, g_zzzz_0_zz_0, wa_z, c_exps, d_exps  : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fbe_0 = 0.5 / a_exp;

        const double fi_acd_0 = 0.5 / (a_exp + c_exps[i] + d_exps[i]);

        const double fz_be_0 = 2.0 * (c_exps[i] + d_exps[i]) * fbe_0 * fi_acd_0;

        g_zzzz_0_xx_0[i] = 3.0 * g_zz_0_xx_0[i] * fbe_0 - 3.0 * g_zz_0_xx_1[i] * fz_be_0 + g_zzz_0_xx_1[i] * wa_z[i];

        g_zzzz_0_xy_0[i] = 3.0 * g_zz_0_xy_0[i] * fbe_0 - 3.0 * g_zz_0_xy_1[i] * fz_be_0 + g_zzz_0_xy_1[i] * wa_z[i];

        g_zzzz_0_xz_0[i] = 3.0 * g_zz_0_xz_0[i] * fbe_0 - 3.0 * g_zz_0_xz_1[i] * fz_be_0 + g_zzz_0_x_1[i] * fi_acd_0 + g_zzz_0_xz_1[i] * wa_z[i];

        g_zzzz_0_yy_0[i] = 3.0 * g_zz_0_yy_0[i] * fbe_0 - 3.0 * g_zz_0_yy_1[i] * fz_be_0 + g_zzz_0_yy_1[i] * wa_z[i];

        g_zzzz_0_yz_0[i] = 3.0 * g_zz_0_yz_0[i] * fbe_0 - 3.0 * g_zz_0_yz_1[i] * fz_be_0 + g_zzz_0_y_1[i] * fi_acd_0 + g_zzz_0_yz_1[i] * wa_z[i];

        g_zzzz_0_zz_0[i] = 3.0 * g_zz_0_zz_0[i] * fbe_0 - 3.0 * g_zz_0_zz_1[i] * fz_be_0 + 2.0 * g_zzz_0_z_1[i] * fi_acd_0 + g_zzz_0_zz_1[i] * wa_z[i];
    }
}

} // t3ceri namespace

