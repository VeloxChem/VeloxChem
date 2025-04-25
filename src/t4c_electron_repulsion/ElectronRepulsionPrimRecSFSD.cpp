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

#include "ElectronRepulsionPrimRecSFSD.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sfsd(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sfsd,
                                  size_t                idx_eri_0_spsd,
                                  size_t                idx_eri_1_spsd,
                                  size_t                idx_eri_1_sdsp,
                                  size_t                idx_eri_0_sdsd,
                                  size_t                idx_eri_1_sdsd,
                                  CSimdArray<double>&   factors,
                                  const size_t          idx_wp,
                                  const TPoint<double>& r_pb,
                                  const double          a_exp,
                                  const double          b_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto c_exps = factors.data(0);

    auto d_exps = factors.data(1);

    // Set up R(WP) distances

    auto wp_x = factors.data(idx_wp);

    auto wp_y = factors.data(idx_wp + 1);

    auto wp_z = factors.data(idx_wp + 2);

    // set up R(PB) distances

    const auto xyz = r_pb.coordinates();

    const auto pb_x = xyz[0];

    const auto pb_y = xyz[1];

    const auto pb_z = xyz[2];

    /// Set up components of auxilary buffer : SPSD

    auto g_0_x_0_xx_0 = pbuffer.data(idx_eri_0_spsd);

    auto g_0_x_0_xy_0 = pbuffer.data(idx_eri_0_spsd + 1);

    auto g_0_x_0_xz_0 = pbuffer.data(idx_eri_0_spsd + 2);

    auto g_0_x_0_yy_0 = pbuffer.data(idx_eri_0_spsd + 3);

    auto g_0_x_0_yz_0 = pbuffer.data(idx_eri_0_spsd + 4);

    auto g_0_x_0_zz_0 = pbuffer.data(idx_eri_0_spsd + 5);

    auto g_0_y_0_xx_0 = pbuffer.data(idx_eri_0_spsd + 6);

    auto g_0_y_0_xy_0 = pbuffer.data(idx_eri_0_spsd + 7);

    auto g_0_y_0_xz_0 = pbuffer.data(idx_eri_0_spsd + 8);

    auto g_0_y_0_yy_0 = pbuffer.data(idx_eri_0_spsd + 9);

    auto g_0_y_0_yz_0 = pbuffer.data(idx_eri_0_spsd + 10);

    auto g_0_y_0_zz_0 = pbuffer.data(idx_eri_0_spsd + 11);

    auto g_0_z_0_xx_0 = pbuffer.data(idx_eri_0_spsd + 12);

    auto g_0_z_0_xy_0 = pbuffer.data(idx_eri_0_spsd + 13);

    auto g_0_z_0_xz_0 = pbuffer.data(idx_eri_0_spsd + 14);

    auto g_0_z_0_yy_0 = pbuffer.data(idx_eri_0_spsd + 15);

    auto g_0_z_0_yz_0 = pbuffer.data(idx_eri_0_spsd + 16);

    auto g_0_z_0_zz_0 = pbuffer.data(idx_eri_0_spsd + 17);

    /// Set up components of auxilary buffer : SPSD

    auto g_0_x_0_xx_1 = pbuffer.data(idx_eri_1_spsd);

    auto g_0_x_0_xy_1 = pbuffer.data(idx_eri_1_spsd + 1);

    auto g_0_x_0_xz_1 = pbuffer.data(idx_eri_1_spsd + 2);

    auto g_0_x_0_yy_1 = pbuffer.data(idx_eri_1_spsd + 3);

    auto g_0_x_0_yz_1 = pbuffer.data(idx_eri_1_spsd + 4);

    auto g_0_x_0_zz_1 = pbuffer.data(idx_eri_1_spsd + 5);

    auto g_0_y_0_xx_1 = pbuffer.data(idx_eri_1_spsd + 6);

    auto g_0_y_0_xy_1 = pbuffer.data(idx_eri_1_spsd + 7);

    auto g_0_y_0_xz_1 = pbuffer.data(idx_eri_1_spsd + 8);

    auto g_0_y_0_yy_1 = pbuffer.data(idx_eri_1_spsd + 9);

    auto g_0_y_0_yz_1 = pbuffer.data(idx_eri_1_spsd + 10);

    auto g_0_y_0_zz_1 = pbuffer.data(idx_eri_1_spsd + 11);

    auto g_0_z_0_xx_1 = pbuffer.data(idx_eri_1_spsd + 12);

    auto g_0_z_0_xy_1 = pbuffer.data(idx_eri_1_spsd + 13);

    auto g_0_z_0_xz_1 = pbuffer.data(idx_eri_1_spsd + 14);

    auto g_0_z_0_yy_1 = pbuffer.data(idx_eri_1_spsd + 15);

    auto g_0_z_0_yz_1 = pbuffer.data(idx_eri_1_spsd + 16);

    auto g_0_z_0_zz_1 = pbuffer.data(idx_eri_1_spsd + 17);

    /// Set up components of auxilary buffer : SDSP

    auto g_0_xx_0_x_1 = pbuffer.data(idx_eri_1_sdsp);

    auto g_0_xx_0_y_1 = pbuffer.data(idx_eri_1_sdsp + 1);

    auto g_0_xx_0_z_1 = pbuffer.data(idx_eri_1_sdsp + 2);

    auto g_0_yy_0_x_1 = pbuffer.data(idx_eri_1_sdsp + 9);

    auto g_0_yy_0_y_1 = pbuffer.data(idx_eri_1_sdsp + 10);

    auto g_0_yy_0_z_1 = pbuffer.data(idx_eri_1_sdsp + 11);

    auto g_0_zz_0_x_1 = pbuffer.data(idx_eri_1_sdsp + 15);

    auto g_0_zz_0_y_1 = pbuffer.data(idx_eri_1_sdsp + 16);

    auto g_0_zz_0_z_1 = pbuffer.data(idx_eri_1_sdsp + 17);

    /// Set up components of auxilary buffer : SDSD

    auto g_0_xx_0_xx_0 = pbuffer.data(idx_eri_0_sdsd);

    auto g_0_xx_0_xy_0 = pbuffer.data(idx_eri_0_sdsd + 1);

    auto g_0_xx_0_xz_0 = pbuffer.data(idx_eri_0_sdsd + 2);

    auto g_0_xx_0_yy_0 = pbuffer.data(idx_eri_0_sdsd + 3);

    auto g_0_xx_0_yz_0 = pbuffer.data(idx_eri_0_sdsd + 4);

    auto g_0_xx_0_zz_0 = pbuffer.data(idx_eri_0_sdsd + 5);

    auto g_0_xy_0_xy_0 = pbuffer.data(idx_eri_0_sdsd + 7);

    auto g_0_xz_0_xx_0 = pbuffer.data(idx_eri_0_sdsd + 12);

    auto g_0_xz_0_xz_0 = pbuffer.data(idx_eri_0_sdsd + 14);

    auto g_0_yy_0_xx_0 = pbuffer.data(idx_eri_0_sdsd + 18);

    auto g_0_yy_0_xy_0 = pbuffer.data(idx_eri_0_sdsd + 19);

    auto g_0_yy_0_xz_0 = pbuffer.data(idx_eri_0_sdsd + 20);

    auto g_0_yy_0_yy_0 = pbuffer.data(idx_eri_0_sdsd + 21);

    auto g_0_yy_0_yz_0 = pbuffer.data(idx_eri_0_sdsd + 22);

    auto g_0_yy_0_zz_0 = pbuffer.data(idx_eri_0_sdsd + 23);

    auto g_0_yz_0_yy_0 = pbuffer.data(idx_eri_0_sdsd + 27);

    auto g_0_yz_0_yz_0 = pbuffer.data(idx_eri_0_sdsd + 28);

    auto g_0_yz_0_zz_0 = pbuffer.data(idx_eri_0_sdsd + 29);

    auto g_0_zz_0_xx_0 = pbuffer.data(idx_eri_0_sdsd + 30);

    auto g_0_zz_0_xy_0 = pbuffer.data(idx_eri_0_sdsd + 31);

    auto g_0_zz_0_xz_0 = pbuffer.data(idx_eri_0_sdsd + 32);

    auto g_0_zz_0_yy_0 = pbuffer.data(idx_eri_0_sdsd + 33);

    auto g_0_zz_0_yz_0 = pbuffer.data(idx_eri_0_sdsd + 34);

    auto g_0_zz_0_zz_0 = pbuffer.data(idx_eri_0_sdsd + 35);

    /// Set up components of auxilary buffer : SDSD

    auto g_0_xx_0_xx_1 = pbuffer.data(idx_eri_1_sdsd);

    auto g_0_xx_0_xy_1 = pbuffer.data(idx_eri_1_sdsd + 1);

    auto g_0_xx_0_xz_1 = pbuffer.data(idx_eri_1_sdsd + 2);

    auto g_0_xx_0_yy_1 = pbuffer.data(idx_eri_1_sdsd + 3);

    auto g_0_xx_0_yz_1 = pbuffer.data(idx_eri_1_sdsd + 4);

    auto g_0_xx_0_zz_1 = pbuffer.data(idx_eri_1_sdsd + 5);

    auto g_0_xy_0_xy_1 = pbuffer.data(idx_eri_1_sdsd + 7);

    auto g_0_xz_0_xx_1 = pbuffer.data(idx_eri_1_sdsd + 12);

    auto g_0_xz_0_xz_1 = pbuffer.data(idx_eri_1_sdsd + 14);

    auto g_0_yy_0_xx_1 = pbuffer.data(idx_eri_1_sdsd + 18);

    auto g_0_yy_0_xy_1 = pbuffer.data(idx_eri_1_sdsd + 19);

    auto g_0_yy_0_xz_1 = pbuffer.data(idx_eri_1_sdsd + 20);

    auto g_0_yy_0_yy_1 = pbuffer.data(idx_eri_1_sdsd + 21);

    auto g_0_yy_0_yz_1 = pbuffer.data(idx_eri_1_sdsd + 22);

    auto g_0_yy_0_zz_1 = pbuffer.data(idx_eri_1_sdsd + 23);

    auto g_0_yz_0_yy_1 = pbuffer.data(idx_eri_1_sdsd + 27);

    auto g_0_yz_0_yz_1 = pbuffer.data(idx_eri_1_sdsd + 28);

    auto g_0_yz_0_zz_1 = pbuffer.data(idx_eri_1_sdsd + 29);

    auto g_0_zz_0_xx_1 = pbuffer.data(idx_eri_1_sdsd + 30);

    auto g_0_zz_0_xy_1 = pbuffer.data(idx_eri_1_sdsd + 31);

    auto g_0_zz_0_xz_1 = pbuffer.data(idx_eri_1_sdsd + 32);

    auto g_0_zz_0_yy_1 = pbuffer.data(idx_eri_1_sdsd + 33);

    auto g_0_zz_0_yz_1 = pbuffer.data(idx_eri_1_sdsd + 34);

    auto g_0_zz_0_zz_1 = pbuffer.data(idx_eri_1_sdsd + 35);

    /// Set up 0-6 components of targeted buffer : SFSD

    auto g_0_xxx_0_xx_0 = pbuffer.data(idx_eri_0_sfsd);

    auto g_0_xxx_0_xy_0 = pbuffer.data(idx_eri_0_sfsd + 1);

    auto g_0_xxx_0_xz_0 = pbuffer.data(idx_eri_0_sfsd + 2);

    auto g_0_xxx_0_yy_0 = pbuffer.data(idx_eri_0_sfsd + 3);

    auto g_0_xxx_0_yz_0 = pbuffer.data(idx_eri_0_sfsd + 4);

    auto g_0_xxx_0_zz_0 = pbuffer.data(idx_eri_0_sfsd + 5);

#pragma omp simd aligned(g_0_x_0_xx_0,       \
                             g_0_x_0_xx_1,   \
                             g_0_x_0_xy_0,   \
                             g_0_x_0_xy_1,   \
                             g_0_x_0_xz_0,   \
                             g_0_x_0_xz_1,   \
                             g_0_x_0_yy_0,   \
                             g_0_x_0_yy_1,   \
                             g_0_x_0_yz_0,   \
                             g_0_x_0_yz_1,   \
                             g_0_x_0_zz_0,   \
                             g_0_x_0_zz_1,   \
                             g_0_xx_0_x_1,   \
                             g_0_xx_0_xx_0,  \
                             g_0_xx_0_xx_1,  \
                             g_0_xx_0_xy_0,  \
                             g_0_xx_0_xy_1,  \
                             g_0_xx_0_xz_0,  \
                             g_0_xx_0_xz_1,  \
                             g_0_xx_0_y_1,   \
                             g_0_xx_0_yy_0,  \
                             g_0_xx_0_yy_1,  \
                             g_0_xx_0_yz_0,  \
                             g_0_xx_0_yz_1,  \
                             g_0_xx_0_z_1,   \
                             g_0_xx_0_zz_0,  \
                             g_0_xx_0_zz_1,  \
                             g_0_xxx_0_xx_0, \
                             g_0_xxx_0_xy_0, \
                             g_0_xxx_0_xz_0, \
                             g_0_xxx_0_yy_0, \
                             g_0_xxx_0_yz_0, \
                             g_0_xxx_0_zz_0, \
                             wp_x,           \
                             c_exps,         \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxx_0_xx_0[i] = 2.0 * g_0_x_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xx_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_x_1[i] * fi_abcd_0 +
                            g_0_xx_0_xx_0[i] * pb_x + g_0_xx_0_xx_1[i] * wp_x[i];

        g_0_xxx_0_xy_0[i] = 2.0 * g_0_x_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xy_1[i] * fti_ab_0 + g_0_xx_0_y_1[i] * fi_abcd_0 +
                            g_0_xx_0_xy_0[i] * pb_x + g_0_xx_0_xy_1[i] * wp_x[i];

        g_0_xxx_0_xz_0[i] = 2.0 * g_0_x_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xz_1[i] * fti_ab_0 + g_0_xx_0_z_1[i] * fi_abcd_0 +
                            g_0_xx_0_xz_0[i] * pb_x + g_0_xx_0_xz_1[i] * wp_x[i];

        g_0_xxx_0_yy_0[i] = 2.0 * g_0_x_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yy_1[i] * fti_ab_0 + g_0_xx_0_yy_0[i] * pb_x + g_0_xx_0_yy_1[i] * wp_x[i];

        g_0_xxx_0_yz_0[i] = 2.0 * g_0_x_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yz_1[i] * fti_ab_0 + g_0_xx_0_yz_0[i] * pb_x + g_0_xx_0_yz_1[i] * wp_x[i];

        g_0_xxx_0_zz_0[i] = 2.0 * g_0_x_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_zz_1[i] * fti_ab_0 + g_0_xx_0_zz_0[i] * pb_x + g_0_xx_0_zz_1[i] * wp_x[i];
    }

    /// Set up 6-12 components of targeted buffer : SFSD

    auto g_0_xxy_0_xx_0 = pbuffer.data(idx_eri_0_sfsd + 6);

    auto g_0_xxy_0_xy_0 = pbuffer.data(idx_eri_0_sfsd + 7);

    auto g_0_xxy_0_xz_0 = pbuffer.data(idx_eri_0_sfsd + 8);

    auto g_0_xxy_0_yy_0 = pbuffer.data(idx_eri_0_sfsd + 9);

    auto g_0_xxy_0_yz_0 = pbuffer.data(idx_eri_0_sfsd + 10);

    auto g_0_xxy_0_zz_0 = pbuffer.data(idx_eri_0_sfsd + 11);

#pragma omp simd aligned(g_0_xx_0_x_1,       \
                             g_0_xx_0_xx_0,  \
                             g_0_xx_0_xx_1,  \
                             g_0_xx_0_xy_0,  \
                             g_0_xx_0_xy_1,  \
                             g_0_xx_0_xz_0,  \
                             g_0_xx_0_xz_1,  \
                             g_0_xx_0_y_1,   \
                             g_0_xx_0_yy_0,  \
                             g_0_xx_0_yy_1,  \
                             g_0_xx_0_yz_0,  \
                             g_0_xx_0_yz_1,  \
                             g_0_xx_0_z_1,   \
                             g_0_xx_0_zz_0,  \
                             g_0_xx_0_zz_1,  \
                             g_0_xxy_0_xx_0, \
                             g_0_xxy_0_xy_0, \
                             g_0_xxy_0_xz_0, \
                             g_0_xxy_0_yy_0, \
                             g_0_xxy_0_yz_0, \
                             g_0_xxy_0_zz_0, \
                             wp_y,           \
                             c_exps,         \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxy_0_xx_0[i] = g_0_xx_0_xx_0[i] * pb_y + g_0_xx_0_xx_1[i] * wp_y[i];

        g_0_xxy_0_xy_0[i] = g_0_xx_0_x_1[i] * fi_abcd_0 + g_0_xx_0_xy_0[i] * pb_y + g_0_xx_0_xy_1[i] * wp_y[i];

        g_0_xxy_0_xz_0[i] = g_0_xx_0_xz_0[i] * pb_y + g_0_xx_0_xz_1[i] * wp_y[i];

        g_0_xxy_0_yy_0[i] = 2.0 * g_0_xx_0_y_1[i] * fi_abcd_0 + g_0_xx_0_yy_0[i] * pb_y + g_0_xx_0_yy_1[i] * wp_y[i];

        g_0_xxy_0_yz_0[i] = g_0_xx_0_z_1[i] * fi_abcd_0 + g_0_xx_0_yz_0[i] * pb_y + g_0_xx_0_yz_1[i] * wp_y[i];

        g_0_xxy_0_zz_0[i] = g_0_xx_0_zz_0[i] * pb_y + g_0_xx_0_zz_1[i] * wp_y[i];
    }

    /// Set up 12-18 components of targeted buffer : SFSD

    auto g_0_xxz_0_xx_0 = pbuffer.data(idx_eri_0_sfsd + 12);

    auto g_0_xxz_0_xy_0 = pbuffer.data(idx_eri_0_sfsd + 13);

    auto g_0_xxz_0_xz_0 = pbuffer.data(idx_eri_0_sfsd + 14);

    auto g_0_xxz_0_yy_0 = pbuffer.data(idx_eri_0_sfsd + 15);

    auto g_0_xxz_0_yz_0 = pbuffer.data(idx_eri_0_sfsd + 16);

    auto g_0_xxz_0_zz_0 = pbuffer.data(idx_eri_0_sfsd + 17);

#pragma omp simd aligned(g_0_xx_0_x_1,       \
                             g_0_xx_0_xx_0,  \
                             g_0_xx_0_xx_1,  \
                             g_0_xx_0_xy_0,  \
                             g_0_xx_0_xy_1,  \
                             g_0_xx_0_xz_0,  \
                             g_0_xx_0_xz_1,  \
                             g_0_xx_0_y_1,   \
                             g_0_xx_0_yy_0,  \
                             g_0_xx_0_yy_1,  \
                             g_0_xx_0_yz_0,  \
                             g_0_xx_0_yz_1,  \
                             g_0_xx_0_z_1,   \
                             g_0_xx_0_zz_0,  \
                             g_0_xx_0_zz_1,  \
                             g_0_xxz_0_xx_0, \
                             g_0_xxz_0_xy_0, \
                             g_0_xxz_0_xz_0, \
                             g_0_xxz_0_yy_0, \
                             g_0_xxz_0_yz_0, \
                             g_0_xxz_0_zz_0, \
                             wp_z,           \
                             c_exps,         \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxz_0_xx_0[i] = g_0_xx_0_xx_0[i] * pb_z + g_0_xx_0_xx_1[i] * wp_z[i];

        g_0_xxz_0_xy_0[i] = g_0_xx_0_xy_0[i] * pb_z + g_0_xx_0_xy_1[i] * wp_z[i];

        g_0_xxz_0_xz_0[i] = g_0_xx_0_x_1[i] * fi_abcd_0 + g_0_xx_0_xz_0[i] * pb_z + g_0_xx_0_xz_1[i] * wp_z[i];

        g_0_xxz_0_yy_0[i] = g_0_xx_0_yy_0[i] * pb_z + g_0_xx_0_yy_1[i] * wp_z[i];

        g_0_xxz_0_yz_0[i] = g_0_xx_0_y_1[i] * fi_abcd_0 + g_0_xx_0_yz_0[i] * pb_z + g_0_xx_0_yz_1[i] * wp_z[i];

        g_0_xxz_0_zz_0[i] = 2.0 * g_0_xx_0_z_1[i] * fi_abcd_0 + g_0_xx_0_zz_0[i] * pb_z + g_0_xx_0_zz_1[i] * wp_z[i];
    }

    /// Set up 18-24 components of targeted buffer : SFSD

    auto g_0_xyy_0_xx_0 = pbuffer.data(idx_eri_0_sfsd + 18);

    auto g_0_xyy_0_xy_0 = pbuffer.data(idx_eri_0_sfsd + 19);

    auto g_0_xyy_0_xz_0 = pbuffer.data(idx_eri_0_sfsd + 20);

    auto g_0_xyy_0_yy_0 = pbuffer.data(idx_eri_0_sfsd + 21);

    auto g_0_xyy_0_yz_0 = pbuffer.data(idx_eri_0_sfsd + 22);

    auto g_0_xyy_0_zz_0 = pbuffer.data(idx_eri_0_sfsd + 23);

#pragma omp simd aligned(g_0_xyy_0_xx_0,     \
                             g_0_xyy_0_xy_0, \
                             g_0_xyy_0_xz_0, \
                             g_0_xyy_0_yy_0, \
                             g_0_xyy_0_yz_0, \
                             g_0_xyy_0_zz_0, \
                             g_0_yy_0_x_1,   \
                             g_0_yy_0_xx_0,  \
                             g_0_yy_0_xx_1,  \
                             g_0_yy_0_xy_0,  \
                             g_0_yy_0_xy_1,  \
                             g_0_yy_0_xz_0,  \
                             g_0_yy_0_xz_1,  \
                             g_0_yy_0_y_1,   \
                             g_0_yy_0_yy_0,  \
                             g_0_yy_0_yy_1,  \
                             g_0_yy_0_yz_0,  \
                             g_0_yy_0_yz_1,  \
                             g_0_yy_0_z_1,   \
                             g_0_yy_0_zz_0,  \
                             g_0_yy_0_zz_1,  \
                             wp_x,           \
                             c_exps,         \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyy_0_xx_0[i] = 2.0 * g_0_yy_0_x_1[i] * fi_abcd_0 + g_0_yy_0_xx_0[i] * pb_x + g_0_yy_0_xx_1[i] * wp_x[i];

        g_0_xyy_0_xy_0[i] = g_0_yy_0_y_1[i] * fi_abcd_0 + g_0_yy_0_xy_0[i] * pb_x + g_0_yy_0_xy_1[i] * wp_x[i];

        g_0_xyy_0_xz_0[i] = g_0_yy_0_z_1[i] * fi_abcd_0 + g_0_yy_0_xz_0[i] * pb_x + g_0_yy_0_xz_1[i] * wp_x[i];

        g_0_xyy_0_yy_0[i] = g_0_yy_0_yy_0[i] * pb_x + g_0_yy_0_yy_1[i] * wp_x[i];

        g_0_xyy_0_yz_0[i] = g_0_yy_0_yz_0[i] * pb_x + g_0_yy_0_yz_1[i] * wp_x[i];

        g_0_xyy_0_zz_0[i] = g_0_yy_0_zz_0[i] * pb_x + g_0_yy_0_zz_1[i] * wp_x[i];
    }

    /// Set up 24-30 components of targeted buffer : SFSD

    auto g_0_xyz_0_xx_0 = pbuffer.data(idx_eri_0_sfsd + 24);

    auto g_0_xyz_0_xy_0 = pbuffer.data(idx_eri_0_sfsd + 25);

    auto g_0_xyz_0_xz_0 = pbuffer.data(idx_eri_0_sfsd + 26);

    auto g_0_xyz_0_yy_0 = pbuffer.data(idx_eri_0_sfsd + 27);

    auto g_0_xyz_0_yz_0 = pbuffer.data(idx_eri_0_sfsd + 28);

    auto g_0_xyz_0_zz_0 = pbuffer.data(idx_eri_0_sfsd + 29);

#pragma omp simd aligned(g_0_xy_0_xy_0,      \
                             g_0_xy_0_xy_1,  \
                             g_0_xyz_0_xx_0, \
                             g_0_xyz_0_xy_0, \
                             g_0_xyz_0_xz_0, \
                             g_0_xyz_0_yy_0, \
                             g_0_xyz_0_yz_0, \
                             g_0_xyz_0_zz_0, \
                             g_0_xz_0_xx_0,  \
                             g_0_xz_0_xx_1,  \
                             g_0_xz_0_xz_0,  \
                             g_0_xz_0_xz_1,  \
                             g_0_yz_0_yy_0,  \
                             g_0_yz_0_yy_1,  \
                             g_0_yz_0_yz_0,  \
                             g_0_yz_0_yz_1,  \
                             g_0_yz_0_zz_0,  \
                             g_0_yz_0_zz_1,  \
                             wp_x,           \
                             wp_y,           \
                             wp_z,           \
                             c_exps,         \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        g_0_xyz_0_xx_0[i] = g_0_xz_0_xx_0[i] * pb_y + g_0_xz_0_xx_1[i] * wp_y[i];

        g_0_xyz_0_xy_0[i] = g_0_xy_0_xy_0[i] * pb_z + g_0_xy_0_xy_1[i] * wp_z[i];

        g_0_xyz_0_xz_0[i] = g_0_xz_0_xz_0[i] * pb_y + g_0_xz_0_xz_1[i] * wp_y[i];

        g_0_xyz_0_yy_0[i] = g_0_yz_0_yy_0[i] * pb_x + g_0_yz_0_yy_1[i] * wp_x[i];

        g_0_xyz_0_yz_0[i] = g_0_yz_0_yz_0[i] * pb_x + g_0_yz_0_yz_1[i] * wp_x[i];

        g_0_xyz_0_zz_0[i] = g_0_yz_0_zz_0[i] * pb_x + g_0_yz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 30-36 components of targeted buffer : SFSD

    auto g_0_xzz_0_xx_0 = pbuffer.data(idx_eri_0_sfsd + 30);

    auto g_0_xzz_0_xy_0 = pbuffer.data(idx_eri_0_sfsd + 31);

    auto g_0_xzz_0_xz_0 = pbuffer.data(idx_eri_0_sfsd + 32);

    auto g_0_xzz_0_yy_0 = pbuffer.data(idx_eri_0_sfsd + 33);

    auto g_0_xzz_0_yz_0 = pbuffer.data(idx_eri_0_sfsd + 34);

    auto g_0_xzz_0_zz_0 = pbuffer.data(idx_eri_0_sfsd + 35);

#pragma omp simd aligned(g_0_xzz_0_xx_0,     \
                             g_0_xzz_0_xy_0, \
                             g_0_xzz_0_xz_0, \
                             g_0_xzz_0_yy_0, \
                             g_0_xzz_0_yz_0, \
                             g_0_xzz_0_zz_0, \
                             g_0_zz_0_x_1,   \
                             g_0_zz_0_xx_0,  \
                             g_0_zz_0_xx_1,  \
                             g_0_zz_0_xy_0,  \
                             g_0_zz_0_xy_1,  \
                             g_0_zz_0_xz_0,  \
                             g_0_zz_0_xz_1,  \
                             g_0_zz_0_y_1,   \
                             g_0_zz_0_yy_0,  \
                             g_0_zz_0_yy_1,  \
                             g_0_zz_0_yz_0,  \
                             g_0_zz_0_yz_1,  \
                             g_0_zz_0_z_1,   \
                             g_0_zz_0_zz_0,  \
                             g_0_zz_0_zz_1,  \
                             wp_x,           \
                             c_exps,         \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzz_0_xx_0[i] = 2.0 * g_0_zz_0_x_1[i] * fi_abcd_0 + g_0_zz_0_xx_0[i] * pb_x + g_0_zz_0_xx_1[i] * wp_x[i];

        g_0_xzz_0_xy_0[i] = g_0_zz_0_y_1[i] * fi_abcd_0 + g_0_zz_0_xy_0[i] * pb_x + g_0_zz_0_xy_1[i] * wp_x[i];

        g_0_xzz_0_xz_0[i] = g_0_zz_0_z_1[i] * fi_abcd_0 + g_0_zz_0_xz_0[i] * pb_x + g_0_zz_0_xz_1[i] * wp_x[i];

        g_0_xzz_0_yy_0[i] = g_0_zz_0_yy_0[i] * pb_x + g_0_zz_0_yy_1[i] * wp_x[i];

        g_0_xzz_0_yz_0[i] = g_0_zz_0_yz_0[i] * pb_x + g_0_zz_0_yz_1[i] * wp_x[i];

        g_0_xzz_0_zz_0[i] = g_0_zz_0_zz_0[i] * pb_x + g_0_zz_0_zz_1[i] * wp_x[i];
    }

    /// Set up 36-42 components of targeted buffer : SFSD

    auto g_0_yyy_0_xx_0 = pbuffer.data(idx_eri_0_sfsd + 36);

    auto g_0_yyy_0_xy_0 = pbuffer.data(idx_eri_0_sfsd + 37);

    auto g_0_yyy_0_xz_0 = pbuffer.data(idx_eri_0_sfsd + 38);

    auto g_0_yyy_0_yy_0 = pbuffer.data(idx_eri_0_sfsd + 39);

    auto g_0_yyy_0_yz_0 = pbuffer.data(idx_eri_0_sfsd + 40);

    auto g_0_yyy_0_zz_0 = pbuffer.data(idx_eri_0_sfsd + 41);

#pragma omp simd aligned(g_0_y_0_xx_0,       \
                             g_0_y_0_xx_1,   \
                             g_0_y_0_xy_0,   \
                             g_0_y_0_xy_1,   \
                             g_0_y_0_xz_0,   \
                             g_0_y_0_xz_1,   \
                             g_0_y_0_yy_0,   \
                             g_0_y_0_yy_1,   \
                             g_0_y_0_yz_0,   \
                             g_0_y_0_yz_1,   \
                             g_0_y_0_zz_0,   \
                             g_0_y_0_zz_1,   \
                             g_0_yy_0_x_1,   \
                             g_0_yy_0_xx_0,  \
                             g_0_yy_0_xx_1,  \
                             g_0_yy_0_xy_0,  \
                             g_0_yy_0_xy_1,  \
                             g_0_yy_0_xz_0,  \
                             g_0_yy_0_xz_1,  \
                             g_0_yy_0_y_1,   \
                             g_0_yy_0_yy_0,  \
                             g_0_yy_0_yy_1,  \
                             g_0_yy_0_yz_0,  \
                             g_0_yy_0_yz_1,  \
                             g_0_yy_0_z_1,   \
                             g_0_yy_0_zz_0,  \
                             g_0_yy_0_zz_1,  \
                             g_0_yyy_0_xx_0, \
                             g_0_yyy_0_xy_0, \
                             g_0_yyy_0_xz_0, \
                             g_0_yyy_0_yy_0, \
                             g_0_yyy_0_yz_0, \
                             g_0_yyy_0_zz_0, \
                             wp_y,           \
                             c_exps,         \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyy_0_xx_0[i] = 2.0 * g_0_y_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xx_1[i] * fti_ab_0 + g_0_yy_0_xx_0[i] * pb_y + g_0_yy_0_xx_1[i] * wp_y[i];

        g_0_yyy_0_xy_0[i] = 2.0 * g_0_y_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xy_1[i] * fti_ab_0 + g_0_yy_0_x_1[i] * fi_abcd_0 +
                            g_0_yy_0_xy_0[i] * pb_y + g_0_yy_0_xy_1[i] * wp_y[i];

        g_0_yyy_0_xz_0[i] = 2.0 * g_0_y_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xz_1[i] * fti_ab_0 + g_0_yy_0_xz_0[i] * pb_y + g_0_yy_0_xz_1[i] * wp_y[i];

        g_0_yyy_0_yy_0[i] = 2.0 * g_0_y_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yy_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_y_1[i] * fi_abcd_0 +
                            g_0_yy_0_yy_0[i] * pb_y + g_0_yy_0_yy_1[i] * wp_y[i];

        g_0_yyy_0_yz_0[i] = 2.0 * g_0_y_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yz_1[i] * fti_ab_0 + g_0_yy_0_z_1[i] * fi_abcd_0 +
                            g_0_yy_0_yz_0[i] * pb_y + g_0_yy_0_yz_1[i] * wp_y[i];

        g_0_yyy_0_zz_0[i] = 2.0 * g_0_y_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_zz_1[i] * fti_ab_0 + g_0_yy_0_zz_0[i] * pb_y + g_0_yy_0_zz_1[i] * wp_y[i];
    }

    /// Set up 42-48 components of targeted buffer : SFSD

    auto g_0_yyz_0_xx_0 = pbuffer.data(idx_eri_0_sfsd + 42);

    auto g_0_yyz_0_xy_0 = pbuffer.data(idx_eri_0_sfsd + 43);

    auto g_0_yyz_0_xz_0 = pbuffer.data(idx_eri_0_sfsd + 44);

    auto g_0_yyz_0_yy_0 = pbuffer.data(idx_eri_0_sfsd + 45);

    auto g_0_yyz_0_yz_0 = pbuffer.data(idx_eri_0_sfsd + 46);

    auto g_0_yyz_0_zz_0 = pbuffer.data(idx_eri_0_sfsd + 47);

#pragma omp simd aligned(g_0_yy_0_x_1,       \
                             g_0_yy_0_xx_0,  \
                             g_0_yy_0_xx_1,  \
                             g_0_yy_0_xy_0,  \
                             g_0_yy_0_xy_1,  \
                             g_0_yy_0_xz_0,  \
                             g_0_yy_0_xz_1,  \
                             g_0_yy_0_y_1,   \
                             g_0_yy_0_yy_0,  \
                             g_0_yy_0_yy_1,  \
                             g_0_yy_0_yz_0,  \
                             g_0_yy_0_yz_1,  \
                             g_0_yy_0_z_1,   \
                             g_0_yy_0_zz_0,  \
                             g_0_yy_0_zz_1,  \
                             g_0_yyz_0_xx_0, \
                             g_0_yyz_0_xy_0, \
                             g_0_yyz_0_xz_0, \
                             g_0_yyz_0_yy_0, \
                             g_0_yyz_0_yz_0, \
                             g_0_yyz_0_zz_0, \
                             wp_z,           \
                             c_exps,         \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyz_0_xx_0[i] = g_0_yy_0_xx_0[i] * pb_z + g_0_yy_0_xx_1[i] * wp_z[i];

        g_0_yyz_0_xy_0[i] = g_0_yy_0_xy_0[i] * pb_z + g_0_yy_0_xy_1[i] * wp_z[i];

        g_0_yyz_0_xz_0[i] = g_0_yy_0_x_1[i] * fi_abcd_0 + g_0_yy_0_xz_0[i] * pb_z + g_0_yy_0_xz_1[i] * wp_z[i];

        g_0_yyz_0_yy_0[i] = g_0_yy_0_yy_0[i] * pb_z + g_0_yy_0_yy_1[i] * wp_z[i];

        g_0_yyz_0_yz_0[i] = g_0_yy_0_y_1[i] * fi_abcd_0 + g_0_yy_0_yz_0[i] * pb_z + g_0_yy_0_yz_1[i] * wp_z[i];

        g_0_yyz_0_zz_0[i] = 2.0 * g_0_yy_0_z_1[i] * fi_abcd_0 + g_0_yy_0_zz_0[i] * pb_z + g_0_yy_0_zz_1[i] * wp_z[i];
    }

    /// Set up 48-54 components of targeted buffer : SFSD

    auto g_0_yzz_0_xx_0 = pbuffer.data(idx_eri_0_sfsd + 48);

    auto g_0_yzz_0_xy_0 = pbuffer.data(idx_eri_0_sfsd + 49);

    auto g_0_yzz_0_xz_0 = pbuffer.data(idx_eri_0_sfsd + 50);

    auto g_0_yzz_0_yy_0 = pbuffer.data(idx_eri_0_sfsd + 51);

    auto g_0_yzz_0_yz_0 = pbuffer.data(idx_eri_0_sfsd + 52);

    auto g_0_yzz_0_zz_0 = pbuffer.data(idx_eri_0_sfsd + 53);

#pragma omp simd aligned(g_0_yzz_0_xx_0,     \
                             g_0_yzz_0_xy_0, \
                             g_0_yzz_0_xz_0, \
                             g_0_yzz_0_yy_0, \
                             g_0_yzz_0_yz_0, \
                             g_0_yzz_0_zz_0, \
                             g_0_zz_0_x_1,   \
                             g_0_zz_0_xx_0,  \
                             g_0_zz_0_xx_1,  \
                             g_0_zz_0_xy_0,  \
                             g_0_zz_0_xy_1,  \
                             g_0_zz_0_xz_0,  \
                             g_0_zz_0_xz_1,  \
                             g_0_zz_0_y_1,   \
                             g_0_zz_0_yy_0,  \
                             g_0_zz_0_yy_1,  \
                             g_0_zz_0_yz_0,  \
                             g_0_zz_0_yz_1,  \
                             g_0_zz_0_z_1,   \
                             g_0_zz_0_zz_0,  \
                             g_0_zz_0_zz_1,  \
                             wp_y,           \
                             c_exps,         \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzz_0_xx_0[i] = g_0_zz_0_xx_0[i] * pb_y + g_0_zz_0_xx_1[i] * wp_y[i];

        g_0_yzz_0_xy_0[i] = g_0_zz_0_x_1[i] * fi_abcd_0 + g_0_zz_0_xy_0[i] * pb_y + g_0_zz_0_xy_1[i] * wp_y[i];

        g_0_yzz_0_xz_0[i] = g_0_zz_0_xz_0[i] * pb_y + g_0_zz_0_xz_1[i] * wp_y[i];

        g_0_yzz_0_yy_0[i] = 2.0 * g_0_zz_0_y_1[i] * fi_abcd_0 + g_0_zz_0_yy_0[i] * pb_y + g_0_zz_0_yy_1[i] * wp_y[i];

        g_0_yzz_0_yz_0[i] = g_0_zz_0_z_1[i] * fi_abcd_0 + g_0_zz_0_yz_0[i] * pb_y + g_0_zz_0_yz_1[i] * wp_y[i];

        g_0_yzz_0_zz_0[i] = g_0_zz_0_zz_0[i] * pb_y + g_0_zz_0_zz_1[i] * wp_y[i];
    }

    /// Set up 54-60 components of targeted buffer : SFSD

    auto g_0_zzz_0_xx_0 = pbuffer.data(idx_eri_0_sfsd + 54);

    auto g_0_zzz_0_xy_0 = pbuffer.data(idx_eri_0_sfsd + 55);

    auto g_0_zzz_0_xz_0 = pbuffer.data(idx_eri_0_sfsd + 56);

    auto g_0_zzz_0_yy_0 = pbuffer.data(idx_eri_0_sfsd + 57);

    auto g_0_zzz_0_yz_0 = pbuffer.data(idx_eri_0_sfsd + 58);

    auto g_0_zzz_0_zz_0 = pbuffer.data(idx_eri_0_sfsd + 59);

#pragma omp simd aligned(g_0_z_0_xx_0,       \
                             g_0_z_0_xx_1,   \
                             g_0_z_0_xy_0,   \
                             g_0_z_0_xy_1,   \
                             g_0_z_0_xz_0,   \
                             g_0_z_0_xz_1,   \
                             g_0_z_0_yy_0,   \
                             g_0_z_0_yy_1,   \
                             g_0_z_0_yz_0,   \
                             g_0_z_0_yz_1,   \
                             g_0_z_0_zz_0,   \
                             g_0_z_0_zz_1,   \
                             g_0_zz_0_x_1,   \
                             g_0_zz_0_xx_0,  \
                             g_0_zz_0_xx_1,  \
                             g_0_zz_0_xy_0,  \
                             g_0_zz_0_xy_1,  \
                             g_0_zz_0_xz_0,  \
                             g_0_zz_0_xz_1,  \
                             g_0_zz_0_y_1,   \
                             g_0_zz_0_yy_0,  \
                             g_0_zz_0_yy_1,  \
                             g_0_zz_0_yz_0,  \
                             g_0_zz_0_yz_1,  \
                             g_0_zz_0_z_1,   \
                             g_0_zz_0_zz_0,  \
                             g_0_zz_0_zz_1,  \
                             g_0_zzz_0_xx_0, \
                             g_0_zzz_0_xy_0, \
                             g_0_zzz_0_xz_0, \
                             g_0_zzz_0_yy_0, \
                             g_0_zzz_0_yz_0, \
                             g_0_zzz_0_zz_0, \
                             wp_z,           \
                             c_exps,         \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzz_0_xx_0[i] = 2.0 * g_0_z_0_xx_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xx_1[i] * fti_ab_0 + g_0_zz_0_xx_0[i] * pb_z + g_0_zz_0_xx_1[i] * wp_z[i];

        g_0_zzz_0_xy_0[i] = 2.0 * g_0_z_0_xy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xy_1[i] * fti_ab_0 + g_0_zz_0_xy_0[i] * pb_z + g_0_zz_0_xy_1[i] * wp_z[i];

        g_0_zzz_0_xz_0[i] = 2.0 * g_0_z_0_xz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xz_1[i] * fti_ab_0 + g_0_zz_0_x_1[i] * fi_abcd_0 +
                            g_0_zz_0_xz_0[i] * pb_z + g_0_zz_0_xz_1[i] * wp_z[i];

        g_0_zzz_0_yy_0[i] = 2.0 * g_0_z_0_yy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yy_1[i] * fti_ab_0 + g_0_zz_0_yy_0[i] * pb_z + g_0_zz_0_yy_1[i] * wp_z[i];

        g_0_zzz_0_yz_0[i] = 2.0 * g_0_z_0_yz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yz_1[i] * fti_ab_0 + g_0_zz_0_y_1[i] * fi_abcd_0 +
                            g_0_zz_0_yz_0[i] * pb_z + g_0_zz_0_yz_1[i] * wp_z[i];

        g_0_zzz_0_zz_0[i] = 2.0 * g_0_z_0_zz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_zz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_z_1[i] * fi_abcd_0 +
                            g_0_zz_0_zz_0[i] * pb_z + g_0_zz_0_zz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
