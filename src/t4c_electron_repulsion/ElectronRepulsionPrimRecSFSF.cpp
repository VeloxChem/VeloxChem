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

#include "ElectronRepulsionPrimRecSFSF.hpp"

namespace erirec {  // erirec namespace

auto
comp_prim_electron_repulsion_sfsf(CSimdArray<double>&   pbuffer,
                                  const size_t          idx_eri_0_sfsf,
                                  size_t                idx_eri_0_spsf,
                                  size_t                idx_eri_1_spsf,
                                  size_t                idx_eri_1_sdsd,
                                  size_t                idx_eri_0_sdsf,
                                  size_t                idx_eri_1_sdsf,
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

    /// Set up components of auxilary buffer : SPSF

    auto g_0_x_0_xxx_0 = pbuffer.data(idx_eri_0_spsf);

    auto g_0_x_0_xxy_0 = pbuffer.data(idx_eri_0_spsf + 1);

    auto g_0_x_0_xxz_0 = pbuffer.data(idx_eri_0_spsf + 2);

    auto g_0_x_0_xyy_0 = pbuffer.data(idx_eri_0_spsf + 3);

    auto g_0_x_0_xyz_0 = pbuffer.data(idx_eri_0_spsf + 4);

    auto g_0_x_0_xzz_0 = pbuffer.data(idx_eri_0_spsf + 5);

    auto g_0_x_0_yyy_0 = pbuffer.data(idx_eri_0_spsf + 6);

    auto g_0_x_0_yyz_0 = pbuffer.data(idx_eri_0_spsf + 7);

    auto g_0_x_0_yzz_0 = pbuffer.data(idx_eri_0_spsf + 8);

    auto g_0_x_0_zzz_0 = pbuffer.data(idx_eri_0_spsf + 9);

    auto g_0_y_0_xxx_0 = pbuffer.data(idx_eri_0_spsf + 10);

    auto g_0_y_0_xxy_0 = pbuffer.data(idx_eri_0_spsf + 11);

    auto g_0_y_0_xxz_0 = pbuffer.data(idx_eri_0_spsf + 12);

    auto g_0_y_0_xyy_0 = pbuffer.data(idx_eri_0_spsf + 13);

    auto g_0_y_0_xyz_0 = pbuffer.data(idx_eri_0_spsf + 14);

    auto g_0_y_0_xzz_0 = pbuffer.data(idx_eri_0_spsf + 15);

    auto g_0_y_0_yyy_0 = pbuffer.data(idx_eri_0_spsf + 16);

    auto g_0_y_0_yyz_0 = pbuffer.data(idx_eri_0_spsf + 17);

    auto g_0_y_0_yzz_0 = pbuffer.data(idx_eri_0_spsf + 18);

    auto g_0_y_0_zzz_0 = pbuffer.data(idx_eri_0_spsf + 19);

    auto g_0_z_0_xxx_0 = pbuffer.data(idx_eri_0_spsf + 20);

    auto g_0_z_0_xxy_0 = pbuffer.data(idx_eri_0_spsf + 21);

    auto g_0_z_0_xxz_0 = pbuffer.data(idx_eri_0_spsf + 22);

    auto g_0_z_0_xyy_0 = pbuffer.data(idx_eri_0_spsf + 23);

    auto g_0_z_0_xyz_0 = pbuffer.data(idx_eri_0_spsf + 24);

    auto g_0_z_0_xzz_0 = pbuffer.data(idx_eri_0_spsf + 25);

    auto g_0_z_0_yyy_0 = pbuffer.data(idx_eri_0_spsf + 26);

    auto g_0_z_0_yyz_0 = pbuffer.data(idx_eri_0_spsf + 27);

    auto g_0_z_0_yzz_0 = pbuffer.data(idx_eri_0_spsf + 28);

    auto g_0_z_0_zzz_0 = pbuffer.data(idx_eri_0_spsf + 29);

    /// Set up components of auxilary buffer : SPSF

    auto g_0_x_0_xxx_1 = pbuffer.data(idx_eri_1_spsf);

    auto g_0_x_0_xxy_1 = pbuffer.data(idx_eri_1_spsf + 1);

    auto g_0_x_0_xxz_1 = pbuffer.data(idx_eri_1_spsf + 2);

    auto g_0_x_0_xyy_1 = pbuffer.data(idx_eri_1_spsf + 3);

    auto g_0_x_0_xyz_1 = pbuffer.data(idx_eri_1_spsf + 4);

    auto g_0_x_0_xzz_1 = pbuffer.data(idx_eri_1_spsf + 5);

    auto g_0_x_0_yyy_1 = pbuffer.data(idx_eri_1_spsf + 6);

    auto g_0_x_0_yyz_1 = pbuffer.data(idx_eri_1_spsf + 7);

    auto g_0_x_0_yzz_1 = pbuffer.data(idx_eri_1_spsf + 8);

    auto g_0_x_0_zzz_1 = pbuffer.data(idx_eri_1_spsf + 9);

    auto g_0_y_0_xxx_1 = pbuffer.data(idx_eri_1_spsf + 10);

    auto g_0_y_0_xxy_1 = pbuffer.data(idx_eri_1_spsf + 11);

    auto g_0_y_0_xxz_1 = pbuffer.data(idx_eri_1_spsf + 12);

    auto g_0_y_0_xyy_1 = pbuffer.data(idx_eri_1_spsf + 13);

    auto g_0_y_0_xyz_1 = pbuffer.data(idx_eri_1_spsf + 14);

    auto g_0_y_0_xzz_1 = pbuffer.data(idx_eri_1_spsf + 15);

    auto g_0_y_0_yyy_1 = pbuffer.data(idx_eri_1_spsf + 16);

    auto g_0_y_0_yyz_1 = pbuffer.data(idx_eri_1_spsf + 17);

    auto g_0_y_0_yzz_1 = pbuffer.data(idx_eri_1_spsf + 18);

    auto g_0_y_0_zzz_1 = pbuffer.data(idx_eri_1_spsf + 19);

    auto g_0_z_0_xxx_1 = pbuffer.data(idx_eri_1_spsf + 20);

    auto g_0_z_0_xxy_1 = pbuffer.data(idx_eri_1_spsf + 21);

    auto g_0_z_0_xxz_1 = pbuffer.data(idx_eri_1_spsf + 22);

    auto g_0_z_0_xyy_1 = pbuffer.data(idx_eri_1_spsf + 23);

    auto g_0_z_0_xyz_1 = pbuffer.data(idx_eri_1_spsf + 24);

    auto g_0_z_0_xzz_1 = pbuffer.data(idx_eri_1_spsf + 25);

    auto g_0_z_0_yyy_1 = pbuffer.data(idx_eri_1_spsf + 26);

    auto g_0_z_0_yyz_1 = pbuffer.data(idx_eri_1_spsf + 27);

    auto g_0_z_0_yzz_1 = pbuffer.data(idx_eri_1_spsf + 28);

    auto g_0_z_0_zzz_1 = pbuffer.data(idx_eri_1_spsf + 29);

    /// Set up components of auxilary buffer : SDSD

    auto g_0_xx_0_xx_1 = pbuffer.data(idx_eri_1_sdsd);

    auto g_0_xx_0_xy_1 = pbuffer.data(idx_eri_1_sdsd + 1);

    auto g_0_xx_0_xz_1 = pbuffer.data(idx_eri_1_sdsd + 2);

    auto g_0_xx_0_yy_1 = pbuffer.data(idx_eri_1_sdsd + 3);

    auto g_0_xx_0_yz_1 = pbuffer.data(idx_eri_1_sdsd + 4);

    auto g_0_xx_0_zz_1 = pbuffer.data(idx_eri_1_sdsd + 5);

    auto g_0_yy_0_xx_1 = pbuffer.data(idx_eri_1_sdsd + 18);

    auto g_0_yy_0_xy_1 = pbuffer.data(idx_eri_1_sdsd + 19);

    auto g_0_yy_0_xz_1 = pbuffer.data(idx_eri_1_sdsd + 20);

    auto g_0_yy_0_yy_1 = pbuffer.data(idx_eri_1_sdsd + 21);

    auto g_0_yy_0_yz_1 = pbuffer.data(idx_eri_1_sdsd + 22);

    auto g_0_yy_0_zz_1 = pbuffer.data(idx_eri_1_sdsd + 23);

    auto g_0_yz_0_yz_1 = pbuffer.data(idx_eri_1_sdsd + 28);

    auto g_0_zz_0_xx_1 = pbuffer.data(idx_eri_1_sdsd + 30);

    auto g_0_zz_0_xy_1 = pbuffer.data(idx_eri_1_sdsd + 31);

    auto g_0_zz_0_xz_1 = pbuffer.data(idx_eri_1_sdsd + 32);

    auto g_0_zz_0_yy_1 = pbuffer.data(idx_eri_1_sdsd + 33);

    auto g_0_zz_0_yz_1 = pbuffer.data(idx_eri_1_sdsd + 34);

    auto g_0_zz_0_zz_1 = pbuffer.data(idx_eri_1_sdsd + 35);

    /// Set up components of auxilary buffer : SDSF

    auto g_0_xx_0_xxx_0 = pbuffer.data(idx_eri_0_sdsf);

    auto g_0_xx_0_xxy_0 = pbuffer.data(idx_eri_0_sdsf + 1);

    auto g_0_xx_0_xxz_0 = pbuffer.data(idx_eri_0_sdsf + 2);

    auto g_0_xx_0_xyy_0 = pbuffer.data(idx_eri_0_sdsf + 3);

    auto g_0_xx_0_xyz_0 = pbuffer.data(idx_eri_0_sdsf + 4);

    auto g_0_xx_0_xzz_0 = pbuffer.data(idx_eri_0_sdsf + 5);

    auto g_0_xx_0_yyy_0 = pbuffer.data(idx_eri_0_sdsf + 6);

    auto g_0_xx_0_yyz_0 = pbuffer.data(idx_eri_0_sdsf + 7);

    auto g_0_xx_0_yzz_0 = pbuffer.data(idx_eri_0_sdsf + 8);

    auto g_0_xx_0_zzz_0 = pbuffer.data(idx_eri_0_sdsf + 9);

    auto g_0_xy_0_xxy_0 = pbuffer.data(idx_eri_0_sdsf + 11);

    auto g_0_xy_0_xyy_0 = pbuffer.data(idx_eri_0_sdsf + 13);

    auto g_0_xz_0_xxx_0 = pbuffer.data(idx_eri_0_sdsf + 20);

    auto g_0_xz_0_xxz_0 = pbuffer.data(idx_eri_0_sdsf + 22);

    auto g_0_xz_0_xzz_0 = pbuffer.data(idx_eri_0_sdsf + 25);

    auto g_0_yy_0_xxx_0 = pbuffer.data(idx_eri_0_sdsf + 30);

    auto g_0_yy_0_xxy_0 = pbuffer.data(idx_eri_0_sdsf + 31);

    auto g_0_yy_0_xxz_0 = pbuffer.data(idx_eri_0_sdsf + 32);

    auto g_0_yy_0_xyy_0 = pbuffer.data(idx_eri_0_sdsf + 33);

    auto g_0_yy_0_xyz_0 = pbuffer.data(idx_eri_0_sdsf + 34);

    auto g_0_yy_0_xzz_0 = pbuffer.data(idx_eri_0_sdsf + 35);

    auto g_0_yy_0_yyy_0 = pbuffer.data(idx_eri_0_sdsf + 36);

    auto g_0_yy_0_yyz_0 = pbuffer.data(idx_eri_0_sdsf + 37);

    auto g_0_yy_0_yzz_0 = pbuffer.data(idx_eri_0_sdsf + 38);

    auto g_0_yy_0_zzz_0 = pbuffer.data(idx_eri_0_sdsf + 39);

    auto g_0_yz_0_xyz_0 = pbuffer.data(idx_eri_0_sdsf + 44);

    auto g_0_yz_0_yyy_0 = pbuffer.data(idx_eri_0_sdsf + 46);

    auto g_0_yz_0_yyz_0 = pbuffer.data(idx_eri_0_sdsf + 47);

    auto g_0_yz_0_yzz_0 = pbuffer.data(idx_eri_0_sdsf + 48);

    auto g_0_yz_0_zzz_0 = pbuffer.data(idx_eri_0_sdsf + 49);

    auto g_0_zz_0_xxx_0 = pbuffer.data(idx_eri_0_sdsf + 50);

    auto g_0_zz_0_xxy_0 = pbuffer.data(idx_eri_0_sdsf + 51);

    auto g_0_zz_0_xxz_0 = pbuffer.data(idx_eri_0_sdsf + 52);

    auto g_0_zz_0_xyy_0 = pbuffer.data(idx_eri_0_sdsf + 53);

    auto g_0_zz_0_xyz_0 = pbuffer.data(idx_eri_0_sdsf + 54);

    auto g_0_zz_0_xzz_0 = pbuffer.data(idx_eri_0_sdsf + 55);

    auto g_0_zz_0_yyy_0 = pbuffer.data(idx_eri_0_sdsf + 56);

    auto g_0_zz_0_yyz_0 = pbuffer.data(idx_eri_0_sdsf + 57);

    auto g_0_zz_0_yzz_0 = pbuffer.data(idx_eri_0_sdsf + 58);

    auto g_0_zz_0_zzz_0 = pbuffer.data(idx_eri_0_sdsf + 59);

    /// Set up components of auxilary buffer : SDSF

    auto g_0_xx_0_xxx_1 = pbuffer.data(idx_eri_1_sdsf);

    auto g_0_xx_0_xxy_1 = pbuffer.data(idx_eri_1_sdsf + 1);

    auto g_0_xx_0_xxz_1 = pbuffer.data(idx_eri_1_sdsf + 2);

    auto g_0_xx_0_xyy_1 = pbuffer.data(idx_eri_1_sdsf + 3);

    auto g_0_xx_0_xyz_1 = pbuffer.data(idx_eri_1_sdsf + 4);

    auto g_0_xx_0_xzz_1 = pbuffer.data(idx_eri_1_sdsf + 5);

    auto g_0_xx_0_yyy_1 = pbuffer.data(idx_eri_1_sdsf + 6);

    auto g_0_xx_0_yyz_1 = pbuffer.data(idx_eri_1_sdsf + 7);

    auto g_0_xx_0_yzz_1 = pbuffer.data(idx_eri_1_sdsf + 8);

    auto g_0_xx_0_zzz_1 = pbuffer.data(idx_eri_1_sdsf + 9);

    auto g_0_xy_0_xxy_1 = pbuffer.data(idx_eri_1_sdsf + 11);

    auto g_0_xy_0_xyy_1 = pbuffer.data(idx_eri_1_sdsf + 13);

    auto g_0_xz_0_xxx_1 = pbuffer.data(idx_eri_1_sdsf + 20);

    auto g_0_xz_0_xxz_1 = pbuffer.data(idx_eri_1_sdsf + 22);

    auto g_0_xz_0_xzz_1 = pbuffer.data(idx_eri_1_sdsf + 25);

    auto g_0_yy_0_xxx_1 = pbuffer.data(idx_eri_1_sdsf + 30);

    auto g_0_yy_0_xxy_1 = pbuffer.data(idx_eri_1_sdsf + 31);

    auto g_0_yy_0_xxz_1 = pbuffer.data(idx_eri_1_sdsf + 32);

    auto g_0_yy_0_xyy_1 = pbuffer.data(idx_eri_1_sdsf + 33);

    auto g_0_yy_0_xyz_1 = pbuffer.data(idx_eri_1_sdsf + 34);

    auto g_0_yy_0_xzz_1 = pbuffer.data(idx_eri_1_sdsf + 35);

    auto g_0_yy_0_yyy_1 = pbuffer.data(idx_eri_1_sdsf + 36);

    auto g_0_yy_0_yyz_1 = pbuffer.data(idx_eri_1_sdsf + 37);

    auto g_0_yy_0_yzz_1 = pbuffer.data(idx_eri_1_sdsf + 38);

    auto g_0_yy_0_zzz_1 = pbuffer.data(idx_eri_1_sdsf + 39);

    auto g_0_yz_0_xyz_1 = pbuffer.data(idx_eri_1_sdsf + 44);

    auto g_0_yz_0_yyy_1 = pbuffer.data(idx_eri_1_sdsf + 46);

    auto g_0_yz_0_yyz_1 = pbuffer.data(idx_eri_1_sdsf + 47);

    auto g_0_yz_0_yzz_1 = pbuffer.data(idx_eri_1_sdsf + 48);

    auto g_0_yz_0_zzz_1 = pbuffer.data(idx_eri_1_sdsf + 49);

    auto g_0_zz_0_xxx_1 = pbuffer.data(idx_eri_1_sdsf + 50);

    auto g_0_zz_0_xxy_1 = pbuffer.data(idx_eri_1_sdsf + 51);

    auto g_0_zz_0_xxz_1 = pbuffer.data(idx_eri_1_sdsf + 52);

    auto g_0_zz_0_xyy_1 = pbuffer.data(idx_eri_1_sdsf + 53);

    auto g_0_zz_0_xyz_1 = pbuffer.data(idx_eri_1_sdsf + 54);

    auto g_0_zz_0_xzz_1 = pbuffer.data(idx_eri_1_sdsf + 55);

    auto g_0_zz_0_yyy_1 = pbuffer.data(idx_eri_1_sdsf + 56);

    auto g_0_zz_0_yyz_1 = pbuffer.data(idx_eri_1_sdsf + 57);

    auto g_0_zz_0_yzz_1 = pbuffer.data(idx_eri_1_sdsf + 58);

    auto g_0_zz_0_zzz_1 = pbuffer.data(idx_eri_1_sdsf + 59);

    /// Set up 0-10 components of targeted buffer : SFSF

    auto g_0_xxx_0_xxx_0 = pbuffer.data(idx_eri_0_sfsf);

    auto g_0_xxx_0_xxy_0 = pbuffer.data(idx_eri_0_sfsf + 1);

    auto g_0_xxx_0_xxz_0 = pbuffer.data(idx_eri_0_sfsf + 2);

    auto g_0_xxx_0_xyy_0 = pbuffer.data(idx_eri_0_sfsf + 3);

    auto g_0_xxx_0_xyz_0 = pbuffer.data(idx_eri_0_sfsf + 4);

    auto g_0_xxx_0_xzz_0 = pbuffer.data(idx_eri_0_sfsf + 5);

    auto g_0_xxx_0_yyy_0 = pbuffer.data(idx_eri_0_sfsf + 6);

    auto g_0_xxx_0_yyz_0 = pbuffer.data(idx_eri_0_sfsf + 7);

    auto g_0_xxx_0_yzz_0 = pbuffer.data(idx_eri_0_sfsf + 8);

    auto g_0_xxx_0_zzz_0 = pbuffer.data(idx_eri_0_sfsf + 9);

#pragma omp simd aligned(g_0_x_0_xxx_0,       \
                             g_0_x_0_xxx_1,   \
                             g_0_x_0_xxy_0,   \
                             g_0_x_0_xxy_1,   \
                             g_0_x_0_xxz_0,   \
                             g_0_x_0_xxz_1,   \
                             g_0_x_0_xyy_0,   \
                             g_0_x_0_xyy_1,   \
                             g_0_x_0_xyz_0,   \
                             g_0_x_0_xyz_1,   \
                             g_0_x_0_xzz_0,   \
                             g_0_x_0_xzz_1,   \
                             g_0_x_0_yyy_0,   \
                             g_0_x_0_yyy_1,   \
                             g_0_x_0_yyz_0,   \
                             g_0_x_0_yyz_1,   \
                             g_0_x_0_yzz_0,   \
                             g_0_x_0_yzz_1,   \
                             g_0_x_0_zzz_0,   \
                             g_0_x_0_zzz_1,   \
                             g_0_xx_0_xx_1,   \
                             g_0_xx_0_xxx_0,  \
                             g_0_xx_0_xxx_1,  \
                             g_0_xx_0_xxy_0,  \
                             g_0_xx_0_xxy_1,  \
                             g_0_xx_0_xxz_0,  \
                             g_0_xx_0_xxz_1,  \
                             g_0_xx_0_xy_1,   \
                             g_0_xx_0_xyy_0,  \
                             g_0_xx_0_xyy_1,  \
                             g_0_xx_0_xyz_0,  \
                             g_0_xx_0_xyz_1,  \
                             g_0_xx_0_xz_1,   \
                             g_0_xx_0_xzz_0,  \
                             g_0_xx_0_xzz_1,  \
                             g_0_xx_0_yy_1,   \
                             g_0_xx_0_yyy_0,  \
                             g_0_xx_0_yyy_1,  \
                             g_0_xx_0_yyz_0,  \
                             g_0_xx_0_yyz_1,  \
                             g_0_xx_0_yz_1,   \
                             g_0_xx_0_yzz_0,  \
                             g_0_xx_0_yzz_1,  \
                             g_0_xx_0_zz_1,   \
                             g_0_xx_0_zzz_0,  \
                             g_0_xx_0_zzz_1,  \
                             g_0_xxx_0_xxx_0, \
                             g_0_xxx_0_xxy_0, \
                             g_0_xxx_0_xxz_0, \
                             g_0_xxx_0_xyy_0, \
                             g_0_xxx_0_xyz_0, \
                             g_0_xxx_0_xzz_0, \
                             g_0_xxx_0_yyy_0, \
                             g_0_xxx_0_yyz_0, \
                             g_0_xxx_0_yzz_0, \
                             g_0_xxx_0_zzz_0, \
                             wp_x,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_xxx_0_xxx_0[i] = 2.0 * g_0_x_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxx_1[i] * fti_ab_0 + 3.0 * g_0_xx_0_xx_1[i] * fi_abcd_0 +
                             g_0_xx_0_xxx_0[i] * pb_x + g_0_xx_0_xxx_1[i] * wp_x[i];

        g_0_xxx_0_xxy_0[i] = 2.0 * g_0_x_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxy_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_xy_1[i] * fi_abcd_0 +
                             g_0_xx_0_xxy_0[i] * pb_x + g_0_xx_0_xxy_1[i] * wp_x[i];

        g_0_xxx_0_xxz_0[i] = 2.0 * g_0_x_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xxz_1[i] * fti_ab_0 + 2.0 * g_0_xx_0_xz_1[i] * fi_abcd_0 +
                             g_0_xx_0_xxz_0[i] * pb_x + g_0_xx_0_xxz_1[i] * wp_x[i];

        g_0_xxx_0_xyy_0[i] = 2.0 * g_0_x_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyy_1[i] * fti_ab_0 + g_0_xx_0_yy_1[i] * fi_abcd_0 +
                             g_0_xx_0_xyy_0[i] * pb_x + g_0_xx_0_xyy_1[i] * wp_x[i];

        g_0_xxx_0_xyz_0[i] = 2.0 * g_0_x_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xyz_1[i] * fti_ab_0 + g_0_xx_0_yz_1[i] * fi_abcd_0 +
                             g_0_xx_0_xyz_0[i] * pb_x + g_0_xx_0_xyz_1[i] * wp_x[i];

        g_0_xxx_0_xzz_0[i] = 2.0 * g_0_x_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_xzz_1[i] * fti_ab_0 + g_0_xx_0_zz_1[i] * fi_abcd_0 +
                             g_0_xx_0_xzz_0[i] * pb_x + g_0_xx_0_xzz_1[i] * wp_x[i];

        g_0_xxx_0_yyy_0[i] =
            2.0 * g_0_x_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyy_1[i] * fti_ab_0 + g_0_xx_0_yyy_0[i] * pb_x + g_0_xx_0_yyy_1[i] * wp_x[i];

        g_0_xxx_0_yyz_0[i] =
            2.0 * g_0_x_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yyz_1[i] * fti_ab_0 + g_0_xx_0_yyz_0[i] * pb_x + g_0_xx_0_yyz_1[i] * wp_x[i];

        g_0_xxx_0_yzz_0[i] =
            2.0 * g_0_x_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_yzz_1[i] * fti_ab_0 + g_0_xx_0_yzz_0[i] * pb_x + g_0_xx_0_yzz_1[i] * wp_x[i];

        g_0_xxx_0_zzz_0[i] =
            2.0 * g_0_x_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_x_0_zzz_1[i] * fti_ab_0 + g_0_xx_0_zzz_0[i] * pb_x + g_0_xx_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 10-20 components of targeted buffer : SFSF

    auto g_0_xxy_0_xxx_0 = pbuffer.data(idx_eri_0_sfsf + 10);

    auto g_0_xxy_0_xxy_0 = pbuffer.data(idx_eri_0_sfsf + 11);

    auto g_0_xxy_0_xxz_0 = pbuffer.data(idx_eri_0_sfsf + 12);

    auto g_0_xxy_0_xyy_0 = pbuffer.data(idx_eri_0_sfsf + 13);

    auto g_0_xxy_0_xyz_0 = pbuffer.data(idx_eri_0_sfsf + 14);

    auto g_0_xxy_0_xzz_0 = pbuffer.data(idx_eri_0_sfsf + 15);

    auto g_0_xxy_0_yyy_0 = pbuffer.data(idx_eri_0_sfsf + 16);

    auto g_0_xxy_0_yyz_0 = pbuffer.data(idx_eri_0_sfsf + 17);

    auto g_0_xxy_0_yzz_0 = pbuffer.data(idx_eri_0_sfsf + 18);

    auto g_0_xxy_0_zzz_0 = pbuffer.data(idx_eri_0_sfsf + 19);

#pragma omp simd aligned(g_0_xx_0_xx_1,       \
                             g_0_xx_0_xxx_0,  \
                             g_0_xx_0_xxx_1,  \
                             g_0_xx_0_xxy_0,  \
                             g_0_xx_0_xxy_1,  \
                             g_0_xx_0_xxz_0,  \
                             g_0_xx_0_xxz_1,  \
                             g_0_xx_0_xy_1,   \
                             g_0_xx_0_xyy_0,  \
                             g_0_xx_0_xyy_1,  \
                             g_0_xx_0_xyz_0,  \
                             g_0_xx_0_xyz_1,  \
                             g_0_xx_0_xz_1,   \
                             g_0_xx_0_xzz_0,  \
                             g_0_xx_0_xzz_1,  \
                             g_0_xx_0_yy_1,   \
                             g_0_xx_0_yyy_0,  \
                             g_0_xx_0_yyy_1,  \
                             g_0_xx_0_yyz_0,  \
                             g_0_xx_0_yyz_1,  \
                             g_0_xx_0_yz_1,   \
                             g_0_xx_0_yzz_0,  \
                             g_0_xx_0_yzz_1,  \
                             g_0_xx_0_zz_1,   \
                             g_0_xx_0_zzz_0,  \
                             g_0_xx_0_zzz_1,  \
                             g_0_xxy_0_xxx_0, \
                             g_0_xxy_0_xxy_0, \
                             g_0_xxy_0_xxz_0, \
                             g_0_xxy_0_xyy_0, \
                             g_0_xxy_0_xyz_0, \
                             g_0_xxy_0_xzz_0, \
                             g_0_xxy_0_yyy_0, \
                             g_0_xxy_0_yyz_0, \
                             g_0_xxy_0_yzz_0, \
                             g_0_xxy_0_zzz_0, \
                             wp_y,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxy_0_xxx_0[i] = g_0_xx_0_xxx_0[i] * pb_y + g_0_xx_0_xxx_1[i] * wp_y[i];

        g_0_xxy_0_xxy_0[i] = g_0_xx_0_xx_1[i] * fi_abcd_0 + g_0_xx_0_xxy_0[i] * pb_y + g_0_xx_0_xxy_1[i] * wp_y[i];

        g_0_xxy_0_xxz_0[i] = g_0_xx_0_xxz_0[i] * pb_y + g_0_xx_0_xxz_1[i] * wp_y[i];

        g_0_xxy_0_xyy_0[i] = 2.0 * g_0_xx_0_xy_1[i] * fi_abcd_0 + g_0_xx_0_xyy_0[i] * pb_y + g_0_xx_0_xyy_1[i] * wp_y[i];

        g_0_xxy_0_xyz_0[i] = g_0_xx_0_xz_1[i] * fi_abcd_0 + g_0_xx_0_xyz_0[i] * pb_y + g_0_xx_0_xyz_1[i] * wp_y[i];

        g_0_xxy_0_xzz_0[i] = g_0_xx_0_xzz_0[i] * pb_y + g_0_xx_0_xzz_1[i] * wp_y[i];

        g_0_xxy_0_yyy_0[i] = 3.0 * g_0_xx_0_yy_1[i] * fi_abcd_0 + g_0_xx_0_yyy_0[i] * pb_y + g_0_xx_0_yyy_1[i] * wp_y[i];

        g_0_xxy_0_yyz_0[i] = 2.0 * g_0_xx_0_yz_1[i] * fi_abcd_0 + g_0_xx_0_yyz_0[i] * pb_y + g_0_xx_0_yyz_1[i] * wp_y[i];

        g_0_xxy_0_yzz_0[i] = g_0_xx_0_zz_1[i] * fi_abcd_0 + g_0_xx_0_yzz_0[i] * pb_y + g_0_xx_0_yzz_1[i] * wp_y[i];

        g_0_xxy_0_zzz_0[i] = g_0_xx_0_zzz_0[i] * pb_y + g_0_xx_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 20-30 components of targeted buffer : SFSF

    auto g_0_xxz_0_xxx_0 = pbuffer.data(idx_eri_0_sfsf + 20);

    auto g_0_xxz_0_xxy_0 = pbuffer.data(idx_eri_0_sfsf + 21);

    auto g_0_xxz_0_xxz_0 = pbuffer.data(idx_eri_0_sfsf + 22);

    auto g_0_xxz_0_xyy_0 = pbuffer.data(idx_eri_0_sfsf + 23);

    auto g_0_xxz_0_xyz_0 = pbuffer.data(idx_eri_0_sfsf + 24);

    auto g_0_xxz_0_xzz_0 = pbuffer.data(idx_eri_0_sfsf + 25);

    auto g_0_xxz_0_yyy_0 = pbuffer.data(idx_eri_0_sfsf + 26);

    auto g_0_xxz_0_yyz_0 = pbuffer.data(idx_eri_0_sfsf + 27);

    auto g_0_xxz_0_yzz_0 = pbuffer.data(idx_eri_0_sfsf + 28);

    auto g_0_xxz_0_zzz_0 = pbuffer.data(idx_eri_0_sfsf + 29);

#pragma omp simd aligned(g_0_xx_0_xx_1,       \
                             g_0_xx_0_xxx_0,  \
                             g_0_xx_0_xxx_1,  \
                             g_0_xx_0_xxy_0,  \
                             g_0_xx_0_xxy_1,  \
                             g_0_xx_0_xxz_0,  \
                             g_0_xx_0_xxz_1,  \
                             g_0_xx_0_xy_1,   \
                             g_0_xx_0_xyy_0,  \
                             g_0_xx_0_xyy_1,  \
                             g_0_xx_0_xyz_0,  \
                             g_0_xx_0_xyz_1,  \
                             g_0_xx_0_xz_1,   \
                             g_0_xx_0_xzz_0,  \
                             g_0_xx_0_xzz_1,  \
                             g_0_xx_0_yy_1,   \
                             g_0_xx_0_yyy_0,  \
                             g_0_xx_0_yyy_1,  \
                             g_0_xx_0_yyz_0,  \
                             g_0_xx_0_yyz_1,  \
                             g_0_xx_0_yz_1,   \
                             g_0_xx_0_yzz_0,  \
                             g_0_xx_0_yzz_1,  \
                             g_0_xx_0_zz_1,   \
                             g_0_xx_0_zzz_0,  \
                             g_0_xx_0_zzz_1,  \
                             g_0_xxz_0_xxx_0, \
                             g_0_xxz_0_xxy_0, \
                             g_0_xxz_0_xxz_0, \
                             g_0_xxz_0_xyy_0, \
                             g_0_xxz_0_xyz_0, \
                             g_0_xxz_0_xzz_0, \
                             g_0_xxz_0_yyy_0, \
                             g_0_xxz_0_yyz_0, \
                             g_0_xxz_0_yzz_0, \
                             g_0_xxz_0_zzz_0, \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xxz_0_xxx_0[i] = g_0_xx_0_xxx_0[i] * pb_z + g_0_xx_0_xxx_1[i] * wp_z[i];

        g_0_xxz_0_xxy_0[i] = g_0_xx_0_xxy_0[i] * pb_z + g_0_xx_0_xxy_1[i] * wp_z[i];

        g_0_xxz_0_xxz_0[i] = g_0_xx_0_xx_1[i] * fi_abcd_0 + g_0_xx_0_xxz_0[i] * pb_z + g_0_xx_0_xxz_1[i] * wp_z[i];

        g_0_xxz_0_xyy_0[i] = g_0_xx_0_xyy_0[i] * pb_z + g_0_xx_0_xyy_1[i] * wp_z[i];

        g_0_xxz_0_xyz_0[i] = g_0_xx_0_xy_1[i] * fi_abcd_0 + g_0_xx_0_xyz_0[i] * pb_z + g_0_xx_0_xyz_1[i] * wp_z[i];

        g_0_xxz_0_xzz_0[i] = 2.0 * g_0_xx_0_xz_1[i] * fi_abcd_0 + g_0_xx_0_xzz_0[i] * pb_z + g_0_xx_0_xzz_1[i] * wp_z[i];

        g_0_xxz_0_yyy_0[i] = g_0_xx_0_yyy_0[i] * pb_z + g_0_xx_0_yyy_1[i] * wp_z[i];

        g_0_xxz_0_yyz_0[i] = g_0_xx_0_yy_1[i] * fi_abcd_0 + g_0_xx_0_yyz_0[i] * pb_z + g_0_xx_0_yyz_1[i] * wp_z[i];

        g_0_xxz_0_yzz_0[i] = 2.0 * g_0_xx_0_yz_1[i] * fi_abcd_0 + g_0_xx_0_yzz_0[i] * pb_z + g_0_xx_0_yzz_1[i] * wp_z[i];

        g_0_xxz_0_zzz_0[i] = 3.0 * g_0_xx_0_zz_1[i] * fi_abcd_0 + g_0_xx_0_zzz_0[i] * pb_z + g_0_xx_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 30-40 components of targeted buffer : SFSF

    auto g_0_xyy_0_xxx_0 = pbuffer.data(idx_eri_0_sfsf + 30);

    auto g_0_xyy_0_xxy_0 = pbuffer.data(idx_eri_0_sfsf + 31);

    auto g_0_xyy_0_xxz_0 = pbuffer.data(idx_eri_0_sfsf + 32);

    auto g_0_xyy_0_xyy_0 = pbuffer.data(idx_eri_0_sfsf + 33);

    auto g_0_xyy_0_xyz_0 = pbuffer.data(idx_eri_0_sfsf + 34);

    auto g_0_xyy_0_xzz_0 = pbuffer.data(idx_eri_0_sfsf + 35);

    auto g_0_xyy_0_yyy_0 = pbuffer.data(idx_eri_0_sfsf + 36);

    auto g_0_xyy_0_yyz_0 = pbuffer.data(idx_eri_0_sfsf + 37);

    auto g_0_xyy_0_yzz_0 = pbuffer.data(idx_eri_0_sfsf + 38);

    auto g_0_xyy_0_zzz_0 = pbuffer.data(idx_eri_0_sfsf + 39);

#pragma omp simd aligned(g_0_xyy_0_xxx_0,     \
                             g_0_xyy_0_xxy_0, \
                             g_0_xyy_0_xxz_0, \
                             g_0_xyy_0_xyy_0, \
                             g_0_xyy_0_xyz_0, \
                             g_0_xyy_0_xzz_0, \
                             g_0_xyy_0_yyy_0, \
                             g_0_xyy_0_yyz_0, \
                             g_0_xyy_0_yzz_0, \
                             g_0_xyy_0_zzz_0, \
                             g_0_yy_0_xx_1,   \
                             g_0_yy_0_xxx_0,  \
                             g_0_yy_0_xxx_1,  \
                             g_0_yy_0_xxy_0,  \
                             g_0_yy_0_xxy_1,  \
                             g_0_yy_0_xxz_0,  \
                             g_0_yy_0_xxz_1,  \
                             g_0_yy_0_xy_1,   \
                             g_0_yy_0_xyy_0,  \
                             g_0_yy_0_xyy_1,  \
                             g_0_yy_0_xyz_0,  \
                             g_0_yy_0_xyz_1,  \
                             g_0_yy_0_xz_1,   \
                             g_0_yy_0_xzz_0,  \
                             g_0_yy_0_xzz_1,  \
                             g_0_yy_0_yy_1,   \
                             g_0_yy_0_yyy_0,  \
                             g_0_yy_0_yyy_1,  \
                             g_0_yy_0_yyz_0,  \
                             g_0_yy_0_yyz_1,  \
                             g_0_yy_0_yz_1,   \
                             g_0_yy_0_yzz_0,  \
                             g_0_yy_0_yzz_1,  \
                             g_0_yy_0_zz_1,   \
                             g_0_yy_0_zzz_0,  \
                             g_0_yy_0_zzz_1,  \
                             wp_x,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyy_0_xxx_0[i] = 3.0 * g_0_yy_0_xx_1[i] * fi_abcd_0 + g_0_yy_0_xxx_0[i] * pb_x + g_0_yy_0_xxx_1[i] * wp_x[i];

        g_0_xyy_0_xxy_0[i] = 2.0 * g_0_yy_0_xy_1[i] * fi_abcd_0 + g_0_yy_0_xxy_0[i] * pb_x + g_0_yy_0_xxy_1[i] * wp_x[i];

        g_0_xyy_0_xxz_0[i] = 2.0 * g_0_yy_0_xz_1[i] * fi_abcd_0 + g_0_yy_0_xxz_0[i] * pb_x + g_0_yy_0_xxz_1[i] * wp_x[i];

        g_0_xyy_0_xyy_0[i] = g_0_yy_0_yy_1[i] * fi_abcd_0 + g_0_yy_0_xyy_0[i] * pb_x + g_0_yy_0_xyy_1[i] * wp_x[i];

        g_0_xyy_0_xyz_0[i] = g_0_yy_0_yz_1[i] * fi_abcd_0 + g_0_yy_0_xyz_0[i] * pb_x + g_0_yy_0_xyz_1[i] * wp_x[i];

        g_0_xyy_0_xzz_0[i] = g_0_yy_0_zz_1[i] * fi_abcd_0 + g_0_yy_0_xzz_0[i] * pb_x + g_0_yy_0_xzz_1[i] * wp_x[i];

        g_0_xyy_0_yyy_0[i] = g_0_yy_0_yyy_0[i] * pb_x + g_0_yy_0_yyy_1[i] * wp_x[i];

        g_0_xyy_0_yyz_0[i] = g_0_yy_0_yyz_0[i] * pb_x + g_0_yy_0_yyz_1[i] * wp_x[i];

        g_0_xyy_0_yzz_0[i] = g_0_yy_0_yzz_0[i] * pb_x + g_0_yy_0_yzz_1[i] * wp_x[i];

        g_0_xyy_0_zzz_0[i] = g_0_yy_0_zzz_0[i] * pb_x + g_0_yy_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 40-50 components of targeted buffer : SFSF

    auto g_0_xyz_0_xxx_0 = pbuffer.data(idx_eri_0_sfsf + 40);

    auto g_0_xyz_0_xxy_0 = pbuffer.data(idx_eri_0_sfsf + 41);

    auto g_0_xyz_0_xxz_0 = pbuffer.data(idx_eri_0_sfsf + 42);

    auto g_0_xyz_0_xyy_0 = pbuffer.data(idx_eri_0_sfsf + 43);

    auto g_0_xyz_0_xyz_0 = pbuffer.data(idx_eri_0_sfsf + 44);

    auto g_0_xyz_0_xzz_0 = pbuffer.data(idx_eri_0_sfsf + 45);

    auto g_0_xyz_0_yyy_0 = pbuffer.data(idx_eri_0_sfsf + 46);

    auto g_0_xyz_0_yyz_0 = pbuffer.data(idx_eri_0_sfsf + 47);

    auto g_0_xyz_0_yzz_0 = pbuffer.data(idx_eri_0_sfsf + 48);

    auto g_0_xyz_0_zzz_0 = pbuffer.data(idx_eri_0_sfsf + 49);

#pragma omp simd aligned(g_0_xy_0_xxy_0,      \
                             g_0_xy_0_xxy_1,  \
                             g_0_xy_0_xyy_0,  \
                             g_0_xy_0_xyy_1,  \
                             g_0_xyz_0_xxx_0, \
                             g_0_xyz_0_xxy_0, \
                             g_0_xyz_0_xxz_0, \
                             g_0_xyz_0_xyy_0, \
                             g_0_xyz_0_xyz_0, \
                             g_0_xyz_0_xzz_0, \
                             g_0_xyz_0_yyy_0, \
                             g_0_xyz_0_yyz_0, \
                             g_0_xyz_0_yzz_0, \
                             g_0_xyz_0_zzz_0, \
                             g_0_xz_0_xxx_0,  \
                             g_0_xz_0_xxx_1,  \
                             g_0_xz_0_xxz_0,  \
                             g_0_xz_0_xxz_1,  \
                             g_0_xz_0_xzz_0,  \
                             g_0_xz_0_xzz_1,  \
                             g_0_yz_0_xyz_0,  \
                             g_0_yz_0_xyz_1,  \
                             g_0_yz_0_yyy_0,  \
                             g_0_yz_0_yyy_1,  \
                             g_0_yz_0_yyz_0,  \
                             g_0_yz_0_yyz_1,  \
                             g_0_yz_0_yz_1,   \
                             g_0_yz_0_yzz_0,  \
                             g_0_yz_0_yzz_1,  \
                             g_0_yz_0_zzz_0,  \
                             g_0_yz_0_zzz_1,  \
                             wp_x,            \
                             wp_y,            \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xyz_0_xxx_0[i] = g_0_xz_0_xxx_0[i] * pb_y + g_0_xz_0_xxx_1[i] * wp_y[i];

        g_0_xyz_0_xxy_0[i] = g_0_xy_0_xxy_0[i] * pb_z + g_0_xy_0_xxy_1[i] * wp_z[i];

        g_0_xyz_0_xxz_0[i] = g_0_xz_0_xxz_0[i] * pb_y + g_0_xz_0_xxz_1[i] * wp_y[i];

        g_0_xyz_0_xyy_0[i] = g_0_xy_0_xyy_0[i] * pb_z + g_0_xy_0_xyy_1[i] * wp_z[i];

        g_0_xyz_0_xyz_0[i] = g_0_yz_0_yz_1[i] * fi_abcd_0 + g_0_yz_0_xyz_0[i] * pb_x + g_0_yz_0_xyz_1[i] * wp_x[i];

        g_0_xyz_0_xzz_0[i] = g_0_xz_0_xzz_0[i] * pb_y + g_0_xz_0_xzz_1[i] * wp_y[i];

        g_0_xyz_0_yyy_0[i] = g_0_yz_0_yyy_0[i] * pb_x + g_0_yz_0_yyy_1[i] * wp_x[i];

        g_0_xyz_0_yyz_0[i] = g_0_yz_0_yyz_0[i] * pb_x + g_0_yz_0_yyz_1[i] * wp_x[i];

        g_0_xyz_0_yzz_0[i] = g_0_yz_0_yzz_0[i] * pb_x + g_0_yz_0_yzz_1[i] * wp_x[i];

        g_0_xyz_0_zzz_0[i] = g_0_yz_0_zzz_0[i] * pb_x + g_0_yz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 50-60 components of targeted buffer : SFSF

    auto g_0_xzz_0_xxx_0 = pbuffer.data(idx_eri_0_sfsf + 50);

    auto g_0_xzz_0_xxy_0 = pbuffer.data(idx_eri_0_sfsf + 51);

    auto g_0_xzz_0_xxz_0 = pbuffer.data(idx_eri_0_sfsf + 52);

    auto g_0_xzz_0_xyy_0 = pbuffer.data(idx_eri_0_sfsf + 53);

    auto g_0_xzz_0_xyz_0 = pbuffer.data(idx_eri_0_sfsf + 54);

    auto g_0_xzz_0_xzz_0 = pbuffer.data(idx_eri_0_sfsf + 55);

    auto g_0_xzz_0_yyy_0 = pbuffer.data(idx_eri_0_sfsf + 56);

    auto g_0_xzz_0_yyz_0 = pbuffer.data(idx_eri_0_sfsf + 57);

    auto g_0_xzz_0_yzz_0 = pbuffer.data(idx_eri_0_sfsf + 58);

    auto g_0_xzz_0_zzz_0 = pbuffer.data(idx_eri_0_sfsf + 59);

#pragma omp simd aligned(g_0_xzz_0_xxx_0,     \
                             g_0_xzz_0_xxy_0, \
                             g_0_xzz_0_xxz_0, \
                             g_0_xzz_0_xyy_0, \
                             g_0_xzz_0_xyz_0, \
                             g_0_xzz_0_xzz_0, \
                             g_0_xzz_0_yyy_0, \
                             g_0_xzz_0_yyz_0, \
                             g_0_xzz_0_yzz_0, \
                             g_0_xzz_0_zzz_0, \
                             g_0_zz_0_xx_1,   \
                             g_0_zz_0_xxx_0,  \
                             g_0_zz_0_xxx_1,  \
                             g_0_zz_0_xxy_0,  \
                             g_0_zz_0_xxy_1,  \
                             g_0_zz_0_xxz_0,  \
                             g_0_zz_0_xxz_1,  \
                             g_0_zz_0_xy_1,   \
                             g_0_zz_0_xyy_0,  \
                             g_0_zz_0_xyy_1,  \
                             g_0_zz_0_xyz_0,  \
                             g_0_zz_0_xyz_1,  \
                             g_0_zz_0_xz_1,   \
                             g_0_zz_0_xzz_0,  \
                             g_0_zz_0_xzz_1,  \
                             g_0_zz_0_yy_1,   \
                             g_0_zz_0_yyy_0,  \
                             g_0_zz_0_yyy_1,  \
                             g_0_zz_0_yyz_0,  \
                             g_0_zz_0_yyz_1,  \
                             g_0_zz_0_yz_1,   \
                             g_0_zz_0_yzz_0,  \
                             g_0_zz_0_yzz_1,  \
                             g_0_zz_0_zz_1,   \
                             g_0_zz_0_zzz_0,  \
                             g_0_zz_0_zzz_1,  \
                             wp_x,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_xzz_0_xxx_0[i] = 3.0 * g_0_zz_0_xx_1[i] * fi_abcd_0 + g_0_zz_0_xxx_0[i] * pb_x + g_0_zz_0_xxx_1[i] * wp_x[i];

        g_0_xzz_0_xxy_0[i] = 2.0 * g_0_zz_0_xy_1[i] * fi_abcd_0 + g_0_zz_0_xxy_0[i] * pb_x + g_0_zz_0_xxy_1[i] * wp_x[i];

        g_0_xzz_0_xxz_0[i] = 2.0 * g_0_zz_0_xz_1[i] * fi_abcd_0 + g_0_zz_0_xxz_0[i] * pb_x + g_0_zz_0_xxz_1[i] * wp_x[i];

        g_0_xzz_0_xyy_0[i] = g_0_zz_0_yy_1[i] * fi_abcd_0 + g_0_zz_0_xyy_0[i] * pb_x + g_0_zz_0_xyy_1[i] * wp_x[i];

        g_0_xzz_0_xyz_0[i] = g_0_zz_0_yz_1[i] * fi_abcd_0 + g_0_zz_0_xyz_0[i] * pb_x + g_0_zz_0_xyz_1[i] * wp_x[i];

        g_0_xzz_0_xzz_0[i] = g_0_zz_0_zz_1[i] * fi_abcd_0 + g_0_zz_0_xzz_0[i] * pb_x + g_0_zz_0_xzz_1[i] * wp_x[i];

        g_0_xzz_0_yyy_0[i] = g_0_zz_0_yyy_0[i] * pb_x + g_0_zz_0_yyy_1[i] * wp_x[i];

        g_0_xzz_0_yyz_0[i] = g_0_zz_0_yyz_0[i] * pb_x + g_0_zz_0_yyz_1[i] * wp_x[i];

        g_0_xzz_0_yzz_0[i] = g_0_zz_0_yzz_0[i] * pb_x + g_0_zz_0_yzz_1[i] * wp_x[i];

        g_0_xzz_0_zzz_0[i] = g_0_zz_0_zzz_0[i] * pb_x + g_0_zz_0_zzz_1[i] * wp_x[i];
    }

    /// Set up 60-70 components of targeted buffer : SFSF

    auto g_0_yyy_0_xxx_0 = pbuffer.data(idx_eri_0_sfsf + 60);

    auto g_0_yyy_0_xxy_0 = pbuffer.data(idx_eri_0_sfsf + 61);

    auto g_0_yyy_0_xxz_0 = pbuffer.data(idx_eri_0_sfsf + 62);

    auto g_0_yyy_0_xyy_0 = pbuffer.data(idx_eri_0_sfsf + 63);

    auto g_0_yyy_0_xyz_0 = pbuffer.data(idx_eri_0_sfsf + 64);

    auto g_0_yyy_0_xzz_0 = pbuffer.data(idx_eri_0_sfsf + 65);

    auto g_0_yyy_0_yyy_0 = pbuffer.data(idx_eri_0_sfsf + 66);

    auto g_0_yyy_0_yyz_0 = pbuffer.data(idx_eri_0_sfsf + 67);

    auto g_0_yyy_0_yzz_0 = pbuffer.data(idx_eri_0_sfsf + 68);

    auto g_0_yyy_0_zzz_0 = pbuffer.data(idx_eri_0_sfsf + 69);

#pragma omp simd aligned(g_0_y_0_xxx_0,       \
                             g_0_y_0_xxx_1,   \
                             g_0_y_0_xxy_0,   \
                             g_0_y_0_xxy_1,   \
                             g_0_y_0_xxz_0,   \
                             g_0_y_0_xxz_1,   \
                             g_0_y_0_xyy_0,   \
                             g_0_y_0_xyy_1,   \
                             g_0_y_0_xyz_0,   \
                             g_0_y_0_xyz_1,   \
                             g_0_y_0_xzz_0,   \
                             g_0_y_0_xzz_1,   \
                             g_0_y_0_yyy_0,   \
                             g_0_y_0_yyy_1,   \
                             g_0_y_0_yyz_0,   \
                             g_0_y_0_yyz_1,   \
                             g_0_y_0_yzz_0,   \
                             g_0_y_0_yzz_1,   \
                             g_0_y_0_zzz_0,   \
                             g_0_y_0_zzz_1,   \
                             g_0_yy_0_xx_1,   \
                             g_0_yy_0_xxx_0,  \
                             g_0_yy_0_xxx_1,  \
                             g_0_yy_0_xxy_0,  \
                             g_0_yy_0_xxy_1,  \
                             g_0_yy_0_xxz_0,  \
                             g_0_yy_0_xxz_1,  \
                             g_0_yy_0_xy_1,   \
                             g_0_yy_0_xyy_0,  \
                             g_0_yy_0_xyy_1,  \
                             g_0_yy_0_xyz_0,  \
                             g_0_yy_0_xyz_1,  \
                             g_0_yy_0_xz_1,   \
                             g_0_yy_0_xzz_0,  \
                             g_0_yy_0_xzz_1,  \
                             g_0_yy_0_yy_1,   \
                             g_0_yy_0_yyy_0,  \
                             g_0_yy_0_yyy_1,  \
                             g_0_yy_0_yyz_0,  \
                             g_0_yy_0_yyz_1,  \
                             g_0_yy_0_yz_1,   \
                             g_0_yy_0_yzz_0,  \
                             g_0_yy_0_yzz_1,  \
                             g_0_yy_0_zz_1,   \
                             g_0_yy_0_zzz_0,  \
                             g_0_yy_0_zzz_1,  \
                             g_0_yyy_0_xxx_0, \
                             g_0_yyy_0_xxy_0, \
                             g_0_yyy_0_xxz_0, \
                             g_0_yyy_0_xyy_0, \
                             g_0_yyy_0_xyz_0, \
                             g_0_yyy_0_xzz_0, \
                             g_0_yyy_0_yyy_0, \
                             g_0_yyy_0_yyz_0, \
                             g_0_yyy_0_yzz_0, \
                             g_0_yyy_0_zzz_0, \
                             wp_y,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_yyy_0_xxx_0[i] =
            2.0 * g_0_y_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxx_1[i] * fti_ab_0 + g_0_yy_0_xxx_0[i] * pb_y + g_0_yy_0_xxx_1[i] * wp_y[i];

        g_0_yyy_0_xxy_0[i] = 2.0 * g_0_y_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxy_1[i] * fti_ab_0 + g_0_yy_0_xx_1[i] * fi_abcd_0 +
                             g_0_yy_0_xxy_0[i] * pb_y + g_0_yy_0_xxy_1[i] * wp_y[i];

        g_0_yyy_0_xxz_0[i] =
            2.0 * g_0_y_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xxz_1[i] * fti_ab_0 + g_0_yy_0_xxz_0[i] * pb_y + g_0_yy_0_xxz_1[i] * wp_y[i];

        g_0_yyy_0_xyy_0[i] = 2.0 * g_0_y_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyy_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_xy_1[i] * fi_abcd_0 +
                             g_0_yy_0_xyy_0[i] * pb_y + g_0_yy_0_xyy_1[i] * wp_y[i];

        g_0_yyy_0_xyz_0[i] = 2.0 * g_0_y_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xyz_1[i] * fti_ab_0 + g_0_yy_0_xz_1[i] * fi_abcd_0 +
                             g_0_yy_0_xyz_0[i] * pb_y + g_0_yy_0_xyz_1[i] * wp_y[i];

        g_0_yyy_0_xzz_0[i] =
            2.0 * g_0_y_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_xzz_1[i] * fti_ab_0 + g_0_yy_0_xzz_0[i] * pb_y + g_0_yy_0_xzz_1[i] * wp_y[i];

        g_0_yyy_0_yyy_0[i] = 2.0 * g_0_y_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyy_1[i] * fti_ab_0 + 3.0 * g_0_yy_0_yy_1[i] * fi_abcd_0 +
                             g_0_yy_0_yyy_0[i] * pb_y + g_0_yy_0_yyy_1[i] * wp_y[i];

        g_0_yyy_0_yyz_0[i] = 2.0 * g_0_y_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yyz_1[i] * fti_ab_0 + 2.0 * g_0_yy_0_yz_1[i] * fi_abcd_0 +
                             g_0_yy_0_yyz_0[i] * pb_y + g_0_yy_0_yyz_1[i] * wp_y[i];

        g_0_yyy_0_yzz_0[i] = 2.0 * g_0_y_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_yzz_1[i] * fti_ab_0 + g_0_yy_0_zz_1[i] * fi_abcd_0 +
                             g_0_yy_0_yzz_0[i] * pb_y + g_0_yy_0_yzz_1[i] * wp_y[i];

        g_0_yyy_0_zzz_0[i] =
            2.0 * g_0_y_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_y_0_zzz_1[i] * fti_ab_0 + g_0_yy_0_zzz_0[i] * pb_y + g_0_yy_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 70-80 components of targeted buffer : SFSF

    auto g_0_yyz_0_xxx_0 = pbuffer.data(idx_eri_0_sfsf + 70);

    auto g_0_yyz_0_xxy_0 = pbuffer.data(idx_eri_0_sfsf + 71);

    auto g_0_yyz_0_xxz_0 = pbuffer.data(idx_eri_0_sfsf + 72);

    auto g_0_yyz_0_xyy_0 = pbuffer.data(idx_eri_0_sfsf + 73);

    auto g_0_yyz_0_xyz_0 = pbuffer.data(idx_eri_0_sfsf + 74);

    auto g_0_yyz_0_xzz_0 = pbuffer.data(idx_eri_0_sfsf + 75);

    auto g_0_yyz_0_yyy_0 = pbuffer.data(idx_eri_0_sfsf + 76);

    auto g_0_yyz_0_yyz_0 = pbuffer.data(idx_eri_0_sfsf + 77);

    auto g_0_yyz_0_yzz_0 = pbuffer.data(idx_eri_0_sfsf + 78);

    auto g_0_yyz_0_zzz_0 = pbuffer.data(idx_eri_0_sfsf + 79);

#pragma omp simd aligned(g_0_yy_0_xx_1,       \
                             g_0_yy_0_xxx_0,  \
                             g_0_yy_0_xxx_1,  \
                             g_0_yy_0_xxy_0,  \
                             g_0_yy_0_xxy_1,  \
                             g_0_yy_0_xxz_0,  \
                             g_0_yy_0_xxz_1,  \
                             g_0_yy_0_xy_1,   \
                             g_0_yy_0_xyy_0,  \
                             g_0_yy_0_xyy_1,  \
                             g_0_yy_0_xyz_0,  \
                             g_0_yy_0_xyz_1,  \
                             g_0_yy_0_xz_1,   \
                             g_0_yy_0_xzz_0,  \
                             g_0_yy_0_xzz_1,  \
                             g_0_yy_0_yy_1,   \
                             g_0_yy_0_yyy_0,  \
                             g_0_yy_0_yyy_1,  \
                             g_0_yy_0_yyz_0,  \
                             g_0_yy_0_yyz_1,  \
                             g_0_yy_0_yz_1,   \
                             g_0_yy_0_yzz_0,  \
                             g_0_yy_0_yzz_1,  \
                             g_0_yy_0_zz_1,   \
                             g_0_yy_0_zzz_0,  \
                             g_0_yy_0_zzz_1,  \
                             g_0_yyz_0_xxx_0, \
                             g_0_yyz_0_xxy_0, \
                             g_0_yyz_0_xxz_0, \
                             g_0_yyz_0_xyy_0, \
                             g_0_yyz_0_xyz_0, \
                             g_0_yyz_0_xzz_0, \
                             g_0_yyz_0_yyy_0, \
                             g_0_yyz_0_yyz_0, \
                             g_0_yyz_0_yzz_0, \
                             g_0_yyz_0_zzz_0, \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yyz_0_xxx_0[i] = g_0_yy_0_xxx_0[i] * pb_z + g_0_yy_0_xxx_1[i] * wp_z[i];

        g_0_yyz_0_xxy_0[i] = g_0_yy_0_xxy_0[i] * pb_z + g_0_yy_0_xxy_1[i] * wp_z[i];

        g_0_yyz_0_xxz_0[i] = g_0_yy_0_xx_1[i] * fi_abcd_0 + g_0_yy_0_xxz_0[i] * pb_z + g_0_yy_0_xxz_1[i] * wp_z[i];

        g_0_yyz_0_xyy_0[i] = g_0_yy_0_xyy_0[i] * pb_z + g_0_yy_0_xyy_1[i] * wp_z[i];

        g_0_yyz_0_xyz_0[i] = g_0_yy_0_xy_1[i] * fi_abcd_0 + g_0_yy_0_xyz_0[i] * pb_z + g_0_yy_0_xyz_1[i] * wp_z[i];

        g_0_yyz_0_xzz_0[i] = 2.0 * g_0_yy_0_xz_1[i] * fi_abcd_0 + g_0_yy_0_xzz_0[i] * pb_z + g_0_yy_0_xzz_1[i] * wp_z[i];

        g_0_yyz_0_yyy_0[i] = g_0_yy_0_yyy_0[i] * pb_z + g_0_yy_0_yyy_1[i] * wp_z[i];

        g_0_yyz_0_yyz_0[i] = g_0_yy_0_yy_1[i] * fi_abcd_0 + g_0_yy_0_yyz_0[i] * pb_z + g_0_yy_0_yyz_1[i] * wp_z[i];

        g_0_yyz_0_yzz_0[i] = 2.0 * g_0_yy_0_yz_1[i] * fi_abcd_0 + g_0_yy_0_yzz_0[i] * pb_z + g_0_yy_0_yzz_1[i] * wp_z[i];

        g_0_yyz_0_zzz_0[i] = 3.0 * g_0_yy_0_zz_1[i] * fi_abcd_0 + g_0_yy_0_zzz_0[i] * pb_z + g_0_yy_0_zzz_1[i] * wp_z[i];
    }

    /// Set up 80-90 components of targeted buffer : SFSF

    auto g_0_yzz_0_xxx_0 = pbuffer.data(idx_eri_0_sfsf + 80);

    auto g_0_yzz_0_xxy_0 = pbuffer.data(idx_eri_0_sfsf + 81);

    auto g_0_yzz_0_xxz_0 = pbuffer.data(idx_eri_0_sfsf + 82);

    auto g_0_yzz_0_xyy_0 = pbuffer.data(idx_eri_0_sfsf + 83);

    auto g_0_yzz_0_xyz_0 = pbuffer.data(idx_eri_0_sfsf + 84);

    auto g_0_yzz_0_xzz_0 = pbuffer.data(idx_eri_0_sfsf + 85);

    auto g_0_yzz_0_yyy_0 = pbuffer.data(idx_eri_0_sfsf + 86);

    auto g_0_yzz_0_yyz_0 = pbuffer.data(idx_eri_0_sfsf + 87);

    auto g_0_yzz_0_yzz_0 = pbuffer.data(idx_eri_0_sfsf + 88);

    auto g_0_yzz_0_zzz_0 = pbuffer.data(idx_eri_0_sfsf + 89);

#pragma omp simd aligned(g_0_yzz_0_xxx_0,     \
                             g_0_yzz_0_xxy_0, \
                             g_0_yzz_0_xxz_0, \
                             g_0_yzz_0_xyy_0, \
                             g_0_yzz_0_xyz_0, \
                             g_0_yzz_0_xzz_0, \
                             g_0_yzz_0_yyy_0, \
                             g_0_yzz_0_yyz_0, \
                             g_0_yzz_0_yzz_0, \
                             g_0_yzz_0_zzz_0, \
                             g_0_zz_0_xx_1,   \
                             g_0_zz_0_xxx_0,  \
                             g_0_zz_0_xxx_1,  \
                             g_0_zz_0_xxy_0,  \
                             g_0_zz_0_xxy_1,  \
                             g_0_zz_0_xxz_0,  \
                             g_0_zz_0_xxz_1,  \
                             g_0_zz_0_xy_1,   \
                             g_0_zz_0_xyy_0,  \
                             g_0_zz_0_xyy_1,  \
                             g_0_zz_0_xyz_0,  \
                             g_0_zz_0_xyz_1,  \
                             g_0_zz_0_xz_1,   \
                             g_0_zz_0_xzz_0,  \
                             g_0_zz_0_xzz_1,  \
                             g_0_zz_0_yy_1,   \
                             g_0_zz_0_yyy_0,  \
                             g_0_zz_0_yyy_1,  \
                             g_0_zz_0_yyz_0,  \
                             g_0_zz_0_yyz_1,  \
                             g_0_zz_0_yz_1,   \
                             g_0_zz_0_yzz_0,  \
                             g_0_zz_0_yzz_1,  \
                             g_0_zz_0_zz_1,   \
                             g_0_zz_0_zzz_0,  \
                             g_0_zz_0_zzz_1,  \
                             wp_y,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        g_0_yzz_0_xxx_0[i] = g_0_zz_0_xxx_0[i] * pb_y + g_0_zz_0_xxx_1[i] * wp_y[i];

        g_0_yzz_0_xxy_0[i] = g_0_zz_0_xx_1[i] * fi_abcd_0 + g_0_zz_0_xxy_0[i] * pb_y + g_0_zz_0_xxy_1[i] * wp_y[i];

        g_0_yzz_0_xxz_0[i] = g_0_zz_0_xxz_0[i] * pb_y + g_0_zz_0_xxz_1[i] * wp_y[i];

        g_0_yzz_0_xyy_0[i] = 2.0 * g_0_zz_0_xy_1[i] * fi_abcd_0 + g_0_zz_0_xyy_0[i] * pb_y + g_0_zz_0_xyy_1[i] * wp_y[i];

        g_0_yzz_0_xyz_0[i] = g_0_zz_0_xz_1[i] * fi_abcd_0 + g_0_zz_0_xyz_0[i] * pb_y + g_0_zz_0_xyz_1[i] * wp_y[i];

        g_0_yzz_0_xzz_0[i] = g_0_zz_0_xzz_0[i] * pb_y + g_0_zz_0_xzz_1[i] * wp_y[i];

        g_0_yzz_0_yyy_0[i] = 3.0 * g_0_zz_0_yy_1[i] * fi_abcd_0 + g_0_zz_0_yyy_0[i] * pb_y + g_0_zz_0_yyy_1[i] * wp_y[i];

        g_0_yzz_0_yyz_0[i] = 2.0 * g_0_zz_0_yz_1[i] * fi_abcd_0 + g_0_zz_0_yyz_0[i] * pb_y + g_0_zz_0_yyz_1[i] * wp_y[i];

        g_0_yzz_0_yzz_0[i] = g_0_zz_0_zz_1[i] * fi_abcd_0 + g_0_zz_0_yzz_0[i] * pb_y + g_0_zz_0_yzz_1[i] * wp_y[i];

        g_0_yzz_0_zzz_0[i] = g_0_zz_0_zzz_0[i] * pb_y + g_0_zz_0_zzz_1[i] * wp_y[i];
    }

    /// Set up 90-100 components of targeted buffer : SFSF

    auto g_0_zzz_0_xxx_0 = pbuffer.data(idx_eri_0_sfsf + 90);

    auto g_0_zzz_0_xxy_0 = pbuffer.data(idx_eri_0_sfsf + 91);

    auto g_0_zzz_0_xxz_0 = pbuffer.data(idx_eri_0_sfsf + 92);

    auto g_0_zzz_0_xyy_0 = pbuffer.data(idx_eri_0_sfsf + 93);

    auto g_0_zzz_0_xyz_0 = pbuffer.data(idx_eri_0_sfsf + 94);

    auto g_0_zzz_0_xzz_0 = pbuffer.data(idx_eri_0_sfsf + 95);

    auto g_0_zzz_0_yyy_0 = pbuffer.data(idx_eri_0_sfsf + 96);

    auto g_0_zzz_0_yyz_0 = pbuffer.data(idx_eri_0_sfsf + 97);

    auto g_0_zzz_0_yzz_0 = pbuffer.data(idx_eri_0_sfsf + 98);

    auto g_0_zzz_0_zzz_0 = pbuffer.data(idx_eri_0_sfsf + 99);

#pragma omp simd aligned(g_0_z_0_xxx_0,       \
                             g_0_z_0_xxx_1,   \
                             g_0_z_0_xxy_0,   \
                             g_0_z_0_xxy_1,   \
                             g_0_z_0_xxz_0,   \
                             g_0_z_0_xxz_1,   \
                             g_0_z_0_xyy_0,   \
                             g_0_z_0_xyy_1,   \
                             g_0_z_0_xyz_0,   \
                             g_0_z_0_xyz_1,   \
                             g_0_z_0_xzz_0,   \
                             g_0_z_0_xzz_1,   \
                             g_0_z_0_yyy_0,   \
                             g_0_z_0_yyy_1,   \
                             g_0_z_0_yyz_0,   \
                             g_0_z_0_yyz_1,   \
                             g_0_z_0_yzz_0,   \
                             g_0_z_0_yzz_1,   \
                             g_0_z_0_zzz_0,   \
                             g_0_z_0_zzz_1,   \
                             g_0_zz_0_xx_1,   \
                             g_0_zz_0_xxx_0,  \
                             g_0_zz_0_xxx_1,  \
                             g_0_zz_0_xxy_0,  \
                             g_0_zz_0_xxy_1,  \
                             g_0_zz_0_xxz_0,  \
                             g_0_zz_0_xxz_1,  \
                             g_0_zz_0_xy_1,   \
                             g_0_zz_0_xyy_0,  \
                             g_0_zz_0_xyy_1,  \
                             g_0_zz_0_xyz_0,  \
                             g_0_zz_0_xyz_1,  \
                             g_0_zz_0_xz_1,   \
                             g_0_zz_0_xzz_0,  \
                             g_0_zz_0_xzz_1,  \
                             g_0_zz_0_yy_1,   \
                             g_0_zz_0_yyy_0,  \
                             g_0_zz_0_yyy_1,  \
                             g_0_zz_0_yyz_0,  \
                             g_0_zz_0_yyz_1,  \
                             g_0_zz_0_yz_1,   \
                             g_0_zz_0_yzz_0,  \
                             g_0_zz_0_yzz_1,  \
                             g_0_zz_0_zz_1,   \
                             g_0_zz_0_zzz_0,  \
                             g_0_zz_0_zzz_1,  \
                             g_0_zzz_0_xxx_0, \
                             g_0_zzz_0_xxy_0, \
                             g_0_zzz_0_xxz_0, \
                             g_0_zzz_0_xyy_0, \
                             g_0_zzz_0_xyz_0, \
                             g_0_zzz_0_xzz_0, \
                             g_0_zzz_0_yyy_0, \
                             g_0_zzz_0_yyz_0, \
                             g_0_zzz_0_yzz_0, \
                             g_0_zzz_0_zzz_0, \
                             wp_z,            \
                             c_exps,          \
                             d_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fi_ab_0 = 0.5 / (a_exp + b_exp);

        const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);

        const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);

        g_0_zzz_0_xxx_0[i] =
            2.0 * g_0_z_0_xxx_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxx_1[i] * fti_ab_0 + g_0_zz_0_xxx_0[i] * pb_z + g_0_zz_0_xxx_1[i] * wp_z[i];

        g_0_zzz_0_xxy_0[i] =
            2.0 * g_0_z_0_xxy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxy_1[i] * fti_ab_0 + g_0_zz_0_xxy_0[i] * pb_z + g_0_zz_0_xxy_1[i] * wp_z[i];

        g_0_zzz_0_xxz_0[i] = 2.0 * g_0_z_0_xxz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xxz_1[i] * fti_ab_0 + g_0_zz_0_xx_1[i] * fi_abcd_0 +
                             g_0_zz_0_xxz_0[i] * pb_z + g_0_zz_0_xxz_1[i] * wp_z[i];

        g_0_zzz_0_xyy_0[i] =
            2.0 * g_0_z_0_xyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyy_1[i] * fti_ab_0 + g_0_zz_0_xyy_0[i] * pb_z + g_0_zz_0_xyy_1[i] * wp_z[i];

        g_0_zzz_0_xyz_0[i] = 2.0 * g_0_z_0_xyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xyz_1[i] * fti_ab_0 + g_0_zz_0_xy_1[i] * fi_abcd_0 +
                             g_0_zz_0_xyz_0[i] * pb_z + g_0_zz_0_xyz_1[i] * wp_z[i];

        g_0_zzz_0_xzz_0[i] = 2.0 * g_0_z_0_xzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_xzz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_xz_1[i] * fi_abcd_0 +
                             g_0_zz_0_xzz_0[i] * pb_z + g_0_zz_0_xzz_1[i] * wp_z[i];

        g_0_zzz_0_yyy_0[i] =
            2.0 * g_0_z_0_yyy_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyy_1[i] * fti_ab_0 + g_0_zz_0_yyy_0[i] * pb_z + g_0_zz_0_yyy_1[i] * wp_z[i];

        g_0_zzz_0_yyz_0[i] = 2.0 * g_0_z_0_yyz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yyz_1[i] * fti_ab_0 + g_0_zz_0_yy_1[i] * fi_abcd_0 +
                             g_0_zz_0_yyz_0[i] * pb_z + g_0_zz_0_yyz_1[i] * wp_z[i];

        g_0_zzz_0_yzz_0[i] = 2.0 * g_0_z_0_yzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_yzz_1[i] * fti_ab_0 + 2.0 * g_0_zz_0_yz_1[i] * fi_abcd_0 +
                             g_0_zz_0_yzz_0[i] * pb_z + g_0_zz_0_yzz_1[i] * wp_z[i];

        g_0_zzz_0_zzz_0[i] = 2.0 * g_0_z_0_zzz_0[i] * fi_ab_0 - 2.0 * g_0_z_0_zzz_1[i] * fti_ab_0 + 3.0 * g_0_zz_0_zz_1[i] * fi_abcd_0 +
                             g_0_zz_0_zzz_0[i] * pb_z + g_0_zz_0_zzz_1[i] * wp_z[i];
    }
}

}  // namespace erirec
