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

#include "NuclearPotentialGeom020PrimRecDF.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_df(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_020_0_df,
                                        const size_t              idx_npot_geom_020_0_sf,
                                        const size_t              idx_npot_geom_020_1_sf,
                                        const size_t              idx_npot_geom_020_0_pd,
                                        const size_t              idx_npot_geom_020_1_pd,
                                        const size_t              idx_npot_geom_010_1_pf,
                                        const size_t              idx_npot_geom_020_0_pf,
                                        const size_t              idx_npot_geom_020_1_pf,
                                        const CSimdArray<double>& factors,
                                        const size_t              idx_rpa,
                                        const size_t              idx_rpc,
                                        const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up R(PC) distances

    auto pc_x = factors.data(idx_rpc);

    auto pc_y = factors.data(idx_rpc + 1);

    auto pc_z = factors.data(idx_rpc + 2);

    // Set up components of auxiliary buffer : SF

    auto ta2_xx_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf);

    auto ta2_xx_0_xxy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 1);

    auto ta2_xx_0_xxz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 2);

    auto ta2_xx_0_xyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 3);

    auto ta2_xx_0_xyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 4);

    auto ta2_xx_0_xzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 5);

    auto ta2_xx_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 6);

    auto ta2_xx_0_yyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 7);

    auto ta2_xx_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 8);

    auto ta2_xx_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 9);

    auto ta2_xy_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf + 10);

    auto ta2_xy_0_xxy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 11);

    auto ta2_xy_0_xxz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 12);

    auto ta2_xy_0_xyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 13);

    auto ta2_xy_0_xyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 14);

    auto ta2_xy_0_xzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 15);

    auto ta2_xy_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 16);

    auto ta2_xy_0_yyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 17);

    auto ta2_xy_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 18);

    auto ta2_xy_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 19);

    auto ta2_xz_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf + 20);

    auto ta2_xz_0_xxy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 21);

    auto ta2_xz_0_xxz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 22);

    auto ta2_xz_0_xyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 23);

    auto ta2_xz_0_xyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 24);

    auto ta2_xz_0_xzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 25);

    auto ta2_xz_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 26);

    auto ta2_xz_0_yyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 27);

    auto ta2_xz_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 28);

    auto ta2_xz_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 29);

    auto ta2_yy_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf + 30);

    auto ta2_yy_0_xxy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 31);

    auto ta2_yy_0_xxz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 32);

    auto ta2_yy_0_xyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 33);

    auto ta2_yy_0_xyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 34);

    auto ta2_yy_0_xzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 35);

    auto ta2_yy_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 36);

    auto ta2_yy_0_yyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 37);

    auto ta2_yy_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 38);

    auto ta2_yy_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 39);

    auto ta2_yz_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf + 40);

    auto ta2_yz_0_xxy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 41);

    auto ta2_yz_0_xxz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 42);

    auto ta2_yz_0_xyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 43);

    auto ta2_yz_0_xyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 44);

    auto ta2_yz_0_xzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 45);

    auto ta2_yz_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 46);

    auto ta2_yz_0_yyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 47);

    auto ta2_yz_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 48);

    auto ta2_yz_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 49);

    auto ta2_zz_0_xxx_0 = pbuffer.data(idx_npot_geom_020_0_sf + 50);

    auto ta2_zz_0_xxy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 51);

    auto ta2_zz_0_xxz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 52);

    auto ta2_zz_0_xyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 53);

    auto ta2_zz_0_xyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 54);

    auto ta2_zz_0_xzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 55);

    auto ta2_zz_0_yyy_0 = pbuffer.data(idx_npot_geom_020_0_sf + 56);

    auto ta2_zz_0_yyz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 57);

    auto ta2_zz_0_yzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 58);

    auto ta2_zz_0_zzz_0 = pbuffer.data(idx_npot_geom_020_0_sf + 59);

    // Set up components of auxiliary buffer : SF

    auto ta2_xx_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf);

    auto ta2_xx_0_xxy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 1);

    auto ta2_xx_0_xxz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 2);

    auto ta2_xx_0_xyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 3);

    auto ta2_xx_0_xyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 4);

    auto ta2_xx_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 5);

    auto ta2_xx_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 6);

    auto ta2_xx_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 7);

    auto ta2_xx_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 8);

    auto ta2_xx_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 9);

    auto ta2_xy_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf + 10);

    auto ta2_xy_0_xxy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 11);

    auto ta2_xy_0_xxz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 12);

    auto ta2_xy_0_xyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 13);

    auto ta2_xy_0_xyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 14);

    auto ta2_xy_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 15);

    auto ta2_xy_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 16);

    auto ta2_xy_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 17);

    auto ta2_xy_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 18);

    auto ta2_xy_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 19);

    auto ta2_xz_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf + 20);

    auto ta2_xz_0_xxy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 21);

    auto ta2_xz_0_xxz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 22);

    auto ta2_xz_0_xyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 23);

    auto ta2_xz_0_xyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 24);

    auto ta2_xz_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 25);

    auto ta2_xz_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 26);

    auto ta2_xz_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 27);

    auto ta2_xz_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 28);

    auto ta2_xz_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 29);

    auto ta2_yy_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf + 30);

    auto ta2_yy_0_xxy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 31);

    auto ta2_yy_0_xxz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 32);

    auto ta2_yy_0_xyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 33);

    auto ta2_yy_0_xyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 34);

    auto ta2_yy_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 35);

    auto ta2_yy_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 36);

    auto ta2_yy_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 37);

    auto ta2_yy_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 38);

    auto ta2_yy_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 39);

    auto ta2_yz_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf + 40);

    auto ta2_yz_0_xxy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 41);

    auto ta2_yz_0_xxz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 42);

    auto ta2_yz_0_xyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 43);

    auto ta2_yz_0_xyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 44);

    auto ta2_yz_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 45);

    auto ta2_yz_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 46);

    auto ta2_yz_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 47);

    auto ta2_yz_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 48);

    auto ta2_yz_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 49);

    auto ta2_zz_0_xxx_1 = pbuffer.data(idx_npot_geom_020_1_sf + 50);

    auto ta2_zz_0_xxy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 51);

    auto ta2_zz_0_xxz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 52);

    auto ta2_zz_0_xyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 53);

    auto ta2_zz_0_xyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 54);

    auto ta2_zz_0_xzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 55);

    auto ta2_zz_0_yyy_1 = pbuffer.data(idx_npot_geom_020_1_sf + 56);

    auto ta2_zz_0_yyz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 57);

    auto ta2_zz_0_yzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 58);

    auto ta2_zz_0_zzz_1 = pbuffer.data(idx_npot_geom_020_1_sf + 59);

    // Set up components of auxiliary buffer : PD

    auto ta2_xx_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd);

    auto ta2_xx_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 1);

    auto ta2_xx_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 2);

    auto ta2_xx_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 3);

    auto ta2_xx_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 4);

    auto ta2_xx_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 5);

    auto ta2_xx_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 6);

    auto ta2_xx_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 7);

    auto ta2_xx_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 8);

    auto ta2_xx_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 9);

    auto ta2_xx_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 10);

    auto ta2_xx_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 11);

    auto ta2_xx_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 12);

    auto ta2_xx_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 13);

    auto ta2_xx_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 14);

    auto ta2_xx_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 15);

    auto ta2_xx_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 16);

    auto ta2_xx_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 17);

    auto ta2_xy_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 18);

    auto ta2_xy_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 19);

    auto ta2_xy_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 20);

    auto ta2_xy_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 21);

    auto ta2_xy_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 22);

    auto ta2_xy_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 23);

    auto ta2_xy_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 24);

    auto ta2_xy_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 25);

    auto ta2_xy_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 26);

    auto ta2_xy_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 27);

    auto ta2_xy_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 28);

    auto ta2_xy_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 29);

    auto ta2_xy_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 30);

    auto ta2_xy_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 31);

    auto ta2_xy_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 32);

    auto ta2_xy_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 33);

    auto ta2_xy_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 34);

    auto ta2_xy_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 35);

    auto ta2_xz_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 36);

    auto ta2_xz_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 37);

    auto ta2_xz_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 38);

    auto ta2_xz_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 39);

    auto ta2_xz_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 40);

    auto ta2_xz_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 41);

    auto ta2_xz_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 42);

    auto ta2_xz_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 43);

    auto ta2_xz_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 44);

    auto ta2_xz_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 45);

    auto ta2_xz_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 46);

    auto ta2_xz_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 47);

    auto ta2_xz_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 48);

    auto ta2_xz_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 49);

    auto ta2_xz_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 50);

    auto ta2_xz_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 51);

    auto ta2_xz_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 52);

    auto ta2_xz_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 53);

    auto ta2_yy_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 54);

    auto ta2_yy_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 55);

    auto ta2_yy_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 56);

    auto ta2_yy_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 57);

    auto ta2_yy_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 58);

    auto ta2_yy_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 59);

    auto ta2_yy_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 60);

    auto ta2_yy_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 61);

    auto ta2_yy_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 62);

    auto ta2_yy_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 63);

    auto ta2_yy_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 64);

    auto ta2_yy_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 65);

    auto ta2_yy_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 66);

    auto ta2_yy_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 67);

    auto ta2_yy_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 68);

    auto ta2_yy_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 69);

    auto ta2_yy_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 70);

    auto ta2_yy_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 71);

    auto ta2_yz_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 72);

    auto ta2_yz_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 73);

    auto ta2_yz_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 74);

    auto ta2_yz_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 75);

    auto ta2_yz_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 76);

    auto ta2_yz_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 77);

    auto ta2_yz_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 78);

    auto ta2_yz_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 79);

    auto ta2_yz_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 80);

    auto ta2_yz_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 81);

    auto ta2_yz_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 82);

    auto ta2_yz_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 83);

    auto ta2_yz_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 84);

    auto ta2_yz_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 85);

    auto ta2_yz_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 86);

    auto ta2_yz_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 87);

    auto ta2_yz_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 88);

    auto ta2_yz_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 89);

    auto ta2_zz_x_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 90);

    auto ta2_zz_x_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 91);

    auto ta2_zz_x_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 92);

    auto ta2_zz_x_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 93);

    auto ta2_zz_x_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 94);

    auto ta2_zz_x_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 95);

    auto ta2_zz_y_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 96);

    auto ta2_zz_y_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 97);

    auto ta2_zz_y_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 98);

    auto ta2_zz_y_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 99);

    auto ta2_zz_y_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 100);

    auto ta2_zz_y_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 101);

    auto ta2_zz_z_xx_0 = pbuffer.data(idx_npot_geom_020_0_pd + 102);

    auto ta2_zz_z_xy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 103);

    auto ta2_zz_z_xz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 104);

    auto ta2_zz_z_yy_0 = pbuffer.data(idx_npot_geom_020_0_pd + 105);

    auto ta2_zz_z_yz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 106);

    auto ta2_zz_z_zz_0 = pbuffer.data(idx_npot_geom_020_0_pd + 107);

    // Set up components of auxiliary buffer : PD

    auto ta2_xx_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd);

    auto ta2_xx_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 1);

    auto ta2_xx_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 2);

    auto ta2_xx_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 3);

    auto ta2_xx_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 4);

    auto ta2_xx_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 5);

    auto ta2_xx_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 6);

    auto ta2_xx_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 7);

    auto ta2_xx_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 8);

    auto ta2_xx_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 9);

    auto ta2_xx_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 10);

    auto ta2_xx_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 11);

    auto ta2_xx_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 12);

    auto ta2_xx_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 13);

    auto ta2_xx_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 14);

    auto ta2_xx_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 15);

    auto ta2_xx_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 16);

    auto ta2_xx_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 17);

    auto ta2_xy_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 18);

    auto ta2_xy_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 19);

    auto ta2_xy_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 20);

    auto ta2_xy_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 21);

    auto ta2_xy_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 22);

    auto ta2_xy_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 23);

    auto ta2_xy_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 24);

    auto ta2_xy_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 25);

    auto ta2_xy_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 26);

    auto ta2_xy_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 27);

    auto ta2_xy_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 28);

    auto ta2_xy_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 29);

    auto ta2_xy_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 30);

    auto ta2_xy_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 31);

    auto ta2_xy_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 32);

    auto ta2_xy_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 33);

    auto ta2_xy_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 34);

    auto ta2_xy_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 35);

    auto ta2_xz_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 36);

    auto ta2_xz_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 37);

    auto ta2_xz_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 38);

    auto ta2_xz_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 39);

    auto ta2_xz_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 40);

    auto ta2_xz_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 41);

    auto ta2_xz_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 42);

    auto ta2_xz_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 43);

    auto ta2_xz_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 44);

    auto ta2_xz_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 45);

    auto ta2_xz_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 46);

    auto ta2_xz_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 47);

    auto ta2_xz_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 48);

    auto ta2_xz_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 49);

    auto ta2_xz_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 50);

    auto ta2_xz_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 51);

    auto ta2_xz_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 52);

    auto ta2_xz_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 53);

    auto ta2_yy_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 54);

    auto ta2_yy_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 55);

    auto ta2_yy_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 56);

    auto ta2_yy_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 57);

    auto ta2_yy_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 58);

    auto ta2_yy_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 59);

    auto ta2_yy_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 60);

    auto ta2_yy_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 61);

    auto ta2_yy_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 62);

    auto ta2_yy_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 63);

    auto ta2_yy_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 64);

    auto ta2_yy_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 65);

    auto ta2_yy_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 66);

    auto ta2_yy_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 67);

    auto ta2_yy_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 68);

    auto ta2_yy_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 69);

    auto ta2_yy_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 70);

    auto ta2_yy_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 71);

    auto ta2_yz_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 72);

    auto ta2_yz_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 73);

    auto ta2_yz_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 74);

    auto ta2_yz_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 75);

    auto ta2_yz_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 76);

    auto ta2_yz_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 77);

    auto ta2_yz_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 78);

    auto ta2_yz_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 79);

    auto ta2_yz_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 80);

    auto ta2_yz_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 81);

    auto ta2_yz_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 82);

    auto ta2_yz_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 83);

    auto ta2_yz_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 84);

    auto ta2_yz_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 85);

    auto ta2_yz_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 86);

    auto ta2_yz_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 87);

    auto ta2_yz_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 88);

    auto ta2_yz_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 89);

    auto ta2_zz_x_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 90);

    auto ta2_zz_x_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 91);

    auto ta2_zz_x_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 92);

    auto ta2_zz_x_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 93);

    auto ta2_zz_x_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 94);

    auto ta2_zz_x_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 95);

    auto ta2_zz_y_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 96);

    auto ta2_zz_y_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 97);

    auto ta2_zz_y_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 98);

    auto ta2_zz_y_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 99);

    auto ta2_zz_y_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 100);

    auto ta2_zz_y_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 101);

    auto ta2_zz_z_xx_1 = pbuffer.data(idx_npot_geom_020_1_pd + 102);

    auto ta2_zz_z_xy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 103);

    auto ta2_zz_z_xz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 104);

    auto ta2_zz_z_yy_1 = pbuffer.data(idx_npot_geom_020_1_pd + 105);

    auto ta2_zz_z_yz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 106);

    auto ta2_zz_z_zz_1 = pbuffer.data(idx_npot_geom_020_1_pd + 107);

    // Set up components of auxiliary buffer : PF

    auto ta1_x_x_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf);

    auto ta1_x_x_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 1);

    auto ta1_x_x_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 2);

    auto ta1_x_x_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 3);

    auto ta1_x_x_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 4);

    auto ta1_x_x_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 5);

    auto ta1_x_x_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 6);

    auto ta1_x_x_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 7);

    auto ta1_x_x_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 8);

    auto ta1_x_x_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 9);

    auto ta1_x_y_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 10);

    auto ta1_x_y_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 11);

    auto ta1_x_y_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 12);

    auto ta1_x_y_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 13);

    auto ta1_x_y_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 14);

    auto ta1_x_y_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 15);

    auto ta1_x_y_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 16);

    auto ta1_x_y_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 17);

    auto ta1_x_y_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 18);

    auto ta1_x_y_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 19);

    auto ta1_x_z_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 20);

    auto ta1_x_z_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 21);

    auto ta1_x_z_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 22);

    auto ta1_x_z_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 23);

    auto ta1_x_z_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 24);

    auto ta1_x_z_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 25);

    auto ta1_x_z_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 26);

    auto ta1_x_z_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 27);

    auto ta1_x_z_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 28);

    auto ta1_x_z_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 29);

    auto ta1_y_x_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 30);

    auto ta1_y_x_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 31);

    auto ta1_y_x_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 32);

    auto ta1_y_x_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 33);

    auto ta1_y_x_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 34);

    auto ta1_y_x_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 35);

    auto ta1_y_x_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 36);

    auto ta1_y_x_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 37);

    auto ta1_y_x_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 38);

    auto ta1_y_x_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 39);

    auto ta1_y_y_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 40);

    auto ta1_y_y_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 41);

    auto ta1_y_y_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 42);

    auto ta1_y_y_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 43);

    auto ta1_y_y_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 44);

    auto ta1_y_y_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 45);

    auto ta1_y_y_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 46);

    auto ta1_y_y_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 47);

    auto ta1_y_y_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 48);

    auto ta1_y_y_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 49);

    auto ta1_y_z_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 50);

    auto ta1_y_z_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 51);

    auto ta1_y_z_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 52);

    auto ta1_y_z_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 53);

    auto ta1_y_z_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 54);

    auto ta1_y_z_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 55);

    auto ta1_y_z_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 56);

    auto ta1_y_z_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 57);

    auto ta1_y_z_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 58);

    auto ta1_y_z_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 59);

    auto ta1_z_x_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 60);

    auto ta1_z_x_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 61);

    auto ta1_z_x_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 62);

    auto ta1_z_x_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 63);

    auto ta1_z_x_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 64);

    auto ta1_z_x_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 65);

    auto ta1_z_x_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 66);

    auto ta1_z_x_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 67);

    auto ta1_z_x_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 68);

    auto ta1_z_x_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 69);

    auto ta1_z_y_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 70);

    auto ta1_z_y_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 71);

    auto ta1_z_y_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 72);

    auto ta1_z_y_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 73);

    auto ta1_z_y_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 74);

    auto ta1_z_y_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 75);

    auto ta1_z_y_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 76);

    auto ta1_z_y_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 77);

    auto ta1_z_y_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 78);

    auto ta1_z_y_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 79);

    auto ta1_z_z_xxx_1 = pbuffer.data(idx_npot_geom_010_1_pf + 80);

    auto ta1_z_z_xxy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 81);

    auto ta1_z_z_xxz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 82);

    auto ta1_z_z_xyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 83);

    auto ta1_z_z_xyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 84);

    auto ta1_z_z_xzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 85);

    auto ta1_z_z_yyy_1 = pbuffer.data(idx_npot_geom_010_1_pf + 86);

    auto ta1_z_z_yyz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 87);

    auto ta1_z_z_yzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 88);

    auto ta1_z_z_zzz_1 = pbuffer.data(idx_npot_geom_010_1_pf + 89);

    // Set up components of auxiliary buffer : PF

    auto ta2_xx_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf);

    auto ta2_xx_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 1);

    auto ta2_xx_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 2);

    auto ta2_xx_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 3);

    auto ta2_xx_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 4);

    auto ta2_xx_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 5);

    auto ta2_xx_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 6);

    auto ta2_xx_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 7);

    auto ta2_xx_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 8);

    auto ta2_xx_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 9);

    auto ta2_xx_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 10);

    auto ta2_xx_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 11);

    auto ta2_xx_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 12);

    auto ta2_xx_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 13);

    auto ta2_xx_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 14);

    auto ta2_xx_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 15);

    auto ta2_xx_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 16);

    auto ta2_xx_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 17);

    auto ta2_xx_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 18);

    auto ta2_xx_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 19);

    auto ta2_xx_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 20);

    auto ta2_xx_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 21);

    auto ta2_xx_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 22);

    auto ta2_xx_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 23);

    auto ta2_xx_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 24);

    auto ta2_xx_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 25);

    auto ta2_xx_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 26);

    auto ta2_xx_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 27);

    auto ta2_xx_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 28);

    auto ta2_xx_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 29);

    auto ta2_xy_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 30);

    auto ta2_xy_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 31);

    auto ta2_xy_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 32);

    auto ta2_xy_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 33);

    auto ta2_xy_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 34);

    auto ta2_xy_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 35);

    auto ta2_xy_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 36);

    auto ta2_xy_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 37);

    auto ta2_xy_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 38);

    auto ta2_xy_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 39);

    auto ta2_xy_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 40);

    auto ta2_xy_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 41);

    auto ta2_xy_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 42);

    auto ta2_xy_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 43);

    auto ta2_xy_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 44);

    auto ta2_xy_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 45);

    auto ta2_xy_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 46);

    auto ta2_xy_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 47);

    auto ta2_xy_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 48);

    auto ta2_xy_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 49);

    auto ta2_xy_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 50);

    auto ta2_xy_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 51);

    auto ta2_xy_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 52);

    auto ta2_xy_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 53);

    auto ta2_xy_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 54);

    auto ta2_xy_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 55);

    auto ta2_xy_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 56);

    auto ta2_xy_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 57);

    auto ta2_xy_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 58);

    auto ta2_xy_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 59);

    auto ta2_xz_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 60);

    auto ta2_xz_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 61);

    auto ta2_xz_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 62);

    auto ta2_xz_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 63);

    auto ta2_xz_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 64);

    auto ta2_xz_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 65);

    auto ta2_xz_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 66);

    auto ta2_xz_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 67);

    auto ta2_xz_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 68);

    auto ta2_xz_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 69);

    auto ta2_xz_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 70);

    auto ta2_xz_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 71);

    auto ta2_xz_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 72);

    auto ta2_xz_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 73);

    auto ta2_xz_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 74);

    auto ta2_xz_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 75);

    auto ta2_xz_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 76);

    auto ta2_xz_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 77);

    auto ta2_xz_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 78);

    auto ta2_xz_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 79);

    auto ta2_xz_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 80);

    auto ta2_xz_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 81);

    auto ta2_xz_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 82);

    auto ta2_xz_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 83);

    auto ta2_xz_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 84);

    auto ta2_xz_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 85);

    auto ta2_xz_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 86);

    auto ta2_xz_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 87);

    auto ta2_xz_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 88);

    auto ta2_xz_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 89);

    auto ta2_yy_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 90);

    auto ta2_yy_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 91);

    auto ta2_yy_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 92);

    auto ta2_yy_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 93);

    auto ta2_yy_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 94);

    auto ta2_yy_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 95);

    auto ta2_yy_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 96);

    auto ta2_yy_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 97);

    auto ta2_yy_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 98);

    auto ta2_yy_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 99);

    auto ta2_yy_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 100);

    auto ta2_yy_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 101);

    auto ta2_yy_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 102);

    auto ta2_yy_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 103);

    auto ta2_yy_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 104);

    auto ta2_yy_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 105);

    auto ta2_yy_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 106);

    auto ta2_yy_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 107);

    auto ta2_yy_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 108);

    auto ta2_yy_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 109);

    auto ta2_yy_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 110);

    auto ta2_yy_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 111);

    auto ta2_yy_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 112);

    auto ta2_yy_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 113);

    auto ta2_yy_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 114);

    auto ta2_yy_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 115);

    auto ta2_yy_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 116);

    auto ta2_yy_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 117);

    auto ta2_yy_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 118);

    auto ta2_yy_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 119);

    auto ta2_yz_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 120);

    auto ta2_yz_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 121);

    auto ta2_yz_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 122);

    auto ta2_yz_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 123);

    auto ta2_yz_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 124);

    auto ta2_yz_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 125);

    auto ta2_yz_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 126);

    auto ta2_yz_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 127);

    auto ta2_yz_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 128);

    auto ta2_yz_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 129);

    auto ta2_yz_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 130);

    auto ta2_yz_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 131);

    auto ta2_yz_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 132);

    auto ta2_yz_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 133);

    auto ta2_yz_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 134);

    auto ta2_yz_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 135);

    auto ta2_yz_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 136);

    auto ta2_yz_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 137);

    auto ta2_yz_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 138);

    auto ta2_yz_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 139);

    auto ta2_yz_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 140);

    auto ta2_yz_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 141);

    auto ta2_yz_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 142);

    auto ta2_yz_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 143);

    auto ta2_yz_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 144);

    auto ta2_yz_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 145);

    auto ta2_yz_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 146);

    auto ta2_yz_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 147);

    auto ta2_yz_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 148);

    auto ta2_yz_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 149);

    auto ta2_zz_x_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 150);

    auto ta2_zz_x_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 151);

    auto ta2_zz_x_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 152);

    auto ta2_zz_x_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 153);

    auto ta2_zz_x_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 154);

    auto ta2_zz_x_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 155);

    auto ta2_zz_x_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 156);

    auto ta2_zz_x_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 157);

    auto ta2_zz_x_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 158);

    auto ta2_zz_x_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 159);

    auto ta2_zz_y_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 160);

    auto ta2_zz_y_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 161);

    auto ta2_zz_y_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 162);

    auto ta2_zz_y_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 163);

    auto ta2_zz_y_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 164);

    auto ta2_zz_y_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 165);

    auto ta2_zz_y_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 166);

    auto ta2_zz_y_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 167);

    auto ta2_zz_y_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 168);

    auto ta2_zz_y_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 169);

    auto ta2_zz_z_xxx_0 = pbuffer.data(idx_npot_geom_020_0_pf + 170);

    auto ta2_zz_z_xxy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 171);

    auto ta2_zz_z_xxz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 172);

    auto ta2_zz_z_xyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 173);

    auto ta2_zz_z_xyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 174);

    auto ta2_zz_z_xzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 175);

    auto ta2_zz_z_yyy_0 = pbuffer.data(idx_npot_geom_020_0_pf + 176);

    auto ta2_zz_z_yyz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 177);

    auto ta2_zz_z_yzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 178);

    auto ta2_zz_z_zzz_0 = pbuffer.data(idx_npot_geom_020_0_pf + 179);

    // Set up components of auxiliary buffer : PF

    auto ta2_xx_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf);

    auto ta2_xx_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 1);

    auto ta2_xx_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 2);

    auto ta2_xx_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 3);

    auto ta2_xx_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 4);

    auto ta2_xx_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 5);

    auto ta2_xx_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 6);

    auto ta2_xx_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 7);

    auto ta2_xx_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 8);

    auto ta2_xx_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 9);

    auto ta2_xx_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 10);

    auto ta2_xx_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 11);

    auto ta2_xx_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 12);

    auto ta2_xx_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 13);

    auto ta2_xx_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 14);

    auto ta2_xx_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 15);

    auto ta2_xx_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 16);

    auto ta2_xx_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 17);

    auto ta2_xx_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 18);

    auto ta2_xx_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 19);

    auto ta2_xx_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 20);

    auto ta2_xx_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 21);

    auto ta2_xx_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 22);

    auto ta2_xx_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 23);

    auto ta2_xx_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 24);

    auto ta2_xx_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 25);

    auto ta2_xx_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 26);

    auto ta2_xx_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 27);

    auto ta2_xx_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 28);

    auto ta2_xx_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 29);

    auto ta2_xy_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 30);

    auto ta2_xy_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 31);

    auto ta2_xy_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 32);

    auto ta2_xy_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 33);

    auto ta2_xy_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 34);

    auto ta2_xy_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 35);

    auto ta2_xy_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 36);

    auto ta2_xy_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 37);

    auto ta2_xy_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 38);

    auto ta2_xy_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 39);

    auto ta2_xy_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 40);

    auto ta2_xy_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 41);

    auto ta2_xy_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 42);

    auto ta2_xy_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 43);

    auto ta2_xy_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 44);

    auto ta2_xy_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 45);

    auto ta2_xy_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 46);

    auto ta2_xy_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 47);

    auto ta2_xy_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 48);

    auto ta2_xy_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 49);

    auto ta2_xy_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 50);

    auto ta2_xy_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 51);

    auto ta2_xy_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 52);

    auto ta2_xy_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 53);

    auto ta2_xy_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 54);

    auto ta2_xy_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 55);

    auto ta2_xy_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 56);

    auto ta2_xy_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 57);

    auto ta2_xy_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 58);

    auto ta2_xy_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 59);

    auto ta2_xz_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 60);

    auto ta2_xz_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 61);

    auto ta2_xz_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 62);

    auto ta2_xz_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 63);

    auto ta2_xz_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 64);

    auto ta2_xz_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 65);

    auto ta2_xz_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 66);

    auto ta2_xz_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 67);

    auto ta2_xz_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 68);

    auto ta2_xz_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 69);

    auto ta2_xz_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 70);

    auto ta2_xz_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 71);

    auto ta2_xz_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 72);

    auto ta2_xz_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 73);

    auto ta2_xz_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 74);

    auto ta2_xz_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 75);

    auto ta2_xz_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 76);

    auto ta2_xz_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 77);

    auto ta2_xz_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 78);

    auto ta2_xz_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 79);

    auto ta2_xz_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 80);

    auto ta2_xz_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 81);

    auto ta2_xz_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 82);

    auto ta2_xz_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 83);

    auto ta2_xz_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 84);

    auto ta2_xz_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 85);

    auto ta2_xz_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 86);

    auto ta2_xz_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 87);

    auto ta2_xz_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 88);

    auto ta2_xz_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 89);

    auto ta2_yy_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 90);

    auto ta2_yy_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 91);

    auto ta2_yy_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 92);

    auto ta2_yy_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 93);

    auto ta2_yy_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 94);

    auto ta2_yy_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 95);

    auto ta2_yy_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 96);

    auto ta2_yy_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 97);

    auto ta2_yy_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 98);

    auto ta2_yy_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 99);

    auto ta2_yy_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 100);

    auto ta2_yy_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 101);

    auto ta2_yy_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 102);

    auto ta2_yy_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 103);

    auto ta2_yy_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 104);

    auto ta2_yy_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 105);

    auto ta2_yy_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 106);

    auto ta2_yy_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 107);

    auto ta2_yy_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 108);

    auto ta2_yy_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 109);

    auto ta2_yy_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 110);

    auto ta2_yy_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 111);

    auto ta2_yy_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 112);

    auto ta2_yy_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 113);

    auto ta2_yy_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 114);

    auto ta2_yy_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 115);

    auto ta2_yy_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 116);

    auto ta2_yy_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 117);

    auto ta2_yy_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 118);

    auto ta2_yy_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 119);

    auto ta2_yz_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 120);

    auto ta2_yz_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 121);

    auto ta2_yz_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 122);

    auto ta2_yz_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 123);

    auto ta2_yz_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 124);

    auto ta2_yz_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 125);

    auto ta2_yz_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 126);

    auto ta2_yz_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 127);

    auto ta2_yz_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 128);

    auto ta2_yz_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 129);

    auto ta2_yz_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 130);

    auto ta2_yz_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 131);

    auto ta2_yz_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 132);

    auto ta2_yz_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 133);

    auto ta2_yz_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 134);

    auto ta2_yz_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 135);

    auto ta2_yz_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 136);

    auto ta2_yz_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 137);

    auto ta2_yz_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 138);

    auto ta2_yz_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 139);

    auto ta2_yz_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 140);

    auto ta2_yz_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 141);

    auto ta2_yz_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 142);

    auto ta2_yz_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 143);

    auto ta2_yz_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 144);

    auto ta2_yz_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 145);

    auto ta2_yz_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 146);

    auto ta2_yz_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 147);

    auto ta2_yz_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 148);

    auto ta2_yz_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 149);

    auto ta2_zz_x_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 150);

    auto ta2_zz_x_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 151);

    auto ta2_zz_x_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 152);

    auto ta2_zz_x_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 153);

    auto ta2_zz_x_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 154);

    auto ta2_zz_x_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 155);

    auto ta2_zz_x_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 156);

    auto ta2_zz_x_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 157);

    auto ta2_zz_x_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 158);

    auto ta2_zz_x_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 159);

    auto ta2_zz_y_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 160);

    auto ta2_zz_y_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 161);

    auto ta2_zz_y_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 162);

    auto ta2_zz_y_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 163);

    auto ta2_zz_y_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 164);

    auto ta2_zz_y_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 165);

    auto ta2_zz_y_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 166);

    auto ta2_zz_y_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 167);

    auto ta2_zz_y_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 168);

    auto ta2_zz_y_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 169);

    auto ta2_zz_z_xxx_1 = pbuffer.data(idx_npot_geom_020_1_pf + 170);

    auto ta2_zz_z_xxy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 171);

    auto ta2_zz_z_xxz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 172);

    auto ta2_zz_z_xyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 173);

    auto ta2_zz_z_xyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 174);

    auto ta2_zz_z_xzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 175);

    auto ta2_zz_z_yyy_1 = pbuffer.data(idx_npot_geom_020_1_pf + 176);

    auto ta2_zz_z_yyz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 177);

    auto ta2_zz_z_yzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 178);

    auto ta2_zz_z_zzz_1 = pbuffer.data(idx_npot_geom_020_1_pf + 179);

    // Set up 0-10 components of targeted buffer : DF

    auto ta2_xx_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df);

    auto ta2_xx_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 1);

    auto ta2_xx_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 2);

    auto ta2_xx_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 3);

    auto ta2_xx_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 4);

    auto ta2_xx_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 5);

    auto ta2_xx_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 6);

    auto ta2_xx_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 7);

    auto ta2_xx_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 8);

    auto ta2_xx_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 9);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_x_x_xxx_1,   \
                             ta1_x_x_xxy_1,   \
                             ta1_x_x_xxz_1,   \
                             ta1_x_x_xyy_1,   \
                             ta1_x_x_xyz_1,   \
                             ta1_x_x_xzz_1,   \
                             ta1_x_x_yyy_1,   \
                             ta1_x_x_yyz_1,   \
                             ta1_x_x_yzz_1,   \
                             ta1_x_x_zzz_1,   \
                             ta2_xx_0_xxx_0,  \
                             ta2_xx_0_xxx_1,  \
                             ta2_xx_0_xxy_0,  \
                             ta2_xx_0_xxy_1,  \
                             ta2_xx_0_xxz_0,  \
                             ta2_xx_0_xxz_1,  \
                             ta2_xx_0_xyy_0,  \
                             ta2_xx_0_xyy_1,  \
                             ta2_xx_0_xyz_0,  \
                             ta2_xx_0_xyz_1,  \
                             ta2_xx_0_xzz_0,  \
                             ta2_xx_0_xzz_1,  \
                             ta2_xx_0_yyy_0,  \
                             ta2_xx_0_yyy_1,  \
                             ta2_xx_0_yyz_0,  \
                             ta2_xx_0_yyz_1,  \
                             ta2_xx_0_yzz_0,  \
                             ta2_xx_0_yzz_1,  \
                             ta2_xx_0_zzz_0,  \
                             ta2_xx_0_zzz_1,  \
                             ta2_xx_x_xx_0,   \
                             ta2_xx_x_xx_1,   \
                             ta2_xx_x_xxx_0,  \
                             ta2_xx_x_xxx_1,  \
                             ta2_xx_x_xxy_0,  \
                             ta2_xx_x_xxy_1,  \
                             ta2_xx_x_xxz_0,  \
                             ta2_xx_x_xxz_1,  \
                             ta2_xx_x_xy_0,   \
                             ta2_xx_x_xy_1,   \
                             ta2_xx_x_xyy_0,  \
                             ta2_xx_x_xyy_1,  \
                             ta2_xx_x_xyz_0,  \
                             ta2_xx_x_xyz_1,  \
                             ta2_xx_x_xz_0,   \
                             ta2_xx_x_xz_1,   \
                             ta2_xx_x_xzz_0,  \
                             ta2_xx_x_xzz_1,  \
                             ta2_xx_x_yy_0,   \
                             ta2_xx_x_yy_1,   \
                             ta2_xx_x_yyy_0,  \
                             ta2_xx_x_yyy_1,  \
                             ta2_xx_x_yyz_0,  \
                             ta2_xx_x_yyz_1,  \
                             ta2_xx_x_yz_0,   \
                             ta2_xx_x_yz_1,   \
                             ta2_xx_x_yzz_0,  \
                             ta2_xx_x_yzz_1,  \
                             ta2_xx_x_zz_0,   \
                             ta2_xx_x_zz_1,   \
                             ta2_xx_x_zzz_0,  \
                             ta2_xx_x_zzz_1,  \
                             ta2_xx_xx_xxx_0, \
                             ta2_xx_xx_xxy_0, \
                             ta2_xx_xx_xxz_0, \
                             ta2_xx_xx_xyy_0, \
                             ta2_xx_xx_xyz_0, \
                             ta2_xx_xx_xzz_0, \
                             ta2_xx_xx_yyy_0, \
                             ta2_xx_xx_yyz_0, \
                             ta2_xx_xx_yzz_0, \
                             ta2_xx_xx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xx_xxx_0[i] = ta2_xx_0_xxx_0[i] * fe_0 - ta2_xx_0_xxx_1[i] * fe_0 + 3.0 * ta2_xx_x_xx_0[i] * fe_0 - 3.0 * ta2_xx_x_xx_1[i] * fe_0 +
                             2.0 * ta1_x_x_xxx_1[i] + ta2_xx_x_xxx_0[i] * pa_x[i] - ta2_xx_x_xxx_1[i] * pc_x[i];

        ta2_xx_xx_xxy_0[i] = ta2_xx_0_xxy_0[i] * fe_0 - ta2_xx_0_xxy_1[i] * fe_0 + 2.0 * ta2_xx_x_xy_0[i] * fe_0 - 2.0 * ta2_xx_x_xy_1[i] * fe_0 +
                             2.0 * ta1_x_x_xxy_1[i] + ta2_xx_x_xxy_0[i] * pa_x[i] - ta2_xx_x_xxy_1[i] * pc_x[i];

        ta2_xx_xx_xxz_0[i] = ta2_xx_0_xxz_0[i] * fe_0 - ta2_xx_0_xxz_1[i] * fe_0 + 2.0 * ta2_xx_x_xz_0[i] * fe_0 - 2.0 * ta2_xx_x_xz_1[i] * fe_0 +
                             2.0 * ta1_x_x_xxz_1[i] + ta2_xx_x_xxz_0[i] * pa_x[i] - ta2_xx_x_xxz_1[i] * pc_x[i];

        ta2_xx_xx_xyy_0[i] = ta2_xx_0_xyy_0[i] * fe_0 - ta2_xx_0_xyy_1[i] * fe_0 + ta2_xx_x_yy_0[i] * fe_0 - ta2_xx_x_yy_1[i] * fe_0 +
                             2.0 * ta1_x_x_xyy_1[i] + ta2_xx_x_xyy_0[i] * pa_x[i] - ta2_xx_x_xyy_1[i] * pc_x[i];

        ta2_xx_xx_xyz_0[i] = ta2_xx_0_xyz_0[i] * fe_0 - ta2_xx_0_xyz_1[i] * fe_0 + ta2_xx_x_yz_0[i] * fe_0 - ta2_xx_x_yz_1[i] * fe_0 +
                             2.0 * ta1_x_x_xyz_1[i] + ta2_xx_x_xyz_0[i] * pa_x[i] - ta2_xx_x_xyz_1[i] * pc_x[i];

        ta2_xx_xx_xzz_0[i] = ta2_xx_0_xzz_0[i] * fe_0 - ta2_xx_0_xzz_1[i] * fe_0 + ta2_xx_x_zz_0[i] * fe_0 - ta2_xx_x_zz_1[i] * fe_0 +
                             2.0 * ta1_x_x_xzz_1[i] + ta2_xx_x_xzz_0[i] * pa_x[i] - ta2_xx_x_xzz_1[i] * pc_x[i];

        ta2_xx_xx_yyy_0[i] =
            ta2_xx_0_yyy_0[i] * fe_0 - ta2_xx_0_yyy_1[i] * fe_0 + 2.0 * ta1_x_x_yyy_1[i] + ta2_xx_x_yyy_0[i] * pa_x[i] - ta2_xx_x_yyy_1[i] * pc_x[i];

        ta2_xx_xx_yyz_0[i] =
            ta2_xx_0_yyz_0[i] * fe_0 - ta2_xx_0_yyz_1[i] * fe_0 + 2.0 * ta1_x_x_yyz_1[i] + ta2_xx_x_yyz_0[i] * pa_x[i] - ta2_xx_x_yyz_1[i] * pc_x[i];

        ta2_xx_xx_yzz_0[i] =
            ta2_xx_0_yzz_0[i] * fe_0 - ta2_xx_0_yzz_1[i] * fe_0 + 2.0 * ta1_x_x_yzz_1[i] + ta2_xx_x_yzz_0[i] * pa_x[i] - ta2_xx_x_yzz_1[i] * pc_x[i];

        ta2_xx_xx_zzz_0[i] =
            ta2_xx_0_zzz_0[i] * fe_0 - ta2_xx_0_zzz_1[i] * fe_0 + 2.0 * ta1_x_x_zzz_1[i] + ta2_xx_x_zzz_0[i] * pa_x[i] - ta2_xx_x_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : DF

    auto ta2_xx_xy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 10);

    auto ta2_xx_xy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 11);

    auto ta2_xx_xy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 12);

    auto ta2_xx_xy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 13);

    auto ta2_xx_xy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 14);

    auto ta2_xx_xy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 15);

    auto ta2_xx_xy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 16);

    auto ta2_xx_xy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 17);

    auto ta2_xx_xy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 18);

    auto ta2_xx_xy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 19);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_y_yyy_1,   \
                             ta1_x_y_yyz_1,   \
                             ta1_x_y_yzz_1,   \
                             ta2_xx_x_xx_0,   \
                             ta2_xx_x_xx_1,   \
                             ta2_xx_x_xxx_0,  \
                             ta2_xx_x_xxx_1,  \
                             ta2_xx_x_xxy_0,  \
                             ta2_xx_x_xxy_1,  \
                             ta2_xx_x_xxz_0,  \
                             ta2_xx_x_xxz_1,  \
                             ta2_xx_x_xy_0,   \
                             ta2_xx_x_xy_1,   \
                             ta2_xx_x_xyy_0,  \
                             ta2_xx_x_xyy_1,  \
                             ta2_xx_x_xyz_0,  \
                             ta2_xx_x_xyz_1,  \
                             ta2_xx_x_xz_0,   \
                             ta2_xx_x_xz_1,   \
                             ta2_xx_x_xzz_0,  \
                             ta2_xx_x_xzz_1,  \
                             ta2_xx_x_zzz_0,  \
                             ta2_xx_x_zzz_1,  \
                             ta2_xx_xy_xxx_0, \
                             ta2_xx_xy_xxy_0, \
                             ta2_xx_xy_xxz_0, \
                             ta2_xx_xy_xyy_0, \
                             ta2_xx_xy_xyz_0, \
                             ta2_xx_xy_xzz_0, \
                             ta2_xx_xy_yyy_0, \
                             ta2_xx_xy_yyz_0, \
                             ta2_xx_xy_yzz_0, \
                             ta2_xx_xy_zzz_0, \
                             ta2_xx_y_yyy_0,  \
                             ta2_xx_y_yyy_1,  \
                             ta2_xx_y_yyz_0,  \
                             ta2_xx_y_yyz_1,  \
                             ta2_xx_y_yzz_0,  \
                             ta2_xx_y_yzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xy_xxx_0[i] = ta2_xx_x_xxx_0[i] * pa_y[i] - ta2_xx_x_xxx_1[i] * pc_y[i];

        ta2_xx_xy_xxy_0[i] = ta2_xx_x_xx_0[i] * fe_0 - ta2_xx_x_xx_1[i] * fe_0 + ta2_xx_x_xxy_0[i] * pa_y[i] - ta2_xx_x_xxy_1[i] * pc_y[i];

        ta2_xx_xy_xxz_0[i] = ta2_xx_x_xxz_0[i] * pa_y[i] - ta2_xx_x_xxz_1[i] * pc_y[i];

        ta2_xx_xy_xyy_0[i] =
            2.0 * ta2_xx_x_xy_0[i] * fe_0 - 2.0 * ta2_xx_x_xy_1[i] * fe_0 + ta2_xx_x_xyy_0[i] * pa_y[i] - ta2_xx_x_xyy_1[i] * pc_y[i];

        ta2_xx_xy_xyz_0[i] = ta2_xx_x_xz_0[i] * fe_0 - ta2_xx_x_xz_1[i] * fe_0 + ta2_xx_x_xyz_0[i] * pa_y[i] - ta2_xx_x_xyz_1[i] * pc_y[i];

        ta2_xx_xy_xzz_0[i] = ta2_xx_x_xzz_0[i] * pa_y[i] - ta2_xx_x_xzz_1[i] * pc_y[i];

        ta2_xx_xy_yyy_0[i] = 2.0 * ta1_x_y_yyy_1[i] + ta2_xx_y_yyy_0[i] * pa_x[i] - ta2_xx_y_yyy_1[i] * pc_x[i];

        ta2_xx_xy_yyz_0[i] = 2.0 * ta1_x_y_yyz_1[i] + ta2_xx_y_yyz_0[i] * pa_x[i] - ta2_xx_y_yyz_1[i] * pc_x[i];

        ta2_xx_xy_yzz_0[i] = 2.0 * ta1_x_y_yzz_1[i] + ta2_xx_y_yzz_0[i] * pa_x[i] - ta2_xx_y_yzz_1[i] * pc_x[i];

        ta2_xx_xy_zzz_0[i] = ta2_xx_x_zzz_0[i] * pa_y[i] - ta2_xx_x_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : DF

    auto ta2_xx_xz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 20);

    auto ta2_xx_xz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 21);

    auto ta2_xx_xz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 22);

    auto ta2_xx_xz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 23);

    auto ta2_xx_xz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 24);

    auto ta2_xx_xz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 25);

    auto ta2_xx_xz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 26);

    auto ta2_xx_xz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 27);

    auto ta2_xx_xz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 28);

    auto ta2_xx_xz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 29);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_x_z_yyz_1,   \
                             ta1_x_z_yzz_1,   \
                             ta1_x_z_zzz_1,   \
                             ta2_xx_x_xx_0,   \
                             ta2_xx_x_xx_1,   \
                             ta2_xx_x_xxx_0,  \
                             ta2_xx_x_xxx_1,  \
                             ta2_xx_x_xxy_0,  \
                             ta2_xx_x_xxy_1,  \
                             ta2_xx_x_xxz_0,  \
                             ta2_xx_x_xxz_1,  \
                             ta2_xx_x_xy_0,   \
                             ta2_xx_x_xy_1,   \
                             ta2_xx_x_xyy_0,  \
                             ta2_xx_x_xyy_1,  \
                             ta2_xx_x_xyz_0,  \
                             ta2_xx_x_xyz_1,  \
                             ta2_xx_x_xz_0,   \
                             ta2_xx_x_xz_1,   \
                             ta2_xx_x_xzz_0,  \
                             ta2_xx_x_xzz_1,  \
                             ta2_xx_x_yyy_0,  \
                             ta2_xx_x_yyy_1,  \
                             ta2_xx_xz_xxx_0, \
                             ta2_xx_xz_xxy_0, \
                             ta2_xx_xz_xxz_0, \
                             ta2_xx_xz_xyy_0, \
                             ta2_xx_xz_xyz_0, \
                             ta2_xx_xz_xzz_0, \
                             ta2_xx_xz_yyy_0, \
                             ta2_xx_xz_yyz_0, \
                             ta2_xx_xz_yzz_0, \
                             ta2_xx_xz_zzz_0, \
                             ta2_xx_z_yyz_0,  \
                             ta2_xx_z_yyz_1,  \
                             ta2_xx_z_yzz_0,  \
                             ta2_xx_z_yzz_1,  \
                             ta2_xx_z_zzz_0,  \
                             ta2_xx_z_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xz_xxx_0[i] = ta2_xx_x_xxx_0[i] * pa_z[i] - ta2_xx_x_xxx_1[i] * pc_z[i];

        ta2_xx_xz_xxy_0[i] = ta2_xx_x_xxy_0[i] * pa_z[i] - ta2_xx_x_xxy_1[i] * pc_z[i];

        ta2_xx_xz_xxz_0[i] = ta2_xx_x_xx_0[i] * fe_0 - ta2_xx_x_xx_1[i] * fe_0 + ta2_xx_x_xxz_0[i] * pa_z[i] - ta2_xx_x_xxz_1[i] * pc_z[i];

        ta2_xx_xz_xyy_0[i] = ta2_xx_x_xyy_0[i] * pa_z[i] - ta2_xx_x_xyy_1[i] * pc_z[i];

        ta2_xx_xz_xyz_0[i] = ta2_xx_x_xy_0[i] * fe_0 - ta2_xx_x_xy_1[i] * fe_0 + ta2_xx_x_xyz_0[i] * pa_z[i] - ta2_xx_x_xyz_1[i] * pc_z[i];

        ta2_xx_xz_xzz_0[i] =
            2.0 * ta2_xx_x_xz_0[i] * fe_0 - 2.0 * ta2_xx_x_xz_1[i] * fe_0 + ta2_xx_x_xzz_0[i] * pa_z[i] - ta2_xx_x_xzz_1[i] * pc_z[i];

        ta2_xx_xz_yyy_0[i] = ta2_xx_x_yyy_0[i] * pa_z[i] - ta2_xx_x_yyy_1[i] * pc_z[i];

        ta2_xx_xz_yyz_0[i] = 2.0 * ta1_x_z_yyz_1[i] + ta2_xx_z_yyz_0[i] * pa_x[i] - ta2_xx_z_yyz_1[i] * pc_x[i];

        ta2_xx_xz_yzz_0[i] = 2.0 * ta1_x_z_yzz_1[i] + ta2_xx_z_yzz_0[i] * pa_x[i] - ta2_xx_z_yzz_1[i] * pc_x[i];

        ta2_xx_xz_zzz_0[i] = 2.0 * ta1_x_z_zzz_1[i] + ta2_xx_z_zzz_0[i] * pa_x[i] - ta2_xx_z_zzz_1[i] * pc_x[i];
    }

    // Set up 30-40 components of targeted buffer : DF

    auto ta2_xx_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 30);

    auto ta2_xx_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 31);

    auto ta2_xx_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 32);

    auto ta2_xx_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 33);

    auto ta2_xx_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 34);

    auto ta2_xx_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 35);

    auto ta2_xx_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 36);

    auto ta2_xx_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 37);

    auto ta2_xx_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 38);

    auto ta2_xx_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 39);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta2_xx_0_xxx_0,  \
                             ta2_xx_0_xxx_1,  \
                             ta2_xx_0_xxy_0,  \
                             ta2_xx_0_xxy_1,  \
                             ta2_xx_0_xxz_0,  \
                             ta2_xx_0_xxz_1,  \
                             ta2_xx_0_xyy_0,  \
                             ta2_xx_0_xyy_1,  \
                             ta2_xx_0_xyz_0,  \
                             ta2_xx_0_xyz_1,  \
                             ta2_xx_0_xzz_0,  \
                             ta2_xx_0_xzz_1,  \
                             ta2_xx_0_yyy_0,  \
                             ta2_xx_0_yyy_1,  \
                             ta2_xx_0_yyz_0,  \
                             ta2_xx_0_yyz_1,  \
                             ta2_xx_0_yzz_0,  \
                             ta2_xx_0_yzz_1,  \
                             ta2_xx_0_zzz_0,  \
                             ta2_xx_0_zzz_1,  \
                             ta2_xx_y_xx_0,   \
                             ta2_xx_y_xx_1,   \
                             ta2_xx_y_xxx_0,  \
                             ta2_xx_y_xxx_1,  \
                             ta2_xx_y_xxy_0,  \
                             ta2_xx_y_xxy_1,  \
                             ta2_xx_y_xxz_0,  \
                             ta2_xx_y_xxz_1,  \
                             ta2_xx_y_xy_0,   \
                             ta2_xx_y_xy_1,   \
                             ta2_xx_y_xyy_0,  \
                             ta2_xx_y_xyy_1,  \
                             ta2_xx_y_xyz_0,  \
                             ta2_xx_y_xyz_1,  \
                             ta2_xx_y_xz_0,   \
                             ta2_xx_y_xz_1,   \
                             ta2_xx_y_xzz_0,  \
                             ta2_xx_y_xzz_1,  \
                             ta2_xx_y_yy_0,   \
                             ta2_xx_y_yy_1,   \
                             ta2_xx_y_yyy_0,  \
                             ta2_xx_y_yyy_1,  \
                             ta2_xx_y_yyz_0,  \
                             ta2_xx_y_yyz_1,  \
                             ta2_xx_y_yz_0,   \
                             ta2_xx_y_yz_1,   \
                             ta2_xx_y_yzz_0,  \
                             ta2_xx_y_yzz_1,  \
                             ta2_xx_y_zz_0,   \
                             ta2_xx_y_zz_1,   \
                             ta2_xx_y_zzz_0,  \
                             ta2_xx_y_zzz_1,  \
                             ta2_xx_yy_xxx_0, \
                             ta2_xx_yy_xxy_0, \
                             ta2_xx_yy_xxz_0, \
                             ta2_xx_yy_xyy_0, \
                             ta2_xx_yy_xyz_0, \
                             ta2_xx_yy_xzz_0, \
                             ta2_xx_yy_yyy_0, \
                             ta2_xx_yy_yyz_0, \
                             ta2_xx_yy_yzz_0, \
                             ta2_xx_yy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yy_xxx_0[i] = ta2_xx_0_xxx_0[i] * fe_0 - ta2_xx_0_xxx_1[i] * fe_0 + ta2_xx_y_xxx_0[i] * pa_y[i] - ta2_xx_y_xxx_1[i] * pc_y[i];

        ta2_xx_yy_xxy_0[i] = ta2_xx_0_xxy_0[i] * fe_0 - ta2_xx_0_xxy_1[i] * fe_0 + ta2_xx_y_xx_0[i] * fe_0 - ta2_xx_y_xx_1[i] * fe_0 +
                             ta2_xx_y_xxy_0[i] * pa_y[i] - ta2_xx_y_xxy_1[i] * pc_y[i];

        ta2_xx_yy_xxz_0[i] = ta2_xx_0_xxz_0[i] * fe_0 - ta2_xx_0_xxz_1[i] * fe_0 + ta2_xx_y_xxz_0[i] * pa_y[i] - ta2_xx_y_xxz_1[i] * pc_y[i];

        ta2_xx_yy_xyy_0[i] = ta2_xx_0_xyy_0[i] * fe_0 - ta2_xx_0_xyy_1[i] * fe_0 + 2.0 * ta2_xx_y_xy_0[i] * fe_0 - 2.0 * ta2_xx_y_xy_1[i] * fe_0 +
                             ta2_xx_y_xyy_0[i] * pa_y[i] - ta2_xx_y_xyy_1[i] * pc_y[i];

        ta2_xx_yy_xyz_0[i] = ta2_xx_0_xyz_0[i] * fe_0 - ta2_xx_0_xyz_1[i] * fe_0 + ta2_xx_y_xz_0[i] * fe_0 - ta2_xx_y_xz_1[i] * fe_0 +
                             ta2_xx_y_xyz_0[i] * pa_y[i] - ta2_xx_y_xyz_1[i] * pc_y[i];

        ta2_xx_yy_xzz_0[i] = ta2_xx_0_xzz_0[i] * fe_0 - ta2_xx_0_xzz_1[i] * fe_0 + ta2_xx_y_xzz_0[i] * pa_y[i] - ta2_xx_y_xzz_1[i] * pc_y[i];

        ta2_xx_yy_yyy_0[i] = ta2_xx_0_yyy_0[i] * fe_0 - ta2_xx_0_yyy_1[i] * fe_0 + 3.0 * ta2_xx_y_yy_0[i] * fe_0 - 3.0 * ta2_xx_y_yy_1[i] * fe_0 +
                             ta2_xx_y_yyy_0[i] * pa_y[i] - ta2_xx_y_yyy_1[i] * pc_y[i];

        ta2_xx_yy_yyz_0[i] = ta2_xx_0_yyz_0[i] * fe_0 - ta2_xx_0_yyz_1[i] * fe_0 + 2.0 * ta2_xx_y_yz_0[i] * fe_0 - 2.0 * ta2_xx_y_yz_1[i] * fe_0 +
                             ta2_xx_y_yyz_0[i] * pa_y[i] - ta2_xx_y_yyz_1[i] * pc_y[i];

        ta2_xx_yy_yzz_0[i] = ta2_xx_0_yzz_0[i] * fe_0 - ta2_xx_0_yzz_1[i] * fe_0 + ta2_xx_y_zz_0[i] * fe_0 - ta2_xx_y_zz_1[i] * fe_0 +
                             ta2_xx_y_yzz_0[i] * pa_y[i] - ta2_xx_y_yzz_1[i] * pc_y[i];

        ta2_xx_yy_zzz_0[i] = ta2_xx_0_zzz_0[i] * fe_0 - ta2_xx_0_zzz_1[i] * fe_0 + ta2_xx_y_zzz_0[i] * pa_y[i] - ta2_xx_y_zzz_1[i] * pc_y[i];
    }

    // Set up 40-50 components of targeted buffer : DF

    auto ta2_xx_yz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 40);

    auto ta2_xx_yz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 41);

    auto ta2_xx_yz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 42);

    auto ta2_xx_yz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 43);

    auto ta2_xx_yz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 44);

    auto ta2_xx_yz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 45);

    auto ta2_xx_yz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 46);

    auto ta2_xx_yz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 47);

    auto ta2_xx_yz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 48);

    auto ta2_xx_yz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 49);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta2_xx_y_xxy_0,  \
                             ta2_xx_y_xxy_1,  \
                             ta2_xx_y_xyy_0,  \
                             ta2_xx_y_xyy_1,  \
                             ta2_xx_y_yyy_0,  \
                             ta2_xx_y_yyy_1,  \
                             ta2_xx_yz_xxx_0, \
                             ta2_xx_yz_xxy_0, \
                             ta2_xx_yz_xxz_0, \
                             ta2_xx_yz_xyy_0, \
                             ta2_xx_yz_xyz_0, \
                             ta2_xx_yz_xzz_0, \
                             ta2_xx_yz_yyy_0, \
                             ta2_xx_yz_yyz_0, \
                             ta2_xx_yz_yzz_0, \
                             ta2_xx_yz_zzz_0, \
                             ta2_xx_z_xxx_0,  \
                             ta2_xx_z_xxx_1,  \
                             ta2_xx_z_xxz_0,  \
                             ta2_xx_z_xxz_1,  \
                             ta2_xx_z_xyz_0,  \
                             ta2_xx_z_xyz_1,  \
                             ta2_xx_z_xz_0,   \
                             ta2_xx_z_xz_1,   \
                             ta2_xx_z_xzz_0,  \
                             ta2_xx_z_xzz_1,  \
                             ta2_xx_z_yyz_0,  \
                             ta2_xx_z_yyz_1,  \
                             ta2_xx_z_yz_0,   \
                             ta2_xx_z_yz_1,   \
                             ta2_xx_z_yzz_0,  \
                             ta2_xx_z_yzz_1,  \
                             ta2_xx_z_zz_0,   \
                             ta2_xx_z_zz_1,   \
                             ta2_xx_z_zzz_0,  \
                             ta2_xx_z_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yz_xxx_0[i] = ta2_xx_z_xxx_0[i] * pa_y[i] - ta2_xx_z_xxx_1[i] * pc_y[i];

        ta2_xx_yz_xxy_0[i] = ta2_xx_y_xxy_0[i] * pa_z[i] - ta2_xx_y_xxy_1[i] * pc_z[i];

        ta2_xx_yz_xxz_0[i] = ta2_xx_z_xxz_0[i] * pa_y[i] - ta2_xx_z_xxz_1[i] * pc_y[i];

        ta2_xx_yz_xyy_0[i] = ta2_xx_y_xyy_0[i] * pa_z[i] - ta2_xx_y_xyy_1[i] * pc_z[i];

        ta2_xx_yz_xyz_0[i] = ta2_xx_z_xz_0[i] * fe_0 - ta2_xx_z_xz_1[i] * fe_0 + ta2_xx_z_xyz_0[i] * pa_y[i] - ta2_xx_z_xyz_1[i] * pc_y[i];

        ta2_xx_yz_xzz_0[i] = ta2_xx_z_xzz_0[i] * pa_y[i] - ta2_xx_z_xzz_1[i] * pc_y[i];

        ta2_xx_yz_yyy_0[i] = ta2_xx_y_yyy_0[i] * pa_z[i] - ta2_xx_y_yyy_1[i] * pc_z[i];

        ta2_xx_yz_yyz_0[i] =
            2.0 * ta2_xx_z_yz_0[i] * fe_0 - 2.0 * ta2_xx_z_yz_1[i] * fe_0 + ta2_xx_z_yyz_0[i] * pa_y[i] - ta2_xx_z_yyz_1[i] * pc_y[i];

        ta2_xx_yz_yzz_0[i] = ta2_xx_z_zz_0[i] * fe_0 - ta2_xx_z_zz_1[i] * fe_0 + ta2_xx_z_yzz_0[i] * pa_y[i] - ta2_xx_z_yzz_1[i] * pc_y[i];

        ta2_xx_yz_zzz_0[i] = ta2_xx_z_zzz_0[i] * pa_y[i] - ta2_xx_z_zzz_1[i] * pc_y[i];
    }

    // Set up 50-60 components of targeted buffer : DF

    auto ta2_xx_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 50);

    auto ta2_xx_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 51);

    auto ta2_xx_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 52);

    auto ta2_xx_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 53);

    auto ta2_xx_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 54);

    auto ta2_xx_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 55);

    auto ta2_xx_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 56);

    auto ta2_xx_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 57);

    auto ta2_xx_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 58);

    auto ta2_xx_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 59);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta2_xx_0_xxx_0,  \
                             ta2_xx_0_xxx_1,  \
                             ta2_xx_0_xxy_0,  \
                             ta2_xx_0_xxy_1,  \
                             ta2_xx_0_xxz_0,  \
                             ta2_xx_0_xxz_1,  \
                             ta2_xx_0_xyy_0,  \
                             ta2_xx_0_xyy_1,  \
                             ta2_xx_0_xyz_0,  \
                             ta2_xx_0_xyz_1,  \
                             ta2_xx_0_xzz_0,  \
                             ta2_xx_0_xzz_1,  \
                             ta2_xx_0_yyy_0,  \
                             ta2_xx_0_yyy_1,  \
                             ta2_xx_0_yyz_0,  \
                             ta2_xx_0_yyz_1,  \
                             ta2_xx_0_yzz_0,  \
                             ta2_xx_0_yzz_1,  \
                             ta2_xx_0_zzz_0,  \
                             ta2_xx_0_zzz_1,  \
                             ta2_xx_z_xx_0,   \
                             ta2_xx_z_xx_1,   \
                             ta2_xx_z_xxx_0,  \
                             ta2_xx_z_xxx_1,  \
                             ta2_xx_z_xxy_0,  \
                             ta2_xx_z_xxy_1,  \
                             ta2_xx_z_xxz_0,  \
                             ta2_xx_z_xxz_1,  \
                             ta2_xx_z_xy_0,   \
                             ta2_xx_z_xy_1,   \
                             ta2_xx_z_xyy_0,  \
                             ta2_xx_z_xyy_1,  \
                             ta2_xx_z_xyz_0,  \
                             ta2_xx_z_xyz_1,  \
                             ta2_xx_z_xz_0,   \
                             ta2_xx_z_xz_1,   \
                             ta2_xx_z_xzz_0,  \
                             ta2_xx_z_xzz_1,  \
                             ta2_xx_z_yy_0,   \
                             ta2_xx_z_yy_1,   \
                             ta2_xx_z_yyy_0,  \
                             ta2_xx_z_yyy_1,  \
                             ta2_xx_z_yyz_0,  \
                             ta2_xx_z_yyz_1,  \
                             ta2_xx_z_yz_0,   \
                             ta2_xx_z_yz_1,   \
                             ta2_xx_z_yzz_0,  \
                             ta2_xx_z_yzz_1,  \
                             ta2_xx_z_zz_0,   \
                             ta2_xx_z_zz_1,   \
                             ta2_xx_z_zzz_0,  \
                             ta2_xx_z_zzz_1,  \
                             ta2_xx_zz_xxx_0, \
                             ta2_xx_zz_xxy_0, \
                             ta2_xx_zz_xxz_0, \
                             ta2_xx_zz_xyy_0, \
                             ta2_xx_zz_xyz_0, \
                             ta2_xx_zz_xzz_0, \
                             ta2_xx_zz_yyy_0, \
                             ta2_xx_zz_yyz_0, \
                             ta2_xx_zz_yzz_0, \
                             ta2_xx_zz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_zz_xxx_0[i] = ta2_xx_0_xxx_0[i] * fe_0 - ta2_xx_0_xxx_1[i] * fe_0 + ta2_xx_z_xxx_0[i] * pa_z[i] - ta2_xx_z_xxx_1[i] * pc_z[i];

        ta2_xx_zz_xxy_0[i] = ta2_xx_0_xxy_0[i] * fe_0 - ta2_xx_0_xxy_1[i] * fe_0 + ta2_xx_z_xxy_0[i] * pa_z[i] - ta2_xx_z_xxy_1[i] * pc_z[i];

        ta2_xx_zz_xxz_0[i] = ta2_xx_0_xxz_0[i] * fe_0 - ta2_xx_0_xxz_1[i] * fe_0 + ta2_xx_z_xx_0[i] * fe_0 - ta2_xx_z_xx_1[i] * fe_0 +
                             ta2_xx_z_xxz_0[i] * pa_z[i] - ta2_xx_z_xxz_1[i] * pc_z[i];

        ta2_xx_zz_xyy_0[i] = ta2_xx_0_xyy_0[i] * fe_0 - ta2_xx_0_xyy_1[i] * fe_0 + ta2_xx_z_xyy_0[i] * pa_z[i] - ta2_xx_z_xyy_1[i] * pc_z[i];

        ta2_xx_zz_xyz_0[i] = ta2_xx_0_xyz_0[i] * fe_0 - ta2_xx_0_xyz_1[i] * fe_0 + ta2_xx_z_xy_0[i] * fe_0 - ta2_xx_z_xy_1[i] * fe_0 +
                             ta2_xx_z_xyz_0[i] * pa_z[i] - ta2_xx_z_xyz_1[i] * pc_z[i];

        ta2_xx_zz_xzz_0[i] = ta2_xx_0_xzz_0[i] * fe_0 - ta2_xx_0_xzz_1[i] * fe_0 + 2.0 * ta2_xx_z_xz_0[i] * fe_0 - 2.0 * ta2_xx_z_xz_1[i] * fe_0 +
                             ta2_xx_z_xzz_0[i] * pa_z[i] - ta2_xx_z_xzz_1[i] * pc_z[i];

        ta2_xx_zz_yyy_0[i] = ta2_xx_0_yyy_0[i] * fe_0 - ta2_xx_0_yyy_1[i] * fe_0 + ta2_xx_z_yyy_0[i] * pa_z[i] - ta2_xx_z_yyy_1[i] * pc_z[i];

        ta2_xx_zz_yyz_0[i] = ta2_xx_0_yyz_0[i] * fe_0 - ta2_xx_0_yyz_1[i] * fe_0 + ta2_xx_z_yy_0[i] * fe_0 - ta2_xx_z_yy_1[i] * fe_0 +
                             ta2_xx_z_yyz_0[i] * pa_z[i] - ta2_xx_z_yyz_1[i] * pc_z[i];

        ta2_xx_zz_yzz_0[i] = ta2_xx_0_yzz_0[i] * fe_0 - ta2_xx_0_yzz_1[i] * fe_0 + 2.0 * ta2_xx_z_yz_0[i] * fe_0 - 2.0 * ta2_xx_z_yz_1[i] * fe_0 +
                             ta2_xx_z_yzz_0[i] * pa_z[i] - ta2_xx_z_yzz_1[i] * pc_z[i];

        ta2_xx_zz_zzz_0[i] = ta2_xx_0_zzz_0[i] * fe_0 - ta2_xx_0_zzz_1[i] * fe_0 + 3.0 * ta2_xx_z_zz_0[i] * fe_0 - 3.0 * ta2_xx_z_zz_1[i] * fe_0 +
                             ta2_xx_z_zzz_0[i] * pa_z[i] - ta2_xx_z_zzz_1[i] * pc_z[i];
    }

    // Set up 60-70 components of targeted buffer : DF

    auto ta2_xy_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 60);

    auto ta2_xy_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 61);

    auto ta2_xy_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 62);

    auto ta2_xy_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 63);

    auto ta2_xy_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 64);

    auto ta2_xy_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 65);

    auto ta2_xy_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 66);

    auto ta2_xy_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 67);

    auto ta2_xy_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 68);

    auto ta2_xy_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 69);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_y_x_xxx_1,   \
                             ta1_y_x_xxy_1,   \
                             ta1_y_x_xxz_1,   \
                             ta1_y_x_xyy_1,   \
                             ta1_y_x_xyz_1,   \
                             ta1_y_x_xzz_1,   \
                             ta1_y_x_yyy_1,   \
                             ta1_y_x_yyz_1,   \
                             ta1_y_x_yzz_1,   \
                             ta1_y_x_zzz_1,   \
                             ta2_xy_0_xxx_0,  \
                             ta2_xy_0_xxx_1,  \
                             ta2_xy_0_xxy_0,  \
                             ta2_xy_0_xxy_1,  \
                             ta2_xy_0_xxz_0,  \
                             ta2_xy_0_xxz_1,  \
                             ta2_xy_0_xyy_0,  \
                             ta2_xy_0_xyy_1,  \
                             ta2_xy_0_xyz_0,  \
                             ta2_xy_0_xyz_1,  \
                             ta2_xy_0_xzz_0,  \
                             ta2_xy_0_xzz_1,  \
                             ta2_xy_0_yyy_0,  \
                             ta2_xy_0_yyy_1,  \
                             ta2_xy_0_yyz_0,  \
                             ta2_xy_0_yyz_1,  \
                             ta2_xy_0_yzz_0,  \
                             ta2_xy_0_yzz_1,  \
                             ta2_xy_0_zzz_0,  \
                             ta2_xy_0_zzz_1,  \
                             ta2_xy_x_xx_0,   \
                             ta2_xy_x_xx_1,   \
                             ta2_xy_x_xxx_0,  \
                             ta2_xy_x_xxx_1,  \
                             ta2_xy_x_xxy_0,  \
                             ta2_xy_x_xxy_1,  \
                             ta2_xy_x_xxz_0,  \
                             ta2_xy_x_xxz_1,  \
                             ta2_xy_x_xy_0,   \
                             ta2_xy_x_xy_1,   \
                             ta2_xy_x_xyy_0,  \
                             ta2_xy_x_xyy_1,  \
                             ta2_xy_x_xyz_0,  \
                             ta2_xy_x_xyz_1,  \
                             ta2_xy_x_xz_0,   \
                             ta2_xy_x_xz_1,   \
                             ta2_xy_x_xzz_0,  \
                             ta2_xy_x_xzz_1,  \
                             ta2_xy_x_yy_0,   \
                             ta2_xy_x_yy_1,   \
                             ta2_xy_x_yyy_0,  \
                             ta2_xy_x_yyy_1,  \
                             ta2_xy_x_yyz_0,  \
                             ta2_xy_x_yyz_1,  \
                             ta2_xy_x_yz_0,   \
                             ta2_xy_x_yz_1,   \
                             ta2_xy_x_yzz_0,  \
                             ta2_xy_x_yzz_1,  \
                             ta2_xy_x_zz_0,   \
                             ta2_xy_x_zz_1,   \
                             ta2_xy_x_zzz_0,  \
                             ta2_xy_x_zzz_1,  \
                             ta2_xy_xx_xxx_0, \
                             ta2_xy_xx_xxy_0, \
                             ta2_xy_xx_xxz_0, \
                             ta2_xy_xx_xyy_0, \
                             ta2_xy_xx_xyz_0, \
                             ta2_xy_xx_xzz_0, \
                             ta2_xy_xx_yyy_0, \
                             ta2_xy_xx_yyz_0, \
                             ta2_xy_xx_yzz_0, \
                             ta2_xy_xx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xx_xxx_0[i] = ta2_xy_0_xxx_0[i] * fe_0 - ta2_xy_0_xxx_1[i] * fe_0 + 3.0 * ta2_xy_x_xx_0[i] * fe_0 - 3.0 * ta2_xy_x_xx_1[i] * fe_0 +
                             ta1_y_x_xxx_1[i] + ta2_xy_x_xxx_0[i] * pa_x[i] - ta2_xy_x_xxx_1[i] * pc_x[i];

        ta2_xy_xx_xxy_0[i] = ta2_xy_0_xxy_0[i] * fe_0 - ta2_xy_0_xxy_1[i] * fe_0 + 2.0 * ta2_xy_x_xy_0[i] * fe_0 - 2.0 * ta2_xy_x_xy_1[i] * fe_0 +
                             ta1_y_x_xxy_1[i] + ta2_xy_x_xxy_0[i] * pa_x[i] - ta2_xy_x_xxy_1[i] * pc_x[i];

        ta2_xy_xx_xxz_0[i] = ta2_xy_0_xxz_0[i] * fe_0 - ta2_xy_0_xxz_1[i] * fe_0 + 2.0 * ta2_xy_x_xz_0[i] * fe_0 - 2.0 * ta2_xy_x_xz_1[i] * fe_0 +
                             ta1_y_x_xxz_1[i] + ta2_xy_x_xxz_0[i] * pa_x[i] - ta2_xy_x_xxz_1[i] * pc_x[i];

        ta2_xy_xx_xyy_0[i] = ta2_xy_0_xyy_0[i] * fe_0 - ta2_xy_0_xyy_1[i] * fe_0 + ta2_xy_x_yy_0[i] * fe_0 - ta2_xy_x_yy_1[i] * fe_0 +
                             ta1_y_x_xyy_1[i] + ta2_xy_x_xyy_0[i] * pa_x[i] - ta2_xy_x_xyy_1[i] * pc_x[i];

        ta2_xy_xx_xyz_0[i] = ta2_xy_0_xyz_0[i] * fe_0 - ta2_xy_0_xyz_1[i] * fe_0 + ta2_xy_x_yz_0[i] * fe_0 - ta2_xy_x_yz_1[i] * fe_0 +
                             ta1_y_x_xyz_1[i] + ta2_xy_x_xyz_0[i] * pa_x[i] - ta2_xy_x_xyz_1[i] * pc_x[i];

        ta2_xy_xx_xzz_0[i] = ta2_xy_0_xzz_0[i] * fe_0 - ta2_xy_0_xzz_1[i] * fe_0 + ta2_xy_x_zz_0[i] * fe_0 - ta2_xy_x_zz_1[i] * fe_0 +
                             ta1_y_x_xzz_1[i] + ta2_xy_x_xzz_0[i] * pa_x[i] - ta2_xy_x_xzz_1[i] * pc_x[i];

        ta2_xy_xx_yyy_0[i] =
            ta2_xy_0_yyy_0[i] * fe_0 - ta2_xy_0_yyy_1[i] * fe_0 + ta1_y_x_yyy_1[i] + ta2_xy_x_yyy_0[i] * pa_x[i] - ta2_xy_x_yyy_1[i] * pc_x[i];

        ta2_xy_xx_yyz_0[i] =
            ta2_xy_0_yyz_0[i] * fe_0 - ta2_xy_0_yyz_1[i] * fe_0 + ta1_y_x_yyz_1[i] + ta2_xy_x_yyz_0[i] * pa_x[i] - ta2_xy_x_yyz_1[i] * pc_x[i];

        ta2_xy_xx_yzz_0[i] =
            ta2_xy_0_yzz_0[i] * fe_0 - ta2_xy_0_yzz_1[i] * fe_0 + ta1_y_x_yzz_1[i] + ta2_xy_x_yzz_0[i] * pa_x[i] - ta2_xy_x_yzz_1[i] * pc_x[i];

        ta2_xy_xx_zzz_0[i] =
            ta2_xy_0_zzz_0[i] * fe_0 - ta2_xy_0_zzz_1[i] * fe_0 + ta1_y_x_zzz_1[i] + ta2_xy_x_zzz_0[i] * pa_x[i] - ta2_xy_x_zzz_1[i] * pc_x[i];
    }

    // Set up 70-80 components of targeted buffer : DF

    auto ta2_xy_xy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 70);

    auto ta2_xy_xy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 71);

    auto ta2_xy_xy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 72);

    auto ta2_xy_xy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 73);

    auto ta2_xy_xy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 74);

    auto ta2_xy_xy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 75);

    auto ta2_xy_xy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 76);

    auto ta2_xy_xy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 77);

    auto ta2_xy_xy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 78);

    auto ta2_xy_xy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 79);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_x_x_xxx_1,   \
                             ta1_x_x_xxz_1,   \
                             ta1_x_x_xzz_1,   \
                             ta1_y_y_xxy_1,   \
                             ta1_y_y_xyy_1,   \
                             ta1_y_y_xyz_1,   \
                             ta1_y_y_yyy_1,   \
                             ta1_y_y_yyz_1,   \
                             ta1_y_y_yzz_1,   \
                             ta1_y_y_zzz_1,   \
                             ta2_xy_x_xxx_0,  \
                             ta2_xy_x_xxx_1,  \
                             ta2_xy_x_xxz_0,  \
                             ta2_xy_x_xxz_1,  \
                             ta2_xy_x_xzz_0,  \
                             ta2_xy_x_xzz_1,  \
                             ta2_xy_xy_xxx_0, \
                             ta2_xy_xy_xxy_0, \
                             ta2_xy_xy_xxz_0, \
                             ta2_xy_xy_xyy_0, \
                             ta2_xy_xy_xyz_0, \
                             ta2_xy_xy_xzz_0, \
                             ta2_xy_xy_yyy_0, \
                             ta2_xy_xy_yyz_0, \
                             ta2_xy_xy_yzz_0, \
                             ta2_xy_xy_zzz_0, \
                             ta2_xy_y_xxy_0,  \
                             ta2_xy_y_xxy_1,  \
                             ta2_xy_y_xy_0,   \
                             ta2_xy_y_xy_1,   \
                             ta2_xy_y_xyy_0,  \
                             ta2_xy_y_xyy_1,  \
                             ta2_xy_y_xyz_0,  \
                             ta2_xy_y_xyz_1,  \
                             ta2_xy_y_yy_0,   \
                             ta2_xy_y_yy_1,   \
                             ta2_xy_y_yyy_0,  \
                             ta2_xy_y_yyy_1,  \
                             ta2_xy_y_yyz_0,  \
                             ta2_xy_y_yyz_1,  \
                             ta2_xy_y_yz_0,   \
                             ta2_xy_y_yz_1,   \
                             ta2_xy_y_yzz_0,  \
                             ta2_xy_y_yzz_1,  \
                             ta2_xy_y_zzz_0,  \
                             ta2_xy_y_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xy_xxx_0[i] = ta1_x_x_xxx_1[i] + ta2_xy_x_xxx_0[i] * pa_y[i] - ta2_xy_x_xxx_1[i] * pc_y[i];

        ta2_xy_xy_xxy_0[i] = 2.0 * ta2_xy_y_xy_0[i] * fe_0 - 2.0 * ta2_xy_y_xy_1[i] * fe_0 + ta1_y_y_xxy_1[i] + ta2_xy_y_xxy_0[i] * pa_x[i] -
                             ta2_xy_y_xxy_1[i] * pc_x[i];

        ta2_xy_xy_xxz_0[i] = ta1_x_x_xxz_1[i] + ta2_xy_x_xxz_0[i] * pa_y[i] - ta2_xy_x_xxz_1[i] * pc_y[i];

        ta2_xy_xy_xyy_0[i] =
            ta2_xy_y_yy_0[i] * fe_0 - ta2_xy_y_yy_1[i] * fe_0 + ta1_y_y_xyy_1[i] + ta2_xy_y_xyy_0[i] * pa_x[i] - ta2_xy_y_xyy_1[i] * pc_x[i];

        ta2_xy_xy_xyz_0[i] =
            ta2_xy_y_yz_0[i] * fe_0 - ta2_xy_y_yz_1[i] * fe_0 + ta1_y_y_xyz_1[i] + ta2_xy_y_xyz_0[i] * pa_x[i] - ta2_xy_y_xyz_1[i] * pc_x[i];

        ta2_xy_xy_xzz_0[i] = ta1_x_x_xzz_1[i] + ta2_xy_x_xzz_0[i] * pa_y[i] - ta2_xy_x_xzz_1[i] * pc_y[i];

        ta2_xy_xy_yyy_0[i] = ta1_y_y_yyy_1[i] + ta2_xy_y_yyy_0[i] * pa_x[i] - ta2_xy_y_yyy_1[i] * pc_x[i];

        ta2_xy_xy_yyz_0[i] = ta1_y_y_yyz_1[i] + ta2_xy_y_yyz_0[i] * pa_x[i] - ta2_xy_y_yyz_1[i] * pc_x[i];

        ta2_xy_xy_yzz_0[i] = ta1_y_y_yzz_1[i] + ta2_xy_y_yzz_0[i] * pa_x[i] - ta2_xy_y_yzz_1[i] * pc_x[i];

        ta2_xy_xy_zzz_0[i] = ta1_y_y_zzz_1[i] + ta2_xy_y_zzz_0[i] * pa_x[i] - ta2_xy_y_zzz_1[i] * pc_x[i];
    }

    // Set up 80-90 components of targeted buffer : DF

    auto ta2_xy_xz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 80);

    auto ta2_xy_xz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 81);

    auto ta2_xy_xz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 82);

    auto ta2_xy_xz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 83);

    auto ta2_xy_xz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 84);

    auto ta2_xy_xz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 85);

    auto ta2_xy_xz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 86);

    auto ta2_xy_xz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 87);

    auto ta2_xy_xz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 88);

    auto ta2_xy_xz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 89);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_z_yyz_1,   \
                             ta1_y_z_yzz_1,   \
                             ta1_y_z_zzz_1,   \
                             ta2_xy_x_xx_0,   \
                             ta2_xy_x_xx_1,   \
                             ta2_xy_x_xxx_0,  \
                             ta2_xy_x_xxx_1,  \
                             ta2_xy_x_xxy_0,  \
                             ta2_xy_x_xxy_1,  \
                             ta2_xy_x_xxz_0,  \
                             ta2_xy_x_xxz_1,  \
                             ta2_xy_x_xy_0,   \
                             ta2_xy_x_xy_1,   \
                             ta2_xy_x_xyy_0,  \
                             ta2_xy_x_xyy_1,  \
                             ta2_xy_x_xyz_0,  \
                             ta2_xy_x_xyz_1,  \
                             ta2_xy_x_xz_0,   \
                             ta2_xy_x_xz_1,   \
                             ta2_xy_x_xzz_0,  \
                             ta2_xy_x_xzz_1,  \
                             ta2_xy_x_yyy_0,  \
                             ta2_xy_x_yyy_1,  \
                             ta2_xy_xz_xxx_0, \
                             ta2_xy_xz_xxy_0, \
                             ta2_xy_xz_xxz_0, \
                             ta2_xy_xz_xyy_0, \
                             ta2_xy_xz_xyz_0, \
                             ta2_xy_xz_xzz_0, \
                             ta2_xy_xz_yyy_0, \
                             ta2_xy_xz_yyz_0, \
                             ta2_xy_xz_yzz_0, \
                             ta2_xy_xz_zzz_0, \
                             ta2_xy_z_yyz_0,  \
                             ta2_xy_z_yyz_1,  \
                             ta2_xy_z_yzz_0,  \
                             ta2_xy_z_yzz_1,  \
                             ta2_xy_z_zzz_0,  \
                             ta2_xy_z_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xz_xxx_0[i] = ta2_xy_x_xxx_0[i] * pa_z[i] - ta2_xy_x_xxx_1[i] * pc_z[i];

        ta2_xy_xz_xxy_0[i] = ta2_xy_x_xxy_0[i] * pa_z[i] - ta2_xy_x_xxy_1[i] * pc_z[i];

        ta2_xy_xz_xxz_0[i] = ta2_xy_x_xx_0[i] * fe_0 - ta2_xy_x_xx_1[i] * fe_0 + ta2_xy_x_xxz_0[i] * pa_z[i] - ta2_xy_x_xxz_1[i] * pc_z[i];

        ta2_xy_xz_xyy_0[i] = ta2_xy_x_xyy_0[i] * pa_z[i] - ta2_xy_x_xyy_1[i] * pc_z[i];

        ta2_xy_xz_xyz_0[i] = ta2_xy_x_xy_0[i] * fe_0 - ta2_xy_x_xy_1[i] * fe_0 + ta2_xy_x_xyz_0[i] * pa_z[i] - ta2_xy_x_xyz_1[i] * pc_z[i];

        ta2_xy_xz_xzz_0[i] =
            2.0 * ta2_xy_x_xz_0[i] * fe_0 - 2.0 * ta2_xy_x_xz_1[i] * fe_0 + ta2_xy_x_xzz_0[i] * pa_z[i] - ta2_xy_x_xzz_1[i] * pc_z[i];

        ta2_xy_xz_yyy_0[i] = ta2_xy_x_yyy_0[i] * pa_z[i] - ta2_xy_x_yyy_1[i] * pc_z[i];

        ta2_xy_xz_yyz_0[i] = ta1_y_z_yyz_1[i] + ta2_xy_z_yyz_0[i] * pa_x[i] - ta2_xy_z_yyz_1[i] * pc_x[i];

        ta2_xy_xz_yzz_0[i] = ta1_y_z_yzz_1[i] + ta2_xy_z_yzz_0[i] * pa_x[i] - ta2_xy_z_yzz_1[i] * pc_x[i];

        ta2_xy_xz_zzz_0[i] = ta1_y_z_zzz_1[i] + ta2_xy_z_zzz_0[i] * pa_x[i] - ta2_xy_z_zzz_1[i] * pc_x[i];
    }

    // Set up 90-100 components of targeted buffer : DF

    auto ta2_xy_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 90);

    auto ta2_xy_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 91);

    auto ta2_xy_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 92);

    auto ta2_xy_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 93);

    auto ta2_xy_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 94);

    auto ta2_xy_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 95);

    auto ta2_xy_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 96);

    auto ta2_xy_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 97);

    auto ta2_xy_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 98);

    auto ta2_xy_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 99);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_x_y_xxx_1,   \
                             ta1_x_y_xxy_1,   \
                             ta1_x_y_xxz_1,   \
                             ta1_x_y_xyy_1,   \
                             ta1_x_y_xyz_1,   \
                             ta1_x_y_xzz_1,   \
                             ta1_x_y_yyy_1,   \
                             ta1_x_y_yyz_1,   \
                             ta1_x_y_yzz_1,   \
                             ta1_x_y_zzz_1,   \
                             ta2_xy_0_xxx_0,  \
                             ta2_xy_0_xxx_1,  \
                             ta2_xy_0_xxy_0,  \
                             ta2_xy_0_xxy_1,  \
                             ta2_xy_0_xxz_0,  \
                             ta2_xy_0_xxz_1,  \
                             ta2_xy_0_xyy_0,  \
                             ta2_xy_0_xyy_1,  \
                             ta2_xy_0_xyz_0,  \
                             ta2_xy_0_xyz_1,  \
                             ta2_xy_0_xzz_0,  \
                             ta2_xy_0_xzz_1,  \
                             ta2_xy_0_yyy_0,  \
                             ta2_xy_0_yyy_1,  \
                             ta2_xy_0_yyz_0,  \
                             ta2_xy_0_yyz_1,  \
                             ta2_xy_0_yzz_0,  \
                             ta2_xy_0_yzz_1,  \
                             ta2_xy_0_zzz_0,  \
                             ta2_xy_0_zzz_1,  \
                             ta2_xy_y_xx_0,   \
                             ta2_xy_y_xx_1,   \
                             ta2_xy_y_xxx_0,  \
                             ta2_xy_y_xxx_1,  \
                             ta2_xy_y_xxy_0,  \
                             ta2_xy_y_xxy_1,  \
                             ta2_xy_y_xxz_0,  \
                             ta2_xy_y_xxz_1,  \
                             ta2_xy_y_xy_0,   \
                             ta2_xy_y_xy_1,   \
                             ta2_xy_y_xyy_0,  \
                             ta2_xy_y_xyy_1,  \
                             ta2_xy_y_xyz_0,  \
                             ta2_xy_y_xyz_1,  \
                             ta2_xy_y_xz_0,   \
                             ta2_xy_y_xz_1,   \
                             ta2_xy_y_xzz_0,  \
                             ta2_xy_y_xzz_1,  \
                             ta2_xy_y_yy_0,   \
                             ta2_xy_y_yy_1,   \
                             ta2_xy_y_yyy_0,  \
                             ta2_xy_y_yyy_1,  \
                             ta2_xy_y_yyz_0,  \
                             ta2_xy_y_yyz_1,  \
                             ta2_xy_y_yz_0,   \
                             ta2_xy_y_yz_1,   \
                             ta2_xy_y_yzz_0,  \
                             ta2_xy_y_yzz_1,  \
                             ta2_xy_y_zz_0,   \
                             ta2_xy_y_zz_1,   \
                             ta2_xy_y_zzz_0,  \
                             ta2_xy_y_zzz_1,  \
                             ta2_xy_yy_xxx_0, \
                             ta2_xy_yy_xxy_0, \
                             ta2_xy_yy_xxz_0, \
                             ta2_xy_yy_xyy_0, \
                             ta2_xy_yy_xyz_0, \
                             ta2_xy_yy_xzz_0, \
                             ta2_xy_yy_yyy_0, \
                             ta2_xy_yy_yyz_0, \
                             ta2_xy_yy_yzz_0, \
                             ta2_xy_yy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yy_xxx_0[i] =
            ta2_xy_0_xxx_0[i] * fe_0 - ta2_xy_0_xxx_1[i] * fe_0 + ta1_x_y_xxx_1[i] + ta2_xy_y_xxx_0[i] * pa_y[i] - ta2_xy_y_xxx_1[i] * pc_y[i];

        ta2_xy_yy_xxy_0[i] = ta2_xy_0_xxy_0[i] * fe_0 - ta2_xy_0_xxy_1[i] * fe_0 + ta2_xy_y_xx_0[i] * fe_0 - ta2_xy_y_xx_1[i] * fe_0 +
                             ta1_x_y_xxy_1[i] + ta2_xy_y_xxy_0[i] * pa_y[i] - ta2_xy_y_xxy_1[i] * pc_y[i];

        ta2_xy_yy_xxz_0[i] =
            ta2_xy_0_xxz_0[i] * fe_0 - ta2_xy_0_xxz_1[i] * fe_0 + ta1_x_y_xxz_1[i] + ta2_xy_y_xxz_0[i] * pa_y[i] - ta2_xy_y_xxz_1[i] * pc_y[i];

        ta2_xy_yy_xyy_0[i] = ta2_xy_0_xyy_0[i] * fe_0 - ta2_xy_0_xyy_1[i] * fe_0 + 2.0 * ta2_xy_y_xy_0[i] * fe_0 - 2.0 * ta2_xy_y_xy_1[i] * fe_0 +
                             ta1_x_y_xyy_1[i] + ta2_xy_y_xyy_0[i] * pa_y[i] - ta2_xy_y_xyy_1[i] * pc_y[i];

        ta2_xy_yy_xyz_0[i] = ta2_xy_0_xyz_0[i] * fe_0 - ta2_xy_0_xyz_1[i] * fe_0 + ta2_xy_y_xz_0[i] * fe_0 - ta2_xy_y_xz_1[i] * fe_0 +
                             ta1_x_y_xyz_1[i] + ta2_xy_y_xyz_0[i] * pa_y[i] - ta2_xy_y_xyz_1[i] * pc_y[i];

        ta2_xy_yy_xzz_0[i] =
            ta2_xy_0_xzz_0[i] * fe_0 - ta2_xy_0_xzz_1[i] * fe_0 + ta1_x_y_xzz_1[i] + ta2_xy_y_xzz_0[i] * pa_y[i] - ta2_xy_y_xzz_1[i] * pc_y[i];

        ta2_xy_yy_yyy_0[i] = ta2_xy_0_yyy_0[i] * fe_0 - ta2_xy_0_yyy_1[i] * fe_0 + 3.0 * ta2_xy_y_yy_0[i] * fe_0 - 3.0 * ta2_xy_y_yy_1[i] * fe_0 +
                             ta1_x_y_yyy_1[i] + ta2_xy_y_yyy_0[i] * pa_y[i] - ta2_xy_y_yyy_1[i] * pc_y[i];

        ta2_xy_yy_yyz_0[i] = ta2_xy_0_yyz_0[i] * fe_0 - ta2_xy_0_yyz_1[i] * fe_0 + 2.0 * ta2_xy_y_yz_0[i] * fe_0 - 2.0 * ta2_xy_y_yz_1[i] * fe_0 +
                             ta1_x_y_yyz_1[i] + ta2_xy_y_yyz_0[i] * pa_y[i] - ta2_xy_y_yyz_1[i] * pc_y[i];

        ta2_xy_yy_yzz_0[i] = ta2_xy_0_yzz_0[i] * fe_0 - ta2_xy_0_yzz_1[i] * fe_0 + ta2_xy_y_zz_0[i] * fe_0 - ta2_xy_y_zz_1[i] * fe_0 +
                             ta1_x_y_yzz_1[i] + ta2_xy_y_yzz_0[i] * pa_y[i] - ta2_xy_y_yzz_1[i] * pc_y[i];

        ta2_xy_yy_zzz_0[i] =
            ta2_xy_0_zzz_0[i] * fe_0 - ta2_xy_0_zzz_1[i] * fe_0 + ta1_x_y_zzz_1[i] + ta2_xy_y_zzz_0[i] * pa_y[i] - ta2_xy_y_zzz_1[i] * pc_y[i];
    }

    // Set up 100-110 components of targeted buffer : DF

    auto ta2_xy_yz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 100);

    auto ta2_xy_yz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 101);

    auto ta2_xy_yz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 102);

    auto ta2_xy_yz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 103);

    auto ta2_xy_yz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 104);

    auto ta2_xy_yz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 105);

    auto ta2_xy_yz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 106);

    auto ta2_xy_yz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 107);

    auto ta2_xy_yz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 108);

    auto ta2_xy_yz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 109);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_z_xxz_1,   \
                             ta1_x_z_xzz_1,   \
                             ta1_x_z_zzz_1,   \
                             ta2_xy_y_xxx_0,  \
                             ta2_xy_y_xxx_1,  \
                             ta2_xy_y_xxy_0,  \
                             ta2_xy_y_xxy_1,  \
                             ta2_xy_y_xy_0,   \
                             ta2_xy_y_xy_1,   \
                             ta2_xy_y_xyy_0,  \
                             ta2_xy_y_xyy_1,  \
                             ta2_xy_y_xyz_0,  \
                             ta2_xy_y_xyz_1,  \
                             ta2_xy_y_yy_0,   \
                             ta2_xy_y_yy_1,   \
                             ta2_xy_y_yyy_0,  \
                             ta2_xy_y_yyy_1,  \
                             ta2_xy_y_yyz_0,  \
                             ta2_xy_y_yyz_1,  \
                             ta2_xy_y_yz_0,   \
                             ta2_xy_y_yz_1,   \
                             ta2_xy_y_yzz_0,  \
                             ta2_xy_y_yzz_1,  \
                             ta2_xy_yz_xxx_0, \
                             ta2_xy_yz_xxy_0, \
                             ta2_xy_yz_xxz_0, \
                             ta2_xy_yz_xyy_0, \
                             ta2_xy_yz_xyz_0, \
                             ta2_xy_yz_xzz_0, \
                             ta2_xy_yz_yyy_0, \
                             ta2_xy_yz_yyz_0, \
                             ta2_xy_yz_yzz_0, \
                             ta2_xy_yz_zzz_0, \
                             ta2_xy_z_xxz_0,  \
                             ta2_xy_z_xxz_1,  \
                             ta2_xy_z_xzz_0,  \
                             ta2_xy_z_xzz_1,  \
                             ta2_xy_z_zzz_0,  \
                             ta2_xy_z_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yz_xxx_0[i] = ta2_xy_y_xxx_0[i] * pa_z[i] - ta2_xy_y_xxx_1[i] * pc_z[i];

        ta2_xy_yz_xxy_0[i] = ta2_xy_y_xxy_0[i] * pa_z[i] - ta2_xy_y_xxy_1[i] * pc_z[i];

        ta2_xy_yz_xxz_0[i] = ta1_x_z_xxz_1[i] + ta2_xy_z_xxz_0[i] * pa_y[i] - ta2_xy_z_xxz_1[i] * pc_y[i];

        ta2_xy_yz_xyy_0[i] = ta2_xy_y_xyy_0[i] * pa_z[i] - ta2_xy_y_xyy_1[i] * pc_z[i];

        ta2_xy_yz_xyz_0[i] = ta2_xy_y_xy_0[i] * fe_0 - ta2_xy_y_xy_1[i] * fe_0 + ta2_xy_y_xyz_0[i] * pa_z[i] - ta2_xy_y_xyz_1[i] * pc_z[i];

        ta2_xy_yz_xzz_0[i] = ta1_x_z_xzz_1[i] + ta2_xy_z_xzz_0[i] * pa_y[i] - ta2_xy_z_xzz_1[i] * pc_y[i];

        ta2_xy_yz_yyy_0[i] = ta2_xy_y_yyy_0[i] * pa_z[i] - ta2_xy_y_yyy_1[i] * pc_z[i];

        ta2_xy_yz_yyz_0[i] = ta2_xy_y_yy_0[i] * fe_0 - ta2_xy_y_yy_1[i] * fe_0 + ta2_xy_y_yyz_0[i] * pa_z[i] - ta2_xy_y_yyz_1[i] * pc_z[i];

        ta2_xy_yz_yzz_0[i] =
            2.0 * ta2_xy_y_yz_0[i] * fe_0 - 2.0 * ta2_xy_y_yz_1[i] * fe_0 + ta2_xy_y_yzz_0[i] * pa_z[i] - ta2_xy_y_yzz_1[i] * pc_z[i];

        ta2_xy_yz_zzz_0[i] = ta1_x_z_zzz_1[i] + ta2_xy_z_zzz_0[i] * pa_y[i] - ta2_xy_z_zzz_1[i] * pc_y[i];
    }

    // Set up 110-120 components of targeted buffer : DF

    auto ta2_xy_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 110);

    auto ta2_xy_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 111);

    auto ta2_xy_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 112);

    auto ta2_xy_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 113);

    auto ta2_xy_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 114);

    auto ta2_xy_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 115);

    auto ta2_xy_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 116);

    auto ta2_xy_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 117);

    auto ta2_xy_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 118);

    auto ta2_xy_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 119);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta2_xy_0_xxx_0,  \
                             ta2_xy_0_xxx_1,  \
                             ta2_xy_0_xxy_0,  \
                             ta2_xy_0_xxy_1,  \
                             ta2_xy_0_xxz_0,  \
                             ta2_xy_0_xxz_1,  \
                             ta2_xy_0_xyy_0,  \
                             ta2_xy_0_xyy_1,  \
                             ta2_xy_0_xyz_0,  \
                             ta2_xy_0_xyz_1,  \
                             ta2_xy_0_xzz_0,  \
                             ta2_xy_0_xzz_1,  \
                             ta2_xy_0_yyy_0,  \
                             ta2_xy_0_yyy_1,  \
                             ta2_xy_0_yyz_0,  \
                             ta2_xy_0_yyz_1,  \
                             ta2_xy_0_yzz_0,  \
                             ta2_xy_0_yzz_1,  \
                             ta2_xy_0_zzz_0,  \
                             ta2_xy_0_zzz_1,  \
                             ta2_xy_z_xx_0,   \
                             ta2_xy_z_xx_1,   \
                             ta2_xy_z_xxx_0,  \
                             ta2_xy_z_xxx_1,  \
                             ta2_xy_z_xxy_0,  \
                             ta2_xy_z_xxy_1,  \
                             ta2_xy_z_xxz_0,  \
                             ta2_xy_z_xxz_1,  \
                             ta2_xy_z_xy_0,   \
                             ta2_xy_z_xy_1,   \
                             ta2_xy_z_xyy_0,  \
                             ta2_xy_z_xyy_1,  \
                             ta2_xy_z_xyz_0,  \
                             ta2_xy_z_xyz_1,  \
                             ta2_xy_z_xz_0,   \
                             ta2_xy_z_xz_1,   \
                             ta2_xy_z_xzz_0,  \
                             ta2_xy_z_xzz_1,  \
                             ta2_xy_z_yy_0,   \
                             ta2_xy_z_yy_1,   \
                             ta2_xy_z_yyy_0,  \
                             ta2_xy_z_yyy_1,  \
                             ta2_xy_z_yyz_0,  \
                             ta2_xy_z_yyz_1,  \
                             ta2_xy_z_yz_0,   \
                             ta2_xy_z_yz_1,   \
                             ta2_xy_z_yzz_0,  \
                             ta2_xy_z_yzz_1,  \
                             ta2_xy_z_zz_0,   \
                             ta2_xy_z_zz_1,   \
                             ta2_xy_z_zzz_0,  \
                             ta2_xy_z_zzz_1,  \
                             ta2_xy_zz_xxx_0, \
                             ta2_xy_zz_xxy_0, \
                             ta2_xy_zz_xxz_0, \
                             ta2_xy_zz_xyy_0, \
                             ta2_xy_zz_xyz_0, \
                             ta2_xy_zz_xzz_0, \
                             ta2_xy_zz_yyy_0, \
                             ta2_xy_zz_yyz_0, \
                             ta2_xy_zz_yzz_0, \
                             ta2_xy_zz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_zz_xxx_0[i] = ta2_xy_0_xxx_0[i] * fe_0 - ta2_xy_0_xxx_1[i] * fe_0 + ta2_xy_z_xxx_0[i] * pa_z[i] - ta2_xy_z_xxx_1[i] * pc_z[i];

        ta2_xy_zz_xxy_0[i] = ta2_xy_0_xxy_0[i] * fe_0 - ta2_xy_0_xxy_1[i] * fe_0 + ta2_xy_z_xxy_0[i] * pa_z[i] - ta2_xy_z_xxy_1[i] * pc_z[i];

        ta2_xy_zz_xxz_0[i] = ta2_xy_0_xxz_0[i] * fe_0 - ta2_xy_0_xxz_1[i] * fe_0 + ta2_xy_z_xx_0[i] * fe_0 - ta2_xy_z_xx_1[i] * fe_0 +
                             ta2_xy_z_xxz_0[i] * pa_z[i] - ta2_xy_z_xxz_1[i] * pc_z[i];

        ta2_xy_zz_xyy_0[i] = ta2_xy_0_xyy_0[i] * fe_0 - ta2_xy_0_xyy_1[i] * fe_0 + ta2_xy_z_xyy_0[i] * pa_z[i] - ta2_xy_z_xyy_1[i] * pc_z[i];

        ta2_xy_zz_xyz_0[i] = ta2_xy_0_xyz_0[i] * fe_0 - ta2_xy_0_xyz_1[i] * fe_0 + ta2_xy_z_xy_0[i] * fe_0 - ta2_xy_z_xy_1[i] * fe_0 +
                             ta2_xy_z_xyz_0[i] * pa_z[i] - ta2_xy_z_xyz_1[i] * pc_z[i];

        ta2_xy_zz_xzz_0[i] = ta2_xy_0_xzz_0[i] * fe_0 - ta2_xy_0_xzz_1[i] * fe_0 + 2.0 * ta2_xy_z_xz_0[i] * fe_0 - 2.0 * ta2_xy_z_xz_1[i] * fe_0 +
                             ta2_xy_z_xzz_0[i] * pa_z[i] - ta2_xy_z_xzz_1[i] * pc_z[i];

        ta2_xy_zz_yyy_0[i] = ta2_xy_0_yyy_0[i] * fe_0 - ta2_xy_0_yyy_1[i] * fe_0 + ta2_xy_z_yyy_0[i] * pa_z[i] - ta2_xy_z_yyy_1[i] * pc_z[i];

        ta2_xy_zz_yyz_0[i] = ta2_xy_0_yyz_0[i] * fe_0 - ta2_xy_0_yyz_1[i] * fe_0 + ta2_xy_z_yy_0[i] * fe_0 - ta2_xy_z_yy_1[i] * fe_0 +
                             ta2_xy_z_yyz_0[i] * pa_z[i] - ta2_xy_z_yyz_1[i] * pc_z[i];

        ta2_xy_zz_yzz_0[i] = ta2_xy_0_yzz_0[i] * fe_0 - ta2_xy_0_yzz_1[i] * fe_0 + 2.0 * ta2_xy_z_yz_0[i] * fe_0 - 2.0 * ta2_xy_z_yz_1[i] * fe_0 +
                             ta2_xy_z_yzz_0[i] * pa_z[i] - ta2_xy_z_yzz_1[i] * pc_z[i];

        ta2_xy_zz_zzz_0[i] = ta2_xy_0_zzz_0[i] * fe_0 - ta2_xy_0_zzz_1[i] * fe_0 + 3.0 * ta2_xy_z_zz_0[i] * fe_0 - 3.0 * ta2_xy_z_zz_1[i] * fe_0 +
                             ta2_xy_z_zzz_0[i] * pa_z[i] - ta2_xy_z_zzz_1[i] * pc_z[i];
    }

    // Set up 120-130 components of targeted buffer : DF

    auto ta2_xz_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 120);

    auto ta2_xz_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 121);

    auto ta2_xz_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 122);

    auto ta2_xz_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 123);

    auto ta2_xz_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 124);

    auto ta2_xz_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 125);

    auto ta2_xz_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 126);

    auto ta2_xz_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 127);

    auto ta2_xz_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 128);

    auto ta2_xz_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 129);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta1_z_x_xxx_1,   \
                             ta1_z_x_xxy_1,   \
                             ta1_z_x_xxz_1,   \
                             ta1_z_x_xyy_1,   \
                             ta1_z_x_xyz_1,   \
                             ta1_z_x_xzz_1,   \
                             ta1_z_x_yyy_1,   \
                             ta1_z_x_yyz_1,   \
                             ta1_z_x_yzz_1,   \
                             ta1_z_x_zzz_1,   \
                             ta2_xz_0_xxx_0,  \
                             ta2_xz_0_xxx_1,  \
                             ta2_xz_0_xxy_0,  \
                             ta2_xz_0_xxy_1,  \
                             ta2_xz_0_xxz_0,  \
                             ta2_xz_0_xxz_1,  \
                             ta2_xz_0_xyy_0,  \
                             ta2_xz_0_xyy_1,  \
                             ta2_xz_0_xyz_0,  \
                             ta2_xz_0_xyz_1,  \
                             ta2_xz_0_xzz_0,  \
                             ta2_xz_0_xzz_1,  \
                             ta2_xz_0_yyy_0,  \
                             ta2_xz_0_yyy_1,  \
                             ta2_xz_0_yyz_0,  \
                             ta2_xz_0_yyz_1,  \
                             ta2_xz_0_yzz_0,  \
                             ta2_xz_0_yzz_1,  \
                             ta2_xz_0_zzz_0,  \
                             ta2_xz_0_zzz_1,  \
                             ta2_xz_x_xx_0,   \
                             ta2_xz_x_xx_1,   \
                             ta2_xz_x_xxx_0,  \
                             ta2_xz_x_xxx_1,  \
                             ta2_xz_x_xxy_0,  \
                             ta2_xz_x_xxy_1,  \
                             ta2_xz_x_xxz_0,  \
                             ta2_xz_x_xxz_1,  \
                             ta2_xz_x_xy_0,   \
                             ta2_xz_x_xy_1,   \
                             ta2_xz_x_xyy_0,  \
                             ta2_xz_x_xyy_1,  \
                             ta2_xz_x_xyz_0,  \
                             ta2_xz_x_xyz_1,  \
                             ta2_xz_x_xz_0,   \
                             ta2_xz_x_xz_1,   \
                             ta2_xz_x_xzz_0,  \
                             ta2_xz_x_xzz_1,  \
                             ta2_xz_x_yy_0,   \
                             ta2_xz_x_yy_1,   \
                             ta2_xz_x_yyy_0,  \
                             ta2_xz_x_yyy_1,  \
                             ta2_xz_x_yyz_0,  \
                             ta2_xz_x_yyz_1,  \
                             ta2_xz_x_yz_0,   \
                             ta2_xz_x_yz_1,   \
                             ta2_xz_x_yzz_0,  \
                             ta2_xz_x_yzz_1,  \
                             ta2_xz_x_zz_0,   \
                             ta2_xz_x_zz_1,   \
                             ta2_xz_x_zzz_0,  \
                             ta2_xz_x_zzz_1,  \
                             ta2_xz_xx_xxx_0, \
                             ta2_xz_xx_xxy_0, \
                             ta2_xz_xx_xxz_0, \
                             ta2_xz_xx_xyy_0, \
                             ta2_xz_xx_xyz_0, \
                             ta2_xz_xx_xzz_0, \
                             ta2_xz_xx_yyy_0, \
                             ta2_xz_xx_yyz_0, \
                             ta2_xz_xx_yzz_0, \
                             ta2_xz_xx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xx_xxx_0[i] = ta2_xz_0_xxx_0[i] * fe_0 - ta2_xz_0_xxx_1[i] * fe_0 + 3.0 * ta2_xz_x_xx_0[i] * fe_0 - 3.0 * ta2_xz_x_xx_1[i] * fe_0 +
                             ta1_z_x_xxx_1[i] + ta2_xz_x_xxx_0[i] * pa_x[i] - ta2_xz_x_xxx_1[i] * pc_x[i];

        ta2_xz_xx_xxy_0[i] = ta2_xz_0_xxy_0[i] * fe_0 - ta2_xz_0_xxy_1[i] * fe_0 + 2.0 * ta2_xz_x_xy_0[i] * fe_0 - 2.0 * ta2_xz_x_xy_1[i] * fe_0 +
                             ta1_z_x_xxy_1[i] + ta2_xz_x_xxy_0[i] * pa_x[i] - ta2_xz_x_xxy_1[i] * pc_x[i];

        ta2_xz_xx_xxz_0[i] = ta2_xz_0_xxz_0[i] * fe_0 - ta2_xz_0_xxz_1[i] * fe_0 + 2.0 * ta2_xz_x_xz_0[i] * fe_0 - 2.0 * ta2_xz_x_xz_1[i] * fe_0 +
                             ta1_z_x_xxz_1[i] + ta2_xz_x_xxz_0[i] * pa_x[i] - ta2_xz_x_xxz_1[i] * pc_x[i];

        ta2_xz_xx_xyy_0[i] = ta2_xz_0_xyy_0[i] * fe_0 - ta2_xz_0_xyy_1[i] * fe_0 + ta2_xz_x_yy_0[i] * fe_0 - ta2_xz_x_yy_1[i] * fe_0 +
                             ta1_z_x_xyy_1[i] + ta2_xz_x_xyy_0[i] * pa_x[i] - ta2_xz_x_xyy_1[i] * pc_x[i];

        ta2_xz_xx_xyz_0[i] = ta2_xz_0_xyz_0[i] * fe_0 - ta2_xz_0_xyz_1[i] * fe_0 + ta2_xz_x_yz_0[i] * fe_0 - ta2_xz_x_yz_1[i] * fe_0 +
                             ta1_z_x_xyz_1[i] + ta2_xz_x_xyz_0[i] * pa_x[i] - ta2_xz_x_xyz_1[i] * pc_x[i];

        ta2_xz_xx_xzz_0[i] = ta2_xz_0_xzz_0[i] * fe_0 - ta2_xz_0_xzz_1[i] * fe_0 + ta2_xz_x_zz_0[i] * fe_0 - ta2_xz_x_zz_1[i] * fe_0 +
                             ta1_z_x_xzz_1[i] + ta2_xz_x_xzz_0[i] * pa_x[i] - ta2_xz_x_xzz_1[i] * pc_x[i];

        ta2_xz_xx_yyy_0[i] =
            ta2_xz_0_yyy_0[i] * fe_0 - ta2_xz_0_yyy_1[i] * fe_0 + ta1_z_x_yyy_1[i] + ta2_xz_x_yyy_0[i] * pa_x[i] - ta2_xz_x_yyy_1[i] * pc_x[i];

        ta2_xz_xx_yyz_0[i] =
            ta2_xz_0_yyz_0[i] * fe_0 - ta2_xz_0_yyz_1[i] * fe_0 + ta1_z_x_yyz_1[i] + ta2_xz_x_yyz_0[i] * pa_x[i] - ta2_xz_x_yyz_1[i] * pc_x[i];

        ta2_xz_xx_yzz_0[i] =
            ta2_xz_0_yzz_0[i] * fe_0 - ta2_xz_0_yzz_1[i] * fe_0 + ta1_z_x_yzz_1[i] + ta2_xz_x_yzz_0[i] * pa_x[i] - ta2_xz_x_yzz_1[i] * pc_x[i];

        ta2_xz_xx_zzz_0[i] =
            ta2_xz_0_zzz_0[i] * fe_0 - ta2_xz_0_zzz_1[i] * fe_0 + ta1_z_x_zzz_1[i] + ta2_xz_x_zzz_0[i] * pa_x[i] - ta2_xz_x_zzz_1[i] * pc_x[i];
    }

    // Set up 130-140 components of targeted buffer : DF

    auto ta2_xz_xy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 130);

    auto ta2_xz_xy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 131);

    auto ta2_xz_xy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 132);

    auto ta2_xz_xy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 133);

    auto ta2_xz_xy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 134);

    auto ta2_xz_xy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 135);

    auto ta2_xz_xy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 136);

    auto ta2_xz_xy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 137);

    auto ta2_xz_xy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 138);

    auto ta2_xz_xy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 139);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_y_yyy_1,   \
                             ta1_z_y_yyz_1,   \
                             ta1_z_y_yzz_1,   \
                             ta2_xz_x_xx_0,   \
                             ta2_xz_x_xx_1,   \
                             ta2_xz_x_xxx_0,  \
                             ta2_xz_x_xxx_1,  \
                             ta2_xz_x_xxy_0,  \
                             ta2_xz_x_xxy_1,  \
                             ta2_xz_x_xxz_0,  \
                             ta2_xz_x_xxz_1,  \
                             ta2_xz_x_xy_0,   \
                             ta2_xz_x_xy_1,   \
                             ta2_xz_x_xyy_0,  \
                             ta2_xz_x_xyy_1,  \
                             ta2_xz_x_xyz_0,  \
                             ta2_xz_x_xyz_1,  \
                             ta2_xz_x_xz_0,   \
                             ta2_xz_x_xz_1,   \
                             ta2_xz_x_xzz_0,  \
                             ta2_xz_x_xzz_1,  \
                             ta2_xz_x_zzz_0,  \
                             ta2_xz_x_zzz_1,  \
                             ta2_xz_xy_xxx_0, \
                             ta2_xz_xy_xxy_0, \
                             ta2_xz_xy_xxz_0, \
                             ta2_xz_xy_xyy_0, \
                             ta2_xz_xy_xyz_0, \
                             ta2_xz_xy_xzz_0, \
                             ta2_xz_xy_yyy_0, \
                             ta2_xz_xy_yyz_0, \
                             ta2_xz_xy_yzz_0, \
                             ta2_xz_xy_zzz_0, \
                             ta2_xz_y_yyy_0,  \
                             ta2_xz_y_yyy_1,  \
                             ta2_xz_y_yyz_0,  \
                             ta2_xz_y_yyz_1,  \
                             ta2_xz_y_yzz_0,  \
                             ta2_xz_y_yzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xy_xxx_0[i] = ta2_xz_x_xxx_0[i] * pa_y[i] - ta2_xz_x_xxx_1[i] * pc_y[i];

        ta2_xz_xy_xxy_0[i] = ta2_xz_x_xx_0[i] * fe_0 - ta2_xz_x_xx_1[i] * fe_0 + ta2_xz_x_xxy_0[i] * pa_y[i] - ta2_xz_x_xxy_1[i] * pc_y[i];

        ta2_xz_xy_xxz_0[i] = ta2_xz_x_xxz_0[i] * pa_y[i] - ta2_xz_x_xxz_1[i] * pc_y[i];

        ta2_xz_xy_xyy_0[i] =
            2.0 * ta2_xz_x_xy_0[i] * fe_0 - 2.0 * ta2_xz_x_xy_1[i] * fe_0 + ta2_xz_x_xyy_0[i] * pa_y[i] - ta2_xz_x_xyy_1[i] * pc_y[i];

        ta2_xz_xy_xyz_0[i] = ta2_xz_x_xz_0[i] * fe_0 - ta2_xz_x_xz_1[i] * fe_0 + ta2_xz_x_xyz_0[i] * pa_y[i] - ta2_xz_x_xyz_1[i] * pc_y[i];

        ta2_xz_xy_xzz_0[i] = ta2_xz_x_xzz_0[i] * pa_y[i] - ta2_xz_x_xzz_1[i] * pc_y[i];

        ta2_xz_xy_yyy_0[i] = ta1_z_y_yyy_1[i] + ta2_xz_y_yyy_0[i] * pa_x[i] - ta2_xz_y_yyy_1[i] * pc_x[i];

        ta2_xz_xy_yyz_0[i] = ta1_z_y_yyz_1[i] + ta2_xz_y_yyz_0[i] * pa_x[i] - ta2_xz_y_yyz_1[i] * pc_x[i];

        ta2_xz_xy_yzz_0[i] = ta1_z_y_yzz_1[i] + ta2_xz_y_yzz_0[i] * pa_x[i] - ta2_xz_y_yzz_1[i] * pc_x[i];

        ta2_xz_xy_zzz_0[i] = ta2_xz_x_zzz_0[i] * pa_y[i] - ta2_xz_x_zzz_1[i] * pc_y[i];
    }

    // Set up 140-150 components of targeted buffer : DF

    auto ta2_xz_xz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 140);

    auto ta2_xz_xz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 141);

    auto ta2_xz_xz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 142);

    auto ta2_xz_xz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 143);

    auto ta2_xz_xz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 144);

    auto ta2_xz_xz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 145);

    auto ta2_xz_xz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 146);

    auto ta2_xz_xz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 147);

    auto ta2_xz_xz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 148);

    auto ta2_xz_xz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 149);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_x_x_xxx_1,   \
                             ta1_x_x_xxy_1,   \
                             ta1_x_x_xyy_1,   \
                             ta1_z_z_xxz_1,   \
                             ta1_z_z_xyz_1,   \
                             ta1_z_z_xzz_1,   \
                             ta1_z_z_yyy_1,   \
                             ta1_z_z_yyz_1,   \
                             ta1_z_z_yzz_1,   \
                             ta1_z_z_zzz_1,   \
                             ta2_xz_x_xxx_0,  \
                             ta2_xz_x_xxx_1,  \
                             ta2_xz_x_xxy_0,  \
                             ta2_xz_x_xxy_1,  \
                             ta2_xz_x_xyy_0,  \
                             ta2_xz_x_xyy_1,  \
                             ta2_xz_xz_xxx_0, \
                             ta2_xz_xz_xxy_0, \
                             ta2_xz_xz_xxz_0, \
                             ta2_xz_xz_xyy_0, \
                             ta2_xz_xz_xyz_0, \
                             ta2_xz_xz_xzz_0, \
                             ta2_xz_xz_yyy_0, \
                             ta2_xz_xz_yyz_0, \
                             ta2_xz_xz_yzz_0, \
                             ta2_xz_xz_zzz_0, \
                             ta2_xz_z_xxz_0,  \
                             ta2_xz_z_xxz_1,  \
                             ta2_xz_z_xyz_0,  \
                             ta2_xz_z_xyz_1,  \
                             ta2_xz_z_xz_0,   \
                             ta2_xz_z_xz_1,   \
                             ta2_xz_z_xzz_0,  \
                             ta2_xz_z_xzz_1,  \
                             ta2_xz_z_yyy_0,  \
                             ta2_xz_z_yyy_1,  \
                             ta2_xz_z_yyz_0,  \
                             ta2_xz_z_yyz_1,  \
                             ta2_xz_z_yz_0,   \
                             ta2_xz_z_yz_1,   \
                             ta2_xz_z_yzz_0,  \
                             ta2_xz_z_yzz_1,  \
                             ta2_xz_z_zz_0,   \
                             ta2_xz_z_zz_1,   \
                             ta2_xz_z_zzz_0,  \
                             ta2_xz_z_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xz_xxx_0[i] = ta1_x_x_xxx_1[i] + ta2_xz_x_xxx_0[i] * pa_z[i] - ta2_xz_x_xxx_1[i] * pc_z[i];

        ta2_xz_xz_xxy_0[i] = ta1_x_x_xxy_1[i] + ta2_xz_x_xxy_0[i] * pa_z[i] - ta2_xz_x_xxy_1[i] * pc_z[i];

        ta2_xz_xz_xxz_0[i] = 2.0 * ta2_xz_z_xz_0[i] * fe_0 - 2.0 * ta2_xz_z_xz_1[i] * fe_0 + ta1_z_z_xxz_1[i] + ta2_xz_z_xxz_0[i] * pa_x[i] -
                             ta2_xz_z_xxz_1[i] * pc_x[i];

        ta2_xz_xz_xyy_0[i] = ta1_x_x_xyy_1[i] + ta2_xz_x_xyy_0[i] * pa_z[i] - ta2_xz_x_xyy_1[i] * pc_z[i];

        ta2_xz_xz_xyz_0[i] =
            ta2_xz_z_yz_0[i] * fe_0 - ta2_xz_z_yz_1[i] * fe_0 + ta1_z_z_xyz_1[i] + ta2_xz_z_xyz_0[i] * pa_x[i] - ta2_xz_z_xyz_1[i] * pc_x[i];

        ta2_xz_xz_xzz_0[i] =
            ta2_xz_z_zz_0[i] * fe_0 - ta2_xz_z_zz_1[i] * fe_0 + ta1_z_z_xzz_1[i] + ta2_xz_z_xzz_0[i] * pa_x[i] - ta2_xz_z_xzz_1[i] * pc_x[i];

        ta2_xz_xz_yyy_0[i] = ta1_z_z_yyy_1[i] + ta2_xz_z_yyy_0[i] * pa_x[i] - ta2_xz_z_yyy_1[i] * pc_x[i];

        ta2_xz_xz_yyz_0[i] = ta1_z_z_yyz_1[i] + ta2_xz_z_yyz_0[i] * pa_x[i] - ta2_xz_z_yyz_1[i] * pc_x[i];

        ta2_xz_xz_yzz_0[i] = ta1_z_z_yzz_1[i] + ta2_xz_z_yzz_0[i] * pa_x[i] - ta2_xz_z_yzz_1[i] * pc_x[i];

        ta2_xz_xz_zzz_0[i] = ta1_z_z_zzz_1[i] + ta2_xz_z_zzz_0[i] * pa_x[i] - ta2_xz_z_zzz_1[i] * pc_x[i];
    }

    // Set up 150-160 components of targeted buffer : DF

    auto ta2_xz_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 150);

    auto ta2_xz_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 151);

    auto ta2_xz_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 152);

    auto ta2_xz_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 153);

    auto ta2_xz_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 154);

    auto ta2_xz_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 155);

    auto ta2_xz_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 156);

    auto ta2_xz_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 157);

    auto ta2_xz_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 158);

    auto ta2_xz_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 159);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta2_xz_0_xxx_0,  \
                             ta2_xz_0_xxx_1,  \
                             ta2_xz_0_xxy_0,  \
                             ta2_xz_0_xxy_1,  \
                             ta2_xz_0_xxz_0,  \
                             ta2_xz_0_xxz_1,  \
                             ta2_xz_0_xyy_0,  \
                             ta2_xz_0_xyy_1,  \
                             ta2_xz_0_xyz_0,  \
                             ta2_xz_0_xyz_1,  \
                             ta2_xz_0_xzz_0,  \
                             ta2_xz_0_xzz_1,  \
                             ta2_xz_0_yyy_0,  \
                             ta2_xz_0_yyy_1,  \
                             ta2_xz_0_yyz_0,  \
                             ta2_xz_0_yyz_1,  \
                             ta2_xz_0_yzz_0,  \
                             ta2_xz_0_yzz_1,  \
                             ta2_xz_0_zzz_0,  \
                             ta2_xz_0_zzz_1,  \
                             ta2_xz_y_xx_0,   \
                             ta2_xz_y_xx_1,   \
                             ta2_xz_y_xxx_0,  \
                             ta2_xz_y_xxx_1,  \
                             ta2_xz_y_xxy_0,  \
                             ta2_xz_y_xxy_1,  \
                             ta2_xz_y_xxz_0,  \
                             ta2_xz_y_xxz_1,  \
                             ta2_xz_y_xy_0,   \
                             ta2_xz_y_xy_1,   \
                             ta2_xz_y_xyy_0,  \
                             ta2_xz_y_xyy_1,  \
                             ta2_xz_y_xyz_0,  \
                             ta2_xz_y_xyz_1,  \
                             ta2_xz_y_xz_0,   \
                             ta2_xz_y_xz_1,   \
                             ta2_xz_y_xzz_0,  \
                             ta2_xz_y_xzz_1,  \
                             ta2_xz_y_yy_0,   \
                             ta2_xz_y_yy_1,   \
                             ta2_xz_y_yyy_0,  \
                             ta2_xz_y_yyy_1,  \
                             ta2_xz_y_yyz_0,  \
                             ta2_xz_y_yyz_1,  \
                             ta2_xz_y_yz_0,   \
                             ta2_xz_y_yz_1,   \
                             ta2_xz_y_yzz_0,  \
                             ta2_xz_y_yzz_1,  \
                             ta2_xz_y_zz_0,   \
                             ta2_xz_y_zz_1,   \
                             ta2_xz_y_zzz_0,  \
                             ta2_xz_y_zzz_1,  \
                             ta2_xz_yy_xxx_0, \
                             ta2_xz_yy_xxy_0, \
                             ta2_xz_yy_xxz_0, \
                             ta2_xz_yy_xyy_0, \
                             ta2_xz_yy_xyz_0, \
                             ta2_xz_yy_xzz_0, \
                             ta2_xz_yy_yyy_0, \
                             ta2_xz_yy_yyz_0, \
                             ta2_xz_yy_yzz_0, \
                             ta2_xz_yy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yy_xxx_0[i] = ta2_xz_0_xxx_0[i] * fe_0 - ta2_xz_0_xxx_1[i] * fe_0 + ta2_xz_y_xxx_0[i] * pa_y[i] - ta2_xz_y_xxx_1[i] * pc_y[i];

        ta2_xz_yy_xxy_0[i] = ta2_xz_0_xxy_0[i] * fe_0 - ta2_xz_0_xxy_1[i] * fe_0 + ta2_xz_y_xx_0[i] * fe_0 - ta2_xz_y_xx_1[i] * fe_0 +
                             ta2_xz_y_xxy_0[i] * pa_y[i] - ta2_xz_y_xxy_1[i] * pc_y[i];

        ta2_xz_yy_xxz_0[i] = ta2_xz_0_xxz_0[i] * fe_0 - ta2_xz_0_xxz_1[i] * fe_0 + ta2_xz_y_xxz_0[i] * pa_y[i] - ta2_xz_y_xxz_1[i] * pc_y[i];

        ta2_xz_yy_xyy_0[i] = ta2_xz_0_xyy_0[i] * fe_0 - ta2_xz_0_xyy_1[i] * fe_0 + 2.0 * ta2_xz_y_xy_0[i] * fe_0 - 2.0 * ta2_xz_y_xy_1[i] * fe_0 +
                             ta2_xz_y_xyy_0[i] * pa_y[i] - ta2_xz_y_xyy_1[i] * pc_y[i];

        ta2_xz_yy_xyz_0[i] = ta2_xz_0_xyz_0[i] * fe_0 - ta2_xz_0_xyz_1[i] * fe_0 + ta2_xz_y_xz_0[i] * fe_0 - ta2_xz_y_xz_1[i] * fe_0 +
                             ta2_xz_y_xyz_0[i] * pa_y[i] - ta2_xz_y_xyz_1[i] * pc_y[i];

        ta2_xz_yy_xzz_0[i] = ta2_xz_0_xzz_0[i] * fe_0 - ta2_xz_0_xzz_1[i] * fe_0 + ta2_xz_y_xzz_0[i] * pa_y[i] - ta2_xz_y_xzz_1[i] * pc_y[i];

        ta2_xz_yy_yyy_0[i] = ta2_xz_0_yyy_0[i] * fe_0 - ta2_xz_0_yyy_1[i] * fe_0 + 3.0 * ta2_xz_y_yy_0[i] * fe_0 - 3.0 * ta2_xz_y_yy_1[i] * fe_0 +
                             ta2_xz_y_yyy_0[i] * pa_y[i] - ta2_xz_y_yyy_1[i] * pc_y[i];

        ta2_xz_yy_yyz_0[i] = ta2_xz_0_yyz_0[i] * fe_0 - ta2_xz_0_yyz_1[i] * fe_0 + 2.0 * ta2_xz_y_yz_0[i] * fe_0 - 2.0 * ta2_xz_y_yz_1[i] * fe_0 +
                             ta2_xz_y_yyz_0[i] * pa_y[i] - ta2_xz_y_yyz_1[i] * pc_y[i];

        ta2_xz_yy_yzz_0[i] = ta2_xz_0_yzz_0[i] * fe_0 - ta2_xz_0_yzz_1[i] * fe_0 + ta2_xz_y_zz_0[i] * fe_0 - ta2_xz_y_zz_1[i] * fe_0 +
                             ta2_xz_y_yzz_0[i] * pa_y[i] - ta2_xz_y_yzz_1[i] * pc_y[i];

        ta2_xz_yy_zzz_0[i] = ta2_xz_0_zzz_0[i] * fe_0 - ta2_xz_0_zzz_1[i] * fe_0 + ta2_xz_y_zzz_0[i] * pa_y[i] - ta2_xz_y_zzz_1[i] * pc_y[i];
    }

    // Set up 160-170 components of targeted buffer : DF

    auto ta2_xz_yz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 160);

    auto ta2_xz_yz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 161);

    auto ta2_xz_yz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 162);

    auto ta2_xz_yz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 163);

    auto ta2_xz_yz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 164);

    auto ta2_xz_yz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 165);

    auto ta2_xz_yz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 166);

    auto ta2_xz_yz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 167);

    auto ta2_xz_yz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 168);

    auto ta2_xz_yz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 169);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_x_y_xxy_1,   \
                             ta1_x_y_xyy_1,   \
                             ta1_x_y_yyy_1,   \
                             ta2_xz_y_xxy_0,  \
                             ta2_xz_y_xxy_1,  \
                             ta2_xz_y_xyy_0,  \
                             ta2_xz_y_xyy_1,  \
                             ta2_xz_y_yyy_0,  \
                             ta2_xz_y_yyy_1,  \
                             ta2_xz_yz_xxx_0, \
                             ta2_xz_yz_xxy_0, \
                             ta2_xz_yz_xxz_0, \
                             ta2_xz_yz_xyy_0, \
                             ta2_xz_yz_xyz_0, \
                             ta2_xz_yz_xzz_0, \
                             ta2_xz_yz_yyy_0, \
                             ta2_xz_yz_yyz_0, \
                             ta2_xz_yz_yzz_0, \
                             ta2_xz_yz_zzz_0, \
                             ta2_xz_z_xxx_0,  \
                             ta2_xz_z_xxx_1,  \
                             ta2_xz_z_xxz_0,  \
                             ta2_xz_z_xxz_1,  \
                             ta2_xz_z_xyz_0,  \
                             ta2_xz_z_xyz_1,  \
                             ta2_xz_z_xz_0,   \
                             ta2_xz_z_xz_1,   \
                             ta2_xz_z_xzz_0,  \
                             ta2_xz_z_xzz_1,  \
                             ta2_xz_z_yyz_0,  \
                             ta2_xz_z_yyz_1,  \
                             ta2_xz_z_yz_0,   \
                             ta2_xz_z_yz_1,   \
                             ta2_xz_z_yzz_0,  \
                             ta2_xz_z_yzz_1,  \
                             ta2_xz_z_zz_0,   \
                             ta2_xz_z_zz_1,   \
                             ta2_xz_z_zzz_0,  \
                             ta2_xz_z_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yz_xxx_0[i] = ta2_xz_z_xxx_0[i] * pa_y[i] - ta2_xz_z_xxx_1[i] * pc_y[i];

        ta2_xz_yz_xxy_0[i] = ta1_x_y_xxy_1[i] + ta2_xz_y_xxy_0[i] * pa_z[i] - ta2_xz_y_xxy_1[i] * pc_z[i];

        ta2_xz_yz_xxz_0[i] = ta2_xz_z_xxz_0[i] * pa_y[i] - ta2_xz_z_xxz_1[i] * pc_y[i];

        ta2_xz_yz_xyy_0[i] = ta1_x_y_xyy_1[i] + ta2_xz_y_xyy_0[i] * pa_z[i] - ta2_xz_y_xyy_1[i] * pc_z[i];

        ta2_xz_yz_xyz_0[i] = ta2_xz_z_xz_0[i] * fe_0 - ta2_xz_z_xz_1[i] * fe_0 + ta2_xz_z_xyz_0[i] * pa_y[i] - ta2_xz_z_xyz_1[i] * pc_y[i];

        ta2_xz_yz_xzz_0[i] = ta2_xz_z_xzz_0[i] * pa_y[i] - ta2_xz_z_xzz_1[i] * pc_y[i];

        ta2_xz_yz_yyy_0[i] = ta1_x_y_yyy_1[i] + ta2_xz_y_yyy_0[i] * pa_z[i] - ta2_xz_y_yyy_1[i] * pc_z[i];

        ta2_xz_yz_yyz_0[i] =
            2.0 * ta2_xz_z_yz_0[i] * fe_0 - 2.0 * ta2_xz_z_yz_1[i] * fe_0 + ta2_xz_z_yyz_0[i] * pa_y[i] - ta2_xz_z_yyz_1[i] * pc_y[i];

        ta2_xz_yz_yzz_0[i] = ta2_xz_z_zz_0[i] * fe_0 - ta2_xz_z_zz_1[i] * fe_0 + ta2_xz_z_yzz_0[i] * pa_y[i] - ta2_xz_z_yzz_1[i] * pc_y[i];

        ta2_xz_yz_zzz_0[i] = ta2_xz_z_zzz_0[i] * pa_y[i] - ta2_xz_z_zzz_1[i] * pc_y[i];
    }

    // Set up 170-180 components of targeted buffer : DF

    auto ta2_xz_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 170);

    auto ta2_xz_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 171);

    auto ta2_xz_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 172);

    auto ta2_xz_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 173);

    auto ta2_xz_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 174);

    auto ta2_xz_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 175);

    auto ta2_xz_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 176);

    auto ta2_xz_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 177);

    auto ta2_xz_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 178);

    auto ta2_xz_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 179);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_x_z_xxx_1,   \
                             ta1_x_z_xxy_1,   \
                             ta1_x_z_xxz_1,   \
                             ta1_x_z_xyy_1,   \
                             ta1_x_z_xyz_1,   \
                             ta1_x_z_xzz_1,   \
                             ta1_x_z_yyy_1,   \
                             ta1_x_z_yyz_1,   \
                             ta1_x_z_yzz_1,   \
                             ta1_x_z_zzz_1,   \
                             ta2_xz_0_xxx_0,  \
                             ta2_xz_0_xxx_1,  \
                             ta2_xz_0_xxy_0,  \
                             ta2_xz_0_xxy_1,  \
                             ta2_xz_0_xxz_0,  \
                             ta2_xz_0_xxz_1,  \
                             ta2_xz_0_xyy_0,  \
                             ta2_xz_0_xyy_1,  \
                             ta2_xz_0_xyz_0,  \
                             ta2_xz_0_xyz_1,  \
                             ta2_xz_0_xzz_0,  \
                             ta2_xz_0_xzz_1,  \
                             ta2_xz_0_yyy_0,  \
                             ta2_xz_0_yyy_1,  \
                             ta2_xz_0_yyz_0,  \
                             ta2_xz_0_yyz_1,  \
                             ta2_xz_0_yzz_0,  \
                             ta2_xz_0_yzz_1,  \
                             ta2_xz_0_zzz_0,  \
                             ta2_xz_0_zzz_1,  \
                             ta2_xz_z_xx_0,   \
                             ta2_xz_z_xx_1,   \
                             ta2_xz_z_xxx_0,  \
                             ta2_xz_z_xxx_1,  \
                             ta2_xz_z_xxy_0,  \
                             ta2_xz_z_xxy_1,  \
                             ta2_xz_z_xxz_0,  \
                             ta2_xz_z_xxz_1,  \
                             ta2_xz_z_xy_0,   \
                             ta2_xz_z_xy_1,   \
                             ta2_xz_z_xyy_0,  \
                             ta2_xz_z_xyy_1,  \
                             ta2_xz_z_xyz_0,  \
                             ta2_xz_z_xyz_1,  \
                             ta2_xz_z_xz_0,   \
                             ta2_xz_z_xz_1,   \
                             ta2_xz_z_xzz_0,  \
                             ta2_xz_z_xzz_1,  \
                             ta2_xz_z_yy_0,   \
                             ta2_xz_z_yy_1,   \
                             ta2_xz_z_yyy_0,  \
                             ta2_xz_z_yyy_1,  \
                             ta2_xz_z_yyz_0,  \
                             ta2_xz_z_yyz_1,  \
                             ta2_xz_z_yz_0,   \
                             ta2_xz_z_yz_1,   \
                             ta2_xz_z_yzz_0,  \
                             ta2_xz_z_yzz_1,  \
                             ta2_xz_z_zz_0,   \
                             ta2_xz_z_zz_1,   \
                             ta2_xz_z_zzz_0,  \
                             ta2_xz_z_zzz_1,  \
                             ta2_xz_zz_xxx_0, \
                             ta2_xz_zz_xxy_0, \
                             ta2_xz_zz_xxz_0, \
                             ta2_xz_zz_xyy_0, \
                             ta2_xz_zz_xyz_0, \
                             ta2_xz_zz_xzz_0, \
                             ta2_xz_zz_yyy_0, \
                             ta2_xz_zz_yyz_0, \
                             ta2_xz_zz_yzz_0, \
                             ta2_xz_zz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_zz_xxx_0[i] =
            ta2_xz_0_xxx_0[i] * fe_0 - ta2_xz_0_xxx_1[i] * fe_0 + ta1_x_z_xxx_1[i] + ta2_xz_z_xxx_0[i] * pa_z[i] - ta2_xz_z_xxx_1[i] * pc_z[i];

        ta2_xz_zz_xxy_0[i] =
            ta2_xz_0_xxy_0[i] * fe_0 - ta2_xz_0_xxy_1[i] * fe_0 + ta1_x_z_xxy_1[i] + ta2_xz_z_xxy_0[i] * pa_z[i] - ta2_xz_z_xxy_1[i] * pc_z[i];

        ta2_xz_zz_xxz_0[i] = ta2_xz_0_xxz_0[i] * fe_0 - ta2_xz_0_xxz_1[i] * fe_0 + ta2_xz_z_xx_0[i] * fe_0 - ta2_xz_z_xx_1[i] * fe_0 +
                             ta1_x_z_xxz_1[i] + ta2_xz_z_xxz_0[i] * pa_z[i] - ta2_xz_z_xxz_1[i] * pc_z[i];

        ta2_xz_zz_xyy_0[i] =
            ta2_xz_0_xyy_0[i] * fe_0 - ta2_xz_0_xyy_1[i] * fe_0 + ta1_x_z_xyy_1[i] + ta2_xz_z_xyy_0[i] * pa_z[i] - ta2_xz_z_xyy_1[i] * pc_z[i];

        ta2_xz_zz_xyz_0[i] = ta2_xz_0_xyz_0[i] * fe_0 - ta2_xz_0_xyz_1[i] * fe_0 + ta2_xz_z_xy_0[i] * fe_0 - ta2_xz_z_xy_1[i] * fe_0 +
                             ta1_x_z_xyz_1[i] + ta2_xz_z_xyz_0[i] * pa_z[i] - ta2_xz_z_xyz_1[i] * pc_z[i];

        ta2_xz_zz_xzz_0[i] = ta2_xz_0_xzz_0[i] * fe_0 - ta2_xz_0_xzz_1[i] * fe_0 + 2.0 * ta2_xz_z_xz_0[i] * fe_0 - 2.0 * ta2_xz_z_xz_1[i] * fe_0 +
                             ta1_x_z_xzz_1[i] + ta2_xz_z_xzz_0[i] * pa_z[i] - ta2_xz_z_xzz_1[i] * pc_z[i];

        ta2_xz_zz_yyy_0[i] =
            ta2_xz_0_yyy_0[i] * fe_0 - ta2_xz_0_yyy_1[i] * fe_0 + ta1_x_z_yyy_1[i] + ta2_xz_z_yyy_0[i] * pa_z[i] - ta2_xz_z_yyy_1[i] * pc_z[i];

        ta2_xz_zz_yyz_0[i] = ta2_xz_0_yyz_0[i] * fe_0 - ta2_xz_0_yyz_1[i] * fe_0 + ta2_xz_z_yy_0[i] * fe_0 - ta2_xz_z_yy_1[i] * fe_0 +
                             ta1_x_z_yyz_1[i] + ta2_xz_z_yyz_0[i] * pa_z[i] - ta2_xz_z_yyz_1[i] * pc_z[i];

        ta2_xz_zz_yzz_0[i] = ta2_xz_0_yzz_0[i] * fe_0 - ta2_xz_0_yzz_1[i] * fe_0 + 2.0 * ta2_xz_z_yz_0[i] * fe_0 - 2.0 * ta2_xz_z_yz_1[i] * fe_0 +
                             ta1_x_z_yzz_1[i] + ta2_xz_z_yzz_0[i] * pa_z[i] - ta2_xz_z_yzz_1[i] * pc_z[i];

        ta2_xz_zz_zzz_0[i] = ta2_xz_0_zzz_0[i] * fe_0 - ta2_xz_0_zzz_1[i] * fe_0 + 3.0 * ta2_xz_z_zz_0[i] * fe_0 - 3.0 * ta2_xz_z_zz_1[i] * fe_0 +
                             ta1_x_z_zzz_1[i] + ta2_xz_z_zzz_0[i] * pa_z[i] - ta2_xz_z_zzz_1[i] * pc_z[i];
    }

    // Set up 180-190 components of targeted buffer : DF

    auto ta2_yy_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 180);

    auto ta2_yy_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 181);

    auto ta2_yy_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 182);

    auto ta2_yy_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 183);

    auto ta2_yy_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 184);

    auto ta2_yy_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 185);

    auto ta2_yy_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 186);

    auto ta2_yy_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 187);

    auto ta2_yy_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 188);

    auto ta2_yy_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 189);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta2_yy_0_xxx_0,  \
                             ta2_yy_0_xxx_1,  \
                             ta2_yy_0_xxy_0,  \
                             ta2_yy_0_xxy_1,  \
                             ta2_yy_0_xxz_0,  \
                             ta2_yy_0_xxz_1,  \
                             ta2_yy_0_xyy_0,  \
                             ta2_yy_0_xyy_1,  \
                             ta2_yy_0_xyz_0,  \
                             ta2_yy_0_xyz_1,  \
                             ta2_yy_0_xzz_0,  \
                             ta2_yy_0_xzz_1,  \
                             ta2_yy_0_yyy_0,  \
                             ta2_yy_0_yyy_1,  \
                             ta2_yy_0_yyz_0,  \
                             ta2_yy_0_yyz_1,  \
                             ta2_yy_0_yzz_0,  \
                             ta2_yy_0_yzz_1,  \
                             ta2_yy_0_zzz_0,  \
                             ta2_yy_0_zzz_1,  \
                             ta2_yy_x_xx_0,   \
                             ta2_yy_x_xx_1,   \
                             ta2_yy_x_xxx_0,  \
                             ta2_yy_x_xxx_1,  \
                             ta2_yy_x_xxy_0,  \
                             ta2_yy_x_xxy_1,  \
                             ta2_yy_x_xxz_0,  \
                             ta2_yy_x_xxz_1,  \
                             ta2_yy_x_xy_0,   \
                             ta2_yy_x_xy_1,   \
                             ta2_yy_x_xyy_0,  \
                             ta2_yy_x_xyy_1,  \
                             ta2_yy_x_xyz_0,  \
                             ta2_yy_x_xyz_1,  \
                             ta2_yy_x_xz_0,   \
                             ta2_yy_x_xz_1,   \
                             ta2_yy_x_xzz_0,  \
                             ta2_yy_x_xzz_1,  \
                             ta2_yy_x_yy_0,   \
                             ta2_yy_x_yy_1,   \
                             ta2_yy_x_yyy_0,  \
                             ta2_yy_x_yyy_1,  \
                             ta2_yy_x_yyz_0,  \
                             ta2_yy_x_yyz_1,  \
                             ta2_yy_x_yz_0,   \
                             ta2_yy_x_yz_1,   \
                             ta2_yy_x_yzz_0,  \
                             ta2_yy_x_yzz_1,  \
                             ta2_yy_x_zz_0,   \
                             ta2_yy_x_zz_1,   \
                             ta2_yy_x_zzz_0,  \
                             ta2_yy_x_zzz_1,  \
                             ta2_yy_xx_xxx_0, \
                             ta2_yy_xx_xxy_0, \
                             ta2_yy_xx_xxz_0, \
                             ta2_yy_xx_xyy_0, \
                             ta2_yy_xx_xyz_0, \
                             ta2_yy_xx_xzz_0, \
                             ta2_yy_xx_yyy_0, \
                             ta2_yy_xx_yyz_0, \
                             ta2_yy_xx_yzz_0, \
                             ta2_yy_xx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xx_xxx_0[i] = ta2_yy_0_xxx_0[i] * fe_0 - ta2_yy_0_xxx_1[i] * fe_0 + 3.0 * ta2_yy_x_xx_0[i] * fe_0 - 3.0 * ta2_yy_x_xx_1[i] * fe_0 +
                             ta2_yy_x_xxx_0[i] * pa_x[i] - ta2_yy_x_xxx_1[i] * pc_x[i];

        ta2_yy_xx_xxy_0[i] = ta2_yy_0_xxy_0[i] * fe_0 - ta2_yy_0_xxy_1[i] * fe_0 + 2.0 * ta2_yy_x_xy_0[i] * fe_0 - 2.0 * ta2_yy_x_xy_1[i] * fe_0 +
                             ta2_yy_x_xxy_0[i] * pa_x[i] - ta2_yy_x_xxy_1[i] * pc_x[i];

        ta2_yy_xx_xxz_0[i] = ta2_yy_0_xxz_0[i] * fe_0 - ta2_yy_0_xxz_1[i] * fe_0 + 2.0 * ta2_yy_x_xz_0[i] * fe_0 - 2.0 * ta2_yy_x_xz_1[i] * fe_0 +
                             ta2_yy_x_xxz_0[i] * pa_x[i] - ta2_yy_x_xxz_1[i] * pc_x[i];

        ta2_yy_xx_xyy_0[i] = ta2_yy_0_xyy_0[i] * fe_0 - ta2_yy_0_xyy_1[i] * fe_0 + ta2_yy_x_yy_0[i] * fe_0 - ta2_yy_x_yy_1[i] * fe_0 +
                             ta2_yy_x_xyy_0[i] * pa_x[i] - ta2_yy_x_xyy_1[i] * pc_x[i];

        ta2_yy_xx_xyz_0[i] = ta2_yy_0_xyz_0[i] * fe_0 - ta2_yy_0_xyz_1[i] * fe_0 + ta2_yy_x_yz_0[i] * fe_0 - ta2_yy_x_yz_1[i] * fe_0 +
                             ta2_yy_x_xyz_0[i] * pa_x[i] - ta2_yy_x_xyz_1[i] * pc_x[i];

        ta2_yy_xx_xzz_0[i] = ta2_yy_0_xzz_0[i] * fe_0 - ta2_yy_0_xzz_1[i] * fe_0 + ta2_yy_x_zz_0[i] * fe_0 - ta2_yy_x_zz_1[i] * fe_0 +
                             ta2_yy_x_xzz_0[i] * pa_x[i] - ta2_yy_x_xzz_1[i] * pc_x[i];

        ta2_yy_xx_yyy_0[i] = ta2_yy_0_yyy_0[i] * fe_0 - ta2_yy_0_yyy_1[i] * fe_0 + ta2_yy_x_yyy_0[i] * pa_x[i] - ta2_yy_x_yyy_1[i] * pc_x[i];

        ta2_yy_xx_yyz_0[i] = ta2_yy_0_yyz_0[i] * fe_0 - ta2_yy_0_yyz_1[i] * fe_0 + ta2_yy_x_yyz_0[i] * pa_x[i] - ta2_yy_x_yyz_1[i] * pc_x[i];

        ta2_yy_xx_yzz_0[i] = ta2_yy_0_yzz_0[i] * fe_0 - ta2_yy_0_yzz_1[i] * fe_0 + ta2_yy_x_yzz_0[i] * pa_x[i] - ta2_yy_x_yzz_1[i] * pc_x[i];

        ta2_yy_xx_zzz_0[i] = ta2_yy_0_zzz_0[i] * fe_0 - ta2_yy_0_zzz_1[i] * fe_0 + ta2_yy_x_zzz_0[i] * pa_x[i] - ta2_yy_x_zzz_1[i] * pc_x[i];
    }

    // Set up 190-200 components of targeted buffer : DF

    auto ta2_yy_xy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 190);

    auto ta2_yy_xy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 191);

    auto ta2_yy_xy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 192);

    auto ta2_yy_xy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 193);

    auto ta2_yy_xy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 194);

    auto ta2_yy_xy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 195);

    auto ta2_yy_xy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 196);

    auto ta2_yy_xy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 197);

    auto ta2_yy_xy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 198);

    auto ta2_yy_xy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 199);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_y_x_xxx_1,   \
                             ta1_y_x_xxz_1,   \
                             ta1_y_x_xzz_1,   \
                             ta2_yy_x_xxx_0,  \
                             ta2_yy_x_xxx_1,  \
                             ta2_yy_x_xxz_0,  \
                             ta2_yy_x_xxz_1,  \
                             ta2_yy_x_xzz_0,  \
                             ta2_yy_x_xzz_1,  \
                             ta2_yy_xy_xxx_0, \
                             ta2_yy_xy_xxy_0, \
                             ta2_yy_xy_xxz_0, \
                             ta2_yy_xy_xyy_0, \
                             ta2_yy_xy_xyz_0, \
                             ta2_yy_xy_xzz_0, \
                             ta2_yy_xy_yyy_0, \
                             ta2_yy_xy_yyz_0, \
                             ta2_yy_xy_yzz_0, \
                             ta2_yy_xy_zzz_0, \
                             ta2_yy_y_xxy_0,  \
                             ta2_yy_y_xxy_1,  \
                             ta2_yy_y_xy_0,   \
                             ta2_yy_y_xy_1,   \
                             ta2_yy_y_xyy_0,  \
                             ta2_yy_y_xyy_1,  \
                             ta2_yy_y_xyz_0,  \
                             ta2_yy_y_xyz_1,  \
                             ta2_yy_y_yy_0,   \
                             ta2_yy_y_yy_1,   \
                             ta2_yy_y_yyy_0,  \
                             ta2_yy_y_yyy_1,  \
                             ta2_yy_y_yyz_0,  \
                             ta2_yy_y_yyz_1,  \
                             ta2_yy_y_yz_0,   \
                             ta2_yy_y_yz_1,   \
                             ta2_yy_y_yzz_0,  \
                             ta2_yy_y_yzz_1,  \
                             ta2_yy_y_zzz_0,  \
                             ta2_yy_y_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xy_xxx_0[i] = 2.0 * ta1_y_x_xxx_1[i] + ta2_yy_x_xxx_0[i] * pa_y[i] - ta2_yy_x_xxx_1[i] * pc_y[i];

        ta2_yy_xy_xxy_0[i] =
            2.0 * ta2_yy_y_xy_0[i] * fe_0 - 2.0 * ta2_yy_y_xy_1[i] * fe_0 + ta2_yy_y_xxy_0[i] * pa_x[i] - ta2_yy_y_xxy_1[i] * pc_x[i];

        ta2_yy_xy_xxz_0[i] = 2.0 * ta1_y_x_xxz_1[i] + ta2_yy_x_xxz_0[i] * pa_y[i] - ta2_yy_x_xxz_1[i] * pc_y[i];

        ta2_yy_xy_xyy_0[i] = ta2_yy_y_yy_0[i] * fe_0 - ta2_yy_y_yy_1[i] * fe_0 + ta2_yy_y_xyy_0[i] * pa_x[i] - ta2_yy_y_xyy_1[i] * pc_x[i];

        ta2_yy_xy_xyz_0[i] = ta2_yy_y_yz_0[i] * fe_0 - ta2_yy_y_yz_1[i] * fe_0 + ta2_yy_y_xyz_0[i] * pa_x[i] - ta2_yy_y_xyz_1[i] * pc_x[i];

        ta2_yy_xy_xzz_0[i] = 2.0 * ta1_y_x_xzz_1[i] + ta2_yy_x_xzz_0[i] * pa_y[i] - ta2_yy_x_xzz_1[i] * pc_y[i];

        ta2_yy_xy_yyy_0[i] = ta2_yy_y_yyy_0[i] * pa_x[i] - ta2_yy_y_yyy_1[i] * pc_x[i];

        ta2_yy_xy_yyz_0[i] = ta2_yy_y_yyz_0[i] * pa_x[i] - ta2_yy_y_yyz_1[i] * pc_x[i];

        ta2_yy_xy_yzz_0[i] = ta2_yy_y_yzz_0[i] * pa_x[i] - ta2_yy_y_yzz_1[i] * pc_x[i];

        ta2_yy_xy_zzz_0[i] = ta2_yy_y_zzz_0[i] * pa_x[i] - ta2_yy_y_zzz_1[i] * pc_x[i];
    }

    // Set up 200-210 components of targeted buffer : DF

    auto ta2_yy_xz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 200);

    auto ta2_yy_xz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 201);

    auto ta2_yy_xz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 202);

    auto ta2_yy_xz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 203);

    auto ta2_yy_xz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 204);

    auto ta2_yy_xz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 205);

    auto ta2_yy_xz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 206);

    auto ta2_yy_xz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 207);

    auto ta2_yy_xz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 208);

    auto ta2_yy_xz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 209);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta2_yy_x_xxx_0,  \
                             ta2_yy_x_xxx_1,  \
                             ta2_yy_x_xxy_0,  \
                             ta2_yy_x_xxy_1,  \
                             ta2_yy_x_xyy_0,  \
                             ta2_yy_x_xyy_1,  \
                             ta2_yy_xz_xxx_0, \
                             ta2_yy_xz_xxy_0, \
                             ta2_yy_xz_xxz_0, \
                             ta2_yy_xz_xyy_0, \
                             ta2_yy_xz_xyz_0, \
                             ta2_yy_xz_xzz_0, \
                             ta2_yy_xz_yyy_0, \
                             ta2_yy_xz_yyz_0, \
                             ta2_yy_xz_yzz_0, \
                             ta2_yy_xz_zzz_0, \
                             ta2_yy_z_xxz_0,  \
                             ta2_yy_z_xxz_1,  \
                             ta2_yy_z_xyz_0,  \
                             ta2_yy_z_xyz_1,  \
                             ta2_yy_z_xz_0,   \
                             ta2_yy_z_xz_1,   \
                             ta2_yy_z_xzz_0,  \
                             ta2_yy_z_xzz_1,  \
                             ta2_yy_z_yyy_0,  \
                             ta2_yy_z_yyy_1,  \
                             ta2_yy_z_yyz_0,  \
                             ta2_yy_z_yyz_1,  \
                             ta2_yy_z_yz_0,   \
                             ta2_yy_z_yz_1,   \
                             ta2_yy_z_yzz_0,  \
                             ta2_yy_z_yzz_1,  \
                             ta2_yy_z_zz_0,   \
                             ta2_yy_z_zz_1,   \
                             ta2_yy_z_zzz_0,  \
                             ta2_yy_z_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xz_xxx_0[i] = ta2_yy_x_xxx_0[i] * pa_z[i] - ta2_yy_x_xxx_1[i] * pc_z[i];

        ta2_yy_xz_xxy_0[i] = ta2_yy_x_xxy_0[i] * pa_z[i] - ta2_yy_x_xxy_1[i] * pc_z[i];

        ta2_yy_xz_xxz_0[i] =
            2.0 * ta2_yy_z_xz_0[i] * fe_0 - 2.0 * ta2_yy_z_xz_1[i] * fe_0 + ta2_yy_z_xxz_0[i] * pa_x[i] - ta2_yy_z_xxz_1[i] * pc_x[i];

        ta2_yy_xz_xyy_0[i] = ta2_yy_x_xyy_0[i] * pa_z[i] - ta2_yy_x_xyy_1[i] * pc_z[i];

        ta2_yy_xz_xyz_0[i] = ta2_yy_z_yz_0[i] * fe_0 - ta2_yy_z_yz_1[i] * fe_0 + ta2_yy_z_xyz_0[i] * pa_x[i] - ta2_yy_z_xyz_1[i] * pc_x[i];

        ta2_yy_xz_xzz_0[i] = ta2_yy_z_zz_0[i] * fe_0 - ta2_yy_z_zz_1[i] * fe_0 + ta2_yy_z_xzz_0[i] * pa_x[i] - ta2_yy_z_xzz_1[i] * pc_x[i];

        ta2_yy_xz_yyy_0[i] = ta2_yy_z_yyy_0[i] * pa_x[i] - ta2_yy_z_yyy_1[i] * pc_x[i];

        ta2_yy_xz_yyz_0[i] = ta2_yy_z_yyz_0[i] * pa_x[i] - ta2_yy_z_yyz_1[i] * pc_x[i];

        ta2_yy_xz_yzz_0[i] = ta2_yy_z_yzz_0[i] * pa_x[i] - ta2_yy_z_yzz_1[i] * pc_x[i];

        ta2_yy_xz_zzz_0[i] = ta2_yy_z_zzz_0[i] * pa_x[i] - ta2_yy_z_zzz_1[i] * pc_x[i];
    }

    // Set up 210-220 components of targeted buffer : DF

    auto ta2_yy_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 210);

    auto ta2_yy_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 211);

    auto ta2_yy_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 212);

    auto ta2_yy_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 213);

    auto ta2_yy_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 214);

    auto ta2_yy_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 215);

    auto ta2_yy_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 216);

    auto ta2_yy_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 217);

    auto ta2_yy_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 218);

    auto ta2_yy_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 219);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_y_y_xxx_1,   \
                             ta1_y_y_xxy_1,   \
                             ta1_y_y_xxz_1,   \
                             ta1_y_y_xyy_1,   \
                             ta1_y_y_xyz_1,   \
                             ta1_y_y_xzz_1,   \
                             ta1_y_y_yyy_1,   \
                             ta1_y_y_yyz_1,   \
                             ta1_y_y_yzz_1,   \
                             ta1_y_y_zzz_1,   \
                             ta2_yy_0_xxx_0,  \
                             ta2_yy_0_xxx_1,  \
                             ta2_yy_0_xxy_0,  \
                             ta2_yy_0_xxy_1,  \
                             ta2_yy_0_xxz_0,  \
                             ta2_yy_0_xxz_1,  \
                             ta2_yy_0_xyy_0,  \
                             ta2_yy_0_xyy_1,  \
                             ta2_yy_0_xyz_0,  \
                             ta2_yy_0_xyz_1,  \
                             ta2_yy_0_xzz_0,  \
                             ta2_yy_0_xzz_1,  \
                             ta2_yy_0_yyy_0,  \
                             ta2_yy_0_yyy_1,  \
                             ta2_yy_0_yyz_0,  \
                             ta2_yy_0_yyz_1,  \
                             ta2_yy_0_yzz_0,  \
                             ta2_yy_0_yzz_1,  \
                             ta2_yy_0_zzz_0,  \
                             ta2_yy_0_zzz_1,  \
                             ta2_yy_y_xx_0,   \
                             ta2_yy_y_xx_1,   \
                             ta2_yy_y_xxx_0,  \
                             ta2_yy_y_xxx_1,  \
                             ta2_yy_y_xxy_0,  \
                             ta2_yy_y_xxy_1,  \
                             ta2_yy_y_xxz_0,  \
                             ta2_yy_y_xxz_1,  \
                             ta2_yy_y_xy_0,   \
                             ta2_yy_y_xy_1,   \
                             ta2_yy_y_xyy_0,  \
                             ta2_yy_y_xyy_1,  \
                             ta2_yy_y_xyz_0,  \
                             ta2_yy_y_xyz_1,  \
                             ta2_yy_y_xz_0,   \
                             ta2_yy_y_xz_1,   \
                             ta2_yy_y_xzz_0,  \
                             ta2_yy_y_xzz_1,  \
                             ta2_yy_y_yy_0,   \
                             ta2_yy_y_yy_1,   \
                             ta2_yy_y_yyy_0,  \
                             ta2_yy_y_yyy_1,  \
                             ta2_yy_y_yyz_0,  \
                             ta2_yy_y_yyz_1,  \
                             ta2_yy_y_yz_0,   \
                             ta2_yy_y_yz_1,   \
                             ta2_yy_y_yzz_0,  \
                             ta2_yy_y_yzz_1,  \
                             ta2_yy_y_zz_0,   \
                             ta2_yy_y_zz_1,   \
                             ta2_yy_y_zzz_0,  \
                             ta2_yy_y_zzz_1,  \
                             ta2_yy_yy_xxx_0, \
                             ta2_yy_yy_xxy_0, \
                             ta2_yy_yy_xxz_0, \
                             ta2_yy_yy_xyy_0, \
                             ta2_yy_yy_xyz_0, \
                             ta2_yy_yy_xzz_0, \
                             ta2_yy_yy_yyy_0, \
                             ta2_yy_yy_yyz_0, \
                             ta2_yy_yy_yzz_0, \
                             ta2_yy_yy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yy_xxx_0[i] =
            ta2_yy_0_xxx_0[i] * fe_0 - ta2_yy_0_xxx_1[i] * fe_0 + 2.0 * ta1_y_y_xxx_1[i] + ta2_yy_y_xxx_0[i] * pa_y[i] - ta2_yy_y_xxx_1[i] * pc_y[i];

        ta2_yy_yy_xxy_0[i] = ta2_yy_0_xxy_0[i] * fe_0 - ta2_yy_0_xxy_1[i] * fe_0 + ta2_yy_y_xx_0[i] * fe_0 - ta2_yy_y_xx_1[i] * fe_0 +
                             2.0 * ta1_y_y_xxy_1[i] + ta2_yy_y_xxy_0[i] * pa_y[i] - ta2_yy_y_xxy_1[i] * pc_y[i];

        ta2_yy_yy_xxz_0[i] =
            ta2_yy_0_xxz_0[i] * fe_0 - ta2_yy_0_xxz_1[i] * fe_0 + 2.0 * ta1_y_y_xxz_1[i] + ta2_yy_y_xxz_0[i] * pa_y[i] - ta2_yy_y_xxz_1[i] * pc_y[i];

        ta2_yy_yy_xyy_0[i] = ta2_yy_0_xyy_0[i] * fe_0 - ta2_yy_0_xyy_1[i] * fe_0 + 2.0 * ta2_yy_y_xy_0[i] * fe_0 - 2.0 * ta2_yy_y_xy_1[i] * fe_0 +
                             2.0 * ta1_y_y_xyy_1[i] + ta2_yy_y_xyy_0[i] * pa_y[i] - ta2_yy_y_xyy_1[i] * pc_y[i];

        ta2_yy_yy_xyz_0[i] = ta2_yy_0_xyz_0[i] * fe_0 - ta2_yy_0_xyz_1[i] * fe_0 + ta2_yy_y_xz_0[i] * fe_0 - ta2_yy_y_xz_1[i] * fe_0 +
                             2.0 * ta1_y_y_xyz_1[i] + ta2_yy_y_xyz_0[i] * pa_y[i] - ta2_yy_y_xyz_1[i] * pc_y[i];

        ta2_yy_yy_xzz_0[i] =
            ta2_yy_0_xzz_0[i] * fe_0 - ta2_yy_0_xzz_1[i] * fe_0 + 2.0 * ta1_y_y_xzz_1[i] + ta2_yy_y_xzz_0[i] * pa_y[i] - ta2_yy_y_xzz_1[i] * pc_y[i];

        ta2_yy_yy_yyy_0[i] = ta2_yy_0_yyy_0[i] * fe_0 - ta2_yy_0_yyy_1[i] * fe_0 + 3.0 * ta2_yy_y_yy_0[i] * fe_0 - 3.0 * ta2_yy_y_yy_1[i] * fe_0 +
                             2.0 * ta1_y_y_yyy_1[i] + ta2_yy_y_yyy_0[i] * pa_y[i] - ta2_yy_y_yyy_1[i] * pc_y[i];

        ta2_yy_yy_yyz_0[i] = ta2_yy_0_yyz_0[i] * fe_0 - ta2_yy_0_yyz_1[i] * fe_0 + 2.0 * ta2_yy_y_yz_0[i] * fe_0 - 2.0 * ta2_yy_y_yz_1[i] * fe_0 +
                             2.0 * ta1_y_y_yyz_1[i] + ta2_yy_y_yyz_0[i] * pa_y[i] - ta2_yy_y_yyz_1[i] * pc_y[i];

        ta2_yy_yy_yzz_0[i] = ta2_yy_0_yzz_0[i] * fe_0 - ta2_yy_0_yzz_1[i] * fe_0 + ta2_yy_y_zz_0[i] * fe_0 - ta2_yy_y_zz_1[i] * fe_0 +
                             2.0 * ta1_y_y_yzz_1[i] + ta2_yy_y_yzz_0[i] * pa_y[i] - ta2_yy_y_yzz_1[i] * pc_y[i];

        ta2_yy_yy_zzz_0[i] =
            ta2_yy_0_zzz_0[i] * fe_0 - ta2_yy_0_zzz_1[i] * fe_0 + 2.0 * ta1_y_y_zzz_1[i] + ta2_yy_y_zzz_0[i] * pa_y[i] - ta2_yy_y_zzz_1[i] * pc_y[i];
    }

    // Set up 220-230 components of targeted buffer : DF

    auto ta2_yy_yz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 220);

    auto ta2_yy_yz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 221);

    auto ta2_yy_yz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 222);

    auto ta2_yy_yz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 223);

    auto ta2_yy_yz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 224);

    auto ta2_yy_yz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 225);

    auto ta2_yy_yz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 226);

    auto ta2_yy_yz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 227);

    auto ta2_yy_yz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 228);

    auto ta2_yy_yz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 229);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_z_xxz_1,   \
                             ta1_y_z_xzz_1,   \
                             ta1_y_z_zzz_1,   \
                             ta2_yy_y_xxx_0,  \
                             ta2_yy_y_xxx_1,  \
                             ta2_yy_y_xxy_0,  \
                             ta2_yy_y_xxy_1,  \
                             ta2_yy_y_xy_0,   \
                             ta2_yy_y_xy_1,   \
                             ta2_yy_y_xyy_0,  \
                             ta2_yy_y_xyy_1,  \
                             ta2_yy_y_xyz_0,  \
                             ta2_yy_y_xyz_1,  \
                             ta2_yy_y_yy_0,   \
                             ta2_yy_y_yy_1,   \
                             ta2_yy_y_yyy_0,  \
                             ta2_yy_y_yyy_1,  \
                             ta2_yy_y_yyz_0,  \
                             ta2_yy_y_yyz_1,  \
                             ta2_yy_y_yz_0,   \
                             ta2_yy_y_yz_1,   \
                             ta2_yy_y_yzz_0,  \
                             ta2_yy_y_yzz_1,  \
                             ta2_yy_yz_xxx_0, \
                             ta2_yy_yz_xxy_0, \
                             ta2_yy_yz_xxz_0, \
                             ta2_yy_yz_xyy_0, \
                             ta2_yy_yz_xyz_0, \
                             ta2_yy_yz_xzz_0, \
                             ta2_yy_yz_yyy_0, \
                             ta2_yy_yz_yyz_0, \
                             ta2_yy_yz_yzz_0, \
                             ta2_yy_yz_zzz_0, \
                             ta2_yy_z_xxz_0,  \
                             ta2_yy_z_xxz_1,  \
                             ta2_yy_z_xzz_0,  \
                             ta2_yy_z_xzz_1,  \
                             ta2_yy_z_zzz_0,  \
                             ta2_yy_z_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yz_xxx_0[i] = ta2_yy_y_xxx_0[i] * pa_z[i] - ta2_yy_y_xxx_1[i] * pc_z[i];

        ta2_yy_yz_xxy_0[i] = ta2_yy_y_xxy_0[i] * pa_z[i] - ta2_yy_y_xxy_1[i] * pc_z[i];

        ta2_yy_yz_xxz_0[i] = 2.0 * ta1_y_z_xxz_1[i] + ta2_yy_z_xxz_0[i] * pa_y[i] - ta2_yy_z_xxz_1[i] * pc_y[i];

        ta2_yy_yz_xyy_0[i] = ta2_yy_y_xyy_0[i] * pa_z[i] - ta2_yy_y_xyy_1[i] * pc_z[i];

        ta2_yy_yz_xyz_0[i] = ta2_yy_y_xy_0[i] * fe_0 - ta2_yy_y_xy_1[i] * fe_0 + ta2_yy_y_xyz_0[i] * pa_z[i] - ta2_yy_y_xyz_1[i] * pc_z[i];

        ta2_yy_yz_xzz_0[i] = 2.0 * ta1_y_z_xzz_1[i] + ta2_yy_z_xzz_0[i] * pa_y[i] - ta2_yy_z_xzz_1[i] * pc_y[i];

        ta2_yy_yz_yyy_0[i] = ta2_yy_y_yyy_0[i] * pa_z[i] - ta2_yy_y_yyy_1[i] * pc_z[i];

        ta2_yy_yz_yyz_0[i] = ta2_yy_y_yy_0[i] * fe_0 - ta2_yy_y_yy_1[i] * fe_0 + ta2_yy_y_yyz_0[i] * pa_z[i] - ta2_yy_y_yyz_1[i] * pc_z[i];

        ta2_yy_yz_yzz_0[i] =
            2.0 * ta2_yy_y_yz_0[i] * fe_0 - 2.0 * ta2_yy_y_yz_1[i] * fe_0 + ta2_yy_y_yzz_0[i] * pa_z[i] - ta2_yy_y_yzz_1[i] * pc_z[i];

        ta2_yy_yz_zzz_0[i] = 2.0 * ta1_y_z_zzz_1[i] + ta2_yy_z_zzz_0[i] * pa_y[i] - ta2_yy_z_zzz_1[i] * pc_y[i];
    }

    // Set up 230-240 components of targeted buffer : DF

    auto ta2_yy_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 230);

    auto ta2_yy_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 231);

    auto ta2_yy_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 232);

    auto ta2_yy_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 233);

    auto ta2_yy_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 234);

    auto ta2_yy_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 235);

    auto ta2_yy_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 236);

    auto ta2_yy_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 237);

    auto ta2_yy_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 238);

    auto ta2_yy_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 239);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta2_yy_0_xxx_0,  \
                             ta2_yy_0_xxx_1,  \
                             ta2_yy_0_xxy_0,  \
                             ta2_yy_0_xxy_1,  \
                             ta2_yy_0_xxz_0,  \
                             ta2_yy_0_xxz_1,  \
                             ta2_yy_0_xyy_0,  \
                             ta2_yy_0_xyy_1,  \
                             ta2_yy_0_xyz_0,  \
                             ta2_yy_0_xyz_1,  \
                             ta2_yy_0_xzz_0,  \
                             ta2_yy_0_xzz_1,  \
                             ta2_yy_0_yyy_0,  \
                             ta2_yy_0_yyy_1,  \
                             ta2_yy_0_yyz_0,  \
                             ta2_yy_0_yyz_1,  \
                             ta2_yy_0_yzz_0,  \
                             ta2_yy_0_yzz_1,  \
                             ta2_yy_0_zzz_0,  \
                             ta2_yy_0_zzz_1,  \
                             ta2_yy_z_xx_0,   \
                             ta2_yy_z_xx_1,   \
                             ta2_yy_z_xxx_0,  \
                             ta2_yy_z_xxx_1,  \
                             ta2_yy_z_xxy_0,  \
                             ta2_yy_z_xxy_1,  \
                             ta2_yy_z_xxz_0,  \
                             ta2_yy_z_xxz_1,  \
                             ta2_yy_z_xy_0,   \
                             ta2_yy_z_xy_1,   \
                             ta2_yy_z_xyy_0,  \
                             ta2_yy_z_xyy_1,  \
                             ta2_yy_z_xyz_0,  \
                             ta2_yy_z_xyz_1,  \
                             ta2_yy_z_xz_0,   \
                             ta2_yy_z_xz_1,   \
                             ta2_yy_z_xzz_0,  \
                             ta2_yy_z_xzz_1,  \
                             ta2_yy_z_yy_0,   \
                             ta2_yy_z_yy_1,   \
                             ta2_yy_z_yyy_0,  \
                             ta2_yy_z_yyy_1,  \
                             ta2_yy_z_yyz_0,  \
                             ta2_yy_z_yyz_1,  \
                             ta2_yy_z_yz_0,   \
                             ta2_yy_z_yz_1,   \
                             ta2_yy_z_yzz_0,  \
                             ta2_yy_z_yzz_1,  \
                             ta2_yy_z_zz_0,   \
                             ta2_yy_z_zz_1,   \
                             ta2_yy_z_zzz_0,  \
                             ta2_yy_z_zzz_1,  \
                             ta2_yy_zz_xxx_0, \
                             ta2_yy_zz_xxy_0, \
                             ta2_yy_zz_xxz_0, \
                             ta2_yy_zz_xyy_0, \
                             ta2_yy_zz_xyz_0, \
                             ta2_yy_zz_xzz_0, \
                             ta2_yy_zz_yyy_0, \
                             ta2_yy_zz_yyz_0, \
                             ta2_yy_zz_yzz_0, \
                             ta2_yy_zz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_zz_xxx_0[i] = ta2_yy_0_xxx_0[i] * fe_0 - ta2_yy_0_xxx_1[i] * fe_0 + ta2_yy_z_xxx_0[i] * pa_z[i] - ta2_yy_z_xxx_1[i] * pc_z[i];

        ta2_yy_zz_xxy_0[i] = ta2_yy_0_xxy_0[i] * fe_0 - ta2_yy_0_xxy_1[i] * fe_0 + ta2_yy_z_xxy_0[i] * pa_z[i] - ta2_yy_z_xxy_1[i] * pc_z[i];

        ta2_yy_zz_xxz_0[i] = ta2_yy_0_xxz_0[i] * fe_0 - ta2_yy_0_xxz_1[i] * fe_0 + ta2_yy_z_xx_0[i] * fe_0 - ta2_yy_z_xx_1[i] * fe_0 +
                             ta2_yy_z_xxz_0[i] * pa_z[i] - ta2_yy_z_xxz_1[i] * pc_z[i];

        ta2_yy_zz_xyy_0[i] = ta2_yy_0_xyy_0[i] * fe_0 - ta2_yy_0_xyy_1[i] * fe_0 + ta2_yy_z_xyy_0[i] * pa_z[i] - ta2_yy_z_xyy_1[i] * pc_z[i];

        ta2_yy_zz_xyz_0[i] = ta2_yy_0_xyz_0[i] * fe_0 - ta2_yy_0_xyz_1[i] * fe_0 + ta2_yy_z_xy_0[i] * fe_0 - ta2_yy_z_xy_1[i] * fe_0 +
                             ta2_yy_z_xyz_0[i] * pa_z[i] - ta2_yy_z_xyz_1[i] * pc_z[i];

        ta2_yy_zz_xzz_0[i] = ta2_yy_0_xzz_0[i] * fe_0 - ta2_yy_0_xzz_1[i] * fe_0 + 2.0 * ta2_yy_z_xz_0[i] * fe_0 - 2.0 * ta2_yy_z_xz_1[i] * fe_0 +
                             ta2_yy_z_xzz_0[i] * pa_z[i] - ta2_yy_z_xzz_1[i] * pc_z[i];

        ta2_yy_zz_yyy_0[i] = ta2_yy_0_yyy_0[i] * fe_0 - ta2_yy_0_yyy_1[i] * fe_0 + ta2_yy_z_yyy_0[i] * pa_z[i] - ta2_yy_z_yyy_1[i] * pc_z[i];

        ta2_yy_zz_yyz_0[i] = ta2_yy_0_yyz_0[i] * fe_0 - ta2_yy_0_yyz_1[i] * fe_0 + ta2_yy_z_yy_0[i] * fe_0 - ta2_yy_z_yy_1[i] * fe_0 +
                             ta2_yy_z_yyz_0[i] * pa_z[i] - ta2_yy_z_yyz_1[i] * pc_z[i];

        ta2_yy_zz_yzz_0[i] = ta2_yy_0_yzz_0[i] * fe_0 - ta2_yy_0_yzz_1[i] * fe_0 + 2.0 * ta2_yy_z_yz_0[i] * fe_0 - 2.0 * ta2_yy_z_yz_1[i] * fe_0 +
                             ta2_yy_z_yzz_0[i] * pa_z[i] - ta2_yy_z_yzz_1[i] * pc_z[i];

        ta2_yy_zz_zzz_0[i] = ta2_yy_0_zzz_0[i] * fe_0 - ta2_yy_0_zzz_1[i] * fe_0 + 3.0 * ta2_yy_z_zz_0[i] * fe_0 - 3.0 * ta2_yy_z_zz_1[i] * fe_0 +
                             ta2_yy_z_zzz_0[i] * pa_z[i] - ta2_yy_z_zzz_1[i] * pc_z[i];
    }

    // Set up 240-250 components of targeted buffer : DF

    auto ta2_yz_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 240);

    auto ta2_yz_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 241);

    auto ta2_yz_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 242);

    auto ta2_yz_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 243);

    auto ta2_yz_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 244);

    auto ta2_yz_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 245);

    auto ta2_yz_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 246);

    auto ta2_yz_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 247);

    auto ta2_yz_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 248);

    auto ta2_yz_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 249);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta2_yz_0_xxx_0,  \
                             ta2_yz_0_xxx_1,  \
                             ta2_yz_0_xxy_0,  \
                             ta2_yz_0_xxy_1,  \
                             ta2_yz_0_xxz_0,  \
                             ta2_yz_0_xxz_1,  \
                             ta2_yz_0_xyy_0,  \
                             ta2_yz_0_xyy_1,  \
                             ta2_yz_0_xyz_0,  \
                             ta2_yz_0_xyz_1,  \
                             ta2_yz_0_xzz_0,  \
                             ta2_yz_0_xzz_1,  \
                             ta2_yz_0_yyy_0,  \
                             ta2_yz_0_yyy_1,  \
                             ta2_yz_0_yyz_0,  \
                             ta2_yz_0_yyz_1,  \
                             ta2_yz_0_yzz_0,  \
                             ta2_yz_0_yzz_1,  \
                             ta2_yz_0_zzz_0,  \
                             ta2_yz_0_zzz_1,  \
                             ta2_yz_x_xx_0,   \
                             ta2_yz_x_xx_1,   \
                             ta2_yz_x_xxx_0,  \
                             ta2_yz_x_xxx_1,  \
                             ta2_yz_x_xxy_0,  \
                             ta2_yz_x_xxy_1,  \
                             ta2_yz_x_xxz_0,  \
                             ta2_yz_x_xxz_1,  \
                             ta2_yz_x_xy_0,   \
                             ta2_yz_x_xy_1,   \
                             ta2_yz_x_xyy_0,  \
                             ta2_yz_x_xyy_1,  \
                             ta2_yz_x_xyz_0,  \
                             ta2_yz_x_xyz_1,  \
                             ta2_yz_x_xz_0,   \
                             ta2_yz_x_xz_1,   \
                             ta2_yz_x_xzz_0,  \
                             ta2_yz_x_xzz_1,  \
                             ta2_yz_x_yy_0,   \
                             ta2_yz_x_yy_1,   \
                             ta2_yz_x_yyy_0,  \
                             ta2_yz_x_yyy_1,  \
                             ta2_yz_x_yyz_0,  \
                             ta2_yz_x_yyz_1,  \
                             ta2_yz_x_yz_0,   \
                             ta2_yz_x_yz_1,   \
                             ta2_yz_x_yzz_0,  \
                             ta2_yz_x_yzz_1,  \
                             ta2_yz_x_zz_0,   \
                             ta2_yz_x_zz_1,   \
                             ta2_yz_x_zzz_0,  \
                             ta2_yz_x_zzz_1,  \
                             ta2_yz_xx_xxx_0, \
                             ta2_yz_xx_xxy_0, \
                             ta2_yz_xx_xxz_0, \
                             ta2_yz_xx_xyy_0, \
                             ta2_yz_xx_xyz_0, \
                             ta2_yz_xx_xzz_0, \
                             ta2_yz_xx_yyy_0, \
                             ta2_yz_xx_yyz_0, \
                             ta2_yz_xx_yzz_0, \
                             ta2_yz_xx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xx_xxx_0[i] = ta2_yz_0_xxx_0[i] * fe_0 - ta2_yz_0_xxx_1[i] * fe_0 + 3.0 * ta2_yz_x_xx_0[i] * fe_0 - 3.0 * ta2_yz_x_xx_1[i] * fe_0 +
                             ta2_yz_x_xxx_0[i] * pa_x[i] - ta2_yz_x_xxx_1[i] * pc_x[i];

        ta2_yz_xx_xxy_0[i] = ta2_yz_0_xxy_0[i] * fe_0 - ta2_yz_0_xxy_1[i] * fe_0 + 2.0 * ta2_yz_x_xy_0[i] * fe_0 - 2.0 * ta2_yz_x_xy_1[i] * fe_0 +
                             ta2_yz_x_xxy_0[i] * pa_x[i] - ta2_yz_x_xxy_1[i] * pc_x[i];

        ta2_yz_xx_xxz_0[i] = ta2_yz_0_xxz_0[i] * fe_0 - ta2_yz_0_xxz_1[i] * fe_0 + 2.0 * ta2_yz_x_xz_0[i] * fe_0 - 2.0 * ta2_yz_x_xz_1[i] * fe_0 +
                             ta2_yz_x_xxz_0[i] * pa_x[i] - ta2_yz_x_xxz_1[i] * pc_x[i];

        ta2_yz_xx_xyy_0[i] = ta2_yz_0_xyy_0[i] * fe_0 - ta2_yz_0_xyy_1[i] * fe_0 + ta2_yz_x_yy_0[i] * fe_0 - ta2_yz_x_yy_1[i] * fe_0 +
                             ta2_yz_x_xyy_0[i] * pa_x[i] - ta2_yz_x_xyy_1[i] * pc_x[i];

        ta2_yz_xx_xyz_0[i] = ta2_yz_0_xyz_0[i] * fe_0 - ta2_yz_0_xyz_1[i] * fe_0 + ta2_yz_x_yz_0[i] * fe_0 - ta2_yz_x_yz_1[i] * fe_0 +
                             ta2_yz_x_xyz_0[i] * pa_x[i] - ta2_yz_x_xyz_1[i] * pc_x[i];

        ta2_yz_xx_xzz_0[i] = ta2_yz_0_xzz_0[i] * fe_0 - ta2_yz_0_xzz_1[i] * fe_0 + ta2_yz_x_zz_0[i] * fe_0 - ta2_yz_x_zz_1[i] * fe_0 +
                             ta2_yz_x_xzz_0[i] * pa_x[i] - ta2_yz_x_xzz_1[i] * pc_x[i];

        ta2_yz_xx_yyy_0[i] = ta2_yz_0_yyy_0[i] * fe_0 - ta2_yz_0_yyy_1[i] * fe_0 + ta2_yz_x_yyy_0[i] * pa_x[i] - ta2_yz_x_yyy_1[i] * pc_x[i];

        ta2_yz_xx_yyz_0[i] = ta2_yz_0_yyz_0[i] * fe_0 - ta2_yz_0_yyz_1[i] * fe_0 + ta2_yz_x_yyz_0[i] * pa_x[i] - ta2_yz_x_yyz_1[i] * pc_x[i];

        ta2_yz_xx_yzz_0[i] = ta2_yz_0_yzz_0[i] * fe_0 - ta2_yz_0_yzz_1[i] * fe_0 + ta2_yz_x_yzz_0[i] * pa_x[i] - ta2_yz_x_yzz_1[i] * pc_x[i];

        ta2_yz_xx_zzz_0[i] = ta2_yz_0_zzz_0[i] * fe_0 - ta2_yz_0_zzz_1[i] * fe_0 + ta2_yz_x_zzz_0[i] * pa_x[i] - ta2_yz_x_zzz_1[i] * pc_x[i];
    }

    // Set up 250-260 components of targeted buffer : DF

    auto ta2_yz_xy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 250);

    auto ta2_yz_xy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 251);

    auto ta2_yz_xy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 252);

    auto ta2_yz_xy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 253);

    auto ta2_yz_xy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 254);

    auto ta2_yz_xy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 255);

    auto ta2_yz_xy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 256);

    auto ta2_yz_xy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 257);

    auto ta2_yz_xy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 258);

    auto ta2_yz_xy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 259);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta1_z_x_xxx_1,   \
                             ta1_z_x_xxz_1,   \
                             ta1_z_x_xzz_1,   \
                             ta2_yz_x_xxx_0,  \
                             ta2_yz_x_xxx_1,  \
                             ta2_yz_x_xxz_0,  \
                             ta2_yz_x_xxz_1,  \
                             ta2_yz_x_xzz_0,  \
                             ta2_yz_x_xzz_1,  \
                             ta2_yz_xy_xxx_0, \
                             ta2_yz_xy_xxy_0, \
                             ta2_yz_xy_xxz_0, \
                             ta2_yz_xy_xyy_0, \
                             ta2_yz_xy_xyz_0, \
                             ta2_yz_xy_xzz_0, \
                             ta2_yz_xy_yyy_0, \
                             ta2_yz_xy_yyz_0, \
                             ta2_yz_xy_yzz_0, \
                             ta2_yz_xy_zzz_0, \
                             ta2_yz_y_xxy_0,  \
                             ta2_yz_y_xxy_1,  \
                             ta2_yz_y_xy_0,   \
                             ta2_yz_y_xy_1,   \
                             ta2_yz_y_xyy_0,  \
                             ta2_yz_y_xyy_1,  \
                             ta2_yz_y_xyz_0,  \
                             ta2_yz_y_xyz_1,  \
                             ta2_yz_y_yy_0,   \
                             ta2_yz_y_yy_1,   \
                             ta2_yz_y_yyy_0,  \
                             ta2_yz_y_yyy_1,  \
                             ta2_yz_y_yyz_0,  \
                             ta2_yz_y_yyz_1,  \
                             ta2_yz_y_yz_0,   \
                             ta2_yz_y_yz_1,   \
                             ta2_yz_y_yzz_0,  \
                             ta2_yz_y_yzz_1,  \
                             ta2_yz_y_zzz_0,  \
                             ta2_yz_y_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xy_xxx_0[i] = ta1_z_x_xxx_1[i] + ta2_yz_x_xxx_0[i] * pa_y[i] - ta2_yz_x_xxx_1[i] * pc_y[i];

        ta2_yz_xy_xxy_0[i] =
            2.0 * ta2_yz_y_xy_0[i] * fe_0 - 2.0 * ta2_yz_y_xy_1[i] * fe_0 + ta2_yz_y_xxy_0[i] * pa_x[i] - ta2_yz_y_xxy_1[i] * pc_x[i];

        ta2_yz_xy_xxz_0[i] = ta1_z_x_xxz_1[i] + ta2_yz_x_xxz_0[i] * pa_y[i] - ta2_yz_x_xxz_1[i] * pc_y[i];

        ta2_yz_xy_xyy_0[i] = ta2_yz_y_yy_0[i] * fe_0 - ta2_yz_y_yy_1[i] * fe_0 + ta2_yz_y_xyy_0[i] * pa_x[i] - ta2_yz_y_xyy_1[i] * pc_x[i];

        ta2_yz_xy_xyz_0[i] = ta2_yz_y_yz_0[i] * fe_0 - ta2_yz_y_yz_1[i] * fe_0 + ta2_yz_y_xyz_0[i] * pa_x[i] - ta2_yz_y_xyz_1[i] * pc_x[i];

        ta2_yz_xy_xzz_0[i] = ta1_z_x_xzz_1[i] + ta2_yz_x_xzz_0[i] * pa_y[i] - ta2_yz_x_xzz_1[i] * pc_y[i];

        ta2_yz_xy_yyy_0[i] = ta2_yz_y_yyy_0[i] * pa_x[i] - ta2_yz_y_yyy_1[i] * pc_x[i];

        ta2_yz_xy_yyz_0[i] = ta2_yz_y_yyz_0[i] * pa_x[i] - ta2_yz_y_yyz_1[i] * pc_x[i];

        ta2_yz_xy_yzz_0[i] = ta2_yz_y_yzz_0[i] * pa_x[i] - ta2_yz_y_yzz_1[i] * pc_x[i];

        ta2_yz_xy_zzz_0[i] = ta2_yz_y_zzz_0[i] * pa_x[i] - ta2_yz_y_zzz_1[i] * pc_x[i];
    }

    // Set up 260-270 components of targeted buffer : DF

    auto ta2_yz_xz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 260);

    auto ta2_yz_xz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 261);

    auto ta2_yz_xz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 262);

    auto ta2_yz_xz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 263);

    auto ta2_yz_xz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 264);

    auto ta2_yz_xz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 265);

    auto ta2_yz_xz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 266);

    auto ta2_yz_xz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 267);

    auto ta2_yz_xz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 268);

    auto ta2_yz_xz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 269);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_y_x_xxx_1,   \
                             ta1_y_x_xxy_1,   \
                             ta1_y_x_xyy_1,   \
                             ta2_yz_x_xxx_0,  \
                             ta2_yz_x_xxx_1,  \
                             ta2_yz_x_xxy_0,  \
                             ta2_yz_x_xxy_1,  \
                             ta2_yz_x_xyy_0,  \
                             ta2_yz_x_xyy_1,  \
                             ta2_yz_xz_xxx_0, \
                             ta2_yz_xz_xxy_0, \
                             ta2_yz_xz_xxz_0, \
                             ta2_yz_xz_xyy_0, \
                             ta2_yz_xz_xyz_0, \
                             ta2_yz_xz_xzz_0, \
                             ta2_yz_xz_yyy_0, \
                             ta2_yz_xz_yyz_0, \
                             ta2_yz_xz_yzz_0, \
                             ta2_yz_xz_zzz_0, \
                             ta2_yz_z_xxz_0,  \
                             ta2_yz_z_xxz_1,  \
                             ta2_yz_z_xyz_0,  \
                             ta2_yz_z_xyz_1,  \
                             ta2_yz_z_xz_0,   \
                             ta2_yz_z_xz_1,   \
                             ta2_yz_z_xzz_0,  \
                             ta2_yz_z_xzz_1,  \
                             ta2_yz_z_yyy_0,  \
                             ta2_yz_z_yyy_1,  \
                             ta2_yz_z_yyz_0,  \
                             ta2_yz_z_yyz_1,  \
                             ta2_yz_z_yz_0,   \
                             ta2_yz_z_yz_1,   \
                             ta2_yz_z_yzz_0,  \
                             ta2_yz_z_yzz_1,  \
                             ta2_yz_z_zz_0,   \
                             ta2_yz_z_zz_1,   \
                             ta2_yz_z_zzz_0,  \
                             ta2_yz_z_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xz_xxx_0[i] = ta1_y_x_xxx_1[i] + ta2_yz_x_xxx_0[i] * pa_z[i] - ta2_yz_x_xxx_1[i] * pc_z[i];

        ta2_yz_xz_xxy_0[i] = ta1_y_x_xxy_1[i] + ta2_yz_x_xxy_0[i] * pa_z[i] - ta2_yz_x_xxy_1[i] * pc_z[i];

        ta2_yz_xz_xxz_0[i] =
            2.0 * ta2_yz_z_xz_0[i] * fe_0 - 2.0 * ta2_yz_z_xz_1[i] * fe_0 + ta2_yz_z_xxz_0[i] * pa_x[i] - ta2_yz_z_xxz_1[i] * pc_x[i];

        ta2_yz_xz_xyy_0[i] = ta1_y_x_xyy_1[i] + ta2_yz_x_xyy_0[i] * pa_z[i] - ta2_yz_x_xyy_1[i] * pc_z[i];

        ta2_yz_xz_xyz_0[i] = ta2_yz_z_yz_0[i] * fe_0 - ta2_yz_z_yz_1[i] * fe_0 + ta2_yz_z_xyz_0[i] * pa_x[i] - ta2_yz_z_xyz_1[i] * pc_x[i];

        ta2_yz_xz_xzz_0[i] = ta2_yz_z_zz_0[i] * fe_0 - ta2_yz_z_zz_1[i] * fe_0 + ta2_yz_z_xzz_0[i] * pa_x[i] - ta2_yz_z_xzz_1[i] * pc_x[i];

        ta2_yz_xz_yyy_0[i] = ta2_yz_z_yyy_0[i] * pa_x[i] - ta2_yz_z_yyy_1[i] * pc_x[i];

        ta2_yz_xz_yyz_0[i] = ta2_yz_z_yyz_0[i] * pa_x[i] - ta2_yz_z_yyz_1[i] * pc_x[i];

        ta2_yz_xz_yzz_0[i] = ta2_yz_z_yzz_0[i] * pa_x[i] - ta2_yz_z_yzz_1[i] * pc_x[i];

        ta2_yz_xz_zzz_0[i] = ta2_yz_z_zzz_0[i] * pa_x[i] - ta2_yz_z_zzz_1[i] * pc_x[i];
    }

    // Set up 270-280 components of targeted buffer : DF

    auto ta2_yz_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 270);

    auto ta2_yz_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 271);

    auto ta2_yz_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 272);

    auto ta2_yz_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 273);

    auto ta2_yz_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 274);

    auto ta2_yz_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 275);

    auto ta2_yz_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 276);

    auto ta2_yz_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 277);

    auto ta2_yz_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 278);

    auto ta2_yz_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 279);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta1_z_y_xxx_1,   \
                             ta1_z_y_xxy_1,   \
                             ta1_z_y_xxz_1,   \
                             ta1_z_y_xyy_1,   \
                             ta1_z_y_xyz_1,   \
                             ta1_z_y_xzz_1,   \
                             ta1_z_y_yyy_1,   \
                             ta1_z_y_yyz_1,   \
                             ta1_z_y_yzz_1,   \
                             ta1_z_y_zzz_1,   \
                             ta2_yz_0_xxx_0,  \
                             ta2_yz_0_xxx_1,  \
                             ta2_yz_0_xxy_0,  \
                             ta2_yz_0_xxy_1,  \
                             ta2_yz_0_xxz_0,  \
                             ta2_yz_0_xxz_1,  \
                             ta2_yz_0_xyy_0,  \
                             ta2_yz_0_xyy_1,  \
                             ta2_yz_0_xyz_0,  \
                             ta2_yz_0_xyz_1,  \
                             ta2_yz_0_xzz_0,  \
                             ta2_yz_0_xzz_1,  \
                             ta2_yz_0_yyy_0,  \
                             ta2_yz_0_yyy_1,  \
                             ta2_yz_0_yyz_0,  \
                             ta2_yz_0_yyz_1,  \
                             ta2_yz_0_yzz_0,  \
                             ta2_yz_0_yzz_1,  \
                             ta2_yz_0_zzz_0,  \
                             ta2_yz_0_zzz_1,  \
                             ta2_yz_y_xx_0,   \
                             ta2_yz_y_xx_1,   \
                             ta2_yz_y_xxx_0,  \
                             ta2_yz_y_xxx_1,  \
                             ta2_yz_y_xxy_0,  \
                             ta2_yz_y_xxy_1,  \
                             ta2_yz_y_xxz_0,  \
                             ta2_yz_y_xxz_1,  \
                             ta2_yz_y_xy_0,   \
                             ta2_yz_y_xy_1,   \
                             ta2_yz_y_xyy_0,  \
                             ta2_yz_y_xyy_1,  \
                             ta2_yz_y_xyz_0,  \
                             ta2_yz_y_xyz_1,  \
                             ta2_yz_y_xz_0,   \
                             ta2_yz_y_xz_1,   \
                             ta2_yz_y_xzz_0,  \
                             ta2_yz_y_xzz_1,  \
                             ta2_yz_y_yy_0,   \
                             ta2_yz_y_yy_1,   \
                             ta2_yz_y_yyy_0,  \
                             ta2_yz_y_yyy_1,  \
                             ta2_yz_y_yyz_0,  \
                             ta2_yz_y_yyz_1,  \
                             ta2_yz_y_yz_0,   \
                             ta2_yz_y_yz_1,   \
                             ta2_yz_y_yzz_0,  \
                             ta2_yz_y_yzz_1,  \
                             ta2_yz_y_zz_0,   \
                             ta2_yz_y_zz_1,   \
                             ta2_yz_y_zzz_0,  \
                             ta2_yz_y_zzz_1,  \
                             ta2_yz_yy_xxx_0, \
                             ta2_yz_yy_xxy_0, \
                             ta2_yz_yy_xxz_0, \
                             ta2_yz_yy_xyy_0, \
                             ta2_yz_yy_xyz_0, \
                             ta2_yz_yy_xzz_0, \
                             ta2_yz_yy_yyy_0, \
                             ta2_yz_yy_yyz_0, \
                             ta2_yz_yy_yzz_0, \
                             ta2_yz_yy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yy_xxx_0[i] =
            ta2_yz_0_xxx_0[i] * fe_0 - ta2_yz_0_xxx_1[i] * fe_0 + ta1_z_y_xxx_1[i] + ta2_yz_y_xxx_0[i] * pa_y[i] - ta2_yz_y_xxx_1[i] * pc_y[i];

        ta2_yz_yy_xxy_0[i] = ta2_yz_0_xxy_0[i] * fe_0 - ta2_yz_0_xxy_1[i] * fe_0 + ta2_yz_y_xx_0[i] * fe_0 - ta2_yz_y_xx_1[i] * fe_0 +
                             ta1_z_y_xxy_1[i] + ta2_yz_y_xxy_0[i] * pa_y[i] - ta2_yz_y_xxy_1[i] * pc_y[i];

        ta2_yz_yy_xxz_0[i] =
            ta2_yz_0_xxz_0[i] * fe_0 - ta2_yz_0_xxz_1[i] * fe_0 + ta1_z_y_xxz_1[i] + ta2_yz_y_xxz_0[i] * pa_y[i] - ta2_yz_y_xxz_1[i] * pc_y[i];

        ta2_yz_yy_xyy_0[i] = ta2_yz_0_xyy_0[i] * fe_0 - ta2_yz_0_xyy_1[i] * fe_0 + 2.0 * ta2_yz_y_xy_0[i] * fe_0 - 2.0 * ta2_yz_y_xy_1[i] * fe_0 +
                             ta1_z_y_xyy_1[i] + ta2_yz_y_xyy_0[i] * pa_y[i] - ta2_yz_y_xyy_1[i] * pc_y[i];

        ta2_yz_yy_xyz_0[i] = ta2_yz_0_xyz_0[i] * fe_0 - ta2_yz_0_xyz_1[i] * fe_0 + ta2_yz_y_xz_0[i] * fe_0 - ta2_yz_y_xz_1[i] * fe_0 +
                             ta1_z_y_xyz_1[i] + ta2_yz_y_xyz_0[i] * pa_y[i] - ta2_yz_y_xyz_1[i] * pc_y[i];

        ta2_yz_yy_xzz_0[i] =
            ta2_yz_0_xzz_0[i] * fe_0 - ta2_yz_0_xzz_1[i] * fe_0 + ta1_z_y_xzz_1[i] + ta2_yz_y_xzz_0[i] * pa_y[i] - ta2_yz_y_xzz_1[i] * pc_y[i];

        ta2_yz_yy_yyy_0[i] = ta2_yz_0_yyy_0[i] * fe_0 - ta2_yz_0_yyy_1[i] * fe_0 + 3.0 * ta2_yz_y_yy_0[i] * fe_0 - 3.0 * ta2_yz_y_yy_1[i] * fe_0 +
                             ta1_z_y_yyy_1[i] + ta2_yz_y_yyy_0[i] * pa_y[i] - ta2_yz_y_yyy_1[i] * pc_y[i];

        ta2_yz_yy_yyz_0[i] = ta2_yz_0_yyz_0[i] * fe_0 - ta2_yz_0_yyz_1[i] * fe_0 + 2.0 * ta2_yz_y_yz_0[i] * fe_0 - 2.0 * ta2_yz_y_yz_1[i] * fe_0 +
                             ta1_z_y_yyz_1[i] + ta2_yz_y_yyz_0[i] * pa_y[i] - ta2_yz_y_yyz_1[i] * pc_y[i];

        ta2_yz_yy_yzz_0[i] = ta2_yz_0_yzz_0[i] * fe_0 - ta2_yz_0_yzz_1[i] * fe_0 + ta2_yz_y_zz_0[i] * fe_0 - ta2_yz_y_zz_1[i] * fe_0 +
                             ta1_z_y_yzz_1[i] + ta2_yz_y_yzz_0[i] * pa_y[i] - ta2_yz_y_yzz_1[i] * pc_y[i];

        ta2_yz_yy_zzz_0[i] =
            ta2_yz_0_zzz_0[i] * fe_0 - ta2_yz_0_zzz_1[i] * fe_0 + ta1_z_y_zzz_1[i] + ta2_yz_y_zzz_0[i] * pa_y[i] - ta2_yz_y_zzz_1[i] * pc_y[i];
    }

    // Set up 280-290 components of targeted buffer : DF

    auto ta2_yz_yz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 280);

    auto ta2_yz_yz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 281);

    auto ta2_yz_yz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 282);

    auto ta2_yz_yz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 283);

    auto ta2_yz_yz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 284);

    auto ta2_yz_yz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 285);

    auto ta2_yz_yz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 286);

    auto ta2_yz_yz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 287);

    auto ta2_yz_yz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 288);

    auto ta2_yz_yz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 289);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_y_y_xxy_1,   \
                             ta1_y_y_xyy_1,   \
                             ta1_y_y_yyy_1,   \
                             ta1_z_z_xxx_1,   \
                             ta1_z_z_xxz_1,   \
                             ta1_z_z_xyz_1,   \
                             ta1_z_z_xzz_1,   \
                             ta1_z_z_yyz_1,   \
                             ta1_z_z_yzz_1,   \
                             ta1_z_z_zzz_1,   \
                             ta2_yz_y_xxy_0,  \
                             ta2_yz_y_xxy_1,  \
                             ta2_yz_y_xyy_0,  \
                             ta2_yz_y_xyy_1,  \
                             ta2_yz_y_yyy_0,  \
                             ta2_yz_y_yyy_1,  \
                             ta2_yz_yz_xxx_0, \
                             ta2_yz_yz_xxy_0, \
                             ta2_yz_yz_xxz_0, \
                             ta2_yz_yz_xyy_0, \
                             ta2_yz_yz_xyz_0, \
                             ta2_yz_yz_xzz_0, \
                             ta2_yz_yz_yyy_0, \
                             ta2_yz_yz_yyz_0, \
                             ta2_yz_yz_yzz_0, \
                             ta2_yz_yz_zzz_0, \
                             ta2_yz_z_xxx_0,  \
                             ta2_yz_z_xxx_1,  \
                             ta2_yz_z_xxz_0,  \
                             ta2_yz_z_xxz_1,  \
                             ta2_yz_z_xyz_0,  \
                             ta2_yz_z_xyz_1,  \
                             ta2_yz_z_xz_0,   \
                             ta2_yz_z_xz_1,   \
                             ta2_yz_z_xzz_0,  \
                             ta2_yz_z_xzz_1,  \
                             ta2_yz_z_yyz_0,  \
                             ta2_yz_z_yyz_1,  \
                             ta2_yz_z_yz_0,   \
                             ta2_yz_z_yz_1,   \
                             ta2_yz_z_yzz_0,  \
                             ta2_yz_z_yzz_1,  \
                             ta2_yz_z_zz_0,   \
                             ta2_yz_z_zz_1,   \
                             ta2_yz_z_zzz_0,  \
                             ta2_yz_z_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yz_xxx_0[i] = ta1_z_z_xxx_1[i] + ta2_yz_z_xxx_0[i] * pa_y[i] - ta2_yz_z_xxx_1[i] * pc_y[i];

        ta2_yz_yz_xxy_0[i] = ta1_y_y_xxy_1[i] + ta2_yz_y_xxy_0[i] * pa_z[i] - ta2_yz_y_xxy_1[i] * pc_z[i];

        ta2_yz_yz_xxz_0[i] = ta1_z_z_xxz_1[i] + ta2_yz_z_xxz_0[i] * pa_y[i] - ta2_yz_z_xxz_1[i] * pc_y[i];

        ta2_yz_yz_xyy_0[i] = ta1_y_y_xyy_1[i] + ta2_yz_y_xyy_0[i] * pa_z[i] - ta2_yz_y_xyy_1[i] * pc_z[i];

        ta2_yz_yz_xyz_0[i] =
            ta2_yz_z_xz_0[i] * fe_0 - ta2_yz_z_xz_1[i] * fe_0 + ta1_z_z_xyz_1[i] + ta2_yz_z_xyz_0[i] * pa_y[i] - ta2_yz_z_xyz_1[i] * pc_y[i];

        ta2_yz_yz_xzz_0[i] = ta1_z_z_xzz_1[i] + ta2_yz_z_xzz_0[i] * pa_y[i] - ta2_yz_z_xzz_1[i] * pc_y[i];

        ta2_yz_yz_yyy_0[i] = ta1_y_y_yyy_1[i] + ta2_yz_y_yyy_0[i] * pa_z[i] - ta2_yz_y_yyy_1[i] * pc_z[i];

        ta2_yz_yz_yyz_0[i] = 2.0 * ta2_yz_z_yz_0[i] * fe_0 - 2.0 * ta2_yz_z_yz_1[i] * fe_0 + ta1_z_z_yyz_1[i] + ta2_yz_z_yyz_0[i] * pa_y[i] -
                             ta2_yz_z_yyz_1[i] * pc_y[i];

        ta2_yz_yz_yzz_0[i] =
            ta2_yz_z_zz_0[i] * fe_0 - ta2_yz_z_zz_1[i] * fe_0 + ta1_z_z_yzz_1[i] + ta2_yz_z_yzz_0[i] * pa_y[i] - ta2_yz_z_yzz_1[i] * pc_y[i];

        ta2_yz_yz_zzz_0[i] = ta1_z_z_zzz_1[i] + ta2_yz_z_zzz_0[i] * pa_y[i] - ta2_yz_z_zzz_1[i] * pc_y[i];
    }

    // Set up 290-300 components of targeted buffer : DF

    auto ta2_yz_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 290);

    auto ta2_yz_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 291);

    auto ta2_yz_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 292);

    auto ta2_yz_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 293);

    auto ta2_yz_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 294);

    auto ta2_yz_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 295);

    auto ta2_yz_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 296);

    auto ta2_yz_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 297);

    auto ta2_yz_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 298);

    auto ta2_yz_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 299);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_y_z_xxx_1,   \
                             ta1_y_z_xxy_1,   \
                             ta1_y_z_xxz_1,   \
                             ta1_y_z_xyy_1,   \
                             ta1_y_z_xyz_1,   \
                             ta1_y_z_xzz_1,   \
                             ta1_y_z_yyy_1,   \
                             ta1_y_z_yyz_1,   \
                             ta1_y_z_yzz_1,   \
                             ta1_y_z_zzz_1,   \
                             ta2_yz_0_xxx_0,  \
                             ta2_yz_0_xxx_1,  \
                             ta2_yz_0_xxy_0,  \
                             ta2_yz_0_xxy_1,  \
                             ta2_yz_0_xxz_0,  \
                             ta2_yz_0_xxz_1,  \
                             ta2_yz_0_xyy_0,  \
                             ta2_yz_0_xyy_1,  \
                             ta2_yz_0_xyz_0,  \
                             ta2_yz_0_xyz_1,  \
                             ta2_yz_0_xzz_0,  \
                             ta2_yz_0_xzz_1,  \
                             ta2_yz_0_yyy_0,  \
                             ta2_yz_0_yyy_1,  \
                             ta2_yz_0_yyz_0,  \
                             ta2_yz_0_yyz_1,  \
                             ta2_yz_0_yzz_0,  \
                             ta2_yz_0_yzz_1,  \
                             ta2_yz_0_zzz_0,  \
                             ta2_yz_0_zzz_1,  \
                             ta2_yz_z_xx_0,   \
                             ta2_yz_z_xx_1,   \
                             ta2_yz_z_xxx_0,  \
                             ta2_yz_z_xxx_1,  \
                             ta2_yz_z_xxy_0,  \
                             ta2_yz_z_xxy_1,  \
                             ta2_yz_z_xxz_0,  \
                             ta2_yz_z_xxz_1,  \
                             ta2_yz_z_xy_0,   \
                             ta2_yz_z_xy_1,   \
                             ta2_yz_z_xyy_0,  \
                             ta2_yz_z_xyy_1,  \
                             ta2_yz_z_xyz_0,  \
                             ta2_yz_z_xyz_1,  \
                             ta2_yz_z_xz_0,   \
                             ta2_yz_z_xz_1,   \
                             ta2_yz_z_xzz_0,  \
                             ta2_yz_z_xzz_1,  \
                             ta2_yz_z_yy_0,   \
                             ta2_yz_z_yy_1,   \
                             ta2_yz_z_yyy_0,  \
                             ta2_yz_z_yyy_1,  \
                             ta2_yz_z_yyz_0,  \
                             ta2_yz_z_yyz_1,  \
                             ta2_yz_z_yz_0,   \
                             ta2_yz_z_yz_1,   \
                             ta2_yz_z_yzz_0,  \
                             ta2_yz_z_yzz_1,  \
                             ta2_yz_z_zz_0,   \
                             ta2_yz_z_zz_1,   \
                             ta2_yz_z_zzz_0,  \
                             ta2_yz_z_zzz_1,  \
                             ta2_yz_zz_xxx_0, \
                             ta2_yz_zz_xxy_0, \
                             ta2_yz_zz_xxz_0, \
                             ta2_yz_zz_xyy_0, \
                             ta2_yz_zz_xyz_0, \
                             ta2_yz_zz_xzz_0, \
                             ta2_yz_zz_yyy_0, \
                             ta2_yz_zz_yyz_0, \
                             ta2_yz_zz_yzz_0, \
                             ta2_yz_zz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_zz_xxx_0[i] =
            ta2_yz_0_xxx_0[i] * fe_0 - ta2_yz_0_xxx_1[i] * fe_0 + ta1_y_z_xxx_1[i] + ta2_yz_z_xxx_0[i] * pa_z[i] - ta2_yz_z_xxx_1[i] * pc_z[i];

        ta2_yz_zz_xxy_0[i] =
            ta2_yz_0_xxy_0[i] * fe_0 - ta2_yz_0_xxy_1[i] * fe_0 + ta1_y_z_xxy_1[i] + ta2_yz_z_xxy_0[i] * pa_z[i] - ta2_yz_z_xxy_1[i] * pc_z[i];

        ta2_yz_zz_xxz_0[i] = ta2_yz_0_xxz_0[i] * fe_0 - ta2_yz_0_xxz_1[i] * fe_0 + ta2_yz_z_xx_0[i] * fe_0 - ta2_yz_z_xx_1[i] * fe_0 +
                             ta1_y_z_xxz_1[i] + ta2_yz_z_xxz_0[i] * pa_z[i] - ta2_yz_z_xxz_1[i] * pc_z[i];

        ta2_yz_zz_xyy_0[i] =
            ta2_yz_0_xyy_0[i] * fe_0 - ta2_yz_0_xyy_1[i] * fe_0 + ta1_y_z_xyy_1[i] + ta2_yz_z_xyy_0[i] * pa_z[i] - ta2_yz_z_xyy_1[i] * pc_z[i];

        ta2_yz_zz_xyz_0[i] = ta2_yz_0_xyz_0[i] * fe_0 - ta2_yz_0_xyz_1[i] * fe_0 + ta2_yz_z_xy_0[i] * fe_0 - ta2_yz_z_xy_1[i] * fe_0 +
                             ta1_y_z_xyz_1[i] + ta2_yz_z_xyz_0[i] * pa_z[i] - ta2_yz_z_xyz_1[i] * pc_z[i];

        ta2_yz_zz_xzz_0[i] = ta2_yz_0_xzz_0[i] * fe_0 - ta2_yz_0_xzz_1[i] * fe_0 + 2.0 * ta2_yz_z_xz_0[i] * fe_0 - 2.0 * ta2_yz_z_xz_1[i] * fe_0 +
                             ta1_y_z_xzz_1[i] + ta2_yz_z_xzz_0[i] * pa_z[i] - ta2_yz_z_xzz_1[i] * pc_z[i];

        ta2_yz_zz_yyy_0[i] =
            ta2_yz_0_yyy_0[i] * fe_0 - ta2_yz_0_yyy_1[i] * fe_0 + ta1_y_z_yyy_1[i] + ta2_yz_z_yyy_0[i] * pa_z[i] - ta2_yz_z_yyy_1[i] * pc_z[i];

        ta2_yz_zz_yyz_0[i] = ta2_yz_0_yyz_0[i] * fe_0 - ta2_yz_0_yyz_1[i] * fe_0 + ta2_yz_z_yy_0[i] * fe_0 - ta2_yz_z_yy_1[i] * fe_0 +
                             ta1_y_z_yyz_1[i] + ta2_yz_z_yyz_0[i] * pa_z[i] - ta2_yz_z_yyz_1[i] * pc_z[i];

        ta2_yz_zz_yzz_0[i] = ta2_yz_0_yzz_0[i] * fe_0 - ta2_yz_0_yzz_1[i] * fe_0 + 2.0 * ta2_yz_z_yz_0[i] * fe_0 - 2.0 * ta2_yz_z_yz_1[i] * fe_0 +
                             ta1_y_z_yzz_1[i] + ta2_yz_z_yzz_0[i] * pa_z[i] - ta2_yz_z_yzz_1[i] * pc_z[i];

        ta2_yz_zz_zzz_0[i] = ta2_yz_0_zzz_0[i] * fe_0 - ta2_yz_0_zzz_1[i] * fe_0 + 3.0 * ta2_yz_z_zz_0[i] * fe_0 - 3.0 * ta2_yz_z_zz_1[i] * fe_0 +
                             ta1_y_z_zzz_1[i] + ta2_yz_z_zzz_0[i] * pa_z[i] - ta2_yz_z_zzz_1[i] * pc_z[i];
    }

    // Set up 300-310 components of targeted buffer : DF

    auto ta2_zz_xx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 300);

    auto ta2_zz_xx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 301);

    auto ta2_zz_xx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 302);

    auto ta2_zz_xx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 303);

    auto ta2_zz_xx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 304);

    auto ta2_zz_xx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 305);

    auto ta2_zz_xx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 306);

    auto ta2_zz_xx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 307);

    auto ta2_zz_xx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 308);

    auto ta2_zz_xx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 309);

#pragma omp simd aligned(pa_x,                \
                             pc_x,            \
                             ta2_zz_0_xxx_0,  \
                             ta2_zz_0_xxx_1,  \
                             ta2_zz_0_xxy_0,  \
                             ta2_zz_0_xxy_1,  \
                             ta2_zz_0_xxz_0,  \
                             ta2_zz_0_xxz_1,  \
                             ta2_zz_0_xyy_0,  \
                             ta2_zz_0_xyy_1,  \
                             ta2_zz_0_xyz_0,  \
                             ta2_zz_0_xyz_1,  \
                             ta2_zz_0_xzz_0,  \
                             ta2_zz_0_xzz_1,  \
                             ta2_zz_0_yyy_0,  \
                             ta2_zz_0_yyy_1,  \
                             ta2_zz_0_yyz_0,  \
                             ta2_zz_0_yyz_1,  \
                             ta2_zz_0_yzz_0,  \
                             ta2_zz_0_yzz_1,  \
                             ta2_zz_0_zzz_0,  \
                             ta2_zz_0_zzz_1,  \
                             ta2_zz_x_xx_0,   \
                             ta2_zz_x_xx_1,   \
                             ta2_zz_x_xxx_0,  \
                             ta2_zz_x_xxx_1,  \
                             ta2_zz_x_xxy_0,  \
                             ta2_zz_x_xxy_1,  \
                             ta2_zz_x_xxz_0,  \
                             ta2_zz_x_xxz_1,  \
                             ta2_zz_x_xy_0,   \
                             ta2_zz_x_xy_1,   \
                             ta2_zz_x_xyy_0,  \
                             ta2_zz_x_xyy_1,  \
                             ta2_zz_x_xyz_0,  \
                             ta2_zz_x_xyz_1,  \
                             ta2_zz_x_xz_0,   \
                             ta2_zz_x_xz_1,   \
                             ta2_zz_x_xzz_0,  \
                             ta2_zz_x_xzz_1,  \
                             ta2_zz_x_yy_0,   \
                             ta2_zz_x_yy_1,   \
                             ta2_zz_x_yyy_0,  \
                             ta2_zz_x_yyy_1,  \
                             ta2_zz_x_yyz_0,  \
                             ta2_zz_x_yyz_1,  \
                             ta2_zz_x_yz_0,   \
                             ta2_zz_x_yz_1,   \
                             ta2_zz_x_yzz_0,  \
                             ta2_zz_x_yzz_1,  \
                             ta2_zz_x_zz_0,   \
                             ta2_zz_x_zz_1,   \
                             ta2_zz_x_zzz_0,  \
                             ta2_zz_x_zzz_1,  \
                             ta2_zz_xx_xxx_0, \
                             ta2_zz_xx_xxy_0, \
                             ta2_zz_xx_xxz_0, \
                             ta2_zz_xx_xyy_0, \
                             ta2_zz_xx_xyz_0, \
                             ta2_zz_xx_xzz_0, \
                             ta2_zz_xx_yyy_0, \
                             ta2_zz_xx_yyz_0, \
                             ta2_zz_xx_yzz_0, \
                             ta2_zz_xx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xx_xxx_0[i] = ta2_zz_0_xxx_0[i] * fe_0 - ta2_zz_0_xxx_1[i] * fe_0 + 3.0 * ta2_zz_x_xx_0[i] * fe_0 - 3.0 * ta2_zz_x_xx_1[i] * fe_0 +
                             ta2_zz_x_xxx_0[i] * pa_x[i] - ta2_zz_x_xxx_1[i] * pc_x[i];

        ta2_zz_xx_xxy_0[i] = ta2_zz_0_xxy_0[i] * fe_0 - ta2_zz_0_xxy_1[i] * fe_0 + 2.0 * ta2_zz_x_xy_0[i] * fe_0 - 2.0 * ta2_zz_x_xy_1[i] * fe_0 +
                             ta2_zz_x_xxy_0[i] * pa_x[i] - ta2_zz_x_xxy_1[i] * pc_x[i];

        ta2_zz_xx_xxz_0[i] = ta2_zz_0_xxz_0[i] * fe_0 - ta2_zz_0_xxz_1[i] * fe_0 + 2.0 * ta2_zz_x_xz_0[i] * fe_0 - 2.0 * ta2_zz_x_xz_1[i] * fe_0 +
                             ta2_zz_x_xxz_0[i] * pa_x[i] - ta2_zz_x_xxz_1[i] * pc_x[i];

        ta2_zz_xx_xyy_0[i] = ta2_zz_0_xyy_0[i] * fe_0 - ta2_zz_0_xyy_1[i] * fe_0 + ta2_zz_x_yy_0[i] * fe_0 - ta2_zz_x_yy_1[i] * fe_0 +
                             ta2_zz_x_xyy_0[i] * pa_x[i] - ta2_zz_x_xyy_1[i] * pc_x[i];

        ta2_zz_xx_xyz_0[i] = ta2_zz_0_xyz_0[i] * fe_0 - ta2_zz_0_xyz_1[i] * fe_0 + ta2_zz_x_yz_0[i] * fe_0 - ta2_zz_x_yz_1[i] * fe_0 +
                             ta2_zz_x_xyz_0[i] * pa_x[i] - ta2_zz_x_xyz_1[i] * pc_x[i];

        ta2_zz_xx_xzz_0[i] = ta2_zz_0_xzz_0[i] * fe_0 - ta2_zz_0_xzz_1[i] * fe_0 + ta2_zz_x_zz_0[i] * fe_0 - ta2_zz_x_zz_1[i] * fe_0 +
                             ta2_zz_x_xzz_0[i] * pa_x[i] - ta2_zz_x_xzz_1[i] * pc_x[i];

        ta2_zz_xx_yyy_0[i] = ta2_zz_0_yyy_0[i] * fe_0 - ta2_zz_0_yyy_1[i] * fe_0 + ta2_zz_x_yyy_0[i] * pa_x[i] - ta2_zz_x_yyy_1[i] * pc_x[i];

        ta2_zz_xx_yyz_0[i] = ta2_zz_0_yyz_0[i] * fe_0 - ta2_zz_0_yyz_1[i] * fe_0 + ta2_zz_x_yyz_0[i] * pa_x[i] - ta2_zz_x_yyz_1[i] * pc_x[i];

        ta2_zz_xx_yzz_0[i] = ta2_zz_0_yzz_0[i] * fe_0 - ta2_zz_0_yzz_1[i] * fe_0 + ta2_zz_x_yzz_0[i] * pa_x[i] - ta2_zz_x_yzz_1[i] * pc_x[i];

        ta2_zz_xx_zzz_0[i] = ta2_zz_0_zzz_0[i] * fe_0 - ta2_zz_0_zzz_1[i] * fe_0 + ta2_zz_x_zzz_0[i] * pa_x[i] - ta2_zz_x_zzz_1[i] * pc_x[i];
    }

    // Set up 310-320 components of targeted buffer : DF

    auto ta2_zz_xy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 310);

    auto ta2_zz_xy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 311);

    auto ta2_zz_xy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 312);

    auto ta2_zz_xy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 313);

    auto ta2_zz_xy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 314);

    auto ta2_zz_xy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 315);

    auto ta2_zz_xy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 316);

    auto ta2_zz_xy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 317);

    auto ta2_zz_xy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 318);

    auto ta2_zz_xy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 319);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pc_x,            \
                             pc_y,            \
                             ta2_zz_x_xxx_0,  \
                             ta2_zz_x_xxx_1,  \
                             ta2_zz_x_xxz_0,  \
                             ta2_zz_x_xxz_1,  \
                             ta2_zz_x_xzz_0,  \
                             ta2_zz_x_xzz_1,  \
                             ta2_zz_xy_xxx_0, \
                             ta2_zz_xy_xxy_0, \
                             ta2_zz_xy_xxz_0, \
                             ta2_zz_xy_xyy_0, \
                             ta2_zz_xy_xyz_0, \
                             ta2_zz_xy_xzz_0, \
                             ta2_zz_xy_yyy_0, \
                             ta2_zz_xy_yyz_0, \
                             ta2_zz_xy_yzz_0, \
                             ta2_zz_xy_zzz_0, \
                             ta2_zz_y_xxy_0,  \
                             ta2_zz_y_xxy_1,  \
                             ta2_zz_y_xy_0,   \
                             ta2_zz_y_xy_1,   \
                             ta2_zz_y_xyy_0,  \
                             ta2_zz_y_xyy_1,  \
                             ta2_zz_y_xyz_0,  \
                             ta2_zz_y_xyz_1,  \
                             ta2_zz_y_yy_0,   \
                             ta2_zz_y_yy_1,   \
                             ta2_zz_y_yyy_0,  \
                             ta2_zz_y_yyy_1,  \
                             ta2_zz_y_yyz_0,  \
                             ta2_zz_y_yyz_1,  \
                             ta2_zz_y_yz_0,   \
                             ta2_zz_y_yz_1,   \
                             ta2_zz_y_yzz_0,  \
                             ta2_zz_y_yzz_1,  \
                             ta2_zz_y_zzz_0,  \
                             ta2_zz_y_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xy_xxx_0[i] = ta2_zz_x_xxx_0[i] * pa_y[i] - ta2_zz_x_xxx_1[i] * pc_y[i];

        ta2_zz_xy_xxy_0[i] =
            2.0 * ta2_zz_y_xy_0[i] * fe_0 - 2.0 * ta2_zz_y_xy_1[i] * fe_0 + ta2_zz_y_xxy_0[i] * pa_x[i] - ta2_zz_y_xxy_1[i] * pc_x[i];

        ta2_zz_xy_xxz_0[i] = ta2_zz_x_xxz_0[i] * pa_y[i] - ta2_zz_x_xxz_1[i] * pc_y[i];

        ta2_zz_xy_xyy_0[i] = ta2_zz_y_yy_0[i] * fe_0 - ta2_zz_y_yy_1[i] * fe_0 + ta2_zz_y_xyy_0[i] * pa_x[i] - ta2_zz_y_xyy_1[i] * pc_x[i];

        ta2_zz_xy_xyz_0[i] = ta2_zz_y_yz_0[i] * fe_0 - ta2_zz_y_yz_1[i] * fe_0 + ta2_zz_y_xyz_0[i] * pa_x[i] - ta2_zz_y_xyz_1[i] * pc_x[i];

        ta2_zz_xy_xzz_0[i] = ta2_zz_x_xzz_0[i] * pa_y[i] - ta2_zz_x_xzz_1[i] * pc_y[i];

        ta2_zz_xy_yyy_0[i] = ta2_zz_y_yyy_0[i] * pa_x[i] - ta2_zz_y_yyy_1[i] * pc_x[i];

        ta2_zz_xy_yyz_0[i] = ta2_zz_y_yyz_0[i] * pa_x[i] - ta2_zz_y_yyz_1[i] * pc_x[i];

        ta2_zz_xy_yzz_0[i] = ta2_zz_y_yzz_0[i] * pa_x[i] - ta2_zz_y_yzz_1[i] * pc_x[i];

        ta2_zz_xy_zzz_0[i] = ta2_zz_y_zzz_0[i] * pa_x[i] - ta2_zz_y_zzz_1[i] * pc_x[i];
    }

    // Set up 320-330 components of targeted buffer : DF

    auto ta2_zz_xz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 320);

    auto ta2_zz_xz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 321);

    auto ta2_zz_xz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 322);

    auto ta2_zz_xz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 323);

    auto ta2_zz_xz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 324);

    auto ta2_zz_xz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 325);

    auto ta2_zz_xz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 326);

    auto ta2_zz_xz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 327);

    auto ta2_zz_xz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 328);

    auto ta2_zz_xz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 329);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             pc_x,            \
                             pc_z,            \
                             ta1_z_x_xxx_1,   \
                             ta1_z_x_xxy_1,   \
                             ta1_z_x_xyy_1,   \
                             ta2_zz_x_xxx_0,  \
                             ta2_zz_x_xxx_1,  \
                             ta2_zz_x_xxy_0,  \
                             ta2_zz_x_xxy_1,  \
                             ta2_zz_x_xyy_0,  \
                             ta2_zz_x_xyy_1,  \
                             ta2_zz_xz_xxx_0, \
                             ta2_zz_xz_xxy_0, \
                             ta2_zz_xz_xxz_0, \
                             ta2_zz_xz_xyy_0, \
                             ta2_zz_xz_xyz_0, \
                             ta2_zz_xz_xzz_0, \
                             ta2_zz_xz_yyy_0, \
                             ta2_zz_xz_yyz_0, \
                             ta2_zz_xz_yzz_0, \
                             ta2_zz_xz_zzz_0, \
                             ta2_zz_z_xxz_0,  \
                             ta2_zz_z_xxz_1,  \
                             ta2_zz_z_xyz_0,  \
                             ta2_zz_z_xyz_1,  \
                             ta2_zz_z_xz_0,   \
                             ta2_zz_z_xz_1,   \
                             ta2_zz_z_xzz_0,  \
                             ta2_zz_z_xzz_1,  \
                             ta2_zz_z_yyy_0,  \
                             ta2_zz_z_yyy_1,  \
                             ta2_zz_z_yyz_0,  \
                             ta2_zz_z_yyz_1,  \
                             ta2_zz_z_yz_0,   \
                             ta2_zz_z_yz_1,   \
                             ta2_zz_z_yzz_0,  \
                             ta2_zz_z_yzz_1,  \
                             ta2_zz_z_zz_0,   \
                             ta2_zz_z_zz_1,   \
                             ta2_zz_z_zzz_0,  \
                             ta2_zz_z_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xz_xxx_0[i] = 2.0 * ta1_z_x_xxx_1[i] + ta2_zz_x_xxx_0[i] * pa_z[i] - ta2_zz_x_xxx_1[i] * pc_z[i];

        ta2_zz_xz_xxy_0[i] = 2.0 * ta1_z_x_xxy_1[i] + ta2_zz_x_xxy_0[i] * pa_z[i] - ta2_zz_x_xxy_1[i] * pc_z[i];

        ta2_zz_xz_xxz_0[i] =
            2.0 * ta2_zz_z_xz_0[i] * fe_0 - 2.0 * ta2_zz_z_xz_1[i] * fe_0 + ta2_zz_z_xxz_0[i] * pa_x[i] - ta2_zz_z_xxz_1[i] * pc_x[i];

        ta2_zz_xz_xyy_0[i] = 2.0 * ta1_z_x_xyy_1[i] + ta2_zz_x_xyy_0[i] * pa_z[i] - ta2_zz_x_xyy_1[i] * pc_z[i];

        ta2_zz_xz_xyz_0[i] = ta2_zz_z_yz_0[i] * fe_0 - ta2_zz_z_yz_1[i] * fe_0 + ta2_zz_z_xyz_0[i] * pa_x[i] - ta2_zz_z_xyz_1[i] * pc_x[i];

        ta2_zz_xz_xzz_0[i] = ta2_zz_z_zz_0[i] * fe_0 - ta2_zz_z_zz_1[i] * fe_0 + ta2_zz_z_xzz_0[i] * pa_x[i] - ta2_zz_z_xzz_1[i] * pc_x[i];

        ta2_zz_xz_yyy_0[i] = ta2_zz_z_yyy_0[i] * pa_x[i] - ta2_zz_z_yyy_1[i] * pc_x[i];

        ta2_zz_xz_yyz_0[i] = ta2_zz_z_yyz_0[i] * pa_x[i] - ta2_zz_z_yyz_1[i] * pc_x[i];

        ta2_zz_xz_yzz_0[i] = ta2_zz_z_yzz_0[i] * pa_x[i] - ta2_zz_z_yzz_1[i] * pc_x[i];

        ta2_zz_xz_zzz_0[i] = ta2_zz_z_zzz_0[i] * pa_x[i] - ta2_zz_z_zzz_1[i] * pc_x[i];
    }

    // Set up 330-340 components of targeted buffer : DF

    auto ta2_zz_yy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 330);

    auto ta2_zz_yy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 331);

    auto ta2_zz_yy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 332);

    auto ta2_zz_yy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 333);

    auto ta2_zz_yy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 334);

    auto ta2_zz_yy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 335);

    auto ta2_zz_yy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 336);

    auto ta2_zz_yy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 337);

    auto ta2_zz_yy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 338);

    auto ta2_zz_yy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 339);

#pragma omp simd aligned(pa_y,                \
                             pc_y,            \
                             ta2_zz_0_xxx_0,  \
                             ta2_zz_0_xxx_1,  \
                             ta2_zz_0_xxy_0,  \
                             ta2_zz_0_xxy_1,  \
                             ta2_zz_0_xxz_0,  \
                             ta2_zz_0_xxz_1,  \
                             ta2_zz_0_xyy_0,  \
                             ta2_zz_0_xyy_1,  \
                             ta2_zz_0_xyz_0,  \
                             ta2_zz_0_xyz_1,  \
                             ta2_zz_0_xzz_0,  \
                             ta2_zz_0_xzz_1,  \
                             ta2_zz_0_yyy_0,  \
                             ta2_zz_0_yyy_1,  \
                             ta2_zz_0_yyz_0,  \
                             ta2_zz_0_yyz_1,  \
                             ta2_zz_0_yzz_0,  \
                             ta2_zz_0_yzz_1,  \
                             ta2_zz_0_zzz_0,  \
                             ta2_zz_0_zzz_1,  \
                             ta2_zz_y_xx_0,   \
                             ta2_zz_y_xx_1,   \
                             ta2_zz_y_xxx_0,  \
                             ta2_zz_y_xxx_1,  \
                             ta2_zz_y_xxy_0,  \
                             ta2_zz_y_xxy_1,  \
                             ta2_zz_y_xxz_0,  \
                             ta2_zz_y_xxz_1,  \
                             ta2_zz_y_xy_0,   \
                             ta2_zz_y_xy_1,   \
                             ta2_zz_y_xyy_0,  \
                             ta2_zz_y_xyy_1,  \
                             ta2_zz_y_xyz_0,  \
                             ta2_zz_y_xyz_1,  \
                             ta2_zz_y_xz_0,   \
                             ta2_zz_y_xz_1,   \
                             ta2_zz_y_xzz_0,  \
                             ta2_zz_y_xzz_1,  \
                             ta2_zz_y_yy_0,   \
                             ta2_zz_y_yy_1,   \
                             ta2_zz_y_yyy_0,  \
                             ta2_zz_y_yyy_1,  \
                             ta2_zz_y_yyz_0,  \
                             ta2_zz_y_yyz_1,  \
                             ta2_zz_y_yz_0,   \
                             ta2_zz_y_yz_1,   \
                             ta2_zz_y_yzz_0,  \
                             ta2_zz_y_yzz_1,  \
                             ta2_zz_y_zz_0,   \
                             ta2_zz_y_zz_1,   \
                             ta2_zz_y_zzz_0,  \
                             ta2_zz_y_zzz_1,  \
                             ta2_zz_yy_xxx_0, \
                             ta2_zz_yy_xxy_0, \
                             ta2_zz_yy_xxz_0, \
                             ta2_zz_yy_xyy_0, \
                             ta2_zz_yy_xyz_0, \
                             ta2_zz_yy_xzz_0, \
                             ta2_zz_yy_yyy_0, \
                             ta2_zz_yy_yyz_0, \
                             ta2_zz_yy_yzz_0, \
                             ta2_zz_yy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yy_xxx_0[i] = ta2_zz_0_xxx_0[i] * fe_0 - ta2_zz_0_xxx_1[i] * fe_0 + ta2_zz_y_xxx_0[i] * pa_y[i] - ta2_zz_y_xxx_1[i] * pc_y[i];

        ta2_zz_yy_xxy_0[i] = ta2_zz_0_xxy_0[i] * fe_0 - ta2_zz_0_xxy_1[i] * fe_0 + ta2_zz_y_xx_0[i] * fe_0 - ta2_zz_y_xx_1[i] * fe_0 +
                             ta2_zz_y_xxy_0[i] * pa_y[i] - ta2_zz_y_xxy_1[i] * pc_y[i];

        ta2_zz_yy_xxz_0[i] = ta2_zz_0_xxz_0[i] * fe_0 - ta2_zz_0_xxz_1[i] * fe_0 + ta2_zz_y_xxz_0[i] * pa_y[i] - ta2_zz_y_xxz_1[i] * pc_y[i];

        ta2_zz_yy_xyy_0[i] = ta2_zz_0_xyy_0[i] * fe_0 - ta2_zz_0_xyy_1[i] * fe_0 + 2.0 * ta2_zz_y_xy_0[i] * fe_0 - 2.0 * ta2_zz_y_xy_1[i] * fe_0 +
                             ta2_zz_y_xyy_0[i] * pa_y[i] - ta2_zz_y_xyy_1[i] * pc_y[i];

        ta2_zz_yy_xyz_0[i] = ta2_zz_0_xyz_0[i] * fe_0 - ta2_zz_0_xyz_1[i] * fe_0 + ta2_zz_y_xz_0[i] * fe_0 - ta2_zz_y_xz_1[i] * fe_0 +
                             ta2_zz_y_xyz_0[i] * pa_y[i] - ta2_zz_y_xyz_1[i] * pc_y[i];

        ta2_zz_yy_xzz_0[i] = ta2_zz_0_xzz_0[i] * fe_0 - ta2_zz_0_xzz_1[i] * fe_0 + ta2_zz_y_xzz_0[i] * pa_y[i] - ta2_zz_y_xzz_1[i] * pc_y[i];

        ta2_zz_yy_yyy_0[i] = ta2_zz_0_yyy_0[i] * fe_0 - ta2_zz_0_yyy_1[i] * fe_0 + 3.0 * ta2_zz_y_yy_0[i] * fe_0 - 3.0 * ta2_zz_y_yy_1[i] * fe_0 +
                             ta2_zz_y_yyy_0[i] * pa_y[i] - ta2_zz_y_yyy_1[i] * pc_y[i];

        ta2_zz_yy_yyz_0[i] = ta2_zz_0_yyz_0[i] * fe_0 - ta2_zz_0_yyz_1[i] * fe_0 + 2.0 * ta2_zz_y_yz_0[i] * fe_0 - 2.0 * ta2_zz_y_yz_1[i] * fe_0 +
                             ta2_zz_y_yyz_0[i] * pa_y[i] - ta2_zz_y_yyz_1[i] * pc_y[i];

        ta2_zz_yy_yzz_0[i] = ta2_zz_0_yzz_0[i] * fe_0 - ta2_zz_0_yzz_1[i] * fe_0 + ta2_zz_y_zz_0[i] * fe_0 - ta2_zz_y_zz_1[i] * fe_0 +
                             ta2_zz_y_yzz_0[i] * pa_y[i] - ta2_zz_y_yzz_1[i] * pc_y[i];

        ta2_zz_yy_zzz_0[i] = ta2_zz_0_zzz_0[i] * fe_0 - ta2_zz_0_zzz_1[i] * fe_0 + ta2_zz_y_zzz_0[i] * pa_y[i] - ta2_zz_y_zzz_1[i] * pc_y[i];
    }

    // Set up 340-350 components of targeted buffer : DF

    auto ta2_zz_yz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 340);

    auto ta2_zz_yz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 341);

    auto ta2_zz_yz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 342);

    auto ta2_zz_yz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 343);

    auto ta2_zz_yz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 344);

    auto ta2_zz_yz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 345);

    auto ta2_zz_yz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 346);

    auto ta2_zz_yz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 347);

    auto ta2_zz_yz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 348);

    auto ta2_zz_yz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 349);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             pc_y,            \
                             pc_z,            \
                             ta1_z_y_xxy_1,   \
                             ta1_z_y_xyy_1,   \
                             ta1_z_y_yyy_1,   \
                             ta2_zz_y_xxy_0,  \
                             ta2_zz_y_xxy_1,  \
                             ta2_zz_y_xyy_0,  \
                             ta2_zz_y_xyy_1,  \
                             ta2_zz_y_yyy_0,  \
                             ta2_zz_y_yyy_1,  \
                             ta2_zz_yz_xxx_0, \
                             ta2_zz_yz_xxy_0, \
                             ta2_zz_yz_xxz_0, \
                             ta2_zz_yz_xyy_0, \
                             ta2_zz_yz_xyz_0, \
                             ta2_zz_yz_xzz_0, \
                             ta2_zz_yz_yyy_0, \
                             ta2_zz_yz_yyz_0, \
                             ta2_zz_yz_yzz_0, \
                             ta2_zz_yz_zzz_0, \
                             ta2_zz_z_xxx_0,  \
                             ta2_zz_z_xxx_1,  \
                             ta2_zz_z_xxz_0,  \
                             ta2_zz_z_xxz_1,  \
                             ta2_zz_z_xyz_0,  \
                             ta2_zz_z_xyz_1,  \
                             ta2_zz_z_xz_0,   \
                             ta2_zz_z_xz_1,   \
                             ta2_zz_z_xzz_0,  \
                             ta2_zz_z_xzz_1,  \
                             ta2_zz_z_yyz_0,  \
                             ta2_zz_z_yyz_1,  \
                             ta2_zz_z_yz_0,   \
                             ta2_zz_z_yz_1,   \
                             ta2_zz_z_yzz_0,  \
                             ta2_zz_z_yzz_1,  \
                             ta2_zz_z_zz_0,   \
                             ta2_zz_z_zz_1,   \
                             ta2_zz_z_zzz_0,  \
                             ta2_zz_z_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yz_xxx_0[i] = ta2_zz_z_xxx_0[i] * pa_y[i] - ta2_zz_z_xxx_1[i] * pc_y[i];

        ta2_zz_yz_xxy_0[i] = 2.0 * ta1_z_y_xxy_1[i] + ta2_zz_y_xxy_0[i] * pa_z[i] - ta2_zz_y_xxy_1[i] * pc_z[i];

        ta2_zz_yz_xxz_0[i] = ta2_zz_z_xxz_0[i] * pa_y[i] - ta2_zz_z_xxz_1[i] * pc_y[i];

        ta2_zz_yz_xyy_0[i] = 2.0 * ta1_z_y_xyy_1[i] + ta2_zz_y_xyy_0[i] * pa_z[i] - ta2_zz_y_xyy_1[i] * pc_z[i];

        ta2_zz_yz_xyz_0[i] = ta2_zz_z_xz_0[i] * fe_0 - ta2_zz_z_xz_1[i] * fe_0 + ta2_zz_z_xyz_0[i] * pa_y[i] - ta2_zz_z_xyz_1[i] * pc_y[i];

        ta2_zz_yz_xzz_0[i] = ta2_zz_z_xzz_0[i] * pa_y[i] - ta2_zz_z_xzz_1[i] * pc_y[i];

        ta2_zz_yz_yyy_0[i] = 2.0 * ta1_z_y_yyy_1[i] + ta2_zz_y_yyy_0[i] * pa_z[i] - ta2_zz_y_yyy_1[i] * pc_z[i];

        ta2_zz_yz_yyz_0[i] =
            2.0 * ta2_zz_z_yz_0[i] * fe_0 - 2.0 * ta2_zz_z_yz_1[i] * fe_0 + ta2_zz_z_yyz_0[i] * pa_y[i] - ta2_zz_z_yyz_1[i] * pc_y[i];

        ta2_zz_yz_yzz_0[i] = ta2_zz_z_zz_0[i] * fe_0 - ta2_zz_z_zz_1[i] * fe_0 + ta2_zz_z_yzz_0[i] * pa_y[i] - ta2_zz_z_yzz_1[i] * pc_y[i];

        ta2_zz_yz_zzz_0[i] = ta2_zz_z_zzz_0[i] * pa_y[i] - ta2_zz_z_zzz_1[i] * pc_y[i];
    }

    // Set up 350-360 components of targeted buffer : DF

    auto ta2_zz_zz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_df + 350);

    auto ta2_zz_zz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_df + 351);

    auto ta2_zz_zz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_df + 352);

    auto ta2_zz_zz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 353);

    auto ta2_zz_zz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 354);

    auto ta2_zz_zz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 355);

    auto ta2_zz_zz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_df + 356);

    auto ta2_zz_zz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_df + 357);

    auto ta2_zz_zz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 358);

    auto ta2_zz_zz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_df + 359);

#pragma omp simd aligned(pa_z,                \
                             pc_z,            \
                             ta1_z_z_xxx_1,   \
                             ta1_z_z_xxy_1,   \
                             ta1_z_z_xxz_1,   \
                             ta1_z_z_xyy_1,   \
                             ta1_z_z_xyz_1,   \
                             ta1_z_z_xzz_1,   \
                             ta1_z_z_yyy_1,   \
                             ta1_z_z_yyz_1,   \
                             ta1_z_z_yzz_1,   \
                             ta1_z_z_zzz_1,   \
                             ta2_zz_0_xxx_0,  \
                             ta2_zz_0_xxx_1,  \
                             ta2_zz_0_xxy_0,  \
                             ta2_zz_0_xxy_1,  \
                             ta2_zz_0_xxz_0,  \
                             ta2_zz_0_xxz_1,  \
                             ta2_zz_0_xyy_0,  \
                             ta2_zz_0_xyy_1,  \
                             ta2_zz_0_xyz_0,  \
                             ta2_zz_0_xyz_1,  \
                             ta2_zz_0_xzz_0,  \
                             ta2_zz_0_xzz_1,  \
                             ta2_zz_0_yyy_0,  \
                             ta2_zz_0_yyy_1,  \
                             ta2_zz_0_yyz_0,  \
                             ta2_zz_0_yyz_1,  \
                             ta2_zz_0_yzz_0,  \
                             ta2_zz_0_yzz_1,  \
                             ta2_zz_0_zzz_0,  \
                             ta2_zz_0_zzz_1,  \
                             ta2_zz_z_xx_0,   \
                             ta2_zz_z_xx_1,   \
                             ta2_zz_z_xxx_0,  \
                             ta2_zz_z_xxx_1,  \
                             ta2_zz_z_xxy_0,  \
                             ta2_zz_z_xxy_1,  \
                             ta2_zz_z_xxz_0,  \
                             ta2_zz_z_xxz_1,  \
                             ta2_zz_z_xy_0,   \
                             ta2_zz_z_xy_1,   \
                             ta2_zz_z_xyy_0,  \
                             ta2_zz_z_xyy_1,  \
                             ta2_zz_z_xyz_0,  \
                             ta2_zz_z_xyz_1,  \
                             ta2_zz_z_xz_0,   \
                             ta2_zz_z_xz_1,   \
                             ta2_zz_z_xzz_0,  \
                             ta2_zz_z_xzz_1,  \
                             ta2_zz_z_yy_0,   \
                             ta2_zz_z_yy_1,   \
                             ta2_zz_z_yyy_0,  \
                             ta2_zz_z_yyy_1,  \
                             ta2_zz_z_yyz_0,  \
                             ta2_zz_z_yyz_1,  \
                             ta2_zz_z_yz_0,   \
                             ta2_zz_z_yz_1,   \
                             ta2_zz_z_yzz_0,  \
                             ta2_zz_z_yzz_1,  \
                             ta2_zz_z_zz_0,   \
                             ta2_zz_z_zz_1,   \
                             ta2_zz_z_zzz_0,  \
                             ta2_zz_z_zzz_1,  \
                             ta2_zz_zz_xxx_0, \
                             ta2_zz_zz_xxy_0, \
                             ta2_zz_zz_xxz_0, \
                             ta2_zz_zz_xyy_0, \
                             ta2_zz_zz_xyz_0, \
                             ta2_zz_zz_xzz_0, \
                             ta2_zz_zz_yyy_0, \
                             ta2_zz_zz_yyz_0, \
                             ta2_zz_zz_yzz_0, \
                             ta2_zz_zz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_zz_xxx_0[i] =
            ta2_zz_0_xxx_0[i] * fe_0 - ta2_zz_0_xxx_1[i] * fe_0 + 2.0 * ta1_z_z_xxx_1[i] + ta2_zz_z_xxx_0[i] * pa_z[i] - ta2_zz_z_xxx_1[i] * pc_z[i];

        ta2_zz_zz_xxy_0[i] =
            ta2_zz_0_xxy_0[i] * fe_0 - ta2_zz_0_xxy_1[i] * fe_0 + 2.0 * ta1_z_z_xxy_1[i] + ta2_zz_z_xxy_0[i] * pa_z[i] - ta2_zz_z_xxy_1[i] * pc_z[i];

        ta2_zz_zz_xxz_0[i] = ta2_zz_0_xxz_0[i] * fe_0 - ta2_zz_0_xxz_1[i] * fe_0 + ta2_zz_z_xx_0[i] * fe_0 - ta2_zz_z_xx_1[i] * fe_0 +
                             2.0 * ta1_z_z_xxz_1[i] + ta2_zz_z_xxz_0[i] * pa_z[i] - ta2_zz_z_xxz_1[i] * pc_z[i];

        ta2_zz_zz_xyy_0[i] =
            ta2_zz_0_xyy_0[i] * fe_0 - ta2_zz_0_xyy_1[i] * fe_0 + 2.0 * ta1_z_z_xyy_1[i] + ta2_zz_z_xyy_0[i] * pa_z[i] - ta2_zz_z_xyy_1[i] * pc_z[i];

        ta2_zz_zz_xyz_0[i] = ta2_zz_0_xyz_0[i] * fe_0 - ta2_zz_0_xyz_1[i] * fe_0 + ta2_zz_z_xy_0[i] * fe_0 - ta2_zz_z_xy_1[i] * fe_0 +
                             2.0 * ta1_z_z_xyz_1[i] + ta2_zz_z_xyz_0[i] * pa_z[i] - ta2_zz_z_xyz_1[i] * pc_z[i];

        ta2_zz_zz_xzz_0[i] = ta2_zz_0_xzz_0[i] * fe_0 - ta2_zz_0_xzz_1[i] * fe_0 + 2.0 * ta2_zz_z_xz_0[i] * fe_0 - 2.0 * ta2_zz_z_xz_1[i] * fe_0 +
                             2.0 * ta1_z_z_xzz_1[i] + ta2_zz_z_xzz_0[i] * pa_z[i] - ta2_zz_z_xzz_1[i] * pc_z[i];

        ta2_zz_zz_yyy_0[i] =
            ta2_zz_0_yyy_0[i] * fe_0 - ta2_zz_0_yyy_1[i] * fe_0 + 2.0 * ta1_z_z_yyy_1[i] + ta2_zz_z_yyy_0[i] * pa_z[i] - ta2_zz_z_yyy_1[i] * pc_z[i];

        ta2_zz_zz_yyz_0[i] = ta2_zz_0_yyz_0[i] * fe_0 - ta2_zz_0_yyz_1[i] * fe_0 + ta2_zz_z_yy_0[i] * fe_0 - ta2_zz_z_yy_1[i] * fe_0 +
                             2.0 * ta1_z_z_yyz_1[i] + ta2_zz_z_yyz_0[i] * pa_z[i] - ta2_zz_z_yyz_1[i] * pc_z[i];

        ta2_zz_zz_yzz_0[i] = ta2_zz_0_yzz_0[i] * fe_0 - ta2_zz_0_yzz_1[i] * fe_0 + 2.0 * ta2_zz_z_yz_0[i] * fe_0 - 2.0 * ta2_zz_z_yz_1[i] * fe_0 +
                             2.0 * ta1_z_z_yzz_1[i] + ta2_zz_z_yzz_0[i] * pa_z[i] - ta2_zz_z_yzz_1[i] * pc_z[i];

        ta2_zz_zz_zzz_0[i] = ta2_zz_0_zzz_0[i] * fe_0 - ta2_zz_0_zzz_1[i] * fe_0 + 3.0 * ta2_zz_z_zz_0[i] * fe_0 - 3.0 * ta2_zz_z_zz_1[i] * fe_0 +
                             2.0 * ta1_z_z_zzz_1[i] + ta2_zz_z_zzz_0[i] * pa_z[i] - ta2_zz_z_zzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
