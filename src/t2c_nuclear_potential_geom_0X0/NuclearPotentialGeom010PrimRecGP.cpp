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

#include "NuclearPotentialGeom010PrimRecGP.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_gp(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_gp,
                                        const size_t              idx_npot_geom_010_0_dp,
                                        const size_t              idx_npot_geom_010_1_dp,
                                        const size_t              idx_npot_geom_010_0_fs,
                                        const size_t              idx_npot_geom_010_1_fs,
                                        const size_t              idx_npot_1_fp,
                                        const size_t              idx_npot_geom_010_0_fp,
                                        const size_t              idx_npot_geom_010_1_fp,
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

    // Set up components of auxiliary buffer : DP

    auto ta1_x_xx_x_0 = pbuffer.data(idx_npot_geom_010_0_dp);

    auto ta1_x_xx_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 1);

    auto ta1_x_xx_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 2);

    auto ta1_x_xy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 3);

    auto ta1_x_xz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 6);

    auto ta1_x_yy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 9);

    auto ta1_x_yy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 10);

    auto ta1_x_yy_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 11);

    auto ta1_x_yz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 14);

    auto ta1_x_zz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 15);

    auto ta1_x_zz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 16);

    auto ta1_x_zz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 17);

    auto ta1_y_xx_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 18);

    auto ta1_y_xx_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 19);

    auto ta1_y_xx_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 20);

    auto ta1_y_xy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 22);

    auto ta1_y_xz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 26);

    auto ta1_y_yy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 27);

    auto ta1_y_yy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 28);

    auto ta1_y_yy_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 29);

    auto ta1_y_yz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 31);

    auto ta1_y_zz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 33);

    auto ta1_y_zz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 34);

    auto ta1_y_zz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 35);

    auto ta1_z_xx_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 36);

    auto ta1_z_xx_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 37);

    auto ta1_z_xx_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 38);

    auto ta1_z_xy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 40);

    auto ta1_z_xz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 44);

    auto ta1_z_yy_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 45);

    auto ta1_z_yy_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 46);

    auto ta1_z_yy_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 47);

    auto ta1_z_yz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 50);

    auto ta1_z_zz_x_0 = pbuffer.data(idx_npot_geom_010_0_dp + 51);

    auto ta1_z_zz_y_0 = pbuffer.data(idx_npot_geom_010_0_dp + 52);

    auto ta1_z_zz_z_0 = pbuffer.data(idx_npot_geom_010_0_dp + 53);

    // Set up components of auxiliary buffer : DP

    auto ta1_x_xx_x_1 = pbuffer.data(idx_npot_geom_010_1_dp);

    auto ta1_x_xx_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 1);

    auto ta1_x_xx_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 2);

    auto ta1_x_xy_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 3);

    auto ta1_x_xz_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 6);

    auto ta1_x_yy_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 9);

    auto ta1_x_yy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 10);

    auto ta1_x_yy_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 11);

    auto ta1_x_yz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 14);

    auto ta1_x_zz_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 15);

    auto ta1_x_zz_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 16);

    auto ta1_x_zz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 17);

    auto ta1_y_xx_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 18);

    auto ta1_y_xx_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 19);

    auto ta1_y_xx_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 20);

    auto ta1_y_xy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 22);

    auto ta1_y_xz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 26);

    auto ta1_y_yy_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 27);

    auto ta1_y_yy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 28);

    auto ta1_y_yy_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 29);

    auto ta1_y_yz_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 31);

    auto ta1_y_zz_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 33);

    auto ta1_y_zz_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 34);

    auto ta1_y_zz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 35);

    auto ta1_z_xx_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 36);

    auto ta1_z_xx_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 37);

    auto ta1_z_xx_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 38);

    auto ta1_z_xy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 40);

    auto ta1_z_xz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 44);

    auto ta1_z_yy_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 45);

    auto ta1_z_yy_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 46);

    auto ta1_z_yy_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 47);

    auto ta1_z_yz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 50);

    auto ta1_z_zz_x_1 = pbuffer.data(idx_npot_geom_010_1_dp + 51);

    auto ta1_z_zz_y_1 = pbuffer.data(idx_npot_geom_010_1_dp + 52);

    auto ta1_z_zz_z_1 = pbuffer.data(idx_npot_geom_010_1_dp + 53);

    // Set up components of auxiliary buffer : FS

    auto ta1_x_xxx_0_0 = pbuffer.data(idx_npot_geom_010_0_fs);

    auto ta1_x_yyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 6);

    auto ta1_x_zzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 9);

    auto ta1_y_xxx_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 10);

    auto ta1_y_yyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 16);

    auto ta1_y_zzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 19);

    auto ta1_z_xxx_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 20);

    auto ta1_z_yyy_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 26);

    auto ta1_z_zzz_0_0 = pbuffer.data(idx_npot_geom_010_0_fs + 29);

    // Set up components of auxiliary buffer : FS

    auto ta1_x_xxx_0_1 = pbuffer.data(idx_npot_geom_010_1_fs);

    auto ta1_x_yyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 6);

    auto ta1_x_zzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 9);

    auto ta1_y_xxx_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 10);

    auto ta1_y_yyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 16);

    auto ta1_y_zzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 19);

    auto ta1_z_xxx_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 20);

    auto ta1_z_yyy_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 26);

    auto ta1_z_zzz_0_1 = pbuffer.data(idx_npot_geom_010_1_fs + 29);

    // Set up components of auxiliary buffer : FP

    auto ta_xxx_x_1 = pbuffer.data(idx_npot_1_fp);

    auto ta_xxx_y_1 = pbuffer.data(idx_npot_1_fp + 1);

    auto ta_xxx_z_1 = pbuffer.data(idx_npot_1_fp + 2);

    auto ta_xxy_x_1 = pbuffer.data(idx_npot_1_fp + 3);

    auto ta_xxy_y_1 = pbuffer.data(idx_npot_1_fp + 4);

    auto ta_xxz_x_1 = pbuffer.data(idx_npot_1_fp + 6);

    auto ta_xxz_z_1 = pbuffer.data(idx_npot_1_fp + 8);

    auto ta_xyy_x_1 = pbuffer.data(idx_npot_1_fp + 9);

    auto ta_xyy_y_1 = pbuffer.data(idx_npot_1_fp + 10);

    auto ta_xzz_x_1 = pbuffer.data(idx_npot_1_fp + 15);

    auto ta_xzz_z_1 = pbuffer.data(idx_npot_1_fp + 17);

    auto ta_yyy_x_1 = pbuffer.data(idx_npot_1_fp + 18);

    auto ta_yyy_y_1 = pbuffer.data(idx_npot_1_fp + 19);

    auto ta_yyy_z_1 = pbuffer.data(idx_npot_1_fp + 20);

    auto ta_yyz_y_1 = pbuffer.data(idx_npot_1_fp + 22);

    auto ta_yyz_z_1 = pbuffer.data(idx_npot_1_fp + 23);

    auto ta_yzz_y_1 = pbuffer.data(idx_npot_1_fp + 25);

    auto ta_yzz_z_1 = pbuffer.data(idx_npot_1_fp + 26);

    auto ta_zzz_x_1 = pbuffer.data(idx_npot_1_fp + 27);

    auto ta_zzz_y_1 = pbuffer.data(idx_npot_1_fp + 28);

    auto ta_zzz_z_1 = pbuffer.data(idx_npot_1_fp + 29);

    // Set up components of auxiliary buffer : FP

    auto ta1_x_xxx_x_0 = pbuffer.data(idx_npot_geom_010_0_fp);

    auto ta1_x_xxx_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 1);

    auto ta1_x_xxx_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 2);

    auto ta1_x_xxy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 3);

    auto ta1_x_xxy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 4);

    auto ta1_x_xxy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 5);

    auto ta1_x_xxz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 6);

    auto ta1_x_xxz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 7);

    auto ta1_x_xxz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 8);

    auto ta1_x_xyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 9);

    auto ta1_x_xyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 10);

    auto ta1_x_xzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 15);

    auto ta1_x_xzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 17);

    auto ta1_x_yyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 18);

    auto ta1_x_yyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 19);

    auto ta1_x_yyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 20);

    auto ta1_x_yyz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 22);

    auto ta1_x_yyz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 23);

    auto ta1_x_yzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 24);

    auto ta1_x_yzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 25);

    auto ta1_x_yzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 26);

    auto ta1_x_zzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 27);

    auto ta1_x_zzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 28);

    auto ta1_x_zzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 29);

    auto ta1_y_xxx_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 30);

    auto ta1_y_xxx_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 31);

    auto ta1_y_xxx_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 32);

    auto ta1_y_xxy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 33);

    auto ta1_y_xxy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 34);

    auto ta1_y_xxz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 36);

    auto ta1_y_xxz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 38);

    auto ta1_y_xyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 39);

    auto ta1_y_xyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 40);

    auto ta1_y_xyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 41);

    auto ta1_y_xzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 45);

    auto ta1_y_xzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 46);

    auto ta1_y_xzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 47);

    auto ta1_y_yyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 48);

    auto ta1_y_yyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 49);

    auto ta1_y_yyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 50);

    auto ta1_y_yyz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 51);

    auto ta1_y_yyz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 52);

    auto ta1_y_yyz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 53);

    auto ta1_y_yzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 55);

    auto ta1_y_yzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 56);

    auto ta1_y_zzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 57);

    auto ta1_y_zzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 58);

    auto ta1_y_zzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 59);

    auto ta1_z_xxx_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 60);

    auto ta1_z_xxx_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 61);

    auto ta1_z_xxx_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 62);

    auto ta1_z_xxy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 63);

    auto ta1_z_xxy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 64);

    auto ta1_z_xxz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 66);

    auto ta1_z_xxz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 68);

    auto ta1_z_xyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 69);

    auto ta1_z_xyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 70);

    auto ta1_z_xyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 71);

    auto ta1_z_xzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 75);

    auto ta1_z_xzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 76);

    auto ta1_z_xzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 77);

    auto ta1_z_yyy_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 78);

    auto ta1_z_yyy_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 79);

    auto ta1_z_yyy_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 80);

    auto ta1_z_yyz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 82);

    auto ta1_z_yyz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 83);

    auto ta1_z_yzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 84);

    auto ta1_z_yzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 85);

    auto ta1_z_yzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 86);

    auto ta1_z_zzz_x_0 = pbuffer.data(idx_npot_geom_010_0_fp + 87);

    auto ta1_z_zzz_y_0 = pbuffer.data(idx_npot_geom_010_0_fp + 88);

    auto ta1_z_zzz_z_0 = pbuffer.data(idx_npot_geom_010_0_fp + 89);

    // Set up components of auxiliary buffer : FP

    auto ta1_x_xxx_x_1 = pbuffer.data(idx_npot_geom_010_1_fp);

    auto ta1_x_xxx_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 1);

    auto ta1_x_xxx_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 2);

    auto ta1_x_xxy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 3);

    auto ta1_x_xxy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 4);

    auto ta1_x_xxy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 5);

    auto ta1_x_xxz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 6);

    auto ta1_x_xxz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 7);

    auto ta1_x_xxz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 8);

    auto ta1_x_xyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 9);

    auto ta1_x_xyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 10);

    auto ta1_x_xzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 15);

    auto ta1_x_xzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 17);

    auto ta1_x_yyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 18);

    auto ta1_x_yyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 19);

    auto ta1_x_yyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 20);

    auto ta1_x_yyz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 22);

    auto ta1_x_yyz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 23);

    auto ta1_x_yzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 24);

    auto ta1_x_yzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 25);

    auto ta1_x_yzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 26);

    auto ta1_x_zzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 27);

    auto ta1_x_zzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 28);

    auto ta1_x_zzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 29);

    auto ta1_y_xxx_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 30);

    auto ta1_y_xxx_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 31);

    auto ta1_y_xxx_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 32);

    auto ta1_y_xxy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 33);

    auto ta1_y_xxy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 34);

    auto ta1_y_xxz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 36);

    auto ta1_y_xxz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 38);

    auto ta1_y_xyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 39);

    auto ta1_y_xyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 40);

    auto ta1_y_xyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 41);

    auto ta1_y_xzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 45);

    auto ta1_y_xzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 46);

    auto ta1_y_xzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 47);

    auto ta1_y_yyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 48);

    auto ta1_y_yyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 49);

    auto ta1_y_yyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 50);

    auto ta1_y_yyz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 51);

    auto ta1_y_yyz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 52);

    auto ta1_y_yyz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 53);

    auto ta1_y_yzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 55);

    auto ta1_y_yzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 56);

    auto ta1_y_zzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 57);

    auto ta1_y_zzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 58);

    auto ta1_y_zzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 59);

    auto ta1_z_xxx_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 60);

    auto ta1_z_xxx_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 61);

    auto ta1_z_xxx_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 62);

    auto ta1_z_xxy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 63);

    auto ta1_z_xxy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 64);

    auto ta1_z_xxz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 66);

    auto ta1_z_xxz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 68);

    auto ta1_z_xyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 69);

    auto ta1_z_xyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 70);

    auto ta1_z_xyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 71);

    auto ta1_z_xzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 75);

    auto ta1_z_xzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 76);

    auto ta1_z_xzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 77);

    auto ta1_z_yyy_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 78);

    auto ta1_z_yyy_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 79);

    auto ta1_z_yyy_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 80);

    auto ta1_z_yyz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 82);

    auto ta1_z_yyz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 83);

    auto ta1_z_yzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 84);

    auto ta1_z_yzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 85);

    auto ta1_z_yzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 86);

    auto ta1_z_zzz_x_1 = pbuffer.data(idx_npot_geom_010_1_fp + 87);

    auto ta1_z_zzz_y_1 = pbuffer.data(idx_npot_geom_010_1_fp + 88);

    auto ta1_z_zzz_z_1 = pbuffer.data(idx_npot_geom_010_1_fp + 89);

    // Set up 0-3 components of targeted buffer : GP

    auto ta1_x_xxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_gp);

    auto ta1_x_xxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 1);

    auto ta1_x_xxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 2);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_x_xx_x_0,   \
                             ta1_x_xx_x_1,   \
                             ta1_x_xx_y_0,   \
                             ta1_x_xx_y_1,   \
                             ta1_x_xx_z_0,   \
                             ta1_x_xx_z_1,   \
                             ta1_x_xxx_0_0,  \
                             ta1_x_xxx_0_1,  \
                             ta1_x_xxx_x_0,  \
                             ta1_x_xxx_x_1,  \
                             ta1_x_xxx_y_0,  \
                             ta1_x_xxx_y_1,  \
                             ta1_x_xxx_z_0,  \
                             ta1_x_xxx_z_1,  \
                             ta1_x_xxxx_x_0, \
                             ta1_x_xxxx_y_0, \
                             ta1_x_xxxx_z_0, \
                             ta_xxx_x_1,     \
                             ta_xxx_y_1,     \
                             ta_xxx_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxx_x_0[i] = 3.0 * ta1_x_xx_x_0[i] * fe_0 - 3.0 * ta1_x_xx_x_1[i] * fe_0 + ta1_x_xxx_0_0[i] * fe_0 - ta1_x_xxx_0_1[i] * fe_0 +
                            ta_xxx_x_1[i] + ta1_x_xxx_x_0[i] * pa_x[i] - ta1_x_xxx_x_1[i] * pc_x[i];

        ta1_x_xxxx_y_0[i] =
            3.0 * ta1_x_xx_y_0[i] * fe_0 - 3.0 * ta1_x_xx_y_1[i] * fe_0 + ta_xxx_y_1[i] + ta1_x_xxx_y_0[i] * pa_x[i] - ta1_x_xxx_y_1[i] * pc_x[i];

        ta1_x_xxxx_z_0[i] =
            3.0 * ta1_x_xx_z_0[i] * fe_0 - 3.0 * ta1_x_xx_z_1[i] * fe_0 + ta_xxx_z_1[i] + ta1_x_xxx_z_0[i] * pa_x[i] - ta1_x_xxx_z_1[i] * pc_x[i];
    }

    // Set up 3-6 components of targeted buffer : GP

    auto ta1_x_xxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 3);

    auto ta1_x_xxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 4);

    auto ta1_x_xxxy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 5);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_x_xxx_0_0,  \
                             ta1_x_xxx_0_1,  \
                             ta1_x_xxx_x_0,  \
                             ta1_x_xxx_x_1,  \
                             ta1_x_xxx_y_0,  \
                             ta1_x_xxx_y_1,  \
                             ta1_x_xxx_z_0,  \
                             ta1_x_xxx_z_1,  \
                             ta1_x_xxxy_x_0, \
                             ta1_x_xxxy_y_0, \
                             ta1_x_xxxy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxy_x_0[i] = ta1_x_xxx_x_0[i] * pa_y[i] - ta1_x_xxx_x_1[i] * pc_y[i];

        ta1_x_xxxy_y_0[i] = ta1_x_xxx_0_0[i] * fe_0 - ta1_x_xxx_0_1[i] * fe_0 + ta1_x_xxx_y_0[i] * pa_y[i] - ta1_x_xxx_y_1[i] * pc_y[i];

        ta1_x_xxxy_z_0[i] = ta1_x_xxx_z_0[i] * pa_y[i] - ta1_x_xxx_z_1[i] * pc_y[i];
    }

    // Set up 6-9 components of targeted buffer : GP

    auto ta1_x_xxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 6);

    auto ta1_x_xxxz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 7);

    auto ta1_x_xxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 8);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_x_xxx_0_0,  \
                             ta1_x_xxx_0_1,  \
                             ta1_x_xxx_x_0,  \
                             ta1_x_xxx_x_1,  \
                             ta1_x_xxx_y_0,  \
                             ta1_x_xxx_y_1,  \
                             ta1_x_xxx_z_0,  \
                             ta1_x_xxx_z_1,  \
                             ta1_x_xxxz_x_0, \
                             ta1_x_xxxz_y_0, \
                             ta1_x_xxxz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxz_x_0[i] = ta1_x_xxx_x_0[i] * pa_z[i] - ta1_x_xxx_x_1[i] * pc_z[i];

        ta1_x_xxxz_y_0[i] = ta1_x_xxx_y_0[i] * pa_z[i] - ta1_x_xxx_y_1[i] * pc_z[i];

        ta1_x_xxxz_z_0[i] = ta1_x_xxx_0_0[i] * fe_0 - ta1_x_xxx_0_1[i] * fe_0 + ta1_x_xxx_z_0[i] * pa_z[i] - ta1_x_xxx_z_1[i] * pc_z[i];
    }

    // Set up 9-12 components of targeted buffer : GP

    auto ta1_x_xxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 9);

    auto ta1_x_xxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 10);

    auto ta1_x_xxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 11);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_x_xx_x_0,   \
                             ta1_x_xx_x_1,   \
                             ta1_x_xx_z_0,   \
                             ta1_x_xx_z_1,   \
                             ta1_x_xxy_x_0,  \
                             ta1_x_xxy_x_1,  \
                             ta1_x_xxy_z_0,  \
                             ta1_x_xxy_z_1,  \
                             ta1_x_xxyy_x_0, \
                             ta1_x_xxyy_y_0, \
                             ta1_x_xxyy_z_0, \
                             ta1_x_xyy_y_0,  \
                             ta1_x_xyy_y_1,  \
                             ta1_x_yy_y_0,   \
                             ta1_x_yy_y_1,   \
                             ta_xyy_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyy_x_0[i] = ta1_x_xx_x_0[i] * fe_0 - ta1_x_xx_x_1[i] * fe_0 + ta1_x_xxy_x_0[i] * pa_y[i] - ta1_x_xxy_x_1[i] * pc_y[i];

        ta1_x_xxyy_y_0[i] = ta1_x_yy_y_0[i] * fe_0 - ta1_x_yy_y_1[i] * fe_0 + ta_xyy_y_1[i] + ta1_x_xyy_y_0[i] * pa_x[i] - ta1_x_xyy_y_1[i] * pc_x[i];

        ta1_x_xxyy_z_0[i] = ta1_x_xx_z_0[i] * fe_0 - ta1_x_xx_z_1[i] * fe_0 + ta1_x_xxy_z_0[i] * pa_y[i] - ta1_x_xxy_z_1[i] * pc_y[i];
    }

    // Set up 12-15 components of targeted buffer : GP

    auto ta1_x_xxyz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 12);

    auto ta1_x_xxyz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 13);

    auto ta1_x_xxyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 14);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_x_xxy_y_0,  \
                             ta1_x_xxy_y_1,  \
                             ta1_x_xxyz_x_0, \
                             ta1_x_xxyz_y_0, \
                             ta1_x_xxyz_z_0, \
                             ta1_x_xxz_x_0,  \
                             ta1_x_xxz_x_1,  \
                             ta1_x_xxz_z_0,  \
                             ta1_x_xxz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_x_xxyz_x_0[i] = ta1_x_xxz_x_0[i] * pa_y[i] - ta1_x_xxz_x_1[i] * pc_y[i];

        ta1_x_xxyz_y_0[i] = ta1_x_xxy_y_0[i] * pa_z[i] - ta1_x_xxy_y_1[i] * pc_z[i];

        ta1_x_xxyz_z_0[i] = ta1_x_xxz_z_0[i] * pa_y[i] - ta1_x_xxz_z_1[i] * pc_y[i];
    }

    // Set up 15-18 components of targeted buffer : GP

    auto ta1_x_xxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 15);

    auto ta1_x_xxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 16);

    auto ta1_x_xxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 17);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_x_xx_x_0,   \
                             ta1_x_xx_x_1,   \
                             ta1_x_xx_y_0,   \
                             ta1_x_xx_y_1,   \
                             ta1_x_xxz_x_0,  \
                             ta1_x_xxz_x_1,  \
                             ta1_x_xxz_y_0,  \
                             ta1_x_xxz_y_1,  \
                             ta1_x_xxzz_x_0, \
                             ta1_x_xxzz_y_0, \
                             ta1_x_xxzz_z_0, \
                             ta1_x_xzz_z_0,  \
                             ta1_x_xzz_z_1,  \
                             ta1_x_zz_z_0,   \
                             ta1_x_zz_z_1,   \
                             ta_xzz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxzz_x_0[i] = ta1_x_xx_x_0[i] * fe_0 - ta1_x_xx_x_1[i] * fe_0 + ta1_x_xxz_x_0[i] * pa_z[i] - ta1_x_xxz_x_1[i] * pc_z[i];

        ta1_x_xxzz_y_0[i] = ta1_x_xx_y_0[i] * fe_0 - ta1_x_xx_y_1[i] * fe_0 + ta1_x_xxz_y_0[i] * pa_z[i] - ta1_x_xxz_y_1[i] * pc_z[i];

        ta1_x_xxzz_z_0[i] = ta1_x_zz_z_0[i] * fe_0 - ta1_x_zz_z_1[i] * fe_0 + ta_xzz_z_1[i] + ta1_x_xzz_z_0[i] * pa_x[i] - ta1_x_xzz_z_1[i] * pc_x[i];
    }

    // Set up 18-21 components of targeted buffer : GP

    auto ta1_x_xyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 18);

    auto ta1_x_xyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 19);

    auto ta1_x_xyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 20);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_x_xy_x_0,   \
                             ta1_x_xy_x_1,   \
                             ta1_x_xyy_x_0,  \
                             ta1_x_xyy_x_1,  \
                             ta1_x_xyyy_x_0, \
                             ta1_x_xyyy_y_0, \
                             ta1_x_xyyy_z_0, \
                             ta1_x_yyy_y_0,  \
                             ta1_x_yyy_y_1,  \
                             ta1_x_yyy_z_0,  \
                             ta1_x_yyy_z_1,  \
                             ta_yyy_y_1,     \
                             ta_yyy_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyy_x_0[i] = 2.0 * ta1_x_xy_x_0[i] * fe_0 - 2.0 * ta1_x_xy_x_1[i] * fe_0 + ta1_x_xyy_x_0[i] * pa_y[i] - ta1_x_xyy_x_1[i] * pc_y[i];

        ta1_x_xyyy_y_0[i] = ta_yyy_y_1[i] + ta1_x_yyy_y_0[i] * pa_x[i] - ta1_x_yyy_y_1[i] * pc_x[i];

        ta1_x_xyyy_z_0[i] = ta_yyy_z_1[i] + ta1_x_yyy_z_0[i] * pa_x[i] - ta1_x_yyy_z_1[i] * pc_x[i];
    }

    // Set up 21-24 components of targeted buffer : GP

    auto ta1_x_xyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 21);

    auto ta1_x_xyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 22);

    auto ta1_x_xyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 23);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_x_xyy_x_0,  \
                             ta1_x_xyy_x_1,  \
                             ta1_x_xyy_y_0,  \
                             ta1_x_xyy_y_1,  \
                             ta1_x_xyyz_x_0, \
                             ta1_x_xyyz_y_0, \
                             ta1_x_xyyz_z_0, \
                             ta1_x_yyz_z_0,  \
                             ta1_x_yyz_z_1,  \
                             ta_yyz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_x_xyyz_x_0[i] = ta1_x_xyy_x_0[i] * pa_z[i] - ta1_x_xyy_x_1[i] * pc_z[i];

        ta1_x_xyyz_y_0[i] = ta1_x_xyy_y_0[i] * pa_z[i] - ta1_x_xyy_y_1[i] * pc_z[i];

        ta1_x_xyyz_z_0[i] = ta_yyz_z_1[i] + ta1_x_yyz_z_0[i] * pa_x[i] - ta1_x_yyz_z_1[i] * pc_x[i];
    }

    // Set up 24-27 components of targeted buffer : GP

    auto ta1_x_xyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 24);

    auto ta1_x_xyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 25);

    auto ta1_x_xyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 26);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_x_xyzz_x_0, \
                             ta1_x_xyzz_y_0, \
                             ta1_x_xyzz_z_0, \
                             ta1_x_xzz_x_0,  \
                             ta1_x_xzz_x_1,  \
                             ta1_x_xzz_z_0,  \
                             ta1_x_xzz_z_1,  \
                             ta1_x_yzz_y_0,  \
                             ta1_x_yzz_y_1,  \
                             ta_yzz_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_x_xyzz_x_0[i] = ta1_x_xzz_x_0[i] * pa_y[i] - ta1_x_xzz_x_1[i] * pc_y[i];

        ta1_x_xyzz_y_0[i] = ta_yzz_y_1[i] + ta1_x_yzz_y_0[i] * pa_x[i] - ta1_x_yzz_y_1[i] * pc_x[i];

        ta1_x_xyzz_z_0[i] = ta1_x_xzz_z_0[i] * pa_y[i] - ta1_x_xzz_z_1[i] * pc_y[i];
    }

    // Set up 27-30 components of targeted buffer : GP

    auto ta1_x_xzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 27);

    auto ta1_x_xzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 28);

    auto ta1_x_xzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 29);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_x_xz_x_0,   \
                             ta1_x_xz_x_1,   \
                             ta1_x_xzz_x_0,  \
                             ta1_x_xzz_x_1,  \
                             ta1_x_xzzz_x_0, \
                             ta1_x_xzzz_y_0, \
                             ta1_x_xzzz_z_0, \
                             ta1_x_zzz_y_0,  \
                             ta1_x_zzz_y_1,  \
                             ta1_x_zzz_z_0,  \
                             ta1_x_zzz_z_1,  \
                             ta_zzz_y_1,     \
                             ta_zzz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzzz_x_0[i] = 2.0 * ta1_x_xz_x_0[i] * fe_0 - 2.0 * ta1_x_xz_x_1[i] * fe_0 + ta1_x_xzz_x_0[i] * pa_z[i] - ta1_x_xzz_x_1[i] * pc_z[i];

        ta1_x_xzzz_y_0[i] = ta_zzz_y_1[i] + ta1_x_zzz_y_0[i] * pa_x[i] - ta1_x_zzz_y_1[i] * pc_x[i];

        ta1_x_xzzz_z_0[i] = ta_zzz_z_1[i] + ta1_x_zzz_z_0[i] * pa_x[i] - ta1_x_zzz_z_1[i] * pc_x[i];
    }

    // Set up 30-33 components of targeted buffer : GP

    auto ta1_x_yyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 30);

    auto ta1_x_yyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 31);

    auto ta1_x_yyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 32);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_x_yy_x_0,   \
                             ta1_x_yy_x_1,   \
                             ta1_x_yy_y_0,   \
                             ta1_x_yy_y_1,   \
                             ta1_x_yy_z_0,   \
                             ta1_x_yy_z_1,   \
                             ta1_x_yyy_0_0,  \
                             ta1_x_yyy_0_1,  \
                             ta1_x_yyy_x_0,  \
                             ta1_x_yyy_x_1,  \
                             ta1_x_yyy_y_0,  \
                             ta1_x_yyy_y_1,  \
                             ta1_x_yyy_z_0,  \
                             ta1_x_yyy_z_1,  \
                             ta1_x_yyyy_x_0, \
                             ta1_x_yyyy_y_0, \
                             ta1_x_yyyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyy_x_0[i] = 3.0 * ta1_x_yy_x_0[i] * fe_0 - 3.0 * ta1_x_yy_x_1[i] * fe_0 + ta1_x_yyy_x_0[i] * pa_y[i] - ta1_x_yyy_x_1[i] * pc_y[i];

        ta1_x_yyyy_y_0[i] = 3.0 * ta1_x_yy_y_0[i] * fe_0 - 3.0 * ta1_x_yy_y_1[i] * fe_0 + ta1_x_yyy_0_0[i] * fe_0 - ta1_x_yyy_0_1[i] * fe_0 +
                            ta1_x_yyy_y_0[i] * pa_y[i] - ta1_x_yyy_y_1[i] * pc_y[i];

        ta1_x_yyyy_z_0[i] = 3.0 * ta1_x_yy_z_0[i] * fe_0 - 3.0 * ta1_x_yy_z_1[i] * fe_0 + ta1_x_yyy_z_0[i] * pa_y[i] - ta1_x_yyy_z_1[i] * pc_y[i];
    }

    // Set up 33-36 components of targeted buffer : GP

    auto ta1_x_yyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 33);

    auto ta1_x_yyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 34);

    auto ta1_x_yyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 35);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_x_yyy_x_0,  \
                             ta1_x_yyy_x_1,  \
                             ta1_x_yyy_y_0,  \
                             ta1_x_yyy_y_1,  \
                             ta1_x_yyyz_x_0, \
                             ta1_x_yyyz_y_0, \
                             ta1_x_yyyz_z_0, \
                             ta1_x_yyz_z_0,  \
                             ta1_x_yyz_z_1,  \
                             ta1_x_yz_z_0,   \
                             ta1_x_yz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyz_x_0[i] = ta1_x_yyy_x_0[i] * pa_z[i] - ta1_x_yyy_x_1[i] * pc_z[i];

        ta1_x_yyyz_y_0[i] = ta1_x_yyy_y_0[i] * pa_z[i] - ta1_x_yyy_y_1[i] * pc_z[i];

        ta1_x_yyyz_z_0[i] = 2.0 * ta1_x_yz_z_0[i] * fe_0 - 2.0 * ta1_x_yz_z_1[i] * fe_0 + ta1_x_yyz_z_0[i] * pa_y[i] - ta1_x_yyz_z_1[i] * pc_y[i];
    }

    // Set up 36-39 components of targeted buffer : GP

    auto ta1_x_yyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 36);

    auto ta1_x_yyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 37);

    auto ta1_x_yyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 38);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_x_yy_y_0,   \
                             ta1_x_yy_y_1,   \
                             ta1_x_yyz_y_0,  \
                             ta1_x_yyz_y_1,  \
                             ta1_x_yyzz_x_0, \
                             ta1_x_yyzz_y_0, \
                             ta1_x_yyzz_z_0, \
                             ta1_x_yzz_x_0,  \
                             ta1_x_yzz_x_1,  \
                             ta1_x_yzz_z_0,  \
                             ta1_x_yzz_z_1,  \
                             ta1_x_zz_x_0,   \
                             ta1_x_zz_x_1,   \
                             ta1_x_zz_z_0,   \
                             ta1_x_zz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyzz_x_0[i] = ta1_x_zz_x_0[i] * fe_0 - ta1_x_zz_x_1[i] * fe_0 + ta1_x_yzz_x_0[i] * pa_y[i] - ta1_x_yzz_x_1[i] * pc_y[i];

        ta1_x_yyzz_y_0[i] = ta1_x_yy_y_0[i] * fe_0 - ta1_x_yy_y_1[i] * fe_0 + ta1_x_yyz_y_0[i] * pa_z[i] - ta1_x_yyz_y_1[i] * pc_z[i];

        ta1_x_yyzz_z_0[i] = ta1_x_zz_z_0[i] * fe_0 - ta1_x_zz_z_1[i] * fe_0 + ta1_x_yzz_z_0[i] * pa_y[i] - ta1_x_yzz_z_1[i] * pc_y[i];
    }

    // Set up 39-42 components of targeted buffer : GP

    auto ta1_x_yzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 39);

    auto ta1_x_yzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 40);

    auto ta1_x_yzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 41);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_x_yzzz_x_0, \
                             ta1_x_yzzz_y_0, \
                             ta1_x_yzzz_z_0, \
                             ta1_x_zzz_0_0,  \
                             ta1_x_zzz_0_1,  \
                             ta1_x_zzz_x_0,  \
                             ta1_x_zzz_x_1,  \
                             ta1_x_zzz_y_0,  \
                             ta1_x_zzz_y_1,  \
                             ta1_x_zzz_z_0,  \
                             ta1_x_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzzz_x_0[i] = ta1_x_zzz_x_0[i] * pa_y[i] - ta1_x_zzz_x_1[i] * pc_y[i];

        ta1_x_yzzz_y_0[i] = ta1_x_zzz_0_0[i] * fe_0 - ta1_x_zzz_0_1[i] * fe_0 + ta1_x_zzz_y_0[i] * pa_y[i] - ta1_x_zzz_y_1[i] * pc_y[i];

        ta1_x_yzzz_z_0[i] = ta1_x_zzz_z_0[i] * pa_y[i] - ta1_x_zzz_z_1[i] * pc_y[i];
    }

    // Set up 42-45 components of targeted buffer : GP

    auto ta1_x_zzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 42);

    auto ta1_x_zzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 43);

    auto ta1_x_zzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 44);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_x_zz_x_0,   \
                             ta1_x_zz_x_1,   \
                             ta1_x_zz_y_0,   \
                             ta1_x_zz_y_1,   \
                             ta1_x_zz_z_0,   \
                             ta1_x_zz_z_1,   \
                             ta1_x_zzz_0_0,  \
                             ta1_x_zzz_0_1,  \
                             ta1_x_zzz_x_0,  \
                             ta1_x_zzz_x_1,  \
                             ta1_x_zzz_y_0,  \
                             ta1_x_zzz_y_1,  \
                             ta1_x_zzz_z_0,  \
                             ta1_x_zzz_z_1,  \
                             ta1_x_zzzz_x_0, \
                             ta1_x_zzzz_y_0, \
                             ta1_x_zzzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzzz_x_0[i] = 3.0 * ta1_x_zz_x_0[i] * fe_0 - 3.0 * ta1_x_zz_x_1[i] * fe_0 + ta1_x_zzz_x_0[i] * pa_z[i] - ta1_x_zzz_x_1[i] * pc_z[i];

        ta1_x_zzzz_y_0[i] = 3.0 * ta1_x_zz_y_0[i] * fe_0 - 3.0 * ta1_x_zz_y_1[i] * fe_0 + ta1_x_zzz_y_0[i] * pa_z[i] - ta1_x_zzz_y_1[i] * pc_z[i];

        ta1_x_zzzz_z_0[i] = 3.0 * ta1_x_zz_z_0[i] * fe_0 - 3.0 * ta1_x_zz_z_1[i] * fe_0 + ta1_x_zzz_0_0[i] * fe_0 - ta1_x_zzz_0_1[i] * fe_0 +
                            ta1_x_zzz_z_0[i] * pa_z[i] - ta1_x_zzz_z_1[i] * pc_z[i];
    }

    // Set up 45-48 components of targeted buffer : GP

    auto ta1_y_xxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 45);

    auto ta1_y_xxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 46);

    auto ta1_y_xxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 47);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_y_xx_x_0,   \
                             ta1_y_xx_x_1,   \
                             ta1_y_xx_y_0,   \
                             ta1_y_xx_y_1,   \
                             ta1_y_xx_z_0,   \
                             ta1_y_xx_z_1,   \
                             ta1_y_xxx_0_0,  \
                             ta1_y_xxx_0_1,  \
                             ta1_y_xxx_x_0,  \
                             ta1_y_xxx_x_1,  \
                             ta1_y_xxx_y_0,  \
                             ta1_y_xxx_y_1,  \
                             ta1_y_xxx_z_0,  \
                             ta1_y_xxx_z_1,  \
                             ta1_y_xxxx_x_0, \
                             ta1_y_xxxx_y_0, \
                             ta1_y_xxxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxx_x_0[i] = 3.0 * ta1_y_xx_x_0[i] * fe_0 - 3.0 * ta1_y_xx_x_1[i] * fe_0 + ta1_y_xxx_0_0[i] * fe_0 - ta1_y_xxx_0_1[i] * fe_0 +
                            ta1_y_xxx_x_0[i] * pa_x[i] - ta1_y_xxx_x_1[i] * pc_x[i];

        ta1_y_xxxx_y_0[i] = 3.0 * ta1_y_xx_y_0[i] * fe_0 - 3.0 * ta1_y_xx_y_1[i] * fe_0 + ta1_y_xxx_y_0[i] * pa_x[i] - ta1_y_xxx_y_1[i] * pc_x[i];

        ta1_y_xxxx_z_0[i] = 3.0 * ta1_y_xx_z_0[i] * fe_0 - 3.0 * ta1_y_xx_z_1[i] * fe_0 + ta1_y_xxx_z_0[i] * pa_x[i] - ta1_y_xxx_z_1[i] * pc_x[i];
    }

    // Set up 48-51 components of targeted buffer : GP

    auto ta1_y_xxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 48);

    auto ta1_y_xxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 49);

    auto ta1_y_xxxy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 50);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_y_xxx_x_0,  \
                             ta1_y_xxx_x_1,  \
                             ta1_y_xxx_z_0,  \
                             ta1_y_xxx_z_1,  \
                             ta1_y_xxxy_x_0, \
                             ta1_y_xxxy_y_0, \
                             ta1_y_xxxy_z_0, \
                             ta1_y_xxy_y_0,  \
                             ta1_y_xxy_y_1,  \
                             ta1_y_xy_y_0,   \
                             ta1_y_xy_y_1,   \
                             ta_xxx_x_1,     \
                             ta_xxx_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxy_x_0[i] = ta_xxx_x_1[i] + ta1_y_xxx_x_0[i] * pa_y[i] - ta1_y_xxx_x_1[i] * pc_y[i];

        ta1_y_xxxy_y_0[i] = 2.0 * ta1_y_xy_y_0[i] * fe_0 - 2.0 * ta1_y_xy_y_1[i] * fe_0 + ta1_y_xxy_y_0[i] * pa_x[i] - ta1_y_xxy_y_1[i] * pc_x[i];

        ta1_y_xxxy_z_0[i] = ta_xxx_z_1[i] + ta1_y_xxx_z_0[i] * pa_y[i] - ta1_y_xxx_z_1[i] * pc_y[i];
    }

    // Set up 51-54 components of targeted buffer : GP

    auto ta1_y_xxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 51);

    auto ta1_y_xxxz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 52);

    auto ta1_y_xxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 53);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_y_xxx_x_0,  \
                             ta1_y_xxx_x_1,  \
                             ta1_y_xxx_y_0,  \
                             ta1_y_xxx_y_1,  \
                             ta1_y_xxxz_x_0, \
                             ta1_y_xxxz_y_0, \
                             ta1_y_xxxz_z_0, \
                             ta1_y_xxz_z_0,  \
                             ta1_y_xxz_z_1,  \
                             ta1_y_xz_z_0,   \
                             ta1_y_xz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxz_x_0[i] = ta1_y_xxx_x_0[i] * pa_z[i] - ta1_y_xxx_x_1[i] * pc_z[i];

        ta1_y_xxxz_y_0[i] = ta1_y_xxx_y_0[i] * pa_z[i] - ta1_y_xxx_y_1[i] * pc_z[i];

        ta1_y_xxxz_z_0[i] = 2.0 * ta1_y_xz_z_0[i] * fe_0 - 2.0 * ta1_y_xz_z_1[i] * fe_0 + ta1_y_xxz_z_0[i] * pa_x[i] - ta1_y_xxz_z_1[i] * pc_x[i];
    }

    // Set up 54-57 components of targeted buffer : GP

    auto ta1_y_xxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 54);

    auto ta1_y_xxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 55);

    auto ta1_y_xxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 56);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_y_xx_x_0,   \
                             ta1_y_xx_x_1,   \
                             ta1_y_xxy_x_0,  \
                             ta1_y_xxy_x_1,  \
                             ta1_y_xxyy_x_0, \
                             ta1_y_xxyy_y_0, \
                             ta1_y_xxyy_z_0, \
                             ta1_y_xyy_y_0,  \
                             ta1_y_xyy_y_1,  \
                             ta1_y_xyy_z_0,  \
                             ta1_y_xyy_z_1,  \
                             ta1_y_yy_y_0,   \
                             ta1_y_yy_y_1,   \
                             ta1_y_yy_z_0,   \
                             ta1_y_yy_z_1,   \
                             ta_xxy_x_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyy_x_0[i] = ta1_y_xx_x_0[i] * fe_0 - ta1_y_xx_x_1[i] * fe_0 + ta_xxy_x_1[i] + ta1_y_xxy_x_0[i] * pa_y[i] - ta1_y_xxy_x_1[i] * pc_y[i];

        ta1_y_xxyy_y_0[i] = ta1_y_yy_y_0[i] * fe_0 - ta1_y_yy_y_1[i] * fe_0 + ta1_y_xyy_y_0[i] * pa_x[i] - ta1_y_xyy_y_1[i] * pc_x[i];

        ta1_y_xxyy_z_0[i] = ta1_y_yy_z_0[i] * fe_0 - ta1_y_yy_z_1[i] * fe_0 + ta1_y_xyy_z_0[i] * pa_x[i] - ta1_y_xyy_z_1[i] * pc_x[i];
    }

    // Set up 57-60 components of targeted buffer : GP

    auto ta1_y_xxyz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 57);

    auto ta1_y_xxyz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 58);

    auto ta1_y_xxyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 59);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_y_xxy_x_0,  \
                             ta1_y_xxy_x_1,  \
                             ta1_y_xxy_y_0,  \
                             ta1_y_xxy_y_1,  \
                             ta1_y_xxyz_x_0, \
                             ta1_y_xxyz_y_0, \
                             ta1_y_xxyz_z_0, \
                             ta1_y_xxz_z_0,  \
                             ta1_y_xxz_z_1,  \
                             ta_xxz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_y_xxyz_x_0[i] = ta1_y_xxy_x_0[i] * pa_z[i] - ta1_y_xxy_x_1[i] * pc_z[i];

        ta1_y_xxyz_y_0[i] = ta1_y_xxy_y_0[i] * pa_z[i] - ta1_y_xxy_y_1[i] * pc_z[i];

        ta1_y_xxyz_z_0[i] = ta_xxz_z_1[i] + ta1_y_xxz_z_0[i] * pa_y[i] - ta1_y_xxz_z_1[i] * pc_y[i];
    }

    // Set up 60-63 components of targeted buffer : GP

    auto ta1_y_xxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 60);

    auto ta1_y_xxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 61);

    auto ta1_y_xxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 62);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_y_xx_x_0,   \
                             ta1_y_xx_x_1,   \
                             ta1_y_xxz_x_0,  \
                             ta1_y_xxz_x_1,  \
                             ta1_y_xxzz_x_0, \
                             ta1_y_xxzz_y_0, \
                             ta1_y_xxzz_z_0, \
                             ta1_y_xzz_y_0,  \
                             ta1_y_xzz_y_1,  \
                             ta1_y_xzz_z_0,  \
                             ta1_y_xzz_z_1,  \
                             ta1_y_zz_y_0,   \
                             ta1_y_zz_y_1,   \
                             ta1_y_zz_z_0,   \
                             ta1_y_zz_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxzz_x_0[i] = ta1_y_xx_x_0[i] * fe_0 - ta1_y_xx_x_1[i] * fe_0 + ta1_y_xxz_x_0[i] * pa_z[i] - ta1_y_xxz_x_1[i] * pc_z[i];

        ta1_y_xxzz_y_0[i] = ta1_y_zz_y_0[i] * fe_0 - ta1_y_zz_y_1[i] * fe_0 + ta1_y_xzz_y_0[i] * pa_x[i] - ta1_y_xzz_y_1[i] * pc_x[i];

        ta1_y_xxzz_z_0[i] = ta1_y_zz_z_0[i] * fe_0 - ta1_y_zz_z_1[i] * fe_0 + ta1_y_xzz_z_0[i] * pa_x[i] - ta1_y_xzz_z_1[i] * pc_x[i];
    }

    // Set up 63-66 components of targeted buffer : GP

    auto ta1_y_xyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 63);

    auto ta1_y_xyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 64);

    auto ta1_y_xyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 65);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_y_xyyy_x_0, \
                             ta1_y_xyyy_y_0, \
                             ta1_y_xyyy_z_0, \
                             ta1_y_yyy_0_0,  \
                             ta1_y_yyy_0_1,  \
                             ta1_y_yyy_x_0,  \
                             ta1_y_yyy_x_1,  \
                             ta1_y_yyy_y_0,  \
                             ta1_y_yyy_y_1,  \
                             ta1_y_yyy_z_0,  \
                             ta1_y_yyy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyy_x_0[i] = ta1_y_yyy_0_0[i] * fe_0 - ta1_y_yyy_0_1[i] * fe_0 + ta1_y_yyy_x_0[i] * pa_x[i] - ta1_y_yyy_x_1[i] * pc_x[i];

        ta1_y_xyyy_y_0[i] = ta1_y_yyy_y_0[i] * pa_x[i] - ta1_y_yyy_y_1[i] * pc_x[i];

        ta1_y_xyyy_z_0[i] = ta1_y_yyy_z_0[i] * pa_x[i] - ta1_y_yyy_z_1[i] * pc_x[i];
    }

    // Set up 66-69 components of targeted buffer : GP

    auto ta1_y_xyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 66);

    auto ta1_y_xyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 67);

    auto ta1_y_xyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 68);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_y_xyy_x_0,  \
                             ta1_y_xyy_x_1,  \
                             ta1_y_xyyz_x_0, \
                             ta1_y_xyyz_y_0, \
                             ta1_y_xyyz_z_0, \
                             ta1_y_yyz_y_0,  \
                             ta1_y_yyz_y_1,  \
                             ta1_y_yyz_z_0,  \
                             ta1_y_yyz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_y_xyyz_x_0[i] = ta1_y_xyy_x_0[i] * pa_z[i] - ta1_y_xyy_x_1[i] * pc_z[i];

        ta1_y_xyyz_y_0[i] = ta1_y_yyz_y_0[i] * pa_x[i] - ta1_y_yyz_y_1[i] * pc_x[i];

        ta1_y_xyyz_z_0[i] = ta1_y_yyz_z_0[i] * pa_x[i] - ta1_y_yyz_z_1[i] * pc_x[i];
    }

    // Set up 69-72 components of targeted buffer : GP

    auto ta1_y_xyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 69);

    auto ta1_y_xyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 70);

    auto ta1_y_xyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 71);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_y_xyzz_x_0, \
                             ta1_y_xyzz_y_0, \
                             ta1_y_xyzz_z_0, \
                             ta1_y_xzz_x_0,  \
                             ta1_y_xzz_x_1,  \
                             ta1_y_yzz_y_0,  \
                             ta1_y_yzz_y_1,  \
                             ta1_y_yzz_z_0,  \
                             ta1_y_yzz_z_1,  \
                             ta_xzz_x_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_y_xyzz_x_0[i] = ta_xzz_x_1[i] + ta1_y_xzz_x_0[i] * pa_y[i] - ta1_y_xzz_x_1[i] * pc_y[i];

        ta1_y_xyzz_y_0[i] = ta1_y_yzz_y_0[i] * pa_x[i] - ta1_y_yzz_y_1[i] * pc_x[i];

        ta1_y_xyzz_z_0[i] = ta1_y_yzz_z_0[i] * pa_x[i] - ta1_y_yzz_z_1[i] * pc_x[i];
    }

    // Set up 72-75 components of targeted buffer : GP

    auto ta1_y_xzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 72);

    auto ta1_y_xzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 73);

    auto ta1_y_xzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 74);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_y_xzzz_x_0, \
                             ta1_y_xzzz_y_0, \
                             ta1_y_xzzz_z_0, \
                             ta1_y_zzz_0_0,  \
                             ta1_y_zzz_0_1,  \
                             ta1_y_zzz_x_0,  \
                             ta1_y_zzz_x_1,  \
                             ta1_y_zzz_y_0,  \
                             ta1_y_zzz_y_1,  \
                             ta1_y_zzz_z_0,  \
                             ta1_y_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzzz_x_0[i] = ta1_y_zzz_0_0[i] * fe_0 - ta1_y_zzz_0_1[i] * fe_0 + ta1_y_zzz_x_0[i] * pa_x[i] - ta1_y_zzz_x_1[i] * pc_x[i];

        ta1_y_xzzz_y_0[i] = ta1_y_zzz_y_0[i] * pa_x[i] - ta1_y_zzz_y_1[i] * pc_x[i];

        ta1_y_xzzz_z_0[i] = ta1_y_zzz_z_0[i] * pa_x[i] - ta1_y_zzz_z_1[i] * pc_x[i];
    }

    // Set up 75-78 components of targeted buffer : GP

    auto ta1_y_yyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 75);

    auto ta1_y_yyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 76);

    auto ta1_y_yyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 77);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_y_yy_x_0,   \
                             ta1_y_yy_x_1,   \
                             ta1_y_yy_y_0,   \
                             ta1_y_yy_y_1,   \
                             ta1_y_yy_z_0,   \
                             ta1_y_yy_z_1,   \
                             ta1_y_yyy_0_0,  \
                             ta1_y_yyy_0_1,  \
                             ta1_y_yyy_x_0,  \
                             ta1_y_yyy_x_1,  \
                             ta1_y_yyy_y_0,  \
                             ta1_y_yyy_y_1,  \
                             ta1_y_yyy_z_0,  \
                             ta1_y_yyy_z_1,  \
                             ta1_y_yyyy_x_0, \
                             ta1_y_yyyy_y_0, \
                             ta1_y_yyyy_z_0, \
                             ta_yyy_x_1,     \
                             ta_yyy_y_1,     \
                             ta_yyy_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyy_x_0[i] =
            3.0 * ta1_y_yy_x_0[i] * fe_0 - 3.0 * ta1_y_yy_x_1[i] * fe_0 + ta_yyy_x_1[i] + ta1_y_yyy_x_0[i] * pa_y[i] - ta1_y_yyy_x_1[i] * pc_y[i];

        ta1_y_yyyy_y_0[i] = 3.0 * ta1_y_yy_y_0[i] * fe_0 - 3.0 * ta1_y_yy_y_1[i] * fe_0 + ta1_y_yyy_0_0[i] * fe_0 - ta1_y_yyy_0_1[i] * fe_0 +
                            ta_yyy_y_1[i] + ta1_y_yyy_y_0[i] * pa_y[i] - ta1_y_yyy_y_1[i] * pc_y[i];

        ta1_y_yyyy_z_0[i] =
            3.0 * ta1_y_yy_z_0[i] * fe_0 - 3.0 * ta1_y_yy_z_1[i] * fe_0 + ta_yyy_z_1[i] + ta1_y_yyy_z_0[i] * pa_y[i] - ta1_y_yyy_z_1[i] * pc_y[i];
    }

    // Set up 78-81 components of targeted buffer : GP

    auto ta1_y_yyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 78);

    auto ta1_y_yyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 79);

    auto ta1_y_yyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 80);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_y_yyy_0_0,  \
                             ta1_y_yyy_0_1,  \
                             ta1_y_yyy_x_0,  \
                             ta1_y_yyy_x_1,  \
                             ta1_y_yyy_y_0,  \
                             ta1_y_yyy_y_1,  \
                             ta1_y_yyy_z_0,  \
                             ta1_y_yyy_z_1,  \
                             ta1_y_yyyz_x_0, \
                             ta1_y_yyyz_y_0, \
                             ta1_y_yyyz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyz_x_0[i] = ta1_y_yyy_x_0[i] * pa_z[i] - ta1_y_yyy_x_1[i] * pc_z[i];

        ta1_y_yyyz_y_0[i] = ta1_y_yyy_y_0[i] * pa_z[i] - ta1_y_yyy_y_1[i] * pc_z[i];

        ta1_y_yyyz_z_0[i] = ta1_y_yyy_0_0[i] * fe_0 - ta1_y_yyy_0_1[i] * fe_0 + ta1_y_yyy_z_0[i] * pa_z[i] - ta1_y_yyy_z_1[i] * pc_z[i];
    }

    // Set up 81-84 components of targeted buffer : GP

    auto ta1_y_yyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 81);

    auto ta1_y_yyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 82);

    auto ta1_y_yyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 83);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_y_yy_x_0,   \
                             ta1_y_yy_x_1,   \
                             ta1_y_yy_y_0,   \
                             ta1_y_yy_y_1,   \
                             ta1_y_yyz_x_0,  \
                             ta1_y_yyz_x_1,  \
                             ta1_y_yyz_y_0,  \
                             ta1_y_yyz_y_1,  \
                             ta1_y_yyzz_x_0, \
                             ta1_y_yyzz_y_0, \
                             ta1_y_yyzz_z_0, \
                             ta1_y_yzz_z_0,  \
                             ta1_y_yzz_z_1,  \
                             ta1_y_zz_z_0,   \
                             ta1_y_zz_z_1,   \
                             ta_yzz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyzz_x_0[i] = ta1_y_yy_x_0[i] * fe_0 - ta1_y_yy_x_1[i] * fe_0 + ta1_y_yyz_x_0[i] * pa_z[i] - ta1_y_yyz_x_1[i] * pc_z[i];

        ta1_y_yyzz_y_0[i] = ta1_y_yy_y_0[i] * fe_0 - ta1_y_yy_y_1[i] * fe_0 + ta1_y_yyz_y_0[i] * pa_z[i] - ta1_y_yyz_y_1[i] * pc_z[i];

        ta1_y_yyzz_z_0[i] = ta1_y_zz_z_0[i] * fe_0 - ta1_y_zz_z_1[i] * fe_0 + ta_yzz_z_1[i] + ta1_y_yzz_z_0[i] * pa_y[i] - ta1_y_yzz_z_1[i] * pc_y[i];
    }

    // Set up 84-87 components of targeted buffer : GP

    auto ta1_y_yzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 84);

    auto ta1_y_yzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 85);

    auto ta1_y_yzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 86);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_y_yz_y_0,   \
                             ta1_y_yz_y_1,   \
                             ta1_y_yzz_y_0,  \
                             ta1_y_yzz_y_1,  \
                             ta1_y_yzzz_x_0, \
                             ta1_y_yzzz_y_0, \
                             ta1_y_yzzz_z_0, \
                             ta1_y_zzz_x_0,  \
                             ta1_y_zzz_x_1,  \
                             ta1_y_zzz_z_0,  \
                             ta1_y_zzz_z_1,  \
                             ta_zzz_x_1,     \
                             ta_zzz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzzz_x_0[i] = ta_zzz_x_1[i] + ta1_y_zzz_x_0[i] * pa_y[i] - ta1_y_zzz_x_1[i] * pc_y[i];

        ta1_y_yzzz_y_0[i] = 2.0 * ta1_y_yz_y_0[i] * fe_0 - 2.0 * ta1_y_yz_y_1[i] * fe_0 + ta1_y_yzz_y_0[i] * pa_z[i] - ta1_y_yzz_y_1[i] * pc_z[i];

        ta1_y_yzzz_z_0[i] = ta_zzz_z_1[i] + ta1_y_zzz_z_0[i] * pa_y[i] - ta1_y_zzz_z_1[i] * pc_y[i];
    }

    // Set up 87-90 components of targeted buffer : GP

    auto ta1_y_zzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 87);

    auto ta1_y_zzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 88);

    auto ta1_y_zzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 89);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_y_zz_x_0,   \
                             ta1_y_zz_x_1,   \
                             ta1_y_zz_y_0,   \
                             ta1_y_zz_y_1,   \
                             ta1_y_zz_z_0,   \
                             ta1_y_zz_z_1,   \
                             ta1_y_zzz_0_0,  \
                             ta1_y_zzz_0_1,  \
                             ta1_y_zzz_x_0,  \
                             ta1_y_zzz_x_1,  \
                             ta1_y_zzz_y_0,  \
                             ta1_y_zzz_y_1,  \
                             ta1_y_zzz_z_0,  \
                             ta1_y_zzz_z_1,  \
                             ta1_y_zzzz_x_0, \
                             ta1_y_zzzz_y_0, \
                             ta1_y_zzzz_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzzz_x_0[i] = 3.0 * ta1_y_zz_x_0[i] * fe_0 - 3.0 * ta1_y_zz_x_1[i] * fe_0 + ta1_y_zzz_x_0[i] * pa_z[i] - ta1_y_zzz_x_1[i] * pc_z[i];

        ta1_y_zzzz_y_0[i] = 3.0 * ta1_y_zz_y_0[i] * fe_0 - 3.0 * ta1_y_zz_y_1[i] * fe_0 + ta1_y_zzz_y_0[i] * pa_z[i] - ta1_y_zzz_y_1[i] * pc_z[i];

        ta1_y_zzzz_z_0[i] = 3.0 * ta1_y_zz_z_0[i] * fe_0 - 3.0 * ta1_y_zz_z_1[i] * fe_0 + ta1_y_zzz_0_0[i] * fe_0 - ta1_y_zzz_0_1[i] * fe_0 +
                            ta1_y_zzz_z_0[i] * pa_z[i] - ta1_y_zzz_z_1[i] * pc_z[i];
    }

    // Set up 90-93 components of targeted buffer : GP

    auto ta1_z_xxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 90);

    auto ta1_z_xxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 91);

    auto ta1_z_xxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 92);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_z_xx_x_0,   \
                             ta1_z_xx_x_1,   \
                             ta1_z_xx_y_0,   \
                             ta1_z_xx_y_1,   \
                             ta1_z_xx_z_0,   \
                             ta1_z_xx_z_1,   \
                             ta1_z_xxx_0_0,  \
                             ta1_z_xxx_0_1,  \
                             ta1_z_xxx_x_0,  \
                             ta1_z_xxx_x_1,  \
                             ta1_z_xxx_y_0,  \
                             ta1_z_xxx_y_1,  \
                             ta1_z_xxx_z_0,  \
                             ta1_z_xxx_z_1,  \
                             ta1_z_xxxx_x_0, \
                             ta1_z_xxxx_y_0, \
                             ta1_z_xxxx_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxx_x_0[i] = 3.0 * ta1_z_xx_x_0[i] * fe_0 - 3.0 * ta1_z_xx_x_1[i] * fe_0 + ta1_z_xxx_0_0[i] * fe_0 - ta1_z_xxx_0_1[i] * fe_0 +
                            ta1_z_xxx_x_0[i] * pa_x[i] - ta1_z_xxx_x_1[i] * pc_x[i];

        ta1_z_xxxx_y_0[i] = 3.0 * ta1_z_xx_y_0[i] * fe_0 - 3.0 * ta1_z_xx_y_1[i] * fe_0 + ta1_z_xxx_y_0[i] * pa_x[i] - ta1_z_xxx_y_1[i] * pc_x[i];

        ta1_z_xxxx_z_0[i] = 3.0 * ta1_z_xx_z_0[i] * fe_0 - 3.0 * ta1_z_xx_z_1[i] * fe_0 + ta1_z_xxx_z_0[i] * pa_x[i] - ta1_z_xxx_z_1[i] * pc_x[i];
    }

    // Set up 93-96 components of targeted buffer : GP

    auto ta1_z_xxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 93);

    auto ta1_z_xxxy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 94);

    auto ta1_z_xxxy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 95);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_z_xxx_x_0,  \
                             ta1_z_xxx_x_1,  \
                             ta1_z_xxx_z_0,  \
                             ta1_z_xxx_z_1,  \
                             ta1_z_xxxy_x_0, \
                             ta1_z_xxxy_y_0, \
                             ta1_z_xxxy_z_0, \
                             ta1_z_xxy_y_0,  \
                             ta1_z_xxy_y_1,  \
                             ta1_z_xy_y_0,   \
                             ta1_z_xy_y_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxy_x_0[i] = ta1_z_xxx_x_0[i] * pa_y[i] - ta1_z_xxx_x_1[i] * pc_y[i];

        ta1_z_xxxy_y_0[i] = 2.0 * ta1_z_xy_y_0[i] * fe_0 - 2.0 * ta1_z_xy_y_1[i] * fe_0 + ta1_z_xxy_y_0[i] * pa_x[i] - ta1_z_xxy_y_1[i] * pc_x[i];

        ta1_z_xxxy_z_0[i] = ta1_z_xxx_z_0[i] * pa_y[i] - ta1_z_xxx_z_1[i] * pc_y[i];
    }

    // Set up 96-99 components of targeted buffer : GP

    auto ta1_z_xxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 96);

    auto ta1_z_xxxz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 97);

    auto ta1_z_xxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 98);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_z_xxx_x_0,  \
                             ta1_z_xxx_x_1,  \
                             ta1_z_xxx_y_0,  \
                             ta1_z_xxx_y_1,  \
                             ta1_z_xxxz_x_0, \
                             ta1_z_xxxz_y_0, \
                             ta1_z_xxxz_z_0, \
                             ta1_z_xxz_z_0,  \
                             ta1_z_xxz_z_1,  \
                             ta1_z_xz_z_0,   \
                             ta1_z_xz_z_1,   \
                             ta_xxx_x_1,     \
                             ta_xxx_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxz_x_0[i] = ta_xxx_x_1[i] + ta1_z_xxx_x_0[i] * pa_z[i] - ta1_z_xxx_x_1[i] * pc_z[i];

        ta1_z_xxxz_y_0[i] = ta_xxx_y_1[i] + ta1_z_xxx_y_0[i] * pa_z[i] - ta1_z_xxx_y_1[i] * pc_z[i];

        ta1_z_xxxz_z_0[i] = 2.0 * ta1_z_xz_z_0[i] * fe_0 - 2.0 * ta1_z_xz_z_1[i] * fe_0 + ta1_z_xxz_z_0[i] * pa_x[i] - ta1_z_xxz_z_1[i] * pc_x[i];
    }

    // Set up 99-102 components of targeted buffer : GP

    auto ta1_z_xxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 99);

    auto ta1_z_xxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 100);

    auto ta1_z_xxyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 101);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_z_xx_x_0,   \
                             ta1_z_xx_x_1,   \
                             ta1_z_xxy_x_0,  \
                             ta1_z_xxy_x_1,  \
                             ta1_z_xxyy_x_0, \
                             ta1_z_xxyy_y_0, \
                             ta1_z_xxyy_z_0, \
                             ta1_z_xyy_y_0,  \
                             ta1_z_xyy_y_1,  \
                             ta1_z_xyy_z_0,  \
                             ta1_z_xyy_z_1,  \
                             ta1_z_yy_y_0,   \
                             ta1_z_yy_y_1,   \
                             ta1_z_yy_z_0,   \
                             ta1_z_yy_z_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyy_x_0[i] = ta1_z_xx_x_0[i] * fe_0 - ta1_z_xx_x_1[i] * fe_0 + ta1_z_xxy_x_0[i] * pa_y[i] - ta1_z_xxy_x_1[i] * pc_y[i];

        ta1_z_xxyy_y_0[i] = ta1_z_yy_y_0[i] * fe_0 - ta1_z_yy_y_1[i] * fe_0 + ta1_z_xyy_y_0[i] * pa_x[i] - ta1_z_xyy_y_1[i] * pc_x[i];

        ta1_z_xxyy_z_0[i] = ta1_z_yy_z_0[i] * fe_0 - ta1_z_yy_z_1[i] * fe_0 + ta1_z_xyy_z_0[i] * pa_x[i] - ta1_z_xyy_z_1[i] * pc_x[i];
    }

    // Set up 102-105 components of targeted buffer : GP

    auto ta1_z_xxyz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 102);

    auto ta1_z_xxyz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 103);

    auto ta1_z_xxyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 104);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_z_xxy_y_0,  \
                             ta1_z_xxy_y_1,  \
                             ta1_z_xxyz_x_0, \
                             ta1_z_xxyz_y_0, \
                             ta1_z_xxyz_z_0, \
                             ta1_z_xxz_x_0,  \
                             ta1_z_xxz_x_1,  \
                             ta1_z_xxz_z_0,  \
                             ta1_z_xxz_z_1,  \
                             ta_xxy_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_z_xxyz_x_0[i] = ta1_z_xxz_x_0[i] * pa_y[i] - ta1_z_xxz_x_1[i] * pc_y[i];

        ta1_z_xxyz_y_0[i] = ta_xxy_y_1[i] + ta1_z_xxy_y_0[i] * pa_z[i] - ta1_z_xxy_y_1[i] * pc_z[i];

        ta1_z_xxyz_z_0[i] = ta1_z_xxz_z_0[i] * pa_y[i] - ta1_z_xxz_z_1[i] * pc_y[i];
    }

    // Set up 105-108 components of targeted buffer : GP

    auto ta1_z_xxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 105);

    auto ta1_z_xxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 106);

    auto ta1_z_xxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 107);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_z_xx_x_0,   \
                             ta1_z_xx_x_1,   \
                             ta1_z_xxz_x_0,  \
                             ta1_z_xxz_x_1,  \
                             ta1_z_xxzz_x_0, \
                             ta1_z_xxzz_y_0, \
                             ta1_z_xxzz_z_0, \
                             ta1_z_xzz_y_0,  \
                             ta1_z_xzz_y_1,  \
                             ta1_z_xzz_z_0,  \
                             ta1_z_xzz_z_1,  \
                             ta1_z_zz_y_0,   \
                             ta1_z_zz_y_1,   \
                             ta1_z_zz_z_0,   \
                             ta1_z_zz_z_1,   \
                             ta_xxz_x_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxzz_x_0[i] = ta1_z_xx_x_0[i] * fe_0 - ta1_z_xx_x_1[i] * fe_0 + ta_xxz_x_1[i] + ta1_z_xxz_x_0[i] * pa_z[i] - ta1_z_xxz_x_1[i] * pc_z[i];

        ta1_z_xxzz_y_0[i] = ta1_z_zz_y_0[i] * fe_0 - ta1_z_zz_y_1[i] * fe_0 + ta1_z_xzz_y_0[i] * pa_x[i] - ta1_z_xzz_y_1[i] * pc_x[i];

        ta1_z_xxzz_z_0[i] = ta1_z_zz_z_0[i] * fe_0 - ta1_z_zz_z_1[i] * fe_0 + ta1_z_xzz_z_0[i] * pa_x[i] - ta1_z_xzz_z_1[i] * pc_x[i];
    }

    // Set up 108-111 components of targeted buffer : GP

    auto ta1_z_xyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 108);

    auto ta1_z_xyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 109);

    auto ta1_z_xyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 110);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_z_xyyy_x_0, \
                             ta1_z_xyyy_y_0, \
                             ta1_z_xyyy_z_0, \
                             ta1_z_yyy_0_0,  \
                             ta1_z_yyy_0_1,  \
                             ta1_z_yyy_x_0,  \
                             ta1_z_yyy_x_1,  \
                             ta1_z_yyy_y_0,  \
                             ta1_z_yyy_y_1,  \
                             ta1_z_yyy_z_0,  \
                             ta1_z_yyy_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyy_x_0[i] = ta1_z_yyy_0_0[i] * fe_0 - ta1_z_yyy_0_1[i] * fe_0 + ta1_z_yyy_x_0[i] * pa_x[i] - ta1_z_yyy_x_1[i] * pc_x[i];

        ta1_z_xyyy_y_0[i] = ta1_z_yyy_y_0[i] * pa_x[i] - ta1_z_yyy_y_1[i] * pc_x[i];

        ta1_z_xyyy_z_0[i] = ta1_z_yyy_z_0[i] * pa_x[i] - ta1_z_yyy_z_1[i] * pc_x[i];
    }

    // Set up 111-114 components of targeted buffer : GP

    auto ta1_z_xyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 111);

    auto ta1_z_xyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 112);

    auto ta1_z_xyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 113);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta1_z_xyy_x_0,  \
                             ta1_z_xyy_x_1,  \
                             ta1_z_xyyz_x_0, \
                             ta1_z_xyyz_y_0, \
                             ta1_z_xyyz_z_0, \
                             ta1_z_yyz_y_0,  \
                             ta1_z_yyz_y_1,  \
                             ta1_z_yyz_z_0,  \
                             ta1_z_yyz_z_1,  \
                             ta_xyy_x_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_z_xyyz_x_0[i] = ta_xyy_x_1[i] + ta1_z_xyy_x_0[i] * pa_z[i] - ta1_z_xyy_x_1[i] * pc_z[i];

        ta1_z_xyyz_y_0[i] = ta1_z_yyz_y_0[i] * pa_x[i] - ta1_z_yyz_y_1[i] * pc_x[i];

        ta1_z_xyyz_z_0[i] = ta1_z_yyz_z_0[i] * pa_x[i] - ta1_z_yyz_z_1[i] * pc_x[i];
    }

    // Set up 114-117 components of targeted buffer : GP

    auto ta1_z_xyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 114);

    auto ta1_z_xyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 115);

    auto ta1_z_xyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 116);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta1_z_xyzz_x_0, \
                             ta1_z_xyzz_y_0, \
                             ta1_z_xyzz_z_0, \
                             ta1_z_xzz_x_0,  \
                             ta1_z_xzz_x_1,  \
                             ta1_z_yzz_y_0,  \
                             ta1_z_yzz_y_1,  \
                             ta1_z_yzz_z_0,  \
                             ta1_z_yzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ta1_z_xyzz_x_0[i] = ta1_z_xzz_x_0[i] * pa_y[i] - ta1_z_xzz_x_1[i] * pc_y[i];

        ta1_z_xyzz_y_0[i] = ta1_z_yzz_y_0[i] * pa_x[i] - ta1_z_yzz_y_1[i] * pc_x[i];

        ta1_z_xyzz_z_0[i] = ta1_z_yzz_z_0[i] * pa_x[i] - ta1_z_yzz_z_1[i] * pc_x[i];
    }

    // Set up 117-120 components of targeted buffer : GP

    auto ta1_z_xzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 117);

    auto ta1_z_xzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 118);

    auto ta1_z_xzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 119);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta1_z_xzzz_x_0, \
                             ta1_z_xzzz_y_0, \
                             ta1_z_xzzz_z_0, \
                             ta1_z_zzz_0_0,  \
                             ta1_z_zzz_0_1,  \
                             ta1_z_zzz_x_0,  \
                             ta1_z_zzz_x_1,  \
                             ta1_z_zzz_y_0,  \
                             ta1_z_zzz_y_1,  \
                             ta1_z_zzz_z_0,  \
                             ta1_z_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzzz_x_0[i] = ta1_z_zzz_0_0[i] * fe_0 - ta1_z_zzz_0_1[i] * fe_0 + ta1_z_zzz_x_0[i] * pa_x[i] - ta1_z_zzz_x_1[i] * pc_x[i];

        ta1_z_xzzz_y_0[i] = ta1_z_zzz_y_0[i] * pa_x[i] - ta1_z_zzz_y_1[i] * pc_x[i];

        ta1_z_xzzz_z_0[i] = ta1_z_zzz_z_0[i] * pa_x[i] - ta1_z_zzz_z_1[i] * pc_x[i];
    }

    // Set up 120-123 components of targeted buffer : GP

    auto ta1_z_yyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 120);

    auto ta1_z_yyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 121);

    auto ta1_z_yyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 122);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_z_yy_x_0,   \
                             ta1_z_yy_x_1,   \
                             ta1_z_yy_y_0,   \
                             ta1_z_yy_y_1,   \
                             ta1_z_yy_z_0,   \
                             ta1_z_yy_z_1,   \
                             ta1_z_yyy_0_0,  \
                             ta1_z_yyy_0_1,  \
                             ta1_z_yyy_x_0,  \
                             ta1_z_yyy_x_1,  \
                             ta1_z_yyy_y_0,  \
                             ta1_z_yyy_y_1,  \
                             ta1_z_yyy_z_0,  \
                             ta1_z_yyy_z_1,  \
                             ta1_z_yyyy_x_0, \
                             ta1_z_yyyy_y_0, \
                             ta1_z_yyyy_z_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyy_x_0[i] = 3.0 * ta1_z_yy_x_0[i] * fe_0 - 3.0 * ta1_z_yy_x_1[i] * fe_0 + ta1_z_yyy_x_0[i] * pa_y[i] - ta1_z_yyy_x_1[i] * pc_y[i];

        ta1_z_yyyy_y_0[i] = 3.0 * ta1_z_yy_y_0[i] * fe_0 - 3.0 * ta1_z_yy_y_1[i] * fe_0 + ta1_z_yyy_0_0[i] * fe_0 - ta1_z_yyy_0_1[i] * fe_0 +
                            ta1_z_yyy_y_0[i] * pa_y[i] - ta1_z_yyy_y_1[i] * pc_y[i];

        ta1_z_yyyy_z_0[i] = 3.0 * ta1_z_yy_z_0[i] * fe_0 - 3.0 * ta1_z_yy_z_1[i] * fe_0 + ta1_z_yyy_z_0[i] * pa_y[i] - ta1_z_yyy_z_1[i] * pc_y[i];
    }

    // Set up 123-126 components of targeted buffer : GP

    auto ta1_z_yyyz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 123);

    auto ta1_z_yyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 124);

    auto ta1_z_yyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 125);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_z_yyy_x_0,  \
                             ta1_z_yyy_x_1,  \
                             ta1_z_yyy_y_0,  \
                             ta1_z_yyy_y_1,  \
                             ta1_z_yyyz_x_0, \
                             ta1_z_yyyz_y_0, \
                             ta1_z_yyyz_z_0, \
                             ta1_z_yyz_z_0,  \
                             ta1_z_yyz_z_1,  \
                             ta1_z_yz_z_0,   \
                             ta1_z_yz_z_1,   \
                             ta_yyy_x_1,     \
                             ta_yyy_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyz_x_0[i] = ta_yyy_x_1[i] + ta1_z_yyy_x_0[i] * pa_z[i] - ta1_z_yyy_x_1[i] * pc_z[i];

        ta1_z_yyyz_y_0[i] = ta_yyy_y_1[i] + ta1_z_yyy_y_0[i] * pa_z[i] - ta1_z_yyy_y_1[i] * pc_z[i];

        ta1_z_yyyz_z_0[i] = 2.0 * ta1_z_yz_z_0[i] * fe_0 - 2.0 * ta1_z_yz_z_1[i] * fe_0 + ta1_z_yyz_z_0[i] * pa_y[i] - ta1_z_yyz_z_1[i] * pc_y[i];
    }

    // Set up 126-129 components of targeted buffer : GP

    auto ta1_z_yyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 126);

    auto ta1_z_yyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 127);

    auto ta1_z_yyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 128);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta1_z_yy_y_0,   \
                             ta1_z_yy_y_1,   \
                             ta1_z_yyz_y_0,  \
                             ta1_z_yyz_y_1,  \
                             ta1_z_yyzz_x_0, \
                             ta1_z_yyzz_y_0, \
                             ta1_z_yyzz_z_0, \
                             ta1_z_yzz_x_0,  \
                             ta1_z_yzz_x_1,  \
                             ta1_z_yzz_z_0,  \
                             ta1_z_yzz_z_1,  \
                             ta1_z_zz_x_0,   \
                             ta1_z_zz_x_1,   \
                             ta1_z_zz_z_0,   \
                             ta1_z_zz_z_1,   \
                             ta_yyz_y_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyzz_x_0[i] = ta1_z_zz_x_0[i] * fe_0 - ta1_z_zz_x_1[i] * fe_0 + ta1_z_yzz_x_0[i] * pa_y[i] - ta1_z_yzz_x_1[i] * pc_y[i];

        ta1_z_yyzz_y_0[i] = ta1_z_yy_y_0[i] * fe_0 - ta1_z_yy_y_1[i] * fe_0 + ta_yyz_y_1[i] + ta1_z_yyz_y_0[i] * pa_z[i] - ta1_z_yyz_y_1[i] * pc_z[i];

        ta1_z_yyzz_z_0[i] = ta1_z_zz_z_0[i] * fe_0 - ta1_z_zz_z_1[i] * fe_0 + ta1_z_yzz_z_0[i] * pa_y[i] - ta1_z_yzz_z_1[i] * pc_y[i];
    }

    // Set up 129-132 components of targeted buffer : GP

    auto ta1_z_yzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 129);

    auto ta1_z_yzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 130);

    auto ta1_z_yzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 131);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta1_z_yzzz_x_0, \
                             ta1_z_yzzz_y_0, \
                             ta1_z_yzzz_z_0, \
                             ta1_z_zzz_0_0,  \
                             ta1_z_zzz_0_1,  \
                             ta1_z_zzz_x_0,  \
                             ta1_z_zzz_x_1,  \
                             ta1_z_zzz_y_0,  \
                             ta1_z_zzz_y_1,  \
                             ta1_z_zzz_z_0,  \
                             ta1_z_zzz_z_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzzz_x_0[i] = ta1_z_zzz_x_0[i] * pa_y[i] - ta1_z_zzz_x_1[i] * pc_y[i];

        ta1_z_yzzz_y_0[i] = ta1_z_zzz_0_0[i] * fe_0 - ta1_z_zzz_0_1[i] * fe_0 + ta1_z_zzz_y_0[i] * pa_y[i] - ta1_z_zzz_y_1[i] * pc_y[i];

        ta1_z_yzzz_z_0[i] = ta1_z_zzz_z_0[i] * pa_y[i] - ta1_z_zzz_z_1[i] * pc_y[i];
    }

    // Set up 132-135 components of targeted buffer : GP

    auto ta1_z_zzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_gp + 132);

    auto ta1_z_zzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_gp + 133);

    auto ta1_z_zzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_gp + 134);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta1_z_zz_x_0,   \
                             ta1_z_zz_x_1,   \
                             ta1_z_zz_y_0,   \
                             ta1_z_zz_y_1,   \
                             ta1_z_zz_z_0,   \
                             ta1_z_zz_z_1,   \
                             ta1_z_zzz_0_0,  \
                             ta1_z_zzz_0_1,  \
                             ta1_z_zzz_x_0,  \
                             ta1_z_zzz_x_1,  \
                             ta1_z_zzz_y_0,  \
                             ta1_z_zzz_y_1,  \
                             ta1_z_zzz_z_0,  \
                             ta1_z_zzz_z_1,  \
                             ta1_z_zzzz_x_0, \
                             ta1_z_zzzz_y_0, \
                             ta1_z_zzzz_z_0, \
                             ta_zzz_x_1,     \
                             ta_zzz_y_1,     \
                             ta_zzz_z_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzzz_x_0[i] =
            3.0 * ta1_z_zz_x_0[i] * fe_0 - 3.0 * ta1_z_zz_x_1[i] * fe_0 + ta_zzz_x_1[i] + ta1_z_zzz_x_0[i] * pa_z[i] - ta1_z_zzz_x_1[i] * pc_z[i];

        ta1_z_zzzz_y_0[i] =
            3.0 * ta1_z_zz_y_0[i] * fe_0 - 3.0 * ta1_z_zz_y_1[i] * fe_0 + ta_zzz_y_1[i] + ta1_z_zzz_y_0[i] * pa_z[i] - ta1_z_zzz_y_1[i] * pc_z[i];

        ta1_z_zzzz_z_0[i] = 3.0 * ta1_z_zz_z_0[i] * fe_0 - 3.0 * ta1_z_zz_z_1[i] * fe_0 + ta1_z_zzz_0_0[i] * fe_0 - ta1_z_zzz_0_1[i] * fe_0 +
                            ta_zzz_z_1[i] + ta1_z_zzz_z_0[i] * pa_z[i] - ta1_z_zzz_z_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
