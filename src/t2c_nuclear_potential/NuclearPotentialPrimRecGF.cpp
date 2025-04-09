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

#include "NuclearPotentialPrimRecGF.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_gf(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_gf,
                               const size_t              idx_npot_0_df,
                               const size_t              idx_npot_1_df,
                               const size_t              idx_npot_0_fd,
                               const size_t              idx_npot_1_fd,
                               const size_t              idx_npot_0_ff,
                               const size_t              idx_npot_1_ff,
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

    // Set up components of auxiliary buffer : DF

    auto ta_xx_xxx_0 = pbuffer.data(idx_npot_0_df);

    auto ta_xx_xxy_0 = pbuffer.data(idx_npot_0_df + 1);

    auto ta_xx_xxz_0 = pbuffer.data(idx_npot_0_df + 2);

    auto ta_xx_xyy_0 = pbuffer.data(idx_npot_0_df + 3);

    auto ta_xx_xyz_0 = pbuffer.data(idx_npot_0_df + 4);

    auto ta_xx_xzz_0 = pbuffer.data(idx_npot_0_df + 5);

    auto ta_xx_yyy_0 = pbuffer.data(idx_npot_0_df + 6);

    auto ta_xx_yyz_0 = pbuffer.data(idx_npot_0_df + 7);

    auto ta_xx_yzz_0 = pbuffer.data(idx_npot_0_df + 8);

    auto ta_xx_zzz_0 = pbuffer.data(idx_npot_0_df + 9);

    auto ta_xy_yyy_0 = pbuffer.data(idx_npot_0_df + 16);

    auto ta_xy_yyz_0 = pbuffer.data(idx_npot_0_df + 17);

    auto ta_xy_yzz_0 = pbuffer.data(idx_npot_0_df + 18);

    auto ta_xz_yyz_0 = pbuffer.data(idx_npot_0_df + 27);

    auto ta_xz_yzz_0 = pbuffer.data(idx_npot_0_df + 28);

    auto ta_xz_zzz_0 = pbuffer.data(idx_npot_0_df + 29);

    auto ta_yy_xxx_0 = pbuffer.data(idx_npot_0_df + 30);

    auto ta_yy_xxy_0 = pbuffer.data(idx_npot_0_df + 31);

    auto ta_yy_xxz_0 = pbuffer.data(idx_npot_0_df + 32);

    auto ta_yy_xyy_0 = pbuffer.data(idx_npot_0_df + 33);

    auto ta_yy_xyz_0 = pbuffer.data(idx_npot_0_df + 34);

    auto ta_yy_xzz_0 = pbuffer.data(idx_npot_0_df + 35);

    auto ta_yy_yyy_0 = pbuffer.data(idx_npot_0_df + 36);

    auto ta_yy_yyz_0 = pbuffer.data(idx_npot_0_df + 37);

    auto ta_yy_yzz_0 = pbuffer.data(idx_npot_0_df + 38);

    auto ta_yy_zzz_0 = pbuffer.data(idx_npot_0_df + 39);

    auto ta_yz_xxz_0 = pbuffer.data(idx_npot_0_df + 42);

    auto ta_yz_xzz_0 = pbuffer.data(idx_npot_0_df + 45);

    auto ta_yz_yyz_0 = pbuffer.data(idx_npot_0_df + 47);

    auto ta_yz_yzz_0 = pbuffer.data(idx_npot_0_df + 48);

    auto ta_yz_zzz_0 = pbuffer.data(idx_npot_0_df + 49);

    auto ta_zz_xxx_0 = pbuffer.data(idx_npot_0_df + 50);

    auto ta_zz_xxy_0 = pbuffer.data(idx_npot_0_df + 51);

    auto ta_zz_xxz_0 = pbuffer.data(idx_npot_0_df + 52);

    auto ta_zz_xyy_0 = pbuffer.data(idx_npot_0_df + 53);

    auto ta_zz_xyz_0 = pbuffer.data(idx_npot_0_df + 54);

    auto ta_zz_xzz_0 = pbuffer.data(idx_npot_0_df + 55);

    auto ta_zz_yyy_0 = pbuffer.data(idx_npot_0_df + 56);

    auto ta_zz_yyz_0 = pbuffer.data(idx_npot_0_df + 57);

    auto ta_zz_yzz_0 = pbuffer.data(idx_npot_0_df + 58);

    auto ta_zz_zzz_0 = pbuffer.data(idx_npot_0_df + 59);

    // Set up components of auxiliary buffer : DF

    auto ta_xx_xxx_1 = pbuffer.data(idx_npot_1_df);

    auto ta_xx_xxy_1 = pbuffer.data(idx_npot_1_df + 1);

    auto ta_xx_xxz_1 = pbuffer.data(idx_npot_1_df + 2);

    auto ta_xx_xyy_1 = pbuffer.data(idx_npot_1_df + 3);

    auto ta_xx_xyz_1 = pbuffer.data(idx_npot_1_df + 4);

    auto ta_xx_xzz_1 = pbuffer.data(idx_npot_1_df + 5);

    auto ta_xx_yyy_1 = pbuffer.data(idx_npot_1_df + 6);

    auto ta_xx_yyz_1 = pbuffer.data(idx_npot_1_df + 7);

    auto ta_xx_yzz_1 = pbuffer.data(idx_npot_1_df + 8);

    auto ta_xx_zzz_1 = pbuffer.data(idx_npot_1_df + 9);

    auto ta_xy_yyy_1 = pbuffer.data(idx_npot_1_df + 16);

    auto ta_xy_yyz_1 = pbuffer.data(idx_npot_1_df + 17);

    auto ta_xy_yzz_1 = pbuffer.data(idx_npot_1_df + 18);

    auto ta_xz_yyz_1 = pbuffer.data(idx_npot_1_df + 27);

    auto ta_xz_yzz_1 = pbuffer.data(idx_npot_1_df + 28);

    auto ta_xz_zzz_1 = pbuffer.data(idx_npot_1_df + 29);

    auto ta_yy_xxx_1 = pbuffer.data(idx_npot_1_df + 30);

    auto ta_yy_xxy_1 = pbuffer.data(idx_npot_1_df + 31);

    auto ta_yy_xxz_1 = pbuffer.data(idx_npot_1_df + 32);

    auto ta_yy_xyy_1 = pbuffer.data(idx_npot_1_df + 33);

    auto ta_yy_xyz_1 = pbuffer.data(idx_npot_1_df + 34);

    auto ta_yy_xzz_1 = pbuffer.data(idx_npot_1_df + 35);

    auto ta_yy_yyy_1 = pbuffer.data(idx_npot_1_df + 36);

    auto ta_yy_yyz_1 = pbuffer.data(idx_npot_1_df + 37);

    auto ta_yy_yzz_1 = pbuffer.data(idx_npot_1_df + 38);

    auto ta_yy_zzz_1 = pbuffer.data(idx_npot_1_df + 39);

    auto ta_yz_xxz_1 = pbuffer.data(idx_npot_1_df + 42);

    auto ta_yz_xzz_1 = pbuffer.data(idx_npot_1_df + 45);

    auto ta_yz_yyz_1 = pbuffer.data(idx_npot_1_df + 47);

    auto ta_yz_yzz_1 = pbuffer.data(idx_npot_1_df + 48);

    auto ta_yz_zzz_1 = pbuffer.data(idx_npot_1_df + 49);

    auto ta_zz_xxx_1 = pbuffer.data(idx_npot_1_df + 50);

    auto ta_zz_xxy_1 = pbuffer.data(idx_npot_1_df + 51);

    auto ta_zz_xxz_1 = pbuffer.data(idx_npot_1_df + 52);

    auto ta_zz_xyy_1 = pbuffer.data(idx_npot_1_df + 53);

    auto ta_zz_xyz_1 = pbuffer.data(idx_npot_1_df + 54);

    auto ta_zz_xzz_1 = pbuffer.data(idx_npot_1_df + 55);

    auto ta_zz_yyy_1 = pbuffer.data(idx_npot_1_df + 56);

    auto ta_zz_yyz_1 = pbuffer.data(idx_npot_1_df + 57);

    auto ta_zz_yzz_1 = pbuffer.data(idx_npot_1_df + 58);

    auto ta_zz_zzz_1 = pbuffer.data(idx_npot_1_df + 59);

    // Set up components of auxiliary buffer : FD

    auto ta_xxx_xx_0 = pbuffer.data(idx_npot_0_fd);

    auto ta_xxx_xy_0 = pbuffer.data(idx_npot_0_fd + 1);

    auto ta_xxx_xz_0 = pbuffer.data(idx_npot_0_fd + 2);

    auto ta_xxx_yy_0 = pbuffer.data(idx_npot_0_fd + 3);

    auto ta_xxx_yz_0 = pbuffer.data(idx_npot_0_fd + 4);

    auto ta_xxx_zz_0 = pbuffer.data(idx_npot_0_fd + 5);

    auto ta_xxz_xz_0 = pbuffer.data(idx_npot_0_fd + 14);

    auto ta_xyy_xy_0 = pbuffer.data(idx_npot_0_fd + 19);

    auto ta_xyy_yy_0 = pbuffer.data(idx_npot_0_fd + 21);

    auto ta_xyy_yz_0 = pbuffer.data(idx_npot_0_fd + 22);

    auto ta_xzz_xz_0 = pbuffer.data(idx_npot_0_fd + 32);

    auto ta_xzz_yz_0 = pbuffer.data(idx_npot_0_fd + 34);

    auto ta_xzz_zz_0 = pbuffer.data(idx_npot_0_fd + 35);

    auto ta_yyy_xx_0 = pbuffer.data(idx_npot_0_fd + 36);

    auto ta_yyy_xy_0 = pbuffer.data(idx_npot_0_fd + 37);

    auto ta_yyy_xz_0 = pbuffer.data(idx_npot_0_fd + 38);

    auto ta_yyy_yy_0 = pbuffer.data(idx_npot_0_fd + 39);

    auto ta_yyy_yz_0 = pbuffer.data(idx_npot_0_fd + 40);

    auto ta_yyy_zz_0 = pbuffer.data(idx_npot_0_fd + 41);

    auto ta_yyz_xz_0 = pbuffer.data(idx_npot_0_fd + 44);

    auto ta_yyz_yz_0 = pbuffer.data(idx_npot_0_fd + 46);

    auto ta_yyz_zz_0 = pbuffer.data(idx_npot_0_fd + 47);

    auto ta_yzz_xy_0 = pbuffer.data(idx_npot_0_fd + 49);

    auto ta_yzz_xz_0 = pbuffer.data(idx_npot_0_fd + 50);

    auto ta_yzz_yy_0 = pbuffer.data(idx_npot_0_fd + 51);

    auto ta_yzz_yz_0 = pbuffer.data(idx_npot_0_fd + 52);

    auto ta_yzz_zz_0 = pbuffer.data(idx_npot_0_fd + 53);

    auto ta_zzz_xx_0 = pbuffer.data(idx_npot_0_fd + 54);

    auto ta_zzz_xy_0 = pbuffer.data(idx_npot_0_fd + 55);

    auto ta_zzz_xz_0 = pbuffer.data(idx_npot_0_fd + 56);

    auto ta_zzz_yy_0 = pbuffer.data(idx_npot_0_fd + 57);

    auto ta_zzz_yz_0 = pbuffer.data(idx_npot_0_fd + 58);

    auto ta_zzz_zz_0 = pbuffer.data(idx_npot_0_fd + 59);

    // Set up components of auxiliary buffer : FD

    auto ta_xxx_xx_1 = pbuffer.data(idx_npot_1_fd);

    auto ta_xxx_xy_1 = pbuffer.data(idx_npot_1_fd + 1);

    auto ta_xxx_xz_1 = pbuffer.data(idx_npot_1_fd + 2);

    auto ta_xxx_yy_1 = pbuffer.data(idx_npot_1_fd + 3);

    auto ta_xxx_yz_1 = pbuffer.data(idx_npot_1_fd + 4);

    auto ta_xxx_zz_1 = pbuffer.data(idx_npot_1_fd + 5);

    auto ta_xxz_xz_1 = pbuffer.data(idx_npot_1_fd + 14);

    auto ta_xyy_xy_1 = pbuffer.data(idx_npot_1_fd + 19);

    auto ta_xyy_yy_1 = pbuffer.data(idx_npot_1_fd + 21);

    auto ta_xyy_yz_1 = pbuffer.data(idx_npot_1_fd + 22);

    auto ta_xzz_xz_1 = pbuffer.data(idx_npot_1_fd + 32);

    auto ta_xzz_yz_1 = pbuffer.data(idx_npot_1_fd + 34);

    auto ta_xzz_zz_1 = pbuffer.data(idx_npot_1_fd + 35);

    auto ta_yyy_xx_1 = pbuffer.data(idx_npot_1_fd + 36);

    auto ta_yyy_xy_1 = pbuffer.data(idx_npot_1_fd + 37);

    auto ta_yyy_xz_1 = pbuffer.data(idx_npot_1_fd + 38);

    auto ta_yyy_yy_1 = pbuffer.data(idx_npot_1_fd + 39);

    auto ta_yyy_yz_1 = pbuffer.data(idx_npot_1_fd + 40);

    auto ta_yyy_zz_1 = pbuffer.data(idx_npot_1_fd + 41);

    auto ta_yyz_xz_1 = pbuffer.data(idx_npot_1_fd + 44);

    auto ta_yyz_yz_1 = pbuffer.data(idx_npot_1_fd + 46);

    auto ta_yyz_zz_1 = pbuffer.data(idx_npot_1_fd + 47);

    auto ta_yzz_xy_1 = pbuffer.data(idx_npot_1_fd + 49);

    auto ta_yzz_xz_1 = pbuffer.data(idx_npot_1_fd + 50);

    auto ta_yzz_yy_1 = pbuffer.data(idx_npot_1_fd + 51);

    auto ta_yzz_yz_1 = pbuffer.data(idx_npot_1_fd + 52);

    auto ta_yzz_zz_1 = pbuffer.data(idx_npot_1_fd + 53);

    auto ta_zzz_xx_1 = pbuffer.data(idx_npot_1_fd + 54);

    auto ta_zzz_xy_1 = pbuffer.data(idx_npot_1_fd + 55);

    auto ta_zzz_xz_1 = pbuffer.data(idx_npot_1_fd + 56);

    auto ta_zzz_yy_1 = pbuffer.data(idx_npot_1_fd + 57);

    auto ta_zzz_yz_1 = pbuffer.data(idx_npot_1_fd + 58);

    auto ta_zzz_zz_1 = pbuffer.data(idx_npot_1_fd + 59);

    // Set up components of auxiliary buffer : FF

    auto ta_xxx_xxx_0 = pbuffer.data(idx_npot_0_ff);

    auto ta_xxx_xxy_0 = pbuffer.data(idx_npot_0_ff + 1);

    auto ta_xxx_xxz_0 = pbuffer.data(idx_npot_0_ff + 2);

    auto ta_xxx_xyy_0 = pbuffer.data(idx_npot_0_ff + 3);

    auto ta_xxx_xyz_0 = pbuffer.data(idx_npot_0_ff + 4);

    auto ta_xxx_xzz_0 = pbuffer.data(idx_npot_0_ff + 5);

    auto ta_xxx_yyy_0 = pbuffer.data(idx_npot_0_ff + 6);

    auto ta_xxx_yyz_0 = pbuffer.data(idx_npot_0_ff + 7);

    auto ta_xxx_yzz_0 = pbuffer.data(idx_npot_0_ff + 8);

    auto ta_xxx_zzz_0 = pbuffer.data(idx_npot_0_ff + 9);

    auto ta_xxy_xxx_0 = pbuffer.data(idx_npot_0_ff + 10);

    auto ta_xxy_xxy_0 = pbuffer.data(idx_npot_0_ff + 11);

    auto ta_xxy_xxz_0 = pbuffer.data(idx_npot_0_ff + 12);

    auto ta_xxy_xyy_0 = pbuffer.data(idx_npot_0_ff + 13);

    auto ta_xxy_xzz_0 = pbuffer.data(idx_npot_0_ff + 15);

    auto ta_xxy_yyy_0 = pbuffer.data(idx_npot_0_ff + 16);

    auto ta_xxy_yyz_0 = pbuffer.data(idx_npot_0_ff + 17);

    auto ta_xxy_yzz_0 = pbuffer.data(idx_npot_0_ff + 18);

    auto ta_xxz_xxx_0 = pbuffer.data(idx_npot_0_ff + 20);

    auto ta_xxz_xxy_0 = pbuffer.data(idx_npot_0_ff + 21);

    auto ta_xxz_xxz_0 = pbuffer.data(idx_npot_0_ff + 22);

    auto ta_xxz_xyy_0 = pbuffer.data(idx_npot_0_ff + 23);

    auto ta_xxz_xyz_0 = pbuffer.data(idx_npot_0_ff + 24);

    auto ta_xxz_xzz_0 = pbuffer.data(idx_npot_0_ff + 25);

    auto ta_xxz_yyz_0 = pbuffer.data(idx_npot_0_ff + 27);

    auto ta_xxz_yzz_0 = pbuffer.data(idx_npot_0_ff + 28);

    auto ta_xxz_zzz_0 = pbuffer.data(idx_npot_0_ff + 29);

    auto ta_xyy_xxx_0 = pbuffer.data(idx_npot_0_ff + 30);

    auto ta_xyy_xxy_0 = pbuffer.data(idx_npot_0_ff + 31);

    auto ta_xyy_xyy_0 = pbuffer.data(idx_npot_0_ff + 33);

    auto ta_xyy_xyz_0 = pbuffer.data(idx_npot_0_ff + 34);

    auto ta_xyy_yyy_0 = pbuffer.data(idx_npot_0_ff + 36);

    auto ta_xyy_yyz_0 = pbuffer.data(idx_npot_0_ff + 37);

    auto ta_xyy_yzz_0 = pbuffer.data(idx_npot_0_ff + 38);

    auto ta_xyy_zzz_0 = pbuffer.data(idx_npot_0_ff + 39);

    auto ta_xyz_yyz_0 = pbuffer.data(idx_npot_0_ff + 47);

    auto ta_xyz_yzz_0 = pbuffer.data(idx_npot_0_ff + 48);

    auto ta_xzz_xxx_0 = pbuffer.data(idx_npot_0_ff + 50);

    auto ta_xzz_xxz_0 = pbuffer.data(idx_npot_0_ff + 52);

    auto ta_xzz_xyz_0 = pbuffer.data(idx_npot_0_ff + 54);

    auto ta_xzz_xzz_0 = pbuffer.data(idx_npot_0_ff + 55);

    auto ta_xzz_yyy_0 = pbuffer.data(idx_npot_0_ff + 56);

    auto ta_xzz_yyz_0 = pbuffer.data(idx_npot_0_ff + 57);

    auto ta_xzz_yzz_0 = pbuffer.data(idx_npot_0_ff + 58);

    auto ta_xzz_zzz_0 = pbuffer.data(idx_npot_0_ff + 59);

    auto ta_yyy_xxx_0 = pbuffer.data(idx_npot_0_ff + 60);

    auto ta_yyy_xxy_0 = pbuffer.data(idx_npot_0_ff + 61);

    auto ta_yyy_xxz_0 = pbuffer.data(idx_npot_0_ff + 62);

    auto ta_yyy_xyy_0 = pbuffer.data(idx_npot_0_ff + 63);

    auto ta_yyy_xyz_0 = pbuffer.data(idx_npot_0_ff + 64);

    auto ta_yyy_xzz_0 = pbuffer.data(idx_npot_0_ff + 65);

    auto ta_yyy_yyy_0 = pbuffer.data(idx_npot_0_ff + 66);

    auto ta_yyy_yyz_0 = pbuffer.data(idx_npot_0_ff + 67);

    auto ta_yyy_yzz_0 = pbuffer.data(idx_npot_0_ff + 68);

    auto ta_yyy_zzz_0 = pbuffer.data(idx_npot_0_ff + 69);

    auto ta_yyz_xxy_0 = pbuffer.data(idx_npot_0_ff + 71);

    auto ta_yyz_xxz_0 = pbuffer.data(idx_npot_0_ff + 72);

    auto ta_yyz_xyy_0 = pbuffer.data(idx_npot_0_ff + 73);

    auto ta_yyz_xyz_0 = pbuffer.data(idx_npot_0_ff + 74);

    auto ta_yyz_xzz_0 = pbuffer.data(idx_npot_0_ff + 75);

    auto ta_yyz_yyy_0 = pbuffer.data(idx_npot_0_ff + 76);

    auto ta_yyz_yyz_0 = pbuffer.data(idx_npot_0_ff + 77);

    auto ta_yyz_yzz_0 = pbuffer.data(idx_npot_0_ff + 78);

    auto ta_yyz_zzz_0 = pbuffer.data(idx_npot_0_ff + 79);

    auto ta_yzz_xxx_0 = pbuffer.data(idx_npot_0_ff + 80);

    auto ta_yzz_xxy_0 = pbuffer.data(idx_npot_0_ff + 81);

    auto ta_yzz_xxz_0 = pbuffer.data(idx_npot_0_ff + 82);

    auto ta_yzz_xyy_0 = pbuffer.data(idx_npot_0_ff + 83);

    auto ta_yzz_xyz_0 = pbuffer.data(idx_npot_0_ff + 84);

    auto ta_yzz_xzz_0 = pbuffer.data(idx_npot_0_ff + 85);

    auto ta_yzz_yyy_0 = pbuffer.data(idx_npot_0_ff + 86);

    auto ta_yzz_yyz_0 = pbuffer.data(idx_npot_0_ff + 87);

    auto ta_yzz_yzz_0 = pbuffer.data(idx_npot_0_ff + 88);

    auto ta_yzz_zzz_0 = pbuffer.data(idx_npot_0_ff + 89);

    auto ta_zzz_xxx_0 = pbuffer.data(idx_npot_0_ff + 90);

    auto ta_zzz_xxy_0 = pbuffer.data(idx_npot_0_ff + 91);

    auto ta_zzz_xxz_0 = pbuffer.data(idx_npot_0_ff + 92);

    auto ta_zzz_xyy_0 = pbuffer.data(idx_npot_0_ff + 93);

    auto ta_zzz_xyz_0 = pbuffer.data(idx_npot_0_ff + 94);

    auto ta_zzz_xzz_0 = pbuffer.data(idx_npot_0_ff + 95);

    auto ta_zzz_yyy_0 = pbuffer.data(idx_npot_0_ff + 96);

    auto ta_zzz_yyz_0 = pbuffer.data(idx_npot_0_ff + 97);

    auto ta_zzz_yzz_0 = pbuffer.data(idx_npot_0_ff + 98);

    auto ta_zzz_zzz_0 = pbuffer.data(idx_npot_0_ff + 99);

    // Set up components of auxiliary buffer : FF

    auto ta_xxx_xxx_1 = pbuffer.data(idx_npot_1_ff);

    auto ta_xxx_xxy_1 = pbuffer.data(idx_npot_1_ff + 1);

    auto ta_xxx_xxz_1 = pbuffer.data(idx_npot_1_ff + 2);

    auto ta_xxx_xyy_1 = pbuffer.data(idx_npot_1_ff + 3);

    auto ta_xxx_xyz_1 = pbuffer.data(idx_npot_1_ff + 4);

    auto ta_xxx_xzz_1 = pbuffer.data(idx_npot_1_ff + 5);

    auto ta_xxx_yyy_1 = pbuffer.data(idx_npot_1_ff + 6);

    auto ta_xxx_yyz_1 = pbuffer.data(idx_npot_1_ff + 7);

    auto ta_xxx_yzz_1 = pbuffer.data(idx_npot_1_ff + 8);

    auto ta_xxx_zzz_1 = pbuffer.data(idx_npot_1_ff + 9);

    auto ta_xxy_xxx_1 = pbuffer.data(idx_npot_1_ff + 10);

    auto ta_xxy_xxy_1 = pbuffer.data(idx_npot_1_ff + 11);

    auto ta_xxy_xxz_1 = pbuffer.data(idx_npot_1_ff + 12);

    auto ta_xxy_xyy_1 = pbuffer.data(idx_npot_1_ff + 13);

    auto ta_xxy_xzz_1 = pbuffer.data(idx_npot_1_ff + 15);

    auto ta_xxy_yyy_1 = pbuffer.data(idx_npot_1_ff + 16);

    auto ta_xxy_yyz_1 = pbuffer.data(idx_npot_1_ff + 17);

    auto ta_xxy_yzz_1 = pbuffer.data(idx_npot_1_ff + 18);

    auto ta_xxz_xxx_1 = pbuffer.data(idx_npot_1_ff + 20);

    auto ta_xxz_xxy_1 = pbuffer.data(idx_npot_1_ff + 21);

    auto ta_xxz_xxz_1 = pbuffer.data(idx_npot_1_ff + 22);

    auto ta_xxz_xyy_1 = pbuffer.data(idx_npot_1_ff + 23);

    auto ta_xxz_xyz_1 = pbuffer.data(idx_npot_1_ff + 24);

    auto ta_xxz_xzz_1 = pbuffer.data(idx_npot_1_ff + 25);

    auto ta_xxz_yyz_1 = pbuffer.data(idx_npot_1_ff + 27);

    auto ta_xxz_yzz_1 = pbuffer.data(idx_npot_1_ff + 28);

    auto ta_xxz_zzz_1 = pbuffer.data(idx_npot_1_ff + 29);

    auto ta_xyy_xxx_1 = pbuffer.data(idx_npot_1_ff + 30);

    auto ta_xyy_xxy_1 = pbuffer.data(idx_npot_1_ff + 31);

    auto ta_xyy_xyy_1 = pbuffer.data(idx_npot_1_ff + 33);

    auto ta_xyy_xyz_1 = pbuffer.data(idx_npot_1_ff + 34);

    auto ta_xyy_yyy_1 = pbuffer.data(idx_npot_1_ff + 36);

    auto ta_xyy_yyz_1 = pbuffer.data(idx_npot_1_ff + 37);

    auto ta_xyy_yzz_1 = pbuffer.data(idx_npot_1_ff + 38);

    auto ta_xyy_zzz_1 = pbuffer.data(idx_npot_1_ff + 39);

    auto ta_xyz_yyz_1 = pbuffer.data(idx_npot_1_ff + 47);

    auto ta_xyz_yzz_1 = pbuffer.data(idx_npot_1_ff + 48);

    auto ta_xzz_xxx_1 = pbuffer.data(idx_npot_1_ff + 50);

    auto ta_xzz_xxz_1 = pbuffer.data(idx_npot_1_ff + 52);

    auto ta_xzz_xyz_1 = pbuffer.data(idx_npot_1_ff + 54);

    auto ta_xzz_xzz_1 = pbuffer.data(idx_npot_1_ff + 55);

    auto ta_xzz_yyy_1 = pbuffer.data(idx_npot_1_ff + 56);

    auto ta_xzz_yyz_1 = pbuffer.data(idx_npot_1_ff + 57);

    auto ta_xzz_yzz_1 = pbuffer.data(idx_npot_1_ff + 58);

    auto ta_xzz_zzz_1 = pbuffer.data(idx_npot_1_ff + 59);

    auto ta_yyy_xxx_1 = pbuffer.data(idx_npot_1_ff + 60);

    auto ta_yyy_xxy_1 = pbuffer.data(idx_npot_1_ff + 61);

    auto ta_yyy_xxz_1 = pbuffer.data(idx_npot_1_ff + 62);

    auto ta_yyy_xyy_1 = pbuffer.data(idx_npot_1_ff + 63);

    auto ta_yyy_xyz_1 = pbuffer.data(idx_npot_1_ff + 64);

    auto ta_yyy_xzz_1 = pbuffer.data(idx_npot_1_ff + 65);

    auto ta_yyy_yyy_1 = pbuffer.data(idx_npot_1_ff + 66);

    auto ta_yyy_yyz_1 = pbuffer.data(idx_npot_1_ff + 67);

    auto ta_yyy_yzz_1 = pbuffer.data(idx_npot_1_ff + 68);

    auto ta_yyy_zzz_1 = pbuffer.data(idx_npot_1_ff + 69);

    auto ta_yyz_xxy_1 = pbuffer.data(idx_npot_1_ff + 71);

    auto ta_yyz_xxz_1 = pbuffer.data(idx_npot_1_ff + 72);

    auto ta_yyz_xyy_1 = pbuffer.data(idx_npot_1_ff + 73);

    auto ta_yyz_xyz_1 = pbuffer.data(idx_npot_1_ff + 74);

    auto ta_yyz_xzz_1 = pbuffer.data(idx_npot_1_ff + 75);

    auto ta_yyz_yyy_1 = pbuffer.data(idx_npot_1_ff + 76);

    auto ta_yyz_yyz_1 = pbuffer.data(idx_npot_1_ff + 77);

    auto ta_yyz_yzz_1 = pbuffer.data(idx_npot_1_ff + 78);

    auto ta_yyz_zzz_1 = pbuffer.data(idx_npot_1_ff + 79);

    auto ta_yzz_xxx_1 = pbuffer.data(idx_npot_1_ff + 80);

    auto ta_yzz_xxy_1 = pbuffer.data(idx_npot_1_ff + 81);

    auto ta_yzz_xxz_1 = pbuffer.data(idx_npot_1_ff + 82);

    auto ta_yzz_xyy_1 = pbuffer.data(idx_npot_1_ff + 83);

    auto ta_yzz_xyz_1 = pbuffer.data(idx_npot_1_ff + 84);

    auto ta_yzz_xzz_1 = pbuffer.data(idx_npot_1_ff + 85);

    auto ta_yzz_yyy_1 = pbuffer.data(idx_npot_1_ff + 86);

    auto ta_yzz_yyz_1 = pbuffer.data(idx_npot_1_ff + 87);

    auto ta_yzz_yzz_1 = pbuffer.data(idx_npot_1_ff + 88);

    auto ta_yzz_zzz_1 = pbuffer.data(idx_npot_1_ff + 89);

    auto ta_zzz_xxx_1 = pbuffer.data(idx_npot_1_ff + 90);

    auto ta_zzz_xxy_1 = pbuffer.data(idx_npot_1_ff + 91);

    auto ta_zzz_xxz_1 = pbuffer.data(idx_npot_1_ff + 92);

    auto ta_zzz_xyy_1 = pbuffer.data(idx_npot_1_ff + 93);

    auto ta_zzz_xyz_1 = pbuffer.data(idx_npot_1_ff + 94);

    auto ta_zzz_xzz_1 = pbuffer.data(idx_npot_1_ff + 95);

    auto ta_zzz_yyy_1 = pbuffer.data(idx_npot_1_ff + 96);

    auto ta_zzz_yyz_1 = pbuffer.data(idx_npot_1_ff + 97);

    auto ta_zzz_yzz_1 = pbuffer.data(idx_npot_1_ff + 98);

    auto ta_zzz_zzz_1 = pbuffer.data(idx_npot_1_ff + 99);

    // Set up 0-10 components of targeted buffer : GF

    auto ta_xxxx_xxx_0 = pbuffer.data(idx_npot_0_gf);

    auto ta_xxxx_xxy_0 = pbuffer.data(idx_npot_0_gf + 1);

    auto ta_xxxx_xxz_0 = pbuffer.data(idx_npot_0_gf + 2);

    auto ta_xxxx_xyy_0 = pbuffer.data(idx_npot_0_gf + 3);

    auto ta_xxxx_xyz_0 = pbuffer.data(idx_npot_0_gf + 4);

    auto ta_xxxx_xzz_0 = pbuffer.data(idx_npot_0_gf + 5);

    auto ta_xxxx_yyy_0 = pbuffer.data(idx_npot_0_gf + 6);

    auto ta_xxxx_yyz_0 = pbuffer.data(idx_npot_0_gf + 7);

    auto ta_xxxx_yzz_0 = pbuffer.data(idx_npot_0_gf + 8);

    auto ta_xxxx_zzz_0 = pbuffer.data(idx_npot_0_gf + 9);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta_xx_xxx_0,   \
                             ta_xx_xxx_1,   \
                             ta_xx_xxy_0,   \
                             ta_xx_xxy_1,   \
                             ta_xx_xxz_0,   \
                             ta_xx_xxz_1,   \
                             ta_xx_xyy_0,   \
                             ta_xx_xyy_1,   \
                             ta_xx_xyz_0,   \
                             ta_xx_xyz_1,   \
                             ta_xx_xzz_0,   \
                             ta_xx_xzz_1,   \
                             ta_xx_yyy_0,   \
                             ta_xx_yyy_1,   \
                             ta_xx_yyz_0,   \
                             ta_xx_yyz_1,   \
                             ta_xx_yzz_0,   \
                             ta_xx_yzz_1,   \
                             ta_xx_zzz_0,   \
                             ta_xx_zzz_1,   \
                             ta_xxx_xx_0,   \
                             ta_xxx_xx_1,   \
                             ta_xxx_xxx_0,  \
                             ta_xxx_xxx_1,  \
                             ta_xxx_xxy_0,  \
                             ta_xxx_xxy_1,  \
                             ta_xxx_xxz_0,  \
                             ta_xxx_xxz_1,  \
                             ta_xxx_xy_0,   \
                             ta_xxx_xy_1,   \
                             ta_xxx_xyy_0,  \
                             ta_xxx_xyy_1,  \
                             ta_xxx_xyz_0,  \
                             ta_xxx_xyz_1,  \
                             ta_xxx_xz_0,   \
                             ta_xxx_xz_1,   \
                             ta_xxx_xzz_0,  \
                             ta_xxx_xzz_1,  \
                             ta_xxx_yy_0,   \
                             ta_xxx_yy_1,   \
                             ta_xxx_yyy_0,  \
                             ta_xxx_yyy_1,  \
                             ta_xxx_yyz_0,  \
                             ta_xxx_yyz_1,  \
                             ta_xxx_yz_0,   \
                             ta_xxx_yz_1,   \
                             ta_xxx_yzz_0,  \
                             ta_xxx_yzz_1,  \
                             ta_xxx_zz_0,   \
                             ta_xxx_zz_1,   \
                             ta_xxx_zzz_0,  \
                             ta_xxx_zzz_1,  \
                             ta_xxxx_xxx_0, \
                             ta_xxxx_xxy_0, \
                             ta_xxxx_xxz_0, \
                             ta_xxxx_xyy_0, \
                             ta_xxxx_xyz_0, \
                             ta_xxxx_xzz_0, \
                             ta_xxxx_yyy_0, \
                             ta_xxxx_yyz_0, \
                             ta_xxxx_yzz_0, \
                             ta_xxxx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxx_xxx_0[i] = 3.0 * ta_xx_xxx_0[i] * fe_0 - 3.0 * ta_xx_xxx_1[i] * fe_0 + 3.0 * ta_xxx_xx_0[i] * fe_0 - 3.0 * ta_xxx_xx_1[i] * fe_0 +
                           ta_xxx_xxx_0[i] * pa_x[i] - ta_xxx_xxx_1[i] * pc_x[i];

        ta_xxxx_xxy_0[i] = 3.0 * ta_xx_xxy_0[i] * fe_0 - 3.0 * ta_xx_xxy_1[i] * fe_0 + 2.0 * ta_xxx_xy_0[i] * fe_0 - 2.0 * ta_xxx_xy_1[i] * fe_0 +
                           ta_xxx_xxy_0[i] * pa_x[i] - ta_xxx_xxy_1[i] * pc_x[i];

        ta_xxxx_xxz_0[i] = 3.0 * ta_xx_xxz_0[i] * fe_0 - 3.0 * ta_xx_xxz_1[i] * fe_0 + 2.0 * ta_xxx_xz_0[i] * fe_0 - 2.0 * ta_xxx_xz_1[i] * fe_0 +
                           ta_xxx_xxz_0[i] * pa_x[i] - ta_xxx_xxz_1[i] * pc_x[i];

        ta_xxxx_xyy_0[i] = 3.0 * ta_xx_xyy_0[i] * fe_0 - 3.0 * ta_xx_xyy_1[i] * fe_0 + ta_xxx_yy_0[i] * fe_0 - ta_xxx_yy_1[i] * fe_0 +
                           ta_xxx_xyy_0[i] * pa_x[i] - ta_xxx_xyy_1[i] * pc_x[i];

        ta_xxxx_xyz_0[i] = 3.0 * ta_xx_xyz_0[i] * fe_0 - 3.0 * ta_xx_xyz_1[i] * fe_0 + ta_xxx_yz_0[i] * fe_0 - ta_xxx_yz_1[i] * fe_0 +
                           ta_xxx_xyz_0[i] * pa_x[i] - ta_xxx_xyz_1[i] * pc_x[i];

        ta_xxxx_xzz_0[i] = 3.0 * ta_xx_xzz_0[i] * fe_0 - 3.0 * ta_xx_xzz_1[i] * fe_0 + ta_xxx_zz_0[i] * fe_0 - ta_xxx_zz_1[i] * fe_0 +
                           ta_xxx_xzz_0[i] * pa_x[i] - ta_xxx_xzz_1[i] * pc_x[i];

        ta_xxxx_yyy_0[i] = 3.0 * ta_xx_yyy_0[i] * fe_0 - 3.0 * ta_xx_yyy_1[i] * fe_0 + ta_xxx_yyy_0[i] * pa_x[i] - ta_xxx_yyy_1[i] * pc_x[i];

        ta_xxxx_yyz_0[i] = 3.0 * ta_xx_yyz_0[i] * fe_0 - 3.0 * ta_xx_yyz_1[i] * fe_0 + ta_xxx_yyz_0[i] * pa_x[i] - ta_xxx_yyz_1[i] * pc_x[i];

        ta_xxxx_yzz_0[i] = 3.0 * ta_xx_yzz_0[i] * fe_0 - 3.0 * ta_xx_yzz_1[i] * fe_0 + ta_xxx_yzz_0[i] * pa_x[i] - ta_xxx_yzz_1[i] * pc_x[i];

        ta_xxxx_zzz_0[i] = 3.0 * ta_xx_zzz_0[i] * fe_0 - 3.0 * ta_xx_zzz_1[i] * fe_0 + ta_xxx_zzz_0[i] * pa_x[i] - ta_xxx_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : GF

    auto ta_xxxy_xxx_0 = pbuffer.data(idx_npot_0_gf + 10);

    auto ta_xxxy_xxy_0 = pbuffer.data(idx_npot_0_gf + 11);

    auto ta_xxxy_xxz_0 = pbuffer.data(idx_npot_0_gf + 12);

    auto ta_xxxy_xyy_0 = pbuffer.data(idx_npot_0_gf + 13);

    auto ta_xxxy_xyz_0 = pbuffer.data(idx_npot_0_gf + 14);

    auto ta_xxxy_xzz_0 = pbuffer.data(idx_npot_0_gf + 15);

    auto ta_xxxy_yyy_0 = pbuffer.data(idx_npot_0_gf + 16);

    auto ta_xxxy_yyz_0 = pbuffer.data(idx_npot_0_gf + 17);

    auto ta_xxxy_yzz_0 = pbuffer.data(idx_npot_0_gf + 18);

    auto ta_xxxy_zzz_0 = pbuffer.data(idx_npot_0_gf + 19);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xxx_xx_0,   \
                             ta_xxx_xx_1,   \
                             ta_xxx_xxx_0,  \
                             ta_xxx_xxx_1,  \
                             ta_xxx_xxy_0,  \
                             ta_xxx_xxy_1,  \
                             ta_xxx_xxz_0,  \
                             ta_xxx_xxz_1,  \
                             ta_xxx_xy_0,   \
                             ta_xxx_xy_1,   \
                             ta_xxx_xyy_0,  \
                             ta_xxx_xyy_1,  \
                             ta_xxx_xyz_0,  \
                             ta_xxx_xyz_1,  \
                             ta_xxx_xz_0,   \
                             ta_xxx_xz_1,   \
                             ta_xxx_xzz_0,  \
                             ta_xxx_xzz_1,  \
                             ta_xxx_zzz_0,  \
                             ta_xxx_zzz_1,  \
                             ta_xxxy_xxx_0, \
                             ta_xxxy_xxy_0, \
                             ta_xxxy_xxz_0, \
                             ta_xxxy_xyy_0, \
                             ta_xxxy_xyz_0, \
                             ta_xxxy_xzz_0, \
                             ta_xxxy_yyy_0, \
                             ta_xxxy_yyz_0, \
                             ta_xxxy_yzz_0, \
                             ta_xxxy_zzz_0, \
                             ta_xxy_yyy_0,  \
                             ta_xxy_yyy_1,  \
                             ta_xxy_yyz_0,  \
                             ta_xxy_yyz_1,  \
                             ta_xxy_yzz_0,  \
                             ta_xxy_yzz_1,  \
                             ta_xy_yyy_0,   \
                             ta_xy_yyy_1,   \
                             ta_xy_yyz_0,   \
                             ta_xy_yyz_1,   \
                             ta_xy_yzz_0,   \
                             ta_xy_yzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxy_xxx_0[i] = ta_xxx_xxx_0[i] * pa_y[i] - ta_xxx_xxx_1[i] * pc_y[i];

        ta_xxxy_xxy_0[i] = ta_xxx_xx_0[i] * fe_0 - ta_xxx_xx_1[i] * fe_0 + ta_xxx_xxy_0[i] * pa_y[i] - ta_xxx_xxy_1[i] * pc_y[i];

        ta_xxxy_xxz_0[i] = ta_xxx_xxz_0[i] * pa_y[i] - ta_xxx_xxz_1[i] * pc_y[i];

        ta_xxxy_xyy_0[i] = 2.0 * ta_xxx_xy_0[i] * fe_0 - 2.0 * ta_xxx_xy_1[i] * fe_0 + ta_xxx_xyy_0[i] * pa_y[i] - ta_xxx_xyy_1[i] * pc_y[i];

        ta_xxxy_xyz_0[i] = ta_xxx_xz_0[i] * fe_0 - ta_xxx_xz_1[i] * fe_0 + ta_xxx_xyz_0[i] * pa_y[i] - ta_xxx_xyz_1[i] * pc_y[i];

        ta_xxxy_xzz_0[i] = ta_xxx_xzz_0[i] * pa_y[i] - ta_xxx_xzz_1[i] * pc_y[i];

        ta_xxxy_yyy_0[i] = 2.0 * ta_xy_yyy_0[i] * fe_0 - 2.0 * ta_xy_yyy_1[i] * fe_0 + ta_xxy_yyy_0[i] * pa_x[i] - ta_xxy_yyy_1[i] * pc_x[i];

        ta_xxxy_yyz_0[i] = 2.0 * ta_xy_yyz_0[i] * fe_0 - 2.0 * ta_xy_yyz_1[i] * fe_0 + ta_xxy_yyz_0[i] * pa_x[i] - ta_xxy_yyz_1[i] * pc_x[i];

        ta_xxxy_yzz_0[i] = 2.0 * ta_xy_yzz_0[i] * fe_0 - 2.0 * ta_xy_yzz_1[i] * fe_0 + ta_xxy_yzz_0[i] * pa_x[i] - ta_xxy_yzz_1[i] * pc_x[i];

        ta_xxxy_zzz_0[i] = ta_xxx_zzz_0[i] * pa_y[i] - ta_xxx_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : GF

    auto ta_xxxz_xxx_0 = pbuffer.data(idx_npot_0_gf + 20);

    auto ta_xxxz_xxy_0 = pbuffer.data(idx_npot_0_gf + 21);

    auto ta_xxxz_xxz_0 = pbuffer.data(idx_npot_0_gf + 22);

    auto ta_xxxz_xyy_0 = pbuffer.data(idx_npot_0_gf + 23);

    auto ta_xxxz_xyz_0 = pbuffer.data(idx_npot_0_gf + 24);

    auto ta_xxxz_xzz_0 = pbuffer.data(idx_npot_0_gf + 25);

    auto ta_xxxz_yyy_0 = pbuffer.data(idx_npot_0_gf + 26);

    auto ta_xxxz_yyz_0 = pbuffer.data(idx_npot_0_gf + 27);

    auto ta_xxxz_yzz_0 = pbuffer.data(idx_npot_0_gf + 28);

    auto ta_xxxz_zzz_0 = pbuffer.data(idx_npot_0_gf + 29);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta_xxx_xx_0,   \
                             ta_xxx_xx_1,   \
                             ta_xxx_xxx_0,  \
                             ta_xxx_xxx_1,  \
                             ta_xxx_xxy_0,  \
                             ta_xxx_xxy_1,  \
                             ta_xxx_xxz_0,  \
                             ta_xxx_xxz_1,  \
                             ta_xxx_xy_0,   \
                             ta_xxx_xy_1,   \
                             ta_xxx_xyy_0,  \
                             ta_xxx_xyy_1,  \
                             ta_xxx_xyz_0,  \
                             ta_xxx_xyz_1,  \
                             ta_xxx_xz_0,   \
                             ta_xxx_xz_1,   \
                             ta_xxx_xzz_0,  \
                             ta_xxx_xzz_1,  \
                             ta_xxx_yyy_0,  \
                             ta_xxx_yyy_1,  \
                             ta_xxxz_xxx_0, \
                             ta_xxxz_xxy_0, \
                             ta_xxxz_xxz_0, \
                             ta_xxxz_xyy_0, \
                             ta_xxxz_xyz_0, \
                             ta_xxxz_xzz_0, \
                             ta_xxxz_yyy_0, \
                             ta_xxxz_yyz_0, \
                             ta_xxxz_yzz_0, \
                             ta_xxxz_zzz_0, \
                             ta_xxz_yyz_0,  \
                             ta_xxz_yyz_1,  \
                             ta_xxz_yzz_0,  \
                             ta_xxz_yzz_1,  \
                             ta_xxz_zzz_0,  \
                             ta_xxz_zzz_1,  \
                             ta_xz_yyz_0,   \
                             ta_xz_yyz_1,   \
                             ta_xz_yzz_0,   \
                             ta_xz_yzz_1,   \
                             ta_xz_zzz_0,   \
                             ta_xz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxz_xxx_0[i] = ta_xxx_xxx_0[i] * pa_z[i] - ta_xxx_xxx_1[i] * pc_z[i];

        ta_xxxz_xxy_0[i] = ta_xxx_xxy_0[i] * pa_z[i] - ta_xxx_xxy_1[i] * pc_z[i];

        ta_xxxz_xxz_0[i] = ta_xxx_xx_0[i] * fe_0 - ta_xxx_xx_1[i] * fe_0 + ta_xxx_xxz_0[i] * pa_z[i] - ta_xxx_xxz_1[i] * pc_z[i];

        ta_xxxz_xyy_0[i] = ta_xxx_xyy_0[i] * pa_z[i] - ta_xxx_xyy_1[i] * pc_z[i];

        ta_xxxz_xyz_0[i] = ta_xxx_xy_0[i] * fe_0 - ta_xxx_xy_1[i] * fe_0 + ta_xxx_xyz_0[i] * pa_z[i] - ta_xxx_xyz_1[i] * pc_z[i];

        ta_xxxz_xzz_0[i] = 2.0 * ta_xxx_xz_0[i] * fe_0 - 2.0 * ta_xxx_xz_1[i] * fe_0 + ta_xxx_xzz_0[i] * pa_z[i] - ta_xxx_xzz_1[i] * pc_z[i];

        ta_xxxz_yyy_0[i] = ta_xxx_yyy_0[i] * pa_z[i] - ta_xxx_yyy_1[i] * pc_z[i];

        ta_xxxz_yyz_0[i] = 2.0 * ta_xz_yyz_0[i] * fe_0 - 2.0 * ta_xz_yyz_1[i] * fe_0 + ta_xxz_yyz_0[i] * pa_x[i] - ta_xxz_yyz_1[i] * pc_x[i];

        ta_xxxz_yzz_0[i] = 2.0 * ta_xz_yzz_0[i] * fe_0 - 2.0 * ta_xz_yzz_1[i] * fe_0 + ta_xxz_yzz_0[i] * pa_x[i] - ta_xxz_yzz_1[i] * pc_x[i];

        ta_xxxz_zzz_0[i] = 2.0 * ta_xz_zzz_0[i] * fe_0 - 2.0 * ta_xz_zzz_1[i] * fe_0 + ta_xxz_zzz_0[i] * pa_x[i] - ta_xxz_zzz_1[i] * pc_x[i];
    }

    // Set up 30-40 components of targeted buffer : GF

    auto ta_xxyy_xxx_0 = pbuffer.data(idx_npot_0_gf + 30);

    auto ta_xxyy_xxy_0 = pbuffer.data(idx_npot_0_gf + 31);

    auto ta_xxyy_xxz_0 = pbuffer.data(idx_npot_0_gf + 32);

    auto ta_xxyy_xyy_0 = pbuffer.data(idx_npot_0_gf + 33);

    auto ta_xxyy_xyz_0 = pbuffer.data(idx_npot_0_gf + 34);

    auto ta_xxyy_xzz_0 = pbuffer.data(idx_npot_0_gf + 35);

    auto ta_xxyy_yyy_0 = pbuffer.data(idx_npot_0_gf + 36);

    auto ta_xxyy_yyz_0 = pbuffer.data(idx_npot_0_gf + 37);

    auto ta_xxyy_yzz_0 = pbuffer.data(idx_npot_0_gf + 38);

    auto ta_xxyy_zzz_0 = pbuffer.data(idx_npot_0_gf + 39);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xx_xxx_0,   \
                             ta_xx_xxx_1,   \
                             ta_xx_xxz_0,   \
                             ta_xx_xxz_1,   \
                             ta_xx_xzz_0,   \
                             ta_xx_xzz_1,   \
                             ta_xxy_xxx_0,  \
                             ta_xxy_xxx_1,  \
                             ta_xxy_xxz_0,  \
                             ta_xxy_xxz_1,  \
                             ta_xxy_xzz_0,  \
                             ta_xxy_xzz_1,  \
                             ta_xxyy_xxx_0, \
                             ta_xxyy_xxy_0, \
                             ta_xxyy_xxz_0, \
                             ta_xxyy_xyy_0, \
                             ta_xxyy_xyz_0, \
                             ta_xxyy_xzz_0, \
                             ta_xxyy_yyy_0, \
                             ta_xxyy_yyz_0, \
                             ta_xxyy_yzz_0, \
                             ta_xxyy_zzz_0, \
                             ta_xyy_xxy_0,  \
                             ta_xyy_xxy_1,  \
                             ta_xyy_xy_0,   \
                             ta_xyy_xy_1,   \
                             ta_xyy_xyy_0,  \
                             ta_xyy_xyy_1,  \
                             ta_xyy_xyz_0,  \
                             ta_xyy_xyz_1,  \
                             ta_xyy_yy_0,   \
                             ta_xyy_yy_1,   \
                             ta_xyy_yyy_0,  \
                             ta_xyy_yyy_1,  \
                             ta_xyy_yyz_0,  \
                             ta_xyy_yyz_1,  \
                             ta_xyy_yz_0,   \
                             ta_xyy_yz_1,   \
                             ta_xyy_yzz_0,  \
                             ta_xyy_yzz_1,  \
                             ta_xyy_zzz_0,  \
                             ta_xyy_zzz_1,  \
                             ta_yy_xxy_0,   \
                             ta_yy_xxy_1,   \
                             ta_yy_xyy_0,   \
                             ta_yy_xyy_1,   \
                             ta_yy_xyz_0,   \
                             ta_yy_xyz_1,   \
                             ta_yy_yyy_0,   \
                             ta_yy_yyy_1,   \
                             ta_yy_yyz_0,   \
                             ta_yy_yyz_1,   \
                             ta_yy_yzz_0,   \
                             ta_yy_yzz_1,   \
                             ta_yy_zzz_0,   \
                             ta_yy_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyy_xxx_0[i] = ta_xx_xxx_0[i] * fe_0 - ta_xx_xxx_1[i] * fe_0 + ta_xxy_xxx_0[i] * pa_y[i] - ta_xxy_xxx_1[i] * pc_y[i];

        ta_xxyy_xxy_0[i] = ta_yy_xxy_0[i] * fe_0 - ta_yy_xxy_1[i] * fe_0 + 2.0 * ta_xyy_xy_0[i] * fe_0 - 2.0 * ta_xyy_xy_1[i] * fe_0 +
                           ta_xyy_xxy_0[i] * pa_x[i] - ta_xyy_xxy_1[i] * pc_x[i];

        ta_xxyy_xxz_0[i] = ta_xx_xxz_0[i] * fe_0 - ta_xx_xxz_1[i] * fe_0 + ta_xxy_xxz_0[i] * pa_y[i] - ta_xxy_xxz_1[i] * pc_y[i];

        ta_xxyy_xyy_0[i] = ta_yy_xyy_0[i] * fe_0 - ta_yy_xyy_1[i] * fe_0 + ta_xyy_yy_0[i] * fe_0 - ta_xyy_yy_1[i] * fe_0 + ta_xyy_xyy_0[i] * pa_x[i] -
                           ta_xyy_xyy_1[i] * pc_x[i];

        ta_xxyy_xyz_0[i] = ta_yy_xyz_0[i] * fe_0 - ta_yy_xyz_1[i] * fe_0 + ta_xyy_yz_0[i] * fe_0 - ta_xyy_yz_1[i] * fe_0 + ta_xyy_xyz_0[i] * pa_x[i] -
                           ta_xyy_xyz_1[i] * pc_x[i];

        ta_xxyy_xzz_0[i] = ta_xx_xzz_0[i] * fe_0 - ta_xx_xzz_1[i] * fe_0 + ta_xxy_xzz_0[i] * pa_y[i] - ta_xxy_xzz_1[i] * pc_y[i];

        ta_xxyy_yyy_0[i] = ta_yy_yyy_0[i] * fe_0 - ta_yy_yyy_1[i] * fe_0 + ta_xyy_yyy_0[i] * pa_x[i] - ta_xyy_yyy_1[i] * pc_x[i];

        ta_xxyy_yyz_0[i] = ta_yy_yyz_0[i] * fe_0 - ta_yy_yyz_1[i] * fe_0 + ta_xyy_yyz_0[i] * pa_x[i] - ta_xyy_yyz_1[i] * pc_x[i];

        ta_xxyy_yzz_0[i] = ta_yy_yzz_0[i] * fe_0 - ta_yy_yzz_1[i] * fe_0 + ta_xyy_yzz_0[i] * pa_x[i] - ta_xyy_yzz_1[i] * pc_x[i];

        ta_xxyy_zzz_0[i] = ta_yy_zzz_0[i] * fe_0 - ta_yy_zzz_1[i] * fe_0 + ta_xyy_zzz_0[i] * pa_x[i] - ta_xyy_zzz_1[i] * pc_x[i];
    }

    // Set up 40-50 components of targeted buffer : GF

    auto ta_xxyz_xxx_0 = pbuffer.data(idx_npot_0_gf + 40);

    auto ta_xxyz_xxy_0 = pbuffer.data(idx_npot_0_gf + 41);

    auto ta_xxyz_xxz_0 = pbuffer.data(idx_npot_0_gf + 42);

    auto ta_xxyz_xyy_0 = pbuffer.data(idx_npot_0_gf + 43);

    auto ta_xxyz_xyz_0 = pbuffer.data(idx_npot_0_gf + 44);

    auto ta_xxyz_xzz_0 = pbuffer.data(idx_npot_0_gf + 45);

    auto ta_xxyz_yyy_0 = pbuffer.data(idx_npot_0_gf + 46);

    auto ta_xxyz_yyz_0 = pbuffer.data(idx_npot_0_gf + 47);

    auto ta_xxyz_yzz_0 = pbuffer.data(idx_npot_0_gf + 48);

    auto ta_xxyz_zzz_0 = pbuffer.data(idx_npot_0_gf + 49);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pa_z,          \
                             pc_x,          \
                             pc_y,          \
                             pc_z,          \
                             ta_xxy_xxy_0,  \
                             ta_xxy_xxy_1,  \
                             ta_xxy_xyy_0,  \
                             ta_xxy_xyy_1,  \
                             ta_xxy_yyy_0,  \
                             ta_xxy_yyy_1,  \
                             ta_xxyz_xxx_0, \
                             ta_xxyz_xxy_0, \
                             ta_xxyz_xxz_0, \
                             ta_xxyz_xyy_0, \
                             ta_xxyz_xyz_0, \
                             ta_xxyz_xzz_0, \
                             ta_xxyz_yyy_0, \
                             ta_xxyz_yyz_0, \
                             ta_xxyz_yzz_0, \
                             ta_xxyz_zzz_0, \
                             ta_xxz_xxx_0,  \
                             ta_xxz_xxx_1,  \
                             ta_xxz_xxz_0,  \
                             ta_xxz_xxz_1,  \
                             ta_xxz_xyz_0,  \
                             ta_xxz_xyz_1,  \
                             ta_xxz_xz_0,   \
                             ta_xxz_xz_1,   \
                             ta_xxz_xzz_0,  \
                             ta_xxz_xzz_1,  \
                             ta_xxz_zzz_0,  \
                             ta_xxz_zzz_1,  \
                             ta_xyz_yyz_0,  \
                             ta_xyz_yyz_1,  \
                             ta_xyz_yzz_0,  \
                             ta_xyz_yzz_1,  \
                             ta_yz_yyz_0,   \
                             ta_yz_yyz_1,   \
                             ta_yz_yzz_0,   \
                             ta_yz_yzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyz_xxx_0[i] = ta_xxz_xxx_0[i] * pa_y[i] - ta_xxz_xxx_1[i] * pc_y[i];

        ta_xxyz_xxy_0[i] = ta_xxy_xxy_0[i] * pa_z[i] - ta_xxy_xxy_1[i] * pc_z[i];

        ta_xxyz_xxz_0[i] = ta_xxz_xxz_0[i] * pa_y[i] - ta_xxz_xxz_1[i] * pc_y[i];

        ta_xxyz_xyy_0[i] = ta_xxy_xyy_0[i] * pa_z[i] - ta_xxy_xyy_1[i] * pc_z[i];

        ta_xxyz_xyz_0[i] = ta_xxz_xz_0[i] * fe_0 - ta_xxz_xz_1[i] * fe_0 + ta_xxz_xyz_0[i] * pa_y[i] - ta_xxz_xyz_1[i] * pc_y[i];

        ta_xxyz_xzz_0[i] = ta_xxz_xzz_0[i] * pa_y[i] - ta_xxz_xzz_1[i] * pc_y[i];

        ta_xxyz_yyy_0[i] = ta_xxy_yyy_0[i] * pa_z[i] - ta_xxy_yyy_1[i] * pc_z[i];

        ta_xxyz_yyz_0[i] = ta_yz_yyz_0[i] * fe_0 - ta_yz_yyz_1[i] * fe_0 + ta_xyz_yyz_0[i] * pa_x[i] - ta_xyz_yyz_1[i] * pc_x[i];

        ta_xxyz_yzz_0[i] = ta_yz_yzz_0[i] * fe_0 - ta_yz_yzz_1[i] * fe_0 + ta_xyz_yzz_0[i] * pa_x[i] - ta_xyz_yzz_1[i] * pc_x[i];

        ta_xxyz_zzz_0[i] = ta_xxz_zzz_0[i] * pa_y[i] - ta_xxz_zzz_1[i] * pc_y[i];
    }

    // Set up 50-60 components of targeted buffer : GF

    auto ta_xxzz_xxx_0 = pbuffer.data(idx_npot_0_gf + 50);

    auto ta_xxzz_xxy_0 = pbuffer.data(idx_npot_0_gf + 51);

    auto ta_xxzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 52);

    auto ta_xxzz_xyy_0 = pbuffer.data(idx_npot_0_gf + 53);

    auto ta_xxzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 54);

    auto ta_xxzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 55);

    auto ta_xxzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 56);

    auto ta_xxzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 57);

    auto ta_xxzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 58);

    auto ta_xxzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 59);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta_xx_xxx_0,   \
                             ta_xx_xxx_1,   \
                             ta_xx_xxy_0,   \
                             ta_xx_xxy_1,   \
                             ta_xx_xyy_0,   \
                             ta_xx_xyy_1,   \
                             ta_xxz_xxx_0,  \
                             ta_xxz_xxx_1,  \
                             ta_xxz_xxy_0,  \
                             ta_xxz_xxy_1,  \
                             ta_xxz_xyy_0,  \
                             ta_xxz_xyy_1,  \
                             ta_xxzz_xxx_0, \
                             ta_xxzz_xxy_0, \
                             ta_xxzz_xxz_0, \
                             ta_xxzz_xyy_0, \
                             ta_xxzz_xyz_0, \
                             ta_xxzz_xzz_0, \
                             ta_xxzz_yyy_0, \
                             ta_xxzz_yyz_0, \
                             ta_xxzz_yzz_0, \
                             ta_xxzz_zzz_0, \
                             ta_xzz_xxz_0,  \
                             ta_xzz_xxz_1,  \
                             ta_xzz_xyz_0,  \
                             ta_xzz_xyz_1,  \
                             ta_xzz_xz_0,   \
                             ta_xzz_xz_1,   \
                             ta_xzz_xzz_0,  \
                             ta_xzz_xzz_1,  \
                             ta_xzz_yyy_0,  \
                             ta_xzz_yyy_1,  \
                             ta_xzz_yyz_0,  \
                             ta_xzz_yyz_1,  \
                             ta_xzz_yz_0,   \
                             ta_xzz_yz_1,   \
                             ta_xzz_yzz_0,  \
                             ta_xzz_yzz_1,  \
                             ta_xzz_zz_0,   \
                             ta_xzz_zz_1,   \
                             ta_xzz_zzz_0,  \
                             ta_xzz_zzz_1,  \
                             ta_zz_xxz_0,   \
                             ta_zz_xxz_1,   \
                             ta_zz_xyz_0,   \
                             ta_zz_xyz_1,   \
                             ta_zz_xzz_0,   \
                             ta_zz_xzz_1,   \
                             ta_zz_yyy_0,   \
                             ta_zz_yyy_1,   \
                             ta_zz_yyz_0,   \
                             ta_zz_yyz_1,   \
                             ta_zz_yzz_0,   \
                             ta_zz_yzz_1,   \
                             ta_zz_zzz_0,   \
                             ta_zz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzz_xxx_0[i] = ta_xx_xxx_0[i] * fe_0 - ta_xx_xxx_1[i] * fe_0 + ta_xxz_xxx_0[i] * pa_z[i] - ta_xxz_xxx_1[i] * pc_z[i];

        ta_xxzz_xxy_0[i] = ta_xx_xxy_0[i] * fe_0 - ta_xx_xxy_1[i] * fe_0 + ta_xxz_xxy_0[i] * pa_z[i] - ta_xxz_xxy_1[i] * pc_z[i];

        ta_xxzz_xxz_0[i] = ta_zz_xxz_0[i] * fe_0 - ta_zz_xxz_1[i] * fe_0 + 2.0 * ta_xzz_xz_0[i] * fe_0 - 2.0 * ta_xzz_xz_1[i] * fe_0 +
                           ta_xzz_xxz_0[i] * pa_x[i] - ta_xzz_xxz_1[i] * pc_x[i];

        ta_xxzz_xyy_0[i] = ta_xx_xyy_0[i] * fe_0 - ta_xx_xyy_1[i] * fe_0 + ta_xxz_xyy_0[i] * pa_z[i] - ta_xxz_xyy_1[i] * pc_z[i];

        ta_xxzz_xyz_0[i] = ta_zz_xyz_0[i] * fe_0 - ta_zz_xyz_1[i] * fe_0 + ta_xzz_yz_0[i] * fe_0 - ta_xzz_yz_1[i] * fe_0 + ta_xzz_xyz_0[i] * pa_x[i] -
                           ta_xzz_xyz_1[i] * pc_x[i];

        ta_xxzz_xzz_0[i] = ta_zz_xzz_0[i] * fe_0 - ta_zz_xzz_1[i] * fe_0 + ta_xzz_zz_0[i] * fe_0 - ta_xzz_zz_1[i] * fe_0 + ta_xzz_xzz_0[i] * pa_x[i] -
                           ta_xzz_xzz_1[i] * pc_x[i];

        ta_xxzz_yyy_0[i] = ta_zz_yyy_0[i] * fe_0 - ta_zz_yyy_1[i] * fe_0 + ta_xzz_yyy_0[i] * pa_x[i] - ta_xzz_yyy_1[i] * pc_x[i];

        ta_xxzz_yyz_0[i] = ta_zz_yyz_0[i] * fe_0 - ta_zz_yyz_1[i] * fe_0 + ta_xzz_yyz_0[i] * pa_x[i] - ta_xzz_yyz_1[i] * pc_x[i];

        ta_xxzz_yzz_0[i] = ta_zz_yzz_0[i] * fe_0 - ta_zz_yzz_1[i] * fe_0 + ta_xzz_yzz_0[i] * pa_x[i] - ta_xzz_yzz_1[i] * pc_x[i];

        ta_xxzz_zzz_0[i] = ta_zz_zzz_0[i] * fe_0 - ta_zz_zzz_1[i] * fe_0 + ta_xzz_zzz_0[i] * pa_x[i] - ta_xzz_zzz_1[i] * pc_x[i];
    }

    // Set up 60-70 components of targeted buffer : GF

    auto ta_xyyy_xxx_0 = pbuffer.data(idx_npot_0_gf + 60);

    auto ta_xyyy_xxy_0 = pbuffer.data(idx_npot_0_gf + 61);

    auto ta_xyyy_xxz_0 = pbuffer.data(idx_npot_0_gf + 62);

    auto ta_xyyy_xyy_0 = pbuffer.data(idx_npot_0_gf + 63);

    auto ta_xyyy_xyz_0 = pbuffer.data(idx_npot_0_gf + 64);

    auto ta_xyyy_xzz_0 = pbuffer.data(idx_npot_0_gf + 65);

    auto ta_xyyy_yyy_0 = pbuffer.data(idx_npot_0_gf + 66);

    auto ta_xyyy_yyz_0 = pbuffer.data(idx_npot_0_gf + 67);

    auto ta_xyyy_yzz_0 = pbuffer.data(idx_npot_0_gf + 68);

    auto ta_xyyy_zzz_0 = pbuffer.data(idx_npot_0_gf + 69);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta_xyyy_xxx_0, \
                             ta_xyyy_xxy_0, \
                             ta_xyyy_xxz_0, \
                             ta_xyyy_xyy_0, \
                             ta_xyyy_xyz_0, \
                             ta_xyyy_xzz_0, \
                             ta_xyyy_yyy_0, \
                             ta_xyyy_yyz_0, \
                             ta_xyyy_yzz_0, \
                             ta_xyyy_zzz_0, \
                             ta_yyy_xx_0,   \
                             ta_yyy_xx_1,   \
                             ta_yyy_xxx_0,  \
                             ta_yyy_xxx_1,  \
                             ta_yyy_xxy_0,  \
                             ta_yyy_xxy_1,  \
                             ta_yyy_xxz_0,  \
                             ta_yyy_xxz_1,  \
                             ta_yyy_xy_0,   \
                             ta_yyy_xy_1,   \
                             ta_yyy_xyy_0,  \
                             ta_yyy_xyy_1,  \
                             ta_yyy_xyz_0,  \
                             ta_yyy_xyz_1,  \
                             ta_yyy_xz_0,   \
                             ta_yyy_xz_1,   \
                             ta_yyy_xzz_0,  \
                             ta_yyy_xzz_1,  \
                             ta_yyy_yy_0,   \
                             ta_yyy_yy_1,   \
                             ta_yyy_yyy_0,  \
                             ta_yyy_yyy_1,  \
                             ta_yyy_yyz_0,  \
                             ta_yyy_yyz_1,  \
                             ta_yyy_yz_0,   \
                             ta_yyy_yz_1,   \
                             ta_yyy_yzz_0,  \
                             ta_yyy_yzz_1,  \
                             ta_yyy_zz_0,   \
                             ta_yyy_zz_1,   \
                             ta_yyy_zzz_0,  \
                             ta_yyy_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyy_xxx_0[i] = 3.0 * ta_yyy_xx_0[i] * fe_0 - 3.0 * ta_yyy_xx_1[i] * fe_0 + ta_yyy_xxx_0[i] * pa_x[i] - ta_yyy_xxx_1[i] * pc_x[i];

        ta_xyyy_xxy_0[i] = 2.0 * ta_yyy_xy_0[i] * fe_0 - 2.0 * ta_yyy_xy_1[i] * fe_0 + ta_yyy_xxy_0[i] * pa_x[i] - ta_yyy_xxy_1[i] * pc_x[i];

        ta_xyyy_xxz_0[i] = 2.0 * ta_yyy_xz_0[i] * fe_0 - 2.0 * ta_yyy_xz_1[i] * fe_0 + ta_yyy_xxz_0[i] * pa_x[i] - ta_yyy_xxz_1[i] * pc_x[i];

        ta_xyyy_xyy_0[i] = ta_yyy_yy_0[i] * fe_0 - ta_yyy_yy_1[i] * fe_0 + ta_yyy_xyy_0[i] * pa_x[i] - ta_yyy_xyy_1[i] * pc_x[i];

        ta_xyyy_xyz_0[i] = ta_yyy_yz_0[i] * fe_0 - ta_yyy_yz_1[i] * fe_0 + ta_yyy_xyz_0[i] * pa_x[i] - ta_yyy_xyz_1[i] * pc_x[i];

        ta_xyyy_xzz_0[i] = ta_yyy_zz_0[i] * fe_0 - ta_yyy_zz_1[i] * fe_0 + ta_yyy_xzz_0[i] * pa_x[i] - ta_yyy_xzz_1[i] * pc_x[i];

        ta_xyyy_yyy_0[i] = ta_yyy_yyy_0[i] * pa_x[i] - ta_yyy_yyy_1[i] * pc_x[i];

        ta_xyyy_yyz_0[i] = ta_yyy_yyz_0[i] * pa_x[i] - ta_yyy_yyz_1[i] * pc_x[i];

        ta_xyyy_yzz_0[i] = ta_yyy_yzz_0[i] * pa_x[i] - ta_yyy_yzz_1[i] * pc_x[i];

        ta_xyyy_zzz_0[i] = ta_yyy_zzz_0[i] * pa_x[i] - ta_yyy_zzz_1[i] * pc_x[i];
    }

    // Set up 70-80 components of targeted buffer : GF

    auto ta_xyyz_xxx_0 = pbuffer.data(idx_npot_0_gf + 70);

    auto ta_xyyz_xxy_0 = pbuffer.data(idx_npot_0_gf + 71);

    auto ta_xyyz_xxz_0 = pbuffer.data(idx_npot_0_gf + 72);

    auto ta_xyyz_xyy_0 = pbuffer.data(idx_npot_0_gf + 73);

    auto ta_xyyz_xyz_0 = pbuffer.data(idx_npot_0_gf + 74);

    auto ta_xyyz_xzz_0 = pbuffer.data(idx_npot_0_gf + 75);

    auto ta_xyyz_yyy_0 = pbuffer.data(idx_npot_0_gf + 76);

    auto ta_xyyz_yyz_0 = pbuffer.data(idx_npot_0_gf + 77);

    auto ta_xyyz_yzz_0 = pbuffer.data(idx_npot_0_gf + 78);

    auto ta_xyyz_zzz_0 = pbuffer.data(idx_npot_0_gf + 79);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta_xyy_xxx_0,  \
                             ta_xyy_xxx_1,  \
                             ta_xyy_xxy_0,  \
                             ta_xyy_xxy_1,  \
                             ta_xyy_xyy_0,  \
                             ta_xyy_xyy_1,  \
                             ta_xyyz_xxx_0, \
                             ta_xyyz_xxy_0, \
                             ta_xyyz_xxz_0, \
                             ta_xyyz_xyy_0, \
                             ta_xyyz_xyz_0, \
                             ta_xyyz_xzz_0, \
                             ta_xyyz_yyy_0, \
                             ta_xyyz_yyz_0, \
                             ta_xyyz_yzz_0, \
                             ta_xyyz_zzz_0, \
                             ta_yyz_xxz_0,  \
                             ta_yyz_xxz_1,  \
                             ta_yyz_xyz_0,  \
                             ta_yyz_xyz_1,  \
                             ta_yyz_xz_0,   \
                             ta_yyz_xz_1,   \
                             ta_yyz_xzz_0,  \
                             ta_yyz_xzz_1,  \
                             ta_yyz_yyy_0,  \
                             ta_yyz_yyy_1,  \
                             ta_yyz_yyz_0,  \
                             ta_yyz_yyz_1,  \
                             ta_yyz_yz_0,   \
                             ta_yyz_yz_1,   \
                             ta_yyz_yzz_0,  \
                             ta_yyz_yzz_1,  \
                             ta_yyz_zz_0,   \
                             ta_yyz_zz_1,   \
                             ta_yyz_zzz_0,  \
                             ta_yyz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyz_xxx_0[i] = ta_xyy_xxx_0[i] * pa_z[i] - ta_xyy_xxx_1[i] * pc_z[i];

        ta_xyyz_xxy_0[i] = ta_xyy_xxy_0[i] * pa_z[i] - ta_xyy_xxy_1[i] * pc_z[i];

        ta_xyyz_xxz_0[i] = 2.0 * ta_yyz_xz_0[i] * fe_0 - 2.0 * ta_yyz_xz_1[i] * fe_0 + ta_yyz_xxz_0[i] * pa_x[i] - ta_yyz_xxz_1[i] * pc_x[i];

        ta_xyyz_xyy_0[i] = ta_xyy_xyy_0[i] * pa_z[i] - ta_xyy_xyy_1[i] * pc_z[i];

        ta_xyyz_xyz_0[i] = ta_yyz_yz_0[i] * fe_0 - ta_yyz_yz_1[i] * fe_0 + ta_yyz_xyz_0[i] * pa_x[i] - ta_yyz_xyz_1[i] * pc_x[i];

        ta_xyyz_xzz_0[i] = ta_yyz_zz_0[i] * fe_0 - ta_yyz_zz_1[i] * fe_0 + ta_yyz_xzz_0[i] * pa_x[i] - ta_yyz_xzz_1[i] * pc_x[i];

        ta_xyyz_yyy_0[i] = ta_yyz_yyy_0[i] * pa_x[i] - ta_yyz_yyy_1[i] * pc_x[i];

        ta_xyyz_yyz_0[i] = ta_yyz_yyz_0[i] * pa_x[i] - ta_yyz_yyz_1[i] * pc_x[i];

        ta_xyyz_yzz_0[i] = ta_yyz_yzz_0[i] * pa_x[i] - ta_yyz_yzz_1[i] * pc_x[i];

        ta_xyyz_zzz_0[i] = ta_yyz_zzz_0[i] * pa_x[i] - ta_yyz_zzz_1[i] * pc_x[i];
    }

    // Set up 80-90 components of targeted buffer : GF

    auto ta_xyzz_xxx_0 = pbuffer.data(idx_npot_0_gf + 80);

    auto ta_xyzz_xxy_0 = pbuffer.data(idx_npot_0_gf + 81);

    auto ta_xyzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 82);

    auto ta_xyzz_xyy_0 = pbuffer.data(idx_npot_0_gf + 83);

    auto ta_xyzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 84);

    auto ta_xyzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 85);

    auto ta_xyzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 86);

    auto ta_xyzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 87);

    auto ta_xyzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 88);

    auto ta_xyzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 89);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xyzz_xxx_0, \
                             ta_xyzz_xxy_0, \
                             ta_xyzz_xxz_0, \
                             ta_xyzz_xyy_0, \
                             ta_xyzz_xyz_0, \
                             ta_xyzz_xzz_0, \
                             ta_xyzz_yyy_0, \
                             ta_xyzz_yyz_0, \
                             ta_xyzz_yzz_0, \
                             ta_xyzz_zzz_0, \
                             ta_xzz_xxx_0,  \
                             ta_xzz_xxx_1,  \
                             ta_xzz_xxz_0,  \
                             ta_xzz_xxz_1,  \
                             ta_xzz_xzz_0,  \
                             ta_xzz_xzz_1,  \
                             ta_yzz_xxy_0,  \
                             ta_yzz_xxy_1,  \
                             ta_yzz_xy_0,   \
                             ta_yzz_xy_1,   \
                             ta_yzz_xyy_0,  \
                             ta_yzz_xyy_1,  \
                             ta_yzz_xyz_0,  \
                             ta_yzz_xyz_1,  \
                             ta_yzz_yy_0,   \
                             ta_yzz_yy_1,   \
                             ta_yzz_yyy_0,  \
                             ta_yzz_yyy_1,  \
                             ta_yzz_yyz_0,  \
                             ta_yzz_yyz_1,  \
                             ta_yzz_yz_0,   \
                             ta_yzz_yz_1,   \
                             ta_yzz_yzz_0,  \
                             ta_yzz_yzz_1,  \
                             ta_yzz_zzz_0,  \
                             ta_yzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyzz_xxx_0[i] = ta_xzz_xxx_0[i] * pa_y[i] - ta_xzz_xxx_1[i] * pc_y[i];

        ta_xyzz_xxy_0[i] = 2.0 * ta_yzz_xy_0[i] * fe_0 - 2.0 * ta_yzz_xy_1[i] * fe_0 + ta_yzz_xxy_0[i] * pa_x[i] - ta_yzz_xxy_1[i] * pc_x[i];

        ta_xyzz_xxz_0[i] = ta_xzz_xxz_0[i] * pa_y[i] - ta_xzz_xxz_1[i] * pc_y[i];

        ta_xyzz_xyy_0[i] = ta_yzz_yy_0[i] * fe_0 - ta_yzz_yy_1[i] * fe_0 + ta_yzz_xyy_0[i] * pa_x[i] - ta_yzz_xyy_1[i] * pc_x[i];

        ta_xyzz_xyz_0[i] = ta_yzz_yz_0[i] * fe_0 - ta_yzz_yz_1[i] * fe_0 + ta_yzz_xyz_0[i] * pa_x[i] - ta_yzz_xyz_1[i] * pc_x[i];

        ta_xyzz_xzz_0[i] = ta_xzz_xzz_0[i] * pa_y[i] - ta_xzz_xzz_1[i] * pc_y[i];

        ta_xyzz_yyy_0[i] = ta_yzz_yyy_0[i] * pa_x[i] - ta_yzz_yyy_1[i] * pc_x[i];

        ta_xyzz_yyz_0[i] = ta_yzz_yyz_0[i] * pa_x[i] - ta_yzz_yyz_1[i] * pc_x[i];

        ta_xyzz_yzz_0[i] = ta_yzz_yzz_0[i] * pa_x[i] - ta_yzz_yzz_1[i] * pc_x[i];

        ta_xyzz_zzz_0[i] = ta_yzz_zzz_0[i] * pa_x[i] - ta_yzz_zzz_1[i] * pc_x[i];
    }

    // Set up 90-100 components of targeted buffer : GF

    auto ta_xzzz_xxx_0 = pbuffer.data(idx_npot_0_gf + 90);

    auto ta_xzzz_xxy_0 = pbuffer.data(idx_npot_0_gf + 91);

    auto ta_xzzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 92);

    auto ta_xzzz_xyy_0 = pbuffer.data(idx_npot_0_gf + 93);

    auto ta_xzzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 94);

    auto ta_xzzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 95);

    auto ta_xzzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 96);

    auto ta_xzzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 97);

    auto ta_xzzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 98);

    auto ta_xzzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 99);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta_xzzz_xxx_0, \
                             ta_xzzz_xxy_0, \
                             ta_xzzz_xxz_0, \
                             ta_xzzz_xyy_0, \
                             ta_xzzz_xyz_0, \
                             ta_xzzz_xzz_0, \
                             ta_xzzz_yyy_0, \
                             ta_xzzz_yyz_0, \
                             ta_xzzz_yzz_0, \
                             ta_xzzz_zzz_0, \
                             ta_zzz_xx_0,   \
                             ta_zzz_xx_1,   \
                             ta_zzz_xxx_0,  \
                             ta_zzz_xxx_1,  \
                             ta_zzz_xxy_0,  \
                             ta_zzz_xxy_1,  \
                             ta_zzz_xxz_0,  \
                             ta_zzz_xxz_1,  \
                             ta_zzz_xy_0,   \
                             ta_zzz_xy_1,   \
                             ta_zzz_xyy_0,  \
                             ta_zzz_xyy_1,  \
                             ta_zzz_xyz_0,  \
                             ta_zzz_xyz_1,  \
                             ta_zzz_xz_0,   \
                             ta_zzz_xz_1,   \
                             ta_zzz_xzz_0,  \
                             ta_zzz_xzz_1,  \
                             ta_zzz_yy_0,   \
                             ta_zzz_yy_1,   \
                             ta_zzz_yyy_0,  \
                             ta_zzz_yyy_1,  \
                             ta_zzz_yyz_0,  \
                             ta_zzz_yyz_1,  \
                             ta_zzz_yz_0,   \
                             ta_zzz_yz_1,   \
                             ta_zzz_yzz_0,  \
                             ta_zzz_yzz_1,  \
                             ta_zzz_zz_0,   \
                             ta_zzz_zz_1,   \
                             ta_zzz_zzz_0,  \
                             ta_zzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzz_xxx_0[i] = 3.0 * ta_zzz_xx_0[i] * fe_0 - 3.0 * ta_zzz_xx_1[i] * fe_0 + ta_zzz_xxx_0[i] * pa_x[i] - ta_zzz_xxx_1[i] * pc_x[i];

        ta_xzzz_xxy_0[i] = 2.0 * ta_zzz_xy_0[i] * fe_0 - 2.0 * ta_zzz_xy_1[i] * fe_0 + ta_zzz_xxy_0[i] * pa_x[i] - ta_zzz_xxy_1[i] * pc_x[i];

        ta_xzzz_xxz_0[i] = 2.0 * ta_zzz_xz_0[i] * fe_0 - 2.0 * ta_zzz_xz_1[i] * fe_0 + ta_zzz_xxz_0[i] * pa_x[i] - ta_zzz_xxz_1[i] * pc_x[i];

        ta_xzzz_xyy_0[i] = ta_zzz_yy_0[i] * fe_0 - ta_zzz_yy_1[i] * fe_0 + ta_zzz_xyy_0[i] * pa_x[i] - ta_zzz_xyy_1[i] * pc_x[i];

        ta_xzzz_xyz_0[i] = ta_zzz_yz_0[i] * fe_0 - ta_zzz_yz_1[i] * fe_0 + ta_zzz_xyz_0[i] * pa_x[i] - ta_zzz_xyz_1[i] * pc_x[i];

        ta_xzzz_xzz_0[i] = ta_zzz_zz_0[i] * fe_0 - ta_zzz_zz_1[i] * fe_0 + ta_zzz_xzz_0[i] * pa_x[i] - ta_zzz_xzz_1[i] * pc_x[i];

        ta_xzzz_yyy_0[i] = ta_zzz_yyy_0[i] * pa_x[i] - ta_zzz_yyy_1[i] * pc_x[i];

        ta_xzzz_yyz_0[i] = ta_zzz_yyz_0[i] * pa_x[i] - ta_zzz_yyz_1[i] * pc_x[i];

        ta_xzzz_yzz_0[i] = ta_zzz_yzz_0[i] * pa_x[i] - ta_zzz_yzz_1[i] * pc_x[i];

        ta_xzzz_zzz_0[i] = ta_zzz_zzz_0[i] * pa_x[i] - ta_zzz_zzz_1[i] * pc_x[i];
    }

    // Set up 100-110 components of targeted buffer : GF

    auto ta_yyyy_xxx_0 = pbuffer.data(idx_npot_0_gf + 100);

    auto ta_yyyy_xxy_0 = pbuffer.data(idx_npot_0_gf + 101);

    auto ta_yyyy_xxz_0 = pbuffer.data(idx_npot_0_gf + 102);

    auto ta_yyyy_xyy_0 = pbuffer.data(idx_npot_0_gf + 103);

    auto ta_yyyy_xyz_0 = pbuffer.data(idx_npot_0_gf + 104);

    auto ta_yyyy_xzz_0 = pbuffer.data(idx_npot_0_gf + 105);

    auto ta_yyyy_yyy_0 = pbuffer.data(idx_npot_0_gf + 106);

    auto ta_yyyy_yyz_0 = pbuffer.data(idx_npot_0_gf + 107);

    auto ta_yyyy_yzz_0 = pbuffer.data(idx_npot_0_gf + 108);

    auto ta_yyyy_zzz_0 = pbuffer.data(idx_npot_0_gf + 109);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta_yy_xxx_0,   \
                             ta_yy_xxx_1,   \
                             ta_yy_xxy_0,   \
                             ta_yy_xxy_1,   \
                             ta_yy_xxz_0,   \
                             ta_yy_xxz_1,   \
                             ta_yy_xyy_0,   \
                             ta_yy_xyy_1,   \
                             ta_yy_xyz_0,   \
                             ta_yy_xyz_1,   \
                             ta_yy_xzz_0,   \
                             ta_yy_xzz_1,   \
                             ta_yy_yyy_0,   \
                             ta_yy_yyy_1,   \
                             ta_yy_yyz_0,   \
                             ta_yy_yyz_1,   \
                             ta_yy_yzz_0,   \
                             ta_yy_yzz_1,   \
                             ta_yy_zzz_0,   \
                             ta_yy_zzz_1,   \
                             ta_yyy_xx_0,   \
                             ta_yyy_xx_1,   \
                             ta_yyy_xxx_0,  \
                             ta_yyy_xxx_1,  \
                             ta_yyy_xxy_0,  \
                             ta_yyy_xxy_1,  \
                             ta_yyy_xxz_0,  \
                             ta_yyy_xxz_1,  \
                             ta_yyy_xy_0,   \
                             ta_yyy_xy_1,   \
                             ta_yyy_xyy_0,  \
                             ta_yyy_xyy_1,  \
                             ta_yyy_xyz_0,  \
                             ta_yyy_xyz_1,  \
                             ta_yyy_xz_0,   \
                             ta_yyy_xz_1,   \
                             ta_yyy_xzz_0,  \
                             ta_yyy_xzz_1,  \
                             ta_yyy_yy_0,   \
                             ta_yyy_yy_1,   \
                             ta_yyy_yyy_0,  \
                             ta_yyy_yyy_1,  \
                             ta_yyy_yyz_0,  \
                             ta_yyy_yyz_1,  \
                             ta_yyy_yz_0,   \
                             ta_yyy_yz_1,   \
                             ta_yyy_yzz_0,  \
                             ta_yyy_yzz_1,  \
                             ta_yyy_zz_0,   \
                             ta_yyy_zz_1,   \
                             ta_yyy_zzz_0,  \
                             ta_yyy_zzz_1,  \
                             ta_yyyy_xxx_0, \
                             ta_yyyy_xxy_0, \
                             ta_yyyy_xxz_0, \
                             ta_yyyy_xyy_0, \
                             ta_yyyy_xyz_0, \
                             ta_yyyy_xzz_0, \
                             ta_yyyy_yyy_0, \
                             ta_yyyy_yyz_0, \
                             ta_yyyy_yzz_0, \
                             ta_yyyy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyy_xxx_0[i] = 3.0 * ta_yy_xxx_0[i] * fe_0 - 3.0 * ta_yy_xxx_1[i] * fe_0 + ta_yyy_xxx_0[i] * pa_y[i] - ta_yyy_xxx_1[i] * pc_y[i];

        ta_yyyy_xxy_0[i] = 3.0 * ta_yy_xxy_0[i] * fe_0 - 3.0 * ta_yy_xxy_1[i] * fe_0 + ta_yyy_xx_0[i] * fe_0 - ta_yyy_xx_1[i] * fe_0 +
                           ta_yyy_xxy_0[i] * pa_y[i] - ta_yyy_xxy_1[i] * pc_y[i];

        ta_yyyy_xxz_0[i] = 3.0 * ta_yy_xxz_0[i] * fe_0 - 3.0 * ta_yy_xxz_1[i] * fe_0 + ta_yyy_xxz_0[i] * pa_y[i] - ta_yyy_xxz_1[i] * pc_y[i];

        ta_yyyy_xyy_0[i] = 3.0 * ta_yy_xyy_0[i] * fe_0 - 3.0 * ta_yy_xyy_1[i] * fe_0 + 2.0 * ta_yyy_xy_0[i] * fe_0 - 2.0 * ta_yyy_xy_1[i] * fe_0 +
                           ta_yyy_xyy_0[i] * pa_y[i] - ta_yyy_xyy_1[i] * pc_y[i];

        ta_yyyy_xyz_0[i] = 3.0 * ta_yy_xyz_0[i] * fe_0 - 3.0 * ta_yy_xyz_1[i] * fe_0 + ta_yyy_xz_0[i] * fe_0 - ta_yyy_xz_1[i] * fe_0 +
                           ta_yyy_xyz_0[i] * pa_y[i] - ta_yyy_xyz_1[i] * pc_y[i];

        ta_yyyy_xzz_0[i] = 3.0 * ta_yy_xzz_0[i] * fe_0 - 3.0 * ta_yy_xzz_1[i] * fe_0 + ta_yyy_xzz_0[i] * pa_y[i] - ta_yyy_xzz_1[i] * pc_y[i];

        ta_yyyy_yyy_0[i] = 3.0 * ta_yy_yyy_0[i] * fe_0 - 3.0 * ta_yy_yyy_1[i] * fe_0 + 3.0 * ta_yyy_yy_0[i] * fe_0 - 3.0 * ta_yyy_yy_1[i] * fe_0 +
                           ta_yyy_yyy_0[i] * pa_y[i] - ta_yyy_yyy_1[i] * pc_y[i];

        ta_yyyy_yyz_0[i] = 3.0 * ta_yy_yyz_0[i] * fe_0 - 3.0 * ta_yy_yyz_1[i] * fe_0 + 2.0 * ta_yyy_yz_0[i] * fe_0 - 2.0 * ta_yyy_yz_1[i] * fe_0 +
                           ta_yyy_yyz_0[i] * pa_y[i] - ta_yyy_yyz_1[i] * pc_y[i];

        ta_yyyy_yzz_0[i] = 3.0 * ta_yy_yzz_0[i] * fe_0 - 3.0 * ta_yy_yzz_1[i] * fe_0 + ta_yyy_zz_0[i] * fe_0 - ta_yyy_zz_1[i] * fe_0 +
                           ta_yyy_yzz_0[i] * pa_y[i] - ta_yyy_yzz_1[i] * pc_y[i];

        ta_yyyy_zzz_0[i] = 3.0 * ta_yy_zzz_0[i] * fe_0 - 3.0 * ta_yy_zzz_1[i] * fe_0 + ta_yyy_zzz_0[i] * pa_y[i] - ta_yyy_zzz_1[i] * pc_y[i];
    }

    // Set up 110-120 components of targeted buffer : GF

    auto ta_yyyz_xxx_0 = pbuffer.data(idx_npot_0_gf + 110);

    auto ta_yyyz_xxy_0 = pbuffer.data(idx_npot_0_gf + 111);

    auto ta_yyyz_xxz_0 = pbuffer.data(idx_npot_0_gf + 112);

    auto ta_yyyz_xyy_0 = pbuffer.data(idx_npot_0_gf + 113);

    auto ta_yyyz_xyz_0 = pbuffer.data(idx_npot_0_gf + 114);

    auto ta_yyyz_xzz_0 = pbuffer.data(idx_npot_0_gf + 115);

    auto ta_yyyz_yyy_0 = pbuffer.data(idx_npot_0_gf + 116);

    auto ta_yyyz_yyz_0 = pbuffer.data(idx_npot_0_gf + 117);

    auto ta_yyyz_yzz_0 = pbuffer.data(idx_npot_0_gf + 118);

    auto ta_yyyz_zzz_0 = pbuffer.data(idx_npot_0_gf + 119);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta_yyy_xxx_0,  \
                             ta_yyy_xxx_1,  \
                             ta_yyy_xxy_0,  \
                             ta_yyy_xxy_1,  \
                             ta_yyy_xy_0,   \
                             ta_yyy_xy_1,   \
                             ta_yyy_xyy_0,  \
                             ta_yyy_xyy_1,  \
                             ta_yyy_xyz_0,  \
                             ta_yyy_xyz_1,  \
                             ta_yyy_yy_0,   \
                             ta_yyy_yy_1,   \
                             ta_yyy_yyy_0,  \
                             ta_yyy_yyy_1,  \
                             ta_yyy_yyz_0,  \
                             ta_yyy_yyz_1,  \
                             ta_yyy_yz_0,   \
                             ta_yyy_yz_1,   \
                             ta_yyy_yzz_0,  \
                             ta_yyy_yzz_1,  \
                             ta_yyyz_xxx_0, \
                             ta_yyyz_xxy_0, \
                             ta_yyyz_xxz_0, \
                             ta_yyyz_xyy_0, \
                             ta_yyyz_xyz_0, \
                             ta_yyyz_xzz_0, \
                             ta_yyyz_yyy_0, \
                             ta_yyyz_yyz_0, \
                             ta_yyyz_yzz_0, \
                             ta_yyyz_zzz_0, \
                             ta_yyz_xxz_0,  \
                             ta_yyz_xxz_1,  \
                             ta_yyz_xzz_0,  \
                             ta_yyz_xzz_1,  \
                             ta_yyz_zzz_0,  \
                             ta_yyz_zzz_1,  \
                             ta_yz_xxz_0,   \
                             ta_yz_xxz_1,   \
                             ta_yz_xzz_0,   \
                             ta_yz_xzz_1,   \
                             ta_yz_zzz_0,   \
                             ta_yz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyz_xxx_0[i] = ta_yyy_xxx_0[i] * pa_z[i] - ta_yyy_xxx_1[i] * pc_z[i];

        ta_yyyz_xxy_0[i] = ta_yyy_xxy_0[i] * pa_z[i] - ta_yyy_xxy_1[i] * pc_z[i];

        ta_yyyz_xxz_0[i] = 2.0 * ta_yz_xxz_0[i] * fe_0 - 2.0 * ta_yz_xxz_1[i] * fe_0 + ta_yyz_xxz_0[i] * pa_y[i] - ta_yyz_xxz_1[i] * pc_y[i];

        ta_yyyz_xyy_0[i] = ta_yyy_xyy_0[i] * pa_z[i] - ta_yyy_xyy_1[i] * pc_z[i];

        ta_yyyz_xyz_0[i] = ta_yyy_xy_0[i] * fe_0 - ta_yyy_xy_1[i] * fe_0 + ta_yyy_xyz_0[i] * pa_z[i] - ta_yyy_xyz_1[i] * pc_z[i];

        ta_yyyz_xzz_0[i] = 2.0 * ta_yz_xzz_0[i] * fe_0 - 2.0 * ta_yz_xzz_1[i] * fe_0 + ta_yyz_xzz_0[i] * pa_y[i] - ta_yyz_xzz_1[i] * pc_y[i];

        ta_yyyz_yyy_0[i] = ta_yyy_yyy_0[i] * pa_z[i] - ta_yyy_yyy_1[i] * pc_z[i];

        ta_yyyz_yyz_0[i] = ta_yyy_yy_0[i] * fe_0 - ta_yyy_yy_1[i] * fe_0 + ta_yyy_yyz_0[i] * pa_z[i] - ta_yyy_yyz_1[i] * pc_z[i];

        ta_yyyz_yzz_0[i] = 2.0 * ta_yyy_yz_0[i] * fe_0 - 2.0 * ta_yyy_yz_1[i] * fe_0 + ta_yyy_yzz_0[i] * pa_z[i] - ta_yyy_yzz_1[i] * pc_z[i];

        ta_yyyz_zzz_0[i] = 2.0 * ta_yz_zzz_0[i] * fe_0 - 2.0 * ta_yz_zzz_1[i] * fe_0 + ta_yyz_zzz_0[i] * pa_y[i] - ta_yyz_zzz_1[i] * pc_y[i];
    }

    // Set up 120-130 components of targeted buffer : GF

    auto ta_yyzz_xxx_0 = pbuffer.data(idx_npot_0_gf + 120);

    auto ta_yyzz_xxy_0 = pbuffer.data(idx_npot_0_gf + 121);

    auto ta_yyzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 122);

    auto ta_yyzz_xyy_0 = pbuffer.data(idx_npot_0_gf + 123);

    auto ta_yyzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 124);

    auto ta_yyzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 125);

    auto ta_yyzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 126);

    auto ta_yyzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 127);

    auto ta_yyzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 128);

    auto ta_yyzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 129);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta_yy_xxy_0,   \
                             ta_yy_xxy_1,   \
                             ta_yy_xyy_0,   \
                             ta_yy_xyy_1,   \
                             ta_yy_yyy_0,   \
                             ta_yy_yyy_1,   \
                             ta_yyz_xxy_0,  \
                             ta_yyz_xxy_1,  \
                             ta_yyz_xyy_0,  \
                             ta_yyz_xyy_1,  \
                             ta_yyz_yyy_0,  \
                             ta_yyz_yyy_1,  \
                             ta_yyzz_xxx_0, \
                             ta_yyzz_xxy_0, \
                             ta_yyzz_xxz_0, \
                             ta_yyzz_xyy_0, \
                             ta_yyzz_xyz_0, \
                             ta_yyzz_xzz_0, \
                             ta_yyzz_yyy_0, \
                             ta_yyzz_yyz_0, \
                             ta_yyzz_yzz_0, \
                             ta_yyzz_zzz_0, \
                             ta_yzz_xxx_0,  \
                             ta_yzz_xxx_1,  \
                             ta_yzz_xxz_0,  \
                             ta_yzz_xxz_1,  \
                             ta_yzz_xyz_0,  \
                             ta_yzz_xyz_1,  \
                             ta_yzz_xz_0,   \
                             ta_yzz_xz_1,   \
                             ta_yzz_xzz_0,  \
                             ta_yzz_xzz_1,  \
                             ta_yzz_yyz_0,  \
                             ta_yzz_yyz_1,  \
                             ta_yzz_yz_0,   \
                             ta_yzz_yz_1,   \
                             ta_yzz_yzz_0,  \
                             ta_yzz_yzz_1,  \
                             ta_yzz_zz_0,   \
                             ta_yzz_zz_1,   \
                             ta_yzz_zzz_0,  \
                             ta_yzz_zzz_1,  \
                             ta_zz_xxx_0,   \
                             ta_zz_xxx_1,   \
                             ta_zz_xxz_0,   \
                             ta_zz_xxz_1,   \
                             ta_zz_xyz_0,   \
                             ta_zz_xyz_1,   \
                             ta_zz_xzz_0,   \
                             ta_zz_xzz_1,   \
                             ta_zz_yyz_0,   \
                             ta_zz_yyz_1,   \
                             ta_zz_yzz_0,   \
                             ta_zz_yzz_1,   \
                             ta_zz_zzz_0,   \
                             ta_zz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzz_xxx_0[i] = ta_zz_xxx_0[i] * fe_0 - ta_zz_xxx_1[i] * fe_0 + ta_yzz_xxx_0[i] * pa_y[i] - ta_yzz_xxx_1[i] * pc_y[i];

        ta_yyzz_xxy_0[i] = ta_yy_xxy_0[i] * fe_0 - ta_yy_xxy_1[i] * fe_0 + ta_yyz_xxy_0[i] * pa_z[i] - ta_yyz_xxy_1[i] * pc_z[i];

        ta_yyzz_xxz_0[i] = ta_zz_xxz_0[i] * fe_0 - ta_zz_xxz_1[i] * fe_0 + ta_yzz_xxz_0[i] * pa_y[i] - ta_yzz_xxz_1[i] * pc_y[i];

        ta_yyzz_xyy_0[i] = ta_yy_xyy_0[i] * fe_0 - ta_yy_xyy_1[i] * fe_0 + ta_yyz_xyy_0[i] * pa_z[i] - ta_yyz_xyy_1[i] * pc_z[i];

        ta_yyzz_xyz_0[i] = ta_zz_xyz_0[i] * fe_0 - ta_zz_xyz_1[i] * fe_0 + ta_yzz_xz_0[i] * fe_0 - ta_yzz_xz_1[i] * fe_0 + ta_yzz_xyz_0[i] * pa_y[i] -
                           ta_yzz_xyz_1[i] * pc_y[i];

        ta_yyzz_xzz_0[i] = ta_zz_xzz_0[i] * fe_0 - ta_zz_xzz_1[i] * fe_0 + ta_yzz_xzz_0[i] * pa_y[i] - ta_yzz_xzz_1[i] * pc_y[i];

        ta_yyzz_yyy_0[i] = ta_yy_yyy_0[i] * fe_0 - ta_yy_yyy_1[i] * fe_0 + ta_yyz_yyy_0[i] * pa_z[i] - ta_yyz_yyy_1[i] * pc_z[i];

        ta_yyzz_yyz_0[i] = ta_zz_yyz_0[i] * fe_0 - ta_zz_yyz_1[i] * fe_0 + 2.0 * ta_yzz_yz_0[i] * fe_0 - 2.0 * ta_yzz_yz_1[i] * fe_0 +
                           ta_yzz_yyz_0[i] * pa_y[i] - ta_yzz_yyz_1[i] * pc_y[i];

        ta_yyzz_yzz_0[i] = ta_zz_yzz_0[i] * fe_0 - ta_zz_yzz_1[i] * fe_0 + ta_yzz_zz_0[i] * fe_0 - ta_yzz_zz_1[i] * fe_0 + ta_yzz_yzz_0[i] * pa_y[i] -
                           ta_yzz_yzz_1[i] * pc_y[i];

        ta_yyzz_zzz_0[i] = ta_zz_zzz_0[i] * fe_0 - ta_zz_zzz_1[i] * fe_0 + ta_yzz_zzz_0[i] * pa_y[i] - ta_yzz_zzz_1[i] * pc_y[i];
    }

    // Set up 130-140 components of targeted buffer : GF

    auto ta_yzzz_xxx_0 = pbuffer.data(idx_npot_0_gf + 130);

    auto ta_yzzz_xxy_0 = pbuffer.data(idx_npot_0_gf + 131);

    auto ta_yzzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 132);

    auto ta_yzzz_xyy_0 = pbuffer.data(idx_npot_0_gf + 133);

    auto ta_yzzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 134);

    auto ta_yzzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 135);

    auto ta_yzzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 136);

    auto ta_yzzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 137);

    auto ta_yzzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 138);

    auto ta_yzzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 139);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta_yzzz_xxx_0, \
                             ta_yzzz_xxy_0, \
                             ta_yzzz_xxz_0, \
                             ta_yzzz_xyy_0, \
                             ta_yzzz_xyz_0, \
                             ta_yzzz_xzz_0, \
                             ta_yzzz_yyy_0, \
                             ta_yzzz_yyz_0, \
                             ta_yzzz_yzz_0, \
                             ta_yzzz_zzz_0, \
                             ta_zzz_xx_0,   \
                             ta_zzz_xx_1,   \
                             ta_zzz_xxx_0,  \
                             ta_zzz_xxx_1,  \
                             ta_zzz_xxy_0,  \
                             ta_zzz_xxy_1,  \
                             ta_zzz_xxz_0,  \
                             ta_zzz_xxz_1,  \
                             ta_zzz_xy_0,   \
                             ta_zzz_xy_1,   \
                             ta_zzz_xyy_0,  \
                             ta_zzz_xyy_1,  \
                             ta_zzz_xyz_0,  \
                             ta_zzz_xyz_1,  \
                             ta_zzz_xz_0,   \
                             ta_zzz_xz_1,   \
                             ta_zzz_xzz_0,  \
                             ta_zzz_xzz_1,  \
                             ta_zzz_yy_0,   \
                             ta_zzz_yy_1,   \
                             ta_zzz_yyy_0,  \
                             ta_zzz_yyy_1,  \
                             ta_zzz_yyz_0,  \
                             ta_zzz_yyz_1,  \
                             ta_zzz_yz_0,   \
                             ta_zzz_yz_1,   \
                             ta_zzz_yzz_0,  \
                             ta_zzz_yzz_1,  \
                             ta_zzz_zz_0,   \
                             ta_zzz_zz_1,   \
                             ta_zzz_zzz_0,  \
                             ta_zzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzz_xxx_0[i] = ta_zzz_xxx_0[i] * pa_y[i] - ta_zzz_xxx_1[i] * pc_y[i];

        ta_yzzz_xxy_0[i] = ta_zzz_xx_0[i] * fe_0 - ta_zzz_xx_1[i] * fe_0 + ta_zzz_xxy_0[i] * pa_y[i] - ta_zzz_xxy_1[i] * pc_y[i];

        ta_yzzz_xxz_0[i] = ta_zzz_xxz_0[i] * pa_y[i] - ta_zzz_xxz_1[i] * pc_y[i];

        ta_yzzz_xyy_0[i] = 2.0 * ta_zzz_xy_0[i] * fe_0 - 2.0 * ta_zzz_xy_1[i] * fe_0 + ta_zzz_xyy_0[i] * pa_y[i] - ta_zzz_xyy_1[i] * pc_y[i];

        ta_yzzz_xyz_0[i] = ta_zzz_xz_0[i] * fe_0 - ta_zzz_xz_1[i] * fe_0 + ta_zzz_xyz_0[i] * pa_y[i] - ta_zzz_xyz_1[i] * pc_y[i];

        ta_yzzz_xzz_0[i] = ta_zzz_xzz_0[i] * pa_y[i] - ta_zzz_xzz_1[i] * pc_y[i];

        ta_yzzz_yyy_0[i] = 3.0 * ta_zzz_yy_0[i] * fe_0 - 3.0 * ta_zzz_yy_1[i] * fe_0 + ta_zzz_yyy_0[i] * pa_y[i] - ta_zzz_yyy_1[i] * pc_y[i];

        ta_yzzz_yyz_0[i] = 2.0 * ta_zzz_yz_0[i] * fe_0 - 2.0 * ta_zzz_yz_1[i] * fe_0 + ta_zzz_yyz_0[i] * pa_y[i] - ta_zzz_yyz_1[i] * pc_y[i];

        ta_yzzz_yzz_0[i] = ta_zzz_zz_0[i] * fe_0 - ta_zzz_zz_1[i] * fe_0 + ta_zzz_yzz_0[i] * pa_y[i] - ta_zzz_yzz_1[i] * pc_y[i];

        ta_yzzz_zzz_0[i] = ta_zzz_zzz_0[i] * pa_y[i] - ta_zzz_zzz_1[i] * pc_y[i];
    }

    // Set up 140-150 components of targeted buffer : GF

    auto ta_zzzz_xxx_0 = pbuffer.data(idx_npot_0_gf + 140);

    auto ta_zzzz_xxy_0 = pbuffer.data(idx_npot_0_gf + 141);

    auto ta_zzzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 142);

    auto ta_zzzz_xyy_0 = pbuffer.data(idx_npot_0_gf + 143);

    auto ta_zzzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 144);

    auto ta_zzzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 145);

    auto ta_zzzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 146);

    auto ta_zzzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 147);

    auto ta_zzzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 148);

    auto ta_zzzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 149);

#pragma omp simd aligned(pa_z,              \
                             pc_z,          \
                             ta_zz_xxx_0,   \
                             ta_zz_xxx_1,   \
                             ta_zz_xxy_0,   \
                             ta_zz_xxy_1,   \
                             ta_zz_xxz_0,   \
                             ta_zz_xxz_1,   \
                             ta_zz_xyy_0,   \
                             ta_zz_xyy_1,   \
                             ta_zz_xyz_0,   \
                             ta_zz_xyz_1,   \
                             ta_zz_xzz_0,   \
                             ta_zz_xzz_1,   \
                             ta_zz_yyy_0,   \
                             ta_zz_yyy_1,   \
                             ta_zz_yyz_0,   \
                             ta_zz_yyz_1,   \
                             ta_zz_yzz_0,   \
                             ta_zz_yzz_1,   \
                             ta_zz_zzz_0,   \
                             ta_zz_zzz_1,   \
                             ta_zzz_xx_0,   \
                             ta_zzz_xx_1,   \
                             ta_zzz_xxx_0,  \
                             ta_zzz_xxx_1,  \
                             ta_zzz_xxy_0,  \
                             ta_zzz_xxy_1,  \
                             ta_zzz_xxz_0,  \
                             ta_zzz_xxz_1,  \
                             ta_zzz_xy_0,   \
                             ta_zzz_xy_1,   \
                             ta_zzz_xyy_0,  \
                             ta_zzz_xyy_1,  \
                             ta_zzz_xyz_0,  \
                             ta_zzz_xyz_1,  \
                             ta_zzz_xz_0,   \
                             ta_zzz_xz_1,   \
                             ta_zzz_xzz_0,  \
                             ta_zzz_xzz_1,  \
                             ta_zzz_yy_0,   \
                             ta_zzz_yy_1,   \
                             ta_zzz_yyy_0,  \
                             ta_zzz_yyy_1,  \
                             ta_zzz_yyz_0,  \
                             ta_zzz_yyz_1,  \
                             ta_zzz_yz_0,   \
                             ta_zzz_yz_1,   \
                             ta_zzz_yzz_0,  \
                             ta_zzz_yzz_1,  \
                             ta_zzz_zz_0,   \
                             ta_zzz_zz_1,   \
                             ta_zzz_zzz_0,  \
                             ta_zzz_zzz_1,  \
                             ta_zzzz_xxx_0, \
                             ta_zzzz_xxy_0, \
                             ta_zzzz_xxz_0, \
                             ta_zzzz_xyy_0, \
                             ta_zzzz_xyz_0, \
                             ta_zzzz_xzz_0, \
                             ta_zzzz_yyy_0, \
                             ta_zzzz_yyz_0, \
                             ta_zzzz_yzz_0, \
                             ta_zzzz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzz_xxx_0[i] = 3.0 * ta_zz_xxx_0[i] * fe_0 - 3.0 * ta_zz_xxx_1[i] * fe_0 + ta_zzz_xxx_0[i] * pa_z[i] - ta_zzz_xxx_1[i] * pc_z[i];

        ta_zzzz_xxy_0[i] = 3.0 * ta_zz_xxy_0[i] * fe_0 - 3.0 * ta_zz_xxy_1[i] * fe_0 + ta_zzz_xxy_0[i] * pa_z[i] - ta_zzz_xxy_1[i] * pc_z[i];

        ta_zzzz_xxz_0[i] = 3.0 * ta_zz_xxz_0[i] * fe_0 - 3.0 * ta_zz_xxz_1[i] * fe_0 + ta_zzz_xx_0[i] * fe_0 - ta_zzz_xx_1[i] * fe_0 +
                           ta_zzz_xxz_0[i] * pa_z[i] - ta_zzz_xxz_1[i] * pc_z[i];

        ta_zzzz_xyy_0[i] = 3.0 * ta_zz_xyy_0[i] * fe_0 - 3.0 * ta_zz_xyy_1[i] * fe_0 + ta_zzz_xyy_0[i] * pa_z[i] - ta_zzz_xyy_1[i] * pc_z[i];

        ta_zzzz_xyz_0[i] = 3.0 * ta_zz_xyz_0[i] * fe_0 - 3.0 * ta_zz_xyz_1[i] * fe_0 + ta_zzz_xy_0[i] * fe_0 - ta_zzz_xy_1[i] * fe_0 +
                           ta_zzz_xyz_0[i] * pa_z[i] - ta_zzz_xyz_1[i] * pc_z[i];

        ta_zzzz_xzz_0[i] = 3.0 * ta_zz_xzz_0[i] * fe_0 - 3.0 * ta_zz_xzz_1[i] * fe_0 + 2.0 * ta_zzz_xz_0[i] * fe_0 - 2.0 * ta_zzz_xz_1[i] * fe_0 +
                           ta_zzz_xzz_0[i] * pa_z[i] - ta_zzz_xzz_1[i] * pc_z[i];

        ta_zzzz_yyy_0[i] = 3.0 * ta_zz_yyy_0[i] * fe_0 - 3.0 * ta_zz_yyy_1[i] * fe_0 + ta_zzz_yyy_0[i] * pa_z[i] - ta_zzz_yyy_1[i] * pc_z[i];

        ta_zzzz_yyz_0[i] = 3.0 * ta_zz_yyz_0[i] * fe_0 - 3.0 * ta_zz_yyz_1[i] * fe_0 + ta_zzz_yy_0[i] * fe_0 - ta_zzz_yy_1[i] * fe_0 +
                           ta_zzz_yyz_0[i] * pa_z[i] - ta_zzz_yyz_1[i] * pc_z[i];

        ta_zzzz_yzz_0[i] = 3.0 * ta_zz_yzz_0[i] * fe_0 - 3.0 * ta_zz_yzz_1[i] * fe_0 + 2.0 * ta_zzz_yz_0[i] * fe_0 - 2.0 * ta_zzz_yz_1[i] * fe_0 +
                           ta_zzz_yzz_0[i] * pa_z[i] - ta_zzz_yzz_1[i] * pc_z[i];

        ta_zzzz_zzz_0[i] = 3.0 * ta_zz_zzz_0[i] * fe_0 - 3.0 * ta_zz_zzz_1[i] * fe_0 + 3.0 * ta_zzz_zz_0[i] * fe_0 - 3.0 * ta_zzz_zz_1[i] * fe_0 +
                           ta_zzz_zzz_0[i] * pa_z[i] - ta_zzz_zzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
