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

#include "NuclearPotentialPrimRecHF.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_hf(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_hf,
                               const size_t              idx_npot_0_ff,
                               const size_t              idx_npot_1_ff,
                               const size_t              idx_npot_0_gd,
                               const size_t              idx_npot_1_gd,
                               const size_t              idx_npot_0_gf,
                               const size_t              idx_npot_1_gf,
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

    auto ta_xxy_xxz_0 = pbuffer.data(idx_npot_0_ff + 12);

    auto ta_xxy_xzz_0 = pbuffer.data(idx_npot_0_ff + 15);

    auto ta_xxy_yyy_0 = pbuffer.data(idx_npot_0_ff + 16);

    auto ta_xxy_yyz_0 = pbuffer.data(idx_npot_0_ff + 17);

    auto ta_xxy_yzz_0 = pbuffer.data(idx_npot_0_ff + 18);

    auto ta_xxz_xxx_0 = pbuffer.data(idx_npot_0_ff + 20);

    auto ta_xxz_xxy_0 = pbuffer.data(idx_npot_0_ff + 21);

    auto ta_xxz_xxz_0 = pbuffer.data(idx_npot_0_ff + 22);

    auto ta_xxz_xyy_0 = pbuffer.data(idx_npot_0_ff + 23);

    auto ta_xxz_xzz_0 = pbuffer.data(idx_npot_0_ff + 25);

    auto ta_xxz_yyz_0 = pbuffer.data(idx_npot_0_ff + 27);

    auto ta_xxz_yzz_0 = pbuffer.data(idx_npot_0_ff + 28);

    auto ta_xxz_zzz_0 = pbuffer.data(idx_npot_0_ff + 29);

    auto ta_xyy_xxy_0 = pbuffer.data(idx_npot_0_ff + 31);

    auto ta_xyy_xyy_0 = pbuffer.data(idx_npot_0_ff + 33);

    auto ta_xyy_xyz_0 = pbuffer.data(idx_npot_0_ff + 34);

    auto ta_xyy_yyy_0 = pbuffer.data(idx_npot_0_ff + 36);

    auto ta_xyy_yyz_0 = pbuffer.data(idx_npot_0_ff + 37);

    auto ta_xyy_yzz_0 = pbuffer.data(idx_npot_0_ff + 38);

    auto ta_xyy_zzz_0 = pbuffer.data(idx_npot_0_ff + 39);

    auto ta_xyz_yyz_0 = pbuffer.data(idx_npot_0_ff + 47);

    auto ta_xyz_yzz_0 = pbuffer.data(idx_npot_0_ff + 48);

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

    auto ta_yyz_xzz_0 = pbuffer.data(idx_npot_0_ff + 75);

    auto ta_yyz_yyy_0 = pbuffer.data(idx_npot_0_ff + 76);

    auto ta_yyz_yyz_0 = pbuffer.data(idx_npot_0_ff + 77);

    auto ta_yyz_yzz_0 = pbuffer.data(idx_npot_0_ff + 78);

    auto ta_yyz_zzz_0 = pbuffer.data(idx_npot_0_ff + 79);

    auto ta_yzz_xxx_0 = pbuffer.data(idx_npot_0_ff + 80);

    auto ta_yzz_xxz_0 = pbuffer.data(idx_npot_0_ff + 82);

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

    auto ta_xxy_xxz_1 = pbuffer.data(idx_npot_1_ff + 12);

    auto ta_xxy_xzz_1 = pbuffer.data(idx_npot_1_ff + 15);

    auto ta_xxy_yyy_1 = pbuffer.data(idx_npot_1_ff + 16);

    auto ta_xxy_yyz_1 = pbuffer.data(idx_npot_1_ff + 17);

    auto ta_xxy_yzz_1 = pbuffer.data(idx_npot_1_ff + 18);

    auto ta_xxz_xxx_1 = pbuffer.data(idx_npot_1_ff + 20);

    auto ta_xxz_xxy_1 = pbuffer.data(idx_npot_1_ff + 21);

    auto ta_xxz_xxz_1 = pbuffer.data(idx_npot_1_ff + 22);

    auto ta_xxz_xyy_1 = pbuffer.data(idx_npot_1_ff + 23);

    auto ta_xxz_xzz_1 = pbuffer.data(idx_npot_1_ff + 25);

    auto ta_xxz_yyz_1 = pbuffer.data(idx_npot_1_ff + 27);

    auto ta_xxz_yzz_1 = pbuffer.data(idx_npot_1_ff + 28);

    auto ta_xxz_zzz_1 = pbuffer.data(idx_npot_1_ff + 29);

    auto ta_xyy_xxy_1 = pbuffer.data(idx_npot_1_ff + 31);

    auto ta_xyy_xyy_1 = pbuffer.data(idx_npot_1_ff + 33);

    auto ta_xyy_xyz_1 = pbuffer.data(idx_npot_1_ff + 34);

    auto ta_xyy_yyy_1 = pbuffer.data(idx_npot_1_ff + 36);

    auto ta_xyy_yyz_1 = pbuffer.data(idx_npot_1_ff + 37);

    auto ta_xyy_yzz_1 = pbuffer.data(idx_npot_1_ff + 38);

    auto ta_xyy_zzz_1 = pbuffer.data(idx_npot_1_ff + 39);

    auto ta_xyz_yyz_1 = pbuffer.data(idx_npot_1_ff + 47);

    auto ta_xyz_yzz_1 = pbuffer.data(idx_npot_1_ff + 48);

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

    auto ta_yyz_xzz_1 = pbuffer.data(idx_npot_1_ff + 75);

    auto ta_yyz_yyy_1 = pbuffer.data(idx_npot_1_ff + 76);

    auto ta_yyz_yyz_1 = pbuffer.data(idx_npot_1_ff + 77);

    auto ta_yyz_yzz_1 = pbuffer.data(idx_npot_1_ff + 78);

    auto ta_yyz_zzz_1 = pbuffer.data(idx_npot_1_ff + 79);

    auto ta_yzz_xxx_1 = pbuffer.data(idx_npot_1_ff + 80);

    auto ta_yzz_xxz_1 = pbuffer.data(idx_npot_1_ff + 82);

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

    // Set up components of auxiliary buffer : GD

    auto ta_xxxx_xx_0 = pbuffer.data(idx_npot_0_gd);

    auto ta_xxxx_xy_0 = pbuffer.data(idx_npot_0_gd + 1);

    auto ta_xxxx_xz_0 = pbuffer.data(idx_npot_0_gd + 2);

    auto ta_xxxx_yy_0 = pbuffer.data(idx_npot_0_gd + 3);

    auto ta_xxxx_yz_0 = pbuffer.data(idx_npot_0_gd + 4);

    auto ta_xxxx_zz_0 = pbuffer.data(idx_npot_0_gd + 5);

    auto ta_xxxz_xz_0 = pbuffer.data(idx_npot_0_gd + 14);

    auto ta_xxyy_xy_0 = pbuffer.data(idx_npot_0_gd + 19);

    auto ta_xxyy_yy_0 = pbuffer.data(idx_npot_0_gd + 21);

    auto ta_xxyy_yz_0 = pbuffer.data(idx_npot_0_gd + 22);

    auto ta_xxzz_xx_0 = pbuffer.data(idx_npot_0_gd + 30);

    auto ta_xxzz_xy_0 = pbuffer.data(idx_npot_0_gd + 31);

    auto ta_xxzz_xz_0 = pbuffer.data(idx_npot_0_gd + 32);

    auto ta_xxzz_yz_0 = pbuffer.data(idx_npot_0_gd + 34);

    auto ta_xxzz_zz_0 = pbuffer.data(idx_npot_0_gd + 35);

    auto ta_xyyy_xy_0 = pbuffer.data(idx_npot_0_gd + 37);

    auto ta_xyyy_yy_0 = pbuffer.data(idx_npot_0_gd + 39);

    auto ta_xyyy_yz_0 = pbuffer.data(idx_npot_0_gd + 40);

    auto ta_xzzz_xz_0 = pbuffer.data(idx_npot_0_gd + 56);

    auto ta_xzzz_yz_0 = pbuffer.data(idx_npot_0_gd + 58);

    auto ta_xzzz_zz_0 = pbuffer.data(idx_npot_0_gd + 59);

    auto ta_yyyy_xx_0 = pbuffer.data(idx_npot_0_gd + 60);

    auto ta_yyyy_xy_0 = pbuffer.data(idx_npot_0_gd + 61);

    auto ta_yyyy_xz_0 = pbuffer.data(idx_npot_0_gd + 62);

    auto ta_yyyy_yy_0 = pbuffer.data(idx_npot_0_gd + 63);

    auto ta_yyyy_yz_0 = pbuffer.data(idx_npot_0_gd + 64);

    auto ta_yyyy_zz_0 = pbuffer.data(idx_npot_0_gd + 65);

    auto ta_yyyz_xz_0 = pbuffer.data(idx_npot_0_gd + 68);

    auto ta_yyyz_yz_0 = pbuffer.data(idx_npot_0_gd + 70);

    auto ta_yyyz_zz_0 = pbuffer.data(idx_npot_0_gd + 71);

    auto ta_yyzz_xx_0 = pbuffer.data(idx_npot_0_gd + 72);

    auto ta_yyzz_xy_0 = pbuffer.data(idx_npot_0_gd + 73);

    auto ta_yyzz_xz_0 = pbuffer.data(idx_npot_0_gd + 74);

    auto ta_yyzz_yy_0 = pbuffer.data(idx_npot_0_gd + 75);

    auto ta_yyzz_yz_0 = pbuffer.data(idx_npot_0_gd + 76);

    auto ta_yyzz_zz_0 = pbuffer.data(idx_npot_0_gd + 77);

    auto ta_yzzz_xy_0 = pbuffer.data(idx_npot_0_gd + 79);

    auto ta_yzzz_xz_0 = pbuffer.data(idx_npot_0_gd + 80);

    auto ta_yzzz_yy_0 = pbuffer.data(idx_npot_0_gd + 81);

    auto ta_yzzz_yz_0 = pbuffer.data(idx_npot_0_gd + 82);

    auto ta_yzzz_zz_0 = pbuffer.data(idx_npot_0_gd + 83);

    auto ta_zzzz_xx_0 = pbuffer.data(idx_npot_0_gd + 84);

    auto ta_zzzz_xy_0 = pbuffer.data(idx_npot_0_gd + 85);

    auto ta_zzzz_xz_0 = pbuffer.data(idx_npot_0_gd + 86);

    auto ta_zzzz_yy_0 = pbuffer.data(idx_npot_0_gd + 87);

    auto ta_zzzz_yz_0 = pbuffer.data(idx_npot_0_gd + 88);

    auto ta_zzzz_zz_0 = pbuffer.data(idx_npot_0_gd + 89);

    // Set up components of auxiliary buffer : GD

    auto ta_xxxx_xx_1 = pbuffer.data(idx_npot_1_gd);

    auto ta_xxxx_xy_1 = pbuffer.data(idx_npot_1_gd + 1);

    auto ta_xxxx_xz_1 = pbuffer.data(idx_npot_1_gd + 2);

    auto ta_xxxx_yy_1 = pbuffer.data(idx_npot_1_gd + 3);

    auto ta_xxxx_yz_1 = pbuffer.data(idx_npot_1_gd + 4);

    auto ta_xxxx_zz_1 = pbuffer.data(idx_npot_1_gd + 5);

    auto ta_xxxz_xz_1 = pbuffer.data(idx_npot_1_gd + 14);

    auto ta_xxyy_xy_1 = pbuffer.data(idx_npot_1_gd + 19);

    auto ta_xxyy_yy_1 = pbuffer.data(idx_npot_1_gd + 21);

    auto ta_xxyy_yz_1 = pbuffer.data(idx_npot_1_gd + 22);

    auto ta_xxzz_xx_1 = pbuffer.data(idx_npot_1_gd + 30);

    auto ta_xxzz_xy_1 = pbuffer.data(idx_npot_1_gd + 31);

    auto ta_xxzz_xz_1 = pbuffer.data(idx_npot_1_gd + 32);

    auto ta_xxzz_yz_1 = pbuffer.data(idx_npot_1_gd + 34);

    auto ta_xxzz_zz_1 = pbuffer.data(idx_npot_1_gd + 35);

    auto ta_xyyy_xy_1 = pbuffer.data(idx_npot_1_gd + 37);

    auto ta_xyyy_yy_1 = pbuffer.data(idx_npot_1_gd + 39);

    auto ta_xyyy_yz_1 = pbuffer.data(idx_npot_1_gd + 40);

    auto ta_xzzz_xz_1 = pbuffer.data(idx_npot_1_gd + 56);

    auto ta_xzzz_yz_1 = pbuffer.data(idx_npot_1_gd + 58);

    auto ta_xzzz_zz_1 = pbuffer.data(idx_npot_1_gd + 59);

    auto ta_yyyy_xx_1 = pbuffer.data(idx_npot_1_gd + 60);

    auto ta_yyyy_xy_1 = pbuffer.data(idx_npot_1_gd + 61);

    auto ta_yyyy_xz_1 = pbuffer.data(idx_npot_1_gd + 62);

    auto ta_yyyy_yy_1 = pbuffer.data(idx_npot_1_gd + 63);

    auto ta_yyyy_yz_1 = pbuffer.data(idx_npot_1_gd + 64);

    auto ta_yyyy_zz_1 = pbuffer.data(idx_npot_1_gd + 65);

    auto ta_yyyz_xz_1 = pbuffer.data(idx_npot_1_gd + 68);

    auto ta_yyyz_yz_1 = pbuffer.data(idx_npot_1_gd + 70);

    auto ta_yyyz_zz_1 = pbuffer.data(idx_npot_1_gd + 71);

    auto ta_yyzz_xx_1 = pbuffer.data(idx_npot_1_gd + 72);

    auto ta_yyzz_xy_1 = pbuffer.data(idx_npot_1_gd + 73);

    auto ta_yyzz_xz_1 = pbuffer.data(idx_npot_1_gd + 74);

    auto ta_yyzz_yy_1 = pbuffer.data(idx_npot_1_gd + 75);

    auto ta_yyzz_yz_1 = pbuffer.data(idx_npot_1_gd + 76);

    auto ta_yyzz_zz_1 = pbuffer.data(idx_npot_1_gd + 77);

    auto ta_yzzz_xy_1 = pbuffer.data(idx_npot_1_gd + 79);

    auto ta_yzzz_xz_1 = pbuffer.data(idx_npot_1_gd + 80);

    auto ta_yzzz_yy_1 = pbuffer.data(idx_npot_1_gd + 81);

    auto ta_yzzz_yz_1 = pbuffer.data(idx_npot_1_gd + 82);

    auto ta_yzzz_zz_1 = pbuffer.data(idx_npot_1_gd + 83);

    auto ta_zzzz_xx_1 = pbuffer.data(idx_npot_1_gd + 84);

    auto ta_zzzz_xy_1 = pbuffer.data(idx_npot_1_gd + 85);

    auto ta_zzzz_xz_1 = pbuffer.data(idx_npot_1_gd + 86);

    auto ta_zzzz_yy_1 = pbuffer.data(idx_npot_1_gd + 87);

    auto ta_zzzz_yz_1 = pbuffer.data(idx_npot_1_gd + 88);

    auto ta_zzzz_zz_1 = pbuffer.data(idx_npot_1_gd + 89);

    // Set up components of auxiliary buffer : GF

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

    auto ta_xxxy_xxx_0 = pbuffer.data(idx_npot_0_gf + 10);

    auto ta_xxxy_xxy_0 = pbuffer.data(idx_npot_0_gf + 11);

    auto ta_xxxy_xxz_0 = pbuffer.data(idx_npot_0_gf + 12);

    auto ta_xxxy_xyy_0 = pbuffer.data(idx_npot_0_gf + 13);

    auto ta_xxxy_xzz_0 = pbuffer.data(idx_npot_0_gf + 15);

    auto ta_xxxy_yyy_0 = pbuffer.data(idx_npot_0_gf + 16);

    auto ta_xxxy_yyz_0 = pbuffer.data(idx_npot_0_gf + 17);

    auto ta_xxxy_yzz_0 = pbuffer.data(idx_npot_0_gf + 18);

    auto ta_xxxz_xxx_0 = pbuffer.data(idx_npot_0_gf + 20);

    auto ta_xxxz_xxy_0 = pbuffer.data(idx_npot_0_gf + 21);

    auto ta_xxxz_xxz_0 = pbuffer.data(idx_npot_0_gf + 22);

    auto ta_xxxz_xyy_0 = pbuffer.data(idx_npot_0_gf + 23);

    auto ta_xxxz_xyz_0 = pbuffer.data(idx_npot_0_gf + 24);

    auto ta_xxxz_xzz_0 = pbuffer.data(idx_npot_0_gf + 25);

    auto ta_xxxz_yyz_0 = pbuffer.data(idx_npot_0_gf + 27);

    auto ta_xxxz_yzz_0 = pbuffer.data(idx_npot_0_gf + 28);

    auto ta_xxxz_zzz_0 = pbuffer.data(idx_npot_0_gf + 29);

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

    auto ta_xxyz_xxz_0 = pbuffer.data(idx_npot_0_gf + 42);

    auto ta_xxyz_xzz_0 = pbuffer.data(idx_npot_0_gf + 45);

    auto ta_xxyz_yyz_0 = pbuffer.data(idx_npot_0_gf + 47);

    auto ta_xxyz_yzz_0 = pbuffer.data(idx_npot_0_gf + 48);

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

    auto ta_xyyy_xxx_0 = pbuffer.data(idx_npot_0_gf + 60);

    auto ta_xyyy_xxy_0 = pbuffer.data(idx_npot_0_gf + 61);

    auto ta_xyyy_xyy_0 = pbuffer.data(idx_npot_0_gf + 63);

    auto ta_xyyy_xyz_0 = pbuffer.data(idx_npot_0_gf + 64);

    auto ta_xyyy_yyy_0 = pbuffer.data(idx_npot_0_gf + 66);

    auto ta_xyyy_yyz_0 = pbuffer.data(idx_npot_0_gf + 67);

    auto ta_xyyy_yzz_0 = pbuffer.data(idx_npot_0_gf + 68);

    auto ta_xyyy_zzz_0 = pbuffer.data(idx_npot_0_gf + 69);

    auto ta_xyyz_yyz_0 = pbuffer.data(idx_npot_0_gf + 77);

    auto ta_xyyz_yzz_0 = pbuffer.data(idx_npot_0_gf + 78);

    auto ta_xyyz_zzz_0 = pbuffer.data(idx_npot_0_gf + 79);

    auto ta_xyzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 86);

    auto ta_xyzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 87);

    auto ta_xyzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 88);

    auto ta_xzzz_xxx_0 = pbuffer.data(idx_npot_0_gf + 90);

    auto ta_xzzz_xxz_0 = pbuffer.data(idx_npot_0_gf + 92);

    auto ta_xzzz_xyz_0 = pbuffer.data(idx_npot_0_gf + 94);

    auto ta_xzzz_xzz_0 = pbuffer.data(idx_npot_0_gf + 95);

    auto ta_xzzz_yyy_0 = pbuffer.data(idx_npot_0_gf + 96);

    auto ta_xzzz_yyz_0 = pbuffer.data(idx_npot_0_gf + 97);

    auto ta_xzzz_yzz_0 = pbuffer.data(idx_npot_0_gf + 98);

    auto ta_xzzz_zzz_0 = pbuffer.data(idx_npot_0_gf + 99);

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

    auto ta_yyyz_xxy_0 = pbuffer.data(idx_npot_0_gf + 111);

    auto ta_yyyz_xxz_0 = pbuffer.data(idx_npot_0_gf + 112);

    auto ta_yyyz_xyy_0 = pbuffer.data(idx_npot_0_gf + 113);

    auto ta_yyyz_xyz_0 = pbuffer.data(idx_npot_0_gf + 114);

    auto ta_yyyz_xzz_0 = pbuffer.data(idx_npot_0_gf + 115);

    auto ta_yyyz_yyy_0 = pbuffer.data(idx_npot_0_gf + 116);

    auto ta_yyyz_yyz_0 = pbuffer.data(idx_npot_0_gf + 117);

    auto ta_yyyz_yzz_0 = pbuffer.data(idx_npot_0_gf + 118);

    auto ta_yyyz_zzz_0 = pbuffer.data(idx_npot_0_gf + 119);

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

    // Set up components of auxiliary buffer : GF

    auto ta_xxxx_xxx_1 = pbuffer.data(idx_npot_1_gf);

    auto ta_xxxx_xxy_1 = pbuffer.data(idx_npot_1_gf + 1);

    auto ta_xxxx_xxz_1 = pbuffer.data(idx_npot_1_gf + 2);

    auto ta_xxxx_xyy_1 = pbuffer.data(idx_npot_1_gf + 3);

    auto ta_xxxx_xyz_1 = pbuffer.data(idx_npot_1_gf + 4);

    auto ta_xxxx_xzz_1 = pbuffer.data(idx_npot_1_gf + 5);

    auto ta_xxxx_yyy_1 = pbuffer.data(idx_npot_1_gf + 6);

    auto ta_xxxx_yyz_1 = pbuffer.data(idx_npot_1_gf + 7);

    auto ta_xxxx_yzz_1 = pbuffer.data(idx_npot_1_gf + 8);

    auto ta_xxxx_zzz_1 = pbuffer.data(idx_npot_1_gf + 9);

    auto ta_xxxy_xxx_1 = pbuffer.data(idx_npot_1_gf + 10);

    auto ta_xxxy_xxy_1 = pbuffer.data(idx_npot_1_gf + 11);

    auto ta_xxxy_xxz_1 = pbuffer.data(idx_npot_1_gf + 12);

    auto ta_xxxy_xyy_1 = pbuffer.data(idx_npot_1_gf + 13);

    auto ta_xxxy_xzz_1 = pbuffer.data(idx_npot_1_gf + 15);

    auto ta_xxxy_yyy_1 = pbuffer.data(idx_npot_1_gf + 16);

    auto ta_xxxy_yyz_1 = pbuffer.data(idx_npot_1_gf + 17);

    auto ta_xxxy_yzz_1 = pbuffer.data(idx_npot_1_gf + 18);

    auto ta_xxxz_xxx_1 = pbuffer.data(idx_npot_1_gf + 20);

    auto ta_xxxz_xxy_1 = pbuffer.data(idx_npot_1_gf + 21);

    auto ta_xxxz_xxz_1 = pbuffer.data(idx_npot_1_gf + 22);

    auto ta_xxxz_xyy_1 = pbuffer.data(idx_npot_1_gf + 23);

    auto ta_xxxz_xyz_1 = pbuffer.data(idx_npot_1_gf + 24);

    auto ta_xxxz_xzz_1 = pbuffer.data(idx_npot_1_gf + 25);

    auto ta_xxxz_yyz_1 = pbuffer.data(idx_npot_1_gf + 27);

    auto ta_xxxz_yzz_1 = pbuffer.data(idx_npot_1_gf + 28);

    auto ta_xxxz_zzz_1 = pbuffer.data(idx_npot_1_gf + 29);

    auto ta_xxyy_xxx_1 = pbuffer.data(idx_npot_1_gf + 30);

    auto ta_xxyy_xxy_1 = pbuffer.data(idx_npot_1_gf + 31);

    auto ta_xxyy_xxz_1 = pbuffer.data(idx_npot_1_gf + 32);

    auto ta_xxyy_xyy_1 = pbuffer.data(idx_npot_1_gf + 33);

    auto ta_xxyy_xyz_1 = pbuffer.data(idx_npot_1_gf + 34);

    auto ta_xxyy_xzz_1 = pbuffer.data(idx_npot_1_gf + 35);

    auto ta_xxyy_yyy_1 = pbuffer.data(idx_npot_1_gf + 36);

    auto ta_xxyy_yyz_1 = pbuffer.data(idx_npot_1_gf + 37);

    auto ta_xxyy_yzz_1 = pbuffer.data(idx_npot_1_gf + 38);

    auto ta_xxyy_zzz_1 = pbuffer.data(idx_npot_1_gf + 39);

    auto ta_xxyz_xxz_1 = pbuffer.data(idx_npot_1_gf + 42);

    auto ta_xxyz_xzz_1 = pbuffer.data(idx_npot_1_gf + 45);

    auto ta_xxyz_yyz_1 = pbuffer.data(idx_npot_1_gf + 47);

    auto ta_xxyz_yzz_1 = pbuffer.data(idx_npot_1_gf + 48);

    auto ta_xxzz_xxx_1 = pbuffer.data(idx_npot_1_gf + 50);

    auto ta_xxzz_xxy_1 = pbuffer.data(idx_npot_1_gf + 51);

    auto ta_xxzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 52);

    auto ta_xxzz_xyy_1 = pbuffer.data(idx_npot_1_gf + 53);

    auto ta_xxzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 54);

    auto ta_xxzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 55);

    auto ta_xxzz_yyy_1 = pbuffer.data(idx_npot_1_gf + 56);

    auto ta_xxzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 57);

    auto ta_xxzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 58);

    auto ta_xxzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 59);

    auto ta_xyyy_xxx_1 = pbuffer.data(idx_npot_1_gf + 60);

    auto ta_xyyy_xxy_1 = pbuffer.data(idx_npot_1_gf + 61);

    auto ta_xyyy_xyy_1 = pbuffer.data(idx_npot_1_gf + 63);

    auto ta_xyyy_xyz_1 = pbuffer.data(idx_npot_1_gf + 64);

    auto ta_xyyy_yyy_1 = pbuffer.data(idx_npot_1_gf + 66);

    auto ta_xyyy_yyz_1 = pbuffer.data(idx_npot_1_gf + 67);

    auto ta_xyyy_yzz_1 = pbuffer.data(idx_npot_1_gf + 68);

    auto ta_xyyy_zzz_1 = pbuffer.data(idx_npot_1_gf + 69);

    auto ta_xyyz_yyz_1 = pbuffer.data(idx_npot_1_gf + 77);

    auto ta_xyyz_yzz_1 = pbuffer.data(idx_npot_1_gf + 78);

    auto ta_xyyz_zzz_1 = pbuffer.data(idx_npot_1_gf + 79);

    auto ta_xyzz_yyy_1 = pbuffer.data(idx_npot_1_gf + 86);

    auto ta_xyzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 87);

    auto ta_xyzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 88);

    auto ta_xzzz_xxx_1 = pbuffer.data(idx_npot_1_gf + 90);

    auto ta_xzzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 92);

    auto ta_xzzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 94);

    auto ta_xzzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 95);

    auto ta_xzzz_yyy_1 = pbuffer.data(idx_npot_1_gf + 96);

    auto ta_xzzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 97);

    auto ta_xzzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 98);

    auto ta_xzzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 99);

    auto ta_yyyy_xxx_1 = pbuffer.data(idx_npot_1_gf + 100);

    auto ta_yyyy_xxy_1 = pbuffer.data(idx_npot_1_gf + 101);

    auto ta_yyyy_xxz_1 = pbuffer.data(idx_npot_1_gf + 102);

    auto ta_yyyy_xyy_1 = pbuffer.data(idx_npot_1_gf + 103);

    auto ta_yyyy_xyz_1 = pbuffer.data(idx_npot_1_gf + 104);

    auto ta_yyyy_xzz_1 = pbuffer.data(idx_npot_1_gf + 105);

    auto ta_yyyy_yyy_1 = pbuffer.data(idx_npot_1_gf + 106);

    auto ta_yyyy_yyz_1 = pbuffer.data(idx_npot_1_gf + 107);

    auto ta_yyyy_yzz_1 = pbuffer.data(idx_npot_1_gf + 108);

    auto ta_yyyy_zzz_1 = pbuffer.data(idx_npot_1_gf + 109);

    auto ta_yyyz_xxy_1 = pbuffer.data(idx_npot_1_gf + 111);

    auto ta_yyyz_xxz_1 = pbuffer.data(idx_npot_1_gf + 112);

    auto ta_yyyz_xyy_1 = pbuffer.data(idx_npot_1_gf + 113);

    auto ta_yyyz_xyz_1 = pbuffer.data(idx_npot_1_gf + 114);

    auto ta_yyyz_xzz_1 = pbuffer.data(idx_npot_1_gf + 115);

    auto ta_yyyz_yyy_1 = pbuffer.data(idx_npot_1_gf + 116);

    auto ta_yyyz_yyz_1 = pbuffer.data(idx_npot_1_gf + 117);

    auto ta_yyyz_yzz_1 = pbuffer.data(idx_npot_1_gf + 118);

    auto ta_yyyz_zzz_1 = pbuffer.data(idx_npot_1_gf + 119);

    auto ta_yyzz_xxx_1 = pbuffer.data(idx_npot_1_gf + 120);

    auto ta_yyzz_xxy_1 = pbuffer.data(idx_npot_1_gf + 121);

    auto ta_yyzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 122);

    auto ta_yyzz_xyy_1 = pbuffer.data(idx_npot_1_gf + 123);

    auto ta_yyzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 124);

    auto ta_yyzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 125);

    auto ta_yyzz_yyy_1 = pbuffer.data(idx_npot_1_gf + 126);

    auto ta_yyzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 127);

    auto ta_yyzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 128);

    auto ta_yyzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 129);

    auto ta_yzzz_xxx_1 = pbuffer.data(idx_npot_1_gf + 130);

    auto ta_yzzz_xxy_1 = pbuffer.data(idx_npot_1_gf + 131);

    auto ta_yzzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 132);

    auto ta_yzzz_xyy_1 = pbuffer.data(idx_npot_1_gf + 133);

    auto ta_yzzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 134);

    auto ta_yzzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 135);

    auto ta_yzzz_yyy_1 = pbuffer.data(idx_npot_1_gf + 136);

    auto ta_yzzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 137);

    auto ta_yzzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 138);

    auto ta_yzzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 139);

    auto ta_zzzz_xxx_1 = pbuffer.data(idx_npot_1_gf + 140);

    auto ta_zzzz_xxy_1 = pbuffer.data(idx_npot_1_gf + 141);

    auto ta_zzzz_xxz_1 = pbuffer.data(idx_npot_1_gf + 142);

    auto ta_zzzz_xyy_1 = pbuffer.data(idx_npot_1_gf + 143);

    auto ta_zzzz_xyz_1 = pbuffer.data(idx_npot_1_gf + 144);

    auto ta_zzzz_xzz_1 = pbuffer.data(idx_npot_1_gf + 145);

    auto ta_zzzz_yyy_1 = pbuffer.data(idx_npot_1_gf + 146);

    auto ta_zzzz_yyz_1 = pbuffer.data(idx_npot_1_gf + 147);

    auto ta_zzzz_yzz_1 = pbuffer.data(idx_npot_1_gf + 148);

    auto ta_zzzz_zzz_1 = pbuffer.data(idx_npot_1_gf + 149);

    // Set up 0-10 components of targeted buffer : HF

    auto ta_xxxxx_xxx_0 = pbuffer.data(idx_npot_0_hf);

    auto ta_xxxxx_xxy_0 = pbuffer.data(idx_npot_0_hf + 1);

    auto ta_xxxxx_xxz_0 = pbuffer.data(idx_npot_0_hf + 2);

    auto ta_xxxxx_xyy_0 = pbuffer.data(idx_npot_0_hf + 3);

    auto ta_xxxxx_xyz_0 = pbuffer.data(idx_npot_0_hf + 4);

    auto ta_xxxxx_xzz_0 = pbuffer.data(idx_npot_0_hf + 5);

    auto ta_xxxxx_yyy_0 = pbuffer.data(idx_npot_0_hf + 6);

    auto ta_xxxxx_yyz_0 = pbuffer.data(idx_npot_0_hf + 7);

    auto ta_xxxxx_yzz_0 = pbuffer.data(idx_npot_0_hf + 8);

    auto ta_xxxxx_zzz_0 = pbuffer.data(idx_npot_0_hf + 9);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta_xxx_xxx_0,   \
                             ta_xxx_xxx_1,   \
                             ta_xxx_xxy_0,   \
                             ta_xxx_xxy_1,   \
                             ta_xxx_xxz_0,   \
                             ta_xxx_xxz_1,   \
                             ta_xxx_xyy_0,   \
                             ta_xxx_xyy_1,   \
                             ta_xxx_xyz_0,   \
                             ta_xxx_xyz_1,   \
                             ta_xxx_xzz_0,   \
                             ta_xxx_xzz_1,   \
                             ta_xxx_yyy_0,   \
                             ta_xxx_yyy_1,   \
                             ta_xxx_yyz_0,   \
                             ta_xxx_yyz_1,   \
                             ta_xxx_yzz_0,   \
                             ta_xxx_yzz_1,   \
                             ta_xxx_zzz_0,   \
                             ta_xxx_zzz_1,   \
                             ta_xxxx_xx_0,   \
                             ta_xxxx_xx_1,   \
                             ta_xxxx_xxx_0,  \
                             ta_xxxx_xxx_1,  \
                             ta_xxxx_xxy_0,  \
                             ta_xxxx_xxy_1,  \
                             ta_xxxx_xxz_0,  \
                             ta_xxxx_xxz_1,  \
                             ta_xxxx_xy_0,   \
                             ta_xxxx_xy_1,   \
                             ta_xxxx_xyy_0,  \
                             ta_xxxx_xyy_1,  \
                             ta_xxxx_xyz_0,  \
                             ta_xxxx_xyz_1,  \
                             ta_xxxx_xz_0,   \
                             ta_xxxx_xz_1,   \
                             ta_xxxx_xzz_0,  \
                             ta_xxxx_xzz_1,  \
                             ta_xxxx_yy_0,   \
                             ta_xxxx_yy_1,   \
                             ta_xxxx_yyy_0,  \
                             ta_xxxx_yyy_1,  \
                             ta_xxxx_yyz_0,  \
                             ta_xxxx_yyz_1,  \
                             ta_xxxx_yz_0,   \
                             ta_xxxx_yz_1,   \
                             ta_xxxx_yzz_0,  \
                             ta_xxxx_yzz_1,  \
                             ta_xxxx_zz_0,   \
                             ta_xxxx_zz_1,   \
                             ta_xxxx_zzz_0,  \
                             ta_xxxx_zzz_1,  \
                             ta_xxxxx_xxx_0, \
                             ta_xxxxx_xxy_0, \
                             ta_xxxxx_xxz_0, \
                             ta_xxxxx_xyy_0, \
                             ta_xxxxx_xyz_0, \
                             ta_xxxxx_xzz_0, \
                             ta_xxxxx_yyy_0, \
                             ta_xxxxx_yyz_0, \
                             ta_xxxxx_yzz_0, \
                             ta_xxxxx_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxx_xxx_0[i] = 4.0 * ta_xxx_xxx_0[i] * fe_0 - 4.0 * ta_xxx_xxx_1[i] * fe_0 + 3.0 * ta_xxxx_xx_0[i] * fe_0 -
                            3.0 * ta_xxxx_xx_1[i] * fe_0 + ta_xxxx_xxx_0[i] * pa_x[i] - ta_xxxx_xxx_1[i] * pc_x[i];

        ta_xxxxx_xxy_0[i] = 4.0 * ta_xxx_xxy_0[i] * fe_0 - 4.0 * ta_xxx_xxy_1[i] * fe_0 + 2.0 * ta_xxxx_xy_0[i] * fe_0 -
                            2.0 * ta_xxxx_xy_1[i] * fe_0 + ta_xxxx_xxy_0[i] * pa_x[i] - ta_xxxx_xxy_1[i] * pc_x[i];

        ta_xxxxx_xxz_0[i] = 4.0 * ta_xxx_xxz_0[i] * fe_0 - 4.0 * ta_xxx_xxz_1[i] * fe_0 + 2.0 * ta_xxxx_xz_0[i] * fe_0 -
                            2.0 * ta_xxxx_xz_1[i] * fe_0 + ta_xxxx_xxz_0[i] * pa_x[i] - ta_xxxx_xxz_1[i] * pc_x[i];

        ta_xxxxx_xyy_0[i] = 4.0 * ta_xxx_xyy_0[i] * fe_0 - 4.0 * ta_xxx_xyy_1[i] * fe_0 + ta_xxxx_yy_0[i] * fe_0 - ta_xxxx_yy_1[i] * fe_0 +
                            ta_xxxx_xyy_0[i] * pa_x[i] - ta_xxxx_xyy_1[i] * pc_x[i];

        ta_xxxxx_xyz_0[i] = 4.0 * ta_xxx_xyz_0[i] * fe_0 - 4.0 * ta_xxx_xyz_1[i] * fe_0 + ta_xxxx_yz_0[i] * fe_0 - ta_xxxx_yz_1[i] * fe_0 +
                            ta_xxxx_xyz_0[i] * pa_x[i] - ta_xxxx_xyz_1[i] * pc_x[i];

        ta_xxxxx_xzz_0[i] = 4.0 * ta_xxx_xzz_0[i] * fe_0 - 4.0 * ta_xxx_xzz_1[i] * fe_0 + ta_xxxx_zz_0[i] * fe_0 - ta_xxxx_zz_1[i] * fe_0 +
                            ta_xxxx_xzz_0[i] * pa_x[i] - ta_xxxx_xzz_1[i] * pc_x[i];

        ta_xxxxx_yyy_0[i] = 4.0 * ta_xxx_yyy_0[i] * fe_0 - 4.0 * ta_xxx_yyy_1[i] * fe_0 + ta_xxxx_yyy_0[i] * pa_x[i] - ta_xxxx_yyy_1[i] * pc_x[i];

        ta_xxxxx_yyz_0[i] = 4.0 * ta_xxx_yyz_0[i] * fe_0 - 4.0 * ta_xxx_yyz_1[i] * fe_0 + ta_xxxx_yyz_0[i] * pa_x[i] - ta_xxxx_yyz_1[i] * pc_x[i];

        ta_xxxxx_yzz_0[i] = 4.0 * ta_xxx_yzz_0[i] * fe_0 - 4.0 * ta_xxx_yzz_1[i] * fe_0 + ta_xxxx_yzz_0[i] * pa_x[i] - ta_xxxx_yzz_1[i] * pc_x[i];

        ta_xxxxx_zzz_0[i] = 4.0 * ta_xxx_zzz_0[i] * fe_0 - 4.0 * ta_xxx_zzz_1[i] * fe_0 + ta_xxxx_zzz_0[i] * pa_x[i] - ta_xxxx_zzz_1[i] * pc_x[i];
    }

    // Set up 10-20 components of targeted buffer : HF

    auto ta_xxxxy_xxx_0 = pbuffer.data(idx_npot_0_hf + 10);

    auto ta_xxxxy_xxy_0 = pbuffer.data(idx_npot_0_hf + 11);

    auto ta_xxxxy_xxz_0 = pbuffer.data(idx_npot_0_hf + 12);

    auto ta_xxxxy_xyy_0 = pbuffer.data(idx_npot_0_hf + 13);

    auto ta_xxxxy_xyz_0 = pbuffer.data(idx_npot_0_hf + 14);

    auto ta_xxxxy_xzz_0 = pbuffer.data(idx_npot_0_hf + 15);

    auto ta_xxxxy_yyy_0 = pbuffer.data(idx_npot_0_hf + 16);

    auto ta_xxxxy_yyz_0 = pbuffer.data(idx_npot_0_hf + 17);

    auto ta_xxxxy_yzz_0 = pbuffer.data(idx_npot_0_hf + 18);

    auto ta_xxxxy_zzz_0 = pbuffer.data(idx_npot_0_hf + 19);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta_xxxx_xx_0,   \
                             ta_xxxx_xx_1,   \
                             ta_xxxx_xxx_0,  \
                             ta_xxxx_xxx_1,  \
                             ta_xxxx_xxy_0,  \
                             ta_xxxx_xxy_1,  \
                             ta_xxxx_xxz_0,  \
                             ta_xxxx_xxz_1,  \
                             ta_xxxx_xy_0,   \
                             ta_xxxx_xy_1,   \
                             ta_xxxx_xyy_0,  \
                             ta_xxxx_xyy_1,  \
                             ta_xxxx_xyz_0,  \
                             ta_xxxx_xyz_1,  \
                             ta_xxxx_xz_0,   \
                             ta_xxxx_xz_1,   \
                             ta_xxxx_xzz_0,  \
                             ta_xxxx_xzz_1,  \
                             ta_xxxx_zzz_0,  \
                             ta_xxxx_zzz_1,  \
                             ta_xxxxy_xxx_0, \
                             ta_xxxxy_xxy_0, \
                             ta_xxxxy_xxz_0, \
                             ta_xxxxy_xyy_0, \
                             ta_xxxxy_xyz_0, \
                             ta_xxxxy_xzz_0, \
                             ta_xxxxy_yyy_0, \
                             ta_xxxxy_yyz_0, \
                             ta_xxxxy_yzz_0, \
                             ta_xxxxy_zzz_0, \
                             ta_xxxy_yyy_0,  \
                             ta_xxxy_yyy_1,  \
                             ta_xxxy_yyz_0,  \
                             ta_xxxy_yyz_1,  \
                             ta_xxxy_yzz_0,  \
                             ta_xxxy_yzz_1,  \
                             ta_xxy_yyy_0,   \
                             ta_xxy_yyy_1,   \
                             ta_xxy_yyz_0,   \
                             ta_xxy_yyz_1,   \
                             ta_xxy_yzz_0,   \
                             ta_xxy_yzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxy_xxx_0[i] = ta_xxxx_xxx_0[i] * pa_y[i] - ta_xxxx_xxx_1[i] * pc_y[i];

        ta_xxxxy_xxy_0[i] = ta_xxxx_xx_0[i] * fe_0 - ta_xxxx_xx_1[i] * fe_0 + ta_xxxx_xxy_0[i] * pa_y[i] - ta_xxxx_xxy_1[i] * pc_y[i];

        ta_xxxxy_xxz_0[i] = ta_xxxx_xxz_0[i] * pa_y[i] - ta_xxxx_xxz_1[i] * pc_y[i];

        ta_xxxxy_xyy_0[i] = 2.0 * ta_xxxx_xy_0[i] * fe_0 - 2.0 * ta_xxxx_xy_1[i] * fe_0 + ta_xxxx_xyy_0[i] * pa_y[i] - ta_xxxx_xyy_1[i] * pc_y[i];

        ta_xxxxy_xyz_0[i] = ta_xxxx_xz_0[i] * fe_0 - ta_xxxx_xz_1[i] * fe_0 + ta_xxxx_xyz_0[i] * pa_y[i] - ta_xxxx_xyz_1[i] * pc_y[i];

        ta_xxxxy_xzz_0[i] = ta_xxxx_xzz_0[i] * pa_y[i] - ta_xxxx_xzz_1[i] * pc_y[i];

        ta_xxxxy_yyy_0[i] = 3.0 * ta_xxy_yyy_0[i] * fe_0 - 3.0 * ta_xxy_yyy_1[i] * fe_0 + ta_xxxy_yyy_0[i] * pa_x[i] - ta_xxxy_yyy_1[i] * pc_x[i];

        ta_xxxxy_yyz_0[i] = 3.0 * ta_xxy_yyz_0[i] * fe_0 - 3.0 * ta_xxy_yyz_1[i] * fe_0 + ta_xxxy_yyz_0[i] * pa_x[i] - ta_xxxy_yyz_1[i] * pc_x[i];

        ta_xxxxy_yzz_0[i] = 3.0 * ta_xxy_yzz_0[i] * fe_0 - 3.0 * ta_xxy_yzz_1[i] * fe_0 + ta_xxxy_yzz_0[i] * pa_x[i] - ta_xxxy_yzz_1[i] * pc_x[i];

        ta_xxxxy_zzz_0[i] = ta_xxxx_zzz_0[i] * pa_y[i] - ta_xxxx_zzz_1[i] * pc_y[i];
    }

    // Set up 20-30 components of targeted buffer : HF

    auto ta_xxxxz_xxx_0 = pbuffer.data(idx_npot_0_hf + 20);

    auto ta_xxxxz_xxy_0 = pbuffer.data(idx_npot_0_hf + 21);

    auto ta_xxxxz_xxz_0 = pbuffer.data(idx_npot_0_hf + 22);

    auto ta_xxxxz_xyy_0 = pbuffer.data(idx_npot_0_hf + 23);

    auto ta_xxxxz_xyz_0 = pbuffer.data(idx_npot_0_hf + 24);

    auto ta_xxxxz_xzz_0 = pbuffer.data(idx_npot_0_hf + 25);

    auto ta_xxxxz_yyy_0 = pbuffer.data(idx_npot_0_hf + 26);

    auto ta_xxxxz_yyz_0 = pbuffer.data(idx_npot_0_hf + 27);

    auto ta_xxxxz_yzz_0 = pbuffer.data(idx_npot_0_hf + 28);

    auto ta_xxxxz_zzz_0 = pbuffer.data(idx_npot_0_hf + 29);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta_xxxx_xx_0,   \
                             ta_xxxx_xx_1,   \
                             ta_xxxx_xxx_0,  \
                             ta_xxxx_xxx_1,  \
                             ta_xxxx_xxy_0,  \
                             ta_xxxx_xxy_1,  \
                             ta_xxxx_xxz_0,  \
                             ta_xxxx_xxz_1,  \
                             ta_xxxx_xy_0,   \
                             ta_xxxx_xy_1,   \
                             ta_xxxx_xyy_0,  \
                             ta_xxxx_xyy_1,  \
                             ta_xxxx_xyz_0,  \
                             ta_xxxx_xyz_1,  \
                             ta_xxxx_xz_0,   \
                             ta_xxxx_xz_1,   \
                             ta_xxxx_xzz_0,  \
                             ta_xxxx_xzz_1,  \
                             ta_xxxx_yyy_0,  \
                             ta_xxxx_yyy_1,  \
                             ta_xxxxz_xxx_0, \
                             ta_xxxxz_xxy_0, \
                             ta_xxxxz_xxz_0, \
                             ta_xxxxz_xyy_0, \
                             ta_xxxxz_xyz_0, \
                             ta_xxxxz_xzz_0, \
                             ta_xxxxz_yyy_0, \
                             ta_xxxxz_yyz_0, \
                             ta_xxxxz_yzz_0, \
                             ta_xxxxz_zzz_0, \
                             ta_xxxz_yyz_0,  \
                             ta_xxxz_yyz_1,  \
                             ta_xxxz_yzz_0,  \
                             ta_xxxz_yzz_1,  \
                             ta_xxxz_zzz_0,  \
                             ta_xxxz_zzz_1,  \
                             ta_xxz_yyz_0,   \
                             ta_xxz_yyz_1,   \
                             ta_xxz_yzz_0,   \
                             ta_xxz_yzz_1,   \
                             ta_xxz_zzz_0,   \
                             ta_xxz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxz_xxx_0[i] = ta_xxxx_xxx_0[i] * pa_z[i] - ta_xxxx_xxx_1[i] * pc_z[i];

        ta_xxxxz_xxy_0[i] = ta_xxxx_xxy_0[i] * pa_z[i] - ta_xxxx_xxy_1[i] * pc_z[i];

        ta_xxxxz_xxz_0[i] = ta_xxxx_xx_0[i] * fe_0 - ta_xxxx_xx_1[i] * fe_0 + ta_xxxx_xxz_0[i] * pa_z[i] - ta_xxxx_xxz_1[i] * pc_z[i];

        ta_xxxxz_xyy_0[i] = ta_xxxx_xyy_0[i] * pa_z[i] - ta_xxxx_xyy_1[i] * pc_z[i];

        ta_xxxxz_xyz_0[i] = ta_xxxx_xy_0[i] * fe_0 - ta_xxxx_xy_1[i] * fe_0 + ta_xxxx_xyz_0[i] * pa_z[i] - ta_xxxx_xyz_1[i] * pc_z[i];

        ta_xxxxz_xzz_0[i] = 2.0 * ta_xxxx_xz_0[i] * fe_0 - 2.0 * ta_xxxx_xz_1[i] * fe_0 + ta_xxxx_xzz_0[i] * pa_z[i] - ta_xxxx_xzz_1[i] * pc_z[i];

        ta_xxxxz_yyy_0[i] = ta_xxxx_yyy_0[i] * pa_z[i] - ta_xxxx_yyy_1[i] * pc_z[i];

        ta_xxxxz_yyz_0[i] = 3.0 * ta_xxz_yyz_0[i] * fe_0 - 3.0 * ta_xxz_yyz_1[i] * fe_0 + ta_xxxz_yyz_0[i] * pa_x[i] - ta_xxxz_yyz_1[i] * pc_x[i];

        ta_xxxxz_yzz_0[i] = 3.0 * ta_xxz_yzz_0[i] * fe_0 - 3.0 * ta_xxz_yzz_1[i] * fe_0 + ta_xxxz_yzz_0[i] * pa_x[i] - ta_xxxz_yzz_1[i] * pc_x[i];

        ta_xxxxz_zzz_0[i] = 3.0 * ta_xxz_zzz_0[i] * fe_0 - 3.0 * ta_xxz_zzz_1[i] * fe_0 + ta_xxxz_zzz_0[i] * pa_x[i] - ta_xxxz_zzz_1[i] * pc_x[i];
    }

    // Set up 30-40 components of targeted buffer : HF

    auto ta_xxxyy_xxx_0 = pbuffer.data(idx_npot_0_hf + 30);

    auto ta_xxxyy_xxy_0 = pbuffer.data(idx_npot_0_hf + 31);

    auto ta_xxxyy_xxz_0 = pbuffer.data(idx_npot_0_hf + 32);

    auto ta_xxxyy_xyy_0 = pbuffer.data(idx_npot_0_hf + 33);

    auto ta_xxxyy_xyz_0 = pbuffer.data(idx_npot_0_hf + 34);

    auto ta_xxxyy_xzz_0 = pbuffer.data(idx_npot_0_hf + 35);

    auto ta_xxxyy_yyy_0 = pbuffer.data(idx_npot_0_hf + 36);

    auto ta_xxxyy_yyz_0 = pbuffer.data(idx_npot_0_hf + 37);

    auto ta_xxxyy_yzz_0 = pbuffer.data(idx_npot_0_hf + 38);

    auto ta_xxxyy_zzz_0 = pbuffer.data(idx_npot_0_hf + 39);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta_xxx_xxx_0,   \
                             ta_xxx_xxx_1,   \
                             ta_xxx_xxz_0,   \
                             ta_xxx_xxz_1,   \
                             ta_xxx_xzz_0,   \
                             ta_xxx_xzz_1,   \
                             ta_xxxy_xxx_0,  \
                             ta_xxxy_xxx_1,  \
                             ta_xxxy_xxz_0,  \
                             ta_xxxy_xxz_1,  \
                             ta_xxxy_xzz_0,  \
                             ta_xxxy_xzz_1,  \
                             ta_xxxyy_xxx_0, \
                             ta_xxxyy_xxy_0, \
                             ta_xxxyy_xxz_0, \
                             ta_xxxyy_xyy_0, \
                             ta_xxxyy_xyz_0, \
                             ta_xxxyy_xzz_0, \
                             ta_xxxyy_yyy_0, \
                             ta_xxxyy_yyz_0, \
                             ta_xxxyy_yzz_0, \
                             ta_xxxyy_zzz_0, \
                             ta_xxyy_xxy_0,  \
                             ta_xxyy_xxy_1,  \
                             ta_xxyy_xy_0,   \
                             ta_xxyy_xy_1,   \
                             ta_xxyy_xyy_0,  \
                             ta_xxyy_xyy_1,  \
                             ta_xxyy_xyz_0,  \
                             ta_xxyy_xyz_1,  \
                             ta_xxyy_yy_0,   \
                             ta_xxyy_yy_1,   \
                             ta_xxyy_yyy_0,  \
                             ta_xxyy_yyy_1,  \
                             ta_xxyy_yyz_0,  \
                             ta_xxyy_yyz_1,  \
                             ta_xxyy_yz_0,   \
                             ta_xxyy_yz_1,   \
                             ta_xxyy_yzz_0,  \
                             ta_xxyy_yzz_1,  \
                             ta_xxyy_zzz_0,  \
                             ta_xxyy_zzz_1,  \
                             ta_xyy_xxy_0,   \
                             ta_xyy_xxy_1,   \
                             ta_xyy_xyy_0,   \
                             ta_xyy_xyy_1,   \
                             ta_xyy_xyz_0,   \
                             ta_xyy_xyz_1,   \
                             ta_xyy_yyy_0,   \
                             ta_xyy_yyy_1,   \
                             ta_xyy_yyz_0,   \
                             ta_xyy_yyz_1,   \
                             ta_xyy_yzz_0,   \
                             ta_xyy_yzz_1,   \
                             ta_xyy_zzz_0,   \
                             ta_xyy_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyy_xxx_0[i] = ta_xxx_xxx_0[i] * fe_0 - ta_xxx_xxx_1[i] * fe_0 + ta_xxxy_xxx_0[i] * pa_y[i] - ta_xxxy_xxx_1[i] * pc_y[i];

        ta_xxxyy_xxy_0[i] = 2.0 * ta_xyy_xxy_0[i] * fe_0 - 2.0 * ta_xyy_xxy_1[i] * fe_0 + 2.0 * ta_xxyy_xy_0[i] * fe_0 -
                            2.0 * ta_xxyy_xy_1[i] * fe_0 + ta_xxyy_xxy_0[i] * pa_x[i] - ta_xxyy_xxy_1[i] * pc_x[i];

        ta_xxxyy_xxz_0[i] = ta_xxx_xxz_0[i] * fe_0 - ta_xxx_xxz_1[i] * fe_0 + ta_xxxy_xxz_0[i] * pa_y[i] - ta_xxxy_xxz_1[i] * pc_y[i];

        ta_xxxyy_xyy_0[i] = 2.0 * ta_xyy_xyy_0[i] * fe_0 - 2.0 * ta_xyy_xyy_1[i] * fe_0 + ta_xxyy_yy_0[i] * fe_0 - ta_xxyy_yy_1[i] * fe_0 +
                            ta_xxyy_xyy_0[i] * pa_x[i] - ta_xxyy_xyy_1[i] * pc_x[i];

        ta_xxxyy_xyz_0[i] = 2.0 * ta_xyy_xyz_0[i] * fe_0 - 2.0 * ta_xyy_xyz_1[i] * fe_0 + ta_xxyy_yz_0[i] * fe_0 - ta_xxyy_yz_1[i] * fe_0 +
                            ta_xxyy_xyz_0[i] * pa_x[i] - ta_xxyy_xyz_1[i] * pc_x[i];

        ta_xxxyy_xzz_0[i] = ta_xxx_xzz_0[i] * fe_0 - ta_xxx_xzz_1[i] * fe_0 + ta_xxxy_xzz_0[i] * pa_y[i] - ta_xxxy_xzz_1[i] * pc_y[i];

        ta_xxxyy_yyy_0[i] = 2.0 * ta_xyy_yyy_0[i] * fe_0 - 2.0 * ta_xyy_yyy_1[i] * fe_0 + ta_xxyy_yyy_0[i] * pa_x[i] - ta_xxyy_yyy_1[i] * pc_x[i];

        ta_xxxyy_yyz_0[i] = 2.0 * ta_xyy_yyz_0[i] * fe_0 - 2.0 * ta_xyy_yyz_1[i] * fe_0 + ta_xxyy_yyz_0[i] * pa_x[i] - ta_xxyy_yyz_1[i] * pc_x[i];

        ta_xxxyy_yzz_0[i] = 2.0 * ta_xyy_yzz_0[i] * fe_0 - 2.0 * ta_xyy_yzz_1[i] * fe_0 + ta_xxyy_yzz_0[i] * pa_x[i] - ta_xxyy_yzz_1[i] * pc_x[i];

        ta_xxxyy_zzz_0[i] = 2.0 * ta_xyy_zzz_0[i] * fe_0 - 2.0 * ta_xyy_zzz_1[i] * fe_0 + ta_xxyy_zzz_0[i] * pa_x[i] - ta_xxyy_zzz_1[i] * pc_x[i];
    }

    // Set up 40-50 components of targeted buffer : HF

    auto ta_xxxyz_xxx_0 = pbuffer.data(idx_npot_0_hf + 40);

    auto ta_xxxyz_xxy_0 = pbuffer.data(idx_npot_0_hf + 41);

    auto ta_xxxyz_xxz_0 = pbuffer.data(idx_npot_0_hf + 42);

    auto ta_xxxyz_xyy_0 = pbuffer.data(idx_npot_0_hf + 43);

    auto ta_xxxyz_xyz_0 = pbuffer.data(idx_npot_0_hf + 44);

    auto ta_xxxyz_xzz_0 = pbuffer.data(idx_npot_0_hf + 45);

    auto ta_xxxyz_yyy_0 = pbuffer.data(idx_npot_0_hf + 46);

    auto ta_xxxyz_yyz_0 = pbuffer.data(idx_npot_0_hf + 47);

    auto ta_xxxyz_yzz_0 = pbuffer.data(idx_npot_0_hf + 48);

    auto ta_xxxyz_zzz_0 = pbuffer.data(idx_npot_0_hf + 49);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             pc_x,           \
                             pc_y,           \
                             pc_z,           \
                             ta_xxxy_xxy_0,  \
                             ta_xxxy_xxy_1,  \
                             ta_xxxy_xyy_0,  \
                             ta_xxxy_xyy_1,  \
                             ta_xxxy_yyy_0,  \
                             ta_xxxy_yyy_1,  \
                             ta_xxxyz_xxx_0, \
                             ta_xxxyz_xxy_0, \
                             ta_xxxyz_xxz_0, \
                             ta_xxxyz_xyy_0, \
                             ta_xxxyz_xyz_0, \
                             ta_xxxyz_xzz_0, \
                             ta_xxxyz_yyy_0, \
                             ta_xxxyz_yyz_0, \
                             ta_xxxyz_yzz_0, \
                             ta_xxxyz_zzz_0, \
                             ta_xxxz_xxx_0,  \
                             ta_xxxz_xxx_1,  \
                             ta_xxxz_xxz_0,  \
                             ta_xxxz_xxz_1,  \
                             ta_xxxz_xyz_0,  \
                             ta_xxxz_xyz_1,  \
                             ta_xxxz_xz_0,   \
                             ta_xxxz_xz_1,   \
                             ta_xxxz_xzz_0,  \
                             ta_xxxz_xzz_1,  \
                             ta_xxxz_zzz_0,  \
                             ta_xxxz_zzz_1,  \
                             ta_xxyz_yyz_0,  \
                             ta_xxyz_yyz_1,  \
                             ta_xxyz_yzz_0,  \
                             ta_xxyz_yzz_1,  \
                             ta_xyz_yyz_0,   \
                             ta_xyz_yyz_1,   \
                             ta_xyz_yzz_0,   \
                             ta_xyz_yzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyz_xxx_0[i] = ta_xxxz_xxx_0[i] * pa_y[i] - ta_xxxz_xxx_1[i] * pc_y[i];

        ta_xxxyz_xxy_0[i] = ta_xxxy_xxy_0[i] * pa_z[i] - ta_xxxy_xxy_1[i] * pc_z[i];

        ta_xxxyz_xxz_0[i] = ta_xxxz_xxz_0[i] * pa_y[i] - ta_xxxz_xxz_1[i] * pc_y[i];

        ta_xxxyz_xyy_0[i] = ta_xxxy_xyy_0[i] * pa_z[i] - ta_xxxy_xyy_1[i] * pc_z[i];

        ta_xxxyz_xyz_0[i] = ta_xxxz_xz_0[i] * fe_0 - ta_xxxz_xz_1[i] * fe_0 + ta_xxxz_xyz_0[i] * pa_y[i] - ta_xxxz_xyz_1[i] * pc_y[i];

        ta_xxxyz_xzz_0[i] = ta_xxxz_xzz_0[i] * pa_y[i] - ta_xxxz_xzz_1[i] * pc_y[i];

        ta_xxxyz_yyy_0[i] = ta_xxxy_yyy_0[i] * pa_z[i] - ta_xxxy_yyy_1[i] * pc_z[i];

        ta_xxxyz_yyz_0[i] = 2.0 * ta_xyz_yyz_0[i] * fe_0 - 2.0 * ta_xyz_yyz_1[i] * fe_0 + ta_xxyz_yyz_0[i] * pa_x[i] - ta_xxyz_yyz_1[i] * pc_x[i];

        ta_xxxyz_yzz_0[i] = 2.0 * ta_xyz_yzz_0[i] * fe_0 - 2.0 * ta_xyz_yzz_1[i] * fe_0 + ta_xxyz_yzz_0[i] * pa_x[i] - ta_xxyz_yzz_1[i] * pc_x[i];

        ta_xxxyz_zzz_0[i] = ta_xxxz_zzz_0[i] * pa_y[i] - ta_xxxz_zzz_1[i] * pc_y[i];
    }

    // Set up 50-60 components of targeted buffer : HF

    auto ta_xxxzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 50);

    auto ta_xxxzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 51);

    auto ta_xxxzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 52);

    auto ta_xxxzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 53);

    auto ta_xxxzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 54);

    auto ta_xxxzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 55);

    auto ta_xxxzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 56);

    auto ta_xxxzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 57);

    auto ta_xxxzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 58);

    auto ta_xxxzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 59);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta_xxx_xxx_0,   \
                             ta_xxx_xxx_1,   \
                             ta_xxx_xxy_0,   \
                             ta_xxx_xxy_1,   \
                             ta_xxx_xyy_0,   \
                             ta_xxx_xyy_1,   \
                             ta_xxxz_xxx_0,  \
                             ta_xxxz_xxx_1,  \
                             ta_xxxz_xxy_0,  \
                             ta_xxxz_xxy_1,  \
                             ta_xxxz_xyy_0,  \
                             ta_xxxz_xyy_1,  \
                             ta_xxxzz_xxx_0, \
                             ta_xxxzz_xxy_0, \
                             ta_xxxzz_xxz_0, \
                             ta_xxxzz_xyy_0, \
                             ta_xxxzz_xyz_0, \
                             ta_xxxzz_xzz_0, \
                             ta_xxxzz_yyy_0, \
                             ta_xxxzz_yyz_0, \
                             ta_xxxzz_yzz_0, \
                             ta_xxxzz_zzz_0, \
                             ta_xxzz_xxz_0,  \
                             ta_xxzz_xxz_1,  \
                             ta_xxzz_xyz_0,  \
                             ta_xxzz_xyz_1,  \
                             ta_xxzz_xz_0,   \
                             ta_xxzz_xz_1,   \
                             ta_xxzz_xzz_0,  \
                             ta_xxzz_xzz_1,  \
                             ta_xxzz_yyy_0,  \
                             ta_xxzz_yyy_1,  \
                             ta_xxzz_yyz_0,  \
                             ta_xxzz_yyz_1,  \
                             ta_xxzz_yz_0,   \
                             ta_xxzz_yz_1,   \
                             ta_xxzz_yzz_0,  \
                             ta_xxzz_yzz_1,  \
                             ta_xxzz_zz_0,   \
                             ta_xxzz_zz_1,   \
                             ta_xxzz_zzz_0,  \
                             ta_xxzz_zzz_1,  \
                             ta_xzz_xxz_0,   \
                             ta_xzz_xxz_1,   \
                             ta_xzz_xyz_0,   \
                             ta_xzz_xyz_1,   \
                             ta_xzz_xzz_0,   \
                             ta_xzz_xzz_1,   \
                             ta_xzz_yyy_0,   \
                             ta_xzz_yyy_1,   \
                             ta_xzz_yyz_0,   \
                             ta_xzz_yyz_1,   \
                             ta_xzz_yzz_0,   \
                             ta_xzz_yzz_1,   \
                             ta_xzz_zzz_0,   \
                             ta_xzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxzz_xxx_0[i] = ta_xxx_xxx_0[i] * fe_0 - ta_xxx_xxx_1[i] * fe_0 + ta_xxxz_xxx_0[i] * pa_z[i] - ta_xxxz_xxx_1[i] * pc_z[i];

        ta_xxxzz_xxy_0[i] = ta_xxx_xxy_0[i] * fe_0 - ta_xxx_xxy_1[i] * fe_0 + ta_xxxz_xxy_0[i] * pa_z[i] - ta_xxxz_xxy_1[i] * pc_z[i];

        ta_xxxzz_xxz_0[i] = 2.0 * ta_xzz_xxz_0[i] * fe_0 - 2.0 * ta_xzz_xxz_1[i] * fe_0 + 2.0 * ta_xxzz_xz_0[i] * fe_0 -
                            2.0 * ta_xxzz_xz_1[i] * fe_0 + ta_xxzz_xxz_0[i] * pa_x[i] - ta_xxzz_xxz_1[i] * pc_x[i];

        ta_xxxzz_xyy_0[i] = ta_xxx_xyy_0[i] * fe_0 - ta_xxx_xyy_1[i] * fe_0 + ta_xxxz_xyy_0[i] * pa_z[i] - ta_xxxz_xyy_1[i] * pc_z[i];

        ta_xxxzz_xyz_0[i] = 2.0 * ta_xzz_xyz_0[i] * fe_0 - 2.0 * ta_xzz_xyz_1[i] * fe_0 + ta_xxzz_yz_0[i] * fe_0 - ta_xxzz_yz_1[i] * fe_0 +
                            ta_xxzz_xyz_0[i] * pa_x[i] - ta_xxzz_xyz_1[i] * pc_x[i];

        ta_xxxzz_xzz_0[i] = 2.0 * ta_xzz_xzz_0[i] * fe_0 - 2.0 * ta_xzz_xzz_1[i] * fe_0 + ta_xxzz_zz_0[i] * fe_0 - ta_xxzz_zz_1[i] * fe_0 +
                            ta_xxzz_xzz_0[i] * pa_x[i] - ta_xxzz_xzz_1[i] * pc_x[i];

        ta_xxxzz_yyy_0[i] = 2.0 * ta_xzz_yyy_0[i] * fe_0 - 2.0 * ta_xzz_yyy_1[i] * fe_0 + ta_xxzz_yyy_0[i] * pa_x[i] - ta_xxzz_yyy_1[i] * pc_x[i];

        ta_xxxzz_yyz_0[i] = 2.0 * ta_xzz_yyz_0[i] * fe_0 - 2.0 * ta_xzz_yyz_1[i] * fe_0 + ta_xxzz_yyz_0[i] * pa_x[i] - ta_xxzz_yyz_1[i] * pc_x[i];

        ta_xxxzz_yzz_0[i] = 2.0 * ta_xzz_yzz_0[i] * fe_0 - 2.0 * ta_xzz_yzz_1[i] * fe_0 + ta_xxzz_yzz_0[i] * pa_x[i] - ta_xxzz_yzz_1[i] * pc_x[i];

        ta_xxxzz_zzz_0[i] = 2.0 * ta_xzz_zzz_0[i] * fe_0 - 2.0 * ta_xzz_zzz_1[i] * fe_0 + ta_xxzz_zzz_0[i] * pa_x[i] - ta_xxzz_zzz_1[i] * pc_x[i];
    }

    // Set up 60-70 components of targeted buffer : HF

    auto ta_xxyyy_xxx_0 = pbuffer.data(idx_npot_0_hf + 60);

    auto ta_xxyyy_xxy_0 = pbuffer.data(idx_npot_0_hf + 61);

    auto ta_xxyyy_xxz_0 = pbuffer.data(idx_npot_0_hf + 62);

    auto ta_xxyyy_xyy_0 = pbuffer.data(idx_npot_0_hf + 63);

    auto ta_xxyyy_xyz_0 = pbuffer.data(idx_npot_0_hf + 64);

    auto ta_xxyyy_xzz_0 = pbuffer.data(idx_npot_0_hf + 65);

    auto ta_xxyyy_yyy_0 = pbuffer.data(idx_npot_0_hf + 66);

    auto ta_xxyyy_yyz_0 = pbuffer.data(idx_npot_0_hf + 67);

    auto ta_xxyyy_yzz_0 = pbuffer.data(idx_npot_0_hf + 68);

    auto ta_xxyyy_zzz_0 = pbuffer.data(idx_npot_0_hf + 69);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta_xxy_xxx_0,   \
                             ta_xxy_xxx_1,   \
                             ta_xxy_xxz_0,   \
                             ta_xxy_xxz_1,   \
                             ta_xxy_xzz_0,   \
                             ta_xxy_xzz_1,   \
                             ta_xxyy_xxx_0,  \
                             ta_xxyy_xxx_1,  \
                             ta_xxyy_xxz_0,  \
                             ta_xxyy_xxz_1,  \
                             ta_xxyy_xzz_0,  \
                             ta_xxyy_xzz_1,  \
                             ta_xxyyy_xxx_0, \
                             ta_xxyyy_xxy_0, \
                             ta_xxyyy_xxz_0, \
                             ta_xxyyy_xyy_0, \
                             ta_xxyyy_xyz_0, \
                             ta_xxyyy_xzz_0, \
                             ta_xxyyy_yyy_0, \
                             ta_xxyyy_yyz_0, \
                             ta_xxyyy_yzz_0, \
                             ta_xxyyy_zzz_0, \
                             ta_xyyy_xxy_0,  \
                             ta_xyyy_xxy_1,  \
                             ta_xyyy_xy_0,   \
                             ta_xyyy_xy_1,   \
                             ta_xyyy_xyy_0,  \
                             ta_xyyy_xyy_1,  \
                             ta_xyyy_xyz_0,  \
                             ta_xyyy_xyz_1,  \
                             ta_xyyy_yy_0,   \
                             ta_xyyy_yy_1,   \
                             ta_xyyy_yyy_0,  \
                             ta_xyyy_yyy_1,  \
                             ta_xyyy_yyz_0,  \
                             ta_xyyy_yyz_1,  \
                             ta_xyyy_yz_0,   \
                             ta_xyyy_yz_1,   \
                             ta_xyyy_yzz_0,  \
                             ta_xyyy_yzz_1,  \
                             ta_xyyy_zzz_0,  \
                             ta_xyyy_zzz_1,  \
                             ta_yyy_xxy_0,   \
                             ta_yyy_xxy_1,   \
                             ta_yyy_xyy_0,   \
                             ta_yyy_xyy_1,   \
                             ta_yyy_xyz_0,   \
                             ta_yyy_xyz_1,   \
                             ta_yyy_yyy_0,   \
                             ta_yyy_yyy_1,   \
                             ta_yyy_yyz_0,   \
                             ta_yyy_yyz_1,   \
                             ta_yyy_yzz_0,   \
                             ta_yyy_yzz_1,   \
                             ta_yyy_zzz_0,   \
                             ta_yyy_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyy_xxx_0[i] = 2.0 * ta_xxy_xxx_0[i] * fe_0 - 2.0 * ta_xxy_xxx_1[i] * fe_0 + ta_xxyy_xxx_0[i] * pa_y[i] - ta_xxyy_xxx_1[i] * pc_y[i];

        ta_xxyyy_xxy_0[i] = ta_yyy_xxy_0[i] * fe_0 - ta_yyy_xxy_1[i] * fe_0 + 2.0 * ta_xyyy_xy_0[i] * fe_0 - 2.0 * ta_xyyy_xy_1[i] * fe_0 +
                            ta_xyyy_xxy_0[i] * pa_x[i] - ta_xyyy_xxy_1[i] * pc_x[i];

        ta_xxyyy_xxz_0[i] = 2.0 * ta_xxy_xxz_0[i] * fe_0 - 2.0 * ta_xxy_xxz_1[i] * fe_0 + ta_xxyy_xxz_0[i] * pa_y[i] - ta_xxyy_xxz_1[i] * pc_y[i];

        ta_xxyyy_xyy_0[i] = ta_yyy_xyy_0[i] * fe_0 - ta_yyy_xyy_1[i] * fe_0 + ta_xyyy_yy_0[i] * fe_0 - ta_xyyy_yy_1[i] * fe_0 +
                            ta_xyyy_xyy_0[i] * pa_x[i] - ta_xyyy_xyy_1[i] * pc_x[i];

        ta_xxyyy_xyz_0[i] = ta_yyy_xyz_0[i] * fe_0 - ta_yyy_xyz_1[i] * fe_0 + ta_xyyy_yz_0[i] * fe_0 - ta_xyyy_yz_1[i] * fe_0 +
                            ta_xyyy_xyz_0[i] * pa_x[i] - ta_xyyy_xyz_1[i] * pc_x[i];

        ta_xxyyy_xzz_0[i] = 2.0 * ta_xxy_xzz_0[i] * fe_0 - 2.0 * ta_xxy_xzz_1[i] * fe_0 + ta_xxyy_xzz_0[i] * pa_y[i] - ta_xxyy_xzz_1[i] * pc_y[i];

        ta_xxyyy_yyy_0[i] = ta_yyy_yyy_0[i] * fe_0 - ta_yyy_yyy_1[i] * fe_0 + ta_xyyy_yyy_0[i] * pa_x[i] - ta_xyyy_yyy_1[i] * pc_x[i];

        ta_xxyyy_yyz_0[i] = ta_yyy_yyz_0[i] * fe_0 - ta_yyy_yyz_1[i] * fe_0 + ta_xyyy_yyz_0[i] * pa_x[i] - ta_xyyy_yyz_1[i] * pc_x[i];

        ta_xxyyy_yzz_0[i] = ta_yyy_yzz_0[i] * fe_0 - ta_yyy_yzz_1[i] * fe_0 + ta_xyyy_yzz_0[i] * pa_x[i] - ta_xyyy_yzz_1[i] * pc_x[i];

        ta_xxyyy_zzz_0[i] = ta_yyy_zzz_0[i] * fe_0 - ta_yyy_zzz_1[i] * fe_0 + ta_xyyy_zzz_0[i] * pa_x[i] - ta_xyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 70-80 components of targeted buffer : HF

    auto ta_xxyyz_xxx_0 = pbuffer.data(idx_npot_0_hf + 70);

    auto ta_xxyyz_xxy_0 = pbuffer.data(idx_npot_0_hf + 71);

    auto ta_xxyyz_xxz_0 = pbuffer.data(idx_npot_0_hf + 72);

    auto ta_xxyyz_xyy_0 = pbuffer.data(idx_npot_0_hf + 73);

    auto ta_xxyyz_xyz_0 = pbuffer.data(idx_npot_0_hf + 74);

    auto ta_xxyyz_xzz_0 = pbuffer.data(idx_npot_0_hf + 75);

    auto ta_xxyyz_yyy_0 = pbuffer.data(idx_npot_0_hf + 76);

    auto ta_xxyyz_yyz_0 = pbuffer.data(idx_npot_0_hf + 77);

    auto ta_xxyyz_yzz_0 = pbuffer.data(idx_npot_0_hf + 78);

    auto ta_xxyyz_zzz_0 = pbuffer.data(idx_npot_0_hf + 79);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             pc_x,           \
                             pc_y,           \
                             pc_z,           \
                             ta_xxyy_xxx_0,  \
                             ta_xxyy_xxx_1,  \
                             ta_xxyy_xxy_0,  \
                             ta_xxyy_xxy_1,  \
                             ta_xxyy_xy_0,   \
                             ta_xxyy_xy_1,   \
                             ta_xxyy_xyy_0,  \
                             ta_xxyy_xyy_1,  \
                             ta_xxyy_xyz_0,  \
                             ta_xxyy_xyz_1,  \
                             ta_xxyy_yyy_0,  \
                             ta_xxyy_yyy_1,  \
                             ta_xxyyz_xxx_0, \
                             ta_xxyyz_xxy_0, \
                             ta_xxyyz_xxz_0, \
                             ta_xxyyz_xyy_0, \
                             ta_xxyyz_xyz_0, \
                             ta_xxyyz_xzz_0, \
                             ta_xxyyz_yyy_0, \
                             ta_xxyyz_yyz_0, \
                             ta_xxyyz_yzz_0, \
                             ta_xxyyz_zzz_0, \
                             ta_xxyz_xxz_0,  \
                             ta_xxyz_xxz_1,  \
                             ta_xxyz_xzz_0,  \
                             ta_xxyz_xzz_1,  \
                             ta_xxz_xxz_0,   \
                             ta_xxz_xxz_1,   \
                             ta_xxz_xzz_0,   \
                             ta_xxz_xzz_1,   \
                             ta_xyyz_yyz_0,  \
                             ta_xyyz_yyz_1,  \
                             ta_xyyz_yzz_0,  \
                             ta_xyyz_yzz_1,  \
                             ta_xyyz_zzz_0,  \
                             ta_xyyz_zzz_1,  \
                             ta_yyz_yyz_0,   \
                             ta_yyz_yyz_1,   \
                             ta_yyz_yzz_0,   \
                             ta_yyz_yzz_1,   \
                             ta_yyz_zzz_0,   \
                             ta_yyz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyz_xxx_0[i] = ta_xxyy_xxx_0[i] * pa_z[i] - ta_xxyy_xxx_1[i] * pc_z[i];

        ta_xxyyz_xxy_0[i] = ta_xxyy_xxy_0[i] * pa_z[i] - ta_xxyy_xxy_1[i] * pc_z[i];

        ta_xxyyz_xxz_0[i] = ta_xxz_xxz_0[i] * fe_0 - ta_xxz_xxz_1[i] * fe_0 + ta_xxyz_xxz_0[i] * pa_y[i] - ta_xxyz_xxz_1[i] * pc_y[i];

        ta_xxyyz_xyy_0[i] = ta_xxyy_xyy_0[i] * pa_z[i] - ta_xxyy_xyy_1[i] * pc_z[i];

        ta_xxyyz_xyz_0[i] = ta_xxyy_xy_0[i] * fe_0 - ta_xxyy_xy_1[i] * fe_0 + ta_xxyy_xyz_0[i] * pa_z[i] - ta_xxyy_xyz_1[i] * pc_z[i];

        ta_xxyyz_xzz_0[i] = ta_xxz_xzz_0[i] * fe_0 - ta_xxz_xzz_1[i] * fe_0 + ta_xxyz_xzz_0[i] * pa_y[i] - ta_xxyz_xzz_1[i] * pc_y[i];

        ta_xxyyz_yyy_0[i] = ta_xxyy_yyy_0[i] * pa_z[i] - ta_xxyy_yyy_1[i] * pc_z[i];

        ta_xxyyz_yyz_0[i] = ta_yyz_yyz_0[i] * fe_0 - ta_yyz_yyz_1[i] * fe_0 + ta_xyyz_yyz_0[i] * pa_x[i] - ta_xyyz_yyz_1[i] * pc_x[i];

        ta_xxyyz_yzz_0[i] = ta_yyz_yzz_0[i] * fe_0 - ta_yyz_yzz_1[i] * fe_0 + ta_xyyz_yzz_0[i] * pa_x[i] - ta_xyyz_yzz_1[i] * pc_x[i];

        ta_xxyyz_zzz_0[i] = ta_yyz_zzz_0[i] * fe_0 - ta_yyz_zzz_1[i] * fe_0 + ta_xyyz_zzz_0[i] * pa_x[i] - ta_xyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 80-90 components of targeted buffer : HF

    auto ta_xxyzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 80);

    auto ta_xxyzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 81);

    auto ta_xxyzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 82);

    auto ta_xxyzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 83);

    auto ta_xxyzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 84);

    auto ta_xxyzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 85);

    auto ta_xxyzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 86);

    auto ta_xxyzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 87);

    auto ta_xxyzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 88);

    auto ta_xxyzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 89);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta_xxyzz_xxx_0, \
                             ta_xxyzz_xxy_0, \
                             ta_xxyzz_xxz_0, \
                             ta_xxyzz_xyy_0, \
                             ta_xxyzz_xyz_0, \
                             ta_xxyzz_xzz_0, \
                             ta_xxyzz_yyy_0, \
                             ta_xxyzz_yyz_0, \
                             ta_xxyzz_yzz_0, \
                             ta_xxyzz_zzz_0, \
                             ta_xxzz_xx_0,   \
                             ta_xxzz_xx_1,   \
                             ta_xxzz_xxx_0,  \
                             ta_xxzz_xxx_1,  \
                             ta_xxzz_xxy_0,  \
                             ta_xxzz_xxy_1,  \
                             ta_xxzz_xxz_0,  \
                             ta_xxzz_xxz_1,  \
                             ta_xxzz_xy_0,   \
                             ta_xxzz_xy_1,   \
                             ta_xxzz_xyy_0,  \
                             ta_xxzz_xyy_1,  \
                             ta_xxzz_xyz_0,  \
                             ta_xxzz_xyz_1,  \
                             ta_xxzz_xz_0,   \
                             ta_xxzz_xz_1,   \
                             ta_xxzz_xzz_0,  \
                             ta_xxzz_xzz_1,  \
                             ta_xxzz_zzz_0,  \
                             ta_xxzz_zzz_1,  \
                             ta_xyzz_yyy_0,  \
                             ta_xyzz_yyy_1,  \
                             ta_xyzz_yyz_0,  \
                             ta_xyzz_yyz_1,  \
                             ta_xyzz_yzz_0,  \
                             ta_xyzz_yzz_1,  \
                             ta_yzz_yyy_0,   \
                             ta_yzz_yyy_1,   \
                             ta_yzz_yyz_0,   \
                             ta_yzz_yyz_1,   \
                             ta_yzz_yzz_0,   \
                             ta_yzz_yzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyzz_xxx_0[i] = ta_xxzz_xxx_0[i] * pa_y[i] - ta_xxzz_xxx_1[i] * pc_y[i];

        ta_xxyzz_xxy_0[i] = ta_xxzz_xx_0[i] * fe_0 - ta_xxzz_xx_1[i] * fe_0 + ta_xxzz_xxy_0[i] * pa_y[i] - ta_xxzz_xxy_1[i] * pc_y[i];

        ta_xxyzz_xxz_0[i] = ta_xxzz_xxz_0[i] * pa_y[i] - ta_xxzz_xxz_1[i] * pc_y[i];

        ta_xxyzz_xyy_0[i] = 2.0 * ta_xxzz_xy_0[i] * fe_0 - 2.0 * ta_xxzz_xy_1[i] * fe_0 + ta_xxzz_xyy_0[i] * pa_y[i] - ta_xxzz_xyy_1[i] * pc_y[i];

        ta_xxyzz_xyz_0[i] = ta_xxzz_xz_0[i] * fe_0 - ta_xxzz_xz_1[i] * fe_0 + ta_xxzz_xyz_0[i] * pa_y[i] - ta_xxzz_xyz_1[i] * pc_y[i];

        ta_xxyzz_xzz_0[i] = ta_xxzz_xzz_0[i] * pa_y[i] - ta_xxzz_xzz_1[i] * pc_y[i];

        ta_xxyzz_yyy_0[i] = ta_yzz_yyy_0[i] * fe_0 - ta_yzz_yyy_1[i] * fe_0 + ta_xyzz_yyy_0[i] * pa_x[i] - ta_xyzz_yyy_1[i] * pc_x[i];

        ta_xxyzz_yyz_0[i] = ta_yzz_yyz_0[i] * fe_0 - ta_yzz_yyz_1[i] * fe_0 + ta_xyzz_yyz_0[i] * pa_x[i] - ta_xyzz_yyz_1[i] * pc_x[i];

        ta_xxyzz_yzz_0[i] = ta_yzz_yzz_0[i] * fe_0 - ta_yzz_yzz_1[i] * fe_0 + ta_xyzz_yzz_0[i] * pa_x[i] - ta_xyzz_yzz_1[i] * pc_x[i];

        ta_xxyzz_zzz_0[i] = ta_xxzz_zzz_0[i] * pa_y[i] - ta_xxzz_zzz_1[i] * pc_y[i];
    }

    // Set up 90-100 components of targeted buffer : HF

    auto ta_xxzzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 90);

    auto ta_xxzzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 91);

    auto ta_xxzzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 92);

    auto ta_xxzzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 93);

    auto ta_xxzzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 94);

    auto ta_xxzzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 95);

    auto ta_xxzzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 96);

    auto ta_xxzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 97);

    auto ta_xxzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 98);

    auto ta_xxzzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 99);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta_xxz_xxx_0,   \
                             ta_xxz_xxx_1,   \
                             ta_xxz_xxy_0,   \
                             ta_xxz_xxy_1,   \
                             ta_xxz_xyy_0,   \
                             ta_xxz_xyy_1,   \
                             ta_xxzz_xxx_0,  \
                             ta_xxzz_xxx_1,  \
                             ta_xxzz_xxy_0,  \
                             ta_xxzz_xxy_1,  \
                             ta_xxzz_xyy_0,  \
                             ta_xxzz_xyy_1,  \
                             ta_xxzzz_xxx_0, \
                             ta_xxzzz_xxy_0, \
                             ta_xxzzz_xxz_0, \
                             ta_xxzzz_xyy_0, \
                             ta_xxzzz_xyz_0, \
                             ta_xxzzz_xzz_0, \
                             ta_xxzzz_yyy_0, \
                             ta_xxzzz_yyz_0, \
                             ta_xxzzz_yzz_0, \
                             ta_xxzzz_zzz_0, \
                             ta_xzzz_xxz_0,  \
                             ta_xzzz_xxz_1,  \
                             ta_xzzz_xyz_0,  \
                             ta_xzzz_xyz_1,  \
                             ta_xzzz_xz_0,   \
                             ta_xzzz_xz_1,   \
                             ta_xzzz_xzz_0,  \
                             ta_xzzz_xzz_1,  \
                             ta_xzzz_yyy_0,  \
                             ta_xzzz_yyy_1,  \
                             ta_xzzz_yyz_0,  \
                             ta_xzzz_yyz_1,  \
                             ta_xzzz_yz_0,   \
                             ta_xzzz_yz_1,   \
                             ta_xzzz_yzz_0,  \
                             ta_xzzz_yzz_1,  \
                             ta_xzzz_zz_0,   \
                             ta_xzzz_zz_1,   \
                             ta_xzzz_zzz_0,  \
                             ta_xzzz_zzz_1,  \
                             ta_zzz_xxz_0,   \
                             ta_zzz_xxz_1,   \
                             ta_zzz_xyz_0,   \
                             ta_zzz_xyz_1,   \
                             ta_zzz_xzz_0,   \
                             ta_zzz_xzz_1,   \
                             ta_zzz_yyy_0,   \
                             ta_zzz_yyy_1,   \
                             ta_zzz_yyz_0,   \
                             ta_zzz_yyz_1,   \
                             ta_zzz_yzz_0,   \
                             ta_zzz_yzz_1,   \
                             ta_zzz_zzz_0,   \
                             ta_zzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzzz_xxx_0[i] = 2.0 * ta_xxz_xxx_0[i] * fe_0 - 2.0 * ta_xxz_xxx_1[i] * fe_0 + ta_xxzz_xxx_0[i] * pa_z[i] - ta_xxzz_xxx_1[i] * pc_z[i];

        ta_xxzzz_xxy_0[i] = 2.0 * ta_xxz_xxy_0[i] * fe_0 - 2.0 * ta_xxz_xxy_1[i] * fe_0 + ta_xxzz_xxy_0[i] * pa_z[i] - ta_xxzz_xxy_1[i] * pc_z[i];

        ta_xxzzz_xxz_0[i] = ta_zzz_xxz_0[i] * fe_0 - ta_zzz_xxz_1[i] * fe_0 + 2.0 * ta_xzzz_xz_0[i] * fe_0 - 2.0 * ta_xzzz_xz_1[i] * fe_0 +
                            ta_xzzz_xxz_0[i] * pa_x[i] - ta_xzzz_xxz_1[i] * pc_x[i];

        ta_xxzzz_xyy_0[i] = 2.0 * ta_xxz_xyy_0[i] * fe_0 - 2.0 * ta_xxz_xyy_1[i] * fe_0 + ta_xxzz_xyy_0[i] * pa_z[i] - ta_xxzz_xyy_1[i] * pc_z[i];

        ta_xxzzz_xyz_0[i] = ta_zzz_xyz_0[i] * fe_0 - ta_zzz_xyz_1[i] * fe_0 + ta_xzzz_yz_0[i] * fe_0 - ta_xzzz_yz_1[i] * fe_0 +
                            ta_xzzz_xyz_0[i] * pa_x[i] - ta_xzzz_xyz_1[i] * pc_x[i];

        ta_xxzzz_xzz_0[i] = ta_zzz_xzz_0[i] * fe_0 - ta_zzz_xzz_1[i] * fe_0 + ta_xzzz_zz_0[i] * fe_0 - ta_xzzz_zz_1[i] * fe_0 +
                            ta_xzzz_xzz_0[i] * pa_x[i] - ta_xzzz_xzz_1[i] * pc_x[i];

        ta_xxzzz_yyy_0[i] = ta_zzz_yyy_0[i] * fe_0 - ta_zzz_yyy_1[i] * fe_0 + ta_xzzz_yyy_0[i] * pa_x[i] - ta_xzzz_yyy_1[i] * pc_x[i];

        ta_xxzzz_yyz_0[i] = ta_zzz_yyz_0[i] * fe_0 - ta_zzz_yyz_1[i] * fe_0 + ta_xzzz_yyz_0[i] * pa_x[i] - ta_xzzz_yyz_1[i] * pc_x[i];

        ta_xxzzz_yzz_0[i] = ta_zzz_yzz_0[i] * fe_0 - ta_zzz_yzz_1[i] * fe_0 + ta_xzzz_yzz_0[i] * pa_x[i] - ta_xzzz_yzz_1[i] * pc_x[i];

        ta_xxzzz_zzz_0[i] = ta_zzz_zzz_0[i] * fe_0 - ta_zzz_zzz_1[i] * fe_0 + ta_xzzz_zzz_0[i] * pa_x[i] - ta_xzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 100-110 components of targeted buffer : HF

    auto ta_xyyyy_xxx_0 = pbuffer.data(idx_npot_0_hf + 100);

    auto ta_xyyyy_xxy_0 = pbuffer.data(idx_npot_0_hf + 101);

    auto ta_xyyyy_xxz_0 = pbuffer.data(idx_npot_0_hf + 102);

    auto ta_xyyyy_xyy_0 = pbuffer.data(idx_npot_0_hf + 103);

    auto ta_xyyyy_xyz_0 = pbuffer.data(idx_npot_0_hf + 104);

    auto ta_xyyyy_xzz_0 = pbuffer.data(idx_npot_0_hf + 105);

    auto ta_xyyyy_yyy_0 = pbuffer.data(idx_npot_0_hf + 106);

    auto ta_xyyyy_yyz_0 = pbuffer.data(idx_npot_0_hf + 107);

    auto ta_xyyyy_yzz_0 = pbuffer.data(idx_npot_0_hf + 108);

    auto ta_xyyyy_zzz_0 = pbuffer.data(idx_npot_0_hf + 109);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta_xyyyy_xxx_0, \
                             ta_xyyyy_xxy_0, \
                             ta_xyyyy_xxz_0, \
                             ta_xyyyy_xyy_0, \
                             ta_xyyyy_xyz_0, \
                             ta_xyyyy_xzz_0, \
                             ta_xyyyy_yyy_0, \
                             ta_xyyyy_yyz_0, \
                             ta_xyyyy_yzz_0, \
                             ta_xyyyy_zzz_0, \
                             ta_yyyy_xx_0,   \
                             ta_yyyy_xx_1,   \
                             ta_yyyy_xxx_0,  \
                             ta_yyyy_xxx_1,  \
                             ta_yyyy_xxy_0,  \
                             ta_yyyy_xxy_1,  \
                             ta_yyyy_xxz_0,  \
                             ta_yyyy_xxz_1,  \
                             ta_yyyy_xy_0,   \
                             ta_yyyy_xy_1,   \
                             ta_yyyy_xyy_0,  \
                             ta_yyyy_xyy_1,  \
                             ta_yyyy_xyz_0,  \
                             ta_yyyy_xyz_1,  \
                             ta_yyyy_xz_0,   \
                             ta_yyyy_xz_1,   \
                             ta_yyyy_xzz_0,  \
                             ta_yyyy_xzz_1,  \
                             ta_yyyy_yy_0,   \
                             ta_yyyy_yy_1,   \
                             ta_yyyy_yyy_0,  \
                             ta_yyyy_yyy_1,  \
                             ta_yyyy_yyz_0,  \
                             ta_yyyy_yyz_1,  \
                             ta_yyyy_yz_0,   \
                             ta_yyyy_yz_1,   \
                             ta_yyyy_yzz_0,  \
                             ta_yyyy_yzz_1,  \
                             ta_yyyy_zz_0,   \
                             ta_yyyy_zz_1,   \
                             ta_yyyy_zzz_0,  \
                             ta_yyyy_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyy_xxx_0[i] = 3.0 * ta_yyyy_xx_0[i] * fe_0 - 3.0 * ta_yyyy_xx_1[i] * fe_0 + ta_yyyy_xxx_0[i] * pa_x[i] - ta_yyyy_xxx_1[i] * pc_x[i];

        ta_xyyyy_xxy_0[i] = 2.0 * ta_yyyy_xy_0[i] * fe_0 - 2.0 * ta_yyyy_xy_1[i] * fe_0 + ta_yyyy_xxy_0[i] * pa_x[i] - ta_yyyy_xxy_1[i] * pc_x[i];

        ta_xyyyy_xxz_0[i] = 2.0 * ta_yyyy_xz_0[i] * fe_0 - 2.0 * ta_yyyy_xz_1[i] * fe_0 + ta_yyyy_xxz_0[i] * pa_x[i] - ta_yyyy_xxz_1[i] * pc_x[i];

        ta_xyyyy_xyy_0[i] = ta_yyyy_yy_0[i] * fe_0 - ta_yyyy_yy_1[i] * fe_0 + ta_yyyy_xyy_0[i] * pa_x[i] - ta_yyyy_xyy_1[i] * pc_x[i];

        ta_xyyyy_xyz_0[i] = ta_yyyy_yz_0[i] * fe_0 - ta_yyyy_yz_1[i] * fe_0 + ta_yyyy_xyz_0[i] * pa_x[i] - ta_yyyy_xyz_1[i] * pc_x[i];

        ta_xyyyy_xzz_0[i] = ta_yyyy_zz_0[i] * fe_0 - ta_yyyy_zz_1[i] * fe_0 + ta_yyyy_xzz_0[i] * pa_x[i] - ta_yyyy_xzz_1[i] * pc_x[i];

        ta_xyyyy_yyy_0[i] = ta_yyyy_yyy_0[i] * pa_x[i] - ta_yyyy_yyy_1[i] * pc_x[i];

        ta_xyyyy_yyz_0[i] = ta_yyyy_yyz_0[i] * pa_x[i] - ta_yyyy_yyz_1[i] * pc_x[i];

        ta_xyyyy_yzz_0[i] = ta_yyyy_yzz_0[i] * pa_x[i] - ta_yyyy_yzz_1[i] * pc_x[i];

        ta_xyyyy_zzz_0[i] = ta_yyyy_zzz_0[i] * pa_x[i] - ta_yyyy_zzz_1[i] * pc_x[i];
    }

    // Set up 110-120 components of targeted buffer : HF

    auto ta_xyyyz_xxx_0 = pbuffer.data(idx_npot_0_hf + 110);

    auto ta_xyyyz_xxy_0 = pbuffer.data(idx_npot_0_hf + 111);

    auto ta_xyyyz_xxz_0 = pbuffer.data(idx_npot_0_hf + 112);

    auto ta_xyyyz_xyy_0 = pbuffer.data(idx_npot_0_hf + 113);

    auto ta_xyyyz_xyz_0 = pbuffer.data(idx_npot_0_hf + 114);

    auto ta_xyyyz_xzz_0 = pbuffer.data(idx_npot_0_hf + 115);

    auto ta_xyyyz_yyy_0 = pbuffer.data(idx_npot_0_hf + 116);

    auto ta_xyyyz_yyz_0 = pbuffer.data(idx_npot_0_hf + 117);

    auto ta_xyyyz_yzz_0 = pbuffer.data(idx_npot_0_hf + 118);

    auto ta_xyyyz_zzz_0 = pbuffer.data(idx_npot_0_hf + 119);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta_xyyy_xxx_0,  \
                             ta_xyyy_xxx_1,  \
                             ta_xyyy_xxy_0,  \
                             ta_xyyy_xxy_1,  \
                             ta_xyyy_xyy_0,  \
                             ta_xyyy_xyy_1,  \
                             ta_xyyyz_xxx_0, \
                             ta_xyyyz_xxy_0, \
                             ta_xyyyz_xxz_0, \
                             ta_xyyyz_xyy_0, \
                             ta_xyyyz_xyz_0, \
                             ta_xyyyz_xzz_0, \
                             ta_xyyyz_yyy_0, \
                             ta_xyyyz_yyz_0, \
                             ta_xyyyz_yzz_0, \
                             ta_xyyyz_zzz_0, \
                             ta_yyyz_xxz_0,  \
                             ta_yyyz_xxz_1,  \
                             ta_yyyz_xyz_0,  \
                             ta_yyyz_xyz_1,  \
                             ta_yyyz_xz_0,   \
                             ta_yyyz_xz_1,   \
                             ta_yyyz_xzz_0,  \
                             ta_yyyz_xzz_1,  \
                             ta_yyyz_yyy_0,  \
                             ta_yyyz_yyy_1,  \
                             ta_yyyz_yyz_0,  \
                             ta_yyyz_yyz_1,  \
                             ta_yyyz_yz_0,   \
                             ta_yyyz_yz_1,   \
                             ta_yyyz_yzz_0,  \
                             ta_yyyz_yzz_1,  \
                             ta_yyyz_zz_0,   \
                             ta_yyyz_zz_1,   \
                             ta_yyyz_zzz_0,  \
                             ta_yyyz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyz_xxx_0[i] = ta_xyyy_xxx_0[i] * pa_z[i] - ta_xyyy_xxx_1[i] * pc_z[i];

        ta_xyyyz_xxy_0[i] = ta_xyyy_xxy_0[i] * pa_z[i] - ta_xyyy_xxy_1[i] * pc_z[i];

        ta_xyyyz_xxz_0[i] = 2.0 * ta_yyyz_xz_0[i] * fe_0 - 2.0 * ta_yyyz_xz_1[i] * fe_0 + ta_yyyz_xxz_0[i] * pa_x[i] - ta_yyyz_xxz_1[i] * pc_x[i];

        ta_xyyyz_xyy_0[i] = ta_xyyy_xyy_0[i] * pa_z[i] - ta_xyyy_xyy_1[i] * pc_z[i];

        ta_xyyyz_xyz_0[i] = ta_yyyz_yz_0[i] * fe_0 - ta_yyyz_yz_1[i] * fe_0 + ta_yyyz_xyz_0[i] * pa_x[i] - ta_yyyz_xyz_1[i] * pc_x[i];

        ta_xyyyz_xzz_0[i] = ta_yyyz_zz_0[i] * fe_0 - ta_yyyz_zz_1[i] * fe_0 + ta_yyyz_xzz_0[i] * pa_x[i] - ta_yyyz_xzz_1[i] * pc_x[i];

        ta_xyyyz_yyy_0[i] = ta_yyyz_yyy_0[i] * pa_x[i] - ta_yyyz_yyy_1[i] * pc_x[i];

        ta_xyyyz_yyz_0[i] = ta_yyyz_yyz_0[i] * pa_x[i] - ta_yyyz_yyz_1[i] * pc_x[i];

        ta_xyyyz_yzz_0[i] = ta_yyyz_yzz_0[i] * pa_x[i] - ta_yyyz_yzz_1[i] * pc_x[i];

        ta_xyyyz_zzz_0[i] = ta_yyyz_zzz_0[i] * pa_x[i] - ta_yyyz_zzz_1[i] * pc_x[i];
    }

    // Set up 120-130 components of targeted buffer : HF

    auto ta_xyyzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 120);

    auto ta_xyyzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 121);

    auto ta_xyyzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 122);

    auto ta_xyyzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 123);

    auto ta_xyyzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 124);

    auto ta_xyyzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 125);

    auto ta_xyyzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 126);

    auto ta_xyyzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 127);

    auto ta_xyyzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 128);

    auto ta_xyyzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 129);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta_xyyzz_xxx_0, \
                             ta_xyyzz_xxy_0, \
                             ta_xyyzz_xxz_0, \
                             ta_xyyzz_xyy_0, \
                             ta_xyyzz_xyz_0, \
                             ta_xyyzz_xzz_0, \
                             ta_xyyzz_yyy_0, \
                             ta_xyyzz_yyz_0, \
                             ta_xyyzz_yzz_0, \
                             ta_xyyzz_zzz_0, \
                             ta_yyzz_xx_0,   \
                             ta_yyzz_xx_1,   \
                             ta_yyzz_xxx_0,  \
                             ta_yyzz_xxx_1,  \
                             ta_yyzz_xxy_0,  \
                             ta_yyzz_xxy_1,  \
                             ta_yyzz_xxz_0,  \
                             ta_yyzz_xxz_1,  \
                             ta_yyzz_xy_0,   \
                             ta_yyzz_xy_1,   \
                             ta_yyzz_xyy_0,  \
                             ta_yyzz_xyy_1,  \
                             ta_yyzz_xyz_0,  \
                             ta_yyzz_xyz_1,  \
                             ta_yyzz_xz_0,   \
                             ta_yyzz_xz_1,   \
                             ta_yyzz_xzz_0,  \
                             ta_yyzz_xzz_1,  \
                             ta_yyzz_yy_0,   \
                             ta_yyzz_yy_1,   \
                             ta_yyzz_yyy_0,  \
                             ta_yyzz_yyy_1,  \
                             ta_yyzz_yyz_0,  \
                             ta_yyzz_yyz_1,  \
                             ta_yyzz_yz_0,   \
                             ta_yyzz_yz_1,   \
                             ta_yyzz_yzz_0,  \
                             ta_yyzz_yzz_1,  \
                             ta_yyzz_zz_0,   \
                             ta_yyzz_zz_1,   \
                             ta_yyzz_zzz_0,  \
                             ta_yyzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyzz_xxx_0[i] = 3.0 * ta_yyzz_xx_0[i] * fe_0 - 3.0 * ta_yyzz_xx_1[i] * fe_0 + ta_yyzz_xxx_0[i] * pa_x[i] - ta_yyzz_xxx_1[i] * pc_x[i];

        ta_xyyzz_xxy_0[i] = 2.0 * ta_yyzz_xy_0[i] * fe_0 - 2.0 * ta_yyzz_xy_1[i] * fe_0 + ta_yyzz_xxy_0[i] * pa_x[i] - ta_yyzz_xxy_1[i] * pc_x[i];

        ta_xyyzz_xxz_0[i] = 2.0 * ta_yyzz_xz_0[i] * fe_0 - 2.0 * ta_yyzz_xz_1[i] * fe_0 + ta_yyzz_xxz_0[i] * pa_x[i] - ta_yyzz_xxz_1[i] * pc_x[i];

        ta_xyyzz_xyy_0[i] = ta_yyzz_yy_0[i] * fe_0 - ta_yyzz_yy_1[i] * fe_0 + ta_yyzz_xyy_0[i] * pa_x[i] - ta_yyzz_xyy_1[i] * pc_x[i];

        ta_xyyzz_xyz_0[i] = ta_yyzz_yz_0[i] * fe_0 - ta_yyzz_yz_1[i] * fe_0 + ta_yyzz_xyz_0[i] * pa_x[i] - ta_yyzz_xyz_1[i] * pc_x[i];

        ta_xyyzz_xzz_0[i] = ta_yyzz_zz_0[i] * fe_0 - ta_yyzz_zz_1[i] * fe_0 + ta_yyzz_xzz_0[i] * pa_x[i] - ta_yyzz_xzz_1[i] * pc_x[i];

        ta_xyyzz_yyy_0[i] = ta_yyzz_yyy_0[i] * pa_x[i] - ta_yyzz_yyy_1[i] * pc_x[i];

        ta_xyyzz_yyz_0[i] = ta_yyzz_yyz_0[i] * pa_x[i] - ta_yyzz_yyz_1[i] * pc_x[i];

        ta_xyyzz_yzz_0[i] = ta_yyzz_yzz_0[i] * pa_x[i] - ta_yyzz_yzz_1[i] * pc_x[i];

        ta_xyyzz_zzz_0[i] = ta_yyzz_zzz_0[i] * pa_x[i] - ta_yyzz_zzz_1[i] * pc_x[i];
    }

    // Set up 130-140 components of targeted buffer : HF

    auto ta_xyzzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 130);

    auto ta_xyzzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 131);

    auto ta_xyzzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 132);

    auto ta_xyzzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 133);

    auto ta_xyzzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 134);

    auto ta_xyzzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 135);

    auto ta_xyzzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 136);

    auto ta_xyzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 137);

    auto ta_xyzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 138);

    auto ta_xyzzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 139);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta_xyzzz_xxx_0, \
                             ta_xyzzz_xxy_0, \
                             ta_xyzzz_xxz_0, \
                             ta_xyzzz_xyy_0, \
                             ta_xyzzz_xyz_0, \
                             ta_xyzzz_xzz_0, \
                             ta_xyzzz_yyy_0, \
                             ta_xyzzz_yyz_0, \
                             ta_xyzzz_yzz_0, \
                             ta_xyzzz_zzz_0, \
                             ta_xzzz_xxx_0,  \
                             ta_xzzz_xxx_1,  \
                             ta_xzzz_xxz_0,  \
                             ta_xzzz_xxz_1,  \
                             ta_xzzz_xzz_0,  \
                             ta_xzzz_xzz_1,  \
                             ta_yzzz_xxy_0,  \
                             ta_yzzz_xxy_1,  \
                             ta_yzzz_xy_0,   \
                             ta_yzzz_xy_1,   \
                             ta_yzzz_xyy_0,  \
                             ta_yzzz_xyy_1,  \
                             ta_yzzz_xyz_0,  \
                             ta_yzzz_xyz_1,  \
                             ta_yzzz_yy_0,   \
                             ta_yzzz_yy_1,   \
                             ta_yzzz_yyy_0,  \
                             ta_yzzz_yyy_1,  \
                             ta_yzzz_yyz_0,  \
                             ta_yzzz_yyz_1,  \
                             ta_yzzz_yz_0,   \
                             ta_yzzz_yz_1,   \
                             ta_yzzz_yzz_0,  \
                             ta_yzzz_yzz_1,  \
                             ta_yzzz_zzz_0,  \
                             ta_yzzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyzzz_xxx_0[i] = ta_xzzz_xxx_0[i] * pa_y[i] - ta_xzzz_xxx_1[i] * pc_y[i];

        ta_xyzzz_xxy_0[i] = 2.0 * ta_yzzz_xy_0[i] * fe_0 - 2.0 * ta_yzzz_xy_1[i] * fe_0 + ta_yzzz_xxy_0[i] * pa_x[i] - ta_yzzz_xxy_1[i] * pc_x[i];

        ta_xyzzz_xxz_0[i] = ta_xzzz_xxz_0[i] * pa_y[i] - ta_xzzz_xxz_1[i] * pc_y[i];

        ta_xyzzz_xyy_0[i] = ta_yzzz_yy_0[i] * fe_0 - ta_yzzz_yy_1[i] * fe_0 + ta_yzzz_xyy_0[i] * pa_x[i] - ta_yzzz_xyy_1[i] * pc_x[i];

        ta_xyzzz_xyz_0[i] = ta_yzzz_yz_0[i] * fe_0 - ta_yzzz_yz_1[i] * fe_0 + ta_yzzz_xyz_0[i] * pa_x[i] - ta_yzzz_xyz_1[i] * pc_x[i];

        ta_xyzzz_xzz_0[i] = ta_xzzz_xzz_0[i] * pa_y[i] - ta_xzzz_xzz_1[i] * pc_y[i];

        ta_xyzzz_yyy_0[i] = ta_yzzz_yyy_0[i] * pa_x[i] - ta_yzzz_yyy_1[i] * pc_x[i];

        ta_xyzzz_yyz_0[i] = ta_yzzz_yyz_0[i] * pa_x[i] - ta_yzzz_yyz_1[i] * pc_x[i];

        ta_xyzzz_yzz_0[i] = ta_yzzz_yzz_0[i] * pa_x[i] - ta_yzzz_yzz_1[i] * pc_x[i];

        ta_xyzzz_zzz_0[i] = ta_yzzz_zzz_0[i] * pa_x[i] - ta_yzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 140-150 components of targeted buffer : HF

    auto ta_xzzzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 140);

    auto ta_xzzzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 141);

    auto ta_xzzzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 142);

    auto ta_xzzzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 143);

    auto ta_xzzzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 144);

    auto ta_xzzzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 145);

    auto ta_xzzzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 146);

    auto ta_xzzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 147);

    auto ta_xzzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 148);

    auto ta_xzzzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 149);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta_xzzzz_xxx_0, \
                             ta_xzzzz_xxy_0, \
                             ta_xzzzz_xxz_0, \
                             ta_xzzzz_xyy_0, \
                             ta_xzzzz_xyz_0, \
                             ta_xzzzz_xzz_0, \
                             ta_xzzzz_yyy_0, \
                             ta_xzzzz_yyz_0, \
                             ta_xzzzz_yzz_0, \
                             ta_xzzzz_zzz_0, \
                             ta_zzzz_xx_0,   \
                             ta_zzzz_xx_1,   \
                             ta_zzzz_xxx_0,  \
                             ta_zzzz_xxx_1,  \
                             ta_zzzz_xxy_0,  \
                             ta_zzzz_xxy_1,  \
                             ta_zzzz_xxz_0,  \
                             ta_zzzz_xxz_1,  \
                             ta_zzzz_xy_0,   \
                             ta_zzzz_xy_1,   \
                             ta_zzzz_xyy_0,  \
                             ta_zzzz_xyy_1,  \
                             ta_zzzz_xyz_0,  \
                             ta_zzzz_xyz_1,  \
                             ta_zzzz_xz_0,   \
                             ta_zzzz_xz_1,   \
                             ta_zzzz_xzz_0,  \
                             ta_zzzz_xzz_1,  \
                             ta_zzzz_yy_0,   \
                             ta_zzzz_yy_1,   \
                             ta_zzzz_yyy_0,  \
                             ta_zzzz_yyy_1,  \
                             ta_zzzz_yyz_0,  \
                             ta_zzzz_yyz_1,  \
                             ta_zzzz_yz_0,   \
                             ta_zzzz_yz_1,   \
                             ta_zzzz_yzz_0,  \
                             ta_zzzz_yzz_1,  \
                             ta_zzzz_zz_0,   \
                             ta_zzzz_zz_1,   \
                             ta_zzzz_zzz_0,  \
                             ta_zzzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzzz_xxx_0[i] = 3.0 * ta_zzzz_xx_0[i] * fe_0 - 3.0 * ta_zzzz_xx_1[i] * fe_0 + ta_zzzz_xxx_0[i] * pa_x[i] - ta_zzzz_xxx_1[i] * pc_x[i];

        ta_xzzzz_xxy_0[i] = 2.0 * ta_zzzz_xy_0[i] * fe_0 - 2.0 * ta_zzzz_xy_1[i] * fe_0 + ta_zzzz_xxy_0[i] * pa_x[i] - ta_zzzz_xxy_1[i] * pc_x[i];

        ta_xzzzz_xxz_0[i] = 2.0 * ta_zzzz_xz_0[i] * fe_0 - 2.0 * ta_zzzz_xz_1[i] * fe_0 + ta_zzzz_xxz_0[i] * pa_x[i] - ta_zzzz_xxz_1[i] * pc_x[i];

        ta_xzzzz_xyy_0[i] = ta_zzzz_yy_0[i] * fe_0 - ta_zzzz_yy_1[i] * fe_0 + ta_zzzz_xyy_0[i] * pa_x[i] - ta_zzzz_xyy_1[i] * pc_x[i];

        ta_xzzzz_xyz_0[i] = ta_zzzz_yz_0[i] * fe_0 - ta_zzzz_yz_1[i] * fe_0 + ta_zzzz_xyz_0[i] * pa_x[i] - ta_zzzz_xyz_1[i] * pc_x[i];

        ta_xzzzz_xzz_0[i] = ta_zzzz_zz_0[i] * fe_0 - ta_zzzz_zz_1[i] * fe_0 + ta_zzzz_xzz_0[i] * pa_x[i] - ta_zzzz_xzz_1[i] * pc_x[i];

        ta_xzzzz_yyy_0[i] = ta_zzzz_yyy_0[i] * pa_x[i] - ta_zzzz_yyy_1[i] * pc_x[i];

        ta_xzzzz_yyz_0[i] = ta_zzzz_yyz_0[i] * pa_x[i] - ta_zzzz_yyz_1[i] * pc_x[i];

        ta_xzzzz_yzz_0[i] = ta_zzzz_yzz_0[i] * pa_x[i] - ta_zzzz_yzz_1[i] * pc_x[i];

        ta_xzzzz_zzz_0[i] = ta_zzzz_zzz_0[i] * pa_x[i] - ta_zzzz_zzz_1[i] * pc_x[i];
    }

    // Set up 150-160 components of targeted buffer : HF

    auto ta_yyyyy_xxx_0 = pbuffer.data(idx_npot_0_hf + 150);

    auto ta_yyyyy_xxy_0 = pbuffer.data(idx_npot_0_hf + 151);

    auto ta_yyyyy_xxz_0 = pbuffer.data(idx_npot_0_hf + 152);

    auto ta_yyyyy_xyy_0 = pbuffer.data(idx_npot_0_hf + 153);

    auto ta_yyyyy_xyz_0 = pbuffer.data(idx_npot_0_hf + 154);

    auto ta_yyyyy_xzz_0 = pbuffer.data(idx_npot_0_hf + 155);

    auto ta_yyyyy_yyy_0 = pbuffer.data(idx_npot_0_hf + 156);

    auto ta_yyyyy_yyz_0 = pbuffer.data(idx_npot_0_hf + 157);

    auto ta_yyyyy_yzz_0 = pbuffer.data(idx_npot_0_hf + 158);

    auto ta_yyyyy_zzz_0 = pbuffer.data(idx_npot_0_hf + 159);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta_yyy_xxx_0,   \
                             ta_yyy_xxx_1,   \
                             ta_yyy_xxy_0,   \
                             ta_yyy_xxy_1,   \
                             ta_yyy_xxz_0,   \
                             ta_yyy_xxz_1,   \
                             ta_yyy_xyy_0,   \
                             ta_yyy_xyy_1,   \
                             ta_yyy_xyz_0,   \
                             ta_yyy_xyz_1,   \
                             ta_yyy_xzz_0,   \
                             ta_yyy_xzz_1,   \
                             ta_yyy_yyy_0,   \
                             ta_yyy_yyy_1,   \
                             ta_yyy_yyz_0,   \
                             ta_yyy_yyz_1,   \
                             ta_yyy_yzz_0,   \
                             ta_yyy_yzz_1,   \
                             ta_yyy_zzz_0,   \
                             ta_yyy_zzz_1,   \
                             ta_yyyy_xx_0,   \
                             ta_yyyy_xx_1,   \
                             ta_yyyy_xxx_0,  \
                             ta_yyyy_xxx_1,  \
                             ta_yyyy_xxy_0,  \
                             ta_yyyy_xxy_1,  \
                             ta_yyyy_xxz_0,  \
                             ta_yyyy_xxz_1,  \
                             ta_yyyy_xy_0,   \
                             ta_yyyy_xy_1,   \
                             ta_yyyy_xyy_0,  \
                             ta_yyyy_xyy_1,  \
                             ta_yyyy_xyz_0,  \
                             ta_yyyy_xyz_1,  \
                             ta_yyyy_xz_0,   \
                             ta_yyyy_xz_1,   \
                             ta_yyyy_xzz_0,  \
                             ta_yyyy_xzz_1,  \
                             ta_yyyy_yy_0,   \
                             ta_yyyy_yy_1,   \
                             ta_yyyy_yyy_0,  \
                             ta_yyyy_yyy_1,  \
                             ta_yyyy_yyz_0,  \
                             ta_yyyy_yyz_1,  \
                             ta_yyyy_yz_0,   \
                             ta_yyyy_yz_1,   \
                             ta_yyyy_yzz_0,  \
                             ta_yyyy_yzz_1,  \
                             ta_yyyy_zz_0,   \
                             ta_yyyy_zz_1,   \
                             ta_yyyy_zzz_0,  \
                             ta_yyyy_zzz_1,  \
                             ta_yyyyy_xxx_0, \
                             ta_yyyyy_xxy_0, \
                             ta_yyyyy_xxz_0, \
                             ta_yyyyy_xyy_0, \
                             ta_yyyyy_xyz_0, \
                             ta_yyyyy_xzz_0, \
                             ta_yyyyy_yyy_0, \
                             ta_yyyyy_yyz_0, \
                             ta_yyyyy_yzz_0, \
                             ta_yyyyy_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyy_xxx_0[i] = 4.0 * ta_yyy_xxx_0[i] * fe_0 - 4.0 * ta_yyy_xxx_1[i] * fe_0 + ta_yyyy_xxx_0[i] * pa_y[i] - ta_yyyy_xxx_1[i] * pc_y[i];

        ta_yyyyy_xxy_0[i] = 4.0 * ta_yyy_xxy_0[i] * fe_0 - 4.0 * ta_yyy_xxy_1[i] * fe_0 + ta_yyyy_xx_0[i] * fe_0 - ta_yyyy_xx_1[i] * fe_0 +
                            ta_yyyy_xxy_0[i] * pa_y[i] - ta_yyyy_xxy_1[i] * pc_y[i];

        ta_yyyyy_xxz_0[i] = 4.0 * ta_yyy_xxz_0[i] * fe_0 - 4.0 * ta_yyy_xxz_1[i] * fe_0 + ta_yyyy_xxz_0[i] * pa_y[i] - ta_yyyy_xxz_1[i] * pc_y[i];

        ta_yyyyy_xyy_0[i] = 4.0 * ta_yyy_xyy_0[i] * fe_0 - 4.0 * ta_yyy_xyy_1[i] * fe_0 + 2.0 * ta_yyyy_xy_0[i] * fe_0 -
                            2.0 * ta_yyyy_xy_1[i] * fe_0 + ta_yyyy_xyy_0[i] * pa_y[i] - ta_yyyy_xyy_1[i] * pc_y[i];

        ta_yyyyy_xyz_0[i] = 4.0 * ta_yyy_xyz_0[i] * fe_0 - 4.0 * ta_yyy_xyz_1[i] * fe_0 + ta_yyyy_xz_0[i] * fe_0 - ta_yyyy_xz_1[i] * fe_0 +
                            ta_yyyy_xyz_0[i] * pa_y[i] - ta_yyyy_xyz_1[i] * pc_y[i];

        ta_yyyyy_xzz_0[i] = 4.0 * ta_yyy_xzz_0[i] * fe_0 - 4.0 * ta_yyy_xzz_1[i] * fe_0 + ta_yyyy_xzz_0[i] * pa_y[i] - ta_yyyy_xzz_1[i] * pc_y[i];

        ta_yyyyy_yyy_0[i] = 4.0 * ta_yyy_yyy_0[i] * fe_0 - 4.0 * ta_yyy_yyy_1[i] * fe_0 + 3.0 * ta_yyyy_yy_0[i] * fe_0 -
                            3.0 * ta_yyyy_yy_1[i] * fe_0 + ta_yyyy_yyy_0[i] * pa_y[i] - ta_yyyy_yyy_1[i] * pc_y[i];

        ta_yyyyy_yyz_0[i] = 4.0 * ta_yyy_yyz_0[i] * fe_0 - 4.0 * ta_yyy_yyz_1[i] * fe_0 + 2.0 * ta_yyyy_yz_0[i] * fe_0 -
                            2.0 * ta_yyyy_yz_1[i] * fe_0 + ta_yyyy_yyz_0[i] * pa_y[i] - ta_yyyy_yyz_1[i] * pc_y[i];

        ta_yyyyy_yzz_0[i] = 4.0 * ta_yyy_yzz_0[i] * fe_0 - 4.0 * ta_yyy_yzz_1[i] * fe_0 + ta_yyyy_zz_0[i] * fe_0 - ta_yyyy_zz_1[i] * fe_0 +
                            ta_yyyy_yzz_0[i] * pa_y[i] - ta_yyyy_yzz_1[i] * pc_y[i];

        ta_yyyyy_zzz_0[i] = 4.0 * ta_yyy_zzz_0[i] * fe_0 - 4.0 * ta_yyy_zzz_1[i] * fe_0 + ta_yyyy_zzz_0[i] * pa_y[i] - ta_yyyy_zzz_1[i] * pc_y[i];
    }

    // Set up 160-170 components of targeted buffer : HF

    auto ta_yyyyz_xxx_0 = pbuffer.data(idx_npot_0_hf + 160);

    auto ta_yyyyz_xxy_0 = pbuffer.data(idx_npot_0_hf + 161);

    auto ta_yyyyz_xxz_0 = pbuffer.data(idx_npot_0_hf + 162);

    auto ta_yyyyz_xyy_0 = pbuffer.data(idx_npot_0_hf + 163);

    auto ta_yyyyz_xyz_0 = pbuffer.data(idx_npot_0_hf + 164);

    auto ta_yyyyz_xzz_0 = pbuffer.data(idx_npot_0_hf + 165);

    auto ta_yyyyz_yyy_0 = pbuffer.data(idx_npot_0_hf + 166);

    auto ta_yyyyz_yyz_0 = pbuffer.data(idx_npot_0_hf + 167);

    auto ta_yyyyz_yzz_0 = pbuffer.data(idx_npot_0_hf + 168);

    auto ta_yyyyz_zzz_0 = pbuffer.data(idx_npot_0_hf + 169);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta_yyyy_xxx_0,  \
                             ta_yyyy_xxx_1,  \
                             ta_yyyy_xxy_0,  \
                             ta_yyyy_xxy_1,  \
                             ta_yyyy_xy_0,   \
                             ta_yyyy_xy_1,   \
                             ta_yyyy_xyy_0,  \
                             ta_yyyy_xyy_1,  \
                             ta_yyyy_xyz_0,  \
                             ta_yyyy_xyz_1,  \
                             ta_yyyy_yy_0,   \
                             ta_yyyy_yy_1,   \
                             ta_yyyy_yyy_0,  \
                             ta_yyyy_yyy_1,  \
                             ta_yyyy_yyz_0,  \
                             ta_yyyy_yyz_1,  \
                             ta_yyyy_yz_0,   \
                             ta_yyyy_yz_1,   \
                             ta_yyyy_yzz_0,  \
                             ta_yyyy_yzz_1,  \
                             ta_yyyyz_xxx_0, \
                             ta_yyyyz_xxy_0, \
                             ta_yyyyz_xxz_0, \
                             ta_yyyyz_xyy_0, \
                             ta_yyyyz_xyz_0, \
                             ta_yyyyz_xzz_0, \
                             ta_yyyyz_yyy_0, \
                             ta_yyyyz_yyz_0, \
                             ta_yyyyz_yzz_0, \
                             ta_yyyyz_zzz_0, \
                             ta_yyyz_xxz_0,  \
                             ta_yyyz_xxz_1,  \
                             ta_yyyz_xzz_0,  \
                             ta_yyyz_xzz_1,  \
                             ta_yyyz_zzz_0,  \
                             ta_yyyz_zzz_1,  \
                             ta_yyz_xxz_0,   \
                             ta_yyz_xxz_1,   \
                             ta_yyz_xzz_0,   \
                             ta_yyz_xzz_1,   \
                             ta_yyz_zzz_0,   \
                             ta_yyz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyz_xxx_0[i] = ta_yyyy_xxx_0[i] * pa_z[i] - ta_yyyy_xxx_1[i] * pc_z[i];

        ta_yyyyz_xxy_0[i] = ta_yyyy_xxy_0[i] * pa_z[i] - ta_yyyy_xxy_1[i] * pc_z[i];

        ta_yyyyz_xxz_0[i] = 3.0 * ta_yyz_xxz_0[i] * fe_0 - 3.0 * ta_yyz_xxz_1[i] * fe_0 + ta_yyyz_xxz_0[i] * pa_y[i] - ta_yyyz_xxz_1[i] * pc_y[i];

        ta_yyyyz_xyy_0[i] = ta_yyyy_xyy_0[i] * pa_z[i] - ta_yyyy_xyy_1[i] * pc_z[i];

        ta_yyyyz_xyz_0[i] = ta_yyyy_xy_0[i] * fe_0 - ta_yyyy_xy_1[i] * fe_0 + ta_yyyy_xyz_0[i] * pa_z[i] - ta_yyyy_xyz_1[i] * pc_z[i];

        ta_yyyyz_xzz_0[i] = 3.0 * ta_yyz_xzz_0[i] * fe_0 - 3.0 * ta_yyz_xzz_1[i] * fe_0 + ta_yyyz_xzz_0[i] * pa_y[i] - ta_yyyz_xzz_1[i] * pc_y[i];

        ta_yyyyz_yyy_0[i] = ta_yyyy_yyy_0[i] * pa_z[i] - ta_yyyy_yyy_1[i] * pc_z[i];

        ta_yyyyz_yyz_0[i] = ta_yyyy_yy_0[i] * fe_0 - ta_yyyy_yy_1[i] * fe_0 + ta_yyyy_yyz_0[i] * pa_z[i] - ta_yyyy_yyz_1[i] * pc_z[i];

        ta_yyyyz_yzz_0[i] = 2.0 * ta_yyyy_yz_0[i] * fe_0 - 2.0 * ta_yyyy_yz_1[i] * fe_0 + ta_yyyy_yzz_0[i] * pa_z[i] - ta_yyyy_yzz_1[i] * pc_z[i];

        ta_yyyyz_zzz_0[i] = 3.0 * ta_yyz_zzz_0[i] * fe_0 - 3.0 * ta_yyz_zzz_1[i] * fe_0 + ta_yyyz_zzz_0[i] * pa_y[i] - ta_yyyz_zzz_1[i] * pc_y[i];
    }

    // Set up 170-180 components of targeted buffer : HF

    auto ta_yyyzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 170);

    auto ta_yyyzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 171);

    auto ta_yyyzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 172);

    auto ta_yyyzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 173);

    auto ta_yyyzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 174);

    auto ta_yyyzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 175);

    auto ta_yyyzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 176);

    auto ta_yyyzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 177);

    auto ta_yyyzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 178);

    auto ta_yyyzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 179);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta_yyy_xxy_0,   \
                             ta_yyy_xxy_1,   \
                             ta_yyy_xyy_0,   \
                             ta_yyy_xyy_1,   \
                             ta_yyy_yyy_0,   \
                             ta_yyy_yyy_1,   \
                             ta_yyyz_xxy_0,  \
                             ta_yyyz_xxy_1,  \
                             ta_yyyz_xyy_0,  \
                             ta_yyyz_xyy_1,  \
                             ta_yyyz_yyy_0,  \
                             ta_yyyz_yyy_1,  \
                             ta_yyyzz_xxx_0, \
                             ta_yyyzz_xxy_0, \
                             ta_yyyzz_xxz_0, \
                             ta_yyyzz_xyy_0, \
                             ta_yyyzz_xyz_0, \
                             ta_yyyzz_xzz_0, \
                             ta_yyyzz_yyy_0, \
                             ta_yyyzz_yyz_0, \
                             ta_yyyzz_yzz_0, \
                             ta_yyyzz_zzz_0, \
                             ta_yyzz_xxx_0,  \
                             ta_yyzz_xxx_1,  \
                             ta_yyzz_xxz_0,  \
                             ta_yyzz_xxz_1,  \
                             ta_yyzz_xyz_0,  \
                             ta_yyzz_xyz_1,  \
                             ta_yyzz_xz_0,   \
                             ta_yyzz_xz_1,   \
                             ta_yyzz_xzz_0,  \
                             ta_yyzz_xzz_1,  \
                             ta_yyzz_yyz_0,  \
                             ta_yyzz_yyz_1,  \
                             ta_yyzz_yz_0,   \
                             ta_yyzz_yz_1,   \
                             ta_yyzz_yzz_0,  \
                             ta_yyzz_yzz_1,  \
                             ta_yyzz_zz_0,   \
                             ta_yyzz_zz_1,   \
                             ta_yyzz_zzz_0,  \
                             ta_yyzz_zzz_1,  \
                             ta_yzz_xxx_0,   \
                             ta_yzz_xxx_1,   \
                             ta_yzz_xxz_0,   \
                             ta_yzz_xxz_1,   \
                             ta_yzz_xyz_0,   \
                             ta_yzz_xyz_1,   \
                             ta_yzz_xzz_0,   \
                             ta_yzz_xzz_1,   \
                             ta_yzz_yyz_0,   \
                             ta_yzz_yyz_1,   \
                             ta_yzz_yzz_0,   \
                             ta_yzz_yzz_1,   \
                             ta_yzz_zzz_0,   \
                             ta_yzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyzz_xxx_0[i] = 2.0 * ta_yzz_xxx_0[i] * fe_0 - 2.0 * ta_yzz_xxx_1[i] * fe_0 + ta_yyzz_xxx_0[i] * pa_y[i] - ta_yyzz_xxx_1[i] * pc_y[i];

        ta_yyyzz_xxy_0[i] = ta_yyy_xxy_0[i] * fe_0 - ta_yyy_xxy_1[i] * fe_0 + ta_yyyz_xxy_0[i] * pa_z[i] - ta_yyyz_xxy_1[i] * pc_z[i];

        ta_yyyzz_xxz_0[i] = 2.0 * ta_yzz_xxz_0[i] * fe_0 - 2.0 * ta_yzz_xxz_1[i] * fe_0 + ta_yyzz_xxz_0[i] * pa_y[i] - ta_yyzz_xxz_1[i] * pc_y[i];

        ta_yyyzz_xyy_0[i] = ta_yyy_xyy_0[i] * fe_0 - ta_yyy_xyy_1[i] * fe_0 + ta_yyyz_xyy_0[i] * pa_z[i] - ta_yyyz_xyy_1[i] * pc_z[i];

        ta_yyyzz_xyz_0[i] = 2.0 * ta_yzz_xyz_0[i] * fe_0 - 2.0 * ta_yzz_xyz_1[i] * fe_0 + ta_yyzz_xz_0[i] * fe_0 - ta_yyzz_xz_1[i] * fe_0 +
                            ta_yyzz_xyz_0[i] * pa_y[i] - ta_yyzz_xyz_1[i] * pc_y[i];

        ta_yyyzz_xzz_0[i] = 2.0 * ta_yzz_xzz_0[i] * fe_0 - 2.0 * ta_yzz_xzz_1[i] * fe_0 + ta_yyzz_xzz_0[i] * pa_y[i] - ta_yyzz_xzz_1[i] * pc_y[i];

        ta_yyyzz_yyy_0[i] = ta_yyy_yyy_0[i] * fe_0 - ta_yyy_yyy_1[i] * fe_0 + ta_yyyz_yyy_0[i] * pa_z[i] - ta_yyyz_yyy_1[i] * pc_z[i];

        ta_yyyzz_yyz_0[i] = 2.0 * ta_yzz_yyz_0[i] * fe_0 - 2.0 * ta_yzz_yyz_1[i] * fe_0 + 2.0 * ta_yyzz_yz_0[i] * fe_0 -
                            2.0 * ta_yyzz_yz_1[i] * fe_0 + ta_yyzz_yyz_0[i] * pa_y[i] - ta_yyzz_yyz_1[i] * pc_y[i];

        ta_yyyzz_yzz_0[i] = 2.0 * ta_yzz_yzz_0[i] * fe_0 - 2.0 * ta_yzz_yzz_1[i] * fe_0 + ta_yyzz_zz_0[i] * fe_0 - ta_yyzz_zz_1[i] * fe_0 +
                            ta_yyzz_yzz_0[i] * pa_y[i] - ta_yyzz_yzz_1[i] * pc_y[i];

        ta_yyyzz_zzz_0[i] = 2.0 * ta_yzz_zzz_0[i] * fe_0 - 2.0 * ta_yzz_zzz_1[i] * fe_0 + ta_yyzz_zzz_0[i] * pa_y[i] - ta_yyzz_zzz_1[i] * pc_y[i];
    }

    // Set up 180-190 components of targeted buffer : HF

    auto ta_yyzzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 180);

    auto ta_yyzzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 181);

    auto ta_yyzzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 182);

    auto ta_yyzzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 183);

    auto ta_yyzzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 184);

    auto ta_yyzzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 185);

    auto ta_yyzzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 186);

    auto ta_yyzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 187);

    auto ta_yyzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 188);

    auto ta_yyzzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 189);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta_yyz_xxy_0,   \
                             ta_yyz_xxy_1,   \
                             ta_yyz_xyy_0,   \
                             ta_yyz_xyy_1,   \
                             ta_yyz_yyy_0,   \
                             ta_yyz_yyy_1,   \
                             ta_yyzz_xxy_0,  \
                             ta_yyzz_xxy_1,  \
                             ta_yyzz_xyy_0,  \
                             ta_yyzz_xyy_1,  \
                             ta_yyzz_yyy_0,  \
                             ta_yyzz_yyy_1,  \
                             ta_yyzzz_xxx_0, \
                             ta_yyzzz_xxy_0, \
                             ta_yyzzz_xxz_0, \
                             ta_yyzzz_xyy_0, \
                             ta_yyzzz_xyz_0, \
                             ta_yyzzz_xzz_0, \
                             ta_yyzzz_yyy_0, \
                             ta_yyzzz_yyz_0, \
                             ta_yyzzz_yzz_0, \
                             ta_yyzzz_zzz_0, \
                             ta_yzzz_xxx_0,  \
                             ta_yzzz_xxx_1,  \
                             ta_yzzz_xxz_0,  \
                             ta_yzzz_xxz_1,  \
                             ta_yzzz_xyz_0,  \
                             ta_yzzz_xyz_1,  \
                             ta_yzzz_xz_0,   \
                             ta_yzzz_xz_1,   \
                             ta_yzzz_xzz_0,  \
                             ta_yzzz_xzz_1,  \
                             ta_yzzz_yyz_0,  \
                             ta_yzzz_yyz_1,  \
                             ta_yzzz_yz_0,   \
                             ta_yzzz_yz_1,   \
                             ta_yzzz_yzz_0,  \
                             ta_yzzz_yzz_1,  \
                             ta_yzzz_zz_0,   \
                             ta_yzzz_zz_1,   \
                             ta_yzzz_zzz_0,  \
                             ta_yzzz_zzz_1,  \
                             ta_zzz_xxx_0,   \
                             ta_zzz_xxx_1,   \
                             ta_zzz_xxz_0,   \
                             ta_zzz_xxz_1,   \
                             ta_zzz_xyz_0,   \
                             ta_zzz_xyz_1,   \
                             ta_zzz_xzz_0,   \
                             ta_zzz_xzz_1,   \
                             ta_zzz_yyz_0,   \
                             ta_zzz_yyz_1,   \
                             ta_zzz_yzz_0,   \
                             ta_zzz_yzz_1,   \
                             ta_zzz_zzz_0,   \
                             ta_zzz_zzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzzz_xxx_0[i] = ta_zzz_xxx_0[i] * fe_0 - ta_zzz_xxx_1[i] * fe_0 + ta_yzzz_xxx_0[i] * pa_y[i] - ta_yzzz_xxx_1[i] * pc_y[i];

        ta_yyzzz_xxy_0[i] = 2.0 * ta_yyz_xxy_0[i] * fe_0 - 2.0 * ta_yyz_xxy_1[i] * fe_0 + ta_yyzz_xxy_0[i] * pa_z[i] - ta_yyzz_xxy_1[i] * pc_z[i];

        ta_yyzzz_xxz_0[i] = ta_zzz_xxz_0[i] * fe_0 - ta_zzz_xxz_1[i] * fe_0 + ta_yzzz_xxz_0[i] * pa_y[i] - ta_yzzz_xxz_1[i] * pc_y[i];

        ta_yyzzz_xyy_0[i] = 2.0 * ta_yyz_xyy_0[i] * fe_0 - 2.0 * ta_yyz_xyy_1[i] * fe_0 + ta_yyzz_xyy_0[i] * pa_z[i] - ta_yyzz_xyy_1[i] * pc_z[i];

        ta_yyzzz_xyz_0[i] = ta_zzz_xyz_0[i] * fe_0 - ta_zzz_xyz_1[i] * fe_0 + ta_yzzz_xz_0[i] * fe_0 - ta_yzzz_xz_1[i] * fe_0 +
                            ta_yzzz_xyz_0[i] * pa_y[i] - ta_yzzz_xyz_1[i] * pc_y[i];

        ta_yyzzz_xzz_0[i] = ta_zzz_xzz_0[i] * fe_0 - ta_zzz_xzz_1[i] * fe_0 + ta_yzzz_xzz_0[i] * pa_y[i] - ta_yzzz_xzz_1[i] * pc_y[i];

        ta_yyzzz_yyy_0[i] = 2.0 * ta_yyz_yyy_0[i] * fe_0 - 2.0 * ta_yyz_yyy_1[i] * fe_0 + ta_yyzz_yyy_0[i] * pa_z[i] - ta_yyzz_yyy_1[i] * pc_z[i];

        ta_yyzzz_yyz_0[i] = ta_zzz_yyz_0[i] * fe_0 - ta_zzz_yyz_1[i] * fe_0 + 2.0 * ta_yzzz_yz_0[i] * fe_0 - 2.0 * ta_yzzz_yz_1[i] * fe_0 +
                            ta_yzzz_yyz_0[i] * pa_y[i] - ta_yzzz_yyz_1[i] * pc_y[i];

        ta_yyzzz_yzz_0[i] = ta_zzz_yzz_0[i] * fe_0 - ta_zzz_yzz_1[i] * fe_0 + ta_yzzz_zz_0[i] * fe_0 - ta_yzzz_zz_1[i] * fe_0 +
                            ta_yzzz_yzz_0[i] * pa_y[i] - ta_yzzz_yzz_1[i] * pc_y[i];

        ta_yyzzz_zzz_0[i] = ta_zzz_zzz_0[i] * fe_0 - ta_zzz_zzz_1[i] * fe_0 + ta_yzzz_zzz_0[i] * pa_y[i] - ta_yzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 190-200 components of targeted buffer : HF

    auto ta_yzzzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 190);

    auto ta_yzzzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 191);

    auto ta_yzzzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 192);

    auto ta_yzzzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 193);

    auto ta_yzzzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 194);

    auto ta_yzzzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 195);

    auto ta_yzzzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 196);

    auto ta_yzzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 197);

    auto ta_yzzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 198);

    auto ta_yzzzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 199);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta_yzzzz_xxx_0, \
                             ta_yzzzz_xxy_0, \
                             ta_yzzzz_xxz_0, \
                             ta_yzzzz_xyy_0, \
                             ta_yzzzz_xyz_0, \
                             ta_yzzzz_xzz_0, \
                             ta_yzzzz_yyy_0, \
                             ta_yzzzz_yyz_0, \
                             ta_yzzzz_yzz_0, \
                             ta_yzzzz_zzz_0, \
                             ta_zzzz_xx_0,   \
                             ta_zzzz_xx_1,   \
                             ta_zzzz_xxx_0,  \
                             ta_zzzz_xxx_1,  \
                             ta_zzzz_xxy_0,  \
                             ta_zzzz_xxy_1,  \
                             ta_zzzz_xxz_0,  \
                             ta_zzzz_xxz_1,  \
                             ta_zzzz_xy_0,   \
                             ta_zzzz_xy_1,   \
                             ta_zzzz_xyy_0,  \
                             ta_zzzz_xyy_1,  \
                             ta_zzzz_xyz_0,  \
                             ta_zzzz_xyz_1,  \
                             ta_zzzz_xz_0,   \
                             ta_zzzz_xz_1,   \
                             ta_zzzz_xzz_0,  \
                             ta_zzzz_xzz_1,  \
                             ta_zzzz_yy_0,   \
                             ta_zzzz_yy_1,   \
                             ta_zzzz_yyy_0,  \
                             ta_zzzz_yyy_1,  \
                             ta_zzzz_yyz_0,  \
                             ta_zzzz_yyz_1,  \
                             ta_zzzz_yz_0,   \
                             ta_zzzz_yz_1,   \
                             ta_zzzz_yzz_0,  \
                             ta_zzzz_yzz_1,  \
                             ta_zzzz_zz_0,   \
                             ta_zzzz_zz_1,   \
                             ta_zzzz_zzz_0,  \
                             ta_zzzz_zzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzzz_xxx_0[i] = ta_zzzz_xxx_0[i] * pa_y[i] - ta_zzzz_xxx_1[i] * pc_y[i];

        ta_yzzzz_xxy_0[i] = ta_zzzz_xx_0[i] * fe_0 - ta_zzzz_xx_1[i] * fe_0 + ta_zzzz_xxy_0[i] * pa_y[i] - ta_zzzz_xxy_1[i] * pc_y[i];

        ta_yzzzz_xxz_0[i] = ta_zzzz_xxz_0[i] * pa_y[i] - ta_zzzz_xxz_1[i] * pc_y[i];

        ta_yzzzz_xyy_0[i] = 2.0 * ta_zzzz_xy_0[i] * fe_0 - 2.0 * ta_zzzz_xy_1[i] * fe_0 + ta_zzzz_xyy_0[i] * pa_y[i] - ta_zzzz_xyy_1[i] * pc_y[i];

        ta_yzzzz_xyz_0[i] = ta_zzzz_xz_0[i] * fe_0 - ta_zzzz_xz_1[i] * fe_0 + ta_zzzz_xyz_0[i] * pa_y[i] - ta_zzzz_xyz_1[i] * pc_y[i];

        ta_yzzzz_xzz_0[i] = ta_zzzz_xzz_0[i] * pa_y[i] - ta_zzzz_xzz_1[i] * pc_y[i];

        ta_yzzzz_yyy_0[i] = 3.0 * ta_zzzz_yy_0[i] * fe_0 - 3.0 * ta_zzzz_yy_1[i] * fe_0 + ta_zzzz_yyy_0[i] * pa_y[i] - ta_zzzz_yyy_1[i] * pc_y[i];

        ta_yzzzz_yyz_0[i] = 2.0 * ta_zzzz_yz_0[i] * fe_0 - 2.0 * ta_zzzz_yz_1[i] * fe_0 + ta_zzzz_yyz_0[i] * pa_y[i] - ta_zzzz_yyz_1[i] * pc_y[i];

        ta_yzzzz_yzz_0[i] = ta_zzzz_zz_0[i] * fe_0 - ta_zzzz_zz_1[i] * fe_0 + ta_zzzz_yzz_0[i] * pa_y[i] - ta_zzzz_yzz_1[i] * pc_y[i];

        ta_yzzzz_zzz_0[i] = ta_zzzz_zzz_0[i] * pa_y[i] - ta_zzzz_zzz_1[i] * pc_y[i];
    }

    // Set up 200-210 components of targeted buffer : HF

    auto ta_zzzzz_xxx_0 = pbuffer.data(idx_npot_0_hf + 200);

    auto ta_zzzzz_xxy_0 = pbuffer.data(idx_npot_0_hf + 201);

    auto ta_zzzzz_xxz_0 = pbuffer.data(idx_npot_0_hf + 202);

    auto ta_zzzzz_xyy_0 = pbuffer.data(idx_npot_0_hf + 203);

    auto ta_zzzzz_xyz_0 = pbuffer.data(idx_npot_0_hf + 204);

    auto ta_zzzzz_xzz_0 = pbuffer.data(idx_npot_0_hf + 205);

    auto ta_zzzzz_yyy_0 = pbuffer.data(idx_npot_0_hf + 206);

    auto ta_zzzzz_yyz_0 = pbuffer.data(idx_npot_0_hf + 207);

    auto ta_zzzzz_yzz_0 = pbuffer.data(idx_npot_0_hf + 208);

    auto ta_zzzzz_zzz_0 = pbuffer.data(idx_npot_0_hf + 209);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta_zzz_xxx_0,   \
                             ta_zzz_xxx_1,   \
                             ta_zzz_xxy_0,   \
                             ta_zzz_xxy_1,   \
                             ta_zzz_xxz_0,   \
                             ta_zzz_xxz_1,   \
                             ta_zzz_xyy_0,   \
                             ta_zzz_xyy_1,   \
                             ta_zzz_xyz_0,   \
                             ta_zzz_xyz_1,   \
                             ta_zzz_xzz_0,   \
                             ta_zzz_xzz_1,   \
                             ta_zzz_yyy_0,   \
                             ta_zzz_yyy_1,   \
                             ta_zzz_yyz_0,   \
                             ta_zzz_yyz_1,   \
                             ta_zzz_yzz_0,   \
                             ta_zzz_yzz_1,   \
                             ta_zzz_zzz_0,   \
                             ta_zzz_zzz_1,   \
                             ta_zzzz_xx_0,   \
                             ta_zzzz_xx_1,   \
                             ta_zzzz_xxx_0,  \
                             ta_zzzz_xxx_1,  \
                             ta_zzzz_xxy_0,  \
                             ta_zzzz_xxy_1,  \
                             ta_zzzz_xxz_0,  \
                             ta_zzzz_xxz_1,  \
                             ta_zzzz_xy_0,   \
                             ta_zzzz_xy_1,   \
                             ta_zzzz_xyy_0,  \
                             ta_zzzz_xyy_1,  \
                             ta_zzzz_xyz_0,  \
                             ta_zzzz_xyz_1,  \
                             ta_zzzz_xz_0,   \
                             ta_zzzz_xz_1,   \
                             ta_zzzz_xzz_0,  \
                             ta_zzzz_xzz_1,  \
                             ta_zzzz_yy_0,   \
                             ta_zzzz_yy_1,   \
                             ta_zzzz_yyy_0,  \
                             ta_zzzz_yyy_1,  \
                             ta_zzzz_yyz_0,  \
                             ta_zzzz_yyz_1,  \
                             ta_zzzz_yz_0,   \
                             ta_zzzz_yz_1,   \
                             ta_zzzz_yzz_0,  \
                             ta_zzzz_yzz_1,  \
                             ta_zzzz_zz_0,   \
                             ta_zzzz_zz_1,   \
                             ta_zzzz_zzz_0,  \
                             ta_zzzz_zzz_1,  \
                             ta_zzzzz_xxx_0, \
                             ta_zzzzz_xxy_0, \
                             ta_zzzzz_xxz_0, \
                             ta_zzzzz_xyy_0, \
                             ta_zzzzz_xyz_0, \
                             ta_zzzzz_xzz_0, \
                             ta_zzzzz_yyy_0, \
                             ta_zzzzz_yyz_0, \
                             ta_zzzzz_yzz_0, \
                             ta_zzzzz_zzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzzz_xxx_0[i] = 4.0 * ta_zzz_xxx_0[i] * fe_0 - 4.0 * ta_zzz_xxx_1[i] * fe_0 + ta_zzzz_xxx_0[i] * pa_z[i] - ta_zzzz_xxx_1[i] * pc_z[i];

        ta_zzzzz_xxy_0[i] = 4.0 * ta_zzz_xxy_0[i] * fe_0 - 4.0 * ta_zzz_xxy_1[i] * fe_0 + ta_zzzz_xxy_0[i] * pa_z[i] - ta_zzzz_xxy_1[i] * pc_z[i];

        ta_zzzzz_xxz_0[i] = 4.0 * ta_zzz_xxz_0[i] * fe_0 - 4.0 * ta_zzz_xxz_1[i] * fe_0 + ta_zzzz_xx_0[i] * fe_0 - ta_zzzz_xx_1[i] * fe_0 +
                            ta_zzzz_xxz_0[i] * pa_z[i] - ta_zzzz_xxz_1[i] * pc_z[i];

        ta_zzzzz_xyy_0[i] = 4.0 * ta_zzz_xyy_0[i] * fe_0 - 4.0 * ta_zzz_xyy_1[i] * fe_0 + ta_zzzz_xyy_0[i] * pa_z[i] - ta_zzzz_xyy_1[i] * pc_z[i];

        ta_zzzzz_xyz_0[i] = 4.0 * ta_zzz_xyz_0[i] * fe_0 - 4.0 * ta_zzz_xyz_1[i] * fe_0 + ta_zzzz_xy_0[i] * fe_0 - ta_zzzz_xy_1[i] * fe_0 +
                            ta_zzzz_xyz_0[i] * pa_z[i] - ta_zzzz_xyz_1[i] * pc_z[i];

        ta_zzzzz_xzz_0[i] = 4.0 * ta_zzz_xzz_0[i] * fe_0 - 4.0 * ta_zzz_xzz_1[i] * fe_0 + 2.0 * ta_zzzz_xz_0[i] * fe_0 -
                            2.0 * ta_zzzz_xz_1[i] * fe_0 + ta_zzzz_xzz_0[i] * pa_z[i] - ta_zzzz_xzz_1[i] * pc_z[i];

        ta_zzzzz_yyy_0[i] = 4.0 * ta_zzz_yyy_0[i] * fe_0 - 4.0 * ta_zzz_yyy_1[i] * fe_0 + ta_zzzz_yyy_0[i] * pa_z[i] - ta_zzzz_yyy_1[i] * pc_z[i];

        ta_zzzzz_yyz_0[i] = 4.0 * ta_zzz_yyz_0[i] * fe_0 - 4.0 * ta_zzz_yyz_1[i] * fe_0 + ta_zzzz_yy_0[i] * fe_0 - ta_zzzz_yy_1[i] * fe_0 +
                            ta_zzzz_yyz_0[i] * pa_z[i] - ta_zzzz_yyz_1[i] * pc_z[i];

        ta_zzzzz_yzz_0[i] = 4.0 * ta_zzz_yzz_0[i] * fe_0 - 4.0 * ta_zzz_yzz_1[i] * fe_0 + 2.0 * ta_zzzz_yz_0[i] * fe_0 -
                            2.0 * ta_zzzz_yz_1[i] * fe_0 + ta_zzzz_yzz_0[i] * pa_z[i] - ta_zzzz_yzz_1[i] * pc_z[i];

        ta_zzzzz_zzz_0[i] = 4.0 * ta_zzz_zzz_0[i] * fe_0 - 4.0 * ta_zzz_zzz_1[i] * fe_0 + 3.0 * ta_zzzz_zz_0[i] * fe_0 -
                            3.0 * ta_zzzz_zz_1[i] * fe_0 + ta_zzzz_zzz_0[i] * pa_z[i] - ta_zzzz_zzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
