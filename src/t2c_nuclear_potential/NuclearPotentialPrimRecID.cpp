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

#include "NuclearPotentialPrimRecID.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_id(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_id,
                               const size_t              idx_npot_0_gd,
                               const size_t              idx_npot_1_gd,
                               const size_t              idx_npot_0_hp,
                               const size_t              idx_npot_1_hp,
                               const size_t              idx_npot_0_hd,
                               const size_t              idx_npot_1_hd,
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

    // Set up components of auxiliary buffer : GD

    auto ta_xxxx_xx_0 = pbuffer.data(idx_npot_0_gd);

    auto ta_xxxx_xy_0 = pbuffer.data(idx_npot_0_gd + 1);

    auto ta_xxxx_xz_0 = pbuffer.data(idx_npot_0_gd + 2);

    auto ta_xxxx_yy_0 = pbuffer.data(idx_npot_0_gd + 3);

    auto ta_xxxx_yz_0 = pbuffer.data(idx_npot_0_gd + 4);

    auto ta_xxxx_zz_0 = pbuffer.data(idx_npot_0_gd + 5);

    auto ta_xxxy_xx_0 = pbuffer.data(idx_npot_0_gd + 6);

    auto ta_xxxy_xz_0 = pbuffer.data(idx_npot_0_gd + 8);

    auto ta_xxxy_yy_0 = pbuffer.data(idx_npot_0_gd + 9);

    auto ta_xxxy_yz_0 = pbuffer.data(idx_npot_0_gd + 10);

    auto ta_xxxz_xx_0 = pbuffer.data(idx_npot_0_gd + 12);

    auto ta_xxxz_xy_0 = pbuffer.data(idx_npot_0_gd + 13);

    auto ta_xxxz_xz_0 = pbuffer.data(idx_npot_0_gd + 14);

    auto ta_xxxz_yz_0 = pbuffer.data(idx_npot_0_gd + 16);

    auto ta_xxxz_zz_0 = pbuffer.data(idx_npot_0_gd + 17);

    auto ta_xxyy_xx_0 = pbuffer.data(idx_npot_0_gd + 18);

    auto ta_xxyy_xy_0 = pbuffer.data(idx_npot_0_gd + 19);

    auto ta_xxyy_xz_0 = pbuffer.data(idx_npot_0_gd + 20);

    auto ta_xxyy_yy_0 = pbuffer.data(idx_npot_0_gd + 21);

    auto ta_xxyy_yz_0 = pbuffer.data(idx_npot_0_gd + 22);

    auto ta_xxyy_zz_0 = pbuffer.data(idx_npot_0_gd + 23);

    auto ta_xxyz_xz_0 = pbuffer.data(idx_npot_0_gd + 26);

    auto ta_xxyz_yz_0 = pbuffer.data(idx_npot_0_gd + 28);

    auto ta_xxzz_xx_0 = pbuffer.data(idx_npot_0_gd + 30);

    auto ta_xxzz_xy_0 = pbuffer.data(idx_npot_0_gd + 31);

    auto ta_xxzz_xz_0 = pbuffer.data(idx_npot_0_gd + 32);

    auto ta_xxzz_yy_0 = pbuffer.data(idx_npot_0_gd + 33);

    auto ta_xxzz_yz_0 = pbuffer.data(idx_npot_0_gd + 34);

    auto ta_xxzz_zz_0 = pbuffer.data(idx_npot_0_gd + 35);

    auto ta_xyyy_xy_0 = pbuffer.data(idx_npot_0_gd + 37);

    auto ta_xyyy_yy_0 = pbuffer.data(idx_npot_0_gd + 39);

    auto ta_xyyy_yz_0 = pbuffer.data(idx_npot_0_gd + 40);

    auto ta_xyyy_zz_0 = pbuffer.data(idx_npot_0_gd + 41);

    auto ta_xyyz_yz_0 = pbuffer.data(idx_npot_0_gd + 46);

    auto ta_xyyz_zz_0 = pbuffer.data(idx_npot_0_gd + 47);

    auto ta_xyzz_yy_0 = pbuffer.data(idx_npot_0_gd + 51);

    auto ta_xyzz_yz_0 = pbuffer.data(idx_npot_0_gd + 52);

    auto ta_xzzz_xz_0 = pbuffer.data(idx_npot_0_gd + 56);

    auto ta_xzzz_yy_0 = pbuffer.data(idx_npot_0_gd + 57);

    auto ta_xzzz_yz_0 = pbuffer.data(idx_npot_0_gd + 58);

    auto ta_xzzz_zz_0 = pbuffer.data(idx_npot_0_gd + 59);

    auto ta_yyyy_xx_0 = pbuffer.data(idx_npot_0_gd + 60);

    auto ta_yyyy_xy_0 = pbuffer.data(idx_npot_0_gd + 61);

    auto ta_yyyy_xz_0 = pbuffer.data(idx_npot_0_gd + 62);

    auto ta_yyyy_yy_0 = pbuffer.data(idx_npot_0_gd + 63);

    auto ta_yyyy_yz_0 = pbuffer.data(idx_npot_0_gd + 64);

    auto ta_yyyy_zz_0 = pbuffer.data(idx_npot_0_gd + 65);

    auto ta_yyyz_xy_0 = pbuffer.data(idx_npot_0_gd + 67);

    auto ta_yyyz_xz_0 = pbuffer.data(idx_npot_0_gd + 68);

    auto ta_yyyz_yy_0 = pbuffer.data(idx_npot_0_gd + 69);

    auto ta_yyyz_yz_0 = pbuffer.data(idx_npot_0_gd + 70);

    auto ta_yyyz_zz_0 = pbuffer.data(idx_npot_0_gd + 71);

    auto ta_yyzz_xx_0 = pbuffer.data(idx_npot_0_gd + 72);

    auto ta_yyzz_xy_0 = pbuffer.data(idx_npot_0_gd + 73);

    auto ta_yyzz_xz_0 = pbuffer.data(idx_npot_0_gd + 74);

    auto ta_yyzz_yy_0 = pbuffer.data(idx_npot_0_gd + 75);

    auto ta_yyzz_yz_0 = pbuffer.data(idx_npot_0_gd + 76);

    auto ta_yyzz_zz_0 = pbuffer.data(idx_npot_0_gd + 77);

    auto ta_yzzz_xx_0 = pbuffer.data(idx_npot_0_gd + 78);

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

    auto ta_xxxy_xx_1 = pbuffer.data(idx_npot_1_gd + 6);

    auto ta_xxxy_xz_1 = pbuffer.data(idx_npot_1_gd + 8);

    auto ta_xxxy_yy_1 = pbuffer.data(idx_npot_1_gd + 9);

    auto ta_xxxy_yz_1 = pbuffer.data(idx_npot_1_gd + 10);

    auto ta_xxxz_xx_1 = pbuffer.data(idx_npot_1_gd + 12);

    auto ta_xxxz_xy_1 = pbuffer.data(idx_npot_1_gd + 13);

    auto ta_xxxz_xz_1 = pbuffer.data(idx_npot_1_gd + 14);

    auto ta_xxxz_yz_1 = pbuffer.data(idx_npot_1_gd + 16);

    auto ta_xxxz_zz_1 = pbuffer.data(idx_npot_1_gd + 17);

    auto ta_xxyy_xx_1 = pbuffer.data(idx_npot_1_gd + 18);

    auto ta_xxyy_xy_1 = pbuffer.data(idx_npot_1_gd + 19);

    auto ta_xxyy_xz_1 = pbuffer.data(idx_npot_1_gd + 20);

    auto ta_xxyy_yy_1 = pbuffer.data(idx_npot_1_gd + 21);

    auto ta_xxyy_yz_1 = pbuffer.data(idx_npot_1_gd + 22);

    auto ta_xxyy_zz_1 = pbuffer.data(idx_npot_1_gd + 23);

    auto ta_xxyz_xz_1 = pbuffer.data(idx_npot_1_gd + 26);

    auto ta_xxyz_yz_1 = pbuffer.data(idx_npot_1_gd + 28);

    auto ta_xxzz_xx_1 = pbuffer.data(idx_npot_1_gd + 30);

    auto ta_xxzz_xy_1 = pbuffer.data(idx_npot_1_gd + 31);

    auto ta_xxzz_xz_1 = pbuffer.data(idx_npot_1_gd + 32);

    auto ta_xxzz_yy_1 = pbuffer.data(idx_npot_1_gd + 33);

    auto ta_xxzz_yz_1 = pbuffer.data(idx_npot_1_gd + 34);

    auto ta_xxzz_zz_1 = pbuffer.data(idx_npot_1_gd + 35);

    auto ta_xyyy_xy_1 = pbuffer.data(idx_npot_1_gd + 37);

    auto ta_xyyy_yy_1 = pbuffer.data(idx_npot_1_gd + 39);

    auto ta_xyyy_yz_1 = pbuffer.data(idx_npot_1_gd + 40);

    auto ta_xyyy_zz_1 = pbuffer.data(idx_npot_1_gd + 41);

    auto ta_xyyz_yz_1 = pbuffer.data(idx_npot_1_gd + 46);

    auto ta_xyyz_zz_1 = pbuffer.data(idx_npot_1_gd + 47);

    auto ta_xyzz_yy_1 = pbuffer.data(idx_npot_1_gd + 51);

    auto ta_xyzz_yz_1 = pbuffer.data(idx_npot_1_gd + 52);

    auto ta_xzzz_xz_1 = pbuffer.data(idx_npot_1_gd + 56);

    auto ta_xzzz_yy_1 = pbuffer.data(idx_npot_1_gd + 57);

    auto ta_xzzz_yz_1 = pbuffer.data(idx_npot_1_gd + 58);

    auto ta_xzzz_zz_1 = pbuffer.data(idx_npot_1_gd + 59);

    auto ta_yyyy_xx_1 = pbuffer.data(idx_npot_1_gd + 60);

    auto ta_yyyy_xy_1 = pbuffer.data(idx_npot_1_gd + 61);

    auto ta_yyyy_xz_1 = pbuffer.data(idx_npot_1_gd + 62);

    auto ta_yyyy_yy_1 = pbuffer.data(idx_npot_1_gd + 63);

    auto ta_yyyy_yz_1 = pbuffer.data(idx_npot_1_gd + 64);

    auto ta_yyyy_zz_1 = pbuffer.data(idx_npot_1_gd + 65);

    auto ta_yyyz_xy_1 = pbuffer.data(idx_npot_1_gd + 67);

    auto ta_yyyz_xz_1 = pbuffer.data(idx_npot_1_gd + 68);

    auto ta_yyyz_yy_1 = pbuffer.data(idx_npot_1_gd + 69);

    auto ta_yyyz_yz_1 = pbuffer.data(idx_npot_1_gd + 70);

    auto ta_yyyz_zz_1 = pbuffer.data(idx_npot_1_gd + 71);

    auto ta_yyzz_xx_1 = pbuffer.data(idx_npot_1_gd + 72);

    auto ta_yyzz_xy_1 = pbuffer.data(idx_npot_1_gd + 73);

    auto ta_yyzz_xz_1 = pbuffer.data(idx_npot_1_gd + 74);

    auto ta_yyzz_yy_1 = pbuffer.data(idx_npot_1_gd + 75);

    auto ta_yyzz_yz_1 = pbuffer.data(idx_npot_1_gd + 76);

    auto ta_yyzz_zz_1 = pbuffer.data(idx_npot_1_gd + 77);

    auto ta_yzzz_xx_1 = pbuffer.data(idx_npot_1_gd + 78);

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

    // Set up components of auxiliary buffer : HP

    auto ta_xxxxx_x_0 = pbuffer.data(idx_npot_0_hp);

    auto ta_xxxxx_y_0 = pbuffer.data(idx_npot_0_hp + 1);

    auto ta_xxxxx_z_0 = pbuffer.data(idx_npot_0_hp + 2);

    auto ta_xxxyy_y_0 = pbuffer.data(idx_npot_0_hp + 10);

    auto ta_xxxzz_x_0 = pbuffer.data(idx_npot_0_hp + 15);

    auto ta_xxxzz_z_0 = pbuffer.data(idx_npot_0_hp + 17);

    auto ta_xxyyy_y_0 = pbuffer.data(idx_npot_0_hp + 19);

    auto ta_xxzzz_x_0 = pbuffer.data(idx_npot_0_hp + 27);

    auto ta_xxzzz_z_0 = pbuffer.data(idx_npot_0_hp + 29);

    auto ta_xyyyy_y_0 = pbuffer.data(idx_npot_0_hp + 31);

    auto ta_xzzzz_z_0 = pbuffer.data(idx_npot_0_hp + 44);

    auto ta_yyyyy_x_0 = pbuffer.data(idx_npot_0_hp + 45);

    auto ta_yyyyy_y_0 = pbuffer.data(idx_npot_0_hp + 46);

    auto ta_yyyyy_z_0 = pbuffer.data(idx_npot_0_hp + 47);

    auto ta_yyyyz_z_0 = pbuffer.data(idx_npot_0_hp + 50);

    auto ta_yyyzz_x_0 = pbuffer.data(idx_npot_0_hp + 51);

    auto ta_yyyzz_y_0 = pbuffer.data(idx_npot_0_hp + 52);

    auto ta_yyyzz_z_0 = pbuffer.data(idx_npot_0_hp + 53);

    auto ta_yyzzz_x_0 = pbuffer.data(idx_npot_0_hp + 54);

    auto ta_yyzzz_y_0 = pbuffer.data(idx_npot_0_hp + 55);

    auto ta_yyzzz_z_0 = pbuffer.data(idx_npot_0_hp + 56);

    auto ta_yzzzz_y_0 = pbuffer.data(idx_npot_0_hp + 58);

    auto ta_yzzzz_z_0 = pbuffer.data(idx_npot_0_hp + 59);

    auto ta_zzzzz_x_0 = pbuffer.data(idx_npot_0_hp + 60);

    auto ta_zzzzz_y_0 = pbuffer.data(idx_npot_0_hp + 61);

    auto ta_zzzzz_z_0 = pbuffer.data(idx_npot_0_hp + 62);

    // Set up components of auxiliary buffer : HP

    auto ta_xxxxx_x_1 = pbuffer.data(idx_npot_1_hp);

    auto ta_xxxxx_y_1 = pbuffer.data(idx_npot_1_hp + 1);

    auto ta_xxxxx_z_1 = pbuffer.data(idx_npot_1_hp + 2);

    auto ta_xxxyy_y_1 = pbuffer.data(idx_npot_1_hp + 10);

    auto ta_xxxzz_x_1 = pbuffer.data(idx_npot_1_hp + 15);

    auto ta_xxxzz_z_1 = pbuffer.data(idx_npot_1_hp + 17);

    auto ta_xxyyy_y_1 = pbuffer.data(idx_npot_1_hp + 19);

    auto ta_xxzzz_x_1 = pbuffer.data(idx_npot_1_hp + 27);

    auto ta_xxzzz_z_1 = pbuffer.data(idx_npot_1_hp + 29);

    auto ta_xyyyy_y_1 = pbuffer.data(idx_npot_1_hp + 31);

    auto ta_xzzzz_z_1 = pbuffer.data(idx_npot_1_hp + 44);

    auto ta_yyyyy_x_1 = pbuffer.data(idx_npot_1_hp + 45);

    auto ta_yyyyy_y_1 = pbuffer.data(idx_npot_1_hp + 46);

    auto ta_yyyyy_z_1 = pbuffer.data(idx_npot_1_hp + 47);

    auto ta_yyyyz_z_1 = pbuffer.data(idx_npot_1_hp + 50);

    auto ta_yyyzz_x_1 = pbuffer.data(idx_npot_1_hp + 51);

    auto ta_yyyzz_y_1 = pbuffer.data(idx_npot_1_hp + 52);

    auto ta_yyyzz_z_1 = pbuffer.data(idx_npot_1_hp + 53);

    auto ta_yyzzz_x_1 = pbuffer.data(idx_npot_1_hp + 54);

    auto ta_yyzzz_y_1 = pbuffer.data(idx_npot_1_hp + 55);

    auto ta_yyzzz_z_1 = pbuffer.data(idx_npot_1_hp + 56);

    auto ta_yzzzz_y_1 = pbuffer.data(idx_npot_1_hp + 58);

    auto ta_yzzzz_z_1 = pbuffer.data(idx_npot_1_hp + 59);

    auto ta_zzzzz_x_1 = pbuffer.data(idx_npot_1_hp + 60);

    auto ta_zzzzz_y_1 = pbuffer.data(idx_npot_1_hp + 61);

    auto ta_zzzzz_z_1 = pbuffer.data(idx_npot_1_hp + 62);

    // Set up components of auxiliary buffer : HD

    auto ta_xxxxx_xx_0 = pbuffer.data(idx_npot_0_hd);

    auto ta_xxxxx_xy_0 = pbuffer.data(idx_npot_0_hd + 1);

    auto ta_xxxxx_xz_0 = pbuffer.data(idx_npot_0_hd + 2);

    auto ta_xxxxx_yy_0 = pbuffer.data(idx_npot_0_hd + 3);

    auto ta_xxxxx_yz_0 = pbuffer.data(idx_npot_0_hd + 4);

    auto ta_xxxxx_zz_0 = pbuffer.data(idx_npot_0_hd + 5);

    auto ta_xxxxy_xx_0 = pbuffer.data(idx_npot_0_hd + 6);

    auto ta_xxxxy_xy_0 = pbuffer.data(idx_npot_0_hd + 7);

    auto ta_xxxxy_xz_0 = pbuffer.data(idx_npot_0_hd + 8);

    auto ta_xxxxy_yy_0 = pbuffer.data(idx_npot_0_hd + 9);

    auto ta_xxxxy_yz_0 = pbuffer.data(idx_npot_0_hd + 10);

    auto ta_xxxxz_xx_0 = pbuffer.data(idx_npot_0_hd + 12);

    auto ta_xxxxz_xy_0 = pbuffer.data(idx_npot_0_hd + 13);

    auto ta_xxxxz_xz_0 = pbuffer.data(idx_npot_0_hd + 14);

    auto ta_xxxxz_yz_0 = pbuffer.data(idx_npot_0_hd + 16);

    auto ta_xxxxz_zz_0 = pbuffer.data(idx_npot_0_hd + 17);

    auto ta_xxxyy_xx_0 = pbuffer.data(idx_npot_0_hd + 18);

    auto ta_xxxyy_xy_0 = pbuffer.data(idx_npot_0_hd + 19);

    auto ta_xxxyy_xz_0 = pbuffer.data(idx_npot_0_hd + 20);

    auto ta_xxxyy_yy_0 = pbuffer.data(idx_npot_0_hd + 21);

    auto ta_xxxyy_yz_0 = pbuffer.data(idx_npot_0_hd + 22);

    auto ta_xxxyy_zz_0 = pbuffer.data(idx_npot_0_hd + 23);

    auto ta_xxxyz_xz_0 = pbuffer.data(idx_npot_0_hd + 26);

    auto ta_xxxyz_yz_0 = pbuffer.data(idx_npot_0_hd + 28);

    auto ta_xxxzz_xx_0 = pbuffer.data(idx_npot_0_hd + 30);

    auto ta_xxxzz_xy_0 = pbuffer.data(idx_npot_0_hd + 31);

    auto ta_xxxzz_xz_0 = pbuffer.data(idx_npot_0_hd + 32);

    auto ta_xxxzz_yy_0 = pbuffer.data(idx_npot_0_hd + 33);

    auto ta_xxxzz_yz_0 = pbuffer.data(idx_npot_0_hd + 34);

    auto ta_xxxzz_zz_0 = pbuffer.data(idx_npot_0_hd + 35);

    auto ta_xxyyy_xx_0 = pbuffer.data(idx_npot_0_hd + 36);

    auto ta_xxyyy_xy_0 = pbuffer.data(idx_npot_0_hd + 37);

    auto ta_xxyyy_xz_0 = pbuffer.data(idx_npot_0_hd + 38);

    auto ta_xxyyy_yy_0 = pbuffer.data(idx_npot_0_hd + 39);

    auto ta_xxyyy_yz_0 = pbuffer.data(idx_npot_0_hd + 40);

    auto ta_xxyyy_zz_0 = pbuffer.data(idx_npot_0_hd + 41);

    auto ta_xxyyz_xy_0 = pbuffer.data(idx_npot_0_hd + 43);

    auto ta_xxyyz_xz_0 = pbuffer.data(idx_npot_0_hd + 44);

    auto ta_xxyyz_yz_0 = pbuffer.data(idx_npot_0_hd + 46);

    auto ta_xxyyz_zz_0 = pbuffer.data(idx_npot_0_hd + 47);

    auto ta_xxyzz_xx_0 = pbuffer.data(idx_npot_0_hd + 48);

    auto ta_xxyzz_xz_0 = pbuffer.data(idx_npot_0_hd + 50);

    auto ta_xxyzz_yy_0 = pbuffer.data(idx_npot_0_hd + 51);

    auto ta_xxyzz_yz_0 = pbuffer.data(idx_npot_0_hd + 52);

    auto ta_xxzzz_xx_0 = pbuffer.data(idx_npot_0_hd + 54);

    auto ta_xxzzz_xy_0 = pbuffer.data(idx_npot_0_hd + 55);

    auto ta_xxzzz_xz_0 = pbuffer.data(idx_npot_0_hd + 56);

    auto ta_xxzzz_yy_0 = pbuffer.data(idx_npot_0_hd + 57);

    auto ta_xxzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 58);

    auto ta_xxzzz_zz_0 = pbuffer.data(idx_npot_0_hd + 59);

    auto ta_xyyyy_xx_0 = pbuffer.data(idx_npot_0_hd + 60);

    auto ta_xyyyy_xy_0 = pbuffer.data(idx_npot_0_hd + 61);

    auto ta_xyyyy_yy_0 = pbuffer.data(idx_npot_0_hd + 63);

    auto ta_xyyyy_yz_0 = pbuffer.data(idx_npot_0_hd + 64);

    auto ta_xyyyy_zz_0 = pbuffer.data(idx_npot_0_hd + 65);

    auto ta_xyyyz_yz_0 = pbuffer.data(idx_npot_0_hd + 70);

    auto ta_xyyyz_zz_0 = pbuffer.data(idx_npot_0_hd + 71);

    auto ta_xyyzz_yy_0 = pbuffer.data(idx_npot_0_hd + 75);

    auto ta_xyyzz_yz_0 = pbuffer.data(idx_npot_0_hd + 76);

    auto ta_xyyzz_zz_0 = pbuffer.data(idx_npot_0_hd + 77);

    auto ta_xyzzz_yy_0 = pbuffer.data(idx_npot_0_hd + 81);

    auto ta_xyzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 82);

    auto ta_xzzzz_xx_0 = pbuffer.data(idx_npot_0_hd + 84);

    auto ta_xzzzz_xz_0 = pbuffer.data(idx_npot_0_hd + 86);

    auto ta_xzzzz_yy_0 = pbuffer.data(idx_npot_0_hd + 87);

    auto ta_xzzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 88);

    auto ta_xzzzz_zz_0 = pbuffer.data(idx_npot_0_hd + 89);

    auto ta_yyyyy_xx_0 = pbuffer.data(idx_npot_0_hd + 90);

    auto ta_yyyyy_xy_0 = pbuffer.data(idx_npot_0_hd + 91);

    auto ta_yyyyy_xz_0 = pbuffer.data(idx_npot_0_hd + 92);

    auto ta_yyyyy_yy_0 = pbuffer.data(idx_npot_0_hd + 93);

    auto ta_yyyyy_yz_0 = pbuffer.data(idx_npot_0_hd + 94);

    auto ta_yyyyy_zz_0 = pbuffer.data(idx_npot_0_hd + 95);

    auto ta_yyyyz_xy_0 = pbuffer.data(idx_npot_0_hd + 97);

    auto ta_yyyyz_xz_0 = pbuffer.data(idx_npot_0_hd + 98);

    auto ta_yyyyz_yy_0 = pbuffer.data(idx_npot_0_hd + 99);

    auto ta_yyyyz_yz_0 = pbuffer.data(idx_npot_0_hd + 100);

    auto ta_yyyyz_zz_0 = pbuffer.data(idx_npot_0_hd + 101);

    auto ta_yyyzz_xx_0 = pbuffer.data(idx_npot_0_hd + 102);

    auto ta_yyyzz_xy_0 = pbuffer.data(idx_npot_0_hd + 103);

    auto ta_yyyzz_xz_0 = pbuffer.data(idx_npot_0_hd + 104);

    auto ta_yyyzz_yy_0 = pbuffer.data(idx_npot_0_hd + 105);

    auto ta_yyyzz_yz_0 = pbuffer.data(idx_npot_0_hd + 106);

    auto ta_yyyzz_zz_0 = pbuffer.data(idx_npot_0_hd + 107);

    auto ta_yyzzz_xx_0 = pbuffer.data(idx_npot_0_hd + 108);

    auto ta_yyzzz_xy_0 = pbuffer.data(idx_npot_0_hd + 109);

    auto ta_yyzzz_xz_0 = pbuffer.data(idx_npot_0_hd + 110);

    auto ta_yyzzz_yy_0 = pbuffer.data(idx_npot_0_hd + 111);

    auto ta_yyzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 112);

    auto ta_yyzzz_zz_0 = pbuffer.data(idx_npot_0_hd + 113);

    auto ta_yzzzz_xx_0 = pbuffer.data(idx_npot_0_hd + 114);

    auto ta_yzzzz_xy_0 = pbuffer.data(idx_npot_0_hd + 115);

    auto ta_yzzzz_xz_0 = pbuffer.data(idx_npot_0_hd + 116);

    auto ta_yzzzz_yy_0 = pbuffer.data(idx_npot_0_hd + 117);

    auto ta_yzzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 118);

    auto ta_yzzzz_zz_0 = pbuffer.data(idx_npot_0_hd + 119);

    auto ta_zzzzz_xx_0 = pbuffer.data(idx_npot_0_hd + 120);

    auto ta_zzzzz_xy_0 = pbuffer.data(idx_npot_0_hd + 121);

    auto ta_zzzzz_xz_0 = pbuffer.data(idx_npot_0_hd + 122);

    auto ta_zzzzz_yy_0 = pbuffer.data(idx_npot_0_hd + 123);

    auto ta_zzzzz_yz_0 = pbuffer.data(idx_npot_0_hd + 124);

    auto ta_zzzzz_zz_0 = pbuffer.data(idx_npot_0_hd + 125);

    // Set up components of auxiliary buffer : HD

    auto ta_xxxxx_xx_1 = pbuffer.data(idx_npot_1_hd);

    auto ta_xxxxx_xy_1 = pbuffer.data(idx_npot_1_hd + 1);

    auto ta_xxxxx_xz_1 = pbuffer.data(idx_npot_1_hd + 2);

    auto ta_xxxxx_yy_1 = pbuffer.data(idx_npot_1_hd + 3);

    auto ta_xxxxx_yz_1 = pbuffer.data(idx_npot_1_hd + 4);

    auto ta_xxxxx_zz_1 = pbuffer.data(idx_npot_1_hd + 5);

    auto ta_xxxxy_xx_1 = pbuffer.data(idx_npot_1_hd + 6);

    auto ta_xxxxy_xy_1 = pbuffer.data(idx_npot_1_hd + 7);

    auto ta_xxxxy_xz_1 = pbuffer.data(idx_npot_1_hd + 8);

    auto ta_xxxxy_yy_1 = pbuffer.data(idx_npot_1_hd + 9);

    auto ta_xxxxy_yz_1 = pbuffer.data(idx_npot_1_hd + 10);

    auto ta_xxxxz_xx_1 = pbuffer.data(idx_npot_1_hd + 12);

    auto ta_xxxxz_xy_1 = pbuffer.data(idx_npot_1_hd + 13);

    auto ta_xxxxz_xz_1 = pbuffer.data(idx_npot_1_hd + 14);

    auto ta_xxxxz_yz_1 = pbuffer.data(idx_npot_1_hd + 16);

    auto ta_xxxxz_zz_1 = pbuffer.data(idx_npot_1_hd + 17);

    auto ta_xxxyy_xx_1 = pbuffer.data(idx_npot_1_hd + 18);

    auto ta_xxxyy_xy_1 = pbuffer.data(idx_npot_1_hd + 19);

    auto ta_xxxyy_xz_1 = pbuffer.data(idx_npot_1_hd + 20);

    auto ta_xxxyy_yy_1 = pbuffer.data(idx_npot_1_hd + 21);

    auto ta_xxxyy_yz_1 = pbuffer.data(idx_npot_1_hd + 22);

    auto ta_xxxyy_zz_1 = pbuffer.data(idx_npot_1_hd + 23);

    auto ta_xxxyz_xz_1 = pbuffer.data(idx_npot_1_hd + 26);

    auto ta_xxxyz_yz_1 = pbuffer.data(idx_npot_1_hd + 28);

    auto ta_xxxzz_xx_1 = pbuffer.data(idx_npot_1_hd + 30);

    auto ta_xxxzz_xy_1 = pbuffer.data(idx_npot_1_hd + 31);

    auto ta_xxxzz_xz_1 = pbuffer.data(idx_npot_1_hd + 32);

    auto ta_xxxzz_yy_1 = pbuffer.data(idx_npot_1_hd + 33);

    auto ta_xxxzz_yz_1 = pbuffer.data(idx_npot_1_hd + 34);

    auto ta_xxxzz_zz_1 = pbuffer.data(idx_npot_1_hd + 35);

    auto ta_xxyyy_xx_1 = pbuffer.data(idx_npot_1_hd + 36);

    auto ta_xxyyy_xy_1 = pbuffer.data(idx_npot_1_hd + 37);

    auto ta_xxyyy_xz_1 = pbuffer.data(idx_npot_1_hd + 38);

    auto ta_xxyyy_yy_1 = pbuffer.data(idx_npot_1_hd + 39);

    auto ta_xxyyy_yz_1 = pbuffer.data(idx_npot_1_hd + 40);

    auto ta_xxyyy_zz_1 = pbuffer.data(idx_npot_1_hd + 41);

    auto ta_xxyyz_xy_1 = pbuffer.data(idx_npot_1_hd + 43);

    auto ta_xxyyz_xz_1 = pbuffer.data(idx_npot_1_hd + 44);

    auto ta_xxyyz_yz_1 = pbuffer.data(idx_npot_1_hd + 46);

    auto ta_xxyyz_zz_1 = pbuffer.data(idx_npot_1_hd + 47);

    auto ta_xxyzz_xx_1 = pbuffer.data(idx_npot_1_hd + 48);

    auto ta_xxyzz_xz_1 = pbuffer.data(idx_npot_1_hd + 50);

    auto ta_xxyzz_yy_1 = pbuffer.data(idx_npot_1_hd + 51);

    auto ta_xxyzz_yz_1 = pbuffer.data(idx_npot_1_hd + 52);

    auto ta_xxzzz_xx_1 = pbuffer.data(idx_npot_1_hd + 54);

    auto ta_xxzzz_xy_1 = pbuffer.data(idx_npot_1_hd + 55);

    auto ta_xxzzz_xz_1 = pbuffer.data(idx_npot_1_hd + 56);

    auto ta_xxzzz_yy_1 = pbuffer.data(idx_npot_1_hd + 57);

    auto ta_xxzzz_yz_1 = pbuffer.data(idx_npot_1_hd + 58);

    auto ta_xxzzz_zz_1 = pbuffer.data(idx_npot_1_hd + 59);

    auto ta_xyyyy_xx_1 = pbuffer.data(idx_npot_1_hd + 60);

    auto ta_xyyyy_xy_1 = pbuffer.data(idx_npot_1_hd + 61);

    auto ta_xyyyy_yy_1 = pbuffer.data(idx_npot_1_hd + 63);

    auto ta_xyyyy_yz_1 = pbuffer.data(idx_npot_1_hd + 64);

    auto ta_xyyyy_zz_1 = pbuffer.data(idx_npot_1_hd + 65);

    auto ta_xyyyz_yz_1 = pbuffer.data(idx_npot_1_hd + 70);

    auto ta_xyyyz_zz_1 = pbuffer.data(idx_npot_1_hd + 71);

    auto ta_xyyzz_yy_1 = pbuffer.data(idx_npot_1_hd + 75);

    auto ta_xyyzz_yz_1 = pbuffer.data(idx_npot_1_hd + 76);

    auto ta_xyyzz_zz_1 = pbuffer.data(idx_npot_1_hd + 77);

    auto ta_xyzzz_yy_1 = pbuffer.data(idx_npot_1_hd + 81);

    auto ta_xyzzz_yz_1 = pbuffer.data(idx_npot_1_hd + 82);

    auto ta_xzzzz_xx_1 = pbuffer.data(idx_npot_1_hd + 84);

    auto ta_xzzzz_xz_1 = pbuffer.data(idx_npot_1_hd + 86);

    auto ta_xzzzz_yy_1 = pbuffer.data(idx_npot_1_hd + 87);

    auto ta_xzzzz_yz_1 = pbuffer.data(idx_npot_1_hd + 88);

    auto ta_xzzzz_zz_1 = pbuffer.data(idx_npot_1_hd + 89);

    auto ta_yyyyy_xx_1 = pbuffer.data(idx_npot_1_hd + 90);

    auto ta_yyyyy_xy_1 = pbuffer.data(idx_npot_1_hd + 91);

    auto ta_yyyyy_xz_1 = pbuffer.data(idx_npot_1_hd + 92);

    auto ta_yyyyy_yy_1 = pbuffer.data(idx_npot_1_hd + 93);

    auto ta_yyyyy_yz_1 = pbuffer.data(idx_npot_1_hd + 94);

    auto ta_yyyyy_zz_1 = pbuffer.data(idx_npot_1_hd + 95);

    auto ta_yyyyz_xy_1 = pbuffer.data(idx_npot_1_hd + 97);

    auto ta_yyyyz_xz_1 = pbuffer.data(idx_npot_1_hd + 98);

    auto ta_yyyyz_yy_1 = pbuffer.data(idx_npot_1_hd + 99);

    auto ta_yyyyz_yz_1 = pbuffer.data(idx_npot_1_hd + 100);

    auto ta_yyyyz_zz_1 = pbuffer.data(idx_npot_1_hd + 101);

    auto ta_yyyzz_xx_1 = pbuffer.data(idx_npot_1_hd + 102);

    auto ta_yyyzz_xy_1 = pbuffer.data(idx_npot_1_hd + 103);

    auto ta_yyyzz_xz_1 = pbuffer.data(idx_npot_1_hd + 104);

    auto ta_yyyzz_yy_1 = pbuffer.data(idx_npot_1_hd + 105);

    auto ta_yyyzz_yz_1 = pbuffer.data(idx_npot_1_hd + 106);

    auto ta_yyyzz_zz_1 = pbuffer.data(idx_npot_1_hd + 107);

    auto ta_yyzzz_xx_1 = pbuffer.data(idx_npot_1_hd + 108);

    auto ta_yyzzz_xy_1 = pbuffer.data(idx_npot_1_hd + 109);

    auto ta_yyzzz_xz_1 = pbuffer.data(idx_npot_1_hd + 110);

    auto ta_yyzzz_yy_1 = pbuffer.data(idx_npot_1_hd + 111);

    auto ta_yyzzz_yz_1 = pbuffer.data(idx_npot_1_hd + 112);

    auto ta_yyzzz_zz_1 = pbuffer.data(idx_npot_1_hd + 113);

    auto ta_yzzzz_xx_1 = pbuffer.data(idx_npot_1_hd + 114);

    auto ta_yzzzz_xy_1 = pbuffer.data(idx_npot_1_hd + 115);

    auto ta_yzzzz_xz_1 = pbuffer.data(idx_npot_1_hd + 116);

    auto ta_yzzzz_yy_1 = pbuffer.data(idx_npot_1_hd + 117);

    auto ta_yzzzz_yz_1 = pbuffer.data(idx_npot_1_hd + 118);

    auto ta_yzzzz_zz_1 = pbuffer.data(idx_npot_1_hd + 119);

    auto ta_zzzzz_xx_1 = pbuffer.data(idx_npot_1_hd + 120);

    auto ta_zzzzz_xy_1 = pbuffer.data(idx_npot_1_hd + 121);

    auto ta_zzzzz_xz_1 = pbuffer.data(idx_npot_1_hd + 122);

    auto ta_zzzzz_yy_1 = pbuffer.data(idx_npot_1_hd + 123);

    auto ta_zzzzz_yz_1 = pbuffer.data(idx_npot_1_hd + 124);

    auto ta_zzzzz_zz_1 = pbuffer.data(idx_npot_1_hd + 125);

    // Set up 0-6 components of targeted buffer : ID

    auto ta_xxxxxx_xx_0 = pbuffer.data(idx_npot_0_id);

    auto ta_xxxxxx_xy_0 = pbuffer.data(idx_npot_0_id + 1);

    auto ta_xxxxxx_xz_0 = pbuffer.data(idx_npot_0_id + 2);

    auto ta_xxxxxx_yy_0 = pbuffer.data(idx_npot_0_id + 3);

    auto ta_xxxxxx_yz_0 = pbuffer.data(idx_npot_0_id + 4);

    auto ta_xxxxxx_zz_0 = pbuffer.data(idx_npot_0_id + 5);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta_xxxx_xx_0,   \
                             ta_xxxx_xx_1,   \
                             ta_xxxx_xy_0,   \
                             ta_xxxx_xy_1,   \
                             ta_xxxx_xz_0,   \
                             ta_xxxx_xz_1,   \
                             ta_xxxx_yy_0,   \
                             ta_xxxx_yy_1,   \
                             ta_xxxx_yz_0,   \
                             ta_xxxx_yz_1,   \
                             ta_xxxx_zz_0,   \
                             ta_xxxx_zz_1,   \
                             ta_xxxxx_x_0,   \
                             ta_xxxxx_x_1,   \
                             ta_xxxxx_xx_0,  \
                             ta_xxxxx_xx_1,  \
                             ta_xxxxx_xy_0,  \
                             ta_xxxxx_xy_1,  \
                             ta_xxxxx_xz_0,  \
                             ta_xxxxx_xz_1,  \
                             ta_xxxxx_y_0,   \
                             ta_xxxxx_y_1,   \
                             ta_xxxxx_yy_0,  \
                             ta_xxxxx_yy_1,  \
                             ta_xxxxx_yz_0,  \
                             ta_xxxxx_yz_1,  \
                             ta_xxxxx_z_0,   \
                             ta_xxxxx_z_1,   \
                             ta_xxxxx_zz_0,  \
                             ta_xxxxx_zz_1,  \
                             ta_xxxxxx_xx_0, \
                             ta_xxxxxx_xy_0, \
                             ta_xxxxxx_xz_0, \
                             ta_xxxxxx_yy_0, \
                             ta_xxxxxx_yz_0, \
                             ta_xxxxxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxx_xx_0[i] = 5.0 * ta_xxxx_xx_0[i] * fe_0 - 5.0 * ta_xxxx_xx_1[i] * fe_0 + 2.0 * ta_xxxxx_x_0[i] * fe_0 -
                            2.0 * ta_xxxxx_x_1[i] * fe_0 + ta_xxxxx_xx_0[i] * pa_x[i] - ta_xxxxx_xx_1[i] * pc_x[i];

        ta_xxxxxx_xy_0[i] = 5.0 * ta_xxxx_xy_0[i] * fe_0 - 5.0 * ta_xxxx_xy_1[i] * fe_0 + ta_xxxxx_y_0[i] * fe_0 - ta_xxxxx_y_1[i] * fe_0 +
                            ta_xxxxx_xy_0[i] * pa_x[i] - ta_xxxxx_xy_1[i] * pc_x[i];

        ta_xxxxxx_xz_0[i] = 5.0 * ta_xxxx_xz_0[i] * fe_0 - 5.0 * ta_xxxx_xz_1[i] * fe_0 + ta_xxxxx_z_0[i] * fe_0 - ta_xxxxx_z_1[i] * fe_0 +
                            ta_xxxxx_xz_0[i] * pa_x[i] - ta_xxxxx_xz_1[i] * pc_x[i];

        ta_xxxxxx_yy_0[i] = 5.0 * ta_xxxx_yy_0[i] * fe_0 - 5.0 * ta_xxxx_yy_1[i] * fe_0 + ta_xxxxx_yy_0[i] * pa_x[i] - ta_xxxxx_yy_1[i] * pc_x[i];

        ta_xxxxxx_yz_0[i] = 5.0 * ta_xxxx_yz_0[i] * fe_0 - 5.0 * ta_xxxx_yz_1[i] * fe_0 + ta_xxxxx_yz_0[i] * pa_x[i] - ta_xxxxx_yz_1[i] * pc_x[i];

        ta_xxxxxx_zz_0[i] = 5.0 * ta_xxxx_zz_0[i] * fe_0 - 5.0 * ta_xxxx_zz_1[i] * fe_0 + ta_xxxxx_zz_0[i] * pa_x[i] - ta_xxxxx_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : ID

    auto ta_xxxxxy_xx_0 = pbuffer.data(idx_npot_0_id + 6);

    auto ta_xxxxxy_xy_0 = pbuffer.data(idx_npot_0_id + 7);

    auto ta_xxxxxy_xz_0 = pbuffer.data(idx_npot_0_id + 8);

    auto ta_xxxxxy_yy_0 = pbuffer.data(idx_npot_0_id + 9);

    auto ta_xxxxxy_yz_0 = pbuffer.data(idx_npot_0_id + 10);

    auto ta_xxxxxy_zz_0 = pbuffer.data(idx_npot_0_id + 11);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta_xxxxx_x_0,   \
                             ta_xxxxx_x_1,   \
                             ta_xxxxx_xx_0,  \
                             ta_xxxxx_xx_1,  \
                             ta_xxxxx_xy_0,  \
                             ta_xxxxx_xy_1,  \
                             ta_xxxxx_xz_0,  \
                             ta_xxxxx_xz_1,  \
                             ta_xxxxx_zz_0,  \
                             ta_xxxxx_zz_1,  \
                             ta_xxxxxy_xx_0, \
                             ta_xxxxxy_xy_0, \
                             ta_xxxxxy_xz_0, \
                             ta_xxxxxy_yy_0, \
                             ta_xxxxxy_yz_0, \
                             ta_xxxxxy_zz_0, \
                             ta_xxxxy_yy_0,  \
                             ta_xxxxy_yy_1,  \
                             ta_xxxxy_yz_0,  \
                             ta_xxxxy_yz_1,  \
                             ta_xxxy_yy_0,   \
                             ta_xxxy_yy_1,   \
                             ta_xxxy_yz_0,   \
                             ta_xxxy_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxy_xx_0[i] = ta_xxxxx_xx_0[i] * pa_y[i] - ta_xxxxx_xx_1[i] * pc_y[i];

        ta_xxxxxy_xy_0[i] = ta_xxxxx_x_0[i] * fe_0 - ta_xxxxx_x_1[i] * fe_0 + ta_xxxxx_xy_0[i] * pa_y[i] - ta_xxxxx_xy_1[i] * pc_y[i];

        ta_xxxxxy_xz_0[i] = ta_xxxxx_xz_0[i] * pa_y[i] - ta_xxxxx_xz_1[i] * pc_y[i];

        ta_xxxxxy_yy_0[i] = 4.0 * ta_xxxy_yy_0[i] * fe_0 - 4.0 * ta_xxxy_yy_1[i] * fe_0 + ta_xxxxy_yy_0[i] * pa_x[i] - ta_xxxxy_yy_1[i] * pc_x[i];

        ta_xxxxxy_yz_0[i] = 4.0 * ta_xxxy_yz_0[i] * fe_0 - 4.0 * ta_xxxy_yz_1[i] * fe_0 + ta_xxxxy_yz_0[i] * pa_x[i] - ta_xxxxy_yz_1[i] * pc_x[i];

        ta_xxxxxy_zz_0[i] = ta_xxxxx_zz_0[i] * pa_y[i] - ta_xxxxx_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : ID

    auto ta_xxxxxz_xx_0 = pbuffer.data(idx_npot_0_id + 12);

    auto ta_xxxxxz_xy_0 = pbuffer.data(idx_npot_0_id + 13);

    auto ta_xxxxxz_xz_0 = pbuffer.data(idx_npot_0_id + 14);

    auto ta_xxxxxz_yy_0 = pbuffer.data(idx_npot_0_id + 15);

    auto ta_xxxxxz_yz_0 = pbuffer.data(idx_npot_0_id + 16);

    auto ta_xxxxxz_zz_0 = pbuffer.data(idx_npot_0_id + 17);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta_xxxxx_x_0,   \
                             ta_xxxxx_x_1,   \
                             ta_xxxxx_xx_0,  \
                             ta_xxxxx_xx_1,  \
                             ta_xxxxx_xy_0,  \
                             ta_xxxxx_xy_1,  \
                             ta_xxxxx_xz_0,  \
                             ta_xxxxx_xz_1,  \
                             ta_xxxxx_yy_0,  \
                             ta_xxxxx_yy_1,  \
                             ta_xxxxxz_xx_0, \
                             ta_xxxxxz_xy_0, \
                             ta_xxxxxz_xz_0, \
                             ta_xxxxxz_yy_0, \
                             ta_xxxxxz_yz_0, \
                             ta_xxxxxz_zz_0, \
                             ta_xxxxz_yz_0,  \
                             ta_xxxxz_yz_1,  \
                             ta_xxxxz_zz_0,  \
                             ta_xxxxz_zz_1,  \
                             ta_xxxz_yz_0,   \
                             ta_xxxz_yz_1,   \
                             ta_xxxz_zz_0,   \
                             ta_xxxz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxxz_xx_0[i] = ta_xxxxx_xx_0[i] * pa_z[i] - ta_xxxxx_xx_1[i] * pc_z[i];

        ta_xxxxxz_xy_0[i] = ta_xxxxx_xy_0[i] * pa_z[i] - ta_xxxxx_xy_1[i] * pc_z[i];

        ta_xxxxxz_xz_0[i] = ta_xxxxx_x_0[i] * fe_0 - ta_xxxxx_x_1[i] * fe_0 + ta_xxxxx_xz_0[i] * pa_z[i] - ta_xxxxx_xz_1[i] * pc_z[i];

        ta_xxxxxz_yy_0[i] = ta_xxxxx_yy_0[i] * pa_z[i] - ta_xxxxx_yy_1[i] * pc_z[i];

        ta_xxxxxz_yz_0[i] = 4.0 * ta_xxxz_yz_0[i] * fe_0 - 4.0 * ta_xxxz_yz_1[i] * fe_0 + ta_xxxxz_yz_0[i] * pa_x[i] - ta_xxxxz_yz_1[i] * pc_x[i];

        ta_xxxxxz_zz_0[i] = 4.0 * ta_xxxz_zz_0[i] * fe_0 - 4.0 * ta_xxxz_zz_1[i] * fe_0 + ta_xxxxz_zz_0[i] * pa_x[i] - ta_xxxxz_zz_1[i] * pc_x[i];
    }

    // Set up 18-24 components of targeted buffer : ID

    auto ta_xxxxyy_xx_0 = pbuffer.data(idx_npot_0_id + 18);

    auto ta_xxxxyy_xy_0 = pbuffer.data(idx_npot_0_id + 19);

    auto ta_xxxxyy_xz_0 = pbuffer.data(idx_npot_0_id + 20);

    auto ta_xxxxyy_yy_0 = pbuffer.data(idx_npot_0_id + 21);

    auto ta_xxxxyy_yz_0 = pbuffer.data(idx_npot_0_id + 22);

    auto ta_xxxxyy_zz_0 = pbuffer.data(idx_npot_0_id + 23);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta_xxxx_xx_0,   \
                             ta_xxxx_xx_1,   \
                             ta_xxxx_xz_0,   \
                             ta_xxxx_xz_1,   \
                             ta_xxxxy_xx_0,  \
                             ta_xxxxy_xx_1,  \
                             ta_xxxxy_xz_0,  \
                             ta_xxxxy_xz_1,  \
                             ta_xxxxyy_xx_0, \
                             ta_xxxxyy_xy_0, \
                             ta_xxxxyy_xz_0, \
                             ta_xxxxyy_yy_0, \
                             ta_xxxxyy_yz_0, \
                             ta_xxxxyy_zz_0, \
                             ta_xxxyy_xy_0,  \
                             ta_xxxyy_xy_1,  \
                             ta_xxxyy_y_0,   \
                             ta_xxxyy_y_1,   \
                             ta_xxxyy_yy_0,  \
                             ta_xxxyy_yy_1,  \
                             ta_xxxyy_yz_0,  \
                             ta_xxxyy_yz_1,  \
                             ta_xxxyy_zz_0,  \
                             ta_xxxyy_zz_1,  \
                             ta_xxyy_xy_0,   \
                             ta_xxyy_xy_1,   \
                             ta_xxyy_yy_0,   \
                             ta_xxyy_yy_1,   \
                             ta_xxyy_yz_0,   \
                             ta_xxyy_yz_1,   \
                             ta_xxyy_zz_0,   \
                             ta_xxyy_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxyy_xx_0[i] = ta_xxxx_xx_0[i] * fe_0 - ta_xxxx_xx_1[i] * fe_0 + ta_xxxxy_xx_0[i] * pa_y[i] - ta_xxxxy_xx_1[i] * pc_y[i];

        ta_xxxxyy_xy_0[i] = 3.0 * ta_xxyy_xy_0[i] * fe_0 - 3.0 * ta_xxyy_xy_1[i] * fe_0 + ta_xxxyy_y_0[i] * fe_0 - ta_xxxyy_y_1[i] * fe_0 +
                            ta_xxxyy_xy_0[i] * pa_x[i] - ta_xxxyy_xy_1[i] * pc_x[i];

        ta_xxxxyy_xz_0[i] = ta_xxxx_xz_0[i] * fe_0 - ta_xxxx_xz_1[i] * fe_0 + ta_xxxxy_xz_0[i] * pa_y[i] - ta_xxxxy_xz_1[i] * pc_y[i];

        ta_xxxxyy_yy_0[i] = 3.0 * ta_xxyy_yy_0[i] * fe_0 - 3.0 * ta_xxyy_yy_1[i] * fe_0 + ta_xxxyy_yy_0[i] * pa_x[i] - ta_xxxyy_yy_1[i] * pc_x[i];

        ta_xxxxyy_yz_0[i] = 3.0 * ta_xxyy_yz_0[i] * fe_0 - 3.0 * ta_xxyy_yz_1[i] * fe_0 + ta_xxxyy_yz_0[i] * pa_x[i] - ta_xxxyy_yz_1[i] * pc_x[i];

        ta_xxxxyy_zz_0[i] = 3.0 * ta_xxyy_zz_0[i] * fe_0 - 3.0 * ta_xxyy_zz_1[i] * fe_0 + ta_xxxyy_zz_0[i] * pa_x[i] - ta_xxxyy_zz_1[i] * pc_x[i];
    }

    // Set up 24-30 components of targeted buffer : ID

    auto ta_xxxxyz_xx_0 = pbuffer.data(idx_npot_0_id + 24);

    auto ta_xxxxyz_xy_0 = pbuffer.data(idx_npot_0_id + 25);

    auto ta_xxxxyz_xz_0 = pbuffer.data(idx_npot_0_id + 26);

    auto ta_xxxxyz_yy_0 = pbuffer.data(idx_npot_0_id + 27);

    auto ta_xxxxyz_yz_0 = pbuffer.data(idx_npot_0_id + 28);

    auto ta_xxxxyz_zz_0 = pbuffer.data(idx_npot_0_id + 29);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             pc_x,           \
                             pc_y,           \
                             pc_z,           \
                             ta_xxxxy_xy_0,  \
                             ta_xxxxy_xy_1,  \
                             ta_xxxxy_yy_0,  \
                             ta_xxxxy_yy_1,  \
                             ta_xxxxyz_xx_0, \
                             ta_xxxxyz_xy_0, \
                             ta_xxxxyz_xz_0, \
                             ta_xxxxyz_yy_0, \
                             ta_xxxxyz_yz_0, \
                             ta_xxxxyz_zz_0, \
                             ta_xxxxz_xx_0,  \
                             ta_xxxxz_xx_1,  \
                             ta_xxxxz_xz_0,  \
                             ta_xxxxz_xz_1,  \
                             ta_xxxxz_zz_0,  \
                             ta_xxxxz_zz_1,  \
                             ta_xxxyz_yz_0,  \
                             ta_xxxyz_yz_1,  \
                             ta_xxyz_yz_0,   \
                             ta_xxyz_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxyz_xx_0[i] = ta_xxxxz_xx_0[i] * pa_y[i] - ta_xxxxz_xx_1[i] * pc_y[i];

        ta_xxxxyz_xy_0[i] = ta_xxxxy_xy_0[i] * pa_z[i] - ta_xxxxy_xy_1[i] * pc_z[i];

        ta_xxxxyz_xz_0[i] = ta_xxxxz_xz_0[i] * pa_y[i] - ta_xxxxz_xz_1[i] * pc_y[i];

        ta_xxxxyz_yy_0[i] = ta_xxxxy_yy_0[i] * pa_z[i] - ta_xxxxy_yy_1[i] * pc_z[i];

        ta_xxxxyz_yz_0[i] = 3.0 * ta_xxyz_yz_0[i] * fe_0 - 3.0 * ta_xxyz_yz_1[i] * fe_0 + ta_xxxyz_yz_0[i] * pa_x[i] - ta_xxxyz_yz_1[i] * pc_x[i];

        ta_xxxxyz_zz_0[i] = ta_xxxxz_zz_0[i] * pa_y[i] - ta_xxxxz_zz_1[i] * pc_y[i];
    }

    // Set up 30-36 components of targeted buffer : ID

    auto ta_xxxxzz_xx_0 = pbuffer.data(idx_npot_0_id + 30);

    auto ta_xxxxzz_xy_0 = pbuffer.data(idx_npot_0_id + 31);

    auto ta_xxxxzz_xz_0 = pbuffer.data(idx_npot_0_id + 32);

    auto ta_xxxxzz_yy_0 = pbuffer.data(idx_npot_0_id + 33);

    auto ta_xxxxzz_yz_0 = pbuffer.data(idx_npot_0_id + 34);

    auto ta_xxxxzz_zz_0 = pbuffer.data(idx_npot_0_id + 35);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta_xxxx_xx_0,   \
                             ta_xxxx_xx_1,   \
                             ta_xxxx_xy_0,   \
                             ta_xxxx_xy_1,   \
                             ta_xxxxz_xx_0,  \
                             ta_xxxxz_xx_1,  \
                             ta_xxxxz_xy_0,  \
                             ta_xxxxz_xy_1,  \
                             ta_xxxxzz_xx_0, \
                             ta_xxxxzz_xy_0, \
                             ta_xxxxzz_xz_0, \
                             ta_xxxxzz_yy_0, \
                             ta_xxxxzz_yz_0, \
                             ta_xxxxzz_zz_0, \
                             ta_xxxzz_xz_0,  \
                             ta_xxxzz_xz_1,  \
                             ta_xxxzz_yy_0,  \
                             ta_xxxzz_yy_1,  \
                             ta_xxxzz_yz_0,  \
                             ta_xxxzz_yz_1,  \
                             ta_xxxzz_z_0,   \
                             ta_xxxzz_z_1,   \
                             ta_xxxzz_zz_0,  \
                             ta_xxxzz_zz_1,  \
                             ta_xxzz_xz_0,   \
                             ta_xxzz_xz_1,   \
                             ta_xxzz_yy_0,   \
                             ta_xxzz_yy_1,   \
                             ta_xxzz_yz_0,   \
                             ta_xxzz_yz_1,   \
                             ta_xxzz_zz_0,   \
                             ta_xxzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxxzz_xx_0[i] = ta_xxxx_xx_0[i] * fe_0 - ta_xxxx_xx_1[i] * fe_0 + ta_xxxxz_xx_0[i] * pa_z[i] - ta_xxxxz_xx_1[i] * pc_z[i];

        ta_xxxxzz_xy_0[i] = ta_xxxx_xy_0[i] * fe_0 - ta_xxxx_xy_1[i] * fe_0 + ta_xxxxz_xy_0[i] * pa_z[i] - ta_xxxxz_xy_1[i] * pc_z[i];

        ta_xxxxzz_xz_0[i] = 3.0 * ta_xxzz_xz_0[i] * fe_0 - 3.0 * ta_xxzz_xz_1[i] * fe_0 + ta_xxxzz_z_0[i] * fe_0 - ta_xxxzz_z_1[i] * fe_0 +
                            ta_xxxzz_xz_0[i] * pa_x[i] - ta_xxxzz_xz_1[i] * pc_x[i];

        ta_xxxxzz_yy_0[i] = 3.0 * ta_xxzz_yy_0[i] * fe_0 - 3.0 * ta_xxzz_yy_1[i] * fe_0 + ta_xxxzz_yy_0[i] * pa_x[i] - ta_xxxzz_yy_1[i] * pc_x[i];

        ta_xxxxzz_yz_0[i] = 3.0 * ta_xxzz_yz_0[i] * fe_0 - 3.0 * ta_xxzz_yz_1[i] * fe_0 + ta_xxxzz_yz_0[i] * pa_x[i] - ta_xxxzz_yz_1[i] * pc_x[i];

        ta_xxxxzz_zz_0[i] = 3.0 * ta_xxzz_zz_0[i] * fe_0 - 3.0 * ta_xxzz_zz_1[i] * fe_0 + ta_xxxzz_zz_0[i] * pa_x[i] - ta_xxxzz_zz_1[i] * pc_x[i];
    }

    // Set up 36-42 components of targeted buffer : ID

    auto ta_xxxyyy_xx_0 = pbuffer.data(idx_npot_0_id + 36);

    auto ta_xxxyyy_xy_0 = pbuffer.data(idx_npot_0_id + 37);

    auto ta_xxxyyy_xz_0 = pbuffer.data(idx_npot_0_id + 38);

    auto ta_xxxyyy_yy_0 = pbuffer.data(idx_npot_0_id + 39);

    auto ta_xxxyyy_yz_0 = pbuffer.data(idx_npot_0_id + 40);

    auto ta_xxxyyy_zz_0 = pbuffer.data(idx_npot_0_id + 41);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta_xxxy_xx_0,   \
                             ta_xxxy_xx_1,   \
                             ta_xxxy_xz_0,   \
                             ta_xxxy_xz_1,   \
                             ta_xxxyy_xx_0,  \
                             ta_xxxyy_xx_1,  \
                             ta_xxxyy_xz_0,  \
                             ta_xxxyy_xz_1,  \
                             ta_xxxyyy_xx_0, \
                             ta_xxxyyy_xy_0, \
                             ta_xxxyyy_xz_0, \
                             ta_xxxyyy_yy_0, \
                             ta_xxxyyy_yz_0, \
                             ta_xxxyyy_zz_0, \
                             ta_xxyyy_xy_0,  \
                             ta_xxyyy_xy_1,  \
                             ta_xxyyy_y_0,   \
                             ta_xxyyy_y_1,   \
                             ta_xxyyy_yy_0,  \
                             ta_xxyyy_yy_1,  \
                             ta_xxyyy_yz_0,  \
                             ta_xxyyy_yz_1,  \
                             ta_xxyyy_zz_0,  \
                             ta_xxyyy_zz_1,  \
                             ta_xyyy_xy_0,   \
                             ta_xyyy_xy_1,   \
                             ta_xyyy_yy_0,   \
                             ta_xyyy_yy_1,   \
                             ta_xyyy_yz_0,   \
                             ta_xyyy_yz_1,   \
                             ta_xyyy_zz_0,   \
                             ta_xyyy_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyyy_xx_0[i] = 2.0 * ta_xxxy_xx_0[i] * fe_0 - 2.0 * ta_xxxy_xx_1[i] * fe_0 + ta_xxxyy_xx_0[i] * pa_y[i] - ta_xxxyy_xx_1[i] * pc_y[i];

        ta_xxxyyy_xy_0[i] = 2.0 * ta_xyyy_xy_0[i] * fe_0 - 2.0 * ta_xyyy_xy_1[i] * fe_0 + ta_xxyyy_y_0[i] * fe_0 - ta_xxyyy_y_1[i] * fe_0 +
                            ta_xxyyy_xy_0[i] * pa_x[i] - ta_xxyyy_xy_1[i] * pc_x[i];

        ta_xxxyyy_xz_0[i] = 2.0 * ta_xxxy_xz_0[i] * fe_0 - 2.0 * ta_xxxy_xz_1[i] * fe_0 + ta_xxxyy_xz_0[i] * pa_y[i] - ta_xxxyy_xz_1[i] * pc_y[i];

        ta_xxxyyy_yy_0[i] = 2.0 * ta_xyyy_yy_0[i] * fe_0 - 2.0 * ta_xyyy_yy_1[i] * fe_0 + ta_xxyyy_yy_0[i] * pa_x[i] - ta_xxyyy_yy_1[i] * pc_x[i];

        ta_xxxyyy_yz_0[i] = 2.0 * ta_xyyy_yz_0[i] * fe_0 - 2.0 * ta_xyyy_yz_1[i] * fe_0 + ta_xxyyy_yz_0[i] * pa_x[i] - ta_xxyyy_yz_1[i] * pc_x[i];

        ta_xxxyyy_zz_0[i] = 2.0 * ta_xyyy_zz_0[i] * fe_0 - 2.0 * ta_xyyy_zz_1[i] * fe_0 + ta_xxyyy_zz_0[i] * pa_x[i] - ta_xxyyy_zz_1[i] * pc_x[i];
    }

    // Set up 42-48 components of targeted buffer : ID

    auto ta_xxxyyz_xx_0 = pbuffer.data(idx_npot_0_id + 42);

    auto ta_xxxyyz_xy_0 = pbuffer.data(idx_npot_0_id + 43);

    auto ta_xxxyyz_xz_0 = pbuffer.data(idx_npot_0_id + 44);

    auto ta_xxxyyz_yy_0 = pbuffer.data(idx_npot_0_id + 45);

    auto ta_xxxyyz_yz_0 = pbuffer.data(idx_npot_0_id + 46);

    auto ta_xxxyyz_zz_0 = pbuffer.data(idx_npot_0_id + 47);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             pc_x,           \
                             pc_y,           \
                             pc_z,           \
                             ta_xxxyy_xx_0,  \
                             ta_xxxyy_xx_1,  \
                             ta_xxxyy_xy_0,  \
                             ta_xxxyy_xy_1,  \
                             ta_xxxyy_yy_0,  \
                             ta_xxxyy_yy_1,  \
                             ta_xxxyyz_xx_0, \
                             ta_xxxyyz_xy_0, \
                             ta_xxxyyz_xz_0, \
                             ta_xxxyyz_yy_0, \
                             ta_xxxyyz_yz_0, \
                             ta_xxxyyz_zz_0, \
                             ta_xxxyz_xz_0,  \
                             ta_xxxyz_xz_1,  \
                             ta_xxxz_xz_0,   \
                             ta_xxxz_xz_1,   \
                             ta_xxyyz_yz_0,  \
                             ta_xxyyz_yz_1,  \
                             ta_xxyyz_zz_0,  \
                             ta_xxyyz_zz_1,  \
                             ta_xyyz_yz_0,   \
                             ta_xyyz_yz_1,   \
                             ta_xyyz_zz_0,   \
                             ta_xyyz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyyz_xx_0[i] = ta_xxxyy_xx_0[i] * pa_z[i] - ta_xxxyy_xx_1[i] * pc_z[i];

        ta_xxxyyz_xy_0[i] = ta_xxxyy_xy_0[i] * pa_z[i] - ta_xxxyy_xy_1[i] * pc_z[i];

        ta_xxxyyz_xz_0[i] = ta_xxxz_xz_0[i] * fe_0 - ta_xxxz_xz_1[i] * fe_0 + ta_xxxyz_xz_0[i] * pa_y[i] - ta_xxxyz_xz_1[i] * pc_y[i];

        ta_xxxyyz_yy_0[i] = ta_xxxyy_yy_0[i] * pa_z[i] - ta_xxxyy_yy_1[i] * pc_z[i];

        ta_xxxyyz_yz_0[i] = 2.0 * ta_xyyz_yz_0[i] * fe_0 - 2.0 * ta_xyyz_yz_1[i] * fe_0 + ta_xxyyz_yz_0[i] * pa_x[i] - ta_xxyyz_yz_1[i] * pc_x[i];

        ta_xxxyyz_zz_0[i] = 2.0 * ta_xyyz_zz_0[i] * fe_0 - 2.0 * ta_xyyz_zz_1[i] * fe_0 + ta_xxyyz_zz_0[i] * pa_x[i] - ta_xxyyz_zz_1[i] * pc_x[i];
    }

    // Set up 48-54 components of targeted buffer : ID

    auto ta_xxxyzz_xx_0 = pbuffer.data(idx_npot_0_id + 48);

    auto ta_xxxyzz_xy_0 = pbuffer.data(idx_npot_0_id + 49);

    auto ta_xxxyzz_xz_0 = pbuffer.data(idx_npot_0_id + 50);

    auto ta_xxxyzz_yy_0 = pbuffer.data(idx_npot_0_id + 51);

    auto ta_xxxyzz_yz_0 = pbuffer.data(idx_npot_0_id + 52);

    auto ta_xxxyzz_zz_0 = pbuffer.data(idx_npot_0_id + 53);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta_xxxyzz_xx_0, \
                             ta_xxxyzz_xy_0, \
                             ta_xxxyzz_xz_0, \
                             ta_xxxyzz_yy_0, \
                             ta_xxxyzz_yz_0, \
                             ta_xxxyzz_zz_0, \
                             ta_xxxzz_x_0,   \
                             ta_xxxzz_x_1,   \
                             ta_xxxzz_xx_0,  \
                             ta_xxxzz_xx_1,  \
                             ta_xxxzz_xy_0,  \
                             ta_xxxzz_xy_1,  \
                             ta_xxxzz_xz_0,  \
                             ta_xxxzz_xz_1,  \
                             ta_xxxzz_zz_0,  \
                             ta_xxxzz_zz_1,  \
                             ta_xxyzz_yy_0,  \
                             ta_xxyzz_yy_1,  \
                             ta_xxyzz_yz_0,  \
                             ta_xxyzz_yz_1,  \
                             ta_xyzz_yy_0,   \
                             ta_xyzz_yy_1,   \
                             ta_xyzz_yz_0,   \
                             ta_xyzz_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxyzz_xx_0[i] = ta_xxxzz_xx_0[i] * pa_y[i] - ta_xxxzz_xx_1[i] * pc_y[i];

        ta_xxxyzz_xy_0[i] = ta_xxxzz_x_0[i] * fe_0 - ta_xxxzz_x_1[i] * fe_0 + ta_xxxzz_xy_0[i] * pa_y[i] - ta_xxxzz_xy_1[i] * pc_y[i];

        ta_xxxyzz_xz_0[i] = ta_xxxzz_xz_0[i] * pa_y[i] - ta_xxxzz_xz_1[i] * pc_y[i];

        ta_xxxyzz_yy_0[i] = 2.0 * ta_xyzz_yy_0[i] * fe_0 - 2.0 * ta_xyzz_yy_1[i] * fe_0 + ta_xxyzz_yy_0[i] * pa_x[i] - ta_xxyzz_yy_1[i] * pc_x[i];

        ta_xxxyzz_yz_0[i] = 2.0 * ta_xyzz_yz_0[i] * fe_0 - 2.0 * ta_xyzz_yz_1[i] * fe_0 + ta_xxyzz_yz_0[i] * pa_x[i] - ta_xxyzz_yz_1[i] * pc_x[i];

        ta_xxxyzz_zz_0[i] = ta_xxxzz_zz_0[i] * pa_y[i] - ta_xxxzz_zz_1[i] * pc_y[i];
    }

    // Set up 54-60 components of targeted buffer : ID

    auto ta_xxxzzz_xx_0 = pbuffer.data(idx_npot_0_id + 54);

    auto ta_xxxzzz_xy_0 = pbuffer.data(idx_npot_0_id + 55);

    auto ta_xxxzzz_xz_0 = pbuffer.data(idx_npot_0_id + 56);

    auto ta_xxxzzz_yy_0 = pbuffer.data(idx_npot_0_id + 57);

    auto ta_xxxzzz_yz_0 = pbuffer.data(idx_npot_0_id + 58);

    auto ta_xxxzzz_zz_0 = pbuffer.data(idx_npot_0_id + 59);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta_xxxz_xx_0,   \
                             ta_xxxz_xx_1,   \
                             ta_xxxz_xy_0,   \
                             ta_xxxz_xy_1,   \
                             ta_xxxzz_xx_0,  \
                             ta_xxxzz_xx_1,  \
                             ta_xxxzz_xy_0,  \
                             ta_xxxzz_xy_1,  \
                             ta_xxxzzz_xx_0, \
                             ta_xxxzzz_xy_0, \
                             ta_xxxzzz_xz_0, \
                             ta_xxxzzz_yy_0, \
                             ta_xxxzzz_yz_0, \
                             ta_xxxzzz_zz_0, \
                             ta_xxzzz_xz_0,  \
                             ta_xxzzz_xz_1,  \
                             ta_xxzzz_yy_0,  \
                             ta_xxzzz_yy_1,  \
                             ta_xxzzz_yz_0,  \
                             ta_xxzzz_yz_1,  \
                             ta_xxzzz_z_0,   \
                             ta_xxzzz_z_1,   \
                             ta_xxzzz_zz_0,  \
                             ta_xxzzz_zz_1,  \
                             ta_xzzz_xz_0,   \
                             ta_xzzz_xz_1,   \
                             ta_xzzz_yy_0,   \
                             ta_xzzz_yy_1,   \
                             ta_xzzz_yz_0,   \
                             ta_xzzz_yz_1,   \
                             ta_xzzz_zz_0,   \
                             ta_xzzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxxzzz_xx_0[i] = 2.0 * ta_xxxz_xx_0[i] * fe_0 - 2.0 * ta_xxxz_xx_1[i] * fe_0 + ta_xxxzz_xx_0[i] * pa_z[i] - ta_xxxzz_xx_1[i] * pc_z[i];

        ta_xxxzzz_xy_0[i] = 2.0 * ta_xxxz_xy_0[i] * fe_0 - 2.0 * ta_xxxz_xy_1[i] * fe_0 + ta_xxxzz_xy_0[i] * pa_z[i] - ta_xxxzz_xy_1[i] * pc_z[i];

        ta_xxxzzz_xz_0[i] = 2.0 * ta_xzzz_xz_0[i] * fe_0 - 2.0 * ta_xzzz_xz_1[i] * fe_0 + ta_xxzzz_z_0[i] * fe_0 - ta_xxzzz_z_1[i] * fe_0 +
                            ta_xxzzz_xz_0[i] * pa_x[i] - ta_xxzzz_xz_1[i] * pc_x[i];

        ta_xxxzzz_yy_0[i] = 2.0 * ta_xzzz_yy_0[i] * fe_0 - 2.0 * ta_xzzz_yy_1[i] * fe_0 + ta_xxzzz_yy_0[i] * pa_x[i] - ta_xxzzz_yy_1[i] * pc_x[i];

        ta_xxxzzz_yz_0[i] = 2.0 * ta_xzzz_yz_0[i] * fe_0 - 2.0 * ta_xzzz_yz_1[i] * fe_0 + ta_xxzzz_yz_0[i] * pa_x[i] - ta_xxzzz_yz_1[i] * pc_x[i];

        ta_xxxzzz_zz_0[i] = 2.0 * ta_xzzz_zz_0[i] * fe_0 - 2.0 * ta_xzzz_zz_1[i] * fe_0 + ta_xxzzz_zz_0[i] * pa_x[i] - ta_xxzzz_zz_1[i] * pc_x[i];
    }

    // Set up 60-66 components of targeted buffer : ID

    auto ta_xxyyyy_xx_0 = pbuffer.data(idx_npot_0_id + 60);

    auto ta_xxyyyy_xy_0 = pbuffer.data(idx_npot_0_id + 61);

    auto ta_xxyyyy_xz_0 = pbuffer.data(idx_npot_0_id + 62);

    auto ta_xxyyyy_yy_0 = pbuffer.data(idx_npot_0_id + 63);

    auto ta_xxyyyy_yz_0 = pbuffer.data(idx_npot_0_id + 64);

    auto ta_xxyyyy_zz_0 = pbuffer.data(idx_npot_0_id + 65);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta_xxyy_xx_0,   \
                             ta_xxyy_xx_1,   \
                             ta_xxyy_xz_0,   \
                             ta_xxyy_xz_1,   \
                             ta_xxyyy_xx_0,  \
                             ta_xxyyy_xx_1,  \
                             ta_xxyyy_xz_0,  \
                             ta_xxyyy_xz_1,  \
                             ta_xxyyyy_xx_0, \
                             ta_xxyyyy_xy_0, \
                             ta_xxyyyy_xz_0, \
                             ta_xxyyyy_yy_0, \
                             ta_xxyyyy_yz_0, \
                             ta_xxyyyy_zz_0, \
                             ta_xyyyy_xy_0,  \
                             ta_xyyyy_xy_1,  \
                             ta_xyyyy_y_0,   \
                             ta_xyyyy_y_1,   \
                             ta_xyyyy_yy_0,  \
                             ta_xyyyy_yy_1,  \
                             ta_xyyyy_yz_0,  \
                             ta_xyyyy_yz_1,  \
                             ta_xyyyy_zz_0,  \
                             ta_xyyyy_zz_1,  \
                             ta_yyyy_xy_0,   \
                             ta_yyyy_xy_1,   \
                             ta_yyyy_yy_0,   \
                             ta_yyyy_yy_1,   \
                             ta_yyyy_yz_0,   \
                             ta_yyyy_yz_1,   \
                             ta_yyyy_zz_0,   \
                             ta_yyyy_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyyy_xx_0[i] = 3.0 * ta_xxyy_xx_0[i] * fe_0 - 3.0 * ta_xxyy_xx_1[i] * fe_0 + ta_xxyyy_xx_0[i] * pa_y[i] - ta_xxyyy_xx_1[i] * pc_y[i];

        ta_xxyyyy_xy_0[i] = ta_yyyy_xy_0[i] * fe_0 - ta_yyyy_xy_1[i] * fe_0 + ta_xyyyy_y_0[i] * fe_0 - ta_xyyyy_y_1[i] * fe_0 +
                            ta_xyyyy_xy_0[i] * pa_x[i] - ta_xyyyy_xy_1[i] * pc_x[i];

        ta_xxyyyy_xz_0[i] = 3.0 * ta_xxyy_xz_0[i] * fe_0 - 3.0 * ta_xxyy_xz_1[i] * fe_0 + ta_xxyyy_xz_0[i] * pa_y[i] - ta_xxyyy_xz_1[i] * pc_y[i];

        ta_xxyyyy_yy_0[i] = ta_yyyy_yy_0[i] * fe_0 - ta_yyyy_yy_1[i] * fe_0 + ta_xyyyy_yy_0[i] * pa_x[i] - ta_xyyyy_yy_1[i] * pc_x[i];

        ta_xxyyyy_yz_0[i] = ta_yyyy_yz_0[i] * fe_0 - ta_yyyy_yz_1[i] * fe_0 + ta_xyyyy_yz_0[i] * pa_x[i] - ta_xyyyy_yz_1[i] * pc_x[i];

        ta_xxyyyy_zz_0[i] = ta_yyyy_zz_0[i] * fe_0 - ta_yyyy_zz_1[i] * fe_0 + ta_xyyyy_zz_0[i] * pa_x[i] - ta_xyyyy_zz_1[i] * pc_x[i];
    }

    // Set up 66-72 components of targeted buffer : ID

    auto ta_xxyyyz_xx_0 = pbuffer.data(idx_npot_0_id + 66);

    auto ta_xxyyyz_xy_0 = pbuffer.data(idx_npot_0_id + 67);

    auto ta_xxyyyz_xz_0 = pbuffer.data(idx_npot_0_id + 68);

    auto ta_xxyyyz_yy_0 = pbuffer.data(idx_npot_0_id + 69);

    auto ta_xxyyyz_yz_0 = pbuffer.data(idx_npot_0_id + 70);

    auto ta_xxyyyz_zz_0 = pbuffer.data(idx_npot_0_id + 71);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             pc_x,           \
                             pc_y,           \
                             pc_z,           \
                             ta_xxyyy_xx_0,  \
                             ta_xxyyy_xx_1,  \
                             ta_xxyyy_xy_0,  \
                             ta_xxyyy_xy_1,  \
                             ta_xxyyy_yy_0,  \
                             ta_xxyyy_yy_1,  \
                             ta_xxyyyz_xx_0, \
                             ta_xxyyyz_xy_0, \
                             ta_xxyyyz_xz_0, \
                             ta_xxyyyz_yy_0, \
                             ta_xxyyyz_yz_0, \
                             ta_xxyyyz_zz_0, \
                             ta_xxyyz_xz_0,  \
                             ta_xxyyz_xz_1,  \
                             ta_xxyz_xz_0,   \
                             ta_xxyz_xz_1,   \
                             ta_xyyyz_yz_0,  \
                             ta_xyyyz_yz_1,  \
                             ta_xyyyz_zz_0,  \
                             ta_xyyyz_zz_1,  \
                             ta_yyyz_yz_0,   \
                             ta_yyyz_yz_1,   \
                             ta_yyyz_zz_0,   \
                             ta_yyyz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyyz_xx_0[i] = ta_xxyyy_xx_0[i] * pa_z[i] - ta_xxyyy_xx_1[i] * pc_z[i];

        ta_xxyyyz_xy_0[i] = ta_xxyyy_xy_0[i] * pa_z[i] - ta_xxyyy_xy_1[i] * pc_z[i];

        ta_xxyyyz_xz_0[i] = 2.0 * ta_xxyz_xz_0[i] * fe_0 - 2.0 * ta_xxyz_xz_1[i] * fe_0 + ta_xxyyz_xz_0[i] * pa_y[i] - ta_xxyyz_xz_1[i] * pc_y[i];

        ta_xxyyyz_yy_0[i] = ta_xxyyy_yy_0[i] * pa_z[i] - ta_xxyyy_yy_1[i] * pc_z[i];

        ta_xxyyyz_yz_0[i] = ta_yyyz_yz_0[i] * fe_0 - ta_yyyz_yz_1[i] * fe_0 + ta_xyyyz_yz_0[i] * pa_x[i] - ta_xyyyz_yz_1[i] * pc_x[i];

        ta_xxyyyz_zz_0[i] = ta_yyyz_zz_0[i] * fe_0 - ta_yyyz_zz_1[i] * fe_0 + ta_xyyyz_zz_0[i] * pa_x[i] - ta_xyyyz_zz_1[i] * pc_x[i];
    }

    // Set up 72-78 components of targeted buffer : ID

    auto ta_xxyyzz_xx_0 = pbuffer.data(idx_npot_0_id + 72);

    auto ta_xxyyzz_xy_0 = pbuffer.data(idx_npot_0_id + 73);

    auto ta_xxyyzz_xz_0 = pbuffer.data(idx_npot_0_id + 74);

    auto ta_xxyyzz_yy_0 = pbuffer.data(idx_npot_0_id + 75);

    auto ta_xxyyzz_yz_0 = pbuffer.data(idx_npot_0_id + 76);

    auto ta_xxyyzz_zz_0 = pbuffer.data(idx_npot_0_id + 77);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             pc_x,           \
                             pc_y,           \
                             pc_z,           \
                             ta_xxyy_xy_0,   \
                             ta_xxyy_xy_1,   \
                             ta_xxyyz_xy_0,  \
                             ta_xxyyz_xy_1,  \
                             ta_xxyyzz_xx_0, \
                             ta_xxyyzz_xy_0, \
                             ta_xxyyzz_xz_0, \
                             ta_xxyyzz_yy_0, \
                             ta_xxyyzz_yz_0, \
                             ta_xxyyzz_zz_0, \
                             ta_xxyzz_xx_0,  \
                             ta_xxyzz_xx_1,  \
                             ta_xxyzz_xz_0,  \
                             ta_xxyzz_xz_1,  \
                             ta_xxzz_xx_0,   \
                             ta_xxzz_xx_1,   \
                             ta_xxzz_xz_0,   \
                             ta_xxzz_xz_1,   \
                             ta_xyyzz_yy_0,  \
                             ta_xyyzz_yy_1,  \
                             ta_xyyzz_yz_0,  \
                             ta_xyyzz_yz_1,  \
                             ta_xyyzz_zz_0,  \
                             ta_xyyzz_zz_1,  \
                             ta_yyzz_yy_0,   \
                             ta_yyzz_yy_1,   \
                             ta_yyzz_yz_0,   \
                             ta_yyzz_yz_1,   \
                             ta_yyzz_zz_0,   \
                             ta_yyzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyyzz_xx_0[i] = ta_xxzz_xx_0[i] * fe_0 - ta_xxzz_xx_1[i] * fe_0 + ta_xxyzz_xx_0[i] * pa_y[i] - ta_xxyzz_xx_1[i] * pc_y[i];

        ta_xxyyzz_xy_0[i] = ta_xxyy_xy_0[i] * fe_0 - ta_xxyy_xy_1[i] * fe_0 + ta_xxyyz_xy_0[i] * pa_z[i] - ta_xxyyz_xy_1[i] * pc_z[i];

        ta_xxyyzz_xz_0[i] = ta_xxzz_xz_0[i] * fe_0 - ta_xxzz_xz_1[i] * fe_0 + ta_xxyzz_xz_0[i] * pa_y[i] - ta_xxyzz_xz_1[i] * pc_y[i];

        ta_xxyyzz_yy_0[i] = ta_yyzz_yy_0[i] * fe_0 - ta_yyzz_yy_1[i] * fe_0 + ta_xyyzz_yy_0[i] * pa_x[i] - ta_xyyzz_yy_1[i] * pc_x[i];

        ta_xxyyzz_yz_0[i] = ta_yyzz_yz_0[i] * fe_0 - ta_yyzz_yz_1[i] * fe_0 + ta_xyyzz_yz_0[i] * pa_x[i] - ta_xyyzz_yz_1[i] * pc_x[i];

        ta_xxyyzz_zz_0[i] = ta_yyzz_zz_0[i] * fe_0 - ta_yyzz_zz_1[i] * fe_0 + ta_xyyzz_zz_0[i] * pa_x[i] - ta_xyyzz_zz_1[i] * pc_x[i];
    }

    // Set up 78-84 components of targeted buffer : ID

    auto ta_xxyzzz_xx_0 = pbuffer.data(idx_npot_0_id + 78);

    auto ta_xxyzzz_xy_0 = pbuffer.data(idx_npot_0_id + 79);

    auto ta_xxyzzz_xz_0 = pbuffer.data(idx_npot_0_id + 80);

    auto ta_xxyzzz_yy_0 = pbuffer.data(idx_npot_0_id + 81);

    auto ta_xxyzzz_yz_0 = pbuffer.data(idx_npot_0_id + 82);

    auto ta_xxyzzz_zz_0 = pbuffer.data(idx_npot_0_id + 83);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta_xxyzzz_xx_0, \
                             ta_xxyzzz_xy_0, \
                             ta_xxyzzz_xz_0, \
                             ta_xxyzzz_yy_0, \
                             ta_xxyzzz_yz_0, \
                             ta_xxyzzz_zz_0, \
                             ta_xxzzz_x_0,   \
                             ta_xxzzz_x_1,   \
                             ta_xxzzz_xx_0,  \
                             ta_xxzzz_xx_1,  \
                             ta_xxzzz_xy_0,  \
                             ta_xxzzz_xy_1,  \
                             ta_xxzzz_xz_0,  \
                             ta_xxzzz_xz_1,  \
                             ta_xxzzz_zz_0,  \
                             ta_xxzzz_zz_1,  \
                             ta_xyzzz_yy_0,  \
                             ta_xyzzz_yy_1,  \
                             ta_xyzzz_yz_0,  \
                             ta_xyzzz_yz_1,  \
                             ta_yzzz_yy_0,   \
                             ta_yzzz_yy_1,   \
                             ta_yzzz_yz_0,   \
                             ta_yzzz_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxyzzz_xx_0[i] = ta_xxzzz_xx_0[i] * pa_y[i] - ta_xxzzz_xx_1[i] * pc_y[i];

        ta_xxyzzz_xy_0[i] = ta_xxzzz_x_0[i] * fe_0 - ta_xxzzz_x_1[i] * fe_0 + ta_xxzzz_xy_0[i] * pa_y[i] - ta_xxzzz_xy_1[i] * pc_y[i];

        ta_xxyzzz_xz_0[i] = ta_xxzzz_xz_0[i] * pa_y[i] - ta_xxzzz_xz_1[i] * pc_y[i];

        ta_xxyzzz_yy_0[i] = ta_yzzz_yy_0[i] * fe_0 - ta_yzzz_yy_1[i] * fe_0 + ta_xyzzz_yy_0[i] * pa_x[i] - ta_xyzzz_yy_1[i] * pc_x[i];

        ta_xxyzzz_yz_0[i] = ta_yzzz_yz_0[i] * fe_0 - ta_yzzz_yz_1[i] * fe_0 + ta_xyzzz_yz_0[i] * pa_x[i] - ta_xyzzz_yz_1[i] * pc_x[i];

        ta_xxyzzz_zz_0[i] = ta_xxzzz_zz_0[i] * pa_y[i] - ta_xxzzz_zz_1[i] * pc_y[i];
    }

    // Set up 84-90 components of targeted buffer : ID

    auto ta_xxzzzz_xx_0 = pbuffer.data(idx_npot_0_id + 84);

    auto ta_xxzzzz_xy_0 = pbuffer.data(idx_npot_0_id + 85);

    auto ta_xxzzzz_xz_0 = pbuffer.data(idx_npot_0_id + 86);

    auto ta_xxzzzz_yy_0 = pbuffer.data(idx_npot_0_id + 87);

    auto ta_xxzzzz_yz_0 = pbuffer.data(idx_npot_0_id + 88);

    auto ta_xxzzzz_zz_0 = pbuffer.data(idx_npot_0_id + 89);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta_xxzz_xx_0,   \
                             ta_xxzz_xx_1,   \
                             ta_xxzz_xy_0,   \
                             ta_xxzz_xy_1,   \
                             ta_xxzzz_xx_0,  \
                             ta_xxzzz_xx_1,  \
                             ta_xxzzz_xy_0,  \
                             ta_xxzzz_xy_1,  \
                             ta_xxzzzz_xx_0, \
                             ta_xxzzzz_xy_0, \
                             ta_xxzzzz_xz_0, \
                             ta_xxzzzz_yy_0, \
                             ta_xxzzzz_yz_0, \
                             ta_xxzzzz_zz_0, \
                             ta_xzzzz_xz_0,  \
                             ta_xzzzz_xz_1,  \
                             ta_xzzzz_yy_0,  \
                             ta_xzzzz_yy_1,  \
                             ta_xzzzz_yz_0,  \
                             ta_xzzzz_yz_1,  \
                             ta_xzzzz_z_0,   \
                             ta_xzzzz_z_1,   \
                             ta_xzzzz_zz_0,  \
                             ta_xzzzz_zz_1,  \
                             ta_zzzz_xz_0,   \
                             ta_zzzz_xz_1,   \
                             ta_zzzz_yy_0,   \
                             ta_zzzz_yy_1,   \
                             ta_zzzz_yz_0,   \
                             ta_zzzz_yz_1,   \
                             ta_zzzz_zz_0,   \
                             ta_zzzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxzzzz_xx_0[i] = 3.0 * ta_xxzz_xx_0[i] * fe_0 - 3.0 * ta_xxzz_xx_1[i] * fe_0 + ta_xxzzz_xx_0[i] * pa_z[i] - ta_xxzzz_xx_1[i] * pc_z[i];

        ta_xxzzzz_xy_0[i] = 3.0 * ta_xxzz_xy_0[i] * fe_0 - 3.0 * ta_xxzz_xy_1[i] * fe_0 + ta_xxzzz_xy_0[i] * pa_z[i] - ta_xxzzz_xy_1[i] * pc_z[i];

        ta_xxzzzz_xz_0[i] = ta_zzzz_xz_0[i] * fe_0 - ta_zzzz_xz_1[i] * fe_0 + ta_xzzzz_z_0[i] * fe_0 - ta_xzzzz_z_1[i] * fe_0 +
                            ta_xzzzz_xz_0[i] * pa_x[i] - ta_xzzzz_xz_1[i] * pc_x[i];

        ta_xxzzzz_yy_0[i] = ta_zzzz_yy_0[i] * fe_0 - ta_zzzz_yy_1[i] * fe_0 + ta_xzzzz_yy_0[i] * pa_x[i] - ta_xzzzz_yy_1[i] * pc_x[i];

        ta_xxzzzz_yz_0[i] = ta_zzzz_yz_0[i] * fe_0 - ta_zzzz_yz_1[i] * fe_0 + ta_xzzzz_yz_0[i] * pa_x[i] - ta_xzzzz_yz_1[i] * pc_x[i];

        ta_xxzzzz_zz_0[i] = ta_zzzz_zz_0[i] * fe_0 - ta_zzzz_zz_1[i] * fe_0 + ta_xzzzz_zz_0[i] * pa_x[i] - ta_xzzzz_zz_1[i] * pc_x[i];
    }

    // Set up 90-96 components of targeted buffer : ID

    auto ta_xyyyyy_xx_0 = pbuffer.data(idx_npot_0_id + 90);

    auto ta_xyyyyy_xy_0 = pbuffer.data(idx_npot_0_id + 91);

    auto ta_xyyyyy_xz_0 = pbuffer.data(idx_npot_0_id + 92);

    auto ta_xyyyyy_yy_0 = pbuffer.data(idx_npot_0_id + 93);

    auto ta_xyyyyy_yz_0 = pbuffer.data(idx_npot_0_id + 94);

    auto ta_xyyyyy_zz_0 = pbuffer.data(idx_npot_0_id + 95);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta_xyyyyy_xx_0, \
                             ta_xyyyyy_xy_0, \
                             ta_xyyyyy_xz_0, \
                             ta_xyyyyy_yy_0, \
                             ta_xyyyyy_yz_0, \
                             ta_xyyyyy_zz_0, \
                             ta_yyyyy_x_0,   \
                             ta_yyyyy_x_1,   \
                             ta_yyyyy_xx_0,  \
                             ta_yyyyy_xx_1,  \
                             ta_yyyyy_xy_0,  \
                             ta_yyyyy_xy_1,  \
                             ta_yyyyy_xz_0,  \
                             ta_yyyyy_xz_1,  \
                             ta_yyyyy_y_0,   \
                             ta_yyyyy_y_1,   \
                             ta_yyyyy_yy_0,  \
                             ta_yyyyy_yy_1,  \
                             ta_yyyyy_yz_0,  \
                             ta_yyyyy_yz_1,  \
                             ta_yyyyy_z_0,   \
                             ta_yyyyy_z_1,   \
                             ta_yyyyy_zz_0,  \
                             ta_yyyyy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyyy_xx_0[i] = 2.0 * ta_yyyyy_x_0[i] * fe_0 - 2.0 * ta_yyyyy_x_1[i] * fe_0 + ta_yyyyy_xx_0[i] * pa_x[i] - ta_yyyyy_xx_1[i] * pc_x[i];

        ta_xyyyyy_xy_0[i] = ta_yyyyy_y_0[i] * fe_0 - ta_yyyyy_y_1[i] * fe_0 + ta_yyyyy_xy_0[i] * pa_x[i] - ta_yyyyy_xy_1[i] * pc_x[i];

        ta_xyyyyy_xz_0[i] = ta_yyyyy_z_0[i] * fe_0 - ta_yyyyy_z_1[i] * fe_0 + ta_yyyyy_xz_0[i] * pa_x[i] - ta_yyyyy_xz_1[i] * pc_x[i];

        ta_xyyyyy_yy_0[i] = ta_yyyyy_yy_0[i] * pa_x[i] - ta_yyyyy_yy_1[i] * pc_x[i];

        ta_xyyyyy_yz_0[i] = ta_yyyyy_yz_0[i] * pa_x[i] - ta_yyyyy_yz_1[i] * pc_x[i];

        ta_xyyyyy_zz_0[i] = ta_yyyyy_zz_0[i] * pa_x[i] - ta_yyyyy_zz_1[i] * pc_x[i];
    }

    // Set up 96-102 components of targeted buffer : ID

    auto ta_xyyyyz_xx_0 = pbuffer.data(idx_npot_0_id + 96);

    auto ta_xyyyyz_xy_0 = pbuffer.data(idx_npot_0_id + 97);

    auto ta_xyyyyz_xz_0 = pbuffer.data(idx_npot_0_id + 98);

    auto ta_xyyyyz_yy_0 = pbuffer.data(idx_npot_0_id + 99);

    auto ta_xyyyyz_yz_0 = pbuffer.data(idx_npot_0_id + 100);

    auto ta_xyyyyz_zz_0 = pbuffer.data(idx_npot_0_id + 101);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             pc_x,           \
                             pc_z,           \
                             ta_xyyyy_xx_0,  \
                             ta_xyyyy_xx_1,  \
                             ta_xyyyy_xy_0,  \
                             ta_xyyyy_xy_1,  \
                             ta_xyyyyz_xx_0, \
                             ta_xyyyyz_xy_0, \
                             ta_xyyyyz_xz_0, \
                             ta_xyyyyz_yy_0, \
                             ta_xyyyyz_yz_0, \
                             ta_xyyyyz_zz_0, \
                             ta_yyyyz_xz_0,  \
                             ta_yyyyz_xz_1,  \
                             ta_yyyyz_yy_0,  \
                             ta_yyyyz_yy_1,  \
                             ta_yyyyz_yz_0,  \
                             ta_yyyyz_yz_1,  \
                             ta_yyyyz_z_0,   \
                             ta_yyyyz_z_1,   \
                             ta_yyyyz_zz_0,  \
                             ta_yyyyz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyyz_xx_0[i] = ta_xyyyy_xx_0[i] * pa_z[i] - ta_xyyyy_xx_1[i] * pc_z[i];

        ta_xyyyyz_xy_0[i] = ta_xyyyy_xy_0[i] * pa_z[i] - ta_xyyyy_xy_1[i] * pc_z[i];

        ta_xyyyyz_xz_0[i] = ta_yyyyz_z_0[i] * fe_0 - ta_yyyyz_z_1[i] * fe_0 + ta_yyyyz_xz_0[i] * pa_x[i] - ta_yyyyz_xz_1[i] * pc_x[i];

        ta_xyyyyz_yy_0[i] = ta_yyyyz_yy_0[i] * pa_x[i] - ta_yyyyz_yy_1[i] * pc_x[i];

        ta_xyyyyz_yz_0[i] = ta_yyyyz_yz_0[i] * pa_x[i] - ta_yyyyz_yz_1[i] * pc_x[i];

        ta_xyyyyz_zz_0[i] = ta_yyyyz_zz_0[i] * pa_x[i] - ta_yyyyz_zz_1[i] * pc_x[i];
    }

    // Set up 102-108 components of targeted buffer : ID

    auto ta_xyyyzz_xx_0 = pbuffer.data(idx_npot_0_id + 102);

    auto ta_xyyyzz_xy_0 = pbuffer.data(idx_npot_0_id + 103);

    auto ta_xyyyzz_xz_0 = pbuffer.data(idx_npot_0_id + 104);

    auto ta_xyyyzz_yy_0 = pbuffer.data(idx_npot_0_id + 105);

    auto ta_xyyyzz_yz_0 = pbuffer.data(idx_npot_0_id + 106);

    auto ta_xyyyzz_zz_0 = pbuffer.data(idx_npot_0_id + 107);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta_xyyyzz_xx_0, \
                             ta_xyyyzz_xy_0, \
                             ta_xyyyzz_xz_0, \
                             ta_xyyyzz_yy_0, \
                             ta_xyyyzz_yz_0, \
                             ta_xyyyzz_zz_0, \
                             ta_yyyzz_x_0,   \
                             ta_yyyzz_x_1,   \
                             ta_yyyzz_xx_0,  \
                             ta_yyyzz_xx_1,  \
                             ta_yyyzz_xy_0,  \
                             ta_yyyzz_xy_1,  \
                             ta_yyyzz_xz_0,  \
                             ta_yyyzz_xz_1,  \
                             ta_yyyzz_y_0,   \
                             ta_yyyzz_y_1,   \
                             ta_yyyzz_yy_0,  \
                             ta_yyyzz_yy_1,  \
                             ta_yyyzz_yz_0,  \
                             ta_yyyzz_yz_1,  \
                             ta_yyyzz_z_0,   \
                             ta_yyyzz_z_1,   \
                             ta_yyyzz_zz_0,  \
                             ta_yyyzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyyzz_xx_0[i] = 2.0 * ta_yyyzz_x_0[i] * fe_0 - 2.0 * ta_yyyzz_x_1[i] * fe_0 + ta_yyyzz_xx_0[i] * pa_x[i] - ta_yyyzz_xx_1[i] * pc_x[i];

        ta_xyyyzz_xy_0[i] = ta_yyyzz_y_0[i] * fe_0 - ta_yyyzz_y_1[i] * fe_0 + ta_yyyzz_xy_0[i] * pa_x[i] - ta_yyyzz_xy_1[i] * pc_x[i];

        ta_xyyyzz_xz_0[i] = ta_yyyzz_z_0[i] * fe_0 - ta_yyyzz_z_1[i] * fe_0 + ta_yyyzz_xz_0[i] * pa_x[i] - ta_yyyzz_xz_1[i] * pc_x[i];

        ta_xyyyzz_yy_0[i] = ta_yyyzz_yy_0[i] * pa_x[i] - ta_yyyzz_yy_1[i] * pc_x[i];

        ta_xyyyzz_yz_0[i] = ta_yyyzz_yz_0[i] * pa_x[i] - ta_yyyzz_yz_1[i] * pc_x[i];

        ta_xyyyzz_zz_0[i] = ta_yyyzz_zz_0[i] * pa_x[i] - ta_yyyzz_zz_1[i] * pc_x[i];
    }

    // Set up 108-114 components of targeted buffer : ID

    auto ta_xyyzzz_xx_0 = pbuffer.data(idx_npot_0_id + 108);

    auto ta_xyyzzz_xy_0 = pbuffer.data(idx_npot_0_id + 109);

    auto ta_xyyzzz_xz_0 = pbuffer.data(idx_npot_0_id + 110);

    auto ta_xyyzzz_yy_0 = pbuffer.data(idx_npot_0_id + 111);

    auto ta_xyyzzz_yz_0 = pbuffer.data(idx_npot_0_id + 112);

    auto ta_xyyzzz_zz_0 = pbuffer.data(idx_npot_0_id + 113);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta_xyyzzz_xx_0, \
                             ta_xyyzzz_xy_0, \
                             ta_xyyzzz_xz_0, \
                             ta_xyyzzz_yy_0, \
                             ta_xyyzzz_yz_0, \
                             ta_xyyzzz_zz_0, \
                             ta_yyzzz_x_0,   \
                             ta_yyzzz_x_1,   \
                             ta_yyzzz_xx_0,  \
                             ta_yyzzz_xx_1,  \
                             ta_yyzzz_xy_0,  \
                             ta_yyzzz_xy_1,  \
                             ta_yyzzz_xz_0,  \
                             ta_yyzzz_xz_1,  \
                             ta_yyzzz_y_0,   \
                             ta_yyzzz_y_1,   \
                             ta_yyzzz_yy_0,  \
                             ta_yyzzz_yy_1,  \
                             ta_yyzzz_yz_0,  \
                             ta_yyzzz_yz_1,  \
                             ta_yyzzz_z_0,   \
                             ta_yyzzz_z_1,   \
                             ta_yyzzz_zz_0,  \
                             ta_yyzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyyzzz_xx_0[i] = 2.0 * ta_yyzzz_x_0[i] * fe_0 - 2.0 * ta_yyzzz_x_1[i] * fe_0 + ta_yyzzz_xx_0[i] * pa_x[i] - ta_yyzzz_xx_1[i] * pc_x[i];

        ta_xyyzzz_xy_0[i] = ta_yyzzz_y_0[i] * fe_0 - ta_yyzzz_y_1[i] * fe_0 + ta_yyzzz_xy_0[i] * pa_x[i] - ta_yyzzz_xy_1[i] * pc_x[i];

        ta_xyyzzz_xz_0[i] = ta_yyzzz_z_0[i] * fe_0 - ta_yyzzz_z_1[i] * fe_0 + ta_yyzzz_xz_0[i] * pa_x[i] - ta_yyzzz_xz_1[i] * pc_x[i];

        ta_xyyzzz_yy_0[i] = ta_yyzzz_yy_0[i] * pa_x[i] - ta_yyzzz_yy_1[i] * pc_x[i];

        ta_xyyzzz_yz_0[i] = ta_yyzzz_yz_0[i] * pa_x[i] - ta_yyzzz_yz_1[i] * pc_x[i];

        ta_xyyzzz_zz_0[i] = ta_yyzzz_zz_0[i] * pa_x[i] - ta_yyzzz_zz_1[i] * pc_x[i];
    }

    // Set up 114-120 components of targeted buffer : ID

    auto ta_xyzzzz_xx_0 = pbuffer.data(idx_npot_0_id + 114);

    auto ta_xyzzzz_xy_0 = pbuffer.data(idx_npot_0_id + 115);

    auto ta_xyzzzz_xz_0 = pbuffer.data(idx_npot_0_id + 116);

    auto ta_xyzzzz_yy_0 = pbuffer.data(idx_npot_0_id + 117);

    auto ta_xyzzzz_yz_0 = pbuffer.data(idx_npot_0_id + 118);

    auto ta_xyzzzz_zz_0 = pbuffer.data(idx_npot_0_id + 119);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pc_x,           \
                             pc_y,           \
                             ta_xyzzzz_xx_0, \
                             ta_xyzzzz_xy_0, \
                             ta_xyzzzz_xz_0, \
                             ta_xyzzzz_yy_0, \
                             ta_xyzzzz_yz_0, \
                             ta_xyzzzz_zz_0, \
                             ta_xzzzz_xx_0,  \
                             ta_xzzzz_xx_1,  \
                             ta_xzzzz_xz_0,  \
                             ta_xzzzz_xz_1,  \
                             ta_yzzzz_xy_0,  \
                             ta_yzzzz_xy_1,  \
                             ta_yzzzz_y_0,   \
                             ta_yzzzz_y_1,   \
                             ta_yzzzz_yy_0,  \
                             ta_yzzzz_yy_1,  \
                             ta_yzzzz_yz_0,  \
                             ta_yzzzz_yz_1,  \
                             ta_yzzzz_zz_0,  \
                             ta_yzzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyzzzz_xx_0[i] = ta_xzzzz_xx_0[i] * pa_y[i] - ta_xzzzz_xx_1[i] * pc_y[i];

        ta_xyzzzz_xy_0[i] = ta_yzzzz_y_0[i] * fe_0 - ta_yzzzz_y_1[i] * fe_0 + ta_yzzzz_xy_0[i] * pa_x[i] - ta_yzzzz_xy_1[i] * pc_x[i];

        ta_xyzzzz_xz_0[i] = ta_xzzzz_xz_0[i] * pa_y[i] - ta_xzzzz_xz_1[i] * pc_y[i];

        ta_xyzzzz_yy_0[i] = ta_yzzzz_yy_0[i] * pa_x[i] - ta_yzzzz_yy_1[i] * pc_x[i];

        ta_xyzzzz_yz_0[i] = ta_yzzzz_yz_0[i] * pa_x[i] - ta_yzzzz_yz_1[i] * pc_x[i];

        ta_xyzzzz_zz_0[i] = ta_yzzzz_zz_0[i] * pa_x[i] - ta_yzzzz_zz_1[i] * pc_x[i];
    }

    // Set up 120-126 components of targeted buffer : ID

    auto ta_xzzzzz_xx_0 = pbuffer.data(idx_npot_0_id + 120);

    auto ta_xzzzzz_xy_0 = pbuffer.data(idx_npot_0_id + 121);

    auto ta_xzzzzz_xz_0 = pbuffer.data(idx_npot_0_id + 122);

    auto ta_xzzzzz_yy_0 = pbuffer.data(idx_npot_0_id + 123);

    auto ta_xzzzzz_yz_0 = pbuffer.data(idx_npot_0_id + 124);

    auto ta_xzzzzz_zz_0 = pbuffer.data(idx_npot_0_id + 125);

#pragma omp simd aligned(pa_x,               \
                             pc_x,           \
                             ta_xzzzzz_xx_0, \
                             ta_xzzzzz_xy_0, \
                             ta_xzzzzz_xz_0, \
                             ta_xzzzzz_yy_0, \
                             ta_xzzzzz_yz_0, \
                             ta_xzzzzz_zz_0, \
                             ta_zzzzz_x_0,   \
                             ta_zzzzz_x_1,   \
                             ta_zzzzz_xx_0,  \
                             ta_zzzzz_xx_1,  \
                             ta_zzzzz_xy_0,  \
                             ta_zzzzz_xy_1,  \
                             ta_zzzzz_xz_0,  \
                             ta_zzzzz_xz_1,  \
                             ta_zzzzz_y_0,   \
                             ta_zzzzz_y_1,   \
                             ta_zzzzz_yy_0,  \
                             ta_zzzzz_yy_1,  \
                             ta_zzzzz_yz_0,  \
                             ta_zzzzz_yz_1,  \
                             ta_zzzzz_z_0,   \
                             ta_zzzzz_z_1,   \
                             ta_zzzzz_zz_0,  \
                             ta_zzzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzzzzz_xx_0[i] = 2.0 * ta_zzzzz_x_0[i] * fe_0 - 2.0 * ta_zzzzz_x_1[i] * fe_0 + ta_zzzzz_xx_0[i] * pa_x[i] - ta_zzzzz_xx_1[i] * pc_x[i];

        ta_xzzzzz_xy_0[i] = ta_zzzzz_y_0[i] * fe_0 - ta_zzzzz_y_1[i] * fe_0 + ta_zzzzz_xy_0[i] * pa_x[i] - ta_zzzzz_xy_1[i] * pc_x[i];

        ta_xzzzzz_xz_0[i] = ta_zzzzz_z_0[i] * fe_0 - ta_zzzzz_z_1[i] * fe_0 + ta_zzzzz_xz_0[i] * pa_x[i] - ta_zzzzz_xz_1[i] * pc_x[i];

        ta_xzzzzz_yy_0[i] = ta_zzzzz_yy_0[i] * pa_x[i] - ta_zzzzz_yy_1[i] * pc_x[i];

        ta_xzzzzz_yz_0[i] = ta_zzzzz_yz_0[i] * pa_x[i] - ta_zzzzz_yz_1[i] * pc_x[i];

        ta_xzzzzz_zz_0[i] = ta_zzzzz_zz_0[i] * pa_x[i] - ta_zzzzz_zz_1[i] * pc_x[i];
    }

    // Set up 126-132 components of targeted buffer : ID

    auto ta_yyyyyy_xx_0 = pbuffer.data(idx_npot_0_id + 126);

    auto ta_yyyyyy_xy_0 = pbuffer.data(idx_npot_0_id + 127);

    auto ta_yyyyyy_xz_0 = pbuffer.data(idx_npot_0_id + 128);

    auto ta_yyyyyy_yy_0 = pbuffer.data(idx_npot_0_id + 129);

    auto ta_yyyyyy_yz_0 = pbuffer.data(idx_npot_0_id + 130);

    auto ta_yyyyyy_zz_0 = pbuffer.data(idx_npot_0_id + 131);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta_yyyy_xx_0,   \
                             ta_yyyy_xx_1,   \
                             ta_yyyy_xy_0,   \
                             ta_yyyy_xy_1,   \
                             ta_yyyy_xz_0,   \
                             ta_yyyy_xz_1,   \
                             ta_yyyy_yy_0,   \
                             ta_yyyy_yy_1,   \
                             ta_yyyy_yz_0,   \
                             ta_yyyy_yz_1,   \
                             ta_yyyy_zz_0,   \
                             ta_yyyy_zz_1,   \
                             ta_yyyyy_x_0,   \
                             ta_yyyyy_x_1,   \
                             ta_yyyyy_xx_0,  \
                             ta_yyyyy_xx_1,  \
                             ta_yyyyy_xy_0,  \
                             ta_yyyyy_xy_1,  \
                             ta_yyyyy_xz_0,  \
                             ta_yyyyy_xz_1,  \
                             ta_yyyyy_y_0,   \
                             ta_yyyyy_y_1,   \
                             ta_yyyyy_yy_0,  \
                             ta_yyyyy_yy_1,  \
                             ta_yyyyy_yz_0,  \
                             ta_yyyyy_yz_1,  \
                             ta_yyyyy_z_0,   \
                             ta_yyyyy_z_1,   \
                             ta_yyyyy_zz_0,  \
                             ta_yyyyy_zz_1,  \
                             ta_yyyyyy_xx_0, \
                             ta_yyyyyy_xy_0, \
                             ta_yyyyyy_xz_0, \
                             ta_yyyyyy_yy_0, \
                             ta_yyyyyy_yz_0, \
                             ta_yyyyyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyyy_xx_0[i] = 5.0 * ta_yyyy_xx_0[i] * fe_0 - 5.0 * ta_yyyy_xx_1[i] * fe_0 + ta_yyyyy_xx_0[i] * pa_y[i] - ta_yyyyy_xx_1[i] * pc_y[i];

        ta_yyyyyy_xy_0[i] = 5.0 * ta_yyyy_xy_0[i] * fe_0 - 5.0 * ta_yyyy_xy_1[i] * fe_0 + ta_yyyyy_x_0[i] * fe_0 - ta_yyyyy_x_1[i] * fe_0 +
                            ta_yyyyy_xy_0[i] * pa_y[i] - ta_yyyyy_xy_1[i] * pc_y[i];

        ta_yyyyyy_xz_0[i] = 5.0 * ta_yyyy_xz_0[i] * fe_0 - 5.0 * ta_yyyy_xz_1[i] * fe_0 + ta_yyyyy_xz_0[i] * pa_y[i] - ta_yyyyy_xz_1[i] * pc_y[i];

        ta_yyyyyy_yy_0[i] = 5.0 * ta_yyyy_yy_0[i] * fe_0 - 5.0 * ta_yyyy_yy_1[i] * fe_0 + 2.0 * ta_yyyyy_y_0[i] * fe_0 -
                            2.0 * ta_yyyyy_y_1[i] * fe_0 + ta_yyyyy_yy_0[i] * pa_y[i] - ta_yyyyy_yy_1[i] * pc_y[i];

        ta_yyyyyy_yz_0[i] = 5.0 * ta_yyyy_yz_0[i] * fe_0 - 5.0 * ta_yyyy_yz_1[i] * fe_0 + ta_yyyyy_z_0[i] * fe_0 - ta_yyyyy_z_1[i] * fe_0 +
                            ta_yyyyy_yz_0[i] * pa_y[i] - ta_yyyyy_yz_1[i] * pc_y[i];

        ta_yyyyyy_zz_0[i] = 5.0 * ta_yyyy_zz_0[i] * fe_0 - 5.0 * ta_yyyy_zz_1[i] * fe_0 + ta_yyyyy_zz_0[i] * pa_y[i] - ta_yyyyy_zz_1[i] * pc_y[i];
    }

    // Set up 132-138 components of targeted buffer : ID

    auto ta_yyyyyz_xx_0 = pbuffer.data(idx_npot_0_id + 132);

    auto ta_yyyyyz_xy_0 = pbuffer.data(idx_npot_0_id + 133);

    auto ta_yyyyyz_xz_0 = pbuffer.data(idx_npot_0_id + 134);

    auto ta_yyyyyz_yy_0 = pbuffer.data(idx_npot_0_id + 135);

    auto ta_yyyyyz_yz_0 = pbuffer.data(idx_npot_0_id + 136);

    auto ta_yyyyyz_zz_0 = pbuffer.data(idx_npot_0_id + 137);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta_yyyyy_xx_0,  \
                             ta_yyyyy_xx_1,  \
                             ta_yyyyy_xy_0,  \
                             ta_yyyyy_xy_1,  \
                             ta_yyyyy_y_0,   \
                             ta_yyyyy_y_1,   \
                             ta_yyyyy_yy_0,  \
                             ta_yyyyy_yy_1,  \
                             ta_yyyyy_yz_0,  \
                             ta_yyyyy_yz_1,  \
                             ta_yyyyyz_xx_0, \
                             ta_yyyyyz_xy_0, \
                             ta_yyyyyz_xz_0, \
                             ta_yyyyyz_yy_0, \
                             ta_yyyyyz_yz_0, \
                             ta_yyyyyz_zz_0, \
                             ta_yyyyz_xz_0,  \
                             ta_yyyyz_xz_1,  \
                             ta_yyyyz_zz_0,  \
                             ta_yyyyz_zz_1,  \
                             ta_yyyz_xz_0,   \
                             ta_yyyz_xz_1,   \
                             ta_yyyz_zz_0,   \
                             ta_yyyz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyyz_xx_0[i] = ta_yyyyy_xx_0[i] * pa_z[i] - ta_yyyyy_xx_1[i] * pc_z[i];

        ta_yyyyyz_xy_0[i] = ta_yyyyy_xy_0[i] * pa_z[i] - ta_yyyyy_xy_1[i] * pc_z[i];

        ta_yyyyyz_xz_0[i] = 4.0 * ta_yyyz_xz_0[i] * fe_0 - 4.0 * ta_yyyz_xz_1[i] * fe_0 + ta_yyyyz_xz_0[i] * pa_y[i] - ta_yyyyz_xz_1[i] * pc_y[i];

        ta_yyyyyz_yy_0[i] = ta_yyyyy_yy_0[i] * pa_z[i] - ta_yyyyy_yy_1[i] * pc_z[i];

        ta_yyyyyz_yz_0[i] = ta_yyyyy_y_0[i] * fe_0 - ta_yyyyy_y_1[i] * fe_0 + ta_yyyyy_yz_0[i] * pa_z[i] - ta_yyyyy_yz_1[i] * pc_z[i];

        ta_yyyyyz_zz_0[i] = 4.0 * ta_yyyz_zz_0[i] * fe_0 - 4.0 * ta_yyyz_zz_1[i] * fe_0 + ta_yyyyz_zz_0[i] * pa_y[i] - ta_yyyyz_zz_1[i] * pc_y[i];
    }

    // Set up 138-144 components of targeted buffer : ID

    auto ta_yyyyzz_xx_0 = pbuffer.data(idx_npot_0_id + 138);

    auto ta_yyyyzz_xy_0 = pbuffer.data(idx_npot_0_id + 139);

    auto ta_yyyyzz_xz_0 = pbuffer.data(idx_npot_0_id + 140);

    auto ta_yyyyzz_yy_0 = pbuffer.data(idx_npot_0_id + 141);

    auto ta_yyyyzz_yz_0 = pbuffer.data(idx_npot_0_id + 142);

    auto ta_yyyyzz_zz_0 = pbuffer.data(idx_npot_0_id + 143);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta_yyyy_xy_0,   \
                             ta_yyyy_xy_1,   \
                             ta_yyyy_yy_0,   \
                             ta_yyyy_yy_1,   \
                             ta_yyyyz_xy_0,  \
                             ta_yyyyz_xy_1,  \
                             ta_yyyyz_yy_0,  \
                             ta_yyyyz_yy_1,  \
                             ta_yyyyzz_xx_0, \
                             ta_yyyyzz_xy_0, \
                             ta_yyyyzz_xz_0, \
                             ta_yyyyzz_yy_0, \
                             ta_yyyyzz_yz_0, \
                             ta_yyyyzz_zz_0, \
                             ta_yyyzz_xx_0,  \
                             ta_yyyzz_xx_1,  \
                             ta_yyyzz_xz_0,  \
                             ta_yyyzz_xz_1,  \
                             ta_yyyzz_yz_0,  \
                             ta_yyyzz_yz_1,  \
                             ta_yyyzz_z_0,   \
                             ta_yyyzz_z_1,   \
                             ta_yyyzz_zz_0,  \
                             ta_yyyzz_zz_1,  \
                             ta_yyzz_xx_0,   \
                             ta_yyzz_xx_1,   \
                             ta_yyzz_xz_0,   \
                             ta_yyzz_xz_1,   \
                             ta_yyzz_yz_0,   \
                             ta_yyzz_yz_1,   \
                             ta_yyzz_zz_0,   \
                             ta_yyzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyyzz_xx_0[i] = 3.0 * ta_yyzz_xx_0[i] * fe_0 - 3.0 * ta_yyzz_xx_1[i] * fe_0 + ta_yyyzz_xx_0[i] * pa_y[i] - ta_yyyzz_xx_1[i] * pc_y[i];

        ta_yyyyzz_xy_0[i] = ta_yyyy_xy_0[i] * fe_0 - ta_yyyy_xy_1[i] * fe_0 + ta_yyyyz_xy_0[i] * pa_z[i] - ta_yyyyz_xy_1[i] * pc_z[i];

        ta_yyyyzz_xz_0[i] = 3.0 * ta_yyzz_xz_0[i] * fe_0 - 3.0 * ta_yyzz_xz_1[i] * fe_0 + ta_yyyzz_xz_0[i] * pa_y[i] - ta_yyyzz_xz_1[i] * pc_y[i];

        ta_yyyyzz_yy_0[i] = ta_yyyy_yy_0[i] * fe_0 - ta_yyyy_yy_1[i] * fe_0 + ta_yyyyz_yy_0[i] * pa_z[i] - ta_yyyyz_yy_1[i] * pc_z[i];

        ta_yyyyzz_yz_0[i] = 3.0 * ta_yyzz_yz_0[i] * fe_0 - 3.0 * ta_yyzz_yz_1[i] * fe_0 + ta_yyyzz_z_0[i] * fe_0 - ta_yyyzz_z_1[i] * fe_0 +
                            ta_yyyzz_yz_0[i] * pa_y[i] - ta_yyyzz_yz_1[i] * pc_y[i];

        ta_yyyyzz_zz_0[i] = 3.0 * ta_yyzz_zz_0[i] * fe_0 - 3.0 * ta_yyzz_zz_1[i] * fe_0 + ta_yyyzz_zz_0[i] * pa_y[i] - ta_yyyzz_zz_1[i] * pc_y[i];
    }

    // Set up 144-150 components of targeted buffer : ID

    auto ta_yyyzzz_xx_0 = pbuffer.data(idx_npot_0_id + 144);

    auto ta_yyyzzz_xy_0 = pbuffer.data(idx_npot_0_id + 145);

    auto ta_yyyzzz_xz_0 = pbuffer.data(idx_npot_0_id + 146);

    auto ta_yyyzzz_yy_0 = pbuffer.data(idx_npot_0_id + 147);

    auto ta_yyyzzz_yz_0 = pbuffer.data(idx_npot_0_id + 148);

    auto ta_yyyzzz_zz_0 = pbuffer.data(idx_npot_0_id + 149);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta_yyyz_xy_0,   \
                             ta_yyyz_xy_1,   \
                             ta_yyyz_yy_0,   \
                             ta_yyyz_yy_1,   \
                             ta_yyyzz_xy_0,  \
                             ta_yyyzz_xy_1,  \
                             ta_yyyzz_yy_0,  \
                             ta_yyyzz_yy_1,  \
                             ta_yyyzzz_xx_0, \
                             ta_yyyzzz_xy_0, \
                             ta_yyyzzz_xz_0, \
                             ta_yyyzzz_yy_0, \
                             ta_yyyzzz_yz_0, \
                             ta_yyyzzz_zz_0, \
                             ta_yyzzz_xx_0,  \
                             ta_yyzzz_xx_1,  \
                             ta_yyzzz_xz_0,  \
                             ta_yyzzz_xz_1,  \
                             ta_yyzzz_yz_0,  \
                             ta_yyzzz_yz_1,  \
                             ta_yyzzz_z_0,   \
                             ta_yyzzz_z_1,   \
                             ta_yyzzz_zz_0,  \
                             ta_yyzzz_zz_1,  \
                             ta_yzzz_xx_0,   \
                             ta_yzzz_xx_1,   \
                             ta_yzzz_xz_0,   \
                             ta_yzzz_xz_1,   \
                             ta_yzzz_yz_0,   \
                             ta_yzzz_yz_1,   \
                             ta_yzzz_zz_0,   \
                             ta_yzzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyyzzz_xx_0[i] = 2.0 * ta_yzzz_xx_0[i] * fe_0 - 2.0 * ta_yzzz_xx_1[i] * fe_0 + ta_yyzzz_xx_0[i] * pa_y[i] - ta_yyzzz_xx_1[i] * pc_y[i];

        ta_yyyzzz_xy_0[i] = 2.0 * ta_yyyz_xy_0[i] * fe_0 - 2.0 * ta_yyyz_xy_1[i] * fe_0 + ta_yyyzz_xy_0[i] * pa_z[i] - ta_yyyzz_xy_1[i] * pc_z[i];

        ta_yyyzzz_xz_0[i] = 2.0 * ta_yzzz_xz_0[i] * fe_0 - 2.0 * ta_yzzz_xz_1[i] * fe_0 + ta_yyzzz_xz_0[i] * pa_y[i] - ta_yyzzz_xz_1[i] * pc_y[i];

        ta_yyyzzz_yy_0[i] = 2.0 * ta_yyyz_yy_0[i] * fe_0 - 2.0 * ta_yyyz_yy_1[i] * fe_0 + ta_yyyzz_yy_0[i] * pa_z[i] - ta_yyyzz_yy_1[i] * pc_z[i];

        ta_yyyzzz_yz_0[i] = 2.0 * ta_yzzz_yz_0[i] * fe_0 - 2.0 * ta_yzzz_yz_1[i] * fe_0 + ta_yyzzz_z_0[i] * fe_0 - ta_yyzzz_z_1[i] * fe_0 +
                            ta_yyzzz_yz_0[i] * pa_y[i] - ta_yyzzz_yz_1[i] * pc_y[i];

        ta_yyyzzz_zz_0[i] = 2.0 * ta_yzzz_zz_0[i] * fe_0 - 2.0 * ta_yzzz_zz_1[i] * fe_0 + ta_yyzzz_zz_0[i] * pa_y[i] - ta_yyzzz_zz_1[i] * pc_y[i];
    }

    // Set up 150-156 components of targeted buffer : ID

    auto ta_yyzzzz_xx_0 = pbuffer.data(idx_npot_0_id + 150);

    auto ta_yyzzzz_xy_0 = pbuffer.data(idx_npot_0_id + 151);

    auto ta_yyzzzz_xz_0 = pbuffer.data(idx_npot_0_id + 152);

    auto ta_yyzzzz_yy_0 = pbuffer.data(idx_npot_0_id + 153);

    auto ta_yyzzzz_yz_0 = pbuffer.data(idx_npot_0_id + 154);

    auto ta_yyzzzz_zz_0 = pbuffer.data(idx_npot_0_id + 155);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             pc_y,           \
                             pc_z,           \
                             ta_yyzz_xy_0,   \
                             ta_yyzz_xy_1,   \
                             ta_yyzz_yy_0,   \
                             ta_yyzz_yy_1,   \
                             ta_yyzzz_xy_0,  \
                             ta_yyzzz_xy_1,  \
                             ta_yyzzz_yy_0,  \
                             ta_yyzzz_yy_1,  \
                             ta_yyzzzz_xx_0, \
                             ta_yyzzzz_xy_0, \
                             ta_yyzzzz_xz_0, \
                             ta_yyzzzz_yy_0, \
                             ta_yyzzzz_yz_0, \
                             ta_yyzzzz_zz_0, \
                             ta_yzzzz_xx_0,  \
                             ta_yzzzz_xx_1,  \
                             ta_yzzzz_xz_0,  \
                             ta_yzzzz_xz_1,  \
                             ta_yzzzz_yz_0,  \
                             ta_yzzzz_yz_1,  \
                             ta_yzzzz_z_0,   \
                             ta_yzzzz_z_1,   \
                             ta_yzzzz_zz_0,  \
                             ta_yzzzz_zz_1,  \
                             ta_zzzz_xx_0,   \
                             ta_zzzz_xx_1,   \
                             ta_zzzz_xz_0,   \
                             ta_zzzz_xz_1,   \
                             ta_zzzz_yz_0,   \
                             ta_zzzz_yz_1,   \
                             ta_zzzz_zz_0,   \
                             ta_zzzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyzzzz_xx_0[i] = ta_zzzz_xx_0[i] * fe_0 - ta_zzzz_xx_1[i] * fe_0 + ta_yzzzz_xx_0[i] * pa_y[i] - ta_yzzzz_xx_1[i] * pc_y[i];

        ta_yyzzzz_xy_0[i] = 3.0 * ta_yyzz_xy_0[i] * fe_0 - 3.0 * ta_yyzz_xy_1[i] * fe_0 + ta_yyzzz_xy_0[i] * pa_z[i] - ta_yyzzz_xy_1[i] * pc_z[i];

        ta_yyzzzz_xz_0[i] = ta_zzzz_xz_0[i] * fe_0 - ta_zzzz_xz_1[i] * fe_0 + ta_yzzzz_xz_0[i] * pa_y[i] - ta_yzzzz_xz_1[i] * pc_y[i];

        ta_yyzzzz_yy_0[i] = 3.0 * ta_yyzz_yy_0[i] * fe_0 - 3.0 * ta_yyzz_yy_1[i] * fe_0 + ta_yyzzz_yy_0[i] * pa_z[i] - ta_yyzzz_yy_1[i] * pc_z[i];

        ta_yyzzzz_yz_0[i] = ta_zzzz_yz_0[i] * fe_0 - ta_zzzz_yz_1[i] * fe_0 + ta_yzzzz_z_0[i] * fe_0 - ta_yzzzz_z_1[i] * fe_0 +
                            ta_yzzzz_yz_0[i] * pa_y[i] - ta_yzzzz_yz_1[i] * pc_y[i];

        ta_yyzzzz_zz_0[i] = ta_zzzz_zz_0[i] * fe_0 - ta_zzzz_zz_1[i] * fe_0 + ta_yzzzz_zz_0[i] * pa_y[i] - ta_yzzzz_zz_1[i] * pc_y[i];
    }

    // Set up 156-162 components of targeted buffer : ID

    auto ta_yzzzzz_xx_0 = pbuffer.data(idx_npot_0_id + 156);

    auto ta_yzzzzz_xy_0 = pbuffer.data(idx_npot_0_id + 157);

    auto ta_yzzzzz_xz_0 = pbuffer.data(idx_npot_0_id + 158);

    auto ta_yzzzzz_yy_0 = pbuffer.data(idx_npot_0_id + 159);

    auto ta_yzzzzz_yz_0 = pbuffer.data(idx_npot_0_id + 160);

    auto ta_yzzzzz_zz_0 = pbuffer.data(idx_npot_0_id + 161);

#pragma omp simd aligned(pa_y,               \
                             pc_y,           \
                             ta_yzzzzz_xx_0, \
                             ta_yzzzzz_xy_0, \
                             ta_yzzzzz_xz_0, \
                             ta_yzzzzz_yy_0, \
                             ta_yzzzzz_yz_0, \
                             ta_yzzzzz_zz_0, \
                             ta_zzzzz_x_0,   \
                             ta_zzzzz_x_1,   \
                             ta_zzzzz_xx_0,  \
                             ta_zzzzz_xx_1,  \
                             ta_zzzzz_xy_0,  \
                             ta_zzzzz_xy_1,  \
                             ta_zzzzz_xz_0,  \
                             ta_zzzzz_xz_1,  \
                             ta_zzzzz_y_0,   \
                             ta_zzzzz_y_1,   \
                             ta_zzzzz_yy_0,  \
                             ta_zzzzz_yy_1,  \
                             ta_zzzzz_yz_0,  \
                             ta_zzzzz_yz_1,  \
                             ta_zzzzz_z_0,   \
                             ta_zzzzz_z_1,   \
                             ta_zzzzz_zz_0,  \
                             ta_zzzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzzzzz_xx_0[i] = ta_zzzzz_xx_0[i] * pa_y[i] - ta_zzzzz_xx_1[i] * pc_y[i];

        ta_yzzzzz_xy_0[i] = ta_zzzzz_x_0[i] * fe_0 - ta_zzzzz_x_1[i] * fe_0 + ta_zzzzz_xy_0[i] * pa_y[i] - ta_zzzzz_xy_1[i] * pc_y[i];

        ta_yzzzzz_xz_0[i] = ta_zzzzz_xz_0[i] * pa_y[i] - ta_zzzzz_xz_1[i] * pc_y[i];

        ta_yzzzzz_yy_0[i] = 2.0 * ta_zzzzz_y_0[i] * fe_0 - 2.0 * ta_zzzzz_y_1[i] * fe_0 + ta_zzzzz_yy_0[i] * pa_y[i] - ta_zzzzz_yy_1[i] * pc_y[i];

        ta_yzzzzz_yz_0[i] = ta_zzzzz_z_0[i] * fe_0 - ta_zzzzz_z_1[i] * fe_0 + ta_zzzzz_yz_0[i] * pa_y[i] - ta_zzzzz_yz_1[i] * pc_y[i];

        ta_yzzzzz_zz_0[i] = ta_zzzzz_zz_0[i] * pa_y[i] - ta_zzzzz_zz_1[i] * pc_y[i];
    }

    // Set up 162-168 components of targeted buffer : ID

    auto ta_zzzzzz_xx_0 = pbuffer.data(idx_npot_0_id + 162);

    auto ta_zzzzzz_xy_0 = pbuffer.data(idx_npot_0_id + 163);

    auto ta_zzzzzz_xz_0 = pbuffer.data(idx_npot_0_id + 164);

    auto ta_zzzzzz_yy_0 = pbuffer.data(idx_npot_0_id + 165);

    auto ta_zzzzzz_yz_0 = pbuffer.data(idx_npot_0_id + 166);

    auto ta_zzzzzz_zz_0 = pbuffer.data(idx_npot_0_id + 167);

#pragma omp simd aligned(pa_z,               \
                             pc_z,           \
                             ta_zzzz_xx_0,   \
                             ta_zzzz_xx_1,   \
                             ta_zzzz_xy_0,   \
                             ta_zzzz_xy_1,   \
                             ta_zzzz_xz_0,   \
                             ta_zzzz_xz_1,   \
                             ta_zzzz_yy_0,   \
                             ta_zzzz_yy_1,   \
                             ta_zzzz_yz_0,   \
                             ta_zzzz_yz_1,   \
                             ta_zzzz_zz_0,   \
                             ta_zzzz_zz_1,   \
                             ta_zzzzz_x_0,   \
                             ta_zzzzz_x_1,   \
                             ta_zzzzz_xx_0,  \
                             ta_zzzzz_xx_1,  \
                             ta_zzzzz_xy_0,  \
                             ta_zzzzz_xy_1,  \
                             ta_zzzzz_xz_0,  \
                             ta_zzzzz_xz_1,  \
                             ta_zzzzz_y_0,   \
                             ta_zzzzz_y_1,   \
                             ta_zzzzz_yy_0,  \
                             ta_zzzzz_yy_1,  \
                             ta_zzzzz_yz_0,  \
                             ta_zzzzz_yz_1,  \
                             ta_zzzzz_z_0,   \
                             ta_zzzzz_z_1,   \
                             ta_zzzzz_zz_0,  \
                             ta_zzzzz_zz_1,  \
                             ta_zzzzzz_xx_0, \
                             ta_zzzzzz_xy_0, \
                             ta_zzzzzz_xz_0, \
                             ta_zzzzzz_yy_0, \
                             ta_zzzzzz_yz_0, \
                             ta_zzzzzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzzzzz_xx_0[i] = 5.0 * ta_zzzz_xx_0[i] * fe_0 - 5.0 * ta_zzzz_xx_1[i] * fe_0 + ta_zzzzz_xx_0[i] * pa_z[i] - ta_zzzzz_xx_1[i] * pc_z[i];

        ta_zzzzzz_xy_0[i] = 5.0 * ta_zzzz_xy_0[i] * fe_0 - 5.0 * ta_zzzz_xy_1[i] * fe_0 + ta_zzzzz_xy_0[i] * pa_z[i] - ta_zzzzz_xy_1[i] * pc_z[i];

        ta_zzzzzz_xz_0[i] = 5.0 * ta_zzzz_xz_0[i] * fe_0 - 5.0 * ta_zzzz_xz_1[i] * fe_0 + ta_zzzzz_x_0[i] * fe_0 - ta_zzzzz_x_1[i] * fe_0 +
                            ta_zzzzz_xz_0[i] * pa_z[i] - ta_zzzzz_xz_1[i] * pc_z[i];

        ta_zzzzzz_yy_0[i] = 5.0 * ta_zzzz_yy_0[i] * fe_0 - 5.0 * ta_zzzz_yy_1[i] * fe_0 + ta_zzzzz_yy_0[i] * pa_z[i] - ta_zzzzz_yy_1[i] * pc_z[i];

        ta_zzzzzz_yz_0[i] = 5.0 * ta_zzzz_yz_0[i] * fe_0 - 5.0 * ta_zzzz_yz_1[i] * fe_0 + ta_zzzzz_y_0[i] * fe_0 - ta_zzzzz_y_1[i] * fe_0 +
                            ta_zzzzz_yz_0[i] * pa_z[i] - ta_zzzzz_yz_1[i] * pc_z[i];

        ta_zzzzzz_zz_0[i] = 5.0 * ta_zzzz_zz_0[i] * fe_0 - 5.0 * ta_zzzz_zz_1[i] * fe_0 + 2.0 * ta_zzzzz_z_0[i] * fe_0 -
                            2.0 * ta_zzzzz_z_1[i] * fe_0 + ta_zzzzz_zz_0[i] * pa_z[i] - ta_zzzzz_zz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
