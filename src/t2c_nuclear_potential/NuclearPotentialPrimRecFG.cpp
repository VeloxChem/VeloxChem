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

#include "NuclearPotentialPrimRecFG.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_fg(CSimdArray<double>&       pbuffer,
                               const size_t              idx_npot_0_fg,
                               const size_t              idx_npot_0_pg,
                               const size_t              idx_npot_1_pg,
                               const size_t              idx_npot_0_df,
                               const size_t              idx_npot_1_df,
                               const size_t              idx_npot_0_dg,
                               const size_t              idx_npot_1_dg,
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

    // Set up components of auxiliary buffer : PG

    auto ta_x_xxxx_0 = pbuffer.data(idx_npot_0_pg);

    auto ta_x_xxxy_0 = pbuffer.data(idx_npot_0_pg + 1);

    auto ta_x_xxxz_0 = pbuffer.data(idx_npot_0_pg + 2);

    auto ta_x_xxyy_0 = pbuffer.data(idx_npot_0_pg + 3);

    auto ta_x_xxyz_0 = pbuffer.data(idx_npot_0_pg + 4);

    auto ta_x_xxzz_0 = pbuffer.data(idx_npot_0_pg + 5);

    auto ta_x_xyyy_0 = pbuffer.data(idx_npot_0_pg + 6);

    auto ta_x_xyyz_0 = pbuffer.data(idx_npot_0_pg + 7);

    auto ta_x_xyzz_0 = pbuffer.data(idx_npot_0_pg + 8);

    auto ta_x_xzzz_0 = pbuffer.data(idx_npot_0_pg + 9);

    auto ta_x_yyyy_0 = pbuffer.data(idx_npot_0_pg + 10);

    auto ta_x_yyyz_0 = pbuffer.data(idx_npot_0_pg + 11);

    auto ta_x_yyzz_0 = pbuffer.data(idx_npot_0_pg + 12);

    auto ta_x_yzzz_0 = pbuffer.data(idx_npot_0_pg + 13);

    auto ta_x_zzzz_0 = pbuffer.data(idx_npot_0_pg + 14);

    auto ta_y_xxxx_0 = pbuffer.data(idx_npot_0_pg + 15);

    auto ta_y_xxxy_0 = pbuffer.data(idx_npot_0_pg + 16);

    auto ta_y_xxxz_0 = pbuffer.data(idx_npot_0_pg + 17);

    auto ta_y_xxyy_0 = pbuffer.data(idx_npot_0_pg + 18);

    auto ta_y_xxyz_0 = pbuffer.data(idx_npot_0_pg + 19);

    auto ta_y_xxzz_0 = pbuffer.data(idx_npot_0_pg + 20);

    auto ta_y_xyyy_0 = pbuffer.data(idx_npot_0_pg + 21);

    auto ta_y_xyyz_0 = pbuffer.data(idx_npot_0_pg + 22);

    auto ta_y_xyzz_0 = pbuffer.data(idx_npot_0_pg + 23);

    auto ta_y_xzzz_0 = pbuffer.data(idx_npot_0_pg + 24);

    auto ta_y_yyyy_0 = pbuffer.data(idx_npot_0_pg + 25);

    auto ta_y_yyyz_0 = pbuffer.data(idx_npot_0_pg + 26);

    auto ta_y_yyzz_0 = pbuffer.data(idx_npot_0_pg + 27);

    auto ta_y_yzzz_0 = pbuffer.data(idx_npot_0_pg + 28);

    auto ta_y_zzzz_0 = pbuffer.data(idx_npot_0_pg + 29);

    auto ta_z_xxxx_0 = pbuffer.data(idx_npot_0_pg + 30);

    auto ta_z_xxxy_0 = pbuffer.data(idx_npot_0_pg + 31);

    auto ta_z_xxxz_0 = pbuffer.data(idx_npot_0_pg + 32);

    auto ta_z_xxyy_0 = pbuffer.data(idx_npot_0_pg + 33);

    auto ta_z_xxyz_0 = pbuffer.data(idx_npot_0_pg + 34);

    auto ta_z_xxzz_0 = pbuffer.data(idx_npot_0_pg + 35);

    auto ta_z_xyyy_0 = pbuffer.data(idx_npot_0_pg + 36);

    auto ta_z_xyyz_0 = pbuffer.data(idx_npot_0_pg + 37);

    auto ta_z_xyzz_0 = pbuffer.data(idx_npot_0_pg + 38);

    auto ta_z_xzzz_0 = pbuffer.data(idx_npot_0_pg + 39);

    auto ta_z_yyyy_0 = pbuffer.data(idx_npot_0_pg + 40);

    auto ta_z_yyyz_0 = pbuffer.data(idx_npot_0_pg + 41);

    auto ta_z_yyzz_0 = pbuffer.data(idx_npot_0_pg + 42);

    auto ta_z_yzzz_0 = pbuffer.data(idx_npot_0_pg + 43);

    auto ta_z_zzzz_0 = pbuffer.data(idx_npot_0_pg + 44);

    // Set up components of auxiliary buffer : PG

    auto ta_x_xxxx_1 = pbuffer.data(idx_npot_1_pg);

    auto ta_x_xxxy_1 = pbuffer.data(idx_npot_1_pg + 1);

    auto ta_x_xxxz_1 = pbuffer.data(idx_npot_1_pg + 2);

    auto ta_x_xxyy_1 = pbuffer.data(idx_npot_1_pg + 3);

    auto ta_x_xxyz_1 = pbuffer.data(idx_npot_1_pg + 4);

    auto ta_x_xxzz_1 = pbuffer.data(idx_npot_1_pg + 5);

    auto ta_x_xyyy_1 = pbuffer.data(idx_npot_1_pg + 6);

    auto ta_x_xyyz_1 = pbuffer.data(idx_npot_1_pg + 7);

    auto ta_x_xyzz_1 = pbuffer.data(idx_npot_1_pg + 8);

    auto ta_x_xzzz_1 = pbuffer.data(idx_npot_1_pg + 9);

    auto ta_x_yyyy_1 = pbuffer.data(idx_npot_1_pg + 10);

    auto ta_x_yyyz_1 = pbuffer.data(idx_npot_1_pg + 11);

    auto ta_x_yyzz_1 = pbuffer.data(idx_npot_1_pg + 12);

    auto ta_x_yzzz_1 = pbuffer.data(idx_npot_1_pg + 13);

    auto ta_x_zzzz_1 = pbuffer.data(idx_npot_1_pg + 14);

    auto ta_y_xxxx_1 = pbuffer.data(idx_npot_1_pg + 15);

    auto ta_y_xxxy_1 = pbuffer.data(idx_npot_1_pg + 16);

    auto ta_y_xxxz_1 = pbuffer.data(idx_npot_1_pg + 17);

    auto ta_y_xxyy_1 = pbuffer.data(idx_npot_1_pg + 18);

    auto ta_y_xxyz_1 = pbuffer.data(idx_npot_1_pg + 19);

    auto ta_y_xxzz_1 = pbuffer.data(idx_npot_1_pg + 20);

    auto ta_y_xyyy_1 = pbuffer.data(idx_npot_1_pg + 21);

    auto ta_y_xyyz_1 = pbuffer.data(idx_npot_1_pg + 22);

    auto ta_y_xyzz_1 = pbuffer.data(idx_npot_1_pg + 23);

    auto ta_y_xzzz_1 = pbuffer.data(idx_npot_1_pg + 24);

    auto ta_y_yyyy_1 = pbuffer.data(idx_npot_1_pg + 25);

    auto ta_y_yyyz_1 = pbuffer.data(idx_npot_1_pg + 26);

    auto ta_y_yyzz_1 = pbuffer.data(idx_npot_1_pg + 27);

    auto ta_y_yzzz_1 = pbuffer.data(idx_npot_1_pg + 28);

    auto ta_y_zzzz_1 = pbuffer.data(idx_npot_1_pg + 29);

    auto ta_z_xxxx_1 = pbuffer.data(idx_npot_1_pg + 30);

    auto ta_z_xxxy_1 = pbuffer.data(idx_npot_1_pg + 31);

    auto ta_z_xxxz_1 = pbuffer.data(idx_npot_1_pg + 32);

    auto ta_z_xxyy_1 = pbuffer.data(idx_npot_1_pg + 33);

    auto ta_z_xxyz_1 = pbuffer.data(idx_npot_1_pg + 34);

    auto ta_z_xxzz_1 = pbuffer.data(idx_npot_1_pg + 35);

    auto ta_z_xyyy_1 = pbuffer.data(idx_npot_1_pg + 36);

    auto ta_z_xyyz_1 = pbuffer.data(idx_npot_1_pg + 37);

    auto ta_z_xyzz_1 = pbuffer.data(idx_npot_1_pg + 38);

    auto ta_z_xzzz_1 = pbuffer.data(idx_npot_1_pg + 39);

    auto ta_z_yyyy_1 = pbuffer.data(idx_npot_1_pg + 40);

    auto ta_z_yyyz_1 = pbuffer.data(idx_npot_1_pg + 41);

    auto ta_z_yyzz_1 = pbuffer.data(idx_npot_1_pg + 42);

    auto ta_z_yzzz_1 = pbuffer.data(idx_npot_1_pg + 43);

    auto ta_z_zzzz_1 = pbuffer.data(idx_npot_1_pg + 44);

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

    auto ta_yz_xyz_0 = pbuffer.data(idx_npot_0_df + 44);

    auto ta_yz_yyz_0 = pbuffer.data(idx_npot_0_df + 47);

    auto ta_yz_yzz_0 = pbuffer.data(idx_npot_0_df + 48);

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

    auto ta_yz_xyz_1 = pbuffer.data(idx_npot_1_df + 44);

    auto ta_yz_yyz_1 = pbuffer.data(idx_npot_1_df + 47);

    auto ta_yz_yzz_1 = pbuffer.data(idx_npot_1_df + 48);

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

    // Set up components of auxiliary buffer : DG

    auto ta_xx_xxxx_0 = pbuffer.data(idx_npot_0_dg);

    auto ta_xx_xxxy_0 = pbuffer.data(idx_npot_0_dg + 1);

    auto ta_xx_xxxz_0 = pbuffer.data(idx_npot_0_dg + 2);

    auto ta_xx_xxyy_0 = pbuffer.data(idx_npot_0_dg + 3);

    auto ta_xx_xxyz_0 = pbuffer.data(idx_npot_0_dg + 4);

    auto ta_xx_xxzz_0 = pbuffer.data(idx_npot_0_dg + 5);

    auto ta_xx_xyyy_0 = pbuffer.data(idx_npot_0_dg + 6);

    auto ta_xx_xyyz_0 = pbuffer.data(idx_npot_0_dg + 7);

    auto ta_xx_xyzz_0 = pbuffer.data(idx_npot_0_dg + 8);

    auto ta_xx_xzzz_0 = pbuffer.data(idx_npot_0_dg + 9);

    auto ta_xx_yyyy_0 = pbuffer.data(idx_npot_0_dg + 10);

    auto ta_xx_yyyz_0 = pbuffer.data(idx_npot_0_dg + 11);

    auto ta_xx_yyzz_0 = pbuffer.data(idx_npot_0_dg + 12);

    auto ta_xx_yzzz_0 = pbuffer.data(idx_npot_0_dg + 13);

    auto ta_xx_zzzz_0 = pbuffer.data(idx_npot_0_dg + 14);

    auto ta_xy_xxxy_0 = pbuffer.data(idx_npot_0_dg + 16);

    auto ta_xy_xxyy_0 = pbuffer.data(idx_npot_0_dg + 18);

    auto ta_xy_xyyy_0 = pbuffer.data(idx_npot_0_dg + 21);

    auto ta_xy_yyyy_0 = pbuffer.data(idx_npot_0_dg + 25);

    auto ta_xy_yyyz_0 = pbuffer.data(idx_npot_0_dg + 26);

    auto ta_xy_yyzz_0 = pbuffer.data(idx_npot_0_dg + 27);

    auto ta_xy_yzzz_0 = pbuffer.data(idx_npot_0_dg + 28);

    auto ta_xz_xxxx_0 = pbuffer.data(idx_npot_0_dg + 30);

    auto ta_xz_xxxz_0 = pbuffer.data(idx_npot_0_dg + 32);

    auto ta_xz_xxzz_0 = pbuffer.data(idx_npot_0_dg + 35);

    auto ta_xz_xzzz_0 = pbuffer.data(idx_npot_0_dg + 39);

    auto ta_xz_yyyz_0 = pbuffer.data(idx_npot_0_dg + 41);

    auto ta_xz_yyzz_0 = pbuffer.data(idx_npot_0_dg + 42);

    auto ta_xz_yzzz_0 = pbuffer.data(idx_npot_0_dg + 43);

    auto ta_xz_zzzz_0 = pbuffer.data(idx_npot_0_dg + 44);

    auto ta_yy_xxxx_0 = pbuffer.data(idx_npot_0_dg + 45);

    auto ta_yy_xxxy_0 = pbuffer.data(idx_npot_0_dg + 46);

    auto ta_yy_xxxz_0 = pbuffer.data(idx_npot_0_dg + 47);

    auto ta_yy_xxyy_0 = pbuffer.data(idx_npot_0_dg + 48);

    auto ta_yy_xxyz_0 = pbuffer.data(idx_npot_0_dg + 49);

    auto ta_yy_xxzz_0 = pbuffer.data(idx_npot_0_dg + 50);

    auto ta_yy_xyyy_0 = pbuffer.data(idx_npot_0_dg + 51);

    auto ta_yy_xyyz_0 = pbuffer.data(idx_npot_0_dg + 52);

    auto ta_yy_xyzz_0 = pbuffer.data(idx_npot_0_dg + 53);

    auto ta_yy_xzzz_0 = pbuffer.data(idx_npot_0_dg + 54);

    auto ta_yy_yyyy_0 = pbuffer.data(idx_npot_0_dg + 55);

    auto ta_yy_yyyz_0 = pbuffer.data(idx_npot_0_dg + 56);

    auto ta_yy_yyzz_0 = pbuffer.data(idx_npot_0_dg + 57);

    auto ta_yy_yzzz_0 = pbuffer.data(idx_npot_0_dg + 58);

    auto ta_yy_zzzz_0 = pbuffer.data(idx_npot_0_dg + 59);

    auto ta_yz_xxxz_0 = pbuffer.data(idx_npot_0_dg + 62);

    auto ta_yz_xxyz_0 = pbuffer.data(idx_npot_0_dg + 64);

    auto ta_yz_xxzz_0 = pbuffer.data(idx_npot_0_dg + 65);

    auto ta_yz_xyyz_0 = pbuffer.data(idx_npot_0_dg + 67);

    auto ta_yz_xyzz_0 = pbuffer.data(idx_npot_0_dg + 68);

    auto ta_yz_xzzz_0 = pbuffer.data(idx_npot_0_dg + 69);

    auto ta_yz_yyyy_0 = pbuffer.data(idx_npot_0_dg + 70);

    auto ta_yz_yyyz_0 = pbuffer.data(idx_npot_0_dg + 71);

    auto ta_yz_yyzz_0 = pbuffer.data(idx_npot_0_dg + 72);

    auto ta_yz_yzzz_0 = pbuffer.data(idx_npot_0_dg + 73);

    auto ta_yz_zzzz_0 = pbuffer.data(idx_npot_0_dg + 74);

    auto ta_zz_xxxx_0 = pbuffer.data(idx_npot_0_dg + 75);

    auto ta_zz_xxxy_0 = pbuffer.data(idx_npot_0_dg + 76);

    auto ta_zz_xxxz_0 = pbuffer.data(idx_npot_0_dg + 77);

    auto ta_zz_xxyy_0 = pbuffer.data(idx_npot_0_dg + 78);

    auto ta_zz_xxyz_0 = pbuffer.data(idx_npot_0_dg + 79);

    auto ta_zz_xxzz_0 = pbuffer.data(idx_npot_0_dg + 80);

    auto ta_zz_xyyy_0 = pbuffer.data(idx_npot_0_dg + 81);

    auto ta_zz_xyyz_0 = pbuffer.data(idx_npot_0_dg + 82);

    auto ta_zz_xyzz_0 = pbuffer.data(idx_npot_0_dg + 83);

    auto ta_zz_xzzz_0 = pbuffer.data(idx_npot_0_dg + 84);

    auto ta_zz_yyyy_0 = pbuffer.data(idx_npot_0_dg + 85);

    auto ta_zz_yyyz_0 = pbuffer.data(idx_npot_0_dg + 86);

    auto ta_zz_yyzz_0 = pbuffer.data(idx_npot_0_dg + 87);

    auto ta_zz_yzzz_0 = pbuffer.data(idx_npot_0_dg + 88);

    auto ta_zz_zzzz_0 = pbuffer.data(idx_npot_0_dg + 89);

    // Set up components of auxiliary buffer : DG

    auto ta_xx_xxxx_1 = pbuffer.data(idx_npot_1_dg);

    auto ta_xx_xxxy_1 = pbuffer.data(idx_npot_1_dg + 1);

    auto ta_xx_xxxz_1 = pbuffer.data(idx_npot_1_dg + 2);

    auto ta_xx_xxyy_1 = pbuffer.data(idx_npot_1_dg + 3);

    auto ta_xx_xxyz_1 = pbuffer.data(idx_npot_1_dg + 4);

    auto ta_xx_xxzz_1 = pbuffer.data(idx_npot_1_dg + 5);

    auto ta_xx_xyyy_1 = pbuffer.data(idx_npot_1_dg + 6);

    auto ta_xx_xyyz_1 = pbuffer.data(idx_npot_1_dg + 7);

    auto ta_xx_xyzz_1 = pbuffer.data(idx_npot_1_dg + 8);

    auto ta_xx_xzzz_1 = pbuffer.data(idx_npot_1_dg + 9);

    auto ta_xx_yyyy_1 = pbuffer.data(idx_npot_1_dg + 10);

    auto ta_xx_yyyz_1 = pbuffer.data(idx_npot_1_dg + 11);

    auto ta_xx_yyzz_1 = pbuffer.data(idx_npot_1_dg + 12);

    auto ta_xx_yzzz_1 = pbuffer.data(idx_npot_1_dg + 13);

    auto ta_xx_zzzz_1 = pbuffer.data(idx_npot_1_dg + 14);

    auto ta_xy_xxxy_1 = pbuffer.data(idx_npot_1_dg + 16);

    auto ta_xy_xxyy_1 = pbuffer.data(idx_npot_1_dg + 18);

    auto ta_xy_xyyy_1 = pbuffer.data(idx_npot_1_dg + 21);

    auto ta_xy_yyyy_1 = pbuffer.data(idx_npot_1_dg + 25);

    auto ta_xy_yyyz_1 = pbuffer.data(idx_npot_1_dg + 26);

    auto ta_xy_yyzz_1 = pbuffer.data(idx_npot_1_dg + 27);

    auto ta_xy_yzzz_1 = pbuffer.data(idx_npot_1_dg + 28);

    auto ta_xz_xxxx_1 = pbuffer.data(idx_npot_1_dg + 30);

    auto ta_xz_xxxz_1 = pbuffer.data(idx_npot_1_dg + 32);

    auto ta_xz_xxzz_1 = pbuffer.data(idx_npot_1_dg + 35);

    auto ta_xz_xzzz_1 = pbuffer.data(idx_npot_1_dg + 39);

    auto ta_xz_yyyz_1 = pbuffer.data(idx_npot_1_dg + 41);

    auto ta_xz_yyzz_1 = pbuffer.data(idx_npot_1_dg + 42);

    auto ta_xz_yzzz_1 = pbuffer.data(idx_npot_1_dg + 43);

    auto ta_xz_zzzz_1 = pbuffer.data(idx_npot_1_dg + 44);

    auto ta_yy_xxxx_1 = pbuffer.data(idx_npot_1_dg + 45);

    auto ta_yy_xxxy_1 = pbuffer.data(idx_npot_1_dg + 46);

    auto ta_yy_xxxz_1 = pbuffer.data(idx_npot_1_dg + 47);

    auto ta_yy_xxyy_1 = pbuffer.data(idx_npot_1_dg + 48);

    auto ta_yy_xxyz_1 = pbuffer.data(idx_npot_1_dg + 49);

    auto ta_yy_xxzz_1 = pbuffer.data(idx_npot_1_dg + 50);

    auto ta_yy_xyyy_1 = pbuffer.data(idx_npot_1_dg + 51);

    auto ta_yy_xyyz_1 = pbuffer.data(idx_npot_1_dg + 52);

    auto ta_yy_xyzz_1 = pbuffer.data(idx_npot_1_dg + 53);

    auto ta_yy_xzzz_1 = pbuffer.data(idx_npot_1_dg + 54);

    auto ta_yy_yyyy_1 = pbuffer.data(idx_npot_1_dg + 55);

    auto ta_yy_yyyz_1 = pbuffer.data(idx_npot_1_dg + 56);

    auto ta_yy_yyzz_1 = pbuffer.data(idx_npot_1_dg + 57);

    auto ta_yy_yzzz_1 = pbuffer.data(idx_npot_1_dg + 58);

    auto ta_yy_zzzz_1 = pbuffer.data(idx_npot_1_dg + 59);

    auto ta_yz_xxxz_1 = pbuffer.data(idx_npot_1_dg + 62);

    auto ta_yz_xxyz_1 = pbuffer.data(idx_npot_1_dg + 64);

    auto ta_yz_xxzz_1 = pbuffer.data(idx_npot_1_dg + 65);

    auto ta_yz_xyyz_1 = pbuffer.data(idx_npot_1_dg + 67);

    auto ta_yz_xyzz_1 = pbuffer.data(idx_npot_1_dg + 68);

    auto ta_yz_xzzz_1 = pbuffer.data(idx_npot_1_dg + 69);

    auto ta_yz_yyyy_1 = pbuffer.data(idx_npot_1_dg + 70);

    auto ta_yz_yyyz_1 = pbuffer.data(idx_npot_1_dg + 71);

    auto ta_yz_yyzz_1 = pbuffer.data(idx_npot_1_dg + 72);

    auto ta_yz_yzzz_1 = pbuffer.data(idx_npot_1_dg + 73);

    auto ta_yz_zzzz_1 = pbuffer.data(idx_npot_1_dg + 74);

    auto ta_zz_xxxx_1 = pbuffer.data(idx_npot_1_dg + 75);

    auto ta_zz_xxxy_1 = pbuffer.data(idx_npot_1_dg + 76);

    auto ta_zz_xxxz_1 = pbuffer.data(idx_npot_1_dg + 77);

    auto ta_zz_xxyy_1 = pbuffer.data(idx_npot_1_dg + 78);

    auto ta_zz_xxyz_1 = pbuffer.data(idx_npot_1_dg + 79);

    auto ta_zz_xxzz_1 = pbuffer.data(idx_npot_1_dg + 80);

    auto ta_zz_xyyy_1 = pbuffer.data(idx_npot_1_dg + 81);

    auto ta_zz_xyyz_1 = pbuffer.data(idx_npot_1_dg + 82);

    auto ta_zz_xyzz_1 = pbuffer.data(idx_npot_1_dg + 83);

    auto ta_zz_xzzz_1 = pbuffer.data(idx_npot_1_dg + 84);

    auto ta_zz_yyyy_1 = pbuffer.data(idx_npot_1_dg + 85);

    auto ta_zz_yyyz_1 = pbuffer.data(idx_npot_1_dg + 86);

    auto ta_zz_yyzz_1 = pbuffer.data(idx_npot_1_dg + 87);

    auto ta_zz_yzzz_1 = pbuffer.data(idx_npot_1_dg + 88);

    auto ta_zz_zzzz_1 = pbuffer.data(idx_npot_1_dg + 89);

    // Set up 0-15 components of targeted buffer : FG

    auto ta_xxx_xxxx_0 = pbuffer.data(idx_npot_0_fg);

    auto ta_xxx_xxxy_0 = pbuffer.data(idx_npot_0_fg + 1);

    auto ta_xxx_xxxz_0 = pbuffer.data(idx_npot_0_fg + 2);

    auto ta_xxx_xxyy_0 = pbuffer.data(idx_npot_0_fg + 3);

    auto ta_xxx_xxyz_0 = pbuffer.data(idx_npot_0_fg + 4);

    auto ta_xxx_xxzz_0 = pbuffer.data(idx_npot_0_fg + 5);

    auto ta_xxx_xyyy_0 = pbuffer.data(idx_npot_0_fg + 6);

    auto ta_xxx_xyyz_0 = pbuffer.data(idx_npot_0_fg + 7);

    auto ta_xxx_xyzz_0 = pbuffer.data(idx_npot_0_fg + 8);

    auto ta_xxx_xzzz_0 = pbuffer.data(idx_npot_0_fg + 9);

    auto ta_xxx_yyyy_0 = pbuffer.data(idx_npot_0_fg + 10);

    auto ta_xxx_yyyz_0 = pbuffer.data(idx_npot_0_fg + 11);

    auto ta_xxx_yyzz_0 = pbuffer.data(idx_npot_0_fg + 12);

    auto ta_xxx_yzzz_0 = pbuffer.data(idx_npot_0_fg + 13);

    auto ta_xxx_zzzz_0 = pbuffer.data(idx_npot_0_fg + 14);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta_x_xxxx_0,   \
                             ta_x_xxxx_1,   \
                             ta_x_xxxy_0,   \
                             ta_x_xxxy_1,   \
                             ta_x_xxxz_0,   \
                             ta_x_xxxz_1,   \
                             ta_x_xxyy_0,   \
                             ta_x_xxyy_1,   \
                             ta_x_xxyz_0,   \
                             ta_x_xxyz_1,   \
                             ta_x_xxzz_0,   \
                             ta_x_xxzz_1,   \
                             ta_x_xyyy_0,   \
                             ta_x_xyyy_1,   \
                             ta_x_xyyz_0,   \
                             ta_x_xyyz_1,   \
                             ta_x_xyzz_0,   \
                             ta_x_xyzz_1,   \
                             ta_x_xzzz_0,   \
                             ta_x_xzzz_1,   \
                             ta_x_yyyy_0,   \
                             ta_x_yyyy_1,   \
                             ta_x_yyyz_0,   \
                             ta_x_yyyz_1,   \
                             ta_x_yyzz_0,   \
                             ta_x_yyzz_1,   \
                             ta_x_yzzz_0,   \
                             ta_x_yzzz_1,   \
                             ta_x_zzzz_0,   \
                             ta_x_zzzz_1,   \
                             ta_xx_xxx_0,   \
                             ta_xx_xxx_1,   \
                             ta_xx_xxxx_0,  \
                             ta_xx_xxxx_1,  \
                             ta_xx_xxxy_0,  \
                             ta_xx_xxxy_1,  \
                             ta_xx_xxxz_0,  \
                             ta_xx_xxxz_1,  \
                             ta_xx_xxy_0,   \
                             ta_xx_xxy_1,   \
                             ta_xx_xxyy_0,  \
                             ta_xx_xxyy_1,  \
                             ta_xx_xxyz_0,  \
                             ta_xx_xxyz_1,  \
                             ta_xx_xxz_0,   \
                             ta_xx_xxz_1,   \
                             ta_xx_xxzz_0,  \
                             ta_xx_xxzz_1,  \
                             ta_xx_xyy_0,   \
                             ta_xx_xyy_1,   \
                             ta_xx_xyyy_0,  \
                             ta_xx_xyyy_1,  \
                             ta_xx_xyyz_0,  \
                             ta_xx_xyyz_1,  \
                             ta_xx_xyz_0,   \
                             ta_xx_xyz_1,   \
                             ta_xx_xyzz_0,  \
                             ta_xx_xyzz_1,  \
                             ta_xx_xzz_0,   \
                             ta_xx_xzz_1,   \
                             ta_xx_xzzz_0,  \
                             ta_xx_xzzz_1,  \
                             ta_xx_yyy_0,   \
                             ta_xx_yyy_1,   \
                             ta_xx_yyyy_0,  \
                             ta_xx_yyyy_1,  \
                             ta_xx_yyyz_0,  \
                             ta_xx_yyyz_1,  \
                             ta_xx_yyz_0,   \
                             ta_xx_yyz_1,   \
                             ta_xx_yyzz_0,  \
                             ta_xx_yyzz_1,  \
                             ta_xx_yzz_0,   \
                             ta_xx_yzz_1,   \
                             ta_xx_yzzz_0,  \
                             ta_xx_yzzz_1,  \
                             ta_xx_zzz_0,   \
                             ta_xx_zzz_1,   \
                             ta_xx_zzzz_0,  \
                             ta_xx_zzzz_1,  \
                             ta_xxx_xxxx_0, \
                             ta_xxx_xxxy_0, \
                             ta_xxx_xxxz_0, \
                             ta_xxx_xxyy_0, \
                             ta_xxx_xxyz_0, \
                             ta_xxx_xxzz_0, \
                             ta_xxx_xyyy_0, \
                             ta_xxx_xyyz_0, \
                             ta_xxx_xyzz_0, \
                             ta_xxx_xzzz_0, \
                             ta_xxx_yyyy_0, \
                             ta_xxx_yyyz_0, \
                             ta_xxx_yyzz_0, \
                             ta_xxx_yzzz_0, \
                             ta_xxx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxx_xxxx_0[i] = 2.0 * ta_x_xxxx_0[i] * fe_0 - 2.0 * ta_x_xxxx_1[i] * fe_0 + 4.0 * ta_xx_xxx_0[i] * fe_0 - 4.0 * ta_xx_xxx_1[i] * fe_0 +
                           ta_xx_xxxx_0[i] * pa_x[i] - ta_xx_xxxx_1[i] * pc_x[i];

        ta_xxx_xxxy_0[i] = 2.0 * ta_x_xxxy_0[i] * fe_0 - 2.0 * ta_x_xxxy_1[i] * fe_0 + 3.0 * ta_xx_xxy_0[i] * fe_0 - 3.0 * ta_xx_xxy_1[i] * fe_0 +
                           ta_xx_xxxy_0[i] * pa_x[i] - ta_xx_xxxy_1[i] * pc_x[i];

        ta_xxx_xxxz_0[i] = 2.0 * ta_x_xxxz_0[i] * fe_0 - 2.0 * ta_x_xxxz_1[i] * fe_0 + 3.0 * ta_xx_xxz_0[i] * fe_0 - 3.0 * ta_xx_xxz_1[i] * fe_0 +
                           ta_xx_xxxz_0[i] * pa_x[i] - ta_xx_xxxz_1[i] * pc_x[i];

        ta_xxx_xxyy_0[i] = 2.0 * ta_x_xxyy_0[i] * fe_0 - 2.0 * ta_x_xxyy_1[i] * fe_0 + 2.0 * ta_xx_xyy_0[i] * fe_0 - 2.0 * ta_xx_xyy_1[i] * fe_0 +
                           ta_xx_xxyy_0[i] * pa_x[i] - ta_xx_xxyy_1[i] * pc_x[i];

        ta_xxx_xxyz_0[i] = 2.0 * ta_x_xxyz_0[i] * fe_0 - 2.0 * ta_x_xxyz_1[i] * fe_0 + 2.0 * ta_xx_xyz_0[i] * fe_0 - 2.0 * ta_xx_xyz_1[i] * fe_0 +
                           ta_xx_xxyz_0[i] * pa_x[i] - ta_xx_xxyz_1[i] * pc_x[i];

        ta_xxx_xxzz_0[i] = 2.0 * ta_x_xxzz_0[i] * fe_0 - 2.0 * ta_x_xxzz_1[i] * fe_0 + 2.0 * ta_xx_xzz_0[i] * fe_0 - 2.0 * ta_xx_xzz_1[i] * fe_0 +
                           ta_xx_xxzz_0[i] * pa_x[i] - ta_xx_xxzz_1[i] * pc_x[i];

        ta_xxx_xyyy_0[i] = 2.0 * ta_x_xyyy_0[i] * fe_0 - 2.0 * ta_x_xyyy_1[i] * fe_0 + ta_xx_yyy_0[i] * fe_0 - ta_xx_yyy_1[i] * fe_0 +
                           ta_xx_xyyy_0[i] * pa_x[i] - ta_xx_xyyy_1[i] * pc_x[i];

        ta_xxx_xyyz_0[i] = 2.0 * ta_x_xyyz_0[i] * fe_0 - 2.0 * ta_x_xyyz_1[i] * fe_0 + ta_xx_yyz_0[i] * fe_0 - ta_xx_yyz_1[i] * fe_0 +
                           ta_xx_xyyz_0[i] * pa_x[i] - ta_xx_xyyz_1[i] * pc_x[i];

        ta_xxx_xyzz_0[i] = 2.0 * ta_x_xyzz_0[i] * fe_0 - 2.0 * ta_x_xyzz_1[i] * fe_0 + ta_xx_yzz_0[i] * fe_0 - ta_xx_yzz_1[i] * fe_0 +
                           ta_xx_xyzz_0[i] * pa_x[i] - ta_xx_xyzz_1[i] * pc_x[i];

        ta_xxx_xzzz_0[i] = 2.0 * ta_x_xzzz_0[i] * fe_0 - 2.0 * ta_x_xzzz_1[i] * fe_0 + ta_xx_zzz_0[i] * fe_0 - ta_xx_zzz_1[i] * fe_0 +
                           ta_xx_xzzz_0[i] * pa_x[i] - ta_xx_xzzz_1[i] * pc_x[i];

        ta_xxx_yyyy_0[i] = 2.0 * ta_x_yyyy_0[i] * fe_0 - 2.0 * ta_x_yyyy_1[i] * fe_0 + ta_xx_yyyy_0[i] * pa_x[i] - ta_xx_yyyy_1[i] * pc_x[i];

        ta_xxx_yyyz_0[i] = 2.0 * ta_x_yyyz_0[i] * fe_0 - 2.0 * ta_x_yyyz_1[i] * fe_0 + ta_xx_yyyz_0[i] * pa_x[i] - ta_xx_yyyz_1[i] * pc_x[i];

        ta_xxx_yyzz_0[i] = 2.0 * ta_x_yyzz_0[i] * fe_0 - 2.0 * ta_x_yyzz_1[i] * fe_0 + ta_xx_yyzz_0[i] * pa_x[i] - ta_xx_yyzz_1[i] * pc_x[i];

        ta_xxx_yzzz_0[i] = 2.0 * ta_x_yzzz_0[i] * fe_0 - 2.0 * ta_x_yzzz_1[i] * fe_0 + ta_xx_yzzz_0[i] * pa_x[i] - ta_xx_yzzz_1[i] * pc_x[i];

        ta_xxx_zzzz_0[i] = 2.0 * ta_x_zzzz_0[i] * fe_0 - 2.0 * ta_x_zzzz_1[i] * fe_0 + ta_xx_zzzz_0[i] * pa_x[i] - ta_xx_zzzz_1[i] * pc_x[i];
    }

    // Set up 15-30 components of targeted buffer : FG

    auto ta_xxy_xxxx_0 = pbuffer.data(idx_npot_0_fg + 15);

    auto ta_xxy_xxxy_0 = pbuffer.data(idx_npot_0_fg + 16);

    auto ta_xxy_xxxz_0 = pbuffer.data(idx_npot_0_fg + 17);

    auto ta_xxy_xxyy_0 = pbuffer.data(idx_npot_0_fg + 18);

    auto ta_xxy_xxyz_0 = pbuffer.data(idx_npot_0_fg + 19);

    auto ta_xxy_xxzz_0 = pbuffer.data(idx_npot_0_fg + 20);

    auto ta_xxy_xyyy_0 = pbuffer.data(idx_npot_0_fg + 21);

    auto ta_xxy_xyyz_0 = pbuffer.data(idx_npot_0_fg + 22);

    auto ta_xxy_xyzz_0 = pbuffer.data(idx_npot_0_fg + 23);

    auto ta_xxy_xzzz_0 = pbuffer.data(idx_npot_0_fg + 24);

    auto ta_xxy_yyyy_0 = pbuffer.data(idx_npot_0_fg + 25);

    auto ta_xxy_yyyz_0 = pbuffer.data(idx_npot_0_fg + 26);

    auto ta_xxy_yyzz_0 = pbuffer.data(idx_npot_0_fg + 27);

    auto ta_xxy_yzzz_0 = pbuffer.data(idx_npot_0_fg + 28);

    auto ta_xxy_zzzz_0 = pbuffer.data(idx_npot_0_fg + 29);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pc_x,          \
                             pc_y,          \
                             ta_xx_xxx_0,   \
                             ta_xx_xxx_1,   \
                             ta_xx_xxxx_0,  \
                             ta_xx_xxxx_1,  \
                             ta_xx_xxxy_0,  \
                             ta_xx_xxxy_1,  \
                             ta_xx_xxxz_0,  \
                             ta_xx_xxxz_1,  \
                             ta_xx_xxy_0,   \
                             ta_xx_xxy_1,   \
                             ta_xx_xxyy_0,  \
                             ta_xx_xxyy_1,  \
                             ta_xx_xxyz_0,  \
                             ta_xx_xxyz_1,  \
                             ta_xx_xxz_0,   \
                             ta_xx_xxz_1,   \
                             ta_xx_xxzz_0,  \
                             ta_xx_xxzz_1,  \
                             ta_xx_xyy_0,   \
                             ta_xx_xyy_1,   \
                             ta_xx_xyyy_0,  \
                             ta_xx_xyyy_1,  \
                             ta_xx_xyyz_0,  \
                             ta_xx_xyyz_1,  \
                             ta_xx_xyz_0,   \
                             ta_xx_xyz_1,   \
                             ta_xx_xyzz_0,  \
                             ta_xx_xyzz_1,  \
                             ta_xx_xzz_0,   \
                             ta_xx_xzz_1,   \
                             ta_xx_xzzz_0,  \
                             ta_xx_xzzz_1,  \
                             ta_xx_zzzz_0,  \
                             ta_xx_zzzz_1,  \
                             ta_xxy_xxxx_0, \
                             ta_xxy_xxxy_0, \
                             ta_xxy_xxxz_0, \
                             ta_xxy_xxyy_0, \
                             ta_xxy_xxyz_0, \
                             ta_xxy_xxzz_0, \
                             ta_xxy_xyyy_0, \
                             ta_xxy_xyyz_0, \
                             ta_xxy_xyzz_0, \
                             ta_xxy_xzzz_0, \
                             ta_xxy_yyyy_0, \
                             ta_xxy_yyyz_0, \
                             ta_xxy_yyzz_0, \
                             ta_xxy_yzzz_0, \
                             ta_xxy_zzzz_0, \
                             ta_xy_yyyy_0,  \
                             ta_xy_yyyy_1,  \
                             ta_xy_yyyz_0,  \
                             ta_xy_yyyz_1,  \
                             ta_xy_yyzz_0,  \
                             ta_xy_yyzz_1,  \
                             ta_xy_yzzz_0,  \
                             ta_xy_yzzz_1,  \
                             ta_y_yyyy_0,   \
                             ta_y_yyyy_1,   \
                             ta_y_yyyz_0,   \
                             ta_y_yyyz_1,   \
                             ta_y_yyzz_0,   \
                             ta_y_yyzz_1,   \
                             ta_y_yzzz_0,   \
                             ta_y_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxy_xxxx_0[i] = ta_xx_xxxx_0[i] * pa_y[i] - ta_xx_xxxx_1[i] * pc_y[i];

        ta_xxy_xxxy_0[i] = ta_xx_xxx_0[i] * fe_0 - ta_xx_xxx_1[i] * fe_0 + ta_xx_xxxy_0[i] * pa_y[i] - ta_xx_xxxy_1[i] * pc_y[i];

        ta_xxy_xxxz_0[i] = ta_xx_xxxz_0[i] * pa_y[i] - ta_xx_xxxz_1[i] * pc_y[i];

        ta_xxy_xxyy_0[i] = 2.0 * ta_xx_xxy_0[i] * fe_0 - 2.0 * ta_xx_xxy_1[i] * fe_0 + ta_xx_xxyy_0[i] * pa_y[i] - ta_xx_xxyy_1[i] * pc_y[i];

        ta_xxy_xxyz_0[i] = ta_xx_xxz_0[i] * fe_0 - ta_xx_xxz_1[i] * fe_0 + ta_xx_xxyz_0[i] * pa_y[i] - ta_xx_xxyz_1[i] * pc_y[i];

        ta_xxy_xxzz_0[i] = ta_xx_xxzz_0[i] * pa_y[i] - ta_xx_xxzz_1[i] * pc_y[i];

        ta_xxy_xyyy_0[i] = 3.0 * ta_xx_xyy_0[i] * fe_0 - 3.0 * ta_xx_xyy_1[i] * fe_0 + ta_xx_xyyy_0[i] * pa_y[i] - ta_xx_xyyy_1[i] * pc_y[i];

        ta_xxy_xyyz_0[i] = 2.0 * ta_xx_xyz_0[i] * fe_0 - 2.0 * ta_xx_xyz_1[i] * fe_0 + ta_xx_xyyz_0[i] * pa_y[i] - ta_xx_xyyz_1[i] * pc_y[i];

        ta_xxy_xyzz_0[i] = ta_xx_xzz_0[i] * fe_0 - ta_xx_xzz_1[i] * fe_0 + ta_xx_xyzz_0[i] * pa_y[i] - ta_xx_xyzz_1[i] * pc_y[i];

        ta_xxy_xzzz_0[i] = ta_xx_xzzz_0[i] * pa_y[i] - ta_xx_xzzz_1[i] * pc_y[i];

        ta_xxy_yyyy_0[i] = ta_y_yyyy_0[i] * fe_0 - ta_y_yyyy_1[i] * fe_0 + ta_xy_yyyy_0[i] * pa_x[i] - ta_xy_yyyy_1[i] * pc_x[i];

        ta_xxy_yyyz_0[i] = ta_y_yyyz_0[i] * fe_0 - ta_y_yyyz_1[i] * fe_0 + ta_xy_yyyz_0[i] * pa_x[i] - ta_xy_yyyz_1[i] * pc_x[i];

        ta_xxy_yyzz_0[i] = ta_y_yyzz_0[i] * fe_0 - ta_y_yyzz_1[i] * fe_0 + ta_xy_yyzz_0[i] * pa_x[i] - ta_xy_yyzz_1[i] * pc_x[i];

        ta_xxy_yzzz_0[i] = ta_y_yzzz_0[i] * fe_0 - ta_y_yzzz_1[i] * fe_0 + ta_xy_yzzz_0[i] * pa_x[i] - ta_xy_yzzz_1[i] * pc_x[i];

        ta_xxy_zzzz_0[i] = ta_xx_zzzz_0[i] * pa_y[i] - ta_xx_zzzz_1[i] * pc_y[i];
    }

    // Set up 30-45 components of targeted buffer : FG

    auto ta_xxz_xxxx_0 = pbuffer.data(idx_npot_0_fg + 30);

    auto ta_xxz_xxxy_0 = pbuffer.data(idx_npot_0_fg + 31);

    auto ta_xxz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 32);

    auto ta_xxz_xxyy_0 = pbuffer.data(idx_npot_0_fg + 33);

    auto ta_xxz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 34);

    auto ta_xxz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 35);

    auto ta_xxz_xyyy_0 = pbuffer.data(idx_npot_0_fg + 36);

    auto ta_xxz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 37);

    auto ta_xxz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 38);

    auto ta_xxz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 39);

    auto ta_xxz_yyyy_0 = pbuffer.data(idx_npot_0_fg + 40);

    auto ta_xxz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 41);

    auto ta_xxz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 42);

    auto ta_xxz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 43);

    auto ta_xxz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 44);

#pragma omp simd aligned(pa_x,              \
                             pa_z,          \
                             pc_x,          \
                             pc_z,          \
                             ta_xx_xxx_0,   \
                             ta_xx_xxx_1,   \
                             ta_xx_xxxx_0,  \
                             ta_xx_xxxx_1,  \
                             ta_xx_xxxy_0,  \
                             ta_xx_xxxy_1,  \
                             ta_xx_xxxz_0,  \
                             ta_xx_xxxz_1,  \
                             ta_xx_xxy_0,   \
                             ta_xx_xxy_1,   \
                             ta_xx_xxyy_0,  \
                             ta_xx_xxyy_1,  \
                             ta_xx_xxyz_0,  \
                             ta_xx_xxyz_1,  \
                             ta_xx_xxz_0,   \
                             ta_xx_xxz_1,   \
                             ta_xx_xxzz_0,  \
                             ta_xx_xxzz_1,  \
                             ta_xx_xyy_0,   \
                             ta_xx_xyy_1,   \
                             ta_xx_xyyy_0,  \
                             ta_xx_xyyy_1,  \
                             ta_xx_xyyz_0,  \
                             ta_xx_xyyz_1,  \
                             ta_xx_xyz_0,   \
                             ta_xx_xyz_1,   \
                             ta_xx_xyzz_0,  \
                             ta_xx_xyzz_1,  \
                             ta_xx_xzz_0,   \
                             ta_xx_xzz_1,   \
                             ta_xx_xzzz_0,  \
                             ta_xx_xzzz_1,  \
                             ta_xx_yyyy_0,  \
                             ta_xx_yyyy_1,  \
                             ta_xxz_xxxx_0, \
                             ta_xxz_xxxy_0, \
                             ta_xxz_xxxz_0, \
                             ta_xxz_xxyy_0, \
                             ta_xxz_xxyz_0, \
                             ta_xxz_xxzz_0, \
                             ta_xxz_xyyy_0, \
                             ta_xxz_xyyz_0, \
                             ta_xxz_xyzz_0, \
                             ta_xxz_xzzz_0, \
                             ta_xxz_yyyy_0, \
                             ta_xxz_yyyz_0, \
                             ta_xxz_yyzz_0, \
                             ta_xxz_yzzz_0, \
                             ta_xxz_zzzz_0, \
                             ta_xz_yyyz_0,  \
                             ta_xz_yyyz_1,  \
                             ta_xz_yyzz_0,  \
                             ta_xz_yyzz_1,  \
                             ta_xz_yzzz_0,  \
                             ta_xz_yzzz_1,  \
                             ta_xz_zzzz_0,  \
                             ta_xz_zzzz_1,  \
                             ta_z_yyyz_0,   \
                             ta_z_yyyz_1,   \
                             ta_z_yyzz_0,   \
                             ta_z_yyzz_1,   \
                             ta_z_yzzz_0,   \
                             ta_z_yzzz_1,   \
                             ta_z_zzzz_0,   \
                             ta_z_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xxz_xxxx_0[i] = ta_xx_xxxx_0[i] * pa_z[i] - ta_xx_xxxx_1[i] * pc_z[i];

        ta_xxz_xxxy_0[i] = ta_xx_xxxy_0[i] * pa_z[i] - ta_xx_xxxy_1[i] * pc_z[i];

        ta_xxz_xxxz_0[i] = ta_xx_xxx_0[i] * fe_0 - ta_xx_xxx_1[i] * fe_0 + ta_xx_xxxz_0[i] * pa_z[i] - ta_xx_xxxz_1[i] * pc_z[i];

        ta_xxz_xxyy_0[i] = ta_xx_xxyy_0[i] * pa_z[i] - ta_xx_xxyy_1[i] * pc_z[i];

        ta_xxz_xxyz_0[i] = ta_xx_xxy_0[i] * fe_0 - ta_xx_xxy_1[i] * fe_0 + ta_xx_xxyz_0[i] * pa_z[i] - ta_xx_xxyz_1[i] * pc_z[i];

        ta_xxz_xxzz_0[i] = 2.0 * ta_xx_xxz_0[i] * fe_0 - 2.0 * ta_xx_xxz_1[i] * fe_0 + ta_xx_xxzz_0[i] * pa_z[i] - ta_xx_xxzz_1[i] * pc_z[i];

        ta_xxz_xyyy_0[i] = ta_xx_xyyy_0[i] * pa_z[i] - ta_xx_xyyy_1[i] * pc_z[i];

        ta_xxz_xyyz_0[i] = ta_xx_xyy_0[i] * fe_0 - ta_xx_xyy_1[i] * fe_0 + ta_xx_xyyz_0[i] * pa_z[i] - ta_xx_xyyz_1[i] * pc_z[i];

        ta_xxz_xyzz_0[i] = 2.0 * ta_xx_xyz_0[i] * fe_0 - 2.0 * ta_xx_xyz_1[i] * fe_0 + ta_xx_xyzz_0[i] * pa_z[i] - ta_xx_xyzz_1[i] * pc_z[i];

        ta_xxz_xzzz_0[i] = 3.0 * ta_xx_xzz_0[i] * fe_0 - 3.0 * ta_xx_xzz_1[i] * fe_0 + ta_xx_xzzz_0[i] * pa_z[i] - ta_xx_xzzz_1[i] * pc_z[i];

        ta_xxz_yyyy_0[i] = ta_xx_yyyy_0[i] * pa_z[i] - ta_xx_yyyy_1[i] * pc_z[i];

        ta_xxz_yyyz_0[i] = ta_z_yyyz_0[i] * fe_0 - ta_z_yyyz_1[i] * fe_0 + ta_xz_yyyz_0[i] * pa_x[i] - ta_xz_yyyz_1[i] * pc_x[i];

        ta_xxz_yyzz_0[i] = ta_z_yyzz_0[i] * fe_0 - ta_z_yyzz_1[i] * fe_0 + ta_xz_yyzz_0[i] * pa_x[i] - ta_xz_yyzz_1[i] * pc_x[i];

        ta_xxz_yzzz_0[i] = ta_z_yzzz_0[i] * fe_0 - ta_z_yzzz_1[i] * fe_0 + ta_xz_yzzz_0[i] * pa_x[i] - ta_xz_yzzz_1[i] * pc_x[i];

        ta_xxz_zzzz_0[i] = ta_z_zzzz_0[i] * fe_0 - ta_z_zzzz_1[i] * fe_0 + ta_xz_zzzz_0[i] * pa_x[i] - ta_xz_zzzz_1[i] * pc_x[i];
    }

    // Set up 45-60 components of targeted buffer : FG

    auto ta_xyy_xxxx_0 = pbuffer.data(idx_npot_0_fg + 45);

    auto ta_xyy_xxxy_0 = pbuffer.data(idx_npot_0_fg + 46);

    auto ta_xyy_xxxz_0 = pbuffer.data(idx_npot_0_fg + 47);

    auto ta_xyy_xxyy_0 = pbuffer.data(idx_npot_0_fg + 48);

    auto ta_xyy_xxyz_0 = pbuffer.data(idx_npot_0_fg + 49);

    auto ta_xyy_xxzz_0 = pbuffer.data(idx_npot_0_fg + 50);

    auto ta_xyy_xyyy_0 = pbuffer.data(idx_npot_0_fg + 51);

    auto ta_xyy_xyyz_0 = pbuffer.data(idx_npot_0_fg + 52);

    auto ta_xyy_xyzz_0 = pbuffer.data(idx_npot_0_fg + 53);

    auto ta_xyy_xzzz_0 = pbuffer.data(idx_npot_0_fg + 54);

    auto ta_xyy_yyyy_0 = pbuffer.data(idx_npot_0_fg + 55);

    auto ta_xyy_yyyz_0 = pbuffer.data(idx_npot_0_fg + 56);

    auto ta_xyy_yyzz_0 = pbuffer.data(idx_npot_0_fg + 57);

    auto ta_xyy_yzzz_0 = pbuffer.data(idx_npot_0_fg + 58);

    auto ta_xyy_zzzz_0 = pbuffer.data(idx_npot_0_fg + 59);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta_xyy_xxxx_0, \
                             ta_xyy_xxxy_0, \
                             ta_xyy_xxxz_0, \
                             ta_xyy_xxyy_0, \
                             ta_xyy_xxyz_0, \
                             ta_xyy_xxzz_0, \
                             ta_xyy_xyyy_0, \
                             ta_xyy_xyyz_0, \
                             ta_xyy_xyzz_0, \
                             ta_xyy_xzzz_0, \
                             ta_xyy_yyyy_0, \
                             ta_xyy_yyyz_0, \
                             ta_xyy_yyzz_0, \
                             ta_xyy_yzzz_0, \
                             ta_xyy_zzzz_0, \
                             ta_yy_xxx_0,   \
                             ta_yy_xxx_1,   \
                             ta_yy_xxxx_0,  \
                             ta_yy_xxxx_1,  \
                             ta_yy_xxxy_0,  \
                             ta_yy_xxxy_1,  \
                             ta_yy_xxxz_0,  \
                             ta_yy_xxxz_1,  \
                             ta_yy_xxy_0,   \
                             ta_yy_xxy_1,   \
                             ta_yy_xxyy_0,  \
                             ta_yy_xxyy_1,  \
                             ta_yy_xxyz_0,  \
                             ta_yy_xxyz_1,  \
                             ta_yy_xxz_0,   \
                             ta_yy_xxz_1,   \
                             ta_yy_xxzz_0,  \
                             ta_yy_xxzz_1,  \
                             ta_yy_xyy_0,   \
                             ta_yy_xyy_1,   \
                             ta_yy_xyyy_0,  \
                             ta_yy_xyyy_1,  \
                             ta_yy_xyyz_0,  \
                             ta_yy_xyyz_1,  \
                             ta_yy_xyz_0,   \
                             ta_yy_xyz_1,   \
                             ta_yy_xyzz_0,  \
                             ta_yy_xyzz_1,  \
                             ta_yy_xzz_0,   \
                             ta_yy_xzz_1,   \
                             ta_yy_xzzz_0,  \
                             ta_yy_xzzz_1,  \
                             ta_yy_yyy_0,   \
                             ta_yy_yyy_1,   \
                             ta_yy_yyyy_0,  \
                             ta_yy_yyyy_1,  \
                             ta_yy_yyyz_0,  \
                             ta_yy_yyyz_1,  \
                             ta_yy_yyz_0,   \
                             ta_yy_yyz_1,   \
                             ta_yy_yyzz_0,  \
                             ta_yy_yyzz_1,  \
                             ta_yy_yzz_0,   \
                             ta_yy_yzz_1,   \
                             ta_yy_yzzz_0,  \
                             ta_yy_yzzz_1,  \
                             ta_yy_zzz_0,   \
                             ta_yy_zzz_1,   \
                             ta_yy_zzzz_0,  \
                             ta_yy_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyy_xxxx_0[i] = 4.0 * ta_yy_xxx_0[i] * fe_0 - 4.0 * ta_yy_xxx_1[i] * fe_0 + ta_yy_xxxx_0[i] * pa_x[i] - ta_yy_xxxx_1[i] * pc_x[i];

        ta_xyy_xxxy_0[i] = 3.0 * ta_yy_xxy_0[i] * fe_0 - 3.0 * ta_yy_xxy_1[i] * fe_0 + ta_yy_xxxy_0[i] * pa_x[i] - ta_yy_xxxy_1[i] * pc_x[i];

        ta_xyy_xxxz_0[i] = 3.0 * ta_yy_xxz_0[i] * fe_0 - 3.0 * ta_yy_xxz_1[i] * fe_0 + ta_yy_xxxz_0[i] * pa_x[i] - ta_yy_xxxz_1[i] * pc_x[i];

        ta_xyy_xxyy_0[i] = 2.0 * ta_yy_xyy_0[i] * fe_0 - 2.0 * ta_yy_xyy_1[i] * fe_0 + ta_yy_xxyy_0[i] * pa_x[i] - ta_yy_xxyy_1[i] * pc_x[i];

        ta_xyy_xxyz_0[i] = 2.0 * ta_yy_xyz_0[i] * fe_0 - 2.0 * ta_yy_xyz_1[i] * fe_0 + ta_yy_xxyz_0[i] * pa_x[i] - ta_yy_xxyz_1[i] * pc_x[i];

        ta_xyy_xxzz_0[i] = 2.0 * ta_yy_xzz_0[i] * fe_0 - 2.0 * ta_yy_xzz_1[i] * fe_0 + ta_yy_xxzz_0[i] * pa_x[i] - ta_yy_xxzz_1[i] * pc_x[i];

        ta_xyy_xyyy_0[i] = ta_yy_yyy_0[i] * fe_0 - ta_yy_yyy_1[i] * fe_0 + ta_yy_xyyy_0[i] * pa_x[i] - ta_yy_xyyy_1[i] * pc_x[i];

        ta_xyy_xyyz_0[i] = ta_yy_yyz_0[i] * fe_0 - ta_yy_yyz_1[i] * fe_0 + ta_yy_xyyz_0[i] * pa_x[i] - ta_yy_xyyz_1[i] * pc_x[i];

        ta_xyy_xyzz_0[i] = ta_yy_yzz_0[i] * fe_0 - ta_yy_yzz_1[i] * fe_0 + ta_yy_xyzz_0[i] * pa_x[i] - ta_yy_xyzz_1[i] * pc_x[i];

        ta_xyy_xzzz_0[i] = ta_yy_zzz_0[i] * fe_0 - ta_yy_zzz_1[i] * fe_0 + ta_yy_xzzz_0[i] * pa_x[i] - ta_yy_xzzz_1[i] * pc_x[i];

        ta_xyy_yyyy_0[i] = ta_yy_yyyy_0[i] * pa_x[i] - ta_yy_yyyy_1[i] * pc_x[i];

        ta_xyy_yyyz_0[i] = ta_yy_yyyz_0[i] * pa_x[i] - ta_yy_yyyz_1[i] * pc_x[i];

        ta_xyy_yyzz_0[i] = ta_yy_yyzz_0[i] * pa_x[i] - ta_yy_yyzz_1[i] * pc_x[i];

        ta_xyy_yzzz_0[i] = ta_yy_yzzz_0[i] * pa_x[i] - ta_yy_yzzz_1[i] * pc_x[i];

        ta_xyy_zzzz_0[i] = ta_yy_zzzz_0[i] * pa_x[i] - ta_yy_zzzz_1[i] * pc_x[i];
    }

    // Set up 60-75 components of targeted buffer : FG

    auto ta_xyz_xxxx_0 = pbuffer.data(idx_npot_0_fg + 60);

    auto ta_xyz_xxxy_0 = pbuffer.data(idx_npot_0_fg + 61);

    auto ta_xyz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 62);

    auto ta_xyz_xxyy_0 = pbuffer.data(idx_npot_0_fg + 63);

    auto ta_xyz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 64);

    auto ta_xyz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 65);

    auto ta_xyz_xyyy_0 = pbuffer.data(idx_npot_0_fg + 66);

    auto ta_xyz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 67);

    auto ta_xyz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 68);

    auto ta_xyz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 69);

    auto ta_xyz_yyyy_0 = pbuffer.data(idx_npot_0_fg + 70);

    auto ta_xyz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 71);

    auto ta_xyz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 72);

    auto ta_xyz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 73);

    auto ta_xyz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 74);

#pragma omp simd aligned(pa_x,              \
                             pa_y,          \
                             pa_z,          \
                             pc_x,          \
                             pc_y,          \
                             pc_z,          \
                             ta_xy_xxxy_0,  \
                             ta_xy_xxxy_1,  \
                             ta_xy_xxyy_0,  \
                             ta_xy_xxyy_1,  \
                             ta_xy_xyyy_0,  \
                             ta_xy_xyyy_1,  \
                             ta_xyz_xxxx_0, \
                             ta_xyz_xxxy_0, \
                             ta_xyz_xxxz_0, \
                             ta_xyz_xxyy_0, \
                             ta_xyz_xxyz_0, \
                             ta_xyz_xxzz_0, \
                             ta_xyz_xyyy_0, \
                             ta_xyz_xyyz_0, \
                             ta_xyz_xyzz_0, \
                             ta_xyz_xzzz_0, \
                             ta_xyz_yyyy_0, \
                             ta_xyz_yyyz_0, \
                             ta_xyz_yyzz_0, \
                             ta_xyz_yzzz_0, \
                             ta_xyz_zzzz_0, \
                             ta_xz_xxxx_0,  \
                             ta_xz_xxxx_1,  \
                             ta_xz_xxxz_0,  \
                             ta_xz_xxxz_1,  \
                             ta_xz_xxzz_0,  \
                             ta_xz_xxzz_1,  \
                             ta_xz_xzzz_0,  \
                             ta_xz_xzzz_1,  \
                             ta_yz_xxyz_0,  \
                             ta_yz_xxyz_1,  \
                             ta_yz_xyyz_0,  \
                             ta_yz_xyyz_1,  \
                             ta_yz_xyz_0,   \
                             ta_yz_xyz_1,   \
                             ta_yz_xyzz_0,  \
                             ta_yz_xyzz_1,  \
                             ta_yz_yyyy_0,  \
                             ta_yz_yyyy_1,  \
                             ta_yz_yyyz_0,  \
                             ta_yz_yyyz_1,  \
                             ta_yz_yyz_0,   \
                             ta_yz_yyz_1,   \
                             ta_yz_yyzz_0,  \
                             ta_yz_yyzz_1,  \
                             ta_yz_yzz_0,   \
                             ta_yz_yzz_1,   \
                             ta_yz_yzzz_0,  \
                             ta_yz_yzzz_1,  \
                             ta_yz_zzzz_0,  \
                             ta_yz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xyz_xxxx_0[i] = ta_xz_xxxx_0[i] * pa_y[i] - ta_xz_xxxx_1[i] * pc_y[i];

        ta_xyz_xxxy_0[i] = ta_xy_xxxy_0[i] * pa_z[i] - ta_xy_xxxy_1[i] * pc_z[i];

        ta_xyz_xxxz_0[i] = ta_xz_xxxz_0[i] * pa_y[i] - ta_xz_xxxz_1[i] * pc_y[i];

        ta_xyz_xxyy_0[i] = ta_xy_xxyy_0[i] * pa_z[i] - ta_xy_xxyy_1[i] * pc_z[i];

        ta_xyz_xxyz_0[i] = 2.0 * ta_yz_xyz_0[i] * fe_0 - 2.0 * ta_yz_xyz_1[i] * fe_0 + ta_yz_xxyz_0[i] * pa_x[i] - ta_yz_xxyz_1[i] * pc_x[i];

        ta_xyz_xxzz_0[i] = ta_xz_xxzz_0[i] * pa_y[i] - ta_xz_xxzz_1[i] * pc_y[i];

        ta_xyz_xyyy_0[i] = ta_xy_xyyy_0[i] * pa_z[i] - ta_xy_xyyy_1[i] * pc_z[i];

        ta_xyz_xyyz_0[i] = ta_yz_yyz_0[i] * fe_0 - ta_yz_yyz_1[i] * fe_0 + ta_yz_xyyz_0[i] * pa_x[i] - ta_yz_xyyz_1[i] * pc_x[i];

        ta_xyz_xyzz_0[i] = ta_yz_yzz_0[i] * fe_0 - ta_yz_yzz_1[i] * fe_0 + ta_yz_xyzz_0[i] * pa_x[i] - ta_yz_xyzz_1[i] * pc_x[i];

        ta_xyz_xzzz_0[i] = ta_xz_xzzz_0[i] * pa_y[i] - ta_xz_xzzz_1[i] * pc_y[i];

        ta_xyz_yyyy_0[i] = ta_yz_yyyy_0[i] * pa_x[i] - ta_yz_yyyy_1[i] * pc_x[i];

        ta_xyz_yyyz_0[i] = ta_yz_yyyz_0[i] * pa_x[i] - ta_yz_yyyz_1[i] * pc_x[i];

        ta_xyz_yyzz_0[i] = ta_yz_yyzz_0[i] * pa_x[i] - ta_yz_yyzz_1[i] * pc_x[i];

        ta_xyz_yzzz_0[i] = ta_yz_yzzz_0[i] * pa_x[i] - ta_yz_yzzz_1[i] * pc_x[i];

        ta_xyz_zzzz_0[i] = ta_yz_zzzz_0[i] * pa_x[i] - ta_yz_zzzz_1[i] * pc_x[i];
    }

    // Set up 75-90 components of targeted buffer : FG

    auto ta_xzz_xxxx_0 = pbuffer.data(idx_npot_0_fg + 75);

    auto ta_xzz_xxxy_0 = pbuffer.data(idx_npot_0_fg + 76);

    auto ta_xzz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 77);

    auto ta_xzz_xxyy_0 = pbuffer.data(idx_npot_0_fg + 78);

    auto ta_xzz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 79);

    auto ta_xzz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 80);

    auto ta_xzz_xyyy_0 = pbuffer.data(idx_npot_0_fg + 81);

    auto ta_xzz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 82);

    auto ta_xzz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 83);

    auto ta_xzz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 84);

    auto ta_xzz_yyyy_0 = pbuffer.data(idx_npot_0_fg + 85);

    auto ta_xzz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 86);

    auto ta_xzz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 87);

    auto ta_xzz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 88);

    auto ta_xzz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 89);

#pragma omp simd aligned(pa_x,              \
                             pc_x,          \
                             ta_xzz_xxxx_0, \
                             ta_xzz_xxxy_0, \
                             ta_xzz_xxxz_0, \
                             ta_xzz_xxyy_0, \
                             ta_xzz_xxyz_0, \
                             ta_xzz_xxzz_0, \
                             ta_xzz_xyyy_0, \
                             ta_xzz_xyyz_0, \
                             ta_xzz_xyzz_0, \
                             ta_xzz_xzzz_0, \
                             ta_xzz_yyyy_0, \
                             ta_xzz_yyyz_0, \
                             ta_xzz_yyzz_0, \
                             ta_xzz_yzzz_0, \
                             ta_xzz_zzzz_0, \
                             ta_zz_xxx_0,   \
                             ta_zz_xxx_1,   \
                             ta_zz_xxxx_0,  \
                             ta_zz_xxxx_1,  \
                             ta_zz_xxxy_0,  \
                             ta_zz_xxxy_1,  \
                             ta_zz_xxxz_0,  \
                             ta_zz_xxxz_1,  \
                             ta_zz_xxy_0,   \
                             ta_zz_xxy_1,   \
                             ta_zz_xxyy_0,  \
                             ta_zz_xxyy_1,  \
                             ta_zz_xxyz_0,  \
                             ta_zz_xxyz_1,  \
                             ta_zz_xxz_0,   \
                             ta_zz_xxz_1,   \
                             ta_zz_xxzz_0,  \
                             ta_zz_xxzz_1,  \
                             ta_zz_xyy_0,   \
                             ta_zz_xyy_1,   \
                             ta_zz_xyyy_0,  \
                             ta_zz_xyyy_1,  \
                             ta_zz_xyyz_0,  \
                             ta_zz_xyyz_1,  \
                             ta_zz_xyz_0,   \
                             ta_zz_xyz_1,   \
                             ta_zz_xyzz_0,  \
                             ta_zz_xyzz_1,  \
                             ta_zz_xzz_0,   \
                             ta_zz_xzz_1,   \
                             ta_zz_xzzz_0,  \
                             ta_zz_xzzz_1,  \
                             ta_zz_yyy_0,   \
                             ta_zz_yyy_1,   \
                             ta_zz_yyyy_0,  \
                             ta_zz_yyyy_1,  \
                             ta_zz_yyyz_0,  \
                             ta_zz_yyyz_1,  \
                             ta_zz_yyz_0,   \
                             ta_zz_yyz_1,   \
                             ta_zz_yyzz_0,  \
                             ta_zz_yyzz_1,  \
                             ta_zz_yzz_0,   \
                             ta_zz_yzz_1,   \
                             ta_zz_yzzz_0,  \
                             ta_zz_yzzz_1,  \
                             ta_zz_zzz_0,   \
                             ta_zz_zzz_1,   \
                             ta_zz_zzzz_0,  \
                             ta_zz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_xzz_xxxx_0[i] = 4.0 * ta_zz_xxx_0[i] * fe_0 - 4.0 * ta_zz_xxx_1[i] * fe_0 + ta_zz_xxxx_0[i] * pa_x[i] - ta_zz_xxxx_1[i] * pc_x[i];

        ta_xzz_xxxy_0[i] = 3.0 * ta_zz_xxy_0[i] * fe_0 - 3.0 * ta_zz_xxy_1[i] * fe_0 + ta_zz_xxxy_0[i] * pa_x[i] - ta_zz_xxxy_1[i] * pc_x[i];

        ta_xzz_xxxz_0[i] = 3.0 * ta_zz_xxz_0[i] * fe_0 - 3.0 * ta_zz_xxz_1[i] * fe_0 + ta_zz_xxxz_0[i] * pa_x[i] - ta_zz_xxxz_1[i] * pc_x[i];

        ta_xzz_xxyy_0[i] = 2.0 * ta_zz_xyy_0[i] * fe_0 - 2.0 * ta_zz_xyy_1[i] * fe_0 + ta_zz_xxyy_0[i] * pa_x[i] - ta_zz_xxyy_1[i] * pc_x[i];

        ta_xzz_xxyz_0[i] = 2.0 * ta_zz_xyz_0[i] * fe_0 - 2.0 * ta_zz_xyz_1[i] * fe_0 + ta_zz_xxyz_0[i] * pa_x[i] - ta_zz_xxyz_1[i] * pc_x[i];

        ta_xzz_xxzz_0[i] = 2.0 * ta_zz_xzz_0[i] * fe_0 - 2.0 * ta_zz_xzz_1[i] * fe_0 + ta_zz_xxzz_0[i] * pa_x[i] - ta_zz_xxzz_1[i] * pc_x[i];

        ta_xzz_xyyy_0[i] = ta_zz_yyy_0[i] * fe_0 - ta_zz_yyy_1[i] * fe_0 + ta_zz_xyyy_0[i] * pa_x[i] - ta_zz_xyyy_1[i] * pc_x[i];

        ta_xzz_xyyz_0[i] = ta_zz_yyz_0[i] * fe_0 - ta_zz_yyz_1[i] * fe_0 + ta_zz_xyyz_0[i] * pa_x[i] - ta_zz_xyyz_1[i] * pc_x[i];

        ta_xzz_xyzz_0[i] = ta_zz_yzz_0[i] * fe_0 - ta_zz_yzz_1[i] * fe_0 + ta_zz_xyzz_0[i] * pa_x[i] - ta_zz_xyzz_1[i] * pc_x[i];

        ta_xzz_xzzz_0[i] = ta_zz_zzz_0[i] * fe_0 - ta_zz_zzz_1[i] * fe_0 + ta_zz_xzzz_0[i] * pa_x[i] - ta_zz_xzzz_1[i] * pc_x[i];

        ta_xzz_yyyy_0[i] = ta_zz_yyyy_0[i] * pa_x[i] - ta_zz_yyyy_1[i] * pc_x[i];

        ta_xzz_yyyz_0[i] = ta_zz_yyyz_0[i] * pa_x[i] - ta_zz_yyyz_1[i] * pc_x[i];

        ta_xzz_yyzz_0[i] = ta_zz_yyzz_0[i] * pa_x[i] - ta_zz_yyzz_1[i] * pc_x[i];

        ta_xzz_yzzz_0[i] = ta_zz_yzzz_0[i] * pa_x[i] - ta_zz_yzzz_1[i] * pc_x[i];

        ta_xzz_zzzz_0[i] = ta_zz_zzzz_0[i] * pa_x[i] - ta_zz_zzzz_1[i] * pc_x[i];
    }

    // Set up 90-105 components of targeted buffer : FG

    auto ta_yyy_xxxx_0 = pbuffer.data(idx_npot_0_fg + 90);

    auto ta_yyy_xxxy_0 = pbuffer.data(idx_npot_0_fg + 91);

    auto ta_yyy_xxxz_0 = pbuffer.data(idx_npot_0_fg + 92);

    auto ta_yyy_xxyy_0 = pbuffer.data(idx_npot_0_fg + 93);

    auto ta_yyy_xxyz_0 = pbuffer.data(idx_npot_0_fg + 94);

    auto ta_yyy_xxzz_0 = pbuffer.data(idx_npot_0_fg + 95);

    auto ta_yyy_xyyy_0 = pbuffer.data(idx_npot_0_fg + 96);

    auto ta_yyy_xyyz_0 = pbuffer.data(idx_npot_0_fg + 97);

    auto ta_yyy_xyzz_0 = pbuffer.data(idx_npot_0_fg + 98);

    auto ta_yyy_xzzz_0 = pbuffer.data(idx_npot_0_fg + 99);

    auto ta_yyy_yyyy_0 = pbuffer.data(idx_npot_0_fg + 100);

    auto ta_yyy_yyyz_0 = pbuffer.data(idx_npot_0_fg + 101);

    auto ta_yyy_yyzz_0 = pbuffer.data(idx_npot_0_fg + 102);

    auto ta_yyy_yzzz_0 = pbuffer.data(idx_npot_0_fg + 103);

    auto ta_yyy_zzzz_0 = pbuffer.data(idx_npot_0_fg + 104);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta_y_xxxx_0,   \
                             ta_y_xxxx_1,   \
                             ta_y_xxxy_0,   \
                             ta_y_xxxy_1,   \
                             ta_y_xxxz_0,   \
                             ta_y_xxxz_1,   \
                             ta_y_xxyy_0,   \
                             ta_y_xxyy_1,   \
                             ta_y_xxyz_0,   \
                             ta_y_xxyz_1,   \
                             ta_y_xxzz_0,   \
                             ta_y_xxzz_1,   \
                             ta_y_xyyy_0,   \
                             ta_y_xyyy_1,   \
                             ta_y_xyyz_0,   \
                             ta_y_xyyz_1,   \
                             ta_y_xyzz_0,   \
                             ta_y_xyzz_1,   \
                             ta_y_xzzz_0,   \
                             ta_y_xzzz_1,   \
                             ta_y_yyyy_0,   \
                             ta_y_yyyy_1,   \
                             ta_y_yyyz_0,   \
                             ta_y_yyyz_1,   \
                             ta_y_yyzz_0,   \
                             ta_y_yyzz_1,   \
                             ta_y_yzzz_0,   \
                             ta_y_yzzz_1,   \
                             ta_y_zzzz_0,   \
                             ta_y_zzzz_1,   \
                             ta_yy_xxx_0,   \
                             ta_yy_xxx_1,   \
                             ta_yy_xxxx_0,  \
                             ta_yy_xxxx_1,  \
                             ta_yy_xxxy_0,  \
                             ta_yy_xxxy_1,  \
                             ta_yy_xxxz_0,  \
                             ta_yy_xxxz_1,  \
                             ta_yy_xxy_0,   \
                             ta_yy_xxy_1,   \
                             ta_yy_xxyy_0,  \
                             ta_yy_xxyy_1,  \
                             ta_yy_xxyz_0,  \
                             ta_yy_xxyz_1,  \
                             ta_yy_xxz_0,   \
                             ta_yy_xxz_1,   \
                             ta_yy_xxzz_0,  \
                             ta_yy_xxzz_1,  \
                             ta_yy_xyy_0,   \
                             ta_yy_xyy_1,   \
                             ta_yy_xyyy_0,  \
                             ta_yy_xyyy_1,  \
                             ta_yy_xyyz_0,  \
                             ta_yy_xyyz_1,  \
                             ta_yy_xyz_0,   \
                             ta_yy_xyz_1,   \
                             ta_yy_xyzz_0,  \
                             ta_yy_xyzz_1,  \
                             ta_yy_xzz_0,   \
                             ta_yy_xzz_1,   \
                             ta_yy_xzzz_0,  \
                             ta_yy_xzzz_1,  \
                             ta_yy_yyy_0,   \
                             ta_yy_yyy_1,   \
                             ta_yy_yyyy_0,  \
                             ta_yy_yyyy_1,  \
                             ta_yy_yyyz_0,  \
                             ta_yy_yyyz_1,  \
                             ta_yy_yyz_0,   \
                             ta_yy_yyz_1,   \
                             ta_yy_yyzz_0,  \
                             ta_yy_yyzz_1,  \
                             ta_yy_yzz_0,   \
                             ta_yy_yzz_1,   \
                             ta_yy_yzzz_0,  \
                             ta_yy_yzzz_1,  \
                             ta_yy_zzz_0,   \
                             ta_yy_zzz_1,   \
                             ta_yy_zzzz_0,  \
                             ta_yy_zzzz_1,  \
                             ta_yyy_xxxx_0, \
                             ta_yyy_xxxy_0, \
                             ta_yyy_xxxz_0, \
                             ta_yyy_xxyy_0, \
                             ta_yyy_xxyz_0, \
                             ta_yyy_xxzz_0, \
                             ta_yyy_xyyy_0, \
                             ta_yyy_xyyz_0, \
                             ta_yyy_xyzz_0, \
                             ta_yyy_xzzz_0, \
                             ta_yyy_yyyy_0, \
                             ta_yyy_yyyz_0, \
                             ta_yyy_yyzz_0, \
                             ta_yyy_yzzz_0, \
                             ta_yyy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyy_xxxx_0[i] = 2.0 * ta_y_xxxx_0[i] * fe_0 - 2.0 * ta_y_xxxx_1[i] * fe_0 + ta_yy_xxxx_0[i] * pa_y[i] - ta_yy_xxxx_1[i] * pc_y[i];

        ta_yyy_xxxy_0[i] = 2.0 * ta_y_xxxy_0[i] * fe_0 - 2.0 * ta_y_xxxy_1[i] * fe_0 + ta_yy_xxx_0[i] * fe_0 - ta_yy_xxx_1[i] * fe_0 +
                           ta_yy_xxxy_0[i] * pa_y[i] - ta_yy_xxxy_1[i] * pc_y[i];

        ta_yyy_xxxz_0[i] = 2.0 * ta_y_xxxz_0[i] * fe_0 - 2.0 * ta_y_xxxz_1[i] * fe_0 + ta_yy_xxxz_0[i] * pa_y[i] - ta_yy_xxxz_1[i] * pc_y[i];

        ta_yyy_xxyy_0[i] = 2.0 * ta_y_xxyy_0[i] * fe_0 - 2.0 * ta_y_xxyy_1[i] * fe_0 + 2.0 * ta_yy_xxy_0[i] * fe_0 - 2.0 * ta_yy_xxy_1[i] * fe_0 +
                           ta_yy_xxyy_0[i] * pa_y[i] - ta_yy_xxyy_1[i] * pc_y[i];

        ta_yyy_xxyz_0[i] = 2.0 * ta_y_xxyz_0[i] * fe_0 - 2.0 * ta_y_xxyz_1[i] * fe_0 + ta_yy_xxz_0[i] * fe_0 - ta_yy_xxz_1[i] * fe_0 +
                           ta_yy_xxyz_0[i] * pa_y[i] - ta_yy_xxyz_1[i] * pc_y[i];

        ta_yyy_xxzz_0[i] = 2.0 * ta_y_xxzz_0[i] * fe_0 - 2.0 * ta_y_xxzz_1[i] * fe_0 + ta_yy_xxzz_0[i] * pa_y[i] - ta_yy_xxzz_1[i] * pc_y[i];

        ta_yyy_xyyy_0[i] = 2.0 * ta_y_xyyy_0[i] * fe_0 - 2.0 * ta_y_xyyy_1[i] * fe_0 + 3.0 * ta_yy_xyy_0[i] * fe_0 - 3.0 * ta_yy_xyy_1[i] * fe_0 +
                           ta_yy_xyyy_0[i] * pa_y[i] - ta_yy_xyyy_1[i] * pc_y[i];

        ta_yyy_xyyz_0[i] = 2.0 * ta_y_xyyz_0[i] * fe_0 - 2.0 * ta_y_xyyz_1[i] * fe_0 + 2.0 * ta_yy_xyz_0[i] * fe_0 - 2.0 * ta_yy_xyz_1[i] * fe_0 +
                           ta_yy_xyyz_0[i] * pa_y[i] - ta_yy_xyyz_1[i] * pc_y[i];

        ta_yyy_xyzz_0[i] = 2.0 * ta_y_xyzz_0[i] * fe_0 - 2.0 * ta_y_xyzz_1[i] * fe_0 + ta_yy_xzz_0[i] * fe_0 - ta_yy_xzz_1[i] * fe_0 +
                           ta_yy_xyzz_0[i] * pa_y[i] - ta_yy_xyzz_1[i] * pc_y[i];

        ta_yyy_xzzz_0[i] = 2.0 * ta_y_xzzz_0[i] * fe_0 - 2.0 * ta_y_xzzz_1[i] * fe_0 + ta_yy_xzzz_0[i] * pa_y[i] - ta_yy_xzzz_1[i] * pc_y[i];

        ta_yyy_yyyy_0[i] = 2.0 * ta_y_yyyy_0[i] * fe_0 - 2.0 * ta_y_yyyy_1[i] * fe_0 + 4.0 * ta_yy_yyy_0[i] * fe_0 - 4.0 * ta_yy_yyy_1[i] * fe_0 +
                           ta_yy_yyyy_0[i] * pa_y[i] - ta_yy_yyyy_1[i] * pc_y[i];

        ta_yyy_yyyz_0[i] = 2.0 * ta_y_yyyz_0[i] * fe_0 - 2.0 * ta_y_yyyz_1[i] * fe_0 + 3.0 * ta_yy_yyz_0[i] * fe_0 - 3.0 * ta_yy_yyz_1[i] * fe_0 +
                           ta_yy_yyyz_0[i] * pa_y[i] - ta_yy_yyyz_1[i] * pc_y[i];

        ta_yyy_yyzz_0[i] = 2.0 * ta_y_yyzz_0[i] * fe_0 - 2.0 * ta_y_yyzz_1[i] * fe_0 + 2.0 * ta_yy_yzz_0[i] * fe_0 - 2.0 * ta_yy_yzz_1[i] * fe_0 +
                           ta_yy_yyzz_0[i] * pa_y[i] - ta_yy_yyzz_1[i] * pc_y[i];

        ta_yyy_yzzz_0[i] = 2.0 * ta_y_yzzz_0[i] * fe_0 - 2.0 * ta_y_yzzz_1[i] * fe_0 + ta_yy_zzz_0[i] * fe_0 - ta_yy_zzz_1[i] * fe_0 +
                           ta_yy_yzzz_0[i] * pa_y[i] - ta_yy_yzzz_1[i] * pc_y[i];

        ta_yyy_zzzz_0[i] = 2.0 * ta_y_zzzz_0[i] * fe_0 - 2.0 * ta_y_zzzz_1[i] * fe_0 + ta_yy_zzzz_0[i] * pa_y[i] - ta_yy_zzzz_1[i] * pc_y[i];
    }

    // Set up 105-120 components of targeted buffer : FG

    auto ta_yyz_xxxx_0 = pbuffer.data(idx_npot_0_fg + 105);

    auto ta_yyz_xxxy_0 = pbuffer.data(idx_npot_0_fg + 106);

    auto ta_yyz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 107);

    auto ta_yyz_xxyy_0 = pbuffer.data(idx_npot_0_fg + 108);

    auto ta_yyz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 109);

    auto ta_yyz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 110);

    auto ta_yyz_xyyy_0 = pbuffer.data(idx_npot_0_fg + 111);

    auto ta_yyz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 112);

    auto ta_yyz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 113);

    auto ta_yyz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 114);

    auto ta_yyz_yyyy_0 = pbuffer.data(idx_npot_0_fg + 115);

    auto ta_yyz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 116);

    auto ta_yyz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 117);

    auto ta_yyz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 118);

    auto ta_yyz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 119);

#pragma omp simd aligned(pa_y,              \
                             pa_z,          \
                             pc_y,          \
                             pc_z,          \
                             ta_yy_xxxx_0,  \
                             ta_yy_xxxx_1,  \
                             ta_yy_xxxy_0,  \
                             ta_yy_xxxy_1,  \
                             ta_yy_xxy_0,   \
                             ta_yy_xxy_1,   \
                             ta_yy_xxyy_0,  \
                             ta_yy_xxyy_1,  \
                             ta_yy_xxyz_0,  \
                             ta_yy_xxyz_1,  \
                             ta_yy_xyy_0,   \
                             ta_yy_xyy_1,   \
                             ta_yy_xyyy_0,  \
                             ta_yy_xyyy_1,  \
                             ta_yy_xyyz_0,  \
                             ta_yy_xyyz_1,  \
                             ta_yy_xyz_0,   \
                             ta_yy_xyz_1,   \
                             ta_yy_xyzz_0,  \
                             ta_yy_xyzz_1,  \
                             ta_yy_yyy_0,   \
                             ta_yy_yyy_1,   \
                             ta_yy_yyyy_0,  \
                             ta_yy_yyyy_1,  \
                             ta_yy_yyyz_0,  \
                             ta_yy_yyyz_1,  \
                             ta_yy_yyz_0,   \
                             ta_yy_yyz_1,   \
                             ta_yy_yyzz_0,  \
                             ta_yy_yyzz_1,  \
                             ta_yy_yzz_0,   \
                             ta_yy_yzz_1,   \
                             ta_yy_yzzz_0,  \
                             ta_yy_yzzz_1,  \
                             ta_yyz_xxxx_0, \
                             ta_yyz_xxxy_0, \
                             ta_yyz_xxxz_0, \
                             ta_yyz_xxyy_0, \
                             ta_yyz_xxyz_0, \
                             ta_yyz_xxzz_0, \
                             ta_yyz_xyyy_0, \
                             ta_yyz_xyyz_0, \
                             ta_yyz_xyzz_0, \
                             ta_yyz_xzzz_0, \
                             ta_yyz_yyyy_0, \
                             ta_yyz_yyyz_0, \
                             ta_yyz_yyzz_0, \
                             ta_yyz_yzzz_0, \
                             ta_yyz_zzzz_0, \
                             ta_yz_xxxz_0,  \
                             ta_yz_xxxz_1,  \
                             ta_yz_xxzz_0,  \
                             ta_yz_xxzz_1,  \
                             ta_yz_xzzz_0,  \
                             ta_yz_xzzz_1,  \
                             ta_yz_zzzz_0,  \
                             ta_yz_zzzz_1,  \
                             ta_z_xxxz_0,   \
                             ta_z_xxxz_1,   \
                             ta_z_xxzz_0,   \
                             ta_z_xxzz_1,   \
                             ta_z_xzzz_0,   \
                             ta_z_xzzz_1,   \
                             ta_z_zzzz_0,   \
                             ta_z_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yyz_xxxx_0[i] = ta_yy_xxxx_0[i] * pa_z[i] - ta_yy_xxxx_1[i] * pc_z[i];

        ta_yyz_xxxy_0[i] = ta_yy_xxxy_0[i] * pa_z[i] - ta_yy_xxxy_1[i] * pc_z[i];

        ta_yyz_xxxz_0[i] = ta_z_xxxz_0[i] * fe_0 - ta_z_xxxz_1[i] * fe_0 + ta_yz_xxxz_0[i] * pa_y[i] - ta_yz_xxxz_1[i] * pc_y[i];

        ta_yyz_xxyy_0[i] = ta_yy_xxyy_0[i] * pa_z[i] - ta_yy_xxyy_1[i] * pc_z[i];

        ta_yyz_xxyz_0[i] = ta_yy_xxy_0[i] * fe_0 - ta_yy_xxy_1[i] * fe_0 + ta_yy_xxyz_0[i] * pa_z[i] - ta_yy_xxyz_1[i] * pc_z[i];

        ta_yyz_xxzz_0[i] = ta_z_xxzz_0[i] * fe_0 - ta_z_xxzz_1[i] * fe_0 + ta_yz_xxzz_0[i] * pa_y[i] - ta_yz_xxzz_1[i] * pc_y[i];

        ta_yyz_xyyy_0[i] = ta_yy_xyyy_0[i] * pa_z[i] - ta_yy_xyyy_1[i] * pc_z[i];

        ta_yyz_xyyz_0[i] = ta_yy_xyy_0[i] * fe_0 - ta_yy_xyy_1[i] * fe_0 + ta_yy_xyyz_0[i] * pa_z[i] - ta_yy_xyyz_1[i] * pc_z[i];

        ta_yyz_xyzz_0[i] = 2.0 * ta_yy_xyz_0[i] * fe_0 - 2.0 * ta_yy_xyz_1[i] * fe_0 + ta_yy_xyzz_0[i] * pa_z[i] - ta_yy_xyzz_1[i] * pc_z[i];

        ta_yyz_xzzz_0[i] = ta_z_xzzz_0[i] * fe_0 - ta_z_xzzz_1[i] * fe_0 + ta_yz_xzzz_0[i] * pa_y[i] - ta_yz_xzzz_1[i] * pc_y[i];

        ta_yyz_yyyy_0[i] = ta_yy_yyyy_0[i] * pa_z[i] - ta_yy_yyyy_1[i] * pc_z[i];

        ta_yyz_yyyz_0[i] = ta_yy_yyy_0[i] * fe_0 - ta_yy_yyy_1[i] * fe_0 + ta_yy_yyyz_0[i] * pa_z[i] - ta_yy_yyyz_1[i] * pc_z[i];

        ta_yyz_yyzz_0[i] = 2.0 * ta_yy_yyz_0[i] * fe_0 - 2.0 * ta_yy_yyz_1[i] * fe_0 + ta_yy_yyzz_0[i] * pa_z[i] - ta_yy_yyzz_1[i] * pc_z[i];

        ta_yyz_yzzz_0[i] = 3.0 * ta_yy_yzz_0[i] * fe_0 - 3.0 * ta_yy_yzz_1[i] * fe_0 + ta_yy_yzzz_0[i] * pa_z[i] - ta_yy_yzzz_1[i] * pc_z[i];

        ta_yyz_zzzz_0[i] = ta_z_zzzz_0[i] * fe_0 - ta_z_zzzz_1[i] * fe_0 + ta_yz_zzzz_0[i] * pa_y[i] - ta_yz_zzzz_1[i] * pc_y[i];
    }

    // Set up 120-135 components of targeted buffer : FG

    auto ta_yzz_xxxx_0 = pbuffer.data(idx_npot_0_fg + 120);

    auto ta_yzz_xxxy_0 = pbuffer.data(idx_npot_0_fg + 121);

    auto ta_yzz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 122);

    auto ta_yzz_xxyy_0 = pbuffer.data(idx_npot_0_fg + 123);

    auto ta_yzz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 124);

    auto ta_yzz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 125);

    auto ta_yzz_xyyy_0 = pbuffer.data(idx_npot_0_fg + 126);

    auto ta_yzz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 127);

    auto ta_yzz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 128);

    auto ta_yzz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 129);

    auto ta_yzz_yyyy_0 = pbuffer.data(idx_npot_0_fg + 130);

    auto ta_yzz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 131);

    auto ta_yzz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 132);

    auto ta_yzz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 133);

    auto ta_yzz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 134);

#pragma omp simd aligned(pa_y,              \
                             pc_y,          \
                             ta_yzz_xxxx_0, \
                             ta_yzz_xxxy_0, \
                             ta_yzz_xxxz_0, \
                             ta_yzz_xxyy_0, \
                             ta_yzz_xxyz_0, \
                             ta_yzz_xxzz_0, \
                             ta_yzz_xyyy_0, \
                             ta_yzz_xyyz_0, \
                             ta_yzz_xyzz_0, \
                             ta_yzz_xzzz_0, \
                             ta_yzz_yyyy_0, \
                             ta_yzz_yyyz_0, \
                             ta_yzz_yyzz_0, \
                             ta_yzz_yzzz_0, \
                             ta_yzz_zzzz_0, \
                             ta_zz_xxx_0,   \
                             ta_zz_xxx_1,   \
                             ta_zz_xxxx_0,  \
                             ta_zz_xxxx_1,  \
                             ta_zz_xxxy_0,  \
                             ta_zz_xxxy_1,  \
                             ta_zz_xxxz_0,  \
                             ta_zz_xxxz_1,  \
                             ta_zz_xxy_0,   \
                             ta_zz_xxy_1,   \
                             ta_zz_xxyy_0,  \
                             ta_zz_xxyy_1,  \
                             ta_zz_xxyz_0,  \
                             ta_zz_xxyz_1,  \
                             ta_zz_xxz_0,   \
                             ta_zz_xxz_1,   \
                             ta_zz_xxzz_0,  \
                             ta_zz_xxzz_1,  \
                             ta_zz_xyy_0,   \
                             ta_zz_xyy_1,   \
                             ta_zz_xyyy_0,  \
                             ta_zz_xyyy_1,  \
                             ta_zz_xyyz_0,  \
                             ta_zz_xyyz_1,  \
                             ta_zz_xyz_0,   \
                             ta_zz_xyz_1,   \
                             ta_zz_xyzz_0,  \
                             ta_zz_xyzz_1,  \
                             ta_zz_xzz_0,   \
                             ta_zz_xzz_1,   \
                             ta_zz_xzzz_0,  \
                             ta_zz_xzzz_1,  \
                             ta_zz_yyy_0,   \
                             ta_zz_yyy_1,   \
                             ta_zz_yyyy_0,  \
                             ta_zz_yyyy_1,  \
                             ta_zz_yyyz_0,  \
                             ta_zz_yyyz_1,  \
                             ta_zz_yyz_0,   \
                             ta_zz_yyz_1,   \
                             ta_zz_yyzz_0,  \
                             ta_zz_yyzz_1,  \
                             ta_zz_yzz_0,   \
                             ta_zz_yzz_1,   \
                             ta_zz_yzzz_0,  \
                             ta_zz_yzzz_1,  \
                             ta_zz_zzz_0,   \
                             ta_zz_zzz_1,   \
                             ta_zz_zzzz_0,  \
                             ta_zz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_yzz_xxxx_0[i] = ta_zz_xxxx_0[i] * pa_y[i] - ta_zz_xxxx_1[i] * pc_y[i];

        ta_yzz_xxxy_0[i] = ta_zz_xxx_0[i] * fe_0 - ta_zz_xxx_1[i] * fe_0 + ta_zz_xxxy_0[i] * pa_y[i] - ta_zz_xxxy_1[i] * pc_y[i];

        ta_yzz_xxxz_0[i] = ta_zz_xxxz_0[i] * pa_y[i] - ta_zz_xxxz_1[i] * pc_y[i];

        ta_yzz_xxyy_0[i] = 2.0 * ta_zz_xxy_0[i] * fe_0 - 2.0 * ta_zz_xxy_1[i] * fe_0 + ta_zz_xxyy_0[i] * pa_y[i] - ta_zz_xxyy_1[i] * pc_y[i];

        ta_yzz_xxyz_0[i] = ta_zz_xxz_0[i] * fe_0 - ta_zz_xxz_1[i] * fe_0 + ta_zz_xxyz_0[i] * pa_y[i] - ta_zz_xxyz_1[i] * pc_y[i];

        ta_yzz_xxzz_0[i] = ta_zz_xxzz_0[i] * pa_y[i] - ta_zz_xxzz_1[i] * pc_y[i];

        ta_yzz_xyyy_0[i] = 3.0 * ta_zz_xyy_0[i] * fe_0 - 3.0 * ta_zz_xyy_1[i] * fe_0 + ta_zz_xyyy_0[i] * pa_y[i] - ta_zz_xyyy_1[i] * pc_y[i];

        ta_yzz_xyyz_0[i] = 2.0 * ta_zz_xyz_0[i] * fe_0 - 2.0 * ta_zz_xyz_1[i] * fe_0 + ta_zz_xyyz_0[i] * pa_y[i] - ta_zz_xyyz_1[i] * pc_y[i];

        ta_yzz_xyzz_0[i] = ta_zz_xzz_0[i] * fe_0 - ta_zz_xzz_1[i] * fe_0 + ta_zz_xyzz_0[i] * pa_y[i] - ta_zz_xyzz_1[i] * pc_y[i];

        ta_yzz_xzzz_0[i] = ta_zz_xzzz_0[i] * pa_y[i] - ta_zz_xzzz_1[i] * pc_y[i];

        ta_yzz_yyyy_0[i] = 4.0 * ta_zz_yyy_0[i] * fe_0 - 4.0 * ta_zz_yyy_1[i] * fe_0 + ta_zz_yyyy_0[i] * pa_y[i] - ta_zz_yyyy_1[i] * pc_y[i];

        ta_yzz_yyyz_0[i] = 3.0 * ta_zz_yyz_0[i] * fe_0 - 3.0 * ta_zz_yyz_1[i] * fe_0 + ta_zz_yyyz_0[i] * pa_y[i] - ta_zz_yyyz_1[i] * pc_y[i];

        ta_yzz_yyzz_0[i] = 2.0 * ta_zz_yzz_0[i] * fe_0 - 2.0 * ta_zz_yzz_1[i] * fe_0 + ta_zz_yyzz_0[i] * pa_y[i] - ta_zz_yyzz_1[i] * pc_y[i];

        ta_yzz_yzzz_0[i] = ta_zz_zzz_0[i] * fe_0 - ta_zz_zzz_1[i] * fe_0 + ta_zz_yzzz_0[i] * pa_y[i] - ta_zz_yzzz_1[i] * pc_y[i];

        ta_yzz_zzzz_0[i] = ta_zz_zzzz_0[i] * pa_y[i] - ta_zz_zzzz_1[i] * pc_y[i];
    }

    // Set up 135-150 components of targeted buffer : FG

    auto ta_zzz_xxxx_0 = pbuffer.data(idx_npot_0_fg + 135);

    auto ta_zzz_xxxy_0 = pbuffer.data(idx_npot_0_fg + 136);

    auto ta_zzz_xxxz_0 = pbuffer.data(idx_npot_0_fg + 137);

    auto ta_zzz_xxyy_0 = pbuffer.data(idx_npot_0_fg + 138);

    auto ta_zzz_xxyz_0 = pbuffer.data(idx_npot_0_fg + 139);

    auto ta_zzz_xxzz_0 = pbuffer.data(idx_npot_0_fg + 140);

    auto ta_zzz_xyyy_0 = pbuffer.data(idx_npot_0_fg + 141);

    auto ta_zzz_xyyz_0 = pbuffer.data(idx_npot_0_fg + 142);

    auto ta_zzz_xyzz_0 = pbuffer.data(idx_npot_0_fg + 143);

    auto ta_zzz_xzzz_0 = pbuffer.data(idx_npot_0_fg + 144);

    auto ta_zzz_yyyy_0 = pbuffer.data(idx_npot_0_fg + 145);

    auto ta_zzz_yyyz_0 = pbuffer.data(idx_npot_0_fg + 146);

    auto ta_zzz_yyzz_0 = pbuffer.data(idx_npot_0_fg + 147);

    auto ta_zzz_yzzz_0 = pbuffer.data(idx_npot_0_fg + 148);

    auto ta_zzz_zzzz_0 = pbuffer.data(idx_npot_0_fg + 149);

#pragma omp simd aligned(pa_z,              \
                             pc_z,          \
                             ta_z_xxxx_0,   \
                             ta_z_xxxx_1,   \
                             ta_z_xxxy_0,   \
                             ta_z_xxxy_1,   \
                             ta_z_xxxz_0,   \
                             ta_z_xxxz_1,   \
                             ta_z_xxyy_0,   \
                             ta_z_xxyy_1,   \
                             ta_z_xxyz_0,   \
                             ta_z_xxyz_1,   \
                             ta_z_xxzz_0,   \
                             ta_z_xxzz_1,   \
                             ta_z_xyyy_0,   \
                             ta_z_xyyy_1,   \
                             ta_z_xyyz_0,   \
                             ta_z_xyyz_1,   \
                             ta_z_xyzz_0,   \
                             ta_z_xyzz_1,   \
                             ta_z_xzzz_0,   \
                             ta_z_xzzz_1,   \
                             ta_z_yyyy_0,   \
                             ta_z_yyyy_1,   \
                             ta_z_yyyz_0,   \
                             ta_z_yyyz_1,   \
                             ta_z_yyzz_0,   \
                             ta_z_yyzz_1,   \
                             ta_z_yzzz_0,   \
                             ta_z_yzzz_1,   \
                             ta_z_zzzz_0,   \
                             ta_z_zzzz_1,   \
                             ta_zz_xxx_0,   \
                             ta_zz_xxx_1,   \
                             ta_zz_xxxx_0,  \
                             ta_zz_xxxx_1,  \
                             ta_zz_xxxy_0,  \
                             ta_zz_xxxy_1,  \
                             ta_zz_xxxz_0,  \
                             ta_zz_xxxz_1,  \
                             ta_zz_xxy_0,   \
                             ta_zz_xxy_1,   \
                             ta_zz_xxyy_0,  \
                             ta_zz_xxyy_1,  \
                             ta_zz_xxyz_0,  \
                             ta_zz_xxyz_1,  \
                             ta_zz_xxz_0,   \
                             ta_zz_xxz_1,   \
                             ta_zz_xxzz_0,  \
                             ta_zz_xxzz_1,  \
                             ta_zz_xyy_0,   \
                             ta_zz_xyy_1,   \
                             ta_zz_xyyy_0,  \
                             ta_zz_xyyy_1,  \
                             ta_zz_xyyz_0,  \
                             ta_zz_xyyz_1,  \
                             ta_zz_xyz_0,   \
                             ta_zz_xyz_1,   \
                             ta_zz_xyzz_0,  \
                             ta_zz_xyzz_1,  \
                             ta_zz_xzz_0,   \
                             ta_zz_xzz_1,   \
                             ta_zz_xzzz_0,  \
                             ta_zz_xzzz_1,  \
                             ta_zz_yyy_0,   \
                             ta_zz_yyy_1,   \
                             ta_zz_yyyy_0,  \
                             ta_zz_yyyy_1,  \
                             ta_zz_yyyz_0,  \
                             ta_zz_yyyz_1,  \
                             ta_zz_yyz_0,   \
                             ta_zz_yyz_1,   \
                             ta_zz_yyzz_0,  \
                             ta_zz_yyzz_1,  \
                             ta_zz_yzz_0,   \
                             ta_zz_yzz_1,   \
                             ta_zz_yzzz_0,  \
                             ta_zz_yzzz_1,  \
                             ta_zz_zzz_0,   \
                             ta_zz_zzz_1,   \
                             ta_zz_zzzz_0,  \
                             ta_zz_zzzz_1,  \
                             ta_zzz_xxxx_0, \
                             ta_zzz_xxxy_0, \
                             ta_zzz_xxxz_0, \
                             ta_zzz_xxyy_0, \
                             ta_zzz_xxyz_0, \
                             ta_zzz_xxzz_0, \
                             ta_zzz_xyyy_0, \
                             ta_zzz_xyyz_0, \
                             ta_zzz_xyzz_0, \
                             ta_zzz_xzzz_0, \
                             ta_zzz_yyyy_0, \
                             ta_zzz_yyyz_0, \
                             ta_zzz_yyzz_0, \
                             ta_zzz_yzzz_0, \
                             ta_zzz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta_zzz_xxxx_0[i] = 2.0 * ta_z_xxxx_0[i] * fe_0 - 2.0 * ta_z_xxxx_1[i] * fe_0 + ta_zz_xxxx_0[i] * pa_z[i] - ta_zz_xxxx_1[i] * pc_z[i];

        ta_zzz_xxxy_0[i] = 2.0 * ta_z_xxxy_0[i] * fe_0 - 2.0 * ta_z_xxxy_1[i] * fe_0 + ta_zz_xxxy_0[i] * pa_z[i] - ta_zz_xxxy_1[i] * pc_z[i];

        ta_zzz_xxxz_0[i] = 2.0 * ta_z_xxxz_0[i] * fe_0 - 2.0 * ta_z_xxxz_1[i] * fe_0 + ta_zz_xxx_0[i] * fe_0 - ta_zz_xxx_1[i] * fe_0 +
                           ta_zz_xxxz_0[i] * pa_z[i] - ta_zz_xxxz_1[i] * pc_z[i];

        ta_zzz_xxyy_0[i] = 2.0 * ta_z_xxyy_0[i] * fe_0 - 2.0 * ta_z_xxyy_1[i] * fe_0 + ta_zz_xxyy_0[i] * pa_z[i] - ta_zz_xxyy_1[i] * pc_z[i];

        ta_zzz_xxyz_0[i] = 2.0 * ta_z_xxyz_0[i] * fe_0 - 2.0 * ta_z_xxyz_1[i] * fe_0 + ta_zz_xxy_0[i] * fe_0 - ta_zz_xxy_1[i] * fe_0 +
                           ta_zz_xxyz_0[i] * pa_z[i] - ta_zz_xxyz_1[i] * pc_z[i];

        ta_zzz_xxzz_0[i] = 2.0 * ta_z_xxzz_0[i] * fe_0 - 2.0 * ta_z_xxzz_1[i] * fe_0 + 2.0 * ta_zz_xxz_0[i] * fe_0 - 2.0 * ta_zz_xxz_1[i] * fe_0 +
                           ta_zz_xxzz_0[i] * pa_z[i] - ta_zz_xxzz_1[i] * pc_z[i];

        ta_zzz_xyyy_0[i] = 2.0 * ta_z_xyyy_0[i] * fe_0 - 2.0 * ta_z_xyyy_1[i] * fe_0 + ta_zz_xyyy_0[i] * pa_z[i] - ta_zz_xyyy_1[i] * pc_z[i];

        ta_zzz_xyyz_0[i] = 2.0 * ta_z_xyyz_0[i] * fe_0 - 2.0 * ta_z_xyyz_1[i] * fe_0 + ta_zz_xyy_0[i] * fe_0 - ta_zz_xyy_1[i] * fe_0 +
                           ta_zz_xyyz_0[i] * pa_z[i] - ta_zz_xyyz_1[i] * pc_z[i];

        ta_zzz_xyzz_0[i] = 2.0 * ta_z_xyzz_0[i] * fe_0 - 2.0 * ta_z_xyzz_1[i] * fe_0 + 2.0 * ta_zz_xyz_0[i] * fe_0 - 2.0 * ta_zz_xyz_1[i] * fe_0 +
                           ta_zz_xyzz_0[i] * pa_z[i] - ta_zz_xyzz_1[i] * pc_z[i];

        ta_zzz_xzzz_0[i] = 2.0 * ta_z_xzzz_0[i] * fe_0 - 2.0 * ta_z_xzzz_1[i] * fe_0 + 3.0 * ta_zz_xzz_0[i] * fe_0 - 3.0 * ta_zz_xzz_1[i] * fe_0 +
                           ta_zz_xzzz_0[i] * pa_z[i] - ta_zz_xzzz_1[i] * pc_z[i];

        ta_zzz_yyyy_0[i] = 2.0 * ta_z_yyyy_0[i] * fe_0 - 2.0 * ta_z_yyyy_1[i] * fe_0 + ta_zz_yyyy_0[i] * pa_z[i] - ta_zz_yyyy_1[i] * pc_z[i];

        ta_zzz_yyyz_0[i] = 2.0 * ta_z_yyyz_0[i] * fe_0 - 2.0 * ta_z_yyyz_1[i] * fe_0 + ta_zz_yyy_0[i] * fe_0 - ta_zz_yyy_1[i] * fe_0 +
                           ta_zz_yyyz_0[i] * pa_z[i] - ta_zz_yyyz_1[i] * pc_z[i];

        ta_zzz_yyzz_0[i] = 2.0 * ta_z_yyzz_0[i] * fe_0 - 2.0 * ta_z_yyzz_1[i] * fe_0 + 2.0 * ta_zz_yyz_0[i] * fe_0 - 2.0 * ta_zz_yyz_1[i] * fe_0 +
                           ta_zz_yyzz_0[i] * pa_z[i] - ta_zz_yyzz_1[i] * pc_z[i];

        ta_zzz_yzzz_0[i] = 2.0 * ta_z_yzzz_0[i] * fe_0 - 2.0 * ta_z_yzzz_1[i] * fe_0 + 3.0 * ta_zz_yzz_0[i] * fe_0 - 3.0 * ta_zz_yzz_1[i] * fe_0 +
                           ta_zz_yzzz_0[i] * pa_z[i] - ta_zz_yzzz_1[i] * pc_z[i];

        ta_zzz_zzzz_0[i] = 2.0 * ta_z_zzzz_0[i] * fe_0 - 2.0 * ta_z_zzzz_1[i] * fe_0 + 4.0 * ta_zz_zzz_0[i] * fe_0 - 4.0 * ta_zz_zzz_1[i] * fe_0 +
                           ta_zz_zzzz_0[i] * pa_z[i] - ta_zz_zzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
