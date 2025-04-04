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

#include "NuclearPotentialGeom010PrimRecID.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_id(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_id,
                                        const size_t              idx_npot_geom_010_0_gd,
                                        const size_t              idx_npot_geom_010_1_gd,
                                        const size_t              idx_npot_geom_010_0_hp,
                                        const size_t              idx_npot_geom_010_1_hp,
                                        const size_t              idx_npot_1_hd,
                                        const size_t              idx_npot_geom_010_0_hd,
                                        const size_t              idx_npot_geom_010_1_hd,
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

    auto ta1_x_xxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd);

    auto ta1_x_xxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 1);

    auto ta1_x_xxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 2);

    auto ta1_x_xxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 3);

    auto ta1_x_xxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 4);

    auto ta1_x_xxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 5);

    auto ta1_x_xxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 6);

    auto ta1_x_xxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 7);

    auto ta1_x_xxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 8);

    auto ta1_x_xxxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 11);

    auto ta1_x_xxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 12);

    auto ta1_x_xxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 13);

    auto ta1_x_xxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 14);

    auto ta1_x_xxxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 15);

    auto ta1_x_xxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 17);

    auto ta1_x_xxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 18);

    auto ta1_x_xxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 19);

    auto ta1_x_xxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 20);

    auto ta1_x_xxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 21);

    auto ta1_x_xxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 22);

    auto ta1_x_xxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 23);

    auto ta1_x_xxyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 26);

    auto ta1_x_xxyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 29);

    auto ta1_x_xxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 30);

    auto ta1_x_xxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 31);

    auto ta1_x_xxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 32);

    auto ta1_x_xxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 33);

    auto ta1_x_xxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 34);

    auto ta1_x_xxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 35);

    auto ta1_x_xyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 36);

    auto ta1_x_xyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 37);

    auto ta1_x_xyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 38);

    auto ta1_x_xyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 39);

    auto ta1_x_xyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 40);

    auto ta1_x_xyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 43);

    auto ta1_x_xyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 44);

    auto ta1_x_xyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 48);

    auto ta1_x_xyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 50);

    auto ta1_x_xzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 54);

    auto ta1_x_xzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 55);

    auto ta1_x_xzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 56);

    auto ta1_x_xzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 58);

    auto ta1_x_xzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 59);

    auto ta1_x_yyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 60);

    auto ta1_x_yyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 61);

    auto ta1_x_yyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 62);

    auto ta1_x_yyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 63);

    auto ta1_x_yyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 64);

    auto ta1_x_yyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 65);

    auto ta1_x_yyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 67);

    auto ta1_x_yyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 68);

    auto ta1_x_yyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 69);

    auto ta1_x_yyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 71);

    auto ta1_x_yyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 72);

    auto ta1_x_yyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 73);

    auto ta1_x_yyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 74);

    auto ta1_x_yyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 75);

    auto ta1_x_yyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 76);

    auto ta1_x_yyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 77);

    auto ta1_x_yzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 78);

    auto ta1_x_yzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 80);

    auto ta1_x_yzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 82);

    auto ta1_x_yzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 83);

    auto ta1_x_zzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 84);

    auto ta1_x_zzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 85);

    auto ta1_x_zzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 86);

    auto ta1_x_zzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 87);

    auto ta1_x_zzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 88);

    auto ta1_x_zzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 89);

    auto ta1_y_xxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 90);

    auto ta1_y_xxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 91);

    auto ta1_y_xxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 92);

    auto ta1_y_xxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 93);

    auto ta1_y_xxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 94);

    auto ta1_y_xxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 95);

    auto ta1_y_xxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 96);

    auto ta1_y_xxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 97);

    auto ta1_y_xxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 98);

    auto ta1_y_xxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 99);

    auto ta1_y_xxxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 100);

    auto ta1_y_xxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 102);

    auto ta1_y_xxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 103);

    auto ta1_y_xxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 106);

    auto ta1_y_xxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 107);

    auto ta1_y_xxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 108);

    auto ta1_y_xxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 109);

    auto ta1_y_xxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 110);

    auto ta1_y_xxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 111);

    auto ta1_y_xxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 112);

    auto ta1_y_xxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 113);

    auto ta1_y_xxyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 115);

    auto ta1_y_xxyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 118);

    auto ta1_y_xxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 120);

    auto ta1_y_xxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 121);

    auto ta1_y_xxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 122);

    auto ta1_y_xxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 123);

    auto ta1_y_xxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 124);

    auto ta1_y_xxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 125);

    auto ta1_y_xyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 127);

    auto ta1_y_xyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 129);

    auto ta1_y_xyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 130);

    auto ta1_y_xyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 131);

    auto ta1_y_xyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 136);

    auto ta1_y_xyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 137);

    auto ta1_y_xyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 141);

    auto ta1_y_xyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 142);

    auto ta1_y_xzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 146);

    auto ta1_y_xzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 147);

    auto ta1_y_xzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 148);

    auto ta1_y_xzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 149);

    auto ta1_y_yyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 150);

    auto ta1_y_yyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 151);

    auto ta1_y_yyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 152);

    auto ta1_y_yyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 153);

    auto ta1_y_yyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 154);

    auto ta1_y_yyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 155);

    auto ta1_y_yyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 156);

    auto ta1_y_yyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 157);

    auto ta1_y_yyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 159);

    auto ta1_y_yyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 160);

    auto ta1_y_yyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 161);

    auto ta1_y_yyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 162);

    auto ta1_y_yyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 163);

    auto ta1_y_yyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 164);

    auto ta1_y_yyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 165);

    auto ta1_y_yyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 166);

    auto ta1_y_yyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 167);

    auto ta1_y_yzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 169);

    auto ta1_y_yzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 170);

    auto ta1_y_yzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 171);

    auto ta1_y_yzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 172);

    auto ta1_y_yzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 173);

    auto ta1_y_zzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 174);

    auto ta1_y_zzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 175);

    auto ta1_y_zzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 176);

    auto ta1_y_zzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 177);

    auto ta1_y_zzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 178);

    auto ta1_y_zzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 179);

    auto ta1_z_xxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 180);

    auto ta1_z_xxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 181);

    auto ta1_z_xxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 182);

    auto ta1_z_xxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 183);

    auto ta1_z_xxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 184);

    auto ta1_z_xxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 185);

    auto ta1_z_xxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 186);

    auto ta1_z_xxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 188);

    auto ta1_z_xxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 189);

    auto ta1_z_xxxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 190);

    auto ta1_z_xxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 192);

    auto ta1_z_xxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 193);

    auto ta1_z_xxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 194);

    auto ta1_z_xxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 196);

    auto ta1_z_xxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 197);

    auto ta1_z_xxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 198);

    auto ta1_z_xxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 199);

    auto ta1_z_xxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 200);

    auto ta1_z_xxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 201);

    auto ta1_z_xxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 202);

    auto ta1_z_xxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 203);

    auto ta1_z_xxyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 206);

    auto ta1_z_xxyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 208);

    auto ta1_z_xxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 210);

    auto ta1_z_xxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 211);

    auto ta1_z_xxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 212);

    auto ta1_z_xxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 213);

    auto ta1_z_xxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 214);

    auto ta1_z_xxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 215);

    auto ta1_z_xyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 217);

    auto ta1_z_xyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 219);

    auto ta1_z_xyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 220);

    auto ta1_z_xyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 221);

    auto ta1_z_xyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 226);

    auto ta1_z_xyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 227);

    auto ta1_z_xyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 231);

    auto ta1_z_xyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 232);

    auto ta1_z_xzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 236);

    auto ta1_z_xzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 237);

    auto ta1_z_xzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 238);

    auto ta1_z_xzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 239);

    auto ta1_z_yyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 240);

    auto ta1_z_yyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 241);

    auto ta1_z_yyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 242);

    auto ta1_z_yyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 243);

    auto ta1_z_yyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 244);

    auto ta1_z_yyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 245);

    auto ta1_z_yyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 247);

    auto ta1_z_yyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 248);

    auto ta1_z_yyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 249);

    auto ta1_z_yyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 250);

    auto ta1_z_yyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 251);

    auto ta1_z_yyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 252);

    auto ta1_z_yyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 253);

    auto ta1_z_yyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 254);

    auto ta1_z_yyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 255);

    auto ta1_z_yyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 256);

    auto ta1_z_yyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 257);

    auto ta1_z_yzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 258);

    auto ta1_z_yzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 260);

    auto ta1_z_yzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 261);

    auto ta1_z_yzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 262);

    auto ta1_z_yzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 263);

    auto ta1_z_zzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_gd + 264);

    auto ta1_z_zzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 265);

    auto ta1_z_zzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 266);

    auto ta1_z_zzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_gd + 267);

    auto ta1_z_zzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 268);

    auto ta1_z_zzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_gd + 269);

    // Set up components of auxiliary buffer : GD

    auto ta1_x_xxxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd);

    auto ta1_x_xxxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 1);

    auto ta1_x_xxxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 2);

    auto ta1_x_xxxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 3);

    auto ta1_x_xxxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 4);

    auto ta1_x_xxxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 5);

    auto ta1_x_xxxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 6);

    auto ta1_x_xxxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 7);

    auto ta1_x_xxxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 8);

    auto ta1_x_xxxy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 11);

    auto ta1_x_xxxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 12);

    auto ta1_x_xxxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 13);

    auto ta1_x_xxxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 14);

    auto ta1_x_xxxz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 15);

    auto ta1_x_xxxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 17);

    auto ta1_x_xxyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 18);

    auto ta1_x_xxyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 19);

    auto ta1_x_xxyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 20);

    auto ta1_x_xxyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 21);

    auto ta1_x_xxyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 22);

    auto ta1_x_xxyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 23);

    auto ta1_x_xxyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 26);

    auto ta1_x_xxyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 29);

    auto ta1_x_xxzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 30);

    auto ta1_x_xxzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 31);

    auto ta1_x_xxzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 32);

    auto ta1_x_xxzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 33);

    auto ta1_x_xxzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 34);

    auto ta1_x_xxzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 35);

    auto ta1_x_xyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 36);

    auto ta1_x_xyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 37);

    auto ta1_x_xyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 38);

    auto ta1_x_xyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 39);

    auto ta1_x_xyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 40);

    auto ta1_x_xyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 43);

    auto ta1_x_xyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 44);

    auto ta1_x_xyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 48);

    auto ta1_x_xyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 50);

    auto ta1_x_xzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 54);

    auto ta1_x_xzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 55);

    auto ta1_x_xzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 56);

    auto ta1_x_xzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 58);

    auto ta1_x_xzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 59);

    auto ta1_x_yyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 60);

    auto ta1_x_yyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 61);

    auto ta1_x_yyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 62);

    auto ta1_x_yyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 63);

    auto ta1_x_yyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 64);

    auto ta1_x_yyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 65);

    auto ta1_x_yyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 67);

    auto ta1_x_yyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 68);

    auto ta1_x_yyyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 69);

    auto ta1_x_yyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 71);

    auto ta1_x_yyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 72);

    auto ta1_x_yyzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 73);

    auto ta1_x_yyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 74);

    auto ta1_x_yyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 75);

    auto ta1_x_yyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 76);

    auto ta1_x_yyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 77);

    auto ta1_x_yzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 78);

    auto ta1_x_yzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 80);

    auto ta1_x_yzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 82);

    auto ta1_x_yzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 83);

    auto ta1_x_zzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 84);

    auto ta1_x_zzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 85);

    auto ta1_x_zzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 86);

    auto ta1_x_zzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 87);

    auto ta1_x_zzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 88);

    auto ta1_x_zzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 89);

    auto ta1_y_xxxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 90);

    auto ta1_y_xxxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 91);

    auto ta1_y_xxxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 92);

    auto ta1_y_xxxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 93);

    auto ta1_y_xxxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 94);

    auto ta1_y_xxxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 95);

    auto ta1_y_xxxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 96);

    auto ta1_y_xxxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 97);

    auto ta1_y_xxxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 98);

    auto ta1_y_xxxy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 99);

    auto ta1_y_xxxy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 100);

    auto ta1_y_xxxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 102);

    auto ta1_y_xxxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 103);

    auto ta1_y_xxxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 106);

    auto ta1_y_xxxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 107);

    auto ta1_y_xxyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 108);

    auto ta1_y_xxyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 109);

    auto ta1_y_xxyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 110);

    auto ta1_y_xxyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 111);

    auto ta1_y_xxyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 112);

    auto ta1_y_xxyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 113);

    auto ta1_y_xxyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 115);

    auto ta1_y_xxyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 118);

    auto ta1_y_xxzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 120);

    auto ta1_y_xxzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 121);

    auto ta1_y_xxzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 122);

    auto ta1_y_xxzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 123);

    auto ta1_y_xxzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 124);

    auto ta1_y_xxzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 125);

    auto ta1_y_xyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 127);

    auto ta1_y_xyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 129);

    auto ta1_y_xyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 130);

    auto ta1_y_xyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 131);

    auto ta1_y_xyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 136);

    auto ta1_y_xyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 137);

    auto ta1_y_xyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 141);

    auto ta1_y_xyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 142);

    auto ta1_y_xzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 146);

    auto ta1_y_xzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 147);

    auto ta1_y_xzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 148);

    auto ta1_y_xzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 149);

    auto ta1_y_yyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 150);

    auto ta1_y_yyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 151);

    auto ta1_y_yyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 152);

    auto ta1_y_yyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 153);

    auto ta1_y_yyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 154);

    auto ta1_y_yyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 155);

    auto ta1_y_yyyz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 156);

    auto ta1_y_yyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 157);

    auto ta1_y_yyyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 159);

    auto ta1_y_yyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 160);

    auto ta1_y_yyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 161);

    auto ta1_y_yyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 162);

    auto ta1_y_yyzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 163);

    auto ta1_y_yyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 164);

    auto ta1_y_yyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 165);

    auto ta1_y_yyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 166);

    auto ta1_y_yyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 167);

    auto ta1_y_yzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 169);

    auto ta1_y_yzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 170);

    auto ta1_y_yzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 171);

    auto ta1_y_yzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 172);

    auto ta1_y_yzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 173);

    auto ta1_y_zzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 174);

    auto ta1_y_zzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 175);

    auto ta1_y_zzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 176);

    auto ta1_y_zzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 177);

    auto ta1_y_zzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 178);

    auto ta1_y_zzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 179);

    auto ta1_z_xxxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 180);

    auto ta1_z_xxxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 181);

    auto ta1_z_xxxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 182);

    auto ta1_z_xxxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 183);

    auto ta1_z_xxxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 184);

    auto ta1_z_xxxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 185);

    auto ta1_z_xxxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 186);

    auto ta1_z_xxxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 188);

    auto ta1_z_xxxy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 189);

    auto ta1_z_xxxy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 190);

    auto ta1_z_xxxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 192);

    auto ta1_z_xxxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 193);

    auto ta1_z_xxxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 194);

    auto ta1_z_xxxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 196);

    auto ta1_z_xxxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 197);

    auto ta1_z_xxyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 198);

    auto ta1_z_xxyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 199);

    auto ta1_z_xxyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 200);

    auto ta1_z_xxyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 201);

    auto ta1_z_xxyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 202);

    auto ta1_z_xxyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 203);

    auto ta1_z_xxyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 206);

    auto ta1_z_xxyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 208);

    auto ta1_z_xxzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 210);

    auto ta1_z_xxzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 211);

    auto ta1_z_xxzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 212);

    auto ta1_z_xxzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 213);

    auto ta1_z_xxzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 214);

    auto ta1_z_xxzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 215);

    auto ta1_z_xyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 217);

    auto ta1_z_xyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 219);

    auto ta1_z_xyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 220);

    auto ta1_z_xyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 221);

    auto ta1_z_xyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 226);

    auto ta1_z_xyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 227);

    auto ta1_z_xyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 231);

    auto ta1_z_xyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 232);

    auto ta1_z_xzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 236);

    auto ta1_z_xzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 237);

    auto ta1_z_xzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 238);

    auto ta1_z_xzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 239);

    auto ta1_z_yyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 240);

    auto ta1_z_yyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 241);

    auto ta1_z_yyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 242);

    auto ta1_z_yyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 243);

    auto ta1_z_yyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 244);

    auto ta1_z_yyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 245);

    auto ta1_z_yyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 247);

    auto ta1_z_yyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 248);

    auto ta1_z_yyyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 249);

    auto ta1_z_yyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 250);

    auto ta1_z_yyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 251);

    auto ta1_z_yyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 252);

    auto ta1_z_yyzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 253);

    auto ta1_z_yyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 254);

    auto ta1_z_yyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 255);

    auto ta1_z_yyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 256);

    auto ta1_z_yyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 257);

    auto ta1_z_yzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 258);

    auto ta1_z_yzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 260);

    auto ta1_z_yzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 261);

    auto ta1_z_yzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 262);

    auto ta1_z_yzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 263);

    auto ta1_z_zzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_gd + 264);

    auto ta1_z_zzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 265);

    auto ta1_z_zzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 266);

    auto ta1_z_zzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_gd + 267);

    auto ta1_z_zzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 268);

    auto ta1_z_zzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_gd + 269);

    // Set up components of auxiliary buffer : HP

    auto ta1_x_xxxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_hp);

    auto ta1_x_xxxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 1);

    auto ta1_x_xxxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 2);

    auto ta1_x_xxxxy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 3);

    auto ta1_x_xxxxz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 6);

    auto ta1_x_xxxxz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 8);

    auto ta1_x_xxxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 9);

    auto ta1_x_xxxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 10);

    auto ta1_x_xxxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 15);

    auto ta1_x_xxxzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 16);

    auto ta1_x_xxxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 17);

    auto ta1_x_xxyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 18);

    auto ta1_x_xxyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 19);

    auto ta1_x_xxzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 27);

    auto ta1_x_xxzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 28);

    auto ta1_x_xxzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 29);

    auto ta1_x_xzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 42);

    auto ta1_x_yyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 45);

    auto ta1_x_yyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 46);

    auto ta1_x_yyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 47);

    auto ta1_x_yyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 53);

    auto ta1_x_yyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 56);

    auto ta1_x_yzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 59);

    auto ta1_x_zzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 60);

    auto ta1_x_zzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 61);

    auto ta1_x_zzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 62);

    auto ta1_y_xxxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 63);

    auto ta1_y_xxxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 64);

    auto ta1_y_xxxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 65);

    auto ta1_y_xxxyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 72);

    auto ta1_y_xxxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 73);

    auto ta1_y_xxxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 80);

    auto ta1_y_xxyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 81);

    auto ta1_y_xxyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 82);

    auto ta1_y_xxzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 92);

    auto ta1_y_xyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 94);

    auto ta1_y_xzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 107);

    auto ta1_y_yyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 108);

    auto ta1_y_yyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 109);

    auto ta1_y_yyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 110);

    auto ta1_y_yyyyz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 112);

    auto ta1_y_yyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 113);

    auto ta1_y_yyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 114);

    auto ta1_y_yyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 115);

    auto ta1_y_yyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 116);

    auto ta1_y_yyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 117);

    auto ta1_y_yyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 118);

    auto ta1_y_yyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 119);

    auto ta1_y_yzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 121);

    auto ta1_y_zzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 123);

    auto ta1_y_zzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 124);

    auto ta1_y_zzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 125);

    auto ta1_z_xxxxx_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 126);

    auto ta1_z_xxxxx_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 127);

    auto ta1_z_xxxxx_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 128);

    auto ta1_z_xxxyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 136);

    auto ta1_z_xxxzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 141);

    auto ta1_z_xxxzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 143);

    auto ta1_z_xxyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 145);

    auto ta1_z_xxzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 153);

    auto ta1_z_xxzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 155);

    auto ta1_z_xyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 157);

    auto ta1_z_xzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 170);

    auto ta1_z_yyyyy_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 171);

    auto ta1_z_yyyyy_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 172);

    auto ta1_z_yyyyy_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 173);

    auto ta1_z_yyyyz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 176);

    auto ta1_z_yyyzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 177);

    auto ta1_z_yyyzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 178);

    auto ta1_z_yyyzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 179);

    auto ta1_z_yyzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 180);

    auto ta1_z_yyzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 181);

    auto ta1_z_yyzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 182);

    auto ta1_z_yzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 184);

    auto ta1_z_yzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 185);

    auto ta1_z_zzzzz_x_0 = pbuffer.data(idx_npot_geom_010_0_hp + 186);

    auto ta1_z_zzzzz_y_0 = pbuffer.data(idx_npot_geom_010_0_hp + 187);

    auto ta1_z_zzzzz_z_0 = pbuffer.data(idx_npot_geom_010_0_hp + 188);

    // Set up components of auxiliary buffer : HP

    auto ta1_x_xxxxx_x_1 = pbuffer.data(idx_npot_geom_010_1_hp);

    auto ta1_x_xxxxx_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 1);

    auto ta1_x_xxxxx_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 2);

    auto ta1_x_xxxxy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 3);

    auto ta1_x_xxxxz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 6);

    auto ta1_x_xxxxz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 8);

    auto ta1_x_xxxyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 9);

    auto ta1_x_xxxyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 10);

    auto ta1_x_xxxzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 15);

    auto ta1_x_xxxzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 16);

    auto ta1_x_xxxzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 17);

    auto ta1_x_xxyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 18);

    auto ta1_x_xxyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 19);

    auto ta1_x_xxzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 27);

    auto ta1_x_xxzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 28);

    auto ta1_x_xxzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 29);

    auto ta1_x_xzzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 42);

    auto ta1_x_yyyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 45);

    auto ta1_x_yyyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 46);

    auto ta1_x_yyyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 47);

    auto ta1_x_yyyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 53);

    auto ta1_x_yyzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 56);

    auto ta1_x_yzzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 59);

    auto ta1_x_zzzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 60);

    auto ta1_x_zzzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 61);

    auto ta1_x_zzzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 62);

    auto ta1_y_xxxxx_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 63);

    auto ta1_y_xxxxx_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 64);

    auto ta1_y_xxxxx_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 65);

    auto ta1_y_xxxyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 72);

    auto ta1_y_xxxyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 73);

    auto ta1_y_xxxzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 80);

    auto ta1_y_xxyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 81);

    auto ta1_y_xxyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 82);

    auto ta1_y_xxzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 92);

    auto ta1_y_xyyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 94);

    auto ta1_y_xzzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 107);

    auto ta1_y_yyyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 108);

    auto ta1_y_yyyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 109);

    auto ta1_y_yyyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 110);

    auto ta1_y_yyyyz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 112);

    auto ta1_y_yyyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 113);

    auto ta1_y_yyyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 114);

    auto ta1_y_yyyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 115);

    auto ta1_y_yyyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 116);

    auto ta1_y_yyzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 117);

    auto ta1_y_yyzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 118);

    auto ta1_y_yyzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 119);

    auto ta1_y_yzzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 121);

    auto ta1_y_zzzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 123);

    auto ta1_y_zzzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 124);

    auto ta1_y_zzzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 125);

    auto ta1_z_xxxxx_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 126);

    auto ta1_z_xxxxx_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 127);

    auto ta1_z_xxxxx_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 128);

    auto ta1_z_xxxyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 136);

    auto ta1_z_xxxzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 141);

    auto ta1_z_xxxzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 143);

    auto ta1_z_xxyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 145);

    auto ta1_z_xxzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 153);

    auto ta1_z_xxzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 155);

    auto ta1_z_xyyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 157);

    auto ta1_z_xzzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 170);

    auto ta1_z_yyyyy_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 171);

    auto ta1_z_yyyyy_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 172);

    auto ta1_z_yyyyy_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 173);

    auto ta1_z_yyyyz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 176);

    auto ta1_z_yyyzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 177);

    auto ta1_z_yyyzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 178);

    auto ta1_z_yyyzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 179);

    auto ta1_z_yyzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 180);

    auto ta1_z_yyzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 181);

    auto ta1_z_yyzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 182);

    auto ta1_z_yzzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 184);

    auto ta1_z_yzzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 185);

    auto ta1_z_zzzzz_x_1 = pbuffer.data(idx_npot_geom_010_1_hp + 186);

    auto ta1_z_zzzzz_y_1 = pbuffer.data(idx_npot_geom_010_1_hp + 187);

    auto ta1_z_zzzzz_z_1 = pbuffer.data(idx_npot_geom_010_1_hp + 188);

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

    auto ta_xxxxz_xx_1 = pbuffer.data(idx_npot_1_hd + 12);

    auto ta_xxxxz_xy_1 = pbuffer.data(idx_npot_1_hd + 13);

    auto ta_xxxxz_xz_1 = pbuffer.data(idx_npot_1_hd + 14);

    auto ta_xxxxz_zz_1 = pbuffer.data(idx_npot_1_hd + 17);

    auto ta_xxxyy_xx_1 = pbuffer.data(idx_npot_1_hd + 18);

    auto ta_xxxyy_xy_1 = pbuffer.data(idx_npot_1_hd + 19);

    auto ta_xxxyy_xz_1 = pbuffer.data(idx_npot_1_hd + 20);

    auto ta_xxxyy_yy_1 = pbuffer.data(idx_npot_1_hd + 21);

    auto ta_xxxyy_yz_1 = pbuffer.data(idx_npot_1_hd + 22);

    auto ta_xxxzz_xx_1 = pbuffer.data(idx_npot_1_hd + 30);

    auto ta_xxxzz_xy_1 = pbuffer.data(idx_npot_1_hd + 31);

    auto ta_xxxzz_xz_1 = pbuffer.data(idx_npot_1_hd + 32);

    auto ta_xxxzz_yz_1 = pbuffer.data(idx_npot_1_hd + 34);

    auto ta_xxxzz_zz_1 = pbuffer.data(idx_npot_1_hd + 35);

    auto ta_xxyyy_xx_1 = pbuffer.data(idx_npot_1_hd + 36);

    auto ta_xxyyy_xy_1 = pbuffer.data(idx_npot_1_hd + 37);

    auto ta_xxyyy_xz_1 = pbuffer.data(idx_npot_1_hd + 38);

    auto ta_xxyyy_yy_1 = pbuffer.data(idx_npot_1_hd + 39);

    auto ta_xxyyy_yz_1 = pbuffer.data(idx_npot_1_hd + 40);

    auto ta_xxyyz_xy_1 = pbuffer.data(idx_npot_1_hd + 43);

    auto ta_xxyzz_xz_1 = pbuffer.data(idx_npot_1_hd + 50);

    auto ta_xxzzz_xx_1 = pbuffer.data(idx_npot_1_hd + 54);

    auto ta_xxzzz_xy_1 = pbuffer.data(idx_npot_1_hd + 55);

    auto ta_xxzzz_xz_1 = pbuffer.data(idx_npot_1_hd + 56);

    auto ta_xxzzz_yz_1 = pbuffer.data(idx_npot_1_hd + 58);

    auto ta_xxzzz_zz_1 = pbuffer.data(idx_npot_1_hd + 59);

    auto ta_xyyyy_xx_1 = pbuffer.data(idx_npot_1_hd + 60);

    auto ta_xyyyy_xy_1 = pbuffer.data(idx_npot_1_hd + 61);

    auto ta_xyyyy_yy_1 = pbuffer.data(idx_npot_1_hd + 63);

    auto ta_xyyyy_yz_1 = pbuffer.data(idx_npot_1_hd + 64);

    auto ta_xyyzz_yz_1 = pbuffer.data(idx_npot_1_hd + 76);

    auto ta_xzzzz_xx_1 = pbuffer.data(idx_npot_1_hd + 84);

    auto ta_xzzzz_xz_1 = pbuffer.data(idx_npot_1_hd + 86);

    auto ta_xzzzz_yz_1 = pbuffer.data(idx_npot_1_hd + 88);

    auto ta_xzzzz_zz_1 = pbuffer.data(idx_npot_1_hd + 89);

    auto ta_yyyyy_xx_1 = pbuffer.data(idx_npot_1_hd + 90);

    auto ta_yyyyy_xy_1 = pbuffer.data(idx_npot_1_hd + 91);

    auto ta_yyyyy_xz_1 = pbuffer.data(idx_npot_1_hd + 92);

    auto ta_yyyyy_yy_1 = pbuffer.data(idx_npot_1_hd + 93);

    auto ta_yyyyy_yz_1 = pbuffer.data(idx_npot_1_hd + 94);

    auto ta_yyyyy_zz_1 = pbuffer.data(idx_npot_1_hd + 95);

    auto ta_yyyyz_xy_1 = pbuffer.data(idx_npot_1_hd + 97);

    auto ta_yyyyz_yy_1 = pbuffer.data(idx_npot_1_hd + 99);

    auto ta_yyyyz_yz_1 = pbuffer.data(idx_npot_1_hd + 100);

    auto ta_yyyyz_zz_1 = pbuffer.data(idx_npot_1_hd + 101);

    auto ta_yyyzz_xy_1 = pbuffer.data(idx_npot_1_hd + 103);

    auto ta_yyyzz_xz_1 = pbuffer.data(idx_npot_1_hd + 104);

    auto ta_yyyzz_yy_1 = pbuffer.data(idx_npot_1_hd + 105);

    auto ta_yyyzz_yz_1 = pbuffer.data(idx_npot_1_hd + 106);

    auto ta_yyyzz_zz_1 = pbuffer.data(idx_npot_1_hd + 107);

    auto ta_yyzzz_xy_1 = pbuffer.data(idx_npot_1_hd + 109);

    auto ta_yyzzz_xz_1 = pbuffer.data(idx_npot_1_hd + 110);

    auto ta_yyzzz_yy_1 = pbuffer.data(idx_npot_1_hd + 111);

    auto ta_yyzzz_yz_1 = pbuffer.data(idx_npot_1_hd + 112);

    auto ta_yyzzz_zz_1 = pbuffer.data(idx_npot_1_hd + 113);

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

    // Set up components of auxiliary buffer : HD

    auto ta1_x_xxxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd);

    auto ta1_x_xxxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 1);

    auto ta1_x_xxxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 2);

    auto ta1_x_xxxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 3);

    auto ta1_x_xxxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 4);

    auto ta1_x_xxxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 5);

    auto ta1_x_xxxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 6);

    auto ta1_x_xxxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 7);

    auto ta1_x_xxxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 8);

    auto ta1_x_xxxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 9);

    auto ta1_x_xxxxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 11);

    auto ta1_x_xxxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 12);

    auto ta1_x_xxxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 13);

    auto ta1_x_xxxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 14);

    auto ta1_x_xxxxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 15);

    auto ta1_x_xxxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 16);

    auto ta1_x_xxxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 17);

    auto ta1_x_xxxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 18);

    auto ta1_x_xxxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 19);

    auto ta1_x_xxxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 20);

    auto ta1_x_xxxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 21);

    auto ta1_x_xxxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 22);

    auto ta1_x_xxxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 23);

    auto ta1_x_xxxyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 26);

    auto ta1_x_xxxyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 29);

    auto ta1_x_xxxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 30);

    auto ta1_x_xxxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 31);

    auto ta1_x_xxxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 32);

    auto ta1_x_xxxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 33);

    auto ta1_x_xxxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 34);

    auto ta1_x_xxxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 35);

    auto ta1_x_xxyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 36);

    auto ta1_x_xxyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 37);

    auto ta1_x_xxyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 38);

    auto ta1_x_xxyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 39);

    auto ta1_x_xxyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 40);

    auto ta1_x_xxyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 41);

    auto ta1_x_xxyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 43);

    auto ta1_x_xxyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 44);

    auto ta1_x_xxyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 45);

    auto ta1_x_xxyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 47);

    auto ta1_x_xxyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 48);

    auto ta1_x_xxyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 50);

    auto ta1_x_xxyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 53);

    auto ta1_x_xxzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 54);

    auto ta1_x_xxzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 55);

    auto ta1_x_xxzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 56);

    auto ta1_x_xxzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 57);

    auto ta1_x_xxzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 58);

    auto ta1_x_xxzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 59);

    auto ta1_x_xyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 60);

    auto ta1_x_xyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 61);

    auto ta1_x_xyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 62);

    auto ta1_x_xyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 63);

    auto ta1_x_xyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 64);

    auto ta1_x_xyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 67);

    auto ta1_x_xyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 68);

    auto ta1_x_xyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 72);

    auto ta1_x_xyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 73);

    auto ta1_x_xyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 74);

    auto ta1_x_xyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 76);

    auto ta1_x_xyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 78);

    auto ta1_x_xyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 80);

    auto ta1_x_xzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 84);

    auto ta1_x_xzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 85);

    auto ta1_x_xzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 86);

    auto ta1_x_xzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 88);

    auto ta1_x_xzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 89);

    auto ta1_x_yyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 90);

    auto ta1_x_yyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 91);

    auto ta1_x_yyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 92);

    auto ta1_x_yyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 93);

    auto ta1_x_yyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 94);

    auto ta1_x_yyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 95);

    auto ta1_x_yyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 97);

    auto ta1_x_yyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 98);

    auto ta1_x_yyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 99);

    auto ta1_x_yyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 100);

    auto ta1_x_yyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 101);

    auto ta1_x_yyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 102);

    auto ta1_x_yyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 103);

    auto ta1_x_yyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 104);

    auto ta1_x_yyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 105);

    auto ta1_x_yyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 106);

    auto ta1_x_yyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 107);

    auto ta1_x_yyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 108);

    auto ta1_x_yyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 109);

    auto ta1_x_yyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 110);

    auto ta1_x_yyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 111);

    auto ta1_x_yyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 112);

    auto ta1_x_yyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 113);

    auto ta1_x_yzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 114);

    auto ta1_x_yzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 116);

    auto ta1_x_yzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 117);

    auto ta1_x_yzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 118);

    auto ta1_x_yzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 119);

    auto ta1_x_zzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 120);

    auto ta1_x_zzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 121);

    auto ta1_x_zzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 122);

    auto ta1_x_zzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 123);

    auto ta1_x_zzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 124);

    auto ta1_x_zzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 125);

    auto ta1_y_xxxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 126);

    auto ta1_y_xxxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 127);

    auto ta1_y_xxxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 128);

    auto ta1_y_xxxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 129);

    auto ta1_y_xxxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 130);

    auto ta1_y_xxxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 131);

    auto ta1_y_xxxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 132);

    auto ta1_y_xxxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 133);

    auto ta1_y_xxxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 134);

    auto ta1_y_xxxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 135);

    auto ta1_y_xxxxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 136);

    auto ta1_y_xxxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 138);

    auto ta1_y_xxxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 139);

    auto ta1_y_xxxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 140);

    auto ta1_y_xxxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 142);

    auto ta1_y_xxxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 143);

    auto ta1_y_xxxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 144);

    auto ta1_y_xxxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 145);

    auto ta1_y_xxxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 146);

    auto ta1_y_xxxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 147);

    auto ta1_y_xxxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 148);

    auto ta1_y_xxxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 149);

    auto ta1_y_xxxyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 151);

    auto ta1_y_xxxyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 154);

    auto ta1_y_xxxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 156);

    auto ta1_y_xxxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 157);

    auto ta1_y_xxxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 158);

    auto ta1_y_xxxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 159);

    auto ta1_y_xxxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 160);

    auto ta1_y_xxxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 161);

    auto ta1_y_xxyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 162);

    auto ta1_y_xxyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 163);

    auto ta1_y_xxyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 164);

    auto ta1_y_xxyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 165);

    auto ta1_y_xxyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 166);

    auto ta1_y_xxyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 167);

    auto ta1_y_xxyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 168);

    auto ta1_y_xxyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 169);

    auto ta1_y_xxyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 172);

    auto ta1_y_xxyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 173);

    auto ta1_y_xxyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 175);

    auto ta1_y_xxyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 176);

    auto ta1_y_xxyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 177);

    auto ta1_y_xxyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 178);

    auto ta1_y_xxzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 180);

    auto ta1_y_xxzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 181);

    auto ta1_y_xxzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 182);

    auto ta1_y_xxzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 183);

    auto ta1_y_xxzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 184);

    auto ta1_y_xxzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 185);

    auto ta1_y_xyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 186);

    auto ta1_y_xyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 187);

    auto ta1_y_xyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 189);

    auto ta1_y_xyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 190);

    auto ta1_y_xyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 191);

    auto ta1_y_xyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 196);

    auto ta1_y_xyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 197);

    auto ta1_y_xyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 201);

    auto ta1_y_xyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 202);

    auto ta1_y_xyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 203);

    auto ta1_y_xyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 207);

    auto ta1_y_xyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 208);

    auto ta1_y_xzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 210);

    auto ta1_y_xzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 212);

    auto ta1_y_xzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 213);

    auto ta1_y_xzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 214);

    auto ta1_y_xzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 215);

    auto ta1_y_yyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 216);

    auto ta1_y_yyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 217);

    auto ta1_y_yyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 218);

    auto ta1_y_yyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 219);

    auto ta1_y_yyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 220);

    auto ta1_y_yyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 221);

    auto ta1_y_yyyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 222);

    auto ta1_y_yyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 223);

    auto ta1_y_yyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 224);

    auto ta1_y_yyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 225);

    auto ta1_y_yyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 226);

    auto ta1_y_yyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 227);

    auto ta1_y_yyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 228);

    auto ta1_y_yyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 229);

    auto ta1_y_yyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 230);

    auto ta1_y_yyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 231);

    auto ta1_y_yyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 232);

    auto ta1_y_yyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 233);

    auto ta1_y_yyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 234);

    auto ta1_y_yyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 235);

    auto ta1_y_yyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 236);

    auto ta1_y_yyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 237);

    auto ta1_y_yyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 238);

    auto ta1_y_yyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 239);

    auto ta1_y_yzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 241);

    auto ta1_y_yzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 242);

    auto ta1_y_yzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 243);

    auto ta1_y_yzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 244);

    auto ta1_y_yzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 245);

    auto ta1_y_zzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 246);

    auto ta1_y_zzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 247);

    auto ta1_y_zzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 248);

    auto ta1_y_zzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 249);

    auto ta1_y_zzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 250);

    auto ta1_y_zzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 251);

    auto ta1_z_xxxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 252);

    auto ta1_z_xxxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 253);

    auto ta1_z_xxxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 254);

    auto ta1_z_xxxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 255);

    auto ta1_z_xxxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 256);

    auto ta1_z_xxxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 257);

    auto ta1_z_xxxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 258);

    auto ta1_z_xxxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 259);

    auto ta1_z_xxxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 260);

    auto ta1_z_xxxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 261);

    auto ta1_z_xxxxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 262);

    auto ta1_z_xxxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 264);

    auto ta1_z_xxxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 265);

    auto ta1_z_xxxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 266);

    auto ta1_z_xxxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 268);

    auto ta1_z_xxxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 269);

    auto ta1_z_xxxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 270);

    auto ta1_z_xxxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 271);

    auto ta1_z_xxxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 272);

    auto ta1_z_xxxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 273);

    auto ta1_z_xxxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 274);

    auto ta1_z_xxxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 275);

    auto ta1_z_xxxyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 278);

    auto ta1_z_xxxyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 280);

    auto ta1_z_xxxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 282);

    auto ta1_z_xxxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 283);

    auto ta1_z_xxxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 284);

    auto ta1_z_xxxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 285);

    auto ta1_z_xxxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 286);

    auto ta1_z_xxxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 287);

    auto ta1_z_xxyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 288);

    auto ta1_z_xxyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 289);

    auto ta1_z_xxyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 290);

    auto ta1_z_xxyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 291);

    auto ta1_z_xxyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 292);

    auto ta1_z_xxyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 293);

    auto ta1_z_xxyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 295);

    auto ta1_z_xxyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 296);

    auto ta1_z_xxyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 298);

    auto ta1_z_xxyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 299);

    auto ta1_z_xxyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 300);

    auto ta1_z_xxyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 302);

    auto ta1_z_xxyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 303);

    auto ta1_z_xxyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 304);

    auto ta1_z_xxzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 306);

    auto ta1_z_xxzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 307);

    auto ta1_z_xxzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 308);

    auto ta1_z_xxzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 309);

    auto ta1_z_xxzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 310);

    auto ta1_z_xxzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 311);

    auto ta1_z_xyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 312);

    auto ta1_z_xyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 313);

    auto ta1_z_xyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 315);

    auto ta1_z_xyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 316);

    auto ta1_z_xyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 317);

    auto ta1_z_xyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 322);

    auto ta1_z_xyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 323);

    auto ta1_z_xyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 327);

    auto ta1_z_xyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 328);

    auto ta1_z_xyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 329);

    auto ta1_z_xyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 333);

    auto ta1_z_xyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 334);

    auto ta1_z_xzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 336);

    auto ta1_z_xzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 338);

    auto ta1_z_xzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 339);

    auto ta1_z_xzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 340);

    auto ta1_z_xzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 341);

    auto ta1_z_yyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 342);

    auto ta1_z_yyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 343);

    auto ta1_z_yyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 344);

    auto ta1_z_yyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 345);

    auto ta1_z_yyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 346);

    auto ta1_z_yyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 347);

    auto ta1_z_yyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 349);

    auto ta1_z_yyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 350);

    auto ta1_z_yyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 351);

    auto ta1_z_yyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 352);

    auto ta1_z_yyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 353);

    auto ta1_z_yyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 354);

    auto ta1_z_yyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 355);

    auto ta1_z_yyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 356);

    auto ta1_z_yyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 357);

    auto ta1_z_yyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 358);

    auto ta1_z_yyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 359);

    auto ta1_z_yyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 360);

    auto ta1_z_yyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 361);

    auto ta1_z_yyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 362);

    auto ta1_z_yyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 363);

    auto ta1_z_yyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 364);

    auto ta1_z_yyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 365);

    auto ta1_z_yzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 366);

    auto ta1_z_yzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 367);

    auto ta1_z_yzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 368);

    auto ta1_z_yzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 369);

    auto ta1_z_yzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 370);

    auto ta1_z_yzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 371);

    auto ta1_z_zzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_hd + 372);

    auto ta1_z_zzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 373);

    auto ta1_z_zzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 374);

    auto ta1_z_zzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_hd + 375);

    auto ta1_z_zzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 376);

    auto ta1_z_zzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_hd + 377);

    // Set up components of auxiliary buffer : HD

    auto ta1_x_xxxxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd);

    auto ta1_x_xxxxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 1);

    auto ta1_x_xxxxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 2);

    auto ta1_x_xxxxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 3);

    auto ta1_x_xxxxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 4);

    auto ta1_x_xxxxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 5);

    auto ta1_x_xxxxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 6);

    auto ta1_x_xxxxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 7);

    auto ta1_x_xxxxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 8);

    auto ta1_x_xxxxy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 9);

    auto ta1_x_xxxxy_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 11);

    auto ta1_x_xxxxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 12);

    auto ta1_x_xxxxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 13);

    auto ta1_x_xxxxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 14);

    auto ta1_x_xxxxz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 15);

    auto ta1_x_xxxxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 16);

    auto ta1_x_xxxxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 17);

    auto ta1_x_xxxyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 18);

    auto ta1_x_xxxyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 19);

    auto ta1_x_xxxyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 20);

    auto ta1_x_xxxyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 21);

    auto ta1_x_xxxyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 22);

    auto ta1_x_xxxyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 23);

    auto ta1_x_xxxyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 26);

    auto ta1_x_xxxyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 29);

    auto ta1_x_xxxzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 30);

    auto ta1_x_xxxzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 31);

    auto ta1_x_xxxzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 32);

    auto ta1_x_xxxzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 33);

    auto ta1_x_xxxzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 34);

    auto ta1_x_xxxzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 35);

    auto ta1_x_xxyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 36);

    auto ta1_x_xxyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 37);

    auto ta1_x_xxyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 38);

    auto ta1_x_xxyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 39);

    auto ta1_x_xxyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 40);

    auto ta1_x_xxyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 41);

    auto ta1_x_xxyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 43);

    auto ta1_x_xxyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 44);

    auto ta1_x_xxyyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 45);

    auto ta1_x_xxyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 47);

    auto ta1_x_xxyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 48);

    auto ta1_x_xxyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 50);

    auto ta1_x_xxyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 53);

    auto ta1_x_xxzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 54);

    auto ta1_x_xxzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 55);

    auto ta1_x_xxzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 56);

    auto ta1_x_xxzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 57);

    auto ta1_x_xxzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 58);

    auto ta1_x_xxzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 59);

    auto ta1_x_xyyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 60);

    auto ta1_x_xyyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 61);

    auto ta1_x_xyyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 62);

    auto ta1_x_xyyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 63);

    auto ta1_x_xyyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 64);

    auto ta1_x_xyyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 67);

    auto ta1_x_xyyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 68);

    auto ta1_x_xyyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 72);

    auto ta1_x_xyyzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 73);

    auto ta1_x_xyyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 74);

    auto ta1_x_xyyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 76);

    auto ta1_x_xyzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 78);

    auto ta1_x_xyzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 80);

    auto ta1_x_xzzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 84);

    auto ta1_x_xzzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 85);

    auto ta1_x_xzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 86);

    auto ta1_x_xzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 88);

    auto ta1_x_xzzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 89);

    auto ta1_x_yyyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 90);

    auto ta1_x_yyyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 91);

    auto ta1_x_yyyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 92);

    auto ta1_x_yyyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 93);

    auto ta1_x_yyyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 94);

    auto ta1_x_yyyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 95);

    auto ta1_x_yyyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 97);

    auto ta1_x_yyyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 98);

    auto ta1_x_yyyyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 99);

    auto ta1_x_yyyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 100);

    auto ta1_x_yyyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 101);

    auto ta1_x_yyyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 102);

    auto ta1_x_yyyzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 103);

    auto ta1_x_yyyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 104);

    auto ta1_x_yyyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 105);

    auto ta1_x_yyyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 106);

    auto ta1_x_yyyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 107);

    auto ta1_x_yyzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 108);

    auto ta1_x_yyzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 109);

    auto ta1_x_yyzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 110);

    auto ta1_x_yyzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 111);

    auto ta1_x_yyzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 112);

    auto ta1_x_yyzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 113);

    auto ta1_x_yzzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 114);

    auto ta1_x_yzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 116);

    auto ta1_x_yzzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 117);

    auto ta1_x_yzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 118);

    auto ta1_x_yzzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 119);

    auto ta1_x_zzzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 120);

    auto ta1_x_zzzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 121);

    auto ta1_x_zzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 122);

    auto ta1_x_zzzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 123);

    auto ta1_x_zzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 124);

    auto ta1_x_zzzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 125);

    auto ta1_y_xxxxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 126);

    auto ta1_y_xxxxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 127);

    auto ta1_y_xxxxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 128);

    auto ta1_y_xxxxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 129);

    auto ta1_y_xxxxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 130);

    auto ta1_y_xxxxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 131);

    auto ta1_y_xxxxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 132);

    auto ta1_y_xxxxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 133);

    auto ta1_y_xxxxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 134);

    auto ta1_y_xxxxy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 135);

    auto ta1_y_xxxxy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 136);

    auto ta1_y_xxxxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 138);

    auto ta1_y_xxxxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 139);

    auto ta1_y_xxxxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 140);

    auto ta1_y_xxxxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 142);

    auto ta1_y_xxxxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 143);

    auto ta1_y_xxxyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 144);

    auto ta1_y_xxxyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 145);

    auto ta1_y_xxxyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 146);

    auto ta1_y_xxxyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 147);

    auto ta1_y_xxxyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 148);

    auto ta1_y_xxxyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 149);

    auto ta1_y_xxxyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 151);

    auto ta1_y_xxxyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 154);

    auto ta1_y_xxxzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 156);

    auto ta1_y_xxxzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 157);

    auto ta1_y_xxxzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 158);

    auto ta1_y_xxxzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 159);

    auto ta1_y_xxxzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 160);

    auto ta1_y_xxxzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 161);

    auto ta1_y_xxyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 162);

    auto ta1_y_xxyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 163);

    auto ta1_y_xxyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 164);

    auto ta1_y_xxyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 165);

    auto ta1_y_xxyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 166);

    auto ta1_y_xxyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 167);

    auto ta1_y_xxyyz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 168);

    auto ta1_y_xxyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 169);

    auto ta1_y_xxyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 172);

    auto ta1_y_xxyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 173);

    auto ta1_y_xxyzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 175);

    auto ta1_y_xxyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 176);

    auto ta1_y_xxyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 177);

    auto ta1_y_xxyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 178);

    auto ta1_y_xxzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 180);

    auto ta1_y_xxzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 181);

    auto ta1_y_xxzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 182);

    auto ta1_y_xxzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 183);

    auto ta1_y_xxzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 184);

    auto ta1_y_xxzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 185);

    auto ta1_y_xyyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 186);

    auto ta1_y_xyyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 187);

    auto ta1_y_xyyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 189);

    auto ta1_y_xyyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 190);

    auto ta1_y_xyyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 191);

    auto ta1_y_xyyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 196);

    auto ta1_y_xyyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 197);

    auto ta1_y_xyyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 201);

    auto ta1_y_xyyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 202);

    auto ta1_y_xyyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 203);

    auto ta1_y_xyzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 207);

    auto ta1_y_xyzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 208);

    auto ta1_y_xzzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 210);

    auto ta1_y_xzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 212);

    auto ta1_y_xzzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 213);

    auto ta1_y_xzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 214);

    auto ta1_y_xzzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 215);

    auto ta1_y_yyyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 216);

    auto ta1_y_yyyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 217);

    auto ta1_y_yyyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 218);

    auto ta1_y_yyyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 219);

    auto ta1_y_yyyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 220);

    auto ta1_y_yyyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 221);

    auto ta1_y_yyyyz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 222);

    auto ta1_y_yyyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 223);

    auto ta1_y_yyyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 224);

    auto ta1_y_yyyyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 225);

    auto ta1_y_yyyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 226);

    auto ta1_y_yyyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 227);

    auto ta1_y_yyyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 228);

    auto ta1_y_yyyzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 229);

    auto ta1_y_yyyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 230);

    auto ta1_y_yyyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 231);

    auto ta1_y_yyyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 232);

    auto ta1_y_yyyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 233);

    auto ta1_y_yyzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 234);

    auto ta1_y_yyzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 235);

    auto ta1_y_yyzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 236);

    auto ta1_y_yyzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 237);

    auto ta1_y_yyzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 238);

    auto ta1_y_yyzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 239);

    auto ta1_y_yzzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 241);

    auto ta1_y_yzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 242);

    auto ta1_y_yzzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 243);

    auto ta1_y_yzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 244);

    auto ta1_y_yzzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 245);

    auto ta1_y_zzzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 246);

    auto ta1_y_zzzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 247);

    auto ta1_y_zzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 248);

    auto ta1_y_zzzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 249);

    auto ta1_y_zzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 250);

    auto ta1_y_zzzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 251);

    auto ta1_z_xxxxx_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 252);

    auto ta1_z_xxxxx_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 253);

    auto ta1_z_xxxxx_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 254);

    auto ta1_z_xxxxx_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 255);

    auto ta1_z_xxxxx_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 256);

    auto ta1_z_xxxxx_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 257);

    auto ta1_z_xxxxy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 258);

    auto ta1_z_xxxxy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 259);

    auto ta1_z_xxxxy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 260);

    auto ta1_z_xxxxy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 261);

    auto ta1_z_xxxxy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 262);

    auto ta1_z_xxxxz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 264);

    auto ta1_z_xxxxz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 265);

    auto ta1_z_xxxxz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 266);

    auto ta1_z_xxxxz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 268);

    auto ta1_z_xxxxz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 269);

    auto ta1_z_xxxyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 270);

    auto ta1_z_xxxyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 271);

    auto ta1_z_xxxyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 272);

    auto ta1_z_xxxyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 273);

    auto ta1_z_xxxyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 274);

    auto ta1_z_xxxyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 275);

    auto ta1_z_xxxyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 278);

    auto ta1_z_xxxyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 280);

    auto ta1_z_xxxzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 282);

    auto ta1_z_xxxzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 283);

    auto ta1_z_xxxzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 284);

    auto ta1_z_xxxzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 285);

    auto ta1_z_xxxzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 286);

    auto ta1_z_xxxzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 287);

    auto ta1_z_xxyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 288);

    auto ta1_z_xxyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 289);

    auto ta1_z_xxyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 290);

    auto ta1_z_xxyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 291);

    auto ta1_z_xxyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 292);

    auto ta1_z_xxyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 293);

    auto ta1_z_xxyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 295);

    auto ta1_z_xxyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 296);

    auto ta1_z_xxyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 298);

    auto ta1_z_xxyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 299);

    auto ta1_z_xxyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 300);

    auto ta1_z_xxyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 302);

    auto ta1_z_xxyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 303);

    auto ta1_z_xxyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 304);

    auto ta1_z_xxzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 306);

    auto ta1_z_xxzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 307);

    auto ta1_z_xxzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 308);

    auto ta1_z_xxzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 309);

    auto ta1_z_xxzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 310);

    auto ta1_z_xxzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 311);

    auto ta1_z_xyyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 312);

    auto ta1_z_xyyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 313);

    auto ta1_z_xyyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 315);

    auto ta1_z_xyyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 316);

    auto ta1_z_xyyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 317);

    auto ta1_z_xyyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 322);

    auto ta1_z_xyyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 323);

    auto ta1_z_xyyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 327);

    auto ta1_z_xyyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 328);

    auto ta1_z_xyyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 329);

    auto ta1_z_xyzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 333);

    auto ta1_z_xyzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 334);

    auto ta1_z_xzzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 336);

    auto ta1_z_xzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 338);

    auto ta1_z_xzzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 339);

    auto ta1_z_xzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 340);

    auto ta1_z_xzzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 341);

    auto ta1_z_yyyyy_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 342);

    auto ta1_z_yyyyy_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 343);

    auto ta1_z_yyyyy_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 344);

    auto ta1_z_yyyyy_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 345);

    auto ta1_z_yyyyy_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 346);

    auto ta1_z_yyyyy_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 347);

    auto ta1_z_yyyyz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 349);

    auto ta1_z_yyyyz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 350);

    auto ta1_z_yyyyz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 351);

    auto ta1_z_yyyyz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 352);

    auto ta1_z_yyyyz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 353);

    auto ta1_z_yyyzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 354);

    auto ta1_z_yyyzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 355);

    auto ta1_z_yyyzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 356);

    auto ta1_z_yyyzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 357);

    auto ta1_z_yyyzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 358);

    auto ta1_z_yyyzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 359);

    auto ta1_z_yyzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 360);

    auto ta1_z_yyzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 361);

    auto ta1_z_yyzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 362);

    auto ta1_z_yyzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 363);

    auto ta1_z_yyzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 364);

    auto ta1_z_yyzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 365);

    auto ta1_z_yzzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 366);

    auto ta1_z_yzzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 367);

    auto ta1_z_yzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 368);

    auto ta1_z_yzzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 369);

    auto ta1_z_yzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 370);

    auto ta1_z_yzzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 371);

    auto ta1_z_zzzzz_xx_1 = pbuffer.data(idx_npot_geom_010_1_hd + 372);

    auto ta1_z_zzzzz_xy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 373);

    auto ta1_z_zzzzz_xz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 374);

    auto ta1_z_zzzzz_yy_1 = pbuffer.data(idx_npot_geom_010_1_hd + 375);

    auto ta1_z_zzzzz_yz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 376);

    auto ta1_z_zzzzz_zz_1 = pbuffer.data(idx_npot_geom_010_1_hd + 377);

    // Set up 0-6 components of targeted buffer : ID

    auto ta1_x_xxxxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_id);

    auto ta1_x_xxxxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 1);

    auto ta1_x_xxxxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 2);

    auto ta1_x_xxxxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 3);

    auto ta1_x_xxxxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 4);

    auto ta1_x_xxxxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 5);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_x_xxxx_xx_0,   \
                             ta1_x_xxxx_xx_1,   \
                             ta1_x_xxxx_xy_0,   \
                             ta1_x_xxxx_xy_1,   \
                             ta1_x_xxxx_xz_0,   \
                             ta1_x_xxxx_xz_1,   \
                             ta1_x_xxxx_yy_0,   \
                             ta1_x_xxxx_yy_1,   \
                             ta1_x_xxxx_yz_0,   \
                             ta1_x_xxxx_yz_1,   \
                             ta1_x_xxxx_zz_0,   \
                             ta1_x_xxxx_zz_1,   \
                             ta1_x_xxxxx_x_0,   \
                             ta1_x_xxxxx_x_1,   \
                             ta1_x_xxxxx_xx_0,  \
                             ta1_x_xxxxx_xx_1,  \
                             ta1_x_xxxxx_xy_0,  \
                             ta1_x_xxxxx_xy_1,  \
                             ta1_x_xxxxx_xz_0,  \
                             ta1_x_xxxxx_xz_1,  \
                             ta1_x_xxxxx_y_0,   \
                             ta1_x_xxxxx_y_1,   \
                             ta1_x_xxxxx_yy_0,  \
                             ta1_x_xxxxx_yy_1,  \
                             ta1_x_xxxxx_yz_0,  \
                             ta1_x_xxxxx_yz_1,  \
                             ta1_x_xxxxx_z_0,   \
                             ta1_x_xxxxx_z_1,   \
                             ta1_x_xxxxx_zz_0,  \
                             ta1_x_xxxxx_zz_1,  \
                             ta1_x_xxxxxx_xx_0, \
                             ta1_x_xxxxxx_xy_0, \
                             ta1_x_xxxxxx_xz_0, \
                             ta1_x_xxxxxx_yy_0, \
                             ta1_x_xxxxxx_yz_0, \
                             ta1_x_xxxxxx_zz_0, \
                             ta_xxxxx_xx_1,     \
                             ta_xxxxx_xy_1,     \
                             ta_xxxxx_xz_1,     \
                             ta_xxxxx_yy_1,     \
                             ta_xxxxx_yz_1,     \
                             ta_xxxxx_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxxx_xx_0[i] = 5.0 * ta1_x_xxxx_xx_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xx_1[i] * fe_0 + 2.0 * ta1_x_xxxxx_x_0[i] * fe_0 -
                               2.0 * ta1_x_xxxxx_x_1[i] * fe_0 + ta_xxxxx_xx_1[i] + ta1_x_xxxxx_xx_0[i] * pa_x[i] - ta1_x_xxxxx_xx_1[i] * pc_x[i];

        ta1_x_xxxxxx_xy_0[i] = 5.0 * ta1_x_xxxx_xy_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xy_1[i] * fe_0 + ta1_x_xxxxx_y_0[i] * fe_0 -
                               ta1_x_xxxxx_y_1[i] * fe_0 + ta_xxxxx_xy_1[i] + ta1_x_xxxxx_xy_0[i] * pa_x[i] - ta1_x_xxxxx_xy_1[i] * pc_x[i];

        ta1_x_xxxxxx_xz_0[i] = 5.0 * ta1_x_xxxx_xz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_xz_1[i] * fe_0 + ta1_x_xxxxx_z_0[i] * fe_0 -
                               ta1_x_xxxxx_z_1[i] * fe_0 + ta_xxxxx_xz_1[i] + ta1_x_xxxxx_xz_0[i] * pa_x[i] - ta1_x_xxxxx_xz_1[i] * pc_x[i];

        ta1_x_xxxxxx_yy_0[i] = 5.0 * ta1_x_xxxx_yy_0[i] * fe_0 - 5.0 * ta1_x_xxxx_yy_1[i] * fe_0 + ta_xxxxx_yy_1[i] + ta1_x_xxxxx_yy_0[i] * pa_x[i] -
                               ta1_x_xxxxx_yy_1[i] * pc_x[i];

        ta1_x_xxxxxx_yz_0[i] = 5.0 * ta1_x_xxxx_yz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_yz_1[i] * fe_0 + ta_xxxxx_yz_1[i] + ta1_x_xxxxx_yz_0[i] * pa_x[i] -
                               ta1_x_xxxxx_yz_1[i] * pc_x[i];

        ta1_x_xxxxxx_zz_0[i] = 5.0 * ta1_x_xxxx_zz_0[i] * fe_0 - 5.0 * ta1_x_xxxx_zz_1[i] * fe_0 + ta_xxxxx_zz_1[i] + ta1_x_xxxxx_zz_0[i] * pa_x[i] -
                               ta1_x_xxxxx_zz_1[i] * pc_x[i];
    }

    // Set up 6-12 components of targeted buffer : ID

    auto ta1_x_xxxxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 6);

    auto ta1_x_xxxxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 7);

    auto ta1_x_xxxxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 8);

    auto ta1_x_xxxxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 9);

    auto ta1_x_xxxxxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 10);

    auto ta1_x_xxxxxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 11);

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_x_xxxxx_x_0,   \
                             ta1_x_xxxxx_x_1,   \
                             ta1_x_xxxxx_xx_0,  \
                             ta1_x_xxxxx_xx_1,  \
                             ta1_x_xxxxx_xy_0,  \
                             ta1_x_xxxxx_xy_1,  \
                             ta1_x_xxxxx_xz_0,  \
                             ta1_x_xxxxx_xz_1,  \
                             ta1_x_xxxxx_y_0,   \
                             ta1_x_xxxxx_y_1,   \
                             ta1_x_xxxxx_yy_0,  \
                             ta1_x_xxxxx_yy_1,  \
                             ta1_x_xxxxx_yz_0,  \
                             ta1_x_xxxxx_yz_1,  \
                             ta1_x_xxxxx_z_0,   \
                             ta1_x_xxxxx_z_1,   \
                             ta1_x_xxxxx_zz_0,  \
                             ta1_x_xxxxx_zz_1,  \
                             ta1_x_xxxxxy_xx_0, \
                             ta1_x_xxxxxy_xy_0, \
                             ta1_x_xxxxxy_xz_0, \
                             ta1_x_xxxxxy_yy_0, \
                             ta1_x_xxxxxy_yz_0, \
                             ta1_x_xxxxxy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxxy_xx_0[i] = ta1_x_xxxxx_xx_0[i] * pa_y[i] - ta1_x_xxxxx_xx_1[i] * pc_y[i];

        ta1_x_xxxxxy_xy_0[i] = ta1_x_xxxxx_x_0[i] * fe_0 - ta1_x_xxxxx_x_1[i] * fe_0 + ta1_x_xxxxx_xy_0[i] * pa_y[i] - ta1_x_xxxxx_xy_1[i] * pc_y[i];

        ta1_x_xxxxxy_xz_0[i] = ta1_x_xxxxx_xz_0[i] * pa_y[i] - ta1_x_xxxxx_xz_1[i] * pc_y[i];

        ta1_x_xxxxxy_yy_0[i] =
            2.0 * ta1_x_xxxxx_y_0[i] * fe_0 - 2.0 * ta1_x_xxxxx_y_1[i] * fe_0 + ta1_x_xxxxx_yy_0[i] * pa_y[i] - ta1_x_xxxxx_yy_1[i] * pc_y[i];

        ta1_x_xxxxxy_yz_0[i] = ta1_x_xxxxx_z_0[i] * fe_0 - ta1_x_xxxxx_z_1[i] * fe_0 + ta1_x_xxxxx_yz_0[i] * pa_y[i] - ta1_x_xxxxx_yz_1[i] * pc_y[i];

        ta1_x_xxxxxy_zz_0[i] = ta1_x_xxxxx_zz_0[i] * pa_y[i] - ta1_x_xxxxx_zz_1[i] * pc_y[i];
    }

    // Set up 12-18 components of targeted buffer : ID

    auto ta1_x_xxxxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 12);

    auto ta1_x_xxxxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 13);

    auto ta1_x_xxxxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 14);

    auto ta1_x_xxxxxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 15);

    auto ta1_x_xxxxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 16);

    auto ta1_x_xxxxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 17);

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta1_x_xxxxx_x_0,   \
                             ta1_x_xxxxx_x_1,   \
                             ta1_x_xxxxx_xx_0,  \
                             ta1_x_xxxxx_xx_1,  \
                             ta1_x_xxxxx_xy_0,  \
                             ta1_x_xxxxx_xy_1,  \
                             ta1_x_xxxxx_xz_0,  \
                             ta1_x_xxxxx_xz_1,  \
                             ta1_x_xxxxx_y_0,   \
                             ta1_x_xxxxx_y_1,   \
                             ta1_x_xxxxx_yy_0,  \
                             ta1_x_xxxxx_yy_1,  \
                             ta1_x_xxxxx_yz_0,  \
                             ta1_x_xxxxx_yz_1,  \
                             ta1_x_xxxxx_z_0,   \
                             ta1_x_xxxxx_z_1,   \
                             ta1_x_xxxxx_zz_0,  \
                             ta1_x_xxxxx_zz_1,  \
                             ta1_x_xxxxxz_xx_0, \
                             ta1_x_xxxxxz_xy_0, \
                             ta1_x_xxxxxz_xz_0, \
                             ta1_x_xxxxxz_yy_0, \
                             ta1_x_xxxxxz_yz_0, \
                             ta1_x_xxxxxz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxxz_xx_0[i] = ta1_x_xxxxx_xx_0[i] * pa_z[i] - ta1_x_xxxxx_xx_1[i] * pc_z[i];

        ta1_x_xxxxxz_xy_0[i] = ta1_x_xxxxx_xy_0[i] * pa_z[i] - ta1_x_xxxxx_xy_1[i] * pc_z[i];

        ta1_x_xxxxxz_xz_0[i] = ta1_x_xxxxx_x_0[i] * fe_0 - ta1_x_xxxxx_x_1[i] * fe_0 + ta1_x_xxxxx_xz_0[i] * pa_z[i] - ta1_x_xxxxx_xz_1[i] * pc_z[i];

        ta1_x_xxxxxz_yy_0[i] = ta1_x_xxxxx_yy_0[i] * pa_z[i] - ta1_x_xxxxx_yy_1[i] * pc_z[i];

        ta1_x_xxxxxz_yz_0[i] = ta1_x_xxxxx_y_0[i] * fe_0 - ta1_x_xxxxx_y_1[i] * fe_0 + ta1_x_xxxxx_yz_0[i] * pa_z[i] - ta1_x_xxxxx_yz_1[i] * pc_z[i];

        ta1_x_xxxxxz_zz_0[i] =
            2.0 * ta1_x_xxxxx_z_0[i] * fe_0 - 2.0 * ta1_x_xxxxx_z_1[i] * fe_0 + ta1_x_xxxxx_zz_0[i] * pa_z[i] - ta1_x_xxxxx_zz_1[i] * pc_z[i];
    }

    // Set up 18-24 components of targeted buffer : ID

    auto ta1_x_xxxxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 18);

    auto ta1_x_xxxxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 19);

    auto ta1_x_xxxxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 20);

    auto ta1_x_xxxxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 21);

    auto ta1_x_xxxxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 22);

    auto ta1_x_xxxxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 23);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_x_xxxx_xx_0,   \
                             ta1_x_xxxx_xx_1,   \
                             ta1_x_xxxx_xy_0,   \
                             ta1_x_xxxx_xy_1,   \
                             ta1_x_xxxx_xz_0,   \
                             ta1_x_xxxx_xz_1,   \
                             ta1_x_xxxx_zz_0,   \
                             ta1_x_xxxx_zz_1,   \
                             ta1_x_xxxxy_x_0,   \
                             ta1_x_xxxxy_x_1,   \
                             ta1_x_xxxxy_xx_0,  \
                             ta1_x_xxxxy_xx_1,  \
                             ta1_x_xxxxy_xy_0,  \
                             ta1_x_xxxxy_xy_1,  \
                             ta1_x_xxxxy_xz_0,  \
                             ta1_x_xxxxy_xz_1,  \
                             ta1_x_xxxxy_zz_0,  \
                             ta1_x_xxxxy_zz_1,  \
                             ta1_x_xxxxyy_xx_0, \
                             ta1_x_xxxxyy_xy_0, \
                             ta1_x_xxxxyy_xz_0, \
                             ta1_x_xxxxyy_yy_0, \
                             ta1_x_xxxxyy_yz_0, \
                             ta1_x_xxxxyy_zz_0, \
                             ta1_x_xxxyy_yy_0,  \
                             ta1_x_xxxyy_yy_1,  \
                             ta1_x_xxxyy_yz_0,  \
                             ta1_x_xxxyy_yz_1,  \
                             ta1_x_xxyy_yy_0,   \
                             ta1_x_xxyy_yy_1,   \
                             ta1_x_xxyy_yz_0,   \
                             ta1_x_xxyy_yz_1,   \
                             ta_xxxyy_yy_1,     \
                             ta_xxxyy_yz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxyy_xx_0[i] = ta1_x_xxxx_xx_0[i] * fe_0 - ta1_x_xxxx_xx_1[i] * fe_0 + ta1_x_xxxxy_xx_0[i] * pa_y[i] - ta1_x_xxxxy_xx_1[i] * pc_y[i];

        ta1_x_xxxxyy_xy_0[i] = ta1_x_xxxx_xy_0[i] * fe_0 - ta1_x_xxxx_xy_1[i] * fe_0 + ta1_x_xxxxy_x_0[i] * fe_0 - ta1_x_xxxxy_x_1[i] * fe_0 +
                               ta1_x_xxxxy_xy_0[i] * pa_y[i] - ta1_x_xxxxy_xy_1[i] * pc_y[i];

        ta1_x_xxxxyy_xz_0[i] = ta1_x_xxxx_xz_0[i] * fe_0 - ta1_x_xxxx_xz_1[i] * fe_0 + ta1_x_xxxxy_xz_0[i] * pa_y[i] - ta1_x_xxxxy_xz_1[i] * pc_y[i];

        ta1_x_xxxxyy_yy_0[i] = 3.0 * ta1_x_xxyy_yy_0[i] * fe_0 - 3.0 * ta1_x_xxyy_yy_1[i] * fe_0 + ta_xxxyy_yy_1[i] + ta1_x_xxxyy_yy_0[i] * pa_x[i] -
                               ta1_x_xxxyy_yy_1[i] * pc_x[i];

        ta1_x_xxxxyy_yz_0[i] = 3.0 * ta1_x_xxyy_yz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_yz_1[i] * fe_0 + ta_xxxyy_yz_1[i] + ta1_x_xxxyy_yz_0[i] * pa_x[i] -
                               ta1_x_xxxyy_yz_1[i] * pc_x[i];

        ta1_x_xxxxyy_zz_0[i] = ta1_x_xxxx_zz_0[i] * fe_0 - ta1_x_xxxx_zz_1[i] * fe_0 + ta1_x_xxxxy_zz_0[i] * pa_y[i] - ta1_x_xxxxy_zz_1[i] * pc_y[i];
    }

    // Set up 24-30 components of targeted buffer : ID

    auto ta1_x_xxxxyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 24);

    auto ta1_x_xxxxyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 25);

    auto ta1_x_xxxxyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 26);

    auto ta1_x_xxxxyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 27);

    auto ta1_x_xxxxyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 28);

    auto ta1_x_xxxxyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 29);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_xxxxy_xy_0,  \
                             ta1_x_xxxxy_xy_1,  \
                             ta1_x_xxxxy_yy_0,  \
                             ta1_x_xxxxy_yy_1,  \
                             ta1_x_xxxxyz_xx_0, \
                             ta1_x_xxxxyz_xy_0, \
                             ta1_x_xxxxyz_xz_0, \
                             ta1_x_xxxxyz_yy_0, \
                             ta1_x_xxxxyz_yz_0, \
                             ta1_x_xxxxyz_zz_0, \
                             ta1_x_xxxxz_xx_0,  \
                             ta1_x_xxxxz_xx_1,  \
                             ta1_x_xxxxz_xz_0,  \
                             ta1_x_xxxxz_xz_1,  \
                             ta1_x_xxxxz_yz_0,  \
                             ta1_x_xxxxz_yz_1,  \
                             ta1_x_xxxxz_z_0,   \
                             ta1_x_xxxxz_z_1,   \
                             ta1_x_xxxxz_zz_0,  \
                             ta1_x_xxxxz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxyz_xx_0[i] = ta1_x_xxxxz_xx_0[i] * pa_y[i] - ta1_x_xxxxz_xx_1[i] * pc_y[i];

        ta1_x_xxxxyz_xy_0[i] = ta1_x_xxxxy_xy_0[i] * pa_z[i] - ta1_x_xxxxy_xy_1[i] * pc_z[i];

        ta1_x_xxxxyz_xz_0[i] = ta1_x_xxxxz_xz_0[i] * pa_y[i] - ta1_x_xxxxz_xz_1[i] * pc_y[i];

        ta1_x_xxxxyz_yy_0[i] = ta1_x_xxxxy_yy_0[i] * pa_z[i] - ta1_x_xxxxy_yy_1[i] * pc_z[i];

        ta1_x_xxxxyz_yz_0[i] = ta1_x_xxxxz_z_0[i] * fe_0 - ta1_x_xxxxz_z_1[i] * fe_0 + ta1_x_xxxxz_yz_0[i] * pa_y[i] - ta1_x_xxxxz_yz_1[i] * pc_y[i];

        ta1_x_xxxxyz_zz_0[i] = ta1_x_xxxxz_zz_0[i] * pa_y[i] - ta1_x_xxxxz_zz_1[i] * pc_y[i];
    }

    // Set up 30-36 components of targeted buffer : ID

    auto ta1_x_xxxxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 30);

    auto ta1_x_xxxxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 31);

    auto ta1_x_xxxxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 32);

    auto ta1_x_xxxxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 33);

    auto ta1_x_xxxxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 34);

    auto ta1_x_xxxxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 35);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_x_xxxx_xx_0,   \
                             ta1_x_xxxx_xx_1,   \
                             ta1_x_xxxx_xy_0,   \
                             ta1_x_xxxx_xy_1,   \
                             ta1_x_xxxx_xz_0,   \
                             ta1_x_xxxx_xz_1,   \
                             ta1_x_xxxx_yy_0,   \
                             ta1_x_xxxx_yy_1,   \
                             ta1_x_xxxxz_x_0,   \
                             ta1_x_xxxxz_x_1,   \
                             ta1_x_xxxxz_xx_0,  \
                             ta1_x_xxxxz_xx_1,  \
                             ta1_x_xxxxz_xy_0,  \
                             ta1_x_xxxxz_xy_1,  \
                             ta1_x_xxxxz_xz_0,  \
                             ta1_x_xxxxz_xz_1,  \
                             ta1_x_xxxxz_yy_0,  \
                             ta1_x_xxxxz_yy_1,  \
                             ta1_x_xxxxzz_xx_0, \
                             ta1_x_xxxxzz_xy_0, \
                             ta1_x_xxxxzz_xz_0, \
                             ta1_x_xxxxzz_yy_0, \
                             ta1_x_xxxxzz_yz_0, \
                             ta1_x_xxxxzz_zz_0, \
                             ta1_x_xxxzz_yz_0,  \
                             ta1_x_xxxzz_yz_1,  \
                             ta1_x_xxxzz_zz_0,  \
                             ta1_x_xxxzz_zz_1,  \
                             ta1_x_xxzz_yz_0,   \
                             ta1_x_xxzz_yz_1,   \
                             ta1_x_xxzz_zz_0,   \
                             ta1_x_xxzz_zz_1,   \
                             ta_xxxzz_yz_1,     \
                             ta_xxxzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxxzz_xx_0[i] = ta1_x_xxxx_xx_0[i] * fe_0 - ta1_x_xxxx_xx_1[i] * fe_0 + ta1_x_xxxxz_xx_0[i] * pa_z[i] - ta1_x_xxxxz_xx_1[i] * pc_z[i];

        ta1_x_xxxxzz_xy_0[i] = ta1_x_xxxx_xy_0[i] * fe_0 - ta1_x_xxxx_xy_1[i] * fe_0 + ta1_x_xxxxz_xy_0[i] * pa_z[i] - ta1_x_xxxxz_xy_1[i] * pc_z[i];

        ta1_x_xxxxzz_xz_0[i] = ta1_x_xxxx_xz_0[i] * fe_0 - ta1_x_xxxx_xz_1[i] * fe_0 + ta1_x_xxxxz_x_0[i] * fe_0 - ta1_x_xxxxz_x_1[i] * fe_0 +
                               ta1_x_xxxxz_xz_0[i] * pa_z[i] - ta1_x_xxxxz_xz_1[i] * pc_z[i];

        ta1_x_xxxxzz_yy_0[i] = ta1_x_xxxx_yy_0[i] * fe_0 - ta1_x_xxxx_yy_1[i] * fe_0 + ta1_x_xxxxz_yy_0[i] * pa_z[i] - ta1_x_xxxxz_yy_1[i] * pc_z[i];

        ta1_x_xxxxzz_yz_0[i] = 3.0 * ta1_x_xxzz_yz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_yz_1[i] * fe_0 + ta_xxxzz_yz_1[i] + ta1_x_xxxzz_yz_0[i] * pa_x[i] -
                               ta1_x_xxxzz_yz_1[i] * pc_x[i];

        ta1_x_xxxxzz_zz_0[i] = 3.0 * ta1_x_xxzz_zz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_zz_1[i] * fe_0 + ta_xxxzz_zz_1[i] + ta1_x_xxxzz_zz_0[i] * pa_x[i] -
                               ta1_x_xxxzz_zz_1[i] * pc_x[i];
    }

    // Set up 36-42 components of targeted buffer : ID

    auto ta1_x_xxxyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 36);

    auto ta1_x_xxxyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 37);

    auto ta1_x_xxxyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 38);

    auto ta1_x_xxxyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 39);

    auto ta1_x_xxxyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 40);

    auto ta1_x_xxxyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 41);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_x_xxxy_xx_0,   \
                             ta1_x_xxxy_xx_1,   \
                             ta1_x_xxxy_xy_0,   \
                             ta1_x_xxxy_xy_1,   \
                             ta1_x_xxxy_xz_0,   \
                             ta1_x_xxxy_xz_1,   \
                             ta1_x_xxxy_zz_0,   \
                             ta1_x_xxxy_zz_1,   \
                             ta1_x_xxxyy_x_0,   \
                             ta1_x_xxxyy_x_1,   \
                             ta1_x_xxxyy_xx_0,  \
                             ta1_x_xxxyy_xx_1,  \
                             ta1_x_xxxyy_xy_0,  \
                             ta1_x_xxxyy_xy_1,  \
                             ta1_x_xxxyy_xz_0,  \
                             ta1_x_xxxyy_xz_1,  \
                             ta1_x_xxxyy_zz_0,  \
                             ta1_x_xxxyy_zz_1,  \
                             ta1_x_xxxyyy_xx_0, \
                             ta1_x_xxxyyy_xy_0, \
                             ta1_x_xxxyyy_xz_0, \
                             ta1_x_xxxyyy_yy_0, \
                             ta1_x_xxxyyy_yz_0, \
                             ta1_x_xxxyyy_zz_0, \
                             ta1_x_xxyyy_yy_0,  \
                             ta1_x_xxyyy_yy_1,  \
                             ta1_x_xxyyy_yz_0,  \
                             ta1_x_xxyyy_yz_1,  \
                             ta1_x_xyyy_yy_0,   \
                             ta1_x_xyyy_yy_1,   \
                             ta1_x_xyyy_yz_0,   \
                             ta1_x_xyyy_yz_1,   \
                             ta_xxyyy_yy_1,     \
                             ta_xxyyy_yz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyyy_xx_0[i] =
            2.0 * ta1_x_xxxy_xx_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xx_1[i] * fe_0 + ta1_x_xxxyy_xx_0[i] * pa_y[i] - ta1_x_xxxyy_xx_1[i] * pc_y[i];

        ta1_x_xxxyyy_xy_0[i] = 2.0 * ta1_x_xxxy_xy_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xy_1[i] * fe_0 + ta1_x_xxxyy_x_0[i] * fe_0 -
                               ta1_x_xxxyy_x_1[i] * fe_0 + ta1_x_xxxyy_xy_0[i] * pa_y[i] - ta1_x_xxxyy_xy_1[i] * pc_y[i];

        ta1_x_xxxyyy_xz_0[i] =
            2.0 * ta1_x_xxxy_xz_0[i] * fe_0 - 2.0 * ta1_x_xxxy_xz_1[i] * fe_0 + ta1_x_xxxyy_xz_0[i] * pa_y[i] - ta1_x_xxxyy_xz_1[i] * pc_y[i];

        ta1_x_xxxyyy_yy_0[i] = 2.0 * ta1_x_xyyy_yy_0[i] * fe_0 - 2.0 * ta1_x_xyyy_yy_1[i] * fe_0 + ta_xxyyy_yy_1[i] + ta1_x_xxyyy_yy_0[i] * pa_x[i] -
                               ta1_x_xxyyy_yy_1[i] * pc_x[i];

        ta1_x_xxxyyy_yz_0[i] = 2.0 * ta1_x_xyyy_yz_0[i] * fe_0 - 2.0 * ta1_x_xyyy_yz_1[i] * fe_0 + ta_xxyyy_yz_1[i] + ta1_x_xxyyy_yz_0[i] * pa_x[i] -
                               ta1_x_xxyyy_yz_1[i] * pc_x[i];

        ta1_x_xxxyyy_zz_0[i] =
            2.0 * ta1_x_xxxy_zz_0[i] * fe_0 - 2.0 * ta1_x_xxxy_zz_1[i] * fe_0 + ta1_x_xxxyy_zz_0[i] * pa_y[i] - ta1_x_xxxyy_zz_1[i] * pc_y[i];
    }

    // Set up 42-48 components of targeted buffer : ID

    auto ta1_x_xxxyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 42);

    auto ta1_x_xxxyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 43);

    auto ta1_x_xxxyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 44);

    auto ta1_x_xxxyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 45);

    auto ta1_x_xxxyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 46);

    auto ta1_x_xxxyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 47);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_xxxyy_xx_0,  \
                             ta1_x_xxxyy_xx_1,  \
                             ta1_x_xxxyy_xy_0,  \
                             ta1_x_xxxyy_xy_1,  \
                             ta1_x_xxxyy_y_0,   \
                             ta1_x_xxxyy_y_1,   \
                             ta1_x_xxxyy_yy_0,  \
                             ta1_x_xxxyy_yy_1,  \
                             ta1_x_xxxyy_yz_0,  \
                             ta1_x_xxxyy_yz_1,  \
                             ta1_x_xxxyyz_xx_0, \
                             ta1_x_xxxyyz_xy_0, \
                             ta1_x_xxxyyz_xz_0, \
                             ta1_x_xxxyyz_yy_0, \
                             ta1_x_xxxyyz_yz_0, \
                             ta1_x_xxxyyz_zz_0, \
                             ta1_x_xxxyz_xz_0,  \
                             ta1_x_xxxyz_xz_1,  \
                             ta1_x_xxxyz_zz_0,  \
                             ta1_x_xxxyz_zz_1,  \
                             ta1_x_xxxz_xz_0,   \
                             ta1_x_xxxz_xz_1,   \
                             ta1_x_xxxz_zz_0,   \
                             ta1_x_xxxz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyyz_xx_0[i] = ta1_x_xxxyy_xx_0[i] * pa_z[i] - ta1_x_xxxyy_xx_1[i] * pc_z[i];

        ta1_x_xxxyyz_xy_0[i] = ta1_x_xxxyy_xy_0[i] * pa_z[i] - ta1_x_xxxyy_xy_1[i] * pc_z[i];

        ta1_x_xxxyyz_xz_0[i] = ta1_x_xxxz_xz_0[i] * fe_0 - ta1_x_xxxz_xz_1[i] * fe_0 + ta1_x_xxxyz_xz_0[i] * pa_y[i] - ta1_x_xxxyz_xz_1[i] * pc_y[i];

        ta1_x_xxxyyz_yy_0[i] = ta1_x_xxxyy_yy_0[i] * pa_z[i] - ta1_x_xxxyy_yy_1[i] * pc_z[i];

        ta1_x_xxxyyz_yz_0[i] = ta1_x_xxxyy_y_0[i] * fe_0 - ta1_x_xxxyy_y_1[i] * fe_0 + ta1_x_xxxyy_yz_0[i] * pa_z[i] - ta1_x_xxxyy_yz_1[i] * pc_z[i];

        ta1_x_xxxyyz_zz_0[i] = ta1_x_xxxz_zz_0[i] * fe_0 - ta1_x_xxxz_zz_1[i] * fe_0 + ta1_x_xxxyz_zz_0[i] * pa_y[i] - ta1_x_xxxyz_zz_1[i] * pc_y[i];
    }

    // Set up 48-54 components of targeted buffer : ID

    auto ta1_x_xxxyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 48);

    auto ta1_x_xxxyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 49);

    auto ta1_x_xxxyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 50);

    auto ta1_x_xxxyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 51);

    auto ta1_x_xxxyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 52);

    auto ta1_x_xxxyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 53);

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_x_xxxyzz_xx_0, \
                             ta1_x_xxxyzz_xy_0, \
                             ta1_x_xxxyzz_xz_0, \
                             ta1_x_xxxyzz_yy_0, \
                             ta1_x_xxxyzz_yz_0, \
                             ta1_x_xxxyzz_zz_0, \
                             ta1_x_xxxzz_x_0,   \
                             ta1_x_xxxzz_x_1,   \
                             ta1_x_xxxzz_xx_0,  \
                             ta1_x_xxxzz_xx_1,  \
                             ta1_x_xxxzz_xy_0,  \
                             ta1_x_xxxzz_xy_1,  \
                             ta1_x_xxxzz_xz_0,  \
                             ta1_x_xxxzz_xz_1,  \
                             ta1_x_xxxzz_y_0,   \
                             ta1_x_xxxzz_y_1,   \
                             ta1_x_xxxzz_yy_0,  \
                             ta1_x_xxxzz_yy_1,  \
                             ta1_x_xxxzz_yz_0,  \
                             ta1_x_xxxzz_yz_1,  \
                             ta1_x_xxxzz_z_0,   \
                             ta1_x_xxxzz_z_1,   \
                             ta1_x_xxxzz_zz_0,  \
                             ta1_x_xxxzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxyzz_xx_0[i] = ta1_x_xxxzz_xx_0[i] * pa_y[i] - ta1_x_xxxzz_xx_1[i] * pc_y[i];

        ta1_x_xxxyzz_xy_0[i] = ta1_x_xxxzz_x_0[i] * fe_0 - ta1_x_xxxzz_x_1[i] * fe_0 + ta1_x_xxxzz_xy_0[i] * pa_y[i] - ta1_x_xxxzz_xy_1[i] * pc_y[i];

        ta1_x_xxxyzz_xz_0[i] = ta1_x_xxxzz_xz_0[i] * pa_y[i] - ta1_x_xxxzz_xz_1[i] * pc_y[i];

        ta1_x_xxxyzz_yy_0[i] =
            2.0 * ta1_x_xxxzz_y_0[i] * fe_0 - 2.0 * ta1_x_xxxzz_y_1[i] * fe_0 + ta1_x_xxxzz_yy_0[i] * pa_y[i] - ta1_x_xxxzz_yy_1[i] * pc_y[i];

        ta1_x_xxxyzz_yz_0[i] = ta1_x_xxxzz_z_0[i] * fe_0 - ta1_x_xxxzz_z_1[i] * fe_0 + ta1_x_xxxzz_yz_0[i] * pa_y[i] - ta1_x_xxxzz_yz_1[i] * pc_y[i];

        ta1_x_xxxyzz_zz_0[i] = ta1_x_xxxzz_zz_0[i] * pa_y[i] - ta1_x_xxxzz_zz_1[i] * pc_y[i];
    }

    // Set up 54-60 components of targeted buffer : ID

    auto ta1_x_xxxzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 54);

    auto ta1_x_xxxzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 55);

    auto ta1_x_xxxzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 56);

    auto ta1_x_xxxzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 57);

    auto ta1_x_xxxzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 58);

    auto ta1_x_xxxzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 59);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_x_xxxz_xx_0,   \
                             ta1_x_xxxz_xx_1,   \
                             ta1_x_xxxz_xy_0,   \
                             ta1_x_xxxz_xy_1,   \
                             ta1_x_xxxz_xz_0,   \
                             ta1_x_xxxz_xz_1,   \
                             ta1_x_xxxz_yy_0,   \
                             ta1_x_xxxz_yy_1,   \
                             ta1_x_xxxzz_x_0,   \
                             ta1_x_xxxzz_x_1,   \
                             ta1_x_xxxzz_xx_0,  \
                             ta1_x_xxxzz_xx_1,  \
                             ta1_x_xxxzz_xy_0,  \
                             ta1_x_xxxzz_xy_1,  \
                             ta1_x_xxxzz_xz_0,  \
                             ta1_x_xxxzz_xz_1,  \
                             ta1_x_xxxzz_yy_0,  \
                             ta1_x_xxxzz_yy_1,  \
                             ta1_x_xxxzzz_xx_0, \
                             ta1_x_xxxzzz_xy_0, \
                             ta1_x_xxxzzz_xz_0, \
                             ta1_x_xxxzzz_yy_0, \
                             ta1_x_xxxzzz_yz_0, \
                             ta1_x_xxxzzz_zz_0, \
                             ta1_x_xxzzz_yz_0,  \
                             ta1_x_xxzzz_yz_1,  \
                             ta1_x_xxzzz_zz_0,  \
                             ta1_x_xxzzz_zz_1,  \
                             ta1_x_xzzz_yz_0,   \
                             ta1_x_xzzz_yz_1,   \
                             ta1_x_xzzz_zz_0,   \
                             ta1_x_xzzz_zz_1,   \
                             ta_xxzzz_yz_1,     \
                             ta_xxzzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxzzz_xx_0[i] =
            2.0 * ta1_x_xxxz_xx_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xx_1[i] * fe_0 + ta1_x_xxxzz_xx_0[i] * pa_z[i] - ta1_x_xxxzz_xx_1[i] * pc_z[i];

        ta1_x_xxxzzz_xy_0[i] =
            2.0 * ta1_x_xxxz_xy_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xy_1[i] * fe_0 + ta1_x_xxxzz_xy_0[i] * pa_z[i] - ta1_x_xxxzz_xy_1[i] * pc_z[i];

        ta1_x_xxxzzz_xz_0[i] = 2.0 * ta1_x_xxxz_xz_0[i] * fe_0 - 2.0 * ta1_x_xxxz_xz_1[i] * fe_0 + ta1_x_xxxzz_x_0[i] * fe_0 -
                               ta1_x_xxxzz_x_1[i] * fe_0 + ta1_x_xxxzz_xz_0[i] * pa_z[i] - ta1_x_xxxzz_xz_1[i] * pc_z[i];

        ta1_x_xxxzzz_yy_0[i] =
            2.0 * ta1_x_xxxz_yy_0[i] * fe_0 - 2.0 * ta1_x_xxxz_yy_1[i] * fe_0 + ta1_x_xxxzz_yy_0[i] * pa_z[i] - ta1_x_xxxzz_yy_1[i] * pc_z[i];

        ta1_x_xxxzzz_yz_0[i] = 2.0 * ta1_x_xzzz_yz_0[i] * fe_0 - 2.0 * ta1_x_xzzz_yz_1[i] * fe_0 + ta_xxzzz_yz_1[i] + ta1_x_xxzzz_yz_0[i] * pa_x[i] -
                               ta1_x_xxzzz_yz_1[i] * pc_x[i];

        ta1_x_xxxzzz_zz_0[i] = 2.0 * ta1_x_xzzz_zz_0[i] * fe_0 - 2.0 * ta1_x_xzzz_zz_1[i] * fe_0 + ta_xxzzz_zz_1[i] + ta1_x_xxzzz_zz_0[i] * pa_x[i] -
                               ta1_x_xxzzz_zz_1[i] * pc_x[i];
    }

    // Set up 60-66 components of targeted buffer : ID

    auto ta1_x_xxyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 60);

    auto ta1_x_xxyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 61);

    auto ta1_x_xxyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 62);

    auto ta1_x_xxyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 63);

    auto ta1_x_xxyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 64);

    auto ta1_x_xxyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 65);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_x_xxyy_xx_0,   \
                             ta1_x_xxyy_xx_1,   \
                             ta1_x_xxyy_xy_0,   \
                             ta1_x_xxyy_xy_1,   \
                             ta1_x_xxyy_xz_0,   \
                             ta1_x_xxyy_xz_1,   \
                             ta1_x_xxyy_zz_0,   \
                             ta1_x_xxyy_zz_1,   \
                             ta1_x_xxyyy_x_0,   \
                             ta1_x_xxyyy_x_1,   \
                             ta1_x_xxyyy_xx_0,  \
                             ta1_x_xxyyy_xx_1,  \
                             ta1_x_xxyyy_xy_0,  \
                             ta1_x_xxyyy_xy_1,  \
                             ta1_x_xxyyy_xz_0,  \
                             ta1_x_xxyyy_xz_1,  \
                             ta1_x_xxyyy_zz_0,  \
                             ta1_x_xxyyy_zz_1,  \
                             ta1_x_xxyyyy_xx_0, \
                             ta1_x_xxyyyy_xy_0, \
                             ta1_x_xxyyyy_xz_0, \
                             ta1_x_xxyyyy_yy_0, \
                             ta1_x_xxyyyy_yz_0, \
                             ta1_x_xxyyyy_zz_0, \
                             ta1_x_xyyyy_yy_0,  \
                             ta1_x_xyyyy_yy_1,  \
                             ta1_x_xyyyy_yz_0,  \
                             ta1_x_xyyyy_yz_1,  \
                             ta1_x_yyyy_yy_0,   \
                             ta1_x_yyyy_yy_1,   \
                             ta1_x_yyyy_yz_0,   \
                             ta1_x_yyyy_yz_1,   \
                             ta_xyyyy_yy_1,     \
                             ta_xyyyy_yz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyyy_xx_0[i] =
            3.0 * ta1_x_xxyy_xx_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xx_1[i] * fe_0 + ta1_x_xxyyy_xx_0[i] * pa_y[i] - ta1_x_xxyyy_xx_1[i] * pc_y[i];

        ta1_x_xxyyyy_xy_0[i] = 3.0 * ta1_x_xxyy_xy_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xy_1[i] * fe_0 + ta1_x_xxyyy_x_0[i] * fe_0 -
                               ta1_x_xxyyy_x_1[i] * fe_0 + ta1_x_xxyyy_xy_0[i] * pa_y[i] - ta1_x_xxyyy_xy_1[i] * pc_y[i];

        ta1_x_xxyyyy_xz_0[i] =
            3.0 * ta1_x_xxyy_xz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_xz_1[i] * fe_0 + ta1_x_xxyyy_xz_0[i] * pa_y[i] - ta1_x_xxyyy_xz_1[i] * pc_y[i];

        ta1_x_xxyyyy_yy_0[i] =
            ta1_x_yyyy_yy_0[i] * fe_0 - ta1_x_yyyy_yy_1[i] * fe_0 + ta_xyyyy_yy_1[i] + ta1_x_xyyyy_yy_0[i] * pa_x[i] - ta1_x_xyyyy_yy_1[i] * pc_x[i];

        ta1_x_xxyyyy_yz_0[i] =
            ta1_x_yyyy_yz_0[i] * fe_0 - ta1_x_yyyy_yz_1[i] * fe_0 + ta_xyyyy_yz_1[i] + ta1_x_xyyyy_yz_0[i] * pa_x[i] - ta1_x_xyyyy_yz_1[i] * pc_x[i];

        ta1_x_xxyyyy_zz_0[i] =
            3.0 * ta1_x_xxyy_zz_0[i] * fe_0 - 3.0 * ta1_x_xxyy_zz_1[i] * fe_0 + ta1_x_xxyyy_zz_0[i] * pa_y[i] - ta1_x_xxyyy_zz_1[i] * pc_y[i];
    }

    // Set up 66-72 components of targeted buffer : ID

    auto ta1_x_xxyyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 66);

    auto ta1_x_xxyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 67);

    auto ta1_x_xxyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 68);

    auto ta1_x_xxyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 69);

    auto ta1_x_xxyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 70);

    auto ta1_x_xxyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 71);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_xxyyy_xx_0,  \
                             ta1_x_xxyyy_xx_1,  \
                             ta1_x_xxyyy_xy_0,  \
                             ta1_x_xxyyy_xy_1,  \
                             ta1_x_xxyyy_y_0,   \
                             ta1_x_xxyyy_y_1,   \
                             ta1_x_xxyyy_yy_0,  \
                             ta1_x_xxyyy_yy_1,  \
                             ta1_x_xxyyy_yz_0,  \
                             ta1_x_xxyyy_yz_1,  \
                             ta1_x_xxyyyz_xx_0, \
                             ta1_x_xxyyyz_xy_0, \
                             ta1_x_xxyyyz_xz_0, \
                             ta1_x_xxyyyz_yy_0, \
                             ta1_x_xxyyyz_yz_0, \
                             ta1_x_xxyyyz_zz_0, \
                             ta1_x_xxyyz_xz_0,  \
                             ta1_x_xxyyz_xz_1,  \
                             ta1_x_xxyyz_zz_0,  \
                             ta1_x_xxyyz_zz_1,  \
                             ta1_x_xxyz_xz_0,   \
                             ta1_x_xxyz_xz_1,   \
                             ta1_x_xxyz_zz_0,   \
                             ta1_x_xxyz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyyz_xx_0[i] = ta1_x_xxyyy_xx_0[i] * pa_z[i] - ta1_x_xxyyy_xx_1[i] * pc_z[i];

        ta1_x_xxyyyz_xy_0[i] = ta1_x_xxyyy_xy_0[i] * pa_z[i] - ta1_x_xxyyy_xy_1[i] * pc_z[i];

        ta1_x_xxyyyz_xz_0[i] =
            2.0 * ta1_x_xxyz_xz_0[i] * fe_0 - 2.0 * ta1_x_xxyz_xz_1[i] * fe_0 + ta1_x_xxyyz_xz_0[i] * pa_y[i] - ta1_x_xxyyz_xz_1[i] * pc_y[i];

        ta1_x_xxyyyz_yy_0[i] = ta1_x_xxyyy_yy_0[i] * pa_z[i] - ta1_x_xxyyy_yy_1[i] * pc_z[i];

        ta1_x_xxyyyz_yz_0[i] = ta1_x_xxyyy_y_0[i] * fe_0 - ta1_x_xxyyy_y_1[i] * fe_0 + ta1_x_xxyyy_yz_0[i] * pa_z[i] - ta1_x_xxyyy_yz_1[i] * pc_z[i];

        ta1_x_xxyyyz_zz_0[i] =
            2.0 * ta1_x_xxyz_zz_0[i] * fe_0 - 2.0 * ta1_x_xxyz_zz_1[i] * fe_0 + ta1_x_xxyyz_zz_0[i] * pa_y[i] - ta1_x_xxyyz_zz_1[i] * pc_y[i];
    }

    // Set up 72-78 components of targeted buffer : ID

    auto ta1_x_xxyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 72);

    auto ta1_x_xxyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 73);

    auto ta1_x_xxyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 74);

    auto ta1_x_xxyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 75);

    auto ta1_x_xxyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 76);

    auto ta1_x_xxyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 77);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_xxyy_xy_0,   \
                             ta1_x_xxyy_xy_1,   \
                             ta1_x_xxyy_yy_0,   \
                             ta1_x_xxyy_yy_1,   \
                             ta1_x_xxyyz_xy_0,  \
                             ta1_x_xxyyz_xy_1,  \
                             ta1_x_xxyyz_yy_0,  \
                             ta1_x_xxyyz_yy_1,  \
                             ta1_x_xxyyzz_xx_0, \
                             ta1_x_xxyyzz_xy_0, \
                             ta1_x_xxyyzz_xz_0, \
                             ta1_x_xxyyzz_yy_0, \
                             ta1_x_xxyyzz_yz_0, \
                             ta1_x_xxyyzz_zz_0, \
                             ta1_x_xxyzz_xx_0,  \
                             ta1_x_xxyzz_xx_1,  \
                             ta1_x_xxyzz_xz_0,  \
                             ta1_x_xxyzz_xz_1,  \
                             ta1_x_xxyzz_zz_0,  \
                             ta1_x_xxyzz_zz_1,  \
                             ta1_x_xxzz_xx_0,   \
                             ta1_x_xxzz_xx_1,   \
                             ta1_x_xxzz_xz_0,   \
                             ta1_x_xxzz_xz_1,   \
                             ta1_x_xxzz_zz_0,   \
                             ta1_x_xxzz_zz_1,   \
                             ta1_x_xyyzz_yz_0,  \
                             ta1_x_xyyzz_yz_1,  \
                             ta1_x_yyzz_yz_0,   \
                             ta1_x_yyzz_yz_1,   \
                             ta_xyyzz_yz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyyzz_xx_0[i] = ta1_x_xxzz_xx_0[i] * fe_0 - ta1_x_xxzz_xx_1[i] * fe_0 + ta1_x_xxyzz_xx_0[i] * pa_y[i] - ta1_x_xxyzz_xx_1[i] * pc_y[i];

        ta1_x_xxyyzz_xy_0[i] = ta1_x_xxyy_xy_0[i] * fe_0 - ta1_x_xxyy_xy_1[i] * fe_0 + ta1_x_xxyyz_xy_0[i] * pa_z[i] - ta1_x_xxyyz_xy_1[i] * pc_z[i];

        ta1_x_xxyyzz_xz_0[i] = ta1_x_xxzz_xz_0[i] * fe_0 - ta1_x_xxzz_xz_1[i] * fe_0 + ta1_x_xxyzz_xz_0[i] * pa_y[i] - ta1_x_xxyzz_xz_1[i] * pc_y[i];

        ta1_x_xxyyzz_yy_0[i] = ta1_x_xxyy_yy_0[i] * fe_0 - ta1_x_xxyy_yy_1[i] * fe_0 + ta1_x_xxyyz_yy_0[i] * pa_z[i] - ta1_x_xxyyz_yy_1[i] * pc_z[i];

        ta1_x_xxyyzz_yz_0[i] =
            ta1_x_yyzz_yz_0[i] * fe_0 - ta1_x_yyzz_yz_1[i] * fe_0 + ta_xyyzz_yz_1[i] + ta1_x_xyyzz_yz_0[i] * pa_x[i] - ta1_x_xyyzz_yz_1[i] * pc_x[i];

        ta1_x_xxyyzz_zz_0[i] = ta1_x_xxzz_zz_0[i] * fe_0 - ta1_x_xxzz_zz_1[i] * fe_0 + ta1_x_xxyzz_zz_0[i] * pa_y[i] - ta1_x_xxyzz_zz_1[i] * pc_y[i];
    }

    // Set up 78-84 components of targeted buffer : ID

    auto ta1_x_xxyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 78);

    auto ta1_x_xxyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 79);

    auto ta1_x_xxyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 80);

    auto ta1_x_xxyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 81);

    auto ta1_x_xxyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 82);

    auto ta1_x_xxyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 83);

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_x_xxyzzz_xx_0, \
                             ta1_x_xxyzzz_xy_0, \
                             ta1_x_xxyzzz_xz_0, \
                             ta1_x_xxyzzz_yy_0, \
                             ta1_x_xxyzzz_yz_0, \
                             ta1_x_xxyzzz_zz_0, \
                             ta1_x_xxzzz_x_0,   \
                             ta1_x_xxzzz_x_1,   \
                             ta1_x_xxzzz_xx_0,  \
                             ta1_x_xxzzz_xx_1,  \
                             ta1_x_xxzzz_xy_0,  \
                             ta1_x_xxzzz_xy_1,  \
                             ta1_x_xxzzz_xz_0,  \
                             ta1_x_xxzzz_xz_1,  \
                             ta1_x_xxzzz_y_0,   \
                             ta1_x_xxzzz_y_1,   \
                             ta1_x_xxzzz_yy_0,  \
                             ta1_x_xxzzz_yy_1,  \
                             ta1_x_xxzzz_yz_0,  \
                             ta1_x_xxzzz_yz_1,  \
                             ta1_x_xxzzz_z_0,   \
                             ta1_x_xxzzz_z_1,   \
                             ta1_x_xxzzz_zz_0,  \
                             ta1_x_xxzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyzzz_xx_0[i] = ta1_x_xxzzz_xx_0[i] * pa_y[i] - ta1_x_xxzzz_xx_1[i] * pc_y[i];

        ta1_x_xxyzzz_xy_0[i] = ta1_x_xxzzz_x_0[i] * fe_0 - ta1_x_xxzzz_x_1[i] * fe_0 + ta1_x_xxzzz_xy_0[i] * pa_y[i] - ta1_x_xxzzz_xy_1[i] * pc_y[i];

        ta1_x_xxyzzz_xz_0[i] = ta1_x_xxzzz_xz_0[i] * pa_y[i] - ta1_x_xxzzz_xz_1[i] * pc_y[i];

        ta1_x_xxyzzz_yy_0[i] =
            2.0 * ta1_x_xxzzz_y_0[i] * fe_0 - 2.0 * ta1_x_xxzzz_y_1[i] * fe_0 + ta1_x_xxzzz_yy_0[i] * pa_y[i] - ta1_x_xxzzz_yy_1[i] * pc_y[i];

        ta1_x_xxyzzz_yz_0[i] = ta1_x_xxzzz_z_0[i] * fe_0 - ta1_x_xxzzz_z_1[i] * fe_0 + ta1_x_xxzzz_yz_0[i] * pa_y[i] - ta1_x_xxzzz_yz_1[i] * pc_y[i];

        ta1_x_xxyzzz_zz_0[i] = ta1_x_xxzzz_zz_0[i] * pa_y[i] - ta1_x_xxzzz_zz_1[i] * pc_y[i];
    }

    // Set up 84-90 components of targeted buffer : ID

    auto ta1_x_xxzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 84);

    auto ta1_x_xxzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 85);

    auto ta1_x_xxzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 86);

    auto ta1_x_xxzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 87);

    auto ta1_x_xxzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 88);

    auto ta1_x_xxzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 89);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_x_xxzz_xx_0,   \
                             ta1_x_xxzz_xx_1,   \
                             ta1_x_xxzz_xy_0,   \
                             ta1_x_xxzz_xy_1,   \
                             ta1_x_xxzz_xz_0,   \
                             ta1_x_xxzz_xz_1,   \
                             ta1_x_xxzz_yy_0,   \
                             ta1_x_xxzz_yy_1,   \
                             ta1_x_xxzzz_x_0,   \
                             ta1_x_xxzzz_x_1,   \
                             ta1_x_xxzzz_xx_0,  \
                             ta1_x_xxzzz_xx_1,  \
                             ta1_x_xxzzz_xy_0,  \
                             ta1_x_xxzzz_xy_1,  \
                             ta1_x_xxzzz_xz_0,  \
                             ta1_x_xxzzz_xz_1,  \
                             ta1_x_xxzzz_yy_0,  \
                             ta1_x_xxzzz_yy_1,  \
                             ta1_x_xxzzzz_xx_0, \
                             ta1_x_xxzzzz_xy_0, \
                             ta1_x_xxzzzz_xz_0, \
                             ta1_x_xxzzzz_yy_0, \
                             ta1_x_xxzzzz_yz_0, \
                             ta1_x_xxzzzz_zz_0, \
                             ta1_x_xzzzz_yz_0,  \
                             ta1_x_xzzzz_yz_1,  \
                             ta1_x_xzzzz_zz_0,  \
                             ta1_x_xzzzz_zz_1,  \
                             ta1_x_zzzz_yz_0,   \
                             ta1_x_zzzz_yz_1,   \
                             ta1_x_zzzz_zz_0,   \
                             ta1_x_zzzz_zz_1,   \
                             ta_xzzzz_yz_1,     \
                             ta_xzzzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxzzzz_xx_0[i] =
            3.0 * ta1_x_xxzz_xx_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xx_1[i] * fe_0 + ta1_x_xxzzz_xx_0[i] * pa_z[i] - ta1_x_xxzzz_xx_1[i] * pc_z[i];

        ta1_x_xxzzzz_xy_0[i] =
            3.0 * ta1_x_xxzz_xy_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xy_1[i] * fe_0 + ta1_x_xxzzz_xy_0[i] * pa_z[i] - ta1_x_xxzzz_xy_1[i] * pc_z[i];

        ta1_x_xxzzzz_xz_0[i] = 3.0 * ta1_x_xxzz_xz_0[i] * fe_0 - 3.0 * ta1_x_xxzz_xz_1[i] * fe_0 + ta1_x_xxzzz_x_0[i] * fe_0 -
                               ta1_x_xxzzz_x_1[i] * fe_0 + ta1_x_xxzzz_xz_0[i] * pa_z[i] - ta1_x_xxzzz_xz_1[i] * pc_z[i];

        ta1_x_xxzzzz_yy_0[i] =
            3.0 * ta1_x_xxzz_yy_0[i] * fe_0 - 3.0 * ta1_x_xxzz_yy_1[i] * fe_0 + ta1_x_xxzzz_yy_0[i] * pa_z[i] - ta1_x_xxzzz_yy_1[i] * pc_z[i];

        ta1_x_xxzzzz_yz_0[i] =
            ta1_x_zzzz_yz_0[i] * fe_0 - ta1_x_zzzz_yz_1[i] * fe_0 + ta_xzzzz_yz_1[i] + ta1_x_xzzzz_yz_0[i] * pa_x[i] - ta1_x_xzzzz_yz_1[i] * pc_x[i];

        ta1_x_xxzzzz_zz_0[i] =
            ta1_x_zzzz_zz_0[i] * fe_0 - ta1_x_zzzz_zz_1[i] * fe_0 + ta_xzzzz_zz_1[i] + ta1_x_xzzzz_zz_0[i] * pa_x[i] - ta1_x_xzzzz_zz_1[i] * pc_x[i];
    }

    // Set up 90-96 components of targeted buffer : ID

    auto ta1_x_xyyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 90);

    auto ta1_x_xyyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 91);

    auto ta1_x_xyyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 92);

    auto ta1_x_xyyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 93);

    auto ta1_x_xyyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 94);

    auto ta1_x_xyyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 95);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_x_xyyy_xx_0,   \
                             ta1_x_xyyy_xx_1,   \
                             ta1_x_xyyy_xz_0,   \
                             ta1_x_xyyy_xz_1,   \
                             ta1_x_xyyyy_xx_0,  \
                             ta1_x_xyyyy_xx_1,  \
                             ta1_x_xyyyy_xz_0,  \
                             ta1_x_xyyyy_xz_1,  \
                             ta1_x_xyyyyy_xx_0, \
                             ta1_x_xyyyyy_xy_0, \
                             ta1_x_xyyyyy_xz_0, \
                             ta1_x_xyyyyy_yy_0, \
                             ta1_x_xyyyyy_yz_0, \
                             ta1_x_xyyyyy_zz_0, \
                             ta1_x_yyyyy_xy_0,  \
                             ta1_x_yyyyy_xy_1,  \
                             ta1_x_yyyyy_y_0,   \
                             ta1_x_yyyyy_y_1,   \
                             ta1_x_yyyyy_yy_0,  \
                             ta1_x_yyyyy_yy_1,  \
                             ta1_x_yyyyy_yz_0,  \
                             ta1_x_yyyyy_yz_1,  \
                             ta1_x_yyyyy_zz_0,  \
                             ta1_x_yyyyy_zz_1,  \
                             ta_yyyyy_xy_1,     \
                             ta_yyyyy_yy_1,     \
                             ta_yyyyy_yz_1,     \
                             ta_yyyyy_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyyy_xx_0[i] =
            4.0 * ta1_x_xyyy_xx_0[i] * fe_0 - 4.0 * ta1_x_xyyy_xx_1[i] * fe_0 + ta1_x_xyyyy_xx_0[i] * pa_y[i] - ta1_x_xyyyy_xx_1[i] * pc_y[i];

        ta1_x_xyyyyy_xy_0[i] =
            ta1_x_yyyyy_y_0[i] * fe_0 - ta1_x_yyyyy_y_1[i] * fe_0 + ta_yyyyy_xy_1[i] + ta1_x_yyyyy_xy_0[i] * pa_x[i] - ta1_x_yyyyy_xy_1[i] * pc_x[i];

        ta1_x_xyyyyy_xz_0[i] =
            4.0 * ta1_x_xyyy_xz_0[i] * fe_0 - 4.0 * ta1_x_xyyy_xz_1[i] * fe_0 + ta1_x_xyyyy_xz_0[i] * pa_y[i] - ta1_x_xyyyy_xz_1[i] * pc_y[i];

        ta1_x_xyyyyy_yy_0[i] = ta_yyyyy_yy_1[i] + ta1_x_yyyyy_yy_0[i] * pa_x[i] - ta1_x_yyyyy_yy_1[i] * pc_x[i];

        ta1_x_xyyyyy_yz_0[i] = ta_yyyyy_yz_1[i] + ta1_x_yyyyy_yz_0[i] * pa_x[i] - ta1_x_yyyyy_yz_1[i] * pc_x[i];

        ta1_x_xyyyyy_zz_0[i] = ta_yyyyy_zz_1[i] + ta1_x_yyyyy_zz_0[i] * pa_x[i] - ta1_x_yyyyy_zz_1[i] * pc_x[i];
    }

    // Set up 96-102 components of targeted buffer : ID

    auto ta1_x_xyyyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 96);

    auto ta1_x_xyyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 97);

    auto ta1_x_xyyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 98);

    auto ta1_x_xyyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 99);

    auto ta1_x_xyyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 100);

    auto ta1_x_xyyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 101);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_xyyyy_xx_0,  \
                             ta1_x_xyyyy_xx_1,  \
                             ta1_x_xyyyy_xy_0,  \
                             ta1_x_xyyyy_xy_1,  \
                             ta1_x_xyyyy_yy_0,  \
                             ta1_x_xyyyy_yy_1,  \
                             ta1_x_xyyyyz_xx_0, \
                             ta1_x_xyyyyz_xy_0, \
                             ta1_x_xyyyyz_xz_0, \
                             ta1_x_xyyyyz_yy_0, \
                             ta1_x_xyyyyz_yz_0, \
                             ta1_x_xyyyyz_zz_0, \
                             ta1_x_xyyyz_xz_0,  \
                             ta1_x_xyyyz_xz_1,  \
                             ta1_x_xyyz_xz_0,   \
                             ta1_x_xyyz_xz_1,   \
                             ta1_x_yyyyz_yz_0,  \
                             ta1_x_yyyyz_yz_1,  \
                             ta1_x_yyyyz_zz_0,  \
                             ta1_x_yyyyz_zz_1,  \
                             ta_yyyyz_yz_1,     \
                             ta_yyyyz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyyz_xx_0[i] = ta1_x_xyyyy_xx_0[i] * pa_z[i] - ta1_x_xyyyy_xx_1[i] * pc_z[i];

        ta1_x_xyyyyz_xy_0[i] = ta1_x_xyyyy_xy_0[i] * pa_z[i] - ta1_x_xyyyy_xy_1[i] * pc_z[i];

        ta1_x_xyyyyz_xz_0[i] =
            3.0 * ta1_x_xyyz_xz_0[i] * fe_0 - 3.0 * ta1_x_xyyz_xz_1[i] * fe_0 + ta1_x_xyyyz_xz_0[i] * pa_y[i] - ta1_x_xyyyz_xz_1[i] * pc_y[i];

        ta1_x_xyyyyz_yy_0[i] = ta1_x_xyyyy_yy_0[i] * pa_z[i] - ta1_x_xyyyy_yy_1[i] * pc_z[i];

        ta1_x_xyyyyz_yz_0[i] = ta_yyyyz_yz_1[i] + ta1_x_yyyyz_yz_0[i] * pa_x[i] - ta1_x_yyyyz_yz_1[i] * pc_x[i];

        ta1_x_xyyyyz_zz_0[i] = ta_yyyyz_zz_1[i] + ta1_x_yyyyz_zz_0[i] * pa_x[i] - ta1_x_yyyyz_zz_1[i] * pc_x[i];
    }

    // Set up 102-108 components of targeted buffer : ID

    auto ta1_x_xyyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 102);

    auto ta1_x_xyyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 103);

    auto ta1_x_xyyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 104);

    auto ta1_x_xyyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 105);

    auto ta1_x_xyyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 106);

    auto ta1_x_xyyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 107);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_xyyy_xy_0,   \
                             ta1_x_xyyy_xy_1,   \
                             ta1_x_xyyyz_xy_0,  \
                             ta1_x_xyyyz_xy_1,  \
                             ta1_x_xyyyzz_xx_0, \
                             ta1_x_xyyyzz_xy_0, \
                             ta1_x_xyyyzz_xz_0, \
                             ta1_x_xyyyzz_yy_0, \
                             ta1_x_xyyyzz_yz_0, \
                             ta1_x_xyyyzz_zz_0, \
                             ta1_x_xyyzz_xx_0,  \
                             ta1_x_xyyzz_xx_1,  \
                             ta1_x_xyyzz_xz_0,  \
                             ta1_x_xyyzz_xz_1,  \
                             ta1_x_xyzz_xx_0,   \
                             ta1_x_xyzz_xx_1,   \
                             ta1_x_xyzz_xz_0,   \
                             ta1_x_xyzz_xz_1,   \
                             ta1_x_yyyzz_yy_0,  \
                             ta1_x_yyyzz_yy_1,  \
                             ta1_x_yyyzz_yz_0,  \
                             ta1_x_yyyzz_yz_1,  \
                             ta1_x_yyyzz_zz_0,  \
                             ta1_x_yyyzz_zz_1,  \
                             ta_yyyzz_yy_1,     \
                             ta_yyyzz_yz_1,     \
                             ta_yyyzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyyzz_xx_0[i] =
            2.0 * ta1_x_xyzz_xx_0[i] * fe_0 - 2.0 * ta1_x_xyzz_xx_1[i] * fe_0 + ta1_x_xyyzz_xx_0[i] * pa_y[i] - ta1_x_xyyzz_xx_1[i] * pc_y[i];

        ta1_x_xyyyzz_xy_0[i] = ta1_x_xyyy_xy_0[i] * fe_0 - ta1_x_xyyy_xy_1[i] * fe_0 + ta1_x_xyyyz_xy_0[i] * pa_z[i] - ta1_x_xyyyz_xy_1[i] * pc_z[i];

        ta1_x_xyyyzz_xz_0[i] =
            2.0 * ta1_x_xyzz_xz_0[i] * fe_0 - 2.0 * ta1_x_xyzz_xz_1[i] * fe_0 + ta1_x_xyyzz_xz_0[i] * pa_y[i] - ta1_x_xyyzz_xz_1[i] * pc_y[i];

        ta1_x_xyyyzz_yy_0[i] = ta_yyyzz_yy_1[i] + ta1_x_yyyzz_yy_0[i] * pa_x[i] - ta1_x_yyyzz_yy_1[i] * pc_x[i];

        ta1_x_xyyyzz_yz_0[i] = ta_yyyzz_yz_1[i] + ta1_x_yyyzz_yz_0[i] * pa_x[i] - ta1_x_yyyzz_yz_1[i] * pc_x[i];

        ta1_x_xyyyzz_zz_0[i] = ta_yyyzz_zz_1[i] + ta1_x_yyyzz_zz_0[i] * pa_x[i] - ta1_x_yyyzz_zz_1[i] * pc_x[i];
    }

    // Set up 108-114 components of targeted buffer : ID

    auto ta1_x_xyyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 108);

    auto ta1_x_xyyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 109);

    auto ta1_x_xyyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 110);

    auto ta1_x_xyyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 111);

    auto ta1_x_xyyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 112);

    auto ta1_x_xyyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 113);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_xyyz_xy_0,   \
                             ta1_x_xyyz_xy_1,   \
                             ta1_x_xyyzz_xy_0,  \
                             ta1_x_xyyzz_xy_1,  \
                             ta1_x_xyyzzz_xx_0, \
                             ta1_x_xyyzzz_xy_0, \
                             ta1_x_xyyzzz_xz_0, \
                             ta1_x_xyyzzz_yy_0, \
                             ta1_x_xyyzzz_yz_0, \
                             ta1_x_xyyzzz_zz_0, \
                             ta1_x_xyzzz_xx_0,  \
                             ta1_x_xyzzz_xx_1,  \
                             ta1_x_xyzzz_xz_0,  \
                             ta1_x_xyzzz_xz_1,  \
                             ta1_x_xzzz_xx_0,   \
                             ta1_x_xzzz_xx_1,   \
                             ta1_x_xzzz_xz_0,   \
                             ta1_x_xzzz_xz_1,   \
                             ta1_x_yyzzz_yy_0,  \
                             ta1_x_yyzzz_yy_1,  \
                             ta1_x_yyzzz_yz_0,  \
                             ta1_x_yyzzz_yz_1,  \
                             ta1_x_yyzzz_zz_0,  \
                             ta1_x_yyzzz_zz_1,  \
                             ta_yyzzz_yy_1,     \
                             ta_yyzzz_yz_1,     \
                             ta_yyzzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyzzz_xx_0[i] = ta1_x_xzzz_xx_0[i] * fe_0 - ta1_x_xzzz_xx_1[i] * fe_0 + ta1_x_xyzzz_xx_0[i] * pa_y[i] - ta1_x_xyzzz_xx_1[i] * pc_y[i];

        ta1_x_xyyzzz_xy_0[i] =
            2.0 * ta1_x_xyyz_xy_0[i] * fe_0 - 2.0 * ta1_x_xyyz_xy_1[i] * fe_0 + ta1_x_xyyzz_xy_0[i] * pa_z[i] - ta1_x_xyyzz_xy_1[i] * pc_z[i];

        ta1_x_xyyzzz_xz_0[i] = ta1_x_xzzz_xz_0[i] * fe_0 - ta1_x_xzzz_xz_1[i] * fe_0 + ta1_x_xyzzz_xz_0[i] * pa_y[i] - ta1_x_xyzzz_xz_1[i] * pc_y[i];

        ta1_x_xyyzzz_yy_0[i] = ta_yyzzz_yy_1[i] + ta1_x_yyzzz_yy_0[i] * pa_x[i] - ta1_x_yyzzz_yy_1[i] * pc_x[i];

        ta1_x_xyyzzz_yz_0[i] = ta_yyzzz_yz_1[i] + ta1_x_yyzzz_yz_0[i] * pa_x[i] - ta1_x_yyzzz_yz_1[i] * pc_x[i];

        ta1_x_xyyzzz_zz_0[i] = ta_yyzzz_zz_1[i] + ta1_x_yyzzz_zz_0[i] * pa_x[i] - ta1_x_yyzzz_zz_1[i] * pc_x[i];
    }

    // Set up 114-120 components of targeted buffer : ID

    auto ta1_x_xyzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 114);

    auto ta1_x_xyzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 115);

    auto ta1_x_xyzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 116);

    auto ta1_x_xyzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 117);

    auto ta1_x_xyzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 118);

    auto ta1_x_xyzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 119);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_x_xyzzzz_xx_0, \
                             ta1_x_xyzzzz_xy_0, \
                             ta1_x_xyzzzz_xz_0, \
                             ta1_x_xyzzzz_yy_0, \
                             ta1_x_xyzzzz_yz_0, \
                             ta1_x_xyzzzz_zz_0, \
                             ta1_x_xzzzz_x_0,   \
                             ta1_x_xzzzz_x_1,   \
                             ta1_x_xzzzz_xx_0,  \
                             ta1_x_xzzzz_xx_1,  \
                             ta1_x_xzzzz_xy_0,  \
                             ta1_x_xzzzz_xy_1,  \
                             ta1_x_xzzzz_xz_0,  \
                             ta1_x_xzzzz_xz_1,  \
                             ta1_x_xzzzz_zz_0,  \
                             ta1_x_xzzzz_zz_1,  \
                             ta1_x_yzzzz_yy_0,  \
                             ta1_x_yzzzz_yy_1,  \
                             ta1_x_yzzzz_yz_0,  \
                             ta1_x_yzzzz_yz_1,  \
                             ta_yzzzz_yy_1,     \
                             ta_yzzzz_yz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyzzzz_xx_0[i] = ta1_x_xzzzz_xx_0[i] * pa_y[i] - ta1_x_xzzzz_xx_1[i] * pc_y[i];

        ta1_x_xyzzzz_xy_0[i] = ta1_x_xzzzz_x_0[i] * fe_0 - ta1_x_xzzzz_x_1[i] * fe_0 + ta1_x_xzzzz_xy_0[i] * pa_y[i] - ta1_x_xzzzz_xy_1[i] * pc_y[i];

        ta1_x_xyzzzz_xz_0[i] = ta1_x_xzzzz_xz_0[i] * pa_y[i] - ta1_x_xzzzz_xz_1[i] * pc_y[i];

        ta1_x_xyzzzz_yy_0[i] = ta_yzzzz_yy_1[i] + ta1_x_yzzzz_yy_0[i] * pa_x[i] - ta1_x_yzzzz_yy_1[i] * pc_x[i];

        ta1_x_xyzzzz_yz_0[i] = ta_yzzzz_yz_1[i] + ta1_x_yzzzz_yz_0[i] * pa_x[i] - ta1_x_yzzzz_yz_1[i] * pc_x[i];

        ta1_x_xyzzzz_zz_0[i] = ta1_x_xzzzz_zz_0[i] * pa_y[i] - ta1_x_xzzzz_zz_1[i] * pc_y[i];
    }

    // Set up 120-126 components of targeted buffer : ID

    auto ta1_x_xzzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 120);

    auto ta1_x_xzzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 121);

    auto ta1_x_xzzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 122);

    auto ta1_x_xzzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 123);

    auto ta1_x_xzzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 124);

    auto ta1_x_xzzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 125);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_x_xzzz_xx_0,   \
                             ta1_x_xzzz_xx_1,   \
                             ta1_x_xzzz_xy_0,   \
                             ta1_x_xzzz_xy_1,   \
                             ta1_x_xzzzz_xx_0,  \
                             ta1_x_xzzzz_xx_1,  \
                             ta1_x_xzzzz_xy_0,  \
                             ta1_x_xzzzz_xy_1,  \
                             ta1_x_xzzzzz_xx_0, \
                             ta1_x_xzzzzz_xy_0, \
                             ta1_x_xzzzzz_xz_0, \
                             ta1_x_xzzzzz_yy_0, \
                             ta1_x_xzzzzz_yz_0, \
                             ta1_x_xzzzzz_zz_0, \
                             ta1_x_zzzzz_xz_0,  \
                             ta1_x_zzzzz_xz_1,  \
                             ta1_x_zzzzz_yy_0,  \
                             ta1_x_zzzzz_yy_1,  \
                             ta1_x_zzzzz_yz_0,  \
                             ta1_x_zzzzz_yz_1,  \
                             ta1_x_zzzzz_z_0,   \
                             ta1_x_zzzzz_z_1,   \
                             ta1_x_zzzzz_zz_0,  \
                             ta1_x_zzzzz_zz_1,  \
                             ta_zzzzz_xz_1,     \
                             ta_zzzzz_yy_1,     \
                             ta_zzzzz_yz_1,     \
                             ta_zzzzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzzzzz_xx_0[i] =
            4.0 * ta1_x_xzzz_xx_0[i] * fe_0 - 4.0 * ta1_x_xzzz_xx_1[i] * fe_0 + ta1_x_xzzzz_xx_0[i] * pa_z[i] - ta1_x_xzzzz_xx_1[i] * pc_z[i];

        ta1_x_xzzzzz_xy_0[i] =
            4.0 * ta1_x_xzzz_xy_0[i] * fe_0 - 4.0 * ta1_x_xzzz_xy_1[i] * fe_0 + ta1_x_xzzzz_xy_0[i] * pa_z[i] - ta1_x_xzzzz_xy_1[i] * pc_z[i];

        ta1_x_xzzzzz_xz_0[i] =
            ta1_x_zzzzz_z_0[i] * fe_0 - ta1_x_zzzzz_z_1[i] * fe_0 + ta_zzzzz_xz_1[i] + ta1_x_zzzzz_xz_0[i] * pa_x[i] - ta1_x_zzzzz_xz_1[i] * pc_x[i];

        ta1_x_xzzzzz_yy_0[i] = ta_zzzzz_yy_1[i] + ta1_x_zzzzz_yy_0[i] * pa_x[i] - ta1_x_zzzzz_yy_1[i] * pc_x[i];

        ta1_x_xzzzzz_yz_0[i] = ta_zzzzz_yz_1[i] + ta1_x_zzzzz_yz_0[i] * pa_x[i] - ta1_x_zzzzz_yz_1[i] * pc_x[i];

        ta1_x_xzzzzz_zz_0[i] = ta_zzzzz_zz_1[i] + ta1_x_zzzzz_zz_0[i] * pa_x[i] - ta1_x_zzzzz_zz_1[i] * pc_x[i];
    }

    // Set up 126-132 components of targeted buffer : ID

    auto ta1_x_yyyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 126);

    auto ta1_x_yyyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 127);

    auto ta1_x_yyyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 128);

    auto ta1_x_yyyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 129);

    auto ta1_x_yyyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 130);

    auto ta1_x_yyyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 131);

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_x_yyyy_xx_0,   \
                             ta1_x_yyyy_xx_1,   \
                             ta1_x_yyyy_xy_0,   \
                             ta1_x_yyyy_xy_1,   \
                             ta1_x_yyyy_xz_0,   \
                             ta1_x_yyyy_xz_1,   \
                             ta1_x_yyyy_yy_0,   \
                             ta1_x_yyyy_yy_1,   \
                             ta1_x_yyyy_yz_0,   \
                             ta1_x_yyyy_yz_1,   \
                             ta1_x_yyyy_zz_0,   \
                             ta1_x_yyyy_zz_1,   \
                             ta1_x_yyyyy_x_0,   \
                             ta1_x_yyyyy_x_1,   \
                             ta1_x_yyyyy_xx_0,  \
                             ta1_x_yyyyy_xx_1,  \
                             ta1_x_yyyyy_xy_0,  \
                             ta1_x_yyyyy_xy_1,  \
                             ta1_x_yyyyy_xz_0,  \
                             ta1_x_yyyyy_xz_1,  \
                             ta1_x_yyyyy_y_0,   \
                             ta1_x_yyyyy_y_1,   \
                             ta1_x_yyyyy_yy_0,  \
                             ta1_x_yyyyy_yy_1,  \
                             ta1_x_yyyyy_yz_0,  \
                             ta1_x_yyyyy_yz_1,  \
                             ta1_x_yyyyy_z_0,   \
                             ta1_x_yyyyy_z_1,   \
                             ta1_x_yyyyy_zz_0,  \
                             ta1_x_yyyyy_zz_1,  \
                             ta1_x_yyyyyy_xx_0, \
                             ta1_x_yyyyyy_xy_0, \
                             ta1_x_yyyyyy_xz_0, \
                             ta1_x_yyyyyy_yy_0, \
                             ta1_x_yyyyyy_yz_0, \
                             ta1_x_yyyyyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyyy_xx_0[i] =
            5.0 * ta1_x_yyyy_xx_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xx_1[i] * fe_0 + ta1_x_yyyyy_xx_0[i] * pa_y[i] - ta1_x_yyyyy_xx_1[i] * pc_y[i];

        ta1_x_yyyyyy_xy_0[i] = 5.0 * ta1_x_yyyy_xy_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xy_1[i] * fe_0 + ta1_x_yyyyy_x_0[i] * fe_0 -
                               ta1_x_yyyyy_x_1[i] * fe_0 + ta1_x_yyyyy_xy_0[i] * pa_y[i] - ta1_x_yyyyy_xy_1[i] * pc_y[i];

        ta1_x_yyyyyy_xz_0[i] =
            5.0 * ta1_x_yyyy_xz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_xz_1[i] * fe_0 + ta1_x_yyyyy_xz_0[i] * pa_y[i] - ta1_x_yyyyy_xz_1[i] * pc_y[i];

        ta1_x_yyyyyy_yy_0[i] = 5.0 * ta1_x_yyyy_yy_0[i] * fe_0 - 5.0 * ta1_x_yyyy_yy_1[i] * fe_0 + 2.0 * ta1_x_yyyyy_y_0[i] * fe_0 -
                               2.0 * ta1_x_yyyyy_y_1[i] * fe_0 + ta1_x_yyyyy_yy_0[i] * pa_y[i] - ta1_x_yyyyy_yy_1[i] * pc_y[i];

        ta1_x_yyyyyy_yz_0[i] = 5.0 * ta1_x_yyyy_yz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_yz_1[i] * fe_0 + ta1_x_yyyyy_z_0[i] * fe_0 -
                               ta1_x_yyyyy_z_1[i] * fe_0 + ta1_x_yyyyy_yz_0[i] * pa_y[i] - ta1_x_yyyyy_yz_1[i] * pc_y[i];

        ta1_x_yyyyyy_zz_0[i] =
            5.0 * ta1_x_yyyy_zz_0[i] * fe_0 - 5.0 * ta1_x_yyyy_zz_1[i] * fe_0 + ta1_x_yyyyy_zz_0[i] * pa_y[i] - ta1_x_yyyyy_zz_1[i] * pc_y[i];
    }

    // Set up 132-138 components of targeted buffer : ID

    auto ta1_x_yyyyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 132);

    auto ta1_x_yyyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 133);

    auto ta1_x_yyyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 134);

    auto ta1_x_yyyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 135);

    auto ta1_x_yyyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 136);

    auto ta1_x_yyyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 137);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_yyyyy_xx_0,  \
                             ta1_x_yyyyy_xx_1,  \
                             ta1_x_yyyyy_xy_0,  \
                             ta1_x_yyyyy_xy_1,  \
                             ta1_x_yyyyy_y_0,   \
                             ta1_x_yyyyy_y_1,   \
                             ta1_x_yyyyy_yy_0,  \
                             ta1_x_yyyyy_yy_1,  \
                             ta1_x_yyyyy_yz_0,  \
                             ta1_x_yyyyy_yz_1,  \
                             ta1_x_yyyyyz_xx_0, \
                             ta1_x_yyyyyz_xy_0, \
                             ta1_x_yyyyyz_xz_0, \
                             ta1_x_yyyyyz_yy_0, \
                             ta1_x_yyyyyz_yz_0, \
                             ta1_x_yyyyyz_zz_0, \
                             ta1_x_yyyyz_xz_0,  \
                             ta1_x_yyyyz_xz_1,  \
                             ta1_x_yyyyz_zz_0,  \
                             ta1_x_yyyyz_zz_1,  \
                             ta1_x_yyyz_xz_0,   \
                             ta1_x_yyyz_xz_1,   \
                             ta1_x_yyyz_zz_0,   \
                             ta1_x_yyyz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyyz_xx_0[i] = ta1_x_yyyyy_xx_0[i] * pa_z[i] - ta1_x_yyyyy_xx_1[i] * pc_z[i];

        ta1_x_yyyyyz_xy_0[i] = ta1_x_yyyyy_xy_0[i] * pa_z[i] - ta1_x_yyyyy_xy_1[i] * pc_z[i];

        ta1_x_yyyyyz_xz_0[i] =
            4.0 * ta1_x_yyyz_xz_0[i] * fe_0 - 4.0 * ta1_x_yyyz_xz_1[i] * fe_0 + ta1_x_yyyyz_xz_0[i] * pa_y[i] - ta1_x_yyyyz_xz_1[i] * pc_y[i];

        ta1_x_yyyyyz_yy_0[i] = ta1_x_yyyyy_yy_0[i] * pa_z[i] - ta1_x_yyyyy_yy_1[i] * pc_z[i];

        ta1_x_yyyyyz_yz_0[i] = ta1_x_yyyyy_y_0[i] * fe_0 - ta1_x_yyyyy_y_1[i] * fe_0 + ta1_x_yyyyy_yz_0[i] * pa_z[i] - ta1_x_yyyyy_yz_1[i] * pc_z[i];

        ta1_x_yyyyyz_zz_0[i] =
            4.0 * ta1_x_yyyz_zz_0[i] * fe_0 - 4.0 * ta1_x_yyyz_zz_1[i] * fe_0 + ta1_x_yyyyz_zz_0[i] * pa_y[i] - ta1_x_yyyyz_zz_1[i] * pc_y[i];
    }

    // Set up 138-144 components of targeted buffer : ID

    auto ta1_x_yyyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 138);

    auto ta1_x_yyyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 139);

    auto ta1_x_yyyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 140);

    auto ta1_x_yyyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 141);

    auto ta1_x_yyyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 142);

    auto ta1_x_yyyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 143);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_yyyy_xy_0,   \
                             ta1_x_yyyy_xy_1,   \
                             ta1_x_yyyy_yy_0,   \
                             ta1_x_yyyy_yy_1,   \
                             ta1_x_yyyyz_xy_0,  \
                             ta1_x_yyyyz_xy_1,  \
                             ta1_x_yyyyz_yy_0,  \
                             ta1_x_yyyyz_yy_1,  \
                             ta1_x_yyyyzz_xx_0, \
                             ta1_x_yyyyzz_xy_0, \
                             ta1_x_yyyyzz_xz_0, \
                             ta1_x_yyyyzz_yy_0, \
                             ta1_x_yyyyzz_yz_0, \
                             ta1_x_yyyyzz_zz_0, \
                             ta1_x_yyyzz_xx_0,  \
                             ta1_x_yyyzz_xx_1,  \
                             ta1_x_yyyzz_xz_0,  \
                             ta1_x_yyyzz_xz_1,  \
                             ta1_x_yyyzz_yz_0,  \
                             ta1_x_yyyzz_yz_1,  \
                             ta1_x_yyyzz_z_0,   \
                             ta1_x_yyyzz_z_1,   \
                             ta1_x_yyyzz_zz_0,  \
                             ta1_x_yyyzz_zz_1,  \
                             ta1_x_yyzz_xx_0,   \
                             ta1_x_yyzz_xx_1,   \
                             ta1_x_yyzz_xz_0,   \
                             ta1_x_yyzz_xz_1,   \
                             ta1_x_yyzz_yz_0,   \
                             ta1_x_yyzz_yz_1,   \
                             ta1_x_yyzz_zz_0,   \
                             ta1_x_yyzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyyzz_xx_0[i] =
            3.0 * ta1_x_yyzz_xx_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xx_1[i] * fe_0 + ta1_x_yyyzz_xx_0[i] * pa_y[i] - ta1_x_yyyzz_xx_1[i] * pc_y[i];

        ta1_x_yyyyzz_xy_0[i] = ta1_x_yyyy_xy_0[i] * fe_0 - ta1_x_yyyy_xy_1[i] * fe_0 + ta1_x_yyyyz_xy_0[i] * pa_z[i] - ta1_x_yyyyz_xy_1[i] * pc_z[i];

        ta1_x_yyyyzz_xz_0[i] =
            3.0 * ta1_x_yyzz_xz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xz_1[i] * fe_0 + ta1_x_yyyzz_xz_0[i] * pa_y[i] - ta1_x_yyyzz_xz_1[i] * pc_y[i];

        ta1_x_yyyyzz_yy_0[i] = ta1_x_yyyy_yy_0[i] * fe_0 - ta1_x_yyyy_yy_1[i] * fe_0 + ta1_x_yyyyz_yy_0[i] * pa_z[i] - ta1_x_yyyyz_yy_1[i] * pc_z[i];

        ta1_x_yyyyzz_yz_0[i] = 3.0 * ta1_x_yyzz_yz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_yz_1[i] * fe_0 + ta1_x_yyyzz_z_0[i] * fe_0 -
                               ta1_x_yyyzz_z_1[i] * fe_0 + ta1_x_yyyzz_yz_0[i] * pa_y[i] - ta1_x_yyyzz_yz_1[i] * pc_y[i];

        ta1_x_yyyyzz_zz_0[i] =
            3.0 * ta1_x_yyzz_zz_0[i] * fe_0 - 3.0 * ta1_x_yyzz_zz_1[i] * fe_0 + ta1_x_yyyzz_zz_0[i] * pa_y[i] - ta1_x_yyyzz_zz_1[i] * pc_y[i];
    }

    // Set up 144-150 components of targeted buffer : ID

    auto ta1_x_yyyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 144);

    auto ta1_x_yyyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 145);

    auto ta1_x_yyyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 146);

    auto ta1_x_yyyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 147);

    auto ta1_x_yyyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 148);

    auto ta1_x_yyyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 149);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_yyyz_xy_0,   \
                             ta1_x_yyyz_xy_1,   \
                             ta1_x_yyyz_yy_0,   \
                             ta1_x_yyyz_yy_1,   \
                             ta1_x_yyyzz_xy_0,  \
                             ta1_x_yyyzz_xy_1,  \
                             ta1_x_yyyzz_yy_0,  \
                             ta1_x_yyyzz_yy_1,  \
                             ta1_x_yyyzzz_xx_0, \
                             ta1_x_yyyzzz_xy_0, \
                             ta1_x_yyyzzz_xz_0, \
                             ta1_x_yyyzzz_yy_0, \
                             ta1_x_yyyzzz_yz_0, \
                             ta1_x_yyyzzz_zz_0, \
                             ta1_x_yyzzz_xx_0,  \
                             ta1_x_yyzzz_xx_1,  \
                             ta1_x_yyzzz_xz_0,  \
                             ta1_x_yyzzz_xz_1,  \
                             ta1_x_yyzzz_yz_0,  \
                             ta1_x_yyzzz_yz_1,  \
                             ta1_x_yyzzz_z_0,   \
                             ta1_x_yyzzz_z_1,   \
                             ta1_x_yyzzz_zz_0,  \
                             ta1_x_yyzzz_zz_1,  \
                             ta1_x_yzzz_xx_0,   \
                             ta1_x_yzzz_xx_1,   \
                             ta1_x_yzzz_xz_0,   \
                             ta1_x_yzzz_xz_1,   \
                             ta1_x_yzzz_yz_0,   \
                             ta1_x_yzzz_yz_1,   \
                             ta1_x_yzzz_zz_0,   \
                             ta1_x_yzzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyzzz_xx_0[i] =
            2.0 * ta1_x_yzzz_xx_0[i] * fe_0 - 2.0 * ta1_x_yzzz_xx_1[i] * fe_0 + ta1_x_yyzzz_xx_0[i] * pa_y[i] - ta1_x_yyzzz_xx_1[i] * pc_y[i];

        ta1_x_yyyzzz_xy_0[i] =
            2.0 * ta1_x_yyyz_xy_0[i] * fe_0 - 2.0 * ta1_x_yyyz_xy_1[i] * fe_0 + ta1_x_yyyzz_xy_0[i] * pa_z[i] - ta1_x_yyyzz_xy_1[i] * pc_z[i];

        ta1_x_yyyzzz_xz_0[i] =
            2.0 * ta1_x_yzzz_xz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_xz_1[i] * fe_0 + ta1_x_yyzzz_xz_0[i] * pa_y[i] - ta1_x_yyzzz_xz_1[i] * pc_y[i];

        ta1_x_yyyzzz_yy_0[i] =
            2.0 * ta1_x_yyyz_yy_0[i] * fe_0 - 2.0 * ta1_x_yyyz_yy_1[i] * fe_0 + ta1_x_yyyzz_yy_0[i] * pa_z[i] - ta1_x_yyyzz_yy_1[i] * pc_z[i];

        ta1_x_yyyzzz_yz_0[i] = 2.0 * ta1_x_yzzz_yz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_yz_1[i] * fe_0 + ta1_x_yyzzz_z_0[i] * fe_0 -
                               ta1_x_yyzzz_z_1[i] * fe_0 + ta1_x_yyzzz_yz_0[i] * pa_y[i] - ta1_x_yyzzz_yz_1[i] * pc_y[i];

        ta1_x_yyyzzz_zz_0[i] =
            2.0 * ta1_x_yzzz_zz_0[i] * fe_0 - 2.0 * ta1_x_yzzz_zz_1[i] * fe_0 + ta1_x_yyzzz_zz_0[i] * pa_y[i] - ta1_x_yyzzz_zz_1[i] * pc_y[i];
    }

    // Set up 150-156 components of targeted buffer : ID

    auto ta1_x_yyzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 150);

    auto ta1_x_yyzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 151);

    auto ta1_x_yyzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 152);

    auto ta1_x_yyzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 153);

    auto ta1_x_yyzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 154);

    auto ta1_x_yyzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 155);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_x_yyzz_xy_0,   \
                             ta1_x_yyzz_xy_1,   \
                             ta1_x_yyzz_yy_0,   \
                             ta1_x_yyzz_yy_1,   \
                             ta1_x_yyzzz_xy_0,  \
                             ta1_x_yyzzz_xy_1,  \
                             ta1_x_yyzzz_yy_0,  \
                             ta1_x_yyzzz_yy_1,  \
                             ta1_x_yyzzzz_xx_0, \
                             ta1_x_yyzzzz_xy_0, \
                             ta1_x_yyzzzz_xz_0, \
                             ta1_x_yyzzzz_yy_0, \
                             ta1_x_yyzzzz_yz_0, \
                             ta1_x_yyzzzz_zz_0, \
                             ta1_x_yzzzz_xx_0,  \
                             ta1_x_yzzzz_xx_1,  \
                             ta1_x_yzzzz_xz_0,  \
                             ta1_x_yzzzz_xz_1,  \
                             ta1_x_yzzzz_yz_0,  \
                             ta1_x_yzzzz_yz_1,  \
                             ta1_x_yzzzz_z_0,   \
                             ta1_x_yzzzz_z_1,   \
                             ta1_x_yzzzz_zz_0,  \
                             ta1_x_yzzzz_zz_1,  \
                             ta1_x_zzzz_xx_0,   \
                             ta1_x_zzzz_xx_1,   \
                             ta1_x_zzzz_xz_0,   \
                             ta1_x_zzzz_xz_1,   \
                             ta1_x_zzzz_yz_0,   \
                             ta1_x_zzzz_yz_1,   \
                             ta1_x_zzzz_zz_0,   \
                             ta1_x_zzzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyzzzz_xx_0[i] = ta1_x_zzzz_xx_0[i] * fe_0 - ta1_x_zzzz_xx_1[i] * fe_0 + ta1_x_yzzzz_xx_0[i] * pa_y[i] - ta1_x_yzzzz_xx_1[i] * pc_y[i];

        ta1_x_yyzzzz_xy_0[i] =
            3.0 * ta1_x_yyzz_xy_0[i] * fe_0 - 3.0 * ta1_x_yyzz_xy_1[i] * fe_0 + ta1_x_yyzzz_xy_0[i] * pa_z[i] - ta1_x_yyzzz_xy_1[i] * pc_z[i];

        ta1_x_yyzzzz_xz_0[i] = ta1_x_zzzz_xz_0[i] * fe_0 - ta1_x_zzzz_xz_1[i] * fe_0 + ta1_x_yzzzz_xz_0[i] * pa_y[i] - ta1_x_yzzzz_xz_1[i] * pc_y[i];

        ta1_x_yyzzzz_yy_0[i] =
            3.0 * ta1_x_yyzz_yy_0[i] * fe_0 - 3.0 * ta1_x_yyzz_yy_1[i] * fe_0 + ta1_x_yyzzz_yy_0[i] * pa_z[i] - ta1_x_yyzzz_yy_1[i] * pc_z[i];

        ta1_x_yyzzzz_yz_0[i] = ta1_x_zzzz_yz_0[i] * fe_0 - ta1_x_zzzz_yz_1[i] * fe_0 + ta1_x_yzzzz_z_0[i] * fe_0 - ta1_x_yzzzz_z_1[i] * fe_0 +
                               ta1_x_yzzzz_yz_0[i] * pa_y[i] - ta1_x_yzzzz_yz_1[i] * pc_y[i];

        ta1_x_yyzzzz_zz_0[i] = ta1_x_zzzz_zz_0[i] * fe_0 - ta1_x_zzzz_zz_1[i] * fe_0 + ta1_x_yzzzz_zz_0[i] * pa_y[i] - ta1_x_yzzzz_zz_1[i] * pc_y[i];
    }

    // Set up 156-162 components of targeted buffer : ID

    auto ta1_x_yzzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 156);

    auto ta1_x_yzzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 157);

    auto ta1_x_yzzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 158);

    auto ta1_x_yzzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 159);

    auto ta1_x_yzzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 160);

    auto ta1_x_yzzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 161);

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_x_yzzzzz_xx_0, \
                             ta1_x_yzzzzz_xy_0, \
                             ta1_x_yzzzzz_xz_0, \
                             ta1_x_yzzzzz_yy_0, \
                             ta1_x_yzzzzz_yz_0, \
                             ta1_x_yzzzzz_zz_0, \
                             ta1_x_zzzzz_x_0,   \
                             ta1_x_zzzzz_x_1,   \
                             ta1_x_zzzzz_xx_0,  \
                             ta1_x_zzzzz_xx_1,  \
                             ta1_x_zzzzz_xy_0,  \
                             ta1_x_zzzzz_xy_1,  \
                             ta1_x_zzzzz_xz_0,  \
                             ta1_x_zzzzz_xz_1,  \
                             ta1_x_zzzzz_y_0,   \
                             ta1_x_zzzzz_y_1,   \
                             ta1_x_zzzzz_yy_0,  \
                             ta1_x_zzzzz_yy_1,  \
                             ta1_x_zzzzz_yz_0,  \
                             ta1_x_zzzzz_yz_1,  \
                             ta1_x_zzzzz_z_0,   \
                             ta1_x_zzzzz_z_1,   \
                             ta1_x_zzzzz_zz_0,  \
                             ta1_x_zzzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzzzzz_xx_0[i] = ta1_x_zzzzz_xx_0[i] * pa_y[i] - ta1_x_zzzzz_xx_1[i] * pc_y[i];

        ta1_x_yzzzzz_xy_0[i] = ta1_x_zzzzz_x_0[i] * fe_0 - ta1_x_zzzzz_x_1[i] * fe_0 + ta1_x_zzzzz_xy_0[i] * pa_y[i] - ta1_x_zzzzz_xy_1[i] * pc_y[i];

        ta1_x_yzzzzz_xz_0[i] = ta1_x_zzzzz_xz_0[i] * pa_y[i] - ta1_x_zzzzz_xz_1[i] * pc_y[i];

        ta1_x_yzzzzz_yy_0[i] =
            2.0 * ta1_x_zzzzz_y_0[i] * fe_0 - 2.0 * ta1_x_zzzzz_y_1[i] * fe_0 + ta1_x_zzzzz_yy_0[i] * pa_y[i] - ta1_x_zzzzz_yy_1[i] * pc_y[i];

        ta1_x_yzzzzz_yz_0[i] = ta1_x_zzzzz_z_0[i] * fe_0 - ta1_x_zzzzz_z_1[i] * fe_0 + ta1_x_zzzzz_yz_0[i] * pa_y[i] - ta1_x_zzzzz_yz_1[i] * pc_y[i];

        ta1_x_yzzzzz_zz_0[i] = ta1_x_zzzzz_zz_0[i] * pa_y[i] - ta1_x_zzzzz_zz_1[i] * pc_y[i];
    }

    // Set up 162-168 components of targeted buffer : ID

    auto ta1_x_zzzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 162);

    auto ta1_x_zzzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 163);

    auto ta1_x_zzzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 164);

    auto ta1_x_zzzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 165);

    auto ta1_x_zzzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 166);

    auto ta1_x_zzzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 167);

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta1_x_zzzz_xx_0,   \
                             ta1_x_zzzz_xx_1,   \
                             ta1_x_zzzz_xy_0,   \
                             ta1_x_zzzz_xy_1,   \
                             ta1_x_zzzz_xz_0,   \
                             ta1_x_zzzz_xz_1,   \
                             ta1_x_zzzz_yy_0,   \
                             ta1_x_zzzz_yy_1,   \
                             ta1_x_zzzz_yz_0,   \
                             ta1_x_zzzz_yz_1,   \
                             ta1_x_zzzz_zz_0,   \
                             ta1_x_zzzz_zz_1,   \
                             ta1_x_zzzzz_x_0,   \
                             ta1_x_zzzzz_x_1,   \
                             ta1_x_zzzzz_xx_0,  \
                             ta1_x_zzzzz_xx_1,  \
                             ta1_x_zzzzz_xy_0,  \
                             ta1_x_zzzzz_xy_1,  \
                             ta1_x_zzzzz_xz_0,  \
                             ta1_x_zzzzz_xz_1,  \
                             ta1_x_zzzzz_y_0,   \
                             ta1_x_zzzzz_y_1,   \
                             ta1_x_zzzzz_yy_0,  \
                             ta1_x_zzzzz_yy_1,  \
                             ta1_x_zzzzz_yz_0,  \
                             ta1_x_zzzzz_yz_1,  \
                             ta1_x_zzzzz_z_0,   \
                             ta1_x_zzzzz_z_1,   \
                             ta1_x_zzzzz_zz_0,  \
                             ta1_x_zzzzz_zz_1,  \
                             ta1_x_zzzzzz_xx_0, \
                             ta1_x_zzzzzz_xy_0, \
                             ta1_x_zzzzzz_xz_0, \
                             ta1_x_zzzzzz_yy_0, \
                             ta1_x_zzzzzz_yz_0, \
                             ta1_x_zzzzzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzzzzz_xx_0[i] =
            5.0 * ta1_x_zzzz_xx_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xx_1[i] * fe_0 + ta1_x_zzzzz_xx_0[i] * pa_z[i] - ta1_x_zzzzz_xx_1[i] * pc_z[i];

        ta1_x_zzzzzz_xy_0[i] =
            5.0 * ta1_x_zzzz_xy_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xy_1[i] * fe_0 + ta1_x_zzzzz_xy_0[i] * pa_z[i] - ta1_x_zzzzz_xy_1[i] * pc_z[i];

        ta1_x_zzzzzz_xz_0[i] = 5.0 * ta1_x_zzzz_xz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_xz_1[i] * fe_0 + ta1_x_zzzzz_x_0[i] * fe_0 -
                               ta1_x_zzzzz_x_1[i] * fe_0 + ta1_x_zzzzz_xz_0[i] * pa_z[i] - ta1_x_zzzzz_xz_1[i] * pc_z[i];

        ta1_x_zzzzzz_yy_0[i] =
            5.0 * ta1_x_zzzz_yy_0[i] * fe_0 - 5.0 * ta1_x_zzzz_yy_1[i] * fe_0 + ta1_x_zzzzz_yy_0[i] * pa_z[i] - ta1_x_zzzzz_yy_1[i] * pc_z[i];

        ta1_x_zzzzzz_yz_0[i] = 5.0 * ta1_x_zzzz_yz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_yz_1[i] * fe_0 + ta1_x_zzzzz_y_0[i] * fe_0 -
                               ta1_x_zzzzz_y_1[i] * fe_0 + ta1_x_zzzzz_yz_0[i] * pa_z[i] - ta1_x_zzzzz_yz_1[i] * pc_z[i];

        ta1_x_zzzzzz_zz_0[i] = 5.0 * ta1_x_zzzz_zz_0[i] * fe_0 - 5.0 * ta1_x_zzzz_zz_1[i] * fe_0 + 2.0 * ta1_x_zzzzz_z_0[i] * fe_0 -
                               2.0 * ta1_x_zzzzz_z_1[i] * fe_0 + ta1_x_zzzzz_zz_0[i] * pa_z[i] - ta1_x_zzzzz_zz_1[i] * pc_z[i];
    }

    // Set up 168-174 components of targeted buffer : ID

    auto ta1_y_xxxxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 168);

    auto ta1_y_xxxxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 169);

    auto ta1_y_xxxxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 170);

    auto ta1_y_xxxxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 171);

    auto ta1_y_xxxxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 172);

    auto ta1_y_xxxxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 173);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_y_xxxx_xx_0,   \
                             ta1_y_xxxx_xx_1,   \
                             ta1_y_xxxx_xy_0,   \
                             ta1_y_xxxx_xy_1,   \
                             ta1_y_xxxx_xz_0,   \
                             ta1_y_xxxx_xz_1,   \
                             ta1_y_xxxx_yy_0,   \
                             ta1_y_xxxx_yy_1,   \
                             ta1_y_xxxx_yz_0,   \
                             ta1_y_xxxx_yz_1,   \
                             ta1_y_xxxx_zz_0,   \
                             ta1_y_xxxx_zz_1,   \
                             ta1_y_xxxxx_x_0,   \
                             ta1_y_xxxxx_x_1,   \
                             ta1_y_xxxxx_xx_0,  \
                             ta1_y_xxxxx_xx_1,  \
                             ta1_y_xxxxx_xy_0,  \
                             ta1_y_xxxxx_xy_1,  \
                             ta1_y_xxxxx_xz_0,  \
                             ta1_y_xxxxx_xz_1,  \
                             ta1_y_xxxxx_y_0,   \
                             ta1_y_xxxxx_y_1,   \
                             ta1_y_xxxxx_yy_0,  \
                             ta1_y_xxxxx_yy_1,  \
                             ta1_y_xxxxx_yz_0,  \
                             ta1_y_xxxxx_yz_1,  \
                             ta1_y_xxxxx_z_0,   \
                             ta1_y_xxxxx_z_1,   \
                             ta1_y_xxxxx_zz_0,  \
                             ta1_y_xxxxx_zz_1,  \
                             ta1_y_xxxxxx_xx_0, \
                             ta1_y_xxxxxx_xy_0, \
                             ta1_y_xxxxxx_xz_0, \
                             ta1_y_xxxxxx_yy_0, \
                             ta1_y_xxxxxx_yz_0, \
                             ta1_y_xxxxxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxxx_xx_0[i] = 5.0 * ta1_y_xxxx_xx_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xx_1[i] * fe_0 + 2.0 * ta1_y_xxxxx_x_0[i] * fe_0 -
                               2.0 * ta1_y_xxxxx_x_1[i] * fe_0 + ta1_y_xxxxx_xx_0[i] * pa_x[i] - ta1_y_xxxxx_xx_1[i] * pc_x[i];

        ta1_y_xxxxxx_xy_0[i] = 5.0 * ta1_y_xxxx_xy_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xy_1[i] * fe_0 + ta1_y_xxxxx_y_0[i] * fe_0 -
                               ta1_y_xxxxx_y_1[i] * fe_0 + ta1_y_xxxxx_xy_0[i] * pa_x[i] - ta1_y_xxxxx_xy_1[i] * pc_x[i];

        ta1_y_xxxxxx_xz_0[i] = 5.0 * ta1_y_xxxx_xz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_xz_1[i] * fe_0 + ta1_y_xxxxx_z_0[i] * fe_0 -
                               ta1_y_xxxxx_z_1[i] * fe_0 + ta1_y_xxxxx_xz_0[i] * pa_x[i] - ta1_y_xxxxx_xz_1[i] * pc_x[i];

        ta1_y_xxxxxx_yy_0[i] =
            5.0 * ta1_y_xxxx_yy_0[i] * fe_0 - 5.0 * ta1_y_xxxx_yy_1[i] * fe_0 + ta1_y_xxxxx_yy_0[i] * pa_x[i] - ta1_y_xxxxx_yy_1[i] * pc_x[i];

        ta1_y_xxxxxx_yz_0[i] =
            5.0 * ta1_y_xxxx_yz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_yz_1[i] * fe_0 + ta1_y_xxxxx_yz_0[i] * pa_x[i] - ta1_y_xxxxx_yz_1[i] * pc_x[i];

        ta1_y_xxxxxx_zz_0[i] =
            5.0 * ta1_y_xxxx_zz_0[i] * fe_0 - 5.0 * ta1_y_xxxx_zz_1[i] * fe_0 + ta1_y_xxxxx_zz_0[i] * pa_x[i] - ta1_y_xxxxx_zz_1[i] * pc_x[i];
    }

    // Set up 174-180 components of targeted buffer : ID

    auto ta1_y_xxxxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 174);

    auto ta1_y_xxxxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 175);

    auto ta1_y_xxxxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 176);

    auto ta1_y_xxxxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 177);

    auto ta1_y_xxxxxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 178);

    auto ta1_y_xxxxxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 179);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_y_xxxxx_x_0,   \
                             ta1_y_xxxxx_x_1,   \
                             ta1_y_xxxxx_xx_0,  \
                             ta1_y_xxxxx_xx_1,  \
                             ta1_y_xxxxx_xy_0,  \
                             ta1_y_xxxxx_xy_1,  \
                             ta1_y_xxxxx_xz_0,  \
                             ta1_y_xxxxx_xz_1,  \
                             ta1_y_xxxxx_zz_0,  \
                             ta1_y_xxxxx_zz_1,  \
                             ta1_y_xxxxxy_xx_0, \
                             ta1_y_xxxxxy_xy_0, \
                             ta1_y_xxxxxy_xz_0, \
                             ta1_y_xxxxxy_yy_0, \
                             ta1_y_xxxxxy_yz_0, \
                             ta1_y_xxxxxy_zz_0, \
                             ta1_y_xxxxy_yy_0,  \
                             ta1_y_xxxxy_yy_1,  \
                             ta1_y_xxxxy_yz_0,  \
                             ta1_y_xxxxy_yz_1,  \
                             ta1_y_xxxy_yy_0,   \
                             ta1_y_xxxy_yy_1,   \
                             ta1_y_xxxy_yz_0,   \
                             ta1_y_xxxy_yz_1,   \
                             ta_xxxxx_xx_1,     \
                             ta_xxxxx_xy_1,     \
                             ta_xxxxx_xz_1,     \
                             ta_xxxxx_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxxy_xx_0[i] = ta_xxxxx_xx_1[i] + ta1_y_xxxxx_xx_0[i] * pa_y[i] - ta1_y_xxxxx_xx_1[i] * pc_y[i];

        ta1_y_xxxxxy_xy_0[i] =
            ta1_y_xxxxx_x_0[i] * fe_0 - ta1_y_xxxxx_x_1[i] * fe_0 + ta_xxxxx_xy_1[i] + ta1_y_xxxxx_xy_0[i] * pa_y[i] - ta1_y_xxxxx_xy_1[i] * pc_y[i];

        ta1_y_xxxxxy_xz_0[i] = ta_xxxxx_xz_1[i] + ta1_y_xxxxx_xz_0[i] * pa_y[i] - ta1_y_xxxxx_xz_1[i] * pc_y[i];

        ta1_y_xxxxxy_yy_0[i] =
            4.0 * ta1_y_xxxy_yy_0[i] * fe_0 - 4.0 * ta1_y_xxxy_yy_1[i] * fe_0 + ta1_y_xxxxy_yy_0[i] * pa_x[i] - ta1_y_xxxxy_yy_1[i] * pc_x[i];

        ta1_y_xxxxxy_yz_0[i] =
            4.0 * ta1_y_xxxy_yz_0[i] * fe_0 - 4.0 * ta1_y_xxxy_yz_1[i] * fe_0 + ta1_y_xxxxy_yz_0[i] * pa_x[i] - ta1_y_xxxxy_yz_1[i] * pc_x[i];

        ta1_y_xxxxxy_zz_0[i] = ta_xxxxx_zz_1[i] + ta1_y_xxxxx_zz_0[i] * pa_y[i] - ta1_y_xxxxx_zz_1[i] * pc_y[i];
    }

    // Set up 180-186 components of targeted buffer : ID

    auto ta1_y_xxxxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 180);

    auto ta1_y_xxxxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 181);

    auto ta1_y_xxxxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 182);

    auto ta1_y_xxxxxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 183);

    auto ta1_y_xxxxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 184);

    auto ta1_y_xxxxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 185);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_y_xxxxx_x_0,   \
                             ta1_y_xxxxx_x_1,   \
                             ta1_y_xxxxx_xx_0,  \
                             ta1_y_xxxxx_xx_1,  \
                             ta1_y_xxxxx_xy_0,  \
                             ta1_y_xxxxx_xy_1,  \
                             ta1_y_xxxxx_xz_0,  \
                             ta1_y_xxxxx_xz_1,  \
                             ta1_y_xxxxx_yy_0,  \
                             ta1_y_xxxxx_yy_1,  \
                             ta1_y_xxxxxz_xx_0, \
                             ta1_y_xxxxxz_xy_0, \
                             ta1_y_xxxxxz_xz_0, \
                             ta1_y_xxxxxz_yy_0, \
                             ta1_y_xxxxxz_yz_0, \
                             ta1_y_xxxxxz_zz_0, \
                             ta1_y_xxxxz_yz_0,  \
                             ta1_y_xxxxz_yz_1,  \
                             ta1_y_xxxxz_zz_0,  \
                             ta1_y_xxxxz_zz_1,  \
                             ta1_y_xxxz_yz_0,   \
                             ta1_y_xxxz_yz_1,   \
                             ta1_y_xxxz_zz_0,   \
                             ta1_y_xxxz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxxz_xx_0[i] = ta1_y_xxxxx_xx_0[i] * pa_z[i] - ta1_y_xxxxx_xx_1[i] * pc_z[i];

        ta1_y_xxxxxz_xy_0[i] = ta1_y_xxxxx_xy_0[i] * pa_z[i] - ta1_y_xxxxx_xy_1[i] * pc_z[i];

        ta1_y_xxxxxz_xz_0[i] = ta1_y_xxxxx_x_0[i] * fe_0 - ta1_y_xxxxx_x_1[i] * fe_0 + ta1_y_xxxxx_xz_0[i] * pa_z[i] - ta1_y_xxxxx_xz_1[i] * pc_z[i];

        ta1_y_xxxxxz_yy_0[i] = ta1_y_xxxxx_yy_0[i] * pa_z[i] - ta1_y_xxxxx_yy_1[i] * pc_z[i];

        ta1_y_xxxxxz_yz_0[i] =
            4.0 * ta1_y_xxxz_yz_0[i] * fe_0 - 4.0 * ta1_y_xxxz_yz_1[i] * fe_0 + ta1_y_xxxxz_yz_0[i] * pa_x[i] - ta1_y_xxxxz_yz_1[i] * pc_x[i];

        ta1_y_xxxxxz_zz_0[i] =
            4.0 * ta1_y_xxxz_zz_0[i] * fe_0 - 4.0 * ta1_y_xxxz_zz_1[i] * fe_0 + ta1_y_xxxxz_zz_0[i] * pa_x[i] - ta1_y_xxxxz_zz_1[i] * pc_x[i];
    }

    // Set up 186-192 components of targeted buffer : ID

    auto ta1_y_xxxxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 186);

    auto ta1_y_xxxxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 187);

    auto ta1_y_xxxxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 188);

    auto ta1_y_xxxxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 189);

    auto ta1_y_xxxxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 190);

    auto ta1_y_xxxxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 191);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_y_xxxx_xx_0,   \
                             ta1_y_xxxx_xx_1,   \
                             ta1_y_xxxx_xz_0,   \
                             ta1_y_xxxx_xz_1,   \
                             ta1_y_xxxxy_xx_0,  \
                             ta1_y_xxxxy_xx_1,  \
                             ta1_y_xxxxy_xz_0,  \
                             ta1_y_xxxxy_xz_1,  \
                             ta1_y_xxxxyy_xx_0, \
                             ta1_y_xxxxyy_xy_0, \
                             ta1_y_xxxxyy_xz_0, \
                             ta1_y_xxxxyy_yy_0, \
                             ta1_y_xxxxyy_yz_0, \
                             ta1_y_xxxxyy_zz_0, \
                             ta1_y_xxxyy_xy_0,  \
                             ta1_y_xxxyy_xy_1,  \
                             ta1_y_xxxyy_y_0,   \
                             ta1_y_xxxyy_y_1,   \
                             ta1_y_xxxyy_yy_0,  \
                             ta1_y_xxxyy_yy_1,  \
                             ta1_y_xxxyy_yz_0,  \
                             ta1_y_xxxyy_yz_1,  \
                             ta1_y_xxxyy_zz_0,  \
                             ta1_y_xxxyy_zz_1,  \
                             ta1_y_xxyy_xy_0,   \
                             ta1_y_xxyy_xy_1,   \
                             ta1_y_xxyy_yy_0,   \
                             ta1_y_xxyy_yy_1,   \
                             ta1_y_xxyy_yz_0,   \
                             ta1_y_xxyy_yz_1,   \
                             ta1_y_xxyy_zz_0,   \
                             ta1_y_xxyy_zz_1,   \
                             ta_xxxxy_xx_1,     \
                             ta_xxxxy_xz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxyy_xx_0[i] =
            ta1_y_xxxx_xx_0[i] * fe_0 - ta1_y_xxxx_xx_1[i] * fe_0 + ta_xxxxy_xx_1[i] + ta1_y_xxxxy_xx_0[i] * pa_y[i] - ta1_y_xxxxy_xx_1[i] * pc_y[i];

        ta1_y_xxxxyy_xy_0[i] = 3.0 * ta1_y_xxyy_xy_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xy_1[i] * fe_0 + ta1_y_xxxyy_y_0[i] * fe_0 -
                               ta1_y_xxxyy_y_1[i] * fe_0 + ta1_y_xxxyy_xy_0[i] * pa_x[i] - ta1_y_xxxyy_xy_1[i] * pc_x[i];

        ta1_y_xxxxyy_xz_0[i] =
            ta1_y_xxxx_xz_0[i] * fe_0 - ta1_y_xxxx_xz_1[i] * fe_0 + ta_xxxxy_xz_1[i] + ta1_y_xxxxy_xz_0[i] * pa_y[i] - ta1_y_xxxxy_xz_1[i] * pc_y[i];

        ta1_y_xxxxyy_yy_0[i] =
            3.0 * ta1_y_xxyy_yy_0[i] * fe_0 - 3.0 * ta1_y_xxyy_yy_1[i] * fe_0 + ta1_y_xxxyy_yy_0[i] * pa_x[i] - ta1_y_xxxyy_yy_1[i] * pc_x[i];

        ta1_y_xxxxyy_yz_0[i] =
            3.0 * ta1_y_xxyy_yz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_yz_1[i] * fe_0 + ta1_y_xxxyy_yz_0[i] * pa_x[i] - ta1_y_xxxyy_yz_1[i] * pc_x[i];

        ta1_y_xxxxyy_zz_0[i] =
            3.0 * ta1_y_xxyy_zz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_zz_1[i] * fe_0 + ta1_y_xxxyy_zz_0[i] * pa_x[i] - ta1_y_xxxyy_zz_1[i] * pc_x[i];
    }

    // Set up 192-198 components of targeted buffer : ID

    auto ta1_y_xxxxyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 192);

    auto ta1_y_xxxxyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 193);

    auto ta1_y_xxxxyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 194);

    auto ta1_y_xxxxyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 195);

    auto ta1_y_xxxxyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 196);

    auto ta1_y_xxxxyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 197);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_y_xxxxy_xx_0,  \
                             ta1_y_xxxxy_xx_1,  \
                             ta1_y_xxxxy_xy_0,  \
                             ta1_y_xxxxy_xy_1,  \
                             ta1_y_xxxxy_yy_0,  \
                             ta1_y_xxxxy_yy_1,  \
                             ta1_y_xxxxyz_xx_0, \
                             ta1_y_xxxxyz_xy_0, \
                             ta1_y_xxxxyz_xz_0, \
                             ta1_y_xxxxyz_yy_0, \
                             ta1_y_xxxxyz_yz_0, \
                             ta1_y_xxxxyz_zz_0, \
                             ta1_y_xxxxz_xz_0,  \
                             ta1_y_xxxxz_xz_1,  \
                             ta1_y_xxxxz_zz_0,  \
                             ta1_y_xxxxz_zz_1,  \
                             ta1_y_xxxyz_yz_0,  \
                             ta1_y_xxxyz_yz_1,  \
                             ta1_y_xxyz_yz_0,   \
                             ta1_y_xxyz_yz_1,   \
                             ta_xxxxz_xz_1,     \
                             ta_xxxxz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxyz_xx_0[i] = ta1_y_xxxxy_xx_0[i] * pa_z[i] - ta1_y_xxxxy_xx_1[i] * pc_z[i];

        ta1_y_xxxxyz_xy_0[i] = ta1_y_xxxxy_xy_0[i] * pa_z[i] - ta1_y_xxxxy_xy_1[i] * pc_z[i];

        ta1_y_xxxxyz_xz_0[i] = ta_xxxxz_xz_1[i] + ta1_y_xxxxz_xz_0[i] * pa_y[i] - ta1_y_xxxxz_xz_1[i] * pc_y[i];

        ta1_y_xxxxyz_yy_0[i] = ta1_y_xxxxy_yy_0[i] * pa_z[i] - ta1_y_xxxxy_yy_1[i] * pc_z[i];

        ta1_y_xxxxyz_yz_0[i] =
            3.0 * ta1_y_xxyz_yz_0[i] * fe_0 - 3.0 * ta1_y_xxyz_yz_1[i] * fe_0 + ta1_y_xxxyz_yz_0[i] * pa_x[i] - ta1_y_xxxyz_yz_1[i] * pc_x[i];

        ta1_y_xxxxyz_zz_0[i] = ta_xxxxz_zz_1[i] + ta1_y_xxxxz_zz_0[i] * pa_y[i] - ta1_y_xxxxz_zz_1[i] * pc_y[i];
    }

    // Set up 198-204 components of targeted buffer : ID

    auto ta1_y_xxxxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 198);

    auto ta1_y_xxxxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 199);

    auto ta1_y_xxxxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 200);

    auto ta1_y_xxxxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 201);

    auto ta1_y_xxxxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 202);

    auto ta1_y_xxxxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 203);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_y_xxxx_xx_0,   \
                             ta1_y_xxxx_xx_1,   \
                             ta1_y_xxxx_xy_0,   \
                             ta1_y_xxxx_xy_1,   \
                             ta1_y_xxxxz_xx_0,  \
                             ta1_y_xxxxz_xx_1,  \
                             ta1_y_xxxxz_xy_0,  \
                             ta1_y_xxxxz_xy_1,  \
                             ta1_y_xxxxzz_xx_0, \
                             ta1_y_xxxxzz_xy_0, \
                             ta1_y_xxxxzz_xz_0, \
                             ta1_y_xxxxzz_yy_0, \
                             ta1_y_xxxxzz_yz_0, \
                             ta1_y_xxxxzz_zz_0, \
                             ta1_y_xxxzz_xz_0,  \
                             ta1_y_xxxzz_xz_1,  \
                             ta1_y_xxxzz_yy_0,  \
                             ta1_y_xxxzz_yy_1,  \
                             ta1_y_xxxzz_yz_0,  \
                             ta1_y_xxxzz_yz_1,  \
                             ta1_y_xxxzz_z_0,   \
                             ta1_y_xxxzz_z_1,   \
                             ta1_y_xxxzz_zz_0,  \
                             ta1_y_xxxzz_zz_1,  \
                             ta1_y_xxzz_xz_0,   \
                             ta1_y_xxzz_xz_1,   \
                             ta1_y_xxzz_yy_0,   \
                             ta1_y_xxzz_yy_1,   \
                             ta1_y_xxzz_yz_0,   \
                             ta1_y_xxzz_yz_1,   \
                             ta1_y_xxzz_zz_0,   \
                             ta1_y_xxzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxxzz_xx_0[i] = ta1_y_xxxx_xx_0[i] * fe_0 - ta1_y_xxxx_xx_1[i] * fe_0 + ta1_y_xxxxz_xx_0[i] * pa_z[i] - ta1_y_xxxxz_xx_1[i] * pc_z[i];

        ta1_y_xxxxzz_xy_0[i] = ta1_y_xxxx_xy_0[i] * fe_0 - ta1_y_xxxx_xy_1[i] * fe_0 + ta1_y_xxxxz_xy_0[i] * pa_z[i] - ta1_y_xxxxz_xy_1[i] * pc_z[i];

        ta1_y_xxxxzz_xz_0[i] = 3.0 * ta1_y_xxzz_xz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xz_1[i] * fe_0 + ta1_y_xxxzz_z_0[i] * fe_0 -
                               ta1_y_xxxzz_z_1[i] * fe_0 + ta1_y_xxxzz_xz_0[i] * pa_x[i] - ta1_y_xxxzz_xz_1[i] * pc_x[i];

        ta1_y_xxxxzz_yy_0[i] =
            3.0 * ta1_y_xxzz_yy_0[i] * fe_0 - 3.0 * ta1_y_xxzz_yy_1[i] * fe_0 + ta1_y_xxxzz_yy_0[i] * pa_x[i] - ta1_y_xxxzz_yy_1[i] * pc_x[i];

        ta1_y_xxxxzz_yz_0[i] =
            3.0 * ta1_y_xxzz_yz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_yz_1[i] * fe_0 + ta1_y_xxxzz_yz_0[i] * pa_x[i] - ta1_y_xxxzz_yz_1[i] * pc_x[i];

        ta1_y_xxxxzz_zz_0[i] =
            3.0 * ta1_y_xxzz_zz_0[i] * fe_0 - 3.0 * ta1_y_xxzz_zz_1[i] * fe_0 + ta1_y_xxxzz_zz_0[i] * pa_x[i] - ta1_y_xxxzz_zz_1[i] * pc_x[i];
    }

    // Set up 204-210 components of targeted buffer : ID

    auto ta1_y_xxxyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 204);

    auto ta1_y_xxxyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 205);

    auto ta1_y_xxxyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 206);

    auto ta1_y_xxxyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 207);

    auto ta1_y_xxxyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 208);

    auto ta1_y_xxxyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 209);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_y_xxxy_xx_0,   \
                             ta1_y_xxxy_xx_1,   \
                             ta1_y_xxxy_xz_0,   \
                             ta1_y_xxxy_xz_1,   \
                             ta1_y_xxxyy_xx_0,  \
                             ta1_y_xxxyy_xx_1,  \
                             ta1_y_xxxyy_xz_0,  \
                             ta1_y_xxxyy_xz_1,  \
                             ta1_y_xxxyyy_xx_0, \
                             ta1_y_xxxyyy_xy_0, \
                             ta1_y_xxxyyy_xz_0, \
                             ta1_y_xxxyyy_yy_0, \
                             ta1_y_xxxyyy_yz_0, \
                             ta1_y_xxxyyy_zz_0, \
                             ta1_y_xxyyy_xy_0,  \
                             ta1_y_xxyyy_xy_1,  \
                             ta1_y_xxyyy_y_0,   \
                             ta1_y_xxyyy_y_1,   \
                             ta1_y_xxyyy_yy_0,  \
                             ta1_y_xxyyy_yy_1,  \
                             ta1_y_xxyyy_yz_0,  \
                             ta1_y_xxyyy_yz_1,  \
                             ta1_y_xxyyy_zz_0,  \
                             ta1_y_xxyyy_zz_1,  \
                             ta1_y_xyyy_xy_0,   \
                             ta1_y_xyyy_xy_1,   \
                             ta1_y_xyyy_yy_0,   \
                             ta1_y_xyyy_yy_1,   \
                             ta1_y_xyyy_yz_0,   \
                             ta1_y_xyyy_yz_1,   \
                             ta1_y_xyyy_zz_0,   \
                             ta1_y_xyyy_zz_1,   \
                             ta_xxxyy_xx_1,     \
                             ta_xxxyy_xz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyyy_xx_0[i] = 2.0 * ta1_y_xxxy_xx_0[i] * fe_0 - 2.0 * ta1_y_xxxy_xx_1[i] * fe_0 + ta_xxxyy_xx_1[i] + ta1_y_xxxyy_xx_0[i] * pa_y[i] -
                               ta1_y_xxxyy_xx_1[i] * pc_y[i];

        ta1_y_xxxyyy_xy_0[i] = 2.0 * ta1_y_xyyy_xy_0[i] * fe_0 - 2.0 * ta1_y_xyyy_xy_1[i] * fe_0 + ta1_y_xxyyy_y_0[i] * fe_0 -
                               ta1_y_xxyyy_y_1[i] * fe_0 + ta1_y_xxyyy_xy_0[i] * pa_x[i] - ta1_y_xxyyy_xy_1[i] * pc_x[i];

        ta1_y_xxxyyy_xz_0[i] = 2.0 * ta1_y_xxxy_xz_0[i] * fe_0 - 2.0 * ta1_y_xxxy_xz_1[i] * fe_0 + ta_xxxyy_xz_1[i] + ta1_y_xxxyy_xz_0[i] * pa_y[i] -
                               ta1_y_xxxyy_xz_1[i] * pc_y[i];

        ta1_y_xxxyyy_yy_0[i] =
            2.0 * ta1_y_xyyy_yy_0[i] * fe_0 - 2.0 * ta1_y_xyyy_yy_1[i] * fe_0 + ta1_y_xxyyy_yy_0[i] * pa_x[i] - ta1_y_xxyyy_yy_1[i] * pc_x[i];

        ta1_y_xxxyyy_yz_0[i] =
            2.0 * ta1_y_xyyy_yz_0[i] * fe_0 - 2.0 * ta1_y_xyyy_yz_1[i] * fe_0 + ta1_y_xxyyy_yz_0[i] * pa_x[i] - ta1_y_xxyyy_yz_1[i] * pc_x[i];

        ta1_y_xxxyyy_zz_0[i] =
            2.0 * ta1_y_xyyy_zz_0[i] * fe_0 - 2.0 * ta1_y_xyyy_zz_1[i] * fe_0 + ta1_y_xxyyy_zz_0[i] * pa_x[i] - ta1_y_xxyyy_zz_1[i] * pc_x[i];
    }

    // Set up 210-216 components of targeted buffer : ID

    auto ta1_y_xxxyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 210);

    auto ta1_y_xxxyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 211);

    auto ta1_y_xxxyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 212);

    auto ta1_y_xxxyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 213);

    auto ta1_y_xxxyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 214);

    auto ta1_y_xxxyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 215);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_y_xxxyy_x_0,   \
                             ta1_y_xxxyy_x_1,   \
                             ta1_y_xxxyy_xx_0,  \
                             ta1_y_xxxyy_xx_1,  \
                             ta1_y_xxxyy_xy_0,  \
                             ta1_y_xxxyy_xy_1,  \
                             ta1_y_xxxyy_xz_0,  \
                             ta1_y_xxxyy_xz_1,  \
                             ta1_y_xxxyy_yy_0,  \
                             ta1_y_xxxyy_yy_1,  \
                             ta1_y_xxxyyz_xx_0, \
                             ta1_y_xxxyyz_xy_0, \
                             ta1_y_xxxyyz_xz_0, \
                             ta1_y_xxxyyz_yy_0, \
                             ta1_y_xxxyyz_yz_0, \
                             ta1_y_xxxyyz_zz_0, \
                             ta1_y_xxyyz_yz_0,  \
                             ta1_y_xxyyz_yz_1,  \
                             ta1_y_xxyyz_zz_0,  \
                             ta1_y_xxyyz_zz_1,  \
                             ta1_y_xyyz_yz_0,   \
                             ta1_y_xyyz_yz_1,   \
                             ta1_y_xyyz_zz_0,   \
                             ta1_y_xyyz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyyz_xx_0[i] = ta1_y_xxxyy_xx_0[i] * pa_z[i] - ta1_y_xxxyy_xx_1[i] * pc_z[i];

        ta1_y_xxxyyz_xy_0[i] = ta1_y_xxxyy_xy_0[i] * pa_z[i] - ta1_y_xxxyy_xy_1[i] * pc_z[i];

        ta1_y_xxxyyz_xz_0[i] = ta1_y_xxxyy_x_0[i] * fe_0 - ta1_y_xxxyy_x_1[i] * fe_0 + ta1_y_xxxyy_xz_0[i] * pa_z[i] - ta1_y_xxxyy_xz_1[i] * pc_z[i];

        ta1_y_xxxyyz_yy_0[i] = ta1_y_xxxyy_yy_0[i] * pa_z[i] - ta1_y_xxxyy_yy_1[i] * pc_z[i];

        ta1_y_xxxyyz_yz_0[i] =
            2.0 * ta1_y_xyyz_yz_0[i] * fe_0 - 2.0 * ta1_y_xyyz_yz_1[i] * fe_0 + ta1_y_xxyyz_yz_0[i] * pa_x[i] - ta1_y_xxyyz_yz_1[i] * pc_x[i];

        ta1_y_xxxyyz_zz_0[i] =
            2.0 * ta1_y_xyyz_zz_0[i] * fe_0 - 2.0 * ta1_y_xyyz_zz_1[i] * fe_0 + ta1_y_xxyyz_zz_0[i] * pa_x[i] - ta1_y_xxyyz_zz_1[i] * pc_x[i];
    }

    // Set up 216-222 components of targeted buffer : ID

    auto ta1_y_xxxyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 216);

    auto ta1_y_xxxyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 217);

    auto ta1_y_xxxyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 218);

    auto ta1_y_xxxyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 219);

    auto ta1_y_xxxyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 220);

    auto ta1_y_xxxyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 221);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_y_xxxy_xy_0,   \
                             ta1_y_xxxy_xy_1,   \
                             ta1_y_xxxyz_xy_0,  \
                             ta1_y_xxxyz_xy_1,  \
                             ta1_y_xxxyzz_xx_0, \
                             ta1_y_xxxyzz_xy_0, \
                             ta1_y_xxxyzz_xz_0, \
                             ta1_y_xxxyzz_yy_0, \
                             ta1_y_xxxyzz_yz_0, \
                             ta1_y_xxxyzz_zz_0, \
                             ta1_y_xxxzz_xx_0,  \
                             ta1_y_xxxzz_xx_1,  \
                             ta1_y_xxxzz_xz_0,  \
                             ta1_y_xxxzz_xz_1,  \
                             ta1_y_xxxzz_zz_0,  \
                             ta1_y_xxxzz_zz_1,  \
                             ta1_y_xxyzz_yy_0,  \
                             ta1_y_xxyzz_yy_1,  \
                             ta1_y_xxyzz_yz_0,  \
                             ta1_y_xxyzz_yz_1,  \
                             ta1_y_xyzz_yy_0,   \
                             ta1_y_xyzz_yy_1,   \
                             ta1_y_xyzz_yz_0,   \
                             ta1_y_xyzz_yz_1,   \
                             ta_xxxzz_xx_1,     \
                             ta_xxxzz_xz_1,     \
                             ta_xxxzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxyzz_xx_0[i] = ta_xxxzz_xx_1[i] + ta1_y_xxxzz_xx_0[i] * pa_y[i] - ta1_y_xxxzz_xx_1[i] * pc_y[i];

        ta1_y_xxxyzz_xy_0[i] = ta1_y_xxxy_xy_0[i] * fe_0 - ta1_y_xxxy_xy_1[i] * fe_0 + ta1_y_xxxyz_xy_0[i] * pa_z[i] - ta1_y_xxxyz_xy_1[i] * pc_z[i];

        ta1_y_xxxyzz_xz_0[i] = ta_xxxzz_xz_1[i] + ta1_y_xxxzz_xz_0[i] * pa_y[i] - ta1_y_xxxzz_xz_1[i] * pc_y[i];

        ta1_y_xxxyzz_yy_0[i] =
            2.0 * ta1_y_xyzz_yy_0[i] * fe_0 - 2.0 * ta1_y_xyzz_yy_1[i] * fe_0 + ta1_y_xxyzz_yy_0[i] * pa_x[i] - ta1_y_xxyzz_yy_1[i] * pc_x[i];

        ta1_y_xxxyzz_yz_0[i] =
            2.0 * ta1_y_xyzz_yz_0[i] * fe_0 - 2.0 * ta1_y_xyzz_yz_1[i] * fe_0 + ta1_y_xxyzz_yz_0[i] * pa_x[i] - ta1_y_xxyzz_yz_1[i] * pc_x[i];

        ta1_y_xxxyzz_zz_0[i] = ta_xxxzz_zz_1[i] + ta1_y_xxxzz_zz_0[i] * pa_y[i] - ta1_y_xxxzz_zz_1[i] * pc_y[i];
    }

    // Set up 222-228 components of targeted buffer : ID

    auto ta1_y_xxxzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 222);

    auto ta1_y_xxxzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 223);

    auto ta1_y_xxxzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 224);

    auto ta1_y_xxxzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 225);

    auto ta1_y_xxxzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 226);

    auto ta1_y_xxxzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 227);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_y_xxxz_xx_0,   \
                             ta1_y_xxxz_xx_1,   \
                             ta1_y_xxxz_xy_0,   \
                             ta1_y_xxxz_xy_1,   \
                             ta1_y_xxxzz_xx_0,  \
                             ta1_y_xxxzz_xx_1,  \
                             ta1_y_xxxzz_xy_0,  \
                             ta1_y_xxxzz_xy_1,  \
                             ta1_y_xxxzzz_xx_0, \
                             ta1_y_xxxzzz_xy_0, \
                             ta1_y_xxxzzz_xz_0, \
                             ta1_y_xxxzzz_yy_0, \
                             ta1_y_xxxzzz_yz_0, \
                             ta1_y_xxxzzz_zz_0, \
                             ta1_y_xxzzz_xz_0,  \
                             ta1_y_xxzzz_xz_1,  \
                             ta1_y_xxzzz_yy_0,  \
                             ta1_y_xxzzz_yy_1,  \
                             ta1_y_xxzzz_yz_0,  \
                             ta1_y_xxzzz_yz_1,  \
                             ta1_y_xxzzz_z_0,   \
                             ta1_y_xxzzz_z_1,   \
                             ta1_y_xxzzz_zz_0,  \
                             ta1_y_xxzzz_zz_1,  \
                             ta1_y_xzzz_xz_0,   \
                             ta1_y_xzzz_xz_1,   \
                             ta1_y_xzzz_yy_0,   \
                             ta1_y_xzzz_yy_1,   \
                             ta1_y_xzzz_yz_0,   \
                             ta1_y_xzzz_yz_1,   \
                             ta1_y_xzzz_zz_0,   \
                             ta1_y_xzzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxzzz_xx_0[i] =
            2.0 * ta1_y_xxxz_xx_0[i] * fe_0 - 2.0 * ta1_y_xxxz_xx_1[i] * fe_0 + ta1_y_xxxzz_xx_0[i] * pa_z[i] - ta1_y_xxxzz_xx_1[i] * pc_z[i];

        ta1_y_xxxzzz_xy_0[i] =
            2.0 * ta1_y_xxxz_xy_0[i] * fe_0 - 2.0 * ta1_y_xxxz_xy_1[i] * fe_0 + ta1_y_xxxzz_xy_0[i] * pa_z[i] - ta1_y_xxxzz_xy_1[i] * pc_z[i];

        ta1_y_xxxzzz_xz_0[i] = 2.0 * ta1_y_xzzz_xz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_xz_1[i] * fe_0 + ta1_y_xxzzz_z_0[i] * fe_0 -
                               ta1_y_xxzzz_z_1[i] * fe_0 + ta1_y_xxzzz_xz_0[i] * pa_x[i] - ta1_y_xxzzz_xz_1[i] * pc_x[i];

        ta1_y_xxxzzz_yy_0[i] =
            2.0 * ta1_y_xzzz_yy_0[i] * fe_0 - 2.0 * ta1_y_xzzz_yy_1[i] * fe_0 + ta1_y_xxzzz_yy_0[i] * pa_x[i] - ta1_y_xxzzz_yy_1[i] * pc_x[i];

        ta1_y_xxxzzz_yz_0[i] =
            2.0 * ta1_y_xzzz_yz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_yz_1[i] * fe_0 + ta1_y_xxzzz_yz_0[i] * pa_x[i] - ta1_y_xxzzz_yz_1[i] * pc_x[i];

        ta1_y_xxxzzz_zz_0[i] =
            2.0 * ta1_y_xzzz_zz_0[i] * fe_0 - 2.0 * ta1_y_xzzz_zz_1[i] * fe_0 + ta1_y_xxzzz_zz_0[i] * pa_x[i] - ta1_y_xxzzz_zz_1[i] * pc_x[i];
    }

    // Set up 228-234 components of targeted buffer : ID

    auto ta1_y_xxyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 228);

    auto ta1_y_xxyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 229);

    auto ta1_y_xxyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 230);

    auto ta1_y_xxyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 231);

    auto ta1_y_xxyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 232);

    auto ta1_y_xxyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 233);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_y_xxyy_xx_0,   \
                             ta1_y_xxyy_xx_1,   \
                             ta1_y_xxyy_xz_0,   \
                             ta1_y_xxyy_xz_1,   \
                             ta1_y_xxyyy_xx_0,  \
                             ta1_y_xxyyy_xx_1,  \
                             ta1_y_xxyyy_xz_0,  \
                             ta1_y_xxyyy_xz_1,  \
                             ta1_y_xxyyyy_xx_0, \
                             ta1_y_xxyyyy_xy_0, \
                             ta1_y_xxyyyy_xz_0, \
                             ta1_y_xxyyyy_yy_0, \
                             ta1_y_xxyyyy_yz_0, \
                             ta1_y_xxyyyy_zz_0, \
                             ta1_y_xyyyy_xy_0,  \
                             ta1_y_xyyyy_xy_1,  \
                             ta1_y_xyyyy_y_0,   \
                             ta1_y_xyyyy_y_1,   \
                             ta1_y_xyyyy_yy_0,  \
                             ta1_y_xyyyy_yy_1,  \
                             ta1_y_xyyyy_yz_0,  \
                             ta1_y_xyyyy_yz_1,  \
                             ta1_y_xyyyy_zz_0,  \
                             ta1_y_xyyyy_zz_1,  \
                             ta1_y_yyyy_xy_0,   \
                             ta1_y_yyyy_xy_1,   \
                             ta1_y_yyyy_yy_0,   \
                             ta1_y_yyyy_yy_1,   \
                             ta1_y_yyyy_yz_0,   \
                             ta1_y_yyyy_yz_1,   \
                             ta1_y_yyyy_zz_0,   \
                             ta1_y_yyyy_zz_1,   \
                             ta_xxyyy_xx_1,     \
                             ta_xxyyy_xz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyyy_xx_0[i] = 3.0 * ta1_y_xxyy_xx_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xx_1[i] * fe_0 + ta_xxyyy_xx_1[i] + ta1_y_xxyyy_xx_0[i] * pa_y[i] -
                               ta1_y_xxyyy_xx_1[i] * pc_y[i];

        ta1_y_xxyyyy_xy_0[i] = ta1_y_yyyy_xy_0[i] * fe_0 - ta1_y_yyyy_xy_1[i] * fe_0 + ta1_y_xyyyy_y_0[i] * fe_0 - ta1_y_xyyyy_y_1[i] * fe_0 +
                               ta1_y_xyyyy_xy_0[i] * pa_x[i] - ta1_y_xyyyy_xy_1[i] * pc_x[i];

        ta1_y_xxyyyy_xz_0[i] = 3.0 * ta1_y_xxyy_xz_0[i] * fe_0 - 3.0 * ta1_y_xxyy_xz_1[i] * fe_0 + ta_xxyyy_xz_1[i] + ta1_y_xxyyy_xz_0[i] * pa_y[i] -
                               ta1_y_xxyyy_xz_1[i] * pc_y[i];

        ta1_y_xxyyyy_yy_0[i] = ta1_y_yyyy_yy_0[i] * fe_0 - ta1_y_yyyy_yy_1[i] * fe_0 + ta1_y_xyyyy_yy_0[i] * pa_x[i] - ta1_y_xyyyy_yy_1[i] * pc_x[i];

        ta1_y_xxyyyy_yz_0[i] = ta1_y_yyyy_yz_0[i] * fe_0 - ta1_y_yyyy_yz_1[i] * fe_0 + ta1_y_xyyyy_yz_0[i] * pa_x[i] - ta1_y_xyyyy_yz_1[i] * pc_x[i];

        ta1_y_xxyyyy_zz_0[i] = ta1_y_yyyy_zz_0[i] * fe_0 - ta1_y_yyyy_zz_1[i] * fe_0 + ta1_y_xyyyy_zz_0[i] * pa_x[i] - ta1_y_xyyyy_zz_1[i] * pc_x[i];
    }

    // Set up 234-240 components of targeted buffer : ID

    auto ta1_y_xxyyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 234);

    auto ta1_y_xxyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 235);

    auto ta1_y_xxyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 236);

    auto ta1_y_xxyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 237);

    auto ta1_y_xxyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 238);

    auto ta1_y_xxyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 239);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_y_xxyyy_x_0,   \
                             ta1_y_xxyyy_x_1,   \
                             ta1_y_xxyyy_xx_0,  \
                             ta1_y_xxyyy_xx_1,  \
                             ta1_y_xxyyy_xy_0,  \
                             ta1_y_xxyyy_xy_1,  \
                             ta1_y_xxyyy_xz_0,  \
                             ta1_y_xxyyy_xz_1,  \
                             ta1_y_xxyyy_yy_0,  \
                             ta1_y_xxyyy_yy_1,  \
                             ta1_y_xxyyyz_xx_0, \
                             ta1_y_xxyyyz_xy_0, \
                             ta1_y_xxyyyz_xz_0, \
                             ta1_y_xxyyyz_yy_0, \
                             ta1_y_xxyyyz_yz_0, \
                             ta1_y_xxyyyz_zz_0, \
                             ta1_y_xyyyz_yz_0,  \
                             ta1_y_xyyyz_yz_1,  \
                             ta1_y_xyyyz_zz_0,  \
                             ta1_y_xyyyz_zz_1,  \
                             ta1_y_yyyz_yz_0,   \
                             ta1_y_yyyz_yz_1,   \
                             ta1_y_yyyz_zz_0,   \
                             ta1_y_yyyz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyyz_xx_0[i] = ta1_y_xxyyy_xx_0[i] * pa_z[i] - ta1_y_xxyyy_xx_1[i] * pc_z[i];

        ta1_y_xxyyyz_xy_0[i] = ta1_y_xxyyy_xy_0[i] * pa_z[i] - ta1_y_xxyyy_xy_1[i] * pc_z[i];

        ta1_y_xxyyyz_xz_0[i] = ta1_y_xxyyy_x_0[i] * fe_0 - ta1_y_xxyyy_x_1[i] * fe_0 + ta1_y_xxyyy_xz_0[i] * pa_z[i] - ta1_y_xxyyy_xz_1[i] * pc_z[i];

        ta1_y_xxyyyz_yy_0[i] = ta1_y_xxyyy_yy_0[i] * pa_z[i] - ta1_y_xxyyy_yy_1[i] * pc_z[i];

        ta1_y_xxyyyz_yz_0[i] = ta1_y_yyyz_yz_0[i] * fe_0 - ta1_y_yyyz_yz_1[i] * fe_0 + ta1_y_xyyyz_yz_0[i] * pa_x[i] - ta1_y_xyyyz_yz_1[i] * pc_x[i];

        ta1_y_xxyyyz_zz_0[i] = ta1_y_yyyz_zz_0[i] * fe_0 - ta1_y_yyyz_zz_1[i] * fe_0 + ta1_y_xyyyz_zz_0[i] * pa_x[i] - ta1_y_xyyyz_zz_1[i] * pc_x[i];
    }

    // Set up 240-246 components of targeted buffer : ID

    auto ta1_y_xxyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 240);

    auto ta1_y_xxyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 241);

    auto ta1_y_xxyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 242);

    auto ta1_y_xxyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 243);

    auto ta1_y_xxyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 244);

    auto ta1_y_xxyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 245);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_y_xxyy_xx_0,   \
                             ta1_y_xxyy_xx_1,   \
                             ta1_y_xxyy_xy_0,   \
                             ta1_y_xxyy_xy_1,   \
                             ta1_y_xxyyz_xx_0,  \
                             ta1_y_xxyyz_xx_1,  \
                             ta1_y_xxyyz_xy_0,  \
                             ta1_y_xxyyz_xy_1,  \
                             ta1_y_xxyyzz_xx_0, \
                             ta1_y_xxyyzz_xy_0, \
                             ta1_y_xxyyzz_xz_0, \
                             ta1_y_xxyyzz_yy_0, \
                             ta1_y_xxyyzz_yz_0, \
                             ta1_y_xxyyzz_zz_0, \
                             ta1_y_xxyzz_xz_0,  \
                             ta1_y_xxyzz_xz_1,  \
                             ta1_y_xxzz_xz_0,   \
                             ta1_y_xxzz_xz_1,   \
                             ta1_y_xyyzz_yy_0,  \
                             ta1_y_xyyzz_yy_1,  \
                             ta1_y_xyyzz_yz_0,  \
                             ta1_y_xyyzz_yz_1,  \
                             ta1_y_xyyzz_zz_0,  \
                             ta1_y_xyyzz_zz_1,  \
                             ta1_y_yyzz_yy_0,   \
                             ta1_y_yyzz_yy_1,   \
                             ta1_y_yyzz_yz_0,   \
                             ta1_y_yyzz_yz_1,   \
                             ta1_y_yyzz_zz_0,   \
                             ta1_y_yyzz_zz_1,   \
                             ta_xxyzz_xz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyyzz_xx_0[i] = ta1_y_xxyy_xx_0[i] * fe_0 - ta1_y_xxyy_xx_1[i] * fe_0 + ta1_y_xxyyz_xx_0[i] * pa_z[i] - ta1_y_xxyyz_xx_1[i] * pc_z[i];

        ta1_y_xxyyzz_xy_0[i] = ta1_y_xxyy_xy_0[i] * fe_0 - ta1_y_xxyy_xy_1[i] * fe_0 + ta1_y_xxyyz_xy_0[i] * pa_z[i] - ta1_y_xxyyz_xy_1[i] * pc_z[i];

        ta1_y_xxyyzz_xz_0[i] =
            ta1_y_xxzz_xz_0[i] * fe_0 - ta1_y_xxzz_xz_1[i] * fe_0 + ta_xxyzz_xz_1[i] + ta1_y_xxyzz_xz_0[i] * pa_y[i] - ta1_y_xxyzz_xz_1[i] * pc_y[i];

        ta1_y_xxyyzz_yy_0[i] = ta1_y_yyzz_yy_0[i] * fe_0 - ta1_y_yyzz_yy_1[i] * fe_0 + ta1_y_xyyzz_yy_0[i] * pa_x[i] - ta1_y_xyyzz_yy_1[i] * pc_x[i];

        ta1_y_xxyyzz_yz_0[i] = ta1_y_yyzz_yz_0[i] * fe_0 - ta1_y_yyzz_yz_1[i] * fe_0 + ta1_y_xyyzz_yz_0[i] * pa_x[i] - ta1_y_xyyzz_yz_1[i] * pc_x[i];

        ta1_y_xxyyzz_zz_0[i] = ta1_y_yyzz_zz_0[i] * fe_0 - ta1_y_yyzz_zz_1[i] * fe_0 + ta1_y_xyyzz_zz_0[i] * pa_x[i] - ta1_y_xyyzz_zz_1[i] * pc_x[i];
    }

    // Set up 246-252 components of targeted buffer : ID

    auto ta1_y_xxyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 246);

    auto ta1_y_xxyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 247);

    auto ta1_y_xxyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 248);

    auto ta1_y_xxyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 249);

    auto ta1_y_xxyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 250);

    auto ta1_y_xxyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 251);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_y_xxyz_xy_0,   \
                             ta1_y_xxyz_xy_1,   \
                             ta1_y_xxyzz_xy_0,  \
                             ta1_y_xxyzz_xy_1,  \
                             ta1_y_xxyzzz_xx_0, \
                             ta1_y_xxyzzz_xy_0, \
                             ta1_y_xxyzzz_xz_0, \
                             ta1_y_xxyzzz_yy_0, \
                             ta1_y_xxyzzz_yz_0, \
                             ta1_y_xxyzzz_zz_0, \
                             ta1_y_xxzzz_xx_0,  \
                             ta1_y_xxzzz_xx_1,  \
                             ta1_y_xxzzz_xz_0,  \
                             ta1_y_xxzzz_xz_1,  \
                             ta1_y_xxzzz_zz_0,  \
                             ta1_y_xxzzz_zz_1,  \
                             ta1_y_xyzzz_yy_0,  \
                             ta1_y_xyzzz_yy_1,  \
                             ta1_y_xyzzz_yz_0,  \
                             ta1_y_xyzzz_yz_1,  \
                             ta1_y_yzzz_yy_0,   \
                             ta1_y_yzzz_yy_1,   \
                             ta1_y_yzzz_yz_0,   \
                             ta1_y_yzzz_yz_1,   \
                             ta_xxzzz_xx_1,     \
                             ta_xxzzz_xz_1,     \
                             ta_xxzzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyzzz_xx_0[i] = ta_xxzzz_xx_1[i] + ta1_y_xxzzz_xx_0[i] * pa_y[i] - ta1_y_xxzzz_xx_1[i] * pc_y[i];

        ta1_y_xxyzzz_xy_0[i] =
            2.0 * ta1_y_xxyz_xy_0[i] * fe_0 - 2.0 * ta1_y_xxyz_xy_1[i] * fe_0 + ta1_y_xxyzz_xy_0[i] * pa_z[i] - ta1_y_xxyzz_xy_1[i] * pc_z[i];

        ta1_y_xxyzzz_xz_0[i] = ta_xxzzz_xz_1[i] + ta1_y_xxzzz_xz_0[i] * pa_y[i] - ta1_y_xxzzz_xz_1[i] * pc_y[i];

        ta1_y_xxyzzz_yy_0[i] = ta1_y_yzzz_yy_0[i] * fe_0 - ta1_y_yzzz_yy_1[i] * fe_0 + ta1_y_xyzzz_yy_0[i] * pa_x[i] - ta1_y_xyzzz_yy_1[i] * pc_x[i];

        ta1_y_xxyzzz_yz_0[i] = ta1_y_yzzz_yz_0[i] * fe_0 - ta1_y_yzzz_yz_1[i] * fe_0 + ta1_y_xyzzz_yz_0[i] * pa_x[i] - ta1_y_xyzzz_yz_1[i] * pc_x[i];

        ta1_y_xxyzzz_zz_0[i] = ta_xxzzz_zz_1[i] + ta1_y_xxzzz_zz_0[i] * pa_y[i] - ta1_y_xxzzz_zz_1[i] * pc_y[i];
    }

    // Set up 252-258 components of targeted buffer : ID

    auto ta1_y_xxzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 252);

    auto ta1_y_xxzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 253);

    auto ta1_y_xxzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 254);

    auto ta1_y_xxzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 255);

    auto ta1_y_xxzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 256);

    auto ta1_y_xxzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 257);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_y_xxzz_xx_0,   \
                             ta1_y_xxzz_xx_1,   \
                             ta1_y_xxzz_xy_0,   \
                             ta1_y_xxzz_xy_1,   \
                             ta1_y_xxzzz_xx_0,  \
                             ta1_y_xxzzz_xx_1,  \
                             ta1_y_xxzzz_xy_0,  \
                             ta1_y_xxzzz_xy_1,  \
                             ta1_y_xxzzzz_xx_0, \
                             ta1_y_xxzzzz_xy_0, \
                             ta1_y_xxzzzz_xz_0, \
                             ta1_y_xxzzzz_yy_0, \
                             ta1_y_xxzzzz_yz_0, \
                             ta1_y_xxzzzz_zz_0, \
                             ta1_y_xzzzz_xz_0,  \
                             ta1_y_xzzzz_xz_1,  \
                             ta1_y_xzzzz_yy_0,  \
                             ta1_y_xzzzz_yy_1,  \
                             ta1_y_xzzzz_yz_0,  \
                             ta1_y_xzzzz_yz_1,  \
                             ta1_y_xzzzz_z_0,   \
                             ta1_y_xzzzz_z_1,   \
                             ta1_y_xzzzz_zz_0,  \
                             ta1_y_xzzzz_zz_1,  \
                             ta1_y_zzzz_xz_0,   \
                             ta1_y_zzzz_xz_1,   \
                             ta1_y_zzzz_yy_0,   \
                             ta1_y_zzzz_yy_1,   \
                             ta1_y_zzzz_yz_0,   \
                             ta1_y_zzzz_yz_1,   \
                             ta1_y_zzzz_zz_0,   \
                             ta1_y_zzzz_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxzzzz_xx_0[i] =
            3.0 * ta1_y_xxzz_xx_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xx_1[i] * fe_0 + ta1_y_xxzzz_xx_0[i] * pa_z[i] - ta1_y_xxzzz_xx_1[i] * pc_z[i];

        ta1_y_xxzzzz_xy_0[i] =
            3.0 * ta1_y_xxzz_xy_0[i] * fe_0 - 3.0 * ta1_y_xxzz_xy_1[i] * fe_0 + ta1_y_xxzzz_xy_0[i] * pa_z[i] - ta1_y_xxzzz_xy_1[i] * pc_z[i];

        ta1_y_xxzzzz_xz_0[i] = ta1_y_zzzz_xz_0[i] * fe_0 - ta1_y_zzzz_xz_1[i] * fe_0 + ta1_y_xzzzz_z_0[i] * fe_0 - ta1_y_xzzzz_z_1[i] * fe_0 +
                               ta1_y_xzzzz_xz_0[i] * pa_x[i] - ta1_y_xzzzz_xz_1[i] * pc_x[i];

        ta1_y_xxzzzz_yy_0[i] = ta1_y_zzzz_yy_0[i] * fe_0 - ta1_y_zzzz_yy_1[i] * fe_0 + ta1_y_xzzzz_yy_0[i] * pa_x[i] - ta1_y_xzzzz_yy_1[i] * pc_x[i];

        ta1_y_xxzzzz_yz_0[i] = ta1_y_zzzz_yz_0[i] * fe_0 - ta1_y_zzzz_yz_1[i] * fe_0 + ta1_y_xzzzz_yz_0[i] * pa_x[i] - ta1_y_xzzzz_yz_1[i] * pc_x[i];

        ta1_y_xxzzzz_zz_0[i] = ta1_y_zzzz_zz_0[i] * fe_0 - ta1_y_zzzz_zz_1[i] * fe_0 + ta1_y_xzzzz_zz_0[i] * pa_x[i] - ta1_y_xzzzz_zz_1[i] * pc_x[i];
    }

    // Set up 258-264 components of targeted buffer : ID

    auto ta1_y_xyyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 258);

    auto ta1_y_xyyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 259);

    auto ta1_y_xyyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 260);

    auto ta1_y_xyyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 261);

    auto ta1_y_xyyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 262);

    auto ta1_y_xyyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 263);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_y_xyyyyy_xx_0, \
                             ta1_y_xyyyyy_xy_0, \
                             ta1_y_xyyyyy_xz_0, \
                             ta1_y_xyyyyy_yy_0, \
                             ta1_y_xyyyyy_yz_0, \
                             ta1_y_xyyyyy_zz_0, \
                             ta1_y_yyyyy_x_0,   \
                             ta1_y_yyyyy_x_1,   \
                             ta1_y_yyyyy_xx_0,  \
                             ta1_y_yyyyy_xx_1,  \
                             ta1_y_yyyyy_xy_0,  \
                             ta1_y_yyyyy_xy_1,  \
                             ta1_y_yyyyy_xz_0,  \
                             ta1_y_yyyyy_xz_1,  \
                             ta1_y_yyyyy_y_0,   \
                             ta1_y_yyyyy_y_1,   \
                             ta1_y_yyyyy_yy_0,  \
                             ta1_y_yyyyy_yy_1,  \
                             ta1_y_yyyyy_yz_0,  \
                             ta1_y_yyyyy_yz_1,  \
                             ta1_y_yyyyy_z_0,   \
                             ta1_y_yyyyy_z_1,   \
                             ta1_y_yyyyy_zz_0,  \
                             ta1_y_yyyyy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyyy_xx_0[i] =
            2.0 * ta1_y_yyyyy_x_0[i] * fe_0 - 2.0 * ta1_y_yyyyy_x_1[i] * fe_0 + ta1_y_yyyyy_xx_0[i] * pa_x[i] - ta1_y_yyyyy_xx_1[i] * pc_x[i];

        ta1_y_xyyyyy_xy_0[i] = ta1_y_yyyyy_y_0[i] * fe_0 - ta1_y_yyyyy_y_1[i] * fe_0 + ta1_y_yyyyy_xy_0[i] * pa_x[i] - ta1_y_yyyyy_xy_1[i] * pc_x[i];

        ta1_y_xyyyyy_xz_0[i] = ta1_y_yyyyy_z_0[i] * fe_0 - ta1_y_yyyyy_z_1[i] * fe_0 + ta1_y_yyyyy_xz_0[i] * pa_x[i] - ta1_y_yyyyy_xz_1[i] * pc_x[i];

        ta1_y_xyyyyy_yy_0[i] = ta1_y_yyyyy_yy_0[i] * pa_x[i] - ta1_y_yyyyy_yy_1[i] * pc_x[i];

        ta1_y_xyyyyy_yz_0[i] = ta1_y_yyyyy_yz_0[i] * pa_x[i] - ta1_y_yyyyy_yz_1[i] * pc_x[i];

        ta1_y_xyyyyy_zz_0[i] = ta1_y_yyyyy_zz_0[i] * pa_x[i] - ta1_y_yyyyy_zz_1[i] * pc_x[i];
    }

    // Set up 264-270 components of targeted buffer : ID

    auto ta1_y_xyyyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 264);

    auto ta1_y_xyyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 265);

    auto ta1_y_xyyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 266);

    auto ta1_y_xyyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 267);

    auto ta1_y_xyyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 268);

    auto ta1_y_xyyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 269);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_y_xyyyy_xx_0,  \
                             ta1_y_xyyyy_xx_1,  \
                             ta1_y_xyyyy_xy_0,  \
                             ta1_y_xyyyy_xy_1,  \
                             ta1_y_xyyyyz_xx_0, \
                             ta1_y_xyyyyz_xy_0, \
                             ta1_y_xyyyyz_xz_0, \
                             ta1_y_xyyyyz_yy_0, \
                             ta1_y_xyyyyz_yz_0, \
                             ta1_y_xyyyyz_zz_0, \
                             ta1_y_yyyyz_xz_0,  \
                             ta1_y_yyyyz_xz_1,  \
                             ta1_y_yyyyz_yy_0,  \
                             ta1_y_yyyyz_yy_1,  \
                             ta1_y_yyyyz_yz_0,  \
                             ta1_y_yyyyz_yz_1,  \
                             ta1_y_yyyyz_z_0,   \
                             ta1_y_yyyyz_z_1,   \
                             ta1_y_yyyyz_zz_0,  \
                             ta1_y_yyyyz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyyz_xx_0[i] = ta1_y_xyyyy_xx_0[i] * pa_z[i] - ta1_y_xyyyy_xx_1[i] * pc_z[i];

        ta1_y_xyyyyz_xy_0[i] = ta1_y_xyyyy_xy_0[i] * pa_z[i] - ta1_y_xyyyy_xy_1[i] * pc_z[i];

        ta1_y_xyyyyz_xz_0[i] = ta1_y_yyyyz_z_0[i] * fe_0 - ta1_y_yyyyz_z_1[i] * fe_0 + ta1_y_yyyyz_xz_0[i] * pa_x[i] - ta1_y_yyyyz_xz_1[i] * pc_x[i];

        ta1_y_xyyyyz_yy_0[i] = ta1_y_yyyyz_yy_0[i] * pa_x[i] - ta1_y_yyyyz_yy_1[i] * pc_x[i];

        ta1_y_xyyyyz_yz_0[i] = ta1_y_yyyyz_yz_0[i] * pa_x[i] - ta1_y_yyyyz_yz_1[i] * pc_x[i];

        ta1_y_xyyyyz_zz_0[i] = ta1_y_yyyyz_zz_0[i] * pa_x[i] - ta1_y_yyyyz_zz_1[i] * pc_x[i];
    }

    // Set up 270-276 components of targeted buffer : ID

    auto ta1_y_xyyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 270);

    auto ta1_y_xyyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 271);

    auto ta1_y_xyyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 272);

    auto ta1_y_xyyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 273);

    auto ta1_y_xyyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 274);

    auto ta1_y_xyyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 275);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_y_xyyyzz_xx_0, \
                             ta1_y_xyyyzz_xy_0, \
                             ta1_y_xyyyzz_xz_0, \
                             ta1_y_xyyyzz_yy_0, \
                             ta1_y_xyyyzz_yz_0, \
                             ta1_y_xyyyzz_zz_0, \
                             ta1_y_yyyzz_x_0,   \
                             ta1_y_yyyzz_x_1,   \
                             ta1_y_yyyzz_xx_0,  \
                             ta1_y_yyyzz_xx_1,  \
                             ta1_y_yyyzz_xy_0,  \
                             ta1_y_yyyzz_xy_1,  \
                             ta1_y_yyyzz_xz_0,  \
                             ta1_y_yyyzz_xz_1,  \
                             ta1_y_yyyzz_y_0,   \
                             ta1_y_yyyzz_y_1,   \
                             ta1_y_yyyzz_yy_0,  \
                             ta1_y_yyyzz_yy_1,  \
                             ta1_y_yyyzz_yz_0,  \
                             ta1_y_yyyzz_yz_1,  \
                             ta1_y_yyyzz_z_0,   \
                             ta1_y_yyyzz_z_1,   \
                             ta1_y_yyyzz_zz_0,  \
                             ta1_y_yyyzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyyzz_xx_0[i] =
            2.0 * ta1_y_yyyzz_x_0[i] * fe_0 - 2.0 * ta1_y_yyyzz_x_1[i] * fe_0 + ta1_y_yyyzz_xx_0[i] * pa_x[i] - ta1_y_yyyzz_xx_1[i] * pc_x[i];

        ta1_y_xyyyzz_xy_0[i] = ta1_y_yyyzz_y_0[i] * fe_0 - ta1_y_yyyzz_y_1[i] * fe_0 + ta1_y_yyyzz_xy_0[i] * pa_x[i] - ta1_y_yyyzz_xy_1[i] * pc_x[i];

        ta1_y_xyyyzz_xz_0[i] = ta1_y_yyyzz_z_0[i] * fe_0 - ta1_y_yyyzz_z_1[i] * fe_0 + ta1_y_yyyzz_xz_0[i] * pa_x[i] - ta1_y_yyyzz_xz_1[i] * pc_x[i];

        ta1_y_xyyyzz_yy_0[i] = ta1_y_yyyzz_yy_0[i] * pa_x[i] - ta1_y_yyyzz_yy_1[i] * pc_x[i];

        ta1_y_xyyyzz_yz_0[i] = ta1_y_yyyzz_yz_0[i] * pa_x[i] - ta1_y_yyyzz_yz_1[i] * pc_x[i];

        ta1_y_xyyyzz_zz_0[i] = ta1_y_yyyzz_zz_0[i] * pa_x[i] - ta1_y_yyyzz_zz_1[i] * pc_x[i];
    }

    // Set up 276-282 components of targeted buffer : ID

    auto ta1_y_xyyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 276);

    auto ta1_y_xyyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 277);

    auto ta1_y_xyyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 278);

    auto ta1_y_xyyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 279);

    auto ta1_y_xyyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 280);

    auto ta1_y_xyyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 281);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_y_xyyzzz_xx_0, \
                             ta1_y_xyyzzz_xy_0, \
                             ta1_y_xyyzzz_xz_0, \
                             ta1_y_xyyzzz_yy_0, \
                             ta1_y_xyyzzz_yz_0, \
                             ta1_y_xyyzzz_zz_0, \
                             ta1_y_yyzzz_x_0,   \
                             ta1_y_yyzzz_x_1,   \
                             ta1_y_yyzzz_xx_0,  \
                             ta1_y_yyzzz_xx_1,  \
                             ta1_y_yyzzz_xy_0,  \
                             ta1_y_yyzzz_xy_1,  \
                             ta1_y_yyzzz_xz_0,  \
                             ta1_y_yyzzz_xz_1,  \
                             ta1_y_yyzzz_y_0,   \
                             ta1_y_yyzzz_y_1,   \
                             ta1_y_yyzzz_yy_0,  \
                             ta1_y_yyzzz_yy_1,  \
                             ta1_y_yyzzz_yz_0,  \
                             ta1_y_yyzzz_yz_1,  \
                             ta1_y_yyzzz_z_0,   \
                             ta1_y_yyzzz_z_1,   \
                             ta1_y_yyzzz_zz_0,  \
                             ta1_y_yyzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyzzz_xx_0[i] =
            2.0 * ta1_y_yyzzz_x_0[i] * fe_0 - 2.0 * ta1_y_yyzzz_x_1[i] * fe_0 + ta1_y_yyzzz_xx_0[i] * pa_x[i] - ta1_y_yyzzz_xx_1[i] * pc_x[i];

        ta1_y_xyyzzz_xy_0[i] = ta1_y_yyzzz_y_0[i] * fe_0 - ta1_y_yyzzz_y_1[i] * fe_0 + ta1_y_yyzzz_xy_0[i] * pa_x[i] - ta1_y_yyzzz_xy_1[i] * pc_x[i];

        ta1_y_xyyzzz_xz_0[i] = ta1_y_yyzzz_z_0[i] * fe_0 - ta1_y_yyzzz_z_1[i] * fe_0 + ta1_y_yyzzz_xz_0[i] * pa_x[i] - ta1_y_yyzzz_xz_1[i] * pc_x[i];

        ta1_y_xyyzzz_yy_0[i] = ta1_y_yyzzz_yy_0[i] * pa_x[i] - ta1_y_yyzzz_yy_1[i] * pc_x[i];

        ta1_y_xyyzzz_yz_0[i] = ta1_y_yyzzz_yz_0[i] * pa_x[i] - ta1_y_yyzzz_yz_1[i] * pc_x[i];

        ta1_y_xyyzzz_zz_0[i] = ta1_y_yyzzz_zz_0[i] * pa_x[i] - ta1_y_yyzzz_zz_1[i] * pc_x[i];
    }

    // Set up 282-288 components of targeted buffer : ID

    auto ta1_y_xyzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 282);

    auto ta1_y_xyzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 283);

    auto ta1_y_xyzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 284);

    auto ta1_y_xyzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 285);

    auto ta1_y_xyzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 286);

    auto ta1_y_xyzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 287);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_y_xyzzzz_xx_0, \
                             ta1_y_xyzzzz_xy_0, \
                             ta1_y_xyzzzz_xz_0, \
                             ta1_y_xyzzzz_yy_0, \
                             ta1_y_xyzzzz_yz_0, \
                             ta1_y_xyzzzz_zz_0, \
                             ta1_y_xzzzz_xx_0,  \
                             ta1_y_xzzzz_xx_1,  \
                             ta1_y_xzzzz_xz_0,  \
                             ta1_y_xzzzz_xz_1,  \
                             ta1_y_yzzzz_xy_0,  \
                             ta1_y_yzzzz_xy_1,  \
                             ta1_y_yzzzz_y_0,   \
                             ta1_y_yzzzz_y_1,   \
                             ta1_y_yzzzz_yy_0,  \
                             ta1_y_yzzzz_yy_1,  \
                             ta1_y_yzzzz_yz_0,  \
                             ta1_y_yzzzz_yz_1,  \
                             ta1_y_yzzzz_zz_0,  \
                             ta1_y_yzzzz_zz_1,  \
                             ta_xzzzz_xx_1,     \
                             ta_xzzzz_xz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyzzzz_xx_0[i] = ta_xzzzz_xx_1[i] + ta1_y_xzzzz_xx_0[i] * pa_y[i] - ta1_y_xzzzz_xx_1[i] * pc_y[i];

        ta1_y_xyzzzz_xy_0[i] = ta1_y_yzzzz_y_0[i] * fe_0 - ta1_y_yzzzz_y_1[i] * fe_0 + ta1_y_yzzzz_xy_0[i] * pa_x[i] - ta1_y_yzzzz_xy_1[i] * pc_x[i];

        ta1_y_xyzzzz_xz_0[i] = ta_xzzzz_xz_1[i] + ta1_y_xzzzz_xz_0[i] * pa_y[i] - ta1_y_xzzzz_xz_1[i] * pc_y[i];

        ta1_y_xyzzzz_yy_0[i] = ta1_y_yzzzz_yy_0[i] * pa_x[i] - ta1_y_yzzzz_yy_1[i] * pc_x[i];

        ta1_y_xyzzzz_yz_0[i] = ta1_y_yzzzz_yz_0[i] * pa_x[i] - ta1_y_yzzzz_yz_1[i] * pc_x[i];

        ta1_y_xyzzzz_zz_0[i] = ta1_y_yzzzz_zz_0[i] * pa_x[i] - ta1_y_yzzzz_zz_1[i] * pc_x[i];
    }

    // Set up 288-294 components of targeted buffer : ID

    auto ta1_y_xzzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 288);

    auto ta1_y_xzzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 289);

    auto ta1_y_xzzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 290);

    auto ta1_y_xzzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 291);

    auto ta1_y_xzzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 292);

    auto ta1_y_xzzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 293);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_y_xzzzzz_xx_0, \
                             ta1_y_xzzzzz_xy_0, \
                             ta1_y_xzzzzz_xz_0, \
                             ta1_y_xzzzzz_yy_0, \
                             ta1_y_xzzzzz_yz_0, \
                             ta1_y_xzzzzz_zz_0, \
                             ta1_y_zzzzz_x_0,   \
                             ta1_y_zzzzz_x_1,   \
                             ta1_y_zzzzz_xx_0,  \
                             ta1_y_zzzzz_xx_1,  \
                             ta1_y_zzzzz_xy_0,  \
                             ta1_y_zzzzz_xy_1,  \
                             ta1_y_zzzzz_xz_0,  \
                             ta1_y_zzzzz_xz_1,  \
                             ta1_y_zzzzz_y_0,   \
                             ta1_y_zzzzz_y_1,   \
                             ta1_y_zzzzz_yy_0,  \
                             ta1_y_zzzzz_yy_1,  \
                             ta1_y_zzzzz_yz_0,  \
                             ta1_y_zzzzz_yz_1,  \
                             ta1_y_zzzzz_z_0,   \
                             ta1_y_zzzzz_z_1,   \
                             ta1_y_zzzzz_zz_0,  \
                             ta1_y_zzzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzzzzz_xx_0[i] =
            2.0 * ta1_y_zzzzz_x_0[i] * fe_0 - 2.0 * ta1_y_zzzzz_x_1[i] * fe_0 + ta1_y_zzzzz_xx_0[i] * pa_x[i] - ta1_y_zzzzz_xx_1[i] * pc_x[i];

        ta1_y_xzzzzz_xy_0[i] = ta1_y_zzzzz_y_0[i] * fe_0 - ta1_y_zzzzz_y_1[i] * fe_0 + ta1_y_zzzzz_xy_0[i] * pa_x[i] - ta1_y_zzzzz_xy_1[i] * pc_x[i];

        ta1_y_xzzzzz_xz_0[i] = ta1_y_zzzzz_z_0[i] * fe_0 - ta1_y_zzzzz_z_1[i] * fe_0 + ta1_y_zzzzz_xz_0[i] * pa_x[i] - ta1_y_zzzzz_xz_1[i] * pc_x[i];

        ta1_y_xzzzzz_yy_0[i] = ta1_y_zzzzz_yy_0[i] * pa_x[i] - ta1_y_zzzzz_yy_1[i] * pc_x[i];

        ta1_y_xzzzzz_yz_0[i] = ta1_y_zzzzz_yz_0[i] * pa_x[i] - ta1_y_zzzzz_yz_1[i] * pc_x[i];

        ta1_y_xzzzzz_zz_0[i] = ta1_y_zzzzz_zz_0[i] * pa_x[i] - ta1_y_zzzzz_zz_1[i] * pc_x[i];
    }

    // Set up 294-300 components of targeted buffer : ID

    auto ta1_y_yyyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 294);

    auto ta1_y_yyyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 295);

    auto ta1_y_yyyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 296);

    auto ta1_y_yyyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 297);

    auto ta1_y_yyyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 298);

    auto ta1_y_yyyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 299);

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_y_yyyy_xx_0,   \
                             ta1_y_yyyy_xx_1,   \
                             ta1_y_yyyy_xy_0,   \
                             ta1_y_yyyy_xy_1,   \
                             ta1_y_yyyy_xz_0,   \
                             ta1_y_yyyy_xz_1,   \
                             ta1_y_yyyy_yy_0,   \
                             ta1_y_yyyy_yy_1,   \
                             ta1_y_yyyy_yz_0,   \
                             ta1_y_yyyy_yz_1,   \
                             ta1_y_yyyy_zz_0,   \
                             ta1_y_yyyy_zz_1,   \
                             ta1_y_yyyyy_x_0,   \
                             ta1_y_yyyyy_x_1,   \
                             ta1_y_yyyyy_xx_0,  \
                             ta1_y_yyyyy_xx_1,  \
                             ta1_y_yyyyy_xy_0,  \
                             ta1_y_yyyyy_xy_1,  \
                             ta1_y_yyyyy_xz_0,  \
                             ta1_y_yyyyy_xz_1,  \
                             ta1_y_yyyyy_y_0,   \
                             ta1_y_yyyyy_y_1,   \
                             ta1_y_yyyyy_yy_0,  \
                             ta1_y_yyyyy_yy_1,  \
                             ta1_y_yyyyy_yz_0,  \
                             ta1_y_yyyyy_yz_1,  \
                             ta1_y_yyyyy_z_0,   \
                             ta1_y_yyyyy_z_1,   \
                             ta1_y_yyyyy_zz_0,  \
                             ta1_y_yyyyy_zz_1,  \
                             ta1_y_yyyyyy_xx_0, \
                             ta1_y_yyyyyy_xy_0, \
                             ta1_y_yyyyyy_xz_0, \
                             ta1_y_yyyyyy_yy_0, \
                             ta1_y_yyyyyy_yz_0, \
                             ta1_y_yyyyyy_zz_0, \
                             ta_yyyyy_xx_1,     \
                             ta_yyyyy_xy_1,     \
                             ta_yyyyy_xz_1,     \
                             ta_yyyyy_yy_1,     \
                             ta_yyyyy_yz_1,     \
                             ta_yyyyy_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyyy_xx_0[i] = 5.0 * ta1_y_yyyy_xx_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xx_1[i] * fe_0 + ta_yyyyy_xx_1[i] + ta1_y_yyyyy_xx_0[i] * pa_y[i] -
                               ta1_y_yyyyy_xx_1[i] * pc_y[i];

        ta1_y_yyyyyy_xy_0[i] = 5.0 * ta1_y_yyyy_xy_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xy_1[i] * fe_0 + ta1_y_yyyyy_x_0[i] * fe_0 -
                               ta1_y_yyyyy_x_1[i] * fe_0 + ta_yyyyy_xy_1[i] + ta1_y_yyyyy_xy_0[i] * pa_y[i] - ta1_y_yyyyy_xy_1[i] * pc_y[i];

        ta1_y_yyyyyy_xz_0[i] = 5.0 * ta1_y_yyyy_xz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_xz_1[i] * fe_0 + ta_yyyyy_xz_1[i] + ta1_y_yyyyy_xz_0[i] * pa_y[i] -
                               ta1_y_yyyyy_xz_1[i] * pc_y[i];

        ta1_y_yyyyyy_yy_0[i] = 5.0 * ta1_y_yyyy_yy_0[i] * fe_0 - 5.0 * ta1_y_yyyy_yy_1[i] * fe_0 + 2.0 * ta1_y_yyyyy_y_0[i] * fe_0 -
                               2.0 * ta1_y_yyyyy_y_1[i] * fe_0 + ta_yyyyy_yy_1[i] + ta1_y_yyyyy_yy_0[i] * pa_y[i] - ta1_y_yyyyy_yy_1[i] * pc_y[i];

        ta1_y_yyyyyy_yz_0[i] = 5.0 * ta1_y_yyyy_yz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_yz_1[i] * fe_0 + ta1_y_yyyyy_z_0[i] * fe_0 -
                               ta1_y_yyyyy_z_1[i] * fe_0 + ta_yyyyy_yz_1[i] + ta1_y_yyyyy_yz_0[i] * pa_y[i] - ta1_y_yyyyy_yz_1[i] * pc_y[i];

        ta1_y_yyyyyy_zz_0[i] = 5.0 * ta1_y_yyyy_zz_0[i] * fe_0 - 5.0 * ta1_y_yyyy_zz_1[i] * fe_0 + ta_yyyyy_zz_1[i] + ta1_y_yyyyy_zz_0[i] * pa_y[i] -
                               ta1_y_yyyyy_zz_1[i] * pc_y[i];
    }

    // Set up 300-306 components of targeted buffer : ID

    auto ta1_y_yyyyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 300);

    auto ta1_y_yyyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 301);

    auto ta1_y_yyyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 302);

    auto ta1_y_yyyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 303);

    auto ta1_y_yyyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 304);

    auto ta1_y_yyyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 305);

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta1_y_yyyyy_x_0,   \
                             ta1_y_yyyyy_x_1,   \
                             ta1_y_yyyyy_xx_0,  \
                             ta1_y_yyyyy_xx_1,  \
                             ta1_y_yyyyy_xy_0,  \
                             ta1_y_yyyyy_xy_1,  \
                             ta1_y_yyyyy_xz_0,  \
                             ta1_y_yyyyy_xz_1,  \
                             ta1_y_yyyyy_y_0,   \
                             ta1_y_yyyyy_y_1,   \
                             ta1_y_yyyyy_yy_0,  \
                             ta1_y_yyyyy_yy_1,  \
                             ta1_y_yyyyy_yz_0,  \
                             ta1_y_yyyyy_yz_1,  \
                             ta1_y_yyyyy_z_0,   \
                             ta1_y_yyyyy_z_1,   \
                             ta1_y_yyyyy_zz_0,  \
                             ta1_y_yyyyy_zz_1,  \
                             ta1_y_yyyyyz_xx_0, \
                             ta1_y_yyyyyz_xy_0, \
                             ta1_y_yyyyyz_xz_0, \
                             ta1_y_yyyyyz_yy_0, \
                             ta1_y_yyyyyz_yz_0, \
                             ta1_y_yyyyyz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyyz_xx_0[i] = ta1_y_yyyyy_xx_0[i] * pa_z[i] - ta1_y_yyyyy_xx_1[i] * pc_z[i];

        ta1_y_yyyyyz_xy_0[i] = ta1_y_yyyyy_xy_0[i] * pa_z[i] - ta1_y_yyyyy_xy_1[i] * pc_z[i];

        ta1_y_yyyyyz_xz_0[i] = ta1_y_yyyyy_x_0[i] * fe_0 - ta1_y_yyyyy_x_1[i] * fe_0 + ta1_y_yyyyy_xz_0[i] * pa_z[i] - ta1_y_yyyyy_xz_1[i] * pc_z[i];

        ta1_y_yyyyyz_yy_0[i] = ta1_y_yyyyy_yy_0[i] * pa_z[i] - ta1_y_yyyyy_yy_1[i] * pc_z[i];

        ta1_y_yyyyyz_yz_0[i] = ta1_y_yyyyy_y_0[i] * fe_0 - ta1_y_yyyyy_y_1[i] * fe_0 + ta1_y_yyyyy_yz_0[i] * pa_z[i] - ta1_y_yyyyy_yz_1[i] * pc_z[i];

        ta1_y_yyyyyz_zz_0[i] =
            2.0 * ta1_y_yyyyy_z_0[i] * fe_0 - 2.0 * ta1_y_yyyyy_z_1[i] * fe_0 + ta1_y_yyyyy_zz_0[i] * pa_z[i] - ta1_y_yyyyy_zz_1[i] * pc_z[i];
    }

    // Set up 306-312 components of targeted buffer : ID

    auto ta1_y_yyyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 306);

    auto ta1_y_yyyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 307);

    auto ta1_y_yyyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 308);

    auto ta1_y_yyyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 309);

    auto ta1_y_yyyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 310);

    auto ta1_y_yyyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 311);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_y_yyyy_xx_0,   \
                             ta1_y_yyyy_xx_1,   \
                             ta1_y_yyyy_xy_0,   \
                             ta1_y_yyyy_xy_1,   \
                             ta1_y_yyyy_yy_0,   \
                             ta1_y_yyyy_yy_1,   \
                             ta1_y_yyyy_yz_0,   \
                             ta1_y_yyyy_yz_1,   \
                             ta1_y_yyyyz_xx_0,  \
                             ta1_y_yyyyz_xx_1,  \
                             ta1_y_yyyyz_xy_0,  \
                             ta1_y_yyyyz_xy_1,  \
                             ta1_y_yyyyz_y_0,   \
                             ta1_y_yyyyz_y_1,   \
                             ta1_y_yyyyz_yy_0,  \
                             ta1_y_yyyyz_yy_1,  \
                             ta1_y_yyyyz_yz_0,  \
                             ta1_y_yyyyz_yz_1,  \
                             ta1_y_yyyyzz_xx_0, \
                             ta1_y_yyyyzz_xy_0, \
                             ta1_y_yyyyzz_xz_0, \
                             ta1_y_yyyyzz_yy_0, \
                             ta1_y_yyyyzz_yz_0, \
                             ta1_y_yyyyzz_zz_0, \
                             ta1_y_yyyzz_xz_0,  \
                             ta1_y_yyyzz_xz_1,  \
                             ta1_y_yyyzz_zz_0,  \
                             ta1_y_yyyzz_zz_1,  \
                             ta1_y_yyzz_xz_0,   \
                             ta1_y_yyzz_xz_1,   \
                             ta1_y_yyzz_zz_0,   \
                             ta1_y_yyzz_zz_1,   \
                             ta_yyyzz_xz_1,     \
                             ta_yyyzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyyzz_xx_0[i] = ta1_y_yyyy_xx_0[i] * fe_0 - ta1_y_yyyy_xx_1[i] * fe_0 + ta1_y_yyyyz_xx_0[i] * pa_z[i] - ta1_y_yyyyz_xx_1[i] * pc_z[i];

        ta1_y_yyyyzz_xy_0[i] = ta1_y_yyyy_xy_0[i] * fe_0 - ta1_y_yyyy_xy_1[i] * fe_0 + ta1_y_yyyyz_xy_0[i] * pa_z[i] - ta1_y_yyyyz_xy_1[i] * pc_z[i];

        ta1_y_yyyyzz_xz_0[i] = 3.0 * ta1_y_yyzz_xz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xz_1[i] * fe_0 + ta_yyyzz_xz_1[i] + ta1_y_yyyzz_xz_0[i] * pa_y[i] -
                               ta1_y_yyyzz_xz_1[i] * pc_y[i];

        ta1_y_yyyyzz_yy_0[i] = ta1_y_yyyy_yy_0[i] * fe_0 - ta1_y_yyyy_yy_1[i] * fe_0 + ta1_y_yyyyz_yy_0[i] * pa_z[i] - ta1_y_yyyyz_yy_1[i] * pc_z[i];

        ta1_y_yyyyzz_yz_0[i] = ta1_y_yyyy_yz_0[i] * fe_0 - ta1_y_yyyy_yz_1[i] * fe_0 + ta1_y_yyyyz_y_0[i] * fe_0 - ta1_y_yyyyz_y_1[i] * fe_0 +
                               ta1_y_yyyyz_yz_0[i] * pa_z[i] - ta1_y_yyyyz_yz_1[i] * pc_z[i];

        ta1_y_yyyyzz_zz_0[i] = 3.0 * ta1_y_yyzz_zz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_zz_1[i] * fe_0 + ta_yyyzz_zz_1[i] + ta1_y_yyyzz_zz_0[i] * pa_y[i] -
                               ta1_y_yyyzz_zz_1[i] * pc_y[i];
    }

    // Set up 312-318 components of targeted buffer : ID

    auto ta1_y_yyyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 312);

    auto ta1_y_yyyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 313);

    auto ta1_y_yyyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 314);

    auto ta1_y_yyyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 315);

    auto ta1_y_yyyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 316);

    auto ta1_y_yyyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 317);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_y_yyyz_xx_0,   \
                             ta1_y_yyyz_xx_1,   \
                             ta1_y_yyyz_xy_0,   \
                             ta1_y_yyyz_xy_1,   \
                             ta1_y_yyyz_yy_0,   \
                             ta1_y_yyyz_yy_1,   \
                             ta1_y_yyyz_yz_0,   \
                             ta1_y_yyyz_yz_1,   \
                             ta1_y_yyyzz_xx_0,  \
                             ta1_y_yyyzz_xx_1,  \
                             ta1_y_yyyzz_xy_0,  \
                             ta1_y_yyyzz_xy_1,  \
                             ta1_y_yyyzz_y_0,   \
                             ta1_y_yyyzz_y_1,   \
                             ta1_y_yyyzz_yy_0,  \
                             ta1_y_yyyzz_yy_1,  \
                             ta1_y_yyyzz_yz_0,  \
                             ta1_y_yyyzz_yz_1,  \
                             ta1_y_yyyzzz_xx_0, \
                             ta1_y_yyyzzz_xy_0, \
                             ta1_y_yyyzzz_xz_0, \
                             ta1_y_yyyzzz_yy_0, \
                             ta1_y_yyyzzz_yz_0, \
                             ta1_y_yyyzzz_zz_0, \
                             ta1_y_yyzzz_xz_0,  \
                             ta1_y_yyzzz_xz_1,  \
                             ta1_y_yyzzz_zz_0,  \
                             ta1_y_yyzzz_zz_1,  \
                             ta1_y_yzzz_xz_0,   \
                             ta1_y_yzzz_xz_1,   \
                             ta1_y_yzzz_zz_0,   \
                             ta1_y_yzzz_zz_1,   \
                             ta_yyzzz_xz_1,     \
                             ta_yyzzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyzzz_xx_0[i] =
            2.0 * ta1_y_yyyz_xx_0[i] * fe_0 - 2.0 * ta1_y_yyyz_xx_1[i] * fe_0 + ta1_y_yyyzz_xx_0[i] * pa_z[i] - ta1_y_yyyzz_xx_1[i] * pc_z[i];

        ta1_y_yyyzzz_xy_0[i] =
            2.0 * ta1_y_yyyz_xy_0[i] * fe_0 - 2.0 * ta1_y_yyyz_xy_1[i] * fe_0 + ta1_y_yyyzz_xy_0[i] * pa_z[i] - ta1_y_yyyzz_xy_1[i] * pc_z[i];

        ta1_y_yyyzzz_xz_0[i] = 2.0 * ta1_y_yzzz_xz_0[i] * fe_0 - 2.0 * ta1_y_yzzz_xz_1[i] * fe_0 + ta_yyzzz_xz_1[i] + ta1_y_yyzzz_xz_0[i] * pa_y[i] -
                               ta1_y_yyzzz_xz_1[i] * pc_y[i];

        ta1_y_yyyzzz_yy_0[i] =
            2.0 * ta1_y_yyyz_yy_0[i] * fe_0 - 2.0 * ta1_y_yyyz_yy_1[i] * fe_0 + ta1_y_yyyzz_yy_0[i] * pa_z[i] - ta1_y_yyyzz_yy_1[i] * pc_z[i];

        ta1_y_yyyzzz_yz_0[i] = 2.0 * ta1_y_yyyz_yz_0[i] * fe_0 - 2.0 * ta1_y_yyyz_yz_1[i] * fe_0 + ta1_y_yyyzz_y_0[i] * fe_0 -
                               ta1_y_yyyzz_y_1[i] * fe_0 + ta1_y_yyyzz_yz_0[i] * pa_z[i] - ta1_y_yyyzz_yz_1[i] * pc_z[i];

        ta1_y_yyyzzz_zz_0[i] = 2.0 * ta1_y_yzzz_zz_0[i] * fe_0 - 2.0 * ta1_y_yzzz_zz_1[i] * fe_0 + ta_yyzzz_zz_1[i] + ta1_y_yyzzz_zz_0[i] * pa_y[i] -
                               ta1_y_yyzzz_zz_1[i] * pc_y[i];
    }

    // Set up 318-324 components of targeted buffer : ID

    auto ta1_y_yyzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 318);

    auto ta1_y_yyzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 319);

    auto ta1_y_yyzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 320);

    auto ta1_y_yyzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 321);

    auto ta1_y_yyzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 322);

    auto ta1_y_yyzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 323);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_y_yyzz_xx_0,   \
                             ta1_y_yyzz_xx_1,   \
                             ta1_y_yyzz_xy_0,   \
                             ta1_y_yyzz_xy_1,   \
                             ta1_y_yyzz_yy_0,   \
                             ta1_y_yyzz_yy_1,   \
                             ta1_y_yyzz_yz_0,   \
                             ta1_y_yyzz_yz_1,   \
                             ta1_y_yyzzz_xx_0,  \
                             ta1_y_yyzzz_xx_1,  \
                             ta1_y_yyzzz_xy_0,  \
                             ta1_y_yyzzz_xy_1,  \
                             ta1_y_yyzzz_y_0,   \
                             ta1_y_yyzzz_y_1,   \
                             ta1_y_yyzzz_yy_0,  \
                             ta1_y_yyzzz_yy_1,  \
                             ta1_y_yyzzz_yz_0,  \
                             ta1_y_yyzzz_yz_1,  \
                             ta1_y_yyzzzz_xx_0, \
                             ta1_y_yyzzzz_xy_0, \
                             ta1_y_yyzzzz_xz_0, \
                             ta1_y_yyzzzz_yy_0, \
                             ta1_y_yyzzzz_yz_0, \
                             ta1_y_yyzzzz_zz_0, \
                             ta1_y_yzzzz_xz_0,  \
                             ta1_y_yzzzz_xz_1,  \
                             ta1_y_yzzzz_zz_0,  \
                             ta1_y_yzzzz_zz_1,  \
                             ta1_y_zzzz_xz_0,   \
                             ta1_y_zzzz_xz_1,   \
                             ta1_y_zzzz_zz_0,   \
                             ta1_y_zzzz_zz_1,   \
                             ta_yzzzz_xz_1,     \
                             ta_yzzzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyzzzz_xx_0[i] =
            3.0 * ta1_y_yyzz_xx_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xx_1[i] * fe_0 + ta1_y_yyzzz_xx_0[i] * pa_z[i] - ta1_y_yyzzz_xx_1[i] * pc_z[i];

        ta1_y_yyzzzz_xy_0[i] =
            3.0 * ta1_y_yyzz_xy_0[i] * fe_0 - 3.0 * ta1_y_yyzz_xy_1[i] * fe_0 + ta1_y_yyzzz_xy_0[i] * pa_z[i] - ta1_y_yyzzz_xy_1[i] * pc_z[i];

        ta1_y_yyzzzz_xz_0[i] =
            ta1_y_zzzz_xz_0[i] * fe_0 - ta1_y_zzzz_xz_1[i] * fe_0 + ta_yzzzz_xz_1[i] + ta1_y_yzzzz_xz_0[i] * pa_y[i] - ta1_y_yzzzz_xz_1[i] * pc_y[i];

        ta1_y_yyzzzz_yy_0[i] =
            3.0 * ta1_y_yyzz_yy_0[i] * fe_0 - 3.0 * ta1_y_yyzz_yy_1[i] * fe_0 + ta1_y_yyzzz_yy_0[i] * pa_z[i] - ta1_y_yyzzz_yy_1[i] * pc_z[i];

        ta1_y_yyzzzz_yz_0[i] = 3.0 * ta1_y_yyzz_yz_0[i] * fe_0 - 3.0 * ta1_y_yyzz_yz_1[i] * fe_0 + ta1_y_yyzzz_y_0[i] * fe_0 -
                               ta1_y_yyzzz_y_1[i] * fe_0 + ta1_y_yyzzz_yz_0[i] * pa_z[i] - ta1_y_yyzzz_yz_1[i] * pc_z[i];

        ta1_y_yyzzzz_zz_0[i] =
            ta1_y_zzzz_zz_0[i] * fe_0 - ta1_y_zzzz_zz_1[i] * fe_0 + ta_yzzzz_zz_1[i] + ta1_y_yzzzz_zz_0[i] * pa_y[i] - ta1_y_yzzzz_zz_1[i] * pc_y[i];
    }

    // Set up 324-330 components of targeted buffer : ID

    auto ta1_y_yzzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 324);

    auto ta1_y_yzzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 325);

    auto ta1_y_yzzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 326);

    auto ta1_y_yzzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 327);

    auto ta1_y_yzzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 328);

    auto ta1_y_yzzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 329);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_y_yzzz_xy_0,   \
                             ta1_y_yzzz_xy_1,   \
                             ta1_y_yzzz_yy_0,   \
                             ta1_y_yzzz_yy_1,   \
                             ta1_y_yzzzz_xy_0,  \
                             ta1_y_yzzzz_xy_1,  \
                             ta1_y_yzzzz_yy_0,  \
                             ta1_y_yzzzz_yy_1,  \
                             ta1_y_yzzzzz_xx_0, \
                             ta1_y_yzzzzz_xy_0, \
                             ta1_y_yzzzzz_xz_0, \
                             ta1_y_yzzzzz_yy_0, \
                             ta1_y_yzzzzz_yz_0, \
                             ta1_y_yzzzzz_zz_0, \
                             ta1_y_zzzzz_xx_0,  \
                             ta1_y_zzzzz_xx_1,  \
                             ta1_y_zzzzz_xz_0,  \
                             ta1_y_zzzzz_xz_1,  \
                             ta1_y_zzzzz_yz_0,  \
                             ta1_y_zzzzz_yz_1,  \
                             ta1_y_zzzzz_z_0,   \
                             ta1_y_zzzzz_z_1,   \
                             ta1_y_zzzzz_zz_0,  \
                             ta1_y_zzzzz_zz_1,  \
                             ta_zzzzz_xx_1,     \
                             ta_zzzzz_xz_1,     \
                             ta_zzzzz_yz_1,     \
                             ta_zzzzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzzzzz_xx_0[i] = ta_zzzzz_xx_1[i] + ta1_y_zzzzz_xx_0[i] * pa_y[i] - ta1_y_zzzzz_xx_1[i] * pc_y[i];

        ta1_y_yzzzzz_xy_0[i] =
            4.0 * ta1_y_yzzz_xy_0[i] * fe_0 - 4.0 * ta1_y_yzzz_xy_1[i] * fe_0 + ta1_y_yzzzz_xy_0[i] * pa_z[i] - ta1_y_yzzzz_xy_1[i] * pc_z[i];

        ta1_y_yzzzzz_xz_0[i] = ta_zzzzz_xz_1[i] + ta1_y_zzzzz_xz_0[i] * pa_y[i] - ta1_y_zzzzz_xz_1[i] * pc_y[i];

        ta1_y_yzzzzz_yy_0[i] =
            4.0 * ta1_y_yzzz_yy_0[i] * fe_0 - 4.0 * ta1_y_yzzz_yy_1[i] * fe_0 + ta1_y_yzzzz_yy_0[i] * pa_z[i] - ta1_y_yzzzz_yy_1[i] * pc_z[i];

        ta1_y_yzzzzz_yz_0[i] =
            ta1_y_zzzzz_z_0[i] * fe_0 - ta1_y_zzzzz_z_1[i] * fe_0 + ta_zzzzz_yz_1[i] + ta1_y_zzzzz_yz_0[i] * pa_y[i] - ta1_y_zzzzz_yz_1[i] * pc_y[i];

        ta1_y_yzzzzz_zz_0[i] = ta_zzzzz_zz_1[i] + ta1_y_zzzzz_zz_0[i] * pa_y[i] - ta1_y_zzzzz_zz_1[i] * pc_y[i];
    }

    // Set up 330-336 components of targeted buffer : ID

    auto ta1_y_zzzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 330);

    auto ta1_y_zzzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 331);

    auto ta1_y_zzzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 332);

    auto ta1_y_zzzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 333);

    auto ta1_y_zzzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 334);

    auto ta1_y_zzzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 335);

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta1_y_zzzz_xx_0,   \
                             ta1_y_zzzz_xx_1,   \
                             ta1_y_zzzz_xy_0,   \
                             ta1_y_zzzz_xy_1,   \
                             ta1_y_zzzz_xz_0,   \
                             ta1_y_zzzz_xz_1,   \
                             ta1_y_zzzz_yy_0,   \
                             ta1_y_zzzz_yy_1,   \
                             ta1_y_zzzz_yz_0,   \
                             ta1_y_zzzz_yz_1,   \
                             ta1_y_zzzz_zz_0,   \
                             ta1_y_zzzz_zz_1,   \
                             ta1_y_zzzzz_x_0,   \
                             ta1_y_zzzzz_x_1,   \
                             ta1_y_zzzzz_xx_0,  \
                             ta1_y_zzzzz_xx_1,  \
                             ta1_y_zzzzz_xy_0,  \
                             ta1_y_zzzzz_xy_1,  \
                             ta1_y_zzzzz_xz_0,  \
                             ta1_y_zzzzz_xz_1,  \
                             ta1_y_zzzzz_y_0,   \
                             ta1_y_zzzzz_y_1,   \
                             ta1_y_zzzzz_yy_0,  \
                             ta1_y_zzzzz_yy_1,  \
                             ta1_y_zzzzz_yz_0,  \
                             ta1_y_zzzzz_yz_1,  \
                             ta1_y_zzzzz_z_0,   \
                             ta1_y_zzzzz_z_1,   \
                             ta1_y_zzzzz_zz_0,  \
                             ta1_y_zzzzz_zz_1,  \
                             ta1_y_zzzzzz_xx_0, \
                             ta1_y_zzzzzz_xy_0, \
                             ta1_y_zzzzzz_xz_0, \
                             ta1_y_zzzzzz_yy_0, \
                             ta1_y_zzzzzz_yz_0, \
                             ta1_y_zzzzzz_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzzzzz_xx_0[i] =
            5.0 * ta1_y_zzzz_xx_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xx_1[i] * fe_0 + ta1_y_zzzzz_xx_0[i] * pa_z[i] - ta1_y_zzzzz_xx_1[i] * pc_z[i];

        ta1_y_zzzzzz_xy_0[i] =
            5.0 * ta1_y_zzzz_xy_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xy_1[i] * fe_0 + ta1_y_zzzzz_xy_0[i] * pa_z[i] - ta1_y_zzzzz_xy_1[i] * pc_z[i];

        ta1_y_zzzzzz_xz_0[i] = 5.0 * ta1_y_zzzz_xz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_xz_1[i] * fe_0 + ta1_y_zzzzz_x_0[i] * fe_0 -
                               ta1_y_zzzzz_x_1[i] * fe_0 + ta1_y_zzzzz_xz_0[i] * pa_z[i] - ta1_y_zzzzz_xz_1[i] * pc_z[i];

        ta1_y_zzzzzz_yy_0[i] =
            5.0 * ta1_y_zzzz_yy_0[i] * fe_0 - 5.0 * ta1_y_zzzz_yy_1[i] * fe_0 + ta1_y_zzzzz_yy_0[i] * pa_z[i] - ta1_y_zzzzz_yy_1[i] * pc_z[i];

        ta1_y_zzzzzz_yz_0[i] = 5.0 * ta1_y_zzzz_yz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_yz_1[i] * fe_0 + ta1_y_zzzzz_y_0[i] * fe_0 -
                               ta1_y_zzzzz_y_1[i] * fe_0 + ta1_y_zzzzz_yz_0[i] * pa_z[i] - ta1_y_zzzzz_yz_1[i] * pc_z[i];

        ta1_y_zzzzzz_zz_0[i] = 5.0 * ta1_y_zzzz_zz_0[i] * fe_0 - 5.0 * ta1_y_zzzz_zz_1[i] * fe_0 + 2.0 * ta1_y_zzzzz_z_0[i] * fe_0 -
                               2.0 * ta1_y_zzzzz_z_1[i] * fe_0 + ta1_y_zzzzz_zz_0[i] * pa_z[i] - ta1_y_zzzzz_zz_1[i] * pc_z[i];
    }

    // Set up 336-342 components of targeted buffer : ID

    auto ta1_z_xxxxxx_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 336);

    auto ta1_z_xxxxxx_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 337);

    auto ta1_z_xxxxxx_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 338);

    auto ta1_z_xxxxxx_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 339);

    auto ta1_z_xxxxxx_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 340);

    auto ta1_z_xxxxxx_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 341);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_z_xxxx_xx_0,   \
                             ta1_z_xxxx_xx_1,   \
                             ta1_z_xxxx_xy_0,   \
                             ta1_z_xxxx_xy_1,   \
                             ta1_z_xxxx_xz_0,   \
                             ta1_z_xxxx_xz_1,   \
                             ta1_z_xxxx_yy_0,   \
                             ta1_z_xxxx_yy_1,   \
                             ta1_z_xxxx_yz_0,   \
                             ta1_z_xxxx_yz_1,   \
                             ta1_z_xxxx_zz_0,   \
                             ta1_z_xxxx_zz_1,   \
                             ta1_z_xxxxx_x_0,   \
                             ta1_z_xxxxx_x_1,   \
                             ta1_z_xxxxx_xx_0,  \
                             ta1_z_xxxxx_xx_1,  \
                             ta1_z_xxxxx_xy_0,  \
                             ta1_z_xxxxx_xy_1,  \
                             ta1_z_xxxxx_xz_0,  \
                             ta1_z_xxxxx_xz_1,  \
                             ta1_z_xxxxx_y_0,   \
                             ta1_z_xxxxx_y_1,   \
                             ta1_z_xxxxx_yy_0,  \
                             ta1_z_xxxxx_yy_1,  \
                             ta1_z_xxxxx_yz_0,  \
                             ta1_z_xxxxx_yz_1,  \
                             ta1_z_xxxxx_z_0,   \
                             ta1_z_xxxxx_z_1,   \
                             ta1_z_xxxxx_zz_0,  \
                             ta1_z_xxxxx_zz_1,  \
                             ta1_z_xxxxxx_xx_0, \
                             ta1_z_xxxxxx_xy_0, \
                             ta1_z_xxxxxx_xz_0, \
                             ta1_z_xxxxxx_yy_0, \
                             ta1_z_xxxxxx_yz_0, \
                             ta1_z_xxxxxx_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxxx_xx_0[i] = 5.0 * ta1_z_xxxx_xx_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xx_1[i] * fe_0 + 2.0 * ta1_z_xxxxx_x_0[i] * fe_0 -
                               2.0 * ta1_z_xxxxx_x_1[i] * fe_0 + ta1_z_xxxxx_xx_0[i] * pa_x[i] - ta1_z_xxxxx_xx_1[i] * pc_x[i];

        ta1_z_xxxxxx_xy_0[i] = 5.0 * ta1_z_xxxx_xy_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xy_1[i] * fe_0 + ta1_z_xxxxx_y_0[i] * fe_0 -
                               ta1_z_xxxxx_y_1[i] * fe_0 + ta1_z_xxxxx_xy_0[i] * pa_x[i] - ta1_z_xxxxx_xy_1[i] * pc_x[i];

        ta1_z_xxxxxx_xz_0[i] = 5.0 * ta1_z_xxxx_xz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_xz_1[i] * fe_0 + ta1_z_xxxxx_z_0[i] * fe_0 -
                               ta1_z_xxxxx_z_1[i] * fe_0 + ta1_z_xxxxx_xz_0[i] * pa_x[i] - ta1_z_xxxxx_xz_1[i] * pc_x[i];

        ta1_z_xxxxxx_yy_0[i] =
            5.0 * ta1_z_xxxx_yy_0[i] * fe_0 - 5.0 * ta1_z_xxxx_yy_1[i] * fe_0 + ta1_z_xxxxx_yy_0[i] * pa_x[i] - ta1_z_xxxxx_yy_1[i] * pc_x[i];

        ta1_z_xxxxxx_yz_0[i] =
            5.0 * ta1_z_xxxx_yz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_yz_1[i] * fe_0 + ta1_z_xxxxx_yz_0[i] * pa_x[i] - ta1_z_xxxxx_yz_1[i] * pc_x[i];

        ta1_z_xxxxxx_zz_0[i] =
            5.0 * ta1_z_xxxx_zz_0[i] * fe_0 - 5.0 * ta1_z_xxxx_zz_1[i] * fe_0 + ta1_z_xxxxx_zz_0[i] * pa_x[i] - ta1_z_xxxxx_zz_1[i] * pc_x[i];
    }

    // Set up 342-348 components of targeted buffer : ID

    auto ta1_z_xxxxxy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 342);

    auto ta1_z_xxxxxy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 343);

    auto ta1_z_xxxxxy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 344);

    auto ta1_z_xxxxxy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 345);

    auto ta1_z_xxxxxy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 346);

    auto ta1_z_xxxxxy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 347);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_z_xxxxx_x_0,   \
                             ta1_z_xxxxx_x_1,   \
                             ta1_z_xxxxx_xx_0,  \
                             ta1_z_xxxxx_xx_1,  \
                             ta1_z_xxxxx_xy_0,  \
                             ta1_z_xxxxx_xy_1,  \
                             ta1_z_xxxxx_xz_0,  \
                             ta1_z_xxxxx_xz_1,  \
                             ta1_z_xxxxx_zz_0,  \
                             ta1_z_xxxxx_zz_1,  \
                             ta1_z_xxxxxy_xx_0, \
                             ta1_z_xxxxxy_xy_0, \
                             ta1_z_xxxxxy_xz_0, \
                             ta1_z_xxxxxy_yy_0, \
                             ta1_z_xxxxxy_yz_0, \
                             ta1_z_xxxxxy_zz_0, \
                             ta1_z_xxxxy_yy_0,  \
                             ta1_z_xxxxy_yy_1,  \
                             ta1_z_xxxxy_yz_0,  \
                             ta1_z_xxxxy_yz_1,  \
                             ta1_z_xxxy_yy_0,   \
                             ta1_z_xxxy_yy_1,   \
                             ta1_z_xxxy_yz_0,   \
                             ta1_z_xxxy_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxxy_xx_0[i] = ta1_z_xxxxx_xx_0[i] * pa_y[i] - ta1_z_xxxxx_xx_1[i] * pc_y[i];

        ta1_z_xxxxxy_xy_0[i] = ta1_z_xxxxx_x_0[i] * fe_0 - ta1_z_xxxxx_x_1[i] * fe_0 + ta1_z_xxxxx_xy_0[i] * pa_y[i] - ta1_z_xxxxx_xy_1[i] * pc_y[i];

        ta1_z_xxxxxy_xz_0[i] = ta1_z_xxxxx_xz_0[i] * pa_y[i] - ta1_z_xxxxx_xz_1[i] * pc_y[i];

        ta1_z_xxxxxy_yy_0[i] =
            4.0 * ta1_z_xxxy_yy_0[i] * fe_0 - 4.0 * ta1_z_xxxy_yy_1[i] * fe_0 + ta1_z_xxxxy_yy_0[i] * pa_x[i] - ta1_z_xxxxy_yy_1[i] * pc_x[i];

        ta1_z_xxxxxy_yz_0[i] =
            4.0 * ta1_z_xxxy_yz_0[i] * fe_0 - 4.0 * ta1_z_xxxy_yz_1[i] * fe_0 + ta1_z_xxxxy_yz_0[i] * pa_x[i] - ta1_z_xxxxy_yz_1[i] * pc_x[i];

        ta1_z_xxxxxy_zz_0[i] = ta1_z_xxxxx_zz_0[i] * pa_y[i] - ta1_z_xxxxx_zz_1[i] * pc_y[i];
    }

    // Set up 348-354 components of targeted buffer : ID

    auto ta1_z_xxxxxz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 348);

    auto ta1_z_xxxxxz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 349);

    auto ta1_z_xxxxxz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 350);

    auto ta1_z_xxxxxz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 351);

    auto ta1_z_xxxxxz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 352);

    auto ta1_z_xxxxxz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 353);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_z_xxxxx_x_0,   \
                             ta1_z_xxxxx_x_1,   \
                             ta1_z_xxxxx_xx_0,  \
                             ta1_z_xxxxx_xx_1,  \
                             ta1_z_xxxxx_xy_0,  \
                             ta1_z_xxxxx_xy_1,  \
                             ta1_z_xxxxx_xz_0,  \
                             ta1_z_xxxxx_xz_1,  \
                             ta1_z_xxxxx_yy_0,  \
                             ta1_z_xxxxx_yy_1,  \
                             ta1_z_xxxxxz_xx_0, \
                             ta1_z_xxxxxz_xy_0, \
                             ta1_z_xxxxxz_xz_0, \
                             ta1_z_xxxxxz_yy_0, \
                             ta1_z_xxxxxz_yz_0, \
                             ta1_z_xxxxxz_zz_0, \
                             ta1_z_xxxxz_yz_0,  \
                             ta1_z_xxxxz_yz_1,  \
                             ta1_z_xxxxz_zz_0,  \
                             ta1_z_xxxxz_zz_1,  \
                             ta1_z_xxxz_yz_0,   \
                             ta1_z_xxxz_yz_1,   \
                             ta1_z_xxxz_zz_0,   \
                             ta1_z_xxxz_zz_1,   \
                             ta_xxxxx_xx_1,     \
                             ta_xxxxx_xy_1,     \
                             ta_xxxxx_xz_1,     \
                             ta_xxxxx_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxxz_xx_0[i] = ta_xxxxx_xx_1[i] + ta1_z_xxxxx_xx_0[i] * pa_z[i] - ta1_z_xxxxx_xx_1[i] * pc_z[i];

        ta1_z_xxxxxz_xy_0[i] = ta_xxxxx_xy_1[i] + ta1_z_xxxxx_xy_0[i] * pa_z[i] - ta1_z_xxxxx_xy_1[i] * pc_z[i];

        ta1_z_xxxxxz_xz_0[i] =
            ta1_z_xxxxx_x_0[i] * fe_0 - ta1_z_xxxxx_x_1[i] * fe_0 + ta_xxxxx_xz_1[i] + ta1_z_xxxxx_xz_0[i] * pa_z[i] - ta1_z_xxxxx_xz_1[i] * pc_z[i];

        ta1_z_xxxxxz_yy_0[i] = ta_xxxxx_yy_1[i] + ta1_z_xxxxx_yy_0[i] * pa_z[i] - ta1_z_xxxxx_yy_1[i] * pc_z[i];

        ta1_z_xxxxxz_yz_0[i] =
            4.0 * ta1_z_xxxz_yz_0[i] * fe_0 - 4.0 * ta1_z_xxxz_yz_1[i] * fe_0 + ta1_z_xxxxz_yz_0[i] * pa_x[i] - ta1_z_xxxxz_yz_1[i] * pc_x[i];

        ta1_z_xxxxxz_zz_0[i] =
            4.0 * ta1_z_xxxz_zz_0[i] * fe_0 - 4.0 * ta1_z_xxxz_zz_1[i] * fe_0 + ta1_z_xxxxz_zz_0[i] * pa_x[i] - ta1_z_xxxxz_zz_1[i] * pc_x[i];
    }

    // Set up 354-360 components of targeted buffer : ID

    auto ta1_z_xxxxyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 354);

    auto ta1_z_xxxxyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 355);

    auto ta1_z_xxxxyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 356);

    auto ta1_z_xxxxyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 357);

    auto ta1_z_xxxxyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 358);

    auto ta1_z_xxxxyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 359);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_z_xxxx_xx_0,   \
                             ta1_z_xxxx_xx_1,   \
                             ta1_z_xxxx_xz_0,   \
                             ta1_z_xxxx_xz_1,   \
                             ta1_z_xxxxy_xx_0,  \
                             ta1_z_xxxxy_xx_1,  \
                             ta1_z_xxxxy_xz_0,  \
                             ta1_z_xxxxy_xz_1,  \
                             ta1_z_xxxxyy_xx_0, \
                             ta1_z_xxxxyy_xy_0, \
                             ta1_z_xxxxyy_xz_0, \
                             ta1_z_xxxxyy_yy_0, \
                             ta1_z_xxxxyy_yz_0, \
                             ta1_z_xxxxyy_zz_0, \
                             ta1_z_xxxyy_xy_0,  \
                             ta1_z_xxxyy_xy_1,  \
                             ta1_z_xxxyy_y_0,   \
                             ta1_z_xxxyy_y_1,   \
                             ta1_z_xxxyy_yy_0,  \
                             ta1_z_xxxyy_yy_1,  \
                             ta1_z_xxxyy_yz_0,  \
                             ta1_z_xxxyy_yz_1,  \
                             ta1_z_xxxyy_zz_0,  \
                             ta1_z_xxxyy_zz_1,  \
                             ta1_z_xxyy_xy_0,   \
                             ta1_z_xxyy_xy_1,   \
                             ta1_z_xxyy_yy_0,   \
                             ta1_z_xxyy_yy_1,   \
                             ta1_z_xxyy_yz_0,   \
                             ta1_z_xxyy_yz_1,   \
                             ta1_z_xxyy_zz_0,   \
                             ta1_z_xxyy_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxyy_xx_0[i] = ta1_z_xxxx_xx_0[i] * fe_0 - ta1_z_xxxx_xx_1[i] * fe_0 + ta1_z_xxxxy_xx_0[i] * pa_y[i] - ta1_z_xxxxy_xx_1[i] * pc_y[i];

        ta1_z_xxxxyy_xy_0[i] = 3.0 * ta1_z_xxyy_xy_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xy_1[i] * fe_0 + ta1_z_xxxyy_y_0[i] * fe_0 -
                               ta1_z_xxxyy_y_1[i] * fe_0 + ta1_z_xxxyy_xy_0[i] * pa_x[i] - ta1_z_xxxyy_xy_1[i] * pc_x[i];

        ta1_z_xxxxyy_xz_0[i] = ta1_z_xxxx_xz_0[i] * fe_0 - ta1_z_xxxx_xz_1[i] * fe_0 + ta1_z_xxxxy_xz_0[i] * pa_y[i] - ta1_z_xxxxy_xz_1[i] * pc_y[i];

        ta1_z_xxxxyy_yy_0[i] =
            3.0 * ta1_z_xxyy_yy_0[i] * fe_0 - 3.0 * ta1_z_xxyy_yy_1[i] * fe_0 + ta1_z_xxxyy_yy_0[i] * pa_x[i] - ta1_z_xxxyy_yy_1[i] * pc_x[i];

        ta1_z_xxxxyy_yz_0[i] =
            3.0 * ta1_z_xxyy_yz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_yz_1[i] * fe_0 + ta1_z_xxxyy_yz_0[i] * pa_x[i] - ta1_z_xxxyy_yz_1[i] * pc_x[i];

        ta1_z_xxxxyy_zz_0[i] =
            3.0 * ta1_z_xxyy_zz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_zz_1[i] * fe_0 + ta1_z_xxxyy_zz_0[i] * pa_x[i] - ta1_z_xxxyy_zz_1[i] * pc_x[i];
    }

    // Set up 360-366 components of targeted buffer : ID

    auto ta1_z_xxxxyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 360);

    auto ta1_z_xxxxyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 361);

    auto ta1_z_xxxxyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 362);

    auto ta1_z_xxxxyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 363);

    auto ta1_z_xxxxyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 364);

    auto ta1_z_xxxxyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 365);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_z_xxxxy_xy_0,  \
                             ta1_z_xxxxy_xy_1,  \
                             ta1_z_xxxxy_yy_0,  \
                             ta1_z_xxxxy_yy_1,  \
                             ta1_z_xxxxyz_xx_0, \
                             ta1_z_xxxxyz_xy_0, \
                             ta1_z_xxxxyz_xz_0, \
                             ta1_z_xxxxyz_yy_0, \
                             ta1_z_xxxxyz_yz_0, \
                             ta1_z_xxxxyz_zz_0, \
                             ta1_z_xxxxz_xx_0,  \
                             ta1_z_xxxxz_xx_1,  \
                             ta1_z_xxxxz_xz_0,  \
                             ta1_z_xxxxz_xz_1,  \
                             ta1_z_xxxxz_zz_0,  \
                             ta1_z_xxxxz_zz_1,  \
                             ta1_z_xxxyz_yz_0,  \
                             ta1_z_xxxyz_yz_1,  \
                             ta1_z_xxyz_yz_0,   \
                             ta1_z_xxyz_yz_1,   \
                             ta_xxxxy_xy_1,     \
                             ta_xxxxy_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxyz_xx_0[i] = ta1_z_xxxxz_xx_0[i] * pa_y[i] - ta1_z_xxxxz_xx_1[i] * pc_y[i];

        ta1_z_xxxxyz_xy_0[i] = ta_xxxxy_xy_1[i] + ta1_z_xxxxy_xy_0[i] * pa_z[i] - ta1_z_xxxxy_xy_1[i] * pc_z[i];

        ta1_z_xxxxyz_xz_0[i] = ta1_z_xxxxz_xz_0[i] * pa_y[i] - ta1_z_xxxxz_xz_1[i] * pc_y[i];

        ta1_z_xxxxyz_yy_0[i] = ta_xxxxy_yy_1[i] + ta1_z_xxxxy_yy_0[i] * pa_z[i] - ta1_z_xxxxy_yy_1[i] * pc_z[i];

        ta1_z_xxxxyz_yz_0[i] =
            3.0 * ta1_z_xxyz_yz_0[i] * fe_0 - 3.0 * ta1_z_xxyz_yz_1[i] * fe_0 + ta1_z_xxxyz_yz_0[i] * pa_x[i] - ta1_z_xxxyz_yz_1[i] * pc_x[i];

        ta1_z_xxxxyz_zz_0[i] = ta1_z_xxxxz_zz_0[i] * pa_y[i] - ta1_z_xxxxz_zz_1[i] * pc_y[i];
    }

    // Set up 366-372 components of targeted buffer : ID

    auto ta1_z_xxxxzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 366);

    auto ta1_z_xxxxzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 367);

    auto ta1_z_xxxxzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 368);

    auto ta1_z_xxxxzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 369);

    auto ta1_z_xxxxzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 370);

    auto ta1_z_xxxxzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 371);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_z_xxxx_xx_0,   \
                             ta1_z_xxxx_xx_1,   \
                             ta1_z_xxxx_xy_0,   \
                             ta1_z_xxxx_xy_1,   \
                             ta1_z_xxxxz_xx_0,  \
                             ta1_z_xxxxz_xx_1,  \
                             ta1_z_xxxxz_xy_0,  \
                             ta1_z_xxxxz_xy_1,  \
                             ta1_z_xxxxzz_xx_0, \
                             ta1_z_xxxxzz_xy_0, \
                             ta1_z_xxxxzz_xz_0, \
                             ta1_z_xxxxzz_yy_0, \
                             ta1_z_xxxxzz_yz_0, \
                             ta1_z_xxxxzz_zz_0, \
                             ta1_z_xxxzz_xz_0,  \
                             ta1_z_xxxzz_xz_1,  \
                             ta1_z_xxxzz_yy_0,  \
                             ta1_z_xxxzz_yy_1,  \
                             ta1_z_xxxzz_yz_0,  \
                             ta1_z_xxxzz_yz_1,  \
                             ta1_z_xxxzz_z_0,   \
                             ta1_z_xxxzz_z_1,   \
                             ta1_z_xxxzz_zz_0,  \
                             ta1_z_xxxzz_zz_1,  \
                             ta1_z_xxzz_xz_0,   \
                             ta1_z_xxzz_xz_1,   \
                             ta1_z_xxzz_yy_0,   \
                             ta1_z_xxzz_yy_1,   \
                             ta1_z_xxzz_yz_0,   \
                             ta1_z_xxzz_yz_1,   \
                             ta1_z_xxzz_zz_0,   \
                             ta1_z_xxzz_zz_1,   \
                             ta_xxxxz_xx_1,     \
                             ta_xxxxz_xy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxxzz_xx_0[i] =
            ta1_z_xxxx_xx_0[i] * fe_0 - ta1_z_xxxx_xx_1[i] * fe_0 + ta_xxxxz_xx_1[i] + ta1_z_xxxxz_xx_0[i] * pa_z[i] - ta1_z_xxxxz_xx_1[i] * pc_z[i];

        ta1_z_xxxxzz_xy_0[i] =
            ta1_z_xxxx_xy_0[i] * fe_0 - ta1_z_xxxx_xy_1[i] * fe_0 + ta_xxxxz_xy_1[i] + ta1_z_xxxxz_xy_0[i] * pa_z[i] - ta1_z_xxxxz_xy_1[i] * pc_z[i];

        ta1_z_xxxxzz_xz_0[i] = 3.0 * ta1_z_xxzz_xz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xz_1[i] * fe_0 + ta1_z_xxxzz_z_0[i] * fe_0 -
                               ta1_z_xxxzz_z_1[i] * fe_0 + ta1_z_xxxzz_xz_0[i] * pa_x[i] - ta1_z_xxxzz_xz_1[i] * pc_x[i];

        ta1_z_xxxxzz_yy_0[i] =
            3.0 * ta1_z_xxzz_yy_0[i] * fe_0 - 3.0 * ta1_z_xxzz_yy_1[i] * fe_0 + ta1_z_xxxzz_yy_0[i] * pa_x[i] - ta1_z_xxxzz_yy_1[i] * pc_x[i];

        ta1_z_xxxxzz_yz_0[i] =
            3.0 * ta1_z_xxzz_yz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_yz_1[i] * fe_0 + ta1_z_xxxzz_yz_0[i] * pa_x[i] - ta1_z_xxxzz_yz_1[i] * pc_x[i];

        ta1_z_xxxxzz_zz_0[i] =
            3.0 * ta1_z_xxzz_zz_0[i] * fe_0 - 3.0 * ta1_z_xxzz_zz_1[i] * fe_0 + ta1_z_xxxzz_zz_0[i] * pa_x[i] - ta1_z_xxxzz_zz_1[i] * pc_x[i];
    }

    // Set up 372-378 components of targeted buffer : ID

    auto ta1_z_xxxyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 372);

    auto ta1_z_xxxyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 373);

    auto ta1_z_xxxyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 374);

    auto ta1_z_xxxyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 375);

    auto ta1_z_xxxyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 376);

    auto ta1_z_xxxyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 377);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_z_xxxy_xx_0,   \
                             ta1_z_xxxy_xx_1,   \
                             ta1_z_xxxy_xz_0,   \
                             ta1_z_xxxy_xz_1,   \
                             ta1_z_xxxyy_xx_0,  \
                             ta1_z_xxxyy_xx_1,  \
                             ta1_z_xxxyy_xz_0,  \
                             ta1_z_xxxyy_xz_1,  \
                             ta1_z_xxxyyy_xx_0, \
                             ta1_z_xxxyyy_xy_0, \
                             ta1_z_xxxyyy_xz_0, \
                             ta1_z_xxxyyy_yy_0, \
                             ta1_z_xxxyyy_yz_0, \
                             ta1_z_xxxyyy_zz_0, \
                             ta1_z_xxyyy_xy_0,  \
                             ta1_z_xxyyy_xy_1,  \
                             ta1_z_xxyyy_y_0,   \
                             ta1_z_xxyyy_y_1,   \
                             ta1_z_xxyyy_yy_0,  \
                             ta1_z_xxyyy_yy_1,  \
                             ta1_z_xxyyy_yz_0,  \
                             ta1_z_xxyyy_yz_1,  \
                             ta1_z_xxyyy_zz_0,  \
                             ta1_z_xxyyy_zz_1,  \
                             ta1_z_xyyy_xy_0,   \
                             ta1_z_xyyy_xy_1,   \
                             ta1_z_xyyy_yy_0,   \
                             ta1_z_xyyy_yy_1,   \
                             ta1_z_xyyy_yz_0,   \
                             ta1_z_xyyy_yz_1,   \
                             ta1_z_xyyy_zz_0,   \
                             ta1_z_xyyy_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyyy_xx_0[i] =
            2.0 * ta1_z_xxxy_xx_0[i] * fe_0 - 2.0 * ta1_z_xxxy_xx_1[i] * fe_0 + ta1_z_xxxyy_xx_0[i] * pa_y[i] - ta1_z_xxxyy_xx_1[i] * pc_y[i];

        ta1_z_xxxyyy_xy_0[i] = 2.0 * ta1_z_xyyy_xy_0[i] * fe_0 - 2.0 * ta1_z_xyyy_xy_1[i] * fe_0 + ta1_z_xxyyy_y_0[i] * fe_0 -
                               ta1_z_xxyyy_y_1[i] * fe_0 + ta1_z_xxyyy_xy_0[i] * pa_x[i] - ta1_z_xxyyy_xy_1[i] * pc_x[i];

        ta1_z_xxxyyy_xz_0[i] =
            2.0 * ta1_z_xxxy_xz_0[i] * fe_0 - 2.0 * ta1_z_xxxy_xz_1[i] * fe_0 + ta1_z_xxxyy_xz_0[i] * pa_y[i] - ta1_z_xxxyy_xz_1[i] * pc_y[i];

        ta1_z_xxxyyy_yy_0[i] =
            2.0 * ta1_z_xyyy_yy_0[i] * fe_0 - 2.0 * ta1_z_xyyy_yy_1[i] * fe_0 + ta1_z_xxyyy_yy_0[i] * pa_x[i] - ta1_z_xxyyy_yy_1[i] * pc_x[i];

        ta1_z_xxxyyy_yz_0[i] =
            2.0 * ta1_z_xyyy_yz_0[i] * fe_0 - 2.0 * ta1_z_xyyy_yz_1[i] * fe_0 + ta1_z_xxyyy_yz_0[i] * pa_x[i] - ta1_z_xxyyy_yz_1[i] * pc_x[i];

        ta1_z_xxxyyy_zz_0[i] =
            2.0 * ta1_z_xyyy_zz_0[i] * fe_0 - 2.0 * ta1_z_xyyy_zz_1[i] * fe_0 + ta1_z_xxyyy_zz_0[i] * pa_x[i] - ta1_z_xxyyy_zz_1[i] * pc_x[i];
    }

    // Set up 378-384 components of targeted buffer : ID

    auto ta1_z_xxxyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 378);

    auto ta1_z_xxxyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 379);

    auto ta1_z_xxxyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 380);

    auto ta1_z_xxxyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 381);

    auto ta1_z_xxxyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 382);

    auto ta1_z_xxxyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 383);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_z_xxxyy_xx_0,  \
                             ta1_z_xxxyy_xx_1,  \
                             ta1_z_xxxyy_xy_0,  \
                             ta1_z_xxxyy_xy_1,  \
                             ta1_z_xxxyy_yy_0,  \
                             ta1_z_xxxyy_yy_1,  \
                             ta1_z_xxxyyz_xx_0, \
                             ta1_z_xxxyyz_xy_0, \
                             ta1_z_xxxyyz_xz_0, \
                             ta1_z_xxxyyz_yy_0, \
                             ta1_z_xxxyyz_yz_0, \
                             ta1_z_xxxyyz_zz_0, \
                             ta1_z_xxxyz_xz_0,  \
                             ta1_z_xxxyz_xz_1,  \
                             ta1_z_xxxz_xz_0,   \
                             ta1_z_xxxz_xz_1,   \
                             ta1_z_xxyyz_yz_0,  \
                             ta1_z_xxyyz_yz_1,  \
                             ta1_z_xxyyz_zz_0,  \
                             ta1_z_xxyyz_zz_1,  \
                             ta1_z_xyyz_yz_0,   \
                             ta1_z_xyyz_yz_1,   \
                             ta1_z_xyyz_zz_0,   \
                             ta1_z_xyyz_zz_1,   \
                             ta_xxxyy_xx_1,     \
                             ta_xxxyy_xy_1,     \
                             ta_xxxyy_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyyz_xx_0[i] = ta_xxxyy_xx_1[i] + ta1_z_xxxyy_xx_0[i] * pa_z[i] - ta1_z_xxxyy_xx_1[i] * pc_z[i];

        ta1_z_xxxyyz_xy_0[i] = ta_xxxyy_xy_1[i] + ta1_z_xxxyy_xy_0[i] * pa_z[i] - ta1_z_xxxyy_xy_1[i] * pc_z[i];

        ta1_z_xxxyyz_xz_0[i] = ta1_z_xxxz_xz_0[i] * fe_0 - ta1_z_xxxz_xz_1[i] * fe_0 + ta1_z_xxxyz_xz_0[i] * pa_y[i] - ta1_z_xxxyz_xz_1[i] * pc_y[i];

        ta1_z_xxxyyz_yy_0[i] = ta_xxxyy_yy_1[i] + ta1_z_xxxyy_yy_0[i] * pa_z[i] - ta1_z_xxxyy_yy_1[i] * pc_z[i];

        ta1_z_xxxyyz_yz_0[i] =
            2.0 * ta1_z_xyyz_yz_0[i] * fe_0 - 2.0 * ta1_z_xyyz_yz_1[i] * fe_0 + ta1_z_xxyyz_yz_0[i] * pa_x[i] - ta1_z_xxyyz_yz_1[i] * pc_x[i];

        ta1_z_xxxyyz_zz_0[i] =
            2.0 * ta1_z_xyyz_zz_0[i] * fe_0 - 2.0 * ta1_z_xyyz_zz_1[i] * fe_0 + ta1_z_xxyyz_zz_0[i] * pa_x[i] - ta1_z_xxyyz_zz_1[i] * pc_x[i];
    }

    // Set up 384-390 components of targeted buffer : ID

    auto ta1_z_xxxyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 384);

    auto ta1_z_xxxyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 385);

    auto ta1_z_xxxyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 386);

    auto ta1_z_xxxyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 387);

    auto ta1_z_xxxyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 388);

    auto ta1_z_xxxyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 389);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_z_xxxyzz_xx_0, \
                             ta1_z_xxxyzz_xy_0, \
                             ta1_z_xxxyzz_xz_0, \
                             ta1_z_xxxyzz_yy_0, \
                             ta1_z_xxxyzz_yz_0, \
                             ta1_z_xxxyzz_zz_0, \
                             ta1_z_xxxzz_x_0,   \
                             ta1_z_xxxzz_x_1,   \
                             ta1_z_xxxzz_xx_0,  \
                             ta1_z_xxxzz_xx_1,  \
                             ta1_z_xxxzz_xy_0,  \
                             ta1_z_xxxzz_xy_1,  \
                             ta1_z_xxxzz_xz_0,  \
                             ta1_z_xxxzz_xz_1,  \
                             ta1_z_xxxzz_zz_0,  \
                             ta1_z_xxxzz_zz_1,  \
                             ta1_z_xxyzz_yy_0,  \
                             ta1_z_xxyzz_yy_1,  \
                             ta1_z_xxyzz_yz_0,  \
                             ta1_z_xxyzz_yz_1,  \
                             ta1_z_xyzz_yy_0,   \
                             ta1_z_xyzz_yy_1,   \
                             ta1_z_xyzz_yz_0,   \
                             ta1_z_xyzz_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxyzz_xx_0[i] = ta1_z_xxxzz_xx_0[i] * pa_y[i] - ta1_z_xxxzz_xx_1[i] * pc_y[i];

        ta1_z_xxxyzz_xy_0[i] = ta1_z_xxxzz_x_0[i] * fe_0 - ta1_z_xxxzz_x_1[i] * fe_0 + ta1_z_xxxzz_xy_0[i] * pa_y[i] - ta1_z_xxxzz_xy_1[i] * pc_y[i];

        ta1_z_xxxyzz_xz_0[i] = ta1_z_xxxzz_xz_0[i] * pa_y[i] - ta1_z_xxxzz_xz_1[i] * pc_y[i];

        ta1_z_xxxyzz_yy_0[i] =
            2.0 * ta1_z_xyzz_yy_0[i] * fe_0 - 2.0 * ta1_z_xyzz_yy_1[i] * fe_0 + ta1_z_xxyzz_yy_0[i] * pa_x[i] - ta1_z_xxyzz_yy_1[i] * pc_x[i];

        ta1_z_xxxyzz_yz_0[i] =
            2.0 * ta1_z_xyzz_yz_0[i] * fe_0 - 2.0 * ta1_z_xyzz_yz_1[i] * fe_0 + ta1_z_xxyzz_yz_0[i] * pa_x[i] - ta1_z_xxyzz_yz_1[i] * pc_x[i];

        ta1_z_xxxyzz_zz_0[i] = ta1_z_xxxzz_zz_0[i] * pa_y[i] - ta1_z_xxxzz_zz_1[i] * pc_y[i];
    }

    // Set up 390-396 components of targeted buffer : ID

    auto ta1_z_xxxzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 390);

    auto ta1_z_xxxzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 391);

    auto ta1_z_xxxzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 392);

    auto ta1_z_xxxzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 393);

    auto ta1_z_xxxzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 394);

    auto ta1_z_xxxzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 395);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_z_xxxz_xx_0,   \
                             ta1_z_xxxz_xx_1,   \
                             ta1_z_xxxz_xy_0,   \
                             ta1_z_xxxz_xy_1,   \
                             ta1_z_xxxzz_xx_0,  \
                             ta1_z_xxxzz_xx_1,  \
                             ta1_z_xxxzz_xy_0,  \
                             ta1_z_xxxzz_xy_1,  \
                             ta1_z_xxxzzz_xx_0, \
                             ta1_z_xxxzzz_xy_0, \
                             ta1_z_xxxzzz_xz_0, \
                             ta1_z_xxxzzz_yy_0, \
                             ta1_z_xxxzzz_yz_0, \
                             ta1_z_xxxzzz_zz_0, \
                             ta1_z_xxzzz_xz_0,  \
                             ta1_z_xxzzz_xz_1,  \
                             ta1_z_xxzzz_yy_0,  \
                             ta1_z_xxzzz_yy_1,  \
                             ta1_z_xxzzz_yz_0,  \
                             ta1_z_xxzzz_yz_1,  \
                             ta1_z_xxzzz_z_0,   \
                             ta1_z_xxzzz_z_1,   \
                             ta1_z_xxzzz_zz_0,  \
                             ta1_z_xxzzz_zz_1,  \
                             ta1_z_xzzz_xz_0,   \
                             ta1_z_xzzz_xz_1,   \
                             ta1_z_xzzz_yy_0,   \
                             ta1_z_xzzz_yy_1,   \
                             ta1_z_xzzz_yz_0,   \
                             ta1_z_xzzz_yz_1,   \
                             ta1_z_xzzz_zz_0,   \
                             ta1_z_xzzz_zz_1,   \
                             ta_xxxzz_xx_1,     \
                             ta_xxxzz_xy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxzzz_xx_0[i] = 2.0 * ta1_z_xxxz_xx_0[i] * fe_0 - 2.0 * ta1_z_xxxz_xx_1[i] * fe_0 + ta_xxxzz_xx_1[i] + ta1_z_xxxzz_xx_0[i] * pa_z[i] -
                               ta1_z_xxxzz_xx_1[i] * pc_z[i];

        ta1_z_xxxzzz_xy_0[i] = 2.0 * ta1_z_xxxz_xy_0[i] * fe_0 - 2.0 * ta1_z_xxxz_xy_1[i] * fe_0 + ta_xxxzz_xy_1[i] + ta1_z_xxxzz_xy_0[i] * pa_z[i] -
                               ta1_z_xxxzz_xy_1[i] * pc_z[i];

        ta1_z_xxxzzz_xz_0[i] = 2.0 * ta1_z_xzzz_xz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_xz_1[i] * fe_0 + ta1_z_xxzzz_z_0[i] * fe_0 -
                               ta1_z_xxzzz_z_1[i] * fe_0 + ta1_z_xxzzz_xz_0[i] * pa_x[i] - ta1_z_xxzzz_xz_1[i] * pc_x[i];

        ta1_z_xxxzzz_yy_0[i] =
            2.0 * ta1_z_xzzz_yy_0[i] * fe_0 - 2.0 * ta1_z_xzzz_yy_1[i] * fe_0 + ta1_z_xxzzz_yy_0[i] * pa_x[i] - ta1_z_xxzzz_yy_1[i] * pc_x[i];

        ta1_z_xxxzzz_yz_0[i] =
            2.0 * ta1_z_xzzz_yz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_yz_1[i] * fe_0 + ta1_z_xxzzz_yz_0[i] * pa_x[i] - ta1_z_xxzzz_yz_1[i] * pc_x[i];

        ta1_z_xxxzzz_zz_0[i] =
            2.0 * ta1_z_xzzz_zz_0[i] * fe_0 - 2.0 * ta1_z_xzzz_zz_1[i] * fe_0 + ta1_z_xxzzz_zz_0[i] * pa_x[i] - ta1_z_xxzzz_zz_1[i] * pc_x[i];
    }

    // Set up 396-402 components of targeted buffer : ID

    auto ta1_z_xxyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 396);

    auto ta1_z_xxyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 397);

    auto ta1_z_xxyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 398);

    auto ta1_z_xxyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 399);

    auto ta1_z_xxyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 400);

    auto ta1_z_xxyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 401);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_z_xxyy_xx_0,   \
                             ta1_z_xxyy_xx_1,   \
                             ta1_z_xxyy_xz_0,   \
                             ta1_z_xxyy_xz_1,   \
                             ta1_z_xxyyy_xx_0,  \
                             ta1_z_xxyyy_xx_1,  \
                             ta1_z_xxyyy_xz_0,  \
                             ta1_z_xxyyy_xz_1,  \
                             ta1_z_xxyyyy_xx_0, \
                             ta1_z_xxyyyy_xy_0, \
                             ta1_z_xxyyyy_xz_0, \
                             ta1_z_xxyyyy_yy_0, \
                             ta1_z_xxyyyy_yz_0, \
                             ta1_z_xxyyyy_zz_0, \
                             ta1_z_xyyyy_xy_0,  \
                             ta1_z_xyyyy_xy_1,  \
                             ta1_z_xyyyy_y_0,   \
                             ta1_z_xyyyy_y_1,   \
                             ta1_z_xyyyy_yy_0,  \
                             ta1_z_xyyyy_yy_1,  \
                             ta1_z_xyyyy_yz_0,  \
                             ta1_z_xyyyy_yz_1,  \
                             ta1_z_xyyyy_zz_0,  \
                             ta1_z_xyyyy_zz_1,  \
                             ta1_z_yyyy_xy_0,   \
                             ta1_z_yyyy_xy_1,   \
                             ta1_z_yyyy_yy_0,   \
                             ta1_z_yyyy_yy_1,   \
                             ta1_z_yyyy_yz_0,   \
                             ta1_z_yyyy_yz_1,   \
                             ta1_z_yyyy_zz_0,   \
                             ta1_z_yyyy_zz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyyy_xx_0[i] =
            3.0 * ta1_z_xxyy_xx_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xx_1[i] * fe_0 + ta1_z_xxyyy_xx_0[i] * pa_y[i] - ta1_z_xxyyy_xx_1[i] * pc_y[i];

        ta1_z_xxyyyy_xy_0[i] = ta1_z_yyyy_xy_0[i] * fe_0 - ta1_z_yyyy_xy_1[i] * fe_0 + ta1_z_xyyyy_y_0[i] * fe_0 - ta1_z_xyyyy_y_1[i] * fe_0 +
                               ta1_z_xyyyy_xy_0[i] * pa_x[i] - ta1_z_xyyyy_xy_1[i] * pc_x[i];

        ta1_z_xxyyyy_xz_0[i] =
            3.0 * ta1_z_xxyy_xz_0[i] * fe_0 - 3.0 * ta1_z_xxyy_xz_1[i] * fe_0 + ta1_z_xxyyy_xz_0[i] * pa_y[i] - ta1_z_xxyyy_xz_1[i] * pc_y[i];

        ta1_z_xxyyyy_yy_0[i] = ta1_z_yyyy_yy_0[i] * fe_0 - ta1_z_yyyy_yy_1[i] * fe_0 + ta1_z_xyyyy_yy_0[i] * pa_x[i] - ta1_z_xyyyy_yy_1[i] * pc_x[i];

        ta1_z_xxyyyy_yz_0[i] = ta1_z_yyyy_yz_0[i] * fe_0 - ta1_z_yyyy_yz_1[i] * fe_0 + ta1_z_xyyyy_yz_0[i] * pa_x[i] - ta1_z_xyyyy_yz_1[i] * pc_x[i];

        ta1_z_xxyyyy_zz_0[i] = ta1_z_yyyy_zz_0[i] * fe_0 - ta1_z_yyyy_zz_1[i] * fe_0 + ta1_z_xyyyy_zz_0[i] * pa_x[i] - ta1_z_xyyyy_zz_1[i] * pc_x[i];
    }

    // Set up 402-408 components of targeted buffer : ID

    auto ta1_z_xxyyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 402);

    auto ta1_z_xxyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 403);

    auto ta1_z_xxyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 404);

    auto ta1_z_xxyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 405);

    auto ta1_z_xxyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 406);

    auto ta1_z_xxyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 407);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_z_xxyyy_xx_0,  \
                             ta1_z_xxyyy_xx_1,  \
                             ta1_z_xxyyy_xy_0,  \
                             ta1_z_xxyyy_xy_1,  \
                             ta1_z_xxyyy_yy_0,  \
                             ta1_z_xxyyy_yy_1,  \
                             ta1_z_xxyyyz_xx_0, \
                             ta1_z_xxyyyz_xy_0, \
                             ta1_z_xxyyyz_xz_0, \
                             ta1_z_xxyyyz_yy_0, \
                             ta1_z_xxyyyz_yz_0, \
                             ta1_z_xxyyyz_zz_0, \
                             ta1_z_xxyyz_xz_0,  \
                             ta1_z_xxyyz_xz_1,  \
                             ta1_z_xxyz_xz_0,   \
                             ta1_z_xxyz_xz_1,   \
                             ta1_z_xyyyz_yz_0,  \
                             ta1_z_xyyyz_yz_1,  \
                             ta1_z_xyyyz_zz_0,  \
                             ta1_z_xyyyz_zz_1,  \
                             ta1_z_yyyz_yz_0,   \
                             ta1_z_yyyz_yz_1,   \
                             ta1_z_yyyz_zz_0,   \
                             ta1_z_yyyz_zz_1,   \
                             ta_xxyyy_xx_1,     \
                             ta_xxyyy_xy_1,     \
                             ta_xxyyy_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyyz_xx_0[i] = ta_xxyyy_xx_1[i] + ta1_z_xxyyy_xx_0[i] * pa_z[i] - ta1_z_xxyyy_xx_1[i] * pc_z[i];

        ta1_z_xxyyyz_xy_0[i] = ta_xxyyy_xy_1[i] + ta1_z_xxyyy_xy_0[i] * pa_z[i] - ta1_z_xxyyy_xy_1[i] * pc_z[i];

        ta1_z_xxyyyz_xz_0[i] =
            2.0 * ta1_z_xxyz_xz_0[i] * fe_0 - 2.0 * ta1_z_xxyz_xz_1[i] * fe_0 + ta1_z_xxyyz_xz_0[i] * pa_y[i] - ta1_z_xxyyz_xz_1[i] * pc_y[i];

        ta1_z_xxyyyz_yy_0[i] = ta_xxyyy_yy_1[i] + ta1_z_xxyyy_yy_0[i] * pa_z[i] - ta1_z_xxyyy_yy_1[i] * pc_z[i];

        ta1_z_xxyyyz_yz_0[i] = ta1_z_yyyz_yz_0[i] * fe_0 - ta1_z_yyyz_yz_1[i] * fe_0 + ta1_z_xyyyz_yz_0[i] * pa_x[i] - ta1_z_xyyyz_yz_1[i] * pc_x[i];

        ta1_z_xxyyyz_zz_0[i] = ta1_z_yyyz_zz_0[i] * fe_0 - ta1_z_yyyz_zz_1[i] * fe_0 + ta1_z_xyyyz_zz_0[i] * pa_x[i] - ta1_z_xyyyz_zz_1[i] * pc_x[i];
    }

    // Set up 408-414 components of targeted buffer : ID

    auto ta1_z_xxyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 408);

    auto ta1_z_xxyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 409);

    auto ta1_z_xxyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 410);

    auto ta1_z_xxyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 411);

    auto ta1_z_xxyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 412);

    auto ta1_z_xxyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 413);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pa_z,              \
                             pc_x,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_z_xxyy_xy_0,   \
                             ta1_z_xxyy_xy_1,   \
                             ta1_z_xxyyz_xy_0,  \
                             ta1_z_xxyyz_xy_1,  \
                             ta1_z_xxyyzz_xx_0, \
                             ta1_z_xxyyzz_xy_0, \
                             ta1_z_xxyyzz_xz_0, \
                             ta1_z_xxyyzz_yy_0, \
                             ta1_z_xxyyzz_yz_0, \
                             ta1_z_xxyyzz_zz_0, \
                             ta1_z_xxyzz_xx_0,  \
                             ta1_z_xxyzz_xx_1,  \
                             ta1_z_xxyzz_xz_0,  \
                             ta1_z_xxyzz_xz_1,  \
                             ta1_z_xxzz_xx_0,   \
                             ta1_z_xxzz_xx_1,   \
                             ta1_z_xxzz_xz_0,   \
                             ta1_z_xxzz_xz_1,   \
                             ta1_z_xyyzz_yy_0,  \
                             ta1_z_xyyzz_yy_1,  \
                             ta1_z_xyyzz_yz_0,  \
                             ta1_z_xyyzz_yz_1,  \
                             ta1_z_xyyzz_zz_0,  \
                             ta1_z_xyyzz_zz_1,  \
                             ta1_z_yyzz_yy_0,   \
                             ta1_z_yyzz_yy_1,   \
                             ta1_z_yyzz_yz_0,   \
                             ta1_z_yyzz_yz_1,   \
                             ta1_z_yyzz_zz_0,   \
                             ta1_z_yyzz_zz_1,   \
                             ta_xxyyz_xy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyyzz_xx_0[i] = ta1_z_xxzz_xx_0[i] * fe_0 - ta1_z_xxzz_xx_1[i] * fe_0 + ta1_z_xxyzz_xx_0[i] * pa_y[i] - ta1_z_xxyzz_xx_1[i] * pc_y[i];

        ta1_z_xxyyzz_xy_0[i] =
            ta1_z_xxyy_xy_0[i] * fe_0 - ta1_z_xxyy_xy_1[i] * fe_0 + ta_xxyyz_xy_1[i] + ta1_z_xxyyz_xy_0[i] * pa_z[i] - ta1_z_xxyyz_xy_1[i] * pc_z[i];

        ta1_z_xxyyzz_xz_0[i] = ta1_z_xxzz_xz_0[i] * fe_0 - ta1_z_xxzz_xz_1[i] * fe_0 + ta1_z_xxyzz_xz_0[i] * pa_y[i] - ta1_z_xxyzz_xz_1[i] * pc_y[i];

        ta1_z_xxyyzz_yy_0[i] = ta1_z_yyzz_yy_0[i] * fe_0 - ta1_z_yyzz_yy_1[i] * fe_0 + ta1_z_xyyzz_yy_0[i] * pa_x[i] - ta1_z_xyyzz_yy_1[i] * pc_x[i];

        ta1_z_xxyyzz_yz_0[i] = ta1_z_yyzz_yz_0[i] * fe_0 - ta1_z_yyzz_yz_1[i] * fe_0 + ta1_z_xyyzz_yz_0[i] * pa_x[i] - ta1_z_xyyzz_yz_1[i] * pc_x[i];

        ta1_z_xxyyzz_zz_0[i] = ta1_z_yyzz_zz_0[i] * fe_0 - ta1_z_yyzz_zz_1[i] * fe_0 + ta1_z_xyyzz_zz_0[i] * pa_x[i] - ta1_z_xyyzz_zz_1[i] * pc_x[i];
    }

    // Set up 414-420 components of targeted buffer : ID

    auto ta1_z_xxyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 414);

    auto ta1_z_xxyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 415);

    auto ta1_z_xxyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 416);

    auto ta1_z_xxyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 417);

    auto ta1_z_xxyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 418);

    auto ta1_z_xxyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 419);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_z_xxyzzz_xx_0, \
                             ta1_z_xxyzzz_xy_0, \
                             ta1_z_xxyzzz_xz_0, \
                             ta1_z_xxyzzz_yy_0, \
                             ta1_z_xxyzzz_yz_0, \
                             ta1_z_xxyzzz_zz_0, \
                             ta1_z_xxzzz_x_0,   \
                             ta1_z_xxzzz_x_1,   \
                             ta1_z_xxzzz_xx_0,  \
                             ta1_z_xxzzz_xx_1,  \
                             ta1_z_xxzzz_xy_0,  \
                             ta1_z_xxzzz_xy_1,  \
                             ta1_z_xxzzz_xz_0,  \
                             ta1_z_xxzzz_xz_1,  \
                             ta1_z_xxzzz_zz_0,  \
                             ta1_z_xxzzz_zz_1,  \
                             ta1_z_xyzzz_yy_0,  \
                             ta1_z_xyzzz_yy_1,  \
                             ta1_z_xyzzz_yz_0,  \
                             ta1_z_xyzzz_yz_1,  \
                             ta1_z_yzzz_yy_0,   \
                             ta1_z_yzzz_yy_1,   \
                             ta1_z_yzzz_yz_0,   \
                             ta1_z_yzzz_yz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyzzz_xx_0[i] = ta1_z_xxzzz_xx_0[i] * pa_y[i] - ta1_z_xxzzz_xx_1[i] * pc_y[i];

        ta1_z_xxyzzz_xy_0[i] = ta1_z_xxzzz_x_0[i] * fe_0 - ta1_z_xxzzz_x_1[i] * fe_0 + ta1_z_xxzzz_xy_0[i] * pa_y[i] - ta1_z_xxzzz_xy_1[i] * pc_y[i];

        ta1_z_xxyzzz_xz_0[i] = ta1_z_xxzzz_xz_0[i] * pa_y[i] - ta1_z_xxzzz_xz_1[i] * pc_y[i];

        ta1_z_xxyzzz_yy_0[i] = ta1_z_yzzz_yy_0[i] * fe_0 - ta1_z_yzzz_yy_1[i] * fe_0 + ta1_z_xyzzz_yy_0[i] * pa_x[i] - ta1_z_xyzzz_yy_1[i] * pc_x[i];

        ta1_z_xxyzzz_yz_0[i] = ta1_z_yzzz_yz_0[i] * fe_0 - ta1_z_yzzz_yz_1[i] * fe_0 + ta1_z_xyzzz_yz_0[i] * pa_x[i] - ta1_z_xyzzz_yz_1[i] * pc_x[i];

        ta1_z_xxyzzz_zz_0[i] = ta1_z_xxzzz_zz_0[i] * pa_y[i] - ta1_z_xxzzz_zz_1[i] * pc_y[i];
    }

    // Set up 420-426 components of targeted buffer : ID

    auto ta1_z_xxzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 420);

    auto ta1_z_xxzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 421);

    auto ta1_z_xxzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 422);

    auto ta1_z_xxzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 423);

    auto ta1_z_xxzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 424);

    auto ta1_z_xxzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 425);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_z_xxzz_xx_0,   \
                             ta1_z_xxzz_xx_1,   \
                             ta1_z_xxzz_xy_0,   \
                             ta1_z_xxzz_xy_1,   \
                             ta1_z_xxzzz_xx_0,  \
                             ta1_z_xxzzz_xx_1,  \
                             ta1_z_xxzzz_xy_0,  \
                             ta1_z_xxzzz_xy_1,  \
                             ta1_z_xxzzzz_xx_0, \
                             ta1_z_xxzzzz_xy_0, \
                             ta1_z_xxzzzz_xz_0, \
                             ta1_z_xxzzzz_yy_0, \
                             ta1_z_xxzzzz_yz_0, \
                             ta1_z_xxzzzz_zz_0, \
                             ta1_z_xzzzz_xz_0,  \
                             ta1_z_xzzzz_xz_1,  \
                             ta1_z_xzzzz_yy_0,  \
                             ta1_z_xzzzz_yy_1,  \
                             ta1_z_xzzzz_yz_0,  \
                             ta1_z_xzzzz_yz_1,  \
                             ta1_z_xzzzz_z_0,   \
                             ta1_z_xzzzz_z_1,   \
                             ta1_z_xzzzz_zz_0,  \
                             ta1_z_xzzzz_zz_1,  \
                             ta1_z_zzzz_xz_0,   \
                             ta1_z_zzzz_xz_1,   \
                             ta1_z_zzzz_yy_0,   \
                             ta1_z_zzzz_yy_1,   \
                             ta1_z_zzzz_yz_0,   \
                             ta1_z_zzzz_yz_1,   \
                             ta1_z_zzzz_zz_0,   \
                             ta1_z_zzzz_zz_1,   \
                             ta_xxzzz_xx_1,     \
                             ta_xxzzz_xy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxzzzz_xx_0[i] = 3.0 * ta1_z_xxzz_xx_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xx_1[i] * fe_0 + ta_xxzzz_xx_1[i] + ta1_z_xxzzz_xx_0[i] * pa_z[i] -
                               ta1_z_xxzzz_xx_1[i] * pc_z[i];

        ta1_z_xxzzzz_xy_0[i] = 3.0 * ta1_z_xxzz_xy_0[i] * fe_0 - 3.0 * ta1_z_xxzz_xy_1[i] * fe_0 + ta_xxzzz_xy_1[i] + ta1_z_xxzzz_xy_0[i] * pa_z[i] -
                               ta1_z_xxzzz_xy_1[i] * pc_z[i];

        ta1_z_xxzzzz_xz_0[i] = ta1_z_zzzz_xz_0[i] * fe_0 - ta1_z_zzzz_xz_1[i] * fe_0 + ta1_z_xzzzz_z_0[i] * fe_0 - ta1_z_xzzzz_z_1[i] * fe_0 +
                               ta1_z_xzzzz_xz_0[i] * pa_x[i] - ta1_z_xzzzz_xz_1[i] * pc_x[i];

        ta1_z_xxzzzz_yy_0[i] = ta1_z_zzzz_yy_0[i] * fe_0 - ta1_z_zzzz_yy_1[i] * fe_0 + ta1_z_xzzzz_yy_0[i] * pa_x[i] - ta1_z_xzzzz_yy_1[i] * pc_x[i];

        ta1_z_xxzzzz_yz_0[i] = ta1_z_zzzz_yz_0[i] * fe_0 - ta1_z_zzzz_yz_1[i] * fe_0 + ta1_z_xzzzz_yz_0[i] * pa_x[i] - ta1_z_xzzzz_yz_1[i] * pc_x[i];

        ta1_z_xxzzzz_zz_0[i] = ta1_z_zzzz_zz_0[i] * fe_0 - ta1_z_zzzz_zz_1[i] * fe_0 + ta1_z_xzzzz_zz_0[i] * pa_x[i] - ta1_z_xzzzz_zz_1[i] * pc_x[i];
    }

    // Set up 426-432 components of targeted buffer : ID

    auto ta1_z_xyyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 426);

    auto ta1_z_xyyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 427);

    auto ta1_z_xyyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 428);

    auto ta1_z_xyyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 429);

    auto ta1_z_xyyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 430);

    auto ta1_z_xyyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 431);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_z_xyyyyy_xx_0, \
                             ta1_z_xyyyyy_xy_0, \
                             ta1_z_xyyyyy_xz_0, \
                             ta1_z_xyyyyy_yy_0, \
                             ta1_z_xyyyyy_yz_0, \
                             ta1_z_xyyyyy_zz_0, \
                             ta1_z_yyyyy_x_0,   \
                             ta1_z_yyyyy_x_1,   \
                             ta1_z_yyyyy_xx_0,  \
                             ta1_z_yyyyy_xx_1,  \
                             ta1_z_yyyyy_xy_0,  \
                             ta1_z_yyyyy_xy_1,  \
                             ta1_z_yyyyy_xz_0,  \
                             ta1_z_yyyyy_xz_1,  \
                             ta1_z_yyyyy_y_0,   \
                             ta1_z_yyyyy_y_1,   \
                             ta1_z_yyyyy_yy_0,  \
                             ta1_z_yyyyy_yy_1,  \
                             ta1_z_yyyyy_yz_0,  \
                             ta1_z_yyyyy_yz_1,  \
                             ta1_z_yyyyy_z_0,   \
                             ta1_z_yyyyy_z_1,   \
                             ta1_z_yyyyy_zz_0,  \
                             ta1_z_yyyyy_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyyy_xx_0[i] =
            2.0 * ta1_z_yyyyy_x_0[i] * fe_0 - 2.0 * ta1_z_yyyyy_x_1[i] * fe_0 + ta1_z_yyyyy_xx_0[i] * pa_x[i] - ta1_z_yyyyy_xx_1[i] * pc_x[i];

        ta1_z_xyyyyy_xy_0[i] = ta1_z_yyyyy_y_0[i] * fe_0 - ta1_z_yyyyy_y_1[i] * fe_0 + ta1_z_yyyyy_xy_0[i] * pa_x[i] - ta1_z_yyyyy_xy_1[i] * pc_x[i];

        ta1_z_xyyyyy_xz_0[i] = ta1_z_yyyyy_z_0[i] * fe_0 - ta1_z_yyyyy_z_1[i] * fe_0 + ta1_z_yyyyy_xz_0[i] * pa_x[i] - ta1_z_yyyyy_xz_1[i] * pc_x[i];

        ta1_z_xyyyyy_yy_0[i] = ta1_z_yyyyy_yy_0[i] * pa_x[i] - ta1_z_yyyyy_yy_1[i] * pc_x[i];

        ta1_z_xyyyyy_yz_0[i] = ta1_z_yyyyy_yz_0[i] * pa_x[i] - ta1_z_yyyyy_yz_1[i] * pc_x[i];

        ta1_z_xyyyyy_zz_0[i] = ta1_z_yyyyy_zz_0[i] * pa_x[i] - ta1_z_yyyyy_zz_1[i] * pc_x[i];
    }

    // Set up 432-438 components of targeted buffer : ID

    auto ta1_z_xyyyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 432);

    auto ta1_z_xyyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 433);

    auto ta1_z_xyyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 434);

    auto ta1_z_xyyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 435);

    auto ta1_z_xyyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 436);

    auto ta1_z_xyyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 437);

#pragma omp simd aligned(pa_x,                  \
                             pa_z,              \
                             pc_x,              \
                             pc_z,              \
                             ta1_z_xyyyy_xx_0,  \
                             ta1_z_xyyyy_xx_1,  \
                             ta1_z_xyyyy_xy_0,  \
                             ta1_z_xyyyy_xy_1,  \
                             ta1_z_xyyyyz_xx_0, \
                             ta1_z_xyyyyz_xy_0, \
                             ta1_z_xyyyyz_xz_0, \
                             ta1_z_xyyyyz_yy_0, \
                             ta1_z_xyyyyz_yz_0, \
                             ta1_z_xyyyyz_zz_0, \
                             ta1_z_yyyyz_xz_0,  \
                             ta1_z_yyyyz_xz_1,  \
                             ta1_z_yyyyz_yy_0,  \
                             ta1_z_yyyyz_yy_1,  \
                             ta1_z_yyyyz_yz_0,  \
                             ta1_z_yyyyz_yz_1,  \
                             ta1_z_yyyyz_z_0,   \
                             ta1_z_yyyyz_z_1,   \
                             ta1_z_yyyyz_zz_0,  \
                             ta1_z_yyyyz_zz_1,  \
                             ta_xyyyy_xx_1,     \
                             ta_xyyyy_xy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyyz_xx_0[i] = ta_xyyyy_xx_1[i] + ta1_z_xyyyy_xx_0[i] * pa_z[i] - ta1_z_xyyyy_xx_1[i] * pc_z[i];

        ta1_z_xyyyyz_xy_0[i] = ta_xyyyy_xy_1[i] + ta1_z_xyyyy_xy_0[i] * pa_z[i] - ta1_z_xyyyy_xy_1[i] * pc_z[i];

        ta1_z_xyyyyz_xz_0[i] = ta1_z_yyyyz_z_0[i] * fe_0 - ta1_z_yyyyz_z_1[i] * fe_0 + ta1_z_yyyyz_xz_0[i] * pa_x[i] - ta1_z_yyyyz_xz_1[i] * pc_x[i];

        ta1_z_xyyyyz_yy_0[i] = ta1_z_yyyyz_yy_0[i] * pa_x[i] - ta1_z_yyyyz_yy_1[i] * pc_x[i];

        ta1_z_xyyyyz_yz_0[i] = ta1_z_yyyyz_yz_0[i] * pa_x[i] - ta1_z_yyyyz_yz_1[i] * pc_x[i];

        ta1_z_xyyyyz_zz_0[i] = ta1_z_yyyyz_zz_0[i] * pa_x[i] - ta1_z_yyyyz_zz_1[i] * pc_x[i];
    }

    // Set up 438-444 components of targeted buffer : ID

    auto ta1_z_xyyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 438);

    auto ta1_z_xyyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 439);

    auto ta1_z_xyyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 440);

    auto ta1_z_xyyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 441);

    auto ta1_z_xyyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 442);

    auto ta1_z_xyyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 443);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_z_xyyyzz_xx_0, \
                             ta1_z_xyyyzz_xy_0, \
                             ta1_z_xyyyzz_xz_0, \
                             ta1_z_xyyyzz_yy_0, \
                             ta1_z_xyyyzz_yz_0, \
                             ta1_z_xyyyzz_zz_0, \
                             ta1_z_yyyzz_x_0,   \
                             ta1_z_yyyzz_x_1,   \
                             ta1_z_yyyzz_xx_0,  \
                             ta1_z_yyyzz_xx_1,  \
                             ta1_z_yyyzz_xy_0,  \
                             ta1_z_yyyzz_xy_1,  \
                             ta1_z_yyyzz_xz_0,  \
                             ta1_z_yyyzz_xz_1,  \
                             ta1_z_yyyzz_y_0,   \
                             ta1_z_yyyzz_y_1,   \
                             ta1_z_yyyzz_yy_0,  \
                             ta1_z_yyyzz_yy_1,  \
                             ta1_z_yyyzz_yz_0,  \
                             ta1_z_yyyzz_yz_1,  \
                             ta1_z_yyyzz_z_0,   \
                             ta1_z_yyyzz_z_1,   \
                             ta1_z_yyyzz_zz_0,  \
                             ta1_z_yyyzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyyzz_xx_0[i] =
            2.0 * ta1_z_yyyzz_x_0[i] * fe_0 - 2.0 * ta1_z_yyyzz_x_1[i] * fe_0 + ta1_z_yyyzz_xx_0[i] * pa_x[i] - ta1_z_yyyzz_xx_1[i] * pc_x[i];

        ta1_z_xyyyzz_xy_0[i] = ta1_z_yyyzz_y_0[i] * fe_0 - ta1_z_yyyzz_y_1[i] * fe_0 + ta1_z_yyyzz_xy_0[i] * pa_x[i] - ta1_z_yyyzz_xy_1[i] * pc_x[i];

        ta1_z_xyyyzz_xz_0[i] = ta1_z_yyyzz_z_0[i] * fe_0 - ta1_z_yyyzz_z_1[i] * fe_0 + ta1_z_yyyzz_xz_0[i] * pa_x[i] - ta1_z_yyyzz_xz_1[i] * pc_x[i];

        ta1_z_xyyyzz_yy_0[i] = ta1_z_yyyzz_yy_0[i] * pa_x[i] - ta1_z_yyyzz_yy_1[i] * pc_x[i];

        ta1_z_xyyyzz_yz_0[i] = ta1_z_yyyzz_yz_0[i] * pa_x[i] - ta1_z_yyyzz_yz_1[i] * pc_x[i];

        ta1_z_xyyyzz_zz_0[i] = ta1_z_yyyzz_zz_0[i] * pa_x[i] - ta1_z_yyyzz_zz_1[i] * pc_x[i];
    }

    // Set up 444-450 components of targeted buffer : ID

    auto ta1_z_xyyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 444);

    auto ta1_z_xyyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 445);

    auto ta1_z_xyyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 446);

    auto ta1_z_xyyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 447);

    auto ta1_z_xyyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 448);

    auto ta1_z_xyyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 449);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_z_xyyzzz_xx_0, \
                             ta1_z_xyyzzz_xy_0, \
                             ta1_z_xyyzzz_xz_0, \
                             ta1_z_xyyzzz_yy_0, \
                             ta1_z_xyyzzz_yz_0, \
                             ta1_z_xyyzzz_zz_0, \
                             ta1_z_yyzzz_x_0,   \
                             ta1_z_yyzzz_x_1,   \
                             ta1_z_yyzzz_xx_0,  \
                             ta1_z_yyzzz_xx_1,  \
                             ta1_z_yyzzz_xy_0,  \
                             ta1_z_yyzzz_xy_1,  \
                             ta1_z_yyzzz_xz_0,  \
                             ta1_z_yyzzz_xz_1,  \
                             ta1_z_yyzzz_y_0,   \
                             ta1_z_yyzzz_y_1,   \
                             ta1_z_yyzzz_yy_0,  \
                             ta1_z_yyzzz_yy_1,  \
                             ta1_z_yyzzz_yz_0,  \
                             ta1_z_yyzzz_yz_1,  \
                             ta1_z_yyzzz_z_0,   \
                             ta1_z_yyzzz_z_1,   \
                             ta1_z_yyzzz_zz_0,  \
                             ta1_z_yyzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyzzz_xx_0[i] =
            2.0 * ta1_z_yyzzz_x_0[i] * fe_0 - 2.0 * ta1_z_yyzzz_x_1[i] * fe_0 + ta1_z_yyzzz_xx_0[i] * pa_x[i] - ta1_z_yyzzz_xx_1[i] * pc_x[i];

        ta1_z_xyyzzz_xy_0[i] = ta1_z_yyzzz_y_0[i] * fe_0 - ta1_z_yyzzz_y_1[i] * fe_0 + ta1_z_yyzzz_xy_0[i] * pa_x[i] - ta1_z_yyzzz_xy_1[i] * pc_x[i];

        ta1_z_xyyzzz_xz_0[i] = ta1_z_yyzzz_z_0[i] * fe_0 - ta1_z_yyzzz_z_1[i] * fe_0 + ta1_z_yyzzz_xz_0[i] * pa_x[i] - ta1_z_yyzzz_xz_1[i] * pc_x[i];

        ta1_z_xyyzzz_yy_0[i] = ta1_z_yyzzz_yy_0[i] * pa_x[i] - ta1_z_yyzzz_yy_1[i] * pc_x[i];

        ta1_z_xyyzzz_yz_0[i] = ta1_z_yyzzz_yz_0[i] * pa_x[i] - ta1_z_yyzzz_yz_1[i] * pc_x[i];

        ta1_z_xyyzzz_zz_0[i] = ta1_z_yyzzz_zz_0[i] * pa_x[i] - ta1_z_yyzzz_zz_1[i] * pc_x[i];
    }

    // Set up 450-456 components of targeted buffer : ID

    auto ta1_z_xyzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 450);

    auto ta1_z_xyzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 451);

    auto ta1_z_xyzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 452);

    auto ta1_z_xyzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 453);

    auto ta1_z_xyzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 454);

    auto ta1_z_xyzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 455);

#pragma omp simd aligned(pa_x,                  \
                             pa_y,              \
                             pc_x,              \
                             pc_y,              \
                             ta1_z_xyzzzz_xx_0, \
                             ta1_z_xyzzzz_xy_0, \
                             ta1_z_xyzzzz_xz_0, \
                             ta1_z_xyzzzz_yy_0, \
                             ta1_z_xyzzzz_yz_0, \
                             ta1_z_xyzzzz_zz_0, \
                             ta1_z_xzzzz_xx_0,  \
                             ta1_z_xzzzz_xx_1,  \
                             ta1_z_xzzzz_xz_0,  \
                             ta1_z_xzzzz_xz_1,  \
                             ta1_z_yzzzz_xy_0,  \
                             ta1_z_yzzzz_xy_1,  \
                             ta1_z_yzzzz_y_0,   \
                             ta1_z_yzzzz_y_1,   \
                             ta1_z_yzzzz_yy_0,  \
                             ta1_z_yzzzz_yy_1,  \
                             ta1_z_yzzzz_yz_0,  \
                             ta1_z_yzzzz_yz_1,  \
                             ta1_z_yzzzz_zz_0,  \
                             ta1_z_yzzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyzzzz_xx_0[i] = ta1_z_xzzzz_xx_0[i] * pa_y[i] - ta1_z_xzzzz_xx_1[i] * pc_y[i];

        ta1_z_xyzzzz_xy_0[i] = ta1_z_yzzzz_y_0[i] * fe_0 - ta1_z_yzzzz_y_1[i] * fe_0 + ta1_z_yzzzz_xy_0[i] * pa_x[i] - ta1_z_yzzzz_xy_1[i] * pc_x[i];

        ta1_z_xyzzzz_xz_0[i] = ta1_z_xzzzz_xz_0[i] * pa_y[i] - ta1_z_xzzzz_xz_1[i] * pc_y[i];

        ta1_z_xyzzzz_yy_0[i] = ta1_z_yzzzz_yy_0[i] * pa_x[i] - ta1_z_yzzzz_yy_1[i] * pc_x[i];

        ta1_z_xyzzzz_yz_0[i] = ta1_z_yzzzz_yz_0[i] * pa_x[i] - ta1_z_yzzzz_yz_1[i] * pc_x[i];

        ta1_z_xyzzzz_zz_0[i] = ta1_z_yzzzz_zz_0[i] * pa_x[i] - ta1_z_yzzzz_zz_1[i] * pc_x[i];
    }

    // Set up 456-462 components of targeted buffer : ID

    auto ta1_z_xzzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 456);

    auto ta1_z_xzzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 457);

    auto ta1_z_xzzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 458);

    auto ta1_z_xzzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 459);

    auto ta1_z_xzzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 460);

    auto ta1_z_xzzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 461);

#pragma omp simd aligned(pa_x,                  \
                             pc_x,              \
                             ta1_z_xzzzzz_xx_0, \
                             ta1_z_xzzzzz_xy_0, \
                             ta1_z_xzzzzz_xz_0, \
                             ta1_z_xzzzzz_yy_0, \
                             ta1_z_xzzzzz_yz_0, \
                             ta1_z_xzzzzz_zz_0, \
                             ta1_z_zzzzz_x_0,   \
                             ta1_z_zzzzz_x_1,   \
                             ta1_z_zzzzz_xx_0,  \
                             ta1_z_zzzzz_xx_1,  \
                             ta1_z_zzzzz_xy_0,  \
                             ta1_z_zzzzz_xy_1,  \
                             ta1_z_zzzzz_xz_0,  \
                             ta1_z_zzzzz_xz_1,  \
                             ta1_z_zzzzz_y_0,   \
                             ta1_z_zzzzz_y_1,   \
                             ta1_z_zzzzz_yy_0,  \
                             ta1_z_zzzzz_yy_1,  \
                             ta1_z_zzzzz_yz_0,  \
                             ta1_z_zzzzz_yz_1,  \
                             ta1_z_zzzzz_z_0,   \
                             ta1_z_zzzzz_z_1,   \
                             ta1_z_zzzzz_zz_0,  \
                             ta1_z_zzzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzzzzz_xx_0[i] =
            2.0 * ta1_z_zzzzz_x_0[i] * fe_0 - 2.0 * ta1_z_zzzzz_x_1[i] * fe_0 + ta1_z_zzzzz_xx_0[i] * pa_x[i] - ta1_z_zzzzz_xx_1[i] * pc_x[i];

        ta1_z_xzzzzz_xy_0[i] = ta1_z_zzzzz_y_0[i] * fe_0 - ta1_z_zzzzz_y_1[i] * fe_0 + ta1_z_zzzzz_xy_0[i] * pa_x[i] - ta1_z_zzzzz_xy_1[i] * pc_x[i];

        ta1_z_xzzzzz_xz_0[i] = ta1_z_zzzzz_z_0[i] * fe_0 - ta1_z_zzzzz_z_1[i] * fe_0 + ta1_z_zzzzz_xz_0[i] * pa_x[i] - ta1_z_zzzzz_xz_1[i] * pc_x[i];

        ta1_z_xzzzzz_yy_0[i] = ta1_z_zzzzz_yy_0[i] * pa_x[i] - ta1_z_zzzzz_yy_1[i] * pc_x[i];

        ta1_z_xzzzzz_yz_0[i] = ta1_z_zzzzz_yz_0[i] * pa_x[i] - ta1_z_zzzzz_yz_1[i] * pc_x[i];

        ta1_z_xzzzzz_zz_0[i] = ta1_z_zzzzz_zz_0[i] * pa_x[i] - ta1_z_zzzzz_zz_1[i] * pc_x[i];
    }

    // Set up 462-468 components of targeted buffer : ID

    auto ta1_z_yyyyyy_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 462);

    auto ta1_z_yyyyyy_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 463);

    auto ta1_z_yyyyyy_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 464);

    auto ta1_z_yyyyyy_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 465);

    auto ta1_z_yyyyyy_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 466);

    auto ta1_z_yyyyyy_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 467);

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_z_yyyy_xx_0,   \
                             ta1_z_yyyy_xx_1,   \
                             ta1_z_yyyy_xy_0,   \
                             ta1_z_yyyy_xy_1,   \
                             ta1_z_yyyy_xz_0,   \
                             ta1_z_yyyy_xz_1,   \
                             ta1_z_yyyy_yy_0,   \
                             ta1_z_yyyy_yy_1,   \
                             ta1_z_yyyy_yz_0,   \
                             ta1_z_yyyy_yz_1,   \
                             ta1_z_yyyy_zz_0,   \
                             ta1_z_yyyy_zz_1,   \
                             ta1_z_yyyyy_x_0,   \
                             ta1_z_yyyyy_x_1,   \
                             ta1_z_yyyyy_xx_0,  \
                             ta1_z_yyyyy_xx_1,  \
                             ta1_z_yyyyy_xy_0,  \
                             ta1_z_yyyyy_xy_1,  \
                             ta1_z_yyyyy_xz_0,  \
                             ta1_z_yyyyy_xz_1,  \
                             ta1_z_yyyyy_y_0,   \
                             ta1_z_yyyyy_y_1,   \
                             ta1_z_yyyyy_yy_0,  \
                             ta1_z_yyyyy_yy_1,  \
                             ta1_z_yyyyy_yz_0,  \
                             ta1_z_yyyyy_yz_1,  \
                             ta1_z_yyyyy_z_0,   \
                             ta1_z_yyyyy_z_1,   \
                             ta1_z_yyyyy_zz_0,  \
                             ta1_z_yyyyy_zz_1,  \
                             ta1_z_yyyyyy_xx_0, \
                             ta1_z_yyyyyy_xy_0, \
                             ta1_z_yyyyyy_xz_0, \
                             ta1_z_yyyyyy_yy_0, \
                             ta1_z_yyyyyy_yz_0, \
                             ta1_z_yyyyyy_zz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyyy_xx_0[i] =
            5.0 * ta1_z_yyyy_xx_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xx_1[i] * fe_0 + ta1_z_yyyyy_xx_0[i] * pa_y[i] - ta1_z_yyyyy_xx_1[i] * pc_y[i];

        ta1_z_yyyyyy_xy_0[i] = 5.0 * ta1_z_yyyy_xy_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xy_1[i] * fe_0 + ta1_z_yyyyy_x_0[i] * fe_0 -
                               ta1_z_yyyyy_x_1[i] * fe_0 + ta1_z_yyyyy_xy_0[i] * pa_y[i] - ta1_z_yyyyy_xy_1[i] * pc_y[i];

        ta1_z_yyyyyy_xz_0[i] =
            5.0 * ta1_z_yyyy_xz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_xz_1[i] * fe_0 + ta1_z_yyyyy_xz_0[i] * pa_y[i] - ta1_z_yyyyy_xz_1[i] * pc_y[i];

        ta1_z_yyyyyy_yy_0[i] = 5.0 * ta1_z_yyyy_yy_0[i] * fe_0 - 5.0 * ta1_z_yyyy_yy_1[i] * fe_0 + 2.0 * ta1_z_yyyyy_y_0[i] * fe_0 -
                               2.0 * ta1_z_yyyyy_y_1[i] * fe_0 + ta1_z_yyyyy_yy_0[i] * pa_y[i] - ta1_z_yyyyy_yy_1[i] * pc_y[i];

        ta1_z_yyyyyy_yz_0[i] = 5.0 * ta1_z_yyyy_yz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_yz_1[i] * fe_0 + ta1_z_yyyyy_z_0[i] * fe_0 -
                               ta1_z_yyyyy_z_1[i] * fe_0 + ta1_z_yyyyy_yz_0[i] * pa_y[i] - ta1_z_yyyyy_yz_1[i] * pc_y[i];

        ta1_z_yyyyyy_zz_0[i] =
            5.0 * ta1_z_yyyy_zz_0[i] * fe_0 - 5.0 * ta1_z_yyyy_zz_1[i] * fe_0 + ta1_z_yyyyy_zz_0[i] * pa_y[i] - ta1_z_yyyyy_zz_1[i] * pc_y[i];
    }

    // Set up 468-474 components of targeted buffer : ID

    auto ta1_z_yyyyyz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 468);

    auto ta1_z_yyyyyz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 469);

    auto ta1_z_yyyyyz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 470);

    auto ta1_z_yyyyyz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 471);

    auto ta1_z_yyyyyz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 472);

    auto ta1_z_yyyyyz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 473);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_z_yyyyy_xx_0,  \
                             ta1_z_yyyyy_xx_1,  \
                             ta1_z_yyyyy_xy_0,  \
                             ta1_z_yyyyy_xy_1,  \
                             ta1_z_yyyyy_y_0,   \
                             ta1_z_yyyyy_y_1,   \
                             ta1_z_yyyyy_yy_0,  \
                             ta1_z_yyyyy_yy_1,  \
                             ta1_z_yyyyy_yz_0,  \
                             ta1_z_yyyyy_yz_1,  \
                             ta1_z_yyyyyz_xx_0, \
                             ta1_z_yyyyyz_xy_0, \
                             ta1_z_yyyyyz_xz_0, \
                             ta1_z_yyyyyz_yy_0, \
                             ta1_z_yyyyyz_yz_0, \
                             ta1_z_yyyyyz_zz_0, \
                             ta1_z_yyyyz_xz_0,  \
                             ta1_z_yyyyz_xz_1,  \
                             ta1_z_yyyyz_zz_0,  \
                             ta1_z_yyyyz_zz_1,  \
                             ta1_z_yyyz_xz_0,   \
                             ta1_z_yyyz_xz_1,   \
                             ta1_z_yyyz_zz_0,   \
                             ta1_z_yyyz_zz_1,   \
                             ta_yyyyy_xx_1,     \
                             ta_yyyyy_xy_1,     \
                             ta_yyyyy_yy_1,     \
                             ta_yyyyy_yz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyyz_xx_0[i] = ta_yyyyy_xx_1[i] + ta1_z_yyyyy_xx_0[i] * pa_z[i] - ta1_z_yyyyy_xx_1[i] * pc_z[i];

        ta1_z_yyyyyz_xy_0[i] = ta_yyyyy_xy_1[i] + ta1_z_yyyyy_xy_0[i] * pa_z[i] - ta1_z_yyyyy_xy_1[i] * pc_z[i];

        ta1_z_yyyyyz_xz_0[i] =
            4.0 * ta1_z_yyyz_xz_0[i] * fe_0 - 4.0 * ta1_z_yyyz_xz_1[i] * fe_0 + ta1_z_yyyyz_xz_0[i] * pa_y[i] - ta1_z_yyyyz_xz_1[i] * pc_y[i];

        ta1_z_yyyyyz_yy_0[i] = ta_yyyyy_yy_1[i] + ta1_z_yyyyy_yy_0[i] * pa_z[i] - ta1_z_yyyyy_yy_1[i] * pc_z[i];

        ta1_z_yyyyyz_yz_0[i] =
            ta1_z_yyyyy_y_0[i] * fe_0 - ta1_z_yyyyy_y_1[i] * fe_0 + ta_yyyyy_yz_1[i] + ta1_z_yyyyy_yz_0[i] * pa_z[i] - ta1_z_yyyyy_yz_1[i] * pc_z[i];

        ta1_z_yyyyyz_zz_0[i] =
            4.0 * ta1_z_yyyz_zz_0[i] * fe_0 - 4.0 * ta1_z_yyyz_zz_1[i] * fe_0 + ta1_z_yyyyz_zz_0[i] * pa_y[i] - ta1_z_yyyyz_zz_1[i] * pc_y[i];
    }

    // Set up 474-480 components of targeted buffer : ID

    auto ta1_z_yyyyzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 474);

    auto ta1_z_yyyyzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 475);

    auto ta1_z_yyyyzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 476);

    auto ta1_z_yyyyzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 477);

    auto ta1_z_yyyyzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 478);

    auto ta1_z_yyyyzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 479);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_z_yyyy_xy_0,   \
                             ta1_z_yyyy_xy_1,   \
                             ta1_z_yyyy_yy_0,   \
                             ta1_z_yyyy_yy_1,   \
                             ta1_z_yyyyz_xy_0,  \
                             ta1_z_yyyyz_xy_1,  \
                             ta1_z_yyyyz_yy_0,  \
                             ta1_z_yyyyz_yy_1,  \
                             ta1_z_yyyyzz_xx_0, \
                             ta1_z_yyyyzz_xy_0, \
                             ta1_z_yyyyzz_xz_0, \
                             ta1_z_yyyyzz_yy_0, \
                             ta1_z_yyyyzz_yz_0, \
                             ta1_z_yyyyzz_zz_0, \
                             ta1_z_yyyzz_xx_0,  \
                             ta1_z_yyyzz_xx_1,  \
                             ta1_z_yyyzz_xz_0,  \
                             ta1_z_yyyzz_xz_1,  \
                             ta1_z_yyyzz_yz_0,  \
                             ta1_z_yyyzz_yz_1,  \
                             ta1_z_yyyzz_z_0,   \
                             ta1_z_yyyzz_z_1,   \
                             ta1_z_yyyzz_zz_0,  \
                             ta1_z_yyyzz_zz_1,  \
                             ta1_z_yyzz_xx_0,   \
                             ta1_z_yyzz_xx_1,   \
                             ta1_z_yyzz_xz_0,   \
                             ta1_z_yyzz_xz_1,   \
                             ta1_z_yyzz_yz_0,   \
                             ta1_z_yyzz_yz_1,   \
                             ta1_z_yyzz_zz_0,   \
                             ta1_z_yyzz_zz_1,   \
                             ta_yyyyz_xy_1,     \
                             ta_yyyyz_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyyzz_xx_0[i] =
            3.0 * ta1_z_yyzz_xx_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xx_1[i] * fe_0 + ta1_z_yyyzz_xx_0[i] * pa_y[i] - ta1_z_yyyzz_xx_1[i] * pc_y[i];

        ta1_z_yyyyzz_xy_0[i] =
            ta1_z_yyyy_xy_0[i] * fe_0 - ta1_z_yyyy_xy_1[i] * fe_0 + ta_yyyyz_xy_1[i] + ta1_z_yyyyz_xy_0[i] * pa_z[i] - ta1_z_yyyyz_xy_1[i] * pc_z[i];

        ta1_z_yyyyzz_xz_0[i] =
            3.0 * ta1_z_yyzz_xz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xz_1[i] * fe_0 + ta1_z_yyyzz_xz_0[i] * pa_y[i] - ta1_z_yyyzz_xz_1[i] * pc_y[i];

        ta1_z_yyyyzz_yy_0[i] =
            ta1_z_yyyy_yy_0[i] * fe_0 - ta1_z_yyyy_yy_1[i] * fe_0 + ta_yyyyz_yy_1[i] + ta1_z_yyyyz_yy_0[i] * pa_z[i] - ta1_z_yyyyz_yy_1[i] * pc_z[i];

        ta1_z_yyyyzz_yz_0[i] = 3.0 * ta1_z_yyzz_yz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_yz_1[i] * fe_0 + ta1_z_yyyzz_z_0[i] * fe_0 -
                               ta1_z_yyyzz_z_1[i] * fe_0 + ta1_z_yyyzz_yz_0[i] * pa_y[i] - ta1_z_yyyzz_yz_1[i] * pc_y[i];

        ta1_z_yyyyzz_zz_0[i] =
            3.0 * ta1_z_yyzz_zz_0[i] * fe_0 - 3.0 * ta1_z_yyzz_zz_1[i] * fe_0 + ta1_z_yyyzz_zz_0[i] * pa_y[i] - ta1_z_yyyzz_zz_1[i] * pc_y[i];
    }

    // Set up 480-486 components of targeted buffer : ID

    auto ta1_z_yyyzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 480);

    auto ta1_z_yyyzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 481);

    auto ta1_z_yyyzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 482);

    auto ta1_z_yyyzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 483);

    auto ta1_z_yyyzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 484);

    auto ta1_z_yyyzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 485);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_z_yyyz_xy_0,   \
                             ta1_z_yyyz_xy_1,   \
                             ta1_z_yyyz_yy_0,   \
                             ta1_z_yyyz_yy_1,   \
                             ta1_z_yyyzz_xy_0,  \
                             ta1_z_yyyzz_xy_1,  \
                             ta1_z_yyyzz_yy_0,  \
                             ta1_z_yyyzz_yy_1,  \
                             ta1_z_yyyzzz_xx_0, \
                             ta1_z_yyyzzz_xy_0, \
                             ta1_z_yyyzzz_xz_0, \
                             ta1_z_yyyzzz_yy_0, \
                             ta1_z_yyyzzz_yz_0, \
                             ta1_z_yyyzzz_zz_0, \
                             ta1_z_yyzzz_xx_0,  \
                             ta1_z_yyzzz_xx_1,  \
                             ta1_z_yyzzz_xz_0,  \
                             ta1_z_yyzzz_xz_1,  \
                             ta1_z_yyzzz_yz_0,  \
                             ta1_z_yyzzz_yz_1,  \
                             ta1_z_yyzzz_z_0,   \
                             ta1_z_yyzzz_z_1,   \
                             ta1_z_yyzzz_zz_0,  \
                             ta1_z_yyzzz_zz_1,  \
                             ta1_z_yzzz_xx_0,   \
                             ta1_z_yzzz_xx_1,   \
                             ta1_z_yzzz_xz_0,   \
                             ta1_z_yzzz_xz_1,   \
                             ta1_z_yzzz_yz_0,   \
                             ta1_z_yzzz_yz_1,   \
                             ta1_z_yzzz_zz_0,   \
                             ta1_z_yzzz_zz_1,   \
                             ta_yyyzz_xy_1,     \
                             ta_yyyzz_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyzzz_xx_0[i] =
            2.0 * ta1_z_yzzz_xx_0[i] * fe_0 - 2.0 * ta1_z_yzzz_xx_1[i] * fe_0 + ta1_z_yyzzz_xx_0[i] * pa_y[i] - ta1_z_yyzzz_xx_1[i] * pc_y[i];

        ta1_z_yyyzzz_xy_0[i] = 2.0 * ta1_z_yyyz_xy_0[i] * fe_0 - 2.0 * ta1_z_yyyz_xy_1[i] * fe_0 + ta_yyyzz_xy_1[i] + ta1_z_yyyzz_xy_0[i] * pa_z[i] -
                               ta1_z_yyyzz_xy_1[i] * pc_z[i];

        ta1_z_yyyzzz_xz_0[i] =
            2.0 * ta1_z_yzzz_xz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_xz_1[i] * fe_0 + ta1_z_yyzzz_xz_0[i] * pa_y[i] - ta1_z_yyzzz_xz_1[i] * pc_y[i];

        ta1_z_yyyzzz_yy_0[i] = 2.0 * ta1_z_yyyz_yy_0[i] * fe_0 - 2.0 * ta1_z_yyyz_yy_1[i] * fe_0 + ta_yyyzz_yy_1[i] + ta1_z_yyyzz_yy_0[i] * pa_z[i] -
                               ta1_z_yyyzz_yy_1[i] * pc_z[i];

        ta1_z_yyyzzz_yz_0[i] = 2.0 * ta1_z_yzzz_yz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_yz_1[i] * fe_0 + ta1_z_yyzzz_z_0[i] * fe_0 -
                               ta1_z_yyzzz_z_1[i] * fe_0 + ta1_z_yyzzz_yz_0[i] * pa_y[i] - ta1_z_yyzzz_yz_1[i] * pc_y[i];

        ta1_z_yyyzzz_zz_0[i] =
            2.0 * ta1_z_yzzz_zz_0[i] * fe_0 - 2.0 * ta1_z_yzzz_zz_1[i] * fe_0 + ta1_z_yyzzz_zz_0[i] * pa_y[i] - ta1_z_yyzzz_zz_1[i] * pc_y[i];
    }

    // Set up 486-492 components of targeted buffer : ID

    auto ta1_z_yyzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 486);

    auto ta1_z_yyzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 487);

    auto ta1_z_yyzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 488);

    auto ta1_z_yyzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 489);

    auto ta1_z_yyzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 490);

    auto ta1_z_yyzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 491);

#pragma omp simd aligned(pa_y,                  \
                             pa_z,              \
                             pc_y,              \
                             pc_z,              \
                             ta1_z_yyzz_xy_0,   \
                             ta1_z_yyzz_xy_1,   \
                             ta1_z_yyzz_yy_0,   \
                             ta1_z_yyzz_yy_1,   \
                             ta1_z_yyzzz_xy_0,  \
                             ta1_z_yyzzz_xy_1,  \
                             ta1_z_yyzzz_yy_0,  \
                             ta1_z_yyzzz_yy_1,  \
                             ta1_z_yyzzzz_xx_0, \
                             ta1_z_yyzzzz_xy_0, \
                             ta1_z_yyzzzz_xz_0, \
                             ta1_z_yyzzzz_yy_0, \
                             ta1_z_yyzzzz_yz_0, \
                             ta1_z_yyzzzz_zz_0, \
                             ta1_z_yzzzz_xx_0,  \
                             ta1_z_yzzzz_xx_1,  \
                             ta1_z_yzzzz_xz_0,  \
                             ta1_z_yzzzz_xz_1,  \
                             ta1_z_yzzzz_yz_0,  \
                             ta1_z_yzzzz_yz_1,  \
                             ta1_z_yzzzz_z_0,   \
                             ta1_z_yzzzz_z_1,   \
                             ta1_z_yzzzz_zz_0,  \
                             ta1_z_yzzzz_zz_1,  \
                             ta1_z_zzzz_xx_0,   \
                             ta1_z_zzzz_xx_1,   \
                             ta1_z_zzzz_xz_0,   \
                             ta1_z_zzzz_xz_1,   \
                             ta1_z_zzzz_yz_0,   \
                             ta1_z_zzzz_yz_1,   \
                             ta1_z_zzzz_zz_0,   \
                             ta1_z_zzzz_zz_1,   \
                             ta_yyzzz_xy_1,     \
                             ta_yyzzz_yy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyzzzz_xx_0[i] = ta1_z_zzzz_xx_0[i] * fe_0 - ta1_z_zzzz_xx_1[i] * fe_0 + ta1_z_yzzzz_xx_0[i] * pa_y[i] - ta1_z_yzzzz_xx_1[i] * pc_y[i];

        ta1_z_yyzzzz_xy_0[i] = 3.0 * ta1_z_yyzz_xy_0[i] * fe_0 - 3.0 * ta1_z_yyzz_xy_1[i] * fe_0 + ta_yyzzz_xy_1[i] + ta1_z_yyzzz_xy_0[i] * pa_z[i] -
                               ta1_z_yyzzz_xy_1[i] * pc_z[i];

        ta1_z_yyzzzz_xz_0[i] = ta1_z_zzzz_xz_0[i] * fe_0 - ta1_z_zzzz_xz_1[i] * fe_0 + ta1_z_yzzzz_xz_0[i] * pa_y[i] - ta1_z_yzzzz_xz_1[i] * pc_y[i];

        ta1_z_yyzzzz_yy_0[i] = 3.0 * ta1_z_yyzz_yy_0[i] * fe_0 - 3.0 * ta1_z_yyzz_yy_1[i] * fe_0 + ta_yyzzz_yy_1[i] + ta1_z_yyzzz_yy_0[i] * pa_z[i] -
                               ta1_z_yyzzz_yy_1[i] * pc_z[i];

        ta1_z_yyzzzz_yz_0[i] = ta1_z_zzzz_yz_0[i] * fe_0 - ta1_z_zzzz_yz_1[i] * fe_0 + ta1_z_yzzzz_z_0[i] * fe_0 - ta1_z_yzzzz_z_1[i] * fe_0 +
                               ta1_z_yzzzz_yz_0[i] * pa_y[i] - ta1_z_yzzzz_yz_1[i] * pc_y[i];

        ta1_z_yyzzzz_zz_0[i] = ta1_z_zzzz_zz_0[i] * fe_0 - ta1_z_zzzz_zz_1[i] * fe_0 + ta1_z_yzzzz_zz_0[i] * pa_y[i] - ta1_z_yzzzz_zz_1[i] * pc_y[i];
    }

    // Set up 492-498 components of targeted buffer : ID

    auto ta1_z_yzzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 492);

    auto ta1_z_yzzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 493);

    auto ta1_z_yzzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 494);

    auto ta1_z_yzzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 495);

    auto ta1_z_yzzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 496);

    auto ta1_z_yzzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 497);

#pragma omp simd aligned(pa_y,                  \
                             pc_y,              \
                             ta1_z_yzzzzz_xx_0, \
                             ta1_z_yzzzzz_xy_0, \
                             ta1_z_yzzzzz_xz_0, \
                             ta1_z_yzzzzz_yy_0, \
                             ta1_z_yzzzzz_yz_0, \
                             ta1_z_yzzzzz_zz_0, \
                             ta1_z_zzzzz_x_0,   \
                             ta1_z_zzzzz_x_1,   \
                             ta1_z_zzzzz_xx_0,  \
                             ta1_z_zzzzz_xx_1,  \
                             ta1_z_zzzzz_xy_0,  \
                             ta1_z_zzzzz_xy_1,  \
                             ta1_z_zzzzz_xz_0,  \
                             ta1_z_zzzzz_xz_1,  \
                             ta1_z_zzzzz_y_0,   \
                             ta1_z_zzzzz_y_1,   \
                             ta1_z_zzzzz_yy_0,  \
                             ta1_z_zzzzz_yy_1,  \
                             ta1_z_zzzzz_yz_0,  \
                             ta1_z_zzzzz_yz_1,  \
                             ta1_z_zzzzz_z_0,   \
                             ta1_z_zzzzz_z_1,   \
                             ta1_z_zzzzz_zz_0,  \
                             ta1_z_zzzzz_zz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzzzzz_xx_0[i] = ta1_z_zzzzz_xx_0[i] * pa_y[i] - ta1_z_zzzzz_xx_1[i] * pc_y[i];

        ta1_z_yzzzzz_xy_0[i] = ta1_z_zzzzz_x_0[i] * fe_0 - ta1_z_zzzzz_x_1[i] * fe_0 + ta1_z_zzzzz_xy_0[i] * pa_y[i] - ta1_z_zzzzz_xy_1[i] * pc_y[i];

        ta1_z_yzzzzz_xz_0[i] = ta1_z_zzzzz_xz_0[i] * pa_y[i] - ta1_z_zzzzz_xz_1[i] * pc_y[i];

        ta1_z_yzzzzz_yy_0[i] =
            2.0 * ta1_z_zzzzz_y_0[i] * fe_0 - 2.0 * ta1_z_zzzzz_y_1[i] * fe_0 + ta1_z_zzzzz_yy_0[i] * pa_y[i] - ta1_z_zzzzz_yy_1[i] * pc_y[i];

        ta1_z_yzzzzz_yz_0[i] = ta1_z_zzzzz_z_0[i] * fe_0 - ta1_z_zzzzz_z_1[i] * fe_0 + ta1_z_zzzzz_yz_0[i] * pa_y[i] - ta1_z_zzzzz_yz_1[i] * pc_y[i];

        ta1_z_yzzzzz_zz_0[i] = ta1_z_zzzzz_zz_0[i] * pa_y[i] - ta1_z_zzzzz_zz_1[i] * pc_y[i];
    }

    // Set up 498-504 components of targeted buffer : ID

    auto ta1_z_zzzzzz_xx_0 = pbuffer.data(idx_npot_geom_010_0_id + 498);

    auto ta1_z_zzzzzz_xy_0 = pbuffer.data(idx_npot_geom_010_0_id + 499);

    auto ta1_z_zzzzzz_xz_0 = pbuffer.data(idx_npot_geom_010_0_id + 500);

    auto ta1_z_zzzzzz_yy_0 = pbuffer.data(idx_npot_geom_010_0_id + 501);

    auto ta1_z_zzzzzz_yz_0 = pbuffer.data(idx_npot_geom_010_0_id + 502);

    auto ta1_z_zzzzzz_zz_0 = pbuffer.data(idx_npot_geom_010_0_id + 503);

#pragma omp simd aligned(pa_z,                  \
                             pc_z,              \
                             ta1_z_zzzz_xx_0,   \
                             ta1_z_zzzz_xx_1,   \
                             ta1_z_zzzz_xy_0,   \
                             ta1_z_zzzz_xy_1,   \
                             ta1_z_zzzz_xz_0,   \
                             ta1_z_zzzz_xz_1,   \
                             ta1_z_zzzz_yy_0,   \
                             ta1_z_zzzz_yy_1,   \
                             ta1_z_zzzz_yz_0,   \
                             ta1_z_zzzz_yz_1,   \
                             ta1_z_zzzz_zz_0,   \
                             ta1_z_zzzz_zz_1,   \
                             ta1_z_zzzzz_x_0,   \
                             ta1_z_zzzzz_x_1,   \
                             ta1_z_zzzzz_xx_0,  \
                             ta1_z_zzzzz_xx_1,  \
                             ta1_z_zzzzz_xy_0,  \
                             ta1_z_zzzzz_xy_1,  \
                             ta1_z_zzzzz_xz_0,  \
                             ta1_z_zzzzz_xz_1,  \
                             ta1_z_zzzzz_y_0,   \
                             ta1_z_zzzzz_y_1,   \
                             ta1_z_zzzzz_yy_0,  \
                             ta1_z_zzzzz_yy_1,  \
                             ta1_z_zzzzz_yz_0,  \
                             ta1_z_zzzzz_yz_1,  \
                             ta1_z_zzzzz_z_0,   \
                             ta1_z_zzzzz_z_1,   \
                             ta1_z_zzzzz_zz_0,  \
                             ta1_z_zzzzz_zz_1,  \
                             ta1_z_zzzzzz_xx_0, \
                             ta1_z_zzzzzz_xy_0, \
                             ta1_z_zzzzzz_xz_0, \
                             ta1_z_zzzzzz_yy_0, \
                             ta1_z_zzzzzz_yz_0, \
                             ta1_z_zzzzzz_zz_0, \
                             ta_zzzzz_xx_1,     \
                             ta_zzzzz_xy_1,     \
                             ta_zzzzz_xz_1,     \
                             ta_zzzzz_yy_1,     \
                             ta_zzzzz_yz_1,     \
                             ta_zzzzz_zz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzzzzz_xx_0[i] = 5.0 * ta1_z_zzzz_xx_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xx_1[i] * fe_0 + ta_zzzzz_xx_1[i] + ta1_z_zzzzz_xx_0[i] * pa_z[i] -
                               ta1_z_zzzzz_xx_1[i] * pc_z[i];

        ta1_z_zzzzzz_xy_0[i] = 5.0 * ta1_z_zzzz_xy_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xy_1[i] * fe_0 + ta_zzzzz_xy_1[i] + ta1_z_zzzzz_xy_0[i] * pa_z[i] -
                               ta1_z_zzzzz_xy_1[i] * pc_z[i];

        ta1_z_zzzzzz_xz_0[i] = 5.0 * ta1_z_zzzz_xz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_xz_1[i] * fe_0 + ta1_z_zzzzz_x_0[i] * fe_0 -
                               ta1_z_zzzzz_x_1[i] * fe_0 + ta_zzzzz_xz_1[i] + ta1_z_zzzzz_xz_0[i] * pa_z[i] - ta1_z_zzzzz_xz_1[i] * pc_z[i];

        ta1_z_zzzzzz_yy_0[i] = 5.0 * ta1_z_zzzz_yy_0[i] * fe_0 - 5.0 * ta1_z_zzzz_yy_1[i] * fe_0 + ta_zzzzz_yy_1[i] + ta1_z_zzzzz_yy_0[i] * pa_z[i] -
                               ta1_z_zzzzz_yy_1[i] * pc_z[i];

        ta1_z_zzzzzz_yz_0[i] = 5.0 * ta1_z_zzzz_yz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_yz_1[i] * fe_0 + ta1_z_zzzzz_y_0[i] * fe_0 -
                               ta1_z_zzzzz_y_1[i] * fe_0 + ta_zzzzz_yz_1[i] + ta1_z_zzzzz_yz_0[i] * pa_z[i] - ta1_z_zzzzz_yz_1[i] * pc_z[i];

        ta1_z_zzzzzz_zz_0[i] = 5.0 * ta1_z_zzzz_zz_0[i] * fe_0 - 5.0 * ta1_z_zzzz_zz_1[i] * fe_0 + 2.0 * ta1_z_zzzzz_z_0[i] * fe_0 -
                               2.0 * ta1_z_zzzzz_z_1[i] * fe_0 + ta_zzzzz_zz_1[i] + ta1_z_zzzzz_zz_0[i] * pa_z[i] - ta1_z_zzzzz_zz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
