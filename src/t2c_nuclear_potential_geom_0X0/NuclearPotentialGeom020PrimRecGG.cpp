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

#include "NuclearPotentialGeom020PrimRecGG.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_020_gg(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_020_0_gg,
                                        const size_t              idx_npot_geom_020_0_dg,
                                        const size_t              idx_npot_geom_020_1_dg,
                                        const size_t              idx_npot_geom_020_0_ff,
                                        const size_t              idx_npot_geom_020_1_ff,
                                        const size_t              idx_npot_geom_010_1_fg,
                                        const size_t              idx_npot_geom_020_0_fg,
                                        const size_t              idx_npot_geom_020_1_fg,
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

    // Set up components of auxiliary buffer : DG

    auto ta2_xx_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg);

    auto ta2_xx_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 1);

    auto ta2_xx_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 2);

    auto ta2_xx_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 3);

    auto ta2_xx_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 4);

    auto ta2_xx_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 5);

    auto ta2_xx_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 6);

    auto ta2_xx_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 7);

    auto ta2_xx_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 8);

    auto ta2_xx_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 9);

    auto ta2_xx_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 10);

    auto ta2_xx_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 11);

    auto ta2_xx_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 12);

    auto ta2_xx_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 13);

    auto ta2_xx_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 14);

    auto ta2_xx_xy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 15);

    auto ta2_xx_xy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 17);

    auto ta2_xx_xy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 20);

    auto ta2_xx_xy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 24);

    auto ta2_xx_xz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 30);

    auto ta2_xx_xz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 31);

    auto ta2_xx_xz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 32);

    auto ta2_xx_xz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 33);

    auto ta2_xx_xz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 35);

    auto ta2_xx_xz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 36);

    auto ta2_xx_xz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 39);

    auto ta2_xx_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 45);

    auto ta2_xx_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 46);

    auto ta2_xx_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 47);

    auto ta2_xx_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 48);

    auto ta2_xx_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 49);

    auto ta2_xx_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 50);

    auto ta2_xx_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 51);

    auto ta2_xx_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 52);

    auto ta2_xx_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 53);

    auto ta2_xx_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 54);

    auto ta2_xx_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 55);

    auto ta2_xx_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 56);

    auto ta2_xx_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 57);

    auto ta2_xx_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 58);

    auto ta2_xx_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 59);

    auto ta2_xx_yz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 62);

    auto ta2_xx_yz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 65);

    auto ta2_xx_yz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 69);

    auto ta2_xx_yz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 74);

    auto ta2_xx_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 75);

    auto ta2_xx_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 76);

    auto ta2_xx_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 77);

    auto ta2_xx_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 78);

    auto ta2_xx_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 79);

    auto ta2_xx_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 80);

    auto ta2_xx_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 81);

    auto ta2_xx_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 82);

    auto ta2_xx_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 83);

    auto ta2_xx_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 84);

    auto ta2_xx_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 85);

    auto ta2_xx_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 86);

    auto ta2_xx_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 87);

    auto ta2_xx_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 88);

    auto ta2_xx_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 89);

    auto ta2_xy_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 90);

    auto ta2_xy_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 91);

    auto ta2_xy_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 92);

    auto ta2_xy_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 93);

    auto ta2_xy_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 94);

    auto ta2_xy_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 95);

    auto ta2_xy_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 96);

    auto ta2_xy_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 97);

    auto ta2_xy_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 98);

    auto ta2_xy_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 99);

    auto ta2_xy_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 100);

    auto ta2_xy_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 101);

    auto ta2_xy_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 102);

    auto ta2_xy_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 103);

    auto ta2_xy_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 104);

    auto ta2_xy_xy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 106);

    auto ta2_xy_xy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 108);

    auto ta2_xy_xy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 111);

    auto ta2_xy_xy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 115);

    auto ta2_xy_xy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 116);

    auto ta2_xy_xy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 117);

    auto ta2_xy_xy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 118);

    auto ta2_xy_xz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 120);

    auto ta2_xy_xz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 121);

    auto ta2_xy_xz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 123);

    auto ta2_xy_xz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 126);

    auto ta2_xy_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 135);

    auto ta2_xy_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 136);

    auto ta2_xy_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 137);

    auto ta2_xy_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 138);

    auto ta2_xy_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 139);

    auto ta2_xy_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 140);

    auto ta2_xy_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 141);

    auto ta2_xy_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 142);

    auto ta2_xy_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 143);

    auto ta2_xy_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 144);

    auto ta2_xy_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 145);

    auto ta2_xy_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 146);

    auto ta2_xy_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 147);

    auto ta2_xy_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 148);

    auto ta2_xy_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 149);

    auto ta2_xy_yz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 151);

    auto ta2_xy_yz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 153);

    auto ta2_xy_yz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 156);

    auto ta2_xy_yz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 160);

    auto ta2_xy_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 165);

    auto ta2_xy_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 166);

    auto ta2_xy_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 167);

    auto ta2_xy_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 168);

    auto ta2_xy_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 169);

    auto ta2_xy_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 170);

    auto ta2_xy_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 171);

    auto ta2_xy_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 172);

    auto ta2_xy_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 173);

    auto ta2_xy_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 174);

    auto ta2_xy_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 175);

    auto ta2_xy_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 176);

    auto ta2_xy_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 177);

    auto ta2_xy_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 178);

    auto ta2_xy_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 179);

    auto ta2_xz_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 180);

    auto ta2_xz_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 181);

    auto ta2_xz_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 182);

    auto ta2_xz_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 183);

    auto ta2_xz_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 184);

    auto ta2_xz_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 185);

    auto ta2_xz_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 186);

    auto ta2_xz_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 187);

    auto ta2_xz_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 188);

    auto ta2_xz_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 189);

    auto ta2_xz_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 190);

    auto ta2_xz_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 191);

    auto ta2_xz_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 192);

    auto ta2_xz_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 193);

    auto ta2_xz_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 194);

    auto ta2_xz_xy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 195);

    auto ta2_xz_xy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 197);

    auto ta2_xz_xy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 200);

    auto ta2_xz_xy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 204);

    auto ta2_xz_xz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 212);

    auto ta2_xz_xz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 215);

    auto ta2_xz_xz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 219);

    auto ta2_xz_xz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 221);

    auto ta2_xz_xz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 222);

    auto ta2_xz_xz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 223);

    auto ta2_xz_xz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 224);

    auto ta2_xz_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 225);

    auto ta2_xz_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 226);

    auto ta2_xz_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 227);

    auto ta2_xz_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 228);

    auto ta2_xz_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 229);

    auto ta2_xz_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 230);

    auto ta2_xz_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 231);

    auto ta2_xz_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 232);

    auto ta2_xz_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 233);

    auto ta2_xz_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 234);

    auto ta2_xz_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 235);

    auto ta2_xz_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 236);

    auto ta2_xz_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 237);

    auto ta2_xz_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 238);

    auto ta2_xz_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 239);

    auto ta2_xz_yz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 242);

    auto ta2_xz_yz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 245);

    auto ta2_xz_yz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 249);

    auto ta2_xz_yz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 254);

    auto ta2_xz_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 255);

    auto ta2_xz_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 256);

    auto ta2_xz_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 257);

    auto ta2_xz_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 258);

    auto ta2_xz_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 259);

    auto ta2_xz_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 260);

    auto ta2_xz_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 261);

    auto ta2_xz_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 262);

    auto ta2_xz_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 263);

    auto ta2_xz_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 264);

    auto ta2_xz_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 265);

    auto ta2_xz_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 266);

    auto ta2_xz_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 267);

    auto ta2_xz_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 268);

    auto ta2_xz_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 269);

    auto ta2_yy_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 270);

    auto ta2_yy_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 271);

    auto ta2_yy_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 272);

    auto ta2_yy_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 273);

    auto ta2_yy_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 274);

    auto ta2_yy_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 275);

    auto ta2_yy_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 276);

    auto ta2_yy_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 277);

    auto ta2_yy_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 278);

    auto ta2_yy_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 279);

    auto ta2_yy_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 280);

    auto ta2_yy_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 281);

    auto ta2_yy_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 282);

    auto ta2_yy_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 283);

    auto ta2_yy_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 284);

    auto ta2_yy_xy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 295);

    auto ta2_yy_xy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 296);

    auto ta2_yy_xy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 297);

    auto ta2_yy_xy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 298);

    auto ta2_yy_xz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 311);

    auto ta2_yy_xz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 312);

    auto ta2_yy_xz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 313);

    auto ta2_yy_xz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 314);

    auto ta2_yy_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 315);

    auto ta2_yy_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 316);

    auto ta2_yy_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 317);

    auto ta2_yy_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 318);

    auto ta2_yy_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 319);

    auto ta2_yy_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 320);

    auto ta2_yy_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 321);

    auto ta2_yy_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 322);

    auto ta2_yy_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 323);

    auto ta2_yy_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 324);

    auto ta2_yy_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 325);

    auto ta2_yy_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 326);

    auto ta2_yy_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 327);

    auto ta2_yy_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 328);

    auto ta2_yy_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 329);

    auto ta2_yy_yz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 331);

    auto ta2_yy_yz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 333);

    auto ta2_yy_yz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 336);

    auto ta2_yy_yz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 340);

    auto ta2_yy_yz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 341);

    auto ta2_yy_yz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 342);

    auto ta2_yy_yz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 343);

    auto ta2_yy_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 345);

    auto ta2_yy_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 346);

    auto ta2_yy_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 347);

    auto ta2_yy_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 348);

    auto ta2_yy_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 349);

    auto ta2_yy_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 350);

    auto ta2_yy_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 351);

    auto ta2_yy_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 352);

    auto ta2_yy_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 353);

    auto ta2_yy_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 354);

    auto ta2_yy_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 355);

    auto ta2_yy_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 356);

    auto ta2_yy_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 357);

    auto ta2_yy_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 358);

    auto ta2_yy_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 359);

    auto ta2_yz_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 360);

    auto ta2_yz_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 361);

    auto ta2_yz_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 362);

    auto ta2_yz_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 363);

    auto ta2_yz_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 364);

    auto ta2_yz_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 365);

    auto ta2_yz_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 366);

    auto ta2_yz_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 367);

    auto ta2_yz_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 368);

    auto ta2_yz_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 369);

    auto ta2_yz_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 370);

    auto ta2_yz_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 371);

    auto ta2_yz_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 372);

    auto ta2_yz_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 373);

    auto ta2_yz_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 374);

    auto ta2_yz_xy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 385);

    auto ta2_yz_xy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 386);

    auto ta2_yz_xy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 387);

    auto ta2_yz_xy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 388);

    auto ta2_yz_xz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 401);

    auto ta2_yz_xz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 402);

    auto ta2_yz_xz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 403);

    auto ta2_yz_xz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 404);

    auto ta2_yz_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 405);

    auto ta2_yz_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 406);

    auto ta2_yz_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 407);

    auto ta2_yz_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 408);

    auto ta2_yz_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 409);

    auto ta2_yz_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 410);

    auto ta2_yz_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 411);

    auto ta2_yz_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 412);

    auto ta2_yz_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 413);

    auto ta2_yz_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 414);

    auto ta2_yz_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 415);

    auto ta2_yz_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 416);

    auto ta2_yz_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 417);

    auto ta2_yz_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 418);

    auto ta2_yz_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 419);

    auto ta2_yz_yz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 422);

    auto ta2_yz_yz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 425);

    auto ta2_yz_yz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 429);

    auto ta2_yz_yz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 431);

    auto ta2_yz_yz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 432);

    auto ta2_yz_yz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 433);

    auto ta2_yz_yz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 434);

    auto ta2_yz_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 435);

    auto ta2_yz_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 436);

    auto ta2_yz_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 437);

    auto ta2_yz_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 438);

    auto ta2_yz_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 439);

    auto ta2_yz_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 440);

    auto ta2_yz_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 441);

    auto ta2_yz_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 442);

    auto ta2_yz_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 443);

    auto ta2_yz_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 444);

    auto ta2_yz_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 445);

    auto ta2_yz_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 446);

    auto ta2_yz_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 447);

    auto ta2_yz_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 448);

    auto ta2_yz_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 449);

    auto ta2_zz_xx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 450);

    auto ta2_zz_xx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 451);

    auto ta2_zz_xx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 452);

    auto ta2_zz_xx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 453);

    auto ta2_zz_xx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 454);

    auto ta2_zz_xx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 455);

    auto ta2_zz_xx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 456);

    auto ta2_zz_xx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 457);

    auto ta2_zz_xx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 458);

    auto ta2_zz_xx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 459);

    auto ta2_zz_xx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 460);

    auto ta2_zz_xx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 461);

    auto ta2_zz_xx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 462);

    auto ta2_zz_xx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 463);

    auto ta2_zz_xx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 464);

    auto ta2_zz_xy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 475);

    auto ta2_zz_xy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 476);

    auto ta2_zz_xy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 477);

    auto ta2_zz_xy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 478);

    auto ta2_zz_xz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 491);

    auto ta2_zz_xz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 492);

    auto ta2_zz_xz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 493);

    auto ta2_zz_xz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 494);

    auto ta2_zz_yy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 495);

    auto ta2_zz_yy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 496);

    auto ta2_zz_yy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 497);

    auto ta2_zz_yy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 498);

    auto ta2_zz_yy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 499);

    auto ta2_zz_yy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 500);

    auto ta2_zz_yy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 501);

    auto ta2_zz_yy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 502);

    auto ta2_zz_yy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 503);

    auto ta2_zz_yy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 504);

    auto ta2_zz_yy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 505);

    auto ta2_zz_yy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 506);

    auto ta2_zz_yy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 507);

    auto ta2_zz_yy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 508);

    auto ta2_zz_yy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 509);

    auto ta2_zz_yz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 512);

    auto ta2_zz_yz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 515);

    auto ta2_zz_yz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 519);

    auto ta2_zz_yz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 521);

    auto ta2_zz_yz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 522);

    auto ta2_zz_yz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 523);

    auto ta2_zz_yz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 524);

    auto ta2_zz_zz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_dg + 525);

    auto ta2_zz_zz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 526);

    auto ta2_zz_zz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 527);

    auto ta2_zz_zz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 528);

    auto ta2_zz_zz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 529);

    auto ta2_zz_zz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 530);

    auto ta2_zz_zz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 531);

    auto ta2_zz_zz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 532);

    auto ta2_zz_zz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 533);

    auto ta2_zz_zz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 534);

    auto ta2_zz_zz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_dg + 535);

    auto ta2_zz_zz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 536);

    auto ta2_zz_zz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 537);

    auto ta2_zz_zz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 538);

    auto ta2_zz_zz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_dg + 539);

    // Set up components of auxiliary buffer : DG

    auto ta2_xx_xx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg);

    auto ta2_xx_xx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 1);

    auto ta2_xx_xx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 2);

    auto ta2_xx_xx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 3);

    auto ta2_xx_xx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 4);

    auto ta2_xx_xx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 5);

    auto ta2_xx_xx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 6);

    auto ta2_xx_xx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 7);

    auto ta2_xx_xx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 8);

    auto ta2_xx_xx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 9);

    auto ta2_xx_xx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 10);

    auto ta2_xx_xx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 11);

    auto ta2_xx_xx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 12);

    auto ta2_xx_xx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 13);

    auto ta2_xx_xx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 14);

    auto ta2_xx_xy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 15);

    auto ta2_xx_xy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 17);

    auto ta2_xx_xy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 20);

    auto ta2_xx_xy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 24);

    auto ta2_xx_xz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 30);

    auto ta2_xx_xz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 31);

    auto ta2_xx_xz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 32);

    auto ta2_xx_xz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 33);

    auto ta2_xx_xz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 35);

    auto ta2_xx_xz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 36);

    auto ta2_xx_xz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 39);

    auto ta2_xx_yy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 45);

    auto ta2_xx_yy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 46);

    auto ta2_xx_yy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 47);

    auto ta2_xx_yy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 48);

    auto ta2_xx_yy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 49);

    auto ta2_xx_yy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 50);

    auto ta2_xx_yy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 51);

    auto ta2_xx_yy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 52);

    auto ta2_xx_yy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 53);

    auto ta2_xx_yy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 54);

    auto ta2_xx_yy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 55);

    auto ta2_xx_yy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 56);

    auto ta2_xx_yy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 57);

    auto ta2_xx_yy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 58);

    auto ta2_xx_yy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 59);

    auto ta2_xx_yz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 62);

    auto ta2_xx_yz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 65);

    auto ta2_xx_yz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 69);

    auto ta2_xx_yz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 74);

    auto ta2_xx_zz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 75);

    auto ta2_xx_zz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 76);

    auto ta2_xx_zz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 77);

    auto ta2_xx_zz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 78);

    auto ta2_xx_zz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 79);

    auto ta2_xx_zz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 80);

    auto ta2_xx_zz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 81);

    auto ta2_xx_zz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 82);

    auto ta2_xx_zz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 83);

    auto ta2_xx_zz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 84);

    auto ta2_xx_zz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 85);

    auto ta2_xx_zz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 86);

    auto ta2_xx_zz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 87);

    auto ta2_xx_zz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 88);

    auto ta2_xx_zz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 89);

    auto ta2_xy_xx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 90);

    auto ta2_xy_xx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 91);

    auto ta2_xy_xx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 92);

    auto ta2_xy_xx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 93);

    auto ta2_xy_xx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 94);

    auto ta2_xy_xx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 95);

    auto ta2_xy_xx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 96);

    auto ta2_xy_xx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 97);

    auto ta2_xy_xx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 98);

    auto ta2_xy_xx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 99);

    auto ta2_xy_xx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 100);

    auto ta2_xy_xx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 101);

    auto ta2_xy_xx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 102);

    auto ta2_xy_xx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 103);

    auto ta2_xy_xx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 104);

    auto ta2_xy_xy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 106);

    auto ta2_xy_xy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 108);

    auto ta2_xy_xy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 111);

    auto ta2_xy_xy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 115);

    auto ta2_xy_xy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 116);

    auto ta2_xy_xy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 117);

    auto ta2_xy_xy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 118);

    auto ta2_xy_xz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 120);

    auto ta2_xy_xz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 121);

    auto ta2_xy_xz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 123);

    auto ta2_xy_xz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 126);

    auto ta2_xy_yy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 135);

    auto ta2_xy_yy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 136);

    auto ta2_xy_yy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 137);

    auto ta2_xy_yy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 138);

    auto ta2_xy_yy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 139);

    auto ta2_xy_yy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 140);

    auto ta2_xy_yy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 141);

    auto ta2_xy_yy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 142);

    auto ta2_xy_yy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 143);

    auto ta2_xy_yy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 144);

    auto ta2_xy_yy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 145);

    auto ta2_xy_yy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 146);

    auto ta2_xy_yy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 147);

    auto ta2_xy_yy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 148);

    auto ta2_xy_yy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 149);

    auto ta2_xy_yz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 151);

    auto ta2_xy_yz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 153);

    auto ta2_xy_yz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 156);

    auto ta2_xy_yz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 160);

    auto ta2_xy_zz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 165);

    auto ta2_xy_zz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 166);

    auto ta2_xy_zz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 167);

    auto ta2_xy_zz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 168);

    auto ta2_xy_zz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 169);

    auto ta2_xy_zz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 170);

    auto ta2_xy_zz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 171);

    auto ta2_xy_zz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 172);

    auto ta2_xy_zz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 173);

    auto ta2_xy_zz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 174);

    auto ta2_xy_zz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 175);

    auto ta2_xy_zz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 176);

    auto ta2_xy_zz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 177);

    auto ta2_xy_zz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 178);

    auto ta2_xy_zz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 179);

    auto ta2_xz_xx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 180);

    auto ta2_xz_xx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 181);

    auto ta2_xz_xx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 182);

    auto ta2_xz_xx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 183);

    auto ta2_xz_xx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 184);

    auto ta2_xz_xx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 185);

    auto ta2_xz_xx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 186);

    auto ta2_xz_xx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 187);

    auto ta2_xz_xx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 188);

    auto ta2_xz_xx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 189);

    auto ta2_xz_xx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 190);

    auto ta2_xz_xx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 191);

    auto ta2_xz_xx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 192);

    auto ta2_xz_xx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 193);

    auto ta2_xz_xx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 194);

    auto ta2_xz_xy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 195);

    auto ta2_xz_xy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 197);

    auto ta2_xz_xy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 200);

    auto ta2_xz_xy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 204);

    auto ta2_xz_xz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 212);

    auto ta2_xz_xz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 215);

    auto ta2_xz_xz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 219);

    auto ta2_xz_xz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 221);

    auto ta2_xz_xz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 222);

    auto ta2_xz_xz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 223);

    auto ta2_xz_xz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 224);

    auto ta2_xz_yy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 225);

    auto ta2_xz_yy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 226);

    auto ta2_xz_yy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 227);

    auto ta2_xz_yy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 228);

    auto ta2_xz_yy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 229);

    auto ta2_xz_yy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 230);

    auto ta2_xz_yy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 231);

    auto ta2_xz_yy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 232);

    auto ta2_xz_yy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 233);

    auto ta2_xz_yy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 234);

    auto ta2_xz_yy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 235);

    auto ta2_xz_yy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 236);

    auto ta2_xz_yy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 237);

    auto ta2_xz_yy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 238);

    auto ta2_xz_yy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 239);

    auto ta2_xz_yz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 242);

    auto ta2_xz_yz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 245);

    auto ta2_xz_yz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 249);

    auto ta2_xz_yz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 254);

    auto ta2_xz_zz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 255);

    auto ta2_xz_zz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 256);

    auto ta2_xz_zz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 257);

    auto ta2_xz_zz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 258);

    auto ta2_xz_zz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 259);

    auto ta2_xz_zz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 260);

    auto ta2_xz_zz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 261);

    auto ta2_xz_zz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 262);

    auto ta2_xz_zz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 263);

    auto ta2_xz_zz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 264);

    auto ta2_xz_zz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 265);

    auto ta2_xz_zz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 266);

    auto ta2_xz_zz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 267);

    auto ta2_xz_zz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 268);

    auto ta2_xz_zz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 269);

    auto ta2_yy_xx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 270);

    auto ta2_yy_xx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 271);

    auto ta2_yy_xx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 272);

    auto ta2_yy_xx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 273);

    auto ta2_yy_xx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 274);

    auto ta2_yy_xx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 275);

    auto ta2_yy_xx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 276);

    auto ta2_yy_xx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 277);

    auto ta2_yy_xx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 278);

    auto ta2_yy_xx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 279);

    auto ta2_yy_xx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 280);

    auto ta2_yy_xx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 281);

    auto ta2_yy_xx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 282);

    auto ta2_yy_xx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 283);

    auto ta2_yy_xx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 284);

    auto ta2_yy_xy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 295);

    auto ta2_yy_xy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 296);

    auto ta2_yy_xy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 297);

    auto ta2_yy_xy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 298);

    auto ta2_yy_xz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 311);

    auto ta2_yy_xz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 312);

    auto ta2_yy_xz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 313);

    auto ta2_yy_xz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 314);

    auto ta2_yy_yy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 315);

    auto ta2_yy_yy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 316);

    auto ta2_yy_yy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 317);

    auto ta2_yy_yy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 318);

    auto ta2_yy_yy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 319);

    auto ta2_yy_yy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 320);

    auto ta2_yy_yy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 321);

    auto ta2_yy_yy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 322);

    auto ta2_yy_yy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 323);

    auto ta2_yy_yy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 324);

    auto ta2_yy_yy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 325);

    auto ta2_yy_yy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 326);

    auto ta2_yy_yy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 327);

    auto ta2_yy_yy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 328);

    auto ta2_yy_yy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 329);

    auto ta2_yy_yz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 331);

    auto ta2_yy_yz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 333);

    auto ta2_yy_yz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 336);

    auto ta2_yy_yz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 340);

    auto ta2_yy_yz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 341);

    auto ta2_yy_yz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 342);

    auto ta2_yy_yz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 343);

    auto ta2_yy_zz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 345);

    auto ta2_yy_zz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 346);

    auto ta2_yy_zz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 347);

    auto ta2_yy_zz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 348);

    auto ta2_yy_zz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 349);

    auto ta2_yy_zz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 350);

    auto ta2_yy_zz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 351);

    auto ta2_yy_zz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 352);

    auto ta2_yy_zz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 353);

    auto ta2_yy_zz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 354);

    auto ta2_yy_zz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 355);

    auto ta2_yy_zz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 356);

    auto ta2_yy_zz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 357);

    auto ta2_yy_zz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 358);

    auto ta2_yy_zz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 359);

    auto ta2_yz_xx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 360);

    auto ta2_yz_xx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 361);

    auto ta2_yz_xx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 362);

    auto ta2_yz_xx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 363);

    auto ta2_yz_xx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 364);

    auto ta2_yz_xx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 365);

    auto ta2_yz_xx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 366);

    auto ta2_yz_xx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 367);

    auto ta2_yz_xx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 368);

    auto ta2_yz_xx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 369);

    auto ta2_yz_xx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 370);

    auto ta2_yz_xx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 371);

    auto ta2_yz_xx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 372);

    auto ta2_yz_xx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 373);

    auto ta2_yz_xx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 374);

    auto ta2_yz_xy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 385);

    auto ta2_yz_xy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 386);

    auto ta2_yz_xy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 387);

    auto ta2_yz_xy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 388);

    auto ta2_yz_xz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 401);

    auto ta2_yz_xz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 402);

    auto ta2_yz_xz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 403);

    auto ta2_yz_xz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 404);

    auto ta2_yz_yy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 405);

    auto ta2_yz_yy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 406);

    auto ta2_yz_yy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 407);

    auto ta2_yz_yy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 408);

    auto ta2_yz_yy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 409);

    auto ta2_yz_yy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 410);

    auto ta2_yz_yy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 411);

    auto ta2_yz_yy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 412);

    auto ta2_yz_yy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 413);

    auto ta2_yz_yy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 414);

    auto ta2_yz_yy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 415);

    auto ta2_yz_yy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 416);

    auto ta2_yz_yy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 417);

    auto ta2_yz_yy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 418);

    auto ta2_yz_yy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 419);

    auto ta2_yz_yz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 422);

    auto ta2_yz_yz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 425);

    auto ta2_yz_yz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 429);

    auto ta2_yz_yz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 431);

    auto ta2_yz_yz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 432);

    auto ta2_yz_yz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 433);

    auto ta2_yz_yz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 434);

    auto ta2_yz_zz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 435);

    auto ta2_yz_zz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 436);

    auto ta2_yz_zz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 437);

    auto ta2_yz_zz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 438);

    auto ta2_yz_zz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 439);

    auto ta2_yz_zz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 440);

    auto ta2_yz_zz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 441);

    auto ta2_yz_zz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 442);

    auto ta2_yz_zz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 443);

    auto ta2_yz_zz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 444);

    auto ta2_yz_zz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 445);

    auto ta2_yz_zz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 446);

    auto ta2_yz_zz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 447);

    auto ta2_yz_zz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 448);

    auto ta2_yz_zz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 449);

    auto ta2_zz_xx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 450);

    auto ta2_zz_xx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 451);

    auto ta2_zz_xx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 452);

    auto ta2_zz_xx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 453);

    auto ta2_zz_xx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 454);

    auto ta2_zz_xx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 455);

    auto ta2_zz_xx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 456);

    auto ta2_zz_xx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 457);

    auto ta2_zz_xx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 458);

    auto ta2_zz_xx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 459);

    auto ta2_zz_xx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 460);

    auto ta2_zz_xx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 461);

    auto ta2_zz_xx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 462);

    auto ta2_zz_xx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 463);

    auto ta2_zz_xx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 464);

    auto ta2_zz_xy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 475);

    auto ta2_zz_xy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 476);

    auto ta2_zz_xy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 477);

    auto ta2_zz_xy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 478);

    auto ta2_zz_xz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 491);

    auto ta2_zz_xz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 492);

    auto ta2_zz_xz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 493);

    auto ta2_zz_xz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 494);

    auto ta2_zz_yy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 495);

    auto ta2_zz_yy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 496);

    auto ta2_zz_yy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 497);

    auto ta2_zz_yy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 498);

    auto ta2_zz_yy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 499);

    auto ta2_zz_yy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 500);

    auto ta2_zz_yy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 501);

    auto ta2_zz_yy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 502);

    auto ta2_zz_yy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 503);

    auto ta2_zz_yy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 504);

    auto ta2_zz_yy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 505);

    auto ta2_zz_yy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 506);

    auto ta2_zz_yy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 507);

    auto ta2_zz_yy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 508);

    auto ta2_zz_yy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 509);

    auto ta2_zz_yz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 512);

    auto ta2_zz_yz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 515);

    auto ta2_zz_yz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 519);

    auto ta2_zz_yz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 521);

    auto ta2_zz_yz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 522);

    auto ta2_zz_yz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 523);

    auto ta2_zz_yz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 524);

    auto ta2_zz_zz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_dg + 525);

    auto ta2_zz_zz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 526);

    auto ta2_zz_zz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 527);

    auto ta2_zz_zz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 528);

    auto ta2_zz_zz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 529);

    auto ta2_zz_zz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 530);

    auto ta2_zz_zz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 531);

    auto ta2_zz_zz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 532);

    auto ta2_zz_zz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 533);

    auto ta2_zz_zz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 534);

    auto ta2_zz_zz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_dg + 535);

    auto ta2_zz_zz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 536);

    auto ta2_zz_zz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 537);

    auto ta2_zz_zz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 538);

    auto ta2_zz_zz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_dg + 539);

    // Set up components of auxiliary buffer : FF

    auto ta2_xx_xxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff);

    auto ta2_xx_xxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 1);

    auto ta2_xx_xxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 2);

    auto ta2_xx_xxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 3);

    auto ta2_xx_xxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 4);

    auto ta2_xx_xxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 5);

    auto ta2_xx_xxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 6);

    auto ta2_xx_xxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 7);

    auto ta2_xx_xxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 8);

    auto ta2_xx_xxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 9);

    auto ta2_xx_xxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 10);

    auto ta2_xx_xxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 11);

    auto ta2_xx_xxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 12);

    auto ta2_xx_xxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 13);

    auto ta2_xx_xxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 14);

    auto ta2_xx_xxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 15);

    auto ta2_xx_xxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 20);

    auto ta2_xx_xxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 21);

    auto ta2_xx_xxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 22);

    auto ta2_xx_xxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 23);

    auto ta2_xx_xxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 24);

    auto ta2_xx_xxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 25);

    auto ta2_xx_xxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 27);

    auto ta2_xx_xxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 28);

    auto ta2_xx_xxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 29);

    auto ta2_xx_xyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 31);

    auto ta2_xx_xyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 33);

    auto ta2_xx_xyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 34);

    auto ta2_xx_xzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 50);

    auto ta2_xx_xzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 51);

    auto ta2_xx_xzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 52);

    auto ta2_xx_xzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 53);

    auto ta2_xx_xzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 54);

    auto ta2_xx_xzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 55);

    auto ta2_xx_yyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 60);

    auto ta2_xx_yyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 61);

    auto ta2_xx_yyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 62);

    auto ta2_xx_yyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 63);

    auto ta2_xx_yyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 64);

    auto ta2_xx_yyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 65);

    auto ta2_xx_yyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 66);

    auto ta2_xx_yyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 67);

    auto ta2_xx_yyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 68);

    auto ta2_xx_yyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 69);

    auto ta2_xx_yzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 82);

    auto ta2_xx_yzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 84);

    auto ta2_xx_yzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 85);

    auto ta2_xx_yzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 87);

    auto ta2_xx_yzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 88);

    auto ta2_xx_yzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 89);

    auto ta2_xx_zzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 90);

    auto ta2_xx_zzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 91);

    auto ta2_xx_zzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 92);

    auto ta2_xx_zzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 93);

    auto ta2_xx_zzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 94);

    auto ta2_xx_zzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 95);

    auto ta2_xx_zzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 96);

    auto ta2_xx_zzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 97);

    auto ta2_xx_zzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 98);

    auto ta2_xx_zzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 99);

    auto ta2_xy_xxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 100);

    auto ta2_xy_xxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 101);

    auto ta2_xy_xxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 102);

    auto ta2_xy_xxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 103);

    auto ta2_xy_xxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 104);

    auto ta2_xy_xxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 105);

    auto ta2_xy_xxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 106);

    auto ta2_xy_xxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 107);

    auto ta2_xy_xxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 108);

    auto ta2_xy_xxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 109);

    auto ta2_xy_xxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 111);

    auto ta2_xy_xxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 113);

    auto ta2_xy_xxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 114);

    auto ta2_xy_xxy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 116);

    auto ta2_xy_xxy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 117);

    auto ta2_xy_xxy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 118);

    auto ta2_xy_xxz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 120);

    auto ta2_xy_xxz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 121);

    auto ta2_xy_xxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 122);

    auto ta2_xy_xxz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 123);

    auto ta2_xy_xxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 124);

    auto ta2_xy_xxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 125);

    auto ta2_xy_xyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 130);

    auto ta2_xy_xyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 131);

    auto ta2_xy_xyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 132);

    auto ta2_xy_xyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 133);

    auto ta2_xy_xyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 134);

    auto ta2_xy_xyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 135);

    auto ta2_xy_xyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 136);

    auto ta2_xy_xyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 137);

    auto ta2_xy_xyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 138);

    auto ta2_xy_yyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 160);

    auto ta2_xy_yyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 161);

    auto ta2_xy_yyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 162);

    auto ta2_xy_yyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 163);

    auto ta2_xy_yyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 164);

    auto ta2_xy_yyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 165);

    auto ta2_xy_yyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 166);

    auto ta2_xy_yyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 167);

    auto ta2_xy_yyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 168);

    auto ta2_xy_yyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 169);

    auto ta2_xy_yyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 171);

    auto ta2_xy_yyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 173);

    auto ta2_xy_yyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 174);

    auto ta2_xy_yyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 176);

    auto ta2_xy_yyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 177);

    auto ta2_xy_yyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 178);

    auto ta2_xy_yzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 184);

    auto ta2_xy_yzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 187);

    auto ta2_xy_yzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 188);

    auto ta2_xy_zzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 190);

    auto ta2_xy_zzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 191);

    auto ta2_xy_zzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 192);

    auto ta2_xy_zzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 193);

    auto ta2_xy_zzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 194);

    auto ta2_xy_zzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 195);

    auto ta2_xy_zzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 196);

    auto ta2_xy_zzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 197);

    auto ta2_xy_zzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 198);

    auto ta2_xy_zzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 199);

    auto ta2_xz_xxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 200);

    auto ta2_xz_xxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 201);

    auto ta2_xz_xxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 202);

    auto ta2_xz_xxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 203);

    auto ta2_xz_xxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 204);

    auto ta2_xz_xxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 205);

    auto ta2_xz_xxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 206);

    auto ta2_xz_xxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 207);

    auto ta2_xz_xxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 208);

    auto ta2_xz_xxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 209);

    auto ta2_xz_xxy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 210);

    auto ta2_xz_xxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 211);

    auto ta2_xz_xxy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 212);

    auto ta2_xz_xxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 213);

    auto ta2_xz_xxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 214);

    auto ta2_xz_xxy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 215);

    auto ta2_xz_xxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 222);

    auto ta2_xz_xxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 224);

    auto ta2_xz_xxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 225);

    auto ta2_xz_xxz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 227);

    auto ta2_xz_xxz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 228);

    auto ta2_xz_xxz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 229);

    auto ta2_xz_xzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 250);

    auto ta2_xz_xzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 251);

    auto ta2_xz_xzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 252);

    auto ta2_xz_xzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 253);

    auto ta2_xz_xzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 254);

    auto ta2_xz_xzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 255);

    auto ta2_xz_xzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 257);

    auto ta2_xz_xzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 258);

    auto ta2_xz_xzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 259);

    auto ta2_xz_yyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 260);

    auto ta2_xz_yyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 261);

    auto ta2_xz_yyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 262);

    auto ta2_xz_yyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 263);

    auto ta2_xz_yyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 264);

    auto ta2_xz_yyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 265);

    auto ta2_xz_yyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 266);

    auto ta2_xz_yyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 267);

    auto ta2_xz_yyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 268);

    auto ta2_xz_yyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 269);

    auto ta2_xz_yyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 274);

    auto ta2_xz_yyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 277);

    auto ta2_xz_yyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 278);

    auto ta2_xz_yzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 282);

    auto ta2_xz_yzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 284);

    auto ta2_xz_yzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 285);

    auto ta2_xz_yzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 287);

    auto ta2_xz_yzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 288);

    auto ta2_xz_yzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 289);

    auto ta2_xz_zzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 290);

    auto ta2_xz_zzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 291);

    auto ta2_xz_zzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 292);

    auto ta2_xz_zzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 293);

    auto ta2_xz_zzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 294);

    auto ta2_xz_zzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 295);

    auto ta2_xz_zzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 296);

    auto ta2_xz_zzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 297);

    auto ta2_xz_zzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 298);

    auto ta2_xz_zzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 299);

    auto ta2_yy_xxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 300);

    auto ta2_yy_xxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 301);

    auto ta2_yy_xxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 302);

    auto ta2_yy_xxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 303);

    auto ta2_yy_xxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 304);

    auto ta2_yy_xxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 305);

    auto ta2_yy_xxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 306);

    auto ta2_yy_xxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 307);

    auto ta2_yy_xxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 308);

    auto ta2_yy_xxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 309);

    auto ta2_yy_xxy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 311);

    auto ta2_yy_xxy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 313);

    auto ta2_yy_xxy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 314);

    auto ta2_yy_xyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 331);

    auto ta2_yy_xyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 333);

    auto ta2_yy_xyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 334);

    auto ta2_yy_xyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 336);

    auto ta2_yy_xyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 337);

    auto ta2_yy_xyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 338);

    auto ta2_yy_xzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 352);

    auto ta2_yy_xzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 354);

    auto ta2_yy_xzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 355);

    auto ta2_yy_xzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 357);

    auto ta2_yy_xzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 358);

    auto ta2_yy_xzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 359);

    auto ta2_yy_yyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 360);

    auto ta2_yy_yyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 361);

    auto ta2_yy_yyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 362);

    auto ta2_yy_yyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 363);

    auto ta2_yy_yyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 364);

    auto ta2_yy_yyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 365);

    auto ta2_yy_yyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 366);

    auto ta2_yy_yyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 367);

    auto ta2_yy_yyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 368);

    auto ta2_yy_yyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 369);

    auto ta2_yy_yyz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 371);

    auto ta2_yy_yyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 372);

    auto ta2_yy_yyz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 373);

    auto ta2_yy_yyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 374);

    auto ta2_yy_yyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 375);

    auto ta2_yy_yyz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 376);

    auto ta2_yy_yyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 377);

    auto ta2_yy_yyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 378);

    auto ta2_yy_yyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 379);

    auto ta2_yy_yzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 381);

    auto ta2_yy_yzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 383);

    auto ta2_yy_yzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 384);

    auto ta2_yy_yzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 386);

    auto ta2_yy_yzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 387);

    auto ta2_yy_yzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 388);

    auto ta2_yy_zzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 390);

    auto ta2_yy_zzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 391);

    auto ta2_yy_zzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 392);

    auto ta2_yy_zzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 393);

    auto ta2_yy_zzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 394);

    auto ta2_yy_zzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 395);

    auto ta2_yy_zzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 396);

    auto ta2_yy_zzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 397);

    auto ta2_yy_zzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 398);

    auto ta2_yy_zzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 399);

    auto ta2_yz_xxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 400);

    auto ta2_yz_xxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 401);

    auto ta2_yz_xxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 402);

    auto ta2_yz_xxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 403);

    auto ta2_yz_xxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 404);

    auto ta2_yz_xxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 405);

    auto ta2_yz_xxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 406);

    auto ta2_yz_xxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 407);

    auto ta2_yz_xxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 408);

    auto ta2_yz_xxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 409);

    auto ta2_yz_xxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 422);

    auto ta2_yz_xxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 424);

    auto ta2_yz_xxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 425);

    auto ta2_yz_xyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 431);

    auto ta2_yz_xyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 433);

    auto ta2_yz_xyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 434);

    auto ta2_yz_xyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 436);

    auto ta2_yz_xyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 437);

    auto ta2_yz_xyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 438);

    auto ta2_yz_xzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 452);

    auto ta2_yz_xzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 454);

    auto ta2_yz_xzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 455);

    auto ta2_yz_xzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 457);

    auto ta2_yz_xzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 458);

    auto ta2_yz_xzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 459);

    auto ta2_yz_yyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 460);

    auto ta2_yz_yyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 461);

    auto ta2_yz_yyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 462);

    auto ta2_yz_yyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 463);

    auto ta2_yz_yyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 464);

    auto ta2_yz_yyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 465);

    auto ta2_yz_yyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 466);

    auto ta2_yz_yyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 467);

    auto ta2_yz_yyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 468);

    auto ta2_yz_yyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 469);

    auto ta2_yz_yyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 472);

    auto ta2_yz_yyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 474);

    auto ta2_yz_yyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 475);

    auto ta2_yz_yyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 477);

    auto ta2_yz_yyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 478);

    auto ta2_yz_yyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 479);

    auto ta2_yz_yzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 481);

    auto ta2_yz_yzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 482);

    auto ta2_yz_yzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 483);

    auto ta2_yz_yzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 484);

    auto ta2_yz_yzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 485);

    auto ta2_yz_yzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 486);

    auto ta2_yz_yzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 487);

    auto ta2_yz_yzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 488);

    auto ta2_yz_yzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 489);

    auto ta2_yz_zzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 490);

    auto ta2_yz_zzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 491);

    auto ta2_yz_zzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 492);

    auto ta2_yz_zzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 493);

    auto ta2_yz_zzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 494);

    auto ta2_yz_zzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 495);

    auto ta2_yz_zzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 496);

    auto ta2_yz_zzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 497);

    auto ta2_yz_zzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 498);

    auto ta2_yz_zzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 499);

    auto ta2_zz_xxx_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 500);

    auto ta2_zz_xxx_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 501);

    auto ta2_zz_xxx_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 502);

    auto ta2_zz_xxx_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 503);

    auto ta2_zz_xxx_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 504);

    auto ta2_zz_xxx_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 505);

    auto ta2_zz_xxx_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 506);

    auto ta2_zz_xxx_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 507);

    auto ta2_zz_xxx_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 508);

    auto ta2_zz_xxx_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 509);

    auto ta2_zz_xxz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 522);

    auto ta2_zz_xxz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 524);

    auto ta2_zz_xxz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 525);

    auto ta2_zz_xyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 531);

    auto ta2_zz_xyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 533);

    auto ta2_zz_xyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 534);

    auto ta2_zz_xyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 536);

    auto ta2_zz_xyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 537);

    auto ta2_zz_xyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 538);

    auto ta2_zz_xzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 552);

    auto ta2_zz_xzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 554);

    auto ta2_zz_xzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 555);

    auto ta2_zz_xzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 557);

    auto ta2_zz_xzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 558);

    auto ta2_zz_xzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 559);

    auto ta2_zz_yyy_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 560);

    auto ta2_zz_yyy_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 561);

    auto ta2_zz_yyy_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 562);

    auto ta2_zz_yyy_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 563);

    auto ta2_zz_yyy_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 564);

    auto ta2_zz_yyy_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 565);

    auto ta2_zz_yyy_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 566);

    auto ta2_zz_yyy_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 567);

    auto ta2_zz_yyy_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 568);

    auto ta2_zz_yyy_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 569);

    auto ta2_zz_yyz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 572);

    auto ta2_zz_yyz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 574);

    auto ta2_zz_yyz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 575);

    auto ta2_zz_yyz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 577);

    auto ta2_zz_yyz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 578);

    auto ta2_zz_yyz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 579);

    auto ta2_zz_yzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 581);

    auto ta2_zz_yzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 582);

    auto ta2_zz_yzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 583);

    auto ta2_zz_yzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 584);

    auto ta2_zz_yzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 585);

    auto ta2_zz_yzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 586);

    auto ta2_zz_yzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 587);

    auto ta2_zz_yzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 588);

    auto ta2_zz_yzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 589);

    auto ta2_zz_zzz_xxx_0 = pbuffer.data(idx_npot_geom_020_0_ff + 590);

    auto ta2_zz_zzz_xxy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 591);

    auto ta2_zz_zzz_xxz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 592);

    auto ta2_zz_zzz_xyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 593);

    auto ta2_zz_zzz_xyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 594);

    auto ta2_zz_zzz_xzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 595);

    auto ta2_zz_zzz_yyy_0 = pbuffer.data(idx_npot_geom_020_0_ff + 596);

    auto ta2_zz_zzz_yyz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 597);

    auto ta2_zz_zzz_yzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 598);

    auto ta2_zz_zzz_zzz_0 = pbuffer.data(idx_npot_geom_020_0_ff + 599);

    // Set up components of auxiliary buffer : FF

    auto ta2_xx_xxx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff);

    auto ta2_xx_xxx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 1);

    auto ta2_xx_xxx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 2);

    auto ta2_xx_xxx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 3);

    auto ta2_xx_xxx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 4);

    auto ta2_xx_xxx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 5);

    auto ta2_xx_xxx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 6);

    auto ta2_xx_xxx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 7);

    auto ta2_xx_xxx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 8);

    auto ta2_xx_xxx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 9);

    auto ta2_xx_xxy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 10);

    auto ta2_xx_xxy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 11);

    auto ta2_xx_xxy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 12);

    auto ta2_xx_xxy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 13);

    auto ta2_xx_xxy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 14);

    auto ta2_xx_xxy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 15);

    auto ta2_xx_xxz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 20);

    auto ta2_xx_xxz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 21);

    auto ta2_xx_xxz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 22);

    auto ta2_xx_xxz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 23);

    auto ta2_xx_xxz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 24);

    auto ta2_xx_xxz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 25);

    auto ta2_xx_xxz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 27);

    auto ta2_xx_xxz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 28);

    auto ta2_xx_xxz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 29);

    auto ta2_xx_xyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 31);

    auto ta2_xx_xyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 33);

    auto ta2_xx_xyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 34);

    auto ta2_xx_xzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 50);

    auto ta2_xx_xzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 51);

    auto ta2_xx_xzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 52);

    auto ta2_xx_xzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 53);

    auto ta2_xx_xzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 54);

    auto ta2_xx_xzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 55);

    auto ta2_xx_yyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 60);

    auto ta2_xx_yyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 61);

    auto ta2_xx_yyy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 62);

    auto ta2_xx_yyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 63);

    auto ta2_xx_yyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 64);

    auto ta2_xx_yyy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 65);

    auto ta2_xx_yyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 66);

    auto ta2_xx_yyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 67);

    auto ta2_xx_yyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 68);

    auto ta2_xx_yyy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 69);

    auto ta2_xx_yzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 82);

    auto ta2_xx_yzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 84);

    auto ta2_xx_yzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 85);

    auto ta2_xx_yzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 87);

    auto ta2_xx_yzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 88);

    auto ta2_xx_yzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 89);

    auto ta2_xx_zzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 90);

    auto ta2_xx_zzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 91);

    auto ta2_xx_zzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 92);

    auto ta2_xx_zzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 93);

    auto ta2_xx_zzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 94);

    auto ta2_xx_zzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 95);

    auto ta2_xx_zzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 96);

    auto ta2_xx_zzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 97);

    auto ta2_xx_zzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 98);

    auto ta2_xx_zzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 99);

    auto ta2_xy_xxx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 100);

    auto ta2_xy_xxx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 101);

    auto ta2_xy_xxx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 102);

    auto ta2_xy_xxx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 103);

    auto ta2_xy_xxx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 104);

    auto ta2_xy_xxx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 105);

    auto ta2_xy_xxx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 106);

    auto ta2_xy_xxx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 107);

    auto ta2_xy_xxx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 108);

    auto ta2_xy_xxx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 109);

    auto ta2_xy_xxy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 111);

    auto ta2_xy_xxy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 113);

    auto ta2_xy_xxy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 114);

    auto ta2_xy_xxy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 116);

    auto ta2_xy_xxy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 117);

    auto ta2_xy_xxy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 118);

    auto ta2_xy_xxz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 120);

    auto ta2_xy_xxz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 121);

    auto ta2_xy_xxz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 122);

    auto ta2_xy_xxz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 123);

    auto ta2_xy_xxz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 124);

    auto ta2_xy_xxz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 125);

    auto ta2_xy_xyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 130);

    auto ta2_xy_xyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 131);

    auto ta2_xy_xyy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 132);

    auto ta2_xy_xyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 133);

    auto ta2_xy_xyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 134);

    auto ta2_xy_xyy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 135);

    auto ta2_xy_xyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 136);

    auto ta2_xy_xyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 137);

    auto ta2_xy_xyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 138);

    auto ta2_xy_yyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 160);

    auto ta2_xy_yyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 161);

    auto ta2_xy_yyy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 162);

    auto ta2_xy_yyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 163);

    auto ta2_xy_yyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 164);

    auto ta2_xy_yyy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 165);

    auto ta2_xy_yyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 166);

    auto ta2_xy_yyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 167);

    auto ta2_xy_yyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 168);

    auto ta2_xy_yyy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 169);

    auto ta2_xy_yyz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 171);

    auto ta2_xy_yyz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 173);

    auto ta2_xy_yyz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 174);

    auto ta2_xy_yyz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 176);

    auto ta2_xy_yyz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 177);

    auto ta2_xy_yyz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 178);

    auto ta2_xy_yzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 184);

    auto ta2_xy_yzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 187);

    auto ta2_xy_yzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 188);

    auto ta2_xy_zzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 190);

    auto ta2_xy_zzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 191);

    auto ta2_xy_zzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 192);

    auto ta2_xy_zzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 193);

    auto ta2_xy_zzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 194);

    auto ta2_xy_zzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 195);

    auto ta2_xy_zzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 196);

    auto ta2_xy_zzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 197);

    auto ta2_xy_zzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 198);

    auto ta2_xy_zzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 199);

    auto ta2_xz_xxx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 200);

    auto ta2_xz_xxx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 201);

    auto ta2_xz_xxx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 202);

    auto ta2_xz_xxx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 203);

    auto ta2_xz_xxx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 204);

    auto ta2_xz_xxx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 205);

    auto ta2_xz_xxx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 206);

    auto ta2_xz_xxx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 207);

    auto ta2_xz_xxx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 208);

    auto ta2_xz_xxx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 209);

    auto ta2_xz_xxy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 210);

    auto ta2_xz_xxy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 211);

    auto ta2_xz_xxy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 212);

    auto ta2_xz_xxy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 213);

    auto ta2_xz_xxy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 214);

    auto ta2_xz_xxy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 215);

    auto ta2_xz_xxz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 222);

    auto ta2_xz_xxz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 224);

    auto ta2_xz_xxz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 225);

    auto ta2_xz_xxz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 227);

    auto ta2_xz_xxz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 228);

    auto ta2_xz_xxz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 229);

    auto ta2_xz_xzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 250);

    auto ta2_xz_xzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 251);

    auto ta2_xz_xzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 252);

    auto ta2_xz_xzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 253);

    auto ta2_xz_xzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 254);

    auto ta2_xz_xzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 255);

    auto ta2_xz_xzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 257);

    auto ta2_xz_xzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 258);

    auto ta2_xz_xzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 259);

    auto ta2_xz_yyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 260);

    auto ta2_xz_yyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 261);

    auto ta2_xz_yyy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 262);

    auto ta2_xz_yyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 263);

    auto ta2_xz_yyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 264);

    auto ta2_xz_yyy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 265);

    auto ta2_xz_yyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 266);

    auto ta2_xz_yyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 267);

    auto ta2_xz_yyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 268);

    auto ta2_xz_yyy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 269);

    auto ta2_xz_yyz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 274);

    auto ta2_xz_yyz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 277);

    auto ta2_xz_yyz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 278);

    auto ta2_xz_yzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 282);

    auto ta2_xz_yzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 284);

    auto ta2_xz_yzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 285);

    auto ta2_xz_yzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 287);

    auto ta2_xz_yzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 288);

    auto ta2_xz_yzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 289);

    auto ta2_xz_zzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 290);

    auto ta2_xz_zzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 291);

    auto ta2_xz_zzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 292);

    auto ta2_xz_zzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 293);

    auto ta2_xz_zzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 294);

    auto ta2_xz_zzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 295);

    auto ta2_xz_zzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 296);

    auto ta2_xz_zzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 297);

    auto ta2_xz_zzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 298);

    auto ta2_xz_zzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 299);

    auto ta2_yy_xxx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 300);

    auto ta2_yy_xxx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 301);

    auto ta2_yy_xxx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 302);

    auto ta2_yy_xxx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 303);

    auto ta2_yy_xxx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 304);

    auto ta2_yy_xxx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 305);

    auto ta2_yy_xxx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 306);

    auto ta2_yy_xxx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 307);

    auto ta2_yy_xxx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 308);

    auto ta2_yy_xxx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 309);

    auto ta2_yy_xxy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 311);

    auto ta2_yy_xxy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 313);

    auto ta2_yy_xxy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 314);

    auto ta2_yy_xyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 331);

    auto ta2_yy_xyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 333);

    auto ta2_yy_xyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 334);

    auto ta2_yy_xyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 336);

    auto ta2_yy_xyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 337);

    auto ta2_yy_xyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 338);

    auto ta2_yy_xzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 352);

    auto ta2_yy_xzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 354);

    auto ta2_yy_xzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 355);

    auto ta2_yy_xzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 357);

    auto ta2_yy_xzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 358);

    auto ta2_yy_xzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 359);

    auto ta2_yy_yyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 360);

    auto ta2_yy_yyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 361);

    auto ta2_yy_yyy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 362);

    auto ta2_yy_yyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 363);

    auto ta2_yy_yyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 364);

    auto ta2_yy_yyy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 365);

    auto ta2_yy_yyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 366);

    auto ta2_yy_yyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 367);

    auto ta2_yy_yyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 368);

    auto ta2_yy_yyy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 369);

    auto ta2_yy_yyz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 371);

    auto ta2_yy_yyz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 372);

    auto ta2_yy_yyz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 373);

    auto ta2_yy_yyz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 374);

    auto ta2_yy_yyz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 375);

    auto ta2_yy_yyz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 376);

    auto ta2_yy_yyz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 377);

    auto ta2_yy_yyz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 378);

    auto ta2_yy_yyz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 379);

    auto ta2_yy_yzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 381);

    auto ta2_yy_yzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 383);

    auto ta2_yy_yzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 384);

    auto ta2_yy_yzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 386);

    auto ta2_yy_yzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 387);

    auto ta2_yy_yzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 388);

    auto ta2_yy_zzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 390);

    auto ta2_yy_zzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 391);

    auto ta2_yy_zzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 392);

    auto ta2_yy_zzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 393);

    auto ta2_yy_zzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 394);

    auto ta2_yy_zzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 395);

    auto ta2_yy_zzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 396);

    auto ta2_yy_zzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 397);

    auto ta2_yy_zzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 398);

    auto ta2_yy_zzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 399);

    auto ta2_yz_xxx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 400);

    auto ta2_yz_xxx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 401);

    auto ta2_yz_xxx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 402);

    auto ta2_yz_xxx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 403);

    auto ta2_yz_xxx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 404);

    auto ta2_yz_xxx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 405);

    auto ta2_yz_xxx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 406);

    auto ta2_yz_xxx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 407);

    auto ta2_yz_xxx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 408);

    auto ta2_yz_xxx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 409);

    auto ta2_yz_xxz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 422);

    auto ta2_yz_xxz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 424);

    auto ta2_yz_xxz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 425);

    auto ta2_yz_xyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 431);

    auto ta2_yz_xyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 433);

    auto ta2_yz_xyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 434);

    auto ta2_yz_xyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 436);

    auto ta2_yz_xyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 437);

    auto ta2_yz_xyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 438);

    auto ta2_yz_xzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 452);

    auto ta2_yz_xzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 454);

    auto ta2_yz_xzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 455);

    auto ta2_yz_xzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 457);

    auto ta2_yz_xzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 458);

    auto ta2_yz_xzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 459);

    auto ta2_yz_yyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 460);

    auto ta2_yz_yyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 461);

    auto ta2_yz_yyy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 462);

    auto ta2_yz_yyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 463);

    auto ta2_yz_yyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 464);

    auto ta2_yz_yyy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 465);

    auto ta2_yz_yyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 466);

    auto ta2_yz_yyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 467);

    auto ta2_yz_yyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 468);

    auto ta2_yz_yyy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 469);

    auto ta2_yz_yyz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 472);

    auto ta2_yz_yyz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 474);

    auto ta2_yz_yyz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 475);

    auto ta2_yz_yyz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 477);

    auto ta2_yz_yyz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 478);

    auto ta2_yz_yyz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 479);

    auto ta2_yz_yzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 481);

    auto ta2_yz_yzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 482);

    auto ta2_yz_yzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 483);

    auto ta2_yz_yzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 484);

    auto ta2_yz_yzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 485);

    auto ta2_yz_yzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 486);

    auto ta2_yz_yzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 487);

    auto ta2_yz_yzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 488);

    auto ta2_yz_yzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 489);

    auto ta2_yz_zzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 490);

    auto ta2_yz_zzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 491);

    auto ta2_yz_zzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 492);

    auto ta2_yz_zzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 493);

    auto ta2_yz_zzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 494);

    auto ta2_yz_zzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 495);

    auto ta2_yz_zzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 496);

    auto ta2_yz_zzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 497);

    auto ta2_yz_zzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 498);

    auto ta2_yz_zzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 499);

    auto ta2_zz_xxx_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 500);

    auto ta2_zz_xxx_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 501);

    auto ta2_zz_xxx_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 502);

    auto ta2_zz_xxx_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 503);

    auto ta2_zz_xxx_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 504);

    auto ta2_zz_xxx_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 505);

    auto ta2_zz_xxx_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 506);

    auto ta2_zz_xxx_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 507);

    auto ta2_zz_xxx_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 508);

    auto ta2_zz_xxx_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 509);

    auto ta2_zz_xxz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 522);

    auto ta2_zz_xxz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 524);

    auto ta2_zz_xxz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 525);

    auto ta2_zz_xyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 531);

    auto ta2_zz_xyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 533);

    auto ta2_zz_xyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 534);

    auto ta2_zz_xyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 536);

    auto ta2_zz_xyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 537);

    auto ta2_zz_xyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 538);

    auto ta2_zz_xzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 552);

    auto ta2_zz_xzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 554);

    auto ta2_zz_xzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 555);

    auto ta2_zz_xzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 557);

    auto ta2_zz_xzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 558);

    auto ta2_zz_xzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 559);

    auto ta2_zz_yyy_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 560);

    auto ta2_zz_yyy_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 561);

    auto ta2_zz_yyy_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 562);

    auto ta2_zz_yyy_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 563);

    auto ta2_zz_yyy_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 564);

    auto ta2_zz_yyy_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 565);

    auto ta2_zz_yyy_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 566);

    auto ta2_zz_yyy_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 567);

    auto ta2_zz_yyy_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 568);

    auto ta2_zz_yyy_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 569);

    auto ta2_zz_yyz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 572);

    auto ta2_zz_yyz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 574);

    auto ta2_zz_yyz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 575);

    auto ta2_zz_yyz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 577);

    auto ta2_zz_yyz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 578);

    auto ta2_zz_yyz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 579);

    auto ta2_zz_yzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 581);

    auto ta2_zz_yzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 582);

    auto ta2_zz_yzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 583);

    auto ta2_zz_yzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 584);

    auto ta2_zz_yzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 585);

    auto ta2_zz_yzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 586);

    auto ta2_zz_yzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 587);

    auto ta2_zz_yzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 588);

    auto ta2_zz_yzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 589);

    auto ta2_zz_zzz_xxx_1 = pbuffer.data(idx_npot_geom_020_1_ff + 590);

    auto ta2_zz_zzz_xxy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 591);

    auto ta2_zz_zzz_xxz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 592);

    auto ta2_zz_zzz_xyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 593);

    auto ta2_zz_zzz_xyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 594);

    auto ta2_zz_zzz_xzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 595);

    auto ta2_zz_zzz_yyy_1 = pbuffer.data(idx_npot_geom_020_1_ff + 596);

    auto ta2_zz_zzz_yyz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 597);

    auto ta2_zz_zzz_yzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 598);

    auto ta2_zz_zzz_zzz_1 = pbuffer.data(idx_npot_geom_020_1_ff + 599);

    // Set up components of auxiliary buffer : FG

    auto ta1_x_xxx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg);

    auto ta1_x_xxx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 1);

    auto ta1_x_xxx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 2);

    auto ta1_x_xxx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 3);

    auto ta1_x_xxx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 4);

    auto ta1_x_xxx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 5);

    auto ta1_x_xxx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 6);

    auto ta1_x_xxx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 7);

    auto ta1_x_xxx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 8);

    auto ta1_x_xxx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 9);

    auto ta1_x_xxx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 10);

    auto ta1_x_xxx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 11);

    auto ta1_x_xxx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 12);

    auto ta1_x_xxx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 13);

    auto ta1_x_xxx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 14);

    auto ta1_x_xxy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 15);

    auto ta1_x_xxy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 16);

    auto ta1_x_xxy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 17);

    auto ta1_x_xxy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 18);

    auto ta1_x_xxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 20);

    auto ta1_x_xxy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 21);

    auto ta1_x_xxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 24);

    auto ta1_x_xxy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 25);

    auto ta1_x_xxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 30);

    auto ta1_x_xxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 31);

    auto ta1_x_xxz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 32);

    auto ta1_x_xxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 33);

    auto ta1_x_xxz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 35);

    auto ta1_x_xxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 36);

    auto ta1_x_xxz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 39);

    auto ta1_x_xxz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 44);

    auto ta1_x_xyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 45);

    auto ta1_x_xyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 46);

    auto ta1_x_xyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 48);

    auto ta1_x_xyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 51);

    auto ta1_x_xyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 55);

    auto ta1_x_xyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 56);

    auto ta1_x_xyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 57);

    auto ta1_x_xyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 58);

    auto ta1_x_xzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 75);

    auto ta1_x_xzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 77);

    auto ta1_x_xzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 80);

    auto ta1_x_xzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 84);

    auto ta1_x_xzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 86);

    auto ta1_x_xzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 87);

    auto ta1_x_xzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 88);

    auto ta1_x_xzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 89);

    auto ta1_x_yyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 90);

    auto ta1_x_yyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 91);

    auto ta1_x_yyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 92);

    auto ta1_x_yyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 93);

    auto ta1_x_yyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 94);

    auto ta1_x_yyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 95);

    auto ta1_x_yyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 96);

    auto ta1_x_yyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 97);

    auto ta1_x_yyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 98);

    auto ta1_x_yyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 99);

    auto ta1_x_yyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 100);

    auto ta1_x_yyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 101);

    auto ta1_x_yyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 102);

    auto ta1_x_yyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 103);

    auto ta1_x_yyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 104);

    auto ta1_x_yyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 106);

    auto ta1_x_yyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 108);

    auto ta1_x_yyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 111);

    auto ta1_x_yyz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 115);

    auto ta1_x_yyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 116);

    auto ta1_x_yyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 117);

    auto ta1_x_yyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 118);

    auto ta1_x_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 119);

    auto ta1_x_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 122);

    auto ta1_x_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 125);

    auto ta1_x_yzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 129);

    auto ta1_x_yzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 130);

    auto ta1_x_yzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 131);

    auto ta1_x_yzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 132);

    auto ta1_x_yzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 133);

    auto ta1_x_yzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 134);

    auto ta1_x_zzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 135);

    auto ta1_x_zzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 136);

    auto ta1_x_zzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 137);

    auto ta1_x_zzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 138);

    auto ta1_x_zzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 139);

    auto ta1_x_zzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 140);

    auto ta1_x_zzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 141);

    auto ta1_x_zzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 142);

    auto ta1_x_zzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 143);

    auto ta1_x_zzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 144);

    auto ta1_x_zzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 145);

    auto ta1_x_zzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 146);

    auto ta1_x_zzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 147);

    auto ta1_x_zzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 148);

    auto ta1_x_zzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 149);

    auto ta1_y_xxx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 150);

    auto ta1_y_xxx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 151);

    auto ta1_y_xxx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 152);

    auto ta1_y_xxx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 153);

    auto ta1_y_xxx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 154);

    auto ta1_y_xxx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 155);

    auto ta1_y_xxx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 156);

    auto ta1_y_xxx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 157);

    auto ta1_y_xxx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 158);

    auto ta1_y_xxx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 159);

    auto ta1_y_xxx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 160);

    auto ta1_y_xxx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 161);

    auto ta1_y_xxx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 162);

    auto ta1_y_xxx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 163);

    auto ta1_y_xxx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 164);

    auto ta1_y_xxy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 165);

    auto ta1_y_xxy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 166);

    auto ta1_y_xxy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 167);

    auto ta1_y_xxy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 168);

    auto ta1_y_xxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 170);

    auto ta1_y_xxy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 171);

    auto ta1_y_xxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 174);

    auto ta1_y_xxy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 175);

    auto ta1_y_xxy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 176);

    auto ta1_y_xxy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 177);

    auto ta1_y_xxy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 178);

    auto ta1_y_xxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 180);

    auto ta1_y_xxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 181);

    auto ta1_y_xxz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 182);

    auto ta1_y_xxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 183);

    auto ta1_y_xxz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 185);

    auto ta1_y_xxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 186);

    auto ta1_y_xxz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 189);

    auto ta1_y_xxz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 194);

    auto ta1_y_xyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 195);

    auto ta1_y_xyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 196);

    auto ta1_y_xyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 198);

    auto ta1_y_xyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 199);

    auto ta1_y_xyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 201);

    auto ta1_y_xyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 202);

    auto ta1_y_xyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 203);

    auto ta1_y_xyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 205);

    auto ta1_y_xyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 206);

    auto ta1_y_xyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 207);

    auto ta1_y_xyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 208);

    auto ta1_y_xyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 209);

    auto ta1_y_xzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 225);

    auto ta1_y_xzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 227);

    auto ta1_y_xzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 230);

    auto ta1_y_xzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 234);

    auto ta1_y_xzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 236);

    auto ta1_y_xzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 237);

    auto ta1_y_xzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 238);

    auto ta1_y_xzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 239);

    auto ta1_y_yyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 240);

    auto ta1_y_yyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 241);

    auto ta1_y_yyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 242);

    auto ta1_y_yyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 243);

    auto ta1_y_yyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 244);

    auto ta1_y_yyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 245);

    auto ta1_y_yyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 246);

    auto ta1_y_yyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 247);

    auto ta1_y_yyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 248);

    auto ta1_y_yyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 249);

    auto ta1_y_yyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 250);

    auto ta1_y_yyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 251);

    auto ta1_y_yyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 252);

    auto ta1_y_yyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 253);

    auto ta1_y_yyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 254);

    auto ta1_y_yyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 256);

    auto ta1_y_yyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 258);

    auto ta1_y_yyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 261);

    auto ta1_y_yyz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 265);

    auto ta1_y_yyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 266);

    auto ta1_y_yyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 267);

    auto ta1_y_yyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 268);

    auto ta1_y_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 269);

    auto ta1_y_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 272);

    auto ta1_y_yzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 274);

    auto ta1_y_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 275);

    auto ta1_y_yzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 277);

    auto ta1_y_yzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 278);

    auto ta1_y_yzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 279);

    auto ta1_y_yzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 280);

    auto ta1_y_yzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 281);

    auto ta1_y_yzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 282);

    auto ta1_y_yzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 283);

    auto ta1_y_yzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 284);

    auto ta1_y_zzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 285);

    auto ta1_y_zzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 286);

    auto ta1_y_zzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 287);

    auto ta1_y_zzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 288);

    auto ta1_y_zzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 289);

    auto ta1_y_zzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 290);

    auto ta1_y_zzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 291);

    auto ta1_y_zzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 292);

    auto ta1_y_zzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 293);

    auto ta1_y_zzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 294);

    auto ta1_y_zzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 295);

    auto ta1_y_zzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 296);

    auto ta1_y_zzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 297);

    auto ta1_y_zzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 298);

    auto ta1_y_zzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 299);

    auto ta1_z_xxx_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 300);

    auto ta1_z_xxx_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 301);

    auto ta1_z_xxx_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 302);

    auto ta1_z_xxx_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 303);

    auto ta1_z_xxx_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 304);

    auto ta1_z_xxx_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 305);

    auto ta1_z_xxx_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 306);

    auto ta1_z_xxx_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 307);

    auto ta1_z_xxx_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 308);

    auto ta1_z_xxx_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 309);

    auto ta1_z_xxx_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 310);

    auto ta1_z_xxx_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 311);

    auto ta1_z_xxx_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 312);

    auto ta1_z_xxx_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 313);

    auto ta1_z_xxx_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 314);

    auto ta1_z_xxy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 315);

    auto ta1_z_xxy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 316);

    auto ta1_z_xxy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 317);

    auto ta1_z_xxy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 318);

    auto ta1_z_xxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 320);

    auto ta1_z_xxy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 321);

    auto ta1_z_xxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 324);

    auto ta1_z_xxy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 325);

    auto ta1_z_xxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 330);

    auto ta1_z_xxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 331);

    auto ta1_z_xxz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 332);

    auto ta1_z_xxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 333);

    auto ta1_z_xxz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 334);

    auto ta1_z_xxz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 335);

    auto ta1_z_xxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 336);

    auto ta1_z_xxz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 337);

    auto ta1_z_xxz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 338);

    auto ta1_z_xxz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 339);

    auto ta1_z_xxz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 341);

    auto ta1_z_xxz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 342);

    auto ta1_z_xxz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 343);

    auto ta1_z_xxz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 344);

    auto ta1_z_xyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 345);

    auto ta1_z_xyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 346);

    auto ta1_z_xyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 348);

    auto ta1_z_xyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 351);

    auto ta1_z_xyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 355);

    auto ta1_z_xyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 356);

    auto ta1_z_xyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 357);

    auto ta1_z_xyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 358);

    auto ta1_z_xzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 375);

    auto ta1_z_xzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 377);

    auto ta1_z_xzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 379);

    auto ta1_z_xzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 380);

    auto ta1_z_xzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 382);

    auto ta1_z_xzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 383);

    auto ta1_z_xzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 384);

    auto ta1_z_xzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 385);

    auto ta1_z_xzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 386);

    auto ta1_z_xzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 387);

    auto ta1_z_xzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 388);

    auto ta1_z_xzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 389);

    auto ta1_z_yyy_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 390);

    auto ta1_z_yyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 391);

    auto ta1_z_yyy_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 392);

    auto ta1_z_yyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 393);

    auto ta1_z_yyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 394);

    auto ta1_z_yyy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 395);

    auto ta1_z_yyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 396);

    auto ta1_z_yyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 397);

    auto ta1_z_yyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 398);

    auto ta1_z_yyy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 399);

    auto ta1_z_yyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 400);

    auto ta1_z_yyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 401);

    auto ta1_z_yyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 402);

    auto ta1_z_yyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 403);

    auto ta1_z_yyy_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 404);

    auto ta1_z_yyz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 406);

    auto ta1_z_yyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 407);

    auto ta1_z_yyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 408);

    auto ta1_z_yyz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 409);

    auto ta1_z_yyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 410);

    auto ta1_z_yyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 411);

    auto ta1_z_yyz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 412);

    auto ta1_z_yyz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 413);

    auto ta1_z_yyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 414);

    auto ta1_z_yyz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 415);

    auto ta1_z_yyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 416);

    auto ta1_z_yyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 417);

    auto ta1_z_yyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 418);

    auto ta1_z_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 419);

    auto ta1_z_yzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 420);

    auto ta1_z_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 422);

    auto ta1_z_yzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 424);

    auto ta1_z_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 425);

    auto ta1_z_yzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 427);

    auto ta1_z_yzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 428);

    auto ta1_z_yzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 429);

    auto ta1_z_yzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 430);

    auto ta1_z_yzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 431);

    auto ta1_z_yzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 432);

    auto ta1_z_yzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 433);

    auto ta1_z_yzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 434);

    auto ta1_z_zzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 435);

    auto ta1_z_zzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 436);

    auto ta1_z_zzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 437);

    auto ta1_z_zzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 438);

    auto ta1_z_zzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 439);

    auto ta1_z_zzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 440);

    auto ta1_z_zzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 441);

    auto ta1_z_zzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 442);

    auto ta1_z_zzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 443);

    auto ta1_z_zzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 444);

    auto ta1_z_zzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 445);

    auto ta1_z_zzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 446);

    auto ta1_z_zzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 447);

    auto ta1_z_zzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 448);

    auto ta1_z_zzz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 449);

    // Set up components of auxiliary buffer : FG

    auto ta2_xx_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg);

    auto ta2_xx_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 1);

    auto ta2_xx_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 2);

    auto ta2_xx_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 3);

    auto ta2_xx_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 4);

    auto ta2_xx_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 5);

    auto ta2_xx_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 6);

    auto ta2_xx_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 7);

    auto ta2_xx_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 8);

    auto ta2_xx_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 9);

    auto ta2_xx_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 10);

    auto ta2_xx_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 11);

    auto ta2_xx_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 12);

    auto ta2_xx_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 13);

    auto ta2_xx_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 14);

    auto ta2_xx_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 15);

    auto ta2_xx_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 16);

    auto ta2_xx_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 17);

    auto ta2_xx_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 18);

    auto ta2_xx_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 19);

    auto ta2_xx_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 20);

    auto ta2_xx_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 21);

    auto ta2_xx_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 22);

    auto ta2_xx_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 23);

    auto ta2_xx_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 24);

    auto ta2_xx_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 25);

    auto ta2_xx_xxy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 29);

    auto ta2_xx_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 30);

    auto ta2_xx_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 31);

    auto ta2_xx_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 32);

    auto ta2_xx_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 33);

    auto ta2_xx_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 34);

    auto ta2_xx_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 35);

    auto ta2_xx_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 36);

    auto ta2_xx_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 37);

    auto ta2_xx_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 38);

    auto ta2_xx_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 39);

    auto ta2_xx_xxz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 40);

    auto ta2_xx_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 41);

    auto ta2_xx_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 42);

    auto ta2_xx_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 43);

    auto ta2_xx_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 44);

    auto ta2_xx_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 45);

    auto ta2_xx_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 46);

    auto ta2_xx_xyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 47);

    auto ta2_xx_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 48);

    auto ta2_xx_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 49);

    auto ta2_xx_xyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 50);

    auto ta2_xx_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 51);

    auto ta2_xx_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 52);

    auto ta2_xx_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 53);

    auto ta2_xx_xyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 54);

    auto ta2_xx_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 55);

    auto ta2_xx_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 56);

    auto ta2_xx_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 57);

    auto ta2_xx_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 58);

    auto ta2_xx_xyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 62);

    auto ta2_xx_xyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 65);

    auto ta2_xx_xyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 69);

    auto ta2_xx_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 75);

    auto ta2_xx_xzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 76);

    auto ta2_xx_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 77);

    auto ta2_xx_xzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 78);

    auto ta2_xx_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 79);

    auto ta2_xx_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 80);

    auto ta2_xx_xzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 81);

    auto ta2_xx_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 82);

    auto ta2_xx_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 83);

    auto ta2_xx_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 84);

    auto ta2_xx_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 86);

    auto ta2_xx_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 87);

    auto ta2_xx_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 88);

    auto ta2_xx_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 89);

    auto ta2_xx_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 90);

    auto ta2_xx_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 91);

    auto ta2_xx_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 92);

    auto ta2_xx_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 93);

    auto ta2_xx_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 94);

    auto ta2_xx_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 95);

    auto ta2_xx_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 96);

    auto ta2_xx_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 97);

    auto ta2_xx_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 98);

    auto ta2_xx_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 99);

    auto ta2_xx_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 100);

    auto ta2_xx_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 101);

    auto ta2_xx_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 102);

    auto ta2_xx_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 103);

    auto ta2_xx_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 104);

    auto ta2_xx_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 106);

    auto ta2_xx_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 107);

    auto ta2_xx_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 108);

    auto ta2_xx_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 110);

    auto ta2_xx_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 111);

    auto ta2_xx_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 114);

    auto ta2_xx_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 115);

    auto ta2_xx_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 116);

    auto ta2_xx_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 117);

    auto ta2_xx_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 118);

    auto ta2_xx_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 119);

    auto ta2_xx_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 120);

    auto ta2_xx_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 122);

    auto ta2_xx_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 124);

    auto ta2_xx_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 125);

    auto ta2_xx_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 127);

    auto ta2_xx_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 128);

    auto ta2_xx_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 129);

    auto ta2_xx_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 130);

    auto ta2_xx_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 131);

    auto ta2_xx_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 132);

    auto ta2_xx_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 133);

    auto ta2_xx_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 134);

    auto ta2_xx_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 135);

    auto ta2_xx_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 136);

    auto ta2_xx_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 137);

    auto ta2_xx_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 138);

    auto ta2_xx_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 139);

    auto ta2_xx_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 140);

    auto ta2_xx_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 141);

    auto ta2_xx_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 142);

    auto ta2_xx_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 143);

    auto ta2_xx_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 144);

    auto ta2_xx_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 145);

    auto ta2_xx_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 146);

    auto ta2_xx_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 147);

    auto ta2_xx_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 148);

    auto ta2_xx_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 149);

    auto ta2_xy_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 150);

    auto ta2_xy_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 151);

    auto ta2_xy_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 152);

    auto ta2_xy_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 153);

    auto ta2_xy_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 154);

    auto ta2_xy_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 155);

    auto ta2_xy_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 156);

    auto ta2_xy_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 157);

    auto ta2_xy_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 158);

    auto ta2_xy_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 159);

    auto ta2_xy_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 160);

    auto ta2_xy_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 161);

    auto ta2_xy_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 162);

    auto ta2_xy_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 163);

    auto ta2_xy_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 164);

    auto ta2_xy_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 165);

    auto ta2_xy_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 166);

    auto ta2_xy_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 167);

    auto ta2_xy_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 168);

    auto ta2_xy_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 169);

    auto ta2_xy_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 170);

    auto ta2_xy_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 171);

    auto ta2_xy_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 172);

    auto ta2_xy_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 173);

    auto ta2_xy_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 174);

    auto ta2_xy_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 175);

    auto ta2_xy_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 176);

    auto ta2_xy_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 177);

    auto ta2_xy_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 178);

    auto ta2_xy_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 180);

    auto ta2_xy_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 181);

    auto ta2_xy_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 182);

    auto ta2_xy_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 183);

    auto ta2_xy_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 184);

    auto ta2_xy_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 185);

    auto ta2_xy_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 186);

    auto ta2_xy_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 187);

    auto ta2_xy_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 188);

    auto ta2_xy_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 189);

    auto ta2_xy_xxz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 190);

    auto ta2_xy_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 194);

    auto ta2_xy_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 195);

    auto ta2_xy_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 196);

    auto ta2_xy_xyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 197);

    auto ta2_xy_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 198);

    auto ta2_xy_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 199);

    auto ta2_xy_xyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 200);

    auto ta2_xy_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 201);

    auto ta2_xy_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 202);

    auto ta2_xy_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 203);

    auto ta2_xy_xyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 204);

    auto ta2_xy_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 205);

    auto ta2_xy_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 206);

    auto ta2_xy_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 207);

    auto ta2_xy_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 208);

    auto ta2_xy_xyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 209);

    auto ta2_xy_xyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 211);

    auto ta2_xy_xyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 213);

    auto ta2_xy_xyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 216);

    auto ta2_xy_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 225);

    auto ta2_xy_xzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 226);

    auto ta2_xy_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 227);

    auto ta2_xy_xzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 228);

    auto ta2_xy_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 230);

    auto ta2_xy_xzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 231);

    auto ta2_xy_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 234);

    auto ta2_xy_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 236);

    auto ta2_xy_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 237);

    auto ta2_xy_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 238);

    auto ta2_xy_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 239);

    auto ta2_xy_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 240);

    auto ta2_xy_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 241);

    auto ta2_xy_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 242);

    auto ta2_xy_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 243);

    auto ta2_xy_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 244);

    auto ta2_xy_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 245);

    auto ta2_xy_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 246);

    auto ta2_xy_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 247);

    auto ta2_xy_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 248);

    auto ta2_xy_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 249);

    auto ta2_xy_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 250);

    auto ta2_xy_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 251);

    auto ta2_xy_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 252);

    auto ta2_xy_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 253);

    auto ta2_xy_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 254);

    auto ta2_xy_yyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 255);

    auto ta2_xy_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 256);

    auto ta2_xy_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 258);

    auto ta2_xy_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 259);

    auto ta2_xy_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 261);

    auto ta2_xy_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 262);

    auto ta2_xy_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 263);

    auto ta2_xy_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 265);

    auto ta2_xy_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 266);

    auto ta2_xy_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 267);

    auto ta2_xy_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 268);

    auto ta2_xy_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 269);

    auto ta2_xy_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 271);

    auto ta2_xy_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 272);

    auto ta2_xy_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 273);

    auto ta2_xy_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 274);

    auto ta2_xy_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 275);

    auto ta2_xy_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 276);

    auto ta2_xy_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 277);

    auto ta2_xy_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 278);

    auto ta2_xy_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 279);

    auto ta2_xy_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 280);

    auto ta2_xy_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 281);

    auto ta2_xy_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 282);

    auto ta2_xy_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 283);

    auto ta2_xy_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 284);

    auto ta2_xy_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 285);

    auto ta2_xy_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 286);

    auto ta2_xy_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 287);

    auto ta2_xy_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 288);

    auto ta2_xy_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 289);

    auto ta2_xy_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 290);

    auto ta2_xy_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 291);

    auto ta2_xy_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 292);

    auto ta2_xy_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 293);

    auto ta2_xy_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 294);

    auto ta2_xy_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 295);

    auto ta2_xy_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 296);

    auto ta2_xy_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 297);

    auto ta2_xy_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 298);

    auto ta2_xy_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 299);

    auto ta2_xz_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 300);

    auto ta2_xz_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 301);

    auto ta2_xz_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 302);

    auto ta2_xz_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 303);

    auto ta2_xz_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 304);

    auto ta2_xz_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 305);

    auto ta2_xz_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 306);

    auto ta2_xz_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 307);

    auto ta2_xz_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 308);

    auto ta2_xz_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 309);

    auto ta2_xz_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 310);

    auto ta2_xz_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 311);

    auto ta2_xz_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 312);

    auto ta2_xz_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 313);

    auto ta2_xz_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 314);

    auto ta2_xz_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 315);

    auto ta2_xz_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 316);

    auto ta2_xz_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 317);

    auto ta2_xz_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 318);

    auto ta2_xz_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 319);

    auto ta2_xz_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 320);

    auto ta2_xz_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 321);

    auto ta2_xz_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 322);

    auto ta2_xz_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 323);

    auto ta2_xz_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 324);

    auto ta2_xz_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 325);

    auto ta2_xz_xxy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 329);

    auto ta2_xz_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 330);

    auto ta2_xz_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 331);

    auto ta2_xz_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 332);

    auto ta2_xz_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 333);

    auto ta2_xz_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 334);

    auto ta2_xz_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 335);

    auto ta2_xz_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 336);

    auto ta2_xz_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 337);

    auto ta2_xz_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 338);

    auto ta2_xz_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 339);

    auto ta2_xz_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 341);

    auto ta2_xz_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 342);

    auto ta2_xz_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 343);

    auto ta2_xz_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 344);

    auto ta2_xz_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 345);

    auto ta2_xz_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 346);

    auto ta2_xz_xyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 347);

    auto ta2_xz_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 348);

    auto ta2_xz_xyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 350);

    auto ta2_xz_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 351);

    auto ta2_xz_xyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 354);

    auto ta2_xz_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 355);

    auto ta2_xz_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 356);

    auto ta2_xz_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 357);

    auto ta2_xz_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 358);

    auto ta2_xz_xyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 362);

    auto ta2_xz_xyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 365);

    auto ta2_xz_xyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 369);

    auto ta2_xz_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 375);

    auto ta2_xz_xzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 376);

    auto ta2_xz_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 377);

    auto ta2_xz_xzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 378);

    auto ta2_xz_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 379);

    auto ta2_xz_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 380);

    auto ta2_xz_xzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 381);

    auto ta2_xz_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 382);

    auto ta2_xz_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 383);

    auto ta2_xz_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 384);

    auto ta2_xz_xzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 385);

    auto ta2_xz_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 386);

    auto ta2_xz_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 387);

    auto ta2_xz_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 388);

    auto ta2_xz_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 389);

    auto ta2_xz_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 390);

    auto ta2_xz_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 391);

    auto ta2_xz_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 392);

    auto ta2_xz_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 393);

    auto ta2_xz_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 394);

    auto ta2_xz_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 395);

    auto ta2_xz_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 396);

    auto ta2_xz_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 397);

    auto ta2_xz_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 398);

    auto ta2_xz_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 399);

    auto ta2_xz_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 400);

    auto ta2_xz_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 401);

    auto ta2_xz_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 402);

    auto ta2_xz_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 403);

    auto ta2_xz_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 404);

    auto ta2_xz_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 406);

    auto ta2_xz_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 407);

    auto ta2_xz_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 408);

    auto ta2_xz_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 409);

    auto ta2_xz_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 410);

    auto ta2_xz_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 411);

    auto ta2_xz_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 412);

    auto ta2_xz_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 413);

    auto ta2_xz_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 414);

    auto ta2_xz_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 415);

    auto ta2_xz_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 416);

    auto ta2_xz_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 417);

    auto ta2_xz_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 418);

    auto ta2_xz_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 419);

    auto ta2_xz_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 420);

    auto ta2_xz_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 422);

    auto ta2_xz_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 424);

    auto ta2_xz_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 425);

    auto ta2_xz_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 427);

    auto ta2_xz_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 428);

    auto ta2_xz_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 429);

    auto ta2_xz_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 430);

    auto ta2_xz_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 431);

    auto ta2_xz_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 432);

    auto ta2_xz_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 433);

    auto ta2_xz_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 434);

    auto ta2_xz_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 435);

    auto ta2_xz_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 436);

    auto ta2_xz_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 437);

    auto ta2_xz_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 438);

    auto ta2_xz_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 439);

    auto ta2_xz_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 440);

    auto ta2_xz_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 441);

    auto ta2_xz_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 442);

    auto ta2_xz_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 443);

    auto ta2_xz_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 444);

    auto ta2_xz_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 445);

    auto ta2_xz_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 446);

    auto ta2_xz_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 447);

    auto ta2_xz_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 448);

    auto ta2_xz_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 449);

    auto ta2_yy_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 450);

    auto ta2_yy_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 451);

    auto ta2_yy_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 452);

    auto ta2_yy_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 453);

    auto ta2_yy_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 454);

    auto ta2_yy_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 455);

    auto ta2_yy_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 456);

    auto ta2_yy_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 457);

    auto ta2_yy_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 458);

    auto ta2_yy_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 459);

    auto ta2_yy_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 460);

    auto ta2_yy_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 461);

    auto ta2_yy_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 462);

    auto ta2_yy_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 463);

    auto ta2_yy_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 464);

    auto ta2_yy_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 465);

    auto ta2_yy_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 466);

    auto ta2_yy_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 467);

    auto ta2_yy_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 468);

    auto ta2_yy_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 469);

    auto ta2_yy_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 470);

    auto ta2_yy_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 471);

    auto ta2_yy_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 472);

    auto ta2_yy_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 473);

    auto ta2_yy_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 474);

    auto ta2_yy_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 475);

    auto ta2_yy_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 476);

    auto ta2_yy_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 477);

    auto ta2_yy_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 478);

    auto ta2_yy_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 480);

    auto ta2_yy_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 481);

    auto ta2_yy_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 482);

    auto ta2_yy_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 483);

    auto ta2_yy_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 485);

    auto ta2_yy_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 486);

    auto ta2_yy_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 489);

    auto ta2_yy_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 491);

    auto ta2_yy_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 492);

    auto ta2_yy_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 493);

    auto ta2_yy_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 494);

    auto ta2_yy_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 495);

    auto ta2_yy_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 496);

    auto ta2_yy_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 498);

    auto ta2_yy_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 499);

    auto ta2_yy_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 501);

    auto ta2_yy_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 502);

    auto ta2_yy_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 503);

    auto ta2_yy_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 505);

    auto ta2_yy_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 506);

    auto ta2_yy_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 507);

    auto ta2_yy_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 508);

    auto ta2_yy_xyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 509);

    auto ta2_yy_xyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 521);

    auto ta2_yy_xyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 522);

    auto ta2_yy_xyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 523);

    auto ta2_yy_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 525);

    auto ta2_yy_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 527);

    auto ta2_yy_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 529);

    auto ta2_yy_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 530);

    auto ta2_yy_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 532);

    auto ta2_yy_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 533);

    auto ta2_yy_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 534);

    auto ta2_yy_xzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 535);

    auto ta2_yy_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 536);

    auto ta2_yy_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 537);

    auto ta2_yy_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 538);

    auto ta2_yy_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 539);

    auto ta2_yy_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 540);

    auto ta2_yy_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 541);

    auto ta2_yy_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 542);

    auto ta2_yy_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 543);

    auto ta2_yy_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 544);

    auto ta2_yy_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 545);

    auto ta2_yy_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 546);

    auto ta2_yy_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 547);

    auto ta2_yy_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 548);

    auto ta2_yy_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 549);

    auto ta2_yy_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 550);

    auto ta2_yy_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 551);

    auto ta2_yy_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 552);

    auto ta2_yy_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 553);

    auto ta2_yy_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 554);

    auto ta2_yy_yyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 555);

    auto ta2_yy_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 556);

    auto ta2_yy_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 557);

    auto ta2_yy_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 558);

    auto ta2_yy_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 559);

    auto ta2_yy_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 560);

    auto ta2_yy_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 561);

    auto ta2_yy_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 562);

    auto ta2_yy_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 563);

    auto ta2_yy_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 564);

    auto ta2_yy_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 565);

    auto ta2_yy_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 566);

    auto ta2_yy_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 567);

    auto ta2_yy_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 568);

    auto ta2_yy_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 569);

    auto ta2_yy_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 571);

    auto ta2_yy_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 572);

    auto ta2_yy_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 573);

    auto ta2_yy_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 574);

    auto ta2_yy_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 575);

    auto ta2_yy_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 576);

    auto ta2_yy_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 577);

    auto ta2_yy_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 578);

    auto ta2_yy_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 579);

    auto ta2_yy_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 580);

    auto ta2_yy_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 581);

    auto ta2_yy_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 582);

    auto ta2_yy_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 583);

    auto ta2_yy_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 584);

    auto ta2_yy_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 585);

    auto ta2_yy_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 586);

    auto ta2_yy_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 587);

    auto ta2_yy_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 588);

    auto ta2_yy_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 589);

    auto ta2_yy_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 590);

    auto ta2_yy_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 591);

    auto ta2_yy_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 592);

    auto ta2_yy_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 593);

    auto ta2_yy_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 594);

    auto ta2_yy_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 595);

    auto ta2_yy_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 596);

    auto ta2_yy_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 597);

    auto ta2_yy_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 598);

    auto ta2_yy_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 599);

    auto ta2_yz_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 600);

    auto ta2_yz_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 601);

    auto ta2_yz_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 602);

    auto ta2_yz_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 603);

    auto ta2_yz_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 604);

    auto ta2_yz_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 605);

    auto ta2_yz_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 606);

    auto ta2_yz_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 607);

    auto ta2_yz_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 608);

    auto ta2_yz_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 609);

    auto ta2_yz_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 610);

    auto ta2_yz_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 611);

    auto ta2_yz_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 612);

    auto ta2_yz_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 613);

    auto ta2_yz_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 614);

    auto ta2_yz_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 615);

    auto ta2_yz_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 616);

    auto ta2_yz_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 617);

    auto ta2_yz_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 618);

    auto ta2_yz_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 620);

    auto ta2_yz_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 621);

    auto ta2_yz_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 624);

    auto ta2_yz_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 625);

    auto ta2_yz_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 626);

    auto ta2_yz_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 627);

    auto ta2_yz_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 628);

    auto ta2_yz_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 630);

    auto ta2_yz_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 631);

    auto ta2_yz_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 632);

    auto ta2_yz_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 633);

    auto ta2_yz_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 634);

    auto ta2_yz_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 635);

    auto ta2_yz_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 636);

    auto ta2_yz_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 637);

    auto ta2_yz_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 638);

    auto ta2_yz_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 639);

    auto ta2_yz_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 641);

    auto ta2_yz_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 642);

    auto ta2_yz_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 643);

    auto ta2_yz_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 644);

    auto ta2_yz_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 645);

    auto ta2_yz_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 646);

    auto ta2_yz_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 648);

    auto ta2_yz_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 649);

    auto ta2_yz_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 651);

    auto ta2_yz_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 652);

    auto ta2_yz_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 653);

    auto ta2_yz_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 655);

    auto ta2_yz_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 656);

    auto ta2_yz_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 657);

    auto ta2_yz_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 658);

    auto ta2_yz_xyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 659);

    auto ta2_yz_xyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 671);

    auto ta2_yz_xyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 672);

    auto ta2_yz_xyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 673);

    auto ta2_yz_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 675);

    auto ta2_yz_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 677);

    auto ta2_yz_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 679);

    auto ta2_yz_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 680);

    auto ta2_yz_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 682);

    auto ta2_yz_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 683);

    auto ta2_yz_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 684);

    auto ta2_yz_xzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 685);

    auto ta2_yz_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 686);

    auto ta2_yz_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 687);

    auto ta2_yz_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 688);

    auto ta2_yz_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 689);

    auto ta2_yz_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 690);

    auto ta2_yz_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 691);

    auto ta2_yz_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 692);

    auto ta2_yz_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 693);

    auto ta2_yz_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 694);

    auto ta2_yz_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 695);

    auto ta2_yz_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 696);

    auto ta2_yz_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 697);

    auto ta2_yz_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 698);

    auto ta2_yz_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 699);

    auto ta2_yz_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 700);

    auto ta2_yz_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 701);

    auto ta2_yz_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 702);

    auto ta2_yz_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 703);

    auto ta2_yz_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 704);

    auto ta2_yz_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 706);

    auto ta2_yz_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 707);

    auto ta2_yz_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 708);

    auto ta2_yz_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 709);

    auto ta2_yz_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 710);

    auto ta2_yz_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 711);

    auto ta2_yz_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 712);

    auto ta2_yz_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 713);

    auto ta2_yz_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 714);

    auto ta2_yz_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 715);

    auto ta2_yz_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 716);

    auto ta2_yz_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 717);

    auto ta2_yz_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 718);

    auto ta2_yz_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 719);

    auto ta2_yz_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 720);

    auto ta2_yz_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 721);

    auto ta2_yz_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 722);

    auto ta2_yz_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 723);

    auto ta2_yz_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 724);

    auto ta2_yz_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 725);

    auto ta2_yz_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 726);

    auto ta2_yz_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 727);

    auto ta2_yz_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 728);

    auto ta2_yz_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 729);

    auto ta2_yz_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 730);

    auto ta2_yz_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 731);

    auto ta2_yz_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 732);

    auto ta2_yz_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 733);

    auto ta2_yz_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 734);

    auto ta2_yz_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 735);

    auto ta2_yz_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 736);

    auto ta2_yz_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 737);

    auto ta2_yz_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 738);

    auto ta2_yz_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 739);

    auto ta2_yz_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 740);

    auto ta2_yz_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 741);

    auto ta2_yz_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 742);

    auto ta2_yz_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 743);

    auto ta2_yz_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 744);

    auto ta2_yz_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 745);

    auto ta2_yz_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 746);

    auto ta2_yz_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 747);

    auto ta2_yz_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 748);

    auto ta2_yz_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 749);

    auto ta2_zz_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 750);

    auto ta2_zz_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 751);

    auto ta2_zz_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 752);

    auto ta2_zz_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 753);

    auto ta2_zz_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 754);

    auto ta2_zz_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 755);

    auto ta2_zz_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 756);

    auto ta2_zz_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 757);

    auto ta2_zz_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 758);

    auto ta2_zz_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 759);

    auto ta2_zz_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 760);

    auto ta2_zz_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 761);

    auto ta2_zz_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 762);

    auto ta2_zz_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 763);

    auto ta2_zz_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 764);

    auto ta2_zz_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 765);

    auto ta2_zz_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 766);

    auto ta2_zz_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 767);

    auto ta2_zz_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 768);

    auto ta2_zz_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 770);

    auto ta2_zz_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 771);

    auto ta2_zz_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 774);

    auto ta2_zz_xxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 775);

    auto ta2_zz_xxy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 776);

    auto ta2_zz_xxy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 777);

    auto ta2_zz_xxy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 778);

    auto ta2_zz_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 780);

    auto ta2_zz_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 781);

    auto ta2_zz_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 782);

    auto ta2_zz_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 783);

    auto ta2_zz_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 784);

    auto ta2_zz_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 785);

    auto ta2_zz_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 786);

    auto ta2_zz_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 787);

    auto ta2_zz_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 788);

    auto ta2_zz_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 789);

    auto ta2_zz_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 791);

    auto ta2_zz_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 792);

    auto ta2_zz_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 793);

    auto ta2_zz_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 794);

    auto ta2_zz_xyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 795);

    auto ta2_zz_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 796);

    auto ta2_zz_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 798);

    auto ta2_zz_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 799);

    auto ta2_zz_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 801);

    auto ta2_zz_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 802);

    auto ta2_zz_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 803);

    auto ta2_zz_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 805);

    auto ta2_zz_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 806);

    auto ta2_zz_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 807);

    auto ta2_zz_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 808);

    auto ta2_zz_xyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 809);

    auto ta2_zz_xyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 821);

    auto ta2_zz_xyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 822);

    auto ta2_zz_xyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 823);

    auto ta2_zz_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 825);

    auto ta2_zz_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 827);

    auto ta2_zz_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 829);

    auto ta2_zz_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 830);

    auto ta2_zz_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 832);

    auto ta2_zz_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 833);

    auto ta2_zz_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 834);

    auto ta2_zz_xzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 835);

    auto ta2_zz_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 836);

    auto ta2_zz_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 837);

    auto ta2_zz_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 838);

    auto ta2_zz_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 839);

    auto ta2_zz_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 840);

    auto ta2_zz_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 841);

    auto ta2_zz_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 842);

    auto ta2_zz_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 843);

    auto ta2_zz_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 844);

    auto ta2_zz_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 845);

    auto ta2_zz_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 846);

    auto ta2_zz_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 847);

    auto ta2_zz_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 848);

    auto ta2_zz_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 849);

    auto ta2_zz_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 850);

    auto ta2_zz_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 851);

    auto ta2_zz_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 852);

    auto ta2_zz_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 853);

    auto ta2_zz_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 854);

    auto ta2_zz_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 856);

    auto ta2_zz_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 857);

    auto ta2_zz_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 858);

    auto ta2_zz_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 859);

    auto ta2_zz_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 860);

    auto ta2_zz_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 861);

    auto ta2_zz_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 862);

    auto ta2_zz_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 863);

    auto ta2_zz_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 864);

    auto ta2_zz_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 865);

    auto ta2_zz_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 866);

    auto ta2_zz_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 867);

    auto ta2_zz_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 868);

    auto ta2_zz_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 869);

    auto ta2_zz_yzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 870);

    auto ta2_zz_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 871);

    auto ta2_zz_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 872);

    auto ta2_zz_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 873);

    auto ta2_zz_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 874);

    auto ta2_zz_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 875);

    auto ta2_zz_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 876);

    auto ta2_zz_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 877);

    auto ta2_zz_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 878);

    auto ta2_zz_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 879);

    auto ta2_zz_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 880);

    auto ta2_zz_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 881);

    auto ta2_zz_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 882);

    auto ta2_zz_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 883);

    auto ta2_zz_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 884);

    auto ta2_zz_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_fg + 885);

    auto ta2_zz_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 886);

    auto ta2_zz_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 887);

    auto ta2_zz_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 888);

    auto ta2_zz_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 889);

    auto ta2_zz_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 890);

    auto ta2_zz_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 891);

    auto ta2_zz_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 892);

    auto ta2_zz_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 893);

    auto ta2_zz_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 894);

    auto ta2_zz_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_fg + 895);

    auto ta2_zz_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 896);

    auto ta2_zz_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 897);

    auto ta2_zz_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 898);

    auto ta2_zz_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_fg + 899);

    // Set up components of auxiliary buffer : FG

    auto ta2_xx_xxx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg);

    auto ta2_xx_xxx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 1);

    auto ta2_xx_xxx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 2);

    auto ta2_xx_xxx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 3);

    auto ta2_xx_xxx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 4);

    auto ta2_xx_xxx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 5);

    auto ta2_xx_xxx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 6);

    auto ta2_xx_xxx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 7);

    auto ta2_xx_xxx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 8);

    auto ta2_xx_xxx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 9);

    auto ta2_xx_xxx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 10);

    auto ta2_xx_xxx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 11);

    auto ta2_xx_xxx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 12);

    auto ta2_xx_xxx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 13);

    auto ta2_xx_xxx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 14);

    auto ta2_xx_xxy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 15);

    auto ta2_xx_xxy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 16);

    auto ta2_xx_xxy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 17);

    auto ta2_xx_xxy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 18);

    auto ta2_xx_xxy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 19);

    auto ta2_xx_xxy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 20);

    auto ta2_xx_xxy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 21);

    auto ta2_xx_xxy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 22);

    auto ta2_xx_xxy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 23);

    auto ta2_xx_xxy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 24);

    auto ta2_xx_xxy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 25);

    auto ta2_xx_xxy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 29);

    auto ta2_xx_xxz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 30);

    auto ta2_xx_xxz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 31);

    auto ta2_xx_xxz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 32);

    auto ta2_xx_xxz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 33);

    auto ta2_xx_xxz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 34);

    auto ta2_xx_xxz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 35);

    auto ta2_xx_xxz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 36);

    auto ta2_xx_xxz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 37);

    auto ta2_xx_xxz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 38);

    auto ta2_xx_xxz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 39);

    auto ta2_xx_xxz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 40);

    auto ta2_xx_xxz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 41);

    auto ta2_xx_xxz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 42);

    auto ta2_xx_xxz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 43);

    auto ta2_xx_xxz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 44);

    auto ta2_xx_xyy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 45);

    auto ta2_xx_xyy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 46);

    auto ta2_xx_xyy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 47);

    auto ta2_xx_xyy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 48);

    auto ta2_xx_xyy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 49);

    auto ta2_xx_xyy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 50);

    auto ta2_xx_xyy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 51);

    auto ta2_xx_xyy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 52);

    auto ta2_xx_xyy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 53);

    auto ta2_xx_xyy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 54);

    auto ta2_xx_xyy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 55);

    auto ta2_xx_xyy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 56);

    auto ta2_xx_xyy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 57);

    auto ta2_xx_xyy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 58);

    auto ta2_xx_xyz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 62);

    auto ta2_xx_xyz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 65);

    auto ta2_xx_xyz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 69);

    auto ta2_xx_xzz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 75);

    auto ta2_xx_xzz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 76);

    auto ta2_xx_xzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 77);

    auto ta2_xx_xzz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 78);

    auto ta2_xx_xzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 79);

    auto ta2_xx_xzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 80);

    auto ta2_xx_xzz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 81);

    auto ta2_xx_xzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 82);

    auto ta2_xx_xzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 83);

    auto ta2_xx_xzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 84);

    auto ta2_xx_xzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 86);

    auto ta2_xx_xzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 87);

    auto ta2_xx_xzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 88);

    auto ta2_xx_xzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 89);

    auto ta2_xx_yyy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 90);

    auto ta2_xx_yyy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 91);

    auto ta2_xx_yyy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 92);

    auto ta2_xx_yyy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 93);

    auto ta2_xx_yyy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 94);

    auto ta2_xx_yyy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 95);

    auto ta2_xx_yyy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 96);

    auto ta2_xx_yyy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 97);

    auto ta2_xx_yyy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 98);

    auto ta2_xx_yyy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 99);

    auto ta2_xx_yyy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 100);

    auto ta2_xx_yyy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 101);

    auto ta2_xx_yyy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 102);

    auto ta2_xx_yyy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 103);

    auto ta2_xx_yyy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 104);

    auto ta2_xx_yyz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 106);

    auto ta2_xx_yyz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 107);

    auto ta2_xx_yyz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 108);

    auto ta2_xx_yyz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 110);

    auto ta2_xx_yyz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 111);

    auto ta2_xx_yyz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 114);

    auto ta2_xx_yyz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 115);

    auto ta2_xx_yyz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 116);

    auto ta2_xx_yyz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 117);

    auto ta2_xx_yyz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 118);

    auto ta2_xx_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 119);

    auto ta2_xx_yzz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 120);

    auto ta2_xx_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 122);

    auto ta2_xx_yzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 124);

    auto ta2_xx_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 125);

    auto ta2_xx_yzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 127);

    auto ta2_xx_yzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 128);

    auto ta2_xx_yzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 129);

    auto ta2_xx_yzz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 130);

    auto ta2_xx_yzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 131);

    auto ta2_xx_yzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 132);

    auto ta2_xx_yzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 133);

    auto ta2_xx_yzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 134);

    auto ta2_xx_zzz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 135);

    auto ta2_xx_zzz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 136);

    auto ta2_xx_zzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 137);

    auto ta2_xx_zzz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 138);

    auto ta2_xx_zzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 139);

    auto ta2_xx_zzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 140);

    auto ta2_xx_zzz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 141);

    auto ta2_xx_zzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 142);

    auto ta2_xx_zzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 143);

    auto ta2_xx_zzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 144);

    auto ta2_xx_zzz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 145);

    auto ta2_xx_zzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 146);

    auto ta2_xx_zzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 147);

    auto ta2_xx_zzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 148);

    auto ta2_xx_zzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 149);

    auto ta2_xy_xxx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 150);

    auto ta2_xy_xxx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 151);

    auto ta2_xy_xxx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 152);

    auto ta2_xy_xxx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 153);

    auto ta2_xy_xxx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 154);

    auto ta2_xy_xxx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 155);

    auto ta2_xy_xxx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 156);

    auto ta2_xy_xxx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 157);

    auto ta2_xy_xxx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 158);

    auto ta2_xy_xxx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 159);

    auto ta2_xy_xxx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 160);

    auto ta2_xy_xxx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 161);

    auto ta2_xy_xxx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 162);

    auto ta2_xy_xxx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 163);

    auto ta2_xy_xxx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 164);

    auto ta2_xy_xxy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 165);

    auto ta2_xy_xxy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 166);

    auto ta2_xy_xxy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 167);

    auto ta2_xy_xxy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 168);

    auto ta2_xy_xxy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 169);

    auto ta2_xy_xxy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 170);

    auto ta2_xy_xxy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 171);

    auto ta2_xy_xxy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 172);

    auto ta2_xy_xxy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 173);

    auto ta2_xy_xxy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 174);

    auto ta2_xy_xxy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 175);

    auto ta2_xy_xxy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 176);

    auto ta2_xy_xxy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 177);

    auto ta2_xy_xxy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 178);

    auto ta2_xy_xxz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 180);

    auto ta2_xy_xxz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 181);

    auto ta2_xy_xxz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 182);

    auto ta2_xy_xxz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 183);

    auto ta2_xy_xxz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 184);

    auto ta2_xy_xxz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 185);

    auto ta2_xy_xxz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 186);

    auto ta2_xy_xxz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 187);

    auto ta2_xy_xxz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 188);

    auto ta2_xy_xxz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 189);

    auto ta2_xy_xxz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 190);

    auto ta2_xy_xxz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 194);

    auto ta2_xy_xyy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 195);

    auto ta2_xy_xyy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 196);

    auto ta2_xy_xyy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 197);

    auto ta2_xy_xyy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 198);

    auto ta2_xy_xyy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 199);

    auto ta2_xy_xyy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 200);

    auto ta2_xy_xyy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 201);

    auto ta2_xy_xyy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 202);

    auto ta2_xy_xyy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 203);

    auto ta2_xy_xyy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 204);

    auto ta2_xy_xyy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 205);

    auto ta2_xy_xyy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 206);

    auto ta2_xy_xyy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 207);

    auto ta2_xy_xyy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 208);

    auto ta2_xy_xyy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 209);

    auto ta2_xy_xyz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 211);

    auto ta2_xy_xyz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 213);

    auto ta2_xy_xyz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 216);

    auto ta2_xy_xzz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 225);

    auto ta2_xy_xzz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 226);

    auto ta2_xy_xzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 227);

    auto ta2_xy_xzz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 228);

    auto ta2_xy_xzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 230);

    auto ta2_xy_xzz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 231);

    auto ta2_xy_xzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 234);

    auto ta2_xy_xzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 236);

    auto ta2_xy_xzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 237);

    auto ta2_xy_xzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 238);

    auto ta2_xy_xzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 239);

    auto ta2_xy_yyy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 240);

    auto ta2_xy_yyy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 241);

    auto ta2_xy_yyy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 242);

    auto ta2_xy_yyy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 243);

    auto ta2_xy_yyy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 244);

    auto ta2_xy_yyy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 245);

    auto ta2_xy_yyy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 246);

    auto ta2_xy_yyy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 247);

    auto ta2_xy_yyy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 248);

    auto ta2_xy_yyy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 249);

    auto ta2_xy_yyy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 250);

    auto ta2_xy_yyy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 251);

    auto ta2_xy_yyy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 252);

    auto ta2_xy_yyy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 253);

    auto ta2_xy_yyy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 254);

    auto ta2_xy_yyz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 255);

    auto ta2_xy_yyz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 256);

    auto ta2_xy_yyz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 258);

    auto ta2_xy_yyz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 259);

    auto ta2_xy_yyz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 261);

    auto ta2_xy_yyz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 262);

    auto ta2_xy_yyz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 263);

    auto ta2_xy_yyz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 265);

    auto ta2_xy_yyz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 266);

    auto ta2_xy_yyz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 267);

    auto ta2_xy_yyz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 268);

    auto ta2_xy_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 269);

    auto ta2_xy_yzz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 271);

    auto ta2_xy_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 272);

    auto ta2_xy_yzz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 273);

    auto ta2_xy_yzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 274);

    auto ta2_xy_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 275);

    auto ta2_xy_yzz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 276);

    auto ta2_xy_yzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 277);

    auto ta2_xy_yzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 278);

    auto ta2_xy_yzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 279);

    auto ta2_xy_yzz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 280);

    auto ta2_xy_yzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 281);

    auto ta2_xy_yzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 282);

    auto ta2_xy_yzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 283);

    auto ta2_xy_yzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 284);

    auto ta2_xy_zzz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 285);

    auto ta2_xy_zzz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 286);

    auto ta2_xy_zzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 287);

    auto ta2_xy_zzz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 288);

    auto ta2_xy_zzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 289);

    auto ta2_xy_zzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 290);

    auto ta2_xy_zzz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 291);

    auto ta2_xy_zzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 292);

    auto ta2_xy_zzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 293);

    auto ta2_xy_zzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 294);

    auto ta2_xy_zzz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 295);

    auto ta2_xy_zzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 296);

    auto ta2_xy_zzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 297);

    auto ta2_xy_zzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 298);

    auto ta2_xy_zzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 299);

    auto ta2_xz_xxx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 300);

    auto ta2_xz_xxx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 301);

    auto ta2_xz_xxx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 302);

    auto ta2_xz_xxx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 303);

    auto ta2_xz_xxx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 304);

    auto ta2_xz_xxx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 305);

    auto ta2_xz_xxx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 306);

    auto ta2_xz_xxx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 307);

    auto ta2_xz_xxx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 308);

    auto ta2_xz_xxx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 309);

    auto ta2_xz_xxx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 310);

    auto ta2_xz_xxx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 311);

    auto ta2_xz_xxx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 312);

    auto ta2_xz_xxx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 313);

    auto ta2_xz_xxx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 314);

    auto ta2_xz_xxy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 315);

    auto ta2_xz_xxy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 316);

    auto ta2_xz_xxy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 317);

    auto ta2_xz_xxy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 318);

    auto ta2_xz_xxy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 319);

    auto ta2_xz_xxy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 320);

    auto ta2_xz_xxy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 321);

    auto ta2_xz_xxy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 322);

    auto ta2_xz_xxy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 323);

    auto ta2_xz_xxy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 324);

    auto ta2_xz_xxy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 325);

    auto ta2_xz_xxy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 329);

    auto ta2_xz_xxz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 330);

    auto ta2_xz_xxz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 331);

    auto ta2_xz_xxz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 332);

    auto ta2_xz_xxz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 333);

    auto ta2_xz_xxz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 334);

    auto ta2_xz_xxz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 335);

    auto ta2_xz_xxz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 336);

    auto ta2_xz_xxz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 337);

    auto ta2_xz_xxz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 338);

    auto ta2_xz_xxz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 339);

    auto ta2_xz_xxz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 341);

    auto ta2_xz_xxz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 342);

    auto ta2_xz_xxz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 343);

    auto ta2_xz_xxz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 344);

    auto ta2_xz_xyy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 345);

    auto ta2_xz_xyy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 346);

    auto ta2_xz_xyy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 347);

    auto ta2_xz_xyy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 348);

    auto ta2_xz_xyy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 350);

    auto ta2_xz_xyy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 351);

    auto ta2_xz_xyy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 354);

    auto ta2_xz_xyy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 355);

    auto ta2_xz_xyy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 356);

    auto ta2_xz_xyy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 357);

    auto ta2_xz_xyy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 358);

    auto ta2_xz_xyz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 362);

    auto ta2_xz_xyz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 365);

    auto ta2_xz_xyz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 369);

    auto ta2_xz_xzz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 375);

    auto ta2_xz_xzz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 376);

    auto ta2_xz_xzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 377);

    auto ta2_xz_xzz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 378);

    auto ta2_xz_xzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 379);

    auto ta2_xz_xzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 380);

    auto ta2_xz_xzz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 381);

    auto ta2_xz_xzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 382);

    auto ta2_xz_xzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 383);

    auto ta2_xz_xzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 384);

    auto ta2_xz_xzz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 385);

    auto ta2_xz_xzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 386);

    auto ta2_xz_xzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 387);

    auto ta2_xz_xzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 388);

    auto ta2_xz_xzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 389);

    auto ta2_xz_yyy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 390);

    auto ta2_xz_yyy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 391);

    auto ta2_xz_yyy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 392);

    auto ta2_xz_yyy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 393);

    auto ta2_xz_yyy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 394);

    auto ta2_xz_yyy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 395);

    auto ta2_xz_yyy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 396);

    auto ta2_xz_yyy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 397);

    auto ta2_xz_yyy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 398);

    auto ta2_xz_yyy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 399);

    auto ta2_xz_yyy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 400);

    auto ta2_xz_yyy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 401);

    auto ta2_xz_yyy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 402);

    auto ta2_xz_yyy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 403);

    auto ta2_xz_yyy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 404);

    auto ta2_xz_yyz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 406);

    auto ta2_xz_yyz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 407);

    auto ta2_xz_yyz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 408);

    auto ta2_xz_yyz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 409);

    auto ta2_xz_yyz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 410);

    auto ta2_xz_yyz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 411);

    auto ta2_xz_yyz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 412);

    auto ta2_xz_yyz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 413);

    auto ta2_xz_yyz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 414);

    auto ta2_xz_yyz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 415);

    auto ta2_xz_yyz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 416);

    auto ta2_xz_yyz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 417);

    auto ta2_xz_yyz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 418);

    auto ta2_xz_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 419);

    auto ta2_xz_yzz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 420);

    auto ta2_xz_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 422);

    auto ta2_xz_yzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 424);

    auto ta2_xz_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 425);

    auto ta2_xz_yzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 427);

    auto ta2_xz_yzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 428);

    auto ta2_xz_yzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 429);

    auto ta2_xz_yzz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 430);

    auto ta2_xz_yzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 431);

    auto ta2_xz_yzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 432);

    auto ta2_xz_yzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 433);

    auto ta2_xz_yzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 434);

    auto ta2_xz_zzz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 435);

    auto ta2_xz_zzz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 436);

    auto ta2_xz_zzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 437);

    auto ta2_xz_zzz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 438);

    auto ta2_xz_zzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 439);

    auto ta2_xz_zzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 440);

    auto ta2_xz_zzz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 441);

    auto ta2_xz_zzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 442);

    auto ta2_xz_zzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 443);

    auto ta2_xz_zzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 444);

    auto ta2_xz_zzz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 445);

    auto ta2_xz_zzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 446);

    auto ta2_xz_zzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 447);

    auto ta2_xz_zzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 448);

    auto ta2_xz_zzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 449);

    auto ta2_yy_xxx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 450);

    auto ta2_yy_xxx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 451);

    auto ta2_yy_xxx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 452);

    auto ta2_yy_xxx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 453);

    auto ta2_yy_xxx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 454);

    auto ta2_yy_xxx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 455);

    auto ta2_yy_xxx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 456);

    auto ta2_yy_xxx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 457);

    auto ta2_yy_xxx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 458);

    auto ta2_yy_xxx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 459);

    auto ta2_yy_xxx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 460);

    auto ta2_yy_xxx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 461);

    auto ta2_yy_xxx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 462);

    auto ta2_yy_xxx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 463);

    auto ta2_yy_xxx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 464);

    auto ta2_yy_xxy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 465);

    auto ta2_yy_xxy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 466);

    auto ta2_yy_xxy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 467);

    auto ta2_yy_xxy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 468);

    auto ta2_yy_xxy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 469);

    auto ta2_yy_xxy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 470);

    auto ta2_yy_xxy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 471);

    auto ta2_yy_xxy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 472);

    auto ta2_yy_xxy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 473);

    auto ta2_yy_xxy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 474);

    auto ta2_yy_xxy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 475);

    auto ta2_yy_xxy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 476);

    auto ta2_yy_xxy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 477);

    auto ta2_yy_xxy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 478);

    auto ta2_yy_xxz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 480);

    auto ta2_yy_xxz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 481);

    auto ta2_yy_xxz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 482);

    auto ta2_yy_xxz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 483);

    auto ta2_yy_xxz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 485);

    auto ta2_yy_xxz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 486);

    auto ta2_yy_xxz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 489);

    auto ta2_yy_xxz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 491);

    auto ta2_yy_xxz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 492);

    auto ta2_yy_xxz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 493);

    auto ta2_yy_xxz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 494);

    auto ta2_yy_xyy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 495);

    auto ta2_yy_xyy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 496);

    auto ta2_yy_xyy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 498);

    auto ta2_yy_xyy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 499);

    auto ta2_yy_xyy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 501);

    auto ta2_yy_xyy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 502);

    auto ta2_yy_xyy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 503);

    auto ta2_yy_xyy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 505);

    auto ta2_yy_xyy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 506);

    auto ta2_yy_xyy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 507);

    auto ta2_yy_xyy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 508);

    auto ta2_yy_xyy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 509);

    auto ta2_yy_xyz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 521);

    auto ta2_yy_xyz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 522);

    auto ta2_yy_xyz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 523);

    auto ta2_yy_xzz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 525);

    auto ta2_yy_xzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 527);

    auto ta2_yy_xzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 529);

    auto ta2_yy_xzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 530);

    auto ta2_yy_xzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 532);

    auto ta2_yy_xzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 533);

    auto ta2_yy_xzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 534);

    auto ta2_yy_xzz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 535);

    auto ta2_yy_xzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 536);

    auto ta2_yy_xzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 537);

    auto ta2_yy_xzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 538);

    auto ta2_yy_xzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 539);

    auto ta2_yy_yyy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 540);

    auto ta2_yy_yyy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 541);

    auto ta2_yy_yyy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 542);

    auto ta2_yy_yyy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 543);

    auto ta2_yy_yyy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 544);

    auto ta2_yy_yyy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 545);

    auto ta2_yy_yyy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 546);

    auto ta2_yy_yyy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 547);

    auto ta2_yy_yyy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 548);

    auto ta2_yy_yyy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 549);

    auto ta2_yy_yyy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 550);

    auto ta2_yy_yyy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 551);

    auto ta2_yy_yyy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 552);

    auto ta2_yy_yyy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 553);

    auto ta2_yy_yyy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 554);

    auto ta2_yy_yyz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 555);

    auto ta2_yy_yyz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 556);

    auto ta2_yy_yyz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 557);

    auto ta2_yy_yyz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 558);

    auto ta2_yy_yyz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 559);

    auto ta2_yy_yyz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 560);

    auto ta2_yy_yyz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 561);

    auto ta2_yy_yyz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 562);

    auto ta2_yy_yyz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 563);

    auto ta2_yy_yyz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 564);

    auto ta2_yy_yyz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 565);

    auto ta2_yy_yyz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 566);

    auto ta2_yy_yyz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 567);

    auto ta2_yy_yyz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 568);

    auto ta2_yy_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 569);

    auto ta2_yy_yzz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 571);

    auto ta2_yy_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 572);

    auto ta2_yy_yzz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 573);

    auto ta2_yy_yzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 574);

    auto ta2_yy_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 575);

    auto ta2_yy_yzz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 576);

    auto ta2_yy_yzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 577);

    auto ta2_yy_yzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 578);

    auto ta2_yy_yzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 579);

    auto ta2_yy_yzz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 580);

    auto ta2_yy_yzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 581);

    auto ta2_yy_yzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 582);

    auto ta2_yy_yzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 583);

    auto ta2_yy_yzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 584);

    auto ta2_yy_zzz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 585);

    auto ta2_yy_zzz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 586);

    auto ta2_yy_zzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 587);

    auto ta2_yy_zzz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 588);

    auto ta2_yy_zzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 589);

    auto ta2_yy_zzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 590);

    auto ta2_yy_zzz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 591);

    auto ta2_yy_zzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 592);

    auto ta2_yy_zzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 593);

    auto ta2_yy_zzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 594);

    auto ta2_yy_zzz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 595);

    auto ta2_yy_zzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 596);

    auto ta2_yy_zzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 597);

    auto ta2_yy_zzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 598);

    auto ta2_yy_zzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 599);

    auto ta2_yz_xxx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 600);

    auto ta2_yz_xxx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 601);

    auto ta2_yz_xxx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 602);

    auto ta2_yz_xxx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 603);

    auto ta2_yz_xxx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 604);

    auto ta2_yz_xxx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 605);

    auto ta2_yz_xxx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 606);

    auto ta2_yz_xxx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 607);

    auto ta2_yz_xxx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 608);

    auto ta2_yz_xxx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 609);

    auto ta2_yz_xxx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 610);

    auto ta2_yz_xxx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 611);

    auto ta2_yz_xxx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 612);

    auto ta2_yz_xxx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 613);

    auto ta2_yz_xxx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 614);

    auto ta2_yz_xxy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 615);

    auto ta2_yz_xxy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 616);

    auto ta2_yz_xxy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 617);

    auto ta2_yz_xxy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 618);

    auto ta2_yz_xxy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 620);

    auto ta2_yz_xxy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 621);

    auto ta2_yz_xxy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 624);

    auto ta2_yz_xxy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 625);

    auto ta2_yz_xxy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 626);

    auto ta2_yz_xxy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 627);

    auto ta2_yz_xxy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 628);

    auto ta2_yz_xxz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 630);

    auto ta2_yz_xxz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 631);

    auto ta2_yz_xxz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 632);

    auto ta2_yz_xxz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 633);

    auto ta2_yz_xxz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 634);

    auto ta2_yz_xxz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 635);

    auto ta2_yz_xxz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 636);

    auto ta2_yz_xxz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 637);

    auto ta2_yz_xxz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 638);

    auto ta2_yz_xxz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 639);

    auto ta2_yz_xxz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 641);

    auto ta2_yz_xxz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 642);

    auto ta2_yz_xxz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 643);

    auto ta2_yz_xxz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 644);

    auto ta2_yz_xyy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 645);

    auto ta2_yz_xyy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 646);

    auto ta2_yz_xyy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 648);

    auto ta2_yz_xyy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 649);

    auto ta2_yz_xyy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 651);

    auto ta2_yz_xyy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 652);

    auto ta2_yz_xyy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 653);

    auto ta2_yz_xyy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 655);

    auto ta2_yz_xyy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 656);

    auto ta2_yz_xyy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 657);

    auto ta2_yz_xyy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 658);

    auto ta2_yz_xyy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 659);

    auto ta2_yz_xyz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 671);

    auto ta2_yz_xyz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 672);

    auto ta2_yz_xyz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 673);

    auto ta2_yz_xzz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 675);

    auto ta2_yz_xzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 677);

    auto ta2_yz_xzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 679);

    auto ta2_yz_xzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 680);

    auto ta2_yz_xzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 682);

    auto ta2_yz_xzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 683);

    auto ta2_yz_xzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 684);

    auto ta2_yz_xzz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 685);

    auto ta2_yz_xzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 686);

    auto ta2_yz_xzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 687);

    auto ta2_yz_xzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 688);

    auto ta2_yz_xzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 689);

    auto ta2_yz_yyy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 690);

    auto ta2_yz_yyy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 691);

    auto ta2_yz_yyy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 692);

    auto ta2_yz_yyy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 693);

    auto ta2_yz_yyy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 694);

    auto ta2_yz_yyy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 695);

    auto ta2_yz_yyy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 696);

    auto ta2_yz_yyy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 697);

    auto ta2_yz_yyy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 698);

    auto ta2_yz_yyy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 699);

    auto ta2_yz_yyy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 700);

    auto ta2_yz_yyy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 701);

    auto ta2_yz_yyy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 702);

    auto ta2_yz_yyy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 703);

    auto ta2_yz_yyy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 704);

    auto ta2_yz_yyz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 706);

    auto ta2_yz_yyz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 707);

    auto ta2_yz_yyz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 708);

    auto ta2_yz_yyz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 709);

    auto ta2_yz_yyz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 710);

    auto ta2_yz_yyz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 711);

    auto ta2_yz_yyz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 712);

    auto ta2_yz_yyz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 713);

    auto ta2_yz_yyz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 714);

    auto ta2_yz_yyz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 715);

    auto ta2_yz_yyz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 716);

    auto ta2_yz_yyz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 717);

    auto ta2_yz_yyz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 718);

    auto ta2_yz_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 719);

    auto ta2_yz_yzz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 720);

    auto ta2_yz_yzz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 721);

    auto ta2_yz_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 722);

    auto ta2_yz_yzz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 723);

    auto ta2_yz_yzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 724);

    auto ta2_yz_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 725);

    auto ta2_yz_yzz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 726);

    auto ta2_yz_yzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 727);

    auto ta2_yz_yzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 728);

    auto ta2_yz_yzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 729);

    auto ta2_yz_yzz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 730);

    auto ta2_yz_yzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 731);

    auto ta2_yz_yzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 732);

    auto ta2_yz_yzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 733);

    auto ta2_yz_yzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 734);

    auto ta2_yz_zzz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 735);

    auto ta2_yz_zzz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 736);

    auto ta2_yz_zzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 737);

    auto ta2_yz_zzz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 738);

    auto ta2_yz_zzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 739);

    auto ta2_yz_zzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 740);

    auto ta2_yz_zzz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 741);

    auto ta2_yz_zzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 742);

    auto ta2_yz_zzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 743);

    auto ta2_yz_zzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 744);

    auto ta2_yz_zzz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 745);

    auto ta2_yz_zzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 746);

    auto ta2_yz_zzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 747);

    auto ta2_yz_zzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 748);

    auto ta2_yz_zzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 749);

    auto ta2_zz_xxx_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 750);

    auto ta2_zz_xxx_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 751);

    auto ta2_zz_xxx_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 752);

    auto ta2_zz_xxx_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 753);

    auto ta2_zz_xxx_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 754);

    auto ta2_zz_xxx_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 755);

    auto ta2_zz_xxx_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 756);

    auto ta2_zz_xxx_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 757);

    auto ta2_zz_xxx_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 758);

    auto ta2_zz_xxx_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 759);

    auto ta2_zz_xxx_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 760);

    auto ta2_zz_xxx_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 761);

    auto ta2_zz_xxx_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 762);

    auto ta2_zz_xxx_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 763);

    auto ta2_zz_xxx_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 764);

    auto ta2_zz_xxy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 765);

    auto ta2_zz_xxy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 766);

    auto ta2_zz_xxy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 767);

    auto ta2_zz_xxy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 768);

    auto ta2_zz_xxy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 770);

    auto ta2_zz_xxy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 771);

    auto ta2_zz_xxy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 774);

    auto ta2_zz_xxy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 775);

    auto ta2_zz_xxy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 776);

    auto ta2_zz_xxy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 777);

    auto ta2_zz_xxy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 778);

    auto ta2_zz_xxz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 780);

    auto ta2_zz_xxz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 781);

    auto ta2_zz_xxz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 782);

    auto ta2_zz_xxz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 783);

    auto ta2_zz_xxz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 784);

    auto ta2_zz_xxz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 785);

    auto ta2_zz_xxz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 786);

    auto ta2_zz_xxz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 787);

    auto ta2_zz_xxz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 788);

    auto ta2_zz_xxz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 789);

    auto ta2_zz_xxz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 791);

    auto ta2_zz_xxz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 792);

    auto ta2_zz_xxz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 793);

    auto ta2_zz_xxz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 794);

    auto ta2_zz_xyy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 795);

    auto ta2_zz_xyy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 796);

    auto ta2_zz_xyy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 798);

    auto ta2_zz_xyy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 799);

    auto ta2_zz_xyy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 801);

    auto ta2_zz_xyy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 802);

    auto ta2_zz_xyy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 803);

    auto ta2_zz_xyy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 805);

    auto ta2_zz_xyy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 806);

    auto ta2_zz_xyy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 807);

    auto ta2_zz_xyy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 808);

    auto ta2_zz_xyy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 809);

    auto ta2_zz_xyz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 821);

    auto ta2_zz_xyz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 822);

    auto ta2_zz_xyz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 823);

    auto ta2_zz_xzz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 825);

    auto ta2_zz_xzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 827);

    auto ta2_zz_xzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 829);

    auto ta2_zz_xzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 830);

    auto ta2_zz_xzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 832);

    auto ta2_zz_xzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 833);

    auto ta2_zz_xzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 834);

    auto ta2_zz_xzz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 835);

    auto ta2_zz_xzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 836);

    auto ta2_zz_xzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 837);

    auto ta2_zz_xzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 838);

    auto ta2_zz_xzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 839);

    auto ta2_zz_yyy_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 840);

    auto ta2_zz_yyy_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 841);

    auto ta2_zz_yyy_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 842);

    auto ta2_zz_yyy_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 843);

    auto ta2_zz_yyy_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 844);

    auto ta2_zz_yyy_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 845);

    auto ta2_zz_yyy_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 846);

    auto ta2_zz_yyy_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 847);

    auto ta2_zz_yyy_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 848);

    auto ta2_zz_yyy_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 849);

    auto ta2_zz_yyy_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 850);

    auto ta2_zz_yyy_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 851);

    auto ta2_zz_yyy_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 852);

    auto ta2_zz_yyy_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 853);

    auto ta2_zz_yyy_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 854);

    auto ta2_zz_yyz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 856);

    auto ta2_zz_yyz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 857);

    auto ta2_zz_yyz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 858);

    auto ta2_zz_yyz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 859);

    auto ta2_zz_yyz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 860);

    auto ta2_zz_yyz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 861);

    auto ta2_zz_yyz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 862);

    auto ta2_zz_yyz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 863);

    auto ta2_zz_yyz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 864);

    auto ta2_zz_yyz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 865);

    auto ta2_zz_yyz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 866);

    auto ta2_zz_yyz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 867);

    auto ta2_zz_yyz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 868);

    auto ta2_zz_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 869);

    auto ta2_zz_yzz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 870);

    auto ta2_zz_yzz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 871);

    auto ta2_zz_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 872);

    auto ta2_zz_yzz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 873);

    auto ta2_zz_yzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 874);

    auto ta2_zz_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 875);

    auto ta2_zz_yzz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 876);

    auto ta2_zz_yzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 877);

    auto ta2_zz_yzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 878);

    auto ta2_zz_yzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 879);

    auto ta2_zz_yzz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 880);

    auto ta2_zz_yzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 881);

    auto ta2_zz_yzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 882);

    auto ta2_zz_yzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 883);

    auto ta2_zz_yzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 884);

    auto ta2_zz_zzz_xxxx_1 = pbuffer.data(idx_npot_geom_020_1_fg + 885);

    auto ta2_zz_zzz_xxxy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 886);

    auto ta2_zz_zzz_xxxz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 887);

    auto ta2_zz_zzz_xxyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 888);

    auto ta2_zz_zzz_xxyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 889);

    auto ta2_zz_zzz_xxzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 890);

    auto ta2_zz_zzz_xyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 891);

    auto ta2_zz_zzz_xyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 892);

    auto ta2_zz_zzz_xyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 893);

    auto ta2_zz_zzz_xzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 894);

    auto ta2_zz_zzz_yyyy_1 = pbuffer.data(idx_npot_geom_020_1_fg + 895);

    auto ta2_zz_zzz_yyyz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 896);

    auto ta2_zz_zzz_yyzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 897);

    auto ta2_zz_zzz_yzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 898);

    auto ta2_zz_zzz_zzzz_1 = pbuffer.data(idx_npot_geom_020_1_fg + 899);

    // Set up 0-15 components of targeted buffer : GG

    auto ta2_xx_xxxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg);

    auto ta2_xx_xxxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1);

    auto ta2_xx_xxxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 2);

    auto ta2_xx_xxxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 3);

    auto ta2_xx_xxxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 4);

    auto ta2_xx_xxxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 5);

    auto ta2_xx_xxxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 6);

    auto ta2_xx_xxxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 7);

    auto ta2_xx_xxxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 8);

    auto ta2_xx_xxxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 9);

    auto ta2_xx_xxxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 10);

    auto ta2_xx_xxxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 11);

    auto ta2_xx_xxxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 12);

    auto ta2_xx_xxxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 13);

    auto ta2_xx_xxxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 14);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_x_xxx_xxxx_1,   \
                             ta1_x_xxx_xxxy_1,   \
                             ta1_x_xxx_xxxz_1,   \
                             ta1_x_xxx_xxyy_1,   \
                             ta1_x_xxx_xxyz_1,   \
                             ta1_x_xxx_xxzz_1,   \
                             ta1_x_xxx_xyyy_1,   \
                             ta1_x_xxx_xyyz_1,   \
                             ta1_x_xxx_xyzz_1,   \
                             ta1_x_xxx_xzzz_1,   \
                             ta1_x_xxx_yyyy_1,   \
                             ta1_x_xxx_yyyz_1,   \
                             ta1_x_xxx_yyzz_1,   \
                             ta1_x_xxx_yzzz_1,   \
                             ta1_x_xxx_zzzz_1,   \
                             ta2_xx_xx_xxxx_0,   \
                             ta2_xx_xx_xxxx_1,   \
                             ta2_xx_xx_xxxy_0,   \
                             ta2_xx_xx_xxxy_1,   \
                             ta2_xx_xx_xxxz_0,   \
                             ta2_xx_xx_xxxz_1,   \
                             ta2_xx_xx_xxyy_0,   \
                             ta2_xx_xx_xxyy_1,   \
                             ta2_xx_xx_xxyz_0,   \
                             ta2_xx_xx_xxyz_1,   \
                             ta2_xx_xx_xxzz_0,   \
                             ta2_xx_xx_xxzz_1,   \
                             ta2_xx_xx_xyyy_0,   \
                             ta2_xx_xx_xyyy_1,   \
                             ta2_xx_xx_xyyz_0,   \
                             ta2_xx_xx_xyyz_1,   \
                             ta2_xx_xx_xyzz_0,   \
                             ta2_xx_xx_xyzz_1,   \
                             ta2_xx_xx_xzzz_0,   \
                             ta2_xx_xx_xzzz_1,   \
                             ta2_xx_xx_yyyy_0,   \
                             ta2_xx_xx_yyyy_1,   \
                             ta2_xx_xx_yyyz_0,   \
                             ta2_xx_xx_yyyz_1,   \
                             ta2_xx_xx_yyzz_0,   \
                             ta2_xx_xx_yyzz_1,   \
                             ta2_xx_xx_yzzz_0,   \
                             ta2_xx_xx_yzzz_1,   \
                             ta2_xx_xx_zzzz_0,   \
                             ta2_xx_xx_zzzz_1,   \
                             ta2_xx_xxx_xxx_0,   \
                             ta2_xx_xxx_xxx_1,   \
                             ta2_xx_xxx_xxxx_0,  \
                             ta2_xx_xxx_xxxx_1,  \
                             ta2_xx_xxx_xxxy_0,  \
                             ta2_xx_xxx_xxxy_1,  \
                             ta2_xx_xxx_xxxz_0,  \
                             ta2_xx_xxx_xxxz_1,  \
                             ta2_xx_xxx_xxy_0,   \
                             ta2_xx_xxx_xxy_1,   \
                             ta2_xx_xxx_xxyy_0,  \
                             ta2_xx_xxx_xxyy_1,  \
                             ta2_xx_xxx_xxyz_0,  \
                             ta2_xx_xxx_xxyz_1,  \
                             ta2_xx_xxx_xxz_0,   \
                             ta2_xx_xxx_xxz_1,   \
                             ta2_xx_xxx_xxzz_0,  \
                             ta2_xx_xxx_xxzz_1,  \
                             ta2_xx_xxx_xyy_0,   \
                             ta2_xx_xxx_xyy_1,   \
                             ta2_xx_xxx_xyyy_0,  \
                             ta2_xx_xxx_xyyy_1,  \
                             ta2_xx_xxx_xyyz_0,  \
                             ta2_xx_xxx_xyyz_1,  \
                             ta2_xx_xxx_xyz_0,   \
                             ta2_xx_xxx_xyz_1,   \
                             ta2_xx_xxx_xyzz_0,  \
                             ta2_xx_xxx_xyzz_1,  \
                             ta2_xx_xxx_xzz_0,   \
                             ta2_xx_xxx_xzz_1,   \
                             ta2_xx_xxx_xzzz_0,  \
                             ta2_xx_xxx_xzzz_1,  \
                             ta2_xx_xxx_yyy_0,   \
                             ta2_xx_xxx_yyy_1,   \
                             ta2_xx_xxx_yyyy_0,  \
                             ta2_xx_xxx_yyyy_1,  \
                             ta2_xx_xxx_yyyz_0,  \
                             ta2_xx_xxx_yyyz_1,  \
                             ta2_xx_xxx_yyz_0,   \
                             ta2_xx_xxx_yyz_1,   \
                             ta2_xx_xxx_yyzz_0,  \
                             ta2_xx_xxx_yyzz_1,  \
                             ta2_xx_xxx_yzz_0,   \
                             ta2_xx_xxx_yzz_1,   \
                             ta2_xx_xxx_yzzz_0,  \
                             ta2_xx_xxx_yzzz_1,  \
                             ta2_xx_xxx_zzz_0,   \
                             ta2_xx_xxx_zzz_1,   \
                             ta2_xx_xxx_zzzz_0,  \
                             ta2_xx_xxx_zzzz_1,  \
                             ta2_xx_xxxx_xxxx_0, \
                             ta2_xx_xxxx_xxxy_0, \
                             ta2_xx_xxxx_xxxz_0, \
                             ta2_xx_xxxx_xxyy_0, \
                             ta2_xx_xxxx_xxyz_0, \
                             ta2_xx_xxxx_xxzz_0, \
                             ta2_xx_xxxx_xyyy_0, \
                             ta2_xx_xxxx_xyyz_0, \
                             ta2_xx_xxxx_xyzz_0, \
                             ta2_xx_xxxx_xzzz_0, \
                             ta2_xx_xxxx_yyyy_0, \
                             ta2_xx_xxxx_yyyz_0, \
                             ta2_xx_xxxx_yyzz_0, \
                             ta2_xx_xxxx_yzzz_0, \
                             ta2_xx_xxxx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxxx_xxxx_0[i] = 3.0 * ta2_xx_xx_xxxx_0[i] * fe_0 - 3.0 * ta2_xx_xx_xxxx_1[i] * fe_0 + 4.0 * ta2_xx_xxx_xxx_0[i] * fe_0 -
                                4.0 * ta2_xx_xxx_xxx_1[i] * fe_0 + 2.0 * ta1_x_xxx_xxxx_1[i] + ta2_xx_xxx_xxxx_0[i] * pa_x[i] -
                                ta2_xx_xxx_xxxx_1[i] * pc_x[i];

        ta2_xx_xxxx_xxxy_0[i] = 3.0 * ta2_xx_xx_xxxy_0[i] * fe_0 - 3.0 * ta2_xx_xx_xxxy_1[i] * fe_0 + 3.0 * ta2_xx_xxx_xxy_0[i] * fe_0 -
                                3.0 * ta2_xx_xxx_xxy_1[i] * fe_0 + 2.0 * ta1_x_xxx_xxxy_1[i] + ta2_xx_xxx_xxxy_0[i] * pa_x[i] -
                                ta2_xx_xxx_xxxy_1[i] * pc_x[i];

        ta2_xx_xxxx_xxxz_0[i] = 3.0 * ta2_xx_xx_xxxz_0[i] * fe_0 - 3.0 * ta2_xx_xx_xxxz_1[i] * fe_0 + 3.0 * ta2_xx_xxx_xxz_0[i] * fe_0 -
                                3.0 * ta2_xx_xxx_xxz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xxxz_1[i] + ta2_xx_xxx_xxxz_0[i] * pa_x[i] -
                                ta2_xx_xxx_xxxz_1[i] * pc_x[i];

        ta2_xx_xxxx_xxyy_0[i] = 3.0 * ta2_xx_xx_xxyy_0[i] * fe_0 - 3.0 * ta2_xx_xx_xxyy_1[i] * fe_0 + 2.0 * ta2_xx_xxx_xyy_0[i] * fe_0 -
                                2.0 * ta2_xx_xxx_xyy_1[i] * fe_0 + 2.0 * ta1_x_xxx_xxyy_1[i] + ta2_xx_xxx_xxyy_0[i] * pa_x[i] -
                                ta2_xx_xxx_xxyy_1[i] * pc_x[i];

        ta2_xx_xxxx_xxyz_0[i] = 3.0 * ta2_xx_xx_xxyz_0[i] * fe_0 - 3.0 * ta2_xx_xx_xxyz_1[i] * fe_0 + 2.0 * ta2_xx_xxx_xyz_0[i] * fe_0 -
                                2.0 * ta2_xx_xxx_xyz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xxyz_1[i] + ta2_xx_xxx_xxyz_0[i] * pa_x[i] -
                                ta2_xx_xxx_xxyz_1[i] * pc_x[i];

        ta2_xx_xxxx_xxzz_0[i] = 3.0 * ta2_xx_xx_xxzz_0[i] * fe_0 - 3.0 * ta2_xx_xx_xxzz_1[i] * fe_0 + 2.0 * ta2_xx_xxx_xzz_0[i] * fe_0 -
                                2.0 * ta2_xx_xxx_xzz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xxzz_1[i] + ta2_xx_xxx_xxzz_0[i] * pa_x[i] -
                                ta2_xx_xxx_xxzz_1[i] * pc_x[i];

        ta2_xx_xxxx_xyyy_0[i] = 3.0 * ta2_xx_xx_xyyy_0[i] * fe_0 - 3.0 * ta2_xx_xx_xyyy_1[i] * fe_0 + ta2_xx_xxx_yyy_0[i] * fe_0 -
                                ta2_xx_xxx_yyy_1[i] * fe_0 + 2.0 * ta1_x_xxx_xyyy_1[i] + ta2_xx_xxx_xyyy_0[i] * pa_x[i] -
                                ta2_xx_xxx_xyyy_1[i] * pc_x[i];

        ta2_xx_xxxx_xyyz_0[i] = 3.0 * ta2_xx_xx_xyyz_0[i] * fe_0 - 3.0 * ta2_xx_xx_xyyz_1[i] * fe_0 + ta2_xx_xxx_yyz_0[i] * fe_0 -
                                ta2_xx_xxx_yyz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xyyz_1[i] + ta2_xx_xxx_xyyz_0[i] * pa_x[i] -
                                ta2_xx_xxx_xyyz_1[i] * pc_x[i];

        ta2_xx_xxxx_xyzz_0[i] = 3.0 * ta2_xx_xx_xyzz_0[i] * fe_0 - 3.0 * ta2_xx_xx_xyzz_1[i] * fe_0 + ta2_xx_xxx_yzz_0[i] * fe_0 -
                                ta2_xx_xxx_yzz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xyzz_1[i] + ta2_xx_xxx_xyzz_0[i] * pa_x[i] -
                                ta2_xx_xxx_xyzz_1[i] * pc_x[i];

        ta2_xx_xxxx_xzzz_0[i] = 3.0 * ta2_xx_xx_xzzz_0[i] * fe_0 - 3.0 * ta2_xx_xx_xzzz_1[i] * fe_0 + ta2_xx_xxx_zzz_0[i] * fe_0 -
                                ta2_xx_xxx_zzz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xzzz_1[i] + ta2_xx_xxx_xzzz_0[i] * pa_x[i] -
                                ta2_xx_xxx_xzzz_1[i] * pc_x[i];

        ta2_xx_xxxx_yyyy_0[i] = 3.0 * ta2_xx_xx_yyyy_0[i] * fe_0 - 3.0 * ta2_xx_xx_yyyy_1[i] * fe_0 + 2.0 * ta1_x_xxx_yyyy_1[i] +
                                ta2_xx_xxx_yyyy_0[i] * pa_x[i] - ta2_xx_xxx_yyyy_1[i] * pc_x[i];

        ta2_xx_xxxx_yyyz_0[i] = 3.0 * ta2_xx_xx_yyyz_0[i] * fe_0 - 3.0 * ta2_xx_xx_yyyz_1[i] * fe_0 + 2.0 * ta1_x_xxx_yyyz_1[i] +
                                ta2_xx_xxx_yyyz_0[i] * pa_x[i] - ta2_xx_xxx_yyyz_1[i] * pc_x[i];

        ta2_xx_xxxx_yyzz_0[i] = 3.0 * ta2_xx_xx_yyzz_0[i] * fe_0 - 3.0 * ta2_xx_xx_yyzz_1[i] * fe_0 + 2.0 * ta1_x_xxx_yyzz_1[i] +
                                ta2_xx_xxx_yyzz_0[i] * pa_x[i] - ta2_xx_xxx_yyzz_1[i] * pc_x[i];

        ta2_xx_xxxx_yzzz_0[i] = 3.0 * ta2_xx_xx_yzzz_0[i] * fe_0 - 3.0 * ta2_xx_xx_yzzz_1[i] * fe_0 + 2.0 * ta1_x_xxx_yzzz_1[i] +
                                ta2_xx_xxx_yzzz_0[i] * pa_x[i] - ta2_xx_xxx_yzzz_1[i] * pc_x[i];

        ta2_xx_xxxx_zzzz_0[i] = 3.0 * ta2_xx_xx_zzzz_0[i] * fe_0 - 3.0 * ta2_xx_xx_zzzz_1[i] * fe_0 + 2.0 * ta1_x_xxx_zzzz_1[i] +
                                ta2_xx_xxx_zzzz_0[i] * pa_x[i] - ta2_xx_xxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 15-30 components of targeted buffer : GG

    auto ta2_xx_xxxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 15);

    auto ta2_xx_xxxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 16);

    auto ta2_xx_xxxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 17);

    auto ta2_xx_xxxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 18);

    auto ta2_xx_xxxy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 19);

    auto ta2_xx_xxxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 20);

    auto ta2_xx_xxxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 21);

    auto ta2_xx_xxxy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 22);

    auto ta2_xx_xxxy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 23);

    auto ta2_xx_xxxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 24);

    auto ta2_xx_xxxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 25);

    auto ta2_xx_xxxy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 26);

    auto ta2_xx_xxxy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 27);

    auto ta2_xx_xxxy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 28);

    auto ta2_xx_xxxy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 29);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta2_xx_xxx_xxx_0,   \
                             ta2_xx_xxx_xxx_1,   \
                             ta2_xx_xxx_xxxx_0,  \
                             ta2_xx_xxx_xxxx_1,  \
                             ta2_xx_xxx_xxxy_0,  \
                             ta2_xx_xxx_xxxy_1,  \
                             ta2_xx_xxx_xxxz_0,  \
                             ta2_xx_xxx_xxxz_1,  \
                             ta2_xx_xxx_xxy_0,   \
                             ta2_xx_xxx_xxy_1,   \
                             ta2_xx_xxx_xxyy_0,  \
                             ta2_xx_xxx_xxyy_1,  \
                             ta2_xx_xxx_xxyz_0,  \
                             ta2_xx_xxx_xxyz_1,  \
                             ta2_xx_xxx_xxz_0,   \
                             ta2_xx_xxx_xxz_1,   \
                             ta2_xx_xxx_xxzz_0,  \
                             ta2_xx_xxx_xxzz_1,  \
                             ta2_xx_xxx_xyy_0,   \
                             ta2_xx_xxx_xyy_1,   \
                             ta2_xx_xxx_xyyy_0,  \
                             ta2_xx_xxx_xyyy_1,  \
                             ta2_xx_xxx_xyyz_0,  \
                             ta2_xx_xxx_xyyz_1,  \
                             ta2_xx_xxx_xyz_0,   \
                             ta2_xx_xxx_xyz_1,   \
                             ta2_xx_xxx_xyzz_0,  \
                             ta2_xx_xxx_xyzz_1,  \
                             ta2_xx_xxx_xzz_0,   \
                             ta2_xx_xxx_xzz_1,   \
                             ta2_xx_xxx_xzzz_0,  \
                             ta2_xx_xxx_xzzz_1,  \
                             ta2_xx_xxx_yyy_0,   \
                             ta2_xx_xxx_yyy_1,   \
                             ta2_xx_xxx_yyyy_0,  \
                             ta2_xx_xxx_yyyy_1,  \
                             ta2_xx_xxx_yyyz_0,  \
                             ta2_xx_xxx_yyyz_1,  \
                             ta2_xx_xxx_yyz_0,   \
                             ta2_xx_xxx_yyz_1,   \
                             ta2_xx_xxx_yyzz_0,  \
                             ta2_xx_xxx_yyzz_1,  \
                             ta2_xx_xxx_yzz_0,   \
                             ta2_xx_xxx_yzz_1,   \
                             ta2_xx_xxx_yzzz_0,  \
                             ta2_xx_xxx_yzzz_1,  \
                             ta2_xx_xxx_zzz_0,   \
                             ta2_xx_xxx_zzz_1,   \
                             ta2_xx_xxx_zzzz_0,  \
                             ta2_xx_xxx_zzzz_1,  \
                             ta2_xx_xxxy_xxxx_0, \
                             ta2_xx_xxxy_xxxy_0, \
                             ta2_xx_xxxy_xxxz_0, \
                             ta2_xx_xxxy_xxyy_0, \
                             ta2_xx_xxxy_xxyz_0, \
                             ta2_xx_xxxy_xxzz_0, \
                             ta2_xx_xxxy_xyyy_0, \
                             ta2_xx_xxxy_xyyz_0, \
                             ta2_xx_xxxy_xyzz_0, \
                             ta2_xx_xxxy_xzzz_0, \
                             ta2_xx_xxxy_yyyy_0, \
                             ta2_xx_xxxy_yyyz_0, \
                             ta2_xx_xxxy_yyzz_0, \
                             ta2_xx_xxxy_yzzz_0, \
                             ta2_xx_xxxy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxxy_xxxx_0[i] = ta2_xx_xxx_xxxx_0[i] * pa_y[i] - ta2_xx_xxx_xxxx_1[i] * pc_y[i];

        ta2_xx_xxxy_xxxy_0[i] =
            ta2_xx_xxx_xxx_0[i] * fe_0 - ta2_xx_xxx_xxx_1[i] * fe_0 + ta2_xx_xxx_xxxy_0[i] * pa_y[i] - ta2_xx_xxx_xxxy_1[i] * pc_y[i];

        ta2_xx_xxxy_xxxz_0[i] = ta2_xx_xxx_xxxz_0[i] * pa_y[i] - ta2_xx_xxx_xxxz_1[i] * pc_y[i];

        ta2_xx_xxxy_xxyy_0[i] =
            2.0 * ta2_xx_xxx_xxy_0[i] * fe_0 - 2.0 * ta2_xx_xxx_xxy_1[i] * fe_0 + ta2_xx_xxx_xxyy_0[i] * pa_y[i] - ta2_xx_xxx_xxyy_1[i] * pc_y[i];

        ta2_xx_xxxy_xxyz_0[i] =
            ta2_xx_xxx_xxz_0[i] * fe_0 - ta2_xx_xxx_xxz_1[i] * fe_0 + ta2_xx_xxx_xxyz_0[i] * pa_y[i] - ta2_xx_xxx_xxyz_1[i] * pc_y[i];

        ta2_xx_xxxy_xxzz_0[i] = ta2_xx_xxx_xxzz_0[i] * pa_y[i] - ta2_xx_xxx_xxzz_1[i] * pc_y[i];

        ta2_xx_xxxy_xyyy_0[i] =
            3.0 * ta2_xx_xxx_xyy_0[i] * fe_0 - 3.0 * ta2_xx_xxx_xyy_1[i] * fe_0 + ta2_xx_xxx_xyyy_0[i] * pa_y[i] - ta2_xx_xxx_xyyy_1[i] * pc_y[i];

        ta2_xx_xxxy_xyyz_0[i] =
            2.0 * ta2_xx_xxx_xyz_0[i] * fe_0 - 2.0 * ta2_xx_xxx_xyz_1[i] * fe_0 + ta2_xx_xxx_xyyz_0[i] * pa_y[i] - ta2_xx_xxx_xyyz_1[i] * pc_y[i];

        ta2_xx_xxxy_xyzz_0[i] =
            ta2_xx_xxx_xzz_0[i] * fe_0 - ta2_xx_xxx_xzz_1[i] * fe_0 + ta2_xx_xxx_xyzz_0[i] * pa_y[i] - ta2_xx_xxx_xyzz_1[i] * pc_y[i];

        ta2_xx_xxxy_xzzz_0[i] = ta2_xx_xxx_xzzz_0[i] * pa_y[i] - ta2_xx_xxx_xzzz_1[i] * pc_y[i];

        ta2_xx_xxxy_yyyy_0[i] =
            4.0 * ta2_xx_xxx_yyy_0[i] * fe_0 - 4.0 * ta2_xx_xxx_yyy_1[i] * fe_0 + ta2_xx_xxx_yyyy_0[i] * pa_y[i] - ta2_xx_xxx_yyyy_1[i] * pc_y[i];

        ta2_xx_xxxy_yyyz_0[i] =
            3.0 * ta2_xx_xxx_yyz_0[i] * fe_0 - 3.0 * ta2_xx_xxx_yyz_1[i] * fe_0 + ta2_xx_xxx_yyyz_0[i] * pa_y[i] - ta2_xx_xxx_yyyz_1[i] * pc_y[i];

        ta2_xx_xxxy_yyzz_0[i] =
            2.0 * ta2_xx_xxx_yzz_0[i] * fe_0 - 2.0 * ta2_xx_xxx_yzz_1[i] * fe_0 + ta2_xx_xxx_yyzz_0[i] * pa_y[i] - ta2_xx_xxx_yyzz_1[i] * pc_y[i];

        ta2_xx_xxxy_yzzz_0[i] =
            ta2_xx_xxx_zzz_0[i] * fe_0 - ta2_xx_xxx_zzz_1[i] * fe_0 + ta2_xx_xxx_yzzz_0[i] * pa_y[i] - ta2_xx_xxx_yzzz_1[i] * pc_y[i];

        ta2_xx_xxxy_zzzz_0[i] = ta2_xx_xxx_zzzz_0[i] * pa_y[i] - ta2_xx_xxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 30-45 components of targeted buffer : GG

    auto ta2_xx_xxxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 30);

    auto ta2_xx_xxxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 31);

    auto ta2_xx_xxxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 32);

    auto ta2_xx_xxxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 33);

    auto ta2_xx_xxxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 34);

    auto ta2_xx_xxxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 35);

    auto ta2_xx_xxxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 36);

    auto ta2_xx_xxxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 37);

    auto ta2_xx_xxxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 38);

    auto ta2_xx_xxxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 39);

    auto ta2_xx_xxxz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 40);

    auto ta2_xx_xxxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 41);

    auto ta2_xx_xxxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 42);

    auto ta2_xx_xxxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 43);

    auto ta2_xx_xxxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 44);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta2_xx_xxx_xxx_0,   \
                             ta2_xx_xxx_xxx_1,   \
                             ta2_xx_xxx_xxxx_0,  \
                             ta2_xx_xxx_xxxx_1,  \
                             ta2_xx_xxx_xxxy_0,  \
                             ta2_xx_xxx_xxxy_1,  \
                             ta2_xx_xxx_xxxz_0,  \
                             ta2_xx_xxx_xxxz_1,  \
                             ta2_xx_xxx_xxy_0,   \
                             ta2_xx_xxx_xxy_1,   \
                             ta2_xx_xxx_xxyy_0,  \
                             ta2_xx_xxx_xxyy_1,  \
                             ta2_xx_xxx_xxyz_0,  \
                             ta2_xx_xxx_xxyz_1,  \
                             ta2_xx_xxx_xxz_0,   \
                             ta2_xx_xxx_xxz_1,   \
                             ta2_xx_xxx_xxzz_0,  \
                             ta2_xx_xxx_xxzz_1,  \
                             ta2_xx_xxx_xyy_0,   \
                             ta2_xx_xxx_xyy_1,   \
                             ta2_xx_xxx_xyyy_0,  \
                             ta2_xx_xxx_xyyy_1,  \
                             ta2_xx_xxx_xyyz_0,  \
                             ta2_xx_xxx_xyyz_1,  \
                             ta2_xx_xxx_xyz_0,   \
                             ta2_xx_xxx_xyz_1,   \
                             ta2_xx_xxx_xyzz_0,  \
                             ta2_xx_xxx_xyzz_1,  \
                             ta2_xx_xxx_xzz_0,   \
                             ta2_xx_xxx_xzz_1,   \
                             ta2_xx_xxx_xzzz_0,  \
                             ta2_xx_xxx_xzzz_1,  \
                             ta2_xx_xxx_yyy_0,   \
                             ta2_xx_xxx_yyy_1,   \
                             ta2_xx_xxx_yyyy_0,  \
                             ta2_xx_xxx_yyyy_1,  \
                             ta2_xx_xxx_yyyz_0,  \
                             ta2_xx_xxx_yyyz_1,  \
                             ta2_xx_xxx_yyz_0,   \
                             ta2_xx_xxx_yyz_1,   \
                             ta2_xx_xxx_yyzz_0,  \
                             ta2_xx_xxx_yyzz_1,  \
                             ta2_xx_xxx_yzz_0,   \
                             ta2_xx_xxx_yzz_1,   \
                             ta2_xx_xxx_yzzz_0,  \
                             ta2_xx_xxx_yzzz_1,  \
                             ta2_xx_xxx_zzz_0,   \
                             ta2_xx_xxx_zzz_1,   \
                             ta2_xx_xxx_zzzz_0,  \
                             ta2_xx_xxx_zzzz_1,  \
                             ta2_xx_xxxz_xxxx_0, \
                             ta2_xx_xxxz_xxxy_0, \
                             ta2_xx_xxxz_xxxz_0, \
                             ta2_xx_xxxz_xxyy_0, \
                             ta2_xx_xxxz_xxyz_0, \
                             ta2_xx_xxxz_xxzz_0, \
                             ta2_xx_xxxz_xyyy_0, \
                             ta2_xx_xxxz_xyyz_0, \
                             ta2_xx_xxxz_xyzz_0, \
                             ta2_xx_xxxz_xzzz_0, \
                             ta2_xx_xxxz_yyyy_0, \
                             ta2_xx_xxxz_yyyz_0, \
                             ta2_xx_xxxz_yyzz_0, \
                             ta2_xx_xxxz_yzzz_0, \
                             ta2_xx_xxxz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxxz_xxxx_0[i] = ta2_xx_xxx_xxxx_0[i] * pa_z[i] - ta2_xx_xxx_xxxx_1[i] * pc_z[i];

        ta2_xx_xxxz_xxxy_0[i] = ta2_xx_xxx_xxxy_0[i] * pa_z[i] - ta2_xx_xxx_xxxy_1[i] * pc_z[i];

        ta2_xx_xxxz_xxxz_0[i] =
            ta2_xx_xxx_xxx_0[i] * fe_0 - ta2_xx_xxx_xxx_1[i] * fe_0 + ta2_xx_xxx_xxxz_0[i] * pa_z[i] - ta2_xx_xxx_xxxz_1[i] * pc_z[i];

        ta2_xx_xxxz_xxyy_0[i] = ta2_xx_xxx_xxyy_0[i] * pa_z[i] - ta2_xx_xxx_xxyy_1[i] * pc_z[i];

        ta2_xx_xxxz_xxyz_0[i] =
            ta2_xx_xxx_xxy_0[i] * fe_0 - ta2_xx_xxx_xxy_1[i] * fe_0 + ta2_xx_xxx_xxyz_0[i] * pa_z[i] - ta2_xx_xxx_xxyz_1[i] * pc_z[i];

        ta2_xx_xxxz_xxzz_0[i] =
            2.0 * ta2_xx_xxx_xxz_0[i] * fe_0 - 2.0 * ta2_xx_xxx_xxz_1[i] * fe_0 + ta2_xx_xxx_xxzz_0[i] * pa_z[i] - ta2_xx_xxx_xxzz_1[i] * pc_z[i];

        ta2_xx_xxxz_xyyy_0[i] = ta2_xx_xxx_xyyy_0[i] * pa_z[i] - ta2_xx_xxx_xyyy_1[i] * pc_z[i];

        ta2_xx_xxxz_xyyz_0[i] =
            ta2_xx_xxx_xyy_0[i] * fe_0 - ta2_xx_xxx_xyy_1[i] * fe_0 + ta2_xx_xxx_xyyz_0[i] * pa_z[i] - ta2_xx_xxx_xyyz_1[i] * pc_z[i];

        ta2_xx_xxxz_xyzz_0[i] =
            2.0 * ta2_xx_xxx_xyz_0[i] * fe_0 - 2.0 * ta2_xx_xxx_xyz_1[i] * fe_0 + ta2_xx_xxx_xyzz_0[i] * pa_z[i] - ta2_xx_xxx_xyzz_1[i] * pc_z[i];

        ta2_xx_xxxz_xzzz_0[i] =
            3.0 * ta2_xx_xxx_xzz_0[i] * fe_0 - 3.0 * ta2_xx_xxx_xzz_1[i] * fe_0 + ta2_xx_xxx_xzzz_0[i] * pa_z[i] - ta2_xx_xxx_xzzz_1[i] * pc_z[i];

        ta2_xx_xxxz_yyyy_0[i] = ta2_xx_xxx_yyyy_0[i] * pa_z[i] - ta2_xx_xxx_yyyy_1[i] * pc_z[i];

        ta2_xx_xxxz_yyyz_0[i] =
            ta2_xx_xxx_yyy_0[i] * fe_0 - ta2_xx_xxx_yyy_1[i] * fe_0 + ta2_xx_xxx_yyyz_0[i] * pa_z[i] - ta2_xx_xxx_yyyz_1[i] * pc_z[i];

        ta2_xx_xxxz_yyzz_0[i] =
            2.0 * ta2_xx_xxx_yyz_0[i] * fe_0 - 2.0 * ta2_xx_xxx_yyz_1[i] * fe_0 + ta2_xx_xxx_yyzz_0[i] * pa_z[i] - ta2_xx_xxx_yyzz_1[i] * pc_z[i];

        ta2_xx_xxxz_yzzz_0[i] =
            3.0 * ta2_xx_xxx_yzz_0[i] * fe_0 - 3.0 * ta2_xx_xxx_yzz_1[i] * fe_0 + ta2_xx_xxx_yzzz_0[i] * pa_z[i] - ta2_xx_xxx_yzzz_1[i] * pc_z[i];

        ta2_xx_xxxz_zzzz_0[i] =
            4.0 * ta2_xx_xxx_zzz_0[i] * fe_0 - 4.0 * ta2_xx_xxx_zzz_1[i] * fe_0 + ta2_xx_xxx_zzzz_0[i] * pa_z[i] - ta2_xx_xxx_zzzz_1[i] * pc_z[i];
    }

    // Set up 45-60 components of targeted buffer : GG

    auto ta2_xx_xxyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 45);

    auto ta2_xx_xxyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 46);

    auto ta2_xx_xxyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 47);

    auto ta2_xx_xxyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 48);

    auto ta2_xx_xxyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 49);

    auto ta2_xx_xxyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 50);

    auto ta2_xx_xxyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 51);

    auto ta2_xx_xxyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 52);

    auto ta2_xx_xxyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 53);

    auto ta2_xx_xxyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 54);

    auto ta2_xx_xxyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 55);

    auto ta2_xx_xxyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 56);

    auto ta2_xx_xxyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 57);

    auto ta2_xx_xxyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 58);

    auto ta2_xx_xxyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 59);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_xyy_yyyy_1,   \
                             ta1_x_xyy_yyyz_1,   \
                             ta1_x_xyy_yyzz_1,   \
                             ta1_x_xyy_yzzz_1,   \
                             ta2_xx_xx_xxxx_0,   \
                             ta2_xx_xx_xxxx_1,   \
                             ta2_xx_xx_xxxy_0,   \
                             ta2_xx_xx_xxxy_1,   \
                             ta2_xx_xx_xxxz_0,   \
                             ta2_xx_xx_xxxz_1,   \
                             ta2_xx_xx_xxyy_0,   \
                             ta2_xx_xx_xxyy_1,   \
                             ta2_xx_xx_xxyz_0,   \
                             ta2_xx_xx_xxyz_1,   \
                             ta2_xx_xx_xxzz_0,   \
                             ta2_xx_xx_xxzz_1,   \
                             ta2_xx_xx_xyyy_0,   \
                             ta2_xx_xx_xyyy_1,   \
                             ta2_xx_xx_xyyz_0,   \
                             ta2_xx_xx_xyyz_1,   \
                             ta2_xx_xx_xyzz_0,   \
                             ta2_xx_xx_xyzz_1,   \
                             ta2_xx_xx_xzzz_0,   \
                             ta2_xx_xx_xzzz_1,   \
                             ta2_xx_xx_zzzz_0,   \
                             ta2_xx_xx_zzzz_1,   \
                             ta2_xx_xxy_xxx_0,   \
                             ta2_xx_xxy_xxx_1,   \
                             ta2_xx_xxy_xxxx_0,  \
                             ta2_xx_xxy_xxxx_1,  \
                             ta2_xx_xxy_xxxy_0,  \
                             ta2_xx_xxy_xxxy_1,  \
                             ta2_xx_xxy_xxxz_0,  \
                             ta2_xx_xxy_xxxz_1,  \
                             ta2_xx_xxy_xxy_0,   \
                             ta2_xx_xxy_xxy_1,   \
                             ta2_xx_xxy_xxyy_0,  \
                             ta2_xx_xxy_xxyy_1,  \
                             ta2_xx_xxy_xxyz_0,  \
                             ta2_xx_xxy_xxyz_1,  \
                             ta2_xx_xxy_xxz_0,   \
                             ta2_xx_xxy_xxz_1,   \
                             ta2_xx_xxy_xxzz_0,  \
                             ta2_xx_xxy_xxzz_1,  \
                             ta2_xx_xxy_xyy_0,   \
                             ta2_xx_xxy_xyy_1,   \
                             ta2_xx_xxy_xyyy_0,  \
                             ta2_xx_xxy_xyyy_1,  \
                             ta2_xx_xxy_xyyz_0,  \
                             ta2_xx_xxy_xyyz_1,  \
                             ta2_xx_xxy_xyz_0,   \
                             ta2_xx_xxy_xyz_1,   \
                             ta2_xx_xxy_xyzz_0,  \
                             ta2_xx_xxy_xyzz_1,  \
                             ta2_xx_xxy_xzz_0,   \
                             ta2_xx_xxy_xzz_1,   \
                             ta2_xx_xxy_xzzz_0,  \
                             ta2_xx_xxy_xzzz_1,  \
                             ta2_xx_xxy_zzzz_0,  \
                             ta2_xx_xxy_zzzz_1,  \
                             ta2_xx_xxyy_xxxx_0, \
                             ta2_xx_xxyy_xxxy_0, \
                             ta2_xx_xxyy_xxxz_0, \
                             ta2_xx_xxyy_xxyy_0, \
                             ta2_xx_xxyy_xxyz_0, \
                             ta2_xx_xxyy_xxzz_0, \
                             ta2_xx_xxyy_xyyy_0, \
                             ta2_xx_xxyy_xyyz_0, \
                             ta2_xx_xxyy_xyzz_0, \
                             ta2_xx_xxyy_xzzz_0, \
                             ta2_xx_xxyy_yyyy_0, \
                             ta2_xx_xxyy_yyyz_0, \
                             ta2_xx_xxyy_yyzz_0, \
                             ta2_xx_xxyy_yzzz_0, \
                             ta2_xx_xxyy_zzzz_0, \
                             ta2_xx_xyy_yyyy_0,  \
                             ta2_xx_xyy_yyyy_1,  \
                             ta2_xx_xyy_yyyz_0,  \
                             ta2_xx_xyy_yyyz_1,  \
                             ta2_xx_xyy_yyzz_0,  \
                             ta2_xx_xyy_yyzz_1,  \
                             ta2_xx_xyy_yzzz_0,  \
                             ta2_xx_xyy_yzzz_1,  \
                             ta2_xx_yy_yyyy_0,   \
                             ta2_xx_yy_yyyy_1,   \
                             ta2_xx_yy_yyyz_0,   \
                             ta2_xx_yy_yyyz_1,   \
                             ta2_xx_yy_yyzz_0,   \
                             ta2_xx_yy_yyzz_1,   \
                             ta2_xx_yy_yzzz_0,   \
                             ta2_xx_yy_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxyy_xxxx_0[i] =
            ta2_xx_xx_xxxx_0[i] * fe_0 - ta2_xx_xx_xxxx_1[i] * fe_0 + ta2_xx_xxy_xxxx_0[i] * pa_y[i] - ta2_xx_xxy_xxxx_1[i] * pc_y[i];

        ta2_xx_xxyy_xxxy_0[i] = ta2_xx_xx_xxxy_0[i] * fe_0 - ta2_xx_xx_xxxy_1[i] * fe_0 + ta2_xx_xxy_xxx_0[i] * fe_0 - ta2_xx_xxy_xxx_1[i] * fe_0 +
                                ta2_xx_xxy_xxxy_0[i] * pa_y[i] - ta2_xx_xxy_xxxy_1[i] * pc_y[i];

        ta2_xx_xxyy_xxxz_0[i] =
            ta2_xx_xx_xxxz_0[i] * fe_0 - ta2_xx_xx_xxxz_1[i] * fe_0 + ta2_xx_xxy_xxxz_0[i] * pa_y[i] - ta2_xx_xxy_xxxz_1[i] * pc_y[i];

        ta2_xx_xxyy_xxyy_0[i] = ta2_xx_xx_xxyy_0[i] * fe_0 - ta2_xx_xx_xxyy_1[i] * fe_0 + 2.0 * ta2_xx_xxy_xxy_0[i] * fe_0 -
                                2.0 * ta2_xx_xxy_xxy_1[i] * fe_0 + ta2_xx_xxy_xxyy_0[i] * pa_y[i] - ta2_xx_xxy_xxyy_1[i] * pc_y[i];

        ta2_xx_xxyy_xxyz_0[i] = ta2_xx_xx_xxyz_0[i] * fe_0 - ta2_xx_xx_xxyz_1[i] * fe_0 + ta2_xx_xxy_xxz_0[i] * fe_0 - ta2_xx_xxy_xxz_1[i] * fe_0 +
                                ta2_xx_xxy_xxyz_0[i] * pa_y[i] - ta2_xx_xxy_xxyz_1[i] * pc_y[i];

        ta2_xx_xxyy_xxzz_0[i] =
            ta2_xx_xx_xxzz_0[i] * fe_0 - ta2_xx_xx_xxzz_1[i] * fe_0 + ta2_xx_xxy_xxzz_0[i] * pa_y[i] - ta2_xx_xxy_xxzz_1[i] * pc_y[i];

        ta2_xx_xxyy_xyyy_0[i] = ta2_xx_xx_xyyy_0[i] * fe_0 - ta2_xx_xx_xyyy_1[i] * fe_0 + 3.0 * ta2_xx_xxy_xyy_0[i] * fe_0 -
                                3.0 * ta2_xx_xxy_xyy_1[i] * fe_0 + ta2_xx_xxy_xyyy_0[i] * pa_y[i] - ta2_xx_xxy_xyyy_1[i] * pc_y[i];

        ta2_xx_xxyy_xyyz_0[i] = ta2_xx_xx_xyyz_0[i] * fe_0 - ta2_xx_xx_xyyz_1[i] * fe_0 + 2.0 * ta2_xx_xxy_xyz_0[i] * fe_0 -
                                2.0 * ta2_xx_xxy_xyz_1[i] * fe_0 + ta2_xx_xxy_xyyz_0[i] * pa_y[i] - ta2_xx_xxy_xyyz_1[i] * pc_y[i];

        ta2_xx_xxyy_xyzz_0[i] = ta2_xx_xx_xyzz_0[i] * fe_0 - ta2_xx_xx_xyzz_1[i] * fe_0 + ta2_xx_xxy_xzz_0[i] * fe_0 - ta2_xx_xxy_xzz_1[i] * fe_0 +
                                ta2_xx_xxy_xyzz_0[i] * pa_y[i] - ta2_xx_xxy_xyzz_1[i] * pc_y[i];

        ta2_xx_xxyy_xzzz_0[i] =
            ta2_xx_xx_xzzz_0[i] * fe_0 - ta2_xx_xx_xzzz_1[i] * fe_0 + ta2_xx_xxy_xzzz_0[i] * pa_y[i] - ta2_xx_xxy_xzzz_1[i] * pc_y[i];

        ta2_xx_xxyy_yyyy_0[i] = ta2_xx_yy_yyyy_0[i] * fe_0 - ta2_xx_yy_yyyy_1[i] * fe_0 + 2.0 * ta1_x_xyy_yyyy_1[i] + ta2_xx_xyy_yyyy_0[i] * pa_x[i] -
                                ta2_xx_xyy_yyyy_1[i] * pc_x[i];

        ta2_xx_xxyy_yyyz_0[i] = ta2_xx_yy_yyyz_0[i] * fe_0 - ta2_xx_yy_yyyz_1[i] * fe_0 + 2.0 * ta1_x_xyy_yyyz_1[i] + ta2_xx_xyy_yyyz_0[i] * pa_x[i] -
                                ta2_xx_xyy_yyyz_1[i] * pc_x[i];

        ta2_xx_xxyy_yyzz_0[i] = ta2_xx_yy_yyzz_0[i] * fe_0 - ta2_xx_yy_yyzz_1[i] * fe_0 + 2.0 * ta1_x_xyy_yyzz_1[i] + ta2_xx_xyy_yyzz_0[i] * pa_x[i] -
                                ta2_xx_xyy_yyzz_1[i] * pc_x[i];

        ta2_xx_xxyy_yzzz_0[i] = ta2_xx_yy_yzzz_0[i] * fe_0 - ta2_xx_yy_yzzz_1[i] * fe_0 + 2.0 * ta1_x_xyy_yzzz_1[i] + ta2_xx_xyy_yzzz_0[i] * pa_x[i] -
                                ta2_xx_xyy_yzzz_1[i] * pc_x[i];

        ta2_xx_xxyy_zzzz_0[i] =
            ta2_xx_xx_zzzz_0[i] * fe_0 - ta2_xx_xx_zzzz_1[i] * fe_0 + ta2_xx_xxy_zzzz_0[i] * pa_y[i] - ta2_xx_xxy_zzzz_1[i] * pc_y[i];
    }

    // Set up 60-75 components of targeted buffer : GG

    auto ta2_xx_xxyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 60);

    auto ta2_xx_xxyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 61);

    auto ta2_xx_xxyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 62);

    auto ta2_xx_xxyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 63);

    auto ta2_xx_xxyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 64);

    auto ta2_xx_xxyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 65);

    auto ta2_xx_xxyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 66);

    auto ta2_xx_xxyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 67);

    auto ta2_xx_xxyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 68);

    auto ta2_xx_xxyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 69);

    auto ta2_xx_xxyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 70);

    auto ta2_xx_xxyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 71);

    auto ta2_xx_xxyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 72);

    auto ta2_xx_xxyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 73);

    auto ta2_xx_xxyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 74);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta2_xx_xxy_xxxy_0,  \
                             ta2_xx_xxy_xxxy_1,  \
                             ta2_xx_xxy_xxyy_0,  \
                             ta2_xx_xxy_xxyy_1,  \
                             ta2_xx_xxy_xyyy_0,  \
                             ta2_xx_xxy_xyyy_1,  \
                             ta2_xx_xxy_yyyy_0,  \
                             ta2_xx_xxy_yyyy_1,  \
                             ta2_xx_xxyz_xxxx_0, \
                             ta2_xx_xxyz_xxxy_0, \
                             ta2_xx_xxyz_xxxz_0, \
                             ta2_xx_xxyz_xxyy_0, \
                             ta2_xx_xxyz_xxyz_0, \
                             ta2_xx_xxyz_xxzz_0, \
                             ta2_xx_xxyz_xyyy_0, \
                             ta2_xx_xxyz_xyyz_0, \
                             ta2_xx_xxyz_xyzz_0, \
                             ta2_xx_xxyz_xzzz_0, \
                             ta2_xx_xxyz_yyyy_0, \
                             ta2_xx_xxyz_yyyz_0, \
                             ta2_xx_xxyz_yyzz_0, \
                             ta2_xx_xxyz_yzzz_0, \
                             ta2_xx_xxyz_zzzz_0, \
                             ta2_xx_xxz_xxxx_0,  \
                             ta2_xx_xxz_xxxx_1,  \
                             ta2_xx_xxz_xxxz_0,  \
                             ta2_xx_xxz_xxxz_1,  \
                             ta2_xx_xxz_xxyz_0,  \
                             ta2_xx_xxz_xxyz_1,  \
                             ta2_xx_xxz_xxz_0,   \
                             ta2_xx_xxz_xxz_1,   \
                             ta2_xx_xxz_xxzz_0,  \
                             ta2_xx_xxz_xxzz_1,  \
                             ta2_xx_xxz_xyyz_0,  \
                             ta2_xx_xxz_xyyz_1,  \
                             ta2_xx_xxz_xyz_0,   \
                             ta2_xx_xxz_xyz_1,   \
                             ta2_xx_xxz_xyzz_0,  \
                             ta2_xx_xxz_xyzz_1,  \
                             ta2_xx_xxz_xzz_0,   \
                             ta2_xx_xxz_xzz_1,   \
                             ta2_xx_xxz_xzzz_0,  \
                             ta2_xx_xxz_xzzz_1,  \
                             ta2_xx_xxz_yyyz_0,  \
                             ta2_xx_xxz_yyyz_1,  \
                             ta2_xx_xxz_yyz_0,   \
                             ta2_xx_xxz_yyz_1,   \
                             ta2_xx_xxz_yyzz_0,  \
                             ta2_xx_xxz_yyzz_1,  \
                             ta2_xx_xxz_yzz_0,   \
                             ta2_xx_xxz_yzz_1,   \
                             ta2_xx_xxz_yzzz_0,  \
                             ta2_xx_xxz_yzzz_1,  \
                             ta2_xx_xxz_zzz_0,   \
                             ta2_xx_xxz_zzz_1,   \
                             ta2_xx_xxz_zzzz_0,  \
                             ta2_xx_xxz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxyz_xxxx_0[i] = ta2_xx_xxz_xxxx_0[i] * pa_y[i] - ta2_xx_xxz_xxxx_1[i] * pc_y[i];

        ta2_xx_xxyz_xxxy_0[i] = ta2_xx_xxy_xxxy_0[i] * pa_z[i] - ta2_xx_xxy_xxxy_1[i] * pc_z[i];

        ta2_xx_xxyz_xxxz_0[i] = ta2_xx_xxz_xxxz_0[i] * pa_y[i] - ta2_xx_xxz_xxxz_1[i] * pc_y[i];

        ta2_xx_xxyz_xxyy_0[i] = ta2_xx_xxy_xxyy_0[i] * pa_z[i] - ta2_xx_xxy_xxyy_1[i] * pc_z[i];

        ta2_xx_xxyz_xxyz_0[i] =
            ta2_xx_xxz_xxz_0[i] * fe_0 - ta2_xx_xxz_xxz_1[i] * fe_0 + ta2_xx_xxz_xxyz_0[i] * pa_y[i] - ta2_xx_xxz_xxyz_1[i] * pc_y[i];

        ta2_xx_xxyz_xxzz_0[i] = ta2_xx_xxz_xxzz_0[i] * pa_y[i] - ta2_xx_xxz_xxzz_1[i] * pc_y[i];

        ta2_xx_xxyz_xyyy_0[i] = ta2_xx_xxy_xyyy_0[i] * pa_z[i] - ta2_xx_xxy_xyyy_1[i] * pc_z[i];

        ta2_xx_xxyz_xyyz_0[i] =
            2.0 * ta2_xx_xxz_xyz_0[i] * fe_0 - 2.0 * ta2_xx_xxz_xyz_1[i] * fe_0 + ta2_xx_xxz_xyyz_0[i] * pa_y[i] - ta2_xx_xxz_xyyz_1[i] * pc_y[i];

        ta2_xx_xxyz_xyzz_0[i] =
            ta2_xx_xxz_xzz_0[i] * fe_0 - ta2_xx_xxz_xzz_1[i] * fe_0 + ta2_xx_xxz_xyzz_0[i] * pa_y[i] - ta2_xx_xxz_xyzz_1[i] * pc_y[i];

        ta2_xx_xxyz_xzzz_0[i] = ta2_xx_xxz_xzzz_0[i] * pa_y[i] - ta2_xx_xxz_xzzz_1[i] * pc_y[i];

        ta2_xx_xxyz_yyyy_0[i] = ta2_xx_xxy_yyyy_0[i] * pa_z[i] - ta2_xx_xxy_yyyy_1[i] * pc_z[i];

        ta2_xx_xxyz_yyyz_0[i] =
            3.0 * ta2_xx_xxz_yyz_0[i] * fe_0 - 3.0 * ta2_xx_xxz_yyz_1[i] * fe_0 + ta2_xx_xxz_yyyz_0[i] * pa_y[i] - ta2_xx_xxz_yyyz_1[i] * pc_y[i];

        ta2_xx_xxyz_yyzz_0[i] =
            2.0 * ta2_xx_xxz_yzz_0[i] * fe_0 - 2.0 * ta2_xx_xxz_yzz_1[i] * fe_0 + ta2_xx_xxz_yyzz_0[i] * pa_y[i] - ta2_xx_xxz_yyzz_1[i] * pc_y[i];

        ta2_xx_xxyz_yzzz_0[i] =
            ta2_xx_xxz_zzz_0[i] * fe_0 - ta2_xx_xxz_zzz_1[i] * fe_0 + ta2_xx_xxz_yzzz_0[i] * pa_y[i] - ta2_xx_xxz_yzzz_1[i] * pc_y[i];

        ta2_xx_xxyz_zzzz_0[i] = ta2_xx_xxz_zzzz_0[i] * pa_y[i] - ta2_xx_xxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 75-90 components of targeted buffer : GG

    auto ta2_xx_xxzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 75);

    auto ta2_xx_xxzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 76);

    auto ta2_xx_xxzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 77);

    auto ta2_xx_xxzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 78);

    auto ta2_xx_xxzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 79);

    auto ta2_xx_xxzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 80);

    auto ta2_xx_xxzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 81);

    auto ta2_xx_xxzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 82);

    auto ta2_xx_xxzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 83);

    auto ta2_xx_xxzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 84);

    auto ta2_xx_xxzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 85);

    auto ta2_xx_xxzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 86);

    auto ta2_xx_xxzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 87);

    auto ta2_xx_xxzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 88);

    auto ta2_xx_xxzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 89);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_x_xzz_yyyz_1,   \
                             ta1_x_xzz_yyzz_1,   \
                             ta1_x_xzz_yzzz_1,   \
                             ta1_x_xzz_zzzz_1,   \
                             ta2_xx_xx_xxxx_0,   \
                             ta2_xx_xx_xxxx_1,   \
                             ta2_xx_xx_xxxy_0,   \
                             ta2_xx_xx_xxxy_1,   \
                             ta2_xx_xx_xxxz_0,   \
                             ta2_xx_xx_xxxz_1,   \
                             ta2_xx_xx_xxyy_0,   \
                             ta2_xx_xx_xxyy_1,   \
                             ta2_xx_xx_xxyz_0,   \
                             ta2_xx_xx_xxyz_1,   \
                             ta2_xx_xx_xxzz_0,   \
                             ta2_xx_xx_xxzz_1,   \
                             ta2_xx_xx_xyyy_0,   \
                             ta2_xx_xx_xyyy_1,   \
                             ta2_xx_xx_xyyz_0,   \
                             ta2_xx_xx_xyyz_1,   \
                             ta2_xx_xx_xyzz_0,   \
                             ta2_xx_xx_xyzz_1,   \
                             ta2_xx_xx_xzzz_0,   \
                             ta2_xx_xx_xzzz_1,   \
                             ta2_xx_xx_yyyy_0,   \
                             ta2_xx_xx_yyyy_1,   \
                             ta2_xx_xxz_xxx_0,   \
                             ta2_xx_xxz_xxx_1,   \
                             ta2_xx_xxz_xxxx_0,  \
                             ta2_xx_xxz_xxxx_1,  \
                             ta2_xx_xxz_xxxy_0,  \
                             ta2_xx_xxz_xxxy_1,  \
                             ta2_xx_xxz_xxxz_0,  \
                             ta2_xx_xxz_xxxz_1,  \
                             ta2_xx_xxz_xxy_0,   \
                             ta2_xx_xxz_xxy_1,   \
                             ta2_xx_xxz_xxyy_0,  \
                             ta2_xx_xxz_xxyy_1,  \
                             ta2_xx_xxz_xxyz_0,  \
                             ta2_xx_xxz_xxyz_1,  \
                             ta2_xx_xxz_xxz_0,   \
                             ta2_xx_xxz_xxz_1,   \
                             ta2_xx_xxz_xxzz_0,  \
                             ta2_xx_xxz_xxzz_1,  \
                             ta2_xx_xxz_xyy_0,   \
                             ta2_xx_xxz_xyy_1,   \
                             ta2_xx_xxz_xyyy_0,  \
                             ta2_xx_xxz_xyyy_1,  \
                             ta2_xx_xxz_xyyz_0,  \
                             ta2_xx_xxz_xyyz_1,  \
                             ta2_xx_xxz_xyz_0,   \
                             ta2_xx_xxz_xyz_1,   \
                             ta2_xx_xxz_xyzz_0,  \
                             ta2_xx_xxz_xyzz_1,  \
                             ta2_xx_xxz_xzz_0,   \
                             ta2_xx_xxz_xzz_1,   \
                             ta2_xx_xxz_xzzz_0,  \
                             ta2_xx_xxz_xzzz_1,  \
                             ta2_xx_xxz_yyyy_0,  \
                             ta2_xx_xxz_yyyy_1,  \
                             ta2_xx_xxzz_xxxx_0, \
                             ta2_xx_xxzz_xxxy_0, \
                             ta2_xx_xxzz_xxxz_0, \
                             ta2_xx_xxzz_xxyy_0, \
                             ta2_xx_xxzz_xxyz_0, \
                             ta2_xx_xxzz_xxzz_0, \
                             ta2_xx_xxzz_xyyy_0, \
                             ta2_xx_xxzz_xyyz_0, \
                             ta2_xx_xxzz_xyzz_0, \
                             ta2_xx_xxzz_xzzz_0, \
                             ta2_xx_xxzz_yyyy_0, \
                             ta2_xx_xxzz_yyyz_0, \
                             ta2_xx_xxzz_yyzz_0, \
                             ta2_xx_xxzz_yzzz_0, \
                             ta2_xx_xxzz_zzzz_0, \
                             ta2_xx_xzz_yyyz_0,  \
                             ta2_xx_xzz_yyyz_1,  \
                             ta2_xx_xzz_yyzz_0,  \
                             ta2_xx_xzz_yyzz_1,  \
                             ta2_xx_xzz_yzzz_0,  \
                             ta2_xx_xzz_yzzz_1,  \
                             ta2_xx_xzz_zzzz_0,  \
                             ta2_xx_xzz_zzzz_1,  \
                             ta2_xx_zz_yyyz_0,   \
                             ta2_xx_zz_yyyz_1,   \
                             ta2_xx_zz_yyzz_0,   \
                             ta2_xx_zz_yyzz_1,   \
                             ta2_xx_zz_yzzz_0,   \
                             ta2_xx_zz_yzzz_1,   \
                             ta2_xx_zz_zzzz_0,   \
                             ta2_xx_zz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xxzz_xxxx_0[i] =
            ta2_xx_xx_xxxx_0[i] * fe_0 - ta2_xx_xx_xxxx_1[i] * fe_0 + ta2_xx_xxz_xxxx_0[i] * pa_z[i] - ta2_xx_xxz_xxxx_1[i] * pc_z[i];

        ta2_xx_xxzz_xxxy_0[i] =
            ta2_xx_xx_xxxy_0[i] * fe_0 - ta2_xx_xx_xxxy_1[i] * fe_0 + ta2_xx_xxz_xxxy_0[i] * pa_z[i] - ta2_xx_xxz_xxxy_1[i] * pc_z[i];

        ta2_xx_xxzz_xxxz_0[i] = ta2_xx_xx_xxxz_0[i] * fe_0 - ta2_xx_xx_xxxz_1[i] * fe_0 + ta2_xx_xxz_xxx_0[i] * fe_0 - ta2_xx_xxz_xxx_1[i] * fe_0 +
                                ta2_xx_xxz_xxxz_0[i] * pa_z[i] - ta2_xx_xxz_xxxz_1[i] * pc_z[i];

        ta2_xx_xxzz_xxyy_0[i] =
            ta2_xx_xx_xxyy_0[i] * fe_0 - ta2_xx_xx_xxyy_1[i] * fe_0 + ta2_xx_xxz_xxyy_0[i] * pa_z[i] - ta2_xx_xxz_xxyy_1[i] * pc_z[i];

        ta2_xx_xxzz_xxyz_0[i] = ta2_xx_xx_xxyz_0[i] * fe_0 - ta2_xx_xx_xxyz_1[i] * fe_0 + ta2_xx_xxz_xxy_0[i] * fe_0 - ta2_xx_xxz_xxy_1[i] * fe_0 +
                                ta2_xx_xxz_xxyz_0[i] * pa_z[i] - ta2_xx_xxz_xxyz_1[i] * pc_z[i];

        ta2_xx_xxzz_xxzz_0[i] = ta2_xx_xx_xxzz_0[i] * fe_0 - ta2_xx_xx_xxzz_1[i] * fe_0 + 2.0 * ta2_xx_xxz_xxz_0[i] * fe_0 -
                                2.0 * ta2_xx_xxz_xxz_1[i] * fe_0 + ta2_xx_xxz_xxzz_0[i] * pa_z[i] - ta2_xx_xxz_xxzz_1[i] * pc_z[i];

        ta2_xx_xxzz_xyyy_0[i] =
            ta2_xx_xx_xyyy_0[i] * fe_0 - ta2_xx_xx_xyyy_1[i] * fe_0 + ta2_xx_xxz_xyyy_0[i] * pa_z[i] - ta2_xx_xxz_xyyy_1[i] * pc_z[i];

        ta2_xx_xxzz_xyyz_0[i] = ta2_xx_xx_xyyz_0[i] * fe_0 - ta2_xx_xx_xyyz_1[i] * fe_0 + ta2_xx_xxz_xyy_0[i] * fe_0 - ta2_xx_xxz_xyy_1[i] * fe_0 +
                                ta2_xx_xxz_xyyz_0[i] * pa_z[i] - ta2_xx_xxz_xyyz_1[i] * pc_z[i];

        ta2_xx_xxzz_xyzz_0[i] = ta2_xx_xx_xyzz_0[i] * fe_0 - ta2_xx_xx_xyzz_1[i] * fe_0 + 2.0 * ta2_xx_xxz_xyz_0[i] * fe_0 -
                                2.0 * ta2_xx_xxz_xyz_1[i] * fe_0 + ta2_xx_xxz_xyzz_0[i] * pa_z[i] - ta2_xx_xxz_xyzz_1[i] * pc_z[i];

        ta2_xx_xxzz_xzzz_0[i] = ta2_xx_xx_xzzz_0[i] * fe_0 - ta2_xx_xx_xzzz_1[i] * fe_0 + 3.0 * ta2_xx_xxz_xzz_0[i] * fe_0 -
                                3.0 * ta2_xx_xxz_xzz_1[i] * fe_0 + ta2_xx_xxz_xzzz_0[i] * pa_z[i] - ta2_xx_xxz_xzzz_1[i] * pc_z[i];

        ta2_xx_xxzz_yyyy_0[i] =
            ta2_xx_xx_yyyy_0[i] * fe_0 - ta2_xx_xx_yyyy_1[i] * fe_0 + ta2_xx_xxz_yyyy_0[i] * pa_z[i] - ta2_xx_xxz_yyyy_1[i] * pc_z[i];

        ta2_xx_xxzz_yyyz_0[i] = ta2_xx_zz_yyyz_0[i] * fe_0 - ta2_xx_zz_yyyz_1[i] * fe_0 + 2.0 * ta1_x_xzz_yyyz_1[i] + ta2_xx_xzz_yyyz_0[i] * pa_x[i] -
                                ta2_xx_xzz_yyyz_1[i] * pc_x[i];

        ta2_xx_xxzz_yyzz_0[i] = ta2_xx_zz_yyzz_0[i] * fe_0 - ta2_xx_zz_yyzz_1[i] * fe_0 + 2.0 * ta1_x_xzz_yyzz_1[i] + ta2_xx_xzz_yyzz_0[i] * pa_x[i] -
                                ta2_xx_xzz_yyzz_1[i] * pc_x[i];

        ta2_xx_xxzz_yzzz_0[i] = ta2_xx_zz_yzzz_0[i] * fe_0 - ta2_xx_zz_yzzz_1[i] * fe_0 + 2.0 * ta1_x_xzz_yzzz_1[i] + ta2_xx_xzz_yzzz_0[i] * pa_x[i] -
                                ta2_xx_xzz_yzzz_1[i] * pc_x[i];

        ta2_xx_xxzz_zzzz_0[i] = ta2_xx_zz_zzzz_0[i] * fe_0 - ta2_xx_zz_zzzz_1[i] * fe_0 + 2.0 * ta1_x_xzz_zzzz_1[i] + ta2_xx_xzz_zzzz_0[i] * pa_x[i] -
                                ta2_xx_xzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 90-105 components of targeted buffer : GG

    auto ta2_xx_xyyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 90);

    auto ta2_xx_xyyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 91);

    auto ta2_xx_xyyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 92);

    auto ta2_xx_xyyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 93);

    auto ta2_xx_xyyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 94);

    auto ta2_xx_xyyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 95);

    auto ta2_xx_xyyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 96);

    auto ta2_xx_xyyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 97);

    auto ta2_xx_xyyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 98);

    auto ta2_xx_xyyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 99);

    auto ta2_xx_xyyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 100);

    auto ta2_xx_xyyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 101);

    auto ta2_xx_xyyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 102);

    auto ta2_xx_xyyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 103);

    auto ta2_xx_xyyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 104);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_yyy_xxxy_1,   \
                             ta1_x_yyy_xxyy_1,   \
                             ta1_x_yyy_xxyz_1,   \
                             ta1_x_yyy_xyyy_1,   \
                             ta1_x_yyy_xyyz_1,   \
                             ta1_x_yyy_xyzz_1,   \
                             ta1_x_yyy_yyyy_1,   \
                             ta1_x_yyy_yyyz_1,   \
                             ta1_x_yyy_yyzz_1,   \
                             ta1_x_yyy_yzzz_1,   \
                             ta1_x_yyy_zzzz_1,   \
                             ta2_xx_xy_xxxx_0,   \
                             ta2_xx_xy_xxxx_1,   \
                             ta2_xx_xy_xxxz_0,   \
                             ta2_xx_xy_xxxz_1,   \
                             ta2_xx_xy_xxzz_0,   \
                             ta2_xx_xy_xxzz_1,   \
                             ta2_xx_xy_xzzz_0,   \
                             ta2_xx_xy_xzzz_1,   \
                             ta2_xx_xyy_xxxx_0,  \
                             ta2_xx_xyy_xxxx_1,  \
                             ta2_xx_xyy_xxxz_0,  \
                             ta2_xx_xyy_xxxz_1,  \
                             ta2_xx_xyy_xxzz_0,  \
                             ta2_xx_xyy_xxzz_1,  \
                             ta2_xx_xyy_xzzz_0,  \
                             ta2_xx_xyy_xzzz_1,  \
                             ta2_xx_xyyy_xxxx_0, \
                             ta2_xx_xyyy_xxxy_0, \
                             ta2_xx_xyyy_xxxz_0, \
                             ta2_xx_xyyy_xxyy_0, \
                             ta2_xx_xyyy_xxyz_0, \
                             ta2_xx_xyyy_xxzz_0, \
                             ta2_xx_xyyy_xyyy_0, \
                             ta2_xx_xyyy_xyyz_0, \
                             ta2_xx_xyyy_xyzz_0, \
                             ta2_xx_xyyy_xzzz_0, \
                             ta2_xx_xyyy_yyyy_0, \
                             ta2_xx_xyyy_yyyz_0, \
                             ta2_xx_xyyy_yyzz_0, \
                             ta2_xx_xyyy_yzzz_0, \
                             ta2_xx_xyyy_zzzz_0, \
                             ta2_xx_yyy_xxxy_0,  \
                             ta2_xx_yyy_xxxy_1,  \
                             ta2_xx_yyy_xxy_0,   \
                             ta2_xx_yyy_xxy_1,   \
                             ta2_xx_yyy_xxyy_0,  \
                             ta2_xx_yyy_xxyy_1,  \
                             ta2_xx_yyy_xxyz_0,  \
                             ta2_xx_yyy_xxyz_1,  \
                             ta2_xx_yyy_xyy_0,   \
                             ta2_xx_yyy_xyy_1,   \
                             ta2_xx_yyy_xyyy_0,  \
                             ta2_xx_yyy_xyyy_1,  \
                             ta2_xx_yyy_xyyz_0,  \
                             ta2_xx_yyy_xyyz_1,  \
                             ta2_xx_yyy_xyz_0,   \
                             ta2_xx_yyy_xyz_1,   \
                             ta2_xx_yyy_xyzz_0,  \
                             ta2_xx_yyy_xyzz_1,  \
                             ta2_xx_yyy_yyy_0,   \
                             ta2_xx_yyy_yyy_1,   \
                             ta2_xx_yyy_yyyy_0,  \
                             ta2_xx_yyy_yyyy_1,  \
                             ta2_xx_yyy_yyyz_0,  \
                             ta2_xx_yyy_yyyz_1,  \
                             ta2_xx_yyy_yyz_0,   \
                             ta2_xx_yyy_yyz_1,   \
                             ta2_xx_yyy_yyzz_0,  \
                             ta2_xx_yyy_yyzz_1,  \
                             ta2_xx_yyy_yzz_0,   \
                             ta2_xx_yyy_yzz_1,   \
                             ta2_xx_yyy_yzzz_0,  \
                             ta2_xx_yyy_yzzz_1,  \
                             ta2_xx_yyy_zzzz_0,  \
                             ta2_xx_yyy_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xyyy_xxxx_0[i] =
            2.0 * ta2_xx_xy_xxxx_0[i] * fe_0 - 2.0 * ta2_xx_xy_xxxx_1[i] * fe_0 + ta2_xx_xyy_xxxx_0[i] * pa_y[i] - ta2_xx_xyy_xxxx_1[i] * pc_y[i];

        ta2_xx_xyyy_xxxy_0[i] = 3.0 * ta2_xx_yyy_xxy_0[i] * fe_0 - 3.0 * ta2_xx_yyy_xxy_1[i] * fe_0 + 2.0 * ta1_x_yyy_xxxy_1[i] +
                                ta2_xx_yyy_xxxy_0[i] * pa_x[i] - ta2_xx_yyy_xxxy_1[i] * pc_x[i];

        ta2_xx_xyyy_xxxz_0[i] =
            2.0 * ta2_xx_xy_xxxz_0[i] * fe_0 - 2.0 * ta2_xx_xy_xxxz_1[i] * fe_0 + ta2_xx_xyy_xxxz_0[i] * pa_y[i] - ta2_xx_xyy_xxxz_1[i] * pc_y[i];

        ta2_xx_xyyy_xxyy_0[i] = 2.0 * ta2_xx_yyy_xyy_0[i] * fe_0 - 2.0 * ta2_xx_yyy_xyy_1[i] * fe_0 + 2.0 * ta1_x_yyy_xxyy_1[i] +
                                ta2_xx_yyy_xxyy_0[i] * pa_x[i] - ta2_xx_yyy_xxyy_1[i] * pc_x[i];

        ta2_xx_xyyy_xxyz_0[i] = 2.0 * ta2_xx_yyy_xyz_0[i] * fe_0 - 2.0 * ta2_xx_yyy_xyz_1[i] * fe_0 + 2.0 * ta1_x_yyy_xxyz_1[i] +
                                ta2_xx_yyy_xxyz_0[i] * pa_x[i] - ta2_xx_yyy_xxyz_1[i] * pc_x[i];

        ta2_xx_xyyy_xxzz_0[i] =
            2.0 * ta2_xx_xy_xxzz_0[i] * fe_0 - 2.0 * ta2_xx_xy_xxzz_1[i] * fe_0 + ta2_xx_xyy_xxzz_0[i] * pa_y[i] - ta2_xx_xyy_xxzz_1[i] * pc_y[i];

        ta2_xx_xyyy_xyyy_0[i] = ta2_xx_yyy_yyy_0[i] * fe_0 - ta2_xx_yyy_yyy_1[i] * fe_0 + 2.0 * ta1_x_yyy_xyyy_1[i] + ta2_xx_yyy_xyyy_0[i] * pa_x[i] -
                                ta2_xx_yyy_xyyy_1[i] * pc_x[i];

        ta2_xx_xyyy_xyyz_0[i] = ta2_xx_yyy_yyz_0[i] * fe_0 - ta2_xx_yyy_yyz_1[i] * fe_0 + 2.0 * ta1_x_yyy_xyyz_1[i] + ta2_xx_yyy_xyyz_0[i] * pa_x[i] -
                                ta2_xx_yyy_xyyz_1[i] * pc_x[i];

        ta2_xx_xyyy_xyzz_0[i] = ta2_xx_yyy_yzz_0[i] * fe_0 - ta2_xx_yyy_yzz_1[i] * fe_0 + 2.0 * ta1_x_yyy_xyzz_1[i] + ta2_xx_yyy_xyzz_0[i] * pa_x[i] -
                                ta2_xx_yyy_xyzz_1[i] * pc_x[i];

        ta2_xx_xyyy_xzzz_0[i] =
            2.0 * ta2_xx_xy_xzzz_0[i] * fe_0 - 2.0 * ta2_xx_xy_xzzz_1[i] * fe_0 + ta2_xx_xyy_xzzz_0[i] * pa_y[i] - ta2_xx_xyy_xzzz_1[i] * pc_y[i];

        ta2_xx_xyyy_yyyy_0[i] = 2.0 * ta1_x_yyy_yyyy_1[i] + ta2_xx_yyy_yyyy_0[i] * pa_x[i] - ta2_xx_yyy_yyyy_1[i] * pc_x[i];

        ta2_xx_xyyy_yyyz_0[i] = 2.0 * ta1_x_yyy_yyyz_1[i] + ta2_xx_yyy_yyyz_0[i] * pa_x[i] - ta2_xx_yyy_yyyz_1[i] * pc_x[i];

        ta2_xx_xyyy_yyzz_0[i] = 2.0 * ta1_x_yyy_yyzz_1[i] + ta2_xx_yyy_yyzz_0[i] * pa_x[i] - ta2_xx_yyy_yyzz_1[i] * pc_x[i];

        ta2_xx_xyyy_yzzz_0[i] = 2.0 * ta1_x_yyy_yzzz_1[i] + ta2_xx_yyy_yzzz_0[i] * pa_x[i] - ta2_xx_yyy_yzzz_1[i] * pc_x[i];

        ta2_xx_xyyy_zzzz_0[i] = 2.0 * ta1_x_yyy_zzzz_1[i] + ta2_xx_yyy_zzzz_0[i] * pa_x[i] - ta2_xx_yyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 105-120 components of targeted buffer : GG

    auto ta2_xx_xyyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 105);

    auto ta2_xx_xyyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 106);

    auto ta2_xx_xyyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 107);

    auto ta2_xx_xyyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 108);

    auto ta2_xx_xyyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 109);

    auto ta2_xx_xyyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 110);

    auto ta2_xx_xyyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 111);

    auto ta2_xx_xyyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 112);

    auto ta2_xx_xyyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 113);

    auto ta2_xx_xyyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 114);

    auto ta2_xx_xyyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 115);

    auto ta2_xx_xyyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 116);

    auto ta2_xx_xyyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 117);

    auto ta2_xx_xyyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 118);

    auto ta2_xx_xyyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 119);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_yyz_yyyz_1,   \
                             ta1_x_yyz_yyzz_1,   \
                             ta1_x_yyz_yzzz_1,   \
                             ta1_x_yyz_zzzz_1,   \
                             ta2_xx_xyy_xxxx_0,  \
                             ta2_xx_xyy_xxxx_1,  \
                             ta2_xx_xyy_xxxy_0,  \
                             ta2_xx_xyy_xxxy_1,  \
                             ta2_xx_xyy_xxy_0,   \
                             ta2_xx_xyy_xxy_1,   \
                             ta2_xx_xyy_xxyy_0,  \
                             ta2_xx_xyy_xxyy_1,  \
                             ta2_xx_xyy_xxyz_0,  \
                             ta2_xx_xyy_xxyz_1,  \
                             ta2_xx_xyy_xyy_0,   \
                             ta2_xx_xyy_xyy_1,   \
                             ta2_xx_xyy_xyyy_0,  \
                             ta2_xx_xyy_xyyy_1,  \
                             ta2_xx_xyy_xyyz_0,  \
                             ta2_xx_xyy_xyyz_1,  \
                             ta2_xx_xyy_xyz_0,   \
                             ta2_xx_xyy_xyz_1,   \
                             ta2_xx_xyy_xyzz_0,  \
                             ta2_xx_xyy_xyzz_1,  \
                             ta2_xx_xyy_yyyy_0,  \
                             ta2_xx_xyy_yyyy_1,  \
                             ta2_xx_xyyz_xxxx_0, \
                             ta2_xx_xyyz_xxxy_0, \
                             ta2_xx_xyyz_xxxz_0, \
                             ta2_xx_xyyz_xxyy_0, \
                             ta2_xx_xyyz_xxyz_0, \
                             ta2_xx_xyyz_xxzz_0, \
                             ta2_xx_xyyz_xyyy_0, \
                             ta2_xx_xyyz_xyyz_0, \
                             ta2_xx_xyyz_xyzz_0, \
                             ta2_xx_xyyz_xzzz_0, \
                             ta2_xx_xyyz_yyyy_0, \
                             ta2_xx_xyyz_yyyz_0, \
                             ta2_xx_xyyz_yyzz_0, \
                             ta2_xx_xyyz_yzzz_0, \
                             ta2_xx_xyyz_zzzz_0, \
                             ta2_xx_xyz_xxxz_0,  \
                             ta2_xx_xyz_xxxz_1,  \
                             ta2_xx_xyz_xxzz_0,  \
                             ta2_xx_xyz_xxzz_1,  \
                             ta2_xx_xyz_xzzz_0,  \
                             ta2_xx_xyz_xzzz_1,  \
                             ta2_xx_xz_xxxz_0,   \
                             ta2_xx_xz_xxxz_1,   \
                             ta2_xx_xz_xxzz_0,   \
                             ta2_xx_xz_xxzz_1,   \
                             ta2_xx_xz_xzzz_0,   \
                             ta2_xx_xz_xzzz_1,   \
                             ta2_xx_yyz_yyyz_0,  \
                             ta2_xx_yyz_yyyz_1,  \
                             ta2_xx_yyz_yyzz_0,  \
                             ta2_xx_yyz_yyzz_1,  \
                             ta2_xx_yyz_yzzz_0,  \
                             ta2_xx_yyz_yzzz_1,  \
                             ta2_xx_yyz_zzzz_0,  \
                             ta2_xx_yyz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xyyz_xxxx_0[i] = ta2_xx_xyy_xxxx_0[i] * pa_z[i] - ta2_xx_xyy_xxxx_1[i] * pc_z[i];

        ta2_xx_xyyz_xxxy_0[i] = ta2_xx_xyy_xxxy_0[i] * pa_z[i] - ta2_xx_xyy_xxxy_1[i] * pc_z[i];

        ta2_xx_xyyz_xxxz_0[i] =
            ta2_xx_xz_xxxz_0[i] * fe_0 - ta2_xx_xz_xxxz_1[i] * fe_0 + ta2_xx_xyz_xxxz_0[i] * pa_y[i] - ta2_xx_xyz_xxxz_1[i] * pc_y[i];

        ta2_xx_xyyz_xxyy_0[i] = ta2_xx_xyy_xxyy_0[i] * pa_z[i] - ta2_xx_xyy_xxyy_1[i] * pc_z[i];

        ta2_xx_xyyz_xxyz_0[i] =
            ta2_xx_xyy_xxy_0[i] * fe_0 - ta2_xx_xyy_xxy_1[i] * fe_0 + ta2_xx_xyy_xxyz_0[i] * pa_z[i] - ta2_xx_xyy_xxyz_1[i] * pc_z[i];

        ta2_xx_xyyz_xxzz_0[i] =
            ta2_xx_xz_xxzz_0[i] * fe_0 - ta2_xx_xz_xxzz_1[i] * fe_0 + ta2_xx_xyz_xxzz_0[i] * pa_y[i] - ta2_xx_xyz_xxzz_1[i] * pc_y[i];

        ta2_xx_xyyz_xyyy_0[i] = ta2_xx_xyy_xyyy_0[i] * pa_z[i] - ta2_xx_xyy_xyyy_1[i] * pc_z[i];

        ta2_xx_xyyz_xyyz_0[i] =
            ta2_xx_xyy_xyy_0[i] * fe_0 - ta2_xx_xyy_xyy_1[i] * fe_0 + ta2_xx_xyy_xyyz_0[i] * pa_z[i] - ta2_xx_xyy_xyyz_1[i] * pc_z[i];

        ta2_xx_xyyz_xyzz_0[i] =
            2.0 * ta2_xx_xyy_xyz_0[i] * fe_0 - 2.0 * ta2_xx_xyy_xyz_1[i] * fe_0 + ta2_xx_xyy_xyzz_0[i] * pa_z[i] - ta2_xx_xyy_xyzz_1[i] * pc_z[i];

        ta2_xx_xyyz_xzzz_0[i] =
            ta2_xx_xz_xzzz_0[i] * fe_0 - ta2_xx_xz_xzzz_1[i] * fe_0 + ta2_xx_xyz_xzzz_0[i] * pa_y[i] - ta2_xx_xyz_xzzz_1[i] * pc_y[i];

        ta2_xx_xyyz_yyyy_0[i] = ta2_xx_xyy_yyyy_0[i] * pa_z[i] - ta2_xx_xyy_yyyy_1[i] * pc_z[i];

        ta2_xx_xyyz_yyyz_0[i] = 2.0 * ta1_x_yyz_yyyz_1[i] + ta2_xx_yyz_yyyz_0[i] * pa_x[i] - ta2_xx_yyz_yyyz_1[i] * pc_x[i];

        ta2_xx_xyyz_yyzz_0[i] = 2.0 * ta1_x_yyz_yyzz_1[i] + ta2_xx_yyz_yyzz_0[i] * pa_x[i] - ta2_xx_yyz_yyzz_1[i] * pc_x[i];

        ta2_xx_xyyz_yzzz_0[i] = 2.0 * ta1_x_yyz_yzzz_1[i] + ta2_xx_yyz_yzzz_0[i] * pa_x[i] - ta2_xx_yyz_yzzz_1[i] * pc_x[i];

        ta2_xx_xyyz_zzzz_0[i] = 2.0 * ta1_x_yyz_zzzz_1[i] + ta2_xx_yyz_zzzz_0[i] * pa_x[i] - ta2_xx_yyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 120-135 components of targeted buffer : GG

    auto ta2_xx_xyzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 120);

    auto ta2_xx_xyzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 121);

    auto ta2_xx_xyzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 122);

    auto ta2_xx_xyzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 123);

    auto ta2_xx_xyzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 124);

    auto ta2_xx_xyzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 125);

    auto ta2_xx_xyzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 126);

    auto ta2_xx_xyzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 127);

    auto ta2_xx_xyzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 128);

    auto ta2_xx_xyzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 129);

    auto ta2_xx_xyzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 130);

    auto ta2_xx_xyzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 131);

    auto ta2_xx_xyzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 132);

    auto ta2_xx_xyzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 133);

    auto ta2_xx_xyzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 134);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_yzz_yyyy_1,   \
                             ta1_x_yzz_yyyz_1,   \
                             ta1_x_yzz_yyzz_1,   \
                             ta1_x_yzz_yzzz_1,   \
                             ta2_xx_xyzz_xxxx_0, \
                             ta2_xx_xyzz_xxxy_0, \
                             ta2_xx_xyzz_xxxz_0, \
                             ta2_xx_xyzz_xxyy_0, \
                             ta2_xx_xyzz_xxyz_0, \
                             ta2_xx_xyzz_xxzz_0, \
                             ta2_xx_xyzz_xyyy_0, \
                             ta2_xx_xyzz_xyyz_0, \
                             ta2_xx_xyzz_xyzz_0, \
                             ta2_xx_xyzz_xzzz_0, \
                             ta2_xx_xyzz_yyyy_0, \
                             ta2_xx_xyzz_yyyz_0, \
                             ta2_xx_xyzz_yyzz_0, \
                             ta2_xx_xyzz_yzzz_0, \
                             ta2_xx_xyzz_zzzz_0, \
                             ta2_xx_xzz_xxx_0,   \
                             ta2_xx_xzz_xxx_1,   \
                             ta2_xx_xzz_xxxx_0,  \
                             ta2_xx_xzz_xxxx_1,  \
                             ta2_xx_xzz_xxxy_0,  \
                             ta2_xx_xzz_xxxy_1,  \
                             ta2_xx_xzz_xxxz_0,  \
                             ta2_xx_xzz_xxxz_1,  \
                             ta2_xx_xzz_xxy_0,   \
                             ta2_xx_xzz_xxy_1,   \
                             ta2_xx_xzz_xxyy_0,  \
                             ta2_xx_xzz_xxyy_1,  \
                             ta2_xx_xzz_xxyz_0,  \
                             ta2_xx_xzz_xxyz_1,  \
                             ta2_xx_xzz_xxz_0,   \
                             ta2_xx_xzz_xxz_1,   \
                             ta2_xx_xzz_xxzz_0,  \
                             ta2_xx_xzz_xxzz_1,  \
                             ta2_xx_xzz_xyy_0,   \
                             ta2_xx_xzz_xyy_1,   \
                             ta2_xx_xzz_xyyy_0,  \
                             ta2_xx_xzz_xyyy_1,  \
                             ta2_xx_xzz_xyyz_0,  \
                             ta2_xx_xzz_xyyz_1,  \
                             ta2_xx_xzz_xyz_0,   \
                             ta2_xx_xzz_xyz_1,   \
                             ta2_xx_xzz_xyzz_0,  \
                             ta2_xx_xzz_xyzz_1,  \
                             ta2_xx_xzz_xzz_0,   \
                             ta2_xx_xzz_xzz_1,   \
                             ta2_xx_xzz_xzzz_0,  \
                             ta2_xx_xzz_xzzz_1,  \
                             ta2_xx_xzz_zzzz_0,  \
                             ta2_xx_xzz_zzzz_1,  \
                             ta2_xx_yzz_yyyy_0,  \
                             ta2_xx_yzz_yyyy_1,  \
                             ta2_xx_yzz_yyyz_0,  \
                             ta2_xx_yzz_yyyz_1,  \
                             ta2_xx_yzz_yyzz_0,  \
                             ta2_xx_yzz_yyzz_1,  \
                             ta2_xx_yzz_yzzz_0,  \
                             ta2_xx_yzz_yzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xyzz_xxxx_0[i] = ta2_xx_xzz_xxxx_0[i] * pa_y[i] - ta2_xx_xzz_xxxx_1[i] * pc_y[i];

        ta2_xx_xyzz_xxxy_0[i] =
            ta2_xx_xzz_xxx_0[i] * fe_0 - ta2_xx_xzz_xxx_1[i] * fe_0 + ta2_xx_xzz_xxxy_0[i] * pa_y[i] - ta2_xx_xzz_xxxy_1[i] * pc_y[i];

        ta2_xx_xyzz_xxxz_0[i] = ta2_xx_xzz_xxxz_0[i] * pa_y[i] - ta2_xx_xzz_xxxz_1[i] * pc_y[i];

        ta2_xx_xyzz_xxyy_0[i] =
            2.0 * ta2_xx_xzz_xxy_0[i] * fe_0 - 2.0 * ta2_xx_xzz_xxy_1[i] * fe_0 + ta2_xx_xzz_xxyy_0[i] * pa_y[i] - ta2_xx_xzz_xxyy_1[i] * pc_y[i];

        ta2_xx_xyzz_xxyz_0[i] =
            ta2_xx_xzz_xxz_0[i] * fe_0 - ta2_xx_xzz_xxz_1[i] * fe_0 + ta2_xx_xzz_xxyz_0[i] * pa_y[i] - ta2_xx_xzz_xxyz_1[i] * pc_y[i];

        ta2_xx_xyzz_xxzz_0[i] = ta2_xx_xzz_xxzz_0[i] * pa_y[i] - ta2_xx_xzz_xxzz_1[i] * pc_y[i];

        ta2_xx_xyzz_xyyy_0[i] =
            3.0 * ta2_xx_xzz_xyy_0[i] * fe_0 - 3.0 * ta2_xx_xzz_xyy_1[i] * fe_0 + ta2_xx_xzz_xyyy_0[i] * pa_y[i] - ta2_xx_xzz_xyyy_1[i] * pc_y[i];

        ta2_xx_xyzz_xyyz_0[i] =
            2.0 * ta2_xx_xzz_xyz_0[i] * fe_0 - 2.0 * ta2_xx_xzz_xyz_1[i] * fe_0 + ta2_xx_xzz_xyyz_0[i] * pa_y[i] - ta2_xx_xzz_xyyz_1[i] * pc_y[i];

        ta2_xx_xyzz_xyzz_0[i] =
            ta2_xx_xzz_xzz_0[i] * fe_0 - ta2_xx_xzz_xzz_1[i] * fe_0 + ta2_xx_xzz_xyzz_0[i] * pa_y[i] - ta2_xx_xzz_xyzz_1[i] * pc_y[i];

        ta2_xx_xyzz_xzzz_0[i] = ta2_xx_xzz_xzzz_0[i] * pa_y[i] - ta2_xx_xzz_xzzz_1[i] * pc_y[i];

        ta2_xx_xyzz_yyyy_0[i] = 2.0 * ta1_x_yzz_yyyy_1[i] + ta2_xx_yzz_yyyy_0[i] * pa_x[i] - ta2_xx_yzz_yyyy_1[i] * pc_x[i];

        ta2_xx_xyzz_yyyz_0[i] = 2.0 * ta1_x_yzz_yyyz_1[i] + ta2_xx_yzz_yyyz_0[i] * pa_x[i] - ta2_xx_yzz_yyyz_1[i] * pc_x[i];

        ta2_xx_xyzz_yyzz_0[i] = 2.0 * ta1_x_yzz_yyzz_1[i] + ta2_xx_yzz_yyzz_0[i] * pa_x[i] - ta2_xx_yzz_yyzz_1[i] * pc_x[i];

        ta2_xx_xyzz_yzzz_0[i] = 2.0 * ta1_x_yzz_yzzz_1[i] + ta2_xx_yzz_yzzz_0[i] * pa_x[i] - ta2_xx_yzz_yzzz_1[i] * pc_x[i];

        ta2_xx_xyzz_zzzz_0[i] = ta2_xx_xzz_zzzz_0[i] * pa_y[i] - ta2_xx_xzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 135-150 components of targeted buffer : GG

    auto ta2_xx_xzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 135);

    auto ta2_xx_xzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 136);

    auto ta2_xx_xzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 137);

    auto ta2_xx_xzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 138);

    auto ta2_xx_xzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 139);

    auto ta2_xx_xzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 140);

    auto ta2_xx_xzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 141);

    auto ta2_xx_xzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 142);

    auto ta2_xx_xzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 143);

    auto ta2_xx_xzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 144);

    auto ta2_xx_xzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 145);

    auto ta2_xx_xzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 146);

    auto ta2_xx_xzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 147);

    auto ta2_xx_xzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 148);

    auto ta2_xx_xzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 149);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_x_zzz_xxxz_1,   \
                             ta1_x_zzz_xxyz_1,   \
                             ta1_x_zzz_xxzz_1,   \
                             ta1_x_zzz_xyyz_1,   \
                             ta1_x_zzz_xyzz_1,   \
                             ta1_x_zzz_xzzz_1,   \
                             ta1_x_zzz_yyyy_1,   \
                             ta1_x_zzz_yyyz_1,   \
                             ta1_x_zzz_yyzz_1,   \
                             ta1_x_zzz_yzzz_1,   \
                             ta1_x_zzz_zzzz_1,   \
                             ta2_xx_xz_xxxx_0,   \
                             ta2_xx_xz_xxxx_1,   \
                             ta2_xx_xz_xxxy_0,   \
                             ta2_xx_xz_xxxy_1,   \
                             ta2_xx_xz_xxyy_0,   \
                             ta2_xx_xz_xxyy_1,   \
                             ta2_xx_xz_xyyy_0,   \
                             ta2_xx_xz_xyyy_1,   \
                             ta2_xx_xzz_xxxx_0,  \
                             ta2_xx_xzz_xxxx_1,  \
                             ta2_xx_xzz_xxxy_0,  \
                             ta2_xx_xzz_xxxy_1,  \
                             ta2_xx_xzz_xxyy_0,  \
                             ta2_xx_xzz_xxyy_1,  \
                             ta2_xx_xzz_xyyy_0,  \
                             ta2_xx_xzz_xyyy_1,  \
                             ta2_xx_xzzz_xxxx_0, \
                             ta2_xx_xzzz_xxxy_0, \
                             ta2_xx_xzzz_xxxz_0, \
                             ta2_xx_xzzz_xxyy_0, \
                             ta2_xx_xzzz_xxyz_0, \
                             ta2_xx_xzzz_xxzz_0, \
                             ta2_xx_xzzz_xyyy_0, \
                             ta2_xx_xzzz_xyyz_0, \
                             ta2_xx_xzzz_xyzz_0, \
                             ta2_xx_xzzz_xzzz_0, \
                             ta2_xx_xzzz_yyyy_0, \
                             ta2_xx_xzzz_yyyz_0, \
                             ta2_xx_xzzz_yyzz_0, \
                             ta2_xx_xzzz_yzzz_0, \
                             ta2_xx_xzzz_zzzz_0, \
                             ta2_xx_zzz_xxxz_0,  \
                             ta2_xx_zzz_xxxz_1,  \
                             ta2_xx_zzz_xxyz_0,  \
                             ta2_xx_zzz_xxyz_1,  \
                             ta2_xx_zzz_xxz_0,   \
                             ta2_xx_zzz_xxz_1,   \
                             ta2_xx_zzz_xxzz_0,  \
                             ta2_xx_zzz_xxzz_1,  \
                             ta2_xx_zzz_xyyz_0,  \
                             ta2_xx_zzz_xyyz_1,  \
                             ta2_xx_zzz_xyz_0,   \
                             ta2_xx_zzz_xyz_1,   \
                             ta2_xx_zzz_xyzz_0,  \
                             ta2_xx_zzz_xyzz_1,  \
                             ta2_xx_zzz_xzz_0,   \
                             ta2_xx_zzz_xzz_1,   \
                             ta2_xx_zzz_xzzz_0,  \
                             ta2_xx_zzz_xzzz_1,  \
                             ta2_xx_zzz_yyyy_0,  \
                             ta2_xx_zzz_yyyy_1,  \
                             ta2_xx_zzz_yyyz_0,  \
                             ta2_xx_zzz_yyyz_1,  \
                             ta2_xx_zzz_yyz_0,   \
                             ta2_xx_zzz_yyz_1,   \
                             ta2_xx_zzz_yyzz_0,  \
                             ta2_xx_zzz_yyzz_1,  \
                             ta2_xx_zzz_yzz_0,   \
                             ta2_xx_zzz_yzz_1,   \
                             ta2_xx_zzz_yzzz_0,  \
                             ta2_xx_zzz_yzzz_1,  \
                             ta2_xx_zzz_zzz_0,   \
                             ta2_xx_zzz_zzz_1,   \
                             ta2_xx_zzz_zzzz_0,  \
                             ta2_xx_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_xzzz_xxxx_0[i] =
            2.0 * ta2_xx_xz_xxxx_0[i] * fe_0 - 2.0 * ta2_xx_xz_xxxx_1[i] * fe_0 + ta2_xx_xzz_xxxx_0[i] * pa_z[i] - ta2_xx_xzz_xxxx_1[i] * pc_z[i];

        ta2_xx_xzzz_xxxy_0[i] =
            2.0 * ta2_xx_xz_xxxy_0[i] * fe_0 - 2.0 * ta2_xx_xz_xxxy_1[i] * fe_0 + ta2_xx_xzz_xxxy_0[i] * pa_z[i] - ta2_xx_xzz_xxxy_1[i] * pc_z[i];

        ta2_xx_xzzz_xxxz_0[i] = 3.0 * ta2_xx_zzz_xxz_0[i] * fe_0 - 3.0 * ta2_xx_zzz_xxz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xxxz_1[i] +
                                ta2_xx_zzz_xxxz_0[i] * pa_x[i] - ta2_xx_zzz_xxxz_1[i] * pc_x[i];

        ta2_xx_xzzz_xxyy_0[i] =
            2.0 * ta2_xx_xz_xxyy_0[i] * fe_0 - 2.0 * ta2_xx_xz_xxyy_1[i] * fe_0 + ta2_xx_xzz_xxyy_0[i] * pa_z[i] - ta2_xx_xzz_xxyy_1[i] * pc_z[i];

        ta2_xx_xzzz_xxyz_0[i] = 2.0 * ta2_xx_zzz_xyz_0[i] * fe_0 - 2.0 * ta2_xx_zzz_xyz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xxyz_1[i] +
                                ta2_xx_zzz_xxyz_0[i] * pa_x[i] - ta2_xx_zzz_xxyz_1[i] * pc_x[i];

        ta2_xx_xzzz_xxzz_0[i] = 2.0 * ta2_xx_zzz_xzz_0[i] * fe_0 - 2.0 * ta2_xx_zzz_xzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xxzz_1[i] +
                                ta2_xx_zzz_xxzz_0[i] * pa_x[i] - ta2_xx_zzz_xxzz_1[i] * pc_x[i];

        ta2_xx_xzzz_xyyy_0[i] =
            2.0 * ta2_xx_xz_xyyy_0[i] * fe_0 - 2.0 * ta2_xx_xz_xyyy_1[i] * fe_0 + ta2_xx_xzz_xyyy_0[i] * pa_z[i] - ta2_xx_xzz_xyyy_1[i] * pc_z[i];

        ta2_xx_xzzz_xyyz_0[i] = ta2_xx_zzz_yyz_0[i] * fe_0 - ta2_xx_zzz_yyz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xyyz_1[i] + ta2_xx_zzz_xyyz_0[i] * pa_x[i] -
                                ta2_xx_zzz_xyyz_1[i] * pc_x[i];

        ta2_xx_xzzz_xyzz_0[i] = ta2_xx_zzz_yzz_0[i] * fe_0 - ta2_xx_zzz_yzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xyzz_1[i] + ta2_xx_zzz_xyzz_0[i] * pa_x[i] -
                                ta2_xx_zzz_xyzz_1[i] * pc_x[i];

        ta2_xx_xzzz_xzzz_0[i] = ta2_xx_zzz_zzz_0[i] * fe_0 - ta2_xx_zzz_zzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xzzz_1[i] + ta2_xx_zzz_xzzz_0[i] * pa_x[i] -
                                ta2_xx_zzz_xzzz_1[i] * pc_x[i];

        ta2_xx_xzzz_yyyy_0[i] = 2.0 * ta1_x_zzz_yyyy_1[i] + ta2_xx_zzz_yyyy_0[i] * pa_x[i] - ta2_xx_zzz_yyyy_1[i] * pc_x[i];

        ta2_xx_xzzz_yyyz_0[i] = 2.0 * ta1_x_zzz_yyyz_1[i] + ta2_xx_zzz_yyyz_0[i] * pa_x[i] - ta2_xx_zzz_yyyz_1[i] * pc_x[i];

        ta2_xx_xzzz_yyzz_0[i] = 2.0 * ta1_x_zzz_yyzz_1[i] + ta2_xx_zzz_yyzz_0[i] * pa_x[i] - ta2_xx_zzz_yyzz_1[i] * pc_x[i];

        ta2_xx_xzzz_yzzz_0[i] = 2.0 * ta1_x_zzz_yzzz_1[i] + ta2_xx_zzz_yzzz_0[i] * pa_x[i] - ta2_xx_zzz_yzzz_1[i] * pc_x[i];

        ta2_xx_xzzz_zzzz_0[i] = 2.0 * ta1_x_zzz_zzzz_1[i] + ta2_xx_zzz_zzzz_0[i] * pa_x[i] - ta2_xx_zzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 150-165 components of targeted buffer : GG

    auto ta2_xx_yyyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 150);

    auto ta2_xx_yyyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 151);

    auto ta2_xx_yyyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 152);

    auto ta2_xx_yyyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 153);

    auto ta2_xx_yyyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 154);

    auto ta2_xx_yyyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 155);

    auto ta2_xx_yyyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 156);

    auto ta2_xx_yyyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 157);

    auto ta2_xx_yyyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 158);

    auto ta2_xx_yyyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 159);

    auto ta2_xx_yyyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 160);

    auto ta2_xx_yyyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 161);

    auto ta2_xx_yyyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 162);

    auto ta2_xx_yyyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 163);

    auto ta2_xx_yyyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 164);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta2_xx_yy_xxxx_0,   \
                             ta2_xx_yy_xxxx_1,   \
                             ta2_xx_yy_xxxy_0,   \
                             ta2_xx_yy_xxxy_1,   \
                             ta2_xx_yy_xxxz_0,   \
                             ta2_xx_yy_xxxz_1,   \
                             ta2_xx_yy_xxyy_0,   \
                             ta2_xx_yy_xxyy_1,   \
                             ta2_xx_yy_xxyz_0,   \
                             ta2_xx_yy_xxyz_1,   \
                             ta2_xx_yy_xxzz_0,   \
                             ta2_xx_yy_xxzz_1,   \
                             ta2_xx_yy_xyyy_0,   \
                             ta2_xx_yy_xyyy_1,   \
                             ta2_xx_yy_xyyz_0,   \
                             ta2_xx_yy_xyyz_1,   \
                             ta2_xx_yy_xyzz_0,   \
                             ta2_xx_yy_xyzz_1,   \
                             ta2_xx_yy_xzzz_0,   \
                             ta2_xx_yy_xzzz_1,   \
                             ta2_xx_yy_yyyy_0,   \
                             ta2_xx_yy_yyyy_1,   \
                             ta2_xx_yy_yyyz_0,   \
                             ta2_xx_yy_yyyz_1,   \
                             ta2_xx_yy_yyzz_0,   \
                             ta2_xx_yy_yyzz_1,   \
                             ta2_xx_yy_yzzz_0,   \
                             ta2_xx_yy_yzzz_1,   \
                             ta2_xx_yy_zzzz_0,   \
                             ta2_xx_yy_zzzz_1,   \
                             ta2_xx_yyy_xxx_0,   \
                             ta2_xx_yyy_xxx_1,   \
                             ta2_xx_yyy_xxxx_0,  \
                             ta2_xx_yyy_xxxx_1,  \
                             ta2_xx_yyy_xxxy_0,  \
                             ta2_xx_yyy_xxxy_1,  \
                             ta2_xx_yyy_xxxz_0,  \
                             ta2_xx_yyy_xxxz_1,  \
                             ta2_xx_yyy_xxy_0,   \
                             ta2_xx_yyy_xxy_1,   \
                             ta2_xx_yyy_xxyy_0,  \
                             ta2_xx_yyy_xxyy_1,  \
                             ta2_xx_yyy_xxyz_0,  \
                             ta2_xx_yyy_xxyz_1,  \
                             ta2_xx_yyy_xxz_0,   \
                             ta2_xx_yyy_xxz_1,   \
                             ta2_xx_yyy_xxzz_0,  \
                             ta2_xx_yyy_xxzz_1,  \
                             ta2_xx_yyy_xyy_0,   \
                             ta2_xx_yyy_xyy_1,   \
                             ta2_xx_yyy_xyyy_0,  \
                             ta2_xx_yyy_xyyy_1,  \
                             ta2_xx_yyy_xyyz_0,  \
                             ta2_xx_yyy_xyyz_1,  \
                             ta2_xx_yyy_xyz_0,   \
                             ta2_xx_yyy_xyz_1,   \
                             ta2_xx_yyy_xyzz_0,  \
                             ta2_xx_yyy_xyzz_1,  \
                             ta2_xx_yyy_xzz_0,   \
                             ta2_xx_yyy_xzz_1,   \
                             ta2_xx_yyy_xzzz_0,  \
                             ta2_xx_yyy_xzzz_1,  \
                             ta2_xx_yyy_yyy_0,   \
                             ta2_xx_yyy_yyy_1,   \
                             ta2_xx_yyy_yyyy_0,  \
                             ta2_xx_yyy_yyyy_1,  \
                             ta2_xx_yyy_yyyz_0,  \
                             ta2_xx_yyy_yyyz_1,  \
                             ta2_xx_yyy_yyz_0,   \
                             ta2_xx_yyy_yyz_1,   \
                             ta2_xx_yyy_yyzz_0,  \
                             ta2_xx_yyy_yyzz_1,  \
                             ta2_xx_yyy_yzz_0,   \
                             ta2_xx_yyy_yzz_1,   \
                             ta2_xx_yyy_yzzz_0,  \
                             ta2_xx_yyy_yzzz_1,  \
                             ta2_xx_yyy_zzz_0,   \
                             ta2_xx_yyy_zzz_1,   \
                             ta2_xx_yyy_zzzz_0,  \
                             ta2_xx_yyy_zzzz_1,  \
                             ta2_xx_yyyy_xxxx_0, \
                             ta2_xx_yyyy_xxxy_0, \
                             ta2_xx_yyyy_xxxz_0, \
                             ta2_xx_yyyy_xxyy_0, \
                             ta2_xx_yyyy_xxyz_0, \
                             ta2_xx_yyyy_xxzz_0, \
                             ta2_xx_yyyy_xyyy_0, \
                             ta2_xx_yyyy_xyyz_0, \
                             ta2_xx_yyyy_xyzz_0, \
                             ta2_xx_yyyy_xzzz_0, \
                             ta2_xx_yyyy_yyyy_0, \
                             ta2_xx_yyyy_yyyz_0, \
                             ta2_xx_yyyy_yyzz_0, \
                             ta2_xx_yyyy_yzzz_0, \
                             ta2_xx_yyyy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyyy_xxxx_0[i] =
            3.0 * ta2_xx_yy_xxxx_0[i] * fe_0 - 3.0 * ta2_xx_yy_xxxx_1[i] * fe_0 + ta2_xx_yyy_xxxx_0[i] * pa_y[i] - ta2_xx_yyy_xxxx_1[i] * pc_y[i];

        ta2_xx_yyyy_xxxy_0[i] = 3.0 * ta2_xx_yy_xxxy_0[i] * fe_0 - 3.0 * ta2_xx_yy_xxxy_1[i] * fe_0 + ta2_xx_yyy_xxx_0[i] * fe_0 -
                                ta2_xx_yyy_xxx_1[i] * fe_0 + ta2_xx_yyy_xxxy_0[i] * pa_y[i] - ta2_xx_yyy_xxxy_1[i] * pc_y[i];

        ta2_xx_yyyy_xxxz_0[i] =
            3.0 * ta2_xx_yy_xxxz_0[i] * fe_0 - 3.0 * ta2_xx_yy_xxxz_1[i] * fe_0 + ta2_xx_yyy_xxxz_0[i] * pa_y[i] - ta2_xx_yyy_xxxz_1[i] * pc_y[i];

        ta2_xx_yyyy_xxyy_0[i] = 3.0 * ta2_xx_yy_xxyy_0[i] * fe_0 - 3.0 * ta2_xx_yy_xxyy_1[i] * fe_0 + 2.0 * ta2_xx_yyy_xxy_0[i] * fe_0 -
                                2.0 * ta2_xx_yyy_xxy_1[i] * fe_0 + ta2_xx_yyy_xxyy_0[i] * pa_y[i] - ta2_xx_yyy_xxyy_1[i] * pc_y[i];

        ta2_xx_yyyy_xxyz_0[i] = 3.0 * ta2_xx_yy_xxyz_0[i] * fe_0 - 3.0 * ta2_xx_yy_xxyz_1[i] * fe_0 + ta2_xx_yyy_xxz_0[i] * fe_0 -
                                ta2_xx_yyy_xxz_1[i] * fe_0 + ta2_xx_yyy_xxyz_0[i] * pa_y[i] - ta2_xx_yyy_xxyz_1[i] * pc_y[i];

        ta2_xx_yyyy_xxzz_0[i] =
            3.0 * ta2_xx_yy_xxzz_0[i] * fe_0 - 3.0 * ta2_xx_yy_xxzz_1[i] * fe_0 + ta2_xx_yyy_xxzz_0[i] * pa_y[i] - ta2_xx_yyy_xxzz_1[i] * pc_y[i];

        ta2_xx_yyyy_xyyy_0[i] = 3.0 * ta2_xx_yy_xyyy_0[i] * fe_0 - 3.0 * ta2_xx_yy_xyyy_1[i] * fe_0 + 3.0 * ta2_xx_yyy_xyy_0[i] * fe_0 -
                                3.0 * ta2_xx_yyy_xyy_1[i] * fe_0 + ta2_xx_yyy_xyyy_0[i] * pa_y[i] - ta2_xx_yyy_xyyy_1[i] * pc_y[i];

        ta2_xx_yyyy_xyyz_0[i] = 3.0 * ta2_xx_yy_xyyz_0[i] * fe_0 - 3.0 * ta2_xx_yy_xyyz_1[i] * fe_0 + 2.0 * ta2_xx_yyy_xyz_0[i] * fe_0 -
                                2.0 * ta2_xx_yyy_xyz_1[i] * fe_0 + ta2_xx_yyy_xyyz_0[i] * pa_y[i] - ta2_xx_yyy_xyyz_1[i] * pc_y[i];

        ta2_xx_yyyy_xyzz_0[i] = 3.0 * ta2_xx_yy_xyzz_0[i] * fe_0 - 3.0 * ta2_xx_yy_xyzz_1[i] * fe_0 + ta2_xx_yyy_xzz_0[i] * fe_0 -
                                ta2_xx_yyy_xzz_1[i] * fe_0 + ta2_xx_yyy_xyzz_0[i] * pa_y[i] - ta2_xx_yyy_xyzz_1[i] * pc_y[i];

        ta2_xx_yyyy_xzzz_0[i] =
            3.0 * ta2_xx_yy_xzzz_0[i] * fe_0 - 3.0 * ta2_xx_yy_xzzz_1[i] * fe_0 + ta2_xx_yyy_xzzz_0[i] * pa_y[i] - ta2_xx_yyy_xzzz_1[i] * pc_y[i];

        ta2_xx_yyyy_yyyy_0[i] = 3.0 * ta2_xx_yy_yyyy_0[i] * fe_0 - 3.0 * ta2_xx_yy_yyyy_1[i] * fe_0 + 4.0 * ta2_xx_yyy_yyy_0[i] * fe_0 -
                                4.0 * ta2_xx_yyy_yyy_1[i] * fe_0 + ta2_xx_yyy_yyyy_0[i] * pa_y[i] - ta2_xx_yyy_yyyy_1[i] * pc_y[i];

        ta2_xx_yyyy_yyyz_0[i] = 3.0 * ta2_xx_yy_yyyz_0[i] * fe_0 - 3.0 * ta2_xx_yy_yyyz_1[i] * fe_0 + 3.0 * ta2_xx_yyy_yyz_0[i] * fe_0 -
                                3.0 * ta2_xx_yyy_yyz_1[i] * fe_0 + ta2_xx_yyy_yyyz_0[i] * pa_y[i] - ta2_xx_yyy_yyyz_1[i] * pc_y[i];

        ta2_xx_yyyy_yyzz_0[i] = 3.0 * ta2_xx_yy_yyzz_0[i] * fe_0 - 3.0 * ta2_xx_yy_yyzz_1[i] * fe_0 + 2.0 * ta2_xx_yyy_yzz_0[i] * fe_0 -
                                2.0 * ta2_xx_yyy_yzz_1[i] * fe_0 + ta2_xx_yyy_yyzz_0[i] * pa_y[i] - ta2_xx_yyy_yyzz_1[i] * pc_y[i];

        ta2_xx_yyyy_yzzz_0[i] = 3.0 * ta2_xx_yy_yzzz_0[i] * fe_0 - 3.0 * ta2_xx_yy_yzzz_1[i] * fe_0 + ta2_xx_yyy_zzz_0[i] * fe_0 -
                                ta2_xx_yyy_zzz_1[i] * fe_0 + ta2_xx_yyy_yzzz_0[i] * pa_y[i] - ta2_xx_yyy_yzzz_1[i] * pc_y[i];

        ta2_xx_yyyy_zzzz_0[i] =
            3.0 * ta2_xx_yy_zzzz_0[i] * fe_0 - 3.0 * ta2_xx_yy_zzzz_1[i] * fe_0 + ta2_xx_yyy_zzzz_0[i] * pa_y[i] - ta2_xx_yyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 165-180 components of targeted buffer : GG

    auto ta2_xx_yyyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 165);

    auto ta2_xx_yyyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 166);

    auto ta2_xx_yyyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 167);

    auto ta2_xx_yyyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 168);

    auto ta2_xx_yyyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 169);

    auto ta2_xx_yyyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 170);

    auto ta2_xx_yyyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 171);

    auto ta2_xx_yyyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 172);

    auto ta2_xx_yyyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 173);

    auto ta2_xx_yyyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 174);

    auto ta2_xx_yyyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 175);

    auto ta2_xx_yyyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 176);

    auto ta2_xx_yyyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 177);

    auto ta2_xx_yyyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 178);

    auto ta2_xx_yyyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 179);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta2_xx_yyy_xxxx_0,  \
                             ta2_xx_yyy_xxxx_1,  \
                             ta2_xx_yyy_xxxy_0,  \
                             ta2_xx_yyy_xxxy_1,  \
                             ta2_xx_yyy_xxy_0,   \
                             ta2_xx_yyy_xxy_1,   \
                             ta2_xx_yyy_xxyy_0,  \
                             ta2_xx_yyy_xxyy_1,  \
                             ta2_xx_yyy_xxyz_0,  \
                             ta2_xx_yyy_xxyz_1,  \
                             ta2_xx_yyy_xyy_0,   \
                             ta2_xx_yyy_xyy_1,   \
                             ta2_xx_yyy_xyyy_0,  \
                             ta2_xx_yyy_xyyy_1,  \
                             ta2_xx_yyy_xyyz_0,  \
                             ta2_xx_yyy_xyyz_1,  \
                             ta2_xx_yyy_xyz_0,   \
                             ta2_xx_yyy_xyz_1,   \
                             ta2_xx_yyy_xyzz_0,  \
                             ta2_xx_yyy_xyzz_1,  \
                             ta2_xx_yyy_yyy_0,   \
                             ta2_xx_yyy_yyy_1,   \
                             ta2_xx_yyy_yyyy_0,  \
                             ta2_xx_yyy_yyyy_1,  \
                             ta2_xx_yyy_yyyz_0,  \
                             ta2_xx_yyy_yyyz_1,  \
                             ta2_xx_yyy_yyz_0,   \
                             ta2_xx_yyy_yyz_1,   \
                             ta2_xx_yyy_yyzz_0,  \
                             ta2_xx_yyy_yyzz_1,  \
                             ta2_xx_yyy_yzz_0,   \
                             ta2_xx_yyy_yzz_1,   \
                             ta2_xx_yyy_yzzz_0,  \
                             ta2_xx_yyy_yzzz_1,  \
                             ta2_xx_yyyz_xxxx_0, \
                             ta2_xx_yyyz_xxxy_0, \
                             ta2_xx_yyyz_xxxz_0, \
                             ta2_xx_yyyz_xxyy_0, \
                             ta2_xx_yyyz_xxyz_0, \
                             ta2_xx_yyyz_xxzz_0, \
                             ta2_xx_yyyz_xyyy_0, \
                             ta2_xx_yyyz_xyyz_0, \
                             ta2_xx_yyyz_xyzz_0, \
                             ta2_xx_yyyz_xzzz_0, \
                             ta2_xx_yyyz_yyyy_0, \
                             ta2_xx_yyyz_yyyz_0, \
                             ta2_xx_yyyz_yyzz_0, \
                             ta2_xx_yyyz_yzzz_0, \
                             ta2_xx_yyyz_zzzz_0, \
                             ta2_xx_yyz_xxxz_0,  \
                             ta2_xx_yyz_xxxz_1,  \
                             ta2_xx_yyz_xxzz_0,  \
                             ta2_xx_yyz_xxzz_1,  \
                             ta2_xx_yyz_xzzz_0,  \
                             ta2_xx_yyz_xzzz_1,  \
                             ta2_xx_yyz_zzzz_0,  \
                             ta2_xx_yyz_zzzz_1,  \
                             ta2_xx_yz_xxxz_0,   \
                             ta2_xx_yz_xxxz_1,   \
                             ta2_xx_yz_xxzz_0,   \
                             ta2_xx_yz_xxzz_1,   \
                             ta2_xx_yz_xzzz_0,   \
                             ta2_xx_yz_xzzz_1,   \
                             ta2_xx_yz_zzzz_0,   \
                             ta2_xx_yz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyyz_xxxx_0[i] = ta2_xx_yyy_xxxx_0[i] * pa_z[i] - ta2_xx_yyy_xxxx_1[i] * pc_z[i];

        ta2_xx_yyyz_xxxy_0[i] = ta2_xx_yyy_xxxy_0[i] * pa_z[i] - ta2_xx_yyy_xxxy_1[i] * pc_z[i];

        ta2_xx_yyyz_xxxz_0[i] =
            2.0 * ta2_xx_yz_xxxz_0[i] * fe_0 - 2.0 * ta2_xx_yz_xxxz_1[i] * fe_0 + ta2_xx_yyz_xxxz_0[i] * pa_y[i] - ta2_xx_yyz_xxxz_1[i] * pc_y[i];

        ta2_xx_yyyz_xxyy_0[i] = ta2_xx_yyy_xxyy_0[i] * pa_z[i] - ta2_xx_yyy_xxyy_1[i] * pc_z[i];

        ta2_xx_yyyz_xxyz_0[i] =
            ta2_xx_yyy_xxy_0[i] * fe_0 - ta2_xx_yyy_xxy_1[i] * fe_0 + ta2_xx_yyy_xxyz_0[i] * pa_z[i] - ta2_xx_yyy_xxyz_1[i] * pc_z[i];

        ta2_xx_yyyz_xxzz_0[i] =
            2.0 * ta2_xx_yz_xxzz_0[i] * fe_0 - 2.0 * ta2_xx_yz_xxzz_1[i] * fe_0 + ta2_xx_yyz_xxzz_0[i] * pa_y[i] - ta2_xx_yyz_xxzz_1[i] * pc_y[i];

        ta2_xx_yyyz_xyyy_0[i] = ta2_xx_yyy_xyyy_0[i] * pa_z[i] - ta2_xx_yyy_xyyy_1[i] * pc_z[i];

        ta2_xx_yyyz_xyyz_0[i] =
            ta2_xx_yyy_xyy_0[i] * fe_0 - ta2_xx_yyy_xyy_1[i] * fe_0 + ta2_xx_yyy_xyyz_0[i] * pa_z[i] - ta2_xx_yyy_xyyz_1[i] * pc_z[i];

        ta2_xx_yyyz_xyzz_0[i] =
            2.0 * ta2_xx_yyy_xyz_0[i] * fe_0 - 2.0 * ta2_xx_yyy_xyz_1[i] * fe_0 + ta2_xx_yyy_xyzz_0[i] * pa_z[i] - ta2_xx_yyy_xyzz_1[i] * pc_z[i];

        ta2_xx_yyyz_xzzz_0[i] =
            2.0 * ta2_xx_yz_xzzz_0[i] * fe_0 - 2.0 * ta2_xx_yz_xzzz_1[i] * fe_0 + ta2_xx_yyz_xzzz_0[i] * pa_y[i] - ta2_xx_yyz_xzzz_1[i] * pc_y[i];

        ta2_xx_yyyz_yyyy_0[i] = ta2_xx_yyy_yyyy_0[i] * pa_z[i] - ta2_xx_yyy_yyyy_1[i] * pc_z[i];

        ta2_xx_yyyz_yyyz_0[i] =
            ta2_xx_yyy_yyy_0[i] * fe_0 - ta2_xx_yyy_yyy_1[i] * fe_0 + ta2_xx_yyy_yyyz_0[i] * pa_z[i] - ta2_xx_yyy_yyyz_1[i] * pc_z[i];

        ta2_xx_yyyz_yyzz_0[i] =
            2.0 * ta2_xx_yyy_yyz_0[i] * fe_0 - 2.0 * ta2_xx_yyy_yyz_1[i] * fe_0 + ta2_xx_yyy_yyzz_0[i] * pa_z[i] - ta2_xx_yyy_yyzz_1[i] * pc_z[i];

        ta2_xx_yyyz_yzzz_0[i] =
            3.0 * ta2_xx_yyy_yzz_0[i] * fe_0 - 3.0 * ta2_xx_yyy_yzz_1[i] * fe_0 + ta2_xx_yyy_yzzz_0[i] * pa_z[i] - ta2_xx_yyy_yzzz_1[i] * pc_z[i];

        ta2_xx_yyyz_zzzz_0[i] =
            2.0 * ta2_xx_yz_zzzz_0[i] * fe_0 - 2.0 * ta2_xx_yz_zzzz_1[i] * fe_0 + ta2_xx_yyz_zzzz_0[i] * pa_y[i] - ta2_xx_yyz_zzzz_1[i] * pc_y[i];
    }

    // Set up 180-195 components of targeted buffer : GG

    auto ta2_xx_yyzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 180);

    auto ta2_xx_yyzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 181);

    auto ta2_xx_yyzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 182);

    auto ta2_xx_yyzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 183);

    auto ta2_xx_yyzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 184);

    auto ta2_xx_yyzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 185);

    auto ta2_xx_yyzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 186);

    auto ta2_xx_yyzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 187);

    auto ta2_xx_yyzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 188);

    auto ta2_xx_yyzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 189);

    auto ta2_xx_yyzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 190);

    auto ta2_xx_yyzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 191);

    auto ta2_xx_yyzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 192);

    auto ta2_xx_yyzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 193);

    auto ta2_xx_yyzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 194);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta2_xx_yy_xxxy_0,   \
                             ta2_xx_yy_xxxy_1,   \
                             ta2_xx_yy_xxyy_0,   \
                             ta2_xx_yy_xxyy_1,   \
                             ta2_xx_yy_xyyy_0,   \
                             ta2_xx_yy_xyyy_1,   \
                             ta2_xx_yy_yyyy_0,   \
                             ta2_xx_yy_yyyy_1,   \
                             ta2_xx_yyz_xxxy_0,  \
                             ta2_xx_yyz_xxxy_1,  \
                             ta2_xx_yyz_xxyy_0,  \
                             ta2_xx_yyz_xxyy_1,  \
                             ta2_xx_yyz_xyyy_0,  \
                             ta2_xx_yyz_xyyy_1,  \
                             ta2_xx_yyz_yyyy_0,  \
                             ta2_xx_yyz_yyyy_1,  \
                             ta2_xx_yyzz_xxxx_0, \
                             ta2_xx_yyzz_xxxy_0, \
                             ta2_xx_yyzz_xxxz_0, \
                             ta2_xx_yyzz_xxyy_0, \
                             ta2_xx_yyzz_xxyz_0, \
                             ta2_xx_yyzz_xxzz_0, \
                             ta2_xx_yyzz_xyyy_0, \
                             ta2_xx_yyzz_xyyz_0, \
                             ta2_xx_yyzz_xyzz_0, \
                             ta2_xx_yyzz_xzzz_0, \
                             ta2_xx_yyzz_yyyy_0, \
                             ta2_xx_yyzz_yyyz_0, \
                             ta2_xx_yyzz_yyzz_0, \
                             ta2_xx_yyzz_yzzz_0, \
                             ta2_xx_yyzz_zzzz_0, \
                             ta2_xx_yzz_xxxx_0,  \
                             ta2_xx_yzz_xxxx_1,  \
                             ta2_xx_yzz_xxxz_0,  \
                             ta2_xx_yzz_xxxz_1,  \
                             ta2_xx_yzz_xxyz_0,  \
                             ta2_xx_yzz_xxyz_1,  \
                             ta2_xx_yzz_xxz_0,   \
                             ta2_xx_yzz_xxz_1,   \
                             ta2_xx_yzz_xxzz_0,  \
                             ta2_xx_yzz_xxzz_1,  \
                             ta2_xx_yzz_xyyz_0,  \
                             ta2_xx_yzz_xyyz_1,  \
                             ta2_xx_yzz_xyz_0,   \
                             ta2_xx_yzz_xyz_1,   \
                             ta2_xx_yzz_xyzz_0,  \
                             ta2_xx_yzz_xyzz_1,  \
                             ta2_xx_yzz_xzz_0,   \
                             ta2_xx_yzz_xzz_1,   \
                             ta2_xx_yzz_xzzz_0,  \
                             ta2_xx_yzz_xzzz_1,  \
                             ta2_xx_yzz_yyyz_0,  \
                             ta2_xx_yzz_yyyz_1,  \
                             ta2_xx_yzz_yyz_0,   \
                             ta2_xx_yzz_yyz_1,   \
                             ta2_xx_yzz_yyzz_0,  \
                             ta2_xx_yzz_yyzz_1,  \
                             ta2_xx_yzz_yzz_0,   \
                             ta2_xx_yzz_yzz_1,   \
                             ta2_xx_yzz_yzzz_0,  \
                             ta2_xx_yzz_yzzz_1,  \
                             ta2_xx_yzz_zzz_0,   \
                             ta2_xx_yzz_zzz_1,   \
                             ta2_xx_yzz_zzzz_0,  \
                             ta2_xx_yzz_zzzz_1,  \
                             ta2_xx_zz_xxxx_0,   \
                             ta2_xx_zz_xxxx_1,   \
                             ta2_xx_zz_xxxz_0,   \
                             ta2_xx_zz_xxxz_1,   \
                             ta2_xx_zz_xxyz_0,   \
                             ta2_xx_zz_xxyz_1,   \
                             ta2_xx_zz_xxzz_0,   \
                             ta2_xx_zz_xxzz_1,   \
                             ta2_xx_zz_xyyz_0,   \
                             ta2_xx_zz_xyyz_1,   \
                             ta2_xx_zz_xyzz_0,   \
                             ta2_xx_zz_xyzz_1,   \
                             ta2_xx_zz_xzzz_0,   \
                             ta2_xx_zz_xzzz_1,   \
                             ta2_xx_zz_yyyz_0,   \
                             ta2_xx_zz_yyyz_1,   \
                             ta2_xx_zz_yyzz_0,   \
                             ta2_xx_zz_yyzz_1,   \
                             ta2_xx_zz_yzzz_0,   \
                             ta2_xx_zz_yzzz_1,   \
                             ta2_xx_zz_zzzz_0,   \
                             ta2_xx_zz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yyzz_xxxx_0[i] =
            ta2_xx_zz_xxxx_0[i] * fe_0 - ta2_xx_zz_xxxx_1[i] * fe_0 + ta2_xx_yzz_xxxx_0[i] * pa_y[i] - ta2_xx_yzz_xxxx_1[i] * pc_y[i];

        ta2_xx_yyzz_xxxy_0[i] =
            ta2_xx_yy_xxxy_0[i] * fe_0 - ta2_xx_yy_xxxy_1[i] * fe_0 + ta2_xx_yyz_xxxy_0[i] * pa_z[i] - ta2_xx_yyz_xxxy_1[i] * pc_z[i];

        ta2_xx_yyzz_xxxz_0[i] =
            ta2_xx_zz_xxxz_0[i] * fe_0 - ta2_xx_zz_xxxz_1[i] * fe_0 + ta2_xx_yzz_xxxz_0[i] * pa_y[i] - ta2_xx_yzz_xxxz_1[i] * pc_y[i];

        ta2_xx_yyzz_xxyy_0[i] =
            ta2_xx_yy_xxyy_0[i] * fe_0 - ta2_xx_yy_xxyy_1[i] * fe_0 + ta2_xx_yyz_xxyy_0[i] * pa_z[i] - ta2_xx_yyz_xxyy_1[i] * pc_z[i];

        ta2_xx_yyzz_xxyz_0[i] = ta2_xx_zz_xxyz_0[i] * fe_0 - ta2_xx_zz_xxyz_1[i] * fe_0 + ta2_xx_yzz_xxz_0[i] * fe_0 - ta2_xx_yzz_xxz_1[i] * fe_0 +
                                ta2_xx_yzz_xxyz_0[i] * pa_y[i] - ta2_xx_yzz_xxyz_1[i] * pc_y[i];

        ta2_xx_yyzz_xxzz_0[i] =
            ta2_xx_zz_xxzz_0[i] * fe_0 - ta2_xx_zz_xxzz_1[i] * fe_0 + ta2_xx_yzz_xxzz_0[i] * pa_y[i] - ta2_xx_yzz_xxzz_1[i] * pc_y[i];

        ta2_xx_yyzz_xyyy_0[i] =
            ta2_xx_yy_xyyy_0[i] * fe_0 - ta2_xx_yy_xyyy_1[i] * fe_0 + ta2_xx_yyz_xyyy_0[i] * pa_z[i] - ta2_xx_yyz_xyyy_1[i] * pc_z[i];

        ta2_xx_yyzz_xyyz_0[i] = ta2_xx_zz_xyyz_0[i] * fe_0 - ta2_xx_zz_xyyz_1[i] * fe_0 + 2.0 * ta2_xx_yzz_xyz_0[i] * fe_0 -
                                2.0 * ta2_xx_yzz_xyz_1[i] * fe_0 + ta2_xx_yzz_xyyz_0[i] * pa_y[i] - ta2_xx_yzz_xyyz_1[i] * pc_y[i];

        ta2_xx_yyzz_xyzz_0[i] = ta2_xx_zz_xyzz_0[i] * fe_0 - ta2_xx_zz_xyzz_1[i] * fe_0 + ta2_xx_yzz_xzz_0[i] * fe_0 - ta2_xx_yzz_xzz_1[i] * fe_0 +
                                ta2_xx_yzz_xyzz_0[i] * pa_y[i] - ta2_xx_yzz_xyzz_1[i] * pc_y[i];

        ta2_xx_yyzz_xzzz_0[i] =
            ta2_xx_zz_xzzz_0[i] * fe_0 - ta2_xx_zz_xzzz_1[i] * fe_0 + ta2_xx_yzz_xzzz_0[i] * pa_y[i] - ta2_xx_yzz_xzzz_1[i] * pc_y[i];

        ta2_xx_yyzz_yyyy_0[i] =
            ta2_xx_yy_yyyy_0[i] * fe_0 - ta2_xx_yy_yyyy_1[i] * fe_0 + ta2_xx_yyz_yyyy_0[i] * pa_z[i] - ta2_xx_yyz_yyyy_1[i] * pc_z[i];

        ta2_xx_yyzz_yyyz_0[i] = ta2_xx_zz_yyyz_0[i] * fe_0 - ta2_xx_zz_yyyz_1[i] * fe_0 + 3.0 * ta2_xx_yzz_yyz_0[i] * fe_0 -
                                3.0 * ta2_xx_yzz_yyz_1[i] * fe_0 + ta2_xx_yzz_yyyz_0[i] * pa_y[i] - ta2_xx_yzz_yyyz_1[i] * pc_y[i];

        ta2_xx_yyzz_yyzz_0[i] = ta2_xx_zz_yyzz_0[i] * fe_0 - ta2_xx_zz_yyzz_1[i] * fe_0 + 2.0 * ta2_xx_yzz_yzz_0[i] * fe_0 -
                                2.0 * ta2_xx_yzz_yzz_1[i] * fe_0 + ta2_xx_yzz_yyzz_0[i] * pa_y[i] - ta2_xx_yzz_yyzz_1[i] * pc_y[i];

        ta2_xx_yyzz_yzzz_0[i] = ta2_xx_zz_yzzz_0[i] * fe_0 - ta2_xx_zz_yzzz_1[i] * fe_0 + ta2_xx_yzz_zzz_0[i] * fe_0 - ta2_xx_yzz_zzz_1[i] * fe_0 +
                                ta2_xx_yzz_yzzz_0[i] * pa_y[i] - ta2_xx_yzz_yzzz_1[i] * pc_y[i];

        ta2_xx_yyzz_zzzz_0[i] =
            ta2_xx_zz_zzzz_0[i] * fe_0 - ta2_xx_zz_zzzz_1[i] * fe_0 + ta2_xx_yzz_zzzz_0[i] * pa_y[i] - ta2_xx_yzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 195-210 components of targeted buffer : GG

    auto ta2_xx_yzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 195);

    auto ta2_xx_yzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 196);

    auto ta2_xx_yzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 197);

    auto ta2_xx_yzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 198);

    auto ta2_xx_yzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 199);

    auto ta2_xx_yzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 200);

    auto ta2_xx_yzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 201);

    auto ta2_xx_yzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 202);

    auto ta2_xx_yzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 203);

    auto ta2_xx_yzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 204);

    auto ta2_xx_yzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 205);

    auto ta2_xx_yzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 206);

    auto ta2_xx_yzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 207);

    auto ta2_xx_yzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 208);

    auto ta2_xx_yzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 209);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta2_xx_yzzz_xxxx_0, \
                             ta2_xx_yzzz_xxxy_0, \
                             ta2_xx_yzzz_xxxz_0, \
                             ta2_xx_yzzz_xxyy_0, \
                             ta2_xx_yzzz_xxyz_0, \
                             ta2_xx_yzzz_xxzz_0, \
                             ta2_xx_yzzz_xyyy_0, \
                             ta2_xx_yzzz_xyyz_0, \
                             ta2_xx_yzzz_xyzz_0, \
                             ta2_xx_yzzz_xzzz_0, \
                             ta2_xx_yzzz_yyyy_0, \
                             ta2_xx_yzzz_yyyz_0, \
                             ta2_xx_yzzz_yyzz_0, \
                             ta2_xx_yzzz_yzzz_0, \
                             ta2_xx_yzzz_zzzz_0, \
                             ta2_xx_zzz_xxx_0,   \
                             ta2_xx_zzz_xxx_1,   \
                             ta2_xx_zzz_xxxx_0,  \
                             ta2_xx_zzz_xxxx_1,  \
                             ta2_xx_zzz_xxxy_0,  \
                             ta2_xx_zzz_xxxy_1,  \
                             ta2_xx_zzz_xxxz_0,  \
                             ta2_xx_zzz_xxxz_1,  \
                             ta2_xx_zzz_xxy_0,   \
                             ta2_xx_zzz_xxy_1,   \
                             ta2_xx_zzz_xxyy_0,  \
                             ta2_xx_zzz_xxyy_1,  \
                             ta2_xx_zzz_xxyz_0,  \
                             ta2_xx_zzz_xxyz_1,  \
                             ta2_xx_zzz_xxz_0,   \
                             ta2_xx_zzz_xxz_1,   \
                             ta2_xx_zzz_xxzz_0,  \
                             ta2_xx_zzz_xxzz_1,  \
                             ta2_xx_zzz_xyy_0,   \
                             ta2_xx_zzz_xyy_1,   \
                             ta2_xx_zzz_xyyy_0,  \
                             ta2_xx_zzz_xyyy_1,  \
                             ta2_xx_zzz_xyyz_0,  \
                             ta2_xx_zzz_xyyz_1,  \
                             ta2_xx_zzz_xyz_0,   \
                             ta2_xx_zzz_xyz_1,   \
                             ta2_xx_zzz_xyzz_0,  \
                             ta2_xx_zzz_xyzz_1,  \
                             ta2_xx_zzz_xzz_0,   \
                             ta2_xx_zzz_xzz_1,   \
                             ta2_xx_zzz_xzzz_0,  \
                             ta2_xx_zzz_xzzz_1,  \
                             ta2_xx_zzz_yyy_0,   \
                             ta2_xx_zzz_yyy_1,   \
                             ta2_xx_zzz_yyyy_0,  \
                             ta2_xx_zzz_yyyy_1,  \
                             ta2_xx_zzz_yyyz_0,  \
                             ta2_xx_zzz_yyyz_1,  \
                             ta2_xx_zzz_yyz_0,   \
                             ta2_xx_zzz_yyz_1,   \
                             ta2_xx_zzz_yyzz_0,  \
                             ta2_xx_zzz_yyzz_1,  \
                             ta2_xx_zzz_yzz_0,   \
                             ta2_xx_zzz_yzz_1,   \
                             ta2_xx_zzz_yzzz_0,  \
                             ta2_xx_zzz_yzzz_1,  \
                             ta2_xx_zzz_zzz_0,   \
                             ta2_xx_zzz_zzz_1,   \
                             ta2_xx_zzz_zzzz_0,  \
                             ta2_xx_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_yzzz_xxxx_0[i] = ta2_xx_zzz_xxxx_0[i] * pa_y[i] - ta2_xx_zzz_xxxx_1[i] * pc_y[i];

        ta2_xx_yzzz_xxxy_0[i] =
            ta2_xx_zzz_xxx_0[i] * fe_0 - ta2_xx_zzz_xxx_1[i] * fe_0 + ta2_xx_zzz_xxxy_0[i] * pa_y[i] - ta2_xx_zzz_xxxy_1[i] * pc_y[i];

        ta2_xx_yzzz_xxxz_0[i] = ta2_xx_zzz_xxxz_0[i] * pa_y[i] - ta2_xx_zzz_xxxz_1[i] * pc_y[i];

        ta2_xx_yzzz_xxyy_0[i] =
            2.0 * ta2_xx_zzz_xxy_0[i] * fe_0 - 2.0 * ta2_xx_zzz_xxy_1[i] * fe_0 + ta2_xx_zzz_xxyy_0[i] * pa_y[i] - ta2_xx_zzz_xxyy_1[i] * pc_y[i];

        ta2_xx_yzzz_xxyz_0[i] =
            ta2_xx_zzz_xxz_0[i] * fe_0 - ta2_xx_zzz_xxz_1[i] * fe_0 + ta2_xx_zzz_xxyz_0[i] * pa_y[i] - ta2_xx_zzz_xxyz_1[i] * pc_y[i];

        ta2_xx_yzzz_xxzz_0[i] = ta2_xx_zzz_xxzz_0[i] * pa_y[i] - ta2_xx_zzz_xxzz_1[i] * pc_y[i];

        ta2_xx_yzzz_xyyy_0[i] =
            3.0 * ta2_xx_zzz_xyy_0[i] * fe_0 - 3.0 * ta2_xx_zzz_xyy_1[i] * fe_0 + ta2_xx_zzz_xyyy_0[i] * pa_y[i] - ta2_xx_zzz_xyyy_1[i] * pc_y[i];

        ta2_xx_yzzz_xyyz_0[i] =
            2.0 * ta2_xx_zzz_xyz_0[i] * fe_0 - 2.0 * ta2_xx_zzz_xyz_1[i] * fe_0 + ta2_xx_zzz_xyyz_0[i] * pa_y[i] - ta2_xx_zzz_xyyz_1[i] * pc_y[i];

        ta2_xx_yzzz_xyzz_0[i] =
            ta2_xx_zzz_xzz_0[i] * fe_0 - ta2_xx_zzz_xzz_1[i] * fe_0 + ta2_xx_zzz_xyzz_0[i] * pa_y[i] - ta2_xx_zzz_xyzz_1[i] * pc_y[i];

        ta2_xx_yzzz_xzzz_0[i] = ta2_xx_zzz_xzzz_0[i] * pa_y[i] - ta2_xx_zzz_xzzz_1[i] * pc_y[i];

        ta2_xx_yzzz_yyyy_0[i] =
            4.0 * ta2_xx_zzz_yyy_0[i] * fe_0 - 4.0 * ta2_xx_zzz_yyy_1[i] * fe_0 + ta2_xx_zzz_yyyy_0[i] * pa_y[i] - ta2_xx_zzz_yyyy_1[i] * pc_y[i];

        ta2_xx_yzzz_yyyz_0[i] =
            3.0 * ta2_xx_zzz_yyz_0[i] * fe_0 - 3.0 * ta2_xx_zzz_yyz_1[i] * fe_0 + ta2_xx_zzz_yyyz_0[i] * pa_y[i] - ta2_xx_zzz_yyyz_1[i] * pc_y[i];

        ta2_xx_yzzz_yyzz_0[i] =
            2.0 * ta2_xx_zzz_yzz_0[i] * fe_0 - 2.0 * ta2_xx_zzz_yzz_1[i] * fe_0 + ta2_xx_zzz_yyzz_0[i] * pa_y[i] - ta2_xx_zzz_yyzz_1[i] * pc_y[i];

        ta2_xx_yzzz_yzzz_0[i] =
            ta2_xx_zzz_zzz_0[i] * fe_0 - ta2_xx_zzz_zzz_1[i] * fe_0 + ta2_xx_zzz_yzzz_0[i] * pa_y[i] - ta2_xx_zzz_yzzz_1[i] * pc_y[i];

        ta2_xx_yzzz_zzzz_0[i] = ta2_xx_zzz_zzzz_0[i] * pa_y[i] - ta2_xx_zzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 210-225 components of targeted buffer : GG

    auto ta2_xx_zzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 210);

    auto ta2_xx_zzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 211);

    auto ta2_xx_zzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 212);

    auto ta2_xx_zzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 213);

    auto ta2_xx_zzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 214);

    auto ta2_xx_zzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 215);

    auto ta2_xx_zzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 216);

    auto ta2_xx_zzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 217);

    auto ta2_xx_zzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 218);

    auto ta2_xx_zzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 219);

    auto ta2_xx_zzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 220);

    auto ta2_xx_zzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 221);

    auto ta2_xx_zzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 222);

    auto ta2_xx_zzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 223);

    auto ta2_xx_zzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 224);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta2_xx_zz_xxxx_0,   \
                             ta2_xx_zz_xxxx_1,   \
                             ta2_xx_zz_xxxy_0,   \
                             ta2_xx_zz_xxxy_1,   \
                             ta2_xx_zz_xxxz_0,   \
                             ta2_xx_zz_xxxz_1,   \
                             ta2_xx_zz_xxyy_0,   \
                             ta2_xx_zz_xxyy_1,   \
                             ta2_xx_zz_xxyz_0,   \
                             ta2_xx_zz_xxyz_1,   \
                             ta2_xx_zz_xxzz_0,   \
                             ta2_xx_zz_xxzz_1,   \
                             ta2_xx_zz_xyyy_0,   \
                             ta2_xx_zz_xyyy_1,   \
                             ta2_xx_zz_xyyz_0,   \
                             ta2_xx_zz_xyyz_1,   \
                             ta2_xx_zz_xyzz_0,   \
                             ta2_xx_zz_xyzz_1,   \
                             ta2_xx_zz_xzzz_0,   \
                             ta2_xx_zz_xzzz_1,   \
                             ta2_xx_zz_yyyy_0,   \
                             ta2_xx_zz_yyyy_1,   \
                             ta2_xx_zz_yyyz_0,   \
                             ta2_xx_zz_yyyz_1,   \
                             ta2_xx_zz_yyzz_0,   \
                             ta2_xx_zz_yyzz_1,   \
                             ta2_xx_zz_yzzz_0,   \
                             ta2_xx_zz_yzzz_1,   \
                             ta2_xx_zz_zzzz_0,   \
                             ta2_xx_zz_zzzz_1,   \
                             ta2_xx_zzz_xxx_0,   \
                             ta2_xx_zzz_xxx_1,   \
                             ta2_xx_zzz_xxxx_0,  \
                             ta2_xx_zzz_xxxx_1,  \
                             ta2_xx_zzz_xxxy_0,  \
                             ta2_xx_zzz_xxxy_1,  \
                             ta2_xx_zzz_xxxz_0,  \
                             ta2_xx_zzz_xxxz_1,  \
                             ta2_xx_zzz_xxy_0,   \
                             ta2_xx_zzz_xxy_1,   \
                             ta2_xx_zzz_xxyy_0,  \
                             ta2_xx_zzz_xxyy_1,  \
                             ta2_xx_zzz_xxyz_0,  \
                             ta2_xx_zzz_xxyz_1,  \
                             ta2_xx_zzz_xxz_0,   \
                             ta2_xx_zzz_xxz_1,   \
                             ta2_xx_zzz_xxzz_0,  \
                             ta2_xx_zzz_xxzz_1,  \
                             ta2_xx_zzz_xyy_0,   \
                             ta2_xx_zzz_xyy_1,   \
                             ta2_xx_zzz_xyyy_0,  \
                             ta2_xx_zzz_xyyy_1,  \
                             ta2_xx_zzz_xyyz_0,  \
                             ta2_xx_zzz_xyyz_1,  \
                             ta2_xx_zzz_xyz_0,   \
                             ta2_xx_zzz_xyz_1,   \
                             ta2_xx_zzz_xyzz_0,  \
                             ta2_xx_zzz_xyzz_1,  \
                             ta2_xx_zzz_xzz_0,   \
                             ta2_xx_zzz_xzz_1,   \
                             ta2_xx_zzz_xzzz_0,  \
                             ta2_xx_zzz_xzzz_1,  \
                             ta2_xx_zzz_yyy_0,   \
                             ta2_xx_zzz_yyy_1,   \
                             ta2_xx_zzz_yyyy_0,  \
                             ta2_xx_zzz_yyyy_1,  \
                             ta2_xx_zzz_yyyz_0,  \
                             ta2_xx_zzz_yyyz_1,  \
                             ta2_xx_zzz_yyz_0,   \
                             ta2_xx_zzz_yyz_1,   \
                             ta2_xx_zzz_yyzz_0,  \
                             ta2_xx_zzz_yyzz_1,  \
                             ta2_xx_zzz_yzz_0,   \
                             ta2_xx_zzz_yzz_1,   \
                             ta2_xx_zzz_yzzz_0,  \
                             ta2_xx_zzz_yzzz_1,  \
                             ta2_xx_zzz_zzz_0,   \
                             ta2_xx_zzz_zzz_1,   \
                             ta2_xx_zzz_zzzz_0,  \
                             ta2_xx_zzz_zzzz_1,  \
                             ta2_xx_zzzz_xxxx_0, \
                             ta2_xx_zzzz_xxxy_0, \
                             ta2_xx_zzzz_xxxz_0, \
                             ta2_xx_zzzz_xxyy_0, \
                             ta2_xx_zzzz_xxyz_0, \
                             ta2_xx_zzzz_xxzz_0, \
                             ta2_xx_zzzz_xyyy_0, \
                             ta2_xx_zzzz_xyyz_0, \
                             ta2_xx_zzzz_xyzz_0, \
                             ta2_xx_zzzz_xzzz_0, \
                             ta2_xx_zzzz_yyyy_0, \
                             ta2_xx_zzzz_yyyz_0, \
                             ta2_xx_zzzz_yyzz_0, \
                             ta2_xx_zzzz_yzzz_0, \
                             ta2_xx_zzzz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xx_zzzz_xxxx_0[i] =
            3.0 * ta2_xx_zz_xxxx_0[i] * fe_0 - 3.0 * ta2_xx_zz_xxxx_1[i] * fe_0 + ta2_xx_zzz_xxxx_0[i] * pa_z[i] - ta2_xx_zzz_xxxx_1[i] * pc_z[i];

        ta2_xx_zzzz_xxxy_0[i] =
            3.0 * ta2_xx_zz_xxxy_0[i] * fe_0 - 3.0 * ta2_xx_zz_xxxy_1[i] * fe_0 + ta2_xx_zzz_xxxy_0[i] * pa_z[i] - ta2_xx_zzz_xxxy_1[i] * pc_z[i];

        ta2_xx_zzzz_xxxz_0[i] = 3.0 * ta2_xx_zz_xxxz_0[i] * fe_0 - 3.0 * ta2_xx_zz_xxxz_1[i] * fe_0 + ta2_xx_zzz_xxx_0[i] * fe_0 -
                                ta2_xx_zzz_xxx_1[i] * fe_0 + ta2_xx_zzz_xxxz_0[i] * pa_z[i] - ta2_xx_zzz_xxxz_1[i] * pc_z[i];

        ta2_xx_zzzz_xxyy_0[i] =
            3.0 * ta2_xx_zz_xxyy_0[i] * fe_0 - 3.0 * ta2_xx_zz_xxyy_1[i] * fe_0 + ta2_xx_zzz_xxyy_0[i] * pa_z[i] - ta2_xx_zzz_xxyy_1[i] * pc_z[i];

        ta2_xx_zzzz_xxyz_0[i] = 3.0 * ta2_xx_zz_xxyz_0[i] * fe_0 - 3.0 * ta2_xx_zz_xxyz_1[i] * fe_0 + ta2_xx_zzz_xxy_0[i] * fe_0 -
                                ta2_xx_zzz_xxy_1[i] * fe_0 + ta2_xx_zzz_xxyz_0[i] * pa_z[i] - ta2_xx_zzz_xxyz_1[i] * pc_z[i];

        ta2_xx_zzzz_xxzz_0[i] = 3.0 * ta2_xx_zz_xxzz_0[i] * fe_0 - 3.0 * ta2_xx_zz_xxzz_1[i] * fe_0 + 2.0 * ta2_xx_zzz_xxz_0[i] * fe_0 -
                                2.0 * ta2_xx_zzz_xxz_1[i] * fe_0 + ta2_xx_zzz_xxzz_0[i] * pa_z[i] - ta2_xx_zzz_xxzz_1[i] * pc_z[i];

        ta2_xx_zzzz_xyyy_0[i] =
            3.0 * ta2_xx_zz_xyyy_0[i] * fe_0 - 3.0 * ta2_xx_zz_xyyy_1[i] * fe_0 + ta2_xx_zzz_xyyy_0[i] * pa_z[i] - ta2_xx_zzz_xyyy_1[i] * pc_z[i];

        ta2_xx_zzzz_xyyz_0[i] = 3.0 * ta2_xx_zz_xyyz_0[i] * fe_0 - 3.0 * ta2_xx_zz_xyyz_1[i] * fe_0 + ta2_xx_zzz_xyy_0[i] * fe_0 -
                                ta2_xx_zzz_xyy_1[i] * fe_0 + ta2_xx_zzz_xyyz_0[i] * pa_z[i] - ta2_xx_zzz_xyyz_1[i] * pc_z[i];

        ta2_xx_zzzz_xyzz_0[i] = 3.0 * ta2_xx_zz_xyzz_0[i] * fe_0 - 3.0 * ta2_xx_zz_xyzz_1[i] * fe_0 + 2.0 * ta2_xx_zzz_xyz_0[i] * fe_0 -
                                2.0 * ta2_xx_zzz_xyz_1[i] * fe_0 + ta2_xx_zzz_xyzz_0[i] * pa_z[i] - ta2_xx_zzz_xyzz_1[i] * pc_z[i];

        ta2_xx_zzzz_xzzz_0[i] = 3.0 * ta2_xx_zz_xzzz_0[i] * fe_0 - 3.0 * ta2_xx_zz_xzzz_1[i] * fe_0 + 3.0 * ta2_xx_zzz_xzz_0[i] * fe_0 -
                                3.0 * ta2_xx_zzz_xzz_1[i] * fe_0 + ta2_xx_zzz_xzzz_0[i] * pa_z[i] - ta2_xx_zzz_xzzz_1[i] * pc_z[i];

        ta2_xx_zzzz_yyyy_0[i] =
            3.0 * ta2_xx_zz_yyyy_0[i] * fe_0 - 3.0 * ta2_xx_zz_yyyy_1[i] * fe_0 + ta2_xx_zzz_yyyy_0[i] * pa_z[i] - ta2_xx_zzz_yyyy_1[i] * pc_z[i];

        ta2_xx_zzzz_yyyz_0[i] = 3.0 * ta2_xx_zz_yyyz_0[i] * fe_0 - 3.0 * ta2_xx_zz_yyyz_1[i] * fe_0 + ta2_xx_zzz_yyy_0[i] * fe_0 -
                                ta2_xx_zzz_yyy_1[i] * fe_0 + ta2_xx_zzz_yyyz_0[i] * pa_z[i] - ta2_xx_zzz_yyyz_1[i] * pc_z[i];

        ta2_xx_zzzz_yyzz_0[i] = 3.0 * ta2_xx_zz_yyzz_0[i] * fe_0 - 3.0 * ta2_xx_zz_yyzz_1[i] * fe_0 + 2.0 * ta2_xx_zzz_yyz_0[i] * fe_0 -
                                2.0 * ta2_xx_zzz_yyz_1[i] * fe_0 + ta2_xx_zzz_yyzz_0[i] * pa_z[i] - ta2_xx_zzz_yyzz_1[i] * pc_z[i];

        ta2_xx_zzzz_yzzz_0[i] = 3.0 * ta2_xx_zz_yzzz_0[i] * fe_0 - 3.0 * ta2_xx_zz_yzzz_1[i] * fe_0 + 3.0 * ta2_xx_zzz_yzz_0[i] * fe_0 -
                                3.0 * ta2_xx_zzz_yzz_1[i] * fe_0 + ta2_xx_zzz_yzzz_0[i] * pa_z[i] - ta2_xx_zzz_yzzz_1[i] * pc_z[i];

        ta2_xx_zzzz_zzzz_0[i] = 3.0 * ta2_xx_zz_zzzz_0[i] * fe_0 - 3.0 * ta2_xx_zz_zzzz_1[i] * fe_0 + 4.0 * ta2_xx_zzz_zzz_0[i] * fe_0 -
                                4.0 * ta2_xx_zzz_zzz_1[i] * fe_0 + ta2_xx_zzz_zzzz_0[i] * pa_z[i] - ta2_xx_zzz_zzzz_1[i] * pc_z[i];
    }

    // Set up 225-240 components of targeted buffer : GG

    auto ta2_xy_xxxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 225);

    auto ta2_xy_xxxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 226);

    auto ta2_xy_xxxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 227);

    auto ta2_xy_xxxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 228);

    auto ta2_xy_xxxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 229);

    auto ta2_xy_xxxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 230);

    auto ta2_xy_xxxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 231);

    auto ta2_xy_xxxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 232);

    auto ta2_xy_xxxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 233);

    auto ta2_xy_xxxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 234);

    auto ta2_xy_xxxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 235);

    auto ta2_xy_xxxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 236);

    auto ta2_xy_xxxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 237);

    auto ta2_xy_xxxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 238);

    auto ta2_xy_xxxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 239);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_xxx_xxxx_1,   \
                             ta1_y_xxx_xxxy_1,   \
                             ta1_y_xxx_xxxz_1,   \
                             ta1_y_xxx_xxyy_1,   \
                             ta1_y_xxx_xxyz_1,   \
                             ta1_y_xxx_xxzz_1,   \
                             ta1_y_xxx_xyyy_1,   \
                             ta1_y_xxx_xyyz_1,   \
                             ta1_y_xxx_xyzz_1,   \
                             ta1_y_xxx_xzzz_1,   \
                             ta1_y_xxx_yyyy_1,   \
                             ta1_y_xxx_yyyz_1,   \
                             ta1_y_xxx_yyzz_1,   \
                             ta1_y_xxx_yzzz_1,   \
                             ta1_y_xxx_zzzz_1,   \
                             ta2_xy_xx_xxxx_0,   \
                             ta2_xy_xx_xxxx_1,   \
                             ta2_xy_xx_xxxy_0,   \
                             ta2_xy_xx_xxxy_1,   \
                             ta2_xy_xx_xxxz_0,   \
                             ta2_xy_xx_xxxz_1,   \
                             ta2_xy_xx_xxyy_0,   \
                             ta2_xy_xx_xxyy_1,   \
                             ta2_xy_xx_xxyz_0,   \
                             ta2_xy_xx_xxyz_1,   \
                             ta2_xy_xx_xxzz_0,   \
                             ta2_xy_xx_xxzz_1,   \
                             ta2_xy_xx_xyyy_0,   \
                             ta2_xy_xx_xyyy_1,   \
                             ta2_xy_xx_xyyz_0,   \
                             ta2_xy_xx_xyyz_1,   \
                             ta2_xy_xx_xyzz_0,   \
                             ta2_xy_xx_xyzz_1,   \
                             ta2_xy_xx_xzzz_0,   \
                             ta2_xy_xx_xzzz_1,   \
                             ta2_xy_xx_yyyy_0,   \
                             ta2_xy_xx_yyyy_1,   \
                             ta2_xy_xx_yyyz_0,   \
                             ta2_xy_xx_yyyz_1,   \
                             ta2_xy_xx_yyzz_0,   \
                             ta2_xy_xx_yyzz_1,   \
                             ta2_xy_xx_yzzz_0,   \
                             ta2_xy_xx_yzzz_1,   \
                             ta2_xy_xx_zzzz_0,   \
                             ta2_xy_xx_zzzz_1,   \
                             ta2_xy_xxx_xxx_0,   \
                             ta2_xy_xxx_xxx_1,   \
                             ta2_xy_xxx_xxxx_0,  \
                             ta2_xy_xxx_xxxx_1,  \
                             ta2_xy_xxx_xxxy_0,  \
                             ta2_xy_xxx_xxxy_1,  \
                             ta2_xy_xxx_xxxz_0,  \
                             ta2_xy_xxx_xxxz_1,  \
                             ta2_xy_xxx_xxy_0,   \
                             ta2_xy_xxx_xxy_1,   \
                             ta2_xy_xxx_xxyy_0,  \
                             ta2_xy_xxx_xxyy_1,  \
                             ta2_xy_xxx_xxyz_0,  \
                             ta2_xy_xxx_xxyz_1,  \
                             ta2_xy_xxx_xxz_0,   \
                             ta2_xy_xxx_xxz_1,   \
                             ta2_xy_xxx_xxzz_0,  \
                             ta2_xy_xxx_xxzz_1,  \
                             ta2_xy_xxx_xyy_0,   \
                             ta2_xy_xxx_xyy_1,   \
                             ta2_xy_xxx_xyyy_0,  \
                             ta2_xy_xxx_xyyy_1,  \
                             ta2_xy_xxx_xyyz_0,  \
                             ta2_xy_xxx_xyyz_1,  \
                             ta2_xy_xxx_xyz_0,   \
                             ta2_xy_xxx_xyz_1,   \
                             ta2_xy_xxx_xyzz_0,  \
                             ta2_xy_xxx_xyzz_1,  \
                             ta2_xy_xxx_xzz_0,   \
                             ta2_xy_xxx_xzz_1,   \
                             ta2_xy_xxx_xzzz_0,  \
                             ta2_xy_xxx_xzzz_1,  \
                             ta2_xy_xxx_yyy_0,   \
                             ta2_xy_xxx_yyy_1,   \
                             ta2_xy_xxx_yyyy_0,  \
                             ta2_xy_xxx_yyyy_1,  \
                             ta2_xy_xxx_yyyz_0,  \
                             ta2_xy_xxx_yyyz_1,  \
                             ta2_xy_xxx_yyz_0,   \
                             ta2_xy_xxx_yyz_1,   \
                             ta2_xy_xxx_yyzz_0,  \
                             ta2_xy_xxx_yyzz_1,  \
                             ta2_xy_xxx_yzz_0,   \
                             ta2_xy_xxx_yzz_1,   \
                             ta2_xy_xxx_yzzz_0,  \
                             ta2_xy_xxx_yzzz_1,  \
                             ta2_xy_xxx_zzz_0,   \
                             ta2_xy_xxx_zzz_1,   \
                             ta2_xy_xxx_zzzz_0,  \
                             ta2_xy_xxx_zzzz_1,  \
                             ta2_xy_xxxx_xxxx_0, \
                             ta2_xy_xxxx_xxxy_0, \
                             ta2_xy_xxxx_xxxz_0, \
                             ta2_xy_xxxx_xxyy_0, \
                             ta2_xy_xxxx_xxyz_0, \
                             ta2_xy_xxxx_xxzz_0, \
                             ta2_xy_xxxx_xyyy_0, \
                             ta2_xy_xxxx_xyyz_0, \
                             ta2_xy_xxxx_xyzz_0, \
                             ta2_xy_xxxx_xzzz_0, \
                             ta2_xy_xxxx_yyyy_0, \
                             ta2_xy_xxxx_yyyz_0, \
                             ta2_xy_xxxx_yyzz_0, \
                             ta2_xy_xxxx_yzzz_0, \
                             ta2_xy_xxxx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxxx_xxxx_0[i] = 3.0 * ta2_xy_xx_xxxx_0[i] * fe_0 - 3.0 * ta2_xy_xx_xxxx_1[i] * fe_0 + 4.0 * ta2_xy_xxx_xxx_0[i] * fe_0 -
                                4.0 * ta2_xy_xxx_xxx_1[i] * fe_0 + ta1_y_xxx_xxxx_1[i] + ta2_xy_xxx_xxxx_0[i] * pa_x[i] -
                                ta2_xy_xxx_xxxx_1[i] * pc_x[i];

        ta2_xy_xxxx_xxxy_0[i] = 3.0 * ta2_xy_xx_xxxy_0[i] * fe_0 - 3.0 * ta2_xy_xx_xxxy_1[i] * fe_0 + 3.0 * ta2_xy_xxx_xxy_0[i] * fe_0 -
                                3.0 * ta2_xy_xxx_xxy_1[i] * fe_0 + ta1_y_xxx_xxxy_1[i] + ta2_xy_xxx_xxxy_0[i] * pa_x[i] -
                                ta2_xy_xxx_xxxy_1[i] * pc_x[i];

        ta2_xy_xxxx_xxxz_0[i] = 3.0 * ta2_xy_xx_xxxz_0[i] * fe_0 - 3.0 * ta2_xy_xx_xxxz_1[i] * fe_0 + 3.0 * ta2_xy_xxx_xxz_0[i] * fe_0 -
                                3.0 * ta2_xy_xxx_xxz_1[i] * fe_0 + ta1_y_xxx_xxxz_1[i] + ta2_xy_xxx_xxxz_0[i] * pa_x[i] -
                                ta2_xy_xxx_xxxz_1[i] * pc_x[i];

        ta2_xy_xxxx_xxyy_0[i] = 3.0 * ta2_xy_xx_xxyy_0[i] * fe_0 - 3.0 * ta2_xy_xx_xxyy_1[i] * fe_0 + 2.0 * ta2_xy_xxx_xyy_0[i] * fe_0 -
                                2.0 * ta2_xy_xxx_xyy_1[i] * fe_0 + ta1_y_xxx_xxyy_1[i] + ta2_xy_xxx_xxyy_0[i] * pa_x[i] -
                                ta2_xy_xxx_xxyy_1[i] * pc_x[i];

        ta2_xy_xxxx_xxyz_0[i] = 3.0 * ta2_xy_xx_xxyz_0[i] * fe_0 - 3.0 * ta2_xy_xx_xxyz_1[i] * fe_0 + 2.0 * ta2_xy_xxx_xyz_0[i] * fe_0 -
                                2.0 * ta2_xy_xxx_xyz_1[i] * fe_0 + ta1_y_xxx_xxyz_1[i] + ta2_xy_xxx_xxyz_0[i] * pa_x[i] -
                                ta2_xy_xxx_xxyz_1[i] * pc_x[i];

        ta2_xy_xxxx_xxzz_0[i] = 3.0 * ta2_xy_xx_xxzz_0[i] * fe_0 - 3.0 * ta2_xy_xx_xxzz_1[i] * fe_0 + 2.0 * ta2_xy_xxx_xzz_0[i] * fe_0 -
                                2.0 * ta2_xy_xxx_xzz_1[i] * fe_0 + ta1_y_xxx_xxzz_1[i] + ta2_xy_xxx_xxzz_0[i] * pa_x[i] -
                                ta2_xy_xxx_xxzz_1[i] * pc_x[i];

        ta2_xy_xxxx_xyyy_0[i] = 3.0 * ta2_xy_xx_xyyy_0[i] * fe_0 - 3.0 * ta2_xy_xx_xyyy_1[i] * fe_0 + ta2_xy_xxx_yyy_0[i] * fe_0 -
                                ta2_xy_xxx_yyy_1[i] * fe_0 + ta1_y_xxx_xyyy_1[i] + ta2_xy_xxx_xyyy_0[i] * pa_x[i] - ta2_xy_xxx_xyyy_1[i] * pc_x[i];

        ta2_xy_xxxx_xyyz_0[i] = 3.0 * ta2_xy_xx_xyyz_0[i] * fe_0 - 3.0 * ta2_xy_xx_xyyz_1[i] * fe_0 + ta2_xy_xxx_yyz_0[i] * fe_0 -
                                ta2_xy_xxx_yyz_1[i] * fe_0 + ta1_y_xxx_xyyz_1[i] + ta2_xy_xxx_xyyz_0[i] * pa_x[i] - ta2_xy_xxx_xyyz_1[i] * pc_x[i];

        ta2_xy_xxxx_xyzz_0[i] = 3.0 * ta2_xy_xx_xyzz_0[i] * fe_0 - 3.0 * ta2_xy_xx_xyzz_1[i] * fe_0 + ta2_xy_xxx_yzz_0[i] * fe_0 -
                                ta2_xy_xxx_yzz_1[i] * fe_0 + ta1_y_xxx_xyzz_1[i] + ta2_xy_xxx_xyzz_0[i] * pa_x[i] - ta2_xy_xxx_xyzz_1[i] * pc_x[i];

        ta2_xy_xxxx_xzzz_0[i] = 3.0 * ta2_xy_xx_xzzz_0[i] * fe_0 - 3.0 * ta2_xy_xx_xzzz_1[i] * fe_0 + ta2_xy_xxx_zzz_0[i] * fe_0 -
                                ta2_xy_xxx_zzz_1[i] * fe_0 + ta1_y_xxx_xzzz_1[i] + ta2_xy_xxx_xzzz_0[i] * pa_x[i] - ta2_xy_xxx_xzzz_1[i] * pc_x[i];

        ta2_xy_xxxx_yyyy_0[i] = 3.0 * ta2_xy_xx_yyyy_0[i] * fe_0 - 3.0 * ta2_xy_xx_yyyy_1[i] * fe_0 + ta1_y_xxx_yyyy_1[i] +
                                ta2_xy_xxx_yyyy_0[i] * pa_x[i] - ta2_xy_xxx_yyyy_1[i] * pc_x[i];

        ta2_xy_xxxx_yyyz_0[i] = 3.0 * ta2_xy_xx_yyyz_0[i] * fe_0 - 3.0 * ta2_xy_xx_yyyz_1[i] * fe_0 + ta1_y_xxx_yyyz_1[i] +
                                ta2_xy_xxx_yyyz_0[i] * pa_x[i] - ta2_xy_xxx_yyyz_1[i] * pc_x[i];

        ta2_xy_xxxx_yyzz_0[i] = 3.0 * ta2_xy_xx_yyzz_0[i] * fe_0 - 3.0 * ta2_xy_xx_yyzz_1[i] * fe_0 + ta1_y_xxx_yyzz_1[i] +
                                ta2_xy_xxx_yyzz_0[i] * pa_x[i] - ta2_xy_xxx_yyzz_1[i] * pc_x[i];

        ta2_xy_xxxx_yzzz_0[i] = 3.0 * ta2_xy_xx_yzzz_0[i] * fe_0 - 3.0 * ta2_xy_xx_yzzz_1[i] * fe_0 + ta1_y_xxx_yzzz_1[i] +
                                ta2_xy_xxx_yzzz_0[i] * pa_x[i] - ta2_xy_xxx_yzzz_1[i] * pc_x[i];

        ta2_xy_xxxx_zzzz_0[i] = 3.0 * ta2_xy_xx_zzzz_0[i] * fe_0 - 3.0 * ta2_xy_xx_zzzz_1[i] * fe_0 + ta1_y_xxx_zzzz_1[i] +
                                ta2_xy_xxx_zzzz_0[i] * pa_x[i] - ta2_xy_xxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 240-255 components of targeted buffer : GG

    auto ta2_xy_xxxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 240);

    auto ta2_xy_xxxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 241);

    auto ta2_xy_xxxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 242);

    auto ta2_xy_xxxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 243);

    auto ta2_xy_xxxy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 244);

    auto ta2_xy_xxxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 245);

    auto ta2_xy_xxxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 246);

    auto ta2_xy_xxxy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 247);

    auto ta2_xy_xxxy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 248);

    auto ta2_xy_xxxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 249);

    auto ta2_xy_xxxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 250);

    auto ta2_xy_xxxy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 251);

    auto ta2_xy_xxxy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 252);

    auto ta2_xy_xxxy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 253);

    auto ta2_xy_xxxy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 254);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_xxx_xxxx_1,   \
                             ta1_x_xxx_xxxy_1,   \
                             ta1_x_xxx_xxxz_1,   \
                             ta1_x_xxx_xxyy_1,   \
                             ta1_x_xxx_xxyz_1,   \
                             ta1_x_xxx_xxzz_1,   \
                             ta1_x_xxx_xyyy_1,   \
                             ta1_x_xxx_xyyz_1,   \
                             ta1_x_xxx_xyzz_1,   \
                             ta1_x_xxx_xzzz_1,   \
                             ta1_x_xxx_zzzz_1,   \
                             ta1_y_xxy_yyyy_1,   \
                             ta1_y_xxy_yyyz_1,   \
                             ta1_y_xxy_yyzz_1,   \
                             ta1_y_xxy_yzzz_1,   \
                             ta2_xy_xxx_xxx_0,   \
                             ta2_xy_xxx_xxx_1,   \
                             ta2_xy_xxx_xxxx_0,  \
                             ta2_xy_xxx_xxxx_1,  \
                             ta2_xy_xxx_xxxy_0,  \
                             ta2_xy_xxx_xxxy_1,  \
                             ta2_xy_xxx_xxxz_0,  \
                             ta2_xy_xxx_xxxz_1,  \
                             ta2_xy_xxx_xxy_0,   \
                             ta2_xy_xxx_xxy_1,   \
                             ta2_xy_xxx_xxyy_0,  \
                             ta2_xy_xxx_xxyy_1,  \
                             ta2_xy_xxx_xxyz_0,  \
                             ta2_xy_xxx_xxyz_1,  \
                             ta2_xy_xxx_xxz_0,   \
                             ta2_xy_xxx_xxz_1,   \
                             ta2_xy_xxx_xxzz_0,  \
                             ta2_xy_xxx_xxzz_1,  \
                             ta2_xy_xxx_xyy_0,   \
                             ta2_xy_xxx_xyy_1,   \
                             ta2_xy_xxx_xyyy_0,  \
                             ta2_xy_xxx_xyyy_1,  \
                             ta2_xy_xxx_xyyz_0,  \
                             ta2_xy_xxx_xyyz_1,  \
                             ta2_xy_xxx_xyz_0,   \
                             ta2_xy_xxx_xyz_1,   \
                             ta2_xy_xxx_xyzz_0,  \
                             ta2_xy_xxx_xyzz_1,  \
                             ta2_xy_xxx_xzz_0,   \
                             ta2_xy_xxx_xzz_1,   \
                             ta2_xy_xxx_xzzz_0,  \
                             ta2_xy_xxx_xzzz_1,  \
                             ta2_xy_xxx_zzzz_0,  \
                             ta2_xy_xxx_zzzz_1,  \
                             ta2_xy_xxxy_xxxx_0, \
                             ta2_xy_xxxy_xxxy_0, \
                             ta2_xy_xxxy_xxxz_0, \
                             ta2_xy_xxxy_xxyy_0, \
                             ta2_xy_xxxy_xxyz_0, \
                             ta2_xy_xxxy_xxzz_0, \
                             ta2_xy_xxxy_xyyy_0, \
                             ta2_xy_xxxy_xyyz_0, \
                             ta2_xy_xxxy_xyzz_0, \
                             ta2_xy_xxxy_xzzz_0, \
                             ta2_xy_xxxy_yyyy_0, \
                             ta2_xy_xxxy_yyyz_0, \
                             ta2_xy_xxxy_yyzz_0, \
                             ta2_xy_xxxy_yzzz_0, \
                             ta2_xy_xxxy_zzzz_0, \
                             ta2_xy_xxy_yyyy_0,  \
                             ta2_xy_xxy_yyyy_1,  \
                             ta2_xy_xxy_yyyz_0,  \
                             ta2_xy_xxy_yyyz_1,  \
                             ta2_xy_xxy_yyzz_0,  \
                             ta2_xy_xxy_yyzz_1,  \
                             ta2_xy_xxy_yzzz_0,  \
                             ta2_xy_xxy_yzzz_1,  \
                             ta2_xy_xy_yyyy_0,   \
                             ta2_xy_xy_yyyy_1,   \
                             ta2_xy_xy_yyyz_0,   \
                             ta2_xy_xy_yyyz_1,   \
                             ta2_xy_xy_yyzz_0,   \
                             ta2_xy_xy_yyzz_1,   \
                             ta2_xy_xy_yzzz_0,   \
                             ta2_xy_xy_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxxy_xxxx_0[i] = ta1_x_xxx_xxxx_1[i] + ta2_xy_xxx_xxxx_0[i] * pa_y[i] - ta2_xy_xxx_xxxx_1[i] * pc_y[i];

        ta2_xy_xxxy_xxxy_0[i] = ta2_xy_xxx_xxx_0[i] * fe_0 - ta2_xy_xxx_xxx_1[i] * fe_0 + ta1_x_xxx_xxxy_1[i] + ta2_xy_xxx_xxxy_0[i] * pa_y[i] -
                                ta2_xy_xxx_xxxy_1[i] * pc_y[i];

        ta2_xy_xxxy_xxxz_0[i] = ta1_x_xxx_xxxz_1[i] + ta2_xy_xxx_xxxz_0[i] * pa_y[i] - ta2_xy_xxx_xxxz_1[i] * pc_y[i];

        ta2_xy_xxxy_xxyy_0[i] = 2.0 * ta2_xy_xxx_xxy_0[i] * fe_0 - 2.0 * ta2_xy_xxx_xxy_1[i] * fe_0 + ta1_x_xxx_xxyy_1[i] +
                                ta2_xy_xxx_xxyy_0[i] * pa_y[i] - ta2_xy_xxx_xxyy_1[i] * pc_y[i];

        ta2_xy_xxxy_xxyz_0[i] = ta2_xy_xxx_xxz_0[i] * fe_0 - ta2_xy_xxx_xxz_1[i] * fe_0 + ta1_x_xxx_xxyz_1[i] + ta2_xy_xxx_xxyz_0[i] * pa_y[i] -
                                ta2_xy_xxx_xxyz_1[i] * pc_y[i];

        ta2_xy_xxxy_xxzz_0[i] = ta1_x_xxx_xxzz_1[i] + ta2_xy_xxx_xxzz_0[i] * pa_y[i] - ta2_xy_xxx_xxzz_1[i] * pc_y[i];

        ta2_xy_xxxy_xyyy_0[i] = 3.0 * ta2_xy_xxx_xyy_0[i] * fe_0 - 3.0 * ta2_xy_xxx_xyy_1[i] * fe_0 + ta1_x_xxx_xyyy_1[i] +
                                ta2_xy_xxx_xyyy_0[i] * pa_y[i] - ta2_xy_xxx_xyyy_1[i] * pc_y[i];

        ta2_xy_xxxy_xyyz_0[i] = 2.0 * ta2_xy_xxx_xyz_0[i] * fe_0 - 2.0 * ta2_xy_xxx_xyz_1[i] * fe_0 + ta1_x_xxx_xyyz_1[i] +
                                ta2_xy_xxx_xyyz_0[i] * pa_y[i] - ta2_xy_xxx_xyyz_1[i] * pc_y[i];

        ta2_xy_xxxy_xyzz_0[i] = ta2_xy_xxx_xzz_0[i] * fe_0 - ta2_xy_xxx_xzz_1[i] * fe_0 + ta1_x_xxx_xyzz_1[i] + ta2_xy_xxx_xyzz_0[i] * pa_y[i] -
                                ta2_xy_xxx_xyzz_1[i] * pc_y[i];

        ta2_xy_xxxy_xzzz_0[i] = ta1_x_xxx_xzzz_1[i] + ta2_xy_xxx_xzzz_0[i] * pa_y[i] - ta2_xy_xxx_xzzz_1[i] * pc_y[i];

        ta2_xy_xxxy_yyyy_0[i] = 2.0 * ta2_xy_xy_yyyy_0[i] * fe_0 - 2.0 * ta2_xy_xy_yyyy_1[i] * fe_0 + ta1_y_xxy_yyyy_1[i] +
                                ta2_xy_xxy_yyyy_0[i] * pa_x[i] - ta2_xy_xxy_yyyy_1[i] * pc_x[i];

        ta2_xy_xxxy_yyyz_0[i] = 2.0 * ta2_xy_xy_yyyz_0[i] * fe_0 - 2.0 * ta2_xy_xy_yyyz_1[i] * fe_0 + ta1_y_xxy_yyyz_1[i] +
                                ta2_xy_xxy_yyyz_0[i] * pa_x[i] - ta2_xy_xxy_yyyz_1[i] * pc_x[i];

        ta2_xy_xxxy_yyzz_0[i] = 2.0 * ta2_xy_xy_yyzz_0[i] * fe_0 - 2.0 * ta2_xy_xy_yyzz_1[i] * fe_0 + ta1_y_xxy_yyzz_1[i] +
                                ta2_xy_xxy_yyzz_0[i] * pa_x[i] - ta2_xy_xxy_yyzz_1[i] * pc_x[i];

        ta2_xy_xxxy_yzzz_0[i] = 2.0 * ta2_xy_xy_yzzz_0[i] * fe_0 - 2.0 * ta2_xy_xy_yzzz_1[i] * fe_0 + ta1_y_xxy_yzzz_1[i] +
                                ta2_xy_xxy_yzzz_0[i] * pa_x[i] - ta2_xy_xxy_yzzz_1[i] * pc_x[i];

        ta2_xy_xxxy_zzzz_0[i] = ta1_x_xxx_zzzz_1[i] + ta2_xy_xxx_zzzz_0[i] * pa_y[i] - ta2_xy_xxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 255-270 components of targeted buffer : GG

    auto ta2_xy_xxxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 255);

    auto ta2_xy_xxxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 256);

    auto ta2_xy_xxxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 257);

    auto ta2_xy_xxxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 258);

    auto ta2_xy_xxxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 259);

    auto ta2_xy_xxxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 260);

    auto ta2_xy_xxxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 261);

    auto ta2_xy_xxxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 262);

    auto ta2_xy_xxxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 263);

    auto ta2_xy_xxxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 264);

    auto ta2_xy_xxxz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 265);

    auto ta2_xy_xxxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 266);

    auto ta2_xy_xxxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 267);

    auto ta2_xy_xxxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 268);

    auto ta2_xy_xxxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 269);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta2_xy_xxx_xxx_0,   \
                             ta2_xy_xxx_xxx_1,   \
                             ta2_xy_xxx_xxxx_0,  \
                             ta2_xy_xxx_xxxx_1,  \
                             ta2_xy_xxx_xxxy_0,  \
                             ta2_xy_xxx_xxxy_1,  \
                             ta2_xy_xxx_xxxz_0,  \
                             ta2_xy_xxx_xxxz_1,  \
                             ta2_xy_xxx_xxy_0,   \
                             ta2_xy_xxx_xxy_1,   \
                             ta2_xy_xxx_xxyy_0,  \
                             ta2_xy_xxx_xxyy_1,  \
                             ta2_xy_xxx_xxyz_0,  \
                             ta2_xy_xxx_xxyz_1,  \
                             ta2_xy_xxx_xxz_0,   \
                             ta2_xy_xxx_xxz_1,   \
                             ta2_xy_xxx_xxzz_0,  \
                             ta2_xy_xxx_xxzz_1,  \
                             ta2_xy_xxx_xyy_0,   \
                             ta2_xy_xxx_xyy_1,   \
                             ta2_xy_xxx_xyyy_0,  \
                             ta2_xy_xxx_xyyy_1,  \
                             ta2_xy_xxx_xyyz_0,  \
                             ta2_xy_xxx_xyyz_1,  \
                             ta2_xy_xxx_xyz_0,   \
                             ta2_xy_xxx_xyz_1,   \
                             ta2_xy_xxx_xyzz_0,  \
                             ta2_xy_xxx_xyzz_1,  \
                             ta2_xy_xxx_xzz_0,   \
                             ta2_xy_xxx_xzz_1,   \
                             ta2_xy_xxx_xzzz_0,  \
                             ta2_xy_xxx_xzzz_1,  \
                             ta2_xy_xxx_yyy_0,   \
                             ta2_xy_xxx_yyy_1,   \
                             ta2_xy_xxx_yyyy_0,  \
                             ta2_xy_xxx_yyyy_1,  \
                             ta2_xy_xxx_yyyz_0,  \
                             ta2_xy_xxx_yyyz_1,  \
                             ta2_xy_xxx_yyz_0,   \
                             ta2_xy_xxx_yyz_1,   \
                             ta2_xy_xxx_yyzz_0,  \
                             ta2_xy_xxx_yyzz_1,  \
                             ta2_xy_xxx_yzz_0,   \
                             ta2_xy_xxx_yzz_1,   \
                             ta2_xy_xxx_yzzz_0,  \
                             ta2_xy_xxx_yzzz_1,  \
                             ta2_xy_xxx_zzz_0,   \
                             ta2_xy_xxx_zzz_1,   \
                             ta2_xy_xxx_zzzz_0,  \
                             ta2_xy_xxx_zzzz_1,  \
                             ta2_xy_xxxz_xxxx_0, \
                             ta2_xy_xxxz_xxxy_0, \
                             ta2_xy_xxxz_xxxz_0, \
                             ta2_xy_xxxz_xxyy_0, \
                             ta2_xy_xxxz_xxyz_0, \
                             ta2_xy_xxxz_xxzz_0, \
                             ta2_xy_xxxz_xyyy_0, \
                             ta2_xy_xxxz_xyyz_0, \
                             ta2_xy_xxxz_xyzz_0, \
                             ta2_xy_xxxz_xzzz_0, \
                             ta2_xy_xxxz_yyyy_0, \
                             ta2_xy_xxxz_yyyz_0, \
                             ta2_xy_xxxz_yyzz_0, \
                             ta2_xy_xxxz_yzzz_0, \
                             ta2_xy_xxxz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxxz_xxxx_0[i] = ta2_xy_xxx_xxxx_0[i] * pa_z[i] - ta2_xy_xxx_xxxx_1[i] * pc_z[i];

        ta2_xy_xxxz_xxxy_0[i] = ta2_xy_xxx_xxxy_0[i] * pa_z[i] - ta2_xy_xxx_xxxy_1[i] * pc_z[i];

        ta2_xy_xxxz_xxxz_0[i] =
            ta2_xy_xxx_xxx_0[i] * fe_0 - ta2_xy_xxx_xxx_1[i] * fe_0 + ta2_xy_xxx_xxxz_0[i] * pa_z[i] - ta2_xy_xxx_xxxz_1[i] * pc_z[i];

        ta2_xy_xxxz_xxyy_0[i] = ta2_xy_xxx_xxyy_0[i] * pa_z[i] - ta2_xy_xxx_xxyy_1[i] * pc_z[i];

        ta2_xy_xxxz_xxyz_0[i] =
            ta2_xy_xxx_xxy_0[i] * fe_0 - ta2_xy_xxx_xxy_1[i] * fe_0 + ta2_xy_xxx_xxyz_0[i] * pa_z[i] - ta2_xy_xxx_xxyz_1[i] * pc_z[i];

        ta2_xy_xxxz_xxzz_0[i] =
            2.0 * ta2_xy_xxx_xxz_0[i] * fe_0 - 2.0 * ta2_xy_xxx_xxz_1[i] * fe_0 + ta2_xy_xxx_xxzz_0[i] * pa_z[i] - ta2_xy_xxx_xxzz_1[i] * pc_z[i];

        ta2_xy_xxxz_xyyy_0[i] = ta2_xy_xxx_xyyy_0[i] * pa_z[i] - ta2_xy_xxx_xyyy_1[i] * pc_z[i];

        ta2_xy_xxxz_xyyz_0[i] =
            ta2_xy_xxx_xyy_0[i] * fe_0 - ta2_xy_xxx_xyy_1[i] * fe_0 + ta2_xy_xxx_xyyz_0[i] * pa_z[i] - ta2_xy_xxx_xyyz_1[i] * pc_z[i];

        ta2_xy_xxxz_xyzz_0[i] =
            2.0 * ta2_xy_xxx_xyz_0[i] * fe_0 - 2.0 * ta2_xy_xxx_xyz_1[i] * fe_0 + ta2_xy_xxx_xyzz_0[i] * pa_z[i] - ta2_xy_xxx_xyzz_1[i] * pc_z[i];

        ta2_xy_xxxz_xzzz_0[i] =
            3.0 * ta2_xy_xxx_xzz_0[i] * fe_0 - 3.0 * ta2_xy_xxx_xzz_1[i] * fe_0 + ta2_xy_xxx_xzzz_0[i] * pa_z[i] - ta2_xy_xxx_xzzz_1[i] * pc_z[i];

        ta2_xy_xxxz_yyyy_0[i] = ta2_xy_xxx_yyyy_0[i] * pa_z[i] - ta2_xy_xxx_yyyy_1[i] * pc_z[i];

        ta2_xy_xxxz_yyyz_0[i] =
            ta2_xy_xxx_yyy_0[i] * fe_0 - ta2_xy_xxx_yyy_1[i] * fe_0 + ta2_xy_xxx_yyyz_0[i] * pa_z[i] - ta2_xy_xxx_yyyz_1[i] * pc_z[i];

        ta2_xy_xxxz_yyzz_0[i] =
            2.0 * ta2_xy_xxx_yyz_0[i] * fe_0 - 2.0 * ta2_xy_xxx_yyz_1[i] * fe_0 + ta2_xy_xxx_yyzz_0[i] * pa_z[i] - ta2_xy_xxx_yyzz_1[i] * pc_z[i];

        ta2_xy_xxxz_yzzz_0[i] =
            3.0 * ta2_xy_xxx_yzz_0[i] * fe_0 - 3.0 * ta2_xy_xxx_yzz_1[i] * fe_0 + ta2_xy_xxx_yzzz_0[i] * pa_z[i] - ta2_xy_xxx_yzzz_1[i] * pc_z[i];

        ta2_xy_xxxz_zzzz_0[i] =
            4.0 * ta2_xy_xxx_zzz_0[i] * fe_0 - 4.0 * ta2_xy_xxx_zzz_1[i] * fe_0 + ta2_xy_xxx_zzzz_0[i] * pa_z[i] - ta2_xy_xxx_zzzz_1[i] * pc_z[i];
    }

    // Set up 270-285 components of targeted buffer : GG

    auto ta2_xy_xxyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 270);

    auto ta2_xy_xxyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 271);

    auto ta2_xy_xxyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 272);

    auto ta2_xy_xxyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 273);

    auto ta2_xy_xxyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 274);

    auto ta2_xy_xxyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 275);

    auto ta2_xy_xxyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 276);

    auto ta2_xy_xxyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 277);

    auto ta2_xy_xxyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 278);

    auto ta2_xy_xxyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 279);

    auto ta2_xy_xxyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 280);

    auto ta2_xy_xxyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 281);

    auto ta2_xy_xxyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 282);

    auto ta2_xy_xxyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 283);

    auto ta2_xy_xxyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 284);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_xxy_xxxx_1,   \
                             ta1_x_xxy_xxxz_1,   \
                             ta1_x_xxy_xxzz_1,   \
                             ta1_x_xxy_xzzz_1,   \
                             ta1_y_xyy_xxxy_1,   \
                             ta1_y_xyy_xxyy_1,   \
                             ta1_y_xyy_xxyz_1,   \
                             ta1_y_xyy_xyyy_1,   \
                             ta1_y_xyy_xyyz_1,   \
                             ta1_y_xyy_xyzz_1,   \
                             ta1_y_xyy_yyyy_1,   \
                             ta1_y_xyy_yyyz_1,   \
                             ta1_y_xyy_yyzz_1,   \
                             ta1_y_xyy_yzzz_1,   \
                             ta1_y_xyy_zzzz_1,   \
                             ta2_xy_xx_xxxx_0,   \
                             ta2_xy_xx_xxxx_1,   \
                             ta2_xy_xx_xxxz_0,   \
                             ta2_xy_xx_xxxz_1,   \
                             ta2_xy_xx_xxzz_0,   \
                             ta2_xy_xx_xxzz_1,   \
                             ta2_xy_xx_xzzz_0,   \
                             ta2_xy_xx_xzzz_1,   \
                             ta2_xy_xxy_xxxx_0,  \
                             ta2_xy_xxy_xxxx_1,  \
                             ta2_xy_xxy_xxxz_0,  \
                             ta2_xy_xxy_xxxz_1,  \
                             ta2_xy_xxy_xxzz_0,  \
                             ta2_xy_xxy_xxzz_1,  \
                             ta2_xy_xxy_xzzz_0,  \
                             ta2_xy_xxy_xzzz_1,  \
                             ta2_xy_xxyy_xxxx_0, \
                             ta2_xy_xxyy_xxxy_0, \
                             ta2_xy_xxyy_xxxz_0, \
                             ta2_xy_xxyy_xxyy_0, \
                             ta2_xy_xxyy_xxyz_0, \
                             ta2_xy_xxyy_xxzz_0, \
                             ta2_xy_xxyy_xyyy_0, \
                             ta2_xy_xxyy_xyyz_0, \
                             ta2_xy_xxyy_xyzz_0, \
                             ta2_xy_xxyy_xzzz_0, \
                             ta2_xy_xxyy_yyyy_0, \
                             ta2_xy_xxyy_yyyz_0, \
                             ta2_xy_xxyy_yyzz_0, \
                             ta2_xy_xxyy_yzzz_0, \
                             ta2_xy_xxyy_zzzz_0, \
                             ta2_xy_xyy_xxxy_0,  \
                             ta2_xy_xyy_xxxy_1,  \
                             ta2_xy_xyy_xxy_0,   \
                             ta2_xy_xyy_xxy_1,   \
                             ta2_xy_xyy_xxyy_0,  \
                             ta2_xy_xyy_xxyy_1,  \
                             ta2_xy_xyy_xxyz_0,  \
                             ta2_xy_xyy_xxyz_1,  \
                             ta2_xy_xyy_xyy_0,   \
                             ta2_xy_xyy_xyy_1,   \
                             ta2_xy_xyy_xyyy_0,  \
                             ta2_xy_xyy_xyyy_1,  \
                             ta2_xy_xyy_xyyz_0,  \
                             ta2_xy_xyy_xyyz_1,  \
                             ta2_xy_xyy_xyz_0,   \
                             ta2_xy_xyy_xyz_1,   \
                             ta2_xy_xyy_xyzz_0,  \
                             ta2_xy_xyy_xyzz_1,  \
                             ta2_xy_xyy_yyy_0,   \
                             ta2_xy_xyy_yyy_1,   \
                             ta2_xy_xyy_yyyy_0,  \
                             ta2_xy_xyy_yyyy_1,  \
                             ta2_xy_xyy_yyyz_0,  \
                             ta2_xy_xyy_yyyz_1,  \
                             ta2_xy_xyy_yyz_0,   \
                             ta2_xy_xyy_yyz_1,   \
                             ta2_xy_xyy_yyzz_0,  \
                             ta2_xy_xyy_yyzz_1,  \
                             ta2_xy_xyy_yzz_0,   \
                             ta2_xy_xyy_yzz_1,   \
                             ta2_xy_xyy_yzzz_0,  \
                             ta2_xy_xyy_yzzz_1,  \
                             ta2_xy_xyy_zzzz_0,  \
                             ta2_xy_xyy_zzzz_1,  \
                             ta2_xy_yy_xxxy_0,   \
                             ta2_xy_yy_xxxy_1,   \
                             ta2_xy_yy_xxyy_0,   \
                             ta2_xy_yy_xxyy_1,   \
                             ta2_xy_yy_xxyz_0,   \
                             ta2_xy_yy_xxyz_1,   \
                             ta2_xy_yy_xyyy_0,   \
                             ta2_xy_yy_xyyy_1,   \
                             ta2_xy_yy_xyyz_0,   \
                             ta2_xy_yy_xyyz_1,   \
                             ta2_xy_yy_xyzz_0,   \
                             ta2_xy_yy_xyzz_1,   \
                             ta2_xy_yy_yyyy_0,   \
                             ta2_xy_yy_yyyy_1,   \
                             ta2_xy_yy_yyyz_0,   \
                             ta2_xy_yy_yyyz_1,   \
                             ta2_xy_yy_yyzz_0,   \
                             ta2_xy_yy_yyzz_1,   \
                             ta2_xy_yy_yzzz_0,   \
                             ta2_xy_yy_yzzz_1,   \
                             ta2_xy_yy_zzzz_0,   \
                             ta2_xy_yy_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxyy_xxxx_0[i] = ta2_xy_xx_xxxx_0[i] * fe_0 - ta2_xy_xx_xxxx_1[i] * fe_0 + ta1_x_xxy_xxxx_1[i] + ta2_xy_xxy_xxxx_0[i] * pa_y[i] -
                                ta2_xy_xxy_xxxx_1[i] * pc_y[i];

        ta2_xy_xxyy_xxxy_0[i] = ta2_xy_yy_xxxy_0[i] * fe_0 - ta2_xy_yy_xxxy_1[i] * fe_0 + 3.0 * ta2_xy_xyy_xxy_0[i] * fe_0 -
                                3.0 * ta2_xy_xyy_xxy_1[i] * fe_0 + ta1_y_xyy_xxxy_1[i] + ta2_xy_xyy_xxxy_0[i] * pa_x[i] -
                                ta2_xy_xyy_xxxy_1[i] * pc_x[i];

        ta2_xy_xxyy_xxxz_0[i] = ta2_xy_xx_xxxz_0[i] * fe_0 - ta2_xy_xx_xxxz_1[i] * fe_0 + ta1_x_xxy_xxxz_1[i] + ta2_xy_xxy_xxxz_0[i] * pa_y[i] -
                                ta2_xy_xxy_xxxz_1[i] * pc_y[i];

        ta2_xy_xxyy_xxyy_0[i] = ta2_xy_yy_xxyy_0[i] * fe_0 - ta2_xy_yy_xxyy_1[i] * fe_0 + 2.0 * ta2_xy_xyy_xyy_0[i] * fe_0 -
                                2.0 * ta2_xy_xyy_xyy_1[i] * fe_0 + ta1_y_xyy_xxyy_1[i] + ta2_xy_xyy_xxyy_0[i] * pa_x[i] -
                                ta2_xy_xyy_xxyy_1[i] * pc_x[i];

        ta2_xy_xxyy_xxyz_0[i] = ta2_xy_yy_xxyz_0[i] * fe_0 - ta2_xy_yy_xxyz_1[i] * fe_0 + 2.0 * ta2_xy_xyy_xyz_0[i] * fe_0 -
                                2.0 * ta2_xy_xyy_xyz_1[i] * fe_0 + ta1_y_xyy_xxyz_1[i] + ta2_xy_xyy_xxyz_0[i] * pa_x[i] -
                                ta2_xy_xyy_xxyz_1[i] * pc_x[i];

        ta2_xy_xxyy_xxzz_0[i] = ta2_xy_xx_xxzz_0[i] * fe_0 - ta2_xy_xx_xxzz_1[i] * fe_0 + ta1_x_xxy_xxzz_1[i] + ta2_xy_xxy_xxzz_0[i] * pa_y[i] -
                                ta2_xy_xxy_xxzz_1[i] * pc_y[i];

        ta2_xy_xxyy_xyyy_0[i] = ta2_xy_yy_xyyy_0[i] * fe_0 - ta2_xy_yy_xyyy_1[i] * fe_0 + ta2_xy_xyy_yyy_0[i] * fe_0 - ta2_xy_xyy_yyy_1[i] * fe_0 +
                                ta1_y_xyy_xyyy_1[i] + ta2_xy_xyy_xyyy_0[i] * pa_x[i] - ta2_xy_xyy_xyyy_1[i] * pc_x[i];

        ta2_xy_xxyy_xyyz_0[i] = ta2_xy_yy_xyyz_0[i] * fe_0 - ta2_xy_yy_xyyz_1[i] * fe_0 + ta2_xy_xyy_yyz_0[i] * fe_0 - ta2_xy_xyy_yyz_1[i] * fe_0 +
                                ta1_y_xyy_xyyz_1[i] + ta2_xy_xyy_xyyz_0[i] * pa_x[i] - ta2_xy_xyy_xyyz_1[i] * pc_x[i];

        ta2_xy_xxyy_xyzz_0[i] = ta2_xy_yy_xyzz_0[i] * fe_0 - ta2_xy_yy_xyzz_1[i] * fe_0 + ta2_xy_xyy_yzz_0[i] * fe_0 - ta2_xy_xyy_yzz_1[i] * fe_0 +
                                ta1_y_xyy_xyzz_1[i] + ta2_xy_xyy_xyzz_0[i] * pa_x[i] - ta2_xy_xyy_xyzz_1[i] * pc_x[i];

        ta2_xy_xxyy_xzzz_0[i] = ta2_xy_xx_xzzz_0[i] * fe_0 - ta2_xy_xx_xzzz_1[i] * fe_0 + ta1_x_xxy_xzzz_1[i] + ta2_xy_xxy_xzzz_0[i] * pa_y[i] -
                                ta2_xy_xxy_xzzz_1[i] * pc_y[i];

        ta2_xy_xxyy_yyyy_0[i] = ta2_xy_yy_yyyy_0[i] * fe_0 - ta2_xy_yy_yyyy_1[i] * fe_0 + ta1_y_xyy_yyyy_1[i] + ta2_xy_xyy_yyyy_0[i] * pa_x[i] -
                                ta2_xy_xyy_yyyy_1[i] * pc_x[i];

        ta2_xy_xxyy_yyyz_0[i] = ta2_xy_yy_yyyz_0[i] * fe_0 - ta2_xy_yy_yyyz_1[i] * fe_0 + ta1_y_xyy_yyyz_1[i] + ta2_xy_xyy_yyyz_0[i] * pa_x[i] -
                                ta2_xy_xyy_yyyz_1[i] * pc_x[i];

        ta2_xy_xxyy_yyzz_0[i] = ta2_xy_yy_yyzz_0[i] * fe_0 - ta2_xy_yy_yyzz_1[i] * fe_0 + ta1_y_xyy_yyzz_1[i] + ta2_xy_xyy_yyzz_0[i] * pa_x[i] -
                                ta2_xy_xyy_yyzz_1[i] * pc_x[i];

        ta2_xy_xxyy_yzzz_0[i] = ta2_xy_yy_yzzz_0[i] * fe_0 - ta2_xy_yy_yzzz_1[i] * fe_0 + ta1_y_xyy_yzzz_1[i] + ta2_xy_xyy_yzzz_0[i] * pa_x[i] -
                                ta2_xy_xyy_yzzz_1[i] * pc_x[i];

        ta2_xy_xxyy_zzzz_0[i] = ta2_xy_yy_zzzz_0[i] * fe_0 - ta2_xy_yy_zzzz_1[i] * fe_0 + ta1_y_xyy_zzzz_1[i] + ta2_xy_xyy_zzzz_0[i] * pa_x[i] -
                                ta2_xy_xyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 285-300 components of targeted buffer : GG

    auto ta2_xy_xxyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 285);

    auto ta2_xy_xxyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 286);

    auto ta2_xy_xxyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 287);

    auto ta2_xy_xxyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 288);

    auto ta2_xy_xxyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 289);

    auto ta2_xy_xxyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 290);

    auto ta2_xy_xxyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 291);

    auto ta2_xy_xxyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 292);

    auto ta2_xy_xxyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 293);

    auto ta2_xy_xxyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 294);

    auto ta2_xy_xxyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 295);

    auto ta2_xy_xxyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 296);

    auto ta2_xy_xxyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 297);

    auto ta2_xy_xxyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 298);

    auto ta2_xy_xxyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 299);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xxz_xxxz_1,   \
                             ta1_x_xxz_xxzz_1,   \
                             ta1_x_xxz_xzzz_1,   \
                             ta1_x_xxz_zzzz_1,   \
                             ta2_xy_xxy_xxxx_0,  \
                             ta2_xy_xxy_xxxx_1,  \
                             ta2_xy_xxy_xxxy_0,  \
                             ta2_xy_xxy_xxxy_1,  \
                             ta2_xy_xxy_xxy_0,   \
                             ta2_xy_xxy_xxy_1,   \
                             ta2_xy_xxy_xxyy_0,  \
                             ta2_xy_xxy_xxyy_1,  \
                             ta2_xy_xxy_xxyz_0,  \
                             ta2_xy_xxy_xxyz_1,  \
                             ta2_xy_xxy_xyy_0,   \
                             ta2_xy_xxy_xyy_1,   \
                             ta2_xy_xxy_xyyy_0,  \
                             ta2_xy_xxy_xyyy_1,  \
                             ta2_xy_xxy_xyyz_0,  \
                             ta2_xy_xxy_xyyz_1,  \
                             ta2_xy_xxy_xyz_0,   \
                             ta2_xy_xxy_xyz_1,   \
                             ta2_xy_xxy_xyzz_0,  \
                             ta2_xy_xxy_xyzz_1,  \
                             ta2_xy_xxy_yyy_0,   \
                             ta2_xy_xxy_yyy_1,   \
                             ta2_xy_xxy_yyyy_0,  \
                             ta2_xy_xxy_yyyy_1,  \
                             ta2_xy_xxy_yyyz_0,  \
                             ta2_xy_xxy_yyyz_1,  \
                             ta2_xy_xxy_yyz_0,   \
                             ta2_xy_xxy_yyz_1,   \
                             ta2_xy_xxy_yyzz_0,  \
                             ta2_xy_xxy_yyzz_1,  \
                             ta2_xy_xxy_yzz_0,   \
                             ta2_xy_xxy_yzz_1,   \
                             ta2_xy_xxy_yzzz_0,  \
                             ta2_xy_xxy_yzzz_1,  \
                             ta2_xy_xxyz_xxxx_0, \
                             ta2_xy_xxyz_xxxy_0, \
                             ta2_xy_xxyz_xxxz_0, \
                             ta2_xy_xxyz_xxyy_0, \
                             ta2_xy_xxyz_xxyz_0, \
                             ta2_xy_xxyz_xxzz_0, \
                             ta2_xy_xxyz_xyyy_0, \
                             ta2_xy_xxyz_xyyz_0, \
                             ta2_xy_xxyz_xyzz_0, \
                             ta2_xy_xxyz_xzzz_0, \
                             ta2_xy_xxyz_yyyy_0, \
                             ta2_xy_xxyz_yyyz_0, \
                             ta2_xy_xxyz_yyzz_0, \
                             ta2_xy_xxyz_yzzz_0, \
                             ta2_xy_xxyz_zzzz_0, \
                             ta2_xy_xxz_xxxz_0,  \
                             ta2_xy_xxz_xxxz_1,  \
                             ta2_xy_xxz_xxzz_0,  \
                             ta2_xy_xxz_xxzz_1,  \
                             ta2_xy_xxz_xzzz_0,  \
                             ta2_xy_xxz_xzzz_1,  \
                             ta2_xy_xxz_zzzz_0,  \
                             ta2_xy_xxz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxyz_xxxx_0[i] = ta2_xy_xxy_xxxx_0[i] * pa_z[i] - ta2_xy_xxy_xxxx_1[i] * pc_z[i];

        ta2_xy_xxyz_xxxy_0[i] = ta2_xy_xxy_xxxy_0[i] * pa_z[i] - ta2_xy_xxy_xxxy_1[i] * pc_z[i];

        ta2_xy_xxyz_xxxz_0[i] = ta1_x_xxz_xxxz_1[i] + ta2_xy_xxz_xxxz_0[i] * pa_y[i] - ta2_xy_xxz_xxxz_1[i] * pc_y[i];

        ta2_xy_xxyz_xxyy_0[i] = ta2_xy_xxy_xxyy_0[i] * pa_z[i] - ta2_xy_xxy_xxyy_1[i] * pc_z[i];

        ta2_xy_xxyz_xxyz_0[i] =
            ta2_xy_xxy_xxy_0[i] * fe_0 - ta2_xy_xxy_xxy_1[i] * fe_0 + ta2_xy_xxy_xxyz_0[i] * pa_z[i] - ta2_xy_xxy_xxyz_1[i] * pc_z[i];

        ta2_xy_xxyz_xxzz_0[i] = ta1_x_xxz_xxzz_1[i] + ta2_xy_xxz_xxzz_0[i] * pa_y[i] - ta2_xy_xxz_xxzz_1[i] * pc_y[i];

        ta2_xy_xxyz_xyyy_0[i] = ta2_xy_xxy_xyyy_0[i] * pa_z[i] - ta2_xy_xxy_xyyy_1[i] * pc_z[i];

        ta2_xy_xxyz_xyyz_0[i] =
            ta2_xy_xxy_xyy_0[i] * fe_0 - ta2_xy_xxy_xyy_1[i] * fe_0 + ta2_xy_xxy_xyyz_0[i] * pa_z[i] - ta2_xy_xxy_xyyz_1[i] * pc_z[i];

        ta2_xy_xxyz_xyzz_0[i] =
            2.0 * ta2_xy_xxy_xyz_0[i] * fe_0 - 2.0 * ta2_xy_xxy_xyz_1[i] * fe_0 + ta2_xy_xxy_xyzz_0[i] * pa_z[i] - ta2_xy_xxy_xyzz_1[i] * pc_z[i];

        ta2_xy_xxyz_xzzz_0[i] = ta1_x_xxz_xzzz_1[i] + ta2_xy_xxz_xzzz_0[i] * pa_y[i] - ta2_xy_xxz_xzzz_1[i] * pc_y[i];

        ta2_xy_xxyz_yyyy_0[i] = ta2_xy_xxy_yyyy_0[i] * pa_z[i] - ta2_xy_xxy_yyyy_1[i] * pc_z[i];

        ta2_xy_xxyz_yyyz_0[i] =
            ta2_xy_xxy_yyy_0[i] * fe_0 - ta2_xy_xxy_yyy_1[i] * fe_0 + ta2_xy_xxy_yyyz_0[i] * pa_z[i] - ta2_xy_xxy_yyyz_1[i] * pc_z[i];

        ta2_xy_xxyz_yyzz_0[i] =
            2.0 * ta2_xy_xxy_yyz_0[i] * fe_0 - 2.0 * ta2_xy_xxy_yyz_1[i] * fe_0 + ta2_xy_xxy_yyzz_0[i] * pa_z[i] - ta2_xy_xxy_yyzz_1[i] * pc_z[i];

        ta2_xy_xxyz_yzzz_0[i] =
            3.0 * ta2_xy_xxy_yzz_0[i] * fe_0 - 3.0 * ta2_xy_xxy_yzz_1[i] * fe_0 + ta2_xy_xxy_yzzz_0[i] * pa_z[i] - ta2_xy_xxy_yzzz_1[i] * pc_z[i];

        ta2_xy_xxyz_zzzz_0[i] = ta1_x_xxz_zzzz_1[i] + ta2_xy_xxz_zzzz_0[i] * pa_y[i] - ta2_xy_xxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 300-315 components of targeted buffer : GG

    auto ta2_xy_xxzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 300);

    auto ta2_xy_xxzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 301);

    auto ta2_xy_xxzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 302);

    auto ta2_xy_xxzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 303);

    auto ta2_xy_xxzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 304);

    auto ta2_xy_xxzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 305);

    auto ta2_xy_xxzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 306);

    auto ta2_xy_xxzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 307);

    auto ta2_xy_xxzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 308);

    auto ta2_xy_xxzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 309);

    auto ta2_xy_xxzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 310);

    auto ta2_xy_xxzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 311);

    auto ta2_xy_xxzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 312);

    auto ta2_xy_xxzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 313);

    auto ta2_xy_xxzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 314);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xzz_yyyz_1,   \
                             ta1_y_xzz_yyzz_1,   \
                             ta1_y_xzz_yzzz_1,   \
                             ta1_y_xzz_zzzz_1,   \
                             ta2_xy_xx_xxxx_0,   \
                             ta2_xy_xx_xxxx_1,   \
                             ta2_xy_xx_xxxy_0,   \
                             ta2_xy_xx_xxxy_1,   \
                             ta2_xy_xx_xxxz_0,   \
                             ta2_xy_xx_xxxz_1,   \
                             ta2_xy_xx_xxyy_0,   \
                             ta2_xy_xx_xxyy_1,   \
                             ta2_xy_xx_xxyz_0,   \
                             ta2_xy_xx_xxyz_1,   \
                             ta2_xy_xx_xxzz_0,   \
                             ta2_xy_xx_xxzz_1,   \
                             ta2_xy_xx_xyyy_0,   \
                             ta2_xy_xx_xyyy_1,   \
                             ta2_xy_xx_xyyz_0,   \
                             ta2_xy_xx_xyyz_1,   \
                             ta2_xy_xx_xyzz_0,   \
                             ta2_xy_xx_xyzz_1,   \
                             ta2_xy_xx_xzzz_0,   \
                             ta2_xy_xx_xzzz_1,   \
                             ta2_xy_xx_yyyy_0,   \
                             ta2_xy_xx_yyyy_1,   \
                             ta2_xy_xxz_xxx_0,   \
                             ta2_xy_xxz_xxx_1,   \
                             ta2_xy_xxz_xxxx_0,  \
                             ta2_xy_xxz_xxxx_1,  \
                             ta2_xy_xxz_xxxy_0,  \
                             ta2_xy_xxz_xxxy_1,  \
                             ta2_xy_xxz_xxxz_0,  \
                             ta2_xy_xxz_xxxz_1,  \
                             ta2_xy_xxz_xxy_0,   \
                             ta2_xy_xxz_xxy_1,   \
                             ta2_xy_xxz_xxyy_0,  \
                             ta2_xy_xxz_xxyy_1,  \
                             ta2_xy_xxz_xxyz_0,  \
                             ta2_xy_xxz_xxyz_1,  \
                             ta2_xy_xxz_xxz_0,   \
                             ta2_xy_xxz_xxz_1,   \
                             ta2_xy_xxz_xxzz_0,  \
                             ta2_xy_xxz_xxzz_1,  \
                             ta2_xy_xxz_xyy_0,   \
                             ta2_xy_xxz_xyy_1,   \
                             ta2_xy_xxz_xyyy_0,  \
                             ta2_xy_xxz_xyyy_1,  \
                             ta2_xy_xxz_xyyz_0,  \
                             ta2_xy_xxz_xyyz_1,  \
                             ta2_xy_xxz_xyz_0,   \
                             ta2_xy_xxz_xyz_1,   \
                             ta2_xy_xxz_xyzz_0,  \
                             ta2_xy_xxz_xyzz_1,  \
                             ta2_xy_xxz_xzz_0,   \
                             ta2_xy_xxz_xzz_1,   \
                             ta2_xy_xxz_xzzz_0,  \
                             ta2_xy_xxz_xzzz_1,  \
                             ta2_xy_xxz_yyyy_0,  \
                             ta2_xy_xxz_yyyy_1,  \
                             ta2_xy_xxzz_xxxx_0, \
                             ta2_xy_xxzz_xxxy_0, \
                             ta2_xy_xxzz_xxxz_0, \
                             ta2_xy_xxzz_xxyy_0, \
                             ta2_xy_xxzz_xxyz_0, \
                             ta2_xy_xxzz_xxzz_0, \
                             ta2_xy_xxzz_xyyy_0, \
                             ta2_xy_xxzz_xyyz_0, \
                             ta2_xy_xxzz_xyzz_0, \
                             ta2_xy_xxzz_xzzz_0, \
                             ta2_xy_xxzz_yyyy_0, \
                             ta2_xy_xxzz_yyyz_0, \
                             ta2_xy_xxzz_yyzz_0, \
                             ta2_xy_xxzz_yzzz_0, \
                             ta2_xy_xxzz_zzzz_0, \
                             ta2_xy_xzz_yyyz_0,  \
                             ta2_xy_xzz_yyyz_1,  \
                             ta2_xy_xzz_yyzz_0,  \
                             ta2_xy_xzz_yyzz_1,  \
                             ta2_xy_xzz_yzzz_0,  \
                             ta2_xy_xzz_yzzz_1,  \
                             ta2_xy_xzz_zzzz_0,  \
                             ta2_xy_xzz_zzzz_1,  \
                             ta2_xy_zz_yyyz_0,   \
                             ta2_xy_zz_yyyz_1,   \
                             ta2_xy_zz_yyzz_0,   \
                             ta2_xy_zz_yyzz_1,   \
                             ta2_xy_zz_yzzz_0,   \
                             ta2_xy_zz_yzzz_1,   \
                             ta2_xy_zz_zzzz_0,   \
                             ta2_xy_zz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xxzz_xxxx_0[i] =
            ta2_xy_xx_xxxx_0[i] * fe_0 - ta2_xy_xx_xxxx_1[i] * fe_0 + ta2_xy_xxz_xxxx_0[i] * pa_z[i] - ta2_xy_xxz_xxxx_1[i] * pc_z[i];

        ta2_xy_xxzz_xxxy_0[i] =
            ta2_xy_xx_xxxy_0[i] * fe_0 - ta2_xy_xx_xxxy_1[i] * fe_0 + ta2_xy_xxz_xxxy_0[i] * pa_z[i] - ta2_xy_xxz_xxxy_1[i] * pc_z[i];

        ta2_xy_xxzz_xxxz_0[i] = ta2_xy_xx_xxxz_0[i] * fe_0 - ta2_xy_xx_xxxz_1[i] * fe_0 + ta2_xy_xxz_xxx_0[i] * fe_0 - ta2_xy_xxz_xxx_1[i] * fe_0 +
                                ta2_xy_xxz_xxxz_0[i] * pa_z[i] - ta2_xy_xxz_xxxz_1[i] * pc_z[i];

        ta2_xy_xxzz_xxyy_0[i] =
            ta2_xy_xx_xxyy_0[i] * fe_0 - ta2_xy_xx_xxyy_1[i] * fe_0 + ta2_xy_xxz_xxyy_0[i] * pa_z[i] - ta2_xy_xxz_xxyy_1[i] * pc_z[i];

        ta2_xy_xxzz_xxyz_0[i] = ta2_xy_xx_xxyz_0[i] * fe_0 - ta2_xy_xx_xxyz_1[i] * fe_0 + ta2_xy_xxz_xxy_0[i] * fe_0 - ta2_xy_xxz_xxy_1[i] * fe_0 +
                                ta2_xy_xxz_xxyz_0[i] * pa_z[i] - ta2_xy_xxz_xxyz_1[i] * pc_z[i];

        ta2_xy_xxzz_xxzz_0[i] = ta2_xy_xx_xxzz_0[i] * fe_0 - ta2_xy_xx_xxzz_1[i] * fe_0 + 2.0 * ta2_xy_xxz_xxz_0[i] * fe_0 -
                                2.0 * ta2_xy_xxz_xxz_1[i] * fe_0 + ta2_xy_xxz_xxzz_0[i] * pa_z[i] - ta2_xy_xxz_xxzz_1[i] * pc_z[i];

        ta2_xy_xxzz_xyyy_0[i] =
            ta2_xy_xx_xyyy_0[i] * fe_0 - ta2_xy_xx_xyyy_1[i] * fe_0 + ta2_xy_xxz_xyyy_0[i] * pa_z[i] - ta2_xy_xxz_xyyy_1[i] * pc_z[i];

        ta2_xy_xxzz_xyyz_0[i] = ta2_xy_xx_xyyz_0[i] * fe_0 - ta2_xy_xx_xyyz_1[i] * fe_0 + ta2_xy_xxz_xyy_0[i] * fe_0 - ta2_xy_xxz_xyy_1[i] * fe_0 +
                                ta2_xy_xxz_xyyz_0[i] * pa_z[i] - ta2_xy_xxz_xyyz_1[i] * pc_z[i];

        ta2_xy_xxzz_xyzz_0[i] = ta2_xy_xx_xyzz_0[i] * fe_0 - ta2_xy_xx_xyzz_1[i] * fe_0 + 2.0 * ta2_xy_xxz_xyz_0[i] * fe_0 -
                                2.0 * ta2_xy_xxz_xyz_1[i] * fe_0 + ta2_xy_xxz_xyzz_0[i] * pa_z[i] - ta2_xy_xxz_xyzz_1[i] * pc_z[i];

        ta2_xy_xxzz_xzzz_0[i] = ta2_xy_xx_xzzz_0[i] * fe_0 - ta2_xy_xx_xzzz_1[i] * fe_0 + 3.0 * ta2_xy_xxz_xzz_0[i] * fe_0 -
                                3.0 * ta2_xy_xxz_xzz_1[i] * fe_0 + ta2_xy_xxz_xzzz_0[i] * pa_z[i] - ta2_xy_xxz_xzzz_1[i] * pc_z[i];

        ta2_xy_xxzz_yyyy_0[i] =
            ta2_xy_xx_yyyy_0[i] * fe_0 - ta2_xy_xx_yyyy_1[i] * fe_0 + ta2_xy_xxz_yyyy_0[i] * pa_z[i] - ta2_xy_xxz_yyyy_1[i] * pc_z[i];

        ta2_xy_xxzz_yyyz_0[i] = ta2_xy_zz_yyyz_0[i] * fe_0 - ta2_xy_zz_yyyz_1[i] * fe_0 + ta1_y_xzz_yyyz_1[i] + ta2_xy_xzz_yyyz_0[i] * pa_x[i] -
                                ta2_xy_xzz_yyyz_1[i] * pc_x[i];

        ta2_xy_xxzz_yyzz_0[i] = ta2_xy_zz_yyzz_0[i] * fe_0 - ta2_xy_zz_yyzz_1[i] * fe_0 + ta1_y_xzz_yyzz_1[i] + ta2_xy_xzz_yyzz_0[i] * pa_x[i] -
                                ta2_xy_xzz_yyzz_1[i] * pc_x[i];

        ta2_xy_xxzz_yzzz_0[i] = ta2_xy_zz_yzzz_0[i] * fe_0 - ta2_xy_zz_yzzz_1[i] * fe_0 + ta1_y_xzz_yzzz_1[i] + ta2_xy_xzz_yzzz_0[i] * pa_x[i] -
                                ta2_xy_xzz_yzzz_1[i] * pc_x[i];

        ta2_xy_xxzz_zzzz_0[i] = ta2_xy_zz_zzzz_0[i] * fe_0 - ta2_xy_zz_zzzz_1[i] * fe_0 + ta1_y_xzz_zzzz_1[i] + ta2_xy_xzz_zzzz_0[i] * pa_x[i] -
                                ta2_xy_xzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 315-330 components of targeted buffer : GG

    auto ta2_xy_xyyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 315);

    auto ta2_xy_xyyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 316);

    auto ta2_xy_xyyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 317);

    auto ta2_xy_xyyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 318);

    auto ta2_xy_xyyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 319);

    auto ta2_xy_xyyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 320);

    auto ta2_xy_xyyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 321);

    auto ta2_xy_xyyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 322);

    auto ta2_xy_xyyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 323);

    auto ta2_xy_xyyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 324);

    auto ta2_xy_xyyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 325);

    auto ta2_xy_xyyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 326);

    auto ta2_xy_xyyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 327);

    auto ta2_xy_xyyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 328);

    auto ta2_xy_xyyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 329);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_yyy_xxxx_1,   \
                             ta1_y_yyy_xxxy_1,   \
                             ta1_y_yyy_xxxz_1,   \
                             ta1_y_yyy_xxyy_1,   \
                             ta1_y_yyy_xxyz_1,   \
                             ta1_y_yyy_xxzz_1,   \
                             ta1_y_yyy_xyyy_1,   \
                             ta1_y_yyy_xyyz_1,   \
                             ta1_y_yyy_xyzz_1,   \
                             ta1_y_yyy_xzzz_1,   \
                             ta1_y_yyy_yyyy_1,   \
                             ta1_y_yyy_yyyz_1,   \
                             ta1_y_yyy_yyzz_1,   \
                             ta1_y_yyy_yzzz_1,   \
                             ta1_y_yyy_zzzz_1,   \
                             ta2_xy_xyyy_xxxx_0, \
                             ta2_xy_xyyy_xxxy_0, \
                             ta2_xy_xyyy_xxxz_0, \
                             ta2_xy_xyyy_xxyy_0, \
                             ta2_xy_xyyy_xxyz_0, \
                             ta2_xy_xyyy_xxzz_0, \
                             ta2_xy_xyyy_xyyy_0, \
                             ta2_xy_xyyy_xyyz_0, \
                             ta2_xy_xyyy_xyzz_0, \
                             ta2_xy_xyyy_xzzz_0, \
                             ta2_xy_xyyy_yyyy_0, \
                             ta2_xy_xyyy_yyyz_0, \
                             ta2_xy_xyyy_yyzz_0, \
                             ta2_xy_xyyy_yzzz_0, \
                             ta2_xy_xyyy_zzzz_0, \
                             ta2_xy_yyy_xxx_0,   \
                             ta2_xy_yyy_xxx_1,   \
                             ta2_xy_yyy_xxxx_0,  \
                             ta2_xy_yyy_xxxx_1,  \
                             ta2_xy_yyy_xxxy_0,  \
                             ta2_xy_yyy_xxxy_1,  \
                             ta2_xy_yyy_xxxz_0,  \
                             ta2_xy_yyy_xxxz_1,  \
                             ta2_xy_yyy_xxy_0,   \
                             ta2_xy_yyy_xxy_1,   \
                             ta2_xy_yyy_xxyy_0,  \
                             ta2_xy_yyy_xxyy_1,  \
                             ta2_xy_yyy_xxyz_0,  \
                             ta2_xy_yyy_xxyz_1,  \
                             ta2_xy_yyy_xxz_0,   \
                             ta2_xy_yyy_xxz_1,   \
                             ta2_xy_yyy_xxzz_0,  \
                             ta2_xy_yyy_xxzz_1,  \
                             ta2_xy_yyy_xyy_0,   \
                             ta2_xy_yyy_xyy_1,   \
                             ta2_xy_yyy_xyyy_0,  \
                             ta2_xy_yyy_xyyy_1,  \
                             ta2_xy_yyy_xyyz_0,  \
                             ta2_xy_yyy_xyyz_1,  \
                             ta2_xy_yyy_xyz_0,   \
                             ta2_xy_yyy_xyz_1,   \
                             ta2_xy_yyy_xyzz_0,  \
                             ta2_xy_yyy_xyzz_1,  \
                             ta2_xy_yyy_xzz_0,   \
                             ta2_xy_yyy_xzz_1,   \
                             ta2_xy_yyy_xzzz_0,  \
                             ta2_xy_yyy_xzzz_1,  \
                             ta2_xy_yyy_yyy_0,   \
                             ta2_xy_yyy_yyy_1,   \
                             ta2_xy_yyy_yyyy_0,  \
                             ta2_xy_yyy_yyyy_1,  \
                             ta2_xy_yyy_yyyz_0,  \
                             ta2_xy_yyy_yyyz_1,  \
                             ta2_xy_yyy_yyz_0,   \
                             ta2_xy_yyy_yyz_1,   \
                             ta2_xy_yyy_yyzz_0,  \
                             ta2_xy_yyy_yyzz_1,  \
                             ta2_xy_yyy_yzz_0,   \
                             ta2_xy_yyy_yzz_1,   \
                             ta2_xy_yyy_yzzz_0,  \
                             ta2_xy_yyy_yzzz_1,  \
                             ta2_xy_yyy_zzz_0,   \
                             ta2_xy_yyy_zzz_1,   \
                             ta2_xy_yyy_zzzz_0,  \
                             ta2_xy_yyy_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xyyy_xxxx_0[i] = 4.0 * ta2_xy_yyy_xxx_0[i] * fe_0 - 4.0 * ta2_xy_yyy_xxx_1[i] * fe_0 + ta1_y_yyy_xxxx_1[i] +
                                ta2_xy_yyy_xxxx_0[i] * pa_x[i] - ta2_xy_yyy_xxxx_1[i] * pc_x[i];

        ta2_xy_xyyy_xxxy_0[i] = 3.0 * ta2_xy_yyy_xxy_0[i] * fe_0 - 3.0 * ta2_xy_yyy_xxy_1[i] * fe_0 + ta1_y_yyy_xxxy_1[i] +
                                ta2_xy_yyy_xxxy_0[i] * pa_x[i] - ta2_xy_yyy_xxxy_1[i] * pc_x[i];

        ta2_xy_xyyy_xxxz_0[i] = 3.0 * ta2_xy_yyy_xxz_0[i] * fe_0 - 3.0 * ta2_xy_yyy_xxz_1[i] * fe_0 + ta1_y_yyy_xxxz_1[i] +
                                ta2_xy_yyy_xxxz_0[i] * pa_x[i] - ta2_xy_yyy_xxxz_1[i] * pc_x[i];

        ta2_xy_xyyy_xxyy_0[i] = 2.0 * ta2_xy_yyy_xyy_0[i] * fe_0 - 2.0 * ta2_xy_yyy_xyy_1[i] * fe_0 + ta1_y_yyy_xxyy_1[i] +
                                ta2_xy_yyy_xxyy_0[i] * pa_x[i] - ta2_xy_yyy_xxyy_1[i] * pc_x[i];

        ta2_xy_xyyy_xxyz_0[i] = 2.0 * ta2_xy_yyy_xyz_0[i] * fe_0 - 2.0 * ta2_xy_yyy_xyz_1[i] * fe_0 + ta1_y_yyy_xxyz_1[i] +
                                ta2_xy_yyy_xxyz_0[i] * pa_x[i] - ta2_xy_yyy_xxyz_1[i] * pc_x[i];

        ta2_xy_xyyy_xxzz_0[i] = 2.0 * ta2_xy_yyy_xzz_0[i] * fe_0 - 2.0 * ta2_xy_yyy_xzz_1[i] * fe_0 + ta1_y_yyy_xxzz_1[i] +
                                ta2_xy_yyy_xxzz_0[i] * pa_x[i] - ta2_xy_yyy_xxzz_1[i] * pc_x[i];

        ta2_xy_xyyy_xyyy_0[i] = ta2_xy_yyy_yyy_0[i] * fe_0 - ta2_xy_yyy_yyy_1[i] * fe_0 + ta1_y_yyy_xyyy_1[i] + ta2_xy_yyy_xyyy_0[i] * pa_x[i] -
                                ta2_xy_yyy_xyyy_1[i] * pc_x[i];

        ta2_xy_xyyy_xyyz_0[i] = ta2_xy_yyy_yyz_0[i] * fe_0 - ta2_xy_yyy_yyz_1[i] * fe_0 + ta1_y_yyy_xyyz_1[i] + ta2_xy_yyy_xyyz_0[i] * pa_x[i] -
                                ta2_xy_yyy_xyyz_1[i] * pc_x[i];

        ta2_xy_xyyy_xyzz_0[i] = ta2_xy_yyy_yzz_0[i] * fe_0 - ta2_xy_yyy_yzz_1[i] * fe_0 + ta1_y_yyy_xyzz_1[i] + ta2_xy_yyy_xyzz_0[i] * pa_x[i] -
                                ta2_xy_yyy_xyzz_1[i] * pc_x[i];

        ta2_xy_xyyy_xzzz_0[i] = ta2_xy_yyy_zzz_0[i] * fe_0 - ta2_xy_yyy_zzz_1[i] * fe_0 + ta1_y_yyy_xzzz_1[i] + ta2_xy_yyy_xzzz_0[i] * pa_x[i] -
                                ta2_xy_yyy_xzzz_1[i] * pc_x[i];

        ta2_xy_xyyy_yyyy_0[i] = ta1_y_yyy_yyyy_1[i] + ta2_xy_yyy_yyyy_0[i] * pa_x[i] - ta2_xy_yyy_yyyy_1[i] * pc_x[i];

        ta2_xy_xyyy_yyyz_0[i] = ta1_y_yyy_yyyz_1[i] + ta2_xy_yyy_yyyz_0[i] * pa_x[i] - ta2_xy_yyy_yyyz_1[i] * pc_x[i];

        ta2_xy_xyyy_yyzz_0[i] = ta1_y_yyy_yyzz_1[i] + ta2_xy_yyy_yyzz_0[i] * pa_x[i] - ta2_xy_yyy_yyzz_1[i] * pc_x[i];

        ta2_xy_xyyy_yzzz_0[i] = ta1_y_yyy_yzzz_1[i] + ta2_xy_yyy_yzzz_0[i] * pa_x[i] - ta2_xy_yyy_yzzz_1[i] * pc_x[i];

        ta2_xy_xyyy_zzzz_0[i] = ta1_y_yyy_zzzz_1[i] + ta2_xy_yyy_zzzz_0[i] * pa_x[i] - ta2_xy_yyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 330-345 components of targeted buffer : GG

    auto ta2_xy_xyyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 330);

    auto ta2_xy_xyyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 331);

    auto ta2_xy_xyyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 332);

    auto ta2_xy_xyyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 333);

    auto ta2_xy_xyyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 334);

    auto ta2_xy_xyyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 335);

    auto ta2_xy_xyyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 336);

    auto ta2_xy_xyyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 337);

    auto ta2_xy_xyyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 338);

    auto ta2_xy_xyyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 339);

    auto ta2_xy_xyyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 340);

    auto ta2_xy_xyyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 341);

    auto ta2_xy_xyyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 342);

    auto ta2_xy_xyyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 343);

    auto ta2_xy_xyyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 344);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_yyz_yyyz_1,   \
                             ta1_y_yyz_yyzz_1,   \
                             ta1_y_yyz_yzzz_1,   \
                             ta1_y_yyz_zzzz_1,   \
                             ta2_xy_xyy_xxx_0,   \
                             ta2_xy_xyy_xxx_1,   \
                             ta2_xy_xyy_xxxx_0,  \
                             ta2_xy_xyy_xxxx_1,  \
                             ta2_xy_xyy_xxxy_0,  \
                             ta2_xy_xyy_xxxy_1,  \
                             ta2_xy_xyy_xxxz_0,  \
                             ta2_xy_xyy_xxxz_1,  \
                             ta2_xy_xyy_xxy_0,   \
                             ta2_xy_xyy_xxy_1,   \
                             ta2_xy_xyy_xxyy_0,  \
                             ta2_xy_xyy_xxyy_1,  \
                             ta2_xy_xyy_xxyz_0,  \
                             ta2_xy_xyy_xxyz_1,  \
                             ta2_xy_xyy_xxz_0,   \
                             ta2_xy_xyy_xxz_1,   \
                             ta2_xy_xyy_xxzz_0,  \
                             ta2_xy_xyy_xxzz_1,  \
                             ta2_xy_xyy_xyy_0,   \
                             ta2_xy_xyy_xyy_1,   \
                             ta2_xy_xyy_xyyy_0,  \
                             ta2_xy_xyy_xyyy_1,  \
                             ta2_xy_xyy_xyyz_0,  \
                             ta2_xy_xyy_xyyz_1,  \
                             ta2_xy_xyy_xyz_0,   \
                             ta2_xy_xyy_xyz_1,   \
                             ta2_xy_xyy_xyzz_0,  \
                             ta2_xy_xyy_xyzz_1,  \
                             ta2_xy_xyy_xzz_0,   \
                             ta2_xy_xyy_xzz_1,   \
                             ta2_xy_xyy_xzzz_0,  \
                             ta2_xy_xyy_xzzz_1,  \
                             ta2_xy_xyy_yyyy_0,  \
                             ta2_xy_xyy_yyyy_1,  \
                             ta2_xy_xyyz_xxxx_0, \
                             ta2_xy_xyyz_xxxy_0, \
                             ta2_xy_xyyz_xxxz_0, \
                             ta2_xy_xyyz_xxyy_0, \
                             ta2_xy_xyyz_xxyz_0, \
                             ta2_xy_xyyz_xxzz_0, \
                             ta2_xy_xyyz_xyyy_0, \
                             ta2_xy_xyyz_xyyz_0, \
                             ta2_xy_xyyz_xyzz_0, \
                             ta2_xy_xyyz_xzzz_0, \
                             ta2_xy_xyyz_yyyy_0, \
                             ta2_xy_xyyz_yyyz_0, \
                             ta2_xy_xyyz_yyzz_0, \
                             ta2_xy_xyyz_yzzz_0, \
                             ta2_xy_xyyz_zzzz_0, \
                             ta2_xy_yyz_yyyz_0,  \
                             ta2_xy_yyz_yyyz_1,  \
                             ta2_xy_yyz_yyzz_0,  \
                             ta2_xy_yyz_yyzz_1,  \
                             ta2_xy_yyz_yzzz_0,  \
                             ta2_xy_yyz_yzzz_1,  \
                             ta2_xy_yyz_zzzz_0,  \
                             ta2_xy_yyz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xyyz_xxxx_0[i] = ta2_xy_xyy_xxxx_0[i] * pa_z[i] - ta2_xy_xyy_xxxx_1[i] * pc_z[i];

        ta2_xy_xyyz_xxxy_0[i] = ta2_xy_xyy_xxxy_0[i] * pa_z[i] - ta2_xy_xyy_xxxy_1[i] * pc_z[i];

        ta2_xy_xyyz_xxxz_0[i] =
            ta2_xy_xyy_xxx_0[i] * fe_0 - ta2_xy_xyy_xxx_1[i] * fe_0 + ta2_xy_xyy_xxxz_0[i] * pa_z[i] - ta2_xy_xyy_xxxz_1[i] * pc_z[i];

        ta2_xy_xyyz_xxyy_0[i] = ta2_xy_xyy_xxyy_0[i] * pa_z[i] - ta2_xy_xyy_xxyy_1[i] * pc_z[i];

        ta2_xy_xyyz_xxyz_0[i] =
            ta2_xy_xyy_xxy_0[i] * fe_0 - ta2_xy_xyy_xxy_1[i] * fe_0 + ta2_xy_xyy_xxyz_0[i] * pa_z[i] - ta2_xy_xyy_xxyz_1[i] * pc_z[i];

        ta2_xy_xyyz_xxzz_0[i] =
            2.0 * ta2_xy_xyy_xxz_0[i] * fe_0 - 2.0 * ta2_xy_xyy_xxz_1[i] * fe_0 + ta2_xy_xyy_xxzz_0[i] * pa_z[i] - ta2_xy_xyy_xxzz_1[i] * pc_z[i];

        ta2_xy_xyyz_xyyy_0[i] = ta2_xy_xyy_xyyy_0[i] * pa_z[i] - ta2_xy_xyy_xyyy_1[i] * pc_z[i];

        ta2_xy_xyyz_xyyz_0[i] =
            ta2_xy_xyy_xyy_0[i] * fe_0 - ta2_xy_xyy_xyy_1[i] * fe_0 + ta2_xy_xyy_xyyz_0[i] * pa_z[i] - ta2_xy_xyy_xyyz_1[i] * pc_z[i];

        ta2_xy_xyyz_xyzz_0[i] =
            2.0 * ta2_xy_xyy_xyz_0[i] * fe_0 - 2.0 * ta2_xy_xyy_xyz_1[i] * fe_0 + ta2_xy_xyy_xyzz_0[i] * pa_z[i] - ta2_xy_xyy_xyzz_1[i] * pc_z[i];

        ta2_xy_xyyz_xzzz_0[i] =
            3.0 * ta2_xy_xyy_xzz_0[i] * fe_0 - 3.0 * ta2_xy_xyy_xzz_1[i] * fe_0 + ta2_xy_xyy_xzzz_0[i] * pa_z[i] - ta2_xy_xyy_xzzz_1[i] * pc_z[i];

        ta2_xy_xyyz_yyyy_0[i] = ta2_xy_xyy_yyyy_0[i] * pa_z[i] - ta2_xy_xyy_yyyy_1[i] * pc_z[i];

        ta2_xy_xyyz_yyyz_0[i] = ta1_y_yyz_yyyz_1[i] + ta2_xy_yyz_yyyz_0[i] * pa_x[i] - ta2_xy_yyz_yyyz_1[i] * pc_x[i];

        ta2_xy_xyyz_yyzz_0[i] = ta1_y_yyz_yyzz_1[i] + ta2_xy_yyz_yyzz_0[i] * pa_x[i] - ta2_xy_yyz_yyzz_1[i] * pc_x[i];

        ta2_xy_xyyz_yzzz_0[i] = ta1_y_yyz_yzzz_1[i] + ta2_xy_yyz_yzzz_0[i] * pa_x[i] - ta2_xy_yyz_yzzz_1[i] * pc_x[i];

        ta2_xy_xyyz_zzzz_0[i] = ta1_y_yyz_zzzz_1[i] + ta2_xy_yyz_zzzz_0[i] * pa_x[i] - ta2_xy_yyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 345-360 components of targeted buffer : GG

    auto ta2_xy_xyzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 345);

    auto ta2_xy_xyzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 346);

    auto ta2_xy_xyzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 347);

    auto ta2_xy_xyzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 348);

    auto ta2_xy_xyzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 349);

    auto ta2_xy_xyzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 350);

    auto ta2_xy_xyzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 351);

    auto ta2_xy_xyzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 352);

    auto ta2_xy_xyzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 353);

    auto ta2_xy_xyzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 354);

    auto ta2_xy_xyzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 355);

    auto ta2_xy_xyzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 356);

    auto ta2_xy_xyzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 357);

    auto ta2_xy_xyzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 358);

    auto ta2_xy_xyzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 359);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xzz_xxxx_1,   \
                             ta1_x_xzz_xxxz_1,   \
                             ta1_x_xzz_xxzz_1,   \
                             ta1_x_xzz_xzzz_1,   \
                             ta1_y_yzz_xxyz_1,   \
                             ta1_y_yzz_xyyz_1,   \
                             ta1_y_yzz_xyzz_1,   \
                             ta1_y_yzz_yyyy_1,   \
                             ta1_y_yzz_yyyz_1,   \
                             ta1_y_yzz_yyzz_1,   \
                             ta1_y_yzz_yzzz_1,   \
                             ta1_y_yzz_zzzz_1,   \
                             ta2_xy_xy_xxxy_0,   \
                             ta2_xy_xy_xxxy_1,   \
                             ta2_xy_xy_xxyy_0,   \
                             ta2_xy_xy_xxyy_1,   \
                             ta2_xy_xy_xyyy_0,   \
                             ta2_xy_xy_xyyy_1,   \
                             ta2_xy_xyz_xxxy_0,  \
                             ta2_xy_xyz_xxxy_1,  \
                             ta2_xy_xyz_xxyy_0,  \
                             ta2_xy_xyz_xxyy_1,  \
                             ta2_xy_xyz_xyyy_0,  \
                             ta2_xy_xyz_xyyy_1,  \
                             ta2_xy_xyzz_xxxx_0, \
                             ta2_xy_xyzz_xxxy_0, \
                             ta2_xy_xyzz_xxxz_0, \
                             ta2_xy_xyzz_xxyy_0, \
                             ta2_xy_xyzz_xxyz_0, \
                             ta2_xy_xyzz_xxzz_0, \
                             ta2_xy_xyzz_xyyy_0, \
                             ta2_xy_xyzz_xyyz_0, \
                             ta2_xy_xyzz_xyzz_0, \
                             ta2_xy_xyzz_xzzz_0, \
                             ta2_xy_xyzz_yyyy_0, \
                             ta2_xy_xyzz_yyyz_0, \
                             ta2_xy_xyzz_yyzz_0, \
                             ta2_xy_xyzz_yzzz_0, \
                             ta2_xy_xyzz_zzzz_0, \
                             ta2_xy_xzz_xxxx_0,  \
                             ta2_xy_xzz_xxxx_1,  \
                             ta2_xy_xzz_xxxz_0,  \
                             ta2_xy_xzz_xxxz_1,  \
                             ta2_xy_xzz_xxzz_0,  \
                             ta2_xy_xzz_xxzz_1,  \
                             ta2_xy_xzz_xzzz_0,  \
                             ta2_xy_xzz_xzzz_1,  \
                             ta2_xy_yzz_xxyz_0,  \
                             ta2_xy_yzz_xxyz_1,  \
                             ta2_xy_yzz_xyyz_0,  \
                             ta2_xy_yzz_xyyz_1,  \
                             ta2_xy_yzz_xyz_0,   \
                             ta2_xy_yzz_xyz_1,   \
                             ta2_xy_yzz_xyzz_0,  \
                             ta2_xy_yzz_xyzz_1,  \
                             ta2_xy_yzz_yyyy_0,  \
                             ta2_xy_yzz_yyyy_1,  \
                             ta2_xy_yzz_yyyz_0,  \
                             ta2_xy_yzz_yyyz_1,  \
                             ta2_xy_yzz_yyz_0,   \
                             ta2_xy_yzz_yyz_1,   \
                             ta2_xy_yzz_yyzz_0,  \
                             ta2_xy_yzz_yyzz_1,  \
                             ta2_xy_yzz_yzz_0,   \
                             ta2_xy_yzz_yzz_1,   \
                             ta2_xy_yzz_yzzz_0,  \
                             ta2_xy_yzz_yzzz_1,  \
                             ta2_xy_yzz_zzzz_0,  \
                             ta2_xy_yzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xyzz_xxxx_0[i] = ta1_x_xzz_xxxx_1[i] + ta2_xy_xzz_xxxx_0[i] * pa_y[i] - ta2_xy_xzz_xxxx_1[i] * pc_y[i];

        ta2_xy_xyzz_xxxy_0[i] =
            ta2_xy_xy_xxxy_0[i] * fe_0 - ta2_xy_xy_xxxy_1[i] * fe_0 + ta2_xy_xyz_xxxy_0[i] * pa_z[i] - ta2_xy_xyz_xxxy_1[i] * pc_z[i];

        ta2_xy_xyzz_xxxz_0[i] = ta1_x_xzz_xxxz_1[i] + ta2_xy_xzz_xxxz_0[i] * pa_y[i] - ta2_xy_xzz_xxxz_1[i] * pc_y[i];

        ta2_xy_xyzz_xxyy_0[i] =
            ta2_xy_xy_xxyy_0[i] * fe_0 - ta2_xy_xy_xxyy_1[i] * fe_0 + ta2_xy_xyz_xxyy_0[i] * pa_z[i] - ta2_xy_xyz_xxyy_1[i] * pc_z[i];

        ta2_xy_xyzz_xxyz_0[i] = 2.0 * ta2_xy_yzz_xyz_0[i] * fe_0 - 2.0 * ta2_xy_yzz_xyz_1[i] * fe_0 + ta1_y_yzz_xxyz_1[i] +
                                ta2_xy_yzz_xxyz_0[i] * pa_x[i] - ta2_xy_yzz_xxyz_1[i] * pc_x[i];

        ta2_xy_xyzz_xxzz_0[i] = ta1_x_xzz_xxzz_1[i] + ta2_xy_xzz_xxzz_0[i] * pa_y[i] - ta2_xy_xzz_xxzz_1[i] * pc_y[i];

        ta2_xy_xyzz_xyyy_0[i] =
            ta2_xy_xy_xyyy_0[i] * fe_0 - ta2_xy_xy_xyyy_1[i] * fe_0 + ta2_xy_xyz_xyyy_0[i] * pa_z[i] - ta2_xy_xyz_xyyy_1[i] * pc_z[i];

        ta2_xy_xyzz_xyyz_0[i] = ta2_xy_yzz_yyz_0[i] * fe_0 - ta2_xy_yzz_yyz_1[i] * fe_0 + ta1_y_yzz_xyyz_1[i] + ta2_xy_yzz_xyyz_0[i] * pa_x[i] -
                                ta2_xy_yzz_xyyz_1[i] * pc_x[i];

        ta2_xy_xyzz_xyzz_0[i] = ta2_xy_yzz_yzz_0[i] * fe_0 - ta2_xy_yzz_yzz_1[i] * fe_0 + ta1_y_yzz_xyzz_1[i] + ta2_xy_yzz_xyzz_0[i] * pa_x[i] -
                                ta2_xy_yzz_xyzz_1[i] * pc_x[i];

        ta2_xy_xyzz_xzzz_0[i] = ta1_x_xzz_xzzz_1[i] + ta2_xy_xzz_xzzz_0[i] * pa_y[i] - ta2_xy_xzz_xzzz_1[i] * pc_y[i];

        ta2_xy_xyzz_yyyy_0[i] = ta1_y_yzz_yyyy_1[i] + ta2_xy_yzz_yyyy_0[i] * pa_x[i] - ta2_xy_yzz_yyyy_1[i] * pc_x[i];

        ta2_xy_xyzz_yyyz_0[i] = ta1_y_yzz_yyyz_1[i] + ta2_xy_yzz_yyyz_0[i] * pa_x[i] - ta2_xy_yzz_yyyz_1[i] * pc_x[i];

        ta2_xy_xyzz_yyzz_0[i] = ta1_y_yzz_yyzz_1[i] + ta2_xy_yzz_yyzz_0[i] * pa_x[i] - ta2_xy_yzz_yyzz_1[i] * pc_x[i];

        ta2_xy_xyzz_yzzz_0[i] = ta1_y_yzz_yzzz_1[i] + ta2_xy_yzz_yzzz_0[i] * pa_x[i] - ta2_xy_yzz_yzzz_1[i] * pc_x[i];

        ta2_xy_xyzz_zzzz_0[i] = ta1_y_yzz_zzzz_1[i] + ta2_xy_yzz_zzzz_0[i] * pa_x[i] - ta2_xy_yzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 360-375 components of targeted buffer : GG

    auto ta2_xy_xzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 360);

    auto ta2_xy_xzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 361);

    auto ta2_xy_xzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 362);

    auto ta2_xy_xzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 363);

    auto ta2_xy_xzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 364);

    auto ta2_xy_xzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 365);

    auto ta2_xy_xzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 366);

    auto ta2_xy_xzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 367);

    auto ta2_xy_xzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 368);

    auto ta2_xy_xzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 369);

    auto ta2_xy_xzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 370);

    auto ta2_xy_xzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 371);

    auto ta2_xy_xzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 372);

    auto ta2_xy_xzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 373);

    auto ta2_xy_xzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 374);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_zzz_xxxz_1,   \
                             ta1_y_zzz_xxyz_1,   \
                             ta1_y_zzz_xxzz_1,   \
                             ta1_y_zzz_xyyz_1,   \
                             ta1_y_zzz_xyzz_1,   \
                             ta1_y_zzz_xzzz_1,   \
                             ta1_y_zzz_yyyy_1,   \
                             ta1_y_zzz_yyyz_1,   \
                             ta1_y_zzz_yyzz_1,   \
                             ta1_y_zzz_yzzz_1,   \
                             ta1_y_zzz_zzzz_1,   \
                             ta2_xy_xz_xxxx_0,   \
                             ta2_xy_xz_xxxx_1,   \
                             ta2_xy_xz_xxxy_0,   \
                             ta2_xy_xz_xxxy_1,   \
                             ta2_xy_xz_xxyy_0,   \
                             ta2_xy_xz_xxyy_1,   \
                             ta2_xy_xz_xyyy_0,   \
                             ta2_xy_xz_xyyy_1,   \
                             ta2_xy_xzz_xxxx_0,  \
                             ta2_xy_xzz_xxxx_1,  \
                             ta2_xy_xzz_xxxy_0,  \
                             ta2_xy_xzz_xxxy_1,  \
                             ta2_xy_xzz_xxyy_0,  \
                             ta2_xy_xzz_xxyy_1,  \
                             ta2_xy_xzz_xyyy_0,  \
                             ta2_xy_xzz_xyyy_1,  \
                             ta2_xy_xzzz_xxxx_0, \
                             ta2_xy_xzzz_xxxy_0, \
                             ta2_xy_xzzz_xxxz_0, \
                             ta2_xy_xzzz_xxyy_0, \
                             ta2_xy_xzzz_xxyz_0, \
                             ta2_xy_xzzz_xxzz_0, \
                             ta2_xy_xzzz_xyyy_0, \
                             ta2_xy_xzzz_xyyz_0, \
                             ta2_xy_xzzz_xyzz_0, \
                             ta2_xy_xzzz_xzzz_0, \
                             ta2_xy_xzzz_yyyy_0, \
                             ta2_xy_xzzz_yyyz_0, \
                             ta2_xy_xzzz_yyzz_0, \
                             ta2_xy_xzzz_yzzz_0, \
                             ta2_xy_xzzz_zzzz_0, \
                             ta2_xy_zzz_xxxz_0,  \
                             ta2_xy_zzz_xxxz_1,  \
                             ta2_xy_zzz_xxyz_0,  \
                             ta2_xy_zzz_xxyz_1,  \
                             ta2_xy_zzz_xxz_0,   \
                             ta2_xy_zzz_xxz_1,   \
                             ta2_xy_zzz_xxzz_0,  \
                             ta2_xy_zzz_xxzz_1,  \
                             ta2_xy_zzz_xyyz_0,  \
                             ta2_xy_zzz_xyyz_1,  \
                             ta2_xy_zzz_xyz_0,   \
                             ta2_xy_zzz_xyz_1,   \
                             ta2_xy_zzz_xyzz_0,  \
                             ta2_xy_zzz_xyzz_1,  \
                             ta2_xy_zzz_xzz_0,   \
                             ta2_xy_zzz_xzz_1,   \
                             ta2_xy_zzz_xzzz_0,  \
                             ta2_xy_zzz_xzzz_1,  \
                             ta2_xy_zzz_yyyy_0,  \
                             ta2_xy_zzz_yyyy_1,  \
                             ta2_xy_zzz_yyyz_0,  \
                             ta2_xy_zzz_yyyz_1,  \
                             ta2_xy_zzz_yyz_0,   \
                             ta2_xy_zzz_yyz_1,   \
                             ta2_xy_zzz_yyzz_0,  \
                             ta2_xy_zzz_yyzz_1,  \
                             ta2_xy_zzz_yzz_0,   \
                             ta2_xy_zzz_yzz_1,   \
                             ta2_xy_zzz_yzzz_0,  \
                             ta2_xy_zzz_yzzz_1,  \
                             ta2_xy_zzz_zzz_0,   \
                             ta2_xy_zzz_zzz_1,   \
                             ta2_xy_zzz_zzzz_0,  \
                             ta2_xy_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_xzzz_xxxx_0[i] =
            2.0 * ta2_xy_xz_xxxx_0[i] * fe_0 - 2.0 * ta2_xy_xz_xxxx_1[i] * fe_0 + ta2_xy_xzz_xxxx_0[i] * pa_z[i] - ta2_xy_xzz_xxxx_1[i] * pc_z[i];

        ta2_xy_xzzz_xxxy_0[i] =
            2.0 * ta2_xy_xz_xxxy_0[i] * fe_0 - 2.0 * ta2_xy_xz_xxxy_1[i] * fe_0 + ta2_xy_xzz_xxxy_0[i] * pa_z[i] - ta2_xy_xzz_xxxy_1[i] * pc_z[i];

        ta2_xy_xzzz_xxxz_0[i] = 3.0 * ta2_xy_zzz_xxz_0[i] * fe_0 - 3.0 * ta2_xy_zzz_xxz_1[i] * fe_0 + ta1_y_zzz_xxxz_1[i] +
                                ta2_xy_zzz_xxxz_0[i] * pa_x[i] - ta2_xy_zzz_xxxz_1[i] * pc_x[i];

        ta2_xy_xzzz_xxyy_0[i] =
            2.0 * ta2_xy_xz_xxyy_0[i] * fe_0 - 2.0 * ta2_xy_xz_xxyy_1[i] * fe_0 + ta2_xy_xzz_xxyy_0[i] * pa_z[i] - ta2_xy_xzz_xxyy_1[i] * pc_z[i];

        ta2_xy_xzzz_xxyz_0[i] = 2.0 * ta2_xy_zzz_xyz_0[i] * fe_0 - 2.0 * ta2_xy_zzz_xyz_1[i] * fe_0 + ta1_y_zzz_xxyz_1[i] +
                                ta2_xy_zzz_xxyz_0[i] * pa_x[i] - ta2_xy_zzz_xxyz_1[i] * pc_x[i];

        ta2_xy_xzzz_xxzz_0[i] = 2.0 * ta2_xy_zzz_xzz_0[i] * fe_0 - 2.0 * ta2_xy_zzz_xzz_1[i] * fe_0 + ta1_y_zzz_xxzz_1[i] +
                                ta2_xy_zzz_xxzz_0[i] * pa_x[i] - ta2_xy_zzz_xxzz_1[i] * pc_x[i];

        ta2_xy_xzzz_xyyy_0[i] =
            2.0 * ta2_xy_xz_xyyy_0[i] * fe_0 - 2.0 * ta2_xy_xz_xyyy_1[i] * fe_0 + ta2_xy_xzz_xyyy_0[i] * pa_z[i] - ta2_xy_xzz_xyyy_1[i] * pc_z[i];

        ta2_xy_xzzz_xyyz_0[i] = ta2_xy_zzz_yyz_0[i] * fe_0 - ta2_xy_zzz_yyz_1[i] * fe_0 + ta1_y_zzz_xyyz_1[i] + ta2_xy_zzz_xyyz_0[i] * pa_x[i] -
                                ta2_xy_zzz_xyyz_1[i] * pc_x[i];

        ta2_xy_xzzz_xyzz_0[i] = ta2_xy_zzz_yzz_0[i] * fe_0 - ta2_xy_zzz_yzz_1[i] * fe_0 + ta1_y_zzz_xyzz_1[i] + ta2_xy_zzz_xyzz_0[i] * pa_x[i] -
                                ta2_xy_zzz_xyzz_1[i] * pc_x[i];

        ta2_xy_xzzz_xzzz_0[i] = ta2_xy_zzz_zzz_0[i] * fe_0 - ta2_xy_zzz_zzz_1[i] * fe_0 + ta1_y_zzz_xzzz_1[i] + ta2_xy_zzz_xzzz_0[i] * pa_x[i] -
                                ta2_xy_zzz_xzzz_1[i] * pc_x[i];

        ta2_xy_xzzz_yyyy_0[i] = ta1_y_zzz_yyyy_1[i] + ta2_xy_zzz_yyyy_0[i] * pa_x[i] - ta2_xy_zzz_yyyy_1[i] * pc_x[i];

        ta2_xy_xzzz_yyyz_0[i] = ta1_y_zzz_yyyz_1[i] + ta2_xy_zzz_yyyz_0[i] * pa_x[i] - ta2_xy_zzz_yyyz_1[i] * pc_x[i];

        ta2_xy_xzzz_yyzz_0[i] = ta1_y_zzz_yyzz_1[i] + ta2_xy_zzz_yyzz_0[i] * pa_x[i] - ta2_xy_zzz_yyzz_1[i] * pc_x[i];

        ta2_xy_xzzz_yzzz_0[i] = ta1_y_zzz_yzzz_1[i] + ta2_xy_zzz_yzzz_0[i] * pa_x[i] - ta2_xy_zzz_yzzz_1[i] * pc_x[i];

        ta2_xy_xzzz_zzzz_0[i] = ta1_y_zzz_zzzz_1[i] + ta2_xy_zzz_zzzz_0[i] * pa_x[i] - ta2_xy_zzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 375-390 components of targeted buffer : GG

    auto ta2_xy_yyyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 375);

    auto ta2_xy_yyyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 376);

    auto ta2_xy_yyyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 377);

    auto ta2_xy_yyyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 378);

    auto ta2_xy_yyyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 379);

    auto ta2_xy_yyyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 380);

    auto ta2_xy_yyyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 381);

    auto ta2_xy_yyyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 382);

    auto ta2_xy_yyyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 383);

    auto ta2_xy_yyyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 384);

    auto ta2_xy_yyyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 385);

    auto ta2_xy_yyyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 386);

    auto ta2_xy_yyyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 387);

    auto ta2_xy_yyyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 388);

    auto ta2_xy_yyyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 389);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_x_yyy_xxxx_1,   \
                             ta1_x_yyy_xxxy_1,   \
                             ta1_x_yyy_xxxz_1,   \
                             ta1_x_yyy_xxyy_1,   \
                             ta1_x_yyy_xxyz_1,   \
                             ta1_x_yyy_xxzz_1,   \
                             ta1_x_yyy_xyyy_1,   \
                             ta1_x_yyy_xyyz_1,   \
                             ta1_x_yyy_xyzz_1,   \
                             ta1_x_yyy_xzzz_1,   \
                             ta1_x_yyy_yyyy_1,   \
                             ta1_x_yyy_yyyz_1,   \
                             ta1_x_yyy_yyzz_1,   \
                             ta1_x_yyy_yzzz_1,   \
                             ta1_x_yyy_zzzz_1,   \
                             ta2_xy_yy_xxxx_0,   \
                             ta2_xy_yy_xxxx_1,   \
                             ta2_xy_yy_xxxy_0,   \
                             ta2_xy_yy_xxxy_1,   \
                             ta2_xy_yy_xxxz_0,   \
                             ta2_xy_yy_xxxz_1,   \
                             ta2_xy_yy_xxyy_0,   \
                             ta2_xy_yy_xxyy_1,   \
                             ta2_xy_yy_xxyz_0,   \
                             ta2_xy_yy_xxyz_1,   \
                             ta2_xy_yy_xxzz_0,   \
                             ta2_xy_yy_xxzz_1,   \
                             ta2_xy_yy_xyyy_0,   \
                             ta2_xy_yy_xyyy_1,   \
                             ta2_xy_yy_xyyz_0,   \
                             ta2_xy_yy_xyyz_1,   \
                             ta2_xy_yy_xyzz_0,   \
                             ta2_xy_yy_xyzz_1,   \
                             ta2_xy_yy_xzzz_0,   \
                             ta2_xy_yy_xzzz_1,   \
                             ta2_xy_yy_yyyy_0,   \
                             ta2_xy_yy_yyyy_1,   \
                             ta2_xy_yy_yyyz_0,   \
                             ta2_xy_yy_yyyz_1,   \
                             ta2_xy_yy_yyzz_0,   \
                             ta2_xy_yy_yyzz_1,   \
                             ta2_xy_yy_yzzz_0,   \
                             ta2_xy_yy_yzzz_1,   \
                             ta2_xy_yy_zzzz_0,   \
                             ta2_xy_yy_zzzz_1,   \
                             ta2_xy_yyy_xxx_0,   \
                             ta2_xy_yyy_xxx_1,   \
                             ta2_xy_yyy_xxxx_0,  \
                             ta2_xy_yyy_xxxx_1,  \
                             ta2_xy_yyy_xxxy_0,  \
                             ta2_xy_yyy_xxxy_1,  \
                             ta2_xy_yyy_xxxz_0,  \
                             ta2_xy_yyy_xxxz_1,  \
                             ta2_xy_yyy_xxy_0,   \
                             ta2_xy_yyy_xxy_1,   \
                             ta2_xy_yyy_xxyy_0,  \
                             ta2_xy_yyy_xxyy_1,  \
                             ta2_xy_yyy_xxyz_0,  \
                             ta2_xy_yyy_xxyz_1,  \
                             ta2_xy_yyy_xxz_0,   \
                             ta2_xy_yyy_xxz_1,   \
                             ta2_xy_yyy_xxzz_0,  \
                             ta2_xy_yyy_xxzz_1,  \
                             ta2_xy_yyy_xyy_0,   \
                             ta2_xy_yyy_xyy_1,   \
                             ta2_xy_yyy_xyyy_0,  \
                             ta2_xy_yyy_xyyy_1,  \
                             ta2_xy_yyy_xyyz_0,  \
                             ta2_xy_yyy_xyyz_1,  \
                             ta2_xy_yyy_xyz_0,   \
                             ta2_xy_yyy_xyz_1,   \
                             ta2_xy_yyy_xyzz_0,  \
                             ta2_xy_yyy_xyzz_1,  \
                             ta2_xy_yyy_xzz_0,   \
                             ta2_xy_yyy_xzz_1,   \
                             ta2_xy_yyy_xzzz_0,  \
                             ta2_xy_yyy_xzzz_1,  \
                             ta2_xy_yyy_yyy_0,   \
                             ta2_xy_yyy_yyy_1,   \
                             ta2_xy_yyy_yyyy_0,  \
                             ta2_xy_yyy_yyyy_1,  \
                             ta2_xy_yyy_yyyz_0,  \
                             ta2_xy_yyy_yyyz_1,  \
                             ta2_xy_yyy_yyz_0,   \
                             ta2_xy_yyy_yyz_1,   \
                             ta2_xy_yyy_yyzz_0,  \
                             ta2_xy_yyy_yyzz_1,  \
                             ta2_xy_yyy_yzz_0,   \
                             ta2_xy_yyy_yzz_1,   \
                             ta2_xy_yyy_yzzz_0,  \
                             ta2_xy_yyy_yzzz_1,  \
                             ta2_xy_yyy_zzz_0,   \
                             ta2_xy_yyy_zzz_1,   \
                             ta2_xy_yyy_zzzz_0,  \
                             ta2_xy_yyy_zzzz_1,  \
                             ta2_xy_yyyy_xxxx_0, \
                             ta2_xy_yyyy_xxxy_0, \
                             ta2_xy_yyyy_xxxz_0, \
                             ta2_xy_yyyy_xxyy_0, \
                             ta2_xy_yyyy_xxyz_0, \
                             ta2_xy_yyyy_xxzz_0, \
                             ta2_xy_yyyy_xyyy_0, \
                             ta2_xy_yyyy_xyyz_0, \
                             ta2_xy_yyyy_xyzz_0, \
                             ta2_xy_yyyy_xzzz_0, \
                             ta2_xy_yyyy_yyyy_0, \
                             ta2_xy_yyyy_yyyz_0, \
                             ta2_xy_yyyy_yyzz_0, \
                             ta2_xy_yyyy_yzzz_0, \
                             ta2_xy_yyyy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyyy_xxxx_0[i] = 3.0 * ta2_xy_yy_xxxx_0[i] * fe_0 - 3.0 * ta2_xy_yy_xxxx_1[i] * fe_0 + ta1_x_yyy_xxxx_1[i] +
                                ta2_xy_yyy_xxxx_0[i] * pa_y[i] - ta2_xy_yyy_xxxx_1[i] * pc_y[i];

        ta2_xy_yyyy_xxxy_0[i] = 3.0 * ta2_xy_yy_xxxy_0[i] * fe_0 - 3.0 * ta2_xy_yy_xxxy_1[i] * fe_0 + ta2_xy_yyy_xxx_0[i] * fe_0 -
                                ta2_xy_yyy_xxx_1[i] * fe_0 + ta1_x_yyy_xxxy_1[i] + ta2_xy_yyy_xxxy_0[i] * pa_y[i] - ta2_xy_yyy_xxxy_1[i] * pc_y[i];

        ta2_xy_yyyy_xxxz_0[i] = 3.0 * ta2_xy_yy_xxxz_0[i] * fe_0 - 3.0 * ta2_xy_yy_xxxz_1[i] * fe_0 + ta1_x_yyy_xxxz_1[i] +
                                ta2_xy_yyy_xxxz_0[i] * pa_y[i] - ta2_xy_yyy_xxxz_1[i] * pc_y[i];

        ta2_xy_yyyy_xxyy_0[i] = 3.0 * ta2_xy_yy_xxyy_0[i] * fe_0 - 3.0 * ta2_xy_yy_xxyy_1[i] * fe_0 + 2.0 * ta2_xy_yyy_xxy_0[i] * fe_0 -
                                2.0 * ta2_xy_yyy_xxy_1[i] * fe_0 + ta1_x_yyy_xxyy_1[i] + ta2_xy_yyy_xxyy_0[i] * pa_y[i] -
                                ta2_xy_yyy_xxyy_1[i] * pc_y[i];

        ta2_xy_yyyy_xxyz_0[i] = 3.0 * ta2_xy_yy_xxyz_0[i] * fe_0 - 3.0 * ta2_xy_yy_xxyz_1[i] * fe_0 + ta2_xy_yyy_xxz_0[i] * fe_0 -
                                ta2_xy_yyy_xxz_1[i] * fe_0 + ta1_x_yyy_xxyz_1[i] + ta2_xy_yyy_xxyz_0[i] * pa_y[i] - ta2_xy_yyy_xxyz_1[i] * pc_y[i];

        ta2_xy_yyyy_xxzz_0[i] = 3.0 * ta2_xy_yy_xxzz_0[i] * fe_0 - 3.0 * ta2_xy_yy_xxzz_1[i] * fe_0 + ta1_x_yyy_xxzz_1[i] +
                                ta2_xy_yyy_xxzz_0[i] * pa_y[i] - ta2_xy_yyy_xxzz_1[i] * pc_y[i];

        ta2_xy_yyyy_xyyy_0[i] = 3.0 * ta2_xy_yy_xyyy_0[i] * fe_0 - 3.0 * ta2_xy_yy_xyyy_1[i] * fe_0 + 3.0 * ta2_xy_yyy_xyy_0[i] * fe_0 -
                                3.0 * ta2_xy_yyy_xyy_1[i] * fe_0 + ta1_x_yyy_xyyy_1[i] + ta2_xy_yyy_xyyy_0[i] * pa_y[i] -
                                ta2_xy_yyy_xyyy_1[i] * pc_y[i];

        ta2_xy_yyyy_xyyz_0[i] = 3.0 * ta2_xy_yy_xyyz_0[i] * fe_0 - 3.0 * ta2_xy_yy_xyyz_1[i] * fe_0 + 2.0 * ta2_xy_yyy_xyz_0[i] * fe_0 -
                                2.0 * ta2_xy_yyy_xyz_1[i] * fe_0 + ta1_x_yyy_xyyz_1[i] + ta2_xy_yyy_xyyz_0[i] * pa_y[i] -
                                ta2_xy_yyy_xyyz_1[i] * pc_y[i];

        ta2_xy_yyyy_xyzz_0[i] = 3.0 * ta2_xy_yy_xyzz_0[i] * fe_0 - 3.0 * ta2_xy_yy_xyzz_1[i] * fe_0 + ta2_xy_yyy_xzz_0[i] * fe_0 -
                                ta2_xy_yyy_xzz_1[i] * fe_0 + ta1_x_yyy_xyzz_1[i] + ta2_xy_yyy_xyzz_0[i] * pa_y[i] - ta2_xy_yyy_xyzz_1[i] * pc_y[i];

        ta2_xy_yyyy_xzzz_0[i] = 3.0 * ta2_xy_yy_xzzz_0[i] * fe_0 - 3.0 * ta2_xy_yy_xzzz_1[i] * fe_0 + ta1_x_yyy_xzzz_1[i] +
                                ta2_xy_yyy_xzzz_0[i] * pa_y[i] - ta2_xy_yyy_xzzz_1[i] * pc_y[i];

        ta2_xy_yyyy_yyyy_0[i] = 3.0 * ta2_xy_yy_yyyy_0[i] * fe_0 - 3.0 * ta2_xy_yy_yyyy_1[i] * fe_0 + 4.0 * ta2_xy_yyy_yyy_0[i] * fe_0 -
                                4.0 * ta2_xy_yyy_yyy_1[i] * fe_0 + ta1_x_yyy_yyyy_1[i] + ta2_xy_yyy_yyyy_0[i] * pa_y[i] -
                                ta2_xy_yyy_yyyy_1[i] * pc_y[i];

        ta2_xy_yyyy_yyyz_0[i] = 3.0 * ta2_xy_yy_yyyz_0[i] * fe_0 - 3.0 * ta2_xy_yy_yyyz_1[i] * fe_0 + 3.0 * ta2_xy_yyy_yyz_0[i] * fe_0 -
                                3.0 * ta2_xy_yyy_yyz_1[i] * fe_0 + ta1_x_yyy_yyyz_1[i] + ta2_xy_yyy_yyyz_0[i] * pa_y[i] -
                                ta2_xy_yyy_yyyz_1[i] * pc_y[i];

        ta2_xy_yyyy_yyzz_0[i] = 3.0 * ta2_xy_yy_yyzz_0[i] * fe_0 - 3.0 * ta2_xy_yy_yyzz_1[i] * fe_0 + 2.0 * ta2_xy_yyy_yzz_0[i] * fe_0 -
                                2.0 * ta2_xy_yyy_yzz_1[i] * fe_0 + ta1_x_yyy_yyzz_1[i] + ta2_xy_yyy_yyzz_0[i] * pa_y[i] -
                                ta2_xy_yyy_yyzz_1[i] * pc_y[i];

        ta2_xy_yyyy_yzzz_0[i] = 3.0 * ta2_xy_yy_yzzz_0[i] * fe_0 - 3.0 * ta2_xy_yy_yzzz_1[i] * fe_0 + ta2_xy_yyy_zzz_0[i] * fe_0 -
                                ta2_xy_yyy_zzz_1[i] * fe_0 + ta1_x_yyy_yzzz_1[i] + ta2_xy_yyy_yzzz_0[i] * pa_y[i] - ta2_xy_yyy_yzzz_1[i] * pc_y[i];

        ta2_xy_yyyy_zzzz_0[i] = 3.0 * ta2_xy_yy_zzzz_0[i] * fe_0 - 3.0 * ta2_xy_yy_zzzz_1[i] * fe_0 + ta1_x_yyy_zzzz_1[i] +
                                ta2_xy_yyy_zzzz_0[i] * pa_y[i] - ta2_xy_yyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 390-405 components of targeted buffer : GG

    auto ta2_xy_yyyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 390);

    auto ta2_xy_yyyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 391);

    auto ta2_xy_yyyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 392);

    auto ta2_xy_yyyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 393);

    auto ta2_xy_yyyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 394);

    auto ta2_xy_yyyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 395);

    auto ta2_xy_yyyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 396);

    auto ta2_xy_yyyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 397);

    auto ta2_xy_yyyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 398);

    auto ta2_xy_yyyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 399);

    auto ta2_xy_yyyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 400);

    auto ta2_xy_yyyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 401);

    auto ta2_xy_yyyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 402);

    auto ta2_xy_yyyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 403);

    auto ta2_xy_yyyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 404);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta2_xy_yyy_xxx_0,   \
                             ta2_xy_yyy_xxx_1,   \
                             ta2_xy_yyy_xxxx_0,  \
                             ta2_xy_yyy_xxxx_1,  \
                             ta2_xy_yyy_xxxy_0,  \
                             ta2_xy_yyy_xxxy_1,  \
                             ta2_xy_yyy_xxxz_0,  \
                             ta2_xy_yyy_xxxz_1,  \
                             ta2_xy_yyy_xxy_0,   \
                             ta2_xy_yyy_xxy_1,   \
                             ta2_xy_yyy_xxyy_0,  \
                             ta2_xy_yyy_xxyy_1,  \
                             ta2_xy_yyy_xxyz_0,  \
                             ta2_xy_yyy_xxyz_1,  \
                             ta2_xy_yyy_xxz_0,   \
                             ta2_xy_yyy_xxz_1,   \
                             ta2_xy_yyy_xxzz_0,  \
                             ta2_xy_yyy_xxzz_1,  \
                             ta2_xy_yyy_xyy_0,   \
                             ta2_xy_yyy_xyy_1,   \
                             ta2_xy_yyy_xyyy_0,  \
                             ta2_xy_yyy_xyyy_1,  \
                             ta2_xy_yyy_xyyz_0,  \
                             ta2_xy_yyy_xyyz_1,  \
                             ta2_xy_yyy_xyz_0,   \
                             ta2_xy_yyy_xyz_1,   \
                             ta2_xy_yyy_xyzz_0,  \
                             ta2_xy_yyy_xyzz_1,  \
                             ta2_xy_yyy_xzz_0,   \
                             ta2_xy_yyy_xzz_1,   \
                             ta2_xy_yyy_xzzz_0,  \
                             ta2_xy_yyy_xzzz_1,  \
                             ta2_xy_yyy_yyy_0,   \
                             ta2_xy_yyy_yyy_1,   \
                             ta2_xy_yyy_yyyy_0,  \
                             ta2_xy_yyy_yyyy_1,  \
                             ta2_xy_yyy_yyyz_0,  \
                             ta2_xy_yyy_yyyz_1,  \
                             ta2_xy_yyy_yyz_0,   \
                             ta2_xy_yyy_yyz_1,   \
                             ta2_xy_yyy_yyzz_0,  \
                             ta2_xy_yyy_yyzz_1,  \
                             ta2_xy_yyy_yzz_0,   \
                             ta2_xy_yyy_yzz_1,   \
                             ta2_xy_yyy_yzzz_0,  \
                             ta2_xy_yyy_yzzz_1,  \
                             ta2_xy_yyy_zzz_0,   \
                             ta2_xy_yyy_zzz_1,   \
                             ta2_xy_yyy_zzzz_0,  \
                             ta2_xy_yyy_zzzz_1,  \
                             ta2_xy_yyyz_xxxx_0, \
                             ta2_xy_yyyz_xxxy_0, \
                             ta2_xy_yyyz_xxxz_0, \
                             ta2_xy_yyyz_xxyy_0, \
                             ta2_xy_yyyz_xxyz_0, \
                             ta2_xy_yyyz_xxzz_0, \
                             ta2_xy_yyyz_xyyy_0, \
                             ta2_xy_yyyz_xyyz_0, \
                             ta2_xy_yyyz_xyzz_0, \
                             ta2_xy_yyyz_xzzz_0, \
                             ta2_xy_yyyz_yyyy_0, \
                             ta2_xy_yyyz_yyyz_0, \
                             ta2_xy_yyyz_yyzz_0, \
                             ta2_xy_yyyz_yzzz_0, \
                             ta2_xy_yyyz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyyz_xxxx_0[i] = ta2_xy_yyy_xxxx_0[i] * pa_z[i] - ta2_xy_yyy_xxxx_1[i] * pc_z[i];

        ta2_xy_yyyz_xxxy_0[i] = ta2_xy_yyy_xxxy_0[i] * pa_z[i] - ta2_xy_yyy_xxxy_1[i] * pc_z[i];

        ta2_xy_yyyz_xxxz_0[i] =
            ta2_xy_yyy_xxx_0[i] * fe_0 - ta2_xy_yyy_xxx_1[i] * fe_0 + ta2_xy_yyy_xxxz_0[i] * pa_z[i] - ta2_xy_yyy_xxxz_1[i] * pc_z[i];

        ta2_xy_yyyz_xxyy_0[i] = ta2_xy_yyy_xxyy_0[i] * pa_z[i] - ta2_xy_yyy_xxyy_1[i] * pc_z[i];

        ta2_xy_yyyz_xxyz_0[i] =
            ta2_xy_yyy_xxy_0[i] * fe_0 - ta2_xy_yyy_xxy_1[i] * fe_0 + ta2_xy_yyy_xxyz_0[i] * pa_z[i] - ta2_xy_yyy_xxyz_1[i] * pc_z[i];

        ta2_xy_yyyz_xxzz_0[i] =
            2.0 * ta2_xy_yyy_xxz_0[i] * fe_0 - 2.0 * ta2_xy_yyy_xxz_1[i] * fe_0 + ta2_xy_yyy_xxzz_0[i] * pa_z[i] - ta2_xy_yyy_xxzz_1[i] * pc_z[i];

        ta2_xy_yyyz_xyyy_0[i] = ta2_xy_yyy_xyyy_0[i] * pa_z[i] - ta2_xy_yyy_xyyy_1[i] * pc_z[i];

        ta2_xy_yyyz_xyyz_0[i] =
            ta2_xy_yyy_xyy_0[i] * fe_0 - ta2_xy_yyy_xyy_1[i] * fe_0 + ta2_xy_yyy_xyyz_0[i] * pa_z[i] - ta2_xy_yyy_xyyz_1[i] * pc_z[i];

        ta2_xy_yyyz_xyzz_0[i] =
            2.0 * ta2_xy_yyy_xyz_0[i] * fe_0 - 2.0 * ta2_xy_yyy_xyz_1[i] * fe_0 + ta2_xy_yyy_xyzz_0[i] * pa_z[i] - ta2_xy_yyy_xyzz_1[i] * pc_z[i];

        ta2_xy_yyyz_xzzz_0[i] =
            3.0 * ta2_xy_yyy_xzz_0[i] * fe_0 - 3.0 * ta2_xy_yyy_xzz_1[i] * fe_0 + ta2_xy_yyy_xzzz_0[i] * pa_z[i] - ta2_xy_yyy_xzzz_1[i] * pc_z[i];

        ta2_xy_yyyz_yyyy_0[i] = ta2_xy_yyy_yyyy_0[i] * pa_z[i] - ta2_xy_yyy_yyyy_1[i] * pc_z[i];

        ta2_xy_yyyz_yyyz_0[i] =
            ta2_xy_yyy_yyy_0[i] * fe_0 - ta2_xy_yyy_yyy_1[i] * fe_0 + ta2_xy_yyy_yyyz_0[i] * pa_z[i] - ta2_xy_yyy_yyyz_1[i] * pc_z[i];

        ta2_xy_yyyz_yyzz_0[i] =
            2.0 * ta2_xy_yyy_yyz_0[i] * fe_0 - 2.0 * ta2_xy_yyy_yyz_1[i] * fe_0 + ta2_xy_yyy_yyzz_0[i] * pa_z[i] - ta2_xy_yyy_yyzz_1[i] * pc_z[i];

        ta2_xy_yyyz_yzzz_0[i] =
            3.0 * ta2_xy_yyy_yzz_0[i] * fe_0 - 3.0 * ta2_xy_yyy_yzz_1[i] * fe_0 + ta2_xy_yyy_yzzz_0[i] * pa_z[i] - ta2_xy_yyy_yzzz_1[i] * pc_z[i];

        ta2_xy_yyyz_zzzz_0[i] =
            4.0 * ta2_xy_yyy_zzz_0[i] * fe_0 - 4.0 * ta2_xy_yyy_zzz_1[i] * fe_0 + ta2_xy_yyy_zzzz_0[i] * pa_z[i] - ta2_xy_yyy_zzzz_1[i] * pc_z[i];
    }

    // Set up 405-420 components of targeted buffer : GG

    auto ta2_xy_yyzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 405);

    auto ta2_xy_yyzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 406);

    auto ta2_xy_yyzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 407);

    auto ta2_xy_yyzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 408);

    auto ta2_xy_yyzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 409);

    auto ta2_xy_yyzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 410);

    auto ta2_xy_yyzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 411);

    auto ta2_xy_yyzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 412);

    auto ta2_xy_yyzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 413);

    auto ta2_xy_yyzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 414);

    auto ta2_xy_yyzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 415);

    auto ta2_xy_yyzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 416);

    auto ta2_xy_yyzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 417);

    auto ta2_xy_yyzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 418);

    auto ta2_xy_yyzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 419);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_yzz_xxxz_1,   \
                             ta1_x_yzz_xxzz_1,   \
                             ta1_x_yzz_xzzz_1,   \
                             ta1_x_yzz_zzzz_1,   \
                             ta2_xy_yy_xxxx_0,   \
                             ta2_xy_yy_xxxx_1,   \
                             ta2_xy_yy_xxxy_0,   \
                             ta2_xy_yy_xxxy_1,   \
                             ta2_xy_yy_xxyy_0,   \
                             ta2_xy_yy_xxyy_1,   \
                             ta2_xy_yy_xxyz_0,   \
                             ta2_xy_yy_xxyz_1,   \
                             ta2_xy_yy_xyyy_0,   \
                             ta2_xy_yy_xyyy_1,   \
                             ta2_xy_yy_xyyz_0,   \
                             ta2_xy_yy_xyyz_1,   \
                             ta2_xy_yy_xyzz_0,   \
                             ta2_xy_yy_xyzz_1,   \
                             ta2_xy_yy_yyyy_0,   \
                             ta2_xy_yy_yyyy_1,   \
                             ta2_xy_yy_yyyz_0,   \
                             ta2_xy_yy_yyyz_1,   \
                             ta2_xy_yy_yyzz_0,   \
                             ta2_xy_yy_yyzz_1,   \
                             ta2_xy_yy_yzzz_0,   \
                             ta2_xy_yy_yzzz_1,   \
                             ta2_xy_yyz_xxxx_0,  \
                             ta2_xy_yyz_xxxx_1,  \
                             ta2_xy_yyz_xxxy_0,  \
                             ta2_xy_yyz_xxxy_1,  \
                             ta2_xy_yyz_xxy_0,   \
                             ta2_xy_yyz_xxy_1,   \
                             ta2_xy_yyz_xxyy_0,  \
                             ta2_xy_yyz_xxyy_1,  \
                             ta2_xy_yyz_xxyz_0,  \
                             ta2_xy_yyz_xxyz_1,  \
                             ta2_xy_yyz_xyy_0,   \
                             ta2_xy_yyz_xyy_1,   \
                             ta2_xy_yyz_xyyy_0,  \
                             ta2_xy_yyz_xyyy_1,  \
                             ta2_xy_yyz_xyyz_0,  \
                             ta2_xy_yyz_xyyz_1,  \
                             ta2_xy_yyz_xyz_0,   \
                             ta2_xy_yyz_xyz_1,   \
                             ta2_xy_yyz_xyzz_0,  \
                             ta2_xy_yyz_xyzz_1,  \
                             ta2_xy_yyz_yyy_0,   \
                             ta2_xy_yyz_yyy_1,   \
                             ta2_xy_yyz_yyyy_0,  \
                             ta2_xy_yyz_yyyy_1,  \
                             ta2_xy_yyz_yyyz_0,  \
                             ta2_xy_yyz_yyyz_1,  \
                             ta2_xy_yyz_yyz_0,   \
                             ta2_xy_yyz_yyz_1,   \
                             ta2_xy_yyz_yyzz_0,  \
                             ta2_xy_yyz_yyzz_1,  \
                             ta2_xy_yyz_yzz_0,   \
                             ta2_xy_yyz_yzz_1,   \
                             ta2_xy_yyz_yzzz_0,  \
                             ta2_xy_yyz_yzzz_1,  \
                             ta2_xy_yyzz_xxxx_0, \
                             ta2_xy_yyzz_xxxy_0, \
                             ta2_xy_yyzz_xxxz_0, \
                             ta2_xy_yyzz_xxyy_0, \
                             ta2_xy_yyzz_xxyz_0, \
                             ta2_xy_yyzz_xxzz_0, \
                             ta2_xy_yyzz_xyyy_0, \
                             ta2_xy_yyzz_xyyz_0, \
                             ta2_xy_yyzz_xyzz_0, \
                             ta2_xy_yyzz_xzzz_0, \
                             ta2_xy_yyzz_yyyy_0, \
                             ta2_xy_yyzz_yyyz_0, \
                             ta2_xy_yyzz_yyzz_0, \
                             ta2_xy_yyzz_yzzz_0, \
                             ta2_xy_yyzz_zzzz_0, \
                             ta2_xy_yzz_xxxz_0,  \
                             ta2_xy_yzz_xxxz_1,  \
                             ta2_xy_yzz_xxzz_0,  \
                             ta2_xy_yzz_xxzz_1,  \
                             ta2_xy_yzz_xzzz_0,  \
                             ta2_xy_yzz_xzzz_1,  \
                             ta2_xy_yzz_zzzz_0,  \
                             ta2_xy_yzz_zzzz_1,  \
                             ta2_xy_zz_xxxz_0,   \
                             ta2_xy_zz_xxxz_1,   \
                             ta2_xy_zz_xxzz_0,   \
                             ta2_xy_zz_xxzz_1,   \
                             ta2_xy_zz_xzzz_0,   \
                             ta2_xy_zz_xzzz_1,   \
                             ta2_xy_zz_zzzz_0,   \
                             ta2_xy_zz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yyzz_xxxx_0[i] =
            ta2_xy_yy_xxxx_0[i] * fe_0 - ta2_xy_yy_xxxx_1[i] * fe_0 + ta2_xy_yyz_xxxx_0[i] * pa_z[i] - ta2_xy_yyz_xxxx_1[i] * pc_z[i];

        ta2_xy_yyzz_xxxy_0[i] =
            ta2_xy_yy_xxxy_0[i] * fe_0 - ta2_xy_yy_xxxy_1[i] * fe_0 + ta2_xy_yyz_xxxy_0[i] * pa_z[i] - ta2_xy_yyz_xxxy_1[i] * pc_z[i];

        ta2_xy_yyzz_xxxz_0[i] = ta2_xy_zz_xxxz_0[i] * fe_0 - ta2_xy_zz_xxxz_1[i] * fe_0 + ta1_x_yzz_xxxz_1[i] + ta2_xy_yzz_xxxz_0[i] * pa_y[i] -
                                ta2_xy_yzz_xxxz_1[i] * pc_y[i];

        ta2_xy_yyzz_xxyy_0[i] =
            ta2_xy_yy_xxyy_0[i] * fe_0 - ta2_xy_yy_xxyy_1[i] * fe_0 + ta2_xy_yyz_xxyy_0[i] * pa_z[i] - ta2_xy_yyz_xxyy_1[i] * pc_z[i];

        ta2_xy_yyzz_xxyz_0[i] = ta2_xy_yy_xxyz_0[i] * fe_0 - ta2_xy_yy_xxyz_1[i] * fe_0 + ta2_xy_yyz_xxy_0[i] * fe_0 - ta2_xy_yyz_xxy_1[i] * fe_0 +
                                ta2_xy_yyz_xxyz_0[i] * pa_z[i] - ta2_xy_yyz_xxyz_1[i] * pc_z[i];

        ta2_xy_yyzz_xxzz_0[i] = ta2_xy_zz_xxzz_0[i] * fe_0 - ta2_xy_zz_xxzz_1[i] * fe_0 + ta1_x_yzz_xxzz_1[i] + ta2_xy_yzz_xxzz_0[i] * pa_y[i] -
                                ta2_xy_yzz_xxzz_1[i] * pc_y[i];

        ta2_xy_yyzz_xyyy_0[i] =
            ta2_xy_yy_xyyy_0[i] * fe_0 - ta2_xy_yy_xyyy_1[i] * fe_0 + ta2_xy_yyz_xyyy_0[i] * pa_z[i] - ta2_xy_yyz_xyyy_1[i] * pc_z[i];

        ta2_xy_yyzz_xyyz_0[i] = ta2_xy_yy_xyyz_0[i] * fe_0 - ta2_xy_yy_xyyz_1[i] * fe_0 + ta2_xy_yyz_xyy_0[i] * fe_0 - ta2_xy_yyz_xyy_1[i] * fe_0 +
                                ta2_xy_yyz_xyyz_0[i] * pa_z[i] - ta2_xy_yyz_xyyz_1[i] * pc_z[i];

        ta2_xy_yyzz_xyzz_0[i] = ta2_xy_yy_xyzz_0[i] * fe_0 - ta2_xy_yy_xyzz_1[i] * fe_0 + 2.0 * ta2_xy_yyz_xyz_0[i] * fe_0 -
                                2.0 * ta2_xy_yyz_xyz_1[i] * fe_0 + ta2_xy_yyz_xyzz_0[i] * pa_z[i] - ta2_xy_yyz_xyzz_1[i] * pc_z[i];

        ta2_xy_yyzz_xzzz_0[i] = ta2_xy_zz_xzzz_0[i] * fe_0 - ta2_xy_zz_xzzz_1[i] * fe_0 + ta1_x_yzz_xzzz_1[i] + ta2_xy_yzz_xzzz_0[i] * pa_y[i] -
                                ta2_xy_yzz_xzzz_1[i] * pc_y[i];

        ta2_xy_yyzz_yyyy_0[i] =
            ta2_xy_yy_yyyy_0[i] * fe_0 - ta2_xy_yy_yyyy_1[i] * fe_0 + ta2_xy_yyz_yyyy_0[i] * pa_z[i] - ta2_xy_yyz_yyyy_1[i] * pc_z[i];

        ta2_xy_yyzz_yyyz_0[i] = ta2_xy_yy_yyyz_0[i] * fe_0 - ta2_xy_yy_yyyz_1[i] * fe_0 + ta2_xy_yyz_yyy_0[i] * fe_0 - ta2_xy_yyz_yyy_1[i] * fe_0 +
                                ta2_xy_yyz_yyyz_0[i] * pa_z[i] - ta2_xy_yyz_yyyz_1[i] * pc_z[i];

        ta2_xy_yyzz_yyzz_0[i] = ta2_xy_yy_yyzz_0[i] * fe_0 - ta2_xy_yy_yyzz_1[i] * fe_0 + 2.0 * ta2_xy_yyz_yyz_0[i] * fe_0 -
                                2.0 * ta2_xy_yyz_yyz_1[i] * fe_0 + ta2_xy_yyz_yyzz_0[i] * pa_z[i] - ta2_xy_yyz_yyzz_1[i] * pc_z[i];

        ta2_xy_yyzz_yzzz_0[i] = ta2_xy_yy_yzzz_0[i] * fe_0 - ta2_xy_yy_yzzz_1[i] * fe_0 + 3.0 * ta2_xy_yyz_yzz_0[i] * fe_0 -
                                3.0 * ta2_xy_yyz_yzz_1[i] * fe_0 + ta2_xy_yyz_yzzz_0[i] * pa_z[i] - ta2_xy_yyz_yzzz_1[i] * pc_z[i];

        ta2_xy_yyzz_zzzz_0[i] = ta2_xy_zz_zzzz_0[i] * fe_0 - ta2_xy_zz_zzzz_1[i] * fe_0 + ta1_x_yzz_zzzz_1[i] + ta2_xy_yzz_zzzz_0[i] * pa_y[i] -
                                ta2_xy_yzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 420-435 components of targeted buffer : GG

    auto ta2_xy_yzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 420);

    auto ta2_xy_yzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 421);

    auto ta2_xy_yzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 422);

    auto ta2_xy_yzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 423);

    auto ta2_xy_yzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 424);

    auto ta2_xy_yzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 425);

    auto ta2_xy_yzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 426);

    auto ta2_xy_yzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 427);

    auto ta2_xy_yzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 428);

    auto ta2_xy_yzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 429);

    auto ta2_xy_yzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 430);

    auto ta2_xy_yzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 431);

    auto ta2_xy_yzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 432);

    auto ta2_xy_yzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 433);

    auto ta2_xy_yzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 434);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_zzz_xxxx_1,   \
                             ta1_x_zzz_xxxz_1,   \
                             ta1_x_zzz_xxyz_1,   \
                             ta1_x_zzz_xxzz_1,   \
                             ta1_x_zzz_xyyz_1,   \
                             ta1_x_zzz_xyzz_1,   \
                             ta1_x_zzz_xzzz_1,   \
                             ta1_x_zzz_yyyz_1,   \
                             ta1_x_zzz_yyzz_1,   \
                             ta1_x_zzz_yzzz_1,   \
                             ta1_x_zzz_zzzz_1,   \
                             ta2_xy_yz_xxxy_0,   \
                             ta2_xy_yz_xxxy_1,   \
                             ta2_xy_yz_xxyy_0,   \
                             ta2_xy_yz_xxyy_1,   \
                             ta2_xy_yz_xyyy_0,   \
                             ta2_xy_yz_xyyy_1,   \
                             ta2_xy_yz_yyyy_0,   \
                             ta2_xy_yz_yyyy_1,   \
                             ta2_xy_yzz_xxxy_0,  \
                             ta2_xy_yzz_xxxy_1,  \
                             ta2_xy_yzz_xxyy_0,  \
                             ta2_xy_yzz_xxyy_1,  \
                             ta2_xy_yzz_xyyy_0,  \
                             ta2_xy_yzz_xyyy_1,  \
                             ta2_xy_yzz_yyyy_0,  \
                             ta2_xy_yzz_yyyy_1,  \
                             ta2_xy_yzzz_xxxx_0, \
                             ta2_xy_yzzz_xxxy_0, \
                             ta2_xy_yzzz_xxxz_0, \
                             ta2_xy_yzzz_xxyy_0, \
                             ta2_xy_yzzz_xxyz_0, \
                             ta2_xy_yzzz_xxzz_0, \
                             ta2_xy_yzzz_xyyy_0, \
                             ta2_xy_yzzz_xyyz_0, \
                             ta2_xy_yzzz_xyzz_0, \
                             ta2_xy_yzzz_xzzz_0, \
                             ta2_xy_yzzz_yyyy_0, \
                             ta2_xy_yzzz_yyyz_0, \
                             ta2_xy_yzzz_yyzz_0, \
                             ta2_xy_yzzz_yzzz_0, \
                             ta2_xy_yzzz_zzzz_0, \
                             ta2_xy_zzz_xxxx_0,  \
                             ta2_xy_zzz_xxxx_1,  \
                             ta2_xy_zzz_xxxz_0,  \
                             ta2_xy_zzz_xxxz_1,  \
                             ta2_xy_zzz_xxyz_0,  \
                             ta2_xy_zzz_xxyz_1,  \
                             ta2_xy_zzz_xxz_0,   \
                             ta2_xy_zzz_xxz_1,   \
                             ta2_xy_zzz_xxzz_0,  \
                             ta2_xy_zzz_xxzz_1,  \
                             ta2_xy_zzz_xyyz_0,  \
                             ta2_xy_zzz_xyyz_1,  \
                             ta2_xy_zzz_xyz_0,   \
                             ta2_xy_zzz_xyz_1,   \
                             ta2_xy_zzz_xyzz_0,  \
                             ta2_xy_zzz_xyzz_1,  \
                             ta2_xy_zzz_xzz_0,   \
                             ta2_xy_zzz_xzz_1,   \
                             ta2_xy_zzz_xzzz_0,  \
                             ta2_xy_zzz_xzzz_1,  \
                             ta2_xy_zzz_yyyz_0,  \
                             ta2_xy_zzz_yyyz_1,  \
                             ta2_xy_zzz_yyz_0,   \
                             ta2_xy_zzz_yyz_1,   \
                             ta2_xy_zzz_yyzz_0,  \
                             ta2_xy_zzz_yyzz_1,  \
                             ta2_xy_zzz_yzz_0,   \
                             ta2_xy_zzz_yzz_1,   \
                             ta2_xy_zzz_yzzz_0,  \
                             ta2_xy_zzz_yzzz_1,  \
                             ta2_xy_zzz_zzz_0,   \
                             ta2_xy_zzz_zzz_1,   \
                             ta2_xy_zzz_zzzz_0,  \
                             ta2_xy_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_yzzz_xxxx_0[i] = ta1_x_zzz_xxxx_1[i] + ta2_xy_zzz_xxxx_0[i] * pa_y[i] - ta2_xy_zzz_xxxx_1[i] * pc_y[i];

        ta2_xy_yzzz_xxxy_0[i] =
            2.0 * ta2_xy_yz_xxxy_0[i] * fe_0 - 2.0 * ta2_xy_yz_xxxy_1[i] * fe_0 + ta2_xy_yzz_xxxy_0[i] * pa_z[i] - ta2_xy_yzz_xxxy_1[i] * pc_z[i];

        ta2_xy_yzzz_xxxz_0[i] = ta1_x_zzz_xxxz_1[i] + ta2_xy_zzz_xxxz_0[i] * pa_y[i] - ta2_xy_zzz_xxxz_1[i] * pc_y[i];

        ta2_xy_yzzz_xxyy_0[i] =
            2.0 * ta2_xy_yz_xxyy_0[i] * fe_0 - 2.0 * ta2_xy_yz_xxyy_1[i] * fe_0 + ta2_xy_yzz_xxyy_0[i] * pa_z[i] - ta2_xy_yzz_xxyy_1[i] * pc_z[i];

        ta2_xy_yzzz_xxyz_0[i] = ta2_xy_zzz_xxz_0[i] * fe_0 - ta2_xy_zzz_xxz_1[i] * fe_0 + ta1_x_zzz_xxyz_1[i] + ta2_xy_zzz_xxyz_0[i] * pa_y[i] -
                                ta2_xy_zzz_xxyz_1[i] * pc_y[i];

        ta2_xy_yzzz_xxzz_0[i] = ta1_x_zzz_xxzz_1[i] + ta2_xy_zzz_xxzz_0[i] * pa_y[i] - ta2_xy_zzz_xxzz_1[i] * pc_y[i];

        ta2_xy_yzzz_xyyy_0[i] =
            2.0 * ta2_xy_yz_xyyy_0[i] * fe_0 - 2.0 * ta2_xy_yz_xyyy_1[i] * fe_0 + ta2_xy_yzz_xyyy_0[i] * pa_z[i] - ta2_xy_yzz_xyyy_1[i] * pc_z[i];

        ta2_xy_yzzz_xyyz_0[i] = 2.0 * ta2_xy_zzz_xyz_0[i] * fe_0 - 2.0 * ta2_xy_zzz_xyz_1[i] * fe_0 + ta1_x_zzz_xyyz_1[i] +
                                ta2_xy_zzz_xyyz_0[i] * pa_y[i] - ta2_xy_zzz_xyyz_1[i] * pc_y[i];

        ta2_xy_yzzz_xyzz_0[i] = ta2_xy_zzz_xzz_0[i] * fe_0 - ta2_xy_zzz_xzz_1[i] * fe_0 + ta1_x_zzz_xyzz_1[i] + ta2_xy_zzz_xyzz_0[i] * pa_y[i] -
                                ta2_xy_zzz_xyzz_1[i] * pc_y[i];

        ta2_xy_yzzz_xzzz_0[i] = ta1_x_zzz_xzzz_1[i] + ta2_xy_zzz_xzzz_0[i] * pa_y[i] - ta2_xy_zzz_xzzz_1[i] * pc_y[i];

        ta2_xy_yzzz_yyyy_0[i] =
            2.0 * ta2_xy_yz_yyyy_0[i] * fe_0 - 2.0 * ta2_xy_yz_yyyy_1[i] * fe_0 + ta2_xy_yzz_yyyy_0[i] * pa_z[i] - ta2_xy_yzz_yyyy_1[i] * pc_z[i];

        ta2_xy_yzzz_yyyz_0[i] = 3.0 * ta2_xy_zzz_yyz_0[i] * fe_0 - 3.0 * ta2_xy_zzz_yyz_1[i] * fe_0 + ta1_x_zzz_yyyz_1[i] +
                                ta2_xy_zzz_yyyz_0[i] * pa_y[i] - ta2_xy_zzz_yyyz_1[i] * pc_y[i];

        ta2_xy_yzzz_yyzz_0[i] = 2.0 * ta2_xy_zzz_yzz_0[i] * fe_0 - 2.0 * ta2_xy_zzz_yzz_1[i] * fe_0 + ta1_x_zzz_yyzz_1[i] +
                                ta2_xy_zzz_yyzz_0[i] * pa_y[i] - ta2_xy_zzz_yyzz_1[i] * pc_y[i];

        ta2_xy_yzzz_yzzz_0[i] = ta2_xy_zzz_zzz_0[i] * fe_0 - ta2_xy_zzz_zzz_1[i] * fe_0 + ta1_x_zzz_yzzz_1[i] + ta2_xy_zzz_yzzz_0[i] * pa_y[i] -
                                ta2_xy_zzz_yzzz_1[i] * pc_y[i];

        ta2_xy_yzzz_zzzz_0[i] = ta1_x_zzz_zzzz_1[i] + ta2_xy_zzz_zzzz_0[i] * pa_y[i] - ta2_xy_zzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 435-450 components of targeted buffer : GG

    auto ta2_xy_zzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 435);

    auto ta2_xy_zzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 436);

    auto ta2_xy_zzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 437);

    auto ta2_xy_zzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 438);

    auto ta2_xy_zzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 439);

    auto ta2_xy_zzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 440);

    auto ta2_xy_zzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 441);

    auto ta2_xy_zzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 442);

    auto ta2_xy_zzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 443);

    auto ta2_xy_zzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 444);

    auto ta2_xy_zzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 445);

    auto ta2_xy_zzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 446);

    auto ta2_xy_zzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 447);

    auto ta2_xy_zzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 448);

    auto ta2_xy_zzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 449);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta2_xy_zz_xxxx_0,   \
                             ta2_xy_zz_xxxx_1,   \
                             ta2_xy_zz_xxxy_0,   \
                             ta2_xy_zz_xxxy_1,   \
                             ta2_xy_zz_xxxz_0,   \
                             ta2_xy_zz_xxxz_1,   \
                             ta2_xy_zz_xxyy_0,   \
                             ta2_xy_zz_xxyy_1,   \
                             ta2_xy_zz_xxyz_0,   \
                             ta2_xy_zz_xxyz_1,   \
                             ta2_xy_zz_xxzz_0,   \
                             ta2_xy_zz_xxzz_1,   \
                             ta2_xy_zz_xyyy_0,   \
                             ta2_xy_zz_xyyy_1,   \
                             ta2_xy_zz_xyyz_0,   \
                             ta2_xy_zz_xyyz_1,   \
                             ta2_xy_zz_xyzz_0,   \
                             ta2_xy_zz_xyzz_1,   \
                             ta2_xy_zz_xzzz_0,   \
                             ta2_xy_zz_xzzz_1,   \
                             ta2_xy_zz_yyyy_0,   \
                             ta2_xy_zz_yyyy_1,   \
                             ta2_xy_zz_yyyz_0,   \
                             ta2_xy_zz_yyyz_1,   \
                             ta2_xy_zz_yyzz_0,   \
                             ta2_xy_zz_yyzz_1,   \
                             ta2_xy_zz_yzzz_0,   \
                             ta2_xy_zz_yzzz_1,   \
                             ta2_xy_zz_zzzz_0,   \
                             ta2_xy_zz_zzzz_1,   \
                             ta2_xy_zzz_xxx_0,   \
                             ta2_xy_zzz_xxx_1,   \
                             ta2_xy_zzz_xxxx_0,  \
                             ta2_xy_zzz_xxxx_1,  \
                             ta2_xy_zzz_xxxy_0,  \
                             ta2_xy_zzz_xxxy_1,  \
                             ta2_xy_zzz_xxxz_0,  \
                             ta2_xy_zzz_xxxz_1,  \
                             ta2_xy_zzz_xxy_0,   \
                             ta2_xy_zzz_xxy_1,   \
                             ta2_xy_zzz_xxyy_0,  \
                             ta2_xy_zzz_xxyy_1,  \
                             ta2_xy_zzz_xxyz_0,  \
                             ta2_xy_zzz_xxyz_1,  \
                             ta2_xy_zzz_xxz_0,   \
                             ta2_xy_zzz_xxz_1,   \
                             ta2_xy_zzz_xxzz_0,  \
                             ta2_xy_zzz_xxzz_1,  \
                             ta2_xy_zzz_xyy_0,   \
                             ta2_xy_zzz_xyy_1,   \
                             ta2_xy_zzz_xyyy_0,  \
                             ta2_xy_zzz_xyyy_1,  \
                             ta2_xy_zzz_xyyz_0,  \
                             ta2_xy_zzz_xyyz_1,  \
                             ta2_xy_zzz_xyz_0,   \
                             ta2_xy_zzz_xyz_1,   \
                             ta2_xy_zzz_xyzz_0,  \
                             ta2_xy_zzz_xyzz_1,  \
                             ta2_xy_zzz_xzz_0,   \
                             ta2_xy_zzz_xzz_1,   \
                             ta2_xy_zzz_xzzz_0,  \
                             ta2_xy_zzz_xzzz_1,  \
                             ta2_xy_zzz_yyy_0,   \
                             ta2_xy_zzz_yyy_1,   \
                             ta2_xy_zzz_yyyy_0,  \
                             ta2_xy_zzz_yyyy_1,  \
                             ta2_xy_zzz_yyyz_0,  \
                             ta2_xy_zzz_yyyz_1,  \
                             ta2_xy_zzz_yyz_0,   \
                             ta2_xy_zzz_yyz_1,   \
                             ta2_xy_zzz_yyzz_0,  \
                             ta2_xy_zzz_yyzz_1,  \
                             ta2_xy_zzz_yzz_0,   \
                             ta2_xy_zzz_yzz_1,   \
                             ta2_xy_zzz_yzzz_0,  \
                             ta2_xy_zzz_yzzz_1,  \
                             ta2_xy_zzz_zzz_0,   \
                             ta2_xy_zzz_zzz_1,   \
                             ta2_xy_zzz_zzzz_0,  \
                             ta2_xy_zzz_zzzz_1,  \
                             ta2_xy_zzzz_xxxx_0, \
                             ta2_xy_zzzz_xxxy_0, \
                             ta2_xy_zzzz_xxxz_0, \
                             ta2_xy_zzzz_xxyy_0, \
                             ta2_xy_zzzz_xxyz_0, \
                             ta2_xy_zzzz_xxzz_0, \
                             ta2_xy_zzzz_xyyy_0, \
                             ta2_xy_zzzz_xyyz_0, \
                             ta2_xy_zzzz_xyzz_0, \
                             ta2_xy_zzzz_xzzz_0, \
                             ta2_xy_zzzz_yyyy_0, \
                             ta2_xy_zzzz_yyyz_0, \
                             ta2_xy_zzzz_yyzz_0, \
                             ta2_xy_zzzz_yzzz_0, \
                             ta2_xy_zzzz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xy_zzzz_xxxx_0[i] =
            3.0 * ta2_xy_zz_xxxx_0[i] * fe_0 - 3.0 * ta2_xy_zz_xxxx_1[i] * fe_0 + ta2_xy_zzz_xxxx_0[i] * pa_z[i] - ta2_xy_zzz_xxxx_1[i] * pc_z[i];

        ta2_xy_zzzz_xxxy_0[i] =
            3.0 * ta2_xy_zz_xxxy_0[i] * fe_0 - 3.0 * ta2_xy_zz_xxxy_1[i] * fe_0 + ta2_xy_zzz_xxxy_0[i] * pa_z[i] - ta2_xy_zzz_xxxy_1[i] * pc_z[i];

        ta2_xy_zzzz_xxxz_0[i] = 3.0 * ta2_xy_zz_xxxz_0[i] * fe_0 - 3.0 * ta2_xy_zz_xxxz_1[i] * fe_0 + ta2_xy_zzz_xxx_0[i] * fe_0 -
                                ta2_xy_zzz_xxx_1[i] * fe_0 + ta2_xy_zzz_xxxz_0[i] * pa_z[i] - ta2_xy_zzz_xxxz_1[i] * pc_z[i];

        ta2_xy_zzzz_xxyy_0[i] =
            3.0 * ta2_xy_zz_xxyy_0[i] * fe_0 - 3.0 * ta2_xy_zz_xxyy_1[i] * fe_0 + ta2_xy_zzz_xxyy_0[i] * pa_z[i] - ta2_xy_zzz_xxyy_1[i] * pc_z[i];

        ta2_xy_zzzz_xxyz_0[i] = 3.0 * ta2_xy_zz_xxyz_0[i] * fe_0 - 3.0 * ta2_xy_zz_xxyz_1[i] * fe_0 + ta2_xy_zzz_xxy_0[i] * fe_0 -
                                ta2_xy_zzz_xxy_1[i] * fe_0 + ta2_xy_zzz_xxyz_0[i] * pa_z[i] - ta2_xy_zzz_xxyz_1[i] * pc_z[i];

        ta2_xy_zzzz_xxzz_0[i] = 3.0 * ta2_xy_zz_xxzz_0[i] * fe_0 - 3.0 * ta2_xy_zz_xxzz_1[i] * fe_0 + 2.0 * ta2_xy_zzz_xxz_0[i] * fe_0 -
                                2.0 * ta2_xy_zzz_xxz_1[i] * fe_0 + ta2_xy_zzz_xxzz_0[i] * pa_z[i] - ta2_xy_zzz_xxzz_1[i] * pc_z[i];

        ta2_xy_zzzz_xyyy_0[i] =
            3.0 * ta2_xy_zz_xyyy_0[i] * fe_0 - 3.0 * ta2_xy_zz_xyyy_1[i] * fe_0 + ta2_xy_zzz_xyyy_0[i] * pa_z[i] - ta2_xy_zzz_xyyy_1[i] * pc_z[i];

        ta2_xy_zzzz_xyyz_0[i] = 3.0 * ta2_xy_zz_xyyz_0[i] * fe_0 - 3.0 * ta2_xy_zz_xyyz_1[i] * fe_0 + ta2_xy_zzz_xyy_0[i] * fe_0 -
                                ta2_xy_zzz_xyy_1[i] * fe_0 + ta2_xy_zzz_xyyz_0[i] * pa_z[i] - ta2_xy_zzz_xyyz_1[i] * pc_z[i];

        ta2_xy_zzzz_xyzz_0[i] = 3.0 * ta2_xy_zz_xyzz_0[i] * fe_0 - 3.0 * ta2_xy_zz_xyzz_1[i] * fe_0 + 2.0 * ta2_xy_zzz_xyz_0[i] * fe_0 -
                                2.0 * ta2_xy_zzz_xyz_1[i] * fe_0 + ta2_xy_zzz_xyzz_0[i] * pa_z[i] - ta2_xy_zzz_xyzz_1[i] * pc_z[i];

        ta2_xy_zzzz_xzzz_0[i] = 3.0 * ta2_xy_zz_xzzz_0[i] * fe_0 - 3.0 * ta2_xy_zz_xzzz_1[i] * fe_0 + 3.0 * ta2_xy_zzz_xzz_0[i] * fe_0 -
                                3.0 * ta2_xy_zzz_xzz_1[i] * fe_0 + ta2_xy_zzz_xzzz_0[i] * pa_z[i] - ta2_xy_zzz_xzzz_1[i] * pc_z[i];

        ta2_xy_zzzz_yyyy_0[i] =
            3.0 * ta2_xy_zz_yyyy_0[i] * fe_0 - 3.0 * ta2_xy_zz_yyyy_1[i] * fe_0 + ta2_xy_zzz_yyyy_0[i] * pa_z[i] - ta2_xy_zzz_yyyy_1[i] * pc_z[i];

        ta2_xy_zzzz_yyyz_0[i] = 3.0 * ta2_xy_zz_yyyz_0[i] * fe_0 - 3.0 * ta2_xy_zz_yyyz_1[i] * fe_0 + ta2_xy_zzz_yyy_0[i] * fe_0 -
                                ta2_xy_zzz_yyy_1[i] * fe_0 + ta2_xy_zzz_yyyz_0[i] * pa_z[i] - ta2_xy_zzz_yyyz_1[i] * pc_z[i];

        ta2_xy_zzzz_yyzz_0[i] = 3.0 * ta2_xy_zz_yyzz_0[i] * fe_0 - 3.0 * ta2_xy_zz_yyzz_1[i] * fe_0 + 2.0 * ta2_xy_zzz_yyz_0[i] * fe_0 -
                                2.0 * ta2_xy_zzz_yyz_1[i] * fe_0 + ta2_xy_zzz_yyzz_0[i] * pa_z[i] - ta2_xy_zzz_yyzz_1[i] * pc_z[i];

        ta2_xy_zzzz_yzzz_0[i] = 3.0 * ta2_xy_zz_yzzz_0[i] * fe_0 - 3.0 * ta2_xy_zz_yzzz_1[i] * fe_0 + 3.0 * ta2_xy_zzz_yzz_0[i] * fe_0 -
                                3.0 * ta2_xy_zzz_yzz_1[i] * fe_0 + ta2_xy_zzz_yzzz_0[i] * pa_z[i] - ta2_xy_zzz_yzzz_1[i] * pc_z[i];

        ta2_xy_zzzz_zzzz_0[i] = 3.0 * ta2_xy_zz_zzzz_0[i] * fe_0 - 3.0 * ta2_xy_zz_zzzz_1[i] * fe_0 + 4.0 * ta2_xy_zzz_zzz_0[i] * fe_0 -
                                4.0 * ta2_xy_zzz_zzz_1[i] * fe_0 + ta2_xy_zzz_zzzz_0[i] * pa_z[i] - ta2_xy_zzz_zzzz_1[i] * pc_z[i];
    }

    // Set up 450-465 components of targeted buffer : GG

    auto ta2_xz_xxxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 450);

    auto ta2_xz_xxxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 451);

    auto ta2_xz_xxxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 452);

    auto ta2_xz_xxxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 453);

    auto ta2_xz_xxxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 454);

    auto ta2_xz_xxxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 455);

    auto ta2_xz_xxxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 456);

    auto ta2_xz_xxxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 457);

    auto ta2_xz_xxxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 458);

    auto ta2_xz_xxxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 459);

    auto ta2_xz_xxxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 460);

    auto ta2_xz_xxxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 461);

    auto ta2_xz_xxxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 462);

    auto ta2_xz_xxxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 463);

    auto ta2_xz_xxxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 464);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_xxx_xxxx_1,   \
                             ta1_z_xxx_xxxy_1,   \
                             ta1_z_xxx_xxxz_1,   \
                             ta1_z_xxx_xxyy_1,   \
                             ta1_z_xxx_xxyz_1,   \
                             ta1_z_xxx_xxzz_1,   \
                             ta1_z_xxx_xyyy_1,   \
                             ta1_z_xxx_xyyz_1,   \
                             ta1_z_xxx_xyzz_1,   \
                             ta1_z_xxx_xzzz_1,   \
                             ta1_z_xxx_yyyy_1,   \
                             ta1_z_xxx_yyyz_1,   \
                             ta1_z_xxx_yyzz_1,   \
                             ta1_z_xxx_yzzz_1,   \
                             ta1_z_xxx_zzzz_1,   \
                             ta2_xz_xx_xxxx_0,   \
                             ta2_xz_xx_xxxx_1,   \
                             ta2_xz_xx_xxxy_0,   \
                             ta2_xz_xx_xxxy_1,   \
                             ta2_xz_xx_xxxz_0,   \
                             ta2_xz_xx_xxxz_1,   \
                             ta2_xz_xx_xxyy_0,   \
                             ta2_xz_xx_xxyy_1,   \
                             ta2_xz_xx_xxyz_0,   \
                             ta2_xz_xx_xxyz_1,   \
                             ta2_xz_xx_xxzz_0,   \
                             ta2_xz_xx_xxzz_1,   \
                             ta2_xz_xx_xyyy_0,   \
                             ta2_xz_xx_xyyy_1,   \
                             ta2_xz_xx_xyyz_0,   \
                             ta2_xz_xx_xyyz_1,   \
                             ta2_xz_xx_xyzz_0,   \
                             ta2_xz_xx_xyzz_1,   \
                             ta2_xz_xx_xzzz_0,   \
                             ta2_xz_xx_xzzz_1,   \
                             ta2_xz_xx_yyyy_0,   \
                             ta2_xz_xx_yyyy_1,   \
                             ta2_xz_xx_yyyz_0,   \
                             ta2_xz_xx_yyyz_1,   \
                             ta2_xz_xx_yyzz_0,   \
                             ta2_xz_xx_yyzz_1,   \
                             ta2_xz_xx_yzzz_0,   \
                             ta2_xz_xx_yzzz_1,   \
                             ta2_xz_xx_zzzz_0,   \
                             ta2_xz_xx_zzzz_1,   \
                             ta2_xz_xxx_xxx_0,   \
                             ta2_xz_xxx_xxx_1,   \
                             ta2_xz_xxx_xxxx_0,  \
                             ta2_xz_xxx_xxxx_1,  \
                             ta2_xz_xxx_xxxy_0,  \
                             ta2_xz_xxx_xxxy_1,  \
                             ta2_xz_xxx_xxxz_0,  \
                             ta2_xz_xxx_xxxz_1,  \
                             ta2_xz_xxx_xxy_0,   \
                             ta2_xz_xxx_xxy_1,   \
                             ta2_xz_xxx_xxyy_0,  \
                             ta2_xz_xxx_xxyy_1,  \
                             ta2_xz_xxx_xxyz_0,  \
                             ta2_xz_xxx_xxyz_1,  \
                             ta2_xz_xxx_xxz_0,   \
                             ta2_xz_xxx_xxz_1,   \
                             ta2_xz_xxx_xxzz_0,  \
                             ta2_xz_xxx_xxzz_1,  \
                             ta2_xz_xxx_xyy_0,   \
                             ta2_xz_xxx_xyy_1,   \
                             ta2_xz_xxx_xyyy_0,  \
                             ta2_xz_xxx_xyyy_1,  \
                             ta2_xz_xxx_xyyz_0,  \
                             ta2_xz_xxx_xyyz_1,  \
                             ta2_xz_xxx_xyz_0,   \
                             ta2_xz_xxx_xyz_1,   \
                             ta2_xz_xxx_xyzz_0,  \
                             ta2_xz_xxx_xyzz_1,  \
                             ta2_xz_xxx_xzz_0,   \
                             ta2_xz_xxx_xzz_1,   \
                             ta2_xz_xxx_xzzz_0,  \
                             ta2_xz_xxx_xzzz_1,  \
                             ta2_xz_xxx_yyy_0,   \
                             ta2_xz_xxx_yyy_1,   \
                             ta2_xz_xxx_yyyy_0,  \
                             ta2_xz_xxx_yyyy_1,  \
                             ta2_xz_xxx_yyyz_0,  \
                             ta2_xz_xxx_yyyz_1,  \
                             ta2_xz_xxx_yyz_0,   \
                             ta2_xz_xxx_yyz_1,   \
                             ta2_xz_xxx_yyzz_0,  \
                             ta2_xz_xxx_yyzz_1,  \
                             ta2_xz_xxx_yzz_0,   \
                             ta2_xz_xxx_yzz_1,   \
                             ta2_xz_xxx_yzzz_0,  \
                             ta2_xz_xxx_yzzz_1,  \
                             ta2_xz_xxx_zzz_0,   \
                             ta2_xz_xxx_zzz_1,   \
                             ta2_xz_xxx_zzzz_0,  \
                             ta2_xz_xxx_zzzz_1,  \
                             ta2_xz_xxxx_xxxx_0, \
                             ta2_xz_xxxx_xxxy_0, \
                             ta2_xz_xxxx_xxxz_0, \
                             ta2_xz_xxxx_xxyy_0, \
                             ta2_xz_xxxx_xxyz_0, \
                             ta2_xz_xxxx_xxzz_0, \
                             ta2_xz_xxxx_xyyy_0, \
                             ta2_xz_xxxx_xyyz_0, \
                             ta2_xz_xxxx_xyzz_0, \
                             ta2_xz_xxxx_xzzz_0, \
                             ta2_xz_xxxx_yyyy_0, \
                             ta2_xz_xxxx_yyyz_0, \
                             ta2_xz_xxxx_yyzz_0, \
                             ta2_xz_xxxx_yzzz_0, \
                             ta2_xz_xxxx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxxx_xxxx_0[i] = 3.0 * ta2_xz_xx_xxxx_0[i] * fe_0 - 3.0 * ta2_xz_xx_xxxx_1[i] * fe_0 + 4.0 * ta2_xz_xxx_xxx_0[i] * fe_0 -
                                4.0 * ta2_xz_xxx_xxx_1[i] * fe_0 + ta1_z_xxx_xxxx_1[i] + ta2_xz_xxx_xxxx_0[i] * pa_x[i] -
                                ta2_xz_xxx_xxxx_1[i] * pc_x[i];

        ta2_xz_xxxx_xxxy_0[i] = 3.0 * ta2_xz_xx_xxxy_0[i] * fe_0 - 3.0 * ta2_xz_xx_xxxy_1[i] * fe_0 + 3.0 * ta2_xz_xxx_xxy_0[i] * fe_0 -
                                3.0 * ta2_xz_xxx_xxy_1[i] * fe_0 + ta1_z_xxx_xxxy_1[i] + ta2_xz_xxx_xxxy_0[i] * pa_x[i] -
                                ta2_xz_xxx_xxxy_1[i] * pc_x[i];

        ta2_xz_xxxx_xxxz_0[i] = 3.0 * ta2_xz_xx_xxxz_0[i] * fe_0 - 3.0 * ta2_xz_xx_xxxz_1[i] * fe_0 + 3.0 * ta2_xz_xxx_xxz_0[i] * fe_0 -
                                3.0 * ta2_xz_xxx_xxz_1[i] * fe_0 + ta1_z_xxx_xxxz_1[i] + ta2_xz_xxx_xxxz_0[i] * pa_x[i] -
                                ta2_xz_xxx_xxxz_1[i] * pc_x[i];

        ta2_xz_xxxx_xxyy_0[i] = 3.0 * ta2_xz_xx_xxyy_0[i] * fe_0 - 3.0 * ta2_xz_xx_xxyy_1[i] * fe_0 + 2.0 * ta2_xz_xxx_xyy_0[i] * fe_0 -
                                2.0 * ta2_xz_xxx_xyy_1[i] * fe_0 + ta1_z_xxx_xxyy_1[i] + ta2_xz_xxx_xxyy_0[i] * pa_x[i] -
                                ta2_xz_xxx_xxyy_1[i] * pc_x[i];

        ta2_xz_xxxx_xxyz_0[i] = 3.0 * ta2_xz_xx_xxyz_0[i] * fe_0 - 3.0 * ta2_xz_xx_xxyz_1[i] * fe_0 + 2.0 * ta2_xz_xxx_xyz_0[i] * fe_0 -
                                2.0 * ta2_xz_xxx_xyz_1[i] * fe_0 + ta1_z_xxx_xxyz_1[i] + ta2_xz_xxx_xxyz_0[i] * pa_x[i] -
                                ta2_xz_xxx_xxyz_1[i] * pc_x[i];

        ta2_xz_xxxx_xxzz_0[i] = 3.0 * ta2_xz_xx_xxzz_0[i] * fe_0 - 3.0 * ta2_xz_xx_xxzz_1[i] * fe_0 + 2.0 * ta2_xz_xxx_xzz_0[i] * fe_0 -
                                2.0 * ta2_xz_xxx_xzz_1[i] * fe_0 + ta1_z_xxx_xxzz_1[i] + ta2_xz_xxx_xxzz_0[i] * pa_x[i] -
                                ta2_xz_xxx_xxzz_1[i] * pc_x[i];

        ta2_xz_xxxx_xyyy_0[i] = 3.0 * ta2_xz_xx_xyyy_0[i] * fe_0 - 3.0 * ta2_xz_xx_xyyy_1[i] * fe_0 + ta2_xz_xxx_yyy_0[i] * fe_0 -
                                ta2_xz_xxx_yyy_1[i] * fe_0 + ta1_z_xxx_xyyy_1[i] + ta2_xz_xxx_xyyy_0[i] * pa_x[i] - ta2_xz_xxx_xyyy_1[i] * pc_x[i];

        ta2_xz_xxxx_xyyz_0[i] = 3.0 * ta2_xz_xx_xyyz_0[i] * fe_0 - 3.0 * ta2_xz_xx_xyyz_1[i] * fe_0 + ta2_xz_xxx_yyz_0[i] * fe_0 -
                                ta2_xz_xxx_yyz_1[i] * fe_0 + ta1_z_xxx_xyyz_1[i] + ta2_xz_xxx_xyyz_0[i] * pa_x[i] - ta2_xz_xxx_xyyz_1[i] * pc_x[i];

        ta2_xz_xxxx_xyzz_0[i] = 3.0 * ta2_xz_xx_xyzz_0[i] * fe_0 - 3.0 * ta2_xz_xx_xyzz_1[i] * fe_0 + ta2_xz_xxx_yzz_0[i] * fe_0 -
                                ta2_xz_xxx_yzz_1[i] * fe_0 + ta1_z_xxx_xyzz_1[i] + ta2_xz_xxx_xyzz_0[i] * pa_x[i] - ta2_xz_xxx_xyzz_1[i] * pc_x[i];

        ta2_xz_xxxx_xzzz_0[i] = 3.0 * ta2_xz_xx_xzzz_0[i] * fe_0 - 3.0 * ta2_xz_xx_xzzz_1[i] * fe_0 + ta2_xz_xxx_zzz_0[i] * fe_0 -
                                ta2_xz_xxx_zzz_1[i] * fe_0 + ta1_z_xxx_xzzz_1[i] + ta2_xz_xxx_xzzz_0[i] * pa_x[i] - ta2_xz_xxx_xzzz_1[i] * pc_x[i];

        ta2_xz_xxxx_yyyy_0[i] = 3.0 * ta2_xz_xx_yyyy_0[i] * fe_0 - 3.0 * ta2_xz_xx_yyyy_1[i] * fe_0 + ta1_z_xxx_yyyy_1[i] +
                                ta2_xz_xxx_yyyy_0[i] * pa_x[i] - ta2_xz_xxx_yyyy_1[i] * pc_x[i];

        ta2_xz_xxxx_yyyz_0[i] = 3.0 * ta2_xz_xx_yyyz_0[i] * fe_0 - 3.0 * ta2_xz_xx_yyyz_1[i] * fe_0 + ta1_z_xxx_yyyz_1[i] +
                                ta2_xz_xxx_yyyz_0[i] * pa_x[i] - ta2_xz_xxx_yyyz_1[i] * pc_x[i];

        ta2_xz_xxxx_yyzz_0[i] = 3.0 * ta2_xz_xx_yyzz_0[i] * fe_0 - 3.0 * ta2_xz_xx_yyzz_1[i] * fe_0 + ta1_z_xxx_yyzz_1[i] +
                                ta2_xz_xxx_yyzz_0[i] * pa_x[i] - ta2_xz_xxx_yyzz_1[i] * pc_x[i];

        ta2_xz_xxxx_yzzz_0[i] = 3.0 * ta2_xz_xx_yzzz_0[i] * fe_0 - 3.0 * ta2_xz_xx_yzzz_1[i] * fe_0 + ta1_z_xxx_yzzz_1[i] +
                                ta2_xz_xxx_yzzz_0[i] * pa_x[i] - ta2_xz_xxx_yzzz_1[i] * pc_x[i];

        ta2_xz_xxxx_zzzz_0[i] = 3.0 * ta2_xz_xx_zzzz_0[i] * fe_0 - 3.0 * ta2_xz_xx_zzzz_1[i] * fe_0 + ta1_z_xxx_zzzz_1[i] +
                                ta2_xz_xxx_zzzz_0[i] * pa_x[i] - ta2_xz_xxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 465-480 components of targeted buffer : GG

    auto ta2_xz_xxxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 465);

    auto ta2_xz_xxxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 466);

    auto ta2_xz_xxxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 467);

    auto ta2_xz_xxxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 468);

    auto ta2_xz_xxxy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 469);

    auto ta2_xz_xxxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 470);

    auto ta2_xz_xxxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 471);

    auto ta2_xz_xxxy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 472);

    auto ta2_xz_xxxy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 473);

    auto ta2_xz_xxxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 474);

    auto ta2_xz_xxxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 475);

    auto ta2_xz_xxxy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 476);

    auto ta2_xz_xxxy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 477);

    auto ta2_xz_xxxy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 478);

    auto ta2_xz_xxxy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 479);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta2_xz_xxx_xxx_0,   \
                             ta2_xz_xxx_xxx_1,   \
                             ta2_xz_xxx_xxxx_0,  \
                             ta2_xz_xxx_xxxx_1,  \
                             ta2_xz_xxx_xxxy_0,  \
                             ta2_xz_xxx_xxxy_1,  \
                             ta2_xz_xxx_xxxz_0,  \
                             ta2_xz_xxx_xxxz_1,  \
                             ta2_xz_xxx_xxy_0,   \
                             ta2_xz_xxx_xxy_1,   \
                             ta2_xz_xxx_xxyy_0,  \
                             ta2_xz_xxx_xxyy_1,  \
                             ta2_xz_xxx_xxyz_0,  \
                             ta2_xz_xxx_xxyz_1,  \
                             ta2_xz_xxx_xxz_0,   \
                             ta2_xz_xxx_xxz_1,   \
                             ta2_xz_xxx_xxzz_0,  \
                             ta2_xz_xxx_xxzz_1,  \
                             ta2_xz_xxx_xyy_0,   \
                             ta2_xz_xxx_xyy_1,   \
                             ta2_xz_xxx_xyyy_0,  \
                             ta2_xz_xxx_xyyy_1,  \
                             ta2_xz_xxx_xyyz_0,  \
                             ta2_xz_xxx_xyyz_1,  \
                             ta2_xz_xxx_xyz_0,   \
                             ta2_xz_xxx_xyz_1,   \
                             ta2_xz_xxx_xyzz_0,  \
                             ta2_xz_xxx_xyzz_1,  \
                             ta2_xz_xxx_xzz_0,   \
                             ta2_xz_xxx_xzz_1,   \
                             ta2_xz_xxx_xzzz_0,  \
                             ta2_xz_xxx_xzzz_1,  \
                             ta2_xz_xxx_yyy_0,   \
                             ta2_xz_xxx_yyy_1,   \
                             ta2_xz_xxx_yyyy_0,  \
                             ta2_xz_xxx_yyyy_1,  \
                             ta2_xz_xxx_yyyz_0,  \
                             ta2_xz_xxx_yyyz_1,  \
                             ta2_xz_xxx_yyz_0,   \
                             ta2_xz_xxx_yyz_1,   \
                             ta2_xz_xxx_yyzz_0,  \
                             ta2_xz_xxx_yyzz_1,  \
                             ta2_xz_xxx_yzz_0,   \
                             ta2_xz_xxx_yzz_1,   \
                             ta2_xz_xxx_yzzz_0,  \
                             ta2_xz_xxx_yzzz_1,  \
                             ta2_xz_xxx_zzz_0,   \
                             ta2_xz_xxx_zzz_1,   \
                             ta2_xz_xxx_zzzz_0,  \
                             ta2_xz_xxx_zzzz_1,  \
                             ta2_xz_xxxy_xxxx_0, \
                             ta2_xz_xxxy_xxxy_0, \
                             ta2_xz_xxxy_xxxz_0, \
                             ta2_xz_xxxy_xxyy_0, \
                             ta2_xz_xxxy_xxyz_0, \
                             ta2_xz_xxxy_xxzz_0, \
                             ta2_xz_xxxy_xyyy_0, \
                             ta2_xz_xxxy_xyyz_0, \
                             ta2_xz_xxxy_xyzz_0, \
                             ta2_xz_xxxy_xzzz_0, \
                             ta2_xz_xxxy_yyyy_0, \
                             ta2_xz_xxxy_yyyz_0, \
                             ta2_xz_xxxy_yyzz_0, \
                             ta2_xz_xxxy_yzzz_0, \
                             ta2_xz_xxxy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxxy_xxxx_0[i] = ta2_xz_xxx_xxxx_0[i] * pa_y[i] - ta2_xz_xxx_xxxx_1[i] * pc_y[i];

        ta2_xz_xxxy_xxxy_0[i] =
            ta2_xz_xxx_xxx_0[i] * fe_0 - ta2_xz_xxx_xxx_1[i] * fe_0 + ta2_xz_xxx_xxxy_0[i] * pa_y[i] - ta2_xz_xxx_xxxy_1[i] * pc_y[i];

        ta2_xz_xxxy_xxxz_0[i] = ta2_xz_xxx_xxxz_0[i] * pa_y[i] - ta2_xz_xxx_xxxz_1[i] * pc_y[i];

        ta2_xz_xxxy_xxyy_0[i] =
            2.0 * ta2_xz_xxx_xxy_0[i] * fe_0 - 2.0 * ta2_xz_xxx_xxy_1[i] * fe_0 + ta2_xz_xxx_xxyy_0[i] * pa_y[i] - ta2_xz_xxx_xxyy_1[i] * pc_y[i];

        ta2_xz_xxxy_xxyz_0[i] =
            ta2_xz_xxx_xxz_0[i] * fe_0 - ta2_xz_xxx_xxz_1[i] * fe_0 + ta2_xz_xxx_xxyz_0[i] * pa_y[i] - ta2_xz_xxx_xxyz_1[i] * pc_y[i];

        ta2_xz_xxxy_xxzz_0[i] = ta2_xz_xxx_xxzz_0[i] * pa_y[i] - ta2_xz_xxx_xxzz_1[i] * pc_y[i];

        ta2_xz_xxxy_xyyy_0[i] =
            3.0 * ta2_xz_xxx_xyy_0[i] * fe_0 - 3.0 * ta2_xz_xxx_xyy_1[i] * fe_0 + ta2_xz_xxx_xyyy_0[i] * pa_y[i] - ta2_xz_xxx_xyyy_1[i] * pc_y[i];

        ta2_xz_xxxy_xyyz_0[i] =
            2.0 * ta2_xz_xxx_xyz_0[i] * fe_0 - 2.0 * ta2_xz_xxx_xyz_1[i] * fe_0 + ta2_xz_xxx_xyyz_0[i] * pa_y[i] - ta2_xz_xxx_xyyz_1[i] * pc_y[i];

        ta2_xz_xxxy_xyzz_0[i] =
            ta2_xz_xxx_xzz_0[i] * fe_0 - ta2_xz_xxx_xzz_1[i] * fe_0 + ta2_xz_xxx_xyzz_0[i] * pa_y[i] - ta2_xz_xxx_xyzz_1[i] * pc_y[i];

        ta2_xz_xxxy_xzzz_0[i] = ta2_xz_xxx_xzzz_0[i] * pa_y[i] - ta2_xz_xxx_xzzz_1[i] * pc_y[i];

        ta2_xz_xxxy_yyyy_0[i] =
            4.0 * ta2_xz_xxx_yyy_0[i] * fe_0 - 4.0 * ta2_xz_xxx_yyy_1[i] * fe_0 + ta2_xz_xxx_yyyy_0[i] * pa_y[i] - ta2_xz_xxx_yyyy_1[i] * pc_y[i];

        ta2_xz_xxxy_yyyz_0[i] =
            3.0 * ta2_xz_xxx_yyz_0[i] * fe_0 - 3.0 * ta2_xz_xxx_yyz_1[i] * fe_0 + ta2_xz_xxx_yyyz_0[i] * pa_y[i] - ta2_xz_xxx_yyyz_1[i] * pc_y[i];

        ta2_xz_xxxy_yyzz_0[i] =
            2.0 * ta2_xz_xxx_yzz_0[i] * fe_0 - 2.0 * ta2_xz_xxx_yzz_1[i] * fe_0 + ta2_xz_xxx_yyzz_0[i] * pa_y[i] - ta2_xz_xxx_yyzz_1[i] * pc_y[i];

        ta2_xz_xxxy_yzzz_0[i] =
            ta2_xz_xxx_zzz_0[i] * fe_0 - ta2_xz_xxx_zzz_1[i] * fe_0 + ta2_xz_xxx_yzzz_0[i] * pa_y[i] - ta2_xz_xxx_yzzz_1[i] * pc_y[i];

        ta2_xz_xxxy_zzzz_0[i] = ta2_xz_xxx_zzzz_0[i] * pa_y[i] - ta2_xz_xxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 480-495 components of targeted buffer : GG

    auto ta2_xz_xxxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 480);

    auto ta2_xz_xxxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 481);

    auto ta2_xz_xxxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 482);

    auto ta2_xz_xxxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 483);

    auto ta2_xz_xxxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 484);

    auto ta2_xz_xxxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 485);

    auto ta2_xz_xxxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 486);

    auto ta2_xz_xxxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 487);

    auto ta2_xz_xxxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 488);

    auto ta2_xz_xxxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 489);

    auto ta2_xz_xxxz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 490);

    auto ta2_xz_xxxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 491);

    auto ta2_xz_xxxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 492);

    auto ta2_xz_xxxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 493);

    auto ta2_xz_xxxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 494);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_x_xxx_xxxx_1,   \
                             ta1_x_xxx_xxxy_1,   \
                             ta1_x_xxx_xxxz_1,   \
                             ta1_x_xxx_xxyy_1,   \
                             ta1_x_xxx_xxyz_1,   \
                             ta1_x_xxx_xxzz_1,   \
                             ta1_x_xxx_xyyy_1,   \
                             ta1_x_xxx_xyyz_1,   \
                             ta1_x_xxx_xyzz_1,   \
                             ta1_x_xxx_xzzz_1,   \
                             ta1_x_xxx_yyyy_1,   \
                             ta1_z_xxz_yyyz_1,   \
                             ta1_z_xxz_yyzz_1,   \
                             ta1_z_xxz_yzzz_1,   \
                             ta1_z_xxz_zzzz_1,   \
                             ta2_xz_xxx_xxx_0,   \
                             ta2_xz_xxx_xxx_1,   \
                             ta2_xz_xxx_xxxx_0,  \
                             ta2_xz_xxx_xxxx_1,  \
                             ta2_xz_xxx_xxxy_0,  \
                             ta2_xz_xxx_xxxy_1,  \
                             ta2_xz_xxx_xxxz_0,  \
                             ta2_xz_xxx_xxxz_1,  \
                             ta2_xz_xxx_xxy_0,   \
                             ta2_xz_xxx_xxy_1,   \
                             ta2_xz_xxx_xxyy_0,  \
                             ta2_xz_xxx_xxyy_1,  \
                             ta2_xz_xxx_xxyz_0,  \
                             ta2_xz_xxx_xxyz_1,  \
                             ta2_xz_xxx_xxz_0,   \
                             ta2_xz_xxx_xxz_1,   \
                             ta2_xz_xxx_xxzz_0,  \
                             ta2_xz_xxx_xxzz_1,  \
                             ta2_xz_xxx_xyy_0,   \
                             ta2_xz_xxx_xyy_1,   \
                             ta2_xz_xxx_xyyy_0,  \
                             ta2_xz_xxx_xyyy_1,  \
                             ta2_xz_xxx_xyyz_0,  \
                             ta2_xz_xxx_xyyz_1,  \
                             ta2_xz_xxx_xyz_0,   \
                             ta2_xz_xxx_xyz_1,   \
                             ta2_xz_xxx_xyzz_0,  \
                             ta2_xz_xxx_xyzz_1,  \
                             ta2_xz_xxx_xzz_0,   \
                             ta2_xz_xxx_xzz_1,   \
                             ta2_xz_xxx_xzzz_0,  \
                             ta2_xz_xxx_xzzz_1,  \
                             ta2_xz_xxx_yyyy_0,  \
                             ta2_xz_xxx_yyyy_1,  \
                             ta2_xz_xxxz_xxxx_0, \
                             ta2_xz_xxxz_xxxy_0, \
                             ta2_xz_xxxz_xxxz_0, \
                             ta2_xz_xxxz_xxyy_0, \
                             ta2_xz_xxxz_xxyz_0, \
                             ta2_xz_xxxz_xxzz_0, \
                             ta2_xz_xxxz_xyyy_0, \
                             ta2_xz_xxxz_xyyz_0, \
                             ta2_xz_xxxz_xyzz_0, \
                             ta2_xz_xxxz_xzzz_0, \
                             ta2_xz_xxxz_yyyy_0, \
                             ta2_xz_xxxz_yyyz_0, \
                             ta2_xz_xxxz_yyzz_0, \
                             ta2_xz_xxxz_yzzz_0, \
                             ta2_xz_xxxz_zzzz_0, \
                             ta2_xz_xxz_yyyz_0,  \
                             ta2_xz_xxz_yyyz_1,  \
                             ta2_xz_xxz_yyzz_0,  \
                             ta2_xz_xxz_yyzz_1,  \
                             ta2_xz_xxz_yzzz_0,  \
                             ta2_xz_xxz_yzzz_1,  \
                             ta2_xz_xxz_zzzz_0,  \
                             ta2_xz_xxz_zzzz_1,  \
                             ta2_xz_xz_yyyz_0,   \
                             ta2_xz_xz_yyyz_1,   \
                             ta2_xz_xz_yyzz_0,   \
                             ta2_xz_xz_yyzz_1,   \
                             ta2_xz_xz_yzzz_0,   \
                             ta2_xz_xz_yzzz_1,   \
                             ta2_xz_xz_zzzz_0,   \
                             ta2_xz_xz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxxz_xxxx_0[i] = ta1_x_xxx_xxxx_1[i] + ta2_xz_xxx_xxxx_0[i] * pa_z[i] - ta2_xz_xxx_xxxx_1[i] * pc_z[i];

        ta2_xz_xxxz_xxxy_0[i] = ta1_x_xxx_xxxy_1[i] + ta2_xz_xxx_xxxy_0[i] * pa_z[i] - ta2_xz_xxx_xxxy_1[i] * pc_z[i];

        ta2_xz_xxxz_xxxz_0[i] = ta2_xz_xxx_xxx_0[i] * fe_0 - ta2_xz_xxx_xxx_1[i] * fe_0 + ta1_x_xxx_xxxz_1[i] + ta2_xz_xxx_xxxz_0[i] * pa_z[i] -
                                ta2_xz_xxx_xxxz_1[i] * pc_z[i];

        ta2_xz_xxxz_xxyy_0[i] = ta1_x_xxx_xxyy_1[i] + ta2_xz_xxx_xxyy_0[i] * pa_z[i] - ta2_xz_xxx_xxyy_1[i] * pc_z[i];

        ta2_xz_xxxz_xxyz_0[i] = ta2_xz_xxx_xxy_0[i] * fe_0 - ta2_xz_xxx_xxy_1[i] * fe_0 + ta1_x_xxx_xxyz_1[i] + ta2_xz_xxx_xxyz_0[i] * pa_z[i] -
                                ta2_xz_xxx_xxyz_1[i] * pc_z[i];

        ta2_xz_xxxz_xxzz_0[i] = 2.0 * ta2_xz_xxx_xxz_0[i] * fe_0 - 2.0 * ta2_xz_xxx_xxz_1[i] * fe_0 + ta1_x_xxx_xxzz_1[i] +
                                ta2_xz_xxx_xxzz_0[i] * pa_z[i] - ta2_xz_xxx_xxzz_1[i] * pc_z[i];

        ta2_xz_xxxz_xyyy_0[i] = ta1_x_xxx_xyyy_1[i] + ta2_xz_xxx_xyyy_0[i] * pa_z[i] - ta2_xz_xxx_xyyy_1[i] * pc_z[i];

        ta2_xz_xxxz_xyyz_0[i] = ta2_xz_xxx_xyy_0[i] * fe_0 - ta2_xz_xxx_xyy_1[i] * fe_0 + ta1_x_xxx_xyyz_1[i] + ta2_xz_xxx_xyyz_0[i] * pa_z[i] -
                                ta2_xz_xxx_xyyz_1[i] * pc_z[i];

        ta2_xz_xxxz_xyzz_0[i] = 2.0 * ta2_xz_xxx_xyz_0[i] * fe_0 - 2.0 * ta2_xz_xxx_xyz_1[i] * fe_0 + ta1_x_xxx_xyzz_1[i] +
                                ta2_xz_xxx_xyzz_0[i] * pa_z[i] - ta2_xz_xxx_xyzz_1[i] * pc_z[i];

        ta2_xz_xxxz_xzzz_0[i] = 3.0 * ta2_xz_xxx_xzz_0[i] * fe_0 - 3.0 * ta2_xz_xxx_xzz_1[i] * fe_0 + ta1_x_xxx_xzzz_1[i] +
                                ta2_xz_xxx_xzzz_0[i] * pa_z[i] - ta2_xz_xxx_xzzz_1[i] * pc_z[i];

        ta2_xz_xxxz_yyyy_0[i] = ta1_x_xxx_yyyy_1[i] + ta2_xz_xxx_yyyy_0[i] * pa_z[i] - ta2_xz_xxx_yyyy_1[i] * pc_z[i];

        ta2_xz_xxxz_yyyz_0[i] = 2.0 * ta2_xz_xz_yyyz_0[i] * fe_0 - 2.0 * ta2_xz_xz_yyyz_1[i] * fe_0 + ta1_z_xxz_yyyz_1[i] +
                                ta2_xz_xxz_yyyz_0[i] * pa_x[i] - ta2_xz_xxz_yyyz_1[i] * pc_x[i];

        ta2_xz_xxxz_yyzz_0[i] = 2.0 * ta2_xz_xz_yyzz_0[i] * fe_0 - 2.0 * ta2_xz_xz_yyzz_1[i] * fe_0 + ta1_z_xxz_yyzz_1[i] +
                                ta2_xz_xxz_yyzz_0[i] * pa_x[i] - ta2_xz_xxz_yyzz_1[i] * pc_x[i];

        ta2_xz_xxxz_yzzz_0[i] = 2.0 * ta2_xz_xz_yzzz_0[i] * fe_0 - 2.0 * ta2_xz_xz_yzzz_1[i] * fe_0 + ta1_z_xxz_yzzz_1[i] +
                                ta2_xz_xxz_yzzz_0[i] * pa_x[i] - ta2_xz_xxz_yzzz_1[i] * pc_x[i];

        ta2_xz_xxxz_zzzz_0[i] = 2.0 * ta2_xz_xz_zzzz_0[i] * fe_0 - 2.0 * ta2_xz_xz_zzzz_1[i] * fe_0 + ta1_z_xxz_zzzz_1[i] +
                                ta2_xz_xxz_zzzz_0[i] * pa_x[i] - ta2_xz_xxz_zzzz_1[i] * pc_x[i];
    }

    // Set up 495-510 components of targeted buffer : GG

    auto ta2_xz_xxyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 495);

    auto ta2_xz_xxyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 496);

    auto ta2_xz_xxyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 497);

    auto ta2_xz_xxyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 498);

    auto ta2_xz_xxyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 499);

    auto ta2_xz_xxyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 500);

    auto ta2_xz_xxyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 501);

    auto ta2_xz_xxyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 502);

    auto ta2_xz_xxyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 503);

    auto ta2_xz_xxyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 504);

    auto ta2_xz_xxyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 505);

    auto ta2_xz_xxyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 506);

    auto ta2_xz_xxyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 507);

    auto ta2_xz_xxyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 508);

    auto ta2_xz_xxyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 509);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xyy_yyyy_1,   \
                             ta1_z_xyy_yyyz_1,   \
                             ta1_z_xyy_yyzz_1,   \
                             ta1_z_xyy_yzzz_1,   \
                             ta2_xz_xx_xxxx_0,   \
                             ta2_xz_xx_xxxx_1,   \
                             ta2_xz_xx_xxxy_0,   \
                             ta2_xz_xx_xxxy_1,   \
                             ta2_xz_xx_xxxz_0,   \
                             ta2_xz_xx_xxxz_1,   \
                             ta2_xz_xx_xxyy_0,   \
                             ta2_xz_xx_xxyy_1,   \
                             ta2_xz_xx_xxyz_0,   \
                             ta2_xz_xx_xxyz_1,   \
                             ta2_xz_xx_xxzz_0,   \
                             ta2_xz_xx_xxzz_1,   \
                             ta2_xz_xx_xyyy_0,   \
                             ta2_xz_xx_xyyy_1,   \
                             ta2_xz_xx_xyyz_0,   \
                             ta2_xz_xx_xyyz_1,   \
                             ta2_xz_xx_xyzz_0,   \
                             ta2_xz_xx_xyzz_1,   \
                             ta2_xz_xx_xzzz_0,   \
                             ta2_xz_xx_xzzz_1,   \
                             ta2_xz_xx_zzzz_0,   \
                             ta2_xz_xx_zzzz_1,   \
                             ta2_xz_xxy_xxx_0,   \
                             ta2_xz_xxy_xxx_1,   \
                             ta2_xz_xxy_xxxx_0,  \
                             ta2_xz_xxy_xxxx_1,  \
                             ta2_xz_xxy_xxxy_0,  \
                             ta2_xz_xxy_xxxy_1,  \
                             ta2_xz_xxy_xxxz_0,  \
                             ta2_xz_xxy_xxxz_1,  \
                             ta2_xz_xxy_xxy_0,   \
                             ta2_xz_xxy_xxy_1,   \
                             ta2_xz_xxy_xxyy_0,  \
                             ta2_xz_xxy_xxyy_1,  \
                             ta2_xz_xxy_xxyz_0,  \
                             ta2_xz_xxy_xxyz_1,  \
                             ta2_xz_xxy_xxz_0,   \
                             ta2_xz_xxy_xxz_1,   \
                             ta2_xz_xxy_xxzz_0,  \
                             ta2_xz_xxy_xxzz_1,  \
                             ta2_xz_xxy_xyy_0,   \
                             ta2_xz_xxy_xyy_1,   \
                             ta2_xz_xxy_xyyy_0,  \
                             ta2_xz_xxy_xyyy_1,  \
                             ta2_xz_xxy_xyyz_0,  \
                             ta2_xz_xxy_xyyz_1,  \
                             ta2_xz_xxy_xyz_0,   \
                             ta2_xz_xxy_xyz_1,   \
                             ta2_xz_xxy_xyzz_0,  \
                             ta2_xz_xxy_xyzz_1,  \
                             ta2_xz_xxy_xzz_0,   \
                             ta2_xz_xxy_xzz_1,   \
                             ta2_xz_xxy_xzzz_0,  \
                             ta2_xz_xxy_xzzz_1,  \
                             ta2_xz_xxy_zzzz_0,  \
                             ta2_xz_xxy_zzzz_1,  \
                             ta2_xz_xxyy_xxxx_0, \
                             ta2_xz_xxyy_xxxy_0, \
                             ta2_xz_xxyy_xxxz_0, \
                             ta2_xz_xxyy_xxyy_0, \
                             ta2_xz_xxyy_xxyz_0, \
                             ta2_xz_xxyy_xxzz_0, \
                             ta2_xz_xxyy_xyyy_0, \
                             ta2_xz_xxyy_xyyz_0, \
                             ta2_xz_xxyy_xyzz_0, \
                             ta2_xz_xxyy_xzzz_0, \
                             ta2_xz_xxyy_yyyy_0, \
                             ta2_xz_xxyy_yyyz_0, \
                             ta2_xz_xxyy_yyzz_0, \
                             ta2_xz_xxyy_yzzz_0, \
                             ta2_xz_xxyy_zzzz_0, \
                             ta2_xz_xyy_yyyy_0,  \
                             ta2_xz_xyy_yyyy_1,  \
                             ta2_xz_xyy_yyyz_0,  \
                             ta2_xz_xyy_yyyz_1,  \
                             ta2_xz_xyy_yyzz_0,  \
                             ta2_xz_xyy_yyzz_1,  \
                             ta2_xz_xyy_yzzz_0,  \
                             ta2_xz_xyy_yzzz_1,  \
                             ta2_xz_yy_yyyy_0,   \
                             ta2_xz_yy_yyyy_1,   \
                             ta2_xz_yy_yyyz_0,   \
                             ta2_xz_yy_yyyz_1,   \
                             ta2_xz_yy_yyzz_0,   \
                             ta2_xz_yy_yyzz_1,   \
                             ta2_xz_yy_yzzz_0,   \
                             ta2_xz_yy_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxyy_xxxx_0[i] =
            ta2_xz_xx_xxxx_0[i] * fe_0 - ta2_xz_xx_xxxx_1[i] * fe_0 + ta2_xz_xxy_xxxx_0[i] * pa_y[i] - ta2_xz_xxy_xxxx_1[i] * pc_y[i];

        ta2_xz_xxyy_xxxy_0[i] = ta2_xz_xx_xxxy_0[i] * fe_0 - ta2_xz_xx_xxxy_1[i] * fe_0 + ta2_xz_xxy_xxx_0[i] * fe_0 - ta2_xz_xxy_xxx_1[i] * fe_0 +
                                ta2_xz_xxy_xxxy_0[i] * pa_y[i] - ta2_xz_xxy_xxxy_1[i] * pc_y[i];

        ta2_xz_xxyy_xxxz_0[i] =
            ta2_xz_xx_xxxz_0[i] * fe_0 - ta2_xz_xx_xxxz_1[i] * fe_0 + ta2_xz_xxy_xxxz_0[i] * pa_y[i] - ta2_xz_xxy_xxxz_1[i] * pc_y[i];

        ta2_xz_xxyy_xxyy_0[i] = ta2_xz_xx_xxyy_0[i] * fe_0 - ta2_xz_xx_xxyy_1[i] * fe_0 + 2.0 * ta2_xz_xxy_xxy_0[i] * fe_0 -
                                2.0 * ta2_xz_xxy_xxy_1[i] * fe_0 + ta2_xz_xxy_xxyy_0[i] * pa_y[i] - ta2_xz_xxy_xxyy_1[i] * pc_y[i];

        ta2_xz_xxyy_xxyz_0[i] = ta2_xz_xx_xxyz_0[i] * fe_0 - ta2_xz_xx_xxyz_1[i] * fe_0 + ta2_xz_xxy_xxz_0[i] * fe_0 - ta2_xz_xxy_xxz_1[i] * fe_0 +
                                ta2_xz_xxy_xxyz_0[i] * pa_y[i] - ta2_xz_xxy_xxyz_1[i] * pc_y[i];

        ta2_xz_xxyy_xxzz_0[i] =
            ta2_xz_xx_xxzz_0[i] * fe_0 - ta2_xz_xx_xxzz_1[i] * fe_0 + ta2_xz_xxy_xxzz_0[i] * pa_y[i] - ta2_xz_xxy_xxzz_1[i] * pc_y[i];

        ta2_xz_xxyy_xyyy_0[i] = ta2_xz_xx_xyyy_0[i] * fe_0 - ta2_xz_xx_xyyy_1[i] * fe_0 + 3.0 * ta2_xz_xxy_xyy_0[i] * fe_0 -
                                3.0 * ta2_xz_xxy_xyy_1[i] * fe_0 + ta2_xz_xxy_xyyy_0[i] * pa_y[i] - ta2_xz_xxy_xyyy_1[i] * pc_y[i];

        ta2_xz_xxyy_xyyz_0[i] = ta2_xz_xx_xyyz_0[i] * fe_0 - ta2_xz_xx_xyyz_1[i] * fe_0 + 2.0 * ta2_xz_xxy_xyz_0[i] * fe_0 -
                                2.0 * ta2_xz_xxy_xyz_1[i] * fe_0 + ta2_xz_xxy_xyyz_0[i] * pa_y[i] - ta2_xz_xxy_xyyz_1[i] * pc_y[i];

        ta2_xz_xxyy_xyzz_0[i] = ta2_xz_xx_xyzz_0[i] * fe_0 - ta2_xz_xx_xyzz_1[i] * fe_0 + ta2_xz_xxy_xzz_0[i] * fe_0 - ta2_xz_xxy_xzz_1[i] * fe_0 +
                                ta2_xz_xxy_xyzz_0[i] * pa_y[i] - ta2_xz_xxy_xyzz_1[i] * pc_y[i];

        ta2_xz_xxyy_xzzz_0[i] =
            ta2_xz_xx_xzzz_0[i] * fe_0 - ta2_xz_xx_xzzz_1[i] * fe_0 + ta2_xz_xxy_xzzz_0[i] * pa_y[i] - ta2_xz_xxy_xzzz_1[i] * pc_y[i];

        ta2_xz_xxyy_yyyy_0[i] = ta2_xz_yy_yyyy_0[i] * fe_0 - ta2_xz_yy_yyyy_1[i] * fe_0 + ta1_z_xyy_yyyy_1[i] + ta2_xz_xyy_yyyy_0[i] * pa_x[i] -
                                ta2_xz_xyy_yyyy_1[i] * pc_x[i];

        ta2_xz_xxyy_yyyz_0[i] = ta2_xz_yy_yyyz_0[i] * fe_0 - ta2_xz_yy_yyyz_1[i] * fe_0 + ta1_z_xyy_yyyz_1[i] + ta2_xz_xyy_yyyz_0[i] * pa_x[i] -
                                ta2_xz_xyy_yyyz_1[i] * pc_x[i];

        ta2_xz_xxyy_yyzz_0[i] = ta2_xz_yy_yyzz_0[i] * fe_0 - ta2_xz_yy_yyzz_1[i] * fe_0 + ta1_z_xyy_yyzz_1[i] + ta2_xz_xyy_yyzz_0[i] * pa_x[i] -
                                ta2_xz_xyy_yyzz_1[i] * pc_x[i];

        ta2_xz_xxyy_yzzz_0[i] = ta2_xz_yy_yzzz_0[i] * fe_0 - ta2_xz_yy_yzzz_1[i] * fe_0 + ta1_z_xyy_yzzz_1[i] + ta2_xz_xyy_yzzz_0[i] * pa_x[i] -
                                ta2_xz_xyy_yzzz_1[i] * pc_x[i];

        ta2_xz_xxyy_zzzz_0[i] =
            ta2_xz_xx_zzzz_0[i] * fe_0 - ta2_xz_xx_zzzz_1[i] * fe_0 + ta2_xz_xxy_zzzz_0[i] * pa_y[i] - ta2_xz_xxy_zzzz_1[i] * pc_y[i];
    }

    // Set up 510-525 components of targeted buffer : GG

    auto ta2_xz_xxyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 510);

    auto ta2_xz_xxyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 511);

    auto ta2_xz_xxyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 512);

    auto ta2_xz_xxyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 513);

    auto ta2_xz_xxyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 514);

    auto ta2_xz_xxyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 515);

    auto ta2_xz_xxyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 516);

    auto ta2_xz_xxyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 517);

    auto ta2_xz_xxyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 518);

    auto ta2_xz_xxyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 519);

    auto ta2_xz_xxyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 520);

    auto ta2_xz_xxyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 521);

    auto ta2_xz_xxyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 522);

    auto ta2_xz_xxyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 523);

    auto ta2_xz_xxyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 524);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xxy_xxxy_1,   \
                             ta1_x_xxy_xxyy_1,   \
                             ta1_x_xxy_xyyy_1,   \
                             ta1_x_xxy_yyyy_1,   \
                             ta2_xz_xxy_xxxy_0,  \
                             ta2_xz_xxy_xxxy_1,  \
                             ta2_xz_xxy_xxyy_0,  \
                             ta2_xz_xxy_xxyy_1,  \
                             ta2_xz_xxy_xyyy_0,  \
                             ta2_xz_xxy_xyyy_1,  \
                             ta2_xz_xxy_yyyy_0,  \
                             ta2_xz_xxy_yyyy_1,  \
                             ta2_xz_xxyz_xxxx_0, \
                             ta2_xz_xxyz_xxxy_0, \
                             ta2_xz_xxyz_xxxz_0, \
                             ta2_xz_xxyz_xxyy_0, \
                             ta2_xz_xxyz_xxyz_0, \
                             ta2_xz_xxyz_xxzz_0, \
                             ta2_xz_xxyz_xyyy_0, \
                             ta2_xz_xxyz_xyyz_0, \
                             ta2_xz_xxyz_xyzz_0, \
                             ta2_xz_xxyz_xzzz_0, \
                             ta2_xz_xxyz_yyyy_0, \
                             ta2_xz_xxyz_yyyz_0, \
                             ta2_xz_xxyz_yyzz_0, \
                             ta2_xz_xxyz_yzzz_0, \
                             ta2_xz_xxyz_zzzz_0, \
                             ta2_xz_xxz_xxxx_0,  \
                             ta2_xz_xxz_xxxx_1,  \
                             ta2_xz_xxz_xxxz_0,  \
                             ta2_xz_xxz_xxxz_1,  \
                             ta2_xz_xxz_xxyz_0,  \
                             ta2_xz_xxz_xxyz_1,  \
                             ta2_xz_xxz_xxz_0,   \
                             ta2_xz_xxz_xxz_1,   \
                             ta2_xz_xxz_xxzz_0,  \
                             ta2_xz_xxz_xxzz_1,  \
                             ta2_xz_xxz_xyyz_0,  \
                             ta2_xz_xxz_xyyz_1,  \
                             ta2_xz_xxz_xyz_0,   \
                             ta2_xz_xxz_xyz_1,   \
                             ta2_xz_xxz_xyzz_0,  \
                             ta2_xz_xxz_xyzz_1,  \
                             ta2_xz_xxz_xzz_0,   \
                             ta2_xz_xxz_xzz_1,   \
                             ta2_xz_xxz_xzzz_0,  \
                             ta2_xz_xxz_xzzz_1,  \
                             ta2_xz_xxz_yyyz_0,  \
                             ta2_xz_xxz_yyyz_1,  \
                             ta2_xz_xxz_yyz_0,   \
                             ta2_xz_xxz_yyz_1,   \
                             ta2_xz_xxz_yyzz_0,  \
                             ta2_xz_xxz_yyzz_1,  \
                             ta2_xz_xxz_yzz_0,   \
                             ta2_xz_xxz_yzz_1,   \
                             ta2_xz_xxz_yzzz_0,  \
                             ta2_xz_xxz_yzzz_1,  \
                             ta2_xz_xxz_zzz_0,   \
                             ta2_xz_xxz_zzz_1,   \
                             ta2_xz_xxz_zzzz_0,  \
                             ta2_xz_xxz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxyz_xxxx_0[i] = ta2_xz_xxz_xxxx_0[i] * pa_y[i] - ta2_xz_xxz_xxxx_1[i] * pc_y[i];

        ta2_xz_xxyz_xxxy_0[i] = ta1_x_xxy_xxxy_1[i] + ta2_xz_xxy_xxxy_0[i] * pa_z[i] - ta2_xz_xxy_xxxy_1[i] * pc_z[i];

        ta2_xz_xxyz_xxxz_0[i] = ta2_xz_xxz_xxxz_0[i] * pa_y[i] - ta2_xz_xxz_xxxz_1[i] * pc_y[i];

        ta2_xz_xxyz_xxyy_0[i] = ta1_x_xxy_xxyy_1[i] + ta2_xz_xxy_xxyy_0[i] * pa_z[i] - ta2_xz_xxy_xxyy_1[i] * pc_z[i];

        ta2_xz_xxyz_xxyz_0[i] =
            ta2_xz_xxz_xxz_0[i] * fe_0 - ta2_xz_xxz_xxz_1[i] * fe_0 + ta2_xz_xxz_xxyz_0[i] * pa_y[i] - ta2_xz_xxz_xxyz_1[i] * pc_y[i];

        ta2_xz_xxyz_xxzz_0[i] = ta2_xz_xxz_xxzz_0[i] * pa_y[i] - ta2_xz_xxz_xxzz_1[i] * pc_y[i];

        ta2_xz_xxyz_xyyy_0[i] = ta1_x_xxy_xyyy_1[i] + ta2_xz_xxy_xyyy_0[i] * pa_z[i] - ta2_xz_xxy_xyyy_1[i] * pc_z[i];

        ta2_xz_xxyz_xyyz_0[i] =
            2.0 * ta2_xz_xxz_xyz_0[i] * fe_0 - 2.0 * ta2_xz_xxz_xyz_1[i] * fe_0 + ta2_xz_xxz_xyyz_0[i] * pa_y[i] - ta2_xz_xxz_xyyz_1[i] * pc_y[i];

        ta2_xz_xxyz_xyzz_0[i] =
            ta2_xz_xxz_xzz_0[i] * fe_0 - ta2_xz_xxz_xzz_1[i] * fe_0 + ta2_xz_xxz_xyzz_0[i] * pa_y[i] - ta2_xz_xxz_xyzz_1[i] * pc_y[i];

        ta2_xz_xxyz_xzzz_0[i] = ta2_xz_xxz_xzzz_0[i] * pa_y[i] - ta2_xz_xxz_xzzz_1[i] * pc_y[i];

        ta2_xz_xxyz_yyyy_0[i] = ta1_x_xxy_yyyy_1[i] + ta2_xz_xxy_yyyy_0[i] * pa_z[i] - ta2_xz_xxy_yyyy_1[i] * pc_z[i];

        ta2_xz_xxyz_yyyz_0[i] =
            3.0 * ta2_xz_xxz_yyz_0[i] * fe_0 - 3.0 * ta2_xz_xxz_yyz_1[i] * fe_0 + ta2_xz_xxz_yyyz_0[i] * pa_y[i] - ta2_xz_xxz_yyyz_1[i] * pc_y[i];

        ta2_xz_xxyz_yyzz_0[i] =
            2.0 * ta2_xz_xxz_yzz_0[i] * fe_0 - 2.0 * ta2_xz_xxz_yzz_1[i] * fe_0 + ta2_xz_xxz_yyzz_0[i] * pa_y[i] - ta2_xz_xxz_yyzz_1[i] * pc_y[i];

        ta2_xz_xxyz_yzzz_0[i] =
            ta2_xz_xxz_zzz_0[i] * fe_0 - ta2_xz_xxz_zzz_1[i] * fe_0 + ta2_xz_xxz_yzzz_0[i] * pa_y[i] - ta2_xz_xxz_yzzz_1[i] * pc_y[i];

        ta2_xz_xxyz_zzzz_0[i] = ta2_xz_xxz_zzzz_0[i] * pa_y[i] - ta2_xz_xxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 525-540 components of targeted buffer : GG

    auto ta2_xz_xxzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 525);

    auto ta2_xz_xxzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 526);

    auto ta2_xz_xxzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 527);

    auto ta2_xz_xxzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 528);

    auto ta2_xz_xxzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 529);

    auto ta2_xz_xxzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 530);

    auto ta2_xz_xxzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 531);

    auto ta2_xz_xxzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 532);

    auto ta2_xz_xxzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 533);

    auto ta2_xz_xxzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 534);

    auto ta2_xz_xxzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 535);

    auto ta2_xz_xxzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 536);

    auto ta2_xz_xxzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 537);

    auto ta2_xz_xxzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 538);

    auto ta2_xz_xxzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 539);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_x_xxz_xxxx_1,   \
                             ta1_x_xxz_xxxy_1,   \
                             ta1_x_xxz_xxyy_1,   \
                             ta1_x_xxz_xyyy_1,   \
                             ta1_z_xzz_xxxz_1,   \
                             ta1_z_xzz_xxyz_1,   \
                             ta1_z_xzz_xxzz_1,   \
                             ta1_z_xzz_xyyz_1,   \
                             ta1_z_xzz_xyzz_1,   \
                             ta1_z_xzz_xzzz_1,   \
                             ta1_z_xzz_yyyy_1,   \
                             ta1_z_xzz_yyyz_1,   \
                             ta1_z_xzz_yyzz_1,   \
                             ta1_z_xzz_yzzz_1,   \
                             ta1_z_xzz_zzzz_1,   \
                             ta2_xz_xx_xxxx_0,   \
                             ta2_xz_xx_xxxx_1,   \
                             ta2_xz_xx_xxxy_0,   \
                             ta2_xz_xx_xxxy_1,   \
                             ta2_xz_xx_xxyy_0,   \
                             ta2_xz_xx_xxyy_1,   \
                             ta2_xz_xx_xyyy_0,   \
                             ta2_xz_xx_xyyy_1,   \
                             ta2_xz_xxz_xxxx_0,  \
                             ta2_xz_xxz_xxxx_1,  \
                             ta2_xz_xxz_xxxy_0,  \
                             ta2_xz_xxz_xxxy_1,  \
                             ta2_xz_xxz_xxyy_0,  \
                             ta2_xz_xxz_xxyy_1,  \
                             ta2_xz_xxz_xyyy_0,  \
                             ta2_xz_xxz_xyyy_1,  \
                             ta2_xz_xxzz_xxxx_0, \
                             ta2_xz_xxzz_xxxy_0, \
                             ta2_xz_xxzz_xxxz_0, \
                             ta2_xz_xxzz_xxyy_0, \
                             ta2_xz_xxzz_xxyz_0, \
                             ta2_xz_xxzz_xxzz_0, \
                             ta2_xz_xxzz_xyyy_0, \
                             ta2_xz_xxzz_xyyz_0, \
                             ta2_xz_xxzz_xyzz_0, \
                             ta2_xz_xxzz_xzzz_0, \
                             ta2_xz_xxzz_yyyy_0, \
                             ta2_xz_xxzz_yyyz_0, \
                             ta2_xz_xxzz_yyzz_0, \
                             ta2_xz_xxzz_yzzz_0, \
                             ta2_xz_xxzz_zzzz_0, \
                             ta2_xz_xzz_xxxz_0,  \
                             ta2_xz_xzz_xxxz_1,  \
                             ta2_xz_xzz_xxyz_0,  \
                             ta2_xz_xzz_xxyz_1,  \
                             ta2_xz_xzz_xxz_0,   \
                             ta2_xz_xzz_xxz_1,   \
                             ta2_xz_xzz_xxzz_0,  \
                             ta2_xz_xzz_xxzz_1,  \
                             ta2_xz_xzz_xyyz_0,  \
                             ta2_xz_xzz_xyyz_1,  \
                             ta2_xz_xzz_xyz_0,   \
                             ta2_xz_xzz_xyz_1,   \
                             ta2_xz_xzz_xyzz_0,  \
                             ta2_xz_xzz_xyzz_1,  \
                             ta2_xz_xzz_xzz_0,   \
                             ta2_xz_xzz_xzz_1,   \
                             ta2_xz_xzz_xzzz_0,  \
                             ta2_xz_xzz_xzzz_1,  \
                             ta2_xz_xzz_yyyy_0,  \
                             ta2_xz_xzz_yyyy_1,  \
                             ta2_xz_xzz_yyyz_0,  \
                             ta2_xz_xzz_yyyz_1,  \
                             ta2_xz_xzz_yyz_0,   \
                             ta2_xz_xzz_yyz_1,   \
                             ta2_xz_xzz_yyzz_0,  \
                             ta2_xz_xzz_yyzz_1,  \
                             ta2_xz_xzz_yzz_0,   \
                             ta2_xz_xzz_yzz_1,   \
                             ta2_xz_xzz_yzzz_0,  \
                             ta2_xz_xzz_yzzz_1,  \
                             ta2_xz_xzz_zzz_0,   \
                             ta2_xz_xzz_zzz_1,   \
                             ta2_xz_xzz_zzzz_0,  \
                             ta2_xz_xzz_zzzz_1,  \
                             ta2_xz_zz_xxxz_0,   \
                             ta2_xz_zz_xxxz_1,   \
                             ta2_xz_zz_xxyz_0,   \
                             ta2_xz_zz_xxyz_1,   \
                             ta2_xz_zz_xxzz_0,   \
                             ta2_xz_zz_xxzz_1,   \
                             ta2_xz_zz_xyyz_0,   \
                             ta2_xz_zz_xyyz_1,   \
                             ta2_xz_zz_xyzz_0,   \
                             ta2_xz_zz_xyzz_1,   \
                             ta2_xz_zz_xzzz_0,   \
                             ta2_xz_zz_xzzz_1,   \
                             ta2_xz_zz_yyyy_0,   \
                             ta2_xz_zz_yyyy_1,   \
                             ta2_xz_zz_yyyz_0,   \
                             ta2_xz_zz_yyyz_1,   \
                             ta2_xz_zz_yyzz_0,   \
                             ta2_xz_zz_yyzz_1,   \
                             ta2_xz_zz_yzzz_0,   \
                             ta2_xz_zz_yzzz_1,   \
                             ta2_xz_zz_zzzz_0,   \
                             ta2_xz_zz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xxzz_xxxx_0[i] = ta2_xz_xx_xxxx_0[i] * fe_0 - ta2_xz_xx_xxxx_1[i] * fe_0 + ta1_x_xxz_xxxx_1[i] + ta2_xz_xxz_xxxx_0[i] * pa_z[i] -
                                ta2_xz_xxz_xxxx_1[i] * pc_z[i];

        ta2_xz_xxzz_xxxy_0[i] = ta2_xz_xx_xxxy_0[i] * fe_0 - ta2_xz_xx_xxxy_1[i] * fe_0 + ta1_x_xxz_xxxy_1[i] + ta2_xz_xxz_xxxy_0[i] * pa_z[i] -
                                ta2_xz_xxz_xxxy_1[i] * pc_z[i];

        ta2_xz_xxzz_xxxz_0[i] = ta2_xz_zz_xxxz_0[i] * fe_0 - ta2_xz_zz_xxxz_1[i] * fe_0 + 3.0 * ta2_xz_xzz_xxz_0[i] * fe_0 -
                                3.0 * ta2_xz_xzz_xxz_1[i] * fe_0 + ta1_z_xzz_xxxz_1[i] + ta2_xz_xzz_xxxz_0[i] * pa_x[i] -
                                ta2_xz_xzz_xxxz_1[i] * pc_x[i];

        ta2_xz_xxzz_xxyy_0[i] = ta2_xz_xx_xxyy_0[i] * fe_0 - ta2_xz_xx_xxyy_1[i] * fe_0 + ta1_x_xxz_xxyy_1[i] + ta2_xz_xxz_xxyy_0[i] * pa_z[i] -
                                ta2_xz_xxz_xxyy_1[i] * pc_z[i];

        ta2_xz_xxzz_xxyz_0[i] = ta2_xz_zz_xxyz_0[i] * fe_0 - ta2_xz_zz_xxyz_1[i] * fe_0 + 2.0 * ta2_xz_xzz_xyz_0[i] * fe_0 -
                                2.0 * ta2_xz_xzz_xyz_1[i] * fe_0 + ta1_z_xzz_xxyz_1[i] + ta2_xz_xzz_xxyz_0[i] * pa_x[i] -
                                ta2_xz_xzz_xxyz_1[i] * pc_x[i];

        ta2_xz_xxzz_xxzz_0[i] = ta2_xz_zz_xxzz_0[i] * fe_0 - ta2_xz_zz_xxzz_1[i] * fe_0 + 2.0 * ta2_xz_xzz_xzz_0[i] * fe_0 -
                                2.0 * ta2_xz_xzz_xzz_1[i] * fe_0 + ta1_z_xzz_xxzz_1[i] + ta2_xz_xzz_xxzz_0[i] * pa_x[i] -
                                ta2_xz_xzz_xxzz_1[i] * pc_x[i];

        ta2_xz_xxzz_xyyy_0[i] = ta2_xz_xx_xyyy_0[i] * fe_0 - ta2_xz_xx_xyyy_1[i] * fe_0 + ta1_x_xxz_xyyy_1[i] + ta2_xz_xxz_xyyy_0[i] * pa_z[i] -
                                ta2_xz_xxz_xyyy_1[i] * pc_z[i];

        ta2_xz_xxzz_xyyz_0[i] = ta2_xz_zz_xyyz_0[i] * fe_0 - ta2_xz_zz_xyyz_1[i] * fe_0 + ta2_xz_xzz_yyz_0[i] * fe_0 - ta2_xz_xzz_yyz_1[i] * fe_0 +
                                ta1_z_xzz_xyyz_1[i] + ta2_xz_xzz_xyyz_0[i] * pa_x[i] - ta2_xz_xzz_xyyz_1[i] * pc_x[i];

        ta2_xz_xxzz_xyzz_0[i] = ta2_xz_zz_xyzz_0[i] * fe_0 - ta2_xz_zz_xyzz_1[i] * fe_0 + ta2_xz_xzz_yzz_0[i] * fe_0 - ta2_xz_xzz_yzz_1[i] * fe_0 +
                                ta1_z_xzz_xyzz_1[i] + ta2_xz_xzz_xyzz_0[i] * pa_x[i] - ta2_xz_xzz_xyzz_1[i] * pc_x[i];

        ta2_xz_xxzz_xzzz_0[i] = ta2_xz_zz_xzzz_0[i] * fe_0 - ta2_xz_zz_xzzz_1[i] * fe_0 + ta2_xz_xzz_zzz_0[i] * fe_0 - ta2_xz_xzz_zzz_1[i] * fe_0 +
                                ta1_z_xzz_xzzz_1[i] + ta2_xz_xzz_xzzz_0[i] * pa_x[i] - ta2_xz_xzz_xzzz_1[i] * pc_x[i];

        ta2_xz_xxzz_yyyy_0[i] = ta2_xz_zz_yyyy_0[i] * fe_0 - ta2_xz_zz_yyyy_1[i] * fe_0 + ta1_z_xzz_yyyy_1[i] + ta2_xz_xzz_yyyy_0[i] * pa_x[i] -
                                ta2_xz_xzz_yyyy_1[i] * pc_x[i];

        ta2_xz_xxzz_yyyz_0[i] = ta2_xz_zz_yyyz_0[i] * fe_0 - ta2_xz_zz_yyyz_1[i] * fe_0 + ta1_z_xzz_yyyz_1[i] + ta2_xz_xzz_yyyz_0[i] * pa_x[i] -
                                ta2_xz_xzz_yyyz_1[i] * pc_x[i];

        ta2_xz_xxzz_yyzz_0[i] = ta2_xz_zz_yyzz_0[i] * fe_0 - ta2_xz_zz_yyzz_1[i] * fe_0 + ta1_z_xzz_yyzz_1[i] + ta2_xz_xzz_yyzz_0[i] * pa_x[i] -
                                ta2_xz_xzz_yyzz_1[i] * pc_x[i];

        ta2_xz_xxzz_yzzz_0[i] = ta2_xz_zz_yzzz_0[i] * fe_0 - ta2_xz_zz_yzzz_1[i] * fe_0 + ta1_z_xzz_yzzz_1[i] + ta2_xz_xzz_yzzz_0[i] * pa_x[i] -
                                ta2_xz_xzz_yzzz_1[i] * pc_x[i];

        ta2_xz_xxzz_zzzz_0[i] = ta2_xz_zz_zzzz_0[i] * fe_0 - ta2_xz_zz_zzzz_1[i] * fe_0 + ta1_z_xzz_zzzz_1[i] + ta2_xz_xzz_zzzz_0[i] * pa_x[i] -
                                ta2_xz_xzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 540-555 components of targeted buffer : GG

    auto ta2_xz_xyyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 540);

    auto ta2_xz_xyyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 541);

    auto ta2_xz_xyyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 542);

    auto ta2_xz_xyyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 543);

    auto ta2_xz_xyyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 544);

    auto ta2_xz_xyyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 545);

    auto ta2_xz_xyyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 546);

    auto ta2_xz_xyyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 547);

    auto ta2_xz_xyyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 548);

    auto ta2_xz_xyyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 549);

    auto ta2_xz_xyyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 550);

    auto ta2_xz_xyyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 551);

    auto ta2_xz_xyyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 552);

    auto ta2_xz_xyyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 553);

    auto ta2_xz_xyyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 554);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_yyy_xxxy_1,   \
                             ta1_z_yyy_xxyy_1,   \
                             ta1_z_yyy_xxyz_1,   \
                             ta1_z_yyy_xyyy_1,   \
                             ta1_z_yyy_xyyz_1,   \
                             ta1_z_yyy_xyzz_1,   \
                             ta1_z_yyy_yyyy_1,   \
                             ta1_z_yyy_yyyz_1,   \
                             ta1_z_yyy_yyzz_1,   \
                             ta1_z_yyy_yzzz_1,   \
                             ta1_z_yyy_zzzz_1,   \
                             ta2_xz_xy_xxxx_0,   \
                             ta2_xz_xy_xxxx_1,   \
                             ta2_xz_xy_xxxz_0,   \
                             ta2_xz_xy_xxxz_1,   \
                             ta2_xz_xy_xxzz_0,   \
                             ta2_xz_xy_xxzz_1,   \
                             ta2_xz_xy_xzzz_0,   \
                             ta2_xz_xy_xzzz_1,   \
                             ta2_xz_xyy_xxxx_0,  \
                             ta2_xz_xyy_xxxx_1,  \
                             ta2_xz_xyy_xxxz_0,  \
                             ta2_xz_xyy_xxxz_1,  \
                             ta2_xz_xyy_xxzz_0,  \
                             ta2_xz_xyy_xxzz_1,  \
                             ta2_xz_xyy_xzzz_0,  \
                             ta2_xz_xyy_xzzz_1,  \
                             ta2_xz_xyyy_xxxx_0, \
                             ta2_xz_xyyy_xxxy_0, \
                             ta2_xz_xyyy_xxxz_0, \
                             ta2_xz_xyyy_xxyy_0, \
                             ta2_xz_xyyy_xxyz_0, \
                             ta2_xz_xyyy_xxzz_0, \
                             ta2_xz_xyyy_xyyy_0, \
                             ta2_xz_xyyy_xyyz_0, \
                             ta2_xz_xyyy_xyzz_0, \
                             ta2_xz_xyyy_xzzz_0, \
                             ta2_xz_xyyy_yyyy_0, \
                             ta2_xz_xyyy_yyyz_0, \
                             ta2_xz_xyyy_yyzz_0, \
                             ta2_xz_xyyy_yzzz_0, \
                             ta2_xz_xyyy_zzzz_0, \
                             ta2_xz_yyy_xxxy_0,  \
                             ta2_xz_yyy_xxxy_1,  \
                             ta2_xz_yyy_xxy_0,   \
                             ta2_xz_yyy_xxy_1,   \
                             ta2_xz_yyy_xxyy_0,  \
                             ta2_xz_yyy_xxyy_1,  \
                             ta2_xz_yyy_xxyz_0,  \
                             ta2_xz_yyy_xxyz_1,  \
                             ta2_xz_yyy_xyy_0,   \
                             ta2_xz_yyy_xyy_1,   \
                             ta2_xz_yyy_xyyy_0,  \
                             ta2_xz_yyy_xyyy_1,  \
                             ta2_xz_yyy_xyyz_0,  \
                             ta2_xz_yyy_xyyz_1,  \
                             ta2_xz_yyy_xyz_0,   \
                             ta2_xz_yyy_xyz_1,   \
                             ta2_xz_yyy_xyzz_0,  \
                             ta2_xz_yyy_xyzz_1,  \
                             ta2_xz_yyy_yyy_0,   \
                             ta2_xz_yyy_yyy_1,   \
                             ta2_xz_yyy_yyyy_0,  \
                             ta2_xz_yyy_yyyy_1,  \
                             ta2_xz_yyy_yyyz_0,  \
                             ta2_xz_yyy_yyyz_1,  \
                             ta2_xz_yyy_yyz_0,   \
                             ta2_xz_yyy_yyz_1,   \
                             ta2_xz_yyy_yyzz_0,  \
                             ta2_xz_yyy_yyzz_1,  \
                             ta2_xz_yyy_yzz_0,   \
                             ta2_xz_yyy_yzz_1,   \
                             ta2_xz_yyy_yzzz_0,  \
                             ta2_xz_yyy_yzzz_1,  \
                             ta2_xz_yyy_zzzz_0,  \
                             ta2_xz_yyy_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xyyy_xxxx_0[i] =
            2.0 * ta2_xz_xy_xxxx_0[i] * fe_0 - 2.0 * ta2_xz_xy_xxxx_1[i] * fe_0 + ta2_xz_xyy_xxxx_0[i] * pa_y[i] - ta2_xz_xyy_xxxx_1[i] * pc_y[i];

        ta2_xz_xyyy_xxxy_0[i] = 3.0 * ta2_xz_yyy_xxy_0[i] * fe_0 - 3.0 * ta2_xz_yyy_xxy_1[i] * fe_0 + ta1_z_yyy_xxxy_1[i] +
                                ta2_xz_yyy_xxxy_0[i] * pa_x[i] - ta2_xz_yyy_xxxy_1[i] * pc_x[i];

        ta2_xz_xyyy_xxxz_0[i] =
            2.0 * ta2_xz_xy_xxxz_0[i] * fe_0 - 2.0 * ta2_xz_xy_xxxz_1[i] * fe_0 + ta2_xz_xyy_xxxz_0[i] * pa_y[i] - ta2_xz_xyy_xxxz_1[i] * pc_y[i];

        ta2_xz_xyyy_xxyy_0[i] = 2.0 * ta2_xz_yyy_xyy_0[i] * fe_0 - 2.0 * ta2_xz_yyy_xyy_1[i] * fe_0 + ta1_z_yyy_xxyy_1[i] +
                                ta2_xz_yyy_xxyy_0[i] * pa_x[i] - ta2_xz_yyy_xxyy_1[i] * pc_x[i];

        ta2_xz_xyyy_xxyz_0[i] = 2.0 * ta2_xz_yyy_xyz_0[i] * fe_0 - 2.0 * ta2_xz_yyy_xyz_1[i] * fe_0 + ta1_z_yyy_xxyz_1[i] +
                                ta2_xz_yyy_xxyz_0[i] * pa_x[i] - ta2_xz_yyy_xxyz_1[i] * pc_x[i];

        ta2_xz_xyyy_xxzz_0[i] =
            2.0 * ta2_xz_xy_xxzz_0[i] * fe_0 - 2.0 * ta2_xz_xy_xxzz_1[i] * fe_0 + ta2_xz_xyy_xxzz_0[i] * pa_y[i] - ta2_xz_xyy_xxzz_1[i] * pc_y[i];

        ta2_xz_xyyy_xyyy_0[i] = ta2_xz_yyy_yyy_0[i] * fe_0 - ta2_xz_yyy_yyy_1[i] * fe_0 + ta1_z_yyy_xyyy_1[i] + ta2_xz_yyy_xyyy_0[i] * pa_x[i] -
                                ta2_xz_yyy_xyyy_1[i] * pc_x[i];

        ta2_xz_xyyy_xyyz_0[i] = ta2_xz_yyy_yyz_0[i] * fe_0 - ta2_xz_yyy_yyz_1[i] * fe_0 + ta1_z_yyy_xyyz_1[i] + ta2_xz_yyy_xyyz_0[i] * pa_x[i] -
                                ta2_xz_yyy_xyyz_1[i] * pc_x[i];

        ta2_xz_xyyy_xyzz_0[i] = ta2_xz_yyy_yzz_0[i] * fe_0 - ta2_xz_yyy_yzz_1[i] * fe_0 + ta1_z_yyy_xyzz_1[i] + ta2_xz_yyy_xyzz_0[i] * pa_x[i] -
                                ta2_xz_yyy_xyzz_1[i] * pc_x[i];

        ta2_xz_xyyy_xzzz_0[i] =
            2.0 * ta2_xz_xy_xzzz_0[i] * fe_0 - 2.0 * ta2_xz_xy_xzzz_1[i] * fe_0 + ta2_xz_xyy_xzzz_0[i] * pa_y[i] - ta2_xz_xyy_xzzz_1[i] * pc_y[i];

        ta2_xz_xyyy_yyyy_0[i] = ta1_z_yyy_yyyy_1[i] + ta2_xz_yyy_yyyy_0[i] * pa_x[i] - ta2_xz_yyy_yyyy_1[i] * pc_x[i];

        ta2_xz_xyyy_yyyz_0[i] = ta1_z_yyy_yyyz_1[i] + ta2_xz_yyy_yyyz_0[i] * pa_x[i] - ta2_xz_yyy_yyyz_1[i] * pc_x[i];

        ta2_xz_xyyy_yyzz_0[i] = ta1_z_yyy_yyzz_1[i] + ta2_xz_yyy_yyzz_0[i] * pa_x[i] - ta2_xz_yyy_yyzz_1[i] * pc_x[i];

        ta2_xz_xyyy_yzzz_0[i] = ta1_z_yyy_yzzz_1[i] + ta2_xz_yyy_yzzz_0[i] * pa_x[i] - ta2_xz_yyy_yzzz_1[i] * pc_x[i];

        ta2_xz_xyyy_zzzz_0[i] = ta1_z_yyy_zzzz_1[i] + ta2_xz_yyy_zzzz_0[i] * pa_x[i] - ta2_xz_yyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 555-570 components of targeted buffer : GG

    auto ta2_xz_xyyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 555);

    auto ta2_xz_xyyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 556);

    auto ta2_xz_xyyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 557);

    auto ta2_xz_xyyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 558);

    auto ta2_xz_xyyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 559);

    auto ta2_xz_xyyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 560);

    auto ta2_xz_xyyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 561);

    auto ta2_xz_xyyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 562);

    auto ta2_xz_xyyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 563);

    auto ta2_xz_xyyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 564);

    auto ta2_xz_xyyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 565);

    auto ta2_xz_xyyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 566);

    auto ta2_xz_xyyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 567);

    auto ta2_xz_xyyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 568);

    auto ta2_xz_xyyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 569);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xyy_xxxx_1,   \
                             ta1_x_xyy_xxxy_1,   \
                             ta1_x_xyy_xxyy_1,   \
                             ta1_x_xyy_xyyy_1,   \
                             ta1_z_yyz_xxyz_1,   \
                             ta1_z_yyz_xyyz_1,   \
                             ta1_z_yyz_xyzz_1,   \
                             ta1_z_yyz_yyyy_1,   \
                             ta1_z_yyz_yyyz_1,   \
                             ta1_z_yyz_yyzz_1,   \
                             ta1_z_yyz_yzzz_1,   \
                             ta1_z_yyz_zzzz_1,   \
                             ta2_xz_xyy_xxxx_0,  \
                             ta2_xz_xyy_xxxx_1,  \
                             ta2_xz_xyy_xxxy_0,  \
                             ta2_xz_xyy_xxxy_1,  \
                             ta2_xz_xyy_xxyy_0,  \
                             ta2_xz_xyy_xxyy_1,  \
                             ta2_xz_xyy_xyyy_0,  \
                             ta2_xz_xyy_xyyy_1,  \
                             ta2_xz_xyyz_xxxx_0, \
                             ta2_xz_xyyz_xxxy_0, \
                             ta2_xz_xyyz_xxxz_0, \
                             ta2_xz_xyyz_xxyy_0, \
                             ta2_xz_xyyz_xxyz_0, \
                             ta2_xz_xyyz_xxzz_0, \
                             ta2_xz_xyyz_xyyy_0, \
                             ta2_xz_xyyz_xyyz_0, \
                             ta2_xz_xyyz_xyzz_0, \
                             ta2_xz_xyyz_xzzz_0, \
                             ta2_xz_xyyz_yyyy_0, \
                             ta2_xz_xyyz_yyyz_0, \
                             ta2_xz_xyyz_yyzz_0, \
                             ta2_xz_xyyz_yzzz_0, \
                             ta2_xz_xyyz_zzzz_0, \
                             ta2_xz_xyz_xxxz_0,  \
                             ta2_xz_xyz_xxxz_1,  \
                             ta2_xz_xyz_xxzz_0,  \
                             ta2_xz_xyz_xxzz_1,  \
                             ta2_xz_xyz_xzzz_0,  \
                             ta2_xz_xyz_xzzz_1,  \
                             ta2_xz_xz_xxxz_0,   \
                             ta2_xz_xz_xxxz_1,   \
                             ta2_xz_xz_xxzz_0,   \
                             ta2_xz_xz_xxzz_1,   \
                             ta2_xz_xz_xzzz_0,   \
                             ta2_xz_xz_xzzz_1,   \
                             ta2_xz_yyz_xxyz_0,  \
                             ta2_xz_yyz_xxyz_1,  \
                             ta2_xz_yyz_xyyz_0,  \
                             ta2_xz_yyz_xyyz_1,  \
                             ta2_xz_yyz_xyz_0,   \
                             ta2_xz_yyz_xyz_1,   \
                             ta2_xz_yyz_xyzz_0,  \
                             ta2_xz_yyz_xyzz_1,  \
                             ta2_xz_yyz_yyyy_0,  \
                             ta2_xz_yyz_yyyy_1,  \
                             ta2_xz_yyz_yyyz_0,  \
                             ta2_xz_yyz_yyyz_1,  \
                             ta2_xz_yyz_yyz_0,   \
                             ta2_xz_yyz_yyz_1,   \
                             ta2_xz_yyz_yyzz_0,  \
                             ta2_xz_yyz_yyzz_1,  \
                             ta2_xz_yyz_yzz_0,   \
                             ta2_xz_yyz_yzz_1,   \
                             ta2_xz_yyz_yzzz_0,  \
                             ta2_xz_yyz_yzzz_1,  \
                             ta2_xz_yyz_zzzz_0,  \
                             ta2_xz_yyz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xyyz_xxxx_0[i] = ta1_x_xyy_xxxx_1[i] + ta2_xz_xyy_xxxx_0[i] * pa_z[i] - ta2_xz_xyy_xxxx_1[i] * pc_z[i];

        ta2_xz_xyyz_xxxy_0[i] = ta1_x_xyy_xxxy_1[i] + ta2_xz_xyy_xxxy_0[i] * pa_z[i] - ta2_xz_xyy_xxxy_1[i] * pc_z[i];

        ta2_xz_xyyz_xxxz_0[i] =
            ta2_xz_xz_xxxz_0[i] * fe_0 - ta2_xz_xz_xxxz_1[i] * fe_0 + ta2_xz_xyz_xxxz_0[i] * pa_y[i] - ta2_xz_xyz_xxxz_1[i] * pc_y[i];

        ta2_xz_xyyz_xxyy_0[i] = ta1_x_xyy_xxyy_1[i] + ta2_xz_xyy_xxyy_0[i] * pa_z[i] - ta2_xz_xyy_xxyy_1[i] * pc_z[i];

        ta2_xz_xyyz_xxyz_0[i] = 2.0 * ta2_xz_yyz_xyz_0[i] * fe_0 - 2.0 * ta2_xz_yyz_xyz_1[i] * fe_0 + ta1_z_yyz_xxyz_1[i] +
                                ta2_xz_yyz_xxyz_0[i] * pa_x[i] - ta2_xz_yyz_xxyz_1[i] * pc_x[i];

        ta2_xz_xyyz_xxzz_0[i] =
            ta2_xz_xz_xxzz_0[i] * fe_0 - ta2_xz_xz_xxzz_1[i] * fe_0 + ta2_xz_xyz_xxzz_0[i] * pa_y[i] - ta2_xz_xyz_xxzz_1[i] * pc_y[i];

        ta2_xz_xyyz_xyyy_0[i] = ta1_x_xyy_xyyy_1[i] + ta2_xz_xyy_xyyy_0[i] * pa_z[i] - ta2_xz_xyy_xyyy_1[i] * pc_z[i];

        ta2_xz_xyyz_xyyz_0[i] = ta2_xz_yyz_yyz_0[i] * fe_0 - ta2_xz_yyz_yyz_1[i] * fe_0 + ta1_z_yyz_xyyz_1[i] + ta2_xz_yyz_xyyz_0[i] * pa_x[i] -
                                ta2_xz_yyz_xyyz_1[i] * pc_x[i];

        ta2_xz_xyyz_xyzz_0[i] = ta2_xz_yyz_yzz_0[i] * fe_0 - ta2_xz_yyz_yzz_1[i] * fe_0 + ta1_z_yyz_xyzz_1[i] + ta2_xz_yyz_xyzz_0[i] * pa_x[i] -
                                ta2_xz_yyz_xyzz_1[i] * pc_x[i];

        ta2_xz_xyyz_xzzz_0[i] =
            ta2_xz_xz_xzzz_0[i] * fe_0 - ta2_xz_xz_xzzz_1[i] * fe_0 + ta2_xz_xyz_xzzz_0[i] * pa_y[i] - ta2_xz_xyz_xzzz_1[i] * pc_y[i];

        ta2_xz_xyyz_yyyy_0[i] = ta1_z_yyz_yyyy_1[i] + ta2_xz_yyz_yyyy_0[i] * pa_x[i] - ta2_xz_yyz_yyyy_1[i] * pc_x[i];

        ta2_xz_xyyz_yyyz_0[i] = ta1_z_yyz_yyyz_1[i] + ta2_xz_yyz_yyyz_0[i] * pa_x[i] - ta2_xz_yyz_yyyz_1[i] * pc_x[i];

        ta2_xz_xyyz_yyzz_0[i] = ta1_z_yyz_yyzz_1[i] + ta2_xz_yyz_yyzz_0[i] * pa_x[i] - ta2_xz_yyz_yyzz_1[i] * pc_x[i];

        ta2_xz_xyyz_yzzz_0[i] = ta1_z_yyz_yzzz_1[i] + ta2_xz_yyz_yzzz_0[i] * pa_x[i] - ta2_xz_yyz_yzzz_1[i] * pc_x[i];

        ta2_xz_xyyz_zzzz_0[i] = ta1_z_yyz_zzzz_1[i] + ta2_xz_yyz_zzzz_0[i] * pa_x[i] - ta2_xz_yyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 570-585 components of targeted buffer : GG

    auto ta2_xz_xyzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 570);

    auto ta2_xz_xyzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 571);

    auto ta2_xz_xyzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 572);

    auto ta2_xz_xyzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 573);

    auto ta2_xz_xyzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 574);

    auto ta2_xz_xyzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 575);

    auto ta2_xz_xyzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 576);

    auto ta2_xz_xyzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 577);

    auto ta2_xz_xyzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 578);

    auto ta2_xz_xyzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 579);

    auto ta2_xz_xyzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 580);

    auto ta2_xz_xyzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 581);

    auto ta2_xz_xyzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 582);

    auto ta2_xz_xyzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 583);

    auto ta2_xz_xyzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 584);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_yzz_yyyy_1,   \
                             ta1_z_yzz_yyyz_1,   \
                             ta1_z_yzz_yyzz_1,   \
                             ta1_z_yzz_yzzz_1,   \
                             ta2_xz_xyzz_xxxx_0, \
                             ta2_xz_xyzz_xxxy_0, \
                             ta2_xz_xyzz_xxxz_0, \
                             ta2_xz_xyzz_xxyy_0, \
                             ta2_xz_xyzz_xxyz_0, \
                             ta2_xz_xyzz_xxzz_0, \
                             ta2_xz_xyzz_xyyy_0, \
                             ta2_xz_xyzz_xyyz_0, \
                             ta2_xz_xyzz_xyzz_0, \
                             ta2_xz_xyzz_xzzz_0, \
                             ta2_xz_xyzz_yyyy_0, \
                             ta2_xz_xyzz_yyyz_0, \
                             ta2_xz_xyzz_yyzz_0, \
                             ta2_xz_xyzz_yzzz_0, \
                             ta2_xz_xyzz_zzzz_0, \
                             ta2_xz_xzz_xxx_0,   \
                             ta2_xz_xzz_xxx_1,   \
                             ta2_xz_xzz_xxxx_0,  \
                             ta2_xz_xzz_xxxx_1,  \
                             ta2_xz_xzz_xxxy_0,  \
                             ta2_xz_xzz_xxxy_1,  \
                             ta2_xz_xzz_xxxz_0,  \
                             ta2_xz_xzz_xxxz_1,  \
                             ta2_xz_xzz_xxy_0,   \
                             ta2_xz_xzz_xxy_1,   \
                             ta2_xz_xzz_xxyy_0,  \
                             ta2_xz_xzz_xxyy_1,  \
                             ta2_xz_xzz_xxyz_0,  \
                             ta2_xz_xzz_xxyz_1,  \
                             ta2_xz_xzz_xxz_0,   \
                             ta2_xz_xzz_xxz_1,   \
                             ta2_xz_xzz_xxzz_0,  \
                             ta2_xz_xzz_xxzz_1,  \
                             ta2_xz_xzz_xyy_0,   \
                             ta2_xz_xzz_xyy_1,   \
                             ta2_xz_xzz_xyyy_0,  \
                             ta2_xz_xzz_xyyy_1,  \
                             ta2_xz_xzz_xyyz_0,  \
                             ta2_xz_xzz_xyyz_1,  \
                             ta2_xz_xzz_xyz_0,   \
                             ta2_xz_xzz_xyz_1,   \
                             ta2_xz_xzz_xyzz_0,  \
                             ta2_xz_xzz_xyzz_1,  \
                             ta2_xz_xzz_xzz_0,   \
                             ta2_xz_xzz_xzz_1,   \
                             ta2_xz_xzz_xzzz_0,  \
                             ta2_xz_xzz_xzzz_1,  \
                             ta2_xz_xzz_zzzz_0,  \
                             ta2_xz_xzz_zzzz_1,  \
                             ta2_xz_yzz_yyyy_0,  \
                             ta2_xz_yzz_yyyy_1,  \
                             ta2_xz_yzz_yyyz_0,  \
                             ta2_xz_yzz_yyyz_1,  \
                             ta2_xz_yzz_yyzz_0,  \
                             ta2_xz_yzz_yyzz_1,  \
                             ta2_xz_yzz_yzzz_0,  \
                             ta2_xz_yzz_yzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xyzz_xxxx_0[i] = ta2_xz_xzz_xxxx_0[i] * pa_y[i] - ta2_xz_xzz_xxxx_1[i] * pc_y[i];

        ta2_xz_xyzz_xxxy_0[i] =
            ta2_xz_xzz_xxx_0[i] * fe_0 - ta2_xz_xzz_xxx_1[i] * fe_0 + ta2_xz_xzz_xxxy_0[i] * pa_y[i] - ta2_xz_xzz_xxxy_1[i] * pc_y[i];

        ta2_xz_xyzz_xxxz_0[i] = ta2_xz_xzz_xxxz_0[i] * pa_y[i] - ta2_xz_xzz_xxxz_1[i] * pc_y[i];

        ta2_xz_xyzz_xxyy_0[i] =
            2.0 * ta2_xz_xzz_xxy_0[i] * fe_0 - 2.0 * ta2_xz_xzz_xxy_1[i] * fe_0 + ta2_xz_xzz_xxyy_0[i] * pa_y[i] - ta2_xz_xzz_xxyy_1[i] * pc_y[i];

        ta2_xz_xyzz_xxyz_0[i] =
            ta2_xz_xzz_xxz_0[i] * fe_0 - ta2_xz_xzz_xxz_1[i] * fe_0 + ta2_xz_xzz_xxyz_0[i] * pa_y[i] - ta2_xz_xzz_xxyz_1[i] * pc_y[i];

        ta2_xz_xyzz_xxzz_0[i] = ta2_xz_xzz_xxzz_0[i] * pa_y[i] - ta2_xz_xzz_xxzz_1[i] * pc_y[i];

        ta2_xz_xyzz_xyyy_0[i] =
            3.0 * ta2_xz_xzz_xyy_0[i] * fe_0 - 3.0 * ta2_xz_xzz_xyy_1[i] * fe_0 + ta2_xz_xzz_xyyy_0[i] * pa_y[i] - ta2_xz_xzz_xyyy_1[i] * pc_y[i];

        ta2_xz_xyzz_xyyz_0[i] =
            2.0 * ta2_xz_xzz_xyz_0[i] * fe_0 - 2.0 * ta2_xz_xzz_xyz_1[i] * fe_0 + ta2_xz_xzz_xyyz_0[i] * pa_y[i] - ta2_xz_xzz_xyyz_1[i] * pc_y[i];

        ta2_xz_xyzz_xyzz_0[i] =
            ta2_xz_xzz_xzz_0[i] * fe_0 - ta2_xz_xzz_xzz_1[i] * fe_0 + ta2_xz_xzz_xyzz_0[i] * pa_y[i] - ta2_xz_xzz_xyzz_1[i] * pc_y[i];

        ta2_xz_xyzz_xzzz_0[i] = ta2_xz_xzz_xzzz_0[i] * pa_y[i] - ta2_xz_xzz_xzzz_1[i] * pc_y[i];

        ta2_xz_xyzz_yyyy_0[i] = ta1_z_yzz_yyyy_1[i] + ta2_xz_yzz_yyyy_0[i] * pa_x[i] - ta2_xz_yzz_yyyy_1[i] * pc_x[i];

        ta2_xz_xyzz_yyyz_0[i] = ta1_z_yzz_yyyz_1[i] + ta2_xz_yzz_yyyz_0[i] * pa_x[i] - ta2_xz_yzz_yyyz_1[i] * pc_x[i];

        ta2_xz_xyzz_yyzz_0[i] = ta1_z_yzz_yyzz_1[i] + ta2_xz_yzz_yyzz_0[i] * pa_x[i] - ta2_xz_yzz_yyzz_1[i] * pc_x[i];

        ta2_xz_xyzz_yzzz_0[i] = ta1_z_yzz_yzzz_1[i] + ta2_xz_yzz_yzzz_0[i] * pa_x[i] - ta2_xz_yzz_yzzz_1[i] * pc_x[i];

        ta2_xz_xyzz_zzzz_0[i] = ta2_xz_xzz_zzzz_0[i] * pa_y[i] - ta2_xz_xzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 585-600 components of targeted buffer : GG

    auto ta2_xz_xzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 585);

    auto ta2_xz_xzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 586);

    auto ta2_xz_xzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 587);

    auto ta2_xz_xzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 588);

    auto ta2_xz_xzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 589);

    auto ta2_xz_xzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 590);

    auto ta2_xz_xzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 591);

    auto ta2_xz_xzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 592);

    auto ta2_xz_xzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 593);

    auto ta2_xz_xzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 594);

    auto ta2_xz_xzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 595);

    auto ta2_xz_xzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 596);

    auto ta2_xz_xzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 597);

    auto ta2_xz_xzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 598);

    auto ta2_xz_xzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 599);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_zzz_xxxx_1,   \
                             ta1_z_zzz_xxxy_1,   \
                             ta1_z_zzz_xxxz_1,   \
                             ta1_z_zzz_xxyy_1,   \
                             ta1_z_zzz_xxyz_1,   \
                             ta1_z_zzz_xxzz_1,   \
                             ta1_z_zzz_xyyy_1,   \
                             ta1_z_zzz_xyyz_1,   \
                             ta1_z_zzz_xyzz_1,   \
                             ta1_z_zzz_xzzz_1,   \
                             ta1_z_zzz_yyyy_1,   \
                             ta1_z_zzz_yyyz_1,   \
                             ta1_z_zzz_yyzz_1,   \
                             ta1_z_zzz_yzzz_1,   \
                             ta1_z_zzz_zzzz_1,   \
                             ta2_xz_xzzz_xxxx_0, \
                             ta2_xz_xzzz_xxxy_0, \
                             ta2_xz_xzzz_xxxz_0, \
                             ta2_xz_xzzz_xxyy_0, \
                             ta2_xz_xzzz_xxyz_0, \
                             ta2_xz_xzzz_xxzz_0, \
                             ta2_xz_xzzz_xyyy_0, \
                             ta2_xz_xzzz_xyyz_0, \
                             ta2_xz_xzzz_xyzz_0, \
                             ta2_xz_xzzz_xzzz_0, \
                             ta2_xz_xzzz_yyyy_0, \
                             ta2_xz_xzzz_yyyz_0, \
                             ta2_xz_xzzz_yyzz_0, \
                             ta2_xz_xzzz_yzzz_0, \
                             ta2_xz_xzzz_zzzz_0, \
                             ta2_xz_zzz_xxx_0,   \
                             ta2_xz_zzz_xxx_1,   \
                             ta2_xz_zzz_xxxx_0,  \
                             ta2_xz_zzz_xxxx_1,  \
                             ta2_xz_zzz_xxxy_0,  \
                             ta2_xz_zzz_xxxy_1,  \
                             ta2_xz_zzz_xxxz_0,  \
                             ta2_xz_zzz_xxxz_1,  \
                             ta2_xz_zzz_xxy_0,   \
                             ta2_xz_zzz_xxy_1,   \
                             ta2_xz_zzz_xxyy_0,  \
                             ta2_xz_zzz_xxyy_1,  \
                             ta2_xz_zzz_xxyz_0,  \
                             ta2_xz_zzz_xxyz_1,  \
                             ta2_xz_zzz_xxz_0,   \
                             ta2_xz_zzz_xxz_1,   \
                             ta2_xz_zzz_xxzz_0,  \
                             ta2_xz_zzz_xxzz_1,  \
                             ta2_xz_zzz_xyy_0,   \
                             ta2_xz_zzz_xyy_1,   \
                             ta2_xz_zzz_xyyy_0,  \
                             ta2_xz_zzz_xyyy_1,  \
                             ta2_xz_zzz_xyyz_0,  \
                             ta2_xz_zzz_xyyz_1,  \
                             ta2_xz_zzz_xyz_0,   \
                             ta2_xz_zzz_xyz_1,   \
                             ta2_xz_zzz_xyzz_0,  \
                             ta2_xz_zzz_xyzz_1,  \
                             ta2_xz_zzz_xzz_0,   \
                             ta2_xz_zzz_xzz_1,   \
                             ta2_xz_zzz_xzzz_0,  \
                             ta2_xz_zzz_xzzz_1,  \
                             ta2_xz_zzz_yyy_0,   \
                             ta2_xz_zzz_yyy_1,   \
                             ta2_xz_zzz_yyyy_0,  \
                             ta2_xz_zzz_yyyy_1,  \
                             ta2_xz_zzz_yyyz_0,  \
                             ta2_xz_zzz_yyyz_1,  \
                             ta2_xz_zzz_yyz_0,   \
                             ta2_xz_zzz_yyz_1,   \
                             ta2_xz_zzz_yyzz_0,  \
                             ta2_xz_zzz_yyzz_1,  \
                             ta2_xz_zzz_yzz_0,   \
                             ta2_xz_zzz_yzz_1,   \
                             ta2_xz_zzz_yzzz_0,  \
                             ta2_xz_zzz_yzzz_1,  \
                             ta2_xz_zzz_zzz_0,   \
                             ta2_xz_zzz_zzz_1,   \
                             ta2_xz_zzz_zzzz_0,  \
                             ta2_xz_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_xzzz_xxxx_0[i] = 4.0 * ta2_xz_zzz_xxx_0[i] * fe_0 - 4.0 * ta2_xz_zzz_xxx_1[i] * fe_0 + ta1_z_zzz_xxxx_1[i] +
                                ta2_xz_zzz_xxxx_0[i] * pa_x[i] - ta2_xz_zzz_xxxx_1[i] * pc_x[i];

        ta2_xz_xzzz_xxxy_0[i] = 3.0 * ta2_xz_zzz_xxy_0[i] * fe_0 - 3.0 * ta2_xz_zzz_xxy_1[i] * fe_0 + ta1_z_zzz_xxxy_1[i] +
                                ta2_xz_zzz_xxxy_0[i] * pa_x[i] - ta2_xz_zzz_xxxy_1[i] * pc_x[i];

        ta2_xz_xzzz_xxxz_0[i] = 3.0 * ta2_xz_zzz_xxz_0[i] * fe_0 - 3.0 * ta2_xz_zzz_xxz_1[i] * fe_0 + ta1_z_zzz_xxxz_1[i] +
                                ta2_xz_zzz_xxxz_0[i] * pa_x[i] - ta2_xz_zzz_xxxz_1[i] * pc_x[i];

        ta2_xz_xzzz_xxyy_0[i] = 2.0 * ta2_xz_zzz_xyy_0[i] * fe_0 - 2.0 * ta2_xz_zzz_xyy_1[i] * fe_0 + ta1_z_zzz_xxyy_1[i] +
                                ta2_xz_zzz_xxyy_0[i] * pa_x[i] - ta2_xz_zzz_xxyy_1[i] * pc_x[i];

        ta2_xz_xzzz_xxyz_0[i] = 2.0 * ta2_xz_zzz_xyz_0[i] * fe_0 - 2.0 * ta2_xz_zzz_xyz_1[i] * fe_0 + ta1_z_zzz_xxyz_1[i] +
                                ta2_xz_zzz_xxyz_0[i] * pa_x[i] - ta2_xz_zzz_xxyz_1[i] * pc_x[i];

        ta2_xz_xzzz_xxzz_0[i] = 2.0 * ta2_xz_zzz_xzz_0[i] * fe_0 - 2.0 * ta2_xz_zzz_xzz_1[i] * fe_0 + ta1_z_zzz_xxzz_1[i] +
                                ta2_xz_zzz_xxzz_0[i] * pa_x[i] - ta2_xz_zzz_xxzz_1[i] * pc_x[i];

        ta2_xz_xzzz_xyyy_0[i] = ta2_xz_zzz_yyy_0[i] * fe_0 - ta2_xz_zzz_yyy_1[i] * fe_0 + ta1_z_zzz_xyyy_1[i] + ta2_xz_zzz_xyyy_0[i] * pa_x[i] -
                                ta2_xz_zzz_xyyy_1[i] * pc_x[i];

        ta2_xz_xzzz_xyyz_0[i] = ta2_xz_zzz_yyz_0[i] * fe_0 - ta2_xz_zzz_yyz_1[i] * fe_0 + ta1_z_zzz_xyyz_1[i] + ta2_xz_zzz_xyyz_0[i] * pa_x[i] -
                                ta2_xz_zzz_xyyz_1[i] * pc_x[i];

        ta2_xz_xzzz_xyzz_0[i] = ta2_xz_zzz_yzz_0[i] * fe_0 - ta2_xz_zzz_yzz_1[i] * fe_0 + ta1_z_zzz_xyzz_1[i] + ta2_xz_zzz_xyzz_0[i] * pa_x[i] -
                                ta2_xz_zzz_xyzz_1[i] * pc_x[i];

        ta2_xz_xzzz_xzzz_0[i] = ta2_xz_zzz_zzz_0[i] * fe_0 - ta2_xz_zzz_zzz_1[i] * fe_0 + ta1_z_zzz_xzzz_1[i] + ta2_xz_zzz_xzzz_0[i] * pa_x[i] -
                                ta2_xz_zzz_xzzz_1[i] * pc_x[i];

        ta2_xz_xzzz_yyyy_0[i] = ta1_z_zzz_yyyy_1[i] + ta2_xz_zzz_yyyy_0[i] * pa_x[i] - ta2_xz_zzz_yyyy_1[i] * pc_x[i];

        ta2_xz_xzzz_yyyz_0[i] = ta1_z_zzz_yyyz_1[i] + ta2_xz_zzz_yyyz_0[i] * pa_x[i] - ta2_xz_zzz_yyyz_1[i] * pc_x[i];

        ta2_xz_xzzz_yyzz_0[i] = ta1_z_zzz_yyzz_1[i] + ta2_xz_zzz_yyzz_0[i] * pa_x[i] - ta2_xz_zzz_yyzz_1[i] * pc_x[i];

        ta2_xz_xzzz_yzzz_0[i] = ta1_z_zzz_yzzz_1[i] + ta2_xz_zzz_yzzz_0[i] * pa_x[i] - ta2_xz_zzz_yzzz_1[i] * pc_x[i];

        ta2_xz_xzzz_zzzz_0[i] = ta1_z_zzz_zzzz_1[i] + ta2_xz_zzz_zzzz_0[i] * pa_x[i] - ta2_xz_zzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 600-615 components of targeted buffer : GG

    auto ta2_xz_yyyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 600);

    auto ta2_xz_yyyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 601);

    auto ta2_xz_yyyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 602);

    auto ta2_xz_yyyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 603);

    auto ta2_xz_yyyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 604);

    auto ta2_xz_yyyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 605);

    auto ta2_xz_yyyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 606);

    auto ta2_xz_yyyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 607);

    auto ta2_xz_yyyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 608);

    auto ta2_xz_yyyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 609);

    auto ta2_xz_yyyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 610);

    auto ta2_xz_yyyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 611);

    auto ta2_xz_yyyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 612);

    auto ta2_xz_yyyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 613);

    auto ta2_xz_yyyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 614);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta2_xz_yy_xxxx_0,   \
                             ta2_xz_yy_xxxx_1,   \
                             ta2_xz_yy_xxxy_0,   \
                             ta2_xz_yy_xxxy_1,   \
                             ta2_xz_yy_xxxz_0,   \
                             ta2_xz_yy_xxxz_1,   \
                             ta2_xz_yy_xxyy_0,   \
                             ta2_xz_yy_xxyy_1,   \
                             ta2_xz_yy_xxyz_0,   \
                             ta2_xz_yy_xxyz_1,   \
                             ta2_xz_yy_xxzz_0,   \
                             ta2_xz_yy_xxzz_1,   \
                             ta2_xz_yy_xyyy_0,   \
                             ta2_xz_yy_xyyy_1,   \
                             ta2_xz_yy_xyyz_0,   \
                             ta2_xz_yy_xyyz_1,   \
                             ta2_xz_yy_xyzz_0,   \
                             ta2_xz_yy_xyzz_1,   \
                             ta2_xz_yy_xzzz_0,   \
                             ta2_xz_yy_xzzz_1,   \
                             ta2_xz_yy_yyyy_0,   \
                             ta2_xz_yy_yyyy_1,   \
                             ta2_xz_yy_yyyz_0,   \
                             ta2_xz_yy_yyyz_1,   \
                             ta2_xz_yy_yyzz_0,   \
                             ta2_xz_yy_yyzz_1,   \
                             ta2_xz_yy_yzzz_0,   \
                             ta2_xz_yy_yzzz_1,   \
                             ta2_xz_yy_zzzz_0,   \
                             ta2_xz_yy_zzzz_1,   \
                             ta2_xz_yyy_xxx_0,   \
                             ta2_xz_yyy_xxx_1,   \
                             ta2_xz_yyy_xxxx_0,  \
                             ta2_xz_yyy_xxxx_1,  \
                             ta2_xz_yyy_xxxy_0,  \
                             ta2_xz_yyy_xxxy_1,  \
                             ta2_xz_yyy_xxxz_0,  \
                             ta2_xz_yyy_xxxz_1,  \
                             ta2_xz_yyy_xxy_0,   \
                             ta2_xz_yyy_xxy_1,   \
                             ta2_xz_yyy_xxyy_0,  \
                             ta2_xz_yyy_xxyy_1,  \
                             ta2_xz_yyy_xxyz_0,  \
                             ta2_xz_yyy_xxyz_1,  \
                             ta2_xz_yyy_xxz_0,   \
                             ta2_xz_yyy_xxz_1,   \
                             ta2_xz_yyy_xxzz_0,  \
                             ta2_xz_yyy_xxzz_1,  \
                             ta2_xz_yyy_xyy_0,   \
                             ta2_xz_yyy_xyy_1,   \
                             ta2_xz_yyy_xyyy_0,  \
                             ta2_xz_yyy_xyyy_1,  \
                             ta2_xz_yyy_xyyz_0,  \
                             ta2_xz_yyy_xyyz_1,  \
                             ta2_xz_yyy_xyz_0,   \
                             ta2_xz_yyy_xyz_1,   \
                             ta2_xz_yyy_xyzz_0,  \
                             ta2_xz_yyy_xyzz_1,  \
                             ta2_xz_yyy_xzz_0,   \
                             ta2_xz_yyy_xzz_1,   \
                             ta2_xz_yyy_xzzz_0,  \
                             ta2_xz_yyy_xzzz_1,  \
                             ta2_xz_yyy_yyy_0,   \
                             ta2_xz_yyy_yyy_1,   \
                             ta2_xz_yyy_yyyy_0,  \
                             ta2_xz_yyy_yyyy_1,  \
                             ta2_xz_yyy_yyyz_0,  \
                             ta2_xz_yyy_yyyz_1,  \
                             ta2_xz_yyy_yyz_0,   \
                             ta2_xz_yyy_yyz_1,   \
                             ta2_xz_yyy_yyzz_0,  \
                             ta2_xz_yyy_yyzz_1,  \
                             ta2_xz_yyy_yzz_0,   \
                             ta2_xz_yyy_yzz_1,   \
                             ta2_xz_yyy_yzzz_0,  \
                             ta2_xz_yyy_yzzz_1,  \
                             ta2_xz_yyy_zzz_0,   \
                             ta2_xz_yyy_zzz_1,   \
                             ta2_xz_yyy_zzzz_0,  \
                             ta2_xz_yyy_zzzz_1,  \
                             ta2_xz_yyyy_xxxx_0, \
                             ta2_xz_yyyy_xxxy_0, \
                             ta2_xz_yyyy_xxxz_0, \
                             ta2_xz_yyyy_xxyy_0, \
                             ta2_xz_yyyy_xxyz_0, \
                             ta2_xz_yyyy_xxzz_0, \
                             ta2_xz_yyyy_xyyy_0, \
                             ta2_xz_yyyy_xyyz_0, \
                             ta2_xz_yyyy_xyzz_0, \
                             ta2_xz_yyyy_xzzz_0, \
                             ta2_xz_yyyy_yyyy_0, \
                             ta2_xz_yyyy_yyyz_0, \
                             ta2_xz_yyyy_yyzz_0, \
                             ta2_xz_yyyy_yzzz_0, \
                             ta2_xz_yyyy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyyy_xxxx_0[i] =
            3.0 * ta2_xz_yy_xxxx_0[i] * fe_0 - 3.0 * ta2_xz_yy_xxxx_1[i] * fe_0 + ta2_xz_yyy_xxxx_0[i] * pa_y[i] - ta2_xz_yyy_xxxx_1[i] * pc_y[i];

        ta2_xz_yyyy_xxxy_0[i] = 3.0 * ta2_xz_yy_xxxy_0[i] * fe_0 - 3.0 * ta2_xz_yy_xxxy_1[i] * fe_0 + ta2_xz_yyy_xxx_0[i] * fe_0 -
                                ta2_xz_yyy_xxx_1[i] * fe_0 + ta2_xz_yyy_xxxy_0[i] * pa_y[i] - ta2_xz_yyy_xxxy_1[i] * pc_y[i];

        ta2_xz_yyyy_xxxz_0[i] =
            3.0 * ta2_xz_yy_xxxz_0[i] * fe_0 - 3.0 * ta2_xz_yy_xxxz_1[i] * fe_0 + ta2_xz_yyy_xxxz_0[i] * pa_y[i] - ta2_xz_yyy_xxxz_1[i] * pc_y[i];

        ta2_xz_yyyy_xxyy_0[i] = 3.0 * ta2_xz_yy_xxyy_0[i] * fe_0 - 3.0 * ta2_xz_yy_xxyy_1[i] * fe_0 + 2.0 * ta2_xz_yyy_xxy_0[i] * fe_0 -
                                2.0 * ta2_xz_yyy_xxy_1[i] * fe_0 + ta2_xz_yyy_xxyy_0[i] * pa_y[i] - ta2_xz_yyy_xxyy_1[i] * pc_y[i];

        ta2_xz_yyyy_xxyz_0[i] = 3.0 * ta2_xz_yy_xxyz_0[i] * fe_0 - 3.0 * ta2_xz_yy_xxyz_1[i] * fe_0 + ta2_xz_yyy_xxz_0[i] * fe_0 -
                                ta2_xz_yyy_xxz_1[i] * fe_0 + ta2_xz_yyy_xxyz_0[i] * pa_y[i] - ta2_xz_yyy_xxyz_1[i] * pc_y[i];

        ta2_xz_yyyy_xxzz_0[i] =
            3.0 * ta2_xz_yy_xxzz_0[i] * fe_0 - 3.0 * ta2_xz_yy_xxzz_1[i] * fe_0 + ta2_xz_yyy_xxzz_0[i] * pa_y[i] - ta2_xz_yyy_xxzz_1[i] * pc_y[i];

        ta2_xz_yyyy_xyyy_0[i] = 3.0 * ta2_xz_yy_xyyy_0[i] * fe_0 - 3.0 * ta2_xz_yy_xyyy_1[i] * fe_0 + 3.0 * ta2_xz_yyy_xyy_0[i] * fe_0 -
                                3.0 * ta2_xz_yyy_xyy_1[i] * fe_0 + ta2_xz_yyy_xyyy_0[i] * pa_y[i] - ta2_xz_yyy_xyyy_1[i] * pc_y[i];

        ta2_xz_yyyy_xyyz_0[i] = 3.0 * ta2_xz_yy_xyyz_0[i] * fe_0 - 3.0 * ta2_xz_yy_xyyz_1[i] * fe_0 + 2.0 * ta2_xz_yyy_xyz_0[i] * fe_0 -
                                2.0 * ta2_xz_yyy_xyz_1[i] * fe_0 + ta2_xz_yyy_xyyz_0[i] * pa_y[i] - ta2_xz_yyy_xyyz_1[i] * pc_y[i];

        ta2_xz_yyyy_xyzz_0[i] = 3.0 * ta2_xz_yy_xyzz_0[i] * fe_0 - 3.0 * ta2_xz_yy_xyzz_1[i] * fe_0 + ta2_xz_yyy_xzz_0[i] * fe_0 -
                                ta2_xz_yyy_xzz_1[i] * fe_0 + ta2_xz_yyy_xyzz_0[i] * pa_y[i] - ta2_xz_yyy_xyzz_1[i] * pc_y[i];

        ta2_xz_yyyy_xzzz_0[i] =
            3.0 * ta2_xz_yy_xzzz_0[i] * fe_0 - 3.0 * ta2_xz_yy_xzzz_1[i] * fe_0 + ta2_xz_yyy_xzzz_0[i] * pa_y[i] - ta2_xz_yyy_xzzz_1[i] * pc_y[i];

        ta2_xz_yyyy_yyyy_0[i] = 3.0 * ta2_xz_yy_yyyy_0[i] * fe_0 - 3.0 * ta2_xz_yy_yyyy_1[i] * fe_0 + 4.0 * ta2_xz_yyy_yyy_0[i] * fe_0 -
                                4.0 * ta2_xz_yyy_yyy_1[i] * fe_0 + ta2_xz_yyy_yyyy_0[i] * pa_y[i] - ta2_xz_yyy_yyyy_1[i] * pc_y[i];

        ta2_xz_yyyy_yyyz_0[i] = 3.0 * ta2_xz_yy_yyyz_0[i] * fe_0 - 3.0 * ta2_xz_yy_yyyz_1[i] * fe_0 + 3.0 * ta2_xz_yyy_yyz_0[i] * fe_0 -
                                3.0 * ta2_xz_yyy_yyz_1[i] * fe_0 + ta2_xz_yyy_yyyz_0[i] * pa_y[i] - ta2_xz_yyy_yyyz_1[i] * pc_y[i];

        ta2_xz_yyyy_yyzz_0[i] = 3.0 * ta2_xz_yy_yyzz_0[i] * fe_0 - 3.0 * ta2_xz_yy_yyzz_1[i] * fe_0 + 2.0 * ta2_xz_yyy_yzz_0[i] * fe_0 -
                                2.0 * ta2_xz_yyy_yzz_1[i] * fe_0 + ta2_xz_yyy_yyzz_0[i] * pa_y[i] - ta2_xz_yyy_yyzz_1[i] * pc_y[i];

        ta2_xz_yyyy_yzzz_0[i] = 3.0 * ta2_xz_yy_yzzz_0[i] * fe_0 - 3.0 * ta2_xz_yy_yzzz_1[i] * fe_0 + ta2_xz_yyy_zzz_0[i] * fe_0 -
                                ta2_xz_yyy_zzz_1[i] * fe_0 + ta2_xz_yyy_yzzz_0[i] * pa_y[i] - ta2_xz_yyy_yzzz_1[i] * pc_y[i];

        ta2_xz_yyyy_zzzz_0[i] =
            3.0 * ta2_xz_yy_zzzz_0[i] * fe_0 - 3.0 * ta2_xz_yy_zzzz_1[i] * fe_0 + ta2_xz_yyy_zzzz_0[i] * pa_y[i] - ta2_xz_yyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 615-630 components of targeted buffer : GG

    auto ta2_xz_yyyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 615);

    auto ta2_xz_yyyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 616);

    auto ta2_xz_yyyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 617);

    auto ta2_xz_yyyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 618);

    auto ta2_xz_yyyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 619);

    auto ta2_xz_yyyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 620);

    auto ta2_xz_yyyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 621);

    auto ta2_xz_yyyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 622);

    auto ta2_xz_yyyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 623);

    auto ta2_xz_yyyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 624);

    auto ta2_xz_yyyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 625);

    auto ta2_xz_yyyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 626);

    auto ta2_xz_yyyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 627);

    auto ta2_xz_yyyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 628);

    auto ta2_xz_yyyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 629);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_yyy_xxxx_1,   \
                             ta1_x_yyy_xxxy_1,   \
                             ta1_x_yyy_xxyy_1,   \
                             ta1_x_yyy_xxyz_1,   \
                             ta1_x_yyy_xyyy_1,   \
                             ta1_x_yyy_xyyz_1,   \
                             ta1_x_yyy_xyzz_1,   \
                             ta1_x_yyy_yyyy_1,   \
                             ta1_x_yyy_yyyz_1,   \
                             ta1_x_yyy_yyzz_1,   \
                             ta1_x_yyy_yzzz_1,   \
                             ta2_xz_yyy_xxxx_0,  \
                             ta2_xz_yyy_xxxx_1,  \
                             ta2_xz_yyy_xxxy_0,  \
                             ta2_xz_yyy_xxxy_1,  \
                             ta2_xz_yyy_xxy_0,   \
                             ta2_xz_yyy_xxy_1,   \
                             ta2_xz_yyy_xxyy_0,  \
                             ta2_xz_yyy_xxyy_1,  \
                             ta2_xz_yyy_xxyz_0,  \
                             ta2_xz_yyy_xxyz_1,  \
                             ta2_xz_yyy_xyy_0,   \
                             ta2_xz_yyy_xyy_1,   \
                             ta2_xz_yyy_xyyy_0,  \
                             ta2_xz_yyy_xyyy_1,  \
                             ta2_xz_yyy_xyyz_0,  \
                             ta2_xz_yyy_xyyz_1,  \
                             ta2_xz_yyy_xyz_0,   \
                             ta2_xz_yyy_xyz_1,   \
                             ta2_xz_yyy_xyzz_0,  \
                             ta2_xz_yyy_xyzz_1,  \
                             ta2_xz_yyy_yyy_0,   \
                             ta2_xz_yyy_yyy_1,   \
                             ta2_xz_yyy_yyyy_0,  \
                             ta2_xz_yyy_yyyy_1,  \
                             ta2_xz_yyy_yyyz_0,  \
                             ta2_xz_yyy_yyyz_1,  \
                             ta2_xz_yyy_yyz_0,   \
                             ta2_xz_yyy_yyz_1,   \
                             ta2_xz_yyy_yyzz_0,  \
                             ta2_xz_yyy_yyzz_1,  \
                             ta2_xz_yyy_yzz_0,   \
                             ta2_xz_yyy_yzz_1,   \
                             ta2_xz_yyy_yzzz_0,  \
                             ta2_xz_yyy_yzzz_1,  \
                             ta2_xz_yyyz_xxxx_0, \
                             ta2_xz_yyyz_xxxy_0, \
                             ta2_xz_yyyz_xxxz_0, \
                             ta2_xz_yyyz_xxyy_0, \
                             ta2_xz_yyyz_xxyz_0, \
                             ta2_xz_yyyz_xxzz_0, \
                             ta2_xz_yyyz_xyyy_0, \
                             ta2_xz_yyyz_xyyz_0, \
                             ta2_xz_yyyz_xyzz_0, \
                             ta2_xz_yyyz_xzzz_0, \
                             ta2_xz_yyyz_yyyy_0, \
                             ta2_xz_yyyz_yyyz_0, \
                             ta2_xz_yyyz_yyzz_0, \
                             ta2_xz_yyyz_yzzz_0, \
                             ta2_xz_yyyz_zzzz_0, \
                             ta2_xz_yyz_xxxz_0,  \
                             ta2_xz_yyz_xxxz_1,  \
                             ta2_xz_yyz_xxzz_0,  \
                             ta2_xz_yyz_xxzz_1,  \
                             ta2_xz_yyz_xzzz_0,  \
                             ta2_xz_yyz_xzzz_1,  \
                             ta2_xz_yyz_zzzz_0,  \
                             ta2_xz_yyz_zzzz_1,  \
                             ta2_xz_yz_xxxz_0,   \
                             ta2_xz_yz_xxxz_1,   \
                             ta2_xz_yz_xxzz_0,   \
                             ta2_xz_yz_xxzz_1,   \
                             ta2_xz_yz_xzzz_0,   \
                             ta2_xz_yz_xzzz_1,   \
                             ta2_xz_yz_zzzz_0,   \
                             ta2_xz_yz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyyz_xxxx_0[i] = ta1_x_yyy_xxxx_1[i] + ta2_xz_yyy_xxxx_0[i] * pa_z[i] - ta2_xz_yyy_xxxx_1[i] * pc_z[i];

        ta2_xz_yyyz_xxxy_0[i] = ta1_x_yyy_xxxy_1[i] + ta2_xz_yyy_xxxy_0[i] * pa_z[i] - ta2_xz_yyy_xxxy_1[i] * pc_z[i];

        ta2_xz_yyyz_xxxz_0[i] =
            2.0 * ta2_xz_yz_xxxz_0[i] * fe_0 - 2.0 * ta2_xz_yz_xxxz_1[i] * fe_0 + ta2_xz_yyz_xxxz_0[i] * pa_y[i] - ta2_xz_yyz_xxxz_1[i] * pc_y[i];

        ta2_xz_yyyz_xxyy_0[i] = ta1_x_yyy_xxyy_1[i] + ta2_xz_yyy_xxyy_0[i] * pa_z[i] - ta2_xz_yyy_xxyy_1[i] * pc_z[i];

        ta2_xz_yyyz_xxyz_0[i] = ta2_xz_yyy_xxy_0[i] * fe_0 - ta2_xz_yyy_xxy_1[i] * fe_0 + ta1_x_yyy_xxyz_1[i] + ta2_xz_yyy_xxyz_0[i] * pa_z[i] -
                                ta2_xz_yyy_xxyz_1[i] * pc_z[i];

        ta2_xz_yyyz_xxzz_0[i] =
            2.0 * ta2_xz_yz_xxzz_0[i] * fe_0 - 2.0 * ta2_xz_yz_xxzz_1[i] * fe_0 + ta2_xz_yyz_xxzz_0[i] * pa_y[i] - ta2_xz_yyz_xxzz_1[i] * pc_y[i];

        ta2_xz_yyyz_xyyy_0[i] = ta1_x_yyy_xyyy_1[i] + ta2_xz_yyy_xyyy_0[i] * pa_z[i] - ta2_xz_yyy_xyyy_1[i] * pc_z[i];

        ta2_xz_yyyz_xyyz_0[i] = ta2_xz_yyy_xyy_0[i] * fe_0 - ta2_xz_yyy_xyy_1[i] * fe_0 + ta1_x_yyy_xyyz_1[i] + ta2_xz_yyy_xyyz_0[i] * pa_z[i] -
                                ta2_xz_yyy_xyyz_1[i] * pc_z[i];

        ta2_xz_yyyz_xyzz_0[i] = 2.0 * ta2_xz_yyy_xyz_0[i] * fe_0 - 2.0 * ta2_xz_yyy_xyz_1[i] * fe_0 + ta1_x_yyy_xyzz_1[i] +
                                ta2_xz_yyy_xyzz_0[i] * pa_z[i] - ta2_xz_yyy_xyzz_1[i] * pc_z[i];

        ta2_xz_yyyz_xzzz_0[i] =
            2.0 * ta2_xz_yz_xzzz_0[i] * fe_0 - 2.0 * ta2_xz_yz_xzzz_1[i] * fe_0 + ta2_xz_yyz_xzzz_0[i] * pa_y[i] - ta2_xz_yyz_xzzz_1[i] * pc_y[i];

        ta2_xz_yyyz_yyyy_0[i] = ta1_x_yyy_yyyy_1[i] + ta2_xz_yyy_yyyy_0[i] * pa_z[i] - ta2_xz_yyy_yyyy_1[i] * pc_z[i];

        ta2_xz_yyyz_yyyz_0[i] = ta2_xz_yyy_yyy_0[i] * fe_0 - ta2_xz_yyy_yyy_1[i] * fe_0 + ta1_x_yyy_yyyz_1[i] + ta2_xz_yyy_yyyz_0[i] * pa_z[i] -
                                ta2_xz_yyy_yyyz_1[i] * pc_z[i];

        ta2_xz_yyyz_yyzz_0[i] = 2.0 * ta2_xz_yyy_yyz_0[i] * fe_0 - 2.0 * ta2_xz_yyy_yyz_1[i] * fe_0 + ta1_x_yyy_yyzz_1[i] +
                                ta2_xz_yyy_yyzz_0[i] * pa_z[i] - ta2_xz_yyy_yyzz_1[i] * pc_z[i];

        ta2_xz_yyyz_yzzz_0[i] = 3.0 * ta2_xz_yyy_yzz_0[i] * fe_0 - 3.0 * ta2_xz_yyy_yzz_1[i] * fe_0 + ta1_x_yyy_yzzz_1[i] +
                                ta2_xz_yyy_yzzz_0[i] * pa_z[i] - ta2_xz_yyy_yzzz_1[i] * pc_z[i];

        ta2_xz_yyyz_zzzz_0[i] =
            2.0 * ta2_xz_yz_zzzz_0[i] * fe_0 - 2.0 * ta2_xz_yz_zzzz_1[i] * fe_0 + ta2_xz_yyz_zzzz_0[i] * pa_y[i] - ta2_xz_yyz_zzzz_1[i] * pc_y[i];
    }

    // Set up 630-645 components of targeted buffer : GG

    auto ta2_xz_yyzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 630);

    auto ta2_xz_yyzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 631);

    auto ta2_xz_yyzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 632);

    auto ta2_xz_yyzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 633);

    auto ta2_xz_yyzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 634);

    auto ta2_xz_yyzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 635);

    auto ta2_xz_yyzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 636);

    auto ta2_xz_yyzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 637);

    auto ta2_xz_yyzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 638);

    auto ta2_xz_yyzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 639);

    auto ta2_xz_yyzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 640);

    auto ta2_xz_yyzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 641);

    auto ta2_xz_yyzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 642);

    auto ta2_xz_yyzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 643);

    auto ta2_xz_yyzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 644);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_yyz_xxxy_1,   \
                             ta1_x_yyz_xxyy_1,   \
                             ta1_x_yyz_xyyy_1,   \
                             ta1_x_yyz_yyyy_1,   \
                             ta2_xz_yy_xxxy_0,   \
                             ta2_xz_yy_xxxy_1,   \
                             ta2_xz_yy_xxyy_0,   \
                             ta2_xz_yy_xxyy_1,   \
                             ta2_xz_yy_xyyy_0,   \
                             ta2_xz_yy_xyyy_1,   \
                             ta2_xz_yy_yyyy_0,   \
                             ta2_xz_yy_yyyy_1,   \
                             ta2_xz_yyz_xxxy_0,  \
                             ta2_xz_yyz_xxxy_1,  \
                             ta2_xz_yyz_xxyy_0,  \
                             ta2_xz_yyz_xxyy_1,  \
                             ta2_xz_yyz_xyyy_0,  \
                             ta2_xz_yyz_xyyy_1,  \
                             ta2_xz_yyz_yyyy_0,  \
                             ta2_xz_yyz_yyyy_1,  \
                             ta2_xz_yyzz_xxxx_0, \
                             ta2_xz_yyzz_xxxy_0, \
                             ta2_xz_yyzz_xxxz_0, \
                             ta2_xz_yyzz_xxyy_0, \
                             ta2_xz_yyzz_xxyz_0, \
                             ta2_xz_yyzz_xxzz_0, \
                             ta2_xz_yyzz_xyyy_0, \
                             ta2_xz_yyzz_xyyz_0, \
                             ta2_xz_yyzz_xyzz_0, \
                             ta2_xz_yyzz_xzzz_0, \
                             ta2_xz_yyzz_yyyy_0, \
                             ta2_xz_yyzz_yyyz_0, \
                             ta2_xz_yyzz_yyzz_0, \
                             ta2_xz_yyzz_yzzz_0, \
                             ta2_xz_yyzz_zzzz_0, \
                             ta2_xz_yzz_xxxx_0,  \
                             ta2_xz_yzz_xxxx_1,  \
                             ta2_xz_yzz_xxxz_0,  \
                             ta2_xz_yzz_xxxz_1,  \
                             ta2_xz_yzz_xxyz_0,  \
                             ta2_xz_yzz_xxyz_1,  \
                             ta2_xz_yzz_xxz_0,   \
                             ta2_xz_yzz_xxz_1,   \
                             ta2_xz_yzz_xxzz_0,  \
                             ta2_xz_yzz_xxzz_1,  \
                             ta2_xz_yzz_xyyz_0,  \
                             ta2_xz_yzz_xyyz_1,  \
                             ta2_xz_yzz_xyz_0,   \
                             ta2_xz_yzz_xyz_1,   \
                             ta2_xz_yzz_xyzz_0,  \
                             ta2_xz_yzz_xyzz_1,  \
                             ta2_xz_yzz_xzz_0,   \
                             ta2_xz_yzz_xzz_1,   \
                             ta2_xz_yzz_xzzz_0,  \
                             ta2_xz_yzz_xzzz_1,  \
                             ta2_xz_yzz_yyyz_0,  \
                             ta2_xz_yzz_yyyz_1,  \
                             ta2_xz_yzz_yyz_0,   \
                             ta2_xz_yzz_yyz_1,   \
                             ta2_xz_yzz_yyzz_0,  \
                             ta2_xz_yzz_yyzz_1,  \
                             ta2_xz_yzz_yzz_0,   \
                             ta2_xz_yzz_yzz_1,   \
                             ta2_xz_yzz_yzzz_0,  \
                             ta2_xz_yzz_yzzz_1,  \
                             ta2_xz_yzz_zzz_0,   \
                             ta2_xz_yzz_zzz_1,   \
                             ta2_xz_yzz_zzzz_0,  \
                             ta2_xz_yzz_zzzz_1,  \
                             ta2_xz_zz_xxxx_0,   \
                             ta2_xz_zz_xxxx_1,   \
                             ta2_xz_zz_xxxz_0,   \
                             ta2_xz_zz_xxxz_1,   \
                             ta2_xz_zz_xxyz_0,   \
                             ta2_xz_zz_xxyz_1,   \
                             ta2_xz_zz_xxzz_0,   \
                             ta2_xz_zz_xxzz_1,   \
                             ta2_xz_zz_xyyz_0,   \
                             ta2_xz_zz_xyyz_1,   \
                             ta2_xz_zz_xyzz_0,   \
                             ta2_xz_zz_xyzz_1,   \
                             ta2_xz_zz_xzzz_0,   \
                             ta2_xz_zz_xzzz_1,   \
                             ta2_xz_zz_yyyz_0,   \
                             ta2_xz_zz_yyyz_1,   \
                             ta2_xz_zz_yyzz_0,   \
                             ta2_xz_zz_yyzz_1,   \
                             ta2_xz_zz_yzzz_0,   \
                             ta2_xz_zz_yzzz_1,   \
                             ta2_xz_zz_zzzz_0,   \
                             ta2_xz_zz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yyzz_xxxx_0[i] =
            ta2_xz_zz_xxxx_0[i] * fe_0 - ta2_xz_zz_xxxx_1[i] * fe_0 + ta2_xz_yzz_xxxx_0[i] * pa_y[i] - ta2_xz_yzz_xxxx_1[i] * pc_y[i];

        ta2_xz_yyzz_xxxy_0[i] = ta2_xz_yy_xxxy_0[i] * fe_0 - ta2_xz_yy_xxxy_1[i] * fe_0 + ta1_x_yyz_xxxy_1[i] + ta2_xz_yyz_xxxy_0[i] * pa_z[i] -
                                ta2_xz_yyz_xxxy_1[i] * pc_z[i];

        ta2_xz_yyzz_xxxz_0[i] =
            ta2_xz_zz_xxxz_0[i] * fe_0 - ta2_xz_zz_xxxz_1[i] * fe_0 + ta2_xz_yzz_xxxz_0[i] * pa_y[i] - ta2_xz_yzz_xxxz_1[i] * pc_y[i];

        ta2_xz_yyzz_xxyy_0[i] = ta2_xz_yy_xxyy_0[i] * fe_0 - ta2_xz_yy_xxyy_1[i] * fe_0 + ta1_x_yyz_xxyy_1[i] + ta2_xz_yyz_xxyy_0[i] * pa_z[i] -
                                ta2_xz_yyz_xxyy_1[i] * pc_z[i];

        ta2_xz_yyzz_xxyz_0[i] = ta2_xz_zz_xxyz_0[i] * fe_0 - ta2_xz_zz_xxyz_1[i] * fe_0 + ta2_xz_yzz_xxz_0[i] * fe_0 - ta2_xz_yzz_xxz_1[i] * fe_0 +
                                ta2_xz_yzz_xxyz_0[i] * pa_y[i] - ta2_xz_yzz_xxyz_1[i] * pc_y[i];

        ta2_xz_yyzz_xxzz_0[i] =
            ta2_xz_zz_xxzz_0[i] * fe_0 - ta2_xz_zz_xxzz_1[i] * fe_0 + ta2_xz_yzz_xxzz_0[i] * pa_y[i] - ta2_xz_yzz_xxzz_1[i] * pc_y[i];

        ta2_xz_yyzz_xyyy_0[i] = ta2_xz_yy_xyyy_0[i] * fe_0 - ta2_xz_yy_xyyy_1[i] * fe_0 + ta1_x_yyz_xyyy_1[i] + ta2_xz_yyz_xyyy_0[i] * pa_z[i] -
                                ta2_xz_yyz_xyyy_1[i] * pc_z[i];

        ta2_xz_yyzz_xyyz_0[i] = ta2_xz_zz_xyyz_0[i] * fe_0 - ta2_xz_zz_xyyz_1[i] * fe_0 + 2.0 * ta2_xz_yzz_xyz_0[i] * fe_0 -
                                2.0 * ta2_xz_yzz_xyz_1[i] * fe_0 + ta2_xz_yzz_xyyz_0[i] * pa_y[i] - ta2_xz_yzz_xyyz_1[i] * pc_y[i];

        ta2_xz_yyzz_xyzz_0[i] = ta2_xz_zz_xyzz_0[i] * fe_0 - ta2_xz_zz_xyzz_1[i] * fe_0 + ta2_xz_yzz_xzz_0[i] * fe_0 - ta2_xz_yzz_xzz_1[i] * fe_0 +
                                ta2_xz_yzz_xyzz_0[i] * pa_y[i] - ta2_xz_yzz_xyzz_1[i] * pc_y[i];

        ta2_xz_yyzz_xzzz_0[i] =
            ta2_xz_zz_xzzz_0[i] * fe_0 - ta2_xz_zz_xzzz_1[i] * fe_0 + ta2_xz_yzz_xzzz_0[i] * pa_y[i] - ta2_xz_yzz_xzzz_1[i] * pc_y[i];

        ta2_xz_yyzz_yyyy_0[i] = ta2_xz_yy_yyyy_0[i] * fe_0 - ta2_xz_yy_yyyy_1[i] * fe_0 + ta1_x_yyz_yyyy_1[i] + ta2_xz_yyz_yyyy_0[i] * pa_z[i] -
                                ta2_xz_yyz_yyyy_1[i] * pc_z[i];

        ta2_xz_yyzz_yyyz_0[i] = ta2_xz_zz_yyyz_0[i] * fe_0 - ta2_xz_zz_yyyz_1[i] * fe_0 + 3.0 * ta2_xz_yzz_yyz_0[i] * fe_0 -
                                3.0 * ta2_xz_yzz_yyz_1[i] * fe_0 + ta2_xz_yzz_yyyz_0[i] * pa_y[i] - ta2_xz_yzz_yyyz_1[i] * pc_y[i];

        ta2_xz_yyzz_yyzz_0[i] = ta2_xz_zz_yyzz_0[i] * fe_0 - ta2_xz_zz_yyzz_1[i] * fe_0 + 2.0 * ta2_xz_yzz_yzz_0[i] * fe_0 -
                                2.0 * ta2_xz_yzz_yzz_1[i] * fe_0 + ta2_xz_yzz_yyzz_0[i] * pa_y[i] - ta2_xz_yzz_yyzz_1[i] * pc_y[i];

        ta2_xz_yyzz_yzzz_0[i] = ta2_xz_zz_yzzz_0[i] * fe_0 - ta2_xz_zz_yzzz_1[i] * fe_0 + ta2_xz_yzz_zzz_0[i] * fe_0 - ta2_xz_yzz_zzz_1[i] * fe_0 +
                                ta2_xz_yzz_yzzz_0[i] * pa_y[i] - ta2_xz_yzz_yzzz_1[i] * pc_y[i];

        ta2_xz_yyzz_zzzz_0[i] =
            ta2_xz_zz_zzzz_0[i] * fe_0 - ta2_xz_zz_zzzz_1[i] * fe_0 + ta2_xz_yzz_zzzz_0[i] * pa_y[i] - ta2_xz_yzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 645-660 components of targeted buffer : GG

    auto ta2_xz_yzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 645);

    auto ta2_xz_yzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 646);

    auto ta2_xz_yzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 647);

    auto ta2_xz_yzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 648);

    auto ta2_xz_yzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 649);

    auto ta2_xz_yzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 650);

    auto ta2_xz_yzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 651);

    auto ta2_xz_yzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 652);

    auto ta2_xz_yzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 653);

    auto ta2_xz_yzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 654);

    auto ta2_xz_yzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 655);

    auto ta2_xz_yzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 656);

    auto ta2_xz_yzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 657);

    auto ta2_xz_yzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 658);

    auto ta2_xz_yzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 659);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta2_xz_yzzz_xxxx_0, \
                             ta2_xz_yzzz_xxxy_0, \
                             ta2_xz_yzzz_xxxz_0, \
                             ta2_xz_yzzz_xxyy_0, \
                             ta2_xz_yzzz_xxyz_0, \
                             ta2_xz_yzzz_xxzz_0, \
                             ta2_xz_yzzz_xyyy_0, \
                             ta2_xz_yzzz_xyyz_0, \
                             ta2_xz_yzzz_xyzz_0, \
                             ta2_xz_yzzz_xzzz_0, \
                             ta2_xz_yzzz_yyyy_0, \
                             ta2_xz_yzzz_yyyz_0, \
                             ta2_xz_yzzz_yyzz_0, \
                             ta2_xz_yzzz_yzzz_0, \
                             ta2_xz_yzzz_zzzz_0, \
                             ta2_xz_zzz_xxx_0,   \
                             ta2_xz_zzz_xxx_1,   \
                             ta2_xz_zzz_xxxx_0,  \
                             ta2_xz_zzz_xxxx_1,  \
                             ta2_xz_zzz_xxxy_0,  \
                             ta2_xz_zzz_xxxy_1,  \
                             ta2_xz_zzz_xxxz_0,  \
                             ta2_xz_zzz_xxxz_1,  \
                             ta2_xz_zzz_xxy_0,   \
                             ta2_xz_zzz_xxy_1,   \
                             ta2_xz_zzz_xxyy_0,  \
                             ta2_xz_zzz_xxyy_1,  \
                             ta2_xz_zzz_xxyz_0,  \
                             ta2_xz_zzz_xxyz_1,  \
                             ta2_xz_zzz_xxz_0,   \
                             ta2_xz_zzz_xxz_1,   \
                             ta2_xz_zzz_xxzz_0,  \
                             ta2_xz_zzz_xxzz_1,  \
                             ta2_xz_zzz_xyy_0,   \
                             ta2_xz_zzz_xyy_1,   \
                             ta2_xz_zzz_xyyy_0,  \
                             ta2_xz_zzz_xyyy_1,  \
                             ta2_xz_zzz_xyyz_0,  \
                             ta2_xz_zzz_xyyz_1,  \
                             ta2_xz_zzz_xyz_0,   \
                             ta2_xz_zzz_xyz_1,   \
                             ta2_xz_zzz_xyzz_0,  \
                             ta2_xz_zzz_xyzz_1,  \
                             ta2_xz_zzz_xzz_0,   \
                             ta2_xz_zzz_xzz_1,   \
                             ta2_xz_zzz_xzzz_0,  \
                             ta2_xz_zzz_xzzz_1,  \
                             ta2_xz_zzz_yyy_0,   \
                             ta2_xz_zzz_yyy_1,   \
                             ta2_xz_zzz_yyyy_0,  \
                             ta2_xz_zzz_yyyy_1,  \
                             ta2_xz_zzz_yyyz_0,  \
                             ta2_xz_zzz_yyyz_1,  \
                             ta2_xz_zzz_yyz_0,   \
                             ta2_xz_zzz_yyz_1,   \
                             ta2_xz_zzz_yyzz_0,  \
                             ta2_xz_zzz_yyzz_1,  \
                             ta2_xz_zzz_yzz_0,   \
                             ta2_xz_zzz_yzz_1,   \
                             ta2_xz_zzz_yzzz_0,  \
                             ta2_xz_zzz_yzzz_1,  \
                             ta2_xz_zzz_zzz_0,   \
                             ta2_xz_zzz_zzz_1,   \
                             ta2_xz_zzz_zzzz_0,  \
                             ta2_xz_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_yzzz_xxxx_0[i] = ta2_xz_zzz_xxxx_0[i] * pa_y[i] - ta2_xz_zzz_xxxx_1[i] * pc_y[i];

        ta2_xz_yzzz_xxxy_0[i] =
            ta2_xz_zzz_xxx_0[i] * fe_0 - ta2_xz_zzz_xxx_1[i] * fe_0 + ta2_xz_zzz_xxxy_0[i] * pa_y[i] - ta2_xz_zzz_xxxy_1[i] * pc_y[i];

        ta2_xz_yzzz_xxxz_0[i] = ta2_xz_zzz_xxxz_0[i] * pa_y[i] - ta2_xz_zzz_xxxz_1[i] * pc_y[i];

        ta2_xz_yzzz_xxyy_0[i] =
            2.0 * ta2_xz_zzz_xxy_0[i] * fe_0 - 2.0 * ta2_xz_zzz_xxy_1[i] * fe_0 + ta2_xz_zzz_xxyy_0[i] * pa_y[i] - ta2_xz_zzz_xxyy_1[i] * pc_y[i];

        ta2_xz_yzzz_xxyz_0[i] =
            ta2_xz_zzz_xxz_0[i] * fe_0 - ta2_xz_zzz_xxz_1[i] * fe_0 + ta2_xz_zzz_xxyz_0[i] * pa_y[i] - ta2_xz_zzz_xxyz_1[i] * pc_y[i];

        ta2_xz_yzzz_xxzz_0[i] = ta2_xz_zzz_xxzz_0[i] * pa_y[i] - ta2_xz_zzz_xxzz_1[i] * pc_y[i];

        ta2_xz_yzzz_xyyy_0[i] =
            3.0 * ta2_xz_zzz_xyy_0[i] * fe_0 - 3.0 * ta2_xz_zzz_xyy_1[i] * fe_0 + ta2_xz_zzz_xyyy_0[i] * pa_y[i] - ta2_xz_zzz_xyyy_1[i] * pc_y[i];

        ta2_xz_yzzz_xyyz_0[i] =
            2.0 * ta2_xz_zzz_xyz_0[i] * fe_0 - 2.0 * ta2_xz_zzz_xyz_1[i] * fe_0 + ta2_xz_zzz_xyyz_0[i] * pa_y[i] - ta2_xz_zzz_xyyz_1[i] * pc_y[i];

        ta2_xz_yzzz_xyzz_0[i] =
            ta2_xz_zzz_xzz_0[i] * fe_0 - ta2_xz_zzz_xzz_1[i] * fe_0 + ta2_xz_zzz_xyzz_0[i] * pa_y[i] - ta2_xz_zzz_xyzz_1[i] * pc_y[i];

        ta2_xz_yzzz_xzzz_0[i] = ta2_xz_zzz_xzzz_0[i] * pa_y[i] - ta2_xz_zzz_xzzz_1[i] * pc_y[i];

        ta2_xz_yzzz_yyyy_0[i] =
            4.0 * ta2_xz_zzz_yyy_0[i] * fe_0 - 4.0 * ta2_xz_zzz_yyy_1[i] * fe_0 + ta2_xz_zzz_yyyy_0[i] * pa_y[i] - ta2_xz_zzz_yyyy_1[i] * pc_y[i];

        ta2_xz_yzzz_yyyz_0[i] =
            3.0 * ta2_xz_zzz_yyz_0[i] * fe_0 - 3.0 * ta2_xz_zzz_yyz_1[i] * fe_0 + ta2_xz_zzz_yyyz_0[i] * pa_y[i] - ta2_xz_zzz_yyyz_1[i] * pc_y[i];

        ta2_xz_yzzz_yyzz_0[i] =
            2.0 * ta2_xz_zzz_yzz_0[i] * fe_0 - 2.0 * ta2_xz_zzz_yzz_1[i] * fe_0 + ta2_xz_zzz_yyzz_0[i] * pa_y[i] - ta2_xz_zzz_yyzz_1[i] * pc_y[i];

        ta2_xz_yzzz_yzzz_0[i] =
            ta2_xz_zzz_zzz_0[i] * fe_0 - ta2_xz_zzz_zzz_1[i] * fe_0 + ta2_xz_zzz_yzzz_0[i] * pa_y[i] - ta2_xz_zzz_yzzz_1[i] * pc_y[i];

        ta2_xz_yzzz_zzzz_0[i] = ta2_xz_zzz_zzzz_0[i] * pa_y[i] - ta2_xz_zzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 660-675 components of targeted buffer : GG

    auto ta2_xz_zzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 660);

    auto ta2_xz_zzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 661);

    auto ta2_xz_zzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 662);

    auto ta2_xz_zzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 663);

    auto ta2_xz_zzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 664);

    auto ta2_xz_zzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 665);

    auto ta2_xz_zzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 666);

    auto ta2_xz_zzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 667);

    auto ta2_xz_zzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 668);

    auto ta2_xz_zzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 669);

    auto ta2_xz_zzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 670);

    auto ta2_xz_zzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 671);

    auto ta2_xz_zzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 672);

    auto ta2_xz_zzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 673);

    auto ta2_xz_zzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 674);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_x_zzz_xxxx_1,   \
                             ta1_x_zzz_xxxy_1,   \
                             ta1_x_zzz_xxxz_1,   \
                             ta1_x_zzz_xxyy_1,   \
                             ta1_x_zzz_xxyz_1,   \
                             ta1_x_zzz_xxzz_1,   \
                             ta1_x_zzz_xyyy_1,   \
                             ta1_x_zzz_xyyz_1,   \
                             ta1_x_zzz_xyzz_1,   \
                             ta1_x_zzz_xzzz_1,   \
                             ta1_x_zzz_yyyy_1,   \
                             ta1_x_zzz_yyyz_1,   \
                             ta1_x_zzz_yyzz_1,   \
                             ta1_x_zzz_yzzz_1,   \
                             ta1_x_zzz_zzzz_1,   \
                             ta2_xz_zz_xxxx_0,   \
                             ta2_xz_zz_xxxx_1,   \
                             ta2_xz_zz_xxxy_0,   \
                             ta2_xz_zz_xxxy_1,   \
                             ta2_xz_zz_xxxz_0,   \
                             ta2_xz_zz_xxxz_1,   \
                             ta2_xz_zz_xxyy_0,   \
                             ta2_xz_zz_xxyy_1,   \
                             ta2_xz_zz_xxyz_0,   \
                             ta2_xz_zz_xxyz_1,   \
                             ta2_xz_zz_xxzz_0,   \
                             ta2_xz_zz_xxzz_1,   \
                             ta2_xz_zz_xyyy_0,   \
                             ta2_xz_zz_xyyy_1,   \
                             ta2_xz_zz_xyyz_0,   \
                             ta2_xz_zz_xyyz_1,   \
                             ta2_xz_zz_xyzz_0,   \
                             ta2_xz_zz_xyzz_1,   \
                             ta2_xz_zz_xzzz_0,   \
                             ta2_xz_zz_xzzz_1,   \
                             ta2_xz_zz_yyyy_0,   \
                             ta2_xz_zz_yyyy_1,   \
                             ta2_xz_zz_yyyz_0,   \
                             ta2_xz_zz_yyyz_1,   \
                             ta2_xz_zz_yyzz_0,   \
                             ta2_xz_zz_yyzz_1,   \
                             ta2_xz_zz_yzzz_0,   \
                             ta2_xz_zz_yzzz_1,   \
                             ta2_xz_zz_zzzz_0,   \
                             ta2_xz_zz_zzzz_1,   \
                             ta2_xz_zzz_xxx_0,   \
                             ta2_xz_zzz_xxx_1,   \
                             ta2_xz_zzz_xxxx_0,  \
                             ta2_xz_zzz_xxxx_1,  \
                             ta2_xz_zzz_xxxy_0,  \
                             ta2_xz_zzz_xxxy_1,  \
                             ta2_xz_zzz_xxxz_0,  \
                             ta2_xz_zzz_xxxz_1,  \
                             ta2_xz_zzz_xxy_0,   \
                             ta2_xz_zzz_xxy_1,   \
                             ta2_xz_zzz_xxyy_0,  \
                             ta2_xz_zzz_xxyy_1,  \
                             ta2_xz_zzz_xxyz_0,  \
                             ta2_xz_zzz_xxyz_1,  \
                             ta2_xz_zzz_xxz_0,   \
                             ta2_xz_zzz_xxz_1,   \
                             ta2_xz_zzz_xxzz_0,  \
                             ta2_xz_zzz_xxzz_1,  \
                             ta2_xz_zzz_xyy_0,   \
                             ta2_xz_zzz_xyy_1,   \
                             ta2_xz_zzz_xyyy_0,  \
                             ta2_xz_zzz_xyyy_1,  \
                             ta2_xz_zzz_xyyz_0,  \
                             ta2_xz_zzz_xyyz_1,  \
                             ta2_xz_zzz_xyz_0,   \
                             ta2_xz_zzz_xyz_1,   \
                             ta2_xz_zzz_xyzz_0,  \
                             ta2_xz_zzz_xyzz_1,  \
                             ta2_xz_zzz_xzz_0,   \
                             ta2_xz_zzz_xzz_1,   \
                             ta2_xz_zzz_xzzz_0,  \
                             ta2_xz_zzz_xzzz_1,  \
                             ta2_xz_zzz_yyy_0,   \
                             ta2_xz_zzz_yyy_1,   \
                             ta2_xz_zzz_yyyy_0,  \
                             ta2_xz_zzz_yyyy_1,  \
                             ta2_xz_zzz_yyyz_0,  \
                             ta2_xz_zzz_yyyz_1,  \
                             ta2_xz_zzz_yyz_0,   \
                             ta2_xz_zzz_yyz_1,   \
                             ta2_xz_zzz_yyzz_0,  \
                             ta2_xz_zzz_yyzz_1,  \
                             ta2_xz_zzz_yzz_0,   \
                             ta2_xz_zzz_yzz_1,   \
                             ta2_xz_zzz_yzzz_0,  \
                             ta2_xz_zzz_yzzz_1,  \
                             ta2_xz_zzz_zzz_0,   \
                             ta2_xz_zzz_zzz_1,   \
                             ta2_xz_zzz_zzzz_0,  \
                             ta2_xz_zzz_zzzz_1,  \
                             ta2_xz_zzzz_xxxx_0, \
                             ta2_xz_zzzz_xxxy_0, \
                             ta2_xz_zzzz_xxxz_0, \
                             ta2_xz_zzzz_xxyy_0, \
                             ta2_xz_zzzz_xxyz_0, \
                             ta2_xz_zzzz_xxzz_0, \
                             ta2_xz_zzzz_xyyy_0, \
                             ta2_xz_zzzz_xyyz_0, \
                             ta2_xz_zzzz_xyzz_0, \
                             ta2_xz_zzzz_xzzz_0, \
                             ta2_xz_zzzz_yyyy_0, \
                             ta2_xz_zzzz_yyyz_0, \
                             ta2_xz_zzzz_yyzz_0, \
                             ta2_xz_zzzz_yzzz_0, \
                             ta2_xz_zzzz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_xz_zzzz_xxxx_0[i] = 3.0 * ta2_xz_zz_xxxx_0[i] * fe_0 - 3.0 * ta2_xz_zz_xxxx_1[i] * fe_0 + ta1_x_zzz_xxxx_1[i] +
                                ta2_xz_zzz_xxxx_0[i] * pa_z[i] - ta2_xz_zzz_xxxx_1[i] * pc_z[i];

        ta2_xz_zzzz_xxxy_0[i] = 3.0 * ta2_xz_zz_xxxy_0[i] * fe_0 - 3.0 * ta2_xz_zz_xxxy_1[i] * fe_0 + ta1_x_zzz_xxxy_1[i] +
                                ta2_xz_zzz_xxxy_0[i] * pa_z[i] - ta2_xz_zzz_xxxy_1[i] * pc_z[i];

        ta2_xz_zzzz_xxxz_0[i] = 3.0 * ta2_xz_zz_xxxz_0[i] * fe_0 - 3.0 * ta2_xz_zz_xxxz_1[i] * fe_0 + ta2_xz_zzz_xxx_0[i] * fe_0 -
                                ta2_xz_zzz_xxx_1[i] * fe_0 + ta1_x_zzz_xxxz_1[i] + ta2_xz_zzz_xxxz_0[i] * pa_z[i] - ta2_xz_zzz_xxxz_1[i] * pc_z[i];

        ta2_xz_zzzz_xxyy_0[i] = 3.0 * ta2_xz_zz_xxyy_0[i] * fe_0 - 3.0 * ta2_xz_zz_xxyy_1[i] * fe_0 + ta1_x_zzz_xxyy_1[i] +
                                ta2_xz_zzz_xxyy_0[i] * pa_z[i] - ta2_xz_zzz_xxyy_1[i] * pc_z[i];

        ta2_xz_zzzz_xxyz_0[i] = 3.0 * ta2_xz_zz_xxyz_0[i] * fe_0 - 3.0 * ta2_xz_zz_xxyz_1[i] * fe_0 + ta2_xz_zzz_xxy_0[i] * fe_0 -
                                ta2_xz_zzz_xxy_1[i] * fe_0 + ta1_x_zzz_xxyz_1[i] + ta2_xz_zzz_xxyz_0[i] * pa_z[i] - ta2_xz_zzz_xxyz_1[i] * pc_z[i];

        ta2_xz_zzzz_xxzz_0[i] = 3.0 * ta2_xz_zz_xxzz_0[i] * fe_0 - 3.0 * ta2_xz_zz_xxzz_1[i] * fe_0 + 2.0 * ta2_xz_zzz_xxz_0[i] * fe_0 -
                                2.0 * ta2_xz_zzz_xxz_1[i] * fe_0 + ta1_x_zzz_xxzz_1[i] + ta2_xz_zzz_xxzz_0[i] * pa_z[i] -
                                ta2_xz_zzz_xxzz_1[i] * pc_z[i];

        ta2_xz_zzzz_xyyy_0[i] = 3.0 * ta2_xz_zz_xyyy_0[i] * fe_0 - 3.0 * ta2_xz_zz_xyyy_1[i] * fe_0 + ta1_x_zzz_xyyy_1[i] +
                                ta2_xz_zzz_xyyy_0[i] * pa_z[i] - ta2_xz_zzz_xyyy_1[i] * pc_z[i];

        ta2_xz_zzzz_xyyz_0[i] = 3.0 * ta2_xz_zz_xyyz_0[i] * fe_0 - 3.0 * ta2_xz_zz_xyyz_1[i] * fe_0 + ta2_xz_zzz_xyy_0[i] * fe_0 -
                                ta2_xz_zzz_xyy_1[i] * fe_0 + ta1_x_zzz_xyyz_1[i] + ta2_xz_zzz_xyyz_0[i] * pa_z[i] - ta2_xz_zzz_xyyz_1[i] * pc_z[i];

        ta2_xz_zzzz_xyzz_0[i] = 3.0 * ta2_xz_zz_xyzz_0[i] * fe_0 - 3.0 * ta2_xz_zz_xyzz_1[i] * fe_0 + 2.0 * ta2_xz_zzz_xyz_0[i] * fe_0 -
                                2.0 * ta2_xz_zzz_xyz_1[i] * fe_0 + ta1_x_zzz_xyzz_1[i] + ta2_xz_zzz_xyzz_0[i] * pa_z[i] -
                                ta2_xz_zzz_xyzz_1[i] * pc_z[i];

        ta2_xz_zzzz_xzzz_0[i] = 3.0 * ta2_xz_zz_xzzz_0[i] * fe_0 - 3.0 * ta2_xz_zz_xzzz_1[i] * fe_0 + 3.0 * ta2_xz_zzz_xzz_0[i] * fe_0 -
                                3.0 * ta2_xz_zzz_xzz_1[i] * fe_0 + ta1_x_zzz_xzzz_1[i] + ta2_xz_zzz_xzzz_0[i] * pa_z[i] -
                                ta2_xz_zzz_xzzz_1[i] * pc_z[i];

        ta2_xz_zzzz_yyyy_0[i] = 3.0 * ta2_xz_zz_yyyy_0[i] * fe_0 - 3.0 * ta2_xz_zz_yyyy_1[i] * fe_0 + ta1_x_zzz_yyyy_1[i] +
                                ta2_xz_zzz_yyyy_0[i] * pa_z[i] - ta2_xz_zzz_yyyy_1[i] * pc_z[i];

        ta2_xz_zzzz_yyyz_0[i] = 3.0 * ta2_xz_zz_yyyz_0[i] * fe_0 - 3.0 * ta2_xz_zz_yyyz_1[i] * fe_0 + ta2_xz_zzz_yyy_0[i] * fe_0 -
                                ta2_xz_zzz_yyy_1[i] * fe_0 + ta1_x_zzz_yyyz_1[i] + ta2_xz_zzz_yyyz_0[i] * pa_z[i] - ta2_xz_zzz_yyyz_1[i] * pc_z[i];

        ta2_xz_zzzz_yyzz_0[i] = 3.0 * ta2_xz_zz_yyzz_0[i] * fe_0 - 3.0 * ta2_xz_zz_yyzz_1[i] * fe_0 + 2.0 * ta2_xz_zzz_yyz_0[i] * fe_0 -
                                2.0 * ta2_xz_zzz_yyz_1[i] * fe_0 + ta1_x_zzz_yyzz_1[i] + ta2_xz_zzz_yyzz_0[i] * pa_z[i] -
                                ta2_xz_zzz_yyzz_1[i] * pc_z[i];

        ta2_xz_zzzz_yzzz_0[i] = 3.0 * ta2_xz_zz_yzzz_0[i] * fe_0 - 3.0 * ta2_xz_zz_yzzz_1[i] * fe_0 + 3.0 * ta2_xz_zzz_yzz_0[i] * fe_0 -
                                3.0 * ta2_xz_zzz_yzz_1[i] * fe_0 + ta1_x_zzz_yzzz_1[i] + ta2_xz_zzz_yzzz_0[i] * pa_z[i] -
                                ta2_xz_zzz_yzzz_1[i] * pc_z[i];

        ta2_xz_zzzz_zzzz_0[i] = 3.0 * ta2_xz_zz_zzzz_0[i] * fe_0 - 3.0 * ta2_xz_zz_zzzz_1[i] * fe_0 + 4.0 * ta2_xz_zzz_zzz_0[i] * fe_0 -
                                4.0 * ta2_xz_zzz_zzz_1[i] * fe_0 + ta1_x_zzz_zzzz_1[i] + ta2_xz_zzz_zzzz_0[i] * pa_z[i] -
                                ta2_xz_zzz_zzzz_1[i] * pc_z[i];
    }

    // Set up 675-690 components of targeted buffer : GG

    auto ta2_yy_xxxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 675);

    auto ta2_yy_xxxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 676);

    auto ta2_yy_xxxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 677);

    auto ta2_yy_xxxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 678);

    auto ta2_yy_xxxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 679);

    auto ta2_yy_xxxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 680);

    auto ta2_yy_xxxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 681);

    auto ta2_yy_xxxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 682);

    auto ta2_yy_xxxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 683);

    auto ta2_yy_xxxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 684);

    auto ta2_yy_xxxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 685);

    auto ta2_yy_xxxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 686);

    auto ta2_yy_xxxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 687);

    auto ta2_yy_xxxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 688);

    auto ta2_yy_xxxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 689);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta2_yy_xx_xxxx_0,   \
                             ta2_yy_xx_xxxx_1,   \
                             ta2_yy_xx_xxxy_0,   \
                             ta2_yy_xx_xxxy_1,   \
                             ta2_yy_xx_xxxz_0,   \
                             ta2_yy_xx_xxxz_1,   \
                             ta2_yy_xx_xxyy_0,   \
                             ta2_yy_xx_xxyy_1,   \
                             ta2_yy_xx_xxyz_0,   \
                             ta2_yy_xx_xxyz_1,   \
                             ta2_yy_xx_xxzz_0,   \
                             ta2_yy_xx_xxzz_1,   \
                             ta2_yy_xx_xyyy_0,   \
                             ta2_yy_xx_xyyy_1,   \
                             ta2_yy_xx_xyyz_0,   \
                             ta2_yy_xx_xyyz_1,   \
                             ta2_yy_xx_xyzz_0,   \
                             ta2_yy_xx_xyzz_1,   \
                             ta2_yy_xx_xzzz_0,   \
                             ta2_yy_xx_xzzz_1,   \
                             ta2_yy_xx_yyyy_0,   \
                             ta2_yy_xx_yyyy_1,   \
                             ta2_yy_xx_yyyz_0,   \
                             ta2_yy_xx_yyyz_1,   \
                             ta2_yy_xx_yyzz_0,   \
                             ta2_yy_xx_yyzz_1,   \
                             ta2_yy_xx_yzzz_0,   \
                             ta2_yy_xx_yzzz_1,   \
                             ta2_yy_xx_zzzz_0,   \
                             ta2_yy_xx_zzzz_1,   \
                             ta2_yy_xxx_xxx_0,   \
                             ta2_yy_xxx_xxx_1,   \
                             ta2_yy_xxx_xxxx_0,  \
                             ta2_yy_xxx_xxxx_1,  \
                             ta2_yy_xxx_xxxy_0,  \
                             ta2_yy_xxx_xxxy_1,  \
                             ta2_yy_xxx_xxxz_0,  \
                             ta2_yy_xxx_xxxz_1,  \
                             ta2_yy_xxx_xxy_0,   \
                             ta2_yy_xxx_xxy_1,   \
                             ta2_yy_xxx_xxyy_0,  \
                             ta2_yy_xxx_xxyy_1,  \
                             ta2_yy_xxx_xxyz_0,  \
                             ta2_yy_xxx_xxyz_1,  \
                             ta2_yy_xxx_xxz_0,   \
                             ta2_yy_xxx_xxz_1,   \
                             ta2_yy_xxx_xxzz_0,  \
                             ta2_yy_xxx_xxzz_1,  \
                             ta2_yy_xxx_xyy_0,   \
                             ta2_yy_xxx_xyy_1,   \
                             ta2_yy_xxx_xyyy_0,  \
                             ta2_yy_xxx_xyyy_1,  \
                             ta2_yy_xxx_xyyz_0,  \
                             ta2_yy_xxx_xyyz_1,  \
                             ta2_yy_xxx_xyz_0,   \
                             ta2_yy_xxx_xyz_1,   \
                             ta2_yy_xxx_xyzz_0,  \
                             ta2_yy_xxx_xyzz_1,  \
                             ta2_yy_xxx_xzz_0,   \
                             ta2_yy_xxx_xzz_1,   \
                             ta2_yy_xxx_xzzz_0,  \
                             ta2_yy_xxx_xzzz_1,  \
                             ta2_yy_xxx_yyy_0,   \
                             ta2_yy_xxx_yyy_1,   \
                             ta2_yy_xxx_yyyy_0,  \
                             ta2_yy_xxx_yyyy_1,  \
                             ta2_yy_xxx_yyyz_0,  \
                             ta2_yy_xxx_yyyz_1,  \
                             ta2_yy_xxx_yyz_0,   \
                             ta2_yy_xxx_yyz_1,   \
                             ta2_yy_xxx_yyzz_0,  \
                             ta2_yy_xxx_yyzz_1,  \
                             ta2_yy_xxx_yzz_0,   \
                             ta2_yy_xxx_yzz_1,   \
                             ta2_yy_xxx_yzzz_0,  \
                             ta2_yy_xxx_yzzz_1,  \
                             ta2_yy_xxx_zzz_0,   \
                             ta2_yy_xxx_zzz_1,   \
                             ta2_yy_xxx_zzzz_0,  \
                             ta2_yy_xxx_zzzz_1,  \
                             ta2_yy_xxxx_xxxx_0, \
                             ta2_yy_xxxx_xxxy_0, \
                             ta2_yy_xxxx_xxxz_0, \
                             ta2_yy_xxxx_xxyy_0, \
                             ta2_yy_xxxx_xxyz_0, \
                             ta2_yy_xxxx_xxzz_0, \
                             ta2_yy_xxxx_xyyy_0, \
                             ta2_yy_xxxx_xyyz_0, \
                             ta2_yy_xxxx_xyzz_0, \
                             ta2_yy_xxxx_xzzz_0, \
                             ta2_yy_xxxx_yyyy_0, \
                             ta2_yy_xxxx_yyyz_0, \
                             ta2_yy_xxxx_yyzz_0, \
                             ta2_yy_xxxx_yzzz_0, \
                             ta2_yy_xxxx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxxx_xxxx_0[i] = 3.0 * ta2_yy_xx_xxxx_0[i] * fe_0 - 3.0 * ta2_yy_xx_xxxx_1[i] * fe_0 + 4.0 * ta2_yy_xxx_xxx_0[i] * fe_0 -
                                4.0 * ta2_yy_xxx_xxx_1[i] * fe_0 + ta2_yy_xxx_xxxx_0[i] * pa_x[i] - ta2_yy_xxx_xxxx_1[i] * pc_x[i];

        ta2_yy_xxxx_xxxy_0[i] = 3.0 * ta2_yy_xx_xxxy_0[i] * fe_0 - 3.0 * ta2_yy_xx_xxxy_1[i] * fe_0 + 3.0 * ta2_yy_xxx_xxy_0[i] * fe_0 -
                                3.0 * ta2_yy_xxx_xxy_1[i] * fe_0 + ta2_yy_xxx_xxxy_0[i] * pa_x[i] - ta2_yy_xxx_xxxy_1[i] * pc_x[i];

        ta2_yy_xxxx_xxxz_0[i] = 3.0 * ta2_yy_xx_xxxz_0[i] * fe_0 - 3.0 * ta2_yy_xx_xxxz_1[i] * fe_0 + 3.0 * ta2_yy_xxx_xxz_0[i] * fe_0 -
                                3.0 * ta2_yy_xxx_xxz_1[i] * fe_0 + ta2_yy_xxx_xxxz_0[i] * pa_x[i] - ta2_yy_xxx_xxxz_1[i] * pc_x[i];

        ta2_yy_xxxx_xxyy_0[i] = 3.0 * ta2_yy_xx_xxyy_0[i] * fe_0 - 3.0 * ta2_yy_xx_xxyy_1[i] * fe_0 + 2.0 * ta2_yy_xxx_xyy_0[i] * fe_0 -
                                2.0 * ta2_yy_xxx_xyy_1[i] * fe_0 + ta2_yy_xxx_xxyy_0[i] * pa_x[i] - ta2_yy_xxx_xxyy_1[i] * pc_x[i];

        ta2_yy_xxxx_xxyz_0[i] = 3.0 * ta2_yy_xx_xxyz_0[i] * fe_0 - 3.0 * ta2_yy_xx_xxyz_1[i] * fe_0 + 2.0 * ta2_yy_xxx_xyz_0[i] * fe_0 -
                                2.0 * ta2_yy_xxx_xyz_1[i] * fe_0 + ta2_yy_xxx_xxyz_0[i] * pa_x[i] - ta2_yy_xxx_xxyz_1[i] * pc_x[i];

        ta2_yy_xxxx_xxzz_0[i] = 3.0 * ta2_yy_xx_xxzz_0[i] * fe_0 - 3.0 * ta2_yy_xx_xxzz_1[i] * fe_0 + 2.0 * ta2_yy_xxx_xzz_0[i] * fe_0 -
                                2.0 * ta2_yy_xxx_xzz_1[i] * fe_0 + ta2_yy_xxx_xxzz_0[i] * pa_x[i] - ta2_yy_xxx_xxzz_1[i] * pc_x[i];

        ta2_yy_xxxx_xyyy_0[i] = 3.0 * ta2_yy_xx_xyyy_0[i] * fe_0 - 3.0 * ta2_yy_xx_xyyy_1[i] * fe_0 + ta2_yy_xxx_yyy_0[i] * fe_0 -
                                ta2_yy_xxx_yyy_1[i] * fe_0 + ta2_yy_xxx_xyyy_0[i] * pa_x[i] - ta2_yy_xxx_xyyy_1[i] * pc_x[i];

        ta2_yy_xxxx_xyyz_0[i] = 3.0 * ta2_yy_xx_xyyz_0[i] * fe_0 - 3.0 * ta2_yy_xx_xyyz_1[i] * fe_0 + ta2_yy_xxx_yyz_0[i] * fe_0 -
                                ta2_yy_xxx_yyz_1[i] * fe_0 + ta2_yy_xxx_xyyz_0[i] * pa_x[i] - ta2_yy_xxx_xyyz_1[i] * pc_x[i];

        ta2_yy_xxxx_xyzz_0[i] = 3.0 * ta2_yy_xx_xyzz_0[i] * fe_0 - 3.0 * ta2_yy_xx_xyzz_1[i] * fe_0 + ta2_yy_xxx_yzz_0[i] * fe_0 -
                                ta2_yy_xxx_yzz_1[i] * fe_0 + ta2_yy_xxx_xyzz_0[i] * pa_x[i] - ta2_yy_xxx_xyzz_1[i] * pc_x[i];

        ta2_yy_xxxx_xzzz_0[i] = 3.0 * ta2_yy_xx_xzzz_0[i] * fe_0 - 3.0 * ta2_yy_xx_xzzz_1[i] * fe_0 + ta2_yy_xxx_zzz_0[i] * fe_0 -
                                ta2_yy_xxx_zzz_1[i] * fe_0 + ta2_yy_xxx_xzzz_0[i] * pa_x[i] - ta2_yy_xxx_xzzz_1[i] * pc_x[i];

        ta2_yy_xxxx_yyyy_0[i] =
            3.0 * ta2_yy_xx_yyyy_0[i] * fe_0 - 3.0 * ta2_yy_xx_yyyy_1[i] * fe_0 + ta2_yy_xxx_yyyy_0[i] * pa_x[i] - ta2_yy_xxx_yyyy_1[i] * pc_x[i];

        ta2_yy_xxxx_yyyz_0[i] =
            3.0 * ta2_yy_xx_yyyz_0[i] * fe_0 - 3.0 * ta2_yy_xx_yyyz_1[i] * fe_0 + ta2_yy_xxx_yyyz_0[i] * pa_x[i] - ta2_yy_xxx_yyyz_1[i] * pc_x[i];

        ta2_yy_xxxx_yyzz_0[i] =
            3.0 * ta2_yy_xx_yyzz_0[i] * fe_0 - 3.0 * ta2_yy_xx_yyzz_1[i] * fe_0 + ta2_yy_xxx_yyzz_0[i] * pa_x[i] - ta2_yy_xxx_yyzz_1[i] * pc_x[i];

        ta2_yy_xxxx_yzzz_0[i] =
            3.0 * ta2_yy_xx_yzzz_0[i] * fe_0 - 3.0 * ta2_yy_xx_yzzz_1[i] * fe_0 + ta2_yy_xxx_yzzz_0[i] * pa_x[i] - ta2_yy_xxx_yzzz_1[i] * pc_x[i];

        ta2_yy_xxxx_zzzz_0[i] =
            3.0 * ta2_yy_xx_zzzz_0[i] * fe_0 - 3.0 * ta2_yy_xx_zzzz_1[i] * fe_0 + ta2_yy_xxx_zzzz_0[i] * pa_x[i] - ta2_yy_xxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 690-705 components of targeted buffer : GG

    auto ta2_yy_xxxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 690);

    auto ta2_yy_xxxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 691);

    auto ta2_yy_xxxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 692);

    auto ta2_yy_xxxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 693);

    auto ta2_yy_xxxy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 694);

    auto ta2_yy_xxxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 695);

    auto ta2_yy_xxxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 696);

    auto ta2_yy_xxxy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 697);

    auto ta2_yy_xxxy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 698);

    auto ta2_yy_xxxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 699);

    auto ta2_yy_xxxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 700);

    auto ta2_yy_xxxy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 701);

    auto ta2_yy_xxxy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 702);

    auto ta2_yy_xxxy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 703);

    auto ta2_yy_xxxy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 704);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_y_xxx_xxxx_1,   \
                             ta1_y_xxx_xxxy_1,   \
                             ta1_y_xxx_xxxz_1,   \
                             ta1_y_xxx_xxyy_1,   \
                             ta1_y_xxx_xxyz_1,   \
                             ta1_y_xxx_xxzz_1,   \
                             ta1_y_xxx_xyyy_1,   \
                             ta1_y_xxx_xyyz_1,   \
                             ta1_y_xxx_xyzz_1,   \
                             ta1_y_xxx_xzzz_1,   \
                             ta1_y_xxx_zzzz_1,   \
                             ta2_yy_xxx_xxx_0,   \
                             ta2_yy_xxx_xxx_1,   \
                             ta2_yy_xxx_xxxx_0,  \
                             ta2_yy_xxx_xxxx_1,  \
                             ta2_yy_xxx_xxxy_0,  \
                             ta2_yy_xxx_xxxy_1,  \
                             ta2_yy_xxx_xxxz_0,  \
                             ta2_yy_xxx_xxxz_1,  \
                             ta2_yy_xxx_xxy_0,   \
                             ta2_yy_xxx_xxy_1,   \
                             ta2_yy_xxx_xxyy_0,  \
                             ta2_yy_xxx_xxyy_1,  \
                             ta2_yy_xxx_xxyz_0,  \
                             ta2_yy_xxx_xxyz_1,  \
                             ta2_yy_xxx_xxz_0,   \
                             ta2_yy_xxx_xxz_1,   \
                             ta2_yy_xxx_xxzz_0,  \
                             ta2_yy_xxx_xxzz_1,  \
                             ta2_yy_xxx_xyy_0,   \
                             ta2_yy_xxx_xyy_1,   \
                             ta2_yy_xxx_xyyy_0,  \
                             ta2_yy_xxx_xyyy_1,  \
                             ta2_yy_xxx_xyyz_0,  \
                             ta2_yy_xxx_xyyz_1,  \
                             ta2_yy_xxx_xyz_0,   \
                             ta2_yy_xxx_xyz_1,   \
                             ta2_yy_xxx_xyzz_0,  \
                             ta2_yy_xxx_xyzz_1,  \
                             ta2_yy_xxx_xzz_0,   \
                             ta2_yy_xxx_xzz_1,   \
                             ta2_yy_xxx_xzzz_0,  \
                             ta2_yy_xxx_xzzz_1,  \
                             ta2_yy_xxx_zzzz_0,  \
                             ta2_yy_xxx_zzzz_1,  \
                             ta2_yy_xxxy_xxxx_0, \
                             ta2_yy_xxxy_xxxy_0, \
                             ta2_yy_xxxy_xxxz_0, \
                             ta2_yy_xxxy_xxyy_0, \
                             ta2_yy_xxxy_xxyz_0, \
                             ta2_yy_xxxy_xxzz_0, \
                             ta2_yy_xxxy_xyyy_0, \
                             ta2_yy_xxxy_xyyz_0, \
                             ta2_yy_xxxy_xyzz_0, \
                             ta2_yy_xxxy_xzzz_0, \
                             ta2_yy_xxxy_yyyy_0, \
                             ta2_yy_xxxy_yyyz_0, \
                             ta2_yy_xxxy_yyzz_0, \
                             ta2_yy_xxxy_yzzz_0, \
                             ta2_yy_xxxy_zzzz_0, \
                             ta2_yy_xxy_yyyy_0,  \
                             ta2_yy_xxy_yyyy_1,  \
                             ta2_yy_xxy_yyyz_0,  \
                             ta2_yy_xxy_yyyz_1,  \
                             ta2_yy_xxy_yyzz_0,  \
                             ta2_yy_xxy_yyzz_1,  \
                             ta2_yy_xxy_yzzz_0,  \
                             ta2_yy_xxy_yzzz_1,  \
                             ta2_yy_xy_yyyy_0,   \
                             ta2_yy_xy_yyyy_1,   \
                             ta2_yy_xy_yyyz_0,   \
                             ta2_yy_xy_yyyz_1,   \
                             ta2_yy_xy_yyzz_0,   \
                             ta2_yy_xy_yyzz_1,   \
                             ta2_yy_xy_yzzz_0,   \
                             ta2_yy_xy_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxxy_xxxx_0[i] = 2.0 * ta1_y_xxx_xxxx_1[i] + ta2_yy_xxx_xxxx_0[i] * pa_y[i] - ta2_yy_xxx_xxxx_1[i] * pc_y[i];

        ta2_yy_xxxy_xxxy_0[i] = ta2_yy_xxx_xxx_0[i] * fe_0 - ta2_yy_xxx_xxx_1[i] * fe_0 + 2.0 * ta1_y_xxx_xxxy_1[i] + ta2_yy_xxx_xxxy_0[i] * pa_y[i] -
                                ta2_yy_xxx_xxxy_1[i] * pc_y[i];

        ta2_yy_xxxy_xxxz_0[i] = 2.0 * ta1_y_xxx_xxxz_1[i] + ta2_yy_xxx_xxxz_0[i] * pa_y[i] - ta2_yy_xxx_xxxz_1[i] * pc_y[i];

        ta2_yy_xxxy_xxyy_0[i] = 2.0 * ta2_yy_xxx_xxy_0[i] * fe_0 - 2.0 * ta2_yy_xxx_xxy_1[i] * fe_0 + 2.0 * ta1_y_xxx_xxyy_1[i] +
                                ta2_yy_xxx_xxyy_0[i] * pa_y[i] - ta2_yy_xxx_xxyy_1[i] * pc_y[i];

        ta2_yy_xxxy_xxyz_0[i] = ta2_yy_xxx_xxz_0[i] * fe_0 - ta2_yy_xxx_xxz_1[i] * fe_0 + 2.0 * ta1_y_xxx_xxyz_1[i] + ta2_yy_xxx_xxyz_0[i] * pa_y[i] -
                                ta2_yy_xxx_xxyz_1[i] * pc_y[i];

        ta2_yy_xxxy_xxzz_0[i] = 2.0 * ta1_y_xxx_xxzz_1[i] + ta2_yy_xxx_xxzz_0[i] * pa_y[i] - ta2_yy_xxx_xxzz_1[i] * pc_y[i];

        ta2_yy_xxxy_xyyy_0[i] = 3.0 * ta2_yy_xxx_xyy_0[i] * fe_0 - 3.0 * ta2_yy_xxx_xyy_1[i] * fe_0 + 2.0 * ta1_y_xxx_xyyy_1[i] +
                                ta2_yy_xxx_xyyy_0[i] * pa_y[i] - ta2_yy_xxx_xyyy_1[i] * pc_y[i];

        ta2_yy_xxxy_xyyz_0[i] = 2.0 * ta2_yy_xxx_xyz_0[i] * fe_0 - 2.0 * ta2_yy_xxx_xyz_1[i] * fe_0 + 2.0 * ta1_y_xxx_xyyz_1[i] +
                                ta2_yy_xxx_xyyz_0[i] * pa_y[i] - ta2_yy_xxx_xyyz_1[i] * pc_y[i];

        ta2_yy_xxxy_xyzz_0[i] = ta2_yy_xxx_xzz_0[i] * fe_0 - ta2_yy_xxx_xzz_1[i] * fe_0 + 2.0 * ta1_y_xxx_xyzz_1[i] + ta2_yy_xxx_xyzz_0[i] * pa_y[i] -
                                ta2_yy_xxx_xyzz_1[i] * pc_y[i];

        ta2_yy_xxxy_xzzz_0[i] = 2.0 * ta1_y_xxx_xzzz_1[i] + ta2_yy_xxx_xzzz_0[i] * pa_y[i] - ta2_yy_xxx_xzzz_1[i] * pc_y[i];

        ta2_yy_xxxy_yyyy_0[i] =
            2.0 * ta2_yy_xy_yyyy_0[i] * fe_0 - 2.0 * ta2_yy_xy_yyyy_1[i] * fe_0 + ta2_yy_xxy_yyyy_0[i] * pa_x[i] - ta2_yy_xxy_yyyy_1[i] * pc_x[i];

        ta2_yy_xxxy_yyyz_0[i] =
            2.0 * ta2_yy_xy_yyyz_0[i] * fe_0 - 2.0 * ta2_yy_xy_yyyz_1[i] * fe_0 + ta2_yy_xxy_yyyz_0[i] * pa_x[i] - ta2_yy_xxy_yyyz_1[i] * pc_x[i];

        ta2_yy_xxxy_yyzz_0[i] =
            2.0 * ta2_yy_xy_yyzz_0[i] * fe_0 - 2.0 * ta2_yy_xy_yyzz_1[i] * fe_0 + ta2_yy_xxy_yyzz_0[i] * pa_x[i] - ta2_yy_xxy_yyzz_1[i] * pc_x[i];

        ta2_yy_xxxy_yzzz_0[i] =
            2.0 * ta2_yy_xy_yzzz_0[i] * fe_0 - 2.0 * ta2_yy_xy_yzzz_1[i] * fe_0 + ta2_yy_xxy_yzzz_0[i] * pa_x[i] - ta2_yy_xxy_yzzz_1[i] * pc_x[i];

        ta2_yy_xxxy_zzzz_0[i] = 2.0 * ta1_y_xxx_zzzz_1[i] + ta2_yy_xxx_zzzz_0[i] * pa_y[i] - ta2_yy_xxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 705-720 components of targeted buffer : GG

    auto ta2_yy_xxxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 705);

    auto ta2_yy_xxxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 706);

    auto ta2_yy_xxxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 707);

    auto ta2_yy_xxxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 708);

    auto ta2_yy_xxxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 709);

    auto ta2_yy_xxxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 710);

    auto ta2_yy_xxxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 711);

    auto ta2_yy_xxxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 712);

    auto ta2_yy_xxxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 713);

    auto ta2_yy_xxxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 714);

    auto ta2_yy_xxxz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 715);

    auto ta2_yy_xxxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 716);

    auto ta2_yy_xxxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 717);

    auto ta2_yy_xxxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 718);

    auto ta2_yy_xxxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 719);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta2_yy_xxx_xxx_0,   \
                             ta2_yy_xxx_xxx_1,   \
                             ta2_yy_xxx_xxxx_0,  \
                             ta2_yy_xxx_xxxx_1,  \
                             ta2_yy_xxx_xxxy_0,  \
                             ta2_yy_xxx_xxxy_1,  \
                             ta2_yy_xxx_xxxz_0,  \
                             ta2_yy_xxx_xxxz_1,  \
                             ta2_yy_xxx_xxy_0,   \
                             ta2_yy_xxx_xxy_1,   \
                             ta2_yy_xxx_xxyy_0,  \
                             ta2_yy_xxx_xxyy_1,  \
                             ta2_yy_xxx_xxyz_0,  \
                             ta2_yy_xxx_xxyz_1,  \
                             ta2_yy_xxx_xxz_0,   \
                             ta2_yy_xxx_xxz_1,   \
                             ta2_yy_xxx_xxzz_0,  \
                             ta2_yy_xxx_xxzz_1,  \
                             ta2_yy_xxx_xyy_0,   \
                             ta2_yy_xxx_xyy_1,   \
                             ta2_yy_xxx_xyyy_0,  \
                             ta2_yy_xxx_xyyy_1,  \
                             ta2_yy_xxx_xyyz_0,  \
                             ta2_yy_xxx_xyyz_1,  \
                             ta2_yy_xxx_xyz_0,   \
                             ta2_yy_xxx_xyz_1,   \
                             ta2_yy_xxx_xyzz_0,  \
                             ta2_yy_xxx_xyzz_1,  \
                             ta2_yy_xxx_xzz_0,   \
                             ta2_yy_xxx_xzz_1,   \
                             ta2_yy_xxx_xzzz_0,  \
                             ta2_yy_xxx_xzzz_1,  \
                             ta2_yy_xxx_yyyy_0,  \
                             ta2_yy_xxx_yyyy_1,  \
                             ta2_yy_xxxz_xxxx_0, \
                             ta2_yy_xxxz_xxxy_0, \
                             ta2_yy_xxxz_xxxz_0, \
                             ta2_yy_xxxz_xxyy_0, \
                             ta2_yy_xxxz_xxyz_0, \
                             ta2_yy_xxxz_xxzz_0, \
                             ta2_yy_xxxz_xyyy_0, \
                             ta2_yy_xxxz_xyyz_0, \
                             ta2_yy_xxxz_xyzz_0, \
                             ta2_yy_xxxz_xzzz_0, \
                             ta2_yy_xxxz_yyyy_0, \
                             ta2_yy_xxxz_yyyz_0, \
                             ta2_yy_xxxz_yyzz_0, \
                             ta2_yy_xxxz_yzzz_0, \
                             ta2_yy_xxxz_zzzz_0, \
                             ta2_yy_xxz_yyyz_0,  \
                             ta2_yy_xxz_yyyz_1,  \
                             ta2_yy_xxz_yyzz_0,  \
                             ta2_yy_xxz_yyzz_1,  \
                             ta2_yy_xxz_yzzz_0,  \
                             ta2_yy_xxz_yzzz_1,  \
                             ta2_yy_xxz_zzzz_0,  \
                             ta2_yy_xxz_zzzz_1,  \
                             ta2_yy_xz_yyyz_0,   \
                             ta2_yy_xz_yyyz_1,   \
                             ta2_yy_xz_yyzz_0,   \
                             ta2_yy_xz_yyzz_1,   \
                             ta2_yy_xz_yzzz_0,   \
                             ta2_yy_xz_yzzz_1,   \
                             ta2_yy_xz_zzzz_0,   \
                             ta2_yy_xz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxxz_xxxx_0[i] = ta2_yy_xxx_xxxx_0[i] * pa_z[i] - ta2_yy_xxx_xxxx_1[i] * pc_z[i];

        ta2_yy_xxxz_xxxy_0[i] = ta2_yy_xxx_xxxy_0[i] * pa_z[i] - ta2_yy_xxx_xxxy_1[i] * pc_z[i];

        ta2_yy_xxxz_xxxz_0[i] =
            ta2_yy_xxx_xxx_0[i] * fe_0 - ta2_yy_xxx_xxx_1[i] * fe_0 + ta2_yy_xxx_xxxz_0[i] * pa_z[i] - ta2_yy_xxx_xxxz_1[i] * pc_z[i];

        ta2_yy_xxxz_xxyy_0[i] = ta2_yy_xxx_xxyy_0[i] * pa_z[i] - ta2_yy_xxx_xxyy_1[i] * pc_z[i];

        ta2_yy_xxxz_xxyz_0[i] =
            ta2_yy_xxx_xxy_0[i] * fe_0 - ta2_yy_xxx_xxy_1[i] * fe_0 + ta2_yy_xxx_xxyz_0[i] * pa_z[i] - ta2_yy_xxx_xxyz_1[i] * pc_z[i];

        ta2_yy_xxxz_xxzz_0[i] =
            2.0 * ta2_yy_xxx_xxz_0[i] * fe_0 - 2.0 * ta2_yy_xxx_xxz_1[i] * fe_0 + ta2_yy_xxx_xxzz_0[i] * pa_z[i] - ta2_yy_xxx_xxzz_1[i] * pc_z[i];

        ta2_yy_xxxz_xyyy_0[i] = ta2_yy_xxx_xyyy_0[i] * pa_z[i] - ta2_yy_xxx_xyyy_1[i] * pc_z[i];

        ta2_yy_xxxz_xyyz_0[i] =
            ta2_yy_xxx_xyy_0[i] * fe_0 - ta2_yy_xxx_xyy_1[i] * fe_0 + ta2_yy_xxx_xyyz_0[i] * pa_z[i] - ta2_yy_xxx_xyyz_1[i] * pc_z[i];

        ta2_yy_xxxz_xyzz_0[i] =
            2.0 * ta2_yy_xxx_xyz_0[i] * fe_0 - 2.0 * ta2_yy_xxx_xyz_1[i] * fe_0 + ta2_yy_xxx_xyzz_0[i] * pa_z[i] - ta2_yy_xxx_xyzz_1[i] * pc_z[i];

        ta2_yy_xxxz_xzzz_0[i] =
            3.0 * ta2_yy_xxx_xzz_0[i] * fe_0 - 3.0 * ta2_yy_xxx_xzz_1[i] * fe_0 + ta2_yy_xxx_xzzz_0[i] * pa_z[i] - ta2_yy_xxx_xzzz_1[i] * pc_z[i];

        ta2_yy_xxxz_yyyy_0[i] = ta2_yy_xxx_yyyy_0[i] * pa_z[i] - ta2_yy_xxx_yyyy_1[i] * pc_z[i];

        ta2_yy_xxxz_yyyz_0[i] =
            2.0 * ta2_yy_xz_yyyz_0[i] * fe_0 - 2.0 * ta2_yy_xz_yyyz_1[i] * fe_0 + ta2_yy_xxz_yyyz_0[i] * pa_x[i] - ta2_yy_xxz_yyyz_1[i] * pc_x[i];

        ta2_yy_xxxz_yyzz_0[i] =
            2.0 * ta2_yy_xz_yyzz_0[i] * fe_0 - 2.0 * ta2_yy_xz_yyzz_1[i] * fe_0 + ta2_yy_xxz_yyzz_0[i] * pa_x[i] - ta2_yy_xxz_yyzz_1[i] * pc_x[i];

        ta2_yy_xxxz_yzzz_0[i] =
            2.0 * ta2_yy_xz_yzzz_0[i] * fe_0 - 2.0 * ta2_yy_xz_yzzz_1[i] * fe_0 + ta2_yy_xxz_yzzz_0[i] * pa_x[i] - ta2_yy_xxz_yzzz_1[i] * pc_x[i];

        ta2_yy_xxxz_zzzz_0[i] =
            2.0 * ta2_yy_xz_zzzz_0[i] * fe_0 - 2.0 * ta2_yy_xz_zzzz_1[i] * fe_0 + ta2_yy_xxz_zzzz_0[i] * pa_x[i] - ta2_yy_xxz_zzzz_1[i] * pc_x[i];
    }

    // Set up 720-735 components of targeted buffer : GG

    auto ta2_yy_xxyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 720);

    auto ta2_yy_xxyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 721);

    auto ta2_yy_xxyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 722);

    auto ta2_yy_xxyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 723);

    auto ta2_yy_xxyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 724);

    auto ta2_yy_xxyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 725);

    auto ta2_yy_xxyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 726);

    auto ta2_yy_xxyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 727);

    auto ta2_yy_xxyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 728);

    auto ta2_yy_xxyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 729);

    auto ta2_yy_xxyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 730);

    auto ta2_yy_xxyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 731);

    auto ta2_yy_xxyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 732);

    auto ta2_yy_xxyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 733);

    auto ta2_yy_xxyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 734);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_y_xxy_xxxx_1,   \
                             ta1_y_xxy_xxxz_1,   \
                             ta1_y_xxy_xxzz_1,   \
                             ta1_y_xxy_xzzz_1,   \
                             ta2_yy_xx_xxxx_0,   \
                             ta2_yy_xx_xxxx_1,   \
                             ta2_yy_xx_xxxz_0,   \
                             ta2_yy_xx_xxxz_1,   \
                             ta2_yy_xx_xxzz_0,   \
                             ta2_yy_xx_xxzz_1,   \
                             ta2_yy_xx_xzzz_0,   \
                             ta2_yy_xx_xzzz_1,   \
                             ta2_yy_xxy_xxxx_0,  \
                             ta2_yy_xxy_xxxx_1,  \
                             ta2_yy_xxy_xxxz_0,  \
                             ta2_yy_xxy_xxxz_1,  \
                             ta2_yy_xxy_xxzz_0,  \
                             ta2_yy_xxy_xxzz_1,  \
                             ta2_yy_xxy_xzzz_0,  \
                             ta2_yy_xxy_xzzz_1,  \
                             ta2_yy_xxyy_xxxx_0, \
                             ta2_yy_xxyy_xxxy_0, \
                             ta2_yy_xxyy_xxxz_0, \
                             ta2_yy_xxyy_xxyy_0, \
                             ta2_yy_xxyy_xxyz_0, \
                             ta2_yy_xxyy_xxzz_0, \
                             ta2_yy_xxyy_xyyy_0, \
                             ta2_yy_xxyy_xyyz_0, \
                             ta2_yy_xxyy_xyzz_0, \
                             ta2_yy_xxyy_xzzz_0, \
                             ta2_yy_xxyy_yyyy_0, \
                             ta2_yy_xxyy_yyyz_0, \
                             ta2_yy_xxyy_yyzz_0, \
                             ta2_yy_xxyy_yzzz_0, \
                             ta2_yy_xxyy_zzzz_0, \
                             ta2_yy_xyy_xxxy_0,  \
                             ta2_yy_xyy_xxxy_1,  \
                             ta2_yy_xyy_xxy_0,   \
                             ta2_yy_xyy_xxy_1,   \
                             ta2_yy_xyy_xxyy_0,  \
                             ta2_yy_xyy_xxyy_1,  \
                             ta2_yy_xyy_xxyz_0,  \
                             ta2_yy_xyy_xxyz_1,  \
                             ta2_yy_xyy_xyy_0,   \
                             ta2_yy_xyy_xyy_1,   \
                             ta2_yy_xyy_xyyy_0,  \
                             ta2_yy_xyy_xyyy_1,  \
                             ta2_yy_xyy_xyyz_0,  \
                             ta2_yy_xyy_xyyz_1,  \
                             ta2_yy_xyy_xyz_0,   \
                             ta2_yy_xyy_xyz_1,   \
                             ta2_yy_xyy_xyzz_0,  \
                             ta2_yy_xyy_xyzz_1,  \
                             ta2_yy_xyy_yyy_0,   \
                             ta2_yy_xyy_yyy_1,   \
                             ta2_yy_xyy_yyyy_0,  \
                             ta2_yy_xyy_yyyy_1,  \
                             ta2_yy_xyy_yyyz_0,  \
                             ta2_yy_xyy_yyyz_1,  \
                             ta2_yy_xyy_yyz_0,   \
                             ta2_yy_xyy_yyz_1,   \
                             ta2_yy_xyy_yyzz_0,  \
                             ta2_yy_xyy_yyzz_1,  \
                             ta2_yy_xyy_yzz_0,   \
                             ta2_yy_xyy_yzz_1,   \
                             ta2_yy_xyy_yzzz_0,  \
                             ta2_yy_xyy_yzzz_1,  \
                             ta2_yy_xyy_zzzz_0,  \
                             ta2_yy_xyy_zzzz_1,  \
                             ta2_yy_yy_xxxy_0,   \
                             ta2_yy_yy_xxxy_1,   \
                             ta2_yy_yy_xxyy_0,   \
                             ta2_yy_yy_xxyy_1,   \
                             ta2_yy_yy_xxyz_0,   \
                             ta2_yy_yy_xxyz_1,   \
                             ta2_yy_yy_xyyy_0,   \
                             ta2_yy_yy_xyyy_1,   \
                             ta2_yy_yy_xyyz_0,   \
                             ta2_yy_yy_xyyz_1,   \
                             ta2_yy_yy_xyzz_0,   \
                             ta2_yy_yy_xyzz_1,   \
                             ta2_yy_yy_yyyy_0,   \
                             ta2_yy_yy_yyyy_1,   \
                             ta2_yy_yy_yyyz_0,   \
                             ta2_yy_yy_yyyz_1,   \
                             ta2_yy_yy_yyzz_0,   \
                             ta2_yy_yy_yyzz_1,   \
                             ta2_yy_yy_yzzz_0,   \
                             ta2_yy_yy_yzzz_1,   \
                             ta2_yy_yy_zzzz_0,   \
                             ta2_yy_yy_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxyy_xxxx_0[i] = ta2_yy_xx_xxxx_0[i] * fe_0 - ta2_yy_xx_xxxx_1[i] * fe_0 + 2.0 * ta1_y_xxy_xxxx_1[i] + ta2_yy_xxy_xxxx_0[i] * pa_y[i] -
                                ta2_yy_xxy_xxxx_1[i] * pc_y[i];

        ta2_yy_xxyy_xxxy_0[i] = ta2_yy_yy_xxxy_0[i] * fe_0 - ta2_yy_yy_xxxy_1[i] * fe_0 + 3.0 * ta2_yy_xyy_xxy_0[i] * fe_0 -
                                3.0 * ta2_yy_xyy_xxy_1[i] * fe_0 + ta2_yy_xyy_xxxy_0[i] * pa_x[i] - ta2_yy_xyy_xxxy_1[i] * pc_x[i];

        ta2_yy_xxyy_xxxz_0[i] = ta2_yy_xx_xxxz_0[i] * fe_0 - ta2_yy_xx_xxxz_1[i] * fe_0 + 2.0 * ta1_y_xxy_xxxz_1[i] + ta2_yy_xxy_xxxz_0[i] * pa_y[i] -
                                ta2_yy_xxy_xxxz_1[i] * pc_y[i];

        ta2_yy_xxyy_xxyy_0[i] = ta2_yy_yy_xxyy_0[i] * fe_0 - ta2_yy_yy_xxyy_1[i] * fe_0 + 2.0 * ta2_yy_xyy_xyy_0[i] * fe_0 -
                                2.0 * ta2_yy_xyy_xyy_1[i] * fe_0 + ta2_yy_xyy_xxyy_0[i] * pa_x[i] - ta2_yy_xyy_xxyy_1[i] * pc_x[i];

        ta2_yy_xxyy_xxyz_0[i] = ta2_yy_yy_xxyz_0[i] * fe_0 - ta2_yy_yy_xxyz_1[i] * fe_0 + 2.0 * ta2_yy_xyy_xyz_0[i] * fe_0 -
                                2.0 * ta2_yy_xyy_xyz_1[i] * fe_0 + ta2_yy_xyy_xxyz_0[i] * pa_x[i] - ta2_yy_xyy_xxyz_1[i] * pc_x[i];

        ta2_yy_xxyy_xxzz_0[i] = ta2_yy_xx_xxzz_0[i] * fe_0 - ta2_yy_xx_xxzz_1[i] * fe_0 + 2.0 * ta1_y_xxy_xxzz_1[i] + ta2_yy_xxy_xxzz_0[i] * pa_y[i] -
                                ta2_yy_xxy_xxzz_1[i] * pc_y[i];

        ta2_yy_xxyy_xyyy_0[i] = ta2_yy_yy_xyyy_0[i] * fe_0 - ta2_yy_yy_xyyy_1[i] * fe_0 + ta2_yy_xyy_yyy_0[i] * fe_0 - ta2_yy_xyy_yyy_1[i] * fe_0 +
                                ta2_yy_xyy_xyyy_0[i] * pa_x[i] - ta2_yy_xyy_xyyy_1[i] * pc_x[i];

        ta2_yy_xxyy_xyyz_0[i] = ta2_yy_yy_xyyz_0[i] * fe_0 - ta2_yy_yy_xyyz_1[i] * fe_0 + ta2_yy_xyy_yyz_0[i] * fe_0 - ta2_yy_xyy_yyz_1[i] * fe_0 +
                                ta2_yy_xyy_xyyz_0[i] * pa_x[i] - ta2_yy_xyy_xyyz_1[i] * pc_x[i];

        ta2_yy_xxyy_xyzz_0[i] = ta2_yy_yy_xyzz_0[i] * fe_0 - ta2_yy_yy_xyzz_1[i] * fe_0 + ta2_yy_xyy_yzz_0[i] * fe_0 - ta2_yy_xyy_yzz_1[i] * fe_0 +
                                ta2_yy_xyy_xyzz_0[i] * pa_x[i] - ta2_yy_xyy_xyzz_1[i] * pc_x[i];

        ta2_yy_xxyy_xzzz_0[i] = ta2_yy_xx_xzzz_0[i] * fe_0 - ta2_yy_xx_xzzz_1[i] * fe_0 + 2.0 * ta1_y_xxy_xzzz_1[i] + ta2_yy_xxy_xzzz_0[i] * pa_y[i] -
                                ta2_yy_xxy_xzzz_1[i] * pc_y[i];

        ta2_yy_xxyy_yyyy_0[i] =
            ta2_yy_yy_yyyy_0[i] * fe_0 - ta2_yy_yy_yyyy_1[i] * fe_0 + ta2_yy_xyy_yyyy_0[i] * pa_x[i] - ta2_yy_xyy_yyyy_1[i] * pc_x[i];

        ta2_yy_xxyy_yyyz_0[i] =
            ta2_yy_yy_yyyz_0[i] * fe_0 - ta2_yy_yy_yyyz_1[i] * fe_0 + ta2_yy_xyy_yyyz_0[i] * pa_x[i] - ta2_yy_xyy_yyyz_1[i] * pc_x[i];

        ta2_yy_xxyy_yyzz_0[i] =
            ta2_yy_yy_yyzz_0[i] * fe_0 - ta2_yy_yy_yyzz_1[i] * fe_0 + ta2_yy_xyy_yyzz_0[i] * pa_x[i] - ta2_yy_xyy_yyzz_1[i] * pc_x[i];

        ta2_yy_xxyy_yzzz_0[i] =
            ta2_yy_yy_yzzz_0[i] * fe_0 - ta2_yy_yy_yzzz_1[i] * fe_0 + ta2_yy_xyy_yzzz_0[i] * pa_x[i] - ta2_yy_xyy_yzzz_1[i] * pc_x[i];

        ta2_yy_xxyy_zzzz_0[i] =
            ta2_yy_yy_zzzz_0[i] * fe_0 - ta2_yy_yy_zzzz_1[i] * fe_0 + ta2_yy_xyy_zzzz_0[i] * pa_x[i] - ta2_yy_xyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 735-750 components of targeted buffer : GG

    auto ta2_yy_xxyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 735);

    auto ta2_yy_xxyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 736);

    auto ta2_yy_xxyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 737);

    auto ta2_yy_xxyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 738);

    auto ta2_yy_xxyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 739);

    auto ta2_yy_xxyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 740);

    auto ta2_yy_xxyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 741);

    auto ta2_yy_xxyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 742);

    auto ta2_yy_xxyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 743);

    auto ta2_yy_xxyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 744);

    auto ta2_yy_xxyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 745);

    auto ta2_yy_xxyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 746);

    auto ta2_yy_xxyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 747);

    auto ta2_yy_xxyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 748);

    auto ta2_yy_xxyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 749);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_xxz_xxxz_1,   \
                             ta1_y_xxz_xxzz_1,   \
                             ta1_y_xxz_xzzz_1,   \
                             ta1_y_xxz_zzzz_1,   \
                             ta2_yy_xxy_xxxx_0,  \
                             ta2_yy_xxy_xxxx_1,  \
                             ta2_yy_xxy_xxxy_0,  \
                             ta2_yy_xxy_xxxy_1,  \
                             ta2_yy_xxy_xxy_0,   \
                             ta2_yy_xxy_xxy_1,   \
                             ta2_yy_xxy_xxyy_0,  \
                             ta2_yy_xxy_xxyy_1,  \
                             ta2_yy_xxy_xxyz_0,  \
                             ta2_yy_xxy_xxyz_1,  \
                             ta2_yy_xxy_xyy_0,   \
                             ta2_yy_xxy_xyy_1,   \
                             ta2_yy_xxy_xyyy_0,  \
                             ta2_yy_xxy_xyyy_1,  \
                             ta2_yy_xxy_xyyz_0,  \
                             ta2_yy_xxy_xyyz_1,  \
                             ta2_yy_xxy_xyz_0,   \
                             ta2_yy_xxy_xyz_1,   \
                             ta2_yy_xxy_xyzz_0,  \
                             ta2_yy_xxy_xyzz_1,  \
                             ta2_yy_xxy_yyyy_0,  \
                             ta2_yy_xxy_yyyy_1,  \
                             ta2_yy_xxyz_xxxx_0, \
                             ta2_yy_xxyz_xxxy_0, \
                             ta2_yy_xxyz_xxxz_0, \
                             ta2_yy_xxyz_xxyy_0, \
                             ta2_yy_xxyz_xxyz_0, \
                             ta2_yy_xxyz_xxzz_0, \
                             ta2_yy_xxyz_xyyy_0, \
                             ta2_yy_xxyz_xyyz_0, \
                             ta2_yy_xxyz_xyzz_0, \
                             ta2_yy_xxyz_xzzz_0, \
                             ta2_yy_xxyz_yyyy_0, \
                             ta2_yy_xxyz_yyyz_0, \
                             ta2_yy_xxyz_yyzz_0, \
                             ta2_yy_xxyz_yzzz_0, \
                             ta2_yy_xxyz_zzzz_0, \
                             ta2_yy_xxz_xxxz_0,  \
                             ta2_yy_xxz_xxxz_1,  \
                             ta2_yy_xxz_xxzz_0,  \
                             ta2_yy_xxz_xxzz_1,  \
                             ta2_yy_xxz_xzzz_0,  \
                             ta2_yy_xxz_xzzz_1,  \
                             ta2_yy_xxz_zzzz_0,  \
                             ta2_yy_xxz_zzzz_1,  \
                             ta2_yy_xyz_yyyz_0,  \
                             ta2_yy_xyz_yyyz_1,  \
                             ta2_yy_xyz_yyzz_0,  \
                             ta2_yy_xyz_yyzz_1,  \
                             ta2_yy_xyz_yzzz_0,  \
                             ta2_yy_xyz_yzzz_1,  \
                             ta2_yy_yz_yyyz_0,   \
                             ta2_yy_yz_yyyz_1,   \
                             ta2_yy_yz_yyzz_0,   \
                             ta2_yy_yz_yyzz_1,   \
                             ta2_yy_yz_yzzz_0,   \
                             ta2_yy_yz_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxyz_xxxx_0[i] = ta2_yy_xxy_xxxx_0[i] * pa_z[i] - ta2_yy_xxy_xxxx_1[i] * pc_z[i];

        ta2_yy_xxyz_xxxy_0[i] = ta2_yy_xxy_xxxy_0[i] * pa_z[i] - ta2_yy_xxy_xxxy_1[i] * pc_z[i];

        ta2_yy_xxyz_xxxz_0[i] = 2.0 * ta1_y_xxz_xxxz_1[i] + ta2_yy_xxz_xxxz_0[i] * pa_y[i] - ta2_yy_xxz_xxxz_1[i] * pc_y[i];

        ta2_yy_xxyz_xxyy_0[i] = ta2_yy_xxy_xxyy_0[i] * pa_z[i] - ta2_yy_xxy_xxyy_1[i] * pc_z[i];

        ta2_yy_xxyz_xxyz_0[i] =
            ta2_yy_xxy_xxy_0[i] * fe_0 - ta2_yy_xxy_xxy_1[i] * fe_0 + ta2_yy_xxy_xxyz_0[i] * pa_z[i] - ta2_yy_xxy_xxyz_1[i] * pc_z[i];

        ta2_yy_xxyz_xxzz_0[i] = 2.0 * ta1_y_xxz_xxzz_1[i] + ta2_yy_xxz_xxzz_0[i] * pa_y[i] - ta2_yy_xxz_xxzz_1[i] * pc_y[i];

        ta2_yy_xxyz_xyyy_0[i] = ta2_yy_xxy_xyyy_0[i] * pa_z[i] - ta2_yy_xxy_xyyy_1[i] * pc_z[i];

        ta2_yy_xxyz_xyyz_0[i] =
            ta2_yy_xxy_xyy_0[i] * fe_0 - ta2_yy_xxy_xyy_1[i] * fe_0 + ta2_yy_xxy_xyyz_0[i] * pa_z[i] - ta2_yy_xxy_xyyz_1[i] * pc_z[i];

        ta2_yy_xxyz_xyzz_0[i] =
            2.0 * ta2_yy_xxy_xyz_0[i] * fe_0 - 2.0 * ta2_yy_xxy_xyz_1[i] * fe_0 + ta2_yy_xxy_xyzz_0[i] * pa_z[i] - ta2_yy_xxy_xyzz_1[i] * pc_z[i];

        ta2_yy_xxyz_xzzz_0[i] = 2.0 * ta1_y_xxz_xzzz_1[i] + ta2_yy_xxz_xzzz_0[i] * pa_y[i] - ta2_yy_xxz_xzzz_1[i] * pc_y[i];

        ta2_yy_xxyz_yyyy_0[i] = ta2_yy_xxy_yyyy_0[i] * pa_z[i] - ta2_yy_xxy_yyyy_1[i] * pc_z[i];

        ta2_yy_xxyz_yyyz_0[i] =
            ta2_yy_yz_yyyz_0[i] * fe_0 - ta2_yy_yz_yyyz_1[i] * fe_0 + ta2_yy_xyz_yyyz_0[i] * pa_x[i] - ta2_yy_xyz_yyyz_1[i] * pc_x[i];

        ta2_yy_xxyz_yyzz_0[i] =
            ta2_yy_yz_yyzz_0[i] * fe_0 - ta2_yy_yz_yyzz_1[i] * fe_0 + ta2_yy_xyz_yyzz_0[i] * pa_x[i] - ta2_yy_xyz_yyzz_1[i] * pc_x[i];

        ta2_yy_xxyz_yzzz_0[i] =
            ta2_yy_yz_yzzz_0[i] * fe_0 - ta2_yy_yz_yzzz_1[i] * fe_0 + ta2_yy_xyz_yzzz_0[i] * pa_x[i] - ta2_yy_xyz_yzzz_1[i] * pc_x[i];

        ta2_yy_xxyz_zzzz_0[i] = 2.0 * ta1_y_xxz_zzzz_1[i] + ta2_yy_xxz_zzzz_0[i] * pa_y[i] - ta2_yy_xxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 750-765 components of targeted buffer : GG

    auto ta2_yy_xxzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 750);

    auto ta2_yy_xxzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 751);

    auto ta2_yy_xxzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 752);

    auto ta2_yy_xxzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 753);

    auto ta2_yy_xxzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 754);

    auto ta2_yy_xxzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 755);

    auto ta2_yy_xxzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 756);

    auto ta2_yy_xxzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 757);

    auto ta2_yy_xxzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 758);

    auto ta2_yy_xxzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 759);

    auto ta2_yy_xxzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 760);

    auto ta2_yy_xxzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 761);

    auto ta2_yy_xxzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 762);

    auto ta2_yy_xxzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 763);

    auto ta2_yy_xxzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 764);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta2_yy_xx_xxxx_0,   \
                             ta2_yy_xx_xxxx_1,   \
                             ta2_yy_xx_xxxy_0,   \
                             ta2_yy_xx_xxxy_1,   \
                             ta2_yy_xx_xxyy_0,   \
                             ta2_yy_xx_xxyy_1,   \
                             ta2_yy_xx_xyyy_0,   \
                             ta2_yy_xx_xyyy_1,   \
                             ta2_yy_xxz_xxxx_0,  \
                             ta2_yy_xxz_xxxx_1,  \
                             ta2_yy_xxz_xxxy_0,  \
                             ta2_yy_xxz_xxxy_1,  \
                             ta2_yy_xxz_xxyy_0,  \
                             ta2_yy_xxz_xxyy_1,  \
                             ta2_yy_xxz_xyyy_0,  \
                             ta2_yy_xxz_xyyy_1,  \
                             ta2_yy_xxzz_xxxx_0, \
                             ta2_yy_xxzz_xxxy_0, \
                             ta2_yy_xxzz_xxxz_0, \
                             ta2_yy_xxzz_xxyy_0, \
                             ta2_yy_xxzz_xxyz_0, \
                             ta2_yy_xxzz_xxzz_0, \
                             ta2_yy_xxzz_xyyy_0, \
                             ta2_yy_xxzz_xyyz_0, \
                             ta2_yy_xxzz_xyzz_0, \
                             ta2_yy_xxzz_xzzz_0, \
                             ta2_yy_xxzz_yyyy_0, \
                             ta2_yy_xxzz_yyyz_0, \
                             ta2_yy_xxzz_yyzz_0, \
                             ta2_yy_xxzz_yzzz_0, \
                             ta2_yy_xxzz_zzzz_0, \
                             ta2_yy_xzz_xxxz_0,  \
                             ta2_yy_xzz_xxxz_1,  \
                             ta2_yy_xzz_xxyz_0,  \
                             ta2_yy_xzz_xxyz_1,  \
                             ta2_yy_xzz_xxz_0,   \
                             ta2_yy_xzz_xxz_1,   \
                             ta2_yy_xzz_xxzz_0,  \
                             ta2_yy_xzz_xxzz_1,  \
                             ta2_yy_xzz_xyyz_0,  \
                             ta2_yy_xzz_xyyz_1,  \
                             ta2_yy_xzz_xyz_0,   \
                             ta2_yy_xzz_xyz_1,   \
                             ta2_yy_xzz_xyzz_0,  \
                             ta2_yy_xzz_xyzz_1,  \
                             ta2_yy_xzz_xzz_0,   \
                             ta2_yy_xzz_xzz_1,   \
                             ta2_yy_xzz_xzzz_0,  \
                             ta2_yy_xzz_xzzz_1,  \
                             ta2_yy_xzz_yyyy_0,  \
                             ta2_yy_xzz_yyyy_1,  \
                             ta2_yy_xzz_yyyz_0,  \
                             ta2_yy_xzz_yyyz_1,  \
                             ta2_yy_xzz_yyz_0,   \
                             ta2_yy_xzz_yyz_1,   \
                             ta2_yy_xzz_yyzz_0,  \
                             ta2_yy_xzz_yyzz_1,  \
                             ta2_yy_xzz_yzz_0,   \
                             ta2_yy_xzz_yzz_1,   \
                             ta2_yy_xzz_yzzz_0,  \
                             ta2_yy_xzz_yzzz_1,  \
                             ta2_yy_xzz_zzz_0,   \
                             ta2_yy_xzz_zzz_1,   \
                             ta2_yy_xzz_zzzz_0,  \
                             ta2_yy_xzz_zzzz_1,  \
                             ta2_yy_zz_xxxz_0,   \
                             ta2_yy_zz_xxxz_1,   \
                             ta2_yy_zz_xxyz_0,   \
                             ta2_yy_zz_xxyz_1,   \
                             ta2_yy_zz_xxzz_0,   \
                             ta2_yy_zz_xxzz_1,   \
                             ta2_yy_zz_xyyz_0,   \
                             ta2_yy_zz_xyyz_1,   \
                             ta2_yy_zz_xyzz_0,   \
                             ta2_yy_zz_xyzz_1,   \
                             ta2_yy_zz_xzzz_0,   \
                             ta2_yy_zz_xzzz_1,   \
                             ta2_yy_zz_yyyy_0,   \
                             ta2_yy_zz_yyyy_1,   \
                             ta2_yy_zz_yyyz_0,   \
                             ta2_yy_zz_yyyz_1,   \
                             ta2_yy_zz_yyzz_0,   \
                             ta2_yy_zz_yyzz_1,   \
                             ta2_yy_zz_yzzz_0,   \
                             ta2_yy_zz_yzzz_1,   \
                             ta2_yy_zz_zzzz_0,   \
                             ta2_yy_zz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xxzz_xxxx_0[i] =
            ta2_yy_xx_xxxx_0[i] * fe_0 - ta2_yy_xx_xxxx_1[i] * fe_0 + ta2_yy_xxz_xxxx_0[i] * pa_z[i] - ta2_yy_xxz_xxxx_1[i] * pc_z[i];

        ta2_yy_xxzz_xxxy_0[i] =
            ta2_yy_xx_xxxy_0[i] * fe_0 - ta2_yy_xx_xxxy_1[i] * fe_0 + ta2_yy_xxz_xxxy_0[i] * pa_z[i] - ta2_yy_xxz_xxxy_1[i] * pc_z[i];

        ta2_yy_xxzz_xxxz_0[i] = ta2_yy_zz_xxxz_0[i] * fe_0 - ta2_yy_zz_xxxz_1[i] * fe_0 + 3.0 * ta2_yy_xzz_xxz_0[i] * fe_0 -
                                3.0 * ta2_yy_xzz_xxz_1[i] * fe_0 + ta2_yy_xzz_xxxz_0[i] * pa_x[i] - ta2_yy_xzz_xxxz_1[i] * pc_x[i];

        ta2_yy_xxzz_xxyy_0[i] =
            ta2_yy_xx_xxyy_0[i] * fe_0 - ta2_yy_xx_xxyy_1[i] * fe_0 + ta2_yy_xxz_xxyy_0[i] * pa_z[i] - ta2_yy_xxz_xxyy_1[i] * pc_z[i];

        ta2_yy_xxzz_xxyz_0[i] = ta2_yy_zz_xxyz_0[i] * fe_0 - ta2_yy_zz_xxyz_1[i] * fe_0 + 2.0 * ta2_yy_xzz_xyz_0[i] * fe_0 -
                                2.0 * ta2_yy_xzz_xyz_1[i] * fe_0 + ta2_yy_xzz_xxyz_0[i] * pa_x[i] - ta2_yy_xzz_xxyz_1[i] * pc_x[i];

        ta2_yy_xxzz_xxzz_0[i] = ta2_yy_zz_xxzz_0[i] * fe_0 - ta2_yy_zz_xxzz_1[i] * fe_0 + 2.0 * ta2_yy_xzz_xzz_0[i] * fe_0 -
                                2.0 * ta2_yy_xzz_xzz_1[i] * fe_0 + ta2_yy_xzz_xxzz_0[i] * pa_x[i] - ta2_yy_xzz_xxzz_1[i] * pc_x[i];

        ta2_yy_xxzz_xyyy_0[i] =
            ta2_yy_xx_xyyy_0[i] * fe_0 - ta2_yy_xx_xyyy_1[i] * fe_0 + ta2_yy_xxz_xyyy_0[i] * pa_z[i] - ta2_yy_xxz_xyyy_1[i] * pc_z[i];

        ta2_yy_xxzz_xyyz_0[i] = ta2_yy_zz_xyyz_0[i] * fe_0 - ta2_yy_zz_xyyz_1[i] * fe_0 + ta2_yy_xzz_yyz_0[i] * fe_0 - ta2_yy_xzz_yyz_1[i] * fe_0 +
                                ta2_yy_xzz_xyyz_0[i] * pa_x[i] - ta2_yy_xzz_xyyz_1[i] * pc_x[i];

        ta2_yy_xxzz_xyzz_0[i] = ta2_yy_zz_xyzz_0[i] * fe_0 - ta2_yy_zz_xyzz_1[i] * fe_0 + ta2_yy_xzz_yzz_0[i] * fe_0 - ta2_yy_xzz_yzz_1[i] * fe_0 +
                                ta2_yy_xzz_xyzz_0[i] * pa_x[i] - ta2_yy_xzz_xyzz_1[i] * pc_x[i];

        ta2_yy_xxzz_xzzz_0[i] = ta2_yy_zz_xzzz_0[i] * fe_0 - ta2_yy_zz_xzzz_1[i] * fe_0 + ta2_yy_xzz_zzz_0[i] * fe_0 - ta2_yy_xzz_zzz_1[i] * fe_0 +
                                ta2_yy_xzz_xzzz_0[i] * pa_x[i] - ta2_yy_xzz_xzzz_1[i] * pc_x[i];

        ta2_yy_xxzz_yyyy_0[i] =
            ta2_yy_zz_yyyy_0[i] * fe_0 - ta2_yy_zz_yyyy_1[i] * fe_0 + ta2_yy_xzz_yyyy_0[i] * pa_x[i] - ta2_yy_xzz_yyyy_1[i] * pc_x[i];

        ta2_yy_xxzz_yyyz_0[i] =
            ta2_yy_zz_yyyz_0[i] * fe_0 - ta2_yy_zz_yyyz_1[i] * fe_0 + ta2_yy_xzz_yyyz_0[i] * pa_x[i] - ta2_yy_xzz_yyyz_1[i] * pc_x[i];

        ta2_yy_xxzz_yyzz_0[i] =
            ta2_yy_zz_yyzz_0[i] * fe_0 - ta2_yy_zz_yyzz_1[i] * fe_0 + ta2_yy_xzz_yyzz_0[i] * pa_x[i] - ta2_yy_xzz_yyzz_1[i] * pc_x[i];

        ta2_yy_xxzz_yzzz_0[i] =
            ta2_yy_zz_yzzz_0[i] * fe_0 - ta2_yy_zz_yzzz_1[i] * fe_0 + ta2_yy_xzz_yzzz_0[i] * pa_x[i] - ta2_yy_xzz_yzzz_1[i] * pc_x[i];

        ta2_yy_xxzz_zzzz_0[i] =
            ta2_yy_zz_zzzz_0[i] * fe_0 - ta2_yy_zz_zzzz_1[i] * fe_0 + ta2_yy_xzz_zzzz_0[i] * pa_x[i] - ta2_yy_xzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 765-780 components of targeted buffer : GG

    auto ta2_yy_xyyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 765);

    auto ta2_yy_xyyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 766);

    auto ta2_yy_xyyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 767);

    auto ta2_yy_xyyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 768);

    auto ta2_yy_xyyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 769);

    auto ta2_yy_xyyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 770);

    auto ta2_yy_xyyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 771);

    auto ta2_yy_xyyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 772);

    auto ta2_yy_xyyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 773);

    auto ta2_yy_xyyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 774);

    auto ta2_yy_xyyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 775);

    auto ta2_yy_xyyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 776);

    auto ta2_yy_xyyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 777);

    auto ta2_yy_xyyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 778);

    auto ta2_yy_xyyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 779);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta2_yy_xyyy_xxxx_0, \
                             ta2_yy_xyyy_xxxy_0, \
                             ta2_yy_xyyy_xxxz_0, \
                             ta2_yy_xyyy_xxyy_0, \
                             ta2_yy_xyyy_xxyz_0, \
                             ta2_yy_xyyy_xxzz_0, \
                             ta2_yy_xyyy_xyyy_0, \
                             ta2_yy_xyyy_xyyz_0, \
                             ta2_yy_xyyy_xyzz_0, \
                             ta2_yy_xyyy_xzzz_0, \
                             ta2_yy_xyyy_yyyy_0, \
                             ta2_yy_xyyy_yyyz_0, \
                             ta2_yy_xyyy_yyzz_0, \
                             ta2_yy_xyyy_yzzz_0, \
                             ta2_yy_xyyy_zzzz_0, \
                             ta2_yy_yyy_xxx_0,   \
                             ta2_yy_yyy_xxx_1,   \
                             ta2_yy_yyy_xxxx_0,  \
                             ta2_yy_yyy_xxxx_1,  \
                             ta2_yy_yyy_xxxy_0,  \
                             ta2_yy_yyy_xxxy_1,  \
                             ta2_yy_yyy_xxxz_0,  \
                             ta2_yy_yyy_xxxz_1,  \
                             ta2_yy_yyy_xxy_0,   \
                             ta2_yy_yyy_xxy_1,   \
                             ta2_yy_yyy_xxyy_0,  \
                             ta2_yy_yyy_xxyy_1,  \
                             ta2_yy_yyy_xxyz_0,  \
                             ta2_yy_yyy_xxyz_1,  \
                             ta2_yy_yyy_xxz_0,   \
                             ta2_yy_yyy_xxz_1,   \
                             ta2_yy_yyy_xxzz_0,  \
                             ta2_yy_yyy_xxzz_1,  \
                             ta2_yy_yyy_xyy_0,   \
                             ta2_yy_yyy_xyy_1,   \
                             ta2_yy_yyy_xyyy_0,  \
                             ta2_yy_yyy_xyyy_1,  \
                             ta2_yy_yyy_xyyz_0,  \
                             ta2_yy_yyy_xyyz_1,  \
                             ta2_yy_yyy_xyz_0,   \
                             ta2_yy_yyy_xyz_1,   \
                             ta2_yy_yyy_xyzz_0,  \
                             ta2_yy_yyy_xyzz_1,  \
                             ta2_yy_yyy_xzz_0,   \
                             ta2_yy_yyy_xzz_1,   \
                             ta2_yy_yyy_xzzz_0,  \
                             ta2_yy_yyy_xzzz_1,  \
                             ta2_yy_yyy_yyy_0,   \
                             ta2_yy_yyy_yyy_1,   \
                             ta2_yy_yyy_yyyy_0,  \
                             ta2_yy_yyy_yyyy_1,  \
                             ta2_yy_yyy_yyyz_0,  \
                             ta2_yy_yyy_yyyz_1,  \
                             ta2_yy_yyy_yyz_0,   \
                             ta2_yy_yyy_yyz_1,   \
                             ta2_yy_yyy_yyzz_0,  \
                             ta2_yy_yyy_yyzz_1,  \
                             ta2_yy_yyy_yzz_0,   \
                             ta2_yy_yyy_yzz_1,   \
                             ta2_yy_yyy_yzzz_0,  \
                             ta2_yy_yyy_yzzz_1,  \
                             ta2_yy_yyy_zzz_0,   \
                             ta2_yy_yyy_zzz_1,   \
                             ta2_yy_yyy_zzzz_0,  \
                             ta2_yy_yyy_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xyyy_xxxx_0[i] =
            4.0 * ta2_yy_yyy_xxx_0[i] * fe_0 - 4.0 * ta2_yy_yyy_xxx_1[i] * fe_0 + ta2_yy_yyy_xxxx_0[i] * pa_x[i] - ta2_yy_yyy_xxxx_1[i] * pc_x[i];

        ta2_yy_xyyy_xxxy_0[i] =
            3.0 * ta2_yy_yyy_xxy_0[i] * fe_0 - 3.0 * ta2_yy_yyy_xxy_1[i] * fe_0 + ta2_yy_yyy_xxxy_0[i] * pa_x[i] - ta2_yy_yyy_xxxy_1[i] * pc_x[i];

        ta2_yy_xyyy_xxxz_0[i] =
            3.0 * ta2_yy_yyy_xxz_0[i] * fe_0 - 3.0 * ta2_yy_yyy_xxz_1[i] * fe_0 + ta2_yy_yyy_xxxz_0[i] * pa_x[i] - ta2_yy_yyy_xxxz_1[i] * pc_x[i];

        ta2_yy_xyyy_xxyy_0[i] =
            2.0 * ta2_yy_yyy_xyy_0[i] * fe_0 - 2.0 * ta2_yy_yyy_xyy_1[i] * fe_0 + ta2_yy_yyy_xxyy_0[i] * pa_x[i] - ta2_yy_yyy_xxyy_1[i] * pc_x[i];

        ta2_yy_xyyy_xxyz_0[i] =
            2.0 * ta2_yy_yyy_xyz_0[i] * fe_0 - 2.0 * ta2_yy_yyy_xyz_1[i] * fe_0 + ta2_yy_yyy_xxyz_0[i] * pa_x[i] - ta2_yy_yyy_xxyz_1[i] * pc_x[i];

        ta2_yy_xyyy_xxzz_0[i] =
            2.0 * ta2_yy_yyy_xzz_0[i] * fe_0 - 2.0 * ta2_yy_yyy_xzz_1[i] * fe_0 + ta2_yy_yyy_xxzz_0[i] * pa_x[i] - ta2_yy_yyy_xxzz_1[i] * pc_x[i];

        ta2_yy_xyyy_xyyy_0[i] =
            ta2_yy_yyy_yyy_0[i] * fe_0 - ta2_yy_yyy_yyy_1[i] * fe_0 + ta2_yy_yyy_xyyy_0[i] * pa_x[i] - ta2_yy_yyy_xyyy_1[i] * pc_x[i];

        ta2_yy_xyyy_xyyz_0[i] =
            ta2_yy_yyy_yyz_0[i] * fe_0 - ta2_yy_yyy_yyz_1[i] * fe_0 + ta2_yy_yyy_xyyz_0[i] * pa_x[i] - ta2_yy_yyy_xyyz_1[i] * pc_x[i];

        ta2_yy_xyyy_xyzz_0[i] =
            ta2_yy_yyy_yzz_0[i] * fe_0 - ta2_yy_yyy_yzz_1[i] * fe_0 + ta2_yy_yyy_xyzz_0[i] * pa_x[i] - ta2_yy_yyy_xyzz_1[i] * pc_x[i];

        ta2_yy_xyyy_xzzz_0[i] =
            ta2_yy_yyy_zzz_0[i] * fe_0 - ta2_yy_yyy_zzz_1[i] * fe_0 + ta2_yy_yyy_xzzz_0[i] * pa_x[i] - ta2_yy_yyy_xzzz_1[i] * pc_x[i];

        ta2_yy_xyyy_yyyy_0[i] = ta2_yy_yyy_yyyy_0[i] * pa_x[i] - ta2_yy_yyy_yyyy_1[i] * pc_x[i];

        ta2_yy_xyyy_yyyz_0[i] = ta2_yy_yyy_yyyz_0[i] * pa_x[i] - ta2_yy_yyy_yyyz_1[i] * pc_x[i];

        ta2_yy_xyyy_yyzz_0[i] = ta2_yy_yyy_yyzz_0[i] * pa_x[i] - ta2_yy_yyy_yyzz_1[i] * pc_x[i];

        ta2_yy_xyyy_yzzz_0[i] = ta2_yy_yyy_yzzz_0[i] * pa_x[i] - ta2_yy_yyy_yzzz_1[i] * pc_x[i];

        ta2_yy_xyyy_zzzz_0[i] = ta2_yy_yyy_zzzz_0[i] * pa_x[i] - ta2_yy_yyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 780-795 components of targeted buffer : GG

    auto ta2_yy_xyyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 780);

    auto ta2_yy_xyyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 781);

    auto ta2_yy_xyyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 782);

    auto ta2_yy_xyyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 783);

    auto ta2_yy_xyyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 784);

    auto ta2_yy_xyyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 785);

    auto ta2_yy_xyyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 786);

    auto ta2_yy_xyyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 787);

    auto ta2_yy_xyyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 788);

    auto ta2_yy_xyyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 789);

    auto ta2_yy_xyyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 790);

    auto ta2_yy_xyyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 791);

    auto ta2_yy_xyyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 792);

    auto ta2_yy_xyyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 793);

    auto ta2_yy_xyyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 794);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta2_yy_xyy_xxxx_0,  \
                             ta2_yy_xyy_xxxx_1,  \
                             ta2_yy_xyy_xxxy_0,  \
                             ta2_yy_xyy_xxxy_1,  \
                             ta2_yy_xyy_xxyy_0,  \
                             ta2_yy_xyy_xxyy_1,  \
                             ta2_yy_xyy_xyyy_0,  \
                             ta2_yy_xyy_xyyy_1,  \
                             ta2_yy_xyyz_xxxx_0, \
                             ta2_yy_xyyz_xxxy_0, \
                             ta2_yy_xyyz_xxxz_0, \
                             ta2_yy_xyyz_xxyy_0, \
                             ta2_yy_xyyz_xxyz_0, \
                             ta2_yy_xyyz_xxzz_0, \
                             ta2_yy_xyyz_xyyy_0, \
                             ta2_yy_xyyz_xyyz_0, \
                             ta2_yy_xyyz_xyzz_0, \
                             ta2_yy_xyyz_xzzz_0, \
                             ta2_yy_xyyz_yyyy_0, \
                             ta2_yy_xyyz_yyyz_0, \
                             ta2_yy_xyyz_yyzz_0, \
                             ta2_yy_xyyz_yzzz_0, \
                             ta2_yy_xyyz_zzzz_0, \
                             ta2_yy_yyz_xxxz_0,  \
                             ta2_yy_yyz_xxxz_1,  \
                             ta2_yy_yyz_xxyz_0,  \
                             ta2_yy_yyz_xxyz_1,  \
                             ta2_yy_yyz_xxz_0,   \
                             ta2_yy_yyz_xxz_1,   \
                             ta2_yy_yyz_xxzz_0,  \
                             ta2_yy_yyz_xxzz_1,  \
                             ta2_yy_yyz_xyyz_0,  \
                             ta2_yy_yyz_xyyz_1,  \
                             ta2_yy_yyz_xyz_0,   \
                             ta2_yy_yyz_xyz_1,   \
                             ta2_yy_yyz_xyzz_0,  \
                             ta2_yy_yyz_xyzz_1,  \
                             ta2_yy_yyz_xzz_0,   \
                             ta2_yy_yyz_xzz_1,   \
                             ta2_yy_yyz_xzzz_0,  \
                             ta2_yy_yyz_xzzz_1,  \
                             ta2_yy_yyz_yyyy_0,  \
                             ta2_yy_yyz_yyyy_1,  \
                             ta2_yy_yyz_yyyz_0,  \
                             ta2_yy_yyz_yyyz_1,  \
                             ta2_yy_yyz_yyz_0,   \
                             ta2_yy_yyz_yyz_1,   \
                             ta2_yy_yyz_yyzz_0,  \
                             ta2_yy_yyz_yyzz_1,  \
                             ta2_yy_yyz_yzz_0,   \
                             ta2_yy_yyz_yzz_1,   \
                             ta2_yy_yyz_yzzz_0,  \
                             ta2_yy_yyz_yzzz_1,  \
                             ta2_yy_yyz_zzz_0,   \
                             ta2_yy_yyz_zzz_1,   \
                             ta2_yy_yyz_zzzz_0,  \
                             ta2_yy_yyz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xyyz_xxxx_0[i] = ta2_yy_xyy_xxxx_0[i] * pa_z[i] - ta2_yy_xyy_xxxx_1[i] * pc_z[i];

        ta2_yy_xyyz_xxxy_0[i] = ta2_yy_xyy_xxxy_0[i] * pa_z[i] - ta2_yy_xyy_xxxy_1[i] * pc_z[i];

        ta2_yy_xyyz_xxxz_0[i] =
            3.0 * ta2_yy_yyz_xxz_0[i] * fe_0 - 3.0 * ta2_yy_yyz_xxz_1[i] * fe_0 + ta2_yy_yyz_xxxz_0[i] * pa_x[i] - ta2_yy_yyz_xxxz_1[i] * pc_x[i];

        ta2_yy_xyyz_xxyy_0[i] = ta2_yy_xyy_xxyy_0[i] * pa_z[i] - ta2_yy_xyy_xxyy_1[i] * pc_z[i];

        ta2_yy_xyyz_xxyz_0[i] =
            2.0 * ta2_yy_yyz_xyz_0[i] * fe_0 - 2.0 * ta2_yy_yyz_xyz_1[i] * fe_0 + ta2_yy_yyz_xxyz_0[i] * pa_x[i] - ta2_yy_yyz_xxyz_1[i] * pc_x[i];

        ta2_yy_xyyz_xxzz_0[i] =
            2.0 * ta2_yy_yyz_xzz_0[i] * fe_0 - 2.0 * ta2_yy_yyz_xzz_1[i] * fe_0 + ta2_yy_yyz_xxzz_0[i] * pa_x[i] - ta2_yy_yyz_xxzz_1[i] * pc_x[i];

        ta2_yy_xyyz_xyyy_0[i] = ta2_yy_xyy_xyyy_0[i] * pa_z[i] - ta2_yy_xyy_xyyy_1[i] * pc_z[i];

        ta2_yy_xyyz_xyyz_0[i] =
            ta2_yy_yyz_yyz_0[i] * fe_0 - ta2_yy_yyz_yyz_1[i] * fe_0 + ta2_yy_yyz_xyyz_0[i] * pa_x[i] - ta2_yy_yyz_xyyz_1[i] * pc_x[i];

        ta2_yy_xyyz_xyzz_0[i] =
            ta2_yy_yyz_yzz_0[i] * fe_0 - ta2_yy_yyz_yzz_1[i] * fe_0 + ta2_yy_yyz_xyzz_0[i] * pa_x[i] - ta2_yy_yyz_xyzz_1[i] * pc_x[i];

        ta2_yy_xyyz_xzzz_0[i] =
            ta2_yy_yyz_zzz_0[i] * fe_0 - ta2_yy_yyz_zzz_1[i] * fe_0 + ta2_yy_yyz_xzzz_0[i] * pa_x[i] - ta2_yy_yyz_xzzz_1[i] * pc_x[i];

        ta2_yy_xyyz_yyyy_0[i] = ta2_yy_yyz_yyyy_0[i] * pa_x[i] - ta2_yy_yyz_yyyy_1[i] * pc_x[i];

        ta2_yy_xyyz_yyyz_0[i] = ta2_yy_yyz_yyyz_0[i] * pa_x[i] - ta2_yy_yyz_yyyz_1[i] * pc_x[i];

        ta2_yy_xyyz_yyzz_0[i] = ta2_yy_yyz_yyzz_0[i] * pa_x[i] - ta2_yy_yyz_yyzz_1[i] * pc_x[i];

        ta2_yy_xyyz_yzzz_0[i] = ta2_yy_yyz_yzzz_0[i] * pa_x[i] - ta2_yy_yyz_yzzz_1[i] * pc_x[i];

        ta2_yy_xyyz_zzzz_0[i] = ta2_yy_yyz_zzzz_0[i] * pa_x[i] - ta2_yy_yyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 795-810 components of targeted buffer : GG

    auto ta2_yy_xyzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 795);

    auto ta2_yy_xyzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 796);

    auto ta2_yy_xyzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 797);

    auto ta2_yy_xyzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 798);

    auto ta2_yy_xyzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 799);

    auto ta2_yy_xyzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 800);

    auto ta2_yy_xyzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 801);

    auto ta2_yy_xyzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 802);

    auto ta2_yy_xyzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 803);

    auto ta2_yy_xyzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 804);

    auto ta2_yy_xyzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 805);

    auto ta2_yy_xyzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 806);

    auto ta2_yy_xyzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 807);

    auto ta2_yy_xyzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 808);

    auto ta2_yy_xyzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 809);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_y_xzz_xxxx_1,   \
                             ta1_y_xzz_xxxz_1,   \
                             ta1_y_xzz_xxzz_1,   \
                             ta1_y_xzz_xzzz_1,   \
                             ta2_yy_xyzz_xxxx_0, \
                             ta2_yy_xyzz_xxxy_0, \
                             ta2_yy_xyzz_xxxz_0, \
                             ta2_yy_xyzz_xxyy_0, \
                             ta2_yy_xyzz_xxyz_0, \
                             ta2_yy_xyzz_xxzz_0, \
                             ta2_yy_xyzz_xyyy_0, \
                             ta2_yy_xyzz_xyyz_0, \
                             ta2_yy_xyzz_xyzz_0, \
                             ta2_yy_xyzz_xzzz_0, \
                             ta2_yy_xyzz_yyyy_0, \
                             ta2_yy_xyzz_yyyz_0, \
                             ta2_yy_xyzz_yyzz_0, \
                             ta2_yy_xyzz_yzzz_0, \
                             ta2_yy_xyzz_zzzz_0, \
                             ta2_yy_xzz_xxxx_0,  \
                             ta2_yy_xzz_xxxx_1,  \
                             ta2_yy_xzz_xxxz_0,  \
                             ta2_yy_xzz_xxxz_1,  \
                             ta2_yy_xzz_xxzz_0,  \
                             ta2_yy_xzz_xxzz_1,  \
                             ta2_yy_xzz_xzzz_0,  \
                             ta2_yy_xzz_xzzz_1,  \
                             ta2_yy_yzz_xxxy_0,  \
                             ta2_yy_yzz_xxxy_1,  \
                             ta2_yy_yzz_xxy_0,   \
                             ta2_yy_yzz_xxy_1,   \
                             ta2_yy_yzz_xxyy_0,  \
                             ta2_yy_yzz_xxyy_1,  \
                             ta2_yy_yzz_xxyz_0,  \
                             ta2_yy_yzz_xxyz_1,  \
                             ta2_yy_yzz_xyy_0,   \
                             ta2_yy_yzz_xyy_1,   \
                             ta2_yy_yzz_xyyy_0,  \
                             ta2_yy_yzz_xyyy_1,  \
                             ta2_yy_yzz_xyyz_0,  \
                             ta2_yy_yzz_xyyz_1,  \
                             ta2_yy_yzz_xyz_0,   \
                             ta2_yy_yzz_xyz_1,   \
                             ta2_yy_yzz_xyzz_0,  \
                             ta2_yy_yzz_xyzz_1,  \
                             ta2_yy_yzz_yyy_0,   \
                             ta2_yy_yzz_yyy_1,   \
                             ta2_yy_yzz_yyyy_0,  \
                             ta2_yy_yzz_yyyy_1,  \
                             ta2_yy_yzz_yyyz_0,  \
                             ta2_yy_yzz_yyyz_1,  \
                             ta2_yy_yzz_yyz_0,   \
                             ta2_yy_yzz_yyz_1,   \
                             ta2_yy_yzz_yyzz_0,  \
                             ta2_yy_yzz_yyzz_1,  \
                             ta2_yy_yzz_yzz_0,   \
                             ta2_yy_yzz_yzz_1,   \
                             ta2_yy_yzz_yzzz_0,  \
                             ta2_yy_yzz_yzzz_1,  \
                             ta2_yy_yzz_zzzz_0,  \
                             ta2_yy_yzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xyzz_xxxx_0[i] = 2.0 * ta1_y_xzz_xxxx_1[i] + ta2_yy_xzz_xxxx_0[i] * pa_y[i] - ta2_yy_xzz_xxxx_1[i] * pc_y[i];

        ta2_yy_xyzz_xxxy_0[i] =
            3.0 * ta2_yy_yzz_xxy_0[i] * fe_0 - 3.0 * ta2_yy_yzz_xxy_1[i] * fe_0 + ta2_yy_yzz_xxxy_0[i] * pa_x[i] - ta2_yy_yzz_xxxy_1[i] * pc_x[i];

        ta2_yy_xyzz_xxxz_0[i] = 2.0 * ta1_y_xzz_xxxz_1[i] + ta2_yy_xzz_xxxz_0[i] * pa_y[i] - ta2_yy_xzz_xxxz_1[i] * pc_y[i];

        ta2_yy_xyzz_xxyy_0[i] =
            2.0 * ta2_yy_yzz_xyy_0[i] * fe_0 - 2.0 * ta2_yy_yzz_xyy_1[i] * fe_0 + ta2_yy_yzz_xxyy_0[i] * pa_x[i] - ta2_yy_yzz_xxyy_1[i] * pc_x[i];

        ta2_yy_xyzz_xxyz_0[i] =
            2.0 * ta2_yy_yzz_xyz_0[i] * fe_0 - 2.0 * ta2_yy_yzz_xyz_1[i] * fe_0 + ta2_yy_yzz_xxyz_0[i] * pa_x[i] - ta2_yy_yzz_xxyz_1[i] * pc_x[i];

        ta2_yy_xyzz_xxzz_0[i] = 2.0 * ta1_y_xzz_xxzz_1[i] + ta2_yy_xzz_xxzz_0[i] * pa_y[i] - ta2_yy_xzz_xxzz_1[i] * pc_y[i];

        ta2_yy_xyzz_xyyy_0[i] =
            ta2_yy_yzz_yyy_0[i] * fe_0 - ta2_yy_yzz_yyy_1[i] * fe_0 + ta2_yy_yzz_xyyy_0[i] * pa_x[i] - ta2_yy_yzz_xyyy_1[i] * pc_x[i];

        ta2_yy_xyzz_xyyz_0[i] =
            ta2_yy_yzz_yyz_0[i] * fe_0 - ta2_yy_yzz_yyz_1[i] * fe_0 + ta2_yy_yzz_xyyz_0[i] * pa_x[i] - ta2_yy_yzz_xyyz_1[i] * pc_x[i];

        ta2_yy_xyzz_xyzz_0[i] =
            ta2_yy_yzz_yzz_0[i] * fe_0 - ta2_yy_yzz_yzz_1[i] * fe_0 + ta2_yy_yzz_xyzz_0[i] * pa_x[i] - ta2_yy_yzz_xyzz_1[i] * pc_x[i];

        ta2_yy_xyzz_xzzz_0[i] = 2.0 * ta1_y_xzz_xzzz_1[i] + ta2_yy_xzz_xzzz_0[i] * pa_y[i] - ta2_yy_xzz_xzzz_1[i] * pc_y[i];

        ta2_yy_xyzz_yyyy_0[i] = ta2_yy_yzz_yyyy_0[i] * pa_x[i] - ta2_yy_yzz_yyyy_1[i] * pc_x[i];

        ta2_yy_xyzz_yyyz_0[i] = ta2_yy_yzz_yyyz_0[i] * pa_x[i] - ta2_yy_yzz_yyyz_1[i] * pc_x[i];

        ta2_yy_xyzz_yyzz_0[i] = ta2_yy_yzz_yyzz_0[i] * pa_x[i] - ta2_yy_yzz_yyzz_1[i] * pc_x[i];

        ta2_yy_xyzz_yzzz_0[i] = ta2_yy_yzz_yzzz_0[i] * pa_x[i] - ta2_yy_yzz_yzzz_1[i] * pc_x[i];

        ta2_yy_xyzz_zzzz_0[i] = ta2_yy_yzz_zzzz_0[i] * pa_x[i] - ta2_yy_yzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 810-825 components of targeted buffer : GG

    auto ta2_yy_xzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 810);

    auto ta2_yy_xzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 811);

    auto ta2_yy_xzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 812);

    auto ta2_yy_xzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 813);

    auto ta2_yy_xzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 814);

    auto ta2_yy_xzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 815);

    auto ta2_yy_xzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 816);

    auto ta2_yy_xzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 817);

    auto ta2_yy_xzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 818);

    auto ta2_yy_xzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 819);

    auto ta2_yy_xzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 820);

    auto ta2_yy_xzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 821);

    auto ta2_yy_xzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 822);

    auto ta2_yy_xzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 823);

    auto ta2_yy_xzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 824);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta2_yy_xzzz_xxxx_0, \
                             ta2_yy_xzzz_xxxy_0, \
                             ta2_yy_xzzz_xxxz_0, \
                             ta2_yy_xzzz_xxyy_0, \
                             ta2_yy_xzzz_xxyz_0, \
                             ta2_yy_xzzz_xxzz_0, \
                             ta2_yy_xzzz_xyyy_0, \
                             ta2_yy_xzzz_xyyz_0, \
                             ta2_yy_xzzz_xyzz_0, \
                             ta2_yy_xzzz_xzzz_0, \
                             ta2_yy_xzzz_yyyy_0, \
                             ta2_yy_xzzz_yyyz_0, \
                             ta2_yy_xzzz_yyzz_0, \
                             ta2_yy_xzzz_yzzz_0, \
                             ta2_yy_xzzz_zzzz_0, \
                             ta2_yy_zzz_xxx_0,   \
                             ta2_yy_zzz_xxx_1,   \
                             ta2_yy_zzz_xxxx_0,  \
                             ta2_yy_zzz_xxxx_1,  \
                             ta2_yy_zzz_xxxy_0,  \
                             ta2_yy_zzz_xxxy_1,  \
                             ta2_yy_zzz_xxxz_0,  \
                             ta2_yy_zzz_xxxz_1,  \
                             ta2_yy_zzz_xxy_0,   \
                             ta2_yy_zzz_xxy_1,   \
                             ta2_yy_zzz_xxyy_0,  \
                             ta2_yy_zzz_xxyy_1,  \
                             ta2_yy_zzz_xxyz_0,  \
                             ta2_yy_zzz_xxyz_1,  \
                             ta2_yy_zzz_xxz_0,   \
                             ta2_yy_zzz_xxz_1,   \
                             ta2_yy_zzz_xxzz_0,  \
                             ta2_yy_zzz_xxzz_1,  \
                             ta2_yy_zzz_xyy_0,   \
                             ta2_yy_zzz_xyy_1,   \
                             ta2_yy_zzz_xyyy_0,  \
                             ta2_yy_zzz_xyyy_1,  \
                             ta2_yy_zzz_xyyz_0,  \
                             ta2_yy_zzz_xyyz_1,  \
                             ta2_yy_zzz_xyz_0,   \
                             ta2_yy_zzz_xyz_1,   \
                             ta2_yy_zzz_xyzz_0,  \
                             ta2_yy_zzz_xyzz_1,  \
                             ta2_yy_zzz_xzz_0,   \
                             ta2_yy_zzz_xzz_1,   \
                             ta2_yy_zzz_xzzz_0,  \
                             ta2_yy_zzz_xzzz_1,  \
                             ta2_yy_zzz_yyy_0,   \
                             ta2_yy_zzz_yyy_1,   \
                             ta2_yy_zzz_yyyy_0,  \
                             ta2_yy_zzz_yyyy_1,  \
                             ta2_yy_zzz_yyyz_0,  \
                             ta2_yy_zzz_yyyz_1,  \
                             ta2_yy_zzz_yyz_0,   \
                             ta2_yy_zzz_yyz_1,   \
                             ta2_yy_zzz_yyzz_0,  \
                             ta2_yy_zzz_yyzz_1,  \
                             ta2_yy_zzz_yzz_0,   \
                             ta2_yy_zzz_yzz_1,   \
                             ta2_yy_zzz_yzzz_0,  \
                             ta2_yy_zzz_yzzz_1,  \
                             ta2_yy_zzz_zzz_0,   \
                             ta2_yy_zzz_zzz_1,   \
                             ta2_yy_zzz_zzzz_0,  \
                             ta2_yy_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_xzzz_xxxx_0[i] =
            4.0 * ta2_yy_zzz_xxx_0[i] * fe_0 - 4.0 * ta2_yy_zzz_xxx_1[i] * fe_0 + ta2_yy_zzz_xxxx_0[i] * pa_x[i] - ta2_yy_zzz_xxxx_1[i] * pc_x[i];

        ta2_yy_xzzz_xxxy_0[i] =
            3.0 * ta2_yy_zzz_xxy_0[i] * fe_0 - 3.0 * ta2_yy_zzz_xxy_1[i] * fe_0 + ta2_yy_zzz_xxxy_0[i] * pa_x[i] - ta2_yy_zzz_xxxy_1[i] * pc_x[i];

        ta2_yy_xzzz_xxxz_0[i] =
            3.0 * ta2_yy_zzz_xxz_0[i] * fe_0 - 3.0 * ta2_yy_zzz_xxz_1[i] * fe_0 + ta2_yy_zzz_xxxz_0[i] * pa_x[i] - ta2_yy_zzz_xxxz_1[i] * pc_x[i];

        ta2_yy_xzzz_xxyy_0[i] =
            2.0 * ta2_yy_zzz_xyy_0[i] * fe_0 - 2.0 * ta2_yy_zzz_xyy_1[i] * fe_0 + ta2_yy_zzz_xxyy_0[i] * pa_x[i] - ta2_yy_zzz_xxyy_1[i] * pc_x[i];

        ta2_yy_xzzz_xxyz_0[i] =
            2.0 * ta2_yy_zzz_xyz_0[i] * fe_0 - 2.0 * ta2_yy_zzz_xyz_1[i] * fe_0 + ta2_yy_zzz_xxyz_0[i] * pa_x[i] - ta2_yy_zzz_xxyz_1[i] * pc_x[i];

        ta2_yy_xzzz_xxzz_0[i] =
            2.0 * ta2_yy_zzz_xzz_0[i] * fe_0 - 2.0 * ta2_yy_zzz_xzz_1[i] * fe_0 + ta2_yy_zzz_xxzz_0[i] * pa_x[i] - ta2_yy_zzz_xxzz_1[i] * pc_x[i];

        ta2_yy_xzzz_xyyy_0[i] =
            ta2_yy_zzz_yyy_0[i] * fe_0 - ta2_yy_zzz_yyy_1[i] * fe_0 + ta2_yy_zzz_xyyy_0[i] * pa_x[i] - ta2_yy_zzz_xyyy_1[i] * pc_x[i];

        ta2_yy_xzzz_xyyz_0[i] =
            ta2_yy_zzz_yyz_0[i] * fe_0 - ta2_yy_zzz_yyz_1[i] * fe_0 + ta2_yy_zzz_xyyz_0[i] * pa_x[i] - ta2_yy_zzz_xyyz_1[i] * pc_x[i];

        ta2_yy_xzzz_xyzz_0[i] =
            ta2_yy_zzz_yzz_0[i] * fe_0 - ta2_yy_zzz_yzz_1[i] * fe_0 + ta2_yy_zzz_xyzz_0[i] * pa_x[i] - ta2_yy_zzz_xyzz_1[i] * pc_x[i];

        ta2_yy_xzzz_xzzz_0[i] =
            ta2_yy_zzz_zzz_0[i] * fe_0 - ta2_yy_zzz_zzz_1[i] * fe_0 + ta2_yy_zzz_xzzz_0[i] * pa_x[i] - ta2_yy_zzz_xzzz_1[i] * pc_x[i];

        ta2_yy_xzzz_yyyy_0[i] = ta2_yy_zzz_yyyy_0[i] * pa_x[i] - ta2_yy_zzz_yyyy_1[i] * pc_x[i];

        ta2_yy_xzzz_yyyz_0[i] = ta2_yy_zzz_yyyz_0[i] * pa_x[i] - ta2_yy_zzz_yyyz_1[i] * pc_x[i];

        ta2_yy_xzzz_yyzz_0[i] = ta2_yy_zzz_yyzz_0[i] * pa_x[i] - ta2_yy_zzz_yyzz_1[i] * pc_x[i];

        ta2_yy_xzzz_yzzz_0[i] = ta2_yy_zzz_yzzz_0[i] * pa_x[i] - ta2_yy_zzz_yzzz_1[i] * pc_x[i];

        ta2_yy_xzzz_zzzz_0[i] = ta2_yy_zzz_zzzz_0[i] * pa_x[i] - ta2_yy_zzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 825-840 components of targeted buffer : GG

    auto ta2_yy_yyyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 825);

    auto ta2_yy_yyyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 826);

    auto ta2_yy_yyyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 827);

    auto ta2_yy_yyyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 828);

    auto ta2_yy_yyyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 829);

    auto ta2_yy_yyyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 830);

    auto ta2_yy_yyyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 831);

    auto ta2_yy_yyyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 832);

    auto ta2_yy_yyyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 833);

    auto ta2_yy_yyyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 834);

    auto ta2_yy_yyyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 835);

    auto ta2_yy_yyyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 836);

    auto ta2_yy_yyyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 837);

    auto ta2_yy_yyyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 838);

    auto ta2_yy_yyyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 839);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_y_yyy_xxxx_1,   \
                             ta1_y_yyy_xxxy_1,   \
                             ta1_y_yyy_xxxz_1,   \
                             ta1_y_yyy_xxyy_1,   \
                             ta1_y_yyy_xxyz_1,   \
                             ta1_y_yyy_xxzz_1,   \
                             ta1_y_yyy_xyyy_1,   \
                             ta1_y_yyy_xyyz_1,   \
                             ta1_y_yyy_xyzz_1,   \
                             ta1_y_yyy_xzzz_1,   \
                             ta1_y_yyy_yyyy_1,   \
                             ta1_y_yyy_yyyz_1,   \
                             ta1_y_yyy_yyzz_1,   \
                             ta1_y_yyy_yzzz_1,   \
                             ta1_y_yyy_zzzz_1,   \
                             ta2_yy_yy_xxxx_0,   \
                             ta2_yy_yy_xxxx_1,   \
                             ta2_yy_yy_xxxy_0,   \
                             ta2_yy_yy_xxxy_1,   \
                             ta2_yy_yy_xxxz_0,   \
                             ta2_yy_yy_xxxz_1,   \
                             ta2_yy_yy_xxyy_0,   \
                             ta2_yy_yy_xxyy_1,   \
                             ta2_yy_yy_xxyz_0,   \
                             ta2_yy_yy_xxyz_1,   \
                             ta2_yy_yy_xxzz_0,   \
                             ta2_yy_yy_xxzz_1,   \
                             ta2_yy_yy_xyyy_0,   \
                             ta2_yy_yy_xyyy_1,   \
                             ta2_yy_yy_xyyz_0,   \
                             ta2_yy_yy_xyyz_1,   \
                             ta2_yy_yy_xyzz_0,   \
                             ta2_yy_yy_xyzz_1,   \
                             ta2_yy_yy_xzzz_0,   \
                             ta2_yy_yy_xzzz_1,   \
                             ta2_yy_yy_yyyy_0,   \
                             ta2_yy_yy_yyyy_1,   \
                             ta2_yy_yy_yyyz_0,   \
                             ta2_yy_yy_yyyz_1,   \
                             ta2_yy_yy_yyzz_0,   \
                             ta2_yy_yy_yyzz_1,   \
                             ta2_yy_yy_yzzz_0,   \
                             ta2_yy_yy_yzzz_1,   \
                             ta2_yy_yy_zzzz_0,   \
                             ta2_yy_yy_zzzz_1,   \
                             ta2_yy_yyy_xxx_0,   \
                             ta2_yy_yyy_xxx_1,   \
                             ta2_yy_yyy_xxxx_0,  \
                             ta2_yy_yyy_xxxx_1,  \
                             ta2_yy_yyy_xxxy_0,  \
                             ta2_yy_yyy_xxxy_1,  \
                             ta2_yy_yyy_xxxz_0,  \
                             ta2_yy_yyy_xxxz_1,  \
                             ta2_yy_yyy_xxy_0,   \
                             ta2_yy_yyy_xxy_1,   \
                             ta2_yy_yyy_xxyy_0,  \
                             ta2_yy_yyy_xxyy_1,  \
                             ta2_yy_yyy_xxyz_0,  \
                             ta2_yy_yyy_xxyz_1,  \
                             ta2_yy_yyy_xxz_0,   \
                             ta2_yy_yyy_xxz_1,   \
                             ta2_yy_yyy_xxzz_0,  \
                             ta2_yy_yyy_xxzz_1,  \
                             ta2_yy_yyy_xyy_0,   \
                             ta2_yy_yyy_xyy_1,   \
                             ta2_yy_yyy_xyyy_0,  \
                             ta2_yy_yyy_xyyy_1,  \
                             ta2_yy_yyy_xyyz_0,  \
                             ta2_yy_yyy_xyyz_1,  \
                             ta2_yy_yyy_xyz_0,   \
                             ta2_yy_yyy_xyz_1,   \
                             ta2_yy_yyy_xyzz_0,  \
                             ta2_yy_yyy_xyzz_1,  \
                             ta2_yy_yyy_xzz_0,   \
                             ta2_yy_yyy_xzz_1,   \
                             ta2_yy_yyy_xzzz_0,  \
                             ta2_yy_yyy_xzzz_1,  \
                             ta2_yy_yyy_yyy_0,   \
                             ta2_yy_yyy_yyy_1,   \
                             ta2_yy_yyy_yyyy_0,  \
                             ta2_yy_yyy_yyyy_1,  \
                             ta2_yy_yyy_yyyz_0,  \
                             ta2_yy_yyy_yyyz_1,  \
                             ta2_yy_yyy_yyz_0,   \
                             ta2_yy_yyy_yyz_1,   \
                             ta2_yy_yyy_yyzz_0,  \
                             ta2_yy_yyy_yyzz_1,  \
                             ta2_yy_yyy_yzz_0,   \
                             ta2_yy_yyy_yzz_1,   \
                             ta2_yy_yyy_yzzz_0,  \
                             ta2_yy_yyy_yzzz_1,  \
                             ta2_yy_yyy_zzz_0,   \
                             ta2_yy_yyy_zzz_1,   \
                             ta2_yy_yyy_zzzz_0,  \
                             ta2_yy_yyy_zzzz_1,  \
                             ta2_yy_yyyy_xxxx_0, \
                             ta2_yy_yyyy_xxxy_0, \
                             ta2_yy_yyyy_xxxz_0, \
                             ta2_yy_yyyy_xxyy_0, \
                             ta2_yy_yyyy_xxyz_0, \
                             ta2_yy_yyyy_xxzz_0, \
                             ta2_yy_yyyy_xyyy_0, \
                             ta2_yy_yyyy_xyyz_0, \
                             ta2_yy_yyyy_xyzz_0, \
                             ta2_yy_yyyy_xzzz_0, \
                             ta2_yy_yyyy_yyyy_0, \
                             ta2_yy_yyyy_yyyz_0, \
                             ta2_yy_yyyy_yyzz_0, \
                             ta2_yy_yyyy_yzzz_0, \
                             ta2_yy_yyyy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyyy_xxxx_0[i] = 3.0 * ta2_yy_yy_xxxx_0[i] * fe_0 - 3.0 * ta2_yy_yy_xxxx_1[i] * fe_0 + 2.0 * ta1_y_yyy_xxxx_1[i] +
                                ta2_yy_yyy_xxxx_0[i] * pa_y[i] - ta2_yy_yyy_xxxx_1[i] * pc_y[i];

        ta2_yy_yyyy_xxxy_0[i] = 3.0 * ta2_yy_yy_xxxy_0[i] * fe_0 - 3.0 * ta2_yy_yy_xxxy_1[i] * fe_0 + ta2_yy_yyy_xxx_0[i] * fe_0 -
                                ta2_yy_yyy_xxx_1[i] * fe_0 + 2.0 * ta1_y_yyy_xxxy_1[i] + ta2_yy_yyy_xxxy_0[i] * pa_y[i] -
                                ta2_yy_yyy_xxxy_1[i] * pc_y[i];

        ta2_yy_yyyy_xxxz_0[i] = 3.0 * ta2_yy_yy_xxxz_0[i] * fe_0 - 3.0 * ta2_yy_yy_xxxz_1[i] * fe_0 + 2.0 * ta1_y_yyy_xxxz_1[i] +
                                ta2_yy_yyy_xxxz_0[i] * pa_y[i] - ta2_yy_yyy_xxxz_1[i] * pc_y[i];

        ta2_yy_yyyy_xxyy_0[i] = 3.0 * ta2_yy_yy_xxyy_0[i] * fe_0 - 3.0 * ta2_yy_yy_xxyy_1[i] * fe_0 + 2.0 * ta2_yy_yyy_xxy_0[i] * fe_0 -
                                2.0 * ta2_yy_yyy_xxy_1[i] * fe_0 + 2.0 * ta1_y_yyy_xxyy_1[i] + ta2_yy_yyy_xxyy_0[i] * pa_y[i] -
                                ta2_yy_yyy_xxyy_1[i] * pc_y[i];

        ta2_yy_yyyy_xxyz_0[i] = 3.0 * ta2_yy_yy_xxyz_0[i] * fe_0 - 3.0 * ta2_yy_yy_xxyz_1[i] * fe_0 + ta2_yy_yyy_xxz_0[i] * fe_0 -
                                ta2_yy_yyy_xxz_1[i] * fe_0 + 2.0 * ta1_y_yyy_xxyz_1[i] + ta2_yy_yyy_xxyz_0[i] * pa_y[i] -
                                ta2_yy_yyy_xxyz_1[i] * pc_y[i];

        ta2_yy_yyyy_xxzz_0[i] = 3.0 * ta2_yy_yy_xxzz_0[i] * fe_0 - 3.0 * ta2_yy_yy_xxzz_1[i] * fe_0 + 2.0 * ta1_y_yyy_xxzz_1[i] +
                                ta2_yy_yyy_xxzz_0[i] * pa_y[i] - ta2_yy_yyy_xxzz_1[i] * pc_y[i];

        ta2_yy_yyyy_xyyy_0[i] = 3.0 * ta2_yy_yy_xyyy_0[i] * fe_0 - 3.0 * ta2_yy_yy_xyyy_1[i] * fe_0 + 3.0 * ta2_yy_yyy_xyy_0[i] * fe_0 -
                                3.0 * ta2_yy_yyy_xyy_1[i] * fe_0 + 2.0 * ta1_y_yyy_xyyy_1[i] + ta2_yy_yyy_xyyy_0[i] * pa_y[i] -
                                ta2_yy_yyy_xyyy_1[i] * pc_y[i];

        ta2_yy_yyyy_xyyz_0[i] = 3.0 * ta2_yy_yy_xyyz_0[i] * fe_0 - 3.0 * ta2_yy_yy_xyyz_1[i] * fe_0 + 2.0 * ta2_yy_yyy_xyz_0[i] * fe_0 -
                                2.0 * ta2_yy_yyy_xyz_1[i] * fe_0 + 2.0 * ta1_y_yyy_xyyz_1[i] + ta2_yy_yyy_xyyz_0[i] * pa_y[i] -
                                ta2_yy_yyy_xyyz_1[i] * pc_y[i];

        ta2_yy_yyyy_xyzz_0[i] = 3.0 * ta2_yy_yy_xyzz_0[i] * fe_0 - 3.0 * ta2_yy_yy_xyzz_1[i] * fe_0 + ta2_yy_yyy_xzz_0[i] * fe_0 -
                                ta2_yy_yyy_xzz_1[i] * fe_0 + 2.0 * ta1_y_yyy_xyzz_1[i] + ta2_yy_yyy_xyzz_0[i] * pa_y[i] -
                                ta2_yy_yyy_xyzz_1[i] * pc_y[i];

        ta2_yy_yyyy_xzzz_0[i] = 3.0 * ta2_yy_yy_xzzz_0[i] * fe_0 - 3.0 * ta2_yy_yy_xzzz_1[i] * fe_0 + 2.0 * ta1_y_yyy_xzzz_1[i] +
                                ta2_yy_yyy_xzzz_0[i] * pa_y[i] - ta2_yy_yyy_xzzz_1[i] * pc_y[i];

        ta2_yy_yyyy_yyyy_0[i] = 3.0 * ta2_yy_yy_yyyy_0[i] * fe_0 - 3.0 * ta2_yy_yy_yyyy_1[i] * fe_0 + 4.0 * ta2_yy_yyy_yyy_0[i] * fe_0 -
                                4.0 * ta2_yy_yyy_yyy_1[i] * fe_0 + 2.0 * ta1_y_yyy_yyyy_1[i] + ta2_yy_yyy_yyyy_0[i] * pa_y[i] -
                                ta2_yy_yyy_yyyy_1[i] * pc_y[i];

        ta2_yy_yyyy_yyyz_0[i] = 3.0 * ta2_yy_yy_yyyz_0[i] * fe_0 - 3.0 * ta2_yy_yy_yyyz_1[i] * fe_0 + 3.0 * ta2_yy_yyy_yyz_0[i] * fe_0 -
                                3.0 * ta2_yy_yyy_yyz_1[i] * fe_0 + 2.0 * ta1_y_yyy_yyyz_1[i] + ta2_yy_yyy_yyyz_0[i] * pa_y[i] -
                                ta2_yy_yyy_yyyz_1[i] * pc_y[i];

        ta2_yy_yyyy_yyzz_0[i] = 3.0 * ta2_yy_yy_yyzz_0[i] * fe_0 - 3.0 * ta2_yy_yy_yyzz_1[i] * fe_0 + 2.0 * ta2_yy_yyy_yzz_0[i] * fe_0 -
                                2.0 * ta2_yy_yyy_yzz_1[i] * fe_0 + 2.0 * ta1_y_yyy_yyzz_1[i] + ta2_yy_yyy_yyzz_0[i] * pa_y[i] -
                                ta2_yy_yyy_yyzz_1[i] * pc_y[i];

        ta2_yy_yyyy_yzzz_0[i] = 3.0 * ta2_yy_yy_yzzz_0[i] * fe_0 - 3.0 * ta2_yy_yy_yzzz_1[i] * fe_0 + ta2_yy_yyy_zzz_0[i] * fe_0 -
                                ta2_yy_yyy_zzz_1[i] * fe_0 + 2.0 * ta1_y_yyy_yzzz_1[i] + ta2_yy_yyy_yzzz_0[i] * pa_y[i] -
                                ta2_yy_yyy_yzzz_1[i] * pc_y[i];

        ta2_yy_yyyy_zzzz_0[i] = 3.0 * ta2_yy_yy_zzzz_0[i] * fe_0 - 3.0 * ta2_yy_yy_zzzz_1[i] * fe_0 + 2.0 * ta1_y_yyy_zzzz_1[i] +
                                ta2_yy_yyy_zzzz_0[i] * pa_y[i] - ta2_yy_yyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 840-855 components of targeted buffer : GG

    auto ta2_yy_yyyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 840);

    auto ta2_yy_yyyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 841);

    auto ta2_yy_yyyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 842);

    auto ta2_yy_yyyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 843);

    auto ta2_yy_yyyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 844);

    auto ta2_yy_yyyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 845);

    auto ta2_yy_yyyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 846);

    auto ta2_yy_yyyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 847);

    auto ta2_yy_yyyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 848);

    auto ta2_yy_yyyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 849);

    auto ta2_yy_yyyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 850);

    auto ta2_yy_yyyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 851);

    auto ta2_yy_yyyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 852);

    auto ta2_yy_yyyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 853);

    auto ta2_yy_yyyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 854);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta2_yy_yyy_xxx_0,   \
                             ta2_yy_yyy_xxx_1,   \
                             ta2_yy_yyy_xxxx_0,  \
                             ta2_yy_yyy_xxxx_1,  \
                             ta2_yy_yyy_xxxy_0,  \
                             ta2_yy_yyy_xxxy_1,  \
                             ta2_yy_yyy_xxxz_0,  \
                             ta2_yy_yyy_xxxz_1,  \
                             ta2_yy_yyy_xxy_0,   \
                             ta2_yy_yyy_xxy_1,   \
                             ta2_yy_yyy_xxyy_0,  \
                             ta2_yy_yyy_xxyy_1,  \
                             ta2_yy_yyy_xxyz_0,  \
                             ta2_yy_yyy_xxyz_1,  \
                             ta2_yy_yyy_xxz_0,   \
                             ta2_yy_yyy_xxz_1,   \
                             ta2_yy_yyy_xxzz_0,  \
                             ta2_yy_yyy_xxzz_1,  \
                             ta2_yy_yyy_xyy_0,   \
                             ta2_yy_yyy_xyy_1,   \
                             ta2_yy_yyy_xyyy_0,  \
                             ta2_yy_yyy_xyyy_1,  \
                             ta2_yy_yyy_xyyz_0,  \
                             ta2_yy_yyy_xyyz_1,  \
                             ta2_yy_yyy_xyz_0,   \
                             ta2_yy_yyy_xyz_1,   \
                             ta2_yy_yyy_xyzz_0,  \
                             ta2_yy_yyy_xyzz_1,  \
                             ta2_yy_yyy_xzz_0,   \
                             ta2_yy_yyy_xzz_1,   \
                             ta2_yy_yyy_xzzz_0,  \
                             ta2_yy_yyy_xzzz_1,  \
                             ta2_yy_yyy_yyy_0,   \
                             ta2_yy_yyy_yyy_1,   \
                             ta2_yy_yyy_yyyy_0,  \
                             ta2_yy_yyy_yyyy_1,  \
                             ta2_yy_yyy_yyyz_0,  \
                             ta2_yy_yyy_yyyz_1,  \
                             ta2_yy_yyy_yyz_0,   \
                             ta2_yy_yyy_yyz_1,   \
                             ta2_yy_yyy_yyzz_0,  \
                             ta2_yy_yyy_yyzz_1,  \
                             ta2_yy_yyy_yzz_0,   \
                             ta2_yy_yyy_yzz_1,   \
                             ta2_yy_yyy_yzzz_0,  \
                             ta2_yy_yyy_yzzz_1,  \
                             ta2_yy_yyy_zzz_0,   \
                             ta2_yy_yyy_zzz_1,   \
                             ta2_yy_yyy_zzzz_0,  \
                             ta2_yy_yyy_zzzz_1,  \
                             ta2_yy_yyyz_xxxx_0, \
                             ta2_yy_yyyz_xxxy_0, \
                             ta2_yy_yyyz_xxxz_0, \
                             ta2_yy_yyyz_xxyy_0, \
                             ta2_yy_yyyz_xxyz_0, \
                             ta2_yy_yyyz_xxzz_0, \
                             ta2_yy_yyyz_xyyy_0, \
                             ta2_yy_yyyz_xyyz_0, \
                             ta2_yy_yyyz_xyzz_0, \
                             ta2_yy_yyyz_xzzz_0, \
                             ta2_yy_yyyz_yyyy_0, \
                             ta2_yy_yyyz_yyyz_0, \
                             ta2_yy_yyyz_yyzz_0, \
                             ta2_yy_yyyz_yzzz_0, \
                             ta2_yy_yyyz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyyz_xxxx_0[i] = ta2_yy_yyy_xxxx_0[i] * pa_z[i] - ta2_yy_yyy_xxxx_1[i] * pc_z[i];

        ta2_yy_yyyz_xxxy_0[i] = ta2_yy_yyy_xxxy_0[i] * pa_z[i] - ta2_yy_yyy_xxxy_1[i] * pc_z[i];

        ta2_yy_yyyz_xxxz_0[i] =
            ta2_yy_yyy_xxx_0[i] * fe_0 - ta2_yy_yyy_xxx_1[i] * fe_0 + ta2_yy_yyy_xxxz_0[i] * pa_z[i] - ta2_yy_yyy_xxxz_1[i] * pc_z[i];

        ta2_yy_yyyz_xxyy_0[i] = ta2_yy_yyy_xxyy_0[i] * pa_z[i] - ta2_yy_yyy_xxyy_1[i] * pc_z[i];

        ta2_yy_yyyz_xxyz_0[i] =
            ta2_yy_yyy_xxy_0[i] * fe_0 - ta2_yy_yyy_xxy_1[i] * fe_0 + ta2_yy_yyy_xxyz_0[i] * pa_z[i] - ta2_yy_yyy_xxyz_1[i] * pc_z[i];

        ta2_yy_yyyz_xxzz_0[i] =
            2.0 * ta2_yy_yyy_xxz_0[i] * fe_0 - 2.0 * ta2_yy_yyy_xxz_1[i] * fe_0 + ta2_yy_yyy_xxzz_0[i] * pa_z[i] - ta2_yy_yyy_xxzz_1[i] * pc_z[i];

        ta2_yy_yyyz_xyyy_0[i] = ta2_yy_yyy_xyyy_0[i] * pa_z[i] - ta2_yy_yyy_xyyy_1[i] * pc_z[i];

        ta2_yy_yyyz_xyyz_0[i] =
            ta2_yy_yyy_xyy_0[i] * fe_0 - ta2_yy_yyy_xyy_1[i] * fe_0 + ta2_yy_yyy_xyyz_0[i] * pa_z[i] - ta2_yy_yyy_xyyz_1[i] * pc_z[i];

        ta2_yy_yyyz_xyzz_0[i] =
            2.0 * ta2_yy_yyy_xyz_0[i] * fe_0 - 2.0 * ta2_yy_yyy_xyz_1[i] * fe_0 + ta2_yy_yyy_xyzz_0[i] * pa_z[i] - ta2_yy_yyy_xyzz_1[i] * pc_z[i];

        ta2_yy_yyyz_xzzz_0[i] =
            3.0 * ta2_yy_yyy_xzz_0[i] * fe_0 - 3.0 * ta2_yy_yyy_xzz_1[i] * fe_0 + ta2_yy_yyy_xzzz_0[i] * pa_z[i] - ta2_yy_yyy_xzzz_1[i] * pc_z[i];

        ta2_yy_yyyz_yyyy_0[i] = ta2_yy_yyy_yyyy_0[i] * pa_z[i] - ta2_yy_yyy_yyyy_1[i] * pc_z[i];

        ta2_yy_yyyz_yyyz_0[i] =
            ta2_yy_yyy_yyy_0[i] * fe_0 - ta2_yy_yyy_yyy_1[i] * fe_0 + ta2_yy_yyy_yyyz_0[i] * pa_z[i] - ta2_yy_yyy_yyyz_1[i] * pc_z[i];

        ta2_yy_yyyz_yyzz_0[i] =
            2.0 * ta2_yy_yyy_yyz_0[i] * fe_0 - 2.0 * ta2_yy_yyy_yyz_1[i] * fe_0 + ta2_yy_yyy_yyzz_0[i] * pa_z[i] - ta2_yy_yyy_yyzz_1[i] * pc_z[i];

        ta2_yy_yyyz_yzzz_0[i] =
            3.0 * ta2_yy_yyy_yzz_0[i] * fe_0 - 3.0 * ta2_yy_yyy_yzz_1[i] * fe_0 + ta2_yy_yyy_yzzz_0[i] * pa_z[i] - ta2_yy_yyy_yzzz_1[i] * pc_z[i];

        ta2_yy_yyyz_zzzz_0[i] =
            4.0 * ta2_yy_yyy_zzz_0[i] * fe_0 - 4.0 * ta2_yy_yyy_zzz_1[i] * fe_0 + ta2_yy_yyy_zzzz_0[i] * pa_z[i] - ta2_yy_yyy_zzzz_1[i] * pc_z[i];
    }

    // Set up 855-870 components of targeted buffer : GG

    auto ta2_yy_yyzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 855);

    auto ta2_yy_yyzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 856);

    auto ta2_yy_yyzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 857);

    auto ta2_yy_yyzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 858);

    auto ta2_yy_yyzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 859);

    auto ta2_yy_yyzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 860);

    auto ta2_yy_yyzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 861);

    auto ta2_yy_yyzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 862);

    auto ta2_yy_yyzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 863);

    auto ta2_yy_yyzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 864);

    auto ta2_yy_yyzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 865);

    auto ta2_yy_yyzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 866);

    auto ta2_yy_yyzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 867);

    auto ta2_yy_yyzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 868);

    auto ta2_yy_yyzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 869);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_yzz_xxxz_1,   \
                             ta1_y_yzz_xxzz_1,   \
                             ta1_y_yzz_xzzz_1,   \
                             ta1_y_yzz_zzzz_1,   \
                             ta2_yy_yy_xxxx_0,   \
                             ta2_yy_yy_xxxx_1,   \
                             ta2_yy_yy_xxxy_0,   \
                             ta2_yy_yy_xxxy_1,   \
                             ta2_yy_yy_xxyy_0,   \
                             ta2_yy_yy_xxyy_1,   \
                             ta2_yy_yy_xxyz_0,   \
                             ta2_yy_yy_xxyz_1,   \
                             ta2_yy_yy_xyyy_0,   \
                             ta2_yy_yy_xyyy_1,   \
                             ta2_yy_yy_xyyz_0,   \
                             ta2_yy_yy_xyyz_1,   \
                             ta2_yy_yy_xyzz_0,   \
                             ta2_yy_yy_xyzz_1,   \
                             ta2_yy_yy_yyyy_0,   \
                             ta2_yy_yy_yyyy_1,   \
                             ta2_yy_yy_yyyz_0,   \
                             ta2_yy_yy_yyyz_1,   \
                             ta2_yy_yy_yyzz_0,   \
                             ta2_yy_yy_yyzz_1,   \
                             ta2_yy_yy_yzzz_0,   \
                             ta2_yy_yy_yzzz_1,   \
                             ta2_yy_yyz_xxxx_0,  \
                             ta2_yy_yyz_xxxx_1,  \
                             ta2_yy_yyz_xxxy_0,  \
                             ta2_yy_yyz_xxxy_1,  \
                             ta2_yy_yyz_xxy_0,   \
                             ta2_yy_yyz_xxy_1,   \
                             ta2_yy_yyz_xxyy_0,  \
                             ta2_yy_yyz_xxyy_1,  \
                             ta2_yy_yyz_xxyz_0,  \
                             ta2_yy_yyz_xxyz_1,  \
                             ta2_yy_yyz_xyy_0,   \
                             ta2_yy_yyz_xyy_1,   \
                             ta2_yy_yyz_xyyy_0,  \
                             ta2_yy_yyz_xyyy_1,  \
                             ta2_yy_yyz_xyyz_0,  \
                             ta2_yy_yyz_xyyz_1,  \
                             ta2_yy_yyz_xyz_0,   \
                             ta2_yy_yyz_xyz_1,   \
                             ta2_yy_yyz_xyzz_0,  \
                             ta2_yy_yyz_xyzz_1,  \
                             ta2_yy_yyz_yyy_0,   \
                             ta2_yy_yyz_yyy_1,   \
                             ta2_yy_yyz_yyyy_0,  \
                             ta2_yy_yyz_yyyy_1,  \
                             ta2_yy_yyz_yyyz_0,  \
                             ta2_yy_yyz_yyyz_1,  \
                             ta2_yy_yyz_yyz_0,   \
                             ta2_yy_yyz_yyz_1,   \
                             ta2_yy_yyz_yyzz_0,  \
                             ta2_yy_yyz_yyzz_1,  \
                             ta2_yy_yyz_yzz_0,   \
                             ta2_yy_yyz_yzz_1,   \
                             ta2_yy_yyz_yzzz_0,  \
                             ta2_yy_yyz_yzzz_1,  \
                             ta2_yy_yyzz_xxxx_0, \
                             ta2_yy_yyzz_xxxy_0, \
                             ta2_yy_yyzz_xxxz_0, \
                             ta2_yy_yyzz_xxyy_0, \
                             ta2_yy_yyzz_xxyz_0, \
                             ta2_yy_yyzz_xxzz_0, \
                             ta2_yy_yyzz_xyyy_0, \
                             ta2_yy_yyzz_xyyz_0, \
                             ta2_yy_yyzz_xyzz_0, \
                             ta2_yy_yyzz_xzzz_0, \
                             ta2_yy_yyzz_yyyy_0, \
                             ta2_yy_yyzz_yyyz_0, \
                             ta2_yy_yyzz_yyzz_0, \
                             ta2_yy_yyzz_yzzz_0, \
                             ta2_yy_yyzz_zzzz_0, \
                             ta2_yy_yzz_xxxz_0,  \
                             ta2_yy_yzz_xxxz_1,  \
                             ta2_yy_yzz_xxzz_0,  \
                             ta2_yy_yzz_xxzz_1,  \
                             ta2_yy_yzz_xzzz_0,  \
                             ta2_yy_yzz_xzzz_1,  \
                             ta2_yy_yzz_zzzz_0,  \
                             ta2_yy_yzz_zzzz_1,  \
                             ta2_yy_zz_xxxz_0,   \
                             ta2_yy_zz_xxxz_1,   \
                             ta2_yy_zz_xxzz_0,   \
                             ta2_yy_zz_xxzz_1,   \
                             ta2_yy_zz_xzzz_0,   \
                             ta2_yy_zz_xzzz_1,   \
                             ta2_yy_zz_zzzz_0,   \
                             ta2_yy_zz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yyzz_xxxx_0[i] =
            ta2_yy_yy_xxxx_0[i] * fe_0 - ta2_yy_yy_xxxx_1[i] * fe_0 + ta2_yy_yyz_xxxx_0[i] * pa_z[i] - ta2_yy_yyz_xxxx_1[i] * pc_z[i];

        ta2_yy_yyzz_xxxy_0[i] =
            ta2_yy_yy_xxxy_0[i] * fe_0 - ta2_yy_yy_xxxy_1[i] * fe_0 + ta2_yy_yyz_xxxy_0[i] * pa_z[i] - ta2_yy_yyz_xxxy_1[i] * pc_z[i];

        ta2_yy_yyzz_xxxz_0[i] = ta2_yy_zz_xxxz_0[i] * fe_0 - ta2_yy_zz_xxxz_1[i] * fe_0 + 2.0 * ta1_y_yzz_xxxz_1[i] + ta2_yy_yzz_xxxz_0[i] * pa_y[i] -
                                ta2_yy_yzz_xxxz_1[i] * pc_y[i];

        ta2_yy_yyzz_xxyy_0[i] =
            ta2_yy_yy_xxyy_0[i] * fe_0 - ta2_yy_yy_xxyy_1[i] * fe_0 + ta2_yy_yyz_xxyy_0[i] * pa_z[i] - ta2_yy_yyz_xxyy_1[i] * pc_z[i];

        ta2_yy_yyzz_xxyz_0[i] = ta2_yy_yy_xxyz_0[i] * fe_0 - ta2_yy_yy_xxyz_1[i] * fe_0 + ta2_yy_yyz_xxy_0[i] * fe_0 - ta2_yy_yyz_xxy_1[i] * fe_0 +
                                ta2_yy_yyz_xxyz_0[i] * pa_z[i] - ta2_yy_yyz_xxyz_1[i] * pc_z[i];

        ta2_yy_yyzz_xxzz_0[i] = ta2_yy_zz_xxzz_0[i] * fe_0 - ta2_yy_zz_xxzz_1[i] * fe_0 + 2.0 * ta1_y_yzz_xxzz_1[i] + ta2_yy_yzz_xxzz_0[i] * pa_y[i] -
                                ta2_yy_yzz_xxzz_1[i] * pc_y[i];

        ta2_yy_yyzz_xyyy_0[i] =
            ta2_yy_yy_xyyy_0[i] * fe_0 - ta2_yy_yy_xyyy_1[i] * fe_0 + ta2_yy_yyz_xyyy_0[i] * pa_z[i] - ta2_yy_yyz_xyyy_1[i] * pc_z[i];

        ta2_yy_yyzz_xyyz_0[i] = ta2_yy_yy_xyyz_0[i] * fe_0 - ta2_yy_yy_xyyz_1[i] * fe_0 + ta2_yy_yyz_xyy_0[i] * fe_0 - ta2_yy_yyz_xyy_1[i] * fe_0 +
                                ta2_yy_yyz_xyyz_0[i] * pa_z[i] - ta2_yy_yyz_xyyz_1[i] * pc_z[i];

        ta2_yy_yyzz_xyzz_0[i] = ta2_yy_yy_xyzz_0[i] * fe_0 - ta2_yy_yy_xyzz_1[i] * fe_0 + 2.0 * ta2_yy_yyz_xyz_0[i] * fe_0 -
                                2.0 * ta2_yy_yyz_xyz_1[i] * fe_0 + ta2_yy_yyz_xyzz_0[i] * pa_z[i] - ta2_yy_yyz_xyzz_1[i] * pc_z[i];

        ta2_yy_yyzz_xzzz_0[i] = ta2_yy_zz_xzzz_0[i] * fe_0 - ta2_yy_zz_xzzz_1[i] * fe_0 + 2.0 * ta1_y_yzz_xzzz_1[i] + ta2_yy_yzz_xzzz_0[i] * pa_y[i] -
                                ta2_yy_yzz_xzzz_1[i] * pc_y[i];

        ta2_yy_yyzz_yyyy_0[i] =
            ta2_yy_yy_yyyy_0[i] * fe_0 - ta2_yy_yy_yyyy_1[i] * fe_0 + ta2_yy_yyz_yyyy_0[i] * pa_z[i] - ta2_yy_yyz_yyyy_1[i] * pc_z[i];

        ta2_yy_yyzz_yyyz_0[i] = ta2_yy_yy_yyyz_0[i] * fe_0 - ta2_yy_yy_yyyz_1[i] * fe_0 + ta2_yy_yyz_yyy_0[i] * fe_0 - ta2_yy_yyz_yyy_1[i] * fe_0 +
                                ta2_yy_yyz_yyyz_0[i] * pa_z[i] - ta2_yy_yyz_yyyz_1[i] * pc_z[i];

        ta2_yy_yyzz_yyzz_0[i] = ta2_yy_yy_yyzz_0[i] * fe_0 - ta2_yy_yy_yyzz_1[i] * fe_0 + 2.0 * ta2_yy_yyz_yyz_0[i] * fe_0 -
                                2.0 * ta2_yy_yyz_yyz_1[i] * fe_0 + ta2_yy_yyz_yyzz_0[i] * pa_z[i] - ta2_yy_yyz_yyzz_1[i] * pc_z[i];

        ta2_yy_yyzz_yzzz_0[i] = ta2_yy_yy_yzzz_0[i] * fe_0 - ta2_yy_yy_yzzz_1[i] * fe_0 + 3.0 * ta2_yy_yyz_yzz_0[i] * fe_0 -
                                3.0 * ta2_yy_yyz_yzz_1[i] * fe_0 + ta2_yy_yyz_yzzz_0[i] * pa_z[i] - ta2_yy_yyz_yzzz_1[i] * pc_z[i];

        ta2_yy_yyzz_zzzz_0[i] = ta2_yy_zz_zzzz_0[i] * fe_0 - ta2_yy_zz_zzzz_1[i] * fe_0 + 2.0 * ta1_y_yzz_zzzz_1[i] + ta2_yy_yzz_zzzz_0[i] * pa_y[i] -
                                ta2_yy_yzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 870-885 components of targeted buffer : GG

    auto ta2_yy_yzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 870);

    auto ta2_yy_yzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 871);

    auto ta2_yy_yzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 872);

    auto ta2_yy_yzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 873);

    auto ta2_yy_yzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 874);

    auto ta2_yy_yzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 875);

    auto ta2_yy_yzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 876);

    auto ta2_yy_yzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 877);

    auto ta2_yy_yzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 878);

    auto ta2_yy_yzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 879);

    auto ta2_yy_yzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 880);

    auto ta2_yy_yzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 881);

    auto ta2_yy_yzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 882);

    auto ta2_yy_yzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 883);

    auto ta2_yy_yzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 884);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_zzz_xxxx_1,   \
                             ta1_y_zzz_xxxz_1,   \
                             ta1_y_zzz_xxyz_1,   \
                             ta1_y_zzz_xxzz_1,   \
                             ta1_y_zzz_xyyz_1,   \
                             ta1_y_zzz_xyzz_1,   \
                             ta1_y_zzz_xzzz_1,   \
                             ta1_y_zzz_yyyz_1,   \
                             ta1_y_zzz_yyzz_1,   \
                             ta1_y_zzz_yzzz_1,   \
                             ta1_y_zzz_zzzz_1,   \
                             ta2_yy_yz_xxxy_0,   \
                             ta2_yy_yz_xxxy_1,   \
                             ta2_yy_yz_xxyy_0,   \
                             ta2_yy_yz_xxyy_1,   \
                             ta2_yy_yz_xyyy_0,   \
                             ta2_yy_yz_xyyy_1,   \
                             ta2_yy_yz_yyyy_0,   \
                             ta2_yy_yz_yyyy_1,   \
                             ta2_yy_yzz_xxxy_0,  \
                             ta2_yy_yzz_xxxy_1,  \
                             ta2_yy_yzz_xxyy_0,  \
                             ta2_yy_yzz_xxyy_1,  \
                             ta2_yy_yzz_xyyy_0,  \
                             ta2_yy_yzz_xyyy_1,  \
                             ta2_yy_yzz_yyyy_0,  \
                             ta2_yy_yzz_yyyy_1,  \
                             ta2_yy_yzzz_xxxx_0, \
                             ta2_yy_yzzz_xxxy_0, \
                             ta2_yy_yzzz_xxxz_0, \
                             ta2_yy_yzzz_xxyy_0, \
                             ta2_yy_yzzz_xxyz_0, \
                             ta2_yy_yzzz_xxzz_0, \
                             ta2_yy_yzzz_xyyy_0, \
                             ta2_yy_yzzz_xyyz_0, \
                             ta2_yy_yzzz_xyzz_0, \
                             ta2_yy_yzzz_xzzz_0, \
                             ta2_yy_yzzz_yyyy_0, \
                             ta2_yy_yzzz_yyyz_0, \
                             ta2_yy_yzzz_yyzz_0, \
                             ta2_yy_yzzz_yzzz_0, \
                             ta2_yy_yzzz_zzzz_0, \
                             ta2_yy_zzz_xxxx_0,  \
                             ta2_yy_zzz_xxxx_1,  \
                             ta2_yy_zzz_xxxz_0,  \
                             ta2_yy_zzz_xxxz_1,  \
                             ta2_yy_zzz_xxyz_0,  \
                             ta2_yy_zzz_xxyz_1,  \
                             ta2_yy_zzz_xxz_0,   \
                             ta2_yy_zzz_xxz_1,   \
                             ta2_yy_zzz_xxzz_0,  \
                             ta2_yy_zzz_xxzz_1,  \
                             ta2_yy_zzz_xyyz_0,  \
                             ta2_yy_zzz_xyyz_1,  \
                             ta2_yy_zzz_xyz_0,   \
                             ta2_yy_zzz_xyz_1,   \
                             ta2_yy_zzz_xyzz_0,  \
                             ta2_yy_zzz_xyzz_1,  \
                             ta2_yy_zzz_xzz_0,   \
                             ta2_yy_zzz_xzz_1,   \
                             ta2_yy_zzz_xzzz_0,  \
                             ta2_yy_zzz_xzzz_1,  \
                             ta2_yy_zzz_yyyz_0,  \
                             ta2_yy_zzz_yyyz_1,  \
                             ta2_yy_zzz_yyz_0,   \
                             ta2_yy_zzz_yyz_1,   \
                             ta2_yy_zzz_yyzz_0,  \
                             ta2_yy_zzz_yyzz_1,  \
                             ta2_yy_zzz_yzz_0,   \
                             ta2_yy_zzz_yzz_1,   \
                             ta2_yy_zzz_yzzz_0,  \
                             ta2_yy_zzz_yzzz_1,  \
                             ta2_yy_zzz_zzz_0,   \
                             ta2_yy_zzz_zzz_1,   \
                             ta2_yy_zzz_zzzz_0,  \
                             ta2_yy_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_yzzz_xxxx_0[i] = 2.0 * ta1_y_zzz_xxxx_1[i] + ta2_yy_zzz_xxxx_0[i] * pa_y[i] - ta2_yy_zzz_xxxx_1[i] * pc_y[i];

        ta2_yy_yzzz_xxxy_0[i] =
            2.0 * ta2_yy_yz_xxxy_0[i] * fe_0 - 2.0 * ta2_yy_yz_xxxy_1[i] * fe_0 + ta2_yy_yzz_xxxy_0[i] * pa_z[i] - ta2_yy_yzz_xxxy_1[i] * pc_z[i];

        ta2_yy_yzzz_xxxz_0[i] = 2.0 * ta1_y_zzz_xxxz_1[i] + ta2_yy_zzz_xxxz_0[i] * pa_y[i] - ta2_yy_zzz_xxxz_1[i] * pc_y[i];

        ta2_yy_yzzz_xxyy_0[i] =
            2.0 * ta2_yy_yz_xxyy_0[i] * fe_0 - 2.0 * ta2_yy_yz_xxyy_1[i] * fe_0 + ta2_yy_yzz_xxyy_0[i] * pa_z[i] - ta2_yy_yzz_xxyy_1[i] * pc_z[i];

        ta2_yy_yzzz_xxyz_0[i] = ta2_yy_zzz_xxz_0[i] * fe_0 - ta2_yy_zzz_xxz_1[i] * fe_0 + 2.0 * ta1_y_zzz_xxyz_1[i] + ta2_yy_zzz_xxyz_0[i] * pa_y[i] -
                                ta2_yy_zzz_xxyz_1[i] * pc_y[i];

        ta2_yy_yzzz_xxzz_0[i] = 2.0 * ta1_y_zzz_xxzz_1[i] + ta2_yy_zzz_xxzz_0[i] * pa_y[i] - ta2_yy_zzz_xxzz_1[i] * pc_y[i];

        ta2_yy_yzzz_xyyy_0[i] =
            2.0 * ta2_yy_yz_xyyy_0[i] * fe_0 - 2.0 * ta2_yy_yz_xyyy_1[i] * fe_0 + ta2_yy_yzz_xyyy_0[i] * pa_z[i] - ta2_yy_yzz_xyyy_1[i] * pc_z[i];

        ta2_yy_yzzz_xyyz_0[i] = 2.0 * ta2_yy_zzz_xyz_0[i] * fe_0 - 2.0 * ta2_yy_zzz_xyz_1[i] * fe_0 + 2.0 * ta1_y_zzz_xyyz_1[i] +
                                ta2_yy_zzz_xyyz_0[i] * pa_y[i] - ta2_yy_zzz_xyyz_1[i] * pc_y[i];

        ta2_yy_yzzz_xyzz_0[i] = ta2_yy_zzz_xzz_0[i] * fe_0 - ta2_yy_zzz_xzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_xyzz_1[i] + ta2_yy_zzz_xyzz_0[i] * pa_y[i] -
                                ta2_yy_zzz_xyzz_1[i] * pc_y[i];

        ta2_yy_yzzz_xzzz_0[i] = 2.0 * ta1_y_zzz_xzzz_1[i] + ta2_yy_zzz_xzzz_0[i] * pa_y[i] - ta2_yy_zzz_xzzz_1[i] * pc_y[i];

        ta2_yy_yzzz_yyyy_0[i] =
            2.0 * ta2_yy_yz_yyyy_0[i] * fe_0 - 2.0 * ta2_yy_yz_yyyy_1[i] * fe_0 + ta2_yy_yzz_yyyy_0[i] * pa_z[i] - ta2_yy_yzz_yyyy_1[i] * pc_z[i];

        ta2_yy_yzzz_yyyz_0[i] = 3.0 * ta2_yy_zzz_yyz_0[i] * fe_0 - 3.0 * ta2_yy_zzz_yyz_1[i] * fe_0 + 2.0 * ta1_y_zzz_yyyz_1[i] +
                                ta2_yy_zzz_yyyz_0[i] * pa_y[i] - ta2_yy_zzz_yyyz_1[i] * pc_y[i];

        ta2_yy_yzzz_yyzz_0[i] = 2.0 * ta2_yy_zzz_yzz_0[i] * fe_0 - 2.0 * ta2_yy_zzz_yzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_yyzz_1[i] +
                                ta2_yy_zzz_yyzz_0[i] * pa_y[i] - ta2_yy_zzz_yyzz_1[i] * pc_y[i];

        ta2_yy_yzzz_yzzz_0[i] = ta2_yy_zzz_zzz_0[i] * fe_0 - ta2_yy_zzz_zzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_yzzz_1[i] + ta2_yy_zzz_yzzz_0[i] * pa_y[i] -
                                ta2_yy_zzz_yzzz_1[i] * pc_y[i];

        ta2_yy_yzzz_zzzz_0[i] = 2.0 * ta1_y_zzz_zzzz_1[i] + ta2_yy_zzz_zzzz_0[i] * pa_y[i] - ta2_yy_zzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 885-900 components of targeted buffer : GG

    auto ta2_yy_zzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 885);

    auto ta2_yy_zzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 886);

    auto ta2_yy_zzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 887);

    auto ta2_yy_zzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 888);

    auto ta2_yy_zzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 889);

    auto ta2_yy_zzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 890);

    auto ta2_yy_zzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 891);

    auto ta2_yy_zzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 892);

    auto ta2_yy_zzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 893);

    auto ta2_yy_zzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 894);

    auto ta2_yy_zzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 895);

    auto ta2_yy_zzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 896);

    auto ta2_yy_zzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 897);

    auto ta2_yy_zzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 898);

    auto ta2_yy_zzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 899);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta2_yy_zz_xxxx_0,   \
                             ta2_yy_zz_xxxx_1,   \
                             ta2_yy_zz_xxxy_0,   \
                             ta2_yy_zz_xxxy_1,   \
                             ta2_yy_zz_xxxz_0,   \
                             ta2_yy_zz_xxxz_1,   \
                             ta2_yy_zz_xxyy_0,   \
                             ta2_yy_zz_xxyy_1,   \
                             ta2_yy_zz_xxyz_0,   \
                             ta2_yy_zz_xxyz_1,   \
                             ta2_yy_zz_xxzz_0,   \
                             ta2_yy_zz_xxzz_1,   \
                             ta2_yy_zz_xyyy_0,   \
                             ta2_yy_zz_xyyy_1,   \
                             ta2_yy_zz_xyyz_0,   \
                             ta2_yy_zz_xyyz_1,   \
                             ta2_yy_zz_xyzz_0,   \
                             ta2_yy_zz_xyzz_1,   \
                             ta2_yy_zz_xzzz_0,   \
                             ta2_yy_zz_xzzz_1,   \
                             ta2_yy_zz_yyyy_0,   \
                             ta2_yy_zz_yyyy_1,   \
                             ta2_yy_zz_yyyz_0,   \
                             ta2_yy_zz_yyyz_1,   \
                             ta2_yy_zz_yyzz_0,   \
                             ta2_yy_zz_yyzz_1,   \
                             ta2_yy_zz_yzzz_0,   \
                             ta2_yy_zz_yzzz_1,   \
                             ta2_yy_zz_zzzz_0,   \
                             ta2_yy_zz_zzzz_1,   \
                             ta2_yy_zzz_xxx_0,   \
                             ta2_yy_zzz_xxx_1,   \
                             ta2_yy_zzz_xxxx_0,  \
                             ta2_yy_zzz_xxxx_1,  \
                             ta2_yy_zzz_xxxy_0,  \
                             ta2_yy_zzz_xxxy_1,  \
                             ta2_yy_zzz_xxxz_0,  \
                             ta2_yy_zzz_xxxz_1,  \
                             ta2_yy_zzz_xxy_0,   \
                             ta2_yy_zzz_xxy_1,   \
                             ta2_yy_zzz_xxyy_0,  \
                             ta2_yy_zzz_xxyy_1,  \
                             ta2_yy_zzz_xxyz_0,  \
                             ta2_yy_zzz_xxyz_1,  \
                             ta2_yy_zzz_xxz_0,   \
                             ta2_yy_zzz_xxz_1,   \
                             ta2_yy_zzz_xxzz_0,  \
                             ta2_yy_zzz_xxzz_1,  \
                             ta2_yy_zzz_xyy_0,   \
                             ta2_yy_zzz_xyy_1,   \
                             ta2_yy_zzz_xyyy_0,  \
                             ta2_yy_zzz_xyyy_1,  \
                             ta2_yy_zzz_xyyz_0,  \
                             ta2_yy_zzz_xyyz_1,  \
                             ta2_yy_zzz_xyz_0,   \
                             ta2_yy_zzz_xyz_1,   \
                             ta2_yy_zzz_xyzz_0,  \
                             ta2_yy_zzz_xyzz_1,  \
                             ta2_yy_zzz_xzz_0,   \
                             ta2_yy_zzz_xzz_1,   \
                             ta2_yy_zzz_xzzz_0,  \
                             ta2_yy_zzz_xzzz_1,  \
                             ta2_yy_zzz_yyy_0,   \
                             ta2_yy_zzz_yyy_1,   \
                             ta2_yy_zzz_yyyy_0,  \
                             ta2_yy_zzz_yyyy_1,  \
                             ta2_yy_zzz_yyyz_0,  \
                             ta2_yy_zzz_yyyz_1,  \
                             ta2_yy_zzz_yyz_0,   \
                             ta2_yy_zzz_yyz_1,   \
                             ta2_yy_zzz_yyzz_0,  \
                             ta2_yy_zzz_yyzz_1,  \
                             ta2_yy_zzz_yzz_0,   \
                             ta2_yy_zzz_yzz_1,   \
                             ta2_yy_zzz_yzzz_0,  \
                             ta2_yy_zzz_yzzz_1,  \
                             ta2_yy_zzz_zzz_0,   \
                             ta2_yy_zzz_zzz_1,   \
                             ta2_yy_zzz_zzzz_0,  \
                             ta2_yy_zzz_zzzz_1,  \
                             ta2_yy_zzzz_xxxx_0, \
                             ta2_yy_zzzz_xxxy_0, \
                             ta2_yy_zzzz_xxxz_0, \
                             ta2_yy_zzzz_xxyy_0, \
                             ta2_yy_zzzz_xxyz_0, \
                             ta2_yy_zzzz_xxzz_0, \
                             ta2_yy_zzzz_xyyy_0, \
                             ta2_yy_zzzz_xyyz_0, \
                             ta2_yy_zzzz_xyzz_0, \
                             ta2_yy_zzzz_xzzz_0, \
                             ta2_yy_zzzz_yyyy_0, \
                             ta2_yy_zzzz_yyyz_0, \
                             ta2_yy_zzzz_yyzz_0, \
                             ta2_yy_zzzz_yzzz_0, \
                             ta2_yy_zzzz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yy_zzzz_xxxx_0[i] =
            3.0 * ta2_yy_zz_xxxx_0[i] * fe_0 - 3.0 * ta2_yy_zz_xxxx_1[i] * fe_0 + ta2_yy_zzz_xxxx_0[i] * pa_z[i] - ta2_yy_zzz_xxxx_1[i] * pc_z[i];

        ta2_yy_zzzz_xxxy_0[i] =
            3.0 * ta2_yy_zz_xxxy_0[i] * fe_0 - 3.0 * ta2_yy_zz_xxxy_1[i] * fe_0 + ta2_yy_zzz_xxxy_0[i] * pa_z[i] - ta2_yy_zzz_xxxy_1[i] * pc_z[i];

        ta2_yy_zzzz_xxxz_0[i] = 3.0 * ta2_yy_zz_xxxz_0[i] * fe_0 - 3.0 * ta2_yy_zz_xxxz_1[i] * fe_0 + ta2_yy_zzz_xxx_0[i] * fe_0 -
                                ta2_yy_zzz_xxx_1[i] * fe_0 + ta2_yy_zzz_xxxz_0[i] * pa_z[i] - ta2_yy_zzz_xxxz_1[i] * pc_z[i];

        ta2_yy_zzzz_xxyy_0[i] =
            3.0 * ta2_yy_zz_xxyy_0[i] * fe_0 - 3.0 * ta2_yy_zz_xxyy_1[i] * fe_0 + ta2_yy_zzz_xxyy_0[i] * pa_z[i] - ta2_yy_zzz_xxyy_1[i] * pc_z[i];

        ta2_yy_zzzz_xxyz_0[i] = 3.0 * ta2_yy_zz_xxyz_0[i] * fe_0 - 3.0 * ta2_yy_zz_xxyz_1[i] * fe_0 + ta2_yy_zzz_xxy_0[i] * fe_0 -
                                ta2_yy_zzz_xxy_1[i] * fe_0 + ta2_yy_zzz_xxyz_0[i] * pa_z[i] - ta2_yy_zzz_xxyz_1[i] * pc_z[i];

        ta2_yy_zzzz_xxzz_0[i] = 3.0 * ta2_yy_zz_xxzz_0[i] * fe_0 - 3.0 * ta2_yy_zz_xxzz_1[i] * fe_0 + 2.0 * ta2_yy_zzz_xxz_0[i] * fe_0 -
                                2.0 * ta2_yy_zzz_xxz_1[i] * fe_0 + ta2_yy_zzz_xxzz_0[i] * pa_z[i] - ta2_yy_zzz_xxzz_1[i] * pc_z[i];

        ta2_yy_zzzz_xyyy_0[i] =
            3.0 * ta2_yy_zz_xyyy_0[i] * fe_0 - 3.0 * ta2_yy_zz_xyyy_1[i] * fe_0 + ta2_yy_zzz_xyyy_0[i] * pa_z[i] - ta2_yy_zzz_xyyy_1[i] * pc_z[i];

        ta2_yy_zzzz_xyyz_0[i] = 3.0 * ta2_yy_zz_xyyz_0[i] * fe_0 - 3.0 * ta2_yy_zz_xyyz_1[i] * fe_0 + ta2_yy_zzz_xyy_0[i] * fe_0 -
                                ta2_yy_zzz_xyy_1[i] * fe_0 + ta2_yy_zzz_xyyz_0[i] * pa_z[i] - ta2_yy_zzz_xyyz_1[i] * pc_z[i];

        ta2_yy_zzzz_xyzz_0[i] = 3.0 * ta2_yy_zz_xyzz_0[i] * fe_0 - 3.0 * ta2_yy_zz_xyzz_1[i] * fe_0 + 2.0 * ta2_yy_zzz_xyz_0[i] * fe_0 -
                                2.0 * ta2_yy_zzz_xyz_1[i] * fe_0 + ta2_yy_zzz_xyzz_0[i] * pa_z[i] - ta2_yy_zzz_xyzz_1[i] * pc_z[i];

        ta2_yy_zzzz_xzzz_0[i] = 3.0 * ta2_yy_zz_xzzz_0[i] * fe_0 - 3.0 * ta2_yy_zz_xzzz_1[i] * fe_0 + 3.0 * ta2_yy_zzz_xzz_0[i] * fe_0 -
                                3.0 * ta2_yy_zzz_xzz_1[i] * fe_0 + ta2_yy_zzz_xzzz_0[i] * pa_z[i] - ta2_yy_zzz_xzzz_1[i] * pc_z[i];

        ta2_yy_zzzz_yyyy_0[i] =
            3.0 * ta2_yy_zz_yyyy_0[i] * fe_0 - 3.0 * ta2_yy_zz_yyyy_1[i] * fe_0 + ta2_yy_zzz_yyyy_0[i] * pa_z[i] - ta2_yy_zzz_yyyy_1[i] * pc_z[i];

        ta2_yy_zzzz_yyyz_0[i] = 3.0 * ta2_yy_zz_yyyz_0[i] * fe_0 - 3.0 * ta2_yy_zz_yyyz_1[i] * fe_0 + ta2_yy_zzz_yyy_0[i] * fe_0 -
                                ta2_yy_zzz_yyy_1[i] * fe_0 + ta2_yy_zzz_yyyz_0[i] * pa_z[i] - ta2_yy_zzz_yyyz_1[i] * pc_z[i];

        ta2_yy_zzzz_yyzz_0[i] = 3.0 * ta2_yy_zz_yyzz_0[i] * fe_0 - 3.0 * ta2_yy_zz_yyzz_1[i] * fe_0 + 2.0 * ta2_yy_zzz_yyz_0[i] * fe_0 -
                                2.0 * ta2_yy_zzz_yyz_1[i] * fe_0 + ta2_yy_zzz_yyzz_0[i] * pa_z[i] - ta2_yy_zzz_yyzz_1[i] * pc_z[i];

        ta2_yy_zzzz_yzzz_0[i] = 3.0 * ta2_yy_zz_yzzz_0[i] * fe_0 - 3.0 * ta2_yy_zz_yzzz_1[i] * fe_0 + 3.0 * ta2_yy_zzz_yzz_0[i] * fe_0 -
                                3.0 * ta2_yy_zzz_yzz_1[i] * fe_0 + ta2_yy_zzz_yzzz_0[i] * pa_z[i] - ta2_yy_zzz_yzzz_1[i] * pc_z[i];

        ta2_yy_zzzz_zzzz_0[i] = 3.0 * ta2_yy_zz_zzzz_0[i] * fe_0 - 3.0 * ta2_yy_zz_zzzz_1[i] * fe_0 + 4.0 * ta2_yy_zzz_zzz_0[i] * fe_0 -
                                4.0 * ta2_yy_zzz_zzz_1[i] * fe_0 + ta2_yy_zzz_zzzz_0[i] * pa_z[i] - ta2_yy_zzz_zzzz_1[i] * pc_z[i];
    }

    // Set up 900-915 components of targeted buffer : GG

    auto ta2_yz_xxxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 900);

    auto ta2_yz_xxxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 901);

    auto ta2_yz_xxxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 902);

    auto ta2_yz_xxxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 903);

    auto ta2_yz_xxxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 904);

    auto ta2_yz_xxxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 905);

    auto ta2_yz_xxxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 906);

    auto ta2_yz_xxxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 907);

    auto ta2_yz_xxxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 908);

    auto ta2_yz_xxxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 909);

    auto ta2_yz_xxxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 910);

    auto ta2_yz_xxxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 911);

    auto ta2_yz_xxxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 912);

    auto ta2_yz_xxxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 913);

    auto ta2_yz_xxxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 914);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta2_yz_xx_xxxx_0,   \
                             ta2_yz_xx_xxxx_1,   \
                             ta2_yz_xx_xxxy_0,   \
                             ta2_yz_xx_xxxy_1,   \
                             ta2_yz_xx_xxxz_0,   \
                             ta2_yz_xx_xxxz_1,   \
                             ta2_yz_xx_xxyy_0,   \
                             ta2_yz_xx_xxyy_1,   \
                             ta2_yz_xx_xxyz_0,   \
                             ta2_yz_xx_xxyz_1,   \
                             ta2_yz_xx_xxzz_0,   \
                             ta2_yz_xx_xxzz_1,   \
                             ta2_yz_xx_xyyy_0,   \
                             ta2_yz_xx_xyyy_1,   \
                             ta2_yz_xx_xyyz_0,   \
                             ta2_yz_xx_xyyz_1,   \
                             ta2_yz_xx_xyzz_0,   \
                             ta2_yz_xx_xyzz_1,   \
                             ta2_yz_xx_xzzz_0,   \
                             ta2_yz_xx_xzzz_1,   \
                             ta2_yz_xx_yyyy_0,   \
                             ta2_yz_xx_yyyy_1,   \
                             ta2_yz_xx_yyyz_0,   \
                             ta2_yz_xx_yyyz_1,   \
                             ta2_yz_xx_yyzz_0,   \
                             ta2_yz_xx_yyzz_1,   \
                             ta2_yz_xx_yzzz_0,   \
                             ta2_yz_xx_yzzz_1,   \
                             ta2_yz_xx_zzzz_0,   \
                             ta2_yz_xx_zzzz_1,   \
                             ta2_yz_xxx_xxx_0,   \
                             ta2_yz_xxx_xxx_1,   \
                             ta2_yz_xxx_xxxx_0,  \
                             ta2_yz_xxx_xxxx_1,  \
                             ta2_yz_xxx_xxxy_0,  \
                             ta2_yz_xxx_xxxy_1,  \
                             ta2_yz_xxx_xxxz_0,  \
                             ta2_yz_xxx_xxxz_1,  \
                             ta2_yz_xxx_xxy_0,   \
                             ta2_yz_xxx_xxy_1,   \
                             ta2_yz_xxx_xxyy_0,  \
                             ta2_yz_xxx_xxyy_1,  \
                             ta2_yz_xxx_xxyz_0,  \
                             ta2_yz_xxx_xxyz_1,  \
                             ta2_yz_xxx_xxz_0,   \
                             ta2_yz_xxx_xxz_1,   \
                             ta2_yz_xxx_xxzz_0,  \
                             ta2_yz_xxx_xxzz_1,  \
                             ta2_yz_xxx_xyy_0,   \
                             ta2_yz_xxx_xyy_1,   \
                             ta2_yz_xxx_xyyy_0,  \
                             ta2_yz_xxx_xyyy_1,  \
                             ta2_yz_xxx_xyyz_0,  \
                             ta2_yz_xxx_xyyz_1,  \
                             ta2_yz_xxx_xyz_0,   \
                             ta2_yz_xxx_xyz_1,   \
                             ta2_yz_xxx_xyzz_0,  \
                             ta2_yz_xxx_xyzz_1,  \
                             ta2_yz_xxx_xzz_0,   \
                             ta2_yz_xxx_xzz_1,   \
                             ta2_yz_xxx_xzzz_0,  \
                             ta2_yz_xxx_xzzz_1,  \
                             ta2_yz_xxx_yyy_0,   \
                             ta2_yz_xxx_yyy_1,   \
                             ta2_yz_xxx_yyyy_0,  \
                             ta2_yz_xxx_yyyy_1,  \
                             ta2_yz_xxx_yyyz_0,  \
                             ta2_yz_xxx_yyyz_1,  \
                             ta2_yz_xxx_yyz_0,   \
                             ta2_yz_xxx_yyz_1,   \
                             ta2_yz_xxx_yyzz_0,  \
                             ta2_yz_xxx_yyzz_1,  \
                             ta2_yz_xxx_yzz_0,   \
                             ta2_yz_xxx_yzz_1,   \
                             ta2_yz_xxx_yzzz_0,  \
                             ta2_yz_xxx_yzzz_1,  \
                             ta2_yz_xxx_zzz_0,   \
                             ta2_yz_xxx_zzz_1,   \
                             ta2_yz_xxx_zzzz_0,  \
                             ta2_yz_xxx_zzzz_1,  \
                             ta2_yz_xxxx_xxxx_0, \
                             ta2_yz_xxxx_xxxy_0, \
                             ta2_yz_xxxx_xxxz_0, \
                             ta2_yz_xxxx_xxyy_0, \
                             ta2_yz_xxxx_xxyz_0, \
                             ta2_yz_xxxx_xxzz_0, \
                             ta2_yz_xxxx_xyyy_0, \
                             ta2_yz_xxxx_xyyz_0, \
                             ta2_yz_xxxx_xyzz_0, \
                             ta2_yz_xxxx_xzzz_0, \
                             ta2_yz_xxxx_yyyy_0, \
                             ta2_yz_xxxx_yyyz_0, \
                             ta2_yz_xxxx_yyzz_0, \
                             ta2_yz_xxxx_yzzz_0, \
                             ta2_yz_xxxx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxxx_xxxx_0[i] = 3.0 * ta2_yz_xx_xxxx_0[i] * fe_0 - 3.0 * ta2_yz_xx_xxxx_1[i] * fe_0 + 4.0 * ta2_yz_xxx_xxx_0[i] * fe_0 -
                                4.0 * ta2_yz_xxx_xxx_1[i] * fe_0 + ta2_yz_xxx_xxxx_0[i] * pa_x[i] - ta2_yz_xxx_xxxx_1[i] * pc_x[i];

        ta2_yz_xxxx_xxxy_0[i] = 3.0 * ta2_yz_xx_xxxy_0[i] * fe_0 - 3.0 * ta2_yz_xx_xxxy_1[i] * fe_0 + 3.0 * ta2_yz_xxx_xxy_0[i] * fe_0 -
                                3.0 * ta2_yz_xxx_xxy_1[i] * fe_0 + ta2_yz_xxx_xxxy_0[i] * pa_x[i] - ta2_yz_xxx_xxxy_1[i] * pc_x[i];

        ta2_yz_xxxx_xxxz_0[i] = 3.0 * ta2_yz_xx_xxxz_0[i] * fe_0 - 3.0 * ta2_yz_xx_xxxz_1[i] * fe_0 + 3.0 * ta2_yz_xxx_xxz_0[i] * fe_0 -
                                3.0 * ta2_yz_xxx_xxz_1[i] * fe_0 + ta2_yz_xxx_xxxz_0[i] * pa_x[i] - ta2_yz_xxx_xxxz_1[i] * pc_x[i];

        ta2_yz_xxxx_xxyy_0[i] = 3.0 * ta2_yz_xx_xxyy_0[i] * fe_0 - 3.0 * ta2_yz_xx_xxyy_1[i] * fe_0 + 2.0 * ta2_yz_xxx_xyy_0[i] * fe_0 -
                                2.0 * ta2_yz_xxx_xyy_1[i] * fe_0 + ta2_yz_xxx_xxyy_0[i] * pa_x[i] - ta2_yz_xxx_xxyy_1[i] * pc_x[i];

        ta2_yz_xxxx_xxyz_0[i] = 3.0 * ta2_yz_xx_xxyz_0[i] * fe_0 - 3.0 * ta2_yz_xx_xxyz_1[i] * fe_0 + 2.0 * ta2_yz_xxx_xyz_0[i] * fe_0 -
                                2.0 * ta2_yz_xxx_xyz_1[i] * fe_0 + ta2_yz_xxx_xxyz_0[i] * pa_x[i] - ta2_yz_xxx_xxyz_1[i] * pc_x[i];

        ta2_yz_xxxx_xxzz_0[i] = 3.0 * ta2_yz_xx_xxzz_0[i] * fe_0 - 3.0 * ta2_yz_xx_xxzz_1[i] * fe_0 + 2.0 * ta2_yz_xxx_xzz_0[i] * fe_0 -
                                2.0 * ta2_yz_xxx_xzz_1[i] * fe_0 + ta2_yz_xxx_xxzz_0[i] * pa_x[i] - ta2_yz_xxx_xxzz_1[i] * pc_x[i];

        ta2_yz_xxxx_xyyy_0[i] = 3.0 * ta2_yz_xx_xyyy_0[i] * fe_0 - 3.0 * ta2_yz_xx_xyyy_1[i] * fe_0 + ta2_yz_xxx_yyy_0[i] * fe_0 -
                                ta2_yz_xxx_yyy_1[i] * fe_0 + ta2_yz_xxx_xyyy_0[i] * pa_x[i] - ta2_yz_xxx_xyyy_1[i] * pc_x[i];

        ta2_yz_xxxx_xyyz_0[i] = 3.0 * ta2_yz_xx_xyyz_0[i] * fe_0 - 3.0 * ta2_yz_xx_xyyz_1[i] * fe_0 + ta2_yz_xxx_yyz_0[i] * fe_0 -
                                ta2_yz_xxx_yyz_1[i] * fe_0 + ta2_yz_xxx_xyyz_0[i] * pa_x[i] - ta2_yz_xxx_xyyz_1[i] * pc_x[i];

        ta2_yz_xxxx_xyzz_0[i] = 3.0 * ta2_yz_xx_xyzz_0[i] * fe_0 - 3.0 * ta2_yz_xx_xyzz_1[i] * fe_0 + ta2_yz_xxx_yzz_0[i] * fe_0 -
                                ta2_yz_xxx_yzz_1[i] * fe_0 + ta2_yz_xxx_xyzz_0[i] * pa_x[i] - ta2_yz_xxx_xyzz_1[i] * pc_x[i];

        ta2_yz_xxxx_xzzz_0[i] = 3.0 * ta2_yz_xx_xzzz_0[i] * fe_0 - 3.0 * ta2_yz_xx_xzzz_1[i] * fe_0 + ta2_yz_xxx_zzz_0[i] * fe_0 -
                                ta2_yz_xxx_zzz_1[i] * fe_0 + ta2_yz_xxx_xzzz_0[i] * pa_x[i] - ta2_yz_xxx_xzzz_1[i] * pc_x[i];

        ta2_yz_xxxx_yyyy_0[i] =
            3.0 * ta2_yz_xx_yyyy_0[i] * fe_0 - 3.0 * ta2_yz_xx_yyyy_1[i] * fe_0 + ta2_yz_xxx_yyyy_0[i] * pa_x[i] - ta2_yz_xxx_yyyy_1[i] * pc_x[i];

        ta2_yz_xxxx_yyyz_0[i] =
            3.0 * ta2_yz_xx_yyyz_0[i] * fe_0 - 3.0 * ta2_yz_xx_yyyz_1[i] * fe_0 + ta2_yz_xxx_yyyz_0[i] * pa_x[i] - ta2_yz_xxx_yyyz_1[i] * pc_x[i];

        ta2_yz_xxxx_yyzz_0[i] =
            3.0 * ta2_yz_xx_yyzz_0[i] * fe_0 - 3.0 * ta2_yz_xx_yyzz_1[i] * fe_0 + ta2_yz_xxx_yyzz_0[i] * pa_x[i] - ta2_yz_xxx_yyzz_1[i] * pc_x[i];

        ta2_yz_xxxx_yzzz_0[i] =
            3.0 * ta2_yz_xx_yzzz_0[i] * fe_0 - 3.0 * ta2_yz_xx_yzzz_1[i] * fe_0 + ta2_yz_xxx_yzzz_0[i] * pa_x[i] - ta2_yz_xxx_yzzz_1[i] * pc_x[i];

        ta2_yz_xxxx_zzzz_0[i] =
            3.0 * ta2_yz_xx_zzzz_0[i] * fe_0 - 3.0 * ta2_yz_xx_zzzz_1[i] * fe_0 + ta2_yz_xxx_zzzz_0[i] * pa_x[i] - ta2_yz_xxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 915-930 components of targeted buffer : GG

    auto ta2_yz_xxxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 915);

    auto ta2_yz_xxxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 916);

    auto ta2_yz_xxxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 917);

    auto ta2_yz_xxxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 918);

    auto ta2_yz_xxxy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 919);

    auto ta2_yz_xxxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 920);

    auto ta2_yz_xxxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 921);

    auto ta2_yz_xxxy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 922);

    auto ta2_yz_xxxy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 923);

    auto ta2_yz_xxxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 924);

    auto ta2_yz_xxxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 925);

    auto ta2_yz_xxxy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 926);

    auto ta2_yz_xxxy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 927);

    auto ta2_yz_xxxy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 928);

    auto ta2_yz_xxxy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 929);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xxx_xxxx_1,   \
                             ta1_z_xxx_xxxy_1,   \
                             ta1_z_xxx_xxxz_1,   \
                             ta1_z_xxx_xxyy_1,   \
                             ta1_z_xxx_xxyz_1,   \
                             ta1_z_xxx_xxzz_1,   \
                             ta1_z_xxx_xyyy_1,   \
                             ta1_z_xxx_xyyz_1,   \
                             ta1_z_xxx_xyzz_1,   \
                             ta1_z_xxx_xzzz_1,   \
                             ta1_z_xxx_zzzz_1,   \
                             ta2_yz_xxx_xxx_0,   \
                             ta2_yz_xxx_xxx_1,   \
                             ta2_yz_xxx_xxxx_0,  \
                             ta2_yz_xxx_xxxx_1,  \
                             ta2_yz_xxx_xxxy_0,  \
                             ta2_yz_xxx_xxxy_1,  \
                             ta2_yz_xxx_xxxz_0,  \
                             ta2_yz_xxx_xxxz_1,  \
                             ta2_yz_xxx_xxy_0,   \
                             ta2_yz_xxx_xxy_1,   \
                             ta2_yz_xxx_xxyy_0,  \
                             ta2_yz_xxx_xxyy_1,  \
                             ta2_yz_xxx_xxyz_0,  \
                             ta2_yz_xxx_xxyz_1,  \
                             ta2_yz_xxx_xxz_0,   \
                             ta2_yz_xxx_xxz_1,   \
                             ta2_yz_xxx_xxzz_0,  \
                             ta2_yz_xxx_xxzz_1,  \
                             ta2_yz_xxx_xyy_0,   \
                             ta2_yz_xxx_xyy_1,   \
                             ta2_yz_xxx_xyyy_0,  \
                             ta2_yz_xxx_xyyy_1,  \
                             ta2_yz_xxx_xyyz_0,  \
                             ta2_yz_xxx_xyyz_1,  \
                             ta2_yz_xxx_xyz_0,   \
                             ta2_yz_xxx_xyz_1,   \
                             ta2_yz_xxx_xyzz_0,  \
                             ta2_yz_xxx_xyzz_1,  \
                             ta2_yz_xxx_xzz_0,   \
                             ta2_yz_xxx_xzz_1,   \
                             ta2_yz_xxx_xzzz_0,  \
                             ta2_yz_xxx_xzzz_1,  \
                             ta2_yz_xxx_zzzz_0,  \
                             ta2_yz_xxx_zzzz_1,  \
                             ta2_yz_xxxy_xxxx_0, \
                             ta2_yz_xxxy_xxxy_0, \
                             ta2_yz_xxxy_xxxz_0, \
                             ta2_yz_xxxy_xxyy_0, \
                             ta2_yz_xxxy_xxyz_0, \
                             ta2_yz_xxxy_xxzz_0, \
                             ta2_yz_xxxy_xyyy_0, \
                             ta2_yz_xxxy_xyyz_0, \
                             ta2_yz_xxxy_xyzz_0, \
                             ta2_yz_xxxy_xzzz_0, \
                             ta2_yz_xxxy_yyyy_0, \
                             ta2_yz_xxxy_yyyz_0, \
                             ta2_yz_xxxy_yyzz_0, \
                             ta2_yz_xxxy_yzzz_0, \
                             ta2_yz_xxxy_zzzz_0, \
                             ta2_yz_xxy_yyyy_0,  \
                             ta2_yz_xxy_yyyy_1,  \
                             ta2_yz_xxy_yyyz_0,  \
                             ta2_yz_xxy_yyyz_1,  \
                             ta2_yz_xxy_yyzz_0,  \
                             ta2_yz_xxy_yyzz_1,  \
                             ta2_yz_xxy_yzzz_0,  \
                             ta2_yz_xxy_yzzz_1,  \
                             ta2_yz_xy_yyyy_0,   \
                             ta2_yz_xy_yyyy_1,   \
                             ta2_yz_xy_yyyz_0,   \
                             ta2_yz_xy_yyyz_1,   \
                             ta2_yz_xy_yyzz_0,   \
                             ta2_yz_xy_yyzz_1,   \
                             ta2_yz_xy_yzzz_0,   \
                             ta2_yz_xy_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxxy_xxxx_0[i] = ta1_z_xxx_xxxx_1[i] + ta2_yz_xxx_xxxx_0[i] * pa_y[i] - ta2_yz_xxx_xxxx_1[i] * pc_y[i];

        ta2_yz_xxxy_xxxy_0[i] = ta2_yz_xxx_xxx_0[i] * fe_0 - ta2_yz_xxx_xxx_1[i] * fe_0 + ta1_z_xxx_xxxy_1[i] + ta2_yz_xxx_xxxy_0[i] * pa_y[i] -
                                ta2_yz_xxx_xxxy_1[i] * pc_y[i];

        ta2_yz_xxxy_xxxz_0[i] = ta1_z_xxx_xxxz_1[i] + ta2_yz_xxx_xxxz_0[i] * pa_y[i] - ta2_yz_xxx_xxxz_1[i] * pc_y[i];

        ta2_yz_xxxy_xxyy_0[i] = 2.0 * ta2_yz_xxx_xxy_0[i] * fe_0 - 2.0 * ta2_yz_xxx_xxy_1[i] * fe_0 + ta1_z_xxx_xxyy_1[i] +
                                ta2_yz_xxx_xxyy_0[i] * pa_y[i] - ta2_yz_xxx_xxyy_1[i] * pc_y[i];

        ta2_yz_xxxy_xxyz_0[i] = ta2_yz_xxx_xxz_0[i] * fe_0 - ta2_yz_xxx_xxz_1[i] * fe_0 + ta1_z_xxx_xxyz_1[i] + ta2_yz_xxx_xxyz_0[i] * pa_y[i] -
                                ta2_yz_xxx_xxyz_1[i] * pc_y[i];

        ta2_yz_xxxy_xxzz_0[i] = ta1_z_xxx_xxzz_1[i] + ta2_yz_xxx_xxzz_0[i] * pa_y[i] - ta2_yz_xxx_xxzz_1[i] * pc_y[i];

        ta2_yz_xxxy_xyyy_0[i] = 3.0 * ta2_yz_xxx_xyy_0[i] * fe_0 - 3.0 * ta2_yz_xxx_xyy_1[i] * fe_0 + ta1_z_xxx_xyyy_1[i] +
                                ta2_yz_xxx_xyyy_0[i] * pa_y[i] - ta2_yz_xxx_xyyy_1[i] * pc_y[i];

        ta2_yz_xxxy_xyyz_0[i] = 2.0 * ta2_yz_xxx_xyz_0[i] * fe_0 - 2.0 * ta2_yz_xxx_xyz_1[i] * fe_0 + ta1_z_xxx_xyyz_1[i] +
                                ta2_yz_xxx_xyyz_0[i] * pa_y[i] - ta2_yz_xxx_xyyz_1[i] * pc_y[i];

        ta2_yz_xxxy_xyzz_0[i] = ta2_yz_xxx_xzz_0[i] * fe_0 - ta2_yz_xxx_xzz_1[i] * fe_0 + ta1_z_xxx_xyzz_1[i] + ta2_yz_xxx_xyzz_0[i] * pa_y[i] -
                                ta2_yz_xxx_xyzz_1[i] * pc_y[i];

        ta2_yz_xxxy_xzzz_0[i] = ta1_z_xxx_xzzz_1[i] + ta2_yz_xxx_xzzz_0[i] * pa_y[i] - ta2_yz_xxx_xzzz_1[i] * pc_y[i];

        ta2_yz_xxxy_yyyy_0[i] =
            2.0 * ta2_yz_xy_yyyy_0[i] * fe_0 - 2.0 * ta2_yz_xy_yyyy_1[i] * fe_0 + ta2_yz_xxy_yyyy_0[i] * pa_x[i] - ta2_yz_xxy_yyyy_1[i] * pc_x[i];

        ta2_yz_xxxy_yyyz_0[i] =
            2.0 * ta2_yz_xy_yyyz_0[i] * fe_0 - 2.0 * ta2_yz_xy_yyyz_1[i] * fe_0 + ta2_yz_xxy_yyyz_0[i] * pa_x[i] - ta2_yz_xxy_yyyz_1[i] * pc_x[i];

        ta2_yz_xxxy_yyzz_0[i] =
            2.0 * ta2_yz_xy_yyzz_0[i] * fe_0 - 2.0 * ta2_yz_xy_yyzz_1[i] * fe_0 + ta2_yz_xxy_yyzz_0[i] * pa_x[i] - ta2_yz_xxy_yyzz_1[i] * pc_x[i];

        ta2_yz_xxxy_yzzz_0[i] =
            2.0 * ta2_yz_xy_yzzz_0[i] * fe_0 - 2.0 * ta2_yz_xy_yzzz_1[i] * fe_0 + ta2_yz_xxy_yzzz_0[i] * pa_x[i] - ta2_yz_xxy_yzzz_1[i] * pc_x[i];

        ta2_yz_xxxy_zzzz_0[i] = ta1_z_xxx_zzzz_1[i] + ta2_yz_xxx_zzzz_0[i] * pa_y[i] - ta2_yz_xxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 930-945 components of targeted buffer : GG

    auto ta2_yz_xxxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 930);

    auto ta2_yz_xxxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 931);

    auto ta2_yz_xxxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 932);

    auto ta2_yz_xxxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 933);

    auto ta2_yz_xxxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 934);

    auto ta2_yz_xxxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 935);

    auto ta2_yz_xxxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 936);

    auto ta2_yz_xxxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 937);

    auto ta2_yz_xxxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 938);

    auto ta2_yz_xxxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 939);

    auto ta2_yz_xxxz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 940);

    auto ta2_yz_xxxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 941);

    auto ta2_yz_xxxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 942);

    auto ta2_yz_xxxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 943);

    auto ta2_yz_xxxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 944);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xxx_xxxx_1,   \
                             ta1_y_xxx_xxxy_1,   \
                             ta1_y_xxx_xxxz_1,   \
                             ta1_y_xxx_xxyy_1,   \
                             ta1_y_xxx_xxyz_1,   \
                             ta1_y_xxx_xxzz_1,   \
                             ta1_y_xxx_xyyy_1,   \
                             ta1_y_xxx_xyyz_1,   \
                             ta1_y_xxx_xyzz_1,   \
                             ta1_y_xxx_xzzz_1,   \
                             ta1_y_xxx_yyyy_1,   \
                             ta2_yz_xxx_xxx_0,   \
                             ta2_yz_xxx_xxx_1,   \
                             ta2_yz_xxx_xxxx_0,  \
                             ta2_yz_xxx_xxxx_1,  \
                             ta2_yz_xxx_xxxy_0,  \
                             ta2_yz_xxx_xxxy_1,  \
                             ta2_yz_xxx_xxxz_0,  \
                             ta2_yz_xxx_xxxz_1,  \
                             ta2_yz_xxx_xxy_0,   \
                             ta2_yz_xxx_xxy_1,   \
                             ta2_yz_xxx_xxyy_0,  \
                             ta2_yz_xxx_xxyy_1,  \
                             ta2_yz_xxx_xxyz_0,  \
                             ta2_yz_xxx_xxyz_1,  \
                             ta2_yz_xxx_xxz_0,   \
                             ta2_yz_xxx_xxz_1,   \
                             ta2_yz_xxx_xxzz_0,  \
                             ta2_yz_xxx_xxzz_1,  \
                             ta2_yz_xxx_xyy_0,   \
                             ta2_yz_xxx_xyy_1,   \
                             ta2_yz_xxx_xyyy_0,  \
                             ta2_yz_xxx_xyyy_1,  \
                             ta2_yz_xxx_xyyz_0,  \
                             ta2_yz_xxx_xyyz_1,  \
                             ta2_yz_xxx_xyz_0,   \
                             ta2_yz_xxx_xyz_1,   \
                             ta2_yz_xxx_xyzz_0,  \
                             ta2_yz_xxx_xyzz_1,  \
                             ta2_yz_xxx_xzz_0,   \
                             ta2_yz_xxx_xzz_1,   \
                             ta2_yz_xxx_xzzz_0,  \
                             ta2_yz_xxx_xzzz_1,  \
                             ta2_yz_xxx_yyyy_0,  \
                             ta2_yz_xxx_yyyy_1,  \
                             ta2_yz_xxxz_xxxx_0, \
                             ta2_yz_xxxz_xxxy_0, \
                             ta2_yz_xxxz_xxxz_0, \
                             ta2_yz_xxxz_xxyy_0, \
                             ta2_yz_xxxz_xxyz_0, \
                             ta2_yz_xxxz_xxzz_0, \
                             ta2_yz_xxxz_xyyy_0, \
                             ta2_yz_xxxz_xyyz_0, \
                             ta2_yz_xxxz_xyzz_0, \
                             ta2_yz_xxxz_xzzz_0, \
                             ta2_yz_xxxz_yyyy_0, \
                             ta2_yz_xxxz_yyyz_0, \
                             ta2_yz_xxxz_yyzz_0, \
                             ta2_yz_xxxz_yzzz_0, \
                             ta2_yz_xxxz_zzzz_0, \
                             ta2_yz_xxz_yyyz_0,  \
                             ta2_yz_xxz_yyyz_1,  \
                             ta2_yz_xxz_yyzz_0,  \
                             ta2_yz_xxz_yyzz_1,  \
                             ta2_yz_xxz_yzzz_0,  \
                             ta2_yz_xxz_yzzz_1,  \
                             ta2_yz_xxz_zzzz_0,  \
                             ta2_yz_xxz_zzzz_1,  \
                             ta2_yz_xz_yyyz_0,   \
                             ta2_yz_xz_yyyz_1,   \
                             ta2_yz_xz_yyzz_0,   \
                             ta2_yz_xz_yyzz_1,   \
                             ta2_yz_xz_yzzz_0,   \
                             ta2_yz_xz_yzzz_1,   \
                             ta2_yz_xz_zzzz_0,   \
                             ta2_yz_xz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxxz_xxxx_0[i] = ta1_y_xxx_xxxx_1[i] + ta2_yz_xxx_xxxx_0[i] * pa_z[i] - ta2_yz_xxx_xxxx_1[i] * pc_z[i];

        ta2_yz_xxxz_xxxy_0[i] = ta1_y_xxx_xxxy_1[i] + ta2_yz_xxx_xxxy_0[i] * pa_z[i] - ta2_yz_xxx_xxxy_1[i] * pc_z[i];

        ta2_yz_xxxz_xxxz_0[i] = ta2_yz_xxx_xxx_0[i] * fe_0 - ta2_yz_xxx_xxx_1[i] * fe_0 + ta1_y_xxx_xxxz_1[i] + ta2_yz_xxx_xxxz_0[i] * pa_z[i] -
                                ta2_yz_xxx_xxxz_1[i] * pc_z[i];

        ta2_yz_xxxz_xxyy_0[i] = ta1_y_xxx_xxyy_1[i] + ta2_yz_xxx_xxyy_0[i] * pa_z[i] - ta2_yz_xxx_xxyy_1[i] * pc_z[i];

        ta2_yz_xxxz_xxyz_0[i] = ta2_yz_xxx_xxy_0[i] * fe_0 - ta2_yz_xxx_xxy_1[i] * fe_0 + ta1_y_xxx_xxyz_1[i] + ta2_yz_xxx_xxyz_0[i] * pa_z[i] -
                                ta2_yz_xxx_xxyz_1[i] * pc_z[i];

        ta2_yz_xxxz_xxzz_0[i] = 2.0 * ta2_yz_xxx_xxz_0[i] * fe_0 - 2.0 * ta2_yz_xxx_xxz_1[i] * fe_0 + ta1_y_xxx_xxzz_1[i] +
                                ta2_yz_xxx_xxzz_0[i] * pa_z[i] - ta2_yz_xxx_xxzz_1[i] * pc_z[i];

        ta2_yz_xxxz_xyyy_0[i] = ta1_y_xxx_xyyy_1[i] + ta2_yz_xxx_xyyy_0[i] * pa_z[i] - ta2_yz_xxx_xyyy_1[i] * pc_z[i];

        ta2_yz_xxxz_xyyz_0[i] = ta2_yz_xxx_xyy_0[i] * fe_0 - ta2_yz_xxx_xyy_1[i] * fe_0 + ta1_y_xxx_xyyz_1[i] + ta2_yz_xxx_xyyz_0[i] * pa_z[i] -
                                ta2_yz_xxx_xyyz_1[i] * pc_z[i];

        ta2_yz_xxxz_xyzz_0[i] = 2.0 * ta2_yz_xxx_xyz_0[i] * fe_0 - 2.0 * ta2_yz_xxx_xyz_1[i] * fe_0 + ta1_y_xxx_xyzz_1[i] +
                                ta2_yz_xxx_xyzz_0[i] * pa_z[i] - ta2_yz_xxx_xyzz_1[i] * pc_z[i];

        ta2_yz_xxxz_xzzz_0[i] = 3.0 * ta2_yz_xxx_xzz_0[i] * fe_0 - 3.0 * ta2_yz_xxx_xzz_1[i] * fe_0 + ta1_y_xxx_xzzz_1[i] +
                                ta2_yz_xxx_xzzz_0[i] * pa_z[i] - ta2_yz_xxx_xzzz_1[i] * pc_z[i];

        ta2_yz_xxxz_yyyy_0[i] = ta1_y_xxx_yyyy_1[i] + ta2_yz_xxx_yyyy_0[i] * pa_z[i] - ta2_yz_xxx_yyyy_1[i] * pc_z[i];

        ta2_yz_xxxz_yyyz_0[i] =
            2.0 * ta2_yz_xz_yyyz_0[i] * fe_0 - 2.0 * ta2_yz_xz_yyyz_1[i] * fe_0 + ta2_yz_xxz_yyyz_0[i] * pa_x[i] - ta2_yz_xxz_yyyz_1[i] * pc_x[i];

        ta2_yz_xxxz_yyzz_0[i] =
            2.0 * ta2_yz_xz_yyzz_0[i] * fe_0 - 2.0 * ta2_yz_xz_yyzz_1[i] * fe_0 + ta2_yz_xxz_yyzz_0[i] * pa_x[i] - ta2_yz_xxz_yyzz_1[i] * pc_x[i];

        ta2_yz_xxxz_yzzz_0[i] =
            2.0 * ta2_yz_xz_yzzz_0[i] * fe_0 - 2.0 * ta2_yz_xz_yzzz_1[i] * fe_0 + ta2_yz_xxz_yzzz_0[i] * pa_x[i] - ta2_yz_xxz_yzzz_1[i] * pc_x[i];

        ta2_yz_xxxz_zzzz_0[i] =
            2.0 * ta2_yz_xz_zzzz_0[i] * fe_0 - 2.0 * ta2_yz_xz_zzzz_1[i] * fe_0 + ta2_yz_xxz_zzzz_0[i] * pa_x[i] - ta2_yz_xxz_zzzz_1[i] * pc_x[i];
    }

    // Set up 945-960 components of targeted buffer : GG

    auto ta2_yz_xxyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 945);

    auto ta2_yz_xxyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 946);

    auto ta2_yz_xxyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 947);

    auto ta2_yz_xxyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 948);

    auto ta2_yz_xxyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 949);

    auto ta2_yz_xxyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 950);

    auto ta2_yz_xxyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 951);

    auto ta2_yz_xxyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 952);

    auto ta2_yz_xxyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 953);

    auto ta2_yz_xxyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 954);

    auto ta2_yz_xxyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 955);

    auto ta2_yz_xxyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 956);

    auto ta2_yz_xxyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 957);

    auto ta2_yz_xxyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 958);

    auto ta2_yz_xxyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 959);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xxy_xxxx_1,   \
                             ta1_z_xxy_xxxz_1,   \
                             ta1_z_xxy_xxzz_1,   \
                             ta1_z_xxy_xzzz_1,   \
                             ta2_yz_xx_xxxx_0,   \
                             ta2_yz_xx_xxxx_1,   \
                             ta2_yz_xx_xxxz_0,   \
                             ta2_yz_xx_xxxz_1,   \
                             ta2_yz_xx_xxzz_0,   \
                             ta2_yz_xx_xxzz_1,   \
                             ta2_yz_xx_xzzz_0,   \
                             ta2_yz_xx_xzzz_1,   \
                             ta2_yz_xxy_xxxx_0,  \
                             ta2_yz_xxy_xxxx_1,  \
                             ta2_yz_xxy_xxxz_0,  \
                             ta2_yz_xxy_xxxz_1,  \
                             ta2_yz_xxy_xxzz_0,  \
                             ta2_yz_xxy_xxzz_1,  \
                             ta2_yz_xxy_xzzz_0,  \
                             ta2_yz_xxy_xzzz_1,  \
                             ta2_yz_xxyy_xxxx_0, \
                             ta2_yz_xxyy_xxxy_0, \
                             ta2_yz_xxyy_xxxz_0, \
                             ta2_yz_xxyy_xxyy_0, \
                             ta2_yz_xxyy_xxyz_0, \
                             ta2_yz_xxyy_xxzz_0, \
                             ta2_yz_xxyy_xyyy_0, \
                             ta2_yz_xxyy_xyyz_0, \
                             ta2_yz_xxyy_xyzz_0, \
                             ta2_yz_xxyy_xzzz_0, \
                             ta2_yz_xxyy_yyyy_0, \
                             ta2_yz_xxyy_yyyz_0, \
                             ta2_yz_xxyy_yyzz_0, \
                             ta2_yz_xxyy_yzzz_0, \
                             ta2_yz_xxyy_zzzz_0, \
                             ta2_yz_xyy_xxxy_0,  \
                             ta2_yz_xyy_xxxy_1,  \
                             ta2_yz_xyy_xxy_0,   \
                             ta2_yz_xyy_xxy_1,   \
                             ta2_yz_xyy_xxyy_0,  \
                             ta2_yz_xyy_xxyy_1,  \
                             ta2_yz_xyy_xxyz_0,  \
                             ta2_yz_xyy_xxyz_1,  \
                             ta2_yz_xyy_xyy_0,   \
                             ta2_yz_xyy_xyy_1,   \
                             ta2_yz_xyy_xyyy_0,  \
                             ta2_yz_xyy_xyyy_1,  \
                             ta2_yz_xyy_xyyz_0,  \
                             ta2_yz_xyy_xyyz_1,  \
                             ta2_yz_xyy_xyz_0,   \
                             ta2_yz_xyy_xyz_1,   \
                             ta2_yz_xyy_xyzz_0,  \
                             ta2_yz_xyy_xyzz_1,  \
                             ta2_yz_xyy_yyy_0,   \
                             ta2_yz_xyy_yyy_1,   \
                             ta2_yz_xyy_yyyy_0,  \
                             ta2_yz_xyy_yyyy_1,  \
                             ta2_yz_xyy_yyyz_0,  \
                             ta2_yz_xyy_yyyz_1,  \
                             ta2_yz_xyy_yyz_0,   \
                             ta2_yz_xyy_yyz_1,   \
                             ta2_yz_xyy_yyzz_0,  \
                             ta2_yz_xyy_yyzz_1,  \
                             ta2_yz_xyy_yzz_0,   \
                             ta2_yz_xyy_yzz_1,   \
                             ta2_yz_xyy_yzzz_0,  \
                             ta2_yz_xyy_yzzz_1,  \
                             ta2_yz_xyy_zzzz_0,  \
                             ta2_yz_xyy_zzzz_1,  \
                             ta2_yz_yy_xxxy_0,   \
                             ta2_yz_yy_xxxy_1,   \
                             ta2_yz_yy_xxyy_0,   \
                             ta2_yz_yy_xxyy_1,   \
                             ta2_yz_yy_xxyz_0,   \
                             ta2_yz_yy_xxyz_1,   \
                             ta2_yz_yy_xyyy_0,   \
                             ta2_yz_yy_xyyy_1,   \
                             ta2_yz_yy_xyyz_0,   \
                             ta2_yz_yy_xyyz_1,   \
                             ta2_yz_yy_xyzz_0,   \
                             ta2_yz_yy_xyzz_1,   \
                             ta2_yz_yy_yyyy_0,   \
                             ta2_yz_yy_yyyy_1,   \
                             ta2_yz_yy_yyyz_0,   \
                             ta2_yz_yy_yyyz_1,   \
                             ta2_yz_yy_yyzz_0,   \
                             ta2_yz_yy_yyzz_1,   \
                             ta2_yz_yy_yzzz_0,   \
                             ta2_yz_yy_yzzz_1,   \
                             ta2_yz_yy_zzzz_0,   \
                             ta2_yz_yy_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxyy_xxxx_0[i] = ta2_yz_xx_xxxx_0[i] * fe_0 - ta2_yz_xx_xxxx_1[i] * fe_0 + ta1_z_xxy_xxxx_1[i] + ta2_yz_xxy_xxxx_0[i] * pa_y[i] -
                                ta2_yz_xxy_xxxx_1[i] * pc_y[i];

        ta2_yz_xxyy_xxxy_0[i] = ta2_yz_yy_xxxy_0[i] * fe_0 - ta2_yz_yy_xxxy_1[i] * fe_0 + 3.0 * ta2_yz_xyy_xxy_0[i] * fe_0 -
                                3.0 * ta2_yz_xyy_xxy_1[i] * fe_0 + ta2_yz_xyy_xxxy_0[i] * pa_x[i] - ta2_yz_xyy_xxxy_1[i] * pc_x[i];

        ta2_yz_xxyy_xxxz_0[i] = ta2_yz_xx_xxxz_0[i] * fe_0 - ta2_yz_xx_xxxz_1[i] * fe_0 + ta1_z_xxy_xxxz_1[i] + ta2_yz_xxy_xxxz_0[i] * pa_y[i] -
                                ta2_yz_xxy_xxxz_1[i] * pc_y[i];

        ta2_yz_xxyy_xxyy_0[i] = ta2_yz_yy_xxyy_0[i] * fe_0 - ta2_yz_yy_xxyy_1[i] * fe_0 + 2.0 * ta2_yz_xyy_xyy_0[i] * fe_0 -
                                2.0 * ta2_yz_xyy_xyy_1[i] * fe_0 + ta2_yz_xyy_xxyy_0[i] * pa_x[i] - ta2_yz_xyy_xxyy_1[i] * pc_x[i];

        ta2_yz_xxyy_xxyz_0[i] = ta2_yz_yy_xxyz_0[i] * fe_0 - ta2_yz_yy_xxyz_1[i] * fe_0 + 2.0 * ta2_yz_xyy_xyz_0[i] * fe_0 -
                                2.0 * ta2_yz_xyy_xyz_1[i] * fe_0 + ta2_yz_xyy_xxyz_0[i] * pa_x[i] - ta2_yz_xyy_xxyz_1[i] * pc_x[i];

        ta2_yz_xxyy_xxzz_0[i] = ta2_yz_xx_xxzz_0[i] * fe_0 - ta2_yz_xx_xxzz_1[i] * fe_0 + ta1_z_xxy_xxzz_1[i] + ta2_yz_xxy_xxzz_0[i] * pa_y[i] -
                                ta2_yz_xxy_xxzz_1[i] * pc_y[i];

        ta2_yz_xxyy_xyyy_0[i] = ta2_yz_yy_xyyy_0[i] * fe_0 - ta2_yz_yy_xyyy_1[i] * fe_0 + ta2_yz_xyy_yyy_0[i] * fe_0 - ta2_yz_xyy_yyy_1[i] * fe_0 +
                                ta2_yz_xyy_xyyy_0[i] * pa_x[i] - ta2_yz_xyy_xyyy_1[i] * pc_x[i];

        ta2_yz_xxyy_xyyz_0[i] = ta2_yz_yy_xyyz_0[i] * fe_0 - ta2_yz_yy_xyyz_1[i] * fe_0 + ta2_yz_xyy_yyz_0[i] * fe_0 - ta2_yz_xyy_yyz_1[i] * fe_0 +
                                ta2_yz_xyy_xyyz_0[i] * pa_x[i] - ta2_yz_xyy_xyyz_1[i] * pc_x[i];

        ta2_yz_xxyy_xyzz_0[i] = ta2_yz_yy_xyzz_0[i] * fe_0 - ta2_yz_yy_xyzz_1[i] * fe_0 + ta2_yz_xyy_yzz_0[i] * fe_0 - ta2_yz_xyy_yzz_1[i] * fe_0 +
                                ta2_yz_xyy_xyzz_0[i] * pa_x[i] - ta2_yz_xyy_xyzz_1[i] * pc_x[i];

        ta2_yz_xxyy_xzzz_0[i] = ta2_yz_xx_xzzz_0[i] * fe_0 - ta2_yz_xx_xzzz_1[i] * fe_0 + ta1_z_xxy_xzzz_1[i] + ta2_yz_xxy_xzzz_0[i] * pa_y[i] -
                                ta2_yz_xxy_xzzz_1[i] * pc_y[i];

        ta2_yz_xxyy_yyyy_0[i] =
            ta2_yz_yy_yyyy_0[i] * fe_0 - ta2_yz_yy_yyyy_1[i] * fe_0 + ta2_yz_xyy_yyyy_0[i] * pa_x[i] - ta2_yz_xyy_yyyy_1[i] * pc_x[i];

        ta2_yz_xxyy_yyyz_0[i] =
            ta2_yz_yy_yyyz_0[i] * fe_0 - ta2_yz_yy_yyyz_1[i] * fe_0 + ta2_yz_xyy_yyyz_0[i] * pa_x[i] - ta2_yz_xyy_yyyz_1[i] * pc_x[i];

        ta2_yz_xxyy_yyzz_0[i] =
            ta2_yz_yy_yyzz_0[i] * fe_0 - ta2_yz_yy_yyzz_1[i] * fe_0 + ta2_yz_xyy_yyzz_0[i] * pa_x[i] - ta2_yz_xyy_yyzz_1[i] * pc_x[i];

        ta2_yz_xxyy_yzzz_0[i] =
            ta2_yz_yy_yzzz_0[i] * fe_0 - ta2_yz_yy_yzzz_1[i] * fe_0 + ta2_yz_xyy_yzzz_0[i] * pa_x[i] - ta2_yz_xyy_yzzz_1[i] * pc_x[i];

        ta2_yz_xxyy_zzzz_0[i] =
            ta2_yz_yy_zzzz_0[i] * fe_0 - ta2_yz_yy_zzzz_1[i] * fe_0 + ta2_yz_xyy_zzzz_0[i] * pa_x[i] - ta2_yz_xyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 960-975 components of targeted buffer : GG

    auto ta2_yz_xxyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 960);

    auto ta2_yz_xxyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 961);

    auto ta2_yz_xxyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 962);

    auto ta2_yz_xxyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 963);

    auto ta2_yz_xxyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 964);

    auto ta2_yz_xxyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 965);

    auto ta2_yz_xxyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 966);

    auto ta2_yz_xxyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 967);

    auto ta2_yz_xxyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 968);

    auto ta2_yz_xxyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 969);

    auto ta2_yz_xxyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 970);

    auto ta2_yz_xxyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 971);

    auto ta2_yz_xxyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 972);

    auto ta2_yz_xxyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 973);

    auto ta2_yz_xxyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 974);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_xxy_xxxy_1,   \
                             ta1_y_xxy_xxyy_1,   \
                             ta1_y_xxy_xyyy_1,   \
                             ta1_y_xxy_yyyy_1,   \
                             ta1_z_xxz_xxxx_1,   \
                             ta1_z_xxz_xxxz_1,   \
                             ta1_z_xxz_xxyz_1,   \
                             ta1_z_xxz_xxzz_1,   \
                             ta1_z_xxz_xyyz_1,   \
                             ta1_z_xxz_xyzz_1,   \
                             ta1_z_xxz_xzzz_1,   \
                             ta1_z_xxz_zzzz_1,   \
                             ta2_yz_xxy_xxxy_0,  \
                             ta2_yz_xxy_xxxy_1,  \
                             ta2_yz_xxy_xxyy_0,  \
                             ta2_yz_xxy_xxyy_1,  \
                             ta2_yz_xxy_xyyy_0,  \
                             ta2_yz_xxy_xyyy_1,  \
                             ta2_yz_xxy_yyyy_0,  \
                             ta2_yz_xxy_yyyy_1,  \
                             ta2_yz_xxyz_xxxx_0, \
                             ta2_yz_xxyz_xxxy_0, \
                             ta2_yz_xxyz_xxxz_0, \
                             ta2_yz_xxyz_xxyy_0, \
                             ta2_yz_xxyz_xxyz_0, \
                             ta2_yz_xxyz_xxzz_0, \
                             ta2_yz_xxyz_xyyy_0, \
                             ta2_yz_xxyz_xyyz_0, \
                             ta2_yz_xxyz_xyzz_0, \
                             ta2_yz_xxyz_xzzz_0, \
                             ta2_yz_xxyz_yyyy_0, \
                             ta2_yz_xxyz_yyyz_0, \
                             ta2_yz_xxyz_yyzz_0, \
                             ta2_yz_xxyz_yzzz_0, \
                             ta2_yz_xxyz_zzzz_0, \
                             ta2_yz_xxz_xxxx_0,  \
                             ta2_yz_xxz_xxxx_1,  \
                             ta2_yz_xxz_xxxz_0,  \
                             ta2_yz_xxz_xxxz_1,  \
                             ta2_yz_xxz_xxyz_0,  \
                             ta2_yz_xxz_xxyz_1,  \
                             ta2_yz_xxz_xxz_0,   \
                             ta2_yz_xxz_xxz_1,   \
                             ta2_yz_xxz_xxzz_0,  \
                             ta2_yz_xxz_xxzz_1,  \
                             ta2_yz_xxz_xyyz_0,  \
                             ta2_yz_xxz_xyyz_1,  \
                             ta2_yz_xxz_xyz_0,   \
                             ta2_yz_xxz_xyz_1,   \
                             ta2_yz_xxz_xyzz_0,  \
                             ta2_yz_xxz_xyzz_1,  \
                             ta2_yz_xxz_xzz_0,   \
                             ta2_yz_xxz_xzz_1,   \
                             ta2_yz_xxz_xzzz_0,  \
                             ta2_yz_xxz_xzzz_1,  \
                             ta2_yz_xxz_zzzz_0,  \
                             ta2_yz_xxz_zzzz_1,  \
                             ta2_yz_xyz_yyyz_0,  \
                             ta2_yz_xyz_yyyz_1,  \
                             ta2_yz_xyz_yyzz_0,  \
                             ta2_yz_xyz_yyzz_1,  \
                             ta2_yz_xyz_yzzz_0,  \
                             ta2_yz_xyz_yzzz_1,  \
                             ta2_yz_yz_yyyz_0,   \
                             ta2_yz_yz_yyyz_1,   \
                             ta2_yz_yz_yyzz_0,   \
                             ta2_yz_yz_yyzz_1,   \
                             ta2_yz_yz_yzzz_0,   \
                             ta2_yz_yz_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxyz_xxxx_0[i] = ta1_z_xxz_xxxx_1[i] + ta2_yz_xxz_xxxx_0[i] * pa_y[i] - ta2_yz_xxz_xxxx_1[i] * pc_y[i];

        ta2_yz_xxyz_xxxy_0[i] = ta1_y_xxy_xxxy_1[i] + ta2_yz_xxy_xxxy_0[i] * pa_z[i] - ta2_yz_xxy_xxxy_1[i] * pc_z[i];

        ta2_yz_xxyz_xxxz_0[i] = ta1_z_xxz_xxxz_1[i] + ta2_yz_xxz_xxxz_0[i] * pa_y[i] - ta2_yz_xxz_xxxz_1[i] * pc_y[i];

        ta2_yz_xxyz_xxyy_0[i] = ta1_y_xxy_xxyy_1[i] + ta2_yz_xxy_xxyy_0[i] * pa_z[i] - ta2_yz_xxy_xxyy_1[i] * pc_z[i];

        ta2_yz_xxyz_xxyz_0[i] = ta2_yz_xxz_xxz_0[i] * fe_0 - ta2_yz_xxz_xxz_1[i] * fe_0 + ta1_z_xxz_xxyz_1[i] + ta2_yz_xxz_xxyz_0[i] * pa_y[i] -
                                ta2_yz_xxz_xxyz_1[i] * pc_y[i];

        ta2_yz_xxyz_xxzz_0[i] = ta1_z_xxz_xxzz_1[i] + ta2_yz_xxz_xxzz_0[i] * pa_y[i] - ta2_yz_xxz_xxzz_1[i] * pc_y[i];

        ta2_yz_xxyz_xyyy_0[i] = ta1_y_xxy_xyyy_1[i] + ta2_yz_xxy_xyyy_0[i] * pa_z[i] - ta2_yz_xxy_xyyy_1[i] * pc_z[i];

        ta2_yz_xxyz_xyyz_0[i] = 2.0 * ta2_yz_xxz_xyz_0[i] * fe_0 - 2.0 * ta2_yz_xxz_xyz_1[i] * fe_0 + ta1_z_xxz_xyyz_1[i] +
                                ta2_yz_xxz_xyyz_0[i] * pa_y[i] - ta2_yz_xxz_xyyz_1[i] * pc_y[i];

        ta2_yz_xxyz_xyzz_0[i] = ta2_yz_xxz_xzz_0[i] * fe_0 - ta2_yz_xxz_xzz_1[i] * fe_0 + ta1_z_xxz_xyzz_1[i] + ta2_yz_xxz_xyzz_0[i] * pa_y[i] -
                                ta2_yz_xxz_xyzz_1[i] * pc_y[i];

        ta2_yz_xxyz_xzzz_0[i] = ta1_z_xxz_xzzz_1[i] + ta2_yz_xxz_xzzz_0[i] * pa_y[i] - ta2_yz_xxz_xzzz_1[i] * pc_y[i];

        ta2_yz_xxyz_yyyy_0[i] = ta1_y_xxy_yyyy_1[i] + ta2_yz_xxy_yyyy_0[i] * pa_z[i] - ta2_yz_xxy_yyyy_1[i] * pc_z[i];

        ta2_yz_xxyz_yyyz_0[i] =
            ta2_yz_yz_yyyz_0[i] * fe_0 - ta2_yz_yz_yyyz_1[i] * fe_0 + ta2_yz_xyz_yyyz_0[i] * pa_x[i] - ta2_yz_xyz_yyyz_1[i] * pc_x[i];

        ta2_yz_xxyz_yyzz_0[i] =
            ta2_yz_yz_yyzz_0[i] * fe_0 - ta2_yz_yz_yyzz_1[i] * fe_0 + ta2_yz_xyz_yyzz_0[i] * pa_x[i] - ta2_yz_xyz_yyzz_1[i] * pc_x[i];

        ta2_yz_xxyz_yzzz_0[i] =
            ta2_yz_yz_yzzz_0[i] * fe_0 - ta2_yz_yz_yzzz_1[i] * fe_0 + ta2_yz_xyz_yzzz_0[i] * pa_x[i] - ta2_yz_xyz_yzzz_1[i] * pc_x[i];

        ta2_yz_xxyz_zzzz_0[i] = ta1_z_xxz_zzzz_1[i] + ta2_yz_xxz_zzzz_0[i] * pa_y[i] - ta2_yz_xxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 975-990 components of targeted buffer : GG

    auto ta2_yz_xxzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 975);

    auto ta2_yz_xxzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 976);

    auto ta2_yz_xxzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 977);

    auto ta2_yz_xxzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 978);

    auto ta2_yz_xxzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 979);

    auto ta2_yz_xxzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 980);

    auto ta2_yz_xxzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 981);

    auto ta2_yz_xxzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 982);

    auto ta2_yz_xxzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 983);

    auto ta2_yz_xxzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 984);

    auto ta2_yz_xxzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 985);

    auto ta2_yz_xxzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 986);

    auto ta2_yz_xxzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 987);

    auto ta2_yz_xxzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 988);

    auto ta2_yz_xxzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 989);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xxz_xxxx_1,   \
                             ta1_y_xxz_xxxy_1,   \
                             ta1_y_xxz_xxyy_1,   \
                             ta1_y_xxz_xyyy_1,   \
                             ta2_yz_xx_xxxx_0,   \
                             ta2_yz_xx_xxxx_1,   \
                             ta2_yz_xx_xxxy_0,   \
                             ta2_yz_xx_xxxy_1,   \
                             ta2_yz_xx_xxyy_0,   \
                             ta2_yz_xx_xxyy_1,   \
                             ta2_yz_xx_xyyy_0,   \
                             ta2_yz_xx_xyyy_1,   \
                             ta2_yz_xxz_xxxx_0,  \
                             ta2_yz_xxz_xxxx_1,  \
                             ta2_yz_xxz_xxxy_0,  \
                             ta2_yz_xxz_xxxy_1,  \
                             ta2_yz_xxz_xxyy_0,  \
                             ta2_yz_xxz_xxyy_1,  \
                             ta2_yz_xxz_xyyy_0,  \
                             ta2_yz_xxz_xyyy_1,  \
                             ta2_yz_xxzz_xxxx_0, \
                             ta2_yz_xxzz_xxxy_0, \
                             ta2_yz_xxzz_xxxz_0, \
                             ta2_yz_xxzz_xxyy_0, \
                             ta2_yz_xxzz_xxyz_0, \
                             ta2_yz_xxzz_xxzz_0, \
                             ta2_yz_xxzz_xyyy_0, \
                             ta2_yz_xxzz_xyyz_0, \
                             ta2_yz_xxzz_xyzz_0, \
                             ta2_yz_xxzz_xzzz_0, \
                             ta2_yz_xxzz_yyyy_0, \
                             ta2_yz_xxzz_yyyz_0, \
                             ta2_yz_xxzz_yyzz_0, \
                             ta2_yz_xxzz_yzzz_0, \
                             ta2_yz_xxzz_zzzz_0, \
                             ta2_yz_xzz_xxxz_0,  \
                             ta2_yz_xzz_xxxz_1,  \
                             ta2_yz_xzz_xxyz_0,  \
                             ta2_yz_xzz_xxyz_1,  \
                             ta2_yz_xzz_xxz_0,   \
                             ta2_yz_xzz_xxz_1,   \
                             ta2_yz_xzz_xxzz_0,  \
                             ta2_yz_xzz_xxzz_1,  \
                             ta2_yz_xzz_xyyz_0,  \
                             ta2_yz_xzz_xyyz_1,  \
                             ta2_yz_xzz_xyz_0,   \
                             ta2_yz_xzz_xyz_1,   \
                             ta2_yz_xzz_xyzz_0,  \
                             ta2_yz_xzz_xyzz_1,  \
                             ta2_yz_xzz_xzz_0,   \
                             ta2_yz_xzz_xzz_1,   \
                             ta2_yz_xzz_xzzz_0,  \
                             ta2_yz_xzz_xzzz_1,  \
                             ta2_yz_xzz_yyyy_0,  \
                             ta2_yz_xzz_yyyy_1,  \
                             ta2_yz_xzz_yyyz_0,  \
                             ta2_yz_xzz_yyyz_1,  \
                             ta2_yz_xzz_yyz_0,   \
                             ta2_yz_xzz_yyz_1,   \
                             ta2_yz_xzz_yyzz_0,  \
                             ta2_yz_xzz_yyzz_1,  \
                             ta2_yz_xzz_yzz_0,   \
                             ta2_yz_xzz_yzz_1,   \
                             ta2_yz_xzz_yzzz_0,  \
                             ta2_yz_xzz_yzzz_1,  \
                             ta2_yz_xzz_zzz_0,   \
                             ta2_yz_xzz_zzz_1,   \
                             ta2_yz_xzz_zzzz_0,  \
                             ta2_yz_xzz_zzzz_1,  \
                             ta2_yz_zz_xxxz_0,   \
                             ta2_yz_zz_xxxz_1,   \
                             ta2_yz_zz_xxyz_0,   \
                             ta2_yz_zz_xxyz_1,   \
                             ta2_yz_zz_xxzz_0,   \
                             ta2_yz_zz_xxzz_1,   \
                             ta2_yz_zz_xyyz_0,   \
                             ta2_yz_zz_xyyz_1,   \
                             ta2_yz_zz_xyzz_0,   \
                             ta2_yz_zz_xyzz_1,   \
                             ta2_yz_zz_xzzz_0,   \
                             ta2_yz_zz_xzzz_1,   \
                             ta2_yz_zz_yyyy_0,   \
                             ta2_yz_zz_yyyy_1,   \
                             ta2_yz_zz_yyyz_0,   \
                             ta2_yz_zz_yyyz_1,   \
                             ta2_yz_zz_yyzz_0,   \
                             ta2_yz_zz_yyzz_1,   \
                             ta2_yz_zz_yzzz_0,   \
                             ta2_yz_zz_yzzz_1,   \
                             ta2_yz_zz_zzzz_0,   \
                             ta2_yz_zz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xxzz_xxxx_0[i] = ta2_yz_xx_xxxx_0[i] * fe_0 - ta2_yz_xx_xxxx_1[i] * fe_0 + ta1_y_xxz_xxxx_1[i] + ta2_yz_xxz_xxxx_0[i] * pa_z[i] -
                                ta2_yz_xxz_xxxx_1[i] * pc_z[i];

        ta2_yz_xxzz_xxxy_0[i] = ta2_yz_xx_xxxy_0[i] * fe_0 - ta2_yz_xx_xxxy_1[i] * fe_0 + ta1_y_xxz_xxxy_1[i] + ta2_yz_xxz_xxxy_0[i] * pa_z[i] -
                                ta2_yz_xxz_xxxy_1[i] * pc_z[i];

        ta2_yz_xxzz_xxxz_0[i] = ta2_yz_zz_xxxz_0[i] * fe_0 - ta2_yz_zz_xxxz_1[i] * fe_0 + 3.0 * ta2_yz_xzz_xxz_0[i] * fe_0 -
                                3.0 * ta2_yz_xzz_xxz_1[i] * fe_0 + ta2_yz_xzz_xxxz_0[i] * pa_x[i] - ta2_yz_xzz_xxxz_1[i] * pc_x[i];

        ta2_yz_xxzz_xxyy_0[i] = ta2_yz_xx_xxyy_0[i] * fe_0 - ta2_yz_xx_xxyy_1[i] * fe_0 + ta1_y_xxz_xxyy_1[i] + ta2_yz_xxz_xxyy_0[i] * pa_z[i] -
                                ta2_yz_xxz_xxyy_1[i] * pc_z[i];

        ta2_yz_xxzz_xxyz_0[i] = ta2_yz_zz_xxyz_0[i] * fe_0 - ta2_yz_zz_xxyz_1[i] * fe_0 + 2.0 * ta2_yz_xzz_xyz_0[i] * fe_0 -
                                2.0 * ta2_yz_xzz_xyz_1[i] * fe_0 + ta2_yz_xzz_xxyz_0[i] * pa_x[i] - ta2_yz_xzz_xxyz_1[i] * pc_x[i];

        ta2_yz_xxzz_xxzz_0[i] = ta2_yz_zz_xxzz_0[i] * fe_0 - ta2_yz_zz_xxzz_1[i] * fe_0 + 2.0 * ta2_yz_xzz_xzz_0[i] * fe_0 -
                                2.0 * ta2_yz_xzz_xzz_1[i] * fe_0 + ta2_yz_xzz_xxzz_0[i] * pa_x[i] - ta2_yz_xzz_xxzz_1[i] * pc_x[i];

        ta2_yz_xxzz_xyyy_0[i] = ta2_yz_xx_xyyy_0[i] * fe_0 - ta2_yz_xx_xyyy_1[i] * fe_0 + ta1_y_xxz_xyyy_1[i] + ta2_yz_xxz_xyyy_0[i] * pa_z[i] -
                                ta2_yz_xxz_xyyy_1[i] * pc_z[i];

        ta2_yz_xxzz_xyyz_0[i] = ta2_yz_zz_xyyz_0[i] * fe_0 - ta2_yz_zz_xyyz_1[i] * fe_0 + ta2_yz_xzz_yyz_0[i] * fe_0 - ta2_yz_xzz_yyz_1[i] * fe_0 +
                                ta2_yz_xzz_xyyz_0[i] * pa_x[i] - ta2_yz_xzz_xyyz_1[i] * pc_x[i];

        ta2_yz_xxzz_xyzz_0[i] = ta2_yz_zz_xyzz_0[i] * fe_0 - ta2_yz_zz_xyzz_1[i] * fe_0 + ta2_yz_xzz_yzz_0[i] * fe_0 - ta2_yz_xzz_yzz_1[i] * fe_0 +
                                ta2_yz_xzz_xyzz_0[i] * pa_x[i] - ta2_yz_xzz_xyzz_1[i] * pc_x[i];

        ta2_yz_xxzz_xzzz_0[i] = ta2_yz_zz_xzzz_0[i] * fe_0 - ta2_yz_zz_xzzz_1[i] * fe_0 + ta2_yz_xzz_zzz_0[i] * fe_0 - ta2_yz_xzz_zzz_1[i] * fe_0 +
                                ta2_yz_xzz_xzzz_0[i] * pa_x[i] - ta2_yz_xzz_xzzz_1[i] * pc_x[i];

        ta2_yz_xxzz_yyyy_0[i] =
            ta2_yz_zz_yyyy_0[i] * fe_0 - ta2_yz_zz_yyyy_1[i] * fe_0 + ta2_yz_xzz_yyyy_0[i] * pa_x[i] - ta2_yz_xzz_yyyy_1[i] * pc_x[i];

        ta2_yz_xxzz_yyyz_0[i] =
            ta2_yz_zz_yyyz_0[i] * fe_0 - ta2_yz_zz_yyyz_1[i] * fe_0 + ta2_yz_xzz_yyyz_0[i] * pa_x[i] - ta2_yz_xzz_yyyz_1[i] * pc_x[i];

        ta2_yz_xxzz_yyzz_0[i] =
            ta2_yz_zz_yyzz_0[i] * fe_0 - ta2_yz_zz_yyzz_1[i] * fe_0 + ta2_yz_xzz_yyzz_0[i] * pa_x[i] - ta2_yz_xzz_yyzz_1[i] * pc_x[i];

        ta2_yz_xxzz_yzzz_0[i] =
            ta2_yz_zz_yzzz_0[i] * fe_0 - ta2_yz_zz_yzzz_1[i] * fe_0 + ta2_yz_xzz_yzzz_0[i] * pa_x[i] - ta2_yz_xzz_yzzz_1[i] * pc_x[i];

        ta2_yz_xxzz_zzzz_0[i] =
            ta2_yz_zz_zzzz_0[i] * fe_0 - ta2_yz_zz_zzzz_1[i] * fe_0 + ta2_yz_xzz_zzzz_0[i] * pa_x[i] - ta2_yz_xzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 990-1005 components of targeted buffer : GG

    auto ta2_yz_xyyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 990);

    auto ta2_yz_xyyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 991);

    auto ta2_yz_xyyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 992);

    auto ta2_yz_xyyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 993);

    auto ta2_yz_xyyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 994);

    auto ta2_yz_xyyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 995);

    auto ta2_yz_xyyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 996);

    auto ta2_yz_xyyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 997);

    auto ta2_yz_xyyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 998);

    auto ta2_yz_xyyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 999);

    auto ta2_yz_xyyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1000);

    auto ta2_yz_xyyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1001);

    auto ta2_yz_xyyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1002);

    auto ta2_yz_xyyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1003);

    auto ta2_yz_xyyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1004);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta2_yz_xyyy_xxxx_0, \
                             ta2_yz_xyyy_xxxy_0, \
                             ta2_yz_xyyy_xxxz_0, \
                             ta2_yz_xyyy_xxyy_0, \
                             ta2_yz_xyyy_xxyz_0, \
                             ta2_yz_xyyy_xxzz_0, \
                             ta2_yz_xyyy_xyyy_0, \
                             ta2_yz_xyyy_xyyz_0, \
                             ta2_yz_xyyy_xyzz_0, \
                             ta2_yz_xyyy_xzzz_0, \
                             ta2_yz_xyyy_yyyy_0, \
                             ta2_yz_xyyy_yyyz_0, \
                             ta2_yz_xyyy_yyzz_0, \
                             ta2_yz_xyyy_yzzz_0, \
                             ta2_yz_xyyy_zzzz_0, \
                             ta2_yz_yyy_xxx_0,   \
                             ta2_yz_yyy_xxx_1,   \
                             ta2_yz_yyy_xxxx_0,  \
                             ta2_yz_yyy_xxxx_1,  \
                             ta2_yz_yyy_xxxy_0,  \
                             ta2_yz_yyy_xxxy_1,  \
                             ta2_yz_yyy_xxxz_0,  \
                             ta2_yz_yyy_xxxz_1,  \
                             ta2_yz_yyy_xxy_0,   \
                             ta2_yz_yyy_xxy_1,   \
                             ta2_yz_yyy_xxyy_0,  \
                             ta2_yz_yyy_xxyy_1,  \
                             ta2_yz_yyy_xxyz_0,  \
                             ta2_yz_yyy_xxyz_1,  \
                             ta2_yz_yyy_xxz_0,   \
                             ta2_yz_yyy_xxz_1,   \
                             ta2_yz_yyy_xxzz_0,  \
                             ta2_yz_yyy_xxzz_1,  \
                             ta2_yz_yyy_xyy_0,   \
                             ta2_yz_yyy_xyy_1,   \
                             ta2_yz_yyy_xyyy_0,  \
                             ta2_yz_yyy_xyyy_1,  \
                             ta2_yz_yyy_xyyz_0,  \
                             ta2_yz_yyy_xyyz_1,  \
                             ta2_yz_yyy_xyz_0,   \
                             ta2_yz_yyy_xyz_1,   \
                             ta2_yz_yyy_xyzz_0,  \
                             ta2_yz_yyy_xyzz_1,  \
                             ta2_yz_yyy_xzz_0,   \
                             ta2_yz_yyy_xzz_1,   \
                             ta2_yz_yyy_xzzz_0,  \
                             ta2_yz_yyy_xzzz_1,  \
                             ta2_yz_yyy_yyy_0,   \
                             ta2_yz_yyy_yyy_1,   \
                             ta2_yz_yyy_yyyy_0,  \
                             ta2_yz_yyy_yyyy_1,  \
                             ta2_yz_yyy_yyyz_0,  \
                             ta2_yz_yyy_yyyz_1,  \
                             ta2_yz_yyy_yyz_0,   \
                             ta2_yz_yyy_yyz_1,   \
                             ta2_yz_yyy_yyzz_0,  \
                             ta2_yz_yyy_yyzz_1,  \
                             ta2_yz_yyy_yzz_0,   \
                             ta2_yz_yyy_yzz_1,   \
                             ta2_yz_yyy_yzzz_0,  \
                             ta2_yz_yyy_yzzz_1,  \
                             ta2_yz_yyy_zzz_0,   \
                             ta2_yz_yyy_zzz_1,   \
                             ta2_yz_yyy_zzzz_0,  \
                             ta2_yz_yyy_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xyyy_xxxx_0[i] =
            4.0 * ta2_yz_yyy_xxx_0[i] * fe_0 - 4.0 * ta2_yz_yyy_xxx_1[i] * fe_0 + ta2_yz_yyy_xxxx_0[i] * pa_x[i] - ta2_yz_yyy_xxxx_1[i] * pc_x[i];

        ta2_yz_xyyy_xxxy_0[i] =
            3.0 * ta2_yz_yyy_xxy_0[i] * fe_0 - 3.0 * ta2_yz_yyy_xxy_1[i] * fe_0 + ta2_yz_yyy_xxxy_0[i] * pa_x[i] - ta2_yz_yyy_xxxy_1[i] * pc_x[i];

        ta2_yz_xyyy_xxxz_0[i] =
            3.0 * ta2_yz_yyy_xxz_0[i] * fe_0 - 3.0 * ta2_yz_yyy_xxz_1[i] * fe_0 + ta2_yz_yyy_xxxz_0[i] * pa_x[i] - ta2_yz_yyy_xxxz_1[i] * pc_x[i];

        ta2_yz_xyyy_xxyy_0[i] =
            2.0 * ta2_yz_yyy_xyy_0[i] * fe_0 - 2.0 * ta2_yz_yyy_xyy_1[i] * fe_0 + ta2_yz_yyy_xxyy_0[i] * pa_x[i] - ta2_yz_yyy_xxyy_1[i] * pc_x[i];

        ta2_yz_xyyy_xxyz_0[i] =
            2.0 * ta2_yz_yyy_xyz_0[i] * fe_0 - 2.0 * ta2_yz_yyy_xyz_1[i] * fe_0 + ta2_yz_yyy_xxyz_0[i] * pa_x[i] - ta2_yz_yyy_xxyz_1[i] * pc_x[i];

        ta2_yz_xyyy_xxzz_0[i] =
            2.0 * ta2_yz_yyy_xzz_0[i] * fe_0 - 2.0 * ta2_yz_yyy_xzz_1[i] * fe_0 + ta2_yz_yyy_xxzz_0[i] * pa_x[i] - ta2_yz_yyy_xxzz_1[i] * pc_x[i];

        ta2_yz_xyyy_xyyy_0[i] =
            ta2_yz_yyy_yyy_0[i] * fe_0 - ta2_yz_yyy_yyy_1[i] * fe_0 + ta2_yz_yyy_xyyy_0[i] * pa_x[i] - ta2_yz_yyy_xyyy_1[i] * pc_x[i];

        ta2_yz_xyyy_xyyz_0[i] =
            ta2_yz_yyy_yyz_0[i] * fe_0 - ta2_yz_yyy_yyz_1[i] * fe_0 + ta2_yz_yyy_xyyz_0[i] * pa_x[i] - ta2_yz_yyy_xyyz_1[i] * pc_x[i];

        ta2_yz_xyyy_xyzz_0[i] =
            ta2_yz_yyy_yzz_0[i] * fe_0 - ta2_yz_yyy_yzz_1[i] * fe_0 + ta2_yz_yyy_xyzz_0[i] * pa_x[i] - ta2_yz_yyy_xyzz_1[i] * pc_x[i];

        ta2_yz_xyyy_xzzz_0[i] =
            ta2_yz_yyy_zzz_0[i] * fe_0 - ta2_yz_yyy_zzz_1[i] * fe_0 + ta2_yz_yyy_xzzz_0[i] * pa_x[i] - ta2_yz_yyy_xzzz_1[i] * pc_x[i];

        ta2_yz_xyyy_yyyy_0[i] = ta2_yz_yyy_yyyy_0[i] * pa_x[i] - ta2_yz_yyy_yyyy_1[i] * pc_x[i];

        ta2_yz_xyyy_yyyz_0[i] = ta2_yz_yyy_yyyz_0[i] * pa_x[i] - ta2_yz_yyy_yyyz_1[i] * pc_x[i];

        ta2_yz_xyyy_yyzz_0[i] = ta2_yz_yyy_yyzz_0[i] * pa_x[i] - ta2_yz_yyy_yyzz_1[i] * pc_x[i];

        ta2_yz_xyyy_yzzz_0[i] = ta2_yz_yyy_yzzz_0[i] * pa_x[i] - ta2_yz_yyy_yzzz_1[i] * pc_x[i];

        ta2_yz_xyyy_zzzz_0[i] = ta2_yz_yyy_zzzz_0[i] * pa_x[i] - ta2_yz_yyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 1005-1020 components of targeted buffer : GG

    auto ta2_yz_xyyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1005);

    auto ta2_yz_xyyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1006);

    auto ta2_yz_xyyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1007);

    auto ta2_yz_xyyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1008);

    auto ta2_yz_xyyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1009);

    auto ta2_yz_xyyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1010);

    auto ta2_yz_xyyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1011);

    auto ta2_yz_xyyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1012);

    auto ta2_yz_xyyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1013);

    auto ta2_yz_xyyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1014);

    auto ta2_yz_xyyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1015);

    auto ta2_yz_xyyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1016);

    auto ta2_yz_xyyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1017);

    auto ta2_yz_xyyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1018);

    auto ta2_yz_xyyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1019);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xyy_xxxx_1,   \
                             ta1_y_xyy_xxxy_1,   \
                             ta1_y_xyy_xxyy_1,   \
                             ta1_y_xyy_xyyy_1,   \
                             ta2_yz_xyy_xxxx_0,  \
                             ta2_yz_xyy_xxxx_1,  \
                             ta2_yz_xyy_xxxy_0,  \
                             ta2_yz_xyy_xxxy_1,  \
                             ta2_yz_xyy_xxyy_0,  \
                             ta2_yz_xyy_xxyy_1,  \
                             ta2_yz_xyy_xyyy_0,  \
                             ta2_yz_xyy_xyyy_1,  \
                             ta2_yz_xyyz_xxxx_0, \
                             ta2_yz_xyyz_xxxy_0, \
                             ta2_yz_xyyz_xxxz_0, \
                             ta2_yz_xyyz_xxyy_0, \
                             ta2_yz_xyyz_xxyz_0, \
                             ta2_yz_xyyz_xxzz_0, \
                             ta2_yz_xyyz_xyyy_0, \
                             ta2_yz_xyyz_xyyz_0, \
                             ta2_yz_xyyz_xyzz_0, \
                             ta2_yz_xyyz_xzzz_0, \
                             ta2_yz_xyyz_yyyy_0, \
                             ta2_yz_xyyz_yyyz_0, \
                             ta2_yz_xyyz_yyzz_0, \
                             ta2_yz_xyyz_yzzz_0, \
                             ta2_yz_xyyz_zzzz_0, \
                             ta2_yz_yyz_xxxz_0,  \
                             ta2_yz_yyz_xxxz_1,  \
                             ta2_yz_yyz_xxyz_0,  \
                             ta2_yz_yyz_xxyz_1,  \
                             ta2_yz_yyz_xxz_0,   \
                             ta2_yz_yyz_xxz_1,   \
                             ta2_yz_yyz_xxzz_0,  \
                             ta2_yz_yyz_xxzz_1,  \
                             ta2_yz_yyz_xyyz_0,  \
                             ta2_yz_yyz_xyyz_1,  \
                             ta2_yz_yyz_xyz_0,   \
                             ta2_yz_yyz_xyz_1,   \
                             ta2_yz_yyz_xyzz_0,  \
                             ta2_yz_yyz_xyzz_1,  \
                             ta2_yz_yyz_xzz_0,   \
                             ta2_yz_yyz_xzz_1,   \
                             ta2_yz_yyz_xzzz_0,  \
                             ta2_yz_yyz_xzzz_1,  \
                             ta2_yz_yyz_yyyy_0,  \
                             ta2_yz_yyz_yyyy_1,  \
                             ta2_yz_yyz_yyyz_0,  \
                             ta2_yz_yyz_yyyz_1,  \
                             ta2_yz_yyz_yyz_0,   \
                             ta2_yz_yyz_yyz_1,   \
                             ta2_yz_yyz_yyzz_0,  \
                             ta2_yz_yyz_yyzz_1,  \
                             ta2_yz_yyz_yzz_0,   \
                             ta2_yz_yyz_yzz_1,   \
                             ta2_yz_yyz_yzzz_0,  \
                             ta2_yz_yyz_yzzz_1,  \
                             ta2_yz_yyz_zzz_0,   \
                             ta2_yz_yyz_zzz_1,   \
                             ta2_yz_yyz_zzzz_0,  \
                             ta2_yz_yyz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xyyz_xxxx_0[i] = ta1_y_xyy_xxxx_1[i] + ta2_yz_xyy_xxxx_0[i] * pa_z[i] - ta2_yz_xyy_xxxx_1[i] * pc_z[i];

        ta2_yz_xyyz_xxxy_0[i] = ta1_y_xyy_xxxy_1[i] + ta2_yz_xyy_xxxy_0[i] * pa_z[i] - ta2_yz_xyy_xxxy_1[i] * pc_z[i];

        ta2_yz_xyyz_xxxz_0[i] =
            3.0 * ta2_yz_yyz_xxz_0[i] * fe_0 - 3.0 * ta2_yz_yyz_xxz_1[i] * fe_0 + ta2_yz_yyz_xxxz_0[i] * pa_x[i] - ta2_yz_yyz_xxxz_1[i] * pc_x[i];

        ta2_yz_xyyz_xxyy_0[i] = ta1_y_xyy_xxyy_1[i] + ta2_yz_xyy_xxyy_0[i] * pa_z[i] - ta2_yz_xyy_xxyy_1[i] * pc_z[i];

        ta2_yz_xyyz_xxyz_0[i] =
            2.0 * ta2_yz_yyz_xyz_0[i] * fe_0 - 2.0 * ta2_yz_yyz_xyz_1[i] * fe_0 + ta2_yz_yyz_xxyz_0[i] * pa_x[i] - ta2_yz_yyz_xxyz_1[i] * pc_x[i];

        ta2_yz_xyyz_xxzz_0[i] =
            2.0 * ta2_yz_yyz_xzz_0[i] * fe_0 - 2.0 * ta2_yz_yyz_xzz_1[i] * fe_0 + ta2_yz_yyz_xxzz_0[i] * pa_x[i] - ta2_yz_yyz_xxzz_1[i] * pc_x[i];

        ta2_yz_xyyz_xyyy_0[i] = ta1_y_xyy_xyyy_1[i] + ta2_yz_xyy_xyyy_0[i] * pa_z[i] - ta2_yz_xyy_xyyy_1[i] * pc_z[i];

        ta2_yz_xyyz_xyyz_0[i] =
            ta2_yz_yyz_yyz_0[i] * fe_0 - ta2_yz_yyz_yyz_1[i] * fe_0 + ta2_yz_yyz_xyyz_0[i] * pa_x[i] - ta2_yz_yyz_xyyz_1[i] * pc_x[i];

        ta2_yz_xyyz_xyzz_0[i] =
            ta2_yz_yyz_yzz_0[i] * fe_0 - ta2_yz_yyz_yzz_1[i] * fe_0 + ta2_yz_yyz_xyzz_0[i] * pa_x[i] - ta2_yz_yyz_xyzz_1[i] * pc_x[i];

        ta2_yz_xyyz_xzzz_0[i] =
            ta2_yz_yyz_zzz_0[i] * fe_0 - ta2_yz_yyz_zzz_1[i] * fe_0 + ta2_yz_yyz_xzzz_0[i] * pa_x[i] - ta2_yz_yyz_xzzz_1[i] * pc_x[i];

        ta2_yz_xyyz_yyyy_0[i] = ta2_yz_yyz_yyyy_0[i] * pa_x[i] - ta2_yz_yyz_yyyy_1[i] * pc_x[i];

        ta2_yz_xyyz_yyyz_0[i] = ta2_yz_yyz_yyyz_0[i] * pa_x[i] - ta2_yz_yyz_yyyz_1[i] * pc_x[i];

        ta2_yz_xyyz_yyzz_0[i] = ta2_yz_yyz_yyzz_0[i] * pa_x[i] - ta2_yz_yyz_yyzz_1[i] * pc_x[i];

        ta2_yz_xyyz_yzzz_0[i] = ta2_yz_yyz_yzzz_0[i] * pa_x[i] - ta2_yz_yyz_yzzz_1[i] * pc_x[i];

        ta2_yz_xyyz_zzzz_0[i] = ta2_yz_yyz_zzzz_0[i] * pa_x[i] - ta2_yz_yyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 1020-1035 components of targeted buffer : GG

    auto ta2_yz_xyzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1020);

    auto ta2_yz_xyzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1021);

    auto ta2_yz_xyzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1022);

    auto ta2_yz_xyzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1023);

    auto ta2_yz_xyzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1024);

    auto ta2_yz_xyzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1025);

    auto ta2_yz_xyzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1026);

    auto ta2_yz_xyzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1027);

    auto ta2_yz_xyzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1028);

    auto ta2_yz_xyzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1029);

    auto ta2_yz_xyzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1030);

    auto ta2_yz_xyzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1031);

    auto ta2_yz_xyzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1032);

    auto ta2_yz_xyzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1033);

    auto ta2_yz_xyzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1034);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xzz_xxxx_1,   \
                             ta1_z_xzz_xxxz_1,   \
                             ta1_z_xzz_xxzz_1,   \
                             ta1_z_xzz_xzzz_1,   \
                             ta2_yz_xyzz_xxxx_0, \
                             ta2_yz_xyzz_xxxy_0, \
                             ta2_yz_xyzz_xxxz_0, \
                             ta2_yz_xyzz_xxyy_0, \
                             ta2_yz_xyzz_xxyz_0, \
                             ta2_yz_xyzz_xxzz_0, \
                             ta2_yz_xyzz_xyyy_0, \
                             ta2_yz_xyzz_xyyz_0, \
                             ta2_yz_xyzz_xyzz_0, \
                             ta2_yz_xyzz_xzzz_0, \
                             ta2_yz_xyzz_yyyy_0, \
                             ta2_yz_xyzz_yyyz_0, \
                             ta2_yz_xyzz_yyzz_0, \
                             ta2_yz_xyzz_yzzz_0, \
                             ta2_yz_xyzz_zzzz_0, \
                             ta2_yz_xzz_xxxx_0,  \
                             ta2_yz_xzz_xxxx_1,  \
                             ta2_yz_xzz_xxxz_0,  \
                             ta2_yz_xzz_xxxz_1,  \
                             ta2_yz_xzz_xxzz_0,  \
                             ta2_yz_xzz_xxzz_1,  \
                             ta2_yz_xzz_xzzz_0,  \
                             ta2_yz_xzz_xzzz_1,  \
                             ta2_yz_yzz_xxxy_0,  \
                             ta2_yz_yzz_xxxy_1,  \
                             ta2_yz_yzz_xxy_0,   \
                             ta2_yz_yzz_xxy_1,   \
                             ta2_yz_yzz_xxyy_0,  \
                             ta2_yz_yzz_xxyy_1,  \
                             ta2_yz_yzz_xxyz_0,  \
                             ta2_yz_yzz_xxyz_1,  \
                             ta2_yz_yzz_xyy_0,   \
                             ta2_yz_yzz_xyy_1,   \
                             ta2_yz_yzz_xyyy_0,  \
                             ta2_yz_yzz_xyyy_1,  \
                             ta2_yz_yzz_xyyz_0,  \
                             ta2_yz_yzz_xyyz_1,  \
                             ta2_yz_yzz_xyz_0,   \
                             ta2_yz_yzz_xyz_1,   \
                             ta2_yz_yzz_xyzz_0,  \
                             ta2_yz_yzz_xyzz_1,  \
                             ta2_yz_yzz_yyy_0,   \
                             ta2_yz_yzz_yyy_1,   \
                             ta2_yz_yzz_yyyy_0,  \
                             ta2_yz_yzz_yyyy_1,  \
                             ta2_yz_yzz_yyyz_0,  \
                             ta2_yz_yzz_yyyz_1,  \
                             ta2_yz_yzz_yyz_0,   \
                             ta2_yz_yzz_yyz_1,   \
                             ta2_yz_yzz_yyzz_0,  \
                             ta2_yz_yzz_yyzz_1,  \
                             ta2_yz_yzz_yzz_0,   \
                             ta2_yz_yzz_yzz_1,   \
                             ta2_yz_yzz_yzzz_0,  \
                             ta2_yz_yzz_yzzz_1,  \
                             ta2_yz_yzz_zzzz_0,  \
                             ta2_yz_yzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xyzz_xxxx_0[i] = ta1_z_xzz_xxxx_1[i] + ta2_yz_xzz_xxxx_0[i] * pa_y[i] - ta2_yz_xzz_xxxx_1[i] * pc_y[i];

        ta2_yz_xyzz_xxxy_0[i] =
            3.0 * ta2_yz_yzz_xxy_0[i] * fe_0 - 3.0 * ta2_yz_yzz_xxy_1[i] * fe_0 + ta2_yz_yzz_xxxy_0[i] * pa_x[i] - ta2_yz_yzz_xxxy_1[i] * pc_x[i];

        ta2_yz_xyzz_xxxz_0[i] = ta1_z_xzz_xxxz_1[i] + ta2_yz_xzz_xxxz_0[i] * pa_y[i] - ta2_yz_xzz_xxxz_1[i] * pc_y[i];

        ta2_yz_xyzz_xxyy_0[i] =
            2.0 * ta2_yz_yzz_xyy_0[i] * fe_0 - 2.0 * ta2_yz_yzz_xyy_1[i] * fe_0 + ta2_yz_yzz_xxyy_0[i] * pa_x[i] - ta2_yz_yzz_xxyy_1[i] * pc_x[i];

        ta2_yz_xyzz_xxyz_0[i] =
            2.0 * ta2_yz_yzz_xyz_0[i] * fe_0 - 2.0 * ta2_yz_yzz_xyz_1[i] * fe_0 + ta2_yz_yzz_xxyz_0[i] * pa_x[i] - ta2_yz_yzz_xxyz_1[i] * pc_x[i];

        ta2_yz_xyzz_xxzz_0[i] = ta1_z_xzz_xxzz_1[i] + ta2_yz_xzz_xxzz_0[i] * pa_y[i] - ta2_yz_xzz_xxzz_1[i] * pc_y[i];

        ta2_yz_xyzz_xyyy_0[i] =
            ta2_yz_yzz_yyy_0[i] * fe_0 - ta2_yz_yzz_yyy_1[i] * fe_0 + ta2_yz_yzz_xyyy_0[i] * pa_x[i] - ta2_yz_yzz_xyyy_1[i] * pc_x[i];

        ta2_yz_xyzz_xyyz_0[i] =
            ta2_yz_yzz_yyz_0[i] * fe_0 - ta2_yz_yzz_yyz_1[i] * fe_0 + ta2_yz_yzz_xyyz_0[i] * pa_x[i] - ta2_yz_yzz_xyyz_1[i] * pc_x[i];

        ta2_yz_xyzz_xyzz_0[i] =
            ta2_yz_yzz_yzz_0[i] * fe_0 - ta2_yz_yzz_yzz_1[i] * fe_0 + ta2_yz_yzz_xyzz_0[i] * pa_x[i] - ta2_yz_yzz_xyzz_1[i] * pc_x[i];

        ta2_yz_xyzz_xzzz_0[i] = ta1_z_xzz_xzzz_1[i] + ta2_yz_xzz_xzzz_0[i] * pa_y[i] - ta2_yz_xzz_xzzz_1[i] * pc_y[i];

        ta2_yz_xyzz_yyyy_0[i] = ta2_yz_yzz_yyyy_0[i] * pa_x[i] - ta2_yz_yzz_yyyy_1[i] * pc_x[i];

        ta2_yz_xyzz_yyyz_0[i] = ta2_yz_yzz_yyyz_0[i] * pa_x[i] - ta2_yz_yzz_yyyz_1[i] * pc_x[i];

        ta2_yz_xyzz_yyzz_0[i] = ta2_yz_yzz_yyzz_0[i] * pa_x[i] - ta2_yz_yzz_yyzz_1[i] * pc_x[i];

        ta2_yz_xyzz_yzzz_0[i] = ta2_yz_yzz_yzzz_0[i] * pa_x[i] - ta2_yz_yzz_yzzz_1[i] * pc_x[i];

        ta2_yz_xyzz_zzzz_0[i] = ta2_yz_yzz_zzzz_0[i] * pa_x[i] - ta2_yz_yzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 1035-1050 components of targeted buffer : GG

    auto ta2_yz_xzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1035);

    auto ta2_yz_xzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1036);

    auto ta2_yz_xzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1037);

    auto ta2_yz_xzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1038);

    auto ta2_yz_xzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1039);

    auto ta2_yz_xzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1040);

    auto ta2_yz_xzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1041);

    auto ta2_yz_xzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1042);

    auto ta2_yz_xzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1043);

    auto ta2_yz_xzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1044);

    auto ta2_yz_xzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1045);

    auto ta2_yz_xzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1046);

    auto ta2_yz_xzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1047);

    auto ta2_yz_xzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1048);

    auto ta2_yz_xzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1049);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta2_yz_xzzz_xxxx_0, \
                             ta2_yz_xzzz_xxxy_0, \
                             ta2_yz_xzzz_xxxz_0, \
                             ta2_yz_xzzz_xxyy_0, \
                             ta2_yz_xzzz_xxyz_0, \
                             ta2_yz_xzzz_xxzz_0, \
                             ta2_yz_xzzz_xyyy_0, \
                             ta2_yz_xzzz_xyyz_0, \
                             ta2_yz_xzzz_xyzz_0, \
                             ta2_yz_xzzz_xzzz_0, \
                             ta2_yz_xzzz_yyyy_0, \
                             ta2_yz_xzzz_yyyz_0, \
                             ta2_yz_xzzz_yyzz_0, \
                             ta2_yz_xzzz_yzzz_0, \
                             ta2_yz_xzzz_zzzz_0, \
                             ta2_yz_zzz_xxx_0,   \
                             ta2_yz_zzz_xxx_1,   \
                             ta2_yz_zzz_xxxx_0,  \
                             ta2_yz_zzz_xxxx_1,  \
                             ta2_yz_zzz_xxxy_0,  \
                             ta2_yz_zzz_xxxy_1,  \
                             ta2_yz_zzz_xxxz_0,  \
                             ta2_yz_zzz_xxxz_1,  \
                             ta2_yz_zzz_xxy_0,   \
                             ta2_yz_zzz_xxy_1,   \
                             ta2_yz_zzz_xxyy_0,  \
                             ta2_yz_zzz_xxyy_1,  \
                             ta2_yz_zzz_xxyz_0,  \
                             ta2_yz_zzz_xxyz_1,  \
                             ta2_yz_zzz_xxz_0,   \
                             ta2_yz_zzz_xxz_1,   \
                             ta2_yz_zzz_xxzz_0,  \
                             ta2_yz_zzz_xxzz_1,  \
                             ta2_yz_zzz_xyy_0,   \
                             ta2_yz_zzz_xyy_1,   \
                             ta2_yz_zzz_xyyy_0,  \
                             ta2_yz_zzz_xyyy_1,  \
                             ta2_yz_zzz_xyyz_0,  \
                             ta2_yz_zzz_xyyz_1,  \
                             ta2_yz_zzz_xyz_0,   \
                             ta2_yz_zzz_xyz_1,   \
                             ta2_yz_zzz_xyzz_0,  \
                             ta2_yz_zzz_xyzz_1,  \
                             ta2_yz_zzz_xzz_0,   \
                             ta2_yz_zzz_xzz_1,   \
                             ta2_yz_zzz_xzzz_0,  \
                             ta2_yz_zzz_xzzz_1,  \
                             ta2_yz_zzz_yyy_0,   \
                             ta2_yz_zzz_yyy_1,   \
                             ta2_yz_zzz_yyyy_0,  \
                             ta2_yz_zzz_yyyy_1,  \
                             ta2_yz_zzz_yyyz_0,  \
                             ta2_yz_zzz_yyyz_1,  \
                             ta2_yz_zzz_yyz_0,   \
                             ta2_yz_zzz_yyz_1,   \
                             ta2_yz_zzz_yyzz_0,  \
                             ta2_yz_zzz_yyzz_1,  \
                             ta2_yz_zzz_yzz_0,   \
                             ta2_yz_zzz_yzz_1,   \
                             ta2_yz_zzz_yzzz_0,  \
                             ta2_yz_zzz_yzzz_1,  \
                             ta2_yz_zzz_zzz_0,   \
                             ta2_yz_zzz_zzz_1,   \
                             ta2_yz_zzz_zzzz_0,  \
                             ta2_yz_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_xzzz_xxxx_0[i] =
            4.0 * ta2_yz_zzz_xxx_0[i] * fe_0 - 4.0 * ta2_yz_zzz_xxx_1[i] * fe_0 + ta2_yz_zzz_xxxx_0[i] * pa_x[i] - ta2_yz_zzz_xxxx_1[i] * pc_x[i];

        ta2_yz_xzzz_xxxy_0[i] =
            3.0 * ta2_yz_zzz_xxy_0[i] * fe_0 - 3.0 * ta2_yz_zzz_xxy_1[i] * fe_0 + ta2_yz_zzz_xxxy_0[i] * pa_x[i] - ta2_yz_zzz_xxxy_1[i] * pc_x[i];

        ta2_yz_xzzz_xxxz_0[i] =
            3.0 * ta2_yz_zzz_xxz_0[i] * fe_0 - 3.0 * ta2_yz_zzz_xxz_1[i] * fe_0 + ta2_yz_zzz_xxxz_0[i] * pa_x[i] - ta2_yz_zzz_xxxz_1[i] * pc_x[i];

        ta2_yz_xzzz_xxyy_0[i] =
            2.0 * ta2_yz_zzz_xyy_0[i] * fe_0 - 2.0 * ta2_yz_zzz_xyy_1[i] * fe_0 + ta2_yz_zzz_xxyy_0[i] * pa_x[i] - ta2_yz_zzz_xxyy_1[i] * pc_x[i];

        ta2_yz_xzzz_xxyz_0[i] =
            2.0 * ta2_yz_zzz_xyz_0[i] * fe_0 - 2.0 * ta2_yz_zzz_xyz_1[i] * fe_0 + ta2_yz_zzz_xxyz_0[i] * pa_x[i] - ta2_yz_zzz_xxyz_1[i] * pc_x[i];

        ta2_yz_xzzz_xxzz_0[i] =
            2.0 * ta2_yz_zzz_xzz_0[i] * fe_0 - 2.0 * ta2_yz_zzz_xzz_1[i] * fe_0 + ta2_yz_zzz_xxzz_0[i] * pa_x[i] - ta2_yz_zzz_xxzz_1[i] * pc_x[i];

        ta2_yz_xzzz_xyyy_0[i] =
            ta2_yz_zzz_yyy_0[i] * fe_0 - ta2_yz_zzz_yyy_1[i] * fe_0 + ta2_yz_zzz_xyyy_0[i] * pa_x[i] - ta2_yz_zzz_xyyy_1[i] * pc_x[i];

        ta2_yz_xzzz_xyyz_0[i] =
            ta2_yz_zzz_yyz_0[i] * fe_0 - ta2_yz_zzz_yyz_1[i] * fe_0 + ta2_yz_zzz_xyyz_0[i] * pa_x[i] - ta2_yz_zzz_xyyz_1[i] * pc_x[i];

        ta2_yz_xzzz_xyzz_0[i] =
            ta2_yz_zzz_yzz_0[i] * fe_0 - ta2_yz_zzz_yzz_1[i] * fe_0 + ta2_yz_zzz_xyzz_0[i] * pa_x[i] - ta2_yz_zzz_xyzz_1[i] * pc_x[i];

        ta2_yz_xzzz_xzzz_0[i] =
            ta2_yz_zzz_zzz_0[i] * fe_0 - ta2_yz_zzz_zzz_1[i] * fe_0 + ta2_yz_zzz_xzzz_0[i] * pa_x[i] - ta2_yz_zzz_xzzz_1[i] * pc_x[i];

        ta2_yz_xzzz_yyyy_0[i] = ta2_yz_zzz_yyyy_0[i] * pa_x[i] - ta2_yz_zzz_yyyy_1[i] * pc_x[i];

        ta2_yz_xzzz_yyyz_0[i] = ta2_yz_zzz_yyyz_0[i] * pa_x[i] - ta2_yz_zzz_yyyz_1[i] * pc_x[i];

        ta2_yz_xzzz_yyzz_0[i] = ta2_yz_zzz_yyzz_0[i] * pa_x[i] - ta2_yz_zzz_yyzz_1[i] * pc_x[i];

        ta2_yz_xzzz_yzzz_0[i] = ta2_yz_zzz_yzzz_0[i] * pa_x[i] - ta2_yz_zzz_yzzz_1[i] * pc_x[i];

        ta2_yz_xzzz_zzzz_0[i] = ta2_yz_zzz_zzzz_0[i] * pa_x[i] - ta2_yz_zzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 1050-1065 components of targeted buffer : GG

    auto ta2_yz_yyyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1050);

    auto ta2_yz_yyyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1051);

    auto ta2_yz_yyyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1052);

    auto ta2_yz_yyyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1053);

    auto ta2_yz_yyyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1054);

    auto ta2_yz_yyyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1055);

    auto ta2_yz_yyyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1056);

    auto ta2_yz_yyyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1057);

    auto ta2_yz_yyyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1058);

    auto ta2_yz_yyyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1059);

    auto ta2_yz_yyyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1060);

    auto ta2_yz_yyyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1061);

    auto ta2_yz_yyyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1062);

    auto ta2_yz_yyyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1063);

    auto ta2_yz_yyyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1064);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_z_yyy_xxxx_1,   \
                             ta1_z_yyy_xxxy_1,   \
                             ta1_z_yyy_xxxz_1,   \
                             ta1_z_yyy_xxyy_1,   \
                             ta1_z_yyy_xxyz_1,   \
                             ta1_z_yyy_xxzz_1,   \
                             ta1_z_yyy_xyyy_1,   \
                             ta1_z_yyy_xyyz_1,   \
                             ta1_z_yyy_xyzz_1,   \
                             ta1_z_yyy_xzzz_1,   \
                             ta1_z_yyy_yyyy_1,   \
                             ta1_z_yyy_yyyz_1,   \
                             ta1_z_yyy_yyzz_1,   \
                             ta1_z_yyy_yzzz_1,   \
                             ta1_z_yyy_zzzz_1,   \
                             ta2_yz_yy_xxxx_0,   \
                             ta2_yz_yy_xxxx_1,   \
                             ta2_yz_yy_xxxy_0,   \
                             ta2_yz_yy_xxxy_1,   \
                             ta2_yz_yy_xxxz_0,   \
                             ta2_yz_yy_xxxz_1,   \
                             ta2_yz_yy_xxyy_0,   \
                             ta2_yz_yy_xxyy_1,   \
                             ta2_yz_yy_xxyz_0,   \
                             ta2_yz_yy_xxyz_1,   \
                             ta2_yz_yy_xxzz_0,   \
                             ta2_yz_yy_xxzz_1,   \
                             ta2_yz_yy_xyyy_0,   \
                             ta2_yz_yy_xyyy_1,   \
                             ta2_yz_yy_xyyz_0,   \
                             ta2_yz_yy_xyyz_1,   \
                             ta2_yz_yy_xyzz_0,   \
                             ta2_yz_yy_xyzz_1,   \
                             ta2_yz_yy_xzzz_0,   \
                             ta2_yz_yy_xzzz_1,   \
                             ta2_yz_yy_yyyy_0,   \
                             ta2_yz_yy_yyyy_1,   \
                             ta2_yz_yy_yyyz_0,   \
                             ta2_yz_yy_yyyz_1,   \
                             ta2_yz_yy_yyzz_0,   \
                             ta2_yz_yy_yyzz_1,   \
                             ta2_yz_yy_yzzz_0,   \
                             ta2_yz_yy_yzzz_1,   \
                             ta2_yz_yy_zzzz_0,   \
                             ta2_yz_yy_zzzz_1,   \
                             ta2_yz_yyy_xxx_0,   \
                             ta2_yz_yyy_xxx_1,   \
                             ta2_yz_yyy_xxxx_0,  \
                             ta2_yz_yyy_xxxx_1,  \
                             ta2_yz_yyy_xxxy_0,  \
                             ta2_yz_yyy_xxxy_1,  \
                             ta2_yz_yyy_xxxz_0,  \
                             ta2_yz_yyy_xxxz_1,  \
                             ta2_yz_yyy_xxy_0,   \
                             ta2_yz_yyy_xxy_1,   \
                             ta2_yz_yyy_xxyy_0,  \
                             ta2_yz_yyy_xxyy_1,  \
                             ta2_yz_yyy_xxyz_0,  \
                             ta2_yz_yyy_xxyz_1,  \
                             ta2_yz_yyy_xxz_0,   \
                             ta2_yz_yyy_xxz_1,   \
                             ta2_yz_yyy_xxzz_0,  \
                             ta2_yz_yyy_xxzz_1,  \
                             ta2_yz_yyy_xyy_0,   \
                             ta2_yz_yyy_xyy_1,   \
                             ta2_yz_yyy_xyyy_0,  \
                             ta2_yz_yyy_xyyy_1,  \
                             ta2_yz_yyy_xyyz_0,  \
                             ta2_yz_yyy_xyyz_1,  \
                             ta2_yz_yyy_xyz_0,   \
                             ta2_yz_yyy_xyz_1,   \
                             ta2_yz_yyy_xyzz_0,  \
                             ta2_yz_yyy_xyzz_1,  \
                             ta2_yz_yyy_xzz_0,   \
                             ta2_yz_yyy_xzz_1,   \
                             ta2_yz_yyy_xzzz_0,  \
                             ta2_yz_yyy_xzzz_1,  \
                             ta2_yz_yyy_yyy_0,   \
                             ta2_yz_yyy_yyy_1,   \
                             ta2_yz_yyy_yyyy_0,  \
                             ta2_yz_yyy_yyyy_1,  \
                             ta2_yz_yyy_yyyz_0,  \
                             ta2_yz_yyy_yyyz_1,  \
                             ta2_yz_yyy_yyz_0,   \
                             ta2_yz_yyy_yyz_1,   \
                             ta2_yz_yyy_yyzz_0,  \
                             ta2_yz_yyy_yyzz_1,  \
                             ta2_yz_yyy_yzz_0,   \
                             ta2_yz_yyy_yzz_1,   \
                             ta2_yz_yyy_yzzz_0,  \
                             ta2_yz_yyy_yzzz_1,  \
                             ta2_yz_yyy_zzz_0,   \
                             ta2_yz_yyy_zzz_1,   \
                             ta2_yz_yyy_zzzz_0,  \
                             ta2_yz_yyy_zzzz_1,  \
                             ta2_yz_yyyy_xxxx_0, \
                             ta2_yz_yyyy_xxxy_0, \
                             ta2_yz_yyyy_xxxz_0, \
                             ta2_yz_yyyy_xxyy_0, \
                             ta2_yz_yyyy_xxyz_0, \
                             ta2_yz_yyyy_xxzz_0, \
                             ta2_yz_yyyy_xyyy_0, \
                             ta2_yz_yyyy_xyyz_0, \
                             ta2_yz_yyyy_xyzz_0, \
                             ta2_yz_yyyy_xzzz_0, \
                             ta2_yz_yyyy_yyyy_0, \
                             ta2_yz_yyyy_yyyz_0, \
                             ta2_yz_yyyy_yyzz_0, \
                             ta2_yz_yyyy_yzzz_0, \
                             ta2_yz_yyyy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyyy_xxxx_0[i] = 3.0 * ta2_yz_yy_xxxx_0[i] * fe_0 - 3.0 * ta2_yz_yy_xxxx_1[i] * fe_0 + ta1_z_yyy_xxxx_1[i] +
                                ta2_yz_yyy_xxxx_0[i] * pa_y[i] - ta2_yz_yyy_xxxx_1[i] * pc_y[i];

        ta2_yz_yyyy_xxxy_0[i] = 3.0 * ta2_yz_yy_xxxy_0[i] * fe_0 - 3.0 * ta2_yz_yy_xxxy_1[i] * fe_0 + ta2_yz_yyy_xxx_0[i] * fe_0 -
                                ta2_yz_yyy_xxx_1[i] * fe_0 + ta1_z_yyy_xxxy_1[i] + ta2_yz_yyy_xxxy_0[i] * pa_y[i] - ta2_yz_yyy_xxxy_1[i] * pc_y[i];

        ta2_yz_yyyy_xxxz_0[i] = 3.0 * ta2_yz_yy_xxxz_0[i] * fe_0 - 3.0 * ta2_yz_yy_xxxz_1[i] * fe_0 + ta1_z_yyy_xxxz_1[i] +
                                ta2_yz_yyy_xxxz_0[i] * pa_y[i] - ta2_yz_yyy_xxxz_1[i] * pc_y[i];

        ta2_yz_yyyy_xxyy_0[i] = 3.0 * ta2_yz_yy_xxyy_0[i] * fe_0 - 3.0 * ta2_yz_yy_xxyy_1[i] * fe_0 + 2.0 * ta2_yz_yyy_xxy_0[i] * fe_0 -
                                2.0 * ta2_yz_yyy_xxy_1[i] * fe_0 + ta1_z_yyy_xxyy_1[i] + ta2_yz_yyy_xxyy_0[i] * pa_y[i] -
                                ta2_yz_yyy_xxyy_1[i] * pc_y[i];

        ta2_yz_yyyy_xxyz_0[i] = 3.0 * ta2_yz_yy_xxyz_0[i] * fe_0 - 3.0 * ta2_yz_yy_xxyz_1[i] * fe_0 + ta2_yz_yyy_xxz_0[i] * fe_0 -
                                ta2_yz_yyy_xxz_1[i] * fe_0 + ta1_z_yyy_xxyz_1[i] + ta2_yz_yyy_xxyz_0[i] * pa_y[i] - ta2_yz_yyy_xxyz_1[i] * pc_y[i];

        ta2_yz_yyyy_xxzz_0[i] = 3.0 * ta2_yz_yy_xxzz_0[i] * fe_0 - 3.0 * ta2_yz_yy_xxzz_1[i] * fe_0 + ta1_z_yyy_xxzz_1[i] +
                                ta2_yz_yyy_xxzz_0[i] * pa_y[i] - ta2_yz_yyy_xxzz_1[i] * pc_y[i];

        ta2_yz_yyyy_xyyy_0[i] = 3.0 * ta2_yz_yy_xyyy_0[i] * fe_0 - 3.0 * ta2_yz_yy_xyyy_1[i] * fe_0 + 3.0 * ta2_yz_yyy_xyy_0[i] * fe_0 -
                                3.0 * ta2_yz_yyy_xyy_1[i] * fe_0 + ta1_z_yyy_xyyy_1[i] + ta2_yz_yyy_xyyy_0[i] * pa_y[i] -
                                ta2_yz_yyy_xyyy_1[i] * pc_y[i];

        ta2_yz_yyyy_xyyz_0[i] = 3.0 * ta2_yz_yy_xyyz_0[i] * fe_0 - 3.0 * ta2_yz_yy_xyyz_1[i] * fe_0 + 2.0 * ta2_yz_yyy_xyz_0[i] * fe_0 -
                                2.0 * ta2_yz_yyy_xyz_1[i] * fe_0 + ta1_z_yyy_xyyz_1[i] + ta2_yz_yyy_xyyz_0[i] * pa_y[i] -
                                ta2_yz_yyy_xyyz_1[i] * pc_y[i];

        ta2_yz_yyyy_xyzz_0[i] = 3.0 * ta2_yz_yy_xyzz_0[i] * fe_0 - 3.0 * ta2_yz_yy_xyzz_1[i] * fe_0 + ta2_yz_yyy_xzz_0[i] * fe_0 -
                                ta2_yz_yyy_xzz_1[i] * fe_0 + ta1_z_yyy_xyzz_1[i] + ta2_yz_yyy_xyzz_0[i] * pa_y[i] - ta2_yz_yyy_xyzz_1[i] * pc_y[i];

        ta2_yz_yyyy_xzzz_0[i] = 3.0 * ta2_yz_yy_xzzz_0[i] * fe_0 - 3.0 * ta2_yz_yy_xzzz_1[i] * fe_0 + ta1_z_yyy_xzzz_1[i] +
                                ta2_yz_yyy_xzzz_0[i] * pa_y[i] - ta2_yz_yyy_xzzz_1[i] * pc_y[i];

        ta2_yz_yyyy_yyyy_0[i] = 3.0 * ta2_yz_yy_yyyy_0[i] * fe_0 - 3.0 * ta2_yz_yy_yyyy_1[i] * fe_0 + 4.0 * ta2_yz_yyy_yyy_0[i] * fe_0 -
                                4.0 * ta2_yz_yyy_yyy_1[i] * fe_0 + ta1_z_yyy_yyyy_1[i] + ta2_yz_yyy_yyyy_0[i] * pa_y[i] -
                                ta2_yz_yyy_yyyy_1[i] * pc_y[i];

        ta2_yz_yyyy_yyyz_0[i] = 3.0 * ta2_yz_yy_yyyz_0[i] * fe_0 - 3.0 * ta2_yz_yy_yyyz_1[i] * fe_0 + 3.0 * ta2_yz_yyy_yyz_0[i] * fe_0 -
                                3.0 * ta2_yz_yyy_yyz_1[i] * fe_0 + ta1_z_yyy_yyyz_1[i] + ta2_yz_yyy_yyyz_0[i] * pa_y[i] -
                                ta2_yz_yyy_yyyz_1[i] * pc_y[i];

        ta2_yz_yyyy_yyzz_0[i] = 3.0 * ta2_yz_yy_yyzz_0[i] * fe_0 - 3.0 * ta2_yz_yy_yyzz_1[i] * fe_0 + 2.0 * ta2_yz_yyy_yzz_0[i] * fe_0 -
                                2.0 * ta2_yz_yyy_yzz_1[i] * fe_0 + ta1_z_yyy_yyzz_1[i] + ta2_yz_yyy_yyzz_0[i] * pa_y[i] -
                                ta2_yz_yyy_yyzz_1[i] * pc_y[i];

        ta2_yz_yyyy_yzzz_0[i] = 3.0 * ta2_yz_yy_yzzz_0[i] * fe_0 - 3.0 * ta2_yz_yy_yzzz_1[i] * fe_0 + ta2_yz_yyy_zzz_0[i] * fe_0 -
                                ta2_yz_yyy_zzz_1[i] * fe_0 + ta1_z_yyy_yzzz_1[i] + ta2_yz_yyy_yzzz_0[i] * pa_y[i] - ta2_yz_yyy_yzzz_1[i] * pc_y[i];

        ta2_yz_yyyy_zzzz_0[i] = 3.0 * ta2_yz_yy_zzzz_0[i] * fe_0 - 3.0 * ta2_yz_yy_zzzz_1[i] * fe_0 + ta1_z_yyy_zzzz_1[i] +
                                ta2_yz_yyy_zzzz_0[i] * pa_y[i] - ta2_yz_yyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 1065-1080 components of targeted buffer : GG

    auto ta2_yz_yyyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1065);

    auto ta2_yz_yyyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1066);

    auto ta2_yz_yyyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1067);

    auto ta2_yz_yyyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1068);

    auto ta2_yz_yyyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1069);

    auto ta2_yz_yyyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1070);

    auto ta2_yz_yyyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1071);

    auto ta2_yz_yyyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1072);

    auto ta2_yz_yyyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1073);

    auto ta2_yz_yyyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1074);

    auto ta2_yz_yyyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1075);

    auto ta2_yz_yyyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1076);

    auto ta2_yz_yyyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1077);

    auto ta2_yz_yyyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1078);

    auto ta2_yz_yyyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1079);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_yyy_xxxx_1,   \
                             ta1_y_yyy_xxxy_1,   \
                             ta1_y_yyy_xxyy_1,   \
                             ta1_y_yyy_xxyz_1,   \
                             ta1_y_yyy_xyyy_1,   \
                             ta1_y_yyy_xyyz_1,   \
                             ta1_y_yyy_xyzz_1,   \
                             ta1_y_yyy_yyyy_1,   \
                             ta1_y_yyy_yyyz_1,   \
                             ta1_y_yyy_yyzz_1,   \
                             ta1_y_yyy_yzzz_1,   \
                             ta1_z_yyz_xxxz_1,   \
                             ta1_z_yyz_xxzz_1,   \
                             ta1_z_yyz_xzzz_1,   \
                             ta1_z_yyz_zzzz_1,   \
                             ta2_yz_yyy_xxxx_0,  \
                             ta2_yz_yyy_xxxx_1,  \
                             ta2_yz_yyy_xxxy_0,  \
                             ta2_yz_yyy_xxxy_1,  \
                             ta2_yz_yyy_xxy_0,   \
                             ta2_yz_yyy_xxy_1,   \
                             ta2_yz_yyy_xxyy_0,  \
                             ta2_yz_yyy_xxyy_1,  \
                             ta2_yz_yyy_xxyz_0,  \
                             ta2_yz_yyy_xxyz_1,  \
                             ta2_yz_yyy_xyy_0,   \
                             ta2_yz_yyy_xyy_1,   \
                             ta2_yz_yyy_xyyy_0,  \
                             ta2_yz_yyy_xyyy_1,  \
                             ta2_yz_yyy_xyyz_0,  \
                             ta2_yz_yyy_xyyz_1,  \
                             ta2_yz_yyy_xyz_0,   \
                             ta2_yz_yyy_xyz_1,   \
                             ta2_yz_yyy_xyzz_0,  \
                             ta2_yz_yyy_xyzz_1,  \
                             ta2_yz_yyy_yyy_0,   \
                             ta2_yz_yyy_yyy_1,   \
                             ta2_yz_yyy_yyyy_0,  \
                             ta2_yz_yyy_yyyy_1,  \
                             ta2_yz_yyy_yyyz_0,  \
                             ta2_yz_yyy_yyyz_1,  \
                             ta2_yz_yyy_yyz_0,   \
                             ta2_yz_yyy_yyz_1,   \
                             ta2_yz_yyy_yyzz_0,  \
                             ta2_yz_yyy_yyzz_1,  \
                             ta2_yz_yyy_yzz_0,   \
                             ta2_yz_yyy_yzz_1,   \
                             ta2_yz_yyy_yzzz_0,  \
                             ta2_yz_yyy_yzzz_1,  \
                             ta2_yz_yyyz_xxxx_0, \
                             ta2_yz_yyyz_xxxy_0, \
                             ta2_yz_yyyz_xxxz_0, \
                             ta2_yz_yyyz_xxyy_0, \
                             ta2_yz_yyyz_xxyz_0, \
                             ta2_yz_yyyz_xxzz_0, \
                             ta2_yz_yyyz_xyyy_0, \
                             ta2_yz_yyyz_xyyz_0, \
                             ta2_yz_yyyz_xyzz_0, \
                             ta2_yz_yyyz_xzzz_0, \
                             ta2_yz_yyyz_yyyy_0, \
                             ta2_yz_yyyz_yyyz_0, \
                             ta2_yz_yyyz_yyzz_0, \
                             ta2_yz_yyyz_yzzz_0, \
                             ta2_yz_yyyz_zzzz_0, \
                             ta2_yz_yyz_xxxz_0,  \
                             ta2_yz_yyz_xxxz_1,  \
                             ta2_yz_yyz_xxzz_0,  \
                             ta2_yz_yyz_xxzz_1,  \
                             ta2_yz_yyz_xzzz_0,  \
                             ta2_yz_yyz_xzzz_1,  \
                             ta2_yz_yyz_zzzz_0,  \
                             ta2_yz_yyz_zzzz_1,  \
                             ta2_yz_yz_xxxz_0,   \
                             ta2_yz_yz_xxxz_1,   \
                             ta2_yz_yz_xxzz_0,   \
                             ta2_yz_yz_xxzz_1,   \
                             ta2_yz_yz_xzzz_0,   \
                             ta2_yz_yz_xzzz_1,   \
                             ta2_yz_yz_zzzz_0,   \
                             ta2_yz_yz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyyz_xxxx_0[i] = ta1_y_yyy_xxxx_1[i] + ta2_yz_yyy_xxxx_0[i] * pa_z[i] - ta2_yz_yyy_xxxx_1[i] * pc_z[i];

        ta2_yz_yyyz_xxxy_0[i] = ta1_y_yyy_xxxy_1[i] + ta2_yz_yyy_xxxy_0[i] * pa_z[i] - ta2_yz_yyy_xxxy_1[i] * pc_z[i];

        ta2_yz_yyyz_xxxz_0[i] = 2.0 * ta2_yz_yz_xxxz_0[i] * fe_0 - 2.0 * ta2_yz_yz_xxxz_1[i] * fe_0 + ta1_z_yyz_xxxz_1[i] +
                                ta2_yz_yyz_xxxz_0[i] * pa_y[i] - ta2_yz_yyz_xxxz_1[i] * pc_y[i];

        ta2_yz_yyyz_xxyy_0[i] = ta1_y_yyy_xxyy_1[i] + ta2_yz_yyy_xxyy_0[i] * pa_z[i] - ta2_yz_yyy_xxyy_1[i] * pc_z[i];

        ta2_yz_yyyz_xxyz_0[i] = ta2_yz_yyy_xxy_0[i] * fe_0 - ta2_yz_yyy_xxy_1[i] * fe_0 + ta1_y_yyy_xxyz_1[i] + ta2_yz_yyy_xxyz_0[i] * pa_z[i] -
                                ta2_yz_yyy_xxyz_1[i] * pc_z[i];

        ta2_yz_yyyz_xxzz_0[i] = 2.0 * ta2_yz_yz_xxzz_0[i] * fe_0 - 2.0 * ta2_yz_yz_xxzz_1[i] * fe_0 + ta1_z_yyz_xxzz_1[i] +
                                ta2_yz_yyz_xxzz_0[i] * pa_y[i] - ta2_yz_yyz_xxzz_1[i] * pc_y[i];

        ta2_yz_yyyz_xyyy_0[i] = ta1_y_yyy_xyyy_1[i] + ta2_yz_yyy_xyyy_0[i] * pa_z[i] - ta2_yz_yyy_xyyy_1[i] * pc_z[i];

        ta2_yz_yyyz_xyyz_0[i] = ta2_yz_yyy_xyy_0[i] * fe_0 - ta2_yz_yyy_xyy_1[i] * fe_0 + ta1_y_yyy_xyyz_1[i] + ta2_yz_yyy_xyyz_0[i] * pa_z[i] -
                                ta2_yz_yyy_xyyz_1[i] * pc_z[i];

        ta2_yz_yyyz_xyzz_0[i] = 2.0 * ta2_yz_yyy_xyz_0[i] * fe_0 - 2.0 * ta2_yz_yyy_xyz_1[i] * fe_0 + ta1_y_yyy_xyzz_1[i] +
                                ta2_yz_yyy_xyzz_0[i] * pa_z[i] - ta2_yz_yyy_xyzz_1[i] * pc_z[i];

        ta2_yz_yyyz_xzzz_0[i] = 2.0 * ta2_yz_yz_xzzz_0[i] * fe_0 - 2.0 * ta2_yz_yz_xzzz_1[i] * fe_0 + ta1_z_yyz_xzzz_1[i] +
                                ta2_yz_yyz_xzzz_0[i] * pa_y[i] - ta2_yz_yyz_xzzz_1[i] * pc_y[i];

        ta2_yz_yyyz_yyyy_0[i] = ta1_y_yyy_yyyy_1[i] + ta2_yz_yyy_yyyy_0[i] * pa_z[i] - ta2_yz_yyy_yyyy_1[i] * pc_z[i];

        ta2_yz_yyyz_yyyz_0[i] = ta2_yz_yyy_yyy_0[i] * fe_0 - ta2_yz_yyy_yyy_1[i] * fe_0 + ta1_y_yyy_yyyz_1[i] + ta2_yz_yyy_yyyz_0[i] * pa_z[i] -
                                ta2_yz_yyy_yyyz_1[i] * pc_z[i];

        ta2_yz_yyyz_yyzz_0[i] = 2.0 * ta2_yz_yyy_yyz_0[i] * fe_0 - 2.0 * ta2_yz_yyy_yyz_1[i] * fe_0 + ta1_y_yyy_yyzz_1[i] +
                                ta2_yz_yyy_yyzz_0[i] * pa_z[i] - ta2_yz_yyy_yyzz_1[i] * pc_z[i];

        ta2_yz_yyyz_yzzz_0[i] = 3.0 * ta2_yz_yyy_yzz_0[i] * fe_0 - 3.0 * ta2_yz_yyy_yzz_1[i] * fe_0 + ta1_y_yyy_yzzz_1[i] +
                                ta2_yz_yyy_yzzz_0[i] * pa_z[i] - ta2_yz_yyy_yzzz_1[i] * pc_z[i];

        ta2_yz_yyyz_zzzz_0[i] = 2.0 * ta2_yz_yz_zzzz_0[i] * fe_0 - 2.0 * ta2_yz_yz_zzzz_1[i] * fe_0 + ta1_z_yyz_zzzz_1[i] +
                                ta2_yz_yyz_zzzz_0[i] * pa_y[i] - ta2_yz_yyz_zzzz_1[i] * pc_y[i];
    }

    // Set up 1080-1095 components of targeted buffer : GG

    auto ta2_yz_yyzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1080);

    auto ta2_yz_yyzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1081);

    auto ta2_yz_yyzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1082);

    auto ta2_yz_yyzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1083);

    auto ta2_yz_yyzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1084);

    auto ta2_yz_yyzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1085);

    auto ta2_yz_yyzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1086);

    auto ta2_yz_yyzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1087);

    auto ta2_yz_yyzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1088);

    auto ta2_yz_yyzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1089);

    auto ta2_yz_yyzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1090);

    auto ta2_yz_yyzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1091);

    auto ta2_yz_yyzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1092);

    auto ta2_yz_yyzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1093);

    auto ta2_yz_yyzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1094);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_yyz_xxxy_1,   \
                             ta1_y_yyz_xxyy_1,   \
                             ta1_y_yyz_xyyy_1,   \
                             ta1_y_yyz_yyyy_1,   \
                             ta1_z_yzz_xxxx_1,   \
                             ta1_z_yzz_xxxz_1,   \
                             ta1_z_yzz_xxyz_1,   \
                             ta1_z_yzz_xxzz_1,   \
                             ta1_z_yzz_xyyz_1,   \
                             ta1_z_yzz_xyzz_1,   \
                             ta1_z_yzz_xzzz_1,   \
                             ta1_z_yzz_yyyz_1,   \
                             ta1_z_yzz_yyzz_1,   \
                             ta1_z_yzz_yzzz_1,   \
                             ta1_z_yzz_zzzz_1,   \
                             ta2_yz_yy_xxxy_0,   \
                             ta2_yz_yy_xxxy_1,   \
                             ta2_yz_yy_xxyy_0,   \
                             ta2_yz_yy_xxyy_1,   \
                             ta2_yz_yy_xyyy_0,   \
                             ta2_yz_yy_xyyy_1,   \
                             ta2_yz_yy_yyyy_0,   \
                             ta2_yz_yy_yyyy_1,   \
                             ta2_yz_yyz_xxxy_0,  \
                             ta2_yz_yyz_xxxy_1,  \
                             ta2_yz_yyz_xxyy_0,  \
                             ta2_yz_yyz_xxyy_1,  \
                             ta2_yz_yyz_xyyy_0,  \
                             ta2_yz_yyz_xyyy_1,  \
                             ta2_yz_yyz_yyyy_0,  \
                             ta2_yz_yyz_yyyy_1,  \
                             ta2_yz_yyzz_xxxx_0, \
                             ta2_yz_yyzz_xxxy_0, \
                             ta2_yz_yyzz_xxxz_0, \
                             ta2_yz_yyzz_xxyy_0, \
                             ta2_yz_yyzz_xxyz_0, \
                             ta2_yz_yyzz_xxzz_0, \
                             ta2_yz_yyzz_xyyy_0, \
                             ta2_yz_yyzz_xyyz_0, \
                             ta2_yz_yyzz_xyzz_0, \
                             ta2_yz_yyzz_xzzz_0, \
                             ta2_yz_yyzz_yyyy_0, \
                             ta2_yz_yyzz_yyyz_0, \
                             ta2_yz_yyzz_yyzz_0, \
                             ta2_yz_yyzz_yzzz_0, \
                             ta2_yz_yyzz_zzzz_0, \
                             ta2_yz_yzz_xxxx_0,  \
                             ta2_yz_yzz_xxxx_1,  \
                             ta2_yz_yzz_xxxz_0,  \
                             ta2_yz_yzz_xxxz_1,  \
                             ta2_yz_yzz_xxyz_0,  \
                             ta2_yz_yzz_xxyz_1,  \
                             ta2_yz_yzz_xxz_0,   \
                             ta2_yz_yzz_xxz_1,   \
                             ta2_yz_yzz_xxzz_0,  \
                             ta2_yz_yzz_xxzz_1,  \
                             ta2_yz_yzz_xyyz_0,  \
                             ta2_yz_yzz_xyyz_1,  \
                             ta2_yz_yzz_xyz_0,   \
                             ta2_yz_yzz_xyz_1,   \
                             ta2_yz_yzz_xyzz_0,  \
                             ta2_yz_yzz_xyzz_1,  \
                             ta2_yz_yzz_xzz_0,   \
                             ta2_yz_yzz_xzz_1,   \
                             ta2_yz_yzz_xzzz_0,  \
                             ta2_yz_yzz_xzzz_1,  \
                             ta2_yz_yzz_yyyz_0,  \
                             ta2_yz_yzz_yyyz_1,  \
                             ta2_yz_yzz_yyz_0,   \
                             ta2_yz_yzz_yyz_1,   \
                             ta2_yz_yzz_yyzz_0,  \
                             ta2_yz_yzz_yyzz_1,  \
                             ta2_yz_yzz_yzz_0,   \
                             ta2_yz_yzz_yzz_1,   \
                             ta2_yz_yzz_yzzz_0,  \
                             ta2_yz_yzz_yzzz_1,  \
                             ta2_yz_yzz_zzz_0,   \
                             ta2_yz_yzz_zzz_1,   \
                             ta2_yz_yzz_zzzz_0,  \
                             ta2_yz_yzz_zzzz_1,  \
                             ta2_yz_zz_xxxx_0,   \
                             ta2_yz_zz_xxxx_1,   \
                             ta2_yz_zz_xxxz_0,   \
                             ta2_yz_zz_xxxz_1,   \
                             ta2_yz_zz_xxyz_0,   \
                             ta2_yz_zz_xxyz_1,   \
                             ta2_yz_zz_xxzz_0,   \
                             ta2_yz_zz_xxzz_1,   \
                             ta2_yz_zz_xyyz_0,   \
                             ta2_yz_zz_xyyz_1,   \
                             ta2_yz_zz_xyzz_0,   \
                             ta2_yz_zz_xyzz_1,   \
                             ta2_yz_zz_xzzz_0,   \
                             ta2_yz_zz_xzzz_1,   \
                             ta2_yz_zz_yyyz_0,   \
                             ta2_yz_zz_yyyz_1,   \
                             ta2_yz_zz_yyzz_0,   \
                             ta2_yz_zz_yyzz_1,   \
                             ta2_yz_zz_yzzz_0,   \
                             ta2_yz_zz_yzzz_1,   \
                             ta2_yz_zz_zzzz_0,   \
                             ta2_yz_zz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yyzz_xxxx_0[i] = ta2_yz_zz_xxxx_0[i] * fe_0 - ta2_yz_zz_xxxx_1[i] * fe_0 + ta1_z_yzz_xxxx_1[i] + ta2_yz_yzz_xxxx_0[i] * pa_y[i] -
                                ta2_yz_yzz_xxxx_1[i] * pc_y[i];

        ta2_yz_yyzz_xxxy_0[i] = ta2_yz_yy_xxxy_0[i] * fe_0 - ta2_yz_yy_xxxy_1[i] * fe_0 + ta1_y_yyz_xxxy_1[i] + ta2_yz_yyz_xxxy_0[i] * pa_z[i] -
                                ta2_yz_yyz_xxxy_1[i] * pc_z[i];

        ta2_yz_yyzz_xxxz_0[i] = ta2_yz_zz_xxxz_0[i] * fe_0 - ta2_yz_zz_xxxz_1[i] * fe_0 + ta1_z_yzz_xxxz_1[i] + ta2_yz_yzz_xxxz_0[i] * pa_y[i] -
                                ta2_yz_yzz_xxxz_1[i] * pc_y[i];

        ta2_yz_yyzz_xxyy_0[i] = ta2_yz_yy_xxyy_0[i] * fe_0 - ta2_yz_yy_xxyy_1[i] * fe_0 + ta1_y_yyz_xxyy_1[i] + ta2_yz_yyz_xxyy_0[i] * pa_z[i] -
                                ta2_yz_yyz_xxyy_1[i] * pc_z[i];

        ta2_yz_yyzz_xxyz_0[i] = ta2_yz_zz_xxyz_0[i] * fe_0 - ta2_yz_zz_xxyz_1[i] * fe_0 + ta2_yz_yzz_xxz_0[i] * fe_0 - ta2_yz_yzz_xxz_1[i] * fe_0 +
                                ta1_z_yzz_xxyz_1[i] + ta2_yz_yzz_xxyz_0[i] * pa_y[i] - ta2_yz_yzz_xxyz_1[i] * pc_y[i];

        ta2_yz_yyzz_xxzz_0[i] = ta2_yz_zz_xxzz_0[i] * fe_0 - ta2_yz_zz_xxzz_1[i] * fe_0 + ta1_z_yzz_xxzz_1[i] + ta2_yz_yzz_xxzz_0[i] * pa_y[i] -
                                ta2_yz_yzz_xxzz_1[i] * pc_y[i];

        ta2_yz_yyzz_xyyy_0[i] = ta2_yz_yy_xyyy_0[i] * fe_0 - ta2_yz_yy_xyyy_1[i] * fe_0 + ta1_y_yyz_xyyy_1[i] + ta2_yz_yyz_xyyy_0[i] * pa_z[i] -
                                ta2_yz_yyz_xyyy_1[i] * pc_z[i];

        ta2_yz_yyzz_xyyz_0[i] = ta2_yz_zz_xyyz_0[i] * fe_0 - ta2_yz_zz_xyyz_1[i] * fe_0 + 2.0 * ta2_yz_yzz_xyz_0[i] * fe_0 -
                                2.0 * ta2_yz_yzz_xyz_1[i] * fe_0 + ta1_z_yzz_xyyz_1[i] + ta2_yz_yzz_xyyz_0[i] * pa_y[i] -
                                ta2_yz_yzz_xyyz_1[i] * pc_y[i];

        ta2_yz_yyzz_xyzz_0[i] = ta2_yz_zz_xyzz_0[i] * fe_0 - ta2_yz_zz_xyzz_1[i] * fe_0 + ta2_yz_yzz_xzz_0[i] * fe_0 - ta2_yz_yzz_xzz_1[i] * fe_0 +
                                ta1_z_yzz_xyzz_1[i] + ta2_yz_yzz_xyzz_0[i] * pa_y[i] - ta2_yz_yzz_xyzz_1[i] * pc_y[i];

        ta2_yz_yyzz_xzzz_0[i] = ta2_yz_zz_xzzz_0[i] * fe_0 - ta2_yz_zz_xzzz_1[i] * fe_0 + ta1_z_yzz_xzzz_1[i] + ta2_yz_yzz_xzzz_0[i] * pa_y[i] -
                                ta2_yz_yzz_xzzz_1[i] * pc_y[i];

        ta2_yz_yyzz_yyyy_0[i] = ta2_yz_yy_yyyy_0[i] * fe_0 - ta2_yz_yy_yyyy_1[i] * fe_0 + ta1_y_yyz_yyyy_1[i] + ta2_yz_yyz_yyyy_0[i] * pa_z[i] -
                                ta2_yz_yyz_yyyy_1[i] * pc_z[i];

        ta2_yz_yyzz_yyyz_0[i] = ta2_yz_zz_yyyz_0[i] * fe_0 - ta2_yz_zz_yyyz_1[i] * fe_0 + 3.0 * ta2_yz_yzz_yyz_0[i] * fe_0 -
                                3.0 * ta2_yz_yzz_yyz_1[i] * fe_0 + ta1_z_yzz_yyyz_1[i] + ta2_yz_yzz_yyyz_0[i] * pa_y[i] -
                                ta2_yz_yzz_yyyz_1[i] * pc_y[i];

        ta2_yz_yyzz_yyzz_0[i] = ta2_yz_zz_yyzz_0[i] * fe_0 - ta2_yz_zz_yyzz_1[i] * fe_0 + 2.0 * ta2_yz_yzz_yzz_0[i] * fe_0 -
                                2.0 * ta2_yz_yzz_yzz_1[i] * fe_0 + ta1_z_yzz_yyzz_1[i] + ta2_yz_yzz_yyzz_0[i] * pa_y[i] -
                                ta2_yz_yzz_yyzz_1[i] * pc_y[i];

        ta2_yz_yyzz_yzzz_0[i] = ta2_yz_zz_yzzz_0[i] * fe_0 - ta2_yz_zz_yzzz_1[i] * fe_0 + ta2_yz_yzz_zzz_0[i] * fe_0 - ta2_yz_yzz_zzz_1[i] * fe_0 +
                                ta1_z_yzz_yzzz_1[i] + ta2_yz_yzz_yzzz_0[i] * pa_y[i] - ta2_yz_yzz_yzzz_1[i] * pc_y[i];

        ta2_yz_yyzz_zzzz_0[i] = ta2_yz_zz_zzzz_0[i] * fe_0 - ta2_yz_zz_zzzz_1[i] * fe_0 + ta1_z_yzz_zzzz_1[i] + ta2_yz_yzz_zzzz_0[i] * pa_y[i] -
                                ta2_yz_yzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 1095-1110 components of targeted buffer : GG

    auto ta2_yz_yzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1095);

    auto ta2_yz_yzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1096);

    auto ta2_yz_yzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1097);

    auto ta2_yz_yzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1098);

    auto ta2_yz_yzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1099);

    auto ta2_yz_yzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1100);

    auto ta2_yz_yzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1101);

    auto ta2_yz_yzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1102);

    auto ta2_yz_yzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1103);

    auto ta2_yz_yzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1104);

    auto ta2_yz_yzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1105);

    auto ta2_yz_yzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1106);

    auto ta2_yz_yzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1107);

    auto ta2_yz_yzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1108);

    auto ta2_yz_yzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1109);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_z_zzz_xxxx_1,   \
                             ta1_z_zzz_xxxy_1,   \
                             ta1_z_zzz_xxxz_1,   \
                             ta1_z_zzz_xxyy_1,   \
                             ta1_z_zzz_xxyz_1,   \
                             ta1_z_zzz_xxzz_1,   \
                             ta1_z_zzz_xyyy_1,   \
                             ta1_z_zzz_xyyz_1,   \
                             ta1_z_zzz_xyzz_1,   \
                             ta1_z_zzz_xzzz_1,   \
                             ta1_z_zzz_yyyy_1,   \
                             ta1_z_zzz_yyyz_1,   \
                             ta1_z_zzz_yyzz_1,   \
                             ta1_z_zzz_yzzz_1,   \
                             ta1_z_zzz_zzzz_1,   \
                             ta2_yz_yzzz_xxxx_0, \
                             ta2_yz_yzzz_xxxy_0, \
                             ta2_yz_yzzz_xxxz_0, \
                             ta2_yz_yzzz_xxyy_0, \
                             ta2_yz_yzzz_xxyz_0, \
                             ta2_yz_yzzz_xxzz_0, \
                             ta2_yz_yzzz_xyyy_0, \
                             ta2_yz_yzzz_xyyz_0, \
                             ta2_yz_yzzz_xyzz_0, \
                             ta2_yz_yzzz_xzzz_0, \
                             ta2_yz_yzzz_yyyy_0, \
                             ta2_yz_yzzz_yyyz_0, \
                             ta2_yz_yzzz_yyzz_0, \
                             ta2_yz_yzzz_yzzz_0, \
                             ta2_yz_yzzz_zzzz_0, \
                             ta2_yz_zzz_xxx_0,   \
                             ta2_yz_zzz_xxx_1,   \
                             ta2_yz_zzz_xxxx_0,  \
                             ta2_yz_zzz_xxxx_1,  \
                             ta2_yz_zzz_xxxy_0,  \
                             ta2_yz_zzz_xxxy_1,  \
                             ta2_yz_zzz_xxxz_0,  \
                             ta2_yz_zzz_xxxz_1,  \
                             ta2_yz_zzz_xxy_0,   \
                             ta2_yz_zzz_xxy_1,   \
                             ta2_yz_zzz_xxyy_0,  \
                             ta2_yz_zzz_xxyy_1,  \
                             ta2_yz_zzz_xxyz_0,  \
                             ta2_yz_zzz_xxyz_1,  \
                             ta2_yz_zzz_xxz_0,   \
                             ta2_yz_zzz_xxz_1,   \
                             ta2_yz_zzz_xxzz_0,  \
                             ta2_yz_zzz_xxzz_1,  \
                             ta2_yz_zzz_xyy_0,   \
                             ta2_yz_zzz_xyy_1,   \
                             ta2_yz_zzz_xyyy_0,  \
                             ta2_yz_zzz_xyyy_1,  \
                             ta2_yz_zzz_xyyz_0,  \
                             ta2_yz_zzz_xyyz_1,  \
                             ta2_yz_zzz_xyz_0,   \
                             ta2_yz_zzz_xyz_1,   \
                             ta2_yz_zzz_xyzz_0,  \
                             ta2_yz_zzz_xyzz_1,  \
                             ta2_yz_zzz_xzz_0,   \
                             ta2_yz_zzz_xzz_1,   \
                             ta2_yz_zzz_xzzz_0,  \
                             ta2_yz_zzz_xzzz_1,  \
                             ta2_yz_zzz_yyy_0,   \
                             ta2_yz_zzz_yyy_1,   \
                             ta2_yz_zzz_yyyy_0,  \
                             ta2_yz_zzz_yyyy_1,  \
                             ta2_yz_zzz_yyyz_0,  \
                             ta2_yz_zzz_yyyz_1,  \
                             ta2_yz_zzz_yyz_0,   \
                             ta2_yz_zzz_yyz_1,   \
                             ta2_yz_zzz_yyzz_0,  \
                             ta2_yz_zzz_yyzz_1,  \
                             ta2_yz_zzz_yzz_0,   \
                             ta2_yz_zzz_yzz_1,   \
                             ta2_yz_zzz_yzzz_0,  \
                             ta2_yz_zzz_yzzz_1,  \
                             ta2_yz_zzz_zzz_0,   \
                             ta2_yz_zzz_zzz_1,   \
                             ta2_yz_zzz_zzzz_0,  \
                             ta2_yz_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_yzzz_xxxx_0[i] = ta1_z_zzz_xxxx_1[i] + ta2_yz_zzz_xxxx_0[i] * pa_y[i] - ta2_yz_zzz_xxxx_1[i] * pc_y[i];

        ta2_yz_yzzz_xxxy_0[i] = ta2_yz_zzz_xxx_0[i] * fe_0 - ta2_yz_zzz_xxx_1[i] * fe_0 + ta1_z_zzz_xxxy_1[i] + ta2_yz_zzz_xxxy_0[i] * pa_y[i] -
                                ta2_yz_zzz_xxxy_1[i] * pc_y[i];

        ta2_yz_yzzz_xxxz_0[i] = ta1_z_zzz_xxxz_1[i] + ta2_yz_zzz_xxxz_0[i] * pa_y[i] - ta2_yz_zzz_xxxz_1[i] * pc_y[i];

        ta2_yz_yzzz_xxyy_0[i] = 2.0 * ta2_yz_zzz_xxy_0[i] * fe_0 - 2.0 * ta2_yz_zzz_xxy_1[i] * fe_0 + ta1_z_zzz_xxyy_1[i] +
                                ta2_yz_zzz_xxyy_0[i] * pa_y[i] - ta2_yz_zzz_xxyy_1[i] * pc_y[i];

        ta2_yz_yzzz_xxyz_0[i] = ta2_yz_zzz_xxz_0[i] * fe_0 - ta2_yz_zzz_xxz_1[i] * fe_0 + ta1_z_zzz_xxyz_1[i] + ta2_yz_zzz_xxyz_0[i] * pa_y[i] -
                                ta2_yz_zzz_xxyz_1[i] * pc_y[i];

        ta2_yz_yzzz_xxzz_0[i] = ta1_z_zzz_xxzz_1[i] + ta2_yz_zzz_xxzz_0[i] * pa_y[i] - ta2_yz_zzz_xxzz_1[i] * pc_y[i];

        ta2_yz_yzzz_xyyy_0[i] = 3.0 * ta2_yz_zzz_xyy_0[i] * fe_0 - 3.0 * ta2_yz_zzz_xyy_1[i] * fe_0 + ta1_z_zzz_xyyy_1[i] +
                                ta2_yz_zzz_xyyy_0[i] * pa_y[i] - ta2_yz_zzz_xyyy_1[i] * pc_y[i];

        ta2_yz_yzzz_xyyz_0[i] = 2.0 * ta2_yz_zzz_xyz_0[i] * fe_0 - 2.0 * ta2_yz_zzz_xyz_1[i] * fe_0 + ta1_z_zzz_xyyz_1[i] +
                                ta2_yz_zzz_xyyz_0[i] * pa_y[i] - ta2_yz_zzz_xyyz_1[i] * pc_y[i];

        ta2_yz_yzzz_xyzz_0[i] = ta2_yz_zzz_xzz_0[i] * fe_0 - ta2_yz_zzz_xzz_1[i] * fe_0 + ta1_z_zzz_xyzz_1[i] + ta2_yz_zzz_xyzz_0[i] * pa_y[i] -
                                ta2_yz_zzz_xyzz_1[i] * pc_y[i];

        ta2_yz_yzzz_xzzz_0[i] = ta1_z_zzz_xzzz_1[i] + ta2_yz_zzz_xzzz_0[i] * pa_y[i] - ta2_yz_zzz_xzzz_1[i] * pc_y[i];

        ta2_yz_yzzz_yyyy_0[i] = 4.0 * ta2_yz_zzz_yyy_0[i] * fe_0 - 4.0 * ta2_yz_zzz_yyy_1[i] * fe_0 + ta1_z_zzz_yyyy_1[i] +
                                ta2_yz_zzz_yyyy_0[i] * pa_y[i] - ta2_yz_zzz_yyyy_1[i] * pc_y[i];

        ta2_yz_yzzz_yyyz_0[i] = 3.0 * ta2_yz_zzz_yyz_0[i] * fe_0 - 3.0 * ta2_yz_zzz_yyz_1[i] * fe_0 + ta1_z_zzz_yyyz_1[i] +
                                ta2_yz_zzz_yyyz_0[i] * pa_y[i] - ta2_yz_zzz_yyyz_1[i] * pc_y[i];

        ta2_yz_yzzz_yyzz_0[i] = 2.0 * ta2_yz_zzz_yzz_0[i] * fe_0 - 2.0 * ta2_yz_zzz_yzz_1[i] * fe_0 + ta1_z_zzz_yyzz_1[i] +
                                ta2_yz_zzz_yyzz_0[i] * pa_y[i] - ta2_yz_zzz_yyzz_1[i] * pc_y[i];

        ta2_yz_yzzz_yzzz_0[i] = ta2_yz_zzz_zzz_0[i] * fe_0 - ta2_yz_zzz_zzz_1[i] * fe_0 + ta1_z_zzz_yzzz_1[i] + ta2_yz_zzz_yzzz_0[i] * pa_y[i] -
                                ta2_yz_zzz_yzzz_1[i] * pc_y[i];

        ta2_yz_yzzz_zzzz_0[i] = ta1_z_zzz_zzzz_1[i] + ta2_yz_zzz_zzzz_0[i] * pa_y[i] - ta2_yz_zzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 1110-1125 components of targeted buffer : GG

    auto ta2_yz_zzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1110);

    auto ta2_yz_zzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1111);

    auto ta2_yz_zzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1112);

    auto ta2_yz_zzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1113);

    auto ta2_yz_zzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1114);

    auto ta2_yz_zzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1115);

    auto ta2_yz_zzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1116);

    auto ta2_yz_zzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1117);

    auto ta2_yz_zzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1118);

    auto ta2_yz_zzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1119);

    auto ta2_yz_zzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1120);

    auto ta2_yz_zzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1121);

    auto ta2_yz_zzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1122);

    auto ta2_yz_zzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1123);

    auto ta2_yz_zzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1124);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_y_zzz_xxxx_1,   \
                             ta1_y_zzz_xxxy_1,   \
                             ta1_y_zzz_xxxz_1,   \
                             ta1_y_zzz_xxyy_1,   \
                             ta1_y_zzz_xxyz_1,   \
                             ta1_y_zzz_xxzz_1,   \
                             ta1_y_zzz_xyyy_1,   \
                             ta1_y_zzz_xyyz_1,   \
                             ta1_y_zzz_xyzz_1,   \
                             ta1_y_zzz_xzzz_1,   \
                             ta1_y_zzz_yyyy_1,   \
                             ta1_y_zzz_yyyz_1,   \
                             ta1_y_zzz_yyzz_1,   \
                             ta1_y_zzz_yzzz_1,   \
                             ta1_y_zzz_zzzz_1,   \
                             ta2_yz_zz_xxxx_0,   \
                             ta2_yz_zz_xxxx_1,   \
                             ta2_yz_zz_xxxy_0,   \
                             ta2_yz_zz_xxxy_1,   \
                             ta2_yz_zz_xxxz_0,   \
                             ta2_yz_zz_xxxz_1,   \
                             ta2_yz_zz_xxyy_0,   \
                             ta2_yz_zz_xxyy_1,   \
                             ta2_yz_zz_xxyz_0,   \
                             ta2_yz_zz_xxyz_1,   \
                             ta2_yz_zz_xxzz_0,   \
                             ta2_yz_zz_xxzz_1,   \
                             ta2_yz_zz_xyyy_0,   \
                             ta2_yz_zz_xyyy_1,   \
                             ta2_yz_zz_xyyz_0,   \
                             ta2_yz_zz_xyyz_1,   \
                             ta2_yz_zz_xyzz_0,   \
                             ta2_yz_zz_xyzz_1,   \
                             ta2_yz_zz_xzzz_0,   \
                             ta2_yz_zz_xzzz_1,   \
                             ta2_yz_zz_yyyy_0,   \
                             ta2_yz_zz_yyyy_1,   \
                             ta2_yz_zz_yyyz_0,   \
                             ta2_yz_zz_yyyz_1,   \
                             ta2_yz_zz_yyzz_0,   \
                             ta2_yz_zz_yyzz_1,   \
                             ta2_yz_zz_yzzz_0,   \
                             ta2_yz_zz_yzzz_1,   \
                             ta2_yz_zz_zzzz_0,   \
                             ta2_yz_zz_zzzz_1,   \
                             ta2_yz_zzz_xxx_0,   \
                             ta2_yz_zzz_xxx_1,   \
                             ta2_yz_zzz_xxxx_0,  \
                             ta2_yz_zzz_xxxx_1,  \
                             ta2_yz_zzz_xxxy_0,  \
                             ta2_yz_zzz_xxxy_1,  \
                             ta2_yz_zzz_xxxz_0,  \
                             ta2_yz_zzz_xxxz_1,  \
                             ta2_yz_zzz_xxy_0,   \
                             ta2_yz_zzz_xxy_1,   \
                             ta2_yz_zzz_xxyy_0,  \
                             ta2_yz_zzz_xxyy_1,  \
                             ta2_yz_zzz_xxyz_0,  \
                             ta2_yz_zzz_xxyz_1,  \
                             ta2_yz_zzz_xxz_0,   \
                             ta2_yz_zzz_xxz_1,   \
                             ta2_yz_zzz_xxzz_0,  \
                             ta2_yz_zzz_xxzz_1,  \
                             ta2_yz_zzz_xyy_0,   \
                             ta2_yz_zzz_xyy_1,   \
                             ta2_yz_zzz_xyyy_0,  \
                             ta2_yz_zzz_xyyy_1,  \
                             ta2_yz_zzz_xyyz_0,  \
                             ta2_yz_zzz_xyyz_1,  \
                             ta2_yz_zzz_xyz_0,   \
                             ta2_yz_zzz_xyz_1,   \
                             ta2_yz_zzz_xyzz_0,  \
                             ta2_yz_zzz_xyzz_1,  \
                             ta2_yz_zzz_xzz_0,   \
                             ta2_yz_zzz_xzz_1,   \
                             ta2_yz_zzz_xzzz_0,  \
                             ta2_yz_zzz_xzzz_1,  \
                             ta2_yz_zzz_yyy_0,   \
                             ta2_yz_zzz_yyy_1,   \
                             ta2_yz_zzz_yyyy_0,  \
                             ta2_yz_zzz_yyyy_1,  \
                             ta2_yz_zzz_yyyz_0,  \
                             ta2_yz_zzz_yyyz_1,  \
                             ta2_yz_zzz_yyz_0,   \
                             ta2_yz_zzz_yyz_1,   \
                             ta2_yz_zzz_yyzz_0,  \
                             ta2_yz_zzz_yyzz_1,  \
                             ta2_yz_zzz_yzz_0,   \
                             ta2_yz_zzz_yzz_1,   \
                             ta2_yz_zzz_yzzz_0,  \
                             ta2_yz_zzz_yzzz_1,  \
                             ta2_yz_zzz_zzz_0,   \
                             ta2_yz_zzz_zzz_1,   \
                             ta2_yz_zzz_zzzz_0,  \
                             ta2_yz_zzz_zzzz_1,  \
                             ta2_yz_zzzz_xxxx_0, \
                             ta2_yz_zzzz_xxxy_0, \
                             ta2_yz_zzzz_xxxz_0, \
                             ta2_yz_zzzz_xxyy_0, \
                             ta2_yz_zzzz_xxyz_0, \
                             ta2_yz_zzzz_xxzz_0, \
                             ta2_yz_zzzz_xyyy_0, \
                             ta2_yz_zzzz_xyyz_0, \
                             ta2_yz_zzzz_xyzz_0, \
                             ta2_yz_zzzz_xzzz_0, \
                             ta2_yz_zzzz_yyyy_0, \
                             ta2_yz_zzzz_yyyz_0, \
                             ta2_yz_zzzz_yyzz_0, \
                             ta2_yz_zzzz_yzzz_0, \
                             ta2_yz_zzzz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_yz_zzzz_xxxx_0[i] = 3.0 * ta2_yz_zz_xxxx_0[i] * fe_0 - 3.0 * ta2_yz_zz_xxxx_1[i] * fe_0 + ta1_y_zzz_xxxx_1[i] +
                                ta2_yz_zzz_xxxx_0[i] * pa_z[i] - ta2_yz_zzz_xxxx_1[i] * pc_z[i];

        ta2_yz_zzzz_xxxy_0[i] = 3.0 * ta2_yz_zz_xxxy_0[i] * fe_0 - 3.0 * ta2_yz_zz_xxxy_1[i] * fe_0 + ta1_y_zzz_xxxy_1[i] +
                                ta2_yz_zzz_xxxy_0[i] * pa_z[i] - ta2_yz_zzz_xxxy_1[i] * pc_z[i];

        ta2_yz_zzzz_xxxz_0[i] = 3.0 * ta2_yz_zz_xxxz_0[i] * fe_0 - 3.0 * ta2_yz_zz_xxxz_1[i] * fe_0 + ta2_yz_zzz_xxx_0[i] * fe_0 -
                                ta2_yz_zzz_xxx_1[i] * fe_0 + ta1_y_zzz_xxxz_1[i] + ta2_yz_zzz_xxxz_0[i] * pa_z[i] - ta2_yz_zzz_xxxz_1[i] * pc_z[i];

        ta2_yz_zzzz_xxyy_0[i] = 3.0 * ta2_yz_zz_xxyy_0[i] * fe_0 - 3.0 * ta2_yz_zz_xxyy_1[i] * fe_0 + ta1_y_zzz_xxyy_1[i] +
                                ta2_yz_zzz_xxyy_0[i] * pa_z[i] - ta2_yz_zzz_xxyy_1[i] * pc_z[i];

        ta2_yz_zzzz_xxyz_0[i] = 3.0 * ta2_yz_zz_xxyz_0[i] * fe_0 - 3.0 * ta2_yz_zz_xxyz_1[i] * fe_0 + ta2_yz_zzz_xxy_0[i] * fe_0 -
                                ta2_yz_zzz_xxy_1[i] * fe_0 + ta1_y_zzz_xxyz_1[i] + ta2_yz_zzz_xxyz_0[i] * pa_z[i] - ta2_yz_zzz_xxyz_1[i] * pc_z[i];

        ta2_yz_zzzz_xxzz_0[i] = 3.0 * ta2_yz_zz_xxzz_0[i] * fe_0 - 3.0 * ta2_yz_zz_xxzz_1[i] * fe_0 + 2.0 * ta2_yz_zzz_xxz_0[i] * fe_0 -
                                2.0 * ta2_yz_zzz_xxz_1[i] * fe_0 + ta1_y_zzz_xxzz_1[i] + ta2_yz_zzz_xxzz_0[i] * pa_z[i] -
                                ta2_yz_zzz_xxzz_1[i] * pc_z[i];

        ta2_yz_zzzz_xyyy_0[i] = 3.0 * ta2_yz_zz_xyyy_0[i] * fe_0 - 3.0 * ta2_yz_zz_xyyy_1[i] * fe_0 + ta1_y_zzz_xyyy_1[i] +
                                ta2_yz_zzz_xyyy_0[i] * pa_z[i] - ta2_yz_zzz_xyyy_1[i] * pc_z[i];

        ta2_yz_zzzz_xyyz_0[i] = 3.0 * ta2_yz_zz_xyyz_0[i] * fe_0 - 3.0 * ta2_yz_zz_xyyz_1[i] * fe_0 + ta2_yz_zzz_xyy_0[i] * fe_0 -
                                ta2_yz_zzz_xyy_1[i] * fe_0 + ta1_y_zzz_xyyz_1[i] + ta2_yz_zzz_xyyz_0[i] * pa_z[i] - ta2_yz_zzz_xyyz_1[i] * pc_z[i];

        ta2_yz_zzzz_xyzz_0[i] = 3.0 * ta2_yz_zz_xyzz_0[i] * fe_0 - 3.0 * ta2_yz_zz_xyzz_1[i] * fe_0 + 2.0 * ta2_yz_zzz_xyz_0[i] * fe_0 -
                                2.0 * ta2_yz_zzz_xyz_1[i] * fe_0 + ta1_y_zzz_xyzz_1[i] + ta2_yz_zzz_xyzz_0[i] * pa_z[i] -
                                ta2_yz_zzz_xyzz_1[i] * pc_z[i];

        ta2_yz_zzzz_xzzz_0[i] = 3.0 * ta2_yz_zz_xzzz_0[i] * fe_0 - 3.0 * ta2_yz_zz_xzzz_1[i] * fe_0 + 3.0 * ta2_yz_zzz_xzz_0[i] * fe_0 -
                                3.0 * ta2_yz_zzz_xzz_1[i] * fe_0 + ta1_y_zzz_xzzz_1[i] + ta2_yz_zzz_xzzz_0[i] * pa_z[i] -
                                ta2_yz_zzz_xzzz_1[i] * pc_z[i];

        ta2_yz_zzzz_yyyy_0[i] = 3.0 * ta2_yz_zz_yyyy_0[i] * fe_0 - 3.0 * ta2_yz_zz_yyyy_1[i] * fe_0 + ta1_y_zzz_yyyy_1[i] +
                                ta2_yz_zzz_yyyy_0[i] * pa_z[i] - ta2_yz_zzz_yyyy_1[i] * pc_z[i];

        ta2_yz_zzzz_yyyz_0[i] = 3.0 * ta2_yz_zz_yyyz_0[i] * fe_0 - 3.0 * ta2_yz_zz_yyyz_1[i] * fe_0 + ta2_yz_zzz_yyy_0[i] * fe_0 -
                                ta2_yz_zzz_yyy_1[i] * fe_0 + ta1_y_zzz_yyyz_1[i] + ta2_yz_zzz_yyyz_0[i] * pa_z[i] - ta2_yz_zzz_yyyz_1[i] * pc_z[i];

        ta2_yz_zzzz_yyzz_0[i] = 3.0 * ta2_yz_zz_yyzz_0[i] * fe_0 - 3.0 * ta2_yz_zz_yyzz_1[i] * fe_0 + 2.0 * ta2_yz_zzz_yyz_0[i] * fe_0 -
                                2.0 * ta2_yz_zzz_yyz_1[i] * fe_0 + ta1_y_zzz_yyzz_1[i] + ta2_yz_zzz_yyzz_0[i] * pa_z[i] -
                                ta2_yz_zzz_yyzz_1[i] * pc_z[i];

        ta2_yz_zzzz_yzzz_0[i] = 3.0 * ta2_yz_zz_yzzz_0[i] * fe_0 - 3.0 * ta2_yz_zz_yzzz_1[i] * fe_0 + 3.0 * ta2_yz_zzz_yzz_0[i] * fe_0 -
                                3.0 * ta2_yz_zzz_yzz_1[i] * fe_0 + ta1_y_zzz_yzzz_1[i] + ta2_yz_zzz_yzzz_0[i] * pa_z[i] -
                                ta2_yz_zzz_yzzz_1[i] * pc_z[i];

        ta2_yz_zzzz_zzzz_0[i] = 3.0 * ta2_yz_zz_zzzz_0[i] * fe_0 - 3.0 * ta2_yz_zz_zzzz_1[i] * fe_0 + 4.0 * ta2_yz_zzz_zzz_0[i] * fe_0 -
                                4.0 * ta2_yz_zzz_zzz_1[i] * fe_0 + ta1_y_zzz_zzzz_1[i] + ta2_yz_zzz_zzzz_0[i] * pa_z[i] -
                                ta2_yz_zzz_zzzz_1[i] * pc_z[i];
    }

    // Set up 1125-1140 components of targeted buffer : GG

    auto ta2_zz_xxxx_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1125);

    auto ta2_zz_xxxx_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1126);

    auto ta2_zz_xxxx_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1127);

    auto ta2_zz_xxxx_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1128);

    auto ta2_zz_xxxx_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1129);

    auto ta2_zz_xxxx_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1130);

    auto ta2_zz_xxxx_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1131);

    auto ta2_zz_xxxx_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1132);

    auto ta2_zz_xxxx_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1133);

    auto ta2_zz_xxxx_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1134);

    auto ta2_zz_xxxx_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1135);

    auto ta2_zz_xxxx_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1136);

    auto ta2_zz_xxxx_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1137);

    auto ta2_zz_xxxx_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1138);

    auto ta2_zz_xxxx_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1139);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta2_zz_xx_xxxx_0,   \
                             ta2_zz_xx_xxxx_1,   \
                             ta2_zz_xx_xxxy_0,   \
                             ta2_zz_xx_xxxy_1,   \
                             ta2_zz_xx_xxxz_0,   \
                             ta2_zz_xx_xxxz_1,   \
                             ta2_zz_xx_xxyy_0,   \
                             ta2_zz_xx_xxyy_1,   \
                             ta2_zz_xx_xxyz_0,   \
                             ta2_zz_xx_xxyz_1,   \
                             ta2_zz_xx_xxzz_0,   \
                             ta2_zz_xx_xxzz_1,   \
                             ta2_zz_xx_xyyy_0,   \
                             ta2_zz_xx_xyyy_1,   \
                             ta2_zz_xx_xyyz_0,   \
                             ta2_zz_xx_xyyz_1,   \
                             ta2_zz_xx_xyzz_0,   \
                             ta2_zz_xx_xyzz_1,   \
                             ta2_zz_xx_xzzz_0,   \
                             ta2_zz_xx_xzzz_1,   \
                             ta2_zz_xx_yyyy_0,   \
                             ta2_zz_xx_yyyy_1,   \
                             ta2_zz_xx_yyyz_0,   \
                             ta2_zz_xx_yyyz_1,   \
                             ta2_zz_xx_yyzz_0,   \
                             ta2_zz_xx_yyzz_1,   \
                             ta2_zz_xx_yzzz_0,   \
                             ta2_zz_xx_yzzz_1,   \
                             ta2_zz_xx_zzzz_0,   \
                             ta2_zz_xx_zzzz_1,   \
                             ta2_zz_xxx_xxx_0,   \
                             ta2_zz_xxx_xxx_1,   \
                             ta2_zz_xxx_xxxx_0,  \
                             ta2_zz_xxx_xxxx_1,  \
                             ta2_zz_xxx_xxxy_0,  \
                             ta2_zz_xxx_xxxy_1,  \
                             ta2_zz_xxx_xxxz_0,  \
                             ta2_zz_xxx_xxxz_1,  \
                             ta2_zz_xxx_xxy_0,   \
                             ta2_zz_xxx_xxy_1,   \
                             ta2_zz_xxx_xxyy_0,  \
                             ta2_zz_xxx_xxyy_1,  \
                             ta2_zz_xxx_xxyz_0,  \
                             ta2_zz_xxx_xxyz_1,  \
                             ta2_zz_xxx_xxz_0,   \
                             ta2_zz_xxx_xxz_1,   \
                             ta2_zz_xxx_xxzz_0,  \
                             ta2_zz_xxx_xxzz_1,  \
                             ta2_zz_xxx_xyy_0,   \
                             ta2_zz_xxx_xyy_1,   \
                             ta2_zz_xxx_xyyy_0,  \
                             ta2_zz_xxx_xyyy_1,  \
                             ta2_zz_xxx_xyyz_0,  \
                             ta2_zz_xxx_xyyz_1,  \
                             ta2_zz_xxx_xyz_0,   \
                             ta2_zz_xxx_xyz_1,   \
                             ta2_zz_xxx_xyzz_0,  \
                             ta2_zz_xxx_xyzz_1,  \
                             ta2_zz_xxx_xzz_0,   \
                             ta2_zz_xxx_xzz_1,   \
                             ta2_zz_xxx_xzzz_0,  \
                             ta2_zz_xxx_xzzz_1,  \
                             ta2_zz_xxx_yyy_0,   \
                             ta2_zz_xxx_yyy_1,   \
                             ta2_zz_xxx_yyyy_0,  \
                             ta2_zz_xxx_yyyy_1,  \
                             ta2_zz_xxx_yyyz_0,  \
                             ta2_zz_xxx_yyyz_1,  \
                             ta2_zz_xxx_yyz_0,   \
                             ta2_zz_xxx_yyz_1,   \
                             ta2_zz_xxx_yyzz_0,  \
                             ta2_zz_xxx_yyzz_1,  \
                             ta2_zz_xxx_yzz_0,   \
                             ta2_zz_xxx_yzz_1,   \
                             ta2_zz_xxx_yzzz_0,  \
                             ta2_zz_xxx_yzzz_1,  \
                             ta2_zz_xxx_zzz_0,   \
                             ta2_zz_xxx_zzz_1,   \
                             ta2_zz_xxx_zzzz_0,  \
                             ta2_zz_xxx_zzzz_1,  \
                             ta2_zz_xxxx_xxxx_0, \
                             ta2_zz_xxxx_xxxy_0, \
                             ta2_zz_xxxx_xxxz_0, \
                             ta2_zz_xxxx_xxyy_0, \
                             ta2_zz_xxxx_xxyz_0, \
                             ta2_zz_xxxx_xxzz_0, \
                             ta2_zz_xxxx_xyyy_0, \
                             ta2_zz_xxxx_xyyz_0, \
                             ta2_zz_xxxx_xyzz_0, \
                             ta2_zz_xxxx_xzzz_0, \
                             ta2_zz_xxxx_yyyy_0, \
                             ta2_zz_xxxx_yyyz_0, \
                             ta2_zz_xxxx_yyzz_0, \
                             ta2_zz_xxxx_yzzz_0, \
                             ta2_zz_xxxx_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxxx_xxxx_0[i] = 3.0 * ta2_zz_xx_xxxx_0[i] * fe_0 - 3.0 * ta2_zz_xx_xxxx_1[i] * fe_0 + 4.0 * ta2_zz_xxx_xxx_0[i] * fe_0 -
                                4.0 * ta2_zz_xxx_xxx_1[i] * fe_0 + ta2_zz_xxx_xxxx_0[i] * pa_x[i] - ta2_zz_xxx_xxxx_1[i] * pc_x[i];

        ta2_zz_xxxx_xxxy_0[i] = 3.0 * ta2_zz_xx_xxxy_0[i] * fe_0 - 3.0 * ta2_zz_xx_xxxy_1[i] * fe_0 + 3.0 * ta2_zz_xxx_xxy_0[i] * fe_0 -
                                3.0 * ta2_zz_xxx_xxy_1[i] * fe_0 + ta2_zz_xxx_xxxy_0[i] * pa_x[i] - ta2_zz_xxx_xxxy_1[i] * pc_x[i];

        ta2_zz_xxxx_xxxz_0[i] = 3.0 * ta2_zz_xx_xxxz_0[i] * fe_0 - 3.0 * ta2_zz_xx_xxxz_1[i] * fe_0 + 3.0 * ta2_zz_xxx_xxz_0[i] * fe_0 -
                                3.0 * ta2_zz_xxx_xxz_1[i] * fe_0 + ta2_zz_xxx_xxxz_0[i] * pa_x[i] - ta2_zz_xxx_xxxz_1[i] * pc_x[i];

        ta2_zz_xxxx_xxyy_0[i] = 3.0 * ta2_zz_xx_xxyy_0[i] * fe_0 - 3.0 * ta2_zz_xx_xxyy_1[i] * fe_0 + 2.0 * ta2_zz_xxx_xyy_0[i] * fe_0 -
                                2.0 * ta2_zz_xxx_xyy_1[i] * fe_0 + ta2_zz_xxx_xxyy_0[i] * pa_x[i] - ta2_zz_xxx_xxyy_1[i] * pc_x[i];

        ta2_zz_xxxx_xxyz_0[i] = 3.0 * ta2_zz_xx_xxyz_0[i] * fe_0 - 3.0 * ta2_zz_xx_xxyz_1[i] * fe_0 + 2.0 * ta2_zz_xxx_xyz_0[i] * fe_0 -
                                2.0 * ta2_zz_xxx_xyz_1[i] * fe_0 + ta2_zz_xxx_xxyz_0[i] * pa_x[i] - ta2_zz_xxx_xxyz_1[i] * pc_x[i];

        ta2_zz_xxxx_xxzz_0[i] = 3.0 * ta2_zz_xx_xxzz_0[i] * fe_0 - 3.0 * ta2_zz_xx_xxzz_1[i] * fe_0 + 2.0 * ta2_zz_xxx_xzz_0[i] * fe_0 -
                                2.0 * ta2_zz_xxx_xzz_1[i] * fe_0 + ta2_zz_xxx_xxzz_0[i] * pa_x[i] - ta2_zz_xxx_xxzz_1[i] * pc_x[i];

        ta2_zz_xxxx_xyyy_0[i] = 3.0 * ta2_zz_xx_xyyy_0[i] * fe_0 - 3.0 * ta2_zz_xx_xyyy_1[i] * fe_0 + ta2_zz_xxx_yyy_0[i] * fe_0 -
                                ta2_zz_xxx_yyy_1[i] * fe_0 + ta2_zz_xxx_xyyy_0[i] * pa_x[i] - ta2_zz_xxx_xyyy_1[i] * pc_x[i];

        ta2_zz_xxxx_xyyz_0[i] = 3.0 * ta2_zz_xx_xyyz_0[i] * fe_0 - 3.0 * ta2_zz_xx_xyyz_1[i] * fe_0 + ta2_zz_xxx_yyz_0[i] * fe_0 -
                                ta2_zz_xxx_yyz_1[i] * fe_0 + ta2_zz_xxx_xyyz_0[i] * pa_x[i] - ta2_zz_xxx_xyyz_1[i] * pc_x[i];

        ta2_zz_xxxx_xyzz_0[i] = 3.0 * ta2_zz_xx_xyzz_0[i] * fe_0 - 3.0 * ta2_zz_xx_xyzz_1[i] * fe_0 + ta2_zz_xxx_yzz_0[i] * fe_0 -
                                ta2_zz_xxx_yzz_1[i] * fe_0 + ta2_zz_xxx_xyzz_0[i] * pa_x[i] - ta2_zz_xxx_xyzz_1[i] * pc_x[i];

        ta2_zz_xxxx_xzzz_0[i] = 3.0 * ta2_zz_xx_xzzz_0[i] * fe_0 - 3.0 * ta2_zz_xx_xzzz_1[i] * fe_0 + ta2_zz_xxx_zzz_0[i] * fe_0 -
                                ta2_zz_xxx_zzz_1[i] * fe_0 + ta2_zz_xxx_xzzz_0[i] * pa_x[i] - ta2_zz_xxx_xzzz_1[i] * pc_x[i];

        ta2_zz_xxxx_yyyy_0[i] =
            3.0 * ta2_zz_xx_yyyy_0[i] * fe_0 - 3.0 * ta2_zz_xx_yyyy_1[i] * fe_0 + ta2_zz_xxx_yyyy_0[i] * pa_x[i] - ta2_zz_xxx_yyyy_1[i] * pc_x[i];

        ta2_zz_xxxx_yyyz_0[i] =
            3.0 * ta2_zz_xx_yyyz_0[i] * fe_0 - 3.0 * ta2_zz_xx_yyyz_1[i] * fe_0 + ta2_zz_xxx_yyyz_0[i] * pa_x[i] - ta2_zz_xxx_yyyz_1[i] * pc_x[i];

        ta2_zz_xxxx_yyzz_0[i] =
            3.0 * ta2_zz_xx_yyzz_0[i] * fe_0 - 3.0 * ta2_zz_xx_yyzz_1[i] * fe_0 + ta2_zz_xxx_yyzz_0[i] * pa_x[i] - ta2_zz_xxx_yyzz_1[i] * pc_x[i];

        ta2_zz_xxxx_yzzz_0[i] =
            3.0 * ta2_zz_xx_yzzz_0[i] * fe_0 - 3.0 * ta2_zz_xx_yzzz_1[i] * fe_0 + ta2_zz_xxx_yzzz_0[i] * pa_x[i] - ta2_zz_xxx_yzzz_1[i] * pc_x[i];

        ta2_zz_xxxx_zzzz_0[i] =
            3.0 * ta2_zz_xx_zzzz_0[i] * fe_0 - 3.0 * ta2_zz_xx_zzzz_1[i] * fe_0 + ta2_zz_xxx_zzzz_0[i] * pa_x[i] - ta2_zz_xxx_zzzz_1[i] * pc_x[i];
    }

    // Set up 1140-1155 components of targeted buffer : GG

    auto ta2_zz_xxxy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1140);

    auto ta2_zz_xxxy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1141);

    auto ta2_zz_xxxy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1142);

    auto ta2_zz_xxxy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1143);

    auto ta2_zz_xxxy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1144);

    auto ta2_zz_xxxy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1145);

    auto ta2_zz_xxxy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1146);

    auto ta2_zz_xxxy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1147);

    auto ta2_zz_xxxy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1148);

    auto ta2_zz_xxxy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1149);

    auto ta2_zz_xxxy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1150);

    auto ta2_zz_xxxy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1151);

    auto ta2_zz_xxxy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1152);

    auto ta2_zz_xxxy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1153);

    auto ta2_zz_xxxy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1154);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta2_zz_xxx_xxx_0,   \
                             ta2_zz_xxx_xxx_1,   \
                             ta2_zz_xxx_xxxx_0,  \
                             ta2_zz_xxx_xxxx_1,  \
                             ta2_zz_xxx_xxxy_0,  \
                             ta2_zz_xxx_xxxy_1,  \
                             ta2_zz_xxx_xxxz_0,  \
                             ta2_zz_xxx_xxxz_1,  \
                             ta2_zz_xxx_xxy_0,   \
                             ta2_zz_xxx_xxy_1,   \
                             ta2_zz_xxx_xxyy_0,  \
                             ta2_zz_xxx_xxyy_1,  \
                             ta2_zz_xxx_xxyz_0,  \
                             ta2_zz_xxx_xxyz_1,  \
                             ta2_zz_xxx_xxz_0,   \
                             ta2_zz_xxx_xxz_1,   \
                             ta2_zz_xxx_xxzz_0,  \
                             ta2_zz_xxx_xxzz_1,  \
                             ta2_zz_xxx_xyy_0,   \
                             ta2_zz_xxx_xyy_1,   \
                             ta2_zz_xxx_xyyy_0,  \
                             ta2_zz_xxx_xyyy_1,  \
                             ta2_zz_xxx_xyyz_0,  \
                             ta2_zz_xxx_xyyz_1,  \
                             ta2_zz_xxx_xyz_0,   \
                             ta2_zz_xxx_xyz_1,   \
                             ta2_zz_xxx_xyzz_0,  \
                             ta2_zz_xxx_xyzz_1,  \
                             ta2_zz_xxx_xzz_0,   \
                             ta2_zz_xxx_xzz_1,   \
                             ta2_zz_xxx_xzzz_0,  \
                             ta2_zz_xxx_xzzz_1,  \
                             ta2_zz_xxx_zzzz_0,  \
                             ta2_zz_xxx_zzzz_1,  \
                             ta2_zz_xxxy_xxxx_0, \
                             ta2_zz_xxxy_xxxy_0, \
                             ta2_zz_xxxy_xxxz_0, \
                             ta2_zz_xxxy_xxyy_0, \
                             ta2_zz_xxxy_xxyz_0, \
                             ta2_zz_xxxy_xxzz_0, \
                             ta2_zz_xxxy_xyyy_0, \
                             ta2_zz_xxxy_xyyz_0, \
                             ta2_zz_xxxy_xyzz_0, \
                             ta2_zz_xxxy_xzzz_0, \
                             ta2_zz_xxxy_yyyy_0, \
                             ta2_zz_xxxy_yyyz_0, \
                             ta2_zz_xxxy_yyzz_0, \
                             ta2_zz_xxxy_yzzz_0, \
                             ta2_zz_xxxy_zzzz_0, \
                             ta2_zz_xxy_yyyy_0,  \
                             ta2_zz_xxy_yyyy_1,  \
                             ta2_zz_xxy_yyyz_0,  \
                             ta2_zz_xxy_yyyz_1,  \
                             ta2_zz_xxy_yyzz_0,  \
                             ta2_zz_xxy_yyzz_1,  \
                             ta2_zz_xxy_yzzz_0,  \
                             ta2_zz_xxy_yzzz_1,  \
                             ta2_zz_xy_yyyy_0,   \
                             ta2_zz_xy_yyyy_1,   \
                             ta2_zz_xy_yyyz_0,   \
                             ta2_zz_xy_yyyz_1,   \
                             ta2_zz_xy_yyzz_0,   \
                             ta2_zz_xy_yyzz_1,   \
                             ta2_zz_xy_yzzz_0,   \
                             ta2_zz_xy_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxxy_xxxx_0[i] = ta2_zz_xxx_xxxx_0[i] * pa_y[i] - ta2_zz_xxx_xxxx_1[i] * pc_y[i];

        ta2_zz_xxxy_xxxy_0[i] =
            ta2_zz_xxx_xxx_0[i] * fe_0 - ta2_zz_xxx_xxx_1[i] * fe_0 + ta2_zz_xxx_xxxy_0[i] * pa_y[i] - ta2_zz_xxx_xxxy_1[i] * pc_y[i];

        ta2_zz_xxxy_xxxz_0[i] = ta2_zz_xxx_xxxz_0[i] * pa_y[i] - ta2_zz_xxx_xxxz_1[i] * pc_y[i];

        ta2_zz_xxxy_xxyy_0[i] =
            2.0 * ta2_zz_xxx_xxy_0[i] * fe_0 - 2.0 * ta2_zz_xxx_xxy_1[i] * fe_0 + ta2_zz_xxx_xxyy_0[i] * pa_y[i] - ta2_zz_xxx_xxyy_1[i] * pc_y[i];

        ta2_zz_xxxy_xxyz_0[i] =
            ta2_zz_xxx_xxz_0[i] * fe_0 - ta2_zz_xxx_xxz_1[i] * fe_0 + ta2_zz_xxx_xxyz_0[i] * pa_y[i] - ta2_zz_xxx_xxyz_1[i] * pc_y[i];

        ta2_zz_xxxy_xxzz_0[i] = ta2_zz_xxx_xxzz_0[i] * pa_y[i] - ta2_zz_xxx_xxzz_1[i] * pc_y[i];

        ta2_zz_xxxy_xyyy_0[i] =
            3.0 * ta2_zz_xxx_xyy_0[i] * fe_0 - 3.0 * ta2_zz_xxx_xyy_1[i] * fe_0 + ta2_zz_xxx_xyyy_0[i] * pa_y[i] - ta2_zz_xxx_xyyy_1[i] * pc_y[i];

        ta2_zz_xxxy_xyyz_0[i] =
            2.0 * ta2_zz_xxx_xyz_0[i] * fe_0 - 2.0 * ta2_zz_xxx_xyz_1[i] * fe_0 + ta2_zz_xxx_xyyz_0[i] * pa_y[i] - ta2_zz_xxx_xyyz_1[i] * pc_y[i];

        ta2_zz_xxxy_xyzz_0[i] =
            ta2_zz_xxx_xzz_0[i] * fe_0 - ta2_zz_xxx_xzz_1[i] * fe_0 + ta2_zz_xxx_xyzz_0[i] * pa_y[i] - ta2_zz_xxx_xyzz_1[i] * pc_y[i];

        ta2_zz_xxxy_xzzz_0[i] = ta2_zz_xxx_xzzz_0[i] * pa_y[i] - ta2_zz_xxx_xzzz_1[i] * pc_y[i];

        ta2_zz_xxxy_yyyy_0[i] =
            2.0 * ta2_zz_xy_yyyy_0[i] * fe_0 - 2.0 * ta2_zz_xy_yyyy_1[i] * fe_0 + ta2_zz_xxy_yyyy_0[i] * pa_x[i] - ta2_zz_xxy_yyyy_1[i] * pc_x[i];

        ta2_zz_xxxy_yyyz_0[i] =
            2.0 * ta2_zz_xy_yyyz_0[i] * fe_0 - 2.0 * ta2_zz_xy_yyyz_1[i] * fe_0 + ta2_zz_xxy_yyyz_0[i] * pa_x[i] - ta2_zz_xxy_yyyz_1[i] * pc_x[i];

        ta2_zz_xxxy_yyzz_0[i] =
            2.0 * ta2_zz_xy_yyzz_0[i] * fe_0 - 2.0 * ta2_zz_xy_yyzz_1[i] * fe_0 + ta2_zz_xxy_yyzz_0[i] * pa_x[i] - ta2_zz_xxy_yyzz_1[i] * pc_x[i];

        ta2_zz_xxxy_yzzz_0[i] =
            2.0 * ta2_zz_xy_yzzz_0[i] * fe_0 - 2.0 * ta2_zz_xy_yzzz_1[i] * fe_0 + ta2_zz_xxy_yzzz_0[i] * pa_x[i] - ta2_zz_xxy_yzzz_1[i] * pc_x[i];

        ta2_zz_xxxy_zzzz_0[i] = ta2_zz_xxx_zzzz_0[i] * pa_y[i] - ta2_zz_xxx_zzzz_1[i] * pc_y[i];
    }

    // Set up 1155-1170 components of targeted buffer : GG

    auto ta2_zz_xxxz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1155);

    auto ta2_zz_xxxz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1156);

    auto ta2_zz_xxxz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1157);

    auto ta2_zz_xxxz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1158);

    auto ta2_zz_xxxz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1159);

    auto ta2_zz_xxxz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1160);

    auto ta2_zz_xxxz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1161);

    auto ta2_zz_xxxz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1162);

    auto ta2_zz_xxxz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1163);

    auto ta2_zz_xxxz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1164);

    auto ta2_zz_xxxz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1165);

    auto ta2_zz_xxxz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1166);

    auto ta2_zz_xxxz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1167);

    auto ta2_zz_xxxz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1168);

    auto ta2_zz_xxxz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1169);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_z_xxx_xxxx_1,   \
                             ta1_z_xxx_xxxy_1,   \
                             ta1_z_xxx_xxxz_1,   \
                             ta1_z_xxx_xxyy_1,   \
                             ta1_z_xxx_xxyz_1,   \
                             ta1_z_xxx_xxzz_1,   \
                             ta1_z_xxx_xyyy_1,   \
                             ta1_z_xxx_xyyz_1,   \
                             ta1_z_xxx_xyzz_1,   \
                             ta1_z_xxx_xzzz_1,   \
                             ta1_z_xxx_yyyy_1,   \
                             ta2_zz_xxx_xxx_0,   \
                             ta2_zz_xxx_xxx_1,   \
                             ta2_zz_xxx_xxxx_0,  \
                             ta2_zz_xxx_xxxx_1,  \
                             ta2_zz_xxx_xxxy_0,  \
                             ta2_zz_xxx_xxxy_1,  \
                             ta2_zz_xxx_xxxz_0,  \
                             ta2_zz_xxx_xxxz_1,  \
                             ta2_zz_xxx_xxy_0,   \
                             ta2_zz_xxx_xxy_1,   \
                             ta2_zz_xxx_xxyy_0,  \
                             ta2_zz_xxx_xxyy_1,  \
                             ta2_zz_xxx_xxyz_0,  \
                             ta2_zz_xxx_xxyz_1,  \
                             ta2_zz_xxx_xxz_0,   \
                             ta2_zz_xxx_xxz_1,   \
                             ta2_zz_xxx_xxzz_0,  \
                             ta2_zz_xxx_xxzz_1,  \
                             ta2_zz_xxx_xyy_0,   \
                             ta2_zz_xxx_xyy_1,   \
                             ta2_zz_xxx_xyyy_0,  \
                             ta2_zz_xxx_xyyy_1,  \
                             ta2_zz_xxx_xyyz_0,  \
                             ta2_zz_xxx_xyyz_1,  \
                             ta2_zz_xxx_xyz_0,   \
                             ta2_zz_xxx_xyz_1,   \
                             ta2_zz_xxx_xyzz_0,  \
                             ta2_zz_xxx_xyzz_1,  \
                             ta2_zz_xxx_xzz_0,   \
                             ta2_zz_xxx_xzz_1,   \
                             ta2_zz_xxx_xzzz_0,  \
                             ta2_zz_xxx_xzzz_1,  \
                             ta2_zz_xxx_yyyy_0,  \
                             ta2_zz_xxx_yyyy_1,  \
                             ta2_zz_xxxz_xxxx_0, \
                             ta2_zz_xxxz_xxxy_0, \
                             ta2_zz_xxxz_xxxz_0, \
                             ta2_zz_xxxz_xxyy_0, \
                             ta2_zz_xxxz_xxyz_0, \
                             ta2_zz_xxxz_xxzz_0, \
                             ta2_zz_xxxz_xyyy_0, \
                             ta2_zz_xxxz_xyyz_0, \
                             ta2_zz_xxxz_xyzz_0, \
                             ta2_zz_xxxz_xzzz_0, \
                             ta2_zz_xxxz_yyyy_0, \
                             ta2_zz_xxxz_yyyz_0, \
                             ta2_zz_xxxz_yyzz_0, \
                             ta2_zz_xxxz_yzzz_0, \
                             ta2_zz_xxxz_zzzz_0, \
                             ta2_zz_xxz_yyyz_0,  \
                             ta2_zz_xxz_yyyz_1,  \
                             ta2_zz_xxz_yyzz_0,  \
                             ta2_zz_xxz_yyzz_1,  \
                             ta2_zz_xxz_yzzz_0,  \
                             ta2_zz_xxz_yzzz_1,  \
                             ta2_zz_xxz_zzzz_0,  \
                             ta2_zz_xxz_zzzz_1,  \
                             ta2_zz_xz_yyyz_0,   \
                             ta2_zz_xz_yyyz_1,   \
                             ta2_zz_xz_yyzz_0,   \
                             ta2_zz_xz_yyzz_1,   \
                             ta2_zz_xz_yzzz_0,   \
                             ta2_zz_xz_yzzz_1,   \
                             ta2_zz_xz_zzzz_0,   \
                             ta2_zz_xz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxxz_xxxx_0[i] = 2.0 * ta1_z_xxx_xxxx_1[i] + ta2_zz_xxx_xxxx_0[i] * pa_z[i] - ta2_zz_xxx_xxxx_1[i] * pc_z[i];

        ta2_zz_xxxz_xxxy_0[i] = 2.0 * ta1_z_xxx_xxxy_1[i] + ta2_zz_xxx_xxxy_0[i] * pa_z[i] - ta2_zz_xxx_xxxy_1[i] * pc_z[i];

        ta2_zz_xxxz_xxxz_0[i] = ta2_zz_xxx_xxx_0[i] * fe_0 - ta2_zz_xxx_xxx_1[i] * fe_0 + 2.0 * ta1_z_xxx_xxxz_1[i] + ta2_zz_xxx_xxxz_0[i] * pa_z[i] -
                                ta2_zz_xxx_xxxz_1[i] * pc_z[i];

        ta2_zz_xxxz_xxyy_0[i] = 2.0 * ta1_z_xxx_xxyy_1[i] + ta2_zz_xxx_xxyy_0[i] * pa_z[i] - ta2_zz_xxx_xxyy_1[i] * pc_z[i];

        ta2_zz_xxxz_xxyz_0[i] = ta2_zz_xxx_xxy_0[i] * fe_0 - ta2_zz_xxx_xxy_1[i] * fe_0 + 2.0 * ta1_z_xxx_xxyz_1[i] + ta2_zz_xxx_xxyz_0[i] * pa_z[i] -
                                ta2_zz_xxx_xxyz_1[i] * pc_z[i];

        ta2_zz_xxxz_xxzz_0[i] = 2.0 * ta2_zz_xxx_xxz_0[i] * fe_0 - 2.0 * ta2_zz_xxx_xxz_1[i] * fe_0 + 2.0 * ta1_z_xxx_xxzz_1[i] +
                                ta2_zz_xxx_xxzz_0[i] * pa_z[i] - ta2_zz_xxx_xxzz_1[i] * pc_z[i];

        ta2_zz_xxxz_xyyy_0[i] = 2.0 * ta1_z_xxx_xyyy_1[i] + ta2_zz_xxx_xyyy_0[i] * pa_z[i] - ta2_zz_xxx_xyyy_1[i] * pc_z[i];

        ta2_zz_xxxz_xyyz_0[i] = ta2_zz_xxx_xyy_0[i] * fe_0 - ta2_zz_xxx_xyy_1[i] * fe_0 + 2.0 * ta1_z_xxx_xyyz_1[i] + ta2_zz_xxx_xyyz_0[i] * pa_z[i] -
                                ta2_zz_xxx_xyyz_1[i] * pc_z[i];

        ta2_zz_xxxz_xyzz_0[i] = 2.0 * ta2_zz_xxx_xyz_0[i] * fe_0 - 2.0 * ta2_zz_xxx_xyz_1[i] * fe_0 + 2.0 * ta1_z_xxx_xyzz_1[i] +
                                ta2_zz_xxx_xyzz_0[i] * pa_z[i] - ta2_zz_xxx_xyzz_1[i] * pc_z[i];

        ta2_zz_xxxz_xzzz_0[i] = 3.0 * ta2_zz_xxx_xzz_0[i] * fe_0 - 3.0 * ta2_zz_xxx_xzz_1[i] * fe_0 + 2.0 * ta1_z_xxx_xzzz_1[i] +
                                ta2_zz_xxx_xzzz_0[i] * pa_z[i] - ta2_zz_xxx_xzzz_1[i] * pc_z[i];

        ta2_zz_xxxz_yyyy_0[i] = 2.0 * ta1_z_xxx_yyyy_1[i] + ta2_zz_xxx_yyyy_0[i] * pa_z[i] - ta2_zz_xxx_yyyy_1[i] * pc_z[i];

        ta2_zz_xxxz_yyyz_0[i] =
            2.0 * ta2_zz_xz_yyyz_0[i] * fe_0 - 2.0 * ta2_zz_xz_yyyz_1[i] * fe_0 + ta2_zz_xxz_yyyz_0[i] * pa_x[i] - ta2_zz_xxz_yyyz_1[i] * pc_x[i];

        ta2_zz_xxxz_yyzz_0[i] =
            2.0 * ta2_zz_xz_yyzz_0[i] * fe_0 - 2.0 * ta2_zz_xz_yyzz_1[i] * fe_0 + ta2_zz_xxz_yyzz_0[i] * pa_x[i] - ta2_zz_xxz_yyzz_1[i] * pc_x[i];

        ta2_zz_xxxz_yzzz_0[i] =
            2.0 * ta2_zz_xz_yzzz_0[i] * fe_0 - 2.0 * ta2_zz_xz_yzzz_1[i] * fe_0 + ta2_zz_xxz_yzzz_0[i] * pa_x[i] - ta2_zz_xxz_yzzz_1[i] * pc_x[i];

        ta2_zz_xxxz_zzzz_0[i] =
            2.0 * ta2_zz_xz_zzzz_0[i] * fe_0 - 2.0 * ta2_zz_xz_zzzz_1[i] * fe_0 + ta2_zz_xxz_zzzz_0[i] * pa_x[i] - ta2_zz_xxz_zzzz_1[i] * pc_x[i];
    }

    // Set up 1170-1185 components of targeted buffer : GG

    auto ta2_zz_xxyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1170);

    auto ta2_zz_xxyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1171);

    auto ta2_zz_xxyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1172);

    auto ta2_zz_xxyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1173);

    auto ta2_zz_xxyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1174);

    auto ta2_zz_xxyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1175);

    auto ta2_zz_xxyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1176);

    auto ta2_zz_xxyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1177);

    auto ta2_zz_xxyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1178);

    auto ta2_zz_xxyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1179);

    auto ta2_zz_xxyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1180);

    auto ta2_zz_xxyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1181);

    auto ta2_zz_xxyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1182);

    auto ta2_zz_xxyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1183);

    auto ta2_zz_xxyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1184);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta2_zz_xx_xxxx_0,   \
                             ta2_zz_xx_xxxx_1,   \
                             ta2_zz_xx_xxxz_0,   \
                             ta2_zz_xx_xxxz_1,   \
                             ta2_zz_xx_xxzz_0,   \
                             ta2_zz_xx_xxzz_1,   \
                             ta2_zz_xx_xzzz_0,   \
                             ta2_zz_xx_xzzz_1,   \
                             ta2_zz_xxy_xxxx_0,  \
                             ta2_zz_xxy_xxxx_1,  \
                             ta2_zz_xxy_xxxz_0,  \
                             ta2_zz_xxy_xxxz_1,  \
                             ta2_zz_xxy_xxzz_0,  \
                             ta2_zz_xxy_xxzz_1,  \
                             ta2_zz_xxy_xzzz_0,  \
                             ta2_zz_xxy_xzzz_1,  \
                             ta2_zz_xxyy_xxxx_0, \
                             ta2_zz_xxyy_xxxy_0, \
                             ta2_zz_xxyy_xxxz_0, \
                             ta2_zz_xxyy_xxyy_0, \
                             ta2_zz_xxyy_xxyz_0, \
                             ta2_zz_xxyy_xxzz_0, \
                             ta2_zz_xxyy_xyyy_0, \
                             ta2_zz_xxyy_xyyz_0, \
                             ta2_zz_xxyy_xyzz_0, \
                             ta2_zz_xxyy_xzzz_0, \
                             ta2_zz_xxyy_yyyy_0, \
                             ta2_zz_xxyy_yyyz_0, \
                             ta2_zz_xxyy_yyzz_0, \
                             ta2_zz_xxyy_yzzz_0, \
                             ta2_zz_xxyy_zzzz_0, \
                             ta2_zz_xyy_xxxy_0,  \
                             ta2_zz_xyy_xxxy_1,  \
                             ta2_zz_xyy_xxy_0,   \
                             ta2_zz_xyy_xxy_1,   \
                             ta2_zz_xyy_xxyy_0,  \
                             ta2_zz_xyy_xxyy_1,  \
                             ta2_zz_xyy_xxyz_0,  \
                             ta2_zz_xyy_xxyz_1,  \
                             ta2_zz_xyy_xyy_0,   \
                             ta2_zz_xyy_xyy_1,   \
                             ta2_zz_xyy_xyyy_0,  \
                             ta2_zz_xyy_xyyy_1,  \
                             ta2_zz_xyy_xyyz_0,  \
                             ta2_zz_xyy_xyyz_1,  \
                             ta2_zz_xyy_xyz_0,   \
                             ta2_zz_xyy_xyz_1,   \
                             ta2_zz_xyy_xyzz_0,  \
                             ta2_zz_xyy_xyzz_1,  \
                             ta2_zz_xyy_yyy_0,   \
                             ta2_zz_xyy_yyy_1,   \
                             ta2_zz_xyy_yyyy_0,  \
                             ta2_zz_xyy_yyyy_1,  \
                             ta2_zz_xyy_yyyz_0,  \
                             ta2_zz_xyy_yyyz_1,  \
                             ta2_zz_xyy_yyz_0,   \
                             ta2_zz_xyy_yyz_1,   \
                             ta2_zz_xyy_yyzz_0,  \
                             ta2_zz_xyy_yyzz_1,  \
                             ta2_zz_xyy_yzz_0,   \
                             ta2_zz_xyy_yzz_1,   \
                             ta2_zz_xyy_yzzz_0,  \
                             ta2_zz_xyy_yzzz_1,  \
                             ta2_zz_xyy_zzzz_0,  \
                             ta2_zz_xyy_zzzz_1,  \
                             ta2_zz_yy_xxxy_0,   \
                             ta2_zz_yy_xxxy_1,   \
                             ta2_zz_yy_xxyy_0,   \
                             ta2_zz_yy_xxyy_1,   \
                             ta2_zz_yy_xxyz_0,   \
                             ta2_zz_yy_xxyz_1,   \
                             ta2_zz_yy_xyyy_0,   \
                             ta2_zz_yy_xyyy_1,   \
                             ta2_zz_yy_xyyz_0,   \
                             ta2_zz_yy_xyyz_1,   \
                             ta2_zz_yy_xyzz_0,   \
                             ta2_zz_yy_xyzz_1,   \
                             ta2_zz_yy_yyyy_0,   \
                             ta2_zz_yy_yyyy_1,   \
                             ta2_zz_yy_yyyz_0,   \
                             ta2_zz_yy_yyyz_1,   \
                             ta2_zz_yy_yyzz_0,   \
                             ta2_zz_yy_yyzz_1,   \
                             ta2_zz_yy_yzzz_0,   \
                             ta2_zz_yy_yzzz_1,   \
                             ta2_zz_yy_zzzz_0,   \
                             ta2_zz_yy_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxyy_xxxx_0[i] =
            ta2_zz_xx_xxxx_0[i] * fe_0 - ta2_zz_xx_xxxx_1[i] * fe_0 + ta2_zz_xxy_xxxx_0[i] * pa_y[i] - ta2_zz_xxy_xxxx_1[i] * pc_y[i];

        ta2_zz_xxyy_xxxy_0[i] = ta2_zz_yy_xxxy_0[i] * fe_0 - ta2_zz_yy_xxxy_1[i] * fe_0 + 3.0 * ta2_zz_xyy_xxy_0[i] * fe_0 -
                                3.0 * ta2_zz_xyy_xxy_1[i] * fe_0 + ta2_zz_xyy_xxxy_0[i] * pa_x[i] - ta2_zz_xyy_xxxy_1[i] * pc_x[i];

        ta2_zz_xxyy_xxxz_0[i] =
            ta2_zz_xx_xxxz_0[i] * fe_0 - ta2_zz_xx_xxxz_1[i] * fe_0 + ta2_zz_xxy_xxxz_0[i] * pa_y[i] - ta2_zz_xxy_xxxz_1[i] * pc_y[i];

        ta2_zz_xxyy_xxyy_0[i] = ta2_zz_yy_xxyy_0[i] * fe_0 - ta2_zz_yy_xxyy_1[i] * fe_0 + 2.0 * ta2_zz_xyy_xyy_0[i] * fe_0 -
                                2.0 * ta2_zz_xyy_xyy_1[i] * fe_0 + ta2_zz_xyy_xxyy_0[i] * pa_x[i] - ta2_zz_xyy_xxyy_1[i] * pc_x[i];

        ta2_zz_xxyy_xxyz_0[i] = ta2_zz_yy_xxyz_0[i] * fe_0 - ta2_zz_yy_xxyz_1[i] * fe_0 + 2.0 * ta2_zz_xyy_xyz_0[i] * fe_0 -
                                2.0 * ta2_zz_xyy_xyz_1[i] * fe_0 + ta2_zz_xyy_xxyz_0[i] * pa_x[i] - ta2_zz_xyy_xxyz_1[i] * pc_x[i];

        ta2_zz_xxyy_xxzz_0[i] =
            ta2_zz_xx_xxzz_0[i] * fe_0 - ta2_zz_xx_xxzz_1[i] * fe_0 + ta2_zz_xxy_xxzz_0[i] * pa_y[i] - ta2_zz_xxy_xxzz_1[i] * pc_y[i];

        ta2_zz_xxyy_xyyy_0[i] = ta2_zz_yy_xyyy_0[i] * fe_0 - ta2_zz_yy_xyyy_1[i] * fe_0 + ta2_zz_xyy_yyy_0[i] * fe_0 - ta2_zz_xyy_yyy_1[i] * fe_0 +
                                ta2_zz_xyy_xyyy_0[i] * pa_x[i] - ta2_zz_xyy_xyyy_1[i] * pc_x[i];

        ta2_zz_xxyy_xyyz_0[i] = ta2_zz_yy_xyyz_0[i] * fe_0 - ta2_zz_yy_xyyz_1[i] * fe_0 + ta2_zz_xyy_yyz_0[i] * fe_0 - ta2_zz_xyy_yyz_1[i] * fe_0 +
                                ta2_zz_xyy_xyyz_0[i] * pa_x[i] - ta2_zz_xyy_xyyz_1[i] * pc_x[i];

        ta2_zz_xxyy_xyzz_0[i] = ta2_zz_yy_xyzz_0[i] * fe_0 - ta2_zz_yy_xyzz_1[i] * fe_0 + ta2_zz_xyy_yzz_0[i] * fe_0 - ta2_zz_xyy_yzz_1[i] * fe_0 +
                                ta2_zz_xyy_xyzz_0[i] * pa_x[i] - ta2_zz_xyy_xyzz_1[i] * pc_x[i];

        ta2_zz_xxyy_xzzz_0[i] =
            ta2_zz_xx_xzzz_0[i] * fe_0 - ta2_zz_xx_xzzz_1[i] * fe_0 + ta2_zz_xxy_xzzz_0[i] * pa_y[i] - ta2_zz_xxy_xzzz_1[i] * pc_y[i];

        ta2_zz_xxyy_yyyy_0[i] =
            ta2_zz_yy_yyyy_0[i] * fe_0 - ta2_zz_yy_yyyy_1[i] * fe_0 + ta2_zz_xyy_yyyy_0[i] * pa_x[i] - ta2_zz_xyy_yyyy_1[i] * pc_x[i];

        ta2_zz_xxyy_yyyz_0[i] =
            ta2_zz_yy_yyyz_0[i] * fe_0 - ta2_zz_yy_yyyz_1[i] * fe_0 + ta2_zz_xyy_yyyz_0[i] * pa_x[i] - ta2_zz_xyy_yyyz_1[i] * pc_x[i];

        ta2_zz_xxyy_yyzz_0[i] =
            ta2_zz_yy_yyzz_0[i] * fe_0 - ta2_zz_yy_yyzz_1[i] * fe_0 + ta2_zz_xyy_yyzz_0[i] * pa_x[i] - ta2_zz_xyy_yyzz_1[i] * pc_x[i];

        ta2_zz_xxyy_yzzz_0[i] =
            ta2_zz_yy_yzzz_0[i] * fe_0 - ta2_zz_yy_yzzz_1[i] * fe_0 + ta2_zz_xyy_yzzz_0[i] * pa_x[i] - ta2_zz_xyy_yzzz_1[i] * pc_x[i];

        ta2_zz_xxyy_zzzz_0[i] =
            ta2_zz_yy_zzzz_0[i] * fe_0 - ta2_zz_yy_zzzz_1[i] * fe_0 + ta2_zz_xyy_zzzz_0[i] * pa_x[i] - ta2_zz_xyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 1185-1200 components of targeted buffer : GG

    auto ta2_zz_xxyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1185);

    auto ta2_zz_xxyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1186);

    auto ta2_zz_xxyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1187);

    auto ta2_zz_xxyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1188);

    auto ta2_zz_xxyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1189);

    auto ta2_zz_xxyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1190);

    auto ta2_zz_xxyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1191);

    auto ta2_zz_xxyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1192);

    auto ta2_zz_xxyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1193);

    auto ta2_zz_xxyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1194);

    auto ta2_zz_xxyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1195);

    auto ta2_zz_xxyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1196);

    auto ta2_zz_xxyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1197);

    auto ta2_zz_xxyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1198);

    auto ta2_zz_xxyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1199);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_xxy_xxxy_1,   \
                             ta1_z_xxy_xxyy_1,   \
                             ta1_z_xxy_xyyy_1,   \
                             ta1_z_xxy_yyyy_1,   \
                             ta2_zz_xxy_xxxy_0,  \
                             ta2_zz_xxy_xxxy_1,  \
                             ta2_zz_xxy_xxyy_0,  \
                             ta2_zz_xxy_xxyy_1,  \
                             ta2_zz_xxy_xyyy_0,  \
                             ta2_zz_xxy_xyyy_1,  \
                             ta2_zz_xxy_yyyy_0,  \
                             ta2_zz_xxy_yyyy_1,  \
                             ta2_zz_xxyz_xxxx_0, \
                             ta2_zz_xxyz_xxxy_0, \
                             ta2_zz_xxyz_xxxz_0, \
                             ta2_zz_xxyz_xxyy_0, \
                             ta2_zz_xxyz_xxyz_0, \
                             ta2_zz_xxyz_xxzz_0, \
                             ta2_zz_xxyz_xyyy_0, \
                             ta2_zz_xxyz_xyyz_0, \
                             ta2_zz_xxyz_xyzz_0, \
                             ta2_zz_xxyz_xzzz_0, \
                             ta2_zz_xxyz_yyyy_0, \
                             ta2_zz_xxyz_yyyz_0, \
                             ta2_zz_xxyz_yyzz_0, \
                             ta2_zz_xxyz_yzzz_0, \
                             ta2_zz_xxyz_zzzz_0, \
                             ta2_zz_xxz_xxxx_0,  \
                             ta2_zz_xxz_xxxx_1,  \
                             ta2_zz_xxz_xxxz_0,  \
                             ta2_zz_xxz_xxxz_1,  \
                             ta2_zz_xxz_xxyz_0,  \
                             ta2_zz_xxz_xxyz_1,  \
                             ta2_zz_xxz_xxz_0,   \
                             ta2_zz_xxz_xxz_1,   \
                             ta2_zz_xxz_xxzz_0,  \
                             ta2_zz_xxz_xxzz_1,  \
                             ta2_zz_xxz_xyyz_0,  \
                             ta2_zz_xxz_xyyz_1,  \
                             ta2_zz_xxz_xyz_0,   \
                             ta2_zz_xxz_xyz_1,   \
                             ta2_zz_xxz_xyzz_0,  \
                             ta2_zz_xxz_xyzz_1,  \
                             ta2_zz_xxz_xzz_0,   \
                             ta2_zz_xxz_xzz_1,   \
                             ta2_zz_xxz_xzzz_0,  \
                             ta2_zz_xxz_xzzz_1,  \
                             ta2_zz_xxz_zzzz_0,  \
                             ta2_zz_xxz_zzzz_1,  \
                             ta2_zz_xyz_yyyz_0,  \
                             ta2_zz_xyz_yyyz_1,  \
                             ta2_zz_xyz_yyzz_0,  \
                             ta2_zz_xyz_yyzz_1,  \
                             ta2_zz_xyz_yzzz_0,  \
                             ta2_zz_xyz_yzzz_1,  \
                             ta2_zz_yz_yyyz_0,   \
                             ta2_zz_yz_yyyz_1,   \
                             ta2_zz_yz_yyzz_0,   \
                             ta2_zz_yz_yyzz_1,   \
                             ta2_zz_yz_yzzz_0,   \
                             ta2_zz_yz_yzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxyz_xxxx_0[i] = ta2_zz_xxz_xxxx_0[i] * pa_y[i] - ta2_zz_xxz_xxxx_1[i] * pc_y[i];

        ta2_zz_xxyz_xxxy_0[i] = 2.0 * ta1_z_xxy_xxxy_1[i] + ta2_zz_xxy_xxxy_0[i] * pa_z[i] - ta2_zz_xxy_xxxy_1[i] * pc_z[i];

        ta2_zz_xxyz_xxxz_0[i] = ta2_zz_xxz_xxxz_0[i] * pa_y[i] - ta2_zz_xxz_xxxz_1[i] * pc_y[i];

        ta2_zz_xxyz_xxyy_0[i] = 2.0 * ta1_z_xxy_xxyy_1[i] + ta2_zz_xxy_xxyy_0[i] * pa_z[i] - ta2_zz_xxy_xxyy_1[i] * pc_z[i];

        ta2_zz_xxyz_xxyz_0[i] =
            ta2_zz_xxz_xxz_0[i] * fe_0 - ta2_zz_xxz_xxz_1[i] * fe_0 + ta2_zz_xxz_xxyz_0[i] * pa_y[i] - ta2_zz_xxz_xxyz_1[i] * pc_y[i];

        ta2_zz_xxyz_xxzz_0[i] = ta2_zz_xxz_xxzz_0[i] * pa_y[i] - ta2_zz_xxz_xxzz_1[i] * pc_y[i];

        ta2_zz_xxyz_xyyy_0[i] = 2.0 * ta1_z_xxy_xyyy_1[i] + ta2_zz_xxy_xyyy_0[i] * pa_z[i] - ta2_zz_xxy_xyyy_1[i] * pc_z[i];

        ta2_zz_xxyz_xyyz_0[i] =
            2.0 * ta2_zz_xxz_xyz_0[i] * fe_0 - 2.0 * ta2_zz_xxz_xyz_1[i] * fe_0 + ta2_zz_xxz_xyyz_0[i] * pa_y[i] - ta2_zz_xxz_xyyz_1[i] * pc_y[i];

        ta2_zz_xxyz_xyzz_0[i] =
            ta2_zz_xxz_xzz_0[i] * fe_0 - ta2_zz_xxz_xzz_1[i] * fe_0 + ta2_zz_xxz_xyzz_0[i] * pa_y[i] - ta2_zz_xxz_xyzz_1[i] * pc_y[i];

        ta2_zz_xxyz_xzzz_0[i] = ta2_zz_xxz_xzzz_0[i] * pa_y[i] - ta2_zz_xxz_xzzz_1[i] * pc_y[i];

        ta2_zz_xxyz_yyyy_0[i] = 2.0 * ta1_z_xxy_yyyy_1[i] + ta2_zz_xxy_yyyy_0[i] * pa_z[i] - ta2_zz_xxy_yyyy_1[i] * pc_z[i];

        ta2_zz_xxyz_yyyz_0[i] =
            ta2_zz_yz_yyyz_0[i] * fe_0 - ta2_zz_yz_yyyz_1[i] * fe_0 + ta2_zz_xyz_yyyz_0[i] * pa_x[i] - ta2_zz_xyz_yyyz_1[i] * pc_x[i];

        ta2_zz_xxyz_yyzz_0[i] =
            ta2_zz_yz_yyzz_0[i] * fe_0 - ta2_zz_yz_yyzz_1[i] * fe_0 + ta2_zz_xyz_yyzz_0[i] * pa_x[i] - ta2_zz_xyz_yyzz_1[i] * pc_x[i];

        ta2_zz_xxyz_yzzz_0[i] =
            ta2_zz_yz_yzzz_0[i] * fe_0 - ta2_zz_yz_yzzz_1[i] * fe_0 + ta2_zz_xyz_yzzz_0[i] * pa_x[i] - ta2_zz_xyz_yzzz_1[i] * pc_x[i];

        ta2_zz_xxyz_zzzz_0[i] = ta2_zz_xxz_zzzz_0[i] * pa_y[i] - ta2_zz_xxz_zzzz_1[i] * pc_y[i];
    }

    // Set up 1200-1215 components of targeted buffer : GG

    auto ta2_zz_xxzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1200);

    auto ta2_zz_xxzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1201);

    auto ta2_zz_xxzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1202);

    auto ta2_zz_xxzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1203);

    auto ta2_zz_xxzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1204);

    auto ta2_zz_xxzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1205);

    auto ta2_zz_xxzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1206);

    auto ta2_zz_xxzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1207);

    auto ta2_zz_xxzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1208);

    auto ta2_zz_xxzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1209);

    auto ta2_zz_xxzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1210);

    auto ta2_zz_xxzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1211);

    auto ta2_zz_xxzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1212);

    auto ta2_zz_xxzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1213);

    auto ta2_zz_xxzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1214);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_z_xxz_xxxx_1,   \
                             ta1_z_xxz_xxxy_1,   \
                             ta1_z_xxz_xxyy_1,   \
                             ta1_z_xxz_xyyy_1,   \
                             ta2_zz_xx_xxxx_0,   \
                             ta2_zz_xx_xxxx_1,   \
                             ta2_zz_xx_xxxy_0,   \
                             ta2_zz_xx_xxxy_1,   \
                             ta2_zz_xx_xxyy_0,   \
                             ta2_zz_xx_xxyy_1,   \
                             ta2_zz_xx_xyyy_0,   \
                             ta2_zz_xx_xyyy_1,   \
                             ta2_zz_xxz_xxxx_0,  \
                             ta2_zz_xxz_xxxx_1,  \
                             ta2_zz_xxz_xxxy_0,  \
                             ta2_zz_xxz_xxxy_1,  \
                             ta2_zz_xxz_xxyy_0,  \
                             ta2_zz_xxz_xxyy_1,  \
                             ta2_zz_xxz_xyyy_0,  \
                             ta2_zz_xxz_xyyy_1,  \
                             ta2_zz_xxzz_xxxx_0, \
                             ta2_zz_xxzz_xxxy_0, \
                             ta2_zz_xxzz_xxxz_0, \
                             ta2_zz_xxzz_xxyy_0, \
                             ta2_zz_xxzz_xxyz_0, \
                             ta2_zz_xxzz_xxzz_0, \
                             ta2_zz_xxzz_xyyy_0, \
                             ta2_zz_xxzz_xyyz_0, \
                             ta2_zz_xxzz_xyzz_0, \
                             ta2_zz_xxzz_xzzz_0, \
                             ta2_zz_xxzz_yyyy_0, \
                             ta2_zz_xxzz_yyyz_0, \
                             ta2_zz_xxzz_yyzz_0, \
                             ta2_zz_xxzz_yzzz_0, \
                             ta2_zz_xxzz_zzzz_0, \
                             ta2_zz_xzz_xxxz_0,  \
                             ta2_zz_xzz_xxxz_1,  \
                             ta2_zz_xzz_xxyz_0,  \
                             ta2_zz_xzz_xxyz_1,  \
                             ta2_zz_xzz_xxz_0,   \
                             ta2_zz_xzz_xxz_1,   \
                             ta2_zz_xzz_xxzz_0,  \
                             ta2_zz_xzz_xxzz_1,  \
                             ta2_zz_xzz_xyyz_0,  \
                             ta2_zz_xzz_xyyz_1,  \
                             ta2_zz_xzz_xyz_0,   \
                             ta2_zz_xzz_xyz_1,   \
                             ta2_zz_xzz_xyzz_0,  \
                             ta2_zz_xzz_xyzz_1,  \
                             ta2_zz_xzz_xzz_0,   \
                             ta2_zz_xzz_xzz_1,   \
                             ta2_zz_xzz_xzzz_0,  \
                             ta2_zz_xzz_xzzz_1,  \
                             ta2_zz_xzz_yyyy_0,  \
                             ta2_zz_xzz_yyyy_1,  \
                             ta2_zz_xzz_yyyz_0,  \
                             ta2_zz_xzz_yyyz_1,  \
                             ta2_zz_xzz_yyz_0,   \
                             ta2_zz_xzz_yyz_1,   \
                             ta2_zz_xzz_yyzz_0,  \
                             ta2_zz_xzz_yyzz_1,  \
                             ta2_zz_xzz_yzz_0,   \
                             ta2_zz_xzz_yzz_1,   \
                             ta2_zz_xzz_yzzz_0,  \
                             ta2_zz_xzz_yzzz_1,  \
                             ta2_zz_xzz_zzz_0,   \
                             ta2_zz_xzz_zzz_1,   \
                             ta2_zz_xzz_zzzz_0,  \
                             ta2_zz_xzz_zzzz_1,  \
                             ta2_zz_zz_xxxz_0,   \
                             ta2_zz_zz_xxxz_1,   \
                             ta2_zz_zz_xxyz_0,   \
                             ta2_zz_zz_xxyz_1,   \
                             ta2_zz_zz_xxzz_0,   \
                             ta2_zz_zz_xxzz_1,   \
                             ta2_zz_zz_xyyz_0,   \
                             ta2_zz_zz_xyyz_1,   \
                             ta2_zz_zz_xyzz_0,   \
                             ta2_zz_zz_xyzz_1,   \
                             ta2_zz_zz_xzzz_0,   \
                             ta2_zz_zz_xzzz_1,   \
                             ta2_zz_zz_yyyy_0,   \
                             ta2_zz_zz_yyyy_1,   \
                             ta2_zz_zz_yyyz_0,   \
                             ta2_zz_zz_yyyz_1,   \
                             ta2_zz_zz_yyzz_0,   \
                             ta2_zz_zz_yyzz_1,   \
                             ta2_zz_zz_yzzz_0,   \
                             ta2_zz_zz_yzzz_1,   \
                             ta2_zz_zz_zzzz_0,   \
                             ta2_zz_zz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xxzz_xxxx_0[i] = ta2_zz_xx_xxxx_0[i] * fe_0 - ta2_zz_xx_xxxx_1[i] * fe_0 + 2.0 * ta1_z_xxz_xxxx_1[i] + ta2_zz_xxz_xxxx_0[i] * pa_z[i] -
                                ta2_zz_xxz_xxxx_1[i] * pc_z[i];

        ta2_zz_xxzz_xxxy_0[i] = ta2_zz_xx_xxxy_0[i] * fe_0 - ta2_zz_xx_xxxy_1[i] * fe_0 + 2.0 * ta1_z_xxz_xxxy_1[i] + ta2_zz_xxz_xxxy_0[i] * pa_z[i] -
                                ta2_zz_xxz_xxxy_1[i] * pc_z[i];

        ta2_zz_xxzz_xxxz_0[i] = ta2_zz_zz_xxxz_0[i] * fe_0 - ta2_zz_zz_xxxz_1[i] * fe_0 + 3.0 * ta2_zz_xzz_xxz_0[i] * fe_0 -
                                3.0 * ta2_zz_xzz_xxz_1[i] * fe_0 + ta2_zz_xzz_xxxz_0[i] * pa_x[i] - ta2_zz_xzz_xxxz_1[i] * pc_x[i];

        ta2_zz_xxzz_xxyy_0[i] = ta2_zz_xx_xxyy_0[i] * fe_0 - ta2_zz_xx_xxyy_1[i] * fe_0 + 2.0 * ta1_z_xxz_xxyy_1[i] + ta2_zz_xxz_xxyy_0[i] * pa_z[i] -
                                ta2_zz_xxz_xxyy_1[i] * pc_z[i];

        ta2_zz_xxzz_xxyz_0[i] = ta2_zz_zz_xxyz_0[i] * fe_0 - ta2_zz_zz_xxyz_1[i] * fe_0 + 2.0 * ta2_zz_xzz_xyz_0[i] * fe_0 -
                                2.0 * ta2_zz_xzz_xyz_1[i] * fe_0 + ta2_zz_xzz_xxyz_0[i] * pa_x[i] - ta2_zz_xzz_xxyz_1[i] * pc_x[i];

        ta2_zz_xxzz_xxzz_0[i] = ta2_zz_zz_xxzz_0[i] * fe_0 - ta2_zz_zz_xxzz_1[i] * fe_0 + 2.0 * ta2_zz_xzz_xzz_0[i] * fe_0 -
                                2.0 * ta2_zz_xzz_xzz_1[i] * fe_0 + ta2_zz_xzz_xxzz_0[i] * pa_x[i] - ta2_zz_xzz_xxzz_1[i] * pc_x[i];

        ta2_zz_xxzz_xyyy_0[i] = ta2_zz_xx_xyyy_0[i] * fe_0 - ta2_zz_xx_xyyy_1[i] * fe_0 + 2.0 * ta1_z_xxz_xyyy_1[i] + ta2_zz_xxz_xyyy_0[i] * pa_z[i] -
                                ta2_zz_xxz_xyyy_1[i] * pc_z[i];

        ta2_zz_xxzz_xyyz_0[i] = ta2_zz_zz_xyyz_0[i] * fe_0 - ta2_zz_zz_xyyz_1[i] * fe_0 + ta2_zz_xzz_yyz_0[i] * fe_0 - ta2_zz_xzz_yyz_1[i] * fe_0 +
                                ta2_zz_xzz_xyyz_0[i] * pa_x[i] - ta2_zz_xzz_xyyz_1[i] * pc_x[i];

        ta2_zz_xxzz_xyzz_0[i] = ta2_zz_zz_xyzz_0[i] * fe_0 - ta2_zz_zz_xyzz_1[i] * fe_0 + ta2_zz_xzz_yzz_0[i] * fe_0 - ta2_zz_xzz_yzz_1[i] * fe_0 +
                                ta2_zz_xzz_xyzz_0[i] * pa_x[i] - ta2_zz_xzz_xyzz_1[i] * pc_x[i];

        ta2_zz_xxzz_xzzz_0[i] = ta2_zz_zz_xzzz_0[i] * fe_0 - ta2_zz_zz_xzzz_1[i] * fe_0 + ta2_zz_xzz_zzz_0[i] * fe_0 - ta2_zz_xzz_zzz_1[i] * fe_0 +
                                ta2_zz_xzz_xzzz_0[i] * pa_x[i] - ta2_zz_xzz_xzzz_1[i] * pc_x[i];

        ta2_zz_xxzz_yyyy_0[i] =
            ta2_zz_zz_yyyy_0[i] * fe_0 - ta2_zz_zz_yyyy_1[i] * fe_0 + ta2_zz_xzz_yyyy_0[i] * pa_x[i] - ta2_zz_xzz_yyyy_1[i] * pc_x[i];

        ta2_zz_xxzz_yyyz_0[i] =
            ta2_zz_zz_yyyz_0[i] * fe_0 - ta2_zz_zz_yyyz_1[i] * fe_0 + ta2_zz_xzz_yyyz_0[i] * pa_x[i] - ta2_zz_xzz_yyyz_1[i] * pc_x[i];

        ta2_zz_xxzz_yyzz_0[i] =
            ta2_zz_zz_yyzz_0[i] * fe_0 - ta2_zz_zz_yyzz_1[i] * fe_0 + ta2_zz_xzz_yyzz_0[i] * pa_x[i] - ta2_zz_xzz_yyzz_1[i] * pc_x[i];

        ta2_zz_xxzz_yzzz_0[i] =
            ta2_zz_zz_yzzz_0[i] * fe_0 - ta2_zz_zz_yzzz_1[i] * fe_0 + ta2_zz_xzz_yzzz_0[i] * pa_x[i] - ta2_zz_xzz_yzzz_1[i] * pc_x[i];

        ta2_zz_xxzz_zzzz_0[i] =
            ta2_zz_zz_zzzz_0[i] * fe_0 - ta2_zz_zz_zzzz_1[i] * fe_0 + ta2_zz_xzz_zzzz_0[i] * pa_x[i] - ta2_zz_xzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 1215-1230 components of targeted buffer : GG

    auto ta2_zz_xyyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1215);

    auto ta2_zz_xyyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1216);

    auto ta2_zz_xyyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1217);

    auto ta2_zz_xyyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1218);

    auto ta2_zz_xyyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1219);

    auto ta2_zz_xyyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1220);

    auto ta2_zz_xyyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1221);

    auto ta2_zz_xyyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1222);

    auto ta2_zz_xyyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1223);

    auto ta2_zz_xyyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1224);

    auto ta2_zz_xyyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1225);

    auto ta2_zz_xyyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1226);

    auto ta2_zz_xyyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1227);

    auto ta2_zz_xyyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1228);

    auto ta2_zz_xyyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1229);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta2_zz_xyyy_xxxx_0, \
                             ta2_zz_xyyy_xxxy_0, \
                             ta2_zz_xyyy_xxxz_0, \
                             ta2_zz_xyyy_xxyy_0, \
                             ta2_zz_xyyy_xxyz_0, \
                             ta2_zz_xyyy_xxzz_0, \
                             ta2_zz_xyyy_xyyy_0, \
                             ta2_zz_xyyy_xyyz_0, \
                             ta2_zz_xyyy_xyzz_0, \
                             ta2_zz_xyyy_xzzz_0, \
                             ta2_zz_xyyy_yyyy_0, \
                             ta2_zz_xyyy_yyyz_0, \
                             ta2_zz_xyyy_yyzz_0, \
                             ta2_zz_xyyy_yzzz_0, \
                             ta2_zz_xyyy_zzzz_0, \
                             ta2_zz_yyy_xxx_0,   \
                             ta2_zz_yyy_xxx_1,   \
                             ta2_zz_yyy_xxxx_0,  \
                             ta2_zz_yyy_xxxx_1,  \
                             ta2_zz_yyy_xxxy_0,  \
                             ta2_zz_yyy_xxxy_1,  \
                             ta2_zz_yyy_xxxz_0,  \
                             ta2_zz_yyy_xxxz_1,  \
                             ta2_zz_yyy_xxy_0,   \
                             ta2_zz_yyy_xxy_1,   \
                             ta2_zz_yyy_xxyy_0,  \
                             ta2_zz_yyy_xxyy_1,  \
                             ta2_zz_yyy_xxyz_0,  \
                             ta2_zz_yyy_xxyz_1,  \
                             ta2_zz_yyy_xxz_0,   \
                             ta2_zz_yyy_xxz_1,   \
                             ta2_zz_yyy_xxzz_0,  \
                             ta2_zz_yyy_xxzz_1,  \
                             ta2_zz_yyy_xyy_0,   \
                             ta2_zz_yyy_xyy_1,   \
                             ta2_zz_yyy_xyyy_0,  \
                             ta2_zz_yyy_xyyy_1,  \
                             ta2_zz_yyy_xyyz_0,  \
                             ta2_zz_yyy_xyyz_1,  \
                             ta2_zz_yyy_xyz_0,   \
                             ta2_zz_yyy_xyz_1,   \
                             ta2_zz_yyy_xyzz_0,  \
                             ta2_zz_yyy_xyzz_1,  \
                             ta2_zz_yyy_xzz_0,   \
                             ta2_zz_yyy_xzz_1,   \
                             ta2_zz_yyy_xzzz_0,  \
                             ta2_zz_yyy_xzzz_1,  \
                             ta2_zz_yyy_yyy_0,   \
                             ta2_zz_yyy_yyy_1,   \
                             ta2_zz_yyy_yyyy_0,  \
                             ta2_zz_yyy_yyyy_1,  \
                             ta2_zz_yyy_yyyz_0,  \
                             ta2_zz_yyy_yyyz_1,  \
                             ta2_zz_yyy_yyz_0,   \
                             ta2_zz_yyy_yyz_1,   \
                             ta2_zz_yyy_yyzz_0,  \
                             ta2_zz_yyy_yyzz_1,  \
                             ta2_zz_yyy_yzz_0,   \
                             ta2_zz_yyy_yzz_1,   \
                             ta2_zz_yyy_yzzz_0,  \
                             ta2_zz_yyy_yzzz_1,  \
                             ta2_zz_yyy_zzz_0,   \
                             ta2_zz_yyy_zzz_1,   \
                             ta2_zz_yyy_zzzz_0,  \
                             ta2_zz_yyy_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xyyy_xxxx_0[i] =
            4.0 * ta2_zz_yyy_xxx_0[i] * fe_0 - 4.0 * ta2_zz_yyy_xxx_1[i] * fe_0 + ta2_zz_yyy_xxxx_0[i] * pa_x[i] - ta2_zz_yyy_xxxx_1[i] * pc_x[i];

        ta2_zz_xyyy_xxxy_0[i] =
            3.0 * ta2_zz_yyy_xxy_0[i] * fe_0 - 3.0 * ta2_zz_yyy_xxy_1[i] * fe_0 + ta2_zz_yyy_xxxy_0[i] * pa_x[i] - ta2_zz_yyy_xxxy_1[i] * pc_x[i];

        ta2_zz_xyyy_xxxz_0[i] =
            3.0 * ta2_zz_yyy_xxz_0[i] * fe_0 - 3.0 * ta2_zz_yyy_xxz_1[i] * fe_0 + ta2_zz_yyy_xxxz_0[i] * pa_x[i] - ta2_zz_yyy_xxxz_1[i] * pc_x[i];

        ta2_zz_xyyy_xxyy_0[i] =
            2.0 * ta2_zz_yyy_xyy_0[i] * fe_0 - 2.0 * ta2_zz_yyy_xyy_1[i] * fe_0 + ta2_zz_yyy_xxyy_0[i] * pa_x[i] - ta2_zz_yyy_xxyy_1[i] * pc_x[i];

        ta2_zz_xyyy_xxyz_0[i] =
            2.0 * ta2_zz_yyy_xyz_0[i] * fe_0 - 2.0 * ta2_zz_yyy_xyz_1[i] * fe_0 + ta2_zz_yyy_xxyz_0[i] * pa_x[i] - ta2_zz_yyy_xxyz_1[i] * pc_x[i];

        ta2_zz_xyyy_xxzz_0[i] =
            2.0 * ta2_zz_yyy_xzz_0[i] * fe_0 - 2.0 * ta2_zz_yyy_xzz_1[i] * fe_0 + ta2_zz_yyy_xxzz_0[i] * pa_x[i] - ta2_zz_yyy_xxzz_1[i] * pc_x[i];

        ta2_zz_xyyy_xyyy_0[i] =
            ta2_zz_yyy_yyy_0[i] * fe_0 - ta2_zz_yyy_yyy_1[i] * fe_0 + ta2_zz_yyy_xyyy_0[i] * pa_x[i] - ta2_zz_yyy_xyyy_1[i] * pc_x[i];

        ta2_zz_xyyy_xyyz_0[i] =
            ta2_zz_yyy_yyz_0[i] * fe_0 - ta2_zz_yyy_yyz_1[i] * fe_0 + ta2_zz_yyy_xyyz_0[i] * pa_x[i] - ta2_zz_yyy_xyyz_1[i] * pc_x[i];

        ta2_zz_xyyy_xyzz_0[i] =
            ta2_zz_yyy_yzz_0[i] * fe_0 - ta2_zz_yyy_yzz_1[i] * fe_0 + ta2_zz_yyy_xyzz_0[i] * pa_x[i] - ta2_zz_yyy_xyzz_1[i] * pc_x[i];

        ta2_zz_xyyy_xzzz_0[i] =
            ta2_zz_yyy_zzz_0[i] * fe_0 - ta2_zz_yyy_zzz_1[i] * fe_0 + ta2_zz_yyy_xzzz_0[i] * pa_x[i] - ta2_zz_yyy_xzzz_1[i] * pc_x[i];

        ta2_zz_xyyy_yyyy_0[i] = ta2_zz_yyy_yyyy_0[i] * pa_x[i] - ta2_zz_yyy_yyyy_1[i] * pc_x[i];

        ta2_zz_xyyy_yyyz_0[i] = ta2_zz_yyy_yyyz_0[i] * pa_x[i] - ta2_zz_yyy_yyyz_1[i] * pc_x[i];

        ta2_zz_xyyy_yyzz_0[i] = ta2_zz_yyy_yyzz_0[i] * pa_x[i] - ta2_zz_yyy_yyzz_1[i] * pc_x[i];

        ta2_zz_xyyy_yzzz_0[i] = ta2_zz_yyy_yzzz_0[i] * pa_x[i] - ta2_zz_yyy_yzzz_1[i] * pc_x[i];

        ta2_zz_xyyy_zzzz_0[i] = ta2_zz_yyy_zzzz_0[i] * pa_x[i] - ta2_zz_yyy_zzzz_1[i] * pc_x[i];
    }

    // Set up 1230-1245 components of targeted buffer : GG

    auto ta2_zz_xyyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1230);

    auto ta2_zz_xyyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1231);

    auto ta2_zz_xyyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1232);

    auto ta2_zz_xyyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1233);

    auto ta2_zz_xyyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1234);

    auto ta2_zz_xyyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1235);

    auto ta2_zz_xyyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1236);

    auto ta2_zz_xyyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1237);

    auto ta2_zz_xyyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1238);

    auto ta2_zz_xyyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1239);

    auto ta2_zz_xyyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1240);

    auto ta2_zz_xyyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1241);

    auto ta2_zz_xyyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1242);

    auto ta2_zz_xyyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1243);

    auto ta2_zz_xyyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1244);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_z_xyy_xxxx_1,   \
                             ta1_z_xyy_xxxy_1,   \
                             ta1_z_xyy_xxyy_1,   \
                             ta1_z_xyy_xyyy_1,   \
                             ta2_zz_xyy_xxxx_0,  \
                             ta2_zz_xyy_xxxx_1,  \
                             ta2_zz_xyy_xxxy_0,  \
                             ta2_zz_xyy_xxxy_1,  \
                             ta2_zz_xyy_xxyy_0,  \
                             ta2_zz_xyy_xxyy_1,  \
                             ta2_zz_xyy_xyyy_0,  \
                             ta2_zz_xyy_xyyy_1,  \
                             ta2_zz_xyyz_xxxx_0, \
                             ta2_zz_xyyz_xxxy_0, \
                             ta2_zz_xyyz_xxxz_0, \
                             ta2_zz_xyyz_xxyy_0, \
                             ta2_zz_xyyz_xxyz_0, \
                             ta2_zz_xyyz_xxzz_0, \
                             ta2_zz_xyyz_xyyy_0, \
                             ta2_zz_xyyz_xyyz_0, \
                             ta2_zz_xyyz_xyzz_0, \
                             ta2_zz_xyyz_xzzz_0, \
                             ta2_zz_xyyz_yyyy_0, \
                             ta2_zz_xyyz_yyyz_0, \
                             ta2_zz_xyyz_yyzz_0, \
                             ta2_zz_xyyz_yzzz_0, \
                             ta2_zz_xyyz_zzzz_0, \
                             ta2_zz_yyz_xxxz_0,  \
                             ta2_zz_yyz_xxxz_1,  \
                             ta2_zz_yyz_xxyz_0,  \
                             ta2_zz_yyz_xxyz_1,  \
                             ta2_zz_yyz_xxz_0,   \
                             ta2_zz_yyz_xxz_1,   \
                             ta2_zz_yyz_xxzz_0,  \
                             ta2_zz_yyz_xxzz_1,  \
                             ta2_zz_yyz_xyyz_0,  \
                             ta2_zz_yyz_xyyz_1,  \
                             ta2_zz_yyz_xyz_0,   \
                             ta2_zz_yyz_xyz_1,   \
                             ta2_zz_yyz_xyzz_0,  \
                             ta2_zz_yyz_xyzz_1,  \
                             ta2_zz_yyz_xzz_0,   \
                             ta2_zz_yyz_xzz_1,   \
                             ta2_zz_yyz_xzzz_0,  \
                             ta2_zz_yyz_xzzz_1,  \
                             ta2_zz_yyz_yyyy_0,  \
                             ta2_zz_yyz_yyyy_1,  \
                             ta2_zz_yyz_yyyz_0,  \
                             ta2_zz_yyz_yyyz_1,  \
                             ta2_zz_yyz_yyz_0,   \
                             ta2_zz_yyz_yyz_1,   \
                             ta2_zz_yyz_yyzz_0,  \
                             ta2_zz_yyz_yyzz_1,  \
                             ta2_zz_yyz_yzz_0,   \
                             ta2_zz_yyz_yzz_1,   \
                             ta2_zz_yyz_yzzz_0,  \
                             ta2_zz_yyz_yzzz_1,  \
                             ta2_zz_yyz_zzz_0,   \
                             ta2_zz_yyz_zzz_1,   \
                             ta2_zz_yyz_zzzz_0,  \
                             ta2_zz_yyz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xyyz_xxxx_0[i] = 2.0 * ta1_z_xyy_xxxx_1[i] + ta2_zz_xyy_xxxx_0[i] * pa_z[i] - ta2_zz_xyy_xxxx_1[i] * pc_z[i];

        ta2_zz_xyyz_xxxy_0[i] = 2.0 * ta1_z_xyy_xxxy_1[i] + ta2_zz_xyy_xxxy_0[i] * pa_z[i] - ta2_zz_xyy_xxxy_1[i] * pc_z[i];

        ta2_zz_xyyz_xxxz_0[i] =
            3.0 * ta2_zz_yyz_xxz_0[i] * fe_0 - 3.0 * ta2_zz_yyz_xxz_1[i] * fe_0 + ta2_zz_yyz_xxxz_0[i] * pa_x[i] - ta2_zz_yyz_xxxz_1[i] * pc_x[i];

        ta2_zz_xyyz_xxyy_0[i] = 2.0 * ta1_z_xyy_xxyy_1[i] + ta2_zz_xyy_xxyy_0[i] * pa_z[i] - ta2_zz_xyy_xxyy_1[i] * pc_z[i];

        ta2_zz_xyyz_xxyz_0[i] =
            2.0 * ta2_zz_yyz_xyz_0[i] * fe_0 - 2.0 * ta2_zz_yyz_xyz_1[i] * fe_0 + ta2_zz_yyz_xxyz_0[i] * pa_x[i] - ta2_zz_yyz_xxyz_1[i] * pc_x[i];

        ta2_zz_xyyz_xxzz_0[i] =
            2.0 * ta2_zz_yyz_xzz_0[i] * fe_0 - 2.0 * ta2_zz_yyz_xzz_1[i] * fe_0 + ta2_zz_yyz_xxzz_0[i] * pa_x[i] - ta2_zz_yyz_xxzz_1[i] * pc_x[i];

        ta2_zz_xyyz_xyyy_0[i] = 2.0 * ta1_z_xyy_xyyy_1[i] + ta2_zz_xyy_xyyy_0[i] * pa_z[i] - ta2_zz_xyy_xyyy_1[i] * pc_z[i];

        ta2_zz_xyyz_xyyz_0[i] =
            ta2_zz_yyz_yyz_0[i] * fe_0 - ta2_zz_yyz_yyz_1[i] * fe_0 + ta2_zz_yyz_xyyz_0[i] * pa_x[i] - ta2_zz_yyz_xyyz_1[i] * pc_x[i];

        ta2_zz_xyyz_xyzz_0[i] =
            ta2_zz_yyz_yzz_0[i] * fe_0 - ta2_zz_yyz_yzz_1[i] * fe_0 + ta2_zz_yyz_xyzz_0[i] * pa_x[i] - ta2_zz_yyz_xyzz_1[i] * pc_x[i];

        ta2_zz_xyyz_xzzz_0[i] =
            ta2_zz_yyz_zzz_0[i] * fe_0 - ta2_zz_yyz_zzz_1[i] * fe_0 + ta2_zz_yyz_xzzz_0[i] * pa_x[i] - ta2_zz_yyz_xzzz_1[i] * pc_x[i];

        ta2_zz_xyyz_yyyy_0[i] = ta2_zz_yyz_yyyy_0[i] * pa_x[i] - ta2_zz_yyz_yyyy_1[i] * pc_x[i];

        ta2_zz_xyyz_yyyz_0[i] = ta2_zz_yyz_yyyz_0[i] * pa_x[i] - ta2_zz_yyz_yyyz_1[i] * pc_x[i];

        ta2_zz_xyyz_yyzz_0[i] = ta2_zz_yyz_yyzz_0[i] * pa_x[i] - ta2_zz_yyz_yyzz_1[i] * pc_x[i];

        ta2_zz_xyyz_yzzz_0[i] = ta2_zz_yyz_yzzz_0[i] * pa_x[i] - ta2_zz_yyz_yzzz_1[i] * pc_x[i];

        ta2_zz_xyyz_zzzz_0[i] = ta2_zz_yyz_zzzz_0[i] * pa_x[i] - ta2_zz_yyz_zzzz_1[i] * pc_x[i];
    }

    // Set up 1245-1260 components of targeted buffer : GG

    auto ta2_zz_xyzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1245);

    auto ta2_zz_xyzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1246);

    auto ta2_zz_xyzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1247);

    auto ta2_zz_xyzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1248);

    auto ta2_zz_xyzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1249);

    auto ta2_zz_xyzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1250);

    auto ta2_zz_xyzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1251);

    auto ta2_zz_xyzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1252);

    auto ta2_zz_xyzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1253);

    auto ta2_zz_xyzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1254);

    auto ta2_zz_xyzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1255);

    auto ta2_zz_xyzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1256);

    auto ta2_zz_xyzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1257);

    auto ta2_zz_xyzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1258);

    auto ta2_zz_xyzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1259);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta2_zz_xyzz_xxxx_0, \
                             ta2_zz_xyzz_xxxy_0, \
                             ta2_zz_xyzz_xxxz_0, \
                             ta2_zz_xyzz_xxyy_0, \
                             ta2_zz_xyzz_xxyz_0, \
                             ta2_zz_xyzz_xxzz_0, \
                             ta2_zz_xyzz_xyyy_0, \
                             ta2_zz_xyzz_xyyz_0, \
                             ta2_zz_xyzz_xyzz_0, \
                             ta2_zz_xyzz_xzzz_0, \
                             ta2_zz_xyzz_yyyy_0, \
                             ta2_zz_xyzz_yyyz_0, \
                             ta2_zz_xyzz_yyzz_0, \
                             ta2_zz_xyzz_yzzz_0, \
                             ta2_zz_xyzz_zzzz_0, \
                             ta2_zz_xzz_xxxx_0,  \
                             ta2_zz_xzz_xxxx_1,  \
                             ta2_zz_xzz_xxxz_0,  \
                             ta2_zz_xzz_xxxz_1,  \
                             ta2_zz_xzz_xxzz_0,  \
                             ta2_zz_xzz_xxzz_1,  \
                             ta2_zz_xzz_xzzz_0,  \
                             ta2_zz_xzz_xzzz_1,  \
                             ta2_zz_yzz_xxxy_0,  \
                             ta2_zz_yzz_xxxy_1,  \
                             ta2_zz_yzz_xxy_0,   \
                             ta2_zz_yzz_xxy_1,   \
                             ta2_zz_yzz_xxyy_0,  \
                             ta2_zz_yzz_xxyy_1,  \
                             ta2_zz_yzz_xxyz_0,  \
                             ta2_zz_yzz_xxyz_1,  \
                             ta2_zz_yzz_xyy_0,   \
                             ta2_zz_yzz_xyy_1,   \
                             ta2_zz_yzz_xyyy_0,  \
                             ta2_zz_yzz_xyyy_1,  \
                             ta2_zz_yzz_xyyz_0,  \
                             ta2_zz_yzz_xyyz_1,  \
                             ta2_zz_yzz_xyz_0,   \
                             ta2_zz_yzz_xyz_1,   \
                             ta2_zz_yzz_xyzz_0,  \
                             ta2_zz_yzz_xyzz_1,  \
                             ta2_zz_yzz_yyy_0,   \
                             ta2_zz_yzz_yyy_1,   \
                             ta2_zz_yzz_yyyy_0,  \
                             ta2_zz_yzz_yyyy_1,  \
                             ta2_zz_yzz_yyyz_0,  \
                             ta2_zz_yzz_yyyz_1,  \
                             ta2_zz_yzz_yyz_0,   \
                             ta2_zz_yzz_yyz_1,   \
                             ta2_zz_yzz_yyzz_0,  \
                             ta2_zz_yzz_yyzz_1,  \
                             ta2_zz_yzz_yzz_0,   \
                             ta2_zz_yzz_yzz_1,   \
                             ta2_zz_yzz_yzzz_0,  \
                             ta2_zz_yzz_yzzz_1,  \
                             ta2_zz_yzz_zzzz_0,  \
                             ta2_zz_yzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xyzz_xxxx_0[i] = ta2_zz_xzz_xxxx_0[i] * pa_y[i] - ta2_zz_xzz_xxxx_1[i] * pc_y[i];

        ta2_zz_xyzz_xxxy_0[i] =
            3.0 * ta2_zz_yzz_xxy_0[i] * fe_0 - 3.0 * ta2_zz_yzz_xxy_1[i] * fe_0 + ta2_zz_yzz_xxxy_0[i] * pa_x[i] - ta2_zz_yzz_xxxy_1[i] * pc_x[i];

        ta2_zz_xyzz_xxxz_0[i] = ta2_zz_xzz_xxxz_0[i] * pa_y[i] - ta2_zz_xzz_xxxz_1[i] * pc_y[i];

        ta2_zz_xyzz_xxyy_0[i] =
            2.0 * ta2_zz_yzz_xyy_0[i] * fe_0 - 2.0 * ta2_zz_yzz_xyy_1[i] * fe_0 + ta2_zz_yzz_xxyy_0[i] * pa_x[i] - ta2_zz_yzz_xxyy_1[i] * pc_x[i];

        ta2_zz_xyzz_xxyz_0[i] =
            2.0 * ta2_zz_yzz_xyz_0[i] * fe_0 - 2.0 * ta2_zz_yzz_xyz_1[i] * fe_0 + ta2_zz_yzz_xxyz_0[i] * pa_x[i] - ta2_zz_yzz_xxyz_1[i] * pc_x[i];

        ta2_zz_xyzz_xxzz_0[i] = ta2_zz_xzz_xxzz_0[i] * pa_y[i] - ta2_zz_xzz_xxzz_1[i] * pc_y[i];

        ta2_zz_xyzz_xyyy_0[i] =
            ta2_zz_yzz_yyy_0[i] * fe_0 - ta2_zz_yzz_yyy_1[i] * fe_0 + ta2_zz_yzz_xyyy_0[i] * pa_x[i] - ta2_zz_yzz_xyyy_1[i] * pc_x[i];

        ta2_zz_xyzz_xyyz_0[i] =
            ta2_zz_yzz_yyz_0[i] * fe_0 - ta2_zz_yzz_yyz_1[i] * fe_0 + ta2_zz_yzz_xyyz_0[i] * pa_x[i] - ta2_zz_yzz_xyyz_1[i] * pc_x[i];

        ta2_zz_xyzz_xyzz_0[i] =
            ta2_zz_yzz_yzz_0[i] * fe_0 - ta2_zz_yzz_yzz_1[i] * fe_0 + ta2_zz_yzz_xyzz_0[i] * pa_x[i] - ta2_zz_yzz_xyzz_1[i] * pc_x[i];

        ta2_zz_xyzz_xzzz_0[i] = ta2_zz_xzz_xzzz_0[i] * pa_y[i] - ta2_zz_xzz_xzzz_1[i] * pc_y[i];

        ta2_zz_xyzz_yyyy_0[i] = ta2_zz_yzz_yyyy_0[i] * pa_x[i] - ta2_zz_yzz_yyyy_1[i] * pc_x[i];

        ta2_zz_xyzz_yyyz_0[i] = ta2_zz_yzz_yyyz_0[i] * pa_x[i] - ta2_zz_yzz_yyyz_1[i] * pc_x[i];

        ta2_zz_xyzz_yyzz_0[i] = ta2_zz_yzz_yyzz_0[i] * pa_x[i] - ta2_zz_yzz_yyzz_1[i] * pc_x[i];

        ta2_zz_xyzz_yzzz_0[i] = ta2_zz_yzz_yzzz_0[i] * pa_x[i] - ta2_zz_yzz_yzzz_1[i] * pc_x[i];

        ta2_zz_xyzz_zzzz_0[i] = ta2_zz_yzz_zzzz_0[i] * pa_x[i] - ta2_zz_yzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 1260-1275 components of targeted buffer : GG

    auto ta2_zz_xzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1260);

    auto ta2_zz_xzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1261);

    auto ta2_zz_xzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1262);

    auto ta2_zz_xzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1263);

    auto ta2_zz_xzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1264);

    auto ta2_zz_xzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1265);

    auto ta2_zz_xzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1266);

    auto ta2_zz_xzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1267);

    auto ta2_zz_xzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1268);

    auto ta2_zz_xzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1269);

    auto ta2_zz_xzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1270);

    auto ta2_zz_xzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1271);

    auto ta2_zz_xzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1272);

    auto ta2_zz_xzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1273);

    auto ta2_zz_xzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1274);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta2_zz_xzzz_xxxx_0, \
                             ta2_zz_xzzz_xxxy_0, \
                             ta2_zz_xzzz_xxxz_0, \
                             ta2_zz_xzzz_xxyy_0, \
                             ta2_zz_xzzz_xxyz_0, \
                             ta2_zz_xzzz_xxzz_0, \
                             ta2_zz_xzzz_xyyy_0, \
                             ta2_zz_xzzz_xyyz_0, \
                             ta2_zz_xzzz_xyzz_0, \
                             ta2_zz_xzzz_xzzz_0, \
                             ta2_zz_xzzz_yyyy_0, \
                             ta2_zz_xzzz_yyyz_0, \
                             ta2_zz_xzzz_yyzz_0, \
                             ta2_zz_xzzz_yzzz_0, \
                             ta2_zz_xzzz_zzzz_0, \
                             ta2_zz_zzz_xxx_0,   \
                             ta2_zz_zzz_xxx_1,   \
                             ta2_zz_zzz_xxxx_0,  \
                             ta2_zz_zzz_xxxx_1,  \
                             ta2_zz_zzz_xxxy_0,  \
                             ta2_zz_zzz_xxxy_1,  \
                             ta2_zz_zzz_xxxz_0,  \
                             ta2_zz_zzz_xxxz_1,  \
                             ta2_zz_zzz_xxy_0,   \
                             ta2_zz_zzz_xxy_1,   \
                             ta2_zz_zzz_xxyy_0,  \
                             ta2_zz_zzz_xxyy_1,  \
                             ta2_zz_zzz_xxyz_0,  \
                             ta2_zz_zzz_xxyz_1,  \
                             ta2_zz_zzz_xxz_0,   \
                             ta2_zz_zzz_xxz_1,   \
                             ta2_zz_zzz_xxzz_0,  \
                             ta2_zz_zzz_xxzz_1,  \
                             ta2_zz_zzz_xyy_0,   \
                             ta2_zz_zzz_xyy_1,   \
                             ta2_zz_zzz_xyyy_0,  \
                             ta2_zz_zzz_xyyy_1,  \
                             ta2_zz_zzz_xyyz_0,  \
                             ta2_zz_zzz_xyyz_1,  \
                             ta2_zz_zzz_xyz_0,   \
                             ta2_zz_zzz_xyz_1,   \
                             ta2_zz_zzz_xyzz_0,  \
                             ta2_zz_zzz_xyzz_1,  \
                             ta2_zz_zzz_xzz_0,   \
                             ta2_zz_zzz_xzz_1,   \
                             ta2_zz_zzz_xzzz_0,  \
                             ta2_zz_zzz_xzzz_1,  \
                             ta2_zz_zzz_yyy_0,   \
                             ta2_zz_zzz_yyy_1,   \
                             ta2_zz_zzz_yyyy_0,  \
                             ta2_zz_zzz_yyyy_1,  \
                             ta2_zz_zzz_yyyz_0,  \
                             ta2_zz_zzz_yyyz_1,  \
                             ta2_zz_zzz_yyz_0,   \
                             ta2_zz_zzz_yyz_1,   \
                             ta2_zz_zzz_yyzz_0,  \
                             ta2_zz_zzz_yyzz_1,  \
                             ta2_zz_zzz_yzz_0,   \
                             ta2_zz_zzz_yzz_1,   \
                             ta2_zz_zzz_yzzz_0,  \
                             ta2_zz_zzz_yzzz_1,  \
                             ta2_zz_zzz_zzz_0,   \
                             ta2_zz_zzz_zzz_1,   \
                             ta2_zz_zzz_zzzz_0,  \
                             ta2_zz_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_xzzz_xxxx_0[i] =
            4.0 * ta2_zz_zzz_xxx_0[i] * fe_0 - 4.0 * ta2_zz_zzz_xxx_1[i] * fe_0 + ta2_zz_zzz_xxxx_0[i] * pa_x[i] - ta2_zz_zzz_xxxx_1[i] * pc_x[i];

        ta2_zz_xzzz_xxxy_0[i] =
            3.0 * ta2_zz_zzz_xxy_0[i] * fe_0 - 3.0 * ta2_zz_zzz_xxy_1[i] * fe_0 + ta2_zz_zzz_xxxy_0[i] * pa_x[i] - ta2_zz_zzz_xxxy_1[i] * pc_x[i];

        ta2_zz_xzzz_xxxz_0[i] =
            3.0 * ta2_zz_zzz_xxz_0[i] * fe_0 - 3.0 * ta2_zz_zzz_xxz_1[i] * fe_0 + ta2_zz_zzz_xxxz_0[i] * pa_x[i] - ta2_zz_zzz_xxxz_1[i] * pc_x[i];

        ta2_zz_xzzz_xxyy_0[i] =
            2.0 * ta2_zz_zzz_xyy_0[i] * fe_0 - 2.0 * ta2_zz_zzz_xyy_1[i] * fe_0 + ta2_zz_zzz_xxyy_0[i] * pa_x[i] - ta2_zz_zzz_xxyy_1[i] * pc_x[i];

        ta2_zz_xzzz_xxyz_0[i] =
            2.0 * ta2_zz_zzz_xyz_0[i] * fe_0 - 2.0 * ta2_zz_zzz_xyz_1[i] * fe_0 + ta2_zz_zzz_xxyz_0[i] * pa_x[i] - ta2_zz_zzz_xxyz_1[i] * pc_x[i];

        ta2_zz_xzzz_xxzz_0[i] =
            2.0 * ta2_zz_zzz_xzz_0[i] * fe_0 - 2.0 * ta2_zz_zzz_xzz_1[i] * fe_0 + ta2_zz_zzz_xxzz_0[i] * pa_x[i] - ta2_zz_zzz_xxzz_1[i] * pc_x[i];

        ta2_zz_xzzz_xyyy_0[i] =
            ta2_zz_zzz_yyy_0[i] * fe_0 - ta2_zz_zzz_yyy_1[i] * fe_0 + ta2_zz_zzz_xyyy_0[i] * pa_x[i] - ta2_zz_zzz_xyyy_1[i] * pc_x[i];

        ta2_zz_xzzz_xyyz_0[i] =
            ta2_zz_zzz_yyz_0[i] * fe_0 - ta2_zz_zzz_yyz_1[i] * fe_0 + ta2_zz_zzz_xyyz_0[i] * pa_x[i] - ta2_zz_zzz_xyyz_1[i] * pc_x[i];

        ta2_zz_xzzz_xyzz_0[i] =
            ta2_zz_zzz_yzz_0[i] * fe_0 - ta2_zz_zzz_yzz_1[i] * fe_0 + ta2_zz_zzz_xyzz_0[i] * pa_x[i] - ta2_zz_zzz_xyzz_1[i] * pc_x[i];

        ta2_zz_xzzz_xzzz_0[i] =
            ta2_zz_zzz_zzz_0[i] * fe_0 - ta2_zz_zzz_zzz_1[i] * fe_0 + ta2_zz_zzz_xzzz_0[i] * pa_x[i] - ta2_zz_zzz_xzzz_1[i] * pc_x[i];

        ta2_zz_xzzz_yyyy_0[i] = ta2_zz_zzz_yyyy_0[i] * pa_x[i] - ta2_zz_zzz_yyyy_1[i] * pc_x[i];

        ta2_zz_xzzz_yyyz_0[i] = ta2_zz_zzz_yyyz_0[i] * pa_x[i] - ta2_zz_zzz_yyyz_1[i] * pc_x[i];

        ta2_zz_xzzz_yyzz_0[i] = ta2_zz_zzz_yyzz_0[i] * pa_x[i] - ta2_zz_zzz_yyzz_1[i] * pc_x[i];

        ta2_zz_xzzz_yzzz_0[i] = ta2_zz_zzz_yzzz_0[i] * pa_x[i] - ta2_zz_zzz_yzzz_1[i] * pc_x[i];

        ta2_zz_xzzz_zzzz_0[i] = ta2_zz_zzz_zzzz_0[i] * pa_x[i] - ta2_zz_zzz_zzzz_1[i] * pc_x[i];
    }

    // Set up 1275-1290 components of targeted buffer : GG

    auto ta2_zz_yyyy_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1275);

    auto ta2_zz_yyyy_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1276);

    auto ta2_zz_yyyy_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1277);

    auto ta2_zz_yyyy_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1278);

    auto ta2_zz_yyyy_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1279);

    auto ta2_zz_yyyy_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1280);

    auto ta2_zz_yyyy_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1281);

    auto ta2_zz_yyyy_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1282);

    auto ta2_zz_yyyy_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1283);

    auto ta2_zz_yyyy_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1284);

    auto ta2_zz_yyyy_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1285);

    auto ta2_zz_yyyy_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1286);

    auto ta2_zz_yyyy_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1287);

    auto ta2_zz_yyyy_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1288);

    auto ta2_zz_yyyy_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1289);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta2_zz_yy_xxxx_0,   \
                             ta2_zz_yy_xxxx_1,   \
                             ta2_zz_yy_xxxy_0,   \
                             ta2_zz_yy_xxxy_1,   \
                             ta2_zz_yy_xxxz_0,   \
                             ta2_zz_yy_xxxz_1,   \
                             ta2_zz_yy_xxyy_0,   \
                             ta2_zz_yy_xxyy_1,   \
                             ta2_zz_yy_xxyz_0,   \
                             ta2_zz_yy_xxyz_1,   \
                             ta2_zz_yy_xxzz_0,   \
                             ta2_zz_yy_xxzz_1,   \
                             ta2_zz_yy_xyyy_0,   \
                             ta2_zz_yy_xyyy_1,   \
                             ta2_zz_yy_xyyz_0,   \
                             ta2_zz_yy_xyyz_1,   \
                             ta2_zz_yy_xyzz_0,   \
                             ta2_zz_yy_xyzz_1,   \
                             ta2_zz_yy_xzzz_0,   \
                             ta2_zz_yy_xzzz_1,   \
                             ta2_zz_yy_yyyy_0,   \
                             ta2_zz_yy_yyyy_1,   \
                             ta2_zz_yy_yyyz_0,   \
                             ta2_zz_yy_yyyz_1,   \
                             ta2_zz_yy_yyzz_0,   \
                             ta2_zz_yy_yyzz_1,   \
                             ta2_zz_yy_yzzz_0,   \
                             ta2_zz_yy_yzzz_1,   \
                             ta2_zz_yy_zzzz_0,   \
                             ta2_zz_yy_zzzz_1,   \
                             ta2_zz_yyy_xxx_0,   \
                             ta2_zz_yyy_xxx_1,   \
                             ta2_zz_yyy_xxxx_0,  \
                             ta2_zz_yyy_xxxx_1,  \
                             ta2_zz_yyy_xxxy_0,  \
                             ta2_zz_yyy_xxxy_1,  \
                             ta2_zz_yyy_xxxz_0,  \
                             ta2_zz_yyy_xxxz_1,  \
                             ta2_zz_yyy_xxy_0,   \
                             ta2_zz_yyy_xxy_1,   \
                             ta2_zz_yyy_xxyy_0,  \
                             ta2_zz_yyy_xxyy_1,  \
                             ta2_zz_yyy_xxyz_0,  \
                             ta2_zz_yyy_xxyz_1,  \
                             ta2_zz_yyy_xxz_0,   \
                             ta2_zz_yyy_xxz_1,   \
                             ta2_zz_yyy_xxzz_0,  \
                             ta2_zz_yyy_xxzz_1,  \
                             ta2_zz_yyy_xyy_0,   \
                             ta2_zz_yyy_xyy_1,   \
                             ta2_zz_yyy_xyyy_0,  \
                             ta2_zz_yyy_xyyy_1,  \
                             ta2_zz_yyy_xyyz_0,  \
                             ta2_zz_yyy_xyyz_1,  \
                             ta2_zz_yyy_xyz_0,   \
                             ta2_zz_yyy_xyz_1,   \
                             ta2_zz_yyy_xyzz_0,  \
                             ta2_zz_yyy_xyzz_1,  \
                             ta2_zz_yyy_xzz_0,   \
                             ta2_zz_yyy_xzz_1,   \
                             ta2_zz_yyy_xzzz_0,  \
                             ta2_zz_yyy_xzzz_1,  \
                             ta2_zz_yyy_yyy_0,   \
                             ta2_zz_yyy_yyy_1,   \
                             ta2_zz_yyy_yyyy_0,  \
                             ta2_zz_yyy_yyyy_1,  \
                             ta2_zz_yyy_yyyz_0,  \
                             ta2_zz_yyy_yyyz_1,  \
                             ta2_zz_yyy_yyz_0,   \
                             ta2_zz_yyy_yyz_1,   \
                             ta2_zz_yyy_yyzz_0,  \
                             ta2_zz_yyy_yyzz_1,  \
                             ta2_zz_yyy_yzz_0,   \
                             ta2_zz_yyy_yzz_1,   \
                             ta2_zz_yyy_yzzz_0,  \
                             ta2_zz_yyy_yzzz_1,  \
                             ta2_zz_yyy_zzz_0,   \
                             ta2_zz_yyy_zzz_1,   \
                             ta2_zz_yyy_zzzz_0,  \
                             ta2_zz_yyy_zzzz_1,  \
                             ta2_zz_yyyy_xxxx_0, \
                             ta2_zz_yyyy_xxxy_0, \
                             ta2_zz_yyyy_xxxz_0, \
                             ta2_zz_yyyy_xxyy_0, \
                             ta2_zz_yyyy_xxyz_0, \
                             ta2_zz_yyyy_xxzz_0, \
                             ta2_zz_yyyy_xyyy_0, \
                             ta2_zz_yyyy_xyyz_0, \
                             ta2_zz_yyyy_xyzz_0, \
                             ta2_zz_yyyy_xzzz_0, \
                             ta2_zz_yyyy_yyyy_0, \
                             ta2_zz_yyyy_yyyz_0, \
                             ta2_zz_yyyy_yyzz_0, \
                             ta2_zz_yyyy_yzzz_0, \
                             ta2_zz_yyyy_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyyy_xxxx_0[i] =
            3.0 * ta2_zz_yy_xxxx_0[i] * fe_0 - 3.0 * ta2_zz_yy_xxxx_1[i] * fe_0 + ta2_zz_yyy_xxxx_0[i] * pa_y[i] - ta2_zz_yyy_xxxx_1[i] * pc_y[i];

        ta2_zz_yyyy_xxxy_0[i] = 3.0 * ta2_zz_yy_xxxy_0[i] * fe_0 - 3.0 * ta2_zz_yy_xxxy_1[i] * fe_0 + ta2_zz_yyy_xxx_0[i] * fe_0 -
                                ta2_zz_yyy_xxx_1[i] * fe_0 + ta2_zz_yyy_xxxy_0[i] * pa_y[i] - ta2_zz_yyy_xxxy_1[i] * pc_y[i];

        ta2_zz_yyyy_xxxz_0[i] =
            3.0 * ta2_zz_yy_xxxz_0[i] * fe_0 - 3.0 * ta2_zz_yy_xxxz_1[i] * fe_0 + ta2_zz_yyy_xxxz_0[i] * pa_y[i] - ta2_zz_yyy_xxxz_1[i] * pc_y[i];

        ta2_zz_yyyy_xxyy_0[i] = 3.0 * ta2_zz_yy_xxyy_0[i] * fe_0 - 3.0 * ta2_zz_yy_xxyy_1[i] * fe_0 + 2.0 * ta2_zz_yyy_xxy_0[i] * fe_0 -
                                2.0 * ta2_zz_yyy_xxy_1[i] * fe_0 + ta2_zz_yyy_xxyy_0[i] * pa_y[i] - ta2_zz_yyy_xxyy_1[i] * pc_y[i];

        ta2_zz_yyyy_xxyz_0[i] = 3.0 * ta2_zz_yy_xxyz_0[i] * fe_0 - 3.0 * ta2_zz_yy_xxyz_1[i] * fe_0 + ta2_zz_yyy_xxz_0[i] * fe_0 -
                                ta2_zz_yyy_xxz_1[i] * fe_0 + ta2_zz_yyy_xxyz_0[i] * pa_y[i] - ta2_zz_yyy_xxyz_1[i] * pc_y[i];

        ta2_zz_yyyy_xxzz_0[i] =
            3.0 * ta2_zz_yy_xxzz_0[i] * fe_0 - 3.0 * ta2_zz_yy_xxzz_1[i] * fe_0 + ta2_zz_yyy_xxzz_0[i] * pa_y[i] - ta2_zz_yyy_xxzz_1[i] * pc_y[i];

        ta2_zz_yyyy_xyyy_0[i] = 3.0 * ta2_zz_yy_xyyy_0[i] * fe_0 - 3.0 * ta2_zz_yy_xyyy_1[i] * fe_0 + 3.0 * ta2_zz_yyy_xyy_0[i] * fe_0 -
                                3.0 * ta2_zz_yyy_xyy_1[i] * fe_0 + ta2_zz_yyy_xyyy_0[i] * pa_y[i] - ta2_zz_yyy_xyyy_1[i] * pc_y[i];

        ta2_zz_yyyy_xyyz_0[i] = 3.0 * ta2_zz_yy_xyyz_0[i] * fe_0 - 3.0 * ta2_zz_yy_xyyz_1[i] * fe_0 + 2.0 * ta2_zz_yyy_xyz_0[i] * fe_0 -
                                2.0 * ta2_zz_yyy_xyz_1[i] * fe_0 + ta2_zz_yyy_xyyz_0[i] * pa_y[i] - ta2_zz_yyy_xyyz_1[i] * pc_y[i];

        ta2_zz_yyyy_xyzz_0[i] = 3.0 * ta2_zz_yy_xyzz_0[i] * fe_0 - 3.0 * ta2_zz_yy_xyzz_1[i] * fe_0 + ta2_zz_yyy_xzz_0[i] * fe_0 -
                                ta2_zz_yyy_xzz_1[i] * fe_0 + ta2_zz_yyy_xyzz_0[i] * pa_y[i] - ta2_zz_yyy_xyzz_1[i] * pc_y[i];

        ta2_zz_yyyy_xzzz_0[i] =
            3.0 * ta2_zz_yy_xzzz_0[i] * fe_0 - 3.0 * ta2_zz_yy_xzzz_1[i] * fe_0 + ta2_zz_yyy_xzzz_0[i] * pa_y[i] - ta2_zz_yyy_xzzz_1[i] * pc_y[i];

        ta2_zz_yyyy_yyyy_0[i] = 3.0 * ta2_zz_yy_yyyy_0[i] * fe_0 - 3.0 * ta2_zz_yy_yyyy_1[i] * fe_0 + 4.0 * ta2_zz_yyy_yyy_0[i] * fe_0 -
                                4.0 * ta2_zz_yyy_yyy_1[i] * fe_0 + ta2_zz_yyy_yyyy_0[i] * pa_y[i] - ta2_zz_yyy_yyyy_1[i] * pc_y[i];

        ta2_zz_yyyy_yyyz_0[i] = 3.0 * ta2_zz_yy_yyyz_0[i] * fe_0 - 3.0 * ta2_zz_yy_yyyz_1[i] * fe_0 + 3.0 * ta2_zz_yyy_yyz_0[i] * fe_0 -
                                3.0 * ta2_zz_yyy_yyz_1[i] * fe_0 + ta2_zz_yyy_yyyz_0[i] * pa_y[i] - ta2_zz_yyy_yyyz_1[i] * pc_y[i];

        ta2_zz_yyyy_yyzz_0[i] = 3.0 * ta2_zz_yy_yyzz_0[i] * fe_0 - 3.0 * ta2_zz_yy_yyzz_1[i] * fe_0 + 2.0 * ta2_zz_yyy_yzz_0[i] * fe_0 -
                                2.0 * ta2_zz_yyy_yzz_1[i] * fe_0 + ta2_zz_yyy_yyzz_0[i] * pa_y[i] - ta2_zz_yyy_yyzz_1[i] * pc_y[i];

        ta2_zz_yyyy_yzzz_0[i] = 3.0 * ta2_zz_yy_yzzz_0[i] * fe_0 - 3.0 * ta2_zz_yy_yzzz_1[i] * fe_0 + ta2_zz_yyy_zzz_0[i] * fe_0 -
                                ta2_zz_yyy_zzz_1[i] * fe_0 + ta2_zz_yyy_yzzz_0[i] * pa_y[i] - ta2_zz_yyy_yzzz_1[i] * pc_y[i];

        ta2_zz_yyyy_zzzz_0[i] =
            3.0 * ta2_zz_yy_zzzz_0[i] * fe_0 - 3.0 * ta2_zz_yy_zzzz_1[i] * fe_0 + ta2_zz_yyy_zzzz_0[i] * pa_y[i] - ta2_zz_yyy_zzzz_1[i] * pc_y[i];
    }

    // Set up 1290-1305 components of targeted buffer : GG

    auto ta2_zz_yyyz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1290);

    auto ta2_zz_yyyz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1291);

    auto ta2_zz_yyyz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1292);

    auto ta2_zz_yyyz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1293);

    auto ta2_zz_yyyz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1294);

    auto ta2_zz_yyyz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1295);

    auto ta2_zz_yyyz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1296);

    auto ta2_zz_yyyz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1297);

    auto ta2_zz_yyyz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1298);

    auto ta2_zz_yyyz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1299);

    auto ta2_zz_yyyz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1300);

    auto ta2_zz_yyyz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1301);

    auto ta2_zz_yyyz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1302);

    auto ta2_zz_yyyz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1303);

    auto ta2_zz_yyyz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1304);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_yyy_xxxx_1,   \
                             ta1_z_yyy_xxxy_1,   \
                             ta1_z_yyy_xxyy_1,   \
                             ta1_z_yyy_xxyz_1,   \
                             ta1_z_yyy_xyyy_1,   \
                             ta1_z_yyy_xyyz_1,   \
                             ta1_z_yyy_xyzz_1,   \
                             ta1_z_yyy_yyyy_1,   \
                             ta1_z_yyy_yyyz_1,   \
                             ta1_z_yyy_yyzz_1,   \
                             ta1_z_yyy_yzzz_1,   \
                             ta2_zz_yyy_xxxx_0,  \
                             ta2_zz_yyy_xxxx_1,  \
                             ta2_zz_yyy_xxxy_0,  \
                             ta2_zz_yyy_xxxy_1,  \
                             ta2_zz_yyy_xxy_0,   \
                             ta2_zz_yyy_xxy_1,   \
                             ta2_zz_yyy_xxyy_0,  \
                             ta2_zz_yyy_xxyy_1,  \
                             ta2_zz_yyy_xxyz_0,  \
                             ta2_zz_yyy_xxyz_1,  \
                             ta2_zz_yyy_xyy_0,   \
                             ta2_zz_yyy_xyy_1,   \
                             ta2_zz_yyy_xyyy_0,  \
                             ta2_zz_yyy_xyyy_1,  \
                             ta2_zz_yyy_xyyz_0,  \
                             ta2_zz_yyy_xyyz_1,  \
                             ta2_zz_yyy_xyz_0,   \
                             ta2_zz_yyy_xyz_1,   \
                             ta2_zz_yyy_xyzz_0,  \
                             ta2_zz_yyy_xyzz_1,  \
                             ta2_zz_yyy_yyy_0,   \
                             ta2_zz_yyy_yyy_1,   \
                             ta2_zz_yyy_yyyy_0,  \
                             ta2_zz_yyy_yyyy_1,  \
                             ta2_zz_yyy_yyyz_0,  \
                             ta2_zz_yyy_yyyz_1,  \
                             ta2_zz_yyy_yyz_0,   \
                             ta2_zz_yyy_yyz_1,   \
                             ta2_zz_yyy_yyzz_0,  \
                             ta2_zz_yyy_yyzz_1,  \
                             ta2_zz_yyy_yzz_0,   \
                             ta2_zz_yyy_yzz_1,   \
                             ta2_zz_yyy_yzzz_0,  \
                             ta2_zz_yyy_yzzz_1,  \
                             ta2_zz_yyyz_xxxx_0, \
                             ta2_zz_yyyz_xxxy_0, \
                             ta2_zz_yyyz_xxxz_0, \
                             ta2_zz_yyyz_xxyy_0, \
                             ta2_zz_yyyz_xxyz_0, \
                             ta2_zz_yyyz_xxzz_0, \
                             ta2_zz_yyyz_xyyy_0, \
                             ta2_zz_yyyz_xyyz_0, \
                             ta2_zz_yyyz_xyzz_0, \
                             ta2_zz_yyyz_xzzz_0, \
                             ta2_zz_yyyz_yyyy_0, \
                             ta2_zz_yyyz_yyyz_0, \
                             ta2_zz_yyyz_yyzz_0, \
                             ta2_zz_yyyz_yzzz_0, \
                             ta2_zz_yyyz_zzzz_0, \
                             ta2_zz_yyz_xxxz_0,  \
                             ta2_zz_yyz_xxxz_1,  \
                             ta2_zz_yyz_xxzz_0,  \
                             ta2_zz_yyz_xxzz_1,  \
                             ta2_zz_yyz_xzzz_0,  \
                             ta2_zz_yyz_xzzz_1,  \
                             ta2_zz_yyz_zzzz_0,  \
                             ta2_zz_yyz_zzzz_1,  \
                             ta2_zz_yz_xxxz_0,   \
                             ta2_zz_yz_xxxz_1,   \
                             ta2_zz_yz_xxzz_0,   \
                             ta2_zz_yz_xxzz_1,   \
                             ta2_zz_yz_xzzz_0,   \
                             ta2_zz_yz_xzzz_1,   \
                             ta2_zz_yz_zzzz_0,   \
                             ta2_zz_yz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyyz_xxxx_0[i] = 2.0 * ta1_z_yyy_xxxx_1[i] + ta2_zz_yyy_xxxx_0[i] * pa_z[i] - ta2_zz_yyy_xxxx_1[i] * pc_z[i];

        ta2_zz_yyyz_xxxy_0[i] = 2.0 * ta1_z_yyy_xxxy_1[i] + ta2_zz_yyy_xxxy_0[i] * pa_z[i] - ta2_zz_yyy_xxxy_1[i] * pc_z[i];

        ta2_zz_yyyz_xxxz_0[i] =
            2.0 * ta2_zz_yz_xxxz_0[i] * fe_0 - 2.0 * ta2_zz_yz_xxxz_1[i] * fe_0 + ta2_zz_yyz_xxxz_0[i] * pa_y[i] - ta2_zz_yyz_xxxz_1[i] * pc_y[i];

        ta2_zz_yyyz_xxyy_0[i] = 2.0 * ta1_z_yyy_xxyy_1[i] + ta2_zz_yyy_xxyy_0[i] * pa_z[i] - ta2_zz_yyy_xxyy_1[i] * pc_z[i];

        ta2_zz_yyyz_xxyz_0[i] = ta2_zz_yyy_xxy_0[i] * fe_0 - ta2_zz_yyy_xxy_1[i] * fe_0 + 2.0 * ta1_z_yyy_xxyz_1[i] + ta2_zz_yyy_xxyz_0[i] * pa_z[i] -
                                ta2_zz_yyy_xxyz_1[i] * pc_z[i];

        ta2_zz_yyyz_xxzz_0[i] =
            2.0 * ta2_zz_yz_xxzz_0[i] * fe_0 - 2.0 * ta2_zz_yz_xxzz_1[i] * fe_0 + ta2_zz_yyz_xxzz_0[i] * pa_y[i] - ta2_zz_yyz_xxzz_1[i] * pc_y[i];

        ta2_zz_yyyz_xyyy_0[i] = 2.0 * ta1_z_yyy_xyyy_1[i] + ta2_zz_yyy_xyyy_0[i] * pa_z[i] - ta2_zz_yyy_xyyy_1[i] * pc_z[i];

        ta2_zz_yyyz_xyyz_0[i] = ta2_zz_yyy_xyy_0[i] * fe_0 - ta2_zz_yyy_xyy_1[i] * fe_0 + 2.0 * ta1_z_yyy_xyyz_1[i] + ta2_zz_yyy_xyyz_0[i] * pa_z[i] -
                                ta2_zz_yyy_xyyz_1[i] * pc_z[i];

        ta2_zz_yyyz_xyzz_0[i] = 2.0 * ta2_zz_yyy_xyz_0[i] * fe_0 - 2.0 * ta2_zz_yyy_xyz_1[i] * fe_0 + 2.0 * ta1_z_yyy_xyzz_1[i] +
                                ta2_zz_yyy_xyzz_0[i] * pa_z[i] - ta2_zz_yyy_xyzz_1[i] * pc_z[i];

        ta2_zz_yyyz_xzzz_0[i] =
            2.0 * ta2_zz_yz_xzzz_0[i] * fe_0 - 2.0 * ta2_zz_yz_xzzz_1[i] * fe_0 + ta2_zz_yyz_xzzz_0[i] * pa_y[i] - ta2_zz_yyz_xzzz_1[i] * pc_y[i];

        ta2_zz_yyyz_yyyy_0[i] = 2.0 * ta1_z_yyy_yyyy_1[i] + ta2_zz_yyy_yyyy_0[i] * pa_z[i] - ta2_zz_yyy_yyyy_1[i] * pc_z[i];

        ta2_zz_yyyz_yyyz_0[i] = ta2_zz_yyy_yyy_0[i] * fe_0 - ta2_zz_yyy_yyy_1[i] * fe_0 + 2.0 * ta1_z_yyy_yyyz_1[i] + ta2_zz_yyy_yyyz_0[i] * pa_z[i] -
                                ta2_zz_yyy_yyyz_1[i] * pc_z[i];

        ta2_zz_yyyz_yyzz_0[i] = 2.0 * ta2_zz_yyy_yyz_0[i] * fe_0 - 2.0 * ta2_zz_yyy_yyz_1[i] * fe_0 + 2.0 * ta1_z_yyy_yyzz_1[i] +
                                ta2_zz_yyy_yyzz_0[i] * pa_z[i] - ta2_zz_yyy_yyzz_1[i] * pc_z[i];

        ta2_zz_yyyz_yzzz_0[i] = 3.0 * ta2_zz_yyy_yzz_0[i] * fe_0 - 3.0 * ta2_zz_yyy_yzz_1[i] * fe_0 + 2.0 * ta1_z_yyy_yzzz_1[i] +
                                ta2_zz_yyy_yzzz_0[i] * pa_z[i] - ta2_zz_yyy_yzzz_1[i] * pc_z[i];

        ta2_zz_yyyz_zzzz_0[i] =
            2.0 * ta2_zz_yz_zzzz_0[i] * fe_0 - 2.0 * ta2_zz_yz_zzzz_1[i] * fe_0 + ta2_zz_yyz_zzzz_0[i] * pa_y[i] - ta2_zz_yyz_zzzz_1[i] * pc_y[i];
    }

    // Set up 1305-1320 components of targeted buffer : GG

    auto ta2_zz_yyzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1305);

    auto ta2_zz_yyzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1306);

    auto ta2_zz_yyzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1307);

    auto ta2_zz_yyzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1308);

    auto ta2_zz_yyzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1309);

    auto ta2_zz_yyzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1310);

    auto ta2_zz_yyzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1311);

    auto ta2_zz_yyzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1312);

    auto ta2_zz_yyzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1313);

    auto ta2_zz_yyzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1314);

    auto ta2_zz_yyzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1315);

    auto ta2_zz_yyzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1316);

    auto ta2_zz_yyzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1317);

    auto ta2_zz_yyzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1318);

    auto ta2_zz_yyzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1319);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_yyz_xxxy_1,   \
                             ta1_z_yyz_xxyy_1,   \
                             ta1_z_yyz_xyyy_1,   \
                             ta1_z_yyz_yyyy_1,   \
                             ta2_zz_yy_xxxy_0,   \
                             ta2_zz_yy_xxxy_1,   \
                             ta2_zz_yy_xxyy_0,   \
                             ta2_zz_yy_xxyy_1,   \
                             ta2_zz_yy_xyyy_0,   \
                             ta2_zz_yy_xyyy_1,   \
                             ta2_zz_yy_yyyy_0,   \
                             ta2_zz_yy_yyyy_1,   \
                             ta2_zz_yyz_xxxy_0,  \
                             ta2_zz_yyz_xxxy_1,  \
                             ta2_zz_yyz_xxyy_0,  \
                             ta2_zz_yyz_xxyy_1,  \
                             ta2_zz_yyz_xyyy_0,  \
                             ta2_zz_yyz_xyyy_1,  \
                             ta2_zz_yyz_yyyy_0,  \
                             ta2_zz_yyz_yyyy_1,  \
                             ta2_zz_yyzz_xxxx_0, \
                             ta2_zz_yyzz_xxxy_0, \
                             ta2_zz_yyzz_xxxz_0, \
                             ta2_zz_yyzz_xxyy_0, \
                             ta2_zz_yyzz_xxyz_0, \
                             ta2_zz_yyzz_xxzz_0, \
                             ta2_zz_yyzz_xyyy_0, \
                             ta2_zz_yyzz_xyyz_0, \
                             ta2_zz_yyzz_xyzz_0, \
                             ta2_zz_yyzz_xzzz_0, \
                             ta2_zz_yyzz_yyyy_0, \
                             ta2_zz_yyzz_yyyz_0, \
                             ta2_zz_yyzz_yyzz_0, \
                             ta2_zz_yyzz_yzzz_0, \
                             ta2_zz_yyzz_zzzz_0, \
                             ta2_zz_yzz_xxxx_0,  \
                             ta2_zz_yzz_xxxx_1,  \
                             ta2_zz_yzz_xxxz_0,  \
                             ta2_zz_yzz_xxxz_1,  \
                             ta2_zz_yzz_xxyz_0,  \
                             ta2_zz_yzz_xxyz_1,  \
                             ta2_zz_yzz_xxz_0,   \
                             ta2_zz_yzz_xxz_1,   \
                             ta2_zz_yzz_xxzz_0,  \
                             ta2_zz_yzz_xxzz_1,  \
                             ta2_zz_yzz_xyyz_0,  \
                             ta2_zz_yzz_xyyz_1,  \
                             ta2_zz_yzz_xyz_0,   \
                             ta2_zz_yzz_xyz_1,   \
                             ta2_zz_yzz_xyzz_0,  \
                             ta2_zz_yzz_xyzz_1,  \
                             ta2_zz_yzz_xzz_0,   \
                             ta2_zz_yzz_xzz_1,   \
                             ta2_zz_yzz_xzzz_0,  \
                             ta2_zz_yzz_xzzz_1,  \
                             ta2_zz_yzz_yyyz_0,  \
                             ta2_zz_yzz_yyyz_1,  \
                             ta2_zz_yzz_yyz_0,   \
                             ta2_zz_yzz_yyz_1,   \
                             ta2_zz_yzz_yyzz_0,  \
                             ta2_zz_yzz_yyzz_1,  \
                             ta2_zz_yzz_yzz_0,   \
                             ta2_zz_yzz_yzz_1,   \
                             ta2_zz_yzz_yzzz_0,  \
                             ta2_zz_yzz_yzzz_1,  \
                             ta2_zz_yzz_zzz_0,   \
                             ta2_zz_yzz_zzz_1,   \
                             ta2_zz_yzz_zzzz_0,  \
                             ta2_zz_yzz_zzzz_1,  \
                             ta2_zz_zz_xxxx_0,   \
                             ta2_zz_zz_xxxx_1,   \
                             ta2_zz_zz_xxxz_0,   \
                             ta2_zz_zz_xxxz_1,   \
                             ta2_zz_zz_xxyz_0,   \
                             ta2_zz_zz_xxyz_1,   \
                             ta2_zz_zz_xxzz_0,   \
                             ta2_zz_zz_xxzz_1,   \
                             ta2_zz_zz_xyyz_0,   \
                             ta2_zz_zz_xyyz_1,   \
                             ta2_zz_zz_xyzz_0,   \
                             ta2_zz_zz_xyzz_1,   \
                             ta2_zz_zz_xzzz_0,   \
                             ta2_zz_zz_xzzz_1,   \
                             ta2_zz_zz_yyyz_0,   \
                             ta2_zz_zz_yyyz_1,   \
                             ta2_zz_zz_yyzz_0,   \
                             ta2_zz_zz_yyzz_1,   \
                             ta2_zz_zz_yzzz_0,   \
                             ta2_zz_zz_yzzz_1,   \
                             ta2_zz_zz_zzzz_0,   \
                             ta2_zz_zz_zzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yyzz_xxxx_0[i] =
            ta2_zz_zz_xxxx_0[i] * fe_0 - ta2_zz_zz_xxxx_1[i] * fe_0 + ta2_zz_yzz_xxxx_0[i] * pa_y[i] - ta2_zz_yzz_xxxx_1[i] * pc_y[i];

        ta2_zz_yyzz_xxxy_0[i] = ta2_zz_yy_xxxy_0[i] * fe_0 - ta2_zz_yy_xxxy_1[i] * fe_0 + 2.0 * ta1_z_yyz_xxxy_1[i] + ta2_zz_yyz_xxxy_0[i] * pa_z[i] -
                                ta2_zz_yyz_xxxy_1[i] * pc_z[i];

        ta2_zz_yyzz_xxxz_0[i] =
            ta2_zz_zz_xxxz_0[i] * fe_0 - ta2_zz_zz_xxxz_1[i] * fe_0 + ta2_zz_yzz_xxxz_0[i] * pa_y[i] - ta2_zz_yzz_xxxz_1[i] * pc_y[i];

        ta2_zz_yyzz_xxyy_0[i] = ta2_zz_yy_xxyy_0[i] * fe_0 - ta2_zz_yy_xxyy_1[i] * fe_0 + 2.0 * ta1_z_yyz_xxyy_1[i] + ta2_zz_yyz_xxyy_0[i] * pa_z[i] -
                                ta2_zz_yyz_xxyy_1[i] * pc_z[i];

        ta2_zz_yyzz_xxyz_0[i] = ta2_zz_zz_xxyz_0[i] * fe_0 - ta2_zz_zz_xxyz_1[i] * fe_0 + ta2_zz_yzz_xxz_0[i] * fe_0 - ta2_zz_yzz_xxz_1[i] * fe_0 +
                                ta2_zz_yzz_xxyz_0[i] * pa_y[i] - ta2_zz_yzz_xxyz_1[i] * pc_y[i];

        ta2_zz_yyzz_xxzz_0[i] =
            ta2_zz_zz_xxzz_0[i] * fe_0 - ta2_zz_zz_xxzz_1[i] * fe_0 + ta2_zz_yzz_xxzz_0[i] * pa_y[i] - ta2_zz_yzz_xxzz_1[i] * pc_y[i];

        ta2_zz_yyzz_xyyy_0[i] = ta2_zz_yy_xyyy_0[i] * fe_0 - ta2_zz_yy_xyyy_1[i] * fe_0 + 2.0 * ta1_z_yyz_xyyy_1[i] + ta2_zz_yyz_xyyy_0[i] * pa_z[i] -
                                ta2_zz_yyz_xyyy_1[i] * pc_z[i];

        ta2_zz_yyzz_xyyz_0[i] = ta2_zz_zz_xyyz_0[i] * fe_0 - ta2_zz_zz_xyyz_1[i] * fe_0 + 2.0 * ta2_zz_yzz_xyz_0[i] * fe_0 -
                                2.0 * ta2_zz_yzz_xyz_1[i] * fe_0 + ta2_zz_yzz_xyyz_0[i] * pa_y[i] - ta2_zz_yzz_xyyz_1[i] * pc_y[i];

        ta2_zz_yyzz_xyzz_0[i] = ta2_zz_zz_xyzz_0[i] * fe_0 - ta2_zz_zz_xyzz_1[i] * fe_0 + ta2_zz_yzz_xzz_0[i] * fe_0 - ta2_zz_yzz_xzz_1[i] * fe_0 +
                                ta2_zz_yzz_xyzz_0[i] * pa_y[i] - ta2_zz_yzz_xyzz_1[i] * pc_y[i];

        ta2_zz_yyzz_xzzz_0[i] =
            ta2_zz_zz_xzzz_0[i] * fe_0 - ta2_zz_zz_xzzz_1[i] * fe_0 + ta2_zz_yzz_xzzz_0[i] * pa_y[i] - ta2_zz_yzz_xzzz_1[i] * pc_y[i];

        ta2_zz_yyzz_yyyy_0[i] = ta2_zz_yy_yyyy_0[i] * fe_0 - ta2_zz_yy_yyyy_1[i] * fe_0 + 2.0 * ta1_z_yyz_yyyy_1[i] + ta2_zz_yyz_yyyy_0[i] * pa_z[i] -
                                ta2_zz_yyz_yyyy_1[i] * pc_z[i];

        ta2_zz_yyzz_yyyz_0[i] = ta2_zz_zz_yyyz_0[i] * fe_0 - ta2_zz_zz_yyyz_1[i] * fe_0 + 3.0 * ta2_zz_yzz_yyz_0[i] * fe_0 -
                                3.0 * ta2_zz_yzz_yyz_1[i] * fe_0 + ta2_zz_yzz_yyyz_0[i] * pa_y[i] - ta2_zz_yzz_yyyz_1[i] * pc_y[i];

        ta2_zz_yyzz_yyzz_0[i] = ta2_zz_zz_yyzz_0[i] * fe_0 - ta2_zz_zz_yyzz_1[i] * fe_0 + 2.0 * ta2_zz_yzz_yzz_0[i] * fe_0 -
                                2.0 * ta2_zz_yzz_yzz_1[i] * fe_0 + ta2_zz_yzz_yyzz_0[i] * pa_y[i] - ta2_zz_yzz_yyzz_1[i] * pc_y[i];

        ta2_zz_yyzz_yzzz_0[i] = ta2_zz_zz_yzzz_0[i] * fe_0 - ta2_zz_zz_yzzz_1[i] * fe_0 + ta2_zz_yzz_zzz_0[i] * fe_0 - ta2_zz_yzz_zzz_1[i] * fe_0 +
                                ta2_zz_yzz_yzzz_0[i] * pa_y[i] - ta2_zz_yzz_yzzz_1[i] * pc_y[i];

        ta2_zz_yyzz_zzzz_0[i] =
            ta2_zz_zz_zzzz_0[i] * fe_0 - ta2_zz_zz_zzzz_1[i] * fe_0 + ta2_zz_yzz_zzzz_0[i] * pa_y[i] - ta2_zz_yzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 1320-1335 components of targeted buffer : GG

    auto ta2_zz_yzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1320);

    auto ta2_zz_yzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1321);

    auto ta2_zz_yzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1322);

    auto ta2_zz_yzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1323);

    auto ta2_zz_yzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1324);

    auto ta2_zz_yzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1325);

    auto ta2_zz_yzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1326);

    auto ta2_zz_yzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1327);

    auto ta2_zz_yzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1328);

    auto ta2_zz_yzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1329);

    auto ta2_zz_yzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1330);

    auto ta2_zz_yzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1331);

    auto ta2_zz_yzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1332);

    auto ta2_zz_yzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1333);

    auto ta2_zz_yzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1334);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta2_zz_yzzz_xxxx_0, \
                             ta2_zz_yzzz_xxxy_0, \
                             ta2_zz_yzzz_xxxz_0, \
                             ta2_zz_yzzz_xxyy_0, \
                             ta2_zz_yzzz_xxyz_0, \
                             ta2_zz_yzzz_xxzz_0, \
                             ta2_zz_yzzz_xyyy_0, \
                             ta2_zz_yzzz_xyyz_0, \
                             ta2_zz_yzzz_xyzz_0, \
                             ta2_zz_yzzz_xzzz_0, \
                             ta2_zz_yzzz_yyyy_0, \
                             ta2_zz_yzzz_yyyz_0, \
                             ta2_zz_yzzz_yyzz_0, \
                             ta2_zz_yzzz_yzzz_0, \
                             ta2_zz_yzzz_zzzz_0, \
                             ta2_zz_zzz_xxx_0,   \
                             ta2_zz_zzz_xxx_1,   \
                             ta2_zz_zzz_xxxx_0,  \
                             ta2_zz_zzz_xxxx_1,  \
                             ta2_zz_zzz_xxxy_0,  \
                             ta2_zz_zzz_xxxy_1,  \
                             ta2_zz_zzz_xxxz_0,  \
                             ta2_zz_zzz_xxxz_1,  \
                             ta2_zz_zzz_xxy_0,   \
                             ta2_zz_zzz_xxy_1,   \
                             ta2_zz_zzz_xxyy_0,  \
                             ta2_zz_zzz_xxyy_1,  \
                             ta2_zz_zzz_xxyz_0,  \
                             ta2_zz_zzz_xxyz_1,  \
                             ta2_zz_zzz_xxz_0,   \
                             ta2_zz_zzz_xxz_1,   \
                             ta2_zz_zzz_xxzz_0,  \
                             ta2_zz_zzz_xxzz_1,  \
                             ta2_zz_zzz_xyy_0,   \
                             ta2_zz_zzz_xyy_1,   \
                             ta2_zz_zzz_xyyy_0,  \
                             ta2_zz_zzz_xyyy_1,  \
                             ta2_zz_zzz_xyyz_0,  \
                             ta2_zz_zzz_xyyz_1,  \
                             ta2_zz_zzz_xyz_0,   \
                             ta2_zz_zzz_xyz_1,   \
                             ta2_zz_zzz_xyzz_0,  \
                             ta2_zz_zzz_xyzz_1,  \
                             ta2_zz_zzz_xzz_0,   \
                             ta2_zz_zzz_xzz_1,   \
                             ta2_zz_zzz_xzzz_0,  \
                             ta2_zz_zzz_xzzz_1,  \
                             ta2_zz_zzz_yyy_0,   \
                             ta2_zz_zzz_yyy_1,   \
                             ta2_zz_zzz_yyyy_0,  \
                             ta2_zz_zzz_yyyy_1,  \
                             ta2_zz_zzz_yyyz_0,  \
                             ta2_zz_zzz_yyyz_1,  \
                             ta2_zz_zzz_yyz_0,   \
                             ta2_zz_zzz_yyz_1,   \
                             ta2_zz_zzz_yyzz_0,  \
                             ta2_zz_zzz_yyzz_1,  \
                             ta2_zz_zzz_yzz_0,   \
                             ta2_zz_zzz_yzz_1,   \
                             ta2_zz_zzz_yzzz_0,  \
                             ta2_zz_zzz_yzzz_1,  \
                             ta2_zz_zzz_zzz_0,   \
                             ta2_zz_zzz_zzz_1,   \
                             ta2_zz_zzz_zzzz_0,  \
                             ta2_zz_zzz_zzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_yzzz_xxxx_0[i] = ta2_zz_zzz_xxxx_0[i] * pa_y[i] - ta2_zz_zzz_xxxx_1[i] * pc_y[i];

        ta2_zz_yzzz_xxxy_0[i] =
            ta2_zz_zzz_xxx_0[i] * fe_0 - ta2_zz_zzz_xxx_1[i] * fe_0 + ta2_zz_zzz_xxxy_0[i] * pa_y[i] - ta2_zz_zzz_xxxy_1[i] * pc_y[i];

        ta2_zz_yzzz_xxxz_0[i] = ta2_zz_zzz_xxxz_0[i] * pa_y[i] - ta2_zz_zzz_xxxz_1[i] * pc_y[i];

        ta2_zz_yzzz_xxyy_0[i] =
            2.0 * ta2_zz_zzz_xxy_0[i] * fe_0 - 2.0 * ta2_zz_zzz_xxy_1[i] * fe_0 + ta2_zz_zzz_xxyy_0[i] * pa_y[i] - ta2_zz_zzz_xxyy_1[i] * pc_y[i];

        ta2_zz_yzzz_xxyz_0[i] =
            ta2_zz_zzz_xxz_0[i] * fe_0 - ta2_zz_zzz_xxz_1[i] * fe_0 + ta2_zz_zzz_xxyz_0[i] * pa_y[i] - ta2_zz_zzz_xxyz_1[i] * pc_y[i];

        ta2_zz_yzzz_xxzz_0[i] = ta2_zz_zzz_xxzz_0[i] * pa_y[i] - ta2_zz_zzz_xxzz_1[i] * pc_y[i];

        ta2_zz_yzzz_xyyy_0[i] =
            3.0 * ta2_zz_zzz_xyy_0[i] * fe_0 - 3.0 * ta2_zz_zzz_xyy_1[i] * fe_0 + ta2_zz_zzz_xyyy_0[i] * pa_y[i] - ta2_zz_zzz_xyyy_1[i] * pc_y[i];

        ta2_zz_yzzz_xyyz_0[i] =
            2.0 * ta2_zz_zzz_xyz_0[i] * fe_0 - 2.0 * ta2_zz_zzz_xyz_1[i] * fe_0 + ta2_zz_zzz_xyyz_0[i] * pa_y[i] - ta2_zz_zzz_xyyz_1[i] * pc_y[i];

        ta2_zz_yzzz_xyzz_0[i] =
            ta2_zz_zzz_xzz_0[i] * fe_0 - ta2_zz_zzz_xzz_1[i] * fe_0 + ta2_zz_zzz_xyzz_0[i] * pa_y[i] - ta2_zz_zzz_xyzz_1[i] * pc_y[i];

        ta2_zz_yzzz_xzzz_0[i] = ta2_zz_zzz_xzzz_0[i] * pa_y[i] - ta2_zz_zzz_xzzz_1[i] * pc_y[i];

        ta2_zz_yzzz_yyyy_0[i] =
            4.0 * ta2_zz_zzz_yyy_0[i] * fe_0 - 4.0 * ta2_zz_zzz_yyy_1[i] * fe_0 + ta2_zz_zzz_yyyy_0[i] * pa_y[i] - ta2_zz_zzz_yyyy_1[i] * pc_y[i];

        ta2_zz_yzzz_yyyz_0[i] =
            3.0 * ta2_zz_zzz_yyz_0[i] * fe_0 - 3.0 * ta2_zz_zzz_yyz_1[i] * fe_0 + ta2_zz_zzz_yyyz_0[i] * pa_y[i] - ta2_zz_zzz_yyyz_1[i] * pc_y[i];

        ta2_zz_yzzz_yyzz_0[i] =
            2.0 * ta2_zz_zzz_yzz_0[i] * fe_0 - 2.0 * ta2_zz_zzz_yzz_1[i] * fe_0 + ta2_zz_zzz_yyzz_0[i] * pa_y[i] - ta2_zz_zzz_yyzz_1[i] * pc_y[i];

        ta2_zz_yzzz_yzzz_0[i] =
            ta2_zz_zzz_zzz_0[i] * fe_0 - ta2_zz_zzz_zzz_1[i] * fe_0 + ta2_zz_zzz_yzzz_0[i] * pa_y[i] - ta2_zz_zzz_yzzz_1[i] * pc_y[i];

        ta2_zz_yzzz_zzzz_0[i] = ta2_zz_zzz_zzzz_0[i] * pa_y[i] - ta2_zz_zzz_zzzz_1[i] * pc_y[i];
    }

    // Set up 1335-1350 components of targeted buffer : GG

    auto ta2_zz_zzzz_xxxx_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1335);

    auto ta2_zz_zzzz_xxxy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1336);

    auto ta2_zz_zzzz_xxxz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1337);

    auto ta2_zz_zzzz_xxyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1338);

    auto ta2_zz_zzzz_xxyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1339);

    auto ta2_zz_zzzz_xxzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1340);

    auto ta2_zz_zzzz_xyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1341);

    auto ta2_zz_zzzz_xyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1342);

    auto ta2_zz_zzzz_xyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1343);

    auto ta2_zz_zzzz_xzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1344);

    auto ta2_zz_zzzz_yyyy_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1345);

    auto ta2_zz_zzzz_yyyz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1346);

    auto ta2_zz_zzzz_yyzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1347);

    auto ta2_zz_zzzz_yzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1348);

    auto ta2_zz_zzzz_zzzz_0 = pbuffer.data(idx_npot_geom_020_0_gg + 1349);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_z_zzz_xxxx_1,   \
                             ta1_z_zzz_xxxy_1,   \
                             ta1_z_zzz_xxxz_1,   \
                             ta1_z_zzz_xxyy_1,   \
                             ta1_z_zzz_xxyz_1,   \
                             ta1_z_zzz_xxzz_1,   \
                             ta1_z_zzz_xyyy_1,   \
                             ta1_z_zzz_xyyz_1,   \
                             ta1_z_zzz_xyzz_1,   \
                             ta1_z_zzz_xzzz_1,   \
                             ta1_z_zzz_yyyy_1,   \
                             ta1_z_zzz_yyyz_1,   \
                             ta1_z_zzz_yyzz_1,   \
                             ta1_z_zzz_yzzz_1,   \
                             ta1_z_zzz_zzzz_1,   \
                             ta2_zz_zz_xxxx_0,   \
                             ta2_zz_zz_xxxx_1,   \
                             ta2_zz_zz_xxxy_0,   \
                             ta2_zz_zz_xxxy_1,   \
                             ta2_zz_zz_xxxz_0,   \
                             ta2_zz_zz_xxxz_1,   \
                             ta2_zz_zz_xxyy_0,   \
                             ta2_zz_zz_xxyy_1,   \
                             ta2_zz_zz_xxyz_0,   \
                             ta2_zz_zz_xxyz_1,   \
                             ta2_zz_zz_xxzz_0,   \
                             ta2_zz_zz_xxzz_1,   \
                             ta2_zz_zz_xyyy_0,   \
                             ta2_zz_zz_xyyy_1,   \
                             ta2_zz_zz_xyyz_0,   \
                             ta2_zz_zz_xyyz_1,   \
                             ta2_zz_zz_xyzz_0,   \
                             ta2_zz_zz_xyzz_1,   \
                             ta2_zz_zz_xzzz_0,   \
                             ta2_zz_zz_xzzz_1,   \
                             ta2_zz_zz_yyyy_0,   \
                             ta2_zz_zz_yyyy_1,   \
                             ta2_zz_zz_yyyz_0,   \
                             ta2_zz_zz_yyyz_1,   \
                             ta2_zz_zz_yyzz_0,   \
                             ta2_zz_zz_yyzz_1,   \
                             ta2_zz_zz_yzzz_0,   \
                             ta2_zz_zz_yzzz_1,   \
                             ta2_zz_zz_zzzz_0,   \
                             ta2_zz_zz_zzzz_1,   \
                             ta2_zz_zzz_xxx_0,   \
                             ta2_zz_zzz_xxx_1,   \
                             ta2_zz_zzz_xxxx_0,  \
                             ta2_zz_zzz_xxxx_1,  \
                             ta2_zz_zzz_xxxy_0,  \
                             ta2_zz_zzz_xxxy_1,  \
                             ta2_zz_zzz_xxxz_0,  \
                             ta2_zz_zzz_xxxz_1,  \
                             ta2_zz_zzz_xxy_0,   \
                             ta2_zz_zzz_xxy_1,   \
                             ta2_zz_zzz_xxyy_0,  \
                             ta2_zz_zzz_xxyy_1,  \
                             ta2_zz_zzz_xxyz_0,  \
                             ta2_zz_zzz_xxyz_1,  \
                             ta2_zz_zzz_xxz_0,   \
                             ta2_zz_zzz_xxz_1,   \
                             ta2_zz_zzz_xxzz_0,  \
                             ta2_zz_zzz_xxzz_1,  \
                             ta2_zz_zzz_xyy_0,   \
                             ta2_zz_zzz_xyy_1,   \
                             ta2_zz_zzz_xyyy_0,  \
                             ta2_zz_zzz_xyyy_1,  \
                             ta2_zz_zzz_xyyz_0,  \
                             ta2_zz_zzz_xyyz_1,  \
                             ta2_zz_zzz_xyz_0,   \
                             ta2_zz_zzz_xyz_1,   \
                             ta2_zz_zzz_xyzz_0,  \
                             ta2_zz_zzz_xyzz_1,  \
                             ta2_zz_zzz_xzz_0,   \
                             ta2_zz_zzz_xzz_1,   \
                             ta2_zz_zzz_xzzz_0,  \
                             ta2_zz_zzz_xzzz_1,  \
                             ta2_zz_zzz_yyy_0,   \
                             ta2_zz_zzz_yyy_1,   \
                             ta2_zz_zzz_yyyy_0,  \
                             ta2_zz_zzz_yyyy_1,  \
                             ta2_zz_zzz_yyyz_0,  \
                             ta2_zz_zzz_yyyz_1,  \
                             ta2_zz_zzz_yyz_0,   \
                             ta2_zz_zzz_yyz_1,   \
                             ta2_zz_zzz_yyzz_0,  \
                             ta2_zz_zzz_yyzz_1,  \
                             ta2_zz_zzz_yzz_0,   \
                             ta2_zz_zzz_yzz_1,   \
                             ta2_zz_zzz_yzzz_0,  \
                             ta2_zz_zzz_yzzz_1,  \
                             ta2_zz_zzz_zzz_0,   \
                             ta2_zz_zzz_zzz_1,   \
                             ta2_zz_zzz_zzzz_0,  \
                             ta2_zz_zzz_zzzz_1,  \
                             ta2_zz_zzzz_xxxx_0, \
                             ta2_zz_zzzz_xxxy_0, \
                             ta2_zz_zzzz_xxxz_0, \
                             ta2_zz_zzzz_xxyy_0, \
                             ta2_zz_zzzz_xxyz_0, \
                             ta2_zz_zzzz_xxzz_0, \
                             ta2_zz_zzzz_xyyy_0, \
                             ta2_zz_zzzz_xyyz_0, \
                             ta2_zz_zzzz_xyzz_0, \
                             ta2_zz_zzzz_xzzz_0, \
                             ta2_zz_zzzz_yyyy_0, \
                             ta2_zz_zzzz_yyyz_0, \
                             ta2_zz_zzzz_yyzz_0, \
                             ta2_zz_zzzz_yzzz_0, \
                             ta2_zz_zzzz_zzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta2_zz_zzzz_xxxx_0[i] = 3.0 * ta2_zz_zz_xxxx_0[i] * fe_0 - 3.0 * ta2_zz_zz_xxxx_1[i] * fe_0 + 2.0 * ta1_z_zzz_xxxx_1[i] +
                                ta2_zz_zzz_xxxx_0[i] * pa_z[i] - ta2_zz_zzz_xxxx_1[i] * pc_z[i];

        ta2_zz_zzzz_xxxy_0[i] = 3.0 * ta2_zz_zz_xxxy_0[i] * fe_0 - 3.0 * ta2_zz_zz_xxxy_1[i] * fe_0 + 2.0 * ta1_z_zzz_xxxy_1[i] +
                                ta2_zz_zzz_xxxy_0[i] * pa_z[i] - ta2_zz_zzz_xxxy_1[i] * pc_z[i];

        ta2_zz_zzzz_xxxz_0[i] = 3.0 * ta2_zz_zz_xxxz_0[i] * fe_0 - 3.0 * ta2_zz_zz_xxxz_1[i] * fe_0 + ta2_zz_zzz_xxx_0[i] * fe_0 -
                                ta2_zz_zzz_xxx_1[i] * fe_0 + 2.0 * ta1_z_zzz_xxxz_1[i] + ta2_zz_zzz_xxxz_0[i] * pa_z[i] -
                                ta2_zz_zzz_xxxz_1[i] * pc_z[i];

        ta2_zz_zzzz_xxyy_0[i] = 3.0 * ta2_zz_zz_xxyy_0[i] * fe_0 - 3.0 * ta2_zz_zz_xxyy_1[i] * fe_0 + 2.0 * ta1_z_zzz_xxyy_1[i] +
                                ta2_zz_zzz_xxyy_0[i] * pa_z[i] - ta2_zz_zzz_xxyy_1[i] * pc_z[i];

        ta2_zz_zzzz_xxyz_0[i] = 3.0 * ta2_zz_zz_xxyz_0[i] * fe_0 - 3.0 * ta2_zz_zz_xxyz_1[i] * fe_0 + ta2_zz_zzz_xxy_0[i] * fe_0 -
                                ta2_zz_zzz_xxy_1[i] * fe_0 + 2.0 * ta1_z_zzz_xxyz_1[i] + ta2_zz_zzz_xxyz_0[i] * pa_z[i] -
                                ta2_zz_zzz_xxyz_1[i] * pc_z[i];

        ta2_zz_zzzz_xxzz_0[i] = 3.0 * ta2_zz_zz_xxzz_0[i] * fe_0 - 3.0 * ta2_zz_zz_xxzz_1[i] * fe_0 + 2.0 * ta2_zz_zzz_xxz_0[i] * fe_0 -
                                2.0 * ta2_zz_zzz_xxz_1[i] * fe_0 + 2.0 * ta1_z_zzz_xxzz_1[i] + ta2_zz_zzz_xxzz_0[i] * pa_z[i] -
                                ta2_zz_zzz_xxzz_1[i] * pc_z[i];

        ta2_zz_zzzz_xyyy_0[i] = 3.0 * ta2_zz_zz_xyyy_0[i] * fe_0 - 3.0 * ta2_zz_zz_xyyy_1[i] * fe_0 + 2.0 * ta1_z_zzz_xyyy_1[i] +
                                ta2_zz_zzz_xyyy_0[i] * pa_z[i] - ta2_zz_zzz_xyyy_1[i] * pc_z[i];

        ta2_zz_zzzz_xyyz_0[i] = 3.0 * ta2_zz_zz_xyyz_0[i] * fe_0 - 3.0 * ta2_zz_zz_xyyz_1[i] * fe_0 + ta2_zz_zzz_xyy_0[i] * fe_0 -
                                ta2_zz_zzz_xyy_1[i] * fe_0 + 2.0 * ta1_z_zzz_xyyz_1[i] + ta2_zz_zzz_xyyz_0[i] * pa_z[i] -
                                ta2_zz_zzz_xyyz_1[i] * pc_z[i];

        ta2_zz_zzzz_xyzz_0[i] = 3.0 * ta2_zz_zz_xyzz_0[i] * fe_0 - 3.0 * ta2_zz_zz_xyzz_1[i] * fe_0 + 2.0 * ta2_zz_zzz_xyz_0[i] * fe_0 -
                                2.0 * ta2_zz_zzz_xyz_1[i] * fe_0 + 2.0 * ta1_z_zzz_xyzz_1[i] + ta2_zz_zzz_xyzz_0[i] * pa_z[i] -
                                ta2_zz_zzz_xyzz_1[i] * pc_z[i];

        ta2_zz_zzzz_xzzz_0[i] = 3.0 * ta2_zz_zz_xzzz_0[i] * fe_0 - 3.0 * ta2_zz_zz_xzzz_1[i] * fe_0 + 3.0 * ta2_zz_zzz_xzz_0[i] * fe_0 -
                                3.0 * ta2_zz_zzz_xzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_xzzz_1[i] + ta2_zz_zzz_xzzz_0[i] * pa_z[i] -
                                ta2_zz_zzz_xzzz_1[i] * pc_z[i];

        ta2_zz_zzzz_yyyy_0[i] = 3.0 * ta2_zz_zz_yyyy_0[i] * fe_0 - 3.0 * ta2_zz_zz_yyyy_1[i] * fe_0 + 2.0 * ta1_z_zzz_yyyy_1[i] +
                                ta2_zz_zzz_yyyy_0[i] * pa_z[i] - ta2_zz_zzz_yyyy_1[i] * pc_z[i];

        ta2_zz_zzzz_yyyz_0[i] = 3.0 * ta2_zz_zz_yyyz_0[i] * fe_0 - 3.0 * ta2_zz_zz_yyyz_1[i] * fe_0 + ta2_zz_zzz_yyy_0[i] * fe_0 -
                                ta2_zz_zzz_yyy_1[i] * fe_0 + 2.0 * ta1_z_zzz_yyyz_1[i] + ta2_zz_zzz_yyyz_0[i] * pa_z[i] -
                                ta2_zz_zzz_yyyz_1[i] * pc_z[i];

        ta2_zz_zzzz_yyzz_0[i] = 3.0 * ta2_zz_zz_yyzz_0[i] * fe_0 - 3.0 * ta2_zz_zz_yyzz_1[i] * fe_0 + 2.0 * ta2_zz_zzz_yyz_0[i] * fe_0 -
                                2.0 * ta2_zz_zzz_yyz_1[i] * fe_0 + 2.0 * ta1_z_zzz_yyzz_1[i] + ta2_zz_zzz_yyzz_0[i] * pa_z[i] -
                                ta2_zz_zzz_yyzz_1[i] * pc_z[i];

        ta2_zz_zzzz_yzzz_0[i] = 3.0 * ta2_zz_zz_yzzz_0[i] * fe_0 - 3.0 * ta2_zz_zz_yzzz_1[i] * fe_0 + 3.0 * ta2_zz_zzz_yzz_0[i] * fe_0 -
                                3.0 * ta2_zz_zzz_yzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_yzzz_1[i] + ta2_zz_zzz_yzzz_0[i] * pa_z[i] -
                                ta2_zz_zzz_yzzz_1[i] * pc_z[i];

        ta2_zz_zzzz_zzzz_0[i] = 3.0 * ta2_zz_zz_zzzz_0[i] * fe_0 - 3.0 * ta2_zz_zz_zzzz_1[i] * fe_0 + 4.0 * ta2_zz_zzz_zzz_0[i] * fe_0 -
                                4.0 * ta2_zz_zzz_zzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_zzzz_1[i] + ta2_zz_zzz_zzzz_0[i] * pa_z[i] -
                                ta2_zz_zzz_zzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
