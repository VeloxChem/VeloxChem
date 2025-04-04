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

#include "NuclearPotentialGeom010PrimRecGH.hpp"

namespace npotrec {  // npotrec namespace

auto
comp_prim_nuclear_potential_geom_010_gh(CSimdArray<double>&       pbuffer,
                                        const size_t              idx_npot_geom_010_0_gh,
                                        const size_t              idx_npot_geom_010_0_dh,
                                        const size_t              idx_npot_geom_010_1_dh,
                                        const size_t              idx_npot_geom_010_0_fg,
                                        const size_t              idx_npot_geom_010_1_fg,
                                        const size_t              idx_npot_1_fh,
                                        const size_t              idx_npot_geom_010_0_fh,
                                        const size_t              idx_npot_geom_010_1_fh,
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

    // Set up components of auxiliary buffer : DH

    auto ta1_x_xx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh);

    auto ta1_x_xx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 1);

    auto ta1_x_xx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 2);

    auto ta1_x_xx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 3);

    auto ta1_x_xx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 4);

    auto ta1_x_xx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 5);

    auto ta1_x_xx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 6);

    auto ta1_x_xx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 7);

    auto ta1_x_xx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 8);

    auto ta1_x_xx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 9);

    auto ta1_x_xx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 10);

    auto ta1_x_xx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 11);

    auto ta1_x_xx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 12);

    auto ta1_x_xx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 13);

    auto ta1_x_xx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 14);

    auto ta1_x_xx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 15);

    auto ta1_x_xx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 16);

    auto ta1_x_xx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 17);

    auto ta1_x_xx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 18);

    auto ta1_x_xx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 19);

    auto ta1_x_xx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 20);

    auto ta1_x_xy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 21);

    auto ta1_x_xy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 23);

    auto ta1_x_xy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 26);

    auto ta1_x_xy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 30);

    auto ta1_x_xy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 35);

    auto ta1_x_xz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 42);

    auto ta1_x_xz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 43);

    auto ta1_x_xz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 44);

    auto ta1_x_xz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 45);

    auto ta1_x_xz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 47);

    auto ta1_x_xz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 48);

    auto ta1_x_xz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 51);

    auto ta1_x_xz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 52);

    auto ta1_x_xz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 56);

    auto ta1_x_yy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 63);

    auto ta1_x_yy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 64);

    auto ta1_x_yy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 65);

    auto ta1_x_yy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 66);

    auto ta1_x_yy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 67);

    auto ta1_x_yy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 68);

    auto ta1_x_yy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 69);

    auto ta1_x_yy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 70);

    auto ta1_x_yy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 71);

    auto ta1_x_yy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 72);

    auto ta1_x_yy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 73);

    auto ta1_x_yy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 74);

    auto ta1_x_yy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 75);

    auto ta1_x_yy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 76);

    auto ta1_x_yy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 77);

    auto ta1_x_yy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 78);

    auto ta1_x_yy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 79);

    auto ta1_x_yy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 80);

    auto ta1_x_yy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 81);

    auto ta1_x_yy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 82);

    auto ta1_x_yy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 83);

    auto ta1_x_yz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 86);

    auto ta1_x_yz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 89);

    auto ta1_x_yz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 93);

    auto ta1_x_yz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 98);

    auto ta1_x_yz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 104);

    auto ta1_x_zz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 105);

    auto ta1_x_zz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 106);

    auto ta1_x_zz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 107);

    auto ta1_x_zz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 108);

    auto ta1_x_zz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 109);

    auto ta1_x_zz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 110);

    auto ta1_x_zz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 111);

    auto ta1_x_zz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 112);

    auto ta1_x_zz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 113);

    auto ta1_x_zz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 114);

    auto ta1_x_zz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 115);

    auto ta1_x_zz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 116);

    auto ta1_x_zz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 117);

    auto ta1_x_zz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 118);

    auto ta1_x_zz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 119);

    auto ta1_x_zz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 120);

    auto ta1_x_zz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 121);

    auto ta1_x_zz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 122);

    auto ta1_x_zz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 123);

    auto ta1_x_zz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 124);

    auto ta1_x_zz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 125);

    auto ta1_y_xx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 126);

    auto ta1_y_xx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 127);

    auto ta1_y_xx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 128);

    auto ta1_y_xx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 129);

    auto ta1_y_xx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 130);

    auto ta1_y_xx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 131);

    auto ta1_y_xx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 132);

    auto ta1_y_xx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 133);

    auto ta1_y_xx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 134);

    auto ta1_y_xx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 135);

    auto ta1_y_xx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 136);

    auto ta1_y_xx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 137);

    auto ta1_y_xx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 138);

    auto ta1_y_xx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 139);

    auto ta1_y_xx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 140);

    auto ta1_y_xx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 141);

    auto ta1_y_xx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 142);

    auto ta1_y_xx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 143);

    auto ta1_y_xx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 144);

    auto ta1_y_xx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 145);

    auto ta1_y_xx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 146);

    auto ta1_y_xy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 162);

    auto ta1_y_xy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 163);

    auto ta1_y_xy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 164);

    auto ta1_y_xy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 165);

    auto ta1_y_xy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 166);

    auto ta1_y_xz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 184);

    auto ta1_y_xz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 185);

    auto ta1_y_xz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 186);

    auto ta1_y_xz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 187);

    auto ta1_y_xz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 188);

    auto ta1_y_yy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 189);

    auto ta1_y_yy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 190);

    auto ta1_y_yy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 191);

    auto ta1_y_yy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 192);

    auto ta1_y_yy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 193);

    auto ta1_y_yy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 194);

    auto ta1_y_yy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 195);

    auto ta1_y_yy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 196);

    auto ta1_y_yy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 197);

    auto ta1_y_yy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 198);

    auto ta1_y_yy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 199);

    auto ta1_y_yy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 200);

    auto ta1_y_yy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 201);

    auto ta1_y_yy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 202);

    auto ta1_y_yy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 203);

    auto ta1_y_yy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 204);

    auto ta1_y_yy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 205);

    auto ta1_y_yy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 206);

    auto ta1_y_yy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 207);

    auto ta1_y_yy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 208);

    auto ta1_y_yy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 209);

    auto ta1_y_yz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 211);

    auto ta1_y_yz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 213);

    auto ta1_y_yz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 216);

    auto ta1_y_yz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 220);

    auto ta1_y_yz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 225);

    auto ta1_y_yz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 226);

    auto ta1_y_yz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 227);

    auto ta1_y_yz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 228);

    auto ta1_y_yz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 229);

    auto ta1_y_zz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 231);

    auto ta1_y_zz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 232);

    auto ta1_y_zz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 233);

    auto ta1_y_zz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 234);

    auto ta1_y_zz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 235);

    auto ta1_y_zz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 236);

    auto ta1_y_zz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 237);

    auto ta1_y_zz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 238);

    auto ta1_y_zz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 239);

    auto ta1_y_zz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 240);

    auto ta1_y_zz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 241);

    auto ta1_y_zz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 242);

    auto ta1_y_zz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 243);

    auto ta1_y_zz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 244);

    auto ta1_y_zz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 245);

    auto ta1_y_zz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 246);

    auto ta1_y_zz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 247);

    auto ta1_y_zz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 248);

    auto ta1_y_zz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 249);

    auto ta1_y_zz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 250);

    auto ta1_y_zz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 251);

    auto ta1_z_xx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 252);

    auto ta1_z_xx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 253);

    auto ta1_z_xx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 254);

    auto ta1_z_xx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 255);

    auto ta1_z_xx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 256);

    auto ta1_z_xx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 257);

    auto ta1_z_xx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 258);

    auto ta1_z_xx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 259);

    auto ta1_z_xx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 260);

    auto ta1_z_xx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 261);

    auto ta1_z_xx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 262);

    auto ta1_z_xx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 263);

    auto ta1_z_xx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 264);

    auto ta1_z_xx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 265);

    auto ta1_z_xx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 266);

    auto ta1_z_xx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 267);

    auto ta1_z_xx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 268);

    auto ta1_z_xx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 269);

    auto ta1_z_xx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 270);

    auto ta1_z_xx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 271);

    auto ta1_z_xx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 272);

    auto ta1_z_xy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 288);

    auto ta1_z_xy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 289);

    auto ta1_z_xy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 290);

    auto ta1_z_xy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 291);

    auto ta1_z_xy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 292);

    auto ta1_z_xz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 310);

    auto ta1_z_xz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 311);

    auto ta1_z_xz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 312);

    auto ta1_z_xz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 313);

    auto ta1_z_xz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 314);

    auto ta1_z_yy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 315);

    auto ta1_z_yy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 316);

    auto ta1_z_yy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 317);

    auto ta1_z_yy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 318);

    auto ta1_z_yy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 319);

    auto ta1_z_yy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 320);

    auto ta1_z_yy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 321);

    auto ta1_z_yy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 322);

    auto ta1_z_yy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 323);

    auto ta1_z_yy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 324);

    auto ta1_z_yy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 325);

    auto ta1_z_yy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 326);

    auto ta1_z_yy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 327);

    auto ta1_z_yy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 328);

    auto ta1_z_yy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 329);

    auto ta1_z_yy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 330);

    auto ta1_z_yy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 331);

    auto ta1_z_yy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 332);

    auto ta1_z_yy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 333);

    auto ta1_z_yy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 334);

    auto ta1_z_yy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 335);

    auto ta1_z_yz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 338);

    auto ta1_z_yz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 341);

    auto ta1_z_yz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 345);

    auto ta1_z_yz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 350);

    auto ta1_z_yz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 352);

    auto ta1_z_yz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 353);

    auto ta1_z_yz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 354);

    auto ta1_z_yz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 355);

    auto ta1_z_yz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 356);

    auto ta1_z_zz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_dh + 357);

    auto ta1_z_zz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 358);

    auto ta1_z_zz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 359);

    auto ta1_z_zz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 360);

    auto ta1_z_zz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 361);

    auto ta1_z_zz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 362);

    auto ta1_z_zz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 363);

    auto ta1_z_zz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 364);

    auto ta1_z_zz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 365);

    auto ta1_z_zz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 366);

    auto ta1_z_zz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 367);

    auto ta1_z_zz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 368);

    auto ta1_z_zz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 369);

    auto ta1_z_zz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 370);

    auto ta1_z_zz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 371);

    auto ta1_z_zz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_dh + 372);

    auto ta1_z_zz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 373);

    auto ta1_z_zz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 374);

    auto ta1_z_zz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 375);

    auto ta1_z_zz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 376);

    auto ta1_z_zz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_dh + 377);

    // Set up components of auxiliary buffer : DH

    auto ta1_x_xx_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh);

    auto ta1_x_xx_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 1);

    auto ta1_x_xx_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 2);

    auto ta1_x_xx_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 3);

    auto ta1_x_xx_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 4);

    auto ta1_x_xx_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 5);

    auto ta1_x_xx_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 6);

    auto ta1_x_xx_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 7);

    auto ta1_x_xx_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 8);

    auto ta1_x_xx_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 9);

    auto ta1_x_xx_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 10);

    auto ta1_x_xx_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 11);

    auto ta1_x_xx_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 12);

    auto ta1_x_xx_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 13);

    auto ta1_x_xx_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 14);

    auto ta1_x_xx_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 15);

    auto ta1_x_xx_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 16);

    auto ta1_x_xx_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 17);

    auto ta1_x_xx_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 18);

    auto ta1_x_xx_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 19);

    auto ta1_x_xx_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 20);

    auto ta1_x_xy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 21);

    auto ta1_x_xy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 23);

    auto ta1_x_xy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 26);

    auto ta1_x_xy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 30);

    auto ta1_x_xy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 35);

    auto ta1_x_xz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 42);

    auto ta1_x_xz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 43);

    auto ta1_x_xz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 44);

    auto ta1_x_xz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 45);

    auto ta1_x_xz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 47);

    auto ta1_x_xz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 48);

    auto ta1_x_xz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 51);

    auto ta1_x_xz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 52);

    auto ta1_x_xz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 56);

    auto ta1_x_yy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 63);

    auto ta1_x_yy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 64);

    auto ta1_x_yy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 65);

    auto ta1_x_yy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 66);

    auto ta1_x_yy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 67);

    auto ta1_x_yy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 68);

    auto ta1_x_yy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 69);

    auto ta1_x_yy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 70);

    auto ta1_x_yy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 71);

    auto ta1_x_yy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 72);

    auto ta1_x_yy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 73);

    auto ta1_x_yy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 74);

    auto ta1_x_yy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 75);

    auto ta1_x_yy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 76);

    auto ta1_x_yy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 77);

    auto ta1_x_yy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 78);

    auto ta1_x_yy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 79);

    auto ta1_x_yy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 80);

    auto ta1_x_yy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 81);

    auto ta1_x_yy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 82);

    auto ta1_x_yy_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 83);

    auto ta1_x_yz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 86);

    auto ta1_x_yz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 89);

    auto ta1_x_yz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 93);

    auto ta1_x_yz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 98);

    auto ta1_x_yz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 104);

    auto ta1_x_zz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 105);

    auto ta1_x_zz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 106);

    auto ta1_x_zz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 107);

    auto ta1_x_zz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 108);

    auto ta1_x_zz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 109);

    auto ta1_x_zz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 110);

    auto ta1_x_zz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 111);

    auto ta1_x_zz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 112);

    auto ta1_x_zz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 113);

    auto ta1_x_zz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 114);

    auto ta1_x_zz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 115);

    auto ta1_x_zz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 116);

    auto ta1_x_zz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 117);

    auto ta1_x_zz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 118);

    auto ta1_x_zz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 119);

    auto ta1_x_zz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 120);

    auto ta1_x_zz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 121);

    auto ta1_x_zz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 122);

    auto ta1_x_zz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 123);

    auto ta1_x_zz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 124);

    auto ta1_x_zz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 125);

    auto ta1_y_xx_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 126);

    auto ta1_y_xx_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 127);

    auto ta1_y_xx_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 128);

    auto ta1_y_xx_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 129);

    auto ta1_y_xx_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 130);

    auto ta1_y_xx_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 131);

    auto ta1_y_xx_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 132);

    auto ta1_y_xx_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 133);

    auto ta1_y_xx_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 134);

    auto ta1_y_xx_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 135);

    auto ta1_y_xx_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 136);

    auto ta1_y_xx_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 137);

    auto ta1_y_xx_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 138);

    auto ta1_y_xx_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 139);

    auto ta1_y_xx_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 140);

    auto ta1_y_xx_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 141);

    auto ta1_y_xx_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 142);

    auto ta1_y_xx_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 143);

    auto ta1_y_xx_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 144);

    auto ta1_y_xx_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 145);

    auto ta1_y_xx_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 146);

    auto ta1_y_xy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 162);

    auto ta1_y_xy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 163);

    auto ta1_y_xy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 164);

    auto ta1_y_xy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 165);

    auto ta1_y_xy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 166);

    auto ta1_y_xz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 184);

    auto ta1_y_xz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 185);

    auto ta1_y_xz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 186);

    auto ta1_y_xz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 187);

    auto ta1_y_xz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 188);

    auto ta1_y_yy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 189);

    auto ta1_y_yy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 190);

    auto ta1_y_yy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 191);

    auto ta1_y_yy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 192);

    auto ta1_y_yy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 193);

    auto ta1_y_yy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 194);

    auto ta1_y_yy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 195);

    auto ta1_y_yy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 196);

    auto ta1_y_yy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 197);

    auto ta1_y_yy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 198);

    auto ta1_y_yy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 199);

    auto ta1_y_yy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 200);

    auto ta1_y_yy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 201);

    auto ta1_y_yy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 202);

    auto ta1_y_yy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 203);

    auto ta1_y_yy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 204);

    auto ta1_y_yy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 205);

    auto ta1_y_yy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 206);

    auto ta1_y_yy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 207);

    auto ta1_y_yy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 208);

    auto ta1_y_yy_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 209);

    auto ta1_y_yz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 211);

    auto ta1_y_yz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 213);

    auto ta1_y_yz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 216);

    auto ta1_y_yz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 220);

    auto ta1_y_yz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 225);

    auto ta1_y_yz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 226);

    auto ta1_y_yz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 227);

    auto ta1_y_yz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 228);

    auto ta1_y_yz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 229);

    auto ta1_y_zz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 231);

    auto ta1_y_zz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 232);

    auto ta1_y_zz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 233);

    auto ta1_y_zz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 234);

    auto ta1_y_zz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 235);

    auto ta1_y_zz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 236);

    auto ta1_y_zz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 237);

    auto ta1_y_zz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 238);

    auto ta1_y_zz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 239);

    auto ta1_y_zz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 240);

    auto ta1_y_zz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 241);

    auto ta1_y_zz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 242);

    auto ta1_y_zz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 243);

    auto ta1_y_zz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 244);

    auto ta1_y_zz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 245);

    auto ta1_y_zz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 246);

    auto ta1_y_zz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 247);

    auto ta1_y_zz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 248);

    auto ta1_y_zz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 249);

    auto ta1_y_zz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 250);

    auto ta1_y_zz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 251);

    auto ta1_z_xx_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 252);

    auto ta1_z_xx_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 253);

    auto ta1_z_xx_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 254);

    auto ta1_z_xx_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 255);

    auto ta1_z_xx_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 256);

    auto ta1_z_xx_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 257);

    auto ta1_z_xx_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 258);

    auto ta1_z_xx_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 259);

    auto ta1_z_xx_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 260);

    auto ta1_z_xx_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 261);

    auto ta1_z_xx_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 262);

    auto ta1_z_xx_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 263);

    auto ta1_z_xx_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 264);

    auto ta1_z_xx_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 265);

    auto ta1_z_xx_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 266);

    auto ta1_z_xx_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 267);

    auto ta1_z_xx_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 268);

    auto ta1_z_xx_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 269);

    auto ta1_z_xx_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 270);

    auto ta1_z_xx_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 271);

    auto ta1_z_xx_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 272);

    auto ta1_z_xy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 288);

    auto ta1_z_xy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 289);

    auto ta1_z_xy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 290);

    auto ta1_z_xy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 291);

    auto ta1_z_xy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 292);

    auto ta1_z_xz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 310);

    auto ta1_z_xz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 311);

    auto ta1_z_xz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 312);

    auto ta1_z_xz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 313);

    auto ta1_z_xz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 314);

    auto ta1_z_yy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 315);

    auto ta1_z_yy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 316);

    auto ta1_z_yy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 317);

    auto ta1_z_yy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 318);

    auto ta1_z_yy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 319);

    auto ta1_z_yy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 320);

    auto ta1_z_yy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 321);

    auto ta1_z_yy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 322);

    auto ta1_z_yy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 323);

    auto ta1_z_yy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 324);

    auto ta1_z_yy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 325);

    auto ta1_z_yy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 326);

    auto ta1_z_yy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 327);

    auto ta1_z_yy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 328);

    auto ta1_z_yy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 329);

    auto ta1_z_yy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 330);

    auto ta1_z_yy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 331);

    auto ta1_z_yy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 332);

    auto ta1_z_yy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 333);

    auto ta1_z_yy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 334);

    auto ta1_z_yy_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 335);

    auto ta1_z_yz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 338);

    auto ta1_z_yz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 341);

    auto ta1_z_yz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 345);

    auto ta1_z_yz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 350);

    auto ta1_z_yz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 352);

    auto ta1_z_yz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 353);

    auto ta1_z_yz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 354);

    auto ta1_z_yz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 355);

    auto ta1_z_yz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 356);

    auto ta1_z_zz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_dh + 357);

    auto ta1_z_zz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 358);

    auto ta1_z_zz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 359);

    auto ta1_z_zz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 360);

    auto ta1_z_zz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 361);

    auto ta1_z_zz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 362);

    auto ta1_z_zz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 363);

    auto ta1_z_zz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 364);

    auto ta1_z_zz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 365);

    auto ta1_z_zz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 366);

    auto ta1_z_zz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 367);

    auto ta1_z_zz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 368);

    auto ta1_z_zz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 369);

    auto ta1_z_zz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 370);

    auto ta1_z_zz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 371);

    auto ta1_z_zz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_dh + 372);

    auto ta1_z_zz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 373);

    auto ta1_z_zz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 374);

    auto ta1_z_zz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 375);

    auto ta1_z_zz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 376);

    auto ta1_z_zz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_dh + 377);

    // Set up components of auxiliary buffer : FG

    auto ta1_x_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg);

    auto ta1_x_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 1);

    auto ta1_x_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 2);

    auto ta1_x_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 3);

    auto ta1_x_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 4);

    auto ta1_x_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 5);

    auto ta1_x_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 6);

    auto ta1_x_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 7);

    auto ta1_x_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 8);

    auto ta1_x_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 9);

    auto ta1_x_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 10);

    auto ta1_x_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 11);

    auto ta1_x_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 12);

    auto ta1_x_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 13);

    auto ta1_x_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 14);

    auto ta1_x_xxy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 15);

    auto ta1_x_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 16);

    auto ta1_x_xxy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 17);

    auto ta1_x_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 18);

    auto ta1_x_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 19);

    auto ta1_x_xxy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 20);

    auto ta1_x_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 21);

    auto ta1_x_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 22);

    auto ta1_x_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 23);

    auto ta1_x_xxy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 24);

    auto ta1_x_xxz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 30);

    auto ta1_x_xxz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 31);

    auto ta1_x_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 32);

    auto ta1_x_xxz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 33);

    auto ta1_x_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 34);

    auto ta1_x_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 35);

    auto ta1_x_xxz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 36);

    auto ta1_x_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 37);

    auto ta1_x_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 38);

    auto ta1_x_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 39);

    auto ta1_x_xxz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 41);

    auto ta1_x_xxz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 42);

    auto ta1_x_xxz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 43);

    auto ta1_x_xxz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 44);

    auto ta1_x_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 46);

    auto ta1_x_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 48);

    auto ta1_x_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 49);

    auto ta1_x_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 51);

    auto ta1_x_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 52);

    auto ta1_x_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 53);

    auto ta1_x_xzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 75);

    auto ta1_x_xzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 76);

    auto ta1_x_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 77);

    auto ta1_x_xzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 78);

    auto ta1_x_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 79);

    auto ta1_x_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 80);

    auto ta1_x_xzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 81);

    auto ta1_x_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 82);

    auto ta1_x_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 83);

    auto ta1_x_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 84);

    auto ta1_x_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 90);

    auto ta1_x_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 91);

    auto ta1_x_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 92);

    auto ta1_x_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 93);

    auto ta1_x_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 94);

    auto ta1_x_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 95);

    auto ta1_x_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 96);

    auto ta1_x_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 97);

    auto ta1_x_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 98);

    auto ta1_x_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 99);

    auto ta1_x_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 100);

    auto ta1_x_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 101);

    auto ta1_x_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 102);

    auto ta1_x_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 103);

    auto ta1_x_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 104);

    auto ta1_x_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 122);

    auto ta1_x_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 124);

    auto ta1_x_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 125);

    auto ta1_x_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 127);

    auto ta1_x_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 128);

    auto ta1_x_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 129);

    auto ta1_x_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 131);

    auto ta1_x_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 132);

    auto ta1_x_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 133);

    auto ta1_x_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 134);

    auto ta1_x_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 135);

    auto ta1_x_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 136);

    auto ta1_x_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 137);

    auto ta1_x_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 138);

    auto ta1_x_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 139);

    auto ta1_x_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 140);

    auto ta1_x_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 141);

    auto ta1_x_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 142);

    auto ta1_x_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 143);

    auto ta1_x_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 144);

    auto ta1_x_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 145);

    auto ta1_x_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 146);

    auto ta1_x_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 147);

    auto ta1_x_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 148);

    auto ta1_x_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 149);

    auto ta1_y_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 150);

    auto ta1_y_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 151);

    auto ta1_y_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 152);

    auto ta1_y_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 153);

    auto ta1_y_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 154);

    auto ta1_y_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 155);

    auto ta1_y_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 156);

    auto ta1_y_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 157);

    auto ta1_y_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 158);

    auto ta1_y_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 159);

    auto ta1_y_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 160);

    auto ta1_y_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 161);

    auto ta1_y_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 162);

    auto ta1_y_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 163);

    auto ta1_y_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 164);

    auto ta1_y_xxy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 166);

    auto ta1_y_xxy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 168);

    auto ta1_y_xxy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 169);

    auto ta1_y_xxy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 171);

    auto ta1_y_xxy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 172);

    auto ta1_y_xxy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 173);

    auto ta1_y_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 196);

    auto ta1_y_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 198);

    auto ta1_y_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 199);

    auto ta1_y_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 201);

    auto ta1_y_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 202);

    auto ta1_y_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 203);

    auto ta1_y_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 205);

    auto ta1_y_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 206);

    auto ta1_y_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 207);

    auto ta1_y_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 208);

    auto ta1_y_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 227);

    auto ta1_y_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 229);

    auto ta1_y_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 230);

    auto ta1_y_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 232);

    auto ta1_y_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 233);

    auto ta1_y_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 234);

    auto ta1_y_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 236);

    auto ta1_y_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 237);

    auto ta1_y_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 238);

    auto ta1_y_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 239);

    auto ta1_y_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 240);

    auto ta1_y_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 241);

    auto ta1_y_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 242);

    auto ta1_y_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 243);

    auto ta1_y_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 244);

    auto ta1_y_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 245);

    auto ta1_y_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 246);

    auto ta1_y_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 247);

    auto ta1_y_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 248);

    auto ta1_y_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 249);

    auto ta1_y_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 250);

    auto ta1_y_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 251);

    auto ta1_y_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 252);

    auto ta1_y_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 253);

    auto ta1_y_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 254);

    auto ta1_y_yyz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 256);

    auto ta1_y_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 257);

    auto ta1_y_yyz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 258);

    auto ta1_y_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 259);

    auto ta1_y_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 260);

    auto ta1_y_yyz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 261);

    auto ta1_y_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 262);

    auto ta1_y_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 263);

    auto ta1_y_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 264);

    auto ta1_y_yyz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 265);

    auto ta1_y_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 266);

    auto ta1_y_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 267);

    auto ta1_y_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 268);

    auto ta1_y_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 269);

    auto ta1_y_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 271);

    auto ta1_y_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 273);

    auto ta1_y_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 274);

    auto ta1_y_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 276);

    auto ta1_y_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 277);

    auto ta1_y_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 278);

    auto ta1_y_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 280);

    auto ta1_y_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 281);

    auto ta1_y_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 282);

    auto ta1_y_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 283);

    auto ta1_y_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 285);

    auto ta1_y_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 286);

    auto ta1_y_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 287);

    auto ta1_y_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 288);

    auto ta1_y_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 289);

    auto ta1_y_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 290);

    auto ta1_y_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 291);

    auto ta1_y_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 292);

    auto ta1_y_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 293);

    auto ta1_y_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 294);

    auto ta1_y_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 295);

    auto ta1_y_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 296);

    auto ta1_y_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 297);

    auto ta1_y_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 298);

    auto ta1_y_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 299);

    auto ta1_z_xxx_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 300);

    auto ta1_z_xxx_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 301);

    auto ta1_z_xxx_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 302);

    auto ta1_z_xxx_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 303);

    auto ta1_z_xxx_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 304);

    auto ta1_z_xxx_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 305);

    auto ta1_z_xxx_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 306);

    auto ta1_z_xxx_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 307);

    auto ta1_z_xxx_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 308);

    auto ta1_z_xxx_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 309);

    auto ta1_z_xxx_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 310);

    auto ta1_z_xxx_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 311);

    auto ta1_z_xxx_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 312);

    auto ta1_z_xxx_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 313);

    auto ta1_z_xxx_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 314);

    auto ta1_z_xxz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 332);

    auto ta1_z_xxz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 334);

    auto ta1_z_xxz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 335);

    auto ta1_z_xxz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 337);

    auto ta1_z_xxz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 338);

    auto ta1_z_xxz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 339);

    auto ta1_z_xyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 346);

    auto ta1_z_xyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 348);

    auto ta1_z_xyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 349);

    auto ta1_z_xyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 351);

    auto ta1_z_xyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 352);

    auto ta1_z_xyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 353);

    auto ta1_z_xyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 355);

    auto ta1_z_xyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 356);

    auto ta1_z_xyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 357);

    auto ta1_z_xyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 358);

    auto ta1_z_xzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 377);

    auto ta1_z_xzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 379);

    auto ta1_z_xzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 380);

    auto ta1_z_xzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 382);

    auto ta1_z_xzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 383);

    auto ta1_z_xzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 384);

    auto ta1_z_xzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 386);

    auto ta1_z_xzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 387);

    auto ta1_z_xzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 388);

    auto ta1_z_xzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 389);

    auto ta1_z_yyy_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 390);

    auto ta1_z_yyy_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 391);

    auto ta1_z_yyy_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 392);

    auto ta1_z_yyy_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 393);

    auto ta1_z_yyy_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 394);

    auto ta1_z_yyy_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 395);

    auto ta1_z_yyy_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 396);

    auto ta1_z_yyy_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 397);

    auto ta1_z_yyy_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 398);

    auto ta1_z_yyy_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 399);

    auto ta1_z_yyy_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 400);

    auto ta1_z_yyy_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 401);

    auto ta1_z_yyy_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 402);

    auto ta1_z_yyy_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 403);

    auto ta1_z_yyy_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 404);

    auto ta1_z_yyz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 407);

    auto ta1_z_yyz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 409);

    auto ta1_z_yyz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 410);

    auto ta1_z_yyz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 412);

    auto ta1_z_yyz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 413);

    auto ta1_z_yyz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 414);

    auto ta1_z_yyz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 416);

    auto ta1_z_yyz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 417);

    auto ta1_z_yyz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 418);

    auto ta1_z_yyz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 419);

    auto ta1_z_yzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 421);

    auto ta1_z_yzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 422);

    auto ta1_z_yzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 423);

    auto ta1_z_yzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 424);

    auto ta1_z_yzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 425);

    auto ta1_z_yzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 426);

    auto ta1_z_yzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 427);

    auto ta1_z_yzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 428);

    auto ta1_z_yzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 429);

    auto ta1_z_yzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 430);

    auto ta1_z_yzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 431);

    auto ta1_z_yzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 432);

    auto ta1_z_yzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 433);

    auto ta1_z_yzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 434);

    auto ta1_z_zzz_xxxx_0 = pbuffer.data(idx_npot_geom_010_0_fg + 435);

    auto ta1_z_zzz_xxxy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 436);

    auto ta1_z_zzz_xxxz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 437);

    auto ta1_z_zzz_xxyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 438);

    auto ta1_z_zzz_xxyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 439);

    auto ta1_z_zzz_xxzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 440);

    auto ta1_z_zzz_xyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 441);

    auto ta1_z_zzz_xyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 442);

    auto ta1_z_zzz_xyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 443);

    auto ta1_z_zzz_xzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 444);

    auto ta1_z_zzz_yyyy_0 = pbuffer.data(idx_npot_geom_010_0_fg + 445);

    auto ta1_z_zzz_yyyz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 446);

    auto ta1_z_zzz_yyzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 447);

    auto ta1_z_zzz_yzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 448);

    auto ta1_z_zzz_zzzz_0 = pbuffer.data(idx_npot_geom_010_0_fg + 449);

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

    auto ta1_x_xxy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 19);

    auto ta1_x_xxy_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 20);

    auto ta1_x_xxy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 21);

    auto ta1_x_xxy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 22);

    auto ta1_x_xxy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 23);

    auto ta1_x_xxy_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 24);

    auto ta1_x_xxz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 30);

    auto ta1_x_xxz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 31);

    auto ta1_x_xxz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 32);

    auto ta1_x_xxz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 33);

    auto ta1_x_xxz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 34);

    auto ta1_x_xxz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 35);

    auto ta1_x_xxz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 36);

    auto ta1_x_xxz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 37);

    auto ta1_x_xxz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 38);

    auto ta1_x_xxz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 39);

    auto ta1_x_xxz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 41);

    auto ta1_x_xxz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 42);

    auto ta1_x_xxz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 43);

    auto ta1_x_xxz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 44);

    auto ta1_x_xyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 46);

    auto ta1_x_xyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 48);

    auto ta1_x_xyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 49);

    auto ta1_x_xyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 51);

    auto ta1_x_xyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 52);

    auto ta1_x_xyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 53);

    auto ta1_x_xzz_xxxx_1 = pbuffer.data(idx_npot_geom_010_1_fg + 75);

    auto ta1_x_xzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 76);

    auto ta1_x_xzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 77);

    auto ta1_x_xzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 78);

    auto ta1_x_xzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 79);

    auto ta1_x_xzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 80);

    auto ta1_x_xzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 81);

    auto ta1_x_xzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 82);

    auto ta1_x_xzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 83);

    auto ta1_x_xzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 84);

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

    auto ta1_x_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 122);

    auto ta1_x_yzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 124);

    auto ta1_x_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 125);

    auto ta1_x_yzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 127);

    auto ta1_x_yzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 128);

    auto ta1_x_yzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 129);

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

    auto ta1_y_xxy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 166);

    auto ta1_y_xxy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 168);

    auto ta1_y_xxy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 169);

    auto ta1_y_xxy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 171);

    auto ta1_y_xxy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 172);

    auto ta1_y_xxy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 173);

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

    auto ta1_y_xzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 227);

    auto ta1_y_xzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 229);

    auto ta1_y_xzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 230);

    auto ta1_y_xzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 232);

    auto ta1_y_xzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 233);

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

    auto ta1_y_yyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 257);

    auto ta1_y_yyz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 258);

    auto ta1_y_yyz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 259);

    auto ta1_y_yyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 260);

    auto ta1_y_yyz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 261);

    auto ta1_y_yyz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 262);

    auto ta1_y_yyz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 263);

    auto ta1_y_yyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 264);

    auto ta1_y_yyz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 265);

    auto ta1_y_yyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 266);

    auto ta1_y_yyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 267);

    auto ta1_y_yyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 268);

    auto ta1_y_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 269);

    auto ta1_y_yzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 271);

    auto ta1_y_yzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 273);

    auto ta1_y_yzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 274);

    auto ta1_y_yzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 276);

    auto ta1_y_yzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 277);

    auto ta1_y_yzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 278);

    auto ta1_y_yzz_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 280);

    auto ta1_y_yzz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 281);

    auto ta1_y_yzz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 282);

    auto ta1_y_yzz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 283);

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

    auto ta1_z_xxz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 332);

    auto ta1_z_xxz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 334);

    auto ta1_z_xxz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 335);

    auto ta1_z_xxz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 337);

    auto ta1_z_xxz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 338);

    auto ta1_z_xxz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 339);

    auto ta1_z_xyy_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 346);

    auto ta1_z_xyy_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 348);

    auto ta1_z_xyy_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 349);

    auto ta1_z_xyy_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 351);

    auto ta1_z_xyy_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 352);

    auto ta1_z_xyy_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 353);

    auto ta1_z_xyy_yyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 355);

    auto ta1_z_xyy_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 356);

    auto ta1_z_xyy_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 357);

    auto ta1_z_xyy_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 358);

    auto ta1_z_xzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 377);

    auto ta1_z_xzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 379);

    auto ta1_z_xzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 380);

    auto ta1_z_xzz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 382);

    auto ta1_z_xzz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 383);

    auto ta1_z_xzz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 384);

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

    auto ta1_z_yyz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 407);

    auto ta1_z_yyz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 409);

    auto ta1_z_yyz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 410);

    auto ta1_z_yyz_xyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 412);

    auto ta1_z_yyz_xyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 413);

    auto ta1_z_yyz_xzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 414);

    auto ta1_z_yyz_yyyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 416);

    auto ta1_z_yyz_yyzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 417);

    auto ta1_z_yyz_yzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 418);

    auto ta1_z_yyz_zzzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 419);

    auto ta1_z_yzz_xxxy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 421);

    auto ta1_z_yzz_xxxz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 422);

    auto ta1_z_yzz_xxyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 423);

    auto ta1_z_yzz_xxyz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 424);

    auto ta1_z_yzz_xxzz_1 = pbuffer.data(idx_npot_geom_010_1_fg + 425);

    auto ta1_z_yzz_xyyy_1 = pbuffer.data(idx_npot_geom_010_1_fg + 426);

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

    // Set up components of auxiliary buffer : FH

    auto ta_xxx_xxxxx_1 = pbuffer.data(idx_npot_1_fh);

    auto ta_xxx_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 1);

    auto ta_xxx_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 2);

    auto ta_xxx_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 3);

    auto ta_xxx_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 4);

    auto ta_xxx_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 5);

    auto ta_xxx_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 6);

    auto ta_xxx_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 7);

    auto ta_xxx_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 8);

    auto ta_xxx_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 9);

    auto ta_xxx_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 10);

    auto ta_xxx_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 11);

    auto ta_xxx_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 12);

    auto ta_xxx_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 13);

    auto ta_xxx_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 14);

    auto ta_xxx_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 15);

    auto ta_xxx_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 16);

    auto ta_xxx_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 17);

    auto ta_xxx_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 18);

    auto ta_xxx_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 19);

    auto ta_xxx_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 20);

    auto ta_xxy_xxxxx_1 = pbuffer.data(idx_npot_1_fh + 21);

    auto ta_xxy_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 22);

    auto ta_xxy_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 23);

    auto ta_xxy_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 24);

    auto ta_xxy_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 26);

    auto ta_xxy_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 27);

    auto ta_xxy_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 30);

    auto ta_xxy_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 31);

    auto ta_xxy_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 35);

    auto ta_xxy_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 36);

    auto ta_xxz_xxxxx_1 = pbuffer.data(idx_npot_1_fh + 42);

    auto ta_xxz_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 43);

    auto ta_xxz_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 44);

    auto ta_xxz_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 45);

    auto ta_xxz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 47);

    auto ta_xxz_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 48);

    auto ta_xxz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 51);

    auto ta_xxz_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 52);

    auto ta_xxz_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 56);

    auto ta_xxz_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 62);

    auto ta_xyy_xxxxx_1 = pbuffer.data(idx_npot_1_fh + 63);

    auto ta_xyy_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 64);

    auto ta_xyy_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 66);

    auto ta_xyy_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 69);

    auto ta_xyy_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 73);

    auto ta_xyy_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 78);

    auto ta_xyy_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 79);

    auto ta_xyy_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 80);

    auto ta_xyy_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 81);

    auto ta_xyy_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 82);

    auto ta_xzz_xxxxx_1 = pbuffer.data(idx_npot_1_fh + 105);

    auto ta_xzz_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 107);

    auto ta_xzz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 110);

    auto ta_xzz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 114);

    auto ta_xzz_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 119);

    auto ta_xzz_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 121);

    auto ta_xzz_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 122);

    auto ta_xzz_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 123);

    auto ta_xzz_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 124);

    auto ta_xzz_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 125);

    auto ta_yyy_xxxxx_1 = pbuffer.data(idx_npot_1_fh + 126);

    auto ta_yyy_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 127);

    auto ta_yyy_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 128);

    auto ta_yyy_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 129);

    auto ta_yyy_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 130);

    auto ta_yyy_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 131);

    auto ta_yyy_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 132);

    auto ta_yyy_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 133);

    auto ta_yyy_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 134);

    auto ta_yyy_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 135);

    auto ta_yyy_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 136);

    auto ta_yyy_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 137);

    auto ta_yyy_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 138);

    auto ta_yyy_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 139);

    auto ta_yyy_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 140);

    auto ta_yyy_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 141);

    auto ta_yyy_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 142);

    auto ta_yyy_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 143);

    auto ta_yyy_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 144);

    auto ta_yyy_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 145);

    auto ta_yyy_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 146);

    auto ta_yyz_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 148);

    auto ta_yyz_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 150);

    auto ta_yyz_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 153);

    auto ta_yyz_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 157);

    auto ta_yyz_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 162);

    auto ta_yyz_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 163);

    auto ta_yyz_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 164);

    auto ta_yyz_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 165);

    auto ta_yyz_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 166);

    auto ta_yyz_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 167);

    auto ta_yzz_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 170);

    auto ta_yzz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 173);

    auto ta_yzz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 177);

    auto ta_yzz_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 182);

    auto ta_yzz_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 183);

    auto ta_yzz_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 184);

    auto ta_yzz_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 185);

    auto ta_yzz_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 186);

    auto ta_yzz_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 187);

    auto ta_yzz_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 188);

    auto ta_zzz_xxxxx_1 = pbuffer.data(idx_npot_1_fh + 189);

    auto ta_zzz_xxxxy_1 = pbuffer.data(idx_npot_1_fh + 190);

    auto ta_zzz_xxxxz_1 = pbuffer.data(idx_npot_1_fh + 191);

    auto ta_zzz_xxxyy_1 = pbuffer.data(idx_npot_1_fh + 192);

    auto ta_zzz_xxxyz_1 = pbuffer.data(idx_npot_1_fh + 193);

    auto ta_zzz_xxxzz_1 = pbuffer.data(idx_npot_1_fh + 194);

    auto ta_zzz_xxyyy_1 = pbuffer.data(idx_npot_1_fh + 195);

    auto ta_zzz_xxyyz_1 = pbuffer.data(idx_npot_1_fh + 196);

    auto ta_zzz_xxyzz_1 = pbuffer.data(idx_npot_1_fh + 197);

    auto ta_zzz_xxzzz_1 = pbuffer.data(idx_npot_1_fh + 198);

    auto ta_zzz_xyyyy_1 = pbuffer.data(idx_npot_1_fh + 199);

    auto ta_zzz_xyyyz_1 = pbuffer.data(idx_npot_1_fh + 200);

    auto ta_zzz_xyyzz_1 = pbuffer.data(idx_npot_1_fh + 201);

    auto ta_zzz_xyzzz_1 = pbuffer.data(idx_npot_1_fh + 202);

    auto ta_zzz_xzzzz_1 = pbuffer.data(idx_npot_1_fh + 203);

    auto ta_zzz_yyyyy_1 = pbuffer.data(idx_npot_1_fh + 204);

    auto ta_zzz_yyyyz_1 = pbuffer.data(idx_npot_1_fh + 205);

    auto ta_zzz_yyyzz_1 = pbuffer.data(idx_npot_1_fh + 206);

    auto ta_zzz_yyzzz_1 = pbuffer.data(idx_npot_1_fh + 207);

    auto ta_zzz_yzzzz_1 = pbuffer.data(idx_npot_1_fh + 208);

    auto ta_zzz_zzzzz_1 = pbuffer.data(idx_npot_1_fh + 209);

    // Set up components of auxiliary buffer : FH

    auto ta1_x_xxx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh);

    auto ta1_x_xxx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 1);

    auto ta1_x_xxx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 2);

    auto ta1_x_xxx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 3);

    auto ta1_x_xxx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 4);

    auto ta1_x_xxx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 5);

    auto ta1_x_xxx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 6);

    auto ta1_x_xxx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 7);

    auto ta1_x_xxx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 8);

    auto ta1_x_xxx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 9);

    auto ta1_x_xxx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 10);

    auto ta1_x_xxx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 11);

    auto ta1_x_xxx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 12);

    auto ta1_x_xxx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 13);

    auto ta1_x_xxx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 14);

    auto ta1_x_xxx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 15);

    auto ta1_x_xxx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 16);

    auto ta1_x_xxx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 17);

    auto ta1_x_xxx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 18);

    auto ta1_x_xxx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 19);

    auto ta1_x_xxx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 20);

    auto ta1_x_xxy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 21);

    auto ta1_x_xxy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 22);

    auto ta1_x_xxy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 23);

    auto ta1_x_xxy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 24);

    auto ta1_x_xxy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 25);

    auto ta1_x_xxy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 26);

    auto ta1_x_xxy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 27);

    auto ta1_x_xxy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 28);

    auto ta1_x_xxy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 29);

    auto ta1_x_xxy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 30);

    auto ta1_x_xxy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 31);

    auto ta1_x_xxy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 32);

    auto ta1_x_xxy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 33);

    auto ta1_x_xxy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 34);

    auto ta1_x_xxy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 35);

    auto ta1_x_xxy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 36);

    auto ta1_x_xxy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 41);

    auto ta1_x_xxz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 42);

    auto ta1_x_xxz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 43);

    auto ta1_x_xxz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 44);

    auto ta1_x_xxz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 45);

    auto ta1_x_xxz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 46);

    auto ta1_x_xxz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 47);

    auto ta1_x_xxz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 48);

    auto ta1_x_xxz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 49);

    auto ta1_x_xxz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 50);

    auto ta1_x_xxz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 51);

    auto ta1_x_xxz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 52);

    auto ta1_x_xxz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 53);

    auto ta1_x_xxz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 54);

    auto ta1_x_xxz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 55);

    auto ta1_x_xxz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 56);

    auto ta1_x_xxz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 57);

    auto ta1_x_xxz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 58);

    auto ta1_x_xxz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 59);

    auto ta1_x_xxz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 60);

    auto ta1_x_xxz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 61);

    auto ta1_x_xxz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 62);

    auto ta1_x_xyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 63);

    auto ta1_x_xyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 64);

    auto ta1_x_xyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 65);

    auto ta1_x_xyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 66);

    auto ta1_x_xyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 67);

    auto ta1_x_xyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 68);

    auto ta1_x_xyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 69);

    auto ta1_x_xyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 70);

    auto ta1_x_xyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 71);

    auto ta1_x_xyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 72);

    auto ta1_x_xyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 73);

    auto ta1_x_xyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 74);

    auto ta1_x_xyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 75);

    auto ta1_x_xyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 76);

    auto ta1_x_xyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 77);

    auto ta1_x_xyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 78);

    auto ta1_x_xyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 79);

    auto ta1_x_xyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 80);

    auto ta1_x_xyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 81);

    auto ta1_x_xyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 82);

    auto ta1_x_xyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 86);

    auto ta1_x_xyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 89);

    auto ta1_x_xyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 93);

    auto ta1_x_xyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 98);

    auto ta1_x_xzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 105);

    auto ta1_x_xzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 106);

    auto ta1_x_xzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 107);

    auto ta1_x_xzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 108);

    auto ta1_x_xzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 109);

    auto ta1_x_xzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 110);

    auto ta1_x_xzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 111);

    auto ta1_x_xzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 112);

    auto ta1_x_xzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 113);

    auto ta1_x_xzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 114);

    auto ta1_x_xzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 115);

    auto ta1_x_xzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 116);

    auto ta1_x_xzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 117);

    auto ta1_x_xzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 118);

    auto ta1_x_xzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 119);

    auto ta1_x_xzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 121);

    auto ta1_x_xzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 122);

    auto ta1_x_xzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 123);

    auto ta1_x_xzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 124);

    auto ta1_x_xzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 125);

    auto ta1_x_yyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 126);

    auto ta1_x_yyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 127);

    auto ta1_x_yyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 128);

    auto ta1_x_yyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 129);

    auto ta1_x_yyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 130);

    auto ta1_x_yyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 131);

    auto ta1_x_yyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 132);

    auto ta1_x_yyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 133);

    auto ta1_x_yyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 134);

    auto ta1_x_yyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 135);

    auto ta1_x_yyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 136);

    auto ta1_x_yyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 137);

    auto ta1_x_yyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 138);

    auto ta1_x_yyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 139);

    auto ta1_x_yyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 140);

    auto ta1_x_yyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 141);

    auto ta1_x_yyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 142);

    auto ta1_x_yyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 143);

    auto ta1_x_yyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 144);

    auto ta1_x_yyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 145);

    auto ta1_x_yyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 146);

    auto ta1_x_yyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 148);

    auto ta1_x_yyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 149);

    auto ta1_x_yyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 150);

    auto ta1_x_yyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 152);

    auto ta1_x_yyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 153);

    auto ta1_x_yyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 156);

    auto ta1_x_yyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 157);

    auto ta1_x_yyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 161);

    auto ta1_x_yyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 162);

    auto ta1_x_yyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 163);

    auto ta1_x_yyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 164);

    auto ta1_x_yyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 165);

    auto ta1_x_yyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 166);

    auto ta1_x_yyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 167);

    auto ta1_x_yzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 168);

    auto ta1_x_yzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 170);

    auto ta1_x_yzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 172);

    auto ta1_x_yzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 173);

    auto ta1_x_yzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 175);

    auto ta1_x_yzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 176);

    auto ta1_x_yzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 177);

    auto ta1_x_yzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 179);

    auto ta1_x_yzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 180);

    auto ta1_x_yzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 181);

    auto ta1_x_yzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 182);

    auto ta1_x_yzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 183);

    auto ta1_x_yzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 184);

    auto ta1_x_yzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 185);

    auto ta1_x_yzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 186);

    auto ta1_x_yzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 187);

    auto ta1_x_yzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 188);

    auto ta1_x_zzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 189);

    auto ta1_x_zzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 190);

    auto ta1_x_zzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 191);

    auto ta1_x_zzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 192);

    auto ta1_x_zzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 193);

    auto ta1_x_zzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 194);

    auto ta1_x_zzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 195);

    auto ta1_x_zzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 196);

    auto ta1_x_zzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 197);

    auto ta1_x_zzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 198);

    auto ta1_x_zzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 199);

    auto ta1_x_zzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 200);

    auto ta1_x_zzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 201);

    auto ta1_x_zzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 202);

    auto ta1_x_zzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 203);

    auto ta1_x_zzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 204);

    auto ta1_x_zzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 205);

    auto ta1_x_zzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 206);

    auto ta1_x_zzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 207);

    auto ta1_x_zzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 208);

    auto ta1_x_zzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 209);

    auto ta1_y_xxx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 210);

    auto ta1_y_xxx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 211);

    auto ta1_y_xxx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 212);

    auto ta1_y_xxx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 213);

    auto ta1_y_xxx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 214);

    auto ta1_y_xxx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 215);

    auto ta1_y_xxx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 216);

    auto ta1_y_xxx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 217);

    auto ta1_y_xxx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 218);

    auto ta1_y_xxx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 219);

    auto ta1_y_xxx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 220);

    auto ta1_y_xxx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 221);

    auto ta1_y_xxx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 222);

    auto ta1_y_xxx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 223);

    auto ta1_y_xxx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 224);

    auto ta1_y_xxx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 225);

    auto ta1_y_xxx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 226);

    auto ta1_y_xxx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 227);

    auto ta1_y_xxx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 228);

    auto ta1_y_xxx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 229);

    auto ta1_y_xxx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 230);

    auto ta1_y_xxy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 231);

    auto ta1_y_xxy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 232);

    auto ta1_y_xxy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 233);

    auto ta1_y_xxy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 234);

    auto ta1_y_xxy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 235);

    auto ta1_y_xxy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 236);

    auto ta1_y_xxy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 237);

    auto ta1_y_xxy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 238);

    auto ta1_y_xxy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 239);

    auto ta1_y_xxy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 240);

    auto ta1_y_xxy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 241);

    auto ta1_y_xxy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 242);

    auto ta1_y_xxy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 243);

    auto ta1_y_xxy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 244);

    auto ta1_y_xxy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 245);

    auto ta1_y_xxy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 246);

    auto ta1_y_xxy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 247);

    auto ta1_y_xxy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 248);

    auto ta1_y_xxy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 249);

    auto ta1_y_xxy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 250);

    auto ta1_y_xxz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 252);

    auto ta1_y_xxz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 253);

    auto ta1_y_xxz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 254);

    auto ta1_y_xxz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 255);

    auto ta1_y_xxz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 257);

    auto ta1_y_xxz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 258);

    auto ta1_y_xxz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 261);

    auto ta1_y_xxz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 262);

    auto ta1_y_xxz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 266);

    auto ta1_y_xxz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 268);

    auto ta1_y_xxz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 269);

    auto ta1_y_xxz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 270);

    auto ta1_y_xxz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 271);

    auto ta1_y_xxz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 272);

    auto ta1_y_xyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 273);

    auto ta1_y_xyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 274);

    auto ta1_y_xyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 276);

    auto ta1_y_xyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 277);

    auto ta1_y_xyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 279);

    auto ta1_y_xyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 280);

    auto ta1_y_xyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 281);

    auto ta1_y_xyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 283);

    auto ta1_y_xyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 284);

    auto ta1_y_xyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 285);

    auto ta1_y_xyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 286);

    auto ta1_y_xyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 288);

    auto ta1_y_xyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 289);

    auto ta1_y_xyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 290);

    auto ta1_y_xyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 291);

    auto ta1_y_xyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 292);

    auto ta1_y_xyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 293);

    auto ta1_y_xyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 310);

    auto ta1_y_xyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 311);

    auto ta1_y_xyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 312);

    auto ta1_y_xyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 313);

    auto ta1_y_xzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 315);

    auto ta1_y_xzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 317);

    auto ta1_y_xzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 319);

    auto ta1_y_xzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 320);

    auto ta1_y_xzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 322);

    auto ta1_y_xzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 323);

    auto ta1_y_xzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 324);

    auto ta1_y_xzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 326);

    auto ta1_y_xzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 327);

    auto ta1_y_xzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 328);

    auto ta1_y_xzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 329);

    auto ta1_y_xzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 330);

    auto ta1_y_xzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 331);

    auto ta1_y_xzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 332);

    auto ta1_y_xzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 333);

    auto ta1_y_xzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 334);

    auto ta1_y_xzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 335);

    auto ta1_y_yyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 336);

    auto ta1_y_yyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 337);

    auto ta1_y_yyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 338);

    auto ta1_y_yyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 339);

    auto ta1_y_yyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 340);

    auto ta1_y_yyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 341);

    auto ta1_y_yyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 342);

    auto ta1_y_yyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 343);

    auto ta1_y_yyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 344);

    auto ta1_y_yyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 345);

    auto ta1_y_yyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 346);

    auto ta1_y_yyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 347);

    auto ta1_y_yyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 348);

    auto ta1_y_yyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 349);

    auto ta1_y_yyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 350);

    auto ta1_y_yyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 351);

    auto ta1_y_yyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 352);

    auto ta1_y_yyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 353);

    auto ta1_y_yyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 354);

    auto ta1_y_yyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 355);

    auto ta1_y_yyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 356);

    auto ta1_y_yyz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 357);

    auto ta1_y_yyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 358);

    auto ta1_y_yyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 359);

    auto ta1_y_yyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 360);

    auto ta1_y_yyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 361);

    auto ta1_y_yyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 362);

    auto ta1_y_yyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 363);

    auto ta1_y_yyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 364);

    auto ta1_y_yyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 365);

    auto ta1_y_yyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 366);

    auto ta1_y_yyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 367);

    auto ta1_y_yyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 368);

    auto ta1_y_yyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 369);

    auto ta1_y_yyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 370);

    auto ta1_y_yyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 371);

    auto ta1_y_yyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 372);

    auto ta1_y_yyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 373);

    auto ta1_y_yyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 374);

    auto ta1_y_yyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 375);

    auto ta1_y_yyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 376);

    auto ta1_y_yyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 377);

    auto ta1_y_yzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 379);

    auto ta1_y_yzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 380);

    auto ta1_y_yzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 381);

    auto ta1_y_yzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 382);

    auto ta1_y_yzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 383);

    auto ta1_y_yzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 384);

    auto ta1_y_yzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 385);

    auto ta1_y_yzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 386);

    auto ta1_y_yzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 387);

    auto ta1_y_yzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 388);

    auto ta1_y_yzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 389);

    auto ta1_y_yzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 390);

    auto ta1_y_yzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 391);

    auto ta1_y_yzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 392);

    auto ta1_y_yzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 393);

    auto ta1_y_yzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 394);

    auto ta1_y_yzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 395);

    auto ta1_y_yzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 396);

    auto ta1_y_yzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 397);

    auto ta1_y_yzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 398);

    auto ta1_y_zzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 399);

    auto ta1_y_zzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 400);

    auto ta1_y_zzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 401);

    auto ta1_y_zzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 402);

    auto ta1_y_zzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 403);

    auto ta1_y_zzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 404);

    auto ta1_y_zzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 405);

    auto ta1_y_zzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 406);

    auto ta1_y_zzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 407);

    auto ta1_y_zzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 408);

    auto ta1_y_zzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 409);

    auto ta1_y_zzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 410);

    auto ta1_y_zzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 411);

    auto ta1_y_zzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 412);

    auto ta1_y_zzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 413);

    auto ta1_y_zzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 414);

    auto ta1_y_zzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 415);

    auto ta1_y_zzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 416);

    auto ta1_y_zzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 417);

    auto ta1_y_zzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 418);

    auto ta1_y_zzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 419);

    auto ta1_z_xxx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 420);

    auto ta1_z_xxx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 421);

    auto ta1_z_xxx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 422);

    auto ta1_z_xxx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 423);

    auto ta1_z_xxx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 424);

    auto ta1_z_xxx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 425);

    auto ta1_z_xxx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 426);

    auto ta1_z_xxx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 427);

    auto ta1_z_xxx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 428);

    auto ta1_z_xxx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 429);

    auto ta1_z_xxx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 430);

    auto ta1_z_xxx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 431);

    auto ta1_z_xxx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 432);

    auto ta1_z_xxx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 433);

    auto ta1_z_xxx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 434);

    auto ta1_z_xxx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 435);

    auto ta1_z_xxx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 436);

    auto ta1_z_xxx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 437);

    auto ta1_z_xxx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 438);

    auto ta1_z_xxx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 439);

    auto ta1_z_xxx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 440);

    auto ta1_z_xxy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 441);

    auto ta1_z_xxy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 442);

    auto ta1_z_xxy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 443);

    auto ta1_z_xxy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 444);

    auto ta1_z_xxy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 446);

    auto ta1_z_xxy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 447);

    auto ta1_z_xxy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 450);

    auto ta1_z_xxy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 451);

    auto ta1_z_xxy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 455);

    auto ta1_z_xxy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 456);

    auto ta1_z_xxy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 457);

    auto ta1_z_xxy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 458);

    auto ta1_z_xxy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 459);

    auto ta1_z_xxy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 460);

    auto ta1_z_xxz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 462);

    auto ta1_z_xxz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 463);

    auto ta1_z_xxz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 464);

    auto ta1_z_xxz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 465);

    auto ta1_z_xxz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 466);

    auto ta1_z_xxz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 467);

    auto ta1_z_xxz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 468);

    auto ta1_z_xxz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 469);

    auto ta1_z_xxz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 470);

    auto ta1_z_xxz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 471);

    auto ta1_z_xxz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 472);

    auto ta1_z_xxz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 473);

    auto ta1_z_xxz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 474);

    auto ta1_z_xxz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 475);

    auto ta1_z_xxz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 476);

    auto ta1_z_xxz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 478);

    auto ta1_z_xxz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 479);

    auto ta1_z_xxz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 480);

    auto ta1_z_xxz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 481);

    auto ta1_z_xxz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 482);

    auto ta1_z_xyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 483);

    auto ta1_z_xyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 484);

    auto ta1_z_xyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 486);

    auto ta1_z_xyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 487);

    auto ta1_z_xyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 489);

    auto ta1_z_xyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 490);

    auto ta1_z_xyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 491);

    auto ta1_z_xyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 493);

    auto ta1_z_xyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 494);

    auto ta1_z_xyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 495);

    auto ta1_z_xyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 496);

    auto ta1_z_xyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 498);

    auto ta1_z_xyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 499);

    auto ta1_z_xyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 500);

    auto ta1_z_xyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 501);

    auto ta1_z_xyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 502);

    auto ta1_z_xyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 503);

    auto ta1_z_xyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 520);

    auto ta1_z_xyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 521);

    auto ta1_z_xyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 522);

    auto ta1_z_xyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 523);

    auto ta1_z_xzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 525);

    auto ta1_z_xzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 527);

    auto ta1_z_xzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 529);

    auto ta1_z_xzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 530);

    auto ta1_z_xzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 532);

    auto ta1_z_xzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 533);

    auto ta1_z_xzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 534);

    auto ta1_z_xzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 536);

    auto ta1_z_xzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 537);

    auto ta1_z_xzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 538);

    auto ta1_z_xzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 539);

    auto ta1_z_xzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 540);

    auto ta1_z_xzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 541);

    auto ta1_z_xzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 542);

    auto ta1_z_xzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 543);

    auto ta1_z_xzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 544);

    auto ta1_z_xzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 545);

    auto ta1_z_yyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 546);

    auto ta1_z_yyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 547);

    auto ta1_z_yyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 548);

    auto ta1_z_yyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 549);

    auto ta1_z_yyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 550);

    auto ta1_z_yyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 551);

    auto ta1_z_yyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 552);

    auto ta1_z_yyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 553);

    auto ta1_z_yyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 554);

    auto ta1_z_yyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 555);

    auto ta1_z_yyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 556);

    auto ta1_z_yyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 557);

    auto ta1_z_yyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 558);

    auto ta1_z_yyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 559);

    auto ta1_z_yyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 560);

    auto ta1_z_yyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 561);

    auto ta1_z_yyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 562);

    auto ta1_z_yyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 563);

    auto ta1_z_yyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 564);

    auto ta1_z_yyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 565);

    auto ta1_z_yyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 566);

    auto ta1_z_yyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 568);

    auto ta1_z_yyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 569);

    auto ta1_z_yyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 570);

    auto ta1_z_yyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 571);

    auto ta1_z_yyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 572);

    auto ta1_z_yyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 573);

    auto ta1_z_yyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 574);

    auto ta1_z_yyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 575);

    auto ta1_z_yyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 576);

    auto ta1_z_yyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 577);

    auto ta1_z_yyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 578);

    auto ta1_z_yyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 579);

    auto ta1_z_yyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 580);

    auto ta1_z_yyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 581);

    auto ta1_z_yyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 582);

    auto ta1_z_yyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 583);

    auto ta1_z_yyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 584);

    auto ta1_z_yyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 585);

    auto ta1_z_yyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 586);

    auto ta1_z_yyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 587);

    auto ta1_z_yzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 588);

    auto ta1_z_yzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 589);

    auto ta1_z_yzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 590);

    auto ta1_z_yzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 591);

    auto ta1_z_yzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 592);

    auto ta1_z_yzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 593);

    auto ta1_z_yzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 594);

    auto ta1_z_yzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 595);

    auto ta1_z_yzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 596);

    auto ta1_z_yzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 597);

    auto ta1_z_yzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 598);

    auto ta1_z_yzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 599);

    auto ta1_z_yzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 600);

    auto ta1_z_yzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 601);

    auto ta1_z_yzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 602);

    auto ta1_z_yzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 603);

    auto ta1_z_yzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 604);

    auto ta1_z_yzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 605);

    auto ta1_z_yzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 606);

    auto ta1_z_yzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 607);

    auto ta1_z_yzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 608);

    auto ta1_z_zzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_fh + 609);

    auto ta1_z_zzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 610);

    auto ta1_z_zzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 611);

    auto ta1_z_zzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 612);

    auto ta1_z_zzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 613);

    auto ta1_z_zzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 614);

    auto ta1_z_zzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 615);

    auto ta1_z_zzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 616);

    auto ta1_z_zzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 617);

    auto ta1_z_zzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 618);

    auto ta1_z_zzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 619);

    auto ta1_z_zzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 620);

    auto ta1_z_zzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 621);

    auto ta1_z_zzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 622);

    auto ta1_z_zzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 623);

    auto ta1_z_zzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_fh + 624);

    auto ta1_z_zzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 625);

    auto ta1_z_zzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 626);

    auto ta1_z_zzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 627);

    auto ta1_z_zzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 628);

    auto ta1_z_zzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_fh + 629);

    // Set up components of auxiliary buffer : FH

    auto ta1_x_xxx_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh);

    auto ta1_x_xxx_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 1);

    auto ta1_x_xxx_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 2);

    auto ta1_x_xxx_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 3);

    auto ta1_x_xxx_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 4);

    auto ta1_x_xxx_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 5);

    auto ta1_x_xxx_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 6);

    auto ta1_x_xxx_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 7);

    auto ta1_x_xxx_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 8);

    auto ta1_x_xxx_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 9);

    auto ta1_x_xxx_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 10);

    auto ta1_x_xxx_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 11);

    auto ta1_x_xxx_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 12);

    auto ta1_x_xxx_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 13);

    auto ta1_x_xxx_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 14);

    auto ta1_x_xxx_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 15);

    auto ta1_x_xxx_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 16);

    auto ta1_x_xxx_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 17);

    auto ta1_x_xxx_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 18);

    auto ta1_x_xxx_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 19);

    auto ta1_x_xxx_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 20);

    auto ta1_x_xxy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 21);

    auto ta1_x_xxy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 22);

    auto ta1_x_xxy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 23);

    auto ta1_x_xxy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 24);

    auto ta1_x_xxy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 25);

    auto ta1_x_xxy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 26);

    auto ta1_x_xxy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 27);

    auto ta1_x_xxy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 28);

    auto ta1_x_xxy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 29);

    auto ta1_x_xxy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 30);

    auto ta1_x_xxy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 31);

    auto ta1_x_xxy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 32);

    auto ta1_x_xxy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 33);

    auto ta1_x_xxy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 34);

    auto ta1_x_xxy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 35);

    auto ta1_x_xxy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 36);

    auto ta1_x_xxy_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 41);

    auto ta1_x_xxz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 42);

    auto ta1_x_xxz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 43);

    auto ta1_x_xxz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 44);

    auto ta1_x_xxz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 45);

    auto ta1_x_xxz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 46);

    auto ta1_x_xxz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 47);

    auto ta1_x_xxz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 48);

    auto ta1_x_xxz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 49);

    auto ta1_x_xxz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 50);

    auto ta1_x_xxz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 51);

    auto ta1_x_xxz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 52);

    auto ta1_x_xxz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 53);

    auto ta1_x_xxz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 54);

    auto ta1_x_xxz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 55);

    auto ta1_x_xxz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 56);

    auto ta1_x_xxz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 57);

    auto ta1_x_xxz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 58);

    auto ta1_x_xxz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 59);

    auto ta1_x_xxz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 60);

    auto ta1_x_xxz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 61);

    auto ta1_x_xxz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 62);

    auto ta1_x_xyy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 63);

    auto ta1_x_xyy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 64);

    auto ta1_x_xyy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 65);

    auto ta1_x_xyy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 66);

    auto ta1_x_xyy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 67);

    auto ta1_x_xyy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 68);

    auto ta1_x_xyy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 69);

    auto ta1_x_xyy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 70);

    auto ta1_x_xyy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 71);

    auto ta1_x_xyy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 72);

    auto ta1_x_xyy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 73);

    auto ta1_x_xyy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 74);

    auto ta1_x_xyy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 75);

    auto ta1_x_xyy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 76);

    auto ta1_x_xyy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 77);

    auto ta1_x_xyy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 78);

    auto ta1_x_xyy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 79);

    auto ta1_x_xyy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 80);

    auto ta1_x_xyy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 81);

    auto ta1_x_xyy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 82);

    auto ta1_x_xyz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 86);

    auto ta1_x_xyz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 89);

    auto ta1_x_xyz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 93);

    auto ta1_x_xyz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 98);

    auto ta1_x_xzz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 105);

    auto ta1_x_xzz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 106);

    auto ta1_x_xzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 107);

    auto ta1_x_xzz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 108);

    auto ta1_x_xzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 109);

    auto ta1_x_xzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 110);

    auto ta1_x_xzz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 111);

    auto ta1_x_xzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 112);

    auto ta1_x_xzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 113);

    auto ta1_x_xzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 114);

    auto ta1_x_xzz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 115);

    auto ta1_x_xzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 116);

    auto ta1_x_xzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 117);

    auto ta1_x_xzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 118);

    auto ta1_x_xzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 119);

    auto ta1_x_xzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 121);

    auto ta1_x_xzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 122);

    auto ta1_x_xzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 123);

    auto ta1_x_xzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 124);

    auto ta1_x_xzz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 125);

    auto ta1_x_yyy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 126);

    auto ta1_x_yyy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 127);

    auto ta1_x_yyy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 128);

    auto ta1_x_yyy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 129);

    auto ta1_x_yyy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 130);

    auto ta1_x_yyy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 131);

    auto ta1_x_yyy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 132);

    auto ta1_x_yyy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 133);

    auto ta1_x_yyy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 134);

    auto ta1_x_yyy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 135);

    auto ta1_x_yyy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 136);

    auto ta1_x_yyy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 137);

    auto ta1_x_yyy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 138);

    auto ta1_x_yyy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 139);

    auto ta1_x_yyy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 140);

    auto ta1_x_yyy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 141);

    auto ta1_x_yyy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 142);

    auto ta1_x_yyy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 143);

    auto ta1_x_yyy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 144);

    auto ta1_x_yyy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 145);

    auto ta1_x_yyy_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 146);

    auto ta1_x_yyz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 148);

    auto ta1_x_yyz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 149);

    auto ta1_x_yyz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 150);

    auto ta1_x_yyz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 152);

    auto ta1_x_yyz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 153);

    auto ta1_x_yyz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 156);

    auto ta1_x_yyz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 157);

    auto ta1_x_yyz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 161);

    auto ta1_x_yyz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 162);

    auto ta1_x_yyz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 163);

    auto ta1_x_yyz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 164);

    auto ta1_x_yyz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 165);

    auto ta1_x_yyz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 166);

    auto ta1_x_yyz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 167);

    auto ta1_x_yzz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 168);

    auto ta1_x_yzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 170);

    auto ta1_x_yzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 172);

    auto ta1_x_yzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 173);

    auto ta1_x_yzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 175);

    auto ta1_x_yzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 176);

    auto ta1_x_yzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 177);

    auto ta1_x_yzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 179);

    auto ta1_x_yzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 180);

    auto ta1_x_yzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 181);

    auto ta1_x_yzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 182);

    auto ta1_x_yzz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 183);

    auto ta1_x_yzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 184);

    auto ta1_x_yzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 185);

    auto ta1_x_yzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 186);

    auto ta1_x_yzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 187);

    auto ta1_x_yzz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 188);

    auto ta1_x_zzz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 189);

    auto ta1_x_zzz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 190);

    auto ta1_x_zzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 191);

    auto ta1_x_zzz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 192);

    auto ta1_x_zzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 193);

    auto ta1_x_zzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 194);

    auto ta1_x_zzz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 195);

    auto ta1_x_zzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 196);

    auto ta1_x_zzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 197);

    auto ta1_x_zzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 198);

    auto ta1_x_zzz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 199);

    auto ta1_x_zzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 200);

    auto ta1_x_zzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 201);

    auto ta1_x_zzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 202);

    auto ta1_x_zzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 203);

    auto ta1_x_zzz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 204);

    auto ta1_x_zzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 205);

    auto ta1_x_zzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 206);

    auto ta1_x_zzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 207);

    auto ta1_x_zzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 208);

    auto ta1_x_zzz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 209);

    auto ta1_y_xxx_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 210);

    auto ta1_y_xxx_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 211);

    auto ta1_y_xxx_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 212);

    auto ta1_y_xxx_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 213);

    auto ta1_y_xxx_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 214);

    auto ta1_y_xxx_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 215);

    auto ta1_y_xxx_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 216);

    auto ta1_y_xxx_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 217);

    auto ta1_y_xxx_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 218);

    auto ta1_y_xxx_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 219);

    auto ta1_y_xxx_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 220);

    auto ta1_y_xxx_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 221);

    auto ta1_y_xxx_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 222);

    auto ta1_y_xxx_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 223);

    auto ta1_y_xxx_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 224);

    auto ta1_y_xxx_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 225);

    auto ta1_y_xxx_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 226);

    auto ta1_y_xxx_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 227);

    auto ta1_y_xxx_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 228);

    auto ta1_y_xxx_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 229);

    auto ta1_y_xxx_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 230);

    auto ta1_y_xxy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 231);

    auto ta1_y_xxy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 232);

    auto ta1_y_xxy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 233);

    auto ta1_y_xxy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 234);

    auto ta1_y_xxy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 235);

    auto ta1_y_xxy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 236);

    auto ta1_y_xxy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 237);

    auto ta1_y_xxy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 238);

    auto ta1_y_xxy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 239);

    auto ta1_y_xxy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 240);

    auto ta1_y_xxy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 241);

    auto ta1_y_xxy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 242);

    auto ta1_y_xxy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 243);

    auto ta1_y_xxy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 244);

    auto ta1_y_xxy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 245);

    auto ta1_y_xxy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 246);

    auto ta1_y_xxy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 247);

    auto ta1_y_xxy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 248);

    auto ta1_y_xxy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 249);

    auto ta1_y_xxy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 250);

    auto ta1_y_xxz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 252);

    auto ta1_y_xxz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 253);

    auto ta1_y_xxz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 254);

    auto ta1_y_xxz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 255);

    auto ta1_y_xxz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 257);

    auto ta1_y_xxz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 258);

    auto ta1_y_xxz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 261);

    auto ta1_y_xxz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 262);

    auto ta1_y_xxz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 266);

    auto ta1_y_xxz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 268);

    auto ta1_y_xxz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 269);

    auto ta1_y_xxz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 270);

    auto ta1_y_xxz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 271);

    auto ta1_y_xxz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 272);

    auto ta1_y_xyy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 273);

    auto ta1_y_xyy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 274);

    auto ta1_y_xyy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 276);

    auto ta1_y_xyy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 277);

    auto ta1_y_xyy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 279);

    auto ta1_y_xyy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 280);

    auto ta1_y_xyy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 281);

    auto ta1_y_xyy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 283);

    auto ta1_y_xyy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 284);

    auto ta1_y_xyy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 285);

    auto ta1_y_xyy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 286);

    auto ta1_y_xyy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 288);

    auto ta1_y_xyy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 289);

    auto ta1_y_xyy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 290);

    auto ta1_y_xyy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 291);

    auto ta1_y_xyy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 292);

    auto ta1_y_xyy_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 293);

    auto ta1_y_xyz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 310);

    auto ta1_y_xyz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 311);

    auto ta1_y_xyz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 312);

    auto ta1_y_xyz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 313);

    auto ta1_y_xzz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 315);

    auto ta1_y_xzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 317);

    auto ta1_y_xzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 319);

    auto ta1_y_xzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 320);

    auto ta1_y_xzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 322);

    auto ta1_y_xzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 323);

    auto ta1_y_xzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 324);

    auto ta1_y_xzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 326);

    auto ta1_y_xzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 327);

    auto ta1_y_xzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 328);

    auto ta1_y_xzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 329);

    auto ta1_y_xzz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 330);

    auto ta1_y_xzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 331);

    auto ta1_y_xzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 332);

    auto ta1_y_xzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 333);

    auto ta1_y_xzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 334);

    auto ta1_y_xzz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 335);

    auto ta1_y_yyy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 336);

    auto ta1_y_yyy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 337);

    auto ta1_y_yyy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 338);

    auto ta1_y_yyy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 339);

    auto ta1_y_yyy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 340);

    auto ta1_y_yyy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 341);

    auto ta1_y_yyy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 342);

    auto ta1_y_yyy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 343);

    auto ta1_y_yyy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 344);

    auto ta1_y_yyy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 345);

    auto ta1_y_yyy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 346);

    auto ta1_y_yyy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 347);

    auto ta1_y_yyy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 348);

    auto ta1_y_yyy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 349);

    auto ta1_y_yyy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 350);

    auto ta1_y_yyy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 351);

    auto ta1_y_yyy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 352);

    auto ta1_y_yyy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 353);

    auto ta1_y_yyy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 354);

    auto ta1_y_yyy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 355);

    auto ta1_y_yyy_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 356);

    auto ta1_y_yyz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 357);

    auto ta1_y_yyz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 358);

    auto ta1_y_yyz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 359);

    auto ta1_y_yyz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 360);

    auto ta1_y_yyz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 361);

    auto ta1_y_yyz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 362);

    auto ta1_y_yyz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 363);

    auto ta1_y_yyz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 364);

    auto ta1_y_yyz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 365);

    auto ta1_y_yyz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 366);

    auto ta1_y_yyz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 367);

    auto ta1_y_yyz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 368);

    auto ta1_y_yyz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 369);

    auto ta1_y_yyz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 370);

    auto ta1_y_yyz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 371);

    auto ta1_y_yyz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 372);

    auto ta1_y_yyz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 373);

    auto ta1_y_yyz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 374);

    auto ta1_y_yyz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 375);

    auto ta1_y_yyz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 376);

    auto ta1_y_yyz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 377);

    auto ta1_y_yzz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 379);

    auto ta1_y_yzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 380);

    auto ta1_y_yzz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 381);

    auto ta1_y_yzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 382);

    auto ta1_y_yzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 383);

    auto ta1_y_yzz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 384);

    auto ta1_y_yzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 385);

    auto ta1_y_yzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 386);

    auto ta1_y_yzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 387);

    auto ta1_y_yzz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 388);

    auto ta1_y_yzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 389);

    auto ta1_y_yzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 390);

    auto ta1_y_yzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 391);

    auto ta1_y_yzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 392);

    auto ta1_y_yzz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 393);

    auto ta1_y_yzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 394);

    auto ta1_y_yzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 395);

    auto ta1_y_yzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 396);

    auto ta1_y_yzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 397);

    auto ta1_y_yzz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 398);

    auto ta1_y_zzz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 399);

    auto ta1_y_zzz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 400);

    auto ta1_y_zzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 401);

    auto ta1_y_zzz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 402);

    auto ta1_y_zzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 403);

    auto ta1_y_zzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 404);

    auto ta1_y_zzz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 405);

    auto ta1_y_zzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 406);

    auto ta1_y_zzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 407);

    auto ta1_y_zzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 408);

    auto ta1_y_zzz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 409);

    auto ta1_y_zzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 410);

    auto ta1_y_zzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 411);

    auto ta1_y_zzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 412);

    auto ta1_y_zzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 413);

    auto ta1_y_zzz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 414);

    auto ta1_y_zzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 415);

    auto ta1_y_zzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 416);

    auto ta1_y_zzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 417);

    auto ta1_y_zzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 418);

    auto ta1_y_zzz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 419);

    auto ta1_z_xxx_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 420);

    auto ta1_z_xxx_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 421);

    auto ta1_z_xxx_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 422);

    auto ta1_z_xxx_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 423);

    auto ta1_z_xxx_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 424);

    auto ta1_z_xxx_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 425);

    auto ta1_z_xxx_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 426);

    auto ta1_z_xxx_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 427);

    auto ta1_z_xxx_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 428);

    auto ta1_z_xxx_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 429);

    auto ta1_z_xxx_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 430);

    auto ta1_z_xxx_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 431);

    auto ta1_z_xxx_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 432);

    auto ta1_z_xxx_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 433);

    auto ta1_z_xxx_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 434);

    auto ta1_z_xxx_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 435);

    auto ta1_z_xxx_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 436);

    auto ta1_z_xxx_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 437);

    auto ta1_z_xxx_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 438);

    auto ta1_z_xxx_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 439);

    auto ta1_z_xxx_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 440);

    auto ta1_z_xxy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 441);

    auto ta1_z_xxy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 442);

    auto ta1_z_xxy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 443);

    auto ta1_z_xxy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 444);

    auto ta1_z_xxy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 446);

    auto ta1_z_xxy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 447);

    auto ta1_z_xxy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 450);

    auto ta1_z_xxy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 451);

    auto ta1_z_xxy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 455);

    auto ta1_z_xxy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 456);

    auto ta1_z_xxy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 457);

    auto ta1_z_xxy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 458);

    auto ta1_z_xxy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 459);

    auto ta1_z_xxy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 460);

    auto ta1_z_xxz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 462);

    auto ta1_z_xxz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 463);

    auto ta1_z_xxz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 464);

    auto ta1_z_xxz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 465);

    auto ta1_z_xxz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 466);

    auto ta1_z_xxz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 467);

    auto ta1_z_xxz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 468);

    auto ta1_z_xxz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 469);

    auto ta1_z_xxz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 470);

    auto ta1_z_xxz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 471);

    auto ta1_z_xxz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 472);

    auto ta1_z_xxz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 473);

    auto ta1_z_xxz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 474);

    auto ta1_z_xxz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 475);

    auto ta1_z_xxz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 476);

    auto ta1_z_xxz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 478);

    auto ta1_z_xxz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 479);

    auto ta1_z_xxz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 480);

    auto ta1_z_xxz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 481);

    auto ta1_z_xxz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 482);

    auto ta1_z_xyy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 483);

    auto ta1_z_xyy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 484);

    auto ta1_z_xyy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 486);

    auto ta1_z_xyy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 487);

    auto ta1_z_xyy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 489);

    auto ta1_z_xyy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 490);

    auto ta1_z_xyy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 491);

    auto ta1_z_xyy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 493);

    auto ta1_z_xyy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 494);

    auto ta1_z_xyy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 495);

    auto ta1_z_xyy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 496);

    auto ta1_z_xyy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 498);

    auto ta1_z_xyy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 499);

    auto ta1_z_xyy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 500);

    auto ta1_z_xyy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 501);

    auto ta1_z_xyy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 502);

    auto ta1_z_xyy_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 503);

    auto ta1_z_xyz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 520);

    auto ta1_z_xyz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 521);

    auto ta1_z_xyz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 522);

    auto ta1_z_xyz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 523);

    auto ta1_z_xzz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 525);

    auto ta1_z_xzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 527);

    auto ta1_z_xzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 529);

    auto ta1_z_xzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 530);

    auto ta1_z_xzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 532);

    auto ta1_z_xzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 533);

    auto ta1_z_xzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 534);

    auto ta1_z_xzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 536);

    auto ta1_z_xzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 537);

    auto ta1_z_xzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 538);

    auto ta1_z_xzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 539);

    auto ta1_z_xzz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 540);

    auto ta1_z_xzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 541);

    auto ta1_z_xzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 542);

    auto ta1_z_xzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 543);

    auto ta1_z_xzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 544);

    auto ta1_z_xzz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 545);

    auto ta1_z_yyy_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 546);

    auto ta1_z_yyy_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 547);

    auto ta1_z_yyy_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 548);

    auto ta1_z_yyy_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 549);

    auto ta1_z_yyy_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 550);

    auto ta1_z_yyy_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 551);

    auto ta1_z_yyy_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 552);

    auto ta1_z_yyy_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 553);

    auto ta1_z_yyy_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 554);

    auto ta1_z_yyy_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 555);

    auto ta1_z_yyy_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 556);

    auto ta1_z_yyy_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 557);

    auto ta1_z_yyy_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 558);

    auto ta1_z_yyy_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 559);

    auto ta1_z_yyy_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 560);

    auto ta1_z_yyy_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 561);

    auto ta1_z_yyy_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 562);

    auto ta1_z_yyy_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 563);

    auto ta1_z_yyy_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 564);

    auto ta1_z_yyy_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 565);

    auto ta1_z_yyy_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 566);

    auto ta1_z_yyz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 568);

    auto ta1_z_yyz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 569);

    auto ta1_z_yyz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 570);

    auto ta1_z_yyz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 571);

    auto ta1_z_yyz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 572);

    auto ta1_z_yyz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 573);

    auto ta1_z_yyz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 574);

    auto ta1_z_yyz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 575);

    auto ta1_z_yyz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 576);

    auto ta1_z_yyz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 577);

    auto ta1_z_yyz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 578);

    auto ta1_z_yyz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 579);

    auto ta1_z_yyz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 580);

    auto ta1_z_yyz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 581);

    auto ta1_z_yyz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 582);

    auto ta1_z_yyz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 583);

    auto ta1_z_yyz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 584);

    auto ta1_z_yyz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 585);

    auto ta1_z_yyz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 586);

    auto ta1_z_yyz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 587);

    auto ta1_z_yzz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 588);

    auto ta1_z_yzz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 589);

    auto ta1_z_yzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 590);

    auto ta1_z_yzz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 591);

    auto ta1_z_yzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 592);

    auto ta1_z_yzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 593);

    auto ta1_z_yzz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 594);

    auto ta1_z_yzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 595);

    auto ta1_z_yzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 596);

    auto ta1_z_yzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 597);

    auto ta1_z_yzz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 598);

    auto ta1_z_yzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 599);

    auto ta1_z_yzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 600);

    auto ta1_z_yzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 601);

    auto ta1_z_yzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 602);

    auto ta1_z_yzz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 603);

    auto ta1_z_yzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 604);

    auto ta1_z_yzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 605);

    auto ta1_z_yzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 606);

    auto ta1_z_yzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 607);

    auto ta1_z_yzz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 608);

    auto ta1_z_zzz_xxxxx_1 = pbuffer.data(idx_npot_geom_010_1_fh + 609);

    auto ta1_z_zzz_xxxxy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 610);

    auto ta1_z_zzz_xxxxz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 611);

    auto ta1_z_zzz_xxxyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 612);

    auto ta1_z_zzz_xxxyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 613);

    auto ta1_z_zzz_xxxzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 614);

    auto ta1_z_zzz_xxyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 615);

    auto ta1_z_zzz_xxyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 616);

    auto ta1_z_zzz_xxyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 617);

    auto ta1_z_zzz_xxzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 618);

    auto ta1_z_zzz_xyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 619);

    auto ta1_z_zzz_xyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 620);

    auto ta1_z_zzz_xyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 621);

    auto ta1_z_zzz_xyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 622);

    auto ta1_z_zzz_xzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 623);

    auto ta1_z_zzz_yyyyy_1 = pbuffer.data(idx_npot_geom_010_1_fh + 624);

    auto ta1_z_zzz_yyyyz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 625);

    auto ta1_z_zzz_yyyzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 626);

    auto ta1_z_zzz_yyzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 627);

    auto ta1_z_zzz_yzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 628);

    auto ta1_z_zzz_zzzzz_1 = pbuffer.data(idx_npot_geom_010_1_fh + 629);

    // Set up 0-21 components of targeted buffer : GH

    auto ta1_x_xxxx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh);

    auto ta1_x_xxxx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 1);

    auto ta1_x_xxxx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 2);

    auto ta1_x_xxxx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 3);

    auto ta1_x_xxxx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 4);

    auto ta1_x_xxxx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 5);

    auto ta1_x_xxxx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 6);

    auto ta1_x_xxxx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 7);

    auto ta1_x_xxxx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 8);

    auto ta1_x_xxxx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 9);

    auto ta1_x_xxxx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 10);

    auto ta1_x_xxxx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 11);

    auto ta1_x_xxxx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 12);

    auto ta1_x_xxxx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 13);

    auto ta1_x_xxxx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 14);

    auto ta1_x_xxxx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 15);

    auto ta1_x_xxxx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 16);

    auto ta1_x_xxxx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 17);

    auto ta1_x_xxxx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 18);

    auto ta1_x_xxxx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 19);

    auto ta1_x_xxxx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 20);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_x_xx_xxxxx_0,   \
                             ta1_x_xx_xxxxx_1,   \
                             ta1_x_xx_xxxxy_0,   \
                             ta1_x_xx_xxxxy_1,   \
                             ta1_x_xx_xxxxz_0,   \
                             ta1_x_xx_xxxxz_1,   \
                             ta1_x_xx_xxxyy_0,   \
                             ta1_x_xx_xxxyy_1,   \
                             ta1_x_xx_xxxyz_0,   \
                             ta1_x_xx_xxxyz_1,   \
                             ta1_x_xx_xxxzz_0,   \
                             ta1_x_xx_xxxzz_1,   \
                             ta1_x_xx_xxyyy_0,   \
                             ta1_x_xx_xxyyy_1,   \
                             ta1_x_xx_xxyyz_0,   \
                             ta1_x_xx_xxyyz_1,   \
                             ta1_x_xx_xxyzz_0,   \
                             ta1_x_xx_xxyzz_1,   \
                             ta1_x_xx_xxzzz_0,   \
                             ta1_x_xx_xxzzz_1,   \
                             ta1_x_xx_xyyyy_0,   \
                             ta1_x_xx_xyyyy_1,   \
                             ta1_x_xx_xyyyz_0,   \
                             ta1_x_xx_xyyyz_1,   \
                             ta1_x_xx_xyyzz_0,   \
                             ta1_x_xx_xyyzz_1,   \
                             ta1_x_xx_xyzzz_0,   \
                             ta1_x_xx_xyzzz_1,   \
                             ta1_x_xx_xzzzz_0,   \
                             ta1_x_xx_xzzzz_1,   \
                             ta1_x_xx_yyyyy_0,   \
                             ta1_x_xx_yyyyy_1,   \
                             ta1_x_xx_yyyyz_0,   \
                             ta1_x_xx_yyyyz_1,   \
                             ta1_x_xx_yyyzz_0,   \
                             ta1_x_xx_yyyzz_1,   \
                             ta1_x_xx_yyzzz_0,   \
                             ta1_x_xx_yyzzz_1,   \
                             ta1_x_xx_yzzzz_0,   \
                             ta1_x_xx_yzzzz_1,   \
                             ta1_x_xx_zzzzz_0,   \
                             ta1_x_xx_zzzzz_1,   \
                             ta1_x_xxx_xxxx_0,   \
                             ta1_x_xxx_xxxx_1,   \
                             ta1_x_xxx_xxxxx_0,  \
                             ta1_x_xxx_xxxxx_1,  \
                             ta1_x_xxx_xxxxy_0,  \
                             ta1_x_xxx_xxxxy_1,  \
                             ta1_x_xxx_xxxxz_0,  \
                             ta1_x_xxx_xxxxz_1,  \
                             ta1_x_xxx_xxxy_0,   \
                             ta1_x_xxx_xxxy_1,   \
                             ta1_x_xxx_xxxyy_0,  \
                             ta1_x_xxx_xxxyy_1,  \
                             ta1_x_xxx_xxxyz_0,  \
                             ta1_x_xxx_xxxyz_1,  \
                             ta1_x_xxx_xxxz_0,   \
                             ta1_x_xxx_xxxz_1,   \
                             ta1_x_xxx_xxxzz_0,  \
                             ta1_x_xxx_xxxzz_1,  \
                             ta1_x_xxx_xxyy_0,   \
                             ta1_x_xxx_xxyy_1,   \
                             ta1_x_xxx_xxyyy_0,  \
                             ta1_x_xxx_xxyyy_1,  \
                             ta1_x_xxx_xxyyz_0,  \
                             ta1_x_xxx_xxyyz_1,  \
                             ta1_x_xxx_xxyz_0,   \
                             ta1_x_xxx_xxyz_1,   \
                             ta1_x_xxx_xxyzz_0,  \
                             ta1_x_xxx_xxyzz_1,  \
                             ta1_x_xxx_xxzz_0,   \
                             ta1_x_xxx_xxzz_1,   \
                             ta1_x_xxx_xxzzz_0,  \
                             ta1_x_xxx_xxzzz_1,  \
                             ta1_x_xxx_xyyy_0,   \
                             ta1_x_xxx_xyyy_1,   \
                             ta1_x_xxx_xyyyy_0,  \
                             ta1_x_xxx_xyyyy_1,  \
                             ta1_x_xxx_xyyyz_0,  \
                             ta1_x_xxx_xyyyz_1,  \
                             ta1_x_xxx_xyyz_0,   \
                             ta1_x_xxx_xyyz_1,   \
                             ta1_x_xxx_xyyzz_0,  \
                             ta1_x_xxx_xyyzz_1,  \
                             ta1_x_xxx_xyzz_0,   \
                             ta1_x_xxx_xyzz_1,   \
                             ta1_x_xxx_xyzzz_0,  \
                             ta1_x_xxx_xyzzz_1,  \
                             ta1_x_xxx_xzzz_0,   \
                             ta1_x_xxx_xzzz_1,   \
                             ta1_x_xxx_xzzzz_0,  \
                             ta1_x_xxx_xzzzz_1,  \
                             ta1_x_xxx_yyyy_0,   \
                             ta1_x_xxx_yyyy_1,   \
                             ta1_x_xxx_yyyyy_0,  \
                             ta1_x_xxx_yyyyy_1,  \
                             ta1_x_xxx_yyyyz_0,  \
                             ta1_x_xxx_yyyyz_1,  \
                             ta1_x_xxx_yyyz_0,   \
                             ta1_x_xxx_yyyz_1,   \
                             ta1_x_xxx_yyyzz_0,  \
                             ta1_x_xxx_yyyzz_1,  \
                             ta1_x_xxx_yyzz_0,   \
                             ta1_x_xxx_yyzz_1,   \
                             ta1_x_xxx_yyzzz_0,  \
                             ta1_x_xxx_yyzzz_1,  \
                             ta1_x_xxx_yzzz_0,   \
                             ta1_x_xxx_yzzz_1,   \
                             ta1_x_xxx_yzzzz_0,  \
                             ta1_x_xxx_yzzzz_1,  \
                             ta1_x_xxx_zzzz_0,   \
                             ta1_x_xxx_zzzz_1,   \
                             ta1_x_xxx_zzzzz_0,  \
                             ta1_x_xxx_zzzzz_1,  \
                             ta1_x_xxxx_xxxxx_0, \
                             ta1_x_xxxx_xxxxy_0, \
                             ta1_x_xxxx_xxxxz_0, \
                             ta1_x_xxxx_xxxyy_0, \
                             ta1_x_xxxx_xxxyz_0, \
                             ta1_x_xxxx_xxxzz_0, \
                             ta1_x_xxxx_xxyyy_0, \
                             ta1_x_xxxx_xxyyz_0, \
                             ta1_x_xxxx_xxyzz_0, \
                             ta1_x_xxxx_xxzzz_0, \
                             ta1_x_xxxx_xyyyy_0, \
                             ta1_x_xxxx_xyyyz_0, \
                             ta1_x_xxxx_xyyzz_0, \
                             ta1_x_xxxx_xyzzz_0, \
                             ta1_x_xxxx_xzzzz_0, \
                             ta1_x_xxxx_yyyyy_0, \
                             ta1_x_xxxx_yyyyz_0, \
                             ta1_x_xxxx_yyyzz_0, \
                             ta1_x_xxxx_yyzzz_0, \
                             ta1_x_xxxx_yzzzz_0, \
                             ta1_x_xxxx_zzzzz_0, \
                             ta_xxx_xxxxx_1,     \
                             ta_xxx_xxxxy_1,     \
                             ta_xxx_xxxxz_1,     \
                             ta_xxx_xxxyy_1,     \
                             ta_xxx_xxxyz_1,     \
                             ta_xxx_xxxzz_1,     \
                             ta_xxx_xxyyy_1,     \
                             ta_xxx_xxyyz_1,     \
                             ta_xxx_xxyzz_1,     \
                             ta_xxx_xxzzz_1,     \
                             ta_xxx_xyyyy_1,     \
                             ta_xxx_xyyyz_1,     \
                             ta_xxx_xyyzz_1,     \
                             ta_xxx_xyzzz_1,     \
                             ta_xxx_xzzzz_1,     \
                             ta_xxx_yyyyy_1,     \
                             ta_xxx_yyyyz_1,     \
                             ta_xxx_yyyzz_1,     \
                             ta_xxx_yyzzz_1,     \
                             ta_xxx_yzzzz_1,     \
                             ta_xxx_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxx_xxxxx_0[i] = 3.0 * ta1_x_xx_xxxxx_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxxx_1[i] * fe_0 + 5.0 * ta1_x_xxx_xxxx_0[i] * fe_0 -
                                5.0 * ta1_x_xxx_xxxx_1[i] * fe_0 + ta_xxx_xxxxx_1[i] + ta1_x_xxx_xxxxx_0[i] * pa_x[i] -
                                ta1_x_xxx_xxxxx_1[i] * pc_x[i];

        ta1_x_xxxx_xxxxy_0[i] = 3.0 * ta1_x_xx_xxxxy_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxxy_1[i] * fe_0 + 4.0 * ta1_x_xxx_xxxy_0[i] * fe_0 -
                                4.0 * ta1_x_xxx_xxxy_1[i] * fe_0 + ta_xxx_xxxxy_1[i] + ta1_x_xxx_xxxxy_0[i] * pa_x[i] -
                                ta1_x_xxx_xxxxy_1[i] * pc_x[i];

        ta1_x_xxxx_xxxxz_0[i] = 3.0 * ta1_x_xx_xxxxz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxxz_1[i] * fe_0 + 4.0 * ta1_x_xxx_xxxz_0[i] * fe_0 -
                                4.0 * ta1_x_xxx_xxxz_1[i] * fe_0 + ta_xxx_xxxxz_1[i] + ta1_x_xxx_xxxxz_0[i] * pa_x[i] -
                                ta1_x_xxx_xxxxz_1[i] * pc_x[i];

        ta1_x_xxxx_xxxyy_0[i] = 3.0 * ta1_x_xx_xxxyy_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxyy_1[i] * fe_0 + 3.0 * ta1_x_xxx_xxyy_0[i] * fe_0 -
                                3.0 * ta1_x_xxx_xxyy_1[i] * fe_0 + ta_xxx_xxxyy_1[i] + ta1_x_xxx_xxxyy_0[i] * pa_x[i] -
                                ta1_x_xxx_xxxyy_1[i] * pc_x[i];

        ta1_x_xxxx_xxxyz_0[i] = 3.0 * ta1_x_xx_xxxyz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxyz_1[i] * fe_0 + 3.0 * ta1_x_xxx_xxyz_0[i] * fe_0 -
                                3.0 * ta1_x_xxx_xxyz_1[i] * fe_0 + ta_xxx_xxxyz_1[i] + ta1_x_xxx_xxxyz_0[i] * pa_x[i] -
                                ta1_x_xxx_xxxyz_1[i] * pc_x[i];

        ta1_x_xxxx_xxxzz_0[i] = 3.0 * ta1_x_xx_xxxzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxxzz_1[i] * fe_0 + 3.0 * ta1_x_xxx_xxzz_0[i] * fe_0 -
                                3.0 * ta1_x_xxx_xxzz_1[i] * fe_0 + ta_xxx_xxxzz_1[i] + ta1_x_xxx_xxxzz_0[i] * pa_x[i] -
                                ta1_x_xxx_xxxzz_1[i] * pc_x[i];

        ta1_x_xxxx_xxyyy_0[i] = 3.0 * ta1_x_xx_xxyyy_0[i] * fe_0 - 3.0 * ta1_x_xx_xxyyy_1[i] * fe_0 + 2.0 * ta1_x_xxx_xyyy_0[i] * fe_0 -
                                2.0 * ta1_x_xxx_xyyy_1[i] * fe_0 + ta_xxx_xxyyy_1[i] + ta1_x_xxx_xxyyy_0[i] * pa_x[i] -
                                ta1_x_xxx_xxyyy_1[i] * pc_x[i];

        ta1_x_xxxx_xxyyz_0[i] = 3.0 * ta1_x_xx_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxyyz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xyyz_0[i] * fe_0 -
                                2.0 * ta1_x_xxx_xyyz_1[i] * fe_0 + ta_xxx_xxyyz_1[i] + ta1_x_xxx_xxyyz_0[i] * pa_x[i] -
                                ta1_x_xxx_xxyyz_1[i] * pc_x[i];

        ta1_x_xxxx_xxyzz_0[i] = 3.0 * ta1_x_xx_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxyzz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xyzz_0[i] * fe_0 -
                                2.0 * ta1_x_xxx_xyzz_1[i] * fe_0 + ta_xxx_xxyzz_1[i] + ta1_x_xxx_xxyzz_0[i] * pa_x[i] -
                                ta1_x_xxx_xxyzz_1[i] * pc_x[i];

        ta1_x_xxxx_xxzzz_0[i] = 3.0 * ta1_x_xx_xxzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xxzzz_1[i] * fe_0 + 2.0 * ta1_x_xxx_xzzz_0[i] * fe_0 -
                                2.0 * ta1_x_xxx_xzzz_1[i] * fe_0 + ta_xxx_xxzzz_1[i] + ta1_x_xxx_xxzzz_0[i] * pa_x[i] -
                                ta1_x_xxx_xxzzz_1[i] * pc_x[i];

        ta1_x_xxxx_xyyyy_0[i] = 3.0 * ta1_x_xx_xyyyy_0[i] * fe_0 - 3.0 * ta1_x_xx_xyyyy_1[i] * fe_0 + ta1_x_xxx_yyyy_0[i] * fe_0 -
                                ta1_x_xxx_yyyy_1[i] * fe_0 + ta_xxx_xyyyy_1[i] + ta1_x_xxx_xyyyy_0[i] * pa_x[i] - ta1_x_xxx_xyyyy_1[i] * pc_x[i];

        ta1_x_xxxx_xyyyz_0[i] = 3.0 * ta1_x_xx_xyyyz_0[i] * fe_0 - 3.0 * ta1_x_xx_xyyyz_1[i] * fe_0 + ta1_x_xxx_yyyz_0[i] * fe_0 -
                                ta1_x_xxx_yyyz_1[i] * fe_0 + ta_xxx_xyyyz_1[i] + ta1_x_xxx_xyyyz_0[i] * pa_x[i] - ta1_x_xxx_xyyyz_1[i] * pc_x[i];

        ta1_x_xxxx_xyyzz_0[i] = 3.0 * ta1_x_xx_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xyyzz_1[i] * fe_0 + ta1_x_xxx_yyzz_0[i] * fe_0 -
                                ta1_x_xxx_yyzz_1[i] * fe_0 + ta_xxx_xyyzz_1[i] + ta1_x_xxx_xyyzz_0[i] * pa_x[i] - ta1_x_xxx_xyyzz_1[i] * pc_x[i];

        ta1_x_xxxx_xyzzz_0[i] = 3.0 * ta1_x_xx_xyzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xyzzz_1[i] * fe_0 + ta1_x_xxx_yzzz_0[i] * fe_0 -
                                ta1_x_xxx_yzzz_1[i] * fe_0 + ta_xxx_xyzzz_1[i] + ta1_x_xxx_xyzzz_0[i] * pa_x[i] - ta1_x_xxx_xyzzz_1[i] * pc_x[i];

        ta1_x_xxxx_xzzzz_0[i] = 3.0 * ta1_x_xx_xzzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_xzzzz_1[i] * fe_0 + ta1_x_xxx_zzzz_0[i] * fe_0 -
                                ta1_x_xxx_zzzz_1[i] * fe_0 + ta_xxx_xzzzz_1[i] + ta1_x_xxx_xzzzz_0[i] * pa_x[i] - ta1_x_xxx_xzzzz_1[i] * pc_x[i];

        ta1_x_xxxx_yyyyy_0[i] = 3.0 * ta1_x_xx_yyyyy_0[i] * fe_0 - 3.0 * ta1_x_xx_yyyyy_1[i] * fe_0 + ta_xxx_yyyyy_1[i] +
                                ta1_x_xxx_yyyyy_0[i] * pa_x[i] - ta1_x_xxx_yyyyy_1[i] * pc_x[i];

        ta1_x_xxxx_yyyyz_0[i] = 3.0 * ta1_x_xx_yyyyz_0[i] * fe_0 - 3.0 * ta1_x_xx_yyyyz_1[i] * fe_0 + ta_xxx_yyyyz_1[i] +
                                ta1_x_xxx_yyyyz_0[i] * pa_x[i] - ta1_x_xxx_yyyyz_1[i] * pc_x[i];

        ta1_x_xxxx_yyyzz_0[i] = 3.0 * ta1_x_xx_yyyzz_0[i] * fe_0 - 3.0 * ta1_x_xx_yyyzz_1[i] * fe_0 + ta_xxx_yyyzz_1[i] +
                                ta1_x_xxx_yyyzz_0[i] * pa_x[i] - ta1_x_xxx_yyyzz_1[i] * pc_x[i];

        ta1_x_xxxx_yyzzz_0[i] = 3.0 * ta1_x_xx_yyzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_yyzzz_1[i] * fe_0 + ta_xxx_yyzzz_1[i] +
                                ta1_x_xxx_yyzzz_0[i] * pa_x[i] - ta1_x_xxx_yyzzz_1[i] * pc_x[i];

        ta1_x_xxxx_yzzzz_0[i] = 3.0 * ta1_x_xx_yzzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_yzzzz_1[i] * fe_0 + ta_xxx_yzzzz_1[i] +
                                ta1_x_xxx_yzzzz_0[i] * pa_x[i] - ta1_x_xxx_yzzzz_1[i] * pc_x[i];

        ta1_x_xxxx_zzzzz_0[i] = 3.0 * ta1_x_xx_zzzzz_0[i] * fe_0 - 3.0 * ta1_x_xx_zzzzz_1[i] * fe_0 + ta_xxx_zzzzz_1[i] +
                                ta1_x_xxx_zzzzz_0[i] * pa_x[i] - ta1_x_xxx_zzzzz_1[i] * pc_x[i];
    }

    // Set up 21-42 components of targeted buffer : GH

    auto ta1_x_xxxy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 21);

    auto ta1_x_xxxy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 22);

    auto ta1_x_xxxy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 23);

    auto ta1_x_xxxy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 24);

    auto ta1_x_xxxy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 25);

    auto ta1_x_xxxy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 26);

    auto ta1_x_xxxy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 27);

    auto ta1_x_xxxy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 28);

    auto ta1_x_xxxy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 29);

    auto ta1_x_xxxy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 30);

    auto ta1_x_xxxy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 31);

    auto ta1_x_xxxy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 32);

    auto ta1_x_xxxy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 33);

    auto ta1_x_xxxy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 34);

    auto ta1_x_xxxy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 35);

    auto ta1_x_xxxy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 36);

    auto ta1_x_xxxy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 37);

    auto ta1_x_xxxy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 38);

    auto ta1_x_xxxy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 39);

    auto ta1_x_xxxy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 40);

    auto ta1_x_xxxy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 41);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_x_xxx_xxxx_0,   \
                             ta1_x_xxx_xxxx_1,   \
                             ta1_x_xxx_xxxxx_0,  \
                             ta1_x_xxx_xxxxx_1,  \
                             ta1_x_xxx_xxxxy_0,  \
                             ta1_x_xxx_xxxxy_1,  \
                             ta1_x_xxx_xxxxz_0,  \
                             ta1_x_xxx_xxxxz_1,  \
                             ta1_x_xxx_xxxy_0,   \
                             ta1_x_xxx_xxxy_1,   \
                             ta1_x_xxx_xxxyy_0,  \
                             ta1_x_xxx_xxxyy_1,  \
                             ta1_x_xxx_xxxyz_0,  \
                             ta1_x_xxx_xxxyz_1,  \
                             ta1_x_xxx_xxxz_0,   \
                             ta1_x_xxx_xxxz_1,   \
                             ta1_x_xxx_xxxzz_0,  \
                             ta1_x_xxx_xxxzz_1,  \
                             ta1_x_xxx_xxyy_0,   \
                             ta1_x_xxx_xxyy_1,   \
                             ta1_x_xxx_xxyyy_0,  \
                             ta1_x_xxx_xxyyy_1,  \
                             ta1_x_xxx_xxyyz_0,  \
                             ta1_x_xxx_xxyyz_1,  \
                             ta1_x_xxx_xxyz_0,   \
                             ta1_x_xxx_xxyz_1,   \
                             ta1_x_xxx_xxyzz_0,  \
                             ta1_x_xxx_xxyzz_1,  \
                             ta1_x_xxx_xxzz_0,   \
                             ta1_x_xxx_xxzz_1,   \
                             ta1_x_xxx_xxzzz_0,  \
                             ta1_x_xxx_xxzzz_1,  \
                             ta1_x_xxx_xyyy_0,   \
                             ta1_x_xxx_xyyy_1,   \
                             ta1_x_xxx_xyyyy_0,  \
                             ta1_x_xxx_xyyyy_1,  \
                             ta1_x_xxx_xyyyz_0,  \
                             ta1_x_xxx_xyyyz_1,  \
                             ta1_x_xxx_xyyz_0,   \
                             ta1_x_xxx_xyyz_1,   \
                             ta1_x_xxx_xyyzz_0,  \
                             ta1_x_xxx_xyyzz_1,  \
                             ta1_x_xxx_xyzz_0,   \
                             ta1_x_xxx_xyzz_1,   \
                             ta1_x_xxx_xyzzz_0,  \
                             ta1_x_xxx_xyzzz_1,  \
                             ta1_x_xxx_xzzz_0,   \
                             ta1_x_xxx_xzzz_1,   \
                             ta1_x_xxx_xzzzz_0,  \
                             ta1_x_xxx_xzzzz_1,  \
                             ta1_x_xxx_yyyy_0,   \
                             ta1_x_xxx_yyyy_1,   \
                             ta1_x_xxx_yyyyy_0,  \
                             ta1_x_xxx_yyyyy_1,  \
                             ta1_x_xxx_yyyyz_0,  \
                             ta1_x_xxx_yyyyz_1,  \
                             ta1_x_xxx_yyyz_0,   \
                             ta1_x_xxx_yyyz_1,   \
                             ta1_x_xxx_yyyzz_0,  \
                             ta1_x_xxx_yyyzz_1,  \
                             ta1_x_xxx_yyzz_0,   \
                             ta1_x_xxx_yyzz_1,   \
                             ta1_x_xxx_yyzzz_0,  \
                             ta1_x_xxx_yyzzz_1,  \
                             ta1_x_xxx_yzzz_0,   \
                             ta1_x_xxx_yzzz_1,   \
                             ta1_x_xxx_yzzzz_0,  \
                             ta1_x_xxx_yzzzz_1,  \
                             ta1_x_xxx_zzzz_0,   \
                             ta1_x_xxx_zzzz_1,   \
                             ta1_x_xxx_zzzzz_0,  \
                             ta1_x_xxx_zzzzz_1,  \
                             ta1_x_xxxy_xxxxx_0, \
                             ta1_x_xxxy_xxxxy_0, \
                             ta1_x_xxxy_xxxxz_0, \
                             ta1_x_xxxy_xxxyy_0, \
                             ta1_x_xxxy_xxxyz_0, \
                             ta1_x_xxxy_xxxzz_0, \
                             ta1_x_xxxy_xxyyy_0, \
                             ta1_x_xxxy_xxyyz_0, \
                             ta1_x_xxxy_xxyzz_0, \
                             ta1_x_xxxy_xxzzz_0, \
                             ta1_x_xxxy_xyyyy_0, \
                             ta1_x_xxxy_xyyyz_0, \
                             ta1_x_xxxy_xyyzz_0, \
                             ta1_x_xxxy_xyzzz_0, \
                             ta1_x_xxxy_xzzzz_0, \
                             ta1_x_xxxy_yyyyy_0, \
                             ta1_x_xxxy_yyyyz_0, \
                             ta1_x_xxxy_yyyzz_0, \
                             ta1_x_xxxy_yyzzz_0, \
                             ta1_x_xxxy_yzzzz_0, \
                             ta1_x_xxxy_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxy_xxxxx_0[i] = ta1_x_xxx_xxxxx_0[i] * pa_y[i] - ta1_x_xxx_xxxxx_1[i] * pc_y[i];

        ta1_x_xxxy_xxxxy_0[i] =
            ta1_x_xxx_xxxx_0[i] * fe_0 - ta1_x_xxx_xxxx_1[i] * fe_0 + ta1_x_xxx_xxxxy_0[i] * pa_y[i] - ta1_x_xxx_xxxxy_1[i] * pc_y[i];

        ta1_x_xxxy_xxxxz_0[i] = ta1_x_xxx_xxxxz_0[i] * pa_y[i] - ta1_x_xxx_xxxxz_1[i] * pc_y[i];

        ta1_x_xxxy_xxxyy_0[i] =
            2.0 * ta1_x_xxx_xxxy_0[i] * fe_0 - 2.0 * ta1_x_xxx_xxxy_1[i] * fe_0 + ta1_x_xxx_xxxyy_0[i] * pa_y[i] - ta1_x_xxx_xxxyy_1[i] * pc_y[i];

        ta1_x_xxxy_xxxyz_0[i] =
            ta1_x_xxx_xxxz_0[i] * fe_0 - ta1_x_xxx_xxxz_1[i] * fe_0 + ta1_x_xxx_xxxyz_0[i] * pa_y[i] - ta1_x_xxx_xxxyz_1[i] * pc_y[i];

        ta1_x_xxxy_xxxzz_0[i] = ta1_x_xxx_xxxzz_0[i] * pa_y[i] - ta1_x_xxx_xxxzz_1[i] * pc_y[i];

        ta1_x_xxxy_xxyyy_0[i] =
            3.0 * ta1_x_xxx_xxyy_0[i] * fe_0 - 3.0 * ta1_x_xxx_xxyy_1[i] * fe_0 + ta1_x_xxx_xxyyy_0[i] * pa_y[i] - ta1_x_xxx_xxyyy_1[i] * pc_y[i];

        ta1_x_xxxy_xxyyz_0[i] =
            2.0 * ta1_x_xxx_xxyz_0[i] * fe_0 - 2.0 * ta1_x_xxx_xxyz_1[i] * fe_0 + ta1_x_xxx_xxyyz_0[i] * pa_y[i] - ta1_x_xxx_xxyyz_1[i] * pc_y[i];

        ta1_x_xxxy_xxyzz_0[i] =
            ta1_x_xxx_xxzz_0[i] * fe_0 - ta1_x_xxx_xxzz_1[i] * fe_0 + ta1_x_xxx_xxyzz_0[i] * pa_y[i] - ta1_x_xxx_xxyzz_1[i] * pc_y[i];

        ta1_x_xxxy_xxzzz_0[i] = ta1_x_xxx_xxzzz_0[i] * pa_y[i] - ta1_x_xxx_xxzzz_1[i] * pc_y[i];

        ta1_x_xxxy_xyyyy_0[i] =
            4.0 * ta1_x_xxx_xyyy_0[i] * fe_0 - 4.0 * ta1_x_xxx_xyyy_1[i] * fe_0 + ta1_x_xxx_xyyyy_0[i] * pa_y[i] - ta1_x_xxx_xyyyy_1[i] * pc_y[i];

        ta1_x_xxxy_xyyyz_0[i] =
            3.0 * ta1_x_xxx_xyyz_0[i] * fe_0 - 3.0 * ta1_x_xxx_xyyz_1[i] * fe_0 + ta1_x_xxx_xyyyz_0[i] * pa_y[i] - ta1_x_xxx_xyyyz_1[i] * pc_y[i];

        ta1_x_xxxy_xyyzz_0[i] =
            2.0 * ta1_x_xxx_xyzz_0[i] * fe_0 - 2.0 * ta1_x_xxx_xyzz_1[i] * fe_0 + ta1_x_xxx_xyyzz_0[i] * pa_y[i] - ta1_x_xxx_xyyzz_1[i] * pc_y[i];

        ta1_x_xxxy_xyzzz_0[i] =
            ta1_x_xxx_xzzz_0[i] * fe_0 - ta1_x_xxx_xzzz_1[i] * fe_0 + ta1_x_xxx_xyzzz_0[i] * pa_y[i] - ta1_x_xxx_xyzzz_1[i] * pc_y[i];

        ta1_x_xxxy_xzzzz_0[i] = ta1_x_xxx_xzzzz_0[i] * pa_y[i] - ta1_x_xxx_xzzzz_1[i] * pc_y[i];

        ta1_x_xxxy_yyyyy_0[i] =
            5.0 * ta1_x_xxx_yyyy_0[i] * fe_0 - 5.0 * ta1_x_xxx_yyyy_1[i] * fe_0 + ta1_x_xxx_yyyyy_0[i] * pa_y[i] - ta1_x_xxx_yyyyy_1[i] * pc_y[i];

        ta1_x_xxxy_yyyyz_0[i] =
            4.0 * ta1_x_xxx_yyyz_0[i] * fe_0 - 4.0 * ta1_x_xxx_yyyz_1[i] * fe_0 + ta1_x_xxx_yyyyz_0[i] * pa_y[i] - ta1_x_xxx_yyyyz_1[i] * pc_y[i];

        ta1_x_xxxy_yyyzz_0[i] =
            3.0 * ta1_x_xxx_yyzz_0[i] * fe_0 - 3.0 * ta1_x_xxx_yyzz_1[i] * fe_0 + ta1_x_xxx_yyyzz_0[i] * pa_y[i] - ta1_x_xxx_yyyzz_1[i] * pc_y[i];

        ta1_x_xxxy_yyzzz_0[i] =
            2.0 * ta1_x_xxx_yzzz_0[i] * fe_0 - 2.0 * ta1_x_xxx_yzzz_1[i] * fe_0 + ta1_x_xxx_yyzzz_0[i] * pa_y[i] - ta1_x_xxx_yyzzz_1[i] * pc_y[i];

        ta1_x_xxxy_yzzzz_0[i] =
            ta1_x_xxx_zzzz_0[i] * fe_0 - ta1_x_xxx_zzzz_1[i] * fe_0 + ta1_x_xxx_yzzzz_0[i] * pa_y[i] - ta1_x_xxx_yzzzz_1[i] * pc_y[i];

        ta1_x_xxxy_zzzzz_0[i] = ta1_x_xxx_zzzzz_0[i] * pa_y[i] - ta1_x_xxx_zzzzz_1[i] * pc_y[i];
    }

    // Set up 42-63 components of targeted buffer : GH

    auto ta1_x_xxxz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 42);

    auto ta1_x_xxxz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 43);

    auto ta1_x_xxxz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 44);

    auto ta1_x_xxxz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 45);

    auto ta1_x_xxxz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 46);

    auto ta1_x_xxxz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 47);

    auto ta1_x_xxxz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 48);

    auto ta1_x_xxxz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 49);

    auto ta1_x_xxxz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 50);

    auto ta1_x_xxxz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 51);

    auto ta1_x_xxxz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 52);

    auto ta1_x_xxxz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 53);

    auto ta1_x_xxxz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 54);

    auto ta1_x_xxxz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 55);

    auto ta1_x_xxxz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 56);

    auto ta1_x_xxxz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 57);

    auto ta1_x_xxxz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 58);

    auto ta1_x_xxxz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 59);

    auto ta1_x_xxxz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 60);

    auto ta1_x_xxxz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 61);

    auto ta1_x_xxxz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 62);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_x_xxx_xxxx_0,   \
                             ta1_x_xxx_xxxx_1,   \
                             ta1_x_xxx_xxxxx_0,  \
                             ta1_x_xxx_xxxxx_1,  \
                             ta1_x_xxx_xxxxy_0,  \
                             ta1_x_xxx_xxxxy_1,  \
                             ta1_x_xxx_xxxxz_0,  \
                             ta1_x_xxx_xxxxz_1,  \
                             ta1_x_xxx_xxxy_0,   \
                             ta1_x_xxx_xxxy_1,   \
                             ta1_x_xxx_xxxyy_0,  \
                             ta1_x_xxx_xxxyy_1,  \
                             ta1_x_xxx_xxxyz_0,  \
                             ta1_x_xxx_xxxyz_1,  \
                             ta1_x_xxx_xxxz_0,   \
                             ta1_x_xxx_xxxz_1,   \
                             ta1_x_xxx_xxxzz_0,  \
                             ta1_x_xxx_xxxzz_1,  \
                             ta1_x_xxx_xxyy_0,   \
                             ta1_x_xxx_xxyy_1,   \
                             ta1_x_xxx_xxyyy_0,  \
                             ta1_x_xxx_xxyyy_1,  \
                             ta1_x_xxx_xxyyz_0,  \
                             ta1_x_xxx_xxyyz_1,  \
                             ta1_x_xxx_xxyz_0,   \
                             ta1_x_xxx_xxyz_1,   \
                             ta1_x_xxx_xxyzz_0,  \
                             ta1_x_xxx_xxyzz_1,  \
                             ta1_x_xxx_xxzz_0,   \
                             ta1_x_xxx_xxzz_1,   \
                             ta1_x_xxx_xxzzz_0,  \
                             ta1_x_xxx_xxzzz_1,  \
                             ta1_x_xxx_xyyy_0,   \
                             ta1_x_xxx_xyyy_1,   \
                             ta1_x_xxx_xyyyy_0,  \
                             ta1_x_xxx_xyyyy_1,  \
                             ta1_x_xxx_xyyyz_0,  \
                             ta1_x_xxx_xyyyz_1,  \
                             ta1_x_xxx_xyyz_0,   \
                             ta1_x_xxx_xyyz_1,   \
                             ta1_x_xxx_xyyzz_0,  \
                             ta1_x_xxx_xyyzz_1,  \
                             ta1_x_xxx_xyzz_0,   \
                             ta1_x_xxx_xyzz_1,   \
                             ta1_x_xxx_xyzzz_0,  \
                             ta1_x_xxx_xyzzz_1,  \
                             ta1_x_xxx_xzzz_0,   \
                             ta1_x_xxx_xzzz_1,   \
                             ta1_x_xxx_xzzzz_0,  \
                             ta1_x_xxx_xzzzz_1,  \
                             ta1_x_xxx_yyyy_0,   \
                             ta1_x_xxx_yyyy_1,   \
                             ta1_x_xxx_yyyyy_0,  \
                             ta1_x_xxx_yyyyy_1,  \
                             ta1_x_xxx_yyyyz_0,  \
                             ta1_x_xxx_yyyyz_1,  \
                             ta1_x_xxx_yyyz_0,   \
                             ta1_x_xxx_yyyz_1,   \
                             ta1_x_xxx_yyyzz_0,  \
                             ta1_x_xxx_yyyzz_1,  \
                             ta1_x_xxx_yyzz_0,   \
                             ta1_x_xxx_yyzz_1,   \
                             ta1_x_xxx_yyzzz_0,  \
                             ta1_x_xxx_yyzzz_1,  \
                             ta1_x_xxx_yzzz_0,   \
                             ta1_x_xxx_yzzz_1,   \
                             ta1_x_xxx_yzzzz_0,  \
                             ta1_x_xxx_yzzzz_1,  \
                             ta1_x_xxx_zzzz_0,   \
                             ta1_x_xxx_zzzz_1,   \
                             ta1_x_xxx_zzzzz_0,  \
                             ta1_x_xxx_zzzzz_1,  \
                             ta1_x_xxxz_xxxxx_0, \
                             ta1_x_xxxz_xxxxy_0, \
                             ta1_x_xxxz_xxxxz_0, \
                             ta1_x_xxxz_xxxyy_0, \
                             ta1_x_xxxz_xxxyz_0, \
                             ta1_x_xxxz_xxxzz_0, \
                             ta1_x_xxxz_xxyyy_0, \
                             ta1_x_xxxz_xxyyz_0, \
                             ta1_x_xxxz_xxyzz_0, \
                             ta1_x_xxxz_xxzzz_0, \
                             ta1_x_xxxz_xyyyy_0, \
                             ta1_x_xxxz_xyyyz_0, \
                             ta1_x_xxxz_xyyzz_0, \
                             ta1_x_xxxz_xyzzz_0, \
                             ta1_x_xxxz_xzzzz_0, \
                             ta1_x_xxxz_yyyyy_0, \
                             ta1_x_xxxz_yyyyz_0, \
                             ta1_x_xxxz_yyyzz_0, \
                             ta1_x_xxxz_yyzzz_0, \
                             ta1_x_xxxz_yzzzz_0, \
                             ta1_x_xxxz_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxxz_xxxxx_0[i] = ta1_x_xxx_xxxxx_0[i] * pa_z[i] - ta1_x_xxx_xxxxx_1[i] * pc_z[i];

        ta1_x_xxxz_xxxxy_0[i] = ta1_x_xxx_xxxxy_0[i] * pa_z[i] - ta1_x_xxx_xxxxy_1[i] * pc_z[i];

        ta1_x_xxxz_xxxxz_0[i] =
            ta1_x_xxx_xxxx_0[i] * fe_0 - ta1_x_xxx_xxxx_1[i] * fe_0 + ta1_x_xxx_xxxxz_0[i] * pa_z[i] - ta1_x_xxx_xxxxz_1[i] * pc_z[i];

        ta1_x_xxxz_xxxyy_0[i] = ta1_x_xxx_xxxyy_0[i] * pa_z[i] - ta1_x_xxx_xxxyy_1[i] * pc_z[i];

        ta1_x_xxxz_xxxyz_0[i] =
            ta1_x_xxx_xxxy_0[i] * fe_0 - ta1_x_xxx_xxxy_1[i] * fe_0 + ta1_x_xxx_xxxyz_0[i] * pa_z[i] - ta1_x_xxx_xxxyz_1[i] * pc_z[i];

        ta1_x_xxxz_xxxzz_0[i] =
            2.0 * ta1_x_xxx_xxxz_0[i] * fe_0 - 2.0 * ta1_x_xxx_xxxz_1[i] * fe_0 + ta1_x_xxx_xxxzz_0[i] * pa_z[i] - ta1_x_xxx_xxxzz_1[i] * pc_z[i];

        ta1_x_xxxz_xxyyy_0[i] = ta1_x_xxx_xxyyy_0[i] * pa_z[i] - ta1_x_xxx_xxyyy_1[i] * pc_z[i];

        ta1_x_xxxz_xxyyz_0[i] =
            ta1_x_xxx_xxyy_0[i] * fe_0 - ta1_x_xxx_xxyy_1[i] * fe_0 + ta1_x_xxx_xxyyz_0[i] * pa_z[i] - ta1_x_xxx_xxyyz_1[i] * pc_z[i];

        ta1_x_xxxz_xxyzz_0[i] =
            2.0 * ta1_x_xxx_xxyz_0[i] * fe_0 - 2.0 * ta1_x_xxx_xxyz_1[i] * fe_0 + ta1_x_xxx_xxyzz_0[i] * pa_z[i] - ta1_x_xxx_xxyzz_1[i] * pc_z[i];

        ta1_x_xxxz_xxzzz_0[i] =
            3.0 * ta1_x_xxx_xxzz_0[i] * fe_0 - 3.0 * ta1_x_xxx_xxzz_1[i] * fe_0 + ta1_x_xxx_xxzzz_0[i] * pa_z[i] - ta1_x_xxx_xxzzz_1[i] * pc_z[i];

        ta1_x_xxxz_xyyyy_0[i] = ta1_x_xxx_xyyyy_0[i] * pa_z[i] - ta1_x_xxx_xyyyy_1[i] * pc_z[i];

        ta1_x_xxxz_xyyyz_0[i] =
            ta1_x_xxx_xyyy_0[i] * fe_0 - ta1_x_xxx_xyyy_1[i] * fe_0 + ta1_x_xxx_xyyyz_0[i] * pa_z[i] - ta1_x_xxx_xyyyz_1[i] * pc_z[i];

        ta1_x_xxxz_xyyzz_0[i] =
            2.0 * ta1_x_xxx_xyyz_0[i] * fe_0 - 2.0 * ta1_x_xxx_xyyz_1[i] * fe_0 + ta1_x_xxx_xyyzz_0[i] * pa_z[i] - ta1_x_xxx_xyyzz_1[i] * pc_z[i];

        ta1_x_xxxz_xyzzz_0[i] =
            3.0 * ta1_x_xxx_xyzz_0[i] * fe_0 - 3.0 * ta1_x_xxx_xyzz_1[i] * fe_0 + ta1_x_xxx_xyzzz_0[i] * pa_z[i] - ta1_x_xxx_xyzzz_1[i] * pc_z[i];

        ta1_x_xxxz_xzzzz_0[i] =
            4.0 * ta1_x_xxx_xzzz_0[i] * fe_0 - 4.0 * ta1_x_xxx_xzzz_1[i] * fe_0 + ta1_x_xxx_xzzzz_0[i] * pa_z[i] - ta1_x_xxx_xzzzz_1[i] * pc_z[i];

        ta1_x_xxxz_yyyyy_0[i] = ta1_x_xxx_yyyyy_0[i] * pa_z[i] - ta1_x_xxx_yyyyy_1[i] * pc_z[i];

        ta1_x_xxxz_yyyyz_0[i] =
            ta1_x_xxx_yyyy_0[i] * fe_0 - ta1_x_xxx_yyyy_1[i] * fe_0 + ta1_x_xxx_yyyyz_0[i] * pa_z[i] - ta1_x_xxx_yyyyz_1[i] * pc_z[i];

        ta1_x_xxxz_yyyzz_0[i] =
            2.0 * ta1_x_xxx_yyyz_0[i] * fe_0 - 2.0 * ta1_x_xxx_yyyz_1[i] * fe_0 + ta1_x_xxx_yyyzz_0[i] * pa_z[i] - ta1_x_xxx_yyyzz_1[i] * pc_z[i];

        ta1_x_xxxz_yyzzz_0[i] =
            3.0 * ta1_x_xxx_yyzz_0[i] * fe_0 - 3.0 * ta1_x_xxx_yyzz_1[i] * fe_0 + ta1_x_xxx_yyzzz_0[i] * pa_z[i] - ta1_x_xxx_yyzzz_1[i] * pc_z[i];

        ta1_x_xxxz_yzzzz_0[i] =
            4.0 * ta1_x_xxx_yzzz_0[i] * fe_0 - 4.0 * ta1_x_xxx_yzzz_1[i] * fe_0 + ta1_x_xxx_yzzzz_0[i] * pa_z[i] - ta1_x_xxx_yzzzz_1[i] * pc_z[i];

        ta1_x_xxxz_zzzzz_0[i] =
            5.0 * ta1_x_xxx_zzzz_0[i] * fe_0 - 5.0 * ta1_x_xxx_zzzz_1[i] * fe_0 + ta1_x_xxx_zzzzz_0[i] * pa_z[i] - ta1_x_xxx_zzzzz_1[i] * pc_z[i];
    }

    // Set up 63-84 components of targeted buffer : GH

    auto ta1_x_xxyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 63);

    auto ta1_x_xxyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 64);

    auto ta1_x_xxyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 65);

    auto ta1_x_xxyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 66);

    auto ta1_x_xxyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 67);

    auto ta1_x_xxyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 68);

    auto ta1_x_xxyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 69);

    auto ta1_x_xxyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 70);

    auto ta1_x_xxyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 71);

    auto ta1_x_xxyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 72);

    auto ta1_x_xxyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 73);

    auto ta1_x_xxyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 74);

    auto ta1_x_xxyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 75);

    auto ta1_x_xxyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 76);

    auto ta1_x_xxyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 77);

    auto ta1_x_xxyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 78);

    auto ta1_x_xxyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 79);

    auto ta1_x_xxyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 80);

    auto ta1_x_xxyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 81);

    auto ta1_x_xxyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 82);

    auto ta1_x_xxyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 83);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_xx_xxxxx_0,   \
                             ta1_x_xx_xxxxx_1,   \
                             ta1_x_xx_xxxxy_0,   \
                             ta1_x_xx_xxxxy_1,   \
                             ta1_x_xx_xxxxz_0,   \
                             ta1_x_xx_xxxxz_1,   \
                             ta1_x_xx_xxxyy_0,   \
                             ta1_x_xx_xxxyy_1,   \
                             ta1_x_xx_xxxyz_0,   \
                             ta1_x_xx_xxxyz_1,   \
                             ta1_x_xx_xxxzz_0,   \
                             ta1_x_xx_xxxzz_1,   \
                             ta1_x_xx_xxyyy_0,   \
                             ta1_x_xx_xxyyy_1,   \
                             ta1_x_xx_xxyyz_0,   \
                             ta1_x_xx_xxyyz_1,   \
                             ta1_x_xx_xxyzz_0,   \
                             ta1_x_xx_xxyzz_1,   \
                             ta1_x_xx_xxzzz_0,   \
                             ta1_x_xx_xxzzz_1,   \
                             ta1_x_xx_xyyyy_0,   \
                             ta1_x_xx_xyyyy_1,   \
                             ta1_x_xx_xyyyz_0,   \
                             ta1_x_xx_xyyyz_1,   \
                             ta1_x_xx_xyyzz_0,   \
                             ta1_x_xx_xyyzz_1,   \
                             ta1_x_xx_xyzzz_0,   \
                             ta1_x_xx_xyzzz_1,   \
                             ta1_x_xx_xzzzz_0,   \
                             ta1_x_xx_xzzzz_1,   \
                             ta1_x_xx_zzzzz_0,   \
                             ta1_x_xx_zzzzz_1,   \
                             ta1_x_xxy_xxxx_0,   \
                             ta1_x_xxy_xxxx_1,   \
                             ta1_x_xxy_xxxxx_0,  \
                             ta1_x_xxy_xxxxx_1,  \
                             ta1_x_xxy_xxxxy_0,  \
                             ta1_x_xxy_xxxxy_1,  \
                             ta1_x_xxy_xxxxz_0,  \
                             ta1_x_xxy_xxxxz_1,  \
                             ta1_x_xxy_xxxy_0,   \
                             ta1_x_xxy_xxxy_1,   \
                             ta1_x_xxy_xxxyy_0,  \
                             ta1_x_xxy_xxxyy_1,  \
                             ta1_x_xxy_xxxyz_0,  \
                             ta1_x_xxy_xxxyz_1,  \
                             ta1_x_xxy_xxxz_0,   \
                             ta1_x_xxy_xxxz_1,   \
                             ta1_x_xxy_xxxzz_0,  \
                             ta1_x_xxy_xxxzz_1,  \
                             ta1_x_xxy_xxyy_0,   \
                             ta1_x_xxy_xxyy_1,   \
                             ta1_x_xxy_xxyyy_0,  \
                             ta1_x_xxy_xxyyy_1,  \
                             ta1_x_xxy_xxyyz_0,  \
                             ta1_x_xxy_xxyyz_1,  \
                             ta1_x_xxy_xxyz_0,   \
                             ta1_x_xxy_xxyz_1,   \
                             ta1_x_xxy_xxyzz_0,  \
                             ta1_x_xxy_xxyzz_1,  \
                             ta1_x_xxy_xxzz_0,   \
                             ta1_x_xxy_xxzz_1,   \
                             ta1_x_xxy_xxzzz_0,  \
                             ta1_x_xxy_xxzzz_1,  \
                             ta1_x_xxy_xyyy_0,   \
                             ta1_x_xxy_xyyy_1,   \
                             ta1_x_xxy_xyyyy_0,  \
                             ta1_x_xxy_xyyyy_1,  \
                             ta1_x_xxy_xyyyz_0,  \
                             ta1_x_xxy_xyyyz_1,  \
                             ta1_x_xxy_xyyz_0,   \
                             ta1_x_xxy_xyyz_1,   \
                             ta1_x_xxy_xyyzz_0,  \
                             ta1_x_xxy_xyyzz_1,  \
                             ta1_x_xxy_xyzz_0,   \
                             ta1_x_xxy_xyzz_1,   \
                             ta1_x_xxy_xyzzz_0,  \
                             ta1_x_xxy_xyzzz_1,  \
                             ta1_x_xxy_xzzz_0,   \
                             ta1_x_xxy_xzzz_1,   \
                             ta1_x_xxy_xzzzz_0,  \
                             ta1_x_xxy_xzzzz_1,  \
                             ta1_x_xxy_zzzzz_0,  \
                             ta1_x_xxy_zzzzz_1,  \
                             ta1_x_xxyy_xxxxx_0, \
                             ta1_x_xxyy_xxxxy_0, \
                             ta1_x_xxyy_xxxxz_0, \
                             ta1_x_xxyy_xxxyy_0, \
                             ta1_x_xxyy_xxxyz_0, \
                             ta1_x_xxyy_xxxzz_0, \
                             ta1_x_xxyy_xxyyy_0, \
                             ta1_x_xxyy_xxyyz_0, \
                             ta1_x_xxyy_xxyzz_0, \
                             ta1_x_xxyy_xxzzz_0, \
                             ta1_x_xxyy_xyyyy_0, \
                             ta1_x_xxyy_xyyyz_0, \
                             ta1_x_xxyy_xyyzz_0, \
                             ta1_x_xxyy_xyzzz_0, \
                             ta1_x_xxyy_xzzzz_0, \
                             ta1_x_xxyy_yyyyy_0, \
                             ta1_x_xxyy_yyyyz_0, \
                             ta1_x_xxyy_yyyzz_0, \
                             ta1_x_xxyy_yyzzz_0, \
                             ta1_x_xxyy_yzzzz_0, \
                             ta1_x_xxyy_zzzzz_0, \
                             ta1_x_xyy_yyyyy_0,  \
                             ta1_x_xyy_yyyyy_1,  \
                             ta1_x_xyy_yyyyz_0,  \
                             ta1_x_xyy_yyyyz_1,  \
                             ta1_x_xyy_yyyzz_0,  \
                             ta1_x_xyy_yyyzz_1,  \
                             ta1_x_xyy_yyzzz_0,  \
                             ta1_x_xyy_yyzzz_1,  \
                             ta1_x_xyy_yzzzz_0,  \
                             ta1_x_xyy_yzzzz_1,  \
                             ta1_x_yy_yyyyy_0,   \
                             ta1_x_yy_yyyyy_1,   \
                             ta1_x_yy_yyyyz_0,   \
                             ta1_x_yy_yyyyz_1,   \
                             ta1_x_yy_yyyzz_0,   \
                             ta1_x_yy_yyyzz_1,   \
                             ta1_x_yy_yyzzz_0,   \
                             ta1_x_yy_yyzzz_1,   \
                             ta1_x_yy_yzzzz_0,   \
                             ta1_x_yy_yzzzz_1,   \
                             ta_xyy_yyyyy_1,     \
                             ta_xyy_yyyyz_1,     \
                             ta_xyy_yyyzz_1,     \
                             ta_xyy_yyzzz_1,     \
                             ta_xyy_yzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyy_xxxxx_0[i] =
            ta1_x_xx_xxxxx_0[i] * fe_0 - ta1_x_xx_xxxxx_1[i] * fe_0 + ta1_x_xxy_xxxxx_0[i] * pa_y[i] - ta1_x_xxy_xxxxx_1[i] * pc_y[i];

        ta1_x_xxyy_xxxxy_0[i] = ta1_x_xx_xxxxy_0[i] * fe_0 - ta1_x_xx_xxxxy_1[i] * fe_0 + ta1_x_xxy_xxxx_0[i] * fe_0 - ta1_x_xxy_xxxx_1[i] * fe_0 +
                                ta1_x_xxy_xxxxy_0[i] * pa_y[i] - ta1_x_xxy_xxxxy_1[i] * pc_y[i];

        ta1_x_xxyy_xxxxz_0[i] =
            ta1_x_xx_xxxxz_0[i] * fe_0 - ta1_x_xx_xxxxz_1[i] * fe_0 + ta1_x_xxy_xxxxz_0[i] * pa_y[i] - ta1_x_xxy_xxxxz_1[i] * pc_y[i];

        ta1_x_xxyy_xxxyy_0[i] = ta1_x_xx_xxxyy_0[i] * fe_0 - ta1_x_xx_xxxyy_1[i] * fe_0 + 2.0 * ta1_x_xxy_xxxy_0[i] * fe_0 -
                                2.0 * ta1_x_xxy_xxxy_1[i] * fe_0 + ta1_x_xxy_xxxyy_0[i] * pa_y[i] - ta1_x_xxy_xxxyy_1[i] * pc_y[i];

        ta1_x_xxyy_xxxyz_0[i] = ta1_x_xx_xxxyz_0[i] * fe_0 - ta1_x_xx_xxxyz_1[i] * fe_0 + ta1_x_xxy_xxxz_0[i] * fe_0 - ta1_x_xxy_xxxz_1[i] * fe_0 +
                                ta1_x_xxy_xxxyz_0[i] * pa_y[i] - ta1_x_xxy_xxxyz_1[i] * pc_y[i];

        ta1_x_xxyy_xxxzz_0[i] =
            ta1_x_xx_xxxzz_0[i] * fe_0 - ta1_x_xx_xxxzz_1[i] * fe_0 + ta1_x_xxy_xxxzz_0[i] * pa_y[i] - ta1_x_xxy_xxxzz_1[i] * pc_y[i];

        ta1_x_xxyy_xxyyy_0[i] = ta1_x_xx_xxyyy_0[i] * fe_0 - ta1_x_xx_xxyyy_1[i] * fe_0 + 3.0 * ta1_x_xxy_xxyy_0[i] * fe_0 -
                                3.0 * ta1_x_xxy_xxyy_1[i] * fe_0 + ta1_x_xxy_xxyyy_0[i] * pa_y[i] - ta1_x_xxy_xxyyy_1[i] * pc_y[i];

        ta1_x_xxyy_xxyyz_0[i] = ta1_x_xx_xxyyz_0[i] * fe_0 - ta1_x_xx_xxyyz_1[i] * fe_0 + 2.0 * ta1_x_xxy_xxyz_0[i] * fe_0 -
                                2.0 * ta1_x_xxy_xxyz_1[i] * fe_0 + ta1_x_xxy_xxyyz_0[i] * pa_y[i] - ta1_x_xxy_xxyyz_1[i] * pc_y[i];

        ta1_x_xxyy_xxyzz_0[i] = ta1_x_xx_xxyzz_0[i] * fe_0 - ta1_x_xx_xxyzz_1[i] * fe_0 + ta1_x_xxy_xxzz_0[i] * fe_0 - ta1_x_xxy_xxzz_1[i] * fe_0 +
                                ta1_x_xxy_xxyzz_0[i] * pa_y[i] - ta1_x_xxy_xxyzz_1[i] * pc_y[i];

        ta1_x_xxyy_xxzzz_0[i] =
            ta1_x_xx_xxzzz_0[i] * fe_0 - ta1_x_xx_xxzzz_1[i] * fe_0 + ta1_x_xxy_xxzzz_0[i] * pa_y[i] - ta1_x_xxy_xxzzz_1[i] * pc_y[i];

        ta1_x_xxyy_xyyyy_0[i] = ta1_x_xx_xyyyy_0[i] * fe_0 - ta1_x_xx_xyyyy_1[i] * fe_0 + 4.0 * ta1_x_xxy_xyyy_0[i] * fe_0 -
                                4.0 * ta1_x_xxy_xyyy_1[i] * fe_0 + ta1_x_xxy_xyyyy_0[i] * pa_y[i] - ta1_x_xxy_xyyyy_1[i] * pc_y[i];

        ta1_x_xxyy_xyyyz_0[i] = ta1_x_xx_xyyyz_0[i] * fe_0 - ta1_x_xx_xyyyz_1[i] * fe_0 + 3.0 * ta1_x_xxy_xyyz_0[i] * fe_0 -
                                3.0 * ta1_x_xxy_xyyz_1[i] * fe_0 + ta1_x_xxy_xyyyz_0[i] * pa_y[i] - ta1_x_xxy_xyyyz_1[i] * pc_y[i];

        ta1_x_xxyy_xyyzz_0[i] = ta1_x_xx_xyyzz_0[i] * fe_0 - ta1_x_xx_xyyzz_1[i] * fe_0 + 2.0 * ta1_x_xxy_xyzz_0[i] * fe_0 -
                                2.0 * ta1_x_xxy_xyzz_1[i] * fe_0 + ta1_x_xxy_xyyzz_0[i] * pa_y[i] - ta1_x_xxy_xyyzz_1[i] * pc_y[i];

        ta1_x_xxyy_xyzzz_0[i] = ta1_x_xx_xyzzz_0[i] * fe_0 - ta1_x_xx_xyzzz_1[i] * fe_0 + ta1_x_xxy_xzzz_0[i] * fe_0 - ta1_x_xxy_xzzz_1[i] * fe_0 +
                                ta1_x_xxy_xyzzz_0[i] * pa_y[i] - ta1_x_xxy_xyzzz_1[i] * pc_y[i];

        ta1_x_xxyy_xzzzz_0[i] =
            ta1_x_xx_xzzzz_0[i] * fe_0 - ta1_x_xx_xzzzz_1[i] * fe_0 + ta1_x_xxy_xzzzz_0[i] * pa_y[i] - ta1_x_xxy_xzzzz_1[i] * pc_y[i];

        ta1_x_xxyy_yyyyy_0[i] = ta1_x_yy_yyyyy_0[i] * fe_0 - ta1_x_yy_yyyyy_1[i] * fe_0 + ta_xyy_yyyyy_1[i] + ta1_x_xyy_yyyyy_0[i] * pa_x[i] -
                                ta1_x_xyy_yyyyy_1[i] * pc_x[i];

        ta1_x_xxyy_yyyyz_0[i] = ta1_x_yy_yyyyz_0[i] * fe_0 - ta1_x_yy_yyyyz_1[i] * fe_0 + ta_xyy_yyyyz_1[i] + ta1_x_xyy_yyyyz_0[i] * pa_x[i] -
                                ta1_x_xyy_yyyyz_1[i] * pc_x[i];

        ta1_x_xxyy_yyyzz_0[i] = ta1_x_yy_yyyzz_0[i] * fe_0 - ta1_x_yy_yyyzz_1[i] * fe_0 + ta_xyy_yyyzz_1[i] + ta1_x_xyy_yyyzz_0[i] * pa_x[i] -
                                ta1_x_xyy_yyyzz_1[i] * pc_x[i];

        ta1_x_xxyy_yyzzz_0[i] = ta1_x_yy_yyzzz_0[i] * fe_0 - ta1_x_yy_yyzzz_1[i] * fe_0 + ta_xyy_yyzzz_1[i] + ta1_x_xyy_yyzzz_0[i] * pa_x[i] -
                                ta1_x_xyy_yyzzz_1[i] * pc_x[i];

        ta1_x_xxyy_yzzzz_0[i] = ta1_x_yy_yzzzz_0[i] * fe_0 - ta1_x_yy_yzzzz_1[i] * fe_0 + ta_xyy_yzzzz_1[i] + ta1_x_xyy_yzzzz_0[i] * pa_x[i] -
                                ta1_x_xyy_yzzzz_1[i] * pc_x[i];

        ta1_x_xxyy_zzzzz_0[i] =
            ta1_x_xx_zzzzz_0[i] * fe_0 - ta1_x_xx_zzzzz_1[i] * fe_0 + ta1_x_xxy_zzzzz_0[i] * pa_y[i] - ta1_x_xxy_zzzzz_1[i] * pc_y[i];
    }

    // Set up 84-105 components of targeted buffer : GH

    auto ta1_x_xxyz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 84);

    auto ta1_x_xxyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 85);

    auto ta1_x_xxyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 86);

    auto ta1_x_xxyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 87);

    auto ta1_x_xxyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 88);

    auto ta1_x_xxyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 89);

    auto ta1_x_xxyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 90);

    auto ta1_x_xxyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 91);

    auto ta1_x_xxyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 92);

    auto ta1_x_xxyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 93);

    auto ta1_x_xxyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 94);

    auto ta1_x_xxyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 95);

    auto ta1_x_xxyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 96);

    auto ta1_x_xxyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 97);

    auto ta1_x_xxyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 98);

    auto ta1_x_xxyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 99);

    auto ta1_x_xxyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 100);

    auto ta1_x_xxyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 101);

    auto ta1_x_xxyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 102);

    auto ta1_x_xxyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 103);

    auto ta1_x_xxyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 104);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xxy_xxxxy_0,  \
                             ta1_x_xxy_xxxxy_1,  \
                             ta1_x_xxy_xxxyy_0,  \
                             ta1_x_xxy_xxxyy_1,  \
                             ta1_x_xxy_xxyyy_0,  \
                             ta1_x_xxy_xxyyy_1,  \
                             ta1_x_xxy_xyyyy_0,  \
                             ta1_x_xxy_xyyyy_1,  \
                             ta1_x_xxy_yyyyy_0,  \
                             ta1_x_xxy_yyyyy_1,  \
                             ta1_x_xxyz_xxxxx_0, \
                             ta1_x_xxyz_xxxxy_0, \
                             ta1_x_xxyz_xxxxz_0, \
                             ta1_x_xxyz_xxxyy_0, \
                             ta1_x_xxyz_xxxyz_0, \
                             ta1_x_xxyz_xxxzz_0, \
                             ta1_x_xxyz_xxyyy_0, \
                             ta1_x_xxyz_xxyyz_0, \
                             ta1_x_xxyz_xxyzz_0, \
                             ta1_x_xxyz_xxzzz_0, \
                             ta1_x_xxyz_xyyyy_0, \
                             ta1_x_xxyz_xyyyz_0, \
                             ta1_x_xxyz_xyyzz_0, \
                             ta1_x_xxyz_xyzzz_0, \
                             ta1_x_xxyz_xzzzz_0, \
                             ta1_x_xxyz_yyyyy_0, \
                             ta1_x_xxyz_yyyyz_0, \
                             ta1_x_xxyz_yyyzz_0, \
                             ta1_x_xxyz_yyzzz_0, \
                             ta1_x_xxyz_yzzzz_0, \
                             ta1_x_xxyz_zzzzz_0, \
                             ta1_x_xxz_xxxxx_0,  \
                             ta1_x_xxz_xxxxx_1,  \
                             ta1_x_xxz_xxxxz_0,  \
                             ta1_x_xxz_xxxxz_1,  \
                             ta1_x_xxz_xxxyz_0,  \
                             ta1_x_xxz_xxxyz_1,  \
                             ta1_x_xxz_xxxz_0,   \
                             ta1_x_xxz_xxxz_1,   \
                             ta1_x_xxz_xxxzz_0,  \
                             ta1_x_xxz_xxxzz_1,  \
                             ta1_x_xxz_xxyyz_0,  \
                             ta1_x_xxz_xxyyz_1,  \
                             ta1_x_xxz_xxyz_0,   \
                             ta1_x_xxz_xxyz_1,   \
                             ta1_x_xxz_xxyzz_0,  \
                             ta1_x_xxz_xxyzz_1,  \
                             ta1_x_xxz_xxzz_0,   \
                             ta1_x_xxz_xxzz_1,   \
                             ta1_x_xxz_xxzzz_0,  \
                             ta1_x_xxz_xxzzz_1,  \
                             ta1_x_xxz_xyyyz_0,  \
                             ta1_x_xxz_xyyyz_1,  \
                             ta1_x_xxz_xyyz_0,   \
                             ta1_x_xxz_xyyz_1,   \
                             ta1_x_xxz_xyyzz_0,  \
                             ta1_x_xxz_xyyzz_1,  \
                             ta1_x_xxz_xyzz_0,   \
                             ta1_x_xxz_xyzz_1,   \
                             ta1_x_xxz_xyzzz_0,  \
                             ta1_x_xxz_xyzzz_1,  \
                             ta1_x_xxz_xzzz_0,   \
                             ta1_x_xxz_xzzz_1,   \
                             ta1_x_xxz_xzzzz_0,  \
                             ta1_x_xxz_xzzzz_1,  \
                             ta1_x_xxz_yyyyz_0,  \
                             ta1_x_xxz_yyyyz_1,  \
                             ta1_x_xxz_yyyz_0,   \
                             ta1_x_xxz_yyyz_1,   \
                             ta1_x_xxz_yyyzz_0,  \
                             ta1_x_xxz_yyyzz_1,  \
                             ta1_x_xxz_yyzz_0,   \
                             ta1_x_xxz_yyzz_1,   \
                             ta1_x_xxz_yyzzz_0,  \
                             ta1_x_xxz_yyzzz_1,  \
                             ta1_x_xxz_yzzz_0,   \
                             ta1_x_xxz_yzzz_1,   \
                             ta1_x_xxz_yzzzz_0,  \
                             ta1_x_xxz_yzzzz_1,  \
                             ta1_x_xxz_zzzz_0,   \
                             ta1_x_xxz_zzzz_1,   \
                             ta1_x_xxz_zzzzz_0,  \
                             ta1_x_xxz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxyz_xxxxx_0[i] = ta1_x_xxz_xxxxx_0[i] * pa_y[i] - ta1_x_xxz_xxxxx_1[i] * pc_y[i];

        ta1_x_xxyz_xxxxy_0[i] = ta1_x_xxy_xxxxy_0[i] * pa_z[i] - ta1_x_xxy_xxxxy_1[i] * pc_z[i];

        ta1_x_xxyz_xxxxz_0[i] = ta1_x_xxz_xxxxz_0[i] * pa_y[i] - ta1_x_xxz_xxxxz_1[i] * pc_y[i];

        ta1_x_xxyz_xxxyy_0[i] = ta1_x_xxy_xxxyy_0[i] * pa_z[i] - ta1_x_xxy_xxxyy_1[i] * pc_z[i];

        ta1_x_xxyz_xxxyz_0[i] =
            ta1_x_xxz_xxxz_0[i] * fe_0 - ta1_x_xxz_xxxz_1[i] * fe_0 + ta1_x_xxz_xxxyz_0[i] * pa_y[i] - ta1_x_xxz_xxxyz_1[i] * pc_y[i];

        ta1_x_xxyz_xxxzz_0[i] = ta1_x_xxz_xxxzz_0[i] * pa_y[i] - ta1_x_xxz_xxxzz_1[i] * pc_y[i];

        ta1_x_xxyz_xxyyy_0[i] = ta1_x_xxy_xxyyy_0[i] * pa_z[i] - ta1_x_xxy_xxyyy_1[i] * pc_z[i];

        ta1_x_xxyz_xxyyz_0[i] =
            2.0 * ta1_x_xxz_xxyz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xxyz_1[i] * fe_0 + ta1_x_xxz_xxyyz_0[i] * pa_y[i] - ta1_x_xxz_xxyyz_1[i] * pc_y[i];

        ta1_x_xxyz_xxyzz_0[i] =
            ta1_x_xxz_xxzz_0[i] * fe_0 - ta1_x_xxz_xxzz_1[i] * fe_0 + ta1_x_xxz_xxyzz_0[i] * pa_y[i] - ta1_x_xxz_xxyzz_1[i] * pc_y[i];

        ta1_x_xxyz_xxzzz_0[i] = ta1_x_xxz_xxzzz_0[i] * pa_y[i] - ta1_x_xxz_xxzzz_1[i] * pc_y[i];

        ta1_x_xxyz_xyyyy_0[i] = ta1_x_xxy_xyyyy_0[i] * pa_z[i] - ta1_x_xxy_xyyyy_1[i] * pc_z[i];

        ta1_x_xxyz_xyyyz_0[i] =
            3.0 * ta1_x_xxz_xyyz_0[i] * fe_0 - 3.0 * ta1_x_xxz_xyyz_1[i] * fe_0 + ta1_x_xxz_xyyyz_0[i] * pa_y[i] - ta1_x_xxz_xyyyz_1[i] * pc_y[i];

        ta1_x_xxyz_xyyzz_0[i] =
            2.0 * ta1_x_xxz_xyzz_0[i] * fe_0 - 2.0 * ta1_x_xxz_xyzz_1[i] * fe_0 + ta1_x_xxz_xyyzz_0[i] * pa_y[i] - ta1_x_xxz_xyyzz_1[i] * pc_y[i];

        ta1_x_xxyz_xyzzz_0[i] =
            ta1_x_xxz_xzzz_0[i] * fe_0 - ta1_x_xxz_xzzz_1[i] * fe_0 + ta1_x_xxz_xyzzz_0[i] * pa_y[i] - ta1_x_xxz_xyzzz_1[i] * pc_y[i];

        ta1_x_xxyz_xzzzz_0[i] = ta1_x_xxz_xzzzz_0[i] * pa_y[i] - ta1_x_xxz_xzzzz_1[i] * pc_y[i];

        ta1_x_xxyz_yyyyy_0[i] = ta1_x_xxy_yyyyy_0[i] * pa_z[i] - ta1_x_xxy_yyyyy_1[i] * pc_z[i];

        ta1_x_xxyz_yyyyz_0[i] =
            4.0 * ta1_x_xxz_yyyz_0[i] * fe_0 - 4.0 * ta1_x_xxz_yyyz_1[i] * fe_0 + ta1_x_xxz_yyyyz_0[i] * pa_y[i] - ta1_x_xxz_yyyyz_1[i] * pc_y[i];

        ta1_x_xxyz_yyyzz_0[i] =
            3.0 * ta1_x_xxz_yyzz_0[i] * fe_0 - 3.0 * ta1_x_xxz_yyzz_1[i] * fe_0 + ta1_x_xxz_yyyzz_0[i] * pa_y[i] - ta1_x_xxz_yyyzz_1[i] * pc_y[i];

        ta1_x_xxyz_yyzzz_0[i] =
            2.0 * ta1_x_xxz_yzzz_0[i] * fe_0 - 2.0 * ta1_x_xxz_yzzz_1[i] * fe_0 + ta1_x_xxz_yyzzz_0[i] * pa_y[i] - ta1_x_xxz_yyzzz_1[i] * pc_y[i];

        ta1_x_xxyz_yzzzz_0[i] =
            ta1_x_xxz_zzzz_0[i] * fe_0 - ta1_x_xxz_zzzz_1[i] * fe_0 + ta1_x_xxz_yzzzz_0[i] * pa_y[i] - ta1_x_xxz_yzzzz_1[i] * pc_y[i];

        ta1_x_xxyz_zzzzz_0[i] = ta1_x_xxz_zzzzz_0[i] * pa_y[i] - ta1_x_xxz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 105-126 components of targeted buffer : GH

    auto ta1_x_xxzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 105);

    auto ta1_x_xxzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 106);

    auto ta1_x_xxzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 107);

    auto ta1_x_xxzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 108);

    auto ta1_x_xxzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 109);

    auto ta1_x_xxzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 110);

    auto ta1_x_xxzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 111);

    auto ta1_x_xxzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 112);

    auto ta1_x_xxzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 113);

    auto ta1_x_xxzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 114);

    auto ta1_x_xxzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 115);

    auto ta1_x_xxzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 116);

    auto ta1_x_xxzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 117);

    auto ta1_x_xxzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 118);

    auto ta1_x_xxzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 119);

    auto ta1_x_xxzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 120);

    auto ta1_x_xxzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 121);

    auto ta1_x_xxzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 122);

    auto ta1_x_xxzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 123);

    auto ta1_x_xxzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 124);

    auto ta1_x_xxzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 125);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_x_xx_xxxxx_0,   \
                             ta1_x_xx_xxxxx_1,   \
                             ta1_x_xx_xxxxy_0,   \
                             ta1_x_xx_xxxxy_1,   \
                             ta1_x_xx_xxxxz_0,   \
                             ta1_x_xx_xxxxz_1,   \
                             ta1_x_xx_xxxyy_0,   \
                             ta1_x_xx_xxxyy_1,   \
                             ta1_x_xx_xxxyz_0,   \
                             ta1_x_xx_xxxyz_1,   \
                             ta1_x_xx_xxxzz_0,   \
                             ta1_x_xx_xxxzz_1,   \
                             ta1_x_xx_xxyyy_0,   \
                             ta1_x_xx_xxyyy_1,   \
                             ta1_x_xx_xxyyz_0,   \
                             ta1_x_xx_xxyyz_1,   \
                             ta1_x_xx_xxyzz_0,   \
                             ta1_x_xx_xxyzz_1,   \
                             ta1_x_xx_xxzzz_0,   \
                             ta1_x_xx_xxzzz_1,   \
                             ta1_x_xx_xyyyy_0,   \
                             ta1_x_xx_xyyyy_1,   \
                             ta1_x_xx_xyyyz_0,   \
                             ta1_x_xx_xyyyz_1,   \
                             ta1_x_xx_xyyzz_0,   \
                             ta1_x_xx_xyyzz_1,   \
                             ta1_x_xx_xyzzz_0,   \
                             ta1_x_xx_xyzzz_1,   \
                             ta1_x_xx_xzzzz_0,   \
                             ta1_x_xx_xzzzz_1,   \
                             ta1_x_xx_yyyyy_0,   \
                             ta1_x_xx_yyyyy_1,   \
                             ta1_x_xxz_xxxx_0,   \
                             ta1_x_xxz_xxxx_1,   \
                             ta1_x_xxz_xxxxx_0,  \
                             ta1_x_xxz_xxxxx_1,  \
                             ta1_x_xxz_xxxxy_0,  \
                             ta1_x_xxz_xxxxy_1,  \
                             ta1_x_xxz_xxxxz_0,  \
                             ta1_x_xxz_xxxxz_1,  \
                             ta1_x_xxz_xxxy_0,   \
                             ta1_x_xxz_xxxy_1,   \
                             ta1_x_xxz_xxxyy_0,  \
                             ta1_x_xxz_xxxyy_1,  \
                             ta1_x_xxz_xxxyz_0,  \
                             ta1_x_xxz_xxxyz_1,  \
                             ta1_x_xxz_xxxz_0,   \
                             ta1_x_xxz_xxxz_1,   \
                             ta1_x_xxz_xxxzz_0,  \
                             ta1_x_xxz_xxxzz_1,  \
                             ta1_x_xxz_xxyy_0,   \
                             ta1_x_xxz_xxyy_1,   \
                             ta1_x_xxz_xxyyy_0,  \
                             ta1_x_xxz_xxyyy_1,  \
                             ta1_x_xxz_xxyyz_0,  \
                             ta1_x_xxz_xxyyz_1,  \
                             ta1_x_xxz_xxyz_0,   \
                             ta1_x_xxz_xxyz_1,   \
                             ta1_x_xxz_xxyzz_0,  \
                             ta1_x_xxz_xxyzz_1,  \
                             ta1_x_xxz_xxzz_0,   \
                             ta1_x_xxz_xxzz_1,   \
                             ta1_x_xxz_xxzzz_0,  \
                             ta1_x_xxz_xxzzz_1,  \
                             ta1_x_xxz_xyyy_0,   \
                             ta1_x_xxz_xyyy_1,   \
                             ta1_x_xxz_xyyyy_0,  \
                             ta1_x_xxz_xyyyy_1,  \
                             ta1_x_xxz_xyyyz_0,  \
                             ta1_x_xxz_xyyyz_1,  \
                             ta1_x_xxz_xyyz_0,   \
                             ta1_x_xxz_xyyz_1,   \
                             ta1_x_xxz_xyyzz_0,  \
                             ta1_x_xxz_xyyzz_1,  \
                             ta1_x_xxz_xyzz_0,   \
                             ta1_x_xxz_xyzz_1,   \
                             ta1_x_xxz_xyzzz_0,  \
                             ta1_x_xxz_xyzzz_1,  \
                             ta1_x_xxz_xzzz_0,   \
                             ta1_x_xxz_xzzz_1,   \
                             ta1_x_xxz_xzzzz_0,  \
                             ta1_x_xxz_xzzzz_1,  \
                             ta1_x_xxz_yyyyy_0,  \
                             ta1_x_xxz_yyyyy_1,  \
                             ta1_x_xxzz_xxxxx_0, \
                             ta1_x_xxzz_xxxxy_0, \
                             ta1_x_xxzz_xxxxz_0, \
                             ta1_x_xxzz_xxxyy_0, \
                             ta1_x_xxzz_xxxyz_0, \
                             ta1_x_xxzz_xxxzz_0, \
                             ta1_x_xxzz_xxyyy_0, \
                             ta1_x_xxzz_xxyyz_0, \
                             ta1_x_xxzz_xxyzz_0, \
                             ta1_x_xxzz_xxzzz_0, \
                             ta1_x_xxzz_xyyyy_0, \
                             ta1_x_xxzz_xyyyz_0, \
                             ta1_x_xxzz_xyyzz_0, \
                             ta1_x_xxzz_xyzzz_0, \
                             ta1_x_xxzz_xzzzz_0, \
                             ta1_x_xxzz_yyyyy_0, \
                             ta1_x_xxzz_yyyyz_0, \
                             ta1_x_xxzz_yyyzz_0, \
                             ta1_x_xxzz_yyzzz_0, \
                             ta1_x_xxzz_yzzzz_0, \
                             ta1_x_xxzz_zzzzz_0, \
                             ta1_x_xzz_yyyyz_0,  \
                             ta1_x_xzz_yyyyz_1,  \
                             ta1_x_xzz_yyyzz_0,  \
                             ta1_x_xzz_yyyzz_1,  \
                             ta1_x_xzz_yyzzz_0,  \
                             ta1_x_xzz_yyzzz_1,  \
                             ta1_x_xzz_yzzzz_0,  \
                             ta1_x_xzz_yzzzz_1,  \
                             ta1_x_xzz_zzzzz_0,  \
                             ta1_x_xzz_zzzzz_1,  \
                             ta1_x_zz_yyyyz_0,   \
                             ta1_x_zz_yyyyz_1,   \
                             ta1_x_zz_yyyzz_0,   \
                             ta1_x_zz_yyyzz_1,   \
                             ta1_x_zz_yyzzz_0,   \
                             ta1_x_zz_yyzzz_1,   \
                             ta1_x_zz_yzzzz_0,   \
                             ta1_x_zz_yzzzz_1,   \
                             ta1_x_zz_zzzzz_0,   \
                             ta1_x_zz_zzzzz_1,   \
                             ta_xzz_yyyyz_1,     \
                             ta_xzz_yyyzz_1,     \
                             ta_xzz_yyzzz_1,     \
                             ta_xzz_yzzzz_1,     \
                             ta_xzz_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xxzz_xxxxx_0[i] =
            ta1_x_xx_xxxxx_0[i] * fe_0 - ta1_x_xx_xxxxx_1[i] * fe_0 + ta1_x_xxz_xxxxx_0[i] * pa_z[i] - ta1_x_xxz_xxxxx_1[i] * pc_z[i];

        ta1_x_xxzz_xxxxy_0[i] =
            ta1_x_xx_xxxxy_0[i] * fe_0 - ta1_x_xx_xxxxy_1[i] * fe_0 + ta1_x_xxz_xxxxy_0[i] * pa_z[i] - ta1_x_xxz_xxxxy_1[i] * pc_z[i];

        ta1_x_xxzz_xxxxz_0[i] = ta1_x_xx_xxxxz_0[i] * fe_0 - ta1_x_xx_xxxxz_1[i] * fe_0 + ta1_x_xxz_xxxx_0[i] * fe_0 - ta1_x_xxz_xxxx_1[i] * fe_0 +
                                ta1_x_xxz_xxxxz_0[i] * pa_z[i] - ta1_x_xxz_xxxxz_1[i] * pc_z[i];

        ta1_x_xxzz_xxxyy_0[i] =
            ta1_x_xx_xxxyy_0[i] * fe_0 - ta1_x_xx_xxxyy_1[i] * fe_0 + ta1_x_xxz_xxxyy_0[i] * pa_z[i] - ta1_x_xxz_xxxyy_1[i] * pc_z[i];

        ta1_x_xxzz_xxxyz_0[i] = ta1_x_xx_xxxyz_0[i] * fe_0 - ta1_x_xx_xxxyz_1[i] * fe_0 + ta1_x_xxz_xxxy_0[i] * fe_0 - ta1_x_xxz_xxxy_1[i] * fe_0 +
                                ta1_x_xxz_xxxyz_0[i] * pa_z[i] - ta1_x_xxz_xxxyz_1[i] * pc_z[i];

        ta1_x_xxzz_xxxzz_0[i] = ta1_x_xx_xxxzz_0[i] * fe_0 - ta1_x_xx_xxxzz_1[i] * fe_0 + 2.0 * ta1_x_xxz_xxxz_0[i] * fe_0 -
                                2.0 * ta1_x_xxz_xxxz_1[i] * fe_0 + ta1_x_xxz_xxxzz_0[i] * pa_z[i] - ta1_x_xxz_xxxzz_1[i] * pc_z[i];

        ta1_x_xxzz_xxyyy_0[i] =
            ta1_x_xx_xxyyy_0[i] * fe_0 - ta1_x_xx_xxyyy_1[i] * fe_0 + ta1_x_xxz_xxyyy_0[i] * pa_z[i] - ta1_x_xxz_xxyyy_1[i] * pc_z[i];

        ta1_x_xxzz_xxyyz_0[i] = ta1_x_xx_xxyyz_0[i] * fe_0 - ta1_x_xx_xxyyz_1[i] * fe_0 + ta1_x_xxz_xxyy_0[i] * fe_0 - ta1_x_xxz_xxyy_1[i] * fe_0 +
                                ta1_x_xxz_xxyyz_0[i] * pa_z[i] - ta1_x_xxz_xxyyz_1[i] * pc_z[i];

        ta1_x_xxzz_xxyzz_0[i] = ta1_x_xx_xxyzz_0[i] * fe_0 - ta1_x_xx_xxyzz_1[i] * fe_0 + 2.0 * ta1_x_xxz_xxyz_0[i] * fe_0 -
                                2.0 * ta1_x_xxz_xxyz_1[i] * fe_0 + ta1_x_xxz_xxyzz_0[i] * pa_z[i] - ta1_x_xxz_xxyzz_1[i] * pc_z[i];

        ta1_x_xxzz_xxzzz_0[i] = ta1_x_xx_xxzzz_0[i] * fe_0 - ta1_x_xx_xxzzz_1[i] * fe_0 + 3.0 * ta1_x_xxz_xxzz_0[i] * fe_0 -
                                3.0 * ta1_x_xxz_xxzz_1[i] * fe_0 + ta1_x_xxz_xxzzz_0[i] * pa_z[i] - ta1_x_xxz_xxzzz_1[i] * pc_z[i];

        ta1_x_xxzz_xyyyy_0[i] =
            ta1_x_xx_xyyyy_0[i] * fe_0 - ta1_x_xx_xyyyy_1[i] * fe_0 + ta1_x_xxz_xyyyy_0[i] * pa_z[i] - ta1_x_xxz_xyyyy_1[i] * pc_z[i];

        ta1_x_xxzz_xyyyz_0[i] = ta1_x_xx_xyyyz_0[i] * fe_0 - ta1_x_xx_xyyyz_1[i] * fe_0 + ta1_x_xxz_xyyy_0[i] * fe_0 - ta1_x_xxz_xyyy_1[i] * fe_0 +
                                ta1_x_xxz_xyyyz_0[i] * pa_z[i] - ta1_x_xxz_xyyyz_1[i] * pc_z[i];

        ta1_x_xxzz_xyyzz_0[i] = ta1_x_xx_xyyzz_0[i] * fe_0 - ta1_x_xx_xyyzz_1[i] * fe_0 + 2.0 * ta1_x_xxz_xyyz_0[i] * fe_0 -
                                2.0 * ta1_x_xxz_xyyz_1[i] * fe_0 + ta1_x_xxz_xyyzz_0[i] * pa_z[i] - ta1_x_xxz_xyyzz_1[i] * pc_z[i];

        ta1_x_xxzz_xyzzz_0[i] = ta1_x_xx_xyzzz_0[i] * fe_0 - ta1_x_xx_xyzzz_1[i] * fe_0 + 3.0 * ta1_x_xxz_xyzz_0[i] * fe_0 -
                                3.0 * ta1_x_xxz_xyzz_1[i] * fe_0 + ta1_x_xxz_xyzzz_0[i] * pa_z[i] - ta1_x_xxz_xyzzz_1[i] * pc_z[i];

        ta1_x_xxzz_xzzzz_0[i] = ta1_x_xx_xzzzz_0[i] * fe_0 - ta1_x_xx_xzzzz_1[i] * fe_0 + 4.0 * ta1_x_xxz_xzzz_0[i] * fe_0 -
                                4.0 * ta1_x_xxz_xzzz_1[i] * fe_0 + ta1_x_xxz_xzzzz_0[i] * pa_z[i] - ta1_x_xxz_xzzzz_1[i] * pc_z[i];

        ta1_x_xxzz_yyyyy_0[i] =
            ta1_x_xx_yyyyy_0[i] * fe_0 - ta1_x_xx_yyyyy_1[i] * fe_0 + ta1_x_xxz_yyyyy_0[i] * pa_z[i] - ta1_x_xxz_yyyyy_1[i] * pc_z[i];

        ta1_x_xxzz_yyyyz_0[i] = ta1_x_zz_yyyyz_0[i] * fe_0 - ta1_x_zz_yyyyz_1[i] * fe_0 + ta_xzz_yyyyz_1[i] + ta1_x_xzz_yyyyz_0[i] * pa_x[i] -
                                ta1_x_xzz_yyyyz_1[i] * pc_x[i];

        ta1_x_xxzz_yyyzz_0[i] = ta1_x_zz_yyyzz_0[i] * fe_0 - ta1_x_zz_yyyzz_1[i] * fe_0 + ta_xzz_yyyzz_1[i] + ta1_x_xzz_yyyzz_0[i] * pa_x[i] -
                                ta1_x_xzz_yyyzz_1[i] * pc_x[i];

        ta1_x_xxzz_yyzzz_0[i] = ta1_x_zz_yyzzz_0[i] * fe_0 - ta1_x_zz_yyzzz_1[i] * fe_0 + ta_xzz_yyzzz_1[i] + ta1_x_xzz_yyzzz_0[i] * pa_x[i] -
                                ta1_x_xzz_yyzzz_1[i] * pc_x[i];

        ta1_x_xxzz_yzzzz_0[i] = ta1_x_zz_yzzzz_0[i] * fe_0 - ta1_x_zz_yzzzz_1[i] * fe_0 + ta_xzz_yzzzz_1[i] + ta1_x_xzz_yzzzz_0[i] * pa_x[i] -
                                ta1_x_xzz_yzzzz_1[i] * pc_x[i];

        ta1_x_xxzz_zzzzz_0[i] = ta1_x_zz_zzzzz_0[i] * fe_0 - ta1_x_zz_zzzzz_1[i] * fe_0 + ta_xzz_zzzzz_1[i] + ta1_x_xzz_zzzzz_0[i] * pa_x[i] -
                                ta1_x_xzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 126-147 components of targeted buffer : GH

    auto ta1_x_xyyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 126);

    auto ta1_x_xyyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 127);

    auto ta1_x_xyyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 128);

    auto ta1_x_xyyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 129);

    auto ta1_x_xyyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 130);

    auto ta1_x_xyyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 131);

    auto ta1_x_xyyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 132);

    auto ta1_x_xyyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 133);

    auto ta1_x_xyyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 134);

    auto ta1_x_xyyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 135);

    auto ta1_x_xyyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 136);

    auto ta1_x_xyyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 137);

    auto ta1_x_xyyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 138);

    auto ta1_x_xyyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 139);

    auto ta1_x_xyyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 140);

    auto ta1_x_xyyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 141);

    auto ta1_x_xyyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 142);

    auto ta1_x_xyyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 143);

    auto ta1_x_xyyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 144);

    auto ta1_x_xyyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 145);

    auto ta1_x_xyyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 146);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_xy_xxxxx_0,   \
                             ta1_x_xy_xxxxx_1,   \
                             ta1_x_xy_xxxxz_0,   \
                             ta1_x_xy_xxxxz_1,   \
                             ta1_x_xy_xxxzz_0,   \
                             ta1_x_xy_xxxzz_1,   \
                             ta1_x_xy_xxzzz_0,   \
                             ta1_x_xy_xxzzz_1,   \
                             ta1_x_xy_xzzzz_0,   \
                             ta1_x_xy_xzzzz_1,   \
                             ta1_x_xyy_xxxxx_0,  \
                             ta1_x_xyy_xxxxx_1,  \
                             ta1_x_xyy_xxxxz_0,  \
                             ta1_x_xyy_xxxxz_1,  \
                             ta1_x_xyy_xxxzz_0,  \
                             ta1_x_xyy_xxxzz_1,  \
                             ta1_x_xyy_xxzzz_0,  \
                             ta1_x_xyy_xxzzz_1,  \
                             ta1_x_xyy_xzzzz_0,  \
                             ta1_x_xyy_xzzzz_1,  \
                             ta1_x_xyyy_xxxxx_0, \
                             ta1_x_xyyy_xxxxy_0, \
                             ta1_x_xyyy_xxxxz_0, \
                             ta1_x_xyyy_xxxyy_0, \
                             ta1_x_xyyy_xxxyz_0, \
                             ta1_x_xyyy_xxxzz_0, \
                             ta1_x_xyyy_xxyyy_0, \
                             ta1_x_xyyy_xxyyz_0, \
                             ta1_x_xyyy_xxyzz_0, \
                             ta1_x_xyyy_xxzzz_0, \
                             ta1_x_xyyy_xyyyy_0, \
                             ta1_x_xyyy_xyyyz_0, \
                             ta1_x_xyyy_xyyzz_0, \
                             ta1_x_xyyy_xyzzz_0, \
                             ta1_x_xyyy_xzzzz_0, \
                             ta1_x_xyyy_yyyyy_0, \
                             ta1_x_xyyy_yyyyz_0, \
                             ta1_x_xyyy_yyyzz_0, \
                             ta1_x_xyyy_yyzzz_0, \
                             ta1_x_xyyy_yzzzz_0, \
                             ta1_x_xyyy_zzzzz_0, \
                             ta1_x_yyy_xxxxy_0,  \
                             ta1_x_yyy_xxxxy_1,  \
                             ta1_x_yyy_xxxy_0,   \
                             ta1_x_yyy_xxxy_1,   \
                             ta1_x_yyy_xxxyy_0,  \
                             ta1_x_yyy_xxxyy_1,  \
                             ta1_x_yyy_xxxyz_0,  \
                             ta1_x_yyy_xxxyz_1,  \
                             ta1_x_yyy_xxyy_0,   \
                             ta1_x_yyy_xxyy_1,   \
                             ta1_x_yyy_xxyyy_0,  \
                             ta1_x_yyy_xxyyy_1,  \
                             ta1_x_yyy_xxyyz_0,  \
                             ta1_x_yyy_xxyyz_1,  \
                             ta1_x_yyy_xxyz_0,   \
                             ta1_x_yyy_xxyz_1,   \
                             ta1_x_yyy_xxyzz_0,  \
                             ta1_x_yyy_xxyzz_1,  \
                             ta1_x_yyy_xyyy_0,   \
                             ta1_x_yyy_xyyy_1,   \
                             ta1_x_yyy_xyyyy_0,  \
                             ta1_x_yyy_xyyyy_1,  \
                             ta1_x_yyy_xyyyz_0,  \
                             ta1_x_yyy_xyyyz_1,  \
                             ta1_x_yyy_xyyz_0,   \
                             ta1_x_yyy_xyyz_1,   \
                             ta1_x_yyy_xyyzz_0,  \
                             ta1_x_yyy_xyyzz_1,  \
                             ta1_x_yyy_xyzz_0,   \
                             ta1_x_yyy_xyzz_1,   \
                             ta1_x_yyy_xyzzz_0,  \
                             ta1_x_yyy_xyzzz_1,  \
                             ta1_x_yyy_yyyy_0,   \
                             ta1_x_yyy_yyyy_1,   \
                             ta1_x_yyy_yyyyy_0,  \
                             ta1_x_yyy_yyyyy_1,  \
                             ta1_x_yyy_yyyyz_0,  \
                             ta1_x_yyy_yyyyz_1,  \
                             ta1_x_yyy_yyyz_0,   \
                             ta1_x_yyy_yyyz_1,   \
                             ta1_x_yyy_yyyzz_0,  \
                             ta1_x_yyy_yyyzz_1,  \
                             ta1_x_yyy_yyzz_0,   \
                             ta1_x_yyy_yyzz_1,   \
                             ta1_x_yyy_yyzzz_0,  \
                             ta1_x_yyy_yyzzz_1,  \
                             ta1_x_yyy_yzzz_0,   \
                             ta1_x_yyy_yzzz_1,   \
                             ta1_x_yyy_yzzzz_0,  \
                             ta1_x_yyy_yzzzz_1,  \
                             ta1_x_yyy_zzzzz_0,  \
                             ta1_x_yyy_zzzzz_1,  \
                             ta_yyy_xxxxy_1,     \
                             ta_yyy_xxxyy_1,     \
                             ta_yyy_xxxyz_1,     \
                             ta_yyy_xxyyy_1,     \
                             ta_yyy_xxyyz_1,     \
                             ta_yyy_xxyzz_1,     \
                             ta_yyy_xyyyy_1,     \
                             ta_yyy_xyyyz_1,     \
                             ta_yyy_xyyzz_1,     \
                             ta_yyy_xyzzz_1,     \
                             ta_yyy_yyyyy_1,     \
                             ta_yyy_yyyyz_1,     \
                             ta_yyy_yyyzz_1,     \
                             ta_yyy_yyzzz_1,     \
                             ta_yyy_yzzzz_1,     \
                             ta_yyy_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyy_xxxxx_0[i] =
            2.0 * ta1_x_xy_xxxxx_0[i] * fe_0 - 2.0 * ta1_x_xy_xxxxx_1[i] * fe_0 + ta1_x_xyy_xxxxx_0[i] * pa_y[i] - ta1_x_xyy_xxxxx_1[i] * pc_y[i];

        ta1_x_xyyy_xxxxy_0[i] = 4.0 * ta1_x_yyy_xxxy_0[i] * fe_0 - 4.0 * ta1_x_yyy_xxxy_1[i] * fe_0 + ta_yyy_xxxxy_1[i] +
                                ta1_x_yyy_xxxxy_0[i] * pa_x[i] - ta1_x_yyy_xxxxy_1[i] * pc_x[i];

        ta1_x_xyyy_xxxxz_0[i] =
            2.0 * ta1_x_xy_xxxxz_0[i] * fe_0 - 2.0 * ta1_x_xy_xxxxz_1[i] * fe_0 + ta1_x_xyy_xxxxz_0[i] * pa_y[i] - ta1_x_xyy_xxxxz_1[i] * pc_y[i];

        ta1_x_xyyy_xxxyy_0[i] = 3.0 * ta1_x_yyy_xxyy_0[i] * fe_0 - 3.0 * ta1_x_yyy_xxyy_1[i] * fe_0 + ta_yyy_xxxyy_1[i] +
                                ta1_x_yyy_xxxyy_0[i] * pa_x[i] - ta1_x_yyy_xxxyy_1[i] * pc_x[i];

        ta1_x_xyyy_xxxyz_0[i] = 3.0 * ta1_x_yyy_xxyz_0[i] * fe_0 - 3.0 * ta1_x_yyy_xxyz_1[i] * fe_0 + ta_yyy_xxxyz_1[i] +
                                ta1_x_yyy_xxxyz_0[i] * pa_x[i] - ta1_x_yyy_xxxyz_1[i] * pc_x[i];

        ta1_x_xyyy_xxxzz_0[i] =
            2.0 * ta1_x_xy_xxxzz_0[i] * fe_0 - 2.0 * ta1_x_xy_xxxzz_1[i] * fe_0 + ta1_x_xyy_xxxzz_0[i] * pa_y[i] - ta1_x_xyy_xxxzz_1[i] * pc_y[i];

        ta1_x_xyyy_xxyyy_0[i] = 2.0 * ta1_x_yyy_xyyy_0[i] * fe_0 - 2.0 * ta1_x_yyy_xyyy_1[i] * fe_0 + ta_yyy_xxyyy_1[i] +
                                ta1_x_yyy_xxyyy_0[i] * pa_x[i] - ta1_x_yyy_xxyyy_1[i] * pc_x[i];

        ta1_x_xyyy_xxyyz_0[i] = 2.0 * ta1_x_yyy_xyyz_0[i] * fe_0 - 2.0 * ta1_x_yyy_xyyz_1[i] * fe_0 + ta_yyy_xxyyz_1[i] +
                                ta1_x_yyy_xxyyz_0[i] * pa_x[i] - ta1_x_yyy_xxyyz_1[i] * pc_x[i];

        ta1_x_xyyy_xxyzz_0[i] = 2.0 * ta1_x_yyy_xyzz_0[i] * fe_0 - 2.0 * ta1_x_yyy_xyzz_1[i] * fe_0 + ta_yyy_xxyzz_1[i] +
                                ta1_x_yyy_xxyzz_0[i] * pa_x[i] - ta1_x_yyy_xxyzz_1[i] * pc_x[i];

        ta1_x_xyyy_xxzzz_0[i] =
            2.0 * ta1_x_xy_xxzzz_0[i] * fe_0 - 2.0 * ta1_x_xy_xxzzz_1[i] * fe_0 + ta1_x_xyy_xxzzz_0[i] * pa_y[i] - ta1_x_xyy_xxzzz_1[i] * pc_y[i];

        ta1_x_xyyy_xyyyy_0[i] = ta1_x_yyy_yyyy_0[i] * fe_0 - ta1_x_yyy_yyyy_1[i] * fe_0 + ta_yyy_xyyyy_1[i] + ta1_x_yyy_xyyyy_0[i] * pa_x[i] -
                                ta1_x_yyy_xyyyy_1[i] * pc_x[i];

        ta1_x_xyyy_xyyyz_0[i] = ta1_x_yyy_yyyz_0[i] * fe_0 - ta1_x_yyy_yyyz_1[i] * fe_0 + ta_yyy_xyyyz_1[i] + ta1_x_yyy_xyyyz_0[i] * pa_x[i] -
                                ta1_x_yyy_xyyyz_1[i] * pc_x[i];

        ta1_x_xyyy_xyyzz_0[i] = ta1_x_yyy_yyzz_0[i] * fe_0 - ta1_x_yyy_yyzz_1[i] * fe_0 + ta_yyy_xyyzz_1[i] + ta1_x_yyy_xyyzz_0[i] * pa_x[i] -
                                ta1_x_yyy_xyyzz_1[i] * pc_x[i];

        ta1_x_xyyy_xyzzz_0[i] = ta1_x_yyy_yzzz_0[i] * fe_0 - ta1_x_yyy_yzzz_1[i] * fe_0 + ta_yyy_xyzzz_1[i] + ta1_x_yyy_xyzzz_0[i] * pa_x[i] -
                                ta1_x_yyy_xyzzz_1[i] * pc_x[i];

        ta1_x_xyyy_xzzzz_0[i] =
            2.0 * ta1_x_xy_xzzzz_0[i] * fe_0 - 2.0 * ta1_x_xy_xzzzz_1[i] * fe_0 + ta1_x_xyy_xzzzz_0[i] * pa_y[i] - ta1_x_xyy_xzzzz_1[i] * pc_y[i];

        ta1_x_xyyy_yyyyy_0[i] = ta_yyy_yyyyy_1[i] + ta1_x_yyy_yyyyy_0[i] * pa_x[i] - ta1_x_yyy_yyyyy_1[i] * pc_x[i];

        ta1_x_xyyy_yyyyz_0[i] = ta_yyy_yyyyz_1[i] + ta1_x_yyy_yyyyz_0[i] * pa_x[i] - ta1_x_yyy_yyyyz_1[i] * pc_x[i];

        ta1_x_xyyy_yyyzz_0[i] = ta_yyy_yyyzz_1[i] + ta1_x_yyy_yyyzz_0[i] * pa_x[i] - ta1_x_yyy_yyyzz_1[i] * pc_x[i];

        ta1_x_xyyy_yyzzz_0[i] = ta_yyy_yyzzz_1[i] + ta1_x_yyy_yyzzz_0[i] * pa_x[i] - ta1_x_yyy_yyzzz_1[i] * pc_x[i];

        ta1_x_xyyy_yzzzz_0[i] = ta_yyy_yzzzz_1[i] + ta1_x_yyy_yzzzz_0[i] * pa_x[i] - ta1_x_yyy_yzzzz_1[i] * pc_x[i];

        ta1_x_xyyy_zzzzz_0[i] = ta_yyy_zzzzz_1[i] + ta1_x_yyy_zzzzz_0[i] * pa_x[i] - ta1_x_yyy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 147-168 components of targeted buffer : GH

    auto ta1_x_xyyz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 147);

    auto ta1_x_xyyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 148);

    auto ta1_x_xyyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 149);

    auto ta1_x_xyyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 150);

    auto ta1_x_xyyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 151);

    auto ta1_x_xyyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 152);

    auto ta1_x_xyyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 153);

    auto ta1_x_xyyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 154);

    auto ta1_x_xyyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 155);

    auto ta1_x_xyyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 156);

    auto ta1_x_xyyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 157);

    auto ta1_x_xyyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 158);

    auto ta1_x_xyyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 159);

    auto ta1_x_xyyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 160);

    auto ta1_x_xyyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 161);

    auto ta1_x_xyyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 162);

    auto ta1_x_xyyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 163);

    auto ta1_x_xyyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 164);

    auto ta1_x_xyyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 165);

    auto ta1_x_xyyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 166);

    auto ta1_x_xyyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 167);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_xyy_xxxxx_0,  \
                             ta1_x_xyy_xxxxx_1,  \
                             ta1_x_xyy_xxxxy_0,  \
                             ta1_x_xyy_xxxxy_1,  \
                             ta1_x_xyy_xxxy_0,   \
                             ta1_x_xyy_xxxy_1,   \
                             ta1_x_xyy_xxxyy_0,  \
                             ta1_x_xyy_xxxyy_1,  \
                             ta1_x_xyy_xxxyz_0,  \
                             ta1_x_xyy_xxxyz_1,  \
                             ta1_x_xyy_xxyy_0,   \
                             ta1_x_xyy_xxyy_1,   \
                             ta1_x_xyy_xxyyy_0,  \
                             ta1_x_xyy_xxyyy_1,  \
                             ta1_x_xyy_xxyyz_0,  \
                             ta1_x_xyy_xxyyz_1,  \
                             ta1_x_xyy_xxyz_0,   \
                             ta1_x_xyy_xxyz_1,   \
                             ta1_x_xyy_xxyzz_0,  \
                             ta1_x_xyy_xxyzz_1,  \
                             ta1_x_xyy_xyyy_0,   \
                             ta1_x_xyy_xyyy_1,   \
                             ta1_x_xyy_xyyyy_0,  \
                             ta1_x_xyy_xyyyy_1,  \
                             ta1_x_xyy_xyyyz_0,  \
                             ta1_x_xyy_xyyyz_1,  \
                             ta1_x_xyy_xyyz_0,   \
                             ta1_x_xyy_xyyz_1,   \
                             ta1_x_xyy_xyyzz_0,  \
                             ta1_x_xyy_xyyzz_1,  \
                             ta1_x_xyy_xyzz_0,   \
                             ta1_x_xyy_xyzz_1,   \
                             ta1_x_xyy_xyzzz_0,  \
                             ta1_x_xyy_xyzzz_1,  \
                             ta1_x_xyy_yyyyy_0,  \
                             ta1_x_xyy_yyyyy_1,  \
                             ta1_x_xyyz_xxxxx_0, \
                             ta1_x_xyyz_xxxxy_0, \
                             ta1_x_xyyz_xxxxz_0, \
                             ta1_x_xyyz_xxxyy_0, \
                             ta1_x_xyyz_xxxyz_0, \
                             ta1_x_xyyz_xxxzz_0, \
                             ta1_x_xyyz_xxyyy_0, \
                             ta1_x_xyyz_xxyyz_0, \
                             ta1_x_xyyz_xxyzz_0, \
                             ta1_x_xyyz_xxzzz_0, \
                             ta1_x_xyyz_xyyyy_0, \
                             ta1_x_xyyz_xyyyz_0, \
                             ta1_x_xyyz_xyyzz_0, \
                             ta1_x_xyyz_xyzzz_0, \
                             ta1_x_xyyz_xzzzz_0, \
                             ta1_x_xyyz_yyyyy_0, \
                             ta1_x_xyyz_yyyyz_0, \
                             ta1_x_xyyz_yyyzz_0, \
                             ta1_x_xyyz_yyzzz_0, \
                             ta1_x_xyyz_yzzzz_0, \
                             ta1_x_xyyz_zzzzz_0, \
                             ta1_x_xyz_xxxxz_0,  \
                             ta1_x_xyz_xxxxz_1,  \
                             ta1_x_xyz_xxxzz_0,  \
                             ta1_x_xyz_xxxzz_1,  \
                             ta1_x_xyz_xxzzz_0,  \
                             ta1_x_xyz_xxzzz_1,  \
                             ta1_x_xyz_xzzzz_0,  \
                             ta1_x_xyz_xzzzz_1,  \
                             ta1_x_xz_xxxxz_0,   \
                             ta1_x_xz_xxxxz_1,   \
                             ta1_x_xz_xxxzz_0,   \
                             ta1_x_xz_xxxzz_1,   \
                             ta1_x_xz_xxzzz_0,   \
                             ta1_x_xz_xxzzz_1,   \
                             ta1_x_xz_xzzzz_0,   \
                             ta1_x_xz_xzzzz_1,   \
                             ta1_x_yyz_yyyyz_0,  \
                             ta1_x_yyz_yyyyz_1,  \
                             ta1_x_yyz_yyyzz_0,  \
                             ta1_x_yyz_yyyzz_1,  \
                             ta1_x_yyz_yyzzz_0,  \
                             ta1_x_yyz_yyzzz_1,  \
                             ta1_x_yyz_yzzzz_0,  \
                             ta1_x_yyz_yzzzz_1,  \
                             ta1_x_yyz_zzzzz_0,  \
                             ta1_x_yyz_zzzzz_1,  \
                             ta_yyz_yyyyz_1,     \
                             ta_yyz_yyyzz_1,     \
                             ta_yyz_yyzzz_1,     \
                             ta_yyz_yzzzz_1,     \
                             ta_yyz_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyyz_xxxxx_0[i] = ta1_x_xyy_xxxxx_0[i] * pa_z[i] - ta1_x_xyy_xxxxx_1[i] * pc_z[i];

        ta1_x_xyyz_xxxxy_0[i] = ta1_x_xyy_xxxxy_0[i] * pa_z[i] - ta1_x_xyy_xxxxy_1[i] * pc_z[i];

        ta1_x_xyyz_xxxxz_0[i] =
            ta1_x_xz_xxxxz_0[i] * fe_0 - ta1_x_xz_xxxxz_1[i] * fe_0 + ta1_x_xyz_xxxxz_0[i] * pa_y[i] - ta1_x_xyz_xxxxz_1[i] * pc_y[i];

        ta1_x_xyyz_xxxyy_0[i] = ta1_x_xyy_xxxyy_0[i] * pa_z[i] - ta1_x_xyy_xxxyy_1[i] * pc_z[i];

        ta1_x_xyyz_xxxyz_0[i] =
            ta1_x_xyy_xxxy_0[i] * fe_0 - ta1_x_xyy_xxxy_1[i] * fe_0 + ta1_x_xyy_xxxyz_0[i] * pa_z[i] - ta1_x_xyy_xxxyz_1[i] * pc_z[i];

        ta1_x_xyyz_xxxzz_0[i] =
            ta1_x_xz_xxxzz_0[i] * fe_0 - ta1_x_xz_xxxzz_1[i] * fe_0 + ta1_x_xyz_xxxzz_0[i] * pa_y[i] - ta1_x_xyz_xxxzz_1[i] * pc_y[i];

        ta1_x_xyyz_xxyyy_0[i] = ta1_x_xyy_xxyyy_0[i] * pa_z[i] - ta1_x_xyy_xxyyy_1[i] * pc_z[i];

        ta1_x_xyyz_xxyyz_0[i] =
            ta1_x_xyy_xxyy_0[i] * fe_0 - ta1_x_xyy_xxyy_1[i] * fe_0 + ta1_x_xyy_xxyyz_0[i] * pa_z[i] - ta1_x_xyy_xxyyz_1[i] * pc_z[i];

        ta1_x_xyyz_xxyzz_0[i] =
            2.0 * ta1_x_xyy_xxyz_0[i] * fe_0 - 2.0 * ta1_x_xyy_xxyz_1[i] * fe_0 + ta1_x_xyy_xxyzz_0[i] * pa_z[i] - ta1_x_xyy_xxyzz_1[i] * pc_z[i];

        ta1_x_xyyz_xxzzz_0[i] =
            ta1_x_xz_xxzzz_0[i] * fe_0 - ta1_x_xz_xxzzz_1[i] * fe_0 + ta1_x_xyz_xxzzz_0[i] * pa_y[i] - ta1_x_xyz_xxzzz_1[i] * pc_y[i];

        ta1_x_xyyz_xyyyy_0[i] = ta1_x_xyy_xyyyy_0[i] * pa_z[i] - ta1_x_xyy_xyyyy_1[i] * pc_z[i];

        ta1_x_xyyz_xyyyz_0[i] =
            ta1_x_xyy_xyyy_0[i] * fe_0 - ta1_x_xyy_xyyy_1[i] * fe_0 + ta1_x_xyy_xyyyz_0[i] * pa_z[i] - ta1_x_xyy_xyyyz_1[i] * pc_z[i];

        ta1_x_xyyz_xyyzz_0[i] =
            2.0 * ta1_x_xyy_xyyz_0[i] * fe_0 - 2.0 * ta1_x_xyy_xyyz_1[i] * fe_0 + ta1_x_xyy_xyyzz_0[i] * pa_z[i] - ta1_x_xyy_xyyzz_1[i] * pc_z[i];

        ta1_x_xyyz_xyzzz_0[i] =
            3.0 * ta1_x_xyy_xyzz_0[i] * fe_0 - 3.0 * ta1_x_xyy_xyzz_1[i] * fe_0 + ta1_x_xyy_xyzzz_0[i] * pa_z[i] - ta1_x_xyy_xyzzz_1[i] * pc_z[i];

        ta1_x_xyyz_xzzzz_0[i] =
            ta1_x_xz_xzzzz_0[i] * fe_0 - ta1_x_xz_xzzzz_1[i] * fe_0 + ta1_x_xyz_xzzzz_0[i] * pa_y[i] - ta1_x_xyz_xzzzz_1[i] * pc_y[i];

        ta1_x_xyyz_yyyyy_0[i] = ta1_x_xyy_yyyyy_0[i] * pa_z[i] - ta1_x_xyy_yyyyy_1[i] * pc_z[i];

        ta1_x_xyyz_yyyyz_0[i] = ta_yyz_yyyyz_1[i] + ta1_x_yyz_yyyyz_0[i] * pa_x[i] - ta1_x_yyz_yyyyz_1[i] * pc_x[i];

        ta1_x_xyyz_yyyzz_0[i] = ta_yyz_yyyzz_1[i] + ta1_x_yyz_yyyzz_0[i] * pa_x[i] - ta1_x_yyz_yyyzz_1[i] * pc_x[i];

        ta1_x_xyyz_yyzzz_0[i] = ta_yyz_yyzzz_1[i] + ta1_x_yyz_yyzzz_0[i] * pa_x[i] - ta1_x_yyz_yyzzz_1[i] * pc_x[i];

        ta1_x_xyyz_yzzzz_0[i] = ta_yyz_yzzzz_1[i] + ta1_x_yyz_yzzzz_0[i] * pa_x[i] - ta1_x_yyz_yzzzz_1[i] * pc_x[i];

        ta1_x_xyyz_zzzzz_0[i] = ta_yyz_zzzzz_1[i] + ta1_x_yyz_zzzzz_0[i] * pa_x[i] - ta1_x_yyz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 168-189 components of targeted buffer : GH

    auto ta1_x_xyzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 168);

    auto ta1_x_xyzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 169);

    auto ta1_x_xyzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 170);

    auto ta1_x_xyzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 171);

    auto ta1_x_xyzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 172);

    auto ta1_x_xyzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 173);

    auto ta1_x_xyzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 174);

    auto ta1_x_xyzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 175);

    auto ta1_x_xyzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 176);

    auto ta1_x_xyzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 177);

    auto ta1_x_xyzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 178);

    auto ta1_x_xyzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 179);

    auto ta1_x_xyzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 180);

    auto ta1_x_xyzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 181);

    auto ta1_x_xyzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 182);

    auto ta1_x_xyzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 183);

    auto ta1_x_xyzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 184);

    auto ta1_x_xyzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 185);

    auto ta1_x_xyzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 186);

    auto ta1_x_xyzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 187);

    auto ta1_x_xyzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 188);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_x_xyzz_xxxxx_0, \
                             ta1_x_xyzz_xxxxy_0, \
                             ta1_x_xyzz_xxxxz_0, \
                             ta1_x_xyzz_xxxyy_0, \
                             ta1_x_xyzz_xxxyz_0, \
                             ta1_x_xyzz_xxxzz_0, \
                             ta1_x_xyzz_xxyyy_0, \
                             ta1_x_xyzz_xxyyz_0, \
                             ta1_x_xyzz_xxyzz_0, \
                             ta1_x_xyzz_xxzzz_0, \
                             ta1_x_xyzz_xyyyy_0, \
                             ta1_x_xyzz_xyyyz_0, \
                             ta1_x_xyzz_xyyzz_0, \
                             ta1_x_xyzz_xyzzz_0, \
                             ta1_x_xyzz_xzzzz_0, \
                             ta1_x_xyzz_yyyyy_0, \
                             ta1_x_xyzz_yyyyz_0, \
                             ta1_x_xyzz_yyyzz_0, \
                             ta1_x_xyzz_yyzzz_0, \
                             ta1_x_xyzz_yzzzz_0, \
                             ta1_x_xyzz_zzzzz_0, \
                             ta1_x_xzz_xxxx_0,   \
                             ta1_x_xzz_xxxx_1,   \
                             ta1_x_xzz_xxxxx_0,  \
                             ta1_x_xzz_xxxxx_1,  \
                             ta1_x_xzz_xxxxy_0,  \
                             ta1_x_xzz_xxxxy_1,  \
                             ta1_x_xzz_xxxxz_0,  \
                             ta1_x_xzz_xxxxz_1,  \
                             ta1_x_xzz_xxxy_0,   \
                             ta1_x_xzz_xxxy_1,   \
                             ta1_x_xzz_xxxyy_0,  \
                             ta1_x_xzz_xxxyy_1,  \
                             ta1_x_xzz_xxxyz_0,  \
                             ta1_x_xzz_xxxyz_1,  \
                             ta1_x_xzz_xxxz_0,   \
                             ta1_x_xzz_xxxz_1,   \
                             ta1_x_xzz_xxxzz_0,  \
                             ta1_x_xzz_xxxzz_1,  \
                             ta1_x_xzz_xxyy_0,   \
                             ta1_x_xzz_xxyy_1,   \
                             ta1_x_xzz_xxyyy_0,  \
                             ta1_x_xzz_xxyyy_1,  \
                             ta1_x_xzz_xxyyz_0,  \
                             ta1_x_xzz_xxyyz_1,  \
                             ta1_x_xzz_xxyz_0,   \
                             ta1_x_xzz_xxyz_1,   \
                             ta1_x_xzz_xxyzz_0,  \
                             ta1_x_xzz_xxyzz_1,  \
                             ta1_x_xzz_xxzz_0,   \
                             ta1_x_xzz_xxzz_1,   \
                             ta1_x_xzz_xxzzz_0,  \
                             ta1_x_xzz_xxzzz_1,  \
                             ta1_x_xzz_xyyy_0,   \
                             ta1_x_xzz_xyyy_1,   \
                             ta1_x_xzz_xyyyy_0,  \
                             ta1_x_xzz_xyyyy_1,  \
                             ta1_x_xzz_xyyyz_0,  \
                             ta1_x_xzz_xyyyz_1,  \
                             ta1_x_xzz_xyyz_0,   \
                             ta1_x_xzz_xyyz_1,   \
                             ta1_x_xzz_xyyzz_0,  \
                             ta1_x_xzz_xyyzz_1,  \
                             ta1_x_xzz_xyzz_0,   \
                             ta1_x_xzz_xyzz_1,   \
                             ta1_x_xzz_xyzzz_0,  \
                             ta1_x_xzz_xyzzz_1,  \
                             ta1_x_xzz_xzzz_0,   \
                             ta1_x_xzz_xzzz_1,   \
                             ta1_x_xzz_xzzzz_0,  \
                             ta1_x_xzz_xzzzz_1,  \
                             ta1_x_xzz_zzzzz_0,  \
                             ta1_x_xzz_zzzzz_1,  \
                             ta1_x_yzz_yyyyy_0,  \
                             ta1_x_yzz_yyyyy_1,  \
                             ta1_x_yzz_yyyyz_0,  \
                             ta1_x_yzz_yyyyz_1,  \
                             ta1_x_yzz_yyyzz_0,  \
                             ta1_x_yzz_yyyzz_1,  \
                             ta1_x_yzz_yyzzz_0,  \
                             ta1_x_yzz_yyzzz_1,  \
                             ta1_x_yzz_yzzzz_0,  \
                             ta1_x_yzz_yzzzz_1,  \
                             ta_yzz_yyyyy_1,     \
                             ta_yzz_yyyyz_1,     \
                             ta_yzz_yyyzz_1,     \
                             ta_yzz_yyzzz_1,     \
                             ta_yzz_yzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xyzz_xxxxx_0[i] = ta1_x_xzz_xxxxx_0[i] * pa_y[i] - ta1_x_xzz_xxxxx_1[i] * pc_y[i];

        ta1_x_xyzz_xxxxy_0[i] =
            ta1_x_xzz_xxxx_0[i] * fe_0 - ta1_x_xzz_xxxx_1[i] * fe_0 + ta1_x_xzz_xxxxy_0[i] * pa_y[i] - ta1_x_xzz_xxxxy_1[i] * pc_y[i];

        ta1_x_xyzz_xxxxz_0[i] = ta1_x_xzz_xxxxz_0[i] * pa_y[i] - ta1_x_xzz_xxxxz_1[i] * pc_y[i];

        ta1_x_xyzz_xxxyy_0[i] =
            2.0 * ta1_x_xzz_xxxy_0[i] * fe_0 - 2.0 * ta1_x_xzz_xxxy_1[i] * fe_0 + ta1_x_xzz_xxxyy_0[i] * pa_y[i] - ta1_x_xzz_xxxyy_1[i] * pc_y[i];

        ta1_x_xyzz_xxxyz_0[i] =
            ta1_x_xzz_xxxz_0[i] * fe_0 - ta1_x_xzz_xxxz_1[i] * fe_0 + ta1_x_xzz_xxxyz_0[i] * pa_y[i] - ta1_x_xzz_xxxyz_1[i] * pc_y[i];

        ta1_x_xyzz_xxxzz_0[i] = ta1_x_xzz_xxxzz_0[i] * pa_y[i] - ta1_x_xzz_xxxzz_1[i] * pc_y[i];

        ta1_x_xyzz_xxyyy_0[i] =
            3.0 * ta1_x_xzz_xxyy_0[i] * fe_0 - 3.0 * ta1_x_xzz_xxyy_1[i] * fe_0 + ta1_x_xzz_xxyyy_0[i] * pa_y[i] - ta1_x_xzz_xxyyy_1[i] * pc_y[i];

        ta1_x_xyzz_xxyyz_0[i] =
            2.0 * ta1_x_xzz_xxyz_0[i] * fe_0 - 2.0 * ta1_x_xzz_xxyz_1[i] * fe_0 + ta1_x_xzz_xxyyz_0[i] * pa_y[i] - ta1_x_xzz_xxyyz_1[i] * pc_y[i];

        ta1_x_xyzz_xxyzz_0[i] =
            ta1_x_xzz_xxzz_0[i] * fe_0 - ta1_x_xzz_xxzz_1[i] * fe_0 + ta1_x_xzz_xxyzz_0[i] * pa_y[i] - ta1_x_xzz_xxyzz_1[i] * pc_y[i];

        ta1_x_xyzz_xxzzz_0[i] = ta1_x_xzz_xxzzz_0[i] * pa_y[i] - ta1_x_xzz_xxzzz_1[i] * pc_y[i];

        ta1_x_xyzz_xyyyy_0[i] =
            4.0 * ta1_x_xzz_xyyy_0[i] * fe_0 - 4.0 * ta1_x_xzz_xyyy_1[i] * fe_0 + ta1_x_xzz_xyyyy_0[i] * pa_y[i] - ta1_x_xzz_xyyyy_1[i] * pc_y[i];

        ta1_x_xyzz_xyyyz_0[i] =
            3.0 * ta1_x_xzz_xyyz_0[i] * fe_0 - 3.0 * ta1_x_xzz_xyyz_1[i] * fe_0 + ta1_x_xzz_xyyyz_0[i] * pa_y[i] - ta1_x_xzz_xyyyz_1[i] * pc_y[i];

        ta1_x_xyzz_xyyzz_0[i] =
            2.0 * ta1_x_xzz_xyzz_0[i] * fe_0 - 2.0 * ta1_x_xzz_xyzz_1[i] * fe_0 + ta1_x_xzz_xyyzz_0[i] * pa_y[i] - ta1_x_xzz_xyyzz_1[i] * pc_y[i];

        ta1_x_xyzz_xyzzz_0[i] =
            ta1_x_xzz_xzzz_0[i] * fe_0 - ta1_x_xzz_xzzz_1[i] * fe_0 + ta1_x_xzz_xyzzz_0[i] * pa_y[i] - ta1_x_xzz_xyzzz_1[i] * pc_y[i];

        ta1_x_xyzz_xzzzz_0[i] = ta1_x_xzz_xzzzz_0[i] * pa_y[i] - ta1_x_xzz_xzzzz_1[i] * pc_y[i];

        ta1_x_xyzz_yyyyy_0[i] = ta_yzz_yyyyy_1[i] + ta1_x_yzz_yyyyy_0[i] * pa_x[i] - ta1_x_yzz_yyyyy_1[i] * pc_x[i];

        ta1_x_xyzz_yyyyz_0[i] = ta_yzz_yyyyz_1[i] + ta1_x_yzz_yyyyz_0[i] * pa_x[i] - ta1_x_yzz_yyyyz_1[i] * pc_x[i];

        ta1_x_xyzz_yyyzz_0[i] = ta_yzz_yyyzz_1[i] + ta1_x_yzz_yyyzz_0[i] * pa_x[i] - ta1_x_yzz_yyyzz_1[i] * pc_x[i];

        ta1_x_xyzz_yyzzz_0[i] = ta_yzz_yyzzz_1[i] + ta1_x_yzz_yyzzz_0[i] * pa_x[i] - ta1_x_yzz_yyzzz_1[i] * pc_x[i];

        ta1_x_xyzz_yzzzz_0[i] = ta_yzz_yzzzz_1[i] + ta1_x_yzz_yzzzz_0[i] * pa_x[i] - ta1_x_yzz_yzzzz_1[i] * pc_x[i];

        ta1_x_xyzz_zzzzz_0[i] = ta1_x_xzz_zzzzz_0[i] * pa_y[i] - ta1_x_xzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 189-210 components of targeted buffer : GH

    auto ta1_x_xzzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 189);

    auto ta1_x_xzzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 190);

    auto ta1_x_xzzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 191);

    auto ta1_x_xzzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 192);

    auto ta1_x_xzzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 193);

    auto ta1_x_xzzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 194);

    auto ta1_x_xzzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 195);

    auto ta1_x_xzzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 196);

    auto ta1_x_xzzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 197);

    auto ta1_x_xzzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 198);

    auto ta1_x_xzzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 199);

    auto ta1_x_xzzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 200);

    auto ta1_x_xzzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 201);

    auto ta1_x_xzzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 202);

    auto ta1_x_xzzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 203);

    auto ta1_x_xzzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 204);

    auto ta1_x_xzzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 205);

    auto ta1_x_xzzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 206);

    auto ta1_x_xzzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 207);

    auto ta1_x_xzzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 208);

    auto ta1_x_xzzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 209);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_x_xz_xxxxx_0,   \
                             ta1_x_xz_xxxxx_1,   \
                             ta1_x_xz_xxxxy_0,   \
                             ta1_x_xz_xxxxy_1,   \
                             ta1_x_xz_xxxyy_0,   \
                             ta1_x_xz_xxxyy_1,   \
                             ta1_x_xz_xxyyy_0,   \
                             ta1_x_xz_xxyyy_1,   \
                             ta1_x_xz_xyyyy_0,   \
                             ta1_x_xz_xyyyy_1,   \
                             ta1_x_xzz_xxxxx_0,  \
                             ta1_x_xzz_xxxxx_1,  \
                             ta1_x_xzz_xxxxy_0,  \
                             ta1_x_xzz_xxxxy_1,  \
                             ta1_x_xzz_xxxyy_0,  \
                             ta1_x_xzz_xxxyy_1,  \
                             ta1_x_xzz_xxyyy_0,  \
                             ta1_x_xzz_xxyyy_1,  \
                             ta1_x_xzz_xyyyy_0,  \
                             ta1_x_xzz_xyyyy_1,  \
                             ta1_x_xzzz_xxxxx_0, \
                             ta1_x_xzzz_xxxxy_0, \
                             ta1_x_xzzz_xxxxz_0, \
                             ta1_x_xzzz_xxxyy_0, \
                             ta1_x_xzzz_xxxyz_0, \
                             ta1_x_xzzz_xxxzz_0, \
                             ta1_x_xzzz_xxyyy_0, \
                             ta1_x_xzzz_xxyyz_0, \
                             ta1_x_xzzz_xxyzz_0, \
                             ta1_x_xzzz_xxzzz_0, \
                             ta1_x_xzzz_xyyyy_0, \
                             ta1_x_xzzz_xyyyz_0, \
                             ta1_x_xzzz_xyyzz_0, \
                             ta1_x_xzzz_xyzzz_0, \
                             ta1_x_xzzz_xzzzz_0, \
                             ta1_x_xzzz_yyyyy_0, \
                             ta1_x_xzzz_yyyyz_0, \
                             ta1_x_xzzz_yyyzz_0, \
                             ta1_x_xzzz_yyzzz_0, \
                             ta1_x_xzzz_yzzzz_0, \
                             ta1_x_xzzz_zzzzz_0, \
                             ta1_x_zzz_xxxxz_0,  \
                             ta1_x_zzz_xxxxz_1,  \
                             ta1_x_zzz_xxxyz_0,  \
                             ta1_x_zzz_xxxyz_1,  \
                             ta1_x_zzz_xxxz_0,   \
                             ta1_x_zzz_xxxz_1,   \
                             ta1_x_zzz_xxxzz_0,  \
                             ta1_x_zzz_xxxzz_1,  \
                             ta1_x_zzz_xxyyz_0,  \
                             ta1_x_zzz_xxyyz_1,  \
                             ta1_x_zzz_xxyz_0,   \
                             ta1_x_zzz_xxyz_1,   \
                             ta1_x_zzz_xxyzz_0,  \
                             ta1_x_zzz_xxyzz_1,  \
                             ta1_x_zzz_xxzz_0,   \
                             ta1_x_zzz_xxzz_1,   \
                             ta1_x_zzz_xxzzz_0,  \
                             ta1_x_zzz_xxzzz_1,  \
                             ta1_x_zzz_xyyyz_0,  \
                             ta1_x_zzz_xyyyz_1,  \
                             ta1_x_zzz_xyyz_0,   \
                             ta1_x_zzz_xyyz_1,   \
                             ta1_x_zzz_xyyzz_0,  \
                             ta1_x_zzz_xyyzz_1,  \
                             ta1_x_zzz_xyzz_0,   \
                             ta1_x_zzz_xyzz_1,   \
                             ta1_x_zzz_xyzzz_0,  \
                             ta1_x_zzz_xyzzz_1,  \
                             ta1_x_zzz_xzzz_0,   \
                             ta1_x_zzz_xzzz_1,   \
                             ta1_x_zzz_xzzzz_0,  \
                             ta1_x_zzz_xzzzz_1,  \
                             ta1_x_zzz_yyyyy_0,  \
                             ta1_x_zzz_yyyyy_1,  \
                             ta1_x_zzz_yyyyz_0,  \
                             ta1_x_zzz_yyyyz_1,  \
                             ta1_x_zzz_yyyz_0,   \
                             ta1_x_zzz_yyyz_1,   \
                             ta1_x_zzz_yyyzz_0,  \
                             ta1_x_zzz_yyyzz_1,  \
                             ta1_x_zzz_yyzz_0,   \
                             ta1_x_zzz_yyzz_1,   \
                             ta1_x_zzz_yyzzz_0,  \
                             ta1_x_zzz_yyzzz_1,  \
                             ta1_x_zzz_yzzz_0,   \
                             ta1_x_zzz_yzzz_1,   \
                             ta1_x_zzz_yzzzz_0,  \
                             ta1_x_zzz_yzzzz_1,  \
                             ta1_x_zzz_zzzz_0,   \
                             ta1_x_zzz_zzzz_1,   \
                             ta1_x_zzz_zzzzz_0,  \
                             ta1_x_zzz_zzzzz_1,  \
                             ta_zzz_xxxxz_1,     \
                             ta_zzz_xxxyz_1,     \
                             ta_zzz_xxxzz_1,     \
                             ta_zzz_xxyyz_1,     \
                             ta_zzz_xxyzz_1,     \
                             ta_zzz_xxzzz_1,     \
                             ta_zzz_xyyyz_1,     \
                             ta_zzz_xyyzz_1,     \
                             ta_zzz_xyzzz_1,     \
                             ta_zzz_xzzzz_1,     \
                             ta_zzz_yyyyy_1,     \
                             ta_zzz_yyyyz_1,     \
                             ta_zzz_yyyzz_1,     \
                             ta_zzz_yyzzz_1,     \
                             ta_zzz_yzzzz_1,     \
                             ta_zzz_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_xzzz_xxxxx_0[i] =
            2.0 * ta1_x_xz_xxxxx_0[i] * fe_0 - 2.0 * ta1_x_xz_xxxxx_1[i] * fe_0 + ta1_x_xzz_xxxxx_0[i] * pa_z[i] - ta1_x_xzz_xxxxx_1[i] * pc_z[i];

        ta1_x_xzzz_xxxxy_0[i] =
            2.0 * ta1_x_xz_xxxxy_0[i] * fe_0 - 2.0 * ta1_x_xz_xxxxy_1[i] * fe_0 + ta1_x_xzz_xxxxy_0[i] * pa_z[i] - ta1_x_xzz_xxxxy_1[i] * pc_z[i];

        ta1_x_xzzz_xxxxz_0[i] = 4.0 * ta1_x_zzz_xxxz_0[i] * fe_0 - 4.0 * ta1_x_zzz_xxxz_1[i] * fe_0 + ta_zzz_xxxxz_1[i] +
                                ta1_x_zzz_xxxxz_0[i] * pa_x[i] - ta1_x_zzz_xxxxz_1[i] * pc_x[i];

        ta1_x_xzzz_xxxyy_0[i] =
            2.0 * ta1_x_xz_xxxyy_0[i] * fe_0 - 2.0 * ta1_x_xz_xxxyy_1[i] * fe_0 + ta1_x_xzz_xxxyy_0[i] * pa_z[i] - ta1_x_xzz_xxxyy_1[i] * pc_z[i];

        ta1_x_xzzz_xxxyz_0[i] = 3.0 * ta1_x_zzz_xxyz_0[i] * fe_0 - 3.0 * ta1_x_zzz_xxyz_1[i] * fe_0 + ta_zzz_xxxyz_1[i] +
                                ta1_x_zzz_xxxyz_0[i] * pa_x[i] - ta1_x_zzz_xxxyz_1[i] * pc_x[i];

        ta1_x_xzzz_xxxzz_0[i] = 3.0 * ta1_x_zzz_xxzz_0[i] * fe_0 - 3.0 * ta1_x_zzz_xxzz_1[i] * fe_0 + ta_zzz_xxxzz_1[i] +
                                ta1_x_zzz_xxxzz_0[i] * pa_x[i] - ta1_x_zzz_xxxzz_1[i] * pc_x[i];

        ta1_x_xzzz_xxyyy_0[i] =
            2.0 * ta1_x_xz_xxyyy_0[i] * fe_0 - 2.0 * ta1_x_xz_xxyyy_1[i] * fe_0 + ta1_x_xzz_xxyyy_0[i] * pa_z[i] - ta1_x_xzz_xxyyy_1[i] * pc_z[i];

        ta1_x_xzzz_xxyyz_0[i] = 2.0 * ta1_x_zzz_xyyz_0[i] * fe_0 - 2.0 * ta1_x_zzz_xyyz_1[i] * fe_0 + ta_zzz_xxyyz_1[i] +
                                ta1_x_zzz_xxyyz_0[i] * pa_x[i] - ta1_x_zzz_xxyyz_1[i] * pc_x[i];

        ta1_x_xzzz_xxyzz_0[i] = 2.0 * ta1_x_zzz_xyzz_0[i] * fe_0 - 2.0 * ta1_x_zzz_xyzz_1[i] * fe_0 + ta_zzz_xxyzz_1[i] +
                                ta1_x_zzz_xxyzz_0[i] * pa_x[i] - ta1_x_zzz_xxyzz_1[i] * pc_x[i];

        ta1_x_xzzz_xxzzz_0[i] = 2.0 * ta1_x_zzz_xzzz_0[i] * fe_0 - 2.0 * ta1_x_zzz_xzzz_1[i] * fe_0 + ta_zzz_xxzzz_1[i] +
                                ta1_x_zzz_xxzzz_0[i] * pa_x[i] - ta1_x_zzz_xxzzz_1[i] * pc_x[i];

        ta1_x_xzzz_xyyyy_0[i] =
            2.0 * ta1_x_xz_xyyyy_0[i] * fe_0 - 2.0 * ta1_x_xz_xyyyy_1[i] * fe_0 + ta1_x_xzz_xyyyy_0[i] * pa_z[i] - ta1_x_xzz_xyyyy_1[i] * pc_z[i];

        ta1_x_xzzz_xyyyz_0[i] = ta1_x_zzz_yyyz_0[i] * fe_0 - ta1_x_zzz_yyyz_1[i] * fe_0 + ta_zzz_xyyyz_1[i] + ta1_x_zzz_xyyyz_0[i] * pa_x[i] -
                                ta1_x_zzz_xyyyz_1[i] * pc_x[i];

        ta1_x_xzzz_xyyzz_0[i] = ta1_x_zzz_yyzz_0[i] * fe_0 - ta1_x_zzz_yyzz_1[i] * fe_0 + ta_zzz_xyyzz_1[i] + ta1_x_zzz_xyyzz_0[i] * pa_x[i] -
                                ta1_x_zzz_xyyzz_1[i] * pc_x[i];

        ta1_x_xzzz_xyzzz_0[i] = ta1_x_zzz_yzzz_0[i] * fe_0 - ta1_x_zzz_yzzz_1[i] * fe_0 + ta_zzz_xyzzz_1[i] + ta1_x_zzz_xyzzz_0[i] * pa_x[i] -
                                ta1_x_zzz_xyzzz_1[i] * pc_x[i];

        ta1_x_xzzz_xzzzz_0[i] = ta1_x_zzz_zzzz_0[i] * fe_0 - ta1_x_zzz_zzzz_1[i] * fe_0 + ta_zzz_xzzzz_1[i] + ta1_x_zzz_xzzzz_0[i] * pa_x[i] -
                                ta1_x_zzz_xzzzz_1[i] * pc_x[i];

        ta1_x_xzzz_yyyyy_0[i] = ta_zzz_yyyyy_1[i] + ta1_x_zzz_yyyyy_0[i] * pa_x[i] - ta1_x_zzz_yyyyy_1[i] * pc_x[i];

        ta1_x_xzzz_yyyyz_0[i] = ta_zzz_yyyyz_1[i] + ta1_x_zzz_yyyyz_0[i] * pa_x[i] - ta1_x_zzz_yyyyz_1[i] * pc_x[i];

        ta1_x_xzzz_yyyzz_0[i] = ta_zzz_yyyzz_1[i] + ta1_x_zzz_yyyzz_0[i] * pa_x[i] - ta1_x_zzz_yyyzz_1[i] * pc_x[i];

        ta1_x_xzzz_yyzzz_0[i] = ta_zzz_yyzzz_1[i] + ta1_x_zzz_yyzzz_0[i] * pa_x[i] - ta1_x_zzz_yyzzz_1[i] * pc_x[i];

        ta1_x_xzzz_yzzzz_0[i] = ta_zzz_yzzzz_1[i] + ta1_x_zzz_yzzzz_0[i] * pa_x[i] - ta1_x_zzz_yzzzz_1[i] * pc_x[i];

        ta1_x_xzzz_zzzzz_0[i] = ta_zzz_zzzzz_1[i] + ta1_x_zzz_zzzzz_0[i] * pa_x[i] - ta1_x_zzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 210-231 components of targeted buffer : GH

    auto ta1_x_yyyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 210);

    auto ta1_x_yyyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 211);

    auto ta1_x_yyyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 212);

    auto ta1_x_yyyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 213);

    auto ta1_x_yyyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 214);

    auto ta1_x_yyyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 215);

    auto ta1_x_yyyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 216);

    auto ta1_x_yyyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 217);

    auto ta1_x_yyyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 218);

    auto ta1_x_yyyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 219);

    auto ta1_x_yyyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 220);

    auto ta1_x_yyyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 221);

    auto ta1_x_yyyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 222);

    auto ta1_x_yyyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 223);

    auto ta1_x_yyyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 224);

    auto ta1_x_yyyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 225);

    auto ta1_x_yyyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 226);

    auto ta1_x_yyyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 227);

    auto ta1_x_yyyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 228);

    auto ta1_x_yyyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 229);

    auto ta1_x_yyyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 230);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_x_yy_xxxxx_0,   \
                             ta1_x_yy_xxxxx_1,   \
                             ta1_x_yy_xxxxy_0,   \
                             ta1_x_yy_xxxxy_1,   \
                             ta1_x_yy_xxxxz_0,   \
                             ta1_x_yy_xxxxz_1,   \
                             ta1_x_yy_xxxyy_0,   \
                             ta1_x_yy_xxxyy_1,   \
                             ta1_x_yy_xxxyz_0,   \
                             ta1_x_yy_xxxyz_1,   \
                             ta1_x_yy_xxxzz_0,   \
                             ta1_x_yy_xxxzz_1,   \
                             ta1_x_yy_xxyyy_0,   \
                             ta1_x_yy_xxyyy_1,   \
                             ta1_x_yy_xxyyz_0,   \
                             ta1_x_yy_xxyyz_1,   \
                             ta1_x_yy_xxyzz_0,   \
                             ta1_x_yy_xxyzz_1,   \
                             ta1_x_yy_xxzzz_0,   \
                             ta1_x_yy_xxzzz_1,   \
                             ta1_x_yy_xyyyy_0,   \
                             ta1_x_yy_xyyyy_1,   \
                             ta1_x_yy_xyyyz_0,   \
                             ta1_x_yy_xyyyz_1,   \
                             ta1_x_yy_xyyzz_0,   \
                             ta1_x_yy_xyyzz_1,   \
                             ta1_x_yy_xyzzz_0,   \
                             ta1_x_yy_xyzzz_1,   \
                             ta1_x_yy_xzzzz_0,   \
                             ta1_x_yy_xzzzz_1,   \
                             ta1_x_yy_yyyyy_0,   \
                             ta1_x_yy_yyyyy_1,   \
                             ta1_x_yy_yyyyz_0,   \
                             ta1_x_yy_yyyyz_1,   \
                             ta1_x_yy_yyyzz_0,   \
                             ta1_x_yy_yyyzz_1,   \
                             ta1_x_yy_yyzzz_0,   \
                             ta1_x_yy_yyzzz_1,   \
                             ta1_x_yy_yzzzz_0,   \
                             ta1_x_yy_yzzzz_1,   \
                             ta1_x_yy_zzzzz_0,   \
                             ta1_x_yy_zzzzz_1,   \
                             ta1_x_yyy_xxxx_0,   \
                             ta1_x_yyy_xxxx_1,   \
                             ta1_x_yyy_xxxxx_0,  \
                             ta1_x_yyy_xxxxx_1,  \
                             ta1_x_yyy_xxxxy_0,  \
                             ta1_x_yyy_xxxxy_1,  \
                             ta1_x_yyy_xxxxz_0,  \
                             ta1_x_yyy_xxxxz_1,  \
                             ta1_x_yyy_xxxy_0,   \
                             ta1_x_yyy_xxxy_1,   \
                             ta1_x_yyy_xxxyy_0,  \
                             ta1_x_yyy_xxxyy_1,  \
                             ta1_x_yyy_xxxyz_0,  \
                             ta1_x_yyy_xxxyz_1,  \
                             ta1_x_yyy_xxxz_0,   \
                             ta1_x_yyy_xxxz_1,   \
                             ta1_x_yyy_xxxzz_0,  \
                             ta1_x_yyy_xxxzz_1,  \
                             ta1_x_yyy_xxyy_0,   \
                             ta1_x_yyy_xxyy_1,   \
                             ta1_x_yyy_xxyyy_0,  \
                             ta1_x_yyy_xxyyy_1,  \
                             ta1_x_yyy_xxyyz_0,  \
                             ta1_x_yyy_xxyyz_1,  \
                             ta1_x_yyy_xxyz_0,   \
                             ta1_x_yyy_xxyz_1,   \
                             ta1_x_yyy_xxyzz_0,  \
                             ta1_x_yyy_xxyzz_1,  \
                             ta1_x_yyy_xxzz_0,   \
                             ta1_x_yyy_xxzz_1,   \
                             ta1_x_yyy_xxzzz_0,  \
                             ta1_x_yyy_xxzzz_1,  \
                             ta1_x_yyy_xyyy_0,   \
                             ta1_x_yyy_xyyy_1,   \
                             ta1_x_yyy_xyyyy_0,  \
                             ta1_x_yyy_xyyyy_1,  \
                             ta1_x_yyy_xyyyz_0,  \
                             ta1_x_yyy_xyyyz_1,  \
                             ta1_x_yyy_xyyz_0,   \
                             ta1_x_yyy_xyyz_1,   \
                             ta1_x_yyy_xyyzz_0,  \
                             ta1_x_yyy_xyyzz_1,  \
                             ta1_x_yyy_xyzz_0,   \
                             ta1_x_yyy_xyzz_1,   \
                             ta1_x_yyy_xyzzz_0,  \
                             ta1_x_yyy_xyzzz_1,  \
                             ta1_x_yyy_xzzz_0,   \
                             ta1_x_yyy_xzzz_1,   \
                             ta1_x_yyy_xzzzz_0,  \
                             ta1_x_yyy_xzzzz_1,  \
                             ta1_x_yyy_yyyy_0,   \
                             ta1_x_yyy_yyyy_1,   \
                             ta1_x_yyy_yyyyy_0,  \
                             ta1_x_yyy_yyyyy_1,  \
                             ta1_x_yyy_yyyyz_0,  \
                             ta1_x_yyy_yyyyz_1,  \
                             ta1_x_yyy_yyyz_0,   \
                             ta1_x_yyy_yyyz_1,   \
                             ta1_x_yyy_yyyzz_0,  \
                             ta1_x_yyy_yyyzz_1,  \
                             ta1_x_yyy_yyzz_0,   \
                             ta1_x_yyy_yyzz_1,   \
                             ta1_x_yyy_yyzzz_0,  \
                             ta1_x_yyy_yyzzz_1,  \
                             ta1_x_yyy_yzzz_0,   \
                             ta1_x_yyy_yzzz_1,   \
                             ta1_x_yyy_yzzzz_0,  \
                             ta1_x_yyy_yzzzz_1,  \
                             ta1_x_yyy_zzzz_0,   \
                             ta1_x_yyy_zzzz_1,   \
                             ta1_x_yyy_zzzzz_0,  \
                             ta1_x_yyy_zzzzz_1,  \
                             ta1_x_yyyy_xxxxx_0, \
                             ta1_x_yyyy_xxxxy_0, \
                             ta1_x_yyyy_xxxxz_0, \
                             ta1_x_yyyy_xxxyy_0, \
                             ta1_x_yyyy_xxxyz_0, \
                             ta1_x_yyyy_xxxzz_0, \
                             ta1_x_yyyy_xxyyy_0, \
                             ta1_x_yyyy_xxyyz_0, \
                             ta1_x_yyyy_xxyzz_0, \
                             ta1_x_yyyy_xxzzz_0, \
                             ta1_x_yyyy_xyyyy_0, \
                             ta1_x_yyyy_xyyyz_0, \
                             ta1_x_yyyy_xyyzz_0, \
                             ta1_x_yyyy_xyzzz_0, \
                             ta1_x_yyyy_xzzzz_0, \
                             ta1_x_yyyy_yyyyy_0, \
                             ta1_x_yyyy_yyyyz_0, \
                             ta1_x_yyyy_yyyzz_0, \
                             ta1_x_yyyy_yyzzz_0, \
                             ta1_x_yyyy_yzzzz_0, \
                             ta1_x_yyyy_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyy_xxxxx_0[i] =
            3.0 * ta1_x_yy_xxxxx_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxxx_1[i] * fe_0 + ta1_x_yyy_xxxxx_0[i] * pa_y[i] - ta1_x_yyy_xxxxx_1[i] * pc_y[i];

        ta1_x_yyyy_xxxxy_0[i] = 3.0 * ta1_x_yy_xxxxy_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxxy_1[i] * fe_0 + ta1_x_yyy_xxxx_0[i] * fe_0 -
                                ta1_x_yyy_xxxx_1[i] * fe_0 + ta1_x_yyy_xxxxy_0[i] * pa_y[i] - ta1_x_yyy_xxxxy_1[i] * pc_y[i];

        ta1_x_yyyy_xxxxz_0[i] =
            3.0 * ta1_x_yy_xxxxz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxxz_1[i] * fe_0 + ta1_x_yyy_xxxxz_0[i] * pa_y[i] - ta1_x_yyy_xxxxz_1[i] * pc_y[i];

        ta1_x_yyyy_xxxyy_0[i] = 3.0 * ta1_x_yy_xxxyy_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxyy_1[i] * fe_0 + 2.0 * ta1_x_yyy_xxxy_0[i] * fe_0 -
                                2.0 * ta1_x_yyy_xxxy_1[i] * fe_0 + ta1_x_yyy_xxxyy_0[i] * pa_y[i] - ta1_x_yyy_xxxyy_1[i] * pc_y[i];

        ta1_x_yyyy_xxxyz_0[i] = 3.0 * ta1_x_yy_xxxyz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxyz_1[i] * fe_0 + ta1_x_yyy_xxxz_0[i] * fe_0 -
                                ta1_x_yyy_xxxz_1[i] * fe_0 + ta1_x_yyy_xxxyz_0[i] * pa_y[i] - ta1_x_yyy_xxxyz_1[i] * pc_y[i];

        ta1_x_yyyy_xxxzz_0[i] =
            3.0 * ta1_x_yy_xxxzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxxzz_1[i] * fe_0 + ta1_x_yyy_xxxzz_0[i] * pa_y[i] - ta1_x_yyy_xxxzz_1[i] * pc_y[i];

        ta1_x_yyyy_xxyyy_0[i] = 3.0 * ta1_x_yy_xxyyy_0[i] * fe_0 - 3.0 * ta1_x_yy_xxyyy_1[i] * fe_0 + 3.0 * ta1_x_yyy_xxyy_0[i] * fe_0 -
                                3.0 * ta1_x_yyy_xxyy_1[i] * fe_0 + ta1_x_yyy_xxyyy_0[i] * pa_y[i] - ta1_x_yyy_xxyyy_1[i] * pc_y[i];

        ta1_x_yyyy_xxyyz_0[i] = 3.0 * ta1_x_yy_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxyyz_1[i] * fe_0 + 2.0 * ta1_x_yyy_xxyz_0[i] * fe_0 -
                                2.0 * ta1_x_yyy_xxyz_1[i] * fe_0 + ta1_x_yyy_xxyyz_0[i] * pa_y[i] - ta1_x_yyy_xxyyz_1[i] * pc_y[i];

        ta1_x_yyyy_xxyzz_0[i] = 3.0 * ta1_x_yy_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxyzz_1[i] * fe_0 + ta1_x_yyy_xxzz_0[i] * fe_0 -
                                ta1_x_yyy_xxzz_1[i] * fe_0 + ta1_x_yyy_xxyzz_0[i] * pa_y[i] - ta1_x_yyy_xxyzz_1[i] * pc_y[i];

        ta1_x_yyyy_xxzzz_0[i] =
            3.0 * ta1_x_yy_xxzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xxzzz_1[i] * fe_0 + ta1_x_yyy_xxzzz_0[i] * pa_y[i] - ta1_x_yyy_xxzzz_1[i] * pc_y[i];

        ta1_x_yyyy_xyyyy_0[i] = 3.0 * ta1_x_yy_xyyyy_0[i] * fe_0 - 3.0 * ta1_x_yy_xyyyy_1[i] * fe_0 + 4.0 * ta1_x_yyy_xyyy_0[i] * fe_0 -
                                4.0 * ta1_x_yyy_xyyy_1[i] * fe_0 + ta1_x_yyy_xyyyy_0[i] * pa_y[i] - ta1_x_yyy_xyyyy_1[i] * pc_y[i];

        ta1_x_yyyy_xyyyz_0[i] = 3.0 * ta1_x_yy_xyyyz_0[i] * fe_0 - 3.0 * ta1_x_yy_xyyyz_1[i] * fe_0 + 3.0 * ta1_x_yyy_xyyz_0[i] * fe_0 -
                                3.0 * ta1_x_yyy_xyyz_1[i] * fe_0 + ta1_x_yyy_xyyyz_0[i] * pa_y[i] - ta1_x_yyy_xyyyz_1[i] * pc_y[i];

        ta1_x_yyyy_xyyzz_0[i] = 3.0 * ta1_x_yy_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xyyzz_1[i] * fe_0 + 2.0 * ta1_x_yyy_xyzz_0[i] * fe_0 -
                                2.0 * ta1_x_yyy_xyzz_1[i] * fe_0 + ta1_x_yyy_xyyzz_0[i] * pa_y[i] - ta1_x_yyy_xyyzz_1[i] * pc_y[i];

        ta1_x_yyyy_xyzzz_0[i] = 3.0 * ta1_x_yy_xyzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xyzzz_1[i] * fe_0 + ta1_x_yyy_xzzz_0[i] * fe_0 -
                                ta1_x_yyy_xzzz_1[i] * fe_0 + ta1_x_yyy_xyzzz_0[i] * pa_y[i] - ta1_x_yyy_xyzzz_1[i] * pc_y[i];

        ta1_x_yyyy_xzzzz_0[i] =
            3.0 * ta1_x_yy_xzzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_xzzzz_1[i] * fe_0 + ta1_x_yyy_xzzzz_0[i] * pa_y[i] - ta1_x_yyy_xzzzz_1[i] * pc_y[i];

        ta1_x_yyyy_yyyyy_0[i] = 3.0 * ta1_x_yy_yyyyy_0[i] * fe_0 - 3.0 * ta1_x_yy_yyyyy_1[i] * fe_0 + 5.0 * ta1_x_yyy_yyyy_0[i] * fe_0 -
                                5.0 * ta1_x_yyy_yyyy_1[i] * fe_0 + ta1_x_yyy_yyyyy_0[i] * pa_y[i] - ta1_x_yyy_yyyyy_1[i] * pc_y[i];

        ta1_x_yyyy_yyyyz_0[i] = 3.0 * ta1_x_yy_yyyyz_0[i] * fe_0 - 3.0 * ta1_x_yy_yyyyz_1[i] * fe_0 + 4.0 * ta1_x_yyy_yyyz_0[i] * fe_0 -
                                4.0 * ta1_x_yyy_yyyz_1[i] * fe_0 + ta1_x_yyy_yyyyz_0[i] * pa_y[i] - ta1_x_yyy_yyyyz_1[i] * pc_y[i];

        ta1_x_yyyy_yyyzz_0[i] = 3.0 * ta1_x_yy_yyyzz_0[i] * fe_0 - 3.0 * ta1_x_yy_yyyzz_1[i] * fe_0 + 3.0 * ta1_x_yyy_yyzz_0[i] * fe_0 -
                                3.0 * ta1_x_yyy_yyzz_1[i] * fe_0 + ta1_x_yyy_yyyzz_0[i] * pa_y[i] - ta1_x_yyy_yyyzz_1[i] * pc_y[i];

        ta1_x_yyyy_yyzzz_0[i] = 3.0 * ta1_x_yy_yyzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_yyzzz_1[i] * fe_0 + 2.0 * ta1_x_yyy_yzzz_0[i] * fe_0 -
                                2.0 * ta1_x_yyy_yzzz_1[i] * fe_0 + ta1_x_yyy_yyzzz_0[i] * pa_y[i] - ta1_x_yyy_yyzzz_1[i] * pc_y[i];

        ta1_x_yyyy_yzzzz_0[i] = 3.0 * ta1_x_yy_yzzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_yzzzz_1[i] * fe_0 + ta1_x_yyy_zzzz_0[i] * fe_0 -
                                ta1_x_yyy_zzzz_1[i] * fe_0 + ta1_x_yyy_yzzzz_0[i] * pa_y[i] - ta1_x_yyy_yzzzz_1[i] * pc_y[i];

        ta1_x_yyyy_zzzzz_0[i] =
            3.0 * ta1_x_yy_zzzzz_0[i] * fe_0 - 3.0 * ta1_x_yy_zzzzz_1[i] * fe_0 + ta1_x_yyy_zzzzz_0[i] * pa_y[i] - ta1_x_yyy_zzzzz_1[i] * pc_y[i];
    }

    // Set up 231-252 components of targeted buffer : GH

    auto ta1_x_yyyz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 231);

    auto ta1_x_yyyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 232);

    auto ta1_x_yyyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 233);

    auto ta1_x_yyyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 234);

    auto ta1_x_yyyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 235);

    auto ta1_x_yyyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 236);

    auto ta1_x_yyyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 237);

    auto ta1_x_yyyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 238);

    auto ta1_x_yyyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 239);

    auto ta1_x_yyyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 240);

    auto ta1_x_yyyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 241);

    auto ta1_x_yyyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 242);

    auto ta1_x_yyyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 243);

    auto ta1_x_yyyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 244);

    auto ta1_x_yyyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 245);

    auto ta1_x_yyyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 246);

    auto ta1_x_yyyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 247);

    auto ta1_x_yyyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 248);

    auto ta1_x_yyyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 249);

    auto ta1_x_yyyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 250);

    auto ta1_x_yyyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 251);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_yyy_xxxxx_0,  \
                             ta1_x_yyy_xxxxx_1,  \
                             ta1_x_yyy_xxxxy_0,  \
                             ta1_x_yyy_xxxxy_1,  \
                             ta1_x_yyy_xxxy_0,   \
                             ta1_x_yyy_xxxy_1,   \
                             ta1_x_yyy_xxxyy_0,  \
                             ta1_x_yyy_xxxyy_1,  \
                             ta1_x_yyy_xxxyz_0,  \
                             ta1_x_yyy_xxxyz_1,  \
                             ta1_x_yyy_xxyy_0,   \
                             ta1_x_yyy_xxyy_1,   \
                             ta1_x_yyy_xxyyy_0,  \
                             ta1_x_yyy_xxyyy_1,  \
                             ta1_x_yyy_xxyyz_0,  \
                             ta1_x_yyy_xxyyz_1,  \
                             ta1_x_yyy_xxyz_0,   \
                             ta1_x_yyy_xxyz_1,   \
                             ta1_x_yyy_xxyzz_0,  \
                             ta1_x_yyy_xxyzz_1,  \
                             ta1_x_yyy_xyyy_0,   \
                             ta1_x_yyy_xyyy_1,   \
                             ta1_x_yyy_xyyyy_0,  \
                             ta1_x_yyy_xyyyy_1,  \
                             ta1_x_yyy_xyyyz_0,  \
                             ta1_x_yyy_xyyyz_1,  \
                             ta1_x_yyy_xyyz_0,   \
                             ta1_x_yyy_xyyz_1,   \
                             ta1_x_yyy_xyyzz_0,  \
                             ta1_x_yyy_xyyzz_1,  \
                             ta1_x_yyy_xyzz_0,   \
                             ta1_x_yyy_xyzz_1,   \
                             ta1_x_yyy_xyzzz_0,  \
                             ta1_x_yyy_xyzzz_1,  \
                             ta1_x_yyy_yyyy_0,   \
                             ta1_x_yyy_yyyy_1,   \
                             ta1_x_yyy_yyyyy_0,  \
                             ta1_x_yyy_yyyyy_1,  \
                             ta1_x_yyy_yyyyz_0,  \
                             ta1_x_yyy_yyyyz_1,  \
                             ta1_x_yyy_yyyz_0,   \
                             ta1_x_yyy_yyyz_1,   \
                             ta1_x_yyy_yyyzz_0,  \
                             ta1_x_yyy_yyyzz_1,  \
                             ta1_x_yyy_yyzz_0,   \
                             ta1_x_yyy_yyzz_1,   \
                             ta1_x_yyy_yyzzz_0,  \
                             ta1_x_yyy_yyzzz_1,  \
                             ta1_x_yyy_yzzz_0,   \
                             ta1_x_yyy_yzzz_1,   \
                             ta1_x_yyy_yzzzz_0,  \
                             ta1_x_yyy_yzzzz_1,  \
                             ta1_x_yyyz_xxxxx_0, \
                             ta1_x_yyyz_xxxxy_0, \
                             ta1_x_yyyz_xxxxz_0, \
                             ta1_x_yyyz_xxxyy_0, \
                             ta1_x_yyyz_xxxyz_0, \
                             ta1_x_yyyz_xxxzz_0, \
                             ta1_x_yyyz_xxyyy_0, \
                             ta1_x_yyyz_xxyyz_0, \
                             ta1_x_yyyz_xxyzz_0, \
                             ta1_x_yyyz_xxzzz_0, \
                             ta1_x_yyyz_xyyyy_0, \
                             ta1_x_yyyz_xyyyz_0, \
                             ta1_x_yyyz_xyyzz_0, \
                             ta1_x_yyyz_xyzzz_0, \
                             ta1_x_yyyz_xzzzz_0, \
                             ta1_x_yyyz_yyyyy_0, \
                             ta1_x_yyyz_yyyyz_0, \
                             ta1_x_yyyz_yyyzz_0, \
                             ta1_x_yyyz_yyzzz_0, \
                             ta1_x_yyyz_yzzzz_0, \
                             ta1_x_yyyz_zzzzz_0, \
                             ta1_x_yyz_xxxxz_0,  \
                             ta1_x_yyz_xxxxz_1,  \
                             ta1_x_yyz_xxxzz_0,  \
                             ta1_x_yyz_xxxzz_1,  \
                             ta1_x_yyz_xxzzz_0,  \
                             ta1_x_yyz_xxzzz_1,  \
                             ta1_x_yyz_xzzzz_0,  \
                             ta1_x_yyz_xzzzz_1,  \
                             ta1_x_yyz_zzzzz_0,  \
                             ta1_x_yyz_zzzzz_1,  \
                             ta1_x_yz_xxxxz_0,   \
                             ta1_x_yz_xxxxz_1,   \
                             ta1_x_yz_xxxzz_0,   \
                             ta1_x_yz_xxxzz_1,   \
                             ta1_x_yz_xxzzz_0,   \
                             ta1_x_yz_xxzzz_1,   \
                             ta1_x_yz_xzzzz_0,   \
                             ta1_x_yz_xzzzz_1,   \
                             ta1_x_yz_zzzzz_0,   \
                             ta1_x_yz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyyz_xxxxx_0[i] = ta1_x_yyy_xxxxx_0[i] * pa_z[i] - ta1_x_yyy_xxxxx_1[i] * pc_z[i];

        ta1_x_yyyz_xxxxy_0[i] = ta1_x_yyy_xxxxy_0[i] * pa_z[i] - ta1_x_yyy_xxxxy_1[i] * pc_z[i];

        ta1_x_yyyz_xxxxz_0[i] =
            2.0 * ta1_x_yz_xxxxz_0[i] * fe_0 - 2.0 * ta1_x_yz_xxxxz_1[i] * fe_0 + ta1_x_yyz_xxxxz_0[i] * pa_y[i] - ta1_x_yyz_xxxxz_1[i] * pc_y[i];

        ta1_x_yyyz_xxxyy_0[i] = ta1_x_yyy_xxxyy_0[i] * pa_z[i] - ta1_x_yyy_xxxyy_1[i] * pc_z[i];

        ta1_x_yyyz_xxxyz_0[i] =
            ta1_x_yyy_xxxy_0[i] * fe_0 - ta1_x_yyy_xxxy_1[i] * fe_0 + ta1_x_yyy_xxxyz_0[i] * pa_z[i] - ta1_x_yyy_xxxyz_1[i] * pc_z[i];

        ta1_x_yyyz_xxxzz_0[i] =
            2.0 * ta1_x_yz_xxxzz_0[i] * fe_0 - 2.0 * ta1_x_yz_xxxzz_1[i] * fe_0 + ta1_x_yyz_xxxzz_0[i] * pa_y[i] - ta1_x_yyz_xxxzz_1[i] * pc_y[i];

        ta1_x_yyyz_xxyyy_0[i] = ta1_x_yyy_xxyyy_0[i] * pa_z[i] - ta1_x_yyy_xxyyy_1[i] * pc_z[i];

        ta1_x_yyyz_xxyyz_0[i] =
            ta1_x_yyy_xxyy_0[i] * fe_0 - ta1_x_yyy_xxyy_1[i] * fe_0 + ta1_x_yyy_xxyyz_0[i] * pa_z[i] - ta1_x_yyy_xxyyz_1[i] * pc_z[i];

        ta1_x_yyyz_xxyzz_0[i] =
            2.0 * ta1_x_yyy_xxyz_0[i] * fe_0 - 2.0 * ta1_x_yyy_xxyz_1[i] * fe_0 + ta1_x_yyy_xxyzz_0[i] * pa_z[i] - ta1_x_yyy_xxyzz_1[i] * pc_z[i];

        ta1_x_yyyz_xxzzz_0[i] =
            2.0 * ta1_x_yz_xxzzz_0[i] * fe_0 - 2.0 * ta1_x_yz_xxzzz_1[i] * fe_0 + ta1_x_yyz_xxzzz_0[i] * pa_y[i] - ta1_x_yyz_xxzzz_1[i] * pc_y[i];

        ta1_x_yyyz_xyyyy_0[i] = ta1_x_yyy_xyyyy_0[i] * pa_z[i] - ta1_x_yyy_xyyyy_1[i] * pc_z[i];

        ta1_x_yyyz_xyyyz_0[i] =
            ta1_x_yyy_xyyy_0[i] * fe_0 - ta1_x_yyy_xyyy_1[i] * fe_0 + ta1_x_yyy_xyyyz_0[i] * pa_z[i] - ta1_x_yyy_xyyyz_1[i] * pc_z[i];

        ta1_x_yyyz_xyyzz_0[i] =
            2.0 * ta1_x_yyy_xyyz_0[i] * fe_0 - 2.0 * ta1_x_yyy_xyyz_1[i] * fe_0 + ta1_x_yyy_xyyzz_0[i] * pa_z[i] - ta1_x_yyy_xyyzz_1[i] * pc_z[i];

        ta1_x_yyyz_xyzzz_0[i] =
            3.0 * ta1_x_yyy_xyzz_0[i] * fe_0 - 3.0 * ta1_x_yyy_xyzz_1[i] * fe_0 + ta1_x_yyy_xyzzz_0[i] * pa_z[i] - ta1_x_yyy_xyzzz_1[i] * pc_z[i];

        ta1_x_yyyz_xzzzz_0[i] =
            2.0 * ta1_x_yz_xzzzz_0[i] * fe_0 - 2.0 * ta1_x_yz_xzzzz_1[i] * fe_0 + ta1_x_yyz_xzzzz_0[i] * pa_y[i] - ta1_x_yyz_xzzzz_1[i] * pc_y[i];

        ta1_x_yyyz_yyyyy_0[i] = ta1_x_yyy_yyyyy_0[i] * pa_z[i] - ta1_x_yyy_yyyyy_1[i] * pc_z[i];

        ta1_x_yyyz_yyyyz_0[i] =
            ta1_x_yyy_yyyy_0[i] * fe_0 - ta1_x_yyy_yyyy_1[i] * fe_0 + ta1_x_yyy_yyyyz_0[i] * pa_z[i] - ta1_x_yyy_yyyyz_1[i] * pc_z[i];

        ta1_x_yyyz_yyyzz_0[i] =
            2.0 * ta1_x_yyy_yyyz_0[i] * fe_0 - 2.0 * ta1_x_yyy_yyyz_1[i] * fe_0 + ta1_x_yyy_yyyzz_0[i] * pa_z[i] - ta1_x_yyy_yyyzz_1[i] * pc_z[i];

        ta1_x_yyyz_yyzzz_0[i] =
            3.0 * ta1_x_yyy_yyzz_0[i] * fe_0 - 3.0 * ta1_x_yyy_yyzz_1[i] * fe_0 + ta1_x_yyy_yyzzz_0[i] * pa_z[i] - ta1_x_yyy_yyzzz_1[i] * pc_z[i];

        ta1_x_yyyz_yzzzz_0[i] =
            4.0 * ta1_x_yyy_yzzz_0[i] * fe_0 - 4.0 * ta1_x_yyy_yzzz_1[i] * fe_0 + ta1_x_yyy_yzzzz_0[i] * pa_z[i] - ta1_x_yyy_yzzzz_1[i] * pc_z[i];

        ta1_x_yyyz_zzzzz_0[i] =
            2.0 * ta1_x_yz_zzzzz_0[i] * fe_0 - 2.0 * ta1_x_yz_zzzzz_1[i] * fe_0 + ta1_x_yyz_zzzzz_0[i] * pa_y[i] - ta1_x_yyz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 252-273 components of targeted buffer : GH

    auto ta1_x_yyzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 252);

    auto ta1_x_yyzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 253);

    auto ta1_x_yyzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 254);

    auto ta1_x_yyzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 255);

    auto ta1_x_yyzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 256);

    auto ta1_x_yyzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 257);

    auto ta1_x_yyzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 258);

    auto ta1_x_yyzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 259);

    auto ta1_x_yyzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 260);

    auto ta1_x_yyzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 261);

    auto ta1_x_yyzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 262);

    auto ta1_x_yyzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 263);

    auto ta1_x_yyzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 264);

    auto ta1_x_yyzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 265);

    auto ta1_x_yyzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 266);

    auto ta1_x_yyzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 267);

    auto ta1_x_yyzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 268);

    auto ta1_x_yyzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 269);

    auto ta1_x_yyzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 270);

    auto ta1_x_yyzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 271);

    auto ta1_x_yyzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 272);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_x_yy_xxxxy_0,   \
                             ta1_x_yy_xxxxy_1,   \
                             ta1_x_yy_xxxyy_0,   \
                             ta1_x_yy_xxxyy_1,   \
                             ta1_x_yy_xxyyy_0,   \
                             ta1_x_yy_xxyyy_1,   \
                             ta1_x_yy_xyyyy_0,   \
                             ta1_x_yy_xyyyy_1,   \
                             ta1_x_yy_yyyyy_0,   \
                             ta1_x_yy_yyyyy_1,   \
                             ta1_x_yyz_xxxxy_0,  \
                             ta1_x_yyz_xxxxy_1,  \
                             ta1_x_yyz_xxxyy_0,  \
                             ta1_x_yyz_xxxyy_1,  \
                             ta1_x_yyz_xxyyy_0,  \
                             ta1_x_yyz_xxyyy_1,  \
                             ta1_x_yyz_xyyyy_0,  \
                             ta1_x_yyz_xyyyy_1,  \
                             ta1_x_yyz_yyyyy_0,  \
                             ta1_x_yyz_yyyyy_1,  \
                             ta1_x_yyzz_xxxxx_0, \
                             ta1_x_yyzz_xxxxy_0, \
                             ta1_x_yyzz_xxxxz_0, \
                             ta1_x_yyzz_xxxyy_0, \
                             ta1_x_yyzz_xxxyz_0, \
                             ta1_x_yyzz_xxxzz_0, \
                             ta1_x_yyzz_xxyyy_0, \
                             ta1_x_yyzz_xxyyz_0, \
                             ta1_x_yyzz_xxyzz_0, \
                             ta1_x_yyzz_xxzzz_0, \
                             ta1_x_yyzz_xyyyy_0, \
                             ta1_x_yyzz_xyyyz_0, \
                             ta1_x_yyzz_xyyzz_0, \
                             ta1_x_yyzz_xyzzz_0, \
                             ta1_x_yyzz_xzzzz_0, \
                             ta1_x_yyzz_yyyyy_0, \
                             ta1_x_yyzz_yyyyz_0, \
                             ta1_x_yyzz_yyyzz_0, \
                             ta1_x_yyzz_yyzzz_0, \
                             ta1_x_yyzz_yzzzz_0, \
                             ta1_x_yyzz_zzzzz_0, \
                             ta1_x_yzz_xxxxx_0,  \
                             ta1_x_yzz_xxxxx_1,  \
                             ta1_x_yzz_xxxxz_0,  \
                             ta1_x_yzz_xxxxz_1,  \
                             ta1_x_yzz_xxxyz_0,  \
                             ta1_x_yzz_xxxyz_1,  \
                             ta1_x_yzz_xxxz_0,   \
                             ta1_x_yzz_xxxz_1,   \
                             ta1_x_yzz_xxxzz_0,  \
                             ta1_x_yzz_xxxzz_1,  \
                             ta1_x_yzz_xxyyz_0,  \
                             ta1_x_yzz_xxyyz_1,  \
                             ta1_x_yzz_xxyz_0,   \
                             ta1_x_yzz_xxyz_1,   \
                             ta1_x_yzz_xxyzz_0,  \
                             ta1_x_yzz_xxyzz_1,  \
                             ta1_x_yzz_xxzz_0,   \
                             ta1_x_yzz_xxzz_1,   \
                             ta1_x_yzz_xxzzz_0,  \
                             ta1_x_yzz_xxzzz_1,  \
                             ta1_x_yzz_xyyyz_0,  \
                             ta1_x_yzz_xyyyz_1,  \
                             ta1_x_yzz_xyyz_0,   \
                             ta1_x_yzz_xyyz_1,   \
                             ta1_x_yzz_xyyzz_0,  \
                             ta1_x_yzz_xyyzz_1,  \
                             ta1_x_yzz_xyzz_0,   \
                             ta1_x_yzz_xyzz_1,   \
                             ta1_x_yzz_xyzzz_0,  \
                             ta1_x_yzz_xyzzz_1,  \
                             ta1_x_yzz_xzzz_0,   \
                             ta1_x_yzz_xzzz_1,   \
                             ta1_x_yzz_xzzzz_0,  \
                             ta1_x_yzz_xzzzz_1,  \
                             ta1_x_yzz_yyyyz_0,  \
                             ta1_x_yzz_yyyyz_1,  \
                             ta1_x_yzz_yyyz_0,   \
                             ta1_x_yzz_yyyz_1,   \
                             ta1_x_yzz_yyyzz_0,  \
                             ta1_x_yzz_yyyzz_1,  \
                             ta1_x_yzz_yyzz_0,   \
                             ta1_x_yzz_yyzz_1,   \
                             ta1_x_yzz_yyzzz_0,  \
                             ta1_x_yzz_yyzzz_1,  \
                             ta1_x_yzz_yzzz_0,   \
                             ta1_x_yzz_yzzz_1,   \
                             ta1_x_yzz_yzzzz_0,  \
                             ta1_x_yzz_yzzzz_1,  \
                             ta1_x_yzz_zzzz_0,   \
                             ta1_x_yzz_zzzz_1,   \
                             ta1_x_yzz_zzzzz_0,  \
                             ta1_x_yzz_zzzzz_1,  \
                             ta1_x_zz_xxxxx_0,   \
                             ta1_x_zz_xxxxx_1,   \
                             ta1_x_zz_xxxxz_0,   \
                             ta1_x_zz_xxxxz_1,   \
                             ta1_x_zz_xxxyz_0,   \
                             ta1_x_zz_xxxyz_1,   \
                             ta1_x_zz_xxxzz_0,   \
                             ta1_x_zz_xxxzz_1,   \
                             ta1_x_zz_xxyyz_0,   \
                             ta1_x_zz_xxyyz_1,   \
                             ta1_x_zz_xxyzz_0,   \
                             ta1_x_zz_xxyzz_1,   \
                             ta1_x_zz_xxzzz_0,   \
                             ta1_x_zz_xxzzz_1,   \
                             ta1_x_zz_xyyyz_0,   \
                             ta1_x_zz_xyyyz_1,   \
                             ta1_x_zz_xyyzz_0,   \
                             ta1_x_zz_xyyzz_1,   \
                             ta1_x_zz_xyzzz_0,   \
                             ta1_x_zz_xyzzz_1,   \
                             ta1_x_zz_xzzzz_0,   \
                             ta1_x_zz_xzzzz_1,   \
                             ta1_x_zz_yyyyz_0,   \
                             ta1_x_zz_yyyyz_1,   \
                             ta1_x_zz_yyyzz_0,   \
                             ta1_x_zz_yyyzz_1,   \
                             ta1_x_zz_yyzzz_0,   \
                             ta1_x_zz_yyzzz_1,   \
                             ta1_x_zz_yzzzz_0,   \
                             ta1_x_zz_yzzzz_1,   \
                             ta1_x_zz_zzzzz_0,   \
                             ta1_x_zz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yyzz_xxxxx_0[i] =
            ta1_x_zz_xxxxx_0[i] * fe_0 - ta1_x_zz_xxxxx_1[i] * fe_0 + ta1_x_yzz_xxxxx_0[i] * pa_y[i] - ta1_x_yzz_xxxxx_1[i] * pc_y[i];

        ta1_x_yyzz_xxxxy_0[i] =
            ta1_x_yy_xxxxy_0[i] * fe_0 - ta1_x_yy_xxxxy_1[i] * fe_0 + ta1_x_yyz_xxxxy_0[i] * pa_z[i] - ta1_x_yyz_xxxxy_1[i] * pc_z[i];

        ta1_x_yyzz_xxxxz_0[i] =
            ta1_x_zz_xxxxz_0[i] * fe_0 - ta1_x_zz_xxxxz_1[i] * fe_0 + ta1_x_yzz_xxxxz_0[i] * pa_y[i] - ta1_x_yzz_xxxxz_1[i] * pc_y[i];

        ta1_x_yyzz_xxxyy_0[i] =
            ta1_x_yy_xxxyy_0[i] * fe_0 - ta1_x_yy_xxxyy_1[i] * fe_0 + ta1_x_yyz_xxxyy_0[i] * pa_z[i] - ta1_x_yyz_xxxyy_1[i] * pc_z[i];

        ta1_x_yyzz_xxxyz_0[i] = ta1_x_zz_xxxyz_0[i] * fe_0 - ta1_x_zz_xxxyz_1[i] * fe_0 + ta1_x_yzz_xxxz_0[i] * fe_0 - ta1_x_yzz_xxxz_1[i] * fe_0 +
                                ta1_x_yzz_xxxyz_0[i] * pa_y[i] - ta1_x_yzz_xxxyz_1[i] * pc_y[i];

        ta1_x_yyzz_xxxzz_0[i] =
            ta1_x_zz_xxxzz_0[i] * fe_0 - ta1_x_zz_xxxzz_1[i] * fe_0 + ta1_x_yzz_xxxzz_0[i] * pa_y[i] - ta1_x_yzz_xxxzz_1[i] * pc_y[i];

        ta1_x_yyzz_xxyyy_0[i] =
            ta1_x_yy_xxyyy_0[i] * fe_0 - ta1_x_yy_xxyyy_1[i] * fe_0 + ta1_x_yyz_xxyyy_0[i] * pa_z[i] - ta1_x_yyz_xxyyy_1[i] * pc_z[i];

        ta1_x_yyzz_xxyyz_0[i] = ta1_x_zz_xxyyz_0[i] * fe_0 - ta1_x_zz_xxyyz_1[i] * fe_0 + 2.0 * ta1_x_yzz_xxyz_0[i] * fe_0 -
                                2.0 * ta1_x_yzz_xxyz_1[i] * fe_0 + ta1_x_yzz_xxyyz_0[i] * pa_y[i] - ta1_x_yzz_xxyyz_1[i] * pc_y[i];

        ta1_x_yyzz_xxyzz_0[i] = ta1_x_zz_xxyzz_0[i] * fe_0 - ta1_x_zz_xxyzz_1[i] * fe_0 + ta1_x_yzz_xxzz_0[i] * fe_0 - ta1_x_yzz_xxzz_1[i] * fe_0 +
                                ta1_x_yzz_xxyzz_0[i] * pa_y[i] - ta1_x_yzz_xxyzz_1[i] * pc_y[i];

        ta1_x_yyzz_xxzzz_0[i] =
            ta1_x_zz_xxzzz_0[i] * fe_0 - ta1_x_zz_xxzzz_1[i] * fe_0 + ta1_x_yzz_xxzzz_0[i] * pa_y[i] - ta1_x_yzz_xxzzz_1[i] * pc_y[i];

        ta1_x_yyzz_xyyyy_0[i] =
            ta1_x_yy_xyyyy_0[i] * fe_0 - ta1_x_yy_xyyyy_1[i] * fe_0 + ta1_x_yyz_xyyyy_0[i] * pa_z[i] - ta1_x_yyz_xyyyy_1[i] * pc_z[i];

        ta1_x_yyzz_xyyyz_0[i] = ta1_x_zz_xyyyz_0[i] * fe_0 - ta1_x_zz_xyyyz_1[i] * fe_0 + 3.0 * ta1_x_yzz_xyyz_0[i] * fe_0 -
                                3.0 * ta1_x_yzz_xyyz_1[i] * fe_0 + ta1_x_yzz_xyyyz_0[i] * pa_y[i] - ta1_x_yzz_xyyyz_1[i] * pc_y[i];

        ta1_x_yyzz_xyyzz_0[i] = ta1_x_zz_xyyzz_0[i] * fe_0 - ta1_x_zz_xyyzz_1[i] * fe_0 + 2.0 * ta1_x_yzz_xyzz_0[i] * fe_0 -
                                2.0 * ta1_x_yzz_xyzz_1[i] * fe_0 + ta1_x_yzz_xyyzz_0[i] * pa_y[i] - ta1_x_yzz_xyyzz_1[i] * pc_y[i];

        ta1_x_yyzz_xyzzz_0[i] = ta1_x_zz_xyzzz_0[i] * fe_0 - ta1_x_zz_xyzzz_1[i] * fe_0 + ta1_x_yzz_xzzz_0[i] * fe_0 - ta1_x_yzz_xzzz_1[i] * fe_0 +
                                ta1_x_yzz_xyzzz_0[i] * pa_y[i] - ta1_x_yzz_xyzzz_1[i] * pc_y[i];

        ta1_x_yyzz_xzzzz_0[i] =
            ta1_x_zz_xzzzz_0[i] * fe_0 - ta1_x_zz_xzzzz_1[i] * fe_0 + ta1_x_yzz_xzzzz_0[i] * pa_y[i] - ta1_x_yzz_xzzzz_1[i] * pc_y[i];

        ta1_x_yyzz_yyyyy_0[i] =
            ta1_x_yy_yyyyy_0[i] * fe_0 - ta1_x_yy_yyyyy_1[i] * fe_0 + ta1_x_yyz_yyyyy_0[i] * pa_z[i] - ta1_x_yyz_yyyyy_1[i] * pc_z[i];

        ta1_x_yyzz_yyyyz_0[i] = ta1_x_zz_yyyyz_0[i] * fe_0 - ta1_x_zz_yyyyz_1[i] * fe_0 + 4.0 * ta1_x_yzz_yyyz_0[i] * fe_0 -
                                4.0 * ta1_x_yzz_yyyz_1[i] * fe_0 + ta1_x_yzz_yyyyz_0[i] * pa_y[i] - ta1_x_yzz_yyyyz_1[i] * pc_y[i];

        ta1_x_yyzz_yyyzz_0[i] = ta1_x_zz_yyyzz_0[i] * fe_0 - ta1_x_zz_yyyzz_1[i] * fe_0 + 3.0 * ta1_x_yzz_yyzz_0[i] * fe_0 -
                                3.0 * ta1_x_yzz_yyzz_1[i] * fe_0 + ta1_x_yzz_yyyzz_0[i] * pa_y[i] - ta1_x_yzz_yyyzz_1[i] * pc_y[i];

        ta1_x_yyzz_yyzzz_0[i] = ta1_x_zz_yyzzz_0[i] * fe_0 - ta1_x_zz_yyzzz_1[i] * fe_0 + 2.0 * ta1_x_yzz_yzzz_0[i] * fe_0 -
                                2.0 * ta1_x_yzz_yzzz_1[i] * fe_0 + ta1_x_yzz_yyzzz_0[i] * pa_y[i] - ta1_x_yzz_yyzzz_1[i] * pc_y[i];

        ta1_x_yyzz_yzzzz_0[i] = ta1_x_zz_yzzzz_0[i] * fe_0 - ta1_x_zz_yzzzz_1[i] * fe_0 + ta1_x_yzz_zzzz_0[i] * fe_0 - ta1_x_yzz_zzzz_1[i] * fe_0 +
                                ta1_x_yzz_yzzzz_0[i] * pa_y[i] - ta1_x_yzz_yzzzz_1[i] * pc_y[i];

        ta1_x_yyzz_zzzzz_0[i] =
            ta1_x_zz_zzzzz_0[i] * fe_0 - ta1_x_zz_zzzzz_1[i] * fe_0 + ta1_x_yzz_zzzzz_0[i] * pa_y[i] - ta1_x_yzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 273-294 components of targeted buffer : GH

    auto ta1_x_yzzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 273);

    auto ta1_x_yzzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 274);

    auto ta1_x_yzzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 275);

    auto ta1_x_yzzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 276);

    auto ta1_x_yzzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 277);

    auto ta1_x_yzzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 278);

    auto ta1_x_yzzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 279);

    auto ta1_x_yzzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 280);

    auto ta1_x_yzzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 281);

    auto ta1_x_yzzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 282);

    auto ta1_x_yzzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 283);

    auto ta1_x_yzzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 284);

    auto ta1_x_yzzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 285);

    auto ta1_x_yzzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 286);

    auto ta1_x_yzzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 287);

    auto ta1_x_yzzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 288);

    auto ta1_x_yzzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 289);

    auto ta1_x_yzzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 290);

    auto ta1_x_yzzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 291);

    auto ta1_x_yzzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 292);

    auto ta1_x_yzzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 293);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_x_yzzz_xxxxx_0, \
                             ta1_x_yzzz_xxxxy_0, \
                             ta1_x_yzzz_xxxxz_0, \
                             ta1_x_yzzz_xxxyy_0, \
                             ta1_x_yzzz_xxxyz_0, \
                             ta1_x_yzzz_xxxzz_0, \
                             ta1_x_yzzz_xxyyy_0, \
                             ta1_x_yzzz_xxyyz_0, \
                             ta1_x_yzzz_xxyzz_0, \
                             ta1_x_yzzz_xxzzz_0, \
                             ta1_x_yzzz_xyyyy_0, \
                             ta1_x_yzzz_xyyyz_0, \
                             ta1_x_yzzz_xyyzz_0, \
                             ta1_x_yzzz_xyzzz_0, \
                             ta1_x_yzzz_xzzzz_0, \
                             ta1_x_yzzz_yyyyy_0, \
                             ta1_x_yzzz_yyyyz_0, \
                             ta1_x_yzzz_yyyzz_0, \
                             ta1_x_yzzz_yyzzz_0, \
                             ta1_x_yzzz_yzzzz_0, \
                             ta1_x_yzzz_zzzzz_0, \
                             ta1_x_zzz_xxxx_0,   \
                             ta1_x_zzz_xxxx_1,   \
                             ta1_x_zzz_xxxxx_0,  \
                             ta1_x_zzz_xxxxx_1,  \
                             ta1_x_zzz_xxxxy_0,  \
                             ta1_x_zzz_xxxxy_1,  \
                             ta1_x_zzz_xxxxz_0,  \
                             ta1_x_zzz_xxxxz_1,  \
                             ta1_x_zzz_xxxy_0,   \
                             ta1_x_zzz_xxxy_1,   \
                             ta1_x_zzz_xxxyy_0,  \
                             ta1_x_zzz_xxxyy_1,  \
                             ta1_x_zzz_xxxyz_0,  \
                             ta1_x_zzz_xxxyz_1,  \
                             ta1_x_zzz_xxxz_0,   \
                             ta1_x_zzz_xxxz_1,   \
                             ta1_x_zzz_xxxzz_0,  \
                             ta1_x_zzz_xxxzz_1,  \
                             ta1_x_zzz_xxyy_0,   \
                             ta1_x_zzz_xxyy_1,   \
                             ta1_x_zzz_xxyyy_0,  \
                             ta1_x_zzz_xxyyy_1,  \
                             ta1_x_zzz_xxyyz_0,  \
                             ta1_x_zzz_xxyyz_1,  \
                             ta1_x_zzz_xxyz_0,   \
                             ta1_x_zzz_xxyz_1,   \
                             ta1_x_zzz_xxyzz_0,  \
                             ta1_x_zzz_xxyzz_1,  \
                             ta1_x_zzz_xxzz_0,   \
                             ta1_x_zzz_xxzz_1,   \
                             ta1_x_zzz_xxzzz_0,  \
                             ta1_x_zzz_xxzzz_1,  \
                             ta1_x_zzz_xyyy_0,   \
                             ta1_x_zzz_xyyy_1,   \
                             ta1_x_zzz_xyyyy_0,  \
                             ta1_x_zzz_xyyyy_1,  \
                             ta1_x_zzz_xyyyz_0,  \
                             ta1_x_zzz_xyyyz_1,  \
                             ta1_x_zzz_xyyz_0,   \
                             ta1_x_zzz_xyyz_1,   \
                             ta1_x_zzz_xyyzz_0,  \
                             ta1_x_zzz_xyyzz_1,  \
                             ta1_x_zzz_xyzz_0,   \
                             ta1_x_zzz_xyzz_1,   \
                             ta1_x_zzz_xyzzz_0,  \
                             ta1_x_zzz_xyzzz_1,  \
                             ta1_x_zzz_xzzz_0,   \
                             ta1_x_zzz_xzzz_1,   \
                             ta1_x_zzz_xzzzz_0,  \
                             ta1_x_zzz_xzzzz_1,  \
                             ta1_x_zzz_yyyy_0,   \
                             ta1_x_zzz_yyyy_1,   \
                             ta1_x_zzz_yyyyy_0,  \
                             ta1_x_zzz_yyyyy_1,  \
                             ta1_x_zzz_yyyyz_0,  \
                             ta1_x_zzz_yyyyz_1,  \
                             ta1_x_zzz_yyyz_0,   \
                             ta1_x_zzz_yyyz_1,   \
                             ta1_x_zzz_yyyzz_0,  \
                             ta1_x_zzz_yyyzz_1,  \
                             ta1_x_zzz_yyzz_0,   \
                             ta1_x_zzz_yyzz_1,   \
                             ta1_x_zzz_yyzzz_0,  \
                             ta1_x_zzz_yyzzz_1,  \
                             ta1_x_zzz_yzzz_0,   \
                             ta1_x_zzz_yzzz_1,   \
                             ta1_x_zzz_yzzzz_0,  \
                             ta1_x_zzz_yzzzz_1,  \
                             ta1_x_zzz_zzzz_0,   \
                             ta1_x_zzz_zzzz_1,   \
                             ta1_x_zzz_zzzzz_0,  \
                             ta1_x_zzz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_yzzz_xxxxx_0[i] = ta1_x_zzz_xxxxx_0[i] * pa_y[i] - ta1_x_zzz_xxxxx_1[i] * pc_y[i];

        ta1_x_yzzz_xxxxy_0[i] =
            ta1_x_zzz_xxxx_0[i] * fe_0 - ta1_x_zzz_xxxx_1[i] * fe_0 + ta1_x_zzz_xxxxy_0[i] * pa_y[i] - ta1_x_zzz_xxxxy_1[i] * pc_y[i];

        ta1_x_yzzz_xxxxz_0[i] = ta1_x_zzz_xxxxz_0[i] * pa_y[i] - ta1_x_zzz_xxxxz_1[i] * pc_y[i];

        ta1_x_yzzz_xxxyy_0[i] =
            2.0 * ta1_x_zzz_xxxy_0[i] * fe_0 - 2.0 * ta1_x_zzz_xxxy_1[i] * fe_0 + ta1_x_zzz_xxxyy_0[i] * pa_y[i] - ta1_x_zzz_xxxyy_1[i] * pc_y[i];

        ta1_x_yzzz_xxxyz_0[i] =
            ta1_x_zzz_xxxz_0[i] * fe_0 - ta1_x_zzz_xxxz_1[i] * fe_0 + ta1_x_zzz_xxxyz_0[i] * pa_y[i] - ta1_x_zzz_xxxyz_1[i] * pc_y[i];

        ta1_x_yzzz_xxxzz_0[i] = ta1_x_zzz_xxxzz_0[i] * pa_y[i] - ta1_x_zzz_xxxzz_1[i] * pc_y[i];

        ta1_x_yzzz_xxyyy_0[i] =
            3.0 * ta1_x_zzz_xxyy_0[i] * fe_0 - 3.0 * ta1_x_zzz_xxyy_1[i] * fe_0 + ta1_x_zzz_xxyyy_0[i] * pa_y[i] - ta1_x_zzz_xxyyy_1[i] * pc_y[i];

        ta1_x_yzzz_xxyyz_0[i] =
            2.0 * ta1_x_zzz_xxyz_0[i] * fe_0 - 2.0 * ta1_x_zzz_xxyz_1[i] * fe_0 + ta1_x_zzz_xxyyz_0[i] * pa_y[i] - ta1_x_zzz_xxyyz_1[i] * pc_y[i];

        ta1_x_yzzz_xxyzz_0[i] =
            ta1_x_zzz_xxzz_0[i] * fe_0 - ta1_x_zzz_xxzz_1[i] * fe_0 + ta1_x_zzz_xxyzz_0[i] * pa_y[i] - ta1_x_zzz_xxyzz_1[i] * pc_y[i];

        ta1_x_yzzz_xxzzz_0[i] = ta1_x_zzz_xxzzz_0[i] * pa_y[i] - ta1_x_zzz_xxzzz_1[i] * pc_y[i];

        ta1_x_yzzz_xyyyy_0[i] =
            4.0 * ta1_x_zzz_xyyy_0[i] * fe_0 - 4.0 * ta1_x_zzz_xyyy_1[i] * fe_0 + ta1_x_zzz_xyyyy_0[i] * pa_y[i] - ta1_x_zzz_xyyyy_1[i] * pc_y[i];

        ta1_x_yzzz_xyyyz_0[i] =
            3.0 * ta1_x_zzz_xyyz_0[i] * fe_0 - 3.0 * ta1_x_zzz_xyyz_1[i] * fe_0 + ta1_x_zzz_xyyyz_0[i] * pa_y[i] - ta1_x_zzz_xyyyz_1[i] * pc_y[i];

        ta1_x_yzzz_xyyzz_0[i] =
            2.0 * ta1_x_zzz_xyzz_0[i] * fe_0 - 2.0 * ta1_x_zzz_xyzz_1[i] * fe_0 + ta1_x_zzz_xyyzz_0[i] * pa_y[i] - ta1_x_zzz_xyyzz_1[i] * pc_y[i];

        ta1_x_yzzz_xyzzz_0[i] =
            ta1_x_zzz_xzzz_0[i] * fe_0 - ta1_x_zzz_xzzz_1[i] * fe_0 + ta1_x_zzz_xyzzz_0[i] * pa_y[i] - ta1_x_zzz_xyzzz_1[i] * pc_y[i];

        ta1_x_yzzz_xzzzz_0[i] = ta1_x_zzz_xzzzz_0[i] * pa_y[i] - ta1_x_zzz_xzzzz_1[i] * pc_y[i];

        ta1_x_yzzz_yyyyy_0[i] =
            5.0 * ta1_x_zzz_yyyy_0[i] * fe_0 - 5.0 * ta1_x_zzz_yyyy_1[i] * fe_0 + ta1_x_zzz_yyyyy_0[i] * pa_y[i] - ta1_x_zzz_yyyyy_1[i] * pc_y[i];

        ta1_x_yzzz_yyyyz_0[i] =
            4.0 * ta1_x_zzz_yyyz_0[i] * fe_0 - 4.0 * ta1_x_zzz_yyyz_1[i] * fe_0 + ta1_x_zzz_yyyyz_0[i] * pa_y[i] - ta1_x_zzz_yyyyz_1[i] * pc_y[i];

        ta1_x_yzzz_yyyzz_0[i] =
            3.0 * ta1_x_zzz_yyzz_0[i] * fe_0 - 3.0 * ta1_x_zzz_yyzz_1[i] * fe_0 + ta1_x_zzz_yyyzz_0[i] * pa_y[i] - ta1_x_zzz_yyyzz_1[i] * pc_y[i];

        ta1_x_yzzz_yyzzz_0[i] =
            2.0 * ta1_x_zzz_yzzz_0[i] * fe_0 - 2.0 * ta1_x_zzz_yzzz_1[i] * fe_0 + ta1_x_zzz_yyzzz_0[i] * pa_y[i] - ta1_x_zzz_yyzzz_1[i] * pc_y[i];

        ta1_x_yzzz_yzzzz_0[i] =
            ta1_x_zzz_zzzz_0[i] * fe_0 - ta1_x_zzz_zzzz_1[i] * fe_0 + ta1_x_zzz_yzzzz_0[i] * pa_y[i] - ta1_x_zzz_yzzzz_1[i] * pc_y[i];

        ta1_x_yzzz_zzzzz_0[i] = ta1_x_zzz_zzzzz_0[i] * pa_y[i] - ta1_x_zzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 294-315 components of targeted buffer : GH

    auto ta1_x_zzzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 294);

    auto ta1_x_zzzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 295);

    auto ta1_x_zzzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 296);

    auto ta1_x_zzzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 297);

    auto ta1_x_zzzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 298);

    auto ta1_x_zzzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 299);

    auto ta1_x_zzzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 300);

    auto ta1_x_zzzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 301);

    auto ta1_x_zzzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 302);

    auto ta1_x_zzzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 303);

    auto ta1_x_zzzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 304);

    auto ta1_x_zzzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 305);

    auto ta1_x_zzzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 306);

    auto ta1_x_zzzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 307);

    auto ta1_x_zzzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 308);

    auto ta1_x_zzzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 309);

    auto ta1_x_zzzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 310);

    auto ta1_x_zzzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 311);

    auto ta1_x_zzzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 312);

    auto ta1_x_zzzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 313);

    auto ta1_x_zzzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 314);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_x_zz_xxxxx_0,   \
                             ta1_x_zz_xxxxx_1,   \
                             ta1_x_zz_xxxxy_0,   \
                             ta1_x_zz_xxxxy_1,   \
                             ta1_x_zz_xxxxz_0,   \
                             ta1_x_zz_xxxxz_1,   \
                             ta1_x_zz_xxxyy_0,   \
                             ta1_x_zz_xxxyy_1,   \
                             ta1_x_zz_xxxyz_0,   \
                             ta1_x_zz_xxxyz_1,   \
                             ta1_x_zz_xxxzz_0,   \
                             ta1_x_zz_xxxzz_1,   \
                             ta1_x_zz_xxyyy_0,   \
                             ta1_x_zz_xxyyy_1,   \
                             ta1_x_zz_xxyyz_0,   \
                             ta1_x_zz_xxyyz_1,   \
                             ta1_x_zz_xxyzz_0,   \
                             ta1_x_zz_xxyzz_1,   \
                             ta1_x_zz_xxzzz_0,   \
                             ta1_x_zz_xxzzz_1,   \
                             ta1_x_zz_xyyyy_0,   \
                             ta1_x_zz_xyyyy_1,   \
                             ta1_x_zz_xyyyz_0,   \
                             ta1_x_zz_xyyyz_1,   \
                             ta1_x_zz_xyyzz_0,   \
                             ta1_x_zz_xyyzz_1,   \
                             ta1_x_zz_xyzzz_0,   \
                             ta1_x_zz_xyzzz_1,   \
                             ta1_x_zz_xzzzz_0,   \
                             ta1_x_zz_xzzzz_1,   \
                             ta1_x_zz_yyyyy_0,   \
                             ta1_x_zz_yyyyy_1,   \
                             ta1_x_zz_yyyyz_0,   \
                             ta1_x_zz_yyyyz_1,   \
                             ta1_x_zz_yyyzz_0,   \
                             ta1_x_zz_yyyzz_1,   \
                             ta1_x_zz_yyzzz_0,   \
                             ta1_x_zz_yyzzz_1,   \
                             ta1_x_zz_yzzzz_0,   \
                             ta1_x_zz_yzzzz_1,   \
                             ta1_x_zz_zzzzz_0,   \
                             ta1_x_zz_zzzzz_1,   \
                             ta1_x_zzz_xxxx_0,   \
                             ta1_x_zzz_xxxx_1,   \
                             ta1_x_zzz_xxxxx_0,  \
                             ta1_x_zzz_xxxxx_1,  \
                             ta1_x_zzz_xxxxy_0,  \
                             ta1_x_zzz_xxxxy_1,  \
                             ta1_x_zzz_xxxxz_0,  \
                             ta1_x_zzz_xxxxz_1,  \
                             ta1_x_zzz_xxxy_0,   \
                             ta1_x_zzz_xxxy_1,   \
                             ta1_x_zzz_xxxyy_0,  \
                             ta1_x_zzz_xxxyy_1,  \
                             ta1_x_zzz_xxxyz_0,  \
                             ta1_x_zzz_xxxyz_1,  \
                             ta1_x_zzz_xxxz_0,   \
                             ta1_x_zzz_xxxz_1,   \
                             ta1_x_zzz_xxxzz_0,  \
                             ta1_x_zzz_xxxzz_1,  \
                             ta1_x_zzz_xxyy_0,   \
                             ta1_x_zzz_xxyy_1,   \
                             ta1_x_zzz_xxyyy_0,  \
                             ta1_x_zzz_xxyyy_1,  \
                             ta1_x_zzz_xxyyz_0,  \
                             ta1_x_zzz_xxyyz_1,  \
                             ta1_x_zzz_xxyz_0,   \
                             ta1_x_zzz_xxyz_1,   \
                             ta1_x_zzz_xxyzz_0,  \
                             ta1_x_zzz_xxyzz_1,  \
                             ta1_x_zzz_xxzz_0,   \
                             ta1_x_zzz_xxzz_1,   \
                             ta1_x_zzz_xxzzz_0,  \
                             ta1_x_zzz_xxzzz_1,  \
                             ta1_x_zzz_xyyy_0,   \
                             ta1_x_zzz_xyyy_1,   \
                             ta1_x_zzz_xyyyy_0,  \
                             ta1_x_zzz_xyyyy_1,  \
                             ta1_x_zzz_xyyyz_0,  \
                             ta1_x_zzz_xyyyz_1,  \
                             ta1_x_zzz_xyyz_0,   \
                             ta1_x_zzz_xyyz_1,   \
                             ta1_x_zzz_xyyzz_0,  \
                             ta1_x_zzz_xyyzz_1,  \
                             ta1_x_zzz_xyzz_0,   \
                             ta1_x_zzz_xyzz_1,   \
                             ta1_x_zzz_xyzzz_0,  \
                             ta1_x_zzz_xyzzz_1,  \
                             ta1_x_zzz_xzzz_0,   \
                             ta1_x_zzz_xzzz_1,   \
                             ta1_x_zzz_xzzzz_0,  \
                             ta1_x_zzz_xzzzz_1,  \
                             ta1_x_zzz_yyyy_0,   \
                             ta1_x_zzz_yyyy_1,   \
                             ta1_x_zzz_yyyyy_0,  \
                             ta1_x_zzz_yyyyy_1,  \
                             ta1_x_zzz_yyyyz_0,  \
                             ta1_x_zzz_yyyyz_1,  \
                             ta1_x_zzz_yyyz_0,   \
                             ta1_x_zzz_yyyz_1,   \
                             ta1_x_zzz_yyyzz_0,  \
                             ta1_x_zzz_yyyzz_1,  \
                             ta1_x_zzz_yyzz_0,   \
                             ta1_x_zzz_yyzz_1,   \
                             ta1_x_zzz_yyzzz_0,  \
                             ta1_x_zzz_yyzzz_1,  \
                             ta1_x_zzz_yzzz_0,   \
                             ta1_x_zzz_yzzz_1,   \
                             ta1_x_zzz_yzzzz_0,  \
                             ta1_x_zzz_yzzzz_1,  \
                             ta1_x_zzz_zzzz_0,   \
                             ta1_x_zzz_zzzz_1,   \
                             ta1_x_zzz_zzzzz_0,  \
                             ta1_x_zzz_zzzzz_1,  \
                             ta1_x_zzzz_xxxxx_0, \
                             ta1_x_zzzz_xxxxy_0, \
                             ta1_x_zzzz_xxxxz_0, \
                             ta1_x_zzzz_xxxyy_0, \
                             ta1_x_zzzz_xxxyz_0, \
                             ta1_x_zzzz_xxxzz_0, \
                             ta1_x_zzzz_xxyyy_0, \
                             ta1_x_zzzz_xxyyz_0, \
                             ta1_x_zzzz_xxyzz_0, \
                             ta1_x_zzzz_xxzzz_0, \
                             ta1_x_zzzz_xyyyy_0, \
                             ta1_x_zzzz_xyyyz_0, \
                             ta1_x_zzzz_xyyzz_0, \
                             ta1_x_zzzz_xyzzz_0, \
                             ta1_x_zzzz_xzzzz_0, \
                             ta1_x_zzzz_yyyyy_0, \
                             ta1_x_zzzz_yyyyz_0, \
                             ta1_x_zzzz_yyyzz_0, \
                             ta1_x_zzzz_yyzzz_0, \
                             ta1_x_zzzz_yzzzz_0, \
                             ta1_x_zzzz_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_x_zzzz_xxxxx_0[i] =
            3.0 * ta1_x_zz_xxxxx_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxxx_1[i] * fe_0 + ta1_x_zzz_xxxxx_0[i] * pa_z[i] - ta1_x_zzz_xxxxx_1[i] * pc_z[i];

        ta1_x_zzzz_xxxxy_0[i] =
            3.0 * ta1_x_zz_xxxxy_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxxy_1[i] * fe_0 + ta1_x_zzz_xxxxy_0[i] * pa_z[i] - ta1_x_zzz_xxxxy_1[i] * pc_z[i];

        ta1_x_zzzz_xxxxz_0[i] = 3.0 * ta1_x_zz_xxxxz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxxz_1[i] * fe_0 + ta1_x_zzz_xxxx_0[i] * fe_0 -
                                ta1_x_zzz_xxxx_1[i] * fe_0 + ta1_x_zzz_xxxxz_0[i] * pa_z[i] - ta1_x_zzz_xxxxz_1[i] * pc_z[i];

        ta1_x_zzzz_xxxyy_0[i] =
            3.0 * ta1_x_zz_xxxyy_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxyy_1[i] * fe_0 + ta1_x_zzz_xxxyy_0[i] * pa_z[i] - ta1_x_zzz_xxxyy_1[i] * pc_z[i];

        ta1_x_zzzz_xxxyz_0[i] = 3.0 * ta1_x_zz_xxxyz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxyz_1[i] * fe_0 + ta1_x_zzz_xxxy_0[i] * fe_0 -
                                ta1_x_zzz_xxxy_1[i] * fe_0 + ta1_x_zzz_xxxyz_0[i] * pa_z[i] - ta1_x_zzz_xxxyz_1[i] * pc_z[i];

        ta1_x_zzzz_xxxzz_0[i] = 3.0 * ta1_x_zz_xxxzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxxzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xxxz_0[i] * fe_0 -
                                2.0 * ta1_x_zzz_xxxz_1[i] * fe_0 + ta1_x_zzz_xxxzz_0[i] * pa_z[i] - ta1_x_zzz_xxxzz_1[i] * pc_z[i];

        ta1_x_zzzz_xxyyy_0[i] =
            3.0 * ta1_x_zz_xxyyy_0[i] * fe_0 - 3.0 * ta1_x_zz_xxyyy_1[i] * fe_0 + ta1_x_zzz_xxyyy_0[i] * pa_z[i] - ta1_x_zzz_xxyyy_1[i] * pc_z[i];

        ta1_x_zzzz_xxyyz_0[i] = 3.0 * ta1_x_zz_xxyyz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxyyz_1[i] * fe_0 + ta1_x_zzz_xxyy_0[i] * fe_0 -
                                ta1_x_zzz_xxyy_1[i] * fe_0 + ta1_x_zzz_xxyyz_0[i] * pa_z[i] - ta1_x_zzz_xxyyz_1[i] * pc_z[i];

        ta1_x_zzzz_xxyzz_0[i] = 3.0 * ta1_x_zz_xxyzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxyzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xxyz_0[i] * fe_0 -
                                2.0 * ta1_x_zzz_xxyz_1[i] * fe_0 + ta1_x_zzz_xxyzz_0[i] * pa_z[i] - ta1_x_zzz_xxyzz_1[i] * pc_z[i];

        ta1_x_zzzz_xxzzz_0[i] = 3.0 * ta1_x_zz_xxzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xxzzz_1[i] * fe_0 + 3.0 * ta1_x_zzz_xxzz_0[i] * fe_0 -
                                3.0 * ta1_x_zzz_xxzz_1[i] * fe_0 + ta1_x_zzz_xxzzz_0[i] * pa_z[i] - ta1_x_zzz_xxzzz_1[i] * pc_z[i];

        ta1_x_zzzz_xyyyy_0[i] =
            3.0 * ta1_x_zz_xyyyy_0[i] * fe_0 - 3.0 * ta1_x_zz_xyyyy_1[i] * fe_0 + ta1_x_zzz_xyyyy_0[i] * pa_z[i] - ta1_x_zzz_xyyyy_1[i] * pc_z[i];

        ta1_x_zzzz_xyyyz_0[i] = 3.0 * ta1_x_zz_xyyyz_0[i] * fe_0 - 3.0 * ta1_x_zz_xyyyz_1[i] * fe_0 + ta1_x_zzz_xyyy_0[i] * fe_0 -
                                ta1_x_zzz_xyyy_1[i] * fe_0 + ta1_x_zzz_xyyyz_0[i] * pa_z[i] - ta1_x_zzz_xyyyz_1[i] * pc_z[i];

        ta1_x_zzzz_xyyzz_0[i] = 3.0 * ta1_x_zz_xyyzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xyyzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_xyyz_0[i] * fe_0 -
                                2.0 * ta1_x_zzz_xyyz_1[i] * fe_0 + ta1_x_zzz_xyyzz_0[i] * pa_z[i] - ta1_x_zzz_xyyzz_1[i] * pc_z[i];

        ta1_x_zzzz_xyzzz_0[i] = 3.0 * ta1_x_zz_xyzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xyzzz_1[i] * fe_0 + 3.0 * ta1_x_zzz_xyzz_0[i] * fe_0 -
                                3.0 * ta1_x_zzz_xyzz_1[i] * fe_0 + ta1_x_zzz_xyzzz_0[i] * pa_z[i] - ta1_x_zzz_xyzzz_1[i] * pc_z[i];

        ta1_x_zzzz_xzzzz_0[i] = 3.0 * ta1_x_zz_xzzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_xzzzz_1[i] * fe_0 + 4.0 * ta1_x_zzz_xzzz_0[i] * fe_0 -
                                4.0 * ta1_x_zzz_xzzz_1[i] * fe_0 + ta1_x_zzz_xzzzz_0[i] * pa_z[i] - ta1_x_zzz_xzzzz_1[i] * pc_z[i];

        ta1_x_zzzz_yyyyy_0[i] =
            3.0 * ta1_x_zz_yyyyy_0[i] * fe_0 - 3.0 * ta1_x_zz_yyyyy_1[i] * fe_0 + ta1_x_zzz_yyyyy_0[i] * pa_z[i] - ta1_x_zzz_yyyyy_1[i] * pc_z[i];

        ta1_x_zzzz_yyyyz_0[i] = 3.0 * ta1_x_zz_yyyyz_0[i] * fe_0 - 3.0 * ta1_x_zz_yyyyz_1[i] * fe_0 + ta1_x_zzz_yyyy_0[i] * fe_0 -
                                ta1_x_zzz_yyyy_1[i] * fe_0 + ta1_x_zzz_yyyyz_0[i] * pa_z[i] - ta1_x_zzz_yyyyz_1[i] * pc_z[i];

        ta1_x_zzzz_yyyzz_0[i] = 3.0 * ta1_x_zz_yyyzz_0[i] * fe_0 - 3.0 * ta1_x_zz_yyyzz_1[i] * fe_0 + 2.0 * ta1_x_zzz_yyyz_0[i] * fe_0 -
                                2.0 * ta1_x_zzz_yyyz_1[i] * fe_0 + ta1_x_zzz_yyyzz_0[i] * pa_z[i] - ta1_x_zzz_yyyzz_1[i] * pc_z[i];

        ta1_x_zzzz_yyzzz_0[i] = 3.0 * ta1_x_zz_yyzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_yyzzz_1[i] * fe_0 + 3.0 * ta1_x_zzz_yyzz_0[i] * fe_0 -
                                3.0 * ta1_x_zzz_yyzz_1[i] * fe_0 + ta1_x_zzz_yyzzz_0[i] * pa_z[i] - ta1_x_zzz_yyzzz_1[i] * pc_z[i];

        ta1_x_zzzz_yzzzz_0[i] = 3.0 * ta1_x_zz_yzzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_yzzzz_1[i] * fe_0 + 4.0 * ta1_x_zzz_yzzz_0[i] * fe_0 -
                                4.0 * ta1_x_zzz_yzzz_1[i] * fe_0 + ta1_x_zzz_yzzzz_0[i] * pa_z[i] - ta1_x_zzz_yzzzz_1[i] * pc_z[i];

        ta1_x_zzzz_zzzzz_0[i] = 3.0 * ta1_x_zz_zzzzz_0[i] * fe_0 - 3.0 * ta1_x_zz_zzzzz_1[i] * fe_0 + 5.0 * ta1_x_zzz_zzzz_0[i] * fe_0 -
                                5.0 * ta1_x_zzz_zzzz_1[i] * fe_0 + ta1_x_zzz_zzzzz_0[i] * pa_z[i] - ta1_x_zzz_zzzzz_1[i] * pc_z[i];
    }

    // Set up 315-336 components of targeted buffer : GH

    auto ta1_y_xxxx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 315);

    auto ta1_y_xxxx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 316);

    auto ta1_y_xxxx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 317);

    auto ta1_y_xxxx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 318);

    auto ta1_y_xxxx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 319);

    auto ta1_y_xxxx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 320);

    auto ta1_y_xxxx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 321);

    auto ta1_y_xxxx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 322);

    auto ta1_y_xxxx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 323);

    auto ta1_y_xxxx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 324);

    auto ta1_y_xxxx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 325);

    auto ta1_y_xxxx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 326);

    auto ta1_y_xxxx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 327);

    auto ta1_y_xxxx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 328);

    auto ta1_y_xxxx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 329);

    auto ta1_y_xxxx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 330);

    auto ta1_y_xxxx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 331);

    auto ta1_y_xxxx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 332);

    auto ta1_y_xxxx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 333);

    auto ta1_y_xxxx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 334);

    auto ta1_y_xxxx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 335);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_xx_xxxxx_0,   \
                             ta1_y_xx_xxxxx_1,   \
                             ta1_y_xx_xxxxy_0,   \
                             ta1_y_xx_xxxxy_1,   \
                             ta1_y_xx_xxxxz_0,   \
                             ta1_y_xx_xxxxz_1,   \
                             ta1_y_xx_xxxyy_0,   \
                             ta1_y_xx_xxxyy_1,   \
                             ta1_y_xx_xxxyz_0,   \
                             ta1_y_xx_xxxyz_1,   \
                             ta1_y_xx_xxxzz_0,   \
                             ta1_y_xx_xxxzz_1,   \
                             ta1_y_xx_xxyyy_0,   \
                             ta1_y_xx_xxyyy_1,   \
                             ta1_y_xx_xxyyz_0,   \
                             ta1_y_xx_xxyyz_1,   \
                             ta1_y_xx_xxyzz_0,   \
                             ta1_y_xx_xxyzz_1,   \
                             ta1_y_xx_xxzzz_0,   \
                             ta1_y_xx_xxzzz_1,   \
                             ta1_y_xx_xyyyy_0,   \
                             ta1_y_xx_xyyyy_1,   \
                             ta1_y_xx_xyyyz_0,   \
                             ta1_y_xx_xyyyz_1,   \
                             ta1_y_xx_xyyzz_0,   \
                             ta1_y_xx_xyyzz_1,   \
                             ta1_y_xx_xyzzz_0,   \
                             ta1_y_xx_xyzzz_1,   \
                             ta1_y_xx_xzzzz_0,   \
                             ta1_y_xx_xzzzz_1,   \
                             ta1_y_xx_yyyyy_0,   \
                             ta1_y_xx_yyyyy_1,   \
                             ta1_y_xx_yyyyz_0,   \
                             ta1_y_xx_yyyyz_1,   \
                             ta1_y_xx_yyyzz_0,   \
                             ta1_y_xx_yyyzz_1,   \
                             ta1_y_xx_yyzzz_0,   \
                             ta1_y_xx_yyzzz_1,   \
                             ta1_y_xx_yzzzz_0,   \
                             ta1_y_xx_yzzzz_1,   \
                             ta1_y_xx_zzzzz_0,   \
                             ta1_y_xx_zzzzz_1,   \
                             ta1_y_xxx_xxxx_0,   \
                             ta1_y_xxx_xxxx_1,   \
                             ta1_y_xxx_xxxxx_0,  \
                             ta1_y_xxx_xxxxx_1,  \
                             ta1_y_xxx_xxxxy_0,  \
                             ta1_y_xxx_xxxxy_1,  \
                             ta1_y_xxx_xxxxz_0,  \
                             ta1_y_xxx_xxxxz_1,  \
                             ta1_y_xxx_xxxy_0,   \
                             ta1_y_xxx_xxxy_1,   \
                             ta1_y_xxx_xxxyy_0,  \
                             ta1_y_xxx_xxxyy_1,  \
                             ta1_y_xxx_xxxyz_0,  \
                             ta1_y_xxx_xxxyz_1,  \
                             ta1_y_xxx_xxxz_0,   \
                             ta1_y_xxx_xxxz_1,   \
                             ta1_y_xxx_xxxzz_0,  \
                             ta1_y_xxx_xxxzz_1,  \
                             ta1_y_xxx_xxyy_0,   \
                             ta1_y_xxx_xxyy_1,   \
                             ta1_y_xxx_xxyyy_0,  \
                             ta1_y_xxx_xxyyy_1,  \
                             ta1_y_xxx_xxyyz_0,  \
                             ta1_y_xxx_xxyyz_1,  \
                             ta1_y_xxx_xxyz_0,   \
                             ta1_y_xxx_xxyz_1,   \
                             ta1_y_xxx_xxyzz_0,  \
                             ta1_y_xxx_xxyzz_1,  \
                             ta1_y_xxx_xxzz_0,   \
                             ta1_y_xxx_xxzz_1,   \
                             ta1_y_xxx_xxzzz_0,  \
                             ta1_y_xxx_xxzzz_1,  \
                             ta1_y_xxx_xyyy_0,   \
                             ta1_y_xxx_xyyy_1,   \
                             ta1_y_xxx_xyyyy_0,  \
                             ta1_y_xxx_xyyyy_1,  \
                             ta1_y_xxx_xyyyz_0,  \
                             ta1_y_xxx_xyyyz_1,  \
                             ta1_y_xxx_xyyz_0,   \
                             ta1_y_xxx_xyyz_1,   \
                             ta1_y_xxx_xyyzz_0,  \
                             ta1_y_xxx_xyyzz_1,  \
                             ta1_y_xxx_xyzz_0,   \
                             ta1_y_xxx_xyzz_1,   \
                             ta1_y_xxx_xyzzz_0,  \
                             ta1_y_xxx_xyzzz_1,  \
                             ta1_y_xxx_xzzz_0,   \
                             ta1_y_xxx_xzzz_1,   \
                             ta1_y_xxx_xzzzz_0,  \
                             ta1_y_xxx_xzzzz_1,  \
                             ta1_y_xxx_yyyy_0,   \
                             ta1_y_xxx_yyyy_1,   \
                             ta1_y_xxx_yyyyy_0,  \
                             ta1_y_xxx_yyyyy_1,  \
                             ta1_y_xxx_yyyyz_0,  \
                             ta1_y_xxx_yyyyz_1,  \
                             ta1_y_xxx_yyyz_0,   \
                             ta1_y_xxx_yyyz_1,   \
                             ta1_y_xxx_yyyzz_0,  \
                             ta1_y_xxx_yyyzz_1,  \
                             ta1_y_xxx_yyzz_0,   \
                             ta1_y_xxx_yyzz_1,   \
                             ta1_y_xxx_yyzzz_0,  \
                             ta1_y_xxx_yyzzz_1,  \
                             ta1_y_xxx_yzzz_0,   \
                             ta1_y_xxx_yzzz_1,   \
                             ta1_y_xxx_yzzzz_0,  \
                             ta1_y_xxx_yzzzz_1,  \
                             ta1_y_xxx_zzzz_0,   \
                             ta1_y_xxx_zzzz_1,   \
                             ta1_y_xxx_zzzzz_0,  \
                             ta1_y_xxx_zzzzz_1,  \
                             ta1_y_xxxx_xxxxx_0, \
                             ta1_y_xxxx_xxxxy_0, \
                             ta1_y_xxxx_xxxxz_0, \
                             ta1_y_xxxx_xxxyy_0, \
                             ta1_y_xxxx_xxxyz_0, \
                             ta1_y_xxxx_xxxzz_0, \
                             ta1_y_xxxx_xxyyy_0, \
                             ta1_y_xxxx_xxyyz_0, \
                             ta1_y_xxxx_xxyzz_0, \
                             ta1_y_xxxx_xxzzz_0, \
                             ta1_y_xxxx_xyyyy_0, \
                             ta1_y_xxxx_xyyyz_0, \
                             ta1_y_xxxx_xyyzz_0, \
                             ta1_y_xxxx_xyzzz_0, \
                             ta1_y_xxxx_xzzzz_0, \
                             ta1_y_xxxx_yyyyy_0, \
                             ta1_y_xxxx_yyyyz_0, \
                             ta1_y_xxxx_yyyzz_0, \
                             ta1_y_xxxx_yyzzz_0, \
                             ta1_y_xxxx_yzzzz_0, \
                             ta1_y_xxxx_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxx_xxxxx_0[i] = 3.0 * ta1_y_xx_xxxxx_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxxx_1[i] * fe_0 + 5.0 * ta1_y_xxx_xxxx_0[i] * fe_0 -
                                5.0 * ta1_y_xxx_xxxx_1[i] * fe_0 + ta1_y_xxx_xxxxx_0[i] * pa_x[i] - ta1_y_xxx_xxxxx_1[i] * pc_x[i];

        ta1_y_xxxx_xxxxy_0[i] = 3.0 * ta1_y_xx_xxxxy_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxxy_1[i] * fe_0 + 4.0 * ta1_y_xxx_xxxy_0[i] * fe_0 -
                                4.0 * ta1_y_xxx_xxxy_1[i] * fe_0 + ta1_y_xxx_xxxxy_0[i] * pa_x[i] - ta1_y_xxx_xxxxy_1[i] * pc_x[i];

        ta1_y_xxxx_xxxxz_0[i] = 3.0 * ta1_y_xx_xxxxz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxxz_1[i] * fe_0 + 4.0 * ta1_y_xxx_xxxz_0[i] * fe_0 -
                                4.0 * ta1_y_xxx_xxxz_1[i] * fe_0 + ta1_y_xxx_xxxxz_0[i] * pa_x[i] - ta1_y_xxx_xxxxz_1[i] * pc_x[i];

        ta1_y_xxxx_xxxyy_0[i] = 3.0 * ta1_y_xx_xxxyy_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxyy_1[i] * fe_0 + 3.0 * ta1_y_xxx_xxyy_0[i] * fe_0 -
                                3.0 * ta1_y_xxx_xxyy_1[i] * fe_0 + ta1_y_xxx_xxxyy_0[i] * pa_x[i] - ta1_y_xxx_xxxyy_1[i] * pc_x[i];

        ta1_y_xxxx_xxxyz_0[i] = 3.0 * ta1_y_xx_xxxyz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxyz_1[i] * fe_0 + 3.0 * ta1_y_xxx_xxyz_0[i] * fe_0 -
                                3.0 * ta1_y_xxx_xxyz_1[i] * fe_0 + ta1_y_xxx_xxxyz_0[i] * pa_x[i] - ta1_y_xxx_xxxyz_1[i] * pc_x[i];

        ta1_y_xxxx_xxxzz_0[i] = 3.0 * ta1_y_xx_xxxzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxxzz_1[i] * fe_0 + 3.0 * ta1_y_xxx_xxzz_0[i] * fe_0 -
                                3.0 * ta1_y_xxx_xxzz_1[i] * fe_0 + ta1_y_xxx_xxxzz_0[i] * pa_x[i] - ta1_y_xxx_xxxzz_1[i] * pc_x[i];

        ta1_y_xxxx_xxyyy_0[i] = 3.0 * ta1_y_xx_xxyyy_0[i] * fe_0 - 3.0 * ta1_y_xx_xxyyy_1[i] * fe_0 + 2.0 * ta1_y_xxx_xyyy_0[i] * fe_0 -
                                2.0 * ta1_y_xxx_xyyy_1[i] * fe_0 + ta1_y_xxx_xxyyy_0[i] * pa_x[i] - ta1_y_xxx_xxyyy_1[i] * pc_x[i];

        ta1_y_xxxx_xxyyz_0[i] = 3.0 * ta1_y_xx_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxyyz_1[i] * fe_0 + 2.0 * ta1_y_xxx_xyyz_0[i] * fe_0 -
                                2.0 * ta1_y_xxx_xyyz_1[i] * fe_0 + ta1_y_xxx_xxyyz_0[i] * pa_x[i] - ta1_y_xxx_xxyyz_1[i] * pc_x[i];

        ta1_y_xxxx_xxyzz_0[i] = 3.0 * ta1_y_xx_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxyzz_1[i] * fe_0 + 2.0 * ta1_y_xxx_xyzz_0[i] * fe_0 -
                                2.0 * ta1_y_xxx_xyzz_1[i] * fe_0 + ta1_y_xxx_xxyzz_0[i] * pa_x[i] - ta1_y_xxx_xxyzz_1[i] * pc_x[i];

        ta1_y_xxxx_xxzzz_0[i] = 3.0 * ta1_y_xx_xxzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xxzzz_1[i] * fe_0 + 2.0 * ta1_y_xxx_xzzz_0[i] * fe_0 -
                                2.0 * ta1_y_xxx_xzzz_1[i] * fe_0 + ta1_y_xxx_xxzzz_0[i] * pa_x[i] - ta1_y_xxx_xxzzz_1[i] * pc_x[i];

        ta1_y_xxxx_xyyyy_0[i] = 3.0 * ta1_y_xx_xyyyy_0[i] * fe_0 - 3.0 * ta1_y_xx_xyyyy_1[i] * fe_0 + ta1_y_xxx_yyyy_0[i] * fe_0 -
                                ta1_y_xxx_yyyy_1[i] * fe_0 + ta1_y_xxx_xyyyy_0[i] * pa_x[i] - ta1_y_xxx_xyyyy_1[i] * pc_x[i];

        ta1_y_xxxx_xyyyz_0[i] = 3.0 * ta1_y_xx_xyyyz_0[i] * fe_0 - 3.0 * ta1_y_xx_xyyyz_1[i] * fe_0 + ta1_y_xxx_yyyz_0[i] * fe_0 -
                                ta1_y_xxx_yyyz_1[i] * fe_0 + ta1_y_xxx_xyyyz_0[i] * pa_x[i] - ta1_y_xxx_xyyyz_1[i] * pc_x[i];

        ta1_y_xxxx_xyyzz_0[i] = 3.0 * ta1_y_xx_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xyyzz_1[i] * fe_0 + ta1_y_xxx_yyzz_0[i] * fe_0 -
                                ta1_y_xxx_yyzz_1[i] * fe_0 + ta1_y_xxx_xyyzz_0[i] * pa_x[i] - ta1_y_xxx_xyyzz_1[i] * pc_x[i];

        ta1_y_xxxx_xyzzz_0[i] = 3.0 * ta1_y_xx_xyzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xyzzz_1[i] * fe_0 + ta1_y_xxx_yzzz_0[i] * fe_0 -
                                ta1_y_xxx_yzzz_1[i] * fe_0 + ta1_y_xxx_xyzzz_0[i] * pa_x[i] - ta1_y_xxx_xyzzz_1[i] * pc_x[i];

        ta1_y_xxxx_xzzzz_0[i] = 3.0 * ta1_y_xx_xzzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_xzzzz_1[i] * fe_0 + ta1_y_xxx_zzzz_0[i] * fe_0 -
                                ta1_y_xxx_zzzz_1[i] * fe_0 + ta1_y_xxx_xzzzz_0[i] * pa_x[i] - ta1_y_xxx_xzzzz_1[i] * pc_x[i];

        ta1_y_xxxx_yyyyy_0[i] =
            3.0 * ta1_y_xx_yyyyy_0[i] * fe_0 - 3.0 * ta1_y_xx_yyyyy_1[i] * fe_0 + ta1_y_xxx_yyyyy_0[i] * pa_x[i] - ta1_y_xxx_yyyyy_1[i] * pc_x[i];

        ta1_y_xxxx_yyyyz_0[i] =
            3.0 * ta1_y_xx_yyyyz_0[i] * fe_0 - 3.0 * ta1_y_xx_yyyyz_1[i] * fe_0 + ta1_y_xxx_yyyyz_0[i] * pa_x[i] - ta1_y_xxx_yyyyz_1[i] * pc_x[i];

        ta1_y_xxxx_yyyzz_0[i] =
            3.0 * ta1_y_xx_yyyzz_0[i] * fe_0 - 3.0 * ta1_y_xx_yyyzz_1[i] * fe_0 + ta1_y_xxx_yyyzz_0[i] * pa_x[i] - ta1_y_xxx_yyyzz_1[i] * pc_x[i];

        ta1_y_xxxx_yyzzz_0[i] =
            3.0 * ta1_y_xx_yyzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_yyzzz_1[i] * fe_0 + ta1_y_xxx_yyzzz_0[i] * pa_x[i] - ta1_y_xxx_yyzzz_1[i] * pc_x[i];

        ta1_y_xxxx_yzzzz_0[i] =
            3.0 * ta1_y_xx_yzzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_yzzzz_1[i] * fe_0 + ta1_y_xxx_yzzzz_0[i] * pa_x[i] - ta1_y_xxx_yzzzz_1[i] * pc_x[i];

        ta1_y_xxxx_zzzzz_0[i] =
            3.0 * ta1_y_xx_zzzzz_0[i] * fe_0 - 3.0 * ta1_y_xx_zzzzz_1[i] * fe_0 + ta1_y_xxx_zzzzz_0[i] * pa_x[i] - ta1_y_xxx_zzzzz_1[i] * pc_x[i];
    }

    // Set up 336-357 components of targeted buffer : GH

    auto ta1_y_xxxy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 336);

    auto ta1_y_xxxy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 337);

    auto ta1_y_xxxy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 338);

    auto ta1_y_xxxy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 339);

    auto ta1_y_xxxy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 340);

    auto ta1_y_xxxy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 341);

    auto ta1_y_xxxy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 342);

    auto ta1_y_xxxy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 343);

    auto ta1_y_xxxy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 344);

    auto ta1_y_xxxy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 345);

    auto ta1_y_xxxy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 346);

    auto ta1_y_xxxy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 347);

    auto ta1_y_xxxy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 348);

    auto ta1_y_xxxy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 349);

    auto ta1_y_xxxy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 350);

    auto ta1_y_xxxy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 351);

    auto ta1_y_xxxy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 352);

    auto ta1_y_xxxy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 353);

    auto ta1_y_xxxy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 354);

    auto ta1_y_xxxy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 355);

    auto ta1_y_xxxy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 356);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_y_xxx_xxxx_0,   \
                             ta1_y_xxx_xxxx_1,   \
                             ta1_y_xxx_xxxxx_0,  \
                             ta1_y_xxx_xxxxx_1,  \
                             ta1_y_xxx_xxxxy_0,  \
                             ta1_y_xxx_xxxxy_1,  \
                             ta1_y_xxx_xxxxz_0,  \
                             ta1_y_xxx_xxxxz_1,  \
                             ta1_y_xxx_xxxy_0,   \
                             ta1_y_xxx_xxxy_1,   \
                             ta1_y_xxx_xxxyy_0,  \
                             ta1_y_xxx_xxxyy_1,  \
                             ta1_y_xxx_xxxyz_0,  \
                             ta1_y_xxx_xxxyz_1,  \
                             ta1_y_xxx_xxxz_0,   \
                             ta1_y_xxx_xxxz_1,   \
                             ta1_y_xxx_xxxzz_0,  \
                             ta1_y_xxx_xxxzz_1,  \
                             ta1_y_xxx_xxyy_0,   \
                             ta1_y_xxx_xxyy_1,   \
                             ta1_y_xxx_xxyyy_0,  \
                             ta1_y_xxx_xxyyy_1,  \
                             ta1_y_xxx_xxyyz_0,  \
                             ta1_y_xxx_xxyyz_1,  \
                             ta1_y_xxx_xxyz_0,   \
                             ta1_y_xxx_xxyz_1,   \
                             ta1_y_xxx_xxyzz_0,  \
                             ta1_y_xxx_xxyzz_1,  \
                             ta1_y_xxx_xxzz_0,   \
                             ta1_y_xxx_xxzz_1,   \
                             ta1_y_xxx_xxzzz_0,  \
                             ta1_y_xxx_xxzzz_1,  \
                             ta1_y_xxx_xyyy_0,   \
                             ta1_y_xxx_xyyy_1,   \
                             ta1_y_xxx_xyyyy_0,  \
                             ta1_y_xxx_xyyyy_1,  \
                             ta1_y_xxx_xyyyz_0,  \
                             ta1_y_xxx_xyyyz_1,  \
                             ta1_y_xxx_xyyz_0,   \
                             ta1_y_xxx_xyyz_1,   \
                             ta1_y_xxx_xyyzz_0,  \
                             ta1_y_xxx_xyyzz_1,  \
                             ta1_y_xxx_xyzz_0,   \
                             ta1_y_xxx_xyzz_1,   \
                             ta1_y_xxx_xyzzz_0,  \
                             ta1_y_xxx_xyzzz_1,  \
                             ta1_y_xxx_xzzz_0,   \
                             ta1_y_xxx_xzzz_1,   \
                             ta1_y_xxx_xzzzz_0,  \
                             ta1_y_xxx_xzzzz_1,  \
                             ta1_y_xxx_zzzzz_0,  \
                             ta1_y_xxx_zzzzz_1,  \
                             ta1_y_xxxy_xxxxx_0, \
                             ta1_y_xxxy_xxxxy_0, \
                             ta1_y_xxxy_xxxxz_0, \
                             ta1_y_xxxy_xxxyy_0, \
                             ta1_y_xxxy_xxxyz_0, \
                             ta1_y_xxxy_xxxzz_0, \
                             ta1_y_xxxy_xxyyy_0, \
                             ta1_y_xxxy_xxyyz_0, \
                             ta1_y_xxxy_xxyzz_0, \
                             ta1_y_xxxy_xxzzz_0, \
                             ta1_y_xxxy_xyyyy_0, \
                             ta1_y_xxxy_xyyyz_0, \
                             ta1_y_xxxy_xyyzz_0, \
                             ta1_y_xxxy_xyzzz_0, \
                             ta1_y_xxxy_xzzzz_0, \
                             ta1_y_xxxy_yyyyy_0, \
                             ta1_y_xxxy_yyyyz_0, \
                             ta1_y_xxxy_yyyzz_0, \
                             ta1_y_xxxy_yyzzz_0, \
                             ta1_y_xxxy_yzzzz_0, \
                             ta1_y_xxxy_zzzzz_0, \
                             ta1_y_xxy_yyyyy_0,  \
                             ta1_y_xxy_yyyyy_1,  \
                             ta1_y_xxy_yyyyz_0,  \
                             ta1_y_xxy_yyyyz_1,  \
                             ta1_y_xxy_yyyzz_0,  \
                             ta1_y_xxy_yyyzz_1,  \
                             ta1_y_xxy_yyzzz_0,  \
                             ta1_y_xxy_yyzzz_1,  \
                             ta1_y_xxy_yzzzz_0,  \
                             ta1_y_xxy_yzzzz_1,  \
                             ta1_y_xy_yyyyy_0,   \
                             ta1_y_xy_yyyyy_1,   \
                             ta1_y_xy_yyyyz_0,   \
                             ta1_y_xy_yyyyz_1,   \
                             ta1_y_xy_yyyzz_0,   \
                             ta1_y_xy_yyyzz_1,   \
                             ta1_y_xy_yyzzz_0,   \
                             ta1_y_xy_yyzzz_1,   \
                             ta1_y_xy_yzzzz_0,   \
                             ta1_y_xy_yzzzz_1,   \
                             ta_xxx_xxxxx_1,     \
                             ta_xxx_xxxxy_1,     \
                             ta_xxx_xxxxz_1,     \
                             ta_xxx_xxxyy_1,     \
                             ta_xxx_xxxyz_1,     \
                             ta_xxx_xxxzz_1,     \
                             ta_xxx_xxyyy_1,     \
                             ta_xxx_xxyyz_1,     \
                             ta_xxx_xxyzz_1,     \
                             ta_xxx_xxzzz_1,     \
                             ta_xxx_xyyyy_1,     \
                             ta_xxx_xyyyz_1,     \
                             ta_xxx_xyyzz_1,     \
                             ta_xxx_xyzzz_1,     \
                             ta_xxx_xzzzz_1,     \
                             ta_xxx_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxy_xxxxx_0[i] = ta_xxx_xxxxx_1[i] + ta1_y_xxx_xxxxx_0[i] * pa_y[i] - ta1_y_xxx_xxxxx_1[i] * pc_y[i];

        ta1_y_xxxy_xxxxy_0[i] = ta1_y_xxx_xxxx_0[i] * fe_0 - ta1_y_xxx_xxxx_1[i] * fe_0 + ta_xxx_xxxxy_1[i] + ta1_y_xxx_xxxxy_0[i] * pa_y[i] -
                                ta1_y_xxx_xxxxy_1[i] * pc_y[i];

        ta1_y_xxxy_xxxxz_0[i] = ta_xxx_xxxxz_1[i] + ta1_y_xxx_xxxxz_0[i] * pa_y[i] - ta1_y_xxx_xxxxz_1[i] * pc_y[i];

        ta1_y_xxxy_xxxyy_0[i] = 2.0 * ta1_y_xxx_xxxy_0[i] * fe_0 - 2.0 * ta1_y_xxx_xxxy_1[i] * fe_0 + ta_xxx_xxxyy_1[i] +
                                ta1_y_xxx_xxxyy_0[i] * pa_y[i] - ta1_y_xxx_xxxyy_1[i] * pc_y[i];

        ta1_y_xxxy_xxxyz_0[i] = ta1_y_xxx_xxxz_0[i] * fe_0 - ta1_y_xxx_xxxz_1[i] * fe_0 + ta_xxx_xxxyz_1[i] + ta1_y_xxx_xxxyz_0[i] * pa_y[i] -
                                ta1_y_xxx_xxxyz_1[i] * pc_y[i];

        ta1_y_xxxy_xxxzz_0[i] = ta_xxx_xxxzz_1[i] + ta1_y_xxx_xxxzz_0[i] * pa_y[i] - ta1_y_xxx_xxxzz_1[i] * pc_y[i];

        ta1_y_xxxy_xxyyy_0[i] = 3.0 * ta1_y_xxx_xxyy_0[i] * fe_0 - 3.0 * ta1_y_xxx_xxyy_1[i] * fe_0 + ta_xxx_xxyyy_1[i] +
                                ta1_y_xxx_xxyyy_0[i] * pa_y[i] - ta1_y_xxx_xxyyy_1[i] * pc_y[i];

        ta1_y_xxxy_xxyyz_0[i] = 2.0 * ta1_y_xxx_xxyz_0[i] * fe_0 - 2.0 * ta1_y_xxx_xxyz_1[i] * fe_0 + ta_xxx_xxyyz_1[i] +
                                ta1_y_xxx_xxyyz_0[i] * pa_y[i] - ta1_y_xxx_xxyyz_1[i] * pc_y[i];

        ta1_y_xxxy_xxyzz_0[i] = ta1_y_xxx_xxzz_0[i] * fe_0 - ta1_y_xxx_xxzz_1[i] * fe_0 + ta_xxx_xxyzz_1[i] + ta1_y_xxx_xxyzz_0[i] * pa_y[i] -
                                ta1_y_xxx_xxyzz_1[i] * pc_y[i];

        ta1_y_xxxy_xxzzz_0[i] = ta_xxx_xxzzz_1[i] + ta1_y_xxx_xxzzz_0[i] * pa_y[i] - ta1_y_xxx_xxzzz_1[i] * pc_y[i];

        ta1_y_xxxy_xyyyy_0[i] = 4.0 * ta1_y_xxx_xyyy_0[i] * fe_0 - 4.0 * ta1_y_xxx_xyyy_1[i] * fe_0 + ta_xxx_xyyyy_1[i] +
                                ta1_y_xxx_xyyyy_0[i] * pa_y[i] - ta1_y_xxx_xyyyy_1[i] * pc_y[i];

        ta1_y_xxxy_xyyyz_0[i] = 3.0 * ta1_y_xxx_xyyz_0[i] * fe_0 - 3.0 * ta1_y_xxx_xyyz_1[i] * fe_0 + ta_xxx_xyyyz_1[i] +
                                ta1_y_xxx_xyyyz_0[i] * pa_y[i] - ta1_y_xxx_xyyyz_1[i] * pc_y[i];

        ta1_y_xxxy_xyyzz_0[i] = 2.0 * ta1_y_xxx_xyzz_0[i] * fe_0 - 2.0 * ta1_y_xxx_xyzz_1[i] * fe_0 + ta_xxx_xyyzz_1[i] +
                                ta1_y_xxx_xyyzz_0[i] * pa_y[i] - ta1_y_xxx_xyyzz_1[i] * pc_y[i];

        ta1_y_xxxy_xyzzz_0[i] = ta1_y_xxx_xzzz_0[i] * fe_0 - ta1_y_xxx_xzzz_1[i] * fe_0 + ta_xxx_xyzzz_1[i] + ta1_y_xxx_xyzzz_0[i] * pa_y[i] -
                                ta1_y_xxx_xyzzz_1[i] * pc_y[i];

        ta1_y_xxxy_xzzzz_0[i] = ta_xxx_xzzzz_1[i] + ta1_y_xxx_xzzzz_0[i] * pa_y[i] - ta1_y_xxx_xzzzz_1[i] * pc_y[i];

        ta1_y_xxxy_yyyyy_0[i] =
            2.0 * ta1_y_xy_yyyyy_0[i] * fe_0 - 2.0 * ta1_y_xy_yyyyy_1[i] * fe_0 + ta1_y_xxy_yyyyy_0[i] * pa_x[i] - ta1_y_xxy_yyyyy_1[i] * pc_x[i];

        ta1_y_xxxy_yyyyz_0[i] =
            2.0 * ta1_y_xy_yyyyz_0[i] * fe_0 - 2.0 * ta1_y_xy_yyyyz_1[i] * fe_0 + ta1_y_xxy_yyyyz_0[i] * pa_x[i] - ta1_y_xxy_yyyyz_1[i] * pc_x[i];

        ta1_y_xxxy_yyyzz_0[i] =
            2.0 * ta1_y_xy_yyyzz_0[i] * fe_0 - 2.0 * ta1_y_xy_yyyzz_1[i] * fe_0 + ta1_y_xxy_yyyzz_0[i] * pa_x[i] - ta1_y_xxy_yyyzz_1[i] * pc_x[i];

        ta1_y_xxxy_yyzzz_0[i] =
            2.0 * ta1_y_xy_yyzzz_0[i] * fe_0 - 2.0 * ta1_y_xy_yyzzz_1[i] * fe_0 + ta1_y_xxy_yyzzz_0[i] * pa_x[i] - ta1_y_xxy_yyzzz_1[i] * pc_x[i];

        ta1_y_xxxy_yzzzz_0[i] =
            2.0 * ta1_y_xy_yzzzz_0[i] * fe_0 - 2.0 * ta1_y_xy_yzzzz_1[i] * fe_0 + ta1_y_xxy_yzzzz_0[i] * pa_x[i] - ta1_y_xxy_yzzzz_1[i] * pc_x[i];

        ta1_y_xxxy_zzzzz_0[i] = ta_xxx_zzzzz_1[i] + ta1_y_xxx_zzzzz_0[i] * pa_y[i] - ta1_y_xxx_zzzzz_1[i] * pc_y[i];
    }

    // Set up 357-378 components of targeted buffer : GH

    auto ta1_y_xxxz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 357);

    auto ta1_y_xxxz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 358);

    auto ta1_y_xxxz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 359);

    auto ta1_y_xxxz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 360);

    auto ta1_y_xxxz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 361);

    auto ta1_y_xxxz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 362);

    auto ta1_y_xxxz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 363);

    auto ta1_y_xxxz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 364);

    auto ta1_y_xxxz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 365);

    auto ta1_y_xxxz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 366);

    auto ta1_y_xxxz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 367);

    auto ta1_y_xxxz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 368);

    auto ta1_y_xxxz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 369);

    auto ta1_y_xxxz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 370);

    auto ta1_y_xxxz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 371);

    auto ta1_y_xxxz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 372);

    auto ta1_y_xxxz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 373);

    auto ta1_y_xxxz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 374);

    auto ta1_y_xxxz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 375);

    auto ta1_y_xxxz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 376);

    auto ta1_y_xxxz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 377);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xxx_xxxx_0,   \
                             ta1_y_xxx_xxxx_1,   \
                             ta1_y_xxx_xxxxx_0,  \
                             ta1_y_xxx_xxxxx_1,  \
                             ta1_y_xxx_xxxxy_0,  \
                             ta1_y_xxx_xxxxy_1,  \
                             ta1_y_xxx_xxxxz_0,  \
                             ta1_y_xxx_xxxxz_1,  \
                             ta1_y_xxx_xxxy_0,   \
                             ta1_y_xxx_xxxy_1,   \
                             ta1_y_xxx_xxxyy_0,  \
                             ta1_y_xxx_xxxyy_1,  \
                             ta1_y_xxx_xxxyz_0,  \
                             ta1_y_xxx_xxxyz_1,  \
                             ta1_y_xxx_xxxz_0,   \
                             ta1_y_xxx_xxxz_1,   \
                             ta1_y_xxx_xxxzz_0,  \
                             ta1_y_xxx_xxxzz_1,  \
                             ta1_y_xxx_xxyy_0,   \
                             ta1_y_xxx_xxyy_1,   \
                             ta1_y_xxx_xxyyy_0,  \
                             ta1_y_xxx_xxyyy_1,  \
                             ta1_y_xxx_xxyyz_0,  \
                             ta1_y_xxx_xxyyz_1,  \
                             ta1_y_xxx_xxyz_0,   \
                             ta1_y_xxx_xxyz_1,   \
                             ta1_y_xxx_xxyzz_0,  \
                             ta1_y_xxx_xxyzz_1,  \
                             ta1_y_xxx_xxzz_0,   \
                             ta1_y_xxx_xxzz_1,   \
                             ta1_y_xxx_xxzzz_0,  \
                             ta1_y_xxx_xxzzz_1,  \
                             ta1_y_xxx_xyyy_0,   \
                             ta1_y_xxx_xyyy_1,   \
                             ta1_y_xxx_xyyyy_0,  \
                             ta1_y_xxx_xyyyy_1,  \
                             ta1_y_xxx_xyyyz_0,  \
                             ta1_y_xxx_xyyyz_1,  \
                             ta1_y_xxx_xyyz_0,   \
                             ta1_y_xxx_xyyz_1,   \
                             ta1_y_xxx_xyyzz_0,  \
                             ta1_y_xxx_xyyzz_1,  \
                             ta1_y_xxx_xyzz_0,   \
                             ta1_y_xxx_xyzz_1,   \
                             ta1_y_xxx_xyzzz_0,  \
                             ta1_y_xxx_xyzzz_1,  \
                             ta1_y_xxx_xzzz_0,   \
                             ta1_y_xxx_xzzz_1,   \
                             ta1_y_xxx_xzzzz_0,  \
                             ta1_y_xxx_xzzzz_1,  \
                             ta1_y_xxx_yyyyy_0,  \
                             ta1_y_xxx_yyyyy_1,  \
                             ta1_y_xxxz_xxxxx_0, \
                             ta1_y_xxxz_xxxxy_0, \
                             ta1_y_xxxz_xxxxz_0, \
                             ta1_y_xxxz_xxxyy_0, \
                             ta1_y_xxxz_xxxyz_0, \
                             ta1_y_xxxz_xxxzz_0, \
                             ta1_y_xxxz_xxyyy_0, \
                             ta1_y_xxxz_xxyyz_0, \
                             ta1_y_xxxz_xxyzz_0, \
                             ta1_y_xxxz_xxzzz_0, \
                             ta1_y_xxxz_xyyyy_0, \
                             ta1_y_xxxz_xyyyz_0, \
                             ta1_y_xxxz_xyyzz_0, \
                             ta1_y_xxxz_xyzzz_0, \
                             ta1_y_xxxz_xzzzz_0, \
                             ta1_y_xxxz_yyyyy_0, \
                             ta1_y_xxxz_yyyyz_0, \
                             ta1_y_xxxz_yyyzz_0, \
                             ta1_y_xxxz_yyzzz_0, \
                             ta1_y_xxxz_yzzzz_0, \
                             ta1_y_xxxz_zzzzz_0, \
                             ta1_y_xxz_yyyyz_0,  \
                             ta1_y_xxz_yyyyz_1,  \
                             ta1_y_xxz_yyyzz_0,  \
                             ta1_y_xxz_yyyzz_1,  \
                             ta1_y_xxz_yyzzz_0,  \
                             ta1_y_xxz_yyzzz_1,  \
                             ta1_y_xxz_yzzzz_0,  \
                             ta1_y_xxz_yzzzz_1,  \
                             ta1_y_xxz_zzzzz_0,  \
                             ta1_y_xxz_zzzzz_1,  \
                             ta1_y_xz_yyyyz_0,   \
                             ta1_y_xz_yyyyz_1,   \
                             ta1_y_xz_yyyzz_0,   \
                             ta1_y_xz_yyyzz_1,   \
                             ta1_y_xz_yyzzz_0,   \
                             ta1_y_xz_yyzzz_1,   \
                             ta1_y_xz_yzzzz_0,   \
                             ta1_y_xz_yzzzz_1,   \
                             ta1_y_xz_zzzzz_0,   \
                             ta1_y_xz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxxz_xxxxx_0[i] = ta1_y_xxx_xxxxx_0[i] * pa_z[i] - ta1_y_xxx_xxxxx_1[i] * pc_z[i];

        ta1_y_xxxz_xxxxy_0[i] = ta1_y_xxx_xxxxy_0[i] * pa_z[i] - ta1_y_xxx_xxxxy_1[i] * pc_z[i];

        ta1_y_xxxz_xxxxz_0[i] =
            ta1_y_xxx_xxxx_0[i] * fe_0 - ta1_y_xxx_xxxx_1[i] * fe_0 + ta1_y_xxx_xxxxz_0[i] * pa_z[i] - ta1_y_xxx_xxxxz_1[i] * pc_z[i];

        ta1_y_xxxz_xxxyy_0[i] = ta1_y_xxx_xxxyy_0[i] * pa_z[i] - ta1_y_xxx_xxxyy_1[i] * pc_z[i];

        ta1_y_xxxz_xxxyz_0[i] =
            ta1_y_xxx_xxxy_0[i] * fe_0 - ta1_y_xxx_xxxy_1[i] * fe_0 + ta1_y_xxx_xxxyz_0[i] * pa_z[i] - ta1_y_xxx_xxxyz_1[i] * pc_z[i];

        ta1_y_xxxz_xxxzz_0[i] =
            2.0 * ta1_y_xxx_xxxz_0[i] * fe_0 - 2.0 * ta1_y_xxx_xxxz_1[i] * fe_0 + ta1_y_xxx_xxxzz_0[i] * pa_z[i] - ta1_y_xxx_xxxzz_1[i] * pc_z[i];

        ta1_y_xxxz_xxyyy_0[i] = ta1_y_xxx_xxyyy_0[i] * pa_z[i] - ta1_y_xxx_xxyyy_1[i] * pc_z[i];

        ta1_y_xxxz_xxyyz_0[i] =
            ta1_y_xxx_xxyy_0[i] * fe_0 - ta1_y_xxx_xxyy_1[i] * fe_0 + ta1_y_xxx_xxyyz_0[i] * pa_z[i] - ta1_y_xxx_xxyyz_1[i] * pc_z[i];

        ta1_y_xxxz_xxyzz_0[i] =
            2.0 * ta1_y_xxx_xxyz_0[i] * fe_0 - 2.0 * ta1_y_xxx_xxyz_1[i] * fe_0 + ta1_y_xxx_xxyzz_0[i] * pa_z[i] - ta1_y_xxx_xxyzz_1[i] * pc_z[i];

        ta1_y_xxxz_xxzzz_0[i] =
            3.0 * ta1_y_xxx_xxzz_0[i] * fe_0 - 3.0 * ta1_y_xxx_xxzz_1[i] * fe_0 + ta1_y_xxx_xxzzz_0[i] * pa_z[i] - ta1_y_xxx_xxzzz_1[i] * pc_z[i];

        ta1_y_xxxz_xyyyy_0[i] = ta1_y_xxx_xyyyy_0[i] * pa_z[i] - ta1_y_xxx_xyyyy_1[i] * pc_z[i];

        ta1_y_xxxz_xyyyz_0[i] =
            ta1_y_xxx_xyyy_0[i] * fe_0 - ta1_y_xxx_xyyy_1[i] * fe_0 + ta1_y_xxx_xyyyz_0[i] * pa_z[i] - ta1_y_xxx_xyyyz_1[i] * pc_z[i];

        ta1_y_xxxz_xyyzz_0[i] =
            2.0 * ta1_y_xxx_xyyz_0[i] * fe_0 - 2.0 * ta1_y_xxx_xyyz_1[i] * fe_0 + ta1_y_xxx_xyyzz_0[i] * pa_z[i] - ta1_y_xxx_xyyzz_1[i] * pc_z[i];

        ta1_y_xxxz_xyzzz_0[i] =
            3.0 * ta1_y_xxx_xyzz_0[i] * fe_0 - 3.0 * ta1_y_xxx_xyzz_1[i] * fe_0 + ta1_y_xxx_xyzzz_0[i] * pa_z[i] - ta1_y_xxx_xyzzz_1[i] * pc_z[i];

        ta1_y_xxxz_xzzzz_0[i] =
            4.0 * ta1_y_xxx_xzzz_0[i] * fe_0 - 4.0 * ta1_y_xxx_xzzz_1[i] * fe_0 + ta1_y_xxx_xzzzz_0[i] * pa_z[i] - ta1_y_xxx_xzzzz_1[i] * pc_z[i];

        ta1_y_xxxz_yyyyy_0[i] = ta1_y_xxx_yyyyy_0[i] * pa_z[i] - ta1_y_xxx_yyyyy_1[i] * pc_z[i];

        ta1_y_xxxz_yyyyz_0[i] =
            2.0 * ta1_y_xz_yyyyz_0[i] * fe_0 - 2.0 * ta1_y_xz_yyyyz_1[i] * fe_0 + ta1_y_xxz_yyyyz_0[i] * pa_x[i] - ta1_y_xxz_yyyyz_1[i] * pc_x[i];

        ta1_y_xxxz_yyyzz_0[i] =
            2.0 * ta1_y_xz_yyyzz_0[i] * fe_0 - 2.0 * ta1_y_xz_yyyzz_1[i] * fe_0 + ta1_y_xxz_yyyzz_0[i] * pa_x[i] - ta1_y_xxz_yyyzz_1[i] * pc_x[i];

        ta1_y_xxxz_yyzzz_0[i] =
            2.0 * ta1_y_xz_yyzzz_0[i] * fe_0 - 2.0 * ta1_y_xz_yyzzz_1[i] * fe_0 + ta1_y_xxz_yyzzz_0[i] * pa_x[i] - ta1_y_xxz_yyzzz_1[i] * pc_x[i];

        ta1_y_xxxz_yzzzz_0[i] =
            2.0 * ta1_y_xz_yzzzz_0[i] * fe_0 - 2.0 * ta1_y_xz_yzzzz_1[i] * fe_0 + ta1_y_xxz_yzzzz_0[i] * pa_x[i] - ta1_y_xxz_yzzzz_1[i] * pc_x[i];

        ta1_y_xxxz_zzzzz_0[i] =
            2.0 * ta1_y_xz_zzzzz_0[i] * fe_0 - 2.0 * ta1_y_xz_zzzzz_1[i] * fe_0 + ta1_y_xxz_zzzzz_0[i] * pa_x[i] - ta1_y_xxz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 378-399 components of targeted buffer : GH

    auto ta1_y_xxyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 378);

    auto ta1_y_xxyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 379);

    auto ta1_y_xxyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 380);

    auto ta1_y_xxyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 381);

    auto ta1_y_xxyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 382);

    auto ta1_y_xxyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 383);

    auto ta1_y_xxyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 384);

    auto ta1_y_xxyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 385);

    auto ta1_y_xxyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 386);

    auto ta1_y_xxyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 387);

    auto ta1_y_xxyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 388);

    auto ta1_y_xxyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 389);

    auto ta1_y_xxyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 390);

    auto ta1_y_xxyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 391);

    auto ta1_y_xxyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 392);

    auto ta1_y_xxyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 393);

    auto ta1_y_xxyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 394);

    auto ta1_y_xxyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 395);

    auto ta1_y_xxyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 396);

    auto ta1_y_xxyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 397);

    auto ta1_y_xxyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 398);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_y_xx_xxxxx_0,   \
                             ta1_y_xx_xxxxx_1,   \
                             ta1_y_xx_xxxxz_0,   \
                             ta1_y_xx_xxxxz_1,   \
                             ta1_y_xx_xxxzz_0,   \
                             ta1_y_xx_xxxzz_1,   \
                             ta1_y_xx_xxzzz_0,   \
                             ta1_y_xx_xxzzz_1,   \
                             ta1_y_xx_xzzzz_0,   \
                             ta1_y_xx_xzzzz_1,   \
                             ta1_y_xxy_xxxxx_0,  \
                             ta1_y_xxy_xxxxx_1,  \
                             ta1_y_xxy_xxxxz_0,  \
                             ta1_y_xxy_xxxxz_1,  \
                             ta1_y_xxy_xxxzz_0,  \
                             ta1_y_xxy_xxxzz_1,  \
                             ta1_y_xxy_xxzzz_0,  \
                             ta1_y_xxy_xxzzz_1,  \
                             ta1_y_xxy_xzzzz_0,  \
                             ta1_y_xxy_xzzzz_1,  \
                             ta1_y_xxyy_xxxxx_0, \
                             ta1_y_xxyy_xxxxy_0, \
                             ta1_y_xxyy_xxxxz_0, \
                             ta1_y_xxyy_xxxyy_0, \
                             ta1_y_xxyy_xxxyz_0, \
                             ta1_y_xxyy_xxxzz_0, \
                             ta1_y_xxyy_xxyyy_0, \
                             ta1_y_xxyy_xxyyz_0, \
                             ta1_y_xxyy_xxyzz_0, \
                             ta1_y_xxyy_xxzzz_0, \
                             ta1_y_xxyy_xyyyy_0, \
                             ta1_y_xxyy_xyyyz_0, \
                             ta1_y_xxyy_xyyzz_0, \
                             ta1_y_xxyy_xyzzz_0, \
                             ta1_y_xxyy_xzzzz_0, \
                             ta1_y_xxyy_yyyyy_0, \
                             ta1_y_xxyy_yyyyz_0, \
                             ta1_y_xxyy_yyyzz_0, \
                             ta1_y_xxyy_yyzzz_0, \
                             ta1_y_xxyy_yzzzz_0, \
                             ta1_y_xxyy_zzzzz_0, \
                             ta1_y_xyy_xxxxy_0,  \
                             ta1_y_xyy_xxxxy_1,  \
                             ta1_y_xyy_xxxy_0,   \
                             ta1_y_xyy_xxxy_1,   \
                             ta1_y_xyy_xxxyy_0,  \
                             ta1_y_xyy_xxxyy_1,  \
                             ta1_y_xyy_xxxyz_0,  \
                             ta1_y_xyy_xxxyz_1,  \
                             ta1_y_xyy_xxyy_0,   \
                             ta1_y_xyy_xxyy_1,   \
                             ta1_y_xyy_xxyyy_0,  \
                             ta1_y_xyy_xxyyy_1,  \
                             ta1_y_xyy_xxyyz_0,  \
                             ta1_y_xyy_xxyyz_1,  \
                             ta1_y_xyy_xxyz_0,   \
                             ta1_y_xyy_xxyz_1,   \
                             ta1_y_xyy_xxyzz_0,  \
                             ta1_y_xyy_xxyzz_1,  \
                             ta1_y_xyy_xyyy_0,   \
                             ta1_y_xyy_xyyy_1,   \
                             ta1_y_xyy_xyyyy_0,  \
                             ta1_y_xyy_xyyyy_1,  \
                             ta1_y_xyy_xyyyz_0,  \
                             ta1_y_xyy_xyyyz_1,  \
                             ta1_y_xyy_xyyz_0,   \
                             ta1_y_xyy_xyyz_1,   \
                             ta1_y_xyy_xyyzz_0,  \
                             ta1_y_xyy_xyyzz_1,  \
                             ta1_y_xyy_xyzz_0,   \
                             ta1_y_xyy_xyzz_1,   \
                             ta1_y_xyy_xyzzz_0,  \
                             ta1_y_xyy_xyzzz_1,  \
                             ta1_y_xyy_yyyy_0,   \
                             ta1_y_xyy_yyyy_1,   \
                             ta1_y_xyy_yyyyy_0,  \
                             ta1_y_xyy_yyyyy_1,  \
                             ta1_y_xyy_yyyyz_0,  \
                             ta1_y_xyy_yyyyz_1,  \
                             ta1_y_xyy_yyyz_0,   \
                             ta1_y_xyy_yyyz_1,   \
                             ta1_y_xyy_yyyzz_0,  \
                             ta1_y_xyy_yyyzz_1,  \
                             ta1_y_xyy_yyzz_0,   \
                             ta1_y_xyy_yyzz_1,   \
                             ta1_y_xyy_yyzzz_0,  \
                             ta1_y_xyy_yyzzz_1,  \
                             ta1_y_xyy_yzzz_0,   \
                             ta1_y_xyy_yzzz_1,   \
                             ta1_y_xyy_yzzzz_0,  \
                             ta1_y_xyy_yzzzz_1,  \
                             ta1_y_xyy_zzzzz_0,  \
                             ta1_y_xyy_zzzzz_1,  \
                             ta1_y_yy_xxxxy_0,   \
                             ta1_y_yy_xxxxy_1,   \
                             ta1_y_yy_xxxyy_0,   \
                             ta1_y_yy_xxxyy_1,   \
                             ta1_y_yy_xxxyz_0,   \
                             ta1_y_yy_xxxyz_1,   \
                             ta1_y_yy_xxyyy_0,   \
                             ta1_y_yy_xxyyy_1,   \
                             ta1_y_yy_xxyyz_0,   \
                             ta1_y_yy_xxyyz_1,   \
                             ta1_y_yy_xxyzz_0,   \
                             ta1_y_yy_xxyzz_1,   \
                             ta1_y_yy_xyyyy_0,   \
                             ta1_y_yy_xyyyy_1,   \
                             ta1_y_yy_xyyyz_0,   \
                             ta1_y_yy_xyyyz_1,   \
                             ta1_y_yy_xyyzz_0,   \
                             ta1_y_yy_xyyzz_1,   \
                             ta1_y_yy_xyzzz_0,   \
                             ta1_y_yy_xyzzz_1,   \
                             ta1_y_yy_yyyyy_0,   \
                             ta1_y_yy_yyyyy_1,   \
                             ta1_y_yy_yyyyz_0,   \
                             ta1_y_yy_yyyyz_1,   \
                             ta1_y_yy_yyyzz_0,   \
                             ta1_y_yy_yyyzz_1,   \
                             ta1_y_yy_yyzzz_0,   \
                             ta1_y_yy_yyzzz_1,   \
                             ta1_y_yy_yzzzz_0,   \
                             ta1_y_yy_yzzzz_1,   \
                             ta1_y_yy_zzzzz_0,   \
                             ta1_y_yy_zzzzz_1,   \
                             ta_xxy_xxxxx_1,     \
                             ta_xxy_xxxxz_1,     \
                             ta_xxy_xxxzz_1,     \
                             ta_xxy_xxzzz_1,     \
                             ta_xxy_xzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyy_xxxxx_0[i] = ta1_y_xx_xxxxx_0[i] * fe_0 - ta1_y_xx_xxxxx_1[i] * fe_0 + ta_xxy_xxxxx_1[i] + ta1_y_xxy_xxxxx_0[i] * pa_y[i] -
                                ta1_y_xxy_xxxxx_1[i] * pc_y[i];

        ta1_y_xxyy_xxxxy_0[i] = ta1_y_yy_xxxxy_0[i] * fe_0 - ta1_y_yy_xxxxy_1[i] * fe_0 + 4.0 * ta1_y_xyy_xxxy_0[i] * fe_0 -
                                4.0 * ta1_y_xyy_xxxy_1[i] * fe_0 + ta1_y_xyy_xxxxy_0[i] * pa_x[i] - ta1_y_xyy_xxxxy_1[i] * pc_x[i];

        ta1_y_xxyy_xxxxz_0[i] = ta1_y_xx_xxxxz_0[i] * fe_0 - ta1_y_xx_xxxxz_1[i] * fe_0 + ta_xxy_xxxxz_1[i] + ta1_y_xxy_xxxxz_0[i] * pa_y[i] -
                                ta1_y_xxy_xxxxz_1[i] * pc_y[i];

        ta1_y_xxyy_xxxyy_0[i] = ta1_y_yy_xxxyy_0[i] * fe_0 - ta1_y_yy_xxxyy_1[i] * fe_0 + 3.0 * ta1_y_xyy_xxyy_0[i] * fe_0 -
                                3.0 * ta1_y_xyy_xxyy_1[i] * fe_0 + ta1_y_xyy_xxxyy_0[i] * pa_x[i] - ta1_y_xyy_xxxyy_1[i] * pc_x[i];

        ta1_y_xxyy_xxxyz_0[i] = ta1_y_yy_xxxyz_0[i] * fe_0 - ta1_y_yy_xxxyz_1[i] * fe_0 + 3.0 * ta1_y_xyy_xxyz_0[i] * fe_0 -
                                3.0 * ta1_y_xyy_xxyz_1[i] * fe_0 + ta1_y_xyy_xxxyz_0[i] * pa_x[i] - ta1_y_xyy_xxxyz_1[i] * pc_x[i];

        ta1_y_xxyy_xxxzz_0[i] = ta1_y_xx_xxxzz_0[i] * fe_0 - ta1_y_xx_xxxzz_1[i] * fe_0 + ta_xxy_xxxzz_1[i] + ta1_y_xxy_xxxzz_0[i] * pa_y[i] -
                                ta1_y_xxy_xxxzz_1[i] * pc_y[i];

        ta1_y_xxyy_xxyyy_0[i] = ta1_y_yy_xxyyy_0[i] * fe_0 - ta1_y_yy_xxyyy_1[i] * fe_0 + 2.0 * ta1_y_xyy_xyyy_0[i] * fe_0 -
                                2.0 * ta1_y_xyy_xyyy_1[i] * fe_0 + ta1_y_xyy_xxyyy_0[i] * pa_x[i] - ta1_y_xyy_xxyyy_1[i] * pc_x[i];

        ta1_y_xxyy_xxyyz_0[i] = ta1_y_yy_xxyyz_0[i] * fe_0 - ta1_y_yy_xxyyz_1[i] * fe_0 + 2.0 * ta1_y_xyy_xyyz_0[i] * fe_0 -
                                2.0 * ta1_y_xyy_xyyz_1[i] * fe_0 + ta1_y_xyy_xxyyz_0[i] * pa_x[i] - ta1_y_xyy_xxyyz_1[i] * pc_x[i];

        ta1_y_xxyy_xxyzz_0[i] = ta1_y_yy_xxyzz_0[i] * fe_0 - ta1_y_yy_xxyzz_1[i] * fe_0 + 2.0 * ta1_y_xyy_xyzz_0[i] * fe_0 -
                                2.0 * ta1_y_xyy_xyzz_1[i] * fe_0 + ta1_y_xyy_xxyzz_0[i] * pa_x[i] - ta1_y_xyy_xxyzz_1[i] * pc_x[i];

        ta1_y_xxyy_xxzzz_0[i] = ta1_y_xx_xxzzz_0[i] * fe_0 - ta1_y_xx_xxzzz_1[i] * fe_0 + ta_xxy_xxzzz_1[i] + ta1_y_xxy_xxzzz_0[i] * pa_y[i] -
                                ta1_y_xxy_xxzzz_1[i] * pc_y[i];

        ta1_y_xxyy_xyyyy_0[i] = ta1_y_yy_xyyyy_0[i] * fe_0 - ta1_y_yy_xyyyy_1[i] * fe_0 + ta1_y_xyy_yyyy_0[i] * fe_0 - ta1_y_xyy_yyyy_1[i] * fe_0 +
                                ta1_y_xyy_xyyyy_0[i] * pa_x[i] - ta1_y_xyy_xyyyy_1[i] * pc_x[i];

        ta1_y_xxyy_xyyyz_0[i] = ta1_y_yy_xyyyz_0[i] * fe_0 - ta1_y_yy_xyyyz_1[i] * fe_0 + ta1_y_xyy_yyyz_0[i] * fe_0 - ta1_y_xyy_yyyz_1[i] * fe_0 +
                                ta1_y_xyy_xyyyz_0[i] * pa_x[i] - ta1_y_xyy_xyyyz_1[i] * pc_x[i];

        ta1_y_xxyy_xyyzz_0[i] = ta1_y_yy_xyyzz_0[i] * fe_0 - ta1_y_yy_xyyzz_1[i] * fe_0 + ta1_y_xyy_yyzz_0[i] * fe_0 - ta1_y_xyy_yyzz_1[i] * fe_0 +
                                ta1_y_xyy_xyyzz_0[i] * pa_x[i] - ta1_y_xyy_xyyzz_1[i] * pc_x[i];

        ta1_y_xxyy_xyzzz_0[i] = ta1_y_yy_xyzzz_0[i] * fe_0 - ta1_y_yy_xyzzz_1[i] * fe_0 + ta1_y_xyy_yzzz_0[i] * fe_0 - ta1_y_xyy_yzzz_1[i] * fe_0 +
                                ta1_y_xyy_xyzzz_0[i] * pa_x[i] - ta1_y_xyy_xyzzz_1[i] * pc_x[i];

        ta1_y_xxyy_xzzzz_0[i] = ta1_y_xx_xzzzz_0[i] * fe_0 - ta1_y_xx_xzzzz_1[i] * fe_0 + ta_xxy_xzzzz_1[i] + ta1_y_xxy_xzzzz_0[i] * pa_y[i] -
                                ta1_y_xxy_xzzzz_1[i] * pc_y[i];

        ta1_y_xxyy_yyyyy_0[i] =
            ta1_y_yy_yyyyy_0[i] * fe_0 - ta1_y_yy_yyyyy_1[i] * fe_0 + ta1_y_xyy_yyyyy_0[i] * pa_x[i] - ta1_y_xyy_yyyyy_1[i] * pc_x[i];

        ta1_y_xxyy_yyyyz_0[i] =
            ta1_y_yy_yyyyz_0[i] * fe_0 - ta1_y_yy_yyyyz_1[i] * fe_0 + ta1_y_xyy_yyyyz_0[i] * pa_x[i] - ta1_y_xyy_yyyyz_1[i] * pc_x[i];

        ta1_y_xxyy_yyyzz_0[i] =
            ta1_y_yy_yyyzz_0[i] * fe_0 - ta1_y_yy_yyyzz_1[i] * fe_0 + ta1_y_xyy_yyyzz_0[i] * pa_x[i] - ta1_y_xyy_yyyzz_1[i] * pc_x[i];

        ta1_y_xxyy_yyzzz_0[i] =
            ta1_y_yy_yyzzz_0[i] * fe_0 - ta1_y_yy_yyzzz_1[i] * fe_0 + ta1_y_xyy_yyzzz_0[i] * pa_x[i] - ta1_y_xyy_yyzzz_1[i] * pc_x[i];

        ta1_y_xxyy_yzzzz_0[i] =
            ta1_y_yy_yzzzz_0[i] * fe_0 - ta1_y_yy_yzzzz_1[i] * fe_0 + ta1_y_xyy_yzzzz_0[i] * pa_x[i] - ta1_y_xyy_yzzzz_1[i] * pc_x[i];

        ta1_y_xxyy_zzzzz_0[i] =
            ta1_y_yy_zzzzz_0[i] * fe_0 - ta1_y_yy_zzzzz_1[i] * fe_0 + ta1_y_xyy_zzzzz_0[i] * pa_x[i] - ta1_y_xyy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 399-420 components of targeted buffer : GH

    auto ta1_y_xxyz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 399);

    auto ta1_y_xxyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 400);

    auto ta1_y_xxyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 401);

    auto ta1_y_xxyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 402);

    auto ta1_y_xxyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 403);

    auto ta1_y_xxyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 404);

    auto ta1_y_xxyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 405);

    auto ta1_y_xxyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 406);

    auto ta1_y_xxyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 407);

    auto ta1_y_xxyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 408);

    auto ta1_y_xxyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 409);

    auto ta1_y_xxyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 410);

    auto ta1_y_xxyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 411);

    auto ta1_y_xxyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 412);

    auto ta1_y_xxyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 413);

    auto ta1_y_xxyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 414);

    auto ta1_y_xxyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 415);

    auto ta1_y_xxyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 416);

    auto ta1_y_xxyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 417);

    auto ta1_y_xxyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 418);

    auto ta1_y_xxyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 419);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_xxy_xxxxx_0,  \
                             ta1_y_xxy_xxxxx_1,  \
                             ta1_y_xxy_xxxxy_0,  \
                             ta1_y_xxy_xxxxy_1,  \
                             ta1_y_xxy_xxxy_0,   \
                             ta1_y_xxy_xxxy_1,   \
                             ta1_y_xxy_xxxyy_0,  \
                             ta1_y_xxy_xxxyy_1,  \
                             ta1_y_xxy_xxxyz_0,  \
                             ta1_y_xxy_xxxyz_1,  \
                             ta1_y_xxy_xxyy_0,   \
                             ta1_y_xxy_xxyy_1,   \
                             ta1_y_xxy_xxyyy_0,  \
                             ta1_y_xxy_xxyyy_1,  \
                             ta1_y_xxy_xxyyz_0,  \
                             ta1_y_xxy_xxyyz_1,  \
                             ta1_y_xxy_xxyz_0,   \
                             ta1_y_xxy_xxyz_1,   \
                             ta1_y_xxy_xxyzz_0,  \
                             ta1_y_xxy_xxyzz_1,  \
                             ta1_y_xxy_xyyy_0,   \
                             ta1_y_xxy_xyyy_1,   \
                             ta1_y_xxy_xyyyy_0,  \
                             ta1_y_xxy_xyyyy_1,  \
                             ta1_y_xxy_xyyyz_0,  \
                             ta1_y_xxy_xyyyz_1,  \
                             ta1_y_xxy_xyyz_0,   \
                             ta1_y_xxy_xyyz_1,   \
                             ta1_y_xxy_xyyzz_0,  \
                             ta1_y_xxy_xyyzz_1,  \
                             ta1_y_xxy_xyzz_0,   \
                             ta1_y_xxy_xyzz_1,   \
                             ta1_y_xxy_xyzzz_0,  \
                             ta1_y_xxy_xyzzz_1,  \
                             ta1_y_xxy_yyyyy_0,  \
                             ta1_y_xxy_yyyyy_1,  \
                             ta1_y_xxyz_xxxxx_0, \
                             ta1_y_xxyz_xxxxy_0, \
                             ta1_y_xxyz_xxxxz_0, \
                             ta1_y_xxyz_xxxyy_0, \
                             ta1_y_xxyz_xxxyz_0, \
                             ta1_y_xxyz_xxxzz_0, \
                             ta1_y_xxyz_xxyyy_0, \
                             ta1_y_xxyz_xxyyz_0, \
                             ta1_y_xxyz_xxyzz_0, \
                             ta1_y_xxyz_xxzzz_0, \
                             ta1_y_xxyz_xyyyy_0, \
                             ta1_y_xxyz_xyyyz_0, \
                             ta1_y_xxyz_xyyzz_0, \
                             ta1_y_xxyz_xyzzz_0, \
                             ta1_y_xxyz_xzzzz_0, \
                             ta1_y_xxyz_yyyyy_0, \
                             ta1_y_xxyz_yyyyz_0, \
                             ta1_y_xxyz_yyyzz_0, \
                             ta1_y_xxyz_yyzzz_0, \
                             ta1_y_xxyz_yzzzz_0, \
                             ta1_y_xxyz_zzzzz_0, \
                             ta1_y_xxz_xxxxz_0,  \
                             ta1_y_xxz_xxxxz_1,  \
                             ta1_y_xxz_xxxzz_0,  \
                             ta1_y_xxz_xxxzz_1,  \
                             ta1_y_xxz_xxzzz_0,  \
                             ta1_y_xxz_xxzzz_1,  \
                             ta1_y_xxz_xzzzz_0,  \
                             ta1_y_xxz_xzzzz_1,  \
                             ta1_y_xxz_zzzzz_0,  \
                             ta1_y_xxz_zzzzz_1,  \
                             ta1_y_xyz_yyyyz_0,  \
                             ta1_y_xyz_yyyyz_1,  \
                             ta1_y_xyz_yyyzz_0,  \
                             ta1_y_xyz_yyyzz_1,  \
                             ta1_y_xyz_yyzzz_0,  \
                             ta1_y_xyz_yyzzz_1,  \
                             ta1_y_xyz_yzzzz_0,  \
                             ta1_y_xyz_yzzzz_1,  \
                             ta1_y_yz_yyyyz_0,   \
                             ta1_y_yz_yyyyz_1,   \
                             ta1_y_yz_yyyzz_0,   \
                             ta1_y_yz_yyyzz_1,   \
                             ta1_y_yz_yyzzz_0,   \
                             ta1_y_yz_yyzzz_1,   \
                             ta1_y_yz_yzzzz_0,   \
                             ta1_y_yz_yzzzz_1,   \
                             ta_xxz_xxxxz_1,     \
                             ta_xxz_xxxzz_1,     \
                             ta_xxz_xxzzz_1,     \
                             ta_xxz_xzzzz_1,     \
                             ta_xxz_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxyz_xxxxx_0[i] = ta1_y_xxy_xxxxx_0[i] * pa_z[i] - ta1_y_xxy_xxxxx_1[i] * pc_z[i];

        ta1_y_xxyz_xxxxy_0[i] = ta1_y_xxy_xxxxy_0[i] * pa_z[i] - ta1_y_xxy_xxxxy_1[i] * pc_z[i];

        ta1_y_xxyz_xxxxz_0[i] = ta_xxz_xxxxz_1[i] + ta1_y_xxz_xxxxz_0[i] * pa_y[i] - ta1_y_xxz_xxxxz_1[i] * pc_y[i];

        ta1_y_xxyz_xxxyy_0[i] = ta1_y_xxy_xxxyy_0[i] * pa_z[i] - ta1_y_xxy_xxxyy_1[i] * pc_z[i];

        ta1_y_xxyz_xxxyz_0[i] =
            ta1_y_xxy_xxxy_0[i] * fe_0 - ta1_y_xxy_xxxy_1[i] * fe_0 + ta1_y_xxy_xxxyz_0[i] * pa_z[i] - ta1_y_xxy_xxxyz_1[i] * pc_z[i];

        ta1_y_xxyz_xxxzz_0[i] = ta_xxz_xxxzz_1[i] + ta1_y_xxz_xxxzz_0[i] * pa_y[i] - ta1_y_xxz_xxxzz_1[i] * pc_y[i];

        ta1_y_xxyz_xxyyy_0[i] = ta1_y_xxy_xxyyy_0[i] * pa_z[i] - ta1_y_xxy_xxyyy_1[i] * pc_z[i];

        ta1_y_xxyz_xxyyz_0[i] =
            ta1_y_xxy_xxyy_0[i] * fe_0 - ta1_y_xxy_xxyy_1[i] * fe_0 + ta1_y_xxy_xxyyz_0[i] * pa_z[i] - ta1_y_xxy_xxyyz_1[i] * pc_z[i];

        ta1_y_xxyz_xxyzz_0[i] =
            2.0 * ta1_y_xxy_xxyz_0[i] * fe_0 - 2.0 * ta1_y_xxy_xxyz_1[i] * fe_0 + ta1_y_xxy_xxyzz_0[i] * pa_z[i] - ta1_y_xxy_xxyzz_1[i] * pc_z[i];

        ta1_y_xxyz_xxzzz_0[i] = ta_xxz_xxzzz_1[i] + ta1_y_xxz_xxzzz_0[i] * pa_y[i] - ta1_y_xxz_xxzzz_1[i] * pc_y[i];

        ta1_y_xxyz_xyyyy_0[i] = ta1_y_xxy_xyyyy_0[i] * pa_z[i] - ta1_y_xxy_xyyyy_1[i] * pc_z[i];

        ta1_y_xxyz_xyyyz_0[i] =
            ta1_y_xxy_xyyy_0[i] * fe_0 - ta1_y_xxy_xyyy_1[i] * fe_0 + ta1_y_xxy_xyyyz_0[i] * pa_z[i] - ta1_y_xxy_xyyyz_1[i] * pc_z[i];

        ta1_y_xxyz_xyyzz_0[i] =
            2.0 * ta1_y_xxy_xyyz_0[i] * fe_0 - 2.0 * ta1_y_xxy_xyyz_1[i] * fe_0 + ta1_y_xxy_xyyzz_0[i] * pa_z[i] - ta1_y_xxy_xyyzz_1[i] * pc_z[i];

        ta1_y_xxyz_xyzzz_0[i] =
            3.0 * ta1_y_xxy_xyzz_0[i] * fe_0 - 3.0 * ta1_y_xxy_xyzz_1[i] * fe_0 + ta1_y_xxy_xyzzz_0[i] * pa_z[i] - ta1_y_xxy_xyzzz_1[i] * pc_z[i];

        ta1_y_xxyz_xzzzz_0[i] = ta_xxz_xzzzz_1[i] + ta1_y_xxz_xzzzz_0[i] * pa_y[i] - ta1_y_xxz_xzzzz_1[i] * pc_y[i];

        ta1_y_xxyz_yyyyy_0[i] = ta1_y_xxy_yyyyy_0[i] * pa_z[i] - ta1_y_xxy_yyyyy_1[i] * pc_z[i];

        ta1_y_xxyz_yyyyz_0[i] =
            ta1_y_yz_yyyyz_0[i] * fe_0 - ta1_y_yz_yyyyz_1[i] * fe_0 + ta1_y_xyz_yyyyz_0[i] * pa_x[i] - ta1_y_xyz_yyyyz_1[i] * pc_x[i];

        ta1_y_xxyz_yyyzz_0[i] =
            ta1_y_yz_yyyzz_0[i] * fe_0 - ta1_y_yz_yyyzz_1[i] * fe_0 + ta1_y_xyz_yyyzz_0[i] * pa_x[i] - ta1_y_xyz_yyyzz_1[i] * pc_x[i];

        ta1_y_xxyz_yyzzz_0[i] =
            ta1_y_yz_yyzzz_0[i] * fe_0 - ta1_y_yz_yyzzz_1[i] * fe_0 + ta1_y_xyz_yyzzz_0[i] * pa_x[i] - ta1_y_xyz_yyzzz_1[i] * pc_x[i];

        ta1_y_xxyz_yzzzz_0[i] =
            ta1_y_yz_yzzzz_0[i] * fe_0 - ta1_y_yz_yzzzz_1[i] * fe_0 + ta1_y_xyz_yzzzz_0[i] * pa_x[i] - ta1_y_xyz_yzzzz_1[i] * pc_x[i];

        ta1_y_xxyz_zzzzz_0[i] = ta_xxz_zzzzz_1[i] + ta1_y_xxz_zzzzz_0[i] * pa_y[i] - ta1_y_xxz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 420-441 components of targeted buffer : GH

    auto ta1_y_xxzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 420);

    auto ta1_y_xxzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 421);

    auto ta1_y_xxzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 422);

    auto ta1_y_xxzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 423);

    auto ta1_y_xxzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 424);

    auto ta1_y_xxzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 425);

    auto ta1_y_xxzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 426);

    auto ta1_y_xxzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 427);

    auto ta1_y_xxzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 428);

    auto ta1_y_xxzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 429);

    auto ta1_y_xxzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 430);

    auto ta1_y_xxzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 431);

    auto ta1_y_xxzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 432);

    auto ta1_y_xxzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 433);

    auto ta1_y_xxzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 434);

    auto ta1_y_xxzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 435);

    auto ta1_y_xxzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 436);

    auto ta1_y_xxzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 437);

    auto ta1_y_xxzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 438);

    auto ta1_y_xxzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 439);

    auto ta1_y_xxzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 440);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xx_xxxxx_0,   \
                             ta1_y_xx_xxxxx_1,   \
                             ta1_y_xx_xxxxy_0,   \
                             ta1_y_xx_xxxxy_1,   \
                             ta1_y_xx_xxxyy_0,   \
                             ta1_y_xx_xxxyy_1,   \
                             ta1_y_xx_xxyyy_0,   \
                             ta1_y_xx_xxyyy_1,   \
                             ta1_y_xx_xyyyy_0,   \
                             ta1_y_xx_xyyyy_1,   \
                             ta1_y_xxz_xxxxx_0,  \
                             ta1_y_xxz_xxxxx_1,  \
                             ta1_y_xxz_xxxxy_0,  \
                             ta1_y_xxz_xxxxy_1,  \
                             ta1_y_xxz_xxxyy_0,  \
                             ta1_y_xxz_xxxyy_1,  \
                             ta1_y_xxz_xxyyy_0,  \
                             ta1_y_xxz_xxyyy_1,  \
                             ta1_y_xxz_xyyyy_0,  \
                             ta1_y_xxz_xyyyy_1,  \
                             ta1_y_xxzz_xxxxx_0, \
                             ta1_y_xxzz_xxxxy_0, \
                             ta1_y_xxzz_xxxxz_0, \
                             ta1_y_xxzz_xxxyy_0, \
                             ta1_y_xxzz_xxxyz_0, \
                             ta1_y_xxzz_xxxzz_0, \
                             ta1_y_xxzz_xxyyy_0, \
                             ta1_y_xxzz_xxyyz_0, \
                             ta1_y_xxzz_xxyzz_0, \
                             ta1_y_xxzz_xxzzz_0, \
                             ta1_y_xxzz_xyyyy_0, \
                             ta1_y_xxzz_xyyyz_0, \
                             ta1_y_xxzz_xyyzz_0, \
                             ta1_y_xxzz_xyzzz_0, \
                             ta1_y_xxzz_xzzzz_0, \
                             ta1_y_xxzz_yyyyy_0, \
                             ta1_y_xxzz_yyyyz_0, \
                             ta1_y_xxzz_yyyzz_0, \
                             ta1_y_xxzz_yyzzz_0, \
                             ta1_y_xxzz_yzzzz_0, \
                             ta1_y_xxzz_zzzzz_0, \
                             ta1_y_xzz_xxxxz_0,  \
                             ta1_y_xzz_xxxxz_1,  \
                             ta1_y_xzz_xxxyz_0,  \
                             ta1_y_xzz_xxxyz_1,  \
                             ta1_y_xzz_xxxz_0,   \
                             ta1_y_xzz_xxxz_1,   \
                             ta1_y_xzz_xxxzz_0,  \
                             ta1_y_xzz_xxxzz_1,  \
                             ta1_y_xzz_xxyyz_0,  \
                             ta1_y_xzz_xxyyz_1,  \
                             ta1_y_xzz_xxyz_0,   \
                             ta1_y_xzz_xxyz_1,   \
                             ta1_y_xzz_xxyzz_0,  \
                             ta1_y_xzz_xxyzz_1,  \
                             ta1_y_xzz_xxzz_0,   \
                             ta1_y_xzz_xxzz_1,   \
                             ta1_y_xzz_xxzzz_0,  \
                             ta1_y_xzz_xxzzz_1,  \
                             ta1_y_xzz_xyyyz_0,  \
                             ta1_y_xzz_xyyyz_1,  \
                             ta1_y_xzz_xyyz_0,   \
                             ta1_y_xzz_xyyz_1,   \
                             ta1_y_xzz_xyyzz_0,  \
                             ta1_y_xzz_xyyzz_1,  \
                             ta1_y_xzz_xyzz_0,   \
                             ta1_y_xzz_xyzz_1,   \
                             ta1_y_xzz_xyzzz_0,  \
                             ta1_y_xzz_xyzzz_1,  \
                             ta1_y_xzz_xzzz_0,   \
                             ta1_y_xzz_xzzz_1,   \
                             ta1_y_xzz_xzzzz_0,  \
                             ta1_y_xzz_xzzzz_1,  \
                             ta1_y_xzz_yyyyy_0,  \
                             ta1_y_xzz_yyyyy_1,  \
                             ta1_y_xzz_yyyyz_0,  \
                             ta1_y_xzz_yyyyz_1,  \
                             ta1_y_xzz_yyyz_0,   \
                             ta1_y_xzz_yyyz_1,   \
                             ta1_y_xzz_yyyzz_0,  \
                             ta1_y_xzz_yyyzz_1,  \
                             ta1_y_xzz_yyzz_0,   \
                             ta1_y_xzz_yyzz_1,   \
                             ta1_y_xzz_yyzzz_0,  \
                             ta1_y_xzz_yyzzz_1,  \
                             ta1_y_xzz_yzzz_0,   \
                             ta1_y_xzz_yzzz_1,   \
                             ta1_y_xzz_yzzzz_0,  \
                             ta1_y_xzz_yzzzz_1,  \
                             ta1_y_xzz_zzzz_0,   \
                             ta1_y_xzz_zzzz_1,   \
                             ta1_y_xzz_zzzzz_0,  \
                             ta1_y_xzz_zzzzz_1,  \
                             ta1_y_zz_xxxxz_0,   \
                             ta1_y_zz_xxxxz_1,   \
                             ta1_y_zz_xxxyz_0,   \
                             ta1_y_zz_xxxyz_1,   \
                             ta1_y_zz_xxxzz_0,   \
                             ta1_y_zz_xxxzz_1,   \
                             ta1_y_zz_xxyyz_0,   \
                             ta1_y_zz_xxyyz_1,   \
                             ta1_y_zz_xxyzz_0,   \
                             ta1_y_zz_xxyzz_1,   \
                             ta1_y_zz_xxzzz_0,   \
                             ta1_y_zz_xxzzz_1,   \
                             ta1_y_zz_xyyyz_0,   \
                             ta1_y_zz_xyyyz_1,   \
                             ta1_y_zz_xyyzz_0,   \
                             ta1_y_zz_xyyzz_1,   \
                             ta1_y_zz_xyzzz_0,   \
                             ta1_y_zz_xyzzz_1,   \
                             ta1_y_zz_xzzzz_0,   \
                             ta1_y_zz_xzzzz_1,   \
                             ta1_y_zz_yyyyy_0,   \
                             ta1_y_zz_yyyyy_1,   \
                             ta1_y_zz_yyyyz_0,   \
                             ta1_y_zz_yyyyz_1,   \
                             ta1_y_zz_yyyzz_0,   \
                             ta1_y_zz_yyyzz_1,   \
                             ta1_y_zz_yyzzz_0,   \
                             ta1_y_zz_yyzzz_1,   \
                             ta1_y_zz_yzzzz_0,   \
                             ta1_y_zz_yzzzz_1,   \
                             ta1_y_zz_zzzzz_0,   \
                             ta1_y_zz_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xxzz_xxxxx_0[i] =
            ta1_y_xx_xxxxx_0[i] * fe_0 - ta1_y_xx_xxxxx_1[i] * fe_0 + ta1_y_xxz_xxxxx_0[i] * pa_z[i] - ta1_y_xxz_xxxxx_1[i] * pc_z[i];

        ta1_y_xxzz_xxxxy_0[i] =
            ta1_y_xx_xxxxy_0[i] * fe_0 - ta1_y_xx_xxxxy_1[i] * fe_0 + ta1_y_xxz_xxxxy_0[i] * pa_z[i] - ta1_y_xxz_xxxxy_1[i] * pc_z[i];

        ta1_y_xxzz_xxxxz_0[i] = ta1_y_zz_xxxxz_0[i] * fe_0 - ta1_y_zz_xxxxz_1[i] * fe_0 + 4.0 * ta1_y_xzz_xxxz_0[i] * fe_0 -
                                4.0 * ta1_y_xzz_xxxz_1[i] * fe_0 + ta1_y_xzz_xxxxz_0[i] * pa_x[i] - ta1_y_xzz_xxxxz_1[i] * pc_x[i];

        ta1_y_xxzz_xxxyy_0[i] =
            ta1_y_xx_xxxyy_0[i] * fe_0 - ta1_y_xx_xxxyy_1[i] * fe_0 + ta1_y_xxz_xxxyy_0[i] * pa_z[i] - ta1_y_xxz_xxxyy_1[i] * pc_z[i];

        ta1_y_xxzz_xxxyz_0[i] = ta1_y_zz_xxxyz_0[i] * fe_0 - ta1_y_zz_xxxyz_1[i] * fe_0 + 3.0 * ta1_y_xzz_xxyz_0[i] * fe_0 -
                                3.0 * ta1_y_xzz_xxyz_1[i] * fe_0 + ta1_y_xzz_xxxyz_0[i] * pa_x[i] - ta1_y_xzz_xxxyz_1[i] * pc_x[i];

        ta1_y_xxzz_xxxzz_0[i] = ta1_y_zz_xxxzz_0[i] * fe_0 - ta1_y_zz_xxxzz_1[i] * fe_0 + 3.0 * ta1_y_xzz_xxzz_0[i] * fe_0 -
                                3.0 * ta1_y_xzz_xxzz_1[i] * fe_0 + ta1_y_xzz_xxxzz_0[i] * pa_x[i] - ta1_y_xzz_xxxzz_1[i] * pc_x[i];

        ta1_y_xxzz_xxyyy_0[i] =
            ta1_y_xx_xxyyy_0[i] * fe_0 - ta1_y_xx_xxyyy_1[i] * fe_0 + ta1_y_xxz_xxyyy_0[i] * pa_z[i] - ta1_y_xxz_xxyyy_1[i] * pc_z[i];

        ta1_y_xxzz_xxyyz_0[i] = ta1_y_zz_xxyyz_0[i] * fe_0 - ta1_y_zz_xxyyz_1[i] * fe_0 + 2.0 * ta1_y_xzz_xyyz_0[i] * fe_0 -
                                2.0 * ta1_y_xzz_xyyz_1[i] * fe_0 + ta1_y_xzz_xxyyz_0[i] * pa_x[i] - ta1_y_xzz_xxyyz_1[i] * pc_x[i];

        ta1_y_xxzz_xxyzz_0[i] = ta1_y_zz_xxyzz_0[i] * fe_0 - ta1_y_zz_xxyzz_1[i] * fe_0 + 2.0 * ta1_y_xzz_xyzz_0[i] * fe_0 -
                                2.0 * ta1_y_xzz_xyzz_1[i] * fe_0 + ta1_y_xzz_xxyzz_0[i] * pa_x[i] - ta1_y_xzz_xxyzz_1[i] * pc_x[i];

        ta1_y_xxzz_xxzzz_0[i] = ta1_y_zz_xxzzz_0[i] * fe_0 - ta1_y_zz_xxzzz_1[i] * fe_0 + 2.0 * ta1_y_xzz_xzzz_0[i] * fe_0 -
                                2.0 * ta1_y_xzz_xzzz_1[i] * fe_0 + ta1_y_xzz_xxzzz_0[i] * pa_x[i] - ta1_y_xzz_xxzzz_1[i] * pc_x[i];

        ta1_y_xxzz_xyyyy_0[i] =
            ta1_y_xx_xyyyy_0[i] * fe_0 - ta1_y_xx_xyyyy_1[i] * fe_0 + ta1_y_xxz_xyyyy_0[i] * pa_z[i] - ta1_y_xxz_xyyyy_1[i] * pc_z[i];

        ta1_y_xxzz_xyyyz_0[i] = ta1_y_zz_xyyyz_0[i] * fe_0 - ta1_y_zz_xyyyz_1[i] * fe_0 + ta1_y_xzz_yyyz_0[i] * fe_0 - ta1_y_xzz_yyyz_1[i] * fe_0 +
                                ta1_y_xzz_xyyyz_0[i] * pa_x[i] - ta1_y_xzz_xyyyz_1[i] * pc_x[i];

        ta1_y_xxzz_xyyzz_0[i] = ta1_y_zz_xyyzz_0[i] * fe_0 - ta1_y_zz_xyyzz_1[i] * fe_0 + ta1_y_xzz_yyzz_0[i] * fe_0 - ta1_y_xzz_yyzz_1[i] * fe_0 +
                                ta1_y_xzz_xyyzz_0[i] * pa_x[i] - ta1_y_xzz_xyyzz_1[i] * pc_x[i];

        ta1_y_xxzz_xyzzz_0[i] = ta1_y_zz_xyzzz_0[i] * fe_0 - ta1_y_zz_xyzzz_1[i] * fe_0 + ta1_y_xzz_yzzz_0[i] * fe_0 - ta1_y_xzz_yzzz_1[i] * fe_0 +
                                ta1_y_xzz_xyzzz_0[i] * pa_x[i] - ta1_y_xzz_xyzzz_1[i] * pc_x[i];

        ta1_y_xxzz_xzzzz_0[i] = ta1_y_zz_xzzzz_0[i] * fe_0 - ta1_y_zz_xzzzz_1[i] * fe_0 + ta1_y_xzz_zzzz_0[i] * fe_0 - ta1_y_xzz_zzzz_1[i] * fe_0 +
                                ta1_y_xzz_xzzzz_0[i] * pa_x[i] - ta1_y_xzz_xzzzz_1[i] * pc_x[i];

        ta1_y_xxzz_yyyyy_0[i] =
            ta1_y_zz_yyyyy_0[i] * fe_0 - ta1_y_zz_yyyyy_1[i] * fe_0 + ta1_y_xzz_yyyyy_0[i] * pa_x[i] - ta1_y_xzz_yyyyy_1[i] * pc_x[i];

        ta1_y_xxzz_yyyyz_0[i] =
            ta1_y_zz_yyyyz_0[i] * fe_0 - ta1_y_zz_yyyyz_1[i] * fe_0 + ta1_y_xzz_yyyyz_0[i] * pa_x[i] - ta1_y_xzz_yyyyz_1[i] * pc_x[i];

        ta1_y_xxzz_yyyzz_0[i] =
            ta1_y_zz_yyyzz_0[i] * fe_0 - ta1_y_zz_yyyzz_1[i] * fe_0 + ta1_y_xzz_yyyzz_0[i] * pa_x[i] - ta1_y_xzz_yyyzz_1[i] * pc_x[i];

        ta1_y_xxzz_yyzzz_0[i] =
            ta1_y_zz_yyzzz_0[i] * fe_0 - ta1_y_zz_yyzzz_1[i] * fe_0 + ta1_y_xzz_yyzzz_0[i] * pa_x[i] - ta1_y_xzz_yyzzz_1[i] * pc_x[i];

        ta1_y_xxzz_yzzzz_0[i] =
            ta1_y_zz_yzzzz_0[i] * fe_0 - ta1_y_zz_yzzzz_1[i] * fe_0 + ta1_y_xzz_yzzzz_0[i] * pa_x[i] - ta1_y_xzz_yzzzz_1[i] * pc_x[i];

        ta1_y_xxzz_zzzzz_0[i] =
            ta1_y_zz_zzzzz_0[i] * fe_0 - ta1_y_zz_zzzzz_1[i] * fe_0 + ta1_y_xzz_zzzzz_0[i] * pa_x[i] - ta1_y_xzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 441-462 components of targeted buffer : GH

    auto ta1_y_xyyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 441);

    auto ta1_y_xyyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 442);

    auto ta1_y_xyyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 443);

    auto ta1_y_xyyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 444);

    auto ta1_y_xyyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 445);

    auto ta1_y_xyyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 446);

    auto ta1_y_xyyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 447);

    auto ta1_y_xyyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 448);

    auto ta1_y_xyyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 449);

    auto ta1_y_xyyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 450);

    auto ta1_y_xyyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 451);

    auto ta1_y_xyyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 452);

    auto ta1_y_xyyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 453);

    auto ta1_y_xyyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 454);

    auto ta1_y_xyyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 455);

    auto ta1_y_xyyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 456);

    auto ta1_y_xyyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 457);

    auto ta1_y_xyyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 458);

    auto ta1_y_xyyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 459);

    auto ta1_y_xyyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 460);

    auto ta1_y_xyyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 461);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_xyyy_xxxxx_0, \
                             ta1_y_xyyy_xxxxy_0, \
                             ta1_y_xyyy_xxxxz_0, \
                             ta1_y_xyyy_xxxyy_0, \
                             ta1_y_xyyy_xxxyz_0, \
                             ta1_y_xyyy_xxxzz_0, \
                             ta1_y_xyyy_xxyyy_0, \
                             ta1_y_xyyy_xxyyz_0, \
                             ta1_y_xyyy_xxyzz_0, \
                             ta1_y_xyyy_xxzzz_0, \
                             ta1_y_xyyy_xyyyy_0, \
                             ta1_y_xyyy_xyyyz_0, \
                             ta1_y_xyyy_xyyzz_0, \
                             ta1_y_xyyy_xyzzz_0, \
                             ta1_y_xyyy_xzzzz_0, \
                             ta1_y_xyyy_yyyyy_0, \
                             ta1_y_xyyy_yyyyz_0, \
                             ta1_y_xyyy_yyyzz_0, \
                             ta1_y_xyyy_yyzzz_0, \
                             ta1_y_xyyy_yzzzz_0, \
                             ta1_y_xyyy_zzzzz_0, \
                             ta1_y_yyy_xxxx_0,   \
                             ta1_y_yyy_xxxx_1,   \
                             ta1_y_yyy_xxxxx_0,  \
                             ta1_y_yyy_xxxxx_1,  \
                             ta1_y_yyy_xxxxy_0,  \
                             ta1_y_yyy_xxxxy_1,  \
                             ta1_y_yyy_xxxxz_0,  \
                             ta1_y_yyy_xxxxz_1,  \
                             ta1_y_yyy_xxxy_0,   \
                             ta1_y_yyy_xxxy_1,   \
                             ta1_y_yyy_xxxyy_0,  \
                             ta1_y_yyy_xxxyy_1,  \
                             ta1_y_yyy_xxxyz_0,  \
                             ta1_y_yyy_xxxyz_1,  \
                             ta1_y_yyy_xxxz_0,   \
                             ta1_y_yyy_xxxz_1,   \
                             ta1_y_yyy_xxxzz_0,  \
                             ta1_y_yyy_xxxzz_1,  \
                             ta1_y_yyy_xxyy_0,   \
                             ta1_y_yyy_xxyy_1,   \
                             ta1_y_yyy_xxyyy_0,  \
                             ta1_y_yyy_xxyyy_1,  \
                             ta1_y_yyy_xxyyz_0,  \
                             ta1_y_yyy_xxyyz_1,  \
                             ta1_y_yyy_xxyz_0,   \
                             ta1_y_yyy_xxyz_1,   \
                             ta1_y_yyy_xxyzz_0,  \
                             ta1_y_yyy_xxyzz_1,  \
                             ta1_y_yyy_xxzz_0,   \
                             ta1_y_yyy_xxzz_1,   \
                             ta1_y_yyy_xxzzz_0,  \
                             ta1_y_yyy_xxzzz_1,  \
                             ta1_y_yyy_xyyy_0,   \
                             ta1_y_yyy_xyyy_1,   \
                             ta1_y_yyy_xyyyy_0,  \
                             ta1_y_yyy_xyyyy_1,  \
                             ta1_y_yyy_xyyyz_0,  \
                             ta1_y_yyy_xyyyz_1,  \
                             ta1_y_yyy_xyyz_0,   \
                             ta1_y_yyy_xyyz_1,   \
                             ta1_y_yyy_xyyzz_0,  \
                             ta1_y_yyy_xyyzz_1,  \
                             ta1_y_yyy_xyzz_0,   \
                             ta1_y_yyy_xyzz_1,   \
                             ta1_y_yyy_xyzzz_0,  \
                             ta1_y_yyy_xyzzz_1,  \
                             ta1_y_yyy_xzzz_0,   \
                             ta1_y_yyy_xzzz_1,   \
                             ta1_y_yyy_xzzzz_0,  \
                             ta1_y_yyy_xzzzz_1,  \
                             ta1_y_yyy_yyyy_0,   \
                             ta1_y_yyy_yyyy_1,   \
                             ta1_y_yyy_yyyyy_0,  \
                             ta1_y_yyy_yyyyy_1,  \
                             ta1_y_yyy_yyyyz_0,  \
                             ta1_y_yyy_yyyyz_1,  \
                             ta1_y_yyy_yyyz_0,   \
                             ta1_y_yyy_yyyz_1,   \
                             ta1_y_yyy_yyyzz_0,  \
                             ta1_y_yyy_yyyzz_1,  \
                             ta1_y_yyy_yyzz_0,   \
                             ta1_y_yyy_yyzz_1,   \
                             ta1_y_yyy_yyzzz_0,  \
                             ta1_y_yyy_yyzzz_1,  \
                             ta1_y_yyy_yzzz_0,   \
                             ta1_y_yyy_yzzz_1,   \
                             ta1_y_yyy_yzzzz_0,  \
                             ta1_y_yyy_yzzzz_1,  \
                             ta1_y_yyy_zzzz_0,   \
                             ta1_y_yyy_zzzz_1,   \
                             ta1_y_yyy_zzzzz_0,  \
                             ta1_y_yyy_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyy_xxxxx_0[i] =
            5.0 * ta1_y_yyy_xxxx_0[i] * fe_0 - 5.0 * ta1_y_yyy_xxxx_1[i] * fe_0 + ta1_y_yyy_xxxxx_0[i] * pa_x[i] - ta1_y_yyy_xxxxx_1[i] * pc_x[i];

        ta1_y_xyyy_xxxxy_0[i] =
            4.0 * ta1_y_yyy_xxxy_0[i] * fe_0 - 4.0 * ta1_y_yyy_xxxy_1[i] * fe_0 + ta1_y_yyy_xxxxy_0[i] * pa_x[i] - ta1_y_yyy_xxxxy_1[i] * pc_x[i];

        ta1_y_xyyy_xxxxz_0[i] =
            4.0 * ta1_y_yyy_xxxz_0[i] * fe_0 - 4.0 * ta1_y_yyy_xxxz_1[i] * fe_0 + ta1_y_yyy_xxxxz_0[i] * pa_x[i] - ta1_y_yyy_xxxxz_1[i] * pc_x[i];

        ta1_y_xyyy_xxxyy_0[i] =
            3.0 * ta1_y_yyy_xxyy_0[i] * fe_0 - 3.0 * ta1_y_yyy_xxyy_1[i] * fe_0 + ta1_y_yyy_xxxyy_0[i] * pa_x[i] - ta1_y_yyy_xxxyy_1[i] * pc_x[i];

        ta1_y_xyyy_xxxyz_0[i] =
            3.0 * ta1_y_yyy_xxyz_0[i] * fe_0 - 3.0 * ta1_y_yyy_xxyz_1[i] * fe_0 + ta1_y_yyy_xxxyz_0[i] * pa_x[i] - ta1_y_yyy_xxxyz_1[i] * pc_x[i];

        ta1_y_xyyy_xxxzz_0[i] =
            3.0 * ta1_y_yyy_xxzz_0[i] * fe_0 - 3.0 * ta1_y_yyy_xxzz_1[i] * fe_0 + ta1_y_yyy_xxxzz_0[i] * pa_x[i] - ta1_y_yyy_xxxzz_1[i] * pc_x[i];

        ta1_y_xyyy_xxyyy_0[i] =
            2.0 * ta1_y_yyy_xyyy_0[i] * fe_0 - 2.0 * ta1_y_yyy_xyyy_1[i] * fe_0 + ta1_y_yyy_xxyyy_0[i] * pa_x[i] - ta1_y_yyy_xxyyy_1[i] * pc_x[i];

        ta1_y_xyyy_xxyyz_0[i] =
            2.0 * ta1_y_yyy_xyyz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xyyz_1[i] * fe_0 + ta1_y_yyy_xxyyz_0[i] * pa_x[i] - ta1_y_yyy_xxyyz_1[i] * pc_x[i];

        ta1_y_xyyy_xxyzz_0[i] =
            2.0 * ta1_y_yyy_xyzz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xyzz_1[i] * fe_0 + ta1_y_yyy_xxyzz_0[i] * pa_x[i] - ta1_y_yyy_xxyzz_1[i] * pc_x[i];

        ta1_y_xyyy_xxzzz_0[i] =
            2.0 * ta1_y_yyy_xzzz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xzzz_1[i] * fe_0 + ta1_y_yyy_xxzzz_0[i] * pa_x[i] - ta1_y_yyy_xxzzz_1[i] * pc_x[i];

        ta1_y_xyyy_xyyyy_0[i] =
            ta1_y_yyy_yyyy_0[i] * fe_0 - ta1_y_yyy_yyyy_1[i] * fe_0 + ta1_y_yyy_xyyyy_0[i] * pa_x[i] - ta1_y_yyy_xyyyy_1[i] * pc_x[i];

        ta1_y_xyyy_xyyyz_0[i] =
            ta1_y_yyy_yyyz_0[i] * fe_0 - ta1_y_yyy_yyyz_1[i] * fe_0 + ta1_y_yyy_xyyyz_0[i] * pa_x[i] - ta1_y_yyy_xyyyz_1[i] * pc_x[i];

        ta1_y_xyyy_xyyzz_0[i] =
            ta1_y_yyy_yyzz_0[i] * fe_0 - ta1_y_yyy_yyzz_1[i] * fe_0 + ta1_y_yyy_xyyzz_0[i] * pa_x[i] - ta1_y_yyy_xyyzz_1[i] * pc_x[i];

        ta1_y_xyyy_xyzzz_0[i] =
            ta1_y_yyy_yzzz_0[i] * fe_0 - ta1_y_yyy_yzzz_1[i] * fe_0 + ta1_y_yyy_xyzzz_0[i] * pa_x[i] - ta1_y_yyy_xyzzz_1[i] * pc_x[i];

        ta1_y_xyyy_xzzzz_0[i] =
            ta1_y_yyy_zzzz_0[i] * fe_0 - ta1_y_yyy_zzzz_1[i] * fe_0 + ta1_y_yyy_xzzzz_0[i] * pa_x[i] - ta1_y_yyy_xzzzz_1[i] * pc_x[i];

        ta1_y_xyyy_yyyyy_0[i] = ta1_y_yyy_yyyyy_0[i] * pa_x[i] - ta1_y_yyy_yyyyy_1[i] * pc_x[i];

        ta1_y_xyyy_yyyyz_0[i] = ta1_y_yyy_yyyyz_0[i] * pa_x[i] - ta1_y_yyy_yyyyz_1[i] * pc_x[i];

        ta1_y_xyyy_yyyzz_0[i] = ta1_y_yyy_yyyzz_0[i] * pa_x[i] - ta1_y_yyy_yyyzz_1[i] * pc_x[i];

        ta1_y_xyyy_yyzzz_0[i] = ta1_y_yyy_yyzzz_0[i] * pa_x[i] - ta1_y_yyy_yyzzz_1[i] * pc_x[i];

        ta1_y_xyyy_yzzzz_0[i] = ta1_y_yyy_yzzzz_0[i] * pa_x[i] - ta1_y_yyy_yzzzz_1[i] * pc_x[i];

        ta1_y_xyyy_zzzzz_0[i] = ta1_y_yyy_zzzzz_0[i] * pa_x[i] - ta1_y_yyy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 462-483 components of targeted buffer : GH

    auto ta1_y_xyyz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 462);

    auto ta1_y_xyyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 463);

    auto ta1_y_xyyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 464);

    auto ta1_y_xyyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 465);

    auto ta1_y_xyyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 466);

    auto ta1_y_xyyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 467);

    auto ta1_y_xyyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 468);

    auto ta1_y_xyyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 469);

    auto ta1_y_xyyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 470);

    auto ta1_y_xyyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 471);

    auto ta1_y_xyyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 472);

    auto ta1_y_xyyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 473);

    auto ta1_y_xyyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 474);

    auto ta1_y_xyyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 475);

    auto ta1_y_xyyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 476);

    auto ta1_y_xyyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 477);

    auto ta1_y_xyyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 478);

    auto ta1_y_xyyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 479);

    auto ta1_y_xyyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 480);

    auto ta1_y_xyyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 481);

    auto ta1_y_xyyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 482);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_y_xyy_xxxxx_0,  \
                             ta1_y_xyy_xxxxx_1,  \
                             ta1_y_xyy_xxxxy_0,  \
                             ta1_y_xyy_xxxxy_1,  \
                             ta1_y_xyy_xxxyy_0,  \
                             ta1_y_xyy_xxxyy_1,  \
                             ta1_y_xyy_xxyyy_0,  \
                             ta1_y_xyy_xxyyy_1,  \
                             ta1_y_xyy_xyyyy_0,  \
                             ta1_y_xyy_xyyyy_1,  \
                             ta1_y_xyyz_xxxxx_0, \
                             ta1_y_xyyz_xxxxy_0, \
                             ta1_y_xyyz_xxxxz_0, \
                             ta1_y_xyyz_xxxyy_0, \
                             ta1_y_xyyz_xxxyz_0, \
                             ta1_y_xyyz_xxxzz_0, \
                             ta1_y_xyyz_xxyyy_0, \
                             ta1_y_xyyz_xxyyz_0, \
                             ta1_y_xyyz_xxyzz_0, \
                             ta1_y_xyyz_xxzzz_0, \
                             ta1_y_xyyz_xyyyy_0, \
                             ta1_y_xyyz_xyyyz_0, \
                             ta1_y_xyyz_xyyzz_0, \
                             ta1_y_xyyz_xyzzz_0, \
                             ta1_y_xyyz_xzzzz_0, \
                             ta1_y_xyyz_yyyyy_0, \
                             ta1_y_xyyz_yyyyz_0, \
                             ta1_y_xyyz_yyyzz_0, \
                             ta1_y_xyyz_yyzzz_0, \
                             ta1_y_xyyz_yzzzz_0, \
                             ta1_y_xyyz_zzzzz_0, \
                             ta1_y_yyz_xxxxz_0,  \
                             ta1_y_yyz_xxxxz_1,  \
                             ta1_y_yyz_xxxyz_0,  \
                             ta1_y_yyz_xxxyz_1,  \
                             ta1_y_yyz_xxxz_0,   \
                             ta1_y_yyz_xxxz_1,   \
                             ta1_y_yyz_xxxzz_0,  \
                             ta1_y_yyz_xxxzz_1,  \
                             ta1_y_yyz_xxyyz_0,  \
                             ta1_y_yyz_xxyyz_1,  \
                             ta1_y_yyz_xxyz_0,   \
                             ta1_y_yyz_xxyz_1,   \
                             ta1_y_yyz_xxyzz_0,  \
                             ta1_y_yyz_xxyzz_1,  \
                             ta1_y_yyz_xxzz_0,   \
                             ta1_y_yyz_xxzz_1,   \
                             ta1_y_yyz_xxzzz_0,  \
                             ta1_y_yyz_xxzzz_1,  \
                             ta1_y_yyz_xyyyz_0,  \
                             ta1_y_yyz_xyyyz_1,  \
                             ta1_y_yyz_xyyz_0,   \
                             ta1_y_yyz_xyyz_1,   \
                             ta1_y_yyz_xyyzz_0,  \
                             ta1_y_yyz_xyyzz_1,  \
                             ta1_y_yyz_xyzz_0,   \
                             ta1_y_yyz_xyzz_1,   \
                             ta1_y_yyz_xyzzz_0,  \
                             ta1_y_yyz_xyzzz_1,  \
                             ta1_y_yyz_xzzz_0,   \
                             ta1_y_yyz_xzzz_1,   \
                             ta1_y_yyz_xzzzz_0,  \
                             ta1_y_yyz_xzzzz_1,  \
                             ta1_y_yyz_yyyyy_0,  \
                             ta1_y_yyz_yyyyy_1,  \
                             ta1_y_yyz_yyyyz_0,  \
                             ta1_y_yyz_yyyyz_1,  \
                             ta1_y_yyz_yyyz_0,   \
                             ta1_y_yyz_yyyz_1,   \
                             ta1_y_yyz_yyyzz_0,  \
                             ta1_y_yyz_yyyzz_1,  \
                             ta1_y_yyz_yyzz_0,   \
                             ta1_y_yyz_yyzz_1,   \
                             ta1_y_yyz_yyzzz_0,  \
                             ta1_y_yyz_yyzzz_1,  \
                             ta1_y_yyz_yzzz_0,   \
                             ta1_y_yyz_yzzz_1,   \
                             ta1_y_yyz_yzzzz_0,  \
                             ta1_y_yyz_yzzzz_1,  \
                             ta1_y_yyz_zzzz_0,   \
                             ta1_y_yyz_zzzz_1,   \
                             ta1_y_yyz_zzzzz_0,  \
                             ta1_y_yyz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyyz_xxxxx_0[i] = ta1_y_xyy_xxxxx_0[i] * pa_z[i] - ta1_y_xyy_xxxxx_1[i] * pc_z[i];

        ta1_y_xyyz_xxxxy_0[i] = ta1_y_xyy_xxxxy_0[i] * pa_z[i] - ta1_y_xyy_xxxxy_1[i] * pc_z[i];

        ta1_y_xyyz_xxxxz_0[i] =
            4.0 * ta1_y_yyz_xxxz_0[i] * fe_0 - 4.0 * ta1_y_yyz_xxxz_1[i] * fe_0 + ta1_y_yyz_xxxxz_0[i] * pa_x[i] - ta1_y_yyz_xxxxz_1[i] * pc_x[i];

        ta1_y_xyyz_xxxyy_0[i] = ta1_y_xyy_xxxyy_0[i] * pa_z[i] - ta1_y_xyy_xxxyy_1[i] * pc_z[i];

        ta1_y_xyyz_xxxyz_0[i] =
            3.0 * ta1_y_yyz_xxyz_0[i] * fe_0 - 3.0 * ta1_y_yyz_xxyz_1[i] * fe_0 + ta1_y_yyz_xxxyz_0[i] * pa_x[i] - ta1_y_yyz_xxxyz_1[i] * pc_x[i];

        ta1_y_xyyz_xxxzz_0[i] =
            3.0 * ta1_y_yyz_xxzz_0[i] * fe_0 - 3.0 * ta1_y_yyz_xxzz_1[i] * fe_0 + ta1_y_yyz_xxxzz_0[i] * pa_x[i] - ta1_y_yyz_xxxzz_1[i] * pc_x[i];

        ta1_y_xyyz_xxyyy_0[i] = ta1_y_xyy_xxyyy_0[i] * pa_z[i] - ta1_y_xyy_xxyyy_1[i] * pc_z[i];

        ta1_y_xyyz_xxyyz_0[i] =
            2.0 * ta1_y_yyz_xyyz_0[i] * fe_0 - 2.0 * ta1_y_yyz_xyyz_1[i] * fe_0 + ta1_y_yyz_xxyyz_0[i] * pa_x[i] - ta1_y_yyz_xxyyz_1[i] * pc_x[i];

        ta1_y_xyyz_xxyzz_0[i] =
            2.0 * ta1_y_yyz_xyzz_0[i] * fe_0 - 2.0 * ta1_y_yyz_xyzz_1[i] * fe_0 + ta1_y_yyz_xxyzz_0[i] * pa_x[i] - ta1_y_yyz_xxyzz_1[i] * pc_x[i];

        ta1_y_xyyz_xxzzz_0[i] =
            2.0 * ta1_y_yyz_xzzz_0[i] * fe_0 - 2.0 * ta1_y_yyz_xzzz_1[i] * fe_0 + ta1_y_yyz_xxzzz_0[i] * pa_x[i] - ta1_y_yyz_xxzzz_1[i] * pc_x[i];

        ta1_y_xyyz_xyyyy_0[i] = ta1_y_xyy_xyyyy_0[i] * pa_z[i] - ta1_y_xyy_xyyyy_1[i] * pc_z[i];

        ta1_y_xyyz_xyyyz_0[i] =
            ta1_y_yyz_yyyz_0[i] * fe_0 - ta1_y_yyz_yyyz_1[i] * fe_0 + ta1_y_yyz_xyyyz_0[i] * pa_x[i] - ta1_y_yyz_xyyyz_1[i] * pc_x[i];

        ta1_y_xyyz_xyyzz_0[i] =
            ta1_y_yyz_yyzz_0[i] * fe_0 - ta1_y_yyz_yyzz_1[i] * fe_0 + ta1_y_yyz_xyyzz_0[i] * pa_x[i] - ta1_y_yyz_xyyzz_1[i] * pc_x[i];

        ta1_y_xyyz_xyzzz_0[i] =
            ta1_y_yyz_yzzz_0[i] * fe_0 - ta1_y_yyz_yzzz_1[i] * fe_0 + ta1_y_yyz_xyzzz_0[i] * pa_x[i] - ta1_y_yyz_xyzzz_1[i] * pc_x[i];

        ta1_y_xyyz_xzzzz_0[i] =
            ta1_y_yyz_zzzz_0[i] * fe_0 - ta1_y_yyz_zzzz_1[i] * fe_0 + ta1_y_yyz_xzzzz_0[i] * pa_x[i] - ta1_y_yyz_xzzzz_1[i] * pc_x[i];

        ta1_y_xyyz_yyyyy_0[i] = ta1_y_yyz_yyyyy_0[i] * pa_x[i] - ta1_y_yyz_yyyyy_1[i] * pc_x[i];

        ta1_y_xyyz_yyyyz_0[i] = ta1_y_yyz_yyyyz_0[i] * pa_x[i] - ta1_y_yyz_yyyyz_1[i] * pc_x[i];

        ta1_y_xyyz_yyyzz_0[i] = ta1_y_yyz_yyyzz_0[i] * pa_x[i] - ta1_y_yyz_yyyzz_1[i] * pc_x[i];

        ta1_y_xyyz_yyzzz_0[i] = ta1_y_yyz_yyzzz_0[i] * pa_x[i] - ta1_y_yyz_yyzzz_1[i] * pc_x[i];

        ta1_y_xyyz_yzzzz_0[i] = ta1_y_yyz_yzzzz_0[i] * pa_x[i] - ta1_y_yyz_yzzzz_1[i] * pc_x[i];

        ta1_y_xyyz_zzzzz_0[i] = ta1_y_yyz_zzzzz_0[i] * pa_x[i] - ta1_y_yyz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 483-504 components of targeted buffer : GH

    auto ta1_y_xyzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 483);

    auto ta1_y_xyzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 484);

    auto ta1_y_xyzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 485);

    auto ta1_y_xyzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 486);

    auto ta1_y_xyzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 487);

    auto ta1_y_xyzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 488);

    auto ta1_y_xyzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 489);

    auto ta1_y_xyzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 490);

    auto ta1_y_xyzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 491);

    auto ta1_y_xyzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 492);

    auto ta1_y_xyzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 493);

    auto ta1_y_xyzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 494);

    auto ta1_y_xyzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 495);

    auto ta1_y_xyzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 496);

    auto ta1_y_xyzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 497);

    auto ta1_y_xyzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 498);

    auto ta1_y_xyzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 499);

    auto ta1_y_xyzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 500);

    auto ta1_y_xyzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 501);

    auto ta1_y_xyzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 502);

    auto ta1_y_xyzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 503);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_y_xyzz_xxxxx_0, \
                             ta1_y_xyzz_xxxxy_0, \
                             ta1_y_xyzz_xxxxz_0, \
                             ta1_y_xyzz_xxxyy_0, \
                             ta1_y_xyzz_xxxyz_0, \
                             ta1_y_xyzz_xxxzz_0, \
                             ta1_y_xyzz_xxyyy_0, \
                             ta1_y_xyzz_xxyyz_0, \
                             ta1_y_xyzz_xxyzz_0, \
                             ta1_y_xyzz_xxzzz_0, \
                             ta1_y_xyzz_xyyyy_0, \
                             ta1_y_xyzz_xyyyz_0, \
                             ta1_y_xyzz_xyyzz_0, \
                             ta1_y_xyzz_xyzzz_0, \
                             ta1_y_xyzz_xzzzz_0, \
                             ta1_y_xyzz_yyyyy_0, \
                             ta1_y_xyzz_yyyyz_0, \
                             ta1_y_xyzz_yyyzz_0, \
                             ta1_y_xyzz_yyzzz_0, \
                             ta1_y_xyzz_yzzzz_0, \
                             ta1_y_xyzz_zzzzz_0, \
                             ta1_y_xzz_xxxxx_0,  \
                             ta1_y_xzz_xxxxx_1,  \
                             ta1_y_xzz_xxxxz_0,  \
                             ta1_y_xzz_xxxxz_1,  \
                             ta1_y_xzz_xxxzz_0,  \
                             ta1_y_xzz_xxxzz_1,  \
                             ta1_y_xzz_xxzzz_0,  \
                             ta1_y_xzz_xxzzz_1,  \
                             ta1_y_xzz_xzzzz_0,  \
                             ta1_y_xzz_xzzzz_1,  \
                             ta1_y_yzz_xxxxy_0,  \
                             ta1_y_yzz_xxxxy_1,  \
                             ta1_y_yzz_xxxy_0,   \
                             ta1_y_yzz_xxxy_1,   \
                             ta1_y_yzz_xxxyy_0,  \
                             ta1_y_yzz_xxxyy_1,  \
                             ta1_y_yzz_xxxyz_0,  \
                             ta1_y_yzz_xxxyz_1,  \
                             ta1_y_yzz_xxyy_0,   \
                             ta1_y_yzz_xxyy_1,   \
                             ta1_y_yzz_xxyyy_0,  \
                             ta1_y_yzz_xxyyy_1,  \
                             ta1_y_yzz_xxyyz_0,  \
                             ta1_y_yzz_xxyyz_1,  \
                             ta1_y_yzz_xxyz_0,   \
                             ta1_y_yzz_xxyz_1,   \
                             ta1_y_yzz_xxyzz_0,  \
                             ta1_y_yzz_xxyzz_1,  \
                             ta1_y_yzz_xyyy_0,   \
                             ta1_y_yzz_xyyy_1,   \
                             ta1_y_yzz_xyyyy_0,  \
                             ta1_y_yzz_xyyyy_1,  \
                             ta1_y_yzz_xyyyz_0,  \
                             ta1_y_yzz_xyyyz_1,  \
                             ta1_y_yzz_xyyz_0,   \
                             ta1_y_yzz_xyyz_1,   \
                             ta1_y_yzz_xyyzz_0,  \
                             ta1_y_yzz_xyyzz_1,  \
                             ta1_y_yzz_xyzz_0,   \
                             ta1_y_yzz_xyzz_1,   \
                             ta1_y_yzz_xyzzz_0,  \
                             ta1_y_yzz_xyzzz_1,  \
                             ta1_y_yzz_yyyy_0,   \
                             ta1_y_yzz_yyyy_1,   \
                             ta1_y_yzz_yyyyy_0,  \
                             ta1_y_yzz_yyyyy_1,  \
                             ta1_y_yzz_yyyyz_0,  \
                             ta1_y_yzz_yyyyz_1,  \
                             ta1_y_yzz_yyyz_0,   \
                             ta1_y_yzz_yyyz_1,   \
                             ta1_y_yzz_yyyzz_0,  \
                             ta1_y_yzz_yyyzz_1,  \
                             ta1_y_yzz_yyzz_0,   \
                             ta1_y_yzz_yyzz_1,   \
                             ta1_y_yzz_yyzzz_0,  \
                             ta1_y_yzz_yyzzz_1,  \
                             ta1_y_yzz_yzzz_0,   \
                             ta1_y_yzz_yzzz_1,   \
                             ta1_y_yzz_yzzzz_0,  \
                             ta1_y_yzz_yzzzz_1,  \
                             ta1_y_yzz_zzzzz_0,  \
                             ta1_y_yzz_zzzzz_1,  \
                             ta_xzz_xxxxx_1,     \
                             ta_xzz_xxxxz_1,     \
                             ta_xzz_xxxzz_1,     \
                             ta_xzz_xxzzz_1,     \
                             ta_xzz_xzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xyzz_xxxxx_0[i] = ta_xzz_xxxxx_1[i] + ta1_y_xzz_xxxxx_0[i] * pa_y[i] - ta1_y_xzz_xxxxx_1[i] * pc_y[i];

        ta1_y_xyzz_xxxxy_0[i] =
            4.0 * ta1_y_yzz_xxxy_0[i] * fe_0 - 4.0 * ta1_y_yzz_xxxy_1[i] * fe_0 + ta1_y_yzz_xxxxy_0[i] * pa_x[i] - ta1_y_yzz_xxxxy_1[i] * pc_x[i];

        ta1_y_xyzz_xxxxz_0[i] = ta_xzz_xxxxz_1[i] + ta1_y_xzz_xxxxz_0[i] * pa_y[i] - ta1_y_xzz_xxxxz_1[i] * pc_y[i];

        ta1_y_xyzz_xxxyy_0[i] =
            3.0 * ta1_y_yzz_xxyy_0[i] * fe_0 - 3.0 * ta1_y_yzz_xxyy_1[i] * fe_0 + ta1_y_yzz_xxxyy_0[i] * pa_x[i] - ta1_y_yzz_xxxyy_1[i] * pc_x[i];

        ta1_y_xyzz_xxxyz_0[i] =
            3.0 * ta1_y_yzz_xxyz_0[i] * fe_0 - 3.0 * ta1_y_yzz_xxyz_1[i] * fe_0 + ta1_y_yzz_xxxyz_0[i] * pa_x[i] - ta1_y_yzz_xxxyz_1[i] * pc_x[i];

        ta1_y_xyzz_xxxzz_0[i] = ta_xzz_xxxzz_1[i] + ta1_y_xzz_xxxzz_0[i] * pa_y[i] - ta1_y_xzz_xxxzz_1[i] * pc_y[i];

        ta1_y_xyzz_xxyyy_0[i] =
            2.0 * ta1_y_yzz_xyyy_0[i] * fe_0 - 2.0 * ta1_y_yzz_xyyy_1[i] * fe_0 + ta1_y_yzz_xxyyy_0[i] * pa_x[i] - ta1_y_yzz_xxyyy_1[i] * pc_x[i];

        ta1_y_xyzz_xxyyz_0[i] =
            2.0 * ta1_y_yzz_xyyz_0[i] * fe_0 - 2.0 * ta1_y_yzz_xyyz_1[i] * fe_0 + ta1_y_yzz_xxyyz_0[i] * pa_x[i] - ta1_y_yzz_xxyyz_1[i] * pc_x[i];

        ta1_y_xyzz_xxyzz_0[i] =
            2.0 * ta1_y_yzz_xyzz_0[i] * fe_0 - 2.0 * ta1_y_yzz_xyzz_1[i] * fe_0 + ta1_y_yzz_xxyzz_0[i] * pa_x[i] - ta1_y_yzz_xxyzz_1[i] * pc_x[i];

        ta1_y_xyzz_xxzzz_0[i] = ta_xzz_xxzzz_1[i] + ta1_y_xzz_xxzzz_0[i] * pa_y[i] - ta1_y_xzz_xxzzz_1[i] * pc_y[i];

        ta1_y_xyzz_xyyyy_0[i] =
            ta1_y_yzz_yyyy_0[i] * fe_0 - ta1_y_yzz_yyyy_1[i] * fe_0 + ta1_y_yzz_xyyyy_0[i] * pa_x[i] - ta1_y_yzz_xyyyy_1[i] * pc_x[i];

        ta1_y_xyzz_xyyyz_0[i] =
            ta1_y_yzz_yyyz_0[i] * fe_0 - ta1_y_yzz_yyyz_1[i] * fe_0 + ta1_y_yzz_xyyyz_0[i] * pa_x[i] - ta1_y_yzz_xyyyz_1[i] * pc_x[i];

        ta1_y_xyzz_xyyzz_0[i] =
            ta1_y_yzz_yyzz_0[i] * fe_0 - ta1_y_yzz_yyzz_1[i] * fe_0 + ta1_y_yzz_xyyzz_0[i] * pa_x[i] - ta1_y_yzz_xyyzz_1[i] * pc_x[i];

        ta1_y_xyzz_xyzzz_0[i] =
            ta1_y_yzz_yzzz_0[i] * fe_0 - ta1_y_yzz_yzzz_1[i] * fe_0 + ta1_y_yzz_xyzzz_0[i] * pa_x[i] - ta1_y_yzz_xyzzz_1[i] * pc_x[i];

        ta1_y_xyzz_xzzzz_0[i] = ta_xzz_xzzzz_1[i] + ta1_y_xzz_xzzzz_0[i] * pa_y[i] - ta1_y_xzz_xzzzz_1[i] * pc_y[i];

        ta1_y_xyzz_yyyyy_0[i] = ta1_y_yzz_yyyyy_0[i] * pa_x[i] - ta1_y_yzz_yyyyy_1[i] * pc_x[i];

        ta1_y_xyzz_yyyyz_0[i] = ta1_y_yzz_yyyyz_0[i] * pa_x[i] - ta1_y_yzz_yyyyz_1[i] * pc_x[i];

        ta1_y_xyzz_yyyzz_0[i] = ta1_y_yzz_yyyzz_0[i] * pa_x[i] - ta1_y_yzz_yyyzz_1[i] * pc_x[i];

        ta1_y_xyzz_yyzzz_0[i] = ta1_y_yzz_yyzzz_0[i] * pa_x[i] - ta1_y_yzz_yyzzz_1[i] * pc_x[i];

        ta1_y_xyzz_yzzzz_0[i] = ta1_y_yzz_yzzzz_0[i] * pa_x[i] - ta1_y_yzz_yzzzz_1[i] * pc_x[i];

        ta1_y_xyzz_zzzzz_0[i] = ta1_y_yzz_zzzzz_0[i] * pa_x[i] - ta1_y_yzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 504-525 components of targeted buffer : GH

    auto ta1_y_xzzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 504);

    auto ta1_y_xzzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 505);

    auto ta1_y_xzzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 506);

    auto ta1_y_xzzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 507);

    auto ta1_y_xzzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 508);

    auto ta1_y_xzzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 509);

    auto ta1_y_xzzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 510);

    auto ta1_y_xzzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 511);

    auto ta1_y_xzzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 512);

    auto ta1_y_xzzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 513);

    auto ta1_y_xzzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 514);

    auto ta1_y_xzzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 515);

    auto ta1_y_xzzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 516);

    auto ta1_y_xzzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 517);

    auto ta1_y_xzzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 518);

    auto ta1_y_xzzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 519);

    auto ta1_y_xzzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 520);

    auto ta1_y_xzzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 521);

    auto ta1_y_xzzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 522);

    auto ta1_y_xzzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 523);

    auto ta1_y_xzzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 524);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_y_xzzz_xxxxx_0, \
                             ta1_y_xzzz_xxxxy_0, \
                             ta1_y_xzzz_xxxxz_0, \
                             ta1_y_xzzz_xxxyy_0, \
                             ta1_y_xzzz_xxxyz_0, \
                             ta1_y_xzzz_xxxzz_0, \
                             ta1_y_xzzz_xxyyy_0, \
                             ta1_y_xzzz_xxyyz_0, \
                             ta1_y_xzzz_xxyzz_0, \
                             ta1_y_xzzz_xxzzz_0, \
                             ta1_y_xzzz_xyyyy_0, \
                             ta1_y_xzzz_xyyyz_0, \
                             ta1_y_xzzz_xyyzz_0, \
                             ta1_y_xzzz_xyzzz_0, \
                             ta1_y_xzzz_xzzzz_0, \
                             ta1_y_xzzz_yyyyy_0, \
                             ta1_y_xzzz_yyyyz_0, \
                             ta1_y_xzzz_yyyzz_0, \
                             ta1_y_xzzz_yyzzz_0, \
                             ta1_y_xzzz_yzzzz_0, \
                             ta1_y_xzzz_zzzzz_0, \
                             ta1_y_zzz_xxxx_0,   \
                             ta1_y_zzz_xxxx_1,   \
                             ta1_y_zzz_xxxxx_0,  \
                             ta1_y_zzz_xxxxx_1,  \
                             ta1_y_zzz_xxxxy_0,  \
                             ta1_y_zzz_xxxxy_1,  \
                             ta1_y_zzz_xxxxz_0,  \
                             ta1_y_zzz_xxxxz_1,  \
                             ta1_y_zzz_xxxy_0,   \
                             ta1_y_zzz_xxxy_1,   \
                             ta1_y_zzz_xxxyy_0,  \
                             ta1_y_zzz_xxxyy_1,  \
                             ta1_y_zzz_xxxyz_0,  \
                             ta1_y_zzz_xxxyz_1,  \
                             ta1_y_zzz_xxxz_0,   \
                             ta1_y_zzz_xxxz_1,   \
                             ta1_y_zzz_xxxzz_0,  \
                             ta1_y_zzz_xxxzz_1,  \
                             ta1_y_zzz_xxyy_0,   \
                             ta1_y_zzz_xxyy_1,   \
                             ta1_y_zzz_xxyyy_0,  \
                             ta1_y_zzz_xxyyy_1,  \
                             ta1_y_zzz_xxyyz_0,  \
                             ta1_y_zzz_xxyyz_1,  \
                             ta1_y_zzz_xxyz_0,   \
                             ta1_y_zzz_xxyz_1,   \
                             ta1_y_zzz_xxyzz_0,  \
                             ta1_y_zzz_xxyzz_1,  \
                             ta1_y_zzz_xxzz_0,   \
                             ta1_y_zzz_xxzz_1,   \
                             ta1_y_zzz_xxzzz_0,  \
                             ta1_y_zzz_xxzzz_1,  \
                             ta1_y_zzz_xyyy_0,   \
                             ta1_y_zzz_xyyy_1,   \
                             ta1_y_zzz_xyyyy_0,  \
                             ta1_y_zzz_xyyyy_1,  \
                             ta1_y_zzz_xyyyz_0,  \
                             ta1_y_zzz_xyyyz_1,  \
                             ta1_y_zzz_xyyz_0,   \
                             ta1_y_zzz_xyyz_1,   \
                             ta1_y_zzz_xyyzz_0,  \
                             ta1_y_zzz_xyyzz_1,  \
                             ta1_y_zzz_xyzz_0,   \
                             ta1_y_zzz_xyzz_1,   \
                             ta1_y_zzz_xyzzz_0,  \
                             ta1_y_zzz_xyzzz_1,  \
                             ta1_y_zzz_xzzz_0,   \
                             ta1_y_zzz_xzzz_1,   \
                             ta1_y_zzz_xzzzz_0,  \
                             ta1_y_zzz_xzzzz_1,  \
                             ta1_y_zzz_yyyy_0,   \
                             ta1_y_zzz_yyyy_1,   \
                             ta1_y_zzz_yyyyy_0,  \
                             ta1_y_zzz_yyyyy_1,  \
                             ta1_y_zzz_yyyyz_0,  \
                             ta1_y_zzz_yyyyz_1,  \
                             ta1_y_zzz_yyyz_0,   \
                             ta1_y_zzz_yyyz_1,   \
                             ta1_y_zzz_yyyzz_0,  \
                             ta1_y_zzz_yyyzz_1,  \
                             ta1_y_zzz_yyzz_0,   \
                             ta1_y_zzz_yyzz_1,   \
                             ta1_y_zzz_yyzzz_0,  \
                             ta1_y_zzz_yyzzz_1,  \
                             ta1_y_zzz_yzzz_0,   \
                             ta1_y_zzz_yzzz_1,   \
                             ta1_y_zzz_yzzzz_0,  \
                             ta1_y_zzz_yzzzz_1,  \
                             ta1_y_zzz_zzzz_0,   \
                             ta1_y_zzz_zzzz_1,   \
                             ta1_y_zzz_zzzzz_0,  \
                             ta1_y_zzz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_xzzz_xxxxx_0[i] =
            5.0 * ta1_y_zzz_xxxx_0[i] * fe_0 - 5.0 * ta1_y_zzz_xxxx_1[i] * fe_0 + ta1_y_zzz_xxxxx_0[i] * pa_x[i] - ta1_y_zzz_xxxxx_1[i] * pc_x[i];

        ta1_y_xzzz_xxxxy_0[i] =
            4.0 * ta1_y_zzz_xxxy_0[i] * fe_0 - 4.0 * ta1_y_zzz_xxxy_1[i] * fe_0 + ta1_y_zzz_xxxxy_0[i] * pa_x[i] - ta1_y_zzz_xxxxy_1[i] * pc_x[i];

        ta1_y_xzzz_xxxxz_0[i] =
            4.0 * ta1_y_zzz_xxxz_0[i] * fe_0 - 4.0 * ta1_y_zzz_xxxz_1[i] * fe_0 + ta1_y_zzz_xxxxz_0[i] * pa_x[i] - ta1_y_zzz_xxxxz_1[i] * pc_x[i];

        ta1_y_xzzz_xxxyy_0[i] =
            3.0 * ta1_y_zzz_xxyy_0[i] * fe_0 - 3.0 * ta1_y_zzz_xxyy_1[i] * fe_0 + ta1_y_zzz_xxxyy_0[i] * pa_x[i] - ta1_y_zzz_xxxyy_1[i] * pc_x[i];

        ta1_y_xzzz_xxxyz_0[i] =
            3.0 * ta1_y_zzz_xxyz_0[i] * fe_0 - 3.0 * ta1_y_zzz_xxyz_1[i] * fe_0 + ta1_y_zzz_xxxyz_0[i] * pa_x[i] - ta1_y_zzz_xxxyz_1[i] * pc_x[i];

        ta1_y_xzzz_xxxzz_0[i] =
            3.0 * ta1_y_zzz_xxzz_0[i] * fe_0 - 3.0 * ta1_y_zzz_xxzz_1[i] * fe_0 + ta1_y_zzz_xxxzz_0[i] * pa_x[i] - ta1_y_zzz_xxxzz_1[i] * pc_x[i];

        ta1_y_xzzz_xxyyy_0[i] =
            2.0 * ta1_y_zzz_xyyy_0[i] * fe_0 - 2.0 * ta1_y_zzz_xyyy_1[i] * fe_0 + ta1_y_zzz_xxyyy_0[i] * pa_x[i] - ta1_y_zzz_xxyyy_1[i] * pc_x[i];

        ta1_y_xzzz_xxyyz_0[i] =
            2.0 * ta1_y_zzz_xyyz_0[i] * fe_0 - 2.0 * ta1_y_zzz_xyyz_1[i] * fe_0 + ta1_y_zzz_xxyyz_0[i] * pa_x[i] - ta1_y_zzz_xxyyz_1[i] * pc_x[i];

        ta1_y_xzzz_xxyzz_0[i] =
            2.0 * ta1_y_zzz_xyzz_0[i] * fe_0 - 2.0 * ta1_y_zzz_xyzz_1[i] * fe_0 + ta1_y_zzz_xxyzz_0[i] * pa_x[i] - ta1_y_zzz_xxyzz_1[i] * pc_x[i];

        ta1_y_xzzz_xxzzz_0[i] =
            2.0 * ta1_y_zzz_xzzz_0[i] * fe_0 - 2.0 * ta1_y_zzz_xzzz_1[i] * fe_0 + ta1_y_zzz_xxzzz_0[i] * pa_x[i] - ta1_y_zzz_xxzzz_1[i] * pc_x[i];

        ta1_y_xzzz_xyyyy_0[i] =
            ta1_y_zzz_yyyy_0[i] * fe_0 - ta1_y_zzz_yyyy_1[i] * fe_0 + ta1_y_zzz_xyyyy_0[i] * pa_x[i] - ta1_y_zzz_xyyyy_1[i] * pc_x[i];

        ta1_y_xzzz_xyyyz_0[i] =
            ta1_y_zzz_yyyz_0[i] * fe_0 - ta1_y_zzz_yyyz_1[i] * fe_0 + ta1_y_zzz_xyyyz_0[i] * pa_x[i] - ta1_y_zzz_xyyyz_1[i] * pc_x[i];

        ta1_y_xzzz_xyyzz_0[i] =
            ta1_y_zzz_yyzz_0[i] * fe_0 - ta1_y_zzz_yyzz_1[i] * fe_0 + ta1_y_zzz_xyyzz_0[i] * pa_x[i] - ta1_y_zzz_xyyzz_1[i] * pc_x[i];

        ta1_y_xzzz_xyzzz_0[i] =
            ta1_y_zzz_yzzz_0[i] * fe_0 - ta1_y_zzz_yzzz_1[i] * fe_0 + ta1_y_zzz_xyzzz_0[i] * pa_x[i] - ta1_y_zzz_xyzzz_1[i] * pc_x[i];

        ta1_y_xzzz_xzzzz_0[i] =
            ta1_y_zzz_zzzz_0[i] * fe_0 - ta1_y_zzz_zzzz_1[i] * fe_0 + ta1_y_zzz_xzzzz_0[i] * pa_x[i] - ta1_y_zzz_xzzzz_1[i] * pc_x[i];

        ta1_y_xzzz_yyyyy_0[i] = ta1_y_zzz_yyyyy_0[i] * pa_x[i] - ta1_y_zzz_yyyyy_1[i] * pc_x[i];

        ta1_y_xzzz_yyyyz_0[i] = ta1_y_zzz_yyyyz_0[i] * pa_x[i] - ta1_y_zzz_yyyyz_1[i] * pc_x[i];

        ta1_y_xzzz_yyyzz_0[i] = ta1_y_zzz_yyyzz_0[i] * pa_x[i] - ta1_y_zzz_yyyzz_1[i] * pc_x[i];

        ta1_y_xzzz_yyzzz_0[i] = ta1_y_zzz_yyzzz_0[i] * pa_x[i] - ta1_y_zzz_yyzzz_1[i] * pc_x[i];

        ta1_y_xzzz_yzzzz_0[i] = ta1_y_zzz_yzzzz_0[i] * pa_x[i] - ta1_y_zzz_yzzzz_1[i] * pc_x[i];

        ta1_y_xzzz_zzzzz_0[i] = ta1_y_zzz_zzzzz_0[i] * pa_x[i] - ta1_y_zzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 525-546 components of targeted buffer : GH

    auto ta1_y_yyyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 525);

    auto ta1_y_yyyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 526);

    auto ta1_y_yyyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 527);

    auto ta1_y_yyyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 528);

    auto ta1_y_yyyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 529);

    auto ta1_y_yyyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 530);

    auto ta1_y_yyyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 531);

    auto ta1_y_yyyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 532);

    auto ta1_y_yyyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 533);

    auto ta1_y_yyyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 534);

    auto ta1_y_yyyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 535);

    auto ta1_y_yyyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 536);

    auto ta1_y_yyyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 537);

    auto ta1_y_yyyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 538);

    auto ta1_y_yyyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 539);

    auto ta1_y_yyyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 540);

    auto ta1_y_yyyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 541);

    auto ta1_y_yyyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 542);

    auto ta1_y_yyyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 543);

    auto ta1_y_yyyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 544);

    auto ta1_y_yyyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 545);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_y_yy_xxxxx_0,   \
                             ta1_y_yy_xxxxx_1,   \
                             ta1_y_yy_xxxxy_0,   \
                             ta1_y_yy_xxxxy_1,   \
                             ta1_y_yy_xxxxz_0,   \
                             ta1_y_yy_xxxxz_1,   \
                             ta1_y_yy_xxxyy_0,   \
                             ta1_y_yy_xxxyy_1,   \
                             ta1_y_yy_xxxyz_0,   \
                             ta1_y_yy_xxxyz_1,   \
                             ta1_y_yy_xxxzz_0,   \
                             ta1_y_yy_xxxzz_1,   \
                             ta1_y_yy_xxyyy_0,   \
                             ta1_y_yy_xxyyy_1,   \
                             ta1_y_yy_xxyyz_0,   \
                             ta1_y_yy_xxyyz_1,   \
                             ta1_y_yy_xxyzz_0,   \
                             ta1_y_yy_xxyzz_1,   \
                             ta1_y_yy_xxzzz_0,   \
                             ta1_y_yy_xxzzz_1,   \
                             ta1_y_yy_xyyyy_0,   \
                             ta1_y_yy_xyyyy_1,   \
                             ta1_y_yy_xyyyz_0,   \
                             ta1_y_yy_xyyyz_1,   \
                             ta1_y_yy_xyyzz_0,   \
                             ta1_y_yy_xyyzz_1,   \
                             ta1_y_yy_xyzzz_0,   \
                             ta1_y_yy_xyzzz_1,   \
                             ta1_y_yy_xzzzz_0,   \
                             ta1_y_yy_xzzzz_1,   \
                             ta1_y_yy_yyyyy_0,   \
                             ta1_y_yy_yyyyy_1,   \
                             ta1_y_yy_yyyyz_0,   \
                             ta1_y_yy_yyyyz_1,   \
                             ta1_y_yy_yyyzz_0,   \
                             ta1_y_yy_yyyzz_1,   \
                             ta1_y_yy_yyzzz_0,   \
                             ta1_y_yy_yyzzz_1,   \
                             ta1_y_yy_yzzzz_0,   \
                             ta1_y_yy_yzzzz_1,   \
                             ta1_y_yy_zzzzz_0,   \
                             ta1_y_yy_zzzzz_1,   \
                             ta1_y_yyy_xxxx_0,   \
                             ta1_y_yyy_xxxx_1,   \
                             ta1_y_yyy_xxxxx_0,  \
                             ta1_y_yyy_xxxxx_1,  \
                             ta1_y_yyy_xxxxy_0,  \
                             ta1_y_yyy_xxxxy_1,  \
                             ta1_y_yyy_xxxxz_0,  \
                             ta1_y_yyy_xxxxz_1,  \
                             ta1_y_yyy_xxxy_0,   \
                             ta1_y_yyy_xxxy_1,   \
                             ta1_y_yyy_xxxyy_0,  \
                             ta1_y_yyy_xxxyy_1,  \
                             ta1_y_yyy_xxxyz_0,  \
                             ta1_y_yyy_xxxyz_1,  \
                             ta1_y_yyy_xxxz_0,   \
                             ta1_y_yyy_xxxz_1,   \
                             ta1_y_yyy_xxxzz_0,  \
                             ta1_y_yyy_xxxzz_1,  \
                             ta1_y_yyy_xxyy_0,   \
                             ta1_y_yyy_xxyy_1,   \
                             ta1_y_yyy_xxyyy_0,  \
                             ta1_y_yyy_xxyyy_1,  \
                             ta1_y_yyy_xxyyz_0,  \
                             ta1_y_yyy_xxyyz_1,  \
                             ta1_y_yyy_xxyz_0,   \
                             ta1_y_yyy_xxyz_1,   \
                             ta1_y_yyy_xxyzz_0,  \
                             ta1_y_yyy_xxyzz_1,  \
                             ta1_y_yyy_xxzz_0,   \
                             ta1_y_yyy_xxzz_1,   \
                             ta1_y_yyy_xxzzz_0,  \
                             ta1_y_yyy_xxzzz_1,  \
                             ta1_y_yyy_xyyy_0,   \
                             ta1_y_yyy_xyyy_1,   \
                             ta1_y_yyy_xyyyy_0,  \
                             ta1_y_yyy_xyyyy_1,  \
                             ta1_y_yyy_xyyyz_0,  \
                             ta1_y_yyy_xyyyz_1,  \
                             ta1_y_yyy_xyyz_0,   \
                             ta1_y_yyy_xyyz_1,   \
                             ta1_y_yyy_xyyzz_0,  \
                             ta1_y_yyy_xyyzz_1,  \
                             ta1_y_yyy_xyzz_0,   \
                             ta1_y_yyy_xyzz_1,   \
                             ta1_y_yyy_xyzzz_0,  \
                             ta1_y_yyy_xyzzz_1,  \
                             ta1_y_yyy_xzzz_0,   \
                             ta1_y_yyy_xzzz_1,   \
                             ta1_y_yyy_xzzzz_0,  \
                             ta1_y_yyy_xzzzz_1,  \
                             ta1_y_yyy_yyyy_0,   \
                             ta1_y_yyy_yyyy_1,   \
                             ta1_y_yyy_yyyyy_0,  \
                             ta1_y_yyy_yyyyy_1,  \
                             ta1_y_yyy_yyyyz_0,  \
                             ta1_y_yyy_yyyyz_1,  \
                             ta1_y_yyy_yyyz_0,   \
                             ta1_y_yyy_yyyz_1,   \
                             ta1_y_yyy_yyyzz_0,  \
                             ta1_y_yyy_yyyzz_1,  \
                             ta1_y_yyy_yyzz_0,   \
                             ta1_y_yyy_yyzz_1,   \
                             ta1_y_yyy_yyzzz_0,  \
                             ta1_y_yyy_yyzzz_1,  \
                             ta1_y_yyy_yzzz_0,   \
                             ta1_y_yyy_yzzz_1,   \
                             ta1_y_yyy_yzzzz_0,  \
                             ta1_y_yyy_yzzzz_1,  \
                             ta1_y_yyy_zzzz_0,   \
                             ta1_y_yyy_zzzz_1,   \
                             ta1_y_yyy_zzzzz_0,  \
                             ta1_y_yyy_zzzzz_1,  \
                             ta1_y_yyyy_xxxxx_0, \
                             ta1_y_yyyy_xxxxy_0, \
                             ta1_y_yyyy_xxxxz_0, \
                             ta1_y_yyyy_xxxyy_0, \
                             ta1_y_yyyy_xxxyz_0, \
                             ta1_y_yyyy_xxxzz_0, \
                             ta1_y_yyyy_xxyyy_0, \
                             ta1_y_yyyy_xxyyz_0, \
                             ta1_y_yyyy_xxyzz_0, \
                             ta1_y_yyyy_xxzzz_0, \
                             ta1_y_yyyy_xyyyy_0, \
                             ta1_y_yyyy_xyyyz_0, \
                             ta1_y_yyyy_xyyzz_0, \
                             ta1_y_yyyy_xyzzz_0, \
                             ta1_y_yyyy_xzzzz_0, \
                             ta1_y_yyyy_yyyyy_0, \
                             ta1_y_yyyy_yyyyz_0, \
                             ta1_y_yyyy_yyyzz_0, \
                             ta1_y_yyyy_yyzzz_0, \
                             ta1_y_yyyy_yzzzz_0, \
                             ta1_y_yyyy_zzzzz_0, \
                             ta_yyy_xxxxx_1,     \
                             ta_yyy_xxxxy_1,     \
                             ta_yyy_xxxxz_1,     \
                             ta_yyy_xxxyy_1,     \
                             ta_yyy_xxxyz_1,     \
                             ta_yyy_xxxzz_1,     \
                             ta_yyy_xxyyy_1,     \
                             ta_yyy_xxyyz_1,     \
                             ta_yyy_xxyzz_1,     \
                             ta_yyy_xxzzz_1,     \
                             ta_yyy_xyyyy_1,     \
                             ta_yyy_xyyyz_1,     \
                             ta_yyy_xyyzz_1,     \
                             ta_yyy_xyzzz_1,     \
                             ta_yyy_xzzzz_1,     \
                             ta_yyy_yyyyy_1,     \
                             ta_yyy_yyyyz_1,     \
                             ta_yyy_yyyzz_1,     \
                             ta_yyy_yyzzz_1,     \
                             ta_yyy_yzzzz_1,     \
                             ta_yyy_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyy_xxxxx_0[i] = 3.0 * ta1_y_yy_xxxxx_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxxx_1[i] * fe_0 + ta_yyy_xxxxx_1[i] +
                                ta1_y_yyy_xxxxx_0[i] * pa_y[i] - ta1_y_yyy_xxxxx_1[i] * pc_y[i];

        ta1_y_yyyy_xxxxy_0[i] = 3.0 * ta1_y_yy_xxxxy_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxxy_1[i] * fe_0 + ta1_y_yyy_xxxx_0[i] * fe_0 -
                                ta1_y_yyy_xxxx_1[i] * fe_0 + ta_yyy_xxxxy_1[i] + ta1_y_yyy_xxxxy_0[i] * pa_y[i] - ta1_y_yyy_xxxxy_1[i] * pc_y[i];

        ta1_y_yyyy_xxxxz_0[i] = 3.0 * ta1_y_yy_xxxxz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxxz_1[i] * fe_0 + ta_yyy_xxxxz_1[i] +
                                ta1_y_yyy_xxxxz_0[i] * pa_y[i] - ta1_y_yyy_xxxxz_1[i] * pc_y[i];

        ta1_y_yyyy_xxxyy_0[i] = 3.0 * ta1_y_yy_xxxyy_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxyy_1[i] * fe_0 + 2.0 * ta1_y_yyy_xxxy_0[i] * fe_0 -
                                2.0 * ta1_y_yyy_xxxy_1[i] * fe_0 + ta_yyy_xxxyy_1[i] + ta1_y_yyy_xxxyy_0[i] * pa_y[i] -
                                ta1_y_yyy_xxxyy_1[i] * pc_y[i];

        ta1_y_yyyy_xxxyz_0[i] = 3.0 * ta1_y_yy_xxxyz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxyz_1[i] * fe_0 + ta1_y_yyy_xxxz_0[i] * fe_0 -
                                ta1_y_yyy_xxxz_1[i] * fe_0 + ta_yyy_xxxyz_1[i] + ta1_y_yyy_xxxyz_0[i] * pa_y[i] - ta1_y_yyy_xxxyz_1[i] * pc_y[i];

        ta1_y_yyyy_xxxzz_0[i] = 3.0 * ta1_y_yy_xxxzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxxzz_1[i] * fe_0 + ta_yyy_xxxzz_1[i] +
                                ta1_y_yyy_xxxzz_0[i] * pa_y[i] - ta1_y_yyy_xxxzz_1[i] * pc_y[i];

        ta1_y_yyyy_xxyyy_0[i] = 3.0 * ta1_y_yy_xxyyy_0[i] * fe_0 - 3.0 * ta1_y_yy_xxyyy_1[i] * fe_0 + 3.0 * ta1_y_yyy_xxyy_0[i] * fe_0 -
                                3.0 * ta1_y_yyy_xxyy_1[i] * fe_0 + ta_yyy_xxyyy_1[i] + ta1_y_yyy_xxyyy_0[i] * pa_y[i] -
                                ta1_y_yyy_xxyyy_1[i] * pc_y[i];

        ta1_y_yyyy_xxyyz_0[i] = 3.0 * ta1_y_yy_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxyyz_1[i] * fe_0 + 2.0 * ta1_y_yyy_xxyz_0[i] * fe_0 -
                                2.0 * ta1_y_yyy_xxyz_1[i] * fe_0 + ta_yyy_xxyyz_1[i] + ta1_y_yyy_xxyyz_0[i] * pa_y[i] -
                                ta1_y_yyy_xxyyz_1[i] * pc_y[i];

        ta1_y_yyyy_xxyzz_0[i] = 3.0 * ta1_y_yy_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxyzz_1[i] * fe_0 + ta1_y_yyy_xxzz_0[i] * fe_0 -
                                ta1_y_yyy_xxzz_1[i] * fe_0 + ta_yyy_xxyzz_1[i] + ta1_y_yyy_xxyzz_0[i] * pa_y[i] - ta1_y_yyy_xxyzz_1[i] * pc_y[i];

        ta1_y_yyyy_xxzzz_0[i] = 3.0 * ta1_y_yy_xxzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xxzzz_1[i] * fe_0 + ta_yyy_xxzzz_1[i] +
                                ta1_y_yyy_xxzzz_0[i] * pa_y[i] - ta1_y_yyy_xxzzz_1[i] * pc_y[i];

        ta1_y_yyyy_xyyyy_0[i] = 3.0 * ta1_y_yy_xyyyy_0[i] * fe_0 - 3.0 * ta1_y_yy_xyyyy_1[i] * fe_0 + 4.0 * ta1_y_yyy_xyyy_0[i] * fe_0 -
                                4.0 * ta1_y_yyy_xyyy_1[i] * fe_0 + ta_yyy_xyyyy_1[i] + ta1_y_yyy_xyyyy_0[i] * pa_y[i] -
                                ta1_y_yyy_xyyyy_1[i] * pc_y[i];

        ta1_y_yyyy_xyyyz_0[i] = 3.0 * ta1_y_yy_xyyyz_0[i] * fe_0 - 3.0 * ta1_y_yy_xyyyz_1[i] * fe_0 + 3.0 * ta1_y_yyy_xyyz_0[i] * fe_0 -
                                3.0 * ta1_y_yyy_xyyz_1[i] * fe_0 + ta_yyy_xyyyz_1[i] + ta1_y_yyy_xyyyz_0[i] * pa_y[i] -
                                ta1_y_yyy_xyyyz_1[i] * pc_y[i];

        ta1_y_yyyy_xyyzz_0[i] = 3.0 * ta1_y_yy_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xyyzz_1[i] * fe_0 + 2.0 * ta1_y_yyy_xyzz_0[i] * fe_0 -
                                2.0 * ta1_y_yyy_xyzz_1[i] * fe_0 + ta_yyy_xyyzz_1[i] + ta1_y_yyy_xyyzz_0[i] * pa_y[i] -
                                ta1_y_yyy_xyyzz_1[i] * pc_y[i];

        ta1_y_yyyy_xyzzz_0[i] = 3.0 * ta1_y_yy_xyzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xyzzz_1[i] * fe_0 + ta1_y_yyy_xzzz_0[i] * fe_0 -
                                ta1_y_yyy_xzzz_1[i] * fe_0 + ta_yyy_xyzzz_1[i] + ta1_y_yyy_xyzzz_0[i] * pa_y[i] - ta1_y_yyy_xyzzz_1[i] * pc_y[i];

        ta1_y_yyyy_xzzzz_0[i] = 3.0 * ta1_y_yy_xzzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_xzzzz_1[i] * fe_0 + ta_yyy_xzzzz_1[i] +
                                ta1_y_yyy_xzzzz_0[i] * pa_y[i] - ta1_y_yyy_xzzzz_1[i] * pc_y[i];

        ta1_y_yyyy_yyyyy_0[i] = 3.0 * ta1_y_yy_yyyyy_0[i] * fe_0 - 3.0 * ta1_y_yy_yyyyy_1[i] * fe_0 + 5.0 * ta1_y_yyy_yyyy_0[i] * fe_0 -
                                5.0 * ta1_y_yyy_yyyy_1[i] * fe_0 + ta_yyy_yyyyy_1[i] + ta1_y_yyy_yyyyy_0[i] * pa_y[i] -
                                ta1_y_yyy_yyyyy_1[i] * pc_y[i];

        ta1_y_yyyy_yyyyz_0[i] = 3.0 * ta1_y_yy_yyyyz_0[i] * fe_0 - 3.0 * ta1_y_yy_yyyyz_1[i] * fe_0 + 4.0 * ta1_y_yyy_yyyz_0[i] * fe_0 -
                                4.0 * ta1_y_yyy_yyyz_1[i] * fe_0 + ta_yyy_yyyyz_1[i] + ta1_y_yyy_yyyyz_0[i] * pa_y[i] -
                                ta1_y_yyy_yyyyz_1[i] * pc_y[i];

        ta1_y_yyyy_yyyzz_0[i] = 3.0 * ta1_y_yy_yyyzz_0[i] * fe_0 - 3.0 * ta1_y_yy_yyyzz_1[i] * fe_0 + 3.0 * ta1_y_yyy_yyzz_0[i] * fe_0 -
                                3.0 * ta1_y_yyy_yyzz_1[i] * fe_0 + ta_yyy_yyyzz_1[i] + ta1_y_yyy_yyyzz_0[i] * pa_y[i] -
                                ta1_y_yyy_yyyzz_1[i] * pc_y[i];

        ta1_y_yyyy_yyzzz_0[i] = 3.0 * ta1_y_yy_yyzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_yyzzz_1[i] * fe_0 + 2.0 * ta1_y_yyy_yzzz_0[i] * fe_0 -
                                2.0 * ta1_y_yyy_yzzz_1[i] * fe_0 + ta_yyy_yyzzz_1[i] + ta1_y_yyy_yyzzz_0[i] * pa_y[i] -
                                ta1_y_yyy_yyzzz_1[i] * pc_y[i];

        ta1_y_yyyy_yzzzz_0[i] = 3.0 * ta1_y_yy_yzzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_yzzzz_1[i] * fe_0 + ta1_y_yyy_zzzz_0[i] * fe_0 -
                                ta1_y_yyy_zzzz_1[i] * fe_0 + ta_yyy_yzzzz_1[i] + ta1_y_yyy_yzzzz_0[i] * pa_y[i] - ta1_y_yyy_yzzzz_1[i] * pc_y[i];

        ta1_y_yyyy_zzzzz_0[i] = 3.0 * ta1_y_yy_zzzzz_0[i] * fe_0 - 3.0 * ta1_y_yy_zzzzz_1[i] * fe_0 + ta_yyy_zzzzz_1[i] +
                                ta1_y_yyy_zzzzz_0[i] * pa_y[i] - ta1_y_yyy_zzzzz_1[i] * pc_y[i];
    }

    // Set up 546-567 components of targeted buffer : GH

    auto ta1_y_yyyz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 546);

    auto ta1_y_yyyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 547);

    auto ta1_y_yyyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 548);

    auto ta1_y_yyyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 549);

    auto ta1_y_yyyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 550);

    auto ta1_y_yyyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 551);

    auto ta1_y_yyyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 552);

    auto ta1_y_yyyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 553);

    auto ta1_y_yyyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 554);

    auto ta1_y_yyyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 555);

    auto ta1_y_yyyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 556);

    auto ta1_y_yyyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 557);

    auto ta1_y_yyyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 558);

    auto ta1_y_yyyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 559);

    auto ta1_y_yyyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 560);

    auto ta1_y_yyyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 561);

    auto ta1_y_yyyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 562);

    auto ta1_y_yyyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 563);

    auto ta1_y_yyyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 564);

    auto ta1_y_yyyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 565);

    auto ta1_y_yyyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 566);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_y_yyy_xxxx_0,   \
                             ta1_y_yyy_xxxx_1,   \
                             ta1_y_yyy_xxxxx_0,  \
                             ta1_y_yyy_xxxxx_1,  \
                             ta1_y_yyy_xxxxy_0,  \
                             ta1_y_yyy_xxxxy_1,  \
                             ta1_y_yyy_xxxxz_0,  \
                             ta1_y_yyy_xxxxz_1,  \
                             ta1_y_yyy_xxxy_0,   \
                             ta1_y_yyy_xxxy_1,   \
                             ta1_y_yyy_xxxyy_0,  \
                             ta1_y_yyy_xxxyy_1,  \
                             ta1_y_yyy_xxxyz_0,  \
                             ta1_y_yyy_xxxyz_1,  \
                             ta1_y_yyy_xxxz_0,   \
                             ta1_y_yyy_xxxz_1,   \
                             ta1_y_yyy_xxxzz_0,  \
                             ta1_y_yyy_xxxzz_1,  \
                             ta1_y_yyy_xxyy_0,   \
                             ta1_y_yyy_xxyy_1,   \
                             ta1_y_yyy_xxyyy_0,  \
                             ta1_y_yyy_xxyyy_1,  \
                             ta1_y_yyy_xxyyz_0,  \
                             ta1_y_yyy_xxyyz_1,  \
                             ta1_y_yyy_xxyz_0,   \
                             ta1_y_yyy_xxyz_1,   \
                             ta1_y_yyy_xxyzz_0,  \
                             ta1_y_yyy_xxyzz_1,  \
                             ta1_y_yyy_xxzz_0,   \
                             ta1_y_yyy_xxzz_1,   \
                             ta1_y_yyy_xxzzz_0,  \
                             ta1_y_yyy_xxzzz_1,  \
                             ta1_y_yyy_xyyy_0,   \
                             ta1_y_yyy_xyyy_1,   \
                             ta1_y_yyy_xyyyy_0,  \
                             ta1_y_yyy_xyyyy_1,  \
                             ta1_y_yyy_xyyyz_0,  \
                             ta1_y_yyy_xyyyz_1,  \
                             ta1_y_yyy_xyyz_0,   \
                             ta1_y_yyy_xyyz_1,   \
                             ta1_y_yyy_xyyzz_0,  \
                             ta1_y_yyy_xyyzz_1,  \
                             ta1_y_yyy_xyzz_0,   \
                             ta1_y_yyy_xyzz_1,   \
                             ta1_y_yyy_xyzzz_0,  \
                             ta1_y_yyy_xyzzz_1,  \
                             ta1_y_yyy_xzzz_0,   \
                             ta1_y_yyy_xzzz_1,   \
                             ta1_y_yyy_xzzzz_0,  \
                             ta1_y_yyy_xzzzz_1,  \
                             ta1_y_yyy_yyyy_0,   \
                             ta1_y_yyy_yyyy_1,   \
                             ta1_y_yyy_yyyyy_0,  \
                             ta1_y_yyy_yyyyy_1,  \
                             ta1_y_yyy_yyyyz_0,  \
                             ta1_y_yyy_yyyyz_1,  \
                             ta1_y_yyy_yyyz_0,   \
                             ta1_y_yyy_yyyz_1,   \
                             ta1_y_yyy_yyyzz_0,  \
                             ta1_y_yyy_yyyzz_1,  \
                             ta1_y_yyy_yyzz_0,   \
                             ta1_y_yyy_yyzz_1,   \
                             ta1_y_yyy_yyzzz_0,  \
                             ta1_y_yyy_yyzzz_1,  \
                             ta1_y_yyy_yzzz_0,   \
                             ta1_y_yyy_yzzz_1,   \
                             ta1_y_yyy_yzzzz_0,  \
                             ta1_y_yyy_yzzzz_1,  \
                             ta1_y_yyy_zzzz_0,   \
                             ta1_y_yyy_zzzz_1,   \
                             ta1_y_yyy_zzzzz_0,  \
                             ta1_y_yyy_zzzzz_1,  \
                             ta1_y_yyyz_xxxxx_0, \
                             ta1_y_yyyz_xxxxy_0, \
                             ta1_y_yyyz_xxxxz_0, \
                             ta1_y_yyyz_xxxyy_0, \
                             ta1_y_yyyz_xxxyz_0, \
                             ta1_y_yyyz_xxxzz_0, \
                             ta1_y_yyyz_xxyyy_0, \
                             ta1_y_yyyz_xxyyz_0, \
                             ta1_y_yyyz_xxyzz_0, \
                             ta1_y_yyyz_xxzzz_0, \
                             ta1_y_yyyz_xyyyy_0, \
                             ta1_y_yyyz_xyyyz_0, \
                             ta1_y_yyyz_xyyzz_0, \
                             ta1_y_yyyz_xyzzz_0, \
                             ta1_y_yyyz_xzzzz_0, \
                             ta1_y_yyyz_yyyyy_0, \
                             ta1_y_yyyz_yyyyz_0, \
                             ta1_y_yyyz_yyyzz_0, \
                             ta1_y_yyyz_yyzzz_0, \
                             ta1_y_yyyz_yzzzz_0, \
                             ta1_y_yyyz_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyyz_xxxxx_0[i] = ta1_y_yyy_xxxxx_0[i] * pa_z[i] - ta1_y_yyy_xxxxx_1[i] * pc_z[i];

        ta1_y_yyyz_xxxxy_0[i] = ta1_y_yyy_xxxxy_0[i] * pa_z[i] - ta1_y_yyy_xxxxy_1[i] * pc_z[i];

        ta1_y_yyyz_xxxxz_0[i] =
            ta1_y_yyy_xxxx_0[i] * fe_0 - ta1_y_yyy_xxxx_1[i] * fe_0 + ta1_y_yyy_xxxxz_0[i] * pa_z[i] - ta1_y_yyy_xxxxz_1[i] * pc_z[i];

        ta1_y_yyyz_xxxyy_0[i] = ta1_y_yyy_xxxyy_0[i] * pa_z[i] - ta1_y_yyy_xxxyy_1[i] * pc_z[i];

        ta1_y_yyyz_xxxyz_0[i] =
            ta1_y_yyy_xxxy_0[i] * fe_0 - ta1_y_yyy_xxxy_1[i] * fe_0 + ta1_y_yyy_xxxyz_0[i] * pa_z[i] - ta1_y_yyy_xxxyz_1[i] * pc_z[i];

        ta1_y_yyyz_xxxzz_0[i] =
            2.0 * ta1_y_yyy_xxxz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xxxz_1[i] * fe_0 + ta1_y_yyy_xxxzz_0[i] * pa_z[i] - ta1_y_yyy_xxxzz_1[i] * pc_z[i];

        ta1_y_yyyz_xxyyy_0[i] = ta1_y_yyy_xxyyy_0[i] * pa_z[i] - ta1_y_yyy_xxyyy_1[i] * pc_z[i];

        ta1_y_yyyz_xxyyz_0[i] =
            ta1_y_yyy_xxyy_0[i] * fe_0 - ta1_y_yyy_xxyy_1[i] * fe_0 + ta1_y_yyy_xxyyz_0[i] * pa_z[i] - ta1_y_yyy_xxyyz_1[i] * pc_z[i];

        ta1_y_yyyz_xxyzz_0[i] =
            2.0 * ta1_y_yyy_xxyz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xxyz_1[i] * fe_0 + ta1_y_yyy_xxyzz_0[i] * pa_z[i] - ta1_y_yyy_xxyzz_1[i] * pc_z[i];

        ta1_y_yyyz_xxzzz_0[i] =
            3.0 * ta1_y_yyy_xxzz_0[i] * fe_0 - 3.0 * ta1_y_yyy_xxzz_1[i] * fe_0 + ta1_y_yyy_xxzzz_0[i] * pa_z[i] - ta1_y_yyy_xxzzz_1[i] * pc_z[i];

        ta1_y_yyyz_xyyyy_0[i] = ta1_y_yyy_xyyyy_0[i] * pa_z[i] - ta1_y_yyy_xyyyy_1[i] * pc_z[i];

        ta1_y_yyyz_xyyyz_0[i] =
            ta1_y_yyy_xyyy_0[i] * fe_0 - ta1_y_yyy_xyyy_1[i] * fe_0 + ta1_y_yyy_xyyyz_0[i] * pa_z[i] - ta1_y_yyy_xyyyz_1[i] * pc_z[i];

        ta1_y_yyyz_xyyzz_0[i] =
            2.0 * ta1_y_yyy_xyyz_0[i] * fe_0 - 2.0 * ta1_y_yyy_xyyz_1[i] * fe_0 + ta1_y_yyy_xyyzz_0[i] * pa_z[i] - ta1_y_yyy_xyyzz_1[i] * pc_z[i];

        ta1_y_yyyz_xyzzz_0[i] =
            3.0 * ta1_y_yyy_xyzz_0[i] * fe_0 - 3.0 * ta1_y_yyy_xyzz_1[i] * fe_0 + ta1_y_yyy_xyzzz_0[i] * pa_z[i] - ta1_y_yyy_xyzzz_1[i] * pc_z[i];

        ta1_y_yyyz_xzzzz_0[i] =
            4.0 * ta1_y_yyy_xzzz_0[i] * fe_0 - 4.0 * ta1_y_yyy_xzzz_1[i] * fe_0 + ta1_y_yyy_xzzzz_0[i] * pa_z[i] - ta1_y_yyy_xzzzz_1[i] * pc_z[i];

        ta1_y_yyyz_yyyyy_0[i] = ta1_y_yyy_yyyyy_0[i] * pa_z[i] - ta1_y_yyy_yyyyy_1[i] * pc_z[i];

        ta1_y_yyyz_yyyyz_0[i] =
            ta1_y_yyy_yyyy_0[i] * fe_0 - ta1_y_yyy_yyyy_1[i] * fe_0 + ta1_y_yyy_yyyyz_0[i] * pa_z[i] - ta1_y_yyy_yyyyz_1[i] * pc_z[i];

        ta1_y_yyyz_yyyzz_0[i] =
            2.0 * ta1_y_yyy_yyyz_0[i] * fe_0 - 2.0 * ta1_y_yyy_yyyz_1[i] * fe_0 + ta1_y_yyy_yyyzz_0[i] * pa_z[i] - ta1_y_yyy_yyyzz_1[i] * pc_z[i];

        ta1_y_yyyz_yyzzz_0[i] =
            3.0 * ta1_y_yyy_yyzz_0[i] * fe_0 - 3.0 * ta1_y_yyy_yyzz_1[i] * fe_0 + ta1_y_yyy_yyzzz_0[i] * pa_z[i] - ta1_y_yyy_yyzzz_1[i] * pc_z[i];

        ta1_y_yyyz_yzzzz_0[i] =
            4.0 * ta1_y_yyy_yzzz_0[i] * fe_0 - 4.0 * ta1_y_yyy_yzzz_1[i] * fe_0 + ta1_y_yyy_yzzzz_0[i] * pa_z[i] - ta1_y_yyy_yzzzz_1[i] * pc_z[i];

        ta1_y_yyyz_zzzzz_0[i] =
            5.0 * ta1_y_yyy_zzzz_0[i] * fe_0 - 5.0 * ta1_y_yyy_zzzz_1[i] * fe_0 + ta1_y_yyy_zzzzz_0[i] * pa_z[i] - ta1_y_yyy_zzzzz_1[i] * pc_z[i];
    }

    // Set up 567-588 components of targeted buffer : GH

    auto ta1_y_yyzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 567);

    auto ta1_y_yyzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 568);

    auto ta1_y_yyzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 569);

    auto ta1_y_yyzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 570);

    auto ta1_y_yyzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 571);

    auto ta1_y_yyzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 572);

    auto ta1_y_yyzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 573);

    auto ta1_y_yyzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 574);

    auto ta1_y_yyzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 575);

    auto ta1_y_yyzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 576);

    auto ta1_y_yyzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 577);

    auto ta1_y_yyzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 578);

    auto ta1_y_yyzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 579);

    auto ta1_y_yyzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 580);

    auto ta1_y_yyzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 581);

    auto ta1_y_yyzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 582);

    auto ta1_y_yyzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 583);

    auto ta1_y_yyzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 584);

    auto ta1_y_yyzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 585);

    auto ta1_y_yyzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 586);

    auto ta1_y_yyzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 587);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_yy_xxxxx_0,   \
                             ta1_y_yy_xxxxx_1,   \
                             ta1_y_yy_xxxxy_0,   \
                             ta1_y_yy_xxxxy_1,   \
                             ta1_y_yy_xxxyy_0,   \
                             ta1_y_yy_xxxyy_1,   \
                             ta1_y_yy_xxxyz_0,   \
                             ta1_y_yy_xxxyz_1,   \
                             ta1_y_yy_xxyyy_0,   \
                             ta1_y_yy_xxyyy_1,   \
                             ta1_y_yy_xxyyz_0,   \
                             ta1_y_yy_xxyyz_1,   \
                             ta1_y_yy_xxyzz_0,   \
                             ta1_y_yy_xxyzz_1,   \
                             ta1_y_yy_xyyyy_0,   \
                             ta1_y_yy_xyyyy_1,   \
                             ta1_y_yy_xyyyz_0,   \
                             ta1_y_yy_xyyyz_1,   \
                             ta1_y_yy_xyyzz_0,   \
                             ta1_y_yy_xyyzz_1,   \
                             ta1_y_yy_xyzzz_0,   \
                             ta1_y_yy_xyzzz_1,   \
                             ta1_y_yy_yyyyy_0,   \
                             ta1_y_yy_yyyyy_1,   \
                             ta1_y_yy_yyyyz_0,   \
                             ta1_y_yy_yyyyz_1,   \
                             ta1_y_yy_yyyzz_0,   \
                             ta1_y_yy_yyyzz_1,   \
                             ta1_y_yy_yyzzz_0,   \
                             ta1_y_yy_yyzzz_1,   \
                             ta1_y_yy_yzzzz_0,   \
                             ta1_y_yy_yzzzz_1,   \
                             ta1_y_yyz_xxxxx_0,  \
                             ta1_y_yyz_xxxxx_1,  \
                             ta1_y_yyz_xxxxy_0,  \
                             ta1_y_yyz_xxxxy_1,  \
                             ta1_y_yyz_xxxy_0,   \
                             ta1_y_yyz_xxxy_1,   \
                             ta1_y_yyz_xxxyy_0,  \
                             ta1_y_yyz_xxxyy_1,  \
                             ta1_y_yyz_xxxyz_0,  \
                             ta1_y_yyz_xxxyz_1,  \
                             ta1_y_yyz_xxyy_0,   \
                             ta1_y_yyz_xxyy_1,   \
                             ta1_y_yyz_xxyyy_0,  \
                             ta1_y_yyz_xxyyy_1,  \
                             ta1_y_yyz_xxyyz_0,  \
                             ta1_y_yyz_xxyyz_1,  \
                             ta1_y_yyz_xxyz_0,   \
                             ta1_y_yyz_xxyz_1,   \
                             ta1_y_yyz_xxyzz_0,  \
                             ta1_y_yyz_xxyzz_1,  \
                             ta1_y_yyz_xyyy_0,   \
                             ta1_y_yyz_xyyy_1,   \
                             ta1_y_yyz_xyyyy_0,  \
                             ta1_y_yyz_xyyyy_1,  \
                             ta1_y_yyz_xyyyz_0,  \
                             ta1_y_yyz_xyyyz_1,  \
                             ta1_y_yyz_xyyz_0,   \
                             ta1_y_yyz_xyyz_1,   \
                             ta1_y_yyz_xyyzz_0,  \
                             ta1_y_yyz_xyyzz_1,  \
                             ta1_y_yyz_xyzz_0,   \
                             ta1_y_yyz_xyzz_1,   \
                             ta1_y_yyz_xyzzz_0,  \
                             ta1_y_yyz_xyzzz_1,  \
                             ta1_y_yyz_yyyy_0,   \
                             ta1_y_yyz_yyyy_1,   \
                             ta1_y_yyz_yyyyy_0,  \
                             ta1_y_yyz_yyyyy_1,  \
                             ta1_y_yyz_yyyyz_0,  \
                             ta1_y_yyz_yyyyz_1,  \
                             ta1_y_yyz_yyyz_0,   \
                             ta1_y_yyz_yyyz_1,   \
                             ta1_y_yyz_yyyzz_0,  \
                             ta1_y_yyz_yyyzz_1,  \
                             ta1_y_yyz_yyzz_0,   \
                             ta1_y_yyz_yyzz_1,   \
                             ta1_y_yyz_yyzzz_0,  \
                             ta1_y_yyz_yyzzz_1,  \
                             ta1_y_yyz_yzzz_0,   \
                             ta1_y_yyz_yzzz_1,   \
                             ta1_y_yyz_yzzzz_0,  \
                             ta1_y_yyz_yzzzz_1,  \
                             ta1_y_yyzz_xxxxx_0, \
                             ta1_y_yyzz_xxxxy_0, \
                             ta1_y_yyzz_xxxxz_0, \
                             ta1_y_yyzz_xxxyy_0, \
                             ta1_y_yyzz_xxxyz_0, \
                             ta1_y_yyzz_xxxzz_0, \
                             ta1_y_yyzz_xxyyy_0, \
                             ta1_y_yyzz_xxyyz_0, \
                             ta1_y_yyzz_xxyzz_0, \
                             ta1_y_yyzz_xxzzz_0, \
                             ta1_y_yyzz_xyyyy_0, \
                             ta1_y_yyzz_xyyyz_0, \
                             ta1_y_yyzz_xyyzz_0, \
                             ta1_y_yyzz_xyzzz_0, \
                             ta1_y_yyzz_xzzzz_0, \
                             ta1_y_yyzz_yyyyy_0, \
                             ta1_y_yyzz_yyyyz_0, \
                             ta1_y_yyzz_yyyzz_0, \
                             ta1_y_yyzz_yyzzz_0, \
                             ta1_y_yyzz_yzzzz_0, \
                             ta1_y_yyzz_zzzzz_0, \
                             ta1_y_yzz_xxxxz_0,  \
                             ta1_y_yzz_xxxxz_1,  \
                             ta1_y_yzz_xxxzz_0,  \
                             ta1_y_yzz_xxxzz_1,  \
                             ta1_y_yzz_xxzzz_0,  \
                             ta1_y_yzz_xxzzz_1,  \
                             ta1_y_yzz_xzzzz_0,  \
                             ta1_y_yzz_xzzzz_1,  \
                             ta1_y_yzz_zzzzz_0,  \
                             ta1_y_yzz_zzzzz_1,  \
                             ta1_y_zz_xxxxz_0,   \
                             ta1_y_zz_xxxxz_1,   \
                             ta1_y_zz_xxxzz_0,   \
                             ta1_y_zz_xxxzz_1,   \
                             ta1_y_zz_xxzzz_0,   \
                             ta1_y_zz_xxzzz_1,   \
                             ta1_y_zz_xzzzz_0,   \
                             ta1_y_zz_xzzzz_1,   \
                             ta1_y_zz_zzzzz_0,   \
                             ta1_y_zz_zzzzz_1,   \
                             ta_yzz_xxxxz_1,     \
                             ta_yzz_xxxzz_1,     \
                             ta_yzz_xxzzz_1,     \
                             ta_yzz_xzzzz_1,     \
                             ta_yzz_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yyzz_xxxxx_0[i] =
            ta1_y_yy_xxxxx_0[i] * fe_0 - ta1_y_yy_xxxxx_1[i] * fe_0 + ta1_y_yyz_xxxxx_0[i] * pa_z[i] - ta1_y_yyz_xxxxx_1[i] * pc_z[i];

        ta1_y_yyzz_xxxxy_0[i] =
            ta1_y_yy_xxxxy_0[i] * fe_0 - ta1_y_yy_xxxxy_1[i] * fe_0 + ta1_y_yyz_xxxxy_0[i] * pa_z[i] - ta1_y_yyz_xxxxy_1[i] * pc_z[i];

        ta1_y_yyzz_xxxxz_0[i] = ta1_y_zz_xxxxz_0[i] * fe_0 - ta1_y_zz_xxxxz_1[i] * fe_0 + ta_yzz_xxxxz_1[i] + ta1_y_yzz_xxxxz_0[i] * pa_y[i] -
                                ta1_y_yzz_xxxxz_1[i] * pc_y[i];

        ta1_y_yyzz_xxxyy_0[i] =
            ta1_y_yy_xxxyy_0[i] * fe_0 - ta1_y_yy_xxxyy_1[i] * fe_0 + ta1_y_yyz_xxxyy_0[i] * pa_z[i] - ta1_y_yyz_xxxyy_1[i] * pc_z[i];

        ta1_y_yyzz_xxxyz_0[i] = ta1_y_yy_xxxyz_0[i] * fe_0 - ta1_y_yy_xxxyz_1[i] * fe_0 + ta1_y_yyz_xxxy_0[i] * fe_0 - ta1_y_yyz_xxxy_1[i] * fe_0 +
                                ta1_y_yyz_xxxyz_0[i] * pa_z[i] - ta1_y_yyz_xxxyz_1[i] * pc_z[i];

        ta1_y_yyzz_xxxzz_0[i] = ta1_y_zz_xxxzz_0[i] * fe_0 - ta1_y_zz_xxxzz_1[i] * fe_0 + ta_yzz_xxxzz_1[i] + ta1_y_yzz_xxxzz_0[i] * pa_y[i] -
                                ta1_y_yzz_xxxzz_1[i] * pc_y[i];

        ta1_y_yyzz_xxyyy_0[i] =
            ta1_y_yy_xxyyy_0[i] * fe_0 - ta1_y_yy_xxyyy_1[i] * fe_0 + ta1_y_yyz_xxyyy_0[i] * pa_z[i] - ta1_y_yyz_xxyyy_1[i] * pc_z[i];

        ta1_y_yyzz_xxyyz_0[i] = ta1_y_yy_xxyyz_0[i] * fe_0 - ta1_y_yy_xxyyz_1[i] * fe_0 + ta1_y_yyz_xxyy_0[i] * fe_0 - ta1_y_yyz_xxyy_1[i] * fe_0 +
                                ta1_y_yyz_xxyyz_0[i] * pa_z[i] - ta1_y_yyz_xxyyz_1[i] * pc_z[i];

        ta1_y_yyzz_xxyzz_0[i] = ta1_y_yy_xxyzz_0[i] * fe_0 - ta1_y_yy_xxyzz_1[i] * fe_0 + 2.0 * ta1_y_yyz_xxyz_0[i] * fe_0 -
                                2.0 * ta1_y_yyz_xxyz_1[i] * fe_0 + ta1_y_yyz_xxyzz_0[i] * pa_z[i] - ta1_y_yyz_xxyzz_1[i] * pc_z[i];

        ta1_y_yyzz_xxzzz_0[i] = ta1_y_zz_xxzzz_0[i] * fe_0 - ta1_y_zz_xxzzz_1[i] * fe_0 + ta_yzz_xxzzz_1[i] + ta1_y_yzz_xxzzz_0[i] * pa_y[i] -
                                ta1_y_yzz_xxzzz_1[i] * pc_y[i];

        ta1_y_yyzz_xyyyy_0[i] =
            ta1_y_yy_xyyyy_0[i] * fe_0 - ta1_y_yy_xyyyy_1[i] * fe_0 + ta1_y_yyz_xyyyy_0[i] * pa_z[i] - ta1_y_yyz_xyyyy_1[i] * pc_z[i];

        ta1_y_yyzz_xyyyz_0[i] = ta1_y_yy_xyyyz_0[i] * fe_0 - ta1_y_yy_xyyyz_1[i] * fe_0 + ta1_y_yyz_xyyy_0[i] * fe_0 - ta1_y_yyz_xyyy_1[i] * fe_0 +
                                ta1_y_yyz_xyyyz_0[i] * pa_z[i] - ta1_y_yyz_xyyyz_1[i] * pc_z[i];

        ta1_y_yyzz_xyyzz_0[i] = ta1_y_yy_xyyzz_0[i] * fe_0 - ta1_y_yy_xyyzz_1[i] * fe_0 + 2.0 * ta1_y_yyz_xyyz_0[i] * fe_0 -
                                2.0 * ta1_y_yyz_xyyz_1[i] * fe_0 + ta1_y_yyz_xyyzz_0[i] * pa_z[i] - ta1_y_yyz_xyyzz_1[i] * pc_z[i];

        ta1_y_yyzz_xyzzz_0[i] = ta1_y_yy_xyzzz_0[i] * fe_0 - ta1_y_yy_xyzzz_1[i] * fe_0 + 3.0 * ta1_y_yyz_xyzz_0[i] * fe_0 -
                                3.0 * ta1_y_yyz_xyzz_1[i] * fe_0 + ta1_y_yyz_xyzzz_0[i] * pa_z[i] - ta1_y_yyz_xyzzz_1[i] * pc_z[i];

        ta1_y_yyzz_xzzzz_0[i] = ta1_y_zz_xzzzz_0[i] * fe_0 - ta1_y_zz_xzzzz_1[i] * fe_0 + ta_yzz_xzzzz_1[i] + ta1_y_yzz_xzzzz_0[i] * pa_y[i] -
                                ta1_y_yzz_xzzzz_1[i] * pc_y[i];

        ta1_y_yyzz_yyyyy_0[i] =
            ta1_y_yy_yyyyy_0[i] * fe_0 - ta1_y_yy_yyyyy_1[i] * fe_0 + ta1_y_yyz_yyyyy_0[i] * pa_z[i] - ta1_y_yyz_yyyyy_1[i] * pc_z[i];

        ta1_y_yyzz_yyyyz_0[i] = ta1_y_yy_yyyyz_0[i] * fe_0 - ta1_y_yy_yyyyz_1[i] * fe_0 + ta1_y_yyz_yyyy_0[i] * fe_0 - ta1_y_yyz_yyyy_1[i] * fe_0 +
                                ta1_y_yyz_yyyyz_0[i] * pa_z[i] - ta1_y_yyz_yyyyz_1[i] * pc_z[i];

        ta1_y_yyzz_yyyzz_0[i] = ta1_y_yy_yyyzz_0[i] * fe_0 - ta1_y_yy_yyyzz_1[i] * fe_0 + 2.0 * ta1_y_yyz_yyyz_0[i] * fe_0 -
                                2.0 * ta1_y_yyz_yyyz_1[i] * fe_0 + ta1_y_yyz_yyyzz_0[i] * pa_z[i] - ta1_y_yyz_yyyzz_1[i] * pc_z[i];

        ta1_y_yyzz_yyzzz_0[i] = ta1_y_yy_yyzzz_0[i] * fe_0 - ta1_y_yy_yyzzz_1[i] * fe_0 + 3.0 * ta1_y_yyz_yyzz_0[i] * fe_0 -
                                3.0 * ta1_y_yyz_yyzz_1[i] * fe_0 + ta1_y_yyz_yyzzz_0[i] * pa_z[i] - ta1_y_yyz_yyzzz_1[i] * pc_z[i];

        ta1_y_yyzz_yzzzz_0[i] = ta1_y_yy_yzzzz_0[i] * fe_0 - ta1_y_yy_yzzzz_1[i] * fe_0 + 4.0 * ta1_y_yyz_yzzz_0[i] * fe_0 -
                                4.0 * ta1_y_yyz_yzzz_1[i] * fe_0 + ta1_y_yyz_yzzzz_0[i] * pa_z[i] - ta1_y_yyz_yzzzz_1[i] * pc_z[i];

        ta1_y_yyzz_zzzzz_0[i] = ta1_y_zz_zzzzz_0[i] * fe_0 - ta1_y_zz_zzzzz_1[i] * fe_0 + ta_yzz_zzzzz_1[i] + ta1_y_yzz_zzzzz_0[i] * pa_y[i] -
                                ta1_y_yzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 588-609 components of targeted buffer : GH

    auto ta1_y_yzzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 588);

    auto ta1_y_yzzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 589);

    auto ta1_y_yzzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 590);

    auto ta1_y_yzzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 591);

    auto ta1_y_yzzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 592);

    auto ta1_y_yzzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 593);

    auto ta1_y_yzzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 594);

    auto ta1_y_yzzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 595);

    auto ta1_y_yzzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 596);

    auto ta1_y_yzzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 597);

    auto ta1_y_yzzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 598);

    auto ta1_y_yzzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 599);

    auto ta1_y_yzzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 600);

    auto ta1_y_yzzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 601);

    auto ta1_y_yzzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 602);

    auto ta1_y_yzzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 603);

    auto ta1_y_yzzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 604);

    auto ta1_y_yzzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 605);

    auto ta1_y_yzzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 606);

    auto ta1_y_yzzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 607);

    auto ta1_y_yzzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 608);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_y_yz_xxxxy_0,   \
                             ta1_y_yz_xxxxy_1,   \
                             ta1_y_yz_xxxyy_0,   \
                             ta1_y_yz_xxxyy_1,   \
                             ta1_y_yz_xxyyy_0,   \
                             ta1_y_yz_xxyyy_1,   \
                             ta1_y_yz_xyyyy_0,   \
                             ta1_y_yz_xyyyy_1,   \
                             ta1_y_yz_yyyyy_0,   \
                             ta1_y_yz_yyyyy_1,   \
                             ta1_y_yzz_xxxxy_0,  \
                             ta1_y_yzz_xxxxy_1,  \
                             ta1_y_yzz_xxxyy_0,  \
                             ta1_y_yzz_xxxyy_1,  \
                             ta1_y_yzz_xxyyy_0,  \
                             ta1_y_yzz_xxyyy_1,  \
                             ta1_y_yzz_xyyyy_0,  \
                             ta1_y_yzz_xyyyy_1,  \
                             ta1_y_yzz_yyyyy_0,  \
                             ta1_y_yzz_yyyyy_1,  \
                             ta1_y_yzzz_xxxxx_0, \
                             ta1_y_yzzz_xxxxy_0, \
                             ta1_y_yzzz_xxxxz_0, \
                             ta1_y_yzzz_xxxyy_0, \
                             ta1_y_yzzz_xxxyz_0, \
                             ta1_y_yzzz_xxxzz_0, \
                             ta1_y_yzzz_xxyyy_0, \
                             ta1_y_yzzz_xxyyz_0, \
                             ta1_y_yzzz_xxyzz_0, \
                             ta1_y_yzzz_xxzzz_0, \
                             ta1_y_yzzz_xyyyy_0, \
                             ta1_y_yzzz_xyyyz_0, \
                             ta1_y_yzzz_xyyzz_0, \
                             ta1_y_yzzz_xyzzz_0, \
                             ta1_y_yzzz_xzzzz_0, \
                             ta1_y_yzzz_yyyyy_0, \
                             ta1_y_yzzz_yyyyz_0, \
                             ta1_y_yzzz_yyyzz_0, \
                             ta1_y_yzzz_yyzzz_0, \
                             ta1_y_yzzz_yzzzz_0, \
                             ta1_y_yzzz_zzzzz_0, \
                             ta1_y_zzz_xxxxx_0,  \
                             ta1_y_zzz_xxxxx_1,  \
                             ta1_y_zzz_xxxxz_0,  \
                             ta1_y_zzz_xxxxz_1,  \
                             ta1_y_zzz_xxxyz_0,  \
                             ta1_y_zzz_xxxyz_1,  \
                             ta1_y_zzz_xxxz_0,   \
                             ta1_y_zzz_xxxz_1,   \
                             ta1_y_zzz_xxxzz_0,  \
                             ta1_y_zzz_xxxzz_1,  \
                             ta1_y_zzz_xxyyz_0,  \
                             ta1_y_zzz_xxyyz_1,  \
                             ta1_y_zzz_xxyz_0,   \
                             ta1_y_zzz_xxyz_1,   \
                             ta1_y_zzz_xxyzz_0,  \
                             ta1_y_zzz_xxyzz_1,  \
                             ta1_y_zzz_xxzz_0,   \
                             ta1_y_zzz_xxzz_1,   \
                             ta1_y_zzz_xxzzz_0,  \
                             ta1_y_zzz_xxzzz_1,  \
                             ta1_y_zzz_xyyyz_0,  \
                             ta1_y_zzz_xyyyz_1,  \
                             ta1_y_zzz_xyyz_0,   \
                             ta1_y_zzz_xyyz_1,   \
                             ta1_y_zzz_xyyzz_0,  \
                             ta1_y_zzz_xyyzz_1,  \
                             ta1_y_zzz_xyzz_0,   \
                             ta1_y_zzz_xyzz_1,   \
                             ta1_y_zzz_xyzzz_0,  \
                             ta1_y_zzz_xyzzz_1,  \
                             ta1_y_zzz_xzzz_0,   \
                             ta1_y_zzz_xzzz_1,   \
                             ta1_y_zzz_xzzzz_0,  \
                             ta1_y_zzz_xzzzz_1,  \
                             ta1_y_zzz_yyyyz_0,  \
                             ta1_y_zzz_yyyyz_1,  \
                             ta1_y_zzz_yyyz_0,   \
                             ta1_y_zzz_yyyz_1,   \
                             ta1_y_zzz_yyyzz_0,  \
                             ta1_y_zzz_yyyzz_1,  \
                             ta1_y_zzz_yyzz_0,   \
                             ta1_y_zzz_yyzz_1,   \
                             ta1_y_zzz_yyzzz_0,  \
                             ta1_y_zzz_yyzzz_1,  \
                             ta1_y_zzz_yzzz_0,   \
                             ta1_y_zzz_yzzz_1,   \
                             ta1_y_zzz_yzzzz_0,  \
                             ta1_y_zzz_yzzzz_1,  \
                             ta1_y_zzz_zzzz_0,   \
                             ta1_y_zzz_zzzz_1,   \
                             ta1_y_zzz_zzzzz_0,  \
                             ta1_y_zzz_zzzzz_1,  \
                             ta_zzz_xxxxx_1,     \
                             ta_zzz_xxxxz_1,     \
                             ta_zzz_xxxyz_1,     \
                             ta_zzz_xxxzz_1,     \
                             ta_zzz_xxyyz_1,     \
                             ta_zzz_xxyzz_1,     \
                             ta_zzz_xxzzz_1,     \
                             ta_zzz_xyyyz_1,     \
                             ta_zzz_xyyzz_1,     \
                             ta_zzz_xyzzz_1,     \
                             ta_zzz_xzzzz_1,     \
                             ta_zzz_yyyyz_1,     \
                             ta_zzz_yyyzz_1,     \
                             ta_zzz_yyzzz_1,     \
                             ta_zzz_yzzzz_1,     \
                             ta_zzz_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_yzzz_xxxxx_0[i] = ta_zzz_xxxxx_1[i] + ta1_y_zzz_xxxxx_0[i] * pa_y[i] - ta1_y_zzz_xxxxx_1[i] * pc_y[i];

        ta1_y_yzzz_xxxxy_0[i] =
            2.0 * ta1_y_yz_xxxxy_0[i] * fe_0 - 2.0 * ta1_y_yz_xxxxy_1[i] * fe_0 + ta1_y_yzz_xxxxy_0[i] * pa_z[i] - ta1_y_yzz_xxxxy_1[i] * pc_z[i];

        ta1_y_yzzz_xxxxz_0[i] = ta_zzz_xxxxz_1[i] + ta1_y_zzz_xxxxz_0[i] * pa_y[i] - ta1_y_zzz_xxxxz_1[i] * pc_y[i];

        ta1_y_yzzz_xxxyy_0[i] =
            2.0 * ta1_y_yz_xxxyy_0[i] * fe_0 - 2.0 * ta1_y_yz_xxxyy_1[i] * fe_0 + ta1_y_yzz_xxxyy_0[i] * pa_z[i] - ta1_y_yzz_xxxyy_1[i] * pc_z[i];

        ta1_y_yzzz_xxxyz_0[i] = ta1_y_zzz_xxxz_0[i] * fe_0 - ta1_y_zzz_xxxz_1[i] * fe_0 + ta_zzz_xxxyz_1[i] + ta1_y_zzz_xxxyz_0[i] * pa_y[i] -
                                ta1_y_zzz_xxxyz_1[i] * pc_y[i];

        ta1_y_yzzz_xxxzz_0[i] = ta_zzz_xxxzz_1[i] + ta1_y_zzz_xxxzz_0[i] * pa_y[i] - ta1_y_zzz_xxxzz_1[i] * pc_y[i];

        ta1_y_yzzz_xxyyy_0[i] =
            2.0 * ta1_y_yz_xxyyy_0[i] * fe_0 - 2.0 * ta1_y_yz_xxyyy_1[i] * fe_0 + ta1_y_yzz_xxyyy_0[i] * pa_z[i] - ta1_y_yzz_xxyyy_1[i] * pc_z[i];

        ta1_y_yzzz_xxyyz_0[i] = 2.0 * ta1_y_zzz_xxyz_0[i] * fe_0 - 2.0 * ta1_y_zzz_xxyz_1[i] * fe_0 + ta_zzz_xxyyz_1[i] +
                                ta1_y_zzz_xxyyz_0[i] * pa_y[i] - ta1_y_zzz_xxyyz_1[i] * pc_y[i];

        ta1_y_yzzz_xxyzz_0[i] = ta1_y_zzz_xxzz_0[i] * fe_0 - ta1_y_zzz_xxzz_1[i] * fe_0 + ta_zzz_xxyzz_1[i] + ta1_y_zzz_xxyzz_0[i] * pa_y[i] -
                                ta1_y_zzz_xxyzz_1[i] * pc_y[i];

        ta1_y_yzzz_xxzzz_0[i] = ta_zzz_xxzzz_1[i] + ta1_y_zzz_xxzzz_0[i] * pa_y[i] - ta1_y_zzz_xxzzz_1[i] * pc_y[i];

        ta1_y_yzzz_xyyyy_0[i] =
            2.0 * ta1_y_yz_xyyyy_0[i] * fe_0 - 2.0 * ta1_y_yz_xyyyy_1[i] * fe_0 + ta1_y_yzz_xyyyy_0[i] * pa_z[i] - ta1_y_yzz_xyyyy_1[i] * pc_z[i];

        ta1_y_yzzz_xyyyz_0[i] = 3.0 * ta1_y_zzz_xyyz_0[i] * fe_0 - 3.0 * ta1_y_zzz_xyyz_1[i] * fe_0 + ta_zzz_xyyyz_1[i] +
                                ta1_y_zzz_xyyyz_0[i] * pa_y[i] - ta1_y_zzz_xyyyz_1[i] * pc_y[i];

        ta1_y_yzzz_xyyzz_0[i] = 2.0 * ta1_y_zzz_xyzz_0[i] * fe_0 - 2.0 * ta1_y_zzz_xyzz_1[i] * fe_0 + ta_zzz_xyyzz_1[i] +
                                ta1_y_zzz_xyyzz_0[i] * pa_y[i] - ta1_y_zzz_xyyzz_1[i] * pc_y[i];

        ta1_y_yzzz_xyzzz_0[i] = ta1_y_zzz_xzzz_0[i] * fe_0 - ta1_y_zzz_xzzz_1[i] * fe_0 + ta_zzz_xyzzz_1[i] + ta1_y_zzz_xyzzz_0[i] * pa_y[i] -
                                ta1_y_zzz_xyzzz_1[i] * pc_y[i];

        ta1_y_yzzz_xzzzz_0[i] = ta_zzz_xzzzz_1[i] + ta1_y_zzz_xzzzz_0[i] * pa_y[i] - ta1_y_zzz_xzzzz_1[i] * pc_y[i];

        ta1_y_yzzz_yyyyy_0[i] =
            2.0 * ta1_y_yz_yyyyy_0[i] * fe_0 - 2.0 * ta1_y_yz_yyyyy_1[i] * fe_0 + ta1_y_yzz_yyyyy_0[i] * pa_z[i] - ta1_y_yzz_yyyyy_1[i] * pc_z[i];

        ta1_y_yzzz_yyyyz_0[i] = 4.0 * ta1_y_zzz_yyyz_0[i] * fe_0 - 4.0 * ta1_y_zzz_yyyz_1[i] * fe_0 + ta_zzz_yyyyz_1[i] +
                                ta1_y_zzz_yyyyz_0[i] * pa_y[i] - ta1_y_zzz_yyyyz_1[i] * pc_y[i];

        ta1_y_yzzz_yyyzz_0[i] = 3.0 * ta1_y_zzz_yyzz_0[i] * fe_0 - 3.0 * ta1_y_zzz_yyzz_1[i] * fe_0 + ta_zzz_yyyzz_1[i] +
                                ta1_y_zzz_yyyzz_0[i] * pa_y[i] - ta1_y_zzz_yyyzz_1[i] * pc_y[i];

        ta1_y_yzzz_yyzzz_0[i] = 2.0 * ta1_y_zzz_yzzz_0[i] * fe_0 - 2.0 * ta1_y_zzz_yzzz_1[i] * fe_0 + ta_zzz_yyzzz_1[i] +
                                ta1_y_zzz_yyzzz_0[i] * pa_y[i] - ta1_y_zzz_yyzzz_1[i] * pc_y[i];

        ta1_y_yzzz_yzzzz_0[i] = ta1_y_zzz_zzzz_0[i] * fe_0 - ta1_y_zzz_zzzz_1[i] * fe_0 + ta_zzz_yzzzz_1[i] + ta1_y_zzz_yzzzz_0[i] * pa_y[i] -
                                ta1_y_zzz_yzzzz_1[i] * pc_y[i];

        ta1_y_yzzz_zzzzz_0[i] = ta_zzz_zzzzz_1[i] + ta1_y_zzz_zzzzz_0[i] * pa_y[i] - ta1_y_zzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 609-630 components of targeted buffer : GH

    auto ta1_y_zzzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 609);

    auto ta1_y_zzzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 610);

    auto ta1_y_zzzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 611);

    auto ta1_y_zzzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 612);

    auto ta1_y_zzzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 613);

    auto ta1_y_zzzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 614);

    auto ta1_y_zzzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 615);

    auto ta1_y_zzzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 616);

    auto ta1_y_zzzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 617);

    auto ta1_y_zzzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 618);

    auto ta1_y_zzzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 619);

    auto ta1_y_zzzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 620);

    auto ta1_y_zzzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 621);

    auto ta1_y_zzzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 622);

    auto ta1_y_zzzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 623);

    auto ta1_y_zzzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 624);

    auto ta1_y_zzzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 625);

    auto ta1_y_zzzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 626);

    auto ta1_y_zzzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 627);

    auto ta1_y_zzzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 628);

    auto ta1_y_zzzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 629);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_y_zz_xxxxx_0,   \
                             ta1_y_zz_xxxxx_1,   \
                             ta1_y_zz_xxxxy_0,   \
                             ta1_y_zz_xxxxy_1,   \
                             ta1_y_zz_xxxxz_0,   \
                             ta1_y_zz_xxxxz_1,   \
                             ta1_y_zz_xxxyy_0,   \
                             ta1_y_zz_xxxyy_1,   \
                             ta1_y_zz_xxxyz_0,   \
                             ta1_y_zz_xxxyz_1,   \
                             ta1_y_zz_xxxzz_0,   \
                             ta1_y_zz_xxxzz_1,   \
                             ta1_y_zz_xxyyy_0,   \
                             ta1_y_zz_xxyyy_1,   \
                             ta1_y_zz_xxyyz_0,   \
                             ta1_y_zz_xxyyz_1,   \
                             ta1_y_zz_xxyzz_0,   \
                             ta1_y_zz_xxyzz_1,   \
                             ta1_y_zz_xxzzz_0,   \
                             ta1_y_zz_xxzzz_1,   \
                             ta1_y_zz_xyyyy_0,   \
                             ta1_y_zz_xyyyy_1,   \
                             ta1_y_zz_xyyyz_0,   \
                             ta1_y_zz_xyyyz_1,   \
                             ta1_y_zz_xyyzz_0,   \
                             ta1_y_zz_xyyzz_1,   \
                             ta1_y_zz_xyzzz_0,   \
                             ta1_y_zz_xyzzz_1,   \
                             ta1_y_zz_xzzzz_0,   \
                             ta1_y_zz_xzzzz_1,   \
                             ta1_y_zz_yyyyy_0,   \
                             ta1_y_zz_yyyyy_1,   \
                             ta1_y_zz_yyyyz_0,   \
                             ta1_y_zz_yyyyz_1,   \
                             ta1_y_zz_yyyzz_0,   \
                             ta1_y_zz_yyyzz_1,   \
                             ta1_y_zz_yyzzz_0,   \
                             ta1_y_zz_yyzzz_1,   \
                             ta1_y_zz_yzzzz_0,   \
                             ta1_y_zz_yzzzz_1,   \
                             ta1_y_zz_zzzzz_0,   \
                             ta1_y_zz_zzzzz_1,   \
                             ta1_y_zzz_xxxx_0,   \
                             ta1_y_zzz_xxxx_1,   \
                             ta1_y_zzz_xxxxx_0,  \
                             ta1_y_zzz_xxxxx_1,  \
                             ta1_y_zzz_xxxxy_0,  \
                             ta1_y_zzz_xxxxy_1,  \
                             ta1_y_zzz_xxxxz_0,  \
                             ta1_y_zzz_xxxxz_1,  \
                             ta1_y_zzz_xxxy_0,   \
                             ta1_y_zzz_xxxy_1,   \
                             ta1_y_zzz_xxxyy_0,  \
                             ta1_y_zzz_xxxyy_1,  \
                             ta1_y_zzz_xxxyz_0,  \
                             ta1_y_zzz_xxxyz_1,  \
                             ta1_y_zzz_xxxz_0,   \
                             ta1_y_zzz_xxxz_1,   \
                             ta1_y_zzz_xxxzz_0,  \
                             ta1_y_zzz_xxxzz_1,  \
                             ta1_y_zzz_xxyy_0,   \
                             ta1_y_zzz_xxyy_1,   \
                             ta1_y_zzz_xxyyy_0,  \
                             ta1_y_zzz_xxyyy_1,  \
                             ta1_y_zzz_xxyyz_0,  \
                             ta1_y_zzz_xxyyz_1,  \
                             ta1_y_zzz_xxyz_0,   \
                             ta1_y_zzz_xxyz_1,   \
                             ta1_y_zzz_xxyzz_0,  \
                             ta1_y_zzz_xxyzz_1,  \
                             ta1_y_zzz_xxzz_0,   \
                             ta1_y_zzz_xxzz_1,   \
                             ta1_y_zzz_xxzzz_0,  \
                             ta1_y_zzz_xxzzz_1,  \
                             ta1_y_zzz_xyyy_0,   \
                             ta1_y_zzz_xyyy_1,   \
                             ta1_y_zzz_xyyyy_0,  \
                             ta1_y_zzz_xyyyy_1,  \
                             ta1_y_zzz_xyyyz_0,  \
                             ta1_y_zzz_xyyyz_1,  \
                             ta1_y_zzz_xyyz_0,   \
                             ta1_y_zzz_xyyz_1,   \
                             ta1_y_zzz_xyyzz_0,  \
                             ta1_y_zzz_xyyzz_1,  \
                             ta1_y_zzz_xyzz_0,   \
                             ta1_y_zzz_xyzz_1,   \
                             ta1_y_zzz_xyzzz_0,  \
                             ta1_y_zzz_xyzzz_1,  \
                             ta1_y_zzz_xzzz_0,   \
                             ta1_y_zzz_xzzz_1,   \
                             ta1_y_zzz_xzzzz_0,  \
                             ta1_y_zzz_xzzzz_1,  \
                             ta1_y_zzz_yyyy_0,   \
                             ta1_y_zzz_yyyy_1,   \
                             ta1_y_zzz_yyyyy_0,  \
                             ta1_y_zzz_yyyyy_1,  \
                             ta1_y_zzz_yyyyz_0,  \
                             ta1_y_zzz_yyyyz_1,  \
                             ta1_y_zzz_yyyz_0,   \
                             ta1_y_zzz_yyyz_1,   \
                             ta1_y_zzz_yyyzz_0,  \
                             ta1_y_zzz_yyyzz_1,  \
                             ta1_y_zzz_yyzz_0,   \
                             ta1_y_zzz_yyzz_1,   \
                             ta1_y_zzz_yyzzz_0,  \
                             ta1_y_zzz_yyzzz_1,  \
                             ta1_y_zzz_yzzz_0,   \
                             ta1_y_zzz_yzzz_1,   \
                             ta1_y_zzz_yzzzz_0,  \
                             ta1_y_zzz_yzzzz_1,  \
                             ta1_y_zzz_zzzz_0,   \
                             ta1_y_zzz_zzzz_1,   \
                             ta1_y_zzz_zzzzz_0,  \
                             ta1_y_zzz_zzzzz_1,  \
                             ta1_y_zzzz_xxxxx_0, \
                             ta1_y_zzzz_xxxxy_0, \
                             ta1_y_zzzz_xxxxz_0, \
                             ta1_y_zzzz_xxxyy_0, \
                             ta1_y_zzzz_xxxyz_0, \
                             ta1_y_zzzz_xxxzz_0, \
                             ta1_y_zzzz_xxyyy_0, \
                             ta1_y_zzzz_xxyyz_0, \
                             ta1_y_zzzz_xxyzz_0, \
                             ta1_y_zzzz_xxzzz_0, \
                             ta1_y_zzzz_xyyyy_0, \
                             ta1_y_zzzz_xyyyz_0, \
                             ta1_y_zzzz_xyyzz_0, \
                             ta1_y_zzzz_xyzzz_0, \
                             ta1_y_zzzz_xzzzz_0, \
                             ta1_y_zzzz_yyyyy_0, \
                             ta1_y_zzzz_yyyyz_0, \
                             ta1_y_zzzz_yyyzz_0, \
                             ta1_y_zzzz_yyzzz_0, \
                             ta1_y_zzzz_yzzzz_0, \
                             ta1_y_zzzz_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_y_zzzz_xxxxx_0[i] =
            3.0 * ta1_y_zz_xxxxx_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxxx_1[i] * fe_0 + ta1_y_zzz_xxxxx_0[i] * pa_z[i] - ta1_y_zzz_xxxxx_1[i] * pc_z[i];

        ta1_y_zzzz_xxxxy_0[i] =
            3.0 * ta1_y_zz_xxxxy_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxxy_1[i] * fe_0 + ta1_y_zzz_xxxxy_0[i] * pa_z[i] - ta1_y_zzz_xxxxy_1[i] * pc_z[i];

        ta1_y_zzzz_xxxxz_0[i] = 3.0 * ta1_y_zz_xxxxz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxxz_1[i] * fe_0 + ta1_y_zzz_xxxx_0[i] * fe_0 -
                                ta1_y_zzz_xxxx_1[i] * fe_0 + ta1_y_zzz_xxxxz_0[i] * pa_z[i] - ta1_y_zzz_xxxxz_1[i] * pc_z[i];

        ta1_y_zzzz_xxxyy_0[i] =
            3.0 * ta1_y_zz_xxxyy_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxyy_1[i] * fe_0 + ta1_y_zzz_xxxyy_0[i] * pa_z[i] - ta1_y_zzz_xxxyy_1[i] * pc_z[i];

        ta1_y_zzzz_xxxyz_0[i] = 3.0 * ta1_y_zz_xxxyz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxyz_1[i] * fe_0 + ta1_y_zzz_xxxy_0[i] * fe_0 -
                                ta1_y_zzz_xxxy_1[i] * fe_0 + ta1_y_zzz_xxxyz_0[i] * pa_z[i] - ta1_y_zzz_xxxyz_1[i] * pc_z[i];

        ta1_y_zzzz_xxxzz_0[i] = 3.0 * ta1_y_zz_xxxzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxxzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_xxxz_0[i] * fe_0 -
                                2.0 * ta1_y_zzz_xxxz_1[i] * fe_0 + ta1_y_zzz_xxxzz_0[i] * pa_z[i] - ta1_y_zzz_xxxzz_1[i] * pc_z[i];

        ta1_y_zzzz_xxyyy_0[i] =
            3.0 * ta1_y_zz_xxyyy_0[i] * fe_0 - 3.0 * ta1_y_zz_xxyyy_1[i] * fe_0 + ta1_y_zzz_xxyyy_0[i] * pa_z[i] - ta1_y_zzz_xxyyy_1[i] * pc_z[i];

        ta1_y_zzzz_xxyyz_0[i] = 3.0 * ta1_y_zz_xxyyz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxyyz_1[i] * fe_0 + ta1_y_zzz_xxyy_0[i] * fe_0 -
                                ta1_y_zzz_xxyy_1[i] * fe_0 + ta1_y_zzz_xxyyz_0[i] * pa_z[i] - ta1_y_zzz_xxyyz_1[i] * pc_z[i];

        ta1_y_zzzz_xxyzz_0[i] = 3.0 * ta1_y_zz_xxyzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxyzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_xxyz_0[i] * fe_0 -
                                2.0 * ta1_y_zzz_xxyz_1[i] * fe_0 + ta1_y_zzz_xxyzz_0[i] * pa_z[i] - ta1_y_zzz_xxyzz_1[i] * pc_z[i];

        ta1_y_zzzz_xxzzz_0[i] = 3.0 * ta1_y_zz_xxzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xxzzz_1[i] * fe_0 + 3.0 * ta1_y_zzz_xxzz_0[i] * fe_0 -
                                3.0 * ta1_y_zzz_xxzz_1[i] * fe_0 + ta1_y_zzz_xxzzz_0[i] * pa_z[i] - ta1_y_zzz_xxzzz_1[i] * pc_z[i];

        ta1_y_zzzz_xyyyy_0[i] =
            3.0 * ta1_y_zz_xyyyy_0[i] * fe_0 - 3.0 * ta1_y_zz_xyyyy_1[i] * fe_0 + ta1_y_zzz_xyyyy_0[i] * pa_z[i] - ta1_y_zzz_xyyyy_1[i] * pc_z[i];

        ta1_y_zzzz_xyyyz_0[i] = 3.0 * ta1_y_zz_xyyyz_0[i] * fe_0 - 3.0 * ta1_y_zz_xyyyz_1[i] * fe_0 + ta1_y_zzz_xyyy_0[i] * fe_0 -
                                ta1_y_zzz_xyyy_1[i] * fe_0 + ta1_y_zzz_xyyyz_0[i] * pa_z[i] - ta1_y_zzz_xyyyz_1[i] * pc_z[i];

        ta1_y_zzzz_xyyzz_0[i] = 3.0 * ta1_y_zz_xyyzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xyyzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_xyyz_0[i] * fe_0 -
                                2.0 * ta1_y_zzz_xyyz_1[i] * fe_0 + ta1_y_zzz_xyyzz_0[i] * pa_z[i] - ta1_y_zzz_xyyzz_1[i] * pc_z[i];

        ta1_y_zzzz_xyzzz_0[i] = 3.0 * ta1_y_zz_xyzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xyzzz_1[i] * fe_0 + 3.0 * ta1_y_zzz_xyzz_0[i] * fe_0 -
                                3.0 * ta1_y_zzz_xyzz_1[i] * fe_0 + ta1_y_zzz_xyzzz_0[i] * pa_z[i] - ta1_y_zzz_xyzzz_1[i] * pc_z[i];

        ta1_y_zzzz_xzzzz_0[i] = 3.0 * ta1_y_zz_xzzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_xzzzz_1[i] * fe_0 + 4.0 * ta1_y_zzz_xzzz_0[i] * fe_0 -
                                4.0 * ta1_y_zzz_xzzz_1[i] * fe_0 + ta1_y_zzz_xzzzz_0[i] * pa_z[i] - ta1_y_zzz_xzzzz_1[i] * pc_z[i];

        ta1_y_zzzz_yyyyy_0[i] =
            3.0 * ta1_y_zz_yyyyy_0[i] * fe_0 - 3.0 * ta1_y_zz_yyyyy_1[i] * fe_0 + ta1_y_zzz_yyyyy_0[i] * pa_z[i] - ta1_y_zzz_yyyyy_1[i] * pc_z[i];

        ta1_y_zzzz_yyyyz_0[i] = 3.0 * ta1_y_zz_yyyyz_0[i] * fe_0 - 3.0 * ta1_y_zz_yyyyz_1[i] * fe_0 + ta1_y_zzz_yyyy_0[i] * fe_0 -
                                ta1_y_zzz_yyyy_1[i] * fe_0 + ta1_y_zzz_yyyyz_0[i] * pa_z[i] - ta1_y_zzz_yyyyz_1[i] * pc_z[i];

        ta1_y_zzzz_yyyzz_0[i] = 3.0 * ta1_y_zz_yyyzz_0[i] * fe_0 - 3.0 * ta1_y_zz_yyyzz_1[i] * fe_0 + 2.0 * ta1_y_zzz_yyyz_0[i] * fe_0 -
                                2.0 * ta1_y_zzz_yyyz_1[i] * fe_0 + ta1_y_zzz_yyyzz_0[i] * pa_z[i] - ta1_y_zzz_yyyzz_1[i] * pc_z[i];

        ta1_y_zzzz_yyzzz_0[i] = 3.0 * ta1_y_zz_yyzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_yyzzz_1[i] * fe_0 + 3.0 * ta1_y_zzz_yyzz_0[i] * fe_0 -
                                3.0 * ta1_y_zzz_yyzz_1[i] * fe_0 + ta1_y_zzz_yyzzz_0[i] * pa_z[i] - ta1_y_zzz_yyzzz_1[i] * pc_z[i];

        ta1_y_zzzz_yzzzz_0[i] = 3.0 * ta1_y_zz_yzzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_yzzzz_1[i] * fe_0 + 4.0 * ta1_y_zzz_yzzz_0[i] * fe_0 -
                                4.0 * ta1_y_zzz_yzzz_1[i] * fe_0 + ta1_y_zzz_yzzzz_0[i] * pa_z[i] - ta1_y_zzz_yzzzz_1[i] * pc_z[i];

        ta1_y_zzzz_zzzzz_0[i] = 3.0 * ta1_y_zz_zzzzz_0[i] * fe_0 - 3.0 * ta1_y_zz_zzzzz_1[i] * fe_0 + 5.0 * ta1_y_zzz_zzzz_0[i] * fe_0 -
                                5.0 * ta1_y_zzz_zzzz_1[i] * fe_0 + ta1_y_zzz_zzzzz_0[i] * pa_z[i] - ta1_y_zzz_zzzzz_1[i] * pc_z[i];
    }

    // Set up 630-651 components of targeted buffer : GH

    auto ta1_z_xxxx_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 630);

    auto ta1_z_xxxx_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 631);

    auto ta1_z_xxxx_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 632);

    auto ta1_z_xxxx_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 633);

    auto ta1_z_xxxx_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 634);

    auto ta1_z_xxxx_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 635);

    auto ta1_z_xxxx_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 636);

    auto ta1_z_xxxx_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 637);

    auto ta1_z_xxxx_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 638);

    auto ta1_z_xxxx_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 639);

    auto ta1_z_xxxx_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 640);

    auto ta1_z_xxxx_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 641);

    auto ta1_z_xxxx_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 642);

    auto ta1_z_xxxx_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 643);

    auto ta1_z_xxxx_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 644);

    auto ta1_z_xxxx_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 645);

    auto ta1_z_xxxx_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 646);

    auto ta1_z_xxxx_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 647);

    auto ta1_z_xxxx_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 648);

    auto ta1_z_xxxx_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 649);

    auto ta1_z_xxxx_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 650);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_xx_xxxxx_0,   \
                             ta1_z_xx_xxxxx_1,   \
                             ta1_z_xx_xxxxy_0,   \
                             ta1_z_xx_xxxxy_1,   \
                             ta1_z_xx_xxxxz_0,   \
                             ta1_z_xx_xxxxz_1,   \
                             ta1_z_xx_xxxyy_0,   \
                             ta1_z_xx_xxxyy_1,   \
                             ta1_z_xx_xxxyz_0,   \
                             ta1_z_xx_xxxyz_1,   \
                             ta1_z_xx_xxxzz_0,   \
                             ta1_z_xx_xxxzz_1,   \
                             ta1_z_xx_xxyyy_0,   \
                             ta1_z_xx_xxyyy_1,   \
                             ta1_z_xx_xxyyz_0,   \
                             ta1_z_xx_xxyyz_1,   \
                             ta1_z_xx_xxyzz_0,   \
                             ta1_z_xx_xxyzz_1,   \
                             ta1_z_xx_xxzzz_0,   \
                             ta1_z_xx_xxzzz_1,   \
                             ta1_z_xx_xyyyy_0,   \
                             ta1_z_xx_xyyyy_1,   \
                             ta1_z_xx_xyyyz_0,   \
                             ta1_z_xx_xyyyz_1,   \
                             ta1_z_xx_xyyzz_0,   \
                             ta1_z_xx_xyyzz_1,   \
                             ta1_z_xx_xyzzz_0,   \
                             ta1_z_xx_xyzzz_1,   \
                             ta1_z_xx_xzzzz_0,   \
                             ta1_z_xx_xzzzz_1,   \
                             ta1_z_xx_yyyyy_0,   \
                             ta1_z_xx_yyyyy_1,   \
                             ta1_z_xx_yyyyz_0,   \
                             ta1_z_xx_yyyyz_1,   \
                             ta1_z_xx_yyyzz_0,   \
                             ta1_z_xx_yyyzz_1,   \
                             ta1_z_xx_yyzzz_0,   \
                             ta1_z_xx_yyzzz_1,   \
                             ta1_z_xx_yzzzz_0,   \
                             ta1_z_xx_yzzzz_1,   \
                             ta1_z_xx_zzzzz_0,   \
                             ta1_z_xx_zzzzz_1,   \
                             ta1_z_xxx_xxxx_0,   \
                             ta1_z_xxx_xxxx_1,   \
                             ta1_z_xxx_xxxxx_0,  \
                             ta1_z_xxx_xxxxx_1,  \
                             ta1_z_xxx_xxxxy_0,  \
                             ta1_z_xxx_xxxxy_1,  \
                             ta1_z_xxx_xxxxz_0,  \
                             ta1_z_xxx_xxxxz_1,  \
                             ta1_z_xxx_xxxy_0,   \
                             ta1_z_xxx_xxxy_1,   \
                             ta1_z_xxx_xxxyy_0,  \
                             ta1_z_xxx_xxxyy_1,  \
                             ta1_z_xxx_xxxyz_0,  \
                             ta1_z_xxx_xxxyz_1,  \
                             ta1_z_xxx_xxxz_0,   \
                             ta1_z_xxx_xxxz_1,   \
                             ta1_z_xxx_xxxzz_0,  \
                             ta1_z_xxx_xxxzz_1,  \
                             ta1_z_xxx_xxyy_0,   \
                             ta1_z_xxx_xxyy_1,   \
                             ta1_z_xxx_xxyyy_0,  \
                             ta1_z_xxx_xxyyy_1,  \
                             ta1_z_xxx_xxyyz_0,  \
                             ta1_z_xxx_xxyyz_1,  \
                             ta1_z_xxx_xxyz_0,   \
                             ta1_z_xxx_xxyz_1,   \
                             ta1_z_xxx_xxyzz_0,  \
                             ta1_z_xxx_xxyzz_1,  \
                             ta1_z_xxx_xxzz_0,   \
                             ta1_z_xxx_xxzz_1,   \
                             ta1_z_xxx_xxzzz_0,  \
                             ta1_z_xxx_xxzzz_1,  \
                             ta1_z_xxx_xyyy_0,   \
                             ta1_z_xxx_xyyy_1,   \
                             ta1_z_xxx_xyyyy_0,  \
                             ta1_z_xxx_xyyyy_1,  \
                             ta1_z_xxx_xyyyz_0,  \
                             ta1_z_xxx_xyyyz_1,  \
                             ta1_z_xxx_xyyz_0,   \
                             ta1_z_xxx_xyyz_1,   \
                             ta1_z_xxx_xyyzz_0,  \
                             ta1_z_xxx_xyyzz_1,  \
                             ta1_z_xxx_xyzz_0,   \
                             ta1_z_xxx_xyzz_1,   \
                             ta1_z_xxx_xyzzz_0,  \
                             ta1_z_xxx_xyzzz_1,  \
                             ta1_z_xxx_xzzz_0,   \
                             ta1_z_xxx_xzzz_1,   \
                             ta1_z_xxx_xzzzz_0,  \
                             ta1_z_xxx_xzzzz_1,  \
                             ta1_z_xxx_yyyy_0,   \
                             ta1_z_xxx_yyyy_1,   \
                             ta1_z_xxx_yyyyy_0,  \
                             ta1_z_xxx_yyyyy_1,  \
                             ta1_z_xxx_yyyyz_0,  \
                             ta1_z_xxx_yyyyz_1,  \
                             ta1_z_xxx_yyyz_0,   \
                             ta1_z_xxx_yyyz_1,   \
                             ta1_z_xxx_yyyzz_0,  \
                             ta1_z_xxx_yyyzz_1,  \
                             ta1_z_xxx_yyzz_0,   \
                             ta1_z_xxx_yyzz_1,   \
                             ta1_z_xxx_yyzzz_0,  \
                             ta1_z_xxx_yyzzz_1,  \
                             ta1_z_xxx_yzzz_0,   \
                             ta1_z_xxx_yzzz_1,   \
                             ta1_z_xxx_yzzzz_0,  \
                             ta1_z_xxx_yzzzz_1,  \
                             ta1_z_xxx_zzzz_0,   \
                             ta1_z_xxx_zzzz_1,   \
                             ta1_z_xxx_zzzzz_0,  \
                             ta1_z_xxx_zzzzz_1,  \
                             ta1_z_xxxx_xxxxx_0, \
                             ta1_z_xxxx_xxxxy_0, \
                             ta1_z_xxxx_xxxxz_0, \
                             ta1_z_xxxx_xxxyy_0, \
                             ta1_z_xxxx_xxxyz_0, \
                             ta1_z_xxxx_xxxzz_0, \
                             ta1_z_xxxx_xxyyy_0, \
                             ta1_z_xxxx_xxyyz_0, \
                             ta1_z_xxxx_xxyzz_0, \
                             ta1_z_xxxx_xxzzz_0, \
                             ta1_z_xxxx_xyyyy_0, \
                             ta1_z_xxxx_xyyyz_0, \
                             ta1_z_xxxx_xyyzz_0, \
                             ta1_z_xxxx_xyzzz_0, \
                             ta1_z_xxxx_xzzzz_0, \
                             ta1_z_xxxx_yyyyy_0, \
                             ta1_z_xxxx_yyyyz_0, \
                             ta1_z_xxxx_yyyzz_0, \
                             ta1_z_xxxx_yyzzz_0, \
                             ta1_z_xxxx_yzzzz_0, \
                             ta1_z_xxxx_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxx_xxxxx_0[i] = 3.0 * ta1_z_xx_xxxxx_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxxx_1[i] * fe_0 + 5.0 * ta1_z_xxx_xxxx_0[i] * fe_0 -
                                5.0 * ta1_z_xxx_xxxx_1[i] * fe_0 + ta1_z_xxx_xxxxx_0[i] * pa_x[i] - ta1_z_xxx_xxxxx_1[i] * pc_x[i];

        ta1_z_xxxx_xxxxy_0[i] = 3.0 * ta1_z_xx_xxxxy_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxxy_1[i] * fe_0 + 4.0 * ta1_z_xxx_xxxy_0[i] * fe_0 -
                                4.0 * ta1_z_xxx_xxxy_1[i] * fe_0 + ta1_z_xxx_xxxxy_0[i] * pa_x[i] - ta1_z_xxx_xxxxy_1[i] * pc_x[i];

        ta1_z_xxxx_xxxxz_0[i] = 3.0 * ta1_z_xx_xxxxz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxxz_1[i] * fe_0 + 4.0 * ta1_z_xxx_xxxz_0[i] * fe_0 -
                                4.0 * ta1_z_xxx_xxxz_1[i] * fe_0 + ta1_z_xxx_xxxxz_0[i] * pa_x[i] - ta1_z_xxx_xxxxz_1[i] * pc_x[i];

        ta1_z_xxxx_xxxyy_0[i] = 3.0 * ta1_z_xx_xxxyy_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxyy_1[i] * fe_0 + 3.0 * ta1_z_xxx_xxyy_0[i] * fe_0 -
                                3.0 * ta1_z_xxx_xxyy_1[i] * fe_0 + ta1_z_xxx_xxxyy_0[i] * pa_x[i] - ta1_z_xxx_xxxyy_1[i] * pc_x[i];

        ta1_z_xxxx_xxxyz_0[i] = 3.0 * ta1_z_xx_xxxyz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxyz_1[i] * fe_0 + 3.0 * ta1_z_xxx_xxyz_0[i] * fe_0 -
                                3.0 * ta1_z_xxx_xxyz_1[i] * fe_0 + ta1_z_xxx_xxxyz_0[i] * pa_x[i] - ta1_z_xxx_xxxyz_1[i] * pc_x[i];

        ta1_z_xxxx_xxxzz_0[i] = 3.0 * ta1_z_xx_xxxzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxxzz_1[i] * fe_0 + 3.0 * ta1_z_xxx_xxzz_0[i] * fe_0 -
                                3.0 * ta1_z_xxx_xxzz_1[i] * fe_0 + ta1_z_xxx_xxxzz_0[i] * pa_x[i] - ta1_z_xxx_xxxzz_1[i] * pc_x[i];

        ta1_z_xxxx_xxyyy_0[i] = 3.0 * ta1_z_xx_xxyyy_0[i] * fe_0 - 3.0 * ta1_z_xx_xxyyy_1[i] * fe_0 + 2.0 * ta1_z_xxx_xyyy_0[i] * fe_0 -
                                2.0 * ta1_z_xxx_xyyy_1[i] * fe_0 + ta1_z_xxx_xxyyy_0[i] * pa_x[i] - ta1_z_xxx_xxyyy_1[i] * pc_x[i];

        ta1_z_xxxx_xxyyz_0[i] = 3.0 * ta1_z_xx_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxyyz_1[i] * fe_0 + 2.0 * ta1_z_xxx_xyyz_0[i] * fe_0 -
                                2.0 * ta1_z_xxx_xyyz_1[i] * fe_0 + ta1_z_xxx_xxyyz_0[i] * pa_x[i] - ta1_z_xxx_xxyyz_1[i] * pc_x[i];

        ta1_z_xxxx_xxyzz_0[i] = 3.0 * ta1_z_xx_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxyzz_1[i] * fe_0 + 2.0 * ta1_z_xxx_xyzz_0[i] * fe_0 -
                                2.0 * ta1_z_xxx_xyzz_1[i] * fe_0 + ta1_z_xxx_xxyzz_0[i] * pa_x[i] - ta1_z_xxx_xxyzz_1[i] * pc_x[i];

        ta1_z_xxxx_xxzzz_0[i] = 3.0 * ta1_z_xx_xxzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xxzzz_1[i] * fe_0 + 2.0 * ta1_z_xxx_xzzz_0[i] * fe_0 -
                                2.0 * ta1_z_xxx_xzzz_1[i] * fe_0 + ta1_z_xxx_xxzzz_0[i] * pa_x[i] - ta1_z_xxx_xxzzz_1[i] * pc_x[i];

        ta1_z_xxxx_xyyyy_0[i] = 3.0 * ta1_z_xx_xyyyy_0[i] * fe_0 - 3.0 * ta1_z_xx_xyyyy_1[i] * fe_0 + ta1_z_xxx_yyyy_0[i] * fe_0 -
                                ta1_z_xxx_yyyy_1[i] * fe_0 + ta1_z_xxx_xyyyy_0[i] * pa_x[i] - ta1_z_xxx_xyyyy_1[i] * pc_x[i];

        ta1_z_xxxx_xyyyz_0[i] = 3.0 * ta1_z_xx_xyyyz_0[i] * fe_0 - 3.0 * ta1_z_xx_xyyyz_1[i] * fe_0 + ta1_z_xxx_yyyz_0[i] * fe_0 -
                                ta1_z_xxx_yyyz_1[i] * fe_0 + ta1_z_xxx_xyyyz_0[i] * pa_x[i] - ta1_z_xxx_xyyyz_1[i] * pc_x[i];

        ta1_z_xxxx_xyyzz_0[i] = 3.0 * ta1_z_xx_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xyyzz_1[i] * fe_0 + ta1_z_xxx_yyzz_0[i] * fe_0 -
                                ta1_z_xxx_yyzz_1[i] * fe_0 + ta1_z_xxx_xyyzz_0[i] * pa_x[i] - ta1_z_xxx_xyyzz_1[i] * pc_x[i];

        ta1_z_xxxx_xyzzz_0[i] = 3.0 * ta1_z_xx_xyzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xyzzz_1[i] * fe_0 + ta1_z_xxx_yzzz_0[i] * fe_0 -
                                ta1_z_xxx_yzzz_1[i] * fe_0 + ta1_z_xxx_xyzzz_0[i] * pa_x[i] - ta1_z_xxx_xyzzz_1[i] * pc_x[i];

        ta1_z_xxxx_xzzzz_0[i] = 3.0 * ta1_z_xx_xzzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_xzzzz_1[i] * fe_0 + ta1_z_xxx_zzzz_0[i] * fe_0 -
                                ta1_z_xxx_zzzz_1[i] * fe_0 + ta1_z_xxx_xzzzz_0[i] * pa_x[i] - ta1_z_xxx_xzzzz_1[i] * pc_x[i];

        ta1_z_xxxx_yyyyy_0[i] =
            3.0 * ta1_z_xx_yyyyy_0[i] * fe_0 - 3.0 * ta1_z_xx_yyyyy_1[i] * fe_0 + ta1_z_xxx_yyyyy_0[i] * pa_x[i] - ta1_z_xxx_yyyyy_1[i] * pc_x[i];

        ta1_z_xxxx_yyyyz_0[i] =
            3.0 * ta1_z_xx_yyyyz_0[i] * fe_0 - 3.0 * ta1_z_xx_yyyyz_1[i] * fe_0 + ta1_z_xxx_yyyyz_0[i] * pa_x[i] - ta1_z_xxx_yyyyz_1[i] * pc_x[i];

        ta1_z_xxxx_yyyzz_0[i] =
            3.0 * ta1_z_xx_yyyzz_0[i] * fe_0 - 3.0 * ta1_z_xx_yyyzz_1[i] * fe_0 + ta1_z_xxx_yyyzz_0[i] * pa_x[i] - ta1_z_xxx_yyyzz_1[i] * pc_x[i];

        ta1_z_xxxx_yyzzz_0[i] =
            3.0 * ta1_z_xx_yyzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_yyzzz_1[i] * fe_0 + ta1_z_xxx_yyzzz_0[i] * pa_x[i] - ta1_z_xxx_yyzzz_1[i] * pc_x[i];

        ta1_z_xxxx_yzzzz_0[i] =
            3.0 * ta1_z_xx_yzzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_yzzzz_1[i] * fe_0 + ta1_z_xxx_yzzzz_0[i] * pa_x[i] - ta1_z_xxx_yzzzz_1[i] * pc_x[i];

        ta1_z_xxxx_zzzzz_0[i] =
            3.0 * ta1_z_xx_zzzzz_0[i] * fe_0 - 3.0 * ta1_z_xx_zzzzz_1[i] * fe_0 + ta1_z_xxx_zzzzz_0[i] * pa_x[i] - ta1_z_xxx_zzzzz_1[i] * pc_x[i];
    }

    // Set up 651-672 components of targeted buffer : GH

    auto ta1_z_xxxy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 651);

    auto ta1_z_xxxy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 652);

    auto ta1_z_xxxy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 653);

    auto ta1_z_xxxy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 654);

    auto ta1_z_xxxy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 655);

    auto ta1_z_xxxy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 656);

    auto ta1_z_xxxy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 657);

    auto ta1_z_xxxy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 658);

    auto ta1_z_xxxy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 659);

    auto ta1_z_xxxy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 660);

    auto ta1_z_xxxy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 661);

    auto ta1_z_xxxy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 662);

    auto ta1_z_xxxy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 663);

    auto ta1_z_xxxy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 664);

    auto ta1_z_xxxy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 665);

    auto ta1_z_xxxy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 666);

    auto ta1_z_xxxy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 667);

    auto ta1_z_xxxy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 668);

    auto ta1_z_xxxy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 669);

    auto ta1_z_xxxy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 670);

    auto ta1_z_xxxy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 671);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xxx_xxxx_0,   \
                             ta1_z_xxx_xxxx_1,   \
                             ta1_z_xxx_xxxxx_0,  \
                             ta1_z_xxx_xxxxx_1,  \
                             ta1_z_xxx_xxxxy_0,  \
                             ta1_z_xxx_xxxxy_1,  \
                             ta1_z_xxx_xxxxz_0,  \
                             ta1_z_xxx_xxxxz_1,  \
                             ta1_z_xxx_xxxy_0,   \
                             ta1_z_xxx_xxxy_1,   \
                             ta1_z_xxx_xxxyy_0,  \
                             ta1_z_xxx_xxxyy_1,  \
                             ta1_z_xxx_xxxyz_0,  \
                             ta1_z_xxx_xxxyz_1,  \
                             ta1_z_xxx_xxxz_0,   \
                             ta1_z_xxx_xxxz_1,   \
                             ta1_z_xxx_xxxzz_0,  \
                             ta1_z_xxx_xxxzz_1,  \
                             ta1_z_xxx_xxyy_0,   \
                             ta1_z_xxx_xxyy_1,   \
                             ta1_z_xxx_xxyyy_0,  \
                             ta1_z_xxx_xxyyy_1,  \
                             ta1_z_xxx_xxyyz_0,  \
                             ta1_z_xxx_xxyyz_1,  \
                             ta1_z_xxx_xxyz_0,   \
                             ta1_z_xxx_xxyz_1,   \
                             ta1_z_xxx_xxyzz_0,  \
                             ta1_z_xxx_xxyzz_1,  \
                             ta1_z_xxx_xxzz_0,   \
                             ta1_z_xxx_xxzz_1,   \
                             ta1_z_xxx_xxzzz_0,  \
                             ta1_z_xxx_xxzzz_1,  \
                             ta1_z_xxx_xyyy_0,   \
                             ta1_z_xxx_xyyy_1,   \
                             ta1_z_xxx_xyyyy_0,  \
                             ta1_z_xxx_xyyyy_1,  \
                             ta1_z_xxx_xyyyz_0,  \
                             ta1_z_xxx_xyyyz_1,  \
                             ta1_z_xxx_xyyz_0,   \
                             ta1_z_xxx_xyyz_1,   \
                             ta1_z_xxx_xyyzz_0,  \
                             ta1_z_xxx_xyyzz_1,  \
                             ta1_z_xxx_xyzz_0,   \
                             ta1_z_xxx_xyzz_1,   \
                             ta1_z_xxx_xyzzz_0,  \
                             ta1_z_xxx_xyzzz_1,  \
                             ta1_z_xxx_xzzz_0,   \
                             ta1_z_xxx_xzzz_1,   \
                             ta1_z_xxx_xzzzz_0,  \
                             ta1_z_xxx_xzzzz_1,  \
                             ta1_z_xxx_zzzzz_0,  \
                             ta1_z_xxx_zzzzz_1,  \
                             ta1_z_xxxy_xxxxx_0, \
                             ta1_z_xxxy_xxxxy_0, \
                             ta1_z_xxxy_xxxxz_0, \
                             ta1_z_xxxy_xxxyy_0, \
                             ta1_z_xxxy_xxxyz_0, \
                             ta1_z_xxxy_xxxzz_0, \
                             ta1_z_xxxy_xxyyy_0, \
                             ta1_z_xxxy_xxyyz_0, \
                             ta1_z_xxxy_xxyzz_0, \
                             ta1_z_xxxy_xxzzz_0, \
                             ta1_z_xxxy_xyyyy_0, \
                             ta1_z_xxxy_xyyyz_0, \
                             ta1_z_xxxy_xyyzz_0, \
                             ta1_z_xxxy_xyzzz_0, \
                             ta1_z_xxxy_xzzzz_0, \
                             ta1_z_xxxy_yyyyy_0, \
                             ta1_z_xxxy_yyyyz_0, \
                             ta1_z_xxxy_yyyzz_0, \
                             ta1_z_xxxy_yyzzz_0, \
                             ta1_z_xxxy_yzzzz_0, \
                             ta1_z_xxxy_zzzzz_0, \
                             ta1_z_xxy_yyyyy_0,  \
                             ta1_z_xxy_yyyyy_1,  \
                             ta1_z_xxy_yyyyz_0,  \
                             ta1_z_xxy_yyyyz_1,  \
                             ta1_z_xxy_yyyzz_0,  \
                             ta1_z_xxy_yyyzz_1,  \
                             ta1_z_xxy_yyzzz_0,  \
                             ta1_z_xxy_yyzzz_1,  \
                             ta1_z_xxy_yzzzz_0,  \
                             ta1_z_xxy_yzzzz_1,  \
                             ta1_z_xy_yyyyy_0,   \
                             ta1_z_xy_yyyyy_1,   \
                             ta1_z_xy_yyyyz_0,   \
                             ta1_z_xy_yyyyz_1,   \
                             ta1_z_xy_yyyzz_0,   \
                             ta1_z_xy_yyyzz_1,   \
                             ta1_z_xy_yyzzz_0,   \
                             ta1_z_xy_yyzzz_1,   \
                             ta1_z_xy_yzzzz_0,   \
                             ta1_z_xy_yzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxy_xxxxx_0[i] = ta1_z_xxx_xxxxx_0[i] * pa_y[i] - ta1_z_xxx_xxxxx_1[i] * pc_y[i];

        ta1_z_xxxy_xxxxy_0[i] =
            ta1_z_xxx_xxxx_0[i] * fe_0 - ta1_z_xxx_xxxx_1[i] * fe_0 + ta1_z_xxx_xxxxy_0[i] * pa_y[i] - ta1_z_xxx_xxxxy_1[i] * pc_y[i];

        ta1_z_xxxy_xxxxz_0[i] = ta1_z_xxx_xxxxz_0[i] * pa_y[i] - ta1_z_xxx_xxxxz_1[i] * pc_y[i];

        ta1_z_xxxy_xxxyy_0[i] =
            2.0 * ta1_z_xxx_xxxy_0[i] * fe_0 - 2.0 * ta1_z_xxx_xxxy_1[i] * fe_0 + ta1_z_xxx_xxxyy_0[i] * pa_y[i] - ta1_z_xxx_xxxyy_1[i] * pc_y[i];

        ta1_z_xxxy_xxxyz_0[i] =
            ta1_z_xxx_xxxz_0[i] * fe_0 - ta1_z_xxx_xxxz_1[i] * fe_0 + ta1_z_xxx_xxxyz_0[i] * pa_y[i] - ta1_z_xxx_xxxyz_1[i] * pc_y[i];

        ta1_z_xxxy_xxxzz_0[i] = ta1_z_xxx_xxxzz_0[i] * pa_y[i] - ta1_z_xxx_xxxzz_1[i] * pc_y[i];

        ta1_z_xxxy_xxyyy_0[i] =
            3.0 * ta1_z_xxx_xxyy_0[i] * fe_0 - 3.0 * ta1_z_xxx_xxyy_1[i] * fe_0 + ta1_z_xxx_xxyyy_0[i] * pa_y[i] - ta1_z_xxx_xxyyy_1[i] * pc_y[i];

        ta1_z_xxxy_xxyyz_0[i] =
            2.0 * ta1_z_xxx_xxyz_0[i] * fe_0 - 2.0 * ta1_z_xxx_xxyz_1[i] * fe_0 + ta1_z_xxx_xxyyz_0[i] * pa_y[i] - ta1_z_xxx_xxyyz_1[i] * pc_y[i];

        ta1_z_xxxy_xxyzz_0[i] =
            ta1_z_xxx_xxzz_0[i] * fe_0 - ta1_z_xxx_xxzz_1[i] * fe_0 + ta1_z_xxx_xxyzz_0[i] * pa_y[i] - ta1_z_xxx_xxyzz_1[i] * pc_y[i];

        ta1_z_xxxy_xxzzz_0[i] = ta1_z_xxx_xxzzz_0[i] * pa_y[i] - ta1_z_xxx_xxzzz_1[i] * pc_y[i];

        ta1_z_xxxy_xyyyy_0[i] =
            4.0 * ta1_z_xxx_xyyy_0[i] * fe_0 - 4.0 * ta1_z_xxx_xyyy_1[i] * fe_0 + ta1_z_xxx_xyyyy_0[i] * pa_y[i] - ta1_z_xxx_xyyyy_1[i] * pc_y[i];

        ta1_z_xxxy_xyyyz_0[i] =
            3.0 * ta1_z_xxx_xyyz_0[i] * fe_0 - 3.0 * ta1_z_xxx_xyyz_1[i] * fe_0 + ta1_z_xxx_xyyyz_0[i] * pa_y[i] - ta1_z_xxx_xyyyz_1[i] * pc_y[i];

        ta1_z_xxxy_xyyzz_0[i] =
            2.0 * ta1_z_xxx_xyzz_0[i] * fe_0 - 2.0 * ta1_z_xxx_xyzz_1[i] * fe_0 + ta1_z_xxx_xyyzz_0[i] * pa_y[i] - ta1_z_xxx_xyyzz_1[i] * pc_y[i];

        ta1_z_xxxy_xyzzz_0[i] =
            ta1_z_xxx_xzzz_0[i] * fe_0 - ta1_z_xxx_xzzz_1[i] * fe_0 + ta1_z_xxx_xyzzz_0[i] * pa_y[i] - ta1_z_xxx_xyzzz_1[i] * pc_y[i];

        ta1_z_xxxy_xzzzz_0[i] = ta1_z_xxx_xzzzz_0[i] * pa_y[i] - ta1_z_xxx_xzzzz_1[i] * pc_y[i];

        ta1_z_xxxy_yyyyy_0[i] =
            2.0 * ta1_z_xy_yyyyy_0[i] * fe_0 - 2.0 * ta1_z_xy_yyyyy_1[i] * fe_0 + ta1_z_xxy_yyyyy_0[i] * pa_x[i] - ta1_z_xxy_yyyyy_1[i] * pc_x[i];

        ta1_z_xxxy_yyyyz_0[i] =
            2.0 * ta1_z_xy_yyyyz_0[i] * fe_0 - 2.0 * ta1_z_xy_yyyyz_1[i] * fe_0 + ta1_z_xxy_yyyyz_0[i] * pa_x[i] - ta1_z_xxy_yyyyz_1[i] * pc_x[i];

        ta1_z_xxxy_yyyzz_0[i] =
            2.0 * ta1_z_xy_yyyzz_0[i] * fe_0 - 2.0 * ta1_z_xy_yyyzz_1[i] * fe_0 + ta1_z_xxy_yyyzz_0[i] * pa_x[i] - ta1_z_xxy_yyyzz_1[i] * pc_x[i];

        ta1_z_xxxy_yyzzz_0[i] =
            2.0 * ta1_z_xy_yyzzz_0[i] * fe_0 - 2.0 * ta1_z_xy_yyzzz_1[i] * fe_0 + ta1_z_xxy_yyzzz_0[i] * pa_x[i] - ta1_z_xxy_yyzzz_1[i] * pc_x[i];

        ta1_z_xxxy_yzzzz_0[i] =
            2.0 * ta1_z_xy_yzzzz_0[i] * fe_0 - 2.0 * ta1_z_xy_yzzzz_1[i] * fe_0 + ta1_z_xxy_yzzzz_0[i] * pa_x[i] - ta1_z_xxy_yzzzz_1[i] * pc_x[i];

        ta1_z_xxxy_zzzzz_0[i] = ta1_z_xxx_zzzzz_0[i] * pa_y[i] - ta1_z_xxx_zzzzz_1[i] * pc_y[i];
    }

    // Set up 672-693 components of targeted buffer : GH

    auto ta1_z_xxxz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 672);

    auto ta1_z_xxxz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 673);

    auto ta1_z_xxxz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 674);

    auto ta1_z_xxxz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 675);

    auto ta1_z_xxxz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 676);

    auto ta1_z_xxxz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 677);

    auto ta1_z_xxxz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 678);

    auto ta1_z_xxxz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 679);

    auto ta1_z_xxxz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 680);

    auto ta1_z_xxxz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 681);

    auto ta1_z_xxxz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 682);

    auto ta1_z_xxxz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 683);

    auto ta1_z_xxxz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 684);

    auto ta1_z_xxxz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 685);

    auto ta1_z_xxxz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 686);

    auto ta1_z_xxxz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 687);

    auto ta1_z_xxxz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 688);

    auto ta1_z_xxxz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 689);

    auto ta1_z_xxxz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 690);

    auto ta1_z_xxxz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 691);

    auto ta1_z_xxxz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 692);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_z_xxx_xxxx_0,   \
                             ta1_z_xxx_xxxx_1,   \
                             ta1_z_xxx_xxxxx_0,  \
                             ta1_z_xxx_xxxxx_1,  \
                             ta1_z_xxx_xxxxy_0,  \
                             ta1_z_xxx_xxxxy_1,  \
                             ta1_z_xxx_xxxxz_0,  \
                             ta1_z_xxx_xxxxz_1,  \
                             ta1_z_xxx_xxxy_0,   \
                             ta1_z_xxx_xxxy_1,   \
                             ta1_z_xxx_xxxyy_0,  \
                             ta1_z_xxx_xxxyy_1,  \
                             ta1_z_xxx_xxxyz_0,  \
                             ta1_z_xxx_xxxyz_1,  \
                             ta1_z_xxx_xxxz_0,   \
                             ta1_z_xxx_xxxz_1,   \
                             ta1_z_xxx_xxxzz_0,  \
                             ta1_z_xxx_xxxzz_1,  \
                             ta1_z_xxx_xxyy_0,   \
                             ta1_z_xxx_xxyy_1,   \
                             ta1_z_xxx_xxyyy_0,  \
                             ta1_z_xxx_xxyyy_1,  \
                             ta1_z_xxx_xxyyz_0,  \
                             ta1_z_xxx_xxyyz_1,  \
                             ta1_z_xxx_xxyz_0,   \
                             ta1_z_xxx_xxyz_1,   \
                             ta1_z_xxx_xxyzz_0,  \
                             ta1_z_xxx_xxyzz_1,  \
                             ta1_z_xxx_xxzz_0,   \
                             ta1_z_xxx_xxzz_1,   \
                             ta1_z_xxx_xxzzz_0,  \
                             ta1_z_xxx_xxzzz_1,  \
                             ta1_z_xxx_xyyy_0,   \
                             ta1_z_xxx_xyyy_1,   \
                             ta1_z_xxx_xyyyy_0,  \
                             ta1_z_xxx_xyyyy_1,  \
                             ta1_z_xxx_xyyyz_0,  \
                             ta1_z_xxx_xyyyz_1,  \
                             ta1_z_xxx_xyyz_0,   \
                             ta1_z_xxx_xyyz_1,   \
                             ta1_z_xxx_xyyzz_0,  \
                             ta1_z_xxx_xyyzz_1,  \
                             ta1_z_xxx_xyzz_0,   \
                             ta1_z_xxx_xyzz_1,   \
                             ta1_z_xxx_xyzzz_0,  \
                             ta1_z_xxx_xyzzz_1,  \
                             ta1_z_xxx_xzzz_0,   \
                             ta1_z_xxx_xzzz_1,   \
                             ta1_z_xxx_xzzzz_0,  \
                             ta1_z_xxx_xzzzz_1,  \
                             ta1_z_xxx_yyyyy_0,  \
                             ta1_z_xxx_yyyyy_1,  \
                             ta1_z_xxxz_xxxxx_0, \
                             ta1_z_xxxz_xxxxy_0, \
                             ta1_z_xxxz_xxxxz_0, \
                             ta1_z_xxxz_xxxyy_0, \
                             ta1_z_xxxz_xxxyz_0, \
                             ta1_z_xxxz_xxxzz_0, \
                             ta1_z_xxxz_xxyyy_0, \
                             ta1_z_xxxz_xxyyz_0, \
                             ta1_z_xxxz_xxyzz_0, \
                             ta1_z_xxxz_xxzzz_0, \
                             ta1_z_xxxz_xyyyy_0, \
                             ta1_z_xxxz_xyyyz_0, \
                             ta1_z_xxxz_xyyzz_0, \
                             ta1_z_xxxz_xyzzz_0, \
                             ta1_z_xxxz_xzzzz_0, \
                             ta1_z_xxxz_yyyyy_0, \
                             ta1_z_xxxz_yyyyz_0, \
                             ta1_z_xxxz_yyyzz_0, \
                             ta1_z_xxxz_yyzzz_0, \
                             ta1_z_xxxz_yzzzz_0, \
                             ta1_z_xxxz_zzzzz_0, \
                             ta1_z_xxz_yyyyz_0,  \
                             ta1_z_xxz_yyyyz_1,  \
                             ta1_z_xxz_yyyzz_0,  \
                             ta1_z_xxz_yyyzz_1,  \
                             ta1_z_xxz_yyzzz_0,  \
                             ta1_z_xxz_yyzzz_1,  \
                             ta1_z_xxz_yzzzz_0,  \
                             ta1_z_xxz_yzzzz_1,  \
                             ta1_z_xxz_zzzzz_0,  \
                             ta1_z_xxz_zzzzz_1,  \
                             ta1_z_xz_yyyyz_0,   \
                             ta1_z_xz_yyyyz_1,   \
                             ta1_z_xz_yyyzz_0,   \
                             ta1_z_xz_yyyzz_1,   \
                             ta1_z_xz_yyzzz_0,   \
                             ta1_z_xz_yyzzz_1,   \
                             ta1_z_xz_yzzzz_0,   \
                             ta1_z_xz_yzzzz_1,   \
                             ta1_z_xz_zzzzz_0,   \
                             ta1_z_xz_zzzzz_1,   \
                             ta_xxx_xxxxx_1,     \
                             ta_xxx_xxxxy_1,     \
                             ta_xxx_xxxxz_1,     \
                             ta_xxx_xxxyy_1,     \
                             ta_xxx_xxxyz_1,     \
                             ta_xxx_xxxzz_1,     \
                             ta_xxx_xxyyy_1,     \
                             ta_xxx_xxyyz_1,     \
                             ta_xxx_xxyzz_1,     \
                             ta_xxx_xxzzz_1,     \
                             ta_xxx_xyyyy_1,     \
                             ta_xxx_xyyyz_1,     \
                             ta_xxx_xyyzz_1,     \
                             ta_xxx_xyzzz_1,     \
                             ta_xxx_xzzzz_1,     \
                             ta_xxx_yyyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxxz_xxxxx_0[i] = ta_xxx_xxxxx_1[i] + ta1_z_xxx_xxxxx_0[i] * pa_z[i] - ta1_z_xxx_xxxxx_1[i] * pc_z[i];

        ta1_z_xxxz_xxxxy_0[i] = ta_xxx_xxxxy_1[i] + ta1_z_xxx_xxxxy_0[i] * pa_z[i] - ta1_z_xxx_xxxxy_1[i] * pc_z[i];

        ta1_z_xxxz_xxxxz_0[i] = ta1_z_xxx_xxxx_0[i] * fe_0 - ta1_z_xxx_xxxx_1[i] * fe_0 + ta_xxx_xxxxz_1[i] + ta1_z_xxx_xxxxz_0[i] * pa_z[i] -
                                ta1_z_xxx_xxxxz_1[i] * pc_z[i];

        ta1_z_xxxz_xxxyy_0[i] = ta_xxx_xxxyy_1[i] + ta1_z_xxx_xxxyy_0[i] * pa_z[i] - ta1_z_xxx_xxxyy_1[i] * pc_z[i];

        ta1_z_xxxz_xxxyz_0[i] = ta1_z_xxx_xxxy_0[i] * fe_0 - ta1_z_xxx_xxxy_1[i] * fe_0 + ta_xxx_xxxyz_1[i] + ta1_z_xxx_xxxyz_0[i] * pa_z[i] -
                                ta1_z_xxx_xxxyz_1[i] * pc_z[i];

        ta1_z_xxxz_xxxzz_0[i] = 2.0 * ta1_z_xxx_xxxz_0[i] * fe_0 - 2.0 * ta1_z_xxx_xxxz_1[i] * fe_0 + ta_xxx_xxxzz_1[i] +
                                ta1_z_xxx_xxxzz_0[i] * pa_z[i] - ta1_z_xxx_xxxzz_1[i] * pc_z[i];

        ta1_z_xxxz_xxyyy_0[i] = ta_xxx_xxyyy_1[i] + ta1_z_xxx_xxyyy_0[i] * pa_z[i] - ta1_z_xxx_xxyyy_1[i] * pc_z[i];

        ta1_z_xxxz_xxyyz_0[i] = ta1_z_xxx_xxyy_0[i] * fe_0 - ta1_z_xxx_xxyy_1[i] * fe_0 + ta_xxx_xxyyz_1[i] + ta1_z_xxx_xxyyz_0[i] * pa_z[i] -
                                ta1_z_xxx_xxyyz_1[i] * pc_z[i];

        ta1_z_xxxz_xxyzz_0[i] = 2.0 * ta1_z_xxx_xxyz_0[i] * fe_0 - 2.0 * ta1_z_xxx_xxyz_1[i] * fe_0 + ta_xxx_xxyzz_1[i] +
                                ta1_z_xxx_xxyzz_0[i] * pa_z[i] - ta1_z_xxx_xxyzz_1[i] * pc_z[i];

        ta1_z_xxxz_xxzzz_0[i] = 3.0 * ta1_z_xxx_xxzz_0[i] * fe_0 - 3.0 * ta1_z_xxx_xxzz_1[i] * fe_0 + ta_xxx_xxzzz_1[i] +
                                ta1_z_xxx_xxzzz_0[i] * pa_z[i] - ta1_z_xxx_xxzzz_1[i] * pc_z[i];

        ta1_z_xxxz_xyyyy_0[i] = ta_xxx_xyyyy_1[i] + ta1_z_xxx_xyyyy_0[i] * pa_z[i] - ta1_z_xxx_xyyyy_1[i] * pc_z[i];

        ta1_z_xxxz_xyyyz_0[i] = ta1_z_xxx_xyyy_0[i] * fe_0 - ta1_z_xxx_xyyy_1[i] * fe_0 + ta_xxx_xyyyz_1[i] + ta1_z_xxx_xyyyz_0[i] * pa_z[i] -
                                ta1_z_xxx_xyyyz_1[i] * pc_z[i];

        ta1_z_xxxz_xyyzz_0[i] = 2.0 * ta1_z_xxx_xyyz_0[i] * fe_0 - 2.0 * ta1_z_xxx_xyyz_1[i] * fe_0 + ta_xxx_xyyzz_1[i] +
                                ta1_z_xxx_xyyzz_0[i] * pa_z[i] - ta1_z_xxx_xyyzz_1[i] * pc_z[i];

        ta1_z_xxxz_xyzzz_0[i] = 3.0 * ta1_z_xxx_xyzz_0[i] * fe_0 - 3.0 * ta1_z_xxx_xyzz_1[i] * fe_0 + ta_xxx_xyzzz_1[i] +
                                ta1_z_xxx_xyzzz_0[i] * pa_z[i] - ta1_z_xxx_xyzzz_1[i] * pc_z[i];

        ta1_z_xxxz_xzzzz_0[i] = 4.0 * ta1_z_xxx_xzzz_0[i] * fe_0 - 4.0 * ta1_z_xxx_xzzz_1[i] * fe_0 + ta_xxx_xzzzz_1[i] +
                                ta1_z_xxx_xzzzz_0[i] * pa_z[i] - ta1_z_xxx_xzzzz_1[i] * pc_z[i];

        ta1_z_xxxz_yyyyy_0[i] = ta_xxx_yyyyy_1[i] + ta1_z_xxx_yyyyy_0[i] * pa_z[i] - ta1_z_xxx_yyyyy_1[i] * pc_z[i];

        ta1_z_xxxz_yyyyz_0[i] =
            2.0 * ta1_z_xz_yyyyz_0[i] * fe_0 - 2.0 * ta1_z_xz_yyyyz_1[i] * fe_0 + ta1_z_xxz_yyyyz_0[i] * pa_x[i] - ta1_z_xxz_yyyyz_1[i] * pc_x[i];

        ta1_z_xxxz_yyyzz_0[i] =
            2.0 * ta1_z_xz_yyyzz_0[i] * fe_0 - 2.0 * ta1_z_xz_yyyzz_1[i] * fe_0 + ta1_z_xxz_yyyzz_0[i] * pa_x[i] - ta1_z_xxz_yyyzz_1[i] * pc_x[i];

        ta1_z_xxxz_yyzzz_0[i] =
            2.0 * ta1_z_xz_yyzzz_0[i] * fe_0 - 2.0 * ta1_z_xz_yyzzz_1[i] * fe_0 + ta1_z_xxz_yyzzz_0[i] * pa_x[i] - ta1_z_xxz_yyzzz_1[i] * pc_x[i];

        ta1_z_xxxz_yzzzz_0[i] =
            2.0 * ta1_z_xz_yzzzz_0[i] * fe_0 - 2.0 * ta1_z_xz_yzzzz_1[i] * fe_0 + ta1_z_xxz_yzzzz_0[i] * pa_x[i] - ta1_z_xxz_yzzzz_1[i] * pc_x[i];

        ta1_z_xxxz_zzzzz_0[i] =
            2.0 * ta1_z_xz_zzzzz_0[i] * fe_0 - 2.0 * ta1_z_xz_zzzzz_1[i] * fe_0 + ta1_z_xxz_zzzzz_0[i] * pa_x[i] - ta1_z_xxz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 693-714 components of targeted buffer : GH

    auto ta1_z_xxyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 693);

    auto ta1_z_xxyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 694);

    auto ta1_z_xxyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 695);

    auto ta1_z_xxyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 696);

    auto ta1_z_xxyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 697);

    auto ta1_z_xxyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 698);

    auto ta1_z_xxyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 699);

    auto ta1_z_xxyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 700);

    auto ta1_z_xxyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 701);

    auto ta1_z_xxyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 702);

    auto ta1_z_xxyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 703);

    auto ta1_z_xxyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 704);

    auto ta1_z_xxyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 705);

    auto ta1_z_xxyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 706);

    auto ta1_z_xxyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 707);

    auto ta1_z_xxyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 708);

    auto ta1_z_xxyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 709);

    auto ta1_z_xxyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 710);

    auto ta1_z_xxyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 711);

    auto ta1_z_xxyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 712);

    auto ta1_z_xxyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 713);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xx_xxxxx_0,   \
                             ta1_z_xx_xxxxx_1,   \
                             ta1_z_xx_xxxxz_0,   \
                             ta1_z_xx_xxxxz_1,   \
                             ta1_z_xx_xxxzz_0,   \
                             ta1_z_xx_xxxzz_1,   \
                             ta1_z_xx_xxzzz_0,   \
                             ta1_z_xx_xxzzz_1,   \
                             ta1_z_xx_xzzzz_0,   \
                             ta1_z_xx_xzzzz_1,   \
                             ta1_z_xxy_xxxxx_0,  \
                             ta1_z_xxy_xxxxx_1,  \
                             ta1_z_xxy_xxxxz_0,  \
                             ta1_z_xxy_xxxxz_1,  \
                             ta1_z_xxy_xxxzz_0,  \
                             ta1_z_xxy_xxxzz_1,  \
                             ta1_z_xxy_xxzzz_0,  \
                             ta1_z_xxy_xxzzz_1,  \
                             ta1_z_xxy_xzzzz_0,  \
                             ta1_z_xxy_xzzzz_1,  \
                             ta1_z_xxyy_xxxxx_0, \
                             ta1_z_xxyy_xxxxy_0, \
                             ta1_z_xxyy_xxxxz_0, \
                             ta1_z_xxyy_xxxyy_0, \
                             ta1_z_xxyy_xxxyz_0, \
                             ta1_z_xxyy_xxxzz_0, \
                             ta1_z_xxyy_xxyyy_0, \
                             ta1_z_xxyy_xxyyz_0, \
                             ta1_z_xxyy_xxyzz_0, \
                             ta1_z_xxyy_xxzzz_0, \
                             ta1_z_xxyy_xyyyy_0, \
                             ta1_z_xxyy_xyyyz_0, \
                             ta1_z_xxyy_xyyzz_0, \
                             ta1_z_xxyy_xyzzz_0, \
                             ta1_z_xxyy_xzzzz_0, \
                             ta1_z_xxyy_yyyyy_0, \
                             ta1_z_xxyy_yyyyz_0, \
                             ta1_z_xxyy_yyyzz_0, \
                             ta1_z_xxyy_yyzzz_0, \
                             ta1_z_xxyy_yzzzz_0, \
                             ta1_z_xxyy_zzzzz_0, \
                             ta1_z_xyy_xxxxy_0,  \
                             ta1_z_xyy_xxxxy_1,  \
                             ta1_z_xyy_xxxy_0,   \
                             ta1_z_xyy_xxxy_1,   \
                             ta1_z_xyy_xxxyy_0,  \
                             ta1_z_xyy_xxxyy_1,  \
                             ta1_z_xyy_xxxyz_0,  \
                             ta1_z_xyy_xxxyz_1,  \
                             ta1_z_xyy_xxyy_0,   \
                             ta1_z_xyy_xxyy_1,   \
                             ta1_z_xyy_xxyyy_0,  \
                             ta1_z_xyy_xxyyy_1,  \
                             ta1_z_xyy_xxyyz_0,  \
                             ta1_z_xyy_xxyyz_1,  \
                             ta1_z_xyy_xxyz_0,   \
                             ta1_z_xyy_xxyz_1,   \
                             ta1_z_xyy_xxyzz_0,  \
                             ta1_z_xyy_xxyzz_1,  \
                             ta1_z_xyy_xyyy_0,   \
                             ta1_z_xyy_xyyy_1,   \
                             ta1_z_xyy_xyyyy_0,  \
                             ta1_z_xyy_xyyyy_1,  \
                             ta1_z_xyy_xyyyz_0,  \
                             ta1_z_xyy_xyyyz_1,  \
                             ta1_z_xyy_xyyz_0,   \
                             ta1_z_xyy_xyyz_1,   \
                             ta1_z_xyy_xyyzz_0,  \
                             ta1_z_xyy_xyyzz_1,  \
                             ta1_z_xyy_xyzz_0,   \
                             ta1_z_xyy_xyzz_1,   \
                             ta1_z_xyy_xyzzz_0,  \
                             ta1_z_xyy_xyzzz_1,  \
                             ta1_z_xyy_yyyy_0,   \
                             ta1_z_xyy_yyyy_1,   \
                             ta1_z_xyy_yyyyy_0,  \
                             ta1_z_xyy_yyyyy_1,  \
                             ta1_z_xyy_yyyyz_0,  \
                             ta1_z_xyy_yyyyz_1,  \
                             ta1_z_xyy_yyyz_0,   \
                             ta1_z_xyy_yyyz_1,   \
                             ta1_z_xyy_yyyzz_0,  \
                             ta1_z_xyy_yyyzz_1,  \
                             ta1_z_xyy_yyzz_0,   \
                             ta1_z_xyy_yyzz_1,   \
                             ta1_z_xyy_yyzzz_0,  \
                             ta1_z_xyy_yyzzz_1,  \
                             ta1_z_xyy_yzzz_0,   \
                             ta1_z_xyy_yzzz_1,   \
                             ta1_z_xyy_yzzzz_0,  \
                             ta1_z_xyy_yzzzz_1,  \
                             ta1_z_xyy_zzzzz_0,  \
                             ta1_z_xyy_zzzzz_1,  \
                             ta1_z_yy_xxxxy_0,   \
                             ta1_z_yy_xxxxy_1,   \
                             ta1_z_yy_xxxyy_0,   \
                             ta1_z_yy_xxxyy_1,   \
                             ta1_z_yy_xxxyz_0,   \
                             ta1_z_yy_xxxyz_1,   \
                             ta1_z_yy_xxyyy_0,   \
                             ta1_z_yy_xxyyy_1,   \
                             ta1_z_yy_xxyyz_0,   \
                             ta1_z_yy_xxyyz_1,   \
                             ta1_z_yy_xxyzz_0,   \
                             ta1_z_yy_xxyzz_1,   \
                             ta1_z_yy_xyyyy_0,   \
                             ta1_z_yy_xyyyy_1,   \
                             ta1_z_yy_xyyyz_0,   \
                             ta1_z_yy_xyyyz_1,   \
                             ta1_z_yy_xyyzz_0,   \
                             ta1_z_yy_xyyzz_1,   \
                             ta1_z_yy_xyzzz_0,   \
                             ta1_z_yy_xyzzz_1,   \
                             ta1_z_yy_yyyyy_0,   \
                             ta1_z_yy_yyyyy_1,   \
                             ta1_z_yy_yyyyz_0,   \
                             ta1_z_yy_yyyyz_1,   \
                             ta1_z_yy_yyyzz_0,   \
                             ta1_z_yy_yyyzz_1,   \
                             ta1_z_yy_yyzzz_0,   \
                             ta1_z_yy_yyzzz_1,   \
                             ta1_z_yy_yzzzz_0,   \
                             ta1_z_yy_yzzzz_1,   \
                             ta1_z_yy_zzzzz_0,   \
                             ta1_z_yy_zzzzz_1,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyy_xxxxx_0[i] =
            ta1_z_xx_xxxxx_0[i] * fe_0 - ta1_z_xx_xxxxx_1[i] * fe_0 + ta1_z_xxy_xxxxx_0[i] * pa_y[i] - ta1_z_xxy_xxxxx_1[i] * pc_y[i];

        ta1_z_xxyy_xxxxy_0[i] = ta1_z_yy_xxxxy_0[i] * fe_0 - ta1_z_yy_xxxxy_1[i] * fe_0 + 4.0 * ta1_z_xyy_xxxy_0[i] * fe_0 -
                                4.0 * ta1_z_xyy_xxxy_1[i] * fe_0 + ta1_z_xyy_xxxxy_0[i] * pa_x[i] - ta1_z_xyy_xxxxy_1[i] * pc_x[i];

        ta1_z_xxyy_xxxxz_0[i] =
            ta1_z_xx_xxxxz_0[i] * fe_0 - ta1_z_xx_xxxxz_1[i] * fe_0 + ta1_z_xxy_xxxxz_0[i] * pa_y[i] - ta1_z_xxy_xxxxz_1[i] * pc_y[i];

        ta1_z_xxyy_xxxyy_0[i] = ta1_z_yy_xxxyy_0[i] * fe_0 - ta1_z_yy_xxxyy_1[i] * fe_0 + 3.0 * ta1_z_xyy_xxyy_0[i] * fe_0 -
                                3.0 * ta1_z_xyy_xxyy_1[i] * fe_0 + ta1_z_xyy_xxxyy_0[i] * pa_x[i] - ta1_z_xyy_xxxyy_1[i] * pc_x[i];

        ta1_z_xxyy_xxxyz_0[i] = ta1_z_yy_xxxyz_0[i] * fe_0 - ta1_z_yy_xxxyz_1[i] * fe_0 + 3.0 * ta1_z_xyy_xxyz_0[i] * fe_0 -
                                3.0 * ta1_z_xyy_xxyz_1[i] * fe_0 + ta1_z_xyy_xxxyz_0[i] * pa_x[i] - ta1_z_xyy_xxxyz_1[i] * pc_x[i];

        ta1_z_xxyy_xxxzz_0[i] =
            ta1_z_xx_xxxzz_0[i] * fe_0 - ta1_z_xx_xxxzz_1[i] * fe_0 + ta1_z_xxy_xxxzz_0[i] * pa_y[i] - ta1_z_xxy_xxxzz_1[i] * pc_y[i];

        ta1_z_xxyy_xxyyy_0[i] = ta1_z_yy_xxyyy_0[i] * fe_0 - ta1_z_yy_xxyyy_1[i] * fe_0 + 2.0 * ta1_z_xyy_xyyy_0[i] * fe_0 -
                                2.0 * ta1_z_xyy_xyyy_1[i] * fe_0 + ta1_z_xyy_xxyyy_0[i] * pa_x[i] - ta1_z_xyy_xxyyy_1[i] * pc_x[i];

        ta1_z_xxyy_xxyyz_0[i] = ta1_z_yy_xxyyz_0[i] * fe_0 - ta1_z_yy_xxyyz_1[i] * fe_0 + 2.0 * ta1_z_xyy_xyyz_0[i] * fe_0 -
                                2.0 * ta1_z_xyy_xyyz_1[i] * fe_0 + ta1_z_xyy_xxyyz_0[i] * pa_x[i] - ta1_z_xyy_xxyyz_1[i] * pc_x[i];

        ta1_z_xxyy_xxyzz_0[i] = ta1_z_yy_xxyzz_0[i] * fe_0 - ta1_z_yy_xxyzz_1[i] * fe_0 + 2.0 * ta1_z_xyy_xyzz_0[i] * fe_0 -
                                2.0 * ta1_z_xyy_xyzz_1[i] * fe_0 + ta1_z_xyy_xxyzz_0[i] * pa_x[i] - ta1_z_xyy_xxyzz_1[i] * pc_x[i];

        ta1_z_xxyy_xxzzz_0[i] =
            ta1_z_xx_xxzzz_0[i] * fe_0 - ta1_z_xx_xxzzz_1[i] * fe_0 + ta1_z_xxy_xxzzz_0[i] * pa_y[i] - ta1_z_xxy_xxzzz_1[i] * pc_y[i];

        ta1_z_xxyy_xyyyy_0[i] = ta1_z_yy_xyyyy_0[i] * fe_0 - ta1_z_yy_xyyyy_1[i] * fe_0 + ta1_z_xyy_yyyy_0[i] * fe_0 - ta1_z_xyy_yyyy_1[i] * fe_0 +
                                ta1_z_xyy_xyyyy_0[i] * pa_x[i] - ta1_z_xyy_xyyyy_1[i] * pc_x[i];

        ta1_z_xxyy_xyyyz_0[i] = ta1_z_yy_xyyyz_0[i] * fe_0 - ta1_z_yy_xyyyz_1[i] * fe_0 + ta1_z_xyy_yyyz_0[i] * fe_0 - ta1_z_xyy_yyyz_1[i] * fe_0 +
                                ta1_z_xyy_xyyyz_0[i] * pa_x[i] - ta1_z_xyy_xyyyz_1[i] * pc_x[i];

        ta1_z_xxyy_xyyzz_0[i] = ta1_z_yy_xyyzz_0[i] * fe_0 - ta1_z_yy_xyyzz_1[i] * fe_0 + ta1_z_xyy_yyzz_0[i] * fe_0 - ta1_z_xyy_yyzz_1[i] * fe_0 +
                                ta1_z_xyy_xyyzz_0[i] * pa_x[i] - ta1_z_xyy_xyyzz_1[i] * pc_x[i];

        ta1_z_xxyy_xyzzz_0[i] = ta1_z_yy_xyzzz_0[i] * fe_0 - ta1_z_yy_xyzzz_1[i] * fe_0 + ta1_z_xyy_yzzz_0[i] * fe_0 - ta1_z_xyy_yzzz_1[i] * fe_0 +
                                ta1_z_xyy_xyzzz_0[i] * pa_x[i] - ta1_z_xyy_xyzzz_1[i] * pc_x[i];

        ta1_z_xxyy_xzzzz_0[i] =
            ta1_z_xx_xzzzz_0[i] * fe_0 - ta1_z_xx_xzzzz_1[i] * fe_0 + ta1_z_xxy_xzzzz_0[i] * pa_y[i] - ta1_z_xxy_xzzzz_1[i] * pc_y[i];

        ta1_z_xxyy_yyyyy_0[i] =
            ta1_z_yy_yyyyy_0[i] * fe_0 - ta1_z_yy_yyyyy_1[i] * fe_0 + ta1_z_xyy_yyyyy_0[i] * pa_x[i] - ta1_z_xyy_yyyyy_1[i] * pc_x[i];

        ta1_z_xxyy_yyyyz_0[i] =
            ta1_z_yy_yyyyz_0[i] * fe_0 - ta1_z_yy_yyyyz_1[i] * fe_0 + ta1_z_xyy_yyyyz_0[i] * pa_x[i] - ta1_z_xyy_yyyyz_1[i] * pc_x[i];

        ta1_z_xxyy_yyyzz_0[i] =
            ta1_z_yy_yyyzz_0[i] * fe_0 - ta1_z_yy_yyyzz_1[i] * fe_0 + ta1_z_xyy_yyyzz_0[i] * pa_x[i] - ta1_z_xyy_yyyzz_1[i] * pc_x[i];

        ta1_z_xxyy_yyzzz_0[i] =
            ta1_z_yy_yyzzz_0[i] * fe_0 - ta1_z_yy_yyzzz_1[i] * fe_0 + ta1_z_xyy_yyzzz_0[i] * pa_x[i] - ta1_z_xyy_yyzzz_1[i] * pc_x[i];

        ta1_z_xxyy_yzzzz_0[i] =
            ta1_z_yy_yzzzz_0[i] * fe_0 - ta1_z_yy_yzzzz_1[i] * fe_0 + ta1_z_xyy_yzzzz_0[i] * pa_x[i] - ta1_z_xyy_yzzzz_1[i] * pc_x[i];

        ta1_z_xxyy_zzzzz_0[i] =
            ta1_z_yy_zzzzz_0[i] * fe_0 - ta1_z_yy_zzzzz_1[i] * fe_0 + ta1_z_xyy_zzzzz_0[i] * pa_x[i] - ta1_z_xyy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 714-735 components of targeted buffer : GH

    auto ta1_z_xxyz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 714);

    auto ta1_z_xxyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 715);

    auto ta1_z_xxyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 716);

    auto ta1_z_xxyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 717);

    auto ta1_z_xxyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 718);

    auto ta1_z_xxyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 719);

    auto ta1_z_xxyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 720);

    auto ta1_z_xxyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 721);

    auto ta1_z_xxyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 722);

    auto ta1_z_xxyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 723);

    auto ta1_z_xxyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 724);

    auto ta1_z_xxyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 725);

    auto ta1_z_xxyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 726);

    auto ta1_z_xxyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 727);

    auto ta1_z_xxyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 728);

    auto ta1_z_xxyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 729);

    auto ta1_z_xxyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 730);

    auto ta1_z_xxyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 731);

    auto ta1_z_xxyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 732);

    auto ta1_z_xxyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 733);

    auto ta1_z_xxyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 734);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pa_z,               \
                             pc_x,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_xxy_xxxxy_0,  \
                             ta1_z_xxy_xxxxy_1,  \
                             ta1_z_xxy_xxxyy_0,  \
                             ta1_z_xxy_xxxyy_1,  \
                             ta1_z_xxy_xxyyy_0,  \
                             ta1_z_xxy_xxyyy_1,  \
                             ta1_z_xxy_xyyyy_0,  \
                             ta1_z_xxy_xyyyy_1,  \
                             ta1_z_xxy_yyyyy_0,  \
                             ta1_z_xxy_yyyyy_1,  \
                             ta1_z_xxyz_xxxxx_0, \
                             ta1_z_xxyz_xxxxy_0, \
                             ta1_z_xxyz_xxxxz_0, \
                             ta1_z_xxyz_xxxyy_0, \
                             ta1_z_xxyz_xxxyz_0, \
                             ta1_z_xxyz_xxxzz_0, \
                             ta1_z_xxyz_xxyyy_0, \
                             ta1_z_xxyz_xxyyz_0, \
                             ta1_z_xxyz_xxyzz_0, \
                             ta1_z_xxyz_xxzzz_0, \
                             ta1_z_xxyz_xyyyy_0, \
                             ta1_z_xxyz_xyyyz_0, \
                             ta1_z_xxyz_xyyzz_0, \
                             ta1_z_xxyz_xyzzz_0, \
                             ta1_z_xxyz_xzzzz_0, \
                             ta1_z_xxyz_yyyyy_0, \
                             ta1_z_xxyz_yyyyz_0, \
                             ta1_z_xxyz_yyyzz_0, \
                             ta1_z_xxyz_yyzzz_0, \
                             ta1_z_xxyz_yzzzz_0, \
                             ta1_z_xxyz_zzzzz_0, \
                             ta1_z_xxz_xxxxx_0,  \
                             ta1_z_xxz_xxxxx_1,  \
                             ta1_z_xxz_xxxxz_0,  \
                             ta1_z_xxz_xxxxz_1,  \
                             ta1_z_xxz_xxxyz_0,  \
                             ta1_z_xxz_xxxyz_1,  \
                             ta1_z_xxz_xxxz_0,   \
                             ta1_z_xxz_xxxz_1,   \
                             ta1_z_xxz_xxxzz_0,  \
                             ta1_z_xxz_xxxzz_1,  \
                             ta1_z_xxz_xxyyz_0,  \
                             ta1_z_xxz_xxyyz_1,  \
                             ta1_z_xxz_xxyz_0,   \
                             ta1_z_xxz_xxyz_1,   \
                             ta1_z_xxz_xxyzz_0,  \
                             ta1_z_xxz_xxyzz_1,  \
                             ta1_z_xxz_xxzz_0,   \
                             ta1_z_xxz_xxzz_1,   \
                             ta1_z_xxz_xxzzz_0,  \
                             ta1_z_xxz_xxzzz_1,  \
                             ta1_z_xxz_xyyyz_0,  \
                             ta1_z_xxz_xyyyz_1,  \
                             ta1_z_xxz_xyyz_0,   \
                             ta1_z_xxz_xyyz_1,   \
                             ta1_z_xxz_xyyzz_0,  \
                             ta1_z_xxz_xyyzz_1,  \
                             ta1_z_xxz_xyzz_0,   \
                             ta1_z_xxz_xyzz_1,   \
                             ta1_z_xxz_xyzzz_0,  \
                             ta1_z_xxz_xyzzz_1,  \
                             ta1_z_xxz_xzzz_0,   \
                             ta1_z_xxz_xzzz_1,   \
                             ta1_z_xxz_xzzzz_0,  \
                             ta1_z_xxz_xzzzz_1,  \
                             ta1_z_xxz_zzzzz_0,  \
                             ta1_z_xxz_zzzzz_1,  \
                             ta1_z_xyz_yyyyz_0,  \
                             ta1_z_xyz_yyyyz_1,  \
                             ta1_z_xyz_yyyzz_0,  \
                             ta1_z_xyz_yyyzz_1,  \
                             ta1_z_xyz_yyzzz_0,  \
                             ta1_z_xyz_yyzzz_1,  \
                             ta1_z_xyz_yzzzz_0,  \
                             ta1_z_xyz_yzzzz_1,  \
                             ta1_z_yz_yyyyz_0,   \
                             ta1_z_yz_yyyyz_1,   \
                             ta1_z_yz_yyyzz_0,   \
                             ta1_z_yz_yyyzz_1,   \
                             ta1_z_yz_yyzzz_0,   \
                             ta1_z_yz_yyzzz_1,   \
                             ta1_z_yz_yzzzz_0,   \
                             ta1_z_yz_yzzzz_1,   \
                             ta_xxy_xxxxy_1,     \
                             ta_xxy_xxxyy_1,     \
                             ta_xxy_xxyyy_1,     \
                             ta_xxy_xyyyy_1,     \
                             ta_xxy_yyyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxyz_xxxxx_0[i] = ta1_z_xxz_xxxxx_0[i] * pa_y[i] - ta1_z_xxz_xxxxx_1[i] * pc_y[i];

        ta1_z_xxyz_xxxxy_0[i] = ta_xxy_xxxxy_1[i] + ta1_z_xxy_xxxxy_0[i] * pa_z[i] - ta1_z_xxy_xxxxy_1[i] * pc_z[i];

        ta1_z_xxyz_xxxxz_0[i] = ta1_z_xxz_xxxxz_0[i] * pa_y[i] - ta1_z_xxz_xxxxz_1[i] * pc_y[i];

        ta1_z_xxyz_xxxyy_0[i] = ta_xxy_xxxyy_1[i] + ta1_z_xxy_xxxyy_0[i] * pa_z[i] - ta1_z_xxy_xxxyy_1[i] * pc_z[i];

        ta1_z_xxyz_xxxyz_0[i] =
            ta1_z_xxz_xxxz_0[i] * fe_0 - ta1_z_xxz_xxxz_1[i] * fe_0 + ta1_z_xxz_xxxyz_0[i] * pa_y[i] - ta1_z_xxz_xxxyz_1[i] * pc_y[i];

        ta1_z_xxyz_xxxzz_0[i] = ta1_z_xxz_xxxzz_0[i] * pa_y[i] - ta1_z_xxz_xxxzz_1[i] * pc_y[i];

        ta1_z_xxyz_xxyyy_0[i] = ta_xxy_xxyyy_1[i] + ta1_z_xxy_xxyyy_0[i] * pa_z[i] - ta1_z_xxy_xxyyy_1[i] * pc_z[i];

        ta1_z_xxyz_xxyyz_0[i] =
            2.0 * ta1_z_xxz_xxyz_0[i] * fe_0 - 2.0 * ta1_z_xxz_xxyz_1[i] * fe_0 + ta1_z_xxz_xxyyz_0[i] * pa_y[i] - ta1_z_xxz_xxyyz_1[i] * pc_y[i];

        ta1_z_xxyz_xxyzz_0[i] =
            ta1_z_xxz_xxzz_0[i] * fe_0 - ta1_z_xxz_xxzz_1[i] * fe_0 + ta1_z_xxz_xxyzz_0[i] * pa_y[i] - ta1_z_xxz_xxyzz_1[i] * pc_y[i];

        ta1_z_xxyz_xxzzz_0[i] = ta1_z_xxz_xxzzz_0[i] * pa_y[i] - ta1_z_xxz_xxzzz_1[i] * pc_y[i];

        ta1_z_xxyz_xyyyy_0[i] = ta_xxy_xyyyy_1[i] + ta1_z_xxy_xyyyy_0[i] * pa_z[i] - ta1_z_xxy_xyyyy_1[i] * pc_z[i];

        ta1_z_xxyz_xyyyz_0[i] =
            3.0 * ta1_z_xxz_xyyz_0[i] * fe_0 - 3.0 * ta1_z_xxz_xyyz_1[i] * fe_0 + ta1_z_xxz_xyyyz_0[i] * pa_y[i] - ta1_z_xxz_xyyyz_1[i] * pc_y[i];

        ta1_z_xxyz_xyyzz_0[i] =
            2.0 * ta1_z_xxz_xyzz_0[i] * fe_0 - 2.0 * ta1_z_xxz_xyzz_1[i] * fe_0 + ta1_z_xxz_xyyzz_0[i] * pa_y[i] - ta1_z_xxz_xyyzz_1[i] * pc_y[i];

        ta1_z_xxyz_xyzzz_0[i] =
            ta1_z_xxz_xzzz_0[i] * fe_0 - ta1_z_xxz_xzzz_1[i] * fe_0 + ta1_z_xxz_xyzzz_0[i] * pa_y[i] - ta1_z_xxz_xyzzz_1[i] * pc_y[i];

        ta1_z_xxyz_xzzzz_0[i] = ta1_z_xxz_xzzzz_0[i] * pa_y[i] - ta1_z_xxz_xzzzz_1[i] * pc_y[i];

        ta1_z_xxyz_yyyyy_0[i] = ta_xxy_yyyyy_1[i] + ta1_z_xxy_yyyyy_0[i] * pa_z[i] - ta1_z_xxy_yyyyy_1[i] * pc_z[i];

        ta1_z_xxyz_yyyyz_0[i] =
            ta1_z_yz_yyyyz_0[i] * fe_0 - ta1_z_yz_yyyyz_1[i] * fe_0 + ta1_z_xyz_yyyyz_0[i] * pa_x[i] - ta1_z_xyz_yyyyz_1[i] * pc_x[i];

        ta1_z_xxyz_yyyzz_0[i] =
            ta1_z_yz_yyyzz_0[i] * fe_0 - ta1_z_yz_yyyzz_1[i] * fe_0 + ta1_z_xyz_yyyzz_0[i] * pa_x[i] - ta1_z_xyz_yyyzz_1[i] * pc_x[i];

        ta1_z_xxyz_yyzzz_0[i] =
            ta1_z_yz_yyzzz_0[i] * fe_0 - ta1_z_yz_yyzzz_1[i] * fe_0 + ta1_z_xyz_yyzzz_0[i] * pa_x[i] - ta1_z_xyz_yyzzz_1[i] * pc_x[i];

        ta1_z_xxyz_yzzzz_0[i] =
            ta1_z_yz_yzzzz_0[i] * fe_0 - ta1_z_yz_yzzzz_1[i] * fe_0 + ta1_z_xyz_yzzzz_0[i] * pa_x[i] - ta1_z_xyz_yzzzz_1[i] * pc_x[i];

        ta1_z_xxyz_zzzzz_0[i] = ta1_z_xxz_zzzzz_0[i] * pa_y[i] - ta1_z_xxz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 735-756 components of targeted buffer : GH

    auto ta1_z_xxzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 735);

    auto ta1_z_xxzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 736);

    auto ta1_z_xxzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 737);

    auto ta1_z_xxzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 738);

    auto ta1_z_xxzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 739);

    auto ta1_z_xxzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 740);

    auto ta1_z_xxzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 741);

    auto ta1_z_xxzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 742);

    auto ta1_z_xxzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 743);

    auto ta1_z_xxzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 744);

    auto ta1_z_xxzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 745);

    auto ta1_z_xxzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 746);

    auto ta1_z_xxzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 747);

    auto ta1_z_xxzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 748);

    auto ta1_z_xxzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 749);

    auto ta1_z_xxzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 750);

    auto ta1_z_xxzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 751);

    auto ta1_z_xxzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 752);

    auto ta1_z_xxzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 753);

    auto ta1_z_xxzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 754);

    auto ta1_z_xxzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 755);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_z_xx_xxxxx_0,   \
                             ta1_z_xx_xxxxx_1,   \
                             ta1_z_xx_xxxxy_0,   \
                             ta1_z_xx_xxxxy_1,   \
                             ta1_z_xx_xxxyy_0,   \
                             ta1_z_xx_xxxyy_1,   \
                             ta1_z_xx_xxyyy_0,   \
                             ta1_z_xx_xxyyy_1,   \
                             ta1_z_xx_xyyyy_0,   \
                             ta1_z_xx_xyyyy_1,   \
                             ta1_z_xxz_xxxxx_0,  \
                             ta1_z_xxz_xxxxx_1,  \
                             ta1_z_xxz_xxxxy_0,  \
                             ta1_z_xxz_xxxxy_1,  \
                             ta1_z_xxz_xxxyy_0,  \
                             ta1_z_xxz_xxxyy_1,  \
                             ta1_z_xxz_xxyyy_0,  \
                             ta1_z_xxz_xxyyy_1,  \
                             ta1_z_xxz_xyyyy_0,  \
                             ta1_z_xxz_xyyyy_1,  \
                             ta1_z_xxzz_xxxxx_0, \
                             ta1_z_xxzz_xxxxy_0, \
                             ta1_z_xxzz_xxxxz_0, \
                             ta1_z_xxzz_xxxyy_0, \
                             ta1_z_xxzz_xxxyz_0, \
                             ta1_z_xxzz_xxxzz_0, \
                             ta1_z_xxzz_xxyyy_0, \
                             ta1_z_xxzz_xxyyz_0, \
                             ta1_z_xxzz_xxyzz_0, \
                             ta1_z_xxzz_xxzzz_0, \
                             ta1_z_xxzz_xyyyy_0, \
                             ta1_z_xxzz_xyyyz_0, \
                             ta1_z_xxzz_xyyzz_0, \
                             ta1_z_xxzz_xyzzz_0, \
                             ta1_z_xxzz_xzzzz_0, \
                             ta1_z_xxzz_yyyyy_0, \
                             ta1_z_xxzz_yyyyz_0, \
                             ta1_z_xxzz_yyyzz_0, \
                             ta1_z_xxzz_yyzzz_0, \
                             ta1_z_xxzz_yzzzz_0, \
                             ta1_z_xxzz_zzzzz_0, \
                             ta1_z_xzz_xxxxz_0,  \
                             ta1_z_xzz_xxxxz_1,  \
                             ta1_z_xzz_xxxyz_0,  \
                             ta1_z_xzz_xxxyz_1,  \
                             ta1_z_xzz_xxxz_0,   \
                             ta1_z_xzz_xxxz_1,   \
                             ta1_z_xzz_xxxzz_0,  \
                             ta1_z_xzz_xxxzz_1,  \
                             ta1_z_xzz_xxyyz_0,  \
                             ta1_z_xzz_xxyyz_1,  \
                             ta1_z_xzz_xxyz_0,   \
                             ta1_z_xzz_xxyz_1,   \
                             ta1_z_xzz_xxyzz_0,  \
                             ta1_z_xzz_xxyzz_1,  \
                             ta1_z_xzz_xxzz_0,   \
                             ta1_z_xzz_xxzz_1,   \
                             ta1_z_xzz_xxzzz_0,  \
                             ta1_z_xzz_xxzzz_1,  \
                             ta1_z_xzz_xyyyz_0,  \
                             ta1_z_xzz_xyyyz_1,  \
                             ta1_z_xzz_xyyz_0,   \
                             ta1_z_xzz_xyyz_1,   \
                             ta1_z_xzz_xyyzz_0,  \
                             ta1_z_xzz_xyyzz_1,  \
                             ta1_z_xzz_xyzz_0,   \
                             ta1_z_xzz_xyzz_1,   \
                             ta1_z_xzz_xyzzz_0,  \
                             ta1_z_xzz_xyzzz_1,  \
                             ta1_z_xzz_xzzz_0,   \
                             ta1_z_xzz_xzzz_1,   \
                             ta1_z_xzz_xzzzz_0,  \
                             ta1_z_xzz_xzzzz_1,  \
                             ta1_z_xzz_yyyyy_0,  \
                             ta1_z_xzz_yyyyy_1,  \
                             ta1_z_xzz_yyyyz_0,  \
                             ta1_z_xzz_yyyyz_1,  \
                             ta1_z_xzz_yyyz_0,   \
                             ta1_z_xzz_yyyz_1,   \
                             ta1_z_xzz_yyyzz_0,  \
                             ta1_z_xzz_yyyzz_1,  \
                             ta1_z_xzz_yyzz_0,   \
                             ta1_z_xzz_yyzz_1,   \
                             ta1_z_xzz_yyzzz_0,  \
                             ta1_z_xzz_yyzzz_1,  \
                             ta1_z_xzz_yzzz_0,   \
                             ta1_z_xzz_yzzz_1,   \
                             ta1_z_xzz_yzzzz_0,  \
                             ta1_z_xzz_yzzzz_1,  \
                             ta1_z_xzz_zzzz_0,   \
                             ta1_z_xzz_zzzz_1,   \
                             ta1_z_xzz_zzzzz_0,  \
                             ta1_z_xzz_zzzzz_1,  \
                             ta1_z_zz_xxxxz_0,   \
                             ta1_z_zz_xxxxz_1,   \
                             ta1_z_zz_xxxyz_0,   \
                             ta1_z_zz_xxxyz_1,   \
                             ta1_z_zz_xxxzz_0,   \
                             ta1_z_zz_xxxzz_1,   \
                             ta1_z_zz_xxyyz_0,   \
                             ta1_z_zz_xxyyz_1,   \
                             ta1_z_zz_xxyzz_0,   \
                             ta1_z_zz_xxyzz_1,   \
                             ta1_z_zz_xxzzz_0,   \
                             ta1_z_zz_xxzzz_1,   \
                             ta1_z_zz_xyyyz_0,   \
                             ta1_z_zz_xyyyz_1,   \
                             ta1_z_zz_xyyzz_0,   \
                             ta1_z_zz_xyyzz_1,   \
                             ta1_z_zz_xyzzz_0,   \
                             ta1_z_zz_xyzzz_1,   \
                             ta1_z_zz_xzzzz_0,   \
                             ta1_z_zz_xzzzz_1,   \
                             ta1_z_zz_yyyyy_0,   \
                             ta1_z_zz_yyyyy_1,   \
                             ta1_z_zz_yyyyz_0,   \
                             ta1_z_zz_yyyyz_1,   \
                             ta1_z_zz_yyyzz_0,   \
                             ta1_z_zz_yyyzz_1,   \
                             ta1_z_zz_yyzzz_0,   \
                             ta1_z_zz_yyzzz_1,   \
                             ta1_z_zz_yzzzz_0,   \
                             ta1_z_zz_yzzzz_1,   \
                             ta1_z_zz_zzzzz_0,   \
                             ta1_z_zz_zzzzz_1,   \
                             ta_xxz_xxxxx_1,     \
                             ta_xxz_xxxxy_1,     \
                             ta_xxz_xxxyy_1,     \
                             ta_xxz_xxyyy_1,     \
                             ta_xxz_xyyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xxzz_xxxxx_0[i] = ta1_z_xx_xxxxx_0[i] * fe_0 - ta1_z_xx_xxxxx_1[i] * fe_0 + ta_xxz_xxxxx_1[i] + ta1_z_xxz_xxxxx_0[i] * pa_z[i] -
                                ta1_z_xxz_xxxxx_1[i] * pc_z[i];

        ta1_z_xxzz_xxxxy_0[i] = ta1_z_xx_xxxxy_0[i] * fe_0 - ta1_z_xx_xxxxy_1[i] * fe_0 + ta_xxz_xxxxy_1[i] + ta1_z_xxz_xxxxy_0[i] * pa_z[i] -
                                ta1_z_xxz_xxxxy_1[i] * pc_z[i];

        ta1_z_xxzz_xxxxz_0[i] = ta1_z_zz_xxxxz_0[i] * fe_0 - ta1_z_zz_xxxxz_1[i] * fe_0 + 4.0 * ta1_z_xzz_xxxz_0[i] * fe_0 -
                                4.0 * ta1_z_xzz_xxxz_1[i] * fe_0 + ta1_z_xzz_xxxxz_0[i] * pa_x[i] - ta1_z_xzz_xxxxz_1[i] * pc_x[i];

        ta1_z_xxzz_xxxyy_0[i] = ta1_z_xx_xxxyy_0[i] * fe_0 - ta1_z_xx_xxxyy_1[i] * fe_0 + ta_xxz_xxxyy_1[i] + ta1_z_xxz_xxxyy_0[i] * pa_z[i] -
                                ta1_z_xxz_xxxyy_1[i] * pc_z[i];

        ta1_z_xxzz_xxxyz_0[i] = ta1_z_zz_xxxyz_0[i] * fe_0 - ta1_z_zz_xxxyz_1[i] * fe_0 + 3.0 * ta1_z_xzz_xxyz_0[i] * fe_0 -
                                3.0 * ta1_z_xzz_xxyz_1[i] * fe_0 + ta1_z_xzz_xxxyz_0[i] * pa_x[i] - ta1_z_xzz_xxxyz_1[i] * pc_x[i];

        ta1_z_xxzz_xxxzz_0[i] = ta1_z_zz_xxxzz_0[i] * fe_0 - ta1_z_zz_xxxzz_1[i] * fe_0 + 3.0 * ta1_z_xzz_xxzz_0[i] * fe_0 -
                                3.0 * ta1_z_xzz_xxzz_1[i] * fe_0 + ta1_z_xzz_xxxzz_0[i] * pa_x[i] - ta1_z_xzz_xxxzz_1[i] * pc_x[i];

        ta1_z_xxzz_xxyyy_0[i] = ta1_z_xx_xxyyy_0[i] * fe_0 - ta1_z_xx_xxyyy_1[i] * fe_0 + ta_xxz_xxyyy_1[i] + ta1_z_xxz_xxyyy_0[i] * pa_z[i] -
                                ta1_z_xxz_xxyyy_1[i] * pc_z[i];

        ta1_z_xxzz_xxyyz_0[i] = ta1_z_zz_xxyyz_0[i] * fe_0 - ta1_z_zz_xxyyz_1[i] * fe_0 + 2.0 * ta1_z_xzz_xyyz_0[i] * fe_0 -
                                2.0 * ta1_z_xzz_xyyz_1[i] * fe_0 + ta1_z_xzz_xxyyz_0[i] * pa_x[i] - ta1_z_xzz_xxyyz_1[i] * pc_x[i];

        ta1_z_xxzz_xxyzz_0[i] = ta1_z_zz_xxyzz_0[i] * fe_0 - ta1_z_zz_xxyzz_1[i] * fe_0 + 2.0 * ta1_z_xzz_xyzz_0[i] * fe_0 -
                                2.0 * ta1_z_xzz_xyzz_1[i] * fe_0 + ta1_z_xzz_xxyzz_0[i] * pa_x[i] - ta1_z_xzz_xxyzz_1[i] * pc_x[i];

        ta1_z_xxzz_xxzzz_0[i] = ta1_z_zz_xxzzz_0[i] * fe_0 - ta1_z_zz_xxzzz_1[i] * fe_0 + 2.0 * ta1_z_xzz_xzzz_0[i] * fe_0 -
                                2.0 * ta1_z_xzz_xzzz_1[i] * fe_0 + ta1_z_xzz_xxzzz_0[i] * pa_x[i] - ta1_z_xzz_xxzzz_1[i] * pc_x[i];

        ta1_z_xxzz_xyyyy_0[i] = ta1_z_xx_xyyyy_0[i] * fe_0 - ta1_z_xx_xyyyy_1[i] * fe_0 + ta_xxz_xyyyy_1[i] + ta1_z_xxz_xyyyy_0[i] * pa_z[i] -
                                ta1_z_xxz_xyyyy_1[i] * pc_z[i];

        ta1_z_xxzz_xyyyz_0[i] = ta1_z_zz_xyyyz_0[i] * fe_0 - ta1_z_zz_xyyyz_1[i] * fe_0 + ta1_z_xzz_yyyz_0[i] * fe_0 - ta1_z_xzz_yyyz_1[i] * fe_0 +
                                ta1_z_xzz_xyyyz_0[i] * pa_x[i] - ta1_z_xzz_xyyyz_1[i] * pc_x[i];

        ta1_z_xxzz_xyyzz_0[i] = ta1_z_zz_xyyzz_0[i] * fe_0 - ta1_z_zz_xyyzz_1[i] * fe_0 + ta1_z_xzz_yyzz_0[i] * fe_0 - ta1_z_xzz_yyzz_1[i] * fe_0 +
                                ta1_z_xzz_xyyzz_0[i] * pa_x[i] - ta1_z_xzz_xyyzz_1[i] * pc_x[i];

        ta1_z_xxzz_xyzzz_0[i] = ta1_z_zz_xyzzz_0[i] * fe_0 - ta1_z_zz_xyzzz_1[i] * fe_0 + ta1_z_xzz_yzzz_0[i] * fe_0 - ta1_z_xzz_yzzz_1[i] * fe_0 +
                                ta1_z_xzz_xyzzz_0[i] * pa_x[i] - ta1_z_xzz_xyzzz_1[i] * pc_x[i];

        ta1_z_xxzz_xzzzz_0[i] = ta1_z_zz_xzzzz_0[i] * fe_0 - ta1_z_zz_xzzzz_1[i] * fe_0 + ta1_z_xzz_zzzz_0[i] * fe_0 - ta1_z_xzz_zzzz_1[i] * fe_0 +
                                ta1_z_xzz_xzzzz_0[i] * pa_x[i] - ta1_z_xzz_xzzzz_1[i] * pc_x[i];

        ta1_z_xxzz_yyyyy_0[i] =
            ta1_z_zz_yyyyy_0[i] * fe_0 - ta1_z_zz_yyyyy_1[i] * fe_0 + ta1_z_xzz_yyyyy_0[i] * pa_x[i] - ta1_z_xzz_yyyyy_1[i] * pc_x[i];

        ta1_z_xxzz_yyyyz_0[i] =
            ta1_z_zz_yyyyz_0[i] * fe_0 - ta1_z_zz_yyyyz_1[i] * fe_0 + ta1_z_xzz_yyyyz_0[i] * pa_x[i] - ta1_z_xzz_yyyyz_1[i] * pc_x[i];

        ta1_z_xxzz_yyyzz_0[i] =
            ta1_z_zz_yyyzz_0[i] * fe_0 - ta1_z_zz_yyyzz_1[i] * fe_0 + ta1_z_xzz_yyyzz_0[i] * pa_x[i] - ta1_z_xzz_yyyzz_1[i] * pc_x[i];

        ta1_z_xxzz_yyzzz_0[i] =
            ta1_z_zz_yyzzz_0[i] * fe_0 - ta1_z_zz_yyzzz_1[i] * fe_0 + ta1_z_xzz_yyzzz_0[i] * pa_x[i] - ta1_z_xzz_yyzzz_1[i] * pc_x[i];

        ta1_z_xxzz_yzzzz_0[i] =
            ta1_z_zz_yzzzz_0[i] * fe_0 - ta1_z_zz_yzzzz_1[i] * fe_0 + ta1_z_xzz_yzzzz_0[i] * pa_x[i] - ta1_z_xzz_yzzzz_1[i] * pc_x[i];

        ta1_z_xxzz_zzzzz_0[i] =
            ta1_z_zz_zzzzz_0[i] * fe_0 - ta1_z_zz_zzzzz_1[i] * fe_0 + ta1_z_xzz_zzzzz_0[i] * pa_x[i] - ta1_z_xzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 756-777 components of targeted buffer : GH

    auto ta1_z_xyyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 756);

    auto ta1_z_xyyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 757);

    auto ta1_z_xyyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 758);

    auto ta1_z_xyyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 759);

    auto ta1_z_xyyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 760);

    auto ta1_z_xyyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 761);

    auto ta1_z_xyyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 762);

    auto ta1_z_xyyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 763);

    auto ta1_z_xyyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 764);

    auto ta1_z_xyyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 765);

    auto ta1_z_xyyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 766);

    auto ta1_z_xyyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 767);

    auto ta1_z_xyyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 768);

    auto ta1_z_xyyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 769);

    auto ta1_z_xyyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 770);

    auto ta1_z_xyyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 771);

    auto ta1_z_xyyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 772);

    auto ta1_z_xyyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 773);

    auto ta1_z_xyyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 774);

    auto ta1_z_xyyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 775);

    auto ta1_z_xyyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 776);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_xyyy_xxxxx_0, \
                             ta1_z_xyyy_xxxxy_0, \
                             ta1_z_xyyy_xxxxz_0, \
                             ta1_z_xyyy_xxxyy_0, \
                             ta1_z_xyyy_xxxyz_0, \
                             ta1_z_xyyy_xxxzz_0, \
                             ta1_z_xyyy_xxyyy_0, \
                             ta1_z_xyyy_xxyyz_0, \
                             ta1_z_xyyy_xxyzz_0, \
                             ta1_z_xyyy_xxzzz_0, \
                             ta1_z_xyyy_xyyyy_0, \
                             ta1_z_xyyy_xyyyz_0, \
                             ta1_z_xyyy_xyyzz_0, \
                             ta1_z_xyyy_xyzzz_0, \
                             ta1_z_xyyy_xzzzz_0, \
                             ta1_z_xyyy_yyyyy_0, \
                             ta1_z_xyyy_yyyyz_0, \
                             ta1_z_xyyy_yyyzz_0, \
                             ta1_z_xyyy_yyzzz_0, \
                             ta1_z_xyyy_yzzzz_0, \
                             ta1_z_xyyy_zzzzz_0, \
                             ta1_z_yyy_xxxx_0,   \
                             ta1_z_yyy_xxxx_1,   \
                             ta1_z_yyy_xxxxx_0,  \
                             ta1_z_yyy_xxxxx_1,  \
                             ta1_z_yyy_xxxxy_0,  \
                             ta1_z_yyy_xxxxy_1,  \
                             ta1_z_yyy_xxxxz_0,  \
                             ta1_z_yyy_xxxxz_1,  \
                             ta1_z_yyy_xxxy_0,   \
                             ta1_z_yyy_xxxy_1,   \
                             ta1_z_yyy_xxxyy_0,  \
                             ta1_z_yyy_xxxyy_1,  \
                             ta1_z_yyy_xxxyz_0,  \
                             ta1_z_yyy_xxxyz_1,  \
                             ta1_z_yyy_xxxz_0,   \
                             ta1_z_yyy_xxxz_1,   \
                             ta1_z_yyy_xxxzz_0,  \
                             ta1_z_yyy_xxxzz_1,  \
                             ta1_z_yyy_xxyy_0,   \
                             ta1_z_yyy_xxyy_1,   \
                             ta1_z_yyy_xxyyy_0,  \
                             ta1_z_yyy_xxyyy_1,  \
                             ta1_z_yyy_xxyyz_0,  \
                             ta1_z_yyy_xxyyz_1,  \
                             ta1_z_yyy_xxyz_0,   \
                             ta1_z_yyy_xxyz_1,   \
                             ta1_z_yyy_xxyzz_0,  \
                             ta1_z_yyy_xxyzz_1,  \
                             ta1_z_yyy_xxzz_0,   \
                             ta1_z_yyy_xxzz_1,   \
                             ta1_z_yyy_xxzzz_0,  \
                             ta1_z_yyy_xxzzz_1,  \
                             ta1_z_yyy_xyyy_0,   \
                             ta1_z_yyy_xyyy_1,   \
                             ta1_z_yyy_xyyyy_0,  \
                             ta1_z_yyy_xyyyy_1,  \
                             ta1_z_yyy_xyyyz_0,  \
                             ta1_z_yyy_xyyyz_1,  \
                             ta1_z_yyy_xyyz_0,   \
                             ta1_z_yyy_xyyz_1,   \
                             ta1_z_yyy_xyyzz_0,  \
                             ta1_z_yyy_xyyzz_1,  \
                             ta1_z_yyy_xyzz_0,   \
                             ta1_z_yyy_xyzz_1,   \
                             ta1_z_yyy_xyzzz_0,  \
                             ta1_z_yyy_xyzzz_1,  \
                             ta1_z_yyy_xzzz_0,   \
                             ta1_z_yyy_xzzz_1,   \
                             ta1_z_yyy_xzzzz_0,  \
                             ta1_z_yyy_xzzzz_1,  \
                             ta1_z_yyy_yyyy_0,   \
                             ta1_z_yyy_yyyy_1,   \
                             ta1_z_yyy_yyyyy_0,  \
                             ta1_z_yyy_yyyyy_1,  \
                             ta1_z_yyy_yyyyz_0,  \
                             ta1_z_yyy_yyyyz_1,  \
                             ta1_z_yyy_yyyz_0,   \
                             ta1_z_yyy_yyyz_1,   \
                             ta1_z_yyy_yyyzz_0,  \
                             ta1_z_yyy_yyyzz_1,  \
                             ta1_z_yyy_yyzz_0,   \
                             ta1_z_yyy_yyzz_1,   \
                             ta1_z_yyy_yyzzz_0,  \
                             ta1_z_yyy_yyzzz_1,  \
                             ta1_z_yyy_yzzz_0,   \
                             ta1_z_yyy_yzzz_1,   \
                             ta1_z_yyy_yzzzz_0,  \
                             ta1_z_yyy_yzzzz_1,  \
                             ta1_z_yyy_zzzz_0,   \
                             ta1_z_yyy_zzzz_1,   \
                             ta1_z_yyy_zzzzz_0,  \
                             ta1_z_yyy_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyy_xxxxx_0[i] =
            5.0 * ta1_z_yyy_xxxx_0[i] * fe_0 - 5.0 * ta1_z_yyy_xxxx_1[i] * fe_0 + ta1_z_yyy_xxxxx_0[i] * pa_x[i] - ta1_z_yyy_xxxxx_1[i] * pc_x[i];

        ta1_z_xyyy_xxxxy_0[i] =
            4.0 * ta1_z_yyy_xxxy_0[i] * fe_0 - 4.0 * ta1_z_yyy_xxxy_1[i] * fe_0 + ta1_z_yyy_xxxxy_0[i] * pa_x[i] - ta1_z_yyy_xxxxy_1[i] * pc_x[i];

        ta1_z_xyyy_xxxxz_0[i] =
            4.0 * ta1_z_yyy_xxxz_0[i] * fe_0 - 4.0 * ta1_z_yyy_xxxz_1[i] * fe_0 + ta1_z_yyy_xxxxz_0[i] * pa_x[i] - ta1_z_yyy_xxxxz_1[i] * pc_x[i];

        ta1_z_xyyy_xxxyy_0[i] =
            3.0 * ta1_z_yyy_xxyy_0[i] * fe_0 - 3.0 * ta1_z_yyy_xxyy_1[i] * fe_0 + ta1_z_yyy_xxxyy_0[i] * pa_x[i] - ta1_z_yyy_xxxyy_1[i] * pc_x[i];

        ta1_z_xyyy_xxxyz_0[i] =
            3.0 * ta1_z_yyy_xxyz_0[i] * fe_0 - 3.0 * ta1_z_yyy_xxyz_1[i] * fe_0 + ta1_z_yyy_xxxyz_0[i] * pa_x[i] - ta1_z_yyy_xxxyz_1[i] * pc_x[i];

        ta1_z_xyyy_xxxzz_0[i] =
            3.0 * ta1_z_yyy_xxzz_0[i] * fe_0 - 3.0 * ta1_z_yyy_xxzz_1[i] * fe_0 + ta1_z_yyy_xxxzz_0[i] * pa_x[i] - ta1_z_yyy_xxxzz_1[i] * pc_x[i];

        ta1_z_xyyy_xxyyy_0[i] =
            2.0 * ta1_z_yyy_xyyy_0[i] * fe_0 - 2.0 * ta1_z_yyy_xyyy_1[i] * fe_0 + ta1_z_yyy_xxyyy_0[i] * pa_x[i] - ta1_z_yyy_xxyyy_1[i] * pc_x[i];

        ta1_z_xyyy_xxyyz_0[i] =
            2.0 * ta1_z_yyy_xyyz_0[i] * fe_0 - 2.0 * ta1_z_yyy_xyyz_1[i] * fe_0 + ta1_z_yyy_xxyyz_0[i] * pa_x[i] - ta1_z_yyy_xxyyz_1[i] * pc_x[i];

        ta1_z_xyyy_xxyzz_0[i] =
            2.0 * ta1_z_yyy_xyzz_0[i] * fe_0 - 2.0 * ta1_z_yyy_xyzz_1[i] * fe_0 + ta1_z_yyy_xxyzz_0[i] * pa_x[i] - ta1_z_yyy_xxyzz_1[i] * pc_x[i];

        ta1_z_xyyy_xxzzz_0[i] =
            2.0 * ta1_z_yyy_xzzz_0[i] * fe_0 - 2.0 * ta1_z_yyy_xzzz_1[i] * fe_0 + ta1_z_yyy_xxzzz_0[i] * pa_x[i] - ta1_z_yyy_xxzzz_1[i] * pc_x[i];

        ta1_z_xyyy_xyyyy_0[i] =
            ta1_z_yyy_yyyy_0[i] * fe_0 - ta1_z_yyy_yyyy_1[i] * fe_0 + ta1_z_yyy_xyyyy_0[i] * pa_x[i] - ta1_z_yyy_xyyyy_1[i] * pc_x[i];

        ta1_z_xyyy_xyyyz_0[i] =
            ta1_z_yyy_yyyz_0[i] * fe_0 - ta1_z_yyy_yyyz_1[i] * fe_0 + ta1_z_yyy_xyyyz_0[i] * pa_x[i] - ta1_z_yyy_xyyyz_1[i] * pc_x[i];

        ta1_z_xyyy_xyyzz_0[i] =
            ta1_z_yyy_yyzz_0[i] * fe_0 - ta1_z_yyy_yyzz_1[i] * fe_0 + ta1_z_yyy_xyyzz_0[i] * pa_x[i] - ta1_z_yyy_xyyzz_1[i] * pc_x[i];

        ta1_z_xyyy_xyzzz_0[i] =
            ta1_z_yyy_yzzz_0[i] * fe_0 - ta1_z_yyy_yzzz_1[i] * fe_0 + ta1_z_yyy_xyzzz_0[i] * pa_x[i] - ta1_z_yyy_xyzzz_1[i] * pc_x[i];

        ta1_z_xyyy_xzzzz_0[i] =
            ta1_z_yyy_zzzz_0[i] * fe_0 - ta1_z_yyy_zzzz_1[i] * fe_0 + ta1_z_yyy_xzzzz_0[i] * pa_x[i] - ta1_z_yyy_xzzzz_1[i] * pc_x[i];

        ta1_z_xyyy_yyyyy_0[i] = ta1_z_yyy_yyyyy_0[i] * pa_x[i] - ta1_z_yyy_yyyyy_1[i] * pc_x[i];

        ta1_z_xyyy_yyyyz_0[i] = ta1_z_yyy_yyyyz_0[i] * pa_x[i] - ta1_z_yyy_yyyyz_1[i] * pc_x[i];

        ta1_z_xyyy_yyyzz_0[i] = ta1_z_yyy_yyyzz_0[i] * pa_x[i] - ta1_z_yyy_yyyzz_1[i] * pc_x[i];

        ta1_z_xyyy_yyzzz_0[i] = ta1_z_yyy_yyzzz_0[i] * pa_x[i] - ta1_z_yyy_yyzzz_1[i] * pc_x[i];

        ta1_z_xyyy_yzzzz_0[i] = ta1_z_yyy_yzzzz_0[i] * pa_x[i] - ta1_z_yyy_yzzzz_1[i] * pc_x[i];

        ta1_z_xyyy_zzzzz_0[i] = ta1_z_yyy_zzzzz_0[i] * pa_x[i] - ta1_z_yyy_zzzzz_1[i] * pc_x[i];
    }

    // Set up 777-798 components of targeted buffer : GH

    auto ta1_z_xyyz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 777);

    auto ta1_z_xyyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 778);

    auto ta1_z_xyyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 779);

    auto ta1_z_xyyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 780);

    auto ta1_z_xyyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 781);

    auto ta1_z_xyyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 782);

    auto ta1_z_xyyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 783);

    auto ta1_z_xyyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 784);

    auto ta1_z_xyyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 785);

    auto ta1_z_xyyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 786);

    auto ta1_z_xyyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 787);

    auto ta1_z_xyyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 788);

    auto ta1_z_xyyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 789);

    auto ta1_z_xyyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 790);

    auto ta1_z_xyyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 791);

    auto ta1_z_xyyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 792);

    auto ta1_z_xyyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 793);

    auto ta1_z_xyyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 794);

    auto ta1_z_xyyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 795);

    auto ta1_z_xyyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 796);

    auto ta1_z_xyyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 797);

#pragma omp simd aligned(pa_x,                   \
                             pa_z,               \
                             pc_x,               \
                             pc_z,               \
                             ta1_z_xyy_xxxxx_0,  \
                             ta1_z_xyy_xxxxx_1,  \
                             ta1_z_xyy_xxxxy_0,  \
                             ta1_z_xyy_xxxxy_1,  \
                             ta1_z_xyy_xxxyy_0,  \
                             ta1_z_xyy_xxxyy_1,  \
                             ta1_z_xyy_xxyyy_0,  \
                             ta1_z_xyy_xxyyy_1,  \
                             ta1_z_xyy_xyyyy_0,  \
                             ta1_z_xyy_xyyyy_1,  \
                             ta1_z_xyyz_xxxxx_0, \
                             ta1_z_xyyz_xxxxy_0, \
                             ta1_z_xyyz_xxxxz_0, \
                             ta1_z_xyyz_xxxyy_0, \
                             ta1_z_xyyz_xxxyz_0, \
                             ta1_z_xyyz_xxxzz_0, \
                             ta1_z_xyyz_xxyyy_0, \
                             ta1_z_xyyz_xxyyz_0, \
                             ta1_z_xyyz_xxyzz_0, \
                             ta1_z_xyyz_xxzzz_0, \
                             ta1_z_xyyz_xyyyy_0, \
                             ta1_z_xyyz_xyyyz_0, \
                             ta1_z_xyyz_xyyzz_0, \
                             ta1_z_xyyz_xyzzz_0, \
                             ta1_z_xyyz_xzzzz_0, \
                             ta1_z_xyyz_yyyyy_0, \
                             ta1_z_xyyz_yyyyz_0, \
                             ta1_z_xyyz_yyyzz_0, \
                             ta1_z_xyyz_yyzzz_0, \
                             ta1_z_xyyz_yzzzz_0, \
                             ta1_z_xyyz_zzzzz_0, \
                             ta1_z_yyz_xxxxz_0,  \
                             ta1_z_yyz_xxxxz_1,  \
                             ta1_z_yyz_xxxyz_0,  \
                             ta1_z_yyz_xxxyz_1,  \
                             ta1_z_yyz_xxxz_0,   \
                             ta1_z_yyz_xxxz_1,   \
                             ta1_z_yyz_xxxzz_0,  \
                             ta1_z_yyz_xxxzz_1,  \
                             ta1_z_yyz_xxyyz_0,  \
                             ta1_z_yyz_xxyyz_1,  \
                             ta1_z_yyz_xxyz_0,   \
                             ta1_z_yyz_xxyz_1,   \
                             ta1_z_yyz_xxyzz_0,  \
                             ta1_z_yyz_xxyzz_1,  \
                             ta1_z_yyz_xxzz_0,   \
                             ta1_z_yyz_xxzz_1,   \
                             ta1_z_yyz_xxzzz_0,  \
                             ta1_z_yyz_xxzzz_1,  \
                             ta1_z_yyz_xyyyz_0,  \
                             ta1_z_yyz_xyyyz_1,  \
                             ta1_z_yyz_xyyz_0,   \
                             ta1_z_yyz_xyyz_1,   \
                             ta1_z_yyz_xyyzz_0,  \
                             ta1_z_yyz_xyyzz_1,  \
                             ta1_z_yyz_xyzz_0,   \
                             ta1_z_yyz_xyzz_1,   \
                             ta1_z_yyz_xyzzz_0,  \
                             ta1_z_yyz_xyzzz_1,  \
                             ta1_z_yyz_xzzz_0,   \
                             ta1_z_yyz_xzzz_1,   \
                             ta1_z_yyz_xzzzz_0,  \
                             ta1_z_yyz_xzzzz_1,  \
                             ta1_z_yyz_yyyyy_0,  \
                             ta1_z_yyz_yyyyy_1,  \
                             ta1_z_yyz_yyyyz_0,  \
                             ta1_z_yyz_yyyyz_1,  \
                             ta1_z_yyz_yyyz_0,   \
                             ta1_z_yyz_yyyz_1,   \
                             ta1_z_yyz_yyyzz_0,  \
                             ta1_z_yyz_yyyzz_1,  \
                             ta1_z_yyz_yyzz_0,   \
                             ta1_z_yyz_yyzz_1,   \
                             ta1_z_yyz_yyzzz_0,  \
                             ta1_z_yyz_yyzzz_1,  \
                             ta1_z_yyz_yzzz_0,   \
                             ta1_z_yyz_yzzz_1,   \
                             ta1_z_yyz_yzzzz_0,  \
                             ta1_z_yyz_yzzzz_1,  \
                             ta1_z_yyz_zzzz_0,   \
                             ta1_z_yyz_zzzz_1,   \
                             ta1_z_yyz_zzzzz_0,  \
                             ta1_z_yyz_zzzzz_1,  \
                             ta_xyy_xxxxx_1,     \
                             ta_xyy_xxxxy_1,     \
                             ta_xyy_xxxyy_1,     \
                             ta_xyy_xxyyy_1,     \
                             ta_xyy_xyyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyyz_xxxxx_0[i] = ta_xyy_xxxxx_1[i] + ta1_z_xyy_xxxxx_0[i] * pa_z[i] - ta1_z_xyy_xxxxx_1[i] * pc_z[i];

        ta1_z_xyyz_xxxxy_0[i] = ta_xyy_xxxxy_1[i] + ta1_z_xyy_xxxxy_0[i] * pa_z[i] - ta1_z_xyy_xxxxy_1[i] * pc_z[i];

        ta1_z_xyyz_xxxxz_0[i] =
            4.0 * ta1_z_yyz_xxxz_0[i] * fe_0 - 4.0 * ta1_z_yyz_xxxz_1[i] * fe_0 + ta1_z_yyz_xxxxz_0[i] * pa_x[i] - ta1_z_yyz_xxxxz_1[i] * pc_x[i];

        ta1_z_xyyz_xxxyy_0[i] = ta_xyy_xxxyy_1[i] + ta1_z_xyy_xxxyy_0[i] * pa_z[i] - ta1_z_xyy_xxxyy_1[i] * pc_z[i];

        ta1_z_xyyz_xxxyz_0[i] =
            3.0 * ta1_z_yyz_xxyz_0[i] * fe_0 - 3.0 * ta1_z_yyz_xxyz_1[i] * fe_0 + ta1_z_yyz_xxxyz_0[i] * pa_x[i] - ta1_z_yyz_xxxyz_1[i] * pc_x[i];

        ta1_z_xyyz_xxxzz_0[i] =
            3.0 * ta1_z_yyz_xxzz_0[i] * fe_0 - 3.0 * ta1_z_yyz_xxzz_1[i] * fe_0 + ta1_z_yyz_xxxzz_0[i] * pa_x[i] - ta1_z_yyz_xxxzz_1[i] * pc_x[i];

        ta1_z_xyyz_xxyyy_0[i] = ta_xyy_xxyyy_1[i] + ta1_z_xyy_xxyyy_0[i] * pa_z[i] - ta1_z_xyy_xxyyy_1[i] * pc_z[i];

        ta1_z_xyyz_xxyyz_0[i] =
            2.0 * ta1_z_yyz_xyyz_0[i] * fe_0 - 2.0 * ta1_z_yyz_xyyz_1[i] * fe_0 + ta1_z_yyz_xxyyz_0[i] * pa_x[i] - ta1_z_yyz_xxyyz_1[i] * pc_x[i];

        ta1_z_xyyz_xxyzz_0[i] =
            2.0 * ta1_z_yyz_xyzz_0[i] * fe_0 - 2.0 * ta1_z_yyz_xyzz_1[i] * fe_0 + ta1_z_yyz_xxyzz_0[i] * pa_x[i] - ta1_z_yyz_xxyzz_1[i] * pc_x[i];

        ta1_z_xyyz_xxzzz_0[i] =
            2.0 * ta1_z_yyz_xzzz_0[i] * fe_0 - 2.0 * ta1_z_yyz_xzzz_1[i] * fe_0 + ta1_z_yyz_xxzzz_0[i] * pa_x[i] - ta1_z_yyz_xxzzz_1[i] * pc_x[i];

        ta1_z_xyyz_xyyyy_0[i] = ta_xyy_xyyyy_1[i] + ta1_z_xyy_xyyyy_0[i] * pa_z[i] - ta1_z_xyy_xyyyy_1[i] * pc_z[i];

        ta1_z_xyyz_xyyyz_0[i] =
            ta1_z_yyz_yyyz_0[i] * fe_0 - ta1_z_yyz_yyyz_1[i] * fe_0 + ta1_z_yyz_xyyyz_0[i] * pa_x[i] - ta1_z_yyz_xyyyz_1[i] * pc_x[i];

        ta1_z_xyyz_xyyzz_0[i] =
            ta1_z_yyz_yyzz_0[i] * fe_0 - ta1_z_yyz_yyzz_1[i] * fe_0 + ta1_z_yyz_xyyzz_0[i] * pa_x[i] - ta1_z_yyz_xyyzz_1[i] * pc_x[i];

        ta1_z_xyyz_xyzzz_0[i] =
            ta1_z_yyz_yzzz_0[i] * fe_0 - ta1_z_yyz_yzzz_1[i] * fe_0 + ta1_z_yyz_xyzzz_0[i] * pa_x[i] - ta1_z_yyz_xyzzz_1[i] * pc_x[i];

        ta1_z_xyyz_xzzzz_0[i] =
            ta1_z_yyz_zzzz_0[i] * fe_0 - ta1_z_yyz_zzzz_1[i] * fe_0 + ta1_z_yyz_xzzzz_0[i] * pa_x[i] - ta1_z_yyz_xzzzz_1[i] * pc_x[i];

        ta1_z_xyyz_yyyyy_0[i] = ta1_z_yyz_yyyyy_0[i] * pa_x[i] - ta1_z_yyz_yyyyy_1[i] * pc_x[i];

        ta1_z_xyyz_yyyyz_0[i] = ta1_z_yyz_yyyyz_0[i] * pa_x[i] - ta1_z_yyz_yyyyz_1[i] * pc_x[i];

        ta1_z_xyyz_yyyzz_0[i] = ta1_z_yyz_yyyzz_0[i] * pa_x[i] - ta1_z_yyz_yyyzz_1[i] * pc_x[i];

        ta1_z_xyyz_yyzzz_0[i] = ta1_z_yyz_yyzzz_0[i] * pa_x[i] - ta1_z_yyz_yyzzz_1[i] * pc_x[i];

        ta1_z_xyyz_yzzzz_0[i] = ta1_z_yyz_yzzzz_0[i] * pa_x[i] - ta1_z_yyz_yzzzz_1[i] * pc_x[i];

        ta1_z_xyyz_zzzzz_0[i] = ta1_z_yyz_zzzzz_0[i] * pa_x[i] - ta1_z_yyz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 798-819 components of targeted buffer : GH

    auto ta1_z_xyzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 798);

    auto ta1_z_xyzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 799);

    auto ta1_z_xyzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 800);

    auto ta1_z_xyzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 801);

    auto ta1_z_xyzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 802);

    auto ta1_z_xyzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 803);

    auto ta1_z_xyzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 804);

    auto ta1_z_xyzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 805);

    auto ta1_z_xyzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 806);

    auto ta1_z_xyzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 807);

    auto ta1_z_xyzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 808);

    auto ta1_z_xyzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 809);

    auto ta1_z_xyzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 810);

    auto ta1_z_xyzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 811);

    auto ta1_z_xyzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 812);

    auto ta1_z_xyzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 813);

    auto ta1_z_xyzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 814);

    auto ta1_z_xyzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 815);

    auto ta1_z_xyzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 816);

    auto ta1_z_xyzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 817);

    auto ta1_z_xyzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 818);

#pragma omp simd aligned(pa_x,                   \
                             pa_y,               \
                             pc_x,               \
                             pc_y,               \
                             ta1_z_xyzz_xxxxx_0, \
                             ta1_z_xyzz_xxxxy_0, \
                             ta1_z_xyzz_xxxxz_0, \
                             ta1_z_xyzz_xxxyy_0, \
                             ta1_z_xyzz_xxxyz_0, \
                             ta1_z_xyzz_xxxzz_0, \
                             ta1_z_xyzz_xxyyy_0, \
                             ta1_z_xyzz_xxyyz_0, \
                             ta1_z_xyzz_xxyzz_0, \
                             ta1_z_xyzz_xxzzz_0, \
                             ta1_z_xyzz_xyyyy_0, \
                             ta1_z_xyzz_xyyyz_0, \
                             ta1_z_xyzz_xyyzz_0, \
                             ta1_z_xyzz_xyzzz_0, \
                             ta1_z_xyzz_xzzzz_0, \
                             ta1_z_xyzz_yyyyy_0, \
                             ta1_z_xyzz_yyyyz_0, \
                             ta1_z_xyzz_yyyzz_0, \
                             ta1_z_xyzz_yyzzz_0, \
                             ta1_z_xyzz_yzzzz_0, \
                             ta1_z_xyzz_zzzzz_0, \
                             ta1_z_xzz_xxxxx_0,  \
                             ta1_z_xzz_xxxxx_1,  \
                             ta1_z_xzz_xxxxz_0,  \
                             ta1_z_xzz_xxxxz_1,  \
                             ta1_z_xzz_xxxzz_0,  \
                             ta1_z_xzz_xxxzz_1,  \
                             ta1_z_xzz_xxzzz_0,  \
                             ta1_z_xzz_xxzzz_1,  \
                             ta1_z_xzz_xzzzz_0,  \
                             ta1_z_xzz_xzzzz_1,  \
                             ta1_z_yzz_xxxxy_0,  \
                             ta1_z_yzz_xxxxy_1,  \
                             ta1_z_yzz_xxxy_0,   \
                             ta1_z_yzz_xxxy_1,   \
                             ta1_z_yzz_xxxyy_0,  \
                             ta1_z_yzz_xxxyy_1,  \
                             ta1_z_yzz_xxxyz_0,  \
                             ta1_z_yzz_xxxyz_1,  \
                             ta1_z_yzz_xxyy_0,   \
                             ta1_z_yzz_xxyy_1,   \
                             ta1_z_yzz_xxyyy_0,  \
                             ta1_z_yzz_xxyyy_1,  \
                             ta1_z_yzz_xxyyz_0,  \
                             ta1_z_yzz_xxyyz_1,  \
                             ta1_z_yzz_xxyz_0,   \
                             ta1_z_yzz_xxyz_1,   \
                             ta1_z_yzz_xxyzz_0,  \
                             ta1_z_yzz_xxyzz_1,  \
                             ta1_z_yzz_xyyy_0,   \
                             ta1_z_yzz_xyyy_1,   \
                             ta1_z_yzz_xyyyy_0,  \
                             ta1_z_yzz_xyyyy_1,  \
                             ta1_z_yzz_xyyyz_0,  \
                             ta1_z_yzz_xyyyz_1,  \
                             ta1_z_yzz_xyyz_0,   \
                             ta1_z_yzz_xyyz_1,   \
                             ta1_z_yzz_xyyzz_0,  \
                             ta1_z_yzz_xyyzz_1,  \
                             ta1_z_yzz_xyzz_0,   \
                             ta1_z_yzz_xyzz_1,   \
                             ta1_z_yzz_xyzzz_0,  \
                             ta1_z_yzz_xyzzz_1,  \
                             ta1_z_yzz_yyyy_0,   \
                             ta1_z_yzz_yyyy_1,   \
                             ta1_z_yzz_yyyyy_0,  \
                             ta1_z_yzz_yyyyy_1,  \
                             ta1_z_yzz_yyyyz_0,  \
                             ta1_z_yzz_yyyyz_1,  \
                             ta1_z_yzz_yyyz_0,   \
                             ta1_z_yzz_yyyz_1,   \
                             ta1_z_yzz_yyyzz_0,  \
                             ta1_z_yzz_yyyzz_1,  \
                             ta1_z_yzz_yyzz_0,   \
                             ta1_z_yzz_yyzz_1,   \
                             ta1_z_yzz_yyzzz_0,  \
                             ta1_z_yzz_yyzzz_1,  \
                             ta1_z_yzz_yzzz_0,   \
                             ta1_z_yzz_yzzz_1,   \
                             ta1_z_yzz_yzzzz_0,  \
                             ta1_z_yzz_yzzzz_1,  \
                             ta1_z_yzz_zzzzz_0,  \
                             ta1_z_yzz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xyzz_xxxxx_0[i] = ta1_z_xzz_xxxxx_0[i] * pa_y[i] - ta1_z_xzz_xxxxx_1[i] * pc_y[i];

        ta1_z_xyzz_xxxxy_0[i] =
            4.0 * ta1_z_yzz_xxxy_0[i] * fe_0 - 4.0 * ta1_z_yzz_xxxy_1[i] * fe_0 + ta1_z_yzz_xxxxy_0[i] * pa_x[i] - ta1_z_yzz_xxxxy_1[i] * pc_x[i];

        ta1_z_xyzz_xxxxz_0[i] = ta1_z_xzz_xxxxz_0[i] * pa_y[i] - ta1_z_xzz_xxxxz_1[i] * pc_y[i];

        ta1_z_xyzz_xxxyy_0[i] =
            3.0 * ta1_z_yzz_xxyy_0[i] * fe_0 - 3.0 * ta1_z_yzz_xxyy_1[i] * fe_0 + ta1_z_yzz_xxxyy_0[i] * pa_x[i] - ta1_z_yzz_xxxyy_1[i] * pc_x[i];

        ta1_z_xyzz_xxxyz_0[i] =
            3.0 * ta1_z_yzz_xxyz_0[i] * fe_0 - 3.0 * ta1_z_yzz_xxyz_1[i] * fe_0 + ta1_z_yzz_xxxyz_0[i] * pa_x[i] - ta1_z_yzz_xxxyz_1[i] * pc_x[i];

        ta1_z_xyzz_xxxzz_0[i] = ta1_z_xzz_xxxzz_0[i] * pa_y[i] - ta1_z_xzz_xxxzz_1[i] * pc_y[i];

        ta1_z_xyzz_xxyyy_0[i] =
            2.0 * ta1_z_yzz_xyyy_0[i] * fe_0 - 2.0 * ta1_z_yzz_xyyy_1[i] * fe_0 + ta1_z_yzz_xxyyy_0[i] * pa_x[i] - ta1_z_yzz_xxyyy_1[i] * pc_x[i];

        ta1_z_xyzz_xxyyz_0[i] =
            2.0 * ta1_z_yzz_xyyz_0[i] * fe_0 - 2.0 * ta1_z_yzz_xyyz_1[i] * fe_0 + ta1_z_yzz_xxyyz_0[i] * pa_x[i] - ta1_z_yzz_xxyyz_1[i] * pc_x[i];

        ta1_z_xyzz_xxyzz_0[i] =
            2.0 * ta1_z_yzz_xyzz_0[i] * fe_0 - 2.0 * ta1_z_yzz_xyzz_1[i] * fe_0 + ta1_z_yzz_xxyzz_0[i] * pa_x[i] - ta1_z_yzz_xxyzz_1[i] * pc_x[i];

        ta1_z_xyzz_xxzzz_0[i] = ta1_z_xzz_xxzzz_0[i] * pa_y[i] - ta1_z_xzz_xxzzz_1[i] * pc_y[i];

        ta1_z_xyzz_xyyyy_0[i] =
            ta1_z_yzz_yyyy_0[i] * fe_0 - ta1_z_yzz_yyyy_1[i] * fe_0 + ta1_z_yzz_xyyyy_0[i] * pa_x[i] - ta1_z_yzz_xyyyy_1[i] * pc_x[i];

        ta1_z_xyzz_xyyyz_0[i] =
            ta1_z_yzz_yyyz_0[i] * fe_0 - ta1_z_yzz_yyyz_1[i] * fe_0 + ta1_z_yzz_xyyyz_0[i] * pa_x[i] - ta1_z_yzz_xyyyz_1[i] * pc_x[i];

        ta1_z_xyzz_xyyzz_0[i] =
            ta1_z_yzz_yyzz_0[i] * fe_0 - ta1_z_yzz_yyzz_1[i] * fe_0 + ta1_z_yzz_xyyzz_0[i] * pa_x[i] - ta1_z_yzz_xyyzz_1[i] * pc_x[i];

        ta1_z_xyzz_xyzzz_0[i] =
            ta1_z_yzz_yzzz_0[i] * fe_0 - ta1_z_yzz_yzzz_1[i] * fe_0 + ta1_z_yzz_xyzzz_0[i] * pa_x[i] - ta1_z_yzz_xyzzz_1[i] * pc_x[i];

        ta1_z_xyzz_xzzzz_0[i] = ta1_z_xzz_xzzzz_0[i] * pa_y[i] - ta1_z_xzz_xzzzz_1[i] * pc_y[i];

        ta1_z_xyzz_yyyyy_0[i] = ta1_z_yzz_yyyyy_0[i] * pa_x[i] - ta1_z_yzz_yyyyy_1[i] * pc_x[i];

        ta1_z_xyzz_yyyyz_0[i] = ta1_z_yzz_yyyyz_0[i] * pa_x[i] - ta1_z_yzz_yyyyz_1[i] * pc_x[i];

        ta1_z_xyzz_yyyzz_0[i] = ta1_z_yzz_yyyzz_0[i] * pa_x[i] - ta1_z_yzz_yyyzz_1[i] * pc_x[i];

        ta1_z_xyzz_yyzzz_0[i] = ta1_z_yzz_yyzzz_0[i] * pa_x[i] - ta1_z_yzz_yyzzz_1[i] * pc_x[i];

        ta1_z_xyzz_yzzzz_0[i] = ta1_z_yzz_yzzzz_0[i] * pa_x[i] - ta1_z_yzz_yzzzz_1[i] * pc_x[i];

        ta1_z_xyzz_zzzzz_0[i] = ta1_z_yzz_zzzzz_0[i] * pa_x[i] - ta1_z_yzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 819-840 components of targeted buffer : GH

    auto ta1_z_xzzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 819);

    auto ta1_z_xzzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 820);

    auto ta1_z_xzzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 821);

    auto ta1_z_xzzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 822);

    auto ta1_z_xzzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 823);

    auto ta1_z_xzzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 824);

    auto ta1_z_xzzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 825);

    auto ta1_z_xzzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 826);

    auto ta1_z_xzzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 827);

    auto ta1_z_xzzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 828);

    auto ta1_z_xzzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 829);

    auto ta1_z_xzzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 830);

    auto ta1_z_xzzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 831);

    auto ta1_z_xzzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 832);

    auto ta1_z_xzzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 833);

    auto ta1_z_xzzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 834);

    auto ta1_z_xzzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 835);

    auto ta1_z_xzzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 836);

    auto ta1_z_xzzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 837);

    auto ta1_z_xzzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 838);

    auto ta1_z_xzzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 839);

#pragma omp simd aligned(pa_x,                   \
                             pc_x,               \
                             ta1_z_xzzz_xxxxx_0, \
                             ta1_z_xzzz_xxxxy_0, \
                             ta1_z_xzzz_xxxxz_0, \
                             ta1_z_xzzz_xxxyy_0, \
                             ta1_z_xzzz_xxxyz_0, \
                             ta1_z_xzzz_xxxzz_0, \
                             ta1_z_xzzz_xxyyy_0, \
                             ta1_z_xzzz_xxyyz_0, \
                             ta1_z_xzzz_xxyzz_0, \
                             ta1_z_xzzz_xxzzz_0, \
                             ta1_z_xzzz_xyyyy_0, \
                             ta1_z_xzzz_xyyyz_0, \
                             ta1_z_xzzz_xyyzz_0, \
                             ta1_z_xzzz_xyzzz_0, \
                             ta1_z_xzzz_xzzzz_0, \
                             ta1_z_xzzz_yyyyy_0, \
                             ta1_z_xzzz_yyyyz_0, \
                             ta1_z_xzzz_yyyzz_0, \
                             ta1_z_xzzz_yyzzz_0, \
                             ta1_z_xzzz_yzzzz_0, \
                             ta1_z_xzzz_zzzzz_0, \
                             ta1_z_zzz_xxxx_0,   \
                             ta1_z_zzz_xxxx_1,   \
                             ta1_z_zzz_xxxxx_0,  \
                             ta1_z_zzz_xxxxx_1,  \
                             ta1_z_zzz_xxxxy_0,  \
                             ta1_z_zzz_xxxxy_1,  \
                             ta1_z_zzz_xxxxz_0,  \
                             ta1_z_zzz_xxxxz_1,  \
                             ta1_z_zzz_xxxy_0,   \
                             ta1_z_zzz_xxxy_1,   \
                             ta1_z_zzz_xxxyy_0,  \
                             ta1_z_zzz_xxxyy_1,  \
                             ta1_z_zzz_xxxyz_0,  \
                             ta1_z_zzz_xxxyz_1,  \
                             ta1_z_zzz_xxxz_0,   \
                             ta1_z_zzz_xxxz_1,   \
                             ta1_z_zzz_xxxzz_0,  \
                             ta1_z_zzz_xxxzz_1,  \
                             ta1_z_zzz_xxyy_0,   \
                             ta1_z_zzz_xxyy_1,   \
                             ta1_z_zzz_xxyyy_0,  \
                             ta1_z_zzz_xxyyy_1,  \
                             ta1_z_zzz_xxyyz_0,  \
                             ta1_z_zzz_xxyyz_1,  \
                             ta1_z_zzz_xxyz_0,   \
                             ta1_z_zzz_xxyz_1,   \
                             ta1_z_zzz_xxyzz_0,  \
                             ta1_z_zzz_xxyzz_1,  \
                             ta1_z_zzz_xxzz_0,   \
                             ta1_z_zzz_xxzz_1,   \
                             ta1_z_zzz_xxzzz_0,  \
                             ta1_z_zzz_xxzzz_1,  \
                             ta1_z_zzz_xyyy_0,   \
                             ta1_z_zzz_xyyy_1,   \
                             ta1_z_zzz_xyyyy_0,  \
                             ta1_z_zzz_xyyyy_1,  \
                             ta1_z_zzz_xyyyz_0,  \
                             ta1_z_zzz_xyyyz_1,  \
                             ta1_z_zzz_xyyz_0,   \
                             ta1_z_zzz_xyyz_1,   \
                             ta1_z_zzz_xyyzz_0,  \
                             ta1_z_zzz_xyyzz_1,  \
                             ta1_z_zzz_xyzz_0,   \
                             ta1_z_zzz_xyzz_1,   \
                             ta1_z_zzz_xyzzz_0,  \
                             ta1_z_zzz_xyzzz_1,  \
                             ta1_z_zzz_xzzz_0,   \
                             ta1_z_zzz_xzzz_1,   \
                             ta1_z_zzz_xzzzz_0,  \
                             ta1_z_zzz_xzzzz_1,  \
                             ta1_z_zzz_yyyy_0,   \
                             ta1_z_zzz_yyyy_1,   \
                             ta1_z_zzz_yyyyy_0,  \
                             ta1_z_zzz_yyyyy_1,  \
                             ta1_z_zzz_yyyyz_0,  \
                             ta1_z_zzz_yyyyz_1,  \
                             ta1_z_zzz_yyyz_0,   \
                             ta1_z_zzz_yyyz_1,   \
                             ta1_z_zzz_yyyzz_0,  \
                             ta1_z_zzz_yyyzz_1,  \
                             ta1_z_zzz_yyzz_0,   \
                             ta1_z_zzz_yyzz_1,   \
                             ta1_z_zzz_yyzzz_0,  \
                             ta1_z_zzz_yyzzz_1,  \
                             ta1_z_zzz_yzzz_0,   \
                             ta1_z_zzz_yzzz_1,   \
                             ta1_z_zzz_yzzzz_0,  \
                             ta1_z_zzz_yzzzz_1,  \
                             ta1_z_zzz_zzzz_0,   \
                             ta1_z_zzz_zzzz_1,   \
                             ta1_z_zzz_zzzzz_0,  \
                             ta1_z_zzz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_xzzz_xxxxx_0[i] =
            5.0 * ta1_z_zzz_xxxx_0[i] * fe_0 - 5.0 * ta1_z_zzz_xxxx_1[i] * fe_0 + ta1_z_zzz_xxxxx_0[i] * pa_x[i] - ta1_z_zzz_xxxxx_1[i] * pc_x[i];

        ta1_z_xzzz_xxxxy_0[i] =
            4.0 * ta1_z_zzz_xxxy_0[i] * fe_0 - 4.0 * ta1_z_zzz_xxxy_1[i] * fe_0 + ta1_z_zzz_xxxxy_0[i] * pa_x[i] - ta1_z_zzz_xxxxy_1[i] * pc_x[i];

        ta1_z_xzzz_xxxxz_0[i] =
            4.0 * ta1_z_zzz_xxxz_0[i] * fe_0 - 4.0 * ta1_z_zzz_xxxz_1[i] * fe_0 + ta1_z_zzz_xxxxz_0[i] * pa_x[i] - ta1_z_zzz_xxxxz_1[i] * pc_x[i];

        ta1_z_xzzz_xxxyy_0[i] =
            3.0 * ta1_z_zzz_xxyy_0[i] * fe_0 - 3.0 * ta1_z_zzz_xxyy_1[i] * fe_0 + ta1_z_zzz_xxxyy_0[i] * pa_x[i] - ta1_z_zzz_xxxyy_1[i] * pc_x[i];

        ta1_z_xzzz_xxxyz_0[i] =
            3.0 * ta1_z_zzz_xxyz_0[i] * fe_0 - 3.0 * ta1_z_zzz_xxyz_1[i] * fe_0 + ta1_z_zzz_xxxyz_0[i] * pa_x[i] - ta1_z_zzz_xxxyz_1[i] * pc_x[i];

        ta1_z_xzzz_xxxzz_0[i] =
            3.0 * ta1_z_zzz_xxzz_0[i] * fe_0 - 3.0 * ta1_z_zzz_xxzz_1[i] * fe_0 + ta1_z_zzz_xxxzz_0[i] * pa_x[i] - ta1_z_zzz_xxxzz_1[i] * pc_x[i];

        ta1_z_xzzz_xxyyy_0[i] =
            2.0 * ta1_z_zzz_xyyy_0[i] * fe_0 - 2.0 * ta1_z_zzz_xyyy_1[i] * fe_0 + ta1_z_zzz_xxyyy_0[i] * pa_x[i] - ta1_z_zzz_xxyyy_1[i] * pc_x[i];

        ta1_z_xzzz_xxyyz_0[i] =
            2.0 * ta1_z_zzz_xyyz_0[i] * fe_0 - 2.0 * ta1_z_zzz_xyyz_1[i] * fe_0 + ta1_z_zzz_xxyyz_0[i] * pa_x[i] - ta1_z_zzz_xxyyz_1[i] * pc_x[i];

        ta1_z_xzzz_xxyzz_0[i] =
            2.0 * ta1_z_zzz_xyzz_0[i] * fe_0 - 2.0 * ta1_z_zzz_xyzz_1[i] * fe_0 + ta1_z_zzz_xxyzz_0[i] * pa_x[i] - ta1_z_zzz_xxyzz_1[i] * pc_x[i];

        ta1_z_xzzz_xxzzz_0[i] =
            2.0 * ta1_z_zzz_xzzz_0[i] * fe_0 - 2.0 * ta1_z_zzz_xzzz_1[i] * fe_0 + ta1_z_zzz_xxzzz_0[i] * pa_x[i] - ta1_z_zzz_xxzzz_1[i] * pc_x[i];

        ta1_z_xzzz_xyyyy_0[i] =
            ta1_z_zzz_yyyy_0[i] * fe_0 - ta1_z_zzz_yyyy_1[i] * fe_0 + ta1_z_zzz_xyyyy_0[i] * pa_x[i] - ta1_z_zzz_xyyyy_1[i] * pc_x[i];

        ta1_z_xzzz_xyyyz_0[i] =
            ta1_z_zzz_yyyz_0[i] * fe_0 - ta1_z_zzz_yyyz_1[i] * fe_0 + ta1_z_zzz_xyyyz_0[i] * pa_x[i] - ta1_z_zzz_xyyyz_1[i] * pc_x[i];

        ta1_z_xzzz_xyyzz_0[i] =
            ta1_z_zzz_yyzz_0[i] * fe_0 - ta1_z_zzz_yyzz_1[i] * fe_0 + ta1_z_zzz_xyyzz_0[i] * pa_x[i] - ta1_z_zzz_xyyzz_1[i] * pc_x[i];

        ta1_z_xzzz_xyzzz_0[i] =
            ta1_z_zzz_yzzz_0[i] * fe_0 - ta1_z_zzz_yzzz_1[i] * fe_0 + ta1_z_zzz_xyzzz_0[i] * pa_x[i] - ta1_z_zzz_xyzzz_1[i] * pc_x[i];

        ta1_z_xzzz_xzzzz_0[i] =
            ta1_z_zzz_zzzz_0[i] * fe_0 - ta1_z_zzz_zzzz_1[i] * fe_0 + ta1_z_zzz_xzzzz_0[i] * pa_x[i] - ta1_z_zzz_xzzzz_1[i] * pc_x[i];

        ta1_z_xzzz_yyyyy_0[i] = ta1_z_zzz_yyyyy_0[i] * pa_x[i] - ta1_z_zzz_yyyyy_1[i] * pc_x[i];

        ta1_z_xzzz_yyyyz_0[i] = ta1_z_zzz_yyyyz_0[i] * pa_x[i] - ta1_z_zzz_yyyyz_1[i] * pc_x[i];

        ta1_z_xzzz_yyyzz_0[i] = ta1_z_zzz_yyyzz_0[i] * pa_x[i] - ta1_z_zzz_yyyzz_1[i] * pc_x[i];

        ta1_z_xzzz_yyzzz_0[i] = ta1_z_zzz_yyzzz_0[i] * pa_x[i] - ta1_z_zzz_yyzzz_1[i] * pc_x[i];

        ta1_z_xzzz_yzzzz_0[i] = ta1_z_zzz_yzzzz_0[i] * pa_x[i] - ta1_z_zzz_yzzzz_1[i] * pc_x[i];

        ta1_z_xzzz_zzzzz_0[i] = ta1_z_zzz_zzzzz_0[i] * pa_x[i] - ta1_z_zzz_zzzzz_1[i] * pc_x[i];
    }

    // Set up 840-861 components of targeted buffer : GH

    auto ta1_z_yyyy_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 840);

    auto ta1_z_yyyy_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 841);

    auto ta1_z_yyyy_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 842);

    auto ta1_z_yyyy_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 843);

    auto ta1_z_yyyy_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 844);

    auto ta1_z_yyyy_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 845);

    auto ta1_z_yyyy_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 846);

    auto ta1_z_yyyy_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 847);

    auto ta1_z_yyyy_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 848);

    auto ta1_z_yyyy_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 849);

    auto ta1_z_yyyy_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 850);

    auto ta1_z_yyyy_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 851);

    auto ta1_z_yyyy_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 852);

    auto ta1_z_yyyy_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 853);

    auto ta1_z_yyyy_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 854);

    auto ta1_z_yyyy_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 855);

    auto ta1_z_yyyy_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 856);

    auto ta1_z_yyyy_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 857);

    auto ta1_z_yyyy_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 858);

    auto ta1_z_yyyy_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 859);

    auto ta1_z_yyyy_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 860);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_z_yy_xxxxx_0,   \
                             ta1_z_yy_xxxxx_1,   \
                             ta1_z_yy_xxxxy_0,   \
                             ta1_z_yy_xxxxy_1,   \
                             ta1_z_yy_xxxxz_0,   \
                             ta1_z_yy_xxxxz_1,   \
                             ta1_z_yy_xxxyy_0,   \
                             ta1_z_yy_xxxyy_1,   \
                             ta1_z_yy_xxxyz_0,   \
                             ta1_z_yy_xxxyz_1,   \
                             ta1_z_yy_xxxzz_0,   \
                             ta1_z_yy_xxxzz_1,   \
                             ta1_z_yy_xxyyy_0,   \
                             ta1_z_yy_xxyyy_1,   \
                             ta1_z_yy_xxyyz_0,   \
                             ta1_z_yy_xxyyz_1,   \
                             ta1_z_yy_xxyzz_0,   \
                             ta1_z_yy_xxyzz_1,   \
                             ta1_z_yy_xxzzz_0,   \
                             ta1_z_yy_xxzzz_1,   \
                             ta1_z_yy_xyyyy_0,   \
                             ta1_z_yy_xyyyy_1,   \
                             ta1_z_yy_xyyyz_0,   \
                             ta1_z_yy_xyyyz_1,   \
                             ta1_z_yy_xyyzz_0,   \
                             ta1_z_yy_xyyzz_1,   \
                             ta1_z_yy_xyzzz_0,   \
                             ta1_z_yy_xyzzz_1,   \
                             ta1_z_yy_xzzzz_0,   \
                             ta1_z_yy_xzzzz_1,   \
                             ta1_z_yy_yyyyy_0,   \
                             ta1_z_yy_yyyyy_1,   \
                             ta1_z_yy_yyyyz_0,   \
                             ta1_z_yy_yyyyz_1,   \
                             ta1_z_yy_yyyzz_0,   \
                             ta1_z_yy_yyyzz_1,   \
                             ta1_z_yy_yyzzz_0,   \
                             ta1_z_yy_yyzzz_1,   \
                             ta1_z_yy_yzzzz_0,   \
                             ta1_z_yy_yzzzz_1,   \
                             ta1_z_yy_zzzzz_0,   \
                             ta1_z_yy_zzzzz_1,   \
                             ta1_z_yyy_xxxx_0,   \
                             ta1_z_yyy_xxxx_1,   \
                             ta1_z_yyy_xxxxx_0,  \
                             ta1_z_yyy_xxxxx_1,  \
                             ta1_z_yyy_xxxxy_0,  \
                             ta1_z_yyy_xxxxy_1,  \
                             ta1_z_yyy_xxxxz_0,  \
                             ta1_z_yyy_xxxxz_1,  \
                             ta1_z_yyy_xxxy_0,   \
                             ta1_z_yyy_xxxy_1,   \
                             ta1_z_yyy_xxxyy_0,  \
                             ta1_z_yyy_xxxyy_1,  \
                             ta1_z_yyy_xxxyz_0,  \
                             ta1_z_yyy_xxxyz_1,  \
                             ta1_z_yyy_xxxz_0,   \
                             ta1_z_yyy_xxxz_1,   \
                             ta1_z_yyy_xxxzz_0,  \
                             ta1_z_yyy_xxxzz_1,  \
                             ta1_z_yyy_xxyy_0,   \
                             ta1_z_yyy_xxyy_1,   \
                             ta1_z_yyy_xxyyy_0,  \
                             ta1_z_yyy_xxyyy_1,  \
                             ta1_z_yyy_xxyyz_0,  \
                             ta1_z_yyy_xxyyz_1,  \
                             ta1_z_yyy_xxyz_0,   \
                             ta1_z_yyy_xxyz_1,   \
                             ta1_z_yyy_xxyzz_0,  \
                             ta1_z_yyy_xxyzz_1,  \
                             ta1_z_yyy_xxzz_0,   \
                             ta1_z_yyy_xxzz_1,   \
                             ta1_z_yyy_xxzzz_0,  \
                             ta1_z_yyy_xxzzz_1,  \
                             ta1_z_yyy_xyyy_0,   \
                             ta1_z_yyy_xyyy_1,   \
                             ta1_z_yyy_xyyyy_0,  \
                             ta1_z_yyy_xyyyy_1,  \
                             ta1_z_yyy_xyyyz_0,  \
                             ta1_z_yyy_xyyyz_1,  \
                             ta1_z_yyy_xyyz_0,   \
                             ta1_z_yyy_xyyz_1,   \
                             ta1_z_yyy_xyyzz_0,  \
                             ta1_z_yyy_xyyzz_1,  \
                             ta1_z_yyy_xyzz_0,   \
                             ta1_z_yyy_xyzz_1,   \
                             ta1_z_yyy_xyzzz_0,  \
                             ta1_z_yyy_xyzzz_1,  \
                             ta1_z_yyy_xzzz_0,   \
                             ta1_z_yyy_xzzz_1,   \
                             ta1_z_yyy_xzzzz_0,  \
                             ta1_z_yyy_xzzzz_1,  \
                             ta1_z_yyy_yyyy_0,   \
                             ta1_z_yyy_yyyy_1,   \
                             ta1_z_yyy_yyyyy_0,  \
                             ta1_z_yyy_yyyyy_1,  \
                             ta1_z_yyy_yyyyz_0,  \
                             ta1_z_yyy_yyyyz_1,  \
                             ta1_z_yyy_yyyz_0,   \
                             ta1_z_yyy_yyyz_1,   \
                             ta1_z_yyy_yyyzz_0,  \
                             ta1_z_yyy_yyyzz_1,  \
                             ta1_z_yyy_yyzz_0,   \
                             ta1_z_yyy_yyzz_1,   \
                             ta1_z_yyy_yyzzz_0,  \
                             ta1_z_yyy_yyzzz_1,  \
                             ta1_z_yyy_yzzz_0,   \
                             ta1_z_yyy_yzzz_1,   \
                             ta1_z_yyy_yzzzz_0,  \
                             ta1_z_yyy_yzzzz_1,  \
                             ta1_z_yyy_zzzz_0,   \
                             ta1_z_yyy_zzzz_1,   \
                             ta1_z_yyy_zzzzz_0,  \
                             ta1_z_yyy_zzzzz_1,  \
                             ta1_z_yyyy_xxxxx_0, \
                             ta1_z_yyyy_xxxxy_0, \
                             ta1_z_yyyy_xxxxz_0, \
                             ta1_z_yyyy_xxxyy_0, \
                             ta1_z_yyyy_xxxyz_0, \
                             ta1_z_yyyy_xxxzz_0, \
                             ta1_z_yyyy_xxyyy_0, \
                             ta1_z_yyyy_xxyyz_0, \
                             ta1_z_yyyy_xxyzz_0, \
                             ta1_z_yyyy_xxzzz_0, \
                             ta1_z_yyyy_xyyyy_0, \
                             ta1_z_yyyy_xyyyz_0, \
                             ta1_z_yyyy_xyyzz_0, \
                             ta1_z_yyyy_xyzzz_0, \
                             ta1_z_yyyy_xzzzz_0, \
                             ta1_z_yyyy_yyyyy_0, \
                             ta1_z_yyyy_yyyyz_0, \
                             ta1_z_yyyy_yyyzz_0, \
                             ta1_z_yyyy_yyzzz_0, \
                             ta1_z_yyyy_yzzzz_0, \
                             ta1_z_yyyy_zzzzz_0, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyy_xxxxx_0[i] =
            3.0 * ta1_z_yy_xxxxx_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxxx_1[i] * fe_0 + ta1_z_yyy_xxxxx_0[i] * pa_y[i] - ta1_z_yyy_xxxxx_1[i] * pc_y[i];

        ta1_z_yyyy_xxxxy_0[i] = 3.0 * ta1_z_yy_xxxxy_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxxy_1[i] * fe_0 + ta1_z_yyy_xxxx_0[i] * fe_0 -
                                ta1_z_yyy_xxxx_1[i] * fe_0 + ta1_z_yyy_xxxxy_0[i] * pa_y[i] - ta1_z_yyy_xxxxy_1[i] * pc_y[i];

        ta1_z_yyyy_xxxxz_0[i] =
            3.0 * ta1_z_yy_xxxxz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxxz_1[i] * fe_0 + ta1_z_yyy_xxxxz_0[i] * pa_y[i] - ta1_z_yyy_xxxxz_1[i] * pc_y[i];

        ta1_z_yyyy_xxxyy_0[i] = 3.0 * ta1_z_yy_xxxyy_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxyy_1[i] * fe_0 + 2.0 * ta1_z_yyy_xxxy_0[i] * fe_0 -
                                2.0 * ta1_z_yyy_xxxy_1[i] * fe_0 + ta1_z_yyy_xxxyy_0[i] * pa_y[i] - ta1_z_yyy_xxxyy_1[i] * pc_y[i];

        ta1_z_yyyy_xxxyz_0[i] = 3.0 * ta1_z_yy_xxxyz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxyz_1[i] * fe_0 + ta1_z_yyy_xxxz_0[i] * fe_0 -
                                ta1_z_yyy_xxxz_1[i] * fe_0 + ta1_z_yyy_xxxyz_0[i] * pa_y[i] - ta1_z_yyy_xxxyz_1[i] * pc_y[i];

        ta1_z_yyyy_xxxzz_0[i] =
            3.0 * ta1_z_yy_xxxzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxxzz_1[i] * fe_0 + ta1_z_yyy_xxxzz_0[i] * pa_y[i] - ta1_z_yyy_xxxzz_1[i] * pc_y[i];

        ta1_z_yyyy_xxyyy_0[i] = 3.0 * ta1_z_yy_xxyyy_0[i] * fe_0 - 3.0 * ta1_z_yy_xxyyy_1[i] * fe_0 + 3.0 * ta1_z_yyy_xxyy_0[i] * fe_0 -
                                3.0 * ta1_z_yyy_xxyy_1[i] * fe_0 + ta1_z_yyy_xxyyy_0[i] * pa_y[i] - ta1_z_yyy_xxyyy_1[i] * pc_y[i];

        ta1_z_yyyy_xxyyz_0[i] = 3.0 * ta1_z_yy_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxyyz_1[i] * fe_0 + 2.0 * ta1_z_yyy_xxyz_0[i] * fe_0 -
                                2.0 * ta1_z_yyy_xxyz_1[i] * fe_0 + ta1_z_yyy_xxyyz_0[i] * pa_y[i] - ta1_z_yyy_xxyyz_1[i] * pc_y[i];

        ta1_z_yyyy_xxyzz_0[i] = 3.0 * ta1_z_yy_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxyzz_1[i] * fe_0 + ta1_z_yyy_xxzz_0[i] * fe_0 -
                                ta1_z_yyy_xxzz_1[i] * fe_0 + ta1_z_yyy_xxyzz_0[i] * pa_y[i] - ta1_z_yyy_xxyzz_1[i] * pc_y[i];

        ta1_z_yyyy_xxzzz_0[i] =
            3.0 * ta1_z_yy_xxzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xxzzz_1[i] * fe_0 + ta1_z_yyy_xxzzz_0[i] * pa_y[i] - ta1_z_yyy_xxzzz_1[i] * pc_y[i];

        ta1_z_yyyy_xyyyy_0[i] = 3.0 * ta1_z_yy_xyyyy_0[i] * fe_0 - 3.0 * ta1_z_yy_xyyyy_1[i] * fe_0 + 4.0 * ta1_z_yyy_xyyy_0[i] * fe_0 -
                                4.0 * ta1_z_yyy_xyyy_1[i] * fe_0 + ta1_z_yyy_xyyyy_0[i] * pa_y[i] - ta1_z_yyy_xyyyy_1[i] * pc_y[i];

        ta1_z_yyyy_xyyyz_0[i] = 3.0 * ta1_z_yy_xyyyz_0[i] * fe_0 - 3.0 * ta1_z_yy_xyyyz_1[i] * fe_0 + 3.0 * ta1_z_yyy_xyyz_0[i] * fe_0 -
                                3.0 * ta1_z_yyy_xyyz_1[i] * fe_0 + ta1_z_yyy_xyyyz_0[i] * pa_y[i] - ta1_z_yyy_xyyyz_1[i] * pc_y[i];

        ta1_z_yyyy_xyyzz_0[i] = 3.0 * ta1_z_yy_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xyyzz_1[i] * fe_0 + 2.0 * ta1_z_yyy_xyzz_0[i] * fe_0 -
                                2.0 * ta1_z_yyy_xyzz_1[i] * fe_0 + ta1_z_yyy_xyyzz_0[i] * pa_y[i] - ta1_z_yyy_xyyzz_1[i] * pc_y[i];

        ta1_z_yyyy_xyzzz_0[i] = 3.0 * ta1_z_yy_xyzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xyzzz_1[i] * fe_0 + ta1_z_yyy_xzzz_0[i] * fe_0 -
                                ta1_z_yyy_xzzz_1[i] * fe_0 + ta1_z_yyy_xyzzz_0[i] * pa_y[i] - ta1_z_yyy_xyzzz_1[i] * pc_y[i];

        ta1_z_yyyy_xzzzz_0[i] =
            3.0 * ta1_z_yy_xzzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_xzzzz_1[i] * fe_0 + ta1_z_yyy_xzzzz_0[i] * pa_y[i] - ta1_z_yyy_xzzzz_1[i] * pc_y[i];

        ta1_z_yyyy_yyyyy_0[i] = 3.0 * ta1_z_yy_yyyyy_0[i] * fe_0 - 3.0 * ta1_z_yy_yyyyy_1[i] * fe_0 + 5.0 * ta1_z_yyy_yyyy_0[i] * fe_0 -
                                5.0 * ta1_z_yyy_yyyy_1[i] * fe_0 + ta1_z_yyy_yyyyy_0[i] * pa_y[i] - ta1_z_yyy_yyyyy_1[i] * pc_y[i];

        ta1_z_yyyy_yyyyz_0[i] = 3.0 * ta1_z_yy_yyyyz_0[i] * fe_0 - 3.0 * ta1_z_yy_yyyyz_1[i] * fe_0 + 4.0 * ta1_z_yyy_yyyz_0[i] * fe_0 -
                                4.0 * ta1_z_yyy_yyyz_1[i] * fe_0 + ta1_z_yyy_yyyyz_0[i] * pa_y[i] - ta1_z_yyy_yyyyz_1[i] * pc_y[i];

        ta1_z_yyyy_yyyzz_0[i] = 3.0 * ta1_z_yy_yyyzz_0[i] * fe_0 - 3.0 * ta1_z_yy_yyyzz_1[i] * fe_0 + 3.0 * ta1_z_yyy_yyzz_0[i] * fe_0 -
                                3.0 * ta1_z_yyy_yyzz_1[i] * fe_0 + ta1_z_yyy_yyyzz_0[i] * pa_y[i] - ta1_z_yyy_yyyzz_1[i] * pc_y[i];

        ta1_z_yyyy_yyzzz_0[i] = 3.0 * ta1_z_yy_yyzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_yyzzz_1[i] * fe_0 + 2.0 * ta1_z_yyy_yzzz_0[i] * fe_0 -
                                2.0 * ta1_z_yyy_yzzz_1[i] * fe_0 + ta1_z_yyy_yyzzz_0[i] * pa_y[i] - ta1_z_yyy_yyzzz_1[i] * pc_y[i];

        ta1_z_yyyy_yzzzz_0[i] = 3.0 * ta1_z_yy_yzzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_yzzzz_1[i] * fe_0 + ta1_z_yyy_zzzz_0[i] * fe_0 -
                                ta1_z_yyy_zzzz_1[i] * fe_0 + ta1_z_yyy_yzzzz_0[i] * pa_y[i] - ta1_z_yyy_yzzzz_1[i] * pc_y[i];

        ta1_z_yyyy_zzzzz_0[i] =
            3.0 * ta1_z_yy_zzzzz_0[i] * fe_0 - 3.0 * ta1_z_yy_zzzzz_1[i] * fe_0 + ta1_z_yyy_zzzzz_0[i] * pa_y[i] - ta1_z_yyy_zzzzz_1[i] * pc_y[i];
    }

    // Set up 861-882 components of targeted buffer : GH

    auto ta1_z_yyyz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 861);

    auto ta1_z_yyyz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 862);

    auto ta1_z_yyyz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 863);

    auto ta1_z_yyyz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 864);

    auto ta1_z_yyyz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 865);

    auto ta1_z_yyyz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 866);

    auto ta1_z_yyyz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 867);

    auto ta1_z_yyyz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 868);

    auto ta1_z_yyyz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 869);

    auto ta1_z_yyyz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 870);

    auto ta1_z_yyyz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 871);

    auto ta1_z_yyyz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 872);

    auto ta1_z_yyyz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 873);

    auto ta1_z_yyyz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 874);

    auto ta1_z_yyyz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 875);

    auto ta1_z_yyyz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 876);

    auto ta1_z_yyyz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 877);

    auto ta1_z_yyyz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 878);

    auto ta1_z_yyyz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 879);

    auto ta1_z_yyyz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 880);

    auto ta1_z_yyyz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 881);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_yyy_xxxxx_0,  \
                             ta1_z_yyy_xxxxx_1,  \
                             ta1_z_yyy_xxxxy_0,  \
                             ta1_z_yyy_xxxxy_1,  \
                             ta1_z_yyy_xxxy_0,   \
                             ta1_z_yyy_xxxy_1,   \
                             ta1_z_yyy_xxxyy_0,  \
                             ta1_z_yyy_xxxyy_1,  \
                             ta1_z_yyy_xxxyz_0,  \
                             ta1_z_yyy_xxxyz_1,  \
                             ta1_z_yyy_xxyy_0,   \
                             ta1_z_yyy_xxyy_1,   \
                             ta1_z_yyy_xxyyy_0,  \
                             ta1_z_yyy_xxyyy_1,  \
                             ta1_z_yyy_xxyyz_0,  \
                             ta1_z_yyy_xxyyz_1,  \
                             ta1_z_yyy_xxyz_0,   \
                             ta1_z_yyy_xxyz_1,   \
                             ta1_z_yyy_xxyzz_0,  \
                             ta1_z_yyy_xxyzz_1,  \
                             ta1_z_yyy_xyyy_0,   \
                             ta1_z_yyy_xyyy_1,   \
                             ta1_z_yyy_xyyyy_0,  \
                             ta1_z_yyy_xyyyy_1,  \
                             ta1_z_yyy_xyyyz_0,  \
                             ta1_z_yyy_xyyyz_1,  \
                             ta1_z_yyy_xyyz_0,   \
                             ta1_z_yyy_xyyz_1,   \
                             ta1_z_yyy_xyyzz_0,  \
                             ta1_z_yyy_xyyzz_1,  \
                             ta1_z_yyy_xyzz_0,   \
                             ta1_z_yyy_xyzz_1,   \
                             ta1_z_yyy_xyzzz_0,  \
                             ta1_z_yyy_xyzzz_1,  \
                             ta1_z_yyy_yyyy_0,   \
                             ta1_z_yyy_yyyy_1,   \
                             ta1_z_yyy_yyyyy_0,  \
                             ta1_z_yyy_yyyyy_1,  \
                             ta1_z_yyy_yyyyz_0,  \
                             ta1_z_yyy_yyyyz_1,  \
                             ta1_z_yyy_yyyz_0,   \
                             ta1_z_yyy_yyyz_1,   \
                             ta1_z_yyy_yyyzz_0,  \
                             ta1_z_yyy_yyyzz_1,  \
                             ta1_z_yyy_yyzz_0,   \
                             ta1_z_yyy_yyzz_1,   \
                             ta1_z_yyy_yyzzz_0,  \
                             ta1_z_yyy_yyzzz_1,  \
                             ta1_z_yyy_yzzz_0,   \
                             ta1_z_yyy_yzzz_1,   \
                             ta1_z_yyy_yzzzz_0,  \
                             ta1_z_yyy_yzzzz_1,  \
                             ta1_z_yyyz_xxxxx_0, \
                             ta1_z_yyyz_xxxxy_0, \
                             ta1_z_yyyz_xxxxz_0, \
                             ta1_z_yyyz_xxxyy_0, \
                             ta1_z_yyyz_xxxyz_0, \
                             ta1_z_yyyz_xxxzz_0, \
                             ta1_z_yyyz_xxyyy_0, \
                             ta1_z_yyyz_xxyyz_0, \
                             ta1_z_yyyz_xxyzz_0, \
                             ta1_z_yyyz_xxzzz_0, \
                             ta1_z_yyyz_xyyyy_0, \
                             ta1_z_yyyz_xyyyz_0, \
                             ta1_z_yyyz_xyyzz_0, \
                             ta1_z_yyyz_xyzzz_0, \
                             ta1_z_yyyz_xzzzz_0, \
                             ta1_z_yyyz_yyyyy_0, \
                             ta1_z_yyyz_yyyyz_0, \
                             ta1_z_yyyz_yyyzz_0, \
                             ta1_z_yyyz_yyzzz_0, \
                             ta1_z_yyyz_yzzzz_0, \
                             ta1_z_yyyz_zzzzz_0, \
                             ta1_z_yyz_xxxxz_0,  \
                             ta1_z_yyz_xxxxz_1,  \
                             ta1_z_yyz_xxxzz_0,  \
                             ta1_z_yyz_xxxzz_1,  \
                             ta1_z_yyz_xxzzz_0,  \
                             ta1_z_yyz_xxzzz_1,  \
                             ta1_z_yyz_xzzzz_0,  \
                             ta1_z_yyz_xzzzz_1,  \
                             ta1_z_yyz_zzzzz_0,  \
                             ta1_z_yyz_zzzzz_1,  \
                             ta1_z_yz_xxxxz_0,   \
                             ta1_z_yz_xxxxz_1,   \
                             ta1_z_yz_xxxzz_0,   \
                             ta1_z_yz_xxxzz_1,   \
                             ta1_z_yz_xxzzz_0,   \
                             ta1_z_yz_xxzzz_1,   \
                             ta1_z_yz_xzzzz_0,   \
                             ta1_z_yz_xzzzz_1,   \
                             ta1_z_yz_zzzzz_0,   \
                             ta1_z_yz_zzzzz_1,   \
                             ta_yyy_xxxxx_1,     \
                             ta_yyy_xxxxy_1,     \
                             ta_yyy_xxxyy_1,     \
                             ta_yyy_xxxyz_1,     \
                             ta_yyy_xxyyy_1,     \
                             ta_yyy_xxyyz_1,     \
                             ta_yyy_xxyzz_1,     \
                             ta_yyy_xyyyy_1,     \
                             ta_yyy_xyyyz_1,     \
                             ta_yyy_xyyzz_1,     \
                             ta_yyy_xyzzz_1,     \
                             ta_yyy_yyyyy_1,     \
                             ta_yyy_yyyyz_1,     \
                             ta_yyy_yyyzz_1,     \
                             ta_yyy_yyzzz_1,     \
                             ta_yyy_yzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyyz_xxxxx_0[i] = ta_yyy_xxxxx_1[i] + ta1_z_yyy_xxxxx_0[i] * pa_z[i] - ta1_z_yyy_xxxxx_1[i] * pc_z[i];

        ta1_z_yyyz_xxxxy_0[i] = ta_yyy_xxxxy_1[i] + ta1_z_yyy_xxxxy_0[i] * pa_z[i] - ta1_z_yyy_xxxxy_1[i] * pc_z[i];

        ta1_z_yyyz_xxxxz_0[i] =
            2.0 * ta1_z_yz_xxxxz_0[i] * fe_0 - 2.0 * ta1_z_yz_xxxxz_1[i] * fe_0 + ta1_z_yyz_xxxxz_0[i] * pa_y[i] - ta1_z_yyz_xxxxz_1[i] * pc_y[i];

        ta1_z_yyyz_xxxyy_0[i] = ta_yyy_xxxyy_1[i] + ta1_z_yyy_xxxyy_0[i] * pa_z[i] - ta1_z_yyy_xxxyy_1[i] * pc_z[i];

        ta1_z_yyyz_xxxyz_0[i] = ta1_z_yyy_xxxy_0[i] * fe_0 - ta1_z_yyy_xxxy_1[i] * fe_0 + ta_yyy_xxxyz_1[i] + ta1_z_yyy_xxxyz_0[i] * pa_z[i] -
                                ta1_z_yyy_xxxyz_1[i] * pc_z[i];

        ta1_z_yyyz_xxxzz_0[i] =
            2.0 * ta1_z_yz_xxxzz_0[i] * fe_0 - 2.0 * ta1_z_yz_xxxzz_1[i] * fe_0 + ta1_z_yyz_xxxzz_0[i] * pa_y[i] - ta1_z_yyz_xxxzz_1[i] * pc_y[i];

        ta1_z_yyyz_xxyyy_0[i] = ta_yyy_xxyyy_1[i] + ta1_z_yyy_xxyyy_0[i] * pa_z[i] - ta1_z_yyy_xxyyy_1[i] * pc_z[i];

        ta1_z_yyyz_xxyyz_0[i] = ta1_z_yyy_xxyy_0[i] * fe_0 - ta1_z_yyy_xxyy_1[i] * fe_0 + ta_yyy_xxyyz_1[i] + ta1_z_yyy_xxyyz_0[i] * pa_z[i] -
                                ta1_z_yyy_xxyyz_1[i] * pc_z[i];

        ta1_z_yyyz_xxyzz_0[i] = 2.0 * ta1_z_yyy_xxyz_0[i] * fe_0 - 2.0 * ta1_z_yyy_xxyz_1[i] * fe_0 + ta_yyy_xxyzz_1[i] +
                                ta1_z_yyy_xxyzz_0[i] * pa_z[i] - ta1_z_yyy_xxyzz_1[i] * pc_z[i];

        ta1_z_yyyz_xxzzz_0[i] =
            2.0 * ta1_z_yz_xxzzz_0[i] * fe_0 - 2.0 * ta1_z_yz_xxzzz_1[i] * fe_0 + ta1_z_yyz_xxzzz_0[i] * pa_y[i] - ta1_z_yyz_xxzzz_1[i] * pc_y[i];

        ta1_z_yyyz_xyyyy_0[i] = ta_yyy_xyyyy_1[i] + ta1_z_yyy_xyyyy_0[i] * pa_z[i] - ta1_z_yyy_xyyyy_1[i] * pc_z[i];

        ta1_z_yyyz_xyyyz_0[i] = ta1_z_yyy_xyyy_0[i] * fe_0 - ta1_z_yyy_xyyy_1[i] * fe_0 + ta_yyy_xyyyz_1[i] + ta1_z_yyy_xyyyz_0[i] * pa_z[i] -
                                ta1_z_yyy_xyyyz_1[i] * pc_z[i];

        ta1_z_yyyz_xyyzz_0[i] = 2.0 * ta1_z_yyy_xyyz_0[i] * fe_0 - 2.0 * ta1_z_yyy_xyyz_1[i] * fe_0 + ta_yyy_xyyzz_1[i] +
                                ta1_z_yyy_xyyzz_0[i] * pa_z[i] - ta1_z_yyy_xyyzz_1[i] * pc_z[i];

        ta1_z_yyyz_xyzzz_0[i] = 3.0 * ta1_z_yyy_xyzz_0[i] * fe_0 - 3.0 * ta1_z_yyy_xyzz_1[i] * fe_0 + ta_yyy_xyzzz_1[i] +
                                ta1_z_yyy_xyzzz_0[i] * pa_z[i] - ta1_z_yyy_xyzzz_1[i] * pc_z[i];

        ta1_z_yyyz_xzzzz_0[i] =
            2.0 * ta1_z_yz_xzzzz_0[i] * fe_0 - 2.0 * ta1_z_yz_xzzzz_1[i] * fe_0 + ta1_z_yyz_xzzzz_0[i] * pa_y[i] - ta1_z_yyz_xzzzz_1[i] * pc_y[i];

        ta1_z_yyyz_yyyyy_0[i] = ta_yyy_yyyyy_1[i] + ta1_z_yyy_yyyyy_0[i] * pa_z[i] - ta1_z_yyy_yyyyy_1[i] * pc_z[i];

        ta1_z_yyyz_yyyyz_0[i] = ta1_z_yyy_yyyy_0[i] * fe_0 - ta1_z_yyy_yyyy_1[i] * fe_0 + ta_yyy_yyyyz_1[i] + ta1_z_yyy_yyyyz_0[i] * pa_z[i] -
                                ta1_z_yyy_yyyyz_1[i] * pc_z[i];

        ta1_z_yyyz_yyyzz_0[i] = 2.0 * ta1_z_yyy_yyyz_0[i] * fe_0 - 2.0 * ta1_z_yyy_yyyz_1[i] * fe_0 + ta_yyy_yyyzz_1[i] +
                                ta1_z_yyy_yyyzz_0[i] * pa_z[i] - ta1_z_yyy_yyyzz_1[i] * pc_z[i];

        ta1_z_yyyz_yyzzz_0[i] = 3.0 * ta1_z_yyy_yyzz_0[i] * fe_0 - 3.0 * ta1_z_yyy_yyzz_1[i] * fe_0 + ta_yyy_yyzzz_1[i] +
                                ta1_z_yyy_yyzzz_0[i] * pa_z[i] - ta1_z_yyy_yyzzz_1[i] * pc_z[i];

        ta1_z_yyyz_yzzzz_0[i] = 4.0 * ta1_z_yyy_yzzz_0[i] * fe_0 - 4.0 * ta1_z_yyy_yzzz_1[i] * fe_0 + ta_yyy_yzzzz_1[i] +
                                ta1_z_yyy_yzzzz_0[i] * pa_z[i] - ta1_z_yyy_yzzzz_1[i] * pc_z[i];

        ta1_z_yyyz_zzzzz_0[i] =
            2.0 * ta1_z_yz_zzzzz_0[i] * fe_0 - 2.0 * ta1_z_yz_zzzzz_1[i] * fe_0 + ta1_z_yyz_zzzzz_0[i] * pa_y[i] - ta1_z_yyz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 882-903 components of targeted buffer : GH

    auto ta1_z_yyzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 882);

    auto ta1_z_yyzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 883);

    auto ta1_z_yyzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 884);

    auto ta1_z_yyzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 885);

    auto ta1_z_yyzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 886);

    auto ta1_z_yyzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 887);

    auto ta1_z_yyzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 888);

    auto ta1_z_yyzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 889);

    auto ta1_z_yyzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 890);

    auto ta1_z_yyzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 891);

    auto ta1_z_yyzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 892);

    auto ta1_z_yyzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 893);

    auto ta1_z_yyzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 894);

    auto ta1_z_yyzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 895);

    auto ta1_z_yyzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 896);

    auto ta1_z_yyzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 897);

    auto ta1_z_yyzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 898);

    auto ta1_z_yyzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 899);

    auto ta1_z_yyzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 900);

    auto ta1_z_yyzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 901);

    auto ta1_z_yyzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 902);

#pragma omp simd aligned(pa_y,                   \
                             pa_z,               \
                             pc_y,               \
                             pc_z,               \
                             ta1_z_yy_xxxxy_0,   \
                             ta1_z_yy_xxxxy_1,   \
                             ta1_z_yy_xxxyy_0,   \
                             ta1_z_yy_xxxyy_1,   \
                             ta1_z_yy_xxyyy_0,   \
                             ta1_z_yy_xxyyy_1,   \
                             ta1_z_yy_xyyyy_0,   \
                             ta1_z_yy_xyyyy_1,   \
                             ta1_z_yy_yyyyy_0,   \
                             ta1_z_yy_yyyyy_1,   \
                             ta1_z_yyz_xxxxy_0,  \
                             ta1_z_yyz_xxxxy_1,  \
                             ta1_z_yyz_xxxyy_0,  \
                             ta1_z_yyz_xxxyy_1,  \
                             ta1_z_yyz_xxyyy_0,  \
                             ta1_z_yyz_xxyyy_1,  \
                             ta1_z_yyz_xyyyy_0,  \
                             ta1_z_yyz_xyyyy_1,  \
                             ta1_z_yyz_yyyyy_0,  \
                             ta1_z_yyz_yyyyy_1,  \
                             ta1_z_yyzz_xxxxx_0, \
                             ta1_z_yyzz_xxxxy_0, \
                             ta1_z_yyzz_xxxxz_0, \
                             ta1_z_yyzz_xxxyy_0, \
                             ta1_z_yyzz_xxxyz_0, \
                             ta1_z_yyzz_xxxzz_0, \
                             ta1_z_yyzz_xxyyy_0, \
                             ta1_z_yyzz_xxyyz_0, \
                             ta1_z_yyzz_xxyzz_0, \
                             ta1_z_yyzz_xxzzz_0, \
                             ta1_z_yyzz_xyyyy_0, \
                             ta1_z_yyzz_xyyyz_0, \
                             ta1_z_yyzz_xyyzz_0, \
                             ta1_z_yyzz_xyzzz_0, \
                             ta1_z_yyzz_xzzzz_0, \
                             ta1_z_yyzz_yyyyy_0, \
                             ta1_z_yyzz_yyyyz_0, \
                             ta1_z_yyzz_yyyzz_0, \
                             ta1_z_yyzz_yyzzz_0, \
                             ta1_z_yyzz_yzzzz_0, \
                             ta1_z_yyzz_zzzzz_0, \
                             ta1_z_yzz_xxxxx_0,  \
                             ta1_z_yzz_xxxxx_1,  \
                             ta1_z_yzz_xxxxz_0,  \
                             ta1_z_yzz_xxxxz_1,  \
                             ta1_z_yzz_xxxyz_0,  \
                             ta1_z_yzz_xxxyz_1,  \
                             ta1_z_yzz_xxxz_0,   \
                             ta1_z_yzz_xxxz_1,   \
                             ta1_z_yzz_xxxzz_0,  \
                             ta1_z_yzz_xxxzz_1,  \
                             ta1_z_yzz_xxyyz_0,  \
                             ta1_z_yzz_xxyyz_1,  \
                             ta1_z_yzz_xxyz_0,   \
                             ta1_z_yzz_xxyz_1,   \
                             ta1_z_yzz_xxyzz_0,  \
                             ta1_z_yzz_xxyzz_1,  \
                             ta1_z_yzz_xxzz_0,   \
                             ta1_z_yzz_xxzz_1,   \
                             ta1_z_yzz_xxzzz_0,  \
                             ta1_z_yzz_xxzzz_1,  \
                             ta1_z_yzz_xyyyz_0,  \
                             ta1_z_yzz_xyyyz_1,  \
                             ta1_z_yzz_xyyz_0,   \
                             ta1_z_yzz_xyyz_1,   \
                             ta1_z_yzz_xyyzz_0,  \
                             ta1_z_yzz_xyyzz_1,  \
                             ta1_z_yzz_xyzz_0,   \
                             ta1_z_yzz_xyzz_1,   \
                             ta1_z_yzz_xyzzz_0,  \
                             ta1_z_yzz_xyzzz_1,  \
                             ta1_z_yzz_xzzz_0,   \
                             ta1_z_yzz_xzzz_1,   \
                             ta1_z_yzz_xzzzz_0,  \
                             ta1_z_yzz_xzzzz_1,  \
                             ta1_z_yzz_yyyyz_0,  \
                             ta1_z_yzz_yyyyz_1,  \
                             ta1_z_yzz_yyyz_0,   \
                             ta1_z_yzz_yyyz_1,   \
                             ta1_z_yzz_yyyzz_0,  \
                             ta1_z_yzz_yyyzz_1,  \
                             ta1_z_yzz_yyzz_0,   \
                             ta1_z_yzz_yyzz_1,   \
                             ta1_z_yzz_yyzzz_0,  \
                             ta1_z_yzz_yyzzz_1,  \
                             ta1_z_yzz_yzzz_0,   \
                             ta1_z_yzz_yzzz_1,   \
                             ta1_z_yzz_yzzzz_0,  \
                             ta1_z_yzz_yzzzz_1,  \
                             ta1_z_yzz_zzzz_0,   \
                             ta1_z_yzz_zzzz_1,   \
                             ta1_z_yzz_zzzzz_0,  \
                             ta1_z_yzz_zzzzz_1,  \
                             ta1_z_zz_xxxxx_0,   \
                             ta1_z_zz_xxxxx_1,   \
                             ta1_z_zz_xxxxz_0,   \
                             ta1_z_zz_xxxxz_1,   \
                             ta1_z_zz_xxxyz_0,   \
                             ta1_z_zz_xxxyz_1,   \
                             ta1_z_zz_xxxzz_0,   \
                             ta1_z_zz_xxxzz_1,   \
                             ta1_z_zz_xxyyz_0,   \
                             ta1_z_zz_xxyyz_1,   \
                             ta1_z_zz_xxyzz_0,   \
                             ta1_z_zz_xxyzz_1,   \
                             ta1_z_zz_xxzzz_0,   \
                             ta1_z_zz_xxzzz_1,   \
                             ta1_z_zz_xyyyz_0,   \
                             ta1_z_zz_xyyyz_1,   \
                             ta1_z_zz_xyyzz_0,   \
                             ta1_z_zz_xyyzz_1,   \
                             ta1_z_zz_xyzzz_0,   \
                             ta1_z_zz_xyzzz_1,   \
                             ta1_z_zz_xzzzz_0,   \
                             ta1_z_zz_xzzzz_1,   \
                             ta1_z_zz_yyyyz_0,   \
                             ta1_z_zz_yyyyz_1,   \
                             ta1_z_zz_yyyzz_0,   \
                             ta1_z_zz_yyyzz_1,   \
                             ta1_z_zz_yyzzz_0,   \
                             ta1_z_zz_yyzzz_1,   \
                             ta1_z_zz_yzzzz_0,   \
                             ta1_z_zz_yzzzz_1,   \
                             ta1_z_zz_zzzzz_0,   \
                             ta1_z_zz_zzzzz_1,   \
                             ta_yyz_xxxxy_1,     \
                             ta_yyz_xxxyy_1,     \
                             ta_yyz_xxyyy_1,     \
                             ta_yyz_xyyyy_1,     \
                             ta_yyz_yyyyy_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yyzz_xxxxx_0[i] =
            ta1_z_zz_xxxxx_0[i] * fe_0 - ta1_z_zz_xxxxx_1[i] * fe_0 + ta1_z_yzz_xxxxx_0[i] * pa_y[i] - ta1_z_yzz_xxxxx_1[i] * pc_y[i];

        ta1_z_yyzz_xxxxy_0[i] = ta1_z_yy_xxxxy_0[i] * fe_0 - ta1_z_yy_xxxxy_1[i] * fe_0 + ta_yyz_xxxxy_1[i] + ta1_z_yyz_xxxxy_0[i] * pa_z[i] -
                                ta1_z_yyz_xxxxy_1[i] * pc_z[i];

        ta1_z_yyzz_xxxxz_0[i] =
            ta1_z_zz_xxxxz_0[i] * fe_0 - ta1_z_zz_xxxxz_1[i] * fe_0 + ta1_z_yzz_xxxxz_0[i] * pa_y[i] - ta1_z_yzz_xxxxz_1[i] * pc_y[i];

        ta1_z_yyzz_xxxyy_0[i] = ta1_z_yy_xxxyy_0[i] * fe_0 - ta1_z_yy_xxxyy_1[i] * fe_0 + ta_yyz_xxxyy_1[i] + ta1_z_yyz_xxxyy_0[i] * pa_z[i] -
                                ta1_z_yyz_xxxyy_1[i] * pc_z[i];

        ta1_z_yyzz_xxxyz_0[i] = ta1_z_zz_xxxyz_0[i] * fe_0 - ta1_z_zz_xxxyz_1[i] * fe_0 + ta1_z_yzz_xxxz_0[i] * fe_0 - ta1_z_yzz_xxxz_1[i] * fe_0 +
                                ta1_z_yzz_xxxyz_0[i] * pa_y[i] - ta1_z_yzz_xxxyz_1[i] * pc_y[i];

        ta1_z_yyzz_xxxzz_0[i] =
            ta1_z_zz_xxxzz_0[i] * fe_0 - ta1_z_zz_xxxzz_1[i] * fe_0 + ta1_z_yzz_xxxzz_0[i] * pa_y[i] - ta1_z_yzz_xxxzz_1[i] * pc_y[i];

        ta1_z_yyzz_xxyyy_0[i] = ta1_z_yy_xxyyy_0[i] * fe_0 - ta1_z_yy_xxyyy_1[i] * fe_0 + ta_yyz_xxyyy_1[i] + ta1_z_yyz_xxyyy_0[i] * pa_z[i] -
                                ta1_z_yyz_xxyyy_1[i] * pc_z[i];

        ta1_z_yyzz_xxyyz_0[i] = ta1_z_zz_xxyyz_0[i] * fe_0 - ta1_z_zz_xxyyz_1[i] * fe_0 + 2.0 * ta1_z_yzz_xxyz_0[i] * fe_0 -
                                2.0 * ta1_z_yzz_xxyz_1[i] * fe_0 + ta1_z_yzz_xxyyz_0[i] * pa_y[i] - ta1_z_yzz_xxyyz_1[i] * pc_y[i];

        ta1_z_yyzz_xxyzz_0[i] = ta1_z_zz_xxyzz_0[i] * fe_0 - ta1_z_zz_xxyzz_1[i] * fe_0 + ta1_z_yzz_xxzz_0[i] * fe_0 - ta1_z_yzz_xxzz_1[i] * fe_0 +
                                ta1_z_yzz_xxyzz_0[i] * pa_y[i] - ta1_z_yzz_xxyzz_1[i] * pc_y[i];

        ta1_z_yyzz_xxzzz_0[i] =
            ta1_z_zz_xxzzz_0[i] * fe_0 - ta1_z_zz_xxzzz_1[i] * fe_0 + ta1_z_yzz_xxzzz_0[i] * pa_y[i] - ta1_z_yzz_xxzzz_1[i] * pc_y[i];

        ta1_z_yyzz_xyyyy_0[i] = ta1_z_yy_xyyyy_0[i] * fe_0 - ta1_z_yy_xyyyy_1[i] * fe_0 + ta_yyz_xyyyy_1[i] + ta1_z_yyz_xyyyy_0[i] * pa_z[i] -
                                ta1_z_yyz_xyyyy_1[i] * pc_z[i];

        ta1_z_yyzz_xyyyz_0[i] = ta1_z_zz_xyyyz_0[i] * fe_0 - ta1_z_zz_xyyyz_1[i] * fe_0 + 3.0 * ta1_z_yzz_xyyz_0[i] * fe_0 -
                                3.0 * ta1_z_yzz_xyyz_1[i] * fe_0 + ta1_z_yzz_xyyyz_0[i] * pa_y[i] - ta1_z_yzz_xyyyz_1[i] * pc_y[i];

        ta1_z_yyzz_xyyzz_0[i] = ta1_z_zz_xyyzz_0[i] * fe_0 - ta1_z_zz_xyyzz_1[i] * fe_0 + 2.0 * ta1_z_yzz_xyzz_0[i] * fe_0 -
                                2.0 * ta1_z_yzz_xyzz_1[i] * fe_0 + ta1_z_yzz_xyyzz_0[i] * pa_y[i] - ta1_z_yzz_xyyzz_1[i] * pc_y[i];

        ta1_z_yyzz_xyzzz_0[i] = ta1_z_zz_xyzzz_0[i] * fe_0 - ta1_z_zz_xyzzz_1[i] * fe_0 + ta1_z_yzz_xzzz_0[i] * fe_0 - ta1_z_yzz_xzzz_1[i] * fe_0 +
                                ta1_z_yzz_xyzzz_0[i] * pa_y[i] - ta1_z_yzz_xyzzz_1[i] * pc_y[i];

        ta1_z_yyzz_xzzzz_0[i] =
            ta1_z_zz_xzzzz_0[i] * fe_0 - ta1_z_zz_xzzzz_1[i] * fe_0 + ta1_z_yzz_xzzzz_0[i] * pa_y[i] - ta1_z_yzz_xzzzz_1[i] * pc_y[i];

        ta1_z_yyzz_yyyyy_0[i] = ta1_z_yy_yyyyy_0[i] * fe_0 - ta1_z_yy_yyyyy_1[i] * fe_0 + ta_yyz_yyyyy_1[i] + ta1_z_yyz_yyyyy_0[i] * pa_z[i] -
                                ta1_z_yyz_yyyyy_1[i] * pc_z[i];

        ta1_z_yyzz_yyyyz_0[i] = ta1_z_zz_yyyyz_0[i] * fe_0 - ta1_z_zz_yyyyz_1[i] * fe_0 + 4.0 * ta1_z_yzz_yyyz_0[i] * fe_0 -
                                4.0 * ta1_z_yzz_yyyz_1[i] * fe_0 + ta1_z_yzz_yyyyz_0[i] * pa_y[i] - ta1_z_yzz_yyyyz_1[i] * pc_y[i];

        ta1_z_yyzz_yyyzz_0[i] = ta1_z_zz_yyyzz_0[i] * fe_0 - ta1_z_zz_yyyzz_1[i] * fe_0 + 3.0 * ta1_z_yzz_yyzz_0[i] * fe_0 -
                                3.0 * ta1_z_yzz_yyzz_1[i] * fe_0 + ta1_z_yzz_yyyzz_0[i] * pa_y[i] - ta1_z_yzz_yyyzz_1[i] * pc_y[i];

        ta1_z_yyzz_yyzzz_0[i] = ta1_z_zz_yyzzz_0[i] * fe_0 - ta1_z_zz_yyzzz_1[i] * fe_0 + 2.0 * ta1_z_yzz_yzzz_0[i] * fe_0 -
                                2.0 * ta1_z_yzz_yzzz_1[i] * fe_0 + ta1_z_yzz_yyzzz_0[i] * pa_y[i] - ta1_z_yzz_yyzzz_1[i] * pc_y[i];

        ta1_z_yyzz_yzzzz_0[i] = ta1_z_zz_yzzzz_0[i] * fe_0 - ta1_z_zz_yzzzz_1[i] * fe_0 + ta1_z_yzz_zzzz_0[i] * fe_0 - ta1_z_yzz_zzzz_1[i] * fe_0 +
                                ta1_z_yzz_yzzzz_0[i] * pa_y[i] - ta1_z_yzz_yzzzz_1[i] * pc_y[i];

        ta1_z_yyzz_zzzzz_0[i] =
            ta1_z_zz_zzzzz_0[i] * fe_0 - ta1_z_zz_zzzzz_1[i] * fe_0 + ta1_z_yzz_zzzzz_0[i] * pa_y[i] - ta1_z_yzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 903-924 components of targeted buffer : GH

    auto ta1_z_yzzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 903);

    auto ta1_z_yzzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 904);

    auto ta1_z_yzzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 905);

    auto ta1_z_yzzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 906);

    auto ta1_z_yzzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 907);

    auto ta1_z_yzzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 908);

    auto ta1_z_yzzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 909);

    auto ta1_z_yzzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 910);

    auto ta1_z_yzzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 911);

    auto ta1_z_yzzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 912);

    auto ta1_z_yzzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 913);

    auto ta1_z_yzzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 914);

    auto ta1_z_yzzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 915);

    auto ta1_z_yzzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 916);

    auto ta1_z_yzzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 917);

    auto ta1_z_yzzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 918);

    auto ta1_z_yzzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 919);

    auto ta1_z_yzzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 920);

    auto ta1_z_yzzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 921);

    auto ta1_z_yzzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 922);

    auto ta1_z_yzzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 923);

#pragma omp simd aligned(pa_y,                   \
                             pc_y,               \
                             ta1_z_yzzz_xxxxx_0, \
                             ta1_z_yzzz_xxxxy_0, \
                             ta1_z_yzzz_xxxxz_0, \
                             ta1_z_yzzz_xxxyy_0, \
                             ta1_z_yzzz_xxxyz_0, \
                             ta1_z_yzzz_xxxzz_0, \
                             ta1_z_yzzz_xxyyy_0, \
                             ta1_z_yzzz_xxyyz_0, \
                             ta1_z_yzzz_xxyzz_0, \
                             ta1_z_yzzz_xxzzz_0, \
                             ta1_z_yzzz_xyyyy_0, \
                             ta1_z_yzzz_xyyyz_0, \
                             ta1_z_yzzz_xyyzz_0, \
                             ta1_z_yzzz_xyzzz_0, \
                             ta1_z_yzzz_xzzzz_0, \
                             ta1_z_yzzz_yyyyy_0, \
                             ta1_z_yzzz_yyyyz_0, \
                             ta1_z_yzzz_yyyzz_0, \
                             ta1_z_yzzz_yyzzz_0, \
                             ta1_z_yzzz_yzzzz_0, \
                             ta1_z_yzzz_zzzzz_0, \
                             ta1_z_zzz_xxxx_0,   \
                             ta1_z_zzz_xxxx_1,   \
                             ta1_z_zzz_xxxxx_0,  \
                             ta1_z_zzz_xxxxx_1,  \
                             ta1_z_zzz_xxxxy_0,  \
                             ta1_z_zzz_xxxxy_1,  \
                             ta1_z_zzz_xxxxz_0,  \
                             ta1_z_zzz_xxxxz_1,  \
                             ta1_z_zzz_xxxy_0,   \
                             ta1_z_zzz_xxxy_1,   \
                             ta1_z_zzz_xxxyy_0,  \
                             ta1_z_zzz_xxxyy_1,  \
                             ta1_z_zzz_xxxyz_0,  \
                             ta1_z_zzz_xxxyz_1,  \
                             ta1_z_zzz_xxxz_0,   \
                             ta1_z_zzz_xxxz_1,   \
                             ta1_z_zzz_xxxzz_0,  \
                             ta1_z_zzz_xxxzz_1,  \
                             ta1_z_zzz_xxyy_0,   \
                             ta1_z_zzz_xxyy_1,   \
                             ta1_z_zzz_xxyyy_0,  \
                             ta1_z_zzz_xxyyy_1,  \
                             ta1_z_zzz_xxyyz_0,  \
                             ta1_z_zzz_xxyyz_1,  \
                             ta1_z_zzz_xxyz_0,   \
                             ta1_z_zzz_xxyz_1,   \
                             ta1_z_zzz_xxyzz_0,  \
                             ta1_z_zzz_xxyzz_1,  \
                             ta1_z_zzz_xxzz_0,   \
                             ta1_z_zzz_xxzz_1,   \
                             ta1_z_zzz_xxzzz_0,  \
                             ta1_z_zzz_xxzzz_1,  \
                             ta1_z_zzz_xyyy_0,   \
                             ta1_z_zzz_xyyy_1,   \
                             ta1_z_zzz_xyyyy_0,  \
                             ta1_z_zzz_xyyyy_1,  \
                             ta1_z_zzz_xyyyz_0,  \
                             ta1_z_zzz_xyyyz_1,  \
                             ta1_z_zzz_xyyz_0,   \
                             ta1_z_zzz_xyyz_1,   \
                             ta1_z_zzz_xyyzz_0,  \
                             ta1_z_zzz_xyyzz_1,  \
                             ta1_z_zzz_xyzz_0,   \
                             ta1_z_zzz_xyzz_1,   \
                             ta1_z_zzz_xyzzz_0,  \
                             ta1_z_zzz_xyzzz_1,  \
                             ta1_z_zzz_xzzz_0,   \
                             ta1_z_zzz_xzzz_1,   \
                             ta1_z_zzz_xzzzz_0,  \
                             ta1_z_zzz_xzzzz_1,  \
                             ta1_z_zzz_yyyy_0,   \
                             ta1_z_zzz_yyyy_1,   \
                             ta1_z_zzz_yyyyy_0,  \
                             ta1_z_zzz_yyyyy_1,  \
                             ta1_z_zzz_yyyyz_0,  \
                             ta1_z_zzz_yyyyz_1,  \
                             ta1_z_zzz_yyyz_0,   \
                             ta1_z_zzz_yyyz_1,   \
                             ta1_z_zzz_yyyzz_0,  \
                             ta1_z_zzz_yyyzz_1,  \
                             ta1_z_zzz_yyzz_0,   \
                             ta1_z_zzz_yyzz_1,   \
                             ta1_z_zzz_yyzzz_0,  \
                             ta1_z_zzz_yyzzz_1,  \
                             ta1_z_zzz_yzzz_0,   \
                             ta1_z_zzz_yzzz_1,   \
                             ta1_z_zzz_yzzzz_0,  \
                             ta1_z_zzz_yzzzz_1,  \
                             ta1_z_zzz_zzzz_0,   \
                             ta1_z_zzz_zzzz_1,   \
                             ta1_z_zzz_zzzzz_0,  \
                             ta1_z_zzz_zzzzz_1,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_yzzz_xxxxx_0[i] = ta1_z_zzz_xxxxx_0[i] * pa_y[i] - ta1_z_zzz_xxxxx_1[i] * pc_y[i];

        ta1_z_yzzz_xxxxy_0[i] =
            ta1_z_zzz_xxxx_0[i] * fe_0 - ta1_z_zzz_xxxx_1[i] * fe_0 + ta1_z_zzz_xxxxy_0[i] * pa_y[i] - ta1_z_zzz_xxxxy_1[i] * pc_y[i];

        ta1_z_yzzz_xxxxz_0[i] = ta1_z_zzz_xxxxz_0[i] * pa_y[i] - ta1_z_zzz_xxxxz_1[i] * pc_y[i];

        ta1_z_yzzz_xxxyy_0[i] =
            2.0 * ta1_z_zzz_xxxy_0[i] * fe_0 - 2.0 * ta1_z_zzz_xxxy_1[i] * fe_0 + ta1_z_zzz_xxxyy_0[i] * pa_y[i] - ta1_z_zzz_xxxyy_1[i] * pc_y[i];

        ta1_z_yzzz_xxxyz_0[i] =
            ta1_z_zzz_xxxz_0[i] * fe_0 - ta1_z_zzz_xxxz_1[i] * fe_0 + ta1_z_zzz_xxxyz_0[i] * pa_y[i] - ta1_z_zzz_xxxyz_1[i] * pc_y[i];

        ta1_z_yzzz_xxxzz_0[i] = ta1_z_zzz_xxxzz_0[i] * pa_y[i] - ta1_z_zzz_xxxzz_1[i] * pc_y[i];

        ta1_z_yzzz_xxyyy_0[i] =
            3.0 * ta1_z_zzz_xxyy_0[i] * fe_0 - 3.0 * ta1_z_zzz_xxyy_1[i] * fe_0 + ta1_z_zzz_xxyyy_0[i] * pa_y[i] - ta1_z_zzz_xxyyy_1[i] * pc_y[i];

        ta1_z_yzzz_xxyyz_0[i] =
            2.0 * ta1_z_zzz_xxyz_0[i] * fe_0 - 2.0 * ta1_z_zzz_xxyz_1[i] * fe_0 + ta1_z_zzz_xxyyz_0[i] * pa_y[i] - ta1_z_zzz_xxyyz_1[i] * pc_y[i];

        ta1_z_yzzz_xxyzz_0[i] =
            ta1_z_zzz_xxzz_0[i] * fe_0 - ta1_z_zzz_xxzz_1[i] * fe_0 + ta1_z_zzz_xxyzz_0[i] * pa_y[i] - ta1_z_zzz_xxyzz_1[i] * pc_y[i];

        ta1_z_yzzz_xxzzz_0[i] = ta1_z_zzz_xxzzz_0[i] * pa_y[i] - ta1_z_zzz_xxzzz_1[i] * pc_y[i];

        ta1_z_yzzz_xyyyy_0[i] =
            4.0 * ta1_z_zzz_xyyy_0[i] * fe_0 - 4.0 * ta1_z_zzz_xyyy_1[i] * fe_0 + ta1_z_zzz_xyyyy_0[i] * pa_y[i] - ta1_z_zzz_xyyyy_1[i] * pc_y[i];

        ta1_z_yzzz_xyyyz_0[i] =
            3.0 * ta1_z_zzz_xyyz_0[i] * fe_0 - 3.0 * ta1_z_zzz_xyyz_1[i] * fe_0 + ta1_z_zzz_xyyyz_0[i] * pa_y[i] - ta1_z_zzz_xyyyz_1[i] * pc_y[i];

        ta1_z_yzzz_xyyzz_0[i] =
            2.0 * ta1_z_zzz_xyzz_0[i] * fe_0 - 2.0 * ta1_z_zzz_xyzz_1[i] * fe_0 + ta1_z_zzz_xyyzz_0[i] * pa_y[i] - ta1_z_zzz_xyyzz_1[i] * pc_y[i];

        ta1_z_yzzz_xyzzz_0[i] =
            ta1_z_zzz_xzzz_0[i] * fe_0 - ta1_z_zzz_xzzz_1[i] * fe_0 + ta1_z_zzz_xyzzz_0[i] * pa_y[i] - ta1_z_zzz_xyzzz_1[i] * pc_y[i];

        ta1_z_yzzz_xzzzz_0[i] = ta1_z_zzz_xzzzz_0[i] * pa_y[i] - ta1_z_zzz_xzzzz_1[i] * pc_y[i];

        ta1_z_yzzz_yyyyy_0[i] =
            5.0 * ta1_z_zzz_yyyy_0[i] * fe_0 - 5.0 * ta1_z_zzz_yyyy_1[i] * fe_0 + ta1_z_zzz_yyyyy_0[i] * pa_y[i] - ta1_z_zzz_yyyyy_1[i] * pc_y[i];

        ta1_z_yzzz_yyyyz_0[i] =
            4.0 * ta1_z_zzz_yyyz_0[i] * fe_0 - 4.0 * ta1_z_zzz_yyyz_1[i] * fe_0 + ta1_z_zzz_yyyyz_0[i] * pa_y[i] - ta1_z_zzz_yyyyz_1[i] * pc_y[i];

        ta1_z_yzzz_yyyzz_0[i] =
            3.0 * ta1_z_zzz_yyzz_0[i] * fe_0 - 3.0 * ta1_z_zzz_yyzz_1[i] * fe_0 + ta1_z_zzz_yyyzz_0[i] * pa_y[i] - ta1_z_zzz_yyyzz_1[i] * pc_y[i];

        ta1_z_yzzz_yyzzz_0[i] =
            2.0 * ta1_z_zzz_yzzz_0[i] * fe_0 - 2.0 * ta1_z_zzz_yzzz_1[i] * fe_0 + ta1_z_zzz_yyzzz_0[i] * pa_y[i] - ta1_z_zzz_yyzzz_1[i] * pc_y[i];

        ta1_z_yzzz_yzzzz_0[i] =
            ta1_z_zzz_zzzz_0[i] * fe_0 - ta1_z_zzz_zzzz_1[i] * fe_0 + ta1_z_zzz_yzzzz_0[i] * pa_y[i] - ta1_z_zzz_yzzzz_1[i] * pc_y[i];

        ta1_z_yzzz_zzzzz_0[i] = ta1_z_zzz_zzzzz_0[i] * pa_y[i] - ta1_z_zzz_zzzzz_1[i] * pc_y[i];
    }

    // Set up 924-945 components of targeted buffer : GH

    auto ta1_z_zzzz_xxxxx_0 = pbuffer.data(idx_npot_geom_010_0_gh + 924);

    auto ta1_z_zzzz_xxxxy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 925);

    auto ta1_z_zzzz_xxxxz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 926);

    auto ta1_z_zzzz_xxxyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 927);

    auto ta1_z_zzzz_xxxyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 928);

    auto ta1_z_zzzz_xxxzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 929);

    auto ta1_z_zzzz_xxyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 930);

    auto ta1_z_zzzz_xxyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 931);

    auto ta1_z_zzzz_xxyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 932);

    auto ta1_z_zzzz_xxzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 933);

    auto ta1_z_zzzz_xyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 934);

    auto ta1_z_zzzz_xyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 935);

    auto ta1_z_zzzz_xyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 936);

    auto ta1_z_zzzz_xyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 937);

    auto ta1_z_zzzz_xzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 938);

    auto ta1_z_zzzz_yyyyy_0 = pbuffer.data(idx_npot_geom_010_0_gh + 939);

    auto ta1_z_zzzz_yyyyz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 940);

    auto ta1_z_zzzz_yyyzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 941);

    auto ta1_z_zzzz_yyzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 942);

    auto ta1_z_zzzz_yzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 943);

    auto ta1_z_zzzz_zzzzz_0 = pbuffer.data(idx_npot_geom_010_0_gh + 944);

#pragma omp simd aligned(pa_z,                   \
                             pc_z,               \
                             ta1_z_zz_xxxxx_0,   \
                             ta1_z_zz_xxxxx_1,   \
                             ta1_z_zz_xxxxy_0,   \
                             ta1_z_zz_xxxxy_1,   \
                             ta1_z_zz_xxxxz_0,   \
                             ta1_z_zz_xxxxz_1,   \
                             ta1_z_zz_xxxyy_0,   \
                             ta1_z_zz_xxxyy_1,   \
                             ta1_z_zz_xxxyz_0,   \
                             ta1_z_zz_xxxyz_1,   \
                             ta1_z_zz_xxxzz_0,   \
                             ta1_z_zz_xxxzz_1,   \
                             ta1_z_zz_xxyyy_0,   \
                             ta1_z_zz_xxyyy_1,   \
                             ta1_z_zz_xxyyz_0,   \
                             ta1_z_zz_xxyyz_1,   \
                             ta1_z_zz_xxyzz_0,   \
                             ta1_z_zz_xxyzz_1,   \
                             ta1_z_zz_xxzzz_0,   \
                             ta1_z_zz_xxzzz_1,   \
                             ta1_z_zz_xyyyy_0,   \
                             ta1_z_zz_xyyyy_1,   \
                             ta1_z_zz_xyyyz_0,   \
                             ta1_z_zz_xyyyz_1,   \
                             ta1_z_zz_xyyzz_0,   \
                             ta1_z_zz_xyyzz_1,   \
                             ta1_z_zz_xyzzz_0,   \
                             ta1_z_zz_xyzzz_1,   \
                             ta1_z_zz_xzzzz_0,   \
                             ta1_z_zz_xzzzz_1,   \
                             ta1_z_zz_yyyyy_0,   \
                             ta1_z_zz_yyyyy_1,   \
                             ta1_z_zz_yyyyz_0,   \
                             ta1_z_zz_yyyyz_1,   \
                             ta1_z_zz_yyyzz_0,   \
                             ta1_z_zz_yyyzz_1,   \
                             ta1_z_zz_yyzzz_0,   \
                             ta1_z_zz_yyzzz_1,   \
                             ta1_z_zz_yzzzz_0,   \
                             ta1_z_zz_yzzzz_1,   \
                             ta1_z_zz_zzzzz_0,   \
                             ta1_z_zz_zzzzz_1,   \
                             ta1_z_zzz_xxxx_0,   \
                             ta1_z_zzz_xxxx_1,   \
                             ta1_z_zzz_xxxxx_0,  \
                             ta1_z_zzz_xxxxx_1,  \
                             ta1_z_zzz_xxxxy_0,  \
                             ta1_z_zzz_xxxxy_1,  \
                             ta1_z_zzz_xxxxz_0,  \
                             ta1_z_zzz_xxxxz_1,  \
                             ta1_z_zzz_xxxy_0,   \
                             ta1_z_zzz_xxxy_1,   \
                             ta1_z_zzz_xxxyy_0,  \
                             ta1_z_zzz_xxxyy_1,  \
                             ta1_z_zzz_xxxyz_0,  \
                             ta1_z_zzz_xxxyz_1,  \
                             ta1_z_zzz_xxxz_0,   \
                             ta1_z_zzz_xxxz_1,   \
                             ta1_z_zzz_xxxzz_0,  \
                             ta1_z_zzz_xxxzz_1,  \
                             ta1_z_zzz_xxyy_0,   \
                             ta1_z_zzz_xxyy_1,   \
                             ta1_z_zzz_xxyyy_0,  \
                             ta1_z_zzz_xxyyy_1,  \
                             ta1_z_zzz_xxyyz_0,  \
                             ta1_z_zzz_xxyyz_1,  \
                             ta1_z_zzz_xxyz_0,   \
                             ta1_z_zzz_xxyz_1,   \
                             ta1_z_zzz_xxyzz_0,  \
                             ta1_z_zzz_xxyzz_1,  \
                             ta1_z_zzz_xxzz_0,   \
                             ta1_z_zzz_xxzz_1,   \
                             ta1_z_zzz_xxzzz_0,  \
                             ta1_z_zzz_xxzzz_1,  \
                             ta1_z_zzz_xyyy_0,   \
                             ta1_z_zzz_xyyy_1,   \
                             ta1_z_zzz_xyyyy_0,  \
                             ta1_z_zzz_xyyyy_1,  \
                             ta1_z_zzz_xyyyz_0,  \
                             ta1_z_zzz_xyyyz_1,  \
                             ta1_z_zzz_xyyz_0,   \
                             ta1_z_zzz_xyyz_1,   \
                             ta1_z_zzz_xyyzz_0,  \
                             ta1_z_zzz_xyyzz_1,  \
                             ta1_z_zzz_xyzz_0,   \
                             ta1_z_zzz_xyzz_1,   \
                             ta1_z_zzz_xyzzz_0,  \
                             ta1_z_zzz_xyzzz_1,  \
                             ta1_z_zzz_xzzz_0,   \
                             ta1_z_zzz_xzzz_1,   \
                             ta1_z_zzz_xzzzz_0,  \
                             ta1_z_zzz_xzzzz_1,  \
                             ta1_z_zzz_yyyy_0,   \
                             ta1_z_zzz_yyyy_1,   \
                             ta1_z_zzz_yyyyy_0,  \
                             ta1_z_zzz_yyyyy_1,  \
                             ta1_z_zzz_yyyyz_0,  \
                             ta1_z_zzz_yyyyz_1,  \
                             ta1_z_zzz_yyyz_0,   \
                             ta1_z_zzz_yyyz_1,   \
                             ta1_z_zzz_yyyzz_0,  \
                             ta1_z_zzz_yyyzz_1,  \
                             ta1_z_zzz_yyzz_0,   \
                             ta1_z_zzz_yyzz_1,   \
                             ta1_z_zzz_yyzzz_0,  \
                             ta1_z_zzz_yyzzz_1,  \
                             ta1_z_zzz_yzzz_0,   \
                             ta1_z_zzz_yzzz_1,   \
                             ta1_z_zzz_yzzzz_0,  \
                             ta1_z_zzz_yzzzz_1,  \
                             ta1_z_zzz_zzzz_0,   \
                             ta1_z_zzz_zzzz_1,   \
                             ta1_z_zzz_zzzzz_0,  \
                             ta1_z_zzz_zzzzz_1,  \
                             ta1_z_zzzz_xxxxx_0, \
                             ta1_z_zzzz_xxxxy_0, \
                             ta1_z_zzzz_xxxxz_0, \
                             ta1_z_zzzz_xxxyy_0, \
                             ta1_z_zzzz_xxxyz_0, \
                             ta1_z_zzzz_xxxzz_0, \
                             ta1_z_zzzz_xxyyy_0, \
                             ta1_z_zzzz_xxyyz_0, \
                             ta1_z_zzzz_xxyzz_0, \
                             ta1_z_zzzz_xxzzz_0, \
                             ta1_z_zzzz_xyyyy_0, \
                             ta1_z_zzzz_xyyyz_0, \
                             ta1_z_zzzz_xyyzz_0, \
                             ta1_z_zzzz_xyzzz_0, \
                             ta1_z_zzzz_xzzzz_0, \
                             ta1_z_zzzz_yyyyy_0, \
                             ta1_z_zzzz_yyyyz_0, \
                             ta1_z_zzzz_yyyzz_0, \
                             ta1_z_zzzz_yyzzz_0, \
                             ta1_z_zzzz_yzzzz_0, \
                             ta1_z_zzzz_zzzzz_0, \
                             ta_zzz_xxxxx_1,     \
                             ta_zzz_xxxxy_1,     \
                             ta_zzz_xxxxz_1,     \
                             ta_zzz_xxxyy_1,     \
                             ta_zzz_xxxyz_1,     \
                             ta_zzz_xxxzz_1,     \
                             ta_zzz_xxyyy_1,     \
                             ta_zzz_xxyyz_1,     \
                             ta_zzz_xxyzz_1,     \
                             ta_zzz_xxzzz_1,     \
                             ta_zzz_xyyyy_1,     \
                             ta_zzz_xyyyz_1,     \
                             ta_zzz_xyyzz_1,     \
                             ta_zzz_xyzzz_1,     \
                             ta_zzz_xzzzz_1,     \
                             ta_zzz_yyyyy_1,     \
                             ta_zzz_yyyyz_1,     \
                             ta_zzz_yyyzz_1,     \
                             ta_zzz_yyzzz_1,     \
                             ta_zzz_yzzzz_1,     \
                             ta_zzz_zzzzz_1,     \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        ta1_z_zzzz_xxxxx_0[i] = 3.0 * ta1_z_zz_xxxxx_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxxx_1[i] * fe_0 + ta_zzz_xxxxx_1[i] +
                                ta1_z_zzz_xxxxx_0[i] * pa_z[i] - ta1_z_zzz_xxxxx_1[i] * pc_z[i];

        ta1_z_zzzz_xxxxy_0[i] = 3.0 * ta1_z_zz_xxxxy_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxxy_1[i] * fe_0 + ta_zzz_xxxxy_1[i] +
                                ta1_z_zzz_xxxxy_0[i] * pa_z[i] - ta1_z_zzz_xxxxy_1[i] * pc_z[i];

        ta1_z_zzzz_xxxxz_0[i] = 3.0 * ta1_z_zz_xxxxz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxxz_1[i] * fe_0 + ta1_z_zzz_xxxx_0[i] * fe_0 -
                                ta1_z_zzz_xxxx_1[i] * fe_0 + ta_zzz_xxxxz_1[i] + ta1_z_zzz_xxxxz_0[i] * pa_z[i] - ta1_z_zzz_xxxxz_1[i] * pc_z[i];

        ta1_z_zzzz_xxxyy_0[i] = 3.0 * ta1_z_zz_xxxyy_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxyy_1[i] * fe_0 + ta_zzz_xxxyy_1[i] +
                                ta1_z_zzz_xxxyy_0[i] * pa_z[i] - ta1_z_zzz_xxxyy_1[i] * pc_z[i];

        ta1_z_zzzz_xxxyz_0[i] = 3.0 * ta1_z_zz_xxxyz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxyz_1[i] * fe_0 + ta1_z_zzz_xxxy_0[i] * fe_0 -
                                ta1_z_zzz_xxxy_1[i] * fe_0 + ta_zzz_xxxyz_1[i] + ta1_z_zzz_xxxyz_0[i] * pa_z[i] - ta1_z_zzz_xxxyz_1[i] * pc_z[i];

        ta1_z_zzzz_xxxzz_0[i] = 3.0 * ta1_z_zz_xxxzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxxzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_xxxz_0[i] * fe_0 -
                                2.0 * ta1_z_zzz_xxxz_1[i] * fe_0 + ta_zzz_xxxzz_1[i] + ta1_z_zzz_xxxzz_0[i] * pa_z[i] -
                                ta1_z_zzz_xxxzz_1[i] * pc_z[i];

        ta1_z_zzzz_xxyyy_0[i] = 3.0 * ta1_z_zz_xxyyy_0[i] * fe_0 - 3.0 * ta1_z_zz_xxyyy_1[i] * fe_0 + ta_zzz_xxyyy_1[i] +
                                ta1_z_zzz_xxyyy_0[i] * pa_z[i] - ta1_z_zzz_xxyyy_1[i] * pc_z[i];

        ta1_z_zzzz_xxyyz_0[i] = 3.0 * ta1_z_zz_xxyyz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxyyz_1[i] * fe_0 + ta1_z_zzz_xxyy_0[i] * fe_0 -
                                ta1_z_zzz_xxyy_1[i] * fe_0 + ta_zzz_xxyyz_1[i] + ta1_z_zzz_xxyyz_0[i] * pa_z[i] - ta1_z_zzz_xxyyz_1[i] * pc_z[i];

        ta1_z_zzzz_xxyzz_0[i] = 3.0 * ta1_z_zz_xxyzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxyzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_xxyz_0[i] * fe_0 -
                                2.0 * ta1_z_zzz_xxyz_1[i] * fe_0 + ta_zzz_xxyzz_1[i] + ta1_z_zzz_xxyzz_0[i] * pa_z[i] -
                                ta1_z_zzz_xxyzz_1[i] * pc_z[i];

        ta1_z_zzzz_xxzzz_0[i] = 3.0 * ta1_z_zz_xxzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xxzzz_1[i] * fe_0 + 3.0 * ta1_z_zzz_xxzz_0[i] * fe_0 -
                                3.0 * ta1_z_zzz_xxzz_1[i] * fe_0 + ta_zzz_xxzzz_1[i] + ta1_z_zzz_xxzzz_0[i] * pa_z[i] -
                                ta1_z_zzz_xxzzz_1[i] * pc_z[i];

        ta1_z_zzzz_xyyyy_0[i] = 3.0 * ta1_z_zz_xyyyy_0[i] * fe_0 - 3.0 * ta1_z_zz_xyyyy_1[i] * fe_0 + ta_zzz_xyyyy_1[i] +
                                ta1_z_zzz_xyyyy_0[i] * pa_z[i] - ta1_z_zzz_xyyyy_1[i] * pc_z[i];

        ta1_z_zzzz_xyyyz_0[i] = 3.0 * ta1_z_zz_xyyyz_0[i] * fe_0 - 3.0 * ta1_z_zz_xyyyz_1[i] * fe_0 + ta1_z_zzz_xyyy_0[i] * fe_0 -
                                ta1_z_zzz_xyyy_1[i] * fe_0 + ta_zzz_xyyyz_1[i] + ta1_z_zzz_xyyyz_0[i] * pa_z[i] - ta1_z_zzz_xyyyz_1[i] * pc_z[i];

        ta1_z_zzzz_xyyzz_0[i] = 3.0 * ta1_z_zz_xyyzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xyyzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_xyyz_0[i] * fe_0 -
                                2.0 * ta1_z_zzz_xyyz_1[i] * fe_0 + ta_zzz_xyyzz_1[i] + ta1_z_zzz_xyyzz_0[i] * pa_z[i] -
                                ta1_z_zzz_xyyzz_1[i] * pc_z[i];

        ta1_z_zzzz_xyzzz_0[i] = 3.0 * ta1_z_zz_xyzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xyzzz_1[i] * fe_0 + 3.0 * ta1_z_zzz_xyzz_0[i] * fe_0 -
                                3.0 * ta1_z_zzz_xyzz_1[i] * fe_0 + ta_zzz_xyzzz_1[i] + ta1_z_zzz_xyzzz_0[i] * pa_z[i] -
                                ta1_z_zzz_xyzzz_1[i] * pc_z[i];

        ta1_z_zzzz_xzzzz_0[i] = 3.0 * ta1_z_zz_xzzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_xzzzz_1[i] * fe_0 + 4.0 * ta1_z_zzz_xzzz_0[i] * fe_0 -
                                4.0 * ta1_z_zzz_xzzz_1[i] * fe_0 + ta_zzz_xzzzz_1[i] + ta1_z_zzz_xzzzz_0[i] * pa_z[i] -
                                ta1_z_zzz_xzzzz_1[i] * pc_z[i];

        ta1_z_zzzz_yyyyy_0[i] = 3.0 * ta1_z_zz_yyyyy_0[i] * fe_0 - 3.0 * ta1_z_zz_yyyyy_1[i] * fe_0 + ta_zzz_yyyyy_1[i] +
                                ta1_z_zzz_yyyyy_0[i] * pa_z[i] - ta1_z_zzz_yyyyy_1[i] * pc_z[i];

        ta1_z_zzzz_yyyyz_0[i] = 3.0 * ta1_z_zz_yyyyz_0[i] * fe_0 - 3.0 * ta1_z_zz_yyyyz_1[i] * fe_0 + ta1_z_zzz_yyyy_0[i] * fe_0 -
                                ta1_z_zzz_yyyy_1[i] * fe_0 + ta_zzz_yyyyz_1[i] + ta1_z_zzz_yyyyz_0[i] * pa_z[i] - ta1_z_zzz_yyyyz_1[i] * pc_z[i];

        ta1_z_zzzz_yyyzz_0[i] = 3.0 * ta1_z_zz_yyyzz_0[i] * fe_0 - 3.0 * ta1_z_zz_yyyzz_1[i] * fe_0 + 2.0 * ta1_z_zzz_yyyz_0[i] * fe_0 -
                                2.0 * ta1_z_zzz_yyyz_1[i] * fe_0 + ta_zzz_yyyzz_1[i] + ta1_z_zzz_yyyzz_0[i] * pa_z[i] -
                                ta1_z_zzz_yyyzz_1[i] * pc_z[i];

        ta1_z_zzzz_yyzzz_0[i] = 3.0 * ta1_z_zz_yyzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_yyzzz_1[i] * fe_0 + 3.0 * ta1_z_zzz_yyzz_0[i] * fe_0 -
                                3.0 * ta1_z_zzz_yyzz_1[i] * fe_0 + ta_zzz_yyzzz_1[i] + ta1_z_zzz_yyzzz_0[i] * pa_z[i] -
                                ta1_z_zzz_yyzzz_1[i] * pc_z[i];

        ta1_z_zzzz_yzzzz_0[i] = 3.0 * ta1_z_zz_yzzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_yzzzz_1[i] * fe_0 + 4.0 * ta1_z_zzz_yzzz_0[i] * fe_0 -
                                4.0 * ta1_z_zzz_yzzz_1[i] * fe_0 + ta_zzz_yzzzz_1[i] + ta1_z_zzz_yzzzz_0[i] * pa_z[i] -
                                ta1_z_zzz_yzzzz_1[i] * pc_z[i];

        ta1_z_zzzz_zzzzz_0[i] = 3.0 * ta1_z_zz_zzzzz_0[i] * fe_0 - 3.0 * ta1_z_zz_zzzzz_1[i] * fe_0 + 5.0 * ta1_z_zzz_zzzz_0[i] * fe_0 -
                                5.0 * ta1_z_zzz_zzzz_1[i] * fe_0 + ta_zzz_zzzzz_1[i] + ta1_z_zzz_zzzzz_0[i] * pa_z[i] -
                                ta1_z_zzz_zzzzz_1[i] * pc_z[i];
    }
}

}  // namespace npotrec
