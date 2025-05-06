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

#include "ElectricDipoleMomentumPrimRecHG.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_hg(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_hg,
                                      const size_t              idx_dip_fg,
                                      const size_t              idx_dip_gf,
                                      const size_t              idx_ovl_gg,
                                      const size_t              idx_dip_gg,
                                      const CSimdArray<double>& factors,
                                      const size_t              idx_rpa,
                                      const double              a_exp) -> void
{
    const auto nelems = pbuffer.number_of_active_elements();

    // Set up exponents

    auto b_exps = factors.data(0);

    // Set up R(PA) distances

    auto pa_x = factors.data(idx_rpa);

    auto pa_y = factors.data(idx_rpa + 1);

    auto pa_z = factors.data(idx_rpa + 2);

    // Set up components of auxiliary buffer : FG

    auto tr_x_xxx_xxxx = pbuffer.data(idx_dip_fg);

    auto tr_x_xxx_xxxy = pbuffer.data(idx_dip_fg + 1);

    auto tr_x_xxx_xxxz = pbuffer.data(idx_dip_fg + 2);

    auto tr_x_xxx_xxyy = pbuffer.data(idx_dip_fg + 3);

    auto tr_x_xxx_xxyz = pbuffer.data(idx_dip_fg + 4);

    auto tr_x_xxx_xxzz = pbuffer.data(idx_dip_fg + 5);

    auto tr_x_xxx_xyyy = pbuffer.data(idx_dip_fg + 6);

    auto tr_x_xxx_xyyz = pbuffer.data(idx_dip_fg + 7);

    auto tr_x_xxx_xyzz = pbuffer.data(idx_dip_fg + 8);

    auto tr_x_xxx_xzzz = pbuffer.data(idx_dip_fg + 9);

    auto tr_x_xxx_yyyy = pbuffer.data(idx_dip_fg + 10);

    auto tr_x_xxx_yyyz = pbuffer.data(idx_dip_fg + 11);

    auto tr_x_xxx_yyzz = pbuffer.data(idx_dip_fg + 12);

    auto tr_x_xxx_yzzz = pbuffer.data(idx_dip_fg + 13);

    auto tr_x_xxx_zzzz = pbuffer.data(idx_dip_fg + 14);

    auto tr_x_xxy_xxxx = pbuffer.data(idx_dip_fg + 15);

    auto tr_x_xxy_xxxy = pbuffer.data(idx_dip_fg + 16);

    auto tr_x_xxy_xxxz = pbuffer.data(idx_dip_fg + 17);

    auto tr_x_xxy_xxyy = pbuffer.data(idx_dip_fg + 18);

    auto tr_x_xxy_xxyz = pbuffer.data(idx_dip_fg + 19);

    auto tr_x_xxy_xxzz = pbuffer.data(idx_dip_fg + 20);

    auto tr_x_xxy_xyyy = pbuffer.data(idx_dip_fg + 21);

    auto tr_x_xxy_xyyz = pbuffer.data(idx_dip_fg + 22);

    auto tr_x_xxy_xyzz = pbuffer.data(idx_dip_fg + 23);

    auto tr_x_xxy_xzzz = pbuffer.data(idx_dip_fg + 24);

    auto tr_x_xxy_zzzz = pbuffer.data(idx_dip_fg + 29);

    auto tr_x_xxz_xxxx = pbuffer.data(idx_dip_fg + 30);

    auto tr_x_xxz_xxxy = pbuffer.data(idx_dip_fg + 31);

    auto tr_x_xxz_xxxz = pbuffer.data(idx_dip_fg + 32);

    auto tr_x_xxz_xxyy = pbuffer.data(idx_dip_fg + 33);

    auto tr_x_xxz_xxyz = pbuffer.data(idx_dip_fg + 34);

    auto tr_x_xxz_xxzz = pbuffer.data(idx_dip_fg + 35);

    auto tr_x_xxz_xyyy = pbuffer.data(idx_dip_fg + 36);

    auto tr_x_xxz_xyyz = pbuffer.data(idx_dip_fg + 37);

    auto tr_x_xxz_xyzz = pbuffer.data(idx_dip_fg + 38);

    auto tr_x_xxz_xzzz = pbuffer.data(idx_dip_fg + 39);

    auto tr_x_xxz_yyyy = pbuffer.data(idx_dip_fg + 40);

    auto tr_x_xxz_zzzz = pbuffer.data(idx_dip_fg + 44);

    auto tr_x_xyy_xxxx = pbuffer.data(idx_dip_fg + 45);

    auto tr_x_xyy_xxxy = pbuffer.data(idx_dip_fg + 46);

    auto tr_x_xyy_xxxz = pbuffer.data(idx_dip_fg + 47);

    auto tr_x_xyy_xxyy = pbuffer.data(idx_dip_fg + 48);

    auto tr_x_xyy_xxzz = pbuffer.data(idx_dip_fg + 50);

    auto tr_x_xyy_xyyy = pbuffer.data(idx_dip_fg + 51);

    auto tr_x_xyy_xzzz = pbuffer.data(idx_dip_fg + 54);

    auto tr_x_xyy_yyyy = pbuffer.data(idx_dip_fg + 55);

    auto tr_x_xyy_yyyz = pbuffer.data(idx_dip_fg + 56);

    auto tr_x_xyy_yyzz = pbuffer.data(idx_dip_fg + 57);

    auto tr_x_xyy_yzzz = pbuffer.data(idx_dip_fg + 58);

    auto tr_x_xyz_xxxz = pbuffer.data(idx_dip_fg + 62);

    auto tr_x_xyz_xxzz = pbuffer.data(idx_dip_fg + 65);

    auto tr_x_xyz_xzzz = pbuffer.data(idx_dip_fg + 69);

    auto tr_x_xzz_xxxx = pbuffer.data(idx_dip_fg + 75);

    auto tr_x_xzz_xxxy = pbuffer.data(idx_dip_fg + 76);

    auto tr_x_xzz_xxxz = pbuffer.data(idx_dip_fg + 77);

    auto tr_x_xzz_xxyy = pbuffer.data(idx_dip_fg + 78);

    auto tr_x_xzz_xxzz = pbuffer.data(idx_dip_fg + 80);

    auto tr_x_xzz_xyyy = pbuffer.data(idx_dip_fg + 81);

    auto tr_x_xzz_xzzz = pbuffer.data(idx_dip_fg + 84);

    auto tr_x_xzz_yyyz = pbuffer.data(idx_dip_fg + 86);

    auto tr_x_xzz_yyzz = pbuffer.data(idx_dip_fg + 87);

    auto tr_x_xzz_yzzz = pbuffer.data(idx_dip_fg + 88);

    auto tr_x_xzz_zzzz = pbuffer.data(idx_dip_fg + 89);

    auto tr_x_yyy_xxxx = pbuffer.data(idx_dip_fg + 90);

    auto tr_x_yyy_xxxy = pbuffer.data(idx_dip_fg + 91);

    auto tr_x_yyy_xxxz = pbuffer.data(idx_dip_fg + 92);

    auto tr_x_yyy_xxyy = pbuffer.data(idx_dip_fg + 93);

    auto tr_x_yyy_xxyz = pbuffer.data(idx_dip_fg + 94);

    auto tr_x_yyy_xxzz = pbuffer.data(idx_dip_fg + 95);

    auto tr_x_yyy_xyyy = pbuffer.data(idx_dip_fg + 96);

    auto tr_x_yyy_xyyz = pbuffer.data(idx_dip_fg + 97);

    auto tr_x_yyy_xyzz = pbuffer.data(idx_dip_fg + 98);

    auto tr_x_yyy_xzzz = pbuffer.data(idx_dip_fg + 99);

    auto tr_x_yyy_yyyy = pbuffer.data(idx_dip_fg + 100);

    auto tr_x_yyy_yyyz = pbuffer.data(idx_dip_fg + 101);

    auto tr_x_yyy_yyzz = pbuffer.data(idx_dip_fg + 102);

    auto tr_x_yyy_yzzz = pbuffer.data(idx_dip_fg + 103);

    auto tr_x_yyy_zzzz = pbuffer.data(idx_dip_fg + 104);

    auto tr_x_yyz_xxxy = pbuffer.data(idx_dip_fg + 106);

    auto tr_x_yyz_xxxz = pbuffer.data(idx_dip_fg + 107);

    auto tr_x_yyz_xxyy = pbuffer.data(idx_dip_fg + 108);

    auto tr_x_yyz_xxzz = pbuffer.data(idx_dip_fg + 110);

    auto tr_x_yyz_xyyy = pbuffer.data(idx_dip_fg + 111);

    auto tr_x_yyz_xzzz = pbuffer.data(idx_dip_fg + 114);

    auto tr_x_yyz_yyyy = pbuffer.data(idx_dip_fg + 115);

    auto tr_x_yyz_zzzz = pbuffer.data(idx_dip_fg + 119);

    auto tr_x_yzz_xxxx = pbuffer.data(idx_dip_fg + 120);

    auto tr_x_yzz_xxxz = pbuffer.data(idx_dip_fg + 122);

    auto tr_x_yzz_xxyz = pbuffer.data(idx_dip_fg + 124);

    auto tr_x_yzz_xxzz = pbuffer.data(idx_dip_fg + 125);

    auto tr_x_yzz_xyyz = pbuffer.data(idx_dip_fg + 127);

    auto tr_x_yzz_xyzz = pbuffer.data(idx_dip_fg + 128);

    auto tr_x_yzz_xzzz = pbuffer.data(idx_dip_fg + 129);

    auto tr_x_yzz_yyyz = pbuffer.data(idx_dip_fg + 131);

    auto tr_x_yzz_yyzz = pbuffer.data(idx_dip_fg + 132);

    auto tr_x_yzz_yzzz = pbuffer.data(idx_dip_fg + 133);

    auto tr_x_yzz_zzzz = pbuffer.data(idx_dip_fg + 134);

    auto tr_x_zzz_xxxx = pbuffer.data(idx_dip_fg + 135);

    auto tr_x_zzz_xxxy = pbuffer.data(idx_dip_fg + 136);

    auto tr_x_zzz_xxxz = pbuffer.data(idx_dip_fg + 137);

    auto tr_x_zzz_xxyy = pbuffer.data(idx_dip_fg + 138);

    auto tr_x_zzz_xxyz = pbuffer.data(idx_dip_fg + 139);

    auto tr_x_zzz_xxzz = pbuffer.data(idx_dip_fg + 140);

    auto tr_x_zzz_xyyy = pbuffer.data(idx_dip_fg + 141);

    auto tr_x_zzz_xyyz = pbuffer.data(idx_dip_fg + 142);

    auto tr_x_zzz_xyzz = pbuffer.data(idx_dip_fg + 143);

    auto tr_x_zzz_xzzz = pbuffer.data(idx_dip_fg + 144);

    auto tr_x_zzz_yyyy = pbuffer.data(idx_dip_fg + 145);

    auto tr_x_zzz_yyyz = pbuffer.data(idx_dip_fg + 146);

    auto tr_x_zzz_yyzz = pbuffer.data(idx_dip_fg + 147);

    auto tr_x_zzz_yzzz = pbuffer.data(idx_dip_fg + 148);

    auto tr_x_zzz_zzzz = pbuffer.data(idx_dip_fg + 149);

    auto tr_y_xxx_xxxx = pbuffer.data(idx_dip_fg + 150);

    auto tr_y_xxx_xxxy = pbuffer.data(idx_dip_fg + 151);

    auto tr_y_xxx_xxxz = pbuffer.data(idx_dip_fg + 152);

    auto tr_y_xxx_xxyy = pbuffer.data(idx_dip_fg + 153);

    auto tr_y_xxx_xxyz = pbuffer.data(idx_dip_fg + 154);

    auto tr_y_xxx_xxzz = pbuffer.data(idx_dip_fg + 155);

    auto tr_y_xxx_xyyy = pbuffer.data(idx_dip_fg + 156);

    auto tr_y_xxx_xyyz = pbuffer.data(idx_dip_fg + 157);

    auto tr_y_xxx_xyzz = pbuffer.data(idx_dip_fg + 158);

    auto tr_y_xxx_xzzz = pbuffer.data(idx_dip_fg + 159);

    auto tr_y_xxx_yyyy = pbuffer.data(idx_dip_fg + 160);

    auto tr_y_xxx_yyyz = pbuffer.data(idx_dip_fg + 161);

    auto tr_y_xxx_yyzz = pbuffer.data(idx_dip_fg + 162);

    auto tr_y_xxx_yzzz = pbuffer.data(idx_dip_fg + 163);

    auto tr_y_xxx_zzzz = pbuffer.data(idx_dip_fg + 164);

    auto tr_y_xxy_xxxy = pbuffer.data(idx_dip_fg + 166);

    auto tr_y_xxy_xxyy = pbuffer.data(idx_dip_fg + 168);

    auto tr_y_xxy_xxyz = pbuffer.data(idx_dip_fg + 169);

    auto tr_y_xxy_xyyy = pbuffer.data(idx_dip_fg + 171);

    auto tr_y_xxy_xyyz = pbuffer.data(idx_dip_fg + 172);

    auto tr_y_xxy_xyzz = pbuffer.data(idx_dip_fg + 173);

    auto tr_y_xxy_yyyy = pbuffer.data(idx_dip_fg + 175);

    auto tr_y_xxy_yyyz = pbuffer.data(idx_dip_fg + 176);

    auto tr_y_xxy_yyzz = pbuffer.data(idx_dip_fg + 177);

    auto tr_y_xxy_yzzz = pbuffer.data(idx_dip_fg + 178);

    auto tr_y_xxy_zzzz = pbuffer.data(idx_dip_fg + 179);

    auto tr_y_xxz_xxxx = pbuffer.data(idx_dip_fg + 180);

    auto tr_y_xxz_xxxy = pbuffer.data(idx_dip_fg + 181);

    auto tr_y_xxz_xxyy = pbuffer.data(idx_dip_fg + 183);

    auto tr_y_xxz_xyyy = pbuffer.data(idx_dip_fg + 186);

    auto tr_y_xxz_yyyz = pbuffer.data(idx_dip_fg + 191);

    auto tr_y_xxz_yyzz = pbuffer.data(idx_dip_fg + 192);

    auto tr_y_xxz_yzzz = pbuffer.data(idx_dip_fg + 193);

    auto tr_y_xxz_zzzz = pbuffer.data(idx_dip_fg + 194);

    auto tr_y_xyy_xxxx = pbuffer.data(idx_dip_fg + 195);

    auto tr_y_xyy_xxxy = pbuffer.data(idx_dip_fg + 196);

    auto tr_y_xyy_xxxz = pbuffer.data(idx_dip_fg + 197);

    auto tr_y_xyy_xxyy = pbuffer.data(idx_dip_fg + 198);

    auto tr_y_xyy_xxyz = pbuffer.data(idx_dip_fg + 199);

    auto tr_y_xyy_xxzz = pbuffer.data(idx_dip_fg + 200);

    auto tr_y_xyy_xyyy = pbuffer.data(idx_dip_fg + 201);

    auto tr_y_xyy_xyyz = pbuffer.data(idx_dip_fg + 202);

    auto tr_y_xyy_xyzz = pbuffer.data(idx_dip_fg + 203);

    auto tr_y_xyy_xzzz = pbuffer.data(idx_dip_fg + 204);

    auto tr_y_xyy_yyyy = pbuffer.data(idx_dip_fg + 205);

    auto tr_y_xyy_yyyz = pbuffer.data(idx_dip_fg + 206);

    auto tr_y_xyy_yyzz = pbuffer.data(idx_dip_fg + 207);

    auto tr_y_xyy_yzzz = pbuffer.data(idx_dip_fg + 208);

    auto tr_y_xyy_zzzz = pbuffer.data(idx_dip_fg + 209);

    auto tr_y_xyz_yyyz = pbuffer.data(idx_dip_fg + 221);

    auto tr_y_xyz_yyzz = pbuffer.data(idx_dip_fg + 222);

    auto tr_y_xyz_yzzz = pbuffer.data(idx_dip_fg + 223);

    auto tr_y_xyz_zzzz = pbuffer.data(idx_dip_fg + 224);

    auto tr_y_xzz_xxxz = pbuffer.data(idx_dip_fg + 227);

    auto tr_y_xzz_xxyz = pbuffer.data(idx_dip_fg + 229);

    auto tr_y_xzz_xxzz = pbuffer.data(idx_dip_fg + 230);

    auto tr_y_xzz_xyyz = pbuffer.data(idx_dip_fg + 232);

    auto tr_y_xzz_xyzz = pbuffer.data(idx_dip_fg + 233);

    auto tr_y_xzz_xzzz = pbuffer.data(idx_dip_fg + 234);

    auto tr_y_xzz_yyyy = pbuffer.data(idx_dip_fg + 235);

    auto tr_y_xzz_yyyz = pbuffer.data(idx_dip_fg + 236);

    auto tr_y_xzz_yyzz = pbuffer.data(idx_dip_fg + 237);

    auto tr_y_xzz_yzzz = pbuffer.data(idx_dip_fg + 238);

    auto tr_y_xzz_zzzz = pbuffer.data(idx_dip_fg + 239);

    auto tr_y_yyy_xxxx = pbuffer.data(idx_dip_fg + 240);

    auto tr_y_yyy_xxxy = pbuffer.data(idx_dip_fg + 241);

    auto tr_y_yyy_xxxz = pbuffer.data(idx_dip_fg + 242);

    auto tr_y_yyy_xxyy = pbuffer.data(idx_dip_fg + 243);

    auto tr_y_yyy_xxyz = pbuffer.data(idx_dip_fg + 244);

    auto tr_y_yyy_xxzz = pbuffer.data(idx_dip_fg + 245);

    auto tr_y_yyy_xyyy = pbuffer.data(idx_dip_fg + 246);

    auto tr_y_yyy_xyyz = pbuffer.data(idx_dip_fg + 247);

    auto tr_y_yyy_xyzz = pbuffer.data(idx_dip_fg + 248);

    auto tr_y_yyy_xzzz = pbuffer.data(idx_dip_fg + 249);

    auto tr_y_yyy_yyyy = pbuffer.data(idx_dip_fg + 250);

    auto tr_y_yyy_yyyz = pbuffer.data(idx_dip_fg + 251);

    auto tr_y_yyy_yyzz = pbuffer.data(idx_dip_fg + 252);

    auto tr_y_yyy_yzzz = pbuffer.data(idx_dip_fg + 253);

    auto tr_y_yyy_zzzz = pbuffer.data(idx_dip_fg + 254);

    auto tr_y_yyz_xxxx = pbuffer.data(idx_dip_fg + 255);

    auto tr_y_yyz_xxxy = pbuffer.data(idx_dip_fg + 256);

    auto tr_y_yyz_xxyy = pbuffer.data(idx_dip_fg + 258);

    auto tr_y_yyz_xxyz = pbuffer.data(idx_dip_fg + 259);

    auto tr_y_yyz_xyyy = pbuffer.data(idx_dip_fg + 261);

    auto tr_y_yyz_xyyz = pbuffer.data(idx_dip_fg + 262);

    auto tr_y_yyz_xyzz = pbuffer.data(idx_dip_fg + 263);

    auto tr_y_yyz_yyyy = pbuffer.data(idx_dip_fg + 265);

    auto tr_y_yyz_yyyz = pbuffer.data(idx_dip_fg + 266);

    auto tr_y_yyz_yyzz = pbuffer.data(idx_dip_fg + 267);

    auto tr_y_yyz_yzzz = pbuffer.data(idx_dip_fg + 268);

    auto tr_y_yyz_zzzz = pbuffer.data(idx_dip_fg + 269);

    auto tr_y_yzz_xxxy = pbuffer.data(idx_dip_fg + 271);

    auto tr_y_yzz_xxxz = pbuffer.data(idx_dip_fg + 272);

    auto tr_y_yzz_xxyy = pbuffer.data(idx_dip_fg + 273);

    auto tr_y_yzz_xxyz = pbuffer.data(idx_dip_fg + 274);

    auto tr_y_yzz_xxzz = pbuffer.data(idx_dip_fg + 275);

    auto tr_y_yzz_xyyy = pbuffer.data(idx_dip_fg + 276);

    auto tr_y_yzz_xyyz = pbuffer.data(idx_dip_fg + 277);

    auto tr_y_yzz_xyzz = pbuffer.data(idx_dip_fg + 278);

    auto tr_y_yzz_xzzz = pbuffer.data(idx_dip_fg + 279);

    auto tr_y_yzz_yyyy = pbuffer.data(idx_dip_fg + 280);

    auto tr_y_yzz_yyyz = pbuffer.data(idx_dip_fg + 281);

    auto tr_y_yzz_yyzz = pbuffer.data(idx_dip_fg + 282);

    auto tr_y_yzz_yzzz = pbuffer.data(idx_dip_fg + 283);

    auto tr_y_yzz_zzzz = pbuffer.data(idx_dip_fg + 284);

    auto tr_y_zzz_xxxx = pbuffer.data(idx_dip_fg + 285);

    auto tr_y_zzz_xxxy = pbuffer.data(idx_dip_fg + 286);

    auto tr_y_zzz_xxxz = pbuffer.data(idx_dip_fg + 287);

    auto tr_y_zzz_xxyy = pbuffer.data(idx_dip_fg + 288);

    auto tr_y_zzz_xxyz = pbuffer.data(idx_dip_fg + 289);

    auto tr_y_zzz_xxzz = pbuffer.data(idx_dip_fg + 290);

    auto tr_y_zzz_xyyy = pbuffer.data(idx_dip_fg + 291);

    auto tr_y_zzz_xyyz = pbuffer.data(idx_dip_fg + 292);

    auto tr_y_zzz_xyzz = pbuffer.data(idx_dip_fg + 293);

    auto tr_y_zzz_xzzz = pbuffer.data(idx_dip_fg + 294);

    auto tr_y_zzz_yyyy = pbuffer.data(idx_dip_fg + 295);

    auto tr_y_zzz_yyyz = pbuffer.data(idx_dip_fg + 296);

    auto tr_y_zzz_yyzz = pbuffer.data(idx_dip_fg + 297);

    auto tr_y_zzz_yzzz = pbuffer.data(idx_dip_fg + 298);

    auto tr_y_zzz_zzzz = pbuffer.data(idx_dip_fg + 299);

    auto tr_z_xxx_xxxx = pbuffer.data(idx_dip_fg + 300);

    auto tr_z_xxx_xxxy = pbuffer.data(idx_dip_fg + 301);

    auto tr_z_xxx_xxxz = pbuffer.data(idx_dip_fg + 302);

    auto tr_z_xxx_xxyy = pbuffer.data(idx_dip_fg + 303);

    auto tr_z_xxx_xxyz = pbuffer.data(idx_dip_fg + 304);

    auto tr_z_xxx_xxzz = pbuffer.data(idx_dip_fg + 305);

    auto tr_z_xxx_xyyy = pbuffer.data(idx_dip_fg + 306);

    auto tr_z_xxx_xyyz = pbuffer.data(idx_dip_fg + 307);

    auto tr_z_xxx_xyzz = pbuffer.data(idx_dip_fg + 308);

    auto tr_z_xxx_xzzz = pbuffer.data(idx_dip_fg + 309);

    auto tr_z_xxx_yyyy = pbuffer.data(idx_dip_fg + 310);

    auto tr_z_xxx_yyyz = pbuffer.data(idx_dip_fg + 311);

    auto tr_z_xxx_yyzz = pbuffer.data(idx_dip_fg + 312);

    auto tr_z_xxx_yzzz = pbuffer.data(idx_dip_fg + 313);

    auto tr_z_xxx_zzzz = pbuffer.data(idx_dip_fg + 314);

    auto tr_z_xxy_xxxx = pbuffer.data(idx_dip_fg + 315);

    auto tr_z_xxy_xxxz = pbuffer.data(idx_dip_fg + 317);

    auto tr_z_xxy_xxzz = pbuffer.data(idx_dip_fg + 320);

    auto tr_z_xxy_xzzz = pbuffer.data(idx_dip_fg + 324);

    auto tr_z_xxy_yyyy = pbuffer.data(idx_dip_fg + 325);

    auto tr_z_xxy_yyyz = pbuffer.data(idx_dip_fg + 326);

    auto tr_z_xxy_yyzz = pbuffer.data(idx_dip_fg + 327);

    auto tr_z_xxy_yzzz = pbuffer.data(idx_dip_fg + 328);

    auto tr_z_xxz_xxxx = pbuffer.data(idx_dip_fg + 330);

    auto tr_z_xxz_xxxz = pbuffer.data(idx_dip_fg + 332);

    auto tr_z_xxz_xxyz = pbuffer.data(idx_dip_fg + 334);

    auto tr_z_xxz_xxzz = pbuffer.data(idx_dip_fg + 335);

    auto tr_z_xxz_xyyz = pbuffer.data(idx_dip_fg + 337);

    auto tr_z_xxz_xyzz = pbuffer.data(idx_dip_fg + 338);

    auto tr_z_xxz_xzzz = pbuffer.data(idx_dip_fg + 339);

    auto tr_z_xxz_yyyy = pbuffer.data(idx_dip_fg + 340);

    auto tr_z_xxz_yyyz = pbuffer.data(idx_dip_fg + 341);

    auto tr_z_xxz_yyzz = pbuffer.data(idx_dip_fg + 342);

    auto tr_z_xxz_yzzz = pbuffer.data(idx_dip_fg + 343);

    auto tr_z_xxz_zzzz = pbuffer.data(idx_dip_fg + 344);

    auto tr_z_xyy_xxxy = pbuffer.data(idx_dip_fg + 346);

    auto tr_z_xyy_xxyy = pbuffer.data(idx_dip_fg + 348);

    auto tr_z_xyy_xxyz = pbuffer.data(idx_dip_fg + 349);

    auto tr_z_xyy_xyyy = pbuffer.data(idx_dip_fg + 351);

    auto tr_z_xyy_xyyz = pbuffer.data(idx_dip_fg + 352);

    auto tr_z_xyy_xyzz = pbuffer.data(idx_dip_fg + 353);

    auto tr_z_xyy_yyyy = pbuffer.data(idx_dip_fg + 355);

    auto tr_z_xyy_yyyz = pbuffer.data(idx_dip_fg + 356);

    auto tr_z_xyy_yyzz = pbuffer.data(idx_dip_fg + 357);

    auto tr_z_xyy_yzzz = pbuffer.data(idx_dip_fg + 358);

    auto tr_z_xyy_zzzz = pbuffer.data(idx_dip_fg + 359);

    auto tr_z_xyz_yyyy = pbuffer.data(idx_dip_fg + 370);

    auto tr_z_xyz_yyyz = pbuffer.data(idx_dip_fg + 371);

    auto tr_z_xyz_yyzz = pbuffer.data(idx_dip_fg + 372);

    auto tr_z_xyz_yzzz = pbuffer.data(idx_dip_fg + 373);

    auto tr_z_xzz_xxxx = pbuffer.data(idx_dip_fg + 375);

    auto tr_z_xzz_xxxy = pbuffer.data(idx_dip_fg + 376);

    auto tr_z_xzz_xxxz = pbuffer.data(idx_dip_fg + 377);

    auto tr_z_xzz_xxyy = pbuffer.data(idx_dip_fg + 378);

    auto tr_z_xzz_xxyz = pbuffer.data(idx_dip_fg + 379);

    auto tr_z_xzz_xxzz = pbuffer.data(idx_dip_fg + 380);

    auto tr_z_xzz_xyyy = pbuffer.data(idx_dip_fg + 381);

    auto tr_z_xzz_xyyz = pbuffer.data(idx_dip_fg + 382);

    auto tr_z_xzz_xyzz = pbuffer.data(idx_dip_fg + 383);

    auto tr_z_xzz_xzzz = pbuffer.data(idx_dip_fg + 384);

    auto tr_z_xzz_yyyy = pbuffer.data(idx_dip_fg + 385);

    auto tr_z_xzz_yyyz = pbuffer.data(idx_dip_fg + 386);

    auto tr_z_xzz_yyzz = pbuffer.data(idx_dip_fg + 387);

    auto tr_z_xzz_yzzz = pbuffer.data(idx_dip_fg + 388);

    auto tr_z_xzz_zzzz = pbuffer.data(idx_dip_fg + 389);

    auto tr_z_yyy_xxxx = pbuffer.data(idx_dip_fg + 390);

    auto tr_z_yyy_xxxy = pbuffer.data(idx_dip_fg + 391);

    auto tr_z_yyy_xxxz = pbuffer.data(idx_dip_fg + 392);

    auto tr_z_yyy_xxyy = pbuffer.data(idx_dip_fg + 393);

    auto tr_z_yyy_xxyz = pbuffer.data(idx_dip_fg + 394);

    auto tr_z_yyy_xxzz = pbuffer.data(idx_dip_fg + 395);

    auto tr_z_yyy_xyyy = pbuffer.data(idx_dip_fg + 396);

    auto tr_z_yyy_xyyz = pbuffer.data(idx_dip_fg + 397);

    auto tr_z_yyy_xyzz = pbuffer.data(idx_dip_fg + 398);

    auto tr_z_yyy_xzzz = pbuffer.data(idx_dip_fg + 399);

    auto tr_z_yyy_yyyy = pbuffer.data(idx_dip_fg + 400);

    auto tr_z_yyy_yyyz = pbuffer.data(idx_dip_fg + 401);

    auto tr_z_yyy_yyzz = pbuffer.data(idx_dip_fg + 402);

    auto tr_z_yyy_yzzz = pbuffer.data(idx_dip_fg + 403);

    auto tr_z_yyy_zzzz = pbuffer.data(idx_dip_fg + 404);

    auto tr_z_yyz_xxxx = pbuffer.data(idx_dip_fg + 405);

    auto tr_z_yyz_xxxz = pbuffer.data(idx_dip_fg + 407);

    auto tr_z_yyz_xxyz = pbuffer.data(idx_dip_fg + 409);

    auto tr_z_yyz_xxzz = pbuffer.data(idx_dip_fg + 410);

    auto tr_z_yyz_xyyz = pbuffer.data(idx_dip_fg + 412);

    auto tr_z_yyz_xyzz = pbuffer.data(idx_dip_fg + 413);

    auto tr_z_yyz_xzzz = pbuffer.data(idx_dip_fg + 414);

    auto tr_z_yyz_yyyy = pbuffer.data(idx_dip_fg + 415);

    auto tr_z_yyz_yyyz = pbuffer.data(idx_dip_fg + 416);

    auto tr_z_yyz_yyzz = pbuffer.data(idx_dip_fg + 417);

    auto tr_z_yyz_yzzz = pbuffer.data(idx_dip_fg + 418);

    auto tr_z_yyz_zzzz = pbuffer.data(idx_dip_fg + 419);

    auto tr_z_yzz_xxxx = pbuffer.data(idx_dip_fg + 420);

    auto tr_z_yzz_xxxy = pbuffer.data(idx_dip_fg + 421);

    auto tr_z_yzz_xxxz = pbuffer.data(idx_dip_fg + 422);

    auto tr_z_yzz_xxyy = pbuffer.data(idx_dip_fg + 423);

    auto tr_z_yzz_xxyz = pbuffer.data(idx_dip_fg + 424);

    auto tr_z_yzz_xxzz = pbuffer.data(idx_dip_fg + 425);

    auto tr_z_yzz_xyyy = pbuffer.data(idx_dip_fg + 426);

    auto tr_z_yzz_xyyz = pbuffer.data(idx_dip_fg + 427);

    auto tr_z_yzz_xyzz = pbuffer.data(idx_dip_fg + 428);

    auto tr_z_yzz_xzzz = pbuffer.data(idx_dip_fg + 429);

    auto tr_z_yzz_yyyy = pbuffer.data(idx_dip_fg + 430);

    auto tr_z_yzz_yyyz = pbuffer.data(idx_dip_fg + 431);

    auto tr_z_yzz_yyzz = pbuffer.data(idx_dip_fg + 432);

    auto tr_z_yzz_yzzz = pbuffer.data(idx_dip_fg + 433);

    auto tr_z_yzz_zzzz = pbuffer.data(idx_dip_fg + 434);

    auto tr_z_zzz_xxxx = pbuffer.data(idx_dip_fg + 435);

    auto tr_z_zzz_xxxy = pbuffer.data(idx_dip_fg + 436);

    auto tr_z_zzz_xxxz = pbuffer.data(idx_dip_fg + 437);

    auto tr_z_zzz_xxyy = pbuffer.data(idx_dip_fg + 438);

    auto tr_z_zzz_xxyz = pbuffer.data(idx_dip_fg + 439);

    auto tr_z_zzz_xxzz = pbuffer.data(idx_dip_fg + 440);

    auto tr_z_zzz_xyyy = pbuffer.data(idx_dip_fg + 441);

    auto tr_z_zzz_xyyz = pbuffer.data(idx_dip_fg + 442);

    auto tr_z_zzz_xyzz = pbuffer.data(idx_dip_fg + 443);

    auto tr_z_zzz_xzzz = pbuffer.data(idx_dip_fg + 444);

    auto tr_z_zzz_yyyy = pbuffer.data(idx_dip_fg + 445);

    auto tr_z_zzz_yyyz = pbuffer.data(idx_dip_fg + 446);

    auto tr_z_zzz_yyzz = pbuffer.data(idx_dip_fg + 447);

    auto tr_z_zzz_yzzz = pbuffer.data(idx_dip_fg + 448);

    auto tr_z_zzz_zzzz = pbuffer.data(idx_dip_fg + 449);

    // Set up components of auxiliary buffer : GF

    auto tr_x_xxxx_xxx = pbuffer.data(idx_dip_gf);

    auto tr_x_xxxx_xxy = pbuffer.data(idx_dip_gf + 1);

    auto tr_x_xxxx_xxz = pbuffer.data(idx_dip_gf + 2);

    auto tr_x_xxxx_xyy = pbuffer.data(idx_dip_gf + 3);

    auto tr_x_xxxx_xyz = pbuffer.data(idx_dip_gf + 4);

    auto tr_x_xxxx_xzz = pbuffer.data(idx_dip_gf + 5);

    auto tr_x_xxxx_yyy = pbuffer.data(idx_dip_gf + 6);

    auto tr_x_xxxx_yyz = pbuffer.data(idx_dip_gf + 7);

    auto tr_x_xxxx_yzz = pbuffer.data(idx_dip_gf + 8);

    auto tr_x_xxxx_zzz = pbuffer.data(idx_dip_gf + 9);

    auto tr_x_xxxy_xxx = pbuffer.data(idx_dip_gf + 10);

    auto tr_x_xxxy_xxy = pbuffer.data(idx_dip_gf + 11);

    auto tr_x_xxxy_xxz = pbuffer.data(idx_dip_gf + 12);

    auto tr_x_xxxy_xyy = pbuffer.data(idx_dip_gf + 13);

    auto tr_x_xxxy_xyz = pbuffer.data(idx_dip_gf + 14);

    auto tr_x_xxxy_xzz = pbuffer.data(idx_dip_gf + 15);

    auto tr_x_xxxz_xxx = pbuffer.data(idx_dip_gf + 20);

    auto tr_x_xxxz_xxy = pbuffer.data(idx_dip_gf + 21);

    auto tr_x_xxxz_xxz = pbuffer.data(idx_dip_gf + 22);

    auto tr_x_xxxz_xyy = pbuffer.data(idx_dip_gf + 23);

    auto tr_x_xxxz_xyz = pbuffer.data(idx_dip_gf + 24);

    auto tr_x_xxxz_xzz = pbuffer.data(idx_dip_gf + 25);

    auto tr_x_xxxz_yyz = pbuffer.data(idx_dip_gf + 27);

    auto tr_x_xxxz_yzz = pbuffer.data(idx_dip_gf + 28);

    auto tr_x_xxxz_zzz = pbuffer.data(idx_dip_gf + 29);

    auto tr_x_xxyy_xxx = pbuffer.data(idx_dip_gf + 30);

    auto tr_x_xxyy_xxy = pbuffer.data(idx_dip_gf + 31);

    auto tr_x_xxyy_xxz = pbuffer.data(idx_dip_gf + 32);

    auto tr_x_xxyy_xyy = pbuffer.data(idx_dip_gf + 33);

    auto tr_x_xxyy_xyz = pbuffer.data(idx_dip_gf + 34);

    auto tr_x_xxyy_xzz = pbuffer.data(idx_dip_gf + 35);

    auto tr_x_xxyy_yyy = pbuffer.data(idx_dip_gf + 36);

    auto tr_x_xxyy_yyz = pbuffer.data(idx_dip_gf + 37);

    auto tr_x_xxyy_yzz = pbuffer.data(idx_dip_gf + 38);

    auto tr_x_xxzz_xxx = pbuffer.data(idx_dip_gf + 50);

    auto tr_x_xxzz_xxy = pbuffer.data(idx_dip_gf + 51);

    auto tr_x_xxzz_xxz = pbuffer.data(idx_dip_gf + 52);

    auto tr_x_xxzz_xyy = pbuffer.data(idx_dip_gf + 53);

    auto tr_x_xxzz_xyz = pbuffer.data(idx_dip_gf + 54);

    auto tr_x_xxzz_xzz = pbuffer.data(idx_dip_gf + 55);

    auto tr_x_xxzz_yyy = pbuffer.data(idx_dip_gf + 56);

    auto tr_x_xxzz_yyz = pbuffer.data(idx_dip_gf + 57);

    auto tr_x_xxzz_yzz = pbuffer.data(idx_dip_gf + 58);

    auto tr_x_xxzz_zzz = pbuffer.data(idx_dip_gf + 59);

    auto tr_x_xyyy_xxy = pbuffer.data(idx_dip_gf + 61);

    auto tr_x_xyyy_xyy = pbuffer.data(idx_dip_gf + 63);

    auto tr_x_xyyy_xyz = pbuffer.data(idx_dip_gf + 64);

    auto tr_x_xzzz_xxx = pbuffer.data(idx_dip_gf + 90);

    auto tr_x_xzzz_xxy = pbuffer.data(idx_dip_gf + 91);

    auto tr_x_xzzz_xxz = pbuffer.data(idx_dip_gf + 92);

    auto tr_x_xzzz_xyy = pbuffer.data(idx_dip_gf + 93);

    auto tr_x_xzzz_xyz = pbuffer.data(idx_dip_gf + 94);

    auto tr_x_xzzz_xzz = pbuffer.data(idx_dip_gf + 95);

    auto tr_x_yyyy_xxx = pbuffer.data(idx_dip_gf + 100);

    auto tr_x_yyyy_xxy = pbuffer.data(idx_dip_gf + 101);

    auto tr_x_yyyy_xxz = pbuffer.data(idx_dip_gf + 102);

    auto tr_x_yyyy_xyy = pbuffer.data(idx_dip_gf + 103);

    auto tr_x_yyyy_xyz = pbuffer.data(idx_dip_gf + 104);

    auto tr_x_yyyy_xzz = pbuffer.data(idx_dip_gf + 105);

    auto tr_x_yyyy_yyy = pbuffer.data(idx_dip_gf + 106);

    auto tr_x_yyyy_yyz = pbuffer.data(idx_dip_gf + 107);

    auto tr_x_yyyy_yzz = pbuffer.data(idx_dip_gf + 108);

    auto tr_x_yyyy_zzz = pbuffer.data(idx_dip_gf + 109);

    auto tr_x_yyzz_xxz = pbuffer.data(idx_dip_gf + 122);

    auto tr_x_yyzz_xyz = pbuffer.data(idx_dip_gf + 124);

    auto tr_x_yyzz_xzz = pbuffer.data(idx_dip_gf + 125);

    auto tr_x_yyzz_yyz = pbuffer.data(idx_dip_gf + 127);

    auto tr_x_yyzz_yzz = pbuffer.data(idx_dip_gf + 128);

    auto tr_x_yyzz_zzz = pbuffer.data(idx_dip_gf + 129);

    auto tr_x_yzzz_xxz = pbuffer.data(idx_dip_gf + 132);

    auto tr_x_yzzz_xyz = pbuffer.data(idx_dip_gf + 134);

    auto tr_x_yzzz_xzz = pbuffer.data(idx_dip_gf + 135);

    auto tr_x_yzzz_yyz = pbuffer.data(idx_dip_gf + 137);

    auto tr_x_yzzz_yzz = pbuffer.data(idx_dip_gf + 138);

    auto tr_x_yzzz_zzz = pbuffer.data(idx_dip_gf + 139);

    auto tr_x_zzzz_xxx = pbuffer.data(idx_dip_gf + 140);

    auto tr_x_zzzz_xxy = pbuffer.data(idx_dip_gf + 141);

    auto tr_x_zzzz_xxz = pbuffer.data(idx_dip_gf + 142);

    auto tr_x_zzzz_xyy = pbuffer.data(idx_dip_gf + 143);

    auto tr_x_zzzz_xyz = pbuffer.data(idx_dip_gf + 144);

    auto tr_x_zzzz_xzz = pbuffer.data(idx_dip_gf + 145);

    auto tr_x_zzzz_yyy = pbuffer.data(idx_dip_gf + 146);

    auto tr_x_zzzz_yyz = pbuffer.data(idx_dip_gf + 147);

    auto tr_x_zzzz_yzz = pbuffer.data(idx_dip_gf + 148);

    auto tr_x_zzzz_zzz = pbuffer.data(idx_dip_gf + 149);

    auto tr_y_xxxx_xxx = pbuffer.data(idx_dip_gf + 150);

    auto tr_y_xxxx_xxy = pbuffer.data(idx_dip_gf + 151);

    auto tr_y_xxxx_xxz = pbuffer.data(idx_dip_gf + 152);

    auto tr_y_xxxx_xyy = pbuffer.data(idx_dip_gf + 153);

    auto tr_y_xxxx_xyz = pbuffer.data(idx_dip_gf + 154);

    auto tr_y_xxxx_xzz = pbuffer.data(idx_dip_gf + 155);

    auto tr_y_xxxx_yyy = pbuffer.data(idx_dip_gf + 156);

    auto tr_y_xxxx_yyz = pbuffer.data(idx_dip_gf + 157);

    auto tr_y_xxxx_yzz = pbuffer.data(idx_dip_gf + 158);

    auto tr_y_xxxx_zzz = pbuffer.data(idx_dip_gf + 159);

    auto tr_y_xxxy_xxy = pbuffer.data(idx_dip_gf + 161);

    auto tr_y_xxxy_xyy = pbuffer.data(idx_dip_gf + 163);

    auto tr_y_xxxy_xyz = pbuffer.data(idx_dip_gf + 164);

    auto tr_y_xxxy_yyy = pbuffer.data(idx_dip_gf + 166);

    auto tr_y_xxxy_yyz = pbuffer.data(idx_dip_gf + 167);

    auto tr_y_xxxy_yzz = pbuffer.data(idx_dip_gf + 168);

    auto tr_y_xxyy_xxx = pbuffer.data(idx_dip_gf + 180);

    auto tr_y_xxyy_xxy = pbuffer.data(idx_dip_gf + 181);

    auto tr_y_xxyy_xxz = pbuffer.data(idx_dip_gf + 182);

    auto tr_y_xxyy_xyy = pbuffer.data(idx_dip_gf + 183);

    auto tr_y_xxyy_xyz = pbuffer.data(idx_dip_gf + 184);

    auto tr_y_xxyy_xzz = pbuffer.data(idx_dip_gf + 185);

    auto tr_y_xxyy_yyy = pbuffer.data(idx_dip_gf + 186);

    auto tr_y_xxyy_yyz = pbuffer.data(idx_dip_gf + 187);

    auto tr_y_xxyy_yzz = pbuffer.data(idx_dip_gf + 188);

    auto tr_y_xxyy_zzz = pbuffer.data(idx_dip_gf + 189);

    auto tr_y_xxzz_xxz = pbuffer.data(idx_dip_gf + 202);

    auto tr_y_xxzz_xyz = pbuffer.data(idx_dip_gf + 204);

    auto tr_y_xxzz_xzz = pbuffer.data(idx_dip_gf + 205);

    auto tr_y_xxzz_yyz = pbuffer.data(idx_dip_gf + 207);

    auto tr_y_xxzz_yzz = pbuffer.data(idx_dip_gf + 208);

    auto tr_y_xxzz_zzz = pbuffer.data(idx_dip_gf + 209);

    auto tr_y_xyyy_xxx = pbuffer.data(idx_dip_gf + 210);

    auto tr_y_xyyy_xxy = pbuffer.data(idx_dip_gf + 211);

    auto tr_y_xyyy_xxz = pbuffer.data(idx_dip_gf + 212);

    auto tr_y_xyyy_xyy = pbuffer.data(idx_dip_gf + 213);

    auto tr_y_xyyy_xyz = pbuffer.data(idx_dip_gf + 214);

    auto tr_y_xyyy_xzz = pbuffer.data(idx_dip_gf + 215);

    auto tr_y_xyyy_yyy = pbuffer.data(idx_dip_gf + 216);

    auto tr_y_xyyy_yyz = pbuffer.data(idx_dip_gf + 217);

    auto tr_y_xyyy_yzz = pbuffer.data(idx_dip_gf + 218);

    auto tr_y_xyyy_zzz = pbuffer.data(idx_dip_gf + 219);

    auto tr_y_xyzz_xyz = pbuffer.data(idx_dip_gf + 234);

    auto tr_y_xyzz_yyz = pbuffer.data(idx_dip_gf + 237);

    auto tr_y_xyzz_yzz = pbuffer.data(idx_dip_gf + 238);

    auto tr_y_xzzz_xxz = pbuffer.data(idx_dip_gf + 242);

    auto tr_y_xzzz_xyz = pbuffer.data(idx_dip_gf + 244);

    auto tr_y_xzzz_xzz = pbuffer.data(idx_dip_gf + 245);

    auto tr_y_xzzz_yyz = pbuffer.data(idx_dip_gf + 247);

    auto tr_y_xzzz_yzz = pbuffer.data(idx_dip_gf + 248);

    auto tr_y_xzzz_zzz = pbuffer.data(idx_dip_gf + 249);

    auto tr_y_yyyy_xxx = pbuffer.data(idx_dip_gf + 250);

    auto tr_y_yyyy_xxy = pbuffer.data(idx_dip_gf + 251);

    auto tr_y_yyyy_xxz = pbuffer.data(idx_dip_gf + 252);

    auto tr_y_yyyy_xyy = pbuffer.data(idx_dip_gf + 253);

    auto tr_y_yyyy_xyz = pbuffer.data(idx_dip_gf + 254);

    auto tr_y_yyyy_xzz = pbuffer.data(idx_dip_gf + 255);

    auto tr_y_yyyy_yyy = pbuffer.data(idx_dip_gf + 256);

    auto tr_y_yyyy_yyz = pbuffer.data(idx_dip_gf + 257);

    auto tr_y_yyyy_yzz = pbuffer.data(idx_dip_gf + 258);

    auto tr_y_yyyy_zzz = pbuffer.data(idx_dip_gf + 259);

    auto tr_y_yyyz_xxy = pbuffer.data(idx_dip_gf + 261);

    auto tr_y_yyyz_xxz = pbuffer.data(idx_dip_gf + 262);

    auto tr_y_yyyz_xyy = pbuffer.data(idx_dip_gf + 263);

    auto tr_y_yyyz_xyz = pbuffer.data(idx_dip_gf + 264);

    auto tr_y_yyyz_xzz = pbuffer.data(idx_dip_gf + 265);

    auto tr_y_yyyz_yyy = pbuffer.data(idx_dip_gf + 266);

    auto tr_y_yyyz_yyz = pbuffer.data(idx_dip_gf + 267);

    auto tr_y_yyyz_yzz = pbuffer.data(idx_dip_gf + 268);

    auto tr_y_yyyz_zzz = pbuffer.data(idx_dip_gf + 269);

    auto tr_y_yyzz_xxx = pbuffer.data(idx_dip_gf + 270);

    auto tr_y_yyzz_xxy = pbuffer.data(idx_dip_gf + 271);

    auto tr_y_yyzz_xxz = pbuffer.data(idx_dip_gf + 272);

    auto tr_y_yyzz_xyy = pbuffer.data(idx_dip_gf + 273);

    auto tr_y_yyzz_xyz = pbuffer.data(idx_dip_gf + 274);

    auto tr_y_yyzz_xzz = pbuffer.data(idx_dip_gf + 275);

    auto tr_y_yyzz_yyy = pbuffer.data(idx_dip_gf + 276);

    auto tr_y_yyzz_yyz = pbuffer.data(idx_dip_gf + 277);

    auto tr_y_yyzz_yzz = pbuffer.data(idx_dip_gf + 278);

    auto tr_y_yyzz_zzz = pbuffer.data(idx_dip_gf + 279);

    auto tr_y_yzzz_xxx = pbuffer.data(idx_dip_gf + 280);

    auto tr_y_yzzz_xxy = pbuffer.data(idx_dip_gf + 281);

    auto tr_y_yzzz_xxz = pbuffer.data(idx_dip_gf + 282);

    auto tr_y_yzzz_xyy = pbuffer.data(idx_dip_gf + 283);

    auto tr_y_yzzz_xyz = pbuffer.data(idx_dip_gf + 284);

    auto tr_y_yzzz_xzz = pbuffer.data(idx_dip_gf + 285);

    auto tr_y_yzzz_yyy = pbuffer.data(idx_dip_gf + 286);

    auto tr_y_yzzz_yyz = pbuffer.data(idx_dip_gf + 287);

    auto tr_y_yzzz_yzz = pbuffer.data(idx_dip_gf + 288);

    auto tr_y_yzzz_zzz = pbuffer.data(idx_dip_gf + 289);

    auto tr_y_zzzz_xxx = pbuffer.data(idx_dip_gf + 290);

    auto tr_y_zzzz_xxy = pbuffer.data(idx_dip_gf + 291);

    auto tr_y_zzzz_xxz = pbuffer.data(idx_dip_gf + 292);

    auto tr_y_zzzz_xyy = pbuffer.data(idx_dip_gf + 293);

    auto tr_y_zzzz_xyz = pbuffer.data(idx_dip_gf + 294);

    auto tr_y_zzzz_xzz = pbuffer.data(idx_dip_gf + 295);

    auto tr_y_zzzz_yyy = pbuffer.data(idx_dip_gf + 296);

    auto tr_y_zzzz_yyz = pbuffer.data(idx_dip_gf + 297);

    auto tr_y_zzzz_yzz = pbuffer.data(idx_dip_gf + 298);

    auto tr_y_zzzz_zzz = pbuffer.data(idx_dip_gf + 299);

    auto tr_z_xxxx_xxx = pbuffer.data(idx_dip_gf + 300);

    auto tr_z_xxxx_xxy = pbuffer.data(idx_dip_gf + 301);

    auto tr_z_xxxx_xxz = pbuffer.data(idx_dip_gf + 302);

    auto tr_z_xxxx_xyy = pbuffer.data(idx_dip_gf + 303);

    auto tr_z_xxxx_xyz = pbuffer.data(idx_dip_gf + 304);

    auto tr_z_xxxx_xzz = pbuffer.data(idx_dip_gf + 305);

    auto tr_z_xxxx_yyy = pbuffer.data(idx_dip_gf + 306);

    auto tr_z_xxxx_yyz = pbuffer.data(idx_dip_gf + 307);

    auto tr_z_xxxx_yzz = pbuffer.data(idx_dip_gf + 308);

    auto tr_z_xxxx_zzz = pbuffer.data(idx_dip_gf + 309);

    auto tr_z_xxxz_xxx = pbuffer.data(idx_dip_gf + 320);

    auto tr_z_xxxz_xxy = pbuffer.data(idx_dip_gf + 321);

    auto tr_z_xxxz_xxz = pbuffer.data(idx_dip_gf + 322);

    auto tr_z_xxxz_xyy = pbuffer.data(idx_dip_gf + 323);

    auto tr_z_xxxz_xyz = pbuffer.data(idx_dip_gf + 324);

    auto tr_z_xxxz_xzz = pbuffer.data(idx_dip_gf + 325);

    auto tr_z_xxxz_yyz = pbuffer.data(idx_dip_gf + 327);

    auto tr_z_xxxz_yzz = pbuffer.data(idx_dip_gf + 328);

    auto tr_z_xxxz_zzz = pbuffer.data(idx_dip_gf + 329);

    auto tr_z_xxyy_xxy = pbuffer.data(idx_dip_gf + 331);

    auto tr_z_xxyy_xyy = pbuffer.data(idx_dip_gf + 333);

    auto tr_z_xxyy_xyz = pbuffer.data(idx_dip_gf + 334);

    auto tr_z_xxyy_yyy = pbuffer.data(idx_dip_gf + 336);

    auto tr_z_xxyy_yyz = pbuffer.data(idx_dip_gf + 337);

    auto tr_z_xxyy_yzz = pbuffer.data(idx_dip_gf + 338);

    auto tr_z_xxzz_xxx = pbuffer.data(idx_dip_gf + 350);

    auto tr_z_xxzz_xxy = pbuffer.data(idx_dip_gf + 351);

    auto tr_z_xxzz_xxz = pbuffer.data(idx_dip_gf + 352);

    auto tr_z_xxzz_xyy = pbuffer.data(idx_dip_gf + 353);

    auto tr_z_xxzz_xyz = pbuffer.data(idx_dip_gf + 354);

    auto tr_z_xxzz_xzz = pbuffer.data(idx_dip_gf + 355);

    auto tr_z_xxzz_yyy = pbuffer.data(idx_dip_gf + 356);

    auto tr_z_xxzz_yyz = pbuffer.data(idx_dip_gf + 357);

    auto tr_z_xxzz_yzz = pbuffer.data(idx_dip_gf + 358);

    auto tr_z_xxzz_zzz = pbuffer.data(idx_dip_gf + 359);

    auto tr_z_xyyy_xxy = pbuffer.data(idx_dip_gf + 361);

    auto tr_z_xyyy_xyy = pbuffer.data(idx_dip_gf + 363);

    auto tr_z_xyyy_xyz = pbuffer.data(idx_dip_gf + 364);

    auto tr_z_xyyy_yyy = pbuffer.data(idx_dip_gf + 366);

    auto tr_z_xyyy_yyz = pbuffer.data(idx_dip_gf + 367);

    auto tr_z_xyyy_yzz = pbuffer.data(idx_dip_gf + 368);

    auto tr_z_xyyz_xyz = pbuffer.data(idx_dip_gf + 374);

    auto tr_z_xyyz_yyz = pbuffer.data(idx_dip_gf + 377);

    auto tr_z_xyyz_yzz = pbuffer.data(idx_dip_gf + 378);

    auto tr_z_xzzz_xxx = pbuffer.data(idx_dip_gf + 390);

    auto tr_z_xzzz_xxy = pbuffer.data(idx_dip_gf + 391);

    auto tr_z_xzzz_xxz = pbuffer.data(idx_dip_gf + 392);

    auto tr_z_xzzz_xyy = pbuffer.data(idx_dip_gf + 393);

    auto tr_z_xzzz_xyz = pbuffer.data(idx_dip_gf + 394);

    auto tr_z_xzzz_xzz = pbuffer.data(idx_dip_gf + 395);

    auto tr_z_xzzz_yyy = pbuffer.data(idx_dip_gf + 396);

    auto tr_z_xzzz_yyz = pbuffer.data(idx_dip_gf + 397);

    auto tr_z_xzzz_yzz = pbuffer.data(idx_dip_gf + 398);

    auto tr_z_xzzz_zzz = pbuffer.data(idx_dip_gf + 399);

    auto tr_z_yyyy_xxx = pbuffer.data(idx_dip_gf + 400);

    auto tr_z_yyyy_xxy = pbuffer.data(idx_dip_gf + 401);

    auto tr_z_yyyy_xxz = pbuffer.data(idx_dip_gf + 402);

    auto tr_z_yyyy_xyy = pbuffer.data(idx_dip_gf + 403);

    auto tr_z_yyyy_xyz = pbuffer.data(idx_dip_gf + 404);

    auto tr_z_yyyy_xzz = pbuffer.data(idx_dip_gf + 405);

    auto tr_z_yyyy_yyy = pbuffer.data(idx_dip_gf + 406);

    auto tr_z_yyyy_yyz = pbuffer.data(idx_dip_gf + 407);

    auto tr_z_yyyy_yzz = pbuffer.data(idx_dip_gf + 408);

    auto tr_z_yyyy_zzz = pbuffer.data(idx_dip_gf + 409);

    auto tr_z_yyyz_xxx = pbuffer.data(idx_dip_gf + 410);

    auto tr_z_yyyz_xxy = pbuffer.data(idx_dip_gf + 411);

    auto tr_z_yyyz_xxz = pbuffer.data(idx_dip_gf + 412);

    auto tr_z_yyyz_xyy = pbuffer.data(idx_dip_gf + 413);

    auto tr_z_yyyz_xyz = pbuffer.data(idx_dip_gf + 414);

    auto tr_z_yyyz_xzz = pbuffer.data(idx_dip_gf + 415);

    auto tr_z_yyyz_yyy = pbuffer.data(idx_dip_gf + 416);

    auto tr_z_yyyz_yyz = pbuffer.data(idx_dip_gf + 417);

    auto tr_z_yyyz_yzz = pbuffer.data(idx_dip_gf + 418);

    auto tr_z_yyyz_zzz = pbuffer.data(idx_dip_gf + 419);

    auto tr_z_yyzz_xxx = pbuffer.data(idx_dip_gf + 420);

    auto tr_z_yyzz_xxy = pbuffer.data(idx_dip_gf + 421);

    auto tr_z_yyzz_xxz = pbuffer.data(idx_dip_gf + 422);

    auto tr_z_yyzz_xyy = pbuffer.data(idx_dip_gf + 423);

    auto tr_z_yyzz_xyz = pbuffer.data(idx_dip_gf + 424);

    auto tr_z_yyzz_xzz = pbuffer.data(idx_dip_gf + 425);

    auto tr_z_yyzz_yyy = pbuffer.data(idx_dip_gf + 426);

    auto tr_z_yyzz_yyz = pbuffer.data(idx_dip_gf + 427);

    auto tr_z_yyzz_yzz = pbuffer.data(idx_dip_gf + 428);

    auto tr_z_yyzz_zzz = pbuffer.data(idx_dip_gf + 429);

    auto tr_z_yzzz_xxx = pbuffer.data(idx_dip_gf + 430);

    auto tr_z_yzzz_xxy = pbuffer.data(idx_dip_gf + 431);

    auto tr_z_yzzz_xxz = pbuffer.data(idx_dip_gf + 432);

    auto tr_z_yzzz_xyy = pbuffer.data(idx_dip_gf + 433);

    auto tr_z_yzzz_xyz = pbuffer.data(idx_dip_gf + 434);

    auto tr_z_yzzz_xzz = pbuffer.data(idx_dip_gf + 435);

    auto tr_z_yzzz_yyy = pbuffer.data(idx_dip_gf + 436);

    auto tr_z_yzzz_yyz = pbuffer.data(idx_dip_gf + 437);

    auto tr_z_yzzz_yzz = pbuffer.data(idx_dip_gf + 438);

    auto tr_z_yzzz_zzz = pbuffer.data(idx_dip_gf + 439);

    auto tr_z_zzzz_xxx = pbuffer.data(idx_dip_gf + 440);

    auto tr_z_zzzz_xxy = pbuffer.data(idx_dip_gf + 441);

    auto tr_z_zzzz_xxz = pbuffer.data(idx_dip_gf + 442);

    auto tr_z_zzzz_xyy = pbuffer.data(idx_dip_gf + 443);

    auto tr_z_zzzz_xyz = pbuffer.data(idx_dip_gf + 444);

    auto tr_z_zzzz_xzz = pbuffer.data(idx_dip_gf + 445);

    auto tr_z_zzzz_yyy = pbuffer.data(idx_dip_gf + 446);

    auto tr_z_zzzz_yyz = pbuffer.data(idx_dip_gf + 447);

    auto tr_z_zzzz_yzz = pbuffer.data(idx_dip_gf + 448);

    auto tr_z_zzzz_zzz = pbuffer.data(idx_dip_gf + 449);

    // Set up components of auxiliary buffer : GG

    auto ts_xxxx_xxxx = pbuffer.data(idx_ovl_gg);

    auto ts_xxxx_xxxy = pbuffer.data(idx_ovl_gg + 1);

    auto ts_xxxx_xxxz = pbuffer.data(idx_ovl_gg + 2);

    auto ts_xxxx_xxyy = pbuffer.data(idx_ovl_gg + 3);

    auto ts_xxxx_xxyz = pbuffer.data(idx_ovl_gg + 4);

    auto ts_xxxx_xxzz = pbuffer.data(idx_ovl_gg + 5);

    auto ts_xxxx_xyyy = pbuffer.data(idx_ovl_gg + 6);

    auto ts_xxxx_xyyz = pbuffer.data(idx_ovl_gg + 7);

    auto ts_xxxx_xyzz = pbuffer.data(idx_ovl_gg + 8);

    auto ts_xxxx_xzzz = pbuffer.data(idx_ovl_gg + 9);

    auto ts_xxxx_yyyy = pbuffer.data(idx_ovl_gg + 10);

    auto ts_xxxx_yyyz = pbuffer.data(idx_ovl_gg + 11);

    auto ts_xxxx_yyzz = pbuffer.data(idx_ovl_gg + 12);

    auto ts_xxxx_yzzz = pbuffer.data(idx_ovl_gg + 13);

    auto ts_xxxx_zzzz = pbuffer.data(idx_ovl_gg + 14);

    auto ts_xxxz_xxxz = pbuffer.data(idx_ovl_gg + 32);

    auto ts_xxxz_xxzz = pbuffer.data(idx_ovl_gg + 35);

    auto ts_xxxz_xzzz = pbuffer.data(idx_ovl_gg + 39);

    auto ts_xxyy_xxxy = pbuffer.data(idx_ovl_gg + 46);

    auto ts_xxyy_xxyy = pbuffer.data(idx_ovl_gg + 48);

    auto ts_xxyy_xyyy = pbuffer.data(idx_ovl_gg + 51);

    auto ts_xxyy_yyyy = pbuffer.data(idx_ovl_gg + 55);

    auto ts_xxyy_yyyz = pbuffer.data(idx_ovl_gg + 56);

    auto ts_xxyy_yyzz = pbuffer.data(idx_ovl_gg + 57);

    auto ts_xxyy_yzzz = pbuffer.data(idx_ovl_gg + 58);

    auto ts_xxzz_xxxx = pbuffer.data(idx_ovl_gg + 75);

    auto ts_xxzz_xxxz = pbuffer.data(idx_ovl_gg + 77);

    auto ts_xxzz_xxzz = pbuffer.data(idx_ovl_gg + 80);

    auto ts_xxzz_xzzz = pbuffer.data(idx_ovl_gg + 84);

    auto ts_xxzz_yyyz = pbuffer.data(idx_ovl_gg + 86);

    auto ts_xxzz_yyzz = pbuffer.data(idx_ovl_gg + 87);

    auto ts_xxzz_yzzz = pbuffer.data(idx_ovl_gg + 88);

    auto ts_xxzz_zzzz = pbuffer.data(idx_ovl_gg + 89);

    auto ts_xyyy_yyyy = pbuffer.data(idx_ovl_gg + 100);

    auto ts_xyyy_yyyz = pbuffer.data(idx_ovl_gg + 101);

    auto ts_xyyy_yyzz = pbuffer.data(idx_ovl_gg + 102);

    auto ts_xyyy_yzzz = pbuffer.data(idx_ovl_gg + 103);

    auto ts_xzzz_yyyz = pbuffer.data(idx_ovl_gg + 146);

    auto ts_xzzz_yyzz = pbuffer.data(idx_ovl_gg + 147);

    auto ts_xzzz_yzzz = pbuffer.data(idx_ovl_gg + 148);

    auto ts_xzzz_zzzz = pbuffer.data(idx_ovl_gg + 149);

    auto ts_yyyy_xxxx = pbuffer.data(idx_ovl_gg + 150);

    auto ts_yyyy_xxxy = pbuffer.data(idx_ovl_gg + 151);

    auto ts_yyyy_xxxz = pbuffer.data(idx_ovl_gg + 152);

    auto ts_yyyy_xxyy = pbuffer.data(idx_ovl_gg + 153);

    auto ts_yyyy_xxyz = pbuffer.data(idx_ovl_gg + 154);

    auto ts_yyyy_xxzz = pbuffer.data(idx_ovl_gg + 155);

    auto ts_yyyy_xyyy = pbuffer.data(idx_ovl_gg + 156);

    auto ts_yyyy_xyyz = pbuffer.data(idx_ovl_gg + 157);

    auto ts_yyyy_xyzz = pbuffer.data(idx_ovl_gg + 158);

    auto ts_yyyy_xzzz = pbuffer.data(idx_ovl_gg + 159);

    auto ts_yyyy_yyyy = pbuffer.data(idx_ovl_gg + 160);

    auto ts_yyyy_yyyz = pbuffer.data(idx_ovl_gg + 161);

    auto ts_yyyy_yyzz = pbuffer.data(idx_ovl_gg + 162);

    auto ts_yyyy_yzzz = pbuffer.data(idx_ovl_gg + 163);

    auto ts_yyyy_zzzz = pbuffer.data(idx_ovl_gg + 164);

    auto ts_yyyz_yyyz = pbuffer.data(idx_ovl_gg + 176);

    auto ts_yyyz_yyzz = pbuffer.data(idx_ovl_gg + 177);

    auto ts_yyyz_yzzz = pbuffer.data(idx_ovl_gg + 178);

    auto ts_yyyz_zzzz = pbuffer.data(idx_ovl_gg + 179);

    auto ts_yyzz_xxxz = pbuffer.data(idx_ovl_gg + 182);

    auto ts_yyzz_xxyz = pbuffer.data(idx_ovl_gg + 184);

    auto ts_yyzz_xxzz = pbuffer.data(idx_ovl_gg + 185);

    auto ts_yyzz_xyyz = pbuffer.data(idx_ovl_gg + 187);

    auto ts_yyzz_xyzz = pbuffer.data(idx_ovl_gg + 188);

    auto ts_yyzz_xzzz = pbuffer.data(idx_ovl_gg + 189);

    auto ts_yyzz_yyyy = pbuffer.data(idx_ovl_gg + 190);

    auto ts_yyzz_yyyz = pbuffer.data(idx_ovl_gg + 191);

    auto ts_yyzz_yyzz = pbuffer.data(idx_ovl_gg + 192);

    auto ts_yyzz_yzzz = pbuffer.data(idx_ovl_gg + 193);

    auto ts_yyzz_zzzz = pbuffer.data(idx_ovl_gg + 194);

    auto ts_yzzz_xxxz = pbuffer.data(idx_ovl_gg + 197);

    auto ts_yzzz_xxzz = pbuffer.data(idx_ovl_gg + 200);

    auto ts_yzzz_xzzz = pbuffer.data(idx_ovl_gg + 204);

    auto ts_yzzz_yyyy = pbuffer.data(idx_ovl_gg + 205);

    auto ts_yzzz_yyyz = pbuffer.data(idx_ovl_gg + 206);

    auto ts_yzzz_yyzz = pbuffer.data(idx_ovl_gg + 207);

    auto ts_yzzz_yzzz = pbuffer.data(idx_ovl_gg + 208);

    auto ts_yzzz_zzzz = pbuffer.data(idx_ovl_gg + 209);

    auto ts_zzzz_xxxx = pbuffer.data(idx_ovl_gg + 210);

    auto ts_zzzz_xxxy = pbuffer.data(idx_ovl_gg + 211);

    auto ts_zzzz_xxxz = pbuffer.data(idx_ovl_gg + 212);

    auto ts_zzzz_xxyy = pbuffer.data(idx_ovl_gg + 213);

    auto ts_zzzz_xxyz = pbuffer.data(idx_ovl_gg + 214);

    auto ts_zzzz_xxzz = pbuffer.data(idx_ovl_gg + 215);

    auto ts_zzzz_xyyy = pbuffer.data(idx_ovl_gg + 216);

    auto ts_zzzz_xyyz = pbuffer.data(idx_ovl_gg + 217);

    auto ts_zzzz_xyzz = pbuffer.data(idx_ovl_gg + 218);

    auto ts_zzzz_xzzz = pbuffer.data(idx_ovl_gg + 219);

    auto ts_zzzz_yyyy = pbuffer.data(idx_ovl_gg + 220);

    auto ts_zzzz_yyyz = pbuffer.data(idx_ovl_gg + 221);

    auto ts_zzzz_yyzz = pbuffer.data(idx_ovl_gg + 222);

    auto ts_zzzz_yzzz = pbuffer.data(idx_ovl_gg + 223);

    auto ts_zzzz_zzzz = pbuffer.data(idx_ovl_gg + 224);

    // Set up components of auxiliary buffer : GG

    auto tr_x_xxxx_xxxx = pbuffer.data(idx_dip_gg);

    auto tr_x_xxxx_xxxy = pbuffer.data(idx_dip_gg + 1);

    auto tr_x_xxxx_xxxz = pbuffer.data(idx_dip_gg + 2);

    auto tr_x_xxxx_xxyy = pbuffer.data(idx_dip_gg + 3);

    auto tr_x_xxxx_xxyz = pbuffer.data(idx_dip_gg + 4);

    auto tr_x_xxxx_xxzz = pbuffer.data(idx_dip_gg + 5);

    auto tr_x_xxxx_xyyy = pbuffer.data(idx_dip_gg + 6);

    auto tr_x_xxxx_xyyz = pbuffer.data(idx_dip_gg + 7);

    auto tr_x_xxxx_xyzz = pbuffer.data(idx_dip_gg + 8);

    auto tr_x_xxxx_xzzz = pbuffer.data(idx_dip_gg + 9);

    auto tr_x_xxxx_yyyy = pbuffer.data(idx_dip_gg + 10);

    auto tr_x_xxxx_yyyz = pbuffer.data(idx_dip_gg + 11);

    auto tr_x_xxxx_yyzz = pbuffer.data(idx_dip_gg + 12);

    auto tr_x_xxxx_yzzz = pbuffer.data(idx_dip_gg + 13);

    auto tr_x_xxxx_zzzz = pbuffer.data(idx_dip_gg + 14);

    auto tr_x_xxxy_xxxx = pbuffer.data(idx_dip_gg + 15);

    auto tr_x_xxxy_xxxy = pbuffer.data(idx_dip_gg + 16);

    auto tr_x_xxxy_xxxz = pbuffer.data(idx_dip_gg + 17);

    auto tr_x_xxxy_xxyy = pbuffer.data(idx_dip_gg + 18);

    auto tr_x_xxxy_xxyz = pbuffer.data(idx_dip_gg + 19);

    auto tr_x_xxxy_xxzz = pbuffer.data(idx_dip_gg + 20);

    auto tr_x_xxxy_xyyy = pbuffer.data(idx_dip_gg + 21);

    auto tr_x_xxxy_xyyz = pbuffer.data(idx_dip_gg + 22);

    auto tr_x_xxxy_xyzz = pbuffer.data(idx_dip_gg + 23);

    auto tr_x_xxxy_xzzz = pbuffer.data(idx_dip_gg + 24);

    auto tr_x_xxxy_yyyy = pbuffer.data(idx_dip_gg + 25);

    auto tr_x_xxxy_zzzz = pbuffer.data(idx_dip_gg + 29);

    auto tr_x_xxxz_xxxx = pbuffer.data(idx_dip_gg + 30);

    auto tr_x_xxxz_xxxy = pbuffer.data(idx_dip_gg + 31);

    auto tr_x_xxxz_xxxz = pbuffer.data(idx_dip_gg + 32);

    auto tr_x_xxxz_xxyy = pbuffer.data(idx_dip_gg + 33);

    auto tr_x_xxxz_xxyz = pbuffer.data(idx_dip_gg + 34);

    auto tr_x_xxxz_xxzz = pbuffer.data(idx_dip_gg + 35);

    auto tr_x_xxxz_xyyy = pbuffer.data(idx_dip_gg + 36);

    auto tr_x_xxxz_xyyz = pbuffer.data(idx_dip_gg + 37);

    auto tr_x_xxxz_xyzz = pbuffer.data(idx_dip_gg + 38);

    auto tr_x_xxxz_xzzz = pbuffer.data(idx_dip_gg + 39);

    auto tr_x_xxxz_yyyy = pbuffer.data(idx_dip_gg + 40);

    auto tr_x_xxxz_yyyz = pbuffer.data(idx_dip_gg + 41);

    auto tr_x_xxxz_yyzz = pbuffer.data(idx_dip_gg + 42);

    auto tr_x_xxxz_yzzz = pbuffer.data(idx_dip_gg + 43);

    auto tr_x_xxxz_zzzz = pbuffer.data(idx_dip_gg + 44);

    auto tr_x_xxyy_xxxx = pbuffer.data(idx_dip_gg + 45);

    auto tr_x_xxyy_xxxy = pbuffer.data(idx_dip_gg + 46);

    auto tr_x_xxyy_xxxz = pbuffer.data(idx_dip_gg + 47);

    auto tr_x_xxyy_xxyy = pbuffer.data(idx_dip_gg + 48);

    auto tr_x_xxyy_xxyz = pbuffer.data(idx_dip_gg + 49);

    auto tr_x_xxyy_xxzz = pbuffer.data(idx_dip_gg + 50);

    auto tr_x_xxyy_xyyy = pbuffer.data(idx_dip_gg + 51);

    auto tr_x_xxyy_xyyz = pbuffer.data(idx_dip_gg + 52);

    auto tr_x_xxyy_xyzz = pbuffer.data(idx_dip_gg + 53);

    auto tr_x_xxyy_xzzz = pbuffer.data(idx_dip_gg + 54);

    auto tr_x_xxyy_yyyy = pbuffer.data(idx_dip_gg + 55);

    auto tr_x_xxyy_yyyz = pbuffer.data(idx_dip_gg + 56);

    auto tr_x_xxyy_yyzz = pbuffer.data(idx_dip_gg + 57);

    auto tr_x_xxyy_yzzz = pbuffer.data(idx_dip_gg + 58);

    auto tr_x_xxyy_zzzz = pbuffer.data(idx_dip_gg + 59);

    auto tr_x_xxyz_xxxz = pbuffer.data(idx_dip_gg + 62);

    auto tr_x_xxyz_xxzz = pbuffer.data(idx_dip_gg + 65);

    auto tr_x_xxyz_xzzz = pbuffer.data(idx_dip_gg + 69);

    auto tr_x_xxyz_zzzz = pbuffer.data(idx_dip_gg + 74);

    auto tr_x_xxzz_xxxx = pbuffer.data(idx_dip_gg + 75);

    auto tr_x_xxzz_xxxy = pbuffer.data(idx_dip_gg + 76);

    auto tr_x_xxzz_xxxz = pbuffer.data(idx_dip_gg + 77);

    auto tr_x_xxzz_xxyy = pbuffer.data(idx_dip_gg + 78);

    auto tr_x_xxzz_xxyz = pbuffer.data(idx_dip_gg + 79);

    auto tr_x_xxzz_xxzz = pbuffer.data(idx_dip_gg + 80);

    auto tr_x_xxzz_xyyy = pbuffer.data(idx_dip_gg + 81);

    auto tr_x_xxzz_xyyz = pbuffer.data(idx_dip_gg + 82);

    auto tr_x_xxzz_xyzz = pbuffer.data(idx_dip_gg + 83);

    auto tr_x_xxzz_xzzz = pbuffer.data(idx_dip_gg + 84);

    auto tr_x_xxzz_yyyy = pbuffer.data(idx_dip_gg + 85);

    auto tr_x_xxzz_yyyz = pbuffer.data(idx_dip_gg + 86);

    auto tr_x_xxzz_yyzz = pbuffer.data(idx_dip_gg + 87);

    auto tr_x_xxzz_yzzz = pbuffer.data(idx_dip_gg + 88);

    auto tr_x_xxzz_zzzz = pbuffer.data(idx_dip_gg + 89);

    auto tr_x_xyyy_xxxx = pbuffer.data(idx_dip_gg + 90);

    auto tr_x_xyyy_xxxy = pbuffer.data(idx_dip_gg + 91);

    auto tr_x_xyyy_xxxz = pbuffer.data(idx_dip_gg + 92);

    auto tr_x_xyyy_xxyy = pbuffer.data(idx_dip_gg + 93);

    auto tr_x_xyyy_xxyz = pbuffer.data(idx_dip_gg + 94);

    auto tr_x_xyyy_xxzz = pbuffer.data(idx_dip_gg + 95);

    auto tr_x_xyyy_xyyy = pbuffer.data(idx_dip_gg + 96);

    auto tr_x_xyyy_xyyz = pbuffer.data(idx_dip_gg + 97);

    auto tr_x_xyyy_xyzz = pbuffer.data(idx_dip_gg + 98);

    auto tr_x_xyyy_xzzz = pbuffer.data(idx_dip_gg + 99);

    auto tr_x_xyyy_yyyy = pbuffer.data(idx_dip_gg + 100);

    auto tr_x_xyyy_yyyz = pbuffer.data(idx_dip_gg + 101);

    auto tr_x_xyyy_yyzz = pbuffer.data(idx_dip_gg + 102);

    auto tr_x_xyyy_yzzz = pbuffer.data(idx_dip_gg + 103);

    auto tr_x_xyyz_xxxy = pbuffer.data(idx_dip_gg + 106);

    auto tr_x_xyyz_xxxz = pbuffer.data(idx_dip_gg + 107);

    auto tr_x_xyyz_xxyy = pbuffer.data(idx_dip_gg + 108);

    auto tr_x_xyyz_xxzz = pbuffer.data(idx_dip_gg + 110);

    auto tr_x_xyyz_xyyy = pbuffer.data(idx_dip_gg + 111);

    auto tr_x_xyyz_xzzz = pbuffer.data(idx_dip_gg + 114);

    auto tr_x_xyzz_xxxx = pbuffer.data(idx_dip_gg + 120);

    auto tr_x_xyzz_xxxz = pbuffer.data(idx_dip_gg + 122);

    auto tr_x_xyzz_xxzz = pbuffer.data(idx_dip_gg + 125);

    auto tr_x_xyzz_xzzz = pbuffer.data(idx_dip_gg + 129);

    auto tr_x_xzzz_xxxx = pbuffer.data(idx_dip_gg + 135);

    auto tr_x_xzzz_xxxy = pbuffer.data(idx_dip_gg + 136);

    auto tr_x_xzzz_xxxz = pbuffer.data(idx_dip_gg + 137);

    auto tr_x_xzzz_xxyy = pbuffer.data(idx_dip_gg + 138);

    auto tr_x_xzzz_xxyz = pbuffer.data(idx_dip_gg + 139);

    auto tr_x_xzzz_xxzz = pbuffer.data(idx_dip_gg + 140);

    auto tr_x_xzzz_xyyy = pbuffer.data(idx_dip_gg + 141);

    auto tr_x_xzzz_xyyz = pbuffer.data(idx_dip_gg + 142);

    auto tr_x_xzzz_xyzz = pbuffer.data(idx_dip_gg + 143);

    auto tr_x_xzzz_xzzz = pbuffer.data(idx_dip_gg + 144);

    auto tr_x_xzzz_yyyz = pbuffer.data(idx_dip_gg + 146);

    auto tr_x_xzzz_yyzz = pbuffer.data(idx_dip_gg + 147);

    auto tr_x_xzzz_yzzz = pbuffer.data(idx_dip_gg + 148);

    auto tr_x_xzzz_zzzz = pbuffer.data(idx_dip_gg + 149);

    auto tr_x_yyyy_xxxx = pbuffer.data(idx_dip_gg + 150);

    auto tr_x_yyyy_xxxy = pbuffer.data(idx_dip_gg + 151);

    auto tr_x_yyyy_xxxz = pbuffer.data(idx_dip_gg + 152);

    auto tr_x_yyyy_xxyy = pbuffer.data(idx_dip_gg + 153);

    auto tr_x_yyyy_xxyz = pbuffer.data(idx_dip_gg + 154);

    auto tr_x_yyyy_xxzz = pbuffer.data(idx_dip_gg + 155);

    auto tr_x_yyyy_xyyy = pbuffer.data(idx_dip_gg + 156);

    auto tr_x_yyyy_xyyz = pbuffer.data(idx_dip_gg + 157);

    auto tr_x_yyyy_xyzz = pbuffer.data(idx_dip_gg + 158);

    auto tr_x_yyyy_xzzz = pbuffer.data(idx_dip_gg + 159);

    auto tr_x_yyyy_yyyy = pbuffer.data(idx_dip_gg + 160);

    auto tr_x_yyyy_yyyz = pbuffer.data(idx_dip_gg + 161);

    auto tr_x_yyyy_yyzz = pbuffer.data(idx_dip_gg + 162);

    auto tr_x_yyyy_yzzz = pbuffer.data(idx_dip_gg + 163);

    auto tr_x_yyyy_zzzz = pbuffer.data(idx_dip_gg + 164);

    auto tr_x_yyyz_xxxy = pbuffer.data(idx_dip_gg + 166);

    auto tr_x_yyyz_xxxz = pbuffer.data(idx_dip_gg + 167);

    auto tr_x_yyyz_xxyy = pbuffer.data(idx_dip_gg + 168);

    auto tr_x_yyyz_xxzz = pbuffer.data(idx_dip_gg + 170);

    auto tr_x_yyyz_xyyy = pbuffer.data(idx_dip_gg + 171);

    auto tr_x_yyyz_xzzz = pbuffer.data(idx_dip_gg + 174);

    auto tr_x_yyyz_yyyy = pbuffer.data(idx_dip_gg + 175);

    auto tr_x_yyyz_yyyz = pbuffer.data(idx_dip_gg + 176);

    auto tr_x_yyyz_yyzz = pbuffer.data(idx_dip_gg + 177);

    auto tr_x_yyyz_yzzz = pbuffer.data(idx_dip_gg + 178);

    auto tr_x_yyyz_zzzz = pbuffer.data(idx_dip_gg + 179);

    auto tr_x_yyzz_xxxx = pbuffer.data(idx_dip_gg + 180);

    auto tr_x_yyzz_xxxy = pbuffer.data(idx_dip_gg + 181);

    auto tr_x_yyzz_xxxz = pbuffer.data(idx_dip_gg + 182);

    auto tr_x_yyzz_xxyy = pbuffer.data(idx_dip_gg + 183);

    auto tr_x_yyzz_xxyz = pbuffer.data(idx_dip_gg + 184);

    auto tr_x_yyzz_xxzz = pbuffer.data(idx_dip_gg + 185);

    auto tr_x_yyzz_xyyy = pbuffer.data(idx_dip_gg + 186);

    auto tr_x_yyzz_xyyz = pbuffer.data(idx_dip_gg + 187);

    auto tr_x_yyzz_xyzz = pbuffer.data(idx_dip_gg + 188);

    auto tr_x_yyzz_xzzz = pbuffer.data(idx_dip_gg + 189);

    auto tr_x_yyzz_yyyy = pbuffer.data(idx_dip_gg + 190);

    auto tr_x_yyzz_yyyz = pbuffer.data(idx_dip_gg + 191);

    auto tr_x_yyzz_yyzz = pbuffer.data(idx_dip_gg + 192);

    auto tr_x_yyzz_yzzz = pbuffer.data(idx_dip_gg + 193);

    auto tr_x_yyzz_zzzz = pbuffer.data(idx_dip_gg + 194);

    auto tr_x_yzzz_xxxx = pbuffer.data(idx_dip_gg + 195);

    auto tr_x_yzzz_xxxz = pbuffer.data(idx_dip_gg + 197);

    auto tr_x_yzzz_xxyz = pbuffer.data(idx_dip_gg + 199);

    auto tr_x_yzzz_xxzz = pbuffer.data(idx_dip_gg + 200);

    auto tr_x_yzzz_xyyz = pbuffer.data(idx_dip_gg + 202);

    auto tr_x_yzzz_xyzz = pbuffer.data(idx_dip_gg + 203);

    auto tr_x_yzzz_xzzz = pbuffer.data(idx_dip_gg + 204);

    auto tr_x_yzzz_yyyy = pbuffer.data(idx_dip_gg + 205);

    auto tr_x_yzzz_yyyz = pbuffer.data(idx_dip_gg + 206);

    auto tr_x_yzzz_yyzz = pbuffer.data(idx_dip_gg + 207);

    auto tr_x_yzzz_yzzz = pbuffer.data(idx_dip_gg + 208);

    auto tr_x_yzzz_zzzz = pbuffer.data(idx_dip_gg + 209);

    auto tr_x_zzzz_xxxx = pbuffer.data(idx_dip_gg + 210);

    auto tr_x_zzzz_xxxy = pbuffer.data(idx_dip_gg + 211);

    auto tr_x_zzzz_xxxz = pbuffer.data(idx_dip_gg + 212);

    auto tr_x_zzzz_xxyy = pbuffer.data(idx_dip_gg + 213);

    auto tr_x_zzzz_xxyz = pbuffer.data(idx_dip_gg + 214);

    auto tr_x_zzzz_xxzz = pbuffer.data(idx_dip_gg + 215);

    auto tr_x_zzzz_xyyy = pbuffer.data(idx_dip_gg + 216);

    auto tr_x_zzzz_xyyz = pbuffer.data(idx_dip_gg + 217);

    auto tr_x_zzzz_xyzz = pbuffer.data(idx_dip_gg + 218);

    auto tr_x_zzzz_xzzz = pbuffer.data(idx_dip_gg + 219);

    auto tr_x_zzzz_yyyy = pbuffer.data(idx_dip_gg + 220);

    auto tr_x_zzzz_yyyz = pbuffer.data(idx_dip_gg + 221);

    auto tr_x_zzzz_yyzz = pbuffer.data(idx_dip_gg + 222);

    auto tr_x_zzzz_yzzz = pbuffer.data(idx_dip_gg + 223);

    auto tr_x_zzzz_zzzz = pbuffer.data(idx_dip_gg + 224);

    auto tr_y_xxxx_xxxx = pbuffer.data(idx_dip_gg + 225);

    auto tr_y_xxxx_xxxy = pbuffer.data(idx_dip_gg + 226);

    auto tr_y_xxxx_xxxz = pbuffer.data(idx_dip_gg + 227);

    auto tr_y_xxxx_xxyy = pbuffer.data(idx_dip_gg + 228);

    auto tr_y_xxxx_xxyz = pbuffer.data(idx_dip_gg + 229);

    auto tr_y_xxxx_xxzz = pbuffer.data(idx_dip_gg + 230);

    auto tr_y_xxxx_xyyy = pbuffer.data(idx_dip_gg + 231);

    auto tr_y_xxxx_xyyz = pbuffer.data(idx_dip_gg + 232);

    auto tr_y_xxxx_xyzz = pbuffer.data(idx_dip_gg + 233);

    auto tr_y_xxxx_xzzz = pbuffer.data(idx_dip_gg + 234);

    auto tr_y_xxxx_yyyy = pbuffer.data(idx_dip_gg + 235);

    auto tr_y_xxxx_yyyz = pbuffer.data(idx_dip_gg + 236);

    auto tr_y_xxxx_yyzz = pbuffer.data(idx_dip_gg + 237);

    auto tr_y_xxxx_yzzz = pbuffer.data(idx_dip_gg + 238);

    auto tr_y_xxxx_zzzz = pbuffer.data(idx_dip_gg + 239);

    auto tr_y_xxxy_xxxx = pbuffer.data(idx_dip_gg + 240);

    auto tr_y_xxxy_xxxy = pbuffer.data(idx_dip_gg + 241);

    auto tr_y_xxxy_xxyy = pbuffer.data(idx_dip_gg + 243);

    auto tr_y_xxxy_xxyz = pbuffer.data(idx_dip_gg + 244);

    auto tr_y_xxxy_xyyy = pbuffer.data(idx_dip_gg + 246);

    auto tr_y_xxxy_xyyz = pbuffer.data(idx_dip_gg + 247);

    auto tr_y_xxxy_xyzz = pbuffer.data(idx_dip_gg + 248);

    auto tr_y_xxxy_yyyy = pbuffer.data(idx_dip_gg + 250);

    auto tr_y_xxxy_yyyz = pbuffer.data(idx_dip_gg + 251);

    auto tr_y_xxxy_yyzz = pbuffer.data(idx_dip_gg + 252);

    auto tr_y_xxxy_yzzz = pbuffer.data(idx_dip_gg + 253);

    auto tr_y_xxxy_zzzz = pbuffer.data(idx_dip_gg + 254);

    auto tr_y_xxxz_xxxx = pbuffer.data(idx_dip_gg + 255);

    auto tr_y_xxxz_xxxy = pbuffer.data(idx_dip_gg + 256);

    auto tr_y_xxxz_xxxz = pbuffer.data(idx_dip_gg + 257);

    auto tr_y_xxxz_xxyy = pbuffer.data(idx_dip_gg + 258);

    auto tr_y_xxxz_xxzz = pbuffer.data(idx_dip_gg + 260);

    auto tr_y_xxxz_xyyy = pbuffer.data(idx_dip_gg + 261);

    auto tr_y_xxxz_xzzz = pbuffer.data(idx_dip_gg + 264);

    auto tr_y_xxxz_yyyz = pbuffer.data(idx_dip_gg + 266);

    auto tr_y_xxxz_yyzz = pbuffer.data(idx_dip_gg + 267);

    auto tr_y_xxxz_yzzz = pbuffer.data(idx_dip_gg + 268);

    auto tr_y_xxxz_zzzz = pbuffer.data(idx_dip_gg + 269);

    auto tr_y_xxyy_xxxx = pbuffer.data(idx_dip_gg + 270);

    auto tr_y_xxyy_xxxy = pbuffer.data(idx_dip_gg + 271);

    auto tr_y_xxyy_xxxz = pbuffer.data(idx_dip_gg + 272);

    auto tr_y_xxyy_xxyy = pbuffer.data(idx_dip_gg + 273);

    auto tr_y_xxyy_xxyz = pbuffer.data(idx_dip_gg + 274);

    auto tr_y_xxyy_xxzz = pbuffer.data(idx_dip_gg + 275);

    auto tr_y_xxyy_xyyy = pbuffer.data(idx_dip_gg + 276);

    auto tr_y_xxyy_xyyz = pbuffer.data(idx_dip_gg + 277);

    auto tr_y_xxyy_xyzz = pbuffer.data(idx_dip_gg + 278);

    auto tr_y_xxyy_xzzz = pbuffer.data(idx_dip_gg + 279);

    auto tr_y_xxyy_yyyy = pbuffer.data(idx_dip_gg + 280);

    auto tr_y_xxyy_yyyz = pbuffer.data(idx_dip_gg + 281);

    auto tr_y_xxyy_yyzz = pbuffer.data(idx_dip_gg + 282);

    auto tr_y_xxyy_yzzz = pbuffer.data(idx_dip_gg + 283);

    auto tr_y_xxyy_zzzz = pbuffer.data(idx_dip_gg + 284);

    auto tr_y_xxyz_xxxy = pbuffer.data(idx_dip_gg + 286);

    auto tr_y_xxyz_xxyy = pbuffer.data(idx_dip_gg + 288);

    auto tr_y_xxyz_xyyy = pbuffer.data(idx_dip_gg + 291);

    auto tr_y_xxyz_yyyz = pbuffer.data(idx_dip_gg + 296);

    auto tr_y_xxyz_yyzz = pbuffer.data(idx_dip_gg + 297);

    auto tr_y_xxyz_yzzz = pbuffer.data(idx_dip_gg + 298);

    auto tr_y_xxyz_zzzz = pbuffer.data(idx_dip_gg + 299);

    auto tr_y_xxzz_xxxx = pbuffer.data(idx_dip_gg + 300);

    auto tr_y_xxzz_xxxy = pbuffer.data(idx_dip_gg + 301);

    auto tr_y_xxzz_xxxz = pbuffer.data(idx_dip_gg + 302);

    auto tr_y_xxzz_xxyy = pbuffer.data(idx_dip_gg + 303);

    auto tr_y_xxzz_xxyz = pbuffer.data(idx_dip_gg + 304);

    auto tr_y_xxzz_xxzz = pbuffer.data(idx_dip_gg + 305);

    auto tr_y_xxzz_xyyy = pbuffer.data(idx_dip_gg + 306);

    auto tr_y_xxzz_xyyz = pbuffer.data(idx_dip_gg + 307);

    auto tr_y_xxzz_xyzz = pbuffer.data(idx_dip_gg + 308);

    auto tr_y_xxzz_xzzz = pbuffer.data(idx_dip_gg + 309);

    auto tr_y_xxzz_yyyy = pbuffer.data(idx_dip_gg + 310);

    auto tr_y_xxzz_yyyz = pbuffer.data(idx_dip_gg + 311);

    auto tr_y_xxzz_yyzz = pbuffer.data(idx_dip_gg + 312);

    auto tr_y_xxzz_yzzz = pbuffer.data(idx_dip_gg + 313);

    auto tr_y_xxzz_zzzz = pbuffer.data(idx_dip_gg + 314);

    auto tr_y_xyyy_xxxx = pbuffer.data(idx_dip_gg + 315);

    auto tr_y_xyyy_xxxy = pbuffer.data(idx_dip_gg + 316);

    auto tr_y_xyyy_xxxz = pbuffer.data(idx_dip_gg + 317);

    auto tr_y_xyyy_xxyy = pbuffer.data(idx_dip_gg + 318);

    auto tr_y_xyyy_xxyz = pbuffer.data(idx_dip_gg + 319);

    auto tr_y_xyyy_xxzz = pbuffer.data(idx_dip_gg + 320);

    auto tr_y_xyyy_xyyy = pbuffer.data(idx_dip_gg + 321);

    auto tr_y_xyyy_xyyz = pbuffer.data(idx_dip_gg + 322);

    auto tr_y_xyyy_xyzz = pbuffer.data(idx_dip_gg + 323);

    auto tr_y_xyyy_xzzz = pbuffer.data(idx_dip_gg + 324);

    auto tr_y_xyyy_yyyy = pbuffer.data(idx_dip_gg + 325);

    auto tr_y_xyyy_yyyz = pbuffer.data(idx_dip_gg + 326);

    auto tr_y_xyyy_yyzz = pbuffer.data(idx_dip_gg + 327);

    auto tr_y_xyyy_yzzz = pbuffer.data(idx_dip_gg + 328);

    auto tr_y_xyyy_zzzz = pbuffer.data(idx_dip_gg + 329);

    auto tr_y_xyyz_yyyz = pbuffer.data(idx_dip_gg + 341);

    auto tr_y_xyyz_yyzz = pbuffer.data(idx_dip_gg + 342);

    auto tr_y_xyyz_yzzz = pbuffer.data(idx_dip_gg + 343);

    auto tr_y_xyyz_zzzz = pbuffer.data(idx_dip_gg + 344);

    auto tr_y_xyzz_xxyz = pbuffer.data(idx_dip_gg + 349);

    auto tr_y_xyzz_xyyz = pbuffer.data(idx_dip_gg + 352);

    auto tr_y_xyzz_xyzz = pbuffer.data(idx_dip_gg + 353);

    auto tr_y_xyzz_yyyy = pbuffer.data(idx_dip_gg + 355);

    auto tr_y_xyzz_yyyz = pbuffer.data(idx_dip_gg + 356);

    auto tr_y_xyzz_yyzz = pbuffer.data(idx_dip_gg + 357);

    auto tr_y_xyzz_yzzz = pbuffer.data(idx_dip_gg + 358);

    auto tr_y_xyzz_zzzz = pbuffer.data(idx_dip_gg + 359);

    auto tr_y_xzzz_xxxz = pbuffer.data(idx_dip_gg + 362);

    auto tr_y_xzzz_xxyz = pbuffer.data(idx_dip_gg + 364);

    auto tr_y_xzzz_xxzz = pbuffer.data(idx_dip_gg + 365);

    auto tr_y_xzzz_xyyz = pbuffer.data(idx_dip_gg + 367);

    auto tr_y_xzzz_xyzz = pbuffer.data(idx_dip_gg + 368);

    auto tr_y_xzzz_xzzz = pbuffer.data(idx_dip_gg + 369);

    auto tr_y_xzzz_yyyy = pbuffer.data(idx_dip_gg + 370);

    auto tr_y_xzzz_yyyz = pbuffer.data(idx_dip_gg + 371);

    auto tr_y_xzzz_yyzz = pbuffer.data(idx_dip_gg + 372);

    auto tr_y_xzzz_yzzz = pbuffer.data(idx_dip_gg + 373);

    auto tr_y_xzzz_zzzz = pbuffer.data(idx_dip_gg + 374);

    auto tr_y_yyyy_xxxx = pbuffer.data(idx_dip_gg + 375);

    auto tr_y_yyyy_xxxy = pbuffer.data(idx_dip_gg + 376);

    auto tr_y_yyyy_xxxz = pbuffer.data(idx_dip_gg + 377);

    auto tr_y_yyyy_xxyy = pbuffer.data(idx_dip_gg + 378);

    auto tr_y_yyyy_xxyz = pbuffer.data(idx_dip_gg + 379);

    auto tr_y_yyyy_xxzz = pbuffer.data(idx_dip_gg + 380);

    auto tr_y_yyyy_xyyy = pbuffer.data(idx_dip_gg + 381);

    auto tr_y_yyyy_xyyz = pbuffer.data(idx_dip_gg + 382);

    auto tr_y_yyyy_xyzz = pbuffer.data(idx_dip_gg + 383);

    auto tr_y_yyyy_xzzz = pbuffer.data(idx_dip_gg + 384);

    auto tr_y_yyyy_yyyy = pbuffer.data(idx_dip_gg + 385);

    auto tr_y_yyyy_yyyz = pbuffer.data(idx_dip_gg + 386);

    auto tr_y_yyyy_yyzz = pbuffer.data(idx_dip_gg + 387);

    auto tr_y_yyyy_yzzz = pbuffer.data(idx_dip_gg + 388);

    auto tr_y_yyyy_zzzz = pbuffer.data(idx_dip_gg + 389);

    auto tr_y_yyyz_xxxx = pbuffer.data(idx_dip_gg + 390);

    auto tr_y_yyyz_xxxy = pbuffer.data(idx_dip_gg + 391);

    auto tr_y_yyyz_xxxz = pbuffer.data(idx_dip_gg + 392);

    auto tr_y_yyyz_xxyy = pbuffer.data(idx_dip_gg + 393);

    auto tr_y_yyyz_xxyz = pbuffer.data(idx_dip_gg + 394);

    auto tr_y_yyyz_xxzz = pbuffer.data(idx_dip_gg + 395);

    auto tr_y_yyyz_xyyy = pbuffer.data(idx_dip_gg + 396);

    auto tr_y_yyyz_xyyz = pbuffer.data(idx_dip_gg + 397);

    auto tr_y_yyyz_xyzz = pbuffer.data(idx_dip_gg + 398);

    auto tr_y_yyyz_xzzz = pbuffer.data(idx_dip_gg + 399);

    auto tr_y_yyyz_yyyy = pbuffer.data(idx_dip_gg + 400);

    auto tr_y_yyyz_yyyz = pbuffer.data(idx_dip_gg + 401);

    auto tr_y_yyyz_yyzz = pbuffer.data(idx_dip_gg + 402);

    auto tr_y_yyyz_yzzz = pbuffer.data(idx_dip_gg + 403);

    auto tr_y_yyyz_zzzz = pbuffer.data(idx_dip_gg + 404);

    auto tr_y_yyzz_xxxx = pbuffer.data(idx_dip_gg + 405);

    auto tr_y_yyzz_xxxy = pbuffer.data(idx_dip_gg + 406);

    auto tr_y_yyzz_xxxz = pbuffer.data(idx_dip_gg + 407);

    auto tr_y_yyzz_xxyy = pbuffer.data(idx_dip_gg + 408);

    auto tr_y_yyzz_xxyz = pbuffer.data(idx_dip_gg + 409);

    auto tr_y_yyzz_xxzz = pbuffer.data(idx_dip_gg + 410);

    auto tr_y_yyzz_xyyy = pbuffer.data(idx_dip_gg + 411);

    auto tr_y_yyzz_xyyz = pbuffer.data(idx_dip_gg + 412);

    auto tr_y_yyzz_xyzz = pbuffer.data(idx_dip_gg + 413);

    auto tr_y_yyzz_xzzz = pbuffer.data(idx_dip_gg + 414);

    auto tr_y_yyzz_yyyy = pbuffer.data(idx_dip_gg + 415);

    auto tr_y_yyzz_yyyz = pbuffer.data(idx_dip_gg + 416);

    auto tr_y_yyzz_yyzz = pbuffer.data(idx_dip_gg + 417);

    auto tr_y_yyzz_yzzz = pbuffer.data(idx_dip_gg + 418);

    auto tr_y_yyzz_zzzz = pbuffer.data(idx_dip_gg + 419);

    auto tr_y_yzzz_xxxx = pbuffer.data(idx_dip_gg + 420);

    auto tr_y_yzzz_xxxy = pbuffer.data(idx_dip_gg + 421);

    auto tr_y_yzzz_xxxz = pbuffer.data(idx_dip_gg + 422);

    auto tr_y_yzzz_xxyy = pbuffer.data(idx_dip_gg + 423);

    auto tr_y_yzzz_xxyz = pbuffer.data(idx_dip_gg + 424);

    auto tr_y_yzzz_xxzz = pbuffer.data(idx_dip_gg + 425);

    auto tr_y_yzzz_xyyy = pbuffer.data(idx_dip_gg + 426);

    auto tr_y_yzzz_xyyz = pbuffer.data(idx_dip_gg + 427);

    auto tr_y_yzzz_xyzz = pbuffer.data(idx_dip_gg + 428);

    auto tr_y_yzzz_xzzz = pbuffer.data(idx_dip_gg + 429);

    auto tr_y_yzzz_yyyy = pbuffer.data(idx_dip_gg + 430);

    auto tr_y_yzzz_yyyz = pbuffer.data(idx_dip_gg + 431);

    auto tr_y_yzzz_yyzz = pbuffer.data(idx_dip_gg + 432);

    auto tr_y_yzzz_yzzz = pbuffer.data(idx_dip_gg + 433);

    auto tr_y_yzzz_zzzz = pbuffer.data(idx_dip_gg + 434);

    auto tr_y_zzzz_xxxx = pbuffer.data(idx_dip_gg + 435);

    auto tr_y_zzzz_xxxy = pbuffer.data(idx_dip_gg + 436);

    auto tr_y_zzzz_xxxz = pbuffer.data(idx_dip_gg + 437);

    auto tr_y_zzzz_xxyy = pbuffer.data(idx_dip_gg + 438);

    auto tr_y_zzzz_xxyz = pbuffer.data(idx_dip_gg + 439);

    auto tr_y_zzzz_xxzz = pbuffer.data(idx_dip_gg + 440);

    auto tr_y_zzzz_xyyy = pbuffer.data(idx_dip_gg + 441);

    auto tr_y_zzzz_xyyz = pbuffer.data(idx_dip_gg + 442);

    auto tr_y_zzzz_xyzz = pbuffer.data(idx_dip_gg + 443);

    auto tr_y_zzzz_xzzz = pbuffer.data(idx_dip_gg + 444);

    auto tr_y_zzzz_yyyy = pbuffer.data(idx_dip_gg + 445);

    auto tr_y_zzzz_yyyz = pbuffer.data(idx_dip_gg + 446);

    auto tr_y_zzzz_yyzz = pbuffer.data(idx_dip_gg + 447);

    auto tr_y_zzzz_yzzz = pbuffer.data(idx_dip_gg + 448);

    auto tr_y_zzzz_zzzz = pbuffer.data(idx_dip_gg + 449);

    auto tr_z_xxxx_xxxx = pbuffer.data(idx_dip_gg + 450);

    auto tr_z_xxxx_xxxy = pbuffer.data(idx_dip_gg + 451);

    auto tr_z_xxxx_xxxz = pbuffer.data(idx_dip_gg + 452);

    auto tr_z_xxxx_xxyy = pbuffer.data(idx_dip_gg + 453);

    auto tr_z_xxxx_xxyz = pbuffer.data(idx_dip_gg + 454);

    auto tr_z_xxxx_xxzz = pbuffer.data(idx_dip_gg + 455);

    auto tr_z_xxxx_xyyy = pbuffer.data(idx_dip_gg + 456);

    auto tr_z_xxxx_xyyz = pbuffer.data(idx_dip_gg + 457);

    auto tr_z_xxxx_xyzz = pbuffer.data(idx_dip_gg + 458);

    auto tr_z_xxxx_xzzz = pbuffer.data(idx_dip_gg + 459);

    auto tr_z_xxxx_yyyy = pbuffer.data(idx_dip_gg + 460);

    auto tr_z_xxxx_yyyz = pbuffer.data(idx_dip_gg + 461);

    auto tr_z_xxxx_yyzz = pbuffer.data(idx_dip_gg + 462);

    auto tr_z_xxxx_yzzz = pbuffer.data(idx_dip_gg + 463);

    auto tr_z_xxxx_zzzz = pbuffer.data(idx_dip_gg + 464);

    auto tr_z_xxxy_xxxx = pbuffer.data(idx_dip_gg + 465);

    auto tr_z_xxxy_xxxz = pbuffer.data(idx_dip_gg + 467);

    auto tr_z_xxxy_xxzz = pbuffer.data(idx_dip_gg + 470);

    auto tr_z_xxxy_xzzz = pbuffer.data(idx_dip_gg + 474);

    auto tr_z_xxxy_yyyy = pbuffer.data(idx_dip_gg + 475);

    auto tr_z_xxxy_yyyz = pbuffer.data(idx_dip_gg + 476);

    auto tr_z_xxxy_yyzz = pbuffer.data(idx_dip_gg + 477);

    auto tr_z_xxxy_yzzz = pbuffer.data(idx_dip_gg + 478);

    auto tr_z_xxxz_xxxx = pbuffer.data(idx_dip_gg + 480);

    auto tr_z_xxxz_xxxy = pbuffer.data(idx_dip_gg + 481);

    auto tr_z_xxxz_xxxz = pbuffer.data(idx_dip_gg + 482);

    auto tr_z_xxxz_xxyy = pbuffer.data(idx_dip_gg + 483);

    auto tr_z_xxxz_xxyz = pbuffer.data(idx_dip_gg + 484);

    auto tr_z_xxxz_xxzz = pbuffer.data(idx_dip_gg + 485);

    auto tr_z_xxxz_xyyy = pbuffer.data(idx_dip_gg + 486);

    auto tr_z_xxxz_xyyz = pbuffer.data(idx_dip_gg + 487);

    auto tr_z_xxxz_xyzz = pbuffer.data(idx_dip_gg + 488);

    auto tr_z_xxxz_xzzz = pbuffer.data(idx_dip_gg + 489);

    auto tr_z_xxxz_yyyy = pbuffer.data(idx_dip_gg + 490);

    auto tr_z_xxxz_yyyz = pbuffer.data(idx_dip_gg + 491);

    auto tr_z_xxxz_yyzz = pbuffer.data(idx_dip_gg + 492);

    auto tr_z_xxxz_yzzz = pbuffer.data(idx_dip_gg + 493);

    auto tr_z_xxxz_zzzz = pbuffer.data(idx_dip_gg + 494);

    auto tr_z_xxyy_xxxx = pbuffer.data(idx_dip_gg + 495);

    auto tr_z_xxyy_xxxy = pbuffer.data(idx_dip_gg + 496);

    auto tr_z_xxyy_xxxz = pbuffer.data(idx_dip_gg + 497);

    auto tr_z_xxyy_xxyy = pbuffer.data(idx_dip_gg + 498);

    auto tr_z_xxyy_xxyz = pbuffer.data(idx_dip_gg + 499);

    auto tr_z_xxyy_xxzz = pbuffer.data(idx_dip_gg + 500);

    auto tr_z_xxyy_xyyy = pbuffer.data(idx_dip_gg + 501);

    auto tr_z_xxyy_xyyz = pbuffer.data(idx_dip_gg + 502);

    auto tr_z_xxyy_xyzz = pbuffer.data(idx_dip_gg + 503);

    auto tr_z_xxyy_xzzz = pbuffer.data(idx_dip_gg + 504);

    auto tr_z_xxyy_yyyy = pbuffer.data(idx_dip_gg + 505);

    auto tr_z_xxyy_yyyz = pbuffer.data(idx_dip_gg + 506);

    auto tr_z_xxyy_yyzz = pbuffer.data(idx_dip_gg + 507);

    auto tr_z_xxyy_yzzz = pbuffer.data(idx_dip_gg + 508);

    auto tr_z_xxyy_zzzz = pbuffer.data(idx_dip_gg + 509);

    auto tr_z_xxyz_xxxx = pbuffer.data(idx_dip_gg + 510);

    auto tr_z_xxyz_xxxz = pbuffer.data(idx_dip_gg + 512);

    auto tr_z_xxyz_xxzz = pbuffer.data(idx_dip_gg + 515);

    auto tr_z_xxyz_xzzz = pbuffer.data(idx_dip_gg + 519);

    auto tr_z_xxyz_yyyy = pbuffer.data(idx_dip_gg + 520);

    auto tr_z_xxyz_yyyz = pbuffer.data(idx_dip_gg + 521);

    auto tr_z_xxyz_yyzz = pbuffer.data(idx_dip_gg + 522);

    auto tr_z_xxyz_yzzz = pbuffer.data(idx_dip_gg + 523);

    auto tr_z_xxzz_xxxx = pbuffer.data(idx_dip_gg + 525);

    auto tr_z_xxzz_xxxy = pbuffer.data(idx_dip_gg + 526);

    auto tr_z_xxzz_xxxz = pbuffer.data(idx_dip_gg + 527);

    auto tr_z_xxzz_xxyy = pbuffer.data(idx_dip_gg + 528);

    auto tr_z_xxzz_xxyz = pbuffer.data(idx_dip_gg + 529);

    auto tr_z_xxzz_xxzz = pbuffer.data(idx_dip_gg + 530);

    auto tr_z_xxzz_xyyy = pbuffer.data(idx_dip_gg + 531);

    auto tr_z_xxzz_xyyz = pbuffer.data(idx_dip_gg + 532);

    auto tr_z_xxzz_xyzz = pbuffer.data(idx_dip_gg + 533);

    auto tr_z_xxzz_xzzz = pbuffer.data(idx_dip_gg + 534);

    auto tr_z_xxzz_yyyy = pbuffer.data(idx_dip_gg + 535);

    auto tr_z_xxzz_yyyz = pbuffer.data(idx_dip_gg + 536);

    auto tr_z_xxzz_yyzz = pbuffer.data(idx_dip_gg + 537);

    auto tr_z_xxzz_yzzz = pbuffer.data(idx_dip_gg + 538);

    auto tr_z_xxzz_zzzz = pbuffer.data(idx_dip_gg + 539);

    auto tr_z_xyyy_xxxy = pbuffer.data(idx_dip_gg + 541);

    auto tr_z_xyyy_xxyy = pbuffer.data(idx_dip_gg + 543);

    auto tr_z_xyyy_xxyz = pbuffer.data(idx_dip_gg + 544);

    auto tr_z_xyyy_xyyy = pbuffer.data(idx_dip_gg + 546);

    auto tr_z_xyyy_xyyz = pbuffer.data(idx_dip_gg + 547);

    auto tr_z_xyyy_xyzz = pbuffer.data(idx_dip_gg + 548);

    auto tr_z_xyyy_yyyy = pbuffer.data(idx_dip_gg + 550);

    auto tr_z_xyyy_yyyz = pbuffer.data(idx_dip_gg + 551);

    auto tr_z_xyyy_yyzz = pbuffer.data(idx_dip_gg + 552);

    auto tr_z_xyyy_yzzz = pbuffer.data(idx_dip_gg + 553);

    auto tr_z_xyyy_zzzz = pbuffer.data(idx_dip_gg + 554);

    auto tr_z_xyyz_xxyz = pbuffer.data(idx_dip_gg + 559);

    auto tr_z_xyyz_xyyz = pbuffer.data(idx_dip_gg + 562);

    auto tr_z_xyyz_xyzz = pbuffer.data(idx_dip_gg + 563);

    auto tr_z_xyyz_yyyy = pbuffer.data(idx_dip_gg + 565);

    auto tr_z_xyyz_yyyz = pbuffer.data(idx_dip_gg + 566);

    auto tr_z_xyyz_yyzz = pbuffer.data(idx_dip_gg + 567);

    auto tr_z_xyyz_yzzz = pbuffer.data(idx_dip_gg + 568);

    auto tr_z_xyyz_zzzz = pbuffer.data(idx_dip_gg + 569);

    auto tr_z_xyzz_yyyy = pbuffer.data(idx_dip_gg + 580);

    auto tr_z_xyzz_yyyz = pbuffer.data(idx_dip_gg + 581);

    auto tr_z_xyzz_yyzz = pbuffer.data(idx_dip_gg + 582);

    auto tr_z_xyzz_yzzz = pbuffer.data(idx_dip_gg + 583);

    auto tr_z_xzzz_xxxx = pbuffer.data(idx_dip_gg + 585);

    auto tr_z_xzzz_xxxy = pbuffer.data(idx_dip_gg + 586);

    auto tr_z_xzzz_xxxz = pbuffer.data(idx_dip_gg + 587);

    auto tr_z_xzzz_xxyy = pbuffer.data(idx_dip_gg + 588);

    auto tr_z_xzzz_xxyz = pbuffer.data(idx_dip_gg + 589);

    auto tr_z_xzzz_xxzz = pbuffer.data(idx_dip_gg + 590);

    auto tr_z_xzzz_xyyy = pbuffer.data(idx_dip_gg + 591);

    auto tr_z_xzzz_xyyz = pbuffer.data(idx_dip_gg + 592);

    auto tr_z_xzzz_xyzz = pbuffer.data(idx_dip_gg + 593);

    auto tr_z_xzzz_xzzz = pbuffer.data(idx_dip_gg + 594);

    auto tr_z_xzzz_yyyy = pbuffer.data(idx_dip_gg + 595);

    auto tr_z_xzzz_yyyz = pbuffer.data(idx_dip_gg + 596);

    auto tr_z_xzzz_yyzz = pbuffer.data(idx_dip_gg + 597);

    auto tr_z_xzzz_yzzz = pbuffer.data(idx_dip_gg + 598);

    auto tr_z_xzzz_zzzz = pbuffer.data(idx_dip_gg + 599);

    auto tr_z_yyyy_xxxx = pbuffer.data(idx_dip_gg + 600);

    auto tr_z_yyyy_xxxy = pbuffer.data(idx_dip_gg + 601);

    auto tr_z_yyyy_xxxz = pbuffer.data(idx_dip_gg + 602);

    auto tr_z_yyyy_xxyy = pbuffer.data(idx_dip_gg + 603);

    auto tr_z_yyyy_xxyz = pbuffer.data(idx_dip_gg + 604);

    auto tr_z_yyyy_xxzz = pbuffer.data(idx_dip_gg + 605);

    auto tr_z_yyyy_xyyy = pbuffer.data(idx_dip_gg + 606);

    auto tr_z_yyyy_xyyz = pbuffer.data(idx_dip_gg + 607);

    auto tr_z_yyyy_xyzz = pbuffer.data(idx_dip_gg + 608);

    auto tr_z_yyyy_xzzz = pbuffer.data(idx_dip_gg + 609);

    auto tr_z_yyyy_yyyy = pbuffer.data(idx_dip_gg + 610);

    auto tr_z_yyyy_yyyz = pbuffer.data(idx_dip_gg + 611);

    auto tr_z_yyyy_yyzz = pbuffer.data(idx_dip_gg + 612);

    auto tr_z_yyyy_yzzz = pbuffer.data(idx_dip_gg + 613);

    auto tr_z_yyyy_zzzz = pbuffer.data(idx_dip_gg + 614);

    auto tr_z_yyyz_xxxx = pbuffer.data(idx_dip_gg + 615);

    auto tr_z_yyyz_xxxy = pbuffer.data(idx_dip_gg + 616);

    auto tr_z_yyyz_xxxz = pbuffer.data(idx_dip_gg + 617);

    auto tr_z_yyyz_xxyy = pbuffer.data(idx_dip_gg + 618);

    auto tr_z_yyyz_xxyz = pbuffer.data(idx_dip_gg + 619);

    auto tr_z_yyyz_xxzz = pbuffer.data(idx_dip_gg + 620);

    auto tr_z_yyyz_xyyy = pbuffer.data(idx_dip_gg + 621);

    auto tr_z_yyyz_xyyz = pbuffer.data(idx_dip_gg + 622);

    auto tr_z_yyyz_xyzz = pbuffer.data(idx_dip_gg + 623);

    auto tr_z_yyyz_xzzz = pbuffer.data(idx_dip_gg + 624);

    auto tr_z_yyyz_yyyy = pbuffer.data(idx_dip_gg + 625);

    auto tr_z_yyyz_yyyz = pbuffer.data(idx_dip_gg + 626);

    auto tr_z_yyyz_yyzz = pbuffer.data(idx_dip_gg + 627);

    auto tr_z_yyyz_yzzz = pbuffer.data(idx_dip_gg + 628);

    auto tr_z_yyyz_zzzz = pbuffer.data(idx_dip_gg + 629);

    auto tr_z_yyzz_xxxx = pbuffer.data(idx_dip_gg + 630);

    auto tr_z_yyzz_xxxy = pbuffer.data(idx_dip_gg + 631);

    auto tr_z_yyzz_xxxz = pbuffer.data(idx_dip_gg + 632);

    auto tr_z_yyzz_xxyy = pbuffer.data(idx_dip_gg + 633);

    auto tr_z_yyzz_xxyz = pbuffer.data(idx_dip_gg + 634);

    auto tr_z_yyzz_xxzz = pbuffer.data(idx_dip_gg + 635);

    auto tr_z_yyzz_xyyy = pbuffer.data(idx_dip_gg + 636);

    auto tr_z_yyzz_xyyz = pbuffer.data(idx_dip_gg + 637);

    auto tr_z_yyzz_xyzz = pbuffer.data(idx_dip_gg + 638);

    auto tr_z_yyzz_xzzz = pbuffer.data(idx_dip_gg + 639);

    auto tr_z_yyzz_yyyy = pbuffer.data(idx_dip_gg + 640);

    auto tr_z_yyzz_yyyz = pbuffer.data(idx_dip_gg + 641);

    auto tr_z_yyzz_yyzz = pbuffer.data(idx_dip_gg + 642);

    auto tr_z_yyzz_yzzz = pbuffer.data(idx_dip_gg + 643);

    auto tr_z_yyzz_zzzz = pbuffer.data(idx_dip_gg + 644);

    auto tr_z_yzzz_xxxx = pbuffer.data(idx_dip_gg + 645);

    auto tr_z_yzzz_xxxy = pbuffer.data(idx_dip_gg + 646);

    auto tr_z_yzzz_xxxz = pbuffer.data(idx_dip_gg + 647);

    auto tr_z_yzzz_xxyy = pbuffer.data(idx_dip_gg + 648);

    auto tr_z_yzzz_xxyz = pbuffer.data(idx_dip_gg + 649);

    auto tr_z_yzzz_xxzz = pbuffer.data(idx_dip_gg + 650);

    auto tr_z_yzzz_xyyy = pbuffer.data(idx_dip_gg + 651);

    auto tr_z_yzzz_xyyz = pbuffer.data(idx_dip_gg + 652);

    auto tr_z_yzzz_xyzz = pbuffer.data(idx_dip_gg + 653);

    auto tr_z_yzzz_xzzz = pbuffer.data(idx_dip_gg + 654);

    auto tr_z_yzzz_yyyy = pbuffer.data(idx_dip_gg + 655);

    auto tr_z_yzzz_yyyz = pbuffer.data(idx_dip_gg + 656);

    auto tr_z_yzzz_yyzz = pbuffer.data(idx_dip_gg + 657);

    auto tr_z_yzzz_yzzz = pbuffer.data(idx_dip_gg + 658);

    auto tr_z_yzzz_zzzz = pbuffer.data(idx_dip_gg + 659);

    auto tr_z_zzzz_xxxx = pbuffer.data(idx_dip_gg + 660);

    auto tr_z_zzzz_xxxy = pbuffer.data(idx_dip_gg + 661);

    auto tr_z_zzzz_xxxz = pbuffer.data(idx_dip_gg + 662);

    auto tr_z_zzzz_xxyy = pbuffer.data(idx_dip_gg + 663);

    auto tr_z_zzzz_xxyz = pbuffer.data(idx_dip_gg + 664);

    auto tr_z_zzzz_xxzz = pbuffer.data(idx_dip_gg + 665);

    auto tr_z_zzzz_xyyy = pbuffer.data(idx_dip_gg + 666);

    auto tr_z_zzzz_xyyz = pbuffer.data(idx_dip_gg + 667);

    auto tr_z_zzzz_xyzz = pbuffer.data(idx_dip_gg + 668);

    auto tr_z_zzzz_xzzz = pbuffer.data(idx_dip_gg + 669);

    auto tr_z_zzzz_yyyy = pbuffer.data(idx_dip_gg + 670);

    auto tr_z_zzzz_yyyz = pbuffer.data(idx_dip_gg + 671);

    auto tr_z_zzzz_yyzz = pbuffer.data(idx_dip_gg + 672);

    auto tr_z_zzzz_yzzz = pbuffer.data(idx_dip_gg + 673);

    auto tr_z_zzzz_zzzz = pbuffer.data(idx_dip_gg + 674);

    // Set up 0-15 components of targeted buffer : HG

    auto tr_x_xxxxx_xxxx = pbuffer.data(idx_dip_hg);

    auto tr_x_xxxxx_xxxy = pbuffer.data(idx_dip_hg + 1);

    auto tr_x_xxxxx_xxxz = pbuffer.data(idx_dip_hg + 2);

    auto tr_x_xxxxx_xxyy = pbuffer.data(idx_dip_hg + 3);

    auto tr_x_xxxxx_xxyz = pbuffer.data(idx_dip_hg + 4);

    auto tr_x_xxxxx_xxzz = pbuffer.data(idx_dip_hg + 5);

    auto tr_x_xxxxx_xyyy = pbuffer.data(idx_dip_hg + 6);

    auto tr_x_xxxxx_xyyz = pbuffer.data(idx_dip_hg + 7);

    auto tr_x_xxxxx_xyzz = pbuffer.data(idx_dip_hg + 8);

    auto tr_x_xxxxx_xzzz = pbuffer.data(idx_dip_hg + 9);

    auto tr_x_xxxxx_yyyy = pbuffer.data(idx_dip_hg + 10);

    auto tr_x_xxxxx_yyyz = pbuffer.data(idx_dip_hg + 11);

    auto tr_x_xxxxx_yyzz = pbuffer.data(idx_dip_hg + 12);

    auto tr_x_xxxxx_yzzz = pbuffer.data(idx_dip_hg + 13);

    auto tr_x_xxxxx_zzzz = pbuffer.data(idx_dip_hg + 14);

#pragma omp simd aligned(pa_x,                \
                             tr_x_xxx_xxxx,   \
                             tr_x_xxx_xxxy,   \
                             tr_x_xxx_xxxz,   \
                             tr_x_xxx_xxyy,   \
                             tr_x_xxx_xxyz,   \
                             tr_x_xxx_xxzz,   \
                             tr_x_xxx_xyyy,   \
                             tr_x_xxx_xyyz,   \
                             tr_x_xxx_xyzz,   \
                             tr_x_xxx_xzzz,   \
                             tr_x_xxx_yyyy,   \
                             tr_x_xxx_yyyz,   \
                             tr_x_xxx_yyzz,   \
                             tr_x_xxx_yzzz,   \
                             tr_x_xxx_zzzz,   \
                             tr_x_xxxx_xxx,   \
                             tr_x_xxxx_xxxx,  \
                             tr_x_xxxx_xxxy,  \
                             tr_x_xxxx_xxxz,  \
                             tr_x_xxxx_xxy,   \
                             tr_x_xxxx_xxyy,  \
                             tr_x_xxxx_xxyz,  \
                             tr_x_xxxx_xxz,   \
                             tr_x_xxxx_xxzz,  \
                             tr_x_xxxx_xyy,   \
                             tr_x_xxxx_xyyy,  \
                             tr_x_xxxx_xyyz,  \
                             tr_x_xxxx_xyz,   \
                             tr_x_xxxx_xyzz,  \
                             tr_x_xxxx_xzz,   \
                             tr_x_xxxx_xzzz,  \
                             tr_x_xxxx_yyy,   \
                             tr_x_xxxx_yyyy,  \
                             tr_x_xxxx_yyyz,  \
                             tr_x_xxxx_yyz,   \
                             tr_x_xxxx_yyzz,  \
                             tr_x_xxxx_yzz,   \
                             tr_x_xxxx_yzzz,  \
                             tr_x_xxxx_zzz,   \
                             tr_x_xxxx_zzzz,  \
                             tr_x_xxxxx_xxxx, \
                             tr_x_xxxxx_xxxy, \
                             tr_x_xxxxx_xxxz, \
                             tr_x_xxxxx_xxyy, \
                             tr_x_xxxxx_xxyz, \
                             tr_x_xxxxx_xxzz, \
                             tr_x_xxxxx_xyyy, \
                             tr_x_xxxxx_xyyz, \
                             tr_x_xxxxx_xyzz, \
                             tr_x_xxxxx_xzzz, \
                             tr_x_xxxxx_yyyy, \
                             tr_x_xxxxx_yyyz, \
                             tr_x_xxxxx_yyzz, \
                             tr_x_xxxxx_yzzz, \
                             tr_x_xxxxx_zzzz, \
                             ts_xxxx_xxxx,    \
                             ts_xxxx_xxxy,    \
                             ts_xxxx_xxxz,    \
                             ts_xxxx_xxyy,    \
                             ts_xxxx_xxyz,    \
                             ts_xxxx_xxzz,    \
                             ts_xxxx_xyyy,    \
                             ts_xxxx_xyyz,    \
                             ts_xxxx_xyzz,    \
                             ts_xxxx_xzzz,    \
                             ts_xxxx_yyyy,    \
                             ts_xxxx_yyyz,    \
                             ts_xxxx_yyzz,    \
                             ts_xxxx_yzzz,    \
                             ts_xxxx_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxx_xxxx[i] = 4.0 * tr_x_xxx_xxxx[i] * fe_0 + 4.0 * tr_x_xxxx_xxx[i] * fe_0 + ts_xxxx_xxxx[i] * fe_0 + tr_x_xxxx_xxxx[i] * pa_x[i];

        tr_x_xxxxx_xxxy[i] = 4.0 * tr_x_xxx_xxxy[i] * fe_0 + 3.0 * tr_x_xxxx_xxy[i] * fe_0 + ts_xxxx_xxxy[i] * fe_0 + tr_x_xxxx_xxxy[i] * pa_x[i];

        tr_x_xxxxx_xxxz[i] = 4.0 * tr_x_xxx_xxxz[i] * fe_0 + 3.0 * tr_x_xxxx_xxz[i] * fe_0 + ts_xxxx_xxxz[i] * fe_0 + tr_x_xxxx_xxxz[i] * pa_x[i];

        tr_x_xxxxx_xxyy[i] = 4.0 * tr_x_xxx_xxyy[i] * fe_0 + 2.0 * tr_x_xxxx_xyy[i] * fe_0 + ts_xxxx_xxyy[i] * fe_0 + tr_x_xxxx_xxyy[i] * pa_x[i];

        tr_x_xxxxx_xxyz[i] = 4.0 * tr_x_xxx_xxyz[i] * fe_0 + 2.0 * tr_x_xxxx_xyz[i] * fe_0 + ts_xxxx_xxyz[i] * fe_0 + tr_x_xxxx_xxyz[i] * pa_x[i];

        tr_x_xxxxx_xxzz[i] = 4.0 * tr_x_xxx_xxzz[i] * fe_0 + 2.0 * tr_x_xxxx_xzz[i] * fe_0 + ts_xxxx_xxzz[i] * fe_0 + tr_x_xxxx_xxzz[i] * pa_x[i];

        tr_x_xxxxx_xyyy[i] = 4.0 * tr_x_xxx_xyyy[i] * fe_0 + tr_x_xxxx_yyy[i] * fe_0 + ts_xxxx_xyyy[i] * fe_0 + tr_x_xxxx_xyyy[i] * pa_x[i];

        tr_x_xxxxx_xyyz[i] = 4.0 * tr_x_xxx_xyyz[i] * fe_0 + tr_x_xxxx_yyz[i] * fe_0 + ts_xxxx_xyyz[i] * fe_0 + tr_x_xxxx_xyyz[i] * pa_x[i];

        tr_x_xxxxx_xyzz[i] = 4.0 * tr_x_xxx_xyzz[i] * fe_0 + tr_x_xxxx_yzz[i] * fe_0 + ts_xxxx_xyzz[i] * fe_0 + tr_x_xxxx_xyzz[i] * pa_x[i];

        tr_x_xxxxx_xzzz[i] = 4.0 * tr_x_xxx_xzzz[i] * fe_0 + tr_x_xxxx_zzz[i] * fe_0 + ts_xxxx_xzzz[i] * fe_0 + tr_x_xxxx_xzzz[i] * pa_x[i];

        tr_x_xxxxx_yyyy[i] = 4.0 * tr_x_xxx_yyyy[i] * fe_0 + ts_xxxx_yyyy[i] * fe_0 + tr_x_xxxx_yyyy[i] * pa_x[i];

        tr_x_xxxxx_yyyz[i] = 4.0 * tr_x_xxx_yyyz[i] * fe_0 + ts_xxxx_yyyz[i] * fe_0 + tr_x_xxxx_yyyz[i] * pa_x[i];

        tr_x_xxxxx_yyzz[i] = 4.0 * tr_x_xxx_yyzz[i] * fe_0 + ts_xxxx_yyzz[i] * fe_0 + tr_x_xxxx_yyzz[i] * pa_x[i];

        tr_x_xxxxx_yzzz[i] = 4.0 * tr_x_xxx_yzzz[i] * fe_0 + ts_xxxx_yzzz[i] * fe_0 + tr_x_xxxx_yzzz[i] * pa_x[i];

        tr_x_xxxxx_zzzz[i] = 4.0 * tr_x_xxx_zzzz[i] * fe_0 + ts_xxxx_zzzz[i] * fe_0 + tr_x_xxxx_zzzz[i] * pa_x[i];
    }

    // Set up 15-30 components of targeted buffer : HG

    auto tr_x_xxxxy_xxxx = pbuffer.data(idx_dip_hg + 15);

    auto tr_x_xxxxy_xxxy = pbuffer.data(idx_dip_hg + 16);

    auto tr_x_xxxxy_xxxz = pbuffer.data(idx_dip_hg + 17);

    auto tr_x_xxxxy_xxyy = pbuffer.data(idx_dip_hg + 18);

    auto tr_x_xxxxy_xxyz = pbuffer.data(idx_dip_hg + 19);

    auto tr_x_xxxxy_xxzz = pbuffer.data(idx_dip_hg + 20);

    auto tr_x_xxxxy_xyyy = pbuffer.data(idx_dip_hg + 21);

    auto tr_x_xxxxy_xyyz = pbuffer.data(idx_dip_hg + 22);

    auto tr_x_xxxxy_xyzz = pbuffer.data(idx_dip_hg + 23);

    auto tr_x_xxxxy_xzzz = pbuffer.data(idx_dip_hg + 24);

    auto tr_x_xxxxy_yyyy = pbuffer.data(idx_dip_hg + 25);

    auto tr_x_xxxxy_yyyz = pbuffer.data(idx_dip_hg + 26);

    auto tr_x_xxxxy_yyzz = pbuffer.data(idx_dip_hg + 27);

    auto tr_x_xxxxy_yzzz = pbuffer.data(idx_dip_hg + 28);

    auto tr_x_xxxxy_zzzz = pbuffer.data(idx_dip_hg + 29);

#pragma omp simd aligned(pa_y,                \
                             tr_x_xxxx_xxx,   \
                             tr_x_xxxx_xxxx,  \
                             tr_x_xxxx_xxxy,  \
                             tr_x_xxxx_xxxz,  \
                             tr_x_xxxx_xxy,   \
                             tr_x_xxxx_xxyy,  \
                             tr_x_xxxx_xxyz,  \
                             tr_x_xxxx_xxz,   \
                             tr_x_xxxx_xxzz,  \
                             tr_x_xxxx_xyy,   \
                             tr_x_xxxx_xyyy,  \
                             tr_x_xxxx_xyyz,  \
                             tr_x_xxxx_xyz,   \
                             tr_x_xxxx_xyzz,  \
                             tr_x_xxxx_xzz,   \
                             tr_x_xxxx_xzzz,  \
                             tr_x_xxxx_yyy,   \
                             tr_x_xxxx_yyyy,  \
                             tr_x_xxxx_yyyz,  \
                             tr_x_xxxx_yyz,   \
                             tr_x_xxxx_yyzz,  \
                             tr_x_xxxx_yzz,   \
                             tr_x_xxxx_yzzz,  \
                             tr_x_xxxx_zzz,   \
                             tr_x_xxxx_zzzz,  \
                             tr_x_xxxxy_xxxx, \
                             tr_x_xxxxy_xxxy, \
                             tr_x_xxxxy_xxxz, \
                             tr_x_xxxxy_xxyy, \
                             tr_x_xxxxy_xxyz, \
                             tr_x_xxxxy_xxzz, \
                             tr_x_xxxxy_xyyy, \
                             tr_x_xxxxy_xyyz, \
                             tr_x_xxxxy_xyzz, \
                             tr_x_xxxxy_xzzz, \
                             tr_x_xxxxy_yyyy, \
                             tr_x_xxxxy_yyyz, \
                             tr_x_xxxxy_yyzz, \
                             tr_x_xxxxy_yzzz, \
                             tr_x_xxxxy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxy_xxxx[i] = tr_x_xxxx_xxxx[i] * pa_y[i];

        tr_x_xxxxy_xxxy[i] = tr_x_xxxx_xxx[i] * fe_0 + tr_x_xxxx_xxxy[i] * pa_y[i];

        tr_x_xxxxy_xxxz[i] = tr_x_xxxx_xxxz[i] * pa_y[i];

        tr_x_xxxxy_xxyy[i] = 2.0 * tr_x_xxxx_xxy[i] * fe_0 + tr_x_xxxx_xxyy[i] * pa_y[i];

        tr_x_xxxxy_xxyz[i] = tr_x_xxxx_xxz[i] * fe_0 + tr_x_xxxx_xxyz[i] * pa_y[i];

        tr_x_xxxxy_xxzz[i] = tr_x_xxxx_xxzz[i] * pa_y[i];

        tr_x_xxxxy_xyyy[i] = 3.0 * tr_x_xxxx_xyy[i] * fe_0 + tr_x_xxxx_xyyy[i] * pa_y[i];

        tr_x_xxxxy_xyyz[i] = 2.0 * tr_x_xxxx_xyz[i] * fe_0 + tr_x_xxxx_xyyz[i] * pa_y[i];

        tr_x_xxxxy_xyzz[i] = tr_x_xxxx_xzz[i] * fe_0 + tr_x_xxxx_xyzz[i] * pa_y[i];

        tr_x_xxxxy_xzzz[i] = tr_x_xxxx_xzzz[i] * pa_y[i];

        tr_x_xxxxy_yyyy[i] = 4.0 * tr_x_xxxx_yyy[i] * fe_0 + tr_x_xxxx_yyyy[i] * pa_y[i];

        tr_x_xxxxy_yyyz[i] = 3.0 * tr_x_xxxx_yyz[i] * fe_0 + tr_x_xxxx_yyyz[i] * pa_y[i];

        tr_x_xxxxy_yyzz[i] = 2.0 * tr_x_xxxx_yzz[i] * fe_0 + tr_x_xxxx_yyzz[i] * pa_y[i];

        tr_x_xxxxy_yzzz[i] = tr_x_xxxx_zzz[i] * fe_0 + tr_x_xxxx_yzzz[i] * pa_y[i];

        tr_x_xxxxy_zzzz[i] = tr_x_xxxx_zzzz[i] * pa_y[i];
    }

    // Set up 30-45 components of targeted buffer : HG

    auto tr_x_xxxxz_xxxx = pbuffer.data(idx_dip_hg + 30);

    auto tr_x_xxxxz_xxxy = pbuffer.data(idx_dip_hg + 31);

    auto tr_x_xxxxz_xxxz = pbuffer.data(idx_dip_hg + 32);

    auto tr_x_xxxxz_xxyy = pbuffer.data(idx_dip_hg + 33);

    auto tr_x_xxxxz_xxyz = pbuffer.data(idx_dip_hg + 34);

    auto tr_x_xxxxz_xxzz = pbuffer.data(idx_dip_hg + 35);

    auto tr_x_xxxxz_xyyy = pbuffer.data(idx_dip_hg + 36);

    auto tr_x_xxxxz_xyyz = pbuffer.data(idx_dip_hg + 37);

    auto tr_x_xxxxz_xyzz = pbuffer.data(idx_dip_hg + 38);

    auto tr_x_xxxxz_xzzz = pbuffer.data(idx_dip_hg + 39);

    auto tr_x_xxxxz_yyyy = pbuffer.data(idx_dip_hg + 40);

    auto tr_x_xxxxz_yyyz = pbuffer.data(idx_dip_hg + 41);

    auto tr_x_xxxxz_yyzz = pbuffer.data(idx_dip_hg + 42);

    auto tr_x_xxxxz_yzzz = pbuffer.data(idx_dip_hg + 43);

    auto tr_x_xxxxz_zzzz = pbuffer.data(idx_dip_hg + 44);

#pragma omp simd aligned(pa_z,                \
                             tr_x_xxxx_xxx,   \
                             tr_x_xxxx_xxxx,  \
                             tr_x_xxxx_xxxy,  \
                             tr_x_xxxx_xxxz,  \
                             tr_x_xxxx_xxy,   \
                             tr_x_xxxx_xxyy,  \
                             tr_x_xxxx_xxyz,  \
                             tr_x_xxxx_xxz,   \
                             tr_x_xxxx_xxzz,  \
                             tr_x_xxxx_xyy,   \
                             tr_x_xxxx_xyyy,  \
                             tr_x_xxxx_xyyz,  \
                             tr_x_xxxx_xyz,   \
                             tr_x_xxxx_xyzz,  \
                             tr_x_xxxx_xzz,   \
                             tr_x_xxxx_xzzz,  \
                             tr_x_xxxx_yyy,   \
                             tr_x_xxxx_yyyy,  \
                             tr_x_xxxx_yyyz,  \
                             tr_x_xxxx_yyz,   \
                             tr_x_xxxx_yyzz,  \
                             tr_x_xxxx_yzz,   \
                             tr_x_xxxx_yzzz,  \
                             tr_x_xxxx_zzz,   \
                             tr_x_xxxx_zzzz,  \
                             tr_x_xxxxz_xxxx, \
                             tr_x_xxxxz_xxxy, \
                             tr_x_xxxxz_xxxz, \
                             tr_x_xxxxz_xxyy, \
                             tr_x_xxxxz_xxyz, \
                             tr_x_xxxxz_xxzz, \
                             tr_x_xxxxz_xyyy, \
                             tr_x_xxxxz_xyyz, \
                             tr_x_xxxxz_xyzz, \
                             tr_x_xxxxz_xzzz, \
                             tr_x_xxxxz_yyyy, \
                             tr_x_xxxxz_yyyz, \
                             tr_x_xxxxz_yyzz, \
                             tr_x_xxxxz_yzzz, \
                             tr_x_xxxxz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxz_xxxx[i] = tr_x_xxxx_xxxx[i] * pa_z[i];

        tr_x_xxxxz_xxxy[i] = tr_x_xxxx_xxxy[i] * pa_z[i];

        tr_x_xxxxz_xxxz[i] = tr_x_xxxx_xxx[i] * fe_0 + tr_x_xxxx_xxxz[i] * pa_z[i];

        tr_x_xxxxz_xxyy[i] = tr_x_xxxx_xxyy[i] * pa_z[i];

        tr_x_xxxxz_xxyz[i] = tr_x_xxxx_xxy[i] * fe_0 + tr_x_xxxx_xxyz[i] * pa_z[i];

        tr_x_xxxxz_xxzz[i] = 2.0 * tr_x_xxxx_xxz[i] * fe_0 + tr_x_xxxx_xxzz[i] * pa_z[i];

        tr_x_xxxxz_xyyy[i] = tr_x_xxxx_xyyy[i] * pa_z[i];

        tr_x_xxxxz_xyyz[i] = tr_x_xxxx_xyy[i] * fe_0 + tr_x_xxxx_xyyz[i] * pa_z[i];

        tr_x_xxxxz_xyzz[i] = 2.0 * tr_x_xxxx_xyz[i] * fe_0 + tr_x_xxxx_xyzz[i] * pa_z[i];

        tr_x_xxxxz_xzzz[i] = 3.0 * tr_x_xxxx_xzz[i] * fe_0 + tr_x_xxxx_xzzz[i] * pa_z[i];

        tr_x_xxxxz_yyyy[i] = tr_x_xxxx_yyyy[i] * pa_z[i];

        tr_x_xxxxz_yyyz[i] = tr_x_xxxx_yyy[i] * fe_0 + tr_x_xxxx_yyyz[i] * pa_z[i];

        tr_x_xxxxz_yyzz[i] = 2.0 * tr_x_xxxx_yyz[i] * fe_0 + tr_x_xxxx_yyzz[i] * pa_z[i];

        tr_x_xxxxz_yzzz[i] = 3.0 * tr_x_xxxx_yzz[i] * fe_0 + tr_x_xxxx_yzzz[i] * pa_z[i];

        tr_x_xxxxz_zzzz[i] = 4.0 * tr_x_xxxx_zzz[i] * fe_0 + tr_x_xxxx_zzzz[i] * pa_z[i];
    }

    // Set up 45-60 components of targeted buffer : HG

    auto tr_x_xxxyy_xxxx = pbuffer.data(idx_dip_hg + 45);

    auto tr_x_xxxyy_xxxy = pbuffer.data(idx_dip_hg + 46);

    auto tr_x_xxxyy_xxxz = pbuffer.data(idx_dip_hg + 47);

    auto tr_x_xxxyy_xxyy = pbuffer.data(idx_dip_hg + 48);

    auto tr_x_xxxyy_xxyz = pbuffer.data(idx_dip_hg + 49);

    auto tr_x_xxxyy_xxzz = pbuffer.data(idx_dip_hg + 50);

    auto tr_x_xxxyy_xyyy = pbuffer.data(idx_dip_hg + 51);

    auto tr_x_xxxyy_xyyz = pbuffer.data(idx_dip_hg + 52);

    auto tr_x_xxxyy_xyzz = pbuffer.data(idx_dip_hg + 53);

    auto tr_x_xxxyy_xzzz = pbuffer.data(idx_dip_hg + 54);

    auto tr_x_xxxyy_yyyy = pbuffer.data(idx_dip_hg + 55);

    auto tr_x_xxxyy_yyyz = pbuffer.data(idx_dip_hg + 56);

    auto tr_x_xxxyy_yyzz = pbuffer.data(idx_dip_hg + 57);

    auto tr_x_xxxyy_yzzz = pbuffer.data(idx_dip_hg + 58);

    auto tr_x_xxxyy_zzzz = pbuffer.data(idx_dip_hg + 59);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tr_x_xxx_xxxx,   \
                             tr_x_xxx_xxxy,   \
                             tr_x_xxx_xxxz,   \
                             tr_x_xxx_xxyy,   \
                             tr_x_xxx_xxyz,   \
                             tr_x_xxx_xxzz,   \
                             tr_x_xxx_xyyy,   \
                             tr_x_xxx_xyyz,   \
                             tr_x_xxx_xyzz,   \
                             tr_x_xxx_xzzz,   \
                             tr_x_xxx_zzzz,   \
                             tr_x_xxxy_xxx,   \
                             tr_x_xxxy_xxxx,  \
                             tr_x_xxxy_xxxy,  \
                             tr_x_xxxy_xxxz,  \
                             tr_x_xxxy_xxy,   \
                             tr_x_xxxy_xxyy,  \
                             tr_x_xxxy_xxyz,  \
                             tr_x_xxxy_xxz,   \
                             tr_x_xxxy_xxzz,  \
                             tr_x_xxxy_xyy,   \
                             tr_x_xxxy_xyyy,  \
                             tr_x_xxxy_xyyz,  \
                             tr_x_xxxy_xyz,   \
                             tr_x_xxxy_xyzz,  \
                             tr_x_xxxy_xzz,   \
                             tr_x_xxxy_xzzz,  \
                             tr_x_xxxy_zzzz,  \
                             tr_x_xxxyy_xxxx, \
                             tr_x_xxxyy_xxxy, \
                             tr_x_xxxyy_xxxz, \
                             tr_x_xxxyy_xxyy, \
                             tr_x_xxxyy_xxyz, \
                             tr_x_xxxyy_xxzz, \
                             tr_x_xxxyy_xyyy, \
                             tr_x_xxxyy_xyyz, \
                             tr_x_xxxyy_xyzz, \
                             tr_x_xxxyy_xzzz, \
                             tr_x_xxxyy_yyyy, \
                             tr_x_xxxyy_yyyz, \
                             tr_x_xxxyy_yyzz, \
                             tr_x_xxxyy_yzzz, \
                             tr_x_xxxyy_zzzz, \
                             tr_x_xxyy_yyyy,  \
                             tr_x_xxyy_yyyz,  \
                             tr_x_xxyy_yyzz,  \
                             tr_x_xxyy_yzzz,  \
                             tr_x_xyy_yyyy,   \
                             tr_x_xyy_yyyz,   \
                             tr_x_xyy_yyzz,   \
                             tr_x_xyy_yzzz,   \
                             ts_xxyy_yyyy,    \
                             ts_xxyy_yyyz,    \
                             ts_xxyy_yyzz,    \
                             ts_xxyy_yzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyy_xxxx[i] = tr_x_xxx_xxxx[i] * fe_0 + tr_x_xxxy_xxxx[i] * pa_y[i];

        tr_x_xxxyy_xxxy[i] = tr_x_xxx_xxxy[i] * fe_0 + tr_x_xxxy_xxx[i] * fe_0 + tr_x_xxxy_xxxy[i] * pa_y[i];

        tr_x_xxxyy_xxxz[i] = tr_x_xxx_xxxz[i] * fe_0 + tr_x_xxxy_xxxz[i] * pa_y[i];

        tr_x_xxxyy_xxyy[i] = tr_x_xxx_xxyy[i] * fe_0 + 2.0 * tr_x_xxxy_xxy[i] * fe_0 + tr_x_xxxy_xxyy[i] * pa_y[i];

        tr_x_xxxyy_xxyz[i] = tr_x_xxx_xxyz[i] * fe_0 + tr_x_xxxy_xxz[i] * fe_0 + tr_x_xxxy_xxyz[i] * pa_y[i];

        tr_x_xxxyy_xxzz[i] = tr_x_xxx_xxzz[i] * fe_0 + tr_x_xxxy_xxzz[i] * pa_y[i];

        tr_x_xxxyy_xyyy[i] = tr_x_xxx_xyyy[i] * fe_0 + 3.0 * tr_x_xxxy_xyy[i] * fe_0 + tr_x_xxxy_xyyy[i] * pa_y[i];

        tr_x_xxxyy_xyyz[i] = tr_x_xxx_xyyz[i] * fe_0 + 2.0 * tr_x_xxxy_xyz[i] * fe_0 + tr_x_xxxy_xyyz[i] * pa_y[i];

        tr_x_xxxyy_xyzz[i] = tr_x_xxx_xyzz[i] * fe_0 + tr_x_xxxy_xzz[i] * fe_0 + tr_x_xxxy_xyzz[i] * pa_y[i];

        tr_x_xxxyy_xzzz[i] = tr_x_xxx_xzzz[i] * fe_0 + tr_x_xxxy_xzzz[i] * pa_y[i];

        tr_x_xxxyy_yyyy[i] = 2.0 * tr_x_xyy_yyyy[i] * fe_0 + ts_xxyy_yyyy[i] * fe_0 + tr_x_xxyy_yyyy[i] * pa_x[i];

        tr_x_xxxyy_yyyz[i] = 2.0 * tr_x_xyy_yyyz[i] * fe_0 + ts_xxyy_yyyz[i] * fe_0 + tr_x_xxyy_yyyz[i] * pa_x[i];

        tr_x_xxxyy_yyzz[i] = 2.0 * tr_x_xyy_yyzz[i] * fe_0 + ts_xxyy_yyzz[i] * fe_0 + tr_x_xxyy_yyzz[i] * pa_x[i];

        tr_x_xxxyy_yzzz[i] = 2.0 * tr_x_xyy_yzzz[i] * fe_0 + ts_xxyy_yzzz[i] * fe_0 + tr_x_xxyy_yzzz[i] * pa_x[i];

        tr_x_xxxyy_zzzz[i] = tr_x_xxx_zzzz[i] * fe_0 + tr_x_xxxy_zzzz[i] * pa_y[i];
    }

    // Set up 60-75 components of targeted buffer : HG

    auto tr_x_xxxyz_xxxx = pbuffer.data(idx_dip_hg + 60);

    auto tr_x_xxxyz_xxxy = pbuffer.data(idx_dip_hg + 61);

    auto tr_x_xxxyz_xxxz = pbuffer.data(idx_dip_hg + 62);

    auto tr_x_xxxyz_xxyy = pbuffer.data(idx_dip_hg + 63);

    auto tr_x_xxxyz_xxyz = pbuffer.data(idx_dip_hg + 64);

    auto tr_x_xxxyz_xxzz = pbuffer.data(idx_dip_hg + 65);

    auto tr_x_xxxyz_xyyy = pbuffer.data(idx_dip_hg + 66);

    auto tr_x_xxxyz_xyyz = pbuffer.data(idx_dip_hg + 67);

    auto tr_x_xxxyz_xyzz = pbuffer.data(idx_dip_hg + 68);

    auto tr_x_xxxyz_xzzz = pbuffer.data(idx_dip_hg + 69);

    auto tr_x_xxxyz_yyyy = pbuffer.data(idx_dip_hg + 70);

    auto tr_x_xxxyz_yyyz = pbuffer.data(idx_dip_hg + 71);

    auto tr_x_xxxyz_yyzz = pbuffer.data(idx_dip_hg + 72);

    auto tr_x_xxxyz_yzzz = pbuffer.data(idx_dip_hg + 73);

    auto tr_x_xxxyz_zzzz = pbuffer.data(idx_dip_hg + 74);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             tr_x_xxxy_xxxy,  \
                             tr_x_xxxy_xxyy,  \
                             tr_x_xxxy_xyyy,  \
                             tr_x_xxxy_yyyy,  \
                             tr_x_xxxyz_xxxx, \
                             tr_x_xxxyz_xxxy, \
                             tr_x_xxxyz_xxxz, \
                             tr_x_xxxyz_xxyy, \
                             tr_x_xxxyz_xxyz, \
                             tr_x_xxxyz_xxzz, \
                             tr_x_xxxyz_xyyy, \
                             tr_x_xxxyz_xyyz, \
                             tr_x_xxxyz_xyzz, \
                             tr_x_xxxyz_xzzz, \
                             tr_x_xxxyz_yyyy, \
                             tr_x_xxxyz_yyyz, \
                             tr_x_xxxyz_yyzz, \
                             tr_x_xxxyz_yzzz, \
                             tr_x_xxxyz_zzzz, \
                             tr_x_xxxz_xxxx,  \
                             tr_x_xxxz_xxxz,  \
                             tr_x_xxxz_xxyz,  \
                             tr_x_xxxz_xxz,   \
                             tr_x_xxxz_xxzz,  \
                             tr_x_xxxz_xyyz,  \
                             tr_x_xxxz_xyz,   \
                             tr_x_xxxz_xyzz,  \
                             tr_x_xxxz_xzz,   \
                             tr_x_xxxz_xzzz,  \
                             tr_x_xxxz_yyyz,  \
                             tr_x_xxxz_yyz,   \
                             tr_x_xxxz_yyzz,  \
                             tr_x_xxxz_yzz,   \
                             tr_x_xxxz_yzzz,  \
                             tr_x_xxxz_zzz,   \
                             tr_x_xxxz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyz_xxxx[i] = tr_x_xxxz_xxxx[i] * pa_y[i];

        tr_x_xxxyz_xxxy[i] = tr_x_xxxy_xxxy[i] * pa_z[i];

        tr_x_xxxyz_xxxz[i] = tr_x_xxxz_xxxz[i] * pa_y[i];

        tr_x_xxxyz_xxyy[i] = tr_x_xxxy_xxyy[i] * pa_z[i];

        tr_x_xxxyz_xxyz[i] = tr_x_xxxz_xxz[i] * fe_0 + tr_x_xxxz_xxyz[i] * pa_y[i];

        tr_x_xxxyz_xxzz[i] = tr_x_xxxz_xxzz[i] * pa_y[i];

        tr_x_xxxyz_xyyy[i] = tr_x_xxxy_xyyy[i] * pa_z[i];

        tr_x_xxxyz_xyyz[i] = 2.0 * tr_x_xxxz_xyz[i] * fe_0 + tr_x_xxxz_xyyz[i] * pa_y[i];

        tr_x_xxxyz_xyzz[i] = tr_x_xxxz_xzz[i] * fe_0 + tr_x_xxxz_xyzz[i] * pa_y[i];

        tr_x_xxxyz_xzzz[i] = tr_x_xxxz_xzzz[i] * pa_y[i];

        tr_x_xxxyz_yyyy[i] = tr_x_xxxy_yyyy[i] * pa_z[i];

        tr_x_xxxyz_yyyz[i] = 3.0 * tr_x_xxxz_yyz[i] * fe_0 + tr_x_xxxz_yyyz[i] * pa_y[i];

        tr_x_xxxyz_yyzz[i] = 2.0 * tr_x_xxxz_yzz[i] * fe_0 + tr_x_xxxz_yyzz[i] * pa_y[i];

        tr_x_xxxyz_yzzz[i] = tr_x_xxxz_zzz[i] * fe_0 + tr_x_xxxz_yzzz[i] * pa_y[i];

        tr_x_xxxyz_zzzz[i] = tr_x_xxxz_zzzz[i] * pa_y[i];
    }

    // Set up 75-90 components of targeted buffer : HG

    auto tr_x_xxxzz_xxxx = pbuffer.data(idx_dip_hg + 75);

    auto tr_x_xxxzz_xxxy = pbuffer.data(idx_dip_hg + 76);

    auto tr_x_xxxzz_xxxz = pbuffer.data(idx_dip_hg + 77);

    auto tr_x_xxxzz_xxyy = pbuffer.data(idx_dip_hg + 78);

    auto tr_x_xxxzz_xxyz = pbuffer.data(idx_dip_hg + 79);

    auto tr_x_xxxzz_xxzz = pbuffer.data(idx_dip_hg + 80);

    auto tr_x_xxxzz_xyyy = pbuffer.data(idx_dip_hg + 81);

    auto tr_x_xxxzz_xyyz = pbuffer.data(idx_dip_hg + 82);

    auto tr_x_xxxzz_xyzz = pbuffer.data(idx_dip_hg + 83);

    auto tr_x_xxxzz_xzzz = pbuffer.data(idx_dip_hg + 84);

    auto tr_x_xxxzz_yyyy = pbuffer.data(idx_dip_hg + 85);

    auto tr_x_xxxzz_yyyz = pbuffer.data(idx_dip_hg + 86);

    auto tr_x_xxxzz_yyzz = pbuffer.data(idx_dip_hg + 87);

    auto tr_x_xxxzz_yzzz = pbuffer.data(idx_dip_hg + 88);

    auto tr_x_xxxzz_zzzz = pbuffer.data(idx_dip_hg + 89);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             tr_x_xxx_xxxx,   \
                             tr_x_xxx_xxxy,   \
                             tr_x_xxx_xxxz,   \
                             tr_x_xxx_xxyy,   \
                             tr_x_xxx_xxyz,   \
                             tr_x_xxx_xxzz,   \
                             tr_x_xxx_xyyy,   \
                             tr_x_xxx_xyyz,   \
                             tr_x_xxx_xyzz,   \
                             tr_x_xxx_xzzz,   \
                             tr_x_xxx_yyyy,   \
                             tr_x_xxxz_xxx,   \
                             tr_x_xxxz_xxxx,  \
                             tr_x_xxxz_xxxy,  \
                             tr_x_xxxz_xxxz,  \
                             tr_x_xxxz_xxy,   \
                             tr_x_xxxz_xxyy,  \
                             tr_x_xxxz_xxyz,  \
                             tr_x_xxxz_xxz,   \
                             tr_x_xxxz_xxzz,  \
                             tr_x_xxxz_xyy,   \
                             tr_x_xxxz_xyyy,  \
                             tr_x_xxxz_xyyz,  \
                             tr_x_xxxz_xyz,   \
                             tr_x_xxxz_xyzz,  \
                             tr_x_xxxz_xzz,   \
                             tr_x_xxxz_xzzz,  \
                             tr_x_xxxz_yyyy,  \
                             tr_x_xxxzz_xxxx, \
                             tr_x_xxxzz_xxxy, \
                             tr_x_xxxzz_xxxz, \
                             tr_x_xxxzz_xxyy, \
                             tr_x_xxxzz_xxyz, \
                             tr_x_xxxzz_xxzz, \
                             tr_x_xxxzz_xyyy, \
                             tr_x_xxxzz_xyyz, \
                             tr_x_xxxzz_xyzz, \
                             tr_x_xxxzz_xzzz, \
                             tr_x_xxxzz_yyyy, \
                             tr_x_xxxzz_yyyz, \
                             tr_x_xxxzz_yyzz, \
                             tr_x_xxxzz_yzzz, \
                             tr_x_xxxzz_zzzz, \
                             tr_x_xxzz_yyyz,  \
                             tr_x_xxzz_yyzz,  \
                             tr_x_xxzz_yzzz,  \
                             tr_x_xxzz_zzzz,  \
                             tr_x_xzz_yyyz,   \
                             tr_x_xzz_yyzz,   \
                             tr_x_xzz_yzzz,   \
                             tr_x_xzz_zzzz,   \
                             ts_xxzz_yyyz,    \
                             ts_xxzz_yyzz,    \
                             ts_xxzz_yzzz,    \
                             ts_xxzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxzz_xxxx[i] = tr_x_xxx_xxxx[i] * fe_0 + tr_x_xxxz_xxxx[i] * pa_z[i];

        tr_x_xxxzz_xxxy[i] = tr_x_xxx_xxxy[i] * fe_0 + tr_x_xxxz_xxxy[i] * pa_z[i];

        tr_x_xxxzz_xxxz[i] = tr_x_xxx_xxxz[i] * fe_0 + tr_x_xxxz_xxx[i] * fe_0 + tr_x_xxxz_xxxz[i] * pa_z[i];

        tr_x_xxxzz_xxyy[i] = tr_x_xxx_xxyy[i] * fe_0 + tr_x_xxxz_xxyy[i] * pa_z[i];

        tr_x_xxxzz_xxyz[i] = tr_x_xxx_xxyz[i] * fe_0 + tr_x_xxxz_xxy[i] * fe_0 + tr_x_xxxz_xxyz[i] * pa_z[i];

        tr_x_xxxzz_xxzz[i] = tr_x_xxx_xxzz[i] * fe_0 + 2.0 * tr_x_xxxz_xxz[i] * fe_0 + tr_x_xxxz_xxzz[i] * pa_z[i];

        tr_x_xxxzz_xyyy[i] = tr_x_xxx_xyyy[i] * fe_0 + tr_x_xxxz_xyyy[i] * pa_z[i];

        tr_x_xxxzz_xyyz[i] = tr_x_xxx_xyyz[i] * fe_0 + tr_x_xxxz_xyy[i] * fe_0 + tr_x_xxxz_xyyz[i] * pa_z[i];

        tr_x_xxxzz_xyzz[i] = tr_x_xxx_xyzz[i] * fe_0 + 2.0 * tr_x_xxxz_xyz[i] * fe_0 + tr_x_xxxz_xyzz[i] * pa_z[i];

        tr_x_xxxzz_xzzz[i] = tr_x_xxx_xzzz[i] * fe_0 + 3.0 * tr_x_xxxz_xzz[i] * fe_0 + tr_x_xxxz_xzzz[i] * pa_z[i];

        tr_x_xxxzz_yyyy[i] = tr_x_xxx_yyyy[i] * fe_0 + tr_x_xxxz_yyyy[i] * pa_z[i];

        tr_x_xxxzz_yyyz[i] = 2.0 * tr_x_xzz_yyyz[i] * fe_0 + ts_xxzz_yyyz[i] * fe_0 + tr_x_xxzz_yyyz[i] * pa_x[i];

        tr_x_xxxzz_yyzz[i] = 2.0 * tr_x_xzz_yyzz[i] * fe_0 + ts_xxzz_yyzz[i] * fe_0 + tr_x_xxzz_yyzz[i] * pa_x[i];

        tr_x_xxxzz_yzzz[i] = 2.0 * tr_x_xzz_yzzz[i] * fe_0 + ts_xxzz_yzzz[i] * fe_0 + tr_x_xxzz_yzzz[i] * pa_x[i];

        tr_x_xxxzz_zzzz[i] = 2.0 * tr_x_xzz_zzzz[i] * fe_0 + ts_xxzz_zzzz[i] * fe_0 + tr_x_xxzz_zzzz[i] * pa_x[i];
    }

    // Set up 90-105 components of targeted buffer : HG

    auto tr_x_xxyyy_xxxx = pbuffer.data(idx_dip_hg + 90);

    auto tr_x_xxyyy_xxxy = pbuffer.data(idx_dip_hg + 91);

    auto tr_x_xxyyy_xxxz = pbuffer.data(idx_dip_hg + 92);

    auto tr_x_xxyyy_xxyy = pbuffer.data(idx_dip_hg + 93);

    auto tr_x_xxyyy_xxyz = pbuffer.data(idx_dip_hg + 94);

    auto tr_x_xxyyy_xxzz = pbuffer.data(idx_dip_hg + 95);

    auto tr_x_xxyyy_xyyy = pbuffer.data(idx_dip_hg + 96);

    auto tr_x_xxyyy_xyyz = pbuffer.data(idx_dip_hg + 97);

    auto tr_x_xxyyy_xyzz = pbuffer.data(idx_dip_hg + 98);

    auto tr_x_xxyyy_xzzz = pbuffer.data(idx_dip_hg + 99);

    auto tr_x_xxyyy_yyyy = pbuffer.data(idx_dip_hg + 100);

    auto tr_x_xxyyy_yyyz = pbuffer.data(idx_dip_hg + 101);

    auto tr_x_xxyyy_yyzz = pbuffer.data(idx_dip_hg + 102);

    auto tr_x_xxyyy_yzzz = pbuffer.data(idx_dip_hg + 103);

    auto tr_x_xxyyy_zzzz = pbuffer.data(idx_dip_hg + 104);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tr_x_xxy_xxxx,   \
                             tr_x_xxy_xxxy,   \
                             tr_x_xxy_xxxz,   \
                             tr_x_xxy_xxyy,   \
                             tr_x_xxy_xxyz,   \
                             tr_x_xxy_xxzz,   \
                             tr_x_xxy_xyyy,   \
                             tr_x_xxy_xyyz,   \
                             tr_x_xxy_xyzz,   \
                             tr_x_xxy_xzzz,   \
                             tr_x_xxy_zzzz,   \
                             tr_x_xxyy_xxx,   \
                             tr_x_xxyy_xxxx,  \
                             tr_x_xxyy_xxxy,  \
                             tr_x_xxyy_xxxz,  \
                             tr_x_xxyy_xxy,   \
                             tr_x_xxyy_xxyy,  \
                             tr_x_xxyy_xxyz,  \
                             tr_x_xxyy_xxz,   \
                             tr_x_xxyy_xxzz,  \
                             tr_x_xxyy_xyy,   \
                             tr_x_xxyy_xyyy,  \
                             tr_x_xxyy_xyyz,  \
                             tr_x_xxyy_xyz,   \
                             tr_x_xxyy_xyzz,  \
                             tr_x_xxyy_xzz,   \
                             tr_x_xxyy_xzzz,  \
                             tr_x_xxyy_zzzz,  \
                             tr_x_xxyyy_xxxx, \
                             tr_x_xxyyy_xxxy, \
                             tr_x_xxyyy_xxxz, \
                             tr_x_xxyyy_xxyy, \
                             tr_x_xxyyy_xxyz, \
                             tr_x_xxyyy_xxzz, \
                             tr_x_xxyyy_xyyy, \
                             tr_x_xxyyy_xyyz, \
                             tr_x_xxyyy_xyzz, \
                             tr_x_xxyyy_xzzz, \
                             tr_x_xxyyy_yyyy, \
                             tr_x_xxyyy_yyyz, \
                             tr_x_xxyyy_yyzz, \
                             tr_x_xxyyy_yzzz, \
                             tr_x_xxyyy_zzzz, \
                             tr_x_xyyy_yyyy,  \
                             tr_x_xyyy_yyyz,  \
                             tr_x_xyyy_yyzz,  \
                             tr_x_xyyy_yzzz,  \
                             tr_x_yyy_yyyy,   \
                             tr_x_yyy_yyyz,   \
                             tr_x_yyy_yyzz,   \
                             tr_x_yyy_yzzz,   \
                             ts_xyyy_yyyy,    \
                             ts_xyyy_yyyz,    \
                             ts_xyyy_yyzz,    \
                             ts_xyyy_yzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyy_xxxx[i] = 2.0 * tr_x_xxy_xxxx[i] * fe_0 + tr_x_xxyy_xxxx[i] * pa_y[i];

        tr_x_xxyyy_xxxy[i] = 2.0 * tr_x_xxy_xxxy[i] * fe_0 + tr_x_xxyy_xxx[i] * fe_0 + tr_x_xxyy_xxxy[i] * pa_y[i];

        tr_x_xxyyy_xxxz[i] = 2.0 * tr_x_xxy_xxxz[i] * fe_0 + tr_x_xxyy_xxxz[i] * pa_y[i];

        tr_x_xxyyy_xxyy[i] = 2.0 * tr_x_xxy_xxyy[i] * fe_0 + 2.0 * tr_x_xxyy_xxy[i] * fe_0 + tr_x_xxyy_xxyy[i] * pa_y[i];

        tr_x_xxyyy_xxyz[i] = 2.0 * tr_x_xxy_xxyz[i] * fe_0 + tr_x_xxyy_xxz[i] * fe_0 + tr_x_xxyy_xxyz[i] * pa_y[i];

        tr_x_xxyyy_xxzz[i] = 2.0 * tr_x_xxy_xxzz[i] * fe_0 + tr_x_xxyy_xxzz[i] * pa_y[i];

        tr_x_xxyyy_xyyy[i] = 2.0 * tr_x_xxy_xyyy[i] * fe_0 + 3.0 * tr_x_xxyy_xyy[i] * fe_0 + tr_x_xxyy_xyyy[i] * pa_y[i];

        tr_x_xxyyy_xyyz[i] = 2.0 * tr_x_xxy_xyyz[i] * fe_0 + 2.0 * tr_x_xxyy_xyz[i] * fe_0 + tr_x_xxyy_xyyz[i] * pa_y[i];

        tr_x_xxyyy_xyzz[i] = 2.0 * tr_x_xxy_xyzz[i] * fe_0 + tr_x_xxyy_xzz[i] * fe_0 + tr_x_xxyy_xyzz[i] * pa_y[i];

        tr_x_xxyyy_xzzz[i] = 2.0 * tr_x_xxy_xzzz[i] * fe_0 + tr_x_xxyy_xzzz[i] * pa_y[i];

        tr_x_xxyyy_yyyy[i] = tr_x_yyy_yyyy[i] * fe_0 + ts_xyyy_yyyy[i] * fe_0 + tr_x_xyyy_yyyy[i] * pa_x[i];

        tr_x_xxyyy_yyyz[i] = tr_x_yyy_yyyz[i] * fe_0 + ts_xyyy_yyyz[i] * fe_0 + tr_x_xyyy_yyyz[i] * pa_x[i];

        tr_x_xxyyy_yyzz[i] = tr_x_yyy_yyzz[i] * fe_0 + ts_xyyy_yyzz[i] * fe_0 + tr_x_xyyy_yyzz[i] * pa_x[i];

        tr_x_xxyyy_yzzz[i] = tr_x_yyy_yzzz[i] * fe_0 + ts_xyyy_yzzz[i] * fe_0 + tr_x_xyyy_yzzz[i] * pa_x[i];

        tr_x_xxyyy_zzzz[i] = 2.0 * tr_x_xxy_zzzz[i] * fe_0 + tr_x_xxyy_zzzz[i] * pa_y[i];
    }

    // Set up 105-120 components of targeted buffer : HG

    auto tr_x_xxyyz_xxxx = pbuffer.data(idx_dip_hg + 105);

    auto tr_x_xxyyz_xxxy = pbuffer.data(idx_dip_hg + 106);

    auto tr_x_xxyyz_xxxz = pbuffer.data(idx_dip_hg + 107);

    auto tr_x_xxyyz_xxyy = pbuffer.data(idx_dip_hg + 108);

    auto tr_x_xxyyz_xxyz = pbuffer.data(idx_dip_hg + 109);

    auto tr_x_xxyyz_xxzz = pbuffer.data(idx_dip_hg + 110);

    auto tr_x_xxyyz_xyyy = pbuffer.data(idx_dip_hg + 111);

    auto tr_x_xxyyz_xyyz = pbuffer.data(idx_dip_hg + 112);

    auto tr_x_xxyyz_xyzz = pbuffer.data(idx_dip_hg + 113);

    auto tr_x_xxyyz_xzzz = pbuffer.data(idx_dip_hg + 114);

    auto tr_x_xxyyz_yyyy = pbuffer.data(idx_dip_hg + 115);

    auto tr_x_xxyyz_yyyz = pbuffer.data(idx_dip_hg + 116);

    auto tr_x_xxyyz_yyzz = pbuffer.data(idx_dip_hg + 117);

    auto tr_x_xxyyz_yzzz = pbuffer.data(idx_dip_hg + 118);

    auto tr_x_xxyyz_zzzz = pbuffer.data(idx_dip_hg + 119);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             tr_x_xxyy_xxxx,  \
                             tr_x_xxyy_xxxy,  \
                             tr_x_xxyy_xxy,   \
                             tr_x_xxyy_xxyy,  \
                             tr_x_xxyy_xxyz,  \
                             tr_x_xxyy_xyy,   \
                             tr_x_xxyy_xyyy,  \
                             tr_x_xxyy_xyyz,  \
                             tr_x_xxyy_xyz,   \
                             tr_x_xxyy_xyzz,  \
                             tr_x_xxyy_yyy,   \
                             tr_x_xxyy_yyyy,  \
                             tr_x_xxyy_yyyz,  \
                             tr_x_xxyy_yyz,   \
                             tr_x_xxyy_yyzz,  \
                             tr_x_xxyy_yzz,   \
                             tr_x_xxyy_yzzz,  \
                             tr_x_xxyyz_xxxx, \
                             tr_x_xxyyz_xxxy, \
                             tr_x_xxyyz_xxxz, \
                             tr_x_xxyyz_xxyy, \
                             tr_x_xxyyz_xxyz, \
                             tr_x_xxyyz_xxzz, \
                             tr_x_xxyyz_xyyy, \
                             tr_x_xxyyz_xyyz, \
                             tr_x_xxyyz_xyzz, \
                             tr_x_xxyyz_xzzz, \
                             tr_x_xxyyz_yyyy, \
                             tr_x_xxyyz_yyyz, \
                             tr_x_xxyyz_yyzz, \
                             tr_x_xxyyz_yzzz, \
                             tr_x_xxyyz_zzzz, \
                             tr_x_xxyz_xxxz,  \
                             tr_x_xxyz_xxzz,  \
                             tr_x_xxyz_xzzz,  \
                             tr_x_xxyz_zzzz,  \
                             tr_x_xxz_xxxz,   \
                             tr_x_xxz_xxzz,   \
                             tr_x_xxz_xzzz,   \
                             tr_x_xxz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyz_xxxx[i] = tr_x_xxyy_xxxx[i] * pa_z[i];

        tr_x_xxyyz_xxxy[i] = tr_x_xxyy_xxxy[i] * pa_z[i];

        tr_x_xxyyz_xxxz[i] = tr_x_xxz_xxxz[i] * fe_0 + tr_x_xxyz_xxxz[i] * pa_y[i];

        tr_x_xxyyz_xxyy[i] = tr_x_xxyy_xxyy[i] * pa_z[i];

        tr_x_xxyyz_xxyz[i] = tr_x_xxyy_xxy[i] * fe_0 + tr_x_xxyy_xxyz[i] * pa_z[i];

        tr_x_xxyyz_xxzz[i] = tr_x_xxz_xxzz[i] * fe_0 + tr_x_xxyz_xxzz[i] * pa_y[i];

        tr_x_xxyyz_xyyy[i] = tr_x_xxyy_xyyy[i] * pa_z[i];

        tr_x_xxyyz_xyyz[i] = tr_x_xxyy_xyy[i] * fe_0 + tr_x_xxyy_xyyz[i] * pa_z[i];

        tr_x_xxyyz_xyzz[i] = 2.0 * tr_x_xxyy_xyz[i] * fe_0 + tr_x_xxyy_xyzz[i] * pa_z[i];

        tr_x_xxyyz_xzzz[i] = tr_x_xxz_xzzz[i] * fe_0 + tr_x_xxyz_xzzz[i] * pa_y[i];

        tr_x_xxyyz_yyyy[i] = tr_x_xxyy_yyyy[i] * pa_z[i];

        tr_x_xxyyz_yyyz[i] = tr_x_xxyy_yyy[i] * fe_0 + tr_x_xxyy_yyyz[i] * pa_z[i];

        tr_x_xxyyz_yyzz[i] = 2.0 * tr_x_xxyy_yyz[i] * fe_0 + tr_x_xxyy_yyzz[i] * pa_z[i];

        tr_x_xxyyz_yzzz[i] = 3.0 * tr_x_xxyy_yzz[i] * fe_0 + tr_x_xxyy_yzzz[i] * pa_z[i];

        tr_x_xxyyz_zzzz[i] = tr_x_xxz_zzzz[i] * fe_0 + tr_x_xxyz_zzzz[i] * pa_y[i];
    }

    // Set up 120-135 components of targeted buffer : HG

    auto tr_x_xxyzz_xxxx = pbuffer.data(idx_dip_hg + 120);

    auto tr_x_xxyzz_xxxy = pbuffer.data(idx_dip_hg + 121);

    auto tr_x_xxyzz_xxxz = pbuffer.data(idx_dip_hg + 122);

    auto tr_x_xxyzz_xxyy = pbuffer.data(idx_dip_hg + 123);

    auto tr_x_xxyzz_xxyz = pbuffer.data(idx_dip_hg + 124);

    auto tr_x_xxyzz_xxzz = pbuffer.data(idx_dip_hg + 125);

    auto tr_x_xxyzz_xyyy = pbuffer.data(idx_dip_hg + 126);

    auto tr_x_xxyzz_xyyz = pbuffer.data(idx_dip_hg + 127);

    auto tr_x_xxyzz_xyzz = pbuffer.data(idx_dip_hg + 128);

    auto tr_x_xxyzz_xzzz = pbuffer.data(idx_dip_hg + 129);

    auto tr_x_xxyzz_yyyy = pbuffer.data(idx_dip_hg + 130);

    auto tr_x_xxyzz_yyyz = pbuffer.data(idx_dip_hg + 131);

    auto tr_x_xxyzz_yyzz = pbuffer.data(idx_dip_hg + 132);

    auto tr_x_xxyzz_yzzz = pbuffer.data(idx_dip_hg + 133);

    auto tr_x_xxyzz_zzzz = pbuffer.data(idx_dip_hg + 134);

#pragma omp simd aligned(pa_y,                \
                             tr_x_xxyzz_xxxx, \
                             tr_x_xxyzz_xxxy, \
                             tr_x_xxyzz_xxxz, \
                             tr_x_xxyzz_xxyy, \
                             tr_x_xxyzz_xxyz, \
                             tr_x_xxyzz_xxzz, \
                             tr_x_xxyzz_xyyy, \
                             tr_x_xxyzz_xyyz, \
                             tr_x_xxyzz_xyzz, \
                             tr_x_xxyzz_xzzz, \
                             tr_x_xxyzz_yyyy, \
                             tr_x_xxyzz_yyyz, \
                             tr_x_xxyzz_yyzz, \
                             tr_x_xxyzz_yzzz, \
                             tr_x_xxyzz_zzzz, \
                             tr_x_xxzz_xxx,   \
                             tr_x_xxzz_xxxx,  \
                             tr_x_xxzz_xxxy,  \
                             tr_x_xxzz_xxxz,  \
                             tr_x_xxzz_xxy,   \
                             tr_x_xxzz_xxyy,  \
                             tr_x_xxzz_xxyz,  \
                             tr_x_xxzz_xxz,   \
                             tr_x_xxzz_xxzz,  \
                             tr_x_xxzz_xyy,   \
                             tr_x_xxzz_xyyy,  \
                             tr_x_xxzz_xyyz,  \
                             tr_x_xxzz_xyz,   \
                             tr_x_xxzz_xyzz,  \
                             tr_x_xxzz_xzz,   \
                             tr_x_xxzz_xzzz,  \
                             tr_x_xxzz_yyy,   \
                             tr_x_xxzz_yyyy,  \
                             tr_x_xxzz_yyyz,  \
                             tr_x_xxzz_yyz,   \
                             tr_x_xxzz_yyzz,  \
                             tr_x_xxzz_yzz,   \
                             tr_x_xxzz_yzzz,  \
                             tr_x_xxzz_zzz,   \
                             tr_x_xxzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyzz_xxxx[i] = tr_x_xxzz_xxxx[i] * pa_y[i];

        tr_x_xxyzz_xxxy[i] = tr_x_xxzz_xxx[i] * fe_0 + tr_x_xxzz_xxxy[i] * pa_y[i];

        tr_x_xxyzz_xxxz[i] = tr_x_xxzz_xxxz[i] * pa_y[i];

        tr_x_xxyzz_xxyy[i] = 2.0 * tr_x_xxzz_xxy[i] * fe_0 + tr_x_xxzz_xxyy[i] * pa_y[i];

        tr_x_xxyzz_xxyz[i] = tr_x_xxzz_xxz[i] * fe_0 + tr_x_xxzz_xxyz[i] * pa_y[i];

        tr_x_xxyzz_xxzz[i] = tr_x_xxzz_xxzz[i] * pa_y[i];

        tr_x_xxyzz_xyyy[i] = 3.0 * tr_x_xxzz_xyy[i] * fe_0 + tr_x_xxzz_xyyy[i] * pa_y[i];

        tr_x_xxyzz_xyyz[i] = 2.0 * tr_x_xxzz_xyz[i] * fe_0 + tr_x_xxzz_xyyz[i] * pa_y[i];

        tr_x_xxyzz_xyzz[i] = tr_x_xxzz_xzz[i] * fe_0 + tr_x_xxzz_xyzz[i] * pa_y[i];

        tr_x_xxyzz_xzzz[i] = tr_x_xxzz_xzzz[i] * pa_y[i];

        tr_x_xxyzz_yyyy[i] = 4.0 * tr_x_xxzz_yyy[i] * fe_0 + tr_x_xxzz_yyyy[i] * pa_y[i];

        tr_x_xxyzz_yyyz[i] = 3.0 * tr_x_xxzz_yyz[i] * fe_0 + tr_x_xxzz_yyyz[i] * pa_y[i];

        tr_x_xxyzz_yyzz[i] = 2.0 * tr_x_xxzz_yzz[i] * fe_0 + tr_x_xxzz_yyzz[i] * pa_y[i];

        tr_x_xxyzz_yzzz[i] = tr_x_xxzz_zzz[i] * fe_0 + tr_x_xxzz_yzzz[i] * pa_y[i];

        tr_x_xxyzz_zzzz[i] = tr_x_xxzz_zzzz[i] * pa_y[i];
    }

    // Set up 135-150 components of targeted buffer : HG

    auto tr_x_xxzzz_xxxx = pbuffer.data(idx_dip_hg + 135);

    auto tr_x_xxzzz_xxxy = pbuffer.data(idx_dip_hg + 136);

    auto tr_x_xxzzz_xxxz = pbuffer.data(idx_dip_hg + 137);

    auto tr_x_xxzzz_xxyy = pbuffer.data(idx_dip_hg + 138);

    auto tr_x_xxzzz_xxyz = pbuffer.data(idx_dip_hg + 139);

    auto tr_x_xxzzz_xxzz = pbuffer.data(idx_dip_hg + 140);

    auto tr_x_xxzzz_xyyy = pbuffer.data(idx_dip_hg + 141);

    auto tr_x_xxzzz_xyyz = pbuffer.data(idx_dip_hg + 142);

    auto tr_x_xxzzz_xyzz = pbuffer.data(idx_dip_hg + 143);

    auto tr_x_xxzzz_xzzz = pbuffer.data(idx_dip_hg + 144);

    auto tr_x_xxzzz_yyyy = pbuffer.data(idx_dip_hg + 145);

    auto tr_x_xxzzz_yyyz = pbuffer.data(idx_dip_hg + 146);

    auto tr_x_xxzzz_yyzz = pbuffer.data(idx_dip_hg + 147);

    auto tr_x_xxzzz_yzzz = pbuffer.data(idx_dip_hg + 148);

    auto tr_x_xxzzz_zzzz = pbuffer.data(idx_dip_hg + 149);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             tr_x_xxz_xxxx,   \
                             tr_x_xxz_xxxy,   \
                             tr_x_xxz_xxxz,   \
                             tr_x_xxz_xxyy,   \
                             tr_x_xxz_xxyz,   \
                             tr_x_xxz_xxzz,   \
                             tr_x_xxz_xyyy,   \
                             tr_x_xxz_xyyz,   \
                             tr_x_xxz_xyzz,   \
                             tr_x_xxz_xzzz,   \
                             tr_x_xxz_yyyy,   \
                             tr_x_xxzz_xxx,   \
                             tr_x_xxzz_xxxx,  \
                             tr_x_xxzz_xxxy,  \
                             tr_x_xxzz_xxxz,  \
                             tr_x_xxzz_xxy,   \
                             tr_x_xxzz_xxyy,  \
                             tr_x_xxzz_xxyz,  \
                             tr_x_xxzz_xxz,   \
                             tr_x_xxzz_xxzz,  \
                             tr_x_xxzz_xyy,   \
                             tr_x_xxzz_xyyy,  \
                             tr_x_xxzz_xyyz,  \
                             tr_x_xxzz_xyz,   \
                             tr_x_xxzz_xyzz,  \
                             tr_x_xxzz_xzz,   \
                             tr_x_xxzz_xzzz,  \
                             tr_x_xxzz_yyyy,  \
                             tr_x_xxzzz_xxxx, \
                             tr_x_xxzzz_xxxy, \
                             tr_x_xxzzz_xxxz, \
                             tr_x_xxzzz_xxyy, \
                             tr_x_xxzzz_xxyz, \
                             tr_x_xxzzz_xxzz, \
                             tr_x_xxzzz_xyyy, \
                             tr_x_xxzzz_xyyz, \
                             tr_x_xxzzz_xyzz, \
                             tr_x_xxzzz_xzzz, \
                             tr_x_xxzzz_yyyy, \
                             tr_x_xxzzz_yyyz, \
                             tr_x_xxzzz_yyzz, \
                             tr_x_xxzzz_yzzz, \
                             tr_x_xxzzz_zzzz, \
                             tr_x_xzzz_yyyz,  \
                             tr_x_xzzz_yyzz,  \
                             tr_x_xzzz_yzzz,  \
                             tr_x_xzzz_zzzz,  \
                             tr_x_zzz_yyyz,   \
                             tr_x_zzz_yyzz,   \
                             tr_x_zzz_yzzz,   \
                             tr_x_zzz_zzzz,   \
                             ts_xzzz_yyyz,    \
                             ts_xzzz_yyzz,    \
                             ts_xzzz_yzzz,    \
                             ts_xzzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzzz_xxxx[i] = 2.0 * tr_x_xxz_xxxx[i] * fe_0 + tr_x_xxzz_xxxx[i] * pa_z[i];

        tr_x_xxzzz_xxxy[i] = 2.0 * tr_x_xxz_xxxy[i] * fe_0 + tr_x_xxzz_xxxy[i] * pa_z[i];

        tr_x_xxzzz_xxxz[i] = 2.0 * tr_x_xxz_xxxz[i] * fe_0 + tr_x_xxzz_xxx[i] * fe_0 + tr_x_xxzz_xxxz[i] * pa_z[i];

        tr_x_xxzzz_xxyy[i] = 2.0 * tr_x_xxz_xxyy[i] * fe_0 + tr_x_xxzz_xxyy[i] * pa_z[i];

        tr_x_xxzzz_xxyz[i] = 2.0 * tr_x_xxz_xxyz[i] * fe_0 + tr_x_xxzz_xxy[i] * fe_0 + tr_x_xxzz_xxyz[i] * pa_z[i];

        tr_x_xxzzz_xxzz[i] = 2.0 * tr_x_xxz_xxzz[i] * fe_0 + 2.0 * tr_x_xxzz_xxz[i] * fe_0 + tr_x_xxzz_xxzz[i] * pa_z[i];

        tr_x_xxzzz_xyyy[i] = 2.0 * tr_x_xxz_xyyy[i] * fe_0 + tr_x_xxzz_xyyy[i] * pa_z[i];

        tr_x_xxzzz_xyyz[i] = 2.0 * tr_x_xxz_xyyz[i] * fe_0 + tr_x_xxzz_xyy[i] * fe_0 + tr_x_xxzz_xyyz[i] * pa_z[i];

        tr_x_xxzzz_xyzz[i] = 2.0 * tr_x_xxz_xyzz[i] * fe_0 + 2.0 * tr_x_xxzz_xyz[i] * fe_0 + tr_x_xxzz_xyzz[i] * pa_z[i];

        tr_x_xxzzz_xzzz[i] = 2.0 * tr_x_xxz_xzzz[i] * fe_0 + 3.0 * tr_x_xxzz_xzz[i] * fe_0 + tr_x_xxzz_xzzz[i] * pa_z[i];

        tr_x_xxzzz_yyyy[i] = 2.0 * tr_x_xxz_yyyy[i] * fe_0 + tr_x_xxzz_yyyy[i] * pa_z[i];

        tr_x_xxzzz_yyyz[i] = tr_x_zzz_yyyz[i] * fe_0 + ts_xzzz_yyyz[i] * fe_0 + tr_x_xzzz_yyyz[i] * pa_x[i];

        tr_x_xxzzz_yyzz[i] = tr_x_zzz_yyzz[i] * fe_0 + ts_xzzz_yyzz[i] * fe_0 + tr_x_xzzz_yyzz[i] * pa_x[i];

        tr_x_xxzzz_yzzz[i] = tr_x_zzz_yzzz[i] * fe_0 + ts_xzzz_yzzz[i] * fe_0 + tr_x_xzzz_yzzz[i] * pa_x[i];

        tr_x_xxzzz_zzzz[i] = tr_x_zzz_zzzz[i] * fe_0 + ts_xzzz_zzzz[i] * fe_0 + tr_x_xzzz_zzzz[i] * pa_x[i];
    }

    // Set up 150-165 components of targeted buffer : HG

    auto tr_x_xyyyy_xxxx = pbuffer.data(idx_dip_hg + 150);

    auto tr_x_xyyyy_xxxy = pbuffer.data(idx_dip_hg + 151);

    auto tr_x_xyyyy_xxxz = pbuffer.data(idx_dip_hg + 152);

    auto tr_x_xyyyy_xxyy = pbuffer.data(idx_dip_hg + 153);

    auto tr_x_xyyyy_xxyz = pbuffer.data(idx_dip_hg + 154);

    auto tr_x_xyyyy_xxzz = pbuffer.data(idx_dip_hg + 155);

    auto tr_x_xyyyy_xyyy = pbuffer.data(idx_dip_hg + 156);

    auto tr_x_xyyyy_xyyz = pbuffer.data(idx_dip_hg + 157);

    auto tr_x_xyyyy_xyzz = pbuffer.data(idx_dip_hg + 158);

    auto tr_x_xyyyy_xzzz = pbuffer.data(idx_dip_hg + 159);

    auto tr_x_xyyyy_yyyy = pbuffer.data(idx_dip_hg + 160);

    auto tr_x_xyyyy_yyyz = pbuffer.data(idx_dip_hg + 161);

    auto tr_x_xyyyy_yyzz = pbuffer.data(idx_dip_hg + 162);

    auto tr_x_xyyyy_yzzz = pbuffer.data(idx_dip_hg + 163);

    auto tr_x_xyyyy_zzzz = pbuffer.data(idx_dip_hg + 164);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tr_x_xyy_xxxx,   \
                             tr_x_xyy_xxxz,   \
                             tr_x_xyy_xxzz,   \
                             tr_x_xyy_xzzz,   \
                             tr_x_xyyy_xxxx,  \
                             tr_x_xyyy_xxxz,  \
                             tr_x_xyyy_xxzz,  \
                             tr_x_xyyy_xzzz,  \
                             tr_x_xyyyy_xxxx, \
                             tr_x_xyyyy_xxxy, \
                             tr_x_xyyyy_xxxz, \
                             tr_x_xyyyy_xxyy, \
                             tr_x_xyyyy_xxyz, \
                             tr_x_xyyyy_xxzz, \
                             tr_x_xyyyy_xyyy, \
                             tr_x_xyyyy_xyyz, \
                             tr_x_xyyyy_xyzz, \
                             tr_x_xyyyy_xzzz, \
                             tr_x_xyyyy_yyyy, \
                             tr_x_xyyyy_yyyz, \
                             tr_x_xyyyy_yyzz, \
                             tr_x_xyyyy_yzzz, \
                             tr_x_xyyyy_zzzz, \
                             tr_x_yyyy_xxxy,  \
                             tr_x_yyyy_xxy,   \
                             tr_x_yyyy_xxyy,  \
                             tr_x_yyyy_xxyz,  \
                             tr_x_yyyy_xyy,   \
                             tr_x_yyyy_xyyy,  \
                             tr_x_yyyy_xyyz,  \
                             tr_x_yyyy_xyz,   \
                             tr_x_yyyy_xyzz,  \
                             tr_x_yyyy_yyy,   \
                             tr_x_yyyy_yyyy,  \
                             tr_x_yyyy_yyyz,  \
                             tr_x_yyyy_yyz,   \
                             tr_x_yyyy_yyzz,  \
                             tr_x_yyyy_yzz,   \
                             tr_x_yyyy_yzzz,  \
                             tr_x_yyyy_zzzz,  \
                             ts_yyyy_xxxy,    \
                             ts_yyyy_xxyy,    \
                             ts_yyyy_xxyz,    \
                             ts_yyyy_xyyy,    \
                             ts_yyyy_xyyz,    \
                             ts_yyyy_xyzz,    \
                             ts_yyyy_yyyy,    \
                             ts_yyyy_yyyz,    \
                             ts_yyyy_yyzz,    \
                             ts_yyyy_yzzz,    \
                             ts_yyyy_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyy_xxxx[i] = 3.0 * tr_x_xyy_xxxx[i] * fe_0 + tr_x_xyyy_xxxx[i] * pa_y[i];

        tr_x_xyyyy_xxxy[i] = 3.0 * tr_x_yyyy_xxy[i] * fe_0 + ts_yyyy_xxxy[i] * fe_0 + tr_x_yyyy_xxxy[i] * pa_x[i];

        tr_x_xyyyy_xxxz[i] = 3.0 * tr_x_xyy_xxxz[i] * fe_0 + tr_x_xyyy_xxxz[i] * pa_y[i];

        tr_x_xyyyy_xxyy[i] = 2.0 * tr_x_yyyy_xyy[i] * fe_0 + ts_yyyy_xxyy[i] * fe_0 + tr_x_yyyy_xxyy[i] * pa_x[i];

        tr_x_xyyyy_xxyz[i] = 2.0 * tr_x_yyyy_xyz[i] * fe_0 + ts_yyyy_xxyz[i] * fe_0 + tr_x_yyyy_xxyz[i] * pa_x[i];

        tr_x_xyyyy_xxzz[i] = 3.0 * tr_x_xyy_xxzz[i] * fe_0 + tr_x_xyyy_xxzz[i] * pa_y[i];

        tr_x_xyyyy_xyyy[i] = tr_x_yyyy_yyy[i] * fe_0 + ts_yyyy_xyyy[i] * fe_0 + tr_x_yyyy_xyyy[i] * pa_x[i];

        tr_x_xyyyy_xyyz[i] = tr_x_yyyy_yyz[i] * fe_0 + ts_yyyy_xyyz[i] * fe_0 + tr_x_yyyy_xyyz[i] * pa_x[i];

        tr_x_xyyyy_xyzz[i] = tr_x_yyyy_yzz[i] * fe_0 + ts_yyyy_xyzz[i] * fe_0 + tr_x_yyyy_xyzz[i] * pa_x[i];

        tr_x_xyyyy_xzzz[i] = 3.0 * tr_x_xyy_xzzz[i] * fe_0 + tr_x_xyyy_xzzz[i] * pa_y[i];

        tr_x_xyyyy_yyyy[i] = ts_yyyy_yyyy[i] * fe_0 + tr_x_yyyy_yyyy[i] * pa_x[i];

        tr_x_xyyyy_yyyz[i] = ts_yyyy_yyyz[i] * fe_0 + tr_x_yyyy_yyyz[i] * pa_x[i];

        tr_x_xyyyy_yyzz[i] = ts_yyyy_yyzz[i] * fe_0 + tr_x_yyyy_yyzz[i] * pa_x[i];

        tr_x_xyyyy_yzzz[i] = ts_yyyy_yzzz[i] * fe_0 + tr_x_yyyy_yzzz[i] * pa_x[i];

        tr_x_xyyyy_zzzz[i] = ts_yyyy_zzzz[i] * fe_0 + tr_x_yyyy_zzzz[i] * pa_x[i];
    }

    // Set up 165-180 components of targeted buffer : HG

    auto tr_x_xyyyz_xxxx = pbuffer.data(idx_dip_hg + 165);

    auto tr_x_xyyyz_xxxy = pbuffer.data(idx_dip_hg + 166);

    auto tr_x_xyyyz_xxxz = pbuffer.data(idx_dip_hg + 167);

    auto tr_x_xyyyz_xxyy = pbuffer.data(idx_dip_hg + 168);

    auto tr_x_xyyyz_xxyz = pbuffer.data(idx_dip_hg + 169);

    auto tr_x_xyyyz_xxzz = pbuffer.data(idx_dip_hg + 170);

    auto tr_x_xyyyz_xyyy = pbuffer.data(idx_dip_hg + 171);

    auto tr_x_xyyyz_xyyz = pbuffer.data(idx_dip_hg + 172);

    auto tr_x_xyyyz_xyzz = pbuffer.data(idx_dip_hg + 173);

    auto tr_x_xyyyz_xzzz = pbuffer.data(idx_dip_hg + 174);

    auto tr_x_xyyyz_yyyy = pbuffer.data(idx_dip_hg + 175);

    auto tr_x_xyyyz_yyyz = pbuffer.data(idx_dip_hg + 176);

    auto tr_x_xyyyz_yyzz = pbuffer.data(idx_dip_hg + 177);

    auto tr_x_xyyyz_yzzz = pbuffer.data(idx_dip_hg + 178);

    auto tr_x_xyyyz_zzzz = pbuffer.data(idx_dip_hg + 179);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             tr_x_xyyy_xxxx,  \
                             tr_x_xyyy_xxxy,  \
                             tr_x_xyyy_xxy,   \
                             tr_x_xyyy_xxyy,  \
                             tr_x_xyyy_xxyz,  \
                             tr_x_xyyy_xyy,   \
                             tr_x_xyyy_xyyy,  \
                             tr_x_xyyy_xyyz,  \
                             tr_x_xyyy_xyz,   \
                             tr_x_xyyy_xyzz,  \
                             tr_x_xyyy_yyyy,  \
                             tr_x_xyyyz_xxxx, \
                             tr_x_xyyyz_xxxy, \
                             tr_x_xyyyz_xxxz, \
                             tr_x_xyyyz_xxyy, \
                             tr_x_xyyyz_xxyz, \
                             tr_x_xyyyz_xxzz, \
                             tr_x_xyyyz_xyyy, \
                             tr_x_xyyyz_xyyz, \
                             tr_x_xyyyz_xyzz, \
                             tr_x_xyyyz_xzzz, \
                             tr_x_xyyyz_yyyy, \
                             tr_x_xyyyz_yyyz, \
                             tr_x_xyyyz_yyzz, \
                             tr_x_xyyyz_yzzz, \
                             tr_x_xyyyz_zzzz, \
                             tr_x_xyyz_xxxz,  \
                             tr_x_xyyz_xxzz,  \
                             tr_x_xyyz_xzzz,  \
                             tr_x_xyz_xxxz,   \
                             tr_x_xyz_xxzz,   \
                             tr_x_xyz_xzzz,   \
                             tr_x_yyyz_yyyz,  \
                             tr_x_yyyz_yyzz,  \
                             tr_x_yyyz_yzzz,  \
                             tr_x_yyyz_zzzz,  \
                             ts_yyyz_yyyz,    \
                             ts_yyyz_yyzz,    \
                             ts_yyyz_yzzz,    \
                             ts_yyyz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyz_xxxx[i] = tr_x_xyyy_xxxx[i] * pa_z[i];

        tr_x_xyyyz_xxxy[i] = tr_x_xyyy_xxxy[i] * pa_z[i];

        tr_x_xyyyz_xxxz[i] = 2.0 * tr_x_xyz_xxxz[i] * fe_0 + tr_x_xyyz_xxxz[i] * pa_y[i];

        tr_x_xyyyz_xxyy[i] = tr_x_xyyy_xxyy[i] * pa_z[i];

        tr_x_xyyyz_xxyz[i] = tr_x_xyyy_xxy[i] * fe_0 + tr_x_xyyy_xxyz[i] * pa_z[i];

        tr_x_xyyyz_xxzz[i] = 2.0 * tr_x_xyz_xxzz[i] * fe_0 + tr_x_xyyz_xxzz[i] * pa_y[i];

        tr_x_xyyyz_xyyy[i] = tr_x_xyyy_xyyy[i] * pa_z[i];

        tr_x_xyyyz_xyyz[i] = tr_x_xyyy_xyy[i] * fe_0 + tr_x_xyyy_xyyz[i] * pa_z[i];

        tr_x_xyyyz_xyzz[i] = 2.0 * tr_x_xyyy_xyz[i] * fe_0 + tr_x_xyyy_xyzz[i] * pa_z[i];

        tr_x_xyyyz_xzzz[i] = 2.0 * tr_x_xyz_xzzz[i] * fe_0 + tr_x_xyyz_xzzz[i] * pa_y[i];

        tr_x_xyyyz_yyyy[i] = tr_x_xyyy_yyyy[i] * pa_z[i];

        tr_x_xyyyz_yyyz[i] = ts_yyyz_yyyz[i] * fe_0 + tr_x_yyyz_yyyz[i] * pa_x[i];

        tr_x_xyyyz_yyzz[i] = ts_yyyz_yyzz[i] * fe_0 + tr_x_yyyz_yyzz[i] * pa_x[i];

        tr_x_xyyyz_yzzz[i] = ts_yyyz_yzzz[i] * fe_0 + tr_x_yyyz_yzzz[i] * pa_x[i];

        tr_x_xyyyz_zzzz[i] = ts_yyyz_zzzz[i] * fe_0 + tr_x_yyyz_zzzz[i] * pa_x[i];
    }

    // Set up 180-195 components of targeted buffer : HG

    auto tr_x_xyyzz_xxxx = pbuffer.data(idx_dip_hg + 180);

    auto tr_x_xyyzz_xxxy = pbuffer.data(idx_dip_hg + 181);

    auto tr_x_xyyzz_xxxz = pbuffer.data(idx_dip_hg + 182);

    auto tr_x_xyyzz_xxyy = pbuffer.data(idx_dip_hg + 183);

    auto tr_x_xyyzz_xxyz = pbuffer.data(idx_dip_hg + 184);

    auto tr_x_xyyzz_xxzz = pbuffer.data(idx_dip_hg + 185);

    auto tr_x_xyyzz_xyyy = pbuffer.data(idx_dip_hg + 186);

    auto tr_x_xyyzz_xyyz = pbuffer.data(idx_dip_hg + 187);

    auto tr_x_xyyzz_xyzz = pbuffer.data(idx_dip_hg + 188);

    auto tr_x_xyyzz_xzzz = pbuffer.data(idx_dip_hg + 189);

    auto tr_x_xyyzz_yyyy = pbuffer.data(idx_dip_hg + 190);

    auto tr_x_xyyzz_yyyz = pbuffer.data(idx_dip_hg + 191);

    auto tr_x_xyyzz_yyzz = pbuffer.data(idx_dip_hg + 192);

    auto tr_x_xyyzz_yzzz = pbuffer.data(idx_dip_hg + 193);

    auto tr_x_xyyzz_zzzz = pbuffer.data(idx_dip_hg + 194);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             tr_x_xyy_xxxy,   \
                             tr_x_xyy_xxyy,   \
                             tr_x_xyy_xyyy,   \
                             tr_x_xyyz_xxxy,  \
                             tr_x_xyyz_xxyy,  \
                             tr_x_xyyz_xyyy,  \
                             tr_x_xyyzz_xxxx, \
                             tr_x_xyyzz_xxxy, \
                             tr_x_xyyzz_xxxz, \
                             tr_x_xyyzz_xxyy, \
                             tr_x_xyyzz_xxyz, \
                             tr_x_xyyzz_xxzz, \
                             tr_x_xyyzz_xyyy, \
                             tr_x_xyyzz_xyyz, \
                             tr_x_xyyzz_xyzz, \
                             tr_x_xyyzz_xzzz, \
                             tr_x_xyyzz_yyyy, \
                             tr_x_xyyzz_yyyz, \
                             tr_x_xyyzz_yyzz, \
                             tr_x_xyyzz_yzzz, \
                             tr_x_xyyzz_zzzz, \
                             tr_x_xyzz_xxxx,  \
                             tr_x_xyzz_xxxz,  \
                             tr_x_xyzz_xxzz,  \
                             tr_x_xyzz_xzzz,  \
                             tr_x_xzz_xxxx,   \
                             tr_x_xzz_xxxz,   \
                             tr_x_xzz_xxzz,   \
                             tr_x_xzz_xzzz,   \
                             tr_x_yyzz_xxyz,  \
                             tr_x_yyzz_xyyz,  \
                             tr_x_yyzz_xyz,   \
                             tr_x_yyzz_xyzz,  \
                             tr_x_yyzz_yyyy,  \
                             tr_x_yyzz_yyyz,  \
                             tr_x_yyzz_yyz,   \
                             tr_x_yyzz_yyzz,  \
                             tr_x_yyzz_yzz,   \
                             tr_x_yyzz_yzzz,  \
                             tr_x_yyzz_zzzz,  \
                             ts_yyzz_xxyz,    \
                             ts_yyzz_xyyz,    \
                             ts_yyzz_xyzz,    \
                             ts_yyzz_yyyy,    \
                             ts_yyzz_yyyz,    \
                             ts_yyzz_yyzz,    \
                             ts_yyzz_yzzz,    \
                             ts_yyzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyzz_xxxx[i] = tr_x_xzz_xxxx[i] * fe_0 + tr_x_xyzz_xxxx[i] * pa_y[i];

        tr_x_xyyzz_xxxy[i] = tr_x_xyy_xxxy[i] * fe_0 + tr_x_xyyz_xxxy[i] * pa_z[i];

        tr_x_xyyzz_xxxz[i] = tr_x_xzz_xxxz[i] * fe_0 + tr_x_xyzz_xxxz[i] * pa_y[i];

        tr_x_xyyzz_xxyy[i] = tr_x_xyy_xxyy[i] * fe_0 + tr_x_xyyz_xxyy[i] * pa_z[i];

        tr_x_xyyzz_xxyz[i] = 2.0 * tr_x_yyzz_xyz[i] * fe_0 + ts_yyzz_xxyz[i] * fe_0 + tr_x_yyzz_xxyz[i] * pa_x[i];

        tr_x_xyyzz_xxzz[i] = tr_x_xzz_xxzz[i] * fe_0 + tr_x_xyzz_xxzz[i] * pa_y[i];

        tr_x_xyyzz_xyyy[i] = tr_x_xyy_xyyy[i] * fe_0 + tr_x_xyyz_xyyy[i] * pa_z[i];

        tr_x_xyyzz_xyyz[i] = tr_x_yyzz_yyz[i] * fe_0 + ts_yyzz_xyyz[i] * fe_0 + tr_x_yyzz_xyyz[i] * pa_x[i];

        tr_x_xyyzz_xyzz[i] = tr_x_yyzz_yzz[i] * fe_0 + ts_yyzz_xyzz[i] * fe_0 + tr_x_yyzz_xyzz[i] * pa_x[i];

        tr_x_xyyzz_xzzz[i] = tr_x_xzz_xzzz[i] * fe_0 + tr_x_xyzz_xzzz[i] * pa_y[i];

        tr_x_xyyzz_yyyy[i] = ts_yyzz_yyyy[i] * fe_0 + tr_x_yyzz_yyyy[i] * pa_x[i];

        tr_x_xyyzz_yyyz[i] = ts_yyzz_yyyz[i] * fe_0 + tr_x_yyzz_yyyz[i] * pa_x[i];

        tr_x_xyyzz_yyzz[i] = ts_yyzz_yyzz[i] * fe_0 + tr_x_yyzz_yyzz[i] * pa_x[i];

        tr_x_xyyzz_yzzz[i] = ts_yyzz_yzzz[i] * fe_0 + tr_x_yyzz_yzzz[i] * pa_x[i];

        tr_x_xyyzz_zzzz[i] = ts_yyzz_zzzz[i] * fe_0 + tr_x_yyzz_zzzz[i] * pa_x[i];
    }

    // Set up 195-210 components of targeted buffer : HG

    auto tr_x_xyzzz_xxxx = pbuffer.data(idx_dip_hg + 195);

    auto tr_x_xyzzz_xxxy = pbuffer.data(idx_dip_hg + 196);

    auto tr_x_xyzzz_xxxz = pbuffer.data(idx_dip_hg + 197);

    auto tr_x_xyzzz_xxyy = pbuffer.data(idx_dip_hg + 198);

    auto tr_x_xyzzz_xxyz = pbuffer.data(idx_dip_hg + 199);

    auto tr_x_xyzzz_xxzz = pbuffer.data(idx_dip_hg + 200);

    auto tr_x_xyzzz_xyyy = pbuffer.data(idx_dip_hg + 201);

    auto tr_x_xyzzz_xyyz = pbuffer.data(idx_dip_hg + 202);

    auto tr_x_xyzzz_xyzz = pbuffer.data(idx_dip_hg + 203);

    auto tr_x_xyzzz_xzzz = pbuffer.data(idx_dip_hg + 204);

    auto tr_x_xyzzz_yyyy = pbuffer.data(idx_dip_hg + 205);

    auto tr_x_xyzzz_yyyz = pbuffer.data(idx_dip_hg + 206);

    auto tr_x_xyzzz_yyzz = pbuffer.data(idx_dip_hg + 207);

    auto tr_x_xyzzz_yzzz = pbuffer.data(idx_dip_hg + 208);

    auto tr_x_xyzzz_zzzz = pbuffer.data(idx_dip_hg + 209);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tr_x_xyzzz_xxxx, \
                             tr_x_xyzzz_xxxy, \
                             tr_x_xyzzz_xxxz, \
                             tr_x_xyzzz_xxyy, \
                             tr_x_xyzzz_xxyz, \
                             tr_x_xyzzz_xxzz, \
                             tr_x_xyzzz_xyyy, \
                             tr_x_xyzzz_xyyz, \
                             tr_x_xyzzz_xyzz, \
                             tr_x_xyzzz_xzzz, \
                             tr_x_xyzzz_yyyy, \
                             tr_x_xyzzz_yyyz, \
                             tr_x_xyzzz_yyzz, \
                             tr_x_xyzzz_yzzz, \
                             tr_x_xyzzz_zzzz, \
                             tr_x_xzzz_xxx,   \
                             tr_x_xzzz_xxxx,  \
                             tr_x_xzzz_xxxy,  \
                             tr_x_xzzz_xxxz,  \
                             tr_x_xzzz_xxy,   \
                             tr_x_xzzz_xxyy,  \
                             tr_x_xzzz_xxyz,  \
                             tr_x_xzzz_xxz,   \
                             tr_x_xzzz_xxzz,  \
                             tr_x_xzzz_xyy,   \
                             tr_x_xzzz_xyyy,  \
                             tr_x_xzzz_xyyz,  \
                             tr_x_xzzz_xyz,   \
                             tr_x_xzzz_xyzz,  \
                             tr_x_xzzz_xzz,   \
                             tr_x_xzzz_xzzz,  \
                             tr_x_xzzz_zzzz,  \
                             tr_x_yzzz_yyyy,  \
                             tr_x_yzzz_yyyz,  \
                             tr_x_yzzz_yyzz,  \
                             tr_x_yzzz_yzzz,  \
                             ts_yzzz_yyyy,    \
                             ts_yzzz_yyyz,    \
                             ts_yzzz_yyzz,    \
                             ts_yzzz_yzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzzz_xxxx[i] = tr_x_xzzz_xxxx[i] * pa_y[i];

        tr_x_xyzzz_xxxy[i] = tr_x_xzzz_xxx[i] * fe_0 + tr_x_xzzz_xxxy[i] * pa_y[i];

        tr_x_xyzzz_xxxz[i] = tr_x_xzzz_xxxz[i] * pa_y[i];

        tr_x_xyzzz_xxyy[i] = 2.0 * tr_x_xzzz_xxy[i] * fe_0 + tr_x_xzzz_xxyy[i] * pa_y[i];

        tr_x_xyzzz_xxyz[i] = tr_x_xzzz_xxz[i] * fe_0 + tr_x_xzzz_xxyz[i] * pa_y[i];

        tr_x_xyzzz_xxzz[i] = tr_x_xzzz_xxzz[i] * pa_y[i];

        tr_x_xyzzz_xyyy[i] = 3.0 * tr_x_xzzz_xyy[i] * fe_0 + tr_x_xzzz_xyyy[i] * pa_y[i];

        tr_x_xyzzz_xyyz[i] = 2.0 * tr_x_xzzz_xyz[i] * fe_0 + tr_x_xzzz_xyyz[i] * pa_y[i];

        tr_x_xyzzz_xyzz[i] = tr_x_xzzz_xzz[i] * fe_0 + tr_x_xzzz_xyzz[i] * pa_y[i];

        tr_x_xyzzz_xzzz[i] = tr_x_xzzz_xzzz[i] * pa_y[i];

        tr_x_xyzzz_yyyy[i] = ts_yzzz_yyyy[i] * fe_0 + tr_x_yzzz_yyyy[i] * pa_x[i];

        tr_x_xyzzz_yyyz[i] = ts_yzzz_yyyz[i] * fe_0 + tr_x_yzzz_yyyz[i] * pa_x[i];

        tr_x_xyzzz_yyzz[i] = ts_yzzz_yyzz[i] * fe_0 + tr_x_yzzz_yyzz[i] * pa_x[i];

        tr_x_xyzzz_yzzz[i] = ts_yzzz_yzzz[i] * fe_0 + tr_x_yzzz_yzzz[i] * pa_x[i];

        tr_x_xyzzz_zzzz[i] = tr_x_xzzz_zzzz[i] * pa_y[i];
    }

    // Set up 210-225 components of targeted buffer : HG

    auto tr_x_xzzzz_xxxx = pbuffer.data(idx_dip_hg + 210);

    auto tr_x_xzzzz_xxxy = pbuffer.data(idx_dip_hg + 211);

    auto tr_x_xzzzz_xxxz = pbuffer.data(idx_dip_hg + 212);

    auto tr_x_xzzzz_xxyy = pbuffer.data(idx_dip_hg + 213);

    auto tr_x_xzzzz_xxyz = pbuffer.data(idx_dip_hg + 214);

    auto tr_x_xzzzz_xxzz = pbuffer.data(idx_dip_hg + 215);

    auto tr_x_xzzzz_xyyy = pbuffer.data(idx_dip_hg + 216);

    auto tr_x_xzzzz_xyyz = pbuffer.data(idx_dip_hg + 217);

    auto tr_x_xzzzz_xyzz = pbuffer.data(idx_dip_hg + 218);

    auto tr_x_xzzzz_xzzz = pbuffer.data(idx_dip_hg + 219);

    auto tr_x_xzzzz_yyyy = pbuffer.data(idx_dip_hg + 220);

    auto tr_x_xzzzz_yyyz = pbuffer.data(idx_dip_hg + 221);

    auto tr_x_xzzzz_yyzz = pbuffer.data(idx_dip_hg + 222);

    auto tr_x_xzzzz_yzzz = pbuffer.data(idx_dip_hg + 223);

    auto tr_x_xzzzz_zzzz = pbuffer.data(idx_dip_hg + 224);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             tr_x_xzz_xxxx,   \
                             tr_x_xzz_xxxy,   \
                             tr_x_xzz_xxyy,   \
                             tr_x_xzz_xyyy,   \
                             tr_x_xzzz_xxxx,  \
                             tr_x_xzzz_xxxy,  \
                             tr_x_xzzz_xxyy,  \
                             tr_x_xzzz_xyyy,  \
                             tr_x_xzzzz_xxxx, \
                             tr_x_xzzzz_xxxy, \
                             tr_x_xzzzz_xxxz, \
                             tr_x_xzzzz_xxyy, \
                             tr_x_xzzzz_xxyz, \
                             tr_x_xzzzz_xxzz, \
                             tr_x_xzzzz_xyyy, \
                             tr_x_xzzzz_xyyz, \
                             tr_x_xzzzz_xyzz, \
                             tr_x_xzzzz_xzzz, \
                             tr_x_xzzzz_yyyy, \
                             tr_x_xzzzz_yyyz, \
                             tr_x_xzzzz_yyzz, \
                             tr_x_xzzzz_yzzz, \
                             tr_x_xzzzz_zzzz, \
                             tr_x_zzzz_xxxz,  \
                             tr_x_zzzz_xxyz,  \
                             tr_x_zzzz_xxz,   \
                             tr_x_zzzz_xxzz,  \
                             tr_x_zzzz_xyyz,  \
                             tr_x_zzzz_xyz,   \
                             tr_x_zzzz_xyzz,  \
                             tr_x_zzzz_xzz,   \
                             tr_x_zzzz_xzzz,  \
                             tr_x_zzzz_yyyy,  \
                             tr_x_zzzz_yyyz,  \
                             tr_x_zzzz_yyz,   \
                             tr_x_zzzz_yyzz,  \
                             tr_x_zzzz_yzz,   \
                             tr_x_zzzz_yzzz,  \
                             tr_x_zzzz_zzz,   \
                             tr_x_zzzz_zzzz,  \
                             ts_zzzz_xxxz,    \
                             ts_zzzz_xxyz,    \
                             ts_zzzz_xxzz,    \
                             ts_zzzz_xyyz,    \
                             ts_zzzz_xyzz,    \
                             ts_zzzz_xzzz,    \
                             ts_zzzz_yyyy,    \
                             ts_zzzz_yyyz,    \
                             ts_zzzz_yyzz,    \
                             ts_zzzz_yzzz,    \
                             ts_zzzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzzz_xxxx[i] = 3.0 * tr_x_xzz_xxxx[i] * fe_0 + tr_x_xzzz_xxxx[i] * pa_z[i];

        tr_x_xzzzz_xxxy[i] = 3.0 * tr_x_xzz_xxxy[i] * fe_0 + tr_x_xzzz_xxxy[i] * pa_z[i];

        tr_x_xzzzz_xxxz[i] = 3.0 * tr_x_zzzz_xxz[i] * fe_0 + ts_zzzz_xxxz[i] * fe_0 + tr_x_zzzz_xxxz[i] * pa_x[i];

        tr_x_xzzzz_xxyy[i] = 3.0 * tr_x_xzz_xxyy[i] * fe_0 + tr_x_xzzz_xxyy[i] * pa_z[i];

        tr_x_xzzzz_xxyz[i] = 2.0 * tr_x_zzzz_xyz[i] * fe_0 + ts_zzzz_xxyz[i] * fe_0 + tr_x_zzzz_xxyz[i] * pa_x[i];

        tr_x_xzzzz_xxzz[i] = 2.0 * tr_x_zzzz_xzz[i] * fe_0 + ts_zzzz_xxzz[i] * fe_0 + tr_x_zzzz_xxzz[i] * pa_x[i];

        tr_x_xzzzz_xyyy[i] = 3.0 * tr_x_xzz_xyyy[i] * fe_0 + tr_x_xzzz_xyyy[i] * pa_z[i];

        tr_x_xzzzz_xyyz[i] = tr_x_zzzz_yyz[i] * fe_0 + ts_zzzz_xyyz[i] * fe_0 + tr_x_zzzz_xyyz[i] * pa_x[i];

        tr_x_xzzzz_xyzz[i] = tr_x_zzzz_yzz[i] * fe_0 + ts_zzzz_xyzz[i] * fe_0 + tr_x_zzzz_xyzz[i] * pa_x[i];

        tr_x_xzzzz_xzzz[i] = tr_x_zzzz_zzz[i] * fe_0 + ts_zzzz_xzzz[i] * fe_0 + tr_x_zzzz_xzzz[i] * pa_x[i];

        tr_x_xzzzz_yyyy[i] = ts_zzzz_yyyy[i] * fe_0 + tr_x_zzzz_yyyy[i] * pa_x[i];

        tr_x_xzzzz_yyyz[i] = ts_zzzz_yyyz[i] * fe_0 + tr_x_zzzz_yyyz[i] * pa_x[i];

        tr_x_xzzzz_yyzz[i] = ts_zzzz_yyzz[i] * fe_0 + tr_x_zzzz_yyzz[i] * pa_x[i];

        tr_x_xzzzz_yzzz[i] = ts_zzzz_yzzz[i] * fe_0 + tr_x_zzzz_yzzz[i] * pa_x[i];

        tr_x_xzzzz_zzzz[i] = ts_zzzz_zzzz[i] * fe_0 + tr_x_zzzz_zzzz[i] * pa_x[i];
    }

    // Set up 225-240 components of targeted buffer : HG

    auto tr_x_yyyyy_xxxx = pbuffer.data(idx_dip_hg + 225);

    auto tr_x_yyyyy_xxxy = pbuffer.data(idx_dip_hg + 226);

    auto tr_x_yyyyy_xxxz = pbuffer.data(idx_dip_hg + 227);

    auto tr_x_yyyyy_xxyy = pbuffer.data(idx_dip_hg + 228);

    auto tr_x_yyyyy_xxyz = pbuffer.data(idx_dip_hg + 229);

    auto tr_x_yyyyy_xxzz = pbuffer.data(idx_dip_hg + 230);

    auto tr_x_yyyyy_xyyy = pbuffer.data(idx_dip_hg + 231);

    auto tr_x_yyyyy_xyyz = pbuffer.data(idx_dip_hg + 232);

    auto tr_x_yyyyy_xyzz = pbuffer.data(idx_dip_hg + 233);

    auto tr_x_yyyyy_xzzz = pbuffer.data(idx_dip_hg + 234);

    auto tr_x_yyyyy_yyyy = pbuffer.data(idx_dip_hg + 235);

    auto tr_x_yyyyy_yyyz = pbuffer.data(idx_dip_hg + 236);

    auto tr_x_yyyyy_yyzz = pbuffer.data(idx_dip_hg + 237);

    auto tr_x_yyyyy_yzzz = pbuffer.data(idx_dip_hg + 238);

    auto tr_x_yyyyy_zzzz = pbuffer.data(idx_dip_hg + 239);

#pragma omp simd aligned(pa_y,                \
                             tr_x_yyy_xxxx,   \
                             tr_x_yyy_xxxy,   \
                             tr_x_yyy_xxxz,   \
                             tr_x_yyy_xxyy,   \
                             tr_x_yyy_xxyz,   \
                             tr_x_yyy_xxzz,   \
                             tr_x_yyy_xyyy,   \
                             tr_x_yyy_xyyz,   \
                             tr_x_yyy_xyzz,   \
                             tr_x_yyy_xzzz,   \
                             tr_x_yyy_yyyy,   \
                             tr_x_yyy_yyyz,   \
                             tr_x_yyy_yyzz,   \
                             tr_x_yyy_yzzz,   \
                             tr_x_yyy_zzzz,   \
                             tr_x_yyyy_xxx,   \
                             tr_x_yyyy_xxxx,  \
                             tr_x_yyyy_xxxy,  \
                             tr_x_yyyy_xxxz,  \
                             tr_x_yyyy_xxy,   \
                             tr_x_yyyy_xxyy,  \
                             tr_x_yyyy_xxyz,  \
                             tr_x_yyyy_xxz,   \
                             tr_x_yyyy_xxzz,  \
                             tr_x_yyyy_xyy,   \
                             tr_x_yyyy_xyyy,  \
                             tr_x_yyyy_xyyz,  \
                             tr_x_yyyy_xyz,   \
                             tr_x_yyyy_xyzz,  \
                             tr_x_yyyy_xzz,   \
                             tr_x_yyyy_xzzz,  \
                             tr_x_yyyy_yyy,   \
                             tr_x_yyyy_yyyy,  \
                             tr_x_yyyy_yyyz,  \
                             tr_x_yyyy_yyz,   \
                             tr_x_yyyy_yyzz,  \
                             tr_x_yyyy_yzz,   \
                             tr_x_yyyy_yzzz,  \
                             tr_x_yyyy_zzz,   \
                             tr_x_yyyy_zzzz,  \
                             tr_x_yyyyy_xxxx, \
                             tr_x_yyyyy_xxxy, \
                             tr_x_yyyyy_xxxz, \
                             tr_x_yyyyy_xxyy, \
                             tr_x_yyyyy_xxyz, \
                             tr_x_yyyyy_xxzz, \
                             tr_x_yyyyy_xyyy, \
                             tr_x_yyyyy_xyyz, \
                             tr_x_yyyyy_xyzz, \
                             tr_x_yyyyy_xzzz, \
                             tr_x_yyyyy_yyyy, \
                             tr_x_yyyyy_yyyz, \
                             tr_x_yyyyy_yyzz, \
                             tr_x_yyyyy_yzzz, \
                             tr_x_yyyyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyy_xxxx[i] = 4.0 * tr_x_yyy_xxxx[i] * fe_0 + tr_x_yyyy_xxxx[i] * pa_y[i];

        tr_x_yyyyy_xxxy[i] = 4.0 * tr_x_yyy_xxxy[i] * fe_0 + tr_x_yyyy_xxx[i] * fe_0 + tr_x_yyyy_xxxy[i] * pa_y[i];

        tr_x_yyyyy_xxxz[i] = 4.0 * tr_x_yyy_xxxz[i] * fe_0 + tr_x_yyyy_xxxz[i] * pa_y[i];

        tr_x_yyyyy_xxyy[i] = 4.0 * tr_x_yyy_xxyy[i] * fe_0 + 2.0 * tr_x_yyyy_xxy[i] * fe_0 + tr_x_yyyy_xxyy[i] * pa_y[i];

        tr_x_yyyyy_xxyz[i] = 4.0 * tr_x_yyy_xxyz[i] * fe_0 + tr_x_yyyy_xxz[i] * fe_0 + tr_x_yyyy_xxyz[i] * pa_y[i];

        tr_x_yyyyy_xxzz[i] = 4.0 * tr_x_yyy_xxzz[i] * fe_0 + tr_x_yyyy_xxzz[i] * pa_y[i];

        tr_x_yyyyy_xyyy[i] = 4.0 * tr_x_yyy_xyyy[i] * fe_0 + 3.0 * tr_x_yyyy_xyy[i] * fe_0 + tr_x_yyyy_xyyy[i] * pa_y[i];

        tr_x_yyyyy_xyyz[i] = 4.0 * tr_x_yyy_xyyz[i] * fe_0 + 2.0 * tr_x_yyyy_xyz[i] * fe_0 + tr_x_yyyy_xyyz[i] * pa_y[i];

        tr_x_yyyyy_xyzz[i] = 4.0 * tr_x_yyy_xyzz[i] * fe_0 + tr_x_yyyy_xzz[i] * fe_0 + tr_x_yyyy_xyzz[i] * pa_y[i];

        tr_x_yyyyy_xzzz[i] = 4.0 * tr_x_yyy_xzzz[i] * fe_0 + tr_x_yyyy_xzzz[i] * pa_y[i];

        tr_x_yyyyy_yyyy[i] = 4.0 * tr_x_yyy_yyyy[i] * fe_0 + 4.0 * tr_x_yyyy_yyy[i] * fe_0 + tr_x_yyyy_yyyy[i] * pa_y[i];

        tr_x_yyyyy_yyyz[i] = 4.0 * tr_x_yyy_yyyz[i] * fe_0 + 3.0 * tr_x_yyyy_yyz[i] * fe_0 + tr_x_yyyy_yyyz[i] * pa_y[i];

        tr_x_yyyyy_yyzz[i] = 4.0 * tr_x_yyy_yyzz[i] * fe_0 + 2.0 * tr_x_yyyy_yzz[i] * fe_0 + tr_x_yyyy_yyzz[i] * pa_y[i];

        tr_x_yyyyy_yzzz[i] = 4.0 * tr_x_yyy_yzzz[i] * fe_0 + tr_x_yyyy_zzz[i] * fe_0 + tr_x_yyyy_yzzz[i] * pa_y[i];

        tr_x_yyyyy_zzzz[i] = 4.0 * tr_x_yyy_zzzz[i] * fe_0 + tr_x_yyyy_zzzz[i] * pa_y[i];
    }

    // Set up 240-255 components of targeted buffer : HG

    auto tr_x_yyyyz_xxxx = pbuffer.data(idx_dip_hg + 240);

    auto tr_x_yyyyz_xxxy = pbuffer.data(idx_dip_hg + 241);

    auto tr_x_yyyyz_xxxz = pbuffer.data(idx_dip_hg + 242);

    auto tr_x_yyyyz_xxyy = pbuffer.data(idx_dip_hg + 243);

    auto tr_x_yyyyz_xxyz = pbuffer.data(idx_dip_hg + 244);

    auto tr_x_yyyyz_xxzz = pbuffer.data(idx_dip_hg + 245);

    auto tr_x_yyyyz_xyyy = pbuffer.data(idx_dip_hg + 246);

    auto tr_x_yyyyz_xyyz = pbuffer.data(idx_dip_hg + 247);

    auto tr_x_yyyyz_xyzz = pbuffer.data(idx_dip_hg + 248);

    auto tr_x_yyyyz_xzzz = pbuffer.data(idx_dip_hg + 249);

    auto tr_x_yyyyz_yyyy = pbuffer.data(idx_dip_hg + 250);

    auto tr_x_yyyyz_yyyz = pbuffer.data(idx_dip_hg + 251);

    auto tr_x_yyyyz_yyzz = pbuffer.data(idx_dip_hg + 252);

    auto tr_x_yyyyz_yzzz = pbuffer.data(idx_dip_hg + 253);

    auto tr_x_yyyyz_zzzz = pbuffer.data(idx_dip_hg + 254);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             tr_x_yyyy_xxxx,  \
                             tr_x_yyyy_xxxy,  \
                             tr_x_yyyy_xxy,   \
                             tr_x_yyyy_xxyy,  \
                             tr_x_yyyy_xxyz,  \
                             tr_x_yyyy_xyy,   \
                             tr_x_yyyy_xyyy,  \
                             tr_x_yyyy_xyyz,  \
                             tr_x_yyyy_xyz,   \
                             tr_x_yyyy_xyzz,  \
                             tr_x_yyyy_yyy,   \
                             tr_x_yyyy_yyyy,  \
                             tr_x_yyyy_yyyz,  \
                             tr_x_yyyy_yyz,   \
                             tr_x_yyyy_yyzz,  \
                             tr_x_yyyy_yzz,   \
                             tr_x_yyyy_yzzz,  \
                             tr_x_yyyyz_xxxx, \
                             tr_x_yyyyz_xxxy, \
                             tr_x_yyyyz_xxxz, \
                             tr_x_yyyyz_xxyy, \
                             tr_x_yyyyz_xxyz, \
                             tr_x_yyyyz_xxzz, \
                             tr_x_yyyyz_xyyy, \
                             tr_x_yyyyz_xyyz, \
                             tr_x_yyyyz_xyzz, \
                             tr_x_yyyyz_xzzz, \
                             tr_x_yyyyz_yyyy, \
                             tr_x_yyyyz_yyyz, \
                             tr_x_yyyyz_yyzz, \
                             tr_x_yyyyz_yzzz, \
                             tr_x_yyyyz_zzzz, \
                             tr_x_yyyz_xxxz,  \
                             tr_x_yyyz_xxzz,  \
                             tr_x_yyyz_xzzz,  \
                             tr_x_yyyz_zzzz,  \
                             tr_x_yyz_xxxz,   \
                             tr_x_yyz_xxzz,   \
                             tr_x_yyz_xzzz,   \
                             tr_x_yyz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyz_xxxx[i] = tr_x_yyyy_xxxx[i] * pa_z[i];

        tr_x_yyyyz_xxxy[i] = tr_x_yyyy_xxxy[i] * pa_z[i];

        tr_x_yyyyz_xxxz[i] = 3.0 * tr_x_yyz_xxxz[i] * fe_0 + tr_x_yyyz_xxxz[i] * pa_y[i];

        tr_x_yyyyz_xxyy[i] = tr_x_yyyy_xxyy[i] * pa_z[i];

        tr_x_yyyyz_xxyz[i] = tr_x_yyyy_xxy[i] * fe_0 + tr_x_yyyy_xxyz[i] * pa_z[i];

        tr_x_yyyyz_xxzz[i] = 3.0 * tr_x_yyz_xxzz[i] * fe_0 + tr_x_yyyz_xxzz[i] * pa_y[i];

        tr_x_yyyyz_xyyy[i] = tr_x_yyyy_xyyy[i] * pa_z[i];

        tr_x_yyyyz_xyyz[i] = tr_x_yyyy_xyy[i] * fe_0 + tr_x_yyyy_xyyz[i] * pa_z[i];

        tr_x_yyyyz_xyzz[i] = 2.0 * tr_x_yyyy_xyz[i] * fe_0 + tr_x_yyyy_xyzz[i] * pa_z[i];

        tr_x_yyyyz_xzzz[i] = 3.0 * tr_x_yyz_xzzz[i] * fe_0 + tr_x_yyyz_xzzz[i] * pa_y[i];

        tr_x_yyyyz_yyyy[i] = tr_x_yyyy_yyyy[i] * pa_z[i];

        tr_x_yyyyz_yyyz[i] = tr_x_yyyy_yyy[i] * fe_0 + tr_x_yyyy_yyyz[i] * pa_z[i];

        tr_x_yyyyz_yyzz[i] = 2.0 * tr_x_yyyy_yyz[i] * fe_0 + tr_x_yyyy_yyzz[i] * pa_z[i];

        tr_x_yyyyz_yzzz[i] = 3.0 * tr_x_yyyy_yzz[i] * fe_0 + tr_x_yyyy_yzzz[i] * pa_z[i];

        tr_x_yyyyz_zzzz[i] = 3.0 * tr_x_yyz_zzzz[i] * fe_0 + tr_x_yyyz_zzzz[i] * pa_y[i];
    }

    // Set up 255-270 components of targeted buffer : HG

    auto tr_x_yyyzz_xxxx = pbuffer.data(idx_dip_hg + 255);

    auto tr_x_yyyzz_xxxy = pbuffer.data(idx_dip_hg + 256);

    auto tr_x_yyyzz_xxxz = pbuffer.data(idx_dip_hg + 257);

    auto tr_x_yyyzz_xxyy = pbuffer.data(idx_dip_hg + 258);

    auto tr_x_yyyzz_xxyz = pbuffer.data(idx_dip_hg + 259);

    auto tr_x_yyyzz_xxzz = pbuffer.data(idx_dip_hg + 260);

    auto tr_x_yyyzz_xyyy = pbuffer.data(idx_dip_hg + 261);

    auto tr_x_yyyzz_xyyz = pbuffer.data(idx_dip_hg + 262);

    auto tr_x_yyyzz_xyzz = pbuffer.data(idx_dip_hg + 263);

    auto tr_x_yyyzz_xzzz = pbuffer.data(idx_dip_hg + 264);

    auto tr_x_yyyzz_yyyy = pbuffer.data(idx_dip_hg + 265);

    auto tr_x_yyyzz_yyyz = pbuffer.data(idx_dip_hg + 266);

    auto tr_x_yyyzz_yyzz = pbuffer.data(idx_dip_hg + 267);

    auto tr_x_yyyzz_yzzz = pbuffer.data(idx_dip_hg + 268);

    auto tr_x_yyyzz_zzzz = pbuffer.data(idx_dip_hg + 269);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             tr_x_yyy_xxxy,   \
                             tr_x_yyy_xxyy,   \
                             tr_x_yyy_xyyy,   \
                             tr_x_yyy_yyyy,   \
                             tr_x_yyyz_xxxy,  \
                             tr_x_yyyz_xxyy,  \
                             tr_x_yyyz_xyyy,  \
                             tr_x_yyyz_yyyy,  \
                             tr_x_yyyzz_xxxx, \
                             tr_x_yyyzz_xxxy, \
                             tr_x_yyyzz_xxxz, \
                             tr_x_yyyzz_xxyy, \
                             tr_x_yyyzz_xxyz, \
                             tr_x_yyyzz_xxzz, \
                             tr_x_yyyzz_xyyy, \
                             tr_x_yyyzz_xyyz, \
                             tr_x_yyyzz_xyzz, \
                             tr_x_yyyzz_xzzz, \
                             tr_x_yyyzz_yyyy, \
                             tr_x_yyyzz_yyyz, \
                             tr_x_yyyzz_yyzz, \
                             tr_x_yyyzz_yzzz, \
                             tr_x_yyyzz_zzzz, \
                             tr_x_yyzz_xxxx,  \
                             tr_x_yyzz_xxxz,  \
                             tr_x_yyzz_xxyz,  \
                             tr_x_yyzz_xxz,   \
                             tr_x_yyzz_xxzz,  \
                             tr_x_yyzz_xyyz,  \
                             tr_x_yyzz_xyz,   \
                             tr_x_yyzz_xyzz,  \
                             tr_x_yyzz_xzz,   \
                             tr_x_yyzz_xzzz,  \
                             tr_x_yyzz_yyyz,  \
                             tr_x_yyzz_yyz,   \
                             tr_x_yyzz_yyzz,  \
                             tr_x_yyzz_yzz,   \
                             tr_x_yyzz_yzzz,  \
                             tr_x_yyzz_zzz,   \
                             tr_x_yyzz_zzzz,  \
                             tr_x_yzz_xxxx,   \
                             tr_x_yzz_xxxz,   \
                             tr_x_yzz_xxyz,   \
                             tr_x_yzz_xxzz,   \
                             tr_x_yzz_xyyz,   \
                             tr_x_yzz_xyzz,   \
                             tr_x_yzz_xzzz,   \
                             tr_x_yzz_yyyz,   \
                             tr_x_yzz_yyzz,   \
                             tr_x_yzz_yzzz,   \
                             tr_x_yzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyzz_xxxx[i] = 2.0 * tr_x_yzz_xxxx[i] * fe_0 + tr_x_yyzz_xxxx[i] * pa_y[i];

        tr_x_yyyzz_xxxy[i] = tr_x_yyy_xxxy[i] * fe_0 + tr_x_yyyz_xxxy[i] * pa_z[i];

        tr_x_yyyzz_xxxz[i] = 2.0 * tr_x_yzz_xxxz[i] * fe_0 + tr_x_yyzz_xxxz[i] * pa_y[i];

        tr_x_yyyzz_xxyy[i] = tr_x_yyy_xxyy[i] * fe_0 + tr_x_yyyz_xxyy[i] * pa_z[i];

        tr_x_yyyzz_xxyz[i] = 2.0 * tr_x_yzz_xxyz[i] * fe_0 + tr_x_yyzz_xxz[i] * fe_0 + tr_x_yyzz_xxyz[i] * pa_y[i];

        tr_x_yyyzz_xxzz[i] = 2.0 * tr_x_yzz_xxzz[i] * fe_0 + tr_x_yyzz_xxzz[i] * pa_y[i];

        tr_x_yyyzz_xyyy[i] = tr_x_yyy_xyyy[i] * fe_0 + tr_x_yyyz_xyyy[i] * pa_z[i];

        tr_x_yyyzz_xyyz[i] = 2.0 * tr_x_yzz_xyyz[i] * fe_0 + 2.0 * tr_x_yyzz_xyz[i] * fe_0 + tr_x_yyzz_xyyz[i] * pa_y[i];

        tr_x_yyyzz_xyzz[i] = 2.0 * tr_x_yzz_xyzz[i] * fe_0 + tr_x_yyzz_xzz[i] * fe_0 + tr_x_yyzz_xyzz[i] * pa_y[i];

        tr_x_yyyzz_xzzz[i] = 2.0 * tr_x_yzz_xzzz[i] * fe_0 + tr_x_yyzz_xzzz[i] * pa_y[i];

        tr_x_yyyzz_yyyy[i] = tr_x_yyy_yyyy[i] * fe_0 + tr_x_yyyz_yyyy[i] * pa_z[i];

        tr_x_yyyzz_yyyz[i] = 2.0 * tr_x_yzz_yyyz[i] * fe_0 + 3.0 * tr_x_yyzz_yyz[i] * fe_0 + tr_x_yyzz_yyyz[i] * pa_y[i];

        tr_x_yyyzz_yyzz[i] = 2.0 * tr_x_yzz_yyzz[i] * fe_0 + 2.0 * tr_x_yyzz_yzz[i] * fe_0 + tr_x_yyzz_yyzz[i] * pa_y[i];

        tr_x_yyyzz_yzzz[i] = 2.0 * tr_x_yzz_yzzz[i] * fe_0 + tr_x_yyzz_zzz[i] * fe_0 + tr_x_yyzz_yzzz[i] * pa_y[i];

        tr_x_yyyzz_zzzz[i] = 2.0 * tr_x_yzz_zzzz[i] * fe_0 + tr_x_yyzz_zzzz[i] * pa_y[i];
    }

    // Set up 270-285 components of targeted buffer : HG

    auto tr_x_yyzzz_xxxx = pbuffer.data(idx_dip_hg + 270);

    auto tr_x_yyzzz_xxxy = pbuffer.data(idx_dip_hg + 271);

    auto tr_x_yyzzz_xxxz = pbuffer.data(idx_dip_hg + 272);

    auto tr_x_yyzzz_xxyy = pbuffer.data(idx_dip_hg + 273);

    auto tr_x_yyzzz_xxyz = pbuffer.data(idx_dip_hg + 274);

    auto tr_x_yyzzz_xxzz = pbuffer.data(idx_dip_hg + 275);

    auto tr_x_yyzzz_xyyy = pbuffer.data(idx_dip_hg + 276);

    auto tr_x_yyzzz_xyyz = pbuffer.data(idx_dip_hg + 277);

    auto tr_x_yyzzz_xyzz = pbuffer.data(idx_dip_hg + 278);

    auto tr_x_yyzzz_xzzz = pbuffer.data(idx_dip_hg + 279);

    auto tr_x_yyzzz_yyyy = pbuffer.data(idx_dip_hg + 280);

    auto tr_x_yyzzz_yyyz = pbuffer.data(idx_dip_hg + 281);

    auto tr_x_yyzzz_yyzz = pbuffer.data(idx_dip_hg + 282);

    auto tr_x_yyzzz_yzzz = pbuffer.data(idx_dip_hg + 283);

    auto tr_x_yyzzz_zzzz = pbuffer.data(idx_dip_hg + 284);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             tr_x_yyz_xxxy,   \
                             tr_x_yyz_xxyy,   \
                             tr_x_yyz_xyyy,   \
                             tr_x_yyz_yyyy,   \
                             tr_x_yyzz_xxxy,  \
                             tr_x_yyzz_xxyy,  \
                             tr_x_yyzz_xyyy,  \
                             tr_x_yyzz_yyyy,  \
                             tr_x_yyzzz_xxxx, \
                             tr_x_yyzzz_xxxy, \
                             tr_x_yyzzz_xxxz, \
                             tr_x_yyzzz_xxyy, \
                             tr_x_yyzzz_xxyz, \
                             tr_x_yyzzz_xxzz, \
                             tr_x_yyzzz_xyyy, \
                             tr_x_yyzzz_xyyz, \
                             tr_x_yyzzz_xyzz, \
                             tr_x_yyzzz_xzzz, \
                             tr_x_yyzzz_yyyy, \
                             tr_x_yyzzz_yyyz, \
                             tr_x_yyzzz_yyzz, \
                             tr_x_yyzzz_yzzz, \
                             tr_x_yyzzz_zzzz, \
                             tr_x_yzzz_xxxx,  \
                             tr_x_yzzz_xxxz,  \
                             tr_x_yzzz_xxyz,  \
                             tr_x_yzzz_xxz,   \
                             tr_x_yzzz_xxzz,  \
                             tr_x_yzzz_xyyz,  \
                             tr_x_yzzz_xyz,   \
                             tr_x_yzzz_xyzz,  \
                             tr_x_yzzz_xzz,   \
                             tr_x_yzzz_xzzz,  \
                             tr_x_yzzz_yyyz,  \
                             tr_x_yzzz_yyz,   \
                             tr_x_yzzz_yyzz,  \
                             tr_x_yzzz_yzz,   \
                             tr_x_yzzz_yzzz,  \
                             tr_x_yzzz_zzz,   \
                             tr_x_yzzz_zzzz,  \
                             tr_x_zzz_xxxx,   \
                             tr_x_zzz_xxxz,   \
                             tr_x_zzz_xxyz,   \
                             tr_x_zzz_xxzz,   \
                             tr_x_zzz_xyyz,   \
                             tr_x_zzz_xyzz,   \
                             tr_x_zzz_xzzz,   \
                             tr_x_zzz_yyyz,   \
                             tr_x_zzz_yyzz,   \
                             tr_x_zzz_yzzz,   \
                             tr_x_zzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzzz_xxxx[i] = tr_x_zzz_xxxx[i] * fe_0 + tr_x_yzzz_xxxx[i] * pa_y[i];

        tr_x_yyzzz_xxxy[i] = 2.0 * tr_x_yyz_xxxy[i] * fe_0 + tr_x_yyzz_xxxy[i] * pa_z[i];

        tr_x_yyzzz_xxxz[i] = tr_x_zzz_xxxz[i] * fe_0 + tr_x_yzzz_xxxz[i] * pa_y[i];

        tr_x_yyzzz_xxyy[i] = 2.0 * tr_x_yyz_xxyy[i] * fe_0 + tr_x_yyzz_xxyy[i] * pa_z[i];

        tr_x_yyzzz_xxyz[i] = tr_x_zzz_xxyz[i] * fe_0 + tr_x_yzzz_xxz[i] * fe_0 + tr_x_yzzz_xxyz[i] * pa_y[i];

        tr_x_yyzzz_xxzz[i] = tr_x_zzz_xxzz[i] * fe_0 + tr_x_yzzz_xxzz[i] * pa_y[i];

        tr_x_yyzzz_xyyy[i] = 2.0 * tr_x_yyz_xyyy[i] * fe_0 + tr_x_yyzz_xyyy[i] * pa_z[i];

        tr_x_yyzzz_xyyz[i] = tr_x_zzz_xyyz[i] * fe_0 + 2.0 * tr_x_yzzz_xyz[i] * fe_0 + tr_x_yzzz_xyyz[i] * pa_y[i];

        tr_x_yyzzz_xyzz[i] = tr_x_zzz_xyzz[i] * fe_0 + tr_x_yzzz_xzz[i] * fe_0 + tr_x_yzzz_xyzz[i] * pa_y[i];

        tr_x_yyzzz_xzzz[i] = tr_x_zzz_xzzz[i] * fe_0 + tr_x_yzzz_xzzz[i] * pa_y[i];

        tr_x_yyzzz_yyyy[i] = 2.0 * tr_x_yyz_yyyy[i] * fe_0 + tr_x_yyzz_yyyy[i] * pa_z[i];

        tr_x_yyzzz_yyyz[i] = tr_x_zzz_yyyz[i] * fe_0 + 3.0 * tr_x_yzzz_yyz[i] * fe_0 + tr_x_yzzz_yyyz[i] * pa_y[i];

        tr_x_yyzzz_yyzz[i] = tr_x_zzz_yyzz[i] * fe_0 + 2.0 * tr_x_yzzz_yzz[i] * fe_0 + tr_x_yzzz_yyzz[i] * pa_y[i];

        tr_x_yyzzz_yzzz[i] = tr_x_zzz_yzzz[i] * fe_0 + tr_x_yzzz_zzz[i] * fe_0 + tr_x_yzzz_yzzz[i] * pa_y[i];

        tr_x_yyzzz_zzzz[i] = tr_x_zzz_zzzz[i] * fe_0 + tr_x_yzzz_zzzz[i] * pa_y[i];
    }

    // Set up 285-300 components of targeted buffer : HG

    auto tr_x_yzzzz_xxxx = pbuffer.data(idx_dip_hg + 285);

    auto tr_x_yzzzz_xxxy = pbuffer.data(idx_dip_hg + 286);

    auto tr_x_yzzzz_xxxz = pbuffer.data(idx_dip_hg + 287);

    auto tr_x_yzzzz_xxyy = pbuffer.data(idx_dip_hg + 288);

    auto tr_x_yzzzz_xxyz = pbuffer.data(idx_dip_hg + 289);

    auto tr_x_yzzzz_xxzz = pbuffer.data(idx_dip_hg + 290);

    auto tr_x_yzzzz_xyyy = pbuffer.data(idx_dip_hg + 291);

    auto tr_x_yzzzz_xyyz = pbuffer.data(idx_dip_hg + 292);

    auto tr_x_yzzzz_xyzz = pbuffer.data(idx_dip_hg + 293);

    auto tr_x_yzzzz_xzzz = pbuffer.data(idx_dip_hg + 294);

    auto tr_x_yzzzz_yyyy = pbuffer.data(idx_dip_hg + 295);

    auto tr_x_yzzzz_yyyz = pbuffer.data(idx_dip_hg + 296);

    auto tr_x_yzzzz_yyzz = pbuffer.data(idx_dip_hg + 297);

    auto tr_x_yzzzz_yzzz = pbuffer.data(idx_dip_hg + 298);

    auto tr_x_yzzzz_zzzz = pbuffer.data(idx_dip_hg + 299);

#pragma omp simd aligned(pa_y,                \
                             tr_x_yzzzz_xxxx, \
                             tr_x_yzzzz_xxxy, \
                             tr_x_yzzzz_xxxz, \
                             tr_x_yzzzz_xxyy, \
                             tr_x_yzzzz_xxyz, \
                             tr_x_yzzzz_xxzz, \
                             tr_x_yzzzz_xyyy, \
                             tr_x_yzzzz_xyyz, \
                             tr_x_yzzzz_xyzz, \
                             tr_x_yzzzz_xzzz, \
                             tr_x_yzzzz_yyyy, \
                             tr_x_yzzzz_yyyz, \
                             tr_x_yzzzz_yyzz, \
                             tr_x_yzzzz_yzzz, \
                             tr_x_yzzzz_zzzz, \
                             tr_x_zzzz_xxx,   \
                             tr_x_zzzz_xxxx,  \
                             tr_x_zzzz_xxxy,  \
                             tr_x_zzzz_xxxz,  \
                             tr_x_zzzz_xxy,   \
                             tr_x_zzzz_xxyy,  \
                             tr_x_zzzz_xxyz,  \
                             tr_x_zzzz_xxz,   \
                             tr_x_zzzz_xxzz,  \
                             tr_x_zzzz_xyy,   \
                             tr_x_zzzz_xyyy,  \
                             tr_x_zzzz_xyyz,  \
                             tr_x_zzzz_xyz,   \
                             tr_x_zzzz_xyzz,  \
                             tr_x_zzzz_xzz,   \
                             tr_x_zzzz_xzzz,  \
                             tr_x_zzzz_yyy,   \
                             tr_x_zzzz_yyyy,  \
                             tr_x_zzzz_yyyz,  \
                             tr_x_zzzz_yyz,   \
                             tr_x_zzzz_yyzz,  \
                             tr_x_zzzz_yzz,   \
                             tr_x_zzzz_yzzz,  \
                             tr_x_zzzz_zzz,   \
                             tr_x_zzzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzzz_xxxx[i] = tr_x_zzzz_xxxx[i] * pa_y[i];

        tr_x_yzzzz_xxxy[i] = tr_x_zzzz_xxx[i] * fe_0 + tr_x_zzzz_xxxy[i] * pa_y[i];

        tr_x_yzzzz_xxxz[i] = tr_x_zzzz_xxxz[i] * pa_y[i];

        tr_x_yzzzz_xxyy[i] = 2.0 * tr_x_zzzz_xxy[i] * fe_0 + tr_x_zzzz_xxyy[i] * pa_y[i];

        tr_x_yzzzz_xxyz[i] = tr_x_zzzz_xxz[i] * fe_0 + tr_x_zzzz_xxyz[i] * pa_y[i];

        tr_x_yzzzz_xxzz[i] = tr_x_zzzz_xxzz[i] * pa_y[i];

        tr_x_yzzzz_xyyy[i] = 3.0 * tr_x_zzzz_xyy[i] * fe_0 + tr_x_zzzz_xyyy[i] * pa_y[i];

        tr_x_yzzzz_xyyz[i] = 2.0 * tr_x_zzzz_xyz[i] * fe_0 + tr_x_zzzz_xyyz[i] * pa_y[i];

        tr_x_yzzzz_xyzz[i] = tr_x_zzzz_xzz[i] * fe_0 + tr_x_zzzz_xyzz[i] * pa_y[i];

        tr_x_yzzzz_xzzz[i] = tr_x_zzzz_xzzz[i] * pa_y[i];

        tr_x_yzzzz_yyyy[i] = 4.0 * tr_x_zzzz_yyy[i] * fe_0 + tr_x_zzzz_yyyy[i] * pa_y[i];

        tr_x_yzzzz_yyyz[i] = 3.0 * tr_x_zzzz_yyz[i] * fe_0 + tr_x_zzzz_yyyz[i] * pa_y[i];

        tr_x_yzzzz_yyzz[i] = 2.0 * tr_x_zzzz_yzz[i] * fe_0 + tr_x_zzzz_yyzz[i] * pa_y[i];

        tr_x_yzzzz_yzzz[i] = tr_x_zzzz_zzz[i] * fe_0 + tr_x_zzzz_yzzz[i] * pa_y[i];

        tr_x_yzzzz_zzzz[i] = tr_x_zzzz_zzzz[i] * pa_y[i];
    }

    // Set up 300-315 components of targeted buffer : HG

    auto tr_x_zzzzz_xxxx = pbuffer.data(idx_dip_hg + 300);

    auto tr_x_zzzzz_xxxy = pbuffer.data(idx_dip_hg + 301);

    auto tr_x_zzzzz_xxxz = pbuffer.data(idx_dip_hg + 302);

    auto tr_x_zzzzz_xxyy = pbuffer.data(idx_dip_hg + 303);

    auto tr_x_zzzzz_xxyz = pbuffer.data(idx_dip_hg + 304);

    auto tr_x_zzzzz_xxzz = pbuffer.data(idx_dip_hg + 305);

    auto tr_x_zzzzz_xyyy = pbuffer.data(idx_dip_hg + 306);

    auto tr_x_zzzzz_xyyz = pbuffer.data(idx_dip_hg + 307);

    auto tr_x_zzzzz_xyzz = pbuffer.data(idx_dip_hg + 308);

    auto tr_x_zzzzz_xzzz = pbuffer.data(idx_dip_hg + 309);

    auto tr_x_zzzzz_yyyy = pbuffer.data(idx_dip_hg + 310);

    auto tr_x_zzzzz_yyyz = pbuffer.data(idx_dip_hg + 311);

    auto tr_x_zzzzz_yyzz = pbuffer.data(idx_dip_hg + 312);

    auto tr_x_zzzzz_yzzz = pbuffer.data(idx_dip_hg + 313);

    auto tr_x_zzzzz_zzzz = pbuffer.data(idx_dip_hg + 314);

#pragma omp simd aligned(pa_z,                \
                             tr_x_zzz_xxxx,   \
                             tr_x_zzz_xxxy,   \
                             tr_x_zzz_xxxz,   \
                             tr_x_zzz_xxyy,   \
                             tr_x_zzz_xxyz,   \
                             tr_x_zzz_xxzz,   \
                             tr_x_zzz_xyyy,   \
                             tr_x_zzz_xyyz,   \
                             tr_x_zzz_xyzz,   \
                             tr_x_zzz_xzzz,   \
                             tr_x_zzz_yyyy,   \
                             tr_x_zzz_yyyz,   \
                             tr_x_zzz_yyzz,   \
                             tr_x_zzz_yzzz,   \
                             tr_x_zzz_zzzz,   \
                             tr_x_zzzz_xxx,   \
                             tr_x_zzzz_xxxx,  \
                             tr_x_zzzz_xxxy,  \
                             tr_x_zzzz_xxxz,  \
                             tr_x_zzzz_xxy,   \
                             tr_x_zzzz_xxyy,  \
                             tr_x_zzzz_xxyz,  \
                             tr_x_zzzz_xxz,   \
                             tr_x_zzzz_xxzz,  \
                             tr_x_zzzz_xyy,   \
                             tr_x_zzzz_xyyy,  \
                             tr_x_zzzz_xyyz,  \
                             tr_x_zzzz_xyz,   \
                             tr_x_zzzz_xyzz,  \
                             tr_x_zzzz_xzz,   \
                             tr_x_zzzz_xzzz,  \
                             tr_x_zzzz_yyy,   \
                             tr_x_zzzz_yyyy,  \
                             tr_x_zzzz_yyyz,  \
                             tr_x_zzzz_yyz,   \
                             tr_x_zzzz_yyzz,  \
                             tr_x_zzzz_yzz,   \
                             tr_x_zzzz_yzzz,  \
                             tr_x_zzzz_zzz,   \
                             tr_x_zzzz_zzzz,  \
                             tr_x_zzzzz_xxxx, \
                             tr_x_zzzzz_xxxy, \
                             tr_x_zzzzz_xxxz, \
                             tr_x_zzzzz_xxyy, \
                             tr_x_zzzzz_xxyz, \
                             tr_x_zzzzz_xxzz, \
                             tr_x_zzzzz_xyyy, \
                             tr_x_zzzzz_xyyz, \
                             tr_x_zzzzz_xyzz, \
                             tr_x_zzzzz_xzzz, \
                             tr_x_zzzzz_yyyy, \
                             tr_x_zzzzz_yyyz, \
                             tr_x_zzzzz_yyzz, \
                             tr_x_zzzzz_yzzz, \
                             tr_x_zzzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzzz_xxxx[i] = 4.0 * tr_x_zzz_xxxx[i] * fe_0 + tr_x_zzzz_xxxx[i] * pa_z[i];

        tr_x_zzzzz_xxxy[i] = 4.0 * tr_x_zzz_xxxy[i] * fe_0 + tr_x_zzzz_xxxy[i] * pa_z[i];

        tr_x_zzzzz_xxxz[i] = 4.0 * tr_x_zzz_xxxz[i] * fe_0 + tr_x_zzzz_xxx[i] * fe_0 + tr_x_zzzz_xxxz[i] * pa_z[i];

        tr_x_zzzzz_xxyy[i] = 4.0 * tr_x_zzz_xxyy[i] * fe_0 + tr_x_zzzz_xxyy[i] * pa_z[i];

        tr_x_zzzzz_xxyz[i] = 4.0 * tr_x_zzz_xxyz[i] * fe_0 + tr_x_zzzz_xxy[i] * fe_0 + tr_x_zzzz_xxyz[i] * pa_z[i];

        tr_x_zzzzz_xxzz[i] = 4.0 * tr_x_zzz_xxzz[i] * fe_0 + 2.0 * tr_x_zzzz_xxz[i] * fe_0 + tr_x_zzzz_xxzz[i] * pa_z[i];

        tr_x_zzzzz_xyyy[i] = 4.0 * tr_x_zzz_xyyy[i] * fe_0 + tr_x_zzzz_xyyy[i] * pa_z[i];

        tr_x_zzzzz_xyyz[i] = 4.0 * tr_x_zzz_xyyz[i] * fe_0 + tr_x_zzzz_xyy[i] * fe_0 + tr_x_zzzz_xyyz[i] * pa_z[i];

        tr_x_zzzzz_xyzz[i] = 4.0 * tr_x_zzz_xyzz[i] * fe_0 + 2.0 * tr_x_zzzz_xyz[i] * fe_0 + tr_x_zzzz_xyzz[i] * pa_z[i];

        tr_x_zzzzz_xzzz[i] = 4.0 * tr_x_zzz_xzzz[i] * fe_0 + 3.0 * tr_x_zzzz_xzz[i] * fe_0 + tr_x_zzzz_xzzz[i] * pa_z[i];

        tr_x_zzzzz_yyyy[i] = 4.0 * tr_x_zzz_yyyy[i] * fe_0 + tr_x_zzzz_yyyy[i] * pa_z[i];

        tr_x_zzzzz_yyyz[i] = 4.0 * tr_x_zzz_yyyz[i] * fe_0 + tr_x_zzzz_yyy[i] * fe_0 + tr_x_zzzz_yyyz[i] * pa_z[i];

        tr_x_zzzzz_yyzz[i] = 4.0 * tr_x_zzz_yyzz[i] * fe_0 + 2.0 * tr_x_zzzz_yyz[i] * fe_0 + tr_x_zzzz_yyzz[i] * pa_z[i];

        tr_x_zzzzz_yzzz[i] = 4.0 * tr_x_zzz_yzzz[i] * fe_0 + 3.0 * tr_x_zzzz_yzz[i] * fe_0 + tr_x_zzzz_yzzz[i] * pa_z[i];

        tr_x_zzzzz_zzzz[i] = 4.0 * tr_x_zzz_zzzz[i] * fe_0 + 4.0 * tr_x_zzzz_zzz[i] * fe_0 + tr_x_zzzz_zzzz[i] * pa_z[i];
    }

    // Set up 315-330 components of targeted buffer : HG

    auto tr_y_xxxxx_xxxx = pbuffer.data(idx_dip_hg + 315);

    auto tr_y_xxxxx_xxxy = pbuffer.data(idx_dip_hg + 316);

    auto tr_y_xxxxx_xxxz = pbuffer.data(idx_dip_hg + 317);

    auto tr_y_xxxxx_xxyy = pbuffer.data(idx_dip_hg + 318);

    auto tr_y_xxxxx_xxyz = pbuffer.data(idx_dip_hg + 319);

    auto tr_y_xxxxx_xxzz = pbuffer.data(idx_dip_hg + 320);

    auto tr_y_xxxxx_xyyy = pbuffer.data(idx_dip_hg + 321);

    auto tr_y_xxxxx_xyyz = pbuffer.data(idx_dip_hg + 322);

    auto tr_y_xxxxx_xyzz = pbuffer.data(idx_dip_hg + 323);

    auto tr_y_xxxxx_xzzz = pbuffer.data(idx_dip_hg + 324);

    auto tr_y_xxxxx_yyyy = pbuffer.data(idx_dip_hg + 325);

    auto tr_y_xxxxx_yyyz = pbuffer.data(idx_dip_hg + 326);

    auto tr_y_xxxxx_yyzz = pbuffer.data(idx_dip_hg + 327);

    auto tr_y_xxxxx_yzzz = pbuffer.data(idx_dip_hg + 328);

    auto tr_y_xxxxx_zzzz = pbuffer.data(idx_dip_hg + 329);

#pragma omp simd aligned(pa_x,                \
                             tr_y_xxx_xxxx,   \
                             tr_y_xxx_xxxy,   \
                             tr_y_xxx_xxxz,   \
                             tr_y_xxx_xxyy,   \
                             tr_y_xxx_xxyz,   \
                             tr_y_xxx_xxzz,   \
                             tr_y_xxx_xyyy,   \
                             tr_y_xxx_xyyz,   \
                             tr_y_xxx_xyzz,   \
                             tr_y_xxx_xzzz,   \
                             tr_y_xxx_yyyy,   \
                             tr_y_xxx_yyyz,   \
                             tr_y_xxx_yyzz,   \
                             tr_y_xxx_yzzz,   \
                             tr_y_xxx_zzzz,   \
                             tr_y_xxxx_xxx,   \
                             tr_y_xxxx_xxxx,  \
                             tr_y_xxxx_xxxy,  \
                             tr_y_xxxx_xxxz,  \
                             tr_y_xxxx_xxy,   \
                             tr_y_xxxx_xxyy,  \
                             tr_y_xxxx_xxyz,  \
                             tr_y_xxxx_xxz,   \
                             tr_y_xxxx_xxzz,  \
                             tr_y_xxxx_xyy,   \
                             tr_y_xxxx_xyyy,  \
                             tr_y_xxxx_xyyz,  \
                             tr_y_xxxx_xyz,   \
                             tr_y_xxxx_xyzz,  \
                             tr_y_xxxx_xzz,   \
                             tr_y_xxxx_xzzz,  \
                             tr_y_xxxx_yyy,   \
                             tr_y_xxxx_yyyy,  \
                             tr_y_xxxx_yyyz,  \
                             tr_y_xxxx_yyz,   \
                             tr_y_xxxx_yyzz,  \
                             tr_y_xxxx_yzz,   \
                             tr_y_xxxx_yzzz,  \
                             tr_y_xxxx_zzz,   \
                             tr_y_xxxx_zzzz,  \
                             tr_y_xxxxx_xxxx, \
                             tr_y_xxxxx_xxxy, \
                             tr_y_xxxxx_xxxz, \
                             tr_y_xxxxx_xxyy, \
                             tr_y_xxxxx_xxyz, \
                             tr_y_xxxxx_xxzz, \
                             tr_y_xxxxx_xyyy, \
                             tr_y_xxxxx_xyyz, \
                             tr_y_xxxxx_xyzz, \
                             tr_y_xxxxx_xzzz, \
                             tr_y_xxxxx_yyyy, \
                             tr_y_xxxxx_yyyz, \
                             tr_y_xxxxx_yyzz, \
                             tr_y_xxxxx_yzzz, \
                             tr_y_xxxxx_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxx_xxxx[i] = 4.0 * tr_y_xxx_xxxx[i] * fe_0 + 4.0 * tr_y_xxxx_xxx[i] * fe_0 + tr_y_xxxx_xxxx[i] * pa_x[i];

        tr_y_xxxxx_xxxy[i] = 4.0 * tr_y_xxx_xxxy[i] * fe_0 + 3.0 * tr_y_xxxx_xxy[i] * fe_0 + tr_y_xxxx_xxxy[i] * pa_x[i];

        tr_y_xxxxx_xxxz[i] = 4.0 * tr_y_xxx_xxxz[i] * fe_0 + 3.0 * tr_y_xxxx_xxz[i] * fe_0 + tr_y_xxxx_xxxz[i] * pa_x[i];

        tr_y_xxxxx_xxyy[i] = 4.0 * tr_y_xxx_xxyy[i] * fe_0 + 2.0 * tr_y_xxxx_xyy[i] * fe_0 + tr_y_xxxx_xxyy[i] * pa_x[i];

        tr_y_xxxxx_xxyz[i] = 4.0 * tr_y_xxx_xxyz[i] * fe_0 + 2.0 * tr_y_xxxx_xyz[i] * fe_0 + tr_y_xxxx_xxyz[i] * pa_x[i];

        tr_y_xxxxx_xxzz[i] = 4.0 * tr_y_xxx_xxzz[i] * fe_0 + 2.0 * tr_y_xxxx_xzz[i] * fe_0 + tr_y_xxxx_xxzz[i] * pa_x[i];

        tr_y_xxxxx_xyyy[i] = 4.0 * tr_y_xxx_xyyy[i] * fe_0 + tr_y_xxxx_yyy[i] * fe_0 + tr_y_xxxx_xyyy[i] * pa_x[i];

        tr_y_xxxxx_xyyz[i] = 4.0 * tr_y_xxx_xyyz[i] * fe_0 + tr_y_xxxx_yyz[i] * fe_0 + tr_y_xxxx_xyyz[i] * pa_x[i];

        tr_y_xxxxx_xyzz[i] = 4.0 * tr_y_xxx_xyzz[i] * fe_0 + tr_y_xxxx_yzz[i] * fe_0 + tr_y_xxxx_xyzz[i] * pa_x[i];

        tr_y_xxxxx_xzzz[i] = 4.0 * tr_y_xxx_xzzz[i] * fe_0 + tr_y_xxxx_zzz[i] * fe_0 + tr_y_xxxx_xzzz[i] * pa_x[i];

        tr_y_xxxxx_yyyy[i] = 4.0 * tr_y_xxx_yyyy[i] * fe_0 + tr_y_xxxx_yyyy[i] * pa_x[i];

        tr_y_xxxxx_yyyz[i] = 4.0 * tr_y_xxx_yyyz[i] * fe_0 + tr_y_xxxx_yyyz[i] * pa_x[i];

        tr_y_xxxxx_yyzz[i] = 4.0 * tr_y_xxx_yyzz[i] * fe_0 + tr_y_xxxx_yyzz[i] * pa_x[i];

        tr_y_xxxxx_yzzz[i] = 4.0 * tr_y_xxx_yzzz[i] * fe_0 + tr_y_xxxx_yzzz[i] * pa_x[i];

        tr_y_xxxxx_zzzz[i] = 4.0 * tr_y_xxx_zzzz[i] * fe_0 + tr_y_xxxx_zzzz[i] * pa_x[i];
    }

    // Set up 330-345 components of targeted buffer : HG

    auto tr_y_xxxxy_xxxx = pbuffer.data(idx_dip_hg + 330);

    auto tr_y_xxxxy_xxxy = pbuffer.data(idx_dip_hg + 331);

    auto tr_y_xxxxy_xxxz = pbuffer.data(idx_dip_hg + 332);

    auto tr_y_xxxxy_xxyy = pbuffer.data(idx_dip_hg + 333);

    auto tr_y_xxxxy_xxyz = pbuffer.data(idx_dip_hg + 334);

    auto tr_y_xxxxy_xxzz = pbuffer.data(idx_dip_hg + 335);

    auto tr_y_xxxxy_xyyy = pbuffer.data(idx_dip_hg + 336);

    auto tr_y_xxxxy_xyyz = pbuffer.data(idx_dip_hg + 337);

    auto tr_y_xxxxy_xyzz = pbuffer.data(idx_dip_hg + 338);

    auto tr_y_xxxxy_xzzz = pbuffer.data(idx_dip_hg + 339);

    auto tr_y_xxxxy_yyyy = pbuffer.data(idx_dip_hg + 340);

    auto tr_y_xxxxy_yyyz = pbuffer.data(idx_dip_hg + 341);

    auto tr_y_xxxxy_yyzz = pbuffer.data(idx_dip_hg + 342);

    auto tr_y_xxxxy_yzzz = pbuffer.data(idx_dip_hg + 343);

    auto tr_y_xxxxy_zzzz = pbuffer.data(idx_dip_hg + 344);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tr_y_xxxx_xxxx,  \
                             tr_y_xxxx_xxxz,  \
                             tr_y_xxxx_xxzz,  \
                             tr_y_xxxx_xzzz,  \
                             tr_y_xxxxy_xxxx, \
                             tr_y_xxxxy_xxxy, \
                             tr_y_xxxxy_xxxz, \
                             tr_y_xxxxy_xxyy, \
                             tr_y_xxxxy_xxyz, \
                             tr_y_xxxxy_xxzz, \
                             tr_y_xxxxy_xyyy, \
                             tr_y_xxxxy_xyyz, \
                             tr_y_xxxxy_xyzz, \
                             tr_y_xxxxy_xzzz, \
                             tr_y_xxxxy_yyyy, \
                             tr_y_xxxxy_yyyz, \
                             tr_y_xxxxy_yyzz, \
                             tr_y_xxxxy_yzzz, \
                             tr_y_xxxxy_zzzz, \
                             tr_y_xxxy_xxxy,  \
                             tr_y_xxxy_xxy,   \
                             tr_y_xxxy_xxyy,  \
                             tr_y_xxxy_xxyz,  \
                             tr_y_xxxy_xyy,   \
                             tr_y_xxxy_xyyy,  \
                             tr_y_xxxy_xyyz,  \
                             tr_y_xxxy_xyz,   \
                             tr_y_xxxy_xyzz,  \
                             tr_y_xxxy_yyy,   \
                             tr_y_xxxy_yyyy,  \
                             tr_y_xxxy_yyyz,  \
                             tr_y_xxxy_yyz,   \
                             tr_y_xxxy_yyzz,  \
                             tr_y_xxxy_yzz,   \
                             tr_y_xxxy_yzzz,  \
                             tr_y_xxxy_zzzz,  \
                             tr_y_xxy_xxxy,   \
                             tr_y_xxy_xxyy,   \
                             tr_y_xxy_xxyz,   \
                             tr_y_xxy_xyyy,   \
                             tr_y_xxy_xyyz,   \
                             tr_y_xxy_xyzz,   \
                             tr_y_xxy_yyyy,   \
                             tr_y_xxy_yyyz,   \
                             tr_y_xxy_yyzz,   \
                             tr_y_xxy_yzzz,   \
                             tr_y_xxy_zzzz,   \
                             ts_xxxx_xxxx,    \
                             ts_xxxx_xxxz,    \
                             ts_xxxx_xxzz,    \
                             ts_xxxx_xzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxy_xxxx[i] = ts_xxxx_xxxx[i] * fe_0 + tr_y_xxxx_xxxx[i] * pa_y[i];

        tr_y_xxxxy_xxxy[i] = 3.0 * tr_y_xxy_xxxy[i] * fe_0 + 3.0 * tr_y_xxxy_xxy[i] * fe_0 + tr_y_xxxy_xxxy[i] * pa_x[i];

        tr_y_xxxxy_xxxz[i] = ts_xxxx_xxxz[i] * fe_0 + tr_y_xxxx_xxxz[i] * pa_y[i];

        tr_y_xxxxy_xxyy[i] = 3.0 * tr_y_xxy_xxyy[i] * fe_0 + 2.0 * tr_y_xxxy_xyy[i] * fe_0 + tr_y_xxxy_xxyy[i] * pa_x[i];

        tr_y_xxxxy_xxyz[i] = 3.0 * tr_y_xxy_xxyz[i] * fe_0 + 2.0 * tr_y_xxxy_xyz[i] * fe_0 + tr_y_xxxy_xxyz[i] * pa_x[i];

        tr_y_xxxxy_xxzz[i] = ts_xxxx_xxzz[i] * fe_0 + tr_y_xxxx_xxzz[i] * pa_y[i];

        tr_y_xxxxy_xyyy[i] = 3.0 * tr_y_xxy_xyyy[i] * fe_0 + tr_y_xxxy_yyy[i] * fe_0 + tr_y_xxxy_xyyy[i] * pa_x[i];

        tr_y_xxxxy_xyyz[i] = 3.0 * tr_y_xxy_xyyz[i] * fe_0 + tr_y_xxxy_yyz[i] * fe_0 + tr_y_xxxy_xyyz[i] * pa_x[i];

        tr_y_xxxxy_xyzz[i] = 3.0 * tr_y_xxy_xyzz[i] * fe_0 + tr_y_xxxy_yzz[i] * fe_0 + tr_y_xxxy_xyzz[i] * pa_x[i];

        tr_y_xxxxy_xzzz[i] = ts_xxxx_xzzz[i] * fe_0 + tr_y_xxxx_xzzz[i] * pa_y[i];

        tr_y_xxxxy_yyyy[i] = 3.0 * tr_y_xxy_yyyy[i] * fe_0 + tr_y_xxxy_yyyy[i] * pa_x[i];

        tr_y_xxxxy_yyyz[i] = 3.0 * tr_y_xxy_yyyz[i] * fe_0 + tr_y_xxxy_yyyz[i] * pa_x[i];

        tr_y_xxxxy_yyzz[i] = 3.0 * tr_y_xxy_yyzz[i] * fe_0 + tr_y_xxxy_yyzz[i] * pa_x[i];

        tr_y_xxxxy_yzzz[i] = 3.0 * tr_y_xxy_yzzz[i] * fe_0 + tr_y_xxxy_yzzz[i] * pa_x[i];

        tr_y_xxxxy_zzzz[i] = 3.0 * tr_y_xxy_zzzz[i] * fe_0 + tr_y_xxxy_zzzz[i] * pa_x[i];
    }

    // Set up 345-360 components of targeted buffer : HG

    auto tr_y_xxxxz_xxxx = pbuffer.data(idx_dip_hg + 345);

    auto tr_y_xxxxz_xxxy = pbuffer.data(idx_dip_hg + 346);

    auto tr_y_xxxxz_xxxz = pbuffer.data(idx_dip_hg + 347);

    auto tr_y_xxxxz_xxyy = pbuffer.data(idx_dip_hg + 348);

    auto tr_y_xxxxz_xxyz = pbuffer.data(idx_dip_hg + 349);

    auto tr_y_xxxxz_xxzz = pbuffer.data(idx_dip_hg + 350);

    auto tr_y_xxxxz_xyyy = pbuffer.data(idx_dip_hg + 351);

    auto tr_y_xxxxz_xyyz = pbuffer.data(idx_dip_hg + 352);

    auto tr_y_xxxxz_xyzz = pbuffer.data(idx_dip_hg + 353);

    auto tr_y_xxxxz_xzzz = pbuffer.data(idx_dip_hg + 354);

    auto tr_y_xxxxz_yyyy = pbuffer.data(idx_dip_hg + 355);

    auto tr_y_xxxxz_yyyz = pbuffer.data(idx_dip_hg + 356);

    auto tr_y_xxxxz_yyzz = pbuffer.data(idx_dip_hg + 357);

    auto tr_y_xxxxz_yzzz = pbuffer.data(idx_dip_hg + 358);

    auto tr_y_xxxxz_zzzz = pbuffer.data(idx_dip_hg + 359);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             tr_y_xxxx_xxx,   \
                             tr_y_xxxx_xxxx,  \
                             tr_y_xxxx_xxxy,  \
                             tr_y_xxxx_xxxz,  \
                             tr_y_xxxx_xxy,   \
                             tr_y_xxxx_xxyy,  \
                             tr_y_xxxx_xxyz,  \
                             tr_y_xxxx_xxz,   \
                             tr_y_xxxx_xxzz,  \
                             tr_y_xxxx_xyy,   \
                             tr_y_xxxx_xyyy,  \
                             tr_y_xxxx_xyyz,  \
                             tr_y_xxxx_xyz,   \
                             tr_y_xxxx_xyzz,  \
                             tr_y_xxxx_xzz,   \
                             tr_y_xxxx_xzzz,  \
                             tr_y_xxxx_yyyy,  \
                             tr_y_xxxxz_xxxx, \
                             tr_y_xxxxz_xxxy, \
                             tr_y_xxxxz_xxxz, \
                             tr_y_xxxxz_xxyy, \
                             tr_y_xxxxz_xxyz, \
                             tr_y_xxxxz_xxzz, \
                             tr_y_xxxxz_xyyy, \
                             tr_y_xxxxz_xyyz, \
                             tr_y_xxxxz_xyzz, \
                             tr_y_xxxxz_xzzz, \
                             tr_y_xxxxz_yyyy, \
                             tr_y_xxxxz_yyyz, \
                             tr_y_xxxxz_yyzz, \
                             tr_y_xxxxz_yzzz, \
                             tr_y_xxxxz_zzzz, \
                             tr_y_xxxz_yyyz,  \
                             tr_y_xxxz_yyzz,  \
                             tr_y_xxxz_yzzz,  \
                             tr_y_xxxz_zzzz,  \
                             tr_y_xxz_yyyz,   \
                             tr_y_xxz_yyzz,   \
                             tr_y_xxz_yzzz,   \
                             tr_y_xxz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxz_xxxx[i] = tr_y_xxxx_xxxx[i] * pa_z[i];

        tr_y_xxxxz_xxxy[i] = tr_y_xxxx_xxxy[i] * pa_z[i];

        tr_y_xxxxz_xxxz[i] = tr_y_xxxx_xxx[i] * fe_0 + tr_y_xxxx_xxxz[i] * pa_z[i];

        tr_y_xxxxz_xxyy[i] = tr_y_xxxx_xxyy[i] * pa_z[i];

        tr_y_xxxxz_xxyz[i] = tr_y_xxxx_xxy[i] * fe_0 + tr_y_xxxx_xxyz[i] * pa_z[i];

        tr_y_xxxxz_xxzz[i] = 2.0 * tr_y_xxxx_xxz[i] * fe_0 + tr_y_xxxx_xxzz[i] * pa_z[i];

        tr_y_xxxxz_xyyy[i] = tr_y_xxxx_xyyy[i] * pa_z[i];

        tr_y_xxxxz_xyyz[i] = tr_y_xxxx_xyy[i] * fe_0 + tr_y_xxxx_xyyz[i] * pa_z[i];

        tr_y_xxxxz_xyzz[i] = 2.0 * tr_y_xxxx_xyz[i] * fe_0 + tr_y_xxxx_xyzz[i] * pa_z[i];

        tr_y_xxxxz_xzzz[i] = 3.0 * tr_y_xxxx_xzz[i] * fe_0 + tr_y_xxxx_xzzz[i] * pa_z[i];

        tr_y_xxxxz_yyyy[i] = tr_y_xxxx_yyyy[i] * pa_z[i];

        tr_y_xxxxz_yyyz[i] = 3.0 * tr_y_xxz_yyyz[i] * fe_0 + tr_y_xxxz_yyyz[i] * pa_x[i];

        tr_y_xxxxz_yyzz[i] = 3.0 * tr_y_xxz_yyzz[i] * fe_0 + tr_y_xxxz_yyzz[i] * pa_x[i];

        tr_y_xxxxz_yzzz[i] = 3.0 * tr_y_xxz_yzzz[i] * fe_0 + tr_y_xxxz_yzzz[i] * pa_x[i];

        tr_y_xxxxz_zzzz[i] = 3.0 * tr_y_xxz_zzzz[i] * fe_0 + tr_y_xxxz_zzzz[i] * pa_x[i];
    }

    // Set up 360-375 components of targeted buffer : HG

    auto tr_y_xxxyy_xxxx = pbuffer.data(idx_dip_hg + 360);

    auto tr_y_xxxyy_xxxy = pbuffer.data(idx_dip_hg + 361);

    auto tr_y_xxxyy_xxxz = pbuffer.data(idx_dip_hg + 362);

    auto tr_y_xxxyy_xxyy = pbuffer.data(idx_dip_hg + 363);

    auto tr_y_xxxyy_xxyz = pbuffer.data(idx_dip_hg + 364);

    auto tr_y_xxxyy_xxzz = pbuffer.data(idx_dip_hg + 365);

    auto tr_y_xxxyy_xyyy = pbuffer.data(idx_dip_hg + 366);

    auto tr_y_xxxyy_xyyz = pbuffer.data(idx_dip_hg + 367);

    auto tr_y_xxxyy_xyzz = pbuffer.data(idx_dip_hg + 368);

    auto tr_y_xxxyy_xzzz = pbuffer.data(idx_dip_hg + 369);

    auto tr_y_xxxyy_yyyy = pbuffer.data(idx_dip_hg + 370);

    auto tr_y_xxxyy_yyyz = pbuffer.data(idx_dip_hg + 371);

    auto tr_y_xxxyy_yyzz = pbuffer.data(idx_dip_hg + 372);

    auto tr_y_xxxyy_yzzz = pbuffer.data(idx_dip_hg + 373);

    auto tr_y_xxxyy_zzzz = pbuffer.data(idx_dip_hg + 374);

#pragma omp simd aligned(pa_x,                \
                             tr_y_xxxyy_xxxx, \
                             tr_y_xxxyy_xxxy, \
                             tr_y_xxxyy_xxxz, \
                             tr_y_xxxyy_xxyy, \
                             tr_y_xxxyy_xxyz, \
                             tr_y_xxxyy_xxzz, \
                             tr_y_xxxyy_xyyy, \
                             tr_y_xxxyy_xyyz, \
                             tr_y_xxxyy_xyzz, \
                             tr_y_xxxyy_xzzz, \
                             tr_y_xxxyy_yyyy, \
                             tr_y_xxxyy_yyyz, \
                             tr_y_xxxyy_yyzz, \
                             tr_y_xxxyy_yzzz, \
                             tr_y_xxxyy_zzzz, \
                             tr_y_xxyy_xxx,   \
                             tr_y_xxyy_xxxx,  \
                             tr_y_xxyy_xxxy,  \
                             tr_y_xxyy_xxxz,  \
                             tr_y_xxyy_xxy,   \
                             tr_y_xxyy_xxyy,  \
                             tr_y_xxyy_xxyz,  \
                             tr_y_xxyy_xxz,   \
                             tr_y_xxyy_xxzz,  \
                             tr_y_xxyy_xyy,   \
                             tr_y_xxyy_xyyy,  \
                             tr_y_xxyy_xyyz,  \
                             tr_y_xxyy_xyz,   \
                             tr_y_xxyy_xyzz,  \
                             tr_y_xxyy_xzz,   \
                             tr_y_xxyy_xzzz,  \
                             tr_y_xxyy_yyy,   \
                             tr_y_xxyy_yyyy,  \
                             tr_y_xxyy_yyyz,  \
                             tr_y_xxyy_yyz,   \
                             tr_y_xxyy_yyzz,  \
                             tr_y_xxyy_yzz,   \
                             tr_y_xxyy_yzzz,  \
                             tr_y_xxyy_zzz,   \
                             tr_y_xxyy_zzzz,  \
                             tr_y_xyy_xxxx,   \
                             tr_y_xyy_xxxy,   \
                             tr_y_xyy_xxxz,   \
                             tr_y_xyy_xxyy,   \
                             tr_y_xyy_xxyz,   \
                             tr_y_xyy_xxzz,   \
                             tr_y_xyy_xyyy,   \
                             tr_y_xyy_xyyz,   \
                             tr_y_xyy_xyzz,   \
                             tr_y_xyy_xzzz,   \
                             tr_y_xyy_yyyy,   \
                             tr_y_xyy_yyyz,   \
                             tr_y_xyy_yyzz,   \
                             tr_y_xyy_yzzz,   \
                             tr_y_xyy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyy_xxxx[i] = 2.0 * tr_y_xyy_xxxx[i] * fe_0 + 4.0 * tr_y_xxyy_xxx[i] * fe_0 + tr_y_xxyy_xxxx[i] * pa_x[i];

        tr_y_xxxyy_xxxy[i] = 2.0 * tr_y_xyy_xxxy[i] * fe_0 + 3.0 * tr_y_xxyy_xxy[i] * fe_0 + tr_y_xxyy_xxxy[i] * pa_x[i];

        tr_y_xxxyy_xxxz[i] = 2.0 * tr_y_xyy_xxxz[i] * fe_0 + 3.0 * tr_y_xxyy_xxz[i] * fe_0 + tr_y_xxyy_xxxz[i] * pa_x[i];

        tr_y_xxxyy_xxyy[i] = 2.0 * tr_y_xyy_xxyy[i] * fe_0 + 2.0 * tr_y_xxyy_xyy[i] * fe_0 + tr_y_xxyy_xxyy[i] * pa_x[i];

        tr_y_xxxyy_xxyz[i] = 2.0 * tr_y_xyy_xxyz[i] * fe_0 + 2.0 * tr_y_xxyy_xyz[i] * fe_0 + tr_y_xxyy_xxyz[i] * pa_x[i];

        tr_y_xxxyy_xxzz[i] = 2.0 * tr_y_xyy_xxzz[i] * fe_0 + 2.0 * tr_y_xxyy_xzz[i] * fe_0 + tr_y_xxyy_xxzz[i] * pa_x[i];

        tr_y_xxxyy_xyyy[i] = 2.0 * tr_y_xyy_xyyy[i] * fe_0 + tr_y_xxyy_yyy[i] * fe_0 + tr_y_xxyy_xyyy[i] * pa_x[i];

        tr_y_xxxyy_xyyz[i] = 2.0 * tr_y_xyy_xyyz[i] * fe_0 + tr_y_xxyy_yyz[i] * fe_0 + tr_y_xxyy_xyyz[i] * pa_x[i];

        tr_y_xxxyy_xyzz[i] = 2.0 * tr_y_xyy_xyzz[i] * fe_0 + tr_y_xxyy_yzz[i] * fe_0 + tr_y_xxyy_xyzz[i] * pa_x[i];

        tr_y_xxxyy_xzzz[i] = 2.0 * tr_y_xyy_xzzz[i] * fe_0 + tr_y_xxyy_zzz[i] * fe_0 + tr_y_xxyy_xzzz[i] * pa_x[i];

        tr_y_xxxyy_yyyy[i] = 2.0 * tr_y_xyy_yyyy[i] * fe_0 + tr_y_xxyy_yyyy[i] * pa_x[i];

        tr_y_xxxyy_yyyz[i] = 2.0 * tr_y_xyy_yyyz[i] * fe_0 + tr_y_xxyy_yyyz[i] * pa_x[i];

        tr_y_xxxyy_yyzz[i] = 2.0 * tr_y_xyy_yyzz[i] * fe_0 + tr_y_xxyy_yyzz[i] * pa_x[i];

        tr_y_xxxyy_yzzz[i] = 2.0 * tr_y_xyy_yzzz[i] * fe_0 + tr_y_xxyy_yzzz[i] * pa_x[i];

        tr_y_xxxyy_zzzz[i] = 2.0 * tr_y_xyy_zzzz[i] * fe_0 + tr_y_xxyy_zzzz[i] * pa_x[i];
    }

    // Set up 375-390 components of targeted buffer : HG

    auto tr_y_xxxyz_xxxx = pbuffer.data(idx_dip_hg + 375);

    auto tr_y_xxxyz_xxxy = pbuffer.data(idx_dip_hg + 376);

    auto tr_y_xxxyz_xxxz = pbuffer.data(idx_dip_hg + 377);

    auto tr_y_xxxyz_xxyy = pbuffer.data(idx_dip_hg + 378);

    auto tr_y_xxxyz_xxyz = pbuffer.data(idx_dip_hg + 379);

    auto tr_y_xxxyz_xxzz = pbuffer.data(idx_dip_hg + 380);

    auto tr_y_xxxyz_xyyy = pbuffer.data(idx_dip_hg + 381);

    auto tr_y_xxxyz_xyyz = pbuffer.data(idx_dip_hg + 382);

    auto tr_y_xxxyz_xyzz = pbuffer.data(idx_dip_hg + 383);

    auto tr_y_xxxyz_xzzz = pbuffer.data(idx_dip_hg + 384);

    auto tr_y_xxxyz_yyyy = pbuffer.data(idx_dip_hg + 385);

    auto tr_y_xxxyz_yyyz = pbuffer.data(idx_dip_hg + 386);

    auto tr_y_xxxyz_yyzz = pbuffer.data(idx_dip_hg + 387);

    auto tr_y_xxxyz_yzzz = pbuffer.data(idx_dip_hg + 388);

    auto tr_y_xxxyz_zzzz = pbuffer.data(idx_dip_hg + 389);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             tr_y_xxxy_xxxx,  \
                             tr_y_xxxy_xxxy,  \
                             tr_y_xxxy_xxy,   \
                             tr_y_xxxy_xxyy,  \
                             tr_y_xxxy_xxyz,  \
                             tr_y_xxxy_xyy,   \
                             tr_y_xxxy_xyyy,  \
                             tr_y_xxxy_xyyz,  \
                             tr_y_xxxy_xyz,   \
                             tr_y_xxxy_xyzz,  \
                             tr_y_xxxy_yyyy,  \
                             tr_y_xxxyz_xxxx, \
                             tr_y_xxxyz_xxxy, \
                             tr_y_xxxyz_xxxz, \
                             tr_y_xxxyz_xxyy, \
                             tr_y_xxxyz_xxyz, \
                             tr_y_xxxyz_xxzz, \
                             tr_y_xxxyz_xyyy, \
                             tr_y_xxxyz_xyyz, \
                             tr_y_xxxyz_xyzz, \
                             tr_y_xxxyz_xzzz, \
                             tr_y_xxxyz_yyyy, \
                             tr_y_xxxyz_yyyz, \
                             tr_y_xxxyz_yyzz, \
                             tr_y_xxxyz_yzzz, \
                             tr_y_xxxyz_zzzz, \
                             tr_y_xxxz_xxxz,  \
                             tr_y_xxxz_xxzz,  \
                             tr_y_xxxz_xzzz,  \
                             tr_y_xxyz_yyyz,  \
                             tr_y_xxyz_yyzz,  \
                             tr_y_xxyz_yzzz,  \
                             tr_y_xxyz_zzzz,  \
                             tr_y_xyz_yyyz,   \
                             tr_y_xyz_yyzz,   \
                             tr_y_xyz_yzzz,   \
                             tr_y_xyz_zzzz,   \
                             ts_xxxz_xxxz,    \
                             ts_xxxz_xxzz,    \
                             ts_xxxz_xzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyz_xxxx[i] = tr_y_xxxy_xxxx[i] * pa_z[i];

        tr_y_xxxyz_xxxy[i] = tr_y_xxxy_xxxy[i] * pa_z[i];

        tr_y_xxxyz_xxxz[i] = ts_xxxz_xxxz[i] * fe_0 + tr_y_xxxz_xxxz[i] * pa_y[i];

        tr_y_xxxyz_xxyy[i] = tr_y_xxxy_xxyy[i] * pa_z[i];

        tr_y_xxxyz_xxyz[i] = tr_y_xxxy_xxy[i] * fe_0 + tr_y_xxxy_xxyz[i] * pa_z[i];

        tr_y_xxxyz_xxzz[i] = ts_xxxz_xxzz[i] * fe_0 + tr_y_xxxz_xxzz[i] * pa_y[i];

        tr_y_xxxyz_xyyy[i] = tr_y_xxxy_xyyy[i] * pa_z[i];

        tr_y_xxxyz_xyyz[i] = tr_y_xxxy_xyy[i] * fe_0 + tr_y_xxxy_xyyz[i] * pa_z[i];

        tr_y_xxxyz_xyzz[i] = 2.0 * tr_y_xxxy_xyz[i] * fe_0 + tr_y_xxxy_xyzz[i] * pa_z[i];

        tr_y_xxxyz_xzzz[i] = ts_xxxz_xzzz[i] * fe_0 + tr_y_xxxz_xzzz[i] * pa_y[i];

        tr_y_xxxyz_yyyy[i] = tr_y_xxxy_yyyy[i] * pa_z[i];

        tr_y_xxxyz_yyyz[i] = 2.0 * tr_y_xyz_yyyz[i] * fe_0 + tr_y_xxyz_yyyz[i] * pa_x[i];

        tr_y_xxxyz_yyzz[i] = 2.0 * tr_y_xyz_yyzz[i] * fe_0 + tr_y_xxyz_yyzz[i] * pa_x[i];

        tr_y_xxxyz_yzzz[i] = 2.0 * tr_y_xyz_yzzz[i] * fe_0 + tr_y_xxyz_yzzz[i] * pa_x[i];

        tr_y_xxxyz_zzzz[i] = 2.0 * tr_y_xyz_zzzz[i] * fe_0 + tr_y_xxyz_zzzz[i] * pa_x[i];
    }

    // Set up 390-405 components of targeted buffer : HG

    auto tr_y_xxxzz_xxxx = pbuffer.data(idx_dip_hg + 390);

    auto tr_y_xxxzz_xxxy = pbuffer.data(idx_dip_hg + 391);

    auto tr_y_xxxzz_xxxz = pbuffer.data(idx_dip_hg + 392);

    auto tr_y_xxxzz_xxyy = pbuffer.data(idx_dip_hg + 393);

    auto tr_y_xxxzz_xxyz = pbuffer.data(idx_dip_hg + 394);

    auto tr_y_xxxzz_xxzz = pbuffer.data(idx_dip_hg + 395);

    auto tr_y_xxxzz_xyyy = pbuffer.data(idx_dip_hg + 396);

    auto tr_y_xxxzz_xyyz = pbuffer.data(idx_dip_hg + 397);

    auto tr_y_xxxzz_xyzz = pbuffer.data(idx_dip_hg + 398);

    auto tr_y_xxxzz_xzzz = pbuffer.data(idx_dip_hg + 399);

    auto tr_y_xxxzz_yyyy = pbuffer.data(idx_dip_hg + 400);

    auto tr_y_xxxzz_yyyz = pbuffer.data(idx_dip_hg + 401);

    auto tr_y_xxxzz_yyzz = pbuffer.data(idx_dip_hg + 402);

    auto tr_y_xxxzz_yzzz = pbuffer.data(idx_dip_hg + 403);

    auto tr_y_xxxzz_zzzz = pbuffer.data(idx_dip_hg + 404);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             tr_y_xxx_xxxx,   \
                             tr_y_xxx_xxxy,   \
                             tr_y_xxx_xxyy,   \
                             tr_y_xxx_xyyy,   \
                             tr_y_xxxz_xxxx,  \
                             tr_y_xxxz_xxxy,  \
                             tr_y_xxxz_xxyy,  \
                             tr_y_xxxz_xyyy,  \
                             tr_y_xxxzz_xxxx, \
                             tr_y_xxxzz_xxxy, \
                             tr_y_xxxzz_xxxz, \
                             tr_y_xxxzz_xxyy, \
                             tr_y_xxxzz_xxyz, \
                             tr_y_xxxzz_xxzz, \
                             tr_y_xxxzz_xyyy, \
                             tr_y_xxxzz_xyyz, \
                             tr_y_xxxzz_xyzz, \
                             tr_y_xxxzz_xzzz, \
                             tr_y_xxxzz_yyyy, \
                             tr_y_xxxzz_yyyz, \
                             tr_y_xxxzz_yyzz, \
                             tr_y_xxxzz_yzzz, \
                             tr_y_xxxzz_zzzz, \
                             tr_y_xxzz_xxxz,  \
                             tr_y_xxzz_xxyz,  \
                             tr_y_xxzz_xxz,   \
                             tr_y_xxzz_xxzz,  \
                             tr_y_xxzz_xyyz,  \
                             tr_y_xxzz_xyz,   \
                             tr_y_xxzz_xyzz,  \
                             tr_y_xxzz_xzz,   \
                             tr_y_xxzz_xzzz,  \
                             tr_y_xxzz_yyyy,  \
                             tr_y_xxzz_yyyz,  \
                             tr_y_xxzz_yyz,   \
                             tr_y_xxzz_yyzz,  \
                             tr_y_xxzz_yzz,   \
                             tr_y_xxzz_yzzz,  \
                             tr_y_xxzz_zzz,   \
                             tr_y_xxzz_zzzz,  \
                             tr_y_xzz_xxxz,   \
                             tr_y_xzz_xxyz,   \
                             tr_y_xzz_xxzz,   \
                             tr_y_xzz_xyyz,   \
                             tr_y_xzz_xyzz,   \
                             tr_y_xzz_xzzz,   \
                             tr_y_xzz_yyyy,   \
                             tr_y_xzz_yyyz,   \
                             tr_y_xzz_yyzz,   \
                             tr_y_xzz_yzzz,   \
                             tr_y_xzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxzz_xxxx[i] = tr_y_xxx_xxxx[i] * fe_0 + tr_y_xxxz_xxxx[i] * pa_z[i];

        tr_y_xxxzz_xxxy[i] = tr_y_xxx_xxxy[i] * fe_0 + tr_y_xxxz_xxxy[i] * pa_z[i];

        tr_y_xxxzz_xxxz[i] = 2.0 * tr_y_xzz_xxxz[i] * fe_0 + 3.0 * tr_y_xxzz_xxz[i] * fe_0 + tr_y_xxzz_xxxz[i] * pa_x[i];

        tr_y_xxxzz_xxyy[i] = tr_y_xxx_xxyy[i] * fe_0 + tr_y_xxxz_xxyy[i] * pa_z[i];

        tr_y_xxxzz_xxyz[i] = 2.0 * tr_y_xzz_xxyz[i] * fe_0 + 2.0 * tr_y_xxzz_xyz[i] * fe_0 + tr_y_xxzz_xxyz[i] * pa_x[i];

        tr_y_xxxzz_xxzz[i] = 2.0 * tr_y_xzz_xxzz[i] * fe_0 + 2.0 * tr_y_xxzz_xzz[i] * fe_0 + tr_y_xxzz_xxzz[i] * pa_x[i];

        tr_y_xxxzz_xyyy[i] = tr_y_xxx_xyyy[i] * fe_0 + tr_y_xxxz_xyyy[i] * pa_z[i];

        tr_y_xxxzz_xyyz[i] = 2.0 * tr_y_xzz_xyyz[i] * fe_0 + tr_y_xxzz_yyz[i] * fe_0 + tr_y_xxzz_xyyz[i] * pa_x[i];

        tr_y_xxxzz_xyzz[i] = 2.0 * tr_y_xzz_xyzz[i] * fe_0 + tr_y_xxzz_yzz[i] * fe_0 + tr_y_xxzz_xyzz[i] * pa_x[i];

        tr_y_xxxzz_xzzz[i] = 2.0 * tr_y_xzz_xzzz[i] * fe_0 + tr_y_xxzz_zzz[i] * fe_0 + tr_y_xxzz_xzzz[i] * pa_x[i];

        tr_y_xxxzz_yyyy[i] = 2.0 * tr_y_xzz_yyyy[i] * fe_0 + tr_y_xxzz_yyyy[i] * pa_x[i];

        tr_y_xxxzz_yyyz[i] = 2.0 * tr_y_xzz_yyyz[i] * fe_0 + tr_y_xxzz_yyyz[i] * pa_x[i];

        tr_y_xxxzz_yyzz[i] = 2.0 * tr_y_xzz_yyzz[i] * fe_0 + tr_y_xxzz_yyzz[i] * pa_x[i];

        tr_y_xxxzz_yzzz[i] = 2.0 * tr_y_xzz_yzzz[i] * fe_0 + tr_y_xxzz_yzzz[i] * pa_x[i];

        tr_y_xxxzz_zzzz[i] = 2.0 * tr_y_xzz_zzzz[i] * fe_0 + tr_y_xxzz_zzzz[i] * pa_x[i];
    }

    // Set up 405-420 components of targeted buffer : HG

    auto tr_y_xxyyy_xxxx = pbuffer.data(idx_dip_hg + 405);

    auto tr_y_xxyyy_xxxy = pbuffer.data(idx_dip_hg + 406);

    auto tr_y_xxyyy_xxxz = pbuffer.data(idx_dip_hg + 407);

    auto tr_y_xxyyy_xxyy = pbuffer.data(idx_dip_hg + 408);

    auto tr_y_xxyyy_xxyz = pbuffer.data(idx_dip_hg + 409);

    auto tr_y_xxyyy_xxzz = pbuffer.data(idx_dip_hg + 410);

    auto tr_y_xxyyy_xyyy = pbuffer.data(idx_dip_hg + 411);

    auto tr_y_xxyyy_xyyz = pbuffer.data(idx_dip_hg + 412);

    auto tr_y_xxyyy_xyzz = pbuffer.data(idx_dip_hg + 413);

    auto tr_y_xxyyy_xzzz = pbuffer.data(idx_dip_hg + 414);

    auto tr_y_xxyyy_yyyy = pbuffer.data(idx_dip_hg + 415);

    auto tr_y_xxyyy_yyyz = pbuffer.data(idx_dip_hg + 416);

    auto tr_y_xxyyy_yyzz = pbuffer.data(idx_dip_hg + 417);

    auto tr_y_xxyyy_yzzz = pbuffer.data(idx_dip_hg + 418);

    auto tr_y_xxyyy_zzzz = pbuffer.data(idx_dip_hg + 419);

#pragma omp simd aligned(pa_x,                \
                             tr_y_xxyyy_xxxx, \
                             tr_y_xxyyy_xxxy, \
                             tr_y_xxyyy_xxxz, \
                             tr_y_xxyyy_xxyy, \
                             tr_y_xxyyy_xxyz, \
                             tr_y_xxyyy_xxzz, \
                             tr_y_xxyyy_xyyy, \
                             tr_y_xxyyy_xyyz, \
                             tr_y_xxyyy_xyzz, \
                             tr_y_xxyyy_xzzz, \
                             tr_y_xxyyy_yyyy, \
                             tr_y_xxyyy_yyyz, \
                             tr_y_xxyyy_yyzz, \
                             tr_y_xxyyy_yzzz, \
                             tr_y_xxyyy_zzzz, \
                             tr_y_xyyy_xxx,   \
                             tr_y_xyyy_xxxx,  \
                             tr_y_xyyy_xxxy,  \
                             tr_y_xyyy_xxxz,  \
                             tr_y_xyyy_xxy,   \
                             tr_y_xyyy_xxyy,  \
                             tr_y_xyyy_xxyz,  \
                             tr_y_xyyy_xxz,   \
                             tr_y_xyyy_xxzz,  \
                             tr_y_xyyy_xyy,   \
                             tr_y_xyyy_xyyy,  \
                             tr_y_xyyy_xyyz,  \
                             tr_y_xyyy_xyz,   \
                             tr_y_xyyy_xyzz,  \
                             tr_y_xyyy_xzz,   \
                             tr_y_xyyy_xzzz,  \
                             tr_y_xyyy_yyy,   \
                             tr_y_xyyy_yyyy,  \
                             tr_y_xyyy_yyyz,  \
                             tr_y_xyyy_yyz,   \
                             tr_y_xyyy_yyzz,  \
                             tr_y_xyyy_yzz,   \
                             tr_y_xyyy_yzzz,  \
                             tr_y_xyyy_zzz,   \
                             tr_y_xyyy_zzzz,  \
                             tr_y_yyy_xxxx,   \
                             tr_y_yyy_xxxy,   \
                             tr_y_yyy_xxxz,   \
                             tr_y_yyy_xxyy,   \
                             tr_y_yyy_xxyz,   \
                             tr_y_yyy_xxzz,   \
                             tr_y_yyy_xyyy,   \
                             tr_y_yyy_xyyz,   \
                             tr_y_yyy_xyzz,   \
                             tr_y_yyy_xzzz,   \
                             tr_y_yyy_yyyy,   \
                             tr_y_yyy_yyyz,   \
                             tr_y_yyy_yyzz,   \
                             tr_y_yyy_yzzz,   \
                             tr_y_yyy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyy_xxxx[i] = tr_y_yyy_xxxx[i] * fe_0 + 4.0 * tr_y_xyyy_xxx[i] * fe_0 + tr_y_xyyy_xxxx[i] * pa_x[i];

        tr_y_xxyyy_xxxy[i] = tr_y_yyy_xxxy[i] * fe_0 + 3.0 * tr_y_xyyy_xxy[i] * fe_0 + tr_y_xyyy_xxxy[i] * pa_x[i];

        tr_y_xxyyy_xxxz[i] = tr_y_yyy_xxxz[i] * fe_0 + 3.0 * tr_y_xyyy_xxz[i] * fe_0 + tr_y_xyyy_xxxz[i] * pa_x[i];

        tr_y_xxyyy_xxyy[i] = tr_y_yyy_xxyy[i] * fe_0 + 2.0 * tr_y_xyyy_xyy[i] * fe_0 + tr_y_xyyy_xxyy[i] * pa_x[i];

        tr_y_xxyyy_xxyz[i] = tr_y_yyy_xxyz[i] * fe_0 + 2.0 * tr_y_xyyy_xyz[i] * fe_0 + tr_y_xyyy_xxyz[i] * pa_x[i];

        tr_y_xxyyy_xxzz[i] = tr_y_yyy_xxzz[i] * fe_0 + 2.0 * tr_y_xyyy_xzz[i] * fe_0 + tr_y_xyyy_xxzz[i] * pa_x[i];

        tr_y_xxyyy_xyyy[i] = tr_y_yyy_xyyy[i] * fe_0 + tr_y_xyyy_yyy[i] * fe_0 + tr_y_xyyy_xyyy[i] * pa_x[i];

        tr_y_xxyyy_xyyz[i] = tr_y_yyy_xyyz[i] * fe_0 + tr_y_xyyy_yyz[i] * fe_0 + tr_y_xyyy_xyyz[i] * pa_x[i];

        tr_y_xxyyy_xyzz[i] = tr_y_yyy_xyzz[i] * fe_0 + tr_y_xyyy_yzz[i] * fe_0 + tr_y_xyyy_xyzz[i] * pa_x[i];

        tr_y_xxyyy_xzzz[i] = tr_y_yyy_xzzz[i] * fe_0 + tr_y_xyyy_zzz[i] * fe_0 + tr_y_xyyy_xzzz[i] * pa_x[i];

        tr_y_xxyyy_yyyy[i] = tr_y_yyy_yyyy[i] * fe_0 + tr_y_xyyy_yyyy[i] * pa_x[i];

        tr_y_xxyyy_yyyz[i] = tr_y_yyy_yyyz[i] * fe_0 + tr_y_xyyy_yyyz[i] * pa_x[i];

        tr_y_xxyyy_yyzz[i] = tr_y_yyy_yyzz[i] * fe_0 + tr_y_xyyy_yyzz[i] * pa_x[i];

        tr_y_xxyyy_yzzz[i] = tr_y_yyy_yzzz[i] * fe_0 + tr_y_xyyy_yzzz[i] * pa_x[i];

        tr_y_xxyyy_zzzz[i] = tr_y_yyy_zzzz[i] * fe_0 + tr_y_xyyy_zzzz[i] * pa_x[i];
    }

    // Set up 420-435 components of targeted buffer : HG

    auto tr_y_xxyyz_xxxx = pbuffer.data(idx_dip_hg + 420);

    auto tr_y_xxyyz_xxxy = pbuffer.data(idx_dip_hg + 421);

    auto tr_y_xxyyz_xxxz = pbuffer.data(idx_dip_hg + 422);

    auto tr_y_xxyyz_xxyy = pbuffer.data(idx_dip_hg + 423);

    auto tr_y_xxyyz_xxyz = pbuffer.data(idx_dip_hg + 424);

    auto tr_y_xxyyz_xxzz = pbuffer.data(idx_dip_hg + 425);

    auto tr_y_xxyyz_xyyy = pbuffer.data(idx_dip_hg + 426);

    auto tr_y_xxyyz_xyyz = pbuffer.data(idx_dip_hg + 427);

    auto tr_y_xxyyz_xyzz = pbuffer.data(idx_dip_hg + 428);

    auto tr_y_xxyyz_xzzz = pbuffer.data(idx_dip_hg + 429);

    auto tr_y_xxyyz_yyyy = pbuffer.data(idx_dip_hg + 430);

    auto tr_y_xxyyz_yyyz = pbuffer.data(idx_dip_hg + 431);

    auto tr_y_xxyyz_yyzz = pbuffer.data(idx_dip_hg + 432);

    auto tr_y_xxyyz_yzzz = pbuffer.data(idx_dip_hg + 433);

    auto tr_y_xxyyz_zzzz = pbuffer.data(idx_dip_hg + 434);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             tr_y_xxyy_xxx,   \
                             tr_y_xxyy_xxxx,  \
                             tr_y_xxyy_xxxy,  \
                             tr_y_xxyy_xxxz,  \
                             tr_y_xxyy_xxy,   \
                             tr_y_xxyy_xxyy,  \
                             tr_y_xxyy_xxyz,  \
                             tr_y_xxyy_xxz,   \
                             tr_y_xxyy_xxzz,  \
                             tr_y_xxyy_xyy,   \
                             tr_y_xxyy_xyyy,  \
                             tr_y_xxyy_xyyz,  \
                             tr_y_xxyy_xyz,   \
                             tr_y_xxyy_xyzz,  \
                             tr_y_xxyy_xzz,   \
                             tr_y_xxyy_xzzz,  \
                             tr_y_xxyy_yyyy,  \
                             tr_y_xxyyz_xxxx, \
                             tr_y_xxyyz_xxxy, \
                             tr_y_xxyyz_xxxz, \
                             tr_y_xxyyz_xxyy, \
                             tr_y_xxyyz_xxyz, \
                             tr_y_xxyyz_xxzz, \
                             tr_y_xxyyz_xyyy, \
                             tr_y_xxyyz_xyyz, \
                             tr_y_xxyyz_xyzz, \
                             tr_y_xxyyz_xzzz, \
                             tr_y_xxyyz_yyyy, \
                             tr_y_xxyyz_yyyz, \
                             tr_y_xxyyz_yyzz, \
                             tr_y_xxyyz_yzzz, \
                             tr_y_xxyyz_zzzz, \
                             tr_y_xyyz_yyyz,  \
                             tr_y_xyyz_yyzz,  \
                             tr_y_xyyz_yzzz,  \
                             tr_y_xyyz_zzzz,  \
                             tr_y_yyz_yyyz,   \
                             tr_y_yyz_yyzz,   \
                             tr_y_yyz_yzzz,   \
                             tr_y_yyz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyz_xxxx[i] = tr_y_xxyy_xxxx[i] * pa_z[i];

        tr_y_xxyyz_xxxy[i] = tr_y_xxyy_xxxy[i] * pa_z[i];

        tr_y_xxyyz_xxxz[i] = tr_y_xxyy_xxx[i] * fe_0 + tr_y_xxyy_xxxz[i] * pa_z[i];

        tr_y_xxyyz_xxyy[i] = tr_y_xxyy_xxyy[i] * pa_z[i];

        tr_y_xxyyz_xxyz[i] = tr_y_xxyy_xxy[i] * fe_0 + tr_y_xxyy_xxyz[i] * pa_z[i];

        tr_y_xxyyz_xxzz[i] = 2.0 * tr_y_xxyy_xxz[i] * fe_0 + tr_y_xxyy_xxzz[i] * pa_z[i];

        tr_y_xxyyz_xyyy[i] = tr_y_xxyy_xyyy[i] * pa_z[i];

        tr_y_xxyyz_xyyz[i] = tr_y_xxyy_xyy[i] * fe_0 + tr_y_xxyy_xyyz[i] * pa_z[i];

        tr_y_xxyyz_xyzz[i] = 2.0 * tr_y_xxyy_xyz[i] * fe_0 + tr_y_xxyy_xyzz[i] * pa_z[i];

        tr_y_xxyyz_xzzz[i] = 3.0 * tr_y_xxyy_xzz[i] * fe_0 + tr_y_xxyy_xzzz[i] * pa_z[i];

        tr_y_xxyyz_yyyy[i] = tr_y_xxyy_yyyy[i] * pa_z[i];

        tr_y_xxyyz_yyyz[i] = tr_y_yyz_yyyz[i] * fe_0 + tr_y_xyyz_yyyz[i] * pa_x[i];

        tr_y_xxyyz_yyzz[i] = tr_y_yyz_yyzz[i] * fe_0 + tr_y_xyyz_yyzz[i] * pa_x[i];

        tr_y_xxyyz_yzzz[i] = tr_y_yyz_yzzz[i] * fe_0 + tr_y_xyyz_yzzz[i] * pa_x[i];

        tr_y_xxyyz_zzzz[i] = tr_y_yyz_zzzz[i] * fe_0 + tr_y_xyyz_zzzz[i] * pa_x[i];
    }

    // Set up 435-450 components of targeted buffer : HG

    auto tr_y_xxyzz_xxxx = pbuffer.data(idx_dip_hg + 435);

    auto tr_y_xxyzz_xxxy = pbuffer.data(idx_dip_hg + 436);

    auto tr_y_xxyzz_xxxz = pbuffer.data(idx_dip_hg + 437);

    auto tr_y_xxyzz_xxyy = pbuffer.data(idx_dip_hg + 438);

    auto tr_y_xxyzz_xxyz = pbuffer.data(idx_dip_hg + 439);

    auto tr_y_xxyzz_xxzz = pbuffer.data(idx_dip_hg + 440);

    auto tr_y_xxyzz_xyyy = pbuffer.data(idx_dip_hg + 441);

    auto tr_y_xxyzz_xyyz = pbuffer.data(idx_dip_hg + 442);

    auto tr_y_xxyzz_xyzz = pbuffer.data(idx_dip_hg + 443);

    auto tr_y_xxyzz_xzzz = pbuffer.data(idx_dip_hg + 444);

    auto tr_y_xxyzz_yyyy = pbuffer.data(idx_dip_hg + 445);

    auto tr_y_xxyzz_yyyz = pbuffer.data(idx_dip_hg + 446);

    auto tr_y_xxyzz_yyzz = pbuffer.data(idx_dip_hg + 447);

    auto tr_y_xxyzz_yzzz = pbuffer.data(idx_dip_hg + 448);

    auto tr_y_xxyzz_zzzz = pbuffer.data(idx_dip_hg + 449);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             tr_y_xxy_xxxy,   \
                             tr_y_xxy_xxyy,   \
                             tr_y_xxy_xyyy,   \
                             tr_y_xxyz_xxxy,  \
                             tr_y_xxyz_xxyy,  \
                             tr_y_xxyz_xyyy,  \
                             tr_y_xxyzz_xxxx, \
                             tr_y_xxyzz_xxxy, \
                             tr_y_xxyzz_xxxz, \
                             tr_y_xxyzz_xxyy, \
                             tr_y_xxyzz_xxyz, \
                             tr_y_xxyzz_xxzz, \
                             tr_y_xxyzz_xyyy, \
                             tr_y_xxyzz_xyyz, \
                             tr_y_xxyzz_xyzz, \
                             tr_y_xxyzz_xzzz, \
                             tr_y_xxyzz_yyyy, \
                             tr_y_xxyzz_yyyz, \
                             tr_y_xxyzz_yyzz, \
                             tr_y_xxyzz_yzzz, \
                             tr_y_xxyzz_zzzz, \
                             tr_y_xxzz_xxxx,  \
                             tr_y_xxzz_xxxz,  \
                             tr_y_xxzz_xxzz,  \
                             tr_y_xxzz_xzzz,  \
                             tr_y_xyzz_xxyz,  \
                             tr_y_xyzz_xyyz,  \
                             tr_y_xyzz_xyz,   \
                             tr_y_xyzz_xyzz,  \
                             tr_y_xyzz_yyyy,  \
                             tr_y_xyzz_yyyz,  \
                             tr_y_xyzz_yyz,   \
                             tr_y_xyzz_yyzz,  \
                             tr_y_xyzz_yzz,   \
                             tr_y_xyzz_yzzz,  \
                             tr_y_xyzz_zzzz,  \
                             tr_y_yzz_xxyz,   \
                             tr_y_yzz_xyyz,   \
                             tr_y_yzz_xyzz,   \
                             tr_y_yzz_yyyy,   \
                             tr_y_yzz_yyyz,   \
                             tr_y_yzz_yyzz,   \
                             tr_y_yzz_yzzz,   \
                             tr_y_yzz_zzzz,   \
                             ts_xxzz_xxxx,    \
                             ts_xxzz_xxxz,    \
                             ts_xxzz_xxzz,    \
                             ts_xxzz_xzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyzz_xxxx[i] = ts_xxzz_xxxx[i] * fe_0 + tr_y_xxzz_xxxx[i] * pa_y[i];

        tr_y_xxyzz_xxxy[i] = tr_y_xxy_xxxy[i] * fe_0 + tr_y_xxyz_xxxy[i] * pa_z[i];

        tr_y_xxyzz_xxxz[i] = ts_xxzz_xxxz[i] * fe_0 + tr_y_xxzz_xxxz[i] * pa_y[i];

        tr_y_xxyzz_xxyy[i] = tr_y_xxy_xxyy[i] * fe_0 + tr_y_xxyz_xxyy[i] * pa_z[i];

        tr_y_xxyzz_xxyz[i] = tr_y_yzz_xxyz[i] * fe_0 + 2.0 * tr_y_xyzz_xyz[i] * fe_0 + tr_y_xyzz_xxyz[i] * pa_x[i];

        tr_y_xxyzz_xxzz[i] = ts_xxzz_xxzz[i] * fe_0 + tr_y_xxzz_xxzz[i] * pa_y[i];

        tr_y_xxyzz_xyyy[i] = tr_y_xxy_xyyy[i] * fe_0 + tr_y_xxyz_xyyy[i] * pa_z[i];

        tr_y_xxyzz_xyyz[i] = tr_y_yzz_xyyz[i] * fe_0 + tr_y_xyzz_yyz[i] * fe_0 + tr_y_xyzz_xyyz[i] * pa_x[i];

        tr_y_xxyzz_xyzz[i] = tr_y_yzz_xyzz[i] * fe_0 + tr_y_xyzz_yzz[i] * fe_0 + tr_y_xyzz_xyzz[i] * pa_x[i];

        tr_y_xxyzz_xzzz[i] = ts_xxzz_xzzz[i] * fe_0 + tr_y_xxzz_xzzz[i] * pa_y[i];

        tr_y_xxyzz_yyyy[i] = tr_y_yzz_yyyy[i] * fe_0 + tr_y_xyzz_yyyy[i] * pa_x[i];

        tr_y_xxyzz_yyyz[i] = tr_y_yzz_yyyz[i] * fe_0 + tr_y_xyzz_yyyz[i] * pa_x[i];

        tr_y_xxyzz_yyzz[i] = tr_y_yzz_yyzz[i] * fe_0 + tr_y_xyzz_yyzz[i] * pa_x[i];

        tr_y_xxyzz_yzzz[i] = tr_y_yzz_yzzz[i] * fe_0 + tr_y_xyzz_yzzz[i] * pa_x[i];

        tr_y_xxyzz_zzzz[i] = tr_y_yzz_zzzz[i] * fe_0 + tr_y_xyzz_zzzz[i] * pa_x[i];
    }

    // Set up 450-465 components of targeted buffer : HG

    auto tr_y_xxzzz_xxxx = pbuffer.data(idx_dip_hg + 450);

    auto tr_y_xxzzz_xxxy = pbuffer.data(idx_dip_hg + 451);

    auto tr_y_xxzzz_xxxz = pbuffer.data(idx_dip_hg + 452);

    auto tr_y_xxzzz_xxyy = pbuffer.data(idx_dip_hg + 453);

    auto tr_y_xxzzz_xxyz = pbuffer.data(idx_dip_hg + 454);

    auto tr_y_xxzzz_xxzz = pbuffer.data(idx_dip_hg + 455);

    auto tr_y_xxzzz_xyyy = pbuffer.data(idx_dip_hg + 456);

    auto tr_y_xxzzz_xyyz = pbuffer.data(idx_dip_hg + 457);

    auto tr_y_xxzzz_xyzz = pbuffer.data(idx_dip_hg + 458);

    auto tr_y_xxzzz_xzzz = pbuffer.data(idx_dip_hg + 459);

    auto tr_y_xxzzz_yyyy = pbuffer.data(idx_dip_hg + 460);

    auto tr_y_xxzzz_yyyz = pbuffer.data(idx_dip_hg + 461);

    auto tr_y_xxzzz_yyzz = pbuffer.data(idx_dip_hg + 462);

    auto tr_y_xxzzz_yzzz = pbuffer.data(idx_dip_hg + 463);

    auto tr_y_xxzzz_zzzz = pbuffer.data(idx_dip_hg + 464);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             tr_y_xxz_xxxx,   \
                             tr_y_xxz_xxxy,   \
                             tr_y_xxz_xxyy,   \
                             tr_y_xxz_xyyy,   \
                             tr_y_xxzz_xxxx,  \
                             tr_y_xxzz_xxxy,  \
                             tr_y_xxzz_xxyy,  \
                             tr_y_xxzz_xyyy,  \
                             tr_y_xxzzz_xxxx, \
                             tr_y_xxzzz_xxxy, \
                             tr_y_xxzzz_xxxz, \
                             tr_y_xxzzz_xxyy, \
                             tr_y_xxzzz_xxyz, \
                             tr_y_xxzzz_xxzz, \
                             tr_y_xxzzz_xyyy, \
                             tr_y_xxzzz_xyyz, \
                             tr_y_xxzzz_xyzz, \
                             tr_y_xxzzz_xzzz, \
                             tr_y_xxzzz_yyyy, \
                             tr_y_xxzzz_yyyz, \
                             tr_y_xxzzz_yyzz, \
                             tr_y_xxzzz_yzzz, \
                             tr_y_xxzzz_zzzz, \
                             tr_y_xzzz_xxxz,  \
                             tr_y_xzzz_xxyz,  \
                             tr_y_xzzz_xxz,   \
                             tr_y_xzzz_xxzz,  \
                             tr_y_xzzz_xyyz,  \
                             tr_y_xzzz_xyz,   \
                             tr_y_xzzz_xyzz,  \
                             tr_y_xzzz_xzz,   \
                             tr_y_xzzz_xzzz,  \
                             tr_y_xzzz_yyyy,  \
                             tr_y_xzzz_yyyz,  \
                             tr_y_xzzz_yyz,   \
                             tr_y_xzzz_yyzz,  \
                             tr_y_xzzz_yzz,   \
                             tr_y_xzzz_yzzz,  \
                             tr_y_xzzz_zzz,   \
                             tr_y_xzzz_zzzz,  \
                             tr_y_zzz_xxxz,   \
                             tr_y_zzz_xxyz,   \
                             tr_y_zzz_xxzz,   \
                             tr_y_zzz_xyyz,   \
                             tr_y_zzz_xyzz,   \
                             tr_y_zzz_xzzz,   \
                             tr_y_zzz_yyyy,   \
                             tr_y_zzz_yyyz,   \
                             tr_y_zzz_yyzz,   \
                             tr_y_zzz_yzzz,   \
                             tr_y_zzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzzz_xxxx[i] = 2.0 * tr_y_xxz_xxxx[i] * fe_0 + tr_y_xxzz_xxxx[i] * pa_z[i];

        tr_y_xxzzz_xxxy[i] = 2.0 * tr_y_xxz_xxxy[i] * fe_0 + tr_y_xxzz_xxxy[i] * pa_z[i];

        tr_y_xxzzz_xxxz[i] = tr_y_zzz_xxxz[i] * fe_0 + 3.0 * tr_y_xzzz_xxz[i] * fe_0 + tr_y_xzzz_xxxz[i] * pa_x[i];

        tr_y_xxzzz_xxyy[i] = 2.0 * tr_y_xxz_xxyy[i] * fe_0 + tr_y_xxzz_xxyy[i] * pa_z[i];

        tr_y_xxzzz_xxyz[i] = tr_y_zzz_xxyz[i] * fe_0 + 2.0 * tr_y_xzzz_xyz[i] * fe_0 + tr_y_xzzz_xxyz[i] * pa_x[i];

        tr_y_xxzzz_xxzz[i] = tr_y_zzz_xxzz[i] * fe_0 + 2.0 * tr_y_xzzz_xzz[i] * fe_0 + tr_y_xzzz_xxzz[i] * pa_x[i];

        tr_y_xxzzz_xyyy[i] = 2.0 * tr_y_xxz_xyyy[i] * fe_0 + tr_y_xxzz_xyyy[i] * pa_z[i];

        tr_y_xxzzz_xyyz[i] = tr_y_zzz_xyyz[i] * fe_0 + tr_y_xzzz_yyz[i] * fe_0 + tr_y_xzzz_xyyz[i] * pa_x[i];

        tr_y_xxzzz_xyzz[i] = tr_y_zzz_xyzz[i] * fe_0 + tr_y_xzzz_yzz[i] * fe_0 + tr_y_xzzz_xyzz[i] * pa_x[i];

        tr_y_xxzzz_xzzz[i] = tr_y_zzz_xzzz[i] * fe_0 + tr_y_xzzz_zzz[i] * fe_0 + tr_y_xzzz_xzzz[i] * pa_x[i];

        tr_y_xxzzz_yyyy[i] = tr_y_zzz_yyyy[i] * fe_0 + tr_y_xzzz_yyyy[i] * pa_x[i];

        tr_y_xxzzz_yyyz[i] = tr_y_zzz_yyyz[i] * fe_0 + tr_y_xzzz_yyyz[i] * pa_x[i];

        tr_y_xxzzz_yyzz[i] = tr_y_zzz_yyzz[i] * fe_0 + tr_y_xzzz_yyzz[i] * pa_x[i];

        tr_y_xxzzz_yzzz[i] = tr_y_zzz_yzzz[i] * fe_0 + tr_y_xzzz_yzzz[i] * pa_x[i];

        tr_y_xxzzz_zzzz[i] = tr_y_zzz_zzzz[i] * fe_0 + tr_y_xzzz_zzzz[i] * pa_x[i];
    }

    // Set up 465-480 components of targeted buffer : HG

    auto tr_y_xyyyy_xxxx = pbuffer.data(idx_dip_hg + 465);

    auto tr_y_xyyyy_xxxy = pbuffer.data(idx_dip_hg + 466);

    auto tr_y_xyyyy_xxxz = pbuffer.data(idx_dip_hg + 467);

    auto tr_y_xyyyy_xxyy = pbuffer.data(idx_dip_hg + 468);

    auto tr_y_xyyyy_xxyz = pbuffer.data(idx_dip_hg + 469);

    auto tr_y_xyyyy_xxzz = pbuffer.data(idx_dip_hg + 470);

    auto tr_y_xyyyy_xyyy = pbuffer.data(idx_dip_hg + 471);

    auto tr_y_xyyyy_xyyz = pbuffer.data(idx_dip_hg + 472);

    auto tr_y_xyyyy_xyzz = pbuffer.data(idx_dip_hg + 473);

    auto tr_y_xyyyy_xzzz = pbuffer.data(idx_dip_hg + 474);

    auto tr_y_xyyyy_yyyy = pbuffer.data(idx_dip_hg + 475);

    auto tr_y_xyyyy_yyyz = pbuffer.data(idx_dip_hg + 476);

    auto tr_y_xyyyy_yyzz = pbuffer.data(idx_dip_hg + 477);

    auto tr_y_xyyyy_yzzz = pbuffer.data(idx_dip_hg + 478);

    auto tr_y_xyyyy_zzzz = pbuffer.data(idx_dip_hg + 479);

#pragma omp simd aligned(pa_x,                \
                             tr_y_xyyyy_xxxx, \
                             tr_y_xyyyy_xxxy, \
                             tr_y_xyyyy_xxxz, \
                             tr_y_xyyyy_xxyy, \
                             tr_y_xyyyy_xxyz, \
                             tr_y_xyyyy_xxzz, \
                             tr_y_xyyyy_xyyy, \
                             tr_y_xyyyy_xyyz, \
                             tr_y_xyyyy_xyzz, \
                             tr_y_xyyyy_xzzz, \
                             tr_y_xyyyy_yyyy, \
                             tr_y_xyyyy_yyyz, \
                             tr_y_xyyyy_yyzz, \
                             tr_y_xyyyy_yzzz, \
                             tr_y_xyyyy_zzzz, \
                             tr_y_yyyy_xxx,   \
                             tr_y_yyyy_xxxx,  \
                             tr_y_yyyy_xxxy,  \
                             tr_y_yyyy_xxxz,  \
                             tr_y_yyyy_xxy,   \
                             tr_y_yyyy_xxyy,  \
                             tr_y_yyyy_xxyz,  \
                             tr_y_yyyy_xxz,   \
                             tr_y_yyyy_xxzz,  \
                             tr_y_yyyy_xyy,   \
                             tr_y_yyyy_xyyy,  \
                             tr_y_yyyy_xyyz,  \
                             tr_y_yyyy_xyz,   \
                             tr_y_yyyy_xyzz,  \
                             tr_y_yyyy_xzz,   \
                             tr_y_yyyy_xzzz,  \
                             tr_y_yyyy_yyy,   \
                             tr_y_yyyy_yyyy,  \
                             tr_y_yyyy_yyyz,  \
                             tr_y_yyyy_yyz,   \
                             tr_y_yyyy_yyzz,  \
                             tr_y_yyyy_yzz,   \
                             tr_y_yyyy_yzzz,  \
                             tr_y_yyyy_zzz,   \
                             tr_y_yyyy_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyy_xxxx[i] = 4.0 * tr_y_yyyy_xxx[i] * fe_0 + tr_y_yyyy_xxxx[i] * pa_x[i];

        tr_y_xyyyy_xxxy[i] = 3.0 * tr_y_yyyy_xxy[i] * fe_0 + tr_y_yyyy_xxxy[i] * pa_x[i];

        tr_y_xyyyy_xxxz[i] = 3.0 * tr_y_yyyy_xxz[i] * fe_0 + tr_y_yyyy_xxxz[i] * pa_x[i];

        tr_y_xyyyy_xxyy[i] = 2.0 * tr_y_yyyy_xyy[i] * fe_0 + tr_y_yyyy_xxyy[i] * pa_x[i];

        tr_y_xyyyy_xxyz[i] = 2.0 * tr_y_yyyy_xyz[i] * fe_0 + tr_y_yyyy_xxyz[i] * pa_x[i];

        tr_y_xyyyy_xxzz[i] = 2.0 * tr_y_yyyy_xzz[i] * fe_0 + tr_y_yyyy_xxzz[i] * pa_x[i];

        tr_y_xyyyy_xyyy[i] = tr_y_yyyy_yyy[i] * fe_0 + tr_y_yyyy_xyyy[i] * pa_x[i];

        tr_y_xyyyy_xyyz[i] = tr_y_yyyy_yyz[i] * fe_0 + tr_y_yyyy_xyyz[i] * pa_x[i];

        tr_y_xyyyy_xyzz[i] = tr_y_yyyy_yzz[i] * fe_0 + tr_y_yyyy_xyzz[i] * pa_x[i];

        tr_y_xyyyy_xzzz[i] = tr_y_yyyy_zzz[i] * fe_0 + tr_y_yyyy_xzzz[i] * pa_x[i];

        tr_y_xyyyy_yyyy[i] = tr_y_yyyy_yyyy[i] * pa_x[i];

        tr_y_xyyyy_yyyz[i] = tr_y_yyyy_yyyz[i] * pa_x[i];

        tr_y_xyyyy_yyzz[i] = tr_y_yyyy_yyzz[i] * pa_x[i];

        tr_y_xyyyy_yzzz[i] = tr_y_yyyy_yzzz[i] * pa_x[i];

        tr_y_xyyyy_zzzz[i] = tr_y_yyyy_zzzz[i] * pa_x[i];
    }

    // Set up 480-495 components of targeted buffer : HG

    auto tr_y_xyyyz_xxxx = pbuffer.data(idx_dip_hg + 480);

    auto tr_y_xyyyz_xxxy = pbuffer.data(idx_dip_hg + 481);

    auto tr_y_xyyyz_xxxz = pbuffer.data(idx_dip_hg + 482);

    auto tr_y_xyyyz_xxyy = pbuffer.data(idx_dip_hg + 483);

    auto tr_y_xyyyz_xxyz = pbuffer.data(idx_dip_hg + 484);

    auto tr_y_xyyyz_xxzz = pbuffer.data(idx_dip_hg + 485);

    auto tr_y_xyyyz_xyyy = pbuffer.data(idx_dip_hg + 486);

    auto tr_y_xyyyz_xyyz = pbuffer.data(idx_dip_hg + 487);

    auto tr_y_xyyyz_xyzz = pbuffer.data(idx_dip_hg + 488);

    auto tr_y_xyyyz_xzzz = pbuffer.data(idx_dip_hg + 489);

    auto tr_y_xyyyz_yyyy = pbuffer.data(idx_dip_hg + 490);

    auto tr_y_xyyyz_yyyz = pbuffer.data(idx_dip_hg + 491);

    auto tr_y_xyyyz_yyzz = pbuffer.data(idx_dip_hg + 492);

    auto tr_y_xyyyz_yzzz = pbuffer.data(idx_dip_hg + 493);

    auto tr_y_xyyyz_zzzz = pbuffer.data(idx_dip_hg + 494);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             tr_y_xyyy_xxxx,  \
                             tr_y_xyyy_xxxy,  \
                             tr_y_xyyy_xxyy,  \
                             tr_y_xyyy_xyyy,  \
                             tr_y_xyyyz_xxxx, \
                             tr_y_xyyyz_xxxy, \
                             tr_y_xyyyz_xxxz, \
                             tr_y_xyyyz_xxyy, \
                             tr_y_xyyyz_xxyz, \
                             tr_y_xyyyz_xxzz, \
                             tr_y_xyyyz_xyyy, \
                             tr_y_xyyyz_xyyz, \
                             tr_y_xyyyz_xyzz, \
                             tr_y_xyyyz_xzzz, \
                             tr_y_xyyyz_yyyy, \
                             tr_y_xyyyz_yyyz, \
                             tr_y_xyyyz_yyzz, \
                             tr_y_xyyyz_yzzz, \
                             tr_y_xyyyz_zzzz, \
                             tr_y_yyyz_xxxz,  \
                             tr_y_yyyz_xxyz,  \
                             tr_y_yyyz_xxz,   \
                             tr_y_yyyz_xxzz,  \
                             tr_y_yyyz_xyyz,  \
                             tr_y_yyyz_xyz,   \
                             tr_y_yyyz_xyzz,  \
                             tr_y_yyyz_xzz,   \
                             tr_y_yyyz_xzzz,  \
                             tr_y_yyyz_yyyy,  \
                             tr_y_yyyz_yyyz,  \
                             tr_y_yyyz_yyz,   \
                             tr_y_yyyz_yyzz,  \
                             tr_y_yyyz_yzz,   \
                             tr_y_yyyz_yzzz,  \
                             tr_y_yyyz_zzz,   \
                             tr_y_yyyz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyz_xxxx[i] = tr_y_xyyy_xxxx[i] * pa_z[i];

        tr_y_xyyyz_xxxy[i] = tr_y_xyyy_xxxy[i] * pa_z[i];

        tr_y_xyyyz_xxxz[i] = 3.0 * tr_y_yyyz_xxz[i] * fe_0 + tr_y_yyyz_xxxz[i] * pa_x[i];

        tr_y_xyyyz_xxyy[i] = tr_y_xyyy_xxyy[i] * pa_z[i];

        tr_y_xyyyz_xxyz[i] = 2.0 * tr_y_yyyz_xyz[i] * fe_0 + tr_y_yyyz_xxyz[i] * pa_x[i];

        tr_y_xyyyz_xxzz[i] = 2.0 * tr_y_yyyz_xzz[i] * fe_0 + tr_y_yyyz_xxzz[i] * pa_x[i];

        tr_y_xyyyz_xyyy[i] = tr_y_xyyy_xyyy[i] * pa_z[i];

        tr_y_xyyyz_xyyz[i] = tr_y_yyyz_yyz[i] * fe_0 + tr_y_yyyz_xyyz[i] * pa_x[i];

        tr_y_xyyyz_xyzz[i] = tr_y_yyyz_yzz[i] * fe_0 + tr_y_yyyz_xyzz[i] * pa_x[i];

        tr_y_xyyyz_xzzz[i] = tr_y_yyyz_zzz[i] * fe_0 + tr_y_yyyz_xzzz[i] * pa_x[i];

        tr_y_xyyyz_yyyy[i] = tr_y_yyyz_yyyy[i] * pa_x[i];

        tr_y_xyyyz_yyyz[i] = tr_y_yyyz_yyyz[i] * pa_x[i];

        tr_y_xyyyz_yyzz[i] = tr_y_yyyz_yyzz[i] * pa_x[i];

        tr_y_xyyyz_yzzz[i] = tr_y_yyyz_yzzz[i] * pa_x[i];

        tr_y_xyyyz_zzzz[i] = tr_y_yyyz_zzzz[i] * pa_x[i];
    }

    // Set up 495-510 components of targeted buffer : HG

    auto tr_y_xyyzz_xxxx = pbuffer.data(idx_dip_hg + 495);

    auto tr_y_xyyzz_xxxy = pbuffer.data(idx_dip_hg + 496);

    auto tr_y_xyyzz_xxxz = pbuffer.data(idx_dip_hg + 497);

    auto tr_y_xyyzz_xxyy = pbuffer.data(idx_dip_hg + 498);

    auto tr_y_xyyzz_xxyz = pbuffer.data(idx_dip_hg + 499);

    auto tr_y_xyyzz_xxzz = pbuffer.data(idx_dip_hg + 500);

    auto tr_y_xyyzz_xyyy = pbuffer.data(idx_dip_hg + 501);

    auto tr_y_xyyzz_xyyz = pbuffer.data(idx_dip_hg + 502);

    auto tr_y_xyyzz_xyzz = pbuffer.data(idx_dip_hg + 503);

    auto tr_y_xyyzz_xzzz = pbuffer.data(idx_dip_hg + 504);

    auto tr_y_xyyzz_yyyy = pbuffer.data(idx_dip_hg + 505);

    auto tr_y_xyyzz_yyyz = pbuffer.data(idx_dip_hg + 506);

    auto tr_y_xyyzz_yyzz = pbuffer.data(idx_dip_hg + 507);

    auto tr_y_xyyzz_yzzz = pbuffer.data(idx_dip_hg + 508);

    auto tr_y_xyyzz_zzzz = pbuffer.data(idx_dip_hg + 509);

#pragma omp simd aligned(pa_x,                \
                             tr_y_xyyzz_xxxx, \
                             tr_y_xyyzz_xxxy, \
                             tr_y_xyyzz_xxxz, \
                             tr_y_xyyzz_xxyy, \
                             tr_y_xyyzz_xxyz, \
                             tr_y_xyyzz_xxzz, \
                             tr_y_xyyzz_xyyy, \
                             tr_y_xyyzz_xyyz, \
                             tr_y_xyyzz_xyzz, \
                             tr_y_xyyzz_xzzz, \
                             tr_y_xyyzz_yyyy, \
                             tr_y_xyyzz_yyyz, \
                             tr_y_xyyzz_yyzz, \
                             tr_y_xyyzz_yzzz, \
                             tr_y_xyyzz_zzzz, \
                             tr_y_yyzz_xxx,   \
                             tr_y_yyzz_xxxx,  \
                             tr_y_yyzz_xxxy,  \
                             tr_y_yyzz_xxxz,  \
                             tr_y_yyzz_xxy,   \
                             tr_y_yyzz_xxyy,  \
                             tr_y_yyzz_xxyz,  \
                             tr_y_yyzz_xxz,   \
                             tr_y_yyzz_xxzz,  \
                             tr_y_yyzz_xyy,   \
                             tr_y_yyzz_xyyy,  \
                             tr_y_yyzz_xyyz,  \
                             tr_y_yyzz_xyz,   \
                             tr_y_yyzz_xyzz,  \
                             tr_y_yyzz_xzz,   \
                             tr_y_yyzz_xzzz,  \
                             tr_y_yyzz_yyy,   \
                             tr_y_yyzz_yyyy,  \
                             tr_y_yyzz_yyyz,  \
                             tr_y_yyzz_yyz,   \
                             tr_y_yyzz_yyzz,  \
                             tr_y_yyzz_yzz,   \
                             tr_y_yyzz_yzzz,  \
                             tr_y_yyzz_zzz,   \
                             tr_y_yyzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyzz_xxxx[i] = 4.0 * tr_y_yyzz_xxx[i] * fe_0 + tr_y_yyzz_xxxx[i] * pa_x[i];

        tr_y_xyyzz_xxxy[i] = 3.0 * tr_y_yyzz_xxy[i] * fe_0 + tr_y_yyzz_xxxy[i] * pa_x[i];

        tr_y_xyyzz_xxxz[i] = 3.0 * tr_y_yyzz_xxz[i] * fe_0 + tr_y_yyzz_xxxz[i] * pa_x[i];

        tr_y_xyyzz_xxyy[i] = 2.0 * tr_y_yyzz_xyy[i] * fe_0 + tr_y_yyzz_xxyy[i] * pa_x[i];

        tr_y_xyyzz_xxyz[i] = 2.0 * tr_y_yyzz_xyz[i] * fe_0 + tr_y_yyzz_xxyz[i] * pa_x[i];

        tr_y_xyyzz_xxzz[i] = 2.0 * tr_y_yyzz_xzz[i] * fe_0 + tr_y_yyzz_xxzz[i] * pa_x[i];

        tr_y_xyyzz_xyyy[i] = tr_y_yyzz_yyy[i] * fe_0 + tr_y_yyzz_xyyy[i] * pa_x[i];

        tr_y_xyyzz_xyyz[i] = tr_y_yyzz_yyz[i] * fe_0 + tr_y_yyzz_xyyz[i] * pa_x[i];

        tr_y_xyyzz_xyzz[i] = tr_y_yyzz_yzz[i] * fe_0 + tr_y_yyzz_xyzz[i] * pa_x[i];

        tr_y_xyyzz_xzzz[i] = tr_y_yyzz_zzz[i] * fe_0 + tr_y_yyzz_xzzz[i] * pa_x[i];

        tr_y_xyyzz_yyyy[i] = tr_y_yyzz_yyyy[i] * pa_x[i];

        tr_y_xyyzz_yyyz[i] = tr_y_yyzz_yyyz[i] * pa_x[i];

        tr_y_xyyzz_yyzz[i] = tr_y_yyzz_yyzz[i] * pa_x[i];

        tr_y_xyyzz_yzzz[i] = tr_y_yyzz_yzzz[i] * pa_x[i];

        tr_y_xyyzz_zzzz[i] = tr_y_yyzz_zzzz[i] * pa_x[i];
    }

    // Set up 510-525 components of targeted buffer : HG

    auto tr_y_xyzzz_xxxx = pbuffer.data(idx_dip_hg + 510);

    auto tr_y_xyzzz_xxxy = pbuffer.data(idx_dip_hg + 511);

    auto tr_y_xyzzz_xxxz = pbuffer.data(idx_dip_hg + 512);

    auto tr_y_xyzzz_xxyy = pbuffer.data(idx_dip_hg + 513);

    auto tr_y_xyzzz_xxyz = pbuffer.data(idx_dip_hg + 514);

    auto tr_y_xyzzz_xxzz = pbuffer.data(idx_dip_hg + 515);

    auto tr_y_xyzzz_xyyy = pbuffer.data(idx_dip_hg + 516);

    auto tr_y_xyzzz_xyyz = pbuffer.data(idx_dip_hg + 517);

    auto tr_y_xyzzz_xyzz = pbuffer.data(idx_dip_hg + 518);

    auto tr_y_xyzzz_xzzz = pbuffer.data(idx_dip_hg + 519);

    auto tr_y_xyzzz_yyyy = pbuffer.data(idx_dip_hg + 520);

    auto tr_y_xyzzz_yyyz = pbuffer.data(idx_dip_hg + 521);

    auto tr_y_xyzzz_yyzz = pbuffer.data(idx_dip_hg + 522);

    auto tr_y_xyzzz_yzzz = pbuffer.data(idx_dip_hg + 523);

    auto tr_y_xyzzz_zzzz = pbuffer.data(idx_dip_hg + 524);

#pragma omp simd aligned(pa_x,                \
                             tr_y_xyzzz_xxxx, \
                             tr_y_xyzzz_xxxy, \
                             tr_y_xyzzz_xxxz, \
                             tr_y_xyzzz_xxyy, \
                             tr_y_xyzzz_xxyz, \
                             tr_y_xyzzz_xxzz, \
                             tr_y_xyzzz_xyyy, \
                             tr_y_xyzzz_xyyz, \
                             tr_y_xyzzz_xyzz, \
                             tr_y_xyzzz_xzzz, \
                             tr_y_xyzzz_yyyy, \
                             tr_y_xyzzz_yyyz, \
                             tr_y_xyzzz_yyzz, \
                             tr_y_xyzzz_yzzz, \
                             tr_y_xyzzz_zzzz, \
                             tr_y_yzzz_xxx,   \
                             tr_y_yzzz_xxxx,  \
                             tr_y_yzzz_xxxy,  \
                             tr_y_yzzz_xxxz,  \
                             tr_y_yzzz_xxy,   \
                             tr_y_yzzz_xxyy,  \
                             tr_y_yzzz_xxyz,  \
                             tr_y_yzzz_xxz,   \
                             tr_y_yzzz_xxzz,  \
                             tr_y_yzzz_xyy,   \
                             tr_y_yzzz_xyyy,  \
                             tr_y_yzzz_xyyz,  \
                             tr_y_yzzz_xyz,   \
                             tr_y_yzzz_xyzz,  \
                             tr_y_yzzz_xzz,   \
                             tr_y_yzzz_xzzz,  \
                             tr_y_yzzz_yyy,   \
                             tr_y_yzzz_yyyy,  \
                             tr_y_yzzz_yyyz,  \
                             tr_y_yzzz_yyz,   \
                             tr_y_yzzz_yyzz,  \
                             tr_y_yzzz_yzz,   \
                             tr_y_yzzz_yzzz,  \
                             tr_y_yzzz_zzz,   \
                             tr_y_yzzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzzz_xxxx[i] = 4.0 * tr_y_yzzz_xxx[i] * fe_0 + tr_y_yzzz_xxxx[i] * pa_x[i];

        tr_y_xyzzz_xxxy[i] = 3.0 * tr_y_yzzz_xxy[i] * fe_0 + tr_y_yzzz_xxxy[i] * pa_x[i];

        tr_y_xyzzz_xxxz[i] = 3.0 * tr_y_yzzz_xxz[i] * fe_0 + tr_y_yzzz_xxxz[i] * pa_x[i];

        tr_y_xyzzz_xxyy[i] = 2.0 * tr_y_yzzz_xyy[i] * fe_0 + tr_y_yzzz_xxyy[i] * pa_x[i];

        tr_y_xyzzz_xxyz[i] = 2.0 * tr_y_yzzz_xyz[i] * fe_0 + tr_y_yzzz_xxyz[i] * pa_x[i];

        tr_y_xyzzz_xxzz[i] = 2.0 * tr_y_yzzz_xzz[i] * fe_0 + tr_y_yzzz_xxzz[i] * pa_x[i];

        tr_y_xyzzz_xyyy[i] = tr_y_yzzz_yyy[i] * fe_0 + tr_y_yzzz_xyyy[i] * pa_x[i];

        tr_y_xyzzz_xyyz[i] = tr_y_yzzz_yyz[i] * fe_0 + tr_y_yzzz_xyyz[i] * pa_x[i];

        tr_y_xyzzz_xyzz[i] = tr_y_yzzz_yzz[i] * fe_0 + tr_y_yzzz_xyzz[i] * pa_x[i];

        tr_y_xyzzz_xzzz[i] = tr_y_yzzz_zzz[i] * fe_0 + tr_y_yzzz_xzzz[i] * pa_x[i];

        tr_y_xyzzz_yyyy[i] = tr_y_yzzz_yyyy[i] * pa_x[i];

        tr_y_xyzzz_yyyz[i] = tr_y_yzzz_yyyz[i] * pa_x[i];

        tr_y_xyzzz_yyzz[i] = tr_y_yzzz_yyzz[i] * pa_x[i];

        tr_y_xyzzz_yzzz[i] = tr_y_yzzz_yzzz[i] * pa_x[i];

        tr_y_xyzzz_zzzz[i] = tr_y_yzzz_zzzz[i] * pa_x[i];
    }

    // Set up 525-540 components of targeted buffer : HG

    auto tr_y_xzzzz_xxxx = pbuffer.data(idx_dip_hg + 525);

    auto tr_y_xzzzz_xxxy = pbuffer.data(idx_dip_hg + 526);

    auto tr_y_xzzzz_xxxz = pbuffer.data(idx_dip_hg + 527);

    auto tr_y_xzzzz_xxyy = pbuffer.data(idx_dip_hg + 528);

    auto tr_y_xzzzz_xxyz = pbuffer.data(idx_dip_hg + 529);

    auto tr_y_xzzzz_xxzz = pbuffer.data(idx_dip_hg + 530);

    auto tr_y_xzzzz_xyyy = pbuffer.data(idx_dip_hg + 531);

    auto tr_y_xzzzz_xyyz = pbuffer.data(idx_dip_hg + 532);

    auto tr_y_xzzzz_xyzz = pbuffer.data(idx_dip_hg + 533);

    auto tr_y_xzzzz_xzzz = pbuffer.data(idx_dip_hg + 534);

    auto tr_y_xzzzz_yyyy = pbuffer.data(idx_dip_hg + 535);

    auto tr_y_xzzzz_yyyz = pbuffer.data(idx_dip_hg + 536);

    auto tr_y_xzzzz_yyzz = pbuffer.data(idx_dip_hg + 537);

    auto tr_y_xzzzz_yzzz = pbuffer.data(idx_dip_hg + 538);

    auto tr_y_xzzzz_zzzz = pbuffer.data(idx_dip_hg + 539);

#pragma omp simd aligned(pa_x,                \
                             tr_y_xzzzz_xxxx, \
                             tr_y_xzzzz_xxxy, \
                             tr_y_xzzzz_xxxz, \
                             tr_y_xzzzz_xxyy, \
                             tr_y_xzzzz_xxyz, \
                             tr_y_xzzzz_xxzz, \
                             tr_y_xzzzz_xyyy, \
                             tr_y_xzzzz_xyyz, \
                             tr_y_xzzzz_xyzz, \
                             tr_y_xzzzz_xzzz, \
                             tr_y_xzzzz_yyyy, \
                             tr_y_xzzzz_yyyz, \
                             tr_y_xzzzz_yyzz, \
                             tr_y_xzzzz_yzzz, \
                             tr_y_xzzzz_zzzz, \
                             tr_y_zzzz_xxx,   \
                             tr_y_zzzz_xxxx,  \
                             tr_y_zzzz_xxxy,  \
                             tr_y_zzzz_xxxz,  \
                             tr_y_zzzz_xxy,   \
                             tr_y_zzzz_xxyy,  \
                             tr_y_zzzz_xxyz,  \
                             tr_y_zzzz_xxz,   \
                             tr_y_zzzz_xxzz,  \
                             tr_y_zzzz_xyy,   \
                             tr_y_zzzz_xyyy,  \
                             tr_y_zzzz_xyyz,  \
                             tr_y_zzzz_xyz,   \
                             tr_y_zzzz_xyzz,  \
                             tr_y_zzzz_xzz,   \
                             tr_y_zzzz_xzzz,  \
                             tr_y_zzzz_yyy,   \
                             tr_y_zzzz_yyyy,  \
                             tr_y_zzzz_yyyz,  \
                             tr_y_zzzz_yyz,   \
                             tr_y_zzzz_yyzz,  \
                             tr_y_zzzz_yzz,   \
                             tr_y_zzzz_yzzz,  \
                             tr_y_zzzz_zzz,   \
                             tr_y_zzzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzzz_xxxx[i] = 4.0 * tr_y_zzzz_xxx[i] * fe_0 + tr_y_zzzz_xxxx[i] * pa_x[i];

        tr_y_xzzzz_xxxy[i] = 3.0 * tr_y_zzzz_xxy[i] * fe_0 + tr_y_zzzz_xxxy[i] * pa_x[i];

        tr_y_xzzzz_xxxz[i] = 3.0 * tr_y_zzzz_xxz[i] * fe_0 + tr_y_zzzz_xxxz[i] * pa_x[i];

        tr_y_xzzzz_xxyy[i] = 2.0 * tr_y_zzzz_xyy[i] * fe_0 + tr_y_zzzz_xxyy[i] * pa_x[i];

        tr_y_xzzzz_xxyz[i] = 2.0 * tr_y_zzzz_xyz[i] * fe_0 + tr_y_zzzz_xxyz[i] * pa_x[i];

        tr_y_xzzzz_xxzz[i] = 2.0 * tr_y_zzzz_xzz[i] * fe_0 + tr_y_zzzz_xxzz[i] * pa_x[i];

        tr_y_xzzzz_xyyy[i] = tr_y_zzzz_yyy[i] * fe_0 + tr_y_zzzz_xyyy[i] * pa_x[i];

        tr_y_xzzzz_xyyz[i] = tr_y_zzzz_yyz[i] * fe_0 + tr_y_zzzz_xyyz[i] * pa_x[i];

        tr_y_xzzzz_xyzz[i] = tr_y_zzzz_yzz[i] * fe_0 + tr_y_zzzz_xyzz[i] * pa_x[i];

        tr_y_xzzzz_xzzz[i] = tr_y_zzzz_zzz[i] * fe_0 + tr_y_zzzz_xzzz[i] * pa_x[i];

        tr_y_xzzzz_yyyy[i] = tr_y_zzzz_yyyy[i] * pa_x[i];

        tr_y_xzzzz_yyyz[i] = tr_y_zzzz_yyyz[i] * pa_x[i];

        tr_y_xzzzz_yyzz[i] = tr_y_zzzz_yyzz[i] * pa_x[i];

        tr_y_xzzzz_yzzz[i] = tr_y_zzzz_yzzz[i] * pa_x[i];

        tr_y_xzzzz_zzzz[i] = tr_y_zzzz_zzzz[i] * pa_x[i];
    }

    // Set up 540-555 components of targeted buffer : HG

    auto tr_y_yyyyy_xxxx = pbuffer.data(idx_dip_hg + 540);

    auto tr_y_yyyyy_xxxy = pbuffer.data(idx_dip_hg + 541);

    auto tr_y_yyyyy_xxxz = pbuffer.data(idx_dip_hg + 542);

    auto tr_y_yyyyy_xxyy = pbuffer.data(idx_dip_hg + 543);

    auto tr_y_yyyyy_xxyz = pbuffer.data(idx_dip_hg + 544);

    auto tr_y_yyyyy_xxzz = pbuffer.data(idx_dip_hg + 545);

    auto tr_y_yyyyy_xyyy = pbuffer.data(idx_dip_hg + 546);

    auto tr_y_yyyyy_xyyz = pbuffer.data(idx_dip_hg + 547);

    auto tr_y_yyyyy_xyzz = pbuffer.data(idx_dip_hg + 548);

    auto tr_y_yyyyy_xzzz = pbuffer.data(idx_dip_hg + 549);

    auto tr_y_yyyyy_yyyy = pbuffer.data(idx_dip_hg + 550);

    auto tr_y_yyyyy_yyyz = pbuffer.data(idx_dip_hg + 551);

    auto tr_y_yyyyy_yyzz = pbuffer.data(idx_dip_hg + 552);

    auto tr_y_yyyyy_yzzz = pbuffer.data(idx_dip_hg + 553);

    auto tr_y_yyyyy_zzzz = pbuffer.data(idx_dip_hg + 554);

#pragma omp simd aligned(pa_y,                \
                             tr_y_yyy_xxxx,   \
                             tr_y_yyy_xxxy,   \
                             tr_y_yyy_xxxz,   \
                             tr_y_yyy_xxyy,   \
                             tr_y_yyy_xxyz,   \
                             tr_y_yyy_xxzz,   \
                             tr_y_yyy_xyyy,   \
                             tr_y_yyy_xyyz,   \
                             tr_y_yyy_xyzz,   \
                             tr_y_yyy_xzzz,   \
                             tr_y_yyy_yyyy,   \
                             tr_y_yyy_yyyz,   \
                             tr_y_yyy_yyzz,   \
                             tr_y_yyy_yzzz,   \
                             tr_y_yyy_zzzz,   \
                             tr_y_yyyy_xxx,   \
                             tr_y_yyyy_xxxx,  \
                             tr_y_yyyy_xxxy,  \
                             tr_y_yyyy_xxxz,  \
                             tr_y_yyyy_xxy,   \
                             tr_y_yyyy_xxyy,  \
                             tr_y_yyyy_xxyz,  \
                             tr_y_yyyy_xxz,   \
                             tr_y_yyyy_xxzz,  \
                             tr_y_yyyy_xyy,   \
                             tr_y_yyyy_xyyy,  \
                             tr_y_yyyy_xyyz,  \
                             tr_y_yyyy_xyz,   \
                             tr_y_yyyy_xyzz,  \
                             tr_y_yyyy_xzz,   \
                             tr_y_yyyy_xzzz,  \
                             tr_y_yyyy_yyy,   \
                             tr_y_yyyy_yyyy,  \
                             tr_y_yyyy_yyyz,  \
                             tr_y_yyyy_yyz,   \
                             tr_y_yyyy_yyzz,  \
                             tr_y_yyyy_yzz,   \
                             tr_y_yyyy_yzzz,  \
                             tr_y_yyyy_zzz,   \
                             tr_y_yyyy_zzzz,  \
                             tr_y_yyyyy_xxxx, \
                             tr_y_yyyyy_xxxy, \
                             tr_y_yyyyy_xxxz, \
                             tr_y_yyyyy_xxyy, \
                             tr_y_yyyyy_xxyz, \
                             tr_y_yyyyy_xxzz, \
                             tr_y_yyyyy_xyyy, \
                             tr_y_yyyyy_xyyz, \
                             tr_y_yyyyy_xyzz, \
                             tr_y_yyyyy_xzzz, \
                             tr_y_yyyyy_yyyy, \
                             tr_y_yyyyy_yyyz, \
                             tr_y_yyyyy_yyzz, \
                             tr_y_yyyyy_yzzz, \
                             tr_y_yyyyy_zzzz, \
                             ts_yyyy_xxxx,    \
                             ts_yyyy_xxxy,    \
                             ts_yyyy_xxxz,    \
                             ts_yyyy_xxyy,    \
                             ts_yyyy_xxyz,    \
                             ts_yyyy_xxzz,    \
                             ts_yyyy_xyyy,    \
                             ts_yyyy_xyyz,    \
                             ts_yyyy_xyzz,    \
                             ts_yyyy_xzzz,    \
                             ts_yyyy_yyyy,    \
                             ts_yyyy_yyyz,    \
                             ts_yyyy_yyzz,    \
                             ts_yyyy_yzzz,    \
                             ts_yyyy_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyy_xxxx[i] = 4.0 * tr_y_yyy_xxxx[i] * fe_0 + ts_yyyy_xxxx[i] * fe_0 + tr_y_yyyy_xxxx[i] * pa_y[i];

        tr_y_yyyyy_xxxy[i] = 4.0 * tr_y_yyy_xxxy[i] * fe_0 + tr_y_yyyy_xxx[i] * fe_0 + ts_yyyy_xxxy[i] * fe_0 + tr_y_yyyy_xxxy[i] * pa_y[i];

        tr_y_yyyyy_xxxz[i] = 4.0 * tr_y_yyy_xxxz[i] * fe_0 + ts_yyyy_xxxz[i] * fe_0 + tr_y_yyyy_xxxz[i] * pa_y[i];

        tr_y_yyyyy_xxyy[i] = 4.0 * tr_y_yyy_xxyy[i] * fe_0 + 2.0 * tr_y_yyyy_xxy[i] * fe_0 + ts_yyyy_xxyy[i] * fe_0 + tr_y_yyyy_xxyy[i] * pa_y[i];

        tr_y_yyyyy_xxyz[i] = 4.0 * tr_y_yyy_xxyz[i] * fe_0 + tr_y_yyyy_xxz[i] * fe_0 + ts_yyyy_xxyz[i] * fe_0 + tr_y_yyyy_xxyz[i] * pa_y[i];

        tr_y_yyyyy_xxzz[i] = 4.0 * tr_y_yyy_xxzz[i] * fe_0 + ts_yyyy_xxzz[i] * fe_0 + tr_y_yyyy_xxzz[i] * pa_y[i];

        tr_y_yyyyy_xyyy[i] = 4.0 * tr_y_yyy_xyyy[i] * fe_0 + 3.0 * tr_y_yyyy_xyy[i] * fe_0 + ts_yyyy_xyyy[i] * fe_0 + tr_y_yyyy_xyyy[i] * pa_y[i];

        tr_y_yyyyy_xyyz[i] = 4.0 * tr_y_yyy_xyyz[i] * fe_0 + 2.0 * tr_y_yyyy_xyz[i] * fe_0 + ts_yyyy_xyyz[i] * fe_0 + tr_y_yyyy_xyyz[i] * pa_y[i];

        tr_y_yyyyy_xyzz[i] = 4.0 * tr_y_yyy_xyzz[i] * fe_0 + tr_y_yyyy_xzz[i] * fe_0 + ts_yyyy_xyzz[i] * fe_0 + tr_y_yyyy_xyzz[i] * pa_y[i];

        tr_y_yyyyy_xzzz[i] = 4.0 * tr_y_yyy_xzzz[i] * fe_0 + ts_yyyy_xzzz[i] * fe_0 + tr_y_yyyy_xzzz[i] * pa_y[i];

        tr_y_yyyyy_yyyy[i] = 4.0 * tr_y_yyy_yyyy[i] * fe_0 + 4.0 * tr_y_yyyy_yyy[i] * fe_0 + ts_yyyy_yyyy[i] * fe_0 + tr_y_yyyy_yyyy[i] * pa_y[i];

        tr_y_yyyyy_yyyz[i] = 4.0 * tr_y_yyy_yyyz[i] * fe_0 + 3.0 * tr_y_yyyy_yyz[i] * fe_0 + ts_yyyy_yyyz[i] * fe_0 + tr_y_yyyy_yyyz[i] * pa_y[i];

        tr_y_yyyyy_yyzz[i] = 4.0 * tr_y_yyy_yyzz[i] * fe_0 + 2.0 * tr_y_yyyy_yzz[i] * fe_0 + ts_yyyy_yyzz[i] * fe_0 + tr_y_yyyy_yyzz[i] * pa_y[i];

        tr_y_yyyyy_yzzz[i] = 4.0 * tr_y_yyy_yzzz[i] * fe_0 + tr_y_yyyy_zzz[i] * fe_0 + ts_yyyy_yzzz[i] * fe_0 + tr_y_yyyy_yzzz[i] * pa_y[i];

        tr_y_yyyyy_zzzz[i] = 4.0 * tr_y_yyy_zzzz[i] * fe_0 + ts_yyyy_zzzz[i] * fe_0 + tr_y_yyyy_zzzz[i] * pa_y[i];
    }

    // Set up 555-570 components of targeted buffer : HG

    auto tr_y_yyyyz_xxxx = pbuffer.data(idx_dip_hg + 555);

    auto tr_y_yyyyz_xxxy = pbuffer.data(idx_dip_hg + 556);

    auto tr_y_yyyyz_xxxz = pbuffer.data(idx_dip_hg + 557);

    auto tr_y_yyyyz_xxyy = pbuffer.data(idx_dip_hg + 558);

    auto tr_y_yyyyz_xxyz = pbuffer.data(idx_dip_hg + 559);

    auto tr_y_yyyyz_xxzz = pbuffer.data(idx_dip_hg + 560);

    auto tr_y_yyyyz_xyyy = pbuffer.data(idx_dip_hg + 561);

    auto tr_y_yyyyz_xyyz = pbuffer.data(idx_dip_hg + 562);

    auto tr_y_yyyyz_xyzz = pbuffer.data(idx_dip_hg + 563);

    auto tr_y_yyyyz_xzzz = pbuffer.data(idx_dip_hg + 564);

    auto tr_y_yyyyz_yyyy = pbuffer.data(idx_dip_hg + 565);

    auto tr_y_yyyyz_yyyz = pbuffer.data(idx_dip_hg + 566);

    auto tr_y_yyyyz_yyzz = pbuffer.data(idx_dip_hg + 567);

    auto tr_y_yyyyz_yzzz = pbuffer.data(idx_dip_hg + 568);

    auto tr_y_yyyyz_zzzz = pbuffer.data(idx_dip_hg + 569);

#pragma omp simd aligned(pa_z,                \
                             tr_y_yyyy_xxx,   \
                             tr_y_yyyy_xxxx,  \
                             tr_y_yyyy_xxxy,  \
                             tr_y_yyyy_xxxz,  \
                             tr_y_yyyy_xxy,   \
                             tr_y_yyyy_xxyy,  \
                             tr_y_yyyy_xxyz,  \
                             tr_y_yyyy_xxz,   \
                             tr_y_yyyy_xxzz,  \
                             tr_y_yyyy_xyy,   \
                             tr_y_yyyy_xyyy,  \
                             tr_y_yyyy_xyyz,  \
                             tr_y_yyyy_xyz,   \
                             tr_y_yyyy_xyzz,  \
                             tr_y_yyyy_xzz,   \
                             tr_y_yyyy_xzzz,  \
                             tr_y_yyyy_yyy,   \
                             tr_y_yyyy_yyyy,  \
                             tr_y_yyyy_yyyz,  \
                             tr_y_yyyy_yyz,   \
                             tr_y_yyyy_yyzz,  \
                             tr_y_yyyy_yzz,   \
                             tr_y_yyyy_yzzz,  \
                             tr_y_yyyy_zzz,   \
                             tr_y_yyyy_zzzz,  \
                             tr_y_yyyyz_xxxx, \
                             tr_y_yyyyz_xxxy, \
                             tr_y_yyyyz_xxxz, \
                             tr_y_yyyyz_xxyy, \
                             tr_y_yyyyz_xxyz, \
                             tr_y_yyyyz_xxzz, \
                             tr_y_yyyyz_xyyy, \
                             tr_y_yyyyz_xyyz, \
                             tr_y_yyyyz_xyzz, \
                             tr_y_yyyyz_xzzz, \
                             tr_y_yyyyz_yyyy, \
                             tr_y_yyyyz_yyyz, \
                             tr_y_yyyyz_yyzz, \
                             tr_y_yyyyz_yzzz, \
                             tr_y_yyyyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyz_xxxx[i] = tr_y_yyyy_xxxx[i] * pa_z[i];

        tr_y_yyyyz_xxxy[i] = tr_y_yyyy_xxxy[i] * pa_z[i];

        tr_y_yyyyz_xxxz[i] = tr_y_yyyy_xxx[i] * fe_0 + tr_y_yyyy_xxxz[i] * pa_z[i];

        tr_y_yyyyz_xxyy[i] = tr_y_yyyy_xxyy[i] * pa_z[i];

        tr_y_yyyyz_xxyz[i] = tr_y_yyyy_xxy[i] * fe_0 + tr_y_yyyy_xxyz[i] * pa_z[i];

        tr_y_yyyyz_xxzz[i] = 2.0 * tr_y_yyyy_xxz[i] * fe_0 + tr_y_yyyy_xxzz[i] * pa_z[i];

        tr_y_yyyyz_xyyy[i] = tr_y_yyyy_xyyy[i] * pa_z[i];

        tr_y_yyyyz_xyyz[i] = tr_y_yyyy_xyy[i] * fe_0 + tr_y_yyyy_xyyz[i] * pa_z[i];

        tr_y_yyyyz_xyzz[i] = 2.0 * tr_y_yyyy_xyz[i] * fe_0 + tr_y_yyyy_xyzz[i] * pa_z[i];

        tr_y_yyyyz_xzzz[i] = 3.0 * tr_y_yyyy_xzz[i] * fe_0 + tr_y_yyyy_xzzz[i] * pa_z[i];

        tr_y_yyyyz_yyyy[i] = tr_y_yyyy_yyyy[i] * pa_z[i];

        tr_y_yyyyz_yyyz[i] = tr_y_yyyy_yyy[i] * fe_0 + tr_y_yyyy_yyyz[i] * pa_z[i];

        tr_y_yyyyz_yyzz[i] = 2.0 * tr_y_yyyy_yyz[i] * fe_0 + tr_y_yyyy_yyzz[i] * pa_z[i];

        tr_y_yyyyz_yzzz[i] = 3.0 * tr_y_yyyy_yzz[i] * fe_0 + tr_y_yyyy_yzzz[i] * pa_z[i];

        tr_y_yyyyz_zzzz[i] = 4.0 * tr_y_yyyy_zzz[i] * fe_0 + tr_y_yyyy_zzzz[i] * pa_z[i];
    }

    // Set up 570-585 components of targeted buffer : HG

    auto tr_y_yyyzz_xxxx = pbuffer.data(idx_dip_hg + 570);

    auto tr_y_yyyzz_xxxy = pbuffer.data(idx_dip_hg + 571);

    auto tr_y_yyyzz_xxxz = pbuffer.data(idx_dip_hg + 572);

    auto tr_y_yyyzz_xxyy = pbuffer.data(idx_dip_hg + 573);

    auto tr_y_yyyzz_xxyz = pbuffer.data(idx_dip_hg + 574);

    auto tr_y_yyyzz_xxzz = pbuffer.data(idx_dip_hg + 575);

    auto tr_y_yyyzz_xyyy = pbuffer.data(idx_dip_hg + 576);

    auto tr_y_yyyzz_xyyz = pbuffer.data(idx_dip_hg + 577);

    auto tr_y_yyyzz_xyzz = pbuffer.data(idx_dip_hg + 578);

    auto tr_y_yyyzz_xzzz = pbuffer.data(idx_dip_hg + 579);

    auto tr_y_yyyzz_yyyy = pbuffer.data(idx_dip_hg + 580);

    auto tr_y_yyyzz_yyyz = pbuffer.data(idx_dip_hg + 581);

    auto tr_y_yyyzz_yyzz = pbuffer.data(idx_dip_hg + 582);

    auto tr_y_yyyzz_yzzz = pbuffer.data(idx_dip_hg + 583);

    auto tr_y_yyyzz_zzzz = pbuffer.data(idx_dip_hg + 584);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             tr_y_yyy_xxxx,   \
                             tr_y_yyy_xxxy,   \
                             tr_y_yyy_xxyy,   \
                             tr_y_yyy_xxyz,   \
                             tr_y_yyy_xyyy,   \
                             tr_y_yyy_xyyz,   \
                             tr_y_yyy_xyzz,   \
                             tr_y_yyy_yyyy,   \
                             tr_y_yyy_yyyz,   \
                             tr_y_yyy_yyzz,   \
                             tr_y_yyy_yzzz,   \
                             tr_y_yyyz_xxxx,  \
                             tr_y_yyyz_xxxy,  \
                             tr_y_yyyz_xxy,   \
                             tr_y_yyyz_xxyy,  \
                             tr_y_yyyz_xxyz,  \
                             tr_y_yyyz_xyy,   \
                             tr_y_yyyz_xyyy,  \
                             tr_y_yyyz_xyyz,  \
                             tr_y_yyyz_xyz,   \
                             tr_y_yyyz_xyzz,  \
                             tr_y_yyyz_yyy,   \
                             tr_y_yyyz_yyyy,  \
                             tr_y_yyyz_yyyz,  \
                             tr_y_yyyz_yyz,   \
                             tr_y_yyyz_yyzz,  \
                             tr_y_yyyz_yzz,   \
                             tr_y_yyyz_yzzz,  \
                             tr_y_yyyzz_xxxx, \
                             tr_y_yyyzz_xxxy, \
                             tr_y_yyyzz_xxxz, \
                             tr_y_yyyzz_xxyy, \
                             tr_y_yyyzz_xxyz, \
                             tr_y_yyyzz_xxzz, \
                             tr_y_yyyzz_xyyy, \
                             tr_y_yyyzz_xyyz, \
                             tr_y_yyyzz_xyzz, \
                             tr_y_yyyzz_xzzz, \
                             tr_y_yyyzz_yyyy, \
                             tr_y_yyyzz_yyyz, \
                             tr_y_yyyzz_yyzz, \
                             tr_y_yyyzz_yzzz, \
                             tr_y_yyyzz_zzzz, \
                             tr_y_yyzz_xxxz,  \
                             tr_y_yyzz_xxzz,  \
                             tr_y_yyzz_xzzz,  \
                             tr_y_yyzz_zzzz,  \
                             tr_y_yzz_xxxz,   \
                             tr_y_yzz_xxzz,   \
                             tr_y_yzz_xzzz,   \
                             tr_y_yzz_zzzz,   \
                             ts_yyzz_xxxz,    \
                             ts_yyzz_xxzz,    \
                             ts_yyzz_xzzz,    \
                             ts_yyzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyzz_xxxx[i] = tr_y_yyy_xxxx[i] * fe_0 + tr_y_yyyz_xxxx[i] * pa_z[i];

        tr_y_yyyzz_xxxy[i] = tr_y_yyy_xxxy[i] * fe_0 + tr_y_yyyz_xxxy[i] * pa_z[i];

        tr_y_yyyzz_xxxz[i] = 2.0 * tr_y_yzz_xxxz[i] * fe_0 + ts_yyzz_xxxz[i] * fe_0 + tr_y_yyzz_xxxz[i] * pa_y[i];

        tr_y_yyyzz_xxyy[i] = tr_y_yyy_xxyy[i] * fe_0 + tr_y_yyyz_xxyy[i] * pa_z[i];

        tr_y_yyyzz_xxyz[i] = tr_y_yyy_xxyz[i] * fe_0 + tr_y_yyyz_xxy[i] * fe_0 + tr_y_yyyz_xxyz[i] * pa_z[i];

        tr_y_yyyzz_xxzz[i] = 2.0 * tr_y_yzz_xxzz[i] * fe_0 + ts_yyzz_xxzz[i] * fe_0 + tr_y_yyzz_xxzz[i] * pa_y[i];

        tr_y_yyyzz_xyyy[i] = tr_y_yyy_xyyy[i] * fe_0 + tr_y_yyyz_xyyy[i] * pa_z[i];

        tr_y_yyyzz_xyyz[i] = tr_y_yyy_xyyz[i] * fe_0 + tr_y_yyyz_xyy[i] * fe_0 + tr_y_yyyz_xyyz[i] * pa_z[i];

        tr_y_yyyzz_xyzz[i] = tr_y_yyy_xyzz[i] * fe_0 + 2.0 * tr_y_yyyz_xyz[i] * fe_0 + tr_y_yyyz_xyzz[i] * pa_z[i];

        tr_y_yyyzz_xzzz[i] = 2.0 * tr_y_yzz_xzzz[i] * fe_0 + ts_yyzz_xzzz[i] * fe_0 + tr_y_yyzz_xzzz[i] * pa_y[i];

        tr_y_yyyzz_yyyy[i] = tr_y_yyy_yyyy[i] * fe_0 + tr_y_yyyz_yyyy[i] * pa_z[i];

        tr_y_yyyzz_yyyz[i] = tr_y_yyy_yyyz[i] * fe_0 + tr_y_yyyz_yyy[i] * fe_0 + tr_y_yyyz_yyyz[i] * pa_z[i];

        tr_y_yyyzz_yyzz[i] = tr_y_yyy_yyzz[i] * fe_0 + 2.0 * tr_y_yyyz_yyz[i] * fe_0 + tr_y_yyyz_yyzz[i] * pa_z[i];

        tr_y_yyyzz_yzzz[i] = tr_y_yyy_yzzz[i] * fe_0 + 3.0 * tr_y_yyyz_yzz[i] * fe_0 + tr_y_yyyz_yzzz[i] * pa_z[i];

        tr_y_yyyzz_zzzz[i] = 2.0 * tr_y_yzz_zzzz[i] * fe_0 + ts_yyzz_zzzz[i] * fe_0 + tr_y_yyzz_zzzz[i] * pa_y[i];
    }

    // Set up 585-600 components of targeted buffer : HG

    auto tr_y_yyzzz_xxxx = pbuffer.data(idx_dip_hg + 585);

    auto tr_y_yyzzz_xxxy = pbuffer.data(idx_dip_hg + 586);

    auto tr_y_yyzzz_xxxz = pbuffer.data(idx_dip_hg + 587);

    auto tr_y_yyzzz_xxyy = pbuffer.data(idx_dip_hg + 588);

    auto tr_y_yyzzz_xxyz = pbuffer.data(idx_dip_hg + 589);

    auto tr_y_yyzzz_xxzz = pbuffer.data(idx_dip_hg + 590);

    auto tr_y_yyzzz_xyyy = pbuffer.data(idx_dip_hg + 591);

    auto tr_y_yyzzz_xyyz = pbuffer.data(idx_dip_hg + 592);

    auto tr_y_yyzzz_xyzz = pbuffer.data(idx_dip_hg + 593);

    auto tr_y_yyzzz_xzzz = pbuffer.data(idx_dip_hg + 594);

    auto tr_y_yyzzz_yyyy = pbuffer.data(idx_dip_hg + 595);

    auto tr_y_yyzzz_yyyz = pbuffer.data(idx_dip_hg + 596);

    auto tr_y_yyzzz_yyzz = pbuffer.data(idx_dip_hg + 597);

    auto tr_y_yyzzz_yzzz = pbuffer.data(idx_dip_hg + 598);

    auto tr_y_yyzzz_zzzz = pbuffer.data(idx_dip_hg + 599);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             tr_y_yyz_xxxx,   \
                             tr_y_yyz_xxxy,   \
                             tr_y_yyz_xxyy,   \
                             tr_y_yyz_xxyz,   \
                             tr_y_yyz_xyyy,   \
                             tr_y_yyz_xyyz,   \
                             tr_y_yyz_xyzz,   \
                             tr_y_yyz_yyyy,   \
                             tr_y_yyz_yyyz,   \
                             tr_y_yyz_yyzz,   \
                             tr_y_yyz_yzzz,   \
                             tr_y_yyzz_xxxx,  \
                             tr_y_yyzz_xxxy,  \
                             tr_y_yyzz_xxy,   \
                             tr_y_yyzz_xxyy,  \
                             tr_y_yyzz_xxyz,  \
                             tr_y_yyzz_xyy,   \
                             tr_y_yyzz_xyyy,  \
                             tr_y_yyzz_xyyz,  \
                             tr_y_yyzz_xyz,   \
                             tr_y_yyzz_xyzz,  \
                             tr_y_yyzz_yyy,   \
                             tr_y_yyzz_yyyy,  \
                             tr_y_yyzz_yyyz,  \
                             tr_y_yyzz_yyz,   \
                             tr_y_yyzz_yyzz,  \
                             tr_y_yyzz_yzz,   \
                             tr_y_yyzz_yzzz,  \
                             tr_y_yyzzz_xxxx, \
                             tr_y_yyzzz_xxxy, \
                             tr_y_yyzzz_xxxz, \
                             tr_y_yyzzz_xxyy, \
                             tr_y_yyzzz_xxyz, \
                             tr_y_yyzzz_xxzz, \
                             tr_y_yyzzz_xyyy, \
                             tr_y_yyzzz_xyyz, \
                             tr_y_yyzzz_xyzz, \
                             tr_y_yyzzz_xzzz, \
                             tr_y_yyzzz_yyyy, \
                             tr_y_yyzzz_yyyz, \
                             tr_y_yyzzz_yyzz, \
                             tr_y_yyzzz_yzzz, \
                             tr_y_yyzzz_zzzz, \
                             tr_y_yzzz_xxxz,  \
                             tr_y_yzzz_xxzz,  \
                             tr_y_yzzz_xzzz,  \
                             tr_y_yzzz_zzzz,  \
                             tr_y_zzz_xxxz,   \
                             tr_y_zzz_xxzz,   \
                             tr_y_zzz_xzzz,   \
                             tr_y_zzz_zzzz,   \
                             ts_yzzz_xxxz,    \
                             ts_yzzz_xxzz,    \
                             ts_yzzz_xzzz,    \
                             ts_yzzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzzz_xxxx[i] = 2.0 * tr_y_yyz_xxxx[i] * fe_0 + tr_y_yyzz_xxxx[i] * pa_z[i];

        tr_y_yyzzz_xxxy[i] = 2.0 * tr_y_yyz_xxxy[i] * fe_0 + tr_y_yyzz_xxxy[i] * pa_z[i];

        tr_y_yyzzz_xxxz[i] = tr_y_zzz_xxxz[i] * fe_0 + ts_yzzz_xxxz[i] * fe_0 + tr_y_yzzz_xxxz[i] * pa_y[i];

        tr_y_yyzzz_xxyy[i] = 2.0 * tr_y_yyz_xxyy[i] * fe_0 + tr_y_yyzz_xxyy[i] * pa_z[i];

        tr_y_yyzzz_xxyz[i] = 2.0 * tr_y_yyz_xxyz[i] * fe_0 + tr_y_yyzz_xxy[i] * fe_0 + tr_y_yyzz_xxyz[i] * pa_z[i];

        tr_y_yyzzz_xxzz[i] = tr_y_zzz_xxzz[i] * fe_0 + ts_yzzz_xxzz[i] * fe_0 + tr_y_yzzz_xxzz[i] * pa_y[i];

        tr_y_yyzzz_xyyy[i] = 2.0 * tr_y_yyz_xyyy[i] * fe_0 + tr_y_yyzz_xyyy[i] * pa_z[i];

        tr_y_yyzzz_xyyz[i] = 2.0 * tr_y_yyz_xyyz[i] * fe_0 + tr_y_yyzz_xyy[i] * fe_0 + tr_y_yyzz_xyyz[i] * pa_z[i];

        tr_y_yyzzz_xyzz[i] = 2.0 * tr_y_yyz_xyzz[i] * fe_0 + 2.0 * tr_y_yyzz_xyz[i] * fe_0 + tr_y_yyzz_xyzz[i] * pa_z[i];

        tr_y_yyzzz_xzzz[i] = tr_y_zzz_xzzz[i] * fe_0 + ts_yzzz_xzzz[i] * fe_0 + tr_y_yzzz_xzzz[i] * pa_y[i];

        tr_y_yyzzz_yyyy[i] = 2.0 * tr_y_yyz_yyyy[i] * fe_0 + tr_y_yyzz_yyyy[i] * pa_z[i];

        tr_y_yyzzz_yyyz[i] = 2.0 * tr_y_yyz_yyyz[i] * fe_0 + tr_y_yyzz_yyy[i] * fe_0 + tr_y_yyzz_yyyz[i] * pa_z[i];

        tr_y_yyzzz_yyzz[i] = 2.0 * tr_y_yyz_yyzz[i] * fe_0 + 2.0 * tr_y_yyzz_yyz[i] * fe_0 + tr_y_yyzz_yyzz[i] * pa_z[i];

        tr_y_yyzzz_yzzz[i] = 2.0 * tr_y_yyz_yzzz[i] * fe_0 + 3.0 * tr_y_yyzz_yzz[i] * fe_0 + tr_y_yyzz_yzzz[i] * pa_z[i];

        tr_y_yyzzz_zzzz[i] = tr_y_zzz_zzzz[i] * fe_0 + ts_yzzz_zzzz[i] * fe_0 + tr_y_yzzz_zzzz[i] * pa_y[i];
    }

    // Set up 600-615 components of targeted buffer : HG

    auto tr_y_yzzzz_xxxx = pbuffer.data(idx_dip_hg + 600);

    auto tr_y_yzzzz_xxxy = pbuffer.data(idx_dip_hg + 601);

    auto tr_y_yzzzz_xxxz = pbuffer.data(idx_dip_hg + 602);

    auto tr_y_yzzzz_xxyy = pbuffer.data(idx_dip_hg + 603);

    auto tr_y_yzzzz_xxyz = pbuffer.data(idx_dip_hg + 604);

    auto tr_y_yzzzz_xxzz = pbuffer.data(idx_dip_hg + 605);

    auto tr_y_yzzzz_xyyy = pbuffer.data(idx_dip_hg + 606);

    auto tr_y_yzzzz_xyyz = pbuffer.data(idx_dip_hg + 607);

    auto tr_y_yzzzz_xyzz = pbuffer.data(idx_dip_hg + 608);

    auto tr_y_yzzzz_xzzz = pbuffer.data(idx_dip_hg + 609);

    auto tr_y_yzzzz_yyyy = pbuffer.data(idx_dip_hg + 610);

    auto tr_y_yzzzz_yyyz = pbuffer.data(idx_dip_hg + 611);

    auto tr_y_yzzzz_yyzz = pbuffer.data(idx_dip_hg + 612);

    auto tr_y_yzzzz_yzzz = pbuffer.data(idx_dip_hg + 613);

    auto tr_y_yzzzz_zzzz = pbuffer.data(idx_dip_hg + 614);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             tr_y_yzz_xxxy,   \
                             tr_y_yzz_xxyy,   \
                             tr_y_yzz_xyyy,   \
                             tr_y_yzz_yyyy,   \
                             tr_y_yzzz_xxxy,  \
                             tr_y_yzzz_xxyy,  \
                             tr_y_yzzz_xyyy,  \
                             tr_y_yzzz_yyyy,  \
                             tr_y_yzzzz_xxxx, \
                             tr_y_yzzzz_xxxy, \
                             tr_y_yzzzz_xxxz, \
                             tr_y_yzzzz_xxyy, \
                             tr_y_yzzzz_xxyz, \
                             tr_y_yzzzz_xxzz, \
                             tr_y_yzzzz_xyyy, \
                             tr_y_yzzzz_xyyz, \
                             tr_y_yzzzz_xyzz, \
                             tr_y_yzzzz_xzzz, \
                             tr_y_yzzzz_yyyy, \
                             tr_y_yzzzz_yyyz, \
                             tr_y_yzzzz_yyzz, \
                             tr_y_yzzzz_yzzz, \
                             tr_y_yzzzz_zzzz, \
                             tr_y_zzzz_xxxx,  \
                             tr_y_zzzz_xxxz,  \
                             tr_y_zzzz_xxyz,  \
                             tr_y_zzzz_xxz,   \
                             tr_y_zzzz_xxzz,  \
                             tr_y_zzzz_xyyz,  \
                             tr_y_zzzz_xyz,   \
                             tr_y_zzzz_xyzz,  \
                             tr_y_zzzz_xzz,   \
                             tr_y_zzzz_xzzz,  \
                             tr_y_zzzz_yyyz,  \
                             tr_y_zzzz_yyz,   \
                             tr_y_zzzz_yyzz,  \
                             tr_y_zzzz_yzz,   \
                             tr_y_zzzz_yzzz,  \
                             tr_y_zzzz_zzz,   \
                             tr_y_zzzz_zzzz,  \
                             ts_zzzz_xxxx,    \
                             ts_zzzz_xxxz,    \
                             ts_zzzz_xxyz,    \
                             ts_zzzz_xxzz,    \
                             ts_zzzz_xyyz,    \
                             ts_zzzz_xyzz,    \
                             ts_zzzz_xzzz,    \
                             ts_zzzz_yyyz,    \
                             ts_zzzz_yyzz,    \
                             ts_zzzz_yzzz,    \
                             ts_zzzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzzz_xxxx[i] = ts_zzzz_xxxx[i] * fe_0 + tr_y_zzzz_xxxx[i] * pa_y[i];

        tr_y_yzzzz_xxxy[i] = 3.0 * tr_y_yzz_xxxy[i] * fe_0 + tr_y_yzzz_xxxy[i] * pa_z[i];

        tr_y_yzzzz_xxxz[i] = ts_zzzz_xxxz[i] * fe_0 + tr_y_zzzz_xxxz[i] * pa_y[i];

        tr_y_yzzzz_xxyy[i] = 3.0 * tr_y_yzz_xxyy[i] * fe_0 + tr_y_yzzz_xxyy[i] * pa_z[i];

        tr_y_yzzzz_xxyz[i] = tr_y_zzzz_xxz[i] * fe_0 + ts_zzzz_xxyz[i] * fe_0 + tr_y_zzzz_xxyz[i] * pa_y[i];

        tr_y_yzzzz_xxzz[i] = ts_zzzz_xxzz[i] * fe_0 + tr_y_zzzz_xxzz[i] * pa_y[i];

        tr_y_yzzzz_xyyy[i] = 3.0 * tr_y_yzz_xyyy[i] * fe_0 + tr_y_yzzz_xyyy[i] * pa_z[i];

        tr_y_yzzzz_xyyz[i] = 2.0 * tr_y_zzzz_xyz[i] * fe_0 + ts_zzzz_xyyz[i] * fe_0 + tr_y_zzzz_xyyz[i] * pa_y[i];

        tr_y_yzzzz_xyzz[i] = tr_y_zzzz_xzz[i] * fe_0 + ts_zzzz_xyzz[i] * fe_0 + tr_y_zzzz_xyzz[i] * pa_y[i];

        tr_y_yzzzz_xzzz[i] = ts_zzzz_xzzz[i] * fe_0 + tr_y_zzzz_xzzz[i] * pa_y[i];

        tr_y_yzzzz_yyyy[i] = 3.0 * tr_y_yzz_yyyy[i] * fe_0 + tr_y_yzzz_yyyy[i] * pa_z[i];

        tr_y_yzzzz_yyyz[i] = 3.0 * tr_y_zzzz_yyz[i] * fe_0 + ts_zzzz_yyyz[i] * fe_0 + tr_y_zzzz_yyyz[i] * pa_y[i];

        tr_y_yzzzz_yyzz[i] = 2.0 * tr_y_zzzz_yzz[i] * fe_0 + ts_zzzz_yyzz[i] * fe_0 + tr_y_zzzz_yyzz[i] * pa_y[i];

        tr_y_yzzzz_yzzz[i] = tr_y_zzzz_zzz[i] * fe_0 + ts_zzzz_yzzz[i] * fe_0 + tr_y_zzzz_yzzz[i] * pa_y[i];

        tr_y_yzzzz_zzzz[i] = ts_zzzz_zzzz[i] * fe_0 + tr_y_zzzz_zzzz[i] * pa_y[i];
    }

    // Set up 615-630 components of targeted buffer : HG

    auto tr_y_zzzzz_xxxx = pbuffer.data(idx_dip_hg + 615);

    auto tr_y_zzzzz_xxxy = pbuffer.data(idx_dip_hg + 616);

    auto tr_y_zzzzz_xxxz = pbuffer.data(idx_dip_hg + 617);

    auto tr_y_zzzzz_xxyy = pbuffer.data(idx_dip_hg + 618);

    auto tr_y_zzzzz_xxyz = pbuffer.data(idx_dip_hg + 619);

    auto tr_y_zzzzz_xxzz = pbuffer.data(idx_dip_hg + 620);

    auto tr_y_zzzzz_xyyy = pbuffer.data(idx_dip_hg + 621);

    auto tr_y_zzzzz_xyyz = pbuffer.data(idx_dip_hg + 622);

    auto tr_y_zzzzz_xyzz = pbuffer.data(idx_dip_hg + 623);

    auto tr_y_zzzzz_xzzz = pbuffer.data(idx_dip_hg + 624);

    auto tr_y_zzzzz_yyyy = pbuffer.data(idx_dip_hg + 625);

    auto tr_y_zzzzz_yyyz = pbuffer.data(idx_dip_hg + 626);

    auto tr_y_zzzzz_yyzz = pbuffer.data(idx_dip_hg + 627);

    auto tr_y_zzzzz_yzzz = pbuffer.data(idx_dip_hg + 628);

    auto tr_y_zzzzz_zzzz = pbuffer.data(idx_dip_hg + 629);

#pragma omp simd aligned(pa_z,                \
                             tr_y_zzz_xxxx,   \
                             tr_y_zzz_xxxy,   \
                             tr_y_zzz_xxxz,   \
                             tr_y_zzz_xxyy,   \
                             tr_y_zzz_xxyz,   \
                             tr_y_zzz_xxzz,   \
                             tr_y_zzz_xyyy,   \
                             tr_y_zzz_xyyz,   \
                             tr_y_zzz_xyzz,   \
                             tr_y_zzz_xzzz,   \
                             tr_y_zzz_yyyy,   \
                             tr_y_zzz_yyyz,   \
                             tr_y_zzz_yyzz,   \
                             tr_y_zzz_yzzz,   \
                             tr_y_zzz_zzzz,   \
                             tr_y_zzzz_xxx,   \
                             tr_y_zzzz_xxxx,  \
                             tr_y_zzzz_xxxy,  \
                             tr_y_zzzz_xxxz,  \
                             tr_y_zzzz_xxy,   \
                             tr_y_zzzz_xxyy,  \
                             tr_y_zzzz_xxyz,  \
                             tr_y_zzzz_xxz,   \
                             tr_y_zzzz_xxzz,  \
                             tr_y_zzzz_xyy,   \
                             tr_y_zzzz_xyyy,  \
                             tr_y_zzzz_xyyz,  \
                             tr_y_zzzz_xyz,   \
                             tr_y_zzzz_xyzz,  \
                             tr_y_zzzz_xzz,   \
                             tr_y_zzzz_xzzz,  \
                             tr_y_zzzz_yyy,   \
                             tr_y_zzzz_yyyy,  \
                             tr_y_zzzz_yyyz,  \
                             tr_y_zzzz_yyz,   \
                             tr_y_zzzz_yyzz,  \
                             tr_y_zzzz_yzz,   \
                             tr_y_zzzz_yzzz,  \
                             tr_y_zzzz_zzz,   \
                             tr_y_zzzz_zzzz,  \
                             tr_y_zzzzz_xxxx, \
                             tr_y_zzzzz_xxxy, \
                             tr_y_zzzzz_xxxz, \
                             tr_y_zzzzz_xxyy, \
                             tr_y_zzzzz_xxyz, \
                             tr_y_zzzzz_xxzz, \
                             tr_y_zzzzz_xyyy, \
                             tr_y_zzzzz_xyyz, \
                             tr_y_zzzzz_xyzz, \
                             tr_y_zzzzz_xzzz, \
                             tr_y_zzzzz_yyyy, \
                             tr_y_zzzzz_yyyz, \
                             tr_y_zzzzz_yyzz, \
                             tr_y_zzzzz_yzzz, \
                             tr_y_zzzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzzz_xxxx[i] = 4.0 * tr_y_zzz_xxxx[i] * fe_0 + tr_y_zzzz_xxxx[i] * pa_z[i];

        tr_y_zzzzz_xxxy[i] = 4.0 * tr_y_zzz_xxxy[i] * fe_0 + tr_y_zzzz_xxxy[i] * pa_z[i];

        tr_y_zzzzz_xxxz[i] = 4.0 * tr_y_zzz_xxxz[i] * fe_0 + tr_y_zzzz_xxx[i] * fe_0 + tr_y_zzzz_xxxz[i] * pa_z[i];

        tr_y_zzzzz_xxyy[i] = 4.0 * tr_y_zzz_xxyy[i] * fe_0 + tr_y_zzzz_xxyy[i] * pa_z[i];

        tr_y_zzzzz_xxyz[i] = 4.0 * tr_y_zzz_xxyz[i] * fe_0 + tr_y_zzzz_xxy[i] * fe_0 + tr_y_zzzz_xxyz[i] * pa_z[i];

        tr_y_zzzzz_xxzz[i] = 4.0 * tr_y_zzz_xxzz[i] * fe_0 + 2.0 * tr_y_zzzz_xxz[i] * fe_0 + tr_y_zzzz_xxzz[i] * pa_z[i];

        tr_y_zzzzz_xyyy[i] = 4.0 * tr_y_zzz_xyyy[i] * fe_0 + tr_y_zzzz_xyyy[i] * pa_z[i];

        tr_y_zzzzz_xyyz[i] = 4.0 * tr_y_zzz_xyyz[i] * fe_0 + tr_y_zzzz_xyy[i] * fe_0 + tr_y_zzzz_xyyz[i] * pa_z[i];

        tr_y_zzzzz_xyzz[i] = 4.0 * tr_y_zzz_xyzz[i] * fe_0 + 2.0 * tr_y_zzzz_xyz[i] * fe_0 + tr_y_zzzz_xyzz[i] * pa_z[i];

        tr_y_zzzzz_xzzz[i] = 4.0 * tr_y_zzz_xzzz[i] * fe_0 + 3.0 * tr_y_zzzz_xzz[i] * fe_0 + tr_y_zzzz_xzzz[i] * pa_z[i];

        tr_y_zzzzz_yyyy[i] = 4.0 * tr_y_zzz_yyyy[i] * fe_0 + tr_y_zzzz_yyyy[i] * pa_z[i];

        tr_y_zzzzz_yyyz[i] = 4.0 * tr_y_zzz_yyyz[i] * fe_0 + tr_y_zzzz_yyy[i] * fe_0 + tr_y_zzzz_yyyz[i] * pa_z[i];

        tr_y_zzzzz_yyzz[i] = 4.0 * tr_y_zzz_yyzz[i] * fe_0 + 2.0 * tr_y_zzzz_yyz[i] * fe_0 + tr_y_zzzz_yyzz[i] * pa_z[i];

        tr_y_zzzzz_yzzz[i] = 4.0 * tr_y_zzz_yzzz[i] * fe_0 + 3.0 * tr_y_zzzz_yzz[i] * fe_0 + tr_y_zzzz_yzzz[i] * pa_z[i];

        tr_y_zzzzz_zzzz[i] = 4.0 * tr_y_zzz_zzzz[i] * fe_0 + 4.0 * tr_y_zzzz_zzz[i] * fe_0 + tr_y_zzzz_zzzz[i] * pa_z[i];
    }

    // Set up 630-645 components of targeted buffer : HG

    auto tr_z_xxxxx_xxxx = pbuffer.data(idx_dip_hg + 630);

    auto tr_z_xxxxx_xxxy = pbuffer.data(idx_dip_hg + 631);

    auto tr_z_xxxxx_xxxz = pbuffer.data(idx_dip_hg + 632);

    auto tr_z_xxxxx_xxyy = pbuffer.data(idx_dip_hg + 633);

    auto tr_z_xxxxx_xxyz = pbuffer.data(idx_dip_hg + 634);

    auto tr_z_xxxxx_xxzz = pbuffer.data(idx_dip_hg + 635);

    auto tr_z_xxxxx_xyyy = pbuffer.data(idx_dip_hg + 636);

    auto tr_z_xxxxx_xyyz = pbuffer.data(idx_dip_hg + 637);

    auto tr_z_xxxxx_xyzz = pbuffer.data(idx_dip_hg + 638);

    auto tr_z_xxxxx_xzzz = pbuffer.data(idx_dip_hg + 639);

    auto tr_z_xxxxx_yyyy = pbuffer.data(idx_dip_hg + 640);

    auto tr_z_xxxxx_yyyz = pbuffer.data(idx_dip_hg + 641);

    auto tr_z_xxxxx_yyzz = pbuffer.data(idx_dip_hg + 642);

    auto tr_z_xxxxx_yzzz = pbuffer.data(idx_dip_hg + 643);

    auto tr_z_xxxxx_zzzz = pbuffer.data(idx_dip_hg + 644);

#pragma omp simd aligned(pa_x,                \
                             tr_z_xxx_xxxx,   \
                             tr_z_xxx_xxxy,   \
                             tr_z_xxx_xxxz,   \
                             tr_z_xxx_xxyy,   \
                             tr_z_xxx_xxyz,   \
                             tr_z_xxx_xxzz,   \
                             tr_z_xxx_xyyy,   \
                             tr_z_xxx_xyyz,   \
                             tr_z_xxx_xyzz,   \
                             tr_z_xxx_xzzz,   \
                             tr_z_xxx_yyyy,   \
                             tr_z_xxx_yyyz,   \
                             tr_z_xxx_yyzz,   \
                             tr_z_xxx_yzzz,   \
                             tr_z_xxx_zzzz,   \
                             tr_z_xxxx_xxx,   \
                             tr_z_xxxx_xxxx,  \
                             tr_z_xxxx_xxxy,  \
                             tr_z_xxxx_xxxz,  \
                             tr_z_xxxx_xxy,   \
                             tr_z_xxxx_xxyy,  \
                             tr_z_xxxx_xxyz,  \
                             tr_z_xxxx_xxz,   \
                             tr_z_xxxx_xxzz,  \
                             tr_z_xxxx_xyy,   \
                             tr_z_xxxx_xyyy,  \
                             tr_z_xxxx_xyyz,  \
                             tr_z_xxxx_xyz,   \
                             tr_z_xxxx_xyzz,  \
                             tr_z_xxxx_xzz,   \
                             tr_z_xxxx_xzzz,  \
                             tr_z_xxxx_yyy,   \
                             tr_z_xxxx_yyyy,  \
                             tr_z_xxxx_yyyz,  \
                             tr_z_xxxx_yyz,   \
                             tr_z_xxxx_yyzz,  \
                             tr_z_xxxx_yzz,   \
                             tr_z_xxxx_yzzz,  \
                             tr_z_xxxx_zzz,   \
                             tr_z_xxxx_zzzz,  \
                             tr_z_xxxxx_xxxx, \
                             tr_z_xxxxx_xxxy, \
                             tr_z_xxxxx_xxxz, \
                             tr_z_xxxxx_xxyy, \
                             tr_z_xxxxx_xxyz, \
                             tr_z_xxxxx_xxzz, \
                             tr_z_xxxxx_xyyy, \
                             tr_z_xxxxx_xyyz, \
                             tr_z_xxxxx_xyzz, \
                             tr_z_xxxxx_xzzz, \
                             tr_z_xxxxx_yyyy, \
                             tr_z_xxxxx_yyyz, \
                             tr_z_xxxxx_yyzz, \
                             tr_z_xxxxx_yzzz, \
                             tr_z_xxxxx_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxx_xxxx[i] = 4.0 * tr_z_xxx_xxxx[i] * fe_0 + 4.0 * tr_z_xxxx_xxx[i] * fe_0 + tr_z_xxxx_xxxx[i] * pa_x[i];

        tr_z_xxxxx_xxxy[i] = 4.0 * tr_z_xxx_xxxy[i] * fe_0 + 3.0 * tr_z_xxxx_xxy[i] * fe_0 + tr_z_xxxx_xxxy[i] * pa_x[i];

        tr_z_xxxxx_xxxz[i] = 4.0 * tr_z_xxx_xxxz[i] * fe_0 + 3.0 * tr_z_xxxx_xxz[i] * fe_0 + tr_z_xxxx_xxxz[i] * pa_x[i];

        tr_z_xxxxx_xxyy[i] = 4.0 * tr_z_xxx_xxyy[i] * fe_0 + 2.0 * tr_z_xxxx_xyy[i] * fe_0 + tr_z_xxxx_xxyy[i] * pa_x[i];

        tr_z_xxxxx_xxyz[i] = 4.0 * tr_z_xxx_xxyz[i] * fe_0 + 2.0 * tr_z_xxxx_xyz[i] * fe_0 + tr_z_xxxx_xxyz[i] * pa_x[i];

        tr_z_xxxxx_xxzz[i] = 4.0 * tr_z_xxx_xxzz[i] * fe_0 + 2.0 * tr_z_xxxx_xzz[i] * fe_0 + tr_z_xxxx_xxzz[i] * pa_x[i];

        tr_z_xxxxx_xyyy[i] = 4.0 * tr_z_xxx_xyyy[i] * fe_0 + tr_z_xxxx_yyy[i] * fe_0 + tr_z_xxxx_xyyy[i] * pa_x[i];

        tr_z_xxxxx_xyyz[i] = 4.0 * tr_z_xxx_xyyz[i] * fe_0 + tr_z_xxxx_yyz[i] * fe_0 + tr_z_xxxx_xyyz[i] * pa_x[i];

        tr_z_xxxxx_xyzz[i] = 4.0 * tr_z_xxx_xyzz[i] * fe_0 + tr_z_xxxx_yzz[i] * fe_0 + tr_z_xxxx_xyzz[i] * pa_x[i];

        tr_z_xxxxx_xzzz[i] = 4.0 * tr_z_xxx_xzzz[i] * fe_0 + tr_z_xxxx_zzz[i] * fe_0 + tr_z_xxxx_xzzz[i] * pa_x[i];

        tr_z_xxxxx_yyyy[i] = 4.0 * tr_z_xxx_yyyy[i] * fe_0 + tr_z_xxxx_yyyy[i] * pa_x[i];

        tr_z_xxxxx_yyyz[i] = 4.0 * tr_z_xxx_yyyz[i] * fe_0 + tr_z_xxxx_yyyz[i] * pa_x[i];

        tr_z_xxxxx_yyzz[i] = 4.0 * tr_z_xxx_yyzz[i] * fe_0 + tr_z_xxxx_yyzz[i] * pa_x[i];

        tr_z_xxxxx_yzzz[i] = 4.0 * tr_z_xxx_yzzz[i] * fe_0 + tr_z_xxxx_yzzz[i] * pa_x[i];

        tr_z_xxxxx_zzzz[i] = 4.0 * tr_z_xxx_zzzz[i] * fe_0 + tr_z_xxxx_zzzz[i] * pa_x[i];
    }

    // Set up 645-660 components of targeted buffer : HG

    auto tr_z_xxxxy_xxxx = pbuffer.data(idx_dip_hg + 645);

    auto tr_z_xxxxy_xxxy = pbuffer.data(idx_dip_hg + 646);

    auto tr_z_xxxxy_xxxz = pbuffer.data(idx_dip_hg + 647);

    auto tr_z_xxxxy_xxyy = pbuffer.data(idx_dip_hg + 648);

    auto tr_z_xxxxy_xxyz = pbuffer.data(idx_dip_hg + 649);

    auto tr_z_xxxxy_xxzz = pbuffer.data(idx_dip_hg + 650);

    auto tr_z_xxxxy_xyyy = pbuffer.data(idx_dip_hg + 651);

    auto tr_z_xxxxy_xyyz = pbuffer.data(idx_dip_hg + 652);

    auto tr_z_xxxxy_xyzz = pbuffer.data(idx_dip_hg + 653);

    auto tr_z_xxxxy_xzzz = pbuffer.data(idx_dip_hg + 654);

    auto tr_z_xxxxy_yyyy = pbuffer.data(idx_dip_hg + 655);

    auto tr_z_xxxxy_yyyz = pbuffer.data(idx_dip_hg + 656);

    auto tr_z_xxxxy_yyzz = pbuffer.data(idx_dip_hg + 657);

    auto tr_z_xxxxy_yzzz = pbuffer.data(idx_dip_hg + 658);

    auto tr_z_xxxxy_zzzz = pbuffer.data(idx_dip_hg + 659);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tr_z_xxxx_xxx,   \
                             tr_z_xxxx_xxxx,  \
                             tr_z_xxxx_xxxy,  \
                             tr_z_xxxx_xxxz,  \
                             tr_z_xxxx_xxy,   \
                             tr_z_xxxx_xxyy,  \
                             tr_z_xxxx_xxyz,  \
                             tr_z_xxxx_xxz,   \
                             tr_z_xxxx_xxzz,  \
                             tr_z_xxxx_xyy,   \
                             tr_z_xxxx_xyyy,  \
                             tr_z_xxxx_xyyz,  \
                             tr_z_xxxx_xyz,   \
                             tr_z_xxxx_xyzz,  \
                             tr_z_xxxx_xzz,   \
                             tr_z_xxxx_xzzz,  \
                             tr_z_xxxx_zzzz,  \
                             tr_z_xxxxy_xxxx, \
                             tr_z_xxxxy_xxxy, \
                             tr_z_xxxxy_xxxz, \
                             tr_z_xxxxy_xxyy, \
                             tr_z_xxxxy_xxyz, \
                             tr_z_xxxxy_xxzz, \
                             tr_z_xxxxy_xyyy, \
                             tr_z_xxxxy_xyyz, \
                             tr_z_xxxxy_xyzz, \
                             tr_z_xxxxy_xzzz, \
                             tr_z_xxxxy_yyyy, \
                             tr_z_xxxxy_yyyz, \
                             tr_z_xxxxy_yyzz, \
                             tr_z_xxxxy_yzzz, \
                             tr_z_xxxxy_zzzz, \
                             tr_z_xxxy_yyyy,  \
                             tr_z_xxxy_yyyz,  \
                             tr_z_xxxy_yyzz,  \
                             tr_z_xxxy_yzzz,  \
                             tr_z_xxy_yyyy,   \
                             tr_z_xxy_yyyz,   \
                             tr_z_xxy_yyzz,   \
                             tr_z_xxy_yzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxy_xxxx[i] = tr_z_xxxx_xxxx[i] * pa_y[i];

        tr_z_xxxxy_xxxy[i] = tr_z_xxxx_xxx[i] * fe_0 + tr_z_xxxx_xxxy[i] * pa_y[i];

        tr_z_xxxxy_xxxz[i] = tr_z_xxxx_xxxz[i] * pa_y[i];

        tr_z_xxxxy_xxyy[i] = 2.0 * tr_z_xxxx_xxy[i] * fe_0 + tr_z_xxxx_xxyy[i] * pa_y[i];

        tr_z_xxxxy_xxyz[i] = tr_z_xxxx_xxz[i] * fe_0 + tr_z_xxxx_xxyz[i] * pa_y[i];

        tr_z_xxxxy_xxzz[i] = tr_z_xxxx_xxzz[i] * pa_y[i];

        tr_z_xxxxy_xyyy[i] = 3.0 * tr_z_xxxx_xyy[i] * fe_0 + tr_z_xxxx_xyyy[i] * pa_y[i];

        tr_z_xxxxy_xyyz[i] = 2.0 * tr_z_xxxx_xyz[i] * fe_0 + tr_z_xxxx_xyyz[i] * pa_y[i];

        tr_z_xxxxy_xyzz[i] = tr_z_xxxx_xzz[i] * fe_0 + tr_z_xxxx_xyzz[i] * pa_y[i];

        tr_z_xxxxy_xzzz[i] = tr_z_xxxx_xzzz[i] * pa_y[i];

        tr_z_xxxxy_yyyy[i] = 3.0 * tr_z_xxy_yyyy[i] * fe_0 + tr_z_xxxy_yyyy[i] * pa_x[i];

        tr_z_xxxxy_yyyz[i] = 3.0 * tr_z_xxy_yyyz[i] * fe_0 + tr_z_xxxy_yyyz[i] * pa_x[i];

        tr_z_xxxxy_yyzz[i] = 3.0 * tr_z_xxy_yyzz[i] * fe_0 + tr_z_xxxy_yyzz[i] * pa_x[i];

        tr_z_xxxxy_yzzz[i] = 3.0 * tr_z_xxy_yzzz[i] * fe_0 + tr_z_xxxy_yzzz[i] * pa_x[i];

        tr_z_xxxxy_zzzz[i] = tr_z_xxxx_zzzz[i] * pa_y[i];
    }

    // Set up 660-675 components of targeted buffer : HG

    auto tr_z_xxxxz_xxxx = pbuffer.data(idx_dip_hg + 660);

    auto tr_z_xxxxz_xxxy = pbuffer.data(idx_dip_hg + 661);

    auto tr_z_xxxxz_xxxz = pbuffer.data(idx_dip_hg + 662);

    auto tr_z_xxxxz_xxyy = pbuffer.data(idx_dip_hg + 663);

    auto tr_z_xxxxz_xxyz = pbuffer.data(idx_dip_hg + 664);

    auto tr_z_xxxxz_xxzz = pbuffer.data(idx_dip_hg + 665);

    auto tr_z_xxxxz_xyyy = pbuffer.data(idx_dip_hg + 666);

    auto tr_z_xxxxz_xyyz = pbuffer.data(idx_dip_hg + 667);

    auto tr_z_xxxxz_xyzz = pbuffer.data(idx_dip_hg + 668);

    auto tr_z_xxxxz_xzzz = pbuffer.data(idx_dip_hg + 669);

    auto tr_z_xxxxz_yyyy = pbuffer.data(idx_dip_hg + 670);

    auto tr_z_xxxxz_yyyz = pbuffer.data(idx_dip_hg + 671);

    auto tr_z_xxxxz_yyzz = pbuffer.data(idx_dip_hg + 672);

    auto tr_z_xxxxz_yzzz = pbuffer.data(idx_dip_hg + 673);

    auto tr_z_xxxxz_zzzz = pbuffer.data(idx_dip_hg + 674);

#pragma omp simd aligned(pa_x,                \
                             pa_z,            \
                             tr_z_xxxx_xxxx,  \
                             tr_z_xxxx_xxxy,  \
                             tr_z_xxxx_xxyy,  \
                             tr_z_xxxx_xyyy,  \
                             tr_z_xxxxz_xxxx, \
                             tr_z_xxxxz_xxxy, \
                             tr_z_xxxxz_xxxz, \
                             tr_z_xxxxz_xxyy, \
                             tr_z_xxxxz_xxyz, \
                             tr_z_xxxxz_xxzz, \
                             tr_z_xxxxz_xyyy, \
                             tr_z_xxxxz_xyyz, \
                             tr_z_xxxxz_xyzz, \
                             tr_z_xxxxz_xzzz, \
                             tr_z_xxxxz_yyyy, \
                             tr_z_xxxxz_yyyz, \
                             tr_z_xxxxz_yyzz, \
                             tr_z_xxxxz_yzzz, \
                             tr_z_xxxxz_zzzz, \
                             tr_z_xxxz_xxxz,  \
                             tr_z_xxxz_xxyz,  \
                             tr_z_xxxz_xxz,   \
                             tr_z_xxxz_xxzz,  \
                             tr_z_xxxz_xyyz,  \
                             tr_z_xxxz_xyz,   \
                             tr_z_xxxz_xyzz,  \
                             tr_z_xxxz_xzz,   \
                             tr_z_xxxz_xzzz,  \
                             tr_z_xxxz_yyyy,  \
                             tr_z_xxxz_yyyz,  \
                             tr_z_xxxz_yyz,   \
                             tr_z_xxxz_yyzz,  \
                             tr_z_xxxz_yzz,   \
                             tr_z_xxxz_yzzz,  \
                             tr_z_xxxz_zzz,   \
                             tr_z_xxxz_zzzz,  \
                             tr_z_xxz_xxxz,   \
                             tr_z_xxz_xxyz,   \
                             tr_z_xxz_xxzz,   \
                             tr_z_xxz_xyyz,   \
                             tr_z_xxz_xyzz,   \
                             tr_z_xxz_xzzz,   \
                             tr_z_xxz_yyyy,   \
                             tr_z_xxz_yyyz,   \
                             tr_z_xxz_yyzz,   \
                             tr_z_xxz_yzzz,   \
                             tr_z_xxz_zzzz,   \
                             ts_xxxx_xxxx,    \
                             ts_xxxx_xxxy,    \
                             ts_xxxx_xxyy,    \
                             ts_xxxx_xyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxz_xxxx[i] = ts_xxxx_xxxx[i] * fe_0 + tr_z_xxxx_xxxx[i] * pa_z[i];

        tr_z_xxxxz_xxxy[i] = ts_xxxx_xxxy[i] * fe_0 + tr_z_xxxx_xxxy[i] * pa_z[i];

        tr_z_xxxxz_xxxz[i] = 3.0 * tr_z_xxz_xxxz[i] * fe_0 + 3.0 * tr_z_xxxz_xxz[i] * fe_0 + tr_z_xxxz_xxxz[i] * pa_x[i];

        tr_z_xxxxz_xxyy[i] = ts_xxxx_xxyy[i] * fe_0 + tr_z_xxxx_xxyy[i] * pa_z[i];

        tr_z_xxxxz_xxyz[i] = 3.0 * tr_z_xxz_xxyz[i] * fe_0 + 2.0 * tr_z_xxxz_xyz[i] * fe_0 + tr_z_xxxz_xxyz[i] * pa_x[i];

        tr_z_xxxxz_xxzz[i] = 3.0 * tr_z_xxz_xxzz[i] * fe_0 + 2.0 * tr_z_xxxz_xzz[i] * fe_0 + tr_z_xxxz_xxzz[i] * pa_x[i];

        tr_z_xxxxz_xyyy[i] = ts_xxxx_xyyy[i] * fe_0 + tr_z_xxxx_xyyy[i] * pa_z[i];

        tr_z_xxxxz_xyyz[i] = 3.0 * tr_z_xxz_xyyz[i] * fe_0 + tr_z_xxxz_yyz[i] * fe_0 + tr_z_xxxz_xyyz[i] * pa_x[i];

        tr_z_xxxxz_xyzz[i] = 3.0 * tr_z_xxz_xyzz[i] * fe_0 + tr_z_xxxz_yzz[i] * fe_0 + tr_z_xxxz_xyzz[i] * pa_x[i];

        tr_z_xxxxz_xzzz[i] = 3.0 * tr_z_xxz_xzzz[i] * fe_0 + tr_z_xxxz_zzz[i] * fe_0 + tr_z_xxxz_xzzz[i] * pa_x[i];

        tr_z_xxxxz_yyyy[i] = 3.0 * tr_z_xxz_yyyy[i] * fe_0 + tr_z_xxxz_yyyy[i] * pa_x[i];

        tr_z_xxxxz_yyyz[i] = 3.0 * tr_z_xxz_yyyz[i] * fe_0 + tr_z_xxxz_yyyz[i] * pa_x[i];

        tr_z_xxxxz_yyzz[i] = 3.0 * tr_z_xxz_yyzz[i] * fe_0 + tr_z_xxxz_yyzz[i] * pa_x[i];

        tr_z_xxxxz_yzzz[i] = 3.0 * tr_z_xxz_yzzz[i] * fe_0 + tr_z_xxxz_yzzz[i] * pa_x[i];

        tr_z_xxxxz_zzzz[i] = 3.0 * tr_z_xxz_zzzz[i] * fe_0 + tr_z_xxxz_zzzz[i] * pa_x[i];
    }

    // Set up 675-690 components of targeted buffer : HG

    auto tr_z_xxxyy_xxxx = pbuffer.data(idx_dip_hg + 675);

    auto tr_z_xxxyy_xxxy = pbuffer.data(idx_dip_hg + 676);

    auto tr_z_xxxyy_xxxz = pbuffer.data(idx_dip_hg + 677);

    auto tr_z_xxxyy_xxyy = pbuffer.data(idx_dip_hg + 678);

    auto tr_z_xxxyy_xxyz = pbuffer.data(idx_dip_hg + 679);

    auto tr_z_xxxyy_xxzz = pbuffer.data(idx_dip_hg + 680);

    auto tr_z_xxxyy_xyyy = pbuffer.data(idx_dip_hg + 681);

    auto tr_z_xxxyy_xyyz = pbuffer.data(idx_dip_hg + 682);

    auto tr_z_xxxyy_xyzz = pbuffer.data(idx_dip_hg + 683);

    auto tr_z_xxxyy_xzzz = pbuffer.data(idx_dip_hg + 684);

    auto tr_z_xxxyy_yyyy = pbuffer.data(idx_dip_hg + 685);

    auto tr_z_xxxyy_yyyz = pbuffer.data(idx_dip_hg + 686);

    auto tr_z_xxxyy_yyzz = pbuffer.data(idx_dip_hg + 687);

    auto tr_z_xxxyy_yzzz = pbuffer.data(idx_dip_hg + 688);

    auto tr_z_xxxyy_zzzz = pbuffer.data(idx_dip_hg + 689);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tr_z_xxx_xxxx,   \
                             tr_z_xxx_xxxz,   \
                             tr_z_xxx_xxzz,   \
                             tr_z_xxx_xzzz,   \
                             tr_z_xxxy_xxxx,  \
                             tr_z_xxxy_xxxz,  \
                             tr_z_xxxy_xxzz,  \
                             tr_z_xxxy_xzzz,  \
                             tr_z_xxxyy_xxxx, \
                             tr_z_xxxyy_xxxy, \
                             tr_z_xxxyy_xxxz, \
                             tr_z_xxxyy_xxyy, \
                             tr_z_xxxyy_xxyz, \
                             tr_z_xxxyy_xxzz, \
                             tr_z_xxxyy_xyyy, \
                             tr_z_xxxyy_xyyz, \
                             tr_z_xxxyy_xyzz, \
                             tr_z_xxxyy_xzzz, \
                             tr_z_xxxyy_yyyy, \
                             tr_z_xxxyy_yyyz, \
                             tr_z_xxxyy_yyzz, \
                             tr_z_xxxyy_yzzz, \
                             tr_z_xxxyy_zzzz, \
                             tr_z_xxyy_xxxy,  \
                             tr_z_xxyy_xxy,   \
                             tr_z_xxyy_xxyy,  \
                             tr_z_xxyy_xxyz,  \
                             tr_z_xxyy_xyy,   \
                             tr_z_xxyy_xyyy,  \
                             tr_z_xxyy_xyyz,  \
                             tr_z_xxyy_xyz,   \
                             tr_z_xxyy_xyzz,  \
                             tr_z_xxyy_yyy,   \
                             tr_z_xxyy_yyyy,  \
                             tr_z_xxyy_yyyz,  \
                             tr_z_xxyy_yyz,   \
                             tr_z_xxyy_yyzz,  \
                             tr_z_xxyy_yzz,   \
                             tr_z_xxyy_yzzz,  \
                             tr_z_xxyy_zzzz,  \
                             tr_z_xyy_xxxy,   \
                             tr_z_xyy_xxyy,   \
                             tr_z_xyy_xxyz,   \
                             tr_z_xyy_xyyy,   \
                             tr_z_xyy_xyyz,   \
                             tr_z_xyy_xyzz,   \
                             tr_z_xyy_yyyy,   \
                             tr_z_xyy_yyyz,   \
                             tr_z_xyy_yyzz,   \
                             tr_z_xyy_yzzz,   \
                             tr_z_xyy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyy_xxxx[i] = tr_z_xxx_xxxx[i] * fe_0 + tr_z_xxxy_xxxx[i] * pa_y[i];

        tr_z_xxxyy_xxxy[i] = 2.0 * tr_z_xyy_xxxy[i] * fe_0 + 3.0 * tr_z_xxyy_xxy[i] * fe_0 + tr_z_xxyy_xxxy[i] * pa_x[i];

        tr_z_xxxyy_xxxz[i] = tr_z_xxx_xxxz[i] * fe_0 + tr_z_xxxy_xxxz[i] * pa_y[i];

        tr_z_xxxyy_xxyy[i] = 2.0 * tr_z_xyy_xxyy[i] * fe_0 + 2.0 * tr_z_xxyy_xyy[i] * fe_0 + tr_z_xxyy_xxyy[i] * pa_x[i];

        tr_z_xxxyy_xxyz[i] = 2.0 * tr_z_xyy_xxyz[i] * fe_0 + 2.0 * tr_z_xxyy_xyz[i] * fe_0 + tr_z_xxyy_xxyz[i] * pa_x[i];

        tr_z_xxxyy_xxzz[i] = tr_z_xxx_xxzz[i] * fe_0 + tr_z_xxxy_xxzz[i] * pa_y[i];

        tr_z_xxxyy_xyyy[i] = 2.0 * tr_z_xyy_xyyy[i] * fe_0 + tr_z_xxyy_yyy[i] * fe_0 + tr_z_xxyy_xyyy[i] * pa_x[i];

        tr_z_xxxyy_xyyz[i] = 2.0 * tr_z_xyy_xyyz[i] * fe_0 + tr_z_xxyy_yyz[i] * fe_0 + tr_z_xxyy_xyyz[i] * pa_x[i];

        tr_z_xxxyy_xyzz[i] = 2.0 * tr_z_xyy_xyzz[i] * fe_0 + tr_z_xxyy_yzz[i] * fe_0 + tr_z_xxyy_xyzz[i] * pa_x[i];

        tr_z_xxxyy_xzzz[i] = tr_z_xxx_xzzz[i] * fe_0 + tr_z_xxxy_xzzz[i] * pa_y[i];

        tr_z_xxxyy_yyyy[i] = 2.0 * tr_z_xyy_yyyy[i] * fe_0 + tr_z_xxyy_yyyy[i] * pa_x[i];

        tr_z_xxxyy_yyyz[i] = 2.0 * tr_z_xyy_yyyz[i] * fe_0 + tr_z_xxyy_yyyz[i] * pa_x[i];

        tr_z_xxxyy_yyzz[i] = 2.0 * tr_z_xyy_yyzz[i] * fe_0 + tr_z_xxyy_yyzz[i] * pa_x[i];

        tr_z_xxxyy_yzzz[i] = 2.0 * tr_z_xyy_yzzz[i] * fe_0 + tr_z_xxyy_yzzz[i] * pa_x[i];

        tr_z_xxxyy_zzzz[i] = 2.0 * tr_z_xyy_zzzz[i] * fe_0 + tr_z_xxyy_zzzz[i] * pa_x[i];
    }

    // Set up 690-705 components of targeted buffer : HG

    auto tr_z_xxxyz_xxxx = pbuffer.data(idx_dip_hg + 690);

    auto tr_z_xxxyz_xxxy = pbuffer.data(idx_dip_hg + 691);

    auto tr_z_xxxyz_xxxz = pbuffer.data(idx_dip_hg + 692);

    auto tr_z_xxxyz_xxyy = pbuffer.data(idx_dip_hg + 693);

    auto tr_z_xxxyz_xxyz = pbuffer.data(idx_dip_hg + 694);

    auto tr_z_xxxyz_xxzz = pbuffer.data(idx_dip_hg + 695);

    auto tr_z_xxxyz_xyyy = pbuffer.data(idx_dip_hg + 696);

    auto tr_z_xxxyz_xyyz = pbuffer.data(idx_dip_hg + 697);

    auto tr_z_xxxyz_xyzz = pbuffer.data(idx_dip_hg + 698);

    auto tr_z_xxxyz_xzzz = pbuffer.data(idx_dip_hg + 699);

    auto tr_z_xxxyz_yyyy = pbuffer.data(idx_dip_hg + 700);

    auto tr_z_xxxyz_yyyz = pbuffer.data(idx_dip_hg + 701);

    auto tr_z_xxxyz_yyzz = pbuffer.data(idx_dip_hg + 702);

    auto tr_z_xxxyz_yzzz = pbuffer.data(idx_dip_hg + 703);

    auto tr_z_xxxyz_zzzz = pbuffer.data(idx_dip_hg + 704);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tr_z_xxxyz_xxxx, \
                             tr_z_xxxyz_xxxy, \
                             tr_z_xxxyz_xxxz, \
                             tr_z_xxxyz_xxyy, \
                             tr_z_xxxyz_xxyz, \
                             tr_z_xxxyz_xxzz, \
                             tr_z_xxxyz_xyyy, \
                             tr_z_xxxyz_xyyz, \
                             tr_z_xxxyz_xyzz, \
                             tr_z_xxxyz_xzzz, \
                             tr_z_xxxyz_yyyy, \
                             tr_z_xxxyz_yyyz, \
                             tr_z_xxxyz_yyzz, \
                             tr_z_xxxyz_yzzz, \
                             tr_z_xxxyz_zzzz, \
                             tr_z_xxxz_xxx,   \
                             tr_z_xxxz_xxxx,  \
                             tr_z_xxxz_xxxy,  \
                             tr_z_xxxz_xxxz,  \
                             tr_z_xxxz_xxy,   \
                             tr_z_xxxz_xxyy,  \
                             tr_z_xxxz_xxyz,  \
                             tr_z_xxxz_xxz,   \
                             tr_z_xxxz_xxzz,  \
                             tr_z_xxxz_xyy,   \
                             tr_z_xxxz_xyyy,  \
                             tr_z_xxxz_xyyz,  \
                             tr_z_xxxz_xyz,   \
                             tr_z_xxxz_xyzz,  \
                             tr_z_xxxz_xzz,   \
                             tr_z_xxxz_xzzz,  \
                             tr_z_xxxz_zzzz,  \
                             tr_z_xxyz_yyyy,  \
                             tr_z_xxyz_yyyz,  \
                             tr_z_xxyz_yyzz,  \
                             tr_z_xxyz_yzzz,  \
                             tr_z_xyz_yyyy,   \
                             tr_z_xyz_yyyz,   \
                             tr_z_xyz_yyzz,   \
                             tr_z_xyz_yzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyz_xxxx[i] = tr_z_xxxz_xxxx[i] * pa_y[i];

        tr_z_xxxyz_xxxy[i] = tr_z_xxxz_xxx[i] * fe_0 + tr_z_xxxz_xxxy[i] * pa_y[i];

        tr_z_xxxyz_xxxz[i] = tr_z_xxxz_xxxz[i] * pa_y[i];

        tr_z_xxxyz_xxyy[i] = 2.0 * tr_z_xxxz_xxy[i] * fe_0 + tr_z_xxxz_xxyy[i] * pa_y[i];

        tr_z_xxxyz_xxyz[i] = tr_z_xxxz_xxz[i] * fe_0 + tr_z_xxxz_xxyz[i] * pa_y[i];

        tr_z_xxxyz_xxzz[i] = tr_z_xxxz_xxzz[i] * pa_y[i];

        tr_z_xxxyz_xyyy[i] = 3.0 * tr_z_xxxz_xyy[i] * fe_0 + tr_z_xxxz_xyyy[i] * pa_y[i];

        tr_z_xxxyz_xyyz[i] = 2.0 * tr_z_xxxz_xyz[i] * fe_0 + tr_z_xxxz_xyyz[i] * pa_y[i];

        tr_z_xxxyz_xyzz[i] = tr_z_xxxz_xzz[i] * fe_0 + tr_z_xxxz_xyzz[i] * pa_y[i];

        tr_z_xxxyz_xzzz[i] = tr_z_xxxz_xzzz[i] * pa_y[i];

        tr_z_xxxyz_yyyy[i] = 2.0 * tr_z_xyz_yyyy[i] * fe_0 + tr_z_xxyz_yyyy[i] * pa_x[i];

        tr_z_xxxyz_yyyz[i] = 2.0 * tr_z_xyz_yyyz[i] * fe_0 + tr_z_xxyz_yyyz[i] * pa_x[i];

        tr_z_xxxyz_yyzz[i] = 2.0 * tr_z_xyz_yyzz[i] * fe_0 + tr_z_xxyz_yyzz[i] * pa_x[i];

        tr_z_xxxyz_yzzz[i] = 2.0 * tr_z_xyz_yzzz[i] * fe_0 + tr_z_xxyz_yzzz[i] * pa_x[i];

        tr_z_xxxyz_zzzz[i] = tr_z_xxxz_zzzz[i] * pa_y[i];
    }

    // Set up 705-720 components of targeted buffer : HG

    auto tr_z_xxxzz_xxxx = pbuffer.data(idx_dip_hg + 705);

    auto tr_z_xxxzz_xxxy = pbuffer.data(idx_dip_hg + 706);

    auto tr_z_xxxzz_xxxz = pbuffer.data(idx_dip_hg + 707);

    auto tr_z_xxxzz_xxyy = pbuffer.data(idx_dip_hg + 708);

    auto tr_z_xxxzz_xxyz = pbuffer.data(idx_dip_hg + 709);

    auto tr_z_xxxzz_xxzz = pbuffer.data(idx_dip_hg + 710);

    auto tr_z_xxxzz_xyyy = pbuffer.data(idx_dip_hg + 711);

    auto tr_z_xxxzz_xyyz = pbuffer.data(idx_dip_hg + 712);

    auto tr_z_xxxzz_xyzz = pbuffer.data(idx_dip_hg + 713);

    auto tr_z_xxxzz_xzzz = pbuffer.data(idx_dip_hg + 714);

    auto tr_z_xxxzz_yyyy = pbuffer.data(idx_dip_hg + 715);

    auto tr_z_xxxzz_yyyz = pbuffer.data(idx_dip_hg + 716);

    auto tr_z_xxxzz_yyzz = pbuffer.data(idx_dip_hg + 717);

    auto tr_z_xxxzz_yzzz = pbuffer.data(idx_dip_hg + 718);

    auto tr_z_xxxzz_zzzz = pbuffer.data(idx_dip_hg + 719);

#pragma omp simd aligned(pa_x,                \
                             tr_z_xxxzz_xxxx, \
                             tr_z_xxxzz_xxxy, \
                             tr_z_xxxzz_xxxz, \
                             tr_z_xxxzz_xxyy, \
                             tr_z_xxxzz_xxyz, \
                             tr_z_xxxzz_xxzz, \
                             tr_z_xxxzz_xyyy, \
                             tr_z_xxxzz_xyyz, \
                             tr_z_xxxzz_xyzz, \
                             tr_z_xxxzz_xzzz, \
                             tr_z_xxxzz_yyyy, \
                             tr_z_xxxzz_yyyz, \
                             tr_z_xxxzz_yyzz, \
                             tr_z_xxxzz_yzzz, \
                             tr_z_xxxzz_zzzz, \
                             tr_z_xxzz_xxx,   \
                             tr_z_xxzz_xxxx,  \
                             tr_z_xxzz_xxxy,  \
                             tr_z_xxzz_xxxz,  \
                             tr_z_xxzz_xxy,   \
                             tr_z_xxzz_xxyy,  \
                             tr_z_xxzz_xxyz,  \
                             tr_z_xxzz_xxz,   \
                             tr_z_xxzz_xxzz,  \
                             tr_z_xxzz_xyy,   \
                             tr_z_xxzz_xyyy,  \
                             tr_z_xxzz_xyyz,  \
                             tr_z_xxzz_xyz,   \
                             tr_z_xxzz_xyzz,  \
                             tr_z_xxzz_xzz,   \
                             tr_z_xxzz_xzzz,  \
                             tr_z_xxzz_yyy,   \
                             tr_z_xxzz_yyyy,  \
                             tr_z_xxzz_yyyz,  \
                             tr_z_xxzz_yyz,   \
                             tr_z_xxzz_yyzz,  \
                             tr_z_xxzz_yzz,   \
                             tr_z_xxzz_yzzz,  \
                             tr_z_xxzz_zzz,   \
                             tr_z_xxzz_zzzz,  \
                             tr_z_xzz_xxxx,   \
                             tr_z_xzz_xxxy,   \
                             tr_z_xzz_xxxz,   \
                             tr_z_xzz_xxyy,   \
                             tr_z_xzz_xxyz,   \
                             tr_z_xzz_xxzz,   \
                             tr_z_xzz_xyyy,   \
                             tr_z_xzz_xyyz,   \
                             tr_z_xzz_xyzz,   \
                             tr_z_xzz_xzzz,   \
                             tr_z_xzz_yyyy,   \
                             tr_z_xzz_yyyz,   \
                             tr_z_xzz_yyzz,   \
                             tr_z_xzz_yzzz,   \
                             tr_z_xzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxzz_xxxx[i] = 2.0 * tr_z_xzz_xxxx[i] * fe_0 + 4.0 * tr_z_xxzz_xxx[i] * fe_0 + tr_z_xxzz_xxxx[i] * pa_x[i];

        tr_z_xxxzz_xxxy[i] = 2.0 * tr_z_xzz_xxxy[i] * fe_0 + 3.0 * tr_z_xxzz_xxy[i] * fe_0 + tr_z_xxzz_xxxy[i] * pa_x[i];

        tr_z_xxxzz_xxxz[i] = 2.0 * tr_z_xzz_xxxz[i] * fe_0 + 3.0 * tr_z_xxzz_xxz[i] * fe_0 + tr_z_xxzz_xxxz[i] * pa_x[i];

        tr_z_xxxzz_xxyy[i] = 2.0 * tr_z_xzz_xxyy[i] * fe_0 + 2.0 * tr_z_xxzz_xyy[i] * fe_0 + tr_z_xxzz_xxyy[i] * pa_x[i];

        tr_z_xxxzz_xxyz[i] = 2.0 * tr_z_xzz_xxyz[i] * fe_0 + 2.0 * tr_z_xxzz_xyz[i] * fe_0 + tr_z_xxzz_xxyz[i] * pa_x[i];

        tr_z_xxxzz_xxzz[i] = 2.0 * tr_z_xzz_xxzz[i] * fe_0 + 2.0 * tr_z_xxzz_xzz[i] * fe_0 + tr_z_xxzz_xxzz[i] * pa_x[i];

        tr_z_xxxzz_xyyy[i] = 2.0 * tr_z_xzz_xyyy[i] * fe_0 + tr_z_xxzz_yyy[i] * fe_0 + tr_z_xxzz_xyyy[i] * pa_x[i];

        tr_z_xxxzz_xyyz[i] = 2.0 * tr_z_xzz_xyyz[i] * fe_0 + tr_z_xxzz_yyz[i] * fe_0 + tr_z_xxzz_xyyz[i] * pa_x[i];

        tr_z_xxxzz_xyzz[i] = 2.0 * tr_z_xzz_xyzz[i] * fe_0 + tr_z_xxzz_yzz[i] * fe_0 + tr_z_xxzz_xyzz[i] * pa_x[i];

        tr_z_xxxzz_xzzz[i] = 2.0 * tr_z_xzz_xzzz[i] * fe_0 + tr_z_xxzz_zzz[i] * fe_0 + tr_z_xxzz_xzzz[i] * pa_x[i];

        tr_z_xxxzz_yyyy[i] = 2.0 * tr_z_xzz_yyyy[i] * fe_0 + tr_z_xxzz_yyyy[i] * pa_x[i];

        tr_z_xxxzz_yyyz[i] = 2.0 * tr_z_xzz_yyyz[i] * fe_0 + tr_z_xxzz_yyyz[i] * pa_x[i];

        tr_z_xxxzz_yyzz[i] = 2.0 * tr_z_xzz_yyzz[i] * fe_0 + tr_z_xxzz_yyzz[i] * pa_x[i];

        tr_z_xxxzz_yzzz[i] = 2.0 * tr_z_xzz_yzzz[i] * fe_0 + tr_z_xxzz_yzzz[i] * pa_x[i];

        tr_z_xxxzz_zzzz[i] = 2.0 * tr_z_xzz_zzzz[i] * fe_0 + tr_z_xxzz_zzzz[i] * pa_x[i];
    }

    // Set up 720-735 components of targeted buffer : HG

    auto tr_z_xxyyy_xxxx = pbuffer.data(idx_dip_hg + 720);

    auto tr_z_xxyyy_xxxy = pbuffer.data(idx_dip_hg + 721);

    auto tr_z_xxyyy_xxxz = pbuffer.data(idx_dip_hg + 722);

    auto tr_z_xxyyy_xxyy = pbuffer.data(idx_dip_hg + 723);

    auto tr_z_xxyyy_xxyz = pbuffer.data(idx_dip_hg + 724);

    auto tr_z_xxyyy_xxzz = pbuffer.data(idx_dip_hg + 725);

    auto tr_z_xxyyy_xyyy = pbuffer.data(idx_dip_hg + 726);

    auto tr_z_xxyyy_xyyz = pbuffer.data(idx_dip_hg + 727);

    auto tr_z_xxyyy_xyzz = pbuffer.data(idx_dip_hg + 728);

    auto tr_z_xxyyy_xzzz = pbuffer.data(idx_dip_hg + 729);

    auto tr_z_xxyyy_yyyy = pbuffer.data(idx_dip_hg + 730);

    auto tr_z_xxyyy_yyyz = pbuffer.data(idx_dip_hg + 731);

    auto tr_z_xxyyy_yyzz = pbuffer.data(idx_dip_hg + 732);

    auto tr_z_xxyyy_yzzz = pbuffer.data(idx_dip_hg + 733);

    auto tr_z_xxyyy_zzzz = pbuffer.data(idx_dip_hg + 734);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tr_z_xxy_xxxx,   \
                             tr_z_xxy_xxxz,   \
                             tr_z_xxy_xxzz,   \
                             tr_z_xxy_xzzz,   \
                             tr_z_xxyy_xxxx,  \
                             tr_z_xxyy_xxxz,  \
                             tr_z_xxyy_xxzz,  \
                             tr_z_xxyy_xzzz,  \
                             tr_z_xxyyy_xxxx, \
                             tr_z_xxyyy_xxxy, \
                             tr_z_xxyyy_xxxz, \
                             tr_z_xxyyy_xxyy, \
                             tr_z_xxyyy_xxyz, \
                             tr_z_xxyyy_xxzz, \
                             tr_z_xxyyy_xyyy, \
                             tr_z_xxyyy_xyyz, \
                             tr_z_xxyyy_xyzz, \
                             tr_z_xxyyy_xzzz, \
                             tr_z_xxyyy_yyyy, \
                             tr_z_xxyyy_yyyz, \
                             tr_z_xxyyy_yyzz, \
                             tr_z_xxyyy_yzzz, \
                             tr_z_xxyyy_zzzz, \
                             tr_z_xyyy_xxxy,  \
                             tr_z_xyyy_xxy,   \
                             tr_z_xyyy_xxyy,  \
                             tr_z_xyyy_xxyz,  \
                             tr_z_xyyy_xyy,   \
                             tr_z_xyyy_xyyy,  \
                             tr_z_xyyy_xyyz,  \
                             tr_z_xyyy_xyz,   \
                             tr_z_xyyy_xyzz,  \
                             tr_z_xyyy_yyy,   \
                             tr_z_xyyy_yyyy,  \
                             tr_z_xyyy_yyyz,  \
                             tr_z_xyyy_yyz,   \
                             tr_z_xyyy_yyzz,  \
                             tr_z_xyyy_yzz,   \
                             tr_z_xyyy_yzzz,  \
                             tr_z_xyyy_zzzz,  \
                             tr_z_yyy_xxxy,   \
                             tr_z_yyy_xxyy,   \
                             tr_z_yyy_xxyz,   \
                             tr_z_yyy_xyyy,   \
                             tr_z_yyy_xyyz,   \
                             tr_z_yyy_xyzz,   \
                             tr_z_yyy_yyyy,   \
                             tr_z_yyy_yyyz,   \
                             tr_z_yyy_yyzz,   \
                             tr_z_yyy_yzzz,   \
                             tr_z_yyy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyy_xxxx[i] = 2.0 * tr_z_xxy_xxxx[i] * fe_0 + tr_z_xxyy_xxxx[i] * pa_y[i];

        tr_z_xxyyy_xxxy[i] = tr_z_yyy_xxxy[i] * fe_0 + 3.0 * tr_z_xyyy_xxy[i] * fe_0 + tr_z_xyyy_xxxy[i] * pa_x[i];

        tr_z_xxyyy_xxxz[i] = 2.0 * tr_z_xxy_xxxz[i] * fe_0 + tr_z_xxyy_xxxz[i] * pa_y[i];

        tr_z_xxyyy_xxyy[i] = tr_z_yyy_xxyy[i] * fe_0 + 2.0 * tr_z_xyyy_xyy[i] * fe_0 + tr_z_xyyy_xxyy[i] * pa_x[i];

        tr_z_xxyyy_xxyz[i] = tr_z_yyy_xxyz[i] * fe_0 + 2.0 * tr_z_xyyy_xyz[i] * fe_0 + tr_z_xyyy_xxyz[i] * pa_x[i];

        tr_z_xxyyy_xxzz[i] = 2.0 * tr_z_xxy_xxzz[i] * fe_0 + tr_z_xxyy_xxzz[i] * pa_y[i];

        tr_z_xxyyy_xyyy[i] = tr_z_yyy_xyyy[i] * fe_0 + tr_z_xyyy_yyy[i] * fe_0 + tr_z_xyyy_xyyy[i] * pa_x[i];

        tr_z_xxyyy_xyyz[i] = tr_z_yyy_xyyz[i] * fe_0 + tr_z_xyyy_yyz[i] * fe_0 + tr_z_xyyy_xyyz[i] * pa_x[i];

        tr_z_xxyyy_xyzz[i] = tr_z_yyy_xyzz[i] * fe_0 + tr_z_xyyy_yzz[i] * fe_0 + tr_z_xyyy_xyzz[i] * pa_x[i];

        tr_z_xxyyy_xzzz[i] = 2.0 * tr_z_xxy_xzzz[i] * fe_0 + tr_z_xxyy_xzzz[i] * pa_y[i];

        tr_z_xxyyy_yyyy[i] = tr_z_yyy_yyyy[i] * fe_0 + tr_z_xyyy_yyyy[i] * pa_x[i];

        tr_z_xxyyy_yyyz[i] = tr_z_yyy_yyyz[i] * fe_0 + tr_z_xyyy_yyyz[i] * pa_x[i];

        tr_z_xxyyy_yyzz[i] = tr_z_yyy_yyzz[i] * fe_0 + tr_z_xyyy_yyzz[i] * pa_x[i];

        tr_z_xxyyy_yzzz[i] = tr_z_yyy_yzzz[i] * fe_0 + tr_z_xyyy_yzzz[i] * pa_x[i];

        tr_z_xxyyy_zzzz[i] = tr_z_yyy_zzzz[i] * fe_0 + tr_z_xyyy_zzzz[i] * pa_x[i];
    }

    // Set up 735-750 components of targeted buffer : HG

    auto tr_z_xxyyz_xxxx = pbuffer.data(idx_dip_hg + 735);

    auto tr_z_xxyyz_xxxy = pbuffer.data(idx_dip_hg + 736);

    auto tr_z_xxyyz_xxxz = pbuffer.data(idx_dip_hg + 737);

    auto tr_z_xxyyz_xxyy = pbuffer.data(idx_dip_hg + 738);

    auto tr_z_xxyyz_xxyz = pbuffer.data(idx_dip_hg + 739);

    auto tr_z_xxyyz_xxzz = pbuffer.data(idx_dip_hg + 740);

    auto tr_z_xxyyz_xyyy = pbuffer.data(idx_dip_hg + 741);

    auto tr_z_xxyyz_xyyz = pbuffer.data(idx_dip_hg + 742);

    auto tr_z_xxyyz_xyzz = pbuffer.data(idx_dip_hg + 743);

    auto tr_z_xxyyz_xzzz = pbuffer.data(idx_dip_hg + 744);

    auto tr_z_xxyyz_yyyy = pbuffer.data(idx_dip_hg + 745);

    auto tr_z_xxyyz_yyyz = pbuffer.data(idx_dip_hg + 746);

    auto tr_z_xxyyz_yyzz = pbuffer.data(idx_dip_hg + 747);

    auto tr_z_xxyyz_yzzz = pbuffer.data(idx_dip_hg + 748);

    auto tr_z_xxyyz_zzzz = pbuffer.data(idx_dip_hg + 749);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             pa_z,            \
                             tr_z_xxyy_xxxy,  \
                             tr_z_xxyy_xxyy,  \
                             tr_z_xxyy_xyyy,  \
                             tr_z_xxyyz_xxxx, \
                             tr_z_xxyyz_xxxy, \
                             tr_z_xxyyz_xxxz, \
                             tr_z_xxyyz_xxyy, \
                             tr_z_xxyyz_xxyz, \
                             tr_z_xxyyz_xxzz, \
                             tr_z_xxyyz_xyyy, \
                             tr_z_xxyyz_xyyz, \
                             tr_z_xxyyz_xyzz, \
                             tr_z_xxyyz_xzzz, \
                             tr_z_xxyyz_yyyy, \
                             tr_z_xxyyz_yyyz, \
                             tr_z_xxyyz_yyzz, \
                             tr_z_xxyyz_yzzz, \
                             tr_z_xxyyz_zzzz, \
                             tr_z_xxyz_xxxx,  \
                             tr_z_xxyz_xxxz,  \
                             tr_z_xxyz_xxzz,  \
                             tr_z_xxyz_xzzz,  \
                             tr_z_xxz_xxxx,   \
                             tr_z_xxz_xxxz,   \
                             tr_z_xxz_xxzz,   \
                             tr_z_xxz_xzzz,   \
                             tr_z_xyyz_xxyz,  \
                             tr_z_xyyz_xyyz,  \
                             tr_z_xyyz_xyz,   \
                             tr_z_xyyz_xyzz,  \
                             tr_z_xyyz_yyyy,  \
                             tr_z_xyyz_yyyz,  \
                             tr_z_xyyz_yyz,   \
                             tr_z_xyyz_yyzz,  \
                             tr_z_xyyz_yzz,   \
                             tr_z_xyyz_yzzz,  \
                             tr_z_xyyz_zzzz,  \
                             tr_z_yyz_xxyz,   \
                             tr_z_yyz_xyyz,   \
                             tr_z_yyz_xyzz,   \
                             tr_z_yyz_yyyy,   \
                             tr_z_yyz_yyyz,   \
                             tr_z_yyz_yyzz,   \
                             tr_z_yyz_yzzz,   \
                             tr_z_yyz_zzzz,   \
                             ts_xxyy_xxxy,    \
                             ts_xxyy_xxyy,    \
                             ts_xxyy_xyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyz_xxxx[i] = tr_z_xxz_xxxx[i] * fe_0 + tr_z_xxyz_xxxx[i] * pa_y[i];

        tr_z_xxyyz_xxxy[i] = ts_xxyy_xxxy[i] * fe_0 + tr_z_xxyy_xxxy[i] * pa_z[i];

        tr_z_xxyyz_xxxz[i] = tr_z_xxz_xxxz[i] * fe_0 + tr_z_xxyz_xxxz[i] * pa_y[i];

        tr_z_xxyyz_xxyy[i] = ts_xxyy_xxyy[i] * fe_0 + tr_z_xxyy_xxyy[i] * pa_z[i];

        tr_z_xxyyz_xxyz[i] = tr_z_yyz_xxyz[i] * fe_0 + 2.0 * tr_z_xyyz_xyz[i] * fe_0 + tr_z_xyyz_xxyz[i] * pa_x[i];

        tr_z_xxyyz_xxzz[i] = tr_z_xxz_xxzz[i] * fe_0 + tr_z_xxyz_xxzz[i] * pa_y[i];

        tr_z_xxyyz_xyyy[i] = ts_xxyy_xyyy[i] * fe_0 + tr_z_xxyy_xyyy[i] * pa_z[i];

        tr_z_xxyyz_xyyz[i] = tr_z_yyz_xyyz[i] * fe_0 + tr_z_xyyz_yyz[i] * fe_0 + tr_z_xyyz_xyyz[i] * pa_x[i];

        tr_z_xxyyz_xyzz[i] = tr_z_yyz_xyzz[i] * fe_0 + tr_z_xyyz_yzz[i] * fe_0 + tr_z_xyyz_xyzz[i] * pa_x[i];

        tr_z_xxyyz_xzzz[i] = tr_z_xxz_xzzz[i] * fe_0 + tr_z_xxyz_xzzz[i] * pa_y[i];

        tr_z_xxyyz_yyyy[i] = tr_z_yyz_yyyy[i] * fe_0 + tr_z_xyyz_yyyy[i] * pa_x[i];

        tr_z_xxyyz_yyyz[i] = tr_z_yyz_yyyz[i] * fe_0 + tr_z_xyyz_yyyz[i] * pa_x[i];

        tr_z_xxyyz_yyzz[i] = tr_z_yyz_yyzz[i] * fe_0 + tr_z_xyyz_yyzz[i] * pa_x[i];

        tr_z_xxyyz_yzzz[i] = tr_z_yyz_yzzz[i] * fe_0 + tr_z_xyyz_yzzz[i] * pa_x[i];

        tr_z_xxyyz_zzzz[i] = tr_z_yyz_zzzz[i] * fe_0 + tr_z_xyyz_zzzz[i] * pa_x[i];
    }

    // Set up 750-765 components of targeted buffer : HG

    auto tr_z_xxyzz_xxxx = pbuffer.data(idx_dip_hg + 750);

    auto tr_z_xxyzz_xxxy = pbuffer.data(idx_dip_hg + 751);

    auto tr_z_xxyzz_xxxz = pbuffer.data(idx_dip_hg + 752);

    auto tr_z_xxyzz_xxyy = pbuffer.data(idx_dip_hg + 753);

    auto tr_z_xxyzz_xxyz = pbuffer.data(idx_dip_hg + 754);

    auto tr_z_xxyzz_xxzz = pbuffer.data(idx_dip_hg + 755);

    auto tr_z_xxyzz_xyyy = pbuffer.data(idx_dip_hg + 756);

    auto tr_z_xxyzz_xyyz = pbuffer.data(idx_dip_hg + 757);

    auto tr_z_xxyzz_xyzz = pbuffer.data(idx_dip_hg + 758);

    auto tr_z_xxyzz_xzzz = pbuffer.data(idx_dip_hg + 759);

    auto tr_z_xxyzz_yyyy = pbuffer.data(idx_dip_hg + 760);

    auto tr_z_xxyzz_yyyz = pbuffer.data(idx_dip_hg + 761);

    auto tr_z_xxyzz_yyzz = pbuffer.data(idx_dip_hg + 762);

    auto tr_z_xxyzz_yzzz = pbuffer.data(idx_dip_hg + 763);

    auto tr_z_xxyzz_zzzz = pbuffer.data(idx_dip_hg + 764);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tr_z_xxyzz_xxxx, \
                             tr_z_xxyzz_xxxy, \
                             tr_z_xxyzz_xxxz, \
                             tr_z_xxyzz_xxyy, \
                             tr_z_xxyzz_xxyz, \
                             tr_z_xxyzz_xxzz, \
                             tr_z_xxyzz_xyyy, \
                             tr_z_xxyzz_xyyz, \
                             tr_z_xxyzz_xyzz, \
                             tr_z_xxyzz_xzzz, \
                             tr_z_xxyzz_yyyy, \
                             tr_z_xxyzz_yyyz, \
                             tr_z_xxyzz_yyzz, \
                             tr_z_xxyzz_yzzz, \
                             tr_z_xxyzz_zzzz, \
                             tr_z_xxzz_xxx,   \
                             tr_z_xxzz_xxxx,  \
                             tr_z_xxzz_xxxy,  \
                             tr_z_xxzz_xxxz,  \
                             tr_z_xxzz_xxy,   \
                             tr_z_xxzz_xxyy,  \
                             tr_z_xxzz_xxyz,  \
                             tr_z_xxzz_xxz,   \
                             tr_z_xxzz_xxzz,  \
                             tr_z_xxzz_xyy,   \
                             tr_z_xxzz_xyyy,  \
                             tr_z_xxzz_xyyz,  \
                             tr_z_xxzz_xyz,   \
                             tr_z_xxzz_xyzz,  \
                             tr_z_xxzz_xzz,   \
                             tr_z_xxzz_xzzz,  \
                             tr_z_xxzz_zzzz,  \
                             tr_z_xyzz_yyyy,  \
                             tr_z_xyzz_yyyz,  \
                             tr_z_xyzz_yyzz,  \
                             tr_z_xyzz_yzzz,  \
                             tr_z_yzz_yyyy,   \
                             tr_z_yzz_yyyz,   \
                             tr_z_yzz_yyzz,   \
                             tr_z_yzz_yzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyzz_xxxx[i] = tr_z_xxzz_xxxx[i] * pa_y[i];

        tr_z_xxyzz_xxxy[i] = tr_z_xxzz_xxx[i] * fe_0 + tr_z_xxzz_xxxy[i] * pa_y[i];

        tr_z_xxyzz_xxxz[i] = tr_z_xxzz_xxxz[i] * pa_y[i];

        tr_z_xxyzz_xxyy[i] = 2.0 * tr_z_xxzz_xxy[i] * fe_0 + tr_z_xxzz_xxyy[i] * pa_y[i];

        tr_z_xxyzz_xxyz[i] = tr_z_xxzz_xxz[i] * fe_0 + tr_z_xxzz_xxyz[i] * pa_y[i];

        tr_z_xxyzz_xxzz[i] = tr_z_xxzz_xxzz[i] * pa_y[i];

        tr_z_xxyzz_xyyy[i] = 3.0 * tr_z_xxzz_xyy[i] * fe_0 + tr_z_xxzz_xyyy[i] * pa_y[i];

        tr_z_xxyzz_xyyz[i] = 2.0 * tr_z_xxzz_xyz[i] * fe_0 + tr_z_xxzz_xyyz[i] * pa_y[i];

        tr_z_xxyzz_xyzz[i] = tr_z_xxzz_xzz[i] * fe_0 + tr_z_xxzz_xyzz[i] * pa_y[i];

        tr_z_xxyzz_xzzz[i] = tr_z_xxzz_xzzz[i] * pa_y[i];

        tr_z_xxyzz_yyyy[i] = tr_z_yzz_yyyy[i] * fe_0 + tr_z_xyzz_yyyy[i] * pa_x[i];

        tr_z_xxyzz_yyyz[i] = tr_z_yzz_yyyz[i] * fe_0 + tr_z_xyzz_yyyz[i] * pa_x[i];

        tr_z_xxyzz_yyzz[i] = tr_z_yzz_yyzz[i] * fe_0 + tr_z_xyzz_yyzz[i] * pa_x[i];

        tr_z_xxyzz_yzzz[i] = tr_z_yzz_yzzz[i] * fe_0 + tr_z_xyzz_yzzz[i] * pa_x[i];

        tr_z_xxyzz_zzzz[i] = tr_z_xxzz_zzzz[i] * pa_y[i];
    }

    // Set up 765-780 components of targeted buffer : HG

    auto tr_z_xxzzz_xxxx = pbuffer.data(idx_dip_hg + 765);

    auto tr_z_xxzzz_xxxy = pbuffer.data(idx_dip_hg + 766);

    auto tr_z_xxzzz_xxxz = pbuffer.data(idx_dip_hg + 767);

    auto tr_z_xxzzz_xxyy = pbuffer.data(idx_dip_hg + 768);

    auto tr_z_xxzzz_xxyz = pbuffer.data(idx_dip_hg + 769);

    auto tr_z_xxzzz_xxzz = pbuffer.data(idx_dip_hg + 770);

    auto tr_z_xxzzz_xyyy = pbuffer.data(idx_dip_hg + 771);

    auto tr_z_xxzzz_xyyz = pbuffer.data(idx_dip_hg + 772);

    auto tr_z_xxzzz_xyzz = pbuffer.data(idx_dip_hg + 773);

    auto tr_z_xxzzz_xzzz = pbuffer.data(idx_dip_hg + 774);

    auto tr_z_xxzzz_yyyy = pbuffer.data(idx_dip_hg + 775);

    auto tr_z_xxzzz_yyyz = pbuffer.data(idx_dip_hg + 776);

    auto tr_z_xxzzz_yyzz = pbuffer.data(idx_dip_hg + 777);

    auto tr_z_xxzzz_yzzz = pbuffer.data(idx_dip_hg + 778);

    auto tr_z_xxzzz_zzzz = pbuffer.data(idx_dip_hg + 779);

#pragma omp simd aligned(pa_x,                \
                             tr_z_xxzzz_xxxx, \
                             tr_z_xxzzz_xxxy, \
                             tr_z_xxzzz_xxxz, \
                             tr_z_xxzzz_xxyy, \
                             tr_z_xxzzz_xxyz, \
                             tr_z_xxzzz_xxzz, \
                             tr_z_xxzzz_xyyy, \
                             tr_z_xxzzz_xyyz, \
                             tr_z_xxzzz_xyzz, \
                             tr_z_xxzzz_xzzz, \
                             tr_z_xxzzz_yyyy, \
                             tr_z_xxzzz_yyyz, \
                             tr_z_xxzzz_yyzz, \
                             tr_z_xxzzz_yzzz, \
                             tr_z_xxzzz_zzzz, \
                             tr_z_xzzz_xxx,   \
                             tr_z_xzzz_xxxx,  \
                             tr_z_xzzz_xxxy,  \
                             tr_z_xzzz_xxxz,  \
                             tr_z_xzzz_xxy,   \
                             tr_z_xzzz_xxyy,  \
                             tr_z_xzzz_xxyz,  \
                             tr_z_xzzz_xxz,   \
                             tr_z_xzzz_xxzz,  \
                             tr_z_xzzz_xyy,   \
                             tr_z_xzzz_xyyy,  \
                             tr_z_xzzz_xyyz,  \
                             tr_z_xzzz_xyz,   \
                             tr_z_xzzz_xyzz,  \
                             tr_z_xzzz_xzz,   \
                             tr_z_xzzz_xzzz,  \
                             tr_z_xzzz_yyy,   \
                             tr_z_xzzz_yyyy,  \
                             tr_z_xzzz_yyyz,  \
                             tr_z_xzzz_yyz,   \
                             tr_z_xzzz_yyzz,  \
                             tr_z_xzzz_yzz,   \
                             tr_z_xzzz_yzzz,  \
                             tr_z_xzzz_zzz,   \
                             tr_z_xzzz_zzzz,  \
                             tr_z_zzz_xxxx,   \
                             tr_z_zzz_xxxy,   \
                             tr_z_zzz_xxxz,   \
                             tr_z_zzz_xxyy,   \
                             tr_z_zzz_xxyz,   \
                             tr_z_zzz_xxzz,   \
                             tr_z_zzz_xyyy,   \
                             tr_z_zzz_xyyz,   \
                             tr_z_zzz_xyzz,   \
                             tr_z_zzz_xzzz,   \
                             tr_z_zzz_yyyy,   \
                             tr_z_zzz_yyyz,   \
                             tr_z_zzz_yyzz,   \
                             tr_z_zzz_yzzz,   \
                             tr_z_zzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzzz_xxxx[i] = tr_z_zzz_xxxx[i] * fe_0 + 4.0 * tr_z_xzzz_xxx[i] * fe_0 + tr_z_xzzz_xxxx[i] * pa_x[i];

        tr_z_xxzzz_xxxy[i] = tr_z_zzz_xxxy[i] * fe_0 + 3.0 * tr_z_xzzz_xxy[i] * fe_0 + tr_z_xzzz_xxxy[i] * pa_x[i];

        tr_z_xxzzz_xxxz[i] = tr_z_zzz_xxxz[i] * fe_0 + 3.0 * tr_z_xzzz_xxz[i] * fe_0 + tr_z_xzzz_xxxz[i] * pa_x[i];

        tr_z_xxzzz_xxyy[i] = tr_z_zzz_xxyy[i] * fe_0 + 2.0 * tr_z_xzzz_xyy[i] * fe_0 + tr_z_xzzz_xxyy[i] * pa_x[i];

        tr_z_xxzzz_xxyz[i] = tr_z_zzz_xxyz[i] * fe_0 + 2.0 * tr_z_xzzz_xyz[i] * fe_0 + tr_z_xzzz_xxyz[i] * pa_x[i];

        tr_z_xxzzz_xxzz[i] = tr_z_zzz_xxzz[i] * fe_0 + 2.0 * tr_z_xzzz_xzz[i] * fe_0 + tr_z_xzzz_xxzz[i] * pa_x[i];

        tr_z_xxzzz_xyyy[i] = tr_z_zzz_xyyy[i] * fe_0 + tr_z_xzzz_yyy[i] * fe_0 + tr_z_xzzz_xyyy[i] * pa_x[i];

        tr_z_xxzzz_xyyz[i] = tr_z_zzz_xyyz[i] * fe_0 + tr_z_xzzz_yyz[i] * fe_0 + tr_z_xzzz_xyyz[i] * pa_x[i];

        tr_z_xxzzz_xyzz[i] = tr_z_zzz_xyzz[i] * fe_0 + tr_z_xzzz_yzz[i] * fe_0 + tr_z_xzzz_xyzz[i] * pa_x[i];

        tr_z_xxzzz_xzzz[i] = tr_z_zzz_xzzz[i] * fe_0 + tr_z_xzzz_zzz[i] * fe_0 + tr_z_xzzz_xzzz[i] * pa_x[i];

        tr_z_xxzzz_yyyy[i] = tr_z_zzz_yyyy[i] * fe_0 + tr_z_xzzz_yyyy[i] * pa_x[i];

        tr_z_xxzzz_yyyz[i] = tr_z_zzz_yyyz[i] * fe_0 + tr_z_xzzz_yyyz[i] * pa_x[i];

        tr_z_xxzzz_yyzz[i] = tr_z_zzz_yyzz[i] * fe_0 + tr_z_xzzz_yyzz[i] * pa_x[i];

        tr_z_xxzzz_yzzz[i] = tr_z_zzz_yzzz[i] * fe_0 + tr_z_xzzz_yzzz[i] * pa_x[i];

        tr_z_xxzzz_zzzz[i] = tr_z_zzz_zzzz[i] * fe_0 + tr_z_xzzz_zzzz[i] * pa_x[i];
    }

    // Set up 780-795 components of targeted buffer : HG

    auto tr_z_xyyyy_xxxx = pbuffer.data(idx_dip_hg + 780);

    auto tr_z_xyyyy_xxxy = pbuffer.data(idx_dip_hg + 781);

    auto tr_z_xyyyy_xxxz = pbuffer.data(idx_dip_hg + 782);

    auto tr_z_xyyyy_xxyy = pbuffer.data(idx_dip_hg + 783);

    auto tr_z_xyyyy_xxyz = pbuffer.data(idx_dip_hg + 784);

    auto tr_z_xyyyy_xxzz = pbuffer.data(idx_dip_hg + 785);

    auto tr_z_xyyyy_xyyy = pbuffer.data(idx_dip_hg + 786);

    auto tr_z_xyyyy_xyyz = pbuffer.data(idx_dip_hg + 787);

    auto tr_z_xyyyy_xyzz = pbuffer.data(idx_dip_hg + 788);

    auto tr_z_xyyyy_xzzz = pbuffer.data(idx_dip_hg + 789);

    auto tr_z_xyyyy_yyyy = pbuffer.data(idx_dip_hg + 790);

    auto tr_z_xyyyy_yyyz = pbuffer.data(idx_dip_hg + 791);

    auto tr_z_xyyyy_yyzz = pbuffer.data(idx_dip_hg + 792);

    auto tr_z_xyyyy_yzzz = pbuffer.data(idx_dip_hg + 793);

    auto tr_z_xyyyy_zzzz = pbuffer.data(idx_dip_hg + 794);

#pragma omp simd aligned(pa_x,                \
                             tr_z_xyyyy_xxxx, \
                             tr_z_xyyyy_xxxy, \
                             tr_z_xyyyy_xxxz, \
                             tr_z_xyyyy_xxyy, \
                             tr_z_xyyyy_xxyz, \
                             tr_z_xyyyy_xxzz, \
                             tr_z_xyyyy_xyyy, \
                             tr_z_xyyyy_xyyz, \
                             tr_z_xyyyy_xyzz, \
                             tr_z_xyyyy_xzzz, \
                             tr_z_xyyyy_yyyy, \
                             tr_z_xyyyy_yyyz, \
                             tr_z_xyyyy_yyzz, \
                             tr_z_xyyyy_yzzz, \
                             tr_z_xyyyy_zzzz, \
                             tr_z_yyyy_xxx,   \
                             tr_z_yyyy_xxxx,  \
                             tr_z_yyyy_xxxy,  \
                             tr_z_yyyy_xxxz,  \
                             tr_z_yyyy_xxy,   \
                             tr_z_yyyy_xxyy,  \
                             tr_z_yyyy_xxyz,  \
                             tr_z_yyyy_xxz,   \
                             tr_z_yyyy_xxzz,  \
                             tr_z_yyyy_xyy,   \
                             tr_z_yyyy_xyyy,  \
                             tr_z_yyyy_xyyz,  \
                             tr_z_yyyy_xyz,   \
                             tr_z_yyyy_xyzz,  \
                             tr_z_yyyy_xzz,   \
                             tr_z_yyyy_xzzz,  \
                             tr_z_yyyy_yyy,   \
                             tr_z_yyyy_yyyy,  \
                             tr_z_yyyy_yyyz,  \
                             tr_z_yyyy_yyz,   \
                             tr_z_yyyy_yyzz,  \
                             tr_z_yyyy_yzz,   \
                             tr_z_yyyy_yzzz,  \
                             tr_z_yyyy_zzz,   \
                             tr_z_yyyy_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyy_xxxx[i] = 4.0 * tr_z_yyyy_xxx[i] * fe_0 + tr_z_yyyy_xxxx[i] * pa_x[i];

        tr_z_xyyyy_xxxy[i] = 3.0 * tr_z_yyyy_xxy[i] * fe_0 + tr_z_yyyy_xxxy[i] * pa_x[i];

        tr_z_xyyyy_xxxz[i] = 3.0 * tr_z_yyyy_xxz[i] * fe_0 + tr_z_yyyy_xxxz[i] * pa_x[i];

        tr_z_xyyyy_xxyy[i] = 2.0 * tr_z_yyyy_xyy[i] * fe_0 + tr_z_yyyy_xxyy[i] * pa_x[i];

        tr_z_xyyyy_xxyz[i] = 2.0 * tr_z_yyyy_xyz[i] * fe_0 + tr_z_yyyy_xxyz[i] * pa_x[i];

        tr_z_xyyyy_xxzz[i] = 2.0 * tr_z_yyyy_xzz[i] * fe_0 + tr_z_yyyy_xxzz[i] * pa_x[i];

        tr_z_xyyyy_xyyy[i] = tr_z_yyyy_yyy[i] * fe_0 + tr_z_yyyy_xyyy[i] * pa_x[i];

        tr_z_xyyyy_xyyz[i] = tr_z_yyyy_yyz[i] * fe_0 + tr_z_yyyy_xyyz[i] * pa_x[i];

        tr_z_xyyyy_xyzz[i] = tr_z_yyyy_yzz[i] * fe_0 + tr_z_yyyy_xyzz[i] * pa_x[i];

        tr_z_xyyyy_xzzz[i] = tr_z_yyyy_zzz[i] * fe_0 + tr_z_yyyy_xzzz[i] * pa_x[i];

        tr_z_xyyyy_yyyy[i] = tr_z_yyyy_yyyy[i] * pa_x[i];

        tr_z_xyyyy_yyyz[i] = tr_z_yyyy_yyyz[i] * pa_x[i];

        tr_z_xyyyy_yyzz[i] = tr_z_yyyy_yyzz[i] * pa_x[i];

        tr_z_xyyyy_yzzz[i] = tr_z_yyyy_yzzz[i] * pa_x[i];

        tr_z_xyyyy_zzzz[i] = tr_z_yyyy_zzzz[i] * pa_x[i];
    }

    // Set up 795-810 components of targeted buffer : HG

    auto tr_z_xyyyz_xxxx = pbuffer.data(idx_dip_hg + 795);

    auto tr_z_xyyyz_xxxy = pbuffer.data(idx_dip_hg + 796);

    auto tr_z_xyyyz_xxxz = pbuffer.data(idx_dip_hg + 797);

    auto tr_z_xyyyz_xxyy = pbuffer.data(idx_dip_hg + 798);

    auto tr_z_xyyyz_xxyz = pbuffer.data(idx_dip_hg + 799);

    auto tr_z_xyyyz_xxzz = pbuffer.data(idx_dip_hg + 800);

    auto tr_z_xyyyz_xyyy = pbuffer.data(idx_dip_hg + 801);

    auto tr_z_xyyyz_xyyz = pbuffer.data(idx_dip_hg + 802);

    auto tr_z_xyyyz_xyzz = pbuffer.data(idx_dip_hg + 803);

    auto tr_z_xyyyz_xzzz = pbuffer.data(idx_dip_hg + 804);

    auto tr_z_xyyyz_yyyy = pbuffer.data(idx_dip_hg + 805);

    auto tr_z_xyyyz_yyyz = pbuffer.data(idx_dip_hg + 806);

    auto tr_z_xyyyz_yyzz = pbuffer.data(idx_dip_hg + 807);

    auto tr_z_xyyyz_yzzz = pbuffer.data(idx_dip_hg + 808);

    auto tr_z_xyyyz_zzzz = pbuffer.data(idx_dip_hg + 809);

#pragma omp simd aligned(pa_x,                \
                             tr_z_xyyyz_xxxx, \
                             tr_z_xyyyz_xxxy, \
                             tr_z_xyyyz_xxxz, \
                             tr_z_xyyyz_xxyy, \
                             tr_z_xyyyz_xxyz, \
                             tr_z_xyyyz_xxzz, \
                             tr_z_xyyyz_xyyy, \
                             tr_z_xyyyz_xyyz, \
                             tr_z_xyyyz_xyzz, \
                             tr_z_xyyyz_xzzz, \
                             tr_z_xyyyz_yyyy, \
                             tr_z_xyyyz_yyyz, \
                             tr_z_xyyyz_yyzz, \
                             tr_z_xyyyz_yzzz, \
                             tr_z_xyyyz_zzzz, \
                             tr_z_yyyz_xxx,   \
                             tr_z_yyyz_xxxx,  \
                             tr_z_yyyz_xxxy,  \
                             tr_z_yyyz_xxxz,  \
                             tr_z_yyyz_xxy,   \
                             tr_z_yyyz_xxyy,  \
                             tr_z_yyyz_xxyz,  \
                             tr_z_yyyz_xxz,   \
                             tr_z_yyyz_xxzz,  \
                             tr_z_yyyz_xyy,   \
                             tr_z_yyyz_xyyy,  \
                             tr_z_yyyz_xyyz,  \
                             tr_z_yyyz_xyz,   \
                             tr_z_yyyz_xyzz,  \
                             tr_z_yyyz_xzz,   \
                             tr_z_yyyz_xzzz,  \
                             tr_z_yyyz_yyy,   \
                             tr_z_yyyz_yyyy,  \
                             tr_z_yyyz_yyyz,  \
                             tr_z_yyyz_yyz,   \
                             tr_z_yyyz_yyzz,  \
                             tr_z_yyyz_yzz,   \
                             tr_z_yyyz_yzzz,  \
                             tr_z_yyyz_zzz,   \
                             tr_z_yyyz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyz_xxxx[i] = 4.0 * tr_z_yyyz_xxx[i] * fe_0 + tr_z_yyyz_xxxx[i] * pa_x[i];

        tr_z_xyyyz_xxxy[i] = 3.0 * tr_z_yyyz_xxy[i] * fe_0 + tr_z_yyyz_xxxy[i] * pa_x[i];

        tr_z_xyyyz_xxxz[i] = 3.0 * tr_z_yyyz_xxz[i] * fe_0 + tr_z_yyyz_xxxz[i] * pa_x[i];

        tr_z_xyyyz_xxyy[i] = 2.0 * tr_z_yyyz_xyy[i] * fe_0 + tr_z_yyyz_xxyy[i] * pa_x[i];

        tr_z_xyyyz_xxyz[i] = 2.0 * tr_z_yyyz_xyz[i] * fe_0 + tr_z_yyyz_xxyz[i] * pa_x[i];

        tr_z_xyyyz_xxzz[i] = 2.0 * tr_z_yyyz_xzz[i] * fe_0 + tr_z_yyyz_xxzz[i] * pa_x[i];

        tr_z_xyyyz_xyyy[i] = tr_z_yyyz_yyy[i] * fe_0 + tr_z_yyyz_xyyy[i] * pa_x[i];

        tr_z_xyyyz_xyyz[i] = tr_z_yyyz_yyz[i] * fe_0 + tr_z_yyyz_xyyz[i] * pa_x[i];

        tr_z_xyyyz_xyzz[i] = tr_z_yyyz_yzz[i] * fe_0 + tr_z_yyyz_xyzz[i] * pa_x[i];

        tr_z_xyyyz_xzzz[i] = tr_z_yyyz_zzz[i] * fe_0 + tr_z_yyyz_xzzz[i] * pa_x[i];

        tr_z_xyyyz_yyyy[i] = tr_z_yyyz_yyyy[i] * pa_x[i];

        tr_z_xyyyz_yyyz[i] = tr_z_yyyz_yyyz[i] * pa_x[i];

        tr_z_xyyyz_yyzz[i] = tr_z_yyyz_yyzz[i] * pa_x[i];

        tr_z_xyyyz_yzzz[i] = tr_z_yyyz_yzzz[i] * pa_x[i];

        tr_z_xyyyz_zzzz[i] = tr_z_yyyz_zzzz[i] * pa_x[i];
    }

    // Set up 810-825 components of targeted buffer : HG

    auto tr_z_xyyzz_xxxx = pbuffer.data(idx_dip_hg + 810);

    auto tr_z_xyyzz_xxxy = pbuffer.data(idx_dip_hg + 811);

    auto tr_z_xyyzz_xxxz = pbuffer.data(idx_dip_hg + 812);

    auto tr_z_xyyzz_xxyy = pbuffer.data(idx_dip_hg + 813);

    auto tr_z_xyyzz_xxyz = pbuffer.data(idx_dip_hg + 814);

    auto tr_z_xyyzz_xxzz = pbuffer.data(idx_dip_hg + 815);

    auto tr_z_xyyzz_xyyy = pbuffer.data(idx_dip_hg + 816);

    auto tr_z_xyyzz_xyyz = pbuffer.data(idx_dip_hg + 817);

    auto tr_z_xyyzz_xyzz = pbuffer.data(idx_dip_hg + 818);

    auto tr_z_xyyzz_xzzz = pbuffer.data(idx_dip_hg + 819);

    auto tr_z_xyyzz_yyyy = pbuffer.data(idx_dip_hg + 820);

    auto tr_z_xyyzz_yyyz = pbuffer.data(idx_dip_hg + 821);

    auto tr_z_xyyzz_yyzz = pbuffer.data(idx_dip_hg + 822);

    auto tr_z_xyyzz_yzzz = pbuffer.data(idx_dip_hg + 823);

    auto tr_z_xyyzz_zzzz = pbuffer.data(idx_dip_hg + 824);

#pragma omp simd aligned(pa_x,                \
                             tr_z_xyyzz_xxxx, \
                             tr_z_xyyzz_xxxy, \
                             tr_z_xyyzz_xxxz, \
                             tr_z_xyyzz_xxyy, \
                             tr_z_xyyzz_xxyz, \
                             tr_z_xyyzz_xxzz, \
                             tr_z_xyyzz_xyyy, \
                             tr_z_xyyzz_xyyz, \
                             tr_z_xyyzz_xyzz, \
                             tr_z_xyyzz_xzzz, \
                             tr_z_xyyzz_yyyy, \
                             tr_z_xyyzz_yyyz, \
                             tr_z_xyyzz_yyzz, \
                             tr_z_xyyzz_yzzz, \
                             tr_z_xyyzz_zzzz, \
                             tr_z_yyzz_xxx,   \
                             tr_z_yyzz_xxxx,  \
                             tr_z_yyzz_xxxy,  \
                             tr_z_yyzz_xxxz,  \
                             tr_z_yyzz_xxy,   \
                             tr_z_yyzz_xxyy,  \
                             tr_z_yyzz_xxyz,  \
                             tr_z_yyzz_xxz,   \
                             tr_z_yyzz_xxzz,  \
                             tr_z_yyzz_xyy,   \
                             tr_z_yyzz_xyyy,  \
                             tr_z_yyzz_xyyz,  \
                             tr_z_yyzz_xyz,   \
                             tr_z_yyzz_xyzz,  \
                             tr_z_yyzz_xzz,   \
                             tr_z_yyzz_xzzz,  \
                             tr_z_yyzz_yyy,   \
                             tr_z_yyzz_yyyy,  \
                             tr_z_yyzz_yyyz,  \
                             tr_z_yyzz_yyz,   \
                             tr_z_yyzz_yyzz,  \
                             tr_z_yyzz_yzz,   \
                             tr_z_yyzz_yzzz,  \
                             tr_z_yyzz_zzz,   \
                             tr_z_yyzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyzz_xxxx[i] = 4.0 * tr_z_yyzz_xxx[i] * fe_0 + tr_z_yyzz_xxxx[i] * pa_x[i];

        tr_z_xyyzz_xxxy[i] = 3.0 * tr_z_yyzz_xxy[i] * fe_0 + tr_z_yyzz_xxxy[i] * pa_x[i];

        tr_z_xyyzz_xxxz[i] = 3.0 * tr_z_yyzz_xxz[i] * fe_0 + tr_z_yyzz_xxxz[i] * pa_x[i];

        tr_z_xyyzz_xxyy[i] = 2.0 * tr_z_yyzz_xyy[i] * fe_0 + tr_z_yyzz_xxyy[i] * pa_x[i];

        tr_z_xyyzz_xxyz[i] = 2.0 * tr_z_yyzz_xyz[i] * fe_0 + tr_z_yyzz_xxyz[i] * pa_x[i];

        tr_z_xyyzz_xxzz[i] = 2.0 * tr_z_yyzz_xzz[i] * fe_0 + tr_z_yyzz_xxzz[i] * pa_x[i];

        tr_z_xyyzz_xyyy[i] = tr_z_yyzz_yyy[i] * fe_0 + tr_z_yyzz_xyyy[i] * pa_x[i];

        tr_z_xyyzz_xyyz[i] = tr_z_yyzz_yyz[i] * fe_0 + tr_z_yyzz_xyyz[i] * pa_x[i];

        tr_z_xyyzz_xyzz[i] = tr_z_yyzz_yzz[i] * fe_0 + tr_z_yyzz_xyzz[i] * pa_x[i];

        tr_z_xyyzz_xzzz[i] = tr_z_yyzz_zzz[i] * fe_0 + tr_z_yyzz_xzzz[i] * pa_x[i];

        tr_z_xyyzz_yyyy[i] = tr_z_yyzz_yyyy[i] * pa_x[i];

        tr_z_xyyzz_yyyz[i] = tr_z_yyzz_yyyz[i] * pa_x[i];

        tr_z_xyyzz_yyzz[i] = tr_z_yyzz_yyzz[i] * pa_x[i];

        tr_z_xyyzz_yzzz[i] = tr_z_yyzz_yzzz[i] * pa_x[i];

        tr_z_xyyzz_zzzz[i] = tr_z_yyzz_zzzz[i] * pa_x[i];
    }

    // Set up 825-840 components of targeted buffer : HG

    auto tr_z_xyzzz_xxxx = pbuffer.data(idx_dip_hg + 825);

    auto tr_z_xyzzz_xxxy = pbuffer.data(idx_dip_hg + 826);

    auto tr_z_xyzzz_xxxz = pbuffer.data(idx_dip_hg + 827);

    auto tr_z_xyzzz_xxyy = pbuffer.data(idx_dip_hg + 828);

    auto tr_z_xyzzz_xxyz = pbuffer.data(idx_dip_hg + 829);

    auto tr_z_xyzzz_xxzz = pbuffer.data(idx_dip_hg + 830);

    auto tr_z_xyzzz_xyyy = pbuffer.data(idx_dip_hg + 831);

    auto tr_z_xyzzz_xyyz = pbuffer.data(idx_dip_hg + 832);

    auto tr_z_xyzzz_xyzz = pbuffer.data(idx_dip_hg + 833);

    auto tr_z_xyzzz_xzzz = pbuffer.data(idx_dip_hg + 834);

    auto tr_z_xyzzz_yyyy = pbuffer.data(idx_dip_hg + 835);

    auto tr_z_xyzzz_yyyz = pbuffer.data(idx_dip_hg + 836);

    auto tr_z_xyzzz_yyzz = pbuffer.data(idx_dip_hg + 837);

    auto tr_z_xyzzz_yzzz = pbuffer.data(idx_dip_hg + 838);

    auto tr_z_xyzzz_zzzz = pbuffer.data(idx_dip_hg + 839);

#pragma omp simd aligned(pa_x,                \
                             pa_y,            \
                             tr_z_xyzzz_xxxx, \
                             tr_z_xyzzz_xxxy, \
                             tr_z_xyzzz_xxxz, \
                             tr_z_xyzzz_xxyy, \
                             tr_z_xyzzz_xxyz, \
                             tr_z_xyzzz_xxzz, \
                             tr_z_xyzzz_xyyy, \
                             tr_z_xyzzz_xyyz, \
                             tr_z_xyzzz_xyzz, \
                             tr_z_xyzzz_xzzz, \
                             tr_z_xyzzz_yyyy, \
                             tr_z_xyzzz_yyyz, \
                             tr_z_xyzzz_yyzz, \
                             tr_z_xyzzz_yzzz, \
                             tr_z_xyzzz_zzzz, \
                             tr_z_xzzz_xxxx,  \
                             tr_z_xzzz_xxxz,  \
                             tr_z_xzzz_xxzz,  \
                             tr_z_xzzz_xzzz,  \
                             tr_z_yzzz_xxxy,  \
                             tr_z_yzzz_xxy,   \
                             tr_z_yzzz_xxyy,  \
                             tr_z_yzzz_xxyz,  \
                             tr_z_yzzz_xyy,   \
                             tr_z_yzzz_xyyy,  \
                             tr_z_yzzz_xyyz,  \
                             tr_z_yzzz_xyz,   \
                             tr_z_yzzz_xyzz,  \
                             tr_z_yzzz_yyy,   \
                             tr_z_yzzz_yyyy,  \
                             tr_z_yzzz_yyyz,  \
                             tr_z_yzzz_yyz,   \
                             tr_z_yzzz_yyzz,  \
                             tr_z_yzzz_yzz,   \
                             tr_z_yzzz_yzzz,  \
                             tr_z_yzzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyzzz_xxxx[i] = tr_z_xzzz_xxxx[i] * pa_y[i];

        tr_z_xyzzz_xxxy[i] = 3.0 * tr_z_yzzz_xxy[i] * fe_0 + tr_z_yzzz_xxxy[i] * pa_x[i];

        tr_z_xyzzz_xxxz[i] = tr_z_xzzz_xxxz[i] * pa_y[i];

        tr_z_xyzzz_xxyy[i] = 2.0 * tr_z_yzzz_xyy[i] * fe_0 + tr_z_yzzz_xxyy[i] * pa_x[i];

        tr_z_xyzzz_xxyz[i] = 2.0 * tr_z_yzzz_xyz[i] * fe_0 + tr_z_yzzz_xxyz[i] * pa_x[i];

        tr_z_xyzzz_xxzz[i] = tr_z_xzzz_xxzz[i] * pa_y[i];

        tr_z_xyzzz_xyyy[i] = tr_z_yzzz_yyy[i] * fe_0 + tr_z_yzzz_xyyy[i] * pa_x[i];

        tr_z_xyzzz_xyyz[i] = tr_z_yzzz_yyz[i] * fe_0 + tr_z_yzzz_xyyz[i] * pa_x[i];

        tr_z_xyzzz_xyzz[i] = tr_z_yzzz_yzz[i] * fe_0 + tr_z_yzzz_xyzz[i] * pa_x[i];

        tr_z_xyzzz_xzzz[i] = tr_z_xzzz_xzzz[i] * pa_y[i];

        tr_z_xyzzz_yyyy[i] = tr_z_yzzz_yyyy[i] * pa_x[i];

        tr_z_xyzzz_yyyz[i] = tr_z_yzzz_yyyz[i] * pa_x[i];

        tr_z_xyzzz_yyzz[i] = tr_z_yzzz_yyzz[i] * pa_x[i];

        tr_z_xyzzz_yzzz[i] = tr_z_yzzz_yzzz[i] * pa_x[i];

        tr_z_xyzzz_zzzz[i] = tr_z_yzzz_zzzz[i] * pa_x[i];
    }

    // Set up 840-855 components of targeted buffer : HG

    auto tr_z_xzzzz_xxxx = pbuffer.data(idx_dip_hg + 840);

    auto tr_z_xzzzz_xxxy = pbuffer.data(idx_dip_hg + 841);

    auto tr_z_xzzzz_xxxz = pbuffer.data(idx_dip_hg + 842);

    auto tr_z_xzzzz_xxyy = pbuffer.data(idx_dip_hg + 843);

    auto tr_z_xzzzz_xxyz = pbuffer.data(idx_dip_hg + 844);

    auto tr_z_xzzzz_xxzz = pbuffer.data(idx_dip_hg + 845);

    auto tr_z_xzzzz_xyyy = pbuffer.data(idx_dip_hg + 846);

    auto tr_z_xzzzz_xyyz = pbuffer.data(idx_dip_hg + 847);

    auto tr_z_xzzzz_xyzz = pbuffer.data(idx_dip_hg + 848);

    auto tr_z_xzzzz_xzzz = pbuffer.data(idx_dip_hg + 849);

    auto tr_z_xzzzz_yyyy = pbuffer.data(idx_dip_hg + 850);

    auto tr_z_xzzzz_yyyz = pbuffer.data(idx_dip_hg + 851);

    auto tr_z_xzzzz_yyzz = pbuffer.data(idx_dip_hg + 852);

    auto tr_z_xzzzz_yzzz = pbuffer.data(idx_dip_hg + 853);

    auto tr_z_xzzzz_zzzz = pbuffer.data(idx_dip_hg + 854);

#pragma omp simd aligned(pa_x,                \
                             tr_z_xzzzz_xxxx, \
                             tr_z_xzzzz_xxxy, \
                             tr_z_xzzzz_xxxz, \
                             tr_z_xzzzz_xxyy, \
                             tr_z_xzzzz_xxyz, \
                             tr_z_xzzzz_xxzz, \
                             tr_z_xzzzz_xyyy, \
                             tr_z_xzzzz_xyyz, \
                             tr_z_xzzzz_xyzz, \
                             tr_z_xzzzz_xzzz, \
                             tr_z_xzzzz_yyyy, \
                             tr_z_xzzzz_yyyz, \
                             tr_z_xzzzz_yyzz, \
                             tr_z_xzzzz_yzzz, \
                             tr_z_xzzzz_zzzz, \
                             tr_z_zzzz_xxx,   \
                             tr_z_zzzz_xxxx,  \
                             tr_z_zzzz_xxxy,  \
                             tr_z_zzzz_xxxz,  \
                             tr_z_zzzz_xxy,   \
                             tr_z_zzzz_xxyy,  \
                             tr_z_zzzz_xxyz,  \
                             tr_z_zzzz_xxz,   \
                             tr_z_zzzz_xxzz,  \
                             tr_z_zzzz_xyy,   \
                             tr_z_zzzz_xyyy,  \
                             tr_z_zzzz_xyyz,  \
                             tr_z_zzzz_xyz,   \
                             tr_z_zzzz_xyzz,  \
                             tr_z_zzzz_xzz,   \
                             tr_z_zzzz_xzzz,  \
                             tr_z_zzzz_yyy,   \
                             tr_z_zzzz_yyyy,  \
                             tr_z_zzzz_yyyz,  \
                             tr_z_zzzz_yyz,   \
                             tr_z_zzzz_yyzz,  \
                             tr_z_zzzz_yzz,   \
                             tr_z_zzzz_yzzz,  \
                             tr_z_zzzz_zzz,   \
                             tr_z_zzzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzzz_xxxx[i] = 4.0 * tr_z_zzzz_xxx[i] * fe_0 + tr_z_zzzz_xxxx[i] * pa_x[i];

        tr_z_xzzzz_xxxy[i] = 3.0 * tr_z_zzzz_xxy[i] * fe_0 + tr_z_zzzz_xxxy[i] * pa_x[i];

        tr_z_xzzzz_xxxz[i] = 3.0 * tr_z_zzzz_xxz[i] * fe_0 + tr_z_zzzz_xxxz[i] * pa_x[i];

        tr_z_xzzzz_xxyy[i] = 2.0 * tr_z_zzzz_xyy[i] * fe_0 + tr_z_zzzz_xxyy[i] * pa_x[i];

        tr_z_xzzzz_xxyz[i] = 2.0 * tr_z_zzzz_xyz[i] * fe_0 + tr_z_zzzz_xxyz[i] * pa_x[i];

        tr_z_xzzzz_xxzz[i] = 2.0 * tr_z_zzzz_xzz[i] * fe_0 + tr_z_zzzz_xxzz[i] * pa_x[i];

        tr_z_xzzzz_xyyy[i] = tr_z_zzzz_yyy[i] * fe_0 + tr_z_zzzz_xyyy[i] * pa_x[i];

        tr_z_xzzzz_xyyz[i] = tr_z_zzzz_yyz[i] * fe_0 + tr_z_zzzz_xyyz[i] * pa_x[i];

        tr_z_xzzzz_xyzz[i] = tr_z_zzzz_yzz[i] * fe_0 + tr_z_zzzz_xyzz[i] * pa_x[i];

        tr_z_xzzzz_xzzz[i] = tr_z_zzzz_zzz[i] * fe_0 + tr_z_zzzz_xzzz[i] * pa_x[i];

        tr_z_xzzzz_yyyy[i] = tr_z_zzzz_yyyy[i] * pa_x[i];

        tr_z_xzzzz_yyyz[i] = tr_z_zzzz_yyyz[i] * pa_x[i];

        tr_z_xzzzz_yyzz[i] = tr_z_zzzz_yyzz[i] * pa_x[i];

        tr_z_xzzzz_yzzz[i] = tr_z_zzzz_yzzz[i] * pa_x[i];

        tr_z_xzzzz_zzzz[i] = tr_z_zzzz_zzzz[i] * pa_x[i];
    }

    // Set up 855-870 components of targeted buffer : HG

    auto tr_z_yyyyy_xxxx = pbuffer.data(idx_dip_hg + 855);

    auto tr_z_yyyyy_xxxy = pbuffer.data(idx_dip_hg + 856);

    auto tr_z_yyyyy_xxxz = pbuffer.data(idx_dip_hg + 857);

    auto tr_z_yyyyy_xxyy = pbuffer.data(idx_dip_hg + 858);

    auto tr_z_yyyyy_xxyz = pbuffer.data(idx_dip_hg + 859);

    auto tr_z_yyyyy_xxzz = pbuffer.data(idx_dip_hg + 860);

    auto tr_z_yyyyy_xyyy = pbuffer.data(idx_dip_hg + 861);

    auto tr_z_yyyyy_xyyz = pbuffer.data(idx_dip_hg + 862);

    auto tr_z_yyyyy_xyzz = pbuffer.data(idx_dip_hg + 863);

    auto tr_z_yyyyy_xzzz = pbuffer.data(idx_dip_hg + 864);

    auto tr_z_yyyyy_yyyy = pbuffer.data(idx_dip_hg + 865);

    auto tr_z_yyyyy_yyyz = pbuffer.data(idx_dip_hg + 866);

    auto tr_z_yyyyy_yyzz = pbuffer.data(idx_dip_hg + 867);

    auto tr_z_yyyyy_yzzz = pbuffer.data(idx_dip_hg + 868);

    auto tr_z_yyyyy_zzzz = pbuffer.data(idx_dip_hg + 869);

#pragma omp simd aligned(pa_y,                \
                             tr_z_yyy_xxxx,   \
                             tr_z_yyy_xxxy,   \
                             tr_z_yyy_xxxz,   \
                             tr_z_yyy_xxyy,   \
                             tr_z_yyy_xxyz,   \
                             tr_z_yyy_xxzz,   \
                             tr_z_yyy_xyyy,   \
                             tr_z_yyy_xyyz,   \
                             tr_z_yyy_xyzz,   \
                             tr_z_yyy_xzzz,   \
                             tr_z_yyy_yyyy,   \
                             tr_z_yyy_yyyz,   \
                             tr_z_yyy_yyzz,   \
                             tr_z_yyy_yzzz,   \
                             tr_z_yyy_zzzz,   \
                             tr_z_yyyy_xxx,   \
                             tr_z_yyyy_xxxx,  \
                             tr_z_yyyy_xxxy,  \
                             tr_z_yyyy_xxxz,  \
                             tr_z_yyyy_xxy,   \
                             tr_z_yyyy_xxyy,  \
                             tr_z_yyyy_xxyz,  \
                             tr_z_yyyy_xxz,   \
                             tr_z_yyyy_xxzz,  \
                             tr_z_yyyy_xyy,   \
                             tr_z_yyyy_xyyy,  \
                             tr_z_yyyy_xyyz,  \
                             tr_z_yyyy_xyz,   \
                             tr_z_yyyy_xyzz,  \
                             tr_z_yyyy_xzz,   \
                             tr_z_yyyy_xzzz,  \
                             tr_z_yyyy_yyy,   \
                             tr_z_yyyy_yyyy,  \
                             tr_z_yyyy_yyyz,  \
                             tr_z_yyyy_yyz,   \
                             tr_z_yyyy_yyzz,  \
                             tr_z_yyyy_yzz,   \
                             tr_z_yyyy_yzzz,  \
                             tr_z_yyyy_zzz,   \
                             tr_z_yyyy_zzzz,  \
                             tr_z_yyyyy_xxxx, \
                             tr_z_yyyyy_xxxy, \
                             tr_z_yyyyy_xxxz, \
                             tr_z_yyyyy_xxyy, \
                             tr_z_yyyyy_xxyz, \
                             tr_z_yyyyy_xxzz, \
                             tr_z_yyyyy_xyyy, \
                             tr_z_yyyyy_xyyz, \
                             tr_z_yyyyy_xyzz, \
                             tr_z_yyyyy_xzzz, \
                             tr_z_yyyyy_yyyy, \
                             tr_z_yyyyy_yyyz, \
                             tr_z_yyyyy_yyzz, \
                             tr_z_yyyyy_yzzz, \
                             tr_z_yyyyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyy_xxxx[i] = 4.0 * tr_z_yyy_xxxx[i] * fe_0 + tr_z_yyyy_xxxx[i] * pa_y[i];

        tr_z_yyyyy_xxxy[i] = 4.0 * tr_z_yyy_xxxy[i] * fe_0 + tr_z_yyyy_xxx[i] * fe_0 + tr_z_yyyy_xxxy[i] * pa_y[i];

        tr_z_yyyyy_xxxz[i] = 4.0 * tr_z_yyy_xxxz[i] * fe_0 + tr_z_yyyy_xxxz[i] * pa_y[i];

        tr_z_yyyyy_xxyy[i] = 4.0 * tr_z_yyy_xxyy[i] * fe_0 + 2.0 * tr_z_yyyy_xxy[i] * fe_0 + tr_z_yyyy_xxyy[i] * pa_y[i];

        tr_z_yyyyy_xxyz[i] = 4.0 * tr_z_yyy_xxyz[i] * fe_0 + tr_z_yyyy_xxz[i] * fe_0 + tr_z_yyyy_xxyz[i] * pa_y[i];

        tr_z_yyyyy_xxzz[i] = 4.0 * tr_z_yyy_xxzz[i] * fe_0 + tr_z_yyyy_xxzz[i] * pa_y[i];

        tr_z_yyyyy_xyyy[i] = 4.0 * tr_z_yyy_xyyy[i] * fe_0 + 3.0 * tr_z_yyyy_xyy[i] * fe_0 + tr_z_yyyy_xyyy[i] * pa_y[i];

        tr_z_yyyyy_xyyz[i] = 4.0 * tr_z_yyy_xyyz[i] * fe_0 + 2.0 * tr_z_yyyy_xyz[i] * fe_0 + tr_z_yyyy_xyyz[i] * pa_y[i];

        tr_z_yyyyy_xyzz[i] = 4.0 * tr_z_yyy_xyzz[i] * fe_0 + tr_z_yyyy_xzz[i] * fe_0 + tr_z_yyyy_xyzz[i] * pa_y[i];

        tr_z_yyyyy_xzzz[i] = 4.0 * tr_z_yyy_xzzz[i] * fe_0 + tr_z_yyyy_xzzz[i] * pa_y[i];

        tr_z_yyyyy_yyyy[i] = 4.0 * tr_z_yyy_yyyy[i] * fe_0 + 4.0 * tr_z_yyyy_yyy[i] * fe_0 + tr_z_yyyy_yyyy[i] * pa_y[i];

        tr_z_yyyyy_yyyz[i] = 4.0 * tr_z_yyy_yyyz[i] * fe_0 + 3.0 * tr_z_yyyy_yyz[i] * fe_0 + tr_z_yyyy_yyyz[i] * pa_y[i];

        tr_z_yyyyy_yyzz[i] = 4.0 * tr_z_yyy_yyzz[i] * fe_0 + 2.0 * tr_z_yyyy_yzz[i] * fe_0 + tr_z_yyyy_yyzz[i] * pa_y[i];

        tr_z_yyyyy_yzzz[i] = 4.0 * tr_z_yyy_yzzz[i] * fe_0 + tr_z_yyyy_zzz[i] * fe_0 + tr_z_yyyy_yzzz[i] * pa_y[i];

        tr_z_yyyyy_zzzz[i] = 4.0 * tr_z_yyy_zzzz[i] * fe_0 + tr_z_yyyy_zzzz[i] * pa_y[i];
    }

    // Set up 870-885 components of targeted buffer : HG

    auto tr_z_yyyyz_xxxx = pbuffer.data(idx_dip_hg + 870);

    auto tr_z_yyyyz_xxxy = pbuffer.data(idx_dip_hg + 871);

    auto tr_z_yyyyz_xxxz = pbuffer.data(idx_dip_hg + 872);

    auto tr_z_yyyyz_xxyy = pbuffer.data(idx_dip_hg + 873);

    auto tr_z_yyyyz_xxyz = pbuffer.data(idx_dip_hg + 874);

    auto tr_z_yyyyz_xxzz = pbuffer.data(idx_dip_hg + 875);

    auto tr_z_yyyyz_xyyy = pbuffer.data(idx_dip_hg + 876);

    auto tr_z_yyyyz_xyyz = pbuffer.data(idx_dip_hg + 877);

    auto tr_z_yyyyz_xyzz = pbuffer.data(idx_dip_hg + 878);

    auto tr_z_yyyyz_xzzz = pbuffer.data(idx_dip_hg + 879);

    auto tr_z_yyyyz_yyyy = pbuffer.data(idx_dip_hg + 880);

    auto tr_z_yyyyz_yyyz = pbuffer.data(idx_dip_hg + 881);

    auto tr_z_yyyyz_yyzz = pbuffer.data(idx_dip_hg + 882);

    auto tr_z_yyyyz_yzzz = pbuffer.data(idx_dip_hg + 883);

    auto tr_z_yyyyz_zzzz = pbuffer.data(idx_dip_hg + 884);

#pragma omp simd aligned(pa_y,                \
                             pa_z,            \
                             tr_z_yyyy_xxxy,  \
                             tr_z_yyyy_xxyy,  \
                             tr_z_yyyy_xyyy,  \
                             tr_z_yyyy_yyyy,  \
                             tr_z_yyyyz_xxxx, \
                             tr_z_yyyyz_xxxy, \
                             tr_z_yyyyz_xxxz, \
                             tr_z_yyyyz_xxyy, \
                             tr_z_yyyyz_xxyz, \
                             tr_z_yyyyz_xxzz, \
                             tr_z_yyyyz_xyyy, \
                             tr_z_yyyyz_xyyz, \
                             tr_z_yyyyz_xyzz, \
                             tr_z_yyyyz_xzzz, \
                             tr_z_yyyyz_yyyy, \
                             tr_z_yyyyz_yyyz, \
                             tr_z_yyyyz_yyzz, \
                             tr_z_yyyyz_yzzz, \
                             tr_z_yyyyz_zzzz, \
                             tr_z_yyyz_xxxx,  \
                             tr_z_yyyz_xxxz,  \
                             tr_z_yyyz_xxyz,  \
                             tr_z_yyyz_xxz,   \
                             tr_z_yyyz_xxzz,  \
                             tr_z_yyyz_xyyz,  \
                             tr_z_yyyz_xyz,   \
                             tr_z_yyyz_xyzz,  \
                             tr_z_yyyz_xzz,   \
                             tr_z_yyyz_xzzz,  \
                             tr_z_yyyz_yyyz,  \
                             tr_z_yyyz_yyz,   \
                             tr_z_yyyz_yyzz,  \
                             tr_z_yyyz_yzz,   \
                             tr_z_yyyz_yzzz,  \
                             tr_z_yyyz_zzz,   \
                             tr_z_yyyz_zzzz,  \
                             tr_z_yyz_xxxx,   \
                             tr_z_yyz_xxxz,   \
                             tr_z_yyz_xxyz,   \
                             tr_z_yyz_xxzz,   \
                             tr_z_yyz_xyyz,   \
                             tr_z_yyz_xyzz,   \
                             tr_z_yyz_xzzz,   \
                             tr_z_yyz_yyyz,   \
                             tr_z_yyz_yyzz,   \
                             tr_z_yyz_yzzz,   \
                             tr_z_yyz_zzzz,   \
                             ts_yyyy_xxxy,    \
                             ts_yyyy_xxyy,    \
                             ts_yyyy_xyyy,    \
                             ts_yyyy_yyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyz_xxxx[i] = 3.0 * tr_z_yyz_xxxx[i] * fe_0 + tr_z_yyyz_xxxx[i] * pa_y[i];

        tr_z_yyyyz_xxxy[i] = ts_yyyy_xxxy[i] * fe_0 + tr_z_yyyy_xxxy[i] * pa_z[i];

        tr_z_yyyyz_xxxz[i] = 3.0 * tr_z_yyz_xxxz[i] * fe_0 + tr_z_yyyz_xxxz[i] * pa_y[i];

        tr_z_yyyyz_xxyy[i] = ts_yyyy_xxyy[i] * fe_0 + tr_z_yyyy_xxyy[i] * pa_z[i];

        tr_z_yyyyz_xxyz[i] = 3.0 * tr_z_yyz_xxyz[i] * fe_0 + tr_z_yyyz_xxz[i] * fe_0 + tr_z_yyyz_xxyz[i] * pa_y[i];

        tr_z_yyyyz_xxzz[i] = 3.0 * tr_z_yyz_xxzz[i] * fe_0 + tr_z_yyyz_xxzz[i] * pa_y[i];

        tr_z_yyyyz_xyyy[i] = ts_yyyy_xyyy[i] * fe_0 + tr_z_yyyy_xyyy[i] * pa_z[i];

        tr_z_yyyyz_xyyz[i] = 3.0 * tr_z_yyz_xyyz[i] * fe_0 + 2.0 * tr_z_yyyz_xyz[i] * fe_0 + tr_z_yyyz_xyyz[i] * pa_y[i];

        tr_z_yyyyz_xyzz[i] = 3.0 * tr_z_yyz_xyzz[i] * fe_0 + tr_z_yyyz_xzz[i] * fe_0 + tr_z_yyyz_xyzz[i] * pa_y[i];

        tr_z_yyyyz_xzzz[i] = 3.0 * tr_z_yyz_xzzz[i] * fe_0 + tr_z_yyyz_xzzz[i] * pa_y[i];

        tr_z_yyyyz_yyyy[i] = ts_yyyy_yyyy[i] * fe_0 + tr_z_yyyy_yyyy[i] * pa_z[i];

        tr_z_yyyyz_yyyz[i] = 3.0 * tr_z_yyz_yyyz[i] * fe_0 + 3.0 * tr_z_yyyz_yyz[i] * fe_0 + tr_z_yyyz_yyyz[i] * pa_y[i];

        tr_z_yyyyz_yyzz[i] = 3.0 * tr_z_yyz_yyzz[i] * fe_0 + 2.0 * tr_z_yyyz_yzz[i] * fe_0 + tr_z_yyyz_yyzz[i] * pa_y[i];

        tr_z_yyyyz_yzzz[i] = 3.0 * tr_z_yyz_yzzz[i] * fe_0 + tr_z_yyyz_zzz[i] * fe_0 + tr_z_yyyz_yzzz[i] * pa_y[i];

        tr_z_yyyyz_zzzz[i] = 3.0 * tr_z_yyz_zzzz[i] * fe_0 + tr_z_yyyz_zzzz[i] * pa_y[i];
    }

    // Set up 885-900 components of targeted buffer : HG

    auto tr_z_yyyzz_xxxx = pbuffer.data(idx_dip_hg + 885);

    auto tr_z_yyyzz_xxxy = pbuffer.data(idx_dip_hg + 886);

    auto tr_z_yyyzz_xxxz = pbuffer.data(idx_dip_hg + 887);

    auto tr_z_yyyzz_xxyy = pbuffer.data(idx_dip_hg + 888);

    auto tr_z_yyyzz_xxyz = pbuffer.data(idx_dip_hg + 889);

    auto tr_z_yyyzz_xxzz = pbuffer.data(idx_dip_hg + 890);

    auto tr_z_yyyzz_xyyy = pbuffer.data(idx_dip_hg + 891);

    auto tr_z_yyyzz_xyyz = pbuffer.data(idx_dip_hg + 892);

    auto tr_z_yyyzz_xyzz = pbuffer.data(idx_dip_hg + 893);

    auto tr_z_yyyzz_xzzz = pbuffer.data(idx_dip_hg + 894);

    auto tr_z_yyyzz_yyyy = pbuffer.data(idx_dip_hg + 895);

    auto tr_z_yyyzz_yyyz = pbuffer.data(idx_dip_hg + 896);

    auto tr_z_yyyzz_yyzz = pbuffer.data(idx_dip_hg + 897);

    auto tr_z_yyyzz_yzzz = pbuffer.data(idx_dip_hg + 898);

    auto tr_z_yyyzz_zzzz = pbuffer.data(idx_dip_hg + 899);

#pragma omp simd aligned(pa_y,                \
                             tr_z_yyyzz_xxxx, \
                             tr_z_yyyzz_xxxy, \
                             tr_z_yyyzz_xxxz, \
                             tr_z_yyyzz_xxyy, \
                             tr_z_yyyzz_xxyz, \
                             tr_z_yyyzz_xxzz, \
                             tr_z_yyyzz_xyyy, \
                             tr_z_yyyzz_xyyz, \
                             tr_z_yyyzz_xyzz, \
                             tr_z_yyyzz_xzzz, \
                             tr_z_yyyzz_yyyy, \
                             tr_z_yyyzz_yyyz, \
                             tr_z_yyyzz_yyzz, \
                             tr_z_yyyzz_yzzz, \
                             tr_z_yyyzz_zzzz, \
                             tr_z_yyzz_xxx,   \
                             tr_z_yyzz_xxxx,  \
                             tr_z_yyzz_xxxy,  \
                             tr_z_yyzz_xxxz,  \
                             tr_z_yyzz_xxy,   \
                             tr_z_yyzz_xxyy,  \
                             tr_z_yyzz_xxyz,  \
                             tr_z_yyzz_xxz,   \
                             tr_z_yyzz_xxzz,  \
                             tr_z_yyzz_xyy,   \
                             tr_z_yyzz_xyyy,  \
                             tr_z_yyzz_xyyz,  \
                             tr_z_yyzz_xyz,   \
                             tr_z_yyzz_xyzz,  \
                             tr_z_yyzz_xzz,   \
                             tr_z_yyzz_xzzz,  \
                             tr_z_yyzz_yyy,   \
                             tr_z_yyzz_yyyy,  \
                             tr_z_yyzz_yyyz,  \
                             tr_z_yyzz_yyz,   \
                             tr_z_yyzz_yyzz,  \
                             tr_z_yyzz_yzz,   \
                             tr_z_yyzz_yzzz,  \
                             tr_z_yyzz_zzz,   \
                             tr_z_yyzz_zzzz,  \
                             tr_z_yzz_xxxx,   \
                             tr_z_yzz_xxxy,   \
                             tr_z_yzz_xxxz,   \
                             tr_z_yzz_xxyy,   \
                             tr_z_yzz_xxyz,   \
                             tr_z_yzz_xxzz,   \
                             tr_z_yzz_xyyy,   \
                             tr_z_yzz_xyyz,   \
                             tr_z_yzz_xyzz,   \
                             tr_z_yzz_xzzz,   \
                             tr_z_yzz_yyyy,   \
                             tr_z_yzz_yyyz,   \
                             tr_z_yzz_yyzz,   \
                             tr_z_yzz_yzzz,   \
                             tr_z_yzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyzz_xxxx[i] = 2.0 * tr_z_yzz_xxxx[i] * fe_0 + tr_z_yyzz_xxxx[i] * pa_y[i];

        tr_z_yyyzz_xxxy[i] = 2.0 * tr_z_yzz_xxxy[i] * fe_0 + tr_z_yyzz_xxx[i] * fe_0 + tr_z_yyzz_xxxy[i] * pa_y[i];

        tr_z_yyyzz_xxxz[i] = 2.0 * tr_z_yzz_xxxz[i] * fe_0 + tr_z_yyzz_xxxz[i] * pa_y[i];

        tr_z_yyyzz_xxyy[i] = 2.0 * tr_z_yzz_xxyy[i] * fe_0 + 2.0 * tr_z_yyzz_xxy[i] * fe_0 + tr_z_yyzz_xxyy[i] * pa_y[i];

        tr_z_yyyzz_xxyz[i] = 2.0 * tr_z_yzz_xxyz[i] * fe_0 + tr_z_yyzz_xxz[i] * fe_0 + tr_z_yyzz_xxyz[i] * pa_y[i];

        tr_z_yyyzz_xxzz[i] = 2.0 * tr_z_yzz_xxzz[i] * fe_0 + tr_z_yyzz_xxzz[i] * pa_y[i];

        tr_z_yyyzz_xyyy[i] = 2.0 * tr_z_yzz_xyyy[i] * fe_0 + 3.0 * tr_z_yyzz_xyy[i] * fe_0 + tr_z_yyzz_xyyy[i] * pa_y[i];

        tr_z_yyyzz_xyyz[i] = 2.0 * tr_z_yzz_xyyz[i] * fe_0 + 2.0 * tr_z_yyzz_xyz[i] * fe_0 + tr_z_yyzz_xyyz[i] * pa_y[i];

        tr_z_yyyzz_xyzz[i] = 2.0 * tr_z_yzz_xyzz[i] * fe_0 + tr_z_yyzz_xzz[i] * fe_0 + tr_z_yyzz_xyzz[i] * pa_y[i];

        tr_z_yyyzz_xzzz[i] = 2.0 * tr_z_yzz_xzzz[i] * fe_0 + tr_z_yyzz_xzzz[i] * pa_y[i];

        tr_z_yyyzz_yyyy[i] = 2.0 * tr_z_yzz_yyyy[i] * fe_0 + 4.0 * tr_z_yyzz_yyy[i] * fe_0 + tr_z_yyzz_yyyy[i] * pa_y[i];

        tr_z_yyyzz_yyyz[i] = 2.0 * tr_z_yzz_yyyz[i] * fe_0 + 3.0 * tr_z_yyzz_yyz[i] * fe_0 + tr_z_yyzz_yyyz[i] * pa_y[i];

        tr_z_yyyzz_yyzz[i] = 2.0 * tr_z_yzz_yyzz[i] * fe_0 + 2.0 * tr_z_yyzz_yzz[i] * fe_0 + tr_z_yyzz_yyzz[i] * pa_y[i];

        tr_z_yyyzz_yzzz[i] = 2.0 * tr_z_yzz_yzzz[i] * fe_0 + tr_z_yyzz_zzz[i] * fe_0 + tr_z_yyzz_yzzz[i] * pa_y[i];

        tr_z_yyyzz_zzzz[i] = 2.0 * tr_z_yzz_zzzz[i] * fe_0 + tr_z_yyzz_zzzz[i] * pa_y[i];
    }

    // Set up 900-915 components of targeted buffer : HG

    auto tr_z_yyzzz_xxxx = pbuffer.data(idx_dip_hg + 900);

    auto tr_z_yyzzz_xxxy = pbuffer.data(idx_dip_hg + 901);

    auto tr_z_yyzzz_xxxz = pbuffer.data(idx_dip_hg + 902);

    auto tr_z_yyzzz_xxyy = pbuffer.data(idx_dip_hg + 903);

    auto tr_z_yyzzz_xxyz = pbuffer.data(idx_dip_hg + 904);

    auto tr_z_yyzzz_xxzz = pbuffer.data(idx_dip_hg + 905);

    auto tr_z_yyzzz_xyyy = pbuffer.data(idx_dip_hg + 906);

    auto tr_z_yyzzz_xyyz = pbuffer.data(idx_dip_hg + 907);

    auto tr_z_yyzzz_xyzz = pbuffer.data(idx_dip_hg + 908);

    auto tr_z_yyzzz_xzzz = pbuffer.data(idx_dip_hg + 909);

    auto tr_z_yyzzz_yyyy = pbuffer.data(idx_dip_hg + 910);

    auto tr_z_yyzzz_yyyz = pbuffer.data(idx_dip_hg + 911);

    auto tr_z_yyzzz_yyzz = pbuffer.data(idx_dip_hg + 912);

    auto tr_z_yyzzz_yzzz = pbuffer.data(idx_dip_hg + 913);

    auto tr_z_yyzzz_zzzz = pbuffer.data(idx_dip_hg + 914);

#pragma omp simd aligned(pa_y,                \
                             tr_z_yyzzz_xxxx, \
                             tr_z_yyzzz_xxxy, \
                             tr_z_yyzzz_xxxz, \
                             tr_z_yyzzz_xxyy, \
                             tr_z_yyzzz_xxyz, \
                             tr_z_yyzzz_xxzz, \
                             tr_z_yyzzz_xyyy, \
                             tr_z_yyzzz_xyyz, \
                             tr_z_yyzzz_xyzz, \
                             tr_z_yyzzz_xzzz, \
                             tr_z_yyzzz_yyyy, \
                             tr_z_yyzzz_yyyz, \
                             tr_z_yyzzz_yyzz, \
                             tr_z_yyzzz_yzzz, \
                             tr_z_yyzzz_zzzz, \
                             tr_z_yzzz_xxx,   \
                             tr_z_yzzz_xxxx,  \
                             tr_z_yzzz_xxxy,  \
                             tr_z_yzzz_xxxz,  \
                             tr_z_yzzz_xxy,   \
                             tr_z_yzzz_xxyy,  \
                             tr_z_yzzz_xxyz,  \
                             tr_z_yzzz_xxz,   \
                             tr_z_yzzz_xxzz,  \
                             tr_z_yzzz_xyy,   \
                             tr_z_yzzz_xyyy,  \
                             tr_z_yzzz_xyyz,  \
                             tr_z_yzzz_xyz,   \
                             tr_z_yzzz_xyzz,  \
                             tr_z_yzzz_xzz,   \
                             tr_z_yzzz_xzzz,  \
                             tr_z_yzzz_yyy,   \
                             tr_z_yzzz_yyyy,  \
                             tr_z_yzzz_yyyz,  \
                             tr_z_yzzz_yyz,   \
                             tr_z_yzzz_yyzz,  \
                             tr_z_yzzz_yzz,   \
                             tr_z_yzzz_yzzz,  \
                             tr_z_yzzz_zzz,   \
                             tr_z_yzzz_zzzz,  \
                             tr_z_zzz_xxxx,   \
                             tr_z_zzz_xxxy,   \
                             tr_z_zzz_xxxz,   \
                             tr_z_zzz_xxyy,   \
                             tr_z_zzz_xxyz,   \
                             tr_z_zzz_xxzz,   \
                             tr_z_zzz_xyyy,   \
                             tr_z_zzz_xyyz,   \
                             tr_z_zzz_xyzz,   \
                             tr_z_zzz_xzzz,   \
                             tr_z_zzz_yyyy,   \
                             tr_z_zzz_yyyz,   \
                             tr_z_zzz_yyzz,   \
                             tr_z_zzz_yzzz,   \
                             tr_z_zzz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzzz_xxxx[i] = tr_z_zzz_xxxx[i] * fe_0 + tr_z_yzzz_xxxx[i] * pa_y[i];

        tr_z_yyzzz_xxxy[i] = tr_z_zzz_xxxy[i] * fe_0 + tr_z_yzzz_xxx[i] * fe_0 + tr_z_yzzz_xxxy[i] * pa_y[i];

        tr_z_yyzzz_xxxz[i] = tr_z_zzz_xxxz[i] * fe_0 + tr_z_yzzz_xxxz[i] * pa_y[i];

        tr_z_yyzzz_xxyy[i] = tr_z_zzz_xxyy[i] * fe_0 + 2.0 * tr_z_yzzz_xxy[i] * fe_0 + tr_z_yzzz_xxyy[i] * pa_y[i];

        tr_z_yyzzz_xxyz[i] = tr_z_zzz_xxyz[i] * fe_0 + tr_z_yzzz_xxz[i] * fe_0 + tr_z_yzzz_xxyz[i] * pa_y[i];

        tr_z_yyzzz_xxzz[i] = tr_z_zzz_xxzz[i] * fe_0 + tr_z_yzzz_xxzz[i] * pa_y[i];

        tr_z_yyzzz_xyyy[i] = tr_z_zzz_xyyy[i] * fe_0 + 3.0 * tr_z_yzzz_xyy[i] * fe_0 + tr_z_yzzz_xyyy[i] * pa_y[i];

        tr_z_yyzzz_xyyz[i] = tr_z_zzz_xyyz[i] * fe_0 + 2.0 * tr_z_yzzz_xyz[i] * fe_0 + tr_z_yzzz_xyyz[i] * pa_y[i];

        tr_z_yyzzz_xyzz[i] = tr_z_zzz_xyzz[i] * fe_0 + tr_z_yzzz_xzz[i] * fe_0 + tr_z_yzzz_xyzz[i] * pa_y[i];

        tr_z_yyzzz_xzzz[i] = tr_z_zzz_xzzz[i] * fe_0 + tr_z_yzzz_xzzz[i] * pa_y[i];

        tr_z_yyzzz_yyyy[i] = tr_z_zzz_yyyy[i] * fe_0 + 4.0 * tr_z_yzzz_yyy[i] * fe_0 + tr_z_yzzz_yyyy[i] * pa_y[i];

        tr_z_yyzzz_yyyz[i] = tr_z_zzz_yyyz[i] * fe_0 + 3.0 * tr_z_yzzz_yyz[i] * fe_0 + tr_z_yzzz_yyyz[i] * pa_y[i];

        tr_z_yyzzz_yyzz[i] = tr_z_zzz_yyzz[i] * fe_0 + 2.0 * tr_z_yzzz_yzz[i] * fe_0 + tr_z_yzzz_yyzz[i] * pa_y[i];

        tr_z_yyzzz_yzzz[i] = tr_z_zzz_yzzz[i] * fe_0 + tr_z_yzzz_zzz[i] * fe_0 + tr_z_yzzz_yzzz[i] * pa_y[i];

        tr_z_yyzzz_zzzz[i] = tr_z_zzz_zzzz[i] * fe_0 + tr_z_yzzz_zzzz[i] * pa_y[i];
    }

    // Set up 915-930 components of targeted buffer : HG

    auto tr_z_yzzzz_xxxx = pbuffer.data(idx_dip_hg + 915);

    auto tr_z_yzzzz_xxxy = pbuffer.data(idx_dip_hg + 916);

    auto tr_z_yzzzz_xxxz = pbuffer.data(idx_dip_hg + 917);

    auto tr_z_yzzzz_xxyy = pbuffer.data(idx_dip_hg + 918);

    auto tr_z_yzzzz_xxyz = pbuffer.data(idx_dip_hg + 919);

    auto tr_z_yzzzz_xxzz = pbuffer.data(idx_dip_hg + 920);

    auto tr_z_yzzzz_xyyy = pbuffer.data(idx_dip_hg + 921);

    auto tr_z_yzzzz_xyyz = pbuffer.data(idx_dip_hg + 922);

    auto tr_z_yzzzz_xyzz = pbuffer.data(idx_dip_hg + 923);

    auto tr_z_yzzzz_xzzz = pbuffer.data(idx_dip_hg + 924);

    auto tr_z_yzzzz_yyyy = pbuffer.data(idx_dip_hg + 925);

    auto tr_z_yzzzz_yyyz = pbuffer.data(idx_dip_hg + 926);

    auto tr_z_yzzzz_yyzz = pbuffer.data(idx_dip_hg + 927);

    auto tr_z_yzzzz_yzzz = pbuffer.data(idx_dip_hg + 928);

    auto tr_z_yzzzz_zzzz = pbuffer.data(idx_dip_hg + 929);

#pragma omp simd aligned(pa_y,                \
                             tr_z_yzzzz_xxxx, \
                             tr_z_yzzzz_xxxy, \
                             tr_z_yzzzz_xxxz, \
                             tr_z_yzzzz_xxyy, \
                             tr_z_yzzzz_xxyz, \
                             tr_z_yzzzz_xxzz, \
                             tr_z_yzzzz_xyyy, \
                             tr_z_yzzzz_xyyz, \
                             tr_z_yzzzz_xyzz, \
                             tr_z_yzzzz_xzzz, \
                             tr_z_yzzzz_yyyy, \
                             tr_z_yzzzz_yyyz, \
                             tr_z_yzzzz_yyzz, \
                             tr_z_yzzzz_yzzz, \
                             tr_z_yzzzz_zzzz, \
                             tr_z_zzzz_xxx,   \
                             tr_z_zzzz_xxxx,  \
                             tr_z_zzzz_xxxy,  \
                             tr_z_zzzz_xxxz,  \
                             tr_z_zzzz_xxy,   \
                             tr_z_zzzz_xxyy,  \
                             tr_z_zzzz_xxyz,  \
                             tr_z_zzzz_xxz,   \
                             tr_z_zzzz_xxzz,  \
                             tr_z_zzzz_xyy,   \
                             tr_z_zzzz_xyyy,  \
                             tr_z_zzzz_xyyz,  \
                             tr_z_zzzz_xyz,   \
                             tr_z_zzzz_xyzz,  \
                             tr_z_zzzz_xzz,   \
                             tr_z_zzzz_xzzz,  \
                             tr_z_zzzz_yyy,   \
                             tr_z_zzzz_yyyy,  \
                             tr_z_zzzz_yyyz,  \
                             tr_z_zzzz_yyz,   \
                             tr_z_zzzz_yyzz,  \
                             tr_z_zzzz_yzz,   \
                             tr_z_zzzz_yzzz,  \
                             tr_z_zzzz_zzz,   \
                             tr_z_zzzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzzz_xxxx[i] = tr_z_zzzz_xxxx[i] * pa_y[i];

        tr_z_yzzzz_xxxy[i] = tr_z_zzzz_xxx[i] * fe_0 + tr_z_zzzz_xxxy[i] * pa_y[i];

        tr_z_yzzzz_xxxz[i] = tr_z_zzzz_xxxz[i] * pa_y[i];

        tr_z_yzzzz_xxyy[i] = 2.0 * tr_z_zzzz_xxy[i] * fe_0 + tr_z_zzzz_xxyy[i] * pa_y[i];

        tr_z_yzzzz_xxyz[i] = tr_z_zzzz_xxz[i] * fe_0 + tr_z_zzzz_xxyz[i] * pa_y[i];

        tr_z_yzzzz_xxzz[i] = tr_z_zzzz_xxzz[i] * pa_y[i];

        tr_z_yzzzz_xyyy[i] = 3.0 * tr_z_zzzz_xyy[i] * fe_0 + tr_z_zzzz_xyyy[i] * pa_y[i];

        tr_z_yzzzz_xyyz[i] = 2.0 * tr_z_zzzz_xyz[i] * fe_0 + tr_z_zzzz_xyyz[i] * pa_y[i];

        tr_z_yzzzz_xyzz[i] = tr_z_zzzz_xzz[i] * fe_0 + tr_z_zzzz_xyzz[i] * pa_y[i];

        tr_z_yzzzz_xzzz[i] = tr_z_zzzz_xzzz[i] * pa_y[i];

        tr_z_yzzzz_yyyy[i] = 4.0 * tr_z_zzzz_yyy[i] * fe_0 + tr_z_zzzz_yyyy[i] * pa_y[i];

        tr_z_yzzzz_yyyz[i] = 3.0 * tr_z_zzzz_yyz[i] * fe_0 + tr_z_zzzz_yyyz[i] * pa_y[i];

        tr_z_yzzzz_yyzz[i] = 2.0 * tr_z_zzzz_yzz[i] * fe_0 + tr_z_zzzz_yyzz[i] * pa_y[i];

        tr_z_yzzzz_yzzz[i] = tr_z_zzzz_zzz[i] * fe_0 + tr_z_zzzz_yzzz[i] * pa_y[i];

        tr_z_yzzzz_zzzz[i] = tr_z_zzzz_zzzz[i] * pa_y[i];
    }

    // Set up 930-945 components of targeted buffer : HG

    auto tr_z_zzzzz_xxxx = pbuffer.data(idx_dip_hg + 930);

    auto tr_z_zzzzz_xxxy = pbuffer.data(idx_dip_hg + 931);

    auto tr_z_zzzzz_xxxz = pbuffer.data(idx_dip_hg + 932);

    auto tr_z_zzzzz_xxyy = pbuffer.data(idx_dip_hg + 933);

    auto tr_z_zzzzz_xxyz = pbuffer.data(idx_dip_hg + 934);

    auto tr_z_zzzzz_xxzz = pbuffer.data(idx_dip_hg + 935);

    auto tr_z_zzzzz_xyyy = pbuffer.data(idx_dip_hg + 936);

    auto tr_z_zzzzz_xyyz = pbuffer.data(idx_dip_hg + 937);

    auto tr_z_zzzzz_xyzz = pbuffer.data(idx_dip_hg + 938);

    auto tr_z_zzzzz_xzzz = pbuffer.data(idx_dip_hg + 939);

    auto tr_z_zzzzz_yyyy = pbuffer.data(idx_dip_hg + 940);

    auto tr_z_zzzzz_yyyz = pbuffer.data(idx_dip_hg + 941);

    auto tr_z_zzzzz_yyzz = pbuffer.data(idx_dip_hg + 942);

    auto tr_z_zzzzz_yzzz = pbuffer.data(idx_dip_hg + 943);

    auto tr_z_zzzzz_zzzz = pbuffer.data(idx_dip_hg + 944);

#pragma omp simd aligned(pa_z,                \
                             tr_z_zzz_xxxx,   \
                             tr_z_zzz_xxxy,   \
                             tr_z_zzz_xxxz,   \
                             tr_z_zzz_xxyy,   \
                             tr_z_zzz_xxyz,   \
                             tr_z_zzz_xxzz,   \
                             tr_z_zzz_xyyy,   \
                             tr_z_zzz_xyyz,   \
                             tr_z_zzz_xyzz,   \
                             tr_z_zzz_xzzz,   \
                             tr_z_zzz_yyyy,   \
                             tr_z_zzz_yyyz,   \
                             tr_z_zzz_yyzz,   \
                             tr_z_zzz_yzzz,   \
                             tr_z_zzz_zzzz,   \
                             tr_z_zzzz_xxx,   \
                             tr_z_zzzz_xxxx,  \
                             tr_z_zzzz_xxxy,  \
                             tr_z_zzzz_xxxz,  \
                             tr_z_zzzz_xxy,   \
                             tr_z_zzzz_xxyy,  \
                             tr_z_zzzz_xxyz,  \
                             tr_z_zzzz_xxz,   \
                             tr_z_zzzz_xxzz,  \
                             tr_z_zzzz_xyy,   \
                             tr_z_zzzz_xyyy,  \
                             tr_z_zzzz_xyyz,  \
                             tr_z_zzzz_xyz,   \
                             tr_z_zzzz_xyzz,  \
                             tr_z_zzzz_xzz,   \
                             tr_z_zzzz_xzzz,  \
                             tr_z_zzzz_yyy,   \
                             tr_z_zzzz_yyyy,  \
                             tr_z_zzzz_yyyz,  \
                             tr_z_zzzz_yyz,   \
                             tr_z_zzzz_yyzz,  \
                             tr_z_zzzz_yzz,   \
                             tr_z_zzzz_yzzz,  \
                             tr_z_zzzz_zzz,   \
                             tr_z_zzzz_zzzz,  \
                             tr_z_zzzzz_xxxx, \
                             tr_z_zzzzz_xxxy, \
                             tr_z_zzzzz_xxxz, \
                             tr_z_zzzzz_xxyy, \
                             tr_z_zzzzz_xxyz, \
                             tr_z_zzzzz_xxzz, \
                             tr_z_zzzzz_xyyy, \
                             tr_z_zzzzz_xyyz, \
                             tr_z_zzzzz_xyzz, \
                             tr_z_zzzzz_xzzz, \
                             tr_z_zzzzz_yyyy, \
                             tr_z_zzzzz_yyyz, \
                             tr_z_zzzzz_yyzz, \
                             tr_z_zzzzz_yzzz, \
                             tr_z_zzzzz_zzzz, \
                             ts_zzzz_xxxx,    \
                             ts_zzzz_xxxy,    \
                             ts_zzzz_xxxz,    \
                             ts_zzzz_xxyy,    \
                             ts_zzzz_xxyz,    \
                             ts_zzzz_xxzz,    \
                             ts_zzzz_xyyy,    \
                             ts_zzzz_xyyz,    \
                             ts_zzzz_xyzz,    \
                             ts_zzzz_xzzz,    \
                             ts_zzzz_yyyy,    \
                             ts_zzzz_yyyz,    \
                             ts_zzzz_yyzz,    \
                             ts_zzzz_yzzz,    \
                             ts_zzzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzzz_xxxx[i] = 4.0 * tr_z_zzz_xxxx[i] * fe_0 + ts_zzzz_xxxx[i] * fe_0 + tr_z_zzzz_xxxx[i] * pa_z[i];

        tr_z_zzzzz_xxxy[i] = 4.0 * tr_z_zzz_xxxy[i] * fe_0 + ts_zzzz_xxxy[i] * fe_0 + tr_z_zzzz_xxxy[i] * pa_z[i];

        tr_z_zzzzz_xxxz[i] = 4.0 * tr_z_zzz_xxxz[i] * fe_0 + tr_z_zzzz_xxx[i] * fe_0 + ts_zzzz_xxxz[i] * fe_0 + tr_z_zzzz_xxxz[i] * pa_z[i];

        tr_z_zzzzz_xxyy[i] = 4.0 * tr_z_zzz_xxyy[i] * fe_0 + ts_zzzz_xxyy[i] * fe_0 + tr_z_zzzz_xxyy[i] * pa_z[i];

        tr_z_zzzzz_xxyz[i] = 4.0 * tr_z_zzz_xxyz[i] * fe_0 + tr_z_zzzz_xxy[i] * fe_0 + ts_zzzz_xxyz[i] * fe_0 + tr_z_zzzz_xxyz[i] * pa_z[i];

        tr_z_zzzzz_xxzz[i] = 4.0 * tr_z_zzz_xxzz[i] * fe_0 + 2.0 * tr_z_zzzz_xxz[i] * fe_0 + ts_zzzz_xxzz[i] * fe_0 + tr_z_zzzz_xxzz[i] * pa_z[i];

        tr_z_zzzzz_xyyy[i] = 4.0 * tr_z_zzz_xyyy[i] * fe_0 + ts_zzzz_xyyy[i] * fe_0 + tr_z_zzzz_xyyy[i] * pa_z[i];

        tr_z_zzzzz_xyyz[i] = 4.0 * tr_z_zzz_xyyz[i] * fe_0 + tr_z_zzzz_xyy[i] * fe_0 + ts_zzzz_xyyz[i] * fe_0 + tr_z_zzzz_xyyz[i] * pa_z[i];

        tr_z_zzzzz_xyzz[i] = 4.0 * tr_z_zzz_xyzz[i] * fe_0 + 2.0 * tr_z_zzzz_xyz[i] * fe_0 + ts_zzzz_xyzz[i] * fe_0 + tr_z_zzzz_xyzz[i] * pa_z[i];

        tr_z_zzzzz_xzzz[i] = 4.0 * tr_z_zzz_xzzz[i] * fe_0 + 3.0 * tr_z_zzzz_xzz[i] * fe_0 + ts_zzzz_xzzz[i] * fe_0 + tr_z_zzzz_xzzz[i] * pa_z[i];

        tr_z_zzzzz_yyyy[i] = 4.0 * tr_z_zzz_yyyy[i] * fe_0 + ts_zzzz_yyyy[i] * fe_0 + tr_z_zzzz_yyyy[i] * pa_z[i];

        tr_z_zzzzz_yyyz[i] = 4.0 * tr_z_zzz_yyyz[i] * fe_0 + tr_z_zzzz_yyy[i] * fe_0 + ts_zzzz_yyyz[i] * fe_0 + tr_z_zzzz_yyyz[i] * pa_z[i];

        tr_z_zzzzz_yyzz[i] = 4.0 * tr_z_zzz_yyzz[i] * fe_0 + 2.0 * tr_z_zzzz_yyz[i] * fe_0 + ts_zzzz_yyzz[i] * fe_0 + tr_z_zzzz_yyzz[i] * pa_z[i];

        tr_z_zzzzz_yzzz[i] = 4.0 * tr_z_zzz_yzzz[i] * fe_0 + 3.0 * tr_z_zzzz_yzz[i] * fe_0 + ts_zzzz_yzzz[i] * fe_0 + tr_z_zzzz_yzzz[i] * pa_z[i];

        tr_z_zzzzz_zzzz[i] = 4.0 * tr_z_zzz_zzzz[i] * fe_0 + 4.0 * tr_z_zzzz_zzz[i] * fe_0 + ts_zzzz_zzzz[i] * fe_0 + tr_z_zzzz_zzzz[i] * pa_z[i];
    }
}

}  // namespace diprec
