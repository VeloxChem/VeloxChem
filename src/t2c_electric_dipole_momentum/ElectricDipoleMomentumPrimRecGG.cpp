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

#include "ElectricDipoleMomentumPrimRecGG.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_gg(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_gg,
                                      const size_t              idx_dip_dg,
                                      const size_t              idx_dip_ff,
                                      const size_t              idx_ovl_fg,
                                      const size_t              idx_dip_fg,
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

    // Set up components of auxiliary buffer : DG

    auto tr_x_xx_xxxx = pbuffer.data(idx_dip_dg);

    auto tr_x_xx_xxxy = pbuffer.data(idx_dip_dg + 1);

    auto tr_x_xx_xxxz = pbuffer.data(idx_dip_dg + 2);

    auto tr_x_xx_xxyy = pbuffer.data(idx_dip_dg + 3);

    auto tr_x_xx_xxyz = pbuffer.data(idx_dip_dg + 4);

    auto tr_x_xx_xxzz = pbuffer.data(idx_dip_dg + 5);

    auto tr_x_xx_xyyy = pbuffer.data(idx_dip_dg + 6);

    auto tr_x_xx_xyyz = pbuffer.data(idx_dip_dg + 7);

    auto tr_x_xx_xyzz = pbuffer.data(idx_dip_dg + 8);

    auto tr_x_xx_xzzz = pbuffer.data(idx_dip_dg + 9);

    auto tr_x_xx_yyyy = pbuffer.data(idx_dip_dg + 10);

    auto tr_x_xx_yyyz = pbuffer.data(idx_dip_dg + 11);

    auto tr_x_xx_yyzz = pbuffer.data(idx_dip_dg + 12);

    auto tr_x_xx_yzzz = pbuffer.data(idx_dip_dg + 13);

    auto tr_x_xx_zzzz = pbuffer.data(idx_dip_dg + 14);

    auto tr_x_xy_xxxx = pbuffer.data(idx_dip_dg + 15);

    auto tr_x_xy_xxxz = pbuffer.data(idx_dip_dg + 17);

    auto tr_x_xy_xxzz = pbuffer.data(idx_dip_dg + 20);

    auto tr_x_xy_xzzz = pbuffer.data(idx_dip_dg + 24);

    auto tr_x_xz_xxxx = pbuffer.data(idx_dip_dg + 30);

    auto tr_x_xz_xxxy = pbuffer.data(idx_dip_dg + 31);

    auto tr_x_xz_xxxz = pbuffer.data(idx_dip_dg + 32);

    auto tr_x_xz_xxyy = pbuffer.data(idx_dip_dg + 33);

    auto tr_x_xz_xxzz = pbuffer.data(idx_dip_dg + 35);

    auto tr_x_xz_xyyy = pbuffer.data(idx_dip_dg + 36);

    auto tr_x_xz_xzzz = pbuffer.data(idx_dip_dg + 39);

    auto tr_x_yy_xxxx = pbuffer.data(idx_dip_dg + 45);

    auto tr_x_yy_xxxy = pbuffer.data(idx_dip_dg + 46);

    auto tr_x_yy_xxxz = pbuffer.data(idx_dip_dg + 47);

    auto tr_x_yy_xxyy = pbuffer.data(idx_dip_dg + 48);

    auto tr_x_yy_xxyz = pbuffer.data(idx_dip_dg + 49);

    auto tr_x_yy_xxzz = pbuffer.data(idx_dip_dg + 50);

    auto tr_x_yy_xyyy = pbuffer.data(idx_dip_dg + 51);

    auto tr_x_yy_xyyz = pbuffer.data(idx_dip_dg + 52);

    auto tr_x_yy_xyzz = pbuffer.data(idx_dip_dg + 53);

    auto tr_x_yy_xzzz = pbuffer.data(idx_dip_dg + 54);

    auto tr_x_yy_yyyy = pbuffer.data(idx_dip_dg + 55);

    auto tr_x_yy_yyyz = pbuffer.data(idx_dip_dg + 56);

    auto tr_x_yy_yyzz = pbuffer.data(idx_dip_dg + 57);

    auto tr_x_yy_yzzz = pbuffer.data(idx_dip_dg + 58);

    auto tr_x_yy_zzzz = pbuffer.data(idx_dip_dg + 59);

    auto tr_x_yz_xxxz = pbuffer.data(idx_dip_dg + 62);

    auto tr_x_yz_xxzz = pbuffer.data(idx_dip_dg + 65);

    auto tr_x_yz_xzzz = pbuffer.data(idx_dip_dg + 69);

    auto tr_x_yz_zzzz = pbuffer.data(idx_dip_dg + 74);

    auto tr_x_zz_xxxx = pbuffer.data(idx_dip_dg + 75);

    auto tr_x_zz_xxxy = pbuffer.data(idx_dip_dg + 76);

    auto tr_x_zz_xxxz = pbuffer.data(idx_dip_dg + 77);

    auto tr_x_zz_xxyy = pbuffer.data(idx_dip_dg + 78);

    auto tr_x_zz_xxyz = pbuffer.data(idx_dip_dg + 79);

    auto tr_x_zz_xxzz = pbuffer.data(idx_dip_dg + 80);

    auto tr_x_zz_xyyy = pbuffer.data(idx_dip_dg + 81);

    auto tr_x_zz_xyyz = pbuffer.data(idx_dip_dg + 82);

    auto tr_x_zz_xyzz = pbuffer.data(idx_dip_dg + 83);

    auto tr_x_zz_xzzz = pbuffer.data(idx_dip_dg + 84);

    auto tr_x_zz_yyyy = pbuffer.data(idx_dip_dg + 85);

    auto tr_x_zz_yyyz = pbuffer.data(idx_dip_dg + 86);

    auto tr_x_zz_yyzz = pbuffer.data(idx_dip_dg + 87);

    auto tr_x_zz_yzzz = pbuffer.data(idx_dip_dg + 88);

    auto tr_x_zz_zzzz = pbuffer.data(idx_dip_dg + 89);

    auto tr_y_xx_xxxx = pbuffer.data(idx_dip_dg + 90);

    auto tr_y_xx_xxxy = pbuffer.data(idx_dip_dg + 91);

    auto tr_y_xx_xxxz = pbuffer.data(idx_dip_dg + 92);

    auto tr_y_xx_xxyy = pbuffer.data(idx_dip_dg + 93);

    auto tr_y_xx_xxyz = pbuffer.data(idx_dip_dg + 94);

    auto tr_y_xx_xxzz = pbuffer.data(idx_dip_dg + 95);

    auto tr_y_xx_xyyy = pbuffer.data(idx_dip_dg + 96);

    auto tr_y_xx_xyyz = pbuffer.data(idx_dip_dg + 97);

    auto tr_y_xx_xyzz = pbuffer.data(idx_dip_dg + 98);

    auto tr_y_xx_xzzz = pbuffer.data(idx_dip_dg + 99);

    auto tr_y_xx_yyyy = pbuffer.data(idx_dip_dg + 100);

    auto tr_y_xx_yyyz = pbuffer.data(idx_dip_dg + 101);

    auto tr_y_xx_yyzz = pbuffer.data(idx_dip_dg + 102);

    auto tr_y_xx_yzzz = pbuffer.data(idx_dip_dg + 103);

    auto tr_y_xx_zzzz = pbuffer.data(idx_dip_dg + 104);

    auto tr_y_xy_xxxy = pbuffer.data(idx_dip_dg + 106);

    auto tr_y_xy_xxyy = pbuffer.data(idx_dip_dg + 108);

    auto tr_y_xy_xxyz = pbuffer.data(idx_dip_dg + 109);

    auto tr_y_xy_xyyy = pbuffer.data(idx_dip_dg + 111);

    auto tr_y_xy_xyyz = pbuffer.data(idx_dip_dg + 112);

    auto tr_y_xy_xyzz = pbuffer.data(idx_dip_dg + 113);

    auto tr_y_xy_yyyy = pbuffer.data(idx_dip_dg + 115);

    auto tr_y_xy_yyyz = pbuffer.data(idx_dip_dg + 116);

    auto tr_y_xy_yyzz = pbuffer.data(idx_dip_dg + 117);

    auto tr_y_xy_yzzz = pbuffer.data(idx_dip_dg + 118);

    auto tr_y_xy_zzzz = pbuffer.data(idx_dip_dg + 119);

    auto tr_y_xz_yyyz = pbuffer.data(idx_dip_dg + 131);

    auto tr_y_xz_yyzz = pbuffer.data(idx_dip_dg + 132);

    auto tr_y_xz_yzzz = pbuffer.data(idx_dip_dg + 133);

    auto tr_y_xz_zzzz = pbuffer.data(idx_dip_dg + 134);

    auto tr_y_yy_xxxx = pbuffer.data(idx_dip_dg + 135);

    auto tr_y_yy_xxxy = pbuffer.data(idx_dip_dg + 136);

    auto tr_y_yy_xxxz = pbuffer.data(idx_dip_dg + 137);

    auto tr_y_yy_xxyy = pbuffer.data(idx_dip_dg + 138);

    auto tr_y_yy_xxyz = pbuffer.data(idx_dip_dg + 139);

    auto tr_y_yy_xxzz = pbuffer.data(idx_dip_dg + 140);

    auto tr_y_yy_xyyy = pbuffer.data(idx_dip_dg + 141);

    auto tr_y_yy_xyyz = pbuffer.data(idx_dip_dg + 142);

    auto tr_y_yy_xyzz = pbuffer.data(idx_dip_dg + 143);

    auto tr_y_yy_xzzz = pbuffer.data(idx_dip_dg + 144);

    auto tr_y_yy_yyyy = pbuffer.data(idx_dip_dg + 145);

    auto tr_y_yy_yyyz = pbuffer.data(idx_dip_dg + 146);

    auto tr_y_yy_yyzz = pbuffer.data(idx_dip_dg + 147);

    auto tr_y_yy_yzzz = pbuffer.data(idx_dip_dg + 148);

    auto tr_y_yy_zzzz = pbuffer.data(idx_dip_dg + 149);

    auto tr_y_yz_xxxy = pbuffer.data(idx_dip_dg + 151);

    auto tr_y_yz_xxyy = pbuffer.data(idx_dip_dg + 153);

    auto tr_y_yz_xyyy = pbuffer.data(idx_dip_dg + 156);

    auto tr_y_yz_yyyy = pbuffer.data(idx_dip_dg + 160);

    auto tr_y_yz_yyyz = pbuffer.data(idx_dip_dg + 161);

    auto tr_y_yz_yyzz = pbuffer.data(idx_dip_dg + 162);

    auto tr_y_yz_yzzz = pbuffer.data(idx_dip_dg + 163);

    auto tr_y_yz_zzzz = pbuffer.data(idx_dip_dg + 164);

    auto tr_y_zz_xxxx = pbuffer.data(idx_dip_dg + 165);

    auto tr_y_zz_xxxy = pbuffer.data(idx_dip_dg + 166);

    auto tr_y_zz_xxxz = pbuffer.data(idx_dip_dg + 167);

    auto tr_y_zz_xxyy = pbuffer.data(idx_dip_dg + 168);

    auto tr_y_zz_xxyz = pbuffer.data(idx_dip_dg + 169);

    auto tr_y_zz_xxzz = pbuffer.data(idx_dip_dg + 170);

    auto tr_y_zz_xyyy = pbuffer.data(idx_dip_dg + 171);

    auto tr_y_zz_xyyz = pbuffer.data(idx_dip_dg + 172);

    auto tr_y_zz_xyzz = pbuffer.data(idx_dip_dg + 173);

    auto tr_y_zz_xzzz = pbuffer.data(idx_dip_dg + 174);

    auto tr_y_zz_yyyy = pbuffer.data(idx_dip_dg + 175);

    auto tr_y_zz_yyyz = pbuffer.data(idx_dip_dg + 176);

    auto tr_y_zz_yyzz = pbuffer.data(idx_dip_dg + 177);

    auto tr_y_zz_yzzz = pbuffer.data(idx_dip_dg + 178);

    auto tr_y_zz_zzzz = pbuffer.data(idx_dip_dg + 179);

    auto tr_z_xx_xxxx = pbuffer.data(idx_dip_dg + 180);

    auto tr_z_xx_xxxy = pbuffer.data(idx_dip_dg + 181);

    auto tr_z_xx_xxxz = pbuffer.data(idx_dip_dg + 182);

    auto tr_z_xx_xxyy = pbuffer.data(idx_dip_dg + 183);

    auto tr_z_xx_xxyz = pbuffer.data(idx_dip_dg + 184);

    auto tr_z_xx_xxzz = pbuffer.data(idx_dip_dg + 185);

    auto tr_z_xx_xyyy = pbuffer.data(idx_dip_dg + 186);

    auto tr_z_xx_xyyz = pbuffer.data(idx_dip_dg + 187);

    auto tr_z_xx_xyzz = pbuffer.data(idx_dip_dg + 188);

    auto tr_z_xx_xzzz = pbuffer.data(idx_dip_dg + 189);

    auto tr_z_xx_yyyy = pbuffer.data(idx_dip_dg + 190);

    auto tr_z_xx_yyyz = pbuffer.data(idx_dip_dg + 191);

    auto tr_z_xx_yyzz = pbuffer.data(idx_dip_dg + 192);

    auto tr_z_xx_yzzz = pbuffer.data(idx_dip_dg + 193);

    auto tr_z_xx_zzzz = pbuffer.data(idx_dip_dg + 194);

    auto tr_z_xy_yyyy = pbuffer.data(idx_dip_dg + 205);

    auto tr_z_xy_yyyz = pbuffer.data(idx_dip_dg + 206);

    auto tr_z_xy_yyzz = pbuffer.data(idx_dip_dg + 207);

    auto tr_z_xy_yzzz = pbuffer.data(idx_dip_dg + 208);

    auto tr_z_xz_xxxz = pbuffer.data(idx_dip_dg + 212);

    auto tr_z_xz_xxyz = pbuffer.data(idx_dip_dg + 214);

    auto tr_z_xz_xxzz = pbuffer.data(idx_dip_dg + 215);

    auto tr_z_xz_xyyz = pbuffer.data(idx_dip_dg + 217);

    auto tr_z_xz_xyzz = pbuffer.data(idx_dip_dg + 218);

    auto tr_z_xz_xzzz = pbuffer.data(idx_dip_dg + 219);

    auto tr_z_xz_yyyy = pbuffer.data(idx_dip_dg + 220);

    auto tr_z_xz_yyyz = pbuffer.data(idx_dip_dg + 221);

    auto tr_z_xz_yyzz = pbuffer.data(idx_dip_dg + 222);

    auto tr_z_xz_yzzz = pbuffer.data(idx_dip_dg + 223);

    auto tr_z_xz_zzzz = pbuffer.data(idx_dip_dg + 224);

    auto tr_z_yy_xxxx = pbuffer.data(idx_dip_dg + 225);

    auto tr_z_yy_xxxy = pbuffer.data(idx_dip_dg + 226);

    auto tr_z_yy_xxxz = pbuffer.data(idx_dip_dg + 227);

    auto tr_z_yy_xxyy = pbuffer.data(idx_dip_dg + 228);

    auto tr_z_yy_xxyz = pbuffer.data(idx_dip_dg + 229);

    auto tr_z_yy_xxzz = pbuffer.data(idx_dip_dg + 230);

    auto tr_z_yy_xyyy = pbuffer.data(idx_dip_dg + 231);

    auto tr_z_yy_xyyz = pbuffer.data(idx_dip_dg + 232);

    auto tr_z_yy_xyzz = pbuffer.data(idx_dip_dg + 233);

    auto tr_z_yy_xzzz = pbuffer.data(idx_dip_dg + 234);

    auto tr_z_yy_yyyy = pbuffer.data(idx_dip_dg + 235);

    auto tr_z_yy_yyyz = pbuffer.data(idx_dip_dg + 236);

    auto tr_z_yy_yyzz = pbuffer.data(idx_dip_dg + 237);

    auto tr_z_yy_yzzz = pbuffer.data(idx_dip_dg + 238);

    auto tr_z_yy_zzzz = pbuffer.data(idx_dip_dg + 239);

    auto tr_z_yz_xxxx = pbuffer.data(idx_dip_dg + 240);

    auto tr_z_yz_xxxz = pbuffer.data(idx_dip_dg + 242);

    auto tr_z_yz_xxyz = pbuffer.data(idx_dip_dg + 244);

    auto tr_z_yz_xxzz = pbuffer.data(idx_dip_dg + 245);

    auto tr_z_yz_xyyz = pbuffer.data(idx_dip_dg + 247);

    auto tr_z_yz_xyzz = pbuffer.data(idx_dip_dg + 248);

    auto tr_z_yz_xzzz = pbuffer.data(idx_dip_dg + 249);

    auto tr_z_yz_yyyy = pbuffer.data(idx_dip_dg + 250);

    auto tr_z_yz_yyyz = pbuffer.data(idx_dip_dg + 251);

    auto tr_z_yz_yyzz = pbuffer.data(idx_dip_dg + 252);

    auto tr_z_yz_yzzz = pbuffer.data(idx_dip_dg + 253);

    auto tr_z_yz_zzzz = pbuffer.data(idx_dip_dg + 254);

    auto tr_z_zz_xxxx = pbuffer.data(idx_dip_dg + 255);

    auto tr_z_zz_xxxy = pbuffer.data(idx_dip_dg + 256);

    auto tr_z_zz_xxxz = pbuffer.data(idx_dip_dg + 257);

    auto tr_z_zz_xxyy = pbuffer.data(idx_dip_dg + 258);

    auto tr_z_zz_xxyz = pbuffer.data(idx_dip_dg + 259);

    auto tr_z_zz_xxzz = pbuffer.data(idx_dip_dg + 260);

    auto tr_z_zz_xyyy = pbuffer.data(idx_dip_dg + 261);

    auto tr_z_zz_xyyz = pbuffer.data(idx_dip_dg + 262);

    auto tr_z_zz_xyzz = pbuffer.data(idx_dip_dg + 263);

    auto tr_z_zz_xzzz = pbuffer.data(idx_dip_dg + 264);

    auto tr_z_zz_yyyy = pbuffer.data(idx_dip_dg + 265);

    auto tr_z_zz_yyyz = pbuffer.data(idx_dip_dg + 266);

    auto tr_z_zz_yyzz = pbuffer.data(idx_dip_dg + 267);

    auto tr_z_zz_yzzz = pbuffer.data(idx_dip_dg + 268);

    auto tr_z_zz_zzzz = pbuffer.data(idx_dip_dg + 269);

    // Set up components of auxiliary buffer : FF

    auto tr_x_xxx_xxx = pbuffer.data(idx_dip_ff);

    auto tr_x_xxx_xxy = pbuffer.data(idx_dip_ff + 1);

    auto tr_x_xxx_xxz = pbuffer.data(idx_dip_ff + 2);

    auto tr_x_xxx_xyy = pbuffer.data(idx_dip_ff + 3);

    auto tr_x_xxx_xyz = pbuffer.data(idx_dip_ff + 4);

    auto tr_x_xxx_xzz = pbuffer.data(idx_dip_ff + 5);

    auto tr_x_xxx_yyy = pbuffer.data(idx_dip_ff + 6);

    auto tr_x_xxx_yyz = pbuffer.data(idx_dip_ff + 7);

    auto tr_x_xxx_yzz = pbuffer.data(idx_dip_ff + 8);

    auto tr_x_xxx_zzz = pbuffer.data(idx_dip_ff + 9);

    auto tr_x_xxy_xxx = pbuffer.data(idx_dip_ff + 10);

    auto tr_x_xxy_xxy = pbuffer.data(idx_dip_ff + 11);

    auto tr_x_xxy_xxz = pbuffer.data(idx_dip_ff + 12);

    auto tr_x_xxy_xyy = pbuffer.data(idx_dip_ff + 13);

    auto tr_x_xxy_xyz = pbuffer.data(idx_dip_ff + 14);

    auto tr_x_xxy_xzz = pbuffer.data(idx_dip_ff + 15);

    auto tr_x_xxz_xxx = pbuffer.data(idx_dip_ff + 20);

    auto tr_x_xxz_xxy = pbuffer.data(idx_dip_ff + 21);

    auto tr_x_xxz_xxz = pbuffer.data(idx_dip_ff + 22);

    auto tr_x_xxz_xyy = pbuffer.data(idx_dip_ff + 23);

    auto tr_x_xxz_xyz = pbuffer.data(idx_dip_ff + 24);

    auto tr_x_xxz_xzz = pbuffer.data(idx_dip_ff + 25);

    auto tr_x_xxz_yyz = pbuffer.data(idx_dip_ff + 27);

    auto tr_x_xxz_yzz = pbuffer.data(idx_dip_ff + 28);

    auto tr_x_xxz_zzz = pbuffer.data(idx_dip_ff + 29);

    auto tr_x_xyy_xxy = pbuffer.data(idx_dip_ff + 31);

    auto tr_x_xyy_xyy = pbuffer.data(idx_dip_ff + 33);

    auto tr_x_xyy_xyz = pbuffer.data(idx_dip_ff + 34);

    auto tr_x_xzz_xxx = pbuffer.data(idx_dip_ff + 50);

    auto tr_x_xzz_xxy = pbuffer.data(idx_dip_ff + 51);

    auto tr_x_xzz_xxz = pbuffer.data(idx_dip_ff + 52);

    auto tr_x_xzz_xyy = pbuffer.data(idx_dip_ff + 53);

    auto tr_x_xzz_xyz = pbuffer.data(idx_dip_ff + 54);

    auto tr_x_xzz_xzz = pbuffer.data(idx_dip_ff + 55);

    auto tr_x_yyy_xxx = pbuffer.data(idx_dip_ff + 60);

    auto tr_x_yyy_xxy = pbuffer.data(idx_dip_ff + 61);

    auto tr_x_yyy_xxz = pbuffer.data(idx_dip_ff + 62);

    auto tr_x_yyy_xyy = pbuffer.data(idx_dip_ff + 63);

    auto tr_x_yyy_xyz = pbuffer.data(idx_dip_ff + 64);

    auto tr_x_yyy_xzz = pbuffer.data(idx_dip_ff + 65);

    auto tr_x_yyy_yyy = pbuffer.data(idx_dip_ff + 66);

    auto tr_x_yyy_yyz = pbuffer.data(idx_dip_ff + 67);

    auto tr_x_yyy_yzz = pbuffer.data(idx_dip_ff + 68);

    auto tr_x_yyy_zzz = pbuffer.data(idx_dip_ff + 69);

    auto tr_x_yzz_xxz = pbuffer.data(idx_dip_ff + 82);

    auto tr_x_yzz_xyz = pbuffer.data(idx_dip_ff + 84);

    auto tr_x_yzz_xzz = pbuffer.data(idx_dip_ff + 85);

    auto tr_x_yzz_yyz = pbuffer.data(idx_dip_ff + 87);

    auto tr_x_yzz_yzz = pbuffer.data(idx_dip_ff + 88);

    auto tr_x_yzz_zzz = pbuffer.data(idx_dip_ff + 89);

    auto tr_x_zzz_xxx = pbuffer.data(idx_dip_ff + 90);

    auto tr_x_zzz_xxy = pbuffer.data(idx_dip_ff + 91);

    auto tr_x_zzz_xxz = pbuffer.data(idx_dip_ff + 92);

    auto tr_x_zzz_xyy = pbuffer.data(idx_dip_ff + 93);

    auto tr_x_zzz_xyz = pbuffer.data(idx_dip_ff + 94);

    auto tr_x_zzz_xzz = pbuffer.data(idx_dip_ff + 95);

    auto tr_x_zzz_yyy = pbuffer.data(idx_dip_ff + 96);

    auto tr_x_zzz_yyz = pbuffer.data(idx_dip_ff + 97);

    auto tr_x_zzz_yzz = pbuffer.data(idx_dip_ff + 98);

    auto tr_x_zzz_zzz = pbuffer.data(idx_dip_ff + 99);

    auto tr_y_xxx_xxx = pbuffer.data(idx_dip_ff + 100);

    auto tr_y_xxx_xxy = pbuffer.data(idx_dip_ff + 101);

    auto tr_y_xxx_xxz = pbuffer.data(idx_dip_ff + 102);

    auto tr_y_xxx_xyy = pbuffer.data(idx_dip_ff + 103);

    auto tr_y_xxx_xyz = pbuffer.data(idx_dip_ff + 104);

    auto tr_y_xxx_xzz = pbuffer.data(idx_dip_ff + 105);

    auto tr_y_xxx_yyy = pbuffer.data(idx_dip_ff + 106);

    auto tr_y_xxx_yyz = pbuffer.data(idx_dip_ff + 107);

    auto tr_y_xxx_yzz = pbuffer.data(idx_dip_ff + 108);

    auto tr_y_xxx_zzz = pbuffer.data(idx_dip_ff + 109);

    auto tr_y_xxy_xxy = pbuffer.data(idx_dip_ff + 111);

    auto tr_y_xxy_xyy = pbuffer.data(idx_dip_ff + 113);

    auto tr_y_xxy_xyz = pbuffer.data(idx_dip_ff + 114);

    auto tr_y_xxy_yyy = pbuffer.data(idx_dip_ff + 116);

    auto tr_y_xxy_yyz = pbuffer.data(idx_dip_ff + 117);

    auto tr_y_xxy_yzz = pbuffer.data(idx_dip_ff + 118);

    auto tr_y_xyy_xxx = pbuffer.data(idx_dip_ff + 130);

    auto tr_y_xyy_xxy = pbuffer.data(idx_dip_ff + 131);

    auto tr_y_xyy_xxz = pbuffer.data(idx_dip_ff + 132);

    auto tr_y_xyy_xyy = pbuffer.data(idx_dip_ff + 133);

    auto tr_y_xyy_xyz = pbuffer.data(idx_dip_ff + 134);

    auto tr_y_xyy_xzz = pbuffer.data(idx_dip_ff + 135);

    auto tr_y_xyy_yyy = pbuffer.data(idx_dip_ff + 136);

    auto tr_y_xyy_yyz = pbuffer.data(idx_dip_ff + 137);

    auto tr_y_xyy_yzz = pbuffer.data(idx_dip_ff + 138);

    auto tr_y_xyy_zzz = pbuffer.data(idx_dip_ff + 139);

    auto tr_y_xzz_xxz = pbuffer.data(idx_dip_ff + 152);

    auto tr_y_xzz_xyz = pbuffer.data(idx_dip_ff + 154);

    auto tr_y_xzz_xzz = pbuffer.data(idx_dip_ff + 155);

    auto tr_y_xzz_yyz = pbuffer.data(idx_dip_ff + 157);

    auto tr_y_xzz_yzz = pbuffer.data(idx_dip_ff + 158);

    auto tr_y_xzz_zzz = pbuffer.data(idx_dip_ff + 159);

    auto tr_y_yyy_xxx = pbuffer.data(idx_dip_ff + 160);

    auto tr_y_yyy_xxy = pbuffer.data(idx_dip_ff + 161);

    auto tr_y_yyy_xxz = pbuffer.data(idx_dip_ff + 162);

    auto tr_y_yyy_xyy = pbuffer.data(idx_dip_ff + 163);

    auto tr_y_yyy_xyz = pbuffer.data(idx_dip_ff + 164);

    auto tr_y_yyy_xzz = pbuffer.data(idx_dip_ff + 165);

    auto tr_y_yyy_yyy = pbuffer.data(idx_dip_ff + 166);

    auto tr_y_yyy_yyz = pbuffer.data(idx_dip_ff + 167);

    auto tr_y_yyy_yzz = pbuffer.data(idx_dip_ff + 168);

    auto tr_y_yyy_zzz = pbuffer.data(idx_dip_ff + 169);

    auto tr_y_yyz_xxy = pbuffer.data(idx_dip_ff + 171);

    auto tr_y_yyz_xxz = pbuffer.data(idx_dip_ff + 172);

    auto tr_y_yyz_xyy = pbuffer.data(idx_dip_ff + 173);

    auto tr_y_yyz_xyz = pbuffer.data(idx_dip_ff + 174);

    auto tr_y_yyz_xzz = pbuffer.data(idx_dip_ff + 175);

    auto tr_y_yyz_yyy = pbuffer.data(idx_dip_ff + 176);

    auto tr_y_yyz_yyz = pbuffer.data(idx_dip_ff + 177);

    auto tr_y_yyz_yzz = pbuffer.data(idx_dip_ff + 178);

    auto tr_y_yyz_zzz = pbuffer.data(idx_dip_ff + 179);

    auto tr_y_yzz_xxx = pbuffer.data(idx_dip_ff + 180);

    auto tr_y_yzz_xxy = pbuffer.data(idx_dip_ff + 181);

    auto tr_y_yzz_xxz = pbuffer.data(idx_dip_ff + 182);

    auto tr_y_yzz_xyy = pbuffer.data(idx_dip_ff + 183);

    auto tr_y_yzz_xyz = pbuffer.data(idx_dip_ff + 184);

    auto tr_y_yzz_xzz = pbuffer.data(idx_dip_ff + 185);

    auto tr_y_yzz_yyy = pbuffer.data(idx_dip_ff + 186);

    auto tr_y_yzz_yyz = pbuffer.data(idx_dip_ff + 187);

    auto tr_y_yzz_yzz = pbuffer.data(idx_dip_ff + 188);

    auto tr_y_yzz_zzz = pbuffer.data(idx_dip_ff + 189);

    auto tr_y_zzz_xxx = pbuffer.data(idx_dip_ff + 190);

    auto tr_y_zzz_xxy = pbuffer.data(idx_dip_ff + 191);

    auto tr_y_zzz_xxz = pbuffer.data(idx_dip_ff + 192);

    auto tr_y_zzz_xyy = pbuffer.data(idx_dip_ff + 193);

    auto tr_y_zzz_xyz = pbuffer.data(idx_dip_ff + 194);

    auto tr_y_zzz_xzz = pbuffer.data(idx_dip_ff + 195);

    auto tr_y_zzz_yyy = pbuffer.data(idx_dip_ff + 196);

    auto tr_y_zzz_yyz = pbuffer.data(idx_dip_ff + 197);

    auto tr_y_zzz_yzz = pbuffer.data(idx_dip_ff + 198);

    auto tr_y_zzz_zzz = pbuffer.data(idx_dip_ff + 199);

    auto tr_z_xxx_xxx = pbuffer.data(idx_dip_ff + 200);

    auto tr_z_xxx_xxy = pbuffer.data(idx_dip_ff + 201);

    auto tr_z_xxx_xxz = pbuffer.data(idx_dip_ff + 202);

    auto tr_z_xxx_xyy = pbuffer.data(idx_dip_ff + 203);

    auto tr_z_xxx_xyz = pbuffer.data(idx_dip_ff + 204);

    auto tr_z_xxx_xzz = pbuffer.data(idx_dip_ff + 205);

    auto tr_z_xxx_yyy = pbuffer.data(idx_dip_ff + 206);

    auto tr_z_xxx_yyz = pbuffer.data(idx_dip_ff + 207);

    auto tr_z_xxx_yzz = pbuffer.data(idx_dip_ff + 208);

    auto tr_z_xxx_zzz = pbuffer.data(idx_dip_ff + 209);

    auto tr_z_xxz_xxx = pbuffer.data(idx_dip_ff + 220);

    auto tr_z_xxz_xxy = pbuffer.data(idx_dip_ff + 221);

    auto tr_z_xxz_xxz = pbuffer.data(idx_dip_ff + 222);

    auto tr_z_xxz_xyy = pbuffer.data(idx_dip_ff + 223);

    auto tr_z_xxz_xyz = pbuffer.data(idx_dip_ff + 224);

    auto tr_z_xxz_xzz = pbuffer.data(idx_dip_ff + 225);

    auto tr_z_xxz_yyz = pbuffer.data(idx_dip_ff + 227);

    auto tr_z_xxz_yzz = pbuffer.data(idx_dip_ff + 228);

    auto tr_z_xxz_zzz = pbuffer.data(idx_dip_ff + 229);

    auto tr_z_xyy_xxy = pbuffer.data(idx_dip_ff + 231);

    auto tr_z_xyy_xyy = pbuffer.data(idx_dip_ff + 233);

    auto tr_z_xyy_xyz = pbuffer.data(idx_dip_ff + 234);

    auto tr_z_xyy_yyy = pbuffer.data(idx_dip_ff + 236);

    auto tr_z_xyy_yyz = pbuffer.data(idx_dip_ff + 237);

    auto tr_z_xyy_yzz = pbuffer.data(idx_dip_ff + 238);

    auto tr_z_xzz_xxx = pbuffer.data(idx_dip_ff + 250);

    auto tr_z_xzz_xxy = pbuffer.data(idx_dip_ff + 251);

    auto tr_z_xzz_xxz = pbuffer.data(idx_dip_ff + 252);

    auto tr_z_xzz_xyy = pbuffer.data(idx_dip_ff + 253);

    auto tr_z_xzz_xyz = pbuffer.data(idx_dip_ff + 254);

    auto tr_z_xzz_xzz = pbuffer.data(idx_dip_ff + 255);

    auto tr_z_xzz_yyy = pbuffer.data(idx_dip_ff + 256);

    auto tr_z_xzz_yyz = pbuffer.data(idx_dip_ff + 257);

    auto tr_z_xzz_yzz = pbuffer.data(idx_dip_ff + 258);

    auto tr_z_xzz_zzz = pbuffer.data(idx_dip_ff + 259);

    auto tr_z_yyy_xxx = pbuffer.data(idx_dip_ff + 260);

    auto tr_z_yyy_xxy = pbuffer.data(idx_dip_ff + 261);

    auto tr_z_yyy_xxz = pbuffer.data(idx_dip_ff + 262);

    auto tr_z_yyy_xyy = pbuffer.data(idx_dip_ff + 263);

    auto tr_z_yyy_xyz = pbuffer.data(idx_dip_ff + 264);

    auto tr_z_yyy_xzz = pbuffer.data(idx_dip_ff + 265);

    auto tr_z_yyy_yyy = pbuffer.data(idx_dip_ff + 266);

    auto tr_z_yyy_yyz = pbuffer.data(idx_dip_ff + 267);

    auto tr_z_yyy_yzz = pbuffer.data(idx_dip_ff + 268);

    auto tr_z_yyy_zzz = pbuffer.data(idx_dip_ff + 269);

    auto tr_z_yyz_xxx = pbuffer.data(idx_dip_ff + 270);

    auto tr_z_yyz_xxy = pbuffer.data(idx_dip_ff + 271);

    auto tr_z_yyz_xxz = pbuffer.data(idx_dip_ff + 272);

    auto tr_z_yyz_xyy = pbuffer.data(idx_dip_ff + 273);

    auto tr_z_yyz_xyz = pbuffer.data(idx_dip_ff + 274);

    auto tr_z_yyz_xzz = pbuffer.data(idx_dip_ff + 275);

    auto tr_z_yyz_yyy = pbuffer.data(idx_dip_ff + 276);

    auto tr_z_yyz_yyz = pbuffer.data(idx_dip_ff + 277);

    auto tr_z_yyz_yzz = pbuffer.data(idx_dip_ff + 278);

    auto tr_z_yyz_zzz = pbuffer.data(idx_dip_ff + 279);

    auto tr_z_yzz_xxx = pbuffer.data(idx_dip_ff + 280);

    auto tr_z_yzz_xxy = pbuffer.data(idx_dip_ff + 281);

    auto tr_z_yzz_xxz = pbuffer.data(idx_dip_ff + 282);

    auto tr_z_yzz_xyy = pbuffer.data(idx_dip_ff + 283);

    auto tr_z_yzz_xyz = pbuffer.data(idx_dip_ff + 284);

    auto tr_z_yzz_xzz = pbuffer.data(idx_dip_ff + 285);

    auto tr_z_yzz_yyy = pbuffer.data(idx_dip_ff + 286);

    auto tr_z_yzz_yyz = pbuffer.data(idx_dip_ff + 287);

    auto tr_z_yzz_yzz = pbuffer.data(idx_dip_ff + 288);

    auto tr_z_yzz_zzz = pbuffer.data(idx_dip_ff + 289);

    auto tr_z_zzz_xxx = pbuffer.data(idx_dip_ff + 290);

    auto tr_z_zzz_xxy = pbuffer.data(idx_dip_ff + 291);

    auto tr_z_zzz_xxz = pbuffer.data(idx_dip_ff + 292);

    auto tr_z_zzz_xyy = pbuffer.data(idx_dip_ff + 293);

    auto tr_z_zzz_xyz = pbuffer.data(idx_dip_ff + 294);

    auto tr_z_zzz_xzz = pbuffer.data(idx_dip_ff + 295);

    auto tr_z_zzz_yyy = pbuffer.data(idx_dip_ff + 296);

    auto tr_z_zzz_yyz = pbuffer.data(idx_dip_ff + 297);

    auto tr_z_zzz_yzz = pbuffer.data(idx_dip_ff + 298);

    auto tr_z_zzz_zzz = pbuffer.data(idx_dip_ff + 299);

    // Set up components of auxiliary buffer : FG

    auto ts_xxx_xxxx = pbuffer.data(idx_ovl_fg);

    auto ts_xxx_xxxy = pbuffer.data(idx_ovl_fg + 1);

    auto ts_xxx_xxxz = pbuffer.data(idx_ovl_fg + 2);

    auto ts_xxx_xxyy = pbuffer.data(idx_ovl_fg + 3);

    auto ts_xxx_xxyz = pbuffer.data(idx_ovl_fg + 4);

    auto ts_xxx_xxzz = pbuffer.data(idx_ovl_fg + 5);

    auto ts_xxx_xyyy = pbuffer.data(idx_ovl_fg + 6);

    auto ts_xxx_xyyz = pbuffer.data(idx_ovl_fg + 7);

    auto ts_xxx_xyzz = pbuffer.data(idx_ovl_fg + 8);

    auto ts_xxx_xzzz = pbuffer.data(idx_ovl_fg + 9);

    auto ts_xxx_yyyy = pbuffer.data(idx_ovl_fg + 10);

    auto ts_xxx_yyyz = pbuffer.data(idx_ovl_fg + 11);

    auto ts_xxx_yyzz = pbuffer.data(idx_ovl_fg + 12);

    auto ts_xxx_yzzz = pbuffer.data(idx_ovl_fg + 13);

    auto ts_xxx_zzzz = pbuffer.data(idx_ovl_fg + 14);

    auto ts_xxz_xxxz = pbuffer.data(idx_ovl_fg + 32);

    auto ts_xxz_xxzz = pbuffer.data(idx_ovl_fg + 35);

    auto ts_xxz_xzzz = pbuffer.data(idx_ovl_fg + 39);

    auto ts_xyy_yyyy = pbuffer.data(idx_ovl_fg + 55);

    auto ts_xyy_yyyz = pbuffer.data(idx_ovl_fg + 56);

    auto ts_xyy_yyzz = pbuffer.data(idx_ovl_fg + 57);

    auto ts_xyy_yzzz = pbuffer.data(idx_ovl_fg + 58);

    auto ts_xzz_yyyz = pbuffer.data(idx_ovl_fg + 86);

    auto ts_xzz_yyzz = pbuffer.data(idx_ovl_fg + 87);

    auto ts_xzz_yzzz = pbuffer.data(idx_ovl_fg + 88);

    auto ts_xzz_zzzz = pbuffer.data(idx_ovl_fg + 89);

    auto ts_yyy_xxxx = pbuffer.data(idx_ovl_fg + 90);

    auto ts_yyy_xxxy = pbuffer.data(idx_ovl_fg + 91);

    auto ts_yyy_xxxz = pbuffer.data(idx_ovl_fg + 92);

    auto ts_yyy_xxyy = pbuffer.data(idx_ovl_fg + 93);

    auto ts_yyy_xxyz = pbuffer.data(idx_ovl_fg + 94);

    auto ts_yyy_xxzz = pbuffer.data(idx_ovl_fg + 95);

    auto ts_yyy_xyyy = pbuffer.data(idx_ovl_fg + 96);

    auto ts_yyy_xyyz = pbuffer.data(idx_ovl_fg + 97);

    auto ts_yyy_xyzz = pbuffer.data(idx_ovl_fg + 98);

    auto ts_yyy_xzzz = pbuffer.data(idx_ovl_fg + 99);

    auto ts_yyy_yyyy = pbuffer.data(idx_ovl_fg + 100);

    auto ts_yyy_yyyz = pbuffer.data(idx_ovl_fg + 101);

    auto ts_yyy_yyzz = pbuffer.data(idx_ovl_fg + 102);

    auto ts_yyy_yzzz = pbuffer.data(idx_ovl_fg + 103);

    auto ts_yyy_zzzz = pbuffer.data(idx_ovl_fg + 104);

    auto ts_yyz_yyyz = pbuffer.data(idx_ovl_fg + 116);

    auto ts_yyz_yyzz = pbuffer.data(idx_ovl_fg + 117);

    auto ts_yyz_yzzz = pbuffer.data(idx_ovl_fg + 118);

    auto ts_yyz_zzzz = pbuffer.data(idx_ovl_fg + 119);

    auto ts_yzz_xxxz = pbuffer.data(idx_ovl_fg + 122);

    auto ts_yzz_xxzz = pbuffer.data(idx_ovl_fg + 125);

    auto ts_yzz_xzzz = pbuffer.data(idx_ovl_fg + 129);

    auto ts_yzz_yyyy = pbuffer.data(idx_ovl_fg + 130);

    auto ts_yzz_yyyz = pbuffer.data(idx_ovl_fg + 131);

    auto ts_yzz_yyzz = pbuffer.data(idx_ovl_fg + 132);

    auto ts_yzz_yzzz = pbuffer.data(idx_ovl_fg + 133);

    auto ts_yzz_zzzz = pbuffer.data(idx_ovl_fg + 134);

    auto ts_zzz_xxxx = pbuffer.data(idx_ovl_fg + 135);

    auto ts_zzz_xxxy = pbuffer.data(idx_ovl_fg + 136);

    auto ts_zzz_xxxz = pbuffer.data(idx_ovl_fg + 137);

    auto ts_zzz_xxyy = pbuffer.data(idx_ovl_fg + 138);

    auto ts_zzz_xxyz = pbuffer.data(idx_ovl_fg + 139);

    auto ts_zzz_xxzz = pbuffer.data(idx_ovl_fg + 140);

    auto ts_zzz_xyyy = pbuffer.data(idx_ovl_fg + 141);

    auto ts_zzz_xyyz = pbuffer.data(idx_ovl_fg + 142);

    auto ts_zzz_xyzz = pbuffer.data(idx_ovl_fg + 143);

    auto ts_zzz_xzzz = pbuffer.data(idx_ovl_fg + 144);

    auto ts_zzz_yyyy = pbuffer.data(idx_ovl_fg + 145);

    auto ts_zzz_yyyz = pbuffer.data(idx_ovl_fg + 146);

    auto ts_zzz_yyzz = pbuffer.data(idx_ovl_fg + 147);

    auto ts_zzz_yzzz = pbuffer.data(idx_ovl_fg + 148);

    auto ts_zzz_zzzz = pbuffer.data(idx_ovl_fg + 149);

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

    auto tr_x_xxy_yyyy = pbuffer.data(idx_dip_fg + 25);

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

    auto tr_x_xxz_yyyz = pbuffer.data(idx_dip_fg + 41);

    auto tr_x_xxz_yyzz = pbuffer.data(idx_dip_fg + 42);

    auto tr_x_xxz_yzzz = pbuffer.data(idx_dip_fg + 43);

    auto tr_x_xxz_zzzz = pbuffer.data(idx_dip_fg + 44);

    auto tr_x_xyy_xxxx = pbuffer.data(idx_dip_fg + 45);

    auto tr_x_xyy_xxxy = pbuffer.data(idx_dip_fg + 46);

    auto tr_x_xyy_xxxz = pbuffer.data(idx_dip_fg + 47);

    auto tr_x_xyy_xxyy = pbuffer.data(idx_dip_fg + 48);

    auto tr_x_xyy_xxyz = pbuffer.data(idx_dip_fg + 49);

    auto tr_x_xyy_xxzz = pbuffer.data(idx_dip_fg + 50);

    auto tr_x_xyy_xyyy = pbuffer.data(idx_dip_fg + 51);

    auto tr_x_xyy_xyyz = pbuffer.data(idx_dip_fg + 52);

    auto tr_x_xyy_xyzz = pbuffer.data(idx_dip_fg + 53);

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

    auto tr_x_xzz_xxyz = pbuffer.data(idx_dip_fg + 79);

    auto tr_x_xzz_xxzz = pbuffer.data(idx_dip_fg + 80);

    auto tr_x_xzz_xyyy = pbuffer.data(idx_dip_fg + 81);

    auto tr_x_xzz_xyyz = pbuffer.data(idx_dip_fg + 82);

    auto tr_x_xzz_xyzz = pbuffer.data(idx_dip_fg + 83);

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

    auto tr_x_yyz_yyyz = pbuffer.data(idx_dip_fg + 116);

    auto tr_x_yyz_yyzz = pbuffer.data(idx_dip_fg + 117);

    auto tr_x_yyz_yzzz = pbuffer.data(idx_dip_fg + 118);

    auto tr_x_yyz_zzzz = pbuffer.data(idx_dip_fg + 119);

    auto tr_x_yzz_xxxx = pbuffer.data(idx_dip_fg + 120);

    auto tr_x_yzz_xxxz = pbuffer.data(idx_dip_fg + 122);

    auto tr_x_yzz_xxyz = pbuffer.data(idx_dip_fg + 124);

    auto tr_x_yzz_xxzz = pbuffer.data(idx_dip_fg + 125);

    auto tr_x_yzz_xyyz = pbuffer.data(idx_dip_fg + 127);

    auto tr_x_yzz_xyzz = pbuffer.data(idx_dip_fg + 128);

    auto tr_x_yzz_xzzz = pbuffer.data(idx_dip_fg + 129);

    auto tr_x_yzz_yyyy = pbuffer.data(idx_dip_fg + 130);

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

    auto tr_y_xxy_xxxx = pbuffer.data(idx_dip_fg + 165);

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

    auto tr_y_xxz_xxxz = pbuffer.data(idx_dip_fg + 182);

    auto tr_y_xxz_xxyy = pbuffer.data(idx_dip_fg + 183);

    auto tr_y_xxz_xxzz = pbuffer.data(idx_dip_fg + 185);

    auto tr_y_xxz_xyyy = pbuffer.data(idx_dip_fg + 186);

    auto tr_y_xxz_xzzz = pbuffer.data(idx_dip_fg + 189);

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

    auto tr_y_yyz_xxxz = pbuffer.data(idx_dip_fg + 257);

    auto tr_y_yyz_xxyy = pbuffer.data(idx_dip_fg + 258);

    auto tr_y_yyz_xxyz = pbuffer.data(idx_dip_fg + 259);

    auto tr_y_yyz_xxzz = pbuffer.data(idx_dip_fg + 260);

    auto tr_y_yyz_xyyy = pbuffer.data(idx_dip_fg + 261);

    auto tr_y_yyz_xyyz = pbuffer.data(idx_dip_fg + 262);

    auto tr_y_yyz_xyzz = pbuffer.data(idx_dip_fg + 263);

    auto tr_y_yyz_xzzz = pbuffer.data(idx_dip_fg + 264);

    auto tr_y_yyz_yyyy = pbuffer.data(idx_dip_fg + 265);

    auto tr_y_yyz_yyyz = pbuffer.data(idx_dip_fg + 266);

    auto tr_y_yyz_yyzz = pbuffer.data(idx_dip_fg + 267);

    auto tr_y_yyz_yzzz = pbuffer.data(idx_dip_fg + 268);

    auto tr_y_yyz_zzzz = pbuffer.data(idx_dip_fg + 269);

    auto tr_y_yzz_xxxx = pbuffer.data(idx_dip_fg + 270);

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

    auto tr_z_xxz_xxxy = pbuffer.data(idx_dip_fg + 331);

    auto tr_z_xxz_xxxz = pbuffer.data(idx_dip_fg + 332);

    auto tr_z_xxz_xxyy = pbuffer.data(idx_dip_fg + 333);

    auto tr_z_xxz_xxyz = pbuffer.data(idx_dip_fg + 334);

    auto tr_z_xxz_xxzz = pbuffer.data(idx_dip_fg + 335);

    auto tr_z_xxz_xyyy = pbuffer.data(idx_dip_fg + 336);

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

    auto tr_z_yyz_xxxy = pbuffer.data(idx_dip_fg + 406);

    auto tr_z_yyz_xxxz = pbuffer.data(idx_dip_fg + 407);

    auto tr_z_yyz_xxyy = pbuffer.data(idx_dip_fg + 408);

    auto tr_z_yyz_xxyz = pbuffer.data(idx_dip_fg + 409);

    auto tr_z_yyz_xxzz = pbuffer.data(idx_dip_fg + 410);

    auto tr_z_yyz_xyyy = pbuffer.data(idx_dip_fg + 411);

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

    // Set up 0-15 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,               \
                             tr_x_xx_xxxx,   \
                             tr_x_xx_xxxy,   \
                             tr_x_xx_xxxz,   \
                             tr_x_xx_xxyy,   \
                             tr_x_xx_xxyz,   \
                             tr_x_xx_xxzz,   \
                             tr_x_xx_xyyy,   \
                             tr_x_xx_xyyz,   \
                             tr_x_xx_xyzz,   \
                             tr_x_xx_xzzz,   \
                             tr_x_xx_yyyy,   \
                             tr_x_xx_yyyz,   \
                             tr_x_xx_yyzz,   \
                             tr_x_xx_yzzz,   \
                             tr_x_xx_zzzz,   \
                             tr_x_xxx_xxx,   \
                             tr_x_xxx_xxxx,  \
                             tr_x_xxx_xxxy,  \
                             tr_x_xxx_xxxz,  \
                             tr_x_xxx_xxy,   \
                             tr_x_xxx_xxyy,  \
                             tr_x_xxx_xxyz,  \
                             tr_x_xxx_xxz,   \
                             tr_x_xxx_xxzz,  \
                             tr_x_xxx_xyy,   \
                             tr_x_xxx_xyyy,  \
                             tr_x_xxx_xyyz,  \
                             tr_x_xxx_xyz,   \
                             tr_x_xxx_xyzz,  \
                             tr_x_xxx_xzz,   \
                             tr_x_xxx_xzzz,  \
                             tr_x_xxx_yyy,   \
                             tr_x_xxx_yyyy,  \
                             tr_x_xxx_yyyz,  \
                             tr_x_xxx_yyz,   \
                             tr_x_xxx_yyzz,  \
                             tr_x_xxx_yzz,   \
                             tr_x_xxx_yzzz,  \
                             tr_x_xxx_zzz,   \
                             tr_x_xxx_zzzz,  \
                             tr_x_xxxx_xxxx, \
                             tr_x_xxxx_xxxy, \
                             tr_x_xxxx_xxxz, \
                             tr_x_xxxx_xxyy, \
                             tr_x_xxxx_xxyz, \
                             tr_x_xxxx_xxzz, \
                             tr_x_xxxx_xyyy, \
                             tr_x_xxxx_xyyz, \
                             tr_x_xxxx_xyzz, \
                             tr_x_xxxx_xzzz, \
                             tr_x_xxxx_yyyy, \
                             tr_x_xxxx_yyyz, \
                             tr_x_xxxx_yyzz, \
                             tr_x_xxxx_yzzz, \
                             tr_x_xxxx_zzzz, \
                             ts_xxx_xxxx,    \
                             ts_xxx_xxxy,    \
                             ts_xxx_xxxz,    \
                             ts_xxx_xxyy,    \
                             ts_xxx_xxyz,    \
                             ts_xxx_xxzz,    \
                             ts_xxx_xyyy,    \
                             ts_xxx_xyyz,    \
                             ts_xxx_xyzz,    \
                             ts_xxx_xzzz,    \
                             ts_xxx_yyyy,    \
                             ts_xxx_yyyz,    \
                             ts_xxx_yyzz,    \
                             ts_xxx_yzzz,    \
                             ts_xxx_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxx_xxxx[i] = 3.0 * tr_x_xx_xxxx[i] * fe_0 + 4.0 * tr_x_xxx_xxx[i] * fe_0 + ts_xxx_xxxx[i] * fe_0 + tr_x_xxx_xxxx[i] * pa_x[i];

        tr_x_xxxx_xxxy[i] = 3.0 * tr_x_xx_xxxy[i] * fe_0 + 3.0 * tr_x_xxx_xxy[i] * fe_0 + ts_xxx_xxxy[i] * fe_0 + tr_x_xxx_xxxy[i] * pa_x[i];

        tr_x_xxxx_xxxz[i] = 3.0 * tr_x_xx_xxxz[i] * fe_0 + 3.0 * tr_x_xxx_xxz[i] * fe_0 + ts_xxx_xxxz[i] * fe_0 + tr_x_xxx_xxxz[i] * pa_x[i];

        tr_x_xxxx_xxyy[i] = 3.0 * tr_x_xx_xxyy[i] * fe_0 + 2.0 * tr_x_xxx_xyy[i] * fe_0 + ts_xxx_xxyy[i] * fe_0 + tr_x_xxx_xxyy[i] * pa_x[i];

        tr_x_xxxx_xxyz[i] = 3.0 * tr_x_xx_xxyz[i] * fe_0 + 2.0 * tr_x_xxx_xyz[i] * fe_0 + ts_xxx_xxyz[i] * fe_0 + tr_x_xxx_xxyz[i] * pa_x[i];

        tr_x_xxxx_xxzz[i] = 3.0 * tr_x_xx_xxzz[i] * fe_0 + 2.0 * tr_x_xxx_xzz[i] * fe_0 + ts_xxx_xxzz[i] * fe_0 + tr_x_xxx_xxzz[i] * pa_x[i];

        tr_x_xxxx_xyyy[i] = 3.0 * tr_x_xx_xyyy[i] * fe_0 + tr_x_xxx_yyy[i] * fe_0 + ts_xxx_xyyy[i] * fe_0 + tr_x_xxx_xyyy[i] * pa_x[i];

        tr_x_xxxx_xyyz[i] = 3.0 * tr_x_xx_xyyz[i] * fe_0 + tr_x_xxx_yyz[i] * fe_0 + ts_xxx_xyyz[i] * fe_0 + tr_x_xxx_xyyz[i] * pa_x[i];

        tr_x_xxxx_xyzz[i] = 3.0 * tr_x_xx_xyzz[i] * fe_0 + tr_x_xxx_yzz[i] * fe_0 + ts_xxx_xyzz[i] * fe_0 + tr_x_xxx_xyzz[i] * pa_x[i];

        tr_x_xxxx_xzzz[i] = 3.0 * tr_x_xx_xzzz[i] * fe_0 + tr_x_xxx_zzz[i] * fe_0 + ts_xxx_xzzz[i] * fe_0 + tr_x_xxx_xzzz[i] * pa_x[i];

        tr_x_xxxx_yyyy[i] = 3.0 * tr_x_xx_yyyy[i] * fe_0 + ts_xxx_yyyy[i] * fe_0 + tr_x_xxx_yyyy[i] * pa_x[i];

        tr_x_xxxx_yyyz[i] = 3.0 * tr_x_xx_yyyz[i] * fe_0 + ts_xxx_yyyz[i] * fe_0 + tr_x_xxx_yyyz[i] * pa_x[i];

        tr_x_xxxx_yyzz[i] = 3.0 * tr_x_xx_yyzz[i] * fe_0 + ts_xxx_yyzz[i] * fe_0 + tr_x_xxx_yyzz[i] * pa_x[i];

        tr_x_xxxx_yzzz[i] = 3.0 * tr_x_xx_yzzz[i] * fe_0 + ts_xxx_yzzz[i] * fe_0 + tr_x_xxx_yzzz[i] * pa_x[i];

        tr_x_xxxx_zzzz[i] = 3.0 * tr_x_xx_zzzz[i] * fe_0 + ts_xxx_zzzz[i] * fe_0 + tr_x_xxx_zzzz[i] * pa_x[i];
    }

    // Set up 15-30 components of targeted buffer : GG

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

    auto tr_x_xxxy_yyyz = pbuffer.data(idx_dip_gg + 26);

    auto tr_x_xxxy_yyzz = pbuffer.data(idx_dip_gg + 27);

    auto tr_x_xxxy_yzzz = pbuffer.data(idx_dip_gg + 28);

    auto tr_x_xxxy_zzzz = pbuffer.data(idx_dip_gg + 29);

#pragma omp simd aligned(pa_y,               \
                             tr_x_xxx_xxx,   \
                             tr_x_xxx_xxxx,  \
                             tr_x_xxx_xxxy,  \
                             tr_x_xxx_xxxz,  \
                             tr_x_xxx_xxy,   \
                             tr_x_xxx_xxyy,  \
                             tr_x_xxx_xxyz,  \
                             tr_x_xxx_xxz,   \
                             tr_x_xxx_xxzz,  \
                             tr_x_xxx_xyy,   \
                             tr_x_xxx_xyyy,  \
                             tr_x_xxx_xyyz,  \
                             tr_x_xxx_xyz,   \
                             tr_x_xxx_xyzz,  \
                             tr_x_xxx_xzz,   \
                             tr_x_xxx_xzzz,  \
                             tr_x_xxx_yyy,   \
                             tr_x_xxx_yyyy,  \
                             tr_x_xxx_yyyz,  \
                             tr_x_xxx_yyz,   \
                             tr_x_xxx_yyzz,  \
                             tr_x_xxx_yzz,   \
                             tr_x_xxx_yzzz,  \
                             tr_x_xxx_zzz,   \
                             tr_x_xxx_zzzz,  \
                             tr_x_xxxy_xxxx, \
                             tr_x_xxxy_xxxy, \
                             tr_x_xxxy_xxxz, \
                             tr_x_xxxy_xxyy, \
                             tr_x_xxxy_xxyz, \
                             tr_x_xxxy_xxzz, \
                             tr_x_xxxy_xyyy, \
                             tr_x_xxxy_xyyz, \
                             tr_x_xxxy_xyzz, \
                             tr_x_xxxy_xzzz, \
                             tr_x_xxxy_yyyy, \
                             tr_x_xxxy_yyyz, \
                             tr_x_xxxy_yyzz, \
                             tr_x_xxxy_yzzz, \
                             tr_x_xxxy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxy_xxxx[i] = tr_x_xxx_xxxx[i] * pa_y[i];

        tr_x_xxxy_xxxy[i] = tr_x_xxx_xxx[i] * fe_0 + tr_x_xxx_xxxy[i] * pa_y[i];

        tr_x_xxxy_xxxz[i] = tr_x_xxx_xxxz[i] * pa_y[i];

        tr_x_xxxy_xxyy[i] = 2.0 * tr_x_xxx_xxy[i] * fe_0 + tr_x_xxx_xxyy[i] * pa_y[i];

        tr_x_xxxy_xxyz[i] = tr_x_xxx_xxz[i] * fe_0 + tr_x_xxx_xxyz[i] * pa_y[i];

        tr_x_xxxy_xxzz[i] = tr_x_xxx_xxzz[i] * pa_y[i];

        tr_x_xxxy_xyyy[i] = 3.0 * tr_x_xxx_xyy[i] * fe_0 + tr_x_xxx_xyyy[i] * pa_y[i];

        tr_x_xxxy_xyyz[i] = 2.0 * tr_x_xxx_xyz[i] * fe_0 + tr_x_xxx_xyyz[i] * pa_y[i];

        tr_x_xxxy_xyzz[i] = tr_x_xxx_xzz[i] * fe_0 + tr_x_xxx_xyzz[i] * pa_y[i];

        tr_x_xxxy_xzzz[i] = tr_x_xxx_xzzz[i] * pa_y[i];

        tr_x_xxxy_yyyy[i] = 4.0 * tr_x_xxx_yyy[i] * fe_0 + tr_x_xxx_yyyy[i] * pa_y[i];

        tr_x_xxxy_yyyz[i] = 3.0 * tr_x_xxx_yyz[i] * fe_0 + tr_x_xxx_yyyz[i] * pa_y[i];

        tr_x_xxxy_yyzz[i] = 2.0 * tr_x_xxx_yzz[i] * fe_0 + tr_x_xxx_yyzz[i] * pa_y[i];

        tr_x_xxxy_yzzz[i] = tr_x_xxx_zzz[i] * fe_0 + tr_x_xxx_yzzz[i] * pa_y[i];

        tr_x_xxxy_zzzz[i] = tr_x_xxx_zzzz[i] * pa_y[i];
    }

    // Set up 30-45 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_z,               \
                             tr_x_xxx_xxx,   \
                             tr_x_xxx_xxxx,  \
                             tr_x_xxx_xxxy,  \
                             tr_x_xxx_xxxz,  \
                             tr_x_xxx_xxy,   \
                             tr_x_xxx_xxyy,  \
                             tr_x_xxx_xxyz,  \
                             tr_x_xxx_xxz,   \
                             tr_x_xxx_xxzz,  \
                             tr_x_xxx_xyy,   \
                             tr_x_xxx_xyyy,  \
                             tr_x_xxx_xyyz,  \
                             tr_x_xxx_xyz,   \
                             tr_x_xxx_xyzz,  \
                             tr_x_xxx_xzz,   \
                             tr_x_xxx_xzzz,  \
                             tr_x_xxx_yyy,   \
                             tr_x_xxx_yyyy,  \
                             tr_x_xxx_yyyz,  \
                             tr_x_xxx_yyz,   \
                             tr_x_xxx_yyzz,  \
                             tr_x_xxx_yzz,   \
                             tr_x_xxx_yzzz,  \
                             tr_x_xxx_zzz,   \
                             tr_x_xxx_zzzz,  \
                             tr_x_xxxz_xxxx, \
                             tr_x_xxxz_xxxy, \
                             tr_x_xxxz_xxxz, \
                             tr_x_xxxz_xxyy, \
                             tr_x_xxxz_xxyz, \
                             tr_x_xxxz_xxzz, \
                             tr_x_xxxz_xyyy, \
                             tr_x_xxxz_xyyz, \
                             tr_x_xxxz_xyzz, \
                             tr_x_xxxz_xzzz, \
                             tr_x_xxxz_yyyy, \
                             tr_x_xxxz_yyyz, \
                             tr_x_xxxz_yyzz, \
                             tr_x_xxxz_yzzz, \
                             tr_x_xxxz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxz_xxxx[i] = tr_x_xxx_xxxx[i] * pa_z[i];

        tr_x_xxxz_xxxy[i] = tr_x_xxx_xxxy[i] * pa_z[i];

        tr_x_xxxz_xxxz[i] = tr_x_xxx_xxx[i] * fe_0 + tr_x_xxx_xxxz[i] * pa_z[i];

        tr_x_xxxz_xxyy[i] = tr_x_xxx_xxyy[i] * pa_z[i];

        tr_x_xxxz_xxyz[i] = tr_x_xxx_xxy[i] * fe_0 + tr_x_xxx_xxyz[i] * pa_z[i];

        tr_x_xxxz_xxzz[i] = 2.0 * tr_x_xxx_xxz[i] * fe_0 + tr_x_xxx_xxzz[i] * pa_z[i];

        tr_x_xxxz_xyyy[i] = tr_x_xxx_xyyy[i] * pa_z[i];

        tr_x_xxxz_xyyz[i] = tr_x_xxx_xyy[i] * fe_0 + tr_x_xxx_xyyz[i] * pa_z[i];

        tr_x_xxxz_xyzz[i] = 2.0 * tr_x_xxx_xyz[i] * fe_0 + tr_x_xxx_xyzz[i] * pa_z[i];

        tr_x_xxxz_xzzz[i] = 3.0 * tr_x_xxx_xzz[i] * fe_0 + tr_x_xxx_xzzz[i] * pa_z[i];

        tr_x_xxxz_yyyy[i] = tr_x_xxx_yyyy[i] * pa_z[i];

        tr_x_xxxz_yyyz[i] = tr_x_xxx_yyy[i] * fe_0 + tr_x_xxx_yyyz[i] * pa_z[i];

        tr_x_xxxz_yyzz[i] = 2.0 * tr_x_xxx_yyz[i] * fe_0 + tr_x_xxx_yyzz[i] * pa_z[i];

        tr_x_xxxz_yzzz[i] = 3.0 * tr_x_xxx_yzz[i] * fe_0 + tr_x_xxx_yzzz[i] * pa_z[i];

        tr_x_xxxz_zzzz[i] = 4.0 * tr_x_xxx_zzz[i] * fe_0 + tr_x_xxx_zzzz[i] * pa_z[i];
    }

    // Set up 45-60 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_x_xx_xxxx,   \
                             tr_x_xx_xxxy,   \
                             tr_x_xx_xxxz,   \
                             tr_x_xx_xxyy,   \
                             tr_x_xx_xxyz,   \
                             tr_x_xx_xxzz,   \
                             tr_x_xx_xyyy,   \
                             tr_x_xx_xyyz,   \
                             tr_x_xx_xyzz,   \
                             tr_x_xx_xzzz,   \
                             tr_x_xx_zzzz,   \
                             tr_x_xxy_xxx,   \
                             tr_x_xxy_xxxx,  \
                             tr_x_xxy_xxxy,  \
                             tr_x_xxy_xxxz,  \
                             tr_x_xxy_xxy,   \
                             tr_x_xxy_xxyy,  \
                             tr_x_xxy_xxyz,  \
                             tr_x_xxy_xxz,   \
                             tr_x_xxy_xxzz,  \
                             tr_x_xxy_xyy,   \
                             tr_x_xxy_xyyy,  \
                             tr_x_xxy_xyyz,  \
                             tr_x_xxy_xyz,   \
                             tr_x_xxy_xyzz,  \
                             tr_x_xxy_xzz,   \
                             tr_x_xxy_xzzz,  \
                             tr_x_xxy_zzzz,  \
                             tr_x_xxyy_xxxx, \
                             tr_x_xxyy_xxxy, \
                             tr_x_xxyy_xxxz, \
                             tr_x_xxyy_xxyy, \
                             tr_x_xxyy_xxyz, \
                             tr_x_xxyy_xxzz, \
                             tr_x_xxyy_xyyy, \
                             tr_x_xxyy_xyyz, \
                             tr_x_xxyy_xyzz, \
                             tr_x_xxyy_xzzz, \
                             tr_x_xxyy_yyyy, \
                             tr_x_xxyy_yyyz, \
                             tr_x_xxyy_yyzz, \
                             tr_x_xxyy_yzzz, \
                             tr_x_xxyy_zzzz, \
                             tr_x_xyy_yyyy,  \
                             tr_x_xyy_yyyz,  \
                             tr_x_xyy_yyzz,  \
                             tr_x_xyy_yzzz,  \
                             tr_x_yy_yyyy,   \
                             tr_x_yy_yyyz,   \
                             tr_x_yy_yyzz,   \
                             tr_x_yy_yzzz,   \
                             ts_xyy_yyyy,    \
                             ts_xyy_yyyz,    \
                             ts_xyy_yyzz,    \
                             ts_xyy_yzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyy_xxxx[i] = tr_x_xx_xxxx[i] * fe_0 + tr_x_xxy_xxxx[i] * pa_y[i];

        tr_x_xxyy_xxxy[i] = tr_x_xx_xxxy[i] * fe_0 + tr_x_xxy_xxx[i] * fe_0 + tr_x_xxy_xxxy[i] * pa_y[i];

        tr_x_xxyy_xxxz[i] = tr_x_xx_xxxz[i] * fe_0 + tr_x_xxy_xxxz[i] * pa_y[i];

        tr_x_xxyy_xxyy[i] = tr_x_xx_xxyy[i] * fe_0 + 2.0 * tr_x_xxy_xxy[i] * fe_0 + tr_x_xxy_xxyy[i] * pa_y[i];

        tr_x_xxyy_xxyz[i] = tr_x_xx_xxyz[i] * fe_0 + tr_x_xxy_xxz[i] * fe_0 + tr_x_xxy_xxyz[i] * pa_y[i];

        tr_x_xxyy_xxzz[i] = tr_x_xx_xxzz[i] * fe_0 + tr_x_xxy_xxzz[i] * pa_y[i];

        tr_x_xxyy_xyyy[i] = tr_x_xx_xyyy[i] * fe_0 + 3.0 * tr_x_xxy_xyy[i] * fe_0 + tr_x_xxy_xyyy[i] * pa_y[i];

        tr_x_xxyy_xyyz[i] = tr_x_xx_xyyz[i] * fe_0 + 2.0 * tr_x_xxy_xyz[i] * fe_0 + tr_x_xxy_xyyz[i] * pa_y[i];

        tr_x_xxyy_xyzz[i] = tr_x_xx_xyzz[i] * fe_0 + tr_x_xxy_xzz[i] * fe_0 + tr_x_xxy_xyzz[i] * pa_y[i];

        tr_x_xxyy_xzzz[i] = tr_x_xx_xzzz[i] * fe_0 + tr_x_xxy_xzzz[i] * pa_y[i];

        tr_x_xxyy_yyyy[i] = tr_x_yy_yyyy[i] * fe_0 + ts_xyy_yyyy[i] * fe_0 + tr_x_xyy_yyyy[i] * pa_x[i];

        tr_x_xxyy_yyyz[i] = tr_x_yy_yyyz[i] * fe_0 + ts_xyy_yyyz[i] * fe_0 + tr_x_xyy_yyyz[i] * pa_x[i];

        tr_x_xxyy_yyzz[i] = tr_x_yy_yyzz[i] * fe_0 + ts_xyy_yyzz[i] * fe_0 + tr_x_xyy_yyzz[i] * pa_x[i];

        tr_x_xxyy_yzzz[i] = tr_x_yy_yzzz[i] * fe_0 + ts_xyy_yzzz[i] * fe_0 + tr_x_xyy_yzzz[i] * pa_x[i];

        tr_x_xxyy_zzzz[i] = tr_x_xx_zzzz[i] * fe_0 + tr_x_xxy_zzzz[i] * pa_y[i];
    }

    // Set up 60-75 components of targeted buffer : GG

    auto tr_x_xxyz_xxxx = pbuffer.data(idx_dip_gg + 60);

    auto tr_x_xxyz_xxxy = pbuffer.data(idx_dip_gg + 61);

    auto tr_x_xxyz_xxxz = pbuffer.data(idx_dip_gg + 62);

    auto tr_x_xxyz_xxyy = pbuffer.data(idx_dip_gg + 63);

    auto tr_x_xxyz_xxyz = pbuffer.data(idx_dip_gg + 64);

    auto tr_x_xxyz_xxzz = pbuffer.data(idx_dip_gg + 65);

    auto tr_x_xxyz_xyyy = pbuffer.data(idx_dip_gg + 66);

    auto tr_x_xxyz_xyyz = pbuffer.data(idx_dip_gg + 67);

    auto tr_x_xxyz_xyzz = pbuffer.data(idx_dip_gg + 68);

    auto tr_x_xxyz_xzzz = pbuffer.data(idx_dip_gg + 69);

    auto tr_x_xxyz_yyyy = pbuffer.data(idx_dip_gg + 70);

    auto tr_x_xxyz_yyyz = pbuffer.data(idx_dip_gg + 71);

    auto tr_x_xxyz_yyzz = pbuffer.data(idx_dip_gg + 72);

    auto tr_x_xxyz_yzzz = pbuffer.data(idx_dip_gg + 73);

    auto tr_x_xxyz_zzzz = pbuffer.data(idx_dip_gg + 74);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_xxy_xxxy,  \
                             tr_x_xxy_xxyy,  \
                             tr_x_xxy_xyyy,  \
                             tr_x_xxy_yyyy,  \
                             tr_x_xxyz_xxxx, \
                             tr_x_xxyz_xxxy, \
                             tr_x_xxyz_xxxz, \
                             tr_x_xxyz_xxyy, \
                             tr_x_xxyz_xxyz, \
                             tr_x_xxyz_xxzz, \
                             tr_x_xxyz_xyyy, \
                             tr_x_xxyz_xyyz, \
                             tr_x_xxyz_xyzz, \
                             tr_x_xxyz_xzzz, \
                             tr_x_xxyz_yyyy, \
                             tr_x_xxyz_yyyz, \
                             tr_x_xxyz_yyzz, \
                             tr_x_xxyz_yzzz, \
                             tr_x_xxyz_zzzz, \
                             tr_x_xxz_xxxx,  \
                             tr_x_xxz_xxxz,  \
                             tr_x_xxz_xxyz,  \
                             tr_x_xxz_xxz,   \
                             tr_x_xxz_xxzz,  \
                             tr_x_xxz_xyyz,  \
                             tr_x_xxz_xyz,   \
                             tr_x_xxz_xyzz,  \
                             tr_x_xxz_xzz,   \
                             tr_x_xxz_xzzz,  \
                             tr_x_xxz_yyyz,  \
                             tr_x_xxz_yyz,   \
                             tr_x_xxz_yyzz,  \
                             tr_x_xxz_yzz,   \
                             tr_x_xxz_yzzz,  \
                             tr_x_xxz_zzz,   \
                             tr_x_xxz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyz_xxxx[i] = tr_x_xxz_xxxx[i] * pa_y[i];

        tr_x_xxyz_xxxy[i] = tr_x_xxy_xxxy[i] * pa_z[i];

        tr_x_xxyz_xxxz[i] = tr_x_xxz_xxxz[i] * pa_y[i];

        tr_x_xxyz_xxyy[i] = tr_x_xxy_xxyy[i] * pa_z[i];

        tr_x_xxyz_xxyz[i] = tr_x_xxz_xxz[i] * fe_0 + tr_x_xxz_xxyz[i] * pa_y[i];

        tr_x_xxyz_xxzz[i] = tr_x_xxz_xxzz[i] * pa_y[i];

        tr_x_xxyz_xyyy[i] = tr_x_xxy_xyyy[i] * pa_z[i];

        tr_x_xxyz_xyyz[i] = 2.0 * tr_x_xxz_xyz[i] * fe_0 + tr_x_xxz_xyyz[i] * pa_y[i];

        tr_x_xxyz_xyzz[i] = tr_x_xxz_xzz[i] * fe_0 + tr_x_xxz_xyzz[i] * pa_y[i];

        tr_x_xxyz_xzzz[i] = tr_x_xxz_xzzz[i] * pa_y[i];

        tr_x_xxyz_yyyy[i] = tr_x_xxy_yyyy[i] * pa_z[i];

        tr_x_xxyz_yyyz[i] = 3.0 * tr_x_xxz_yyz[i] * fe_0 + tr_x_xxz_yyyz[i] * pa_y[i];

        tr_x_xxyz_yyzz[i] = 2.0 * tr_x_xxz_yzz[i] * fe_0 + tr_x_xxz_yyzz[i] * pa_y[i];

        tr_x_xxyz_yzzz[i] = tr_x_xxz_zzz[i] * fe_0 + tr_x_xxz_yzzz[i] * pa_y[i];

        tr_x_xxyz_zzzz[i] = tr_x_xxz_zzzz[i] * pa_y[i];
    }

    // Set up 75-90 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_x_xx_xxxx,   \
                             tr_x_xx_xxxy,   \
                             tr_x_xx_xxxz,   \
                             tr_x_xx_xxyy,   \
                             tr_x_xx_xxyz,   \
                             tr_x_xx_xxzz,   \
                             tr_x_xx_xyyy,   \
                             tr_x_xx_xyyz,   \
                             tr_x_xx_xyzz,   \
                             tr_x_xx_xzzz,   \
                             tr_x_xx_yyyy,   \
                             tr_x_xxz_xxx,   \
                             tr_x_xxz_xxxx,  \
                             tr_x_xxz_xxxy,  \
                             tr_x_xxz_xxxz,  \
                             tr_x_xxz_xxy,   \
                             tr_x_xxz_xxyy,  \
                             tr_x_xxz_xxyz,  \
                             tr_x_xxz_xxz,   \
                             tr_x_xxz_xxzz,  \
                             tr_x_xxz_xyy,   \
                             tr_x_xxz_xyyy,  \
                             tr_x_xxz_xyyz,  \
                             tr_x_xxz_xyz,   \
                             tr_x_xxz_xyzz,  \
                             tr_x_xxz_xzz,   \
                             tr_x_xxz_xzzz,  \
                             tr_x_xxz_yyyy,  \
                             tr_x_xxzz_xxxx, \
                             tr_x_xxzz_xxxy, \
                             tr_x_xxzz_xxxz, \
                             tr_x_xxzz_xxyy, \
                             tr_x_xxzz_xxyz, \
                             tr_x_xxzz_xxzz, \
                             tr_x_xxzz_xyyy, \
                             tr_x_xxzz_xyyz, \
                             tr_x_xxzz_xyzz, \
                             tr_x_xxzz_xzzz, \
                             tr_x_xxzz_yyyy, \
                             tr_x_xxzz_yyyz, \
                             tr_x_xxzz_yyzz, \
                             tr_x_xxzz_yzzz, \
                             tr_x_xxzz_zzzz, \
                             tr_x_xzz_yyyz,  \
                             tr_x_xzz_yyzz,  \
                             tr_x_xzz_yzzz,  \
                             tr_x_xzz_zzzz,  \
                             tr_x_zz_yyyz,   \
                             tr_x_zz_yyzz,   \
                             tr_x_zz_yzzz,   \
                             tr_x_zz_zzzz,   \
                             ts_xzz_yyyz,    \
                             ts_xzz_yyzz,    \
                             ts_xzz_yzzz,    \
                             ts_xzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzz_xxxx[i] = tr_x_xx_xxxx[i] * fe_0 + tr_x_xxz_xxxx[i] * pa_z[i];

        tr_x_xxzz_xxxy[i] = tr_x_xx_xxxy[i] * fe_0 + tr_x_xxz_xxxy[i] * pa_z[i];

        tr_x_xxzz_xxxz[i] = tr_x_xx_xxxz[i] * fe_0 + tr_x_xxz_xxx[i] * fe_0 + tr_x_xxz_xxxz[i] * pa_z[i];

        tr_x_xxzz_xxyy[i] = tr_x_xx_xxyy[i] * fe_0 + tr_x_xxz_xxyy[i] * pa_z[i];

        tr_x_xxzz_xxyz[i] = tr_x_xx_xxyz[i] * fe_0 + tr_x_xxz_xxy[i] * fe_0 + tr_x_xxz_xxyz[i] * pa_z[i];

        tr_x_xxzz_xxzz[i] = tr_x_xx_xxzz[i] * fe_0 + 2.0 * tr_x_xxz_xxz[i] * fe_0 + tr_x_xxz_xxzz[i] * pa_z[i];

        tr_x_xxzz_xyyy[i] = tr_x_xx_xyyy[i] * fe_0 + tr_x_xxz_xyyy[i] * pa_z[i];

        tr_x_xxzz_xyyz[i] = tr_x_xx_xyyz[i] * fe_0 + tr_x_xxz_xyy[i] * fe_0 + tr_x_xxz_xyyz[i] * pa_z[i];

        tr_x_xxzz_xyzz[i] = tr_x_xx_xyzz[i] * fe_0 + 2.0 * tr_x_xxz_xyz[i] * fe_0 + tr_x_xxz_xyzz[i] * pa_z[i];

        tr_x_xxzz_xzzz[i] = tr_x_xx_xzzz[i] * fe_0 + 3.0 * tr_x_xxz_xzz[i] * fe_0 + tr_x_xxz_xzzz[i] * pa_z[i];

        tr_x_xxzz_yyyy[i] = tr_x_xx_yyyy[i] * fe_0 + tr_x_xxz_yyyy[i] * pa_z[i];

        tr_x_xxzz_yyyz[i] = tr_x_zz_yyyz[i] * fe_0 + ts_xzz_yyyz[i] * fe_0 + tr_x_xzz_yyyz[i] * pa_x[i];

        tr_x_xxzz_yyzz[i] = tr_x_zz_yyzz[i] * fe_0 + ts_xzz_yyzz[i] * fe_0 + tr_x_xzz_yyzz[i] * pa_x[i];

        tr_x_xxzz_yzzz[i] = tr_x_zz_yzzz[i] * fe_0 + ts_xzz_yzzz[i] * fe_0 + tr_x_xzz_yzzz[i] * pa_x[i];

        tr_x_xxzz_zzzz[i] = tr_x_zz_zzzz[i] * fe_0 + ts_xzz_zzzz[i] * fe_0 + tr_x_xzz_zzzz[i] * pa_x[i];
    }

    // Set up 90-105 components of targeted buffer : GG

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

    auto tr_x_xyyy_zzzz = pbuffer.data(idx_dip_gg + 104);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_x_xy_xxxx,   \
                             tr_x_xy_xxxz,   \
                             tr_x_xy_xxzz,   \
                             tr_x_xy_xzzz,   \
                             tr_x_xyy_xxxx,  \
                             tr_x_xyy_xxxz,  \
                             tr_x_xyy_xxzz,  \
                             tr_x_xyy_xzzz,  \
                             tr_x_xyyy_xxxx, \
                             tr_x_xyyy_xxxy, \
                             tr_x_xyyy_xxxz, \
                             tr_x_xyyy_xxyy, \
                             tr_x_xyyy_xxyz, \
                             tr_x_xyyy_xxzz, \
                             tr_x_xyyy_xyyy, \
                             tr_x_xyyy_xyyz, \
                             tr_x_xyyy_xyzz, \
                             tr_x_xyyy_xzzz, \
                             tr_x_xyyy_yyyy, \
                             tr_x_xyyy_yyyz, \
                             tr_x_xyyy_yyzz, \
                             tr_x_xyyy_yzzz, \
                             tr_x_xyyy_zzzz, \
                             tr_x_yyy_xxxy,  \
                             tr_x_yyy_xxy,   \
                             tr_x_yyy_xxyy,  \
                             tr_x_yyy_xxyz,  \
                             tr_x_yyy_xyy,   \
                             tr_x_yyy_xyyy,  \
                             tr_x_yyy_xyyz,  \
                             tr_x_yyy_xyz,   \
                             tr_x_yyy_xyzz,  \
                             tr_x_yyy_yyy,   \
                             tr_x_yyy_yyyy,  \
                             tr_x_yyy_yyyz,  \
                             tr_x_yyy_yyz,   \
                             tr_x_yyy_yyzz,  \
                             tr_x_yyy_yzz,   \
                             tr_x_yyy_yzzz,  \
                             tr_x_yyy_zzzz,  \
                             ts_yyy_xxxy,    \
                             ts_yyy_xxyy,    \
                             ts_yyy_xxyz,    \
                             ts_yyy_xyyy,    \
                             ts_yyy_xyyz,    \
                             ts_yyy_xyzz,    \
                             ts_yyy_yyyy,    \
                             ts_yyy_yyyz,    \
                             ts_yyy_yyzz,    \
                             ts_yyy_yzzz,    \
                             ts_yyy_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyy_xxxx[i] = 2.0 * tr_x_xy_xxxx[i] * fe_0 + tr_x_xyy_xxxx[i] * pa_y[i];

        tr_x_xyyy_xxxy[i] = 3.0 * tr_x_yyy_xxy[i] * fe_0 + ts_yyy_xxxy[i] * fe_0 + tr_x_yyy_xxxy[i] * pa_x[i];

        tr_x_xyyy_xxxz[i] = 2.0 * tr_x_xy_xxxz[i] * fe_0 + tr_x_xyy_xxxz[i] * pa_y[i];

        tr_x_xyyy_xxyy[i] = 2.0 * tr_x_yyy_xyy[i] * fe_0 + ts_yyy_xxyy[i] * fe_0 + tr_x_yyy_xxyy[i] * pa_x[i];

        tr_x_xyyy_xxyz[i] = 2.0 * tr_x_yyy_xyz[i] * fe_0 + ts_yyy_xxyz[i] * fe_0 + tr_x_yyy_xxyz[i] * pa_x[i];

        tr_x_xyyy_xxzz[i] = 2.0 * tr_x_xy_xxzz[i] * fe_0 + tr_x_xyy_xxzz[i] * pa_y[i];

        tr_x_xyyy_xyyy[i] = tr_x_yyy_yyy[i] * fe_0 + ts_yyy_xyyy[i] * fe_0 + tr_x_yyy_xyyy[i] * pa_x[i];

        tr_x_xyyy_xyyz[i] = tr_x_yyy_yyz[i] * fe_0 + ts_yyy_xyyz[i] * fe_0 + tr_x_yyy_xyyz[i] * pa_x[i];

        tr_x_xyyy_xyzz[i] = tr_x_yyy_yzz[i] * fe_0 + ts_yyy_xyzz[i] * fe_0 + tr_x_yyy_xyzz[i] * pa_x[i];

        tr_x_xyyy_xzzz[i] = 2.0 * tr_x_xy_xzzz[i] * fe_0 + tr_x_xyy_xzzz[i] * pa_y[i];

        tr_x_xyyy_yyyy[i] = ts_yyy_yyyy[i] * fe_0 + tr_x_yyy_yyyy[i] * pa_x[i];

        tr_x_xyyy_yyyz[i] = ts_yyy_yyyz[i] * fe_0 + tr_x_yyy_yyyz[i] * pa_x[i];

        tr_x_xyyy_yyzz[i] = ts_yyy_yyzz[i] * fe_0 + tr_x_yyy_yyzz[i] * pa_x[i];

        tr_x_xyyy_yzzz[i] = ts_yyy_yzzz[i] * fe_0 + tr_x_yyy_yzzz[i] * pa_x[i];

        tr_x_xyyy_zzzz[i] = ts_yyy_zzzz[i] * fe_0 + tr_x_yyy_zzzz[i] * pa_x[i];
    }

    // Set up 105-120 components of targeted buffer : GG

    auto tr_x_xyyz_xxxx = pbuffer.data(idx_dip_gg + 105);

    auto tr_x_xyyz_xxxy = pbuffer.data(idx_dip_gg + 106);

    auto tr_x_xyyz_xxxz = pbuffer.data(idx_dip_gg + 107);

    auto tr_x_xyyz_xxyy = pbuffer.data(idx_dip_gg + 108);

    auto tr_x_xyyz_xxyz = pbuffer.data(idx_dip_gg + 109);

    auto tr_x_xyyz_xxzz = pbuffer.data(idx_dip_gg + 110);

    auto tr_x_xyyz_xyyy = pbuffer.data(idx_dip_gg + 111);

    auto tr_x_xyyz_xyyz = pbuffer.data(idx_dip_gg + 112);

    auto tr_x_xyyz_xyzz = pbuffer.data(idx_dip_gg + 113);

    auto tr_x_xyyz_xzzz = pbuffer.data(idx_dip_gg + 114);

    auto tr_x_xyyz_yyyy = pbuffer.data(idx_dip_gg + 115);

    auto tr_x_xyyz_yyyz = pbuffer.data(idx_dip_gg + 116);

    auto tr_x_xyyz_yyzz = pbuffer.data(idx_dip_gg + 117);

    auto tr_x_xyyz_yzzz = pbuffer.data(idx_dip_gg + 118);

    auto tr_x_xyyz_zzzz = pbuffer.data(idx_dip_gg + 119);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_x_xyy_xxxx,  \
                             tr_x_xyy_xxxy,  \
                             tr_x_xyy_xxy,   \
                             tr_x_xyy_xxyy,  \
                             tr_x_xyy_xxyz,  \
                             tr_x_xyy_xyy,   \
                             tr_x_xyy_xyyy,  \
                             tr_x_xyy_xyyz,  \
                             tr_x_xyy_xyz,   \
                             tr_x_xyy_xyzz,  \
                             tr_x_xyy_yyyy,  \
                             tr_x_xyyz_xxxx, \
                             tr_x_xyyz_xxxy, \
                             tr_x_xyyz_xxxz, \
                             tr_x_xyyz_xxyy, \
                             tr_x_xyyz_xxyz, \
                             tr_x_xyyz_xxzz, \
                             tr_x_xyyz_xyyy, \
                             tr_x_xyyz_xyyz, \
                             tr_x_xyyz_xyzz, \
                             tr_x_xyyz_xzzz, \
                             tr_x_xyyz_yyyy, \
                             tr_x_xyyz_yyyz, \
                             tr_x_xyyz_yyzz, \
                             tr_x_xyyz_yzzz, \
                             tr_x_xyyz_zzzz, \
                             tr_x_xyz_xxxz,  \
                             tr_x_xyz_xxzz,  \
                             tr_x_xyz_xzzz,  \
                             tr_x_xz_xxxz,   \
                             tr_x_xz_xxzz,   \
                             tr_x_xz_xzzz,   \
                             tr_x_yyz_yyyz,  \
                             tr_x_yyz_yyzz,  \
                             tr_x_yyz_yzzz,  \
                             tr_x_yyz_zzzz,  \
                             ts_yyz_yyyz,    \
                             ts_yyz_yyzz,    \
                             ts_yyz_yzzz,    \
                             ts_yyz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyz_xxxx[i] = tr_x_xyy_xxxx[i] * pa_z[i];

        tr_x_xyyz_xxxy[i] = tr_x_xyy_xxxy[i] * pa_z[i];

        tr_x_xyyz_xxxz[i] = tr_x_xz_xxxz[i] * fe_0 + tr_x_xyz_xxxz[i] * pa_y[i];

        tr_x_xyyz_xxyy[i] = tr_x_xyy_xxyy[i] * pa_z[i];

        tr_x_xyyz_xxyz[i] = tr_x_xyy_xxy[i] * fe_0 + tr_x_xyy_xxyz[i] * pa_z[i];

        tr_x_xyyz_xxzz[i] = tr_x_xz_xxzz[i] * fe_0 + tr_x_xyz_xxzz[i] * pa_y[i];

        tr_x_xyyz_xyyy[i] = tr_x_xyy_xyyy[i] * pa_z[i];

        tr_x_xyyz_xyyz[i] = tr_x_xyy_xyy[i] * fe_0 + tr_x_xyy_xyyz[i] * pa_z[i];

        tr_x_xyyz_xyzz[i] = 2.0 * tr_x_xyy_xyz[i] * fe_0 + tr_x_xyy_xyzz[i] * pa_z[i];

        tr_x_xyyz_xzzz[i] = tr_x_xz_xzzz[i] * fe_0 + tr_x_xyz_xzzz[i] * pa_y[i];

        tr_x_xyyz_yyyy[i] = tr_x_xyy_yyyy[i] * pa_z[i];

        tr_x_xyyz_yyyz[i] = ts_yyz_yyyz[i] * fe_0 + tr_x_yyz_yyyz[i] * pa_x[i];

        tr_x_xyyz_yyzz[i] = ts_yyz_yyzz[i] * fe_0 + tr_x_yyz_yyzz[i] * pa_x[i];

        tr_x_xyyz_yzzz[i] = ts_yyz_yzzz[i] * fe_0 + tr_x_yyz_yzzz[i] * pa_x[i];

        tr_x_xyyz_zzzz[i] = ts_yyz_zzzz[i] * fe_0 + tr_x_yyz_zzzz[i] * pa_x[i];
    }

    // Set up 120-135 components of targeted buffer : GG

    auto tr_x_xyzz_xxxx = pbuffer.data(idx_dip_gg + 120);

    auto tr_x_xyzz_xxxy = pbuffer.data(idx_dip_gg + 121);

    auto tr_x_xyzz_xxxz = pbuffer.data(idx_dip_gg + 122);

    auto tr_x_xyzz_xxyy = pbuffer.data(idx_dip_gg + 123);

    auto tr_x_xyzz_xxyz = pbuffer.data(idx_dip_gg + 124);

    auto tr_x_xyzz_xxzz = pbuffer.data(idx_dip_gg + 125);

    auto tr_x_xyzz_xyyy = pbuffer.data(idx_dip_gg + 126);

    auto tr_x_xyzz_xyyz = pbuffer.data(idx_dip_gg + 127);

    auto tr_x_xyzz_xyzz = pbuffer.data(idx_dip_gg + 128);

    auto tr_x_xyzz_xzzz = pbuffer.data(idx_dip_gg + 129);

    auto tr_x_xyzz_yyyy = pbuffer.data(idx_dip_gg + 130);

    auto tr_x_xyzz_yyyz = pbuffer.data(idx_dip_gg + 131);

    auto tr_x_xyzz_yyzz = pbuffer.data(idx_dip_gg + 132);

    auto tr_x_xyzz_yzzz = pbuffer.data(idx_dip_gg + 133);

    auto tr_x_xyzz_zzzz = pbuffer.data(idx_dip_gg + 134);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_x_xyzz_xxxx, \
                             tr_x_xyzz_xxxy, \
                             tr_x_xyzz_xxxz, \
                             tr_x_xyzz_xxyy, \
                             tr_x_xyzz_xxyz, \
                             tr_x_xyzz_xxzz, \
                             tr_x_xyzz_xyyy, \
                             tr_x_xyzz_xyyz, \
                             tr_x_xyzz_xyzz, \
                             tr_x_xyzz_xzzz, \
                             tr_x_xyzz_yyyy, \
                             tr_x_xyzz_yyyz, \
                             tr_x_xyzz_yyzz, \
                             tr_x_xyzz_yzzz, \
                             tr_x_xyzz_zzzz, \
                             tr_x_xzz_xxx,   \
                             tr_x_xzz_xxxx,  \
                             tr_x_xzz_xxxy,  \
                             tr_x_xzz_xxxz,  \
                             tr_x_xzz_xxy,   \
                             tr_x_xzz_xxyy,  \
                             tr_x_xzz_xxyz,  \
                             tr_x_xzz_xxz,   \
                             tr_x_xzz_xxzz,  \
                             tr_x_xzz_xyy,   \
                             tr_x_xzz_xyyy,  \
                             tr_x_xzz_xyyz,  \
                             tr_x_xzz_xyz,   \
                             tr_x_xzz_xyzz,  \
                             tr_x_xzz_xzz,   \
                             tr_x_xzz_xzzz,  \
                             tr_x_xzz_zzzz,  \
                             tr_x_yzz_yyyy,  \
                             tr_x_yzz_yyyz,  \
                             tr_x_yzz_yyzz,  \
                             tr_x_yzz_yzzz,  \
                             ts_yzz_yyyy,    \
                             ts_yzz_yyyz,    \
                             ts_yzz_yyzz,    \
                             ts_yzz_yzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzz_xxxx[i] = tr_x_xzz_xxxx[i] * pa_y[i];

        tr_x_xyzz_xxxy[i] = tr_x_xzz_xxx[i] * fe_0 + tr_x_xzz_xxxy[i] * pa_y[i];

        tr_x_xyzz_xxxz[i] = tr_x_xzz_xxxz[i] * pa_y[i];

        tr_x_xyzz_xxyy[i] = 2.0 * tr_x_xzz_xxy[i] * fe_0 + tr_x_xzz_xxyy[i] * pa_y[i];

        tr_x_xyzz_xxyz[i] = tr_x_xzz_xxz[i] * fe_0 + tr_x_xzz_xxyz[i] * pa_y[i];

        tr_x_xyzz_xxzz[i] = tr_x_xzz_xxzz[i] * pa_y[i];

        tr_x_xyzz_xyyy[i] = 3.0 * tr_x_xzz_xyy[i] * fe_0 + tr_x_xzz_xyyy[i] * pa_y[i];

        tr_x_xyzz_xyyz[i] = 2.0 * tr_x_xzz_xyz[i] * fe_0 + tr_x_xzz_xyyz[i] * pa_y[i];

        tr_x_xyzz_xyzz[i] = tr_x_xzz_xzz[i] * fe_0 + tr_x_xzz_xyzz[i] * pa_y[i];

        tr_x_xyzz_xzzz[i] = tr_x_xzz_xzzz[i] * pa_y[i];

        tr_x_xyzz_yyyy[i] = ts_yzz_yyyy[i] * fe_0 + tr_x_yzz_yyyy[i] * pa_x[i];

        tr_x_xyzz_yyyz[i] = ts_yzz_yyyz[i] * fe_0 + tr_x_yzz_yyyz[i] * pa_x[i];

        tr_x_xyzz_yyzz[i] = ts_yzz_yyzz[i] * fe_0 + tr_x_yzz_yyzz[i] * pa_x[i];

        tr_x_xyzz_yzzz[i] = ts_yzz_yzzz[i] * fe_0 + tr_x_yzz_yzzz[i] * pa_x[i];

        tr_x_xyzz_zzzz[i] = tr_x_xzz_zzzz[i] * pa_y[i];
    }

    // Set up 135-150 components of targeted buffer : GG

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

    auto tr_x_xzzz_yyyy = pbuffer.data(idx_dip_gg + 145);

    auto tr_x_xzzz_yyyz = pbuffer.data(idx_dip_gg + 146);

    auto tr_x_xzzz_yyzz = pbuffer.data(idx_dip_gg + 147);

    auto tr_x_xzzz_yzzz = pbuffer.data(idx_dip_gg + 148);

    auto tr_x_xzzz_zzzz = pbuffer.data(idx_dip_gg + 149);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_x_xz_xxxx,   \
                             tr_x_xz_xxxy,   \
                             tr_x_xz_xxyy,   \
                             tr_x_xz_xyyy,   \
                             tr_x_xzz_xxxx,  \
                             tr_x_xzz_xxxy,  \
                             tr_x_xzz_xxyy,  \
                             tr_x_xzz_xyyy,  \
                             tr_x_xzzz_xxxx, \
                             tr_x_xzzz_xxxy, \
                             tr_x_xzzz_xxxz, \
                             tr_x_xzzz_xxyy, \
                             tr_x_xzzz_xxyz, \
                             tr_x_xzzz_xxzz, \
                             tr_x_xzzz_xyyy, \
                             tr_x_xzzz_xyyz, \
                             tr_x_xzzz_xyzz, \
                             tr_x_xzzz_xzzz, \
                             tr_x_xzzz_yyyy, \
                             tr_x_xzzz_yyyz, \
                             tr_x_xzzz_yyzz, \
                             tr_x_xzzz_yzzz, \
                             tr_x_xzzz_zzzz, \
                             tr_x_zzz_xxxz,  \
                             tr_x_zzz_xxyz,  \
                             tr_x_zzz_xxz,   \
                             tr_x_zzz_xxzz,  \
                             tr_x_zzz_xyyz,  \
                             tr_x_zzz_xyz,   \
                             tr_x_zzz_xyzz,  \
                             tr_x_zzz_xzz,   \
                             tr_x_zzz_xzzz,  \
                             tr_x_zzz_yyyy,  \
                             tr_x_zzz_yyyz,  \
                             tr_x_zzz_yyz,   \
                             tr_x_zzz_yyzz,  \
                             tr_x_zzz_yzz,   \
                             tr_x_zzz_yzzz,  \
                             tr_x_zzz_zzz,   \
                             tr_x_zzz_zzzz,  \
                             ts_zzz_xxxz,    \
                             ts_zzz_xxyz,    \
                             ts_zzz_xxzz,    \
                             ts_zzz_xyyz,    \
                             ts_zzz_xyzz,    \
                             ts_zzz_xzzz,    \
                             ts_zzz_yyyy,    \
                             ts_zzz_yyyz,    \
                             ts_zzz_yyzz,    \
                             ts_zzz_yzzz,    \
                             ts_zzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzz_xxxx[i] = 2.0 * tr_x_xz_xxxx[i] * fe_0 + tr_x_xzz_xxxx[i] * pa_z[i];

        tr_x_xzzz_xxxy[i] = 2.0 * tr_x_xz_xxxy[i] * fe_0 + tr_x_xzz_xxxy[i] * pa_z[i];

        tr_x_xzzz_xxxz[i] = 3.0 * tr_x_zzz_xxz[i] * fe_0 + ts_zzz_xxxz[i] * fe_0 + tr_x_zzz_xxxz[i] * pa_x[i];

        tr_x_xzzz_xxyy[i] = 2.0 * tr_x_xz_xxyy[i] * fe_0 + tr_x_xzz_xxyy[i] * pa_z[i];

        tr_x_xzzz_xxyz[i] = 2.0 * tr_x_zzz_xyz[i] * fe_0 + ts_zzz_xxyz[i] * fe_0 + tr_x_zzz_xxyz[i] * pa_x[i];

        tr_x_xzzz_xxzz[i] = 2.0 * tr_x_zzz_xzz[i] * fe_0 + ts_zzz_xxzz[i] * fe_0 + tr_x_zzz_xxzz[i] * pa_x[i];

        tr_x_xzzz_xyyy[i] = 2.0 * tr_x_xz_xyyy[i] * fe_0 + tr_x_xzz_xyyy[i] * pa_z[i];

        tr_x_xzzz_xyyz[i] = tr_x_zzz_yyz[i] * fe_0 + ts_zzz_xyyz[i] * fe_0 + tr_x_zzz_xyyz[i] * pa_x[i];

        tr_x_xzzz_xyzz[i] = tr_x_zzz_yzz[i] * fe_0 + ts_zzz_xyzz[i] * fe_0 + tr_x_zzz_xyzz[i] * pa_x[i];

        tr_x_xzzz_xzzz[i] = tr_x_zzz_zzz[i] * fe_0 + ts_zzz_xzzz[i] * fe_0 + tr_x_zzz_xzzz[i] * pa_x[i];

        tr_x_xzzz_yyyy[i] = ts_zzz_yyyy[i] * fe_0 + tr_x_zzz_yyyy[i] * pa_x[i];

        tr_x_xzzz_yyyz[i] = ts_zzz_yyyz[i] * fe_0 + tr_x_zzz_yyyz[i] * pa_x[i];

        tr_x_xzzz_yyzz[i] = ts_zzz_yyzz[i] * fe_0 + tr_x_zzz_yyzz[i] * pa_x[i];

        tr_x_xzzz_yzzz[i] = ts_zzz_yzzz[i] * fe_0 + tr_x_zzz_yzzz[i] * pa_x[i];

        tr_x_xzzz_zzzz[i] = ts_zzz_zzzz[i] * fe_0 + tr_x_zzz_zzzz[i] * pa_x[i];
    }

    // Set up 150-165 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_y,               \
                             tr_x_yy_xxxx,   \
                             tr_x_yy_xxxy,   \
                             tr_x_yy_xxxz,   \
                             tr_x_yy_xxyy,   \
                             tr_x_yy_xxyz,   \
                             tr_x_yy_xxzz,   \
                             tr_x_yy_xyyy,   \
                             tr_x_yy_xyyz,   \
                             tr_x_yy_xyzz,   \
                             tr_x_yy_xzzz,   \
                             tr_x_yy_yyyy,   \
                             tr_x_yy_yyyz,   \
                             tr_x_yy_yyzz,   \
                             tr_x_yy_yzzz,   \
                             tr_x_yy_zzzz,   \
                             tr_x_yyy_xxx,   \
                             tr_x_yyy_xxxx,  \
                             tr_x_yyy_xxxy,  \
                             tr_x_yyy_xxxz,  \
                             tr_x_yyy_xxy,   \
                             tr_x_yyy_xxyy,  \
                             tr_x_yyy_xxyz,  \
                             tr_x_yyy_xxz,   \
                             tr_x_yyy_xxzz,  \
                             tr_x_yyy_xyy,   \
                             tr_x_yyy_xyyy,  \
                             tr_x_yyy_xyyz,  \
                             tr_x_yyy_xyz,   \
                             tr_x_yyy_xyzz,  \
                             tr_x_yyy_xzz,   \
                             tr_x_yyy_xzzz,  \
                             tr_x_yyy_yyy,   \
                             tr_x_yyy_yyyy,  \
                             tr_x_yyy_yyyz,  \
                             tr_x_yyy_yyz,   \
                             tr_x_yyy_yyzz,  \
                             tr_x_yyy_yzz,   \
                             tr_x_yyy_yzzz,  \
                             tr_x_yyy_zzz,   \
                             tr_x_yyy_zzzz,  \
                             tr_x_yyyy_xxxx, \
                             tr_x_yyyy_xxxy, \
                             tr_x_yyyy_xxxz, \
                             tr_x_yyyy_xxyy, \
                             tr_x_yyyy_xxyz, \
                             tr_x_yyyy_xxzz, \
                             tr_x_yyyy_xyyy, \
                             tr_x_yyyy_xyyz, \
                             tr_x_yyyy_xyzz, \
                             tr_x_yyyy_xzzz, \
                             tr_x_yyyy_yyyy, \
                             tr_x_yyyy_yyyz, \
                             tr_x_yyyy_yyzz, \
                             tr_x_yyyy_yzzz, \
                             tr_x_yyyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyy_xxxx[i] = 3.0 * tr_x_yy_xxxx[i] * fe_0 + tr_x_yyy_xxxx[i] * pa_y[i];

        tr_x_yyyy_xxxy[i] = 3.0 * tr_x_yy_xxxy[i] * fe_0 + tr_x_yyy_xxx[i] * fe_0 + tr_x_yyy_xxxy[i] * pa_y[i];

        tr_x_yyyy_xxxz[i] = 3.0 * tr_x_yy_xxxz[i] * fe_0 + tr_x_yyy_xxxz[i] * pa_y[i];

        tr_x_yyyy_xxyy[i] = 3.0 * tr_x_yy_xxyy[i] * fe_0 + 2.0 * tr_x_yyy_xxy[i] * fe_0 + tr_x_yyy_xxyy[i] * pa_y[i];

        tr_x_yyyy_xxyz[i] = 3.0 * tr_x_yy_xxyz[i] * fe_0 + tr_x_yyy_xxz[i] * fe_0 + tr_x_yyy_xxyz[i] * pa_y[i];

        tr_x_yyyy_xxzz[i] = 3.0 * tr_x_yy_xxzz[i] * fe_0 + tr_x_yyy_xxzz[i] * pa_y[i];

        tr_x_yyyy_xyyy[i] = 3.0 * tr_x_yy_xyyy[i] * fe_0 + 3.0 * tr_x_yyy_xyy[i] * fe_0 + tr_x_yyy_xyyy[i] * pa_y[i];

        tr_x_yyyy_xyyz[i] = 3.0 * tr_x_yy_xyyz[i] * fe_0 + 2.0 * tr_x_yyy_xyz[i] * fe_0 + tr_x_yyy_xyyz[i] * pa_y[i];

        tr_x_yyyy_xyzz[i] = 3.0 * tr_x_yy_xyzz[i] * fe_0 + tr_x_yyy_xzz[i] * fe_0 + tr_x_yyy_xyzz[i] * pa_y[i];

        tr_x_yyyy_xzzz[i] = 3.0 * tr_x_yy_xzzz[i] * fe_0 + tr_x_yyy_xzzz[i] * pa_y[i];

        tr_x_yyyy_yyyy[i] = 3.0 * tr_x_yy_yyyy[i] * fe_0 + 4.0 * tr_x_yyy_yyy[i] * fe_0 + tr_x_yyy_yyyy[i] * pa_y[i];

        tr_x_yyyy_yyyz[i] = 3.0 * tr_x_yy_yyyz[i] * fe_0 + 3.0 * tr_x_yyy_yyz[i] * fe_0 + tr_x_yyy_yyyz[i] * pa_y[i];

        tr_x_yyyy_yyzz[i] = 3.0 * tr_x_yy_yyzz[i] * fe_0 + 2.0 * tr_x_yyy_yzz[i] * fe_0 + tr_x_yyy_yyzz[i] * pa_y[i];

        tr_x_yyyy_yzzz[i] = 3.0 * tr_x_yy_yzzz[i] * fe_0 + tr_x_yyy_zzz[i] * fe_0 + tr_x_yyy_yzzz[i] * pa_y[i];

        tr_x_yyyy_zzzz[i] = 3.0 * tr_x_yy_zzzz[i] * fe_0 + tr_x_yyy_zzzz[i] * pa_y[i];
    }

    // Set up 165-180 components of targeted buffer : GG

    auto tr_x_yyyz_xxxx = pbuffer.data(idx_dip_gg + 165);

    auto tr_x_yyyz_xxxy = pbuffer.data(idx_dip_gg + 166);

    auto tr_x_yyyz_xxxz = pbuffer.data(idx_dip_gg + 167);

    auto tr_x_yyyz_xxyy = pbuffer.data(idx_dip_gg + 168);

    auto tr_x_yyyz_xxyz = pbuffer.data(idx_dip_gg + 169);

    auto tr_x_yyyz_xxzz = pbuffer.data(idx_dip_gg + 170);

    auto tr_x_yyyz_xyyy = pbuffer.data(idx_dip_gg + 171);

    auto tr_x_yyyz_xyyz = pbuffer.data(idx_dip_gg + 172);

    auto tr_x_yyyz_xyzz = pbuffer.data(idx_dip_gg + 173);

    auto tr_x_yyyz_xzzz = pbuffer.data(idx_dip_gg + 174);

    auto tr_x_yyyz_yyyy = pbuffer.data(idx_dip_gg + 175);

    auto tr_x_yyyz_yyyz = pbuffer.data(idx_dip_gg + 176);

    auto tr_x_yyyz_yyzz = pbuffer.data(idx_dip_gg + 177);

    auto tr_x_yyyz_yzzz = pbuffer.data(idx_dip_gg + 178);

    auto tr_x_yyyz_zzzz = pbuffer.data(idx_dip_gg + 179);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_yyy_xxxx,  \
                             tr_x_yyy_xxxy,  \
                             tr_x_yyy_xxy,   \
                             tr_x_yyy_xxyy,  \
                             tr_x_yyy_xxyz,  \
                             tr_x_yyy_xyy,   \
                             tr_x_yyy_xyyy,  \
                             tr_x_yyy_xyyz,  \
                             tr_x_yyy_xyz,   \
                             tr_x_yyy_xyzz,  \
                             tr_x_yyy_yyy,   \
                             tr_x_yyy_yyyy,  \
                             tr_x_yyy_yyyz,  \
                             tr_x_yyy_yyz,   \
                             tr_x_yyy_yyzz,  \
                             tr_x_yyy_yzz,   \
                             tr_x_yyy_yzzz,  \
                             tr_x_yyyz_xxxx, \
                             tr_x_yyyz_xxxy, \
                             tr_x_yyyz_xxxz, \
                             tr_x_yyyz_xxyy, \
                             tr_x_yyyz_xxyz, \
                             tr_x_yyyz_xxzz, \
                             tr_x_yyyz_xyyy, \
                             tr_x_yyyz_xyyz, \
                             tr_x_yyyz_xyzz, \
                             tr_x_yyyz_xzzz, \
                             tr_x_yyyz_yyyy, \
                             tr_x_yyyz_yyyz, \
                             tr_x_yyyz_yyzz, \
                             tr_x_yyyz_yzzz, \
                             tr_x_yyyz_zzzz, \
                             tr_x_yyz_xxxz,  \
                             tr_x_yyz_xxzz,  \
                             tr_x_yyz_xzzz,  \
                             tr_x_yyz_zzzz,  \
                             tr_x_yz_xxxz,   \
                             tr_x_yz_xxzz,   \
                             tr_x_yz_xzzz,   \
                             tr_x_yz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyz_xxxx[i] = tr_x_yyy_xxxx[i] * pa_z[i];

        tr_x_yyyz_xxxy[i] = tr_x_yyy_xxxy[i] * pa_z[i];

        tr_x_yyyz_xxxz[i] = 2.0 * tr_x_yz_xxxz[i] * fe_0 + tr_x_yyz_xxxz[i] * pa_y[i];

        tr_x_yyyz_xxyy[i] = tr_x_yyy_xxyy[i] * pa_z[i];

        tr_x_yyyz_xxyz[i] = tr_x_yyy_xxy[i] * fe_0 + tr_x_yyy_xxyz[i] * pa_z[i];

        tr_x_yyyz_xxzz[i] = 2.0 * tr_x_yz_xxzz[i] * fe_0 + tr_x_yyz_xxzz[i] * pa_y[i];

        tr_x_yyyz_xyyy[i] = tr_x_yyy_xyyy[i] * pa_z[i];

        tr_x_yyyz_xyyz[i] = tr_x_yyy_xyy[i] * fe_0 + tr_x_yyy_xyyz[i] * pa_z[i];

        tr_x_yyyz_xyzz[i] = 2.0 * tr_x_yyy_xyz[i] * fe_0 + tr_x_yyy_xyzz[i] * pa_z[i];

        tr_x_yyyz_xzzz[i] = 2.0 * tr_x_yz_xzzz[i] * fe_0 + tr_x_yyz_xzzz[i] * pa_y[i];

        tr_x_yyyz_yyyy[i] = tr_x_yyy_yyyy[i] * pa_z[i];

        tr_x_yyyz_yyyz[i] = tr_x_yyy_yyy[i] * fe_0 + tr_x_yyy_yyyz[i] * pa_z[i];

        tr_x_yyyz_yyzz[i] = 2.0 * tr_x_yyy_yyz[i] * fe_0 + tr_x_yyy_yyzz[i] * pa_z[i];

        tr_x_yyyz_yzzz[i] = 3.0 * tr_x_yyy_yzz[i] * fe_0 + tr_x_yyy_yzzz[i] * pa_z[i];

        tr_x_yyyz_zzzz[i] = 2.0 * tr_x_yz_zzzz[i] * fe_0 + tr_x_yyz_zzzz[i] * pa_y[i];
    }

    // Set up 180-195 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_yy_xxxy,   \
                             tr_x_yy_xxyy,   \
                             tr_x_yy_xyyy,   \
                             tr_x_yy_yyyy,   \
                             tr_x_yyz_xxxy,  \
                             tr_x_yyz_xxyy,  \
                             tr_x_yyz_xyyy,  \
                             tr_x_yyz_yyyy,  \
                             tr_x_yyzz_xxxx, \
                             tr_x_yyzz_xxxy, \
                             tr_x_yyzz_xxxz, \
                             tr_x_yyzz_xxyy, \
                             tr_x_yyzz_xxyz, \
                             tr_x_yyzz_xxzz, \
                             tr_x_yyzz_xyyy, \
                             tr_x_yyzz_xyyz, \
                             tr_x_yyzz_xyzz, \
                             tr_x_yyzz_xzzz, \
                             tr_x_yyzz_yyyy, \
                             tr_x_yyzz_yyyz, \
                             tr_x_yyzz_yyzz, \
                             tr_x_yyzz_yzzz, \
                             tr_x_yyzz_zzzz, \
                             tr_x_yzz_xxxx,  \
                             tr_x_yzz_xxxz,  \
                             tr_x_yzz_xxyz,  \
                             tr_x_yzz_xxz,   \
                             tr_x_yzz_xxzz,  \
                             tr_x_yzz_xyyz,  \
                             tr_x_yzz_xyz,   \
                             tr_x_yzz_xyzz,  \
                             tr_x_yzz_xzz,   \
                             tr_x_yzz_xzzz,  \
                             tr_x_yzz_yyyz,  \
                             tr_x_yzz_yyz,   \
                             tr_x_yzz_yyzz,  \
                             tr_x_yzz_yzz,   \
                             tr_x_yzz_yzzz,  \
                             tr_x_yzz_zzz,   \
                             tr_x_yzz_zzzz,  \
                             tr_x_zz_xxxx,   \
                             tr_x_zz_xxxz,   \
                             tr_x_zz_xxyz,   \
                             tr_x_zz_xxzz,   \
                             tr_x_zz_xyyz,   \
                             tr_x_zz_xyzz,   \
                             tr_x_zz_xzzz,   \
                             tr_x_zz_yyyz,   \
                             tr_x_zz_yyzz,   \
                             tr_x_zz_yzzz,   \
                             tr_x_zz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzz_xxxx[i] = tr_x_zz_xxxx[i] * fe_0 + tr_x_yzz_xxxx[i] * pa_y[i];

        tr_x_yyzz_xxxy[i] = tr_x_yy_xxxy[i] * fe_0 + tr_x_yyz_xxxy[i] * pa_z[i];

        tr_x_yyzz_xxxz[i] = tr_x_zz_xxxz[i] * fe_0 + tr_x_yzz_xxxz[i] * pa_y[i];

        tr_x_yyzz_xxyy[i] = tr_x_yy_xxyy[i] * fe_0 + tr_x_yyz_xxyy[i] * pa_z[i];

        tr_x_yyzz_xxyz[i] = tr_x_zz_xxyz[i] * fe_0 + tr_x_yzz_xxz[i] * fe_0 + tr_x_yzz_xxyz[i] * pa_y[i];

        tr_x_yyzz_xxzz[i] = tr_x_zz_xxzz[i] * fe_0 + tr_x_yzz_xxzz[i] * pa_y[i];

        tr_x_yyzz_xyyy[i] = tr_x_yy_xyyy[i] * fe_0 + tr_x_yyz_xyyy[i] * pa_z[i];

        tr_x_yyzz_xyyz[i] = tr_x_zz_xyyz[i] * fe_0 + 2.0 * tr_x_yzz_xyz[i] * fe_0 + tr_x_yzz_xyyz[i] * pa_y[i];

        tr_x_yyzz_xyzz[i] = tr_x_zz_xyzz[i] * fe_0 + tr_x_yzz_xzz[i] * fe_0 + tr_x_yzz_xyzz[i] * pa_y[i];

        tr_x_yyzz_xzzz[i] = tr_x_zz_xzzz[i] * fe_0 + tr_x_yzz_xzzz[i] * pa_y[i];

        tr_x_yyzz_yyyy[i] = tr_x_yy_yyyy[i] * fe_0 + tr_x_yyz_yyyy[i] * pa_z[i];

        tr_x_yyzz_yyyz[i] = tr_x_zz_yyyz[i] * fe_0 + 3.0 * tr_x_yzz_yyz[i] * fe_0 + tr_x_yzz_yyyz[i] * pa_y[i];

        tr_x_yyzz_yyzz[i] = tr_x_zz_yyzz[i] * fe_0 + 2.0 * tr_x_yzz_yzz[i] * fe_0 + tr_x_yzz_yyzz[i] * pa_y[i];

        tr_x_yyzz_yzzz[i] = tr_x_zz_yzzz[i] * fe_0 + tr_x_yzz_zzz[i] * fe_0 + tr_x_yzz_yzzz[i] * pa_y[i];

        tr_x_yyzz_zzzz[i] = tr_x_zz_zzzz[i] * fe_0 + tr_x_yzz_zzzz[i] * pa_y[i];
    }

    // Set up 195-210 components of targeted buffer : GG

    auto tr_x_yzzz_xxxx = pbuffer.data(idx_dip_gg + 195);

    auto tr_x_yzzz_xxxy = pbuffer.data(idx_dip_gg + 196);

    auto tr_x_yzzz_xxxz = pbuffer.data(idx_dip_gg + 197);

    auto tr_x_yzzz_xxyy = pbuffer.data(idx_dip_gg + 198);

    auto tr_x_yzzz_xxyz = pbuffer.data(idx_dip_gg + 199);

    auto tr_x_yzzz_xxzz = pbuffer.data(idx_dip_gg + 200);

    auto tr_x_yzzz_xyyy = pbuffer.data(idx_dip_gg + 201);

    auto tr_x_yzzz_xyyz = pbuffer.data(idx_dip_gg + 202);

    auto tr_x_yzzz_xyzz = pbuffer.data(idx_dip_gg + 203);

    auto tr_x_yzzz_xzzz = pbuffer.data(idx_dip_gg + 204);

    auto tr_x_yzzz_yyyy = pbuffer.data(idx_dip_gg + 205);

    auto tr_x_yzzz_yyyz = pbuffer.data(idx_dip_gg + 206);

    auto tr_x_yzzz_yyzz = pbuffer.data(idx_dip_gg + 207);

    auto tr_x_yzzz_yzzz = pbuffer.data(idx_dip_gg + 208);

    auto tr_x_yzzz_zzzz = pbuffer.data(idx_dip_gg + 209);

#pragma omp simd aligned(pa_y,               \
                             tr_x_yzzz_xxxx, \
                             tr_x_yzzz_xxxy, \
                             tr_x_yzzz_xxxz, \
                             tr_x_yzzz_xxyy, \
                             tr_x_yzzz_xxyz, \
                             tr_x_yzzz_xxzz, \
                             tr_x_yzzz_xyyy, \
                             tr_x_yzzz_xyyz, \
                             tr_x_yzzz_xyzz, \
                             tr_x_yzzz_xzzz, \
                             tr_x_yzzz_yyyy, \
                             tr_x_yzzz_yyyz, \
                             tr_x_yzzz_yyzz, \
                             tr_x_yzzz_yzzz, \
                             tr_x_yzzz_zzzz, \
                             tr_x_zzz_xxx,   \
                             tr_x_zzz_xxxx,  \
                             tr_x_zzz_xxxy,  \
                             tr_x_zzz_xxxz,  \
                             tr_x_zzz_xxy,   \
                             tr_x_zzz_xxyy,  \
                             tr_x_zzz_xxyz,  \
                             tr_x_zzz_xxz,   \
                             tr_x_zzz_xxzz,  \
                             tr_x_zzz_xyy,   \
                             tr_x_zzz_xyyy,  \
                             tr_x_zzz_xyyz,  \
                             tr_x_zzz_xyz,   \
                             tr_x_zzz_xyzz,  \
                             tr_x_zzz_xzz,   \
                             tr_x_zzz_xzzz,  \
                             tr_x_zzz_yyy,   \
                             tr_x_zzz_yyyy,  \
                             tr_x_zzz_yyyz,  \
                             tr_x_zzz_yyz,   \
                             tr_x_zzz_yyzz,  \
                             tr_x_zzz_yzz,   \
                             tr_x_zzz_yzzz,  \
                             tr_x_zzz_zzz,   \
                             tr_x_zzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzz_xxxx[i] = tr_x_zzz_xxxx[i] * pa_y[i];

        tr_x_yzzz_xxxy[i] = tr_x_zzz_xxx[i] * fe_0 + tr_x_zzz_xxxy[i] * pa_y[i];

        tr_x_yzzz_xxxz[i] = tr_x_zzz_xxxz[i] * pa_y[i];

        tr_x_yzzz_xxyy[i] = 2.0 * tr_x_zzz_xxy[i] * fe_0 + tr_x_zzz_xxyy[i] * pa_y[i];

        tr_x_yzzz_xxyz[i] = tr_x_zzz_xxz[i] * fe_0 + tr_x_zzz_xxyz[i] * pa_y[i];

        tr_x_yzzz_xxzz[i] = tr_x_zzz_xxzz[i] * pa_y[i];

        tr_x_yzzz_xyyy[i] = 3.0 * tr_x_zzz_xyy[i] * fe_0 + tr_x_zzz_xyyy[i] * pa_y[i];

        tr_x_yzzz_xyyz[i] = 2.0 * tr_x_zzz_xyz[i] * fe_0 + tr_x_zzz_xyyz[i] * pa_y[i];

        tr_x_yzzz_xyzz[i] = tr_x_zzz_xzz[i] * fe_0 + tr_x_zzz_xyzz[i] * pa_y[i];

        tr_x_yzzz_xzzz[i] = tr_x_zzz_xzzz[i] * pa_y[i];

        tr_x_yzzz_yyyy[i] = 4.0 * tr_x_zzz_yyy[i] * fe_0 + tr_x_zzz_yyyy[i] * pa_y[i];

        tr_x_yzzz_yyyz[i] = 3.0 * tr_x_zzz_yyz[i] * fe_0 + tr_x_zzz_yyyz[i] * pa_y[i];

        tr_x_yzzz_yyzz[i] = 2.0 * tr_x_zzz_yzz[i] * fe_0 + tr_x_zzz_yyzz[i] * pa_y[i];

        tr_x_yzzz_yzzz[i] = tr_x_zzz_zzz[i] * fe_0 + tr_x_zzz_yzzz[i] * pa_y[i];

        tr_x_yzzz_zzzz[i] = tr_x_zzz_zzzz[i] * pa_y[i];
    }

    // Set up 210-225 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_z,               \
                             tr_x_zz_xxxx,   \
                             tr_x_zz_xxxy,   \
                             tr_x_zz_xxxz,   \
                             tr_x_zz_xxyy,   \
                             tr_x_zz_xxyz,   \
                             tr_x_zz_xxzz,   \
                             tr_x_zz_xyyy,   \
                             tr_x_zz_xyyz,   \
                             tr_x_zz_xyzz,   \
                             tr_x_zz_xzzz,   \
                             tr_x_zz_yyyy,   \
                             tr_x_zz_yyyz,   \
                             tr_x_zz_yyzz,   \
                             tr_x_zz_yzzz,   \
                             tr_x_zz_zzzz,   \
                             tr_x_zzz_xxx,   \
                             tr_x_zzz_xxxx,  \
                             tr_x_zzz_xxxy,  \
                             tr_x_zzz_xxxz,  \
                             tr_x_zzz_xxy,   \
                             tr_x_zzz_xxyy,  \
                             tr_x_zzz_xxyz,  \
                             tr_x_zzz_xxz,   \
                             tr_x_zzz_xxzz,  \
                             tr_x_zzz_xyy,   \
                             tr_x_zzz_xyyy,  \
                             tr_x_zzz_xyyz,  \
                             tr_x_zzz_xyz,   \
                             tr_x_zzz_xyzz,  \
                             tr_x_zzz_xzz,   \
                             tr_x_zzz_xzzz,  \
                             tr_x_zzz_yyy,   \
                             tr_x_zzz_yyyy,  \
                             tr_x_zzz_yyyz,  \
                             tr_x_zzz_yyz,   \
                             tr_x_zzz_yyzz,  \
                             tr_x_zzz_yzz,   \
                             tr_x_zzz_yzzz,  \
                             tr_x_zzz_zzz,   \
                             tr_x_zzz_zzzz,  \
                             tr_x_zzzz_xxxx, \
                             tr_x_zzzz_xxxy, \
                             tr_x_zzzz_xxxz, \
                             tr_x_zzzz_xxyy, \
                             tr_x_zzzz_xxyz, \
                             tr_x_zzzz_xxzz, \
                             tr_x_zzzz_xyyy, \
                             tr_x_zzzz_xyyz, \
                             tr_x_zzzz_xyzz, \
                             tr_x_zzzz_xzzz, \
                             tr_x_zzzz_yyyy, \
                             tr_x_zzzz_yyyz, \
                             tr_x_zzzz_yyzz, \
                             tr_x_zzzz_yzzz, \
                             tr_x_zzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzz_xxxx[i] = 3.0 * tr_x_zz_xxxx[i] * fe_0 + tr_x_zzz_xxxx[i] * pa_z[i];

        tr_x_zzzz_xxxy[i] = 3.0 * tr_x_zz_xxxy[i] * fe_0 + tr_x_zzz_xxxy[i] * pa_z[i];

        tr_x_zzzz_xxxz[i] = 3.0 * tr_x_zz_xxxz[i] * fe_0 + tr_x_zzz_xxx[i] * fe_0 + tr_x_zzz_xxxz[i] * pa_z[i];

        tr_x_zzzz_xxyy[i] = 3.0 * tr_x_zz_xxyy[i] * fe_0 + tr_x_zzz_xxyy[i] * pa_z[i];

        tr_x_zzzz_xxyz[i] = 3.0 * tr_x_zz_xxyz[i] * fe_0 + tr_x_zzz_xxy[i] * fe_0 + tr_x_zzz_xxyz[i] * pa_z[i];

        tr_x_zzzz_xxzz[i] = 3.0 * tr_x_zz_xxzz[i] * fe_0 + 2.0 * tr_x_zzz_xxz[i] * fe_0 + tr_x_zzz_xxzz[i] * pa_z[i];

        tr_x_zzzz_xyyy[i] = 3.0 * tr_x_zz_xyyy[i] * fe_0 + tr_x_zzz_xyyy[i] * pa_z[i];

        tr_x_zzzz_xyyz[i] = 3.0 * tr_x_zz_xyyz[i] * fe_0 + tr_x_zzz_xyy[i] * fe_0 + tr_x_zzz_xyyz[i] * pa_z[i];

        tr_x_zzzz_xyzz[i] = 3.0 * tr_x_zz_xyzz[i] * fe_0 + 2.0 * tr_x_zzz_xyz[i] * fe_0 + tr_x_zzz_xyzz[i] * pa_z[i];

        tr_x_zzzz_xzzz[i] = 3.0 * tr_x_zz_xzzz[i] * fe_0 + 3.0 * tr_x_zzz_xzz[i] * fe_0 + tr_x_zzz_xzzz[i] * pa_z[i];

        tr_x_zzzz_yyyy[i] = 3.0 * tr_x_zz_yyyy[i] * fe_0 + tr_x_zzz_yyyy[i] * pa_z[i];

        tr_x_zzzz_yyyz[i] = 3.0 * tr_x_zz_yyyz[i] * fe_0 + tr_x_zzz_yyy[i] * fe_0 + tr_x_zzz_yyyz[i] * pa_z[i];

        tr_x_zzzz_yyzz[i] = 3.0 * tr_x_zz_yyzz[i] * fe_0 + 2.0 * tr_x_zzz_yyz[i] * fe_0 + tr_x_zzz_yyzz[i] * pa_z[i];

        tr_x_zzzz_yzzz[i] = 3.0 * tr_x_zz_yzzz[i] * fe_0 + 3.0 * tr_x_zzz_yzz[i] * fe_0 + tr_x_zzz_yzzz[i] * pa_z[i];

        tr_x_zzzz_zzzz[i] = 3.0 * tr_x_zz_zzzz[i] * fe_0 + 4.0 * tr_x_zzz_zzz[i] * fe_0 + tr_x_zzz_zzzz[i] * pa_z[i];
    }

    // Set up 225-240 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,               \
                             tr_y_xx_xxxx,   \
                             tr_y_xx_xxxy,   \
                             tr_y_xx_xxxz,   \
                             tr_y_xx_xxyy,   \
                             tr_y_xx_xxyz,   \
                             tr_y_xx_xxzz,   \
                             tr_y_xx_xyyy,   \
                             tr_y_xx_xyyz,   \
                             tr_y_xx_xyzz,   \
                             tr_y_xx_xzzz,   \
                             tr_y_xx_yyyy,   \
                             tr_y_xx_yyyz,   \
                             tr_y_xx_yyzz,   \
                             tr_y_xx_yzzz,   \
                             tr_y_xx_zzzz,   \
                             tr_y_xxx_xxx,   \
                             tr_y_xxx_xxxx,  \
                             tr_y_xxx_xxxy,  \
                             tr_y_xxx_xxxz,  \
                             tr_y_xxx_xxy,   \
                             tr_y_xxx_xxyy,  \
                             tr_y_xxx_xxyz,  \
                             tr_y_xxx_xxz,   \
                             tr_y_xxx_xxzz,  \
                             tr_y_xxx_xyy,   \
                             tr_y_xxx_xyyy,  \
                             tr_y_xxx_xyyz,  \
                             tr_y_xxx_xyz,   \
                             tr_y_xxx_xyzz,  \
                             tr_y_xxx_xzz,   \
                             tr_y_xxx_xzzz,  \
                             tr_y_xxx_yyy,   \
                             tr_y_xxx_yyyy,  \
                             tr_y_xxx_yyyz,  \
                             tr_y_xxx_yyz,   \
                             tr_y_xxx_yyzz,  \
                             tr_y_xxx_yzz,   \
                             tr_y_xxx_yzzz,  \
                             tr_y_xxx_zzz,   \
                             tr_y_xxx_zzzz,  \
                             tr_y_xxxx_xxxx, \
                             tr_y_xxxx_xxxy, \
                             tr_y_xxxx_xxxz, \
                             tr_y_xxxx_xxyy, \
                             tr_y_xxxx_xxyz, \
                             tr_y_xxxx_xxzz, \
                             tr_y_xxxx_xyyy, \
                             tr_y_xxxx_xyyz, \
                             tr_y_xxxx_xyzz, \
                             tr_y_xxxx_xzzz, \
                             tr_y_xxxx_yyyy, \
                             tr_y_xxxx_yyyz, \
                             tr_y_xxxx_yyzz, \
                             tr_y_xxxx_yzzz, \
                             tr_y_xxxx_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxx_xxxx[i] = 3.0 * tr_y_xx_xxxx[i] * fe_0 + 4.0 * tr_y_xxx_xxx[i] * fe_0 + tr_y_xxx_xxxx[i] * pa_x[i];

        tr_y_xxxx_xxxy[i] = 3.0 * tr_y_xx_xxxy[i] * fe_0 + 3.0 * tr_y_xxx_xxy[i] * fe_0 + tr_y_xxx_xxxy[i] * pa_x[i];

        tr_y_xxxx_xxxz[i] = 3.0 * tr_y_xx_xxxz[i] * fe_0 + 3.0 * tr_y_xxx_xxz[i] * fe_0 + tr_y_xxx_xxxz[i] * pa_x[i];

        tr_y_xxxx_xxyy[i] = 3.0 * tr_y_xx_xxyy[i] * fe_0 + 2.0 * tr_y_xxx_xyy[i] * fe_0 + tr_y_xxx_xxyy[i] * pa_x[i];

        tr_y_xxxx_xxyz[i] = 3.0 * tr_y_xx_xxyz[i] * fe_0 + 2.0 * tr_y_xxx_xyz[i] * fe_0 + tr_y_xxx_xxyz[i] * pa_x[i];

        tr_y_xxxx_xxzz[i] = 3.0 * tr_y_xx_xxzz[i] * fe_0 + 2.0 * tr_y_xxx_xzz[i] * fe_0 + tr_y_xxx_xxzz[i] * pa_x[i];

        tr_y_xxxx_xyyy[i] = 3.0 * tr_y_xx_xyyy[i] * fe_0 + tr_y_xxx_yyy[i] * fe_0 + tr_y_xxx_xyyy[i] * pa_x[i];

        tr_y_xxxx_xyyz[i] = 3.0 * tr_y_xx_xyyz[i] * fe_0 + tr_y_xxx_yyz[i] * fe_0 + tr_y_xxx_xyyz[i] * pa_x[i];

        tr_y_xxxx_xyzz[i] = 3.0 * tr_y_xx_xyzz[i] * fe_0 + tr_y_xxx_yzz[i] * fe_0 + tr_y_xxx_xyzz[i] * pa_x[i];

        tr_y_xxxx_xzzz[i] = 3.0 * tr_y_xx_xzzz[i] * fe_0 + tr_y_xxx_zzz[i] * fe_0 + tr_y_xxx_xzzz[i] * pa_x[i];

        tr_y_xxxx_yyyy[i] = 3.0 * tr_y_xx_yyyy[i] * fe_0 + tr_y_xxx_yyyy[i] * pa_x[i];

        tr_y_xxxx_yyyz[i] = 3.0 * tr_y_xx_yyyz[i] * fe_0 + tr_y_xxx_yyyz[i] * pa_x[i];

        tr_y_xxxx_yyzz[i] = 3.0 * tr_y_xx_yyzz[i] * fe_0 + tr_y_xxx_yyzz[i] * pa_x[i];

        tr_y_xxxx_yzzz[i] = 3.0 * tr_y_xx_yzzz[i] * fe_0 + tr_y_xxx_yzzz[i] * pa_x[i];

        tr_y_xxxx_zzzz[i] = 3.0 * tr_y_xx_zzzz[i] * fe_0 + tr_y_xxx_zzzz[i] * pa_x[i];
    }

    // Set up 240-255 components of targeted buffer : GG

    auto tr_y_xxxy_xxxx = pbuffer.data(idx_dip_gg + 240);

    auto tr_y_xxxy_xxxy = pbuffer.data(idx_dip_gg + 241);

    auto tr_y_xxxy_xxxz = pbuffer.data(idx_dip_gg + 242);

    auto tr_y_xxxy_xxyy = pbuffer.data(idx_dip_gg + 243);

    auto tr_y_xxxy_xxyz = pbuffer.data(idx_dip_gg + 244);

    auto tr_y_xxxy_xxzz = pbuffer.data(idx_dip_gg + 245);

    auto tr_y_xxxy_xyyy = pbuffer.data(idx_dip_gg + 246);

    auto tr_y_xxxy_xyyz = pbuffer.data(idx_dip_gg + 247);

    auto tr_y_xxxy_xyzz = pbuffer.data(idx_dip_gg + 248);

    auto tr_y_xxxy_xzzz = pbuffer.data(idx_dip_gg + 249);

    auto tr_y_xxxy_yyyy = pbuffer.data(idx_dip_gg + 250);

    auto tr_y_xxxy_yyyz = pbuffer.data(idx_dip_gg + 251);

    auto tr_y_xxxy_yyzz = pbuffer.data(idx_dip_gg + 252);

    auto tr_y_xxxy_yzzz = pbuffer.data(idx_dip_gg + 253);

    auto tr_y_xxxy_zzzz = pbuffer.data(idx_dip_gg + 254);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_y_xxx_xxxx,  \
                             tr_y_xxx_xxxz,  \
                             tr_y_xxx_xxzz,  \
                             tr_y_xxx_xzzz,  \
                             tr_y_xxxy_xxxx, \
                             tr_y_xxxy_xxxy, \
                             tr_y_xxxy_xxxz, \
                             tr_y_xxxy_xxyy, \
                             tr_y_xxxy_xxyz, \
                             tr_y_xxxy_xxzz, \
                             tr_y_xxxy_xyyy, \
                             tr_y_xxxy_xyyz, \
                             tr_y_xxxy_xyzz, \
                             tr_y_xxxy_xzzz, \
                             tr_y_xxxy_yyyy, \
                             tr_y_xxxy_yyyz, \
                             tr_y_xxxy_yyzz, \
                             tr_y_xxxy_yzzz, \
                             tr_y_xxxy_zzzz, \
                             tr_y_xxy_xxxy,  \
                             tr_y_xxy_xxy,   \
                             tr_y_xxy_xxyy,  \
                             tr_y_xxy_xxyz,  \
                             tr_y_xxy_xyy,   \
                             tr_y_xxy_xyyy,  \
                             tr_y_xxy_xyyz,  \
                             tr_y_xxy_xyz,   \
                             tr_y_xxy_xyzz,  \
                             tr_y_xxy_yyy,   \
                             tr_y_xxy_yyyy,  \
                             tr_y_xxy_yyyz,  \
                             tr_y_xxy_yyz,   \
                             tr_y_xxy_yyzz,  \
                             tr_y_xxy_yzz,   \
                             tr_y_xxy_yzzz,  \
                             tr_y_xxy_zzzz,  \
                             tr_y_xy_xxxy,   \
                             tr_y_xy_xxyy,   \
                             tr_y_xy_xxyz,   \
                             tr_y_xy_xyyy,   \
                             tr_y_xy_xyyz,   \
                             tr_y_xy_xyzz,   \
                             tr_y_xy_yyyy,   \
                             tr_y_xy_yyyz,   \
                             tr_y_xy_yyzz,   \
                             tr_y_xy_yzzz,   \
                             tr_y_xy_zzzz,   \
                             ts_xxx_xxxx,    \
                             ts_xxx_xxxz,    \
                             ts_xxx_xxzz,    \
                             ts_xxx_xzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxy_xxxx[i] = ts_xxx_xxxx[i] * fe_0 + tr_y_xxx_xxxx[i] * pa_y[i];

        tr_y_xxxy_xxxy[i] = 2.0 * tr_y_xy_xxxy[i] * fe_0 + 3.0 * tr_y_xxy_xxy[i] * fe_0 + tr_y_xxy_xxxy[i] * pa_x[i];

        tr_y_xxxy_xxxz[i] = ts_xxx_xxxz[i] * fe_0 + tr_y_xxx_xxxz[i] * pa_y[i];

        tr_y_xxxy_xxyy[i] = 2.0 * tr_y_xy_xxyy[i] * fe_0 + 2.0 * tr_y_xxy_xyy[i] * fe_0 + tr_y_xxy_xxyy[i] * pa_x[i];

        tr_y_xxxy_xxyz[i] = 2.0 * tr_y_xy_xxyz[i] * fe_0 + 2.0 * tr_y_xxy_xyz[i] * fe_0 + tr_y_xxy_xxyz[i] * pa_x[i];

        tr_y_xxxy_xxzz[i] = ts_xxx_xxzz[i] * fe_0 + tr_y_xxx_xxzz[i] * pa_y[i];

        tr_y_xxxy_xyyy[i] = 2.0 * tr_y_xy_xyyy[i] * fe_0 + tr_y_xxy_yyy[i] * fe_0 + tr_y_xxy_xyyy[i] * pa_x[i];

        tr_y_xxxy_xyyz[i] = 2.0 * tr_y_xy_xyyz[i] * fe_0 + tr_y_xxy_yyz[i] * fe_0 + tr_y_xxy_xyyz[i] * pa_x[i];

        tr_y_xxxy_xyzz[i] = 2.0 * tr_y_xy_xyzz[i] * fe_0 + tr_y_xxy_yzz[i] * fe_0 + tr_y_xxy_xyzz[i] * pa_x[i];

        tr_y_xxxy_xzzz[i] = ts_xxx_xzzz[i] * fe_0 + tr_y_xxx_xzzz[i] * pa_y[i];

        tr_y_xxxy_yyyy[i] = 2.0 * tr_y_xy_yyyy[i] * fe_0 + tr_y_xxy_yyyy[i] * pa_x[i];

        tr_y_xxxy_yyyz[i] = 2.0 * tr_y_xy_yyyz[i] * fe_0 + tr_y_xxy_yyyz[i] * pa_x[i];

        tr_y_xxxy_yyzz[i] = 2.0 * tr_y_xy_yyzz[i] * fe_0 + tr_y_xxy_yyzz[i] * pa_x[i];

        tr_y_xxxy_yzzz[i] = 2.0 * tr_y_xy_yzzz[i] * fe_0 + tr_y_xxy_yzzz[i] * pa_x[i];

        tr_y_xxxy_zzzz[i] = 2.0 * tr_y_xy_zzzz[i] * fe_0 + tr_y_xxy_zzzz[i] * pa_x[i];
    }

    // Set up 255-270 components of targeted buffer : GG

    auto tr_y_xxxz_xxxx = pbuffer.data(idx_dip_gg + 255);

    auto tr_y_xxxz_xxxy = pbuffer.data(idx_dip_gg + 256);

    auto tr_y_xxxz_xxxz = pbuffer.data(idx_dip_gg + 257);

    auto tr_y_xxxz_xxyy = pbuffer.data(idx_dip_gg + 258);

    auto tr_y_xxxz_xxyz = pbuffer.data(idx_dip_gg + 259);

    auto tr_y_xxxz_xxzz = pbuffer.data(idx_dip_gg + 260);

    auto tr_y_xxxz_xyyy = pbuffer.data(idx_dip_gg + 261);

    auto tr_y_xxxz_xyyz = pbuffer.data(idx_dip_gg + 262);

    auto tr_y_xxxz_xyzz = pbuffer.data(idx_dip_gg + 263);

    auto tr_y_xxxz_xzzz = pbuffer.data(idx_dip_gg + 264);

    auto tr_y_xxxz_yyyy = pbuffer.data(idx_dip_gg + 265);

    auto tr_y_xxxz_yyyz = pbuffer.data(idx_dip_gg + 266);

    auto tr_y_xxxz_yyzz = pbuffer.data(idx_dip_gg + 267);

    auto tr_y_xxxz_yzzz = pbuffer.data(idx_dip_gg + 268);

    auto tr_y_xxxz_zzzz = pbuffer.data(idx_dip_gg + 269);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xxx_xxx,   \
                             tr_y_xxx_xxxx,  \
                             tr_y_xxx_xxxy,  \
                             tr_y_xxx_xxxz,  \
                             tr_y_xxx_xxy,   \
                             tr_y_xxx_xxyy,  \
                             tr_y_xxx_xxyz,  \
                             tr_y_xxx_xxz,   \
                             tr_y_xxx_xxzz,  \
                             tr_y_xxx_xyy,   \
                             tr_y_xxx_xyyy,  \
                             tr_y_xxx_xyyz,  \
                             tr_y_xxx_xyz,   \
                             tr_y_xxx_xyzz,  \
                             tr_y_xxx_xzz,   \
                             tr_y_xxx_xzzz,  \
                             tr_y_xxx_yyyy,  \
                             tr_y_xxxz_xxxx, \
                             tr_y_xxxz_xxxy, \
                             tr_y_xxxz_xxxz, \
                             tr_y_xxxz_xxyy, \
                             tr_y_xxxz_xxyz, \
                             tr_y_xxxz_xxzz, \
                             tr_y_xxxz_xyyy, \
                             tr_y_xxxz_xyyz, \
                             tr_y_xxxz_xyzz, \
                             tr_y_xxxz_xzzz, \
                             tr_y_xxxz_yyyy, \
                             tr_y_xxxz_yyyz, \
                             tr_y_xxxz_yyzz, \
                             tr_y_xxxz_yzzz, \
                             tr_y_xxxz_zzzz, \
                             tr_y_xxz_yyyz,  \
                             tr_y_xxz_yyzz,  \
                             tr_y_xxz_yzzz,  \
                             tr_y_xxz_zzzz,  \
                             tr_y_xz_yyyz,   \
                             tr_y_xz_yyzz,   \
                             tr_y_xz_yzzz,   \
                             tr_y_xz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxz_xxxx[i] = tr_y_xxx_xxxx[i] * pa_z[i];

        tr_y_xxxz_xxxy[i] = tr_y_xxx_xxxy[i] * pa_z[i];

        tr_y_xxxz_xxxz[i] = tr_y_xxx_xxx[i] * fe_0 + tr_y_xxx_xxxz[i] * pa_z[i];

        tr_y_xxxz_xxyy[i] = tr_y_xxx_xxyy[i] * pa_z[i];

        tr_y_xxxz_xxyz[i] = tr_y_xxx_xxy[i] * fe_0 + tr_y_xxx_xxyz[i] * pa_z[i];

        tr_y_xxxz_xxzz[i] = 2.0 * tr_y_xxx_xxz[i] * fe_0 + tr_y_xxx_xxzz[i] * pa_z[i];

        tr_y_xxxz_xyyy[i] = tr_y_xxx_xyyy[i] * pa_z[i];

        tr_y_xxxz_xyyz[i] = tr_y_xxx_xyy[i] * fe_0 + tr_y_xxx_xyyz[i] * pa_z[i];

        tr_y_xxxz_xyzz[i] = 2.0 * tr_y_xxx_xyz[i] * fe_0 + tr_y_xxx_xyzz[i] * pa_z[i];

        tr_y_xxxz_xzzz[i] = 3.0 * tr_y_xxx_xzz[i] * fe_0 + tr_y_xxx_xzzz[i] * pa_z[i];

        tr_y_xxxz_yyyy[i] = tr_y_xxx_yyyy[i] * pa_z[i];

        tr_y_xxxz_yyyz[i] = 2.0 * tr_y_xz_yyyz[i] * fe_0 + tr_y_xxz_yyyz[i] * pa_x[i];

        tr_y_xxxz_yyzz[i] = 2.0 * tr_y_xz_yyzz[i] * fe_0 + tr_y_xxz_yyzz[i] * pa_x[i];

        tr_y_xxxz_yzzz[i] = 2.0 * tr_y_xz_yzzz[i] * fe_0 + tr_y_xxz_yzzz[i] * pa_x[i];

        tr_y_xxxz_zzzz[i] = 2.0 * tr_y_xz_zzzz[i] * fe_0 + tr_y_xxz_zzzz[i] * pa_x[i];
    }

    // Set up 270-285 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,               \
                             tr_y_xxyy_xxxx, \
                             tr_y_xxyy_xxxy, \
                             tr_y_xxyy_xxxz, \
                             tr_y_xxyy_xxyy, \
                             tr_y_xxyy_xxyz, \
                             tr_y_xxyy_xxzz, \
                             tr_y_xxyy_xyyy, \
                             tr_y_xxyy_xyyz, \
                             tr_y_xxyy_xyzz, \
                             tr_y_xxyy_xzzz, \
                             tr_y_xxyy_yyyy, \
                             tr_y_xxyy_yyyz, \
                             tr_y_xxyy_yyzz, \
                             tr_y_xxyy_yzzz, \
                             tr_y_xxyy_zzzz, \
                             tr_y_xyy_xxx,   \
                             tr_y_xyy_xxxx,  \
                             tr_y_xyy_xxxy,  \
                             tr_y_xyy_xxxz,  \
                             tr_y_xyy_xxy,   \
                             tr_y_xyy_xxyy,  \
                             tr_y_xyy_xxyz,  \
                             tr_y_xyy_xxz,   \
                             tr_y_xyy_xxzz,  \
                             tr_y_xyy_xyy,   \
                             tr_y_xyy_xyyy,  \
                             tr_y_xyy_xyyz,  \
                             tr_y_xyy_xyz,   \
                             tr_y_xyy_xyzz,  \
                             tr_y_xyy_xzz,   \
                             tr_y_xyy_xzzz,  \
                             tr_y_xyy_yyy,   \
                             tr_y_xyy_yyyy,  \
                             tr_y_xyy_yyyz,  \
                             tr_y_xyy_yyz,   \
                             tr_y_xyy_yyzz,  \
                             tr_y_xyy_yzz,   \
                             tr_y_xyy_yzzz,  \
                             tr_y_xyy_zzz,   \
                             tr_y_xyy_zzzz,  \
                             tr_y_yy_xxxx,   \
                             tr_y_yy_xxxy,   \
                             tr_y_yy_xxxz,   \
                             tr_y_yy_xxyy,   \
                             tr_y_yy_xxyz,   \
                             tr_y_yy_xxzz,   \
                             tr_y_yy_xyyy,   \
                             tr_y_yy_xyyz,   \
                             tr_y_yy_xyzz,   \
                             tr_y_yy_xzzz,   \
                             tr_y_yy_yyyy,   \
                             tr_y_yy_yyyz,   \
                             tr_y_yy_yyzz,   \
                             tr_y_yy_yzzz,   \
                             tr_y_yy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyy_xxxx[i] = tr_y_yy_xxxx[i] * fe_0 + 4.0 * tr_y_xyy_xxx[i] * fe_0 + tr_y_xyy_xxxx[i] * pa_x[i];

        tr_y_xxyy_xxxy[i] = tr_y_yy_xxxy[i] * fe_0 + 3.0 * tr_y_xyy_xxy[i] * fe_0 + tr_y_xyy_xxxy[i] * pa_x[i];

        tr_y_xxyy_xxxz[i] = tr_y_yy_xxxz[i] * fe_0 + 3.0 * tr_y_xyy_xxz[i] * fe_0 + tr_y_xyy_xxxz[i] * pa_x[i];

        tr_y_xxyy_xxyy[i] = tr_y_yy_xxyy[i] * fe_0 + 2.0 * tr_y_xyy_xyy[i] * fe_0 + tr_y_xyy_xxyy[i] * pa_x[i];

        tr_y_xxyy_xxyz[i] = tr_y_yy_xxyz[i] * fe_0 + 2.0 * tr_y_xyy_xyz[i] * fe_0 + tr_y_xyy_xxyz[i] * pa_x[i];

        tr_y_xxyy_xxzz[i] = tr_y_yy_xxzz[i] * fe_0 + 2.0 * tr_y_xyy_xzz[i] * fe_0 + tr_y_xyy_xxzz[i] * pa_x[i];

        tr_y_xxyy_xyyy[i] = tr_y_yy_xyyy[i] * fe_0 + tr_y_xyy_yyy[i] * fe_0 + tr_y_xyy_xyyy[i] * pa_x[i];

        tr_y_xxyy_xyyz[i] = tr_y_yy_xyyz[i] * fe_0 + tr_y_xyy_yyz[i] * fe_0 + tr_y_xyy_xyyz[i] * pa_x[i];

        tr_y_xxyy_xyzz[i] = tr_y_yy_xyzz[i] * fe_0 + tr_y_xyy_yzz[i] * fe_0 + tr_y_xyy_xyzz[i] * pa_x[i];

        tr_y_xxyy_xzzz[i] = tr_y_yy_xzzz[i] * fe_0 + tr_y_xyy_zzz[i] * fe_0 + tr_y_xyy_xzzz[i] * pa_x[i];

        tr_y_xxyy_yyyy[i] = tr_y_yy_yyyy[i] * fe_0 + tr_y_xyy_yyyy[i] * pa_x[i];

        tr_y_xxyy_yyyz[i] = tr_y_yy_yyyz[i] * fe_0 + tr_y_xyy_yyyz[i] * pa_x[i];

        tr_y_xxyy_yyzz[i] = tr_y_yy_yyzz[i] * fe_0 + tr_y_xyy_yyzz[i] * pa_x[i];

        tr_y_xxyy_yzzz[i] = tr_y_yy_yzzz[i] * fe_0 + tr_y_xyy_yzzz[i] * pa_x[i];

        tr_y_xxyy_zzzz[i] = tr_y_yy_zzzz[i] * fe_0 + tr_y_xyy_zzzz[i] * pa_x[i];
    }

    // Set up 285-300 components of targeted buffer : GG

    auto tr_y_xxyz_xxxx = pbuffer.data(idx_dip_gg + 285);

    auto tr_y_xxyz_xxxy = pbuffer.data(idx_dip_gg + 286);

    auto tr_y_xxyz_xxxz = pbuffer.data(idx_dip_gg + 287);

    auto tr_y_xxyz_xxyy = pbuffer.data(idx_dip_gg + 288);

    auto tr_y_xxyz_xxyz = pbuffer.data(idx_dip_gg + 289);

    auto tr_y_xxyz_xxzz = pbuffer.data(idx_dip_gg + 290);

    auto tr_y_xxyz_xyyy = pbuffer.data(idx_dip_gg + 291);

    auto tr_y_xxyz_xyyz = pbuffer.data(idx_dip_gg + 292);

    auto tr_y_xxyz_xyzz = pbuffer.data(idx_dip_gg + 293);

    auto tr_y_xxyz_xzzz = pbuffer.data(idx_dip_gg + 294);

    auto tr_y_xxyz_yyyy = pbuffer.data(idx_dip_gg + 295);

    auto tr_y_xxyz_yyyz = pbuffer.data(idx_dip_gg + 296);

    auto tr_y_xxyz_yyzz = pbuffer.data(idx_dip_gg + 297);

    auto tr_y_xxyz_yzzz = pbuffer.data(idx_dip_gg + 298);

    auto tr_y_xxyz_zzzz = pbuffer.data(idx_dip_gg + 299);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_y_xxy_xxxx,  \
                             tr_y_xxy_xxxy,  \
                             tr_y_xxy_xxy,   \
                             tr_y_xxy_xxyy,  \
                             tr_y_xxy_xxyz,  \
                             tr_y_xxy_xyy,   \
                             tr_y_xxy_xyyy,  \
                             tr_y_xxy_xyyz,  \
                             tr_y_xxy_xyz,   \
                             tr_y_xxy_xyzz,  \
                             tr_y_xxy_yyyy,  \
                             tr_y_xxyz_xxxx, \
                             tr_y_xxyz_xxxy, \
                             tr_y_xxyz_xxxz, \
                             tr_y_xxyz_xxyy, \
                             tr_y_xxyz_xxyz, \
                             tr_y_xxyz_xxzz, \
                             tr_y_xxyz_xyyy, \
                             tr_y_xxyz_xyyz, \
                             tr_y_xxyz_xyzz, \
                             tr_y_xxyz_xzzz, \
                             tr_y_xxyz_yyyy, \
                             tr_y_xxyz_yyyz, \
                             tr_y_xxyz_yyzz, \
                             tr_y_xxyz_yzzz, \
                             tr_y_xxyz_zzzz, \
                             tr_y_xxz_xxxz,  \
                             tr_y_xxz_xxzz,  \
                             tr_y_xxz_xzzz,  \
                             tr_y_xyz_yyyz,  \
                             tr_y_xyz_yyzz,  \
                             tr_y_xyz_yzzz,  \
                             tr_y_xyz_zzzz,  \
                             tr_y_yz_yyyz,   \
                             tr_y_yz_yyzz,   \
                             tr_y_yz_yzzz,   \
                             tr_y_yz_zzzz,   \
                             ts_xxz_xxxz,    \
                             ts_xxz_xxzz,    \
                             ts_xxz_xzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyz_xxxx[i] = tr_y_xxy_xxxx[i] * pa_z[i];

        tr_y_xxyz_xxxy[i] = tr_y_xxy_xxxy[i] * pa_z[i];

        tr_y_xxyz_xxxz[i] = ts_xxz_xxxz[i] * fe_0 + tr_y_xxz_xxxz[i] * pa_y[i];

        tr_y_xxyz_xxyy[i] = tr_y_xxy_xxyy[i] * pa_z[i];

        tr_y_xxyz_xxyz[i] = tr_y_xxy_xxy[i] * fe_0 + tr_y_xxy_xxyz[i] * pa_z[i];

        tr_y_xxyz_xxzz[i] = ts_xxz_xxzz[i] * fe_0 + tr_y_xxz_xxzz[i] * pa_y[i];

        tr_y_xxyz_xyyy[i] = tr_y_xxy_xyyy[i] * pa_z[i];

        tr_y_xxyz_xyyz[i] = tr_y_xxy_xyy[i] * fe_0 + tr_y_xxy_xyyz[i] * pa_z[i];

        tr_y_xxyz_xyzz[i] = 2.0 * tr_y_xxy_xyz[i] * fe_0 + tr_y_xxy_xyzz[i] * pa_z[i];

        tr_y_xxyz_xzzz[i] = ts_xxz_xzzz[i] * fe_0 + tr_y_xxz_xzzz[i] * pa_y[i];

        tr_y_xxyz_yyyy[i] = tr_y_xxy_yyyy[i] * pa_z[i];

        tr_y_xxyz_yyyz[i] = tr_y_yz_yyyz[i] * fe_0 + tr_y_xyz_yyyz[i] * pa_x[i];

        tr_y_xxyz_yyzz[i] = tr_y_yz_yyzz[i] * fe_0 + tr_y_xyz_yyzz[i] * pa_x[i];

        tr_y_xxyz_yzzz[i] = tr_y_yz_yzzz[i] * fe_0 + tr_y_xyz_yzzz[i] * pa_x[i];

        tr_y_xxyz_zzzz[i] = tr_y_yz_zzzz[i] * fe_0 + tr_y_xyz_zzzz[i] * pa_x[i];
    }

    // Set up 300-315 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xx_xxxx,   \
                             tr_y_xx_xxxy,   \
                             tr_y_xx_xxyy,   \
                             tr_y_xx_xyyy,   \
                             tr_y_xxz_xxxx,  \
                             tr_y_xxz_xxxy,  \
                             tr_y_xxz_xxyy,  \
                             tr_y_xxz_xyyy,  \
                             tr_y_xxzz_xxxx, \
                             tr_y_xxzz_xxxy, \
                             tr_y_xxzz_xxxz, \
                             tr_y_xxzz_xxyy, \
                             tr_y_xxzz_xxyz, \
                             tr_y_xxzz_xxzz, \
                             tr_y_xxzz_xyyy, \
                             tr_y_xxzz_xyyz, \
                             tr_y_xxzz_xyzz, \
                             tr_y_xxzz_xzzz, \
                             tr_y_xxzz_yyyy, \
                             tr_y_xxzz_yyyz, \
                             tr_y_xxzz_yyzz, \
                             tr_y_xxzz_yzzz, \
                             tr_y_xxzz_zzzz, \
                             tr_y_xzz_xxxz,  \
                             tr_y_xzz_xxyz,  \
                             tr_y_xzz_xxz,   \
                             tr_y_xzz_xxzz,  \
                             tr_y_xzz_xyyz,  \
                             tr_y_xzz_xyz,   \
                             tr_y_xzz_xyzz,  \
                             tr_y_xzz_xzz,   \
                             tr_y_xzz_xzzz,  \
                             tr_y_xzz_yyyy,  \
                             tr_y_xzz_yyyz,  \
                             tr_y_xzz_yyz,   \
                             tr_y_xzz_yyzz,  \
                             tr_y_xzz_yzz,   \
                             tr_y_xzz_yzzz,  \
                             tr_y_xzz_zzz,   \
                             tr_y_xzz_zzzz,  \
                             tr_y_zz_xxxz,   \
                             tr_y_zz_xxyz,   \
                             tr_y_zz_xxzz,   \
                             tr_y_zz_xyyz,   \
                             tr_y_zz_xyzz,   \
                             tr_y_zz_xzzz,   \
                             tr_y_zz_yyyy,   \
                             tr_y_zz_yyyz,   \
                             tr_y_zz_yyzz,   \
                             tr_y_zz_yzzz,   \
                             tr_y_zz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzz_xxxx[i] = tr_y_xx_xxxx[i] * fe_0 + tr_y_xxz_xxxx[i] * pa_z[i];

        tr_y_xxzz_xxxy[i] = tr_y_xx_xxxy[i] * fe_0 + tr_y_xxz_xxxy[i] * pa_z[i];

        tr_y_xxzz_xxxz[i] = tr_y_zz_xxxz[i] * fe_0 + 3.0 * tr_y_xzz_xxz[i] * fe_0 + tr_y_xzz_xxxz[i] * pa_x[i];

        tr_y_xxzz_xxyy[i] = tr_y_xx_xxyy[i] * fe_0 + tr_y_xxz_xxyy[i] * pa_z[i];

        tr_y_xxzz_xxyz[i] = tr_y_zz_xxyz[i] * fe_0 + 2.0 * tr_y_xzz_xyz[i] * fe_0 + tr_y_xzz_xxyz[i] * pa_x[i];

        tr_y_xxzz_xxzz[i] = tr_y_zz_xxzz[i] * fe_0 + 2.0 * tr_y_xzz_xzz[i] * fe_0 + tr_y_xzz_xxzz[i] * pa_x[i];

        tr_y_xxzz_xyyy[i] = tr_y_xx_xyyy[i] * fe_0 + tr_y_xxz_xyyy[i] * pa_z[i];

        tr_y_xxzz_xyyz[i] = tr_y_zz_xyyz[i] * fe_0 + tr_y_xzz_yyz[i] * fe_0 + tr_y_xzz_xyyz[i] * pa_x[i];

        tr_y_xxzz_xyzz[i] = tr_y_zz_xyzz[i] * fe_0 + tr_y_xzz_yzz[i] * fe_0 + tr_y_xzz_xyzz[i] * pa_x[i];

        tr_y_xxzz_xzzz[i] = tr_y_zz_xzzz[i] * fe_0 + tr_y_xzz_zzz[i] * fe_0 + tr_y_xzz_xzzz[i] * pa_x[i];

        tr_y_xxzz_yyyy[i] = tr_y_zz_yyyy[i] * fe_0 + tr_y_xzz_yyyy[i] * pa_x[i];

        tr_y_xxzz_yyyz[i] = tr_y_zz_yyyz[i] * fe_0 + tr_y_xzz_yyyz[i] * pa_x[i];

        tr_y_xxzz_yyzz[i] = tr_y_zz_yyzz[i] * fe_0 + tr_y_xzz_yyzz[i] * pa_x[i];

        tr_y_xxzz_yzzz[i] = tr_y_zz_yzzz[i] * fe_0 + tr_y_xzz_yzzz[i] * pa_x[i];

        tr_y_xxzz_zzzz[i] = tr_y_zz_zzzz[i] * fe_0 + tr_y_xzz_zzzz[i] * pa_x[i];
    }

    // Set up 315-330 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,               \
                             tr_y_xyyy_xxxx, \
                             tr_y_xyyy_xxxy, \
                             tr_y_xyyy_xxxz, \
                             tr_y_xyyy_xxyy, \
                             tr_y_xyyy_xxyz, \
                             tr_y_xyyy_xxzz, \
                             tr_y_xyyy_xyyy, \
                             tr_y_xyyy_xyyz, \
                             tr_y_xyyy_xyzz, \
                             tr_y_xyyy_xzzz, \
                             tr_y_xyyy_yyyy, \
                             tr_y_xyyy_yyyz, \
                             tr_y_xyyy_yyzz, \
                             tr_y_xyyy_yzzz, \
                             tr_y_xyyy_zzzz, \
                             tr_y_yyy_xxx,   \
                             tr_y_yyy_xxxx,  \
                             tr_y_yyy_xxxy,  \
                             tr_y_yyy_xxxz,  \
                             tr_y_yyy_xxy,   \
                             tr_y_yyy_xxyy,  \
                             tr_y_yyy_xxyz,  \
                             tr_y_yyy_xxz,   \
                             tr_y_yyy_xxzz,  \
                             tr_y_yyy_xyy,   \
                             tr_y_yyy_xyyy,  \
                             tr_y_yyy_xyyz,  \
                             tr_y_yyy_xyz,   \
                             tr_y_yyy_xyzz,  \
                             tr_y_yyy_xzz,   \
                             tr_y_yyy_xzzz,  \
                             tr_y_yyy_yyy,   \
                             tr_y_yyy_yyyy,  \
                             tr_y_yyy_yyyz,  \
                             tr_y_yyy_yyz,   \
                             tr_y_yyy_yyzz,  \
                             tr_y_yyy_yzz,   \
                             tr_y_yyy_yzzz,  \
                             tr_y_yyy_zzz,   \
                             tr_y_yyy_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyy_xxxx[i] = 4.0 * tr_y_yyy_xxx[i] * fe_0 + tr_y_yyy_xxxx[i] * pa_x[i];

        tr_y_xyyy_xxxy[i] = 3.0 * tr_y_yyy_xxy[i] * fe_0 + tr_y_yyy_xxxy[i] * pa_x[i];

        tr_y_xyyy_xxxz[i] = 3.0 * tr_y_yyy_xxz[i] * fe_0 + tr_y_yyy_xxxz[i] * pa_x[i];

        tr_y_xyyy_xxyy[i] = 2.0 * tr_y_yyy_xyy[i] * fe_0 + tr_y_yyy_xxyy[i] * pa_x[i];

        tr_y_xyyy_xxyz[i] = 2.0 * tr_y_yyy_xyz[i] * fe_0 + tr_y_yyy_xxyz[i] * pa_x[i];

        tr_y_xyyy_xxzz[i] = 2.0 * tr_y_yyy_xzz[i] * fe_0 + tr_y_yyy_xxzz[i] * pa_x[i];

        tr_y_xyyy_xyyy[i] = tr_y_yyy_yyy[i] * fe_0 + tr_y_yyy_xyyy[i] * pa_x[i];

        tr_y_xyyy_xyyz[i] = tr_y_yyy_yyz[i] * fe_0 + tr_y_yyy_xyyz[i] * pa_x[i];

        tr_y_xyyy_xyzz[i] = tr_y_yyy_yzz[i] * fe_0 + tr_y_yyy_xyzz[i] * pa_x[i];

        tr_y_xyyy_xzzz[i] = tr_y_yyy_zzz[i] * fe_0 + tr_y_yyy_xzzz[i] * pa_x[i];

        tr_y_xyyy_yyyy[i] = tr_y_yyy_yyyy[i] * pa_x[i];

        tr_y_xyyy_yyyz[i] = tr_y_yyy_yyyz[i] * pa_x[i];

        tr_y_xyyy_yyzz[i] = tr_y_yyy_yyzz[i] * pa_x[i];

        tr_y_xyyy_yzzz[i] = tr_y_yyy_yzzz[i] * pa_x[i];

        tr_y_xyyy_zzzz[i] = tr_y_yyy_zzzz[i] * pa_x[i];
    }

    // Set up 330-345 components of targeted buffer : GG

    auto tr_y_xyyz_xxxx = pbuffer.data(idx_dip_gg + 330);

    auto tr_y_xyyz_xxxy = pbuffer.data(idx_dip_gg + 331);

    auto tr_y_xyyz_xxxz = pbuffer.data(idx_dip_gg + 332);

    auto tr_y_xyyz_xxyy = pbuffer.data(idx_dip_gg + 333);

    auto tr_y_xyyz_xxyz = pbuffer.data(idx_dip_gg + 334);

    auto tr_y_xyyz_xxzz = pbuffer.data(idx_dip_gg + 335);

    auto tr_y_xyyz_xyyy = pbuffer.data(idx_dip_gg + 336);

    auto tr_y_xyyz_xyyz = pbuffer.data(idx_dip_gg + 337);

    auto tr_y_xyyz_xyzz = pbuffer.data(idx_dip_gg + 338);

    auto tr_y_xyyz_xzzz = pbuffer.data(idx_dip_gg + 339);

    auto tr_y_xyyz_yyyy = pbuffer.data(idx_dip_gg + 340);

    auto tr_y_xyyz_yyyz = pbuffer.data(idx_dip_gg + 341);

    auto tr_y_xyyz_yyzz = pbuffer.data(idx_dip_gg + 342);

    auto tr_y_xyyz_yzzz = pbuffer.data(idx_dip_gg + 343);

    auto tr_y_xyyz_zzzz = pbuffer.data(idx_dip_gg + 344);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xyy_xxxx,  \
                             tr_y_xyy_xxxy,  \
                             tr_y_xyy_xxyy,  \
                             tr_y_xyy_xyyy,  \
                             tr_y_xyyz_xxxx, \
                             tr_y_xyyz_xxxy, \
                             tr_y_xyyz_xxxz, \
                             tr_y_xyyz_xxyy, \
                             tr_y_xyyz_xxyz, \
                             tr_y_xyyz_xxzz, \
                             tr_y_xyyz_xyyy, \
                             tr_y_xyyz_xyyz, \
                             tr_y_xyyz_xyzz, \
                             tr_y_xyyz_xzzz, \
                             tr_y_xyyz_yyyy, \
                             tr_y_xyyz_yyyz, \
                             tr_y_xyyz_yyzz, \
                             tr_y_xyyz_yzzz, \
                             tr_y_xyyz_zzzz, \
                             tr_y_yyz_xxxz,  \
                             tr_y_yyz_xxyz,  \
                             tr_y_yyz_xxz,   \
                             tr_y_yyz_xxzz,  \
                             tr_y_yyz_xyyz,  \
                             tr_y_yyz_xyz,   \
                             tr_y_yyz_xyzz,  \
                             tr_y_yyz_xzz,   \
                             tr_y_yyz_xzzz,  \
                             tr_y_yyz_yyyy,  \
                             tr_y_yyz_yyyz,  \
                             tr_y_yyz_yyz,   \
                             tr_y_yyz_yyzz,  \
                             tr_y_yyz_yzz,   \
                             tr_y_yyz_yzzz,  \
                             tr_y_yyz_zzz,   \
                             tr_y_yyz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyz_xxxx[i] = tr_y_xyy_xxxx[i] * pa_z[i];

        tr_y_xyyz_xxxy[i] = tr_y_xyy_xxxy[i] * pa_z[i];

        tr_y_xyyz_xxxz[i] = 3.0 * tr_y_yyz_xxz[i] * fe_0 + tr_y_yyz_xxxz[i] * pa_x[i];

        tr_y_xyyz_xxyy[i] = tr_y_xyy_xxyy[i] * pa_z[i];

        tr_y_xyyz_xxyz[i] = 2.0 * tr_y_yyz_xyz[i] * fe_0 + tr_y_yyz_xxyz[i] * pa_x[i];

        tr_y_xyyz_xxzz[i] = 2.0 * tr_y_yyz_xzz[i] * fe_0 + tr_y_yyz_xxzz[i] * pa_x[i];

        tr_y_xyyz_xyyy[i] = tr_y_xyy_xyyy[i] * pa_z[i];

        tr_y_xyyz_xyyz[i] = tr_y_yyz_yyz[i] * fe_0 + tr_y_yyz_xyyz[i] * pa_x[i];

        tr_y_xyyz_xyzz[i] = tr_y_yyz_yzz[i] * fe_0 + tr_y_yyz_xyzz[i] * pa_x[i];

        tr_y_xyyz_xzzz[i] = tr_y_yyz_zzz[i] * fe_0 + tr_y_yyz_xzzz[i] * pa_x[i];

        tr_y_xyyz_yyyy[i] = tr_y_yyz_yyyy[i] * pa_x[i];

        tr_y_xyyz_yyyz[i] = tr_y_yyz_yyyz[i] * pa_x[i];

        tr_y_xyyz_yyzz[i] = tr_y_yyz_yyzz[i] * pa_x[i];

        tr_y_xyyz_yzzz[i] = tr_y_yyz_yzzz[i] * pa_x[i];

        tr_y_xyyz_zzzz[i] = tr_y_yyz_zzzz[i] * pa_x[i];
    }

    // Set up 345-360 components of targeted buffer : GG

    auto tr_y_xyzz_xxxx = pbuffer.data(idx_dip_gg + 345);

    auto tr_y_xyzz_xxxy = pbuffer.data(idx_dip_gg + 346);

    auto tr_y_xyzz_xxxz = pbuffer.data(idx_dip_gg + 347);

    auto tr_y_xyzz_xxyy = pbuffer.data(idx_dip_gg + 348);

    auto tr_y_xyzz_xxyz = pbuffer.data(idx_dip_gg + 349);

    auto tr_y_xyzz_xxzz = pbuffer.data(idx_dip_gg + 350);

    auto tr_y_xyzz_xyyy = pbuffer.data(idx_dip_gg + 351);

    auto tr_y_xyzz_xyyz = pbuffer.data(idx_dip_gg + 352);

    auto tr_y_xyzz_xyzz = pbuffer.data(idx_dip_gg + 353);

    auto tr_y_xyzz_xzzz = pbuffer.data(idx_dip_gg + 354);

    auto tr_y_xyzz_yyyy = pbuffer.data(idx_dip_gg + 355);

    auto tr_y_xyzz_yyyz = pbuffer.data(idx_dip_gg + 356);

    auto tr_y_xyzz_yyzz = pbuffer.data(idx_dip_gg + 357);

    auto tr_y_xyzz_yzzz = pbuffer.data(idx_dip_gg + 358);

    auto tr_y_xyzz_zzzz = pbuffer.data(idx_dip_gg + 359);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xyzz_xxxx, \
                             tr_y_xyzz_xxxy, \
                             tr_y_xyzz_xxxz, \
                             tr_y_xyzz_xxyy, \
                             tr_y_xyzz_xxyz, \
                             tr_y_xyzz_xxzz, \
                             tr_y_xyzz_xyyy, \
                             tr_y_xyzz_xyyz, \
                             tr_y_xyzz_xyzz, \
                             tr_y_xyzz_xzzz, \
                             tr_y_xyzz_yyyy, \
                             tr_y_xyzz_yyyz, \
                             tr_y_xyzz_yyzz, \
                             tr_y_xyzz_yzzz, \
                             tr_y_xyzz_zzzz, \
                             tr_y_yzz_xxx,   \
                             tr_y_yzz_xxxx,  \
                             tr_y_yzz_xxxy,  \
                             tr_y_yzz_xxxz,  \
                             tr_y_yzz_xxy,   \
                             tr_y_yzz_xxyy,  \
                             tr_y_yzz_xxyz,  \
                             tr_y_yzz_xxz,   \
                             tr_y_yzz_xxzz,  \
                             tr_y_yzz_xyy,   \
                             tr_y_yzz_xyyy,  \
                             tr_y_yzz_xyyz,  \
                             tr_y_yzz_xyz,   \
                             tr_y_yzz_xyzz,  \
                             tr_y_yzz_xzz,   \
                             tr_y_yzz_xzzz,  \
                             tr_y_yzz_yyy,   \
                             tr_y_yzz_yyyy,  \
                             tr_y_yzz_yyyz,  \
                             tr_y_yzz_yyz,   \
                             tr_y_yzz_yyzz,  \
                             tr_y_yzz_yzz,   \
                             tr_y_yzz_yzzz,  \
                             tr_y_yzz_zzz,   \
                             tr_y_yzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzz_xxxx[i] = 4.0 * tr_y_yzz_xxx[i] * fe_0 + tr_y_yzz_xxxx[i] * pa_x[i];

        tr_y_xyzz_xxxy[i] = 3.0 * tr_y_yzz_xxy[i] * fe_0 + tr_y_yzz_xxxy[i] * pa_x[i];

        tr_y_xyzz_xxxz[i] = 3.0 * tr_y_yzz_xxz[i] * fe_0 + tr_y_yzz_xxxz[i] * pa_x[i];

        tr_y_xyzz_xxyy[i] = 2.0 * tr_y_yzz_xyy[i] * fe_0 + tr_y_yzz_xxyy[i] * pa_x[i];

        tr_y_xyzz_xxyz[i] = 2.0 * tr_y_yzz_xyz[i] * fe_0 + tr_y_yzz_xxyz[i] * pa_x[i];

        tr_y_xyzz_xxzz[i] = 2.0 * tr_y_yzz_xzz[i] * fe_0 + tr_y_yzz_xxzz[i] * pa_x[i];

        tr_y_xyzz_xyyy[i] = tr_y_yzz_yyy[i] * fe_0 + tr_y_yzz_xyyy[i] * pa_x[i];

        tr_y_xyzz_xyyz[i] = tr_y_yzz_yyz[i] * fe_0 + tr_y_yzz_xyyz[i] * pa_x[i];

        tr_y_xyzz_xyzz[i] = tr_y_yzz_yzz[i] * fe_0 + tr_y_yzz_xyzz[i] * pa_x[i];

        tr_y_xyzz_xzzz[i] = tr_y_yzz_zzz[i] * fe_0 + tr_y_yzz_xzzz[i] * pa_x[i];

        tr_y_xyzz_yyyy[i] = tr_y_yzz_yyyy[i] * pa_x[i];

        tr_y_xyzz_yyyz[i] = tr_y_yzz_yyyz[i] * pa_x[i];

        tr_y_xyzz_yyzz[i] = tr_y_yzz_yyzz[i] * pa_x[i];

        tr_y_xyzz_yzzz[i] = tr_y_yzz_yzzz[i] * pa_x[i];

        tr_y_xyzz_zzzz[i] = tr_y_yzz_zzzz[i] * pa_x[i];
    }

    // Set up 360-375 components of targeted buffer : GG

    auto tr_y_xzzz_xxxx = pbuffer.data(idx_dip_gg + 360);

    auto tr_y_xzzz_xxxy = pbuffer.data(idx_dip_gg + 361);

    auto tr_y_xzzz_xxxz = pbuffer.data(idx_dip_gg + 362);

    auto tr_y_xzzz_xxyy = pbuffer.data(idx_dip_gg + 363);

    auto tr_y_xzzz_xxyz = pbuffer.data(idx_dip_gg + 364);

    auto tr_y_xzzz_xxzz = pbuffer.data(idx_dip_gg + 365);

    auto tr_y_xzzz_xyyy = pbuffer.data(idx_dip_gg + 366);

    auto tr_y_xzzz_xyyz = pbuffer.data(idx_dip_gg + 367);

    auto tr_y_xzzz_xyzz = pbuffer.data(idx_dip_gg + 368);

    auto tr_y_xzzz_xzzz = pbuffer.data(idx_dip_gg + 369);

    auto tr_y_xzzz_yyyy = pbuffer.data(idx_dip_gg + 370);

    auto tr_y_xzzz_yyyz = pbuffer.data(idx_dip_gg + 371);

    auto tr_y_xzzz_yyzz = pbuffer.data(idx_dip_gg + 372);

    auto tr_y_xzzz_yzzz = pbuffer.data(idx_dip_gg + 373);

    auto tr_y_xzzz_zzzz = pbuffer.data(idx_dip_gg + 374);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xzzz_xxxx, \
                             tr_y_xzzz_xxxy, \
                             tr_y_xzzz_xxxz, \
                             tr_y_xzzz_xxyy, \
                             tr_y_xzzz_xxyz, \
                             tr_y_xzzz_xxzz, \
                             tr_y_xzzz_xyyy, \
                             tr_y_xzzz_xyyz, \
                             tr_y_xzzz_xyzz, \
                             tr_y_xzzz_xzzz, \
                             tr_y_xzzz_yyyy, \
                             tr_y_xzzz_yyyz, \
                             tr_y_xzzz_yyzz, \
                             tr_y_xzzz_yzzz, \
                             tr_y_xzzz_zzzz, \
                             tr_y_zzz_xxx,   \
                             tr_y_zzz_xxxx,  \
                             tr_y_zzz_xxxy,  \
                             tr_y_zzz_xxxz,  \
                             tr_y_zzz_xxy,   \
                             tr_y_zzz_xxyy,  \
                             tr_y_zzz_xxyz,  \
                             tr_y_zzz_xxz,   \
                             tr_y_zzz_xxzz,  \
                             tr_y_zzz_xyy,   \
                             tr_y_zzz_xyyy,  \
                             tr_y_zzz_xyyz,  \
                             tr_y_zzz_xyz,   \
                             tr_y_zzz_xyzz,  \
                             tr_y_zzz_xzz,   \
                             tr_y_zzz_xzzz,  \
                             tr_y_zzz_yyy,   \
                             tr_y_zzz_yyyy,  \
                             tr_y_zzz_yyyz,  \
                             tr_y_zzz_yyz,   \
                             tr_y_zzz_yyzz,  \
                             tr_y_zzz_yzz,   \
                             tr_y_zzz_yzzz,  \
                             tr_y_zzz_zzz,   \
                             tr_y_zzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzz_xxxx[i] = 4.0 * tr_y_zzz_xxx[i] * fe_0 + tr_y_zzz_xxxx[i] * pa_x[i];

        tr_y_xzzz_xxxy[i] = 3.0 * tr_y_zzz_xxy[i] * fe_0 + tr_y_zzz_xxxy[i] * pa_x[i];

        tr_y_xzzz_xxxz[i] = 3.0 * tr_y_zzz_xxz[i] * fe_0 + tr_y_zzz_xxxz[i] * pa_x[i];

        tr_y_xzzz_xxyy[i] = 2.0 * tr_y_zzz_xyy[i] * fe_0 + tr_y_zzz_xxyy[i] * pa_x[i];

        tr_y_xzzz_xxyz[i] = 2.0 * tr_y_zzz_xyz[i] * fe_0 + tr_y_zzz_xxyz[i] * pa_x[i];

        tr_y_xzzz_xxzz[i] = 2.0 * tr_y_zzz_xzz[i] * fe_0 + tr_y_zzz_xxzz[i] * pa_x[i];

        tr_y_xzzz_xyyy[i] = tr_y_zzz_yyy[i] * fe_0 + tr_y_zzz_xyyy[i] * pa_x[i];

        tr_y_xzzz_xyyz[i] = tr_y_zzz_yyz[i] * fe_0 + tr_y_zzz_xyyz[i] * pa_x[i];

        tr_y_xzzz_xyzz[i] = tr_y_zzz_yzz[i] * fe_0 + tr_y_zzz_xyzz[i] * pa_x[i];

        tr_y_xzzz_xzzz[i] = tr_y_zzz_zzz[i] * fe_0 + tr_y_zzz_xzzz[i] * pa_x[i];

        tr_y_xzzz_yyyy[i] = tr_y_zzz_yyyy[i] * pa_x[i];

        tr_y_xzzz_yyyz[i] = tr_y_zzz_yyyz[i] * pa_x[i];

        tr_y_xzzz_yyzz[i] = tr_y_zzz_yyzz[i] * pa_x[i];

        tr_y_xzzz_yzzz[i] = tr_y_zzz_yzzz[i] * pa_x[i];

        tr_y_xzzz_zzzz[i] = tr_y_zzz_zzzz[i] * pa_x[i];
    }

    // Set up 375-390 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_y,               \
                             tr_y_yy_xxxx,   \
                             tr_y_yy_xxxy,   \
                             tr_y_yy_xxxz,   \
                             tr_y_yy_xxyy,   \
                             tr_y_yy_xxyz,   \
                             tr_y_yy_xxzz,   \
                             tr_y_yy_xyyy,   \
                             tr_y_yy_xyyz,   \
                             tr_y_yy_xyzz,   \
                             tr_y_yy_xzzz,   \
                             tr_y_yy_yyyy,   \
                             tr_y_yy_yyyz,   \
                             tr_y_yy_yyzz,   \
                             tr_y_yy_yzzz,   \
                             tr_y_yy_zzzz,   \
                             tr_y_yyy_xxx,   \
                             tr_y_yyy_xxxx,  \
                             tr_y_yyy_xxxy,  \
                             tr_y_yyy_xxxz,  \
                             tr_y_yyy_xxy,   \
                             tr_y_yyy_xxyy,  \
                             tr_y_yyy_xxyz,  \
                             tr_y_yyy_xxz,   \
                             tr_y_yyy_xxzz,  \
                             tr_y_yyy_xyy,   \
                             tr_y_yyy_xyyy,  \
                             tr_y_yyy_xyyz,  \
                             tr_y_yyy_xyz,   \
                             tr_y_yyy_xyzz,  \
                             tr_y_yyy_xzz,   \
                             tr_y_yyy_xzzz,  \
                             tr_y_yyy_yyy,   \
                             tr_y_yyy_yyyy,  \
                             tr_y_yyy_yyyz,  \
                             tr_y_yyy_yyz,   \
                             tr_y_yyy_yyzz,  \
                             tr_y_yyy_yzz,   \
                             tr_y_yyy_yzzz,  \
                             tr_y_yyy_zzz,   \
                             tr_y_yyy_zzzz,  \
                             tr_y_yyyy_xxxx, \
                             tr_y_yyyy_xxxy, \
                             tr_y_yyyy_xxxz, \
                             tr_y_yyyy_xxyy, \
                             tr_y_yyyy_xxyz, \
                             tr_y_yyyy_xxzz, \
                             tr_y_yyyy_xyyy, \
                             tr_y_yyyy_xyyz, \
                             tr_y_yyyy_xyzz, \
                             tr_y_yyyy_xzzz, \
                             tr_y_yyyy_yyyy, \
                             tr_y_yyyy_yyyz, \
                             tr_y_yyyy_yyzz, \
                             tr_y_yyyy_yzzz, \
                             tr_y_yyyy_zzzz, \
                             ts_yyy_xxxx,    \
                             ts_yyy_xxxy,    \
                             ts_yyy_xxxz,    \
                             ts_yyy_xxyy,    \
                             ts_yyy_xxyz,    \
                             ts_yyy_xxzz,    \
                             ts_yyy_xyyy,    \
                             ts_yyy_xyyz,    \
                             ts_yyy_xyzz,    \
                             ts_yyy_xzzz,    \
                             ts_yyy_yyyy,    \
                             ts_yyy_yyyz,    \
                             ts_yyy_yyzz,    \
                             ts_yyy_yzzz,    \
                             ts_yyy_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyy_xxxx[i] = 3.0 * tr_y_yy_xxxx[i] * fe_0 + ts_yyy_xxxx[i] * fe_0 + tr_y_yyy_xxxx[i] * pa_y[i];

        tr_y_yyyy_xxxy[i] = 3.0 * tr_y_yy_xxxy[i] * fe_0 + tr_y_yyy_xxx[i] * fe_0 + ts_yyy_xxxy[i] * fe_0 + tr_y_yyy_xxxy[i] * pa_y[i];

        tr_y_yyyy_xxxz[i] = 3.0 * tr_y_yy_xxxz[i] * fe_0 + ts_yyy_xxxz[i] * fe_0 + tr_y_yyy_xxxz[i] * pa_y[i];

        tr_y_yyyy_xxyy[i] = 3.0 * tr_y_yy_xxyy[i] * fe_0 + 2.0 * tr_y_yyy_xxy[i] * fe_0 + ts_yyy_xxyy[i] * fe_0 + tr_y_yyy_xxyy[i] * pa_y[i];

        tr_y_yyyy_xxyz[i] = 3.0 * tr_y_yy_xxyz[i] * fe_0 + tr_y_yyy_xxz[i] * fe_0 + ts_yyy_xxyz[i] * fe_0 + tr_y_yyy_xxyz[i] * pa_y[i];

        tr_y_yyyy_xxzz[i] = 3.0 * tr_y_yy_xxzz[i] * fe_0 + ts_yyy_xxzz[i] * fe_0 + tr_y_yyy_xxzz[i] * pa_y[i];

        tr_y_yyyy_xyyy[i] = 3.0 * tr_y_yy_xyyy[i] * fe_0 + 3.0 * tr_y_yyy_xyy[i] * fe_0 + ts_yyy_xyyy[i] * fe_0 + tr_y_yyy_xyyy[i] * pa_y[i];

        tr_y_yyyy_xyyz[i] = 3.0 * tr_y_yy_xyyz[i] * fe_0 + 2.0 * tr_y_yyy_xyz[i] * fe_0 + ts_yyy_xyyz[i] * fe_0 + tr_y_yyy_xyyz[i] * pa_y[i];

        tr_y_yyyy_xyzz[i] = 3.0 * tr_y_yy_xyzz[i] * fe_0 + tr_y_yyy_xzz[i] * fe_0 + ts_yyy_xyzz[i] * fe_0 + tr_y_yyy_xyzz[i] * pa_y[i];

        tr_y_yyyy_xzzz[i] = 3.0 * tr_y_yy_xzzz[i] * fe_0 + ts_yyy_xzzz[i] * fe_0 + tr_y_yyy_xzzz[i] * pa_y[i];

        tr_y_yyyy_yyyy[i] = 3.0 * tr_y_yy_yyyy[i] * fe_0 + 4.0 * tr_y_yyy_yyy[i] * fe_0 + ts_yyy_yyyy[i] * fe_0 + tr_y_yyy_yyyy[i] * pa_y[i];

        tr_y_yyyy_yyyz[i] = 3.0 * tr_y_yy_yyyz[i] * fe_0 + 3.0 * tr_y_yyy_yyz[i] * fe_0 + ts_yyy_yyyz[i] * fe_0 + tr_y_yyy_yyyz[i] * pa_y[i];

        tr_y_yyyy_yyzz[i] = 3.0 * tr_y_yy_yyzz[i] * fe_0 + 2.0 * tr_y_yyy_yzz[i] * fe_0 + ts_yyy_yyzz[i] * fe_0 + tr_y_yyy_yyzz[i] * pa_y[i];

        tr_y_yyyy_yzzz[i] = 3.0 * tr_y_yy_yzzz[i] * fe_0 + tr_y_yyy_zzz[i] * fe_0 + ts_yyy_yzzz[i] * fe_0 + tr_y_yyy_yzzz[i] * pa_y[i];

        tr_y_yyyy_zzzz[i] = 3.0 * tr_y_yy_zzzz[i] * fe_0 + ts_yyy_zzzz[i] * fe_0 + tr_y_yyy_zzzz[i] * pa_y[i];
    }

    // Set up 390-405 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_z,               \
                             tr_y_yyy_xxx,   \
                             tr_y_yyy_xxxx,  \
                             tr_y_yyy_xxxy,  \
                             tr_y_yyy_xxxz,  \
                             tr_y_yyy_xxy,   \
                             tr_y_yyy_xxyy,  \
                             tr_y_yyy_xxyz,  \
                             tr_y_yyy_xxz,   \
                             tr_y_yyy_xxzz,  \
                             tr_y_yyy_xyy,   \
                             tr_y_yyy_xyyy,  \
                             tr_y_yyy_xyyz,  \
                             tr_y_yyy_xyz,   \
                             tr_y_yyy_xyzz,  \
                             tr_y_yyy_xzz,   \
                             tr_y_yyy_xzzz,  \
                             tr_y_yyy_yyy,   \
                             tr_y_yyy_yyyy,  \
                             tr_y_yyy_yyyz,  \
                             tr_y_yyy_yyz,   \
                             tr_y_yyy_yyzz,  \
                             tr_y_yyy_yzz,   \
                             tr_y_yyy_yzzz,  \
                             tr_y_yyy_zzz,   \
                             tr_y_yyy_zzzz,  \
                             tr_y_yyyz_xxxx, \
                             tr_y_yyyz_xxxy, \
                             tr_y_yyyz_xxxz, \
                             tr_y_yyyz_xxyy, \
                             tr_y_yyyz_xxyz, \
                             tr_y_yyyz_xxzz, \
                             tr_y_yyyz_xyyy, \
                             tr_y_yyyz_xyyz, \
                             tr_y_yyyz_xyzz, \
                             tr_y_yyyz_xzzz, \
                             tr_y_yyyz_yyyy, \
                             tr_y_yyyz_yyyz, \
                             tr_y_yyyz_yyzz, \
                             tr_y_yyyz_yzzz, \
                             tr_y_yyyz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyz_xxxx[i] = tr_y_yyy_xxxx[i] * pa_z[i];

        tr_y_yyyz_xxxy[i] = tr_y_yyy_xxxy[i] * pa_z[i];

        tr_y_yyyz_xxxz[i] = tr_y_yyy_xxx[i] * fe_0 + tr_y_yyy_xxxz[i] * pa_z[i];

        tr_y_yyyz_xxyy[i] = tr_y_yyy_xxyy[i] * pa_z[i];

        tr_y_yyyz_xxyz[i] = tr_y_yyy_xxy[i] * fe_0 + tr_y_yyy_xxyz[i] * pa_z[i];

        tr_y_yyyz_xxzz[i] = 2.0 * tr_y_yyy_xxz[i] * fe_0 + tr_y_yyy_xxzz[i] * pa_z[i];

        tr_y_yyyz_xyyy[i] = tr_y_yyy_xyyy[i] * pa_z[i];

        tr_y_yyyz_xyyz[i] = tr_y_yyy_xyy[i] * fe_0 + tr_y_yyy_xyyz[i] * pa_z[i];

        tr_y_yyyz_xyzz[i] = 2.0 * tr_y_yyy_xyz[i] * fe_0 + tr_y_yyy_xyzz[i] * pa_z[i];

        tr_y_yyyz_xzzz[i] = 3.0 * tr_y_yyy_xzz[i] * fe_0 + tr_y_yyy_xzzz[i] * pa_z[i];

        tr_y_yyyz_yyyy[i] = tr_y_yyy_yyyy[i] * pa_z[i];

        tr_y_yyyz_yyyz[i] = tr_y_yyy_yyy[i] * fe_0 + tr_y_yyy_yyyz[i] * pa_z[i];

        tr_y_yyyz_yyzz[i] = 2.0 * tr_y_yyy_yyz[i] * fe_0 + tr_y_yyy_yyzz[i] * pa_z[i];

        tr_y_yyyz_yzzz[i] = 3.0 * tr_y_yyy_yzz[i] * fe_0 + tr_y_yyy_yzzz[i] * pa_z[i];

        tr_y_yyyz_zzzz[i] = 4.0 * tr_y_yyy_zzz[i] * fe_0 + tr_y_yyy_zzzz[i] * pa_z[i];
    }

    // Set up 405-420 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_y_yy_xxxx,   \
                             tr_y_yy_xxxy,   \
                             tr_y_yy_xxyy,   \
                             tr_y_yy_xxyz,   \
                             tr_y_yy_xyyy,   \
                             tr_y_yy_xyyz,   \
                             tr_y_yy_xyzz,   \
                             tr_y_yy_yyyy,   \
                             tr_y_yy_yyyz,   \
                             tr_y_yy_yyzz,   \
                             tr_y_yy_yzzz,   \
                             tr_y_yyz_xxxx,  \
                             tr_y_yyz_xxxy,  \
                             tr_y_yyz_xxy,   \
                             tr_y_yyz_xxyy,  \
                             tr_y_yyz_xxyz,  \
                             tr_y_yyz_xyy,   \
                             tr_y_yyz_xyyy,  \
                             tr_y_yyz_xyyz,  \
                             tr_y_yyz_xyz,   \
                             tr_y_yyz_xyzz,  \
                             tr_y_yyz_yyy,   \
                             tr_y_yyz_yyyy,  \
                             tr_y_yyz_yyyz,  \
                             tr_y_yyz_yyz,   \
                             tr_y_yyz_yyzz,  \
                             tr_y_yyz_yzz,   \
                             tr_y_yyz_yzzz,  \
                             tr_y_yyzz_xxxx, \
                             tr_y_yyzz_xxxy, \
                             tr_y_yyzz_xxxz, \
                             tr_y_yyzz_xxyy, \
                             tr_y_yyzz_xxyz, \
                             tr_y_yyzz_xxzz, \
                             tr_y_yyzz_xyyy, \
                             tr_y_yyzz_xyyz, \
                             tr_y_yyzz_xyzz, \
                             tr_y_yyzz_xzzz, \
                             tr_y_yyzz_yyyy, \
                             tr_y_yyzz_yyyz, \
                             tr_y_yyzz_yyzz, \
                             tr_y_yyzz_yzzz, \
                             tr_y_yyzz_zzzz, \
                             tr_y_yzz_xxxz,  \
                             tr_y_yzz_xxzz,  \
                             tr_y_yzz_xzzz,  \
                             tr_y_yzz_zzzz,  \
                             tr_y_zz_xxxz,   \
                             tr_y_zz_xxzz,   \
                             tr_y_zz_xzzz,   \
                             tr_y_zz_zzzz,   \
                             ts_yzz_xxxz,    \
                             ts_yzz_xxzz,    \
                             ts_yzz_xzzz,    \
                             ts_yzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzz_xxxx[i] = tr_y_yy_xxxx[i] * fe_0 + tr_y_yyz_xxxx[i] * pa_z[i];

        tr_y_yyzz_xxxy[i] = tr_y_yy_xxxy[i] * fe_0 + tr_y_yyz_xxxy[i] * pa_z[i];

        tr_y_yyzz_xxxz[i] = tr_y_zz_xxxz[i] * fe_0 + ts_yzz_xxxz[i] * fe_0 + tr_y_yzz_xxxz[i] * pa_y[i];

        tr_y_yyzz_xxyy[i] = tr_y_yy_xxyy[i] * fe_0 + tr_y_yyz_xxyy[i] * pa_z[i];

        tr_y_yyzz_xxyz[i] = tr_y_yy_xxyz[i] * fe_0 + tr_y_yyz_xxy[i] * fe_0 + tr_y_yyz_xxyz[i] * pa_z[i];

        tr_y_yyzz_xxzz[i] = tr_y_zz_xxzz[i] * fe_0 + ts_yzz_xxzz[i] * fe_0 + tr_y_yzz_xxzz[i] * pa_y[i];

        tr_y_yyzz_xyyy[i] = tr_y_yy_xyyy[i] * fe_0 + tr_y_yyz_xyyy[i] * pa_z[i];

        tr_y_yyzz_xyyz[i] = tr_y_yy_xyyz[i] * fe_0 + tr_y_yyz_xyy[i] * fe_0 + tr_y_yyz_xyyz[i] * pa_z[i];

        tr_y_yyzz_xyzz[i] = tr_y_yy_xyzz[i] * fe_0 + 2.0 * tr_y_yyz_xyz[i] * fe_0 + tr_y_yyz_xyzz[i] * pa_z[i];

        tr_y_yyzz_xzzz[i] = tr_y_zz_xzzz[i] * fe_0 + ts_yzz_xzzz[i] * fe_0 + tr_y_yzz_xzzz[i] * pa_y[i];

        tr_y_yyzz_yyyy[i] = tr_y_yy_yyyy[i] * fe_0 + tr_y_yyz_yyyy[i] * pa_z[i];

        tr_y_yyzz_yyyz[i] = tr_y_yy_yyyz[i] * fe_0 + tr_y_yyz_yyy[i] * fe_0 + tr_y_yyz_yyyz[i] * pa_z[i];

        tr_y_yyzz_yyzz[i] = tr_y_yy_yyzz[i] * fe_0 + 2.0 * tr_y_yyz_yyz[i] * fe_0 + tr_y_yyz_yyzz[i] * pa_z[i];

        tr_y_yyzz_yzzz[i] = tr_y_yy_yzzz[i] * fe_0 + 3.0 * tr_y_yyz_yzz[i] * fe_0 + tr_y_yyz_yzzz[i] * pa_z[i];

        tr_y_yyzz_zzzz[i] = tr_y_zz_zzzz[i] * fe_0 + ts_yzz_zzzz[i] * fe_0 + tr_y_yzz_zzzz[i] * pa_y[i];
    }

    // Set up 420-435 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_y_yz_xxxy,   \
                             tr_y_yz_xxyy,   \
                             tr_y_yz_xyyy,   \
                             tr_y_yz_yyyy,   \
                             tr_y_yzz_xxxy,  \
                             tr_y_yzz_xxyy,  \
                             tr_y_yzz_xyyy,  \
                             tr_y_yzz_yyyy,  \
                             tr_y_yzzz_xxxx, \
                             tr_y_yzzz_xxxy, \
                             tr_y_yzzz_xxxz, \
                             tr_y_yzzz_xxyy, \
                             tr_y_yzzz_xxyz, \
                             tr_y_yzzz_xxzz, \
                             tr_y_yzzz_xyyy, \
                             tr_y_yzzz_xyyz, \
                             tr_y_yzzz_xyzz, \
                             tr_y_yzzz_xzzz, \
                             tr_y_yzzz_yyyy, \
                             tr_y_yzzz_yyyz, \
                             tr_y_yzzz_yyzz, \
                             tr_y_yzzz_yzzz, \
                             tr_y_yzzz_zzzz, \
                             tr_y_zzz_xxxx,  \
                             tr_y_zzz_xxxz,  \
                             tr_y_zzz_xxyz,  \
                             tr_y_zzz_xxz,   \
                             tr_y_zzz_xxzz,  \
                             tr_y_zzz_xyyz,  \
                             tr_y_zzz_xyz,   \
                             tr_y_zzz_xyzz,  \
                             tr_y_zzz_xzz,   \
                             tr_y_zzz_xzzz,  \
                             tr_y_zzz_yyyz,  \
                             tr_y_zzz_yyz,   \
                             tr_y_zzz_yyzz,  \
                             tr_y_zzz_yzz,   \
                             tr_y_zzz_yzzz,  \
                             tr_y_zzz_zzz,   \
                             tr_y_zzz_zzzz,  \
                             ts_zzz_xxxx,    \
                             ts_zzz_xxxz,    \
                             ts_zzz_xxyz,    \
                             ts_zzz_xxzz,    \
                             ts_zzz_xyyz,    \
                             ts_zzz_xyzz,    \
                             ts_zzz_xzzz,    \
                             ts_zzz_yyyz,    \
                             ts_zzz_yyzz,    \
                             ts_zzz_yzzz,    \
                             ts_zzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzz_xxxx[i] = ts_zzz_xxxx[i] * fe_0 + tr_y_zzz_xxxx[i] * pa_y[i];

        tr_y_yzzz_xxxy[i] = 2.0 * tr_y_yz_xxxy[i] * fe_0 + tr_y_yzz_xxxy[i] * pa_z[i];

        tr_y_yzzz_xxxz[i] = ts_zzz_xxxz[i] * fe_0 + tr_y_zzz_xxxz[i] * pa_y[i];

        tr_y_yzzz_xxyy[i] = 2.0 * tr_y_yz_xxyy[i] * fe_0 + tr_y_yzz_xxyy[i] * pa_z[i];

        tr_y_yzzz_xxyz[i] = tr_y_zzz_xxz[i] * fe_0 + ts_zzz_xxyz[i] * fe_0 + tr_y_zzz_xxyz[i] * pa_y[i];

        tr_y_yzzz_xxzz[i] = ts_zzz_xxzz[i] * fe_0 + tr_y_zzz_xxzz[i] * pa_y[i];

        tr_y_yzzz_xyyy[i] = 2.0 * tr_y_yz_xyyy[i] * fe_0 + tr_y_yzz_xyyy[i] * pa_z[i];

        tr_y_yzzz_xyyz[i] = 2.0 * tr_y_zzz_xyz[i] * fe_0 + ts_zzz_xyyz[i] * fe_0 + tr_y_zzz_xyyz[i] * pa_y[i];

        tr_y_yzzz_xyzz[i] = tr_y_zzz_xzz[i] * fe_0 + ts_zzz_xyzz[i] * fe_0 + tr_y_zzz_xyzz[i] * pa_y[i];

        tr_y_yzzz_xzzz[i] = ts_zzz_xzzz[i] * fe_0 + tr_y_zzz_xzzz[i] * pa_y[i];

        tr_y_yzzz_yyyy[i] = 2.0 * tr_y_yz_yyyy[i] * fe_0 + tr_y_yzz_yyyy[i] * pa_z[i];

        tr_y_yzzz_yyyz[i] = 3.0 * tr_y_zzz_yyz[i] * fe_0 + ts_zzz_yyyz[i] * fe_0 + tr_y_zzz_yyyz[i] * pa_y[i];

        tr_y_yzzz_yyzz[i] = 2.0 * tr_y_zzz_yzz[i] * fe_0 + ts_zzz_yyzz[i] * fe_0 + tr_y_zzz_yyzz[i] * pa_y[i];

        tr_y_yzzz_yzzz[i] = tr_y_zzz_zzz[i] * fe_0 + ts_zzz_yzzz[i] * fe_0 + tr_y_zzz_yzzz[i] * pa_y[i];

        tr_y_yzzz_zzzz[i] = ts_zzz_zzzz[i] * fe_0 + tr_y_zzz_zzzz[i] * pa_y[i];
    }

    // Set up 435-450 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_z,               \
                             tr_y_zz_xxxx,   \
                             tr_y_zz_xxxy,   \
                             tr_y_zz_xxxz,   \
                             tr_y_zz_xxyy,   \
                             tr_y_zz_xxyz,   \
                             tr_y_zz_xxzz,   \
                             tr_y_zz_xyyy,   \
                             tr_y_zz_xyyz,   \
                             tr_y_zz_xyzz,   \
                             tr_y_zz_xzzz,   \
                             tr_y_zz_yyyy,   \
                             tr_y_zz_yyyz,   \
                             tr_y_zz_yyzz,   \
                             tr_y_zz_yzzz,   \
                             tr_y_zz_zzzz,   \
                             tr_y_zzz_xxx,   \
                             tr_y_zzz_xxxx,  \
                             tr_y_zzz_xxxy,  \
                             tr_y_zzz_xxxz,  \
                             tr_y_zzz_xxy,   \
                             tr_y_zzz_xxyy,  \
                             tr_y_zzz_xxyz,  \
                             tr_y_zzz_xxz,   \
                             tr_y_zzz_xxzz,  \
                             tr_y_zzz_xyy,   \
                             tr_y_zzz_xyyy,  \
                             tr_y_zzz_xyyz,  \
                             tr_y_zzz_xyz,   \
                             tr_y_zzz_xyzz,  \
                             tr_y_zzz_xzz,   \
                             tr_y_zzz_xzzz,  \
                             tr_y_zzz_yyy,   \
                             tr_y_zzz_yyyy,  \
                             tr_y_zzz_yyyz,  \
                             tr_y_zzz_yyz,   \
                             tr_y_zzz_yyzz,  \
                             tr_y_zzz_yzz,   \
                             tr_y_zzz_yzzz,  \
                             tr_y_zzz_zzz,   \
                             tr_y_zzz_zzzz,  \
                             tr_y_zzzz_xxxx, \
                             tr_y_zzzz_xxxy, \
                             tr_y_zzzz_xxxz, \
                             tr_y_zzzz_xxyy, \
                             tr_y_zzzz_xxyz, \
                             tr_y_zzzz_xxzz, \
                             tr_y_zzzz_xyyy, \
                             tr_y_zzzz_xyyz, \
                             tr_y_zzzz_xyzz, \
                             tr_y_zzzz_xzzz, \
                             tr_y_zzzz_yyyy, \
                             tr_y_zzzz_yyyz, \
                             tr_y_zzzz_yyzz, \
                             tr_y_zzzz_yzzz, \
                             tr_y_zzzz_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzz_xxxx[i] = 3.0 * tr_y_zz_xxxx[i] * fe_0 + tr_y_zzz_xxxx[i] * pa_z[i];

        tr_y_zzzz_xxxy[i] = 3.0 * tr_y_zz_xxxy[i] * fe_0 + tr_y_zzz_xxxy[i] * pa_z[i];

        tr_y_zzzz_xxxz[i] = 3.0 * tr_y_zz_xxxz[i] * fe_0 + tr_y_zzz_xxx[i] * fe_0 + tr_y_zzz_xxxz[i] * pa_z[i];

        tr_y_zzzz_xxyy[i] = 3.0 * tr_y_zz_xxyy[i] * fe_0 + tr_y_zzz_xxyy[i] * pa_z[i];

        tr_y_zzzz_xxyz[i] = 3.0 * tr_y_zz_xxyz[i] * fe_0 + tr_y_zzz_xxy[i] * fe_0 + tr_y_zzz_xxyz[i] * pa_z[i];

        tr_y_zzzz_xxzz[i] = 3.0 * tr_y_zz_xxzz[i] * fe_0 + 2.0 * tr_y_zzz_xxz[i] * fe_0 + tr_y_zzz_xxzz[i] * pa_z[i];

        tr_y_zzzz_xyyy[i] = 3.0 * tr_y_zz_xyyy[i] * fe_0 + tr_y_zzz_xyyy[i] * pa_z[i];

        tr_y_zzzz_xyyz[i] = 3.0 * tr_y_zz_xyyz[i] * fe_0 + tr_y_zzz_xyy[i] * fe_0 + tr_y_zzz_xyyz[i] * pa_z[i];

        tr_y_zzzz_xyzz[i] = 3.0 * tr_y_zz_xyzz[i] * fe_0 + 2.0 * tr_y_zzz_xyz[i] * fe_0 + tr_y_zzz_xyzz[i] * pa_z[i];

        tr_y_zzzz_xzzz[i] = 3.0 * tr_y_zz_xzzz[i] * fe_0 + 3.0 * tr_y_zzz_xzz[i] * fe_0 + tr_y_zzz_xzzz[i] * pa_z[i];

        tr_y_zzzz_yyyy[i] = 3.0 * tr_y_zz_yyyy[i] * fe_0 + tr_y_zzz_yyyy[i] * pa_z[i];

        tr_y_zzzz_yyyz[i] = 3.0 * tr_y_zz_yyyz[i] * fe_0 + tr_y_zzz_yyy[i] * fe_0 + tr_y_zzz_yyyz[i] * pa_z[i];

        tr_y_zzzz_yyzz[i] = 3.0 * tr_y_zz_yyzz[i] * fe_0 + 2.0 * tr_y_zzz_yyz[i] * fe_0 + tr_y_zzz_yyzz[i] * pa_z[i];

        tr_y_zzzz_yzzz[i] = 3.0 * tr_y_zz_yzzz[i] * fe_0 + 3.0 * tr_y_zzz_yzz[i] * fe_0 + tr_y_zzz_yzzz[i] * pa_z[i];

        tr_y_zzzz_zzzz[i] = 3.0 * tr_y_zz_zzzz[i] * fe_0 + 4.0 * tr_y_zzz_zzz[i] * fe_0 + tr_y_zzz_zzzz[i] * pa_z[i];
    }

    // Set up 450-465 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,               \
                             tr_z_xx_xxxx,   \
                             tr_z_xx_xxxy,   \
                             tr_z_xx_xxxz,   \
                             tr_z_xx_xxyy,   \
                             tr_z_xx_xxyz,   \
                             tr_z_xx_xxzz,   \
                             tr_z_xx_xyyy,   \
                             tr_z_xx_xyyz,   \
                             tr_z_xx_xyzz,   \
                             tr_z_xx_xzzz,   \
                             tr_z_xx_yyyy,   \
                             tr_z_xx_yyyz,   \
                             tr_z_xx_yyzz,   \
                             tr_z_xx_yzzz,   \
                             tr_z_xx_zzzz,   \
                             tr_z_xxx_xxx,   \
                             tr_z_xxx_xxxx,  \
                             tr_z_xxx_xxxy,  \
                             tr_z_xxx_xxxz,  \
                             tr_z_xxx_xxy,   \
                             tr_z_xxx_xxyy,  \
                             tr_z_xxx_xxyz,  \
                             tr_z_xxx_xxz,   \
                             tr_z_xxx_xxzz,  \
                             tr_z_xxx_xyy,   \
                             tr_z_xxx_xyyy,  \
                             tr_z_xxx_xyyz,  \
                             tr_z_xxx_xyz,   \
                             tr_z_xxx_xyzz,  \
                             tr_z_xxx_xzz,   \
                             tr_z_xxx_xzzz,  \
                             tr_z_xxx_yyy,   \
                             tr_z_xxx_yyyy,  \
                             tr_z_xxx_yyyz,  \
                             tr_z_xxx_yyz,   \
                             tr_z_xxx_yyzz,  \
                             tr_z_xxx_yzz,   \
                             tr_z_xxx_yzzz,  \
                             tr_z_xxx_zzz,   \
                             tr_z_xxx_zzzz,  \
                             tr_z_xxxx_xxxx, \
                             tr_z_xxxx_xxxy, \
                             tr_z_xxxx_xxxz, \
                             tr_z_xxxx_xxyy, \
                             tr_z_xxxx_xxyz, \
                             tr_z_xxxx_xxzz, \
                             tr_z_xxxx_xyyy, \
                             tr_z_xxxx_xyyz, \
                             tr_z_xxxx_xyzz, \
                             tr_z_xxxx_xzzz, \
                             tr_z_xxxx_yyyy, \
                             tr_z_xxxx_yyyz, \
                             tr_z_xxxx_yyzz, \
                             tr_z_xxxx_yzzz, \
                             tr_z_xxxx_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxx_xxxx[i] = 3.0 * tr_z_xx_xxxx[i] * fe_0 + 4.0 * tr_z_xxx_xxx[i] * fe_0 + tr_z_xxx_xxxx[i] * pa_x[i];

        tr_z_xxxx_xxxy[i] = 3.0 * tr_z_xx_xxxy[i] * fe_0 + 3.0 * tr_z_xxx_xxy[i] * fe_0 + tr_z_xxx_xxxy[i] * pa_x[i];

        tr_z_xxxx_xxxz[i] = 3.0 * tr_z_xx_xxxz[i] * fe_0 + 3.0 * tr_z_xxx_xxz[i] * fe_0 + tr_z_xxx_xxxz[i] * pa_x[i];

        tr_z_xxxx_xxyy[i] = 3.0 * tr_z_xx_xxyy[i] * fe_0 + 2.0 * tr_z_xxx_xyy[i] * fe_0 + tr_z_xxx_xxyy[i] * pa_x[i];

        tr_z_xxxx_xxyz[i] = 3.0 * tr_z_xx_xxyz[i] * fe_0 + 2.0 * tr_z_xxx_xyz[i] * fe_0 + tr_z_xxx_xxyz[i] * pa_x[i];

        tr_z_xxxx_xxzz[i] = 3.0 * tr_z_xx_xxzz[i] * fe_0 + 2.0 * tr_z_xxx_xzz[i] * fe_0 + tr_z_xxx_xxzz[i] * pa_x[i];

        tr_z_xxxx_xyyy[i] = 3.0 * tr_z_xx_xyyy[i] * fe_0 + tr_z_xxx_yyy[i] * fe_0 + tr_z_xxx_xyyy[i] * pa_x[i];

        tr_z_xxxx_xyyz[i] = 3.0 * tr_z_xx_xyyz[i] * fe_0 + tr_z_xxx_yyz[i] * fe_0 + tr_z_xxx_xyyz[i] * pa_x[i];

        tr_z_xxxx_xyzz[i] = 3.0 * tr_z_xx_xyzz[i] * fe_0 + tr_z_xxx_yzz[i] * fe_0 + tr_z_xxx_xyzz[i] * pa_x[i];

        tr_z_xxxx_xzzz[i] = 3.0 * tr_z_xx_xzzz[i] * fe_0 + tr_z_xxx_zzz[i] * fe_0 + tr_z_xxx_xzzz[i] * pa_x[i];

        tr_z_xxxx_yyyy[i] = 3.0 * tr_z_xx_yyyy[i] * fe_0 + tr_z_xxx_yyyy[i] * pa_x[i];

        tr_z_xxxx_yyyz[i] = 3.0 * tr_z_xx_yyyz[i] * fe_0 + tr_z_xxx_yyyz[i] * pa_x[i];

        tr_z_xxxx_yyzz[i] = 3.0 * tr_z_xx_yyzz[i] * fe_0 + tr_z_xxx_yyzz[i] * pa_x[i];

        tr_z_xxxx_yzzz[i] = 3.0 * tr_z_xx_yzzz[i] * fe_0 + tr_z_xxx_yzzz[i] * pa_x[i];

        tr_z_xxxx_zzzz[i] = 3.0 * tr_z_xx_zzzz[i] * fe_0 + tr_z_xxx_zzzz[i] * pa_x[i];
    }

    // Set up 465-480 components of targeted buffer : GG

    auto tr_z_xxxy_xxxx = pbuffer.data(idx_dip_gg + 465);

    auto tr_z_xxxy_xxxy = pbuffer.data(idx_dip_gg + 466);

    auto tr_z_xxxy_xxxz = pbuffer.data(idx_dip_gg + 467);

    auto tr_z_xxxy_xxyy = pbuffer.data(idx_dip_gg + 468);

    auto tr_z_xxxy_xxyz = pbuffer.data(idx_dip_gg + 469);

    auto tr_z_xxxy_xxzz = pbuffer.data(idx_dip_gg + 470);

    auto tr_z_xxxy_xyyy = pbuffer.data(idx_dip_gg + 471);

    auto tr_z_xxxy_xyyz = pbuffer.data(idx_dip_gg + 472);

    auto tr_z_xxxy_xyzz = pbuffer.data(idx_dip_gg + 473);

    auto tr_z_xxxy_xzzz = pbuffer.data(idx_dip_gg + 474);

    auto tr_z_xxxy_yyyy = pbuffer.data(idx_dip_gg + 475);

    auto tr_z_xxxy_yyyz = pbuffer.data(idx_dip_gg + 476);

    auto tr_z_xxxy_yyzz = pbuffer.data(idx_dip_gg + 477);

    auto tr_z_xxxy_yzzz = pbuffer.data(idx_dip_gg + 478);

    auto tr_z_xxxy_zzzz = pbuffer.data(idx_dip_gg + 479);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xxx_xxx,   \
                             tr_z_xxx_xxxx,  \
                             tr_z_xxx_xxxy,  \
                             tr_z_xxx_xxxz,  \
                             tr_z_xxx_xxy,   \
                             tr_z_xxx_xxyy,  \
                             tr_z_xxx_xxyz,  \
                             tr_z_xxx_xxz,   \
                             tr_z_xxx_xxzz,  \
                             tr_z_xxx_xyy,   \
                             tr_z_xxx_xyyy,  \
                             tr_z_xxx_xyyz,  \
                             tr_z_xxx_xyz,   \
                             tr_z_xxx_xyzz,  \
                             tr_z_xxx_xzz,   \
                             tr_z_xxx_xzzz,  \
                             tr_z_xxx_zzzz,  \
                             tr_z_xxxy_xxxx, \
                             tr_z_xxxy_xxxy, \
                             tr_z_xxxy_xxxz, \
                             tr_z_xxxy_xxyy, \
                             tr_z_xxxy_xxyz, \
                             tr_z_xxxy_xxzz, \
                             tr_z_xxxy_xyyy, \
                             tr_z_xxxy_xyyz, \
                             tr_z_xxxy_xyzz, \
                             tr_z_xxxy_xzzz, \
                             tr_z_xxxy_yyyy, \
                             tr_z_xxxy_yyyz, \
                             tr_z_xxxy_yyzz, \
                             tr_z_xxxy_yzzz, \
                             tr_z_xxxy_zzzz, \
                             tr_z_xxy_yyyy,  \
                             tr_z_xxy_yyyz,  \
                             tr_z_xxy_yyzz,  \
                             tr_z_xxy_yzzz,  \
                             tr_z_xy_yyyy,   \
                             tr_z_xy_yyyz,   \
                             tr_z_xy_yyzz,   \
                             tr_z_xy_yzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxy_xxxx[i] = tr_z_xxx_xxxx[i] * pa_y[i];

        tr_z_xxxy_xxxy[i] = tr_z_xxx_xxx[i] * fe_0 + tr_z_xxx_xxxy[i] * pa_y[i];

        tr_z_xxxy_xxxz[i] = tr_z_xxx_xxxz[i] * pa_y[i];

        tr_z_xxxy_xxyy[i] = 2.0 * tr_z_xxx_xxy[i] * fe_0 + tr_z_xxx_xxyy[i] * pa_y[i];

        tr_z_xxxy_xxyz[i] = tr_z_xxx_xxz[i] * fe_0 + tr_z_xxx_xxyz[i] * pa_y[i];

        tr_z_xxxy_xxzz[i] = tr_z_xxx_xxzz[i] * pa_y[i];

        tr_z_xxxy_xyyy[i] = 3.0 * tr_z_xxx_xyy[i] * fe_0 + tr_z_xxx_xyyy[i] * pa_y[i];

        tr_z_xxxy_xyyz[i] = 2.0 * tr_z_xxx_xyz[i] * fe_0 + tr_z_xxx_xyyz[i] * pa_y[i];

        tr_z_xxxy_xyzz[i] = tr_z_xxx_xzz[i] * fe_0 + tr_z_xxx_xyzz[i] * pa_y[i];

        tr_z_xxxy_xzzz[i] = tr_z_xxx_xzzz[i] * pa_y[i];

        tr_z_xxxy_yyyy[i] = 2.0 * tr_z_xy_yyyy[i] * fe_0 + tr_z_xxy_yyyy[i] * pa_x[i];

        tr_z_xxxy_yyyz[i] = 2.0 * tr_z_xy_yyyz[i] * fe_0 + tr_z_xxy_yyyz[i] * pa_x[i];

        tr_z_xxxy_yyzz[i] = 2.0 * tr_z_xy_yyzz[i] * fe_0 + tr_z_xxy_yyzz[i] * pa_x[i];

        tr_z_xxxy_yzzz[i] = 2.0 * tr_z_xy_yzzz[i] * fe_0 + tr_z_xxy_yzzz[i] * pa_x[i];

        tr_z_xxxy_zzzz[i] = tr_z_xxx_zzzz[i] * pa_y[i];
    }

    // Set up 480-495 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_z_xxx_xxxx,  \
                             tr_z_xxx_xxxy,  \
                             tr_z_xxx_xxyy,  \
                             tr_z_xxx_xyyy,  \
                             tr_z_xxxz_xxxx, \
                             tr_z_xxxz_xxxy, \
                             tr_z_xxxz_xxxz, \
                             tr_z_xxxz_xxyy, \
                             tr_z_xxxz_xxyz, \
                             tr_z_xxxz_xxzz, \
                             tr_z_xxxz_xyyy, \
                             tr_z_xxxz_xyyz, \
                             tr_z_xxxz_xyzz, \
                             tr_z_xxxz_xzzz, \
                             tr_z_xxxz_yyyy, \
                             tr_z_xxxz_yyyz, \
                             tr_z_xxxz_yyzz, \
                             tr_z_xxxz_yzzz, \
                             tr_z_xxxz_zzzz, \
                             tr_z_xxz_xxxz,  \
                             tr_z_xxz_xxyz,  \
                             tr_z_xxz_xxz,   \
                             tr_z_xxz_xxzz,  \
                             tr_z_xxz_xyyz,  \
                             tr_z_xxz_xyz,   \
                             tr_z_xxz_xyzz,  \
                             tr_z_xxz_xzz,   \
                             tr_z_xxz_xzzz,  \
                             tr_z_xxz_yyyy,  \
                             tr_z_xxz_yyyz,  \
                             tr_z_xxz_yyz,   \
                             tr_z_xxz_yyzz,  \
                             tr_z_xxz_yzz,   \
                             tr_z_xxz_yzzz,  \
                             tr_z_xxz_zzz,   \
                             tr_z_xxz_zzzz,  \
                             tr_z_xz_xxxz,   \
                             tr_z_xz_xxyz,   \
                             tr_z_xz_xxzz,   \
                             tr_z_xz_xyyz,   \
                             tr_z_xz_xyzz,   \
                             tr_z_xz_xzzz,   \
                             tr_z_xz_yyyy,   \
                             tr_z_xz_yyyz,   \
                             tr_z_xz_yyzz,   \
                             tr_z_xz_yzzz,   \
                             tr_z_xz_zzzz,   \
                             ts_xxx_xxxx,    \
                             ts_xxx_xxxy,    \
                             ts_xxx_xxyy,    \
                             ts_xxx_xyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxz_xxxx[i] = ts_xxx_xxxx[i] * fe_0 + tr_z_xxx_xxxx[i] * pa_z[i];

        tr_z_xxxz_xxxy[i] = ts_xxx_xxxy[i] * fe_0 + tr_z_xxx_xxxy[i] * pa_z[i];

        tr_z_xxxz_xxxz[i] = 2.0 * tr_z_xz_xxxz[i] * fe_0 + 3.0 * tr_z_xxz_xxz[i] * fe_0 + tr_z_xxz_xxxz[i] * pa_x[i];

        tr_z_xxxz_xxyy[i] = ts_xxx_xxyy[i] * fe_0 + tr_z_xxx_xxyy[i] * pa_z[i];

        tr_z_xxxz_xxyz[i] = 2.0 * tr_z_xz_xxyz[i] * fe_0 + 2.0 * tr_z_xxz_xyz[i] * fe_0 + tr_z_xxz_xxyz[i] * pa_x[i];

        tr_z_xxxz_xxzz[i] = 2.0 * tr_z_xz_xxzz[i] * fe_0 + 2.0 * tr_z_xxz_xzz[i] * fe_0 + tr_z_xxz_xxzz[i] * pa_x[i];

        tr_z_xxxz_xyyy[i] = ts_xxx_xyyy[i] * fe_0 + tr_z_xxx_xyyy[i] * pa_z[i];

        tr_z_xxxz_xyyz[i] = 2.0 * tr_z_xz_xyyz[i] * fe_0 + tr_z_xxz_yyz[i] * fe_0 + tr_z_xxz_xyyz[i] * pa_x[i];

        tr_z_xxxz_xyzz[i] = 2.0 * tr_z_xz_xyzz[i] * fe_0 + tr_z_xxz_yzz[i] * fe_0 + tr_z_xxz_xyzz[i] * pa_x[i];

        tr_z_xxxz_xzzz[i] = 2.0 * tr_z_xz_xzzz[i] * fe_0 + tr_z_xxz_zzz[i] * fe_0 + tr_z_xxz_xzzz[i] * pa_x[i];

        tr_z_xxxz_yyyy[i] = 2.0 * tr_z_xz_yyyy[i] * fe_0 + tr_z_xxz_yyyy[i] * pa_x[i];

        tr_z_xxxz_yyyz[i] = 2.0 * tr_z_xz_yyyz[i] * fe_0 + tr_z_xxz_yyyz[i] * pa_x[i];

        tr_z_xxxz_yyzz[i] = 2.0 * tr_z_xz_yyzz[i] * fe_0 + tr_z_xxz_yyzz[i] * pa_x[i];

        tr_z_xxxz_yzzz[i] = 2.0 * tr_z_xz_yzzz[i] * fe_0 + tr_z_xxz_yzzz[i] * pa_x[i];

        tr_z_xxxz_zzzz[i] = 2.0 * tr_z_xz_zzzz[i] * fe_0 + tr_z_xxz_zzzz[i] * pa_x[i];
    }

    // Set up 495-510 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xx_xxxx,   \
                             tr_z_xx_xxxz,   \
                             tr_z_xx_xxzz,   \
                             tr_z_xx_xzzz,   \
                             tr_z_xxy_xxxx,  \
                             tr_z_xxy_xxxz,  \
                             tr_z_xxy_xxzz,  \
                             tr_z_xxy_xzzz,  \
                             tr_z_xxyy_xxxx, \
                             tr_z_xxyy_xxxy, \
                             tr_z_xxyy_xxxz, \
                             tr_z_xxyy_xxyy, \
                             tr_z_xxyy_xxyz, \
                             tr_z_xxyy_xxzz, \
                             tr_z_xxyy_xyyy, \
                             tr_z_xxyy_xyyz, \
                             tr_z_xxyy_xyzz, \
                             tr_z_xxyy_xzzz, \
                             tr_z_xxyy_yyyy, \
                             tr_z_xxyy_yyyz, \
                             tr_z_xxyy_yyzz, \
                             tr_z_xxyy_yzzz, \
                             tr_z_xxyy_zzzz, \
                             tr_z_xyy_xxxy,  \
                             tr_z_xyy_xxy,   \
                             tr_z_xyy_xxyy,  \
                             tr_z_xyy_xxyz,  \
                             tr_z_xyy_xyy,   \
                             tr_z_xyy_xyyy,  \
                             tr_z_xyy_xyyz,  \
                             tr_z_xyy_xyz,   \
                             tr_z_xyy_xyzz,  \
                             tr_z_xyy_yyy,   \
                             tr_z_xyy_yyyy,  \
                             tr_z_xyy_yyyz,  \
                             tr_z_xyy_yyz,   \
                             tr_z_xyy_yyzz,  \
                             tr_z_xyy_yzz,   \
                             tr_z_xyy_yzzz,  \
                             tr_z_xyy_zzzz,  \
                             tr_z_yy_xxxy,   \
                             tr_z_yy_xxyy,   \
                             tr_z_yy_xxyz,   \
                             tr_z_yy_xyyy,   \
                             tr_z_yy_xyyz,   \
                             tr_z_yy_xyzz,   \
                             tr_z_yy_yyyy,   \
                             tr_z_yy_yyyz,   \
                             tr_z_yy_yyzz,   \
                             tr_z_yy_yzzz,   \
                             tr_z_yy_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyy_xxxx[i] = tr_z_xx_xxxx[i] * fe_0 + tr_z_xxy_xxxx[i] * pa_y[i];

        tr_z_xxyy_xxxy[i] = tr_z_yy_xxxy[i] * fe_0 + 3.0 * tr_z_xyy_xxy[i] * fe_0 + tr_z_xyy_xxxy[i] * pa_x[i];

        tr_z_xxyy_xxxz[i] = tr_z_xx_xxxz[i] * fe_0 + tr_z_xxy_xxxz[i] * pa_y[i];

        tr_z_xxyy_xxyy[i] = tr_z_yy_xxyy[i] * fe_0 + 2.0 * tr_z_xyy_xyy[i] * fe_0 + tr_z_xyy_xxyy[i] * pa_x[i];

        tr_z_xxyy_xxyz[i] = tr_z_yy_xxyz[i] * fe_0 + 2.0 * tr_z_xyy_xyz[i] * fe_0 + tr_z_xyy_xxyz[i] * pa_x[i];

        tr_z_xxyy_xxzz[i] = tr_z_xx_xxzz[i] * fe_0 + tr_z_xxy_xxzz[i] * pa_y[i];

        tr_z_xxyy_xyyy[i] = tr_z_yy_xyyy[i] * fe_0 + tr_z_xyy_yyy[i] * fe_0 + tr_z_xyy_xyyy[i] * pa_x[i];

        tr_z_xxyy_xyyz[i] = tr_z_yy_xyyz[i] * fe_0 + tr_z_xyy_yyz[i] * fe_0 + tr_z_xyy_xyyz[i] * pa_x[i];

        tr_z_xxyy_xyzz[i] = tr_z_yy_xyzz[i] * fe_0 + tr_z_xyy_yzz[i] * fe_0 + tr_z_xyy_xyzz[i] * pa_x[i];

        tr_z_xxyy_xzzz[i] = tr_z_xx_xzzz[i] * fe_0 + tr_z_xxy_xzzz[i] * pa_y[i];

        tr_z_xxyy_yyyy[i] = tr_z_yy_yyyy[i] * fe_0 + tr_z_xyy_yyyy[i] * pa_x[i];

        tr_z_xxyy_yyyz[i] = tr_z_yy_yyyz[i] * fe_0 + tr_z_xyy_yyyz[i] * pa_x[i];

        tr_z_xxyy_yyzz[i] = tr_z_yy_yyzz[i] * fe_0 + tr_z_xyy_yyzz[i] * pa_x[i];

        tr_z_xxyy_yzzz[i] = tr_z_yy_yzzz[i] * fe_0 + tr_z_xyy_yzzz[i] * pa_x[i];

        tr_z_xxyy_zzzz[i] = tr_z_yy_zzzz[i] * fe_0 + tr_z_xyy_zzzz[i] * pa_x[i];
    }

    // Set up 510-525 components of targeted buffer : GG

    auto tr_z_xxyz_xxxx = pbuffer.data(idx_dip_gg + 510);

    auto tr_z_xxyz_xxxy = pbuffer.data(idx_dip_gg + 511);

    auto tr_z_xxyz_xxxz = pbuffer.data(idx_dip_gg + 512);

    auto tr_z_xxyz_xxyy = pbuffer.data(idx_dip_gg + 513);

    auto tr_z_xxyz_xxyz = pbuffer.data(idx_dip_gg + 514);

    auto tr_z_xxyz_xxzz = pbuffer.data(idx_dip_gg + 515);

    auto tr_z_xxyz_xyyy = pbuffer.data(idx_dip_gg + 516);

    auto tr_z_xxyz_xyyz = pbuffer.data(idx_dip_gg + 517);

    auto tr_z_xxyz_xyzz = pbuffer.data(idx_dip_gg + 518);

    auto tr_z_xxyz_xzzz = pbuffer.data(idx_dip_gg + 519);

    auto tr_z_xxyz_yyyy = pbuffer.data(idx_dip_gg + 520);

    auto tr_z_xxyz_yyyz = pbuffer.data(idx_dip_gg + 521);

    auto tr_z_xxyz_yyzz = pbuffer.data(idx_dip_gg + 522);

    auto tr_z_xxyz_yzzz = pbuffer.data(idx_dip_gg + 523);

    auto tr_z_xxyz_zzzz = pbuffer.data(idx_dip_gg + 524);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xxyz_xxxx, \
                             tr_z_xxyz_xxxy, \
                             tr_z_xxyz_xxxz, \
                             tr_z_xxyz_xxyy, \
                             tr_z_xxyz_xxyz, \
                             tr_z_xxyz_xxzz, \
                             tr_z_xxyz_xyyy, \
                             tr_z_xxyz_xyyz, \
                             tr_z_xxyz_xyzz, \
                             tr_z_xxyz_xzzz, \
                             tr_z_xxyz_yyyy, \
                             tr_z_xxyz_yyyz, \
                             tr_z_xxyz_yyzz, \
                             tr_z_xxyz_yzzz, \
                             tr_z_xxyz_zzzz, \
                             tr_z_xxz_xxx,   \
                             tr_z_xxz_xxxx,  \
                             tr_z_xxz_xxxy,  \
                             tr_z_xxz_xxxz,  \
                             tr_z_xxz_xxy,   \
                             tr_z_xxz_xxyy,  \
                             tr_z_xxz_xxyz,  \
                             tr_z_xxz_xxz,   \
                             tr_z_xxz_xxzz,  \
                             tr_z_xxz_xyy,   \
                             tr_z_xxz_xyyy,  \
                             tr_z_xxz_xyyz,  \
                             tr_z_xxz_xyz,   \
                             tr_z_xxz_xyzz,  \
                             tr_z_xxz_xzz,   \
                             tr_z_xxz_xzzz,  \
                             tr_z_xxz_zzzz,  \
                             tr_z_xyz_yyyy,  \
                             tr_z_xyz_yyyz,  \
                             tr_z_xyz_yyzz,  \
                             tr_z_xyz_yzzz,  \
                             tr_z_yz_yyyy,   \
                             tr_z_yz_yyyz,   \
                             tr_z_yz_yyzz,   \
                             tr_z_yz_yzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyz_xxxx[i] = tr_z_xxz_xxxx[i] * pa_y[i];

        tr_z_xxyz_xxxy[i] = tr_z_xxz_xxx[i] * fe_0 + tr_z_xxz_xxxy[i] * pa_y[i];

        tr_z_xxyz_xxxz[i] = tr_z_xxz_xxxz[i] * pa_y[i];

        tr_z_xxyz_xxyy[i] = 2.0 * tr_z_xxz_xxy[i] * fe_0 + tr_z_xxz_xxyy[i] * pa_y[i];

        tr_z_xxyz_xxyz[i] = tr_z_xxz_xxz[i] * fe_0 + tr_z_xxz_xxyz[i] * pa_y[i];

        tr_z_xxyz_xxzz[i] = tr_z_xxz_xxzz[i] * pa_y[i];

        tr_z_xxyz_xyyy[i] = 3.0 * tr_z_xxz_xyy[i] * fe_0 + tr_z_xxz_xyyy[i] * pa_y[i];

        tr_z_xxyz_xyyz[i] = 2.0 * tr_z_xxz_xyz[i] * fe_0 + tr_z_xxz_xyyz[i] * pa_y[i];

        tr_z_xxyz_xyzz[i] = tr_z_xxz_xzz[i] * fe_0 + tr_z_xxz_xyzz[i] * pa_y[i];

        tr_z_xxyz_xzzz[i] = tr_z_xxz_xzzz[i] * pa_y[i];

        tr_z_xxyz_yyyy[i] = tr_z_yz_yyyy[i] * fe_0 + tr_z_xyz_yyyy[i] * pa_x[i];

        tr_z_xxyz_yyyz[i] = tr_z_yz_yyyz[i] * fe_0 + tr_z_xyz_yyyz[i] * pa_x[i];

        tr_z_xxyz_yyzz[i] = tr_z_yz_yyzz[i] * fe_0 + tr_z_xyz_yyzz[i] * pa_x[i];

        tr_z_xxyz_yzzz[i] = tr_z_yz_yzzz[i] * fe_0 + tr_z_xyz_yzzz[i] * pa_x[i];

        tr_z_xxyz_zzzz[i] = tr_z_xxz_zzzz[i] * pa_y[i];
    }

    // Set up 525-540 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,               \
                             tr_z_xxzz_xxxx, \
                             tr_z_xxzz_xxxy, \
                             tr_z_xxzz_xxxz, \
                             tr_z_xxzz_xxyy, \
                             tr_z_xxzz_xxyz, \
                             tr_z_xxzz_xxzz, \
                             tr_z_xxzz_xyyy, \
                             tr_z_xxzz_xyyz, \
                             tr_z_xxzz_xyzz, \
                             tr_z_xxzz_xzzz, \
                             tr_z_xxzz_yyyy, \
                             tr_z_xxzz_yyyz, \
                             tr_z_xxzz_yyzz, \
                             tr_z_xxzz_yzzz, \
                             tr_z_xxzz_zzzz, \
                             tr_z_xzz_xxx,   \
                             tr_z_xzz_xxxx,  \
                             tr_z_xzz_xxxy,  \
                             tr_z_xzz_xxxz,  \
                             tr_z_xzz_xxy,   \
                             tr_z_xzz_xxyy,  \
                             tr_z_xzz_xxyz,  \
                             tr_z_xzz_xxz,   \
                             tr_z_xzz_xxzz,  \
                             tr_z_xzz_xyy,   \
                             tr_z_xzz_xyyy,  \
                             tr_z_xzz_xyyz,  \
                             tr_z_xzz_xyz,   \
                             tr_z_xzz_xyzz,  \
                             tr_z_xzz_xzz,   \
                             tr_z_xzz_xzzz,  \
                             tr_z_xzz_yyy,   \
                             tr_z_xzz_yyyy,  \
                             tr_z_xzz_yyyz,  \
                             tr_z_xzz_yyz,   \
                             tr_z_xzz_yyzz,  \
                             tr_z_xzz_yzz,   \
                             tr_z_xzz_yzzz,  \
                             tr_z_xzz_zzz,   \
                             tr_z_xzz_zzzz,  \
                             tr_z_zz_xxxx,   \
                             tr_z_zz_xxxy,   \
                             tr_z_zz_xxxz,   \
                             tr_z_zz_xxyy,   \
                             tr_z_zz_xxyz,   \
                             tr_z_zz_xxzz,   \
                             tr_z_zz_xyyy,   \
                             tr_z_zz_xyyz,   \
                             tr_z_zz_xyzz,   \
                             tr_z_zz_xzzz,   \
                             tr_z_zz_yyyy,   \
                             tr_z_zz_yyyz,   \
                             tr_z_zz_yyzz,   \
                             tr_z_zz_yzzz,   \
                             tr_z_zz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzz_xxxx[i] = tr_z_zz_xxxx[i] * fe_0 + 4.0 * tr_z_xzz_xxx[i] * fe_0 + tr_z_xzz_xxxx[i] * pa_x[i];

        tr_z_xxzz_xxxy[i] = tr_z_zz_xxxy[i] * fe_0 + 3.0 * tr_z_xzz_xxy[i] * fe_0 + tr_z_xzz_xxxy[i] * pa_x[i];

        tr_z_xxzz_xxxz[i] = tr_z_zz_xxxz[i] * fe_0 + 3.0 * tr_z_xzz_xxz[i] * fe_0 + tr_z_xzz_xxxz[i] * pa_x[i];

        tr_z_xxzz_xxyy[i] = tr_z_zz_xxyy[i] * fe_0 + 2.0 * tr_z_xzz_xyy[i] * fe_0 + tr_z_xzz_xxyy[i] * pa_x[i];

        tr_z_xxzz_xxyz[i] = tr_z_zz_xxyz[i] * fe_0 + 2.0 * tr_z_xzz_xyz[i] * fe_0 + tr_z_xzz_xxyz[i] * pa_x[i];

        tr_z_xxzz_xxzz[i] = tr_z_zz_xxzz[i] * fe_0 + 2.0 * tr_z_xzz_xzz[i] * fe_0 + tr_z_xzz_xxzz[i] * pa_x[i];

        tr_z_xxzz_xyyy[i] = tr_z_zz_xyyy[i] * fe_0 + tr_z_xzz_yyy[i] * fe_0 + tr_z_xzz_xyyy[i] * pa_x[i];

        tr_z_xxzz_xyyz[i] = tr_z_zz_xyyz[i] * fe_0 + tr_z_xzz_yyz[i] * fe_0 + tr_z_xzz_xyyz[i] * pa_x[i];

        tr_z_xxzz_xyzz[i] = tr_z_zz_xyzz[i] * fe_0 + tr_z_xzz_yzz[i] * fe_0 + tr_z_xzz_xyzz[i] * pa_x[i];

        tr_z_xxzz_xzzz[i] = tr_z_zz_xzzz[i] * fe_0 + tr_z_xzz_zzz[i] * fe_0 + tr_z_xzz_xzzz[i] * pa_x[i];

        tr_z_xxzz_yyyy[i] = tr_z_zz_yyyy[i] * fe_0 + tr_z_xzz_yyyy[i] * pa_x[i];

        tr_z_xxzz_yyyz[i] = tr_z_zz_yyyz[i] * fe_0 + tr_z_xzz_yyyz[i] * pa_x[i];

        tr_z_xxzz_yyzz[i] = tr_z_zz_yyzz[i] * fe_0 + tr_z_xzz_yyzz[i] * pa_x[i];

        tr_z_xxzz_yzzz[i] = tr_z_zz_yzzz[i] * fe_0 + tr_z_xzz_yzzz[i] * pa_x[i];

        tr_z_xxzz_zzzz[i] = tr_z_zz_zzzz[i] * fe_0 + tr_z_xzz_zzzz[i] * pa_x[i];
    }

    // Set up 540-555 components of targeted buffer : GG

    auto tr_z_xyyy_xxxx = pbuffer.data(idx_dip_gg + 540);

    auto tr_z_xyyy_xxxy = pbuffer.data(idx_dip_gg + 541);

    auto tr_z_xyyy_xxxz = pbuffer.data(idx_dip_gg + 542);

    auto tr_z_xyyy_xxyy = pbuffer.data(idx_dip_gg + 543);

    auto tr_z_xyyy_xxyz = pbuffer.data(idx_dip_gg + 544);

    auto tr_z_xyyy_xxzz = pbuffer.data(idx_dip_gg + 545);

    auto tr_z_xyyy_xyyy = pbuffer.data(idx_dip_gg + 546);

    auto tr_z_xyyy_xyyz = pbuffer.data(idx_dip_gg + 547);

    auto tr_z_xyyy_xyzz = pbuffer.data(idx_dip_gg + 548);

    auto tr_z_xyyy_xzzz = pbuffer.data(idx_dip_gg + 549);

    auto tr_z_xyyy_yyyy = pbuffer.data(idx_dip_gg + 550);

    auto tr_z_xyyy_yyyz = pbuffer.data(idx_dip_gg + 551);

    auto tr_z_xyyy_yyzz = pbuffer.data(idx_dip_gg + 552);

    auto tr_z_xyyy_yzzz = pbuffer.data(idx_dip_gg + 553);

    auto tr_z_xyyy_zzzz = pbuffer.data(idx_dip_gg + 554);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xyyy_xxxx, \
                             tr_z_xyyy_xxxy, \
                             tr_z_xyyy_xxxz, \
                             tr_z_xyyy_xxyy, \
                             tr_z_xyyy_xxyz, \
                             tr_z_xyyy_xxzz, \
                             tr_z_xyyy_xyyy, \
                             tr_z_xyyy_xyyz, \
                             tr_z_xyyy_xyzz, \
                             tr_z_xyyy_xzzz, \
                             tr_z_xyyy_yyyy, \
                             tr_z_xyyy_yyyz, \
                             tr_z_xyyy_yyzz, \
                             tr_z_xyyy_yzzz, \
                             tr_z_xyyy_zzzz, \
                             tr_z_yyy_xxx,   \
                             tr_z_yyy_xxxx,  \
                             tr_z_yyy_xxxy,  \
                             tr_z_yyy_xxxz,  \
                             tr_z_yyy_xxy,   \
                             tr_z_yyy_xxyy,  \
                             tr_z_yyy_xxyz,  \
                             tr_z_yyy_xxz,   \
                             tr_z_yyy_xxzz,  \
                             tr_z_yyy_xyy,   \
                             tr_z_yyy_xyyy,  \
                             tr_z_yyy_xyyz,  \
                             tr_z_yyy_xyz,   \
                             tr_z_yyy_xyzz,  \
                             tr_z_yyy_xzz,   \
                             tr_z_yyy_xzzz,  \
                             tr_z_yyy_yyy,   \
                             tr_z_yyy_yyyy,  \
                             tr_z_yyy_yyyz,  \
                             tr_z_yyy_yyz,   \
                             tr_z_yyy_yyzz,  \
                             tr_z_yyy_yzz,   \
                             tr_z_yyy_yzzz,  \
                             tr_z_yyy_zzz,   \
                             tr_z_yyy_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyy_xxxx[i] = 4.0 * tr_z_yyy_xxx[i] * fe_0 + tr_z_yyy_xxxx[i] * pa_x[i];

        tr_z_xyyy_xxxy[i] = 3.0 * tr_z_yyy_xxy[i] * fe_0 + tr_z_yyy_xxxy[i] * pa_x[i];

        tr_z_xyyy_xxxz[i] = 3.0 * tr_z_yyy_xxz[i] * fe_0 + tr_z_yyy_xxxz[i] * pa_x[i];

        tr_z_xyyy_xxyy[i] = 2.0 * tr_z_yyy_xyy[i] * fe_0 + tr_z_yyy_xxyy[i] * pa_x[i];

        tr_z_xyyy_xxyz[i] = 2.0 * tr_z_yyy_xyz[i] * fe_0 + tr_z_yyy_xxyz[i] * pa_x[i];

        tr_z_xyyy_xxzz[i] = 2.0 * tr_z_yyy_xzz[i] * fe_0 + tr_z_yyy_xxzz[i] * pa_x[i];

        tr_z_xyyy_xyyy[i] = tr_z_yyy_yyy[i] * fe_0 + tr_z_yyy_xyyy[i] * pa_x[i];

        tr_z_xyyy_xyyz[i] = tr_z_yyy_yyz[i] * fe_0 + tr_z_yyy_xyyz[i] * pa_x[i];

        tr_z_xyyy_xyzz[i] = tr_z_yyy_yzz[i] * fe_0 + tr_z_yyy_xyzz[i] * pa_x[i];

        tr_z_xyyy_xzzz[i] = tr_z_yyy_zzz[i] * fe_0 + tr_z_yyy_xzzz[i] * pa_x[i];

        tr_z_xyyy_yyyy[i] = tr_z_yyy_yyyy[i] * pa_x[i];

        tr_z_xyyy_yyyz[i] = tr_z_yyy_yyyz[i] * pa_x[i];

        tr_z_xyyy_yyzz[i] = tr_z_yyy_yyzz[i] * pa_x[i];

        tr_z_xyyy_yzzz[i] = tr_z_yyy_yzzz[i] * pa_x[i];

        tr_z_xyyy_zzzz[i] = tr_z_yyy_zzzz[i] * pa_x[i];
    }

    // Set up 555-570 components of targeted buffer : GG

    auto tr_z_xyyz_xxxx = pbuffer.data(idx_dip_gg + 555);

    auto tr_z_xyyz_xxxy = pbuffer.data(idx_dip_gg + 556);

    auto tr_z_xyyz_xxxz = pbuffer.data(idx_dip_gg + 557);

    auto tr_z_xyyz_xxyy = pbuffer.data(idx_dip_gg + 558);

    auto tr_z_xyyz_xxyz = pbuffer.data(idx_dip_gg + 559);

    auto tr_z_xyyz_xxzz = pbuffer.data(idx_dip_gg + 560);

    auto tr_z_xyyz_xyyy = pbuffer.data(idx_dip_gg + 561);

    auto tr_z_xyyz_xyyz = pbuffer.data(idx_dip_gg + 562);

    auto tr_z_xyyz_xyzz = pbuffer.data(idx_dip_gg + 563);

    auto tr_z_xyyz_xzzz = pbuffer.data(idx_dip_gg + 564);

    auto tr_z_xyyz_yyyy = pbuffer.data(idx_dip_gg + 565);

    auto tr_z_xyyz_yyyz = pbuffer.data(idx_dip_gg + 566);

    auto tr_z_xyyz_yyzz = pbuffer.data(idx_dip_gg + 567);

    auto tr_z_xyyz_yzzz = pbuffer.data(idx_dip_gg + 568);

    auto tr_z_xyyz_zzzz = pbuffer.data(idx_dip_gg + 569);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xyyz_xxxx, \
                             tr_z_xyyz_xxxy, \
                             tr_z_xyyz_xxxz, \
                             tr_z_xyyz_xxyy, \
                             tr_z_xyyz_xxyz, \
                             tr_z_xyyz_xxzz, \
                             tr_z_xyyz_xyyy, \
                             tr_z_xyyz_xyyz, \
                             tr_z_xyyz_xyzz, \
                             tr_z_xyyz_xzzz, \
                             tr_z_xyyz_yyyy, \
                             tr_z_xyyz_yyyz, \
                             tr_z_xyyz_yyzz, \
                             tr_z_xyyz_yzzz, \
                             tr_z_xyyz_zzzz, \
                             tr_z_yyz_xxx,   \
                             tr_z_yyz_xxxx,  \
                             tr_z_yyz_xxxy,  \
                             tr_z_yyz_xxxz,  \
                             tr_z_yyz_xxy,   \
                             tr_z_yyz_xxyy,  \
                             tr_z_yyz_xxyz,  \
                             tr_z_yyz_xxz,   \
                             tr_z_yyz_xxzz,  \
                             tr_z_yyz_xyy,   \
                             tr_z_yyz_xyyy,  \
                             tr_z_yyz_xyyz,  \
                             tr_z_yyz_xyz,   \
                             tr_z_yyz_xyzz,  \
                             tr_z_yyz_xzz,   \
                             tr_z_yyz_xzzz,  \
                             tr_z_yyz_yyy,   \
                             tr_z_yyz_yyyy,  \
                             tr_z_yyz_yyyz,  \
                             tr_z_yyz_yyz,   \
                             tr_z_yyz_yyzz,  \
                             tr_z_yyz_yzz,   \
                             tr_z_yyz_yzzz,  \
                             tr_z_yyz_zzz,   \
                             tr_z_yyz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyz_xxxx[i] = 4.0 * tr_z_yyz_xxx[i] * fe_0 + tr_z_yyz_xxxx[i] * pa_x[i];

        tr_z_xyyz_xxxy[i] = 3.0 * tr_z_yyz_xxy[i] * fe_0 + tr_z_yyz_xxxy[i] * pa_x[i];

        tr_z_xyyz_xxxz[i] = 3.0 * tr_z_yyz_xxz[i] * fe_0 + tr_z_yyz_xxxz[i] * pa_x[i];

        tr_z_xyyz_xxyy[i] = 2.0 * tr_z_yyz_xyy[i] * fe_0 + tr_z_yyz_xxyy[i] * pa_x[i];

        tr_z_xyyz_xxyz[i] = 2.0 * tr_z_yyz_xyz[i] * fe_0 + tr_z_yyz_xxyz[i] * pa_x[i];

        tr_z_xyyz_xxzz[i] = 2.0 * tr_z_yyz_xzz[i] * fe_0 + tr_z_yyz_xxzz[i] * pa_x[i];

        tr_z_xyyz_xyyy[i] = tr_z_yyz_yyy[i] * fe_0 + tr_z_yyz_xyyy[i] * pa_x[i];

        tr_z_xyyz_xyyz[i] = tr_z_yyz_yyz[i] * fe_0 + tr_z_yyz_xyyz[i] * pa_x[i];

        tr_z_xyyz_xyzz[i] = tr_z_yyz_yzz[i] * fe_0 + tr_z_yyz_xyzz[i] * pa_x[i];

        tr_z_xyyz_xzzz[i] = tr_z_yyz_zzz[i] * fe_0 + tr_z_yyz_xzzz[i] * pa_x[i];

        tr_z_xyyz_yyyy[i] = tr_z_yyz_yyyy[i] * pa_x[i];

        tr_z_xyyz_yyyz[i] = tr_z_yyz_yyyz[i] * pa_x[i];

        tr_z_xyyz_yyzz[i] = tr_z_yyz_yyzz[i] * pa_x[i];

        tr_z_xyyz_yzzz[i] = tr_z_yyz_yzzz[i] * pa_x[i];

        tr_z_xyyz_zzzz[i] = tr_z_yyz_zzzz[i] * pa_x[i];
    }

    // Set up 570-585 components of targeted buffer : GG

    auto tr_z_xyzz_xxxx = pbuffer.data(idx_dip_gg + 570);

    auto tr_z_xyzz_xxxy = pbuffer.data(idx_dip_gg + 571);

    auto tr_z_xyzz_xxxz = pbuffer.data(idx_dip_gg + 572);

    auto tr_z_xyzz_xxyy = pbuffer.data(idx_dip_gg + 573);

    auto tr_z_xyzz_xxyz = pbuffer.data(idx_dip_gg + 574);

    auto tr_z_xyzz_xxzz = pbuffer.data(idx_dip_gg + 575);

    auto tr_z_xyzz_xyyy = pbuffer.data(idx_dip_gg + 576);

    auto tr_z_xyzz_xyyz = pbuffer.data(idx_dip_gg + 577);

    auto tr_z_xyzz_xyzz = pbuffer.data(idx_dip_gg + 578);

    auto tr_z_xyzz_xzzz = pbuffer.data(idx_dip_gg + 579);

    auto tr_z_xyzz_yyyy = pbuffer.data(idx_dip_gg + 580);

    auto tr_z_xyzz_yyyz = pbuffer.data(idx_dip_gg + 581);

    auto tr_z_xyzz_yyzz = pbuffer.data(idx_dip_gg + 582);

    auto tr_z_xyzz_yzzz = pbuffer.data(idx_dip_gg + 583);

    auto tr_z_xyzz_zzzz = pbuffer.data(idx_dip_gg + 584);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xyzz_xxxx, \
                             tr_z_xyzz_xxxy, \
                             tr_z_xyzz_xxxz, \
                             tr_z_xyzz_xxyy, \
                             tr_z_xyzz_xxyz, \
                             tr_z_xyzz_xxzz, \
                             tr_z_xyzz_xyyy, \
                             tr_z_xyzz_xyyz, \
                             tr_z_xyzz_xyzz, \
                             tr_z_xyzz_xzzz, \
                             tr_z_xyzz_yyyy, \
                             tr_z_xyzz_yyyz, \
                             tr_z_xyzz_yyzz, \
                             tr_z_xyzz_yzzz, \
                             tr_z_xyzz_zzzz, \
                             tr_z_xzz_xxxx,  \
                             tr_z_xzz_xxxz,  \
                             tr_z_xzz_xxzz,  \
                             tr_z_xzz_xzzz,  \
                             tr_z_yzz_xxxy,  \
                             tr_z_yzz_xxy,   \
                             tr_z_yzz_xxyy,  \
                             tr_z_yzz_xxyz,  \
                             tr_z_yzz_xyy,   \
                             tr_z_yzz_xyyy,  \
                             tr_z_yzz_xyyz,  \
                             tr_z_yzz_xyz,   \
                             tr_z_yzz_xyzz,  \
                             tr_z_yzz_yyy,   \
                             tr_z_yzz_yyyy,  \
                             tr_z_yzz_yyyz,  \
                             tr_z_yzz_yyz,   \
                             tr_z_yzz_yyzz,  \
                             tr_z_yzz_yzz,   \
                             tr_z_yzz_yzzz,  \
                             tr_z_yzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyzz_xxxx[i] = tr_z_xzz_xxxx[i] * pa_y[i];

        tr_z_xyzz_xxxy[i] = 3.0 * tr_z_yzz_xxy[i] * fe_0 + tr_z_yzz_xxxy[i] * pa_x[i];

        tr_z_xyzz_xxxz[i] = tr_z_xzz_xxxz[i] * pa_y[i];

        tr_z_xyzz_xxyy[i] = 2.0 * tr_z_yzz_xyy[i] * fe_0 + tr_z_yzz_xxyy[i] * pa_x[i];

        tr_z_xyzz_xxyz[i] = 2.0 * tr_z_yzz_xyz[i] * fe_0 + tr_z_yzz_xxyz[i] * pa_x[i];

        tr_z_xyzz_xxzz[i] = tr_z_xzz_xxzz[i] * pa_y[i];

        tr_z_xyzz_xyyy[i] = tr_z_yzz_yyy[i] * fe_0 + tr_z_yzz_xyyy[i] * pa_x[i];

        tr_z_xyzz_xyyz[i] = tr_z_yzz_yyz[i] * fe_0 + tr_z_yzz_xyyz[i] * pa_x[i];

        tr_z_xyzz_xyzz[i] = tr_z_yzz_yzz[i] * fe_0 + tr_z_yzz_xyzz[i] * pa_x[i];

        tr_z_xyzz_xzzz[i] = tr_z_xzz_xzzz[i] * pa_y[i];

        tr_z_xyzz_yyyy[i] = tr_z_yzz_yyyy[i] * pa_x[i];

        tr_z_xyzz_yyyz[i] = tr_z_yzz_yyyz[i] * pa_x[i];

        tr_z_xyzz_yyzz[i] = tr_z_yzz_yyzz[i] * pa_x[i];

        tr_z_xyzz_yzzz[i] = tr_z_yzz_yzzz[i] * pa_x[i];

        tr_z_xyzz_zzzz[i] = tr_z_yzz_zzzz[i] * pa_x[i];
    }

    // Set up 585-600 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_x,               \
                             tr_z_xzzz_xxxx, \
                             tr_z_xzzz_xxxy, \
                             tr_z_xzzz_xxxz, \
                             tr_z_xzzz_xxyy, \
                             tr_z_xzzz_xxyz, \
                             tr_z_xzzz_xxzz, \
                             tr_z_xzzz_xyyy, \
                             tr_z_xzzz_xyyz, \
                             tr_z_xzzz_xyzz, \
                             tr_z_xzzz_xzzz, \
                             tr_z_xzzz_yyyy, \
                             tr_z_xzzz_yyyz, \
                             tr_z_xzzz_yyzz, \
                             tr_z_xzzz_yzzz, \
                             tr_z_xzzz_zzzz, \
                             tr_z_zzz_xxx,   \
                             tr_z_zzz_xxxx,  \
                             tr_z_zzz_xxxy,  \
                             tr_z_zzz_xxxz,  \
                             tr_z_zzz_xxy,   \
                             tr_z_zzz_xxyy,  \
                             tr_z_zzz_xxyz,  \
                             tr_z_zzz_xxz,   \
                             tr_z_zzz_xxzz,  \
                             tr_z_zzz_xyy,   \
                             tr_z_zzz_xyyy,  \
                             tr_z_zzz_xyyz,  \
                             tr_z_zzz_xyz,   \
                             tr_z_zzz_xyzz,  \
                             tr_z_zzz_xzz,   \
                             tr_z_zzz_xzzz,  \
                             tr_z_zzz_yyy,   \
                             tr_z_zzz_yyyy,  \
                             tr_z_zzz_yyyz,  \
                             tr_z_zzz_yyz,   \
                             tr_z_zzz_yyzz,  \
                             tr_z_zzz_yzz,   \
                             tr_z_zzz_yzzz,  \
                             tr_z_zzz_zzz,   \
                             tr_z_zzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzz_xxxx[i] = 4.0 * tr_z_zzz_xxx[i] * fe_0 + tr_z_zzz_xxxx[i] * pa_x[i];

        tr_z_xzzz_xxxy[i] = 3.0 * tr_z_zzz_xxy[i] * fe_0 + tr_z_zzz_xxxy[i] * pa_x[i];

        tr_z_xzzz_xxxz[i] = 3.0 * tr_z_zzz_xxz[i] * fe_0 + tr_z_zzz_xxxz[i] * pa_x[i];

        tr_z_xzzz_xxyy[i] = 2.0 * tr_z_zzz_xyy[i] * fe_0 + tr_z_zzz_xxyy[i] * pa_x[i];

        tr_z_xzzz_xxyz[i] = 2.0 * tr_z_zzz_xyz[i] * fe_0 + tr_z_zzz_xxyz[i] * pa_x[i];

        tr_z_xzzz_xxzz[i] = 2.0 * tr_z_zzz_xzz[i] * fe_0 + tr_z_zzz_xxzz[i] * pa_x[i];

        tr_z_xzzz_xyyy[i] = tr_z_zzz_yyy[i] * fe_0 + tr_z_zzz_xyyy[i] * pa_x[i];

        tr_z_xzzz_xyyz[i] = tr_z_zzz_yyz[i] * fe_0 + tr_z_zzz_xyyz[i] * pa_x[i];

        tr_z_xzzz_xyzz[i] = tr_z_zzz_yzz[i] * fe_0 + tr_z_zzz_xyzz[i] * pa_x[i];

        tr_z_xzzz_xzzz[i] = tr_z_zzz_zzz[i] * fe_0 + tr_z_zzz_xzzz[i] * pa_x[i];

        tr_z_xzzz_yyyy[i] = tr_z_zzz_yyyy[i] * pa_x[i];

        tr_z_xzzz_yyyz[i] = tr_z_zzz_yyyz[i] * pa_x[i];

        tr_z_xzzz_yyzz[i] = tr_z_zzz_yyzz[i] * pa_x[i];

        tr_z_xzzz_yzzz[i] = tr_z_zzz_yzzz[i] * pa_x[i];

        tr_z_xzzz_zzzz[i] = tr_z_zzz_zzzz[i] * pa_x[i];
    }

    // Set up 600-615 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_y,               \
                             tr_z_yy_xxxx,   \
                             tr_z_yy_xxxy,   \
                             tr_z_yy_xxxz,   \
                             tr_z_yy_xxyy,   \
                             tr_z_yy_xxyz,   \
                             tr_z_yy_xxzz,   \
                             tr_z_yy_xyyy,   \
                             tr_z_yy_xyyz,   \
                             tr_z_yy_xyzz,   \
                             tr_z_yy_xzzz,   \
                             tr_z_yy_yyyy,   \
                             tr_z_yy_yyyz,   \
                             tr_z_yy_yyzz,   \
                             tr_z_yy_yzzz,   \
                             tr_z_yy_zzzz,   \
                             tr_z_yyy_xxx,   \
                             tr_z_yyy_xxxx,  \
                             tr_z_yyy_xxxy,  \
                             tr_z_yyy_xxxz,  \
                             tr_z_yyy_xxy,   \
                             tr_z_yyy_xxyy,  \
                             tr_z_yyy_xxyz,  \
                             tr_z_yyy_xxz,   \
                             tr_z_yyy_xxzz,  \
                             tr_z_yyy_xyy,   \
                             tr_z_yyy_xyyy,  \
                             tr_z_yyy_xyyz,  \
                             tr_z_yyy_xyz,   \
                             tr_z_yyy_xyzz,  \
                             tr_z_yyy_xzz,   \
                             tr_z_yyy_xzzz,  \
                             tr_z_yyy_yyy,   \
                             tr_z_yyy_yyyy,  \
                             tr_z_yyy_yyyz,  \
                             tr_z_yyy_yyz,   \
                             tr_z_yyy_yyzz,  \
                             tr_z_yyy_yzz,   \
                             tr_z_yyy_yzzz,  \
                             tr_z_yyy_zzz,   \
                             tr_z_yyy_zzzz,  \
                             tr_z_yyyy_xxxx, \
                             tr_z_yyyy_xxxy, \
                             tr_z_yyyy_xxxz, \
                             tr_z_yyyy_xxyy, \
                             tr_z_yyyy_xxyz, \
                             tr_z_yyyy_xxzz, \
                             tr_z_yyyy_xyyy, \
                             tr_z_yyyy_xyyz, \
                             tr_z_yyyy_xyzz, \
                             tr_z_yyyy_xzzz, \
                             tr_z_yyyy_yyyy, \
                             tr_z_yyyy_yyyz, \
                             tr_z_yyyy_yyzz, \
                             tr_z_yyyy_yzzz, \
                             tr_z_yyyy_zzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyy_xxxx[i] = 3.0 * tr_z_yy_xxxx[i] * fe_0 + tr_z_yyy_xxxx[i] * pa_y[i];

        tr_z_yyyy_xxxy[i] = 3.0 * tr_z_yy_xxxy[i] * fe_0 + tr_z_yyy_xxx[i] * fe_0 + tr_z_yyy_xxxy[i] * pa_y[i];

        tr_z_yyyy_xxxz[i] = 3.0 * tr_z_yy_xxxz[i] * fe_0 + tr_z_yyy_xxxz[i] * pa_y[i];

        tr_z_yyyy_xxyy[i] = 3.0 * tr_z_yy_xxyy[i] * fe_0 + 2.0 * tr_z_yyy_xxy[i] * fe_0 + tr_z_yyy_xxyy[i] * pa_y[i];

        tr_z_yyyy_xxyz[i] = 3.0 * tr_z_yy_xxyz[i] * fe_0 + tr_z_yyy_xxz[i] * fe_0 + tr_z_yyy_xxyz[i] * pa_y[i];

        tr_z_yyyy_xxzz[i] = 3.0 * tr_z_yy_xxzz[i] * fe_0 + tr_z_yyy_xxzz[i] * pa_y[i];

        tr_z_yyyy_xyyy[i] = 3.0 * tr_z_yy_xyyy[i] * fe_0 + 3.0 * tr_z_yyy_xyy[i] * fe_0 + tr_z_yyy_xyyy[i] * pa_y[i];

        tr_z_yyyy_xyyz[i] = 3.0 * tr_z_yy_xyyz[i] * fe_0 + 2.0 * tr_z_yyy_xyz[i] * fe_0 + tr_z_yyy_xyyz[i] * pa_y[i];

        tr_z_yyyy_xyzz[i] = 3.0 * tr_z_yy_xyzz[i] * fe_0 + tr_z_yyy_xzz[i] * fe_0 + tr_z_yyy_xyzz[i] * pa_y[i];

        tr_z_yyyy_xzzz[i] = 3.0 * tr_z_yy_xzzz[i] * fe_0 + tr_z_yyy_xzzz[i] * pa_y[i];

        tr_z_yyyy_yyyy[i] = 3.0 * tr_z_yy_yyyy[i] * fe_0 + 4.0 * tr_z_yyy_yyy[i] * fe_0 + tr_z_yyy_yyyy[i] * pa_y[i];

        tr_z_yyyy_yyyz[i] = 3.0 * tr_z_yy_yyyz[i] * fe_0 + 3.0 * tr_z_yyy_yyz[i] * fe_0 + tr_z_yyy_yyyz[i] * pa_y[i];

        tr_z_yyyy_yyzz[i] = 3.0 * tr_z_yy_yyzz[i] * fe_0 + 2.0 * tr_z_yyy_yzz[i] * fe_0 + tr_z_yyy_yyzz[i] * pa_y[i];

        tr_z_yyyy_yzzz[i] = 3.0 * tr_z_yy_yzzz[i] * fe_0 + tr_z_yyy_zzz[i] * fe_0 + tr_z_yyy_yzzz[i] * pa_y[i];

        tr_z_yyyy_zzzz[i] = 3.0 * tr_z_yy_zzzz[i] * fe_0 + tr_z_yyy_zzzz[i] * pa_y[i];
    }

    // Set up 615-630 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_z_yyy_xxxy,  \
                             tr_z_yyy_xxyy,  \
                             tr_z_yyy_xyyy,  \
                             tr_z_yyy_yyyy,  \
                             tr_z_yyyz_xxxx, \
                             tr_z_yyyz_xxxy, \
                             tr_z_yyyz_xxxz, \
                             tr_z_yyyz_xxyy, \
                             tr_z_yyyz_xxyz, \
                             tr_z_yyyz_xxzz, \
                             tr_z_yyyz_xyyy, \
                             tr_z_yyyz_xyyz, \
                             tr_z_yyyz_xyzz, \
                             tr_z_yyyz_xzzz, \
                             tr_z_yyyz_yyyy, \
                             tr_z_yyyz_yyyz, \
                             tr_z_yyyz_yyzz, \
                             tr_z_yyyz_yzzz, \
                             tr_z_yyyz_zzzz, \
                             tr_z_yyz_xxxx,  \
                             tr_z_yyz_xxxz,  \
                             tr_z_yyz_xxyz,  \
                             tr_z_yyz_xxz,   \
                             tr_z_yyz_xxzz,  \
                             tr_z_yyz_xyyz,  \
                             tr_z_yyz_xyz,   \
                             tr_z_yyz_xyzz,  \
                             tr_z_yyz_xzz,   \
                             tr_z_yyz_xzzz,  \
                             tr_z_yyz_yyyz,  \
                             tr_z_yyz_yyz,   \
                             tr_z_yyz_yyzz,  \
                             tr_z_yyz_yzz,   \
                             tr_z_yyz_yzzz,  \
                             tr_z_yyz_zzz,   \
                             tr_z_yyz_zzzz,  \
                             tr_z_yz_xxxx,   \
                             tr_z_yz_xxxz,   \
                             tr_z_yz_xxyz,   \
                             tr_z_yz_xxzz,   \
                             tr_z_yz_xyyz,   \
                             tr_z_yz_xyzz,   \
                             tr_z_yz_xzzz,   \
                             tr_z_yz_yyyz,   \
                             tr_z_yz_yyzz,   \
                             tr_z_yz_yzzz,   \
                             tr_z_yz_zzzz,   \
                             ts_yyy_xxxy,    \
                             ts_yyy_xxyy,    \
                             ts_yyy_xyyy,    \
                             ts_yyy_yyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyz_xxxx[i] = 2.0 * tr_z_yz_xxxx[i] * fe_0 + tr_z_yyz_xxxx[i] * pa_y[i];

        tr_z_yyyz_xxxy[i] = ts_yyy_xxxy[i] * fe_0 + tr_z_yyy_xxxy[i] * pa_z[i];

        tr_z_yyyz_xxxz[i] = 2.0 * tr_z_yz_xxxz[i] * fe_0 + tr_z_yyz_xxxz[i] * pa_y[i];

        tr_z_yyyz_xxyy[i] = ts_yyy_xxyy[i] * fe_0 + tr_z_yyy_xxyy[i] * pa_z[i];

        tr_z_yyyz_xxyz[i] = 2.0 * tr_z_yz_xxyz[i] * fe_0 + tr_z_yyz_xxz[i] * fe_0 + tr_z_yyz_xxyz[i] * pa_y[i];

        tr_z_yyyz_xxzz[i] = 2.0 * tr_z_yz_xxzz[i] * fe_0 + tr_z_yyz_xxzz[i] * pa_y[i];

        tr_z_yyyz_xyyy[i] = ts_yyy_xyyy[i] * fe_0 + tr_z_yyy_xyyy[i] * pa_z[i];

        tr_z_yyyz_xyyz[i] = 2.0 * tr_z_yz_xyyz[i] * fe_0 + 2.0 * tr_z_yyz_xyz[i] * fe_0 + tr_z_yyz_xyyz[i] * pa_y[i];

        tr_z_yyyz_xyzz[i] = 2.0 * tr_z_yz_xyzz[i] * fe_0 + tr_z_yyz_xzz[i] * fe_0 + tr_z_yyz_xyzz[i] * pa_y[i];

        tr_z_yyyz_xzzz[i] = 2.0 * tr_z_yz_xzzz[i] * fe_0 + tr_z_yyz_xzzz[i] * pa_y[i];

        tr_z_yyyz_yyyy[i] = ts_yyy_yyyy[i] * fe_0 + tr_z_yyy_yyyy[i] * pa_z[i];

        tr_z_yyyz_yyyz[i] = 2.0 * tr_z_yz_yyyz[i] * fe_0 + 3.0 * tr_z_yyz_yyz[i] * fe_0 + tr_z_yyz_yyyz[i] * pa_y[i];

        tr_z_yyyz_yyzz[i] = 2.0 * tr_z_yz_yyzz[i] * fe_0 + 2.0 * tr_z_yyz_yzz[i] * fe_0 + tr_z_yyz_yyzz[i] * pa_y[i];

        tr_z_yyyz_yzzz[i] = 2.0 * tr_z_yz_yzzz[i] * fe_0 + tr_z_yyz_zzz[i] * fe_0 + tr_z_yyz_yzzz[i] * pa_y[i];

        tr_z_yyyz_zzzz[i] = 2.0 * tr_z_yz_zzzz[i] * fe_0 + tr_z_yyz_zzzz[i] * pa_y[i];
    }

    // Set up 630-645 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_y,               \
                             tr_z_yyzz_xxxx, \
                             tr_z_yyzz_xxxy, \
                             tr_z_yyzz_xxxz, \
                             tr_z_yyzz_xxyy, \
                             tr_z_yyzz_xxyz, \
                             tr_z_yyzz_xxzz, \
                             tr_z_yyzz_xyyy, \
                             tr_z_yyzz_xyyz, \
                             tr_z_yyzz_xyzz, \
                             tr_z_yyzz_xzzz, \
                             tr_z_yyzz_yyyy, \
                             tr_z_yyzz_yyyz, \
                             tr_z_yyzz_yyzz, \
                             tr_z_yyzz_yzzz, \
                             tr_z_yyzz_zzzz, \
                             tr_z_yzz_xxx,   \
                             tr_z_yzz_xxxx,  \
                             tr_z_yzz_xxxy,  \
                             tr_z_yzz_xxxz,  \
                             tr_z_yzz_xxy,   \
                             tr_z_yzz_xxyy,  \
                             tr_z_yzz_xxyz,  \
                             tr_z_yzz_xxz,   \
                             tr_z_yzz_xxzz,  \
                             tr_z_yzz_xyy,   \
                             tr_z_yzz_xyyy,  \
                             tr_z_yzz_xyyz,  \
                             tr_z_yzz_xyz,   \
                             tr_z_yzz_xyzz,  \
                             tr_z_yzz_xzz,   \
                             tr_z_yzz_xzzz,  \
                             tr_z_yzz_yyy,   \
                             tr_z_yzz_yyyy,  \
                             tr_z_yzz_yyyz,  \
                             tr_z_yzz_yyz,   \
                             tr_z_yzz_yyzz,  \
                             tr_z_yzz_yzz,   \
                             tr_z_yzz_yzzz,  \
                             tr_z_yzz_zzz,   \
                             tr_z_yzz_zzzz,  \
                             tr_z_zz_xxxx,   \
                             tr_z_zz_xxxy,   \
                             tr_z_zz_xxxz,   \
                             tr_z_zz_xxyy,   \
                             tr_z_zz_xxyz,   \
                             tr_z_zz_xxzz,   \
                             tr_z_zz_xyyy,   \
                             tr_z_zz_xyyz,   \
                             tr_z_zz_xyzz,   \
                             tr_z_zz_xzzz,   \
                             tr_z_zz_yyyy,   \
                             tr_z_zz_yyyz,   \
                             tr_z_zz_yyzz,   \
                             tr_z_zz_yzzz,   \
                             tr_z_zz_zzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzz_xxxx[i] = tr_z_zz_xxxx[i] * fe_0 + tr_z_yzz_xxxx[i] * pa_y[i];

        tr_z_yyzz_xxxy[i] = tr_z_zz_xxxy[i] * fe_0 + tr_z_yzz_xxx[i] * fe_0 + tr_z_yzz_xxxy[i] * pa_y[i];

        tr_z_yyzz_xxxz[i] = tr_z_zz_xxxz[i] * fe_0 + tr_z_yzz_xxxz[i] * pa_y[i];

        tr_z_yyzz_xxyy[i] = tr_z_zz_xxyy[i] * fe_0 + 2.0 * tr_z_yzz_xxy[i] * fe_0 + tr_z_yzz_xxyy[i] * pa_y[i];

        tr_z_yyzz_xxyz[i] = tr_z_zz_xxyz[i] * fe_0 + tr_z_yzz_xxz[i] * fe_0 + tr_z_yzz_xxyz[i] * pa_y[i];

        tr_z_yyzz_xxzz[i] = tr_z_zz_xxzz[i] * fe_0 + tr_z_yzz_xxzz[i] * pa_y[i];

        tr_z_yyzz_xyyy[i] = tr_z_zz_xyyy[i] * fe_0 + 3.0 * tr_z_yzz_xyy[i] * fe_0 + tr_z_yzz_xyyy[i] * pa_y[i];

        tr_z_yyzz_xyyz[i] = tr_z_zz_xyyz[i] * fe_0 + 2.0 * tr_z_yzz_xyz[i] * fe_0 + tr_z_yzz_xyyz[i] * pa_y[i];

        tr_z_yyzz_xyzz[i] = tr_z_zz_xyzz[i] * fe_0 + tr_z_yzz_xzz[i] * fe_0 + tr_z_yzz_xyzz[i] * pa_y[i];

        tr_z_yyzz_xzzz[i] = tr_z_zz_xzzz[i] * fe_0 + tr_z_yzz_xzzz[i] * pa_y[i];

        tr_z_yyzz_yyyy[i] = tr_z_zz_yyyy[i] * fe_0 + 4.0 * tr_z_yzz_yyy[i] * fe_0 + tr_z_yzz_yyyy[i] * pa_y[i];

        tr_z_yyzz_yyyz[i] = tr_z_zz_yyyz[i] * fe_0 + 3.0 * tr_z_yzz_yyz[i] * fe_0 + tr_z_yzz_yyyz[i] * pa_y[i];

        tr_z_yyzz_yyzz[i] = tr_z_zz_yyzz[i] * fe_0 + 2.0 * tr_z_yzz_yzz[i] * fe_0 + tr_z_yzz_yyzz[i] * pa_y[i];

        tr_z_yyzz_yzzz[i] = tr_z_zz_yzzz[i] * fe_0 + tr_z_yzz_zzz[i] * fe_0 + tr_z_yzz_yzzz[i] * pa_y[i];

        tr_z_yyzz_zzzz[i] = tr_z_zz_zzzz[i] * fe_0 + tr_z_yzz_zzzz[i] * pa_y[i];
    }

    // Set up 645-660 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_y,               \
                             tr_z_yzzz_xxxx, \
                             tr_z_yzzz_xxxy, \
                             tr_z_yzzz_xxxz, \
                             tr_z_yzzz_xxyy, \
                             tr_z_yzzz_xxyz, \
                             tr_z_yzzz_xxzz, \
                             tr_z_yzzz_xyyy, \
                             tr_z_yzzz_xyyz, \
                             tr_z_yzzz_xyzz, \
                             tr_z_yzzz_xzzz, \
                             tr_z_yzzz_yyyy, \
                             tr_z_yzzz_yyyz, \
                             tr_z_yzzz_yyzz, \
                             tr_z_yzzz_yzzz, \
                             tr_z_yzzz_zzzz, \
                             tr_z_zzz_xxx,   \
                             tr_z_zzz_xxxx,  \
                             tr_z_zzz_xxxy,  \
                             tr_z_zzz_xxxz,  \
                             tr_z_zzz_xxy,   \
                             tr_z_zzz_xxyy,  \
                             tr_z_zzz_xxyz,  \
                             tr_z_zzz_xxz,   \
                             tr_z_zzz_xxzz,  \
                             tr_z_zzz_xyy,   \
                             tr_z_zzz_xyyy,  \
                             tr_z_zzz_xyyz,  \
                             tr_z_zzz_xyz,   \
                             tr_z_zzz_xyzz,  \
                             tr_z_zzz_xzz,   \
                             tr_z_zzz_xzzz,  \
                             tr_z_zzz_yyy,   \
                             tr_z_zzz_yyyy,  \
                             tr_z_zzz_yyyz,  \
                             tr_z_zzz_yyz,   \
                             tr_z_zzz_yyzz,  \
                             tr_z_zzz_yzz,   \
                             tr_z_zzz_yzzz,  \
                             tr_z_zzz_zzz,   \
                             tr_z_zzz_zzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzz_xxxx[i] = tr_z_zzz_xxxx[i] * pa_y[i];

        tr_z_yzzz_xxxy[i] = tr_z_zzz_xxx[i] * fe_0 + tr_z_zzz_xxxy[i] * pa_y[i];

        tr_z_yzzz_xxxz[i] = tr_z_zzz_xxxz[i] * pa_y[i];

        tr_z_yzzz_xxyy[i] = 2.0 * tr_z_zzz_xxy[i] * fe_0 + tr_z_zzz_xxyy[i] * pa_y[i];

        tr_z_yzzz_xxyz[i] = tr_z_zzz_xxz[i] * fe_0 + tr_z_zzz_xxyz[i] * pa_y[i];

        tr_z_yzzz_xxzz[i] = tr_z_zzz_xxzz[i] * pa_y[i];

        tr_z_yzzz_xyyy[i] = 3.0 * tr_z_zzz_xyy[i] * fe_0 + tr_z_zzz_xyyy[i] * pa_y[i];

        tr_z_yzzz_xyyz[i] = 2.0 * tr_z_zzz_xyz[i] * fe_0 + tr_z_zzz_xyyz[i] * pa_y[i];

        tr_z_yzzz_xyzz[i] = tr_z_zzz_xzz[i] * fe_0 + tr_z_zzz_xyzz[i] * pa_y[i];

        tr_z_yzzz_xzzz[i] = tr_z_zzz_xzzz[i] * pa_y[i];

        tr_z_yzzz_yyyy[i] = 4.0 * tr_z_zzz_yyy[i] * fe_0 + tr_z_zzz_yyyy[i] * pa_y[i];

        tr_z_yzzz_yyyz[i] = 3.0 * tr_z_zzz_yyz[i] * fe_0 + tr_z_zzz_yyyz[i] * pa_y[i];

        tr_z_yzzz_yyzz[i] = 2.0 * tr_z_zzz_yzz[i] * fe_0 + tr_z_zzz_yyzz[i] * pa_y[i];

        tr_z_yzzz_yzzz[i] = tr_z_zzz_zzz[i] * fe_0 + tr_z_zzz_yzzz[i] * pa_y[i];

        tr_z_yzzz_zzzz[i] = tr_z_zzz_zzzz[i] * pa_y[i];
    }

    // Set up 660-675 components of targeted buffer : GG

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

#pragma omp simd aligned(pa_z,               \
                             tr_z_zz_xxxx,   \
                             tr_z_zz_xxxy,   \
                             tr_z_zz_xxxz,   \
                             tr_z_zz_xxyy,   \
                             tr_z_zz_xxyz,   \
                             tr_z_zz_xxzz,   \
                             tr_z_zz_xyyy,   \
                             tr_z_zz_xyyz,   \
                             tr_z_zz_xyzz,   \
                             tr_z_zz_xzzz,   \
                             tr_z_zz_yyyy,   \
                             tr_z_zz_yyyz,   \
                             tr_z_zz_yyzz,   \
                             tr_z_zz_yzzz,   \
                             tr_z_zz_zzzz,   \
                             tr_z_zzz_xxx,   \
                             tr_z_zzz_xxxx,  \
                             tr_z_zzz_xxxy,  \
                             tr_z_zzz_xxxz,  \
                             tr_z_zzz_xxy,   \
                             tr_z_zzz_xxyy,  \
                             tr_z_zzz_xxyz,  \
                             tr_z_zzz_xxz,   \
                             tr_z_zzz_xxzz,  \
                             tr_z_zzz_xyy,   \
                             tr_z_zzz_xyyy,  \
                             tr_z_zzz_xyyz,  \
                             tr_z_zzz_xyz,   \
                             tr_z_zzz_xyzz,  \
                             tr_z_zzz_xzz,   \
                             tr_z_zzz_xzzz,  \
                             tr_z_zzz_yyy,   \
                             tr_z_zzz_yyyy,  \
                             tr_z_zzz_yyyz,  \
                             tr_z_zzz_yyz,   \
                             tr_z_zzz_yyzz,  \
                             tr_z_zzz_yzz,   \
                             tr_z_zzz_yzzz,  \
                             tr_z_zzz_zzz,   \
                             tr_z_zzz_zzzz,  \
                             tr_z_zzzz_xxxx, \
                             tr_z_zzzz_xxxy, \
                             tr_z_zzzz_xxxz, \
                             tr_z_zzzz_xxyy, \
                             tr_z_zzzz_xxyz, \
                             tr_z_zzzz_xxzz, \
                             tr_z_zzzz_xyyy, \
                             tr_z_zzzz_xyyz, \
                             tr_z_zzzz_xyzz, \
                             tr_z_zzzz_xzzz, \
                             tr_z_zzzz_yyyy, \
                             tr_z_zzzz_yyyz, \
                             tr_z_zzzz_yyzz, \
                             tr_z_zzzz_yzzz, \
                             tr_z_zzzz_zzzz, \
                             ts_zzz_xxxx,    \
                             ts_zzz_xxxy,    \
                             ts_zzz_xxxz,    \
                             ts_zzz_xxyy,    \
                             ts_zzz_xxyz,    \
                             ts_zzz_xxzz,    \
                             ts_zzz_xyyy,    \
                             ts_zzz_xyyz,    \
                             ts_zzz_xyzz,    \
                             ts_zzz_xzzz,    \
                             ts_zzz_yyyy,    \
                             ts_zzz_yyyz,    \
                             ts_zzz_yyzz,    \
                             ts_zzz_yzzz,    \
                             ts_zzz_zzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzz_xxxx[i] = 3.0 * tr_z_zz_xxxx[i] * fe_0 + ts_zzz_xxxx[i] * fe_0 + tr_z_zzz_xxxx[i] * pa_z[i];

        tr_z_zzzz_xxxy[i] = 3.0 * tr_z_zz_xxxy[i] * fe_0 + ts_zzz_xxxy[i] * fe_0 + tr_z_zzz_xxxy[i] * pa_z[i];

        tr_z_zzzz_xxxz[i] = 3.0 * tr_z_zz_xxxz[i] * fe_0 + tr_z_zzz_xxx[i] * fe_0 + ts_zzz_xxxz[i] * fe_0 + tr_z_zzz_xxxz[i] * pa_z[i];

        tr_z_zzzz_xxyy[i] = 3.0 * tr_z_zz_xxyy[i] * fe_0 + ts_zzz_xxyy[i] * fe_0 + tr_z_zzz_xxyy[i] * pa_z[i];

        tr_z_zzzz_xxyz[i] = 3.0 * tr_z_zz_xxyz[i] * fe_0 + tr_z_zzz_xxy[i] * fe_0 + ts_zzz_xxyz[i] * fe_0 + tr_z_zzz_xxyz[i] * pa_z[i];

        tr_z_zzzz_xxzz[i] = 3.0 * tr_z_zz_xxzz[i] * fe_0 + 2.0 * tr_z_zzz_xxz[i] * fe_0 + ts_zzz_xxzz[i] * fe_0 + tr_z_zzz_xxzz[i] * pa_z[i];

        tr_z_zzzz_xyyy[i] = 3.0 * tr_z_zz_xyyy[i] * fe_0 + ts_zzz_xyyy[i] * fe_0 + tr_z_zzz_xyyy[i] * pa_z[i];

        tr_z_zzzz_xyyz[i] = 3.0 * tr_z_zz_xyyz[i] * fe_0 + tr_z_zzz_xyy[i] * fe_0 + ts_zzz_xyyz[i] * fe_0 + tr_z_zzz_xyyz[i] * pa_z[i];

        tr_z_zzzz_xyzz[i] = 3.0 * tr_z_zz_xyzz[i] * fe_0 + 2.0 * tr_z_zzz_xyz[i] * fe_0 + ts_zzz_xyzz[i] * fe_0 + tr_z_zzz_xyzz[i] * pa_z[i];

        tr_z_zzzz_xzzz[i] = 3.0 * tr_z_zz_xzzz[i] * fe_0 + 3.0 * tr_z_zzz_xzz[i] * fe_0 + ts_zzz_xzzz[i] * fe_0 + tr_z_zzz_xzzz[i] * pa_z[i];

        tr_z_zzzz_yyyy[i] = 3.0 * tr_z_zz_yyyy[i] * fe_0 + ts_zzz_yyyy[i] * fe_0 + tr_z_zzz_yyyy[i] * pa_z[i];

        tr_z_zzzz_yyyz[i] = 3.0 * tr_z_zz_yyyz[i] * fe_0 + tr_z_zzz_yyy[i] * fe_0 + ts_zzz_yyyz[i] * fe_0 + tr_z_zzz_yyyz[i] * pa_z[i];

        tr_z_zzzz_yyzz[i] = 3.0 * tr_z_zz_yyzz[i] * fe_0 + 2.0 * tr_z_zzz_yyz[i] * fe_0 + ts_zzz_yyzz[i] * fe_0 + tr_z_zzz_yyzz[i] * pa_z[i];

        tr_z_zzzz_yzzz[i] = 3.0 * tr_z_zz_yzzz[i] * fe_0 + 3.0 * tr_z_zzz_yzz[i] * fe_0 + ts_zzz_yzzz[i] * fe_0 + tr_z_zzz_yzzz[i] * pa_z[i];

        tr_z_zzzz_zzzz[i] = 3.0 * tr_z_zz_zzzz[i] * fe_0 + 4.0 * tr_z_zzz_zzz[i] * fe_0 + ts_zzz_zzzz[i] * fe_0 + tr_z_zzz_zzzz[i] * pa_z[i];
    }
}

}  // namespace diprec
