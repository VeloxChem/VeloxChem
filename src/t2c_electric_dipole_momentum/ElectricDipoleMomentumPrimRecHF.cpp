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

#include "ElectricDipoleMomentumPrimRecHF.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_hf(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_hf,
                                      const size_t              idx_dip_ff,
                                      const size_t              idx_dip_gd,
                                      const size_t              idx_ovl_gf,
                                      const size_t              idx_dip_gf,
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

    auto tr_x_xxy_zzz = pbuffer.data(idx_dip_ff + 19);

    auto tr_x_xxz_xxx = pbuffer.data(idx_dip_ff + 20);

    auto tr_x_xxz_xxy = pbuffer.data(idx_dip_ff + 21);

    auto tr_x_xxz_xxz = pbuffer.data(idx_dip_ff + 22);

    auto tr_x_xxz_xyy = pbuffer.data(idx_dip_ff + 23);

    auto tr_x_xxz_xyz = pbuffer.data(idx_dip_ff + 24);

    auto tr_x_xxz_xzz = pbuffer.data(idx_dip_ff + 25);

    auto tr_x_xxz_yyy = pbuffer.data(idx_dip_ff + 26);

    auto tr_x_xxz_zzz = pbuffer.data(idx_dip_ff + 29);

    auto tr_x_xyy_xxx = pbuffer.data(idx_dip_ff + 30);

    auto tr_x_xyy_xxy = pbuffer.data(idx_dip_ff + 31);

    auto tr_x_xyy_xxz = pbuffer.data(idx_dip_ff + 32);

    auto tr_x_xyy_xyy = pbuffer.data(idx_dip_ff + 33);

    auto tr_x_xyy_xzz = pbuffer.data(idx_dip_ff + 35);

    auto tr_x_xyy_yyy = pbuffer.data(idx_dip_ff + 36);

    auto tr_x_xyy_yyz = pbuffer.data(idx_dip_ff + 37);

    auto tr_x_xyy_yzz = pbuffer.data(idx_dip_ff + 38);

    auto tr_x_xyz_xxz = pbuffer.data(idx_dip_ff + 42);

    auto tr_x_xyz_xzz = pbuffer.data(idx_dip_ff + 45);

    auto tr_x_xzz_xxx = pbuffer.data(idx_dip_ff + 50);

    auto tr_x_xzz_xxy = pbuffer.data(idx_dip_ff + 51);

    auto tr_x_xzz_xxz = pbuffer.data(idx_dip_ff + 52);

    auto tr_x_xzz_xyy = pbuffer.data(idx_dip_ff + 53);

    auto tr_x_xzz_xzz = pbuffer.data(idx_dip_ff + 55);

    auto tr_x_xzz_yyz = pbuffer.data(idx_dip_ff + 57);

    auto tr_x_xzz_yzz = pbuffer.data(idx_dip_ff + 58);

    auto tr_x_xzz_zzz = pbuffer.data(idx_dip_ff + 59);

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

    auto tr_x_yyz_xxy = pbuffer.data(idx_dip_ff + 71);

    auto tr_x_yyz_xxz = pbuffer.data(idx_dip_ff + 72);

    auto tr_x_yyz_xyy = pbuffer.data(idx_dip_ff + 73);

    auto tr_x_yyz_xzz = pbuffer.data(idx_dip_ff + 75);

    auto tr_x_yyz_yyy = pbuffer.data(idx_dip_ff + 76);

    auto tr_x_yyz_zzz = pbuffer.data(idx_dip_ff + 79);

    auto tr_x_yzz_xxx = pbuffer.data(idx_dip_ff + 80);

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

    auto tr_y_xxy_zzz = pbuffer.data(idx_dip_ff + 119);

    auto tr_y_xxz_xxx = pbuffer.data(idx_dip_ff + 120);

    auto tr_y_xxz_xxy = pbuffer.data(idx_dip_ff + 121);

    auto tr_y_xxz_xyy = pbuffer.data(idx_dip_ff + 123);

    auto tr_y_xxz_yyz = pbuffer.data(idx_dip_ff + 127);

    auto tr_y_xxz_yzz = pbuffer.data(idx_dip_ff + 128);

    auto tr_y_xxz_zzz = pbuffer.data(idx_dip_ff + 129);

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

    auto tr_y_xyz_yyz = pbuffer.data(idx_dip_ff + 147);

    auto tr_y_xyz_yzz = pbuffer.data(idx_dip_ff + 148);

    auto tr_y_xyz_zzz = pbuffer.data(idx_dip_ff + 149);

    auto tr_y_xzz_xxz = pbuffer.data(idx_dip_ff + 152);

    auto tr_y_xzz_xyz = pbuffer.data(idx_dip_ff + 154);

    auto tr_y_xzz_xzz = pbuffer.data(idx_dip_ff + 155);

    auto tr_y_xzz_yyy = pbuffer.data(idx_dip_ff + 156);

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

    auto tr_y_yyz_xxx = pbuffer.data(idx_dip_ff + 170);

    auto tr_y_yyz_xxy = pbuffer.data(idx_dip_ff + 171);

    auto tr_y_yyz_xyy = pbuffer.data(idx_dip_ff + 173);

    auto tr_y_yyz_xyz = pbuffer.data(idx_dip_ff + 174);

    auto tr_y_yyz_yyy = pbuffer.data(idx_dip_ff + 176);

    auto tr_y_yyz_yyz = pbuffer.data(idx_dip_ff + 177);

    auto tr_y_yyz_yzz = pbuffer.data(idx_dip_ff + 178);

    auto tr_y_yyz_zzz = pbuffer.data(idx_dip_ff + 179);

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

    auto tr_z_xxy_xxx = pbuffer.data(idx_dip_ff + 210);

    auto tr_z_xxy_xxz = pbuffer.data(idx_dip_ff + 212);

    auto tr_z_xxy_xzz = pbuffer.data(idx_dip_ff + 215);

    auto tr_z_xxy_yyy = pbuffer.data(idx_dip_ff + 216);

    auto tr_z_xxy_yyz = pbuffer.data(idx_dip_ff + 217);

    auto tr_z_xxy_yzz = pbuffer.data(idx_dip_ff + 218);

    auto tr_z_xxz_xxx = pbuffer.data(idx_dip_ff + 220);

    auto tr_z_xxz_xxz = pbuffer.data(idx_dip_ff + 222);

    auto tr_z_xxz_xyz = pbuffer.data(idx_dip_ff + 224);

    auto tr_z_xxz_xzz = pbuffer.data(idx_dip_ff + 225);

    auto tr_z_xxz_yyy = pbuffer.data(idx_dip_ff + 226);

    auto tr_z_xxz_yyz = pbuffer.data(idx_dip_ff + 227);

    auto tr_z_xxz_yzz = pbuffer.data(idx_dip_ff + 228);

    auto tr_z_xxz_zzz = pbuffer.data(idx_dip_ff + 229);

    auto tr_z_xyy_xxy = pbuffer.data(idx_dip_ff + 231);

    auto tr_z_xyy_xyy = pbuffer.data(idx_dip_ff + 233);

    auto tr_z_xyy_xyz = pbuffer.data(idx_dip_ff + 234);

    auto tr_z_xyy_yyy = pbuffer.data(idx_dip_ff + 236);

    auto tr_z_xyy_yyz = pbuffer.data(idx_dip_ff + 237);

    auto tr_z_xyy_yzz = pbuffer.data(idx_dip_ff + 238);

    auto tr_z_xyy_zzz = pbuffer.data(idx_dip_ff + 239);

    auto tr_z_xyz_yyy = pbuffer.data(idx_dip_ff + 246);

    auto tr_z_xyz_yyz = pbuffer.data(idx_dip_ff + 247);

    auto tr_z_xyz_yzz = pbuffer.data(idx_dip_ff + 248);

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

    auto tr_z_yyz_xxz = pbuffer.data(idx_dip_ff + 272);

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

    // Set up components of auxiliary buffer : GD

    auto tr_x_xxxx_xx = pbuffer.data(idx_dip_gd);

    auto tr_x_xxxx_xy = pbuffer.data(idx_dip_gd + 1);

    auto tr_x_xxxx_xz = pbuffer.data(idx_dip_gd + 2);

    auto tr_x_xxxx_yy = pbuffer.data(idx_dip_gd + 3);

    auto tr_x_xxxx_yz = pbuffer.data(idx_dip_gd + 4);

    auto tr_x_xxxx_zz = pbuffer.data(idx_dip_gd + 5);

    auto tr_x_xxxy_xx = pbuffer.data(idx_dip_gd + 6);

    auto tr_x_xxxy_xy = pbuffer.data(idx_dip_gd + 7);

    auto tr_x_xxxy_xz = pbuffer.data(idx_dip_gd + 8);

    auto tr_x_xxxz_xx = pbuffer.data(idx_dip_gd + 12);

    auto tr_x_xxxz_xy = pbuffer.data(idx_dip_gd + 13);

    auto tr_x_xxxz_xz = pbuffer.data(idx_dip_gd + 14);

    auto tr_x_xxxz_yz = pbuffer.data(idx_dip_gd + 16);

    auto tr_x_xxxz_zz = pbuffer.data(idx_dip_gd + 17);

    auto tr_x_xxyy_xx = pbuffer.data(idx_dip_gd + 18);

    auto tr_x_xxyy_xy = pbuffer.data(idx_dip_gd + 19);

    auto tr_x_xxyy_xz = pbuffer.data(idx_dip_gd + 20);

    auto tr_x_xxyy_yy = pbuffer.data(idx_dip_gd + 21);

    auto tr_x_xxyy_yz = pbuffer.data(idx_dip_gd + 22);

    auto tr_x_xxzz_xx = pbuffer.data(idx_dip_gd + 30);

    auto tr_x_xxzz_xy = pbuffer.data(idx_dip_gd + 31);

    auto tr_x_xxzz_xz = pbuffer.data(idx_dip_gd + 32);

    auto tr_x_xxzz_yy = pbuffer.data(idx_dip_gd + 33);

    auto tr_x_xxzz_yz = pbuffer.data(idx_dip_gd + 34);

    auto tr_x_xxzz_zz = pbuffer.data(idx_dip_gd + 35);

    auto tr_x_xyyy_xy = pbuffer.data(idx_dip_gd + 37);

    auto tr_x_xzzz_xx = pbuffer.data(idx_dip_gd + 54);

    auto tr_x_xzzz_xy = pbuffer.data(idx_dip_gd + 55);

    auto tr_x_xzzz_xz = pbuffer.data(idx_dip_gd + 56);

    auto tr_x_yyyy_xx = pbuffer.data(idx_dip_gd + 60);

    auto tr_x_yyyy_xy = pbuffer.data(idx_dip_gd + 61);

    auto tr_x_yyyy_xz = pbuffer.data(idx_dip_gd + 62);

    auto tr_x_yyyy_yy = pbuffer.data(idx_dip_gd + 63);

    auto tr_x_yyyy_yz = pbuffer.data(idx_dip_gd + 64);

    auto tr_x_yyyy_zz = pbuffer.data(idx_dip_gd + 65);

    auto tr_x_yyzz_xz = pbuffer.data(idx_dip_gd + 74);

    auto tr_x_yyzz_yz = pbuffer.data(idx_dip_gd + 76);

    auto tr_x_yyzz_zz = pbuffer.data(idx_dip_gd + 77);

    auto tr_x_yzzz_xz = pbuffer.data(idx_dip_gd + 80);

    auto tr_x_yzzz_yz = pbuffer.data(idx_dip_gd + 82);

    auto tr_x_yzzz_zz = pbuffer.data(idx_dip_gd + 83);

    auto tr_x_zzzz_xx = pbuffer.data(idx_dip_gd + 84);

    auto tr_x_zzzz_xy = pbuffer.data(idx_dip_gd + 85);

    auto tr_x_zzzz_xz = pbuffer.data(idx_dip_gd + 86);

    auto tr_x_zzzz_yy = pbuffer.data(idx_dip_gd + 87);

    auto tr_x_zzzz_yz = pbuffer.data(idx_dip_gd + 88);

    auto tr_x_zzzz_zz = pbuffer.data(idx_dip_gd + 89);

    auto tr_y_xxxx_xx = pbuffer.data(idx_dip_gd + 90);

    auto tr_y_xxxx_xy = pbuffer.data(idx_dip_gd + 91);

    auto tr_y_xxxx_xz = pbuffer.data(idx_dip_gd + 92);

    auto tr_y_xxxx_yy = pbuffer.data(idx_dip_gd + 93);

    auto tr_y_xxxx_yz = pbuffer.data(idx_dip_gd + 94);

    auto tr_y_xxxx_zz = pbuffer.data(idx_dip_gd + 95);

    auto tr_y_xxxy_xy = pbuffer.data(idx_dip_gd + 97);

    auto tr_y_xxxy_yy = pbuffer.data(idx_dip_gd + 99);

    auto tr_y_xxxy_yz = pbuffer.data(idx_dip_gd + 100);

    auto tr_y_xxyy_xx = pbuffer.data(idx_dip_gd + 108);

    auto tr_y_xxyy_xy = pbuffer.data(idx_dip_gd + 109);

    auto tr_y_xxyy_xz = pbuffer.data(idx_dip_gd + 110);

    auto tr_y_xxyy_yy = pbuffer.data(idx_dip_gd + 111);

    auto tr_y_xxyy_yz = pbuffer.data(idx_dip_gd + 112);

    auto tr_y_xxyy_zz = pbuffer.data(idx_dip_gd + 113);

    auto tr_y_xxzz_xz = pbuffer.data(idx_dip_gd + 122);

    auto tr_y_xxzz_yz = pbuffer.data(idx_dip_gd + 124);

    auto tr_y_xxzz_zz = pbuffer.data(idx_dip_gd + 125);

    auto tr_y_xyyy_xx = pbuffer.data(idx_dip_gd + 126);

    auto tr_y_xyyy_xy = pbuffer.data(idx_dip_gd + 127);

    auto tr_y_xyyy_xz = pbuffer.data(idx_dip_gd + 128);

    auto tr_y_xyyy_yy = pbuffer.data(idx_dip_gd + 129);

    auto tr_y_xyyy_yz = pbuffer.data(idx_dip_gd + 130);

    auto tr_y_xyyy_zz = pbuffer.data(idx_dip_gd + 131);

    auto tr_y_xyzz_yz = pbuffer.data(idx_dip_gd + 142);

    auto tr_y_xzzz_xz = pbuffer.data(idx_dip_gd + 146);

    auto tr_y_xzzz_yz = pbuffer.data(idx_dip_gd + 148);

    auto tr_y_xzzz_zz = pbuffer.data(idx_dip_gd + 149);

    auto tr_y_yyyy_xx = pbuffer.data(idx_dip_gd + 150);

    auto tr_y_yyyy_xy = pbuffer.data(idx_dip_gd + 151);

    auto tr_y_yyyy_xz = pbuffer.data(idx_dip_gd + 152);

    auto tr_y_yyyy_yy = pbuffer.data(idx_dip_gd + 153);

    auto tr_y_yyyy_yz = pbuffer.data(idx_dip_gd + 154);

    auto tr_y_yyyy_zz = pbuffer.data(idx_dip_gd + 155);

    auto tr_y_yyyz_xy = pbuffer.data(idx_dip_gd + 157);

    auto tr_y_yyyz_xz = pbuffer.data(idx_dip_gd + 158);

    auto tr_y_yyyz_yy = pbuffer.data(idx_dip_gd + 159);

    auto tr_y_yyyz_yz = pbuffer.data(idx_dip_gd + 160);

    auto tr_y_yyyz_zz = pbuffer.data(idx_dip_gd + 161);

    auto tr_y_yyzz_xx = pbuffer.data(idx_dip_gd + 162);

    auto tr_y_yyzz_xy = pbuffer.data(idx_dip_gd + 163);

    auto tr_y_yyzz_xz = pbuffer.data(idx_dip_gd + 164);

    auto tr_y_yyzz_yy = pbuffer.data(idx_dip_gd + 165);

    auto tr_y_yyzz_yz = pbuffer.data(idx_dip_gd + 166);

    auto tr_y_yyzz_zz = pbuffer.data(idx_dip_gd + 167);

    auto tr_y_yzzz_xx = pbuffer.data(idx_dip_gd + 168);

    auto tr_y_yzzz_xy = pbuffer.data(idx_dip_gd + 169);

    auto tr_y_yzzz_xz = pbuffer.data(idx_dip_gd + 170);

    auto tr_y_yzzz_yy = pbuffer.data(idx_dip_gd + 171);

    auto tr_y_yzzz_yz = pbuffer.data(idx_dip_gd + 172);

    auto tr_y_yzzz_zz = pbuffer.data(idx_dip_gd + 173);

    auto tr_y_zzzz_xx = pbuffer.data(idx_dip_gd + 174);

    auto tr_y_zzzz_xy = pbuffer.data(idx_dip_gd + 175);

    auto tr_y_zzzz_xz = pbuffer.data(idx_dip_gd + 176);

    auto tr_y_zzzz_yy = pbuffer.data(idx_dip_gd + 177);

    auto tr_y_zzzz_yz = pbuffer.data(idx_dip_gd + 178);

    auto tr_y_zzzz_zz = pbuffer.data(idx_dip_gd + 179);

    auto tr_z_xxxx_xx = pbuffer.data(idx_dip_gd + 180);

    auto tr_z_xxxx_xy = pbuffer.data(idx_dip_gd + 181);

    auto tr_z_xxxx_xz = pbuffer.data(idx_dip_gd + 182);

    auto tr_z_xxxx_yy = pbuffer.data(idx_dip_gd + 183);

    auto tr_z_xxxx_yz = pbuffer.data(idx_dip_gd + 184);

    auto tr_z_xxxx_zz = pbuffer.data(idx_dip_gd + 185);

    auto tr_z_xxxz_xx = pbuffer.data(idx_dip_gd + 192);

    auto tr_z_xxxz_xy = pbuffer.data(idx_dip_gd + 193);

    auto tr_z_xxxz_xz = pbuffer.data(idx_dip_gd + 194);

    auto tr_z_xxxz_yz = pbuffer.data(idx_dip_gd + 196);

    auto tr_z_xxxz_zz = pbuffer.data(idx_dip_gd + 197);

    auto tr_z_xxyy_xy = pbuffer.data(idx_dip_gd + 199);

    auto tr_z_xxyy_yy = pbuffer.data(idx_dip_gd + 201);

    auto tr_z_xxyy_yz = pbuffer.data(idx_dip_gd + 202);

    auto tr_z_xxzz_xx = pbuffer.data(idx_dip_gd + 210);

    auto tr_z_xxzz_xy = pbuffer.data(idx_dip_gd + 211);

    auto tr_z_xxzz_xz = pbuffer.data(idx_dip_gd + 212);

    auto tr_z_xxzz_yy = pbuffer.data(idx_dip_gd + 213);

    auto tr_z_xxzz_yz = pbuffer.data(idx_dip_gd + 214);

    auto tr_z_xxzz_zz = pbuffer.data(idx_dip_gd + 215);

    auto tr_z_xyyy_xy = pbuffer.data(idx_dip_gd + 217);

    auto tr_z_xyyy_yy = pbuffer.data(idx_dip_gd + 219);

    auto tr_z_xyyy_yz = pbuffer.data(idx_dip_gd + 220);

    auto tr_z_xyyz_yz = pbuffer.data(idx_dip_gd + 226);

    auto tr_z_xzzz_xx = pbuffer.data(idx_dip_gd + 234);

    auto tr_z_xzzz_xy = pbuffer.data(idx_dip_gd + 235);

    auto tr_z_xzzz_xz = pbuffer.data(idx_dip_gd + 236);

    auto tr_z_xzzz_yy = pbuffer.data(idx_dip_gd + 237);

    auto tr_z_xzzz_yz = pbuffer.data(idx_dip_gd + 238);

    auto tr_z_xzzz_zz = pbuffer.data(idx_dip_gd + 239);

    auto tr_z_yyyy_xx = pbuffer.data(idx_dip_gd + 240);

    auto tr_z_yyyy_xy = pbuffer.data(idx_dip_gd + 241);

    auto tr_z_yyyy_xz = pbuffer.data(idx_dip_gd + 242);

    auto tr_z_yyyy_yy = pbuffer.data(idx_dip_gd + 243);

    auto tr_z_yyyy_yz = pbuffer.data(idx_dip_gd + 244);

    auto tr_z_yyyy_zz = pbuffer.data(idx_dip_gd + 245);

    auto tr_z_yyyz_xx = pbuffer.data(idx_dip_gd + 246);

    auto tr_z_yyyz_xy = pbuffer.data(idx_dip_gd + 247);

    auto tr_z_yyyz_xz = pbuffer.data(idx_dip_gd + 248);

    auto tr_z_yyyz_yy = pbuffer.data(idx_dip_gd + 249);

    auto tr_z_yyyz_yz = pbuffer.data(idx_dip_gd + 250);

    auto tr_z_yyyz_zz = pbuffer.data(idx_dip_gd + 251);

    auto tr_z_yyzz_xx = pbuffer.data(idx_dip_gd + 252);

    auto tr_z_yyzz_xy = pbuffer.data(idx_dip_gd + 253);

    auto tr_z_yyzz_xz = pbuffer.data(idx_dip_gd + 254);

    auto tr_z_yyzz_yy = pbuffer.data(idx_dip_gd + 255);

    auto tr_z_yyzz_yz = pbuffer.data(idx_dip_gd + 256);

    auto tr_z_yyzz_zz = pbuffer.data(idx_dip_gd + 257);

    auto tr_z_yzzz_xx = pbuffer.data(idx_dip_gd + 258);

    auto tr_z_yzzz_xy = pbuffer.data(idx_dip_gd + 259);

    auto tr_z_yzzz_xz = pbuffer.data(idx_dip_gd + 260);

    auto tr_z_yzzz_yy = pbuffer.data(idx_dip_gd + 261);

    auto tr_z_yzzz_yz = pbuffer.data(idx_dip_gd + 262);

    auto tr_z_yzzz_zz = pbuffer.data(idx_dip_gd + 263);

    auto tr_z_zzzz_xx = pbuffer.data(idx_dip_gd + 264);

    auto tr_z_zzzz_xy = pbuffer.data(idx_dip_gd + 265);

    auto tr_z_zzzz_xz = pbuffer.data(idx_dip_gd + 266);

    auto tr_z_zzzz_yy = pbuffer.data(idx_dip_gd + 267);

    auto tr_z_zzzz_yz = pbuffer.data(idx_dip_gd + 268);

    auto tr_z_zzzz_zz = pbuffer.data(idx_dip_gd + 269);

    // Set up components of auxiliary buffer : GF

    auto ts_xxxx_xxx = pbuffer.data(idx_ovl_gf);

    auto ts_xxxx_xxy = pbuffer.data(idx_ovl_gf + 1);

    auto ts_xxxx_xxz = pbuffer.data(idx_ovl_gf + 2);

    auto ts_xxxx_xyy = pbuffer.data(idx_ovl_gf + 3);

    auto ts_xxxx_xyz = pbuffer.data(idx_ovl_gf + 4);

    auto ts_xxxx_xzz = pbuffer.data(idx_ovl_gf + 5);

    auto ts_xxxx_yyy = pbuffer.data(idx_ovl_gf + 6);

    auto ts_xxxx_yyz = pbuffer.data(idx_ovl_gf + 7);

    auto ts_xxxx_yzz = pbuffer.data(idx_ovl_gf + 8);

    auto ts_xxxx_zzz = pbuffer.data(idx_ovl_gf + 9);

    auto ts_xxxz_xxz = pbuffer.data(idx_ovl_gf + 22);

    auto ts_xxxz_xzz = pbuffer.data(idx_ovl_gf + 25);

    auto ts_xxyy_xxy = pbuffer.data(idx_ovl_gf + 31);

    auto ts_xxyy_xyy = pbuffer.data(idx_ovl_gf + 33);

    auto ts_xxyy_yyy = pbuffer.data(idx_ovl_gf + 36);

    auto ts_xxyy_yyz = pbuffer.data(idx_ovl_gf + 37);

    auto ts_xxyy_yzz = pbuffer.data(idx_ovl_gf + 38);

    auto ts_xxzz_xxx = pbuffer.data(idx_ovl_gf + 50);

    auto ts_xxzz_xxz = pbuffer.data(idx_ovl_gf + 52);

    auto ts_xxzz_xzz = pbuffer.data(idx_ovl_gf + 55);

    auto ts_xxzz_yyz = pbuffer.data(idx_ovl_gf + 57);

    auto ts_xxzz_yzz = pbuffer.data(idx_ovl_gf + 58);

    auto ts_xxzz_zzz = pbuffer.data(idx_ovl_gf + 59);

    auto ts_xyyy_yyy = pbuffer.data(idx_ovl_gf + 66);

    auto ts_xyyy_yyz = pbuffer.data(idx_ovl_gf + 67);

    auto ts_xyyy_yzz = pbuffer.data(idx_ovl_gf + 68);

    auto ts_xzzz_yyz = pbuffer.data(idx_ovl_gf + 97);

    auto ts_xzzz_yzz = pbuffer.data(idx_ovl_gf + 98);

    auto ts_xzzz_zzz = pbuffer.data(idx_ovl_gf + 99);

    auto ts_yyyy_xxx = pbuffer.data(idx_ovl_gf + 100);

    auto ts_yyyy_xxy = pbuffer.data(idx_ovl_gf + 101);

    auto ts_yyyy_xxz = pbuffer.data(idx_ovl_gf + 102);

    auto ts_yyyy_xyy = pbuffer.data(idx_ovl_gf + 103);

    auto ts_yyyy_xyz = pbuffer.data(idx_ovl_gf + 104);

    auto ts_yyyy_xzz = pbuffer.data(idx_ovl_gf + 105);

    auto ts_yyyy_yyy = pbuffer.data(idx_ovl_gf + 106);

    auto ts_yyyy_yyz = pbuffer.data(idx_ovl_gf + 107);

    auto ts_yyyy_yzz = pbuffer.data(idx_ovl_gf + 108);

    auto ts_yyyy_zzz = pbuffer.data(idx_ovl_gf + 109);

    auto ts_yyyz_yyz = pbuffer.data(idx_ovl_gf + 117);

    auto ts_yyyz_yzz = pbuffer.data(idx_ovl_gf + 118);

    auto ts_yyyz_zzz = pbuffer.data(idx_ovl_gf + 119);

    auto ts_yyzz_xxz = pbuffer.data(idx_ovl_gf + 122);

    auto ts_yyzz_xyz = pbuffer.data(idx_ovl_gf + 124);

    auto ts_yyzz_xzz = pbuffer.data(idx_ovl_gf + 125);

    auto ts_yyzz_yyy = pbuffer.data(idx_ovl_gf + 126);

    auto ts_yyzz_yyz = pbuffer.data(idx_ovl_gf + 127);

    auto ts_yyzz_yzz = pbuffer.data(idx_ovl_gf + 128);

    auto ts_yyzz_zzz = pbuffer.data(idx_ovl_gf + 129);

    auto ts_yzzz_xxz = pbuffer.data(idx_ovl_gf + 132);

    auto ts_yzzz_xzz = pbuffer.data(idx_ovl_gf + 135);

    auto ts_yzzz_yyy = pbuffer.data(idx_ovl_gf + 136);

    auto ts_yzzz_yyz = pbuffer.data(idx_ovl_gf + 137);

    auto ts_yzzz_yzz = pbuffer.data(idx_ovl_gf + 138);

    auto ts_yzzz_zzz = pbuffer.data(idx_ovl_gf + 139);

    auto ts_zzzz_xxx = pbuffer.data(idx_ovl_gf + 140);

    auto ts_zzzz_xxy = pbuffer.data(idx_ovl_gf + 141);

    auto ts_zzzz_xxz = pbuffer.data(idx_ovl_gf + 142);

    auto ts_zzzz_xyy = pbuffer.data(idx_ovl_gf + 143);

    auto ts_zzzz_xyz = pbuffer.data(idx_ovl_gf + 144);

    auto ts_zzzz_xzz = pbuffer.data(idx_ovl_gf + 145);

    auto ts_zzzz_yyy = pbuffer.data(idx_ovl_gf + 146);

    auto ts_zzzz_yyz = pbuffer.data(idx_ovl_gf + 147);

    auto ts_zzzz_yzz = pbuffer.data(idx_ovl_gf + 148);

    auto ts_zzzz_zzz = pbuffer.data(idx_ovl_gf + 149);

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

    auto tr_x_xxxy_yyy = pbuffer.data(idx_dip_gf + 16);

    auto tr_x_xxxy_zzz = pbuffer.data(idx_dip_gf + 19);

    auto tr_x_xxxz_xxx = pbuffer.data(idx_dip_gf + 20);

    auto tr_x_xxxz_xxy = pbuffer.data(idx_dip_gf + 21);

    auto tr_x_xxxz_xxz = pbuffer.data(idx_dip_gf + 22);

    auto tr_x_xxxz_xyy = pbuffer.data(idx_dip_gf + 23);

    auto tr_x_xxxz_xyz = pbuffer.data(idx_dip_gf + 24);

    auto tr_x_xxxz_xzz = pbuffer.data(idx_dip_gf + 25);

    auto tr_x_xxxz_yyy = pbuffer.data(idx_dip_gf + 26);

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

    auto tr_x_xxyy_zzz = pbuffer.data(idx_dip_gf + 39);

    auto tr_x_xxyz_xxz = pbuffer.data(idx_dip_gf + 42);

    auto tr_x_xxyz_xzz = pbuffer.data(idx_dip_gf + 45);

    auto tr_x_xxyz_zzz = pbuffer.data(idx_dip_gf + 49);

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

    auto tr_x_xyyy_xxx = pbuffer.data(idx_dip_gf + 60);

    auto tr_x_xyyy_xxy = pbuffer.data(idx_dip_gf + 61);

    auto tr_x_xyyy_xxz = pbuffer.data(idx_dip_gf + 62);

    auto tr_x_xyyy_xyy = pbuffer.data(idx_dip_gf + 63);

    auto tr_x_xyyy_xyz = pbuffer.data(idx_dip_gf + 64);

    auto tr_x_xyyy_xzz = pbuffer.data(idx_dip_gf + 65);

    auto tr_x_xyyy_yyy = pbuffer.data(idx_dip_gf + 66);

    auto tr_x_xyyy_yyz = pbuffer.data(idx_dip_gf + 67);

    auto tr_x_xyyy_yzz = pbuffer.data(idx_dip_gf + 68);

    auto tr_x_xyyz_xxy = pbuffer.data(idx_dip_gf + 71);

    auto tr_x_xyyz_xxz = pbuffer.data(idx_dip_gf + 72);

    auto tr_x_xyyz_xyy = pbuffer.data(idx_dip_gf + 73);

    auto tr_x_xyyz_xzz = pbuffer.data(idx_dip_gf + 75);

    auto tr_x_xyzz_xxx = pbuffer.data(idx_dip_gf + 80);

    auto tr_x_xyzz_xxz = pbuffer.data(idx_dip_gf + 82);

    auto tr_x_xyzz_xzz = pbuffer.data(idx_dip_gf + 85);

    auto tr_x_xzzz_xxx = pbuffer.data(idx_dip_gf + 90);

    auto tr_x_xzzz_xxy = pbuffer.data(idx_dip_gf + 91);

    auto tr_x_xzzz_xxz = pbuffer.data(idx_dip_gf + 92);

    auto tr_x_xzzz_xyy = pbuffer.data(idx_dip_gf + 93);

    auto tr_x_xzzz_xyz = pbuffer.data(idx_dip_gf + 94);

    auto tr_x_xzzz_xzz = pbuffer.data(idx_dip_gf + 95);

    auto tr_x_xzzz_yyz = pbuffer.data(idx_dip_gf + 97);

    auto tr_x_xzzz_yzz = pbuffer.data(idx_dip_gf + 98);

    auto tr_x_xzzz_zzz = pbuffer.data(idx_dip_gf + 99);

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

    auto tr_x_yyyz_xxy = pbuffer.data(idx_dip_gf + 111);

    auto tr_x_yyyz_xxz = pbuffer.data(idx_dip_gf + 112);

    auto tr_x_yyyz_xyy = pbuffer.data(idx_dip_gf + 113);

    auto tr_x_yyyz_xzz = pbuffer.data(idx_dip_gf + 115);

    auto tr_x_yyyz_yyy = pbuffer.data(idx_dip_gf + 116);

    auto tr_x_yyyz_yyz = pbuffer.data(idx_dip_gf + 117);

    auto tr_x_yyyz_yzz = pbuffer.data(idx_dip_gf + 118);

    auto tr_x_yyyz_zzz = pbuffer.data(idx_dip_gf + 119);

    auto tr_x_yyzz_xxx = pbuffer.data(idx_dip_gf + 120);

    auto tr_x_yyzz_xxy = pbuffer.data(idx_dip_gf + 121);

    auto tr_x_yyzz_xxz = pbuffer.data(idx_dip_gf + 122);

    auto tr_x_yyzz_xyy = pbuffer.data(idx_dip_gf + 123);

    auto tr_x_yyzz_xyz = pbuffer.data(idx_dip_gf + 124);

    auto tr_x_yyzz_xzz = pbuffer.data(idx_dip_gf + 125);

    auto tr_x_yyzz_yyy = pbuffer.data(idx_dip_gf + 126);

    auto tr_x_yyzz_yyz = pbuffer.data(idx_dip_gf + 127);

    auto tr_x_yyzz_yzz = pbuffer.data(idx_dip_gf + 128);

    auto tr_x_yyzz_zzz = pbuffer.data(idx_dip_gf + 129);

    auto tr_x_yzzz_xxx = pbuffer.data(idx_dip_gf + 130);

    auto tr_x_yzzz_xxz = pbuffer.data(idx_dip_gf + 132);

    auto tr_x_yzzz_xyz = pbuffer.data(idx_dip_gf + 134);

    auto tr_x_yzzz_xzz = pbuffer.data(idx_dip_gf + 135);

    auto tr_x_yzzz_yyy = pbuffer.data(idx_dip_gf + 136);

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

    auto tr_y_xxxy_xxx = pbuffer.data(idx_dip_gf + 160);

    auto tr_y_xxxy_xxy = pbuffer.data(idx_dip_gf + 161);

    auto tr_y_xxxy_xyy = pbuffer.data(idx_dip_gf + 163);

    auto tr_y_xxxy_xyz = pbuffer.data(idx_dip_gf + 164);

    auto tr_y_xxxy_yyy = pbuffer.data(idx_dip_gf + 166);

    auto tr_y_xxxy_yyz = pbuffer.data(idx_dip_gf + 167);

    auto tr_y_xxxy_yzz = pbuffer.data(idx_dip_gf + 168);

    auto tr_y_xxxy_zzz = pbuffer.data(idx_dip_gf + 169);

    auto tr_y_xxxz_xxx = pbuffer.data(idx_dip_gf + 170);

    auto tr_y_xxxz_xxy = pbuffer.data(idx_dip_gf + 171);

    auto tr_y_xxxz_xxz = pbuffer.data(idx_dip_gf + 172);

    auto tr_y_xxxz_xyy = pbuffer.data(idx_dip_gf + 173);

    auto tr_y_xxxz_xzz = pbuffer.data(idx_dip_gf + 175);

    auto tr_y_xxxz_yyz = pbuffer.data(idx_dip_gf + 177);

    auto tr_y_xxxz_yzz = pbuffer.data(idx_dip_gf + 178);

    auto tr_y_xxxz_zzz = pbuffer.data(idx_dip_gf + 179);

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

    auto tr_y_xxyz_xxy = pbuffer.data(idx_dip_gf + 191);

    auto tr_y_xxyz_xyy = pbuffer.data(idx_dip_gf + 193);

    auto tr_y_xxyz_yyz = pbuffer.data(idx_dip_gf + 197);

    auto tr_y_xxyz_yzz = pbuffer.data(idx_dip_gf + 198);

    auto tr_y_xxyz_zzz = pbuffer.data(idx_dip_gf + 199);

    auto tr_y_xxzz_xxx = pbuffer.data(idx_dip_gf + 200);

    auto tr_y_xxzz_xxy = pbuffer.data(idx_dip_gf + 201);

    auto tr_y_xxzz_xxz = pbuffer.data(idx_dip_gf + 202);

    auto tr_y_xxzz_xyy = pbuffer.data(idx_dip_gf + 203);

    auto tr_y_xxzz_xyz = pbuffer.data(idx_dip_gf + 204);

    auto tr_y_xxzz_xzz = pbuffer.data(idx_dip_gf + 205);

    auto tr_y_xxzz_yyy = pbuffer.data(idx_dip_gf + 206);

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

    auto tr_y_xyyz_yyz = pbuffer.data(idx_dip_gf + 227);

    auto tr_y_xyyz_yzz = pbuffer.data(idx_dip_gf + 228);

    auto tr_y_xyyz_zzz = pbuffer.data(idx_dip_gf + 229);

    auto tr_y_xyzz_xyz = pbuffer.data(idx_dip_gf + 234);

    auto tr_y_xyzz_yyy = pbuffer.data(idx_dip_gf + 236);

    auto tr_y_xyzz_yyz = pbuffer.data(idx_dip_gf + 237);

    auto tr_y_xyzz_yzz = pbuffer.data(idx_dip_gf + 238);

    auto tr_y_xyzz_zzz = pbuffer.data(idx_dip_gf + 239);

    auto tr_y_xzzz_xxz = pbuffer.data(idx_dip_gf + 242);

    auto tr_y_xzzz_xyz = pbuffer.data(idx_dip_gf + 244);

    auto tr_y_xzzz_xzz = pbuffer.data(idx_dip_gf + 245);

    auto tr_y_xzzz_yyy = pbuffer.data(idx_dip_gf + 246);

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

    auto tr_y_yyyz_xxx = pbuffer.data(idx_dip_gf + 260);

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

    auto tr_z_xxxy_xxx = pbuffer.data(idx_dip_gf + 310);

    auto tr_z_xxxy_xxz = pbuffer.data(idx_dip_gf + 312);

    auto tr_z_xxxy_xzz = pbuffer.data(idx_dip_gf + 315);

    auto tr_z_xxxy_yyy = pbuffer.data(idx_dip_gf + 316);

    auto tr_z_xxxy_yyz = pbuffer.data(idx_dip_gf + 317);

    auto tr_z_xxxy_yzz = pbuffer.data(idx_dip_gf + 318);

    auto tr_z_xxxz_xxx = pbuffer.data(idx_dip_gf + 320);

    auto tr_z_xxxz_xxy = pbuffer.data(idx_dip_gf + 321);

    auto tr_z_xxxz_xxz = pbuffer.data(idx_dip_gf + 322);

    auto tr_z_xxxz_xyy = pbuffer.data(idx_dip_gf + 323);

    auto tr_z_xxxz_xyz = pbuffer.data(idx_dip_gf + 324);

    auto tr_z_xxxz_xzz = pbuffer.data(idx_dip_gf + 325);

    auto tr_z_xxxz_yyy = pbuffer.data(idx_dip_gf + 326);

    auto tr_z_xxxz_yyz = pbuffer.data(idx_dip_gf + 327);

    auto tr_z_xxxz_yzz = pbuffer.data(idx_dip_gf + 328);

    auto tr_z_xxxz_zzz = pbuffer.data(idx_dip_gf + 329);

    auto tr_z_xxyy_xxx = pbuffer.data(idx_dip_gf + 330);

    auto tr_z_xxyy_xxy = pbuffer.data(idx_dip_gf + 331);

    auto tr_z_xxyy_xxz = pbuffer.data(idx_dip_gf + 332);

    auto tr_z_xxyy_xyy = pbuffer.data(idx_dip_gf + 333);

    auto tr_z_xxyy_xyz = pbuffer.data(idx_dip_gf + 334);

    auto tr_z_xxyy_xzz = pbuffer.data(idx_dip_gf + 335);

    auto tr_z_xxyy_yyy = pbuffer.data(idx_dip_gf + 336);

    auto tr_z_xxyy_yyz = pbuffer.data(idx_dip_gf + 337);

    auto tr_z_xxyy_yzz = pbuffer.data(idx_dip_gf + 338);

    auto tr_z_xxyy_zzz = pbuffer.data(idx_dip_gf + 339);

    auto tr_z_xxyz_xxx = pbuffer.data(idx_dip_gf + 340);

    auto tr_z_xxyz_xxz = pbuffer.data(idx_dip_gf + 342);

    auto tr_z_xxyz_xzz = pbuffer.data(idx_dip_gf + 345);

    auto tr_z_xxyz_yyy = pbuffer.data(idx_dip_gf + 346);

    auto tr_z_xxyz_yyz = pbuffer.data(idx_dip_gf + 347);

    auto tr_z_xxyz_yzz = pbuffer.data(idx_dip_gf + 348);

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

    auto tr_z_xyyy_zzz = pbuffer.data(idx_dip_gf + 369);

    auto tr_z_xyyz_xyz = pbuffer.data(idx_dip_gf + 374);

    auto tr_z_xyyz_yyy = pbuffer.data(idx_dip_gf + 376);

    auto tr_z_xyyz_yyz = pbuffer.data(idx_dip_gf + 377);

    auto tr_z_xyyz_yzz = pbuffer.data(idx_dip_gf + 378);

    auto tr_z_xyyz_zzz = pbuffer.data(idx_dip_gf + 379);

    auto tr_z_xyzz_yyy = pbuffer.data(idx_dip_gf + 386);

    auto tr_z_xyzz_yyz = pbuffer.data(idx_dip_gf + 387);

    auto tr_z_xyzz_yzz = pbuffer.data(idx_dip_gf + 388);

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

    // Set up 0-10 components of targeted buffer : HF

    auto tr_x_xxxxx_xxx = pbuffer.data(idx_dip_hf);

    auto tr_x_xxxxx_xxy = pbuffer.data(idx_dip_hf + 1);

    auto tr_x_xxxxx_xxz = pbuffer.data(idx_dip_hf + 2);

    auto tr_x_xxxxx_xyy = pbuffer.data(idx_dip_hf + 3);

    auto tr_x_xxxxx_xyz = pbuffer.data(idx_dip_hf + 4);

    auto tr_x_xxxxx_xzz = pbuffer.data(idx_dip_hf + 5);

    auto tr_x_xxxxx_yyy = pbuffer.data(idx_dip_hf + 6);

    auto tr_x_xxxxx_yyz = pbuffer.data(idx_dip_hf + 7);

    auto tr_x_xxxxx_yzz = pbuffer.data(idx_dip_hf + 8);

    auto tr_x_xxxxx_zzz = pbuffer.data(idx_dip_hf + 9);

#pragma omp simd aligned(pa_x,               \
                             tr_x_xxx_xxx,   \
                             tr_x_xxx_xxy,   \
                             tr_x_xxx_xxz,   \
                             tr_x_xxx_xyy,   \
                             tr_x_xxx_xyz,   \
                             tr_x_xxx_xzz,   \
                             tr_x_xxx_yyy,   \
                             tr_x_xxx_yyz,   \
                             tr_x_xxx_yzz,   \
                             tr_x_xxx_zzz,   \
                             tr_x_xxxx_xx,   \
                             tr_x_xxxx_xxx,  \
                             tr_x_xxxx_xxy,  \
                             tr_x_xxxx_xxz,  \
                             tr_x_xxxx_xy,   \
                             tr_x_xxxx_xyy,  \
                             tr_x_xxxx_xyz,  \
                             tr_x_xxxx_xz,   \
                             tr_x_xxxx_xzz,  \
                             tr_x_xxxx_yy,   \
                             tr_x_xxxx_yyy,  \
                             tr_x_xxxx_yyz,  \
                             tr_x_xxxx_yz,   \
                             tr_x_xxxx_yzz,  \
                             tr_x_xxxx_zz,   \
                             tr_x_xxxx_zzz,  \
                             tr_x_xxxxx_xxx, \
                             tr_x_xxxxx_xxy, \
                             tr_x_xxxxx_xxz, \
                             tr_x_xxxxx_xyy, \
                             tr_x_xxxxx_xyz, \
                             tr_x_xxxxx_xzz, \
                             tr_x_xxxxx_yyy, \
                             tr_x_xxxxx_yyz, \
                             tr_x_xxxxx_yzz, \
                             tr_x_xxxxx_zzz, \
                             ts_xxxx_xxx,    \
                             ts_xxxx_xxy,    \
                             ts_xxxx_xxz,    \
                             ts_xxxx_xyy,    \
                             ts_xxxx_xyz,    \
                             ts_xxxx_xzz,    \
                             ts_xxxx_yyy,    \
                             ts_xxxx_yyz,    \
                             ts_xxxx_yzz,    \
                             ts_xxxx_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxx_xxx[i] = 4.0 * tr_x_xxx_xxx[i] * fe_0 + 3.0 * tr_x_xxxx_xx[i] * fe_0 + ts_xxxx_xxx[i] * fe_0 + tr_x_xxxx_xxx[i] * pa_x[i];

        tr_x_xxxxx_xxy[i] = 4.0 * tr_x_xxx_xxy[i] * fe_0 + 2.0 * tr_x_xxxx_xy[i] * fe_0 + ts_xxxx_xxy[i] * fe_0 + tr_x_xxxx_xxy[i] * pa_x[i];

        tr_x_xxxxx_xxz[i] = 4.0 * tr_x_xxx_xxz[i] * fe_0 + 2.0 * tr_x_xxxx_xz[i] * fe_0 + ts_xxxx_xxz[i] * fe_0 + tr_x_xxxx_xxz[i] * pa_x[i];

        tr_x_xxxxx_xyy[i] = 4.0 * tr_x_xxx_xyy[i] * fe_0 + tr_x_xxxx_yy[i] * fe_0 + ts_xxxx_xyy[i] * fe_0 + tr_x_xxxx_xyy[i] * pa_x[i];

        tr_x_xxxxx_xyz[i] = 4.0 * tr_x_xxx_xyz[i] * fe_0 + tr_x_xxxx_yz[i] * fe_0 + ts_xxxx_xyz[i] * fe_0 + tr_x_xxxx_xyz[i] * pa_x[i];

        tr_x_xxxxx_xzz[i] = 4.0 * tr_x_xxx_xzz[i] * fe_0 + tr_x_xxxx_zz[i] * fe_0 + ts_xxxx_xzz[i] * fe_0 + tr_x_xxxx_xzz[i] * pa_x[i];

        tr_x_xxxxx_yyy[i] = 4.0 * tr_x_xxx_yyy[i] * fe_0 + ts_xxxx_yyy[i] * fe_0 + tr_x_xxxx_yyy[i] * pa_x[i];

        tr_x_xxxxx_yyz[i] = 4.0 * tr_x_xxx_yyz[i] * fe_0 + ts_xxxx_yyz[i] * fe_0 + tr_x_xxxx_yyz[i] * pa_x[i];

        tr_x_xxxxx_yzz[i] = 4.0 * tr_x_xxx_yzz[i] * fe_0 + ts_xxxx_yzz[i] * fe_0 + tr_x_xxxx_yzz[i] * pa_x[i];

        tr_x_xxxxx_zzz[i] = 4.0 * tr_x_xxx_zzz[i] * fe_0 + ts_xxxx_zzz[i] * fe_0 + tr_x_xxxx_zzz[i] * pa_x[i];
    }

    // Set up 10-20 components of targeted buffer : HF

    auto tr_x_xxxxy_xxx = pbuffer.data(idx_dip_hf + 10);

    auto tr_x_xxxxy_xxy = pbuffer.data(idx_dip_hf + 11);

    auto tr_x_xxxxy_xxz = pbuffer.data(idx_dip_hf + 12);

    auto tr_x_xxxxy_xyy = pbuffer.data(idx_dip_hf + 13);

    auto tr_x_xxxxy_xyz = pbuffer.data(idx_dip_hf + 14);

    auto tr_x_xxxxy_xzz = pbuffer.data(idx_dip_hf + 15);

    auto tr_x_xxxxy_yyy = pbuffer.data(idx_dip_hf + 16);

    auto tr_x_xxxxy_yyz = pbuffer.data(idx_dip_hf + 17);

    auto tr_x_xxxxy_yzz = pbuffer.data(idx_dip_hf + 18);

    auto tr_x_xxxxy_zzz = pbuffer.data(idx_dip_hf + 19);

#pragma omp simd aligned(pa_y,               \
                             tr_x_xxxx_xx,   \
                             tr_x_xxxx_xxx,  \
                             tr_x_xxxx_xxy,  \
                             tr_x_xxxx_xxz,  \
                             tr_x_xxxx_xy,   \
                             tr_x_xxxx_xyy,  \
                             tr_x_xxxx_xyz,  \
                             tr_x_xxxx_xz,   \
                             tr_x_xxxx_xzz,  \
                             tr_x_xxxx_yy,   \
                             tr_x_xxxx_yyy,  \
                             tr_x_xxxx_yyz,  \
                             tr_x_xxxx_yz,   \
                             tr_x_xxxx_yzz,  \
                             tr_x_xxxx_zz,   \
                             tr_x_xxxx_zzz,  \
                             tr_x_xxxxy_xxx, \
                             tr_x_xxxxy_xxy, \
                             tr_x_xxxxy_xxz, \
                             tr_x_xxxxy_xyy, \
                             tr_x_xxxxy_xyz, \
                             tr_x_xxxxy_xzz, \
                             tr_x_xxxxy_yyy, \
                             tr_x_xxxxy_yyz, \
                             tr_x_xxxxy_yzz, \
                             tr_x_xxxxy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxy_xxx[i] = tr_x_xxxx_xxx[i] * pa_y[i];

        tr_x_xxxxy_xxy[i] = tr_x_xxxx_xx[i] * fe_0 + tr_x_xxxx_xxy[i] * pa_y[i];

        tr_x_xxxxy_xxz[i] = tr_x_xxxx_xxz[i] * pa_y[i];

        tr_x_xxxxy_xyy[i] = 2.0 * tr_x_xxxx_xy[i] * fe_0 + tr_x_xxxx_xyy[i] * pa_y[i];

        tr_x_xxxxy_xyz[i] = tr_x_xxxx_xz[i] * fe_0 + tr_x_xxxx_xyz[i] * pa_y[i];

        tr_x_xxxxy_xzz[i] = tr_x_xxxx_xzz[i] * pa_y[i];

        tr_x_xxxxy_yyy[i] = 3.0 * tr_x_xxxx_yy[i] * fe_0 + tr_x_xxxx_yyy[i] * pa_y[i];

        tr_x_xxxxy_yyz[i] = 2.0 * tr_x_xxxx_yz[i] * fe_0 + tr_x_xxxx_yyz[i] * pa_y[i];

        tr_x_xxxxy_yzz[i] = tr_x_xxxx_zz[i] * fe_0 + tr_x_xxxx_yzz[i] * pa_y[i];

        tr_x_xxxxy_zzz[i] = tr_x_xxxx_zzz[i] * pa_y[i];
    }

    // Set up 20-30 components of targeted buffer : HF

    auto tr_x_xxxxz_xxx = pbuffer.data(idx_dip_hf + 20);

    auto tr_x_xxxxz_xxy = pbuffer.data(idx_dip_hf + 21);

    auto tr_x_xxxxz_xxz = pbuffer.data(idx_dip_hf + 22);

    auto tr_x_xxxxz_xyy = pbuffer.data(idx_dip_hf + 23);

    auto tr_x_xxxxz_xyz = pbuffer.data(idx_dip_hf + 24);

    auto tr_x_xxxxz_xzz = pbuffer.data(idx_dip_hf + 25);

    auto tr_x_xxxxz_yyy = pbuffer.data(idx_dip_hf + 26);

    auto tr_x_xxxxz_yyz = pbuffer.data(idx_dip_hf + 27);

    auto tr_x_xxxxz_yzz = pbuffer.data(idx_dip_hf + 28);

    auto tr_x_xxxxz_zzz = pbuffer.data(idx_dip_hf + 29);

#pragma omp simd aligned(pa_z,               \
                             tr_x_xxxx_xx,   \
                             tr_x_xxxx_xxx,  \
                             tr_x_xxxx_xxy,  \
                             tr_x_xxxx_xxz,  \
                             tr_x_xxxx_xy,   \
                             tr_x_xxxx_xyy,  \
                             tr_x_xxxx_xyz,  \
                             tr_x_xxxx_xz,   \
                             tr_x_xxxx_xzz,  \
                             tr_x_xxxx_yy,   \
                             tr_x_xxxx_yyy,  \
                             tr_x_xxxx_yyz,  \
                             tr_x_xxxx_yz,   \
                             tr_x_xxxx_yzz,  \
                             tr_x_xxxx_zz,   \
                             tr_x_xxxx_zzz,  \
                             tr_x_xxxxz_xxx, \
                             tr_x_xxxxz_xxy, \
                             tr_x_xxxxz_xxz, \
                             tr_x_xxxxz_xyy, \
                             tr_x_xxxxz_xyz, \
                             tr_x_xxxxz_xzz, \
                             tr_x_xxxxz_yyy, \
                             tr_x_xxxxz_yyz, \
                             tr_x_xxxxz_yzz, \
                             tr_x_xxxxz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxxz_xxx[i] = tr_x_xxxx_xxx[i] * pa_z[i];

        tr_x_xxxxz_xxy[i] = tr_x_xxxx_xxy[i] * pa_z[i];

        tr_x_xxxxz_xxz[i] = tr_x_xxxx_xx[i] * fe_0 + tr_x_xxxx_xxz[i] * pa_z[i];

        tr_x_xxxxz_xyy[i] = tr_x_xxxx_xyy[i] * pa_z[i];

        tr_x_xxxxz_xyz[i] = tr_x_xxxx_xy[i] * fe_0 + tr_x_xxxx_xyz[i] * pa_z[i];

        tr_x_xxxxz_xzz[i] = 2.0 * tr_x_xxxx_xz[i] * fe_0 + tr_x_xxxx_xzz[i] * pa_z[i];

        tr_x_xxxxz_yyy[i] = tr_x_xxxx_yyy[i] * pa_z[i];

        tr_x_xxxxz_yyz[i] = tr_x_xxxx_yy[i] * fe_0 + tr_x_xxxx_yyz[i] * pa_z[i];

        tr_x_xxxxz_yzz[i] = 2.0 * tr_x_xxxx_yz[i] * fe_0 + tr_x_xxxx_yzz[i] * pa_z[i];

        tr_x_xxxxz_zzz[i] = 3.0 * tr_x_xxxx_zz[i] * fe_0 + tr_x_xxxx_zzz[i] * pa_z[i];
    }

    // Set up 30-40 components of targeted buffer : HF

    auto tr_x_xxxyy_xxx = pbuffer.data(idx_dip_hf + 30);

    auto tr_x_xxxyy_xxy = pbuffer.data(idx_dip_hf + 31);

    auto tr_x_xxxyy_xxz = pbuffer.data(idx_dip_hf + 32);

    auto tr_x_xxxyy_xyy = pbuffer.data(idx_dip_hf + 33);

    auto tr_x_xxxyy_xyz = pbuffer.data(idx_dip_hf + 34);

    auto tr_x_xxxyy_xzz = pbuffer.data(idx_dip_hf + 35);

    auto tr_x_xxxyy_yyy = pbuffer.data(idx_dip_hf + 36);

    auto tr_x_xxxyy_yyz = pbuffer.data(idx_dip_hf + 37);

    auto tr_x_xxxyy_yzz = pbuffer.data(idx_dip_hf + 38);

    auto tr_x_xxxyy_zzz = pbuffer.data(idx_dip_hf + 39);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_x_xxx_xxx,   \
                             tr_x_xxx_xxy,   \
                             tr_x_xxx_xxz,   \
                             tr_x_xxx_xyy,   \
                             tr_x_xxx_xyz,   \
                             tr_x_xxx_xzz,   \
                             tr_x_xxx_zzz,   \
                             tr_x_xxxy_xx,   \
                             tr_x_xxxy_xxx,  \
                             tr_x_xxxy_xxy,  \
                             tr_x_xxxy_xxz,  \
                             tr_x_xxxy_xy,   \
                             tr_x_xxxy_xyy,  \
                             tr_x_xxxy_xyz,  \
                             tr_x_xxxy_xz,   \
                             tr_x_xxxy_xzz,  \
                             tr_x_xxxy_zzz,  \
                             tr_x_xxxyy_xxx, \
                             tr_x_xxxyy_xxy, \
                             tr_x_xxxyy_xxz, \
                             tr_x_xxxyy_xyy, \
                             tr_x_xxxyy_xyz, \
                             tr_x_xxxyy_xzz, \
                             tr_x_xxxyy_yyy, \
                             tr_x_xxxyy_yyz, \
                             tr_x_xxxyy_yzz, \
                             tr_x_xxxyy_zzz, \
                             tr_x_xxyy_yyy,  \
                             tr_x_xxyy_yyz,  \
                             tr_x_xxyy_yzz,  \
                             tr_x_xyy_yyy,   \
                             tr_x_xyy_yyz,   \
                             tr_x_xyy_yzz,   \
                             ts_xxyy_yyy,    \
                             ts_xxyy_yyz,    \
                             ts_xxyy_yzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyy_xxx[i] = tr_x_xxx_xxx[i] * fe_0 + tr_x_xxxy_xxx[i] * pa_y[i];

        tr_x_xxxyy_xxy[i] = tr_x_xxx_xxy[i] * fe_0 + tr_x_xxxy_xx[i] * fe_0 + tr_x_xxxy_xxy[i] * pa_y[i];

        tr_x_xxxyy_xxz[i] = tr_x_xxx_xxz[i] * fe_0 + tr_x_xxxy_xxz[i] * pa_y[i];

        tr_x_xxxyy_xyy[i] = tr_x_xxx_xyy[i] * fe_0 + 2.0 * tr_x_xxxy_xy[i] * fe_0 + tr_x_xxxy_xyy[i] * pa_y[i];

        tr_x_xxxyy_xyz[i] = tr_x_xxx_xyz[i] * fe_0 + tr_x_xxxy_xz[i] * fe_0 + tr_x_xxxy_xyz[i] * pa_y[i];

        tr_x_xxxyy_xzz[i] = tr_x_xxx_xzz[i] * fe_0 + tr_x_xxxy_xzz[i] * pa_y[i];

        tr_x_xxxyy_yyy[i] = 2.0 * tr_x_xyy_yyy[i] * fe_0 + ts_xxyy_yyy[i] * fe_0 + tr_x_xxyy_yyy[i] * pa_x[i];

        tr_x_xxxyy_yyz[i] = 2.0 * tr_x_xyy_yyz[i] * fe_0 + ts_xxyy_yyz[i] * fe_0 + tr_x_xxyy_yyz[i] * pa_x[i];

        tr_x_xxxyy_yzz[i] = 2.0 * tr_x_xyy_yzz[i] * fe_0 + ts_xxyy_yzz[i] * fe_0 + tr_x_xxyy_yzz[i] * pa_x[i];

        tr_x_xxxyy_zzz[i] = tr_x_xxx_zzz[i] * fe_0 + tr_x_xxxy_zzz[i] * pa_y[i];
    }

    // Set up 40-50 components of targeted buffer : HF

    auto tr_x_xxxyz_xxx = pbuffer.data(idx_dip_hf + 40);

    auto tr_x_xxxyz_xxy = pbuffer.data(idx_dip_hf + 41);

    auto tr_x_xxxyz_xxz = pbuffer.data(idx_dip_hf + 42);

    auto tr_x_xxxyz_xyy = pbuffer.data(idx_dip_hf + 43);

    auto tr_x_xxxyz_xyz = pbuffer.data(idx_dip_hf + 44);

    auto tr_x_xxxyz_xzz = pbuffer.data(idx_dip_hf + 45);

    auto tr_x_xxxyz_yyy = pbuffer.data(idx_dip_hf + 46);

    auto tr_x_xxxyz_yyz = pbuffer.data(idx_dip_hf + 47);

    auto tr_x_xxxyz_yzz = pbuffer.data(idx_dip_hf + 48);

    auto tr_x_xxxyz_zzz = pbuffer.data(idx_dip_hf + 49);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_xxxy_xxy,  \
                             tr_x_xxxy_xyy,  \
                             tr_x_xxxy_yyy,  \
                             tr_x_xxxyz_xxx, \
                             tr_x_xxxyz_xxy, \
                             tr_x_xxxyz_xxz, \
                             tr_x_xxxyz_xyy, \
                             tr_x_xxxyz_xyz, \
                             tr_x_xxxyz_xzz, \
                             tr_x_xxxyz_yyy, \
                             tr_x_xxxyz_yyz, \
                             tr_x_xxxyz_yzz, \
                             tr_x_xxxyz_zzz, \
                             tr_x_xxxz_xxx,  \
                             tr_x_xxxz_xxz,  \
                             tr_x_xxxz_xyz,  \
                             tr_x_xxxz_xz,   \
                             tr_x_xxxz_xzz,  \
                             tr_x_xxxz_yyz,  \
                             tr_x_xxxz_yz,   \
                             tr_x_xxxz_yzz,  \
                             tr_x_xxxz_zz,   \
                             tr_x_xxxz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxyz_xxx[i] = tr_x_xxxz_xxx[i] * pa_y[i];

        tr_x_xxxyz_xxy[i] = tr_x_xxxy_xxy[i] * pa_z[i];

        tr_x_xxxyz_xxz[i] = tr_x_xxxz_xxz[i] * pa_y[i];

        tr_x_xxxyz_xyy[i] = tr_x_xxxy_xyy[i] * pa_z[i];

        tr_x_xxxyz_xyz[i] = tr_x_xxxz_xz[i] * fe_0 + tr_x_xxxz_xyz[i] * pa_y[i];

        tr_x_xxxyz_xzz[i] = tr_x_xxxz_xzz[i] * pa_y[i];

        tr_x_xxxyz_yyy[i] = tr_x_xxxy_yyy[i] * pa_z[i];

        tr_x_xxxyz_yyz[i] = 2.0 * tr_x_xxxz_yz[i] * fe_0 + tr_x_xxxz_yyz[i] * pa_y[i];

        tr_x_xxxyz_yzz[i] = tr_x_xxxz_zz[i] * fe_0 + tr_x_xxxz_yzz[i] * pa_y[i];

        tr_x_xxxyz_zzz[i] = tr_x_xxxz_zzz[i] * pa_y[i];
    }

    // Set up 50-60 components of targeted buffer : HF

    auto tr_x_xxxzz_xxx = pbuffer.data(idx_dip_hf + 50);

    auto tr_x_xxxzz_xxy = pbuffer.data(idx_dip_hf + 51);

    auto tr_x_xxxzz_xxz = pbuffer.data(idx_dip_hf + 52);

    auto tr_x_xxxzz_xyy = pbuffer.data(idx_dip_hf + 53);

    auto tr_x_xxxzz_xyz = pbuffer.data(idx_dip_hf + 54);

    auto tr_x_xxxzz_xzz = pbuffer.data(idx_dip_hf + 55);

    auto tr_x_xxxzz_yyy = pbuffer.data(idx_dip_hf + 56);

    auto tr_x_xxxzz_yyz = pbuffer.data(idx_dip_hf + 57);

    auto tr_x_xxxzz_yzz = pbuffer.data(idx_dip_hf + 58);

    auto tr_x_xxxzz_zzz = pbuffer.data(idx_dip_hf + 59);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_x_xxx_xxx,   \
                             tr_x_xxx_xxy,   \
                             tr_x_xxx_xxz,   \
                             tr_x_xxx_xyy,   \
                             tr_x_xxx_xyz,   \
                             tr_x_xxx_xzz,   \
                             tr_x_xxx_yyy,   \
                             tr_x_xxxz_xx,   \
                             tr_x_xxxz_xxx,  \
                             tr_x_xxxz_xxy,  \
                             tr_x_xxxz_xxz,  \
                             tr_x_xxxz_xy,   \
                             tr_x_xxxz_xyy,  \
                             tr_x_xxxz_xyz,  \
                             tr_x_xxxz_xz,   \
                             tr_x_xxxz_xzz,  \
                             tr_x_xxxz_yyy,  \
                             tr_x_xxxzz_xxx, \
                             tr_x_xxxzz_xxy, \
                             tr_x_xxxzz_xxz, \
                             tr_x_xxxzz_xyy, \
                             tr_x_xxxzz_xyz, \
                             tr_x_xxxzz_xzz, \
                             tr_x_xxxzz_yyy, \
                             tr_x_xxxzz_yyz, \
                             tr_x_xxxzz_yzz, \
                             tr_x_xxxzz_zzz, \
                             tr_x_xxzz_yyz,  \
                             tr_x_xxzz_yzz,  \
                             tr_x_xxzz_zzz,  \
                             tr_x_xzz_yyz,   \
                             tr_x_xzz_yzz,   \
                             tr_x_xzz_zzz,   \
                             ts_xxzz_yyz,    \
                             ts_xxzz_yzz,    \
                             ts_xxzz_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxxzz_xxx[i] = tr_x_xxx_xxx[i] * fe_0 + tr_x_xxxz_xxx[i] * pa_z[i];

        tr_x_xxxzz_xxy[i] = tr_x_xxx_xxy[i] * fe_0 + tr_x_xxxz_xxy[i] * pa_z[i];

        tr_x_xxxzz_xxz[i] = tr_x_xxx_xxz[i] * fe_0 + tr_x_xxxz_xx[i] * fe_0 + tr_x_xxxz_xxz[i] * pa_z[i];

        tr_x_xxxzz_xyy[i] = tr_x_xxx_xyy[i] * fe_0 + tr_x_xxxz_xyy[i] * pa_z[i];

        tr_x_xxxzz_xyz[i] = tr_x_xxx_xyz[i] * fe_0 + tr_x_xxxz_xy[i] * fe_0 + tr_x_xxxz_xyz[i] * pa_z[i];

        tr_x_xxxzz_xzz[i] = tr_x_xxx_xzz[i] * fe_0 + 2.0 * tr_x_xxxz_xz[i] * fe_0 + tr_x_xxxz_xzz[i] * pa_z[i];

        tr_x_xxxzz_yyy[i] = tr_x_xxx_yyy[i] * fe_0 + tr_x_xxxz_yyy[i] * pa_z[i];

        tr_x_xxxzz_yyz[i] = 2.0 * tr_x_xzz_yyz[i] * fe_0 + ts_xxzz_yyz[i] * fe_0 + tr_x_xxzz_yyz[i] * pa_x[i];

        tr_x_xxxzz_yzz[i] = 2.0 * tr_x_xzz_yzz[i] * fe_0 + ts_xxzz_yzz[i] * fe_0 + tr_x_xxzz_yzz[i] * pa_x[i];

        tr_x_xxxzz_zzz[i] = 2.0 * tr_x_xzz_zzz[i] * fe_0 + ts_xxzz_zzz[i] * fe_0 + tr_x_xxzz_zzz[i] * pa_x[i];
    }

    // Set up 60-70 components of targeted buffer : HF

    auto tr_x_xxyyy_xxx = pbuffer.data(idx_dip_hf + 60);

    auto tr_x_xxyyy_xxy = pbuffer.data(idx_dip_hf + 61);

    auto tr_x_xxyyy_xxz = pbuffer.data(idx_dip_hf + 62);

    auto tr_x_xxyyy_xyy = pbuffer.data(idx_dip_hf + 63);

    auto tr_x_xxyyy_xyz = pbuffer.data(idx_dip_hf + 64);

    auto tr_x_xxyyy_xzz = pbuffer.data(idx_dip_hf + 65);

    auto tr_x_xxyyy_yyy = pbuffer.data(idx_dip_hf + 66);

    auto tr_x_xxyyy_yyz = pbuffer.data(idx_dip_hf + 67);

    auto tr_x_xxyyy_yzz = pbuffer.data(idx_dip_hf + 68);

    auto tr_x_xxyyy_zzz = pbuffer.data(idx_dip_hf + 69);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_x_xxy_xxx,   \
                             tr_x_xxy_xxy,   \
                             tr_x_xxy_xxz,   \
                             tr_x_xxy_xyy,   \
                             tr_x_xxy_xyz,   \
                             tr_x_xxy_xzz,   \
                             tr_x_xxy_zzz,   \
                             tr_x_xxyy_xx,   \
                             tr_x_xxyy_xxx,  \
                             tr_x_xxyy_xxy,  \
                             tr_x_xxyy_xxz,  \
                             tr_x_xxyy_xy,   \
                             tr_x_xxyy_xyy,  \
                             tr_x_xxyy_xyz,  \
                             tr_x_xxyy_xz,   \
                             tr_x_xxyy_xzz,  \
                             tr_x_xxyy_zzz,  \
                             tr_x_xxyyy_xxx, \
                             tr_x_xxyyy_xxy, \
                             tr_x_xxyyy_xxz, \
                             tr_x_xxyyy_xyy, \
                             tr_x_xxyyy_xyz, \
                             tr_x_xxyyy_xzz, \
                             tr_x_xxyyy_yyy, \
                             tr_x_xxyyy_yyz, \
                             tr_x_xxyyy_yzz, \
                             tr_x_xxyyy_zzz, \
                             tr_x_xyyy_yyy,  \
                             tr_x_xyyy_yyz,  \
                             tr_x_xyyy_yzz,  \
                             tr_x_yyy_yyy,   \
                             tr_x_yyy_yyz,   \
                             tr_x_yyy_yzz,   \
                             ts_xyyy_yyy,    \
                             ts_xyyy_yyz,    \
                             ts_xyyy_yzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyy_xxx[i] = 2.0 * tr_x_xxy_xxx[i] * fe_0 + tr_x_xxyy_xxx[i] * pa_y[i];

        tr_x_xxyyy_xxy[i] = 2.0 * tr_x_xxy_xxy[i] * fe_0 + tr_x_xxyy_xx[i] * fe_0 + tr_x_xxyy_xxy[i] * pa_y[i];

        tr_x_xxyyy_xxz[i] = 2.0 * tr_x_xxy_xxz[i] * fe_0 + tr_x_xxyy_xxz[i] * pa_y[i];

        tr_x_xxyyy_xyy[i] = 2.0 * tr_x_xxy_xyy[i] * fe_0 + 2.0 * tr_x_xxyy_xy[i] * fe_0 + tr_x_xxyy_xyy[i] * pa_y[i];

        tr_x_xxyyy_xyz[i] = 2.0 * tr_x_xxy_xyz[i] * fe_0 + tr_x_xxyy_xz[i] * fe_0 + tr_x_xxyy_xyz[i] * pa_y[i];

        tr_x_xxyyy_xzz[i] = 2.0 * tr_x_xxy_xzz[i] * fe_0 + tr_x_xxyy_xzz[i] * pa_y[i];

        tr_x_xxyyy_yyy[i] = tr_x_yyy_yyy[i] * fe_0 + ts_xyyy_yyy[i] * fe_0 + tr_x_xyyy_yyy[i] * pa_x[i];

        tr_x_xxyyy_yyz[i] = tr_x_yyy_yyz[i] * fe_0 + ts_xyyy_yyz[i] * fe_0 + tr_x_xyyy_yyz[i] * pa_x[i];

        tr_x_xxyyy_yzz[i] = tr_x_yyy_yzz[i] * fe_0 + ts_xyyy_yzz[i] * fe_0 + tr_x_xyyy_yzz[i] * pa_x[i];

        tr_x_xxyyy_zzz[i] = 2.0 * tr_x_xxy_zzz[i] * fe_0 + tr_x_xxyy_zzz[i] * pa_y[i];
    }

    // Set up 70-80 components of targeted buffer : HF

    auto tr_x_xxyyz_xxx = pbuffer.data(idx_dip_hf + 70);

    auto tr_x_xxyyz_xxy = pbuffer.data(idx_dip_hf + 71);

    auto tr_x_xxyyz_xxz = pbuffer.data(idx_dip_hf + 72);

    auto tr_x_xxyyz_xyy = pbuffer.data(idx_dip_hf + 73);

    auto tr_x_xxyyz_xyz = pbuffer.data(idx_dip_hf + 74);

    auto tr_x_xxyyz_xzz = pbuffer.data(idx_dip_hf + 75);

    auto tr_x_xxyyz_yyy = pbuffer.data(idx_dip_hf + 76);

    auto tr_x_xxyyz_yyz = pbuffer.data(idx_dip_hf + 77);

    auto tr_x_xxyyz_yzz = pbuffer.data(idx_dip_hf + 78);

    auto tr_x_xxyyz_zzz = pbuffer.data(idx_dip_hf + 79);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_xxyy_xxx,  \
                             tr_x_xxyy_xxy,  \
                             tr_x_xxyy_xy,   \
                             tr_x_xxyy_xyy,  \
                             tr_x_xxyy_xyz,  \
                             tr_x_xxyy_yy,   \
                             tr_x_xxyy_yyy,  \
                             tr_x_xxyy_yyz,  \
                             tr_x_xxyy_yz,   \
                             tr_x_xxyy_yzz,  \
                             tr_x_xxyyz_xxx, \
                             tr_x_xxyyz_xxy, \
                             tr_x_xxyyz_xxz, \
                             tr_x_xxyyz_xyy, \
                             tr_x_xxyyz_xyz, \
                             tr_x_xxyyz_xzz, \
                             tr_x_xxyyz_yyy, \
                             tr_x_xxyyz_yyz, \
                             tr_x_xxyyz_yzz, \
                             tr_x_xxyyz_zzz, \
                             tr_x_xxyz_xxz,  \
                             tr_x_xxyz_xzz,  \
                             tr_x_xxyz_zzz,  \
                             tr_x_xxz_xxz,   \
                             tr_x_xxz_xzz,   \
                             tr_x_xxz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyyz_xxx[i] = tr_x_xxyy_xxx[i] * pa_z[i];

        tr_x_xxyyz_xxy[i] = tr_x_xxyy_xxy[i] * pa_z[i];

        tr_x_xxyyz_xxz[i] = tr_x_xxz_xxz[i] * fe_0 + tr_x_xxyz_xxz[i] * pa_y[i];

        tr_x_xxyyz_xyy[i] = tr_x_xxyy_xyy[i] * pa_z[i];

        tr_x_xxyyz_xyz[i] = tr_x_xxyy_xy[i] * fe_0 + tr_x_xxyy_xyz[i] * pa_z[i];

        tr_x_xxyyz_xzz[i] = tr_x_xxz_xzz[i] * fe_0 + tr_x_xxyz_xzz[i] * pa_y[i];

        tr_x_xxyyz_yyy[i] = tr_x_xxyy_yyy[i] * pa_z[i];

        tr_x_xxyyz_yyz[i] = tr_x_xxyy_yy[i] * fe_0 + tr_x_xxyy_yyz[i] * pa_z[i];

        tr_x_xxyyz_yzz[i] = 2.0 * tr_x_xxyy_yz[i] * fe_0 + tr_x_xxyy_yzz[i] * pa_z[i];

        tr_x_xxyyz_zzz[i] = tr_x_xxz_zzz[i] * fe_0 + tr_x_xxyz_zzz[i] * pa_y[i];
    }

    // Set up 80-90 components of targeted buffer : HF

    auto tr_x_xxyzz_xxx = pbuffer.data(idx_dip_hf + 80);

    auto tr_x_xxyzz_xxy = pbuffer.data(idx_dip_hf + 81);

    auto tr_x_xxyzz_xxz = pbuffer.data(idx_dip_hf + 82);

    auto tr_x_xxyzz_xyy = pbuffer.data(idx_dip_hf + 83);

    auto tr_x_xxyzz_xyz = pbuffer.data(idx_dip_hf + 84);

    auto tr_x_xxyzz_xzz = pbuffer.data(idx_dip_hf + 85);

    auto tr_x_xxyzz_yyy = pbuffer.data(idx_dip_hf + 86);

    auto tr_x_xxyzz_yyz = pbuffer.data(idx_dip_hf + 87);

    auto tr_x_xxyzz_yzz = pbuffer.data(idx_dip_hf + 88);

    auto tr_x_xxyzz_zzz = pbuffer.data(idx_dip_hf + 89);

#pragma omp simd aligned(pa_y,               \
                             tr_x_xxyzz_xxx, \
                             tr_x_xxyzz_xxy, \
                             tr_x_xxyzz_xxz, \
                             tr_x_xxyzz_xyy, \
                             tr_x_xxyzz_xyz, \
                             tr_x_xxyzz_xzz, \
                             tr_x_xxyzz_yyy, \
                             tr_x_xxyzz_yyz, \
                             tr_x_xxyzz_yzz, \
                             tr_x_xxyzz_zzz, \
                             tr_x_xxzz_xx,   \
                             tr_x_xxzz_xxx,  \
                             tr_x_xxzz_xxy,  \
                             tr_x_xxzz_xxz,  \
                             tr_x_xxzz_xy,   \
                             tr_x_xxzz_xyy,  \
                             tr_x_xxzz_xyz,  \
                             tr_x_xxzz_xz,   \
                             tr_x_xxzz_xzz,  \
                             tr_x_xxzz_yy,   \
                             tr_x_xxzz_yyy,  \
                             tr_x_xxzz_yyz,  \
                             tr_x_xxzz_yz,   \
                             tr_x_xxzz_yzz,  \
                             tr_x_xxzz_zz,   \
                             tr_x_xxzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxyzz_xxx[i] = tr_x_xxzz_xxx[i] * pa_y[i];

        tr_x_xxyzz_xxy[i] = tr_x_xxzz_xx[i] * fe_0 + tr_x_xxzz_xxy[i] * pa_y[i];

        tr_x_xxyzz_xxz[i] = tr_x_xxzz_xxz[i] * pa_y[i];

        tr_x_xxyzz_xyy[i] = 2.0 * tr_x_xxzz_xy[i] * fe_0 + tr_x_xxzz_xyy[i] * pa_y[i];

        tr_x_xxyzz_xyz[i] = tr_x_xxzz_xz[i] * fe_0 + tr_x_xxzz_xyz[i] * pa_y[i];

        tr_x_xxyzz_xzz[i] = tr_x_xxzz_xzz[i] * pa_y[i];

        tr_x_xxyzz_yyy[i] = 3.0 * tr_x_xxzz_yy[i] * fe_0 + tr_x_xxzz_yyy[i] * pa_y[i];

        tr_x_xxyzz_yyz[i] = 2.0 * tr_x_xxzz_yz[i] * fe_0 + tr_x_xxzz_yyz[i] * pa_y[i];

        tr_x_xxyzz_yzz[i] = tr_x_xxzz_zz[i] * fe_0 + tr_x_xxzz_yzz[i] * pa_y[i];

        tr_x_xxyzz_zzz[i] = tr_x_xxzz_zzz[i] * pa_y[i];
    }

    // Set up 90-100 components of targeted buffer : HF

    auto tr_x_xxzzz_xxx = pbuffer.data(idx_dip_hf + 90);

    auto tr_x_xxzzz_xxy = pbuffer.data(idx_dip_hf + 91);

    auto tr_x_xxzzz_xxz = pbuffer.data(idx_dip_hf + 92);

    auto tr_x_xxzzz_xyy = pbuffer.data(idx_dip_hf + 93);

    auto tr_x_xxzzz_xyz = pbuffer.data(idx_dip_hf + 94);

    auto tr_x_xxzzz_xzz = pbuffer.data(idx_dip_hf + 95);

    auto tr_x_xxzzz_yyy = pbuffer.data(idx_dip_hf + 96);

    auto tr_x_xxzzz_yyz = pbuffer.data(idx_dip_hf + 97);

    auto tr_x_xxzzz_yzz = pbuffer.data(idx_dip_hf + 98);

    auto tr_x_xxzzz_zzz = pbuffer.data(idx_dip_hf + 99);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_x_xxz_xxx,   \
                             tr_x_xxz_xxy,   \
                             tr_x_xxz_xxz,   \
                             tr_x_xxz_xyy,   \
                             tr_x_xxz_xyz,   \
                             tr_x_xxz_xzz,   \
                             tr_x_xxz_yyy,   \
                             tr_x_xxzz_xx,   \
                             tr_x_xxzz_xxx,  \
                             tr_x_xxzz_xxy,  \
                             tr_x_xxzz_xxz,  \
                             tr_x_xxzz_xy,   \
                             tr_x_xxzz_xyy,  \
                             tr_x_xxzz_xyz,  \
                             tr_x_xxzz_xz,   \
                             tr_x_xxzz_xzz,  \
                             tr_x_xxzz_yyy,  \
                             tr_x_xxzzz_xxx, \
                             tr_x_xxzzz_xxy, \
                             tr_x_xxzzz_xxz, \
                             tr_x_xxzzz_xyy, \
                             tr_x_xxzzz_xyz, \
                             tr_x_xxzzz_xzz, \
                             tr_x_xxzzz_yyy, \
                             tr_x_xxzzz_yyz, \
                             tr_x_xxzzz_yzz, \
                             tr_x_xxzzz_zzz, \
                             tr_x_xzzz_yyz,  \
                             tr_x_xzzz_yzz,  \
                             tr_x_xzzz_zzz,  \
                             tr_x_zzz_yyz,   \
                             tr_x_zzz_yzz,   \
                             tr_x_zzz_zzz,   \
                             ts_xzzz_yyz,    \
                             ts_xzzz_yzz,    \
                             ts_xzzz_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxzzz_xxx[i] = 2.0 * tr_x_xxz_xxx[i] * fe_0 + tr_x_xxzz_xxx[i] * pa_z[i];

        tr_x_xxzzz_xxy[i] = 2.0 * tr_x_xxz_xxy[i] * fe_0 + tr_x_xxzz_xxy[i] * pa_z[i];

        tr_x_xxzzz_xxz[i] = 2.0 * tr_x_xxz_xxz[i] * fe_0 + tr_x_xxzz_xx[i] * fe_0 + tr_x_xxzz_xxz[i] * pa_z[i];

        tr_x_xxzzz_xyy[i] = 2.0 * tr_x_xxz_xyy[i] * fe_0 + tr_x_xxzz_xyy[i] * pa_z[i];

        tr_x_xxzzz_xyz[i] = 2.0 * tr_x_xxz_xyz[i] * fe_0 + tr_x_xxzz_xy[i] * fe_0 + tr_x_xxzz_xyz[i] * pa_z[i];

        tr_x_xxzzz_xzz[i] = 2.0 * tr_x_xxz_xzz[i] * fe_0 + 2.0 * tr_x_xxzz_xz[i] * fe_0 + tr_x_xxzz_xzz[i] * pa_z[i];

        tr_x_xxzzz_yyy[i] = 2.0 * tr_x_xxz_yyy[i] * fe_0 + tr_x_xxzz_yyy[i] * pa_z[i];

        tr_x_xxzzz_yyz[i] = tr_x_zzz_yyz[i] * fe_0 + ts_xzzz_yyz[i] * fe_0 + tr_x_xzzz_yyz[i] * pa_x[i];

        tr_x_xxzzz_yzz[i] = tr_x_zzz_yzz[i] * fe_0 + ts_xzzz_yzz[i] * fe_0 + tr_x_xzzz_yzz[i] * pa_x[i];

        tr_x_xxzzz_zzz[i] = tr_x_zzz_zzz[i] * fe_0 + ts_xzzz_zzz[i] * fe_0 + tr_x_xzzz_zzz[i] * pa_x[i];
    }

    // Set up 100-110 components of targeted buffer : HF

    auto tr_x_xyyyy_xxx = pbuffer.data(idx_dip_hf + 100);

    auto tr_x_xyyyy_xxy = pbuffer.data(idx_dip_hf + 101);

    auto tr_x_xyyyy_xxz = pbuffer.data(idx_dip_hf + 102);

    auto tr_x_xyyyy_xyy = pbuffer.data(idx_dip_hf + 103);

    auto tr_x_xyyyy_xyz = pbuffer.data(idx_dip_hf + 104);

    auto tr_x_xyyyy_xzz = pbuffer.data(idx_dip_hf + 105);

    auto tr_x_xyyyy_yyy = pbuffer.data(idx_dip_hf + 106);

    auto tr_x_xyyyy_yyz = pbuffer.data(idx_dip_hf + 107);

    auto tr_x_xyyyy_yzz = pbuffer.data(idx_dip_hf + 108);

    auto tr_x_xyyyy_zzz = pbuffer.data(idx_dip_hf + 109);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_x_xyy_xxx,   \
                             tr_x_xyy_xxz,   \
                             tr_x_xyy_xzz,   \
                             tr_x_xyyy_xxx,  \
                             tr_x_xyyy_xxz,  \
                             tr_x_xyyy_xzz,  \
                             tr_x_xyyyy_xxx, \
                             tr_x_xyyyy_xxy, \
                             tr_x_xyyyy_xxz, \
                             tr_x_xyyyy_xyy, \
                             tr_x_xyyyy_xyz, \
                             tr_x_xyyyy_xzz, \
                             tr_x_xyyyy_yyy, \
                             tr_x_xyyyy_yyz, \
                             tr_x_xyyyy_yzz, \
                             tr_x_xyyyy_zzz, \
                             tr_x_yyyy_xxy,  \
                             tr_x_yyyy_xy,   \
                             tr_x_yyyy_xyy,  \
                             tr_x_yyyy_xyz,  \
                             tr_x_yyyy_yy,   \
                             tr_x_yyyy_yyy,  \
                             tr_x_yyyy_yyz,  \
                             tr_x_yyyy_yz,   \
                             tr_x_yyyy_yzz,  \
                             tr_x_yyyy_zzz,  \
                             ts_yyyy_xxy,    \
                             ts_yyyy_xyy,    \
                             ts_yyyy_xyz,    \
                             ts_yyyy_yyy,    \
                             ts_yyyy_yyz,    \
                             ts_yyyy_yzz,    \
                             ts_yyyy_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyy_xxx[i] = 3.0 * tr_x_xyy_xxx[i] * fe_0 + tr_x_xyyy_xxx[i] * pa_y[i];

        tr_x_xyyyy_xxy[i] = 2.0 * tr_x_yyyy_xy[i] * fe_0 + ts_yyyy_xxy[i] * fe_0 + tr_x_yyyy_xxy[i] * pa_x[i];

        tr_x_xyyyy_xxz[i] = 3.0 * tr_x_xyy_xxz[i] * fe_0 + tr_x_xyyy_xxz[i] * pa_y[i];

        tr_x_xyyyy_xyy[i] = tr_x_yyyy_yy[i] * fe_0 + ts_yyyy_xyy[i] * fe_0 + tr_x_yyyy_xyy[i] * pa_x[i];

        tr_x_xyyyy_xyz[i] = tr_x_yyyy_yz[i] * fe_0 + ts_yyyy_xyz[i] * fe_0 + tr_x_yyyy_xyz[i] * pa_x[i];

        tr_x_xyyyy_xzz[i] = 3.0 * tr_x_xyy_xzz[i] * fe_0 + tr_x_xyyy_xzz[i] * pa_y[i];

        tr_x_xyyyy_yyy[i] = ts_yyyy_yyy[i] * fe_0 + tr_x_yyyy_yyy[i] * pa_x[i];

        tr_x_xyyyy_yyz[i] = ts_yyyy_yyz[i] * fe_0 + tr_x_yyyy_yyz[i] * pa_x[i];

        tr_x_xyyyy_yzz[i] = ts_yyyy_yzz[i] * fe_0 + tr_x_yyyy_yzz[i] * pa_x[i];

        tr_x_xyyyy_zzz[i] = ts_yyyy_zzz[i] * fe_0 + tr_x_yyyy_zzz[i] * pa_x[i];
    }

    // Set up 110-120 components of targeted buffer : HF

    auto tr_x_xyyyz_xxx = pbuffer.data(idx_dip_hf + 110);

    auto tr_x_xyyyz_xxy = pbuffer.data(idx_dip_hf + 111);

    auto tr_x_xyyyz_xxz = pbuffer.data(idx_dip_hf + 112);

    auto tr_x_xyyyz_xyy = pbuffer.data(idx_dip_hf + 113);

    auto tr_x_xyyyz_xyz = pbuffer.data(idx_dip_hf + 114);

    auto tr_x_xyyyz_xzz = pbuffer.data(idx_dip_hf + 115);

    auto tr_x_xyyyz_yyy = pbuffer.data(idx_dip_hf + 116);

    auto tr_x_xyyyz_yyz = pbuffer.data(idx_dip_hf + 117);

    auto tr_x_xyyyz_yzz = pbuffer.data(idx_dip_hf + 118);

    auto tr_x_xyyyz_zzz = pbuffer.data(idx_dip_hf + 119);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_x_xyyy_xxx,  \
                             tr_x_xyyy_xxy,  \
                             tr_x_xyyy_xy,   \
                             tr_x_xyyy_xyy,  \
                             tr_x_xyyy_xyz,  \
                             tr_x_xyyy_yyy,  \
                             tr_x_xyyyz_xxx, \
                             tr_x_xyyyz_xxy, \
                             tr_x_xyyyz_xxz, \
                             tr_x_xyyyz_xyy, \
                             tr_x_xyyyz_xyz, \
                             tr_x_xyyyz_xzz, \
                             tr_x_xyyyz_yyy, \
                             tr_x_xyyyz_yyz, \
                             tr_x_xyyyz_yzz, \
                             tr_x_xyyyz_zzz, \
                             tr_x_xyyz_xxz,  \
                             tr_x_xyyz_xzz,  \
                             tr_x_xyz_xxz,   \
                             tr_x_xyz_xzz,   \
                             tr_x_yyyz_yyz,  \
                             tr_x_yyyz_yzz,  \
                             tr_x_yyyz_zzz,  \
                             ts_yyyz_yyz,    \
                             ts_yyyz_yzz,    \
                             ts_yyyz_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyyz_xxx[i] = tr_x_xyyy_xxx[i] * pa_z[i];

        tr_x_xyyyz_xxy[i] = tr_x_xyyy_xxy[i] * pa_z[i];

        tr_x_xyyyz_xxz[i] = 2.0 * tr_x_xyz_xxz[i] * fe_0 + tr_x_xyyz_xxz[i] * pa_y[i];

        tr_x_xyyyz_xyy[i] = tr_x_xyyy_xyy[i] * pa_z[i];

        tr_x_xyyyz_xyz[i] = tr_x_xyyy_xy[i] * fe_0 + tr_x_xyyy_xyz[i] * pa_z[i];

        tr_x_xyyyz_xzz[i] = 2.0 * tr_x_xyz_xzz[i] * fe_0 + tr_x_xyyz_xzz[i] * pa_y[i];

        tr_x_xyyyz_yyy[i] = tr_x_xyyy_yyy[i] * pa_z[i];

        tr_x_xyyyz_yyz[i] = ts_yyyz_yyz[i] * fe_0 + tr_x_yyyz_yyz[i] * pa_x[i];

        tr_x_xyyyz_yzz[i] = ts_yyyz_yzz[i] * fe_0 + tr_x_yyyz_yzz[i] * pa_x[i];

        tr_x_xyyyz_zzz[i] = ts_yyyz_zzz[i] * fe_0 + tr_x_yyyz_zzz[i] * pa_x[i];
    }

    // Set up 120-130 components of targeted buffer : HF

    auto tr_x_xyyzz_xxx = pbuffer.data(idx_dip_hf + 120);

    auto tr_x_xyyzz_xxy = pbuffer.data(idx_dip_hf + 121);

    auto tr_x_xyyzz_xxz = pbuffer.data(idx_dip_hf + 122);

    auto tr_x_xyyzz_xyy = pbuffer.data(idx_dip_hf + 123);

    auto tr_x_xyyzz_xyz = pbuffer.data(idx_dip_hf + 124);

    auto tr_x_xyyzz_xzz = pbuffer.data(idx_dip_hf + 125);

    auto tr_x_xyyzz_yyy = pbuffer.data(idx_dip_hf + 126);

    auto tr_x_xyyzz_yyz = pbuffer.data(idx_dip_hf + 127);

    auto tr_x_xyyzz_yzz = pbuffer.data(idx_dip_hf + 128);

    auto tr_x_xyyzz_zzz = pbuffer.data(idx_dip_hf + 129);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_x_xyy_xxy,   \
                             tr_x_xyy_xyy,   \
                             tr_x_xyyz_xxy,  \
                             tr_x_xyyz_xyy,  \
                             tr_x_xyyzz_xxx, \
                             tr_x_xyyzz_xxy, \
                             tr_x_xyyzz_xxz, \
                             tr_x_xyyzz_xyy, \
                             tr_x_xyyzz_xyz, \
                             tr_x_xyyzz_xzz, \
                             tr_x_xyyzz_yyy, \
                             tr_x_xyyzz_yyz, \
                             tr_x_xyyzz_yzz, \
                             tr_x_xyyzz_zzz, \
                             tr_x_xyzz_xxx,  \
                             tr_x_xyzz_xxz,  \
                             tr_x_xyzz_xzz,  \
                             tr_x_xzz_xxx,   \
                             tr_x_xzz_xxz,   \
                             tr_x_xzz_xzz,   \
                             tr_x_yyzz_xyz,  \
                             tr_x_yyzz_yyy,  \
                             tr_x_yyzz_yyz,  \
                             tr_x_yyzz_yz,   \
                             tr_x_yyzz_yzz,  \
                             tr_x_yyzz_zzz,  \
                             ts_yyzz_xyz,    \
                             ts_yyzz_yyy,    \
                             ts_yyzz_yyz,    \
                             ts_yyzz_yzz,    \
                             ts_yyzz_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyyzz_xxx[i] = tr_x_xzz_xxx[i] * fe_0 + tr_x_xyzz_xxx[i] * pa_y[i];

        tr_x_xyyzz_xxy[i] = tr_x_xyy_xxy[i] * fe_0 + tr_x_xyyz_xxy[i] * pa_z[i];

        tr_x_xyyzz_xxz[i] = tr_x_xzz_xxz[i] * fe_0 + tr_x_xyzz_xxz[i] * pa_y[i];

        tr_x_xyyzz_xyy[i] = tr_x_xyy_xyy[i] * fe_0 + tr_x_xyyz_xyy[i] * pa_z[i];

        tr_x_xyyzz_xyz[i] = tr_x_yyzz_yz[i] * fe_0 + ts_yyzz_xyz[i] * fe_0 + tr_x_yyzz_xyz[i] * pa_x[i];

        tr_x_xyyzz_xzz[i] = tr_x_xzz_xzz[i] * fe_0 + tr_x_xyzz_xzz[i] * pa_y[i];

        tr_x_xyyzz_yyy[i] = ts_yyzz_yyy[i] * fe_0 + tr_x_yyzz_yyy[i] * pa_x[i];

        tr_x_xyyzz_yyz[i] = ts_yyzz_yyz[i] * fe_0 + tr_x_yyzz_yyz[i] * pa_x[i];

        tr_x_xyyzz_yzz[i] = ts_yyzz_yzz[i] * fe_0 + tr_x_yyzz_yzz[i] * pa_x[i];

        tr_x_xyyzz_zzz[i] = ts_yyzz_zzz[i] * fe_0 + tr_x_yyzz_zzz[i] * pa_x[i];
    }

    // Set up 130-140 components of targeted buffer : HF

    auto tr_x_xyzzz_xxx = pbuffer.data(idx_dip_hf + 130);

    auto tr_x_xyzzz_xxy = pbuffer.data(idx_dip_hf + 131);

    auto tr_x_xyzzz_xxz = pbuffer.data(idx_dip_hf + 132);

    auto tr_x_xyzzz_xyy = pbuffer.data(idx_dip_hf + 133);

    auto tr_x_xyzzz_xyz = pbuffer.data(idx_dip_hf + 134);

    auto tr_x_xyzzz_xzz = pbuffer.data(idx_dip_hf + 135);

    auto tr_x_xyzzz_yyy = pbuffer.data(idx_dip_hf + 136);

    auto tr_x_xyzzz_yyz = pbuffer.data(idx_dip_hf + 137);

    auto tr_x_xyzzz_yzz = pbuffer.data(idx_dip_hf + 138);

    auto tr_x_xyzzz_zzz = pbuffer.data(idx_dip_hf + 139);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_x_xyzzz_xxx, \
                             tr_x_xyzzz_xxy, \
                             tr_x_xyzzz_xxz, \
                             tr_x_xyzzz_xyy, \
                             tr_x_xyzzz_xyz, \
                             tr_x_xyzzz_xzz, \
                             tr_x_xyzzz_yyy, \
                             tr_x_xyzzz_yyz, \
                             tr_x_xyzzz_yzz, \
                             tr_x_xyzzz_zzz, \
                             tr_x_xzzz_xx,   \
                             tr_x_xzzz_xxx,  \
                             tr_x_xzzz_xxy,  \
                             tr_x_xzzz_xxz,  \
                             tr_x_xzzz_xy,   \
                             tr_x_xzzz_xyy,  \
                             tr_x_xzzz_xyz,  \
                             tr_x_xzzz_xz,   \
                             tr_x_xzzz_xzz,  \
                             tr_x_xzzz_zzz,  \
                             tr_x_yzzz_yyy,  \
                             tr_x_yzzz_yyz,  \
                             tr_x_yzzz_yzz,  \
                             ts_yzzz_yyy,    \
                             ts_yzzz_yyz,    \
                             ts_yzzz_yzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyzzz_xxx[i] = tr_x_xzzz_xxx[i] * pa_y[i];

        tr_x_xyzzz_xxy[i] = tr_x_xzzz_xx[i] * fe_0 + tr_x_xzzz_xxy[i] * pa_y[i];

        tr_x_xyzzz_xxz[i] = tr_x_xzzz_xxz[i] * pa_y[i];

        tr_x_xyzzz_xyy[i] = 2.0 * tr_x_xzzz_xy[i] * fe_0 + tr_x_xzzz_xyy[i] * pa_y[i];

        tr_x_xyzzz_xyz[i] = tr_x_xzzz_xz[i] * fe_0 + tr_x_xzzz_xyz[i] * pa_y[i];

        tr_x_xyzzz_xzz[i] = tr_x_xzzz_xzz[i] * pa_y[i];

        tr_x_xyzzz_yyy[i] = ts_yzzz_yyy[i] * fe_0 + tr_x_yzzz_yyy[i] * pa_x[i];

        tr_x_xyzzz_yyz[i] = ts_yzzz_yyz[i] * fe_0 + tr_x_yzzz_yyz[i] * pa_x[i];

        tr_x_xyzzz_yzz[i] = ts_yzzz_yzz[i] * fe_0 + tr_x_yzzz_yzz[i] * pa_x[i];

        tr_x_xyzzz_zzz[i] = tr_x_xzzz_zzz[i] * pa_y[i];
    }

    // Set up 140-150 components of targeted buffer : HF

    auto tr_x_xzzzz_xxx = pbuffer.data(idx_dip_hf + 140);

    auto tr_x_xzzzz_xxy = pbuffer.data(idx_dip_hf + 141);

    auto tr_x_xzzzz_xxz = pbuffer.data(idx_dip_hf + 142);

    auto tr_x_xzzzz_xyy = pbuffer.data(idx_dip_hf + 143);

    auto tr_x_xzzzz_xyz = pbuffer.data(idx_dip_hf + 144);

    auto tr_x_xzzzz_xzz = pbuffer.data(idx_dip_hf + 145);

    auto tr_x_xzzzz_yyy = pbuffer.data(idx_dip_hf + 146);

    auto tr_x_xzzzz_yyz = pbuffer.data(idx_dip_hf + 147);

    auto tr_x_xzzzz_yzz = pbuffer.data(idx_dip_hf + 148);

    auto tr_x_xzzzz_zzz = pbuffer.data(idx_dip_hf + 149);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_x_xzz_xxx,   \
                             tr_x_xzz_xxy,   \
                             tr_x_xzz_xyy,   \
                             tr_x_xzzz_xxx,  \
                             tr_x_xzzz_xxy,  \
                             tr_x_xzzz_xyy,  \
                             tr_x_xzzzz_xxx, \
                             tr_x_xzzzz_xxy, \
                             tr_x_xzzzz_xxz, \
                             tr_x_xzzzz_xyy, \
                             tr_x_xzzzz_xyz, \
                             tr_x_xzzzz_xzz, \
                             tr_x_xzzzz_yyy, \
                             tr_x_xzzzz_yyz, \
                             tr_x_xzzzz_yzz, \
                             tr_x_xzzzz_zzz, \
                             tr_x_zzzz_xxz,  \
                             tr_x_zzzz_xyz,  \
                             tr_x_zzzz_xz,   \
                             tr_x_zzzz_xzz,  \
                             tr_x_zzzz_yyy,  \
                             tr_x_zzzz_yyz,  \
                             tr_x_zzzz_yz,   \
                             tr_x_zzzz_yzz,  \
                             tr_x_zzzz_zz,   \
                             tr_x_zzzz_zzz,  \
                             ts_zzzz_xxz,    \
                             ts_zzzz_xyz,    \
                             ts_zzzz_xzz,    \
                             ts_zzzz_yyy,    \
                             ts_zzzz_yyz,    \
                             ts_zzzz_yzz,    \
                             ts_zzzz_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzzzz_xxx[i] = 3.0 * tr_x_xzz_xxx[i] * fe_0 + tr_x_xzzz_xxx[i] * pa_z[i];

        tr_x_xzzzz_xxy[i] = 3.0 * tr_x_xzz_xxy[i] * fe_0 + tr_x_xzzz_xxy[i] * pa_z[i];

        tr_x_xzzzz_xxz[i] = 2.0 * tr_x_zzzz_xz[i] * fe_0 + ts_zzzz_xxz[i] * fe_0 + tr_x_zzzz_xxz[i] * pa_x[i];

        tr_x_xzzzz_xyy[i] = 3.0 * tr_x_xzz_xyy[i] * fe_0 + tr_x_xzzz_xyy[i] * pa_z[i];

        tr_x_xzzzz_xyz[i] = tr_x_zzzz_yz[i] * fe_0 + ts_zzzz_xyz[i] * fe_0 + tr_x_zzzz_xyz[i] * pa_x[i];

        tr_x_xzzzz_xzz[i] = tr_x_zzzz_zz[i] * fe_0 + ts_zzzz_xzz[i] * fe_0 + tr_x_zzzz_xzz[i] * pa_x[i];

        tr_x_xzzzz_yyy[i] = ts_zzzz_yyy[i] * fe_0 + tr_x_zzzz_yyy[i] * pa_x[i];

        tr_x_xzzzz_yyz[i] = ts_zzzz_yyz[i] * fe_0 + tr_x_zzzz_yyz[i] * pa_x[i];

        tr_x_xzzzz_yzz[i] = ts_zzzz_yzz[i] * fe_0 + tr_x_zzzz_yzz[i] * pa_x[i];

        tr_x_xzzzz_zzz[i] = ts_zzzz_zzz[i] * fe_0 + tr_x_zzzz_zzz[i] * pa_x[i];
    }

    // Set up 150-160 components of targeted buffer : HF

    auto tr_x_yyyyy_xxx = pbuffer.data(idx_dip_hf + 150);

    auto tr_x_yyyyy_xxy = pbuffer.data(idx_dip_hf + 151);

    auto tr_x_yyyyy_xxz = pbuffer.data(idx_dip_hf + 152);

    auto tr_x_yyyyy_xyy = pbuffer.data(idx_dip_hf + 153);

    auto tr_x_yyyyy_xyz = pbuffer.data(idx_dip_hf + 154);

    auto tr_x_yyyyy_xzz = pbuffer.data(idx_dip_hf + 155);

    auto tr_x_yyyyy_yyy = pbuffer.data(idx_dip_hf + 156);

    auto tr_x_yyyyy_yyz = pbuffer.data(idx_dip_hf + 157);

    auto tr_x_yyyyy_yzz = pbuffer.data(idx_dip_hf + 158);

    auto tr_x_yyyyy_zzz = pbuffer.data(idx_dip_hf + 159);

#pragma omp simd aligned(pa_y,               \
                             tr_x_yyy_xxx,   \
                             tr_x_yyy_xxy,   \
                             tr_x_yyy_xxz,   \
                             tr_x_yyy_xyy,   \
                             tr_x_yyy_xyz,   \
                             tr_x_yyy_xzz,   \
                             tr_x_yyy_yyy,   \
                             tr_x_yyy_yyz,   \
                             tr_x_yyy_yzz,   \
                             tr_x_yyy_zzz,   \
                             tr_x_yyyy_xx,   \
                             tr_x_yyyy_xxx,  \
                             tr_x_yyyy_xxy,  \
                             tr_x_yyyy_xxz,  \
                             tr_x_yyyy_xy,   \
                             tr_x_yyyy_xyy,  \
                             tr_x_yyyy_xyz,  \
                             tr_x_yyyy_xz,   \
                             tr_x_yyyy_xzz,  \
                             tr_x_yyyy_yy,   \
                             tr_x_yyyy_yyy,  \
                             tr_x_yyyy_yyz,  \
                             tr_x_yyyy_yz,   \
                             tr_x_yyyy_yzz,  \
                             tr_x_yyyy_zz,   \
                             tr_x_yyyy_zzz,  \
                             tr_x_yyyyy_xxx, \
                             tr_x_yyyyy_xxy, \
                             tr_x_yyyyy_xxz, \
                             tr_x_yyyyy_xyy, \
                             tr_x_yyyyy_xyz, \
                             tr_x_yyyyy_xzz, \
                             tr_x_yyyyy_yyy, \
                             tr_x_yyyyy_yyz, \
                             tr_x_yyyyy_yzz, \
                             tr_x_yyyyy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyy_xxx[i] = 4.0 * tr_x_yyy_xxx[i] * fe_0 + tr_x_yyyy_xxx[i] * pa_y[i];

        tr_x_yyyyy_xxy[i] = 4.0 * tr_x_yyy_xxy[i] * fe_0 + tr_x_yyyy_xx[i] * fe_0 + tr_x_yyyy_xxy[i] * pa_y[i];

        tr_x_yyyyy_xxz[i] = 4.0 * tr_x_yyy_xxz[i] * fe_0 + tr_x_yyyy_xxz[i] * pa_y[i];

        tr_x_yyyyy_xyy[i] = 4.0 * tr_x_yyy_xyy[i] * fe_0 + 2.0 * tr_x_yyyy_xy[i] * fe_0 + tr_x_yyyy_xyy[i] * pa_y[i];

        tr_x_yyyyy_xyz[i] = 4.0 * tr_x_yyy_xyz[i] * fe_0 + tr_x_yyyy_xz[i] * fe_0 + tr_x_yyyy_xyz[i] * pa_y[i];

        tr_x_yyyyy_xzz[i] = 4.0 * tr_x_yyy_xzz[i] * fe_0 + tr_x_yyyy_xzz[i] * pa_y[i];

        tr_x_yyyyy_yyy[i] = 4.0 * tr_x_yyy_yyy[i] * fe_0 + 3.0 * tr_x_yyyy_yy[i] * fe_0 + tr_x_yyyy_yyy[i] * pa_y[i];

        tr_x_yyyyy_yyz[i] = 4.0 * tr_x_yyy_yyz[i] * fe_0 + 2.0 * tr_x_yyyy_yz[i] * fe_0 + tr_x_yyyy_yyz[i] * pa_y[i];

        tr_x_yyyyy_yzz[i] = 4.0 * tr_x_yyy_yzz[i] * fe_0 + tr_x_yyyy_zz[i] * fe_0 + tr_x_yyyy_yzz[i] * pa_y[i];

        tr_x_yyyyy_zzz[i] = 4.0 * tr_x_yyy_zzz[i] * fe_0 + tr_x_yyyy_zzz[i] * pa_y[i];
    }

    // Set up 160-170 components of targeted buffer : HF

    auto tr_x_yyyyz_xxx = pbuffer.data(idx_dip_hf + 160);

    auto tr_x_yyyyz_xxy = pbuffer.data(idx_dip_hf + 161);

    auto tr_x_yyyyz_xxz = pbuffer.data(idx_dip_hf + 162);

    auto tr_x_yyyyz_xyy = pbuffer.data(idx_dip_hf + 163);

    auto tr_x_yyyyz_xyz = pbuffer.data(idx_dip_hf + 164);

    auto tr_x_yyyyz_xzz = pbuffer.data(idx_dip_hf + 165);

    auto tr_x_yyyyz_yyy = pbuffer.data(idx_dip_hf + 166);

    auto tr_x_yyyyz_yyz = pbuffer.data(idx_dip_hf + 167);

    auto tr_x_yyyyz_yzz = pbuffer.data(idx_dip_hf + 168);

    auto tr_x_yyyyz_zzz = pbuffer.data(idx_dip_hf + 169);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_yyyy_xxx,  \
                             tr_x_yyyy_xxy,  \
                             tr_x_yyyy_xy,   \
                             tr_x_yyyy_xyy,  \
                             tr_x_yyyy_xyz,  \
                             tr_x_yyyy_yy,   \
                             tr_x_yyyy_yyy,  \
                             tr_x_yyyy_yyz,  \
                             tr_x_yyyy_yz,   \
                             tr_x_yyyy_yzz,  \
                             tr_x_yyyyz_xxx, \
                             tr_x_yyyyz_xxy, \
                             tr_x_yyyyz_xxz, \
                             tr_x_yyyyz_xyy, \
                             tr_x_yyyyz_xyz, \
                             tr_x_yyyyz_xzz, \
                             tr_x_yyyyz_yyy, \
                             tr_x_yyyyz_yyz, \
                             tr_x_yyyyz_yzz, \
                             tr_x_yyyyz_zzz, \
                             tr_x_yyyz_xxz,  \
                             tr_x_yyyz_xzz,  \
                             tr_x_yyyz_zzz,  \
                             tr_x_yyz_xxz,   \
                             tr_x_yyz_xzz,   \
                             tr_x_yyz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyyz_xxx[i] = tr_x_yyyy_xxx[i] * pa_z[i];

        tr_x_yyyyz_xxy[i] = tr_x_yyyy_xxy[i] * pa_z[i];

        tr_x_yyyyz_xxz[i] = 3.0 * tr_x_yyz_xxz[i] * fe_0 + tr_x_yyyz_xxz[i] * pa_y[i];

        tr_x_yyyyz_xyy[i] = tr_x_yyyy_xyy[i] * pa_z[i];

        tr_x_yyyyz_xyz[i] = tr_x_yyyy_xy[i] * fe_0 + tr_x_yyyy_xyz[i] * pa_z[i];

        tr_x_yyyyz_xzz[i] = 3.0 * tr_x_yyz_xzz[i] * fe_0 + tr_x_yyyz_xzz[i] * pa_y[i];

        tr_x_yyyyz_yyy[i] = tr_x_yyyy_yyy[i] * pa_z[i];

        tr_x_yyyyz_yyz[i] = tr_x_yyyy_yy[i] * fe_0 + tr_x_yyyy_yyz[i] * pa_z[i];

        tr_x_yyyyz_yzz[i] = 2.0 * tr_x_yyyy_yz[i] * fe_0 + tr_x_yyyy_yzz[i] * pa_z[i];

        tr_x_yyyyz_zzz[i] = 3.0 * tr_x_yyz_zzz[i] * fe_0 + tr_x_yyyz_zzz[i] * pa_y[i];
    }

    // Set up 170-180 components of targeted buffer : HF

    auto tr_x_yyyzz_xxx = pbuffer.data(idx_dip_hf + 170);

    auto tr_x_yyyzz_xxy = pbuffer.data(idx_dip_hf + 171);

    auto tr_x_yyyzz_xxz = pbuffer.data(idx_dip_hf + 172);

    auto tr_x_yyyzz_xyy = pbuffer.data(idx_dip_hf + 173);

    auto tr_x_yyyzz_xyz = pbuffer.data(idx_dip_hf + 174);

    auto tr_x_yyyzz_xzz = pbuffer.data(idx_dip_hf + 175);

    auto tr_x_yyyzz_yyy = pbuffer.data(idx_dip_hf + 176);

    auto tr_x_yyyzz_yyz = pbuffer.data(idx_dip_hf + 177);

    auto tr_x_yyyzz_yzz = pbuffer.data(idx_dip_hf + 178);

    auto tr_x_yyyzz_zzz = pbuffer.data(idx_dip_hf + 179);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_yyy_xxy,   \
                             tr_x_yyy_xyy,   \
                             tr_x_yyy_yyy,   \
                             tr_x_yyyz_xxy,  \
                             tr_x_yyyz_xyy,  \
                             tr_x_yyyz_yyy,  \
                             tr_x_yyyzz_xxx, \
                             tr_x_yyyzz_xxy, \
                             tr_x_yyyzz_xxz, \
                             tr_x_yyyzz_xyy, \
                             tr_x_yyyzz_xyz, \
                             tr_x_yyyzz_xzz, \
                             tr_x_yyyzz_yyy, \
                             tr_x_yyyzz_yyz, \
                             tr_x_yyyzz_yzz, \
                             tr_x_yyyzz_zzz, \
                             tr_x_yyzz_xxx,  \
                             tr_x_yyzz_xxz,  \
                             tr_x_yyzz_xyz,  \
                             tr_x_yyzz_xz,   \
                             tr_x_yyzz_xzz,  \
                             tr_x_yyzz_yyz,  \
                             tr_x_yyzz_yz,   \
                             tr_x_yyzz_yzz,  \
                             tr_x_yyzz_zz,   \
                             tr_x_yyzz_zzz,  \
                             tr_x_yzz_xxx,   \
                             tr_x_yzz_xxz,   \
                             tr_x_yzz_xyz,   \
                             tr_x_yzz_xzz,   \
                             tr_x_yzz_yyz,   \
                             tr_x_yzz_yzz,   \
                             tr_x_yzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyyzz_xxx[i] = 2.0 * tr_x_yzz_xxx[i] * fe_0 + tr_x_yyzz_xxx[i] * pa_y[i];

        tr_x_yyyzz_xxy[i] = tr_x_yyy_xxy[i] * fe_0 + tr_x_yyyz_xxy[i] * pa_z[i];

        tr_x_yyyzz_xxz[i] = 2.0 * tr_x_yzz_xxz[i] * fe_0 + tr_x_yyzz_xxz[i] * pa_y[i];

        tr_x_yyyzz_xyy[i] = tr_x_yyy_xyy[i] * fe_0 + tr_x_yyyz_xyy[i] * pa_z[i];

        tr_x_yyyzz_xyz[i] = 2.0 * tr_x_yzz_xyz[i] * fe_0 + tr_x_yyzz_xz[i] * fe_0 + tr_x_yyzz_xyz[i] * pa_y[i];

        tr_x_yyyzz_xzz[i] = 2.0 * tr_x_yzz_xzz[i] * fe_0 + tr_x_yyzz_xzz[i] * pa_y[i];

        tr_x_yyyzz_yyy[i] = tr_x_yyy_yyy[i] * fe_0 + tr_x_yyyz_yyy[i] * pa_z[i];

        tr_x_yyyzz_yyz[i] = 2.0 * tr_x_yzz_yyz[i] * fe_0 + 2.0 * tr_x_yyzz_yz[i] * fe_0 + tr_x_yyzz_yyz[i] * pa_y[i];

        tr_x_yyyzz_yzz[i] = 2.0 * tr_x_yzz_yzz[i] * fe_0 + tr_x_yyzz_zz[i] * fe_0 + tr_x_yyzz_yzz[i] * pa_y[i];

        tr_x_yyyzz_zzz[i] = 2.0 * tr_x_yzz_zzz[i] * fe_0 + tr_x_yyzz_zzz[i] * pa_y[i];
    }

    // Set up 180-190 components of targeted buffer : HF

    auto tr_x_yyzzz_xxx = pbuffer.data(idx_dip_hf + 180);

    auto tr_x_yyzzz_xxy = pbuffer.data(idx_dip_hf + 181);

    auto tr_x_yyzzz_xxz = pbuffer.data(idx_dip_hf + 182);

    auto tr_x_yyzzz_xyy = pbuffer.data(idx_dip_hf + 183);

    auto tr_x_yyzzz_xyz = pbuffer.data(idx_dip_hf + 184);

    auto tr_x_yyzzz_xzz = pbuffer.data(idx_dip_hf + 185);

    auto tr_x_yyzzz_yyy = pbuffer.data(idx_dip_hf + 186);

    auto tr_x_yyzzz_yyz = pbuffer.data(idx_dip_hf + 187);

    auto tr_x_yyzzz_yzz = pbuffer.data(idx_dip_hf + 188);

    auto tr_x_yyzzz_zzz = pbuffer.data(idx_dip_hf + 189);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_yyz_xxy,   \
                             tr_x_yyz_xyy,   \
                             tr_x_yyz_yyy,   \
                             tr_x_yyzz_xxy,  \
                             tr_x_yyzz_xyy,  \
                             tr_x_yyzz_yyy,  \
                             tr_x_yyzzz_xxx, \
                             tr_x_yyzzz_xxy, \
                             tr_x_yyzzz_xxz, \
                             tr_x_yyzzz_xyy, \
                             tr_x_yyzzz_xyz, \
                             tr_x_yyzzz_xzz, \
                             tr_x_yyzzz_yyy, \
                             tr_x_yyzzz_yyz, \
                             tr_x_yyzzz_yzz, \
                             tr_x_yyzzz_zzz, \
                             tr_x_yzzz_xxx,  \
                             tr_x_yzzz_xxz,  \
                             tr_x_yzzz_xyz,  \
                             tr_x_yzzz_xz,   \
                             tr_x_yzzz_xzz,  \
                             tr_x_yzzz_yyz,  \
                             tr_x_yzzz_yz,   \
                             tr_x_yzzz_yzz,  \
                             tr_x_yzzz_zz,   \
                             tr_x_yzzz_zzz,  \
                             tr_x_zzz_xxx,   \
                             tr_x_zzz_xxz,   \
                             tr_x_zzz_xyz,   \
                             tr_x_zzz_xzz,   \
                             tr_x_zzz_yyz,   \
                             tr_x_zzz_yzz,   \
                             tr_x_zzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyzzz_xxx[i] = tr_x_zzz_xxx[i] * fe_0 + tr_x_yzzz_xxx[i] * pa_y[i];

        tr_x_yyzzz_xxy[i] = 2.0 * tr_x_yyz_xxy[i] * fe_0 + tr_x_yyzz_xxy[i] * pa_z[i];

        tr_x_yyzzz_xxz[i] = tr_x_zzz_xxz[i] * fe_0 + tr_x_yzzz_xxz[i] * pa_y[i];

        tr_x_yyzzz_xyy[i] = 2.0 * tr_x_yyz_xyy[i] * fe_0 + tr_x_yyzz_xyy[i] * pa_z[i];

        tr_x_yyzzz_xyz[i] = tr_x_zzz_xyz[i] * fe_0 + tr_x_yzzz_xz[i] * fe_0 + tr_x_yzzz_xyz[i] * pa_y[i];

        tr_x_yyzzz_xzz[i] = tr_x_zzz_xzz[i] * fe_0 + tr_x_yzzz_xzz[i] * pa_y[i];

        tr_x_yyzzz_yyy[i] = 2.0 * tr_x_yyz_yyy[i] * fe_0 + tr_x_yyzz_yyy[i] * pa_z[i];

        tr_x_yyzzz_yyz[i] = tr_x_zzz_yyz[i] * fe_0 + 2.0 * tr_x_yzzz_yz[i] * fe_0 + tr_x_yzzz_yyz[i] * pa_y[i];

        tr_x_yyzzz_yzz[i] = tr_x_zzz_yzz[i] * fe_0 + tr_x_yzzz_zz[i] * fe_0 + tr_x_yzzz_yzz[i] * pa_y[i];

        tr_x_yyzzz_zzz[i] = tr_x_zzz_zzz[i] * fe_0 + tr_x_yzzz_zzz[i] * pa_y[i];
    }

    // Set up 190-200 components of targeted buffer : HF

    auto tr_x_yzzzz_xxx = pbuffer.data(idx_dip_hf + 190);

    auto tr_x_yzzzz_xxy = pbuffer.data(idx_dip_hf + 191);

    auto tr_x_yzzzz_xxz = pbuffer.data(idx_dip_hf + 192);

    auto tr_x_yzzzz_xyy = pbuffer.data(idx_dip_hf + 193);

    auto tr_x_yzzzz_xyz = pbuffer.data(idx_dip_hf + 194);

    auto tr_x_yzzzz_xzz = pbuffer.data(idx_dip_hf + 195);

    auto tr_x_yzzzz_yyy = pbuffer.data(idx_dip_hf + 196);

    auto tr_x_yzzzz_yyz = pbuffer.data(idx_dip_hf + 197);

    auto tr_x_yzzzz_yzz = pbuffer.data(idx_dip_hf + 198);

    auto tr_x_yzzzz_zzz = pbuffer.data(idx_dip_hf + 199);

#pragma omp simd aligned(pa_y,               \
                             tr_x_yzzzz_xxx, \
                             tr_x_yzzzz_xxy, \
                             tr_x_yzzzz_xxz, \
                             tr_x_yzzzz_xyy, \
                             tr_x_yzzzz_xyz, \
                             tr_x_yzzzz_xzz, \
                             tr_x_yzzzz_yyy, \
                             tr_x_yzzzz_yyz, \
                             tr_x_yzzzz_yzz, \
                             tr_x_yzzzz_zzz, \
                             tr_x_zzzz_xx,   \
                             tr_x_zzzz_xxx,  \
                             tr_x_zzzz_xxy,  \
                             tr_x_zzzz_xxz,  \
                             tr_x_zzzz_xy,   \
                             tr_x_zzzz_xyy,  \
                             tr_x_zzzz_xyz,  \
                             tr_x_zzzz_xz,   \
                             tr_x_zzzz_xzz,  \
                             tr_x_zzzz_yy,   \
                             tr_x_zzzz_yyy,  \
                             tr_x_zzzz_yyz,  \
                             tr_x_zzzz_yz,   \
                             tr_x_zzzz_yzz,  \
                             tr_x_zzzz_zz,   \
                             tr_x_zzzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzzzz_xxx[i] = tr_x_zzzz_xxx[i] * pa_y[i];

        tr_x_yzzzz_xxy[i] = tr_x_zzzz_xx[i] * fe_0 + tr_x_zzzz_xxy[i] * pa_y[i];

        tr_x_yzzzz_xxz[i] = tr_x_zzzz_xxz[i] * pa_y[i];

        tr_x_yzzzz_xyy[i] = 2.0 * tr_x_zzzz_xy[i] * fe_0 + tr_x_zzzz_xyy[i] * pa_y[i];

        tr_x_yzzzz_xyz[i] = tr_x_zzzz_xz[i] * fe_0 + tr_x_zzzz_xyz[i] * pa_y[i];

        tr_x_yzzzz_xzz[i] = tr_x_zzzz_xzz[i] * pa_y[i];

        tr_x_yzzzz_yyy[i] = 3.0 * tr_x_zzzz_yy[i] * fe_0 + tr_x_zzzz_yyy[i] * pa_y[i];

        tr_x_yzzzz_yyz[i] = 2.0 * tr_x_zzzz_yz[i] * fe_0 + tr_x_zzzz_yyz[i] * pa_y[i];

        tr_x_yzzzz_yzz[i] = tr_x_zzzz_zz[i] * fe_0 + tr_x_zzzz_yzz[i] * pa_y[i];

        tr_x_yzzzz_zzz[i] = tr_x_zzzz_zzz[i] * pa_y[i];
    }

    // Set up 200-210 components of targeted buffer : HF

    auto tr_x_zzzzz_xxx = pbuffer.data(idx_dip_hf + 200);

    auto tr_x_zzzzz_xxy = pbuffer.data(idx_dip_hf + 201);

    auto tr_x_zzzzz_xxz = pbuffer.data(idx_dip_hf + 202);

    auto tr_x_zzzzz_xyy = pbuffer.data(idx_dip_hf + 203);

    auto tr_x_zzzzz_xyz = pbuffer.data(idx_dip_hf + 204);

    auto tr_x_zzzzz_xzz = pbuffer.data(idx_dip_hf + 205);

    auto tr_x_zzzzz_yyy = pbuffer.data(idx_dip_hf + 206);

    auto tr_x_zzzzz_yyz = pbuffer.data(idx_dip_hf + 207);

    auto tr_x_zzzzz_yzz = pbuffer.data(idx_dip_hf + 208);

    auto tr_x_zzzzz_zzz = pbuffer.data(idx_dip_hf + 209);

#pragma omp simd aligned(pa_z,               \
                             tr_x_zzz_xxx,   \
                             tr_x_zzz_xxy,   \
                             tr_x_zzz_xxz,   \
                             tr_x_zzz_xyy,   \
                             tr_x_zzz_xyz,   \
                             tr_x_zzz_xzz,   \
                             tr_x_zzz_yyy,   \
                             tr_x_zzz_yyz,   \
                             tr_x_zzz_yzz,   \
                             tr_x_zzz_zzz,   \
                             tr_x_zzzz_xx,   \
                             tr_x_zzzz_xxx,  \
                             tr_x_zzzz_xxy,  \
                             tr_x_zzzz_xxz,  \
                             tr_x_zzzz_xy,   \
                             tr_x_zzzz_xyy,  \
                             tr_x_zzzz_xyz,  \
                             tr_x_zzzz_xz,   \
                             tr_x_zzzz_xzz,  \
                             tr_x_zzzz_yy,   \
                             tr_x_zzzz_yyy,  \
                             tr_x_zzzz_yyz,  \
                             tr_x_zzzz_yz,   \
                             tr_x_zzzz_yzz,  \
                             tr_x_zzzz_zz,   \
                             tr_x_zzzz_zzz,  \
                             tr_x_zzzzz_xxx, \
                             tr_x_zzzzz_xxy, \
                             tr_x_zzzzz_xxz, \
                             tr_x_zzzzz_xyy, \
                             tr_x_zzzzz_xyz, \
                             tr_x_zzzzz_xzz, \
                             tr_x_zzzzz_yyy, \
                             tr_x_zzzzz_yyz, \
                             tr_x_zzzzz_yzz, \
                             tr_x_zzzzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzzzz_xxx[i] = 4.0 * tr_x_zzz_xxx[i] * fe_0 + tr_x_zzzz_xxx[i] * pa_z[i];

        tr_x_zzzzz_xxy[i] = 4.0 * tr_x_zzz_xxy[i] * fe_0 + tr_x_zzzz_xxy[i] * pa_z[i];

        tr_x_zzzzz_xxz[i] = 4.0 * tr_x_zzz_xxz[i] * fe_0 + tr_x_zzzz_xx[i] * fe_0 + tr_x_zzzz_xxz[i] * pa_z[i];

        tr_x_zzzzz_xyy[i] = 4.0 * tr_x_zzz_xyy[i] * fe_0 + tr_x_zzzz_xyy[i] * pa_z[i];

        tr_x_zzzzz_xyz[i] = 4.0 * tr_x_zzz_xyz[i] * fe_0 + tr_x_zzzz_xy[i] * fe_0 + tr_x_zzzz_xyz[i] * pa_z[i];

        tr_x_zzzzz_xzz[i] = 4.0 * tr_x_zzz_xzz[i] * fe_0 + 2.0 * tr_x_zzzz_xz[i] * fe_0 + tr_x_zzzz_xzz[i] * pa_z[i];

        tr_x_zzzzz_yyy[i] = 4.0 * tr_x_zzz_yyy[i] * fe_0 + tr_x_zzzz_yyy[i] * pa_z[i];

        tr_x_zzzzz_yyz[i] = 4.0 * tr_x_zzz_yyz[i] * fe_0 + tr_x_zzzz_yy[i] * fe_0 + tr_x_zzzz_yyz[i] * pa_z[i];

        tr_x_zzzzz_yzz[i] = 4.0 * tr_x_zzz_yzz[i] * fe_0 + 2.0 * tr_x_zzzz_yz[i] * fe_0 + tr_x_zzzz_yzz[i] * pa_z[i];

        tr_x_zzzzz_zzz[i] = 4.0 * tr_x_zzz_zzz[i] * fe_0 + 3.0 * tr_x_zzzz_zz[i] * fe_0 + tr_x_zzzz_zzz[i] * pa_z[i];
    }

    // Set up 210-220 components of targeted buffer : HF

    auto tr_y_xxxxx_xxx = pbuffer.data(idx_dip_hf + 210);

    auto tr_y_xxxxx_xxy = pbuffer.data(idx_dip_hf + 211);

    auto tr_y_xxxxx_xxz = pbuffer.data(idx_dip_hf + 212);

    auto tr_y_xxxxx_xyy = pbuffer.data(idx_dip_hf + 213);

    auto tr_y_xxxxx_xyz = pbuffer.data(idx_dip_hf + 214);

    auto tr_y_xxxxx_xzz = pbuffer.data(idx_dip_hf + 215);

    auto tr_y_xxxxx_yyy = pbuffer.data(idx_dip_hf + 216);

    auto tr_y_xxxxx_yyz = pbuffer.data(idx_dip_hf + 217);

    auto tr_y_xxxxx_yzz = pbuffer.data(idx_dip_hf + 218);

    auto tr_y_xxxxx_zzz = pbuffer.data(idx_dip_hf + 219);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xxx_xxx,   \
                             tr_y_xxx_xxy,   \
                             tr_y_xxx_xxz,   \
                             tr_y_xxx_xyy,   \
                             tr_y_xxx_xyz,   \
                             tr_y_xxx_xzz,   \
                             tr_y_xxx_yyy,   \
                             tr_y_xxx_yyz,   \
                             tr_y_xxx_yzz,   \
                             tr_y_xxx_zzz,   \
                             tr_y_xxxx_xx,   \
                             tr_y_xxxx_xxx,  \
                             tr_y_xxxx_xxy,  \
                             tr_y_xxxx_xxz,  \
                             tr_y_xxxx_xy,   \
                             tr_y_xxxx_xyy,  \
                             tr_y_xxxx_xyz,  \
                             tr_y_xxxx_xz,   \
                             tr_y_xxxx_xzz,  \
                             tr_y_xxxx_yy,   \
                             tr_y_xxxx_yyy,  \
                             tr_y_xxxx_yyz,  \
                             tr_y_xxxx_yz,   \
                             tr_y_xxxx_yzz,  \
                             tr_y_xxxx_zz,   \
                             tr_y_xxxx_zzz,  \
                             tr_y_xxxxx_xxx, \
                             tr_y_xxxxx_xxy, \
                             tr_y_xxxxx_xxz, \
                             tr_y_xxxxx_xyy, \
                             tr_y_xxxxx_xyz, \
                             tr_y_xxxxx_xzz, \
                             tr_y_xxxxx_yyy, \
                             tr_y_xxxxx_yyz, \
                             tr_y_xxxxx_yzz, \
                             tr_y_xxxxx_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxx_xxx[i] = 4.0 * tr_y_xxx_xxx[i] * fe_0 + 3.0 * tr_y_xxxx_xx[i] * fe_0 + tr_y_xxxx_xxx[i] * pa_x[i];

        tr_y_xxxxx_xxy[i] = 4.0 * tr_y_xxx_xxy[i] * fe_0 + 2.0 * tr_y_xxxx_xy[i] * fe_0 + tr_y_xxxx_xxy[i] * pa_x[i];

        tr_y_xxxxx_xxz[i] = 4.0 * tr_y_xxx_xxz[i] * fe_0 + 2.0 * tr_y_xxxx_xz[i] * fe_0 + tr_y_xxxx_xxz[i] * pa_x[i];

        tr_y_xxxxx_xyy[i] = 4.0 * tr_y_xxx_xyy[i] * fe_0 + tr_y_xxxx_yy[i] * fe_0 + tr_y_xxxx_xyy[i] * pa_x[i];

        tr_y_xxxxx_xyz[i] = 4.0 * tr_y_xxx_xyz[i] * fe_0 + tr_y_xxxx_yz[i] * fe_0 + tr_y_xxxx_xyz[i] * pa_x[i];

        tr_y_xxxxx_xzz[i] = 4.0 * tr_y_xxx_xzz[i] * fe_0 + tr_y_xxxx_zz[i] * fe_0 + tr_y_xxxx_xzz[i] * pa_x[i];

        tr_y_xxxxx_yyy[i] = 4.0 * tr_y_xxx_yyy[i] * fe_0 + tr_y_xxxx_yyy[i] * pa_x[i];

        tr_y_xxxxx_yyz[i] = 4.0 * tr_y_xxx_yyz[i] * fe_0 + tr_y_xxxx_yyz[i] * pa_x[i];

        tr_y_xxxxx_yzz[i] = 4.0 * tr_y_xxx_yzz[i] * fe_0 + tr_y_xxxx_yzz[i] * pa_x[i];

        tr_y_xxxxx_zzz[i] = 4.0 * tr_y_xxx_zzz[i] * fe_0 + tr_y_xxxx_zzz[i] * pa_x[i];
    }

    // Set up 220-230 components of targeted buffer : HF

    auto tr_y_xxxxy_xxx = pbuffer.data(idx_dip_hf + 220);

    auto tr_y_xxxxy_xxy = pbuffer.data(idx_dip_hf + 221);

    auto tr_y_xxxxy_xxz = pbuffer.data(idx_dip_hf + 222);

    auto tr_y_xxxxy_xyy = pbuffer.data(idx_dip_hf + 223);

    auto tr_y_xxxxy_xyz = pbuffer.data(idx_dip_hf + 224);

    auto tr_y_xxxxy_xzz = pbuffer.data(idx_dip_hf + 225);

    auto tr_y_xxxxy_yyy = pbuffer.data(idx_dip_hf + 226);

    auto tr_y_xxxxy_yyz = pbuffer.data(idx_dip_hf + 227);

    auto tr_y_xxxxy_yzz = pbuffer.data(idx_dip_hf + 228);

    auto tr_y_xxxxy_zzz = pbuffer.data(idx_dip_hf + 229);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_y_xxxx_xxx,  \
                             tr_y_xxxx_xxz,  \
                             tr_y_xxxx_xzz,  \
                             tr_y_xxxxy_xxx, \
                             tr_y_xxxxy_xxy, \
                             tr_y_xxxxy_xxz, \
                             tr_y_xxxxy_xyy, \
                             tr_y_xxxxy_xyz, \
                             tr_y_xxxxy_xzz, \
                             tr_y_xxxxy_yyy, \
                             tr_y_xxxxy_yyz, \
                             tr_y_xxxxy_yzz, \
                             tr_y_xxxxy_zzz, \
                             tr_y_xxxy_xxy,  \
                             tr_y_xxxy_xy,   \
                             tr_y_xxxy_xyy,  \
                             tr_y_xxxy_xyz,  \
                             tr_y_xxxy_yy,   \
                             tr_y_xxxy_yyy,  \
                             tr_y_xxxy_yyz,  \
                             tr_y_xxxy_yz,   \
                             tr_y_xxxy_yzz,  \
                             tr_y_xxxy_zzz,  \
                             tr_y_xxy_xxy,   \
                             tr_y_xxy_xyy,   \
                             tr_y_xxy_xyz,   \
                             tr_y_xxy_yyy,   \
                             tr_y_xxy_yyz,   \
                             tr_y_xxy_yzz,   \
                             tr_y_xxy_zzz,   \
                             ts_xxxx_xxx,    \
                             ts_xxxx_xxz,    \
                             ts_xxxx_xzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxy_xxx[i] = ts_xxxx_xxx[i] * fe_0 + tr_y_xxxx_xxx[i] * pa_y[i];

        tr_y_xxxxy_xxy[i] = 3.0 * tr_y_xxy_xxy[i] * fe_0 + 2.0 * tr_y_xxxy_xy[i] * fe_0 + tr_y_xxxy_xxy[i] * pa_x[i];

        tr_y_xxxxy_xxz[i] = ts_xxxx_xxz[i] * fe_0 + tr_y_xxxx_xxz[i] * pa_y[i];

        tr_y_xxxxy_xyy[i] = 3.0 * tr_y_xxy_xyy[i] * fe_0 + tr_y_xxxy_yy[i] * fe_0 + tr_y_xxxy_xyy[i] * pa_x[i];

        tr_y_xxxxy_xyz[i] = 3.0 * tr_y_xxy_xyz[i] * fe_0 + tr_y_xxxy_yz[i] * fe_0 + tr_y_xxxy_xyz[i] * pa_x[i];

        tr_y_xxxxy_xzz[i] = ts_xxxx_xzz[i] * fe_0 + tr_y_xxxx_xzz[i] * pa_y[i];

        tr_y_xxxxy_yyy[i] = 3.0 * tr_y_xxy_yyy[i] * fe_0 + tr_y_xxxy_yyy[i] * pa_x[i];

        tr_y_xxxxy_yyz[i] = 3.0 * tr_y_xxy_yyz[i] * fe_0 + tr_y_xxxy_yyz[i] * pa_x[i];

        tr_y_xxxxy_yzz[i] = 3.0 * tr_y_xxy_yzz[i] * fe_0 + tr_y_xxxy_yzz[i] * pa_x[i];

        tr_y_xxxxy_zzz[i] = 3.0 * tr_y_xxy_zzz[i] * fe_0 + tr_y_xxxy_zzz[i] * pa_x[i];
    }

    // Set up 230-240 components of targeted buffer : HF

    auto tr_y_xxxxz_xxx = pbuffer.data(idx_dip_hf + 230);

    auto tr_y_xxxxz_xxy = pbuffer.data(idx_dip_hf + 231);

    auto tr_y_xxxxz_xxz = pbuffer.data(idx_dip_hf + 232);

    auto tr_y_xxxxz_xyy = pbuffer.data(idx_dip_hf + 233);

    auto tr_y_xxxxz_xyz = pbuffer.data(idx_dip_hf + 234);

    auto tr_y_xxxxz_xzz = pbuffer.data(idx_dip_hf + 235);

    auto tr_y_xxxxz_yyy = pbuffer.data(idx_dip_hf + 236);

    auto tr_y_xxxxz_yyz = pbuffer.data(idx_dip_hf + 237);

    auto tr_y_xxxxz_yzz = pbuffer.data(idx_dip_hf + 238);

    auto tr_y_xxxxz_zzz = pbuffer.data(idx_dip_hf + 239);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xxxx_xx,   \
                             tr_y_xxxx_xxx,  \
                             tr_y_xxxx_xxy,  \
                             tr_y_xxxx_xxz,  \
                             tr_y_xxxx_xy,   \
                             tr_y_xxxx_xyy,  \
                             tr_y_xxxx_xyz,  \
                             tr_y_xxxx_xz,   \
                             tr_y_xxxx_xzz,  \
                             tr_y_xxxx_yyy,  \
                             tr_y_xxxxz_xxx, \
                             tr_y_xxxxz_xxy, \
                             tr_y_xxxxz_xxz, \
                             tr_y_xxxxz_xyy, \
                             tr_y_xxxxz_xyz, \
                             tr_y_xxxxz_xzz, \
                             tr_y_xxxxz_yyy, \
                             tr_y_xxxxz_yyz, \
                             tr_y_xxxxz_yzz, \
                             tr_y_xxxxz_zzz, \
                             tr_y_xxxz_yyz,  \
                             tr_y_xxxz_yzz,  \
                             tr_y_xxxz_zzz,  \
                             tr_y_xxz_yyz,   \
                             tr_y_xxz_yzz,   \
                             tr_y_xxz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxxz_xxx[i] = tr_y_xxxx_xxx[i] * pa_z[i];

        tr_y_xxxxz_xxy[i] = tr_y_xxxx_xxy[i] * pa_z[i];

        tr_y_xxxxz_xxz[i] = tr_y_xxxx_xx[i] * fe_0 + tr_y_xxxx_xxz[i] * pa_z[i];

        tr_y_xxxxz_xyy[i] = tr_y_xxxx_xyy[i] * pa_z[i];

        tr_y_xxxxz_xyz[i] = tr_y_xxxx_xy[i] * fe_0 + tr_y_xxxx_xyz[i] * pa_z[i];

        tr_y_xxxxz_xzz[i] = 2.0 * tr_y_xxxx_xz[i] * fe_0 + tr_y_xxxx_xzz[i] * pa_z[i];

        tr_y_xxxxz_yyy[i] = tr_y_xxxx_yyy[i] * pa_z[i];

        tr_y_xxxxz_yyz[i] = 3.0 * tr_y_xxz_yyz[i] * fe_0 + tr_y_xxxz_yyz[i] * pa_x[i];

        tr_y_xxxxz_yzz[i] = 3.0 * tr_y_xxz_yzz[i] * fe_0 + tr_y_xxxz_yzz[i] * pa_x[i];

        tr_y_xxxxz_zzz[i] = 3.0 * tr_y_xxz_zzz[i] * fe_0 + tr_y_xxxz_zzz[i] * pa_x[i];
    }

    // Set up 240-250 components of targeted buffer : HF

    auto tr_y_xxxyy_xxx = pbuffer.data(idx_dip_hf + 240);

    auto tr_y_xxxyy_xxy = pbuffer.data(idx_dip_hf + 241);

    auto tr_y_xxxyy_xxz = pbuffer.data(idx_dip_hf + 242);

    auto tr_y_xxxyy_xyy = pbuffer.data(idx_dip_hf + 243);

    auto tr_y_xxxyy_xyz = pbuffer.data(idx_dip_hf + 244);

    auto tr_y_xxxyy_xzz = pbuffer.data(idx_dip_hf + 245);

    auto tr_y_xxxyy_yyy = pbuffer.data(idx_dip_hf + 246);

    auto tr_y_xxxyy_yyz = pbuffer.data(idx_dip_hf + 247);

    auto tr_y_xxxyy_yzz = pbuffer.data(idx_dip_hf + 248);

    auto tr_y_xxxyy_zzz = pbuffer.data(idx_dip_hf + 249);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xxxyy_xxx, \
                             tr_y_xxxyy_xxy, \
                             tr_y_xxxyy_xxz, \
                             tr_y_xxxyy_xyy, \
                             tr_y_xxxyy_xyz, \
                             tr_y_xxxyy_xzz, \
                             tr_y_xxxyy_yyy, \
                             tr_y_xxxyy_yyz, \
                             tr_y_xxxyy_yzz, \
                             tr_y_xxxyy_zzz, \
                             tr_y_xxyy_xx,   \
                             tr_y_xxyy_xxx,  \
                             tr_y_xxyy_xxy,  \
                             tr_y_xxyy_xxz,  \
                             tr_y_xxyy_xy,   \
                             tr_y_xxyy_xyy,  \
                             tr_y_xxyy_xyz,  \
                             tr_y_xxyy_xz,   \
                             tr_y_xxyy_xzz,  \
                             tr_y_xxyy_yy,   \
                             tr_y_xxyy_yyy,  \
                             tr_y_xxyy_yyz,  \
                             tr_y_xxyy_yz,   \
                             tr_y_xxyy_yzz,  \
                             tr_y_xxyy_zz,   \
                             tr_y_xxyy_zzz,  \
                             tr_y_xyy_xxx,   \
                             tr_y_xyy_xxy,   \
                             tr_y_xyy_xxz,   \
                             tr_y_xyy_xyy,   \
                             tr_y_xyy_xyz,   \
                             tr_y_xyy_xzz,   \
                             tr_y_xyy_yyy,   \
                             tr_y_xyy_yyz,   \
                             tr_y_xyy_yzz,   \
                             tr_y_xyy_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyy_xxx[i] = 2.0 * tr_y_xyy_xxx[i] * fe_0 + 3.0 * tr_y_xxyy_xx[i] * fe_0 + tr_y_xxyy_xxx[i] * pa_x[i];

        tr_y_xxxyy_xxy[i] = 2.0 * tr_y_xyy_xxy[i] * fe_0 + 2.0 * tr_y_xxyy_xy[i] * fe_0 + tr_y_xxyy_xxy[i] * pa_x[i];

        tr_y_xxxyy_xxz[i] = 2.0 * tr_y_xyy_xxz[i] * fe_0 + 2.0 * tr_y_xxyy_xz[i] * fe_0 + tr_y_xxyy_xxz[i] * pa_x[i];

        tr_y_xxxyy_xyy[i] = 2.0 * tr_y_xyy_xyy[i] * fe_0 + tr_y_xxyy_yy[i] * fe_0 + tr_y_xxyy_xyy[i] * pa_x[i];

        tr_y_xxxyy_xyz[i] = 2.0 * tr_y_xyy_xyz[i] * fe_0 + tr_y_xxyy_yz[i] * fe_0 + tr_y_xxyy_xyz[i] * pa_x[i];

        tr_y_xxxyy_xzz[i] = 2.0 * tr_y_xyy_xzz[i] * fe_0 + tr_y_xxyy_zz[i] * fe_0 + tr_y_xxyy_xzz[i] * pa_x[i];

        tr_y_xxxyy_yyy[i] = 2.0 * tr_y_xyy_yyy[i] * fe_0 + tr_y_xxyy_yyy[i] * pa_x[i];

        tr_y_xxxyy_yyz[i] = 2.0 * tr_y_xyy_yyz[i] * fe_0 + tr_y_xxyy_yyz[i] * pa_x[i];

        tr_y_xxxyy_yzz[i] = 2.0 * tr_y_xyy_yzz[i] * fe_0 + tr_y_xxyy_yzz[i] * pa_x[i];

        tr_y_xxxyy_zzz[i] = 2.0 * tr_y_xyy_zzz[i] * fe_0 + tr_y_xxyy_zzz[i] * pa_x[i];
    }

    // Set up 250-260 components of targeted buffer : HF

    auto tr_y_xxxyz_xxx = pbuffer.data(idx_dip_hf + 250);

    auto tr_y_xxxyz_xxy = pbuffer.data(idx_dip_hf + 251);

    auto tr_y_xxxyz_xxz = pbuffer.data(idx_dip_hf + 252);

    auto tr_y_xxxyz_xyy = pbuffer.data(idx_dip_hf + 253);

    auto tr_y_xxxyz_xyz = pbuffer.data(idx_dip_hf + 254);

    auto tr_y_xxxyz_xzz = pbuffer.data(idx_dip_hf + 255);

    auto tr_y_xxxyz_yyy = pbuffer.data(idx_dip_hf + 256);

    auto tr_y_xxxyz_yyz = pbuffer.data(idx_dip_hf + 257);

    auto tr_y_xxxyz_yzz = pbuffer.data(idx_dip_hf + 258);

    auto tr_y_xxxyz_zzz = pbuffer.data(idx_dip_hf + 259);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_y_xxxy_xxx,  \
                             tr_y_xxxy_xxy,  \
                             tr_y_xxxy_xy,   \
                             tr_y_xxxy_xyy,  \
                             tr_y_xxxy_xyz,  \
                             tr_y_xxxy_yyy,  \
                             tr_y_xxxyz_xxx, \
                             tr_y_xxxyz_xxy, \
                             tr_y_xxxyz_xxz, \
                             tr_y_xxxyz_xyy, \
                             tr_y_xxxyz_xyz, \
                             tr_y_xxxyz_xzz, \
                             tr_y_xxxyz_yyy, \
                             tr_y_xxxyz_yyz, \
                             tr_y_xxxyz_yzz, \
                             tr_y_xxxyz_zzz, \
                             tr_y_xxxz_xxz,  \
                             tr_y_xxxz_xzz,  \
                             tr_y_xxyz_yyz,  \
                             tr_y_xxyz_yzz,  \
                             tr_y_xxyz_zzz,  \
                             tr_y_xyz_yyz,   \
                             tr_y_xyz_yzz,   \
                             tr_y_xyz_zzz,   \
                             ts_xxxz_xxz,    \
                             ts_xxxz_xzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxyz_xxx[i] = tr_y_xxxy_xxx[i] * pa_z[i];

        tr_y_xxxyz_xxy[i] = tr_y_xxxy_xxy[i] * pa_z[i];

        tr_y_xxxyz_xxz[i] = ts_xxxz_xxz[i] * fe_0 + tr_y_xxxz_xxz[i] * pa_y[i];

        tr_y_xxxyz_xyy[i] = tr_y_xxxy_xyy[i] * pa_z[i];

        tr_y_xxxyz_xyz[i] = tr_y_xxxy_xy[i] * fe_0 + tr_y_xxxy_xyz[i] * pa_z[i];

        tr_y_xxxyz_xzz[i] = ts_xxxz_xzz[i] * fe_0 + tr_y_xxxz_xzz[i] * pa_y[i];

        tr_y_xxxyz_yyy[i] = tr_y_xxxy_yyy[i] * pa_z[i];

        tr_y_xxxyz_yyz[i] = 2.0 * tr_y_xyz_yyz[i] * fe_0 + tr_y_xxyz_yyz[i] * pa_x[i];

        tr_y_xxxyz_yzz[i] = 2.0 * tr_y_xyz_yzz[i] * fe_0 + tr_y_xxyz_yzz[i] * pa_x[i];

        tr_y_xxxyz_zzz[i] = 2.0 * tr_y_xyz_zzz[i] * fe_0 + tr_y_xxyz_zzz[i] * pa_x[i];
    }

    // Set up 260-270 components of targeted buffer : HF

    auto tr_y_xxxzz_xxx = pbuffer.data(idx_dip_hf + 260);

    auto tr_y_xxxzz_xxy = pbuffer.data(idx_dip_hf + 261);

    auto tr_y_xxxzz_xxz = pbuffer.data(idx_dip_hf + 262);

    auto tr_y_xxxzz_xyy = pbuffer.data(idx_dip_hf + 263);

    auto tr_y_xxxzz_xyz = pbuffer.data(idx_dip_hf + 264);

    auto tr_y_xxxzz_xzz = pbuffer.data(idx_dip_hf + 265);

    auto tr_y_xxxzz_yyy = pbuffer.data(idx_dip_hf + 266);

    auto tr_y_xxxzz_yyz = pbuffer.data(idx_dip_hf + 267);

    auto tr_y_xxxzz_yzz = pbuffer.data(idx_dip_hf + 268);

    auto tr_y_xxxzz_zzz = pbuffer.data(idx_dip_hf + 269);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xxx_xxx,   \
                             tr_y_xxx_xxy,   \
                             tr_y_xxx_xyy,   \
                             tr_y_xxxz_xxx,  \
                             tr_y_xxxz_xxy,  \
                             tr_y_xxxz_xyy,  \
                             tr_y_xxxzz_xxx, \
                             tr_y_xxxzz_xxy, \
                             tr_y_xxxzz_xxz, \
                             tr_y_xxxzz_xyy, \
                             tr_y_xxxzz_xyz, \
                             tr_y_xxxzz_xzz, \
                             tr_y_xxxzz_yyy, \
                             tr_y_xxxzz_yyz, \
                             tr_y_xxxzz_yzz, \
                             tr_y_xxxzz_zzz, \
                             tr_y_xxzz_xxz,  \
                             tr_y_xxzz_xyz,  \
                             tr_y_xxzz_xz,   \
                             tr_y_xxzz_xzz,  \
                             tr_y_xxzz_yyy,  \
                             tr_y_xxzz_yyz,  \
                             tr_y_xxzz_yz,   \
                             tr_y_xxzz_yzz,  \
                             tr_y_xxzz_zz,   \
                             tr_y_xxzz_zzz,  \
                             tr_y_xzz_xxz,   \
                             tr_y_xzz_xyz,   \
                             tr_y_xzz_xzz,   \
                             tr_y_xzz_yyy,   \
                             tr_y_xzz_yyz,   \
                             tr_y_xzz_yzz,   \
                             tr_y_xzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxxzz_xxx[i] = tr_y_xxx_xxx[i] * fe_0 + tr_y_xxxz_xxx[i] * pa_z[i];

        tr_y_xxxzz_xxy[i] = tr_y_xxx_xxy[i] * fe_0 + tr_y_xxxz_xxy[i] * pa_z[i];

        tr_y_xxxzz_xxz[i] = 2.0 * tr_y_xzz_xxz[i] * fe_0 + 2.0 * tr_y_xxzz_xz[i] * fe_0 + tr_y_xxzz_xxz[i] * pa_x[i];

        tr_y_xxxzz_xyy[i] = tr_y_xxx_xyy[i] * fe_0 + tr_y_xxxz_xyy[i] * pa_z[i];

        tr_y_xxxzz_xyz[i] = 2.0 * tr_y_xzz_xyz[i] * fe_0 + tr_y_xxzz_yz[i] * fe_0 + tr_y_xxzz_xyz[i] * pa_x[i];

        tr_y_xxxzz_xzz[i] = 2.0 * tr_y_xzz_xzz[i] * fe_0 + tr_y_xxzz_zz[i] * fe_0 + tr_y_xxzz_xzz[i] * pa_x[i];

        tr_y_xxxzz_yyy[i] = 2.0 * tr_y_xzz_yyy[i] * fe_0 + tr_y_xxzz_yyy[i] * pa_x[i];

        tr_y_xxxzz_yyz[i] = 2.0 * tr_y_xzz_yyz[i] * fe_0 + tr_y_xxzz_yyz[i] * pa_x[i];

        tr_y_xxxzz_yzz[i] = 2.0 * tr_y_xzz_yzz[i] * fe_0 + tr_y_xxzz_yzz[i] * pa_x[i];

        tr_y_xxxzz_zzz[i] = 2.0 * tr_y_xzz_zzz[i] * fe_0 + tr_y_xxzz_zzz[i] * pa_x[i];
    }

    // Set up 270-280 components of targeted buffer : HF

    auto tr_y_xxyyy_xxx = pbuffer.data(idx_dip_hf + 270);

    auto tr_y_xxyyy_xxy = pbuffer.data(idx_dip_hf + 271);

    auto tr_y_xxyyy_xxz = pbuffer.data(idx_dip_hf + 272);

    auto tr_y_xxyyy_xyy = pbuffer.data(idx_dip_hf + 273);

    auto tr_y_xxyyy_xyz = pbuffer.data(idx_dip_hf + 274);

    auto tr_y_xxyyy_xzz = pbuffer.data(idx_dip_hf + 275);

    auto tr_y_xxyyy_yyy = pbuffer.data(idx_dip_hf + 276);

    auto tr_y_xxyyy_yyz = pbuffer.data(idx_dip_hf + 277);

    auto tr_y_xxyyy_yzz = pbuffer.data(idx_dip_hf + 278);

    auto tr_y_xxyyy_zzz = pbuffer.data(idx_dip_hf + 279);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xxyyy_xxx, \
                             tr_y_xxyyy_xxy, \
                             tr_y_xxyyy_xxz, \
                             tr_y_xxyyy_xyy, \
                             tr_y_xxyyy_xyz, \
                             tr_y_xxyyy_xzz, \
                             tr_y_xxyyy_yyy, \
                             tr_y_xxyyy_yyz, \
                             tr_y_xxyyy_yzz, \
                             tr_y_xxyyy_zzz, \
                             tr_y_xyyy_xx,   \
                             tr_y_xyyy_xxx,  \
                             tr_y_xyyy_xxy,  \
                             tr_y_xyyy_xxz,  \
                             tr_y_xyyy_xy,   \
                             tr_y_xyyy_xyy,  \
                             tr_y_xyyy_xyz,  \
                             tr_y_xyyy_xz,   \
                             tr_y_xyyy_xzz,  \
                             tr_y_xyyy_yy,   \
                             tr_y_xyyy_yyy,  \
                             tr_y_xyyy_yyz,  \
                             tr_y_xyyy_yz,   \
                             tr_y_xyyy_yzz,  \
                             tr_y_xyyy_zz,   \
                             tr_y_xyyy_zzz,  \
                             tr_y_yyy_xxx,   \
                             tr_y_yyy_xxy,   \
                             tr_y_yyy_xxz,   \
                             tr_y_yyy_xyy,   \
                             tr_y_yyy_xyz,   \
                             tr_y_yyy_xzz,   \
                             tr_y_yyy_yyy,   \
                             tr_y_yyy_yyz,   \
                             tr_y_yyy_yzz,   \
                             tr_y_yyy_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyy_xxx[i] = tr_y_yyy_xxx[i] * fe_0 + 3.0 * tr_y_xyyy_xx[i] * fe_0 + tr_y_xyyy_xxx[i] * pa_x[i];

        tr_y_xxyyy_xxy[i] = tr_y_yyy_xxy[i] * fe_0 + 2.0 * tr_y_xyyy_xy[i] * fe_0 + tr_y_xyyy_xxy[i] * pa_x[i];

        tr_y_xxyyy_xxz[i] = tr_y_yyy_xxz[i] * fe_0 + 2.0 * tr_y_xyyy_xz[i] * fe_0 + tr_y_xyyy_xxz[i] * pa_x[i];

        tr_y_xxyyy_xyy[i] = tr_y_yyy_xyy[i] * fe_0 + tr_y_xyyy_yy[i] * fe_0 + tr_y_xyyy_xyy[i] * pa_x[i];

        tr_y_xxyyy_xyz[i] = tr_y_yyy_xyz[i] * fe_0 + tr_y_xyyy_yz[i] * fe_0 + tr_y_xyyy_xyz[i] * pa_x[i];

        tr_y_xxyyy_xzz[i] = tr_y_yyy_xzz[i] * fe_0 + tr_y_xyyy_zz[i] * fe_0 + tr_y_xyyy_xzz[i] * pa_x[i];

        tr_y_xxyyy_yyy[i] = tr_y_yyy_yyy[i] * fe_0 + tr_y_xyyy_yyy[i] * pa_x[i];

        tr_y_xxyyy_yyz[i] = tr_y_yyy_yyz[i] * fe_0 + tr_y_xyyy_yyz[i] * pa_x[i];

        tr_y_xxyyy_yzz[i] = tr_y_yyy_yzz[i] * fe_0 + tr_y_xyyy_yzz[i] * pa_x[i];

        tr_y_xxyyy_zzz[i] = tr_y_yyy_zzz[i] * fe_0 + tr_y_xyyy_zzz[i] * pa_x[i];
    }

    // Set up 280-290 components of targeted buffer : HF

    auto tr_y_xxyyz_xxx = pbuffer.data(idx_dip_hf + 280);

    auto tr_y_xxyyz_xxy = pbuffer.data(idx_dip_hf + 281);

    auto tr_y_xxyyz_xxz = pbuffer.data(idx_dip_hf + 282);

    auto tr_y_xxyyz_xyy = pbuffer.data(idx_dip_hf + 283);

    auto tr_y_xxyyz_xyz = pbuffer.data(idx_dip_hf + 284);

    auto tr_y_xxyyz_xzz = pbuffer.data(idx_dip_hf + 285);

    auto tr_y_xxyyz_yyy = pbuffer.data(idx_dip_hf + 286);

    auto tr_y_xxyyz_yyz = pbuffer.data(idx_dip_hf + 287);

    auto tr_y_xxyyz_yzz = pbuffer.data(idx_dip_hf + 288);

    auto tr_y_xxyyz_zzz = pbuffer.data(idx_dip_hf + 289);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xxyy_xx,   \
                             tr_y_xxyy_xxx,  \
                             tr_y_xxyy_xxy,  \
                             tr_y_xxyy_xxz,  \
                             tr_y_xxyy_xy,   \
                             tr_y_xxyy_xyy,  \
                             tr_y_xxyy_xyz,  \
                             tr_y_xxyy_xz,   \
                             tr_y_xxyy_xzz,  \
                             tr_y_xxyy_yyy,  \
                             tr_y_xxyyz_xxx, \
                             tr_y_xxyyz_xxy, \
                             tr_y_xxyyz_xxz, \
                             tr_y_xxyyz_xyy, \
                             tr_y_xxyyz_xyz, \
                             tr_y_xxyyz_xzz, \
                             tr_y_xxyyz_yyy, \
                             tr_y_xxyyz_yyz, \
                             tr_y_xxyyz_yzz, \
                             tr_y_xxyyz_zzz, \
                             tr_y_xyyz_yyz,  \
                             tr_y_xyyz_yzz,  \
                             tr_y_xyyz_zzz,  \
                             tr_y_yyz_yyz,   \
                             tr_y_yyz_yzz,   \
                             tr_y_yyz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyyz_xxx[i] = tr_y_xxyy_xxx[i] * pa_z[i];

        tr_y_xxyyz_xxy[i] = tr_y_xxyy_xxy[i] * pa_z[i];

        tr_y_xxyyz_xxz[i] = tr_y_xxyy_xx[i] * fe_0 + tr_y_xxyy_xxz[i] * pa_z[i];

        tr_y_xxyyz_xyy[i] = tr_y_xxyy_xyy[i] * pa_z[i];

        tr_y_xxyyz_xyz[i] = tr_y_xxyy_xy[i] * fe_0 + tr_y_xxyy_xyz[i] * pa_z[i];

        tr_y_xxyyz_xzz[i] = 2.0 * tr_y_xxyy_xz[i] * fe_0 + tr_y_xxyy_xzz[i] * pa_z[i];

        tr_y_xxyyz_yyy[i] = tr_y_xxyy_yyy[i] * pa_z[i];

        tr_y_xxyyz_yyz[i] = tr_y_yyz_yyz[i] * fe_0 + tr_y_xyyz_yyz[i] * pa_x[i];

        tr_y_xxyyz_yzz[i] = tr_y_yyz_yzz[i] * fe_0 + tr_y_xyyz_yzz[i] * pa_x[i];

        tr_y_xxyyz_zzz[i] = tr_y_yyz_zzz[i] * fe_0 + tr_y_xyyz_zzz[i] * pa_x[i];
    }

    // Set up 290-300 components of targeted buffer : HF

    auto tr_y_xxyzz_xxx = pbuffer.data(idx_dip_hf + 290);

    auto tr_y_xxyzz_xxy = pbuffer.data(idx_dip_hf + 291);

    auto tr_y_xxyzz_xxz = pbuffer.data(idx_dip_hf + 292);

    auto tr_y_xxyzz_xyy = pbuffer.data(idx_dip_hf + 293);

    auto tr_y_xxyzz_xyz = pbuffer.data(idx_dip_hf + 294);

    auto tr_y_xxyzz_xzz = pbuffer.data(idx_dip_hf + 295);

    auto tr_y_xxyzz_yyy = pbuffer.data(idx_dip_hf + 296);

    auto tr_y_xxyzz_yyz = pbuffer.data(idx_dip_hf + 297);

    auto tr_y_xxyzz_yzz = pbuffer.data(idx_dip_hf + 298);

    auto tr_y_xxyzz_zzz = pbuffer.data(idx_dip_hf + 299);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_y_xxy_xxy,   \
                             tr_y_xxy_xyy,   \
                             tr_y_xxyz_xxy,  \
                             tr_y_xxyz_xyy,  \
                             tr_y_xxyzz_xxx, \
                             tr_y_xxyzz_xxy, \
                             tr_y_xxyzz_xxz, \
                             tr_y_xxyzz_xyy, \
                             tr_y_xxyzz_xyz, \
                             tr_y_xxyzz_xzz, \
                             tr_y_xxyzz_yyy, \
                             tr_y_xxyzz_yyz, \
                             tr_y_xxyzz_yzz, \
                             tr_y_xxyzz_zzz, \
                             tr_y_xxzz_xxx,  \
                             tr_y_xxzz_xxz,  \
                             tr_y_xxzz_xzz,  \
                             tr_y_xyzz_xyz,  \
                             tr_y_xyzz_yyy,  \
                             tr_y_xyzz_yyz,  \
                             tr_y_xyzz_yz,   \
                             tr_y_xyzz_yzz,  \
                             tr_y_xyzz_zzz,  \
                             tr_y_yzz_xyz,   \
                             tr_y_yzz_yyy,   \
                             tr_y_yzz_yyz,   \
                             tr_y_yzz_yzz,   \
                             tr_y_yzz_zzz,   \
                             ts_xxzz_xxx,    \
                             ts_xxzz_xxz,    \
                             ts_xxzz_xzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxyzz_xxx[i] = ts_xxzz_xxx[i] * fe_0 + tr_y_xxzz_xxx[i] * pa_y[i];

        tr_y_xxyzz_xxy[i] = tr_y_xxy_xxy[i] * fe_0 + tr_y_xxyz_xxy[i] * pa_z[i];

        tr_y_xxyzz_xxz[i] = ts_xxzz_xxz[i] * fe_0 + tr_y_xxzz_xxz[i] * pa_y[i];

        tr_y_xxyzz_xyy[i] = tr_y_xxy_xyy[i] * fe_0 + tr_y_xxyz_xyy[i] * pa_z[i];

        tr_y_xxyzz_xyz[i] = tr_y_yzz_xyz[i] * fe_0 + tr_y_xyzz_yz[i] * fe_0 + tr_y_xyzz_xyz[i] * pa_x[i];

        tr_y_xxyzz_xzz[i] = ts_xxzz_xzz[i] * fe_0 + tr_y_xxzz_xzz[i] * pa_y[i];

        tr_y_xxyzz_yyy[i] = tr_y_yzz_yyy[i] * fe_0 + tr_y_xyzz_yyy[i] * pa_x[i];

        tr_y_xxyzz_yyz[i] = tr_y_yzz_yyz[i] * fe_0 + tr_y_xyzz_yyz[i] * pa_x[i];

        tr_y_xxyzz_yzz[i] = tr_y_yzz_yzz[i] * fe_0 + tr_y_xyzz_yzz[i] * pa_x[i];

        tr_y_xxyzz_zzz[i] = tr_y_yzz_zzz[i] * fe_0 + tr_y_xyzz_zzz[i] * pa_x[i];
    }

    // Set up 300-310 components of targeted buffer : HF

    auto tr_y_xxzzz_xxx = pbuffer.data(idx_dip_hf + 300);

    auto tr_y_xxzzz_xxy = pbuffer.data(idx_dip_hf + 301);

    auto tr_y_xxzzz_xxz = pbuffer.data(idx_dip_hf + 302);

    auto tr_y_xxzzz_xyy = pbuffer.data(idx_dip_hf + 303);

    auto tr_y_xxzzz_xyz = pbuffer.data(idx_dip_hf + 304);

    auto tr_y_xxzzz_xzz = pbuffer.data(idx_dip_hf + 305);

    auto tr_y_xxzzz_yyy = pbuffer.data(idx_dip_hf + 306);

    auto tr_y_xxzzz_yyz = pbuffer.data(idx_dip_hf + 307);

    auto tr_y_xxzzz_yzz = pbuffer.data(idx_dip_hf + 308);

    auto tr_y_xxzzz_zzz = pbuffer.data(idx_dip_hf + 309);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xxz_xxx,   \
                             tr_y_xxz_xxy,   \
                             tr_y_xxz_xyy,   \
                             tr_y_xxzz_xxx,  \
                             tr_y_xxzz_xxy,  \
                             tr_y_xxzz_xyy,  \
                             tr_y_xxzzz_xxx, \
                             tr_y_xxzzz_xxy, \
                             tr_y_xxzzz_xxz, \
                             tr_y_xxzzz_xyy, \
                             tr_y_xxzzz_xyz, \
                             tr_y_xxzzz_xzz, \
                             tr_y_xxzzz_yyy, \
                             tr_y_xxzzz_yyz, \
                             tr_y_xxzzz_yzz, \
                             tr_y_xxzzz_zzz, \
                             tr_y_xzzz_xxz,  \
                             tr_y_xzzz_xyz,  \
                             tr_y_xzzz_xz,   \
                             tr_y_xzzz_xzz,  \
                             tr_y_xzzz_yyy,  \
                             tr_y_xzzz_yyz,  \
                             tr_y_xzzz_yz,   \
                             tr_y_xzzz_yzz,  \
                             tr_y_xzzz_zz,   \
                             tr_y_xzzz_zzz,  \
                             tr_y_zzz_xxz,   \
                             tr_y_zzz_xyz,   \
                             tr_y_zzz_xzz,   \
                             tr_y_zzz_yyy,   \
                             tr_y_zzz_yyz,   \
                             tr_y_zzz_yzz,   \
                             tr_y_zzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxzzz_xxx[i] = 2.0 * tr_y_xxz_xxx[i] * fe_0 + tr_y_xxzz_xxx[i] * pa_z[i];

        tr_y_xxzzz_xxy[i] = 2.0 * tr_y_xxz_xxy[i] * fe_0 + tr_y_xxzz_xxy[i] * pa_z[i];

        tr_y_xxzzz_xxz[i] = tr_y_zzz_xxz[i] * fe_0 + 2.0 * tr_y_xzzz_xz[i] * fe_0 + tr_y_xzzz_xxz[i] * pa_x[i];

        tr_y_xxzzz_xyy[i] = 2.0 * tr_y_xxz_xyy[i] * fe_0 + tr_y_xxzz_xyy[i] * pa_z[i];

        tr_y_xxzzz_xyz[i] = tr_y_zzz_xyz[i] * fe_0 + tr_y_xzzz_yz[i] * fe_0 + tr_y_xzzz_xyz[i] * pa_x[i];

        tr_y_xxzzz_xzz[i] = tr_y_zzz_xzz[i] * fe_0 + tr_y_xzzz_zz[i] * fe_0 + tr_y_xzzz_xzz[i] * pa_x[i];

        tr_y_xxzzz_yyy[i] = tr_y_zzz_yyy[i] * fe_0 + tr_y_xzzz_yyy[i] * pa_x[i];

        tr_y_xxzzz_yyz[i] = tr_y_zzz_yyz[i] * fe_0 + tr_y_xzzz_yyz[i] * pa_x[i];

        tr_y_xxzzz_yzz[i] = tr_y_zzz_yzz[i] * fe_0 + tr_y_xzzz_yzz[i] * pa_x[i];

        tr_y_xxzzz_zzz[i] = tr_y_zzz_zzz[i] * fe_0 + tr_y_xzzz_zzz[i] * pa_x[i];
    }

    // Set up 310-320 components of targeted buffer : HF

    auto tr_y_xyyyy_xxx = pbuffer.data(idx_dip_hf + 310);

    auto tr_y_xyyyy_xxy = pbuffer.data(idx_dip_hf + 311);

    auto tr_y_xyyyy_xxz = pbuffer.data(idx_dip_hf + 312);

    auto tr_y_xyyyy_xyy = pbuffer.data(idx_dip_hf + 313);

    auto tr_y_xyyyy_xyz = pbuffer.data(idx_dip_hf + 314);

    auto tr_y_xyyyy_xzz = pbuffer.data(idx_dip_hf + 315);

    auto tr_y_xyyyy_yyy = pbuffer.data(idx_dip_hf + 316);

    auto tr_y_xyyyy_yyz = pbuffer.data(idx_dip_hf + 317);

    auto tr_y_xyyyy_yzz = pbuffer.data(idx_dip_hf + 318);

    auto tr_y_xyyyy_zzz = pbuffer.data(idx_dip_hf + 319);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xyyyy_xxx, \
                             tr_y_xyyyy_xxy, \
                             tr_y_xyyyy_xxz, \
                             tr_y_xyyyy_xyy, \
                             tr_y_xyyyy_xyz, \
                             tr_y_xyyyy_xzz, \
                             tr_y_xyyyy_yyy, \
                             tr_y_xyyyy_yyz, \
                             tr_y_xyyyy_yzz, \
                             tr_y_xyyyy_zzz, \
                             tr_y_yyyy_xx,   \
                             tr_y_yyyy_xxx,  \
                             tr_y_yyyy_xxy,  \
                             tr_y_yyyy_xxz,  \
                             tr_y_yyyy_xy,   \
                             tr_y_yyyy_xyy,  \
                             tr_y_yyyy_xyz,  \
                             tr_y_yyyy_xz,   \
                             tr_y_yyyy_xzz,  \
                             tr_y_yyyy_yy,   \
                             tr_y_yyyy_yyy,  \
                             tr_y_yyyy_yyz,  \
                             tr_y_yyyy_yz,   \
                             tr_y_yyyy_yzz,  \
                             tr_y_yyyy_zz,   \
                             tr_y_yyyy_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyy_xxx[i] = 3.0 * tr_y_yyyy_xx[i] * fe_0 + tr_y_yyyy_xxx[i] * pa_x[i];

        tr_y_xyyyy_xxy[i] = 2.0 * tr_y_yyyy_xy[i] * fe_0 + tr_y_yyyy_xxy[i] * pa_x[i];

        tr_y_xyyyy_xxz[i] = 2.0 * tr_y_yyyy_xz[i] * fe_0 + tr_y_yyyy_xxz[i] * pa_x[i];

        tr_y_xyyyy_xyy[i] = tr_y_yyyy_yy[i] * fe_0 + tr_y_yyyy_xyy[i] * pa_x[i];

        tr_y_xyyyy_xyz[i] = tr_y_yyyy_yz[i] * fe_0 + tr_y_yyyy_xyz[i] * pa_x[i];

        tr_y_xyyyy_xzz[i] = tr_y_yyyy_zz[i] * fe_0 + tr_y_yyyy_xzz[i] * pa_x[i];

        tr_y_xyyyy_yyy[i] = tr_y_yyyy_yyy[i] * pa_x[i];

        tr_y_xyyyy_yyz[i] = tr_y_yyyy_yyz[i] * pa_x[i];

        tr_y_xyyyy_yzz[i] = tr_y_yyyy_yzz[i] * pa_x[i];

        tr_y_xyyyy_zzz[i] = tr_y_yyyy_zzz[i] * pa_x[i];
    }

    // Set up 320-330 components of targeted buffer : HF

    auto tr_y_xyyyz_xxx = pbuffer.data(idx_dip_hf + 320);

    auto tr_y_xyyyz_xxy = pbuffer.data(idx_dip_hf + 321);

    auto tr_y_xyyyz_xxz = pbuffer.data(idx_dip_hf + 322);

    auto tr_y_xyyyz_xyy = pbuffer.data(idx_dip_hf + 323);

    auto tr_y_xyyyz_xyz = pbuffer.data(idx_dip_hf + 324);

    auto tr_y_xyyyz_xzz = pbuffer.data(idx_dip_hf + 325);

    auto tr_y_xyyyz_yyy = pbuffer.data(idx_dip_hf + 326);

    auto tr_y_xyyyz_yyz = pbuffer.data(idx_dip_hf + 327);

    auto tr_y_xyyyz_yzz = pbuffer.data(idx_dip_hf + 328);

    auto tr_y_xyyyz_zzz = pbuffer.data(idx_dip_hf + 329);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xyyy_xxx,  \
                             tr_y_xyyy_xxy,  \
                             tr_y_xyyy_xyy,  \
                             tr_y_xyyyz_xxx, \
                             tr_y_xyyyz_xxy, \
                             tr_y_xyyyz_xxz, \
                             tr_y_xyyyz_xyy, \
                             tr_y_xyyyz_xyz, \
                             tr_y_xyyyz_xzz, \
                             tr_y_xyyyz_yyy, \
                             tr_y_xyyyz_yyz, \
                             tr_y_xyyyz_yzz, \
                             tr_y_xyyyz_zzz, \
                             tr_y_yyyz_xxz,  \
                             tr_y_yyyz_xyz,  \
                             tr_y_yyyz_xz,   \
                             tr_y_yyyz_xzz,  \
                             tr_y_yyyz_yyy,  \
                             tr_y_yyyz_yyz,  \
                             tr_y_yyyz_yz,   \
                             tr_y_yyyz_yzz,  \
                             tr_y_yyyz_zz,   \
                             tr_y_yyyz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyyz_xxx[i] = tr_y_xyyy_xxx[i] * pa_z[i];

        tr_y_xyyyz_xxy[i] = tr_y_xyyy_xxy[i] * pa_z[i];

        tr_y_xyyyz_xxz[i] = 2.0 * tr_y_yyyz_xz[i] * fe_0 + tr_y_yyyz_xxz[i] * pa_x[i];

        tr_y_xyyyz_xyy[i] = tr_y_xyyy_xyy[i] * pa_z[i];

        tr_y_xyyyz_xyz[i] = tr_y_yyyz_yz[i] * fe_0 + tr_y_yyyz_xyz[i] * pa_x[i];

        tr_y_xyyyz_xzz[i] = tr_y_yyyz_zz[i] * fe_0 + tr_y_yyyz_xzz[i] * pa_x[i];

        tr_y_xyyyz_yyy[i] = tr_y_yyyz_yyy[i] * pa_x[i];

        tr_y_xyyyz_yyz[i] = tr_y_yyyz_yyz[i] * pa_x[i];

        tr_y_xyyyz_yzz[i] = tr_y_yyyz_yzz[i] * pa_x[i];

        tr_y_xyyyz_zzz[i] = tr_y_yyyz_zzz[i] * pa_x[i];
    }

    // Set up 330-340 components of targeted buffer : HF

    auto tr_y_xyyzz_xxx = pbuffer.data(idx_dip_hf + 330);

    auto tr_y_xyyzz_xxy = pbuffer.data(idx_dip_hf + 331);

    auto tr_y_xyyzz_xxz = pbuffer.data(idx_dip_hf + 332);

    auto tr_y_xyyzz_xyy = pbuffer.data(idx_dip_hf + 333);

    auto tr_y_xyyzz_xyz = pbuffer.data(idx_dip_hf + 334);

    auto tr_y_xyyzz_xzz = pbuffer.data(idx_dip_hf + 335);

    auto tr_y_xyyzz_yyy = pbuffer.data(idx_dip_hf + 336);

    auto tr_y_xyyzz_yyz = pbuffer.data(idx_dip_hf + 337);

    auto tr_y_xyyzz_yzz = pbuffer.data(idx_dip_hf + 338);

    auto tr_y_xyyzz_zzz = pbuffer.data(idx_dip_hf + 339);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xyyzz_xxx, \
                             tr_y_xyyzz_xxy, \
                             tr_y_xyyzz_xxz, \
                             tr_y_xyyzz_xyy, \
                             tr_y_xyyzz_xyz, \
                             tr_y_xyyzz_xzz, \
                             tr_y_xyyzz_yyy, \
                             tr_y_xyyzz_yyz, \
                             tr_y_xyyzz_yzz, \
                             tr_y_xyyzz_zzz, \
                             tr_y_yyzz_xx,   \
                             tr_y_yyzz_xxx,  \
                             tr_y_yyzz_xxy,  \
                             tr_y_yyzz_xxz,  \
                             tr_y_yyzz_xy,   \
                             tr_y_yyzz_xyy,  \
                             tr_y_yyzz_xyz,  \
                             tr_y_yyzz_xz,   \
                             tr_y_yyzz_xzz,  \
                             tr_y_yyzz_yy,   \
                             tr_y_yyzz_yyy,  \
                             tr_y_yyzz_yyz,  \
                             tr_y_yyzz_yz,   \
                             tr_y_yyzz_yzz,  \
                             tr_y_yyzz_zz,   \
                             tr_y_yyzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyyzz_xxx[i] = 3.0 * tr_y_yyzz_xx[i] * fe_0 + tr_y_yyzz_xxx[i] * pa_x[i];

        tr_y_xyyzz_xxy[i] = 2.0 * tr_y_yyzz_xy[i] * fe_0 + tr_y_yyzz_xxy[i] * pa_x[i];

        tr_y_xyyzz_xxz[i] = 2.0 * tr_y_yyzz_xz[i] * fe_0 + tr_y_yyzz_xxz[i] * pa_x[i];

        tr_y_xyyzz_xyy[i] = tr_y_yyzz_yy[i] * fe_0 + tr_y_yyzz_xyy[i] * pa_x[i];

        tr_y_xyyzz_xyz[i] = tr_y_yyzz_yz[i] * fe_0 + tr_y_yyzz_xyz[i] * pa_x[i];

        tr_y_xyyzz_xzz[i] = tr_y_yyzz_zz[i] * fe_0 + tr_y_yyzz_xzz[i] * pa_x[i];

        tr_y_xyyzz_yyy[i] = tr_y_yyzz_yyy[i] * pa_x[i];

        tr_y_xyyzz_yyz[i] = tr_y_yyzz_yyz[i] * pa_x[i];

        tr_y_xyyzz_yzz[i] = tr_y_yyzz_yzz[i] * pa_x[i];

        tr_y_xyyzz_zzz[i] = tr_y_yyzz_zzz[i] * pa_x[i];
    }

    // Set up 340-350 components of targeted buffer : HF

    auto tr_y_xyzzz_xxx = pbuffer.data(idx_dip_hf + 340);

    auto tr_y_xyzzz_xxy = pbuffer.data(idx_dip_hf + 341);

    auto tr_y_xyzzz_xxz = pbuffer.data(idx_dip_hf + 342);

    auto tr_y_xyzzz_xyy = pbuffer.data(idx_dip_hf + 343);

    auto tr_y_xyzzz_xyz = pbuffer.data(idx_dip_hf + 344);

    auto tr_y_xyzzz_xzz = pbuffer.data(idx_dip_hf + 345);

    auto tr_y_xyzzz_yyy = pbuffer.data(idx_dip_hf + 346);

    auto tr_y_xyzzz_yyz = pbuffer.data(idx_dip_hf + 347);

    auto tr_y_xyzzz_yzz = pbuffer.data(idx_dip_hf + 348);

    auto tr_y_xyzzz_zzz = pbuffer.data(idx_dip_hf + 349);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xyzzz_xxx, \
                             tr_y_xyzzz_xxy, \
                             tr_y_xyzzz_xxz, \
                             tr_y_xyzzz_xyy, \
                             tr_y_xyzzz_xyz, \
                             tr_y_xyzzz_xzz, \
                             tr_y_xyzzz_yyy, \
                             tr_y_xyzzz_yyz, \
                             tr_y_xyzzz_yzz, \
                             tr_y_xyzzz_zzz, \
                             tr_y_yzzz_xx,   \
                             tr_y_yzzz_xxx,  \
                             tr_y_yzzz_xxy,  \
                             tr_y_yzzz_xxz,  \
                             tr_y_yzzz_xy,   \
                             tr_y_yzzz_xyy,  \
                             tr_y_yzzz_xyz,  \
                             tr_y_yzzz_xz,   \
                             tr_y_yzzz_xzz,  \
                             tr_y_yzzz_yy,   \
                             tr_y_yzzz_yyy,  \
                             tr_y_yzzz_yyz,  \
                             tr_y_yzzz_yz,   \
                             tr_y_yzzz_yzz,  \
                             tr_y_yzzz_zz,   \
                             tr_y_yzzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyzzz_xxx[i] = 3.0 * tr_y_yzzz_xx[i] * fe_0 + tr_y_yzzz_xxx[i] * pa_x[i];

        tr_y_xyzzz_xxy[i] = 2.0 * tr_y_yzzz_xy[i] * fe_0 + tr_y_yzzz_xxy[i] * pa_x[i];

        tr_y_xyzzz_xxz[i] = 2.0 * tr_y_yzzz_xz[i] * fe_0 + tr_y_yzzz_xxz[i] * pa_x[i];

        tr_y_xyzzz_xyy[i] = tr_y_yzzz_yy[i] * fe_0 + tr_y_yzzz_xyy[i] * pa_x[i];

        tr_y_xyzzz_xyz[i] = tr_y_yzzz_yz[i] * fe_0 + tr_y_yzzz_xyz[i] * pa_x[i];

        tr_y_xyzzz_xzz[i] = tr_y_yzzz_zz[i] * fe_0 + tr_y_yzzz_xzz[i] * pa_x[i];

        tr_y_xyzzz_yyy[i] = tr_y_yzzz_yyy[i] * pa_x[i];

        tr_y_xyzzz_yyz[i] = tr_y_yzzz_yyz[i] * pa_x[i];

        tr_y_xyzzz_yzz[i] = tr_y_yzzz_yzz[i] * pa_x[i];

        tr_y_xyzzz_zzz[i] = tr_y_yzzz_zzz[i] * pa_x[i];
    }

    // Set up 350-360 components of targeted buffer : HF

    auto tr_y_xzzzz_xxx = pbuffer.data(idx_dip_hf + 350);

    auto tr_y_xzzzz_xxy = pbuffer.data(idx_dip_hf + 351);

    auto tr_y_xzzzz_xxz = pbuffer.data(idx_dip_hf + 352);

    auto tr_y_xzzzz_xyy = pbuffer.data(idx_dip_hf + 353);

    auto tr_y_xzzzz_xyz = pbuffer.data(idx_dip_hf + 354);

    auto tr_y_xzzzz_xzz = pbuffer.data(idx_dip_hf + 355);

    auto tr_y_xzzzz_yyy = pbuffer.data(idx_dip_hf + 356);

    auto tr_y_xzzzz_yyz = pbuffer.data(idx_dip_hf + 357);

    auto tr_y_xzzzz_yzz = pbuffer.data(idx_dip_hf + 358);

    auto tr_y_xzzzz_zzz = pbuffer.data(idx_dip_hf + 359);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xzzzz_xxx, \
                             tr_y_xzzzz_xxy, \
                             tr_y_xzzzz_xxz, \
                             tr_y_xzzzz_xyy, \
                             tr_y_xzzzz_xyz, \
                             tr_y_xzzzz_xzz, \
                             tr_y_xzzzz_yyy, \
                             tr_y_xzzzz_yyz, \
                             tr_y_xzzzz_yzz, \
                             tr_y_xzzzz_zzz, \
                             tr_y_zzzz_xx,   \
                             tr_y_zzzz_xxx,  \
                             tr_y_zzzz_xxy,  \
                             tr_y_zzzz_xxz,  \
                             tr_y_zzzz_xy,   \
                             tr_y_zzzz_xyy,  \
                             tr_y_zzzz_xyz,  \
                             tr_y_zzzz_xz,   \
                             tr_y_zzzz_xzz,  \
                             tr_y_zzzz_yy,   \
                             tr_y_zzzz_yyy,  \
                             tr_y_zzzz_yyz,  \
                             tr_y_zzzz_yz,   \
                             tr_y_zzzz_yzz,  \
                             tr_y_zzzz_zz,   \
                             tr_y_zzzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzzzz_xxx[i] = 3.0 * tr_y_zzzz_xx[i] * fe_0 + tr_y_zzzz_xxx[i] * pa_x[i];

        tr_y_xzzzz_xxy[i] = 2.0 * tr_y_zzzz_xy[i] * fe_0 + tr_y_zzzz_xxy[i] * pa_x[i];

        tr_y_xzzzz_xxz[i] = 2.0 * tr_y_zzzz_xz[i] * fe_0 + tr_y_zzzz_xxz[i] * pa_x[i];

        tr_y_xzzzz_xyy[i] = tr_y_zzzz_yy[i] * fe_0 + tr_y_zzzz_xyy[i] * pa_x[i];

        tr_y_xzzzz_xyz[i] = tr_y_zzzz_yz[i] * fe_0 + tr_y_zzzz_xyz[i] * pa_x[i];

        tr_y_xzzzz_xzz[i] = tr_y_zzzz_zz[i] * fe_0 + tr_y_zzzz_xzz[i] * pa_x[i];

        tr_y_xzzzz_yyy[i] = tr_y_zzzz_yyy[i] * pa_x[i];

        tr_y_xzzzz_yyz[i] = tr_y_zzzz_yyz[i] * pa_x[i];

        tr_y_xzzzz_yzz[i] = tr_y_zzzz_yzz[i] * pa_x[i];

        tr_y_xzzzz_zzz[i] = tr_y_zzzz_zzz[i] * pa_x[i];
    }

    // Set up 360-370 components of targeted buffer : HF

    auto tr_y_yyyyy_xxx = pbuffer.data(idx_dip_hf + 360);

    auto tr_y_yyyyy_xxy = pbuffer.data(idx_dip_hf + 361);

    auto tr_y_yyyyy_xxz = pbuffer.data(idx_dip_hf + 362);

    auto tr_y_yyyyy_xyy = pbuffer.data(idx_dip_hf + 363);

    auto tr_y_yyyyy_xyz = pbuffer.data(idx_dip_hf + 364);

    auto tr_y_yyyyy_xzz = pbuffer.data(idx_dip_hf + 365);

    auto tr_y_yyyyy_yyy = pbuffer.data(idx_dip_hf + 366);

    auto tr_y_yyyyy_yyz = pbuffer.data(idx_dip_hf + 367);

    auto tr_y_yyyyy_yzz = pbuffer.data(idx_dip_hf + 368);

    auto tr_y_yyyyy_zzz = pbuffer.data(idx_dip_hf + 369);

#pragma omp simd aligned(pa_y,               \
                             tr_y_yyy_xxx,   \
                             tr_y_yyy_xxy,   \
                             tr_y_yyy_xxz,   \
                             tr_y_yyy_xyy,   \
                             tr_y_yyy_xyz,   \
                             tr_y_yyy_xzz,   \
                             tr_y_yyy_yyy,   \
                             tr_y_yyy_yyz,   \
                             tr_y_yyy_yzz,   \
                             tr_y_yyy_zzz,   \
                             tr_y_yyyy_xx,   \
                             tr_y_yyyy_xxx,  \
                             tr_y_yyyy_xxy,  \
                             tr_y_yyyy_xxz,  \
                             tr_y_yyyy_xy,   \
                             tr_y_yyyy_xyy,  \
                             tr_y_yyyy_xyz,  \
                             tr_y_yyyy_xz,   \
                             tr_y_yyyy_xzz,  \
                             tr_y_yyyy_yy,   \
                             tr_y_yyyy_yyy,  \
                             tr_y_yyyy_yyz,  \
                             tr_y_yyyy_yz,   \
                             tr_y_yyyy_yzz,  \
                             tr_y_yyyy_zz,   \
                             tr_y_yyyy_zzz,  \
                             tr_y_yyyyy_xxx, \
                             tr_y_yyyyy_xxy, \
                             tr_y_yyyyy_xxz, \
                             tr_y_yyyyy_xyy, \
                             tr_y_yyyyy_xyz, \
                             tr_y_yyyyy_xzz, \
                             tr_y_yyyyy_yyy, \
                             tr_y_yyyyy_yyz, \
                             tr_y_yyyyy_yzz, \
                             tr_y_yyyyy_zzz, \
                             ts_yyyy_xxx,    \
                             ts_yyyy_xxy,    \
                             ts_yyyy_xxz,    \
                             ts_yyyy_xyy,    \
                             ts_yyyy_xyz,    \
                             ts_yyyy_xzz,    \
                             ts_yyyy_yyy,    \
                             ts_yyyy_yyz,    \
                             ts_yyyy_yzz,    \
                             ts_yyyy_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyy_xxx[i] = 4.0 * tr_y_yyy_xxx[i] * fe_0 + ts_yyyy_xxx[i] * fe_0 + tr_y_yyyy_xxx[i] * pa_y[i];

        tr_y_yyyyy_xxy[i] = 4.0 * tr_y_yyy_xxy[i] * fe_0 + tr_y_yyyy_xx[i] * fe_0 + ts_yyyy_xxy[i] * fe_0 + tr_y_yyyy_xxy[i] * pa_y[i];

        tr_y_yyyyy_xxz[i] = 4.0 * tr_y_yyy_xxz[i] * fe_0 + ts_yyyy_xxz[i] * fe_0 + tr_y_yyyy_xxz[i] * pa_y[i];

        tr_y_yyyyy_xyy[i] = 4.0 * tr_y_yyy_xyy[i] * fe_0 + 2.0 * tr_y_yyyy_xy[i] * fe_0 + ts_yyyy_xyy[i] * fe_0 + tr_y_yyyy_xyy[i] * pa_y[i];

        tr_y_yyyyy_xyz[i] = 4.0 * tr_y_yyy_xyz[i] * fe_0 + tr_y_yyyy_xz[i] * fe_0 + ts_yyyy_xyz[i] * fe_0 + tr_y_yyyy_xyz[i] * pa_y[i];

        tr_y_yyyyy_xzz[i] = 4.0 * tr_y_yyy_xzz[i] * fe_0 + ts_yyyy_xzz[i] * fe_0 + tr_y_yyyy_xzz[i] * pa_y[i];

        tr_y_yyyyy_yyy[i] = 4.0 * tr_y_yyy_yyy[i] * fe_0 + 3.0 * tr_y_yyyy_yy[i] * fe_0 + ts_yyyy_yyy[i] * fe_0 + tr_y_yyyy_yyy[i] * pa_y[i];

        tr_y_yyyyy_yyz[i] = 4.0 * tr_y_yyy_yyz[i] * fe_0 + 2.0 * tr_y_yyyy_yz[i] * fe_0 + ts_yyyy_yyz[i] * fe_0 + tr_y_yyyy_yyz[i] * pa_y[i];

        tr_y_yyyyy_yzz[i] = 4.0 * tr_y_yyy_yzz[i] * fe_0 + tr_y_yyyy_zz[i] * fe_0 + ts_yyyy_yzz[i] * fe_0 + tr_y_yyyy_yzz[i] * pa_y[i];

        tr_y_yyyyy_zzz[i] = 4.0 * tr_y_yyy_zzz[i] * fe_0 + ts_yyyy_zzz[i] * fe_0 + tr_y_yyyy_zzz[i] * pa_y[i];
    }

    // Set up 370-380 components of targeted buffer : HF

    auto tr_y_yyyyz_xxx = pbuffer.data(idx_dip_hf + 370);

    auto tr_y_yyyyz_xxy = pbuffer.data(idx_dip_hf + 371);

    auto tr_y_yyyyz_xxz = pbuffer.data(idx_dip_hf + 372);

    auto tr_y_yyyyz_xyy = pbuffer.data(idx_dip_hf + 373);

    auto tr_y_yyyyz_xyz = pbuffer.data(idx_dip_hf + 374);

    auto tr_y_yyyyz_xzz = pbuffer.data(idx_dip_hf + 375);

    auto tr_y_yyyyz_yyy = pbuffer.data(idx_dip_hf + 376);

    auto tr_y_yyyyz_yyz = pbuffer.data(idx_dip_hf + 377);

    auto tr_y_yyyyz_yzz = pbuffer.data(idx_dip_hf + 378);

    auto tr_y_yyyyz_zzz = pbuffer.data(idx_dip_hf + 379);

#pragma omp simd aligned(pa_z,               \
                             tr_y_yyyy_xx,   \
                             tr_y_yyyy_xxx,  \
                             tr_y_yyyy_xxy,  \
                             tr_y_yyyy_xxz,  \
                             tr_y_yyyy_xy,   \
                             tr_y_yyyy_xyy,  \
                             tr_y_yyyy_xyz,  \
                             tr_y_yyyy_xz,   \
                             tr_y_yyyy_xzz,  \
                             tr_y_yyyy_yy,   \
                             tr_y_yyyy_yyy,  \
                             tr_y_yyyy_yyz,  \
                             tr_y_yyyy_yz,   \
                             tr_y_yyyy_yzz,  \
                             tr_y_yyyy_zz,   \
                             tr_y_yyyy_zzz,  \
                             tr_y_yyyyz_xxx, \
                             tr_y_yyyyz_xxy, \
                             tr_y_yyyyz_xxz, \
                             tr_y_yyyyz_xyy, \
                             tr_y_yyyyz_xyz, \
                             tr_y_yyyyz_xzz, \
                             tr_y_yyyyz_yyy, \
                             tr_y_yyyyz_yyz, \
                             tr_y_yyyyz_yzz, \
                             tr_y_yyyyz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyyz_xxx[i] = tr_y_yyyy_xxx[i] * pa_z[i];

        tr_y_yyyyz_xxy[i] = tr_y_yyyy_xxy[i] * pa_z[i];

        tr_y_yyyyz_xxz[i] = tr_y_yyyy_xx[i] * fe_0 + tr_y_yyyy_xxz[i] * pa_z[i];

        tr_y_yyyyz_xyy[i] = tr_y_yyyy_xyy[i] * pa_z[i];

        tr_y_yyyyz_xyz[i] = tr_y_yyyy_xy[i] * fe_0 + tr_y_yyyy_xyz[i] * pa_z[i];

        tr_y_yyyyz_xzz[i] = 2.0 * tr_y_yyyy_xz[i] * fe_0 + tr_y_yyyy_xzz[i] * pa_z[i];

        tr_y_yyyyz_yyy[i] = tr_y_yyyy_yyy[i] * pa_z[i];

        tr_y_yyyyz_yyz[i] = tr_y_yyyy_yy[i] * fe_0 + tr_y_yyyy_yyz[i] * pa_z[i];

        tr_y_yyyyz_yzz[i] = 2.0 * tr_y_yyyy_yz[i] * fe_0 + tr_y_yyyy_yzz[i] * pa_z[i];

        tr_y_yyyyz_zzz[i] = 3.0 * tr_y_yyyy_zz[i] * fe_0 + tr_y_yyyy_zzz[i] * pa_z[i];
    }

    // Set up 380-390 components of targeted buffer : HF

    auto tr_y_yyyzz_xxx = pbuffer.data(idx_dip_hf + 380);

    auto tr_y_yyyzz_xxy = pbuffer.data(idx_dip_hf + 381);

    auto tr_y_yyyzz_xxz = pbuffer.data(idx_dip_hf + 382);

    auto tr_y_yyyzz_xyy = pbuffer.data(idx_dip_hf + 383);

    auto tr_y_yyyzz_xyz = pbuffer.data(idx_dip_hf + 384);

    auto tr_y_yyyzz_xzz = pbuffer.data(idx_dip_hf + 385);

    auto tr_y_yyyzz_yyy = pbuffer.data(idx_dip_hf + 386);

    auto tr_y_yyyzz_yyz = pbuffer.data(idx_dip_hf + 387);

    auto tr_y_yyyzz_yzz = pbuffer.data(idx_dip_hf + 388);

    auto tr_y_yyyzz_zzz = pbuffer.data(idx_dip_hf + 389);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_y_yyy_xxx,   \
                             tr_y_yyy_xxy,   \
                             tr_y_yyy_xyy,   \
                             tr_y_yyy_xyz,   \
                             tr_y_yyy_yyy,   \
                             tr_y_yyy_yyz,   \
                             tr_y_yyy_yzz,   \
                             tr_y_yyyz_xxx,  \
                             tr_y_yyyz_xxy,  \
                             tr_y_yyyz_xy,   \
                             tr_y_yyyz_xyy,  \
                             tr_y_yyyz_xyz,  \
                             tr_y_yyyz_yy,   \
                             tr_y_yyyz_yyy,  \
                             tr_y_yyyz_yyz,  \
                             tr_y_yyyz_yz,   \
                             tr_y_yyyz_yzz,  \
                             tr_y_yyyzz_xxx, \
                             tr_y_yyyzz_xxy, \
                             tr_y_yyyzz_xxz, \
                             tr_y_yyyzz_xyy, \
                             tr_y_yyyzz_xyz, \
                             tr_y_yyyzz_xzz, \
                             tr_y_yyyzz_yyy, \
                             tr_y_yyyzz_yyz, \
                             tr_y_yyyzz_yzz, \
                             tr_y_yyyzz_zzz, \
                             tr_y_yyzz_xxz,  \
                             tr_y_yyzz_xzz,  \
                             tr_y_yyzz_zzz,  \
                             tr_y_yzz_xxz,   \
                             tr_y_yzz_xzz,   \
                             tr_y_yzz_zzz,   \
                             ts_yyzz_xxz,    \
                             ts_yyzz_xzz,    \
                             ts_yyzz_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyyzz_xxx[i] = tr_y_yyy_xxx[i] * fe_0 + tr_y_yyyz_xxx[i] * pa_z[i];

        tr_y_yyyzz_xxy[i] = tr_y_yyy_xxy[i] * fe_0 + tr_y_yyyz_xxy[i] * pa_z[i];

        tr_y_yyyzz_xxz[i] = 2.0 * tr_y_yzz_xxz[i] * fe_0 + ts_yyzz_xxz[i] * fe_0 + tr_y_yyzz_xxz[i] * pa_y[i];

        tr_y_yyyzz_xyy[i] = tr_y_yyy_xyy[i] * fe_0 + tr_y_yyyz_xyy[i] * pa_z[i];

        tr_y_yyyzz_xyz[i] = tr_y_yyy_xyz[i] * fe_0 + tr_y_yyyz_xy[i] * fe_0 + tr_y_yyyz_xyz[i] * pa_z[i];

        tr_y_yyyzz_xzz[i] = 2.0 * tr_y_yzz_xzz[i] * fe_0 + ts_yyzz_xzz[i] * fe_0 + tr_y_yyzz_xzz[i] * pa_y[i];

        tr_y_yyyzz_yyy[i] = tr_y_yyy_yyy[i] * fe_0 + tr_y_yyyz_yyy[i] * pa_z[i];

        tr_y_yyyzz_yyz[i] = tr_y_yyy_yyz[i] * fe_0 + tr_y_yyyz_yy[i] * fe_0 + tr_y_yyyz_yyz[i] * pa_z[i];

        tr_y_yyyzz_yzz[i] = tr_y_yyy_yzz[i] * fe_0 + 2.0 * tr_y_yyyz_yz[i] * fe_0 + tr_y_yyyz_yzz[i] * pa_z[i];

        tr_y_yyyzz_zzz[i] = 2.0 * tr_y_yzz_zzz[i] * fe_0 + ts_yyzz_zzz[i] * fe_0 + tr_y_yyzz_zzz[i] * pa_y[i];
    }

    // Set up 390-400 components of targeted buffer : HF

    auto tr_y_yyzzz_xxx = pbuffer.data(idx_dip_hf + 390);

    auto tr_y_yyzzz_xxy = pbuffer.data(idx_dip_hf + 391);

    auto tr_y_yyzzz_xxz = pbuffer.data(idx_dip_hf + 392);

    auto tr_y_yyzzz_xyy = pbuffer.data(idx_dip_hf + 393);

    auto tr_y_yyzzz_xyz = pbuffer.data(idx_dip_hf + 394);

    auto tr_y_yyzzz_xzz = pbuffer.data(idx_dip_hf + 395);

    auto tr_y_yyzzz_yyy = pbuffer.data(idx_dip_hf + 396);

    auto tr_y_yyzzz_yyz = pbuffer.data(idx_dip_hf + 397);

    auto tr_y_yyzzz_yzz = pbuffer.data(idx_dip_hf + 398);

    auto tr_y_yyzzz_zzz = pbuffer.data(idx_dip_hf + 399);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_y_yyz_xxx,   \
                             tr_y_yyz_xxy,   \
                             tr_y_yyz_xyy,   \
                             tr_y_yyz_xyz,   \
                             tr_y_yyz_yyy,   \
                             tr_y_yyz_yyz,   \
                             tr_y_yyz_yzz,   \
                             tr_y_yyzz_xxx,  \
                             tr_y_yyzz_xxy,  \
                             tr_y_yyzz_xy,   \
                             tr_y_yyzz_xyy,  \
                             tr_y_yyzz_xyz,  \
                             tr_y_yyzz_yy,   \
                             tr_y_yyzz_yyy,  \
                             tr_y_yyzz_yyz,  \
                             tr_y_yyzz_yz,   \
                             tr_y_yyzz_yzz,  \
                             tr_y_yyzzz_xxx, \
                             tr_y_yyzzz_xxy, \
                             tr_y_yyzzz_xxz, \
                             tr_y_yyzzz_xyy, \
                             tr_y_yyzzz_xyz, \
                             tr_y_yyzzz_xzz, \
                             tr_y_yyzzz_yyy, \
                             tr_y_yyzzz_yyz, \
                             tr_y_yyzzz_yzz, \
                             tr_y_yyzzz_zzz, \
                             tr_y_yzzz_xxz,  \
                             tr_y_yzzz_xzz,  \
                             tr_y_yzzz_zzz,  \
                             tr_y_zzz_xxz,   \
                             tr_y_zzz_xzz,   \
                             tr_y_zzz_zzz,   \
                             ts_yzzz_xxz,    \
                             ts_yzzz_xzz,    \
                             ts_yzzz_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyzzz_xxx[i] = 2.0 * tr_y_yyz_xxx[i] * fe_0 + tr_y_yyzz_xxx[i] * pa_z[i];

        tr_y_yyzzz_xxy[i] = 2.0 * tr_y_yyz_xxy[i] * fe_0 + tr_y_yyzz_xxy[i] * pa_z[i];

        tr_y_yyzzz_xxz[i] = tr_y_zzz_xxz[i] * fe_0 + ts_yzzz_xxz[i] * fe_0 + tr_y_yzzz_xxz[i] * pa_y[i];

        tr_y_yyzzz_xyy[i] = 2.0 * tr_y_yyz_xyy[i] * fe_0 + tr_y_yyzz_xyy[i] * pa_z[i];

        tr_y_yyzzz_xyz[i] = 2.0 * tr_y_yyz_xyz[i] * fe_0 + tr_y_yyzz_xy[i] * fe_0 + tr_y_yyzz_xyz[i] * pa_z[i];

        tr_y_yyzzz_xzz[i] = tr_y_zzz_xzz[i] * fe_0 + ts_yzzz_xzz[i] * fe_0 + tr_y_yzzz_xzz[i] * pa_y[i];

        tr_y_yyzzz_yyy[i] = 2.0 * tr_y_yyz_yyy[i] * fe_0 + tr_y_yyzz_yyy[i] * pa_z[i];

        tr_y_yyzzz_yyz[i] = 2.0 * tr_y_yyz_yyz[i] * fe_0 + tr_y_yyzz_yy[i] * fe_0 + tr_y_yyzz_yyz[i] * pa_z[i];

        tr_y_yyzzz_yzz[i] = 2.0 * tr_y_yyz_yzz[i] * fe_0 + 2.0 * tr_y_yyzz_yz[i] * fe_0 + tr_y_yyzz_yzz[i] * pa_z[i];

        tr_y_yyzzz_zzz[i] = tr_y_zzz_zzz[i] * fe_0 + ts_yzzz_zzz[i] * fe_0 + tr_y_yzzz_zzz[i] * pa_y[i];
    }

    // Set up 400-410 components of targeted buffer : HF

    auto tr_y_yzzzz_xxx = pbuffer.data(idx_dip_hf + 400);

    auto tr_y_yzzzz_xxy = pbuffer.data(idx_dip_hf + 401);

    auto tr_y_yzzzz_xxz = pbuffer.data(idx_dip_hf + 402);

    auto tr_y_yzzzz_xyy = pbuffer.data(idx_dip_hf + 403);

    auto tr_y_yzzzz_xyz = pbuffer.data(idx_dip_hf + 404);

    auto tr_y_yzzzz_xzz = pbuffer.data(idx_dip_hf + 405);

    auto tr_y_yzzzz_yyy = pbuffer.data(idx_dip_hf + 406);

    auto tr_y_yzzzz_yyz = pbuffer.data(idx_dip_hf + 407);

    auto tr_y_yzzzz_yzz = pbuffer.data(idx_dip_hf + 408);

    auto tr_y_yzzzz_zzz = pbuffer.data(idx_dip_hf + 409);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_y_yzz_xxy,   \
                             tr_y_yzz_xyy,   \
                             tr_y_yzz_yyy,   \
                             tr_y_yzzz_xxy,  \
                             tr_y_yzzz_xyy,  \
                             tr_y_yzzz_yyy,  \
                             tr_y_yzzzz_xxx, \
                             tr_y_yzzzz_xxy, \
                             tr_y_yzzzz_xxz, \
                             tr_y_yzzzz_xyy, \
                             tr_y_yzzzz_xyz, \
                             tr_y_yzzzz_xzz, \
                             tr_y_yzzzz_yyy, \
                             tr_y_yzzzz_yyz, \
                             tr_y_yzzzz_yzz, \
                             tr_y_yzzzz_zzz, \
                             tr_y_zzzz_xxx,  \
                             tr_y_zzzz_xxz,  \
                             tr_y_zzzz_xyz,  \
                             tr_y_zzzz_xz,   \
                             tr_y_zzzz_xzz,  \
                             tr_y_zzzz_yyz,  \
                             tr_y_zzzz_yz,   \
                             tr_y_zzzz_yzz,  \
                             tr_y_zzzz_zz,   \
                             tr_y_zzzz_zzz,  \
                             ts_zzzz_xxx,    \
                             ts_zzzz_xxz,    \
                             ts_zzzz_xyz,    \
                             ts_zzzz_xzz,    \
                             ts_zzzz_yyz,    \
                             ts_zzzz_yzz,    \
                             ts_zzzz_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzzzz_xxx[i] = ts_zzzz_xxx[i] * fe_0 + tr_y_zzzz_xxx[i] * pa_y[i];

        tr_y_yzzzz_xxy[i] = 3.0 * tr_y_yzz_xxy[i] * fe_0 + tr_y_yzzz_xxy[i] * pa_z[i];

        tr_y_yzzzz_xxz[i] = ts_zzzz_xxz[i] * fe_0 + tr_y_zzzz_xxz[i] * pa_y[i];

        tr_y_yzzzz_xyy[i] = 3.0 * tr_y_yzz_xyy[i] * fe_0 + tr_y_yzzz_xyy[i] * pa_z[i];

        tr_y_yzzzz_xyz[i] = tr_y_zzzz_xz[i] * fe_0 + ts_zzzz_xyz[i] * fe_0 + tr_y_zzzz_xyz[i] * pa_y[i];

        tr_y_yzzzz_xzz[i] = ts_zzzz_xzz[i] * fe_0 + tr_y_zzzz_xzz[i] * pa_y[i];

        tr_y_yzzzz_yyy[i] = 3.0 * tr_y_yzz_yyy[i] * fe_0 + tr_y_yzzz_yyy[i] * pa_z[i];

        tr_y_yzzzz_yyz[i] = 2.0 * tr_y_zzzz_yz[i] * fe_0 + ts_zzzz_yyz[i] * fe_0 + tr_y_zzzz_yyz[i] * pa_y[i];

        tr_y_yzzzz_yzz[i] = tr_y_zzzz_zz[i] * fe_0 + ts_zzzz_yzz[i] * fe_0 + tr_y_zzzz_yzz[i] * pa_y[i];

        tr_y_yzzzz_zzz[i] = ts_zzzz_zzz[i] * fe_0 + tr_y_zzzz_zzz[i] * pa_y[i];
    }

    // Set up 410-420 components of targeted buffer : HF

    auto tr_y_zzzzz_xxx = pbuffer.data(idx_dip_hf + 410);

    auto tr_y_zzzzz_xxy = pbuffer.data(idx_dip_hf + 411);

    auto tr_y_zzzzz_xxz = pbuffer.data(idx_dip_hf + 412);

    auto tr_y_zzzzz_xyy = pbuffer.data(idx_dip_hf + 413);

    auto tr_y_zzzzz_xyz = pbuffer.data(idx_dip_hf + 414);

    auto tr_y_zzzzz_xzz = pbuffer.data(idx_dip_hf + 415);

    auto tr_y_zzzzz_yyy = pbuffer.data(idx_dip_hf + 416);

    auto tr_y_zzzzz_yyz = pbuffer.data(idx_dip_hf + 417);

    auto tr_y_zzzzz_yzz = pbuffer.data(idx_dip_hf + 418);

    auto tr_y_zzzzz_zzz = pbuffer.data(idx_dip_hf + 419);

#pragma omp simd aligned(pa_z,               \
                             tr_y_zzz_xxx,   \
                             tr_y_zzz_xxy,   \
                             tr_y_zzz_xxz,   \
                             tr_y_zzz_xyy,   \
                             tr_y_zzz_xyz,   \
                             tr_y_zzz_xzz,   \
                             tr_y_zzz_yyy,   \
                             tr_y_zzz_yyz,   \
                             tr_y_zzz_yzz,   \
                             tr_y_zzz_zzz,   \
                             tr_y_zzzz_xx,   \
                             tr_y_zzzz_xxx,  \
                             tr_y_zzzz_xxy,  \
                             tr_y_zzzz_xxz,  \
                             tr_y_zzzz_xy,   \
                             tr_y_zzzz_xyy,  \
                             tr_y_zzzz_xyz,  \
                             tr_y_zzzz_xz,   \
                             tr_y_zzzz_xzz,  \
                             tr_y_zzzz_yy,   \
                             tr_y_zzzz_yyy,  \
                             tr_y_zzzz_yyz,  \
                             tr_y_zzzz_yz,   \
                             tr_y_zzzz_yzz,  \
                             tr_y_zzzz_zz,   \
                             tr_y_zzzz_zzz,  \
                             tr_y_zzzzz_xxx, \
                             tr_y_zzzzz_xxy, \
                             tr_y_zzzzz_xxz, \
                             tr_y_zzzzz_xyy, \
                             tr_y_zzzzz_xyz, \
                             tr_y_zzzzz_xzz, \
                             tr_y_zzzzz_yyy, \
                             tr_y_zzzzz_yyz, \
                             tr_y_zzzzz_yzz, \
                             tr_y_zzzzz_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzzzz_xxx[i] = 4.0 * tr_y_zzz_xxx[i] * fe_0 + tr_y_zzzz_xxx[i] * pa_z[i];

        tr_y_zzzzz_xxy[i] = 4.0 * tr_y_zzz_xxy[i] * fe_0 + tr_y_zzzz_xxy[i] * pa_z[i];

        tr_y_zzzzz_xxz[i] = 4.0 * tr_y_zzz_xxz[i] * fe_0 + tr_y_zzzz_xx[i] * fe_0 + tr_y_zzzz_xxz[i] * pa_z[i];

        tr_y_zzzzz_xyy[i] = 4.0 * tr_y_zzz_xyy[i] * fe_0 + tr_y_zzzz_xyy[i] * pa_z[i];

        tr_y_zzzzz_xyz[i] = 4.0 * tr_y_zzz_xyz[i] * fe_0 + tr_y_zzzz_xy[i] * fe_0 + tr_y_zzzz_xyz[i] * pa_z[i];

        tr_y_zzzzz_xzz[i] = 4.0 * tr_y_zzz_xzz[i] * fe_0 + 2.0 * tr_y_zzzz_xz[i] * fe_0 + tr_y_zzzz_xzz[i] * pa_z[i];

        tr_y_zzzzz_yyy[i] = 4.0 * tr_y_zzz_yyy[i] * fe_0 + tr_y_zzzz_yyy[i] * pa_z[i];

        tr_y_zzzzz_yyz[i] = 4.0 * tr_y_zzz_yyz[i] * fe_0 + tr_y_zzzz_yy[i] * fe_0 + tr_y_zzzz_yyz[i] * pa_z[i];

        tr_y_zzzzz_yzz[i] = 4.0 * tr_y_zzz_yzz[i] * fe_0 + 2.0 * tr_y_zzzz_yz[i] * fe_0 + tr_y_zzzz_yzz[i] * pa_z[i];

        tr_y_zzzzz_zzz[i] = 4.0 * tr_y_zzz_zzz[i] * fe_0 + 3.0 * tr_y_zzzz_zz[i] * fe_0 + tr_y_zzzz_zzz[i] * pa_z[i];
    }

    // Set up 420-430 components of targeted buffer : HF

    auto tr_z_xxxxx_xxx = pbuffer.data(idx_dip_hf + 420);

    auto tr_z_xxxxx_xxy = pbuffer.data(idx_dip_hf + 421);

    auto tr_z_xxxxx_xxz = pbuffer.data(idx_dip_hf + 422);

    auto tr_z_xxxxx_xyy = pbuffer.data(idx_dip_hf + 423);

    auto tr_z_xxxxx_xyz = pbuffer.data(idx_dip_hf + 424);

    auto tr_z_xxxxx_xzz = pbuffer.data(idx_dip_hf + 425);

    auto tr_z_xxxxx_yyy = pbuffer.data(idx_dip_hf + 426);

    auto tr_z_xxxxx_yyz = pbuffer.data(idx_dip_hf + 427);

    auto tr_z_xxxxx_yzz = pbuffer.data(idx_dip_hf + 428);

    auto tr_z_xxxxx_zzz = pbuffer.data(idx_dip_hf + 429);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xxx_xxx,   \
                             tr_z_xxx_xxy,   \
                             tr_z_xxx_xxz,   \
                             tr_z_xxx_xyy,   \
                             tr_z_xxx_xyz,   \
                             tr_z_xxx_xzz,   \
                             tr_z_xxx_yyy,   \
                             tr_z_xxx_yyz,   \
                             tr_z_xxx_yzz,   \
                             tr_z_xxx_zzz,   \
                             tr_z_xxxx_xx,   \
                             tr_z_xxxx_xxx,  \
                             tr_z_xxxx_xxy,  \
                             tr_z_xxxx_xxz,  \
                             tr_z_xxxx_xy,   \
                             tr_z_xxxx_xyy,  \
                             tr_z_xxxx_xyz,  \
                             tr_z_xxxx_xz,   \
                             tr_z_xxxx_xzz,  \
                             tr_z_xxxx_yy,   \
                             tr_z_xxxx_yyy,  \
                             tr_z_xxxx_yyz,  \
                             tr_z_xxxx_yz,   \
                             tr_z_xxxx_yzz,  \
                             tr_z_xxxx_zz,   \
                             tr_z_xxxx_zzz,  \
                             tr_z_xxxxx_xxx, \
                             tr_z_xxxxx_xxy, \
                             tr_z_xxxxx_xxz, \
                             tr_z_xxxxx_xyy, \
                             tr_z_xxxxx_xyz, \
                             tr_z_xxxxx_xzz, \
                             tr_z_xxxxx_yyy, \
                             tr_z_xxxxx_yyz, \
                             tr_z_xxxxx_yzz, \
                             tr_z_xxxxx_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxx_xxx[i] = 4.0 * tr_z_xxx_xxx[i] * fe_0 + 3.0 * tr_z_xxxx_xx[i] * fe_0 + tr_z_xxxx_xxx[i] * pa_x[i];

        tr_z_xxxxx_xxy[i] = 4.0 * tr_z_xxx_xxy[i] * fe_0 + 2.0 * tr_z_xxxx_xy[i] * fe_0 + tr_z_xxxx_xxy[i] * pa_x[i];

        tr_z_xxxxx_xxz[i] = 4.0 * tr_z_xxx_xxz[i] * fe_0 + 2.0 * tr_z_xxxx_xz[i] * fe_0 + tr_z_xxxx_xxz[i] * pa_x[i];

        tr_z_xxxxx_xyy[i] = 4.0 * tr_z_xxx_xyy[i] * fe_0 + tr_z_xxxx_yy[i] * fe_0 + tr_z_xxxx_xyy[i] * pa_x[i];

        tr_z_xxxxx_xyz[i] = 4.0 * tr_z_xxx_xyz[i] * fe_0 + tr_z_xxxx_yz[i] * fe_0 + tr_z_xxxx_xyz[i] * pa_x[i];

        tr_z_xxxxx_xzz[i] = 4.0 * tr_z_xxx_xzz[i] * fe_0 + tr_z_xxxx_zz[i] * fe_0 + tr_z_xxxx_xzz[i] * pa_x[i];

        tr_z_xxxxx_yyy[i] = 4.0 * tr_z_xxx_yyy[i] * fe_0 + tr_z_xxxx_yyy[i] * pa_x[i];

        tr_z_xxxxx_yyz[i] = 4.0 * tr_z_xxx_yyz[i] * fe_0 + tr_z_xxxx_yyz[i] * pa_x[i];

        tr_z_xxxxx_yzz[i] = 4.0 * tr_z_xxx_yzz[i] * fe_0 + tr_z_xxxx_yzz[i] * pa_x[i];

        tr_z_xxxxx_zzz[i] = 4.0 * tr_z_xxx_zzz[i] * fe_0 + tr_z_xxxx_zzz[i] * pa_x[i];
    }

    // Set up 430-440 components of targeted buffer : HF

    auto tr_z_xxxxy_xxx = pbuffer.data(idx_dip_hf + 430);

    auto tr_z_xxxxy_xxy = pbuffer.data(idx_dip_hf + 431);

    auto tr_z_xxxxy_xxz = pbuffer.data(idx_dip_hf + 432);

    auto tr_z_xxxxy_xyy = pbuffer.data(idx_dip_hf + 433);

    auto tr_z_xxxxy_xyz = pbuffer.data(idx_dip_hf + 434);

    auto tr_z_xxxxy_xzz = pbuffer.data(idx_dip_hf + 435);

    auto tr_z_xxxxy_yyy = pbuffer.data(idx_dip_hf + 436);

    auto tr_z_xxxxy_yyz = pbuffer.data(idx_dip_hf + 437);

    auto tr_z_xxxxy_yzz = pbuffer.data(idx_dip_hf + 438);

    auto tr_z_xxxxy_zzz = pbuffer.data(idx_dip_hf + 439);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xxxx_xx,   \
                             tr_z_xxxx_xxx,  \
                             tr_z_xxxx_xxy,  \
                             tr_z_xxxx_xxz,  \
                             tr_z_xxxx_xy,   \
                             tr_z_xxxx_xyy,  \
                             tr_z_xxxx_xyz,  \
                             tr_z_xxxx_xz,   \
                             tr_z_xxxx_xzz,  \
                             tr_z_xxxx_zzz,  \
                             tr_z_xxxxy_xxx, \
                             tr_z_xxxxy_xxy, \
                             tr_z_xxxxy_xxz, \
                             tr_z_xxxxy_xyy, \
                             tr_z_xxxxy_xyz, \
                             tr_z_xxxxy_xzz, \
                             tr_z_xxxxy_yyy, \
                             tr_z_xxxxy_yyz, \
                             tr_z_xxxxy_yzz, \
                             tr_z_xxxxy_zzz, \
                             tr_z_xxxy_yyy,  \
                             tr_z_xxxy_yyz,  \
                             tr_z_xxxy_yzz,  \
                             tr_z_xxy_yyy,   \
                             tr_z_xxy_yyz,   \
                             tr_z_xxy_yzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxy_xxx[i] = tr_z_xxxx_xxx[i] * pa_y[i];

        tr_z_xxxxy_xxy[i] = tr_z_xxxx_xx[i] * fe_0 + tr_z_xxxx_xxy[i] * pa_y[i];

        tr_z_xxxxy_xxz[i] = tr_z_xxxx_xxz[i] * pa_y[i];

        tr_z_xxxxy_xyy[i] = 2.0 * tr_z_xxxx_xy[i] * fe_0 + tr_z_xxxx_xyy[i] * pa_y[i];

        tr_z_xxxxy_xyz[i] = tr_z_xxxx_xz[i] * fe_0 + tr_z_xxxx_xyz[i] * pa_y[i];

        tr_z_xxxxy_xzz[i] = tr_z_xxxx_xzz[i] * pa_y[i];

        tr_z_xxxxy_yyy[i] = 3.0 * tr_z_xxy_yyy[i] * fe_0 + tr_z_xxxy_yyy[i] * pa_x[i];

        tr_z_xxxxy_yyz[i] = 3.0 * tr_z_xxy_yyz[i] * fe_0 + tr_z_xxxy_yyz[i] * pa_x[i];

        tr_z_xxxxy_yzz[i] = 3.0 * tr_z_xxy_yzz[i] * fe_0 + tr_z_xxxy_yzz[i] * pa_x[i];

        tr_z_xxxxy_zzz[i] = tr_z_xxxx_zzz[i] * pa_y[i];
    }

    // Set up 440-450 components of targeted buffer : HF

    auto tr_z_xxxxz_xxx = pbuffer.data(idx_dip_hf + 440);

    auto tr_z_xxxxz_xxy = pbuffer.data(idx_dip_hf + 441);

    auto tr_z_xxxxz_xxz = pbuffer.data(idx_dip_hf + 442);

    auto tr_z_xxxxz_xyy = pbuffer.data(idx_dip_hf + 443);

    auto tr_z_xxxxz_xyz = pbuffer.data(idx_dip_hf + 444);

    auto tr_z_xxxxz_xzz = pbuffer.data(idx_dip_hf + 445);

    auto tr_z_xxxxz_yyy = pbuffer.data(idx_dip_hf + 446);

    auto tr_z_xxxxz_yyz = pbuffer.data(idx_dip_hf + 447);

    auto tr_z_xxxxz_yzz = pbuffer.data(idx_dip_hf + 448);

    auto tr_z_xxxxz_zzz = pbuffer.data(idx_dip_hf + 449);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_z_xxxx_xxx,  \
                             tr_z_xxxx_xxy,  \
                             tr_z_xxxx_xyy,  \
                             tr_z_xxxxz_xxx, \
                             tr_z_xxxxz_xxy, \
                             tr_z_xxxxz_xxz, \
                             tr_z_xxxxz_xyy, \
                             tr_z_xxxxz_xyz, \
                             tr_z_xxxxz_xzz, \
                             tr_z_xxxxz_yyy, \
                             tr_z_xxxxz_yyz, \
                             tr_z_xxxxz_yzz, \
                             tr_z_xxxxz_zzz, \
                             tr_z_xxxz_xxz,  \
                             tr_z_xxxz_xyz,  \
                             tr_z_xxxz_xz,   \
                             tr_z_xxxz_xzz,  \
                             tr_z_xxxz_yyy,  \
                             tr_z_xxxz_yyz,  \
                             tr_z_xxxz_yz,   \
                             tr_z_xxxz_yzz,  \
                             tr_z_xxxz_zz,   \
                             tr_z_xxxz_zzz,  \
                             tr_z_xxz_xxz,   \
                             tr_z_xxz_xyz,   \
                             tr_z_xxz_xzz,   \
                             tr_z_xxz_yyy,   \
                             tr_z_xxz_yyz,   \
                             tr_z_xxz_yzz,   \
                             tr_z_xxz_zzz,   \
                             ts_xxxx_xxx,    \
                             ts_xxxx_xxy,    \
                             ts_xxxx_xyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxxz_xxx[i] = ts_xxxx_xxx[i] * fe_0 + tr_z_xxxx_xxx[i] * pa_z[i];

        tr_z_xxxxz_xxy[i] = ts_xxxx_xxy[i] * fe_0 + tr_z_xxxx_xxy[i] * pa_z[i];

        tr_z_xxxxz_xxz[i] = 3.0 * tr_z_xxz_xxz[i] * fe_0 + 2.0 * tr_z_xxxz_xz[i] * fe_0 + tr_z_xxxz_xxz[i] * pa_x[i];

        tr_z_xxxxz_xyy[i] = ts_xxxx_xyy[i] * fe_0 + tr_z_xxxx_xyy[i] * pa_z[i];

        tr_z_xxxxz_xyz[i] = 3.0 * tr_z_xxz_xyz[i] * fe_0 + tr_z_xxxz_yz[i] * fe_0 + tr_z_xxxz_xyz[i] * pa_x[i];

        tr_z_xxxxz_xzz[i] = 3.0 * tr_z_xxz_xzz[i] * fe_0 + tr_z_xxxz_zz[i] * fe_0 + tr_z_xxxz_xzz[i] * pa_x[i];

        tr_z_xxxxz_yyy[i] = 3.0 * tr_z_xxz_yyy[i] * fe_0 + tr_z_xxxz_yyy[i] * pa_x[i];

        tr_z_xxxxz_yyz[i] = 3.0 * tr_z_xxz_yyz[i] * fe_0 + tr_z_xxxz_yyz[i] * pa_x[i];

        tr_z_xxxxz_yzz[i] = 3.0 * tr_z_xxz_yzz[i] * fe_0 + tr_z_xxxz_yzz[i] * pa_x[i];

        tr_z_xxxxz_zzz[i] = 3.0 * tr_z_xxz_zzz[i] * fe_0 + tr_z_xxxz_zzz[i] * pa_x[i];
    }

    // Set up 450-460 components of targeted buffer : HF

    auto tr_z_xxxyy_xxx = pbuffer.data(idx_dip_hf + 450);

    auto tr_z_xxxyy_xxy = pbuffer.data(idx_dip_hf + 451);

    auto tr_z_xxxyy_xxz = pbuffer.data(idx_dip_hf + 452);

    auto tr_z_xxxyy_xyy = pbuffer.data(idx_dip_hf + 453);

    auto tr_z_xxxyy_xyz = pbuffer.data(idx_dip_hf + 454);

    auto tr_z_xxxyy_xzz = pbuffer.data(idx_dip_hf + 455);

    auto tr_z_xxxyy_yyy = pbuffer.data(idx_dip_hf + 456);

    auto tr_z_xxxyy_yyz = pbuffer.data(idx_dip_hf + 457);

    auto tr_z_xxxyy_yzz = pbuffer.data(idx_dip_hf + 458);

    auto tr_z_xxxyy_zzz = pbuffer.data(idx_dip_hf + 459);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xxx_xxx,   \
                             tr_z_xxx_xxz,   \
                             tr_z_xxx_xzz,   \
                             tr_z_xxxy_xxx,  \
                             tr_z_xxxy_xxz,  \
                             tr_z_xxxy_xzz,  \
                             tr_z_xxxyy_xxx, \
                             tr_z_xxxyy_xxy, \
                             tr_z_xxxyy_xxz, \
                             tr_z_xxxyy_xyy, \
                             tr_z_xxxyy_xyz, \
                             tr_z_xxxyy_xzz, \
                             tr_z_xxxyy_yyy, \
                             tr_z_xxxyy_yyz, \
                             tr_z_xxxyy_yzz, \
                             tr_z_xxxyy_zzz, \
                             tr_z_xxyy_xxy,  \
                             tr_z_xxyy_xy,   \
                             tr_z_xxyy_xyy,  \
                             tr_z_xxyy_xyz,  \
                             tr_z_xxyy_yy,   \
                             tr_z_xxyy_yyy,  \
                             tr_z_xxyy_yyz,  \
                             tr_z_xxyy_yz,   \
                             tr_z_xxyy_yzz,  \
                             tr_z_xxyy_zzz,  \
                             tr_z_xyy_xxy,   \
                             tr_z_xyy_xyy,   \
                             tr_z_xyy_xyz,   \
                             tr_z_xyy_yyy,   \
                             tr_z_xyy_yyz,   \
                             tr_z_xyy_yzz,   \
                             tr_z_xyy_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyy_xxx[i] = tr_z_xxx_xxx[i] * fe_0 + tr_z_xxxy_xxx[i] * pa_y[i];

        tr_z_xxxyy_xxy[i] = 2.0 * tr_z_xyy_xxy[i] * fe_0 + 2.0 * tr_z_xxyy_xy[i] * fe_0 + tr_z_xxyy_xxy[i] * pa_x[i];

        tr_z_xxxyy_xxz[i] = tr_z_xxx_xxz[i] * fe_0 + tr_z_xxxy_xxz[i] * pa_y[i];

        tr_z_xxxyy_xyy[i] = 2.0 * tr_z_xyy_xyy[i] * fe_0 + tr_z_xxyy_yy[i] * fe_0 + tr_z_xxyy_xyy[i] * pa_x[i];

        tr_z_xxxyy_xyz[i] = 2.0 * tr_z_xyy_xyz[i] * fe_0 + tr_z_xxyy_yz[i] * fe_0 + tr_z_xxyy_xyz[i] * pa_x[i];

        tr_z_xxxyy_xzz[i] = tr_z_xxx_xzz[i] * fe_0 + tr_z_xxxy_xzz[i] * pa_y[i];

        tr_z_xxxyy_yyy[i] = 2.0 * tr_z_xyy_yyy[i] * fe_0 + tr_z_xxyy_yyy[i] * pa_x[i];

        tr_z_xxxyy_yyz[i] = 2.0 * tr_z_xyy_yyz[i] * fe_0 + tr_z_xxyy_yyz[i] * pa_x[i];

        tr_z_xxxyy_yzz[i] = 2.0 * tr_z_xyy_yzz[i] * fe_0 + tr_z_xxyy_yzz[i] * pa_x[i];

        tr_z_xxxyy_zzz[i] = 2.0 * tr_z_xyy_zzz[i] * fe_0 + tr_z_xxyy_zzz[i] * pa_x[i];
    }

    // Set up 460-470 components of targeted buffer : HF

    auto tr_z_xxxyz_xxx = pbuffer.data(idx_dip_hf + 460);

    auto tr_z_xxxyz_xxy = pbuffer.data(idx_dip_hf + 461);

    auto tr_z_xxxyz_xxz = pbuffer.data(idx_dip_hf + 462);

    auto tr_z_xxxyz_xyy = pbuffer.data(idx_dip_hf + 463);

    auto tr_z_xxxyz_xyz = pbuffer.data(idx_dip_hf + 464);

    auto tr_z_xxxyz_xzz = pbuffer.data(idx_dip_hf + 465);

    auto tr_z_xxxyz_yyy = pbuffer.data(idx_dip_hf + 466);

    auto tr_z_xxxyz_yyz = pbuffer.data(idx_dip_hf + 467);

    auto tr_z_xxxyz_yzz = pbuffer.data(idx_dip_hf + 468);

    auto tr_z_xxxyz_zzz = pbuffer.data(idx_dip_hf + 469);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xxxyz_xxx, \
                             tr_z_xxxyz_xxy, \
                             tr_z_xxxyz_xxz, \
                             tr_z_xxxyz_xyy, \
                             tr_z_xxxyz_xyz, \
                             tr_z_xxxyz_xzz, \
                             tr_z_xxxyz_yyy, \
                             tr_z_xxxyz_yyz, \
                             tr_z_xxxyz_yzz, \
                             tr_z_xxxyz_zzz, \
                             tr_z_xxxz_xx,   \
                             tr_z_xxxz_xxx,  \
                             tr_z_xxxz_xxy,  \
                             tr_z_xxxz_xxz,  \
                             tr_z_xxxz_xy,   \
                             tr_z_xxxz_xyy,  \
                             tr_z_xxxz_xyz,  \
                             tr_z_xxxz_xz,   \
                             tr_z_xxxz_xzz,  \
                             tr_z_xxxz_zzz,  \
                             tr_z_xxyz_yyy,  \
                             tr_z_xxyz_yyz,  \
                             tr_z_xxyz_yzz,  \
                             tr_z_xyz_yyy,   \
                             tr_z_xyz_yyz,   \
                             tr_z_xyz_yzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxyz_xxx[i] = tr_z_xxxz_xxx[i] * pa_y[i];

        tr_z_xxxyz_xxy[i] = tr_z_xxxz_xx[i] * fe_0 + tr_z_xxxz_xxy[i] * pa_y[i];

        tr_z_xxxyz_xxz[i] = tr_z_xxxz_xxz[i] * pa_y[i];

        tr_z_xxxyz_xyy[i] = 2.0 * tr_z_xxxz_xy[i] * fe_0 + tr_z_xxxz_xyy[i] * pa_y[i];

        tr_z_xxxyz_xyz[i] = tr_z_xxxz_xz[i] * fe_0 + tr_z_xxxz_xyz[i] * pa_y[i];

        tr_z_xxxyz_xzz[i] = tr_z_xxxz_xzz[i] * pa_y[i];

        tr_z_xxxyz_yyy[i] = 2.0 * tr_z_xyz_yyy[i] * fe_0 + tr_z_xxyz_yyy[i] * pa_x[i];

        tr_z_xxxyz_yyz[i] = 2.0 * tr_z_xyz_yyz[i] * fe_0 + tr_z_xxyz_yyz[i] * pa_x[i];

        tr_z_xxxyz_yzz[i] = 2.0 * tr_z_xyz_yzz[i] * fe_0 + tr_z_xxyz_yzz[i] * pa_x[i];

        tr_z_xxxyz_zzz[i] = tr_z_xxxz_zzz[i] * pa_y[i];
    }

    // Set up 470-480 components of targeted buffer : HF

    auto tr_z_xxxzz_xxx = pbuffer.data(idx_dip_hf + 470);

    auto tr_z_xxxzz_xxy = pbuffer.data(idx_dip_hf + 471);

    auto tr_z_xxxzz_xxz = pbuffer.data(idx_dip_hf + 472);

    auto tr_z_xxxzz_xyy = pbuffer.data(idx_dip_hf + 473);

    auto tr_z_xxxzz_xyz = pbuffer.data(idx_dip_hf + 474);

    auto tr_z_xxxzz_xzz = pbuffer.data(idx_dip_hf + 475);

    auto tr_z_xxxzz_yyy = pbuffer.data(idx_dip_hf + 476);

    auto tr_z_xxxzz_yyz = pbuffer.data(idx_dip_hf + 477);

    auto tr_z_xxxzz_yzz = pbuffer.data(idx_dip_hf + 478);

    auto tr_z_xxxzz_zzz = pbuffer.data(idx_dip_hf + 479);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xxxzz_xxx, \
                             tr_z_xxxzz_xxy, \
                             tr_z_xxxzz_xxz, \
                             tr_z_xxxzz_xyy, \
                             tr_z_xxxzz_xyz, \
                             tr_z_xxxzz_xzz, \
                             tr_z_xxxzz_yyy, \
                             tr_z_xxxzz_yyz, \
                             tr_z_xxxzz_yzz, \
                             tr_z_xxxzz_zzz, \
                             tr_z_xxzz_xx,   \
                             tr_z_xxzz_xxx,  \
                             tr_z_xxzz_xxy,  \
                             tr_z_xxzz_xxz,  \
                             tr_z_xxzz_xy,   \
                             tr_z_xxzz_xyy,  \
                             tr_z_xxzz_xyz,  \
                             tr_z_xxzz_xz,   \
                             tr_z_xxzz_xzz,  \
                             tr_z_xxzz_yy,   \
                             tr_z_xxzz_yyy,  \
                             tr_z_xxzz_yyz,  \
                             tr_z_xxzz_yz,   \
                             tr_z_xxzz_yzz,  \
                             tr_z_xxzz_zz,   \
                             tr_z_xxzz_zzz,  \
                             tr_z_xzz_xxx,   \
                             tr_z_xzz_xxy,   \
                             tr_z_xzz_xxz,   \
                             tr_z_xzz_xyy,   \
                             tr_z_xzz_xyz,   \
                             tr_z_xzz_xzz,   \
                             tr_z_xzz_yyy,   \
                             tr_z_xzz_yyz,   \
                             tr_z_xzz_yzz,   \
                             tr_z_xzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxxzz_xxx[i] = 2.0 * tr_z_xzz_xxx[i] * fe_0 + 3.0 * tr_z_xxzz_xx[i] * fe_0 + tr_z_xxzz_xxx[i] * pa_x[i];

        tr_z_xxxzz_xxy[i] = 2.0 * tr_z_xzz_xxy[i] * fe_0 + 2.0 * tr_z_xxzz_xy[i] * fe_0 + tr_z_xxzz_xxy[i] * pa_x[i];

        tr_z_xxxzz_xxz[i] = 2.0 * tr_z_xzz_xxz[i] * fe_0 + 2.0 * tr_z_xxzz_xz[i] * fe_0 + tr_z_xxzz_xxz[i] * pa_x[i];

        tr_z_xxxzz_xyy[i] = 2.0 * tr_z_xzz_xyy[i] * fe_0 + tr_z_xxzz_yy[i] * fe_0 + tr_z_xxzz_xyy[i] * pa_x[i];

        tr_z_xxxzz_xyz[i] = 2.0 * tr_z_xzz_xyz[i] * fe_0 + tr_z_xxzz_yz[i] * fe_0 + tr_z_xxzz_xyz[i] * pa_x[i];

        tr_z_xxxzz_xzz[i] = 2.0 * tr_z_xzz_xzz[i] * fe_0 + tr_z_xxzz_zz[i] * fe_0 + tr_z_xxzz_xzz[i] * pa_x[i];

        tr_z_xxxzz_yyy[i] = 2.0 * tr_z_xzz_yyy[i] * fe_0 + tr_z_xxzz_yyy[i] * pa_x[i];

        tr_z_xxxzz_yyz[i] = 2.0 * tr_z_xzz_yyz[i] * fe_0 + tr_z_xxzz_yyz[i] * pa_x[i];

        tr_z_xxxzz_yzz[i] = 2.0 * tr_z_xzz_yzz[i] * fe_0 + tr_z_xxzz_yzz[i] * pa_x[i];

        tr_z_xxxzz_zzz[i] = 2.0 * tr_z_xzz_zzz[i] * fe_0 + tr_z_xxzz_zzz[i] * pa_x[i];
    }

    // Set up 480-490 components of targeted buffer : HF

    auto tr_z_xxyyy_xxx = pbuffer.data(idx_dip_hf + 480);

    auto tr_z_xxyyy_xxy = pbuffer.data(idx_dip_hf + 481);

    auto tr_z_xxyyy_xxz = pbuffer.data(idx_dip_hf + 482);

    auto tr_z_xxyyy_xyy = pbuffer.data(idx_dip_hf + 483);

    auto tr_z_xxyyy_xyz = pbuffer.data(idx_dip_hf + 484);

    auto tr_z_xxyyy_xzz = pbuffer.data(idx_dip_hf + 485);

    auto tr_z_xxyyy_yyy = pbuffer.data(idx_dip_hf + 486);

    auto tr_z_xxyyy_yyz = pbuffer.data(idx_dip_hf + 487);

    auto tr_z_xxyyy_yzz = pbuffer.data(idx_dip_hf + 488);

    auto tr_z_xxyyy_zzz = pbuffer.data(idx_dip_hf + 489);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xxy_xxx,   \
                             tr_z_xxy_xxz,   \
                             tr_z_xxy_xzz,   \
                             tr_z_xxyy_xxx,  \
                             tr_z_xxyy_xxz,  \
                             tr_z_xxyy_xzz,  \
                             tr_z_xxyyy_xxx, \
                             tr_z_xxyyy_xxy, \
                             tr_z_xxyyy_xxz, \
                             tr_z_xxyyy_xyy, \
                             tr_z_xxyyy_xyz, \
                             tr_z_xxyyy_xzz, \
                             tr_z_xxyyy_yyy, \
                             tr_z_xxyyy_yyz, \
                             tr_z_xxyyy_yzz, \
                             tr_z_xxyyy_zzz, \
                             tr_z_xyyy_xxy,  \
                             tr_z_xyyy_xy,   \
                             tr_z_xyyy_xyy,  \
                             tr_z_xyyy_xyz,  \
                             tr_z_xyyy_yy,   \
                             tr_z_xyyy_yyy,  \
                             tr_z_xyyy_yyz,  \
                             tr_z_xyyy_yz,   \
                             tr_z_xyyy_yzz,  \
                             tr_z_xyyy_zzz,  \
                             tr_z_yyy_xxy,   \
                             tr_z_yyy_xyy,   \
                             tr_z_yyy_xyz,   \
                             tr_z_yyy_yyy,   \
                             tr_z_yyy_yyz,   \
                             tr_z_yyy_yzz,   \
                             tr_z_yyy_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyy_xxx[i] = 2.0 * tr_z_xxy_xxx[i] * fe_0 + tr_z_xxyy_xxx[i] * pa_y[i];

        tr_z_xxyyy_xxy[i] = tr_z_yyy_xxy[i] * fe_0 + 2.0 * tr_z_xyyy_xy[i] * fe_0 + tr_z_xyyy_xxy[i] * pa_x[i];

        tr_z_xxyyy_xxz[i] = 2.0 * tr_z_xxy_xxz[i] * fe_0 + tr_z_xxyy_xxz[i] * pa_y[i];

        tr_z_xxyyy_xyy[i] = tr_z_yyy_xyy[i] * fe_0 + tr_z_xyyy_yy[i] * fe_0 + tr_z_xyyy_xyy[i] * pa_x[i];

        tr_z_xxyyy_xyz[i] = tr_z_yyy_xyz[i] * fe_0 + tr_z_xyyy_yz[i] * fe_0 + tr_z_xyyy_xyz[i] * pa_x[i];

        tr_z_xxyyy_xzz[i] = 2.0 * tr_z_xxy_xzz[i] * fe_0 + tr_z_xxyy_xzz[i] * pa_y[i];

        tr_z_xxyyy_yyy[i] = tr_z_yyy_yyy[i] * fe_0 + tr_z_xyyy_yyy[i] * pa_x[i];

        tr_z_xxyyy_yyz[i] = tr_z_yyy_yyz[i] * fe_0 + tr_z_xyyy_yyz[i] * pa_x[i];

        tr_z_xxyyy_yzz[i] = tr_z_yyy_yzz[i] * fe_0 + tr_z_xyyy_yzz[i] * pa_x[i];

        tr_z_xxyyy_zzz[i] = tr_z_yyy_zzz[i] * fe_0 + tr_z_xyyy_zzz[i] * pa_x[i];
    }

    // Set up 490-500 components of targeted buffer : HF

    auto tr_z_xxyyz_xxx = pbuffer.data(idx_dip_hf + 490);

    auto tr_z_xxyyz_xxy = pbuffer.data(idx_dip_hf + 491);

    auto tr_z_xxyyz_xxz = pbuffer.data(idx_dip_hf + 492);

    auto tr_z_xxyyz_xyy = pbuffer.data(idx_dip_hf + 493);

    auto tr_z_xxyyz_xyz = pbuffer.data(idx_dip_hf + 494);

    auto tr_z_xxyyz_xzz = pbuffer.data(idx_dip_hf + 495);

    auto tr_z_xxyyz_yyy = pbuffer.data(idx_dip_hf + 496);

    auto tr_z_xxyyz_yyz = pbuffer.data(idx_dip_hf + 497);

    auto tr_z_xxyyz_yzz = pbuffer.data(idx_dip_hf + 498);

    auto tr_z_xxyyz_zzz = pbuffer.data(idx_dip_hf + 499);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_z_xxyy_xxy,  \
                             tr_z_xxyy_xyy,  \
                             tr_z_xxyyz_xxx, \
                             tr_z_xxyyz_xxy, \
                             tr_z_xxyyz_xxz, \
                             tr_z_xxyyz_xyy, \
                             tr_z_xxyyz_xyz, \
                             tr_z_xxyyz_xzz, \
                             tr_z_xxyyz_yyy, \
                             tr_z_xxyyz_yyz, \
                             tr_z_xxyyz_yzz, \
                             tr_z_xxyyz_zzz, \
                             tr_z_xxyz_xxx,  \
                             tr_z_xxyz_xxz,  \
                             tr_z_xxyz_xzz,  \
                             tr_z_xxz_xxx,   \
                             tr_z_xxz_xxz,   \
                             tr_z_xxz_xzz,   \
                             tr_z_xyyz_xyz,  \
                             tr_z_xyyz_yyy,  \
                             tr_z_xyyz_yyz,  \
                             tr_z_xyyz_yz,   \
                             tr_z_xyyz_yzz,  \
                             tr_z_xyyz_zzz,  \
                             tr_z_yyz_xyz,   \
                             tr_z_yyz_yyy,   \
                             tr_z_yyz_yyz,   \
                             tr_z_yyz_yzz,   \
                             tr_z_yyz_zzz,   \
                             ts_xxyy_xxy,    \
                             ts_xxyy_xyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyyz_xxx[i] = tr_z_xxz_xxx[i] * fe_0 + tr_z_xxyz_xxx[i] * pa_y[i];

        tr_z_xxyyz_xxy[i] = ts_xxyy_xxy[i] * fe_0 + tr_z_xxyy_xxy[i] * pa_z[i];

        tr_z_xxyyz_xxz[i] = tr_z_xxz_xxz[i] * fe_0 + tr_z_xxyz_xxz[i] * pa_y[i];

        tr_z_xxyyz_xyy[i] = ts_xxyy_xyy[i] * fe_0 + tr_z_xxyy_xyy[i] * pa_z[i];

        tr_z_xxyyz_xyz[i] = tr_z_yyz_xyz[i] * fe_0 + tr_z_xyyz_yz[i] * fe_0 + tr_z_xyyz_xyz[i] * pa_x[i];

        tr_z_xxyyz_xzz[i] = tr_z_xxz_xzz[i] * fe_0 + tr_z_xxyz_xzz[i] * pa_y[i];

        tr_z_xxyyz_yyy[i] = tr_z_yyz_yyy[i] * fe_0 + tr_z_xyyz_yyy[i] * pa_x[i];

        tr_z_xxyyz_yyz[i] = tr_z_yyz_yyz[i] * fe_0 + tr_z_xyyz_yyz[i] * pa_x[i];

        tr_z_xxyyz_yzz[i] = tr_z_yyz_yzz[i] * fe_0 + tr_z_xyyz_yzz[i] * pa_x[i];

        tr_z_xxyyz_zzz[i] = tr_z_yyz_zzz[i] * fe_0 + tr_z_xyyz_zzz[i] * pa_x[i];
    }

    // Set up 500-510 components of targeted buffer : HF

    auto tr_z_xxyzz_xxx = pbuffer.data(idx_dip_hf + 500);

    auto tr_z_xxyzz_xxy = pbuffer.data(idx_dip_hf + 501);

    auto tr_z_xxyzz_xxz = pbuffer.data(idx_dip_hf + 502);

    auto tr_z_xxyzz_xyy = pbuffer.data(idx_dip_hf + 503);

    auto tr_z_xxyzz_xyz = pbuffer.data(idx_dip_hf + 504);

    auto tr_z_xxyzz_xzz = pbuffer.data(idx_dip_hf + 505);

    auto tr_z_xxyzz_yyy = pbuffer.data(idx_dip_hf + 506);

    auto tr_z_xxyzz_yyz = pbuffer.data(idx_dip_hf + 507);

    auto tr_z_xxyzz_yzz = pbuffer.data(idx_dip_hf + 508);

    auto tr_z_xxyzz_zzz = pbuffer.data(idx_dip_hf + 509);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xxyzz_xxx, \
                             tr_z_xxyzz_xxy, \
                             tr_z_xxyzz_xxz, \
                             tr_z_xxyzz_xyy, \
                             tr_z_xxyzz_xyz, \
                             tr_z_xxyzz_xzz, \
                             tr_z_xxyzz_yyy, \
                             tr_z_xxyzz_yyz, \
                             tr_z_xxyzz_yzz, \
                             tr_z_xxyzz_zzz, \
                             tr_z_xxzz_xx,   \
                             tr_z_xxzz_xxx,  \
                             tr_z_xxzz_xxy,  \
                             tr_z_xxzz_xxz,  \
                             tr_z_xxzz_xy,   \
                             tr_z_xxzz_xyy,  \
                             tr_z_xxzz_xyz,  \
                             tr_z_xxzz_xz,   \
                             tr_z_xxzz_xzz,  \
                             tr_z_xxzz_zzz,  \
                             tr_z_xyzz_yyy,  \
                             tr_z_xyzz_yyz,  \
                             tr_z_xyzz_yzz,  \
                             tr_z_yzz_yyy,   \
                             tr_z_yzz_yyz,   \
                             tr_z_yzz_yzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxyzz_xxx[i] = tr_z_xxzz_xxx[i] * pa_y[i];

        tr_z_xxyzz_xxy[i] = tr_z_xxzz_xx[i] * fe_0 + tr_z_xxzz_xxy[i] * pa_y[i];

        tr_z_xxyzz_xxz[i] = tr_z_xxzz_xxz[i] * pa_y[i];

        tr_z_xxyzz_xyy[i] = 2.0 * tr_z_xxzz_xy[i] * fe_0 + tr_z_xxzz_xyy[i] * pa_y[i];

        tr_z_xxyzz_xyz[i] = tr_z_xxzz_xz[i] * fe_0 + tr_z_xxzz_xyz[i] * pa_y[i];

        tr_z_xxyzz_xzz[i] = tr_z_xxzz_xzz[i] * pa_y[i];

        tr_z_xxyzz_yyy[i] = tr_z_yzz_yyy[i] * fe_0 + tr_z_xyzz_yyy[i] * pa_x[i];

        tr_z_xxyzz_yyz[i] = tr_z_yzz_yyz[i] * fe_0 + tr_z_xyzz_yyz[i] * pa_x[i];

        tr_z_xxyzz_yzz[i] = tr_z_yzz_yzz[i] * fe_0 + tr_z_xyzz_yzz[i] * pa_x[i];

        tr_z_xxyzz_zzz[i] = tr_z_xxzz_zzz[i] * pa_y[i];
    }

    // Set up 510-520 components of targeted buffer : HF

    auto tr_z_xxzzz_xxx = pbuffer.data(idx_dip_hf + 510);

    auto tr_z_xxzzz_xxy = pbuffer.data(idx_dip_hf + 511);

    auto tr_z_xxzzz_xxz = pbuffer.data(idx_dip_hf + 512);

    auto tr_z_xxzzz_xyy = pbuffer.data(idx_dip_hf + 513);

    auto tr_z_xxzzz_xyz = pbuffer.data(idx_dip_hf + 514);

    auto tr_z_xxzzz_xzz = pbuffer.data(idx_dip_hf + 515);

    auto tr_z_xxzzz_yyy = pbuffer.data(idx_dip_hf + 516);

    auto tr_z_xxzzz_yyz = pbuffer.data(idx_dip_hf + 517);

    auto tr_z_xxzzz_yzz = pbuffer.data(idx_dip_hf + 518);

    auto tr_z_xxzzz_zzz = pbuffer.data(idx_dip_hf + 519);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xxzzz_xxx, \
                             tr_z_xxzzz_xxy, \
                             tr_z_xxzzz_xxz, \
                             tr_z_xxzzz_xyy, \
                             tr_z_xxzzz_xyz, \
                             tr_z_xxzzz_xzz, \
                             tr_z_xxzzz_yyy, \
                             tr_z_xxzzz_yyz, \
                             tr_z_xxzzz_yzz, \
                             tr_z_xxzzz_zzz, \
                             tr_z_xzzz_xx,   \
                             tr_z_xzzz_xxx,  \
                             tr_z_xzzz_xxy,  \
                             tr_z_xzzz_xxz,  \
                             tr_z_xzzz_xy,   \
                             tr_z_xzzz_xyy,  \
                             tr_z_xzzz_xyz,  \
                             tr_z_xzzz_xz,   \
                             tr_z_xzzz_xzz,  \
                             tr_z_xzzz_yy,   \
                             tr_z_xzzz_yyy,  \
                             tr_z_xzzz_yyz,  \
                             tr_z_xzzz_yz,   \
                             tr_z_xzzz_yzz,  \
                             tr_z_xzzz_zz,   \
                             tr_z_xzzz_zzz,  \
                             tr_z_zzz_xxx,   \
                             tr_z_zzz_xxy,   \
                             tr_z_zzz_xxz,   \
                             tr_z_zzz_xyy,   \
                             tr_z_zzz_xyz,   \
                             tr_z_zzz_xzz,   \
                             tr_z_zzz_yyy,   \
                             tr_z_zzz_yyz,   \
                             tr_z_zzz_yzz,   \
                             tr_z_zzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxzzz_xxx[i] = tr_z_zzz_xxx[i] * fe_0 + 3.0 * tr_z_xzzz_xx[i] * fe_0 + tr_z_xzzz_xxx[i] * pa_x[i];

        tr_z_xxzzz_xxy[i] = tr_z_zzz_xxy[i] * fe_0 + 2.0 * tr_z_xzzz_xy[i] * fe_0 + tr_z_xzzz_xxy[i] * pa_x[i];

        tr_z_xxzzz_xxz[i] = tr_z_zzz_xxz[i] * fe_0 + 2.0 * tr_z_xzzz_xz[i] * fe_0 + tr_z_xzzz_xxz[i] * pa_x[i];

        tr_z_xxzzz_xyy[i] = tr_z_zzz_xyy[i] * fe_0 + tr_z_xzzz_yy[i] * fe_0 + tr_z_xzzz_xyy[i] * pa_x[i];

        tr_z_xxzzz_xyz[i] = tr_z_zzz_xyz[i] * fe_0 + tr_z_xzzz_yz[i] * fe_0 + tr_z_xzzz_xyz[i] * pa_x[i];

        tr_z_xxzzz_xzz[i] = tr_z_zzz_xzz[i] * fe_0 + tr_z_xzzz_zz[i] * fe_0 + tr_z_xzzz_xzz[i] * pa_x[i];

        tr_z_xxzzz_yyy[i] = tr_z_zzz_yyy[i] * fe_0 + tr_z_xzzz_yyy[i] * pa_x[i];

        tr_z_xxzzz_yyz[i] = tr_z_zzz_yyz[i] * fe_0 + tr_z_xzzz_yyz[i] * pa_x[i];

        tr_z_xxzzz_yzz[i] = tr_z_zzz_yzz[i] * fe_0 + tr_z_xzzz_yzz[i] * pa_x[i];

        tr_z_xxzzz_zzz[i] = tr_z_zzz_zzz[i] * fe_0 + tr_z_xzzz_zzz[i] * pa_x[i];
    }

    // Set up 520-530 components of targeted buffer : HF

    auto tr_z_xyyyy_xxx = pbuffer.data(idx_dip_hf + 520);

    auto tr_z_xyyyy_xxy = pbuffer.data(idx_dip_hf + 521);

    auto tr_z_xyyyy_xxz = pbuffer.data(idx_dip_hf + 522);

    auto tr_z_xyyyy_xyy = pbuffer.data(idx_dip_hf + 523);

    auto tr_z_xyyyy_xyz = pbuffer.data(idx_dip_hf + 524);

    auto tr_z_xyyyy_xzz = pbuffer.data(idx_dip_hf + 525);

    auto tr_z_xyyyy_yyy = pbuffer.data(idx_dip_hf + 526);

    auto tr_z_xyyyy_yyz = pbuffer.data(idx_dip_hf + 527);

    auto tr_z_xyyyy_yzz = pbuffer.data(idx_dip_hf + 528);

    auto tr_z_xyyyy_zzz = pbuffer.data(idx_dip_hf + 529);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xyyyy_xxx, \
                             tr_z_xyyyy_xxy, \
                             tr_z_xyyyy_xxz, \
                             tr_z_xyyyy_xyy, \
                             tr_z_xyyyy_xyz, \
                             tr_z_xyyyy_xzz, \
                             tr_z_xyyyy_yyy, \
                             tr_z_xyyyy_yyz, \
                             tr_z_xyyyy_yzz, \
                             tr_z_xyyyy_zzz, \
                             tr_z_yyyy_xx,   \
                             tr_z_yyyy_xxx,  \
                             tr_z_yyyy_xxy,  \
                             tr_z_yyyy_xxz,  \
                             tr_z_yyyy_xy,   \
                             tr_z_yyyy_xyy,  \
                             tr_z_yyyy_xyz,  \
                             tr_z_yyyy_xz,   \
                             tr_z_yyyy_xzz,  \
                             tr_z_yyyy_yy,   \
                             tr_z_yyyy_yyy,  \
                             tr_z_yyyy_yyz,  \
                             tr_z_yyyy_yz,   \
                             tr_z_yyyy_yzz,  \
                             tr_z_yyyy_zz,   \
                             tr_z_yyyy_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyy_xxx[i] = 3.0 * tr_z_yyyy_xx[i] * fe_0 + tr_z_yyyy_xxx[i] * pa_x[i];

        tr_z_xyyyy_xxy[i] = 2.0 * tr_z_yyyy_xy[i] * fe_0 + tr_z_yyyy_xxy[i] * pa_x[i];

        tr_z_xyyyy_xxz[i] = 2.0 * tr_z_yyyy_xz[i] * fe_0 + tr_z_yyyy_xxz[i] * pa_x[i];

        tr_z_xyyyy_xyy[i] = tr_z_yyyy_yy[i] * fe_0 + tr_z_yyyy_xyy[i] * pa_x[i];

        tr_z_xyyyy_xyz[i] = tr_z_yyyy_yz[i] * fe_0 + tr_z_yyyy_xyz[i] * pa_x[i];

        tr_z_xyyyy_xzz[i] = tr_z_yyyy_zz[i] * fe_0 + tr_z_yyyy_xzz[i] * pa_x[i];

        tr_z_xyyyy_yyy[i] = tr_z_yyyy_yyy[i] * pa_x[i];

        tr_z_xyyyy_yyz[i] = tr_z_yyyy_yyz[i] * pa_x[i];

        tr_z_xyyyy_yzz[i] = tr_z_yyyy_yzz[i] * pa_x[i];

        tr_z_xyyyy_zzz[i] = tr_z_yyyy_zzz[i] * pa_x[i];
    }

    // Set up 530-540 components of targeted buffer : HF

    auto tr_z_xyyyz_xxx = pbuffer.data(idx_dip_hf + 530);

    auto tr_z_xyyyz_xxy = pbuffer.data(idx_dip_hf + 531);

    auto tr_z_xyyyz_xxz = pbuffer.data(idx_dip_hf + 532);

    auto tr_z_xyyyz_xyy = pbuffer.data(idx_dip_hf + 533);

    auto tr_z_xyyyz_xyz = pbuffer.data(idx_dip_hf + 534);

    auto tr_z_xyyyz_xzz = pbuffer.data(idx_dip_hf + 535);

    auto tr_z_xyyyz_yyy = pbuffer.data(idx_dip_hf + 536);

    auto tr_z_xyyyz_yyz = pbuffer.data(idx_dip_hf + 537);

    auto tr_z_xyyyz_yzz = pbuffer.data(idx_dip_hf + 538);

    auto tr_z_xyyyz_zzz = pbuffer.data(idx_dip_hf + 539);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xyyyz_xxx, \
                             tr_z_xyyyz_xxy, \
                             tr_z_xyyyz_xxz, \
                             tr_z_xyyyz_xyy, \
                             tr_z_xyyyz_xyz, \
                             tr_z_xyyyz_xzz, \
                             tr_z_xyyyz_yyy, \
                             tr_z_xyyyz_yyz, \
                             tr_z_xyyyz_yzz, \
                             tr_z_xyyyz_zzz, \
                             tr_z_yyyz_xx,   \
                             tr_z_yyyz_xxx,  \
                             tr_z_yyyz_xxy,  \
                             tr_z_yyyz_xxz,  \
                             tr_z_yyyz_xy,   \
                             tr_z_yyyz_xyy,  \
                             tr_z_yyyz_xyz,  \
                             tr_z_yyyz_xz,   \
                             tr_z_yyyz_xzz,  \
                             tr_z_yyyz_yy,   \
                             tr_z_yyyz_yyy,  \
                             tr_z_yyyz_yyz,  \
                             tr_z_yyyz_yz,   \
                             tr_z_yyyz_yzz,  \
                             tr_z_yyyz_zz,   \
                             tr_z_yyyz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyyz_xxx[i] = 3.0 * tr_z_yyyz_xx[i] * fe_0 + tr_z_yyyz_xxx[i] * pa_x[i];

        tr_z_xyyyz_xxy[i] = 2.0 * tr_z_yyyz_xy[i] * fe_0 + tr_z_yyyz_xxy[i] * pa_x[i];

        tr_z_xyyyz_xxz[i] = 2.0 * tr_z_yyyz_xz[i] * fe_0 + tr_z_yyyz_xxz[i] * pa_x[i];

        tr_z_xyyyz_xyy[i] = tr_z_yyyz_yy[i] * fe_0 + tr_z_yyyz_xyy[i] * pa_x[i];

        tr_z_xyyyz_xyz[i] = tr_z_yyyz_yz[i] * fe_0 + tr_z_yyyz_xyz[i] * pa_x[i];

        tr_z_xyyyz_xzz[i] = tr_z_yyyz_zz[i] * fe_0 + tr_z_yyyz_xzz[i] * pa_x[i];

        tr_z_xyyyz_yyy[i] = tr_z_yyyz_yyy[i] * pa_x[i];

        tr_z_xyyyz_yyz[i] = tr_z_yyyz_yyz[i] * pa_x[i];

        tr_z_xyyyz_yzz[i] = tr_z_yyyz_yzz[i] * pa_x[i];

        tr_z_xyyyz_zzz[i] = tr_z_yyyz_zzz[i] * pa_x[i];
    }

    // Set up 540-550 components of targeted buffer : HF

    auto tr_z_xyyzz_xxx = pbuffer.data(idx_dip_hf + 540);

    auto tr_z_xyyzz_xxy = pbuffer.data(idx_dip_hf + 541);

    auto tr_z_xyyzz_xxz = pbuffer.data(idx_dip_hf + 542);

    auto tr_z_xyyzz_xyy = pbuffer.data(idx_dip_hf + 543);

    auto tr_z_xyyzz_xyz = pbuffer.data(idx_dip_hf + 544);

    auto tr_z_xyyzz_xzz = pbuffer.data(idx_dip_hf + 545);

    auto tr_z_xyyzz_yyy = pbuffer.data(idx_dip_hf + 546);

    auto tr_z_xyyzz_yyz = pbuffer.data(idx_dip_hf + 547);

    auto tr_z_xyyzz_yzz = pbuffer.data(idx_dip_hf + 548);

    auto tr_z_xyyzz_zzz = pbuffer.data(idx_dip_hf + 549);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xyyzz_xxx, \
                             tr_z_xyyzz_xxy, \
                             tr_z_xyyzz_xxz, \
                             tr_z_xyyzz_xyy, \
                             tr_z_xyyzz_xyz, \
                             tr_z_xyyzz_xzz, \
                             tr_z_xyyzz_yyy, \
                             tr_z_xyyzz_yyz, \
                             tr_z_xyyzz_yzz, \
                             tr_z_xyyzz_zzz, \
                             tr_z_yyzz_xx,   \
                             tr_z_yyzz_xxx,  \
                             tr_z_yyzz_xxy,  \
                             tr_z_yyzz_xxz,  \
                             tr_z_yyzz_xy,   \
                             tr_z_yyzz_xyy,  \
                             tr_z_yyzz_xyz,  \
                             tr_z_yyzz_xz,   \
                             tr_z_yyzz_xzz,  \
                             tr_z_yyzz_yy,   \
                             tr_z_yyzz_yyy,  \
                             tr_z_yyzz_yyz,  \
                             tr_z_yyzz_yz,   \
                             tr_z_yyzz_yzz,  \
                             tr_z_yyzz_zz,   \
                             tr_z_yyzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyyzz_xxx[i] = 3.0 * tr_z_yyzz_xx[i] * fe_0 + tr_z_yyzz_xxx[i] * pa_x[i];

        tr_z_xyyzz_xxy[i] = 2.0 * tr_z_yyzz_xy[i] * fe_0 + tr_z_yyzz_xxy[i] * pa_x[i];

        tr_z_xyyzz_xxz[i] = 2.0 * tr_z_yyzz_xz[i] * fe_0 + tr_z_yyzz_xxz[i] * pa_x[i];

        tr_z_xyyzz_xyy[i] = tr_z_yyzz_yy[i] * fe_0 + tr_z_yyzz_xyy[i] * pa_x[i];

        tr_z_xyyzz_xyz[i] = tr_z_yyzz_yz[i] * fe_0 + tr_z_yyzz_xyz[i] * pa_x[i];

        tr_z_xyyzz_xzz[i] = tr_z_yyzz_zz[i] * fe_0 + tr_z_yyzz_xzz[i] * pa_x[i];

        tr_z_xyyzz_yyy[i] = tr_z_yyzz_yyy[i] * pa_x[i];

        tr_z_xyyzz_yyz[i] = tr_z_yyzz_yyz[i] * pa_x[i];

        tr_z_xyyzz_yzz[i] = tr_z_yyzz_yzz[i] * pa_x[i];

        tr_z_xyyzz_zzz[i] = tr_z_yyzz_zzz[i] * pa_x[i];
    }

    // Set up 550-560 components of targeted buffer : HF

    auto tr_z_xyzzz_xxx = pbuffer.data(idx_dip_hf + 550);

    auto tr_z_xyzzz_xxy = pbuffer.data(idx_dip_hf + 551);

    auto tr_z_xyzzz_xxz = pbuffer.data(idx_dip_hf + 552);

    auto tr_z_xyzzz_xyy = pbuffer.data(idx_dip_hf + 553);

    auto tr_z_xyzzz_xyz = pbuffer.data(idx_dip_hf + 554);

    auto tr_z_xyzzz_xzz = pbuffer.data(idx_dip_hf + 555);

    auto tr_z_xyzzz_yyy = pbuffer.data(idx_dip_hf + 556);

    auto tr_z_xyzzz_yyz = pbuffer.data(idx_dip_hf + 557);

    auto tr_z_xyzzz_yzz = pbuffer.data(idx_dip_hf + 558);

    auto tr_z_xyzzz_zzz = pbuffer.data(idx_dip_hf + 559);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xyzzz_xxx, \
                             tr_z_xyzzz_xxy, \
                             tr_z_xyzzz_xxz, \
                             tr_z_xyzzz_xyy, \
                             tr_z_xyzzz_xyz, \
                             tr_z_xyzzz_xzz, \
                             tr_z_xyzzz_yyy, \
                             tr_z_xyzzz_yyz, \
                             tr_z_xyzzz_yzz, \
                             tr_z_xyzzz_zzz, \
                             tr_z_xzzz_xxx,  \
                             tr_z_xzzz_xxz,  \
                             tr_z_xzzz_xzz,  \
                             tr_z_yzzz_xxy,  \
                             tr_z_yzzz_xy,   \
                             tr_z_yzzz_xyy,  \
                             tr_z_yzzz_xyz,  \
                             tr_z_yzzz_yy,   \
                             tr_z_yzzz_yyy,  \
                             tr_z_yzzz_yyz,  \
                             tr_z_yzzz_yz,   \
                             tr_z_yzzz_yzz,  \
                             tr_z_yzzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyzzz_xxx[i] = tr_z_xzzz_xxx[i] * pa_y[i];

        tr_z_xyzzz_xxy[i] = 2.0 * tr_z_yzzz_xy[i] * fe_0 + tr_z_yzzz_xxy[i] * pa_x[i];

        tr_z_xyzzz_xxz[i] = tr_z_xzzz_xxz[i] * pa_y[i];

        tr_z_xyzzz_xyy[i] = tr_z_yzzz_yy[i] * fe_0 + tr_z_yzzz_xyy[i] * pa_x[i];

        tr_z_xyzzz_xyz[i] = tr_z_yzzz_yz[i] * fe_0 + tr_z_yzzz_xyz[i] * pa_x[i];

        tr_z_xyzzz_xzz[i] = tr_z_xzzz_xzz[i] * pa_y[i];

        tr_z_xyzzz_yyy[i] = tr_z_yzzz_yyy[i] * pa_x[i];

        tr_z_xyzzz_yyz[i] = tr_z_yzzz_yyz[i] * pa_x[i];

        tr_z_xyzzz_yzz[i] = tr_z_yzzz_yzz[i] * pa_x[i];

        tr_z_xyzzz_zzz[i] = tr_z_yzzz_zzz[i] * pa_x[i];
    }

    // Set up 560-570 components of targeted buffer : HF

    auto tr_z_xzzzz_xxx = pbuffer.data(idx_dip_hf + 560);

    auto tr_z_xzzzz_xxy = pbuffer.data(idx_dip_hf + 561);

    auto tr_z_xzzzz_xxz = pbuffer.data(idx_dip_hf + 562);

    auto tr_z_xzzzz_xyy = pbuffer.data(idx_dip_hf + 563);

    auto tr_z_xzzzz_xyz = pbuffer.data(idx_dip_hf + 564);

    auto tr_z_xzzzz_xzz = pbuffer.data(idx_dip_hf + 565);

    auto tr_z_xzzzz_yyy = pbuffer.data(idx_dip_hf + 566);

    auto tr_z_xzzzz_yyz = pbuffer.data(idx_dip_hf + 567);

    auto tr_z_xzzzz_yzz = pbuffer.data(idx_dip_hf + 568);

    auto tr_z_xzzzz_zzz = pbuffer.data(idx_dip_hf + 569);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xzzzz_xxx, \
                             tr_z_xzzzz_xxy, \
                             tr_z_xzzzz_xxz, \
                             tr_z_xzzzz_xyy, \
                             tr_z_xzzzz_xyz, \
                             tr_z_xzzzz_xzz, \
                             tr_z_xzzzz_yyy, \
                             tr_z_xzzzz_yyz, \
                             tr_z_xzzzz_yzz, \
                             tr_z_xzzzz_zzz, \
                             tr_z_zzzz_xx,   \
                             tr_z_zzzz_xxx,  \
                             tr_z_zzzz_xxy,  \
                             tr_z_zzzz_xxz,  \
                             tr_z_zzzz_xy,   \
                             tr_z_zzzz_xyy,  \
                             tr_z_zzzz_xyz,  \
                             tr_z_zzzz_xz,   \
                             tr_z_zzzz_xzz,  \
                             tr_z_zzzz_yy,   \
                             tr_z_zzzz_yyy,  \
                             tr_z_zzzz_yyz,  \
                             tr_z_zzzz_yz,   \
                             tr_z_zzzz_yzz,  \
                             tr_z_zzzz_zz,   \
                             tr_z_zzzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzzzz_xxx[i] = 3.0 * tr_z_zzzz_xx[i] * fe_0 + tr_z_zzzz_xxx[i] * pa_x[i];

        tr_z_xzzzz_xxy[i] = 2.0 * tr_z_zzzz_xy[i] * fe_0 + tr_z_zzzz_xxy[i] * pa_x[i];

        tr_z_xzzzz_xxz[i] = 2.0 * tr_z_zzzz_xz[i] * fe_0 + tr_z_zzzz_xxz[i] * pa_x[i];

        tr_z_xzzzz_xyy[i] = tr_z_zzzz_yy[i] * fe_0 + tr_z_zzzz_xyy[i] * pa_x[i];

        tr_z_xzzzz_xyz[i] = tr_z_zzzz_yz[i] * fe_0 + tr_z_zzzz_xyz[i] * pa_x[i];

        tr_z_xzzzz_xzz[i] = tr_z_zzzz_zz[i] * fe_0 + tr_z_zzzz_xzz[i] * pa_x[i];

        tr_z_xzzzz_yyy[i] = tr_z_zzzz_yyy[i] * pa_x[i];

        tr_z_xzzzz_yyz[i] = tr_z_zzzz_yyz[i] * pa_x[i];

        tr_z_xzzzz_yzz[i] = tr_z_zzzz_yzz[i] * pa_x[i];

        tr_z_xzzzz_zzz[i] = tr_z_zzzz_zzz[i] * pa_x[i];
    }

    // Set up 570-580 components of targeted buffer : HF

    auto tr_z_yyyyy_xxx = pbuffer.data(idx_dip_hf + 570);

    auto tr_z_yyyyy_xxy = pbuffer.data(idx_dip_hf + 571);

    auto tr_z_yyyyy_xxz = pbuffer.data(idx_dip_hf + 572);

    auto tr_z_yyyyy_xyy = pbuffer.data(idx_dip_hf + 573);

    auto tr_z_yyyyy_xyz = pbuffer.data(idx_dip_hf + 574);

    auto tr_z_yyyyy_xzz = pbuffer.data(idx_dip_hf + 575);

    auto tr_z_yyyyy_yyy = pbuffer.data(idx_dip_hf + 576);

    auto tr_z_yyyyy_yyz = pbuffer.data(idx_dip_hf + 577);

    auto tr_z_yyyyy_yzz = pbuffer.data(idx_dip_hf + 578);

    auto tr_z_yyyyy_zzz = pbuffer.data(idx_dip_hf + 579);

#pragma omp simd aligned(pa_y,               \
                             tr_z_yyy_xxx,   \
                             tr_z_yyy_xxy,   \
                             tr_z_yyy_xxz,   \
                             tr_z_yyy_xyy,   \
                             tr_z_yyy_xyz,   \
                             tr_z_yyy_xzz,   \
                             tr_z_yyy_yyy,   \
                             tr_z_yyy_yyz,   \
                             tr_z_yyy_yzz,   \
                             tr_z_yyy_zzz,   \
                             tr_z_yyyy_xx,   \
                             tr_z_yyyy_xxx,  \
                             tr_z_yyyy_xxy,  \
                             tr_z_yyyy_xxz,  \
                             tr_z_yyyy_xy,   \
                             tr_z_yyyy_xyy,  \
                             tr_z_yyyy_xyz,  \
                             tr_z_yyyy_xz,   \
                             tr_z_yyyy_xzz,  \
                             tr_z_yyyy_yy,   \
                             tr_z_yyyy_yyy,  \
                             tr_z_yyyy_yyz,  \
                             tr_z_yyyy_yz,   \
                             tr_z_yyyy_yzz,  \
                             tr_z_yyyy_zz,   \
                             tr_z_yyyy_zzz,  \
                             tr_z_yyyyy_xxx, \
                             tr_z_yyyyy_xxy, \
                             tr_z_yyyyy_xxz, \
                             tr_z_yyyyy_xyy, \
                             tr_z_yyyyy_xyz, \
                             tr_z_yyyyy_xzz, \
                             tr_z_yyyyy_yyy, \
                             tr_z_yyyyy_yyz, \
                             tr_z_yyyyy_yzz, \
                             tr_z_yyyyy_zzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyy_xxx[i] = 4.0 * tr_z_yyy_xxx[i] * fe_0 + tr_z_yyyy_xxx[i] * pa_y[i];

        tr_z_yyyyy_xxy[i] = 4.0 * tr_z_yyy_xxy[i] * fe_0 + tr_z_yyyy_xx[i] * fe_0 + tr_z_yyyy_xxy[i] * pa_y[i];

        tr_z_yyyyy_xxz[i] = 4.0 * tr_z_yyy_xxz[i] * fe_0 + tr_z_yyyy_xxz[i] * pa_y[i];

        tr_z_yyyyy_xyy[i] = 4.0 * tr_z_yyy_xyy[i] * fe_0 + 2.0 * tr_z_yyyy_xy[i] * fe_0 + tr_z_yyyy_xyy[i] * pa_y[i];

        tr_z_yyyyy_xyz[i] = 4.0 * tr_z_yyy_xyz[i] * fe_0 + tr_z_yyyy_xz[i] * fe_0 + tr_z_yyyy_xyz[i] * pa_y[i];

        tr_z_yyyyy_xzz[i] = 4.0 * tr_z_yyy_xzz[i] * fe_0 + tr_z_yyyy_xzz[i] * pa_y[i];

        tr_z_yyyyy_yyy[i] = 4.0 * tr_z_yyy_yyy[i] * fe_0 + 3.0 * tr_z_yyyy_yy[i] * fe_0 + tr_z_yyyy_yyy[i] * pa_y[i];

        tr_z_yyyyy_yyz[i] = 4.0 * tr_z_yyy_yyz[i] * fe_0 + 2.0 * tr_z_yyyy_yz[i] * fe_0 + tr_z_yyyy_yyz[i] * pa_y[i];

        tr_z_yyyyy_yzz[i] = 4.0 * tr_z_yyy_yzz[i] * fe_0 + tr_z_yyyy_zz[i] * fe_0 + tr_z_yyyy_yzz[i] * pa_y[i];

        tr_z_yyyyy_zzz[i] = 4.0 * tr_z_yyy_zzz[i] * fe_0 + tr_z_yyyy_zzz[i] * pa_y[i];
    }

    // Set up 580-590 components of targeted buffer : HF

    auto tr_z_yyyyz_xxx = pbuffer.data(idx_dip_hf + 580);

    auto tr_z_yyyyz_xxy = pbuffer.data(idx_dip_hf + 581);

    auto tr_z_yyyyz_xxz = pbuffer.data(idx_dip_hf + 582);

    auto tr_z_yyyyz_xyy = pbuffer.data(idx_dip_hf + 583);

    auto tr_z_yyyyz_xyz = pbuffer.data(idx_dip_hf + 584);

    auto tr_z_yyyyz_xzz = pbuffer.data(idx_dip_hf + 585);

    auto tr_z_yyyyz_yyy = pbuffer.data(idx_dip_hf + 586);

    auto tr_z_yyyyz_yyz = pbuffer.data(idx_dip_hf + 587);

    auto tr_z_yyyyz_yzz = pbuffer.data(idx_dip_hf + 588);

    auto tr_z_yyyyz_zzz = pbuffer.data(idx_dip_hf + 589);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_z_yyyy_xxy,  \
                             tr_z_yyyy_xyy,  \
                             tr_z_yyyy_yyy,  \
                             tr_z_yyyyz_xxx, \
                             tr_z_yyyyz_xxy, \
                             tr_z_yyyyz_xxz, \
                             tr_z_yyyyz_xyy, \
                             tr_z_yyyyz_xyz, \
                             tr_z_yyyyz_xzz, \
                             tr_z_yyyyz_yyy, \
                             tr_z_yyyyz_yyz, \
                             tr_z_yyyyz_yzz, \
                             tr_z_yyyyz_zzz, \
                             tr_z_yyyz_xxx,  \
                             tr_z_yyyz_xxz,  \
                             tr_z_yyyz_xyz,  \
                             tr_z_yyyz_xz,   \
                             tr_z_yyyz_xzz,  \
                             tr_z_yyyz_yyz,  \
                             tr_z_yyyz_yz,   \
                             tr_z_yyyz_yzz,  \
                             tr_z_yyyz_zz,   \
                             tr_z_yyyz_zzz,  \
                             tr_z_yyz_xxx,   \
                             tr_z_yyz_xxz,   \
                             tr_z_yyz_xyz,   \
                             tr_z_yyz_xzz,   \
                             tr_z_yyz_yyz,   \
                             tr_z_yyz_yzz,   \
                             tr_z_yyz_zzz,   \
                             ts_yyyy_xxy,    \
                             ts_yyyy_xyy,    \
                             ts_yyyy_yyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyyz_xxx[i] = 3.0 * tr_z_yyz_xxx[i] * fe_0 + tr_z_yyyz_xxx[i] * pa_y[i];

        tr_z_yyyyz_xxy[i] = ts_yyyy_xxy[i] * fe_0 + tr_z_yyyy_xxy[i] * pa_z[i];

        tr_z_yyyyz_xxz[i] = 3.0 * tr_z_yyz_xxz[i] * fe_0 + tr_z_yyyz_xxz[i] * pa_y[i];

        tr_z_yyyyz_xyy[i] = ts_yyyy_xyy[i] * fe_0 + tr_z_yyyy_xyy[i] * pa_z[i];

        tr_z_yyyyz_xyz[i] = 3.0 * tr_z_yyz_xyz[i] * fe_0 + tr_z_yyyz_xz[i] * fe_0 + tr_z_yyyz_xyz[i] * pa_y[i];

        tr_z_yyyyz_xzz[i] = 3.0 * tr_z_yyz_xzz[i] * fe_0 + tr_z_yyyz_xzz[i] * pa_y[i];

        tr_z_yyyyz_yyy[i] = ts_yyyy_yyy[i] * fe_0 + tr_z_yyyy_yyy[i] * pa_z[i];

        tr_z_yyyyz_yyz[i] = 3.0 * tr_z_yyz_yyz[i] * fe_0 + 2.0 * tr_z_yyyz_yz[i] * fe_0 + tr_z_yyyz_yyz[i] * pa_y[i];

        tr_z_yyyyz_yzz[i] = 3.0 * tr_z_yyz_yzz[i] * fe_0 + tr_z_yyyz_zz[i] * fe_0 + tr_z_yyyz_yzz[i] * pa_y[i];

        tr_z_yyyyz_zzz[i] = 3.0 * tr_z_yyz_zzz[i] * fe_0 + tr_z_yyyz_zzz[i] * pa_y[i];
    }

    // Set up 590-600 components of targeted buffer : HF

    auto tr_z_yyyzz_xxx = pbuffer.data(idx_dip_hf + 590);

    auto tr_z_yyyzz_xxy = pbuffer.data(idx_dip_hf + 591);

    auto tr_z_yyyzz_xxz = pbuffer.data(idx_dip_hf + 592);

    auto tr_z_yyyzz_xyy = pbuffer.data(idx_dip_hf + 593);

    auto tr_z_yyyzz_xyz = pbuffer.data(idx_dip_hf + 594);

    auto tr_z_yyyzz_xzz = pbuffer.data(idx_dip_hf + 595);

    auto tr_z_yyyzz_yyy = pbuffer.data(idx_dip_hf + 596);

    auto tr_z_yyyzz_yyz = pbuffer.data(idx_dip_hf + 597);

    auto tr_z_yyyzz_yzz = pbuffer.data(idx_dip_hf + 598);

    auto tr_z_yyyzz_zzz = pbuffer.data(idx_dip_hf + 599);

#pragma omp simd aligned(pa_y,               \
                             tr_z_yyyzz_xxx, \
                             tr_z_yyyzz_xxy, \
                             tr_z_yyyzz_xxz, \
                             tr_z_yyyzz_xyy, \
                             tr_z_yyyzz_xyz, \
                             tr_z_yyyzz_xzz, \
                             tr_z_yyyzz_yyy, \
                             tr_z_yyyzz_yyz, \
                             tr_z_yyyzz_yzz, \
                             tr_z_yyyzz_zzz, \
                             tr_z_yyzz_xx,   \
                             tr_z_yyzz_xxx,  \
                             tr_z_yyzz_xxy,  \
                             tr_z_yyzz_xxz,  \
                             tr_z_yyzz_xy,   \
                             tr_z_yyzz_xyy,  \
                             tr_z_yyzz_xyz,  \
                             tr_z_yyzz_xz,   \
                             tr_z_yyzz_xzz,  \
                             tr_z_yyzz_yy,   \
                             tr_z_yyzz_yyy,  \
                             tr_z_yyzz_yyz,  \
                             tr_z_yyzz_yz,   \
                             tr_z_yyzz_yzz,  \
                             tr_z_yyzz_zz,   \
                             tr_z_yyzz_zzz,  \
                             tr_z_yzz_xxx,   \
                             tr_z_yzz_xxy,   \
                             tr_z_yzz_xxz,   \
                             tr_z_yzz_xyy,   \
                             tr_z_yzz_xyz,   \
                             tr_z_yzz_xzz,   \
                             tr_z_yzz_yyy,   \
                             tr_z_yzz_yyz,   \
                             tr_z_yzz_yzz,   \
                             tr_z_yzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyyzz_xxx[i] = 2.0 * tr_z_yzz_xxx[i] * fe_0 + tr_z_yyzz_xxx[i] * pa_y[i];

        tr_z_yyyzz_xxy[i] = 2.0 * tr_z_yzz_xxy[i] * fe_0 + tr_z_yyzz_xx[i] * fe_0 + tr_z_yyzz_xxy[i] * pa_y[i];

        tr_z_yyyzz_xxz[i] = 2.0 * tr_z_yzz_xxz[i] * fe_0 + tr_z_yyzz_xxz[i] * pa_y[i];

        tr_z_yyyzz_xyy[i] = 2.0 * tr_z_yzz_xyy[i] * fe_0 + 2.0 * tr_z_yyzz_xy[i] * fe_0 + tr_z_yyzz_xyy[i] * pa_y[i];

        tr_z_yyyzz_xyz[i] = 2.0 * tr_z_yzz_xyz[i] * fe_0 + tr_z_yyzz_xz[i] * fe_0 + tr_z_yyzz_xyz[i] * pa_y[i];

        tr_z_yyyzz_xzz[i] = 2.0 * tr_z_yzz_xzz[i] * fe_0 + tr_z_yyzz_xzz[i] * pa_y[i];

        tr_z_yyyzz_yyy[i] = 2.0 * tr_z_yzz_yyy[i] * fe_0 + 3.0 * tr_z_yyzz_yy[i] * fe_0 + tr_z_yyzz_yyy[i] * pa_y[i];

        tr_z_yyyzz_yyz[i] = 2.0 * tr_z_yzz_yyz[i] * fe_0 + 2.0 * tr_z_yyzz_yz[i] * fe_0 + tr_z_yyzz_yyz[i] * pa_y[i];

        tr_z_yyyzz_yzz[i] = 2.0 * tr_z_yzz_yzz[i] * fe_0 + tr_z_yyzz_zz[i] * fe_0 + tr_z_yyzz_yzz[i] * pa_y[i];

        tr_z_yyyzz_zzz[i] = 2.0 * tr_z_yzz_zzz[i] * fe_0 + tr_z_yyzz_zzz[i] * pa_y[i];
    }

    // Set up 600-610 components of targeted buffer : HF

    auto tr_z_yyzzz_xxx = pbuffer.data(idx_dip_hf + 600);

    auto tr_z_yyzzz_xxy = pbuffer.data(idx_dip_hf + 601);

    auto tr_z_yyzzz_xxz = pbuffer.data(idx_dip_hf + 602);

    auto tr_z_yyzzz_xyy = pbuffer.data(idx_dip_hf + 603);

    auto tr_z_yyzzz_xyz = pbuffer.data(idx_dip_hf + 604);

    auto tr_z_yyzzz_xzz = pbuffer.data(idx_dip_hf + 605);

    auto tr_z_yyzzz_yyy = pbuffer.data(idx_dip_hf + 606);

    auto tr_z_yyzzz_yyz = pbuffer.data(idx_dip_hf + 607);

    auto tr_z_yyzzz_yzz = pbuffer.data(idx_dip_hf + 608);

    auto tr_z_yyzzz_zzz = pbuffer.data(idx_dip_hf + 609);

#pragma omp simd aligned(pa_y,               \
                             tr_z_yyzzz_xxx, \
                             tr_z_yyzzz_xxy, \
                             tr_z_yyzzz_xxz, \
                             tr_z_yyzzz_xyy, \
                             tr_z_yyzzz_xyz, \
                             tr_z_yyzzz_xzz, \
                             tr_z_yyzzz_yyy, \
                             tr_z_yyzzz_yyz, \
                             tr_z_yyzzz_yzz, \
                             tr_z_yyzzz_zzz, \
                             tr_z_yzzz_xx,   \
                             tr_z_yzzz_xxx,  \
                             tr_z_yzzz_xxy,  \
                             tr_z_yzzz_xxz,  \
                             tr_z_yzzz_xy,   \
                             tr_z_yzzz_xyy,  \
                             tr_z_yzzz_xyz,  \
                             tr_z_yzzz_xz,   \
                             tr_z_yzzz_xzz,  \
                             tr_z_yzzz_yy,   \
                             tr_z_yzzz_yyy,  \
                             tr_z_yzzz_yyz,  \
                             tr_z_yzzz_yz,   \
                             tr_z_yzzz_yzz,  \
                             tr_z_yzzz_zz,   \
                             tr_z_yzzz_zzz,  \
                             tr_z_zzz_xxx,   \
                             tr_z_zzz_xxy,   \
                             tr_z_zzz_xxz,   \
                             tr_z_zzz_xyy,   \
                             tr_z_zzz_xyz,   \
                             tr_z_zzz_xzz,   \
                             tr_z_zzz_yyy,   \
                             tr_z_zzz_yyz,   \
                             tr_z_zzz_yzz,   \
                             tr_z_zzz_zzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyzzz_xxx[i] = tr_z_zzz_xxx[i] * fe_0 + tr_z_yzzz_xxx[i] * pa_y[i];

        tr_z_yyzzz_xxy[i] = tr_z_zzz_xxy[i] * fe_0 + tr_z_yzzz_xx[i] * fe_0 + tr_z_yzzz_xxy[i] * pa_y[i];

        tr_z_yyzzz_xxz[i] = tr_z_zzz_xxz[i] * fe_0 + tr_z_yzzz_xxz[i] * pa_y[i];

        tr_z_yyzzz_xyy[i] = tr_z_zzz_xyy[i] * fe_0 + 2.0 * tr_z_yzzz_xy[i] * fe_0 + tr_z_yzzz_xyy[i] * pa_y[i];

        tr_z_yyzzz_xyz[i] = tr_z_zzz_xyz[i] * fe_0 + tr_z_yzzz_xz[i] * fe_0 + tr_z_yzzz_xyz[i] * pa_y[i];

        tr_z_yyzzz_xzz[i] = tr_z_zzz_xzz[i] * fe_0 + tr_z_yzzz_xzz[i] * pa_y[i];

        tr_z_yyzzz_yyy[i] = tr_z_zzz_yyy[i] * fe_0 + 3.0 * tr_z_yzzz_yy[i] * fe_0 + tr_z_yzzz_yyy[i] * pa_y[i];

        tr_z_yyzzz_yyz[i] = tr_z_zzz_yyz[i] * fe_0 + 2.0 * tr_z_yzzz_yz[i] * fe_0 + tr_z_yzzz_yyz[i] * pa_y[i];

        tr_z_yyzzz_yzz[i] = tr_z_zzz_yzz[i] * fe_0 + tr_z_yzzz_zz[i] * fe_0 + tr_z_yzzz_yzz[i] * pa_y[i];

        tr_z_yyzzz_zzz[i] = tr_z_zzz_zzz[i] * fe_0 + tr_z_yzzz_zzz[i] * pa_y[i];
    }

    // Set up 610-620 components of targeted buffer : HF

    auto tr_z_yzzzz_xxx = pbuffer.data(idx_dip_hf + 610);

    auto tr_z_yzzzz_xxy = pbuffer.data(idx_dip_hf + 611);

    auto tr_z_yzzzz_xxz = pbuffer.data(idx_dip_hf + 612);

    auto tr_z_yzzzz_xyy = pbuffer.data(idx_dip_hf + 613);

    auto tr_z_yzzzz_xyz = pbuffer.data(idx_dip_hf + 614);

    auto tr_z_yzzzz_xzz = pbuffer.data(idx_dip_hf + 615);

    auto tr_z_yzzzz_yyy = pbuffer.data(idx_dip_hf + 616);

    auto tr_z_yzzzz_yyz = pbuffer.data(idx_dip_hf + 617);

    auto tr_z_yzzzz_yzz = pbuffer.data(idx_dip_hf + 618);

    auto tr_z_yzzzz_zzz = pbuffer.data(idx_dip_hf + 619);

#pragma omp simd aligned(pa_y,               \
                             tr_z_yzzzz_xxx, \
                             tr_z_yzzzz_xxy, \
                             tr_z_yzzzz_xxz, \
                             tr_z_yzzzz_xyy, \
                             tr_z_yzzzz_xyz, \
                             tr_z_yzzzz_xzz, \
                             tr_z_yzzzz_yyy, \
                             tr_z_yzzzz_yyz, \
                             tr_z_yzzzz_yzz, \
                             tr_z_yzzzz_zzz, \
                             tr_z_zzzz_xx,   \
                             tr_z_zzzz_xxx,  \
                             tr_z_zzzz_xxy,  \
                             tr_z_zzzz_xxz,  \
                             tr_z_zzzz_xy,   \
                             tr_z_zzzz_xyy,  \
                             tr_z_zzzz_xyz,  \
                             tr_z_zzzz_xz,   \
                             tr_z_zzzz_xzz,  \
                             tr_z_zzzz_yy,   \
                             tr_z_zzzz_yyy,  \
                             tr_z_zzzz_yyz,  \
                             tr_z_zzzz_yz,   \
                             tr_z_zzzz_yzz,  \
                             tr_z_zzzz_zz,   \
                             tr_z_zzzz_zzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzzzz_xxx[i] = tr_z_zzzz_xxx[i] * pa_y[i];

        tr_z_yzzzz_xxy[i] = tr_z_zzzz_xx[i] * fe_0 + tr_z_zzzz_xxy[i] * pa_y[i];

        tr_z_yzzzz_xxz[i] = tr_z_zzzz_xxz[i] * pa_y[i];

        tr_z_yzzzz_xyy[i] = 2.0 * tr_z_zzzz_xy[i] * fe_0 + tr_z_zzzz_xyy[i] * pa_y[i];

        tr_z_yzzzz_xyz[i] = tr_z_zzzz_xz[i] * fe_0 + tr_z_zzzz_xyz[i] * pa_y[i];

        tr_z_yzzzz_xzz[i] = tr_z_zzzz_xzz[i] * pa_y[i];

        tr_z_yzzzz_yyy[i] = 3.0 * tr_z_zzzz_yy[i] * fe_0 + tr_z_zzzz_yyy[i] * pa_y[i];

        tr_z_yzzzz_yyz[i] = 2.0 * tr_z_zzzz_yz[i] * fe_0 + tr_z_zzzz_yyz[i] * pa_y[i];

        tr_z_yzzzz_yzz[i] = tr_z_zzzz_zz[i] * fe_0 + tr_z_zzzz_yzz[i] * pa_y[i];

        tr_z_yzzzz_zzz[i] = tr_z_zzzz_zzz[i] * pa_y[i];
    }

    // Set up 620-630 components of targeted buffer : HF

    auto tr_z_zzzzz_xxx = pbuffer.data(idx_dip_hf + 620);

    auto tr_z_zzzzz_xxy = pbuffer.data(idx_dip_hf + 621);

    auto tr_z_zzzzz_xxz = pbuffer.data(idx_dip_hf + 622);

    auto tr_z_zzzzz_xyy = pbuffer.data(idx_dip_hf + 623);

    auto tr_z_zzzzz_xyz = pbuffer.data(idx_dip_hf + 624);

    auto tr_z_zzzzz_xzz = pbuffer.data(idx_dip_hf + 625);

    auto tr_z_zzzzz_yyy = pbuffer.data(idx_dip_hf + 626);

    auto tr_z_zzzzz_yyz = pbuffer.data(idx_dip_hf + 627);

    auto tr_z_zzzzz_yzz = pbuffer.data(idx_dip_hf + 628);

    auto tr_z_zzzzz_zzz = pbuffer.data(idx_dip_hf + 629);

#pragma omp simd aligned(pa_z,               \
                             tr_z_zzz_xxx,   \
                             tr_z_zzz_xxy,   \
                             tr_z_zzz_xxz,   \
                             tr_z_zzz_xyy,   \
                             tr_z_zzz_xyz,   \
                             tr_z_zzz_xzz,   \
                             tr_z_zzz_yyy,   \
                             tr_z_zzz_yyz,   \
                             tr_z_zzz_yzz,   \
                             tr_z_zzz_zzz,   \
                             tr_z_zzzz_xx,   \
                             tr_z_zzzz_xxx,  \
                             tr_z_zzzz_xxy,  \
                             tr_z_zzzz_xxz,  \
                             tr_z_zzzz_xy,   \
                             tr_z_zzzz_xyy,  \
                             tr_z_zzzz_xyz,  \
                             tr_z_zzzz_xz,   \
                             tr_z_zzzz_xzz,  \
                             tr_z_zzzz_yy,   \
                             tr_z_zzzz_yyy,  \
                             tr_z_zzzz_yyz,  \
                             tr_z_zzzz_yz,   \
                             tr_z_zzzz_yzz,  \
                             tr_z_zzzz_zz,   \
                             tr_z_zzzz_zzz,  \
                             tr_z_zzzzz_xxx, \
                             tr_z_zzzzz_xxy, \
                             tr_z_zzzzz_xxz, \
                             tr_z_zzzzz_xyy, \
                             tr_z_zzzzz_xyz, \
                             tr_z_zzzzz_xzz, \
                             tr_z_zzzzz_yyy, \
                             tr_z_zzzzz_yyz, \
                             tr_z_zzzzz_yzz, \
                             tr_z_zzzzz_zzz, \
                             ts_zzzz_xxx,    \
                             ts_zzzz_xxy,    \
                             ts_zzzz_xxz,    \
                             ts_zzzz_xyy,    \
                             ts_zzzz_xyz,    \
                             ts_zzzz_xzz,    \
                             ts_zzzz_yyy,    \
                             ts_zzzz_yyz,    \
                             ts_zzzz_yzz,    \
                             ts_zzzz_zzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzzzz_xxx[i] = 4.0 * tr_z_zzz_xxx[i] * fe_0 + ts_zzzz_xxx[i] * fe_0 + tr_z_zzzz_xxx[i] * pa_z[i];

        tr_z_zzzzz_xxy[i] = 4.0 * tr_z_zzz_xxy[i] * fe_0 + ts_zzzz_xxy[i] * fe_0 + tr_z_zzzz_xxy[i] * pa_z[i];

        tr_z_zzzzz_xxz[i] = 4.0 * tr_z_zzz_xxz[i] * fe_0 + tr_z_zzzz_xx[i] * fe_0 + ts_zzzz_xxz[i] * fe_0 + tr_z_zzzz_xxz[i] * pa_z[i];

        tr_z_zzzzz_xyy[i] = 4.0 * tr_z_zzz_xyy[i] * fe_0 + ts_zzzz_xyy[i] * fe_0 + tr_z_zzzz_xyy[i] * pa_z[i];

        tr_z_zzzzz_xyz[i] = 4.0 * tr_z_zzz_xyz[i] * fe_0 + tr_z_zzzz_xy[i] * fe_0 + ts_zzzz_xyz[i] * fe_0 + tr_z_zzzz_xyz[i] * pa_z[i];

        tr_z_zzzzz_xzz[i] = 4.0 * tr_z_zzz_xzz[i] * fe_0 + 2.0 * tr_z_zzzz_xz[i] * fe_0 + ts_zzzz_xzz[i] * fe_0 + tr_z_zzzz_xzz[i] * pa_z[i];

        tr_z_zzzzz_yyy[i] = 4.0 * tr_z_zzz_yyy[i] * fe_0 + ts_zzzz_yyy[i] * fe_0 + tr_z_zzzz_yyy[i] * pa_z[i];

        tr_z_zzzzz_yyz[i] = 4.0 * tr_z_zzz_yyz[i] * fe_0 + tr_z_zzzz_yy[i] * fe_0 + ts_zzzz_yyz[i] * fe_0 + tr_z_zzzz_yyz[i] * pa_z[i];

        tr_z_zzzzz_yzz[i] = 4.0 * tr_z_zzz_yzz[i] * fe_0 + 2.0 * tr_z_zzzz_yz[i] * fe_0 + ts_zzzz_yzz[i] * fe_0 + tr_z_zzzz_yzz[i] * pa_z[i];

        tr_z_zzzzz_zzz[i] = 4.0 * tr_z_zzz_zzz[i] * fe_0 + 3.0 * tr_z_zzzz_zz[i] * fe_0 + ts_zzzz_zzz[i] * fe_0 + tr_z_zzzz_zzz[i] * pa_z[i];
    }
}

}  // namespace diprec
