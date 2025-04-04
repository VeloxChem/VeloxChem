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

#include "ElectricDipoleMomentumPrimRecFH.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_fh(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_fh,
                                      const size_t              idx_dip_ph,
                                      const size_t              idx_dip_dg,
                                      const size_t              idx_ovl_dh,
                                      const size_t              idx_dip_dh,
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

    // Set up components of auxiliary buffer : PH

    auto tr_x_x_xxxxx = pbuffer.data(idx_dip_ph);

    auto tr_x_x_xxxxy = pbuffer.data(idx_dip_ph + 1);

    auto tr_x_x_xxxxz = pbuffer.data(idx_dip_ph + 2);

    auto tr_x_x_xxxyy = pbuffer.data(idx_dip_ph + 3);

    auto tr_x_x_xxxyz = pbuffer.data(idx_dip_ph + 4);

    auto tr_x_x_xxxzz = pbuffer.data(idx_dip_ph + 5);

    auto tr_x_x_xxyyy = pbuffer.data(idx_dip_ph + 6);

    auto tr_x_x_xxyyz = pbuffer.data(idx_dip_ph + 7);

    auto tr_x_x_xxyzz = pbuffer.data(idx_dip_ph + 8);

    auto tr_x_x_xxzzz = pbuffer.data(idx_dip_ph + 9);

    auto tr_x_x_xyyyy = pbuffer.data(idx_dip_ph + 10);

    auto tr_x_x_xyyyz = pbuffer.data(idx_dip_ph + 11);

    auto tr_x_x_xyyzz = pbuffer.data(idx_dip_ph + 12);

    auto tr_x_x_xyzzz = pbuffer.data(idx_dip_ph + 13);

    auto tr_x_x_xzzzz = pbuffer.data(idx_dip_ph + 14);

    auto tr_x_x_yyyyy = pbuffer.data(idx_dip_ph + 15);

    auto tr_x_x_yyyyz = pbuffer.data(idx_dip_ph + 16);

    auto tr_x_x_yyyzz = pbuffer.data(idx_dip_ph + 17);

    auto tr_x_x_yyzzz = pbuffer.data(idx_dip_ph + 18);

    auto tr_x_x_yzzzz = pbuffer.data(idx_dip_ph + 19);

    auto tr_x_x_zzzzz = pbuffer.data(idx_dip_ph + 20);

    auto tr_x_y_xxxxx = pbuffer.data(idx_dip_ph + 21);

    auto tr_x_y_xxxxy = pbuffer.data(idx_dip_ph + 22);

    auto tr_x_y_xxxxz = pbuffer.data(idx_dip_ph + 23);

    auto tr_x_y_xxxyy = pbuffer.data(idx_dip_ph + 24);

    auto tr_x_y_xxxyz = pbuffer.data(idx_dip_ph + 25);

    auto tr_x_y_xxxzz = pbuffer.data(idx_dip_ph + 26);

    auto tr_x_y_xxyyy = pbuffer.data(idx_dip_ph + 27);

    auto tr_x_y_xxyyz = pbuffer.data(idx_dip_ph + 28);

    auto tr_x_y_xxyzz = pbuffer.data(idx_dip_ph + 29);

    auto tr_x_y_xxzzz = pbuffer.data(idx_dip_ph + 30);

    auto tr_x_y_xyyyy = pbuffer.data(idx_dip_ph + 31);

    auto tr_x_y_xyyyz = pbuffer.data(idx_dip_ph + 32);

    auto tr_x_y_xyyzz = pbuffer.data(idx_dip_ph + 33);

    auto tr_x_y_xyzzz = pbuffer.data(idx_dip_ph + 34);

    auto tr_x_y_xzzzz = pbuffer.data(idx_dip_ph + 35);

    auto tr_x_y_yyyyy = pbuffer.data(idx_dip_ph + 36);

    auto tr_x_y_yyyyz = pbuffer.data(idx_dip_ph + 37);

    auto tr_x_y_yyyzz = pbuffer.data(idx_dip_ph + 38);

    auto tr_x_y_yyzzz = pbuffer.data(idx_dip_ph + 39);

    auto tr_x_y_yzzzz = pbuffer.data(idx_dip_ph + 40);

    auto tr_x_y_zzzzz = pbuffer.data(idx_dip_ph + 41);

    auto tr_x_z_xxxxx = pbuffer.data(idx_dip_ph + 42);

    auto tr_x_z_xxxxy = pbuffer.data(idx_dip_ph + 43);

    auto tr_x_z_xxxxz = pbuffer.data(idx_dip_ph + 44);

    auto tr_x_z_xxxyy = pbuffer.data(idx_dip_ph + 45);

    auto tr_x_z_xxxyz = pbuffer.data(idx_dip_ph + 46);

    auto tr_x_z_xxxzz = pbuffer.data(idx_dip_ph + 47);

    auto tr_x_z_xxyyy = pbuffer.data(idx_dip_ph + 48);

    auto tr_x_z_xxyyz = pbuffer.data(idx_dip_ph + 49);

    auto tr_x_z_xxyzz = pbuffer.data(idx_dip_ph + 50);

    auto tr_x_z_xxzzz = pbuffer.data(idx_dip_ph + 51);

    auto tr_x_z_xyyyy = pbuffer.data(idx_dip_ph + 52);

    auto tr_x_z_xyyyz = pbuffer.data(idx_dip_ph + 53);

    auto tr_x_z_xyyzz = pbuffer.data(idx_dip_ph + 54);

    auto tr_x_z_xyzzz = pbuffer.data(idx_dip_ph + 55);

    auto tr_x_z_xzzzz = pbuffer.data(idx_dip_ph + 56);

    auto tr_x_z_yyyyy = pbuffer.data(idx_dip_ph + 57);

    auto tr_x_z_yyyyz = pbuffer.data(idx_dip_ph + 58);

    auto tr_x_z_yyyzz = pbuffer.data(idx_dip_ph + 59);

    auto tr_x_z_yyzzz = pbuffer.data(idx_dip_ph + 60);

    auto tr_x_z_yzzzz = pbuffer.data(idx_dip_ph + 61);

    auto tr_x_z_zzzzz = pbuffer.data(idx_dip_ph + 62);

    auto tr_y_x_xxxxx = pbuffer.data(idx_dip_ph + 63);

    auto tr_y_x_xxxxy = pbuffer.data(idx_dip_ph + 64);

    auto tr_y_x_xxxxz = pbuffer.data(idx_dip_ph + 65);

    auto tr_y_x_xxxyy = pbuffer.data(idx_dip_ph + 66);

    auto tr_y_x_xxxyz = pbuffer.data(idx_dip_ph + 67);

    auto tr_y_x_xxxzz = pbuffer.data(idx_dip_ph + 68);

    auto tr_y_x_xxyyy = pbuffer.data(idx_dip_ph + 69);

    auto tr_y_x_xxyyz = pbuffer.data(idx_dip_ph + 70);

    auto tr_y_x_xxyzz = pbuffer.data(idx_dip_ph + 71);

    auto tr_y_x_xxzzz = pbuffer.data(idx_dip_ph + 72);

    auto tr_y_x_xyyyy = pbuffer.data(idx_dip_ph + 73);

    auto tr_y_x_xyyyz = pbuffer.data(idx_dip_ph + 74);

    auto tr_y_x_xyyzz = pbuffer.data(idx_dip_ph + 75);

    auto tr_y_x_xyzzz = pbuffer.data(idx_dip_ph + 76);

    auto tr_y_x_xzzzz = pbuffer.data(idx_dip_ph + 77);

    auto tr_y_x_yyyyy = pbuffer.data(idx_dip_ph + 78);

    auto tr_y_x_yyyyz = pbuffer.data(idx_dip_ph + 79);

    auto tr_y_x_yyyzz = pbuffer.data(idx_dip_ph + 80);

    auto tr_y_x_yyzzz = pbuffer.data(idx_dip_ph + 81);

    auto tr_y_x_yzzzz = pbuffer.data(idx_dip_ph + 82);

    auto tr_y_x_zzzzz = pbuffer.data(idx_dip_ph + 83);

    auto tr_y_y_xxxxx = pbuffer.data(idx_dip_ph + 84);

    auto tr_y_y_xxxxy = pbuffer.data(idx_dip_ph + 85);

    auto tr_y_y_xxxxz = pbuffer.data(idx_dip_ph + 86);

    auto tr_y_y_xxxyy = pbuffer.data(idx_dip_ph + 87);

    auto tr_y_y_xxxyz = pbuffer.data(idx_dip_ph + 88);

    auto tr_y_y_xxxzz = pbuffer.data(idx_dip_ph + 89);

    auto tr_y_y_xxyyy = pbuffer.data(idx_dip_ph + 90);

    auto tr_y_y_xxyyz = pbuffer.data(idx_dip_ph + 91);

    auto tr_y_y_xxyzz = pbuffer.data(idx_dip_ph + 92);

    auto tr_y_y_xxzzz = pbuffer.data(idx_dip_ph + 93);

    auto tr_y_y_xyyyy = pbuffer.data(idx_dip_ph + 94);

    auto tr_y_y_xyyyz = pbuffer.data(idx_dip_ph + 95);

    auto tr_y_y_xyyzz = pbuffer.data(idx_dip_ph + 96);

    auto tr_y_y_xyzzz = pbuffer.data(idx_dip_ph + 97);

    auto tr_y_y_xzzzz = pbuffer.data(idx_dip_ph + 98);

    auto tr_y_y_yyyyy = pbuffer.data(idx_dip_ph + 99);

    auto tr_y_y_yyyyz = pbuffer.data(idx_dip_ph + 100);

    auto tr_y_y_yyyzz = pbuffer.data(idx_dip_ph + 101);

    auto tr_y_y_yyzzz = pbuffer.data(idx_dip_ph + 102);

    auto tr_y_y_yzzzz = pbuffer.data(idx_dip_ph + 103);

    auto tr_y_y_zzzzz = pbuffer.data(idx_dip_ph + 104);

    auto tr_y_z_xxxxx = pbuffer.data(idx_dip_ph + 105);

    auto tr_y_z_xxxxy = pbuffer.data(idx_dip_ph + 106);

    auto tr_y_z_xxxxz = pbuffer.data(idx_dip_ph + 107);

    auto tr_y_z_xxxyy = pbuffer.data(idx_dip_ph + 108);

    auto tr_y_z_xxxyz = pbuffer.data(idx_dip_ph + 109);

    auto tr_y_z_xxxzz = pbuffer.data(idx_dip_ph + 110);

    auto tr_y_z_xxyyy = pbuffer.data(idx_dip_ph + 111);

    auto tr_y_z_xxyyz = pbuffer.data(idx_dip_ph + 112);

    auto tr_y_z_xxyzz = pbuffer.data(idx_dip_ph + 113);

    auto tr_y_z_xxzzz = pbuffer.data(idx_dip_ph + 114);

    auto tr_y_z_xyyyy = pbuffer.data(idx_dip_ph + 115);

    auto tr_y_z_xyyyz = pbuffer.data(idx_dip_ph + 116);

    auto tr_y_z_xyyzz = pbuffer.data(idx_dip_ph + 117);

    auto tr_y_z_xyzzz = pbuffer.data(idx_dip_ph + 118);

    auto tr_y_z_xzzzz = pbuffer.data(idx_dip_ph + 119);

    auto tr_y_z_yyyyy = pbuffer.data(idx_dip_ph + 120);

    auto tr_y_z_yyyyz = pbuffer.data(idx_dip_ph + 121);

    auto tr_y_z_yyyzz = pbuffer.data(idx_dip_ph + 122);

    auto tr_y_z_yyzzz = pbuffer.data(idx_dip_ph + 123);

    auto tr_y_z_yzzzz = pbuffer.data(idx_dip_ph + 124);

    auto tr_y_z_zzzzz = pbuffer.data(idx_dip_ph + 125);

    auto tr_z_x_xxxxx = pbuffer.data(idx_dip_ph + 126);

    auto tr_z_x_xxxxy = pbuffer.data(idx_dip_ph + 127);

    auto tr_z_x_xxxxz = pbuffer.data(idx_dip_ph + 128);

    auto tr_z_x_xxxyy = pbuffer.data(idx_dip_ph + 129);

    auto tr_z_x_xxxyz = pbuffer.data(idx_dip_ph + 130);

    auto tr_z_x_xxxzz = pbuffer.data(idx_dip_ph + 131);

    auto tr_z_x_xxyyy = pbuffer.data(idx_dip_ph + 132);

    auto tr_z_x_xxyyz = pbuffer.data(idx_dip_ph + 133);

    auto tr_z_x_xxyzz = pbuffer.data(idx_dip_ph + 134);

    auto tr_z_x_xxzzz = pbuffer.data(idx_dip_ph + 135);

    auto tr_z_x_xyyyy = pbuffer.data(idx_dip_ph + 136);

    auto tr_z_x_xyyyz = pbuffer.data(idx_dip_ph + 137);

    auto tr_z_x_xyyzz = pbuffer.data(idx_dip_ph + 138);

    auto tr_z_x_xyzzz = pbuffer.data(idx_dip_ph + 139);

    auto tr_z_x_xzzzz = pbuffer.data(idx_dip_ph + 140);

    auto tr_z_x_yyyyy = pbuffer.data(idx_dip_ph + 141);

    auto tr_z_x_yyyyz = pbuffer.data(idx_dip_ph + 142);

    auto tr_z_x_yyyzz = pbuffer.data(idx_dip_ph + 143);

    auto tr_z_x_yyzzz = pbuffer.data(idx_dip_ph + 144);

    auto tr_z_x_yzzzz = pbuffer.data(idx_dip_ph + 145);

    auto tr_z_x_zzzzz = pbuffer.data(idx_dip_ph + 146);

    auto tr_z_y_xxxxx = pbuffer.data(idx_dip_ph + 147);

    auto tr_z_y_xxxxy = pbuffer.data(idx_dip_ph + 148);

    auto tr_z_y_xxxxz = pbuffer.data(idx_dip_ph + 149);

    auto tr_z_y_xxxyy = pbuffer.data(idx_dip_ph + 150);

    auto tr_z_y_xxxyz = pbuffer.data(idx_dip_ph + 151);

    auto tr_z_y_xxxzz = pbuffer.data(idx_dip_ph + 152);

    auto tr_z_y_xxyyy = pbuffer.data(idx_dip_ph + 153);

    auto tr_z_y_xxyyz = pbuffer.data(idx_dip_ph + 154);

    auto tr_z_y_xxyzz = pbuffer.data(idx_dip_ph + 155);

    auto tr_z_y_xxzzz = pbuffer.data(idx_dip_ph + 156);

    auto tr_z_y_xyyyy = pbuffer.data(idx_dip_ph + 157);

    auto tr_z_y_xyyyz = pbuffer.data(idx_dip_ph + 158);

    auto tr_z_y_xyyzz = pbuffer.data(idx_dip_ph + 159);

    auto tr_z_y_xyzzz = pbuffer.data(idx_dip_ph + 160);

    auto tr_z_y_xzzzz = pbuffer.data(idx_dip_ph + 161);

    auto tr_z_y_yyyyy = pbuffer.data(idx_dip_ph + 162);

    auto tr_z_y_yyyyz = pbuffer.data(idx_dip_ph + 163);

    auto tr_z_y_yyyzz = pbuffer.data(idx_dip_ph + 164);

    auto tr_z_y_yyzzz = pbuffer.data(idx_dip_ph + 165);

    auto tr_z_y_yzzzz = pbuffer.data(idx_dip_ph + 166);

    auto tr_z_y_zzzzz = pbuffer.data(idx_dip_ph + 167);

    auto tr_z_z_xxxxx = pbuffer.data(idx_dip_ph + 168);

    auto tr_z_z_xxxxy = pbuffer.data(idx_dip_ph + 169);

    auto tr_z_z_xxxxz = pbuffer.data(idx_dip_ph + 170);

    auto tr_z_z_xxxyy = pbuffer.data(idx_dip_ph + 171);

    auto tr_z_z_xxxyz = pbuffer.data(idx_dip_ph + 172);

    auto tr_z_z_xxxzz = pbuffer.data(idx_dip_ph + 173);

    auto tr_z_z_xxyyy = pbuffer.data(idx_dip_ph + 174);

    auto tr_z_z_xxyyz = pbuffer.data(idx_dip_ph + 175);

    auto tr_z_z_xxyzz = pbuffer.data(idx_dip_ph + 176);

    auto tr_z_z_xxzzz = pbuffer.data(idx_dip_ph + 177);

    auto tr_z_z_xyyyy = pbuffer.data(idx_dip_ph + 178);

    auto tr_z_z_xyyyz = pbuffer.data(idx_dip_ph + 179);

    auto tr_z_z_xyyzz = pbuffer.data(idx_dip_ph + 180);

    auto tr_z_z_xyzzz = pbuffer.data(idx_dip_ph + 181);

    auto tr_z_z_xzzzz = pbuffer.data(idx_dip_ph + 182);

    auto tr_z_z_yyyyy = pbuffer.data(idx_dip_ph + 183);

    auto tr_z_z_yyyyz = pbuffer.data(idx_dip_ph + 184);

    auto tr_z_z_yyyzz = pbuffer.data(idx_dip_ph + 185);

    auto tr_z_z_yyzzz = pbuffer.data(idx_dip_ph + 186);

    auto tr_z_z_yzzzz = pbuffer.data(idx_dip_ph + 187);

    auto tr_z_z_zzzzz = pbuffer.data(idx_dip_ph + 188);

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

    auto tr_x_xz_xxxz = pbuffer.data(idx_dip_dg + 32);

    auto tr_x_xz_xxyz = pbuffer.data(idx_dip_dg + 34);

    auto tr_x_xz_xxzz = pbuffer.data(idx_dip_dg + 35);

    auto tr_x_xz_xyyz = pbuffer.data(idx_dip_dg + 37);

    auto tr_x_xz_xyzz = pbuffer.data(idx_dip_dg + 38);

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

    auto tr_y_yz_xxxz = pbuffer.data(idx_dip_dg + 152);

    auto tr_y_yz_xxyz = pbuffer.data(idx_dip_dg + 154);

    auto tr_y_yz_xxzz = pbuffer.data(idx_dip_dg + 155);

    auto tr_y_yz_xyyz = pbuffer.data(idx_dip_dg + 157);

    auto tr_y_yz_xyzz = pbuffer.data(idx_dip_dg + 158);

    auto tr_y_yz_xzzz = pbuffer.data(idx_dip_dg + 159);

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

    auto tr_z_xz_xxxz = pbuffer.data(idx_dip_dg + 212);

    auto tr_z_xz_xxyz = pbuffer.data(idx_dip_dg + 214);

    auto tr_z_xz_xxzz = pbuffer.data(idx_dip_dg + 215);

    auto tr_z_xz_xyyz = pbuffer.data(idx_dip_dg + 217);

    auto tr_z_xz_xyzz = pbuffer.data(idx_dip_dg + 218);

    auto tr_z_xz_xzzz = pbuffer.data(idx_dip_dg + 219);

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

    auto tr_z_yz_xxxy = pbuffer.data(idx_dip_dg + 241);

    auto tr_z_yz_xxxz = pbuffer.data(idx_dip_dg + 242);

    auto tr_z_yz_xxyy = pbuffer.data(idx_dip_dg + 243);

    auto tr_z_yz_xxyz = pbuffer.data(idx_dip_dg + 244);

    auto tr_z_yz_xxzz = pbuffer.data(idx_dip_dg + 245);

    auto tr_z_yz_xyyy = pbuffer.data(idx_dip_dg + 246);

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

    // Set up components of auxiliary buffer : DH

    auto ts_xx_xxxxx = pbuffer.data(idx_ovl_dh);

    auto ts_xx_xxxxy = pbuffer.data(idx_ovl_dh + 1);

    auto ts_xx_xxxxz = pbuffer.data(idx_ovl_dh + 2);

    auto ts_xx_xxxyy = pbuffer.data(idx_ovl_dh + 3);

    auto ts_xx_xxxyz = pbuffer.data(idx_ovl_dh + 4);

    auto ts_xx_xxxzz = pbuffer.data(idx_ovl_dh + 5);

    auto ts_xx_xxyyy = pbuffer.data(idx_ovl_dh + 6);

    auto ts_xx_xxyyz = pbuffer.data(idx_ovl_dh + 7);

    auto ts_xx_xxyzz = pbuffer.data(idx_ovl_dh + 8);

    auto ts_xx_xxzzz = pbuffer.data(idx_ovl_dh + 9);

    auto ts_xx_xyyyy = pbuffer.data(idx_ovl_dh + 10);

    auto ts_xx_xyyyz = pbuffer.data(idx_ovl_dh + 11);

    auto ts_xx_xyyzz = pbuffer.data(idx_ovl_dh + 12);

    auto ts_xx_xyzzz = pbuffer.data(idx_ovl_dh + 13);

    auto ts_xx_xzzzz = pbuffer.data(idx_ovl_dh + 14);

    auto ts_xx_yyyyy = pbuffer.data(idx_ovl_dh + 15);

    auto ts_xx_yyyyz = pbuffer.data(idx_ovl_dh + 16);

    auto ts_xx_yyyzz = pbuffer.data(idx_ovl_dh + 17);

    auto ts_xx_yyzzz = pbuffer.data(idx_ovl_dh + 18);

    auto ts_xx_yzzzz = pbuffer.data(idx_ovl_dh + 19);

    auto ts_xx_zzzzz = pbuffer.data(idx_ovl_dh + 20);

    auto ts_yy_xxxxx = pbuffer.data(idx_ovl_dh + 63);

    auto ts_yy_xxxxy = pbuffer.data(idx_ovl_dh + 64);

    auto ts_yy_xxxxz = pbuffer.data(idx_ovl_dh + 65);

    auto ts_yy_xxxyy = pbuffer.data(idx_ovl_dh + 66);

    auto ts_yy_xxxyz = pbuffer.data(idx_ovl_dh + 67);

    auto ts_yy_xxxzz = pbuffer.data(idx_ovl_dh + 68);

    auto ts_yy_xxyyy = pbuffer.data(idx_ovl_dh + 69);

    auto ts_yy_xxyyz = pbuffer.data(idx_ovl_dh + 70);

    auto ts_yy_xxyzz = pbuffer.data(idx_ovl_dh + 71);

    auto ts_yy_xxzzz = pbuffer.data(idx_ovl_dh + 72);

    auto ts_yy_xyyyy = pbuffer.data(idx_ovl_dh + 73);

    auto ts_yy_xyyyz = pbuffer.data(idx_ovl_dh + 74);

    auto ts_yy_xyyzz = pbuffer.data(idx_ovl_dh + 75);

    auto ts_yy_xyzzz = pbuffer.data(idx_ovl_dh + 76);

    auto ts_yy_xzzzz = pbuffer.data(idx_ovl_dh + 77);

    auto ts_yy_yyyyy = pbuffer.data(idx_ovl_dh + 78);

    auto ts_yy_yyyyz = pbuffer.data(idx_ovl_dh + 79);

    auto ts_yy_yyyzz = pbuffer.data(idx_ovl_dh + 80);

    auto ts_yy_yyzzz = pbuffer.data(idx_ovl_dh + 81);

    auto ts_yy_yzzzz = pbuffer.data(idx_ovl_dh + 82);

    auto ts_yy_zzzzz = pbuffer.data(idx_ovl_dh + 83);

    auto ts_yz_yyyyz = pbuffer.data(idx_ovl_dh + 100);

    auto ts_yz_yyyzz = pbuffer.data(idx_ovl_dh + 101);

    auto ts_yz_yyzzz = pbuffer.data(idx_ovl_dh + 102);

    auto ts_yz_yzzzz = pbuffer.data(idx_ovl_dh + 103);

    auto ts_zz_xxxxx = pbuffer.data(idx_ovl_dh + 105);

    auto ts_zz_xxxxy = pbuffer.data(idx_ovl_dh + 106);

    auto ts_zz_xxxxz = pbuffer.data(idx_ovl_dh + 107);

    auto ts_zz_xxxyy = pbuffer.data(idx_ovl_dh + 108);

    auto ts_zz_xxxyz = pbuffer.data(idx_ovl_dh + 109);

    auto ts_zz_xxxzz = pbuffer.data(idx_ovl_dh + 110);

    auto ts_zz_xxyyy = pbuffer.data(idx_ovl_dh + 111);

    auto ts_zz_xxyyz = pbuffer.data(idx_ovl_dh + 112);

    auto ts_zz_xxyzz = pbuffer.data(idx_ovl_dh + 113);

    auto ts_zz_xxzzz = pbuffer.data(idx_ovl_dh + 114);

    auto ts_zz_xyyyy = pbuffer.data(idx_ovl_dh + 115);

    auto ts_zz_xyyyz = pbuffer.data(idx_ovl_dh + 116);

    auto ts_zz_xyyzz = pbuffer.data(idx_ovl_dh + 117);

    auto ts_zz_xyzzz = pbuffer.data(idx_ovl_dh + 118);

    auto ts_zz_xzzzz = pbuffer.data(idx_ovl_dh + 119);

    auto ts_zz_yyyyy = pbuffer.data(idx_ovl_dh + 120);

    auto ts_zz_yyyyz = pbuffer.data(idx_ovl_dh + 121);

    auto ts_zz_yyyzz = pbuffer.data(idx_ovl_dh + 122);

    auto ts_zz_yyzzz = pbuffer.data(idx_ovl_dh + 123);

    auto ts_zz_yzzzz = pbuffer.data(idx_ovl_dh + 124);

    auto ts_zz_zzzzz = pbuffer.data(idx_ovl_dh + 125);

    // Set up components of auxiliary buffer : DH

    auto tr_x_xx_xxxxx = pbuffer.data(idx_dip_dh);

    auto tr_x_xx_xxxxy = pbuffer.data(idx_dip_dh + 1);

    auto tr_x_xx_xxxxz = pbuffer.data(idx_dip_dh + 2);

    auto tr_x_xx_xxxyy = pbuffer.data(idx_dip_dh + 3);

    auto tr_x_xx_xxxyz = pbuffer.data(idx_dip_dh + 4);

    auto tr_x_xx_xxxzz = pbuffer.data(idx_dip_dh + 5);

    auto tr_x_xx_xxyyy = pbuffer.data(idx_dip_dh + 6);

    auto tr_x_xx_xxyyz = pbuffer.data(idx_dip_dh + 7);

    auto tr_x_xx_xxyzz = pbuffer.data(idx_dip_dh + 8);

    auto tr_x_xx_xxzzz = pbuffer.data(idx_dip_dh + 9);

    auto tr_x_xx_xyyyy = pbuffer.data(idx_dip_dh + 10);

    auto tr_x_xx_xyyyz = pbuffer.data(idx_dip_dh + 11);

    auto tr_x_xx_xyyzz = pbuffer.data(idx_dip_dh + 12);

    auto tr_x_xx_xyzzz = pbuffer.data(idx_dip_dh + 13);

    auto tr_x_xx_xzzzz = pbuffer.data(idx_dip_dh + 14);

    auto tr_x_xx_yyyyy = pbuffer.data(idx_dip_dh + 15);

    auto tr_x_xx_yyyyz = pbuffer.data(idx_dip_dh + 16);

    auto tr_x_xx_yyyzz = pbuffer.data(idx_dip_dh + 17);

    auto tr_x_xx_yyzzz = pbuffer.data(idx_dip_dh + 18);

    auto tr_x_xx_yzzzz = pbuffer.data(idx_dip_dh + 19);

    auto tr_x_xx_zzzzz = pbuffer.data(idx_dip_dh + 20);

    auto tr_x_xy_xxxxx = pbuffer.data(idx_dip_dh + 21);

    auto tr_x_xy_xxxxy = pbuffer.data(idx_dip_dh + 22);

    auto tr_x_xy_xxxxz = pbuffer.data(idx_dip_dh + 23);

    auto tr_x_xy_xxxyy = pbuffer.data(idx_dip_dh + 24);

    auto tr_x_xy_xxxzz = pbuffer.data(idx_dip_dh + 26);

    auto tr_x_xy_xxyyy = pbuffer.data(idx_dip_dh + 27);

    auto tr_x_xy_xxzzz = pbuffer.data(idx_dip_dh + 30);

    auto tr_x_xy_xyyyy = pbuffer.data(idx_dip_dh + 31);

    auto tr_x_xy_xzzzz = pbuffer.data(idx_dip_dh + 35);

    auto tr_x_xy_yyyyy = pbuffer.data(idx_dip_dh + 36);

    auto tr_x_xz_xxxxx = pbuffer.data(idx_dip_dh + 42);

    auto tr_x_xz_xxxxy = pbuffer.data(idx_dip_dh + 43);

    auto tr_x_xz_xxxxz = pbuffer.data(idx_dip_dh + 44);

    auto tr_x_xz_xxxyy = pbuffer.data(idx_dip_dh + 45);

    auto tr_x_xz_xxxyz = pbuffer.data(idx_dip_dh + 46);

    auto tr_x_xz_xxxzz = pbuffer.data(idx_dip_dh + 47);

    auto tr_x_xz_xxyyy = pbuffer.data(idx_dip_dh + 48);

    auto tr_x_xz_xxyyz = pbuffer.data(idx_dip_dh + 49);

    auto tr_x_xz_xxyzz = pbuffer.data(idx_dip_dh + 50);

    auto tr_x_xz_xxzzz = pbuffer.data(idx_dip_dh + 51);

    auto tr_x_xz_xyyyy = pbuffer.data(idx_dip_dh + 52);

    auto tr_x_xz_xyyyz = pbuffer.data(idx_dip_dh + 53);

    auto tr_x_xz_xyyzz = pbuffer.data(idx_dip_dh + 54);

    auto tr_x_xz_xyzzz = pbuffer.data(idx_dip_dh + 55);

    auto tr_x_xz_xzzzz = pbuffer.data(idx_dip_dh + 56);

    auto tr_x_xz_zzzzz = pbuffer.data(idx_dip_dh + 62);

    auto tr_x_yy_xxxxx = pbuffer.data(idx_dip_dh + 63);

    auto tr_x_yy_xxxxy = pbuffer.data(idx_dip_dh + 64);

    auto tr_x_yy_xxxxz = pbuffer.data(idx_dip_dh + 65);

    auto tr_x_yy_xxxyy = pbuffer.data(idx_dip_dh + 66);

    auto tr_x_yy_xxxyz = pbuffer.data(idx_dip_dh + 67);

    auto tr_x_yy_xxxzz = pbuffer.data(idx_dip_dh + 68);

    auto tr_x_yy_xxyyy = pbuffer.data(idx_dip_dh + 69);

    auto tr_x_yy_xxyyz = pbuffer.data(idx_dip_dh + 70);

    auto tr_x_yy_xxyzz = pbuffer.data(idx_dip_dh + 71);

    auto tr_x_yy_xxzzz = pbuffer.data(idx_dip_dh + 72);

    auto tr_x_yy_xyyyy = pbuffer.data(idx_dip_dh + 73);

    auto tr_x_yy_xyyyz = pbuffer.data(idx_dip_dh + 74);

    auto tr_x_yy_xyyzz = pbuffer.data(idx_dip_dh + 75);

    auto tr_x_yy_xyzzz = pbuffer.data(idx_dip_dh + 76);

    auto tr_x_yy_xzzzz = pbuffer.data(idx_dip_dh + 77);

    auto tr_x_yy_yyyyy = pbuffer.data(idx_dip_dh + 78);

    auto tr_x_yy_yyyyz = pbuffer.data(idx_dip_dh + 79);

    auto tr_x_yy_yyyzz = pbuffer.data(idx_dip_dh + 80);

    auto tr_x_yy_yyzzz = pbuffer.data(idx_dip_dh + 81);

    auto tr_x_yy_yzzzz = pbuffer.data(idx_dip_dh + 82);

    auto tr_x_yy_zzzzz = pbuffer.data(idx_dip_dh + 83);

    auto tr_x_yz_xxxxz = pbuffer.data(idx_dip_dh + 86);

    auto tr_x_yz_xxxzz = pbuffer.data(idx_dip_dh + 89);

    auto tr_x_yz_xxzzz = pbuffer.data(idx_dip_dh + 93);

    auto tr_x_yz_xzzzz = pbuffer.data(idx_dip_dh + 98);

    auto tr_x_yz_yyyyz = pbuffer.data(idx_dip_dh + 100);

    auto tr_x_yz_yyyzz = pbuffer.data(idx_dip_dh + 101);

    auto tr_x_yz_yyzzz = pbuffer.data(idx_dip_dh + 102);

    auto tr_x_yz_yzzzz = pbuffer.data(idx_dip_dh + 103);

    auto tr_x_yz_zzzzz = pbuffer.data(idx_dip_dh + 104);

    auto tr_x_zz_xxxxx = pbuffer.data(idx_dip_dh + 105);

    auto tr_x_zz_xxxxy = pbuffer.data(idx_dip_dh + 106);

    auto tr_x_zz_xxxxz = pbuffer.data(idx_dip_dh + 107);

    auto tr_x_zz_xxxyy = pbuffer.data(idx_dip_dh + 108);

    auto tr_x_zz_xxxyz = pbuffer.data(idx_dip_dh + 109);

    auto tr_x_zz_xxxzz = pbuffer.data(idx_dip_dh + 110);

    auto tr_x_zz_xxyyy = pbuffer.data(idx_dip_dh + 111);

    auto tr_x_zz_xxyyz = pbuffer.data(idx_dip_dh + 112);

    auto tr_x_zz_xxyzz = pbuffer.data(idx_dip_dh + 113);

    auto tr_x_zz_xxzzz = pbuffer.data(idx_dip_dh + 114);

    auto tr_x_zz_xyyyy = pbuffer.data(idx_dip_dh + 115);

    auto tr_x_zz_xyyyz = pbuffer.data(idx_dip_dh + 116);

    auto tr_x_zz_xyyzz = pbuffer.data(idx_dip_dh + 117);

    auto tr_x_zz_xyzzz = pbuffer.data(idx_dip_dh + 118);

    auto tr_x_zz_xzzzz = pbuffer.data(idx_dip_dh + 119);

    auto tr_x_zz_yyyyy = pbuffer.data(idx_dip_dh + 120);

    auto tr_x_zz_yyyyz = pbuffer.data(idx_dip_dh + 121);

    auto tr_x_zz_yyyzz = pbuffer.data(idx_dip_dh + 122);

    auto tr_x_zz_yyzzz = pbuffer.data(idx_dip_dh + 123);

    auto tr_x_zz_yzzzz = pbuffer.data(idx_dip_dh + 124);

    auto tr_x_zz_zzzzz = pbuffer.data(idx_dip_dh + 125);

    auto tr_y_xx_xxxxx = pbuffer.data(idx_dip_dh + 126);

    auto tr_y_xx_xxxxy = pbuffer.data(idx_dip_dh + 127);

    auto tr_y_xx_xxxxz = pbuffer.data(idx_dip_dh + 128);

    auto tr_y_xx_xxxyy = pbuffer.data(idx_dip_dh + 129);

    auto tr_y_xx_xxxyz = pbuffer.data(idx_dip_dh + 130);

    auto tr_y_xx_xxxzz = pbuffer.data(idx_dip_dh + 131);

    auto tr_y_xx_xxyyy = pbuffer.data(idx_dip_dh + 132);

    auto tr_y_xx_xxyyz = pbuffer.data(idx_dip_dh + 133);

    auto tr_y_xx_xxyzz = pbuffer.data(idx_dip_dh + 134);

    auto tr_y_xx_xxzzz = pbuffer.data(idx_dip_dh + 135);

    auto tr_y_xx_xyyyy = pbuffer.data(idx_dip_dh + 136);

    auto tr_y_xx_xyyyz = pbuffer.data(idx_dip_dh + 137);

    auto tr_y_xx_xyyzz = pbuffer.data(idx_dip_dh + 138);

    auto tr_y_xx_xyzzz = pbuffer.data(idx_dip_dh + 139);

    auto tr_y_xx_xzzzz = pbuffer.data(idx_dip_dh + 140);

    auto tr_y_xx_yyyyy = pbuffer.data(idx_dip_dh + 141);

    auto tr_y_xx_yyyyz = pbuffer.data(idx_dip_dh + 142);

    auto tr_y_xx_yyyzz = pbuffer.data(idx_dip_dh + 143);

    auto tr_y_xx_yyzzz = pbuffer.data(idx_dip_dh + 144);

    auto tr_y_xx_yzzzz = pbuffer.data(idx_dip_dh + 145);

    auto tr_y_xx_zzzzz = pbuffer.data(idx_dip_dh + 146);

    auto tr_y_xy_xxxxx = pbuffer.data(idx_dip_dh + 147);

    auto tr_y_xy_xxxxy = pbuffer.data(idx_dip_dh + 148);

    auto tr_y_xy_xxxyy = pbuffer.data(idx_dip_dh + 150);

    auto tr_y_xy_xxxyz = pbuffer.data(idx_dip_dh + 151);

    auto tr_y_xy_xxyyy = pbuffer.data(idx_dip_dh + 153);

    auto tr_y_xy_xxyyz = pbuffer.data(idx_dip_dh + 154);

    auto tr_y_xy_xxyzz = pbuffer.data(idx_dip_dh + 155);

    auto tr_y_xy_xyyyy = pbuffer.data(idx_dip_dh + 157);

    auto tr_y_xy_xyyyz = pbuffer.data(idx_dip_dh + 158);

    auto tr_y_xy_xyyzz = pbuffer.data(idx_dip_dh + 159);

    auto tr_y_xy_xyzzz = pbuffer.data(idx_dip_dh + 160);

    auto tr_y_xy_yyyyy = pbuffer.data(idx_dip_dh + 162);

    auto tr_y_xy_yyyyz = pbuffer.data(idx_dip_dh + 163);

    auto tr_y_xy_yyyzz = pbuffer.data(idx_dip_dh + 164);

    auto tr_y_xy_yyzzz = pbuffer.data(idx_dip_dh + 165);

    auto tr_y_xy_yzzzz = pbuffer.data(idx_dip_dh + 166);

    auto tr_y_xy_zzzzz = pbuffer.data(idx_dip_dh + 167);

    auto tr_y_xz_yyyyz = pbuffer.data(idx_dip_dh + 184);

    auto tr_y_xz_yyyzz = pbuffer.data(idx_dip_dh + 185);

    auto tr_y_xz_yyzzz = pbuffer.data(idx_dip_dh + 186);

    auto tr_y_xz_yzzzz = pbuffer.data(idx_dip_dh + 187);

    auto tr_y_xz_zzzzz = pbuffer.data(idx_dip_dh + 188);

    auto tr_y_yy_xxxxx = pbuffer.data(idx_dip_dh + 189);

    auto tr_y_yy_xxxxy = pbuffer.data(idx_dip_dh + 190);

    auto tr_y_yy_xxxxz = pbuffer.data(idx_dip_dh + 191);

    auto tr_y_yy_xxxyy = pbuffer.data(idx_dip_dh + 192);

    auto tr_y_yy_xxxyz = pbuffer.data(idx_dip_dh + 193);

    auto tr_y_yy_xxxzz = pbuffer.data(idx_dip_dh + 194);

    auto tr_y_yy_xxyyy = pbuffer.data(idx_dip_dh + 195);

    auto tr_y_yy_xxyyz = pbuffer.data(idx_dip_dh + 196);

    auto tr_y_yy_xxyzz = pbuffer.data(idx_dip_dh + 197);

    auto tr_y_yy_xxzzz = pbuffer.data(idx_dip_dh + 198);

    auto tr_y_yy_xyyyy = pbuffer.data(idx_dip_dh + 199);

    auto tr_y_yy_xyyyz = pbuffer.data(idx_dip_dh + 200);

    auto tr_y_yy_xyyzz = pbuffer.data(idx_dip_dh + 201);

    auto tr_y_yy_xyzzz = pbuffer.data(idx_dip_dh + 202);

    auto tr_y_yy_xzzzz = pbuffer.data(idx_dip_dh + 203);

    auto tr_y_yy_yyyyy = pbuffer.data(idx_dip_dh + 204);

    auto tr_y_yy_yyyyz = pbuffer.data(idx_dip_dh + 205);

    auto tr_y_yy_yyyzz = pbuffer.data(idx_dip_dh + 206);

    auto tr_y_yy_yyzzz = pbuffer.data(idx_dip_dh + 207);

    auto tr_y_yy_yzzzz = pbuffer.data(idx_dip_dh + 208);

    auto tr_y_yy_zzzzz = pbuffer.data(idx_dip_dh + 209);

    auto tr_y_yz_xxxxy = pbuffer.data(idx_dip_dh + 211);

    auto tr_y_yz_xxxxz = pbuffer.data(idx_dip_dh + 212);

    auto tr_y_yz_xxxyy = pbuffer.data(idx_dip_dh + 213);

    auto tr_y_yz_xxxyz = pbuffer.data(idx_dip_dh + 214);

    auto tr_y_yz_xxxzz = pbuffer.data(idx_dip_dh + 215);

    auto tr_y_yz_xxyyy = pbuffer.data(idx_dip_dh + 216);

    auto tr_y_yz_xxyyz = pbuffer.data(idx_dip_dh + 217);

    auto tr_y_yz_xxyzz = pbuffer.data(idx_dip_dh + 218);

    auto tr_y_yz_xxzzz = pbuffer.data(idx_dip_dh + 219);

    auto tr_y_yz_xyyyy = pbuffer.data(idx_dip_dh + 220);

    auto tr_y_yz_xyyyz = pbuffer.data(idx_dip_dh + 221);

    auto tr_y_yz_xyyzz = pbuffer.data(idx_dip_dh + 222);

    auto tr_y_yz_xyzzz = pbuffer.data(idx_dip_dh + 223);

    auto tr_y_yz_xzzzz = pbuffer.data(idx_dip_dh + 224);

    auto tr_y_yz_yyyyy = pbuffer.data(idx_dip_dh + 225);

    auto tr_y_yz_yyyyz = pbuffer.data(idx_dip_dh + 226);

    auto tr_y_yz_yyyzz = pbuffer.data(idx_dip_dh + 227);

    auto tr_y_yz_yyzzz = pbuffer.data(idx_dip_dh + 228);

    auto tr_y_yz_yzzzz = pbuffer.data(idx_dip_dh + 229);

    auto tr_y_yz_zzzzz = pbuffer.data(idx_dip_dh + 230);

    auto tr_y_zz_xxxxx = pbuffer.data(idx_dip_dh + 231);

    auto tr_y_zz_xxxxy = pbuffer.data(idx_dip_dh + 232);

    auto tr_y_zz_xxxxz = pbuffer.data(idx_dip_dh + 233);

    auto tr_y_zz_xxxyy = pbuffer.data(idx_dip_dh + 234);

    auto tr_y_zz_xxxyz = pbuffer.data(idx_dip_dh + 235);

    auto tr_y_zz_xxxzz = pbuffer.data(idx_dip_dh + 236);

    auto tr_y_zz_xxyyy = pbuffer.data(idx_dip_dh + 237);

    auto tr_y_zz_xxyyz = pbuffer.data(idx_dip_dh + 238);

    auto tr_y_zz_xxyzz = pbuffer.data(idx_dip_dh + 239);

    auto tr_y_zz_xxzzz = pbuffer.data(idx_dip_dh + 240);

    auto tr_y_zz_xyyyy = pbuffer.data(idx_dip_dh + 241);

    auto tr_y_zz_xyyyz = pbuffer.data(idx_dip_dh + 242);

    auto tr_y_zz_xyyzz = pbuffer.data(idx_dip_dh + 243);

    auto tr_y_zz_xyzzz = pbuffer.data(idx_dip_dh + 244);

    auto tr_y_zz_xzzzz = pbuffer.data(idx_dip_dh + 245);

    auto tr_y_zz_yyyyy = pbuffer.data(idx_dip_dh + 246);

    auto tr_y_zz_yyyyz = pbuffer.data(idx_dip_dh + 247);

    auto tr_y_zz_yyyzz = pbuffer.data(idx_dip_dh + 248);

    auto tr_y_zz_yyzzz = pbuffer.data(idx_dip_dh + 249);

    auto tr_y_zz_yzzzz = pbuffer.data(idx_dip_dh + 250);

    auto tr_y_zz_zzzzz = pbuffer.data(idx_dip_dh + 251);

    auto tr_z_xx_xxxxx = pbuffer.data(idx_dip_dh + 252);

    auto tr_z_xx_xxxxy = pbuffer.data(idx_dip_dh + 253);

    auto tr_z_xx_xxxxz = pbuffer.data(idx_dip_dh + 254);

    auto tr_z_xx_xxxyy = pbuffer.data(idx_dip_dh + 255);

    auto tr_z_xx_xxxyz = pbuffer.data(idx_dip_dh + 256);

    auto tr_z_xx_xxxzz = pbuffer.data(idx_dip_dh + 257);

    auto tr_z_xx_xxyyy = pbuffer.data(idx_dip_dh + 258);

    auto tr_z_xx_xxyyz = pbuffer.data(idx_dip_dh + 259);

    auto tr_z_xx_xxyzz = pbuffer.data(idx_dip_dh + 260);

    auto tr_z_xx_xxzzz = pbuffer.data(idx_dip_dh + 261);

    auto tr_z_xx_xyyyy = pbuffer.data(idx_dip_dh + 262);

    auto tr_z_xx_xyyyz = pbuffer.data(idx_dip_dh + 263);

    auto tr_z_xx_xyyzz = pbuffer.data(idx_dip_dh + 264);

    auto tr_z_xx_xyzzz = pbuffer.data(idx_dip_dh + 265);

    auto tr_z_xx_xzzzz = pbuffer.data(idx_dip_dh + 266);

    auto tr_z_xx_yyyyy = pbuffer.data(idx_dip_dh + 267);

    auto tr_z_xx_yyyyz = pbuffer.data(idx_dip_dh + 268);

    auto tr_z_xx_yyyzz = pbuffer.data(idx_dip_dh + 269);

    auto tr_z_xx_yyzzz = pbuffer.data(idx_dip_dh + 270);

    auto tr_z_xx_yzzzz = pbuffer.data(idx_dip_dh + 271);

    auto tr_z_xx_zzzzz = pbuffer.data(idx_dip_dh + 272);

    auto tr_z_xy_yyyyy = pbuffer.data(idx_dip_dh + 288);

    auto tr_z_xy_yyyyz = pbuffer.data(idx_dip_dh + 289);

    auto tr_z_xy_yyyzz = pbuffer.data(idx_dip_dh + 290);

    auto tr_z_xy_yyzzz = pbuffer.data(idx_dip_dh + 291);

    auto tr_z_xy_yzzzz = pbuffer.data(idx_dip_dh + 292);

    auto tr_z_xz_xxxxx = pbuffer.data(idx_dip_dh + 294);

    auto tr_z_xz_xxxxz = pbuffer.data(idx_dip_dh + 296);

    auto tr_z_xz_xxxyz = pbuffer.data(idx_dip_dh + 298);

    auto tr_z_xz_xxxzz = pbuffer.data(idx_dip_dh + 299);

    auto tr_z_xz_xxyyz = pbuffer.data(idx_dip_dh + 301);

    auto tr_z_xz_xxyzz = pbuffer.data(idx_dip_dh + 302);

    auto tr_z_xz_xxzzz = pbuffer.data(idx_dip_dh + 303);

    auto tr_z_xz_xyyyz = pbuffer.data(idx_dip_dh + 305);

    auto tr_z_xz_xyyzz = pbuffer.data(idx_dip_dh + 306);

    auto tr_z_xz_xyzzz = pbuffer.data(idx_dip_dh + 307);

    auto tr_z_xz_xzzzz = pbuffer.data(idx_dip_dh + 308);

    auto tr_z_xz_yyyyy = pbuffer.data(idx_dip_dh + 309);

    auto tr_z_xz_yyyyz = pbuffer.data(idx_dip_dh + 310);

    auto tr_z_xz_yyyzz = pbuffer.data(idx_dip_dh + 311);

    auto tr_z_xz_yyzzz = pbuffer.data(idx_dip_dh + 312);

    auto tr_z_xz_yzzzz = pbuffer.data(idx_dip_dh + 313);

    auto tr_z_xz_zzzzz = pbuffer.data(idx_dip_dh + 314);

    auto tr_z_yy_xxxxx = pbuffer.data(idx_dip_dh + 315);

    auto tr_z_yy_xxxxy = pbuffer.data(idx_dip_dh + 316);

    auto tr_z_yy_xxxxz = pbuffer.data(idx_dip_dh + 317);

    auto tr_z_yy_xxxyy = pbuffer.data(idx_dip_dh + 318);

    auto tr_z_yy_xxxyz = pbuffer.data(idx_dip_dh + 319);

    auto tr_z_yy_xxxzz = pbuffer.data(idx_dip_dh + 320);

    auto tr_z_yy_xxyyy = pbuffer.data(idx_dip_dh + 321);

    auto tr_z_yy_xxyyz = pbuffer.data(idx_dip_dh + 322);

    auto tr_z_yy_xxyzz = pbuffer.data(idx_dip_dh + 323);

    auto tr_z_yy_xxzzz = pbuffer.data(idx_dip_dh + 324);

    auto tr_z_yy_xyyyy = pbuffer.data(idx_dip_dh + 325);

    auto tr_z_yy_xyyyz = pbuffer.data(idx_dip_dh + 326);

    auto tr_z_yy_xyyzz = pbuffer.data(idx_dip_dh + 327);

    auto tr_z_yy_xyzzz = pbuffer.data(idx_dip_dh + 328);

    auto tr_z_yy_xzzzz = pbuffer.data(idx_dip_dh + 329);

    auto tr_z_yy_yyyyy = pbuffer.data(idx_dip_dh + 330);

    auto tr_z_yy_yyyyz = pbuffer.data(idx_dip_dh + 331);

    auto tr_z_yy_yyyzz = pbuffer.data(idx_dip_dh + 332);

    auto tr_z_yy_yyzzz = pbuffer.data(idx_dip_dh + 333);

    auto tr_z_yy_yzzzz = pbuffer.data(idx_dip_dh + 334);

    auto tr_z_yy_zzzzz = pbuffer.data(idx_dip_dh + 335);

    auto tr_z_yz_xxxxx = pbuffer.data(idx_dip_dh + 336);

    auto tr_z_yz_xxxxy = pbuffer.data(idx_dip_dh + 337);

    auto tr_z_yz_xxxxz = pbuffer.data(idx_dip_dh + 338);

    auto tr_z_yz_xxxyy = pbuffer.data(idx_dip_dh + 339);

    auto tr_z_yz_xxxyz = pbuffer.data(idx_dip_dh + 340);

    auto tr_z_yz_xxxzz = pbuffer.data(idx_dip_dh + 341);

    auto tr_z_yz_xxyyy = pbuffer.data(idx_dip_dh + 342);

    auto tr_z_yz_xxyyz = pbuffer.data(idx_dip_dh + 343);

    auto tr_z_yz_xxyzz = pbuffer.data(idx_dip_dh + 344);

    auto tr_z_yz_xxzzz = pbuffer.data(idx_dip_dh + 345);

    auto tr_z_yz_xyyyy = pbuffer.data(idx_dip_dh + 346);

    auto tr_z_yz_xyyyz = pbuffer.data(idx_dip_dh + 347);

    auto tr_z_yz_xyyzz = pbuffer.data(idx_dip_dh + 348);

    auto tr_z_yz_xyzzz = pbuffer.data(idx_dip_dh + 349);

    auto tr_z_yz_xzzzz = pbuffer.data(idx_dip_dh + 350);

    auto tr_z_yz_yyyyy = pbuffer.data(idx_dip_dh + 351);

    auto tr_z_yz_yyyyz = pbuffer.data(idx_dip_dh + 352);

    auto tr_z_yz_yyyzz = pbuffer.data(idx_dip_dh + 353);

    auto tr_z_yz_yyzzz = pbuffer.data(idx_dip_dh + 354);

    auto tr_z_yz_yzzzz = pbuffer.data(idx_dip_dh + 355);

    auto tr_z_yz_zzzzz = pbuffer.data(idx_dip_dh + 356);

    auto tr_z_zz_xxxxx = pbuffer.data(idx_dip_dh + 357);

    auto tr_z_zz_xxxxy = pbuffer.data(idx_dip_dh + 358);

    auto tr_z_zz_xxxxz = pbuffer.data(idx_dip_dh + 359);

    auto tr_z_zz_xxxyy = pbuffer.data(idx_dip_dh + 360);

    auto tr_z_zz_xxxyz = pbuffer.data(idx_dip_dh + 361);

    auto tr_z_zz_xxxzz = pbuffer.data(idx_dip_dh + 362);

    auto tr_z_zz_xxyyy = pbuffer.data(idx_dip_dh + 363);

    auto tr_z_zz_xxyyz = pbuffer.data(idx_dip_dh + 364);

    auto tr_z_zz_xxyzz = pbuffer.data(idx_dip_dh + 365);

    auto tr_z_zz_xxzzz = pbuffer.data(idx_dip_dh + 366);

    auto tr_z_zz_xyyyy = pbuffer.data(idx_dip_dh + 367);

    auto tr_z_zz_xyyyz = pbuffer.data(idx_dip_dh + 368);

    auto tr_z_zz_xyyzz = pbuffer.data(idx_dip_dh + 369);

    auto tr_z_zz_xyzzz = pbuffer.data(idx_dip_dh + 370);

    auto tr_z_zz_xzzzz = pbuffer.data(idx_dip_dh + 371);

    auto tr_z_zz_yyyyy = pbuffer.data(idx_dip_dh + 372);

    auto tr_z_zz_yyyyz = pbuffer.data(idx_dip_dh + 373);

    auto tr_z_zz_yyyzz = pbuffer.data(idx_dip_dh + 374);

    auto tr_z_zz_yyzzz = pbuffer.data(idx_dip_dh + 375);

    auto tr_z_zz_yzzzz = pbuffer.data(idx_dip_dh + 376);

    auto tr_z_zz_zzzzz = pbuffer.data(idx_dip_dh + 377);

    // Set up 0-21 components of targeted buffer : FH

    auto tr_x_xxx_xxxxx = pbuffer.data(idx_dip_fh);

    auto tr_x_xxx_xxxxy = pbuffer.data(idx_dip_fh + 1);

    auto tr_x_xxx_xxxxz = pbuffer.data(idx_dip_fh + 2);

    auto tr_x_xxx_xxxyy = pbuffer.data(idx_dip_fh + 3);

    auto tr_x_xxx_xxxyz = pbuffer.data(idx_dip_fh + 4);

    auto tr_x_xxx_xxxzz = pbuffer.data(idx_dip_fh + 5);

    auto tr_x_xxx_xxyyy = pbuffer.data(idx_dip_fh + 6);

    auto tr_x_xxx_xxyyz = pbuffer.data(idx_dip_fh + 7);

    auto tr_x_xxx_xxyzz = pbuffer.data(idx_dip_fh + 8);

    auto tr_x_xxx_xxzzz = pbuffer.data(idx_dip_fh + 9);

    auto tr_x_xxx_xyyyy = pbuffer.data(idx_dip_fh + 10);

    auto tr_x_xxx_xyyyz = pbuffer.data(idx_dip_fh + 11);

    auto tr_x_xxx_xyyzz = pbuffer.data(idx_dip_fh + 12);

    auto tr_x_xxx_xyzzz = pbuffer.data(idx_dip_fh + 13);

    auto tr_x_xxx_xzzzz = pbuffer.data(idx_dip_fh + 14);

    auto tr_x_xxx_yyyyy = pbuffer.data(idx_dip_fh + 15);

    auto tr_x_xxx_yyyyz = pbuffer.data(idx_dip_fh + 16);

    auto tr_x_xxx_yyyzz = pbuffer.data(idx_dip_fh + 17);

    auto tr_x_xxx_yyzzz = pbuffer.data(idx_dip_fh + 18);

    auto tr_x_xxx_yzzzz = pbuffer.data(idx_dip_fh + 19);

    auto tr_x_xxx_zzzzz = pbuffer.data(idx_dip_fh + 20);

#pragma omp simd aligned(pa_x,               \
                             tr_x_x_xxxxx,   \
                             tr_x_x_xxxxy,   \
                             tr_x_x_xxxxz,   \
                             tr_x_x_xxxyy,   \
                             tr_x_x_xxxyz,   \
                             tr_x_x_xxxzz,   \
                             tr_x_x_xxyyy,   \
                             tr_x_x_xxyyz,   \
                             tr_x_x_xxyzz,   \
                             tr_x_x_xxzzz,   \
                             tr_x_x_xyyyy,   \
                             tr_x_x_xyyyz,   \
                             tr_x_x_xyyzz,   \
                             tr_x_x_xyzzz,   \
                             tr_x_x_xzzzz,   \
                             tr_x_x_yyyyy,   \
                             tr_x_x_yyyyz,   \
                             tr_x_x_yyyzz,   \
                             tr_x_x_yyzzz,   \
                             tr_x_x_yzzzz,   \
                             tr_x_x_zzzzz,   \
                             tr_x_xx_xxxx,   \
                             tr_x_xx_xxxxx,  \
                             tr_x_xx_xxxxy,  \
                             tr_x_xx_xxxxz,  \
                             tr_x_xx_xxxy,   \
                             tr_x_xx_xxxyy,  \
                             tr_x_xx_xxxyz,  \
                             tr_x_xx_xxxz,   \
                             tr_x_xx_xxxzz,  \
                             tr_x_xx_xxyy,   \
                             tr_x_xx_xxyyy,  \
                             tr_x_xx_xxyyz,  \
                             tr_x_xx_xxyz,   \
                             tr_x_xx_xxyzz,  \
                             tr_x_xx_xxzz,   \
                             tr_x_xx_xxzzz,  \
                             tr_x_xx_xyyy,   \
                             tr_x_xx_xyyyy,  \
                             tr_x_xx_xyyyz,  \
                             tr_x_xx_xyyz,   \
                             tr_x_xx_xyyzz,  \
                             tr_x_xx_xyzz,   \
                             tr_x_xx_xyzzz,  \
                             tr_x_xx_xzzz,   \
                             tr_x_xx_xzzzz,  \
                             tr_x_xx_yyyy,   \
                             tr_x_xx_yyyyy,  \
                             tr_x_xx_yyyyz,  \
                             tr_x_xx_yyyz,   \
                             tr_x_xx_yyyzz,  \
                             tr_x_xx_yyzz,   \
                             tr_x_xx_yyzzz,  \
                             tr_x_xx_yzzz,   \
                             tr_x_xx_yzzzz,  \
                             tr_x_xx_zzzz,   \
                             tr_x_xx_zzzzz,  \
                             tr_x_xxx_xxxxx, \
                             tr_x_xxx_xxxxy, \
                             tr_x_xxx_xxxxz, \
                             tr_x_xxx_xxxyy, \
                             tr_x_xxx_xxxyz, \
                             tr_x_xxx_xxxzz, \
                             tr_x_xxx_xxyyy, \
                             tr_x_xxx_xxyyz, \
                             tr_x_xxx_xxyzz, \
                             tr_x_xxx_xxzzz, \
                             tr_x_xxx_xyyyy, \
                             tr_x_xxx_xyyyz, \
                             tr_x_xxx_xyyzz, \
                             tr_x_xxx_xyzzz, \
                             tr_x_xxx_xzzzz, \
                             tr_x_xxx_yyyyy, \
                             tr_x_xxx_yyyyz, \
                             tr_x_xxx_yyyzz, \
                             tr_x_xxx_yyzzz, \
                             tr_x_xxx_yzzzz, \
                             tr_x_xxx_zzzzz, \
                             ts_xx_xxxxx,    \
                             ts_xx_xxxxy,    \
                             ts_xx_xxxxz,    \
                             ts_xx_xxxyy,    \
                             ts_xx_xxxyz,    \
                             ts_xx_xxxzz,    \
                             ts_xx_xxyyy,    \
                             ts_xx_xxyyz,    \
                             ts_xx_xxyzz,    \
                             ts_xx_xxzzz,    \
                             ts_xx_xyyyy,    \
                             ts_xx_xyyyz,    \
                             ts_xx_xyyzz,    \
                             ts_xx_xyzzz,    \
                             ts_xx_xzzzz,    \
                             ts_xx_yyyyy,    \
                             ts_xx_yyyyz,    \
                             ts_xx_yyyzz,    \
                             ts_xx_yyzzz,    \
                             ts_xx_yzzzz,    \
                             ts_xx_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxx_xxxxx[i] = 2.0 * tr_x_x_xxxxx[i] * fe_0 + 5.0 * tr_x_xx_xxxx[i] * fe_0 + ts_xx_xxxxx[i] * fe_0 + tr_x_xx_xxxxx[i] * pa_x[i];

        tr_x_xxx_xxxxy[i] = 2.0 * tr_x_x_xxxxy[i] * fe_0 + 4.0 * tr_x_xx_xxxy[i] * fe_0 + ts_xx_xxxxy[i] * fe_0 + tr_x_xx_xxxxy[i] * pa_x[i];

        tr_x_xxx_xxxxz[i] = 2.0 * tr_x_x_xxxxz[i] * fe_0 + 4.0 * tr_x_xx_xxxz[i] * fe_0 + ts_xx_xxxxz[i] * fe_0 + tr_x_xx_xxxxz[i] * pa_x[i];

        tr_x_xxx_xxxyy[i] = 2.0 * tr_x_x_xxxyy[i] * fe_0 + 3.0 * tr_x_xx_xxyy[i] * fe_0 + ts_xx_xxxyy[i] * fe_0 + tr_x_xx_xxxyy[i] * pa_x[i];

        tr_x_xxx_xxxyz[i] = 2.0 * tr_x_x_xxxyz[i] * fe_0 + 3.0 * tr_x_xx_xxyz[i] * fe_0 + ts_xx_xxxyz[i] * fe_0 + tr_x_xx_xxxyz[i] * pa_x[i];

        tr_x_xxx_xxxzz[i] = 2.0 * tr_x_x_xxxzz[i] * fe_0 + 3.0 * tr_x_xx_xxzz[i] * fe_0 + ts_xx_xxxzz[i] * fe_0 + tr_x_xx_xxxzz[i] * pa_x[i];

        tr_x_xxx_xxyyy[i] = 2.0 * tr_x_x_xxyyy[i] * fe_0 + 2.0 * tr_x_xx_xyyy[i] * fe_0 + ts_xx_xxyyy[i] * fe_0 + tr_x_xx_xxyyy[i] * pa_x[i];

        tr_x_xxx_xxyyz[i] = 2.0 * tr_x_x_xxyyz[i] * fe_0 + 2.0 * tr_x_xx_xyyz[i] * fe_0 + ts_xx_xxyyz[i] * fe_0 + tr_x_xx_xxyyz[i] * pa_x[i];

        tr_x_xxx_xxyzz[i] = 2.0 * tr_x_x_xxyzz[i] * fe_0 + 2.0 * tr_x_xx_xyzz[i] * fe_0 + ts_xx_xxyzz[i] * fe_0 + tr_x_xx_xxyzz[i] * pa_x[i];

        tr_x_xxx_xxzzz[i] = 2.0 * tr_x_x_xxzzz[i] * fe_0 + 2.0 * tr_x_xx_xzzz[i] * fe_0 + ts_xx_xxzzz[i] * fe_0 + tr_x_xx_xxzzz[i] * pa_x[i];

        tr_x_xxx_xyyyy[i] = 2.0 * tr_x_x_xyyyy[i] * fe_0 + tr_x_xx_yyyy[i] * fe_0 + ts_xx_xyyyy[i] * fe_0 + tr_x_xx_xyyyy[i] * pa_x[i];

        tr_x_xxx_xyyyz[i] = 2.0 * tr_x_x_xyyyz[i] * fe_0 + tr_x_xx_yyyz[i] * fe_0 + ts_xx_xyyyz[i] * fe_0 + tr_x_xx_xyyyz[i] * pa_x[i];

        tr_x_xxx_xyyzz[i] = 2.0 * tr_x_x_xyyzz[i] * fe_0 + tr_x_xx_yyzz[i] * fe_0 + ts_xx_xyyzz[i] * fe_0 + tr_x_xx_xyyzz[i] * pa_x[i];

        tr_x_xxx_xyzzz[i] = 2.0 * tr_x_x_xyzzz[i] * fe_0 + tr_x_xx_yzzz[i] * fe_0 + ts_xx_xyzzz[i] * fe_0 + tr_x_xx_xyzzz[i] * pa_x[i];

        tr_x_xxx_xzzzz[i] = 2.0 * tr_x_x_xzzzz[i] * fe_0 + tr_x_xx_zzzz[i] * fe_0 + ts_xx_xzzzz[i] * fe_0 + tr_x_xx_xzzzz[i] * pa_x[i];

        tr_x_xxx_yyyyy[i] = 2.0 * tr_x_x_yyyyy[i] * fe_0 + ts_xx_yyyyy[i] * fe_0 + tr_x_xx_yyyyy[i] * pa_x[i];

        tr_x_xxx_yyyyz[i] = 2.0 * tr_x_x_yyyyz[i] * fe_0 + ts_xx_yyyyz[i] * fe_0 + tr_x_xx_yyyyz[i] * pa_x[i];

        tr_x_xxx_yyyzz[i] = 2.0 * tr_x_x_yyyzz[i] * fe_0 + ts_xx_yyyzz[i] * fe_0 + tr_x_xx_yyyzz[i] * pa_x[i];

        tr_x_xxx_yyzzz[i] = 2.0 * tr_x_x_yyzzz[i] * fe_0 + ts_xx_yyzzz[i] * fe_0 + tr_x_xx_yyzzz[i] * pa_x[i];

        tr_x_xxx_yzzzz[i] = 2.0 * tr_x_x_yzzzz[i] * fe_0 + ts_xx_yzzzz[i] * fe_0 + tr_x_xx_yzzzz[i] * pa_x[i];

        tr_x_xxx_zzzzz[i] = 2.0 * tr_x_x_zzzzz[i] * fe_0 + ts_xx_zzzzz[i] * fe_0 + tr_x_xx_zzzzz[i] * pa_x[i];
    }

    // Set up 21-42 components of targeted buffer : FH

    auto tr_x_xxy_xxxxx = pbuffer.data(idx_dip_fh + 21);

    auto tr_x_xxy_xxxxy = pbuffer.data(idx_dip_fh + 22);

    auto tr_x_xxy_xxxxz = pbuffer.data(idx_dip_fh + 23);

    auto tr_x_xxy_xxxyy = pbuffer.data(idx_dip_fh + 24);

    auto tr_x_xxy_xxxyz = pbuffer.data(idx_dip_fh + 25);

    auto tr_x_xxy_xxxzz = pbuffer.data(idx_dip_fh + 26);

    auto tr_x_xxy_xxyyy = pbuffer.data(idx_dip_fh + 27);

    auto tr_x_xxy_xxyyz = pbuffer.data(idx_dip_fh + 28);

    auto tr_x_xxy_xxyzz = pbuffer.data(idx_dip_fh + 29);

    auto tr_x_xxy_xxzzz = pbuffer.data(idx_dip_fh + 30);

    auto tr_x_xxy_xyyyy = pbuffer.data(idx_dip_fh + 31);

    auto tr_x_xxy_xyyyz = pbuffer.data(idx_dip_fh + 32);

    auto tr_x_xxy_xyyzz = pbuffer.data(idx_dip_fh + 33);

    auto tr_x_xxy_xyzzz = pbuffer.data(idx_dip_fh + 34);

    auto tr_x_xxy_xzzzz = pbuffer.data(idx_dip_fh + 35);

    auto tr_x_xxy_yyyyy = pbuffer.data(idx_dip_fh + 36);

    auto tr_x_xxy_yyyyz = pbuffer.data(idx_dip_fh + 37);

    auto tr_x_xxy_yyyzz = pbuffer.data(idx_dip_fh + 38);

    auto tr_x_xxy_yyzzz = pbuffer.data(idx_dip_fh + 39);

    auto tr_x_xxy_yzzzz = pbuffer.data(idx_dip_fh + 40);

    auto tr_x_xxy_zzzzz = pbuffer.data(idx_dip_fh + 41);

#pragma omp simd aligned(pa_y,               \
                             tr_x_xx_xxxx,   \
                             tr_x_xx_xxxxx,  \
                             tr_x_xx_xxxxy,  \
                             tr_x_xx_xxxxz,  \
                             tr_x_xx_xxxy,   \
                             tr_x_xx_xxxyy,  \
                             tr_x_xx_xxxyz,  \
                             tr_x_xx_xxxz,   \
                             tr_x_xx_xxxzz,  \
                             tr_x_xx_xxyy,   \
                             tr_x_xx_xxyyy,  \
                             tr_x_xx_xxyyz,  \
                             tr_x_xx_xxyz,   \
                             tr_x_xx_xxyzz,  \
                             tr_x_xx_xxzz,   \
                             tr_x_xx_xxzzz,  \
                             tr_x_xx_xyyy,   \
                             tr_x_xx_xyyyy,  \
                             tr_x_xx_xyyyz,  \
                             tr_x_xx_xyyz,   \
                             tr_x_xx_xyyzz,  \
                             tr_x_xx_xyzz,   \
                             tr_x_xx_xyzzz,  \
                             tr_x_xx_xzzz,   \
                             tr_x_xx_xzzzz,  \
                             tr_x_xx_yyyy,   \
                             tr_x_xx_yyyyy,  \
                             tr_x_xx_yyyyz,  \
                             tr_x_xx_yyyz,   \
                             tr_x_xx_yyyzz,  \
                             tr_x_xx_yyzz,   \
                             tr_x_xx_yyzzz,  \
                             tr_x_xx_yzzz,   \
                             tr_x_xx_yzzzz,  \
                             tr_x_xx_zzzz,   \
                             tr_x_xx_zzzzz,  \
                             tr_x_xxy_xxxxx, \
                             tr_x_xxy_xxxxy, \
                             tr_x_xxy_xxxxz, \
                             tr_x_xxy_xxxyy, \
                             tr_x_xxy_xxxyz, \
                             tr_x_xxy_xxxzz, \
                             tr_x_xxy_xxyyy, \
                             tr_x_xxy_xxyyz, \
                             tr_x_xxy_xxyzz, \
                             tr_x_xxy_xxzzz, \
                             tr_x_xxy_xyyyy, \
                             tr_x_xxy_xyyyz, \
                             tr_x_xxy_xyyzz, \
                             tr_x_xxy_xyzzz, \
                             tr_x_xxy_xzzzz, \
                             tr_x_xxy_yyyyy, \
                             tr_x_xxy_yyyyz, \
                             tr_x_xxy_yyyzz, \
                             tr_x_xxy_yyzzz, \
                             tr_x_xxy_yzzzz, \
                             tr_x_xxy_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxy_xxxxx[i] = tr_x_xx_xxxxx[i] * pa_y[i];

        tr_x_xxy_xxxxy[i] = tr_x_xx_xxxx[i] * fe_0 + tr_x_xx_xxxxy[i] * pa_y[i];

        tr_x_xxy_xxxxz[i] = tr_x_xx_xxxxz[i] * pa_y[i];

        tr_x_xxy_xxxyy[i] = 2.0 * tr_x_xx_xxxy[i] * fe_0 + tr_x_xx_xxxyy[i] * pa_y[i];

        tr_x_xxy_xxxyz[i] = tr_x_xx_xxxz[i] * fe_0 + tr_x_xx_xxxyz[i] * pa_y[i];

        tr_x_xxy_xxxzz[i] = tr_x_xx_xxxzz[i] * pa_y[i];

        tr_x_xxy_xxyyy[i] = 3.0 * tr_x_xx_xxyy[i] * fe_0 + tr_x_xx_xxyyy[i] * pa_y[i];

        tr_x_xxy_xxyyz[i] = 2.0 * tr_x_xx_xxyz[i] * fe_0 + tr_x_xx_xxyyz[i] * pa_y[i];

        tr_x_xxy_xxyzz[i] = tr_x_xx_xxzz[i] * fe_0 + tr_x_xx_xxyzz[i] * pa_y[i];

        tr_x_xxy_xxzzz[i] = tr_x_xx_xxzzz[i] * pa_y[i];

        tr_x_xxy_xyyyy[i] = 4.0 * tr_x_xx_xyyy[i] * fe_0 + tr_x_xx_xyyyy[i] * pa_y[i];

        tr_x_xxy_xyyyz[i] = 3.0 * tr_x_xx_xyyz[i] * fe_0 + tr_x_xx_xyyyz[i] * pa_y[i];

        tr_x_xxy_xyyzz[i] = 2.0 * tr_x_xx_xyzz[i] * fe_0 + tr_x_xx_xyyzz[i] * pa_y[i];

        tr_x_xxy_xyzzz[i] = tr_x_xx_xzzz[i] * fe_0 + tr_x_xx_xyzzz[i] * pa_y[i];

        tr_x_xxy_xzzzz[i] = tr_x_xx_xzzzz[i] * pa_y[i];

        tr_x_xxy_yyyyy[i] = 5.0 * tr_x_xx_yyyy[i] * fe_0 + tr_x_xx_yyyyy[i] * pa_y[i];

        tr_x_xxy_yyyyz[i] = 4.0 * tr_x_xx_yyyz[i] * fe_0 + tr_x_xx_yyyyz[i] * pa_y[i];

        tr_x_xxy_yyyzz[i] = 3.0 * tr_x_xx_yyzz[i] * fe_0 + tr_x_xx_yyyzz[i] * pa_y[i];

        tr_x_xxy_yyzzz[i] = 2.0 * tr_x_xx_yzzz[i] * fe_0 + tr_x_xx_yyzzz[i] * pa_y[i];

        tr_x_xxy_yzzzz[i] = tr_x_xx_zzzz[i] * fe_0 + tr_x_xx_yzzzz[i] * pa_y[i];

        tr_x_xxy_zzzzz[i] = tr_x_xx_zzzzz[i] * pa_y[i];
    }

    // Set up 42-63 components of targeted buffer : FH

    auto tr_x_xxz_xxxxx = pbuffer.data(idx_dip_fh + 42);

    auto tr_x_xxz_xxxxy = pbuffer.data(idx_dip_fh + 43);

    auto tr_x_xxz_xxxxz = pbuffer.data(idx_dip_fh + 44);

    auto tr_x_xxz_xxxyy = pbuffer.data(idx_dip_fh + 45);

    auto tr_x_xxz_xxxyz = pbuffer.data(idx_dip_fh + 46);

    auto tr_x_xxz_xxxzz = pbuffer.data(idx_dip_fh + 47);

    auto tr_x_xxz_xxyyy = pbuffer.data(idx_dip_fh + 48);

    auto tr_x_xxz_xxyyz = pbuffer.data(idx_dip_fh + 49);

    auto tr_x_xxz_xxyzz = pbuffer.data(idx_dip_fh + 50);

    auto tr_x_xxz_xxzzz = pbuffer.data(idx_dip_fh + 51);

    auto tr_x_xxz_xyyyy = pbuffer.data(idx_dip_fh + 52);

    auto tr_x_xxz_xyyyz = pbuffer.data(idx_dip_fh + 53);

    auto tr_x_xxz_xyyzz = pbuffer.data(idx_dip_fh + 54);

    auto tr_x_xxz_xyzzz = pbuffer.data(idx_dip_fh + 55);

    auto tr_x_xxz_xzzzz = pbuffer.data(idx_dip_fh + 56);

    auto tr_x_xxz_yyyyy = pbuffer.data(idx_dip_fh + 57);

    auto tr_x_xxz_yyyyz = pbuffer.data(idx_dip_fh + 58);

    auto tr_x_xxz_yyyzz = pbuffer.data(idx_dip_fh + 59);

    auto tr_x_xxz_yyzzz = pbuffer.data(idx_dip_fh + 60);

    auto tr_x_xxz_yzzzz = pbuffer.data(idx_dip_fh + 61);

    auto tr_x_xxz_zzzzz = pbuffer.data(idx_dip_fh + 62);

#pragma omp simd aligned(pa_z,               \
                             tr_x_xx_xxxx,   \
                             tr_x_xx_xxxxx,  \
                             tr_x_xx_xxxxy,  \
                             tr_x_xx_xxxxz,  \
                             tr_x_xx_xxxy,   \
                             tr_x_xx_xxxyy,  \
                             tr_x_xx_xxxyz,  \
                             tr_x_xx_xxxz,   \
                             tr_x_xx_xxxzz,  \
                             tr_x_xx_xxyy,   \
                             tr_x_xx_xxyyy,  \
                             tr_x_xx_xxyyz,  \
                             tr_x_xx_xxyz,   \
                             tr_x_xx_xxyzz,  \
                             tr_x_xx_xxzz,   \
                             tr_x_xx_xxzzz,  \
                             tr_x_xx_xyyy,   \
                             tr_x_xx_xyyyy,  \
                             tr_x_xx_xyyyz,  \
                             tr_x_xx_xyyz,   \
                             tr_x_xx_xyyzz,  \
                             tr_x_xx_xyzz,   \
                             tr_x_xx_xyzzz,  \
                             tr_x_xx_xzzz,   \
                             tr_x_xx_xzzzz,  \
                             tr_x_xx_yyyy,   \
                             tr_x_xx_yyyyy,  \
                             tr_x_xx_yyyyz,  \
                             tr_x_xx_yyyz,   \
                             tr_x_xx_yyyzz,  \
                             tr_x_xx_yyzz,   \
                             tr_x_xx_yyzzz,  \
                             tr_x_xx_yzzz,   \
                             tr_x_xx_yzzzz,  \
                             tr_x_xx_zzzz,   \
                             tr_x_xx_zzzzz,  \
                             tr_x_xxz_xxxxx, \
                             tr_x_xxz_xxxxy, \
                             tr_x_xxz_xxxxz, \
                             tr_x_xxz_xxxyy, \
                             tr_x_xxz_xxxyz, \
                             tr_x_xxz_xxxzz, \
                             tr_x_xxz_xxyyy, \
                             tr_x_xxz_xxyyz, \
                             tr_x_xxz_xxyzz, \
                             tr_x_xxz_xxzzz, \
                             tr_x_xxz_xyyyy, \
                             tr_x_xxz_xyyyz, \
                             tr_x_xxz_xyyzz, \
                             tr_x_xxz_xyzzz, \
                             tr_x_xxz_xzzzz, \
                             tr_x_xxz_yyyyy, \
                             tr_x_xxz_yyyyz, \
                             tr_x_xxz_yyyzz, \
                             tr_x_xxz_yyzzz, \
                             tr_x_xxz_yzzzz, \
                             tr_x_xxz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xxz_xxxxx[i] = tr_x_xx_xxxxx[i] * pa_z[i];

        tr_x_xxz_xxxxy[i] = tr_x_xx_xxxxy[i] * pa_z[i];

        tr_x_xxz_xxxxz[i] = tr_x_xx_xxxx[i] * fe_0 + tr_x_xx_xxxxz[i] * pa_z[i];

        tr_x_xxz_xxxyy[i] = tr_x_xx_xxxyy[i] * pa_z[i];

        tr_x_xxz_xxxyz[i] = tr_x_xx_xxxy[i] * fe_0 + tr_x_xx_xxxyz[i] * pa_z[i];

        tr_x_xxz_xxxzz[i] = 2.0 * tr_x_xx_xxxz[i] * fe_0 + tr_x_xx_xxxzz[i] * pa_z[i];

        tr_x_xxz_xxyyy[i] = tr_x_xx_xxyyy[i] * pa_z[i];

        tr_x_xxz_xxyyz[i] = tr_x_xx_xxyy[i] * fe_0 + tr_x_xx_xxyyz[i] * pa_z[i];

        tr_x_xxz_xxyzz[i] = 2.0 * tr_x_xx_xxyz[i] * fe_0 + tr_x_xx_xxyzz[i] * pa_z[i];

        tr_x_xxz_xxzzz[i] = 3.0 * tr_x_xx_xxzz[i] * fe_0 + tr_x_xx_xxzzz[i] * pa_z[i];

        tr_x_xxz_xyyyy[i] = tr_x_xx_xyyyy[i] * pa_z[i];

        tr_x_xxz_xyyyz[i] = tr_x_xx_xyyy[i] * fe_0 + tr_x_xx_xyyyz[i] * pa_z[i];

        tr_x_xxz_xyyzz[i] = 2.0 * tr_x_xx_xyyz[i] * fe_0 + tr_x_xx_xyyzz[i] * pa_z[i];

        tr_x_xxz_xyzzz[i] = 3.0 * tr_x_xx_xyzz[i] * fe_0 + tr_x_xx_xyzzz[i] * pa_z[i];

        tr_x_xxz_xzzzz[i] = 4.0 * tr_x_xx_xzzz[i] * fe_0 + tr_x_xx_xzzzz[i] * pa_z[i];

        tr_x_xxz_yyyyy[i] = tr_x_xx_yyyyy[i] * pa_z[i];

        tr_x_xxz_yyyyz[i] = tr_x_xx_yyyy[i] * fe_0 + tr_x_xx_yyyyz[i] * pa_z[i];

        tr_x_xxz_yyyzz[i] = 2.0 * tr_x_xx_yyyz[i] * fe_0 + tr_x_xx_yyyzz[i] * pa_z[i];

        tr_x_xxz_yyzzz[i] = 3.0 * tr_x_xx_yyzz[i] * fe_0 + tr_x_xx_yyzzz[i] * pa_z[i];

        tr_x_xxz_yzzzz[i] = 4.0 * tr_x_xx_yzzz[i] * fe_0 + tr_x_xx_yzzzz[i] * pa_z[i];

        tr_x_xxz_zzzzz[i] = 5.0 * tr_x_xx_zzzz[i] * fe_0 + tr_x_xx_zzzzz[i] * pa_z[i];
    }

    // Set up 63-84 components of targeted buffer : FH

    auto tr_x_xyy_xxxxx = pbuffer.data(idx_dip_fh + 63);

    auto tr_x_xyy_xxxxy = pbuffer.data(idx_dip_fh + 64);

    auto tr_x_xyy_xxxxz = pbuffer.data(idx_dip_fh + 65);

    auto tr_x_xyy_xxxyy = pbuffer.data(idx_dip_fh + 66);

    auto tr_x_xyy_xxxyz = pbuffer.data(idx_dip_fh + 67);

    auto tr_x_xyy_xxxzz = pbuffer.data(idx_dip_fh + 68);

    auto tr_x_xyy_xxyyy = pbuffer.data(idx_dip_fh + 69);

    auto tr_x_xyy_xxyyz = pbuffer.data(idx_dip_fh + 70);

    auto tr_x_xyy_xxyzz = pbuffer.data(idx_dip_fh + 71);

    auto tr_x_xyy_xxzzz = pbuffer.data(idx_dip_fh + 72);

    auto tr_x_xyy_xyyyy = pbuffer.data(idx_dip_fh + 73);

    auto tr_x_xyy_xyyyz = pbuffer.data(idx_dip_fh + 74);

    auto tr_x_xyy_xyyzz = pbuffer.data(idx_dip_fh + 75);

    auto tr_x_xyy_xyzzz = pbuffer.data(idx_dip_fh + 76);

    auto tr_x_xyy_xzzzz = pbuffer.data(idx_dip_fh + 77);

    auto tr_x_xyy_yyyyy = pbuffer.data(idx_dip_fh + 78);

    auto tr_x_xyy_yyyyz = pbuffer.data(idx_dip_fh + 79);

    auto tr_x_xyy_yyyzz = pbuffer.data(idx_dip_fh + 80);

    auto tr_x_xyy_yyzzz = pbuffer.data(idx_dip_fh + 81);

    auto tr_x_xyy_yzzzz = pbuffer.data(idx_dip_fh + 82);

    auto tr_x_xyy_zzzzz = pbuffer.data(idx_dip_fh + 83);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_x_x_xxxxx,   \
                             tr_x_x_xxxxz,   \
                             tr_x_x_xxxzz,   \
                             tr_x_x_xxzzz,   \
                             tr_x_x_xzzzz,   \
                             tr_x_xy_xxxxx,  \
                             tr_x_xy_xxxxz,  \
                             tr_x_xy_xxxzz,  \
                             tr_x_xy_xxzzz,  \
                             tr_x_xy_xzzzz,  \
                             tr_x_xyy_xxxxx, \
                             tr_x_xyy_xxxxy, \
                             tr_x_xyy_xxxxz, \
                             tr_x_xyy_xxxyy, \
                             tr_x_xyy_xxxyz, \
                             tr_x_xyy_xxxzz, \
                             tr_x_xyy_xxyyy, \
                             tr_x_xyy_xxyyz, \
                             tr_x_xyy_xxyzz, \
                             tr_x_xyy_xxzzz, \
                             tr_x_xyy_xyyyy, \
                             tr_x_xyy_xyyyz, \
                             tr_x_xyy_xyyzz, \
                             tr_x_xyy_xyzzz, \
                             tr_x_xyy_xzzzz, \
                             tr_x_xyy_yyyyy, \
                             tr_x_xyy_yyyyz, \
                             tr_x_xyy_yyyzz, \
                             tr_x_xyy_yyzzz, \
                             tr_x_xyy_yzzzz, \
                             tr_x_xyy_zzzzz, \
                             tr_x_yy_xxxxy,  \
                             tr_x_yy_xxxy,   \
                             tr_x_yy_xxxyy,  \
                             tr_x_yy_xxxyz,  \
                             tr_x_yy_xxyy,   \
                             tr_x_yy_xxyyy,  \
                             tr_x_yy_xxyyz,  \
                             tr_x_yy_xxyz,   \
                             tr_x_yy_xxyzz,  \
                             tr_x_yy_xyyy,   \
                             tr_x_yy_xyyyy,  \
                             tr_x_yy_xyyyz,  \
                             tr_x_yy_xyyz,   \
                             tr_x_yy_xyyzz,  \
                             tr_x_yy_xyzz,   \
                             tr_x_yy_xyzzz,  \
                             tr_x_yy_yyyy,   \
                             tr_x_yy_yyyyy,  \
                             tr_x_yy_yyyyz,  \
                             tr_x_yy_yyyz,   \
                             tr_x_yy_yyyzz,  \
                             tr_x_yy_yyzz,   \
                             tr_x_yy_yyzzz,  \
                             tr_x_yy_yzzz,   \
                             tr_x_yy_yzzzz,  \
                             tr_x_yy_zzzzz,  \
                             ts_yy_xxxxy,    \
                             ts_yy_xxxyy,    \
                             ts_yy_xxxyz,    \
                             ts_yy_xxyyy,    \
                             ts_yy_xxyyz,    \
                             ts_yy_xxyzz,    \
                             ts_yy_xyyyy,    \
                             ts_yy_xyyyz,    \
                             ts_yy_xyyzz,    \
                             ts_yy_xyzzz,    \
                             ts_yy_yyyyy,    \
                             ts_yy_yyyyz,    \
                             ts_yy_yyyzz,    \
                             ts_yy_yyzzz,    \
                             ts_yy_yzzzz,    \
                             ts_yy_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyy_xxxxx[i] = tr_x_x_xxxxx[i] * fe_0 + tr_x_xy_xxxxx[i] * pa_y[i];

        tr_x_xyy_xxxxy[i] = 4.0 * tr_x_yy_xxxy[i] * fe_0 + ts_yy_xxxxy[i] * fe_0 + tr_x_yy_xxxxy[i] * pa_x[i];

        tr_x_xyy_xxxxz[i] = tr_x_x_xxxxz[i] * fe_0 + tr_x_xy_xxxxz[i] * pa_y[i];

        tr_x_xyy_xxxyy[i] = 3.0 * tr_x_yy_xxyy[i] * fe_0 + ts_yy_xxxyy[i] * fe_0 + tr_x_yy_xxxyy[i] * pa_x[i];

        tr_x_xyy_xxxyz[i] = 3.0 * tr_x_yy_xxyz[i] * fe_0 + ts_yy_xxxyz[i] * fe_0 + tr_x_yy_xxxyz[i] * pa_x[i];

        tr_x_xyy_xxxzz[i] = tr_x_x_xxxzz[i] * fe_0 + tr_x_xy_xxxzz[i] * pa_y[i];

        tr_x_xyy_xxyyy[i] = 2.0 * tr_x_yy_xyyy[i] * fe_0 + ts_yy_xxyyy[i] * fe_0 + tr_x_yy_xxyyy[i] * pa_x[i];

        tr_x_xyy_xxyyz[i] = 2.0 * tr_x_yy_xyyz[i] * fe_0 + ts_yy_xxyyz[i] * fe_0 + tr_x_yy_xxyyz[i] * pa_x[i];

        tr_x_xyy_xxyzz[i] = 2.0 * tr_x_yy_xyzz[i] * fe_0 + ts_yy_xxyzz[i] * fe_0 + tr_x_yy_xxyzz[i] * pa_x[i];

        tr_x_xyy_xxzzz[i] = tr_x_x_xxzzz[i] * fe_0 + tr_x_xy_xxzzz[i] * pa_y[i];

        tr_x_xyy_xyyyy[i] = tr_x_yy_yyyy[i] * fe_0 + ts_yy_xyyyy[i] * fe_0 + tr_x_yy_xyyyy[i] * pa_x[i];

        tr_x_xyy_xyyyz[i] = tr_x_yy_yyyz[i] * fe_0 + ts_yy_xyyyz[i] * fe_0 + tr_x_yy_xyyyz[i] * pa_x[i];

        tr_x_xyy_xyyzz[i] = tr_x_yy_yyzz[i] * fe_0 + ts_yy_xyyzz[i] * fe_0 + tr_x_yy_xyyzz[i] * pa_x[i];

        tr_x_xyy_xyzzz[i] = tr_x_yy_yzzz[i] * fe_0 + ts_yy_xyzzz[i] * fe_0 + tr_x_yy_xyzzz[i] * pa_x[i];

        tr_x_xyy_xzzzz[i] = tr_x_x_xzzzz[i] * fe_0 + tr_x_xy_xzzzz[i] * pa_y[i];

        tr_x_xyy_yyyyy[i] = ts_yy_yyyyy[i] * fe_0 + tr_x_yy_yyyyy[i] * pa_x[i];

        tr_x_xyy_yyyyz[i] = ts_yy_yyyyz[i] * fe_0 + tr_x_yy_yyyyz[i] * pa_x[i];

        tr_x_xyy_yyyzz[i] = ts_yy_yyyzz[i] * fe_0 + tr_x_yy_yyyzz[i] * pa_x[i];

        tr_x_xyy_yyzzz[i] = ts_yy_yyzzz[i] * fe_0 + tr_x_yy_yyzzz[i] * pa_x[i];

        tr_x_xyy_yzzzz[i] = ts_yy_yzzzz[i] * fe_0 + tr_x_yy_yzzzz[i] * pa_x[i];

        tr_x_xyy_zzzzz[i] = ts_yy_zzzzz[i] * fe_0 + tr_x_yy_zzzzz[i] * pa_x[i];
    }

    // Set up 84-105 components of targeted buffer : FH

    auto tr_x_xyz_xxxxx = pbuffer.data(idx_dip_fh + 84);

    auto tr_x_xyz_xxxxy = pbuffer.data(idx_dip_fh + 85);

    auto tr_x_xyz_xxxxz = pbuffer.data(idx_dip_fh + 86);

    auto tr_x_xyz_xxxyy = pbuffer.data(idx_dip_fh + 87);

    auto tr_x_xyz_xxxyz = pbuffer.data(idx_dip_fh + 88);

    auto tr_x_xyz_xxxzz = pbuffer.data(idx_dip_fh + 89);

    auto tr_x_xyz_xxyyy = pbuffer.data(idx_dip_fh + 90);

    auto tr_x_xyz_xxyyz = pbuffer.data(idx_dip_fh + 91);

    auto tr_x_xyz_xxyzz = pbuffer.data(idx_dip_fh + 92);

    auto tr_x_xyz_xxzzz = pbuffer.data(idx_dip_fh + 93);

    auto tr_x_xyz_xyyyy = pbuffer.data(idx_dip_fh + 94);

    auto tr_x_xyz_xyyyz = pbuffer.data(idx_dip_fh + 95);

    auto tr_x_xyz_xyyzz = pbuffer.data(idx_dip_fh + 96);

    auto tr_x_xyz_xyzzz = pbuffer.data(idx_dip_fh + 97);

    auto tr_x_xyz_xzzzz = pbuffer.data(idx_dip_fh + 98);

    auto tr_x_xyz_yyyyy = pbuffer.data(idx_dip_fh + 99);

    auto tr_x_xyz_yyyyz = pbuffer.data(idx_dip_fh + 100);

    auto tr_x_xyz_yyyzz = pbuffer.data(idx_dip_fh + 101);

    auto tr_x_xyz_yyzzz = pbuffer.data(idx_dip_fh + 102);

    auto tr_x_xyz_yzzzz = pbuffer.data(idx_dip_fh + 103);

    auto tr_x_xyz_zzzzz = pbuffer.data(idx_dip_fh + 104);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             pa_z,           \
                             tr_x_xy_xxxxy,  \
                             tr_x_xy_xxxyy,  \
                             tr_x_xy_xxyyy,  \
                             tr_x_xy_xyyyy,  \
                             tr_x_xy_yyyyy,  \
                             tr_x_xyz_xxxxx, \
                             tr_x_xyz_xxxxy, \
                             tr_x_xyz_xxxxz, \
                             tr_x_xyz_xxxyy, \
                             tr_x_xyz_xxxyz, \
                             tr_x_xyz_xxxzz, \
                             tr_x_xyz_xxyyy, \
                             tr_x_xyz_xxyyz, \
                             tr_x_xyz_xxyzz, \
                             tr_x_xyz_xxzzz, \
                             tr_x_xyz_xyyyy, \
                             tr_x_xyz_xyyyz, \
                             tr_x_xyz_xyyzz, \
                             tr_x_xyz_xyzzz, \
                             tr_x_xyz_xzzzz, \
                             tr_x_xyz_yyyyy, \
                             tr_x_xyz_yyyyz, \
                             tr_x_xyz_yyyzz, \
                             tr_x_xyz_yyzzz, \
                             tr_x_xyz_yzzzz, \
                             tr_x_xyz_zzzzz, \
                             tr_x_xz_xxxxx,  \
                             tr_x_xz_xxxxz,  \
                             tr_x_xz_xxxyz,  \
                             tr_x_xz_xxxz,   \
                             tr_x_xz_xxxzz,  \
                             tr_x_xz_xxyyz,  \
                             tr_x_xz_xxyz,   \
                             tr_x_xz_xxyzz,  \
                             tr_x_xz_xxzz,   \
                             tr_x_xz_xxzzz,  \
                             tr_x_xz_xyyyz,  \
                             tr_x_xz_xyyz,   \
                             tr_x_xz_xyyzz,  \
                             tr_x_xz_xyzz,   \
                             tr_x_xz_xyzzz,  \
                             tr_x_xz_xzzz,   \
                             tr_x_xz_xzzzz,  \
                             tr_x_xz_zzzzz,  \
                             tr_x_yz_yyyyz,  \
                             tr_x_yz_yyyzz,  \
                             tr_x_yz_yyzzz,  \
                             tr_x_yz_yzzzz,  \
                             ts_yz_yyyyz,    \
                             ts_yz_yyyzz,    \
                             ts_yz_yyzzz,    \
                             ts_yz_yzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xyz_xxxxx[i] = tr_x_xz_xxxxx[i] * pa_y[i];

        tr_x_xyz_xxxxy[i] = tr_x_xy_xxxxy[i] * pa_z[i];

        tr_x_xyz_xxxxz[i] = tr_x_xz_xxxxz[i] * pa_y[i];

        tr_x_xyz_xxxyy[i] = tr_x_xy_xxxyy[i] * pa_z[i];

        tr_x_xyz_xxxyz[i] = tr_x_xz_xxxz[i] * fe_0 + tr_x_xz_xxxyz[i] * pa_y[i];

        tr_x_xyz_xxxzz[i] = tr_x_xz_xxxzz[i] * pa_y[i];

        tr_x_xyz_xxyyy[i] = tr_x_xy_xxyyy[i] * pa_z[i];

        tr_x_xyz_xxyyz[i] = 2.0 * tr_x_xz_xxyz[i] * fe_0 + tr_x_xz_xxyyz[i] * pa_y[i];

        tr_x_xyz_xxyzz[i] = tr_x_xz_xxzz[i] * fe_0 + tr_x_xz_xxyzz[i] * pa_y[i];

        tr_x_xyz_xxzzz[i] = tr_x_xz_xxzzz[i] * pa_y[i];

        tr_x_xyz_xyyyy[i] = tr_x_xy_xyyyy[i] * pa_z[i];

        tr_x_xyz_xyyyz[i] = 3.0 * tr_x_xz_xyyz[i] * fe_0 + tr_x_xz_xyyyz[i] * pa_y[i];

        tr_x_xyz_xyyzz[i] = 2.0 * tr_x_xz_xyzz[i] * fe_0 + tr_x_xz_xyyzz[i] * pa_y[i];

        tr_x_xyz_xyzzz[i] = tr_x_xz_xzzz[i] * fe_0 + tr_x_xz_xyzzz[i] * pa_y[i];

        tr_x_xyz_xzzzz[i] = tr_x_xz_xzzzz[i] * pa_y[i];

        tr_x_xyz_yyyyy[i] = tr_x_xy_yyyyy[i] * pa_z[i];

        tr_x_xyz_yyyyz[i] = ts_yz_yyyyz[i] * fe_0 + tr_x_yz_yyyyz[i] * pa_x[i];

        tr_x_xyz_yyyzz[i] = ts_yz_yyyzz[i] * fe_0 + tr_x_yz_yyyzz[i] * pa_x[i];

        tr_x_xyz_yyzzz[i] = ts_yz_yyzzz[i] * fe_0 + tr_x_yz_yyzzz[i] * pa_x[i];

        tr_x_xyz_yzzzz[i] = ts_yz_yzzzz[i] * fe_0 + tr_x_yz_yzzzz[i] * pa_x[i];

        tr_x_xyz_zzzzz[i] = tr_x_xz_zzzzz[i] * pa_y[i];
    }

    // Set up 105-126 components of targeted buffer : FH

    auto tr_x_xzz_xxxxx = pbuffer.data(idx_dip_fh + 105);

    auto tr_x_xzz_xxxxy = pbuffer.data(idx_dip_fh + 106);

    auto tr_x_xzz_xxxxz = pbuffer.data(idx_dip_fh + 107);

    auto tr_x_xzz_xxxyy = pbuffer.data(idx_dip_fh + 108);

    auto tr_x_xzz_xxxyz = pbuffer.data(idx_dip_fh + 109);

    auto tr_x_xzz_xxxzz = pbuffer.data(idx_dip_fh + 110);

    auto tr_x_xzz_xxyyy = pbuffer.data(idx_dip_fh + 111);

    auto tr_x_xzz_xxyyz = pbuffer.data(idx_dip_fh + 112);

    auto tr_x_xzz_xxyzz = pbuffer.data(idx_dip_fh + 113);

    auto tr_x_xzz_xxzzz = pbuffer.data(idx_dip_fh + 114);

    auto tr_x_xzz_xyyyy = pbuffer.data(idx_dip_fh + 115);

    auto tr_x_xzz_xyyyz = pbuffer.data(idx_dip_fh + 116);

    auto tr_x_xzz_xyyzz = pbuffer.data(idx_dip_fh + 117);

    auto tr_x_xzz_xyzzz = pbuffer.data(idx_dip_fh + 118);

    auto tr_x_xzz_xzzzz = pbuffer.data(idx_dip_fh + 119);

    auto tr_x_xzz_yyyyy = pbuffer.data(idx_dip_fh + 120);

    auto tr_x_xzz_yyyyz = pbuffer.data(idx_dip_fh + 121);

    auto tr_x_xzz_yyyzz = pbuffer.data(idx_dip_fh + 122);

    auto tr_x_xzz_yyzzz = pbuffer.data(idx_dip_fh + 123);

    auto tr_x_xzz_yzzzz = pbuffer.data(idx_dip_fh + 124);

    auto tr_x_xzz_zzzzz = pbuffer.data(idx_dip_fh + 125);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_x_x_xxxxx,   \
                             tr_x_x_xxxxy,   \
                             tr_x_x_xxxyy,   \
                             tr_x_x_xxyyy,   \
                             tr_x_x_xyyyy,   \
                             tr_x_xz_xxxxx,  \
                             tr_x_xz_xxxxy,  \
                             tr_x_xz_xxxyy,  \
                             tr_x_xz_xxyyy,  \
                             tr_x_xz_xyyyy,  \
                             tr_x_xzz_xxxxx, \
                             tr_x_xzz_xxxxy, \
                             tr_x_xzz_xxxxz, \
                             tr_x_xzz_xxxyy, \
                             tr_x_xzz_xxxyz, \
                             tr_x_xzz_xxxzz, \
                             tr_x_xzz_xxyyy, \
                             tr_x_xzz_xxyyz, \
                             tr_x_xzz_xxyzz, \
                             tr_x_xzz_xxzzz, \
                             tr_x_xzz_xyyyy, \
                             tr_x_xzz_xyyyz, \
                             tr_x_xzz_xyyzz, \
                             tr_x_xzz_xyzzz, \
                             tr_x_xzz_xzzzz, \
                             tr_x_xzz_yyyyy, \
                             tr_x_xzz_yyyyz, \
                             tr_x_xzz_yyyzz, \
                             tr_x_xzz_yyzzz, \
                             tr_x_xzz_yzzzz, \
                             tr_x_xzz_zzzzz, \
                             tr_x_zz_xxxxz,  \
                             tr_x_zz_xxxyz,  \
                             tr_x_zz_xxxz,   \
                             tr_x_zz_xxxzz,  \
                             tr_x_zz_xxyyz,  \
                             tr_x_zz_xxyz,   \
                             tr_x_zz_xxyzz,  \
                             tr_x_zz_xxzz,   \
                             tr_x_zz_xxzzz,  \
                             tr_x_zz_xyyyz,  \
                             tr_x_zz_xyyz,   \
                             tr_x_zz_xyyzz,  \
                             tr_x_zz_xyzz,   \
                             tr_x_zz_xyzzz,  \
                             tr_x_zz_xzzz,   \
                             tr_x_zz_xzzzz,  \
                             tr_x_zz_yyyyy,  \
                             tr_x_zz_yyyyz,  \
                             tr_x_zz_yyyz,   \
                             tr_x_zz_yyyzz,  \
                             tr_x_zz_yyzz,   \
                             tr_x_zz_yyzzz,  \
                             tr_x_zz_yzzz,   \
                             tr_x_zz_yzzzz,  \
                             tr_x_zz_zzzz,   \
                             tr_x_zz_zzzzz,  \
                             ts_zz_xxxxz,    \
                             ts_zz_xxxyz,    \
                             ts_zz_xxxzz,    \
                             ts_zz_xxyyz,    \
                             ts_zz_xxyzz,    \
                             ts_zz_xxzzz,    \
                             ts_zz_xyyyz,    \
                             ts_zz_xyyzz,    \
                             ts_zz_xyzzz,    \
                             ts_zz_xzzzz,    \
                             ts_zz_yyyyy,    \
                             ts_zz_yyyyz,    \
                             ts_zz_yyyzz,    \
                             ts_zz_yyzzz,    \
                             ts_zz_yzzzz,    \
                             ts_zz_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xzz_xxxxx[i] = tr_x_x_xxxxx[i] * fe_0 + tr_x_xz_xxxxx[i] * pa_z[i];

        tr_x_xzz_xxxxy[i] = tr_x_x_xxxxy[i] * fe_0 + tr_x_xz_xxxxy[i] * pa_z[i];

        tr_x_xzz_xxxxz[i] = 4.0 * tr_x_zz_xxxz[i] * fe_0 + ts_zz_xxxxz[i] * fe_0 + tr_x_zz_xxxxz[i] * pa_x[i];

        tr_x_xzz_xxxyy[i] = tr_x_x_xxxyy[i] * fe_0 + tr_x_xz_xxxyy[i] * pa_z[i];

        tr_x_xzz_xxxyz[i] = 3.0 * tr_x_zz_xxyz[i] * fe_0 + ts_zz_xxxyz[i] * fe_0 + tr_x_zz_xxxyz[i] * pa_x[i];

        tr_x_xzz_xxxzz[i] = 3.0 * tr_x_zz_xxzz[i] * fe_0 + ts_zz_xxxzz[i] * fe_0 + tr_x_zz_xxxzz[i] * pa_x[i];

        tr_x_xzz_xxyyy[i] = tr_x_x_xxyyy[i] * fe_0 + tr_x_xz_xxyyy[i] * pa_z[i];

        tr_x_xzz_xxyyz[i] = 2.0 * tr_x_zz_xyyz[i] * fe_0 + ts_zz_xxyyz[i] * fe_0 + tr_x_zz_xxyyz[i] * pa_x[i];

        tr_x_xzz_xxyzz[i] = 2.0 * tr_x_zz_xyzz[i] * fe_0 + ts_zz_xxyzz[i] * fe_0 + tr_x_zz_xxyzz[i] * pa_x[i];

        tr_x_xzz_xxzzz[i] = 2.0 * tr_x_zz_xzzz[i] * fe_0 + ts_zz_xxzzz[i] * fe_0 + tr_x_zz_xxzzz[i] * pa_x[i];

        tr_x_xzz_xyyyy[i] = tr_x_x_xyyyy[i] * fe_0 + tr_x_xz_xyyyy[i] * pa_z[i];

        tr_x_xzz_xyyyz[i] = tr_x_zz_yyyz[i] * fe_0 + ts_zz_xyyyz[i] * fe_0 + tr_x_zz_xyyyz[i] * pa_x[i];

        tr_x_xzz_xyyzz[i] = tr_x_zz_yyzz[i] * fe_0 + ts_zz_xyyzz[i] * fe_0 + tr_x_zz_xyyzz[i] * pa_x[i];

        tr_x_xzz_xyzzz[i] = tr_x_zz_yzzz[i] * fe_0 + ts_zz_xyzzz[i] * fe_0 + tr_x_zz_xyzzz[i] * pa_x[i];

        tr_x_xzz_xzzzz[i] = tr_x_zz_zzzz[i] * fe_0 + ts_zz_xzzzz[i] * fe_0 + tr_x_zz_xzzzz[i] * pa_x[i];

        tr_x_xzz_yyyyy[i] = ts_zz_yyyyy[i] * fe_0 + tr_x_zz_yyyyy[i] * pa_x[i];

        tr_x_xzz_yyyyz[i] = ts_zz_yyyyz[i] * fe_0 + tr_x_zz_yyyyz[i] * pa_x[i];

        tr_x_xzz_yyyzz[i] = ts_zz_yyyzz[i] * fe_0 + tr_x_zz_yyyzz[i] * pa_x[i];

        tr_x_xzz_yyzzz[i] = ts_zz_yyzzz[i] * fe_0 + tr_x_zz_yyzzz[i] * pa_x[i];

        tr_x_xzz_yzzzz[i] = ts_zz_yzzzz[i] * fe_0 + tr_x_zz_yzzzz[i] * pa_x[i];

        tr_x_xzz_zzzzz[i] = ts_zz_zzzzz[i] * fe_0 + tr_x_zz_zzzzz[i] * pa_x[i];
    }

    // Set up 126-147 components of targeted buffer : FH

    auto tr_x_yyy_xxxxx = pbuffer.data(idx_dip_fh + 126);

    auto tr_x_yyy_xxxxy = pbuffer.data(idx_dip_fh + 127);

    auto tr_x_yyy_xxxxz = pbuffer.data(idx_dip_fh + 128);

    auto tr_x_yyy_xxxyy = pbuffer.data(idx_dip_fh + 129);

    auto tr_x_yyy_xxxyz = pbuffer.data(idx_dip_fh + 130);

    auto tr_x_yyy_xxxzz = pbuffer.data(idx_dip_fh + 131);

    auto tr_x_yyy_xxyyy = pbuffer.data(idx_dip_fh + 132);

    auto tr_x_yyy_xxyyz = pbuffer.data(idx_dip_fh + 133);

    auto tr_x_yyy_xxyzz = pbuffer.data(idx_dip_fh + 134);

    auto tr_x_yyy_xxzzz = pbuffer.data(idx_dip_fh + 135);

    auto tr_x_yyy_xyyyy = pbuffer.data(idx_dip_fh + 136);

    auto tr_x_yyy_xyyyz = pbuffer.data(idx_dip_fh + 137);

    auto tr_x_yyy_xyyzz = pbuffer.data(idx_dip_fh + 138);

    auto tr_x_yyy_xyzzz = pbuffer.data(idx_dip_fh + 139);

    auto tr_x_yyy_xzzzz = pbuffer.data(idx_dip_fh + 140);

    auto tr_x_yyy_yyyyy = pbuffer.data(idx_dip_fh + 141);

    auto tr_x_yyy_yyyyz = pbuffer.data(idx_dip_fh + 142);

    auto tr_x_yyy_yyyzz = pbuffer.data(idx_dip_fh + 143);

    auto tr_x_yyy_yyzzz = pbuffer.data(idx_dip_fh + 144);

    auto tr_x_yyy_yzzzz = pbuffer.data(idx_dip_fh + 145);

    auto tr_x_yyy_zzzzz = pbuffer.data(idx_dip_fh + 146);

#pragma omp simd aligned(pa_y,               \
                             tr_x_y_xxxxx,   \
                             tr_x_y_xxxxy,   \
                             tr_x_y_xxxxz,   \
                             tr_x_y_xxxyy,   \
                             tr_x_y_xxxyz,   \
                             tr_x_y_xxxzz,   \
                             tr_x_y_xxyyy,   \
                             tr_x_y_xxyyz,   \
                             tr_x_y_xxyzz,   \
                             tr_x_y_xxzzz,   \
                             tr_x_y_xyyyy,   \
                             tr_x_y_xyyyz,   \
                             tr_x_y_xyyzz,   \
                             tr_x_y_xyzzz,   \
                             tr_x_y_xzzzz,   \
                             tr_x_y_yyyyy,   \
                             tr_x_y_yyyyz,   \
                             tr_x_y_yyyzz,   \
                             tr_x_y_yyzzz,   \
                             tr_x_y_yzzzz,   \
                             tr_x_y_zzzzz,   \
                             tr_x_yy_xxxx,   \
                             tr_x_yy_xxxxx,  \
                             tr_x_yy_xxxxy,  \
                             tr_x_yy_xxxxz,  \
                             tr_x_yy_xxxy,   \
                             tr_x_yy_xxxyy,  \
                             tr_x_yy_xxxyz,  \
                             tr_x_yy_xxxz,   \
                             tr_x_yy_xxxzz,  \
                             tr_x_yy_xxyy,   \
                             tr_x_yy_xxyyy,  \
                             tr_x_yy_xxyyz,  \
                             tr_x_yy_xxyz,   \
                             tr_x_yy_xxyzz,  \
                             tr_x_yy_xxzz,   \
                             tr_x_yy_xxzzz,  \
                             tr_x_yy_xyyy,   \
                             tr_x_yy_xyyyy,  \
                             tr_x_yy_xyyyz,  \
                             tr_x_yy_xyyz,   \
                             tr_x_yy_xyyzz,  \
                             tr_x_yy_xyzz,   \
                             tr_x_yy_xyzzz,  \
                             tr_x_yy_xzzz,   \
                             tr_x_yy_xzzzz,  \
                             tr_x_yy_yyyy,   \
                             tr_x_yy_yyyyy,  \
                             tr_x_yy_yyyyz,  \
                             tr_x_yy_yyyz,   \
                             tr_x_yy_yyyzz,  \
                             tr_x_yy_yyzz,   \
                             tr_x_yy_yyzzz,  \
                             tr_x_yy_yzzz,   \
                             tr_x_yy_yzzzz,  \
                             tr_x_yy_zzzz,   \
                             tr_x_yy_zzzzz,  \
                             tr_x_yyy_xxxxx, \
                             tr_x_yyy_xxxxy, \
                             tr_x_yyy_xxxxz, \
                             tr_x_yyy_xxxyy, \
                             tr_x_yyy_xxxyz, \
                             tr_x_yyy_xxxzz, \
                             tr_x_yyy_xxyyy, \
                             tr_x_yyy_xxyyz, \
                             tr_x_yyy_xxyzz, \
                             tr_x_yyy_xxzzz, \
                             tr_x_yyy_xyyyy, \
                             tr_x_yyy_xyyyz, \
                             tr_x_yyy_xyyzz, \
                             tr_x_yyy_xyzzz, \
                             tr_x_yyy_xzzzz, \
                             tr_x_yyy_yyyyy, \
                             tr_x_yyy_yyyyz, \
                             tr_x_yyy_yyyzz, \
                             tr_x_yyy_yyzzz, \
                             tr_x_yyy_yzzzz, \
                             tr_x_yyy_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyy_xxxxx[i] = 2.0 * tr_x_y_xxxxx[i] * fe_0 + tr_x_yy_xxxxx[i] * pa_y[i];

        tr_x_yyy_xxxxy[i] = 2.0 * tr_x_y_xxxxy[i] * fe_0 + tr_x_yy_xxxx[i] * fe_0 + tr_x_yy_xxxxy[i] * pa_y[i];

        tr_x_yyy_xxxxz[i] = 2.0 * tr_x_y_xxxxz[i] * fe_0 + tr_x_yy_xxxxz[i] * pa_y[i];

        tr_x_yyy_xxxyy[i] = 2.0 * tr_x_y_xxxyy[i] * fe_0 + 2.0 * tr_x_yy_xxxy[i] * fe_0 + tr_x_yy_xxxyy[i] * pa_y[i];

        tr_x_yyy_xxxyz[i] = 2.0 * tr_x_y_xxxyz[i] * fe_0 + tr_x_yy_xxxz[i] * fe_0 + tr_x_yy_xxxyz[i] * pa_y[i];

        tr_x_yyy_xxxzz[i] = 2.0 * tr_x_y_xxxzz[i] * fe_0 + tr_x_yy_xxxzz[i] * pa_y[i];

        tr_x_yyy_xxyyy[i] = 2.0 * tr_x_y_xxyyy[i] * fe_0 + 3.0 * tr_x_yy_xxyy[i] * fe_0 + tr_x_yy_xxyyy[i] * pa_y[i];

        tr_x_yyy_xxyyz[i] = 2.0 * tr_x_y_xxyyz[i] * fe_0 + 2.0 * tr_x_yy_xxyz[i] * fe_0 + tr_x_yy_xxyyz[i] * pa_y[i];

        tr_x_yyy_xxyzz[i] = 2.0 * tr_x_y_xxyzz[i] * fe_0 + tr_x_yy_xxzz[i] * fe_0 + tr_x_yy_xxyzz[i] * pa_y[i];

        tr_x_yyy_xxzzz[i] = 2.0 * tr_x_y_xxzzz[i] * fe_0 + tr_x_yy_xxzzz[i] * pa_y[i];

        tr_x_yyy_xyyyy[i] = 2.0 * tr_x_y_xyyyy[i] * fe_0 + 4.0 * tr_x_yy_xyyy[i] * fe_0 + tr_x_yy_xyyyy[i] * pa_y[i];

        tr_x_yyy_xyyyz[i] = 2.0 * tr_x_y_xyyyz[i] * fe_0 + 3.0 * tr_x_yy_xyyz[i] * fe_0 + tr_x_yy_xyyyz[i] * pa_y[i];

        tr_x_yyy_xyyzz[i] = 2.0 * tr_x_y_xyyzz[i] * fe_0 + 2.0 * tr_x_yy_xyzz[i] * fe_0 + tr_x_yy_xyyzz[i] * pa_y[i];

        tr_x_yyy_xyzzz[i] = 2.0 * tr_x_y_xyzzz[i] * fe_0 + tr_x_yy_xzzz[i] * fe_0 + tr_x_yy_xyzzz[i] * pa_y[i];

        tr_x_yyy_xzzzz[i] = 2.0 * tr_x_y_xzzzz[i] * fe_0 + tr_x_yy_xzzzz[i] * pa_y[i];

        tr_x_yyy_yyyyy[i] = 2.0 * tr_x_y_yyyyy[i] * fe_0 + 5.0 * tr_x_yy_yyyy[i] * fe_0 + tr_x_yy_yyyyy[i] * pa_y[i];

        tr_x_yyy_yyyyz[i] = 2.0 * tr_x_y_yyyyz[i] * fe_0 + 4.0 * tr_x_yy_yyyz[i] * fe_0 + tr_x_yy_yyyyz[i] * pa_y[i];

        tr_x_yyy_yyyzz[i] = 2.0 * tr_x_y_yyyzz[i] * fe_0 + 3.0 * tr_x_yy_yyzz[i] * fe_0 + tr_x_yy_yyyzz[i] * pa_y[i];

        tr_x_yyy_yyzzz[i] = 2.0 * tr_x_y_yyzzz[i] * fe_0 + 2.0 * tr_x_yy_yzzz[i] * fe_0 + tr_x_yy_yyzzz[i] * pa_y[i];

        tr_x_yyy_yzzzz[i] = 2.0 * tr_x_y_yzzzz[i] * fe_0 + tr_x_yy_zzzz[i] * fe_0 + tr_x_yy_yzzzz[i] * pa_y[i];

        tr_x_yyy_zzzzz[i] = 2.0 * tr_x_y_zzzzz[i] * fe_0 + tr_x_yy_zzzzz[i] * pa_y[i];
    }

    // Set up 147-168 components of targeted buffer : FH

    auto tr_x_yyz_xxxxx = pbuffer.data(idx_dip_fh + 147);

    auto tr_x_yyz_xxxxy = pbuffer.data(idx_dip_fh + 148);

    auto tr_x_yyz_xxxxz = pbuffer.data(idx_dip_fh + 149);

    auto tr_x_yyz_xxxyy = pbuffer.data(idx_dip_fh + 150);

    auto tr_x_yyz_xxxyz = pbuffer.data(idx_dip_fh + 151);

    auto tr_x_yyz_xxxzz = pbuffer.data(idx_dip_fh + 152);

    auto tr_x_yyz_xxyyy = pbuffer.data(idx_dip_fh + 153);

    auto tr_x_yyz_xxyyz = pbuffer.data(idx_dip_fh + 154);

    auto tr_x_yyz_xxyzz = pbuffer.data(idx_dip_fh + 155);

    auto tr_x_yyz_xxzzz = pbuffer.data(idx_dip_fh + 156);

    auto tr_x_yyz_xyyyy = pbuffer.data(idx_dip_fh + 157);

    auto tr_x_yyz_xyyyz = pbuffer.data(idx_dip_fh + 158);

    auto tr_x_yyz_xyyzz = pbuffer.data(idx_dip_fh + 159);

    auto tr_x_yyz_xyzzz = pbuffer.data(idx_dip_fh + 160);

    auto tr_x_yyz_xzzzz = pbuffer.data(idx_dip_fh + 161);

    auto tr_x_yyz_yyyyy = pbuffer.data(idx_dip_fh + 162);

    auto tr_x_yyz_yyyyz = pbuffer.data(idx_dip_fh + 163);

    auto tr_x_yyz_yyyzz = pbuffer.data(idx_dip_fh + 164);

    auto tr_x_yyz_yyzzz = pbuffer.data(idx_dip_fh + 165);

    auto tr_x_yyz_yzzzz = pbuffer.data(idx_dip_fh + 166);

    auto tr_x_yyz_zzzzz = pbuffer.data(idx_dip_fh + 167);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_yy_xxxxx,  \
                             tr_x_yy_xxxxy,  \
                             tr_x_yy_xxxy,   \
                             tr_x_yy_xxxyy,  \
                             tr_x_yy_xxxyz,  \
                             tr_x_yy_xxyy,   \
                             tr_x_yy_xxyyy,  \
                             tr_x_yy_xxyyz,  \
                             tr_x_yy_xxyz,   \
                             tr_x_yy_xxyzz,  \
                             tr_x_yy_xyyy,   \
                             tr_x_yy_xyyyy,  \
                             tr_x_yy_xyyyz,  \
                             tr_x_yy_xyyz,   \
                             tr_x_yy_xyyzz,  \
                             tr_x_yy_xyzz,   \
                             tr_x_yy_xyzzz,  \
                             tr_x_yy_yyyy,   \
                             tr_x_yy_yyyyy,  \
                             tr_x_yy_yyyyz,  \
                             tr_x_yy_yyyz,   \
                             tr_x_yy_yyyzz,  \
                             tr_x_yy_yyzz,   \
                             tr_x_yy_yyzzz,  \
                             tr_x_yy_yzzz,   \
                             tr_x_yy_yzzzz,  \
                             tr_x_yyz_xxxxx, \
                             tr_x_yyz_xxxxy, \
                             tr_x_yyz_xxxxz, \
                             tr_x_yyz_xxxyy, \
                             tr_x_yyz_xxxyz, \
                             tr_x_yyz_xxxzz, \
                             tr_x_yyz_xxyyy, \
                             tr_x_yyz_xxyyz, \
                             tr_x_yyz_xxyzz, \
                             tr_x_yyz_xxzzz, \
                             tr_x_yyz_xyyyy, \
                             tr_x_yyz_xyyyz, \
                             tr_x_yyz_xyyzz, \
                             tr_x_yyz_xyzzz, \
                             tr_x_yyz_xzzzz, \
                             tr_x_yyz_yyyyy, \
                             tr_x_yyz_yyyyz, \
                             tr_x_yyz_yyyzz, \
                             tr_x_yyz_yyzzz, \
                             tr_x_yyz_yzzzz, \
                             tr_x_yyz_zzzzz, \
                             tr_x_yz_xxxxz,  \
                             tr_x_yz_xxxzz,  \
                             tr_x_yz_xxzzz,  \
                             tr_x_yz_xzzzz,  \
                             tr_x_yz_zzzzz,  \
                             tr_x_z_xxxxz,   \
                             tr_x_z_xxxzz,   \
                             tr_x_z_xxzzz,   \
                             tr_x_z_xzzzz,   \
                             tr_x_z_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yyz_xxxxx[i] = tr_x_yy_xxxxx[i] * pa_z[i];

        tr_x_yyz_xxxxy[i] = tr_x_yy_xxxxy[i] * pa_z[i];

        tr_x_yyz_xxxxz[i] = tr_x_z_xxxxz[i] * fe_0 + tr_x_yz_xxxxz[i] * pa_y[i];

        tr_x_yyz_xxxyy[i] = tr_x_yy_xxxyy[i] * pa_z[i];

        tr_x_yyz_xxxyz[i] = tr_x_yy_xxxy[i] * fe_0 + tr_x_yy_xxxyz[i] * pa_z[i];

        tr_x_yyz_xxxzz[i] = tr_x_z_xxxzz[i] * fe_0 + tr_x_yz_xxxzz[i] * pa_y[i];

        tr_x_yyz_xxyyy[i] = tr_x_yy_xxyyy[i] * pa_z[i];

        tr_x_yyz_xxyyz[i] = tr_x_yy_xxyy[i] * fe_0 + tr_x_yy_xxyyz[i] * pa_z[i];

        tr_x_yyz_xxyzz[i] = 2.0 * tr_x_yy_xxyz[i] * fe_0 + tr_x_yy_xxyzz[i] * pa_z[i];

        tr_x_yyz_xxzzz[i] = tr_x_z_xxzzz[i] * fe_0 + tr_x_yz_xxzzz[i] * pa_y[i];

        tr_x_yyz_xyyyy[i] = tr_x_yy_xyyyy[i] * pa_z[i];

        tr_x_yyz_xyyyz[i] = tr_x_yy_xyyy[i] * fe_0 + tr_x_yy_xyyyz[i] * pa_z[i];

        tr_x_yyz_xyyzz[i] = 2.0 * tr_x_yy_xyyz[i] * fe_0 + tr_x_yy_xyyzz[i] * pa_z[i];

        tr_x_yyz_xyzzz[i] = 3.0 * tr_x_yy_xyzz[i] * fe_0 + tr_x_yy_xyzzz[i] * pa_z[i];

        tr_x_yyz_xzzzz[i] = tr_x_z_xzzzz[i] * fe_0 + tr_x_yz_xzzzz[i] * pa_y[i];

        tr_x_yyz_yyyyy[i] = tr_x_yy_yyyyy[i] * pa_z[i];

        tr_x_yyz_yyyyz[i] = tr_x_yy_yyyy[i] * fe_0 + tr_x_yy_yyyyz[i] * pa_z[i];

        tr_x_yyz_yyyzz[i] = 2.0 * tr_x_yy_yyyz[i] * fe_0 + tr_x_yy_yyyzz[i] * pa_z[i];

        tr_x_yyz_yyzzz[i] = 3.0 * tr_x_yy_yyzz[i] * fe_0 + tr_x_yy_yyzzz[i] * pa_z[i];

        tr_x_yyz_yzzzz[i] = 4.0 * tr_x_yy_yzzz[i] * fe_0 + tr_x_yy_yzzzz[i] * pa_z[i];

        tr_x_yyz_zzzzz[i] = tr_x_z_zzzzz[i] * fe_0 + tr_x_yz_zzzzz[i] * pa_y[i];
    }

    // Set up 168-189 components of targeted buffer : FH

    auto tr_x_yzz_xxxxx = pbuffer.data(idx_dip_fh + 168);

    auto tr_x_yzz_xxxxy = pbuffer.data(idx_dip_fh + 169);

    auto tr_x_yzz_xxxxz = pbuffer.data(idx_dip_fh + 170);

    auto tr_x_yzz_xxxyy = pbuffer.data(idx_dip_fh + 171);

    auto tr_x_yzz_xxxyz = pbuffer.data(idx_dip_fh + 172);

    auto tr_x_yzz_xxxzz = pbuffer.data(idx_dip_fh + 173);

    auto tr_x_yzz_xxyyy = pbuffer.data(idx_dip_fh + 174);

    auto tr_x_yzz_xxyyz = pbuffer.data(idx_dip_fh + 175);

    auto tr_x_yzz_xxyzz = pbuffer.data(idx_dip_fh + 176);

    auto tr_x_yzz_xxzzz = pbuffer.data(idx_dip_fh + 177);

    auto tr_x_yzz_xyyyy = pbuffer.data(idx_dip_fh + 178);

    auto tr_x_yzz_xyyyz = pbuffer.data(idx_dip_fh + 179);

    auto tr_x_yzz_xyyzz = pbuffer.data(idx_dip_fh + 180);

    auto tr_x_yzz_xyzzz = pbuffer.data(idx_dip_fh + 181);

    auto tr_x_yzz_xzzzz = pbuffer.data(idx_dip_fh + 182);

    auto tr_x_yzz_yyyyy = pbuffer.data(idx_dip_fh + 183);

    auto tr_x_yzz_yyyyz = pbuffer.data(idx_dip_fh + 184);

    auto tr_x_yzz_yyyzz = pbuffer.data(idx_dip_fh + 185);

    auto tr_x_yzz_yyzzz = pbuffer.data(idx_dip_fh + 186);

    auto tr_x_yzz_yzzzz = pbuffer.data(idx_dip_fh + 187);

    auto tr_x_yzz_zzzzz = pbuffer.data(idx_dip_fh + 188);

#pragma omp simd aligned(pa_y,               \
                             tr_x_yzz_xxxxx, \
                             tr_x_yzz_xxxxy, \
                             tr_x_yzz_xxxxz, \
                             tr_x_yzz_xxxyy, \
                             tr_x_yzz_xxxyz, \
                             tr_x_yzz_xxxzz, \
                             tr_x_yzz_xxyyy, \
                             tr_x_yzz_xxyyz, \
                             tr_x_yzz_xxyzz, \
                             tr_x_yzz_xxzzz, \
                             tr_x_yzz_xyyyy, \
                             tr_x_yzz_xyyyz, \
                             tr_x_yzz_xyyzz, \
                             tr_x_yzz_xyzzz, \
                             tr_x_yzz_xzzzz, \
                             tr_x_yzz_yyyyy, \
                             tr_x_yzz_yyyyz, \
                             tr_x_yzz_yyyzz, \
                             tr_x_yzz_yyzzz, \
                             tr_x_yzz_yzzzz, \
                             tr_x_yzz_zzzzz, \
                             tr_x_zz_xxxx,   \
                             tr_x_zz_xxxxx,  \
                             tr_x_zz_xxxxy,  \
                             tr_x_zz_xxxxz,  \
                             tr_x_zz_xxxy,   \
                             tr_x_zz_xxxyy,  \
                             tr_x_zz_xxxyz,  \
                             tr_x_zz_xxxz,   \
                             tr_x_zz_xxxzz,  \
                             tr_x_zz_xxyy,   \
                             tr_x_zz_xxyyy,  \
                             tr_x_zz_xxyyz,  \
                             tr_x_zz_xxyz,   \
                             tr_x_zz_xxyzz,  \
                             tr_x_zz_xxzz,   \
                             tr_x_zz_xxzzz,  \
                             tr_x_zz_xyyy,   \
                             tr_x_zz_xyyyy,  \
                             tr_x_zz_xyyyz,  \
                             tr_x_zz_xyyz,   \
                             tr_x_zz_xyyzz,  \
                             tr_x_zz_xyzz,   \
                             tr_x_zz_xyzzz,  \
                             tr_x_zz_xzzz,   \
                             tr_x_zz_xzzzz,  \
                             tr_x_zz_yyyy,   \
                             tr_x_zz_yyyyy,  \
                             tr_x_zz_yyyyz,  \
                             tr_x_zz_yyyz,   \
                             tr_x_zz_yyyzz,  \
                             tr_x_zz_yyzz,   \
                             tr_x_zz_yyzzz,  \
                             tr_x_zz_yzzz,   \
                             tr_x_zz_yzzzz,  \
                             tr_x_zz_zzzz,   \
                             tr_x_zz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yzz_xxxxx[i] = tr_x_zz_xxxxx[i] * pa_y[i];

        tr_x_yzz_xxxxy[i] = tr_x_zz_xxxx[i] * fe_0 + tr_x_zz_xxxxy[i] * pa_y[i];

        tr_x_yzz_xxxxz[i] = tr_x_zz_xxxxz[i] * pa_y[i];

        tr_x_yzz_xxxyy[i] = 2.0 * tr_x_zz_xxxy[i] * fe_0 + tr_x_zz_xxxyy[i] * pa_y[i];

        tr_x_yzz_xxxyz[i] = tr_x_zz_xxxz[i] * fe_0 + tr_x_zz_xxxyz[i] * pa_y[i];

        tr_x_yzz_xxxzz[i] = tr_x_zz_xxxzz[i] * pa_y[i];

        tr_x_yzz_xxyyy[i] = 3.0 * tr_x_zz_xxyy[i] * fe_0 + tr_x_zz_xxyyy[i] * pa_y[i];

        tr_x_yzz_xxyyz[i] = 2.0 * tr_x_zz_xxyz[i] * fe_0 + tr_x_zz_xxyyz[i] * pa_y[i];

        tr_x_yzz_xxyzz[i] = tr_x_zz_xxzz[i] * fe_0 + tr_x_zz_xxyzz[i] * pa_y[i];

        tr_x_yzz_xxzzz[i] = tr_x_zz_xxzzz[i] * pa_y[i];

        tr_x_yzz_xyyyy[i] = 4.0 * tr_x_zz_xyyy[i] * fe_0 + tr_x_zz_xyyyy[i] * pa_y[i];

        tr_x_yzz_xyyyz[i] = 3.0 * tr_x_zz_xyyz[i] * fe_0 + tr_x_zz_xyyyz[i] * pa_y[i];

        tr_x_yzz_xyyzz[i] = 2.0 * tr_x_zz_xyzz[i] * fe_0 + tr_x_zz_xyyzz[i] * pa_y[i];

        tr_x_yzz_xyzzz[i] = tr_x_zz_xzzz[i] * fe_0 + tr_x_zz_xyzzz[i] * pa_y[i];

        tr_x_yzz_xzzzz[i] = tr_x_zz_xzzzz[i] * pa_y[i];

        tr_x_yzz_yyyyy[i] = 5.0 * tr_x_zz_yyyy[i] * fe_0 + tr_x_zz_yyyyy[i] * pa_y[i];

        tr_x_yzz_yyyyz[i] = 4.0 * tr_x_zz_yyyz[i] * fe_0 + tr_x_zz_yyyyz[i] * pa_y[i];

        tr_x_yzz_yyyzz[i] = 3.0 * tr_x_zz_yyzz[i] * fe_0 + tr_x_zz_yyyzz[i] * pa_y[i];

        tr_x_yzz_yyzzz[i] = 2.0 * tr_x_zz_yzzz[i] * fe_0 + tr_x_zz_yyzzz[i] * pa_y[i];

        tr_x_yzz_yzzzz[i] = tr_x_zz_zzzz[i] * fe_0 + tr_x_zz_yzzzz[i] * pa_y[i];

        tr_x_yzz_zzzzz[i] = tr_x_zz_zzzzz[i] * pa_y[i];
    }

    // Set up 189-210 components of targeted buffer : FH

    auto tr_x_zzz_xxxxx = pbuffer.data(idx_dip_fh + 189);

    auto tr_x_zzz_xxxxy = pbuffer.data(idx_dip_fh + 190);

    auto tr_x_zzz_xxxxz = pbuffer.data(idx_dip_fh + 191);

    auto tr_x_zzz_xxxyy = pbuffer.data(idx_dip_fh + 192);

    auto tr_x_zzz_xxxyz = pbuffer.data(idx_dip_fh + 193);

    auto tr_x_zzz_xxxzz = pbuffer.data(idx_dip_fh + 194);

    auto tr_x_zzz_xxyyy = pbuffer.data(idx_dip_fh + 195);

    auto tr_x_zzz_xxyyz = pbuffer.data(idx_dip_fh + 196);

    auto tr_x_zzz_xxyzz = pbuffer.data(idx_dip_fh + 197);

    auto tr_x_zzz_xxzzz = pbuffer.data(idx_dip_fh + 198);

    auto tr_x_zzz_xyyyy = pbuffer.data(idx_dip_fh + 199);

    auto tr_x_zzz_xyyyz = pbuffer.data(idx_dip_fh + 200);

    auto tr_x_zzz_xyyzz = pbuffer.data(idx_dip_fh + 201);

    auto tr_x_zzz_xyzzz = pbuffer.data(idx_dip_fh + 202);

    auto tr_x_zzz_xzzzz = pbuffer.data(idx_dip_fh + 203);

    auto tr_x_zzz_yyyyy = pbuffer.data(idx_dip_fh + 204);

    auto tr_x_zzz_yyyyz = pbuffer.data(idx_dip_fh + 205);

    auto tr_x_zzz_yyyzz = pbuffer.data(idx_dip_fh + 206);

    auto tr_x_zzz_yyzzz = pbuffer.data(idx_dip_fh + 207);

    auto tr_x_zzz_yzzzz = pbuffer.data(idx_dip_fh + 208);

    auto tr_x_zzz_zzzzz = pbuffer.data(idx_dip_fh + 209);

#pragma omp simd aligned(pa_z,               \
                             tr_x_z_xxxxx,   \
                             tr_x_z_xxxxy,   \
                             tr_x_z_xxxxz,   \
                             tr_x_z_xxxyy,   \
                             tr_x_z_xxxyz,   \
                             tr_x_z_xxxzz,   \
                             tr_x_z_xxyyy,   \
                             tr_x_z_xxyyz,   \
                             tr_x_z_xxyzz,   \
                             tr_x_z_xxzzz,   \
                             tr_x_z_xyyyy,   \
                             tr_x_z_xyyyz,   \
                             tr_x_z_xyyzz,   \
                             tr_x_z_xyzzz,   \
                             tr_x_z_xzzzz,   \
                             tr_x_z_yyyyy,   \
                             tr_x_z_yyyyz,   \
                             tr_x_z_yyyzz,   \
                             tr_x_z_yyzzz,   \
                             tr_x_z_yzzzz,   \
                             tr_x_z_zzzzz,   \
                             tr_x_zz_xxxx,   \
                             tr_x_zz_xxxxx,  \
                             tr_x_zz_xxxxy,  \
                             tr_x_zz_xxxxz,  \
                             tr_x_zz_xxxy,   \
                             tr_x_zz_xxxyy,  \
                             tr_x_zz_xxxyz,  \
                             tr_x_zz_xxxz,   \
                             tr_x_zz_xxxzz,  \
                             tr_x_zz_xxyy,   \
                             tr_x_zz_xxyyy,  \
                             tr_x_zz_xxyyz,  \
                             tr_x_zz_xxyz,   \
                             tr_x_zz_xxyzz,  \
                             tr_x_zz_xxzz,   \
                             tr_x_zz_xxzzz,  \
                             tr_x_zz_xyyy,   \
                             tr_x_zz_xyyyy,  \
                             tr_x_zz_xyyyz,  \
                             tr_x_zz_xyyz,   \
                             tr_x_zz_xyyzz,  \
                             tr_x_zz_xyzz,   \
                             tr_x_zz_xyzzz,  \
                             tr_x_zz_xzzz,   \
                             tr_x_zz_xzzzz,  \
                             tr_x_zz_yyyy,   \
                             tr_x_zz_yyyyy,  \
                             tr_x_zz_yyyyz,  \
                             tr_x_zz_yyyz,   \
                             tr_x_zz_yyyzz,  \
                             tr_x_zz_yyzz,   \
                             tr_x_zz_yyzzz,  \
                             tr_x_zz_yzzz,   \
                             tr_x_zz_yzzzz,  \
                             tr_x_zz_zzzz,   \
                             tr_x_zz_zzzzz,  \
                             tr_x_zzz_xxxxx, \
                             tr_x_zzz_xxxxy, \
                             tr_x_zzz_xxxxz, \
                             tr_x_zzz_xxxyy, \
                             tr_x_zzz_xxxyz, \
                             tr_x_zzz_xxxzz, \
                             tr_x_zzz_xxyyy, \
                             tr_x_zzz_xxyyz, \
                             tr_x_zzz_xxyzz, \
                             tr_x_zzz_xxzzz, \
                             tr_x_zzz_xyyyy, \
                             tr_x_zzz_xyyyz, \
                             tr_x_zzz_xyyzz, \
                             tr_x_zzz_xyzzz, \
                             tr_x_zzz_xzzzz, \
                             tr_x_zzz_yyyyy, \
                             tr_x_zzz_yyyyz, \
                             tr_x_zzz_yyyzz, \
                             tr_x_zzz_yyzzz, \
                             tr_x_zzz_yzzzz, \
                             tr_x_zzz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zzz_xxxxx[i] = 2.0 * tr_x_z_xxxxx[i] * fe_0 + tr_x_zz_xxxxx[i] * pa_z[i];

        tr_x_zzz_xxxxy[i] = 2.0 * tr_x_z_xxxxy[i] * fe_0 + tr_x_zz_xxxxy[i] * pa_z[i];

        tr_x_zzz_xxxxz[i] = 2.0 * tr_x_z_xxxxz[i] * fe_0 + tr_x_zz_xxxx[i] * fe_0 + tr_x_zz_xxxxz[i] * pa_z[i];

        tr_x_zzz_xxxyy[i] = 2.0 * tr_x_z_xxxyy[i] * fe_0 + tr_x_zz_xxxyy[i] * pa_z[i];

        tr_x_zzz_xxxyz[i] = 2.0 * tr_x_z_xxxyz[i] * fe_0 + tr_x_zz_xxxy[i] * fe_0 + tr_x_zz_xxxyz[i] * pa_z[i];

        tr_x_zzz_xxxzz[i] = 2.0 * tr_x_z_xxxzz[i] * fe_0 + 2.0 * tr_x_zz_xxxz[i] * fe_0 + tr_x_zz_xxxzz[i] * pa_z[i];

        tr_x_zzz_xxyyy[i] = 2.0 * tr_x_z_xxyyy[i] * fe_0 + tr_x_zz_xxyyy[i] * pa_z[i];

        tr_x_zzz_xxyyz[i] = 2.0 * tr_x_z_xxyyz[i] * fe_0 + tr_x_zz_xxyy[i] * fe_0 + tr_x_zz_xxyyz[i] * pa_z[i];

        tr_x_zzz_xxyzz[i] = 2.0 * tr_x_z_xxyzz[i] * fe_0 + 2.0 * tr_x_zz_xxyz[i] * fe_0 + tr_x_zz_xxyzz[i] * pa_z[i];

        tr_x_zzz_xxzzz[i] = 2.0 * tr_x_z_xxzzz[i] * fe_0 + 3.0 * tr_x_zz_xxzz[i] * fe_0 + tr_x_zz_xxzzz[i] * pa_z[i];

        tr_x_zzz_xyyyy[i] = 2.0 * tr_x_z_xyyyy[i] * fe_0 + tr_x_zz_xyyyy[i] * pa_z[i];

        tr_x_zzz_xyyyz[i] = 2.0 * tr_x_z_xyyyz[i] * fe_0 + tr_x_zz_xyyy[i] * fe_0 + tr_x_zz_xyyyz[i] * pa_z[i];

        tr_x_zzz_xyyzz[i] = 2.0 * tr_x_z_xyyzz[i] * fe_0 + 2.0 * tr_x_zz_xyyz[i] * fe_0 + tr_x_zz_xyyzz[i] * pa_z[i];

        tr_x_zzz_xyzzz[i] = 2.0 * tr_x_z_xyzzz[i] * fe_0 + 3.0 * tr_x_zz_xyzz[i] * fe_0 + tr_x_zz_xyzzz[i] * pa_z[i];

        tr_x_zzz_xzzzz[i] = 2.0 * tr_x_z_xzzzz[i] * fe_0 + 4.0 * tr_x_zz_xzzz[i] * fe_0 + tr_x_zz_xzzzz[i] * pa_z[i];

        tr_x_zzz_yyyyy[i] = 2.0 * tr_x_z_yyyyy[i] * fe_0 + tr_x_zz_yyyyy[i] * pa_z[i];

        tr_x_zzz_yyyyz[i] = 2.0 * tr_x_z_yyyyz[i] * fe_0 + tr_x_zz_yyyy[i] * fe_0 + tr_x_zz_yyyyz[i] * pa_z[i];

        tr_x_zzz_yyyzz[i] = 2.0 * tr_x_z_yyyzz[i] * fe_0 + 2.0 * tr_x_zz_yyyz[i] * fe_0 + tr_x_zz_yyyzz[i] * pa_z[i];

        tr_x_zzz_yyzzz[i] = 2.0 * tr_x_z_yyzzz[i] * fe_0 + 3.0 * tr_x_zz_yyzz[i] * fe_0 + tr_x_zz_yyzzz[i] * pa_z[i];

        tr_x_zzz_yzzzz[i] = 2.0 * tr_x_z_yzzzz[i] * fe_0 + 4.0 * tr_x_zz_yzzz[i] * fe_0 + tr_x_zz_yzzzz[i] * pa_z[i];

        tr_x_zzz_zzzzz[i] = 2.0 * tr_x_z_zzzzz[i] * fe_0 + 5.0 * tr_x_zz_zzzz[i] * fe_0 + tr_x_zz_zzzzz[i] * pa_z[i];
    }

    // Set up 210-231 components of targeted buffer : FH

    auto tr_y_xxx_xxxxx = pbuffer.data(idx_dip_fh + 210);

    auto tr_y_xxx_xxxxy = pbuffer.data(idx_dip_fh + 211);

    auto tr_y_xxx_xxxxz = pbuffer.data(idx_dip_fh + 212);

    auto tr_y_xxx_xxxyy = pbuffer.data(idx_dip_fh + 213);

    auto tr_y_xxx_xxxyz = pbuffer.data(idx_dip_fh + 214);

    auto tr_y_xxx_xxxzz = pbuffer.data(idx_dip_fh + 215);

    auto tr_y_xxx_xxyyy = pbuffer.data(idx_dip_fh + 216);

    auto tr_y_xxx_xxyyz = pbuffer.data(idx_dip_fh + 217);

    auto tr_y_xxx_xxyzz = pbuffer.data(idx_dip_fh + 218);

    auto tr_y_xxx_xxzzz = pbuffer.data(idx_dip_fh + 219);

    auto tr_y_xxx_xyyyy = pbuffer.data(idx_dip_fh + 220);

    auto tr_y_xxx_xyyyz = pbuffer.data(idx_dip_fh + 221);

    auto tr_y_xxx_xyyzz = pbuffer.data(idx_dip_fh + 222);

    auto tr_y_xxx_xyzzz = pbuffer.data(idx_dip_fh + 223);

    auto tr_y_xxx_xzzzz = pbuffer.data(idx_dip_fh + 224);

    auto tr_y_xxx_yyyyy = pbuffer.data(idx_dip_fh + 225);

    auto tr_y_xxx_yyyyz = pbuffer.data(idx_dip_fh + 226);

    auto tr_y_xxx_yyyzz = pbuffer.data(idx_dip_fh + 227);

    auto tr_y_xxx_yyzzz = pbuffer.data(idx_dip_fh + 228);

    auto tr_y_xxx_yzzzz = pbuffer.data(idx_dip_fh + 229);

    auto tr_y_xxx_zzzzz = pbuffer.data(idx_dip_fh + 230);

#pragma omp simd aligned(pa_x,               \
                             tr_y_x_xxxxx,   \
                             tr_y_x_xxxxy,   \
                             tr_y_x_xxxxz,   \
                             tr_y_x_xxxyy,   \
                             tr_y_x_xxxyz,   \
                             tr_y_x_xxxzz,   \
                             tr_y_x_xxyyy,   \
                             tr_y_x_xxyyz,   \
                             tr_y_x_xxyzz,   \
                             tr_y_x_xxzzz,   \
                             tr_y_x_xyyyy,   \
                             tr_y_x_xyyyz,   \
                             tr_y_x_xyyzz,   \
                             tr_y_x_xyzzz,   \
                             tr_y_x_xzzzz,   \
                             tr_y_x_yyyyy,   \
                             tr_y_x_yyyyz,   \
                             tr_y_x_yyyzz,   \
                             tr_y_x_yyzzz,   \
                             tr_y_x_yzzzz,   \
                             tr_y_x_zzzzz,   \
                             tr_y_xx_xxxx,   \
                             tr_y_xx_xxxxx,  \
                             tr_y_xx_xxxxy,  \
                             tr_y_xx_xxxxz,  \
                             tr_y_xx_xxxy,   \
                             tr_y_xx_xxxyy,  \
                             tr_y_xx_xxxyz,  \
                             tr_y_xx_xxxz,   \
                             tr_y_xx_xxxzz,  \
                             tr_y_xx_xxyy,   \
                             tr_y_xx_xxyyy,  \
                             tr_y_xx_xxyyz,  \
                             tr_y_xx_xxyz,   \
                             tr_y_xx_xxyzz,  \
                             tr_y_xx_xxzz,   \
                             tr_y_xx_xxzzz,  \
                             tr_y_xx_xyyy,   \
                             tr_y_xx_xyyyy,  \
                             tr_y_xx_xyyyz,  \
                             tr_y_xx_xyyz,   \
                             tr_y_xx_xyyzz,  \
                             tr_y_xx_xyzz,   \
                             tr_y_xx_xyzzz,  \
                             tr_y_xx_xzzz,   \
                             tr_y_xx_xzzzz,  \
                             tr_y_xx_yyyy,   \
                             tr_y_xx_yyyyy,  \
                             tr_y_xx_yyyyz,  \
                             tr_y_xx_yyyz,   \
                             tr_y_xx_yyyzz,  \
                             tr_y_xx_yyzz,   \
                             tr_y_xx_yyzzz,  \
                             tr_y_xx_yzzz,   \
                             tr_y_xx_yzzzz,  \
                             tr_y_xx_zzzz,   \
                             tr_y_xx_zzzzz,  \
                             tr_y_xxx_xxxxx, \
                             tr_y_xxx_xxxxy, \
                             tr_y_xxx_xxxxz, \
                             tr_y_xxx_xxxyy, \
                             tr_y_xxx_xxxyz, \
                             tr_y_xxx_xxxzz, \
                             tr_y_xxx_xxyyy, \
                             tr_y_xxx_xxyyz, \
                             tr_y_xxx_xxyzz, \
                             tr_y_xxx_xxzzz, \
                             tr_y_xxx_xyyyy, \
                             tr_y_xxx_xyyyz, \
                             tr_y_xxx_xyyzz, \
                             tr_y_xxx_xyzzz, \
                             tr_y_xxx_xzzzz, \
                             tr_y_xxx_yyyyy, \
                             tr_y_xxx_yyyyz, \
                             tr_y_xxx_yyyzz, \
                             tr_y_xxx_yyzzz, \
                             tr_y_xxx_yzzzz, \
                             tr_y_xxx_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxx_xxxxx[i] = 2.0 * tr_y_x_xxxxx[i] * fe_0 + 5.0 * tr_y_xx_xxxx[i] * fe_0 + tr_y_xx_xxxxx[i] * pa_x[i];

        tr_y_xxx_xxxxy[i] = 2.0 * tr_y_x_xxxxy[i] * fe_0 + 4.0 * tr_y_xx_xxxy[i] * fe_0 + tr_y_xx_xxxxy[i] * pa_x[i];

        tr_y_xxx_xxxxz[i] = 2.0 * tr_y_x_xxxxz[i] * fe_0 + 4.0 * tr_y_xx_xxxz[i] * fe_0 + tr_y_xx_xxxxz[i] * pa_x[i];

        tr_y_xxx_xxxyy[i] = 2.0 * tr_y_x_xxxyy[i] * fe_0 + 3.0 * tr_y_xx_xxyy[i] * fe_0 + tr_y_xx_xxxyy[i] * pa_x[i];

        tr_y_xxx_xxxyz[i] = 2.0 * tr_y_x_xxxyz[i] * fe_0 + 3.0 * tr_y_xx_xxyz[i] * fe_0 + tr_y_xx_xxxyz[i] * pa_x[i];

        tr_y_xxx_xxxzz[i] = 2.0 * tr_y_x_xxxzz[i] * fe_0 + 3.0 * tr_y_xx_xxzz[i] * fe_0 + tr_y_xx_xxxzz[i] * pa_x[i];

        tr_y_xxx_xxyyy[i] = 2.0 * tr_y_x_xxyyy[i] * fe_0 + 2.0 * tr_y_xx_xyyy[i] * fe_0 + tr_y_xx_xxyyy[i] * pa_x[i];

        tr_y_xxx_xxyyz[i] = 2.0 * tr_y_x_xxyyz[i] * fe_0 + 2.0 * tr_y_xx_xyyz[i] * fe_0 + tr_y_xx_xxyyz[i] * pa_x[i];

        tr_y_xxx_xxyzz[i] = 2.0 * tr_y_x_xxyzz[i] * fe_0 + 2.0 * tr_y_xx_xyzz[i] * fe_0 + tr_y_xx_xxyzz[i] * pa_x[i];

        tr_y_xxx_xxzzz[i] = 2.0 * tr_y_x_xxzzz[i] * fe_0 + 2.0 * tr_y_xx_xzzz[i] * fe_0 + tr_y_xx_xxzzz[i] * pa_x[i];

        tr_y_xxx_xyyyy[i] = 2.0 * tr_y_x_xyyyy[i] * fe_0 + tr_y_xx_yyyy[i] * fe_0 + tr_y_xx_xyyyy[i] * pa_x[i];

        tr_y_xxx_xyyyz[i] = 2.0 * tr_y_x_xyyyz[i] * fe_0 + tr_y_xx_yyyz[i] * fe_0 + tr_y_xx_xyyyz[i] * pa_x[i];

        tr_y_xxx_xyyzz[i] = 2.0 * tr_y_x_xyyzz[i] * fe_0 + tr_y_xx_yyzz[i] * fe_0 + tr_y_xx_xyyzz[i] * pa_x[i];

        tr_y_xxx_xyzzz[i] = 2.0 * tr_y_x_xyzzz[i] * fe_0 + tr_y_xx_yzzz[i] * fe_0 + tr_y_xx_xyzzz[i] * pa_x[i];

        tr_y_xxx_xzzzz[i] = 2.0 * tr_y_x_xzzzz[i] * fe_0 + tr_y_xx_zzzz[i] * fe_0 + tr_y_xx_xzzzz[i] * pa_x[i];

        tr_y_xxx_yyyyy[i] = 2.0 * tr_y_x_yyyyy[i] * fe_0 + tr_y_xx_yyyyy[i] * pa_x[i];

        tr_y_xxx_yyyyz[i] = 2.0 * tr_y_x_yyyyz[i] * fe_0 + tr_y_xx_yyyyz[i] * pa_x[i];

        tr_y_xxx_yyyzz[i] = 2.0 * tr_y_x_yyyzz[i] * fe_0 + tr_y_xx_yyyzz[i] * pa_x[i];

        tr_y_xxx_yyzzz[i] = 2.0 * tr_y_x_yyzzz[i] * fe_0 + tr_y_xx_yyzzz[i] * pa_x[i];

        tr_y_xxx_yzzzz[i] = 2.0 * tr_y_x_yzzzz[i] * fe_0 + tr_y_xx_yzzzz[i] * pa_x[i];

        tr_y_xxx_zzzzz[i] = 2.0 * tr_y_x_zzzzz[i] * fe_0 + tr_y_xx_zzzzz[i] * pa_x[i];
    }

    // Set up 231-252 components of targeted buffer : FH

    auto tr_y_xxy_xxxxx = pbuffer.data(idx_dip_fh + 231);

    auto tr_y_xxy_xxxxy = pbuffer.data(idx_dip_fh + 232);

    auto tr_y_xxy_xxxxz = pbuffer.data(idx_dip_fh + 233);

    auto tr_y_xxy_xxxyy = pbuffer.data(idx_dip_fh + 234);

    auto tr_y_xxy_xxxyz = pbuffer.data(idx_dip_fh + 235);

    auto tr_y_xxy_xxxzz = pbuffer.data(idx_dip_fh + 236);

    auto tr_y_xxy_xxyyy = pbuffer.data(idx_dip_fh + 237);

    auto tr_y_xxy_xxyyz = pbuffer.data(idx_dip_fh + 238);

    auto tr_y_xxy_xxyzz = pbuffer.data(idx_dip_fh + 239);

    auto tr_y_xxy_xxzzz = pbuffer.data(idx_dip_fh + 240);

    auto tr_y_xxy_xyyyy = pbuffer.data(idx_dip_fh + 241);

    auto tr_y_xxy_xyyyz = pbuffer.data(idx_dip_fh + 242);

    auto tr_y_xxy_xyyzz = pbuffer.data(idx_dip_fh + 243);

    auto tr_y_xxy_xyzzz = pbuffer.data(idx_dip_fh + 244);

    auto tr_y_xxy_xzzzz = pbuffer.data(idx_dip_fh + 245);

    auto tr_y_xxy_yyyyy = pbuffer.data(idx_dip_fh + 246);

    auto tr_y_xxy_yyyyz = pbuffer.data(idx_dip_fh + 247);

    auto tr_y_xxy_yyyzz = pbuffer.data(idx_dip_fh + 248);

    auto tr_y_xxy_yyzzz = pbuffer.data(idx_dip_fh + 249);

    auto tr_y_xxy_yzzzz = pbuffer.data(idx_dip_fh + 250);

    auto tr_y_xxy_zzzzz = pbuffer.data(idx_dip_fh + 251);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_y_xx_xxxxx,  \
                             tr_y_xx_xxxxz,  \
                             tr_y_xx_xxxzz,  \
                             tr_y_xx_xxzzz,  \
                             tr_y_xx_xzzzz,  \
                             tr_y_xxy_xxxxx, \
                             tr_y_xxy_xxxxy, \
                             tr_y_xxy_xxxxz, \
                             tr_y_xxy_xxxyy, \
                             tr_y_xxy_xxxyz, \
                             tr_y_xxy_xxxzz, \
                             tr_y_xxy_xxyyy, \
                             tr_y_xxy_xxyyz, \
                             tr_y_xxy_xxyzz, \
                             tr_y_xxy_xxzzz, \
                             tr_y_xxy_xyyyy, \
                             tr_y_xxy_xyyyz, \
                             tr_y_xxy_xyyzz, \
                             tr_y_xxy_xyzzz, \
                             tr_y_xxy_xzzzz, \
                             tr_y_xxy_yyyyy, \
                             tr_y_xxy_yyyyz, \
                             tr_y_xxy_yyyzz, \
                             tr_y_xxy_yyzzz, \
                             tr_y_xxy_yzzzz, \
                             tr_y_xxy_zzzzz, \
                             tr_y_xy_xxxxy,  \
                             tr_y_xy_xxxy,   \
                             tr_y_xy_xxxyy,  \
                             tr_y_xy_xxxyz,  \
                             tr_y_xy_xxyy,   \
                             tr_y_xy_xxyyy,  \
                             tr_y_xy_xxyyz,  \
                             tr_y_xy_xxyz,   \
                             tr_y_xy_xxyzz,  \
                             tr_y_xy_xyyy,   \
                             tr_y_xy_xyyyy,  \
                             tr_y_xy_xyyyz,  \
                             tr_y_xy_xyyz,   \
                             tr_y_xy_xyyzz,  \
                             tr_y_xy_xyzz,   \
                             tr_y_xy_xyzzz,  \
                             tr_y_xy_yyyy,   \
                             tr_y_xy_yyyyy,  \
                             tr_y_xy_yyyyz,  \
                             tr_y_xy_yyyz,   \
                             tr_y_xy_yyyzz,  \
                             tr_y_xy_yyzz,   \
                             tr_y_xy_yyzzz,  \
                             tr_y_xy_yzzz,   \
                             tr_y_xy_yzzzz,  \
                             tr_y_xy_zzzzz,  \
                             tr_y_y_xxxxy,   \
                             tr_y_y_xxxyy,   \
                             tr_y_y_xxxyz,   \
                             tr_y_y_xxyyy,   \
                             tr_y_y_xxyyz,   \
                             tr_y_y_xxyzz,   \
                             tr_y_y_xyyyy,   \
                             tr_y_y_xyyyz,   \
                             tr_y_y_xyyzz,   \
                             tr_y_y_xyzzz,   \
                             tr_y_y_yyyyy,   \
                             tr_y_y_yyyyz,   \
                             tr_y_y_yyyzz,   \
                             tr_y_y_yyzzz,   \
                             tr_y_y_yzzzz,   \
                             tr_y_y_zzzzz,   \
                             ts_xx_xxxxx,    \
                             ts_xx_xxxxz,    \
                             ts_xx_xxxzz,    \
                             ts_xx_xxzzz,    \
                             ts_xx_xzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxy_xxxxx[i] = ts_xx_xxxxx[i] * fe_0 + tr_y_xx_xxxxx[i] * pa_y[i];

        tr_y_xxy_xxxxy[i] = tr_y_y_xxxxy[i] * fe_0 + 4.0 * tr_y_xy_xxxy[i] * fe_0 + tr_y_xy_xxxxy[i] * pa_x[i];

        tr_y_xxy_xxxxz[i] = ts_xx_xxxxz[i] * fe_0 + tr_y_xx_xxxxz[i] * pa_y[i];

        tr_y_xxy_xxxyy[i] = tr_y_y_xxxyy[i] * fe_0 + 3.0 * tr_y_xy_xxyy[i] * fe_0 + tr_y_xy_xxxyy[i] * pa_x[i];

        tr_y_xxy_xxxyz[i] = tr_y_y_xxxyz[i] * fe_0 + 3.0 * tr_y_xy_xxyz[i] * fe_0 + tr_y_xy_xxxyz[i] * pa_x[i];

        tr_y_xxy_xxxzz[i] = ts_xx_xxxzz[i] * fe_0 + tr_y_xx_xxxzz[i] * pa_y[i];

        tr_y_xxy_xxyyy[i] = tr_y_y_xxyyy[i] * fe_0 + 2.0 * tr_y_xy_xyyy[i] * fe_0 + tr_y_xy_xxyyy[i] * pa_x[i];

        tr_y_xxy_xxyyz[i] = tr_y_y_xxyyz[i] * fe_0 + 2.0 * tr_y_xy_xyyz[i] * fe_0 + tr_y_xy_xxyyz[i] * pa_x[i];

        tr_y_xxy_xxyzz[i] = tr_y_y_xxyzz[i] * fe_0 + 2.0 * tr_y_xy_xyzz[i] * fe_0 + tr_y_xy_xxyzz[i] * pa_x[i];

        tr_y_xxy_xxzzz[i] = ts_xx_xxzzz[i] * fe_0 + tr_y_xx_xxzzz[i] * pa_y[i];

        tr_y_xxy_xyyyy[i] = tr_y_y_xyyyy[i] * fe_0 + tr_y_xy_yyyy[i] * fe_0 + tr_y_xy_xyyyy[i] * pa_x[i];

        tr_y_xxy_xyyyz[i] = tr_y_y_xyyyz[i] * fe_0 + tr_y_xy_yyyz[i] * fe_0 + tr_y_xy_xyyyz[i] * pa_x[i];

        tr_y_xxy_xyyzz[i] = tr_y_y_xyyzz[i] * fe_0 + tr_y_xy_yyzz[i] * fe_0 + tr_y_xy_xyyzz[i] * pa_x[i];

        tr_y_xxy_xyzzz[i] = tr_y_y_xyzzz[i] * fe_0 + tr_y_xy_yzzz[i] * fe_0 + tr_y_xy_xyzzz[i] * pa_x[i];

        tr_y_xxy_xzzzz[i] = ts_xx_xzzzz[i] * fe_0 + tr_y_xx_xzzzz[i] * pa_y[i];

        tr_y_xxy_yyyyy[i] = tr_y_y_yyyyy[i] * fe_0 + tr_y_xy_yyyyy[i] * pa_x[i];

        tr_y_xxy_yyyyz[i] = tr_y_y_yyyyz[i] * fe_0 + tr_y_xy_yyyyz[i] * pa_x[i];

        tr_y_xxy_yyyzz[i] = tr_y_y_yyyzz[i] * fe_0 + tr_y_xy_yyyzz[i] * pa_x[i];

        tr_y_xxy_yyzzz[i] = tr_y_y_yyzzz[i] * fe_0 + tr_y_xy_yyzzz[i] * pa_x[i];

        tr_y_xxy_yzzzz[i] = tr_y_y_yzzzz[i] * fe_0 + tr_y_xy_yzzzz[i] * pa_x[i];

        tr_y_xxy_zzzzz[i] = tr_y_y_zzzzz[i] * fe_0 + tr_y_xy_zzzzz[i] * pa_x[i];
    }

    // Set up 252-273 components of targeted buffer : FH

    auto tr_y_xxz_xxxxx = pbuffer.data(idx_dip_fh + 252);

    auto tr_y_xxz_xxxxy = pbuffer.data(idx_dip_fh + 253);

    auto tr_y_xxz_xxxxz = pbuffer.data(idx_dip_fh + 254);

    auto tr_y_xxz_xxxyy = pbuffer.data(idx_dip_fh + 255);

    auto tr_y_xxz_xxxyz = pbuffer.data(idx_dip_fh + 256);

    auto tr_y_xxz_xxxzz = pbuffer.data(idx_dip_fh + 257);

    auto tr_y_xxz_xxyyy = pbuffer.data(idx_dip_fh + 258);

    auto tr_y_xxz_xxyyz = pbuffer.data(idx_dip_fh + 259);

    auto tr_y_xxz_xxyzz = pbuffer.data(idx_dip_fh + 260);

    auto tr_y_xxz_xxzzz = pbuffer.data(idx_dip_fh + 261);

    auto tr_y_xxz_xyyyy = pbuffer.data(idx_dip_fh + 262);

    auto tr_y_xxz_xyyyz = pbuffer.data(idx_dip_fh + 263);

    auto tr_y_xxz_xyyzz = pbuffer.data(idx_dip_fh + 264);

    auto tr_y_xxz_xyzzz = pbuffer.data(idx_dip_fh + 265);

    auto tr_y_xxz_xzzzz = pbuffer.data(idx_dip_fh + 266);

    auto tr_y_xxz_yyyyy = pbuffer.data(idx_dip_fh + 267);

    auto tr_y_xxz_yyyyz = pbuffer.data(idx_dip_fh + 268);

    auto tr_y_xxz_yyyzz = pbuffer.data(idx_dip_fh + 269);

    auto tr_y_xxz_yyzzz = pbuffer.data(idx_dip_fh + 270);

    auto tr_y_xxz_yzzzz = pbuffer.data(idx_dip_fh + 271);

    auto tr_y_xxz_zzzzz = pbuffer.data(idx_dip_fh + 272);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xx_xxxx,   \
                             tr_y_xx_xxxxx,  \
                             tr_y_xx_xxxxy,  \
                             tr_y_xx_xxxxz,  \
                             tr_y_xx_xxxy,   \
                             tr_y_xx_xxxyy,  \
                             tr_y_xx_xxxyz,  \
                             tr_y_xx_xxxz,   \
                             tr_y_xx_xxxzz,  \
                             tr_y_xx_xxyy,   \
                             tr_y_xx_xxyyy,  \
                             tr_y_xx_xxyyz,  \
                             tr_y_xx_xxyz,   \
                             tr_y_xx_xxyzz,  \
                             tr_y_xx_xxzz,   \
                             tr_y_xx_xxzzz,  \
                             tr_y_xx_xyyy,   \
                             tr_y_xx_xyyyy,  \
                             tr_y_xx_xyyyz,  \
                             tr_y_xx_xyyz,   \
                             tr_y_xx_xyyzz,  \
                             tr_y_xx_xyzz,   \
                             tr_y_xx_xyzzz,  \
                             tr_y_xx_xzzz,   \
                             tr_y_xx_xzzzz,  \
                             tr_y_xx_yyyyy,  \
                             tr_y_xxz_xxxxx, \
                             tr_y_xxz_xxxxy, \
                             tr_y_xxz_xxxxz, \
                             tr_y_xxz_xxxyy, \
                             tr_y_xxz_xxxyz, \
                             tr_y_xxz_xxxzz, \
                             tr_y_xxz_xxyyy, \
                             tr_y_xxz_xxyyz, \
                             tr_y_xxz_xxyzz, \
                             tr_y_xxz_xxzzz, \
                             tr_y_xxz_xyyyy, \
                             tr_y_xxz_xyyyz, \
                             tr_y_xxz_xyyzz, \
                             tr_y_xxz_xyzzz, \
                             tr_y_xxz_xzzzz, \
                             tr_y_xxz_yyyyy, \
                             tr_y_xxz_yyyyz, \
                             tr_y_xxz_yyyzz, \
                             tr_y_xxz_yyzzz, \
                             tr_y_xxz_yzzzz, \
                             tr_y_xxz_zzzzz, \
                             tr_y_xz_yyyyz,  \
                             tr_y_xz_yyyzz,  \
                             tr_y_xz_yyzzz,  \
                             tr_y_xz_yzzzz,  \
                             tr_y_xz_zzzzz,  \
                             tr_y_z_yyyyz,   \
                             tr_y_z_yyyzz,   \
                             tr_y_z_yyzzz,   \
                             tr_y_z_yzzzz,   \
                             tr_y_z_zzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xxz_xxxxx[i] = tr_y_xx_xxxxx[i] * pa_z[i];

        tr_y_xxz_xxxxy[i] = tr_y_xx_xxxxy[i] * pa_z[i];

        tr_y_xxz_xxxxz[i] = tr_y_xx_xxxx[i] * fe_0 + tr_y_xx_xxxxz[i] * pa_z[i];

        tr_y_xxz_xxxyy[i] = tr_y_xx_xxxyy[i] * pa_z[i];

        tr_y_xxz_xxxyz[i] = tr_y_xx_xxxy[i] * fe_0 + tr_y_xx_xxxyz[i] * pa_z[i];

        tr_y_xxz_xxxzz[i] = 2.0 * tr_y_xx_xxxz[i] * fe_0 + tr_y_xx_xxxzz[i] * pa_z[i];

        tr_y_xxz_xxyyy[i] = tr_y_xx_xxyyy[i] * pa_z[i];

        tr_y_xxz_xxyyz[i] = tr_y_xx_xxyy[i] * fe_0 + tr_y_xx_xxyyz[i] * pa_z[i];

        tr_y_xxz_xxyzz[i] = 2.0 * tr_y_xx_xxyz[i] * fe_0 + tr_y_xx_xxyzz[i] * pa_z[i];

        tr_y_xxz_xxzzz[i] = 3.0 * tr_y_xx_xxzz[i] * fe_0 + tr_y_xx_xxzzz[i] * pa_z[i];

        tr_y_xxz_xyyyy[i] = tr_y_xx_xyyyy[i] * pa_z[i];

        tr_y_xxz_xyyyz[i] = tr_y_xx_xyyy[i] * fe_0 + tr_y_xx_xyyyz[i] * pa_z[i];

        tr_y_xxz_xyyzz[i] = 2.0 * tr_y_xx_xyyz[i] * fe_0 + tr_y_xx_xyyzz[i] * pa_z[i];

        tr_y_xxz_xyzzz[i] = 3.0 * tr_y_xx_xyzz[i] * fe_0 + tr_y_xx_xyzzz[i] * pa_z[i];

        tr_y_xxz_xzzzz[i] = 4.0 * tr_y_xx_xzzz[i] * fe_0 + tr_y_xx_xzzzz[i] * pa_z[i];

        tr_y_xxz_yyyyy[i] = tr_y_xx_yyyyy[i] * pa_z[i];

        tr_y_xxz_yyyyz[i] = tr_y_z_yyyyz[i] * fe_0 + tr_y_xz_yyyyz[i] * pa_x[i];

        tr_y_xxz_yyyzz[i] = tr_y_z_yyyzz[i] * fe_0 + tr_y_xz_yyyzz[i] * pa_x[i];

        tr_y_xxz_yyzzz[i] = tr_y_z_yyzzz[i] * fe_0 + tr_y_xz_yyzzz[i] * pa_x[i];

        tr_y_xxz_yzzzz[i] = tr_y_z_yzzzz[i] * fe_0 + tr_y_xz_yzzzz[i] * pa_x[i];

        tr_y_xxz_zzzzz[i] = tr_y_z_zzzzz[i] * fe_0 + tr_y_xz_zzzzz[i] * pa_x[i];
    }

    // Set up 273-294 components of targeted buffer : FH

    auto tr_y_xyy_xxxxx = pbuffer.data(idx_dip_fh + 273);

    auto tr_y_xyy_xxxxy = pbuffer.data(idx_dip_fh + 274);

    auto tr_y_xyy_xxxxz = pbuffer.data(idx_dip_fh + 275);

    auto tr_y_xyy_xxxyy = pbuffer.data(idx_dip_fh + 276);

    auto tr_y_xyy_xxxyz = pbuffer.data(idx_dip_fh + 277);

    auto tr_y_xyy_xxxzz = pbuffer.data(idx_dip_fh + 278);

    auto tr_y_xyy_xxyyy = pbuffer.data(idx_dip_fh + 279);

    auto tr_y_xyy_xxyyz = pbuffer.data(idx_dip_fh + 280);

    auto tr_y_xyy_xxyzz = pbuffer.data(idx_dip_fh + 281);

    auto tr_y_xyy_xxzzz = pbuffer.data(idx_dip_fh + 282);

    auto tr_y_xyy_xyyyy = pbuffer.data(idx_dip_fh + 283);

    auto tr_y_xyy_xyyyz = pbuffer.data(idx_dip_fh + 284);

    auto tr_y_xyy_xyyzz = pbuffer.data(idx_dip_fh + 285);

    auto tr_y_xyy_xyzzz = pbuffer.data(idx_dip_fh + 286);

    auto tr_y_xyy_xzzzz = pbuffer.data(idx_dip_fh + 287);

    auto tr_y_xyy_yyyyy = pbuffer.data(idx_dip_fh + 288);

    auto tr_y_xyy_yyyyz = pbuffer.data(idx_dip_fh + 289);

    auto tr_y_xyy_yyyzz = pbuffer.data(idx_dip_fh + 290);

    auto tr_y_xyy_yyzzz = pbuffer.data(idx_dip_fh + 291);

    auto tr_y_xyy_yzzzz = pbuffer.data(idx_dip_fh + 292);

    auto tr_y_xyy_zzzzz = pbuffer.data(idx_dip_fh + 293);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xyy_xxxxx, \
                             tr_y_xyy_xxxxy, \
                             tr_y_xyy_xxxxz, \
                             tr_y_xyy_xxxyy, \
                             tr_y_xyy_xxxyz, \
                             tr_y_xyy_xxxzz, \
                             tr_y_xyy_xxyyy, \
                             tr_y_xyy_xxyyz, \
                             tr_y_xyy_xxyzz, \
                             tr_y_xyy_xxzzz, \
                             tr_y_xyy_xyyyy, \
                             tr_y_xyy_xyyyz, \
                             tr_y_xyy_xyyzz, \
                             tr_y_xyy_xyzzz, \
                             tr_y_xyy_xzzzz, \
                             tr_y_xyy_yyyyy, \
                             tr_y_xyy_yyyyz, \
                             tr_y_xyy_yyyzz, \
                             tr_y_xyy_yyzzz, \
                             tr_y_xyy_yzzzz, \
                             tr_y_xyy_zzzzz, \
                             tr_y_yy_xxxx,   \
                             tr_y_yy_xxxxx,  \
                             tr_y_yy_xxxxy,  \
                             tr_y_yy_xxxxz,  \
                             tr_y_yy_xxxy,   \
                             tr_y_yy_xxxyy,  \
                             tr_y_yy_xxxyz,  \
                             tr_y_yy_xxxz,   \
                             tr_y_yy_xxxzz,  \
                             tr_y_yy_xxyy,   \
                             tr_y_yy_xxyyy,  \
                             tr_y_yy_xxyyz,  \
                             tr_y_yy_xxyz,   \
                             tr_y_yy_xxyzz,  \
                             tr_y_yy_xxzz,   \
                             tr_y_yy_xxzzz,  \
                             tr_y_yy_xyyy,   \
                             tr_y_yy_xyyyy,  \
                             tr_y_yy_xyyyz,  \
                             tr_y_yy_xyyz,   \
                             tr_y_yy_xyyzz,  \
                             tr_y_yy_xyzz,   \
                             tr_y_yy_xyzzz,  \
                             tr_y_yy_xzzz,   \
                             tr_y_yy_xzzzz,  \
                             tr_y_yy_yyyy,   \
                             tr_y_yy_yyyyy,  \
                             tr_y_yy_yyyyz,  \
                             tr_y_yy_yyyz,   \
                             tr_y_yy_yyyzz,  \
                             tr_y_yy_yyzz,   \
                             tr_y_yy_yyzzz,  \
                             tr_y_yy_yzzz,   \
                             tr_y_yy_yzzzz,  \
                             tr_y_yy_zzzz,   \
                             tr_y_yy_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyy_xxxxx[i] = 5.0 * tr_y_yy_xxxx[i] * fe_0 + tr_y_yy_xxxxx[i] * pa_x[i];

        tr_y_xyy_xxxxy[i] = 4.0 * tr_y_yy_xxxy[i] * fe_0 + tr_y_yy_xxxxy[i] * pa_x[i];

        tr_y_xyy_xxxxz[i] = 4.0 * tr_y_yy_xxxz[i] * fe_0 + tr_y_yy_xxxxz[i] * pa_x[i];

        tr_y_xyy_xxxyy[i] = 3.0 * tr_y_yy_xxyy[i] * fe_0 + tr_y_yy_xxxyy[i] * pa_x[i];

        tr_y_xyy_xxxyz[i] = 3.0 * tr_y_yy_xxyz[i] * fe_0 + tr_y_yy_xxxyz[i] * pa_x[i];

        tr_y_xyy_xxxzz[i] = 3.0 * tr_y_yy_xxzz[i] * fe_0 + tr_y_yy_xxxzz[i] * pa_x[i];

        tr_y_xyy_xxyyy[i] = 2.0 * tr_y_yy_xyyy[i] * fe_0 + tr_y_yy_xxyyy[i] * pa_x[i];

        tr_y_xyy_xxyyz[i] = 2.0 * tr_y_yy_xyyz[i] * fe_0 + tr_y_yy_xxyyz[i] * pa_x[i];

        tr_y_xyy_xxyzz[i] = 2.0 * tr_y_yy_xyzz[i] * fe_0 + tr_y_yy_xxyzz[i] * pa_x[i];

        tr_y_xyy_xxzzz[i] = 2.0 * tr_y_yy_xzzz[i] * fe_0 + tr_y_yy_xxzzz[i] * pa_x[i];

        tr_y_xyy_xyyyy[i] = tr_y_yy_yyyy[i] * fe_0 + tr_y_yy_xyyyy[i] * pa_x[i];

        tr_y_xyy_xyyyz[i] = tr_y_yy_yyyz[i] * fe_0 + tr_y_yy_xyyyz[i] * pa_x[i];

        tr_y_xyy_xyyzz[i] = tr_y_yy_yyzz[i] * fe_0 + tr_y_yy_xyyzz[i] * pa_x[i];

        tr_y_xyy_xyzzz[i] = tr_y_yy_yzzz[i] * fe_0 + tr_y_yy_xyzzz[i] * pa_x[i];

        tr_y_xyy_xzzzz[i] = tr_y_yy_zzzz[i] * fe_0 + tr_y_yy_xzzzz[i] * pa_x[i];

        tr_y_xyy_yyyyy[i] = tr_y_yy_yyyyy[i] * pa_x[i];

        tr_y_xyy_yyyyz[i] = tr_y_yy_yyyyz[i] * pa_x[i];

        tr_y_xyy_yyyzz[i] = tr_y_yy_yyyzz[i] * pa_x[i];

        tr_y_xyy_yyzzz[i] = tr_y_yy_yyzzz[i] * pa_x[i];

        tr_y_xyy_yzzzz[i] = tr_y_yy_yzzzz[i] * pa_x[i];

        tr_y_xyy_zzzzz[i] = tr_y_yy_zzzzz[i] * pa_x[i];
    }

    // Set up 294-315 components of targeted buffer : FH

    auto tr_y_xyz_xxxxx = pbuffer.data(idx_dip_fh + 294);

    auto tr_y_xyz_xxxxy = pbuffer.data(idx_dip_fh + 295);

    auto tr_y_xyz_xxxxz = pbuffer.data(idx_dip_fh + 296);

    auto tr_y_xyz_xxxyy = pbuffer.data(idx_dip_fh + 297);

    auto tr_y_xyz_xxxyz = pbuffer.data(idx_dip_fh + 298);

    auto tr_y_xyz_xxxzz = pbuffer.data(idx_dip_fh + 299);

    auto tr_y_xyz_xxyyy = pbuffer.data(idx_dip_fh + 300);

    auto tr_y_xyz_xxyyz = pbuffer.data(idx_dip_fh + 301);

    auto tr_y_xyz_xxyzz = pbuffer.data(idx_dip_fh + 302);

    auto tr_y_xyz_xxzzz = pbuffer.data(idx_dip_fh + 303);

    auto tr_y_xyz_xyyyy = pbuffer.data(idx_dip_fh + 304);

    auto tr_y_xyz_xyyyz = pbuffer.data(idx_dip_fh + 305);

    auto tr_y_xyz_xyyzz = pbuffer.data(idx_dip_fh + 306);

    auto tr_y_xyz_xyzzz = pbuffer.data(idx_dip_fh + 307);

    auto tr_y_xyz_xzzzz = pbuffer.data(idx_dip_fh + 308);

    auto tr_y_xyz_yyyyy = pbuffer.data(idx_dip_fh + 309);

    auto tr_y_xyz_yyyyz = pbuffer.data(idx_dip_fh + 310);

    auto tr_y_xyz_yyyzz = pbuffer.data(idx_dip_fh + 311);

    auto tr_y_xyz_yyzzz = pbuffer.data(idx_dip_fh + 312);

    auto tr_y_xyz_yzzzz = pbuffer.data(idx_dip_fh + 313);

    auto tr_y_xyz_zzzzz = pbuffer.data(idx_dip_fh + 314);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_xy_xxxxx,  \
                             tr_y_xy_xxxxy,  \
                             tr_y_xy_xxxyy,  \
                             tr_y_xy_xxyyy,  \
                             tr_y_xy_xyyyy,  \
                             tr_y_xyz_xxxxx, \
                             tr_y_xyz_xxxxy, \
                             tr_y_xyz_xxxxz, \
                             tr_y_xyz_xxxyy, \
                             tr_y_xyz_xxxyz, \
                             tr_y_xyz_xxxzz, \
                             tr_y_xyz_xxyyy, \
                             tr_y_xyz_xxyyz, \
                             tr_y_xyz_xxyzz, \
                             tr_y_xyz_xxzzz, \
                             tr_y_xyz_xyyyy, \
                             tr_y_xyz_xyyyz, \
                             tr_y_xyz_xyyzz, \
                             tr_y_xyz_xyzzz, \
                             tr_y_xyz_xzzzz, \
                             tr_y_xyz_yyyyy, \
                             tr_y_xyz_yyyyz, \
                             tr_y_xyz_yyyzz, \
                             tr_y_xyz_yyzzz, \
                             tr_y_xyz_yzzzz, \
                             tr_y_xyz_zzzzz, \
                             tr_y_yz_xxxxz,  \
                             tr_y_yz_xxxyz,  \
                             tr_y_yz_xxxz,   \
                             tr_y_yz_xxxzz,  \
                             tr_y_yz_xxyyz,  \
                             tr_y_yz_xxyz,   \
                             tr_y_yz_xxyzz,  \
                             tr_y_yz_xxzz,   \
                             tr_y_yz_xxzzz,  \
                             tr_y_yz_xyyyz,  \
                             tr_y_yz_xyyz,   \
                             tr_y_yz_xyyzz,  \
                             tr_y_yz_xyzz,   \
                             tr_y_yz_xyzzz,  \
                             tr_y_yz_xzzz,   \
                             tr_y_yz_xzzzz,  \
                             tr_y_yz_yyyyy,  \
                             tr_y_yz_yyyyz,  \
                             tr_y_yz_yyyz,   \
                             tr_y_yz_yyyzz,  \
                             tr_y_yz_yyzz,   \
                             tr_y_yz_yyzzz,  \
                             tr_y_yz_yzzz,   \
                             tr_y_yz_yzzzz,  \
                             tr_y_yz_zzzz,   \
                             tr_y_yz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xyz_xxxxx[i] = tr_y_xy_xxxxx[i] * pa_z[i];

        tr_y_xyz_xxxxy[i] = tr_y_xy_xxxxy[i] * pa_z[i];

        tr_y_xyz_xxxxz[i] = 4.0 * tr_y_yz_xxxz[i] * fe_0 + tr_y_yz_xxxxz[i] * pa_x[i];

        tr_y_xyz_xxxyy[i] = tr_y_xy_xxxyy[i] * pa_z[i];

        tr_y_xyz_xxxyz[i] = 3.0 * tr_y_yz_xxyz[i] * fe_0 + tr_y_yz_xxxyz[i] * pa_x[i];

        tr_y_xyz_xxxzz[i] = 3.0 * tr_y_yz_xxzz[i] * fe_0 + tr_y_yz_xxxzz[i] * pa_x[i];

        tr_y_xyz_xxyyy[i] = tr_y_xy_xxyyy[i] * pa_z[i];

        tr_y_xyz_xxyyz[i] = 2.0 * tr_y_yz_xyyz[i] * fe_0 + tr_y_yz_xxyyz[i] * pa_x[i];

        tr_y_xyz_xxyzz[i] = 2.0 * tr_y_yz_xyzz[i] * fe_0 + tr_y_yz_xxyzz[i] * pa_x[i];

        tr_y_xyz_xxzzz[i] = 2.0 * tr_y_yz_xzzz[i] * fe_0 + tr_y_yz_xxzzz[i] * pa_x[i];

        tr_y_xyz_xyyyy[i] = tr_y_xy_xyyyy[i] * pa_z[i];

        tr_y_xyz_xyyyz[i] = tr_y_yz_yyyz[i] * fe_0 + tr_y_yz_xyyyz[i] * pa_x[i];

        tr_y_xyz_xyyzz[i] = tr_y_yz_yyzz[i] * fe_0 + tr_y_yz_xyyzz[i] * pa_x[i];

        tr_y_xyz_xyzzz[i] = tr_y_yz_yzzz[i] * fe_0 + tr_y_yz_xyzzz[i] * pa_x[i];

        tr_y_xyz_xzzzz[i] = tr_y_yz_zzzz[i] * fe_0 + tr_y_yz_xzzzz[i] * pa_x[i];

        tr_y_xyz_yyyyy[i] = tr_y_yz_yyyyy[i] * pa_x[i];

        tr_y_xyz_yyyyz[i] = tr_y_yz_yyyyz[i] * pa_x[i];

        tr_y_xyz_yyyzz[i] = tr_y_yz_yyyzz[i] * pa_x[i];

        tr_y_xyz_yyzzz[i] = tr_y_yz_yyzzz[i] * pa_x[i];

        tr_y_xyz_yzzzz[i] = tr_y_yz_yzzzz[i] * pa_x[i];

        tr_y_xyz_zzzzz[i] = tr_y_yz_zzzzz[i] * pa_x[i];
    }

    // Set up 315-336 components of targeted buffer : FH

    auto tr_y_xzz_xxxxx = pbuffer.data(idx_dip_fh + 315);

    auto tr_y_xzz_xxxxy = pbuffer.data(idx_dip_fh + 316);

    auto tr_y_xzz_xxxxz = pbuffer.data(idx_dip_fh + 317);

    auto tr_y_xzz_xxxyy = pbuffer.data(idx_dip_fh + 318);

    auto tr_y_xzz_xxxyz = pbuffer.data(idx_dip_fh + 319);

    auto tr_y_xzz_xxxzz = pbuffer.data(idx_dip_fh + 320);

    auto tr_y_xzz_xxyyy = pbuffer.data(idx_dip_fh + 321);

    auto tr_y_xzz_xxyyz = pbuffer.data(idx_dip_fh + 322);

    auto tr_y_xzz_xxyzz = pbuffer.data(idx_dip_fh + 323);

    auto tr_y_xzz_xxzzz = pbuffer.data(idx_dip_fh + 324);

    auto tr_y_xzz_xyyyy = pbuffer.data(idx_dip_fh + 325);

    auto tr_y_xzz_xyyyz = pbuffer.data(idx_dip_fh + 326);

    auto tr_y_xzz_xyyzz = pbuffer.data(idx_dip_fh + 327);

    auto tr_y_xzz_xyzzz = pbuffer.data(idx_dip_fh + 328);

    auto tr_y_xzz_xzzzz = pbuffer.data(idx_dip_fh + 329);

    auto tr_y_xzz_yyyyy = pbuffer.data(idx_dip_fh + 330);

    auto tr_y_xzz_yyyyz = pbuffer.data(idx_dip_fh + 331);

    auto tr_y_xzz_yyyzz = pbuffer.data(idx_dip_fh + 332);

    auto tr_y_xzz_yyzzz = pbuffer.data(idx_dip_fh + 333);

    auto tr_y_xzz_yzzzz = pbuffer.data(idx_dip_fh + 334);

    auto tr_y_xzz_zzzzz = pbuffer.data(idx_dip_fh + 335);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xzz_xxxxx, \
                             tr_y_xzz_xxxxy, \
                             tr_y_xzz_xxxxz, \
                             tr_y_xzz_xxxyy, \
                             tr_y_xzz_xxxyz, \
                             tr_y_xzz_xxxzz, \
                             tr_y_xzz_xxyyy, \
                             tr_y_xzz_xxyyz, \
                             tr_y_xzz_xxyzz, \
                             tr_y_xzz_xxzzz, \
                             tr_y_xzz_xyyyy, \
                             tr_y_xzz_xyyyz, \
                             tr_y_xzz_xyyzz, \
                             tr_y_xzz_xyzzz, \
                             tr_y_xzz_xzzzz, \
                             tr_y_xzz_yyyyy, \
                             tr_y_xzz_yyyyz, \
                             tr_y_xzz_yyyzz, \
                             tr_y_xzz_yyzzz, \
                             tr_y_xzz_yzzzz, \
                             tr_y_xzz_zzzzz, \
                             tr_y_zz_xxxx,   \
                             tr_y_zz_xxxxx,  \
                             tr_y_zz_xxxxy,  \
                             tr_y_zz_xxxxz,  \
                             tr_y_zz_xxxy,   \
                             tr_y_zz_xxxyy,  \
                             tr_y_zz_xxxyz,  \
                             tr_y_zz_xxxz,   \
                             tr_y_zz_xxxzz,  \
                             tr_y_zz_xxyy,   \
                             tr_y_zz_xxyyy,  \
                             tr_y_zz_xxyyz,  \
                             tr_y_zz_xxyz,   \
                             tr_y_zz_xxyzz,  \
                             tr_y_zz_xxzz,   \
                             tr_y_zz_xxzzz,  \
                             tr_y_zz_xyyy,   \
                             tr_y_zz_xyyyy,  \
                             tr_y_zz_xyyyz,  \
                             tr_y_zz_xyyz,   \
                             tr_y_zz_xyyzz,  \
                             tr_y_zz_xyzz,   \
                             tr_y_zz_xyzzz,  \
                             tr_y_zz_xzzz,   \
                             tr_y_zz_xzzzz,  \
                             tr_y_zz_yyyy,   \
                             tr_y_zz_yyyyy,  \
                             tr_y_zz_yyyyz,  \
                             tr_y_zz_yyyz,   \
                             tr_y_zz_yyyzz,  \
                             tr_y_zz_yyzz,   \
                             tr_y_zz_yyzzz,  \
                             tr_y_zz_yzzz,   \
                             tr_y_zz_yzzzz,  \
                             tr_y_zz_zzzz,   \
                             tr_y_zz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xzz_xxxxx[i] = 5.0 * tr_y_zz_xxxx[i] * fe_0 + tr_y_zz_xxxxx[i] * pa_x[i];

        tr_y_xzz_xxxxy[i] = 4.0 * tr_y_zz_xxxy[i] * fe_0 + tr_y_zz_xxxxy[i] * pa_x[i];

        tr_y_xzz_xxxxz[i] = 4.0 * tr_y_zz_xxxz[i] * fe_0 + tr_y_zz_xxxxz[i] * pa_x[i];

        tr_y_xzz_xxxyy[i] = 3.0 * tr_y_zz_xxyy[i] * fe_0 + tr_y_zz_xxxyy[i] * pa_x[i];

        tr_y_xzz_xxxyz[i] = 3.0 * tr_y_zz_xxyz[i] * fe_0 + tr_y_zz_xxxyz[i] * pa_x[i];

        tr_y_xzz_xxxzz[i] = 3.0 * tr_y_zz_xxzz[i] * fe_0 + tr_y_zz_xxxzz[i] * pa_x[i];

        tr_y_xzz_xxyyy[i] = 2.0 * tr_y_zz_xyyy[i] * fe_0 + tr_y_zz_xxyyy[i] * pa_x[i];

        tr_y_xzz_xxyyz[i] = 2.0 * tr_y_zz_xyyz[i] * fe_0 + tr_y_zz_xxyyz[i] * pa_x[i];

        tr_y_xzz_xxyzz[i] = 2.0 * tr_y_zz_xyzz[i] * fe_0 + tr_y_zz_xxyzz[i] * pa_x[i];

        tr_y_xzz_xxzzz[i] = 2.0 * tr_y_zz_xzzz[i] * fe_0 + tr_y_zz_xxzzz[i] * pa_x[i];

        tr_y_xzz_xyyyy[i] = tr_y_zz_yyyy[i] * fe_0 + tr_y_zz_xyyyy[i] * pa_x[i];

        tr_y_xzz_xyyyz[i] = tr_y_zz_yyyz[i] * fe_0 + tr_y_zz_xyyyz[i] * pa_x[i];

        tr_y_xzz_xyyzz[i] = tr_y_zz_yyzz[i] * fe_0 + tr_y_zz_xyyzz[i] * pa_x[i];

        tr_y_xzz_xyzzz[i] = tr_y_zz_yzzz[i] * fe_0 + tr_y_zz_xyzzz[i] * pa_x[i];

        tr_y_xzz_xzzzz[i] = tr_y_zz_zzzz[i] * fe_0 + tr_y_zz_xzzzz[i] * pa_x[i];

        tr_y_xzz_yyyyy[i] = tr_y_zz_yyyyy[i] * pa_x[i];

        tr_y_xzz_yyyyz[i] = tr_y_zz_yyyyz[i] * pa_x[i];

        tr_y_xzz_yyyzz[i] = tr_y_zz_yyyzz[i] * pa_x[i];

        tr_y_xzz_yyzzz[i] = tr_y_zz_yyzzz[i] * pa_x[i];

        tr_y_xzz_yzzzz[i] = tr_y_zz_yzzzz[i] * pa_x[i];

        tr_y_xzz_zzzzz[i] = tr_y_zz_zzzzz[i] * pa_x[i];
    }

    // Set up 336-357 components of targeted buffer : FH

    auto tr_y_yyy_xxxxx = pbuffer.data(idx_dip_fh + 336);

    auto tr_y_yyy_xxxxy = pbuffer.data(idx_dip_fh + 337);

    auto tr_y_yyy_xxxxz = pbuffer.data(idx_dip_fh + 338);

    auto tr_y_yyy_xxxyy = pbuffer.data(idx_dip_fh + 339);

    auto tr_y_yyy_xxxyz = pbuffer.data(idx_dip_fh + 340);

    auto tr_y_yyy_xxxzz = pbuffer.data(idx_dip_fh + 341);

    auto tr_y_yyy_xxyyy = pbuffer.data(idx_dip_fh + 342);

    auto tr_y_yyy_xxyyz = pbuffer.data(idx_dip_fh + 343);

    auto tr_y_yyy_xxyzz = pbuffer.data(idx_dip_fh + 344);

    auto tr_y_yyy_xxzzz = pbuffer.data(idx_dip_fh + 345);

    auto tr_y_yyy_xyyyy = pbuffer.data(idx_dip_fh + 346);

    auto tr_y_yyy_xyyyz = pbuffer.data(idx_dip_fh + 347);

    auto tr_y_yyy_xyyzz = pbuffer.data(idx_dip_fh + 348);

    auto tr_y_yyy_xyzzz = pbuffer.data(idx_dip_fh + 349);

    auto tr_y_yyy_xzzzz = pbuffer.data(idx_dip_fh + 350);

    auto tr_y_yyy_yyyyy = pbuffer.data(idx_dip_fh + 351);

    auto tr_y_yyy_yyyyz = pbuffer.data(idx_dip_fh + 352);

    auto tr_y_yyy_yyyzz = pbuffer.data(idx_dip_fh + 353);

    auto tr_y_yyy_yyzzz = pbuffer.data(idx_dip_fh + 354);

    auto tr_y_yyy_yzzzz = pbuffer.data(idx_dip_fh + 355);

    auto tr_y_yyy_zzzzz = pbuffer.data(idx_dip_fh + 356);

#pragma omp simd aligned(pa_y,               \
                             tr_y_y_xxxxx,   \
                             tr_y_y_xxxxy,   \
                             tr_y_y_xxxxz,   \
                             tr_y_y_xxxyy,   \
                             tr_y_y_xxxyz,   \
                             tr_y_y_xxxzz,   \
                             tr_y_y_xxyyy,   \
                             tr_y_y_xxyyz,   \
                             tr_y_y_xxyzz,   \
                             tr_y_y_xxzzz,   \
                             tr_y_y_xyyyy,   \
                             tr_y_y_xyyyz,   \
                             tr_y_y_xyyzz,   \
                             tr_y_y_xyzzz,   \
                             tr_y_y_xzzzz,   \
                             tr_y_y_yyyyy,   \
                             tr_y_y_yyyyz,   \
                             tr_y_y_yyyzz,   \
                             tr_y_y_yyzzz,   \
                             tr_y_y_yzzzz,   \
                             tr_y_y_zzzzz,   \
                             tr_y_yy_xxxx,   \
                             tr_y_yy_xxxxx,  \
                             tr_y_yy_xxxxy,  \
                             tr_y_yy_xxxxz,  \
                             tr_y_yy_xxxy,   \
                             tr_y_yy_xxxyy,  \
                             tr_y_yy_xxxyz,  \
                             tr_y_yy_xxxz,   \
                             tr_y_yy_xxxzz,  \
                             tr_y_yy_xxyy,   \
                             tr_y_yy_xxyyy,  \
                             tr_y_yy_xxyyz,  \
                             tr_y_yy_xxyz,   \
                             tr_y_yy_xxyzz,  \
                             tr_y_yy_xxzz,   \
                             tr_y_yy_xxzzz,  \
                             tr_y_yy_xyyy,   \
                             tr_y_yy_xyyyy,  \
                             tr_y_yy_xyyyz,  \
                             tr_y_yy_xyyz,   \
                             tr_y_yy_xyyzz,  \
                             tr_y_yy_xyzz,   \
                             tr_y_yy_xyzzz,  \
                             tr_y_yy_xzzz,   \
                             tr_y_yy_xzzzz,  \
                             tr_y_yy_yyyy,   \
                             tr_y_yy_yyyyy,  \
                             tr_y_yy_yyyyz,  \
                             tr_y_yy_yyyz,   \
                             tr_y_yy_yyyzz,  \
                             tr_y_yy_yyzz,   \
                             tr_y_yy_yyzzz,  \
                             tr_y_yy_yzzz,   \
                             tr_y_yy_yzzzz,  \
                             tr_y_yy_zzzz,   \
                             tr_y_yy_zzzzz,  \
                             tr_y_yyy_xxxxx, \
                             tr_y_yyy_xxxxy, \
                             tr_y_yyy_xxxxz, \
                             tr_y_yyy_xxxyy, \
                             tr_y_yyy_xxxyz, \
                             tr_y_yyy_xxxzz, \
                             tr_y_yyy_xxyyy, \
                             tr_y_yyy_xxyyz, \
                             tr_y_yyy_xxyzz, \
                             tr_y_yyy_xxzzz, \
                             tr_y_yyy_xyyyy, \
                             tr_y_yyy_xyyyz, \
                             tr_y_yyy_xyyzz, \
                             tr_y_yyy_xyzzz, \
                             tr_y_yyy_xzzzz, \
                             tr_y_yyy_yyyyy, \
                             tr_y_yyy_yyyyz, \
                             tr_y_yyy_yyyzz, \
                             tr_y_yyy_yyzzz, \
                             tr_y_yyy_yzzzz, \
                             tr_y_yyy_zzzzz, \
                             ts_yy_xxxxx,    \
                             ts_yy_xxxxy,    \
                             ts_yy_xxxxz,    \
                             ts_yy_xxxyy,    \
                             ts_yy_xxxyz,    \
                             ts_yy_xxxzz,    \
                             ts_yy_xxyyy,    \
                             ts_yy_xxyyz,    \
                             ts_yy_xxyzz,    \
                             ts_yy_xxzzz,    \
                             ts_yy_xyyyy,    \
                             ts_yy_xyyyz,    \
                             ts_yy_xyyzz,    \
                             ts_yy_xyzzz,    \
                             ts_yy_xzzzz,    \
                             ts_yy_yyyyy,    \
                             ts_yy_yyyyz,    \
                             ts_yy_yyyzz,    \
                             ts_yy_yyzzz,    \
                             ts_yy_yzzzz,    \
                             ts_yy_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyy_xxxxx[i] = 2.0 * tr_y_y_xxxxx[i] * fe_0 + ts_yy_xxxxx[i] * fe_0 + tr_y_yy_xxxxx[i] * pa_y[i];

        tr_y_yyy_xxxxy[i] = 2.0 * tr_y_y_xxxxy[i] * fe_0 + tr_y_yy_xxxx[i] * fe_0 + ts_yy_xxxxy[i] * fe_0 + tr_y_yy_xxxxy[i] * pa_y[i];

        tr_y_yyy_xxxxz[i] = 2.0 * tr_y_y_xxxxz[i] * fe_0 + ts_yy_xxxxz[i] * fe_0 + tr_y_yy_xxxxz[i] * pa_y[i];

        tr_y_yyy_xxxyy[i] = 2.0 * tr_y_y_xxxyy[i] * fe_0 + 2.0 * tr_y_yy_xxxy[i] * fe_0 + ts_yy_xxxyy[i] * fe_0 + tr_y_yy_xxxyy[i] * pa_y[i];

        tr_y_yyy_xxxyz[i] = 2.0 * tr_y_y_xxxyz[i] * fe_0 + tr_y_yy_xxxz[i] * fe_0 + ts_yy_xxxyz[i] * fe_0 + tr_y_yy_xxxyz[i] * pa_y[i];

        tr_y_yyy_xxxzz[i] = 2.0 * tr_y_y_xxxzz[i] * fe_0 + ts_yy_xxxzz[i] * fe_0 + tr_y_yy_xxxzz[i] * pa_y[i];

        tr_y_yyy_xxyyy[i] = 2.0 * tr_y_y_xxyyy[i] * fe_0 + 3.0 * tr_y_yy_xxyy[i] * fe_0 + ts_yy_xxyyy[i] * fe_0 + tr_y_yy_xxyyy[i] * pa_y[i];

        tr_y_yyy_xxyyz[i] = 2.0 * tr_y_y_xxyyz[i] * fe_0 + 2.0 * tr_y_yy_xxyz[i] * fe_0 + ts_yy_xxyyz[i] * fe_0 + tr_y_yy_xxyyz[i] * pa_y[i];

        tr_y_yyy_xxyzz[i] = 2.0 * tr_y_y_xxyzz[i] * fe_0 + tr_y_yy_xxzz[i] * fe_0 + ts_yy_xxyzz[i] * fe_0 + tr_y_yy_xxyzz[i] * pa_y[i];

        tr_y_yyy_xxzzz[i] = 2.0 * tr_y_y_xxzzz[i] * fe_0 + ts_yy_xxzzz[i] * fe_0 + tr_y_yy_xxzzz[i] * pa_y[i];

        tr_y_yyy_xyyyy[i] = 2.0 * tr_y_y_xyyyy[i] * fe_0 + 4.0 * tr_y_yy_xyyy[i] * fe_0 + ts_yy_xyyyy[i] * fe_0 + tr_y_yy_xyyyy[i] * pa_y[i];

        tr_y_yyy_xyyyz[i] = 2.0 * tr_y_y_xyyyz[i] * fe_0 + 3.0 * tr_y_yy_xyyz[i] * fe_0 + ts_yy_xyyyz[i] * fe_0 + tr_y_yy_xyyyz[i] * pa_y[i];

        tr_y_yyy_xyyzz[i] = 2.0 * tr_y_y_xyyzz[i] * fe_0 + 2.0 * tr_y_yy_xyzz[i] * fe_0 + ts_yy_xyyzz[i] * fe_0 + tr_y_yy_xyyzz[i] * pa_y[i];

        tr_y_yyy_xyzzz[i] = 2.0 * tr_y_y_xyzzz[i] * fe_0 + tr_y_yy_xzzz[i] * fe_0 + ts_yy_xyzzz[i] * fe_0 + tr_y_yy_xyzzz[i] * pa_y[i];

        tr_y_yyy_xzzzz[i] = 2.0 * tr_y_y_xzzzz[i] * fe_0 + ts_yy_xzzzz[i] * fe_0 + tr_y_yy_xzzzz[i] * pa_y[i];

        tr_y_yyy_yyyyy[i] = 2.0 * tr_y_y_yyyyy[i] * fe_0 + 5.0 * tr_y_yy_yyyy[i] * fe_0 + ts_yy_yyyyy[i] * fe_0 + tr_y_yy_yyyyy[i] * pa_y[i];

        tr_y_yyy_yyyyz[i] = 2.0 * tr_y_y_yyyyz[i] * fe_0 + 4.0 * tr_y_yy_yyyz[i] * fe_0 + ts_yy_yyyyz[i] * fe_0 + tr_y_yy_yyyyz[i] * pa_y[i];

        tr_y_yyy_yyyzz[i] = 2.0 * tr_y_y_yyyzz[i] * fe_0 + 3.0 * tr_y_yy_yyzz[i] * fe_0 + ts_yy_yyyzz[i] * fe_0 + tr_y_yy_yyyzz[i] * pa_y[i];

        tr_y_yyy_yyzzz[i] = 2.0 * tr_y_y_yyzzz[i] * fe_0 + 2.0 * tr_y_yy_yzzz[i] * fe_0 + ts_yy_yyzzz[i] * fe_0 + tr_y_yy_yyzzz[i] * pa_y[i];

        tr_y_yyy_yzzzz[i] = 2.0 * tr_y_y_yzzzz[i] * fe_0 + tr_y_yy_zzzz[i] * fe_0 + ts_yy_yzzzz[i] * fe_0 + tr_y_yy_yzzzz[i] * pa_y[i];

        tr_y_yyy_zzzzz[i] = 2.0 * tr_y_y_zzzzz[i] * fe_0 + ts_yy_zzzzz[i] * fe_0 + tr_y_yy_zzzzz[i] * pa_y[i];
    }

    // Set up 357-378 components of targeted buffer : FH

    auto tr_y_yyz_xxxxx = pbuffer.data(idx_dip_fh + 357);

    auto tr_y_yyz_xxxxy = pbuffer.data(idx_dip_fh + 358);

    auto tr_y_yyz_xxxxz = pbuffer.data(idx_dip_fh + 359);

    auto tr_y_yyz_xxxyy = pbuffer.data(idx_dip_fh + 360);

    auto tr_y_yyz_xxxyz = pbuffer.data(idx_dip_fh + 361);

    auto tr_y_yyz_xxxzz = pbuffer.data(idx_dip_fh + 362);

    auto tr_y_yyz_xxyyy = pbuffer.data(idx_dip_fh + 363);

    auto tr_y_yyz_xxyyz = pbuffer.data(idx_dip_fh + 364);

    auto tr_y_yyz_xxyzz = pbuffer.data(idx_dip_fh + 365);

    auto tr_y_yyz_xxzzz = pbuffer.data(idx_dip_fh + 366);

    auto tr_y_yyz_xyyyy = pbuffer.data(idx_dip_fh + 367);

    auto tr_y_yyz_xyyyz = pbuffer.data(idx_dip_fh + 368);

    auto tr_y_yyz_xyyzz = pbuffer.data(idx_dip_fh + 369);

    auto tr_y_yyz_xyzzz = pbuffer.data(idx_dip_fh + 370);

    auto tr_y_yyz_xzzzz = pbuffer.data(idx_dip_fh + 371);

    auto tr_y_yyz_yyyyy = pbuffer.data(idx_dip_fh + 372);

    auto tr_y_yyz_yyyyz = pbuffer.data(idx_dip_fh + 373);

    auto tr_y_yyz_yyyzz = pbuffer.data(idx_dip_fh + 374);

    auto tr_y_yyz_yyzzz = pbuffer.data(idx_dip_fh + 375);

    auto tr_y_yyz_yzzzz = pbuffer.data(idx_dip_fh + 376);

    auto tr_y_yyz_zzzzz = pbuffer.data(idx_dip_fh + 377);

#pragma omp simd aligned(pa_z,               \
                             tr_y_yy_xxxx,   \
                             tr_y_yy_xxxxx,  \
                             tr_y_yy_xxxxy,  \
                             tr_y_yy_xxxxz,  \
                             tr_y_yy_xxxy,   \
                             tr_y_yy_xxxyy,  \
                             tr_y_yy_xxxyz,  \
                             tr_y_yy_xxxz,   \
                             tr_y_yy_xxxzz,  \
                             tr_y_yy_xxyy,   \
                             tr_y_yy_xxyyy,  \
                             tr_y_yy_xxyyz,  \
                             tr_y_yy_xxyz,   \
                             tr_y_yy_xxyzz,  \
                             tr_y_yy_xxzz,   \
                             tr_y_yy_xxzzz,  \
                             tr_y_yy_xyyy,   \
                             tr_y_yy_xyyyy,  \
                             tr_y_yy_xyyyz,  \
                             tr_y_yy_xyyz,   \
                             tr_y_yy_xyyzz,  \
                             tr_y_yy_xyzz,   \
                             tr_y_yy_xyzzz,  \
                             tr_y_yy_xzzz,   \
                             tr_y_yy_xzzzz,  \
                             tr_y_yy_yyyy,   \
                             tr_y_yy_yyyyy,  \
                             tr_y_yy_yyyyz,  \
                             tr_y_yy_yyyz,   \
                             tr_y_yy_yyyzz,  \
                             tr_y_yy_yyzz,   \
                             tr_y_yy_yyzzz,  \
                             tr_y_yy_yzzz,   \
                             tr_y_yy_yzzzz,  \
                             tr_y_yy_zzzz,   \
                             tr_y_yy_zzzzz,  \
                             tr_y_yyz_xxxxx, \
                             tr_y_yyz_xxxxy, \
                             tr_y_yyz_xxxxz, \
                             tr_y_yyz_xxxyy, \
                             tr_y_yyz_xxxyz, \
                             tr_y_yyz_xxxzz, \
                             tr_y_yyz_xxyyy, \
                             tr_y_yyz_xxyyz, \
                             tr_y_yyz_xxyzz, \
                             tr_y_yyz_xxzzz, \
                             tr_y_yyz_xyyyy, \
                             tr_y_yyz_xyyyz, \
                             tr_y_yyz_xyyzz, \
                             tr_y_yyz_xyzzz, \
                             tr_y_yyz_xzzzz, \
                             tr_y_yyz_yyyyy, \
                             tr_y_yyz_yyyyz, \
                             tr_y_yyz_yyyzz, \
                             tr_y_yyz_yyzzz, \
                             tr_y_yyz_yzzzz, \
                             tr_y_yyz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yyz_xxxxx[i] = tr_y_yy_xxxxx[i] * pa_z[i];

        tr_y_yyz_xxxxy[i] = tr_y_yy_xxxxy[i] * pa_z[i];

        tr_y_yyz_xxxxz[i] = tr_y_yy_xxxx[i] * fe_0 + tr_y_yy_xxxxz[i] * pa_z[i];

        tr_y_yyz_xxxyy[i] = tr_y_yy_xxxyy[i] * pa_z[i];

        tr_y_yyz_xxxyz[i] = tr_y_yy_xxxy[i] * fe_0 + tr_y_yy_xxxyz[i] * pa_z[i];

        tr_y_yyz_xxxzz[i] = 2.0 * tr_y_yy_xxxz[i] * fe_0 + tr_y_yy_xxxzz[i] * pa_z[i];

        tr_y_yyz_xxyyy[i] = tr_y_yy_xxyyy[i] * pa_z[i];

        tr_y_yyz_xxyyz[i] = tr_y_yy_xxyy[i] * fe_0 + tr_y_yy_xxyyz[i] * pa_z[i];

        tr_y_yyz_xxyzz[i] = 2.0 * tr_y_yy_xxyz[i] * fe_0 + tr_y_yy_xxyzz[i] * pa_z[i];

        tr_y_yyz_xxzzz[i] = 3.0 * tr_y_yy_xxzz[i] * fe_0 + tr_y_yy_xxzzz[i] * pa_z[i];

        tr_y_yyz_xyyyy[i] = tr_y_yy_xyyyy[i] * pa_z[i];

        tr_y_yyz_xyyyz[i] = tr_y_yy_xyyy[i] * fe_0 + tr_y_yy_xyyyz[i] * pa_z[i];

        tr_y_yyz_xyyzz[i] = 2.0 * tr_y_yy_xyyz[i] * fe_0 + tr_y_yy_xyyzz[i] * pa_z[i];

        tr_y_yyz_xyzzz[i] = 3.0 * tr_y_yy_xyzz[i] * fe_0 + tr_y_yy_xyzzz[i] * pa_z[i];

        tr_y_yyz_xzzzz[i] = 4.0 * tr_y_yy_xzzz[i] * fe_0 + tr_y_yy_xzzzz[i] * pa_z[i];

        tr_y_yyz_yyyyy[i] = tr_y_yy_yyyyy[i] * pa_z[i];

        tr_y_yyz_yyyyz[i] = tr_y_yy_yyyy[i] * fe_0 + tr_y_yy_yyyyz[i] * pa_z[i];

        tr_y_yyz_yyyzz[i] = 2.0 * tr_y_yy_yyyz[i] * fe_0 + tr_y_yy_yyyzz[i] * pa_z[i];

        tr_y_yyz_yyzzz[i] = 3.0 * tr_y_yy_yyzz[i] * fe_0 + tr_y_yy_yyzzz[i] * pa_z[i];

        tr_y_yyz_yzzzz[i] = 4.0 * tr_y_yy_yzzz[i] * fe_0 + tr_y_yy_yzzzz[i] * pa_z[i];

        tr_y_yyz_zzzzz[i] = 5.0 * tr_y_yy_zzzz[i] * fe_0 + tr_y_yy_zzzzz[i] * pa_z[i];
    }

    // Set up 378-399 components of targeted buffer : FH

    auto tr_y_yzz_xxxxx = pbuffer.data(idx_dip_fh + 378);

    auto tr_y_yzz_xxxxy = pbuffer.data(idx_dip_fh + 379);

    auto tr_y_yzz_xxxxz = pbuffer.data(idx_dip_fh + 380);

    auto tr_y_yzz_xxxyy = pbuffer.data(idx_dip_fh + 381);

    auto tr_y_yzz_xxxyz = pbuffer.data(idx_dip_fh + 382);

    auto tr_y_yzz_xxxzz = pbuffer.data(idx_dip_fh + 383);

    auto tr_y_yzz_xxyyy = pbuffer.data(idx_dip_fh + 384);

    auto tr_y_yzz_xxyyz = pbuffer.data(idx_dip_fh + 385);

    auto tr_y_yzz_xxyzz = pbuffer.data(idx_dip_fh + 386);

    auto tr_y_yzz_xxzzz = pbuffer.data(idx_dip_fh + 387);

    auto tr_y_yzz_xyyyy = pbuffer.data(idx_dip_fh + 388);

    auto tr_y_yzz_xyyyz = pbuffer.data(idx_dip_fh + 389);

    auto tr_y_yzz_xyyzz = pbuffer.data(idx_dip_fh + 390);

    auto tr_y_yzz_xyzzz = pbuffer.data(idx_dip_fh + 391);

    auto tr_y_yzz_xzzzz = pbuffer.data(idx_dip_fh + 392);

    auto tr_y_yzz_yyyyy = pbuffer.data(idx_dip_fh + 393);

    auto tr_y_yzz_yyyyz = pbuffer.data(idx_dip_fh + 394);

    auto tr_y_yzz_yyyzz = pbuffer.data(idx_dip_fh + 395);

    auto tr_y_yzz_yyzzz = pbuffer.data(idx_dip_fh + 396);

    auto tr_y_yzz_yzzzz = pbuffer.data(idx_dip_fh + 397);

    auto tr_y_yzz_zzzzz = pbuffer.data(idx_dip_fh + 398);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_y_y_xxxxy,   \
                             tr_y_y_xxxyy,   \
                             tr_y_y_xxyyy,   \
                             tr_y_y_xyyyy,   \
                             tr_y_y_yyyyy,   \
                             tr_y_yz_xxxxy,  \
                             tr_y_yz_xxxyy,  \
                             tr_y_yz_xxyyy,  \
                             tr_y_yz_xyyyy,  \
                             tr_y_yz_yyyyy,  \
                             tr_y_yzz_xxxxx, \
                             tr_y_yzz_xxxxy, \
                             tr_y_yzz_xxxxz, \
                             tr_y_yzz_xxxyy, \
                             tr_y_yzz_xxxyz, \
                             tr_y_yzz_xxxzz, \
                             tr_y_yzz_xxyyy, \
                             tr_y_yzz_xxyyz, \
                             tr_y_yzz_xxyzz, \
                             tr_y_yzz_xxzzz, \
                             tr_y_yzz_xyyyy, \
                             tr_y_yzz_xyyyz, \
                             tr_y_yzz_xyyzz, \
                             tr_y_yzz_xyzzz, \
                             tr_y_yzz_xzzzz, \
                             tr_y_yzz_yyyyy, \
                             tr_y_yzz_yyyyz, \
                             tr_y_yzz_yyyzz, \
                             tr_y_yzz_yyzzz, \
                             tr_y_yzz_yzzzz, \
                             tr_y_yzz_zzzzz, \
                             tr_y_zz_xxxxx,  \
                             tr_y_zz_xxxxz,  \
                             tr_y_zz_xxxyz,  \
                             tr_y_zz_xxxz,   \
                             tr_y_zz_xxxzz,  \
                             tr_y_zz_xxyyz,  \
                             tr_y_zz_xxyz,   \
                             tr_y_zz_xxyzz,  \
                             tr_y_zz_xxzz,   \
                             tr_y_zz_xxzzz,  \
                             tr_y_zz_xyyyz,  \
                             tr_y_zz_xyyz,   \
                             tr_y_zz_xyyzz,  \
                             tr_y_zz_xyzz,   \
                             tr_y_zz_xyzzz,  \
                             tr_y_zz_xzzz,   \
                             tr_y_zz_xzzzz,  \
                             tr_y_zz_yyyyz,  \
                             tr_y_zz_yyyz,   \
                             tr_y_zz_yyyzz,  \
                             tr_y_zz_yyzz,   \
                             tr_y_zz_yyzzz,  \
                             tr_y_zz_yzzz,   \
                             tr_y_zz_yzzzz,  \
                             tr_y_zz_zzzz,   \
                             tr_y_zz_zzzzz,  \
                             ts_zz_xxxxx,    \
                             ts_zz_xxxxz,    \
                             ts_zz_xxxyz,    \
                             ts_zz_xxxzz,    \
                             ts_zz_xxyyz,    \
                             ts_zz_xxyzz,    \
                             ts_zz_xxzzz,    \
                             ts_zz_xyyyz,    \
                             ts_zz_xyyzz,    \
                             ts_zz_xyzzz,    \
                             ts_zz_xzzzz,    \
                             ts_zz_yyyyz,    \
                             ts_zz_yyyzz,    \
                             ts_zz_yyzzz,    \
                             ts_zz_yzzzz,    \
                             ts_zz_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yzz_xxxxx[i] = ts_zz_xxxxx[i] * fe_0 + tr_y_zz_xxxxx[i] * pa_y[i];

        tr_y_yzz_xxxxy[i] = tr_y_y_xxxxy[i] * fe_0 + tr_y_yz_xxxxy[i] * pa_z[i];

        tr_y_yzz_xxxxz[i] = ts_zz_xxxxz[i] * fe_0 + tr_y_zz_xxxxz[i] * pa_y[i];

        tr_y_yzz_xxxyy[i] = tr_y_y_xxxyy[i] * fe_0 + tr_y_yz_xxxyy[i] * pa_z[i];

        tr_y_yzz_xxxyz[i] = tr_y_zz_xxxz[i] * fe_0 + ts_zz_xxxyz[i] * fe_0 + tr_y_zz_xxxyz[i] * pa_y[i];

        tr_y_yzz_xxxzz[i] = ts_zz_xxxzz[i] * fe_0 + tr_y_zz_xxxzz[i] * pa_y[i];

        tr_y_yzz_xxyyy[i] = tr_y_y_xxyyy[i] * fe_0 + tr_y_yz_xxyyy[i] * pa_z[i];

        tr_y_yzz_xxyyz[i] = 2.0 * tr_y_zz_xxyz[i] * fe_0 + ts_zz_xxyyz[i] * fe_0 + tr_y_zz_xxyyz[i] * pa_y[i];

        tr_y_yzz_xxyzz[i] = tr_y_zz_xxzz[i] * fe_0 + ts_zz_xxyzz[i] * fe_0 + tr_y_zz_xxyzz[i] * pa_y[i];

        tr_y_yzz_xxzzz[i] = ts_zz_xxzzz[i] * fe_0 + tr_y_zz_xxzzz[i] * pa_y[i];

        tr_y_yzz_xyyyy[i] = tr_y_y_xyyyy[i] * fe_0 + tr_y_yz_xyyyy[i] * pa_z[i];

        tr_y_yzz_xyyyz[i] = 3.0 * tr_y_zz_xyyz[i] * fe_0 + ts_zz_xyyyz[i] * fe_0 + tr_y_zz_xyyyz[i] * pa_y[i];

        tr_y_yzz_xyyzz[i] = 2.0 * tr_y_zz_xyzz[i] * fe_0 + ts_zz_xyyzz[i] * fe_0 + tr_y_zz_xyyzz[i] * pa_y[i];

        tr_y_yzz_xyzzz[i] = tr_y_zz_xzzz[i] * fe_0 + ts_zz_xyzzz[i] * fe_0 + tr_y_zz_xyzzz[i] * pa_y[i];

        tr_y_yzz_xzzzz[i] = ts_zz_xzzzz[i] * fe_0 + tr_y_zz_xzzzz[i] * pa_y[i];

        tr_y_yzz_yyyyy[i] = tr_y_y_yyyyy[i] * fe_0 + tr_y_yz_yyyyy[i] * pa_z[i];

        tr_y_yzz_yyyyz[i] = 4.0 * tr_y_zz_yyyz[i] * fe_0 + ts_zz_yyyyz[i] * fe_0 + tr_y_zz_yyyyz[i] * pa_y[i];

        tr_y_yzz_yyyzz[i] = 3.0 * tr_y_zz_yyzz[i] * fe_0 + ts_zz_yyyzz[i] * fe_0 + tr_y_zz_yyyzz[i] * pa_y[i];

        tr_y_yzz_yyzzz[i] = 2.0 * tr_y_zz_yzzz[i] * fe_0 + ts_zz_yyzzz[i] * fe_0 + tr_y_zz_yyzzz[i] * pa_y[i];

        tr_y_yzz_yzzzz[i] = tr_y_zz_zzzz[i] * fe_0 + ts_zz_yzzzz[i] * fe_0 + tr_y_zz_yzzzz[i] * pa_y[i];

        tr_y_yzz_zzzzz[i] = ts_zz_zzzzz[i] * fe_0 + tr_y_zz_zzzzz[i] * pa_y[i];
    }

    // Set up 399-420 components of targeted buffer : FH

    auto tr_y_zzz_xxxxx = pbuffer.data(idx_dip_fh + 399);

    auto tr_y_zzz_xxxxy = pbuffer.data(idx_dip_fh + 400);

    auto tr_y_zzz_xxxxz = pbuffer.data(idx_dip_fh + 401);

    auto tr_y_zzz_xxxyy = pbuffer.data(idx_dip_fh + 402);

    auto tr_y_zzz_xxxyz = pbuffer.data(idx_dip_fh + 403);

    auto tr_y_zzz_xxxzz = pbuffer.data(idx_dip_fh + 404);

    auto tr_y_zzz_xxyyy = pbuffer.data(idx_dip_fh + 405);

    auto tr_y_zzz_xxyyz = pbuffer.data(idx_dip_fh + 406);

    auto tr_y_zzz_xxyzz = pbuffer.data(idx_dip_fh + 407);

    auto tr_y_zzz_xxzzz = pbuffer.data(idx_dip_fh + 408);

    auto tr_y_zzz_xyyyy = pbuffer.data(idx_dip_fh + 409);

    auto tr_y_zzz_xyyyz = pbuffer.data(idx_dip_fh + 410);

    auto tr_y_zzz_xyyzz = pbuffer.data(idx_dip_fh + 411);

    auto tr_y_zzz_xyzzz = pbuffer.data(idx_dip_fh + 412);

    auto tr_y_zzz_xzzzz = pbuffer.data(idx_dip_fh + 413);

    auto tr_y_zzz_yyyyy = pbuffer.data(idx_dip_fh + 414);

    auto tr_y_zzz_yyyyz = pbuffer.data(idx_dip_fh + 415);

    auto tr_y_zzz_yyyzz = pbuffer.data(idx_dip_fh + 416);

    auto tr_y_zzz_yyzzz = pbuffer.data(idx_dip_fh + 417);

    auto tr_y_zzz_yzzzz = pbuffer.data(idx_dip_fh + 418);

    auto tr_y_zzz_zzzzz = pbuffer.data(idx_dip_fh + 419);

#pragma omp simd aligned(pa_z,               \
                             tr_y_z_xxxxx,   \
                             tr_y_z_xxxxy,   \
                             tr_y_z_xxxxz,   \
                             tr_y_z_xxxyy,   \
                             tr_y_z_xxxyz,   \
                             tr_y_z_xxxzz,   \
                             tr_y_z_xxyyy,   \
                             tr_y_z_xxyyz,   \
                             tr_y_z_xxyzz,   \
                             tr_y_z_xxzzz,   \
                             tr_y_z_xyyyy,   \
                             tr_y_z_xyyyz,   \
                             tr_y_z_xyyzz,   \
                             tr_y_z_xyzzz,   \
                             tr_y_z_xzzzz,   \
                             tr_y_z_yyyyy,   \
                             tr_y_z_yyyyz,   \
                             tr_y_z_yyyzz,   \
                             tr_y_z_yyzzz,   \
                             tr_y_z_yzzzz,   \
                             tr_y_z_zzzzz,   \
                             tr_y_zz_xxxx,   \
                             tr_y_zz_xxxxx,  \
                             tr_y_zz_xxxxy,  \
                             tr_y_zz_xxxxz,  \
                             tr_y_zz_xxxy,   \
                             tr_y_zz_xxxyy,  \
                             tr_y_zz_xxxyz,  \
                             tr_y_zz_xxxz,   \
                             tr_y_zz_xxxzz,  \
                             tr_y_zz_xxyy,   \
                             tr_y_zz_xxyyy,  \
                             tr_y_zz_xxyyz,  \
                             tr_y_zz_xxyz,   \
                             tr_y_zz_xxyzz,  \
                             tr_y_zz_xxzz,   \
                             tr_y_zz_xxzzz,  \
                             tr_y_zz_xyyy,   \
                             tr_y_zz_xyyyy,  \
                             tr_y_zz_xyyyz,  \
                             tr_y_zz_xyyz,   \
                             tr_y_zz_xyyzz,  \
                             tr_y_zz_xyzz,   \
                             tr_y_zz_xyzzz,  \
                             tr_y_zz_xzzz,   \
                             tr_y_zz_xzzzz,  \
                             tr_y_zz_yyyy,   \
                             tr_y_zz_yyyyy,  \
                             tr_y_zz_yyyyz,  \
                             tr_y_zz_yyyz,   \
                             tr_y_zz_yyyzz,  \
                             tr_y_zz_yyzz,   \
                             tr_y_zz_yyzzz,  \
                             tr_y_zz_yzzz,   \
                             tr_y_zz_yzzzz,  \
                             tr_y_zz_zzzz,   \
                             tr_y_zz_zzzzz,  \
                             tr_y_zzz_xxxxx, \
                             tr_y_zzz_xxxxy, \
                             tr_y_zzz_xxxxz, \
                             tr_y_zzz_xxxyy, \
                             tr_y_zzz_xxxyz, \
                             tr_y_zzz_xxxzz, \
                             tr_y_zzz_xxyyy, \
                             tr_y_zzz_xxyyz, \
                             tr_y_zzz_xxyzz, \
                             tr_y_zzz_xxzzz, \
                             tr_y_zzz_xyyyy, \
                             tr_y_zzz_xyyyz, \
                             tr_y_zzz_xyyzz, \
                             tr_y_zzz_xyzzz, \
                             tr_y_zzz_xzzzz, \
                             tr_y_zzz_yyyyy, \
                             tr_y_zzz_yyyyz, \
                             tr_y_zzz_yyyzz, \
                             tr_y_zzz_yyzzz, \
                             tr_y_zzz_yzzzz, \
                             tr_y_zzz_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zzz_xxxxx[i] = 2.0 * tr_y_z_xxxxx[i] * fe_0 + tr_y_zz_xxxxx[i] * pa_z[i];

        tr_y_zzz_xxxxy[i] = 2.0 * tr_y_z_xxxxy[i] * fe_0 + tr_y_zz_xxxxy[i] * pa_z[i];

        tr_y_zzz_xxxxz[i] = 2.0 * tr_y_z_xxxxz[i] * fe_0 + tr_y_zz_xxxx[i] * fe_0 + tr_y_zz_xxxxz[i] * pa_z[i];

        tr_y_zzz_xxxyy[i] = 2.0 * tr_y_z_xxxyy[i] * fe_0 + tr_y_zz_xxxyy[i] * pa_z[i];

        tr_y_zzz_xxxyz[i] = 2.0 * tr_y_z_xxxyz[i] * fe_0 + tr_y_zz_xxxy[i] * fe_0 + tr_y_zz_xxxyz[i] * pa_z[i];

        tr_y_zzz_xxxzz[i] = 2.0 * tr_y_z_xxxzz[i] * fe_0 + 2.0 * tr_y_zz_xxxz[i] * fe_0 + tr_y_zz_xxxzz[i] * pa_z[i];

        tr_y_zzz_xxyyy[i] = 2.0 * tr_y_z_xxyyy[i] * fe_0 + tr_y_zz_xxyyy[i] * pa_z[i];

        tr_y_zzz_xxyyz[i] = 2.0 * tr_y_z_xxyyz[i] * fe_0 + tr_y_zz_xxyy[i] * fe_0 + tr_y_zz_xxyyz[i] * pa_z[i];

        tr_y_zzz_xxyzz[i] = 2.0 * tr_y_z_xxyzz[i] * fe_0 + 2.0 * tr_y_zz_xxyz[i] * fe_0 + tr_y_zz_xxyzz[i] * pa_z[i];

        tr_y_zzz_xxzzz[i] = 2.0 * tr_y_z_xxzzz[i] * fe_0 + 3.0 * tr_y_zz_xxzz[i] * fe_0 + tr_y_zz_xxzzz[i] * pa_z[i];

        tr_y_zzz_xyyyy[i] = 2.0 * tr_y_z_xyyyy[i] * fe_0 + tr_y_zz_xyyyy[i] * pa_z[i];

        tr_y_zzz_xyyyz[i] = 2.0 * tr_y_z_xyyyz[i] * fe_0 + tr_y_zz_xyyy[i] * fe_0 + tr_y_zz_xyyyz[i] * pa_z[i];

        tr_y_zzz_xyyzz[i] = 2.0 * tr_y_z_xyyzz[i] * fe_0 + 2.0 * tr_y_zz_xyyz[i] * fe_0 + tr_y_zz_xyyzz[i] * pa_z[i];

        tr_y_zzz_xyzzz[i] = 2.0 * tr_y_z_xyzzz[i] * fe_0 + 3.0 * tr_y_zz_xyzz[i] * fe_0 + tr_y_zz_xyzzz[i] * pa_z[i];

        tr_y_zzz_xzzzz[i] = 2.0 * tr_y_z_xzzzz[i] * fe_0 + 4.0 * tr_y_zz_xzzz[i] * fe_0 + tr_y_zz_xzzzz[i] * pa_z[i];

        tr_y_zzz_yyyyy[i] = 2.0 * tr_y_z_yyyyy[i] * fe_0 + tr_y_zz_yyyyy[i] * pa_z[i];

        tr_y_zzz_yyyyz[i] = 2.0 * tr_y_z_yyyyz[i] * fe_0 + tr_y_zz_yyyy[i] * fe_0 + tr_y_zz_yyyyz[i] * pa_z[i];

        tr_y_zzz_yyyzz[i] = 2.0 * tr_y_z_yyyzz[i] * fe_0 + 2.0 * tr_y_zz_yyyz[i] * fe_0 + tr_y_zz_yyyzz[i] * pa_z[i];

        tr_y_zzz_yyzzz[i] = 2.0 * tr_y_z_yyzzz[i] * fe_0 + 3.0 * tr_y_zz_yyzz[i] * fe_0 + tr_y_zz_yyzzz[i] * pa_z[i];

        tr_y_zzz_yzzzz[i] = 2.0 * tr_y_z_yzzzz[i] * fe_0 + 4.0 * tr_y_zz_yzzz[i] * fe_0 + tr_y_zz_yzzzz[i] * pa_z[i];

        tr_y_zzz_zzzzz[i] = 2.0 * tr_y_z_zzzzz[i] * fe_0 + 5.0 * tr_y_zz_zzzz[i] * fe_0 + tr_y_zz_zzzzz[i] * pa_z[i];
    }

    // Set up 420-441 components of targeted buffer : FH

    auto tr_z_xxx_xxxxx = pbuffer.data(idx_dip_fh + 420);

    auto tr_z_xxx_xxxxy = pbuffer.data(idx_dip_fh + 421);

    auto tr_z_xxx_xxxxz = pbuffer.data(idx_dip_fh + 422);

    auto tr_z_xxx_xxxyy = pbuffer.data(idx_dip_fh + 423);

    auto tr_z_xxx_xxxyz = pbuffer.data(idx_dip_fh + 424);

    auto tr_z_xxx_xxxzz = pbuffer.data(idx_dip_fh + 425);

    auto tr_z_xxx_xxyyy = pbuffer.data(idx_dip_fh + 426);

    auto tr_z_xxx_xxyyz = pbuffer.data(idx_dip_fh + 427);

    auto tr_z_xxx_xxyzz = pbuffer.data(idx_dip_fh + 428);

    auto tr_z_xxx_xxzzz = pbuffer.data(idx_dip_fh + 429);

    auto tr_z_xxx_xyyyy = pbuffer.data(idx_dip_fh + 430);

    auto tr_z_xxx_xyyyz = pbuffer.data(idx_dip_fh + 431);

    auto tr_z_xxx_xyyzz = pbuffer.data(idx_dip_fh + 432);

    auto tr_z_xxx_xyzzz = pbuffer.data(idx_dip_fh + 433);

    auto tr_z_xxx_xzzzz = pbuffer.data(idx_dip_fh + 434);

    auto tr_z_xxx_yyyyy = pbuffer.data(idx_dip_fh + 435);

    auto tr_z_xxx_yyyyz = pbuffer.data(idx_dip_fh + 436);

    auto tr_z_xxx_yyyzz = pbuffer.data(idx_dip_fh + 437);

    auto tr_z_xxx_yyzzz = pbuffer.data(idx_dip_fh + 438);

    auto tr_z_xxx_yzzzz = pbuffer.data(idx_dip_fh + 439);

    auto tr_z_xxx_zzzzz = pbuffer.data(idx_dip_fh + 440);

#pragma omp simd aligned(pa_x,               \
                             tr_z_x_xxxxx,   \
                             tr_z_x_xxxxy,   \
                             tr_z_x_xxxxz,   \
                             tr_z_x_xxxyy,   \
                             tr_z_x_xxxyz,   \
                             tr_z_x_xxxzz,   \
                             tr_z_x_xxyyy,   \
                             tr_z_x_xxyyz,   \
                             tr_z_x_xxyzz,   \
                             tr_z_x_xxzzz,   \
                             tr_z_x_xyyyy,   \
                             tr_z_x_xyyyz,   \
                             tr_z_x_xyyzz,   \
                             tr_z_x_xyzzz,   \
                             tr_z_x_xzzzz,   \
                             tr_z_x_yyyyy,   \
                             tr_z_x_yyyyz,   \
                             tr_z_x_yyyzz,   \
                             tr_z_x_yyzzz,   \
                             tr_z_x_yzzzz,   \
                             tr_z_x_zzzzz,   \
                             tr_z_xx_xxxx,   \
                             tr_z_xx_xxxxx,  \
                             tr_z_xx_xxxxy,  \
                             tr_z_xx_xxxxz,  \
                             tr_z_xx_xxxy,   \
                             tr_z_xx_xxxyy,  \
                             tr_z_xx_xxxyz,  \
                             tr_z_xx_xxxz,   \
                             tr_z_xx_xxxzz,  \
                             tr_z_xx_xxyy,   \
                             tr_z_xx_xxyyy,  \
                             tr_z_xx_xxyyz,  \
                             tr_z_xx_xxyz,   \
                             tr_z_xx_xxyzz,  \
                             tr_z_xx_xxzz,   \
                             tr_z_xx_xxzzz,  \
                             tr_z_xx_xyyy,   \
                             tr_z_xx_xyyyy,  \
                             tr_z_xx_xyyyz,  \
                             tr_z_xx_xyyz,   \
                             tr_z_xx_xyyzz,  \
                             tr_z_xx_xyzz,   \
                             tr_z_xx_xyzzz,  \
                             tr_z_xx_xzzz,   \
                             tr_z_xx_xzzzz,  \
                             tr_z_xx_yyyy,   \
                             tr_z_xx_yyyyy,  \
                             tr_z_xx_yyyyz,  \
                             tr_z_xx_yyyz,   \
                             tr_z_xx_yyyzz,  \
                             tr_z_xx_yyzz,   \
                             tr_z_xx_yyzzz,  \
                             tr_z_xx_yzzz,   \
                             tr_z_xx_yzzzz,  \
                             tr_z_xx_zzzz,   \
                             tr_z_xx_zzzzz,  \
                             tr_z_xxx_xxxxx, \
                             tr_z_xxx_xxxxy, \
                             tr_z_xxx_xxxxz, \
                             tr_z_xxx_xxxyy, \
                             tr_z_xxx_xxxyz, \
                             tr_z_xxx_xxxzz, \
                             tr_z_xxx_xxyyy, \
                             tr_z_xxx_xxyyz, \
                             tr_z_xxx_xxyzz, \
                             tr_z_xxx_xxzzz, \
                             tr_z_xxx_xyyyy, \
                             tr_z_xxx_xyyyz, \
                             tr_z_xxx_xyyzz, \
                             tr_z_xxx_xyzzz, \
                             tr_z_xxx_xzzzz, \
                             tr_z_xxx_yyyyy, \
                             tr_z_xxx_yyyyz, \
                             tr_z_xxx_yyyzz, \
                             tr_z_xxx_yyzzz, \
                             tr_z_xxx_yzzzz, \
                             tr_z_xxx_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxx_xxxxx[i] = 2.0 * tr_z_x_xxxxx[i] * fe_0 + 5.0 * tr_z_xx_xxxx[i] * fe_0 + tr_z_xx_xxxxx[i] * pa_x[i];

        tr_z_xxx_xxxxy[i] = 2.0 * tr_z_x_xxxxy[i] * fe_0 + 4.0 * tr_z_xx_xxxy[i] * fe_0 + tr_z_xx_xxxxy[i] * pa_x[i];

        tr_z_xxx_xxxxz[i] = 2.0 * tr_z_x_xxxxz[i] * fe_0 + 4.0 * tr_z_xx_xxxz[i] * fe_0 + tr_z_xx_xxxxz[i] * pa_x[i];

        tr_z_xxx_xxxyy[i] = 2.0 * tr_z_x_xxxyy[i] * fe_0 + 3.0 * tr_z_xx_xxyy[i] * fe_0 + tr_z_xx_xxxyy[i] * pa_x[i];

        tr_z_xxx_xxxyz[i] = 2.0 * tr_z_x_xxxyz[i] * fe_0 + 3.0 * tr_z_xx_xxyz[i] * fe_0 + tr_z_xx_xxxyz[i] * pa_x[i];

        tr_z_xxx_xxxzz[i] = 2.0 * tr_z_x_xxxzz[i] * fe_0 + 3.0 * tr_z_xx_xxzz[i] * fe_0 + tr_z_xx_xxxzz[i] * pa_x[i];

        tr_z_xxx_xxyyy[i] = 2.0 * tr_z_x_xxyyy[i] * fe_0 + 2.0 * tr_z_xx_xyyy[i] * fe_0 + tr_z_xx_xxyyy[i] * pa_x[i];

        tr_z_xxx_xxyyz[i] = 2.0 * tr_z_x_xxyyz[i] * fe_0 + 2.0 * tr_z_xx_xyyz[i] * fe_0 + tr_z_xx_xxyyz[i] * pa_x[i];

        tr_z_xxx_xxyzz[i] = 2.0 * tr_z_x_xxyzz[i] * fe_0 + 2.0 * tr_z_xx_xyzz[i] * fe_0 + tr_z_xx_xxyzz[i] * pa_x[i];

        tr_z_xxx_xxzzz[i] = 2.0 * tr_z_x_xxzzz[i] * fe_0 + 2.0 * tr_z_xx_xzzz[i] * fe_0 + tr_z_xx_xxzzz[i] * pa_x[i];

        tr_z_xxx_xyyyy[i] = 2.0 * tr_z_x_xyyyy[i] * fe_0 + tr_z_xx_yyyy[i] * fe_0 + tr_z_xx_xyyyy[i] * pa_x[i];

        tr_z_xxx_xyyyz[i] = 2.0 * tr_z_x_xyyyz[i] * fe_0 + tr_z_xx_yyyz[i] * fe_0 + tr_z_xx_xyyyz[i] * pa_x[i];

        tr_z_xxx_xyyzz[i] = 2.0 * tr_z_x_xyyzz[i] * fe_0 + tr_z_xx_yyzz[i] * fe_0 + tr_z_xx_xyyzz[i] * pa_x[i];

        tr_z_xxx_xyzzz[i] = 2.0 * tr_z_x_xyzzz[i] * fe_0 + tr_z_xx_yzzz[i] * fe_0 + tr_z_xx_xyzzz[i] * pa_x[i];

        tr_z_xxx_xzzzz[i] = 2.0 * tr_z_x_xzzzz[i] * fe_0 + tr_z_xx_zzzz[i] * fe_0 + tr_z_xx_xzzzz[i] * pa_x[i];

        tr_z_xxx_yyyyy[i] = 2.0 * tr_z_x_yyyyy[i] * fe_0 + tr_z_xx_yyyyy[i] * pa_x[i];

        tr_z_xxx_yyyyz[i] = 2.0 * tr_z_x_yyyyz[i] * fe_0 + tr_z_xx_yyyyz[i] * pa_x[i];

        tr_z_xxx_yyyzz[i] = 2.0 * tr_z_x_yyyzz[i] * fe_0 + tr_z_xx_yyyzz[i] * pa_x[i];

        tr_z_xxx_yyzzz[i] = 2.0 * tr_z_x_yyzzz[i] * fe_0 + tr_z_xx_yyzzz[i] * pa_x[i];

        tr_z_xxx_yzzzz[i] = 2.0 * tr_z_x_yzzzz[i] * fe_0 + tr_z_xx_yzzzz[i] * pa_x[i];

        tr_z_xxx_zzzzz[i] = 2.0 * tr_z_x_zzzzz[i] * fe_0 + tr_z_xx_zzzzz[i] * pa_x[i];
    }

    // Set up 441-462 components of targeted buffer : FH

    auto tr_z_xxy_xxxxx = pbuffer.data(idx_dip_fh + 441);

    auto tr_z_xxy_xxxxy = pbuffer.data(idx_dip_fh + 442);

    auto tr_z_xxy_xxxxz = pbuffer.data(idx_dip_fh + 443);

    auto tr_z_xxy_xxxyy = pbuffer.data(idx_dip_fh + 444);

    auto tr_z_xxy_xxxyz = pbuffer.data(idx_dip_fh + 445);

    auto tr_z_xxy_xxxzz = pbuffer.data(idx_dip_fh + 446);

    auto tr_z_xxy_xxyyy = pbuffer.data(idx_dip_fh + 447);

    auto tr_z_xxy_xxyyz = pbuffer.data(idx_dip_fh + 448);

    auto tr_z_xxy_xxyzz = pbuffer.data(idx_dip_fh + 449);

    auto tr_z_xxy_xxzzz = pbuffer.data(idx_dip_fh + 450);

    auto tr_z_xxy_xyyyy = pbuffer.data(idx_dip_fh + 451);

    auto tr_z_xxy_xyyyz = pbuffer.data(idx_dip_fh + 452);

    auto tr_z_xxy_xyyzz = pbuffer.data(idx_dip_fh + 453);

    auto tr_z_xxy_xyzzz = pbuffer.data(idx_dip_fh + 454);

    auto tr_z_xxy_xzzzz = pbuffer.data(idx_dip_fh + 455);

    auto tr_z_xxy_yyyyy = pbuffer.data(idx_dip_fh + 456);

    auto tr_z_xxy_yyyyz = pbuffer.data(idx_dip_fh + 457);

    auto tr_z_xxy_yyyzz = pbuffer.data(idx_dip_fh + 458);

    auto tr_z_xxy_yyzzz = pbuffer.data(idx_dip_fh + 459);

    auto tr_z_xxy_yzzzz = pbuffer.data(idx_dip_fh + 460);

    auto tr_z_xxy_zzzzz = pbuffer.data(idx_dip_fh + 461);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xx_xxxx,   \
                             tr_z_xx_xxxxx,  \
                             tr_z_xx_xxxxy,  \
                             tr_z_xx_xxxxz,  \
                             tr_z_xx_xxxy,   \
                             tr_z_xx_xxxyy,  \
                             tr_z_xx_xxxyz,  \
                             tr_z_xx_xxxz,   \
                             tr_z_xx_xxxzz,  \
                             tr_z_xx_xxyy,   \
                             tr_z_xx_xxyyy,  \
                             tr_z_xx_xxyyz,  \
                             tr_z_xx_xxyz,   \
                             tr_z_xx_xxyzz,  \
                             tr_z_xx_xxzz,   \
                             tr_z_xx_xxzzz,  \
                             tr_z_xx_xyyy,   \
                             tr_z_xx_xyyyy,  \
                             tr_z_xx_xyyyz,  \
                             tr_z_xx_xyyz,   \
                             tr_z_xx_xyyzz,  \
                             tr_z_xx_xyzz,   \
                             tr_z_xx_xyzzz,  \
                             tr_z_xx_xzzz,   \
                             tr_z_xx_xzzzz,  \
                             tr_z_xx_zzzzz,  \
                             tr_z_xxy_xxxxx, \
                             tr_z_xxy_xxxxy, \
                             tr_z_xxy_xxxxz, \
                             tr_z_xxy_xxxyy, \
                             tr_z_xxy_xxxyz, \
                             tr_z_xxy_xxxzz, \
                             tr_z_xxy_xxyyy, \
                             tr_z_xxy_xxyyz, \
                             tr_z_xxy_xxyzz, \
                             tr_z_xxy_xxzzz, \
                             tr_z_xxy_xyyyy, \
                             tr_z_xxy_xyyyz, \
                             tr_z_xxy_xyyzz, \
                             tr_z_xxy_xyzzz, \
                             tr_z_xxy_xzzzz, \
                             tr_z_xxy_yyyyy, \
                             tr_z_xxy_yyyyz, \
                             tr_z_xxy_yyyzz, \
                             tr_z_xxy_yyzzz, \
                             tr_z_xxy_yzzzz, \
                             tr_z_xxy_zzzzz, \
                             tr_z_xy_yyyyy,  \
                             tr_z_xy_yyyyz,  \
                             tr_z_xy_yyyzz,  \
                             tr_z_xy_yyzzz,  \
                             tr_z_xy_yzzzz,  \
                             tr_z_y_yyyyy,   \
                             tr_z_y_yyyyz,   \
                             tr_z_y_yyyzz,   \
                             tr_z_y_yyzzz,   \
                             tr_z_y_yzzzz,   \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxy_xxxxx[i] = tr_z_xx_xxxxx[i] * pa_y[i];

        tr_z_xxy_xxxxy[i] = tr_z_xx_xxxx[i] * fe_0 + tr_z_xx_xxxxy[i] * pa_y[i];

        tr_z_xxy_xxxxz[i] = tr_z_xx_xxxxz[i] * pa_y[i];

        tr_z_xxy_xxxyy[i] = 2.0 * tr_z_xx_xxxy[i] * fe_0 + tr_z_xx_xxxyy[i] * pa_y[i];

        tr_z_xxy_xxxyz[i] = tr_z_xx_xxxz[i] * fe_0 + tr_z_xx_xxxyz[i] * pa_y[i];

        tr_z_xxy_xxxzz[i] = tr_z_xx_xxxzz[i] * pa_y[i];

        tr_z_xxy_xxyyy[i] = 3.0 * tr_z_xx_xxyy[i] * fe_0 + tr_z_xx_xxyyy[i] * pa_y[i];

        tr_z_xxy_xxyyz[i] = 2.0 * tr_z_xx_xxyz[i] * fe_0 + tr_z_xx_xxyyz[i] * pa_y[i];

        tr_z_xxy_xxyzz[i] = tr_z_xx_xxzz[i] * fe_0 + tr_z_xx_xxyzz[i] * pa_y[i];

        tr_z_xxy_xxzzz[i] = tr_z_xx_xxzzz[i] * pa_y[i];

        tr_z_xxy_xyyyy[i] = 4.0 * tr_z_xx_xyyy[i] * fe_0 + tr_z_xx_xyyyy[i] * pa_y[i];

        tr_z_xxy_xyyyz[i] = 3.0 * tr_z_xx_xyyz[i] * fe_0 + tr_z_xx_xyyyz[i] * pa_y[i];

        tr_z_xxy_xyyzz[i] = 2.0 * tr_z_xx_xyzz[i] * fe_0 + tr_z_xx_xyyzz[i] * pa_y[i];

        tr_z_xxy_xyzzz[i] = tr_z_xx_xzzz[i] * fe_0 + tr_z_xx_xyzzz[i] * pa_y[i];

        tr_z_xxy_xzzzz[i] = tr_z_xx_xzzzz[i] * pa_y[i];

        tr_z_xxy_yyyyy[i] = tr_z_y_yyyyy[i] * fe_0 + tr_z_xy_yyyyy[i] * pa_x[i];

        tr_z_xxy_yyyyz[i] = tr_z_y_yyyyz[i] * fe_0 + tr_z_xy_yyyyz[i] * pa_x[i];

        tr_z_xxy_yyyzz[i] = tr_z_y_yyyzz[i] * fe_0 + tr_z_xy_yyyzz[i] * pa_x[i];

        tr_z_xxy_yyzzz[i] = tr_z_y_yyzzz[i] * fe_0 + tr_z_xy_yyzzz[i] * pa_x[i];

        tr_z_xxy_yzzzz[i] = tr_z_y_yzzzz[i] * fe_0 + tr_z_xy_yzzzz[i] * pa_x[i];

        tr_z_xxy_zzzzz[i] = tr_z_xx_zzzzz[i] * pa_y[i];
    }

    // Set up 462-483 components of targeted buffer : FH

    auto tr_z_xxz_xxxxx = pbuffer.data(idx_dip_fh + 462);

    auto tr_z_xxz_xxxxy = pbuffer.data(idx_dip_fh + 463);

    auto tr_z_xxz_xxxxz = pbuffer.data(idx_dip_fh + 464);

    auto tr_z_xxz_xxxyy = pbuffer.data(idx_dip_fh + 465);

    auto tr_z_xxz_xxxyz = pbuffer.data(idx_dip_fh + 466);

    auto tr_z_xxz_xxxzz = pbuffer.data(idx_dip_fh + 467);

    auto tr_z_xxz_xxyyy = pbuffer.data(idx_dip_fh + 468);

    auto tr_z_xxz_xxyyz = pbuffer.data(idx_dip_fh + 469);

    auto tr_z_xxz_xxyzz = pbuffer.data(idx_dip_fh + 470);

    auto tr_z_xxz_xxzzz = pbuffer.data(idx_dip_fh + 471);

    auto tr_z_xxz_xyyyy = pbuffer.data(idx_dip_fh + 472);

    auto tr_z_xxz_xyyyz = pbuffer.data(idx_dip_fh + 473);

    auto tr_z_xxz_xyyzz = pbuffer.data(idx_dip_fh + 474);

    auto tr_z_xxz_xyzzz = pbuffer.data(idx_dip_fh + 475);

    auto tr_z_xxz_xzzzz = pbuffer.data(idx_dip_fh + 476);

    auto tr_z_xxz_yyyyy = pbuffer.data(idx_dip_fh + 477);

    auto tr_z_xxz_yyyyz = pbuffer.data(idx_dip_fh + 478);

    auto tr_z_xxz_yyyzz = pbuffer.data(idx_dip_fh + 479);

    auto tr_z_xxz_yyzzz = pbuffer.data(idx_dip_fh + 480);

    auto tr_z_xxz_yzzzz = pbuffer.data(idx_dip_fh + 481);

    auto tr_z_xxz_zzzzz = pbuffer.data(idx_dip_fh + 482);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_z_xx_xxxxx,  \
                             tr_z_xx_xxxxy,  \
                             tr_z_xx_xxxyy,  \
                             tr_z_xx_xxyyy,  \
                             tr_z_xx_xyyyy,  \
                             tr_z_xxz_xxxxx, \
                             tr_z_xxz_xxxxy, \
                             tr_z_xxz_xxxxz, \
                             tr_z_xxz_xxxyy, \
                             tr_z_xxz_xxxyz, \
                             tr_z_xxz_xxxzz, \
                             tr_z_xxz_xxyyy, \
                             tr_z_xxz_xxyyz, \
                             tr_z_xxz_xxyzz, \
                             tr_z_xxz_xxzzz, \
                             tr_z_xxz_xyyyy, \
                             tr_z_xxz_xyyyz, \
                             tr_z_xxz_xyyzz, \
                             tr_z_xxz_xyzzz, \
                             tr_z_xxz_xzzzz, \
                             tr_z_xxz_yyyyy, \
                             tr_z_xxz_yyyyz, \
                             tr_z_xxz_yyyzz, \
                             tr_z_xxz_yyzzz, \
                             tr_z_xxz_yzzzz, \
                             tr_z_xxz_zzzzz, \
                             tr_z_xz_xxxxz,  \
                             tr_z_xz_xxxyz,  \
                             tr_z_xz_xxxz,   \
                             tr_z_xz_xxxzz,  \
                             tr_z_xz_xxyyz,  \
                             tr_z_xz_xxyz,   \
                             tr_z_xz_xxyzz,  \
                             tr_z_xz_xxzz,   \
                             tr_z_xz_xxzzz,  \
                             tr_z_xz_xyyyz,  \
                             tr_z_xz_xyyz,   \
                             tr_z_xz_xyyzz,  \
                             tr_z_xz_xyzz,   \
                             tr_z_xz_xyzzz,  \
                             tr_z_xz_xzzz,   \
                             tr_z_xz_xzzzz,  \
                             tr_z_xz_yyyyy,  \
                             tr_z_xz_yyyyz,  \
                             tr_z_xz_yyyz,   \
                             tr_z_xz_yyyzz,  \
                             tr_z_xz_yyzz,   \
                             tr_z_xz_yyzzz,  \
                             tr_z_xz_yzzz,   \
                             tr_z_xz_yzzzz,  \
                             tr_z_xz_zzzz,   \
                             tr_z_xz_zzzzz,  \
                             tr_z_z_xxxxz,   \
                             tr_z_z_xxxyz,   \
                             tr_z_z_xxxzz,   \
                             tr_z_z_xxyyz,   \
                             tr_z_z_xxyzz,   \
                             tr_z_z_xxzzz,   \
                             tr_z_z_xyyyz,   \
                             tr_z_z_xyyzz,   \
                             tr_z_z_xyzzz,   \
                             tr_z_z_xzzzz,   \
                             tr_z_z_yyyyy,   \
                             tr_z_z_yyyyz,   \
                             tr_z_z_yyyzz,   \
                             tr_z_z_yyzzz,   \
                             tr_z_z_yzzzz,   \
                             tr_z_z_zzzzz,   \
                             ts_xx_xxxxx,    \
                             ts_xx_xxxxy,    \
                             ts_xx_xxxyy,    \
                             ts_xx_xxyyy,    \
                             ts_xx_xyyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xxz_xxxxx[i] = ts_xx_xxxxx[i] * fe_0 + tr_z_xx_xxxxx[i] * pa_z[i];

        tr_z_xxz_xxxxy[i] = ts_xx_xxxxy[i] * fe_0 + tr_z_xx_xxxxy[i] * pa_z[i];

        tr_z_xxz_xxxxz[i] = tr_z_z_xxxxz[i] * fe_0 + 4.0 * tr_z_xz_xxxz[i] * fe_0 + tr_z_xz_xxxxz[i] * pa_x[i];

        tr_z_xxz_xxxyy[i] = ts_xx_xxxyy[i] * fe_0 + tr_z_xx_xxxyy[i] * pa_z[i];

        tr_z_xxz_xxxyz[i] = tr_z_z_xxxyz[i] * fe_0 + 3.0 * tr_z_xz_xxyz[i] * fe_0 + tr_z_xz_xxxyz[i] * pa_x[i];

        tr_z_xxz_xxxzz[i] = tr_z_z_xxxzz[i] * fe_0 + 3.0 * tr_z_xz_xxzz[i] * fe_0 + tr_z_xz_xxxzz[i] * pa_x[i];

        tr_z_xxz_xxyyy[i] = ts_xx_xxyyy[i] * fe_0 + tr_z_xx_xxyyy[i] * pa_z[i];

        tr_z_xxz_xxyyz[i] = tr_z_z_xxyyz[i] * fe_0 + 2.0 * tr_z_xz_xyyz[i] * fe_0 + tr_z_xz_xxyyz[i] * pa_x[i];

        tr_z_xxz_xxyzz[i] = tr_z_z_xxyzz[i] * fe_0 + 2.0 * tr_z_xz_xyzz[i] * fe_0 + tr_z_xz_xxyzz[i] * pa_x[i];

        tr_z_xxz_xxzzz[i] = tr_z_z_xxzzz[i] * fe_0 + 2.0 * tr_z_xz_xzzz[i] * fe_0 + tr_z_xz_xxzzz[i] * pa_x[i];

        tr_z_xxz_xyyyy[i] = ts_xx_xyyyy[i] * fe_0 + tr_z_xx_xyyyy[i] * pa_z[i];

        tr_z_xxz_xyyyz[i] = tr_z_z_xyyyz[i] * fe_0 + tr_z_xz_yyyz[i] * fe_0 + tr_z_xz_xyyyz[i] * pa_x[i];

        tr_z_xxz_xyyzz[i] = tr_z_z_xyyzz[i] * fe_0 + tr_z_xz_yyzz[i] * fe_0 + tr_z_xz_xyyzz[i] * pa_x[i];

        tr_z_xxz_xyzzz[i] = tr_z_z_xyzzz[i] * fe_0 + tr_z_xz_yzzz[i] * fe_0 + tr_z_xz_xyzzz[i] * pa_x[i];

        tr_z_xxz_xzzzz[i] = tr_z_z_xzzzz[i] * fe_0 + tr_z_xz_zzzz[i] * fe_0 + tr_z_xz_xzzzz[i] * pa_x[i];

        tr_z_xxz_yyyyy[i] = tr_z_z_yyyyy[i] * fe_0 + tr_z_xz_yyyyy[i] * pa_x[i];

        tr_z_xxz_yyyyz[i] = tr_z_z_yyyyz[i] * fe_0 + tr_z_xz_yyyyz[i] * pa_x[i];

        tr_z_xxz_yyyzz[i] = tr_z_z_yyyzz[i] * fe_0 + tr_z_xz_yyyzz[i] * pa_x[i];

        tr_z_xxz_yyzzz[i] = tr_z_z_yyzzz[i] * fe_0 + tr_z_xz_yyzzz[i] * pa_x[i];

        tr_z_xxz_yzzzz[i] = tr_z_z_yzzzz[i] * fe_0 + tr_z_xz_yzzzz[i] * pa_x[i];

        tr_z_xxz_zzzzz[i] = tr_z_z_zzzzz[i] * fe_0 + tr_z_xz_zzzzz[i] * pa_x[i];
    }

    // Set up 483-504 components of targeted buffer : FH

    auto tr_z_xyy_xxxxx = pbuffer.data(idx_dip_fh + 483);

    auto tr_z_xyy_xxxxy = pbuffer.data(idx_dip_fh + 484);

    auto tr_z_xyy_xxxxz = pbuffer.data(idx_dip_fh + 485);

    auto tr_z_xyy_xxxyy = pbuffer.data(idx_dip_fh + 486);

    auto tr_z_xyy_xxxyz = pbuffer.data(idx_dip_fh + 487);

    auto tr_z_xyy_xxxzz = pbuffer.data(idx_dip_fh + 488);

    auto tr_z_xyy_xxyyy = pbuffer.data(idx_dip_fh + 489);

    auto tr_z_xyy_xxyyz = pbuffer.data(idx_dip_fh + 490);

    auto tr_z_xyy_xxyzz = pbuffer.data(idx_dip_fh + 491);

    auto tr_z_xyy_xxzzz = pbuffer.data(idx_dip_fh + 492);

    auto tr_z_xyy_xyyyy = pbuffer.data(idx_dip_fh + 493);

    auto tr_z_xyy_xyyyz = pbuffer.data(idx_dip_fh + 494);

    auto tr_z_xyy_xyyzz = pbuffer.data(idx_dip_fh + 495);

    auto tr_z_xyy_xyzzz = pbuffer.data(idx_dip_fh + 496);

    auto tr_z_xyy_xzzzz = pbuffer.data(idx_dip_fh + 497);

    auto tr_z_xyy_yyyyy = pbuffer.data(idx_dip_fh + 498);

    auto tr_z_xyy_yyyyz = pbuffer.data(idx_dip_fh + 499);

    auto tr_z_xyy_yyyzz = pbuffer.data(idx_dip_fh + 500);

    auto tr_z_xyy_yyzzz = pbuffer.data(idx_dip_fh + 501);

    auto tr_z_xyy_yzzzz = pbuffer.data(idx_dip_fh + 502);

    auto tr_z_xyy_zzzzz = pbuffer.data(idx_dip_fh + 503);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xyy_xxxxx, \
                             tr_z_xyy_xxxxy, \
                             tr_z_xyy_xxxxz, \
                             tr_z_xyy_xxxyy, \
                             tr_z_xyy_xxxyz, \
                             tr_z_xyy_xxxzz, \
                             tr_z_xyy_xxyyy, \
                             tr_z_xyy_xxyyz, \
                             tr_z_xyy_xxyzz, \
                             tr_z_xyy_xxzzz, \
                             tr_z_xyy_xyyyy, \
                             tr_z_xyy_xyyyz, \
                             tr_z_xyy_xyyzz, \
                             tr_z_xyy_xyzzz, \
                             tr_z_xyy_xzzzz, \
                             tr_z_xyy_yyyyy, \
                             tr_z_xyy_yyyyz, \
                             tr_z_xyy_yyyzz, \
                             tr_z_xyy_yyzzz, \
                             tr_z_xyy_yzzzz, \
                             tr_z_xyy_zzzzz, \
                             tr_z_yy_xxxx,   \
                             tr_z_yy_xxxxx,  \
                             tr_z_yy_xxxxy,  \
                             tr_z_yy_xxxxz,  \
                             tr_z_yy_xxxy,   \
                             tr_z_yy_xxxyy,  \
                             tr_z_yy_xxxyz,  \
                             tr_z_yy_xxxz,   \
                             tr_z_yy_xxxzz,  \
                             tr_z_yy_xxyy,   \
                             tr_z_yy_xxyyy,  \
                             tr_z_yy_xxyyz,  \
                             tr_z_yy_xxyz,   \
                             tr_z_yy_xxyzz,  \
                             tr_z_yy_xxzz,   \
                             tr_z_yy_xxzzz,  \
                             tr_z_yy_xyyy,   \
                             tr_z_yy_xyyyy,  \
                             tr_z_yy_xyyyz,  \
                             tr_z_yy_xyyz,   \
                             tr_z_yy_xyyzz,  \
                             tr_z_yy_xyzz,   \
                             tr_z_yy_xyzzz,  \
                             tr_z_yy_xzzz,   \
                             tr_z_yy_xzzzz,  \
                             tr_z_yy_yyyy,   \
                             tr_z_yy_yyyyy,  \
                             tr_z_yy_yyyyz,  \
                             tr_z_yy_yyyz,   \
                             tr_z_yy_yyyzz,  \
                             tr_z_yy_yyzz,   \
                             tr_z_yy_yyzzz,  \
                             tr_z_yy_yzzz,   \
                             tr_z_yy_yzzzz,  \
                             tr_z_yy_zzzz,   \
                             tr_z_yy_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyy_xxxxx[i] = 5.0 * tr_z_yy_xxxx[i] * fe_0 + tr_z_yy_xxxxx[i] * pa_x[i];

        tr_z_xyy_xxxxy[i] = 4.0 * tr_z_yy_xxxy[i] * fe_0 + tr_z_yy_xxxxy[i] * pa_x[i];

        tr_z_xyy_xxxxz[i] = 4.0 * tr_z_yy_xxxz[i] * fe_0 + tr_z_yy_xxxxz[i] * pa_x[i];

        tr_z_xyy_xxxyy[i] = 3.0 * tr_z_yy_xxyy[i] * fe_0 + tr_z_yy_xxxyy[i] * pa_x[i];

        tr_z_xyy_xxxyz[i] = 3.0 * tr_z_yy_xxyz[i] * fe_0 + tr_z_yy_xxxyz[i] * pa_x[i];

        tr_z_xyy_xxxzz[i] = 3.0 * tr_z_yy_xxzz[i] * fe_0 + tr_z_yy_xxxzz[i] * pa_x[i];

        tr_z_xyy_xxyyy[i] = 2.0 * tr_z_yy_xyyy[i] * fe_0 + tr_z_yy_xxyyy[i] * pa_x[i];

        tr_z_xyy_xxyyz[i] = 2.0 * tr_z_yy_xyyz[i] * fe_0 + tr_z_yy_xxyyz[i] * pa_x[i];

        tr_z_xyy_xxyzz[i] = 2.0 * tr_z_yy_xyzz[i] * fe_0 + tr_z_yy_xxyzz[i] * pa_x[i];

        tr_z_xyy_xxzzz[i] = 2.0 * tr_z_yy_xzzz[i] * fe_0 + tr_z_yy_xxzzz[i] * pa_x[i];

        tr_z_xyy_xyyyy[i] = tr_z_yy_yyyy[i] * fe_0 + tr_z_yy_xyyyy[i] * pa_x[i];

        tr_z_xyy_xyyyz[i] = tr_z_yy_yyyz[i] * fe_0 + tr_z_yy_xyyyz[i] * pa_x[i];

        tr_z_xyy_xyyzz[i] = tr_z_yy_yyzz[i] * fe_0 + tr_z_yy_xyyzz[i] * pa_x[i];

        tr_z_xyy_xyzzz[i] = tr_z_yy_yzzz[i] * fe_0 + tr_z_yy_xyzzz[i] * pa_x[i];

        tr_z_xyy_xzzzz[i] = tr_z_yy_zzzz[i] * fe_0 + tr_z_yy_xzzzz[i] * pa_x[i];

        tr_z_xyy_yyyyy[i] = tr_z_yy_yyyyy[i] * pa_x[i];

        tr_z_xyy_yyyyz[i] = tr_z_yy_yyyyz[i] * pa_x[i];

        tr_z_xyy_yyyzz[i] = tr_z_yy_yyyzz[i] * pa_x[i];

        tr_z_xyy_yyzzz[i] = tr_z_yy_yyzzz[i] * pa_x[i];

        tr_z_xyy_yzzzz[i] = tr_z_yy_yzzzz[i] * pa_x[i];

        tr_z_xyy_zzzzz[i] = tr_z_yy_zzzzz[i] * pa_x[i];
    }

    // Set up 504-525 components of targeted buffer : FH

    auto tr_z_xyz_xxxxx = pbuffer.data(idx_dip_fh + 504);

    auto tr_z_xyz_xxxxy = pbuffer.data(idx_dip_fh + 505);

    auto tr_z_xyz_xxxxz = pbuffer.data(idx_dip_fh + 506);

    auto tr_z_xyz_xxxyy = pbuffer.data(idx_dip_fh + 507);

    auto tr_z_xyz_xxxyz = pbuffer.data(idx_dip_fh + 508);

    auto tr_z_xyz_xxxzz = pbuffer.data(idx_dip_fh + 509);

    auto tr_z_xyz_xxyyy = pbuffer.data(idx_dip_fh + 510);

    auto tr_z_xyz_xxyyz = pbuffer.data(idx_dip_fh + 511);

    auto tr_z_xyz_xxyzz = pbuffer.data(idx_dip_fh + 512);

    auto tr_z_xyz_xxzzz = pbuffer.data(idx_dip_fh + 513);

    auto tr_z_xyz_xyyyy = pbuffer.data(idx_dip_fh + 514);

    auto tr_z_xyz_xyyyz = pbuffer.data(idx_dip_fh + 515);

    auto tr_z_xyz_xyyzz = pbuffer.data(idx_dip_fh + 516);

    auto tr_z_xyz_xyzzz = pbuffer.data(idx_dip_fh + 517);

    auto tr_z_xyz_xzzzz = pbuffer.data(idx_dip_fh + 518);

    auto tr_z_xyz_yyyyy = pbuffer.data(idx_dip_fh + 519);

    auto tr_z_xyz_yyyyz = pbuffer.data(idx_dip_fh + 520);

    auto tr_z_xyz_yyyzz = pbuffer.data(idx_dip_fh + 521);

    auto tr_z_xyz_yyzzz = pbuffer.data(idx_dip_fh + 522);

    auto tr_z_xyz_yzzzz = pbuffer.data(idx_dip_fh + 523);

    auto tr_z_xyz_zzzzz = pbuffer.data(idx_dip_fh + 524);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_xyz_xxxxx, \
                             tr_z_xyz_xxxxy, \
                             tr_z_xyz_xxxxz, \
                             tr_z_xyz_xxxyy, \
                             tr_z_xyz_xxxyz, \
                             tr_z_xyz_xxxzz, \
                             tr_z_xyz_xxyyy, \
                             tr_z_xyz_xxyyz, \
                             tr_z_xyz_xxyzz, \
                             tr_z_xyz_xxzzz, \
                             tr_z_xyz_xyyyy, \
                             tr_z_xyz_xyyyz, \
                             tr_z_xyz_xyyzz, \
                             tr_z_xyz_xyzzz, \
                             tr_z_xyz_xzzzz, \
                             tr_z_xyz_yyyyy, \
                             tr_z_xyz_yyyyz, \
                             tr_z_xyz_yyyzz, \
                             tr_z_xyz_yyzzz, \
                             tr_z_xyz_yzzzz, \
                             tr_z_xyz_zzzzz, \
                             tr_z_xz_xxxxx,  \
                             tr_z_xz_xxxxz,  \
                             tr_z_xz_xxxzz,  \
                             tr_z_xz_xxzzz,  \
                             tr_z_xz_xzzzz,  \
                             tr_z_yz_xxxxy,  \
                             tr_z_yz_xxxy,   \
                             tr_z_yz_xxxyy,  \
                             tr_z_yz_xxxyz,  \
                             tr_z_yz_xxyy,   \
                             tr_z_yz_xxyyy,  \
                             tr_z_yz_xxyyz,  \
                             tr_z_yz_xxyz,   \
                             tr_z_yz_xxyzz,  \
                             tr_z_yz_xyyy,   \
                             tr_z_yz_xyyyy,  \
                             tr_z_yz_xyyyz,  \
                             tr_z_yz_xyyz,   \
                             tr_z_yz_xyyzz,  \
                             tr_z_yz_xyzz,   \
                             tr_z_yz_xyzzz,  \
                             tr_z_yz_yyyy,   \
                             tr_z_yz_yyyyy,  \
                             tr_z_yz_yyyyz,  \
                             tr_z_yz_yyyz,   \
                             tr_z_yz_yyyzz,  \
                             tr_z_yz_yyzz,   \
                             tr_z_yz_yyzzz,  \
                             tr_z_yz_yzzz,   \
                             tr_z_yz_yzzzz,  \
                             tr_z_yz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xyz_xxxxx[i] = tr_z_xz_xxxxx[i] * pa_y[i];

        tr_z_xyz_xxxxy[i] = 4.0 * tr_z_yz_xxxy[i] * fe_0 + tr_z_yz_xxxxy[i] * pa_x[i];

        tr_z_xyz_xxxxz[i] = tr_z_xz_xxxxz[i] * pa_y[i];

        tr_z_xyz_xxxyy[i] = 3.0 * tr_z_yz_xxyy[i] * fe_0 + tr_z_yz_xxxyy[i] * pa_x[i];

        tr_z_xyz_xxxyz[i] = 3.0 * tr_z_yz_xxyz[i] * fe_0 + tr_z_yz_xxxyz[i] * pa_x[i];

        tr_z_xyz_xxxzz[i] = tr_z_xz_xxxzz[i] * pa_y[i];

        tr_z_xyz_xxyyy[i] = 2.0 * tr_z_yz_xyyy[i] * fe_0 + tr_z_yz_xxyyy[i] * pa_x[i];

        tr_z_xyz_xxyyz[i] = 2.0 * tr_z_yz_xyyz[i] * fe_0 + tr_z_yz_xxyyz[i] * pa_x[i];

        tr_z_xyz_xxyzz[i] = 2.0 * tr_z_yz_xyzz[i] * fe_0 + tr_z_yz_xxyzz[i] * pa_x[i];

        tr_z_xyz_xxzzz[i] = tr_z_xz_xxzzz[i] * pa_y[i];

        tr_z_xyz_xyyyy[i] = tr_z_yz_yyyy[i] * fe_0 + tr_z_yz_xyyyy[i] * pa_x[i];

        tr_z_xyz_xyyyz[i] = tr_z_yz_yyyz[i] * fe_0 + tr_z_yz_xyyyz[i] * pa_x[i];

        tr_z_xyz_xyyzz[i] = tr_z_yz_yyzz[i] * fe_0 + tr_z_yz_xyyzz[i] * pa_x[i];

        tr_z_xyz_xyzzz[i] = tr_z_yz_yzzz[i] * fe_0 + tr_z_yz_xyzzz[i] * pa_x[i];

        tr_z_xyz_xzzzz[i] = tr_z_xz_xzzzz[i] * pa_y[i];

        tr_z_xyz_yyyyy[i] = tr_z_yz_yyyyy[i] * pa_x[i];

        tr_z_xyz_yyyyz[i] = tr_z_yz_yyyyz[i] * pa_x[i];

        tr_z_xyz_yyyzz[i] = tr_z_yz_yyyzz[i] * pa_x[i];

        tr_z_xyz_yyzzz[i] = tr_z_yz_yyzzz[i] * pa_x[i];

        tr_z_xyz_yzzzz[i] = tr_z_yz_yzzzz[i] * pa_x[i];

        tr_z_xyz_zzzzz[i] = tr_z_yz_zzzzz[i] * pa_x[i];
    }

    // Set up 525-546 components of targeted buffer : FH

    auto tr_z_xzz_xxxxx = pbuffer.data(idx_dip_fh + 525);

    auto tr_z_xzz_xxxxy = pbuffer.data(idx_dip_fh + 526);

    auto tr_z_xzz_xxxxz = pbuffer.data(idx_dip_fh + 527);

    auto tr_z_xzz_xxxyy = pbuffer.data(idx_dip_fh + 528);

    auto tr_z_xzz_xxxyz = pbuffer.data(idx_dip_fh + 529);

    auto tr_z_xzz_xxxzz = pbuffer.data(idx_dip_fh + 530);

    auto tr_z_xzz_xxyyy = pbuffer.data(idx_dip_fh + 531);

    auto tr_z_xzz_xxyyz = pbuffer.data(idx_dip_fh + 532);

    auto tr_z_xzz_xxyzz = pbuffer.data(idx_dip_fh + 533);

    auto tr_z_xzz_xxzzz = pbuffer.data(idx_dip_fh + 534);

    auto tr_z_xzz_xyyyy = pbuffer.data(idx_dip_fh + 535);

    auto tr_z_xzz_xyyyz = pbuffer.data(idx_dip_fh + 536);

    auto tr_z_xzz_xyyzz = pbuffer.data(idx_dip_fh + 537);

    auto tr_z_xzz_xyzzz = pbuffer.data(idx_dip_fh + 538);

    auto tr_z_xzz_xzzzz = pbuffer.data(idx_dip_fh + 539);

    auto tr_z_xzz_yyyyy = pbuffer.data(idx_dip_fh + 540);

    auto tr_z_xzz_yyyyz = pbuffer.data(idx_dip_fh + 541);

    auto tr_z_xzz_yyyzz = pbuffer.data(idx_dip_fh + 542);

    auto tr_z_xzz_yyzzz = pbuffer.data(idx_dip_fh + 543);

    auto tr_z_xzz_yzzzz = pbuffer.data(idx_dip_fh + 544);

    auto tr_z_xzz_zzzzz = pbuffer.data(idx_dip_fh + 545);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xzz_xxxxx, \
                             tr_z_xzz_xxxxy, \
                             tr_z_xzz_xxxxz, \
                             tr_z_xzz_xxxyy, \
                             tr_z_xzz_xxxyz, \
                             tr_z_xzz_xxxzz, \
                             tr_z_xzz_xxyyy, \
                             tr_z_xzz_xxyyz, \
                             tr_z_xzz_xxyzz, \
                             tr_z_xzz_xxzzz, \
                             tr_z_xzz_xyyyy, \
                             tr_z_xzz_xyyyz, \
                             tr_z_xzz_xyyzz, \
                             tr_z_xzz_xyzzz, \
                             tr_z_xzz_xzzzz, \
                             tr_z_xzz_yyyyy, \
                             tr_z_xzz_yyyyz, \
                             tr_z_xzz_yyyzz, \
                             tr_z_xzz_yyzzz, \
                             tr_z_xzz_yzzzz, \
                             tr_z_xzz_zzzzz, \
                             tr_z_zz_xxxx,   \
                             tr_z_zz_xxxxx,  \
                             tr_z_zz_xxxxy,  \
                             tr_z_zz_xxxxz,  \
                             tr_z_zz_xxxy,   \
                             tr_z_zz_xxxyy,  \
                             tr_z_zz_xxxyz,  \
                             tr_z_zz_xxxz,   \
                             tr_z_zz_xxxzz,  \
                             tr_z_zz_xxyy,   \
                             tr_z_zz_xxyyy,  \
                             tr_z_zz_xxyyz,  \
                             tr_z_zz_xxyz,   \
                             tr_z_zz_xxyzz,  \
                             tr_z_zz_xxzz,   \
                             tr_z_zz_xxzzz,  \
                             tr_z_zz_xyyy,   \
                             tr_z_zz_xyyyy,  \
                             tr_z_zz_xyyyz,  \
                             tr_z_zz_xyyz,   \
                             tr_z_zz_xyyzz,  \
                             tr_z_zz_xyzz,   \
                             tr_z_zz_xyzzz,  \
                             tr_z_zz_xzzz,   \
                             tr_z_zz_xzzzz,  \
                             tr_z_zz_yyyy,   \
                             tr_z_zz_yyyyy,  \
                             tr_z_zz_yyyyz,  \
                             tr_z_zz_yyyz,   \
                             tr_z_zz_yyyzz,  \
                             tr_z_zz_yyzz,   \
                             tr_z_zz_yyzzz,  \
                             tr_z_zz_yzzz,   \
                             tr_z_zz_yzzzz,  \
                             tr_z_zz_zzzz,   \
                             tr_z_zz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xzz_xxxxx[i] = 5.0 * tr_z_zz_xxxx[i] * fe_0 + tr_z_zz_xxxxx[i] * pa_x[i];

        tr_z_xzz_xxxxy[i] = 4.0 * tr_z_zz_xxxy[i] * fe_0 + tr_z_zz_xxxxy[i] * pa_x[i];

        tr_z_xzz_xxxxz[i] = 4.0 * tr_z_zz_xxxz[i] * fe_0 + tr_z_zz_xxxxz[i] * pa_x[i];

        tr_z_xzz_xxxyy[i] = 3.0 * tr_z_zz_xxyy[i] * fe_0 + tr_z_zz_xxxyy[i] * pa_x[i];

        tr_z_xzz_xxxyz[i] = 3.0 * tr_z_zz_xxyz[i] * fe_0 + tr_z_zz_xxxyz[i] * pa_x[i];

        tr_z_xzz_xxxzz[i] = 3.0 * tr_z_zz_xxzz[i] * fe_0 + tr_z_zz_xxxzz[i] * pa_x[i];

        tr_z_xzz_xxyyy[i] = 2.0 * tr_z_zz_xyyy[i] * fe_0 + tr_z_zz_xxyyy[i] * pa_x[i];

        tr_z_xzz_xxyyz[i] = 2.0 * tr_z_zz_xyyz[i] * fe_0 + tr_z_zz_xxyyz[i] * pa_x[i];

        tr_z_xzz_xxyzz[i] = 2.0 * tr_z_zz_xyzz[i] * fe_0 + tr_z_zz_xxyzz[i] * pa_x[i];

        tr_z_xzz_xxzzz[i] = 2.0 * tr_z_zz_xzzz[i] * fe_0 + tr_z_zz_xxzzz[i] * pa_x[i];

        tr_z_xzz_xyyyy[i] = tr_z_zz_yyyy[i] * fe_0 + tr_z_zz_xyyyy[i] * pa_x[i];

        tr_z_xzz_xyyyz[i] = tr_z_zz_yyyz[i] * fe_0 + tr_z_zz_xyyyz[i] * pa_x[i];

        tr_z_xzz_xyyzz[i] = tr_z_zz_yyzz[i] * fe_0 + tr_z_zz_xyyzz[i] * pa_x[i];

        tr_z_xzz_xyzzz[i] = tr_z_zz_yzzz[i] * fe_0 + tr_z_zz_xyzzz[i] * pa_x[i];

        tr_z_xzz_xzzzz[i] = tr_z_zz_zzzz[i] * fe_0 + tr_z_zz_xzzzz[i] * pa_x[i];

        tr_z_xzz_yyyyy[i] = tr_z_zz_yyyyy[i] * pa_x[i];

        tr_z_xzz_yyyyz[i] = tr_z_zz_yyyyz[i] * pa_x[i];

        tr_z_xzz_yyyzz[i] = tr_z_zz_yyyzz[i] * pa_x[i];

        tr_z_xzz_yyzzz[i] = tr_z_zz_yyzzz[i] * pa_x[i];

        tr_z_xzz_yzzzz[i] = tr_z_zz_yzzzz[i] * pa_x[i];

        tr_z_xzz_zzzzz[i] = tr_z_zz_zzzzz[i] * pa_x[i];
    }

    // Set up 546-567 components of targeted buffer : FH

    auto tr_z_yyy_xxxxx = pbuffer.data(idx_dip_fh + 546);

    auto tr_z_yyy_xxxxy = pbuffer.data(idx_dip_fh + 547);

    auto tr_z_yyy_xxxxz = pbuffer.data(idx_dip_fh + 548);

    auto tr_z_yyy_xxxyy = pbuffer.data(idx_dip_fh + 549);

    auto tr_z_yyy_xxxyz = pbuffer.data(idx_dip_fh + 550);

    auto tr_z_yyy_xxxzz = pbuffer.data(idx_dip_fh + 551);

    auto tr_z_yyy_xxyyy = pbuffer.data(idx_dip_fh + 552);

    auto tr_z_yyy_xxyyz = pbuffer.data(idx_dip_fh + 553);

    auto tr_z_yyy_xxyzz = pbuffer.data(idx_dip_fh + 554);

    auto tr_z_yyy_xxzzz = pbuffer.data(idx_dip_fh + 555);

    auto tr_z_yyy_xyyyy = pbuffer.data(idx_dip_fh + 556);

    auto tr_z_yyy_xyyyz = pbuffer.data(idx_dip_fh + 557);

    auto tr_z_yyy_xyyzz = pbuffer.data(idx_dip_fh + 558);

    auto tr_z_yyy_xyzzz = pbuffer.data(idx_dip_fh + 559);

    auto tr_z_yyy_xzzzz = pbuffer.data(idx_dip_fh + 560);

    auto tr_z_yyy_yyyyy = pbuffer.data(idx_dip_fh + 561);

    auto tr_z_yyy_yyyyz = pbuffer.data(idx_dip_fh + 562);

    auto tr_z_yyy_yyyzz = pbuffer.data(idx_dip_fh + 563);

    auto tr_z_yyy_yyzzz = pbuffer.data(idx_dip_fh + 564);

    auto tr_z_yyy_yzzzz = pbuffer.data(idx_dip_fh + 565);

    auto tr_z_yyy_zzzzz = pbuffer.data(idx_dip_fh + 566);

#pragma omp simd aligned(pa_y,               \
                             tr_z_y_xxxxx,   \
                             tr_z_y_xxxxy,   \
                             tr_z_y_xxxxz,   \
                             tr_z_y_xxxyy,   \
                             tr_z_y_xxxyz,   \
                             tr_z_y_xxxzz,   \
                             tr_z_y_xxyyy,   \
                             tr_z_y_xxyyz,   \
                             tr_z_y_xxyzz,   \
                             tr_z_y_xxzzz,   \
                             tr_z_y_xyyyy,   \
                             tr_z_y_xyyyz,   \
                             tr_z_y_xyyzz,   \
                             tr_z_y_xyzzz,   \
                             tr_z_y_xzzzz,   \
                             tr_z_y_yyyyy,   \
                             tr_z_y_yyyyz,   \
                             tr_z_y_yyyzz,   \
                             tr_z_y_yyzzz,   \
                             tr_z_y_yzzzz,   \
                             tr_z_y_zzzzz,   \
                             tr_z_yy_xxxx,   \
                             tr_z_yy_xxxxx,  \
                             tr_z_yy_xxxxy,  \
                             tr_z_yy_xxxxz,  \
                             tr_z_yy_xxxy,   \
                             tr_z_yy_xxxyy,  \
                             tr_z_yy_xxxyz,  \
                             tr_z_yy_xxxz,   \
                             tr_z_yy_xxxzz,  \
                             tr_z_yy_xxyy,   \
                             tr_z_yy_xxyyy,  \
                             tr_z_yy_xxyyz,  \
                             tr_z_yy_xxyz,   \
                             tr_z_yy_xxyzz,  \
                             tr_z_yy_xxzz,   \
                             tr_z_yy_xxzzz,  \
                             tr_z_yy_xyyy,   \
                             tr_z_yy_xyyyy,  \
                             tr_z_yy_xyyyz,  \
                             tr_z_yy_xyyz,   \
                             tr_z_yy_xyyzz,  \
                             tr_z_yy_xyzz,   \
                             tr_z_yy_xyzzz,  \
                             tr_z_yy_xzzz,   \
                             tr_z_yy_xzzzz,  \
                             tr_z_yy_yyyy,   \
                             tr_z_yy_yyyyy,  \
                             tr_z_yy_yyyyz,  \
                             tr_z_yy_yyyz,   \
                             tr_z_yy_yyyzz,  \
                             tr_z_yy_yyzz,   \
                             tr_z_yy_yyzzz,  \
                             tr_z_yy_yzzz,   \
                             tr_z_yy_yzzzz,  \
                             tr_z_yy_zzzz,   \
                             tr_z_yy_zzzzz,  \
                             tr_z_yyy_xxxxx, \
                             tr_z_yyy_xxxxy, \
                             tr_z_yyy_xxxxz, \
                             tr_z_yyy_xxxyy, \
                             tr_z_yyy_xxxyz, \
                             tr_z_yyy_xxxzz, \
                             tr_z_yyy_xxyyy, \
                             tr_z_yyy_xxyyz, \
                             tr_z_yyy_xxyzz, \
                             tr_z_yyy_xxzzz, \
                             tr_z_yyy_xyyyy, \
                             tr_z_yyy_xyyyz, \
                             tr_z_yyy_xyyzz, \
                             tr_z_yyy_xyzzz, \
                             tr_z_yyy_xzzzz, \
                             tr_z_yyy_yyyyy, \
                             tr_z_yyy_yyyyz, \
                             tr_z_yyy_yyyzz, \
                             tr_z_yyy_yyzzz, \
                             tr_z_yyy_yzzzz, \
                             tr_z_yyy_zzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyy_xxxxx[i] = 2.0 * tr_z_y_xxxxx[i] * fe_0 + tr_z_yy_xxxxx[i] * pa_y[i];

        tr_z_yyy_xxxxy[i] = 2.0 * tr_z_y_xxxxy[i] * fe_0 + tr_z_yy_xxxx[i] * fe_0 + tr_z_yy_xxxxy[i] * pa_y[i];

        tr_z_yyy_xxxxz[i] = 2.0 * tr_z_y_xxxxz[i] * fe_0 + tr_z_yy_xxxxz[i] * pa_y[i];

        tr_z_yyy_xxxyy[i] = 2.0 * tr_z_y_xxxyy[i] * fe_0 + 2.0 * tr_z_yy_xxxy[i] * fe_0 + tr_z_yy_xxxyy[i] * pa_y[i];

        tr_z_yyy_xxxyz[i] = 2.0 * tr_z_y_xxxyz[i] * fe_0 + tr_z_yy_xxxz[i] * fe_0 + tr_z_yy_xxxyz[i] * pa_y[i];

        tr_z_yyy_xxxzz[i] = 2.0 * tr_z_y_xxxzz[i] * fe_0 + tr_z_yy_xxxzz[i] * pa_y[i];

        tr_z_yyy_xxyyy[i] = 2.0 * tr_z_y_xxyyy[i] * fe_0 + 3.0 * tr_z_yy_xxyy[i] * fe_0 + tr_z_yy_xxyyy[i] * pa_y[i];

        tr_z_yyy_xxyyz[i] = 2.0 * tr_z_y_xxyyz[i] * fe_0 + 2.0 * tr_z_yy_xxyz[i] * fe_0 + tr_z_yy_xxyyz[i] * pa_y[i];

        tr_z_yyy_xxyzz[i] = 2.0 * tr_z_y_xxyzz[i] * fe_0 + tr_z_yy_xxzz[i] * fe_0 + tr_z_yy_xxyzz[i] * pa_y[i];

        tr_z_yyy_xxzzz[i] = 2.0 * tr_z_y_xxzzz[i] * fe_0 + tr_z_yy_xxzzz[i] * pa_y[i];

        tr_z_yyy_xyyyy[i] = 2.0 * tr_z_y_xyyyy[i] * fe_0 + 4.0 * tr_z_yy_xyyy[i] * fe_0 + tr_z_yy_xyyyy[i] * pa_y[i];

        tr_z_yyy_xyyyz[i] = 2.0 * tr_z_y_xyyyz[i] * fe_0 + 3.0 * tr_z_yy_xyyz[i] * fe_0 + tr_z_yy_xyyyz[i] * pa_y[i];

        tr_z_yyy_xyyzz[i] = 2.0 * tr_z_y_xyyzz[i] * fe_0 + 2.0 * tr_z_yy_xyzz[i] * fe_0 + tr_z_yy_xyyzz[i] * pa_y[i];

        tr_z_yyy_xyzzz[i] = 2.0 * tr_z_y_xyzzz[i] * fe_0 + tr_z_yy_xzzz[i] * fe_0 + tr_z_yy_xyzzz[i] * pa_y[i];

        tr_z_yyy_xzzzz[i] = 2.0 * tr_z_y_xzzzz[i] * fe_0 + tr_z_yy_xzzzz[i] * pa_y[i];

        tr_z_yyy_yyyyy[i] = 2.0 * tr_z_y_yyyyy[i] * fe_0 + 5.0 * tr_z_yy_yyyy[i] * fe_0 + tr_z_yy_yyyyy[i] * pa_y[i];

        tr_z_yyy_yyyyz[i] = 2.0 * tr_z_y_yyyyz[i] * fe_0 + 4.0 * tr_z_yy_yyyz[i] * fe_0 + tr_z_yy_yyyyz[i] * pa_y[i];

        tr_z_yyy_yyyzz[i] = 2.0 * tr_z_y_yyyzz[i] * fe_0 + 3.0 * tr_z_yy_yyzz[i] * fe_0 + tr_z_yy_yyyzz[i] * pa_y[i];

        tr_z_yyy_yyzzz[i] = 2.0 * tr_z_y_yyzzz[i] * fe_0 + 2.0 * tr_z_yy_yzzz[i] * fe_0 + tr_z_yy_yyzzz[i] * pa_y[i];

        tr_z_yyy_yzzzz[i] = 2.0 * tr_z_y_yzzzz[i] * fe_0 + tr_z_yy_zzzz[i] * fe_0 + tr_z_yy_yzzzz[i] * pa_y[i];

        tr_z_yyy_zzzzz[i] = 2.0 * tr_z_y_zzzzz[i] * fe_0 + tr_z_yy_zzzzz[i] * pa_y[i];
    }

    // Set up 567-588 components of targeted buffer : FH

    auto tr_z_yyz_xxxxx = pbuffer.data(idx_dip_fh + 567);

    auto tr_z_yyz_xxxxy = pbuffer.data(idx_dip_fh + 568);

    auto tr_z_yyz_xxxxz = pbuffer.data(idx_dip_fh + 569);

    auto tr_z_yyz_xxxyy = pbuffer.data(idx_dip_fh + 570);

    auto tr_z_yyz_xxxyz = pbuffer.data(idx_dip_fh + 571);

    auto tr_z_yyz_xxxzz = pbuffer.data(idx_dip_fh + 572);

    auto tr_z_yyz_xxyyy = pbuffer.data(idx_dip_fh + 573);

    auto tr_z_yyz_xxyyz = pbuffer.data(idx_dip_fh + 574);

    auto tr_z_yyz_xxyzz = pbuffer.data(idx_dip_fh + 575);

    auto tr_z_yyz_xxzzz = pbuffer.data(idx_dip_fh + 576);

    auto tr_z_yyz_xyyyy = pbuffer.data(idx_dip_fh + 577);

    auto tr_z_yyz_xyyyz = pbuffer.data(idx_dip_fh + 578);

    auto tr_z_yyz_xyyzz = pbuffer.data(idx_dip_fh + 579);

    auto tr_z_yyz_xyzzz = pbuffer.data(idx_dip_fh + 580);

    auto tr_z_yyz_xzzzz = pbuffer.data(idx_dip_fh + 581);

    auto tr_z_yyz_yyyyy = pbuffer.data(idx_dip_fh + 582);

    auto tr_z_yyz_yyyyz = pbuffer.data(idx_dip_fh + 583);

    auto tr_z_yyz_yyyzz = pbuffer.data(idx_dip_fh + 584);

    auto tr_z_yyz_yyzzz = pbuffer.data(idx_dip_fh + 585);

    auto tr_z_yyz_yzzzz = pbuffer.data(idx_dip_fh + 586);

    auto tr_z_yyz_zzzzz = pbuffer.data(idx_dip_fh + 587);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_z_yy_xxxxy,  \
                             tr_z_yy_xxxyy,  \
                             tr_z_yy_xxyyy,  \
                             tr_z_yy_xyyyy,  \
                             tr_z_yy_yyyyy,  \
                             tr_z_yyz_xxxxx, \
                             tr_z_yyz_xxxxy, \
                             tr_z_yyz_xxxxz, \
                             tr_z_yyz_xxxyy, \
                             tr_z_yyz_xxxyz, \
                             tr_z_yyz_xxxzz, \
                             tr_z_yyz_xxyyy, \
                             tr_z_yyz_xxyyz, \
                             tr_z_yyz_xxyzz, \
                             tr_z_yyz_xxzzz, \
                             tr_z_yyz_xyyyy, \
                             tr_z_yyz_xyyyz, \
                             tr_z_yyz_xyyzz, \
                             tr_z_yyz_xyzzz, \
                             tr_z_yyz_xzzzz, \
                             tr_z_yyz_yyyyy, \
                             tr_z_yyz_yyyyz, \
                             tr_z_yyz_yyyzz, \
                             tr_z_yyz_yyzzz, \
                             tr_z_yyz_yzzzz, \
                             tr_z_yyz_zzzzz, \
                             tr_z_yz_xxxxx,  \
                             tr_z_yz_xxxxz,  \
                             tr_z_yz_xxxyz,  \
                             tr_z_yz_xxxz,   \
                             tr_z_yz_xxxzz,  \
                             tr_z_yz_xxyyz,  \
                             tr_z_yz_xxyz,   \
                             tr_z_yz_xxyzz,  \
                             tr_z_yz_xxzz,   \
                             tr_z_yz_xxzzz,  \
                             tr_z_yz_xyyyz,  \
                             tr_z_yz_xyyz,   \
                             tr_z_yz_xyyzz,  \
                             tr_z_yz_xyzz,   \
                             tr_z_yz_xyzzz,  \
                             tr_z_yz_xzzz,   \
                             tr_z_yz_xzzzz,  \
                             tr_z_yz_yyyyz,  \
                             tr_z_yz_yyyz,   \
                             tr_z_yz_yyyzz,  \
                             tr_z_yz_yyzz,   \
                             tr_z_yz_yyzzz,  \
                             tr_z_yz_yzzz,   \
                             tr_z_yz_yzzzz,  \
                             tr_z_yz_zzzz,   \
                             tr_z_yz_zzzzz,  \
                             tr_z_z_xxxxx,   \
                             tr_z_z_xxxxz,   \
                             tr_z_z_xxxyz,   \
                             tr_z_z_xxxzz,   \
                             tr_z_z_xxyyz,   \
                             tr_z_z_xxyzz,   \
                             tr_z_z_xxzzz,   \
                             tr_z_z_xyyyz,   \
                             tr_z_z_xyyzz,   \
                             tr_z_z_xyzzz,   \
                             tr_z_z_xzzzz,   \
                             tr_z_z_yyyyz,   \
                             tr_z_z_yyyzz,   \
                             tr_z_z_yyzzz,   \
                             tr_z_z_yzzzz,   \
                             tr_z_z_zzzzz,   \
                             ts_yy_xxxxy,    \
                             ts_yy_xxxyy,    \
                             ts_yy_xxyyy,    \
                             ts_yy_xyyyy,    \
                             ts_yy_yyyyy,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yyz_xxxxx[i] = tr_z_z_xxxxx[i] * fe_0 + tr_z_yz_xxxxx[i] * pa_y[i];

        tr_z_yyz_xxxxy[i] = ts_yy_xxxxy[i] * fe_0 + tr_z_yy_xxxxy[i] * pa_z[i];

        tr_z_yyz_xxxxz[i] = tr_z_z_xxxxz[i] * fe_0 + tr_z_yz_xxxxz[i] * pa_y[i];

        tr_z_yyz_xxxyy[i] = ts_yy_xxxyy[i] * fe_0 + tr_z_yy_xxxyy[i] * pa_z[i];

        tr_z_yyz_xxxyz[i] = tr_z_z_xxxyz[i] * fe_0 + tr_z_yz_xxxz[i] * fe_0 + tr_z_yz_xxxyz[i] * pa_y[i];

        tr_z_yyz_xxxzz[i] = tr_z_z_xxxzz[i] * fe_0 + tr_z_yz_xxxzz[i] * pa_y[i];

        tr_z_yyz_xxyyy[i] = ts_yy_xxyyy[i] * fe_0 + tr_z_yy_xxyyy[i] * pa_z[i];

        tr_z_yyz_xxyyz[i] = tr_z_z_xxyyz[i] * fe_0 + 2.0 * tr_z_yz_xxyz[i] * fe_0 + tr_z_yz_xxyyz[i] * pa_y[i];

        tr_z_yyz_xxyzz[i] = tr_z_z_xxyzz[i] * fe_0 + tr_z_yz_xxzz[i] * fe_0 + tr_z_yz_xxyzz[i] * pa_y[i];

        tr_z_yyz_xxzzz[i] = tr_z_z_xxzzz[i] * fe_0 + tr_z_yz_xxzzz[i] * pa_y[i];

        tr_z_yyz_xyyyy[i] = ts_yy_xyyyy[i] * fe_0 + tr_z_yy_xyyyy[i] * pa_z[i];

        tr_z_yyz_xyyyz[i] = tr_z_z_xyyyz[i] * fe_0 + 3.0 * tr_z_yz_xyyz[i] * fe_0 + tr_z_yz_xyyyz[i] * pa_y[i];

        tr_z_yyz_xyyzz[i] = tr_z_z_xyyzz[i] * fe_0 + 2.0 * tr_z_yz_xyzz[i] * fe_0 + tr_z_yz_xyyzz[i] * pa_y[i];

        tr_z_yyz_xyzzz[i] = tr_z_z_xyzzz[i] * fe_0 + tr_z_yz_xzzz[i] * fe_0 + tr_z_yz_xyzzz[i] * pa_y[i];

        tr_z_yyz_xzzzz[i] = tr_z_z_xzzzz[i] * fe_0 + tr_z_yz_xzzzz[i] * pa_y[i];

        tr_z_yyz_yyyyy[i] = ts_yy_yyyyy[i] * fe_0 + tr_z_yy_yyyyy[i] * pa_z[i];

        tr_z_yyz_yyyyz[i] = tr_z_z_yyyyz[i] * fe_0 + 4.0 * tr_z_yz_yyyz[i] * fe_0 + tr_z_yz_yyyyz[i] * pa_y[i];

        tr_z_yyz_yyyzz[i] = tr_z_z_yyyzz[i] * fe_0 + 3.0 * tr_z_yz_yyzz[i] * fe_0 + tr_z_yz_yyyzz[i] * pa_y[i];

        tr_z_yyz_yyzzz[i] = tr_z_z_yyzzz[i] * fe_0 + 2.0 * tr_z_yz_yzzz[i] * fe_0 + tr_z_yz_yyzzz[i] * pa_y[i];

        tr_z_yyz_yzzzz[i] = tr_z_z_yzzzz[i] * fe_0 + tr_z_yz_zzzz[i] * fe_0 + tr_z_yz_yzzzz[i] * pa_y[i];

        tr_z_yyz_zzzzz[i] = tr_z_z_zzzzz[i] * fe_0 + tr_z_yz_zzzzz[i] * pa_y[i];
    }

    // Set up 588-609 components of targeted buffer : FH

    auto tr_z_yzz_xxxxx = pbuffer.data(idx_dip_fh + 588);

    auto tr_z_yzz_xxxxy = pbuffer.data(idx_dip_fh + 589);

    auto tr_z_yzz_xxxxz = pbuffer.data(idx_dip_fh + 590);

    auto tr_z_yzz_xxxyy = pbuffer.data(idx_dip_fh + 591);

    auto tr_z_yzz_xxxyz = pbuffer.data(idx_dip_fh + 592);

    auto tr_z_yzz_xxxzz = pbuffer.data(idx_dip_fh + 593);

    auto tr_z_yzz_xxyyy = pbuffer.data(idx_dip_fh + 594);

    auto tr_z_yzz_xxyyz = pbuffer.data(idx_dip_fh + 595);

    auto tr_z_yzz_xxyzz = pbuffer.data(idx_dip_fh + 596);

    auto tr_z_yzz_xxzzz = pbuffer.data(idx_dip_fh + 597);

    auto tr_z_yzz_xyyyy = pbuffer.data(idx_dip_fh + 598);

    auto tr_z_yzz_xyyyz = pbuffer.data(idx_dip_fh + 599);

    auto tr_z_yzz_xyyzz = pbuffer.data(idx_dip_fh + 600);

    auto tr_z_yzz_xyzzz = pbuffer.data(idx_dip_fh + 601);

    auto tr_z_yzz_xzzzz = pbuffer.data(idx_dip_fh + 602);

    auto tr_z_yzz_yyyyy = pbuffer.data(idx_dip_fh + 603);

    auto tr_z_yzz_yyyyz = pbuffer.data(idx_dip_fh + 604);

    auto tr_z_yzz_yyyzz = pbuffer.data(idx_dip_fh + 605);

    auto tr_z_yzz_yyzzz = pbuffer.data(idx_dip_fh + 606);

    auto tr_z_yzz_yzzzz = pbuffer.data(idx_dip_fh + 607);

    auto tr_z_yzz_zzzzz = pbuffer.data(idx_dip_fh + 608);

#pragma omp simd aligned(pa_y,               \
                             tr_z_yzz_xxxxx, \
                             tr_z_yzz_xxxxy, \
                             tr_z_yzz_xxxxz, \
                             tr_z_yzz_xxxyy, \
                             tr_z_yzz_xxxyz, \
                             tr_z_yzz_xxxzz, \
                             tr_z_yzz_xxyyy, \
                             tr_z_yzz_xxyyz, \
                             tr_z_yzz_xxyzz, \
                             tr_z_yzz_xxzzz, \
                             tr_z_yzz_xyyyy, \
                             tr_z_yzz_xyyyz, \
                             tr_z_yzz_xyyzz, \
                             tr_z_yzz_xyzzz, \
                             tr_z_yzz_xzzzz, \
                             tr_z_yzz_yyyyy, \
                             tr_z_yzz_yyyyz, \
                             tr_z_yzz_yyyzz, \
                             tr_z_yzz_yyzzz, \
                             tr_z_yzz_yzzzz, \
                             tr_z_yzz_zzzzz, \
                             tr_z_zz_xxxx,   \
                             tr_z_zz_xxxxx,  \
                             tr_z_zz_xxxxy,  \
                             tr_z_zz_xxxxz,  \
                             tr_z_zz_xxxy,   \
                             tr_z_zz_xxxyy,  \
                             tr_z_zz_xxxyz,  \
                             tr_z_zz_xxxz,   \
                             tr_z_zz_xxxzz,  \
                             tr_z_zz_xxyy,   \
                             tr_z_zz_xxyyy,  \
                             tr_z_zz_xxyyz,  \
                             tr_z_zz_xxyz,   \
                             tr_z_zz_xxyzz,  \
                             tr_z_zz_xxzz,   \
                             tr_z_zz_xxzzz,  \
                             tr_z_zz_xyyy,   \
                             tr_z_zz_xyyyy,  \
                             tr_z_zz_xyyyz,  \
                             tr_z_zz_xyyz,   \
                             tr_z_zz_xyyzz,  \
                             tr_z_zz_xyzz,   \
                             tr_z_zz_xyzzz,  \
                             tr_z_zz_xzzz,   \
                             tr_z_zz_xzzzz,  \
                             tr_z_zz_yyyy,   \
                             tr_z_zz_yyyyy,  \
                             tr_z_zz_yyyyz,  \
                             tr_z_zz_yyyz,   \
                             tr_z_zz_yyyzz,  \
                             tr_z_zz_yyzz,   \
                             tr_z_zz_yyzzz,  \
                             tr_z_zz_yzzz,   \
                             tr_z_zz_yzzzz,  \
                             tr_z_zz_zzzz,   \
                             tr_z_zz_zzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yzz_xxxxx[i] = tr_z_zz_xxxxx[i] * pa_y[i];

        tr_z_yzz_xxxxy[i] = tr_z_zz_xxxx[i] * fe_0 + tr_z_zz_xxxxy[i] * pa_y[i];

        tr_z_yzz_xxxxz[i] = tr_z_zz_xxxxz[i] * pa_y[i];

        tr_z_yzz_xxxyy[i] = 2.0 * tr_z_zz_xxxy[i] * fe_0 + tr_z_zz_xxxyy[i] * pa_y[i];

        tr_z_yzz_xxxyz[i] = tr_z_zz_xxxz[i] * fe_0 + tr_z_zz_xxxyz[i] * pa_y[i];

        tr_z_yzz_xxxzz[i] = tr_z_zz_xxxzz[i] * pa_y[i];

        tr_z_yzz_xxyyy[i] = 3.0 * tr_z_zz_xxyy[i] * fe_0 + tr_z_zz_xxyyy[i] * pa_y[i];

        tr_z_yzz_xxyyz[i] = 2.0 * tr_z_zz_xxyz[i] * fe_0 + tr_z_zz_xxyyz[i] * pa_y[i];

        tr_z_yzz_xxyzz[i] = tr_z_zz_xxzz[i] * fe_0 + tr_z_zz_xxyzz[i] * pa_y[i];

        tr_z_yzz_xxzzz[i] = tr_z_zz_xxzzz[i] * pa_y[i];

        tr_z_yzz_xyyyy[i] = 4.0 * tr_z_zz_xyyy[i] * fe_0 + tr_z_zz_xyyyy[i] * pa_y[i];

        tr_z_yzz_xyyyz[i] = 3.0 * tr_z_zz_xyyz[i] * fe_0 + tr_z_zz_xyyyz[i] * pa_y[i];

        tr_z_yzz_xyyzz[i] = 2.0 * tr_z_zz_xyzz[i] * fe_0 + tr_z_zz_xyyzz[i] * pa_y[i];

        tr_z_yzz_xyzzz[i] = tr_z_zz_xzzz[i] * fe_0 + tr_z_zz_xyzzz[i] * pa_y[i];

        tr_z_yzz_xzzzz[i] = tr_z_zz_xzzzz[i] * pa_y[i];

        tr_z_yzz_yyyyy[i] = 5.0 * tr_z_zz_yyyy[i] * fe_0 + tr_z_zz_yyyyy[i] * pa_y[i];

        tr_z_yzz_yyyyz[i] = 4.0 * tr_z_zz_yyyz[i] * fe_0 + tr_z_zz_yyyyz[i] * pa_y[i];

        tr_z_yzz_yyyzz[i] = 3.0 * tr_z_zz_yyzz[i] * fe_0 + tr_z_zz_yyyzz[i] * pa_y[i];

        tr_z_yzz_yyzzz[i] = 2.0 * tr_z_zz_yzzz[i] * fe_0 + tr_z_zz_yyzzz[i] * pa_y[i];

        tr_z_yzz_yzzzz[i] = tr_z_zz_zzzz[i] * fe_0 + tr_z_zz_yzzzz[i] * pa_y[i];

        tr_z_yzz_zzzzz[i] = tr_z_zz_zzzzz[i] * pa_y[i];
    }

    // Set up 609-630 components of targeted buffer : FH

    auto tr_z_zzz_xxxxx = pbuffer.data(idx_dip_fh + 609);

    auto tr_z_zzz_xxxxy = pbuffer.data(idx_dip_fh + 610);

    auto tr_z_zzz_xxxxz = pbuffer.data(idx_dip_fh + 611);

    auto tr_z_zzz_xxxyy = pbuffer.data(idx_dip_fh + 612);

    auto tr_z_zzz_xxxyz = pbuffer.data(idx_dip_fh + 613);

    auto tr_z_zzz_xxxzz = pbuffer.data(idx_dip_fh + 614);

    auto tr_z_zzz_xxyyy = pbuffer.data(idx_dip_fh + 615);

    auto tr_z_zzz_xxyyz = pbuffer.data(idx_dip_fh + 616);

    auto tr_z_zzz_xxyzz = pbuffer.data(idx_dip_fh + 617);

    auto tr_z_zzz_xxzzz = pbuffer.data(idx_dip_fh + 618);

    auto tr_z_zzz_xyyyy = pbuffer.data(idx_dip_fh + 619);

    auto tr_z_zzz_xyyyz = pbuffer.data(idx_dip_fh + 620);

    auto tr_z_zzz_xyyzz = pbuffer.data(idx_dip_fh + 621);

    auto tr_z_zzz_xyzzz = pbuffer.data(idx_dip_fh + 622);

    auto tr_z_zzz_xzzzz = pbuffer.data(idx_dip_fh + 623);

    auto tr_z_zzz_yyyyy = pbuffer.data(idx_dip_fh + 624);

    auto tr_z_zzz_yyyyz = pbuffer.data(idx_dip_fh + 625);

    auto tr_z_zzz_yyyzz = pbuffer.data(idx_dip_fh + 626);

    auto tr_z_zzz_yyzzz = pbuffer.data(idx_dip_fh + 627);

    auto tr_z_zzz_yzzzz = pbuffer.data(idx_dip_fh + 628);

    auto tr_z_zzz_zzzzz = pbuffer.data(idx_dip_fh + 629);

#pragma omp simd aligned(pa_z,               \
                             tr_z_z_xxxxx,   \
                             tr_z_z_xxxxy,   \
                             tr_z_z_xxxxz,   \
                             tr_z_z_xxxyy,   \
                             tr_z_z_xxxyz,   \
                             tr_z_z_xxxzz,   \
                             tr_z_z_xxyyy,   \
                             tr_z_z_xxyyz,   \
                             tr_z_z_xxyzz,   \
                             tr_z_z_xxzzz,   \
                             tr_z_z_xyyyy,   \
                             tr_z_z_xyyyz,   \
                             tr_z_z_xyyzz,   \
                             tr_z_z_xyzzz,   \
                             tr_z_z_xzzzz,   \
                             tr_z_z_yyyyy,   \
                             tr_z_z_yyyyz,   \
                             tr_z_z_yyyzz,   \
                             tr_z_z_yyzzz,   \
                             tr_z_z_yzzzz,   \
                             tr_z_z_zzzzz,   \
                             tr_z_zz_xxxx,   \
                             tr_z_zz_xxxxx,  \
                             tr_z_zz_xxxxy,  \
                             tr_z_zz_xxxxz,  \
                             tr_z_zz_xxxy,   \
                             tr_z_zz_xxxyy,  \
                             tr_z_zz_xxxyz,  \
                             tr_z_zz_xxxz,   \
                             tr_z_zz_xxxzz,  \
                             tr_z_zz_xxyy,   \
                             tr_z_zz_xxyyy,  \
                             tr_z_zz_xxyyz,  \
                             tr_z_zz_xxyz,   \
                             tr_z_zz_xxyzz,  \
                             tr_z_zz_xxzz,   \
                             tr_z_zz_xxzzz,  \
                             tr_z_zz_xyyy,   \
                             tr_z_zz_xyyyy,  \
                             tr_z_zz_xyyyz,  \
                             tr_z_zz_xyyz,   \
                             tr_z_zz_xyyzz,  \
                             tr_z_zz_xyzz,   \
                             tr_z_zz_xyzzz,  \
                             tr_z_zz_xzzz,   \
                             tr_z_zz_xzzzz,  \
                             tr_z_zz_yyyy,   \
                             tr_z_zz_yyyyy,  \
                             tr_z_zz_yyyyz,  \
                             tr_z_zz_yyyz,   \
                             tr_z_zz_yyyzz,  \
                             tr_z_zz_yyzz,   \
                             tr_z_zz_yyzzz,  \
                             tr_z_zz_yzzz,   \
                             tr_z_zz_yzzzz,  \
                             tr_z_zz_zzzz,   \
                             tr_z_zz_zzzzz,  \
                             tr_z_zzz_xxxxx, \
                             tr_z_zzz_xxxxy, \
                             tr_z_zzz_xxxxz, \
                             tr_z_zzz_xxxyy, \
                             tr_z_zzz_xxxyz, \
                             tr_z_zzz_xxxzz, \
                             tr_z_zzz_xxyyy, \
                             tr_z_zzz_xxyyz, \
                             tr_z_zzz_xxyzz, \
                             tr_z_zzz_xxzzz, \
                             tr_z_zzz_xyyyy, \
                             tr_z_zzz_xyyyz, \
                             tr_z_zzz_xyyzz, \
                             tr_z_zzz_xyzzz, \
                             tr_z_zzz_xzzzz, \
                             tr_z_zzz_yyyyy, \
                             tr_z_zzz_yyyyz, \
                             tr_z_zzz_yyyzz, \
                             tr_z_zzz_yyzzz, \
                             tr_z_zzz_yzzzz, \
                             tr_z_zzz_zzzzz, \
                             ts_zz_xxxxx,    \
                             ts_zz_xxxxy,    \
                             ts_zz_xxxxz,    \
                             ts_zz_xxxyy,    \
                             ts_zz_xxxyz,    \
                             ts_zz_xxxzz,    \
                             ts_zz_xxyyy,    \
                             ts_zz_xxyyz,    \
                             ts_zz_xxyzz,    \
                             ts_zz_xxzzz,    \
                             ts_zz_xyyyy,    \
                             ts_zz_xyyyz,    \
                             ts_zz_xyyzz,    \
                             ts_zz_xyzzz,    \
                             ts_zz_xzzzz,    \
                             ts_zz_yyyyy,    \
                             ts_zz_yyyyz,    \
                             ts_zz_yyyzz,    \
                             ts_zz_yyzzz,    \
                             ts_zz_yzzzz,    \
                             ts_zz_zzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zzz_xxxxx[i] = 2.0 * tr_z_z_xxxxx[i] * fe_0 + ts_zz_xxxxx[i] * fe_0 + tr_z_zz_xxxxx[i] * pa_z[i];

        tr_z_zzz_xxxxy[i] = 2.0 * tr_z_z_xxxxy[i] * fe_0 + ts_zz_xxxxy[i] * fe_0 + tr_z_zz_xxxxy[i] * pa_z[i];

        tr_z_zzz_xxxxz[i] = 2.0 * tr_z_z_xxxxz[i] * fe_0 + tr_z_zz_xxxx[i] * fe_0 + ts_zz_xxxxz[i] * fe_0 + tr_z_zz_xxxxz[i] * pa_z[i];

        tr_z_zzz_xxxyy[i] = 2.0 * tr_z_z_xxxyy[i] * fe_0 + ts_zz_xxxyy[i] * fe_0 + tr_z_zz_xxxyy[i] * pa_z[i];

        tr_z_zzz_xxxyz[i] = 2.0 * tr_z_z_xxxyz[i] * fe_0 + tr_z_zz_xxxy[i] * fe_0 + ts_zz_xxxyz[i] * fe_0 + tr_z_zz_xxxyz[i] * pa_z[i];

        tr_z_zzz_xxxzz[i] = 2.0 * tr_z_z_xxxzz[i] * fe_0 + 2.0 * tr_z_zz_xxxz[i] * fe_0 + ts_zz_xxxzz[i] * fe_0 + tr_z_zz_xxxzz[i] * pa_z[i];

        tr_z_zzz_xxyyy[i] = 2.0 * tr_z_z_xxyyy[i] * fe_0 + ts_zz_xxyyy[i] * fe_0 + tr_z_zz_xxyyy[i] * pa_z[i];

        tr_z_zzz_xxyyz[i] = 2.0 * tr_z_z_xxyyz[i] * fe_0 + tr_z_zz_xxyy[i] * fe_0 + ts_zz_xxyyz[i] * fe_0 + tr_z_zz_xxyyz[i] * pa_z[i];

        tr_z_zzz_xxyzz[i] = 2.0 * tr_z_z_xxyzz[i] * fe_0 + 2.0 * tr_z_zz_xxyz[i] * fe_0 + ts_zz_xxyzz[i] * fe_0 + tr_z_zz_xxyzz[i] * pa_z[i];

        tr_z_zzz_xxzzz[i] = 2.0 * tr_z_z_xxzzz[i] * fe_0 + 3.0 * tr_z_zz_xxzz[i] * fe_0 + ts_zz_xxzzz[i] * fe_0 + tr_z_zz_xxzzz[i] * pa_z[i];

        tr_z_zzz_xyyyy[i] = 2.0 * tr_z_z_xyyyy[i] * fe_0 + ts_zz_xyyyy[i] * fe_0 + tr_z_zz_xyyyy[i] * pa_z[i];

        tr_z_zzz_xyyyz[i] = 2.0 * tr_z_z_xyyyz[i] * fe_0 + tr_z_zz_xyyy[i] * fe_0 + ts_zz_xyyyz[i] * fe_0 + tr_z_zz_xyyyz[i] * pa_z[i];

        tr_z_zzz_xyyzz[i] = 2.0 * tr_z_z_xyyzz[i] * fe_0 + 2.0 * tr_z_zz_xyyz[i] * fe_0 + ts_zz_xyyzz[i] * fe_0 + tr_z_zz_xyyzz[i] * pa_z[i];

        tr_z_zzz_xyzzz[i] = 2.0 * tr_z_z_xyzzz[i] * fe_0 + 3.0 * tr_z_zz_xyzz[i] * fe_0 + ts_zz_xyzzz[i] * fe_0 + tr_z_zz_xyzzz[i] * pa_z[i];

        tr_z_zzz_xzzzz[i] = 2.0 * tr_z_z_xzzzz[i] * fe_0 + 4.0 * tr_z_zz_xzzz[i] * fe_0 + ts_zz_xzzzz[i] * fe_0 + tr_z_zz_xzzzz[i] * pa_z[i];

        tr_z_zzz_yyyyy[i] = 2.0 * tr_z_z_yyyyy[i] * fe_0 + ts_zz_yyyyy[i] * fe_0 + tr_z_zz_yyyyy[i] * pa_z[i];

        tr_z_zzz_yyyyz[i] = 2.0 * tr_z_z_yyyyz[i] * fe_0 + tr_z_zz_yyyy[i] * fe_0 + ts_zz_yyyyz[i] * fe_0 + tr_z_zz_yyyyz[i] * pa_z[i];

        tr_z_zzz_yyyzz[i] = 2.0 * tr_z_z_yyyzz[i] * fe_0 + 2.0 * tr_z_zz_yyyz[i] * fe_0 + ts_zz_yyyzz[i] * fe_0 + tr_z_zz_yyyzz[i] * pa_z[i];

        tr_z_zzz_yyzzz[i] = 2.0 * tr_z_z_yyzzz[i] * fe_0 + 3.0 * tr_z_zz_yyzz[i] * fe_0 + ts_zz_yyzzz[i] * fe_0 + tr_z_zz_yyzzz[i] * pa_z[i];

        tr_z_zzz_yzzzz[i] = 2.0 * tr_z_z_yzzzz[i] * fe_0 + 4.0 * tr_z_zz_yzzz[i] * fe_0 + ts_zz_yzzzz[i] * fe_0 + tr_z_zz_yzzzz[i] * pa_z[i];

        tr_z_zzz_zzzzz[i] = 2.0 * tr_z_z_zzzzz[i] * fe_0 + 5.0 * tr_z_zz_zzzz[i] * fe_0 + ts_zz_zzzzz[i] * fe_0 + tr_z_zz_zzzzz[i] * pa_z[i];
    }
}

}  // namespace diprec
