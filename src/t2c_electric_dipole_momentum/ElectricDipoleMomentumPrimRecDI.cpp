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

#include "ElectricDipoleMomentumPrimRecDI.hpp"

namespace diprec {  // diprec namespace

auto
comp_prim_electric_dipole_momentum_di(CSimdArray<double>&       pbuffer,
                                      const size_t              idx_dip_di,
                                      const size_t              idx_dip_si,
                                      const size_t              idx_dip_ph,
                                      const size_t              idx_ovl_pi,
                                      const size_t              idx_dip_pi,
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

    // Set up components of auxiliary buffer : SI

    auto tr_x_0_xxxxxx = pbuffer.data(idx_dip_si);

    auto tr_x_0_xxxxxy = pbuffer.data(idx_dip_si + 1);

    auto tr_x_0_xxxxxz = pbuffer.data(idx_dip_si + 2);

    auto tr_x_0_xxxxyy = pbuffer.data(idx_dip_si + 3);

    auto tr_x_0_xxxxyz = pbuffer.data(idx_dip_si + 4);

    auto tr_x_0_xxxxzz = pbuffer.data(idx_dip_si + 5);

    auto tr_x_0_xxxyyy = pbuffer.data(idx_dip_si + 6);

    auto tr_x_0_xxxyyz = pbuffer.data(idx_dip_si + 7);

    auto tr_x_0_xxxyzz = pbuffer.data(idx_dip_si + 8);

    auto tr_x_0_xxxzzz = pbuffer.data(idx_dip_si + 9);

    auto tr_x_0_xxyyyy = pbuffer.data(idx_dip_si + 10);

    auto tr_x_0_xxyyyz = pbuffer.data(idx_dip_si + 11);

    auto tr_x_0_xxyyzz = pbuffer.data(idx_dip_si + 12);

    auto tr_x_0_xxyzzz = pbuffer.data(idx_dip_si + 13);

    auto tr_x_0_xxzzzz = pbuffer.data(idx_dip_si + 14);

    auto tr_x_0_xyyyyy = pbuffer.data(idx_dip_si + 15);

    auto tr_x_0_xyyyyz = pbuffer.data(idx_dip_si + 16);

    auto tr_x_0_xyyyzz = pbuffer.data(idx_dip_si + 17);

    auto tr_x_0_xyyzzz = pbuffer.data(idx_dip_si + 18);

    auto tr_x_0_xyzzzz = pbuffer.data(idx_dip_si + 19);

    auto tr_x_0_xzzzzz = pbuffer.data(idx_dip_si + 20);

    auto tr_x_0_yyyyyy = pbuffer.data(idx_dip_si + 21);

    auto tr_x_0_yyyyyz = pbuffer.data(idx_dip_si + 22);

    auto tr_x_0_yyyyzz = pbuffer.data(idx_dip_si + 23);

    auto tr_x_0_yyyzzz = pbuffer.data(idx_dip_si + 24);

    auto tr_x_0_yyzzzz = pbuffer.data(idx_dip_si + 25);

    auto tr_x_0_yzzzzz = pbuffer.data(idx_dip_si + 26);

    auto tr_x_0_zzzzzz = pbuffer.data(idx_dip_si + 27);

    auto tr_y_0_xxxxxx = pbuffer.data(idx_dip_si + 28);

    auto tr_y_0_xxxxxy = pbuffer.data(idx_dip_si + 29);

    auto tr_y_0_xxxxxz = pbuffer.data(idx_dip_si + 30);

    auto tr_y_0_xxxxyy = pbuffer.data(idx_dip_si + 31);

    auto tr_y_0_xxxxyz = pbuffer.data(idx_dip_si + 32);

    auto tr_y_0_xxxxzz = pbuffer.data(idx_dip_si + 33);

    auto tr_y_0_xxxyyy = pbuffer.data(idx_dip_si + 34);

    auto tr_y_0_xxxyyz = pbuffer.data(idx_dip_si + 35);

    auto tr_y_0_xxxyzz = pbuffer.data(idx_dip_si + 36);

    auto tr_y_0_xxxzzz = pbuffer.data(idx_dip_si + 37);

    auto tr_y_0_xxyyyy = pbuffer.data(idx_dip_si + 38);

    auto tr_y_0_xxyyyz = pbuffer.data(idx_dip_si + 39);

    auto tr_y_0_xxyyzz = pbuffer.data(idx_dip_si + 40);

    auto tr_y_0_xxyzzz = pbuffer.data(idx_dip_si + 41);

    auto tr_y_0_xxzzzz = pbuffer.data(idx_dip_si + 42);

    auto tr_y_0_xyyyyy = pbuffer.data(idx_dip_si + 43);

    auto tr_y_0_xyyyyz = pbuffer.data(idx_dip_si + 44);

    auto tr_y_0_xyyyzz = pbuffer.data(idx_dip_si + 45);

    auto tr_y_0_xyyzzz = pbuffer.data(idx_dip_si + 46);

    auto tr_y_0_xyzzzz = pbuffer.data(idx_dip_si + 47);

    auto tr_y_0_xzzzzz = pbuffer.data(idx_dip_si + 48);

    auto tr_y_0_yyyyyy = pbuffer.data(idx_dip_si + 49);

    auto tr_y_0_yyyyyz = pbuffer.data(idx_dip_si + 50);

    auto tr_y_0_yyyyzz = pbuffer.data(idx_dip_si + 51);

    auto tr_y_0_yyyzzz = pbuffer.data(idx_dip_si + 52);

    auto tr_y_0_yyzzzz = pbuffer.data(idx_dip_si + 53);

    auto tr_y_0_yzzzzz = pbuffer.data(idx_dip_si + 54);

    auto tr_y_0_zzzzzz = pbuffer.data(idx_dip_si + 55);

    auto tr_z_0_xxxxxx = pbuffer.data(idx_dip_si + 56);

    auto tr_z_0_xxxxxy = pbuffer.data(idx_dip_si + 57);

    auto tr_z_0_xxxxxz = pbuffer.data(idx_dip_si + 58);

    auto tr_z_0_xxxxyy = pbuffer.data(idx_dip_si + 59);

    auto tr_z_0_xxxxyz = pbuffer.data(idx_dip_si + 60);

    auto tr_z_0_xxxxzz = pbuffer.data(idx_dip_si + 61);

    auto tr_z_0_xxxyyy = pbuffer.data(idx_dip_si + 62);

    auto tr_z_0_xxxyyz = pbuffer.data(idx_dip_si + 63);

    auto tr_z_0_xxxyzz = pbuffer.data(idx_dip_si + 64);

    auto tr_z_0_xxxzzz = pbuffer.data(idx_dip_si + 65);

    auto tr_z_0_xxyyyy = pbuffer.data(idx_dip_si + 66);

    auto tr_z_0_xxyyyz = pbuffer.data(idx_dip_si + 67);

    auto tr_z_0_xxyyzz = pbuffer.data(idx_dip_si + 68);

    auto tr_z_0_xxyzzz = pbuffer.data(idx_dip_si + 69);

    auto tr_z_0_xxzzzz = pbuffer.data(idx_dip_si + 70);

    auto tr_z_0_xyyyyy = pbuffer.data(idx_dip_si + 71);

    auto tr_z_0_xyyyyz = pbuffer.data(idx_dip_si + 72);

    auto tr_z_0_xyyyzz = pbuffer.data(idx_dip_si + 73);

    auto tr_z_0_xyyzzz = pbuffer.data(idx_dip_si + 74);

    auto tr_z_0_xyzzzz = pbuffer.data(idx_dip_si + 75);

    auto tr_z_0_xzzzzz = pbuffer.data(idx_dip_si + 76);

    auto tr_z_0_yyyyyy = pbuffer.data(idx_dip_si + 77);

    auto tr_z_0_yyyyyz = pbuffer.data(idx_dip_si + 78);

    auto tr_z_0_yyyyzz = pbuffer.data(idx_dip_si + 79);

    auto tr_z_0_yyyzzz = pbuffer.data(idx_dip_si + 80);

    auto tr_z_0_yyzzzz = pbuffer.data(idx_dip_si + 81);

    auto tr_z_0_yzzzzz = pbuffer.data(idx_dip_si + 82);

    auto tr_z_0_zzzzzz = pbuffer.data(idx_dip_si + 83);

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

    // Set up components of auxiliary buffer : PI

    auto ts_x_xxxxxx = pbuffer.data(idx_ovl_pi);

    auto ts_x_xxxxxy = pbuffer.data(idx_ovl_pi + 1);

    auto ts_x_xxxxxz = pbuffer.data(idx_ovl_pi + 2);

    auto ts_x_xxxxyy = pbuffer.data(idx_ovl_pi + 3);

    auto ts_x_xxxxyz = pbuffer.data(idx_ovl_pi + 4);

    auto ts_x_xxxxzz = pbuffer.data(idx_ovl_pi + 5);

    auto ts_x_xxxyyy = pbuffer.data(idx_ovl_pi + 6);

    auto ts_x_xxxyyz = pbuffer.data(idx_ovl_pi + 7);

    auto ts_x_xxxyzz = pbuffer.data(idx_ovl_pi + 8);

    auto ts_x_xxxzzz = pbuffer.data(idx_ovl_pi + 9);

    auto ts_x_xxyyyy = pbuffer.data(idx_ovl_pi + 10);

    auto ts_x_xxyyyz = pbuffer.data(idx_ovl_pi + 11);

    auto ts_x_xxyyzz = pbuffer.data(idx_ovl_pi + 12);

    auto ts_x_xxyzzz = pbuffer.data(idx_ovl_pi + 13);

    auto ts_x_xxzzzz = pbuffer.data(idx_ovl_pi + 14);

    auto ts_x_xyyyyy = pbuffer.data(idx_ovl_pi + 15);

    auto ts_x_xyyyyz = pbuffer.data(idx_ovl_pi + 16);

    auto ts_x_xyyyzz = pbuffer.data(idx_ovl_pi + 17);

    auto ts_x_xyyzzz = pbuffer.data(idx_ovl_pi + 18);

    auto ts_x_xyzzzz = pbuffer.data(idx_ovl_pi + 19);

    auto ts_x_xzzzzz = pbuffer.data(idx_ovl_pi + 20);

    auto ts_x_yyyyyy = pbuffer.data(idx_ovl_pi + 21);

    auto ts_x_yyyyyz = pbuffer.data(idx_ovl_pi + 22);

    auto ts_x_yyyyzz = pbuffer.data(idx_ovl_pi + 23);

    auto ts_x_yyyzzz = pbuffer.data(idx_ovl_pi + 24);

    auto ts_x_yyzzzz = pbuffer.data(idx_ovl_pi + 25);

    auto ts_x_yzzzzz = pbuffer.data(idx_ovl_pi + 26);

    auto ts_x_zzzzzz = pbuffer.data(idx_ovl_pi + 27);

    auto ts_y_xxxxxx = pbuffer.data(idx_ovl_pi + 28);

    auto ts_y_xxxxxy = pbuffer.data(idx_ovl_pi + 29);

    auto ts_y_xxxxxz = pbuffer.data(idx_ovl_pi + 30);

    auto ts_y_xxxxyy = pbuffer.data(idx_ovl_pi + 31);

    auto ts_y_xxxxyz = pbuffer.data(idx_ovl_pi + 32);

    auto ts_y_xxxxzz = pbuffer.data(idx_ovl_pi + 33);

    auto ts_y_xxxyyy = pbuffer.data(idx_ovl_pi + 34);

    auto ts_y_xxxyyz = pbuffer.data(idx_ovl_pi + 35);

    auto ts_y_xxxyzz = pbuffer.data(idx_ovl_pi + 36);

    auto ts_y_xxxzzz = pbuffer.data(idx_ovl_pi + 37);

    auto ts_y_xxyyyy = pbuffer.data(idx_ovl_pi + 38);

    auto ts_y_xxyyyz = pbuffer.data(idx_ovl_pi + 39);

    auto ts_y_xxyyzz = pbuffer.data(idx_ovl_pi + 40);

    auto ts_y_xxyzzz = pbuffer.data(idx_ovl_pi + 41);

    auto ts_y_xxzzzz = pbuffer.data(idx_ovl_pi + 42);

    auto ts_y_xyyyyy = pbuffer.data(idx_ovl_pi + 43);

    auto ts_y_xyyyyz = pbuffer.data(idx_ovl_pi + 44);

    auto ts_y_xyyyzz = pbuffer.data(idx_ovl_pi + 45);

    auto ts_y_xyyzzz = pbuffer.data(idx_ovl_pi + 46);

    auto ts_y_xyzzzz = pbuffer.data(idx_ovl_pi + 47);

    auto ts_y_xzzzzz = pbuffer.data(idx_ovl_pi + 48);

    auto ts_y_yyyyyy = pbuffer.data(idx_ovl_pi + 49);

    auto ts_y_yyyyyz = pbuffer.data(idx_ovl_pi + 50);

    auto ts_y_yyyyzz = pbuffer.data(idx_ovl_pi + 51);

    auto ts_y_yyyzzz = pbuffer.data(idx_ovl_pi + 52);

    auto ts_y_yyzzzz = pbuffer.data(idx_ovl_pi + 53);

    auto ts_y_yzzzzz = pbuffer.data(idx_ovl_pi + 54);

    auto ts_y_zzzzzz = pbuffer.data(idx_ovl_pi + 55);

    auto ts_z_xxxxxx = pbuffer.data(idx_ovl_pi + 56);

    auto ts_z_xxxxxy = pbuffer.data(idx_ovl_pi + 57);

    auto ts_z_xxxxxz = pbuffer.data(idx_ovl_pi + 58);

    auto ts_z_xxxxyy = pbuffer.data(idx_ovl_pi + 59);

    auto ts_z_xxxxyz = pbuffer.data(idx_ovl_pi + 60);

    auto ts_z_xxxxzz = pbuffer.data(idx_ovl_pi + 61);

    auto ts_z_xxxyyy = pbuffer.data(idx_ovl_pi + 62);

    auto ts_z_xxxyyz = pbuffer.data(idx_ovl_pi + 63);

    auto ts_z_xxxyzz = pbuffer.data(idx_ovl_pi + 64);

    auto ts_z_xxxzzz = pbuffer.data(idx_ovl_pi + 65);

    auto ts_z_xxyyyy = pbuffer.data(idx_ovl_pi + 66);

    auto ts_z_xxyyyz = pbuffer.data(idx_ovl_pi + 67);

    auto ts_z_xxyyzz = pbuffer.data(idx_ovl_pi + 68);

    auto ts_z_xxyzzz = pbuffer.data(idx_ovl_pi + 69);

    auto ts_z_xxzzzz = pbuffer.data(idx_ovl_pi + 70);

    auto ts_z_xyyyyy = pbuffer.data(idx_ovl_pi + 71);

    auto ts_z_xyyyyz = pbuffer.data(idx_ovl_pi + 72);

    auto ts_z_xyyyzz = pbuffer.data(idx_ovl_pi + 73);

    auto ts_z_xyyzzz = pbuffer.data(idx_ovl_pi + 74);

    auto ts_z_xyzzzz = pbuffer.data(idx_ovl_pi + 75);

    auto ts_z_xzzzzz = pbuffer.data(idx_ovl_pi + 76);

    auto ts_z_yyyyyy = pbuffer.data(idx_ovl_pi + 77);

    auto ts_z_yyyyyz = pbuffer.data(idx_ovl_pi + 78);

    auto ts_z_yyyyzz = pbuffer.data(idx_ovl_pi + 79);

    auto ts_z_yyyzzz = pbuffer.data(idx_ovl_pi + 80);

    auto ts_z_yyzzzz = pbuffer.data(idx_ovl_pi + 81);

    auto ts_z_yzzzzz = pbuffer.data(idx_ovl_pi + 82);

    auto ts_z_zzzzzz = pbuffer.data(idx_ovl_pi + 83);

    // Set up components of auxiliary buffer : PI

    auto tr_x_x_xxxxxx = pbuffer.data(idx_dip_pi);

    auto tr_x_x_xxxxxy = pbuffer.data(idx_dip_pi + 1);

    auto tr_x_x_xxxxxz = pbuffer.data(idx_dip_pi + 2);

    auto tr_x_x_xxxxyy = pbuffer.data(idx_dip_pi + 3);

    auto tr_x_x_xxxxyz = pbuffer.data(idx_dip_pi + 4);

    auto tr_x_x_xxxxzz = pbuffer.data(idx_dip_pi + 5);

    auto tr_x_x_xxxyyy = pbuffer.data(idx_dip_pi + 6);

    auto tr_x_x_xxxyyz = pbuffer.data(idx_dip_pi + 7);

    auto tr_x_x_xxxyzz = pbuffer.data(idx_dip_pi + 8);

    auto tr_x_x_xxxzzz = pbuffer.data(idx_dip_pi + 9);

    auto tr_x_x_xxyyyy = pbuffer.data(idx_dip_pi + 10);

    auto tr_x_x_xxyyyz = pbuffer.data(idx_dip_pi + 11);

    auto tr_x_x_xxyyzz = pbuffer.data(idx_dip_pi + 12);

    auto tr_x_x_xxyzzz = pbuffer.data(idx_dip_pi + 13);

    auto tr_x_x_xxzzzz = pbuffer.data(idx_dip_pi + 14);

    auto tr_x_x_xyyyyy = pbuffer.data(idx_dip_pi + 15);

    auto tr_x_x_xyyyyz = pbuffer.data(idx_dip_pi + 16);

    auto tr_x_x_xyyyzz = pbuffer.data(idx_dip_pi + 17);

    auto tr_x_x_xyyzzz = pbuffer.data(idx_dip_pi + 18);

    auto tr_x_x_xyzzzz = pbuffer.data(idx_dip_pi + 19);

    auto tr_x_x_xzzzzz = pbuffer.data(idx_dip_pi + 20);

    auto tr_x_x_yyyyyy = pbuffer.data(idx_dip_pi + 21);

    auto tr_x_x_yyyyyz = pbuffer.data(idx_dip_pi + 22);

    auto tr_x_x_yyyyzz = pbuffer.data(idx_dip_pi + 23);

    auto tr_x_x_yyyzzz = pbuffer.data(idx_dip_pi + 24);

    auto tr_x_x_yyzzzz = pbuffer.data(idx_dip_pi + 25);

    auto tr_x_x_yzzzzz = pbuffer.data(idx_dip_pi + 26);

    auto tr_x_x_zzzzzz = pbuffer.data(idx_dip_pi + 27);

    auto tr_x_y_xxxxxx = pbuffer.data(idx_dip_pi + 28);

    auto tr_x_y_xxxxxy = pbuffer.data(idx_dip_pi + 29);

    auto tr_x_y_xxxxxz = pbuffer.data(idx_dip_pi + 30);

    auto tr_x_y_xxxxyy = pbuffer.data(idx_dip_pi + 31);

    auto tr_x_y_xxxxyz = pbuffer.data(idx_dip_pi + 32);

    auto tr_x_y_xxxxzz = pbuffer.data(idx_dip_pi + 33);

    auto tr_x_y_xxxyyy = pbuffer.data(idx_dip_pi + 34);

    auto tr_x_y_xxxyyz = pbuffer.data(idx_dip_pi + 35);

    auto tr_x_y_xxxyzz = pbuffer.data(idx_dip_pi + 36);

    auto tr_x_y_xxxzzz = pbuffer.data(idx_dip_pi + 37);

    auto tr_x_y_xxyyyy = pbuffer.data(idx_dip_pi + 38);

    auto tr_x_y_xxyyyz = pbuffer.data(idx_dip_pi + 39);

    auto tr_x_y_xxyyzz = pbuffer.data(idx_dip_pi + 40);

    auto tr_x_y_xxyzzz = pbuffer.data(idx_dip_pi + 41);

    auto tr_x_y_xxzzzz = pbuffer.data(idx_dip_pi + 42);

    auto tr_x_y_xyyyyy = pbuffer.data(idx_dip_pi + 43);

    auto tr_x_y_xyyyyz = pbuffer.data(idx_dip_pi + 44);

    auto tr_x_y_xyyyzz = pbuffer.data(idx_dip_pi + 45);

    auto tr_x_y_xyyzzz = pbuffer.data(idx_dip_pi + 46);

    auto tr_x_y_xyzzzz = pbuffer.data(idx_dip_pi + 47);

    auto tr_x_y_xzzzzz = pbuffer.data(idx_dip_pi + 48);

    auto tr_x_y_yyyyyy = pbuffer.data(idx_dip_pi + 49);

    auto tr_x_y_yyyyyz = pbuffer.data(idx_dip_pi + 50);

    auto tr_x_y_yyyyzz = pbuffer.data(idx_dip_pi + 51);

    auto tr_x_y_yyyzzz = pbuffer.data(idx_dip_pi + 52);

    auto tr_x_y_yyzzzz = pbuffer.data(idx_dip_pi + 53);

    auto tr_x_y_yzzzzz = pbuffer.data(idx_dip_pi + 54);

    auto tr_x_y_zzzzzz = pbuffer.data(idx_dip_pi + 55);

    auto tr_x_z_xxxxxx = pbuffer.data(idx_dip_pi + 56);

    auto tr_x_z_xxxxxy = pbuffer.data(idx_dip_pi + 57);

    auto tr_x_z_xxxxxz = pbuffer.data(idx_dip_pi + 58);

    auto tr_x_z_xxxxyy = pbuffer.data(idx_dip_pi + 59);

    auto tr_x_z_xxxxyz = pbuffer.data(idx_dip_pi + 60);

    auto tr_x_z_xxxxzz = pbuffer.data(idx_dip_pi + 61);

    auto tr_x_z_xxxyyy = pbuffer.data(idx_dip_pi + 62);

    auto tr_x_z_xxxyyz = pbuffer.data(idx_dip_pi + 63);

    auto tr_x_z_xxxyzz = pbuffer.data(idx_dip_pi + 64);

    auto tr_x_z_xxxzzz = pbuffer.data(idx_dip_pi + 65);

    auto tr_x_z_xxyyyy = pbuffer.data(idx_dip_pi + 66);

    auto tr_x_z_xxyyyz = pbuffer.data(idx_dip_pi + 67);

    auto tr_x_z_xxyyzz = pbuffer.data(idx_dip_pi + 68);

    auto tr_x_z_xxyzzz = pbuffer.data(idx_dip_pi + 69);

    auto tr_x_z_xxzzzz = pbuffer.data(idx_dip_pi + 70);

    auto tr_x_z_xyyyyy = pbuffer.data(idx_dip_pi + 71);

    auto tr_x_z_xyyyyz = pbuffer.data(idx_dip_pi + 72);

    auto tr_x_z_xyyyzz = pbuffer.data(idx_dip_pi + 73);

    auto tr_x_z_xyyzzz = pbuffer.data(idx_dip_pi + 74);

    auto tr_x_z_xyzzzz = pbuffer.data(idx_dip_pi + 75);

    auto tr_x_z_xzzzzz = pbuffer.data(idx_dip_pi + 76);

    auto tr_x_z_yyyyyy = pbuffer.data(idx_dip_pi + 77);

    auto tr_x_z_yyyyyz = pbuffer.data(idx_dip_pi + 78);

    auto tr_x_z_yyyyzz = pbuffer.data(idx_dip_pi + 79);

    auto tr_x_z_yyyzzz = pbuffer.data(idx_dip_pi + 80);

    auto tr_x_z_yyzzzz = pbuffer.data(idx_dip_pi + 81);

    auto tr_x_z_yzzzzz = pbuffer.data(idx_dip_pi + 82);

    auto tr_x_z_zzzzzz = pbuffer.data(idx_dip_pi + 83);

    auto tr_y_x_xxxxxx = pbuffer.data(idx_dip_pi + 84);

    auto tr_y_x_xxxxxy = pbuffer.data(idx_dip_pi + 85);

    auto tr_y_x_xxxxxz = pbuffer.data(idx_dip_pi + 86);

    auto tr_y_x_xxxxyy = pbuffer.data(idx_dip_pi + 87);

    auto tr_y_x_xxxxyz = pbuffer.data(idx_dip_pi + 88);

    auto tr_y_x_xxxxzz = pbuffer.data(idx_dip_pi + 89);

    auto tr_y_x_xxxyyy = pbuffer.data(idx_dip_pi + 90);

    auto tr_y_x_xxxyyz = pbuffer.data(idx_dip_pi + 91);

    auto tr_y_x_xxxyzz = pbuffer.data(idx_dip_pi + 92);

    auto tr_y_x_xxxzzz = pbuffer.data(idx_dip_pi + 93);

    auto tr_y_x_xxyyyy = pbuffer.data(idx_dip_pi + 94);

    auto tr_y_x_xxyyyz = pbuffer.data(idx_dip_pi + 95);

    auto tr_y_x_xxyyzz = pbuffer.data(idx_dip_pi + 96);

    auto tr_y_x_xxyzzz = pbuffer.data(idx_dip_pi + 97);

    auto tr_y_x_xxzzzz = pbuffer.data(idx_dip_pi + 98);

    auto tr_y_x_xyyyyy = pbuffer.data(idx_dip_pi + 99);

    auto tr_y_x_xyyyyz = pbuffer.data(idx_dip_pi + 100);

    auto tr_y_x_xyyyzz = pbuffer.data(idx_dip_pi + 101);

    auto tr_y_x_xyyzzz = pbuffer.data(idx_dip_pi + 102);

    auto tr_y_x_xyzzzz = pbuffer.data(idx_dip_pi + 103);

    auto tr_y_x_xzzzzz = pbuffer.data(idx_dip_pi + 104);

    auto tr_y_x_yyyyyy = pbuffer.data(idx_dip_pi + 105);

    auto tr_y_x_yyyyyz = pbuffer.data(idx_dip_pi + 106);

    auto tr_y_x_yyyyzz = pbuffer.data(idx_dip_pi + 107);

    auto tr_y_x_yyyzzz = pbuffer.data(idx_dip_pi + 108);

    auto tr_y_x_yyzzzz = pbuffer.data(idx_dip_pi + 109);

    auto tr_y_x_yzzzzz = pbuffer.data(idx_dip_pi + 110);

    auto tr_y_x_zzzzzz = pbuffer.data(idx_dip_pi + 111);

    auto tr_y_y_xxxxxx = pbuffer.data(idx_dip_pi + 112);

    auto tr_y_y_xxxxxy = pbuffer.data(idx_dip_pi + 113);

    auto tr_y_y_xxxxxz = pbuffer.data(idx_dip_pi + 114);

    auto tr_y_y_xxxxyy = pbuffer.data(idx_dip_pi + 115);

    auto tr_y_y_xxxxyz = pbuffer.data(idx_dip_pi + 116);

    auto tr_y_y_xxxxzz = pbuffer.data(idx_dip_pi + 117);

    auto tr_y_y_xxxyyy = pbuffer.data(idx_dip_pi + 118);

    auto tr_y_y_xxxyyz = pbuffer.data(idx_dip_pi + 119);

    auto tr_y_y_xxxyzz = pbuffer.data(idx_dip_pi + 120);

    auto tr_y_y_xxxzzz = pbuffer.data(idx_dip_pi + 121);

    auto tr_y_y_xxyyyy = pbuffer.data(idx_dip_pi + 122);

    auto tr_y_y_xxyyyz = pbuffer.data(idx_dip_pi + 123);

    auto tr_y_y_xxyyzz = pbuffer.data(idx_dip_pi + 124);

    auto tr_y_y_xxyzzz = pbuffer.data(idx_dip_pi + 125);

    auto tr_y_y_xxzzzz = pbuffer.data(idx_dip_pi + 126);

    auto tr_y_y_xyyyyy = pbuffer.data(idx_dip_pi + 127);

    auto tr_y_y_xyyyyz = pbuffer.data(idx_dip_pi + 128);

    auto tr_y_y_xyyyzz = pbuffer.data(idx_dip_pi + 129);

    auto tr_y_y_xyyzzz = pbuffer.data(idx_dip_pi + 130);

    auto tr_y_y_xyzzzz = pbuffer.data(idx_dip_pi + 131);

    auto tr_y_y_xzzzzz = pbuffer.data(idx_dip_pi + 132);

    auto tr_y_y_yyyyyy = pbuffer.data(idx_dip_pi + 133);

    auto tr_y_y_yyyyyz = pbuffer.data(idx_dip_pi + 134);

    auto tr_y_y_yyyyzz = pbuffer.data(idx_dip_pi + 135);

    auto tr_y_y_yyyzzz = pbuffer.data(idx_dip_pi + 136);

    auto tr_y_y_yyzzzz = pbuffer.data(idx_dip_pi + 137);

    auto tr_y_y_yzzzzz = pbuffer.data(idx_dip_pi + 138);

    auto tr_y_y_zzzzzz = pbuffer.data(idx_dip_pi + 139);

    auto tr_y_z_xxxxxx = pbuffer.data(idx_dip_pi + 140);

    auto tr_y_z_xxxxxy = pbuffer.data(idx_dip_pi + 141);

    auto tr_y_z_xxxxxz = pbuffer.data(idx_dip_pi + 142);

    auto tr_y_z_xxxxyy = pbuffer.data(idx_dip_pi + 143);

    auto tr_y_z_xxxxyz = pbuffer.data(idx_dip_pi + 144);

    auto tr_y_z_xxxxzz = pbuffer.data(idx_dip_pi + 145);

    auto tr_y_z_xxxyyy = pbuffer.data(idx_dip_pi + 146);

    auto tr_y_z_xxxyyz = pbuffer.data(idx_dip_pi + 147);

    auto tr_y_z_xxxyzz = pbuffer.data(idx_dip_pi + 148);

    auto tr_y_z_xxxzzz = pbuffer.data(idx_dip_pi + 149);

    auto tr_y_z_xxyyyy = pbuffer.data(idx_dip_pi + 150);

    auto tr_y_z_xxyyyz = pbuffer.data(idx_dip_pi + 151);

    auto tr_y_z_xxyyzz = pbuffer.data(idx_dip_pi + 152);

    auto tr_y_z_xxyzzz = pbuffer.data(idx_dip_pi + 153);

    auto tr_y_z_xxzzzz = pbuffer.data(idx_dip_pi + 154);

    auto tr_y_z_xyyyyy = pbuffer.data(idx_dip_pi + 155);

    auto tr_y_z_xyyyyz = pbuffer.data(idx_dip_pi + 156);

    auto tr_y_z_xyyyzz = pbuffer.data(idx_dip_pi + 157);

    auto tr_y_z_xyyzzz = pbuffer.data(idx_dip_pi + 158);

    auto tr_y_z_xyzzzz = pbuffer.data(idx_dip_pi + 159);

    auto tr_y_z_xzzzzz = pbuffer.data(idx_dip_pi + 160);

    auto tr_y_z_yyyyyy = pbuffer.data(idx_dip_pi + 161);

    auto tr_y_z_yyyyyz = pbuffer.data(idx_dip_pi + 162);

    auto tr_y_z_yyyyzz = pbuffer.data(idx_dip_pi + 163);

    auto tr_y_z_yyyzzz = pbuffer.data(idx_dip_pi + 164);

    auto tr_y_z_yyzzzz = pbuffer.data(idx_dip_pi + 165);

    auto tr_y_z_yzzzzz = pbuffer.data(idx_dip_pi + 166);

    auto tr_y_z_zzzzzz = pbuffer.data(idx_dip_pi + 167);

    auto tr_z_x_xxxxxx = pbuffer.data(idx_dip_pi + 168);

    auto tr_z_x_xxxxxy = pbuffer.data(idx_dip_pi + 169);

    auto tr_z_x_xxxxxz = pbuffer.data(idx_dip_pi + 170);

    auto tr_z_x_xxxxyy = pbuffer.data(idx_dip_pi + 171);

    auto tr_z_x_xxxxyz = pbuffer.data(idx_dip_pi + 172);

    auto tr_z_x_xxxxzz = pbuffer.data(idx_dip_pi + 173);

    auto tr_z_x_xxxyyy = pbuffer.data(idx_dip_pi + 174);

    auto tr_z_x_xxxyyz = pbuffer.data(idx_dip_pi + 175);

    auto tr_z_x_xxxyzz = pbuffer.data(idx_dip_pi + 176);

    auto tr_z_x_xxxzzz = pbuffer.data(idx_dip_pi + 177);

    auto tr_z_x_xxyyyy = pbuffer.data(idx_dip_pi + 178);

    auto tr_z_x_xxyyyz = pbuffer.data(idx_dip_pi + 179);

    auto tr_z_x_xxyyzz = pbuffer.data(idx_dip_pi + 180);

    auto tr_z_x_xxyzzz = pbuffer.data(idx_dip_pi + 181);

    auto tr_z_x_xxzzzz = pbuffer.data(idx_dip_pi + 182);

    auto tr_z_x_xyyyyy = pbuffer.data(idx_dip_pi + 183);

    auto tr_z_x_xyyyyz = pbuffer.data(idx_dip_pi + 184);

    auto tr_z_x_xyyyzz = pbuffer.data(idx_dip_pi + 185);

    auto tr_z_x_xyyzzz = pbuffer.data(idx_dip_pi + 186);

    auto tr_z_x_xyzzzz = pbuffer.data(idx_dip_pi + 187);

    auto tr_z_x_xzzzzz = pbuffer.data(idx_dip_pi + 188);

    auto tr_z_x_yyyyyy = pbuffer.data(idx_dip_pi + 189);

    auto tr_z_x_yyyyyz = pbuffer.data(idx_dip_pi + 190);

    auto tr_z_x_yyyyzz = pbuffer.data(idx_dip_pi + 191);

    auto tr_z_x_yyyzzz = pbuffer.data(idx_dip_pi + 192);

    auto tr_z_x_yyzzzz = pbuffer.data(idx_dip_pi + 193);

    auto tr_z_x_yzzzzz = pbuffer.data(idx_dip_pi + 194);

    auto tr_z_x_zzzzzz = pbuffer.data(idx_dip_pi + 195);

    auto tr_z_y_xxxxxx = pbuffer.data(idx_dip_pi + 196);

    auto tr_z_y_xxxxxy = pbuffer.data(idx_dip_pi + 197);

    auto tr_z_y_xxxxxz = pbuffer.data(idx_dip_pi + 198);

    auto tr_z_y_xxxxyy = pbuffer.data(idx_dip_pi + 199);

    auto tr_z_y_xxxxyz = pbuffer.data(idx_dip_pi + 200);

    auto tr_z_y_xxxxzz = pbuffer.data(idx_dip_pi + 201);

    auto tr_z_y_xxxyyy = pbuffer.data(idx_dip_pi + 202);

    auto tr_z_y_xxxyyz = pbuffer.data(idx_dip_pi + 203);

    auto tr_z_y_xxxyzz = pbuffer.data(idx_dip_pi + 204);

    auto tr_z_y_xxxzzz = pbuffer.data(idx_dip_pi + 205);

    auto tr_z_y_xxyyyy = pbuffer.data(idx_dip_pi + 206);

    auto tr_z_y_xxyyyz = pbuffer.data(idx_dip_pi + 207);

    auto tr_z_y_xxyyzz = pbuffer.data(idx_dip_pi + 208);

    auto tr_z_y_xxyzzz = pbuffer.data(idx_dip_pi + 209);

    auto tr_z_y_xxzzzz = pbuffer.data(idx_dip_pi + 210);

    auto tr_z_y_xyyyyy = pbuffer.data(idx_dip_pi + 211);

    auto tr_z_y_xyyyyz = pbuffer.data(idx_dip_pi + 212);

    auto tr_z_y_xyyyzz = pbuffer.data(idx_dip_pi + 213);

    auto tr_z_y_xyyzzz = pbuffer.data(idx_dip_pi + 214);

    auto tr_z_y_xyzzzz = pbuffer.data(idx_dip_pi + 215);

    auto tr_z_y_xzzzzz = pbuffer.data(idx_dip_pi + 216);

    auto tr_z_y_yyyyyy = pbuffer.data(idx_dip_pi + 217);

    auto tr_z_y_yyyyyz = pbuffer.data(idx_dip_pi + 218);

    auto tr_z_y_yyyyzz = pbuffer.data(idx_dip_pi + 219);

    auto tr_z_y_yyyzzz = pbuffer.data(idx_dip_pi + 220);

    auto tr_z_y_yyzzzz = pbuffer.data(idx_dip_pi + 221);

    auto tr_z_y_yzzzzz = pbuffer.data(idx_dip_pi + 222);

    auto tr_z_y_zzzzzz = pbuffer.data(idx_dip_pi + 223);

    auto tr_z_z_xxxxxx = pbuffer.data(idx_dip_pi + 224);

    auto tr_z_z_xxxxxy = pbuffer.data(idx_dip_pi + 225);

    auto tr_z_z_xxxxxz = pbuffer.data(idx_dip_pi + 226);

    auto tr_z_z_xxxxyy = pbuffer.data(idx_dip_pi + 227);

    auto tr_z_z_xxxxyz = pbuffer.data(idx_dip_pi + 228);

    auto tr_z_z_xxxxzz = pbuffer.data(idx_dip_pi + 229);

    auto tr_z_z_xxxyyy = pbuffer.data(idx_dip_pi + 230);

    auto tr_z_z_xxxyyz = pbuffer.data(idx_dip_pi + 231);

    auto tr_z_z_xxxyzz = pbuffer.data(idx_dip_pi + 232);

    auto tr_z_z_xxxzzz = pbuffer.data(idx_dip_pi + 233);

    auto tr_z_z_xxyyyy = pbuffer.data(idx_dip_pi + 234);

    auto tr_z_z_xxyyyz = pbuffer.data(idx_dip_pi + 235);

    auto tr_z_z_xxyyzz = pbuffer.data(idx_dip_pi + 236);

    auto tr_z_z_xxyzzz = pbuffer.data(idx_dip_pi + 237);

    auto tr_z_z_xxzzzz = pbuffer.data(idx_dip_pi + 238);

    auto tr_z_z_xyyyyy = pbuffer.data(idx_dip_pi + 239);

    auto tr_z_z_xyyyyz = pbuffer.data(idx_dip_pi + 240);

    auto tr_z_z_xyyyzz = pbuffer.data(idx_dip_pi + 241);

    auto tr_z_z_xyyzzz = pbuffer.data(idx_dip_pi + 242);

    auto tr_z_z_xyzzzz = pbuffer.data(idx_dip_pi + 243);

    auto tr_z_z_xzzzzz = pbuffer.data(idx_dip_pi + 244);

    auto tr_z_z_yyyyyy = pbuffer.data(idx_dip_pi + 245);

    auto tr_z_z_yyyyyz = pbuffer.data(idx_dip_pi + 246);

    auto tr_z_z_yyyyzz = pbuffer.data(idx_dip_pi + 247);

    auto tr_z_z_yyyzzz = pbuffer.data(idx_dip_pi + 248);

    auto tr_z_z_yyzzzz = pbuffer.data(idx_dip_pi + 249);

    auto tr_z_z_yzzzzz = pbuffer.data(idx_dip_pi + 250);

    auto tr_z_z_zzzzzz = pbuffer.data(idx_dip_pi + 251);

    // Set up 0-28 components of targeted buffer : DI

    auto tr_x_xx_xxxxxx = pbuffer.data(idx_dip_di);

    auto tr_x_xx_xxxxxy = pbuffer.data(idx_dip_di + 1);

    auto tr_x_xx_xxxxxz = pbuffer.data(idx_dip_di + 2);

    auto tr_x_xx_xxxxyy = pbuffer.data(idx_dip_di + 3);

    auto tr_x_xx_xxxxyz = pbuffer.data(idx_dip_di + 4);

    auto tr_x_xx_xxxxzz = pbuffer.data(idx_dip_di + 5);

    auto tr_x_xx_xxxyyy = pbuffer.data(idx_dip_di + 6);

    auto tr_x_xx_xxxyyz = pbuffer.data(idx_dip_di + 7);

    auto tr_x_xx_xxxyzz = pbuffer.data(idx_dip_di + 8);

    auto tr_x_xx_xxxzzz = pbuffer.data(idx_dip_di + 9);

    auto tr_x_xx_xxyyyy = pbuffer.data(idx_dip_di + 10);

    auto tr_x_xx_xxyyyz = pbuffer.data(idx_dip_di + 11);

    auto tr_x_xx_xxyyzz = pbuffer.data(idx_dip_di + 12);

    auto tr_x_xx_xxyzzz = pbuffer.data(idx_dip_di + 13);

    auto tr_x_xx_xxzzzz = pbuffer.data(idx_dip_di + 14);

    auto tr_x_xx_xyyyyy = pbuffer.data(idx_dip_di + 15);

    auto tr_x_xx_xyyyyz = pbuffer.data(idx_dip_di + 16);

    auto tr_x_xx_xyyyzz = pbuffer.data(idx_dip_di + 17);

    auto tr_x_xx_xyyzzz = pbuffer.data(idx_dip_di + 18);

    auto tr_x_xx_xyzzzz = pbuffer.data(idx_dip_di + 19);

    auto tr_x_xx_xzzzzz = pbuffer.data(idx_dip_di + 20);

    auto tr_x_xx_yyyyyy = pbuffer.data(idx_dip_di + 21);

    auto tr_x_xx_yyyyyz = pbuffer.data(idx_dip_di + 22);

    auto tr_x_xx_yyyyzz = pbuffer.data(idx_dip_di + 23);

    auto tr_x_xx_yyyzzz = pbuffer.data(idx_dip_di + 24);

    auto tr_x_xx_yyzzzz = pbuffer.data(idx_dip_di + 25);

    auto tr_x_xx_yzzzzz = pbuffer.data(idx_dip_di + 26);

    auto tr_x_xx_zzzzzz = pbuffer.data(idx_dip_di + 27);

#pragma omp simd aligned(pa_x,               \
                             tr_x_0_xxxxxx,  \
                             tr_x_0_xxxxxy,  \
                             tr_x_0_xxxxxz,  \
                             tr_x_0_xxxxyy,  \
                             tr_x_0_xxxxyz,  \
                             tr_x_0_xxxxzz,  \
                             tr_x_0_xxxyyy,  \
                             tr_x_0_xxxyyz,  \
                             tr_x_0_xxxyzz,  \
                             tr_x_0_xxxzzz,  \
                             tr_x_0_xxyyyy,  \
                             tr_x_0_xxyyyz,  \
                             tr_x_0_xxyyzz,  \
                             tr_x_0_xxyzzz,  \
                             tr_x_0_xxzzzz,  \
                             tr_x_0_xyyyyy,  \
                             tr_x_0_xyyyyz,  \
                             tr_x_0_xyyyzz,  \
                             tr_x_0_xyyzzz,  \
                             tr_x_0_xyzzzz,  \
                             tr_x_0_xzzzzz,  \
                             tr_x_0_yyyyyy,  \
                             tr_x_0_yyyyyz,  \
                             tr_x_0_yyyyzz,  \
                             tr_x_0_yyyzzz,  \
                             tr_x_0_yyzzzz,  \
                             tr_x_0_yzzzzz,  \
                             tr_x_0_zzzzzz,  \
                             tr_x_x_xxxxx,   \
                             tr_x_x_xxxxxx,  \
                             tr_x_x_xxxxxy,  \
                             tr_x_x_xxxxxz,  \
                             tr_x_x_xxxxy,   \
                             tr_x_x_xxxxyy,  \
                             tr_x_x_xxxxyz,  \
                             tr_x_x_xxxxz,   \
                             tr_x_x_xxxxzz,  \
                             tr_x_x_xxxyy,   \
                             tr_x_x_xxxyyy,  \
                             tr_x_x_xxxyyz,  \
                             tr_x_x_xxxyz,   \
                             tr_x_x_xxxyzz,  \
                             tr_x_x_xxxzz,   \
                             tr_x_x_xxxzzz,  \
                             tr_x_x_xxyyy,   \
                             tr_x_x_xxyyyy,  \
                             tr_x_x_xxyyyz,  \
                             tr_x_x_xxyyz,   \
                             tr_x_x_xxyyzz,  \
                             tr_x_x_xxyzz,   \
                             tr_x_x_xxyzzz,  \
                             tr_x_x_xxzzz,   \
                             tr_x_x_xxzzzz,  \
                             tr_x_x_xyyyy,   \
                             tr_x_x_xyyyyy,  \
                             tr_x_x_xyyyyz,  \
                             tr_x_x_xyyyz,   \
                             tr_x_x_xyyyzz,  \
                             tr_x_x_xyyzz,   \
                             tr_x_x_xyyzzz,  \
                             tr_x_x_xyzzz,   \
                             tr_x_x_xyzzzz,  \
                             tr_x_x_xzzzz,   \
                             tr_x_x_xzzzzz,  \
                             tr_x_x_yyyyy,   \
                             tr_x_x_yyyyyy,  \
                             tr_x_x_yyyyyz,  \
                             tr_x_x_yyyyz,   \
                             tr_x_x_yyyyzz,  \
                             tr_x_x_yyyzz,   \
                             tr_x_x_yyyzzz,  \
                             tr_x_x_yyzzz,   \
                             tr_x_x_yyzzzz,  \
                             tr_x_x_yzzzz,   \
                             tr_x_x_yzzzzz,  \
                             tr_x_x_zzzzz,   \
                             tr_x_x_zzzzzz,  \
                             tr_x_xx_xxxxxx, \
                             tr_x_xx_xxxxxy, \
                             tr_x_xx_xxxxxz, \
                             tr_x_xx_xxxxyy, \
                             tr_x_xx_xxxxyz, \
                             tr_x_xx_xxxxzz, \
                             tr_x_xx_xxxyyy, \
                             tr_x_xx_xxxyyz, \
                             tr_x_xx_xxxyzz, \
                             tr_x_xx_xxxzzz, \
                             tr_x_xx_xxyyyy, \
                             tr_x_xx_xxyyyz, \
                             tr_x_xx_xxyyzz, \
                             tr_x_xx_xxyzzz, \
                             tr_x_xx_xxzzzz, \
                             tr_x_xx_xyyyyy, \
                             tr_x_xx_xyyyyz, \
                             tr_x_xx_xyyyzz, \
                             tr_x_xx_xyyzzz, \
                             tr_x_xx_xyzzzz, \
                             tr_x_xx_xzzzzz, \
                             tr_x_xx_yyyyyy, \
                             tr_x_xx_yyyyyz, \
                             tr_x_xx_yyyyzz, \
                             tr_x_xx_yyyzzz, \
                             tr_x_xx_yyzzzz, \
                             tr_x_xx_yzzzzz, \
                             tr_x_xx_zzzzzz, \
                             ts_x_xxxxxx,    \
                             ts_x_xxxxxy,    \
                             ts_x_xxxxxz,    \
                             ts_x_xxxxyy,    \
                             ts_x_xxxxyz,    \
                             ts_x_xxxxzz,    \
                             ts_x_xxxyyy,    \
                             ts_x_xxxyyz,    \
                             ts_x_xxxyzz,    \
                             ts_x_xxxzzz,    \
                             ts_x_xxyyyy,    \
                             ts_x_xxyyyz,    \
                             ts_x_xxyyzz,    \
                             ts_x_xxyzzz,    \
                             ts_x_xxzzzz,    \
                             ts_x_xyyyyy,    \
                             ts_x_xyyyyz,    \
                             ts_x_xyyyzz,    \
                             ts_x_xyyzzz,    \
                             ts_x_xyzzzz,    \
                             ts_x_xzzzzz,    \
                             ts_x_yyyyyy,    \
                             ts_x_yyyyyz,    \
                             ts_x_yyyyzz,    \
                             ts_x_yyyzzz,    \
                             ts_x_yyzzzz,    \
                             ts_x_yzzzzz,    \
                             ts_x_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xx_xxxxxx[i] = tr_x_0_xxxxxx[i] * fe_0 + 6.0 * tr_x_x_xxxxx[i] * fe_0 + ts_x_xxxxxx[i] * fe_0 + tr_x_x_xxxxxx[i] * pa_x[i];

        tr_x_xx_xxxxxy[i] = tr_x_0_xxxxxy[i] * fe_0 + 5.0 * tr_x_x_xxxxy[i] * fe_0 + ts_x_xxxxxy[i] * fe_0 + tr_x_x_xxxxxy[i] * pa_x[i];

        tr_x_xx_xxxxxz[i] = tr_x_0_xxxxxz[i] * fe_0 + 5.0 * tr_x_x_xxxxz[i] * fe_0 + ts_x_xxxxxz[i] * fe_0 + tr_x_x_xxxxxz[i] * pa_x[i];

        tr_x_xx_xxxxyy[i] = tr_x_0_xxxxyy[i] * fe_0 + 4.0 * tr_x_x_xxxyy[i] * fe_0 + ts_x_xxxxyy[i] * fe_0 + tr_x_x_xxxxyy[i] * pa_x[i];

        tr_x_xx_xxxxyz[i] = tr_x_0_xxxxyz[i] * fe_0 + 4.0 * tr_x_x_xxxyz[i] * fe_0 + ts_x_xxxxyz[i] * fe_0 + tr_x_x_xxxxyz[i] * pa_x[i];

        tr_x_xx_xxxxzz[i] = tr_x_0_xxxxzz[i] * fe_0 + 4.0 * tr_x_x_xxxzz[i] * fe_0 + ts_x_xxxxzz[i] * fe_0 + tr_x_x_xxxxzz[i] * pa_x[i];

        tr_x_xx_xxxyyy[i] = tr_x_0_xxxyyy[i] * fe_0 + 3.0 * tr_x_x_xxyyy[i] * fe_0 + ts_x_xxxyyy[i] * fe_0 + tr_x_x_xxxyyy[i] * pa_x[i];

        tr_x_xx_xxxyyz[i] = tr_x_0_xxxyyz[i] * fe_0 + 3.0 * tr_x_x_xxyyz[i] * fe_0 + ts_x_xxxyyz[i] * fe_0 + tr_x_x_xxxyyz[i] * pa_x[i];

        tr_x_xx_xxxyzz[i] = tr_x_0_xxxyzz[i] * fe_0 + 3.0 * tr_x_x_xxyzz[i] * fe_0 + ts_x_xxxyzz[i] * fe_0 + tr_x_x_xxxyzz[i] * pa_x[i];

        tr_x_xx_xxxzzz[i] = tr_x_0_xxxzzz[i] * fe_0 + 3.0 * tr_x_x_xxzzz[i] * fe_0 + ts_x_xxxzzz[i] * fe_0 + tr_x_x_xxxzzz[i] * pa_x[i];

        tr_x_xx_xxyyyy[i] = tr_x_0_xxyyyy[i] * fe_0 + 2.0 * tr_x_x_xyyyy[i] * fe_0 + ts_x_xxyyyy[i] * fe_0 + tr_x_x_xxyyyy[i] * pa_x[i];

        tr_x_xx_xxyyyz[i] = tr_x_0_xxyyyz[i] * fe_0 + 2.0 * tr_x_x_xyyyz[i] * fe_0 + ts_x_xxyyyz[i] * fe_0 + tr_x_x_xxyyyz[i] * pa_x[i];

        tr_x_xx_xxyyzz[i] = tr_x_0_xxyyzz[i] * fe_0 + 2.0 * tr_x_x_xyyzz[i] * fe_0 + ts_x_xxyyzz[i] * fe_0 + tr_x_x_xxyyzz[i] * pa_x[i];

        tr_x_xx_xxyzzz[i] = tr_x_0_xxyzzz[i] * fe_0 + 2.0 * tr_x_x_xyzzz[i] * fe_0 + ts_x_xxyzzz[i] * fe_0 + tr_x_x_xxyzzz[i] * pa_x[i];

        tr_x_xx_xxzzzz[i] = tr_x_0_xxzzzz[i] * fe_0 + 2.0 * tr_x_x_xzzzz[i] * fe_0 + ts_x_xxzzzz[i] * fe_0 + tr_x_x_xxzzzz[i] * pa_x[i];

        tr_x_xx_xyyyyy[i] = tr_x_0_xyyyyy[i] * fe_0 + tr_x_x_yyyyy[i] * fe_0 + ts_x_xyyyyy[i] * fe_0 + tr_x_x_xyyyyy[i] * pa_x[i];

        tr_x_xx_xyyyyz[i] = tr_x_0_xyyyyz[i] * fe_0 + tr_x_x_yyyyz[i] * fe_0 + ts_x_xyyyyz[i] * fe_0 + tr_x_x_xyyyyz[i] * pa_x[i];

        tr_x_xx_xyyyzz[i] = tr_x_0_xyyyzz[i] * fe_0 + tr_x_x_yyyzz[i] * fe_0 + ts_x_xyyyzz[i] * fe_0 + tr_x_x_xyyyzz[i] * pa_x[i];

        tr_x_xx_xyyzzz[i] = tr_x_0_xyyzzz[i] * fe_0 + tr_x_x_yyzzz[i] * fe_0 + ts_x_xyyzzz[i] * fe_0 + tr_x_x_xyyzzz[i] * pa_x[i];

        tr_x_xx_xyzzzz[i] = tr_x_0_xyzzzz[i] * fe_0 + tr_x_x_yzzzz[i] * fe_0 + ts_x_xyzzzz[i] * fe_0 + tr_x_x_xyzzzz[i] * pa_x[i];

        tr_x_xx_xzzzzz[i] = tr_x_0_xzzzzz[i] * fe_0 + tr_x_x_zzzzz[i] * fe_0 + ts_x_xzzzzz[i] * fe_0 + tr_x_x_xzzzzz[i] * pa_x[i];

        tr_x_xx_yyyyyy[i] = tr_x_0_yyyyyy[i] * fe_0 + ts_x_yyyyyy[i] * fe_0 + tr_x_x_yyyyyy[i] * pa_x[i];

        tr_x_xx_yyyyyz[i] = tr_x_0_yyyyyz[i] * fe_0 + ts_x_yyyyyz[i] * fe_0 + tr_x_x_yyyyyz[i] * pa_x[i];

        tr_x_xx_yyyyzz[i] = tr_x_0_yyyyzz[i] * fe_0 + ts_x_yyyyzz[i] * fe_0 + tr_x_x_yyyyzz[i] * pa_x[i];

        tr_x_xx_yyyzzz[i] = tr_x_0_yyyzzz[i] * fe_0 + ts_x_yyyzzz[i] * fe_0 + tr_x_x_yyyzzz[i] * pa_x[i];

        tr_x_xx_yyzzzz[i] = tr_x_0_yyzzzz[i] * fe_0 + ts_x_yyzzzz[i] * fe_0 + tr_x_x_yyzzzz[i] * pa_x[i];

        tr_x_xx_yzzzzz[i] = tr_x_0_yzzzzz[i] * fe_0 + ts_x_yzzzzz[i] * fe_0 + tr_x_x_yzzzzz[i] * pa_x[i];

        tr_x_xx_zzzzzz[i] = tr_x_0_zzzzzz[i] * fe_0 + ts_x_zzzzzz[i] * fe_0 + tr_x_x_zzzzzz[i] * pa_x[i];
    }

    // Set up 28-56 components of targeted buffer : DI

    auto tr_x_xy_xxxxxx = pbuffer.data(idx_dip_di + 28);

    auto tr_x_xy_xxxxxy = pbuffer.data(idx_dip_di + 29);

    auto tr_x_xy_xxxxxz = pbuffer.data(idx_dip_di + 30);

    auto tr_x_xy_xxxxyy = pbuffer.data(idx_dip_di + 31);

    auto tr_x_xy_xxxxyz = pbuffer.data(idx_dip_di + 32);

    auto tr_x_xy_xxxxzz = pbuffer.data(idx_dip_di + 33);

    auto tr_x_xy_xxxyyy = pbuffer.data(idx_dip_di + 34);

    auto tr_x_xy_xxxyyz = pbuffer.data(idx_dip_di + 35);

    auto tr_x_xy_xxxyzz = pbuffer.data(idx_dip_di + 36);

    auto tr_x_xy_xxxzzz = pbuffer.data(idx_dip_di + 37);

    auto tr_x_xy_xxyyyy = pbuffer.data(idx_dip_di + 38);

    auto tr_x_xy_xxyyyz = pbuffer.data(idx_dip_di + 39);

    auto tr_x_xy_xxyyzz = pbuffer.data(idx_dip_di + 40);

    auto tr_x_xy_xxyzzz = pbuffer.data(idx_dip_di + 41);

    auto tr_x_xy_xxzzzz = pbuffer.data(idx_dip_di + 42);

    auto tr_x_xy_xyyyyy = pbuffer.data(idx_dip_di + 43);

    auto tr_x_xy_xyyyyz = pbuffer.data(idx_dip_di + 44);

    auto tr_x_xy_xyyyzz = pbuffer.data(idx_dip_di + 45);

    auto tr_x_xy_xyyzzz = pbuffer.data(idx_dip_di + 46);

    auto tr_x_xy_xyzzzz = pbuffer.data(idx_dip_di + 47);

    auto tr_x_xy_xzzzzz = pbuffer.data(idx_dip_di + 48);

    auto tr_x_xy_yyyyyy = pbuffer.data(idx_dip_di + 49);

    auto tr_x_xy_yyyyyz = pbuffer.data(idx_dip_di + 50);

    auto tr_x_xy_yyyyzz = pbuffer.data(idx_dip_di + 51);

    auto tr_x_xy_yyyzzz = pbuffer.data(idx_dip_di + 52);

    auto tr_x_xy_yyzzzz = pbuffer.data(idx_dip_di + 53);

    auto tr_x_xy_yzzzzz = pbuffer.data(idx_dip_di + 54);

    auto tr_x_xy_zzzzzz = pbuffer.data(idx_dip_di + 55);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_x_x_xxxxx,   \
                             tr_x_x_xxxxxx,  \
                             tr_x_x_xxxxxy,  \
                             tr_x_x_xxxxxz,  \
                             tr_x_x_xxxxy,   \
                             tr_x_x_xxxxyy,  \
                             tr_x_x_xxxxyz,  \
                             tr_x_x_xxxxz,   \
                             tr_x_x_xxxxzz,  \
                             tr_x_x_xxxyy,   \
                             tr_x_x_xxxyyy,  \
                             tr_x_x_xxxyyz,  \
                             tr_x_x_xxxyz,   \
                             tr_x_x_xxxyzz,  \
                             tr_x_x_xxxzz,   \
                             tr_x_x_xxxzzz,  \
                             tr_x_x_xxyyy,   \
                             tr_x_x_xxyyyy,  \
                             tr_x_x_xxyyyz,  \
                             tr_x_x_xxyyz,   \
                             tr_x_x_xxyyzz,  \
                             tr_x_x_xxyzz,   \
                             tr_x_x_xxyzzz,  \
                             tr_x_x_xxzzz,   \
                             tr_x_x_xxzzzz,  \
                             tr_x_x_xyyyy,   \
                             tr_x_x_xyyyyy,  \
                             tr_x_x_xyyyyz,  \
                             tr_x_x_xyyyz,   \
                             tr_x_x_xyyyzz,  \
                             tr_x_x_xyyzz,   \
                             tr_x_x_xyyzzz,  \
                             tr_x_x_xyzzz,   \
                             tr_x_x_xyzzzz,  \
                             tr_x_x_xzzzz,   \
                             tr_x_x_xzzzzz,  \
                             tr_x_x_zzzzzz,  \
                             tr_x_xy_xxxxxx, \
                             tr_x_xy_xxxxxy, \
                             tr_x_xy_xxxxxz, \
                             tr_x_xy_xxxxyy, \
                             tr_x_xy_xxxxyz, \
                             tr_x_xy_xxxxzz, \
                             tr_x_xy_xxxyyy, \
                             tr_x_xy_xxxyyz, \
                             tr_x_xy_xxxyzz, \
                             tr_x_xy_xxxzzz, \
                             tr_x_xy_xxyyyy, \
                             tr_x_xy_xxyyyz, \
                             tr_x_xy_xxyyzz, \
                             tr_x_xy_xxyzzz, \
                             tr_x_xy_xxzzzz, \
                             tr_x_xy_xyyyyy, \
                             tr_x_xy_xyyyyz, \
                             tr_x_xy_xyyyzz, \
                             tr_x_xy_xyyzzz, \
                             tr_x_xy_xyzzzz, \
                             tr_x_xy_xzzzzz, \
                             tr_x_xy_yyyyyy, \
                             tr_x_xy_yyyyyz, \
                             tr_x_xy_yyyyzz, \
                             tr_x_xy_yyyzzz, \
                             tr_x_xy_yyzzzz, \
                             tr_x_xy_yzzzzz, \
                             tr_x_xy_zzzzzz, \
                             tr_x_y_yyyyyy,  \
                             tr_x_y_yyyyyz,  \
                             tr_x_y_yyyyzz,  \
                             tr_x_y_yyyzzz,  \
                             tr_x_y_yyzzzz,  \
                             tr_x_y_yzzzzz,  \
                             ts_y_yyyyyy,    \
                             ts_y_yyyyyz,    \
                             ts_y_yyyyzz,    \
                             ts_y_yyyzzz,    \
                             ts_y_yyzzzz,    \
                             ts_y_yzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xy_xxxxxx[i] = tr_x_x_xxxxxx[i] * pa_y[i];

        tr_x_xy_xxxxxy[i] = tr_x_x_xxxxx[i] * fe_0 + tr_x_x_xxxxxy[i] * pa_y[i];

        tr_x_xy_xxxxxz[i] = tr_x_x_xxxxxz[i] * pa_y[i];

        tr_x_xy_xxxxyy[i] = 2.0 * tr_x_x_xxxxy[i] * fe_0 + tr_x_x_xxxxyy[i] * pa_y[i];

        tr_x_xy_xxxxyz[i] = tr_x_x_xxxxz[i] * fe_0 + tr_x_x_xxxxyz[i] * pa_y[i];

        tr_x_xy_xxxxzz[i] = tr_x_x_xxxxzz[i] * pa_y[i];

        tr_x_xy_xxxyyy[i] = 3.0 * tr_x_x_xxxyy[i] * fe_0 + tr_x_x_xxxyyy[i] * pa_y[i];

        tr_x_xy_xxxyyz[i] = 2.0 * tr_x_x_xxxyz[i] * fe_0 + tr_x_x_xxxyyz[i] * pa_y[i];

        tr_x_xy_xxxyzz[i] = tr_x_x_xxxzz[i] * fe_0 + tr_x_x_xxxyzz[i] * pa_y[i];

        tr_x_xy_xxxzzz[i] = tr_x_x_xxxzzz[i] * pa_y[i];

        tr_x_xy_xxyyyy[i] = 4.0 * tr_x_x_xxyyy[i] * fe_0 + tr_x_x_xxyyyy[i] * pa_y[i];

        tr_x_xy_xxyyyz[i] = 3.0 * tr_x_x_xxyyz[i] * fe_0 + tr_x_x_xxyyyz[i] * pa_y[i];

        tr_x_xy_xxyyzz[i] = 2.0 * tr_x_x_xxyzz[i] * fe_0 + tr_x_x_xxyyzz[i] * pa_y[i];

        tr_x_xy_xxyzzz[i] = tr_x_x_xxzzz[i] * fe_0 + tr_x_x_xxyzzz[i] * pa_y[i];

        tr_x_xy_xxzzzz[i] = tr_x_x_xxzzzz[i] * pa_y[i];

        tr_x_xy_xyyyyy[i] = 5.0 * tr_x_x_xyyyy[i] * fe_0 + tr_x_x_xyyyyy[i] * pa_y[i];

        tr_x_xy_xyyyyz[i] = 4.0 * tr_x_x_xyyyz[i] * fe_0 + tr_x_x_xyyyyz[i] * pa_y[i];

        tr_x_xy_xyyyzz[i] = 3.0 * tr_x_x_xyyzz[i] * fe_0 + tr_x_x_xyyyzz[i] * pa_y[i];

        tr_x_xy_xyyzzz[i] = 2.0 * tr_x_x_xyzzz[i] * fe_0 + tr_x_x_xyyzzz[i] * pa_y[i];

        tr_x_xy_xyzzzz[i] = tr_x_x_xzzzz[i] * fe_0 + tr_x_x_xyzzzz[i] * pa_y[i];

        tr_x_xy_xzzzzz[i] = tr_x_x_xzzzzz[i] * pa_y[i];

        tr_x_xy_yyyyyy[i] = ts_y_yyyyyy[i] * fe_0 + tr_x_y_yyyyyy[i] * pa_x[i];

        tr_x_xy_yyyyyz[i] = ts_y_yyyyyz[i] * fe_0 + tr_x_y_yyyyyz[i] * pa_x[i];

        tr_x_xy_yyyyzz[i] = ts_y_yyyyzz[i] * fe_0 + tr_x_y_yyyyzz[i] * pa_x[i];

        tr_x_xy_yyyzzz[i] = ts_y_yyyzzz[i] * fe_0 + tr_x_y_yyyzzz[i] * pa_x[i];

        tr_x_xy_yyzzzz[i] = ts_y_yyzzzz[i] * fe_0 + tr_x_y_yyzzzz[i] * pa_x[i];

        tr_x_xy_yzzzzz[i] = ts_y_yzzzzz[i] * fe_0 + tr_x_y_yzzzzz[i] * pa_x[i];

        tr_x_xy_zzzzzz[i] = tr_x_x_zzzzzz[i] * pa_y[i];
    }

    // Set up 56-84 components of targeted buffer : DI

    auto tr_x_xz_xxxxxx = pbuffer.data(idx_dip_di + 56);

    auto tr_x_xz_xxxxxy = pbuffer.data(idx_dip_di + 57);

    auto tr_x_xz_xxxxxz = pbuffer.data(idx_dip_di + 58);

    auto tr_x_xz_xxxxyy = pbuffer.data(idx_dip_di + 59);

    auto tr_x_xz_xxxxyz = pbuffer.data(idx_dip_di + 60);

    auto tr_x_xz_xxxxzz = pbuffer.data(idx_dip_di + 61);

    auto tr_x_xz_xxxyyy = pbuffer.data(idx_dip_di + 62);

    auto tr_x_xz_xxxyyz = pbuffer.data(idx_dip_di + 63);

    auto tr_x_xz_xxxyzz = pbuffer.data(idx_dip_di + 64);

    auto tr_x_xz_xxxzzz = pbuffer.data(idx_dip_di + 65);

    auto tr_x_xz_xxyyyy = pbuffer.data(idx_dip_di + 66);

    auto tr_x_xz_xxyyyz = pbuffer.data(idx_dip_di + 67);

    auto tr_x_xz_xxyyzz = pbuffer.data(idx_dip_di + 68);

    auto tr_x_xz_xxyzzz = pbuffer.data(idx_dip_di + 69);

    auto tr_x_xz_xxzzzz = pbuffer.data(idx_dip_di + 70);

    auto tr_x_xz_xyyyyy = pbuffer.data(idx_dip_di + 71);

    auto tr_x_xz_xyyyyz = pbuffer.data(idx_dip_di + 72);

    auto tr_x_xz_xyyyzz = pbuffer.data(idx_dip_di + 73);

    auto tr_x_xz_xyyzzz = pbuffer.data(idx_dip_di + 74);

    auto tr_x_xz_xyzzzz = pbuffer.data(idx_dip_di + 75);

    auto tr_x_xz_xzzzzz = pbuffer.data(idx_dip_di + 76);

    auto tr_x_xz_yyyyyy = pbuffer.data(idx_dip_di + 77);

    auto tr_x_xz_yyyyyz = pbuffer.data(idx_dip_di + 78);

    auto tr_x_xz_yyyyzz = pbuffer.data(idx_dip_di + 79);

    auto tr_x_xz_yyyzzz = pbuffer.data(idx_dip_di + 80);

    auto tr_x_xz_yyzzzz = pbuffer.data(idx_dip_di + 81);

    auto tr_x_xz_yzzzzz = pbuffer.data(idx_dip_di + 82);

    auto tr_x_xz_zzzzzz = pbuffer.data(idx_dip_di + 83);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_x_x_xxxxx,   \
                             tr_x_x_xxxxxx,  \
                             tr_x_x_xxxxxy,  \
                             tr_x_x_xxxxxz,  \
                             tr_x_x_xxxxy,   \
                             tr_x_x_xxxxyy,  \
                             tr_x_x_xxxxyz,  \
                             tr_x_x_xxxxz,   \
                             tr_x_x_xxxxzz,  \
                             tr_x_x_xxxyy,   \
                             tr_x_x_xxxyyy,  \
                             tr_x_x_xxxyyz,  \
                             tr_x_x_xxxyz,   \
                             tr_x_x_xxxyzz,  \
                             tr_x_x_xxxzz,   \
                             tr_x_x_xxxzzz,  \
                             tr_x_x_xxyyy,   \
                             tr_x_x_xxyyyy,  \
                             tr_x_x_xxyyyz,  \
                             tr_x_x_xxyyz,   \
                             tr_x_x_xxyyzz,  \
                             tr_x_x_xxyzz,   \
                             tr_x_x_xxyzzz,  \
                             tr_x_x_xxzzz,   \
                             tr_x_x_xxzzzz,  \
                             tr_x_x_xyyyy,   \
                             tr_x_x_xyyyyy,  \
                             tr_x_x_xyyyyz,  \
                             tr_x_x_xyyyz,   \
                             tr_x_x_xyyyzz,  \
                             tr_x_x_xyyzz,   \
                             tr_x_x_xyyzzz,  \
                             tr_x_x_xyzzz,   \
                             tr_x_x_xyzzzz,  \
                             tr_x_x_xzzzz,   \
                             tr_x_x_xzzzzz,  \
                             tr_x_x_yyyyyy,  \
                             tr_x_xz_xxxxxx, \
                             tr_x_xz_xxxxxy, \
                             tr_x_xz_xxxxxz, \
                             tr_x_xz_xxxxyy, \
                             tr_x_xz_xxxxyz, \
                             tr_x_xz_xxxxzz, \
                             tr_x_xz_xxxyyy, \
                             tr_x_xz_xxxyyz, \
                             tr_x_xz_xxxyzz, \
                             tr_x_xz_xxxzzz, \
                             tr_x_xz_xxyyyy, \
                             tr_x_xz_xxyyyz, \
                             tr_x_xz_xxyyzz, \
                             tr_x_xz_xxyzzz, \
                             tr_x_xz_xxzzzz, \
                             tr_x_xz_xyyyyy, \
                             tr_x_xz_xyyyyz, \
                             tr_x_xz_xyyyzz, \
                             tr_x_xz_xyyzzz, \
                             tr_x_xz_xyzzzz, \
                             tr_x_xz_xzzzzz, \
                             tr_x_xz_yyyyyy, \
                             tr_x_xz_yyyyyz, \
                             tr_x_xz_yyyyzz, \
                             tr_x_xz_yyyzzz, \
                             tr_x_xz_yyzzzz, \
                             tr_x_xz_yzzzzz, \
                             tr_x_xz_zzzzzz, \
                             tr_x_z_yyyyyz,  \
                             tr_x_z_yyyyzz,  \
                             tr_x_z_yyyzzz,  \
                             tr_x_z_yyzzzz,  \
                             tr_x_z_yzzzzz,  \
                             tr_x_z_zzzzzz,  \
                             ts_z_yyyyyz,    \
                             ts_z_yyyyzz,    \
                             ts_z_yyyzzz,    \
                             ts_z_yyzzzz,    \
                             ts_z_yzzzzz,    \
                             ts_z_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_xz_xxxxxx[i] = tr_x_x_xxxxxx[i] * pa_z[i];

        tr_x_xz_xxxxxy[i] = tr_x_x_xxxxxy[i] * pa_z[i];

        tr_x_xz_xxxxxz[i] = tr_x_x_xxxxx[i] * fe_0 + tr_x_x_xxxxxz[i] * pa_z[i];

        tr_x_xz_xxxxyy[i] = tr_x_x_xxxxyy[i] * pa_z[i];

        tr_x_xz_xxxxyz[i] = tr_x_x_xxxxy[i] * fe_0 + tr_x_x_xxxxyz[i] * pa_z[i];

        tr_x_xz_xxxxzz[i] = 2.0 * tr_x_x_xxxxz[i] * fe_0 + tr_x_x_xxxxzz[i] * pa_z[i];

        tr_x_xz_xxxyyy[i] = tr_x_x_xxxyyy[i] * pa_z[i];

        tr_x_xz_xxxyyz[i] = tr_x_x_xxxyy[i] * fe_0 + tr_x_x_xxxyyz[i] * pa_z[i];

        tr_x_xz_xxxyzz[i] = 2.0 * tr_x_x_xxxyz[i] * fe_0 + tr_x_x_xxxyzz[i] * pa_z[i];

        tr_x_xz_xxxzzz[i] = 3.0 * tr_x_x_xxxzz[i] * fe_0 + tr_x_x_xxxzzz[i] * pa_z[i];

        tr_x_xz_xxyyyy[i] = tr_x_x_xxyyyy[i] * pa_z[i];

        tr_x_xz_xxyyyz[i] = tr_x_x_xxyyy[i] * fe_0 + tr_x_x_xxyyyz[i] * pa_z[i];

        tr_x_xz_xxyyzz[i] = 2.0 * tr_x_x_xxyyz[i] * fe_0 + tr_x_x_xxyyzz[i] * pa_z[i];

        tr_x_xz_xxyzzz[i] = 3.0 * tr_x_x_xxyzz[i] * fe_0 + tr_x_x_xxyzzz[i] * pa_z[i];

        tr_x_xz_xxzzzz[i] = 4.0 * tr_x_x_xxzzz[i] * fe_0 + tr_x_x_xxzzzz[i] * pa_z[i];

        tr_x_xz_xyyyyy[i] = tr_x_x_xyyyyy[i] * pa_z[i];

        tr_x_xz_xyyyyz[i] = tr_x_x_xyyyy[i] * fe_0 + tr_x_x_xyyyyz[i] * pa_z[i];

        tr_x_xz_xyyyzz[i] = 2.0 * tr_x_x_xyyyz[i] * fe_0 + tr_x_x_xyyyzz[i] * pa_z[i];

        tr_x_xz_xyyzzz[i] = 3.0 * tr_x_x_xyyzz[i] * fe_0 + tr_x_x_xyyzzz[i] * pa_z[i];

        tr_x_xz_xyzzzz[i] = 4.0 * tr_x_x_xyzzz[i] * fe_0 + tr_x_x_xyzzzz[i] * pa_z[i];

        tr_x_xz_xzzzzz[i] = 5.0 * tr_x_x_xzzzz[i] * fe_0 + tr_x_x_xzzzzz[i] * pa_z[i];

        tr_x_xz_yyyyyy[i] = tr_x_x_yyyyyy[i] * pa_z[i];

        tr_x_xz_yyyyyz[i] = ts_z_yyyyyz[i] * fe_0 + tr_x_z_yyyyyz[i] * pa_x[i];

        tr_x_xz_yyyyzz[i] = ts_z_yyyyzz[i] * fe_0 + tr_x_z_yyyyzz[i] * pa_x[i];

        tr_x_xz_yyyzzz[i] = ts_z_yyyzzz[i] * fe_0 + tr_x_z_yyyzzz[i] * pa_x[i];

        tr_x_xz_yyzzzz[i] = ts_z_yyzzzz[i] * fe_0 + tr_x_z_yyzzzz[i] * pa_x[i];

        tr_x_xz_yzzzzz[i] = ts_z_yzzzzz[i] * fe_0 + tr_x_z_yzzzzz[i] * pa_x[i];

        tr_x_xz_zzzzzz[i] = ts_z_zzzzzz[i] * fe_0 + tr_x_z_zzzzzz[i] * pa_x[i];
    }

    // Set up 84-112 components of targeted buffer : DI

    auto tr_x_yy_xxxxxx = pbuffer.data(idx_dip_di + 84);

    auto tr_x_yy_xxxxxy = pbuffer.data(idx_dip_di + 85);

    auto tr_x_yy_xxxxxz = pbuffer.data(idx_dip_di + 86);

    auto tr_x_yy_xxxxyy = pbuffer.data(idx_dip_di + 87);

    auto tr_x_yy_xxxxyz = pbuffer.data(idx_dip_di + 88);

    auto tr_x_yy_xxxxzz = pbuffer.data(idx_dip_di + 89);

    auto tr_x_yy_xxxyyy = pbuffer.data(idx_dip_di + 90);

    auto tr_x_yy_xxxyyz = pbuffer.data(idx_dip_di + 91);

    auto tr_x_yy_xxxyzz = pbuffer.data(idx_dip_di + 92);

    auto tr_x_yy_xxxzzz = pbuffer.data(idx_dip_di + 93);

    auto tr_x_yy_xxyyyy = pbuffer.data(idx_dip_di + 94);

    auto tr_x_yy_xxyyyz = pbuffer.data(idx_dip_di + 95);

    auto tr_x_yy_xxyyzz = pbuffer.data(idx_dip_di + 96);

    auto tr_x_yy_xxyzzz = pbuffer.data(idx_dip_di + 97);

    auto tr_x_yy_xxzzzz = pbuffer.data(idx_dip_di + 98);

    auto tr_x_yy_xyyyyy = pbuffer.data(idx_dip_di + 99);

    auto tr_x_yy_xyyyyz = pbuffer.data(idx_dip_di + 100);

    auto tr_x_yy_xyyyzz = pbuffer.data(idx_dip_di + 101);

    auto tr_x_yy_xyyzzz = pbuffer.data(idx_dip_di + 102);

    auto tr_x_yy_xyzzzz = pbuffer.data(idx_dip_di + 103);

    auto tr_x_yy_xzzzzz = pbuffer.data(idx_dip_di + 104);

    auto tr_x_yy_yyyyyy = pbuffer.data(idx_dip_di + 105);

    auto tr_x_yy_yyyyyz = pbuffer.data(idx_dip_di + 106);

    auto tr_x_yy_yyyyzz = pbuffer.data(idx_dip_di + 107);

    auto tr_x_yy_yyyzzz = pbuffer.data(idx_dip_di + 108);

    auto tr_x_yy_yyzzzz = pbuffer.data(idx_dip_di + 109);

    auto tr_x_yy_yzzzzz = pbuffer.data(idx_dip_di + 110);

    auto tr_x_yy_zzzzzz = pbuffer.data(idx_dip_di + 111);

#pragma omp simd aligned(pa_y,               \
                             tr_x_0_xxxxxx,  \
                             tr_x_0_xxxxxy,  \
                             tr_x_0_xxxxxz,  \
                             tr_x_0_xxxxyy,  \
                             tr_x_0_xxxxyz,  \
                             tr_x_0_xxxxzz,  \
                             tr_x_0_xxxyyy,  \
                             tr_x_0_xxxyyz,  \
                             tr_x_0_xxxyzz,  \
                             tr_x_0_xxxzzz,  \
                             tr_x_0_xxyyyy,  \
                             tr_x_0_xxyyyz,  \
                             tr_x_0_xxyyzz,  \
                             tr_x_0_xxyzzz,  \
                             tr_x_0_xxzzzz,  \
                             tr_x_0_xyyyyy,  \
                             tr_x_0_xyyyyz,  \
                             tr_x_0_xyyyzz,  \
                             tr_x_0_xyyzzz,  \
                             tr_x_0_xyzzzz,  \
                             tr_x_0_xzzzzz,  \
                             tr_x_0_yyyyyy,  \
                             tr_x_0_yyyyyz,  \
                             tr_x_0_yyyyzz,  \
                             tr_x_0_yyyzzz,  \
                             tr_x_0_yyzzzz,  \
                             tr_x_0_yzzzzz,  \
                             tr_x_0_zzzzzz,  \
                             tr_x_y_xxxxx,   \
                             tr_x_y_xxxxxx,  \
                             tr_x_y_xxxxxy,  \
                             tr_x_y_xxxxxz,  \
                             tr_x_y_xxxxy,   \
                             tr_x_y_xxxxyy,  \
                             tr_x_y_xxxxyz,  \
                             tr_x_y_xxxxz,   \
                             tr_x_y_xxxxzz,  \
                             tr_x_y_xxxyy,   \
                             tr_x_y_xxxyyy,  \
                             tr_x_y_xxxyyz,  \
                             tr_x_y_xxxyz,   \
                             tr_x_y_xxxyzz,  \
                             tr_x_y_xxxzz,   \
                             tr_x_y_xxxzzz,  \
                             tr_x_y_xxyyy,   \
                             tr_x_y_xxyyyy,  \
                             tr_x_y_xxyyyz,  \
                             tr_x_y_xxyyz,   \
                             tr_x_y_xxyyzz,  \
                             tr_x_y_xxyzz,   \
                             tr_x_y_xxyzzz,  \
                             tr_x_y_xxzzz,   \
                             tr_x_y_xxzzzz,  \
                             tr_x_y_xyyyy,   \
                             tr_x_y_xyyyyy,  \
                             tr_x_y_xyyyyz,  \
                             tr_x_y_xyyyz,   \
                             tr_x_y_xyyyzz,  \
                             tr_x_y_xyyzz,   \
                             tr_x_y_xyyzzz,  \
                             tr_x_y_xyzzz,   \
                             tr_x_y_xyzzzz,  \
                             tr_x_y_xzzzz,   \
                             tr_x_y_xzzzzz,  \
                             tr_x_y_yyyyy,   \
                             tr_x_y_yyyyyy,  \
                             tr_x_y_yyyyyz,  \
                             tr_x_y_yyyyz,   \
                             tr_x_y_yyyyzz,  \
                             tr_x_y_yyyzz,   \
                             tr_x_y_yyyzzz,  \
                             tr_x_y_yyzzz,   \
                             tr_x_y_yyzzzz,  \
                             tr_x_y_yzzzz,   \
                             tr_x_y_yzzzzz,  \
                             tr_x_y_zzzzz,   \
                             tr_x_y_zzzzzz,  \
                             tr_x_yy_xxxxxx, \
                             tr_x_yy_xxxxxy, \
                             tr_x_yy_xxxxxz, \
                             tr_x_yy_xxxxyy, \
                             tr_x_yy_xxxxyz, \
                             tr_x_yy_xxxxzz, \
                             tr_x_yy_xxxyyy, \
                             tr_x_yy_xxxyyz, \
                             tr_x_yy_xxxyzz, \
                             tr_x_yy_xxxzzz, \
                             tr_x_yy_xxyyyy, \
                             tr_x_yy_xxyyyz, \
                             tr_x_yy_xxyyzz, \
                             tr_x_yy_xxyzzz, \
                             tr_x_yy_xxzzzz, \
                             tr_x_yy_xyyyyy, \
                             tr_x_yy_xyyyyz, \
                             tr_x_yy_xyyyzz, \
                             tr_x_yy_xyyzzz, \
                             tr_x_yy_xyzzzz, \
                             tr_x_yy_xzzzzz, \
                             tr_x_yy_yyyyyy, \
                             tr_x_yy_yyyyyz, \
                             tr_x_yy_yyyyzz, \
                             tr_x_yy_yyyzzz, \
                             tr_x_yy_yyzzzz, \
                             tr_x_yy_yzzzzz, \
                             tr_x_yy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yy_xxxxxx[i] = tr_x_0_xxxxxx[i] * fe_0 + tr_x_y_xxxxxx[i] * pa_y[i];

        tr_x_yy_xxxxxy[i] = tr_x_0_xxxxxy[i] * fe_0 + tr_x_y_xxxxx[i] * fe_0 + tr_x_y_xxxxxy[i] * pa_y[i];

        tr_x_yy_xxxxxz[i] = tr_x_0_xxxxxz[i] * fe_0 + tr_x_y_xxxxxz[i] * pa_y[i];

        tr_x_yy_xxxxyy[i] = tr_x_0_xxxxyy[i] * fe_0 + 2.0 * tr_x_y_xxxxy[i] * fe_0 + tr_x_y_xxxxyy[i] * pa_y[i];

        tr_x_yy_xxxxyz[i] = tr_x_0_xxxxyz[i] * fe_0 + tr_x_y_xxxxz[i] * fe_0 + tr_x_y_xxxxyz[i] * pa_y[i];

        tr_x_yy_xxxxzz[i] = tr_x_0_xxxxzz[i] * fe_0 + tr_x_y_xxxxzz[i] * pa_y[i];

        tr_x_yy_xxxyyy[i] = tr_x_0_xxxyyy[i] * fe_0 + 3.0 * tr_x_y_xxxyy[i] * fe_0 + tr_x_y_xxxyyy[i] * pa_y[i];

        tr_x_yy_xxxyyz[i] = tr_x_0_xxxyyz[i] * fe_0 + 2.0 * tr_x_y_xxxyz[i] * fe_0 + tr_x_y_xxxyyz[i] * pa_y[i];

        tr_x_yy_xxxyzz[i] = tr_x_0_xxxyzz[i] * fe_0 + tr_x_y_xxxzz[i] * fe_0 + tr_x_y_xxxyzz[i] * pa_y[i];

        tr_x_yy_xxxzzz[i] = tr_x_0_xxxzzz[i] * fe_0 + tr_x_y_xxxzzz[i] * pa_y[i];

        tr_x_yy_xxyyyy[i] = tr_x_0_xxyyyy[i] * fe_0 + 4.0 * tr_x_y_xxyyy[i] * fe_0 + tr_x_y_xxyyyy[i] * pa_y[i];

        tr_x_yy_xxyyyz[i] = tr_x_0_xxyyyz[i] * fe_0 + 3.0 * tr_x_y_xxyyz[i] * fe_0 + tr_x_y_xxyyyz[i] * pa_y[i];

        tr_x_yy_xxyyzz[i] = tr_x_0_xxyyzz[i] * fe_0 + 2.0 * tr_x_y_xxyzz[i] * fe_0 + tr_x_y_xxyyzz[i] * pa_y[i];

        tr_x_yy_xxyzzz[i] = tr_x_0_xxyzzz[i] * fe_0 + tr_x_y_xxzzz[i] * fe_0 + tr_x_y_xxyzzz[i] * pa_y[i];

        tr_x_yy_xxzzzz[i] = tr_x_0_xxzzzz[i] * fe_0 + tr_x_y_xxzzzz[i] * pa_y[i];

        tr_x_yy_xyyyyy[i] = tr_x_0_xyyyyy[i] * fe_0 + 5.0 * tr_x_y_xyyyy[i] * fe_0 + tr_x_y_xyyyyy[i] * pa_y[i];

        tr_x_yy_xyyyyz[i] = tr_x_0_xyyyyz[i] * fe_0 + 4.0 * tr_x_y_xyyyz[i] * fe_0 + tr_x_y_xyyyyz[i] * pa_y[i];

        tr_x_yy_xyyyzz[i] = tr_x_0_xyyyzz[i] * fe_0 + 3.0 * tr_x_y_xyyzz[i] * fe_0 + tr_x_y_xyyyzz[i] * pa_y[i];

        tr_x_yy_xyyzzz[i] = tr_x_0_xyyzzz[i] * fe_0 + 2.0 * tr_x_y_xyzzz[i] * fe_0 + tr_x_y_xyyzzz[i] * pa_y[i];

        tr_x_yy_xyzzzz[i] = tr_x_0_xyzzzz[i] * fe_0 + tr_x_y_xzzzz[i] * fe_0 + tr_x_y_xyzzzz[i] * pa_y[i];

        tr_x_yy_xzzzzz[i] = tr_x_0_xzzzzz[i] * fe_0 + tr_x_y_xzzzzz[i] * pa_y[i];

        tr_x_yy_yyyyyy[i] = tr_x_0_yyyyyy[i] * fe_0 + 6.0 * tr_x_y_yyyyy[i] * fe_0 + tr_x_y_yyyyyy[i] * pa_y[i];

        tr_x_yy_yyyyyz[i] = tr_x_0_yyyyyz[i] * fe_0 + 5.0 * tr_x_y_yyyyz[i] * fe_0 + tr_x_y_yyyyyz[i] * pa_y[i];

        tr_x_yy_yyyyzz[i] = tr_x_0_yyyyzz[i] * fe_0 + 4.0 * tr_x_y_yyyzz[i] * fe_0 + tr_x_y_yyyyzz[i] * pa_y[i];

        tr_x_yy_yyyzzz[i] = tr_x_0_yyyzzz[i] * fe_0 + 3.0 * tr_x_y_yyzzz[i] * fe_0 + tr_x_y_yyyzzz[i] * pa_y[i];

        tr_x_yy_yyzzzz[i] = tr_x_0_yyzzzz[i] * fe_0 + 2.0 * tr_x_y_yzzzz[i] * fe_0 + tr_x_y_yyzzzz[i] * pa_y[i];

        tr_x_yy_yzzzzz[i] = tr_x_0_yzzzzz[i] * fe_0 + tr_x_y_zzzzz[i] * fe_0 + tr_x_y_yzzzzz[i] * pa_y[i];

        tr_x_yy_zzzzzz[i] = tr_x_0_zzzzzz[i] * fe_0 + tr_x_y_zzzzzz[i] * pa_y[i];
    }

    // Set up 112-140 components of targeted buffer : DI

    auto tr_x_yz_xxxxxx = pbuffer.data(idx_dip_di + 112);

    auto tr_x_yz_xxxxxy = pbuffer.data(idx_dip_di + 113);

    auto tr_x_yz_xxxxxz = pbuffer.data(idx_dip_di + 114);

    auto tr_x_yz_xxxxyy = pbuffer.data(idx_dip_di + 115);

    auto tr_x_yz_xxxxyz = pbuffer.data(idx_dip_di + 116);

    auto tr_x_yz_xxxxzz = pbuffer.data(idx_dip_di + 117);

    auto tr_x_yz_xxxyyy = pbuffer.data(idx_dip_di + 118);

    auto tr_x_yz_xxxyyz = pbuffer.data(idx_dip_di + 119);

    auto tr_x_yz_xxxyzz = pbuffer.data(idx_dip_di + 120);

    auto tr_x_yz_xxxzzz = pbuffer.data(idx_dip_di + 121);

    auto tr_x_yz_xxyyyy = pbuffer.data(idx_dip_di + 122);

    auto tr_x_yz_xxyyyz = pbuffer.data(idx_dip_di + 123);

    auto tr_x_yz_xxyyzz = pbuffer.data(idx_dip_di + 124);

    auto tr_x_yz_xxyzzz = pbuffer.data(idx_dip_di + 125);

    auto tr_x_yz_xxzzzz = pbuffer.data(idx_dip_di + 126);

    auto tr_x_yz_xyyyyy = pbuffer.data(idx_dip_di + 127);

    auto tr_x_yz_xyyyyz = pbuffer.data(idx_dip_di + 128);

    auto tr_x_yz_xyyyzz = pbuffer.data(idx_dip_di + 129);

    auto tr_x_yz_xyyzzz = pbuffer.data(idx_dip_di + 130);

    auto tr_x_yz_xyzzzz = pbuffer.data(idx_dip_di + 131);

    auto tr_x_yz_xzzzzz = pbuffer.data(idx_dip_di + 132);

    auto tr_x_yz_yyyyyy = pbuffer.data(idx_dip_di + 133);

    auto tr_x_yz_yyyyyz = pbuffer.data(idx_dip_di + 134);

    auto tr_x_yz_yyyyzz = pbuffer.data(idx_dip_di + 135);

    auto tr_x_yz_yyyzzz = pbuffer.data(idx_dip_di + 136);

    auto tr_x_yz_yyzzzz = pbuffer.data(idx_dip_di + 137);

    auto tr_x_yz_yzzzzz = pbuffer.data(idx_dip_di + 138);

    auto tr_x_yz_zzzzzz = pbuffer.data(idx_dip_di + 139);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_x_y_xxxxxy,  \
                             tr_x_y_xxxxyy,  \
                             tr_x_y_xxxyyy,  \
                             tr_x_y_xxyyyy,  \
                             tr_x_y_xyyyyy,  \
                             tr_x_y_yyyyyy,  \
                             tr_x_yz_xxxxxx, \
                             tr_x_yz_xxxxxy, \
                             tr_x_yz_xxxxxz, \
                             tr_x_yz_xxxxyy, \
                             tr_x_yz_xxxxyz, \
                             tr_x_yz_xxxxzz, \
                             tr_x_yz_xxxyyy, \
                             tr_x_yz_xxxyyz, \
                             tr_x_yz_xxxyzz, \
                             tr_x_yz_xxxzzz, \
                             tr_x_yz_xxyyyy, \
                             tr_x_yz_xxyyyz, \
                             tr_x_yz_xxyyzz, \
                             tr_x_yz_xxyzzz, \
                             tr_x_yz_xxzzzz, \
                             tr_x_yz_xyyyyy, \
                             tr_x_yz_xyyyyz, \
                             tr_x_yz_xyyyzz, \
                             tr_x_yz_xyyzzz, \
                             tr_x_yz_xyzzzz, \
                             tr_x_yz_xzzzzz, \
                             tr_x_yz_yyyyyy, \
                             tr_x_yz_yyyyyz, \
                             tr_x_yz_yyyyzz, \
                             tr_x_yz_yyyzzz, \
                             tr_x_yz_yyzzzz, \
                             tr_x_yz_yzzzzz, \
                             tr_x_yz_zzzzzz, \
                             tr_x_z_xxxxxx,  \
                             tr_x_z_xxxxxz,  \
                             tr_x_z_xxxxyz,  \
                             tr_x_z_xxxxz,   \
                             tr_x_z_xxxxzz,  \
                             tr_x_z_xxxyyz,  \
                             tr_x_z_xxxyz,   \
                             tr_x_z_xxxyzz,  \
                             tr_x_z_xxxzz,   \
                             tr_x_z_xxxzzz,  \
                             tr_x_z_xxyyyz,  \
                             tr_x_z_xxyyz,   \
                             tr_x_z_xxyyzz,  \
                             tr_x_z_xxyzz,   \
                             tr_x_z_xxyzzz,  \
                             tr_x_z_xxzzz,   \
                             tr_x_z_xxzzzz,  \
                             tr_x_z_xyyyyz,  \
                             tr_x_z_xyyyz,   \
                             tr_x_z_xyyyzz,  \
                             tr_x_z_xyyzz,   \
                             tr_x_z_xyyzzz,  \
                             tr_x_z_xyzzz,   \
                             tr_x_z_xyzzzz,  \
                             tr_x_z_xzzzz,   \
                             tr_x_z_xzzzzz,  \
                             tr_x_z_yyyyyz,  \
                             tr_x_z_yyyyz,   \
                             tr_x_z_yyyyzz,  \
                             tr_x_z_yyyzz,   \
                             tr_x_z_yyyzzz,  \
                             tr_x_z_yyzzz,   \
                             tr_x_z_yyzzzz,  \
                             tr_x_z_yzzzz,   \
                             tr_x_z_yzzzzz,  \
                             tr_x_z_zzzzz,   \
                             tr_x_z_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_yz_xxxxxx[i] = tr_x_z_xxxxxx[i] * pa_y[i];

        tr_x_yz_xxxxxy[i] = tr_x_y_xxxxxy[i] * pa_z[i];

        tr_x_yz_xxxxxz[i] = tr_x_z_xxxxxz[i] * pa_y[i];

        tr_x_yz_xxxxyy[i] = tr_x_y_xxxxyy[i] * pa_z[i];

        tr_x_yz_xxxxyz[i] = tr_x_z_xxxxz[i] * fe_0 + tr_x_z_xxxxyz[i] * pa_y[i];

        tr_x_yz_xxxxzz[i] = tr_x_z_xxxxzz[i] * pa_y[i];

        tr_x_yz_xxxyyy[i] = tr_x_y_xxxyyy[i] * pa_z[i];

        tr_x_yz_xxxyyz[i] = 2.0 * tr_x_z_xxxyz[i] * fe_0 + tr_x_z_xxxyyz[i] * pa_y[i];

        tr_x_yz_xxxyzz[i] = tr_x_z_xxxzz[i] * fe_0 + tr_x_z_xxxyzz[i] * pa_y[i];

        tr_x_yz_xxxzzz[i] = tr_x_z_xxxzzz[i] * pa_y[i];

        tr_x_yz_xxyyyy[i] = tr_x_y_xxyyyy[i] * pa_z[i];

        tr_x_yz_xxyyyz[i] = 3.0 * tr_x_z_xxyyz[i] * fe_0 + tr_x_z_xxyyyz[i] * pa_y[i];

        tr_x_yz_xxyyzz[i] = 2.0 * tr_x_z_xxyzz[i] * fe_0 + tr_x_z_xxyyzz[i] * pa_y[i];

        tr_x_yz_xxyzzz[i] = tr_x_z_xxzzz[i] * fe_0 + tr_x_z_xxyzzz[i] * pa_y[i];

        tr_x_yz_xxzzzz[i] = tr_x_z_xxzzzz[i] * pa_y[i];

        tr_x_yz_xyyyyy[i] = tr_x_y_xyyyyy[i] * pa_z[i];

        tr_x_yz_xyyyyz[i] = 4.0 * tr_x_z_xyyyz[i] * fe_0 + tr_x_z_xyyyyz[i] * pa_y[i];

        tr_x_yz_xyyyzz[i] = 3.0 * tr_x_z_xyyzz[i] * fe_0 + tr_x_z_xyyyzz[i] * pa_y[i];

        tr_x_yz_xyyzzz[i] = 2.0 * tr_x_z_xyzzz[i] * fe_0 + tr_x_z_xyyzzz[i] * pa_y[i];

        tr_x_yz_xyzzzz[i] = tr_x_z_xzzzz[i] * fe_0 + tr_x_z_xyzzzz[i] * pa_y[i];

        tr_x_yz_xzzzzz[i] = tr_x_z_xzzzzz[i] * pa_y[i];

        tr_x_yz_yyyyyy[i] = tr_x_y_yyyyyy[i] * pa_z[i];

        tr_x_yz_yyyyyz[i] = 5.0 * tr_x_z_yyyyz[i] * fe_0 + tr_x_z_yyyyyz[i] * pa_y[i];

        tr_x_yz_yyyyzz[i] = 4.0 * tr_x_z_yyyzz[i] * fe_0 + tr_x_z_yyyyzz[i] * pa_y[i];

        tr_x_yz_yyyzzz[i] = 3.0 * tr_x_z_yyzzz[i] * fe_0 + tr_x_z_yyyzzz[i] * pa_y[i];

        tr_x_yz_yyzzzz[i] = 2.0 * tr_x_z_yzzzz[i] * fe_0 + tr_x_z_yyzzzz[i] * pa_y[i];

        tr_x_yz_yzzzzz[i] = tr_x_z_zzzzz[i] * fe_0 + tr_x_z_yzzzzz[i] * pa_y[i];

        tr_x_yz_zzzzzz[i] = tr_x_z_zzzzzz[i] * pa_y[i];
    }

    // Set up 140-168 components of targeted buffer : DI

    auto tr_x_zz_xxxxxx = pbuffer.data(idx_dip_di + 140);

    auto tr_x_zz_xxxxxy = pbuffer.data(idx_dip_di + 141);

    auto tr_x_zz_xxxxxz = pbuffer.data(idx_dip_di + 142);

    auto tr_x_zz_xxxxyy = pbuffer.data(idx_dip_di + 143);

    auto tr_x_zz_xxxxyz = pbuffer.data(idx_dip_di + 144);

    auto tr_x_zz_xxxxzz = pbuffer.data(idx_dip_di + 145);

    auto tr_x_zz_xxxyyy = pbuffer.data(idx_dip_di + 146);

    auto tr_x_zz_xxxyyz = pbuffer.data(idx_dip_di + 147);

    auto tr_x_zz_xxxyzz = pbuffer.data(idx_dip_di + 148);

    auto tr_x_zz_xxxzzz = pbuffer.data(idx_dip_di + 149);

    auto tr_x_zz_xxyyyy = pbuffer.data(idx_dip_di + 150);

    auto tr_x_zz_xxyyyz = pbuffer.data(idx_dip_di + 151);

    auto tr_x_zz_xxyyzz = pbuffer.data(idx_dip_di + 152);

    auto tr_x_zz_xxyzzz = pbuffer.data(idx_dip_di + 153);

    auto tr_x_zz_xxzzzz = pbuffer.data(idx_dip_di + 154);

    auto tr_x_zz_xyyyyy = pbuffer.data(idx_dip_di + 155);

    auto tr_x_zz_xyyyyz = pbuffer.data(idx_dip_di + 156);

    auto tr_x_zz_xyyyzz = pbuffer.data(idx_dip_di + 157);

    auto tr_x_zz_xyyzzz = pbuffer.data(idx_dip_di + 158);

    auto tr_x_zz_xyzzzz = pbuffer.data(idx_dip_di + 159);

    auto tr_x_zz_xzzzzz = pbuffer.data(idx_dip_di + 160);

    auto tr_x_zz_yyyyyy = pbuffer.data(idx_dip_di + 161);

    auto tr_x_zz_yyyyyz = pbuffer.data(idx_dip_di + 162);

    auto tr_x_zz_yyyyzz = pbuffer.data(idx_dip_di + 163);

    auto tr_x_zz_yyyzzz = pbuffer.data(idx_dip_di + 164);

    auto tr_x_zz_yyzzzz = pbuffer.data(idx_dip_di + 165);

    auto tr_x_zz_yzzzzz = pbuffer.data(idx_dip_di + 166);

    auto tr_x_zz_zzzzzz = pbuffer.data(idx_dip_di + 167);

#pragma omp simd aligned(pa_z,               \
                             tr_x_0_xxxxxx,  \
                             tr_x_0_xxxxxy,  \
                             tr_x_0_xxxxxz,  \
                             tr_x_0_xxxxyy,  \
                             tr_x_0_xxxxyz,  \
                             tr_x_0_xxxxzz,  \
                             tr_x_0_xxxyyy,  \
                             tr_x_0_xxxyyz,  \
                             tr_x_0_xxxyzz,  \
                             tr_x_0_xxxzzz,  \
                             tr_x_0_xxyyyy,  \
                             tr_x_0_xxyyyz,  \
                             tr_x_0_xxyyzz,  \
                             tr_x_0_xxyzzz,  \
                             tr_x_0_xxzzzz,  \
                             tr_x_0_xyyyyy,  \
                             tr_x_0_xyyyyz,  \
                             tr_x_0_xyyyzz,  \
                             tr_x_0_xyyzzz,  \
                             tr_x_0_xyzzzz,  \
                             tr_x_0_xzzzzz,  \
                             tr_x_0_yyyyyy,  \
                             tr_x_0_yyyyyz,  \
                             tr_x_0_yyyyzz,  \
                             tr_x_0_yyyzzz,  \
                             tr_x_0_yyzzzz,  \
                             tr_x_0_yzzzzz,  \
                             tr_x_0_zzzzzz,  \
                             tr_x_z_xxxxx,   \
                             tr_x_z_xxxxxx,  \
                             tr_x_z_xxxxxy,  \
                             tr_x_z_xxxxxz,  \
                             tr_x_z_xxxxy,   \
                             tr_x_z_xxxxyy,  \
                             tr_x_z_xxxxyz,  \
                             tr_x_z_xxxxz,   \
                             tr_x_z_xxxxzz,  \
                             tr_x_z_xxxyy,   \
                             tr_x_z_xxxyyy,  \
                             tr_x_z_xxxyyz,  \
                             tr_x_z_xxxyz,   \
                             tr_x_z_xxxyzz,  \
                             tr_x_z_xxxzz,   \
                             tr_x_z_xxxzzz,  \
                             tr_x_z_xxyyy,   \
                             tr_x_z_xxyyyy,  \
                             tr_x_z_xxyyyz,  \
                             tr_x_z_xxyyz,   \
                             tr_x_z_xxyyzz,  \
                             tr_x_z_xxyzz,   \
                             tr_x_z_xxyzzz,  \
                             tr_x_z_xxzzz,   \
                             tr_x_z_xxzzzz,  \
                             tr_x_z_xyyyy,   \
                             tr_x_z_xyyyyy,  \
                             tr_x_z_xyyyyz,  \
                             tr_x_z_xyyyz,   \
                             tr_x_z_xyyyzz,  \
                             tr_x_z_xyyzz,   \
                             tr_x_z_xyyzzz,  \
                             tr_x_z_xyzzz,   \
                             tr_x_z_xyzzzz,  \
                             tr_x_z_xzzzz,   \
                             tr_x_z_xzzzzz,  \
                             tr_x_z_yyyyy,   \
                             tr_x_z_yyyyyy,  \
                             tr_x_z_yyyyyz,  \
                             tr_x_z_yyyyz,   \
                             tr_x_z_yyyyzz,  \
                             tr_x_z_yyyzz,   \
                             tr_x_z_yyyzzz,  \
                             tr_x_z_yyzzz,   \
                             tr_x_z_yyzzzz,  \
                             tr_x_z_yzzzz,   \
                             tr_x_z_yzzzzz,  \
                             tr_x_z_zzzzz,   \
                             tr_x_z_zzzzzz,  \
                             tr_x_zz_xxxxxx, \
                             tr_x_zz_xxxxxy, \
                             tr_x_zz_xxxxxz, \
                             tr_x_zz_xxxxyy, \
                             tr_x_zz_xxxxyz, \
                             tr_x_zz_xxxxzz, \
                             tr_x_zz_xxxyyy, \
                             tr_x_zz_xxxyyz, \
                             tr_x_zz_xxxyzz, \
                             tr_x_zz_xxxzzz, \
                             tr_x_zz_xxyyyy, \
                             tr_x_zz_xxyyyz, \
                             tr_x_zz_xxyyzz, \
                             tr_x_zz_xxyzzz, \
                             tr_x_zz_xxzzzz, \
                             tr_x_zz_xyyyyy, \
                             tr_x_zz_xyyyyz, \
                             tr_x_zz_xyyyzz, \
                             tr_x_zz_xyyzzz, \
                             tr_x_zz_xyzzzz, \
                             tr_x_zz_xzzzzz, \
                             tr_x_zz_yyyyyy, \
                             tr_x_zz_yyyyyz, \
                             tr_x_zz_yyyyzz, \
                             tr_x_zz_yyyzzz, \
                             tr_x_zz_yyzzzz, \
                             tr_x_zz_yzzzzz, \
                             tr_x_zz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_x_zz_xxxxxx[i] = tr_x_0_xxxxxx[i] * fe_0 + tr_x_z_xxxxxx[i] * pa_z[i];

        tr_x_zz_xxxxxy[i] = tr_x_0_xxxxxy[i] * fe_0 + tr_x_z_xxxxxy[i] * pa_z[i];

        tr_x_zz_xxxxxz[i] = tr_x_0_xxxxxz[i] * fe_0 + tr_x_z_xxxxx[i] * fe_0 + tr_x_z_xxxxxz[i] * pa_z[i];

        tr_x_zz_xxxxyy[i] = tr_x_0_xxxxyy[i] * fe_0 + tr_x_z_xxxxyy[i] * pa_z[i];

        tr_x_zz_xxxxyz[i] = tr_x_0_xxxxyz[i] * fe_0 + tr_x_z_xxxxy[i] * fe_0 + tr_x_z_xxxxyz[i] * pa_z[i];

        tr_x_zz_xxxxzz[i] = tr_x_0_xxxxzz[i] * fe_0 + 2.0 * tr_x_z_xxxxz[i] * fe_0 + tr_x_z_xxxxzz[i] * pa_z[i];

        tr_x_zz_xxxyyy[i] = tr_x_0_xxxyyy[i] * fe_0 + tr_x_z_xxxyyy[i] * pa_z[i];

        tr_x_zz_xxxyyz[i] = tr_x_0_xxxyyz[i] * fe_0 + tr_x_z_xxxyy[i] * fe_0 + tr_x_z_xxxyyz[i] * pa_z[i];

        tr_x_zz_xxxyzz[i] = tr_x_0_xxxyzz[i] * fe_0 + 2.0 * tr_x_z_xxxyz[i] * fe_0 + tr_x_z_xxxyzz[i] * pa_z[i];

        tr_x_zz_xxxzzz[i] = tr_x_0_xxxzzz[i] * fe_0 + 3.0 * tr_x_z_xxxzz[i] * fe_0 + tr_x_z_xxxzzz[i] * pa_z[i];

        tr_x_zz_xxyyyy[i] = tr_x_0_xxyyyy[i] * fe_0 + tr_x_z_xxyyyy[i] * pa_z[i];

        tr_x_zz_xxyyyz[i] = tr_x_0_xxyyyz[i] * fe_0 + tr_x_z_xxyyy[i] * fe_0 + tr_x_z_xxyyyz[i] * pa_z[i];

        tr_x_zz_xxyyzz[i] = tr_x_0_xxyyzz[i] * fe_0 + 2.0 * tr_x_z_xxyyz[i] * fe_0 + tr_x_z_xxyyzz[i] * pa_z[i];

        tr_x_zz_xxyzzz[i] = tr_x_0_xxyzzz[i] * fe_0 + 3.0 * tr_x_z_xxyzz[i] * fe_0 + tr_x_z_xxyzzz[i] * pa_z[i];

        tr_x_zz_xxzzzz[i] = tr_x_0_xxzzzz[i] * fe_0 + 4.0 * tr_x_z_xxzzz[i] * fe_0 + tr_x_z_xxzzzz[i] * pa_z[i];

        tr_x_zz_xyyyyy[i] = tr_x_0_xyyyyy[i] * fe_0 + tr_x_z_xyyyyy[i] * pa_z[i];

        tr_x_zz_xyyyyz[i] = tr_x_0_xyyyyz[i] * fe_0 + tr_x_z_xyyyy[i] * fe_0 + tr_x_z_xyyyyz[i] * pa_z[i];

        tr_x_zz_xyyyzz[i] = tr_x_0_xyyyzz[i] * fe_0 + 2.0 * tr_x_z_xyyyz[i] * fe_0 + tr_x_z_xyyyzz[i] * pa_z[i];

        tr_x_zz_xyyzzz[i] = tr_x_0_xyyzzz[i] * fe_0 + 3.0 * tr_x_z_xyyzz[i] * fe_0 + tr_x_z_xyyzzz[i] * pa_z[i];

        tr_x_zz_xyzzzz[i] = tr_x_0_xyzzzz[i] * fe_0 + 4.0 * tr_x_z_xyzzz[i] * fe_0 + tr_x_z_xyzzzz[i] * pa_z[i];

        tr_x_zz_xzzzzz[i] = tr_x_0_xzzzzz[i] * fe_0 + 5.0 * tr_x_z_xzzzz[i] * fe_0 + tr_x_z_xzzzzz[i] * pa_z[i];

        tr_x_zz_yyyyyy[i] = tr_x_0_yyyyyy[i] * fe_0 + tr_x_z_yyyyyy[i] * pa_z[i];

        tr_x_zz_yyyyyz[i] = tr_x_0_yyyyyz[i] * fe_0 + tr_x_z_yyyyy[i] * fe_0 + tr_x_z_yyyyyz[i] * pa_z[i];

        tr_x_zz_yyyyzz[i] = tr_x_0_yyyyzz[i] * fe_0 + 2.0 * tr_x_z_yyyyz[i] * fe_0 + tr_x_z_yyyyzz[i] * pa_z[i];

        tr_x_zz_yyyzzz[i] = tr_x_0_yyyzzz[i] * fe_0 + 3.0 * tr_x_z_yyyzz[i] * fe_0 + tr_x_z_yyyzzz[i] * pa_z[i];

        tr_x_zz_yyzzzz[i] = tr_x_0_yyzzzz[i] * fe_0 + 4.0 * tr_x_z_yyzzz[i] * fe_0 + tr_x_z_yyzzzz[i] * pa_z[i];

        tr_x_zz_yzzzzz[i] = tr_x_0_yzzzzz[i] * fe_0 + 5.0 * tr_x_z_yzzzz[i] * fe_0 + tr_x_z_yzzzzz[i] * pa_z[i];

        tr_x_zz_zzzzzz[i] = tr_x_0_zzzzzz[i] * fe_0 + 6.0 * tr_x_z_zzzzz[i] * fe_0 + tr_x_z_zzzzzz[i] * pa_z[i];
    }

    // Set up 168-196 components of targeted buffer : DI

    auto tr_y_xx_xxxxxx = pbuffer.data(idx_dip_di + 168);

    auto tr_y_xx_xxxxxy = pbuffer.data(idx_dip_di + 169);

    auto tr_y_xx_xxxxxz = pbuffer.data(idx_dip_di + 170);

    auto tr_y_xx_xxxxyy = pbuffer.data(idx_dip_di + 171);

    auto tr_y_xx_xxxxyz = pbuffer.data(idx_dip_di + 172);

    auto tr_y_xx_xxxxzz = pbuffer.data(idx_dip_di + 173);

    auto tr_y_xx_xxxyyy = pbuffer.data(idx_dip_di + 174);

    auto tr_y_xx_xxxyyz = pbuffer.data(idx_dip_di + 175);

    auto tr_y_xx_xxxyzz = pbuffer.data(idx_dip_di + 176);

    auto tr_y_xx_xxxzzz = pbuffer.data(idx_dip_di + 177);

    auto tr_y_xx_xxyyyy = pbuffer.data(idx_dip_di + 178);

    auto tr_y_xx_xxyyyz = pbuffer.data(idx_dip_di + 179);

    auto tr_y_xx_xxyyzz = pbuffer.data(idx_dip_di + 180);

    auto tr_y_xx_xxyzzz = pbuffer.data(idx_dip_di + 181);

    auto tr_y_xx_xxzzzz = pbuffer.data(idx_dip_di + 182);

    auto tr_y_xx_xyyyyy = pbuffer.data(idx_dip_di + 183);

    auto tr_y_xx_xyyyyz = pbuffer.data(idx_dip_di + 184);

    auto tr_y_xx_xyyyzz = pbuffer.data(idx_dip_di + 185);

    auto tr_y_xx_xyyzzz = pbuffer.data(idx_dip_di + 186);

    auto tr_y_xx_xyzzzz = pbuffer.data(idx_dip_di + 187);

    auto tr_y_xx_xzzzzz = pbuffer.data(idx_dip_di + 188);

    auto tr_y_xx_yyyyyy = pbuffer.data(idx_dip_di + 189);

    auto tr_y_xx_yyyyyz = pbuffer.data(idx_dip_di + 190);

    auto tr_y_xx_yyyyzz = pbuffer.data(idx_dip_di + 191);

    auto tr_y_xx_yyyzzz = pbuffer.data(idx_dip_di + 192);

    auto tr_y_xx_yyzzzz = pbuffer.data(idx_dip_di + 193);

    auto tr_y_xx_yzzzzz = pbuffer.data(idx_dip_di + 194);

    auto tr_y_xx_zzzzzz = pbuffer.data(idx_dip_di + 195);

#pragma omp simd aligned(pa_x,               \
                             tr_y_0_xxxxxx,  \
                             tr_y_0_xxxxxy,  \
                             tr_y_0_xxxxxz,  \
                             tr_y_0_xxxxyy,  \
                             tr_y_0_xxxxyz,  \
                             tr_y_0_xxxxzz,  \
                             tr_y_0_xxxyyy,  \
                             tr_y_0_xxxyyz,  \
                             tr_y_0_xxxyzz,  \
                             tr_y_0_xxxzzz,  \
                             tr_y_0_xxyyyy,  \
                             tr_y_0_xxyyyz,  \
                             tr_y_0_xxyyzz,  \
                             tr_y_0_xxyzzz,  \
                             tr_y_0_xxzzzz,  \
                             tr_y_0_xyyyyy,  \
                             tr_y_0_xyyyyz,  \
                             tr_y_0_xyyyzz,  \
                             tr_y_0_xyyzzz,  \
                             tr_y_0_xyzzzz,  \
                             tr_y_0_xzzzzz,  \
                             tr_y_0_yyyyyy,  \
                             tr_y_0_yyyyyz,  \
                             tr_y_0_yyyyzz,  \
                             tr_y_0_yyyzzz,  \
                             tr_y_0_yyzzzz,  \
                             tr_y_0_yzzzzz,  \
                             tr_y_0_zzzzzz,  \
                             tr_y_x_xxxxx,   \
                             tr_y_x_xxxxxx,  \
                             tr_y_x_xxxxxy,  \
                             tr_y_x_xxxxxz,  \
                             tr_y_x_xxxxy,   \
                             tr_y_x_xxxxyy,  \
                             tr_y_x_xxxxyz,  \
                             tr_y_x_xxxxz,   \
                             tr_y_x_xxxxzz,  \
                             tr_y_x_xxxyy,   \
                             tr_y_x_xxxyyy,  \
                             tr_y_x_xxxyyz,  \
                             tr_y_x_xxxyz,   \
                             tr_y_x_xxxyzz,  \
                             tr_y_x_xxxzz,   \
                             tr_y_x_xxxzzz,  \
                             tr_y_x_xxyyy,   \
                             tr_y_x_xxyyyy,  \
                             tr_y_x_xxyyyz,  \
                             tr_y_x_xxyyz,   \
                             tr_y_x_xxyyzz,  \
                             tr_y_x_xxyzz,   \
                             tr_y_x_xxyzzz,  \
                             tr_y_x_xxzzz,   \
                             tr_y_x_xxzzzz,  \
                             tr_y_x_xyyyy,   \
                             tr_y_x_xyyyyy,  \
                             tr_y_x_xyyyyz,  \
                             tr_y_x_xyyyz,   \
                             tr_y_x_xyyyzz,  \
                             tr_y_x_xyyzz,   \
                             tr_y_x_xyyzzz,  \
                             tr_y_x_xyzzz,   \
                             tr_y_x_xyzzzz,  \
                             tr_y_x_xzzzz,   \
                             tr_y_x_xzzzzz,  \
                             tr_y_x_yyyyy,   \
                             tr_y_x_yyyyyy,  \
                             tr_y_x_yyyyyz,  \
                             tr_y_x_yyyyz,   \
                             tr_y_x_yyyyzz,  \
                             tr_y_x_yyyzz,   \
                             tr_y_x_yyyzzz,  \
                             tr_y_x_yyzzz,   \
                             tr_y_x_yyzzzz,  \
                             tr_y_x_yzzzz,   \
                             tr_y_x_yzzzzz,  \
                             tr_y_x_zzzzz,   \
                             tr_y_x_zzzzzz,  \
                             tr_y_xx_xxxxxx, \
                             tr_y_xx_xxxxxy, \
                             tr_y_xx_xxxxxz, \
                             tr_y_xx_xxxxyy, \
                             tr_y_xx_xxxxyz, \
                             tr_y_xx_xxxxzz, \
                             tr_y_xx_xxxyyy, \
                             tr_y_xx_xxxyyz, \
                             tr_y_xx_xxxyzz, \
                             tr_y_xx_xxxzzz, \
                             tr_y_xx_xxyyyy, \
                             tr_y_xx_xxyyyz, \
                             tr_y_xx_xxyyzz, \
                             tr_y_xx_xxyzzz, \
                             tr_y_xx_xxzzzz, \
                             tr_y_xx_xyyyyy, \
                             tr_y_xx_xyyyyz, \
                             tr_y_xx_xyyyzz, \
                             tr_y_xx_xyyzzz, \
                             tr_y_xx_xyzzzz, \
                             tr_y_xx_xzzzzz, \
                             tr_y_xx_yyyyyy, \
                             tr_y_xx_yyyyyz, \
                             tr_y_xx_yyyyzz, \
                             tr_y_xx_yyyzzz, \
                             tr_y_xx_yyzzzz, \
                             tr_y_xx_yzzzzz, \
                             tr_y_xx_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xx_xxxxxx[i] = tr_y_0_xxxxxx[i] * fe_0 + 6.0 * tr_y_x_xxxxx[i] * fe_0 + tr_y_x_xxxxxx[i] * pa_x[i];

        tr_y_xx_xxxxxy[i] = tr_y_0_xxxxxy[i] * fe_0 + 5.0 * tr_y_x_xxxxy[i] * fe_0 + tr_y_x_xxxxxy[i] * pa_x[i];

        tr_y_xx_xxxxxz[i] = tr_y_0_xxxxxz[i] * fe_0 + 5.0 * tr_y_x_xxxxz[i] * fe_0 + tr_y_x_xxxxxz[i] * pa_x[i];

        tr_y_xx_xxxxyy[i] = tr_y_0_xxxxyy[i] * fe_0 + 4.0 * tr_y_x_xxxyy[i] * fe_0 + tr_y_x_xxxxyy[i] * pa_x[i];

        tr_y_xx_xxxxyz[i] = tr_y_0_xxxxyz[i] * fe_0 + 4.0 * tr_y_x_xxxyz[i] * fe_0 + tr_y_x_xxxxyz[i] * pa_x[i];

        tr_y_xx_xxxxzz[i] = tr_y_0_xxxxzz[i] * fe_0 + 4.0 * tr_y_x_xxxzz[i] * fe_0 + tr_y_x_xxxxzz[i] * pa_x[i];

        tr_y_xx_xxxyyy[i] = tr_y_0_xxxyyy[i] * fe_0 + 3.0 * tr_y_x_xxyyy[i] * fe_0 + tr_y_x_xxxyyy[i] * pa_x[i];

        tr_y_xx_xxxyyz[i] = tr_y_0_xxxyyz[i] * fe_0 + 3.0 * tr_y_x_xxyyz[i] * fe_0 + tr_y_x_xxxyyz[i] * pa_x[i];

        tr_y_xx_xxxyzz[i] = tr_y_0_xxxyzz[i] * fe_0 + 3.0 * tr_y_x_xxyzz[i] * fe_0 + tr_y_x_xxxyzz[i] * pa_x[i];

        tr_y_xx_xxxzzz[i] = tr_y_0_xxxzzz[i] * fe_0 + 3.0 * tr_y_x_xxzzz[i] * fe_0 + tr_y_x_xxxzzz[i] * pa_x[i];

        tr_y_xx_xxyyyy[i] = tr_y_0_xxyyyy[i] * fe_0 + 2.0 * tr_y_x_xyyyy[i] * fe_0 + tr_y_x_xxyyyy[i] * pa_x[i];

        tr_y_xx_xxyyyz[i] = tr_y_0_xxyyyz[i] * fe_0 + 2.0 * tr_y_x_xyyyz[i] * fe_0 + tr_y_x_xxyyyz[i] * pa_x[i];

        tr_y_xx_xxyyzz[i] = tr_y_0_xxyyzz[i] * fe_0 + 2.0 * tr_y_x_xyyzz[i] * fe_0 + tr_y_x_xxyyzz[i] * pa_x[i];

        tr_y_xx_xxyzzz[i] = tr_y_0_xxyzzz[i] * fe_0 + 2.0 * tr_y_x_xyzzz[i] * fe_0 + tr_y_x_xxyzzz[i] * pa_x[i];

        tr_y_xx_xxzzzz[i] = tr_y_0_xxzzzz[i] * fe_0 + 2.0 * tr_y_x_xzzzz[i] * fe_0 + tr_y_x_xxzzzz[i] * pa_x[i];

        tr_y_xx_xyyyyy[i] = tr_y_0_xyyyyy[i] * fe_0 + tr_y_x_yyyyy[i] * fe_0 + tr_y_x_xyyyyy[i] * pa_x[i];

        tr_y_xx_xyyyyz[i] = tr_y_0_xyyyyz[i] * fe_0 + tr_y_x_yyyyz[i] * fe_0 + tr_y_x_xyyyyz[i] * pa_x[i];

        tr_y_xx_xyyyzz[i] = tr_y_0_xyyyzz[i] * fe_0 + tr_y_x_yyyzz[i] * fe_0 + tr_y_x_xyyyzz[i] * pa_x[i];

        tr_y_xx_xyyzzz[i] = tr_y_0_xyyzzz[i] * fe_0 + tr_y_x_yyzzz[i] * fe_0 + tr_y_x_xyyzzz[i] * pa_x[i];

        tr_y_xx_xyzzzz[i] = tr_y_0_xyzzzz[i] * fe_0 + tr_y_x_yzzzz[i] * fe_0 + tr_y_x_xyzzzz[i] * pa_x[i];

        tr_y_xx_xzzzzz[i] = tr_y_0_xzzzzz[i] * fe_0 + tr_y_x_zzzzz[i] * fe_0 + tr_y_x_xzzzzz[i] * pa_x[i];

        tr_y_xx_yyyyyy[i] = tr_y_0_yyyyyy[i] * fe_0 + tr_y_x_yyyyyy[i] * pa_x[i];

        tr_y_xx_yyyyyz[i] = tr_y_0_yyyyyz[i] * fe_0 + tr_y_x_yyyyyz[i] * pa_x[i];

        tr_y_xx_yyyyzz[i] = tr_y_0_yyyyzz[i] * fe_0 + tr_y_x_yyyyzz[i] * pa_x[i];

        tr_y_xx_yyyzzz[i] = tr_y_0_yyyzzz[i] * fe_0 + tr_y_x_yyyzzz[i] * pa_x[i];

        tr_y_xx_yyzzzz[i] = tr_y_0_yyzzzz[i] * fe_0 + tr_y_x_yyzzzz[i] * pa_x[i];

        tr_y_xx_yzzzzz[i] = tr_y_0_yzzzzz[i] * fe_0 + tr_y_x_yzzzzz[i] * pa_x[i];

        tr_y_xx_zzzzzz[i] = tr_y_0_zzzzzz[i] * fe_0 + tr_y_x_zzzzzz[i] * pa_x[i];
    }

    // Set up 196-224 components of targeted buffer : DI

    auto tr_y_xy_xxxxxx = pbuffer.data(idx_dip_di + 196);

    auto tr_y_xy_xxxxxy = pbuffer.data(idx_dip_di + 197);

    auto tr_y_xy_xxxxxz = pbuffer.data(idx_dip_di + 198);

    auto tr_y_xy_xxxxyy = pbuffer.data(idx_dip_di + 199);

    auto tr_y_xy_xxxxyz = pbuffer.data(idx_dip_di + 200);

    auto tr_y_xy_xxxxzz = pbuffer.data(idx_dip_di + 201);

    auto tr_y_xy_xxxyyy = pbuffer.data(idx_dip_di + 202);

    auto tr_y_xy_xxxyyz = pbuffer.data(idx_dip_di + 203);

    auto tr_y_xy_xxxyzz = pbuffer.data(idx_dip_di + 204);

    auto tr_y_xy_xxxzzz = pbuffer.data(idx_dip_di + 205);

    auto tr_y_xy_xxyyyy = pbuffer.data(idx_dip_di + 206);

    auto tr_y_xy_xxyyyz = pbuffer.data(idx_dip_di + 207);

    auto tr_y_xy_xxyyzz = pbuffer.data(idx_dip_di + 208);

    auto tr_y_xy_xxyzzz = pbuffer.data(idx_dip_di + 209);

    auto tr_y_xy_xxzzzz = pbuffer.data(idx_dip_di + 210);

    auto tr_y_xy_xyyyyy = pbuffer.data(idx_dip_di + 211);

    auto tr_y_xy_xyyyyz = pbuffer.data(idx_dip_di + 212);

    auto tr_y_xy_xyyyzz = pbuffer.data(idx_dip_di + 213);

    auto tr_y_xy_xyyzzz = pbuffer.data(idx_dip_di + 214);

    auto tr_y_xy_xyzzzz = pbuffer.data(idx_dip_di + 215);

    auto tr_y_xy_xzzzzz = pbuffer.data(idx_dip_di + 216);

    auto tr_y_xy_yyyyyy = pbuffer.data(idx_dip_di + 217);

    auto tr_y_xy_yyyyyz = pbuffer.data(idx_dip_di + 218);

    auto tr_y_xy_yyyyzz = pbuffer.data(idx_dip_di + 219);

    auto tr_y_xy_yyyzzz = pbuffer.data(idx_dip_di + 220);

    auto tr_y_xy_yyzzzz = pbuffer.data(idx_dip_di + 221);

    auto tr_y_xy_yzzzzz = pbuffer.data(idx_dip_di + 222);

    auto tr_y_xy_zzzzzz = pbuffer.data(idx_dip_di + 223);

#pragma omp simd aligned(pa_x,               \
                             tr_y_xy_xxxxxx, \
                             tr_y_xy_xxxxxy, \
                             tr_y_xy_xxxxxz, \
                             tr_y_xy_xxxxyy, \
                             tr_y_xy_xxxxyz, \
                             tr_y_xy_xxxxzz, \
                             tr_y_xy_xxxyyy, \
                             tr_y_xy_xxxyyz, \
                             tr_y_xy_xxxyzz, \
                             tr_y_xy_xxxzzz, \
                             tr_y_xy_xxyyyy, \
                             tr_y_xy_xxyyyz, \
                             tr_y_xy_xxyyzz, \
                             tr_y_xy_xxyzzz, \
                             tr_y_xy_xxzzzz, \
                             tr_y_xy_xyyyyy, \
                             tr_y_xy_xyyyyz, \
                             tr_y_xy_xyyyzz, \
                             tr_y_xy_xyyzzz, \
                             tr_y_xy_xyzzzz, \
                             tr_y_xy_xzzzzz, \
                             tr_y_xy_yyyyyy, \
                             tr_y_xy_yyyyyz, \
                             tr_y_xy_yyyyzz, \
                             tr_y_xy_yyyzzz, \
                             tr_y_xy_yyzzzz, \
                             tr_y_xy_yzzzzz, \
                             tr_y_xy_zzzzzz, \
                             tr_y_y_xxxxx,   \
                             tr_y_y_xxxxxx,  \
                             tr_y_y_xxxxxy,  \
                             tr_y_y_xxxxxz,  \
                             tr_y_y_xxxxy,   \
                             tr_y_y_xxxxyy,  \
                             tr_y_y_xxxxyz,  \
                             tr_y_y_xxxxz,   \
                             tr_y_y_xxxxzz,  \
                             tr_y_y_xxxyy,   \
                             tr_y_y_xxxyyy,  \
                             tr_y_y_xxxyyz,  \
                             tr_y_y_xxxyz,   \
                             tr_y_y_xxxyzz,  \
                             tr_y_y_xxxzz,   \
                             tr_y_y_xxxzzz,  \
                             tr_y_y_xxyyy,   \
                             tr_y_y_xxyyyy,  \
                             tr_y_y_xxyyyz,  \
                             tr_y_y_xxyyz,   \
                             tr_y_y_xxyyzz,  \
                             tr_y_y_xxyzz,   \
                             tr_y_y_xxyzzz,  \
                             tr_y_y_xxzzz,   \
                             tr_y_y_xxzzzz,  \
                             tr_y_y_xyyyy,   \
                             tr_y_y_xyyyyy,  \
                             tr_y_y_xyyyyz,  \
                             tr_y_y_xyyyz,   \
                             tr_y_y_xyyyzz,  \
                             tr_y_y_xyyzz,   \
                             tr_y_y_xyyzzz,  \
                             tr_y_y_xyzzz,   \
                             tr_y_y_xyzzzz,  \
                             tr_y_y_xzzzz,   \
                             tr_y_y_xzzzzz,  \
                             tr_y_y_yyyyy,   \
                             tr_y_y_yyyyyy,  \
                             tr_y_y_yyyyyz,  \
                             tr_y_y_yyyyz,   \
                             tr_y_y_yyyyzz,  \
                             tr_y_y_yyyzz,   \
                             tr_y_y_yyyzzz,  \
                             tr_y_y_yyzzz,   \
                             tr_y_y_yyzzzz,  \
                             tr_y_y_yzzzz,   \
                             tr_y_y_yzzzzz,  \
                             tr_y_y_zzzzz,   \
                             tr_y_y_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xy_xxxxxx[i] = 6.0 * tr_y_y_xxxxx[i] * fe_0 + tr_y_y_xxxxxx[i] * pa_x[i];

        tr_y_xy_xxxxxy[i] = 5.0 * tr_y_y_xxxxy[i] * fe_0 + tr_y_y_xxxxxy[i] * pa_x[i];

        tr_y_xy_xxxxxz[i] = 5.0 * tr_y_y_xxxxz[i] * fe_0 + tr_y_y_xxxxxz[i] * pa_x[i];

        tr_y_xy_xxxxyy[i] = 4.0 * tr_y_y_xxxyy[i] * fe_0 + tr_y_y_xxxxyy[i] * pa_x[i];

        tr_y_xy_xxxxyz[i] = 4.0 * tr_y_y_xxxyz[i] * fe_0 + tr_y_y_xxxxyz[i] * pa_x[i];

        tr_y_xy_xxxxzz[i] = 4.0 * tr_y_y_xxxzz[i] * fe_0 + tr_y_y_xxxxzz[i] * pa_x[i];

        tr_y_xy_xxxyyy[i] = 3.0 * tr_y_y_xxyyy[i] * fe_0 + tr_y_y_xxxyyy[i] * pa_x[i];

        tr_y_xy_xxxyyz[i] = 3.0 * tr_y_y_xxyyz[i] * fe_0 + tr_y_y_xxxyyz[i] * pa_x[i];

        tr_y_xy_xxxyzz[i] = 3.0 * tr_y_y_xxyzz[i] * fe_0 + tr_y_y_xxxyzz[i] * pa_x[i];

        tr_y_xy_xxxzzz[i] = 3.0 * tr_y_y_xxzzz[i] * fe_0 + tr_y_y_xxxzzz[i] * pa_x[i];

        tr_y_xy_xxyyyy[i] = 2.0 * tr_y_y_xyyyy[i] * fe_0 + tr_y_y_xxyyyy[i] * pa_x[i];

        tr_y_xy_xxyyyz[i] = 2.0 * tr_y_y_xyyyz[i] * fe_0 + tr_y_y_xxyyyz[i] * pa_x[i];

        tr_y_xy_xxyyzz[i] = 2.0 * tr_y_y_xyyzz[i] * fe_0 + tr_y_y_xxyyzz[i] * pa_x[i];

        tr_y_xy_xxyzzz[i] = 2.0 * tr_y_y_xyzzz[i] * fe_0 + tr_y_y_xxyzzz[i] * pa_x[i];

        tr_y_xy_xxzzzz[i] = 2.0 * tr_y_y_xzzzz[i] * fe_0 + tr_y_y_xxzzzz[i] * pa_x[i];

        tr_y_xy_xyyyyy[i] = tr_y_y_yyyyy[i] * fe_0 + tr_y_y_xyyyyy[i] * pa_x[i];

        tr_y_xy_xyyyyz[i] = tr_y_y_yyyyz[i] * fe_0 + tr_y_y_xyyyyz[i] * pa_x[i];

        tr_y_xy_xyyyzz[i] = tr_y_y_yyyzz[i] * fe_0 + tr_y_y_xyyyzz[i] * pa_x[i];

        tr_y_xy_xyyzzz[i] = tr_y_y_yyzzz[i] * fe_0 + tr_y_y_xyyzzz[i] * pa_x[i];

        tr_y_xy_xyzzzz[i] = tr_y_y_yzzzz[i] * fe_0 + tr_y_y_xyzzzz[i] * pa_x[i];

        tr_y_xy_xzzzzz[i] = tr_y_y_zzzzz[i] * fe_0 + tr_y_y_xzzzzz[i] * pa_x[i];

        tr_y_xy_yyyyyy[i] = tr_y_y_yyyyyy[i] * pa_x[i];

        tr_y_xy_yyyyyz[i] = tr_y_y_yyyyyz[i] * pa_x[i];

        tr_y_xy_yyyyzz[i] = tr_y_y_yyyyzz[i] * pa_x[i];

        tr_y_xy_yyyzzz[i] = tr_y_y_yyyzzz[i] * pa_x[i];

        tr_y_xy_yyzzzz[i] = tr_y_y_yyzzzz[i] * pa_x[i];

        tr_y_xy_yzzzzz[i] = tr_y_y_yzzzzz[i] * pa_x[i];

        tr_y_xy_zzzzzz[i] = tr_y_y_zzzzzz[i] * pa_x[i];
    }

    // Set up 224-252 components of targeted buffer : DI

    auto tr_y_xz_xxxxxx = pbuffer.data(idx_dip_di + 224);

    auto tr_y_xz_xxxxxy = pbuffer.data(idx_dip_di + 225);

    auto tr_y_xz_xxxxxz = pbuffer.data(idx_dip_di + 226);

    auto tr_y_xz_xxxxyy = pbuffer.data(idx_dip_di + 227);

    auto tr_y_xz_xxxxyz = pbuffer.data(idx_dip_di + 228);

    auto tr_y_xz_xxxxzz = pbuffer.data(idx_dip_di + 229);

    auto tr_y_xz_xxxyyy = pbuffer.data(idx_dip_di + 230);

    auto tr_y_xz_xxxyyz = pbuffer.data(idx_dip_di + 231);

    auto tr_y_xz_xxxyzz = pbuffer.data(idx_dip_di + 232);

    auto tr_y_xz_xxxzzz = pbuffer.data(idx_dip_di + 233);

    auto tr_y_xz_xxyyyy = pbuffer.data(idx_dip_di + 234);

    auto tr_y_xz_xxyyyz = pbuffer.data(idx_dip_di + 235);

    auto tr_y_xz_xxyyzz = pbuffer.data(idx_dip_di + 236);

    auto tr_y_xz_xxyzzz = pbuffer.data(idx_dip_di + 237);

    auto tr_y_xz_xxzzzz = pbuffer.data(idx_dip_di + 238);

    auto tr_y_xz_xyyyyy = pbuffer.data(idx_dip_di + 239);

    auto tr_y_xz_xyyyyz = pbuffer.data(idx_dip_di + 240);

    auto tr_y_xz_xyyyzz = pbuffer.data(idx_dip_di + 241);

    auto tr_y_xz_xyyzzz = pbuffer.data(idx_dip_di + 242);

    auto tr_y_xz_xyzzzz = pbuffer.data(idx_dip_di + 243);

    auto tr_y_xz_xzzzzz = pbuffer.data(idx_dip_di + 244);

    auto tr_y_xz_yyyyyy = pbuffer.data(idx_dip_di + 245);

    auto tr_y_xz_yyyyyz = pbuffer.data(idx_dip_di + 246);

    auto tr_y_xz_yyyyzz = pbuffer.data(idx_dip_di + 247);

    auto tr_y_xz_yyyzzz = pbuffer.data(idx_dip_di + 248);

    auto tr_y_xz_yyzzzz = pbuffer.data(idx_dip_di + 249);

    auto tr_y_xz_yzzzzz = pbuffer.data(idx_dip_di + 250);

    auto tr_y_xz_zzzzzz = pbuffer.data(idx_dip_di + 251);

#pragma omp simd aligned(pa_x,               \
                             pa_z,           \
                             tr_y_x_xxxxxx,  \
                             tr_y_x_xxxxxy,  \
                             tr_y_x_xxxxyy,  \
                             tr_y_x_xxxyyy,  \
                             tr_y_x_xxyyyy,  \
                             tr_y_x_xyyyyy,  \
                             tr_y_xz_xxxxxx, \
                             tr_y_xz_xxxxxy, \
                             tr_y_xz_xxxxxz, \
                             tr_y_xz_xxxxyy, \
                             tr_y_xz_xxxxyz, \
                             tr_y_xz_xxxxzz, \
                             tr_y_xz_xxxyyy, \
                             tr_y_xz_xxxyyz, \
                             tr_y_xz_xxxyzz, \
                             tr_y_xz_xxxzzz, \
                             tr_y_xz_xxyyyy, \
                             tr_y_xz_xxyyyz, \
                             tr_y_xz_xxyyzz, \
                             tr_y_xz_xxyzzz, \
                             tr_y_xz_xxzzzz, \
                             tr_y_xz_xyyyyy, \
                             tr_y_xz_xyyyyz, \
                             tr_y_xz_xyyyzz, \
                             tr_y_xz_xyyzzz, \
                             tr_y_xz_xyzzzz, \
                             tr_y_xz_xzzzzz, \
                             tr_y_xz_yyyyyy, \
                             tr_y_xz_yyyyyz, \
                             tr_y_xz_yyyyzz, \
                             tr_y_xz_yyyzzz, \
                             tr_y_xz_yyzzzz, \
                             tr_y_xz_yzzzzz, \
                             tr_y_xz_zzzzzz, \
                             tr_y_z_xxxxxz,  \
                             tr_y_z_xxxxyz,  \
                             tr_y_z_xxxxz,   \
                             tr_y_z_xxxxzz,  \
                             tr_y_z_xxxyyz,  \
                             tr_y_z_xxxyz,   \
                             tr_y_z_xxxyzz,  \
                             tr_y_z_xxxzz,   \
                             tr_y_z_xxxzzz,  \
                             tr_y_z_xxyyyz,  \
                             tr_y_z_xxyyz,   \
                             tr_y_z_xxyyzz,  \
                             tr_y_z_xxyzz,   \
                             tr_y_z_xxyzzz,  \
                             tr_y_z_xxzzz,   \
                             tr_y_z_xxzzzz,  \
                             tr_y_z_xyyyyz,  \
                             tr_y_z_xyyyz,   \
                             tr_y_z_xyyyzz,  \
                             tr_y_z_xyyzz,   \
                             tr_y_z_xyyzzz,  \
                             tr_y_z_xyzzz,   \
                             tr_y_z_xyzzzz,  \
                             tr_y_z_xzzzz,   \
                             tr_y_z_xzzzzz,  \
                             tr_y_z_yyyyyy,  \
                             tr_y_z_yyyyyz,  \
                             tr_y_z_yyyyz,   \
                             tr_y_z_yyyyzz,  \
                             tr_y_z_yyyzz,   \
                             tr_y_z_yyyzzz,  \
                             tr_y_z_yyzzz,   \
                             tr_y_z_yyzzzz,  \
                             tr_y_z_yzzzz,   \
                             tr_y_z_yzzzzz,  \
                             tr_y_z_zzzzz,   \
                             tr_y_z_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_xz_xxxxxx[i] = tr_y_x_xxxxxx[i] * pa_z[i];

        tr_y_xz_xxxxxy[i] = tr_y_x_xxxxxy[i] * pa_z[i];

        tr_y_xz_xxxxxz[i] = 5.0 * tr_y_z_xxxxz[i] * fe_0 + tr_y_z_xxxxxz[i] * pa_x[i];

        tr_y_xz_xxxxyy[i] = tr_y_x_xxxxyy[i] * pa_z[i];

        tr_y_xz_xxxxyz[i] = 4.0 * tr_y_z_xxxyz[i] * fe_0 + tr_y_z_xxxxyz[i] * pa_x[i];

        tr_y_xz_xxxxzz[i] = 4.0 * tr_y_z_xxxzz[i] * fe_0 + tr_y_z_xxxxzz[i] * pa_x[i];

        tr_y_xz_xxxyyy[i] = tr_y_x_xxxyyy[i] * pa_z[i];

        tr_y_xz_xxxyyz[i] = 3.0 * tr_y_z_xxyyz[i] * fe_0 + tr_y_z_xxxyyz[i] * pa_x[i];

        tr_y_xz_xxxyzz[i] = 3.0 * tr_y_z_xxyzz[i] * fe_0 + tr_y_z_xxxyzz[i] * pa_x[i];

        tr_y_xz_xxxzzz[i] = 3.0 * tr_y_z_xxzzz[i] * fe_0 + tr_y_z_xxxzzz[i] * pa_x[i];

        tr_y_xz_xxyyyy[i] = tr_y_x_xxyyyy[i] * pa_z[i];

        tr_y_xz_xxyyyz[i] = 2.0 * tr_y_z_xyyyz[i] * fe_0 + tr_y_z_xxyyyz[i] * pa_x[i];

        tr_y_xz_xxyyzz[i] = 2.0 * tr_y_z_xyyzz[i] * fe_0 + tr_y_z_xxyyzz[i] * pa_x[i];

        tr_y_xz_xxyzzz[i] = 2.0 * tr_y_z_xyzzz[i] * fe_0 + tr_y_z_xxyzzz[i] * pa_x[i];

        tr_y_xz_xxzzzz[i] = 2.0 * tr_y_z_xzzzz[i] * fe_0 + tr_y_z_xxzzzz[i] * pa_x[i];

        tr_y_xz_xyyyyy[i] = tr_y_x_xyyyyy[i] * pa_z[i];

        tr_y_xz_xyyyyz[i] = tr_y_z_yyyyz[i] * fe_0 + tr_y_z_xyyyyz[i] * pa_x[i];

        tr_y_xz_xyyyzz[i] = tr_y_z_yyyzz[i] * fe_0 + tr_y_z_xyyyzz[i] * pa_x[i];

        tr_y_xz_xyyzzz[i] = tr_y_z_yyzzz[i] * fe_0 + tr_y_z_xyyzzz[i] * pa_x[i];

        tr_y_xz_xyzzzz[i] = tr_y_z_yzzzz[i] * fe_0 + tr_y_z_xyzzzz[i] * pa_x[i];

        tr_y_xz_xzzzzz[i] = tr_y_z_zzzzz[i] * fe_0 + tr_y_z_xzzzzz[i] * pa_x[i];

        tr_y_xz_yyyyyy[i] = tr_y_z_yyyyyy[i] * pa_x[i];

        tr_y_xz_yyyyyz[i] = tr_y_z_yyyyyz[i] * pa_x[i];

        tr_y_xz_yyyyzz[i] = tr_y_z_yyyyzz[i] * pa_x[i];

        tr_y_xz_yyyzzz[i] = tr_y_z_yyyzzz[i] * pa_x[i];

        tr_y_xz_yyzzzz[i] = tr_y_z_yyzzzz[i] * pa_x[i];

        tr_y_xz_yzzzzz[i] = tr_y_z_yzzzzz[i] * pa_x[i];

        tr_y_xz_zzzzzz[i] = tr_y_z_zzzzzz[i] * pa_x[i];
    }

    // Set up 252-280 components of targeted buffer : DI

    auto tr_y_yy_xxxxxx = pbuffer.data(idx_dip_di + 252);

    auto tr_y_yy_xxxxxy = pbuffer.data(idx_dip_di + 253);

    auto tr_y_yy_xxxxxz = pbuffer.data(idx_dip_di + 254);

    auto tr_y_yy_xxxxyy = pbuffer.data(idx_dip_di + 255);

    auto tr_y_yy_xxxxyz = pbuffer.data(idx_dip_di + 256);

    auto tr_y_yy_xxxxzz = pbuffer.data(idx_dip_di + 257);

    auto tr_y_yy_xxxyyy = pbuffer.data(idx_dip_di + 258);

    auto tr_y_yy_xxxyyz = pbuffer.data(idx_dip_di + 259);

    auto tr_y_yy_xxxyzz = pbuffer.data(idx_dip_di + 260);

    auto tr_y_yy_xxxzzz = pbuffer.data(idx_dip_di + 261);

    auto tr_y_yy_xxyyyy = pbuffer.data(idx_dip_di + 262);

    auto tr_y_yy_xxyyyz = pbuffer.data(idx_dip_di + 263);

    auto tr_y_yy_xxyyzz = pbuffer.data(idx_dip_di + 264);

    auto tr_y_yy_xxyzzz = pbuffer.data(idx_dip_di + 265);

    auto tr_y_yy_xxzzzz = pbuffer.data(idx_dip_di + 266);

    auto tr_y_yy_xyyyyy = pbuffer.data(idx_dip_di + 267);

    auto tr_y_yy_xyyyyz = pbuffer.data(idx_dip_di + 268);

    auto tr_y_yy_xyyyzz = pbuffer.data(idx_dip_di + 269);

    auto tr_y_yy_xyyzzz = pbuffer.data(idx_dip_di + 270);

    auto tr_y_yy_xyzzzz = pbuffer.data(idx_dip_di + 271);

    auto tr_y_yy_xzzzzz = pbuffer.data(idx_dip_di + 272);

    auto tr_y_yy_yyyyyy = pbuffer.data(idx_dip_di + 273);

    auto tr_y_yy_yyyyyz = pbuffer.data(idx_dip_di + 274);

    auto tr_y_yy_yyyyzz = pbuffer.data(idx_dip_di + 275);

    auto tr_y_yy_yyyzzz = pbuffer.data(idx_dip_di + 276);

    auto tr_y_yy_yyzzzz = pbuffer.data(idx_dip_di + 277);

    auto tr_y_yy_yzzzzz = pbuffer.data(idx_dip_di + 278);

    auto tr_y_yy_zzzzzz = pbuffer.data(idx_dip_di + 279);

#pragma omp simd aligned(pa_y,               \
                             tr_y_0_xxxxxx,  \
                             tr_y_0_xxxxxy,  \
                             tr_y_0_xxxxxz,  \
                             tr_y_0_xxxxyy,  \
                             tr_y_0_xxxxyz,  \
                             tr_y_0_xxxxzz,  \
                             tr_y_0_xxxyyy,  \
                             tr_y_0_xxxyyz,  \
                             tr_y_0_xxxyzz,  \
                             tr_y_0_xxxzzz,  \
                             tr_y_0_xxyyyy,  \
                             tr_y_0_xxyyyz,  \
                             tr_y_0_xxyyzz,  \
                             tr_y_0_xxyzzz,  \
                             tr_y_0_xxzzzz,  \
                             tr_y_0_xyyyyy,  \
                             tr_y_0_xyyyyz,  \
                             tr_y_0_xyyyzz,  \
                             tr_y_0_xyyzzz,  \
                             tr_y_0_xyzzzz,  \
                             tr_y_0_xzzzzz,  \
                             tr_y_0_yyyyyy,  \
                             tr_y_0_yyyyyz,  \
                             tr_y_0_yyyyzz,  \
                             tr_y_0_yyyzzz,  \
                             tr_y_0_yyzzzz,  \
                             tr_y_0_yzzzzz,  \
                             tr_y_0_zzzzzz,  \
                             tr_y_y_xxxxx,   \
                             tr_y_y_xxxxxx,  \
                             tr_y_y_xxxxxy,  \
                             tr_y_y_xxxxxz,  \
                             tr_y_y_xxxxy,   \
                             tr_y_y_xxxxyy,  \
                             tr_y_y_xxxxyz,  \
                             tr_y_y_xxxxz,   \
                             tr_y_y_xxxxzz,  \
                             tr_y_y_xxxyy,   \
                             tr_y_y_xxxyyy,  \
                             tr_y_y_xxxyyz,  \
                             tr_y_y_xxxyz,   \
                             tr_y_y_xxxyzz,  \
                             tr_y_y_xxxzz,   \
                             tr_y_y_xxxzzz,  \
                             tr_y_y_xxyyy,   \
                             tr_y_y_xxyyyy,  \
                             tr_y_y_xxyyyz,  \
                             tr_y_y_xxyyz,   \
                             tr_y_y_xxyyzz,  \
                             tr_y_y_xxyzz,   \
                             tr_y_y_xxyzzz,  \
                             tr_y_y_xxzzz,   \
                             tr_y_y_xxzzzz,  \
                             tr_y_y_xyyyy,   \
                             tr_y_y_xyyyyy,  \
                             tr_y_y_xyyyyz,  \
                             tr_y_y_xyyyz,   \
                             tr_y_y_xyyyzz,  \
                             tr_y_y_xyyzz,   \
                             tr_y_y_xyyzzz,  \
                             tr_y_y_xyzzz,   \
                             tr_y_y_xyzzzz,  \
                             tr_y_y_xzzzz,   \
                             tr_y_y_xzzzzz,  \
                             tr_y_y_yyyyy,   \
                             tr_y_y_yyyyyy,  \
                             tr_y_y_yyyyyz,  \
                             tr_y_y_yyyyz,   \
                             tr_y_y_yyyyzz,  \
                             tr_y_y_yyyzz,   \
                             tr_y_y_yyyzzz,  \
                             tr_y_y_yyzzz,   \
                             tr_y_y_yyzzzz,  \
                             tr_y_y_yzzzz,   \
                             tr_y_y_yzzzzz,  \
                             tr_y_y_zzzzz,   \
                             tr_y_y_zzzzzz,  \
                             tr_y_yy_xxxxxx, \
                             tr_y_yy_xxxxxy, \
                             tr_y_yy_xxxxxz, \
                             tr_y_yy_xxxxyy, \
                             tr_y_yy_xxxxyz, \
                             tr_y_yy_xxxxzz, \
                             tr_y_yy_xxxyyy, \
                             tr_y_yy_xxxyyz, \
                             tr_y_yy_xxxyzz, \
                             tr_y_yy_xxxzzz, \
                             tr_y_yy_xxyyyy, \
                             tr_y_yy_xxyyyz, \
                             tr_y_yy_xxyyzz, \
                             tr_y_yy_xxyzzz, \
                             tr_y_yy_xxzzzz, \
                             tr_y_yy_xyyyyy, \
                             tr_y_yy_xyyyyz, \
                             tr_y_yy_xyyyzz, \
                             tr_y_yy_xyyzzz, \
                             tr_y_yy_xyzzzz, \
                             tr_y_yy_xzzzzz, \
                             tr_y_yy_yyyyyy, \
                             tr_y_yy_yyyyyz, \
                             tr_y_yy_yyyyzz, \
                             tr_y_yy_yyyzzz, \
                             tr_y_yy_yyzzzz, \
                             tr_y_yy_yzzzzz, \
                             tr_y_yy_zzzzzz, \
                             ts_y_xxxxxx,    \
                             ts_y_xxxxxy,    \
                             ts_y_xxxxxz,    \
                             ts_y_xxxxyy,    \
                             ts_y_xxxxyz,    \
                             ts_y_xxxxzz,    \
                             ts_y_xxxyyy,    \
                             ts_y_xxxyyz,    \
                             ts_y_xxxyzz,    \
                             ts_y_xxxzzz,    \
                             ts_y_xxyyyy,    \
                             ts_y_xxyyyz,    \
                             ts_y_xxyyzz,    \
                             ts_y_xxyzzz,    \
                             ts_y_xxzzzz,    \
                             ts_y_xyyyyy,    \
                             ts_y_xyyyyz,    \
                             ts_y_xyyyzz,    \
                             ts_y_xyyzzz,    \
                             ts_y_xyzzzz,    \
                             ts_y_xzzzzz,    \
                             ts_y_yyyyyy,    \
                             ts_y_yyyyyz,    \
                             ts_y_yyyyzz,    \
                             ts_y_yyyzzz,    \
                             ts_y_yyzzzz,    \
                             ts_y_yzzzzz,    \
                             ts_y_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yy_xxxxxx[i] = tr_y_0_xxxxxx[i] * fe_0 + ts_y_xxxxxx[i] * fe_0 + tr_y_y_xxxxxx[i] * pa_y[i];

        tr_y_yy_xxxxxy[i] = tr_y_0_xxxxxy[i] * fe_0 + tr_y_y_xxxxx[i] * fe_0 + ts_y_xxxxxy[i] * fe_0 + tr_y_y_xxxxxy[i] * pa_y[i];

        tr_y_yy_xxxxxz[i] = tr_y_0_xxxxxz[i] * fe_0 + ts_y_xxxxxz[i] * fe_0 + tr_y_y_xxxxxz[i] * pa_y[i];

        tr_y_yy_xxxxyy[i] = tr_y_0_xxxxyy[i] * fe_0 + 2.0 * tr_y_y_xxxxy[i] * fe_0 + ts_y_xxxxyy[i] * fe_0 + tr_y_y_xxxxyy[i] * pa_y[i];

        tr_y_yy_xxxxyz[i] = tr_y_0_xxxxyz[i] * fe_0 + tr_y_y_xxxxz[i] * fe_0 + ts_y_xxxxyz[i] * fe_0 + tr_y_y_xxxxyz[i] * pa_y[i];

        tr_y_yy_xxxxzz[i] = tr_y_0_xxxxzz[i] * fe_0 + ts_y_xxxxzz[i] * fe_0 + tr_y_y_xxxxzz[i] * pa_y[i];

        tr_y_yy_xxxyyy[i] = tr_y_0_xxxyyy[i] * fe_0 + 3.0 * tr_y_y_xxxyy[i] * fe_0 + ts_y_xxxyyy[i] * fe_0 + tr_y_y_xxxyyy[i] * pa_y[i];

        tr_y_yy_xxxyyz[i] = tr_y_0_xxxyyz[i] * fe_0 + 2.0 * tr_y_y_xxxyz[i] * fe_0 + ts_y_xxxyyz[i] * fe_0 + tr_y_y_xxxyyz[i] * pa_y[i];

        tr_y_yy_xxxyzz[i] = tr_y_0_xxxyzz[i] * fe_0 + tr_y_y_xxxzz[i] * fe_0 + ts_y_xxxyzz[i] * fe_0 + tr_y_y_xxxyzz[i] * pa_y[i];

        tr_y_yy_xxxzzz[i] = tr_y_0_xxxzzz[i] * fe_0 + ts_y_xxxzzz[i] * fe_0 + tr_y_y_xxxzzz[i] * pa_y[i];

        tr_y_yy_xxyyyy[i] = tr_y_0_xxyyyy[i] * fe_0 + 4.0 * tr_y_y_xxyyy[i] * fe_0 + ts_y_xxyyyy[i] * fe_0 + tr_y_y_xxyyyy[i] * pa_y[i];

        tr_y_yy_xxyyyz[i] = tr_y_0_xxyyyz[i] * fe_0 + 3.0 * tr_y_y_xxyyz[i] * fe_0 + ts_y_xxyyyz[i] * fe_0 + tr_y_y_xxyyyz[i] * pa_y[i];

        tr_y_yy_xxyyzz[i] = tr_y_0_xxyyzz[i] * fe_0 + 2.0 * tr_y_y_xxyzz[i] * fe_0 + ts_y_xxyyzz[i] * fe_0 + tr_y_y_xxyyzz[i] * pa_y[i];

        tr_y_yy_xxyzzz[i] = tr_y_0_xxyzzz[i] * fe_0 + tr_y_y_xxzzz[i] * fe_0 + ts_y_xxyzzz[i] * fe_0 + tr_y_y_xxyzzz[i] * pa_y[i];

        tr_y_yy_xxzzzz[i] = tr_y_0_xxzzzz[i] * fe_0 + ts_y_xxzzzz[i] * fe_0 + tr_y_y_xxzzzz[i] * pa_y[i];

        tr_y_yy_xyyyyy[i] = tr_y_0_xyyyyy[i] * fe_0 + 5.0 * tr_y_y_xyyyy[i] * fe_0 + ts_y_xyyyyy[i] * fe_0 + tr_y_y_xyyyyy[i] * pa_y[i];

        tr_y_yy_xyyyyz[i] = tr_y_0_xyyyyz[i] * fe_0 + 4.0 * tr_y_y_xyyyz[i] * fe_0 + ts_y_xyyyyz[i] * fe_0 + tr_y_y_xyyyyz[i] * pa_y[i];

        tr_y_yy_xyyyzz[i] = tr_y_0_xyyyzz[i] * fe_0 + 3.0 * tr_y_y_xyyzz[i] * fe_0 + ts_y_xyyyzz[i] * fe_0 + tr_y_y_xyyyzz[i] * pa_y[i];

        tr_y_yy_xyyzzz[i] = tr_y_0_xyyzzz[i] * fe_0 + 2.0 * tr_y_y_xyzzz[i] * fe_0 + ts_y_xyyzzz[i] * fe_0 + tr_y_y_xyyzzz[i] * pa_y[i];

        tr_y_yy_xyzzzz[i] = tr_y_0_xyzzzz[i] * fe_0 + tr_y_y_xzzzz[i] * fe_0 + ts_y_xyzzzz[i] * fe_0 + tr_y_y_xyzzzz[i] * pa_y[i];

        tr_y_yy_xzzzzz[i] = tr_y_0_xzzzzz[i] * fe_0 + ts_y_xzzzzz[i] * fe_0 + tr_y_y_xzzzzz[i] * pa_y[i];

        tr_y_yy_yyyyyy[i] = tr_y_0_yyyyyy[i] * fe_0 + 6.0 * tr_y_y_yyyyy[i] * fe_0 + ts_y_yyyyyy[i] * fe_0 + tr_y_y_yyyyyy[i] * pa_y[i];

        tr_y_yy_yyyyyz[i] = tr_y_0_yyyyyz[i] * fe_0 + 5.0 * tr_y_y_yyyyz[i] * fe_0 + ts_y_yyyyyz[i] * fe_0 + tr_y_y_yyyyyz[i] * pa_y[i];

        tr_y_yy_yyyyzz[i] = tr_y_0_yyyyzz[i] * fe_0 + 4.0 * tr_y_y_yyyzz[i] * fe_0 + ts_y_yyyyzz[i] * fe_0 + tr_y_y_yyyyzz[i] * pa_y[i];

        tr_y_yy_yyyzzz[i] = tr_y_0_yyyzzz[i] * fe_0 + 3.0 * tr_y_y_yyzzz[i] * fe_0 + ts_y_yyyzzz[i] * fe_0 + tr_y_y_yyyzzz[i] * pa_y[i];

        tr_y_yy_yyzzzz[i] = tr_y_0_yyzzzz[i] * fe_0 + 2.0 * tr_y_y_yzzzz[i] * fe_0 + ts_y_yyzzzz[i] * fe_0 + tr_y_y_yyzzzz[i] * pa_y[i];

        tr_y_yy_yzzzzz[i] = tr_y_0_yzzzzz[i] * fe_0 + tr_y_y_zzzzz[i] * fe_0 + ts_y_yzzzzz[i] * fe_0 + tr_y_y_yzzzzz[i] * pa_y[i];

        tr_y_yy_zzzzzz[i] = tr_y_0_zzzzzz[i] * fe_0 + ts_y_zzzzzz[i] * fe_0 + tr_y_y_zzzzzz[i] * pa_y[i];
    }

    // Set up 280-308 components of targeted buffer : DI

    auto tr_y_yz_xxxxxx = pbuffer.data(idx_dip_di + 280);

    auto tr_y_yz_xxxxxy = pbuffer.data(idx_dip_di + 281);

    auto tr_y_yz_xxxxxz = pbuffer.data(idx_dip_di + 282);

    auto tr_y_yz_xxxxyy = pbuffer.data(idx_dip_di + 283);

    auto tr_y_yz_xxxxyz = pbuffer.data(idx_dip_di + 284);

    auto tr_y_yz_xxxxzz = pbuffer.data(idx_dip_di + 285);

    auto tr_y_yz_xxxyyy = pbuffer.data(idx_dip_di + 286);

    auto tr_y_yz_xxxyyz = pbuffer.data(idx_dip_di + 287);

    auto tr_y_yz_xxxyzz = pbuffer.data(idx_dip_di + 288);

    auto tr_y_yz_xxxzzz = pbuffer.data(idx_dip_di + 289);

    auto tr_y_yz_xxyyyy = pbuffer.data(idx_dip_di + 290);

    auto tr_y_yz_xxyyyz = pbuffer.data(idx_dip_di + 291);

    auto tr_y_yz_xxyyzz = pbuffer.data(idx_dip_di + 292);

    auto tr_y_yz_xxyzzz = pbuffer.data(idx_dip_di + 293);

    auto tr_y_yz_xxzzzz = pbuffer.data(idx_dip_di + 294);

    auto tr_y_yz_xyyyyy = pbuffer.data(idx_dip_di + 295);

    auto tr_y_yz_xyyyyz = pbuffer.data(idx_dip_di + 296);

    auto tr_y_yz_xyyyzz = pbuffer.data(idx_dip_di + 297);

    auto tr_y_yz_xyyzzz = pbuffer.data(idx_dip_di + 298);

    auto tr_y_yz_xyzzzz = pbuffer.data(idx_dip_di + 299);

    auto tr_y_yz_xzzzzz = pbuffer.data(idx_dip_di + 300);

    auto tr_y_yz_yyyyyy = pbuffer.data(idx_dip_di + 301);

    auto tr_y_yz_yyyyyz = pbuffer.data(idx_dip_di + 302);

    auto tr_y_yz_yyyyzz = pbuffer.data(idx_dip_di + 303);

    auto tr_y_yz_yyyzzz = pbuffer.data(idx_dip_di + 304);

    auto tr_y_yz_yyzzzz = pbuffer.data(idx_dip_di + 305);

    auto tr_y_yz_yzzzzz = pbuffer.data(idx_dip_di + 306);

    auto tr_y_yz_zzzzzz = pbuffer.data(idx_dip_di + 307);

#pragma omp simd aligned(pa_y,               \
                             pa_z,           \
                             tr_y_y_xxxxxx,  \
                             tr_y_y_xxxxxy,  \
                             tr_y_y_xxxxy,   \
                             tr_y_y_xxxxyy,  \
                             tr_y_y_xxxxyz,  \
                             tr_y_y_xxxyy,   \
                             tr_y_y_xxxyyy,  \
                             tr_y_y_xxxyyz,  \
                             tr_y_y_xxxyz,   \
                             tr_y_y_xxxyzz,  \
                             tr_y_y_xxyyy,   \
                             tr_y_y_xxyyyy,  \
                             tr_y_y_xxyyyz,  \
                             tr_y_y_xxyyz,   \
                             tr_y_y_xxyyzz,  \
                             tr_y_y_xxyzz,   \
                             tr_y_y_xxyzzz,  \
                             tr_y_y_xyyyy,   \
                             tr_y_y_xyyyyy,  \
                             tr_y_y_xyyyyz,  \
                             tr_y_y_xyyyz,   \
                             tr_y_y_xyyyzz,  \
                             tr_y_y_xyyzz,   \
                             tr_y_y_xyyzzz,  \
                             tr_y_y_xyzzz,   \
                             tr_y_y_xyzzzz,  \
                             tr_y_y_yyyyy,   \
                             tr_y_y_yyyyyy,  \
                             tr_y_y_yyyyyz,  \
                             tr_y_y_yyyyz,   \
                             tr_y_y_yyyyzz,  \
                             tr_y_y_yyyzz,   \
                             tr_y_y_yyyzzz,  \
                             tr_y_y_yyzzz,   \
                             tr_y_y_yyzzzz,  \
                             tr_y_y_yzzzz,   \
                             tr_y_y_yzzzzz,  \
                             tr_y_yz_xxxxxx, \
                             tr_y_yz_xxxxxy, \
                             tr_y_yz_xxxxxz, \
                             tr_y_yz_xxxxyy, \
                             tr_y_yz_xxxxyz, \
                             tr_y_yz_xxxxzz, \
                             tr_y_yz_xxxyyy, \
                             tr_y_yz_xxxyyz, \
                             tr_y_yz_xxxyzz, \
                             tr_y_yz_xxxzzz, \
                             tr_y_yz_xxyyyy, \
                             tr_y_yz_xxyyyz, \
                             tr_y_yz_xxyyzz, \
                             tr_y_yz_xxyzzz, \
                             tr_y_yz_xxzzzz, \
                             tr_y_yz_xyyyyy, \
                             tr_y_yz_xyyyyz, \
                             tr_y_yz_xyyyzz, \
                             tr_y_yz_xyyzzz, \
                             tr_y_yz_xyzzzz, \
                             tr_y_yz_xzzzzz, \
                             tr_y_yz_yyyyyy, \
                             tr_y_yz_yyyyyz, \
                             tr_y_yz_yyyyzz, \
                             tr_y_yz_yyyzzz, \
                             tr_y_yz_yyzzzz, \
                             tr_y_yz_yzzzzz, \
                             tr_y_yz_zzzzzz, \
                             tr_y_z_xxxxxz,  \
                             tr_y_z_xxxxzz,  \
                             tr_y_z_xxxzzz,  \
                             tr_y_z_xxzzzz,  \
                             tr_y_z_xzzzzz,  \
                             tr_y_z_zzzzzz,  \
                             ts_z_xxxxxz,    \
                             ts_z_xxxxzz,    \
                             ts_z_xxxzzz,    \
                             ts_z_xxzzzz,    \
                             ts_z_xzzzzz,    \
                             ts_z_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_yz_xxxxxx[i] = tr_y_y_xxxxxx[i] * pa_z[i];

        tr_y_yz_xxxxxy[i] = tr_y_y_xxxxxy[i] * pa_z[i];

        tr_y_yz_xxxxxz[i] = ts_z_xxxxxz[i] * fe_0 + tr_y_z_xxxxxz[i] * pa_y[i];

        tr_y_yz_xxxxyy[i] = tr_y_y_xxxxyy[i] * pa_z[i];

        tr_y_yz_xxxxyz[i] = tr_y_y_xxxxy[i] * fe_0 + tr_y_y_xxxxyz[i] * pa_z[i];

        tr_y_yz_xxxxzz[i] = ts_z_xxxxzz[i] * fe_0 + tr_y_z_xxxxzz[i] * pa_y[i];

        tr_y_yz_xxxyyy[i] = tr_y_y_xxxyyy[i] * pa_z[i];

        tr_y_yz_xxxyyz[i] = tr_y_y_xxxyy[i] * fe_0 + tr_y_y_xxxyyz[i] * pa_z[i];

        tr_y_yz_xxxyzz[i] = 2.0 * tr_y_y_xxxyz[i] * fe_0 + tr_y_y_xxxyzz[i] * pa_z[i];

        tr_y_yz_xxxzzz[i] = ts_z_xxxzzz[i] * fe_0 + tr_y_z_xxxzzz[i] * pa_y[i];

        tr_y_yz_xxyyyy[i] = tr_y_y_xxyyyy[i] * pa_z[i];

        tr_y_yz_xxyyyz[i] = tr_y_y_xxyyy[i] * fe_0 + tr_y_y_xxyyyz[i] * pa_z[i];

        tr_y_yz_xxyyzz[i] = 2.0 * tr_y_y_xxyyz[i] * fe_0 + tr_y_y_xxyyzz[i] * pa_z[i];

        tr_y_yz_xxyzzz[i] = 3.0 * tr_y_y_xxyzz[i] * fe_0 + tr_y_y_xxyzzz[i] * pa_z[i];

        tr_y_yz_xxzzzz[i] = ts_z_xxzzzz[i] * fe_0 + tr_y_z_xxzzzz[i] * pa_y[i];

        tr_y_yz_xyyyyy[i] = tr_y_y_xyyyyy[i] * pa_z[i];

        tr_y_yz_xyyyyz[i] = tr_y_y_xyyyy[i] * fe_0 + tr_y_y_xyyyyz[i] * pa_z[i];

        tr_y_yz_xyyyzz[i] = 2.0 * tr_y_y_xyyyz[i] * fe_0 + tr_y_y_xyyyzz[i] * pa_z[i];

        tr_y_yz_xyyzzz[i] = 3.0 * tr_y_y_xyyzz[i] * fe_0 + tr_y_y_xyyzzz[i] * pa_z[i];

        tr_y_yz_xyzzzz[i] = 4.0 * tr_y_y_xyzzz[i] * fe_0 + tr_y_y_xyzzzz[i] * pa_z[i];

        tr_y_yz_xzzzzz[i] = ts_z_xzzzzz[i] * fe_0 + tr_y_z_xzzzzz[i] * pa_y[i];

        tr_y_yz_yyyyyy[i] = tr_y_y_yyyyyy[i] * pa_z[i];

        tr_y_yz_yyyyyz[i] = tr_y_y_yyyyy[i] * fe_0 + tr_y_y_yyyyyz[i] * pa_z[i];

        tr_y_yz_yyyyzz[i] = 2.0 * tr_y_y_yyyyz[i] * fe_0 + tr_y_y_yyyyzz[i] * pa_z[i];

        tr_y_yz_yyyzzz[i] = 3.0 * tr_y_y_yyyzz[i] * fe_0 + tr_y_y_yyyzzz[i] * pa_z[i];

        tr_y_yz_yyzzzz[i] = 4.0 * tr_y_y_yyzzz[i] * fe_0 + tr_y_y_yyzzzz[i] * pa_z[i];

        tr_y_yz_yzzzzz[i] = 5.0 * tr_y_y_yzzzz[i] * fe_0 + tr_y_y_yzzzzz[i] * pa_z[i];

        tr_y_yz_zzzzzz[i] = ts_z_zzzzzz[i] * fe_0 + tr_y_z_zzzzzz[i] * pa_y[i];
    }

    // Set up 308-336 components of targeted buffer : DI

    auto tr_y_zz_xxxxxx = pbuffer.data(idx_dip_di + 308);

    auto tr_y_zz_xxxxxy = pbuffer.data(idx_dip_di + 309);

    auto tr_y_zz_xxxxxz = pbuffer.data(idx_dip_di + 310);

    auto tr_y_zz_xxxxyy = pbuffer.data(idx_dip_di + 311);

    auto tr_y_zz_xxxxyz = pbuffer.data(idx_dip_di + 312);

    auto tr_y_zz_xxxxzz = pbuffer.data(idx_dip_di + 313);

    auto tr_y_zz_xxxyyy = pbuffer.data(idx_dip_di + 314);

    auto tr_y_zz_xxxyyz = pbuffer.data(idx_dip_di + 315);

    auto tr_y_zz_xxxyzz = pbuffer.data(idx_dip_di + 316);

    auto tr_y_zz_xxxzzz = pbuffer.data(idx_dip_di + 317);

    auto tr_y_zz_xxyyyy = pbuffer.data(idx_dip_di + 318);

    auto tr_y_zz_xxyyyz = pbuffer.data(idx_dip_di + 319);

    auto tr_y_zz_xxyyzz = pbuffer.data(idx_dip_di + 320);

    auto tr_y_zz_xxyzzz = pbuffer.data(idx_dip_di + 321);

    auto tr_y_zz_xxzzzz = pbuffer.data(idx_dip_di + 322);

    auto tr_y_zz_xyyyyy = pbuffer.data(idx_dip_di + 323);

    auto tr_y_zz_xyyyyz = pbuffer.data(idx_dip_di + 324);

    auto tr_y_zz_xyyyzz = pbuffer.data(idx_dip_di + 325);

    auto tr_y_zz_xyyzzz = pbuffer.data(idx_dip_di + 326);

    auto tr_y_zz_xyzzzz = pbuffer.data(idx_dip_di + 327);

    auto tr_y_zz_xzzzzz = pbuffer.data(idx_dip_di + 328);

    auto tr_y_zz_yyyyyy = pbuffer.data(idx_dip_di + 329);

    auto tr_y_zz_yyyyyz = pbuffer.data(idx_dip_di + 330);

    auto tr_y_zz_yyyyzz = pbuffer.data(idx_dip_di + 331);

    auto tr_y_zz_yyyzzz = pbuffer.data(idx_dip_di + 332);

    auto tr_y_zz_yyzzzz = pbuffer.data(idx_dip_di + 333);

    auto tr_y_zz_yzzzzz = pbuffer.data(idx_dip_di + 334);

    auto tr_y_zz_zzzzzz = pbuffer.data(idx_dip_di + 335);

#pragma omp simd aligned(pa_z,               \
                             tr_y_0_xxxxxx,  \
                             tr_y_0_xxxxxy,  \
                             tr_y_0_xxxxxz,  \
                             tr_y_0_xxxxyy,  \
                             tr_y_0_xxxxyz,  \
                             tr_y_0_xxxxzz,  \
                             tr_y_0_xxxyyy,  \
                             tr_y_0_xxxyyz,  \
                             tr_y_0_xxxyzz,  \
                             tr_y_0_xxxzzz,  \
                             tr_y_0_xxyyyy,  \
                             tr_y_0_xxyyyz,  \
                             tr_y_0_xxyyzz,  \
                             tr_y_0_xxyzzz,  \
                             tr_y_0_xxzzzz,  \
                             tr_y_0_xyyyyy,  \
                             tr_y_0_xyyyyz,  \
                             tr_y_0_xyyyzz,  \
                             tr_y_0_xyyzzz,  \
                             tr_y_0_xyzzzz,  \
                             tr_y_0_xzzzzz,  \
                             tr_y_0_yyyyyy,  \
                             tr_y_0_yyyyyz,  \
                             tr_y_0_yyyyzz,  \
                             tr_y_0_yyyzzz,  \
                             tr_y_0_yyzzzz,  \
                             tr_y_0_yzzzzz,  \
                             tr_y_0_zzzzzz,  \
                             tr_y_z_xxxxx,   \
                             tr_y_z_xxxxxx,  \
                             tr_y_z_xxxxxy,  \
                             tr_y_z_xxxxxz,  \
                             tr_y_z_xxxxy,   \
                             tr_y_z_xxxxyy,  \
                             tr_y_z_xxxxyz,  \
                             tr_y_z_xxxxz,   \
                             tr_y_z_xxxxzz,  \
                             tr_y_z_xxxyy,   \
                             tr_y_z_xxxyyy,  \
                             tr_y_z_xxxyyz,  \
                             tr_y_z_xxxyz,   \
                             tr_y_z_xxxyzz,  \
                             tr_y_z_xxxzz,   \
                             tr_y_z_xxxzzz,  \
                             tr_y_z_xxyyy,   \
                             tr_y_z_xxyyyy,  \
                             tr_y_z_xxyyyz,  \
                             tr_y_z_xxyyz,   \
                             tr_y_z_xxyyzz,  \
                             tr_y_z_xxyzz,   \
                             tr_y_z_xxyzzz,  \
                             tr_y_z_xxzzz,   \
                             tr_y_z_xxzzzz,  \
                             tr_y_z_xyyyy,   \
                             tr_y_z_xyyyyy,  \
                             tr_y_z_xyyyyz,  \
                             tr_y_z_xyyyz,   \
                             tr_y_z_xyyyzz,  \
                             tr_y_z_xyyzz,   \
                             tr_y_z_xyyzzz,  \
                             tr_y_z_xyzzz,   \
                             tr_y_z_xyzzzz,  \
                             tr_y_z_xzzzz,   \
                             tr_y_z_xzzzzz,  \
                             tr_y_z_yyyyy,   \
                             tr_y_z_yyyyyy,  \
                             tr_y_z_yyyyyz,  \
                             tr_y_z_yyyyz,   \
                             tr_y_z_yyyyzz,  \
                             tr_y_z_yyyzz,   \
                             tr_y_z_yyyzzz,  \
                             tr_y_z_yyzzz,   \
                             tr_y_z_yyzzzz,  \
                             tr_y_z_yzzzz,   \
                             tr_y_z_yzzzzz,  \
                             tr_y_z_zzzzz,   \
                             tr_y_z_zzzzzz,  \
                             tr_y_zz_xxxxxx, \
                             tr_y_zz_xxxxxy, \
                             tr_y_zz_xxxxxz, \
                             tr_y_zz_xxxxyy, \
                             tr_y_zz_xxxxyz, \
                             tr_y_zz_xxxxzz, \
                             tr_y_zz_xxxyyy, \
                             tr_y_zz_xxxyyz, \
                             tr_y_zz_xxxyzz, \
                             tr_y_zz_xxxzzz, \
                             tr_y_zz_xxyyyy, \
                             tr_y_zz_xxyyyz, \
                             tr_y_zz_xxyyzz, \
                             tr_y_zz_xxyzzz, \
                             tr_y_zz_xxzzzz, \
                             tr_y_zz_xyyyyy, \
                             tr_y_zz_xyyyyz, \
                             tr_y_zz_xyyyzz, \
                             tr_y_zz_xyyzzz, \
                             tr_y_zz_xyzzzz, \
                             tr_y_zz_xzzzzz, \
                             tr_y_zz_yyyyyy, \
                             tr_y_zz_yyyyyz, \
                             tr_y_zz_yyyyzz, \
                             tr_y_zz_yyyzzz, \
                             tr_y_zz_yyzzzz, \
                             tr_y_zz_yzzzzz, \
                             tr_y_zz_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_y_zz_xxxxxx[i] = tr_y_0_xxxxxx[i] * fe_0 + tr_y_z_xxxxxx[i] * pa_z[i];

        tr_y_zz_xxxxxy[i] = tr_y_0_xxxxxy[i] * fe_0 + tr_y_z_xxxxxy[i] * pa_z[i];

        tr_y_zz_xxxxxz[i] = tr_y_0_xxxxxz[i] * fe_0 + tr_y_z_xxxxx[i] * fe_0 + tr_y_z_xxxxxz[i] * pa_z[i];

        tr_y_zz_xxxxyy[i] = tr_y_0_xxxxyy[i] * fe_0 + tr_y_z_xxxxyy[i] * pa_z[i];

        tr_y_zz_xxxxyz[i] = tr_y_0_xxxxyz[i] * fe_0 + tr_y_z_xxxxy[i] * fe_0 + tr_y_z_xxxxyz[i] * pa_z[i];

        tr_y_zz_xxxxzz[i] = tr_y_0_xxxxzz[i] * fe_0 + 2.0 * tr_y_z_xxxxz[i] * fe_0 + tr_y_z_xxxxzz[i] * pa_z[i];

        tr_y_zz_xxxyyy[i] = tr_y_0_xxxyyy[i] * fe_0 + tr_y_z_xxxyyy[i] * pa_z[i];

        tr_y_zz_xxxyyz[i] = tr_y_0_xxxyyz[i] * fe_0 + tr_y_z_xxxyy[i] * fe_0 + tr_y_z_xxxyyz[i] * pa_z[i];

        tr_y_zz_xxxyzz[i] = tr_y_0_xxxyzz[i] * fe_0 + 2.0 * tr_y_z_xxxyz[i] * fe_0 + tr_y_z_xxxyzz[i] * pa_z[i];

        tr_y_zz_xxxzzz[i] = tr_y_0_xxxzzz[i] * fe_0 + 3.0 * tr_y_z_xxxzz[i] * fe_0 + tr_y_z_xxxzzz[i] * pa_z[i];

        tr_y_zz_xxyyyy[i] = tr_y_0_xxyyyy[i] * fe_0 + tr_y_z_xxyyyy[i] * pa_z[i];

        tr_y_zz_xxyyyz[i] = tr_y_0_xxyyyz[i] * fe_0 + tr_y_z_xxyyy[i] * fe_0 + tr_y_z_xxyyyz[i] * pa_z[i];

        tr_y_zz_xxyyzz[i] = tr_y_0_xxyyzz[i] * fe_0 + 2.0 * tr_y_z_xxyyz[i] * fe_0 + tr_y_z_xxyyzz[i] * pa_z[i];

        tr_y_zz_xxyzzz[i] = tr_y_0_xxyzzz[i] * fe_0 + 3.0 * tr_y_z_xxyzz[i] * fe_0 + tr_y_z_xxyzzz[i] * pa_z[i];

        tr_y_zz_xxzzzz[i] = tr_y_0_xxzzzz[i] * fe_0 + 4.0 * tr_y_z_xxzzz[i] * fe_0 + tr_y_z_xxzzzz[i] * pa_z[i];

        tr_y_zz_xyyyyy[i] = tr_y_0_xyyyyy[i] * fe_0 + tr_y_z_xyyyyy[i] * pa_z[i];

        tr_y_zz_xyyyyz[i] = tr_y_0_xyyyyz[i] * fe_0 + tr_y_z_xyyyy[i] * fe_0 + tr_y_z_xyyyyz[i] * pa_z[i];

        tr_y_zz_xyyyzz[i] = tr_y_0_xyyyzz[i] * fe_0 + 2.0 * tr_y_z_xyyyz[i] * fe_0 + tr_y_z_xyyyzz[i] * pa_z[i];

        tr_y_zz_xyyzzz[i] = tr_y_0_xyyzzz[i] * fe_0 + 3.0 * tr_y_z_xyyzz[i] * fe_0 + tr_y_z_xyyzzz[i] * pa_z[i];

        tr_y_zz_xyzzzz[i] = tr_y_0_xyzzzz[i] * fe_0 + 4.0 * tr_y_z_xyzzz[i] * fe_0 + tr_y_z_xyzzzz[i] * pa_z[i];

        tr_y_zz_xzzzzz[i] = tr_y_0_xzzzzz[i] * fe_0 + 5.0 * tr_y_z_xzzzz[i] * fe_0 + tr_y_z_xzzzzz[i] * pa_z[i];

        tr_y_zz_yyyyyy[i] = tr_y_0_yyyyyy[i] * fe_0 + tr_y_z_yyyyyy[i] * pa_z[i];

        tr_y_zz_yyyyyz[i] = tr_y_0_yyyyyz[i] * fe_0 + tr_y_z_yyyyy[i] * fe_0 + tr_y_z_yyyyyz[i] * pa_z[i];

        tr_y_zz_yyyyzz[i] = tr_y_0_yyyyzz[i] * fe_0 + 2.0 * tr_y_z_yyyyz[i] * fe_0 + tr_y_z_yyyyzz[i] * pa_z[i];

        tr_y_zz_yyyzzz[i] = tr_y_0_yyyzzz[i] * fe_0 + 3.0 * tr_y_z_yyyzz[i] * fe_0 + tr_y_z_yyyzzz[i] * pa_z[i];

        tr_y_zz_yyzzzz[i] = tr_y_0_yyzzzz[i] * fe_0 + 4.0 * tr_y_z_yyzzz[i] * fe_0 + tr_y_z_yyzzzz[i] * pa_z[i];

        tr_y_zz_yzzzzz[i] = tr_y_0_yzzzzz[i] * fe_0 + 5.0 * tr_y_z_yzzzz[i] * fe_0 + tr_y_z_yzzzzz[i] * pa_z[i];

        tr_y_zz_zzzzzz[i] = tr_y_0_zzzzzz[i] * fe_0 + 6.0 * tr_y_z_zzzzz[i] * fe_0 + tr_y_z_zzzzzz[i] * pa_z[i];
    }

    // Set up 336-364 components of targeted buffer : DI

    auto tr_z_xx_xxxxxx = pbuffer.data(idx_dip_di + 336);

    auto tr_z_xx_xxxxxy = pbuffer.data(idx_dip_di + 337);

    auto tr_z_xx_xxxxxz = pbuffer.data(idx_dip_di + 338);

    auto tr_z_xx_xxxxyy = pbuffer.data(idx_dip_di + 339);

    auto tr_z_xx_xxxxyz = pbuffer.data(idx_dip_di + 340);

    auto tr_z_xx_xxxxzz = pbuffer.data(idx_dip_di + 341);

    auto tr_z_xx_xxxyyy = pbuffer.data(idx_dip_di + 342);

    auto tr_z_xx_xxxyyz = pbuffer.data(idx_dip_di + 343);

    auto tr_z_xx_xxxyzz = pbuffer.data(idx_dip_di + 344);

    auto tr_z_xx_xxxzzz = pbuffer.data(idx_dip_di + 345);

    auto tr_z_xx_xxyyyy = pbuffer.data(idx_dip_di + 346);

    auto tr_z_xx_xxyyyz = pbuffer.data(idx_dip_di + 347);

    auto tr_z_xx_xxyyzz = pbuffer.data(idx_dip_di + 348);

    auto tr_z_xx_xxyzzz = pbuffer.data(idx_dip_di + 349);

    auto tr_z_xx_xxzzzz = pbuffer.data(idx_dip_di + 350);

    auto tr_z_xx_xyyyyy = pbuffer.data(idx_dip_di + 351);

    auto tr_z_xx_xyyyyz = pbuffer.data(idx_dip_di + 352);

    auto tr_z_xx_xyyyzz = pbuffer.data(idx_dip_di + 353);

    auto tr_z_xx_xyyzzz = pbuffer.data(idx_dip_di + 354);

    auto tr_z_xx_xyzzzz = pbuffer.data(idx_dip_di + 355);

    auto tr_z_xx_xzzzzz = pbuffer.data(idx_dip_di + 356);

    auto tr_z_xx_yyyyyy = pbuffer.data(idx_dip_di + 357);

    auto tr_z_xx_yyyyyz = pbuffer.data(idx_dip_di + 358);

    auto tr_z_xx_yyyyzz = pbuffer.data(idx_dip_di + 359);

    auto tr_z_xx_yyyzzz = pbuffer.data(idx_dip_di + 360);

    auto tr_z_xx_yyzzzz = pbuffer.data(idx_dip_di + 361);

    auto tr_z_xx_yzzzzz = pbuffer.data(idx_dip_di + 362);

    auto tr_z_xx_zzzzzz = pbuffer.data(idx_dip_di + 363);

#pragma omp simd aligned(pa_x,               \
                             tr_z_0_xxxxxx,  \
                             tr_z_0_xxxxxy,  \
                             tr_z_0_xxxxxz,  \
                             tr_z_0_xxxxyy,  \
                             tr_z_0_xxxxyz,  \
                             tr_z_0_xxxxzz,  \
                             tr_z_0_xxxyyy,  \
                             tr_z_0_xxxyyz,  \
                             tr_z_0_xxxyzz,  \
                             tr_z_0_xxxzzz,  \
                             tr_z_0_xxyyyy,  \
                             tr_z_0_xxyyyz,  \
                             tr_z_0_xxyyzz,  \
                             tr_z_0_xxyzzz,  \
                             tr_z_0_xxzzzz,  \
                             tr_z_0_xyyyyy,  \
                             tr_z_0_xyyyyz,  \
                             tr_z_0_xyyyzz,  \
                             tr_z_0_xyyzzz,  \
                             tr_z_0_xyzzzz,  \
                             tr_z_0_xzzzzz,  \
                             tr_z_0_yyyyyy,  \
                             tr_z_0_yyyyyz,  \
                             tr_z_0_yyyyzz,  \
                             tr_z_0_yyyzzz,  \
                             tr_z_0_yyzzzz,  \
                             tr_z_0_yzzzzz,  \
                             tr_z_0_zzzzzz,  \
                             tr_z_x_xxxxx,   \
                             tr_z_x_xxxxxx,  \
                             tr_z_x_xxxxxy,  \
                             tr_z_x_xxxxxz,  \
                             tr_z_x_xxxxy,   \
                             tr_z_x_xxxxyy,  \
                             tr_z_x_xxxxyz,  \
                             tr_z_x_xxxxz,   \
                             tr_z_x_xxxxzz,  \
                             tr_z_x_xxxyy,   \
                             tr_z_x_xxxyyy,  \
                             tr_z_x_xxxyyz,  \
                             tr_z_x_xxxyz,   \
                             tr_z_x_xxxyzz,  \
                             tr_z_x_xxxzz,   \
                             tr_z_x_xxxzzz,  \
                             tr_z_x_xxyyy,   \
                             tr_z_x_xxyyyy,  \
                             tr_z_x_xxyyyz,  \
                             tr_z_x_xxyyz,   \
                             tr_z_x_xxyyzz,  \
                             tr_z_x_xxyzz,   \
                             tr_z_x_xxyzzz,  \
                             tr_z_x_xxzzz,   \
                             tr_z_x_xxzzzz,  \
                             tr_z_x_xyyyy,   \
                             tr_z_x_xyyyyy,  \
                             tr_z_x_xyyyyz,  \
                             tr_z_x_xyyyz,   \
                             tr_z_x_xyyyzz,  \
                             tr_z_x_xyyzz,   \
                             tr_z_x_xyyzzz,  \
                             tr_z_x_xyzzz,   \
                             tr_z_x_xyzzzz,  \
                             tr_z_x_xzzzz,   \
                             tr_z_x_xzzzzz,  \
                             tr_z_x_yyyyy,   \
                             tr_z_x_yyyyyy,  \
                             tr_z_x_yyyyyz,  \
                             tr_z_x_yyyyz,   \
                             tr_z_x_yyyyzz,  \
                             tr_z_x_yyyzz,   \
                             tr_z_x_yyyzzz,  \
                             tr_z_x_yyzzz,   \
                             tr_z_x_yyzzzz,  \
                             tr_z_x_yzzzz,   \
                             tr_z_x_yzzzzz,  \
                             tr_z_x_zzzzz,   \
                             tr_z_x_zzzzzz,  \
                             tr_z_xx_xxxxxx, \
                             tr_z_xx_xxxxxy, \
                             tr_z_xx_xxxxxz, \
                             tr_z_xx_xxxxyy, \
                             tr_z_xx_xxxxyz, \
                             tr_z_xx_xxxxzz, \
                             tr_z_xx_xxxyyy, \
                             tr_z_xx_xxxyyz, \
                             tr_z_xx_xxxyzz, \
                             tr_z_xx_xxxzzz, \
                             tr_z_xx_xxyyyy, \
                             tr_z_xx_xxyyyz, \
                             tr_z_xx_xxyyzz, \
                             tr_z_xx_xxyzzz, \
                             tr_z_xx_xxzzzz, \
                             tr_z_xx_xyyyyy, \
                             tr_z_xx_xyyyyz, \
                             tr_z_xx_xyyyzz, \
                             tr_z_xx_xyyzzz, \
                             tr_z_xx_xyzzzz, \
                             tr_z_xx_xzzzzz, \
                             tr_z_xx_yyyyyy, \
                             tr_z_xx_yyyyyz, \
                             tr_z_xx_yyyyzz, \
                             tr_z_xx_yyyzzz, \
                             tr_z_xx_yyzzzz, \
                             tr_z_xx_yzzzzz, \
                             tr_z_xx_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xx_xxxxxx[i] = tr_z_0_xxxxxx[i] * fe_0 + 6.0 * tr_z_x_xxxxx[i] * fe_0 + tr_z_x_xxxxxx[i] * pa_x[i];

        tr_z_xx_xxxxxy[i] = tr_z_0_xxxxxy[i] * fe_0 + 5.0 * tr_z_x_xxxxy[i] * fe_0 + tr_z_x_xxxxxy[i] * pa_x[i];

        tr_z_xx_xxxxxz[i] = tr_z_0_xxxxxz[i] * fe_0 + 5.0 * tr_z_x_xxxxz[i] * fe_0 + tr_z_x_xxxxxz[i] * pa_x[i];

        tr_z_xx_xxxxyy[i] = tr_z_0_xxxxyy[i] * fe_0 + 4.0 * tr_z_x_xxxyy[i] * fe_0 + tr_z_x_xxxxyy[i] * pa_x[i];

        tr_z_xx_xxxxyz[i] = tr_z_0_xxxxyz[i] * fe_0 + 4.0 * tr_z_x_xxxyz[i] * fe_0 + tr_z_x_xxxxyz[i] * pa_x[i];

        tr_z_xx_xxxxzz[i] = tr_z_0_xxxxzz[i] * fe_0 + 4.0 * tr_z_x_xxxzz[i] * fe_0 + tr_z_x_xxxxzz[i] * pa_x[i];

        tr_z_xx_xxxyyy[i] = tr_z_0_xxxyyy[i] * fe_0 + 3.0 * tr_z_x_xxyyy[i] * fe_0 + tr_z_x_xxxyyy[i] * pa_x[i];

        tr_z_xx_xxxyyz[i] = tr_z_0_xxxyyz[i] * fe_0 + 3.0 * tr_z_x_xxyyz[i] * fe_0 + tr_z_x_xxxyyz[i] * pa_x[i];

        tr_z_xx_xxxyzz[i] = tr_z_0_xxxyzz[i] * fe_0 + 3.0 * tr_z_x_xxyzz[i] * fe_0 + tr_z_x_xxxyzz[i] * pa_x[i];

        tr_z_xx_xxxzzz[i] = tr_z_0_xxxzzz[i] * fe_0 + 3.0 * tr_z_x_xxzzz[i] * fe_0 + tr_z_x_xxxzzz[i] * pa_x[i];

        tr_z_xx_xxyyyy[i] = tr_z_0_xxyyyy[i] * fe_0 + 2.0 * tr_z_x_xyyyy[i] * fe_0 + tr_z_x_xxyyyy[i] * pa_x[i];

        tr_z_xx_xxyyyz[i] = tr_z_0_xxyyyz[i] * fe_0 + 2.0 * tr_z_x_xyyyz[i] * fe_0 + tr_z_x_xxyyyz[i] * pa_x[i];

        tr_z_xx_xxyyzz[i] = tr_z_0_xxyyzz[i] * fe_0 + 2.0 * tr_z_x_xyyzz[i] * fe_0 + tr_z_x_xxyyzz[i] * pa_x[i];

        tr_z_xx_xxyzzz[i] = tr_z_0_xxyzzz[i] * fe_0 + 2.0 * tr_z_x_xyzzz[i] * fe_0 + tr_z_x_xxyzzz[i] * pa_x[i];

        tr_z_xx_xxzzzz[i] = tr_z_0_xxzzzz[i] * fe_0 + 2.0 * tr_z_x_xzzzz[i] * fe_0 + tr_z_x_xxzzzz[i] * pa_x[i];

        tr_z_xx_xyyyyy[i] = tr_z_0_xyyyyy[i] * fe_0 + tr_z_x_yyyyy[i] * fe_0 + tr_z_x_xyyyyy[i] * pa_x[i];

        tr_z_xx_xyyyyz[i] = tr_z_0_xyyyyz[i] * fe_0 + tr_z_x_yyyyz[i] * fe_0 + tr_z_x_xyyyyz[i] * pa_x[i];

        tr_z_xx_xyyyzz[i] = tr_z_0_xyyyzz[i] * fe_0 + tr_z_x_yyyzz[i] * fe_0 + tr_z_x_xyyyzz[i] * pa_x[i];

        tr_z_xx_xyyzzz[i] = tr_z_0_xyyzzz[i] * fe_0 + tr_z_x_yyzzz[i] * fe_0 + tr_z_x_xyyzzz[i] * pa_x[i];

        tr_z_xx_xyzzzz[i] = tr_z_0_xyzzzz[i] * fe_0 + tr_z_x_yzzzz[i] * fe_0 + tr_z_x_xyzzzz[i] * pa_x[i];

        tr_z_xx_xzzzzz[i] = tr_z_0_xzzzzz[i] * fe_0 + tr_z_x_zzzzz[i] * fe_0 + tr_z_x_xzzzzz[i] * pa_x[i];

        tr_z_xx_yyyyyy[i] = tr_z_0_yyyyyy[i] * fe_0 + tr_z_x_yyyyyy[i] * pa_x[i];

        tr_z_xx_yyyyyz[i] = tr_z_0_yyyyyz[i] * fe_0 + tr_z_x_yyyyyz[i] * pa_x[i];

        tr_z_xx_yyyyzz[i] = tr_z_0_yyyyzz[i] * fe_0 + tr_z_x_yyyyzz[i] * pa_x[i];

        tr_z_xx_yyyzzz[i] = tr_z_0_yyyzzz[i] * fe_0 + tr_z_x_yyyzzz[i] * pa_x[i];

        tr_z_xx_yyzzzz[i] = tr_z_0_yyzzzz[i] * fe_0 + tr_z_x_yyzzzz[i] * pa_x[i];

        tr_z_xx_yzzzzz[i] = tr_z_0_yzzzzz[i] * fe_0 + tr_z_x_yzzzzz[i] * pa_x[i];

        tr_z_xx_zzzzzz[i] = tr_z_0_zzzzzz[i] * fe_0 + tr_z_x_zzzzzz[i] * pa_x[i];
    }

    // Set up 364-392 components of targeted buffer : DI

    auto tr_z_xy_xxxxxx = pbuffer.data(idx_dip_di + 364);

    auto tr_z_xy_xxxxxy = pbuffer.data(idx_dip_di + 365);

    auto tr_z_xy_xxxxxz = pbuffer.data(idx_dip_di + 366);

    auto tr_z_xy_xxxxyy = pbuffer.data(idx_dip_di + 367);

    auto tr_z_xy_xxxxyz = pbuffer.data(idx_dip_di + 368);

    auto tr_z_xy_xxxxzz = pbuffer.data(idx_dip_di + 369);

    auto tr_z_xy_xxxyyy = pbuffer.data(idx_dip_di + 370);

    auto tr_z_xy_xxxyyz = pbuffer.data(idx_dip_di + 371);

    auto tr_z_xy_xxxyzz = pbuffer.data(idx_dip_di + 372);

    auto tr_z_xy_xxxzzz = pbuffer.data(idx_dip_di + 373);

    auto tr_z_xy_xxyyyy = pbuffer.data(idx_dip_di + 374);

    auto tr_z_xy_xxyyyz = pbuffer.data(idx_dip_di + 375);

    auto tr_z_xy_xxyyzz = pbuffer.data(idx_dip_di + 376);

    auto tr_z_xy_xxyzzz = pbuffer.data(idx_dip_di + 377);

    auto tr_z_xy_xxzzzz = pbuffer.data(idx_dip_di + 378);

    auto tr_z_xy_xyyyyy = pbuffer.data(idx_dip_di + 379);

    auto tr_z_xy_xyyyyz = pbuffer.data(idx_dip_di + 380);

    auto tr_z_xy_xyyyzz = pbuffer.data(idx_dip_di + 381);

    auto tr_z_xy_xyyzzz = pbuffer.data(idx_dip_di + 382);

    auto tr_z_xy_xyzzzz = pbuffer.data(idx_dip_di + 383);

    auto tr_z_xy_xzzzzz = pbuffer.data(idx_dip_di + 384);

    auto tr_z_xy_yyyyyy = pbuffer.data(idx_dip_di + 385);

    auto tr_z_xy_yyyyyz = pbuffer.data(idx_dip_di + 386);

    auto tr_z_xy_yyyyzz = pbuffer.data(idx_dip_di + 387);

    auto tr_z_xy_yyyzzz = pbuffer.data(idx_dip_di + 388);

    auto tr_z_xy_yyzzzz = pbuffer.data(idx_dip_di + 389);

    auto tr_z_xy_yzzzzz = pbuffer.data(idx_dip_di + 390);

    auto tr_z_xy_zzzzzz = pbuffer.data(idx_dip_di + 391);

#pragma omp simd aligned(pa_x,               \
                             pa_y,           \
                             tr_z_x_xxxxxx,  \
                             tr_z_x_xxxxxz,  \
                             tr_z_x_xxxxzz,  \
                             tr_z_x_xxxzzz,  \
                             tr_z_x_xxzzzz,  \
                             tr_z_x_xzzzzz,  \
                             tr_z_xy_xxxxxx, \
                             tr_z_xy_xxxxxy, \
                             tr_z_xy_xxxxxz, \
                             tr_z_xy_xxxxyy, \
                             tr_z_xy_xxxxyz, \
                             tr_z_xy_xxxxzz, \
                             tr_z_xy_xxxyyy, \
                             tr_z_xy_xxxyyz, \
                             tr_z_xy_xxxyzz, \
                             tr_z_xy_xxxzzz, \
                             tr_z_xy_xxyyyy, \
                             tr_z_xy_xxyyyz, \
                             tr_z_xy_xxyyzz, \
                             tr_z_xy_xxyzzz, \
                             tr_z_xy_xxzzzz, \
                             tr_z_xy_xyyyyy, \
                             tr_z_xy_xyyyyz, \
                             tr_z_xy_xyyyzz, \
                             tr_z_xy_xyyzzz, \
                             tr_z_xy_xyzzzz, \
                             tr_z_xy_xzzzzz, \
                             tr_z_xy_yyyyyy, \
                             tr_z_xy_yyyyyz, \
                             tr_z_xy_yyyyzz, \
                             tr_z_xy_yyyzzz, \
                             tr_z_xy_yyzzzz, \
                             tr_z_xy_yzzzzz, \
                             tr_z_xy_zzzzzz, \
                             tr_z_y_xxxxxy,  \
                             tr_z_y_xxxxy,   \
                             tr_z_y_xxxxyy,  \
                             tr_z_y_xxxxyz,  \
                             tr_z_y_xxxyy,   \
                             tr_z_y_xxxyyy,  \
                             tr_z_y_xxxyyz,  \
                             tr_z_y_xxxyz,   \
                             tr_z_y_xxxyzz,  \
                             tr_z_y_xxyyy,   \
                             tr_z_y_xxyyyy,  \
                             tr_z_y_xxyyyz,  \
                             tr_z_y_xxyyz,   \
                             tr_z_y_xxyyzz,  \
                             tr_z_y_xxyzz,   \
                             tr_z_y_xxyzzz,  \
                             tr_z_y_xyyyy,   \
                             tr_z_y_xyyyyy,  \
                             tr_z_y_xyyyyz,  \
                             tr_z_y_xyyyz,   \
                             tr_z_y_xyyyzz,  \
                             tr_z_y_xyyzz,   \
                             tr_z_y_xyyzzz,  \
                             tr_z_y_xyzzz,   \
                             tr_z_y_xyzzzz,  \
                             tr_z_y_yyyyy,   \
                             tr_z_y_yyyyyy,  \
                             tr_z_y_yyyyyz,  \
                             tr_z_y_yyyyz,   \
                             tr_z_y_yyyyzz,  \
                             tr_z_y_yyyzz,   \
                             tr_z_y_yyyzzz,  \
                             tr_z_y_yyzzz,   \
                             tr_z_y_yyzzzz,  \
                             tr_z_y_yzzzz,   \
                             tr_z_y_yzzzzz,  \
                             tr_z_y_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xy_xxxxxx[i] = tr_z_x_xxxxxx[i] * pa_y[i];

        tr_z_xy_xxxxxy[i] = 5.0 * tr_z_y_xxxxy[i] * fe_0 + tr_z_y_xxxxxy[i] * pa_x[i];

        tr_z_xy_xxxxxz[i] = tr_z_x_xxxxxz[i] * pa_y[i];

        tr_z_xy_xxxxyy[i] = 4.0 * tr_z_y_xxxyy[i] * fe_0 + tr_z_y_xxxxyy[i] * pa_x[i];

        tr_z_xy_xxxxyz[i] = 4.0 * tr_z_y_xxxyz[i] * fe_0 + tr_z_y_xxxxyz[i] * pa_x[i];

        tr_z_xy_xxxxzz[i] = tr_z_x_xxxxzz[i] * pa_y[i];

        tr_z_xy_xxxyyy[i] = 3.0 * tr_z_y_xxyyy[i] * fe_0 + tr_z_y_xxxyyy[i] * pa_x[i];

        tr_z_xy_xxxyyz[i] = 3.0 * tr_z_y_xxyyz[i] * fe_0 + tr_z_y_xxxyyz[i] * pa_x[i];

        tr_z_xy_xxxyzz[i] = 3.0 * tr_z_y_xxyzz[i] * fe_0 + tr_z_y_xxxyzz[i] * pa_x[i];

        tr_z_xy_xxxzzz[i] = tr_z_x_xxxzzz[i] * pa_y[i];

        tr_z_xy_xxyyyy[i] = 2.0 * tr_z_y_xyyyy[i] * fe_0 + tr_z_y_xxyyyy[i] * pa_x[i];

        tr_z_xy_xxyyyz[i] = 2.0 * tr_z_y_xyyyz[i] * fe_0 + tr_z_y_xxyyyz[i] * pa_x[i];

        tr_z_xy_xxyyzz[i] = 2.0 * tr_z_y_xyyzz[i] * fe_0 + tr_z_y_xxyyzz[i] * pa_x[i];

        tr_z_xy_xxyzzz[i] = 2.0 * tr_z_y_xyzzz[i] * fe_0 + tr_z_y_xxyzzz[i] * pa_x[i];

        tr_z_xy_xxzzzz[i] = tr_z_x_xxzzzz[i] * pa_y[i];

        tr_z_xy_xyyyyy[i] = tr_z_y_yyyyy[i] * fe_0 + tr_z_y_xyyyyy[i] * pa_x[i];

        tr_z_xy_xyyyyz[i] = tr_z_y_yyyyz[i] * fe_0 + tr_z_y_xyyyyz[i] * pa_x[i];

        tr_z_xy_xyyyzz[i] = tr_z_y_yyyzz[i] * fe_0 + tr_z_y_xyyyzz[i] * pa_x[i];

        tr_z_xy_xyyzzz[i] = tr_z_y_yyzzz[i] * fe_0 + tr_z_y_xyyzzz[i] * pa_x[i];

        tr_z_xy_xyzzzz[i] = tr_z_y_yzzzz[i] * fe_0 + tr_z_y_xyzzzz[i] * pa_x[i];

        tr_z_xy_xzzzzz[i] = tr_z_x_xzzzzz[i] * pa_y[i];

        tr_z_xy_yyyyyy[i] = tr_z_y_yyyyyy[i] * pa_x[i];

        tr_z_xy_yyyyyz[i] = tr_z_y_yyyyyz[i] * pa_x[i];

        tr_z_xy_yyyyzz[i] = tr_z_y_yyyyzz[i] * pa_x[i];

        tr_z_xy_yyyzzz[i] = tr_z_y_yyyzzz[i] * pa_x[i];

        tr_z_xy_yyzzzz[i] = tr_z_y_yyzzzz[i] * pa_x[i];

        tr_z_xy_yzzzzz[i] = tr_z_y_yzzzzz[i] * pa_x[i];

        tr_z_xy_zzzzzz[i] = tr_z_y_zzzzzz[i] * pa_x[i];
    }

    // Set up 392-420 components of targeted buffer : DI

    auto tr_z_xz_xxxxxx = pbuffer.data(idx_dip_di + 392);

    auto tr_z_xz_xxxxxy = pbuffer.data(idx_dip_di + 393);

    auto tr_z_xz_xxxxxz = pbuffer.data(idx_dip_di + 394);

    auto tr_z_xz_xxxxyy = pbuffer.data(idx_dip_di + 395);

    auto tr_z_xz_xxxxyz = pbuffer.data(idx_dip_di + 396);

    auto tr_z_xz_xxxxzz = pbuffer.data(idx_dip_di + 397);

    auto tr_z_xz_xxxyyy = pbuffer.data(idx_dip_di + 398);

    auto tr_z_xz_xxxyyz = pbuffer.data(idx_dip_di + 399);

    auto tr_z_xz_xxxyzz = pbuffer.data(idx_dip_di + 400);

    auto tr_z_xz_xxxzzz = pbuffer.data(idx_dip_di + 401);

    auto tr_z_xz_xxyyyy = pbuffer.data(idx_dip_di + 402);

    auto tr_z_xz_xxyyyz = pbuffer.data(idx_dip_di + 403);

    auto tr_z_xz_xxyyzz = pbuffer.data(idx_dip_di + 404);

    auto tr_z_xz_xxyzzz = pbuffer.data(idx_dip_di + 405);

    auto tr_z_xz_xxzzzz = pbuffer.data(idx_dip_di + 406);

    auto tr_z_xz_xyyyyy = pbuffer.data(idx_dip_di + 407);

    auto tr_z_xz_xyyyyz = pbuffer.data(idx_dip_di + 408);

    auto tr_z_xz_xyyyzz = pbuffer.data(idx_dip_di + 409);

    auto tr_z_xz_xyyzzz = pbuffer.data(idx_dip_di + 410);

    auto tr_z_xz_xyzzzz = pbuffer.data(idx_dip_di + 411);

    auto tr_z_xz_xzzzzz = pbuffer.data(idx_dip_di + 412);

    auto tr_z_xz_yyyyyy = pbuffer.data(idx_dip_di + 413);

    auto tr_z_xz_yyyyyz = pbuffer.data(idx_dip_di + 414);

    auto tr_z_xz_yyyyzz = pbuffer.data(idx_dip_di + 415);

    auto tr_z_xz_yyyzzz = pbuffer.data(idx_dip_di + 416);

    auto tr_z_xz_yyzzzz = pbuffer.data(idx_dip_di + 417);

    auto tr_z_xz_yzzzzz = pbuffer.data(idx_dip_di + 418);

    auto tr_z_xz_zzzzzz = pbuffer.data(idx_dip_di + 419);

#pragma omp simd aligned(pa_x,               \
                             tr_z_xz_xxxxxx, \
                             tr_z_xz_xxxxxy, \
                             tr_z_xz_xxxxxz, \
                             tr_z_xz_xxxxyy, \
                             tr_z_xz_xxxxyz, \
                             tr_z_xz_xxxxzz, \
                             tr_z_xz_xxxyyy, \
                             tr_z_xz_xxxyyz, \
                             tr_z_xz_xxxyzz, \
                             tr_z_xz_xxxzzz, \
                             tr_z_xz_xxyyyy, \
                             tr_z_xz_xxyyyz, \
                             tr_z_xz_xxyyzz, \
                             tr_z_xz_xxyzzz, \
                             tr_z_xz_xxzzzz, \
                             tr_z_xz_xyyyyy, \
                             tr_z_xz_xyyyyz, \
                             tr_z_xz_xyyyzz, \
                             tr_z_xz_xyyzzz, \
                             tr_z_xz_xyzzzz, \
                             tr_z_xz_xzzzzz, \
                             tr_z_xz_yyyyyy, \
                             tr_z_xz_yyyyyz, \
                             tr_z_xz_yyyyzz, \
                             tr_z_xz_yyyzzz, \
                             tr_z_xz_yyzzzz, \
                             tr_z_xz_yzzzzz, \
                             tr_z_xz_zzzzzz, \
                             tr_z_z_xxxxx,   \
                             tr_z_z_xxxxxx,  \
                             tr_z_z_xxxxxy,  \
                             tr_z_z_xxxxxz,  \
                             tr_z_z_xxxxy,   \
                             tr_z_z_xxxxyy,  \
                             tr_z_z_xxxxyz,  \
                             tr_z_z_xxxxz,   \
                             tr_z_z_xxxxzz,  \
                             tr_z_z_xxxyy,   \
                             tr_z_z_xxxyyy,  \
                             tr_z_z_xxxyyz,  \
                             tr_z_z_xxxyz,   \
                             tr_z_z_xxxyzz,  \
                             tr_z_z_xxxzz,   \
                             tr_z_z_xxxzzz,  \
                             tr_z_z_xxyyy,   \
                             tr_z_z_xxyyyy,  \
                             tr_z_z_xxyyyz,  \
                             tr_z_z_xxyyz,   \
                             tr_z_z_xxyyzz,  \
                             tr_z_z_xxyzz,   \
                             tr_z_z_xxyzzz,  \
                             tr_z_z_xxzzz,   \
                             tr_z_z_xxzzzz,  \
                             tr_z_z_xyyyy,   \
                             tr_z_z_xyyyyy,  \
                             tr_z_z_xyyyyz,  \
                             tr_z_z_xyyyz,   \
                             tr_z_z_xyyyzz,  \
                             tr_z_z_xyyzz,   \
                             tr_z_z_xyyzzz,  \
                             tr_z_z_xyzzz,   \
                             tr_z_z_xyzzzz,  \
                             tr_z_z_xzzzz,   \
                             tr_z_z_xzzzzz,  \
                             tr_z_z_yyyyy,   \
                             tr_z_z_yyyyyy,  \
                             tr_z_z_yyyyyz,  \
                             tr_z_z_yyyyz,   \
                             tr_z_z_yyyyzz,  \
                             tr_z_z_yyyzz,   \
                             tr_z_z_yyyzzz,  \
                             tr_z_z_yyzzz,   \
                             tr_z_z_yyzzzz,  \
                             tr_z_z_yzzzz,   \
                             tr_z_z_yzzzzz,  \
                             tr_z_z_zzzzz,   \
                             tr_z_z_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_xz_xxxxxx[i] = 6.0 * tr_z_z_xxxxx[i] * fe_0 + tr_z_z_xxxxxx[i] * pa_x[i];

        tr_z_xz_xxxxxy[i] = 5.0 * tr_z_z_xxxxy[i] * fe_0 + tr_z_z_xxxxxy[i] * pa_x[i];

        tr_z_xz_xxxxxz[i] = 5.0 * tr_z_z_xxxxz[i] * fe_0 + tr_z_z_xxxxxz[i] * pa_x[i];

        tr_z_xz_xxxxyy[i] = 4.0 * tr_z_z_xxxyy[i] * fe_0 + tr_z_z_xxxxyy[i] * pa_x[i];

        tr_z_xz_xxxxyz[i] = 4.0 * tr_z_z_xxxyz[i] * fe_0 + tr_z_z_xxxxyz[i] * pa_x[i];

        tr_z_xz_xxxxzz[i] = 4.0 * tr_z_z_xxxzz[i] * fe_0 + tr_z_z_xxxxzz[i] * pa_x[i];

        tr_z_xz_xxxyyy[i] = 3.0 * tr_z_z_xxyyy[i] * fe_0 + tr_z_z_xxxyyy[i] * pa_x[i];

        tr_z_xz_xxxyyz[i] = 3.0 * tr_z_z_xxyyz[i] * fe_0 + tr_z_z_xxxyyz[i] * pa_x[i];

        tr_z_xz_xxxyzz[i] = 3.0 * tr_z_z_xxyzz[i] * fe_0 + tr_z_z_xxxyzz[i] * pa_x[i];

        tr_z_xz_xxxzzz[i] = 3.0 * tr_z_z_xxzzz[i] * fe_0 + tr_z_z_xxxzzz[i] * pa_x[i];

        tr_z_xz_xxyyyy[i] = 2.0 * tr_z_z_xyyyy[i] * fe_0 + tr_z_z_xxyyyy[i] * pa_x[i];

        tr_z_xz_xxyyyz[i] = 2.0 * tr_z_z_xyyyz[i] * fe_0 + tr_z_z_xxyyyz[i] * pa_x[i];

        tr_z_xz_xxyyzz[i] = 2.0 * tr_z_z_xyyzz[i] * fe_0 + tr_z_z_xxyyzz[i] * pa_x[i];

        tr_z_xz_xxyzzz[i] = 2.0 * tr_z_z_xyzzz[i] * fe_0 + tr_z_z_xxyzzz[i] * pa_x[i];

        tr_z_xz_xxzzzz[i] = 2.0 * tr_z_z_xzzzz[i] * fe_0 + tr_z_z_xxzzzz[i] * pa_x[i];

        tr_z_xz_xyyyyy[i] = tr_z_z_yyyyy[i] * fe_0 + tr_z_z_xyyyyy[i] * pa_x[i];

        tr_z_xz_xyyyyz[i] = tr_z_z_yyyyz[i] * fe_0 + tr_z_z_xyyyyz[i] * pa_x[i];

        tr_z_xz_xyyyzz[i] = tr_z_z_yyyzz[i] * fe_0 + tr_z_z_xyyyzz[i] * pa_x[i];

        tr_z_xz_xyyzzz[i] = tr_z_z_yyzzz[i] * fe_0 + tr_z_z_xyyzzz[i] * pa_x[i];

        tr_z_xz_xyzzzz[i] = tr_z_z_yzzzz[i] * fe_0 + tr_z_z_xyzzzz[i] * pa_x[i];

        tr_z_xz_xzzzzz[i] = tr_z_z_zzzzz[i] * fe_0 + tr_z_z_xzzzzz[i] * pa_x[i];

        tr_z_xz_yyyyyy[i] = tr_z_z_yyyyyy[i] * pa_x[i];

        tr_z_xz_yyyyyz[i] = tr_z_z_yyyyyz[i] * pa_x[i];

        tr_z_xz_yyyyzz[i] = tr_z_z_yyyyzz[i] * pa_x[i];

        tr_z_xz_yyyzzz[i] = tr_z_z_yyyzzz[i] * pa_x[i];

        tr_z_xz_yyzzzz[i] = tr_z_z_yyzzzz[i] * pa_x[i];

        tr_z_xz_yzzzzz[i] = tr_z_z_yzzzzz[i] * pa_x[i];

        tr_z_xz_zzzzzz[i] = tr_z_z_zzzzzz[i] * pa_x[i];
    }

    // Set up 420-448 components of targeted buffer : DI

    auto tr_z_yy_xxxxxx = pbuffer.data(idx_dip_di + 420);

    auto tr_z_yy_xxxxxy = pbuffer.data(idx_dip_di + 421);

    auto tr_z_yy_xxxxxz = pbuffer.data(idx_dip_di + 422);

    auto tr_z_yy_xxxxyy = pbuffer.data(idx_dip_di + 423);

    auto tr_z_yy_xxxxyz = pbuffer.data(idx_dip_di + 424);

    auto tr_z_yy_xxxxzz = pbuffer.data(idx_dip_di + 425);

    auto tr_z_yy_xxxyyy = pbuffer.data(idx_dip_di + 426);

    auto tr_z_yy_xxxyyz = pbuffer.data(idx_dip_di + 427);

    auto tr_z_yy_xxxyzz = pbuffer.data(idx_dip_di + 428);

    auto tr_z_yy_xxxzzz = pbuffer.data(idx_dip_di + 429);

    auto tr_z_yy_xxyyyy = pbuffer.data(idx_dip_di + 430);

    auto tr_z_yy_xxyyyz = pbuffer.data(idx_dip_di + 431);

    auto tr_z_yy_xxyyzz = pbuffer.data(idx_dip_di + 432);

    auto tr_z_yy_xxyzzz = pbuffer.data(idx_dip_di + 433);

    auto tr_z_yy_xxzzzz = pbuffer.data(idx_dip_di + 434);

    auto tr_z_yy_xyyyyy = pbuffer.data(idx_dip_di + 435);

    auto tr_z_yy_xyyyyz = pbuffer.data(idx_dip_di + 436);

    auto tr_z_yy_xyyyzz = pbuffer.data(idx_dip_di + 437);

    auto tr_z_yy_xyyzzz = pbuffer.data(idx_dip_di + 438);

    auto tr_z_yy_xyzzzz = pbuffer.data(idx_dip_di + 439);

    auto tr_z_yy_xzzzzz = pbuffer.data(idx_dip_di + 440);

    auto tr_z_yy_yyyyyy = pbuffer.data(idx_dip_di + 441);

    auto tr_z_yy_yyyyyz = pbuffer.data(idx_dip_di + 442);

    auto tr_z_yy_yyyyzz = pbuffer.data(idx_dip_di + 443);

    auto tr_z_yy_yyyzzz = pbuffer.data(idx_dip_di + 444);

    auto tr_z_yy_yyzzzz = pbuffer.data(idx_dip_di + 445);

    auto tr_z_yy_yzzzzz = pbuffer.data(idx_dip_di + 446);

    auto tr_z_yy_zzzzzz = pbuffer.data(idx_dip_di + 447);

#pragma omp simd aligned(pa_y,               \
                             tr_z_0_xxxxxx,  \
                             tr_z_0_xxxxxy,  \
                             tr_z_0_xxxxxz,  \
                             tr_z_0_xxxxyy,  \
                             tr_z_0_xxxxyz,  \
                             tr_z_0_xxxxzz,  \
                             tr_z_0_xxxyyy,  \
                             tr_z_0_xxxyyz,  \
                             tr_z_0_xxxyzz,  \
                             tr_z_0_xxxzzz,  \
                             tr_z_0_xxyyyy,  \
                             tr_z_0_xxyyyz,  \
                             tr_z_0_xxyyzz,  \
                             tr_z_0_xxyzzz,  \
                             tr_z_0_xxzzzz,  \
                             tr_z_0_xyyyyy,  \
                             tr_z_0_xyyyyz,  \
                             tr_z_0_xyyyzz,  \
                             tr_z_0_xyyzzz,  \
                             tr_z_0_xyzzzz,  \
                             tr_z_0_xzzzzz,  \
                             tr_z_0_yyyyyy,  \
                             tr_z_0_yyyyyz,  \
                             tr_z_0_yyyyzz,  \
                             tr_z_0_yyyzzz,  \
                             tr_z_0_yyzzzz,  \
                             tr_z_0_yzzzzz,  \
                             tr_z_0_zzzzzz,  \
                             tr_z_y_xxxxx,   \
                             tr_z_y_xxxxxx,  \
                             tr_z_y_xxxxxy,  \
                             tr_z_y_xxxxxz,  \
                             tr_z_y_xxxxy,   \
                             tr_z_y_xxxxyy,  \
                             tr_z_y_xxxxyz,  \
                             tr_z_y_xxxxz,   \
                             tr_z_y_xxxxzz,  \
                             tr_z_y_xxxyy,   \
                             tr_z_y_xxxyyy,  \
                             tr_z_y_xxxyyz,  \
                             tr_z_y_xxxyz,   \
                             tr_z_y_xxxyzz,  \
                             tr_z_y_xxxzz,   \
                             tr_z_y_xxxzzz,  \
                             tr_z_y_xxyyy,   \
                             tr_z_y_xxyyyy,  \
                             tr_z_y_xxyyyz,  \
                             tr_z_y_xxyyz,   \
                             tr_z_y_xxyyzz,  \
                             tr_z_y_xxyzz,   \
                             tr_z_y_xxyzzz,  \
                             tr_z_y_xxzzz,   \
                             tr_z_y_xxzzzz,  \
                             tr_z_y_xyyyy,   \
                             tr_z_y_xyyyyy,  \
                             tr_z_y_xyyyyz,  \
                             tr_z_y_xyyyz,   \
                             tr_z_y_xyyyzz,  \
                             tr_z_y_xyyzz,   \
                             tr_z_y_xyyzzz,  \
                             tr_z_y_xyzzz,   \
                             tr_z_y_xyzzzz,  \
                             tr_z_y_xzzzz,   \
                             tr_z_y_xzzzzz,  \
                             tr_z_y_yyyyy,   \
                             tr_z_y_yyyyyy,  \
                             tr_z_y_yyyyyz,  \
                             tr_z_y_yyyyz,   \
                             tr_z_y_yyyyzz,  \
                             tr_z_y_yyyzz,   \
                             tr_z_y_yyyzzz,  \
                             tr_z_y_yyzzz,   \
                             tr_z_y_yyzzzz,  \
                             tr_z_y_yzzzz,   \
                             tr_z_y_yzzzzz,  \
                             tr_z_y_zzzzz,   \
                             tr_z_y_zzzzzz,  \
                             tr_z_yy_xxxxxx, \
                             tr_z_yy_xxxxxy, \
                             tr_z_yy_xxxxxz, \
                             tr_z_yy_xxxxyy, \
                             tr_z_yy_xxxxyz, \
                             tr_z_yy_xxxxzz, \
                             tr_z_yy_xxxyyy, \
                             tr_z_yy_xxxyyz, \
                             tr_z_yy_xxxyzz, \
                             tr_z_yy_xxxzzz, \
                             tr_z_yy_xxyyyy, \
                             tr_z_yy_xxyyyz, \
                             tr_z_yy_xxyyzz, \
                             tr_z_yy_xxyzzz, \
                             tr_z_yy_xxzzzz, \
                             tr_z_yy_xyyyyy, \
                             tr_z_yy_xyyyyz, \
                             tr_z_yy_xyyyzz, \
                             tr_z_yy_xyyzzz, \
                             tr_z_yy_xyzzzz, \
                             tr_z_yy_xzzzzz, \
                             tr_z_yy_yyyyyy, \
                             tr_z_yy_yyyyyz, \
                             tr_z_yy_yyyyzz, \
                             tr_z_yy_yyyzzz, \
                             tr_z_yy_yyzzzz, \
                             tr_z_yy_yzzzzz, \
                             tr_z_yy_zzzzzz, \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yy_xxxxxx[i] = tr_z_0_xxxxxx[i] * fe_0 + tr_z_y_xxxxxx[i] * pa_y[i];

        tr_z_yy_xxxxxy[i] = tr_z_0_xxxxxy[i] * fe_0 + tr_z_y_xxxxx[i] * fe_0 + tr_z_y_xxxxxy[i] * pa_y[i];

        tr_z_yy_xxxxxz[i] = tr_z_0_xxxxxz[i] * fe_0 + tr_z_y_xxxxxz[i] * pa_y[i];

        tr_z_yy_xxxxyy[i] = tr_z_0_xxxxyy[i] * fe_0 + 2.0 * tr_z_y_xxxxy[i] * fe_0 + tr_z_y_xxxxyy[i] * pa_y[i];

        tr_z_yy_xxxxyz[i] = tr_z_0_xxxxyz[i] * fe_0 + tr_z_y_xxxxz[i] * fe_0 + tr_z_y_xxxxyz[i] * pa_y[i];

        tr_z_yy_xxxxzz[i] = tr_z_0_xxxxzz[i] * fe_0 + tr_z_y_xxxxzz[i] * pa_y[i];

        tr_z_yy_xxxyyy[i] = tr_z_0_xxxyyy[i] * fe_0 + 3.0 * tr_z_y_xxxyy[i] * fe_0 + tr_z_y_xxxyyy[i] * pa_y[i];

        tr_z_yy_xxxyyz[i] = tr_z_0_xxxyyz[i] * fe_0 + 2.0 * tr_z_y_xxxyz[i] * fe_0 + tr_z_y_xxxyyz[i] * pa_y[i];

        tr_z_yy_xxxyzz[i] = tr_z_0_xxxyzz[i] * fe_0 + tr_z_y_xxxzz[i] * fe_0 + tr_z_y_xxxyzz[i] * pa_y[i];

        tr_z_yy_xxxzzz[i] = tr_z_0_xxxzzz[i] * fe_0 + tr_z_y_xxxzzz[i] * pa_y[i];

        tr_z_yy_xxyyyy[i] = tr_z_0_xxyyyy[i] * fe_0 + 4.0 * tr_z_y_xxyyy[i] * fe_0 + tr_z_y_xxyyyy[i] * pa_y[i];

        tr_z_yy_xxyyyz[i] = tr_z_0_xxyyyz[i] * fe_0 + 3.0 * tr_z_y_xxyyz[i] * fe_0 + tr_z_y_xxyyyz[i] * pa_y[i];

        tr_z_yy_xxyyzz[i] = tr_z_0_xxyyzz[i] * fe_0 + 2.0 * tr_z_y_xxyzz[i] * fe_0 + tr_z_y_xxyyzz[i] * pa_y[i];

        tr_z_yy_xxyzzz[i] = tr_z_0_xxyzzz[i] * fe_0 + tr_z_y_xxzzz[i] * fe_0 + tr_z_y_xxyzzz[i] * pa_y[i];

        tr_z_yy_xxzzzz[i] = tr_z_0_xxzzzz[i] * fe_0 + tr_z_y_xxzzzz[i] * pa_y[i];

        tr_z_yy_xyyyyy[i] = tr_z_0_xyyyyy[i] * fe_0 + 5.0 * tr_z_y_xyyyy[i] * fe_0 + tr_z_y_xyyyyy[i] * pa_y[i];

        tr_z_yy_xyyyyz[i] = tr_z_0_xyyyyz[i] * fe_0 + 4.0 * tr_z_y_xyyyz[i] * fe_0 + tr_z_y_xyyyyz[i] * pa_y[i];

        tr_z_yy_xyyyzz[i] = tr_z_0_xyyyzz[i] * fe_0 + 3.0 * tr_z_y_xyyzz[i] * fe_0 + tr_z_y_xyyyzz[i] * pa_y[i];

        tr_z_yy_xyyzzz[i] = tr_z_0_xyyzzz[i] * fe_0 + 2.0 * tr_z_y_xyzzz[i] * fe_0 + tr_z_y_xyyzzz[i] * pa_y[i];

        tr_z_yy_xyzzzz[i] = tr_z_0_xyzzzz[i] * fe_0 + tr_z_y_xzzzz[i] * fe_0 + tr_z_y_xyzzzz[i] * pa_y[i];

        tr_z_yy_xzzzzz[i] = tr_z_0_xzzzzz[i] * fe_0 + tr_z_y_xzzzzz[i] * pa_y[i];

        tr_z_yy_yyyyyy[i] = tr_z_0_yyyyyy[i] * fe_0 + 6.0 * tr_z_y_yyyyy[i] * fe_0 + tr_z_y_yyyyyy[i] * pa_y[i];

        tr_z_yy_yyyyyz[i] = tr_z_0_yyyyyz[i] * fe_0 + 5.0 * tr_z_y_yyyyz[i] * fe_0 + tr_z_y_yyyyyz[i] * pa_y[i];

        tr_z_yy_yyyyzz[i] = tr_z_0_yyyyzz[i] * fe_0 + 4.0 * tr_z_y_yyyzz[i] * fe_0 + tr_z_y_yyyyzz[i] * pa_y[i];

        tr_z_yy_yyyzzz[i] = tr_z_0_yyyzzz[i] * fe_0 + 3.0 * tr_z_y_yyzzz[i] * fe_0 + tr_z_y_yyyzzz[i] * pa_y[i];

        tr_z_yy_yyzzzz[i] = tr_z_0_yyzzzz[i] * fe_0 + 2.0 * tr_z_y_yzzzz[i] * fe_0 + tr_z_y_yyzzzz[i] * pa_y[i];

        tr_z_yy_yzzzzz[i] = tr_z_0_yzzzzz[i] * fe_0 + tr_z_y_zzzzz[i] * fe_0 + tr_z_y_yzzzzz[i] * pa_y[i];

        tr_z_yy_zzzzzz[i] = tr_z_0_zzzzzz[i] * fe_0 + tr_z_y_zzzzzz[i] * pa_y[i];
    }

    // Set up 448-476 components of targeted buffer : DI

    auto tr_z_yz_xxxxxx = pbuffer.data(idx_dip_di + 448);

    auto tr_z_yz_xxxxxy = pbuffer.data(idx_dip_di + 449);

    auto tr_z_yz_xxxxxz = pbuffer.data(idx_dip_di + 450);

    auto tr_z_yz_xxxxyy = pbuffer.data(idx_dip_di + 451);

    auto tr_z_yz_xxxxyz = pbuffer.data(idx_dip_di + 452);

    auto tr_z_yz_xxxxzz = pbuffer.data(idx_dip_di + 453);

    auto tr_z_yz_xxxyyy = pbuffer.data(idx_dip_di + 454);

    auto tr_z_yz_xxxyyz = pbuffer.data(idx_dip_di + 455);

    auto tr_z_yz_xxxyzz = pbuffer.data(idx_dip_di + 456);

    auto tr_z_yz_xxxzzz = pbuffer.data(idx_dip_di + 457);

    auto tr_z_yz_xxyyyy = pbuffer.data(idx_dip_di + 458);

    auto tr_z_yz_xxyyyz = pbuffer.data(idx_dip_di + 459);

    auto tr_z_yz_xxyyzz = pbuffer.data(idx_dip_di + 460);

    auto tr_z_yz_xxyzzz = pbuffer.data(idx_dip_di + 461);

    auto tr_z_yz_xxzzzz = pbuffer.data(idx_dip_di + 462);

    auto tr_z_yz_xyyyyy = pbuffer.data(idx_dip_di + 463);

    auto tr_z_yz_xyyyyz = pbuffer.data(idx_dip_di + 464);

    auto tr_z_yz_xyyyzz = pbuffer.data(idx_dip_di + 465);

    auto tr_z_yz_xyyzzz = pbuffer.data(idx_dip_di + 466);

    auto tr_z_yz_xyzzzz = pbuffer.data(idx_dip_di + 467);

    auto tr_z_yz_xzzzzz = pbuffer.data(idx_dip_di + 468);

    auto tr_z_yz_yyyyyy = pbuffer.data(idx_dip_di + 469);

    auto tr_z_yz_yyyyyz = pbuffer.data(idx_dip_di + 470);

    auto tr_z_yz_yyyyzz = pbuffer.data(idx_dip_di + 471);

    auto tr_z_yz_yyyzzz = pbuffer.data(idx_dip_di + 472);

    auto tr_z_yz_yyzzzz = pbuffer.data(idx_dip_di + 473);

    auto tr_z_yz_yzzzzz = pbuffer.data(idx_dip_di + 474);

    auto tr_z_yz_zzzzzz = pbuffer.data(idx_dip_di + 475);

#pragma omp simd aligned(pa_y,               \
                             tr_z_yz_xxxxxx, \
                             tr_z_yz_xxxxxy, \
                             tr_z_yz_xxxxxz, \
                             tr_z_yz_xxxxyy, \
                             tr_z_yz_xxxxyz, \
                             tr_z_yz_xxxxzz, \
                             tr_z_yz_xxxyyy, \
                             tr_z_yz_xxxyyz, \
                             tr_z_yz_xxxyzz, \
                             tr_z_yz_xxxzzz, \
                             tr_z_yz_xxyyyy, \
                             tr_z_yz_xxyyyz, \
                             tr_z_yz_xxyyzz, \
                             tr_z_yz_xxyzzz, \
                             tr_z_yz_xxzzzz, \
                             tr_z_yz_xyyyyy, \
                             tr_z_yz_xyyyyz, \
                             tr_z_yz_xyyyzz, \
                             tr_z_yz_xyyzzz, \
                             tr_z_yz_xyzzzz, \
                             tr_z_yz_xzzzzz, \
                             tr_z_yz_yyyyyy, \
                             tr_z_yz_yyyyyz, \
                             tr_z_yz_yyyyzz, \
                             tr_z_yz_yyyzzz, \
                             tr_z_yz_yyzzzz, \
                             tr_z_yz_yzzzzz, \
                             tr_z_yz_zzzzzz, \
                             tr_z_z_xxxxx,   \
                             tr_z_z_xxxxxx,  \
                             tr_z_z_xxxxxy,  \
                             tr_z_z_xxxxxz,  \
                             tr_z_z_xxxxy,   \
                             tr_z_z_xxxxyy,  \
                             tr_z_z_xxxxyz,  \
                             tr_z_z_xxxxz,   \
                             tr_z_z_xxxxzz,  \
                             tr_z_z_xxxyy,   \
                             tr_z_z_xxxyyy,  \
                             tr_z_z_xxxyyz,  \
                             tr_z_z_xxxyz,   \
                             tr_z_z_xxxyzz,  \
                             tr_z_z_xxxzz,   \
                             tr_z_z_xxxzzz,  \
                             tr_z_z_xxyyy,   \
                             tr_z_z_xxyyyy,  \
                             tr_z_z_xxyyyz,  \
                             tr_z_z_xxyyz,   \
                             tr_z_z_xxyyzz,  \
                             tr_z_z_xxyzz,   \
                             tr_z_z_xxyzzz,  \
                             tr_z_z_xxzzz,   \
                             tr_z_z_xxzzzz,  \
                             tr_z_z_xyyyy,   \
                             tr_z_z_xyyyyy,  \
                             tr_z_z_xyyyyz,  \
                             tr_z_z_xyyyz,   \
                             tr_z_z_xyyyzz,  \
                             tr_z_z_xyyzz,   \
                             tr_z_z_xyyzzz,  \
                             tr_z_z_xyzzz,   \
                             tr_z_z_xyzzzz,  \
                             tr_z_z_xzzzz,   \
                             tr_z_z_xzzzzz,  \
                             tr_z_z_yyyyy,   \
                             tr_z_z_yyyyyy,  \
                             tr_z_z_yyyyyz,  \
                             tr_z_z_yyyyz,   \
                             tr_z_z_yyyyzz,  \
                             tr_z_z_yyyzz,   \
                             tr_z_z_yyyzzz,  \
                             tr_z_z_yyzzz,   \
                             tr_z_z_yyzzzz,  \
                             tr_z_z_yzzzz,   \
                             tr_z_z_yzzzzz,  \
                             tr_z_z_zzzzz,   \
                             tr_z_z_zzzzzz,  \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_yz_xxxxxx[i] = tr_z_z_xxxxxx[i] * pa_y[i];

        tr_z_yz_xxxxxy[i] = tr_z_z_xxxxx[i] * fe_0 + tr_z_z_xxxxxy[i] * pa_y[i];

        tr_z_yz_xxxxxz[i] = tr_z_z_xxxxxz[i] * pa_y[i];

        tr_z_yz_xxxxyy[i] = 2.0 * tr_z_z_xxxxy[i] * fe_0 + tr_z_z_xxxxyy[i] * pa_y[i];

        tr_z_yz_xxxxyz[i] = tr_z_z_xxxxz[i] * fe_0 + tr_z_z_xxxxyz[i] * pa_y[i];

        tr_z_yz_xxxxzz[i] = tr_z_z_xxxxzz[i] * pa_y[i];

        tr_z_yz_xxxyyy[i] = 3.0 * tr_z_z_xxxyy[i] * fe_0 + tr_z_z_xxxyyy[i] * pa_y[i];

        tr_z_yz_xxxyyz[i] = 2.0 * tr_z_z_xxxyz[i] * fe_0 + tr_z_z_xxxyyz[i] * pa_y[i];

        tr_z_yz_xxxyzz[i] = tr_z_z_xxxzz[i] * fe_0 + tr_z_z_xxxyzz[i] * pa_y[i];

        tr_z_yz_xxxzzz[i] = tr_z_z_xxxzzz[i] * pa_y[i];

        tr_z_yz_xxyyyy[i] = 4.0 * tr_z_z_xxyyy[i] * fe_0 + tr_z_z_xxyyyy[i] * pa_y[i];

        tr_z_yz_xxyyyz[i] = 3.0 * tr_z_z_xxyyz[i] * fe_0 + tr_z_z_xxyyyz[i] * pa_y[i];

        tr_z_yz_xxyyzz[i] = 2.0 * tr_z_z_xxyzz[i] * fe_0 + tr_z_z_xxyyzz[i] * pa_y[i];

        tr_z_yz_xxyzzz[i] = tr_z_z_xxzzz[i] * fe_0 + tr_z_z_xxyzzz[i] * pa_y[i];

        tr_z_yz_xxzzzz[i] = tr_z_z_xxzzzz[i] * pa_y[i];

        tr_z_yz_xyyyyy[i] = 5.0 * tr_z_z_xyyyy[i] * fe_0 + tr_z_z_xyyyyy[i] * pa_y[i];

        tr_z_yz_xyyyyz[i] = 4.0 * tr_z_z_xyyyz[i] * fe_0 + tr_z_z_xyyyyz[i] * pa_y[i];

        tr_z_yz_xyyyzz[i] = 3.0 * tr_z_z_xyyzz[i] * fe_0 + tr_z_z_xyyyzz[i] * pa_y[i];

        tr_z_yz_xyyzzz[i] = 2.0 * tr_z_z_xyzzz[i] * fe_0 + tr_z_z_xyyzzz[i] * pa_y[i];

        tr_z_yz_xyzzzz[i] = tr_z_z_xzzzz[i] * fe_0 + tr_z_z_xyzzzz[i] * pa_y[i];

        tr_z_yz_xzzzzz[i] = tr_z_z_xzzzzz[i] * pa_y[i];

        tr_z_yz_yyyyyy[i] = 6.0 * tr_z_z_yyyyy[i] * fe_0 + tr_z_z_yyyyyy[i] * pa_y[i];

        tr_z_yz_yyyyyz[i] = 5.0 * tr_z_z_yyyyz[i] * fe_0 + tr_z_z_yyyyyz[i] * pa_y[i];

        tr_z_yz_yyyyzz[i] = 4.0 * tr_z_z_yyyzz[i] * fe_0 + tr_z_z_yyyyzz[i] * pa_y[i];

        tr_z_yz_yyyzzz[i] = 3.0 * tr_z_z_yyzzz[i] * fe_0 + tr_z_z_yyyzzz[i] * pa_y[i];

        tr_z_yz_yyzzzz[i] = 2.0 * tr_z_z_yzzzz[i] * fe_0 + tr_z_z_yyzzzz[i] * pa_y[i];

        tr_z_yz_yzzzzz[i] = tr_z_z_zzzzz[i] * fe_0 + tr_z_z_yzzzzz[i] * pa_y[i];

        tr_z_yz_zzzzzz[i] = tr_z_z_zzzzzz[i] * pa_y[i];
    }

    // Set up 476-504 components of targeted buffer : DI

    auto tr_z_zz_xxxxxx = pbuffer.data(idx_dip_di + 476);

    auto tr_z_zz_xxxxxy = pbuffer.data(idx_dip_di + 477);

    auto tr_z_zz_xxxxxz = pbuffer.data(idx_dip_di + 478);

    auto tr_z_zz_xxxxyy = pbuffer.data(idx_dip_di + 479);

    auto tr_z_zz_xxxxyz = pbuffer.data(idx_dip_di + 480);

    auto tr_z_zz_xxxxzz = pbuffer.data(idx_dip_di + 481);

    auto tr_z_zz_xxxyyy = pbuffer.data(idx_dip_di + 482);

    auto tr_z_zz_xxxyyz = pbuffer.data(idx_dip_di + 483);

    auto tr_z_zz_xxxyzz = pbuffer.data(idx_dip_di + 484);

    auto tr_z_zz_xxxzzz = pbuffer.data(idx_dip_di + 485);

    auto tr_z_zz_xxyyyy = pbuffer.data(idx_dip_di + 486);

    auto tr_z_zz_xxyyyz = pbuffer.data(idx_dip_di + 487);

    auto tr_z_zz_xxyyzz = pbuffer.data(idx_dip_di + 488);

    auto tr_z_zz_xxyzzz = pbuffer.data(idx_dip_di + 489);

    auto tr_z_zz_xxzzzz = pbuffer.data(idx_dip_di + 490);

    auto tr_z_zz_xyyyyy = pbuffer.data(idx_dip_di + 491);

    auto tr_z_zz_xyyyyz = pbuffer.data(idx_dip_di + 492);

    auto tr_z_zz_xyyyzz = pbuffer.data(idx_dip_di + 493);

    auto tr_z_zz_xyyzzz = pbuffer.data(idx_dip_di + 494);

    auto tr_z_zz_xyzzzz = pbuffer.data(idx_dip_di + 495);

    auto tr_z_zz_xzzzzz = pbuffer.data(idx_dip_di + 496);

    auto tr_z_zz_yyyyyy = pbuffer.data(idx_dip_di + 497);

    auto tr_z_zz_yyyyyz = pbuffer.data(idx_dip_di + 498);

    auto tr_z_zz_yyyyzz = pbuffer.data(idx_dip_di + 499);

    auto tr_z_zz_yyyzzz = pbuffer.data(idx_dip_di + 500);

    auto tr_z_zz_yyzzzz = pbuffer.data(idx_dip_di + 501);

    auto tr_z_zz_yzzzzz = pbuffer.data(idx_dip_di + 502);

    auto tr_z_zz_zzzzzz = pbuffer.data(idx_dip_di + 503);

#pragma omp simd aligned(pa_z,               \
                             tr_z_0_xxxxxx,  \
                             tr_z_0_xxxxxy,  \
                             tr_z_0_xxxxxz,  \
                             tr_z_0_xxxxyy,  \
                             tr_z_0_xxxxyz,  \
                             tr_z_0_xxxxzz,  \
                             tr_z_0_xxxyyy,  \
                             tr_z_0_xxxyyz,  \
                             tr_z_0_xxxyzz,  \
                             tr_z_0_xxxzzz,  \
                             tr_z_0_xxyyyy,  \
                             tr_z_0_xxyyyz,  \
                             tr_z_0_xxyyzz,  \
                             tr_z_0_xxyzzz,  \
                             tr_z_0_xxzzzz,  \
                             tr_z_0_xyyyyy,  \
                             tr_z_0_xyyyyz,  \
                             tr_z_0_xyyyzz,  \
                             tr_z_0_xyyzzz,  \
                             tr_z_0_xyzzzz,  \
                             tr_z_0_xzzzzz,  \
                             tr_z_0_yyyyyy,  \
                             tr_z_0_yyyyyz,  \
                             tr_z_0_yyyyzz,  \
                             tr_z_0_yyyzzz,  \
                             tr_z_0_yyzzzz,  \
                             tr_z_0_yzzzzz,  \
                             tr_z_0_zzzzzz,  \
                             tr_z_z_xxxxx,   \
                             tr_z_z_xxxxxx,  \
                             tr_z_z_xxxxxy,  \
                             tr_z_z_xxxxxz,  \
                             tr_z_z_xxxxy,   \
                             tr_z_z_xxxxyy,  \
                             tr_z_z_xxxxyz,  \
                             tr_z_z_xxxxz,   \
                             tr_z_z_xxxxzz,  \
                             tr_z_z_xxxyy,   \
                             tr_z_z_xxxyyy,  \
                             tr_z_z_xxxyyz,  \
                             tr_z_z_xxxyz,   \
                             tr_z_z_xxxyzz,  \
                             tr_z_z_xxxzz,   \
                             tr_z_z_xxxzzz,  \
                             tr_z_z_xxyyy,   \
                             tr_z_z_xxyyyy,  \
                             tr_z_z_xxyyyz,  \
                             tr_z_z_xxyyz,   \
                             tr_z_z_xxyyzz,  \
                             tr_z_z_xxyzz,   \
                             tr_z_z_xxyzzz,  \
                             tr_z_z_xxzzz,   \
                             tr_z_z_xxzzzz,  \
                             tr_z_z_xyyyy,   \
                             tr_z_z_xyyyyy,  \
                             tr_z_z_xyyyyz,  \
                             tr_z_z_xyyyz,   \
                             tr_z_z_xyyyzz,  \
                             tr_z_z_xyyzz,   \
                             tr_z_z_xyyzzz,  \
                             tr_z_z_xyzzz,   \
                             tr_z_z_xyzzzz,  \
                             tr_z_z_xzzzz,   \
                             tr_z_z_xzzzzz,  \
                             tr_z_z_yyyyy,   \
                             tr_z_z_yyyyyy,  \
                             tr_z_z_yyyyyz,  \
                             tr_z_z_yyyyz,   \
                             tr_z_z_yyyyzz,  \
                             tr_z_z_yyyzz,   \
                             tr_z_z_yyyzzz,  \
                             tr_z_z_yyzzz,   \
                             tr_z_z_yyzzzz,  \
                             tr_z_z_yzzzz,   \
                             tr_z_z_yzzzzz,  \
                             tr_z_z_zzzzz,   \
                             tr_z_z_zzzzzz,  \
                             tr_z_zz_xxxxxx, \
                             tr_z_zz_xxxxxy, \
                             tr_z_zz_xxxxxz, \
                             tr_z_zz_xxxxyy, \
                             tr_z_zz_xxxxyz, \
                             tr_z_zz_xxxxzz, \
                             tr_z_zz_xxxyyy, \
                             tr_z_zz_xxxyyz, \
                             tr_z_zz_xxxyzz, \
                             tr_z_zz_xxxzzz, \
                             tr_z_zz_xxyyyy, \
                             tr_z_zz_xxyyyz, \
                             tr_z_zz_xxyyzz, \
                             tr_z_zz_xxyzzz, \
                             tr_z_zz_xxzzzz, \
                             tr_z_zz_xyyyyy, \
                             tr_z_zz_xyyyyz, \
                             tr_z_zz_xyyyzz, \
                             tr_z_zz_xyyzzz, \
                             tr_z_zz_xyzzzz, \
                             tr_z_zz_xzzzzz, \
                             tr_z_zz_yyyyyy, \
                             tr_z_zz_yyyyyz, \
                             tr_z_zz_yyyyzz, \
                             tr_z_zz_yyyzzz, \
                             tr_z_zz_yyzzzz, \
                             tr_z_zz_yzzzzz, \
                             tr_z_zz_zzzzzz, \
                             ts_z_xxxxxx,    \
                             ts_z_xxxxxy,    \
                             ts_z_xxxxxz,    \
                             ts_z_xxxxyy,    \
                             ts_z_xxxxyz,    \
                             ts_z_xxxxzz,    \
                             ts_z_xxxyyy,    \
                             ts_z_xxxyyz,    \
                             ts_z_xxxyzz,    \
                             ts_z_xxxzzz,    \
                             ts_z_xxyyyy,    \
                             ts_z_xxyyyz,    \
                             ts_z_xxyyzz,    \
                             ts_z_xxyzzz,    \
                             ts_z_xxzzzz,    \
                             ts_z_xyyyyy,    \
                             ts_z_xyyyyz,    \
                             ts_z_xyyyzz,    \
                             ts_z_xyyzzz,    \
                             ts_z_xyzzzz,    \
                             ts_z_xzzzzz,    \
                             ts_z_yyyyyy,    \
                             ts_z_yyyyyz,    \
                             ts_z_yyyyzz,    \
                             ts_z_yyyzzz,    \
                             ts_z_yyzzzz,    \
                             ts_z_yzzzzz,    \
                             ts_z_zzzzzz,    \
                             b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const double fe_0 = 0.5 / (a_exp + b_exps[i]);

        tr_z_zz_xxxxxx[i] = tr_z_0_xxxxxx[i] * fe_0 + ts_z_xxxxxx[i] * fe_0 + tr_z_z_xxxxxx[i] * pa_z[i];

        tr_z_zz_xxxxxy[i] = tr_z_0_xxxxxy[i] * fe_0 + ts_z_xxxxxy[i] * fe_0 + tr_z_z_xxxxxy[i] * pa_z[i];

        tr_z_zz_xxxxxz[i] = tr_z_0_xxxxxz[i] * fe_0 + tr_z_z_xxxxx[i] * fe_0 + ts_z_xxxxxz[i] * fe_0 + tr_z_z_xxxxxz[i] * pa_z[i];

        tr_z_zz_xxxxyy[i] = tr_z_0_xxxxyy[i] * fe_0 + ts_z_xxxxyy[i] * fe_0 + tr_z_z_xxxxyy[i] * pa_z[i];

        tr_z_zz_xxxxyz[i] = tr_z_0_xxxxyz[i] * fe_0 + tr_z_z_xxxxy[i] * fe_0 + ts_z_xxxxyz[i] * fe_0 + tr_z_z_xxxxyz[i] * pa_z[i];

        tr_z_zz_xxxxzz[i] = tr_z_0_xxxxzz[i] * fe_0 + 2.0 * tr_z_z_xxxxz[i] * fe_0 + ts_z_xxxxzz[i] * fe_0 + tr_z_z_xxxxzz[i] * pa_z[i];

        tr_z_zz_xxxyyy[i] = tr_z_0_xxxyyy[i] * fe_0 + ts_z_xxxyyy[i] * fe_0 + tr_z_z_xxxyyy[i] * pa_z[i];

        tr_z_zz_xxxyyz[i] = tr_z_0_xxxyyz[i] * fe_0 + tr_z_z_xxxyy[i] * fe_0 + ts_z_xxxyyz[i] * fe_0 + tr_z_z_xxxyyz[i] * pa_z[i];

        tr_z_zz_xxxyzz[i] = tr_z_0_xxxyzz[i] * fe_0 + 2.0 * tr_z_z_xxxyz[i] * fe_0 + ts_z_xxxyzz[i] * fe_0 + tr_z_z_xxxyzz[i] * pa_z[i];

        tr_z_zz_xxxzzz[i] = tr_z_0_xxxzzz[i] * fe_0 + 3.0 * tr_z_z_xxxzz[i] * fe_0 + ts_z_xxxzzz[i] * fe_0 + tr_z_z_xxxzzz[i] * pa_z[i];

        tr_z_zz_xxyyyy[i] = tr_z_0_xxyyyy[i] * fe_0 + ts_z_xxyyyy[i] * fe_0 + tr_z_z_xxyyyy[i] * pa_z[i];

        tr_z_zz_xxyyyz[i] = tr_z_0_xxyyyz[i] * fe_0 + tr_z_z_xxyyy[i] * fe_0 + ts_z_xxyyyz[i] * fe_0 + tr_z_z_xxyyyz[i] * pa_z[i];

        tr_z_zz_xxyyzz[i] = tr_z_0_xxyyzz[i] * fe_0 + 2.0 * tr_z_z_xxyyz[i] * fe_0 + ts_z_xxyyzz[i] * fe_0 + tr_z_z_xxyyzz[i] * pa_z[i];

        tr_z_zz_xxyzzz[i] = tr_z_0_xxyzzz[i] * fe_0 + 3.0 * tr_z_z_xxyzz[i] * fe_0 + ts_z_xxyzzz[i] * fe_0 + tr_z_z_xxyzzz[i] * pa_z[i];

        tr_z_zz_xxzzzz[i] = tr_z_0_xxzzzz[i] * fe_0 + 4.0 * tr_z_z_xxzzz[i] * fe_0 + ts_z_xxzzzz[i] * fe_0 + tr_z_z_xxzzzz[i] * pa_z[i];

        tr_z_zz_xyyyyy[i] = tr_z_0_xyyyyy[i] * fe_0 + ts_z_xyyyyy[i] * fe_0 + tr_z_z_xyyyyy[i] * pa_z[i];

        tr_z_zz_xyyyyz[i] = tr_z_0_xyyyyz[i] * fe_0 + tr_z_z_xyyyy[i] * fe_0 + ts_z_xyyyyz[i] * fe_0 + tr_z_z_xyyyyz[i] * pa_z[i];

        tr_z_zz_xyyyzz[i] = tr_z_0_xyyyzz[i] * fe_0 + 2.0 * tr_z_z_xyyyz[i] * fe_0 + ts_z_xyyyzz[i] * fe_0 + tr_z_z_xyyyzz[i] * pa_z[i];

        tr_z_zz_xyyzzz[i] = tr_z_0_xyyzzz[i] * fe_0 + 3.0 * tr_z_z_xyyzz[i] * fe_0 + ts_z_xyyzzz[i] * fe_0 + tr_z_z_xyyzzz[i] * pa_z[i];

        tr_z_zz_xyzzzz[i] = tr_z_0_xyzzzz[i] * fe_0 + 4.0 * tr_z_z_xyzzz[i] * fe_0 + ts_z_xyzzzz[i] * fe_0 + tr_z_z_xyzzzz[i] * pa_z[i];

        tr_z_zz_xzzzzz[i] = tr_z_0_xzzzzz[i] * fe_0 + 5.0 * tr_z_z_xzzzz[i] * fe_0 + ts_z_xzzzzz[i] * fe_0 + tr_z_z_xzzzzz[i] * pa_z[i];

        tr_z_zz_yyyyyy[i] = tr_z_0_yyyyyy[i] * fe_0 + ts_z_yyyyyy[i] * fe_0 + tr_z_z_yyyyyy[i] * pa_z[i];

        tr_z_zz_yyyyyz[i] = tr_z_0_yyyyyz[i] * fe_0 + tr_z_z_yyyyy[i] * fe_0 + ts_z_yyyyyz[i] * fe_0 + tr_z_z_yyyyyz[i] * pa_z[i];

        tr_z_zz_yyyyzz[i] = tr_z_0_yyyyzz[i] * fe_0 + 2.0 * tr_z_z_yyyyz[i] * fe_0 + ts_z_yyyyzz[i] * fe_0 + tr_z_z_yyyyzz[i] * pa_z[i];

        tr_z_zz_yyyzzz[i] = tr_z_0_yyyzzz[i] * fe_0 + 3.0 * tr_z_z_yyyzz[i] * fe_0 + ts_z_yyyzzz[i] * fe_0 + tr_z_z_yyyzzz[i] * pa_z[i];

        tr_z_zz_yyzzzz[i] = tr_z_0_yyzzzz[i] * fe_0 + 4.0 * tr_z_z_yyzzz[i] * fe_0 + ts_z_yyzzzz[i] * fe_0 + tr_z_z_yyzzzz[i] * pa_z[i];

        tr_z_zz_yzzzzz[i] = tr_z_0_yzzzzz[i] * fe_0 + 5.0 * tr_z_z_yzzzz[i] * fe_0 + ts_z_yzzzzz[i] * fe_0 + tr_z_z_yzzzzz[i] * pa_z[i];

        tr_z_zz_zzzzzz[i] = tr_z_0_zzzzzz[i] * fe_0 + 6.0 * tr_z_z_zzzzz[i] * fe_0 + ts_z_zzzzzz[i] * fe_0 + tr_z_z_zzzzzz[i] * pa_z[i];
    }
}

}  // namespace diprec
